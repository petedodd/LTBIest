## Analysis code for Houben & Dodd 2016, distributed under CC BY 4.0 license https://creativecommons.org/licenses/by/4.0/
library(data.table)
library(reshape2)
load('data/whokey.Rdata')
load('data/All_3.Rdata')                #ARI data
load('data/POP2014.Rdata')                   #population data
load('data/POP2035.Rdata')
load('data/POP2050.Rdata')
load('data/POP1997.Rdata')
load('data/DSN2.Rdata') #INH data from Dodd,Sismanidis,Seddon

All <- All[All$lari!=-Inf,]
cnz <- unique(as.character(All$iso3))

## need to have run GP regression first to generate this data
RUNZ <- BDZ <- list()
for(i in 1:length(cnz)){
    cn <- cnz[i]

    fn <- paste0('data/',cn,'.Rdata')
    if(file.exists(fn)){
        load(fn)
        BDZ[[i]] <- erw
    }
    
    fn <- paste0('data/zz_',cn,'.Rdata')
    if(file.exists(fn)){
        load(fn)
        RUNZ[[i]] <- runsdf
    }
}

bdzdf <- do.call('rbind',BDZ)
rundata <- do.call('rbind',RUNZ)


## ------------- analyse data -------------

rundata <- data.table(rundata)
rund <- data.table(rundata)

## rundata
rundata[,year:=2014-year]               #now age
rundata[,lari:=exp(lari)]               #now real ari
rundata <- rundata[order(replicate,iso3,year),list(ari=lari,H=cumsum(lari),year=year),by=list(iso3=iso3,replicate=replicate)]

## ## for past 2 years
mask <- rep(1,length(unique(rundata$year)))
mask[1:2] <- 0                          #all except last 2 years
rundata[,dH:=cumsum(ari*mask),by=list(iso3=iso3,replicate=replicate)] #cumhaz!2y
rundata[,acat:=year %/% 5]
rundata[,P:=1-exp(-H)]                  #ever
rundata[,P1:=-exp(-H)+exp(-dH)]                  #1st recent=prob ever - prob not<2

## CHANGE HERE SENSE!
## Andrews: 0.79 .7-.86
pm <- 0.79                              #0.5 #CHANGE HERE!
pv <- (0.86-0.7)^2/3.92^2
apb <- pm*(1-pm)/pv-1
pa <- pm*apb                            #77.88
pb <- (1-pm)*apb                        #20.70
## curve(dbeta(x,shape1 = pa,shape2=pb),from=0,to=1)
## abline(v=pm,col=2);abline(v=.86,col=2,lty=2);abline(v=.7,col=2,lty=2);
## swap
alph <- rbeta(nrow(rundata),shape1=pb,shape2=pa)

rundata[,P2:=alph*(H-dH) + (1-alph)*(exp(-dH)-exp(-H))]                  #anyrecent

rundata <- rundata[,list(P=mean(P),P1=mean(P1),P2=mean(P2)),
                   by=list(iso3,replicate,acat)]
rundata <- merge(rundata,WHOkey[,c('iso3','g_whoregion')],by='iso3',all.x=TRUE)


POP2014$age <- factor(POP2014$age)
POP2014$year <- as.numeric(POP2014$age)-1
colnames(POP2014)[4] <- 'pop'
colnames(POP2014)[5] <- 'acat'

rundata <- merge(rundata,POP2014[,list(iso3,acat,pop)],by=c('iso3','acat'),all.x=TRUE)

## drop NAs
drpd <- as.character(rundata[is.na(pop),unique(iso3)]) #population data mismatch - 24
rundata <- rundata[!is.na(pop),]        #9K record the countries

## country size threshold 
popsizes <- rundata[replicate==1,list(pop=sum(pop)),by=iso3]
popsizes[,sum(pop<5e2)]                 #24
drp2 <- as.character(popsizes[pop<5e2,iso3])
rundata <- rundata[!iso3 %in% drp2,]
rundata[,length(unique(iso3))]          #169

## drop for having fewer than 15 datapoints
tbl <- table(All$iso3)
drp3 <- names(tbl)[which(tbl<15)][-1]   #4 countries
rundata <- rundata[!iso3 %in% drp3,]
rundata[,length(unique(iso3))]          #168

##  INH
DSN$INH <- DSN$INH + DSN$MDR         #include MDR with the INH-MR
DSN <- DSN[,c('iso3','iter','INH')]
DSN <- DSN[DSN$iso3 %in% as.character(unique(rundata$iso3)),]
smp <- sample(1e3,200)
DSN <- DSN[DSN$iter %in% smp,]            #sample 200
DSN$replicate <- as.numeric(factor(DSN$iter))

## merge in INH data here
rundata <- merge( rundata, DSN[,-2], by=c('iso3','replicate'), all.x=TRUE)

## summary by age by region
rundata0 <- rundata[,list(LTBI=sum(pop*P)/sum(pop),
                          LTBI1=sum(pop*P1)/sum(pop),
                          LTBI2=sum(pop*P2)/sum(pop),
                          LTBIH=sum(pop*P2*INH)/sum(pop)),
                    by=list(acat,replicate,g_whoregion)]

uq <- function(x)quantile(x,probs = .75)
lq <- function(x)quantile(x,probs = .25)

rundata0 <- rundata0[,list(LTBI=median(LTBI),LTBI1=median(LTBI1),
                           LTBI2=median(LTBI2), uq=uq(LTBI),lq=lq(LTBI)),
                     by=list(acat,g_whoregion)]

rundata0[,df:=LTBI2-LTBI1]

rundata0 <- melt(rundata0,id=c('acat','g_whoregion','uq','lq'))


acts <- paste0(5*(0:16),'-',5*(1:17))
acts[17] <- c('80-')
rundata0$age <- acts
rundata0$age <- factor(rundata0$age,levels=acts,ordered = TRUE)
## now examine prevalence by age


## look at country varations
rundatac <- rundata[,list(LTBI=sum(pop*P)/sum(pop),
                          LTBI1=sum(pop*P1)/sum(pop),
                          LTBI2=sum(pop*P2)/sum(pop),
                          LTBIH=sum(pop*P2*INH)/sum(pop)),
                    by=list(iso3,replicate)]

## for mapping
rundatac <- rundatac[,list(LTBI=median(LTBI),
                          LTBI1=median(LTBI1),
                          LTBI2=median(LTBI2),
                          LTBIH=median(LTBIH)),
                    by=list(iso3)]


## look at country varations
rundatacT <- rundata[,list(LTBI=sum(pop*P),
                          LTBI1=sum(pop*P1),
                          LTBI2=sum(pop*P2),
                          LTBIH=sum(pop*P2*INH)),
                    by=list(iso3,replicate)]

## for mapping
rundatacT <- rundatacT[,list(LTBI=median(LTBI),
                             LTBI_lo=quantile(LTBI,probs=0.025),
                             LTBI_hi=quantile(LTBI,probs=0.975),
                          LTBI1=median(LTBI1),
                          LTBI2=median(LTBI2),
                          LTBIH=median(LTBIH)),
                    by=list(iso3)]


top20c <- head(rundatacT[order(LTBI,decreasing = TRUE),list(iso3,LTBI,LTBI_lo,LTBI_hi)],n=20)
top20c$iso3 <- factor(top20c$iso3,levels=rev(top20c$iso3),ordered = TRUE)

top20c2 <- merge(top20c,rundatac[,list(iso3,LTBI)],by='iso3',all.x=TRUE,all.y=FALSE)
top20c2 <- as.data.frame(top20c2)
names(top20c2)[c(2,5)] <- c('LTBI','Percent')
top20c2$Percent <- 1e2*top20c2$Percent

## ---------- numbers ------------
## -- by country --
rundataic <- rundata[,list(nLTBI=sum(pop*P),
                          nLTBI1=sum(pop*P1),
                          nLTBI2=sum(pop*P2),
                          nLTBIH=sum(pop*P2*INH)),
                    by=list(replicate,iso3)]

MID <- rundataic[,list(nLTBI=median(nLTBI),
                          nLTBI1=median(nLTBI1),
                          nLTBI2=median(nLTBI2),
                          nLTBIH=median(nLTBIH)),
                    by=list(iso3)]
LO <- rundataic[,list(nLTBI=lq(nLTBI),
                          nLTBI1=lq(nLTBI1),
                          nLTBI2=lq(nLTBI2),
                          nLTBIH=lq(nLTBIH)),
                    by=list(iso3)]
HI <- rundataic[,list(nLTBI=uq(nLTBI),
                          nLTBI1=uq(nLTBI1),
                          nLTBI2=uq(nLTBI2),
                          nLTBIH=uq(nLTBIH)),
               by=list(iso3)]


## kids
rundataick <- rundata[acat<3,list(nLTBI=sum(pop*P),
                          nLTBI1=sum(pop*P1),
                          nLTBI2=sum(pop*P2),
                          nLTBIH=sum(pop*P2*INH)),
                    by=list(replicate,iso3)]

MIDk <- rundataick[,list(nLTBI=median(nLTBI),
                          nLTBI1=median(nLTBI1),
                          nLTBI2=median(nLTBI2),
                          nLTBIH=median(nLTBIH)),
                    by=list(iso3)]
LOk <- rundataick[,list(nLTBI=lq(nLTBI),
                          nLTBI1=lq(nLTBI1),
                          nLTBI2=lq(nLTBI2),
                          nLTBIH=lq(nLTBIH)),
                    by=list(iso3)]
HIk <- rundataick[,list(nLTBI=uq(nLTBI),
                          nLTBI1=uq(nLTBI1),
                          nLTBI2=uq(nLTBI2),
                          nLTBIH=uq(nLTBIH)),
               by=list(iso3)]


MID <- as.data.frame(MID);LO <- as.data.frame(LO);HI <- as.data.frame(HI);
MID$K <- MIDk$nLTBI; LO$K <- LOk$nLTBI; HI$K <- HIk$nLTBI



## -- by region --
## total numbs
rundata1 <- rundata[,list(nLTBI=sum(pop*P),
                          nLTBI1=sum(pop*P1),
                          nLTBI2=sum(pop*P2),
                          nLTBIH=sum(pop*P2*INH)),
                    by=list(replicate,g_whoregion)]
rundata1g <- rundata[,list(nLTBI=sum(pop*P),
                          nLTBI1=sum(pop*P1),
                          nLTBI2=sum(pop*P2),
                          nLTBIH=sum(pop*P2*INH)),
                    by=list(replicate)]
MID <- rundata1[,list(nLTBI=median(nLTBI),
                          nLTBI1=median(nLTBI1),
                          nLTBI2=median(nLTBI2),
                          nLTBIH=median(nLTBIH)),
                    by=list(g_whoregion)]
LO <- rundata1[,list(nLTBI=lb(nLTBI),
                          nLTBI1=lb(nLTBI1),
                          nLTBI2=lb(nLTBI2),
                          nLTBIH=lb(nLTBIH)),
                    by=list(g_whoregion)]
HI <- rundata1[,list(nLTBI=ub(nLTBI),
                          nLTBI1=ub(nLTBI1),
                          nLTBI2=ub(nLTBI2),
                          nLTBIH=ub(nLTBIH)),
               by=list(g_whoregion)]

## kids
rundata1k <- rundata[acat<3,list(nLTBI=sum(pop*P),
                          nLTBI1=sum(pop*P1),
                          nLTBI2=sum(pop*P2),
                          nLTBIH=sum(pop*P2*INH)),
                    by=list(replicate,g_whoregion)]
rundata1gk <- rundata[acat<3,list(nLTBI=sum(pop*P),
                          nLTBI1=sum(pop*P1),
                          nLTBI2=sum(pop*P2),
                          nLTBIH=sum(pop*P2*INH)),
                    by=list(replicate)]
MIDk <- rundata1k[,list(nLTBI=median(nLTBI),
                          nLTBI1=median(nLTBI1),
                          nLTBI2=median(nLTBI2),
                          nLTBIH=median(nLTBIH)),
                    by=list(g_whoregion)]
LOk <- rundata1k[,list(nLTBI=lb(nLTBI),
                          nLTBI1=lb(nLTBI1),
                          nLTBI2=lb(nLTBI2),
                          nLTBIH=lb(nLTBIH)),
                    by=list(g_whoregion)]
HIk <- rundata1k[,list(nLTBI=ub(nLTBI),
                          nLTBI1=ub(nLTBI1),
                          nLTBI2=ub(nLTBI2),
                          nLTBIH=ub(nLTBIH)),
                    by=list(g_whoregion)]

## all
MID <- as.data.frame(MID);LO <- as.data.frame(LO);HI <- as.data.frame(HI);
MIDg <- rundata1g[,list(nLTBI=median(nLTBI),
                          nLTBI1=median(nLTBI1),
                          nLTBI2=median(nLTBI2),
                          nLTBIH=median(nLTBIH))]
LOg <- rundata1g[,list(nLTBI=lb(nLTBI),
                          nLTBI1=lb(nLTBI1),
                          nLTBI2=lb(nLTBI2),
                          nLTBIH=lb(nLTBIH))]
HIg <- rundata1g[,list(nLTBI=ub(nLTBI),
                          nLTBI1=ub(nLTBI1),
                          nLTBI2=ub(nLTBI2),
                          nLTBIH=ub(nLTBIH))]
MIDg <- as.data.frame(MIDg);LOg <- as.data.frame(LOg);HIg <- as.data.frame(HIg);
MIDg <- cbind(g_whoregion='GLOBAL',MIDg);LOg <- cbind(g_whoregion='GLOBAL',LOg);
HIg <- cbind(g_whoregion='GLOBAL',HIg)
MID <- rbind(MID,MIDg);HI <- rbind(HI,HIg);LO <- rbind(LO,LOg)

## kids
MID <- cbind(MID, c(MIDk[,nLTBI],rundata1gk[,median(nLTBI)] ))
LO <- cbind(LO, c(LOk[,nLTBI],rundata1gk[,lb(nLTBI)] ))
HI <- cbind(HI, c(HIk[,nLTBI],rundata1gk[,ub(nLTBI)] ))
names(MID)[6] <- names(LO)[6] <- names(HI)[6] <- 'nkLTBI'


## ---  for INH proportions ---
rundatR <- rundata1[,list(propRH=nLTBIH/nLTBI2),by=list(replicate,g_whoregion)]
rundatRg <- rundata1g[,list(propRH=nLTBIH/nLTBI2),by=list(replicate)]

MIDR <- rundatR[,list(mid=median(propRH)),by=g_whoregion]
LOR <- rundatR[,list(lo=lb(propRH)),by=g_whoregion]
HIR <- rundatR[,list(hi=ub(propRH)),by=g_whoregion]

MIDR <- as.data.frame(MIDR);LOR <- as.data.frame(LOR);HIR <- as.data.frame(HIR);
MIDR <- rbind(MIDR,data.frame(g_whoregion='GLOBAL',mid=median(rundatRg[,propRH])))
LOR <- rbind(LOR,data.frame(g_whoregion='GLOBAL',lo=lb(rundatRg[,propRH])))
HIR <- rbind(HIR,data.frame(g_whoregion='GLOBAL',hi=ub(rundatRg[,propRH])))
 

## ---------- props ------------
## total props
rundata2 <- rundata[,list(LTBI=sum(pop*P)/sum(pop),
                          LTBI1=sum(pop*P1)/sum(pop),
                          LTBI2=sum(pop*P2)/sum(pop),
                          LTBIH=sum(pop*P2*INH)/sum(pop)),
                    by=list(replicate,g_whoregion)]
rundata2g <- rundata[,list(LTBI=sum(pop*P)/sum(pop),
                          LTBI1=sum(pop*P1)/sum(pop),
                          LTBI2=sum(pop*P2)/sum(pop),
                          LTBIH=sum(pop*P2*INH)/sum(pop)),
                    by=list(replicate)]


## prop in kids
rundata2[,propk:=rundata1k[,nLTBI]/rundata1[,nLTBI]]
rundata2g[,propk:=rundata1gk[,nLTBI]/rundata1g[,nLTBI]]

MID <- rundata2[,list(LTBI=median(LTBI),
                      LTBI1=median(LTBI1),
                      LTBI2=median(LTBI2),
                      LTBIH=median(LTBIH),
                      propk=median(propk)),
                    by=list(g_whoregion)]
LO <- rundata2[,list(LTBI=lb(LTBI),
                     LTBI1=lb(LTBI1),
                     LTBI2=lb(LTBI2),
                     LTBIH=lb(LTBIH),
                     propk=lb(propk)),
                    by=list(g_whoregion)]
HI <- rundata2[,list(LTBI=ub(LTBI),
                     LTBI1=ub(LTBI1),
                     LTBI2=ub(LTBI2),
                     LTBIH=ub(LTBIH),
                     propk=ub(propk)),
                    by=list(g_whoregion)]

MID <- as.data.frame(MID);LO <- as.data.frame(LO);HI <- as.data.frame(HI);
MIDg <- rundata2g[,list(LTBI=median(LTBI),
                        LTBI1=median(LTBI1),
                        LTBI2=median(LTBI2),
                        LTBIH=median(LTBIH),
                        propk=median(propk))]
LOg <- rundata2g[,list(LTBI=lb(LTBI),
                       LTBI1=lb(LTBI1),
                       LTBI2=lb(LTBI2),
                       LTBIH=lb(LTBIH),
                       propk=lb(propk))]
HIg <- rundata2g[,list(LTBI=ub(LTBI),
                       LTBI1=ub(LTBI1),
                       LTBI2=ub(LTBI2),
                       LTBIH=ub(LTBIH),
                       propk=ub(propk))]
MIDg <- as.data.frame(MIDg);LOg <- as.data.frame(LOg);HIg <- as.data.frame(HIg);
MIDg <- cbind(g_whoregion='GLOBAL',MIDg);LOg <- cbind(g_whoregion='GLOBAL',LOg);
HIg <- cbind(g_whoregion='GLOBAL',HIg)
MID <- rbind(MID,MIDg);HI <- rbind(HI,HIg);LO <- rbind(LO,LOg)

## per mille
LO$LTBIH <- LO$LTBIH*10; MID$LTBIH <- MID$LTBIH*10; HI$LTBIH <- HI$LTBIH*10


## -------- 2035 & 2050 -----------


N2035 <- POP2035[,sum(value)]
POP2035$acatn <- 0:16
POP2035$acat <- POP2035$acatn-4
POP2035 <- POP2035[acat>=0,]
POP2035 <- POP2035[,list(iso3,acat,pop35=value)]
N2050 <- POP2050[,sum(value)]
POP2050$acatn <- 0:16
POP2050$acat <- POP2050$acatn-7
POP2050 <- POP2050[acat>=0,]
POP2050 <- POP2050[,list(iso3,acat,pop50=value)]

fac <- 1e5*(0.15*1e-2)                  #reactivation rate

## 2035
fruns <- merge(rundata[,list(iso3,replicate,acat,P,g_whoregion)],
               POP2035,by=c('iso3','acat'),all=TRUE)
fruns <- fruns[!is.na(pop35),]       #ditch the dead
tmp <- fruns[,list(nLTBI=sum(pop35*P)),by=replicate]
tmp <- tmp[!is.na(replicate)]
tmp[,list(mid=1e3*median(nLTBI),lo=1e3*lb(nLTBI),hi=1e3*ub(nLTBI))]


fac*tmp[,list(mid=median(nLTBI)/N2035,lo=lb(nLTBI)/N2035,hi=ub(nLTBI)/N2035)]

## 2050
fruns <- merge(rundata[,list(iso3,replicate,acat,P,g_whoregion)],
               POP2050,by=c('iso3','acat'),all=TRUE)
fruns <- fruns[!is.na(pop50),]       #ditch the dead
tmp <- fruns[,list(nLTBI=sum(pop50*P)),by=replicate]
tmp <- tmp[!is.na(replicate)]
tmp[,list(mid=1e3*median(nLTBI),lo=1e3*lb(nLTBI),hi=1e3*ub(nLTBI))]

fac*tmp[,list(mid=median(nLTBI)/N2050,lo=lb(nLTBI)/N2050,hi=ub(nLTBI)/N2050)]

## 1997 --- NB rundata not available after this

rund[,year:=1997-year]               #now age
rund[,lari:=exp(lari)]               #now real ari
rund <- rund[year>=0,]               #ditch future
rund <- rund[iso3 %in% unique(rundata$iso3)] #droppings

## past
past <- data.frame(year=64:80,iso3=rep(unique(rund$iso3),each=length(64:80)),
                   lari=NA, replicate =
                       rep(1:200,each=length(64:80)*length(unique(rund$iso3)) ))
past <- as.data.table(past)
pastlari <- rund[year==63,list(lari),by=list(iso3,replicate)]
past <- merge(past,pastlari,all.x=TRUE,by=c('iso3','replicate'))
past[,lari:=lari.y]
past <- past[,list(iso3,replicate,year,lari)]

rund <- merge(rund,past,all=TRUE,by=c('iso3','year','replicate'))
rund[,lari:=pmax(lari.x,lari.y,na.rm=TRUE)]
rund <- rund[,list(iso3,replicate,year,lari)]
rund <- rund[order(replicate,iso3,year),]
rund <- rund[,list(ari=lari,H=cumsum(lari),year=year),by=list(iso3=iso3,replicate=replicate)]
rund[,acat:=year %/% 5]
rund[,P:=1-exp(-H)]                  #ever
rund <- rund[,list(P=mean(P)),
                   by=list(iso3,replicate,acat)]
rund <- merge(rund,WHOkey[,c('iso3','g_whoregion')],by='iso3',all.x=TRUE)

## demographic data

N1997 <- POP1997[,sum(value)]
POP1997$acat <- 0:16
POP1997 <- POP1997[,list(iso3,acat,pop97=value)]
## 93% younger than 65

pruns <- merge(rund[,list(iso3,replicate,acat,P,g_whoregion)],
               POP1997,by=c('iso3','acat'),all.x=TRUE)
## pruns[,sum(pop97)/200]*1e-6
pruns <- pruns[!is.na(P),]

tmp <- pruns[,list(nLTBI=sum(pop97*P)),by=replicate]
tmp[,list(mid=1e3*median(nLTBI),lo=1e3*lb(nLTBI),hi=1e3*ub(nLTBI))]*1e-6

tmp <- pruns[,list(pLTBI=sum(pop97*P)/sum(pop97)),by=replicate]
tmp[,list(mid=median(pLTBI),lo=lb(pLTBI),hi=ub(pLTBI))]
