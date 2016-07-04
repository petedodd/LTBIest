## Analysis code for Houben & Dodd 2016, distributed under CC BY 4.0 license https://creativecommons.org/licenses/by/4.0/
library(data.table)
Nro <- read.csv('data/inputData_GTB_estimate.csv')
Nro <- as.data.table(Nro)

## special countries -- mergers
## use SCG for MNE and for SRB
tmp <- copy(Nro[iso3=='SCG'])                 #serbia and mne 1990:2004
tmp$iso3 <- 'MNE'                       #use SCG for MNE
Nro <- rbind(Nro,tmp)                   
tmp$iso3 <- 'SRB'                       #use SCG for SRB
Nro <- rbind(Nro,tmp)
Nro <- Nro[iso3!='SCG',]                #drop SCG
## SSD and SDN
tmp <- copy(Nro[iso3=='SDN'])                 #sudan
tmp <- tmp[year<=2010]                        #before split
tmp$iso3 <- 'SSD'                       #use SDN for SSD before 2011
Nro <- rbind(Nro,tmp)                   
Nro <- Nro[order(iso3,year)]            #order
Nro$iso3 <- factor(Nro$iso3)            #refactor

## Cauthen
C <- read.csv('data/inputData_Cauthen_ARI.csv')
C$ari <- 1e-2*C$ARI
## get ages
st <- unlist(lapply(X=as.character(C$AGE),FUN=function(x)regexpr('\\(',x)[[1]]))
sp <- unlist(lapply(X=as.character(C$AGE),FUN=function(x)regexpr('\\)',x)[[1]]))
age <- substr(as.character(C$AGE),start=st+1,stop=sp-1)
C$age <- as.numeric(age)
C$var <- (C$age*C$ari)/C$N
## log(x+dx) = log(x) + dx/x - .5*(dx/x)^2
## Var[log(X)] ~ Var[X]/E[X]^2  = E^2
CE <- C[,c('iso3','year','ari')]
CE$E <- sqrt(C$var)/C$ari
CE$lari <- log(CE$ari)
CE$type <- 'TST survey'
CE <- na.omit(CE)

## new data ---
C2 <- read.csv('data/inputData_Review_ARI.csv')
C2$ari <- 1e-2*C2$ARI
st <- unlist(lapply(X=as.character(C2$AGE),FUN=function(x)regexpr('\\(',x)[[1]]))
sp <- unlist(lapply(X=as.character(C2$AGE),FUN=function(x)regexpr('\\)',x)[[1]]))
age <- substr(as.character(C2$AGE),start=st+1,stop=sp-1)
C2$age <- as.numeric(age)
C2$var <- (C2$age*C2$ari)/C2$N
C2$var[!is.na(C2$ARI_hi)] <- (1.96e-2*(C2$ARI_hi-C2$ARI_lo)[!is.na(C2$ARI_hi)])^2
CE2 <- C2[,c('iso3','year','ari')]
CE2$E <- sqrt(C2$var)/C2$ari
CE2$lari <- log(CE2$ari)
CE2$type <- 'TST survey'
CE2 <- na.omit(CE2)

## HIV data proportion HIV in TB
NH <- Nro[,list(iso3,year,
              cdr = c_cdr/100,vcdr = (c_cdr_hi-c_cdr_lo)^2/(100*3.92)^2,
              p = e_tbhiv_prct/100,
              vp = (e_tbhiv_prct_hi-e_tbhiv_prct_lo)^2/(100*3.92)^2)]
## drop NA from % HIV
NH <- NH[!is.na(p),]
## drop <1% for o
NH <- NH[p>1e-2,]                       #2931
## NH[is.na(cdr),]                         #95 country-years, many during conflict
bcnz <- unique(NH[is.na(cdr),iso3])     #38 countries
for(cn in as.character(bcnz)){          #min CDR, max var
    NH$cdr[NH$iso3==cn & is.na(NH$cdr)] <- min(NH$cdr[NH$iso3==cn & !is.na(NH$cdr)])
    NH$vcdr[NH$iso3==cn & is.na(NH$vcdr)] <- max(NH$vcdr[NH$iso3==cn & !is.na(NH$vcdr)])
}
NH[is.na(cdr),]                         #check: none
NH[cdr>1]$vcdr <- 0                     #cap CDR for this
NH[cdr>1]$cdr <- 1
NH[is.na(vp)]$vp <- 0                   #these are where p=1
NH2 <- copy(NH)                         #for sense below
## T1 = CDR.T1n + (1-CDR).T1u
## T2 = CDR.T2n + (1-CDR).T2u
## var(T) ~= (Tn-Tu)^2.var(CDR)  + CDR^2.var(Tn) + (1-CDR)^2.var(Tu)
## T1, T2 and variances...
NH[,T1:=cdr*(2+.2)/2+(1-cdr)*(4+1)/2]
NH[,T2:=cdr*(1+.01)/2+(1-cdr)*(.2+.01)/2]
NH[,vT1:=vcdr*.25*(2+.2-4-1)^2 + cdr^2*(2-.2)^2/12 + (1-cdr)^2*(4-1)^2/12]
NH[,vT2:=vcdr*.25*(1+.01-.2-.01)^2 + cdr^2*(1-.01)^2/12 + (1-cdr)^2*(.2-.01)^2/12]
## get f...35/45 with var from U[.3,.4]/U[.4,.5]
mf <- 0.780
vf <- 0.00666
## smear adjustment factor: p is hivntb, f redn if HIV+, T2 hiv+ durn, T1 hiv- durn
## S = (p.T2.f + (1-p).T1) / (p.T2 + (1-p).T1) = A / B
## (up to sign)
## dlogS = (dp).p(1-p)(1-f).T1 / (AB) + (dT2).p(1-p)(1-f).T1 / (AB) 
##  + (dT1).p(1-p)(1-f).T2 / (AB) + (df).p.T2 / A
## variance
## var(log S) = var(p).[p(1-p)(1-f).T1 / (AB)]^2
## + [T1^2.var(T2) + T2^2.var(T1)].[p(1-p)(1-f) / (AB) ]^2  + var(f).[p.T2 / A]^2 
## var log S....
NH[,A:=(p*T2*mf + (1-p)*T1)]
NH[,B:=(p*T2 + (1-p)*T1)]
NH[,vlS:= vp*(p*(1-p)*(1-mf)*T1/(A*B))^2 +
       (T1^2*vT2+T2^2*vT1)*(p*(1-p)*(1-mf)/(A*B))^2 +
           vf*(p*T2/A)^2]               # varlogS/logS^2=E^2
NH[,S:= A/B]                            #S - fraction redn
NH3 <- copy(NH)                         #for sense below
NH <- NH[,list(iso3,year,vlS,S)]

## kids-------
## Kunkel review...
## 0.52 0.4-0.64 adults
## .5% 0-1.9% u5
## 14% 8.9-19.4% 5-14
## proportion of incidence in children
K <- read.csv('data/inputData_GTB_kidestimate.csv')
K <- merge(Nro[year==2014,list(iso3,year,e_inc_num,e_inc_num_lo,e_inc_num_hi)],K,by='iso3')
K[,pk:=inc.num/e_inc_num]               #proportion in kids
K[,vpk:=pk^2*(inc.num.sd^2/inc.num^2 + ((e_inc_num_hi-e_inc_num_lo)/3.92)^2/e_inc_num^2)] #variance in prop in kids
K[pk>1,]$pk <- .1                        #correct a few weirdos
K[vpk>1e-2,]$vpk <- 1e-2                   #correct a few weirdos
## age split
load('data/kidsplit.Rdata')             #Dodd,Sismanidis,Seddon
AF[,pu5:=a/(a+b)]
AF[,vpu5:=(1-pu5)*pu5/(a+b+1)]
K <- merge(K[,list(iso3,pk,vpk)],AF[,list(iso3,pu5,vpu5)],by='iso3')
## fraction by age
YK <- .5e-2
OK <- 14e-2
A <- 0.52
vYK <- (1e-2*1.9/3.92)^2
vOK <- ((19.4-8.9)*1e-2/3.92)^2
vA <- ((0.64-0.4)/3.92)^2
K[,mFr:= pk*pu5*YK + pk*(1-pu5)*OK + (1-pk)*A] #mean fraction smr+
K[,vFr:= mFr^2 * ( vpk*(pu5*YK+(1-pu5)*OK-A)^2 + vpu5*(pk*YK-pk*OK)^2 +
                vYK*(pk*pu5)^2 + vOK*(pk*(1-pu5))^2 + vA*(1-pk)^2)] #var frac smr+
KH <- merge(NH[,list(iso3,year,S,vlS)],K[,list(iso3,mFr,vFr)],by='iso3',all=TRUE)
KH[is.na(year),]$year <- 2014


## now prevalence in same fashion
## log-normal: mu=1.678 sig=0.371
mu <- 1.678
sig <- .371
F <- exp(mu + .5*sig^2)                                  #styblo factor!! (lnorm mn)
vF <- (exp(sig^2)-1)*exp(2*mu+sig^2)                     #log normal variance

## create data
Nr <- copy(Nro)
Nr <- merge(Nr,KH,by=c('iso3','year'),all = TRUE) #merge
## dealing with missing mFr from this merge
bad <- Nr[,all(is.na(mFr)),by=iso3]
for(cn in bad[V1==FALSE,as.character(iso3)]){ #countries with at least one mFr
    if(Nr[iso3==cn,any(is.na(mFr))]){         #a problem like AFG
        vlu <- mean(Nr[iso3==cn,mFr],na.rm=TRUE) #should be just the value
        Nr[iso3==cn,]$mFr <- vlu
        vlu <- mean(Nr[iso3==cn,vFr],na.rm=TRUE) #should be just the value
        Nr[iso3==cn,]$vFr <- vlu
    }
}
## dealing with countries with no mFr data
mmFr <- Nr[,median(mFr,na.rm=TRUE)]
mvFr <- Nr[,median(vFr,na.rm=TRUE)]
for(cn in bad[V1==TRUE,as.character(iso3)]){ #countries with no mFr data
    Nr[iso3==cn,]$mFr <- mmFr
    Nr[iso3==cn,]$vFr <- mvFr
}
## dealing with missing HIV
Nr[is.na(S),]$vlS <- 0                            #countries without HIV
Nr[is.na(S),]$S <- 1                    #countries without HIV

## calculations
Nrc <- copy(Nr)
Nrc[,ari:=e_prev_100k*1e-5 * mFr * S * F]
Nrc[,lari:=log(ari)]
Nrc[,E1:=(e_prev_100k_hi-e_prev_100k_lo)/(3.92*e_prev_100k)] #sd/mean from prev
Nrc[,E:=sqrt(E1^2+vF/F^2 + vlS + vFr/mFr^2)]
Nrc <- as.data.frame(Nrc[,list(iso3,year,ari,E,lari)])
Nrc$type <- 'Prevalence estimate'

All <- rbind(CE,Nrc,CE2)
## main uncertainty is from prevalence, then styblo relation

save(All,file='data/All_3.Rdata')

## sensitivity to HIV+ CDR
NH2[,cdr2:=1.0]
NH2[,T1:=cdr*(2+.2)/2+(1-cdr)*(4+1)/2]
NH2[,T2:=cdr2*(1+.01)/2+(1-cdr2)*(.2+.01)/2]
NH2[,vT1:=vcdr*.25*(2+.2-4-1)^2 + cdr^2*(2-.2)^2/12 + (1-cdr)^2*(4-1)^2/12]
NH2[,vT2:=vcdr*.25*(1+.01-.2-.01)^2 + cdr2^2*(1-.01)^2/12 + (1-cdr2)^2*(.2-.01)^2/12]
NH2[,A:=(p*T2*mf + (1-p)*T1)]
NH2[,B:=(p*T2 + (1-p)*T1)]
NH2[,vlS:= vp*(p*(1-p)*(1-mf)*T1/(A*B))^2 +
       (T1^2*vT2+T2^2*vT1)*(p*(1-p)*(1-mf)/(A*B))^2 +
           vf*(p*T2/A)^2]               # varlogS/logS^2=E^2
NH2[,S:= A/B]                            #S - fraction redn
summary(1e2*(1-NH2$S/NH3$S))             #0.14% IQR=[0.04 - 0.55]

