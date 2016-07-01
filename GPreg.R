## Analysis code for Houben & Dodd 2016, distributed under CC BY 4.0 license https://creativecommons.org/licenses/by/4.0/
args <- commandArgs(TRUE)
jj <- as.numeric(args[1])
print(jj)
library(MASS)
library(Matrix)
load('data/All_3.Rdata')
All <- All[All$lari!=-Inf,]
lin <- 0                              #linear or constant  --- CHANGE HERE!
cn <- unique(as.character(All$iso3))[jj]
All <- All[All$iso3 %in% cn,]

getKKtonly <- function(t1,t2,k=function(x,y) exp(-abs(x-y)),Wm){
    K <- outer(t1,t2,FUN=k)
    K
}                                       #make K(X,Y) matrices

## function to do double-backsolve inversion given cholesky decomp
bsi <- function(CH,v) backsolve(CH,backsolve(CH, v, transpose = TRUE))

## this is basically 1,t x I(country==j)
getHtonly <- function(t,n=1){             #n is highest power
    H <- Matrix(0,nrow=(n+1),ncol=length(t)) #3 for 1,t,t^2
    for(k in 0:n)
        H[1 + k,] <- t^k
    as(H,'sparseMatrix')
}

tr <- function(x) sum(diag(x))

## gets s_k^2, L, tscale
PRMconvert2 <- function(x) c( exp(x[1]), exp(x[2]/2))

getMLnGradT <- function(x,grad=TRUE){                    # see Eqn 2.45 in R&W
    ## preliminaries
    ## -- covariances --
    a <- x[1]; b <- x[2]
    K <- outer(tdz,tdz,FUN=function(x,y) exp(a-exp(-b)*(x-y)^2) )
    K <- Matrix(K)                      #kxx
    K2 <- outer(tdz,tdz,FUN=function(x,y) (x-y)^2*exp(a-exp(-b)*(x-y)^2) )
    K2 <- Matrix(K2)                     #dK/da
    ## -- derived matrices (as used above) --
    ## new version
    cvy <- K + sigz
    U <- chol(cvy)
    Uy <- backsolve(U, y, transpose = TRUE)
    ky <- backsolve(U, Uy) 
    hky <- H %*% ky
    AM <- symmpart(H %*% bsi(U,t(H)))
    V <- chol(AM)
    Vy <- backsolve(V, hky, transpose = TRUE)
    ## -- marginal loglikelihood --
    LML <- -sum(Uy^2) + sum(Vy^2) #data likelihood
    LML <- LML - 2*sum(log(diag(U)))
    LML <- LML - 2*sum(log(diag(V)))
    LML <- LML/2
    ## extras for dLML
    VVy <- backsolve(V, Vy)
    ympy <- bsi(U, t(H) %*% VVy)       #K^{-1} x ...
    dLML <- NULL
    if(grad){
        ## -- gradient --
        ## --- gradient helper function ---
        dHelp <- function(dK){              #takes the local vars from up-level
            dML <- t(ky) %*% dK %*% ky + t(ympy) %*% dK %*% ympy
            dML <- dML - 2*t(ympy) %*% dK %*% ky
            tmp <- bsi(U,dK)        #K^{-1}dK
            dML <- dML - tr(tmp)
            tmp <- bsi(V,H%*%tmp)
            tmp <- bsi(U, t(H)%*%tmp)
            dML <- dML + tr(tmp)  
            return(as.numeric(dML/2))
        }
        ## --- get gradient ---
        dLML <- c( dHelp(K), dHelp(K2) )
    }
    if(LML>0){LML <- -1e3; dLML <- -1e3*rep(1,2)}
    ## return
    return(list(LML=LML,dLML=dLML))
}

getPredztonly <- function(x,tdz,tez,y,Vz){
    ## from here
    usek <- function(i,j)x[1]*exp(-abs(i-j)^2/x[2]^2)
    ## H
    H <- getHtonly(tdz,n=lin)
    Hs <- getHtonly(tez,n=lin)
    ## make matrices
    kxx <- getKKtonly(tdz,tdz,k=usek)
    kxxs <- getKKtonly(tdz,tez,k=usek)
    kxsx <- t(kxxs)
    kxsxs <- getKKtonly(tez,tez,k=usek)
    sigz <- Diagonal(x=Vz)                         #noise
    covy <- kxx + sigz
    U <- chol(covy)
    ## regress
    reg <- 1                                #flag
    HKH <- H %*% bsi(U,t(H))
    V <- chol(symmpart(HKH))                       
    R <- Hs - H %*% bsi(U,kxxs)
    mn <- 0
    y <- y-mn
    ## mean/covar
    mf <- kxsx %*% bsi(U,y)          #mean prediction
    cf <- kxsxs  - kxsx %*% bsi(U,kxxs)
    ## mean stuff
    bbar <- bsi(V,(H %*% bsi(U,y)))
    mg <- mf + reg * t(R) %*% bbar
    if(nrow(V)>1)
        cg <- cf + reg * t(R) %*% bsi(V,R)
    else
        cg <- cf + reg * t(R) %*% (R/V[1,1])
    ## return
    return(list(mg=mg,cg=cg))
}


## ============== work ====================

## time/country vectors for data/extrapolation
fyear <- 1933
tdz <- All$year-fyear
tez <- 1934:2014 - fyear                  #extrapolation times

## H as a global
## dcnz, tdz, sigz, y as globals
## H
H <- getHtonly(tdz,n=lin)
Hs <- getHtonly(tez,n=lin)
y <- All$lari
vz <- (All$E)^2
sigz <- Diagonal(x=vz)                         #noise

mz <- c(log(.5),2*1.5*log(2));sz <- c(1,1)*100
LMLfun2 <- function(x) -(getMLnGradT(x,grad=FALSE)$LML
                         - sum(.5*(x-mz)^2/sz^2))
dLMLfun2 <- function(x) -(getMLnGradT(x)$dLML-(x-mz)/sz^2)
x02 <- mz
## optimize
system.time({                           #<1s
    testo2 <- optim(par=x02,fn = LMLfun2,gr = dLMLfun2)
})
pab <- testo2$par
ab <- PRMconvert2(pab)
## ab
print(ab);
## NB trade-off flexibility vs growing noise going back
xx <- ab
tot <- getPredztonly(xx,tdz,tez,y,vz)
scf <- as.numeric(sqrt(diag(tot$cg)))
erw <- data.frame(year=tez+fyear,iso3=as.character(unique(All$iso3)),
                  lari=as.numeric(tot$mg),
                  upper=as.numeric(tot$mg) + 1.96*scf,
                  lower=as.numeric(tot$mg) - 1.96*scf)

save(erw,file=paste0('data/',cn,'.Rdata'))

runs <- mvrnorm(n=2e2,mu=as.numeric(tot$mg),Sigma=as.matrix(symmpart(tot$cg)))
runsdf <- data.frame(year=tez+fyear,iso3=as.character(unique(All$iso3)),
                  lari=c(t(runs)),replicate=rep(1:nrow(runs),each=ncol(runs)))

save(runsdf,file=paste0('data/zz_',cn,'.Rdata'))
