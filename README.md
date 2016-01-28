# kalmanST

# Warning: package under development!!

rm(list=ls())

library(devtools)
library(kalmanST)
library(spatsurv)
library(sp)
load_all("kalmanST")

source("optifix.R")

set.seed(1)

N <- 50
T <- 100

ids <- rep(1:N,T)
uqids <- 1:N
times <- rep(1:T,each=N)
xcoords <- runif(N)
ycoords <- runif(N)

COVMODEL <- ExponentialCovFct()

u <- as.vector(as.matrix(dist(cbind(xcoords,ycoords))))

a <- 0.8
sigma_S <- 5
phi_S <- 0.1
sigma_Y <- 3

Sigma <- matrix(EvalCov(COVMODEL,u=u,parameters=c(sigma_S,phi_S)),N,N)
cholS <- t(chol(Sigma))

S <- matrix(cholS%*%rnorm(N),N,1)
for(i in 2:T){
    S <- cbind(S,a*S[,i-1]+sqrt(1-a^2)*cholS%*%rnorm(N))
}
S <- as.vector(S)

age <- runif(N*T,20,70)
sex <- sample(c("Male","Female"),N*T,replace=TRUE)

outcome <- age*0.5 + 10
outcome[sex=="Female"] <- outcome[sex=="Female"] - 10
outcome <- outcome + S + rnorm(N*T,0,sigma_Y)

dat <- data.frame(ID=ids,times=times,xcoord=rep(xcoords,T),ycoord=rep(ycoords,T),age=age,sex=sex,outcome=outcome)
coordinates(dat) =~ xcoord + ycoord

#dat <- dat[-sample(1:nrow(dat),floor(N*T/20))]
#dat$outcome[(N*T-20*N):(N*T)] <- NA

out <- kffit_stationary_STGP(   form=outcome ~ age + sex,
                                param=paramSpec(a=a,sigma_S=sigma_S,phi_S=phi_S,sigma_Y=sigma_Y),
                                pid="ID",
                                tid="times",
                                tmax=T,
                                data=dat,
                                cov.model=COVMODEL)


out <- kffit_stationary_STGP(   form="outcome",
                                param=paramSpec(a=a,sigma_S=sigma_S,phi_S=phi_S,sigma_Y=sigma_Y),
                                pid="ID",
                                tid="times",
                                tmax=T,
                                data=dat,
                                cov.model=COVMODEL)


#predictPlots(out)

fun <- kffit_stationary_STGP(   form=outcome ~ age + sex,
                                param=paramSpec(a=a,sigma_S=sigma_S,phi_S=phi_S,sigma_Y=sigma_Y),
                                pid="ID",
                                tid="times",
                                tmax=T,
                                data=dat,
                                cov.model=COVMODEL,
                                returnFunction=TRUE)

pars <- KFparest(fun)
cat("\n\n")
print(pars)
