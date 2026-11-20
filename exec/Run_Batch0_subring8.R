#Script for analyzing age-readings from annual Herring spawning survey in February/March
#Only use data from 2022
library(AgeReading)
library(rstan)
library(TruncatedNormal)
logit = function(p){log(p/(1-p))}
invlogit = function(x){1/(1+exp(-x))}

## Choose dataset  (HerringSurvey/HerringSurvey202?)
data(Ring.readErr)
dataRing = Ring.readErr
data(HerringSurvey)
df = HerringSurvey

## Subsample the ring data
R = dim(dataRing)[1]
A = dim(dataRing)[2]
Mind = matrix(NA,nrow=A,ncol=A)
for(a in 1:A)
  for(b in 1:A)
    Mind[a,b]=(a-1)*A+b
dataRing2 = 0*dataRing
for(r in 1:R)
{
  n = sum(dataRing[r,,])
  n2 = as.integer(n/8)
  v = rep(Mind,dataRing[r,,])
  v2 = sample(v,n2,replace=FALSE)
  b = v2%%A
  b[b==0] = A 
  a = (v2-1)%/%A+1
  for(i in 1:length(v2))
  {
    dataRing2[r,a[i],b[i]] = dataRing2[r,a[i],b[i]]+1
  }
}
dataRing = dataRing2

## Analysis using stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Amin=3;Amax=12
A = Amax-Amin+1
Astar = as.integer(mean(c(Amin,Amax)))
Astar = 7
alpha0=1
R1 = dim(dataRing)[1]
R2=length(unique(df$unitreaderindex))
R=R1+R2
Y=length(unique(df$stratayearindex))
S=dim(df$d)[1]
#Number of strata


## Initial values by running only on ring data
d.comb0 = list(R=R1,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
               readErr=dataRing,
               alpharep=rep(alpha0,A),
               taupar=sqrt(c(20,2,5,1,20)),
               deltapar=50,
               eps=0.001)

load("../save_results/fitring0.year50000.RData")

simext = rstan::extract(fitring0.year)
M = length(simext$Palpha0)
init1 = function(id)
{
  ind = sample(1:M,1)
  list(
    Palpha0=simext$Palpha0[ind],Palpha1=simext$Palpha1[ind],Pbeta0=simext$Pbeta0[ind],Pbeta1=simext$Pbeta1[ind],
       Ptrphi=simext$Ptrphi[ind],delta=simext$delta[ind,],
       tau=simext$tau[ind,],
       alpha0=c(simext$alpha0[ind,],rnorm(R2,simext$Palpha0[ind],sqrt(simext$tau[ind,1]))),
       beta0=c(simext$beta0[ind,],rnorm(R2,simext$Pbeta0[ind],sqrt(simext$tau[ind,2]))),
       alpha1=c(simext$alpha1[ind,],rnorm(R2,simext$Palpha1[ind],sqrt(simext$tau[ind,3]))),
       beta1=c(simext$beta1[ind,],rtnorm(R2,simext$Pbeta1[ind],sqrt(simext$tau[ind,4]),ub=0.0)),
       trphi=c(simext$trphi[ind,],rnorm(R2,simext$Ptrphi[ind],sqrt(simext$tau[ind,5]))))
}

init2 = lapply(1:4,function(id) init1(id=id))

d.comb = list(Y=Y,S=df$S,K=df$K,U=df$U,
              R1=R1,R2=R2,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
              readErr=dataRing,
              D=df$d,N=df$N,
              stratayearindex=df$stratayearindex,
              stationstrataindex=df$stationstrataindex,
              unitstationindex=df$unitstationindex,
              unitreaderindex=df$unitreaderindex,
              alpharep=rep(alpha0,A),
              taupar=c(50,2,5,1,20),
              deltapar=50)

fitcomb0subring8.year = stan(file="stan/agereader_comb_year_strat0.stan",data=d.comb,
                    iter=200000,chains=4,thin=100,init=init2)
save(fitcomb0subring8.year,file="../save_results/fitcomb0subring8.year200000.RData")


