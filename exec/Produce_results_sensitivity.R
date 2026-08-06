library(AgeReading)
library(rstan)
library(ggplot2)
require(gridExtra)
library(loo)
library(bayesplot)

Amin=3;Amax=12
A = Amax-Amin+1
R = 14

dir = "../save_results/"
load(paste0(dir,"fitcomb0.year100000.RData"))              #reader dependent age errors
load(paste0(dir,"fitcomb0.year.sens100000.RData"))             #reader independent age errors

#Harmonic mean for marginal likelihood
foo0 = as.vector(extract(fitcomb0.year)$loglik)
ML0 = max(foo0)-log(mean(exp(-foo0+max(foo0))))
foo0sens = as.vector(extract(fitcomb0.year.sens)$loglik)
ML0sens = max(foo0sens)-log(mean(exp(-foo0sens+max(foo0sens))))

##Calculating number at age
data(Herring_D_boot)
data(HerringSurvey)

Bootmean = Herring_D_boot[Herring_D_boot$BootstrapID<2,]
nunit = length(unique(Bootmean$unitstrat))
#Only using mean of bootstrap samples
for(i in 1:nunit)
{
  Bootmean$Abundance[Bootmean$unitstrat==i] = mean(Herring_D_boot$Abundance[Herring_D_boot$unitstrat==i])
  Bootmean$Biomass[Bootmean$unitstrat==i] = mean(Herring_D_boot$Biomass[Herring_D_boot$unitstrat==i])
}


foo = rstan::extract(fitcomb0.year)$xi1
res0 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)

foo = rstan::extract(fitcomb0.year.sens)$xi1
res0sens = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)



##Rename xi0 so ages become right
## Commented out, already renamed
if(1)
{
library(stringr)
ind = grep("xi0",names(fitcomb0.year.sens))
foo = names(fitcomb0.year.sens)[ind]
for(i in 10:1)
{
  foo = str_replace_all(foo,paste0(",",i,"]"),paste0(",",i+2,"]"))
}
names(fitcomb0.year.sens)[ind] = foo
foo2 = list()
foo3 = substring(foo,first=5,last=5)
foo4 = substring(foo,first=7,last=8)
foo4[1:15] = substring(foo4[1:15],first=1,last=1)
for(i in 1:length(foo))
  foo2[[i]] = substitute(expression(x[i,j]^0),list(i=foo3[i],j=foo4[i]))
}

##Rename M1 so ages become right
## Commented out, already renamed
if(1)
{
ind = grep("M1",names(fitcomb0.year.sens))
ind2 = matrix(nrow=length(ind),ncol=3)
i=0
for(a in 1:A)
  for(a2 in 1:A)
    for(r in 1:R)
    {
      i = i+1
      ind2[i,] = c(r,a2+Amin-1,a+Amin-1)
    }
names(fitcomb0.year.sens)[ind] = paste0("M1[",ind2[,1],",",ind2[,2],",",ind2[,3],"]")
}

M1 = extract(fitcomb0.year,pars="M1")$M1
M1mean = apply(M1,c(2,3,4),mean)

mean0 = round(apply(res0[,2,],2,mean)/10^9,2)
mean0sens = round(apply(res0sens[,2,],2,mean)/10^9,2)
sd0 = round(apply(res0[,2,],2,sd)/10^9,2)
sd0sens = round(apply(res0sens[,2,],2,sd)/10^9,2)

l = list()
ind = 0
y = 2
for(a in 3:5)
{
  ind = ind+1
  df = data.frame(y=c(res0[,y,a],res0sens[,y,a]),
                  Method=c(rep(1,dim(res0)[1]),
                       rep(2,dim(res0sens)[1])))
  df$y = df$y/10^9
  df$Method = as.factor(df$Method)
  l[[ind]] = my_ggplot2(df,a,mean0,mean0sens,sd0,sd0sens,Amin)

}
pdf("Nest_sens_2022.pdf",height=5,width=10)
grid.arrange(l[[1]],l[[2]],l[[3]],ncol=3,nrow=1)
dev.off()

pdf("fit0_xi0_sens_2021.pdf",height=4,width=4)
pars = paste0("xi0[1,",2+1:A,"]")
xx =as.matrix(fitcomb0.year,pars=pars)
xx2 = as.matrix(fitcomb0.year.sens,pars=pars)
xxcomb = matrix(NA,nrow=nrow(xx),ncol=2*ncol(xx))
colnames(xxcomb) = 1:(2*ncol(xx))
for(j in 1:ncol(xx))
{
  print(j)
  xxcomb[,2*j-1] = xx[,j]
  xxcomb[,2*j] = xx2[,j]
  colnames(xxcomb)[c(2*j-1,2*j)] = c(paste0("A-",j),paste0("B-",j))
}
mcmc_intervals(xxcomb,prob=0.8,prob_outer=0.95)+labs(y="Age",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
dev.off()

pdf("fit0_xi0_sens_2022.pdf",height=4,width=4)
pars = paste0("xi0[2,",2+1:A,"]")
xx =as.matrix(fitcomb0.year,pars=pars)
xx2 = as.matrix(fitcomb0.year.sens,pars=pars)
xxcomb = matrix(NA,nrow=nrow(xx),ncol=2*ncol(xx))
colnames(xxcomb) = 1:(2*ncol(xx))
for(j in 1:ncol(xx))
{
  print(j)
  xxcomb[,2*j-1] = xx[,j]
  xxcomb[,2*j] = xx2[,j]
  colnames(xxcomb)[c(2*j-1,2*j)] = c(paste0("A-",j),paste0("B-",j))
}
mcmc_intervals(xxcomb,prob=0.8,prob_outer=0.95)+labs(y="Age",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
dev.off()

pdf("fit0_xi0_sens_2023.pdf",height=4,width=4)
pars = paste0("xi0[3,",2+1:A,"]")
xx =as.matrix(fitcomb0.year,pars=pars)
xx2 = as.matrix(fitcomb0.year.sens,pars=pars)
xxcomb = matrix(NA,nrow=nrow(xx),ncol=2*ncol(xx))
colnames(xxcomb) = 1:(2*ncol(xx))
for(j in 1:ncol(xx))
{
  print(j)
  xxcomb[,2*j-1] = xx[,j]
  xxcomb[,2*j] = xx2[,j]
  colnames(xxcomb)[c(2*j-1,2*j)] = c(paste0("A-",j),paste0("B-",j))
}
mcmc_intervals(xxcomb,prob=0.8,prob_outer=0.95)+labs(y="Age",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
dev.off()



pdf("fit0_M1_sens_11_11.pdf",height=4.5,width=4)
ind = 1:R
pars = paste0("M1[",ind,",11,11]")
xx =as.matrix(fitcomb0.year,pars=pars)
xx2 = as.matrix(fitcomb0.year.sens,pars=pars)
xxcomb = matrix(NA,nrow=nrow(xx),ncol=2*ncol(xx))
colnames(xxcomb) = 1:(2*ncol(xx))
for(j in 1:ncol(xx))
{
  print(j)
  xxcomb[,2*j-1] = xx[,j]
  xxcomb[,2*j] = xx2[,j]
  colnames(xxcomb)[c(2*j-1,2*j)] = c(paste0("A-",j),paste0("B-",j))
}
mcmc_intervals(xxcomb,prob=0.8,prob_outer=0.95)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
dev.off()

pdf("fit0_M1_sens_45.pdf",height=4.5,width=4)
ind = 1:R
pars = paste0("M1[",ind,",4,5]")
xx =as.matrix(fitcomb0.year,pars=pars)
xx2 = as.matrix(fitcomb0.year.sens,pars=pars)
xxcomb = matrix(NA,nrow=nrow(xx),ncol=2*ncol(xx))
colnames(xxcomb) = 1:(2*ncol(xx))
for(j in 1:ncol(xx))
{
  print(j)
  xxcomb[,2*j-1] = xx[,j]
  xxcomb[,2*j] = xx2[,j]
  colnames(xxcomb)[c(2*j-1,2*j)] = c(paste0("A-",j),paste0("B-",j))
}
mcmc_intervals(xxcomb,prob=0.8,prob_outer=0.95)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
dev.off()

pdf("fit0_M1_sens_76.pdf",height=4.5,width=4)
ind = 1:R
pars = paste0("M1[",ind,",7,6]")
xx =as.matrix(fitcomb0.year,pars=pars)
xx2 = as.matrix(fitcomb0.year.sens,pars=pars)
xxcomb = matrix(NA,nrow=nrow(xx),ncol=2*ncol(xx))
colnames(xxcomb) = 1:(2*ncol(xx))
for(j in 1:ncol(xx))
{
  print(j)
  xxcomb[,2*j-1] = xx[,j]
  xxcomb[,2*j] = xx2[,j]
  colnames(xxcomb)[c(2*j-1,2*j)] = c(paste0("A-",j),paste0("B-",j))
}
mcmc_intervals(xxcomb,prob=0.8,prob_outer=0.95)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
dev.off()
