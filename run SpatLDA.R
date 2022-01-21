rm(list=ls(all=TRUE))
library('gtools')
library('Rcpp')
library('RcppArmadillo')
set.seed(1)

#get functions
setwd('U:\\GIT_models\\SpatLDA')
sourceCpp('SpatLDA_aux.cpp')
source('SpatLDA aux functions.R')
source('SpatLDA_main_func.R')

#get data
setwd('U:\\GIT_models\\SpatLDA\\fake data')
dat=read.csv('fake data1.csv',as.is=T)

#potential documents
tmp=unique(dat[,c('xbin','ybin')])
coord.plot=data.frame(x=tmp$xbin,y=tmp$ybin)

#useful stuff
nclust=3
ngibbs=3000
nburn=ngibbs/2

#priors
gamma.theta=0.1
gamma.phi=0.1

#run Gibbs sampler
mod=SpatLDA(nclust=nclust,ngibbs=ngibbs,nburn=nburn,
            coord.plot=coord.plot,dat=dat,
            gamma.theta=gamma.theta,gamma.phi=gamma.phi)

#look at convergence
plot(mod$llk,type='l')
plot(mod$sig2,type='l')
seq1=seq(from=1000,to=ngibbs,length.out=1000)
plot(mod$llk[seq1],type='l')
plot(mod$sig2[seq1],type='l')

#export results
setwd('U:\\GIT_models\\SpatLDA\\results')
nomes=paste0(c('theta','phi','sig2','llk'),'.csv')
write.csv(mod$theta[seq1,],nomes[1],row.names=F)
write.csv(mod$phi[seq1,],  nomes[2],row.names=F)
write.csv(mod$sig2,        nomes[3],row.names=F)
write.csv(mod$llk,         nomes[4],row.names=F)