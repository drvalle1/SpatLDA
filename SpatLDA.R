rm(list=ls(all=TRUE))
library('gtools')
# library('Rcpp')
set.seed(1)

#get functions
setwd('U:\\GIT_models\\SpatLDA')
source('SpatLDA aux functions.R')
# sourceCpp('aux1.cpp')

#get data
setwd('U:\\GIT_models\\SpatLDA\\fake data')
dat0=dat=read.csv('fake data.csv',as.is=T)
# ind=which(colnames(dat)=='omega'); colnames(dat)[ind]='plot.id'
# ind=which(colnames(dat)=='psi'); colnames(dat)[ind]='clust.id'

#potential documents
coord.doc=expand.grid(x=seq(from=0,to=1000,by=200),
                      y=seq(from=0,to=1000,by=200))
# coord.doc=data.frame(x=c(200,600,800),
#                      y=c(200,400,800))
ndoc=nrow(coord.doc)

#useful stuff
nclust=3
ngibbs=1000
nburn=ngibbs/2

#priors
gamma.theta=0.1
gamma.phi=0.1

#useful stuff
ntrees=nrow(dat)
nspp=max(dat$spp)

#get distances
dist1=matrix(NA,ntrees,ndoc)
for (i in 1:ndoc){
  x2=(dat$x-coord.doc$x[i])^2
  y2=(dat$y-coord.doc$y[i])^2
  dist1[,i]=sqrt(x2+y2)
}

#initial values
theta=matrix(1/nclust,ndoc,nclust)
phi=matrix(1/nspp,nclust,nspp)
sig2=100
delta=get.delta(sig2=sig2,dist1=dist1)
dat$clust.id=sample(1:nclust,size=ntrees,replace=T)
dat$plot.id=sample(1:ndoc,size=ntrees,replace=T)

#MH stuff
jump.sd=1 
accept1=0
nchange.jump=50

#to store outcomes from gibbs sampler
theta.out=matrix(NA,ngibbs,nclust*ndoc)
phi.out=matrix(NA,ngibbs,nclust*nspp)
sig2.out=matrix(NA,ngibbs,1)
llk=rep(NA,ngibbs)
  
#run gibbs sampler
options(warn=2)
for (i in 1:ngibbs){
  print(i)   

  #sample clust.id
  dat=sample.clust.id(theta=theta,phi=phi,ntrees=ntrees,dat=dat)
  
  #sample plot.id
  dat=sample.plot.id(theta=theta,delta=delta,dat=dat,ntrees=ntrees)
  
  #sample theta
  theta=sample.theta(dat=dat,ndoc=ndoc,nclust=nclust,gamma1.theta=gamma.theta)
  # theta=theta.true
  
  #sample phi
  phi=sample.phi(dat=dat,gamma1.phi=gamma.phi,nclust=nclust,nspp=nspp)
  # phi=phi.true
  
  #sample sig2
  tmp=sample.sig2(dist1=dist1,sig2=sig2,theta=theta,jump.sd=jump.sd)
  accept1=accept1+tmp$accept1
  sig2=tmp$sig2
  # sig2=sig2.true
  delta=get.delta(sig2=sig2,dist1=dist1)
  
  #adapt MH
  if (i%%nchange.jump==0 & i<nburn){
    if (accept1/nchange.jump < 0.1) jump.sd=jump.sd*0.5
    if (accept1/nchange.jump > 0.6) jump.sd=jump.sd*2
    accept1=0
  } 
  
  #llk
  llk[i]=get.llk(dat=dat,delta=delta,theta=theta,phi=phi,ntrees=ntrees)
  
  #store results  
  theta.out[i,]=theta
  phi.out[i,]=phi
  sig2.out[i]=sig2
}

plot(llk,type='l')
plot(sig2.out,type='l')