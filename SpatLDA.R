rm(list=ls(all=TRUE))
library('gtools')
set.seed(1)

#get functions
setwd('U:\\GIT_models\\SpatLDA')
source('SpatLDA aux functions.R')

#get data
setwd('U:\\GIT_models\\SpatLDA\\fake data')
dat0=dat=read.csv('fake data1.csv',as.is=T)

#potential documents
coord.plot=expand.grid(x=seq(from=0,to=1000,by=200),
                       y=seq(from=0,to=1000,by=200))
nplot=nrow(coord.plot)

#useful stuff
nclust=3
ngibbs=10000
nburn=ngibbs/2

#priors
gamma.theta=0.1
gamma.phi=0.1

#useful stuff
ngrid=nrow(dat)
ind=grep('spp',colnames(dat))
nomes.spp=colnames(dat)[ind]
tmp1=as.numeric(gsub('spp','',nomes.spp))
nspp=max(tmp1)

#get distances
dist1=matrix(NA,ngrid,nplot)
for (i in 1:nplot){
  x2=(dat$xbin-coord.plot$x[i])^2
  y2=(dat$ybin-coord.plot$y[i])^2
  dist1[,i]=sqrt(x2+y2)
}

#initial values for parameters
theta=matrix(1/nclust,nplot,nclust)
phi=matrix(1/nspp,nclust,nspp)
sig2=100
delta=get.delta(sig2=sig2,dist1=dist1)

#initial values for cluster and plot membership
array.gskp=array(0,dim=c(ngrid,nspp,nclust,nplot))
nbins=nclust*nplot
for (i in 1:ngrid){
  for (j in 1:nspp){
    tmp=dat[i,paste0('spp',j)]
    if (tmp>0){
      tmp1=rmultinom(1,size=tmp,prob=rep(1/nbins,nbins))
      array.gskp[i,j,,]=matrix(tmp1,nclust,nplot) 
    }
  }
}
# teste=apply(array.gskp,c(1,2),sum)
# unique(teste-dat[,paste0('spp',1:nspp)])

#MH stuff
jump.sd=1 
accept1=0
nchange.jump=50

#to store outcomes from gibbs sampler
theta.out=matrix(NA,ngibbs,nclust*nplot)
phi.out=matrix(NA,ngibbs,nclust*nspp)
sig2.out=matrix(NA,ngibbs,1)
llk=rep(NA,ngibbs)
  
#run gibbs sampler
options(warn=2)
for (i in 1:ngibbs){
  print(i)   

  #sample clust.id
  for (s in 1:nspp){
    array.gkp=array.gskp[,s,,]
    array.gskp[,s,,]=sample.clust.id(theta=theta,phi=phi[,s],nplot=nplot,
                                     array.gkp=array.gkp)
  }
  
  #sample plot.id
  for (s in 1:nspp){
    array.gkp=array.gskp[,s,,]
    array.gskp[,s,,]=sample.plot.id(theta=theta,delta=delta,array.gkp=array.gkp,
                                    ngrid=ngrid,nclust=nclust)
  }
  
  #sample theta
  theta=sample.theta(array.gskp=array.gskp,nplot=nplot,nclust=nclust,gamma1.theta=gamma.theta)
  # theta=theta.true
  
  #sample phi
  phi=sample.phi(array.gskp=array.gskp,gamma1.phi=gamma.phi,nclust=nclust,nspp=nspp)
  # phi=phi.true
  
  #sample sig2
  tmp=sample.sig2(dist1=dist1,sig2=sig2,theta=theta,jump.sd=jump.sd,
                  ngrid=ngrid,nclust=nclust,array.gskp=array.gskp)
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
  llk[i]=get.llk(dat=dat,delta=delta,theta=theta,phi=phi,
                 ngrid=ngrid,nspp=nspp,nomes.spp=nomes.spp)
  
  #store results  
  theta.out[i,]=theta
  phi.out[i,]=phi
  sig2.out[i]=sig2
}

plot(llk,type='l')
plot(sig2.out,type='l')