rm(list=ls())

#get data
setwd('U:\\GIT_models\\SpatLDA\\fake data')
dat=read.csv('fake data1.csv',as.is=T)
ind=grep('spp',colnames(dat))
nspp=length(ind)
ngrid=nrow(dat)

#potential documents
tmp=unique(dat[,c('xbin','ybin')])
coord.plot=data.frame(x=tmp$xbin,y=tmp$ybin)
nplot=nrow(coord.plot)
nclust=3

#get functions
setwd('U:\\GIT_models\\SpatLDA')
source('SpatLDA aux functions.R')

#get estimated values
setwd('U:\\GIT_models\\SpatLDA\\results')
tmp=read.csv('theta.csv')
theta.estim=matrix(unlist(tmp[nrow(tmp),]),nplot,nclust)
tmp=read.csv('phi.csv')
phi.estim=matrix(unlist(tmp[nrow(tmp),]),nclust,nspp)
tmp=read.csv('sig2.csv')
sig2.estim=tmp$V1[nrow(tmp)]

#get phi.true
setwd('U:\\GIT_models\\SpatLDA\\fake data')
phi.true=data.matrix(read.csv('fake data phi.csv'))

#get ordem
fim=matrix(NA,nclust,nclust)
for (i in 1:nclust){
  for (j in 1:nclust){
    tmp=cbind(phi.true[i,],phi.estim[j,])
    fim[i,j]=cor(tmp)[1,2]
  }
}
fim
ordem=c(2,1,3)
fim[,ordem]

#look at phi
phi.estim1=phi.estim[ordem,]
rango=range(c(phi.true,phi.estim1))
plot(phi.true,phi.estim1,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

#---------------------------------
#look at the spatial distribution of each group
#calculate probabilities
coord=expand.grid(x=seq(from=0,to=1000,by=10),
                  y=seq(from=0,to=1000,by=10))
ncoord=nrow(coord)

#get distances
nplot=nrow(coord.plot)
dist1=matrix(NA,ncoord,nplot)
for (i in 1:nplot){
  x2=(coord$x-coord.plot$x[i])^2
  y2=(coord$y-coord.plot$y[i])^2
  dist1[,i]=sqrt(x2+y2)
}

#calculate probabilities
res=matrix(NA,ncoord,nclust)
for (i in 1:ncoord){
  tmp=exp(-(1/(2*sig2.estim))*dist1[i,])
  delta=tmp/sum(tmp)
  for (j in 1:nclust){
    res[i,j]=sum(theta.estim[,j]*delta)
  }
}

#create final data.frame with prob for each cluster
res=res[,ordem]
colnames(res)=paste0('c',1:nclust)
res1=as.data.frame(res)
coord1=cbind(coord,res1)

#export results
setwd('U:\\GIT_models\\SpatLDA\\results')
write.csv(coord1,'estim proportions.csv',row.names=F)
