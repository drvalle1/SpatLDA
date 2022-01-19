rm(list=ls())
library('ggplot2')
set.seed(3)

#trees
ntrees=10000
coord=data.frame(x=runif(ntrees,min=0,max=1000),
                 y=runif(ntrees,min=0,max=1000))

#documents
ndoc=3
coord.doc=data.frame(x=c(200,600,800),
                     y=c(200,400,800))
plot(y~x,data=coord)
points(y~x,data=coord.doc,col='red',pch=19)

#theta for each doc
nclust=3
theta=matrix(c(0.05,0.05,0.9,
               0.05,0.9,0.05,
               0.8,0.1,0.1),ndoc,nclust,byrow=T)
theta.true=theta
# apply(theta,1,sum)

#phi for each cluster
nspp=50
tmp=rbeta(nspp*nclust,0.1,0.1)
ind=sample(1:(nspp*nclust),size=25)
tmp[ind]=4
tmp1=matrix(tmp,nclust,nspp)
phi=tmp1/apply(tmp1,1,sum)
phi.true=phi
# head(round(phi[,1:10],3))
# hist(phi[5,])

#get distance from each document to each tree
dist1=matrix(NA,ntrees,ndoc)
for (i in 1:ndoc){
  x2=(coord$x-coord.doc$x[i])^2
  y2=(coord$y-coord.doc$y[i])^2
  dist1[,i]=sqrt(x2+y2)
}

#get cluster membership for each tree
sig2=sig2.true=20
for (i in 1:ntrees){
  #which document?
  prob=exp(-(1/(2*sig2))*dist1[i,])
  tmp=rmultinom(1,size=1,prob=prob/sum(prob))
  ind=which(tmp==1)
  coord$omega[i]=ind
  
  #which cluster?
  tmp=rmultinom(1,size=1,prob=theta[ind,])
  ind=which(tmp==1)
  coord$psi[i]=ind
  
  #which species?
  tmp=rmultinom(1,size=1,prob=phi[ind,])
  ind=which(tmp==1)
  coord$spp[i]=ind
}

#spatial distribution of communities
plot(y~x,data=coord,col=coord$psi,pch=19)
plot(y~x,data=coord,col=coord$omega,pch=19)

#export results
setwd('U:\\GIT_models\\SpatLDA\\fake data')
# ind=which(colnames(coord)%in%c('psi','omega'))
# write.csv(coord[,-ind],'fake data.csv',row.names=F)
write.csv(coord,'fake data.csv',row.names=F)
#--------------------------------
#based on these settings, 
#what would a map of the probability of each cluster look like?
coord=expand.grid(x=seq(from=0,to=1000,by=10),
                  y=seq(from=0,to=1000,by=10))
ncoord=nrow(coord)

#get distances
dist1=matrix(NA,ncoord,ndoc)
for (i in 1:ndoc){
  x2=(coord$x-coord.doc$x[i])^2
  y2=(coord$y-coord.doc$y[i])^2
  dist1[,i]=sqrt(x2+y2)
}

#calculate probabilities
res=matrix(NA,ncoord,nclust)
for (i in 1:ncoord){
  tmp=exp(-(1/(2*sig2))*dist1[i,])
  delta=tmp/sum(tmp)
  for (j in 1:nclust){
    res[i,j]=sum(theta[,j]*delta)
  }
}

#create final data.frame with prob for each cluster
colnames(res)=paste0('c',1:nclust)
res1=as.data.frame(res)
coord1=cbind(coord,res1)
setwd('U:\\GIT_models\\SpatLDA\\fake data')
write.csv(coord1,'fake data proportions.csv',row.names=F)

#plot results
library(gridExtra)
library(grid)
p1=ggplot() +
  geom_tile(data = coord1, alpha = 0.8,aes(x = x, y = y,fill = c1)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) 
p2=ggplot() +
  geom_tile(data = coord1, alpha = 0.8,aes(x = x, y = y,fill = c2)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) 
p3=ggplot() +
  geom_tile(data = coord1, alpha = 0.8,aes(x = x, y = y,fill = c3)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) 
grid.arrange(p1, p2, p3,nrow = 1)
