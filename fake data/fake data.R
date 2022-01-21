rm(list=ls())
library('ggplot2')
set.seed(3)

#trees
ntrees=10000
coord=data.frame(x=runif(ntrees,min=0,max=1000),
                 y=runif(ntrees,min=0,max=1000))

#center of each group
center=data.frame(x=c(300,500,200),
                  y=c(100,500,800))
plot(y~x,data=center,xlim=c(0,1000),ylim=c(0,1000))

#theta
nclust=3
sd1=300 #originally 100
theta=data.frame(gr1=dnorm(coord$x,mean=center$x[1],sd=sd1)*
                     dnorm(coord$y,mean=center$y[1],sd=sd1),
                 gr2=dnorm(coord$x,mean=center$x[2],sd=sd1)*
                     dnorm(coord$y,mean=center$y[2],sd=sd1),
                 gr3=dnorm(coord$x,mean=center$x[3],sd=sd1)*
                     dnorm(coord$y,mean=center$y[3],sd=sd1))
theta=theta/apply(theta,1,sum)
apply(theta,2,range)
theta.true=theta
# apply(theta,1,sum)

par(mfrow=c(3,1),mar=rep(1,4))
plot(y~x,data=coord,col=grey(theta$gr1),pch=19)
points(y~x,data=center,pch=19,col='red')
plot(y~x,data=coord,col=grey(theta$gr2),pch=19)
points(y~x,data=center,pch=19,col='red')
plot(y~x,data=coord,col=grey(theta$gr3),pch=19)
points(y~x,data=center,pch=19,col='red')

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

#get cluster membership for each tree
for (i in 1:ntrees){
  tmp=rmultinom(1,size=1,prob=theta[i,])
  ind=which(tmp==1)
  coord$psi[i]=ind
  
  #which species?
  tmp=rmultinom(1,size=1,prob=phi[ind,])
  ind=which(tmp==1)
  coord$spp[i]=ind
}

#spatial distribution of communities
plot(y~x,data=coord,col=coord$psi,pch=19)

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

theta=data.frame(gr1=dnorm(coord$x,mean=center$x[1],sd=sd1)*
                   dnorm(coord$y,mean=center$y[1],sd=sd1),
                 gr2=dnorm(coord$x,mean=center$x[2],sd=sd1)*
                   dnorm(coord$y,mean=center$y[2],sd=sd1),
                 gr3=dnorm(coord$x,mean=center$x[3],sd=sd1)*
                   dnorm(coord$y,mean=center$y[3],sd=sd1))
theta=theta/apply(theta,1,sum)
coord1=cbind(coord,theta)

setwd('U:\\GIT_models\\SpatLDA\\fake data')
write.csv(coord1,'fake data proportions.csv',row.names=F)
write.csv(phi.true,'fake data phi.csv',row.names=F)

#plot results
library(gridExtra)
library(grid)
p1=ggplot() +
  geom_tile(data = coord1, alpha = 0.8,aes(x = x, y = y,fill = gr1)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) 
p2=ggplot() +
  geom_tile(data = coord1, alpha = 0.8,aes(x = x, y = y,fill = gr2)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) 
p3=ggplot() +
  geom_tile(data = coord1, alpha = 0.8,aes(x = x, y = y,fill = gr3)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) 
grid.arrange(p1, p2, p3,nrow = 1)
