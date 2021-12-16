#get estimated values
theta.estim=theta
phi.estim=phi
clust.id.estim=dat$clust.id
plot.id.estim=dat$plot.id
delta.estim=delta
sig2.estim=sig2
  
#get ordem
fim=data.frame(true1=dat0$psi,estim1=clust.id.estim)
k=table(fim)
ordem=c(2,3,1)

#look at theta
theta.estim1=theta.estim[,ordem]
rango=range(c(theta.true,theta.estim1))
plot(theta.true,theta.estim1,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

#look at phi
phi.estim1=phi.estim[ordem,]
rango=range(c(phi.true,phi.estim1))
plot(phi.true,phi.estim1,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

#look at cluster id
clust.id.estim1=clust.id.estim
for (i in 1:length(ordem)){
  cond=clust.id.estim==ordem[i]
  clust.id.estim1[cond]=i
}
fim=data.frame(true1=dat0$psi,estim1=clust.id.estim1)
k=table(fim)
k/apply(k,2,sum)

#look at plot id
fim=data.frame(true1=dat0$omega,estim1=dat$plot.id)
k=table(fim)
k/apply(k,2,sum)

#---------------------------------
#look at the spatial distribution of each group
#calculate probabilities
coord=expand.grid(x=seq(from=0,to=1000,by=10),
                  y=seq(from=0,to=1000,by=10))
ncoord=nrow(coord)

#get distances
ndoc=nrow(coord.doc)
dist1=matrix(NA,ncoord,ndoc)
for (i in 1:ndoc){
  x2=(coord$x-coord.doc$x[i])^2
  y2=(coord$y-coord.doc$y[i])^2
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
colnames(res)=paste0('c',1:nclust)
res1=as.data.frame(res)
coord1=cbind(coord,res1)

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
