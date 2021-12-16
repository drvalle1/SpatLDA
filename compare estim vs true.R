theta.estim=theta
phi.estim=phi
clust.id.estim=dat$clust.id
plot.id.estim=dat$plot.id

ordem=c(3,1,2)
theta.estim1=theta.estim[,ordem]
rango=range(c(theta.true,theta.estim1))
plot(theta.true,theta.estim1,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

phi.estim1=phi.estim[ordem,]
rango=range(c(phi.true,phi.estim1))
plot(phi.true,phi.estim1,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

clust.id.estim1=clust.id.estim
for (i in 1:length(ordem)){
  cond=clust.id.estim==ordem[i]
  clust.id.estim1[cond]=i
}
fim=data.frame(true1=dat0$psi,estim1=clust.id.estim1)
k=table(fim)
k/apply(k,2,sum)

fim=data.frame(true1=dat0$omega,estim1=dat$plot.id)
k=table(fim)
k/apply(k,2,sum)