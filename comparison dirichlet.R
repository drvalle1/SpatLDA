rm(list=ls(all=TRUE))
library('gtools')

set.seed(1)

peso=100
alpha=c(0.1,0.6,0.3)*peso

#Dirichlet gtools
n=10000
k=rdirichlet(n,alpha=alpha)
colMeans(k)
apply(k,2,var)

#Dirichlet mine
ntipo=length(alpha)
tmp=matrix(NA,n,ntipo)
for (i in 1:n){
  tmp[i,]=rgamma(ntipo,alpha)
}
tmp1=tmp/rowSums(tmp)
colMeans(tmp1)
apply(tmp1,2,var)

#--------------------------------
n=1000
ntipos=3
alpha=matrix(NA,n,ntipos)
gtools.res=matrix(NA,n,ntipos)
peso=1000
for (i in 1:n){
  tmp=runif(ntipos)
  alpha[i,]=tmp*peso
  gtools.res[i,]=rdirichlet(1,alpha[i,])
}

tmp=matrix(rgamma(n*ntipos,alpha),n,ntipos)
my.res=tmp/rowSums(tmp)

rango=c(0,1)
plot(gtools.res,my.res,xlim=rango,ylim=rango)
lines(rango,rango,col='red')