sample.clust.id=function(theta,phi,nplot,array.gskp,nspp,ngrid,nclust){
  
  for (s in 1:nspp){
    array.gkp=array.gskp[,s,,]
    tmp=SampleClustID(theta=theta, phi=phi[,s],nplot=nplot, 
                      ngrid=ngrid, nclust=nclust, 
                      ArrayGKP=array.gkp)
    array.gskp[,s,,]=tmp$ArrayGKP1
  }
  array.gskp
}

sample.plot.id=function(theta,delta,array.gskp,ngrid,nclust,nplot,nspp){
  SomaGK1=matrix(0,ngrid,nclust) #useful to sample sig2
  for (s in 1:nspp){
    array.gkp=array.gskp[,s,,]
    SomaGK=rowSums(array.gkp,dims=2)
    SomaK=colSums(SomaGK)
    tmp=SamplePlotID(theta=theta, delta=delta, 
                     ngrid=ngrid, nclust=nclust, nplot=nplot,
                     ArrayGKP=array.gkp,
                     SomaGK=SomaGK,
                     SomaK=SomaK)
    array.gskp[,s,,]=tmp$ArrayGKP1
    SomaGK1=SomaGK1+SomaGK
  }
  list(array.gskp=array.gskp,SomaGK=SomaGK1)
}

get.delta=function(sig2,dist1){
  k=1/(2*sig2)
  tmp=exp(-k*dist1)
  tmp/rowSums(tmp)
}

sample.sig2=function(dist1,sig2,theta,jump.sd,ngrid,nclust,soma.gk){
  #propose new values
  sd1.old=sqrt(sig2)
  tmp=rnorm(1,mean=sd1.old,sd=jump.sd)
  sd1.new=ifelse(tmp<0,-tmp,tmp)
  
  #get deltas
  delta.old=get.delta(sig2=sd1.old^2,dist1=dist1)
  delta.new=get.delta(sig2=sd1.new^2,dist1=dist1)
  ldelta.old.theta=log(delta.old%*%theta)
  ldelta.new.theta=log(delta.new%*%theta)
  p.old=sum(soma.gk*ldelta.old.theta)
  p.new=sum(soma.gk*ldelta.new.theta)

  #accept or reject
  cond=is.na(p.new) #this happens when sd is very small and we divide by zero in get.delta()
  if (cond)  tmp=0
  if (!cond) tmp=exp(p.new-p.old)
  sd1=sd1.old
  accept1=0
  if (runif(1)<tmp) {sd1=sd1.new; accept1=1}
  list(sig2=sd1^2,accept1=accept1)
}

sample.theta=function(soma.skp,nplot,nclust,gamma1.theta,nspp){
  soma.pk=t(colSums(soma.skp,dims=1))
  soma=soma.pk+gamma1.theta
  
  #sample dirichlet
  tmp=matrix(rgamma(nplot*nclust,soma),nplot,nclust)
  tmp/rowSums(tmp)
  # for (i in 1:nplot) theta[i,]=rdirichlet(1,soma[i,])
  # theta
}

sample.phi=function(soma.skp,gamma1.phi,nclust,nspp,nplot){
  phi=matrix(0,nclust,nspp)
  soma.ks=t(rowSums(soma.skp,dims=2))
  soma=soma.ks+gamma1.phi
  
  #sample dirichlet
  tmp=matrix(rgamma(nspp*nclust,soma),nclust,nspp)
  tmp/rowSums(tmp)
  # for (i in 1:nclust) phi[i,]=rdirichlet(1,soma[i,])
  # phi
}

get.llk=function(dat,delta,theta,phi,ngrid,nspp,nomes.spp){
  tmp=delta%*%theta%*%phi
  lprob=log(tmp)
  sum(dat[,nomes.spp]*lprob)
}