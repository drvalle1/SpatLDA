sample.clust.id=function(theta,phi,nplot,array.gkp){
  soma.gp=apply(array.gkp,c(1,3),sum)
  soma.p=apply(soma.gp,1,sum)
  for (i in 1:nplot){
    n=soma.p[i]
    if (n>0){
      tmp=phi*theta[i,]  
      prob=tmp/sum(tmp)
      for (g in 1:ngrid){
        n=soma.gp[g,i]
        if (n>0) array.gkp[g,,i]=rmultinom(1,size=n,prob=prob)  
      }
    }
  }
  array.gkp
}

sample.plot.id=function(theta,delta,array.gkp,ngrid,nclust){
  soma.gk=apply(array.gkp,c(1,2),sum)
  soma.k=apply(soma.gk,2,sum)
  for (k in 1:nclust){
    n=soma.k[k]
    if (n>0){
      for (g in 1:ngrid){
        n=soma.gk[g,k]
        if (n>0){
          tmp=theta[,k]*delta[g,]  
          prob=tmp/sum(tmp)
          array.gkp[g,k,]=rmultinom(1,size=n,prob=prob)
        }
      }
    }
  }
  array.gkp
}

get.delta=function(sig2,dist1){
  k=1/(2*sig2)
  tmp=exp(-k*dist1)
  tmp/rowSums(tmp)
}

sample.sig2=function(dist1,sig2,theta,jump.sd,ngrid,nclust,array.gskp){
  #summarize array
  soma.gk=apply(array.gskp,c(1,3),sum)
  
  #propose new values
  sd1.old=sqrt(sig2)
  tmp=rnorm(1,mean=sd1.old,sd=jump.sd)
  sd1.new=ifelse(tmp<0,-tmp,tmp)
  
  #get deltas
  delta.old=get.delta(sig2=sd1.old^2,dist1=dist1)
  delta.new=get.delta(sig2=sd1.new^2,dist1=dist1)
  ldelta.old.theta=log(delta.old%*%theta)
  ldelta.new.theta=log(delta.new%*%theta)
  p.old=p.new=0
  for (g in 1:ngrid){
    for (k in 1:nclust){
      if (soma.gk[g,k]>0){
        p.old=p.old+soma.gk[g,k]*ldelta.old.theta[g,k]
        p.new=p.new+soma.gk[g,k]*ldelta.new.theta[g,k]
      }
    }
  }
  
  #accept or reject
  cond=is.na(p.new) #this happens when sd is very small and we divide by zero in get.delta()
  if (cond)  tmp=0
  if (!cond) tmp=exp(p.new-p.old)
  sd1=sd1.old
  accept1=0
  if (runif(1)<tmp) {sd1=sd1.new; accept1=1}
  list(sig2=sd1^2,accept1=accept1)
}

sample.theta=function(array.gskp,nplot,nclust,gamma1.theta){
  theta=matrix(0,nplot,nclust)
  soma.pk=t(apply(array.gskp,c(3,4),sum))
  soma=soma.pk+gamma1.theta
  for (i in 1:nplot) theta[i,]=rdirichlet(1,soma[i,])
  theta
}

sample.phi=function(array.gskp,gamma1.phi,nclust,nspp){
  phi=matrix(0,nclust,nspp)
  soma.ks=t(apply(array.gskp,c(2,3),sum))
  soma=soma.ks+gamma1.phi
  for (i in 1:nclust) phi[i,]=rdirichlet(1,soma[i,])
  phi
}

get.llk=function(dat,delta,theta,phi,ngrid,nspp,nomes.spp){
  lprob=matrix(NA,ngrid,nspp)
  for (i in 1:ngrid){
    for (j in 1:nspp){
      tmp=delta[i,]%*%theta%*%phi[,j]
      lprob[i,j]=log(tmp)
    }
  }
  sum(dat[,nomes.spp]*lprob)
}