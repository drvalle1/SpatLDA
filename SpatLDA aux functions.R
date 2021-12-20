sample.clust.id=function(theta,phi,nplot,nspp,array.gskp){
  array.soma=apply(array.gskp,c(1,2,4),sum)
  array.teste=apply(array.soma,c(2,3),sum)
  for (i in 1:nplot){
    for (j in 1:nspp){
      n=array.teste[j,i]
      if (n>0){
        tmp=phi[,j]*theta[i,]  
        prob=tmp/sum(tmp)
        for (g in 1:ngrid){
          n=array.soma[g,j,i]
          if (n>0) array.gskp[g,j,,i]=rmultinom(1,size=n,prob=prob)  
        }
      }
    }
  }      
  array.gskp
}

sample.plot.id=function(theta,delta,array.gskp,ngrid,nspp,nclust){
  array.soma=apply(array.gskp,c(1,2,3),sum)
  array.test=apply(array.soma,c(2,3),sum)
  for (i in 1:nspp){
    for (j in 1:nclust){
      n=array.test[i,j]
      if (n>0){
        for (g in 1:ngrid){
          n=array.soma[g,i,j]
          if (n>0){
            tmp=theta[,j]*delta[g,]  
            prob=tmp/sum(tmp)
            array.gskp[g,i,j,]=rmultinom(1,size=n,prob=prob)
          }
        }
      }
    }
  }
  array.gskp
}

get.delta=function(sig2,dist1){
  k=1/(2*sig2)
  tmp=exp(-k*dist1)
  tmp/rowSums(tmp)
}

sample.sig2=function(dist1,sig2,theta,jump.sd,ngrid,nclust,array.gskp){
  #summarize array
  array.soma=apply(array.gskp,c(1,3),sum)
  
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
  for (i in 1:ngrid){
    for (j in 1:nclust){
      if (array.soma[i,j]>0){
        p.old=p.old+array.soma[i,j]*ldelta.old.theta[i,j]
        p.new=p.new+array.soma[i,j]*ldelta.new.theta[i,j]
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
  soma=t(apply(array.gskp,c(3,4),sum))+gamma1.theta
  for (i in 1:nplot) theta[i,]=rdirichlet(1,soma[i,])
  theta
}

sample.phi=function(array.gskp,gamma1.phi,nclust,nspp){
  phi=matrix(0,nclust,nspp)
  soma=t(apply(array.gskp,c(2,3),sum))+gamma1.phi
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