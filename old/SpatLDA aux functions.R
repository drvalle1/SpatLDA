sample.clust.id=function(theta,phi,ntrees,dat){
  theta1=theta[dat$plot.id,]
  for (i in 1:ntrees){
    tmp=phi[,dat$spp[i]]*theta1[i,]
    prob1=tmp/sum(tmp)
    tmp2=rmultinom(1,size=1,prob=prob1)
    dat$clust.id[i]=which(tmp2==1)
  }
  dat
}

sample.plot.id=function(theta,delta,dat,ntrees){
  for (i in 1:ntrees){
    tmp=theta[,dat$clust.id[i]]*delta[i,]
    prob1=tmp/sum(tmp)
    tmp1=rmultinom(1,size=1,prob=prob1)
    dat$plot.id[i]=which(tmp1==1)
  }
  dat
}

get.delta=function(sig2,dist1){
  k=1/(2*sig2)
  tmp=exp(-k*dist1)
  tmp/rowSums(tmp)
}

sample.sig2=function(dist1,sig2,theta,jump.sd){
  sd1.old=sqrt(sig2)
  tmp=rnorm(1,mean=sd1.old,sd=jump.sd)
  sd1.new=ifelse(tmp<0,-tmp,tmp)
  
  #get deltas
  delta.old=get.delta(sig2=sd1.old^2,dist1=dist1)
  delta.new=get.delta(sig2=sd1.new^2,dist1=dist1)
  ldelta.old.theta=log(delta.old%*%theta)
  ldelta.new.theta=log(delta.new%*%theta)
  p.old=p.new=0
  for (i in 1:ntrees){
    p.old=p.old+ldelta.old.theta[i,dat$clust.id[i]]
    p.new=p.new+ldelta.new.theta[i,dat$clust.id[i]]
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

sample.theta=function(dat,ndoc,nclust,gamma1.theta){
  tab1=table(dat[,c('plot.id','clust.id')])
  nomes.col=as.numeric(colnames(tab1))
  nomes.row=as.numeric(rownames(tab1))
  theta=matrix(0,ndoc,nclust)
  for (i in 1:length(nomes.row)){
    tab2=rep(0,nclust)
    tab2[nomes.col]=tab1[i,]
    theta[nomes.row[i],]=rdirichlet(1,tab2+gamma1.theta)
  }
  theta
}

sample.phi=function(dat,gamma1.phi,nclust,nspp){
  tab1=table(dat[,c('clust.id','spp')])
  nums.col=as.numeric(colnames(tab1))
  nums.row=as.numeric(rownames(tab1))
  phi=matrix(0,nclust,nspp)
  for (i in 1:length(nums.row)){
    tab2=rep(0,nspp)
    tab2[nums.col]=tab1[i,]
    phi[nums.row[i],]=rdirichlet(1,tab2+gamma1.phi)
  }
  phi
}

get.llk=function(dat,delta,theta,phi,ntrees){
  llk=rep(NA,ntrees)
  for (i in 1:ntrees){
    tmp=delta[i,]%*%theta%*%phi[,dat$spp[i]]
    llk[i]=log(tmp)
  }
  sum(llk)
}