SpatLDA=function(nclust,ngibbs,nburn,coord.plot,dat,
                 gamma.theta,gamma.phi,
                 theta.init=NULL,phi.init=NULL){
  nplot=nrow(coord.plot)  
  
  #useful stuff
  ngrid=nrow(dat)
  ind=grep('spp',colnames(dat))
  nomes.spp=colnames(dat)[ind]
  nspp=length(nomes.spp)

  #get distances
  dist1=matrix(NA,ngrid,nplot)
  for (i in 1:nplot){
    x2=(dat$xbin-coord.plot$x[i])^2
    y2=(dat$ybin-coord.plot$y[i])^2
    dist1[,i]=sqrt(x2+y2)
  }

  #initial values for parameters
  cond=is.null(theta.init)
  if (cond)  theta=matrix(1/nclust,nplot,nclust)
  if (!cond) theta=theta.init
  cond=is.null(phi.init)
  if (cond)  phi=matrix(1/nspp,nclust,nspp)
  if (!cond) phi=phi.init
  sig2=100
  delta=get.delta(sig2=sig2,dist1=dist1)
  
  #initial values for cluster and plot membership
  array.gskp=array(0,dim=c(ngrid,nspp,nclust,nplot))
  nbins=nclust*nplot
  for (i in 1:ngrid){
    for (j in 1:nspp){
      tmp=dat[i,paste0('spp',j)]
      if (tmp>0){
        tmp1=rmultinom(1,size=tmp,prob=rep(1/nbins,nbins))
        array.gskp[i,j,,]=matrix(tmp1,nclust,nplot) 
      }
    }
  }
  # teste=apply(array.gskp,c(1,2),sum)
  # unique(teste-dat[,paste0('spp',1:nspp)])
  
  #MH stuff
  jump.sd=1 
  accept1=0
  nchange.jump=50
  
  #to store outcomes from gibbs sampler
  theta.out=matrix(NA,ngibbs,nclust*nplot)
  phi.out=matrix(NA,ngibbs,nclust*nspp)
  sig2.out=matrix(NA,ngibbs,1)
  llk=rep(NA,ngibbs)
  
  #run gibbs sampler
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)   
    
    #sample clust.id
    array.gskp=sample.clust.id(theta=theta,phi=phi,nplot=nplot,
                               array.gskp=array.gskp,
                               nspp=nspp,ngrid=ngrid,nclust=nclust)
    
    #sample plot.id
    tmp=sample.plot.id(theta=theta,delta=delta,array.gskp=array.gskp,
                       ngrid=ngrid,nclust=nclust,nplot=nplot,nspp=nspp)
    array.gskp=tmp$array.gskp
    soma.gk=tmp$SomaGK #useful to sample sig2
    
    #get useful summary for sampling theta and phi
    soma.skp=colSums(array.gskp,dims=1)
    
    #sample theta
    theta=sample.theta(soma.skp=soma.skp,nplot=nplot,nclust=nclust,
                       gamma1.theta=gamma.theta,nspp=nspp)
    # theta=theta.true
    
    #sample phi
    phi=sample.phi(soma.skp=soma.skp,gamma1.phi=gamma.phi,
                   nclust=nclust,nspp=nspp,nplot=nplot)
    # phi=phi.true
    
    #sample sig2
    tmp=sample.sig2(dist1=dist1,sig2=sig2,theta=theta,jump.sd=jump.sd,
                    ngrid=ngrid,nclust=nclust,soma.gk=soma.gk)
    accept1=accept1+tmp$accept1
    sig2=tmp$sig2
    # sig2=sig2.true
    delta=get.delta(sig2=sig2,dist1=dist1)
    
    #adapt MH
    if (i%%nchange.jump==0 & i<nburn){
      if (accept1/nchange.jump < 0.1) jump.sd=jump.sd*0.5
      if (accept1/nchange.jump > 0.6) jump.sd=jump.sd*2
      accept1=0
    } 
    
    #llk 
    llk[i]=get.llk(dat=dat,delta=delta,theta=theta,phi=phi,
                   ngrid=ngrid,nspp=nspp,nomes.spp=nomes.spp)
    
    #store results  
    theta.out[i,]=theta
    phi.out[i,]=phi
    sig2.out[i]=sig2
  }
  
  list(llk=llk,sig2=sig2.out,phi=phi.out,theta=theta.out)
}







