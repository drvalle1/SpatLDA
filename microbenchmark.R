library('microbenchmark')
microbenchmark(
  sample.clust.id(theta=theta,phi=phi,nplot=nplot,
                             array.gskp=array.gskp,
                             nspp=nspp,ngrid=ngrid,nclust=nclust),
  
  #sample plot.id
  sample.plot.id(theta=theta,delta=delta,array.gskp=array.gskp,
                            ngrid=ngrid,nclust=nclust,nplot=nplot,nspp=nspp),
  
  apply(array.gskp,c(2,3,4),sum),
  
  #sample theta
  sample.theta(soma.skp=soma.skp,nplot=nplot,nclust=nclust,gamma1.theta=gamma.theta),
  # theta=theta.true
  
  #sample phi
  sample.phi(soma.skp=soma.skp,gamma1.phi=gamma.phi,nclust=nclust,nspp=nspp),
  # phi=phi.true
  
  #sample sig2
  sample.sig2(dist1=dist1,sig2=sig2,theta=theta,jump.sd=jump.sd,
                  ngrid=ngrid,nclust=nclust,soma.gk=soma.gk),
  
  get.delta(sig2=sig2,dist1=dist1),
  
  #llk 
  get.llk(dat=dat,delta=delta,theta=theta,phi=phi,
                 ngrid=ngrid,nspp=nspp,nomes.spp=nomes.spp)
)
