rm(list=ls())
set.seed(3)

setwd('U:\\GIT_models\\SpatLDA\\fake data')
dat=read.csv('fake data.csv')

#create grid cells of 50 x 50 m
grid.size=50
brk=seq(from=0,to=1000,by=grid.size)
mid=brk[-length(brk)]+(grid.size/2)

tmp=cut(dat$x,breaks=brk); 
ind=as.numeric(tmp)
dat$xbin=mid[ind]

tmp=cut(dat$y,breaks=brk); 
ind=as.numeric(tmp)
dat$ybin=mid[ind]

uni.coord=unique(dat[,c('xbin','ybin')])
ngrid=nrow(uni.coord); ngrid

#get count of trees for each species in each grid cell
nspp=max(dat$spp)
res=matrix(0,ngrid,nspp)
for (i in 1:ngrid){
  cond=dat$xbin==uni.coord$xbin[i] &  
       dat$ybin==uni.coord$ybin[i]
  dat1=dat[cond,]
  tmp=table(dat1$spp)
  res[i,as.numeric(names(tmp))]=tmp
}
colnames(res)=paste0('spp',1:nspp)
res1=cbind(uni.coord,res)

#export results
setwd('U:\\GIT_models\\SpatLDA\\fake data')
write.csv(res1,'fake data1.csv',row.names=F)
