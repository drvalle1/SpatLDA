rm(list=ls())
setwd('U:\\independent studies\\SpatLDA\\sim1')
true.dat=read.csv('fake data\\fake data proportions.csv')
ind=which(colnames(true.dat)%in%paste0('gr',1:3))
colnames(true.dat)[ind]=c('c1.true','c2.true','c3.true')

lda.spat=read.csv('results\\estim proportions.csv')
ind=which(colnames(lda.spat)%in%paste0('c',1:3))
colnames(lda.spat)[ind]=paste0('c',1:3,'.spat')

lda.abund=read.csv('results LdaAbund\\estim proportions.csv')
ind=which(colnames(lda.abund)%in%paste0('c',1:3))
colnames(lda.abund)[ind]=paste0('c',1:3,'.abund')

# lda.agg=read.csv('results agg\\estim proportions.csv')
# ind=which(colnames(lda.agg)%in%paste0('c',1:3))
# colnames(lda.agg)[ind]=paste0('c',1:3,'.agg')

fim=merge(true.dat,merge(lda.spat,lda.abund,all=T),all=T); dim(true.dat); dim(fim)

mae.spat=c(mean(abs(fim$c1.true-fim$c1.spat)),
           mean(abs(fim$c2.true-fim$c2.spat)),
           mean(abs(fim$c3.true-fim$c3.spat)))
mae.abund=c(mean(abs(fim$c1.true-fim$c1.abund)),
            mean(abs(fim$c2.true-fim$c2.abund)),
            mean(abs(fim$c3.true-fim$c3.abund)))
# mae.agg=c(mean(abs(fim$c1.true-fim$c1.agg)),
#           mean(abs(fim$c2.true-fim$c2.agg)),
#           mean(abs(fim$c3.true-fim$c3.agg)))

mae.spat
mae.abund
# mae.agg

sum(mae.spat)
sum(mae.abund)
# sum(mae.agg)

par(mfrow=c(1,2))
plot(fim$c1.true,fim$c1.spat)
plot(fim$c1.true,fim$c1.abund)