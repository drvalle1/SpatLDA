rm(list=ls())
library('ggplot2')
library(gridExtra)
library(grid)

setwd('U:\\GIT_models\\SpatLDA\\fake data')
true.dat=read.csv('fake data proportions.csv')

setwd('U:\\GIT_models\\SpatLDA\\results')
lda.spat=read.csv('estim proportions.csv')

# lda.abund=read.csv('results LdaAbund\\estim proportions.csv')
# lda.agg=read.csv('results agg\\estim proportions.csv')

#plot results
p1t=ggplot() +
  geom_tile(data = true.dat, alpha = 0.8,aes(x = x, y = y,fill = gr1)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) +
  theme(legend.position = "none") 
p2t=ggplot() +
  geom_tile(data = true.dat, alpha = 0.8,aes(x = x, y = y,fill = gr2)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) +
  theme(legend.position = "none") 
p3t=ggplot() +
  geom_tile(data = true.dat, alpha = 0.8,aes(x = x, y = y,fill = gr3)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) +
  theme(legend.position = "none") 

p1.spat=ggplot() +
  geom_tile(data = lda.spat, alpha = 0.8,aes(x = x, y = y,fill = c1)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) +
  theme(legend.position = "none") 
p2.spat=ggplot() +
  geom_tile(data = lda.spat, alpha = 0.8,aes(x = x, y = y,fill = c2)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) +
  theme(legend.position = "none") 
p3.spat=ggplot() +
  geom_tile(data = lda.spat, alpha = 0.8,aes(x = x, y = y,fill = c3)) +
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',limits=c(0,1),midpoint=0.5) +
  theme(legend.position = "none") 

grid.arrange(p1t, p2t, p3t,
             p1.spat,p2.spat,p3.spat,nrow = 2)
