setwd('~/Stat572/myOutput_513/')
library(tidyverse)
library(reshape2)
library(gridExtra)

### This stores the true N and then the actual alpha-observable Nalphas
### which change based on the data generating mechanism
N <- 500
alpha = c(0, 0.005,0.010, 0.050, 0.100)
props <- matrix(0, nrow=5, ncol=5)
#### NEED TO DERIVE AND REPRODUCE THESE. 
props[1,] <- c(1,1,1,1,1) ### BECAUSE OF THE STRICT GREATER THAN IN THE DEFINITION
props[2,] <- c(1, 1, 1, 0.7, 0.7) ### AGAIN BC OF THE STRICT GREATER THAN THIS MATTERS
props[3,] <- 1 - pbeta(1-(1-alpha)^(1/6), 1/2, 1/2)
props[4,] <- 1 - pbeta(1-(1-alpha)^(1/6), 1, 10)
props[5,] <- 1 - pbeta(1-(1-alpha)^(1/6), 1,1)
realNs <- as.data.frame(props*N)
names(realNs) = c("0", "0.005", "0.01", "0.05", "0.1")
realNs$sim <- 1:5
realNdat <- melt(realNs, id='sim')

#### PLOT 1 BETA
load("sim_beta_1.RData")
resBeta <- as.data.frame(resBeta[,1:5])
names(resBeta) = c("0", "0.005", "0.01", "0.05", "0.1")
resBeta$sim <- 1:NROW(resBeta)
resBeta2 <- melt(resBeta, id='sim')
p1 <- ggplot(resBeta2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  geom_point(data=realNdat %>% filter(sim==1), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 1, Beta Inference")+ xlab("")+ylab(expression(log[10](N[alpha])))+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),plot.title = element_text(size = 10))
#### PLOT 1 ATOM
load("sim_atom_1.RData")
resAtom <- as.data.frame(resAtom[,1:5])
names(resAtom) = c("0", "0.005", "0.01", "0.05", "0.1")
resAtom$sim <- 1:NROW(resAtom)
resAtom2 <- melt(resAtom, id='sim')
p2 <- ggplot(resAtom2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975)) +
  geom_point(data=realNdat %>% filter(sim==1), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 1, Atom Inference")+ xlab("")+ylab(expression(log[10](N[alpha])))+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),plot.title = element_text(size = 10))
#### PLOT 2 BETA
load("sim_beta_2.RData")
resBeta <- as.data.frame(resBeta[,1:5])
names(resBeta) = c("0", "0.005", "0.01", "0.05", "0.1")
resBeta$sim <- 1:NROW(resBeta)
resBeta2 <- melt(resBeta, id='sim')
p3 <- ggplot(resBeta2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  geom_point(data=realNdat %>% filter(sim==2), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 2, Beta Inference")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),plot.title = element_text(size = 10))
#### PLOT 2 ATOM
load("sim_atom_2.RData")
resAtom <- as.data.frame(resAtom[,1:5])
names(resAtom) = c("0", "0.005", "0.01", "0.05", "0.1")
resAtom$sim <- 1:NROW(resAtom)
resAtom2 <- melt(resAtom, id='sim')
p4 <- ggplot(resAtom2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975)) +
  geom_point(data=realNdat %>% filter(sim==2), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 2, Atom Inference")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),plot.title = element_text(size = 10))


#### PLOT 3 BETA
load("sim_beta_3.RData")
resBeta <- as.data.frame(resBeta[,1:5])
names(resBeta) = c("0", "0.005", "0.01", "0.05", "0.1")
resBeta$sim <- 1:NROW(resBeta)
resBeta2 <- melt(resBeta, id='sim')
p5 <- ggplot(resBeta2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  geom_point(data=realNdat %>% filter(sim==3), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 3, Beta Inference")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),plot.title = element_text(size = 10))


p52 <- ggplot(resBeta2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  geom_point(data=realNdat %>% filter(sim==3), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Beta(1/2, 1/2)")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=18),axis.text.x=element_text(angle=90),plot.title = element_text(size = 18))+ 
  xlab(expression(alpha))+ylab(expression(log[10](N[alpha])))
#### PLOT 3 ATOM
load("sim_atom_3.RData")
resAtom <- as.data.frame(resAtom[,1:5])
names(resAtom) = c("0", "0.005", "0.01", "0.05", "0.1")
resAtom$sim <- 1:NROW(resAtom)
resAtom2 <- melt(resAtom, id='sim')
p6 <- ggplot(resAtom2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975)) +
  geom_point(data=realNdat %>% filter(sim==3), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 3, Atom Inference") + xlab(expression(alpha))+ylab("")+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),plot.title = element_text(size = 10))


#### PLOT 4 BETA
load("sim_beta_4.RData")
resBeta <- as.data.frame(resBeta[,1:5])
names(resBeta) = c("0", "0.005", "0.01", "0.05", "0.1")
resBeta$sim <- 1:NROW(resBeta)
resBeta2 <- melt(resBeta, id='sim')
p7 <- ggplot(resBeta2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  geom_point(data=realNdat %>% filter(sim==4), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 4, Beta Inference")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),plot.title = element_text(size = 10))
p72 <- ggplot(resBeta2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  geom_point(data=realNdat %>% filter(sim==4), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Beta(1,10)")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=18),axis.text.x=element_text(angle=90),plot.title = element_text(size = 18))+
  xlab(expression(alpha))
#### PLOT 5 ATOM
load("sim_atom_4.RData")
resAtom <- as.data.frame(resAtom[,1:5])
names(resAtom) = c("0", "0.005", "0.01", "0.05", "0.1")
resAtom$sim <- 1:NROW(resAtom)
resAtom2 <- melt(resAtom, id='sim')
p8 <- ggplot(resAtom2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975)) +
  geom_point(data=realNdat %>% filter(sim==4), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 4, Atom Inference")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),plot.title = element_text(size = 10))


#### PLOT 5 BETA
load("sim_beta_5.RData")
resBeta <- as.data.frame(resBeta[,1:5])
names(resBeta) = c("0", "0.005", "0.01", "0.05", "0.1")
resBeta$sim <- 1:NROW(resBeta)
resBeta2 <- melt(resBeta, id='sim')
p9 <- ggplot(resBeta2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  geom_point(data=realNdat %>% filter(sim==5), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 5, Beta Inference")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),plot.title = element_text(size = 10))
p92 <- ggplot(resBeta2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  geom_point(data=realNdat %>% filter(sim==5), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Beta(1,1)")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=18),axis.text.x=element_text(angle=90),plot.title = element_text(size = 18))+
  xlab(expression(alpha))

#### PLOT 5 ATOM
load("sim_atom_5.RData")
resAtom <- as.data.frame(resAtom[,1:5])
names(resAtom) = c("0", "0.005", "0.01", "0.05", "0.1")
resAtom$sim <- 1:NROW(resAtom)
resAtom2 <- melt(resAtom, id='sim')
p10 <- ggplot(resAtom2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975)) +
  geom_point(data=realNdat %>% filter(sim==5), aes(x=variable, y=log10(value)), col="red")+
  ggtitle("Scenario 5, Atom Inference")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90), 
          plot.title = element_text(size = 10))

#### MAKE VIOLIN PLOT!!!!!!! FIGURE 4
png('MYVIOLIN.png',width=800,height=550)
grid.arrange(p1,p3,p5,p7,p9,p2,p4,p6,p8,p10, nrow=2)
dev.off()


png('MYVIOLIN_PARTIAL.png',width=800,height=450)
grid.arrange(p52,p72,p92, nrow=1)
dev.off()




### SNOWHOE HARE
load("hare_beta.RData")
load("hare_atom.RData")
resBeta <- as.data.frame(resBeta1[,1:5])
names(resBeta) = c("0", "0.005", "0.01", "0.05", "0.1")
resBeta$sim <- 1:NROW(resBeta)
resBeta2 <- reshape2::melt(resBeta, id='sim')
resAtom <- as.data.frame(resAtom1[,1:5])
names(resAtom) = c("0", "0.005", "0.01", "0.05", "0.1")
resAtom$sim <- 1:NROW(resAtom)
resAtom2 <- reshape2::melt(resAtom, id='sim')
p1 <- ggplot(resBeta2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  ggtitle("Snowshoe Hare, Beta Inference")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=18),axis.text.x=element_text(angle=90),plot.title = element_text(size = 18))+
  xlab(expression(alpha))
p2 <- ggplot(resAtom2,aes(x=variable,y=log10(value))) + 
  geom_violin(draw_quantiles = c(0.025,0.975))+
  ggtitle("Snowshow Hare, Atom Inference")+ xlab("")+ylab("")+ 
  theme(text=element_text(size=18),axis.text.x=element_text(angle=90),plot.title = element_text(size = 18))+
  xlab(expression(alpha))

png('snowshoeHARE.png',width=800,height=450)
grid.arrange(p1,p2, nrow=1)
dev.off()



### MAKE FIGURE 1 THE THREE HISTOGRAMS
### PANEL A IS SNOWSHOE HARE 
### PANELS BC ARE BETA: SIM 4?? I KNOW THIS FROM LOOKING AT THEIR CODE 
#load("sim_beta_4.RData")
#ggplot(data=NULL, aes(x=dat[['Nhat']])) + geom_histogram() + xlab("N_alpha")+geom_vline(xintercept=500, col="red")+
 # geom_vline(xintercept = quantile(dat[['Nhat']], 0.025))+
#  geom_vline(xintercept = quantile(dat[['Nhat']], 0.975))

#ggplot(data=NULL, aes(x=resBeta[,3])) + geom_histogram() + xlab("N_alpha")+geom_vline(xintercept=500, col="red")+
#  geom_vline(xintercept = quantile(resBeta[,3], 0.025))+
#  geom_vline(xintercept = quantile(resBeta[,3], 0.975))

#ggplot(data=NULL, aes(x=resBeta[,1])) + geom_histogram(bins=50) + xlab("N")+scale_x_log10()+geom_vline(xintercept=500, col="red")+
#    geom_vline(xintercept = quantile(resBeta[,1], 0.025))+
#    geom_vline(xintercept = quantile(resBeta[,1], 0.975))

























