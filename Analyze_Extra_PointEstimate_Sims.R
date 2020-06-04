library(tidyverse)
library(reshape2)
setwd("~/Stat572")
res <- read.csv("cluster_RES.csv", sep=" ", header=FALSE)

names(res) <- c("scenario",  "mean0beta",  "mean005beta", "mean01beta", 
                     "mean05beta", "mean1beta", "median0beta", "median005beta","median01beta", 
                     "median05beta", "median1beta",
                     "mean0atom","mean005atom", 
                     "mean01atom","mean05atom","mean1atom","median0atom",
                     "median005atom", 
                     "median01atom","median05atom","median1atom")
res <- res %>% mutate(scenario = scenario%%5+1)

library(dplyr, warn.conflicts = FALSE)
res2 <- (res %>% group_by(scenario) %>% summarize(across(is.numeric, function(u) mean((u-500)^2))))[,c(1:6, 12:16)]
res2[1:2, 2:6] <- res2[1:2, 7:11] 
res2 <- res2[,1:6]
names(res2) <- c("scenario", "0", "0.005", "0.01", "0.05", "0.1")

names <- c("Atom 0.1,0.2,0.5", "Atom 0.01,0.1,0.6", "Beta(1/2, 1/2)", "Beta(1,10)", "Beta(1,1)")
trial <- reshape2::melt(res2, id=c('scenario')) %>% mutate(scenarioName = names[scenario])

ggplot(data=trial %>% filter(scenario > 2), aes(x=variable, y=sqrt(value), group = scenario, col=as.factor(scenarioName))) + 
  geom_point() + geom_line() + scale_y_continuous(trans="log10")+
  labs(col="Scenario") + xlab(expression(alpha)) + ylab(expression(sqrt(MSE)))+
  ggtitle("Average Root-MSE of Posterior Means")+  theme(text = element_text(size=24))
ggsave("MSE-plot-for-prez.png")

