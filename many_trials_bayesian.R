require(snow)
require(doParallel)
require(foreach)
require(iterators)
#require(MASS)

### Packages I will need and global vars
setwd("~/Prelim/capture-het-pub/My_Sims")
source("generate_Data.R")
source("atommix_me.R")		
source("betamix_me.R")
source("supplementary_functions.R")


##Cluster set up stuff
cl <- makeCluster(4) #Number of cores you want 
registerDoParallel(cl)
clusterSetupRNG(cl,seed=1)

set.seed(555)
nIts <- 1000

### More global params
ts <- 0:6
NU <- rbind(sapply(ts, genProbs1, 6), sapply(ts, genProbs2, 6),sapply(ts, genProbsBeta, 6,1/2,1/2),
            sapply(ts, genProbsBeta, 6,1,10), sapply(ts, genProbsBeta, 6,1,1))


###### THIS IS THE IMPORTANT THING: EXPORT VARS
clusterExport(cl, c("betamix_me", "atommix_me", "NU", "loglik_atom",
                    "loglik_beta"))


##### NOW RUN
out <- foreach(j = c(1:(nIts*5))) %dopar% {
  write(paste("Starting ", j, "th job.\n",sep=''),file='extra_sim_log.txt',append=TRUE)
  N = 500
  t = 6
  param_set = j%%5+1
  alpha = c(0.005, 0.01, 0.05, 0.1)
  real_alpha = 1-(1-alpha)^(1/t)
  nmc = 100000
  a0 = 1/4
  b0 = 1/4
  sb = 0.2 
  sa = 0.2 
  d=3
  s = c(0.05, 0.05, 0.05)
  n <- rmultinom(1,N,NU[param_set,])
  
  ## should i replace n1 w a random guess???
  resBeta <- betamix_me(n, real_alpha, sb, sa, a0, b0, nmc)
  resAtom <- atommix_me(n, d, real_alpha, s, nmc)
  
  betaMeans <- apply(resBeta[,1:5], 2, mean)
  betaMedian <- apply(resBeta[,1:5], 2, median)
  atomMeans <- apply(resAtom[,1:5], 2, mean)
  atomMedian <- apply(resAtom[,1:5], 2, median)
  
  results <- c(j, betaMeans, betaMedian , atomMeans,atomMedian)
  write(paste(results, collapse=" "), file='cluster_RES.csv',append=TRUE)
  results
}

new_out = as.data.frame(matrix(unlist(out), ncol=21 , byrow=TRUE))
names(new_out)=names(out[[1]])
L=6
save(new_out, file = paste("ALL_RES",L,".RData",sep=''))
stopCluster(cl)
