####### ALL SIMULATION RESULTS
setwd("~/Prelim/myCode")
source("atom_mcmc.R")		
source("beta_mcmc.R")
source("generate_nus.R")



### Reproduce the main violin plot from their paper. 5 scenarios, each with two fitting mechanisms..  
ts <- 0:6
NU <- rbind(sapply(ts, genProbs1, 6), sapply(ts, genProbs2, 6),sapply(ts, genProbsBeta, 6,1/2,1/2),
            sapply(ts, genProbsBeta, 6,1,10), sapply(ts, genProbsBeta, 6,1,1))
N = 500
t = 6
alpha = c(0.005, 0.01, 0.05, 0.1)
real_alpha = 1-(1-alpha)^(1/t)
nmc = 1000000
a0 = 1/4
b0 = 1/4
sb = 0.2 
sa = 0.2 
s = c(0.05, 0.05, 0.05)
set.seed(111)
for (j in 1:NROW(NU)) {
  n <- rmultinom(1,N,NU[j,]) ## Data realization!!!
  resBeta <- betamix_me(n, real_alpha, sb, sa, a0, b0, nmc)
  resAtom <- atommix_me(n, 3, real_alpha, s, nmc)
  save(resBeta, file = paste0("myOutput_513/sim_beta_", j, ".RData"))
  save(resAtom, file = paste0("myOutput_513/sim_atom_", j, ".RData"))
}


##### ALL HARE RESULTS
require(Rcapture)
data(hare)
n <- table(apply(hare,1,sum))
### I tried out different guesses for n0 to make sure that things converge
guess1 = 25 ### same as ppl observed once
nfull1 <- c(guess1, n)

nmc = 1000000
t = length(n)
a0 = 1/4
b0 = 1/4
sb = 0.2 
sa = 0.2 
alpha = c(0.005, 0.01, 0.05, 0.1)
real_alpha = 1-(1-alpha)^(1/t)
resBeta1 <- betamix_me(nfull1, real_alpha, sb, sa, a0, b0, nmc)
save(resBeta1, file = "myOutput_513/hare_beta.RData")


s = c(0.05, 0.05, 0.05)
resAtom1 <- atommix_me(nfull1, 3, real_alpha, s, nmc)
save(resAtom1, file = "myOutput_513/hare_atom.RData")
