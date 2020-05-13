### For each of 5 scenarios, I need to compute v_K(H).
## With T=6. 

### SCENARIO 1.
genProbs1 <- function(k,T) {
  qs <- c(0.1, 0.2, 0.5)
  dQs <- c(0.3, 0.2, 0.5)
  choose(T,k)*sum(qs^k*(1-qs)^(T-k)*dQs)
}

### SCENARIO 2.
genProbs2 <- function(k,T) {
  qs <- c(0.01, 0.1, 0.6)
  dQs <- c(0.3, 0.2, 0.5)
  choose(T,k)*sum(qs^k*(1-qs)^(T-k)*dQs)
}

### SCENARIO 3,4,5: Beta(1/2, 1/2)
genProbsBeta<- function(k, T, a,b) {
  choose(T,k)*beta(k+a, b+T-k)/beta(a,b)
}



