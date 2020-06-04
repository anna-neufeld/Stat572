genProbsBeta<- function(k,t, a,b) {
  choose(t,k)*beta(k+a, b+t-k)/beta(a,b)
}

### Parameter A is stuck at 1 in this example
## b is a guess for param b
## n is a vector n1 through n6
conditionalLik <- function(b,n) {
  ts <- 1:6
  probs <- sapply(ts, genProbsBeta, 6,1,b)
  p0 <- genProbsBeta(0,6,1,b)
  probs <- probs/(1-p0)
  return(-dmultinom(n, prob=probs, log=TRUE))
}


## Generate true data
## Where H is a beta(1,b) distrubution
genTrue <- function(b) {
  ts <- 0:6
  probs <- sapply(ts, genProbsBeta, 6,1,b)
  dat <- rmultinom(n=1, size=500, prob=probs)
  return(dat[2:7])
}

oneRep <- function(b) {
  observed <- genTrue(b)
  guess <- optim(par=1, fn=conditionalLik, n=observed, lower=0.00001, upper=1000,
                 method="L-BFGS-B")$par
  nu0 <- genProbsBeta(0,6,1,guess)
  Nhat <- sum(observed)/(1-nu0)
  return(Nhat)
}

oneRepALPHA <- function(b) {
  observed <- genTrue(b)
  guess <- optim(par=1, fn=conditionalLik, n=observed, lower=0.00001, upper=1000,
                 method="L-BFGS-B")$par
  nu0 <- genProbsBeta(0,6,1,guess)
  Nhat <- sum(observed)/(1-nu0)
  realAlpha <- 1-(1-0.1)^{1/6}
  prAlpha <- pbeta(realAlpha, 1, guess, lower.tail=FALSE)
  return( Nhat[1]*prAlpha)
}

bs <- c(seq(0.1,10,by=0.2), seq(11,50,by=1))
nits <- 1000
res <- matrix(0, nrow=length(bs), ncol=nits)
i=1
for (b in bs) {
  res[i,] <- replicate(nits, oneRep(b))
  i=i+1
}

res2 <- matrix(0, nrow=length(bs), ncol=nits)
i=1
for (b in bs) {
  res2[i,] <- replicate(nits, oneRepALPHA(b))
  i=i+1
}

### A FIGURE FOR MY PRELIM. A REAL ONE
sds <- apply(res, 1, sd)
sd2s <- apply(res2, 1, sd)
ggplot(data=NULL, aes(x=bs, y=sds, col="N")) + geom_point() + geom_smooth() +
  ylab("Empirical SD")+xlab("b")+ggtitle("Standard deviation of MLEs")+
  geom_point(data=NULL,aes(x=bs, y=sd2s, col="N_0.1")) +
  geom_smooth(data=NULL,aes(x=bs, y=sd2s, col="N_0.1"))+labs(col="Target Parameter")






### MY OWN EXPERIMENTS
conditionalLikFull <- function(param,n) {
  a <- param[1]
  b <- param[2]
  ts <- 1:6
  probs <- sapply(ts, genProbsBeta, 6,a,b)
  p0 <- genProbsBeta(0,6,a,b)
  probs <- probs/(1-p0)
  return(-dmultinom(n, prob=probs, log=TRUE))
}

genTrueFull <- function(a,b) {
  ts <- 0:6
  probs <- sapply(ts, genProbsBeta, 6,a,b)
  dat <- rmultinom(n=1, size=500, prob=probs)
  return(dat[2:7])
}

oneRepFull <- function(a,b) {
  observed <- genTrueFull(a,b)
  guess <- optim(par=c(1,1), fn=conditionalLikFull, n=observed, lower=c(0.00001,0.00001), upper=c(1000,1000),
                 method="L-BFGS-B")$par
  nu0 <- genProbsBeta(0,6,guess[1], guess[2])
  Nhat <- sum(observed)/(1-nu0)
  return(Nhat)
}



### FOR THE THREE BETA DISTRIBUTIONS FROM THEIR PAPER.
### Infinite MSE in the frequentist case. 
Nhats3 <- replicate(10000, oneRepFull(1/2, 1/2))
Nhats4 <- replicate(10000, try(oneRepFull(1, 10)))
Nhats4 <- as.numeric(Nhats4)
Nhats5 <- replicate(10000, oneRepFull(1,1))
mean((Nhats3-500)^2)
mean((Nhats4-500)^2, na.rm=TRUE)
mean((Nhats5-500)^2)

