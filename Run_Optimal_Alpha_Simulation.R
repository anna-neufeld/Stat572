#### Generates the nus. 
#### when H is Beta(a,b)
genProbsBeta<- function(k,t, a,b) {
  choose(t,k)*beta(k+a, b+t-k)/beta(a,b)
}

## Theta = (a,b) for the Beta distribution
## n is vector n1 through nT of observed data
conditionalLik <- function(theta,n) {
  a <- theta[1]
  b <- theta[2]
  ts <- 1:6
  probs <- sapply(ts, genProbsBeta, 6,a,b)
  p0 <- genProbsBeta(0,6,a,b)
  probs <- probs/(1-p0)
  ## returns negative likelihood for the optim function
  return(-dmultinom(n, prob=probs, log=TRUE))
}


## Generate a realization of true data
## Where H is a beta(a,b) distrubution
genTrue <- function(a,b) {
  ts <- 0:6
  probs <- sapply(ts, genProbsBeta, 6,a,b)
  dat <- rmultinom(n=1, size=500, prob=probs)
  return(dat[2:7])
}

## One repetitiion in my study involves generating observed data, finding the (a,b)
## that maximimzes the conditional likelihood, computing a Horvitz Thomspon estimator Nhat,
## then computing Nhat_alpha for all the desired values of alpha
oneRep <- function(a,b, alphas) {
  try1 <- try({
  Nhat <- rep(0, length(alphas)+1)
  observed <- genTrue(a,b)
  guess <- optim(par=c(1,1), fn=conditionalLik, n=observed, lower=0.00001, upper=1000,
                 method="L-BFGS-B")$par
  nu0 <- genProbsBeta(0,6,guess[1],guess[2])
  Nhat[1] <- sum(observed)/(1-nu0)
  
  ctr=2
  for (alpha in alphas) {
    realAlpha <- 1-(1-alpha)^{1/6}
    prAlpha <- pbeta(realAlpha, guess[1], guess[2], lower.tail=FALSE)
    Nhat[ctr] <- Nhat[1]*prAlpha ## expected value of the binomial
    ctr=ctr+1
  }}, silent=TRUE)
  if (class(try1)=="try-error") {
    print("hello")
    return(rep(NA, length(alpha)+1))
  } else {return(Nhat)}
}


#### Run an experiment to make Figure for results section
## The figure about optimal values of alpha

params <- rbind(c(1/2, 1/2), c(1,10), c(1,1)) ## use Johndrow et al's same scenarios
alphas <- seq(0.0001, 0.25, length.out=50) ## Try out a LOT of alphas!!!
nits <- 10000
res <- matrix(0, nrow=nits*NROW(params), ncol=length(alphas)+3)
res[,1] <- rep(params[,1], each=nits)
res[,2] <- rep(params[,2], each=nits)
colnames(res) <- c("a", "b", "N", round(alphas,3))
j = 1
for (i in 1:NROW(params)) {
  set.seed(i)
  thisres <- t(replicate(nits, oneRep(params[i,1], params[i,2], alphas)))
  res[j:(j+nits-1),-c(1,2)] <- thisres
  j <- j+nits
}

res <- as.data.frame(res)

### This is how I find out the theoretcal bias
### Very closely matches the empirical bias, so that's good. 
realAlphas <- c(0,1-(1-alphas)^{1/6})
prAlphas1 <- pbeta(realAlphas, 1/2, 1/2)
prAlphas2 <- pbeta(realAlphas, 1, 10)
prAlphas3 <- pbeta(realAlphas, 1, 1)


### Compute bias and variance for each alpha
summaryTabMean <- res %>% group_by(a,b) %>% summarize(across(is.numeric, mean))
summaryTabVar <- res %>% group_by(a,b) %>% summarize(across(is.numeric, var))


plot1 <- cbind(c(0, alphas), as.numeric(summaryTabMean[1, -c(1,2)]),as.numeric(summaryTabVar[1, -c(1,2)]),
               (prAlphas1*500)^2)
plot1 <- as.data.frame(plot1)
names(plot1) <- c("alpha", "mean","var", "biassq")
p1 <- ggplot(data=plot1, aes(x=alpha, y = biassq, col="Bias^2")) + geom_line() + ## bias squared
  ylim(0,5000)+ ylab(" ") + ggtitle("Scenario 1: Beta(1/2, 1/2)") + xlab(expression(alpha))+ 
  geom_line(aes(x=alpha, y=var, col="Variance")) + ##variance
  geom_line(aes(x=alpha, y=biassq + var, col="MSE"))+xlim(0,0.25)+guides(col=FALSE)+ ##risk
  geom_vline(xintercept = plot1$alpha[which.min(plot1$biassq + plot1$var)])



plot2 <- cbind(c(0, alphas), as.numeric(summaryTabMean[3, -c(1,2)]),as.numeric(summaryTabVar[3, -c(1,2)]),(prAlphas2*500)^2)
plot2 <- as.data.frame(plot2)
names(plot2) <- c("alpha", "mean","var", "biassq")
p2 <- ggplot(data=plot2, aes(x=alpha, y = biassq, col="Bias^2")) + geom_line() +
  ylim(0,25000)+
  geom_line(aes(x=alpha, y=var, col="Variance")) + ylab(" ") + ggtitle("Scenario 2: Beta(1, 10)") +
  xlab(expression(alpha))+ 
  geom_line(aes(x=alpha, y=biassq + var, col="MSE"))+xlim(0,0.25)+guides(col=FALSE)+
  geom_vline(xintercept = plot2$alpha[which.min(plot2$biassq + plot2$var)])


plot3 <- cbind(c(0, alphas), as.numeric(summaryTabMean[2, -c(1,2)]),as.numeric(summaryTabVar[2, -c(1,2)]),
               (prAlphas3*500)^2)
plot3 <- as.data.frame(plot3)
names(plot3) <- c("alpha", "mean","var", "biassq")
p3 <- ggplot(data=plot3, aes(x=alpha, y = biassq, col="Bias^2")) + geom_line() + ylim(0,1000)+
  geom_line(aes(x=alpha, y=var, col="Variance")) + ylab(" ") + ggtitle("Scenario 3: Beta(1, 1)") +
  xlab(expression(alpha))+ 
  geom_line(aes(x=alpha, y= biassq + var, col="MSE"))+xlim(0,0.25)+
  geom_vline(xintercept=plot3$alpha[which.min(plot3$biassq + plot3$var)])

  #geom_line(aes(x=alpha, y=(mean-500)^2, col="empiricalbias"))

png("~/Stat572/alpha_opt.png", width=700, height=400)
gridExtra::grid.arrange(p1,p2,p3, nrow=1, widths=c(2,2,2.5))
dev.off()







