betamix_me <- function(n, alpha, sa, sb, a0, b0, nmc) {
  A = rep(0, nmc)
  B = rep(0, nmc)
  acc = rep(0, nmc)

  ## Initial starting guesses? 
  a = runif(1)
  b = runif(1)
  
  nalpha = length(alpha)
  
  Nhat = matrix(0, nrow=nmc, ncol=nalpha+1)
  
  for (t in 1:nmc) {
    ### I AM NOT LOOPING OVER J
    ### SO J IS NOT NEEDED ANYWHERE
    
    ## a, b are the old guesses for the parameters
    ## sa and sb are parameters that we passed in. 
    proposed_ab = exp(rnorm(2, mean=c(log(a), log(b)), sd = c(sa, sb)))
    prop_a = proposed_ab[1]
    prop_b = proposed_ab[2]
    
    prop_ll = loglik_beta(prop_a,prop_b, n, TRUE)
    curr_ll = loglik_beta(a,b, n, TRUE)
    
    ## DO THE MATH TO FIGURE OUT WHY
    ## THIS COMPUTES A RATIO OF PRIORS. 
    log_pr = (a0-1)*log(prop_a/a)-b0*(prop_a-a) + (a0-1)*log(prop_b/b)-b0*(prop_b-b)
    
    log_q = log(prop_a/a)+log(prop_b/b)
    
    lrr = prop_ll-curr_ll+log_pr+log_q
    
    acc[t] = runif(1) < exp(lrr)
    
    ### If we get stuck somewhere where we never change we are in trouble!!!
    ### I saw that happening!!!
    ### How to avoid?? 
    if(acc[t]) {
      a = prop_a
      b = prop_b
      A[t] = a
      B[t] = b
    }
    
    ## BASED ON CURRENT a,b PARAMETERS, get probabilities of bins
    pr = loglik_beta(a,b, n, FALSE)
    rho = 1 - pr[1]
    
    ### UPDATE GUESS FOR NUMBER OF UNOBSERVED PEOPLE! Draw it from the neg binom
    ### OMG I JUST READ ABOUT NEGATIVE BINOMIAL
    ### AND ALL THE SUDDEN IT MAKES SENSE
    #n[1] <- rnbinom(n=1, size= sum(n[2:length(n)]), p = rho)
    n1 <- try(rnbinom(n=1, size= sum(n[2:length(n)]), p = rho))
    ctr=0
    while(is.na(n1)) {
      n1 <- try(rnbinom(n=1, size= sum(n[2:length(n)]), p = rho))
      if (ctr > 10) {n1 = 1e8}
    } 
    n[1] <- n1
    Nhat[t,1] <- sum(n) ## YAY
    
    for (l in 1:nalpha) {
      prca = pbeta(alpha[l], a, b, lower.tail=FALSE)
      Nhat[t,l+1] = rbinom(n=1, size=sum(n), p=prca)
    }
  }
  
  res <- cbind(Nhat, a, b)
  colnames(res) <- c(paste("Nhat", 0:length(alpha)), "a", "b")
  
  return(res)
}


loglik_beta <- function(a,b, n, lik = TRUE) {
  pk = rep(0, length(n))
  t <- length(n)-1
  for (k in 0:t) {
    partial = beta(k+a, b+t-k)/beta(a,b)
    pk[k+1] = choose(t,k)*partial
  }
  ll = sum(n*log(pk))
  if (lik) {return(ll)
  } else {return(pk)}
}


