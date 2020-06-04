betamix_me <- function(n, alpha, sa, sb, a0, b0, nmc) {
  ## Empty vecs
  A = rep(0, nmc)
  B = rep(0, nmc)
  acc = rep(0, nmc)

  ## Initial starting guesses
  a = runif(1)
  b = runif(1)
  
  nalpha = length(alpha)
  Nhat = matrix(0, nrow=nmc, ncol=nalpha+1)
  
  for (t in 1:nmc) {
    ## a, b are the old guesses for the parameters
    ## This is the lognormal random walk 
    proposed_ab = exp(rnorm(2, mean=c(log(a), log(b)), sd = c(sa, sb)))
    prop_a = proposed_ab[1]
    prop_b = proposed_ab[2]
    
    ### For the likelihood ratio. 
    prop_ll = loglik_beta(prop_a,prop_b, n)
    curr_ll = loglik_beta(a,b, n)
    
    ## Ratio of priors for gamma(a0, b0) prior
    log_pr = (a0-1)*log(prop_a/a)-b0*(prop_a-a) + (a0-1)*log(prop_b/b)-b0*(prop_b-b)
    
    #### This is just ratio of proposed parameters to old parameters
    #### Turns out to be q ratio for the logNoral random walk
    log_q = log(prop_a*prop_b/(a*b))
    
    ### Log lik ratio - a difference
    lrr = prop_ll-curr_ll+log_pr+log_q
    
    ### Decide if we accept
    acc[t] = runif(1) < exp(lrr)
    
    ### Update if accept.  
    if(acc[t]) {
      a = prop_a
      b = prop_b
      A[t] = a
      B[t] = b
    }
    
    ## BASED ON CURRENT a,b PARAMETERS, get probabilities of bins
    pr = get_nus_beta(a,b, n)
    rho = 1 - pr[1]
    
    ### Update guess for number of unobserved people- draw it from the neg binomial
    ### rnegbinom sometimes literally gives infinity. That's why try statement. 
    n1 <- try(rnbinom(n=1, size= sum(n[2:length(n)]), p = rho))
    ctr=0
    while(is.na(n1)) {
      n1 <- try(rnbinom(n=1, size= sum(n[2:length(n)]), p = rho))
      if (ctr > 10) {n1 = 1e8}
    } 
    n[1] <- n1
    Nhat[t,1] <- sum(n) 
    
    ### Draw all the NAlpha vals
    for (l in 1:nalpha) {
      prca = pbeta(alpha[l], a, b, lower.tail=FALSE)
      Nhat[t,l+1] = rbinom(n=1, size=sum(n), p=prca)
    }
  }
  
  res <- cbind(Nhat, a, b)
  colnames(res) <- c(paste("Nhat", 0:length(alpha)), "a", "b")
  
  return(res)
}


loglik_beta <- function(a,b, n) {
  pk = rep(0, length(n))
  t <- length(n)-1
  for (k in 0:t) {
    pk[k+1] = choose(t,k)*beta(k+a, b+t-k)/beta(a,b)
  }
  ll = sum(n*log(pk)) ### for a ratio, only this part of multinomial likelihood matters
  return(ll)
}

get_nus_beta <- function(a,b,n) {
  nuk = rep(0, length(n))
  t <- length(n)-1
  for (k in 0:t) {
    nuk[k+1] = choose(t,k)*beta(k+a, b+t-k)/beta(a,b)
  }
  return(nuk)
}


