# n is the data with a guess made up for 0 observations
# d is number of atoms
# S is number of time points ob observed data
# alpha is a vector of the alpha values you want to predict for
# s is parameter for prior
atommix_me <- function(n,d,alpha,s,nmc) {
  
  S = length(n)-1
  MU = matrix(0, nrow=nmc, ncol=d)
  acc = matrix(0, nrow=nmc, ncol=d)
  mu = runif(d) 
  eta = rep(1/d, d)
  ETA = matrix(0, nrow=nmc, ncol=d)

  nalpha = length(alpha)
  Nhat= matrix(0,nrow=nmc,ncol=nalpha+1)
  
  for (t in 1:nmc) {
    for (j in 1:d) {
      ## s is a parameter. 
      prop_muj = runif(n=1, mu[j]-s[j], mu[j]+s[j])
      if (prop_muj < 0) { prop_muj = -prop_muj
      } 
      if (prop_muj > 1) {prop_muj = 1-(prop_muj-1)}
      prop_mu = mu
      prop_mu[j] = prop_muj
      
      ### THIS SEEMS WRONG WHY IS IT 7 BY 7
      prop_ll = loglik_atom(prop_mu,eta,n,TRUE)
      curr_ll = loglik_atom(mu,eta,n,TRUE)           
      lrr = prop_ll-curr_ll
      acc[t,j] = runif(1) < exp(lrr)
      if (acc[t,j]) {mu[j] = prop_mu[j]}           
    }
    lr = matrix(0, S+1, d)
    for (j in 1:d) {
      wts0 = rep(0, d)
      wts0[j]=1
      
      #### THIS JUST GETS NU_J(H^**)
      pk = loglik_atom(mu, wts0, n, FALSE)
      lr[,j] = log(pk)
    }
    #### DO THESE THINGS WITH LR
    lr = t(apply(lr, 1, function(u) log(eta) + u))
    lr = t(apply(lr, 1, function(u) u - max(u)))
    lr = exp(lr)
    lcprob = lr * 1/rowSums(lr)
    
    z <- matrix(0, nrow=S+1, ncol=d)
    for (j in 1:(S+1)) {
      z[j,] <- rmultinom(1, n[j], lcprob[j,])
    }
    
    ### Are these guesses that of the people we saw, we had
    ### fullZ[1] with p_i = u1, fullZ[2] with p_i = u2, etc??
    fullZ <- colSums(z)
    
    #### THIS SHOULD BE DIRICHLET. BUT IT TURNS OUT GAMMA(a,1) normalized to sum to 1 is a 
    ### DIRICHLET YAY
    tmp <- rgamma(3, 1/d+fullZ, 1)
    eta = tmp/sum(tmp)
    
    pr = loglik_atom(mu, eta, n,FALSE)
    rho=1-pr[1]
    n[1] <- rnbinom(1,sum(n)-n[1], rho)
    Nhat[t,1] = sum(n)
    
    
    for (l in 1:nalpha) {
        wts = eta/sum(eta)     
        prca = sum(wts[mu>alpha[l]])
        prca[prca<0] = 0
        prca[prca>1] = 1
        Nhat[t,l+1] = rbinom(1,sum(n),prca)
    }
    
    ETA[t,] = eta
    MU[t,] = mu
  }
  
  res <- cbind(Nhat, ETA, MU)
  colnames(res) <- c(paste("Nhat", 0:length(alpha)), paste("eta", 1:d), paste("mu", 1:d))

  return(res)
}


loglik_atom <- function(mu, wts, n,lik = TRUE) {
  pk = rep(0, length(n))
  t <- length(n)-1
  for (k in 0:t) {
    pk[k+1] = choose(t,k)*wts%*%(mu^k*(1-mu)^(t-k))
  }
  ll = sum(n*log(pk))
  if (lik) {return(ll)
  } else {return(pk)}
}





