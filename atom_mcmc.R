# n is the data with a guess made up for 0 observations
# d is number of atoms
# S is number of time points of observed data
# alpha is a vector of the alpha values you want to predict for
# s is parameter for prior
# nmc is number of Monte Carlo iterations
atommix_me <- function(n,d,alpha,s,nmc) {
  
  T = length(n)-1
  
  ## To store the results
  MU = matrix(0, nrow=nmc, ncol=d)
  weights = matrix(0, nrow=nmc, ncol=d)
  nalpha = length(alpha)
  Nhat= matrix(0,nrow=nmc,ncol=nalpha+1)
  
  
  ## Initial guesses for atom locations and weights
  mu = runif(d) 
  eta = rep(1/d, d)

  ## Do the Monte Carlo iterations
  for (t in 1:nmc) {
    for (j in 1:d) {
      ## Relecting Uniform Random Walk for each Atom Location
      prop_muj = runif(n=1, mu[j]-s[j], mu[j]+s[j])
      if (prop_muj < 0) { prop_muj = -prop_muj
      } else {
        if (prop_muj > 1) {prop_muj = 1-(prop_muj-1)}
      }
      
      ## Prop_mu is the old mu with just this one location replaced with the proposal
      prop_mu = mu
      prop_mu[j] = prop_muj
      
      ### Compute our acceptance probability!! This is the log, so
      ### its a difference not a ratio. 
      ### Priors all uniform- so we just use a ratioof likelihoods. 
      prop_ll = lik_atom(prop_mu,eta,n,TRUE)
      curr_ll = lik_atom(mu,eta,n,TRUE)           
      lrr = prop_ll-curr_ll
      acc = runif(1) < exp(lrr)
      if (acc) {mu[j] = prop_mu[j]}           
    }
    
    ## Now we need to do the weights
    lr = matrix(0, T+1, d)
    for (j in 1:d) {
      
      ### What if all the weight was at this atom?
      ## What would the likelihood be in that case? 
      wts0 = rep(0, d)
      wts0[j]=1
      
      #### THIS JUST GETS NU_J(H^**)
      pk = lik_atom(mu, wts0, n, FALSE)
      lr[,j] = log(pk)
    }
    #### DO THESE THINGS WITH LR
    lr = t(apply(lr, 1, function(u) log(eta) + u))
    lr = t(apply(lr, 1, function(u) u - max(u)))
    lr = exp(lr)
    lcprob = lr * 1/rowSums(lr)
    
    z <- matrix(0, nrow=T+1, ncol=d)
    for (j in 1:(T+1)) {
      z[j,] <- rmultinom(1, n[j], lcprob[j,])
    }
    
    ### Are these guesses that of the people we saw, we had
    ### fullZ[1] with p_i = u1, fullZ[2] with p_i = u2, etc??
    fullZ <- colSums(z)
    
    #### THIS SHOULD BE DIRICHLET. BUT IT TURNS OUT GAMMA(a,1) normalized to sum to 1 is a 
    ### DIRICHLET YAY
    tmp <- rgamma(3, 1/d+fullZ, 1)
    eta = tmp/sum(tmp)
    
    pr = lik_atom(mu, eta, n,FALSE)
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
    
    weights[t,] = eta
    MU[t,] = mu
  }
  
  res <- cbind(Nhat, weights, MU)
  colnames(res) <- c(paste("Nhat", 0:length(alpha)), paste("eta", 1:d), paste("mu", 1:d))

  return(res)
}


lik_atom <- function(mu, wts, n,lik = TRUE) {
  pk = rep(0, length(n))
  t <- length(n)-1
  for (k in 0:t) {
    pk[k+1] = choose(t,k)*wts%*%(mu^k*(1-mu)^(t-k))
  }
  ll = sum(n*log(pk))
  if (lik) {return(ll)
  } else {return(pk)}
}





