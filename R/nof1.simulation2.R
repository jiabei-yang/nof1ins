
nof1.ordinal.simulation2 <- function(alpha = 0, beta_B = 1, cut = c(-2,-1.5,-1,-0.5,0, 0.5, 1, 1.5, 2, 2.1, 2.3), ncat = 11){
  
  inv_logit <- function(a){
    1/(1+exp(-a))
  }

  Treat <-  rep(c("A", "B", "A", "B", "A", "B"), each = 3)
  nobs <- length(Treat)
  
  Y <- epsilon <- mu <- rep(NA, nobs)
  
  q <- matrix(0, nrow = length(Treat), ncol = ncat - 1)
  p <- matrix(0, nrow = nobs, ncol = ncat)
  
  for(i in 1:nobs){
    
    mu[i] <- alpha
    if(Treat[i] == "B"){
      mu[i] <- mu[i] + beta_B
    } 
    
    for(r in 1:(ncat-1)){
      q[i,r] <- inv_logit(mu[i] - cut[r])
    }
    
    p[i,1] <- 1 - q[i,1]
    for(r in 2:(ncat-1)){
      p[i,r] <- q[i,r-1] - q[i,r]
    }
    p[i,ncat] <- q[i,(ncat-1)]
    Y[i] <- sample(1:ncat, size = 1,  prob = p[i,])
  }  
  list(Y = Y, Treat = Treat)
}
