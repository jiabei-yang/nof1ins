nof1.prior.default <- function(prior.param){

  with(prior.param, {

  if(is.null(alpha.prior)){
    alpha.prior <- list("dnorm", 0, 1e-6)
  }

  if(is.null(beta.prior)){
    if(response == "ordinal"){
      beta.prior <- list("dnorm", 0, 0.25)
    } else{
      beta.prior <- list("dnorm", 0, 1e-6)
    }
  }

  if(is.null(eta.prior)){
    eta.prior <- list("dnorm", 0, 1e-6)
  }

  if(is.null(rho.prior)){
    rho.prior <- list("dunif", -1, 1)
  }

  if(is.null(hy.prior)){
    if(response == "binomial"){
      hy.prior <- list("dunif", 0, 100)
    } else if (response == "normal"){
      hy.prior <- list("dgamma", 0.001, 0.001)
      # sd <- (max(Y, na.rm = T) - min(Y, na.rm = T)) / 2
      # hy.prior <- list("dunif", 0, 50)
    } else{
      hy.prior <- list("dgamma", 1, 0.001)
    }

  }

  prior.param.default <- list(alpha.prior = alpha.prior, beta.prior = beta.prior, eta.prior = eta.prior, rho.prior = rho.prior, hy.prior = hy.prior)

  if(response == "ordinal"){
    if(is.null(dc.prior)){
      dc.prior <- list("dunif", 0, 20)
    }
    if(is.null(c1.prior)){
      c1.prior <- list("dunif", -20, 20)
    }
    prior.param.default$dc.prior <- dc.prior
    prior.param.default$c1.prior <- c1.prior
  }

  return(prior.param.default)
  })
}

