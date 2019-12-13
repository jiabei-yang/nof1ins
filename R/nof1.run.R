#' Run the model using the nof1 object
#' 
#' This is the core function that runs the model in our program. Before running this function, we need to specify data, prior, 
#' JAGS code, etc. using \code{\link{nof1.data}}.
#'
#' @param nof1 nof1 object created from \code{\link{nof1.data}} function
#' @param inits Initial values for the parameters being sampled. If left unspecified, program will generate 
#' reasonable initial values.
#' @param n.chains Number of chains to run
#' @param max.run Maximum number of iterations that user is willing to run. If the algorithm is not 
#' converging, it will run up to \code{max.run} iterations before printing a message that it did not converge
#' @param setsize Number of iterations that are run between convergence checks. If the algorithm converges 
#' fast, user wouldn't need a big setsize. The number that is printed between each convergence checks is the 
#' gelman-rubin diagnostics and we would want that to be below the conv.limit the user specifies.
#' @param n.run Final number of iterations that the user wants to store. If after the algorithm converges, 
#' user wants less number of iterations, we thin the sequence. If the user wants more iterations, we run 
#' extra iterations to reach the specified number of runs
#' @param conv.limit Convergence limit for Gelman and Rubin's convergence diagnostic.
#' @param extra.pars.save Parameters that user wants to save besides the default parameters saved. See code 
#' using \code{cat(nof1$code)} to see which parameters can be saved.
#' @return
#' \item{nof1}{nof1 object}
#' \item{inits}{Initial values that are either specified by the user or generated as a default}
#' \item{pars.save}{Parameters that are saved. Add more parameters in extra.pars.save if other variables are desired}
#' \item{data_rjags}{Data that is put into rjags function \code{jags.model}}
#' \item{burnin}{Half of the converged sequence is thrown out as a burnin}
#' \item{n.thin}{If the number of iterations user wants (n.run) is less than the number of converged sequence after burnin, we thin the sequence and store the thinning interval}
#' \item{samples}{MCMC samples stored using jags. The returned samples have the form of mcmc.list and can be directly applied to coda functions}
#' \item{max.gelman}{Maximum Gelman and Rubin's convergence diagnostic calculated for the final sample}
#' @examples
#' laughter
#' Y <- laughter$Y
#' Treat <- laughter$Treat
#' nof1 <- nof1.data(Y, Treat, ncat = 11, baseline = "Usual Routine", response = "ordinal")
#' result <- nof1.run(nof1)
#' summary(result$samples)
#' @export

nof1.run <- function(nof1, inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                     conv.limit = 1.05, extra.pars.save = NULL){

  if (!inherits(nof1, "nof1.data")) {
    stop('Given object is not nof1.data. Run nof1.data function first')
  }

  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }

  # Originally within object nof1 using with() function
  if (nof1$response == "ordinal"){
    pars.save <- "c"
  } else {
    pars.save <- NULL
  }
  
  for(i in nof1$Treat.name){
    pars.save <- c(pars.save, paste0("beta_", i))
  }
  
  if(nof1$response == "normal"){
    pars.save <- c(pars.save, "sd")
  }
  
  # bs_df <- nof1$bs_df
  if (nof1$bs.trend) {
    for(i in 1:nof1$bs_df){
      pars.save <- c(pars.save, paste0("eta_", i))
    }
  }
  
  if (nof1$corr.y) {
    pars.save <- c(pars.save, "rho")
  }
  
  # Y <- nof1$Y
  data <- list(Y = nof1$Y)
  for(i in nof1$Treat.name){
    data[[paste0("Treat_", i)]] <- nof1[[paste0("Treat_", i)]]
  }
  
  if (nof1$bs.trend) {
    for (i in 1:nof1$bs_df){
      data[[paste0("bs", i)]] <- nof1[[paste0("bs", i)]]
    }
  }
  
  if(is.null(inits)){
    inits <- nof1.inits(nof1, n.chains)
  }
  
  samples <- jags.fit(nof1, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)
  
  result <- list(nof1 = nof1, inits = inits, pars.save = pars.save, data.rjags = data)
  result <- c(result, samples)
  
  class(result) <- "nof1.result"
  return(result)
  
  # with(nof1, {
  # 
  # # response <- nof1$response
  # 
  #   if (response == "ordinal"){
  #     pars.save <- "c"
  #   } else {
  #     pars.save <- NULL
  #   }
  # 
  # # Treat.name <- nof1$Treat.name
  #   for(i in Treat.name){
  #     pars.save <- c(pars.save, paste0("beta_", i))
  #   }
  # 
  #   if(response == "normal"){
  #     pars.save <- c(pars.save, "sd")
  #   }
  #   
  #   # bs_df <- nof1$bs_df
  #   if (bs.trend) {
  #     for(i in 1:bs_df){
  #       pars.save <- c(pars.save, paste0("eta_", i))
  #     }
  #   }
  # 
  #   if (corr.y) {
  #     pars.save <- c(pars.save, "rho")
  #   }
  #   
  #   # Y <- nof1$Y
  #   data <- list(Y = Y)
  #   for(i in Treat.name){
  #     data[[paste0("Treat_", i)]] <- nof1[[paste0("Treat_", i)]]
  #   }
  # 
  #   if (bs.trend) {
  #     for (i in 1:bs_df){
  #       data[[paste0("bs", i)]] <- nof1[[paste0("bs", i)]]
  #     }
  #   }
  #   
  #   if(is.null(inits)){
  #     inits <- nof1.inits(nof1, n.chains)
  #   }
  #   samples <- jags.fit(nof1, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)
  #   
  #   result <- list(nof1 = nof1, inits = inits, pars.save = pars.save, data.rjags = data)
  #   result <- c(result, samples)
  #   
  #   class(result) <- "nof1.result"
  #   return(result)
  # })
}

jags.fit <- function(nof1, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit){

  mod = rjags::jags.model(textConnection(nof1$code), data = data, inits = inits, n.chains = n.chains, n.adapt = setsize)
  
  adapted <- FALSE
  count <- 0
  while(!adapted){
    adapted <- rjags::adapt(mod, setsize, end.adaptation = FALSE)
    count <- count + 1
    if(count == 100){
      stop("algorithm has not adapted")
    }
  }
  
  samples <- rjags::coda.samples(model = mod, variable.names = pars.save, n.iter = setsize)

  max.gelman <- find.max.gelman(samples)
  print(max.gelman)
  check <- max.gelman > conv.limit

  if(check) {
    count <- 1
    while (check & count < max.run/setsize) {
      samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = setsize)
      samples <- add.mcmc(samples, samples2)

      count <- count + 1

      max.gelman <- find.max.gelman(samples)
      check <- max.gelman > conv.limit
      print(max.gelman)
    }
  }

  start <- mcpar(samples[[1]])[1]
  end <- mcpar(samples[[1]])[2]
  mid <- (end + start-1)/2
  burnin <- ceiling(end - mid)
  samples <- window(samples, start = mid+1, end = end, frequency = 1) #keep the last half of the converged sequence
  samples <- new.mcmc(samples)

  n.thin <- 1
  if(check){
    stop("code didn't converge according to gelman-rubin diagnostics")
  } else if(n.run < burnin){
    n.thin <- ceiling(burnin/n.run)
    extra.run <- n.run * n.thin - burnin
    if(extra.run != 0){
      samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = extra.run)
      samples <- add.mcmc(samples, samples2)
    }
    samples <- window(samples, 1, dim(samples[[1]])[1], n.thin)
  } else if(n.run > burnin){
    extra.run <- n.run - burnin
    samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = extra.run)
    samples <- add.mcmc(samples, samples2)
  }
  max.gelman <- find.max.gelman(samples)
  print(max.gelman)

  out <-list(burnin = burnin, n.thin = n.thin, samples = samples, max.gelman = max.gelman)
  return(out)
}


new.mcmc <- function(x){
  n.chains <- length(x)
  n.var <- coda::nvar(x)
  newobjects <- vector("list", length = n.chains)

  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.var, dimnames = list(NULL, dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- x[[i]]
    newobjects[[i]] <- mcmc(newobjects[[i]])
  }
  mcmc.list(newobjects)
}

add.mcmc <- function(x, y){

  n.chains <- length(x)
  n.var <- coda::nvar(x)
  newobjects <- vector("list", length = n.chains)

  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.var, dimnames = list(NULL, dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- rbind(x[[i]], y[[i]])
    newobjects[[i]] <- mcmc(newobjects[[i]])
  }
  mcmc.list(newobjects)
}

find.max.gelman <- function(samples, index = NULL){

  if(!is.null(index)){
    samples <- lapply(samples, function(x){ x[,index]})
  }
  samples <- lapply(samples, function(x) { x[,colSums(abs(x)) != 0] })

  max(gelman.diag(samples, multivariate = FALSE)$psrf[,1]) #look at point estimate instead of 95% C.I.
}

