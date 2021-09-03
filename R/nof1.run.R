#' Run the model using the nof1 object for individual analysis
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

  # inits = NULL
  # n.chains = 3
  # max.run = 100000
  # setsize = 10000
  # n.run = 50000
  # conv.limit = 1.05
  # extra.pars.save = NULL

  if (!inherits(nof1, "nof1.data")) {
    stop('Given object is not nof1.data. Run nof1.data function first')
  }

  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }

  # Originally within object nof1 using with() function
  if (nof1$response == "ordinal") {
    pars.save <- "alpha"

    if(nof1$n.Treat >= 2) {
      for(i in 2:nof1$n.Treat){
        pars.save <- c(pars.save, paste0("beta_", nof1$Treat.name[i]))
      }
    }

  } else {
    pars.save <- NULL
    for(i in nof1$Treat.name){
      pars.save <- c(pars.save, paste0("beta_", i))
    }
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

  if (nof1$response == "ordinal") {

    if(nof1$n.Treat >= 2) {
      for(i in 2:nof1$n.Treat){
        data[[paste0("Treat_", nof1$Treat.name[i])]] <- nof1[[paste0("Treat_", nof1$Treat.name[i])]]
        data[[paste0("Treat_", nof1$Treat.name[i])]][is.na(data[[paste0("Treat_", nof1$Treat.name[i])]])] <- 0
      }
    }

  } else {
    for(i in nof1$Treat.name){
      data[[paste0("Treat_", i)]] <- nof1[[paste0("Treat_", i)]]
      data[[paste0("Treat_", i)]][is.na(data[[paste0("Treat_", i)]])] <- 0
    }
  }

  if (nof1$bs.trend) {
    for (i in 1:nof1$bs_df){
      data[[paste0("bs", i)]] <- nof1[[paste0("bs", i)]]
    }
  }

  if (is.null(inits)) {
    inits <- nof1.inits(nof1, n.chains)
  }

  samples <- jags.fit(nof1, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)
  # head(samples$samples[[1]])

  result <- list(nof1 = nof1, inits = inits, pars.save = pars.save, data.rjags = data)
  result <- c(result, samples)

  class(result) <- "nof1.result"
  return(result)

}

#' @export

nof1.ma.run <- function(nof1, inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                        conv.limit = 1.05, extra.pars.save = NULL) {

  # inits = NULL
  # n.chains = 3
  # max.run = 100000
  # setsize = 10000
  # n.run = 50000
  # conv.limit = 1.05
  # extra.pars.save = NULL

  if (!inherits(nof1, "nof1.data")) {
    stop('Given object is not nof1.data. Run nof1.data function first')
  }

  if (max.run < setsize) {
    stop("setsize should be smaller than max.run")
  }

  # Create pars.save
  if (nof1$response == "ordinal") {
    pars.save <- c("b", "sigmaSq_beta", "rho")
    for(Treat.name.i in 1:length(nof1$Treat.name)){
      pars.save <- c(pars.save, paste0("beta_", nof1$Treat.name[Treat.name.i]))
    }

  } else {
    if (nof1$model.intcpt == "fixed") {
      pars.save <- c("alpha", "d", "sigmaSq_d", "rho")
    } else if (nof1$model.intcpt == "random") {
      pars.save <- c("b", "sigmaSq_beta", "rho", paste0("beta_", nof1$Treat.name[1]))
    } else if (nof1$model.intcpt == "common") {
      pars.save <- paste0("beta_", nof1$Treat.name[1])
    }

    for(Treat.name.i in 2:length(nof1$Treat.name)){
      pars.save <- c(pars.save, paste0("beta_", nof1$Treat.name[Treat.name.i]))
    }
  }

  if(nof1$response == "normal"){
    pars.save <- c(pars.save, "sd_resid")
  }

  # adjust for covariates
  if (!is.null(nof1$names.covariates)) {
    pars.save <- c(pars.save, paste0("beta_", nof1$names.covariates))
  }

  # trend
  if ((nof1$spline.trend) | (nof1$step.trend)) {
    pars.save <- c(pars.save, "eta")
  }

  # extra parameters
  if (!is.null(extra.pars.save)) {
    pars.save <- c(pars.save, extra.pars.save)
  }

  # Create data for fitting jags
  # Y <- nof1$Y
  data <- list(Y = nof1$Y)
  if (nof1$model.intcpt != "fixed") {
    for (Treat.name.i in 1:length(nof1$Treat.name)) {
      data[[paste0("Treat_", nof1$Treat.name[Treat.name.i])]] <- nof1[[paste0("Treat_", nof1$Treat.name[Treat.name.i])]]
    }
  } else {
    for (Treat.name.i in 2:length(nof1$Treat.name)) {
      data[[paste0("Treat_", nof1$Treat.name[Treat.name.i])]] <- nof1[[paste0("Treat_", nof1$Treat.name[Treat.name.i])]]
    }
  }

  data$n.ID    <- nof1$n.ID
  data$nobs.ID <- nof1$nobs.ID
  data$n.Treat <- nof1$n.Treat

  if (nof1$response == "ordinal") {
    data$ord.ncat <- nof1$ord.ncat
  }

  # data on covariates
  if (!is.null(nof1$names.covariates)) {
    for (covariates.i in 1:length(nof1$names.covariates)) {
      data[[nof1$names.covariates[covariates.i]]] <- nof1[[nof1$names.covariates[covariates.i]]]
    }
  }

  # data on trend - splines
  if (nof1$spline.trend) {
    for (i in 1:nof1$spline_df){
      data[[paste0("spline", i)]] <- nof1[[paste0("spline", i)]]
    }
    data$time <- nof1$time
  }

  # data on trend - steps
  if (nof1$step.trend) {
    for (i in nof1$n.steps){
      data[[paste0("step", i)]] <- nof1[[paste0("step", i)]]
    }
  }

  # Initial values
  if (is.null(inits)) {
    inits <- nof1.ma.inits(nof1, n.chains)
  }

  # samples <- jags.fit(nof1, data, pars.save, NULL, n.chains = 3, max.run = 10^5, setsize = 10^4, n.run = 50000, conv.limit = 1.05)
  samples <- jags.fit(nof1, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)

  result <- list(nof1 = nof1, inits = inits, pars.save = pars.save, data.rjags = data)
  result <- c(result, samples)

  class(result) <- "nof1.result"
  return(result)

}

#' Run the model using the nof1 object for (network) meta-analyses
#'
#' This is the core function that runs the model in our program. Before running this function, we need to specify data, prior,
#' JAGS code, etc. using \code{\link{nof1.nma.data}}.
#'
#' @param nof1 nof1 object created from \code{\link{nof1.nma.data}} function
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
#'
#' @export
# Network meta analysis
nof1.nma.run <- function(nof1, inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                         conv.limit = 1.05, extra.pars.save = NULL) {

  # inits = NULL
  # n.chains = 3
  # max.run = 100000
  # setsize = 10000
  # n.run = 50000
  # conv.limit = 1.05
  # extra.pars.save = NULL

  if (!inherits(nof1, "nof1.data")) {
    stop('Given object is not nof1.data. Run nof1.data function first')
  }

  if (max.run < setsize) {
    stop("setsize should be smaller than max.run")
  }

  # Create pars.save
  # if (nof1$response == "ordinal") {
  #   pars.save <- c("b", "sigmaSq_beta", "rho")
  #   for(Treat.name.i in 1:length(nof1$Treat.name)){
  #     pars.save <- c(pars.save, paste0("beta_", nof1$Treat.name[Treat.name.i]))
  #   }
  #
  # } else {
  if (nof1$model.intcpt == "fixed") {
    pars.save <- c("alpha", "d", "prec_beta", paste0("beta_", 2:nrow(nof1$summ.nID.perTreat)))
  } else if (nof1$model.intcpt == "random") {
    pars.save <- c("b", "prec_beta", "rho", paste0("beta_", 1:nrow(nof1$summ.nID.perTreat)))
  }
  # else if (nof1$model.intcpt == "common") {
  #     pars.save <- paste0("beta_", nof1$Treat.name[1])
  #   }
  # }

  if(nof1$response == "normal"){
    pars.save <- c(pars.save, "prec_resid")
  }

  # adjust for level 2 covariates
  if (!is.null(nof1$cov.matrix)) {
    pars.save <- c(pars.save, "eta_cov")
  }

  # trend
  if ((nof1$spline.trend) | (nof1$step.trend)) {
    pars.save <- c(pars.save, "eta")
  }

  # serial correlation
  if (nof1$corr.y) {
    pars.save <- c(pars.save, "rho_resid")
  }

  # stratification covariates used during randomization
  if (!is.null(nof1$fixed.strata.cov.matrix)) {
    pars.save <- c(pars.save, "eta_fixed_strata_cov")
  }

  if (!is.null(nof1$random.strata.cov.matrix)) {
    pars.save <- c(pars.save, "eta_random_strata_cov", "eta_random_strata_cov_lvls", "prec_eta_random_strata_cov")
  }

  # extra parameters
  if (!is.null(extra.pars.save)) {
    pars.save <- c(pars.save, extra.pars.save)
  }

  # Create data for fitting jags
  # Y <- nof1$Y
  data <- list(Y.matrix = nof1$Y.matrix,
               nobs.ID  = nof1$nobs.ID,
               uniq.Treat.matrix = nof1$uniq.Treat.matrix)
  if (nof1$model.intcpt == "random") {
    data$Treat.1       <- nof1$Treat.1
    # data$Treat.order.1 <- nof1$Treat.order.1
  }
  for (Treat.i in 2:nrow(nof1$summ.nID.perTreat)) {
    data[[paste0("Treat.",  Treat.i)]]       <- nof1[[paste0("Treat.", Treat.i)]]
  }

  # for (Treat.i in 2:length(nof1$Treat.name)) {
  #   data[[paste0("Treat.order.",  Treat.i)]] <- nof1[[paste0("Treat.order.", Treat.i)]]
  # }
  # if (nof1$response == "ordinal") {
  #   data$ord.ncat <- nof1$ord.ncat
  # }

  # data on level 2 covariates
  if (!is.null(nof1$cov.matrix)) {
    data$cov.matrix <- nof1$cov.matrix
  }

  # data on trend - splines
  if (nof1$spline.trend) {
    data$spline.matrix <- nof1$spline.matrix
    data$time.matrix   <- nof1$time.matrix
  }

  # data on trend - steps
  if (nof1$step.trend) {
    data$step.matrix   <- nof1$step.matrix
    data$period.matrix <- nof1$period.matrix
  }

  # data on time - correlation
  if (nof1$corr.y) {
    data$time.matrix   <- nof1$time.matrix
  }

  # stratification covariates used during randomization
  if (!is.null(nof1$fixed.strata.cov.matrix)) {
    data$fixed.strata.cov.matrix <- nof1$fixed.strata.cov.matrix
  }

  if (!is.null(nof1$random.strata.cov.matrix)) {
    data$random.strata.cov.matrix <- nof1$random.strata.cov.matrix
  }

  # Initial values
  if (is.null(inits)) {
    inits <- nof1.nma.inits(nof1, n.chains)
  }

  # samples <- jags.fit(nof1, data, pars.save, NULL, n.chains = 3, max.run = 10^5, setsize = 10^4, n.run = 50000, conv.limit = 1.05)
  samples <- jags.fit(nof1, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)

  result <- list(nof1 = nof1, inits = inits, pars.save = pars.save, data.rjags = data)
  result <- c(result, samples)

  class(result) <- "nof1.result"
  return(result)

}



jags.fit <- function(nof1, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit){

  # set.seed(seed)
  mod = rjags::jags.model(textConnection(nof1$code), data = data, inits = inits, n.chains = n.chains, n.adapt = setsize)

  # adapt model
  adapted <- FALSE
  count <- 0
  while(!adapted){
    adapted <- rjags::adapt(mod, setsize, end.adaptation = FALSE)
    count <- count + 1
    if(count == 100){
      stop("algorithm has not adapted")
    }
  }

  # draw samples
  samples <- rjags::coda.samples(model = mod, variable.names = pars.save, n.iter  = setsize)

  # check convergence
  max.gelman <- find.max.gelman(samples)
  print(max.gelman)
  check <- max.gelman > conv.limit

  # draw samples until convergence
  if(check) {
    count <- 1
    while (check & (count < max.run/setsize)) {
      samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = setsize)
      samples <- add.mcmc(samples, samples2)

      count <- count + 1

      max.gelman <- find.max.gelman(samples)
      check <- max.gelman > conv.limit
      print(max.gelman)
    }
  }

  # create the mcmc object
  start <- mcpar(samples[[1]])[1]
  end <- mcpar(samples[[1]])[2]
  mid <- (end + start-1)/2
  burnin <- ceiling(end - mid)
  samples <- window(samples, start = mid+1, end = end, frequency = 1) #keep the last half of the converged sequence
  samples <- new.mcmc(samples)

  # check if the samples reach the number of n.run
  n.thin <- 1
  if (check) {
    stop("code didn't converge according to gelman-rubin diagnostics")
  } else if (n.run < burnin) {
    n.thin <- ceiling(burnin/n.run)
    extra.run <- n.run * n.thin - burnin
    if(extra.run != 0){
      samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = extra.run)
      samples <- add.mcmc(samples, samples2)
    }
    samples <- window(samples, 1, dim(samples[[1]])[1], n.thin)
  } else if (n.run > burnin) {
    extra.run <- n.run - burnin
    samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = extra.run)
    samples <- add.mcmc(samples, samples2)
    # max.gelman <- find.max.gelman(samples)
    # print(max.gelman)
  }

  max.gelman <- find.max.gelman(samples)
  print(max.gelman)
  check <- max.gelman > conv.limit

  # redraw samples if jumped out of convergence range
  if (check) {

    while (check & (count < max.run/setsize)) {
      # last.count <- count
      samples2   <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = setsize)
      samples    <- add.mcmc(samples, samples2)

      count <- count + 1

      samples <- window(samples,
                        start = setsize + 1,
                        end   = n.run + setsize,
                        frequency = 1)
      # for (i in 1:length(samples)) {
      #   tmp.samples[[i]] <- samples[[i]][round(seq(1, dim(samples[[i]])[1], length.out = n.run)), ]
      # }
      # tmp.samples <- new.mcmc(tmp.samples)

      max.gelman <- find.max.gelman(samples)
      check <- max.gelman > conv.limit
      print(max.gelman)
    }
    # samples <- tmp.samples
  }

  # find DIC
  dic <- dic.samples(mod, n.iter = 10^4)

  if(check){
    stop("code didn't converge according to gelman-rubin diagnostics")
  }

  out <-list(burnin = burnin, n.thin = n.thin, samples = samples, max.gelman = max.gelman, dic = dic)
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

