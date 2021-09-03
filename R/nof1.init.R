nof1.inits <- function(nof1, n.chains){

  response <- nof1$response
  inits <- tryCatch({
    value <- if(response == "normal"){
      nof1.inits.normal(nof1, n.chains)
    } else if(response == "ordinal"){
      nof1.inits.ordinal(nof1, n.chains)
    } else if(response == "binomial" || response == "poisson"){
      nof1.inits.binom.poisson(nof1, n.chains)
    }
    if(any(is.nan(unlist(value)))) value <- NULL

    #if initial value generated is too big (i.e. model didn't work because of sparse data), just set inital value to be NULL
    if(!is.null(value)){
      if(any(abs(unlist(value)) > 100)) value <- NULL
    }
    value
  }, error = function(err){
    print(paste("Initial value not working: ",err))
    return(NULL)
  }, warning = function(warning){
    print(paste("Warning: ", warning))
    return(NULL)
  })
  return(inits)
}

nof1.inits.normal <- function(nof1, n.chains){

  with(nof1, {

    Treat.matrix <- NULL
    for(i in Treat.name){
      Treat.matrix <- cbind(Treat.matrix, nof1[[paste0("Treat_", i)]])
    }

    if (exists("bs_df")){

      for (i in 1:bs_df){
        Treat.matrix <- cbind(Treat.matrix, nof1[[paste0("bs", i)]])
      }
    }

    model <- lm(Y ~ Treat.matrix - 1)
    co <- coef(summary(model))

    ############# Generate initial values
    initial.values = list()
    for(i in 1:n.chains){
      initial.values[[i]] = list()
    }
    for(i in 1:n.chains){

      # initial.values[[i]][["alpha"]] <- co[1,1] + rnorm(1) *co[1,2]

      for(j in 1:length(Treat.name)){
        initial.values[[i]][[paste0("beta_", Treat.name[j])]] <- co[j,1] + rnorm(1) * co[j,2]

        # Take care of bounded beta values if uniform prior prespecified
        if (beta.prior[[1]] == "dunif") {
          initial.values[[i]][[paste0("beta_", Treat.name[j])]] <- max(initial.values[[i]][[paste0("beta_", Treat.name[j])]], beta.prior[[2]])
          initial.values[[i]][[paste0("beta_", Treat.name[j])]] <- min(initial.values[[i]][[paste0("beta_", Treat.name[j])]], beta.prior[[3]])
        }

      }

      if(exists("bs_df")){
        for(j in 1:bs_df){
          initial.values[[i]][[paste0("eta_", j)]] <- co[length(Treat.name)+j, 1] + rnorm(1, 0, co[length(Treat.name)+j, 2])
        }
      }
    }

    # if(!is.nan(summary(model)$fstat[1])){
    #   for(i in 1:n.chains){
    #
    #     fit <- summary(model)
    #     df <- fit$df[2]
    #     random.ISigma <- rchisq(1, df)
    #     resid.var <- fit$sigma^2
    #     sigma2 <- resid.var * df/ random.ISigma
    #
    #     if(hy.prior[[1]] == "dunif"){
    #       if(sqrt(sigma2) > hy.prior[[3]]){
    #         stop("data has more variability than your prior does")
    #       }
    #     }
    #
    #     if(hy.prior[[1]] == "dgamma"){
    #       initial.values[[i]][["prec"]] <- 1/sigma2
    #     } else if((hy.prior[[1]] == "dnorm") | (hy.prior[[1]] == "dunif")){
    #       initial.values[[i]][["sd"]] <- sqrt(sigma2)
    #     }
    #   }
    # }
    return(initial.values)
  })

}

nof1.inits.binom.poisson <- function(nof1, n.chains){

  with(nof1, {

    Treat.matrix <- NULL
    for(i in Treat.name){
      Treat.matrix <- cbind(Treat.matrix, nof1[[paste0("Treat_", i)]])
    }

    if(response == "binomial"){

      model <- glm(Y ~ Treat.matrix - 1, family = binomial(link = "logit"))
      co <- coef(summary(model))
      # else{
      #   model <- glm(Y ~ Treat.matrix + BS, family = binomial(link = "logit"))
      #   co <- coef(summary(model))
      # }
    } else if(response == "poisson"){

      model <- glm(Y ~ Treat.matrix - 1, family = "poisson")
      co <- coef(summary(model))
      # else{
      #   model <- glm(Y ~ Treat.matrix + BS, family = "poisson")
      #   co <- coef(summary(model))
      # }
    }

    initial.values = list()
    for(i in 1:n.chains){
      initial.values[[i]] = list()
    }
    for(i in 1:n.chains){

      # initial.values[[i]][["alpha"]] <- co[1,1] + rnorm(1) *co[1,2]

      for(j in 1:length(Treat.name)){
        # initial.values[[i]][[paste0("beta_", Treat.name[j])]] <- co[1+j,1] + rnorm(1) * co[1+j,2]
        initial.values[[i]][[paste0("beta_", Treat.name[j])]] <- co[j, 1] + rnorm(1) * co[j, 2]
      }

      # if(!is.null(knots)){
      #   for(j in 1:ncol(BS)){
      #     initial.values[[i]][[paste0("gamma", j)]] <- co[1+length(Treat.name)+j, 1] + rnorm(1) * co[1+length(Treat.name)+j, 2]
      #   }
      # }
    }
    return(initial.values)
  })
}



nof1.inits.ordinal <- function(nof1, n.chains){

  Y_matrix <- NULL
  for (i in 1:nof1$ncat) {
    Y_matrix <- cbind(Y_matrix,
                      as.numeric(nof1$Y == i))
  }

  Treat.matrix <- NULL
  for(i in nof1$Treat.name[2:length(nof1$Treat.name)]){
    Treat.matrix <- cbind(Treat.matrix, nof1[[paste0("Treat_", i)]])
  }

  if (nof1$ord.model == "cumulative") {
    fit.resAc <- vglm(Y_matrix ~ Treat.matrix, cumulative(parallel = T))
    co <- coef(summary(fit.resAc))
  } else {
    fit.resAc <- vglm(Y_matrix ~ Treat.matrix, acat(parallel = T))
    co <- coef(summary(fit.resAc))
  }

  ############# Generate initial values
  initial.values = list()
  for(i in 1:n.chains){
    initial.values[[i]] = list()
  }

  for(i in 1:n.chains){

    initial.values[[i]][["alpha"]] <- NULL
    # only initilize alpha if all contrasts exist
    # otherwise we do not know which levels are skipped
    n.contrasts <- sum(grepl("Intercept", rownames(co)))
    if (n.contrasts == (nof1$ncat - 1)) {
      for (j in 2:nof1$ncat) {
        initial.values[[i]][["alpha"]] <- c(initial.values[[i]][["alpha"]],
                                            co[j-1, 1] + rnorm(1) * co[j-1, 2])
      }
    }

    for (j in 2:length(nof1$Treat.name)) {
      initial.values[[i]][[paste0("beta_", nof1$Treat.name[j])]] <- co[n.contrasts + (j-1), 1] + rnorm(1) * co[n.contrasts + (j-1), 2]

    }

  }

  # with(nof1, {
  #
  #   p <- rep(NA, ncat)
  #   c <- rep(NA, ncat-1)
  #
  #   for (i in seq(ncat)) {
  #     p[i] = sum(Y[!is.na(Y)]==i)/nobs
  #     if (!is.na(p[i]) & p[i] == 0) { p[i] = 0.05 }
  #   }
  #   if(sum(p[!is.na(p)])> 1) {
  #     p[max(which(p == max(p)))] = p[max(which(p == max(p)))] + 1 - sum(p)
  #   }
  #
  #   initial.values = list()
  #   for(i in 1:n.chains){
  #     initial.values[[i]] = list()
  #   }
  #
  #   for(i in 1:n.chains){
  #     p <- combinat::rmultz2(nobs,p)/nobs
  #     if (any(p == 0)) {
  #       p[which(p == 0)] <- 0.05
  #       p[max(which(p == max(p)))] <- p[max(which(p == max(p)))] + 1 - sum(p)
  #     }
  #
  #     q <- 1 - cumsum(p)
  #     for(j in seq(ncat -1)){
  #       c[j] <- -log(q[j]/(1-q[j]))
  #     }
  #     dc <- c(c[1], c[-1] - c[-(ncat-1)])
  #     initial.values[[i]][["dc"]] <- dc
  #   }
  #
  #
  #   Treat.matrix <- NULL
  #   for(i in Treat.name){
  #     Treat.matrix <- cbind(Treat.matrix, nof1[[paste0("Treat_", i)]])
  #   }
  #
  #   # if(is.na(knots)){
  #   #   model <- polr(as.ordered(Y) ~ Treat.matrix, Hess = TRUE)
  #   #   co = coef(summary(model))
  #   # }
  #   model <- MASS::polr(as.ordered(Y) ~ Treat.matrix, Hess = TRUE)
  #   co = coef(summary(model))
  #
  #   if(!is.null(model)){
  #     co_Treat <- co[grep('Treat.matrix', rownames(coef(summary(model)))),,drop = FALSE]
  #     #co_BS <- co[grep('BS', rownames(coef(summary(model)))), ]
  #
  #     for(i in 1:n.chains){
  #       for(j in 1:length(Treat.name)){
  #         initial.values[[i]][[paste0("beta_", Treat.name[j])]] <- co_Treat[j,1] + rnorm(1) * co_Treat[j,2]
  #       }
  #
  #       # if(!is.na(knots)){
  #       #   for(j in 1:ncol(BS)){
  #       #     initial.values[[i]][[paste0("gamma", j)]] <- co_BS[j, 1] + rnorm(1) * co_BS[j, 2]
  #       #   }
  #       # }
  #     }
  #   }

  return(initial.values)
  # })
}


# Meta analysis
nof1.ma.inits <- function(nof1, n.chains) {

  if (nof1$response == "normal") {
    initial.values <- nof1.ma.normal.inits(nof1, n.chains)
  } else if (nof1$response == "poisson") {
    initial.values <- nof1.ma.poisson.inits(nof1, n.chains)
  } else if (nof1$response == "binomial") {
    initial.values <- nof1.ma.binomial.inits(nof1, n.chains)
  } else if (nof1$response == "ordinal") {
    # if ordinal outcome, no initial value because no equivalent frequentist model can be fit
    initial.values <- nof1.ma.ordinal.inits(nof1, n.chains)
  }

  # if (!is.null(initial.values)) {
  # if initial value generated is too big (i.e. model didn't work because of sparse data), just set inital value to be NULL
  if (any(abs(unlist(initial.values)) > 100, na.rm = T)) {
    initial.values <- NULL
  }
  # }

  return(initial.values)
}



nof1.ma.normal.inits <- function(nof1, n.chains) {

  if (nof1$model.intcpt == "fixed") {
    tmp.formula <- "nof1$Y.long ~ -1 + nof1$ID + nof1$Treat"
  } else if (nof1$model.intcpt == "random") {
    tmp.formula <- "nof1$Y.long ~ -1 + nof1$Treat"
  }

  if (!is.null(nof1$names.covariates)) {
    for (covariates.i in 1:length(nof1$names.long.covariates)) {
      tmp.formula <- paste0(tmp.formula, " + nof1$", nof1$names.long.covariates[covariates.i])
    }
  }

  if (nof1$spline.trend) {
    for (spline_df.i in 1:nof1$spline_df) {
      tmp.formula <- paste0(tmp.formula, " + nof1$spline.long", spline_df.i)
    }
  }

  if (nof1$step.trend) {
    for (step.i in nof1$n.steps) {
      tmp.formula <- paste0(tmp.formula, " + (nof1$y.step == ", step.i, ")")
    }
  }

  model <- glm(as.formula(tmp.formula), family = gaussian(link = nof1$model.linkfunc))

  summ.model    <- summary(model)
  coefMat.model <- coef(summ.model)
  names.coef    <- rownames(coefMat.model)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]]   <- list()

    if (nof1$model.intcpt == "common") {

      for (n.Treat.i in 1:nof1$n.Treat) {

        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[n.Treat.i]), names.coef)
        initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[n.Treat.i])]] <-
          rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2])
      }

    } else {
      if (nof1$model.intcpt == "fixed") {
        # d
        initial.values[[n.chains.i]]$d <- NULL
        for (Treat.name.i in 2:nof1$n.Treat) {

          ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.name.i]), names.coef)
          initial.values[[n.chains.i]]$d <- c(initial.values[[n.chains.i]]$d,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

          # initialize beta
          # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
        }

        # alpha
        initial.values[[n.chains.i]]$alpha <- NULL
        for (id.i in 1:nof1$n.ID) {

          ind.coefMat <- grepl(nof1$uniq.ID[id.i], names.coef)
          initial.values[[n.chains.i]]$alpha <- c(initial.values[[n.chains.i]]$alpha,
                                                  rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        }

      } else if (nof1$model.intcpt == "random") { # if (nof1$model.intcpt == "fixed") {
        # b
        initial.values[[n.chains.i]]$b <- NULL
        for (n.Treat.i in 1:nof1$n.Treat) {

          ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[n.Treat.i]), names.coef)
          initial.values[[n.chains.i]]$b <- c(initial.values[[n.chains.i]]$b,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        }
      }

      # beta's
      # beta and variance of beta are not settable
      for (id.i in 1:nof1$n.ID) {

        # beta
        ind.id        <- which((nof1$ID == nof1$uniq.ID[id.i]) & (!is.na(nof1$Y.long)))
        if (length(unique(nof1$Treat[ind.id])) == 1) {
          model.id.orig <- lm(nof1$Y.long[ind.id] ~ 1)
        } else {
          if (nof1$model.intcpt == "fixed") {
            model.id.orig <- lm(nof1$Y.long[ind.id] ~ nof1$Treat[ind.id])
          } else {
            model.id.orig <- lm(nof1$Y.long[ind.id] ~ -1 + nof1$Treat[ind.id])
          }
        }
        coefMat.model.id.orig <- coef(summary(model.id.orig))

        if (nof1$model.intcpt == "random") {
          tmp.ind  <- grepl(paste0("]", nof1$Treat.name[1]), rownames(coefMat.model.id.orig))
          # if there is no this treatment on this participant, generate grand mean
          # minus 1 because of there is a reference treatment
          if (sum(tmp.ind) == 0) {
            tmp.coef <- initial.values[[n.chains.i]]$b[1]
          } else if (is.nan(coefMat.model.id.orig[tmp.ind, 2])) { # no standard error
            tmp.coef <- initial.values[[n.chains.i]]$b[1]

          } else {
            tmp.coef <- rnorm(1, mean = coefMat.model.id.orig[tmp.ind, 1], sd = coefMat.model.id.orig[tmp.ind, 2])
          }

          initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]] <-
            c(initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]], tmp.coef)
        }

        for (Treat.name.i in 2:nof1$n.Treat) {

          tmp.ind  <- grepl(paste0("]", nof1$Treat.name[Treat.name.i]), rownames(coefMat.model.id.orig))
          # if there is no this treatment on this participant, generate grand mean
          # minus 1 because of there is a reference treatment
          if (sum(tmp.ind) == 0) {
            if (nof1$model.intcpt == "fixed") {
              tmp.coef <- initial.values[[n.chains.i]]$d[Treat.name.i - 1]
            } else {
              tmp.coef <- initial.values[[n.chains.i]]$b[Treat.name.i]
            }
          } else if (is.nan(coefMat.model.id.orig[tmp.ind, 2])) { # no standard error
            if (nof1$model.intcpt == "fixed") {
              tmp.coef <- initial.values[[n.chains.i]]$d[Treat.name.i - 1]
            } else {
              tmp.coef <- initial.values[[n.chains.i]]$b[Treat.name.i]
            }

          } else {
            tmp.coef <- rnorm(1, mean = coefMat.model.id.orig[tmp.ind, 1], sd = coefMat.model.id.orig[tmp.ind, 2])
          }
          # if (length(tmp.coef) == 0) {
          #   tmp.coef <- coef.model.id.orig[grepl(paste0("]", nof1$Treat.name[Treat.name.i]), names(coef.model.id.orig))]
          # }

          initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <-
            c(initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]], tmp.coef)
        }

      } # for id.i

      # sigmaSq_d
      if (nof1$model.intcpt == "fixed") {
        tmp.betas <- NULL
        initial.values[[n.chains.i]]$beta <- NULL
      } else {
        tmp.betas <- initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]]
        initial.values[[n.chains.i]]$beta <- initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]]
        initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]] <- NULL
      }

      for (Treat.name.i in 2:nof1$n.Treat) {
        tmp.betas <- c(tmp.betas, initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])

        # assign values to beta matrix and remove those treatment specific beta's
        initial.values[[n.chains.i]]$beta <- cbind(initial.values[[n.chains.i]][["beta", exact = TRUE]], initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])
        initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      if (nof1$model.intcpt == "fixed") {
        initial.values[[n.chains.i]]$prec_delta <- 1 / var(tmp.betas)
      } else {
        initial.values[[n.chains.i]]$prec_beta <- 1 / var(tmp.betas)
      }
    } #  if (nof1$model.intcpt == "common") { else


    # adjust for trend - splines
    if (nof1$spline.trend) {
      # eta
      initial.values[[n.chains.i]]$eta <- NULL
      for (eta.i in 1:nof1$spline_df) {

        ind.coefMat <- grepl(paste0("spline.long", eta.i), names.coef)
        initial.values[[n.chains.i]]$eta <- c(initial.values[[n.chains.i]]$eta,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta > 100] <- 100
      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta < -100] <- -100
    }

    # adjust for trend - steps
    if (nof1$step.trend) {
      # eta
      initial.values[[n.chains.i]]$eta <- NULL
      for (eta.i in nof1$n.steps) {

        ind.coefMat <- grepl(paste0("y.step == ", eta.i), names.coef)
        initial.values[[n.chains.i]]$eta <- c(initial.values[[n.chains.i]]$eta,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta > 100] <- 100
      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta < -100] <- -100
    }

    # se_resid
    # glm residuals not easily converted or extracted
    # skip now
    # initial.values[[n.chains.i]]$prec_resid <- 1 / (summ.model$sigma)^2 * summ.model$df[2]/rchisq(1, df = summ.model$df[2])

  } # n.chains.i

  return(initial.values)

}



nof1.ma.poisson.inits <- function(nof1, n.chains) {

  if (nof1$model.intcpt == "fixed") {
    tmp.formula <- "nof1$Y.long ~ -1 + nof1$ID + nof1$Treat"
  } else if (nof1$model.intcpt == "random") {
    tmp.formula <- "nof1$Y.long ~ -1 + nof1$Treat"
  }

  if (!is.null(nof1$names.covariates)) {
    for (covariates.i in 1:length(nof1$names.long.covariates)) {
      tmp.formula <- paste0(tmp.formula, " + nof1$", nof1$names.long.covariates[covariates.i])
    }
  }

  if (nof1$spline.trend) {
    for (bs_df.i in 1:nof1$spline_df) {
      tmp.formula <- paste0(tmp.formula, " + nof1$spline.long", bs_df.i)
    }
  }

  model <- glm(as.formula(tmp.formula), family = poisson(link = nof1$model.linkfunc))

  summ.model    <- summary(model)
  coefMat.model <- coef(summ.model)
  names.coef    <- rownames(coefMat.model)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]]   <- list()

    if (nof1$model.intcpt == "fixed") {
      # d
      initial.values[[n.chains.i]]$d <- NULL
      for (Treat.name.i in 2:nof1$n.Treat) {

        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.name.i]), names.coef)
        initial.values[[n.chains.i]]$d <- c(initial.values[[n.chains.i]]$d,
                                            rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      # alpha
      initial.values[[n.chains.i]]$alpha <- NULL
      for (id.i in 1:nof1$n.ID) {

        ind.coefMat <- grepl(nof1$uniq.ID[id.i], names.coef)
        initial.values[[n.chains.i]]$alpha <- c(initial.values[[n.chains.i]]$alpha,
                                                rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
      }

    } else if (nof1$model.intcpt == "random") { # if (nof1$model.intcpt == "fixed") {
      # b
      initial.values[[n.chains.i]]$b <- NULL
      for (n.Treat.i in 1:nof1$n.Treat) {
        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[n.Treat.i]), names.coef)
        initial.values[[n.chains.i]]$b <- c(initial.values[[n.chains.i]]$b,
                                            rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
      }
    }

    # beta's
    # beta and variance of beta are not settable
    for (id.i in 1:nof1$n.ID) {

      # beta
      ind.id        <- which((nof1$ID == nof1$uniq.ID[id.i]) & (!is.na(nof1$Y.long)))
      if (length(unique(nof1$Treat[ind.id])) == 1) {
        model.id.orig <- glm(nof1$Y.long[ind.id] ~ 1, family = poisson(link = "log"))
      } else {
        if (nof1$model.intcpt == "fixed") {
          model.id.orig <- glm(nof1$Y.long[ind.id] ~ nof1$Treat[ind.id], family = poisson(link = "log"))
        } else {
          model.id.orig <- glm(nof1$Y.long[ind.id] ~ -1 + nof1$Treat[ind.id], family = poisson(link = "log"))
        }
      }
      coefMat.model.id.orig <- coef(summary(model.id.orig))

      if (nof1$model.intcpt == "random") {
        tmp.ind  <- grepl(paste0("]", nof1$Treat.name[1]), rownames(coefMat.model.id.orig))
        # if there is no this treatment on this participant, generate grand mean
        # minus 1 because of there is a reference treatment
        if (sum(tmp.ind) == 0) {
          tmp.coef <- initial.values[[n.chains.i]]$b[1]
        } else if (is.nan(coefMat.model.id.orig[tmp.ind, 2])) { # no standard error
          tmp.coef <- initial.values[[n.chains.i]]$b[1]
        } else {
          tmp.coef <- rnorm(1, mean = coefMat.model.id.orig[tmp.ind, 1], sd = coefMat.model.id.orig[tmp.ind, 2])
        }

        initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]] <-
          c(initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]], tmp.coef)
      }

      for (Treat.name.i in 2:nof1$n.Treat) {

        tmp.ind  <- grepl(paste0("]", nof1$Treat.name[Treat.name.i]), rownames(coefMat.model.id.orig))
        # if there is no this treatment on this participant, generate grand mean
        # minus 1 because of there is a reference treatment
        if (sum(tmp.ind) == 0) {
          if (nof1$model.intcpt == "fixed") {
            tmp.coef <- initial.values[[n.chains.i]]$d[Treat.name.i - 1]
          } else {
            tmp.coef <- initial.values[[n.chains.i]]$b[Treat.name.i]
          }
        }  else if (is.nan(coefMat.model.id.orig[tmp.ind, 2])) { # no standard error
          if (nof1$model.intcpt == "fixed") {
            tmp.coef <- initial.values[[n.chains.i]]$d[Treat.name.i - 1]
          } else {
            tmp.coef <- initial.values[[n.chains.i]]$b[Treat.name.i]
          }

        }else {
          tmp.coef <- rnorm(1, mean = coefMat.model.id.orig[tmp.ind, 1], sd = coefMat.model.id.orig[tmp.ind, 2])
        }
        # if (length(tmp.coef) == 0) {
        #   tmp.coef <- coef.model.id.orig[grepl(paste0("]", nof1$Treat.name[Treat.name.i]), names(coef.model.id.orig))]
        # }

        initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <-
          c(initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]], tmp.coef)
      }

    } # for (id.i in 1:nof1$n.ID) {

    # sigmaSq_d
    if (nof1$model.intcpt == "fixed") {
      tmp.betas <- NULL
      initial.values[[n.chains.i]]$beta <- NULL
    } else if (nof1$model.intcpt == "random") {
      tmp.betas <- initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]]
      initial.values[[n.chains.i]]$beta <- initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]]
      initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]] <- NULL
    }

    for (Treat.name.i in 2:nof1$n.Treat) {
      tmp.betas <- c(tmp.betas, initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])

      # assign values to beta matrix and remove those treatment specific beta's
      initial.values[[n.chains.i]]$beta <- cbind(initial.values[[n.chains.i]][["beta", exact = T]], initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])
      initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
    }

    if (nof1$model.intcpt == "fixed") {
      initial.values[[n.chains.i]]$prec_delta <- 1 / var(tmp.betas)
    } else if (nof1$model.intcpt == "random") {
      initial.values[[n.chains.i]]$prec_beta <- 1 / var(tmp.betas)
    }

    # adjust for trend - splines
    if (nof1$spline.trend) {
      # eta
      initial.values[[n.chains.i]]$eta <- NULL
      for (eta.i in 1:nof1$spline_df) {

        ind.coefMat <- grepl(paste0("spline.long", eta.i), names.coef)
        initial.values[[n.chains.i]]$eta <- c(initial.values[[n.chains.i]]$eta,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta > 100] <- 100
      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta < -100] <- -100
    } # if (nof1$bs.trend) {

  } # n.chains.i

  return(initial.values)

}



nof1.ma.binomial.inits <- function(nof1, n.chains) {

  if (nof1$model.intcpt == "fixed") {
    tmp.formula <- "nof1$Y.long ~ -1 + nof1$ID + nof1$Treat"
  } else if (nof1$model.intcpt == "random") {
    tmp.formula <- "nof1$Y.long ~ -1 + nof1$Treat"
  }

  if (!is.null(nof1$names.covariates)) {
    for (covariates.i in 1:length(nof1$names.long.covariates)) {
      tmp.formula <- paste0(tmp.formula, " + nof1$", nof1$names.long.covariates[covariates.i])
    }
  }

  if (nof1$spline.trend) {
    for (bs_df.i in 1:nof1$spline_df) {
      tmp.formula <- paste0(tmp.formula, " + nof1$spline.long", bs_df.i)
    }
  }

  model <- glm(as.formula(tmp.formula), family = binomial(link = nof1$model.linkfunc))

  summ.model    <- summary(model)
  coefMat.model <- coef(summ.model)
  names.coef    <- rownames(coefMat.model)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]]   <- list()

    if (nof1$model.intcpt == "fixed") {
      # d
      initial.values[[n.chains.i]]$d <- NULL
      for (Treat.name.i in 2:nof1$n.Treat) {

        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.name.i]), names.coef)
        initial.values[[n.chains.i]]$d <- c(initial.values[[n.chains.i]]$d,
                                            rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      # alpha
      initial.values[[n.chains.i]]$alpha <- NULL
      for (id.i in 1:nof1$n.ID) {

        ind.coefMat <- grepl(nof1$uniq.ID[id.i], names.coef)
        initial.values[[n.chains.i]]$alpha <- c(initial.values[[n.chains.i]]$alpha,
                                                rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

      }

      initial.values[[n.chains.i]]$alpha[initial.values[[n.chains.i]]$alpha > 100] <- 100
      initial.values[[n.chains.i]]$alpha[initial.values[[n.chains.i]]$alpha < -100] <- -100

    } else if (nof1$model.intcpt == "random") { # if (nof1$model.intcpt == "fixed") {
      # b
      initial.values[[n.chains.i]]$b <- NULL
      for (n.Treat.i in 1:nof1$n.Treat) {
        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[n.Treat.i]), names.coef)
        initial.values[[n.chains.i]]$b <- c(initial.values[[n.chains.i]]$b,
                                            rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
      }
    }

    # do not generate alpha and beta's
    # values will become really unstable for individual participants and the model will not fit
    # beta's
    # beta and variance of beta are not settable
    for (id.i in 1:nof1$n.ID) {

      # beta
      ind.id        <- which((nof1$ID == nof1$uniq.ID[id.i]) & (!is.na(nof1$Y.long)))
      if (length(unique(nof1$Treat[ind.id])) == 1) {
        model.id.orig <- glm(nof1$Y.long[ind.id] ~ 1, family = binomial(link = nof1$model.linkfunc))
      } else {
        if (nof1$model.intcpt == "fixed") {
          model.id.orig <- glm(nof1$Y.long[ind.id] ~ nof1$Treat[ind.id], family = binomial(link = nof1$model.linkfunc))
        } else {
          model.id.orig <- glm(nof1$Y.long[ind.id] ~ -1 + nof1$Treat[ind.id], family = binomial(link = nof1$model.linkfunc))
        }
      }
      coefMat.model.id.orig <- coef(summary(model.id.orig))

      if (nof1$model.intcpt == "random") {
        tmp.ind  <- grepl(paste0("]", nof1$Treat.name[1]), rownames(coefMat.model.id.orig))
        # if there is no this treatment on this participant, generate grand mean
        # minus 1 because of there is a reference treatment
        if (sum(tmp.ind) == 0) {
          tmp.coef <- initial.values[[n.chains.i]]$b[1]
        } else if (is.nan(coefMat.model.id.orig[tmp.ind, 2])) { # no standard error
          tmp.coef <- initial.values[[n.chains.i]]$b[1]
        } else {
          tmp.coef <- rnorm(1, mean = coefMat.model.id.orig[tmp.ind, 1], sd = coefMat.model.id.orig[tmp.ind, 2])
        }

        initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]] <-
          c(initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]], tmp.coef)
      }

      for (Treat.name.i in 2:nof1$n.Treat) {

        tmp.ind  <- grepl(paste0("]", nof1$Treat.name[Treat.name.i]), rownames(coefMat.model.id.orig))
        # if there is no this treatment on this participant, generate grand mean
        # minus 1 because of there is a reference treatment
        if (sum(tmp.ind) == 0) {
          if (nof1$model.intcpt == "fixed") {
            tmp.coef <- initial.values[[n.chains.i]]$d[Treat.name.i - 1]
          } else {
            tmp.coef <- initial.values[[n.chains.i]]$b[Treat.name.i]
          }
        } else if (is.nan(coefMat.model.id.orig[tmp.ind, 2])) { # no standard error
          if (nof1$model.intcpt == "fixed") {
            tmp.coef <- initial.values[[n.chains.i]]$d[Treat.name.i - 1]
          } else {
            tmp.coef <- initial.values[[n.chains.i]]$b[Treat.name.i]
          }

        } else {
          tmp.coef <- rnorm(1, mean = coefMat.model.id.orig[tmp.ind, 1], sd = coefMat.model.id.orig[tmp.ind, 2])
        }
        # if (length(tmp.coef) == 0) {
        #   tmp.coef <- coef.model.id.orig[grepl(paste0("]", nof1$Treat.name[Treat.name.i]), names(coef.model.id.orig))]
        # }

        initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <-
          c(initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]], tmp.coef)
      }

    } # for (id.i in 1:nof1$n.ID) {

    # sigmaSq_d
    if (nof1$model.intcpt == "fixed") {
      tmp.betas <- NULL
      initial.values[[n.chains.i]]$beta <- NULL
    } else if (nof1$model.intcpt == "random") {
      tmp.betas <- initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]]
      initial.values[[n.chains.i]]$beta <- initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]]
      initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[1])]] <- NULL
    }

    for (Treat.name.i in 2:nof1$n.Treat) {
      tmp.betas <- c(tmp.betas, initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])

      # assign values to beta matrix and remove those treatment specific beta's
      initial.values[[n.chains.i]]$beta <- cbind(initial.values[[n.chains.i]][["beta", exact = T]], initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])
      initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
    }

    if (nof1$model.intcpt == "fixed") {
      initial.values[[n.chains.i]]$prec_delta <- 1 / var(tmp.betas)
    } else if (nof1$model.intcpt == "random") {
      initial.values[[n.chains.i]]$prec_beta <- 1 / var(tmp.betas)
    }

    # adjust for trend - splines
    if (nof1$spline.trend) {
      # eta
      initial.values[[n.chains.i]]$eta <- NULL
      for (eta.i in 1:nof1$spline_df) {

        ind.coefMat <- grepl(paste0("spline.long", eta.i), names.coef)
        initial.values[[n.chains.i]]$eta <- c(initial.values[[n.chains.i]]$eta,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta > 100] <- 100
      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta < -100] <- -100
    } # if (nof1$bs.trend) {

  } # n.chains.i

  return(initial.values)

}



nof1.ma.ordinal.inits <- function(nof1, n.chains) {

  Y_matrix <- NULL
  for (i in 1:nof1$ord.ncat) {
    Y_matrix <- cbind(Y_matrix,
                      as.numeric(nof1$Y.long == i))
  }

  # Treat.matrix <- NULL
  # for(i in nof1$Treat.name[1:length(nof1$Treat.name)]){
  #   Treat.matrix <- cbind(Treat.matrix, as.numeric(nof1$Treat == i))
  # }

  if (nof1$ord.model == "cumulative") {
    fit <- vglm(Y_matrix ~ nof1$Treat, cumulative(parallel = nof1$ord.parallel))
  } else {
    fit <- vglm(Y_matrix ~ nof1$Treat, acat(parallel = nof1$ord.parallel))
  }

  co <- coef(summary(fit))
  names.coef <- rownames(co)
  vcov.coef  <- vcov(fit)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]] <- list()

    # b
    initial.values[[n.chains.i]]$b <- rnorm(1, mean = co[names.coef == "(Intercept):1", 1], sd = co[names.coef == "(Intercept):1", 2])
    for (n.Treat.i in 2:nof1$n.Treat) {
      if (nof1$ord.parallel) {
        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[n.Treat.i]), names.coef)
      } else {
        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[n.Treat.i], ":1"), names.coef)
      }
      initial.values[[n.chains.i]]$b <- c(initial.values[[n.chains.i]]$b,
                                          rnorm(1,
                                                mean = co[ind.coefMat, 1] + co[names.coef == "(Intercept):1", 1],
                                                sd   = sqrt(vcov.coef[names.coef == "(Intercept):1", names.coef == "(Intercept):1"] +
                                                              vcov.coef[ind.coefMat, ind.coefMat] +
                                                              2 * vcov.coef[names.coef == "(Intercept):1", ind.coefMat])))
    }

    for (ncat.i in 2:(nof1$ord.ncat-1)) {
      tmp.b <- rnorm(1,
                     mean = co[names.coef == paste0("(Intercept):", ncat.i), 1],
                     sd = co[names.coef == paste0("(Intercept):", ncat.i), 2])

      if (nof1$ord.parallel) {
        tmp.b <- c(tmp.b, rep(NA, nof1$n.Treat - 1))
      } else {
        for (n.Treat.i in 2:nof1$n.Treat) {
          ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[n.Treat.i], ":", ncat.i), names.coef)
          tmp.b <- c(tmp.b,
                     rnorm(1,
                           mean = co[ind.coefMat, 1] + co[names.coef == paste0("(Intercept):", ncat.i), 1],
                           sd   = sqrt(vcov.coef[names.coef == paste0("(Intercept):", ncat.i), names.coef == paste0("(Intercept):", ncat.i)] +
                                         vcov.coef[ind.coefMat, ind.coefMat] +
                                         2 * vcov.coef[names.coef == paste0("(Intercept):", ncat.i), ind.coefMat])))
        }
      }
      initial.values[[n.chains.i]]$b <- rbind(initial.values[[n.chains.i]]$b, tmp.b)
    }

  } # n.chains.i

  return(initial.values)

}


# Network meta analysis
nof1.nma.inits <- function(nof1, n.chains) {

  if (nof1$response == "normal") {
    initial.values <- nof1.nma.normal.inits(nof1, n.chains)
  } else if (nof1$response == "poisson") {
    initial.values <- nof1.nma.poisson.inits(nof1, n.chains)
  } else if (nof1$response == "binomial") {
    initial.values <- nof1.nma.binomial.inits(nof1, n.chains)
  }
  # else if (nof1$response == "ordinal") {
  #   # if ordinal outcome, no initial value because no equivalent frequentist model can be fit
  #   initial.values <- nof1.ma.ordinal.inits(nof1, n.chains)
  # }

  # if (!is.null(initial.values)) {
  # if initial value generated is too big (i.e. model didn't work because of sparse data), just set inital value to be NULL
  if (any(abs(unlist(initial.values)) > 100, na.rm = T)) {
    initial.values <- NULL
  }
  # }

  return(initial.values)
}

nof1.nma.normal.inits <- function(nof1, n.chains) {

  if (nof1$model.intcpt == "fixed") {
    tmp.formula <- "nof1$data.long$Y ~ nof1$data.long$Treat"
  } else if (nof1$model.intcpt == "random") {
    tmp.formula <- "nof1$data.long$Y ~ -1 + nof1$data.long$Treat"
  }

  if (!is.null(nof1$cov.matrix)) {
    for (cov.i in 1:nof1$n.cov) {
      tmp.formula <- paste0(tmp.formula, " + nof1$data.long$lvl2.cov", cov.i, ":nof1$data.long$Treat")
    }
  }

  if (nof1$spline.trend) {
    tmp.formula <- paste0(tmp.formula, " + nof1$spline.long.matrix")
  }

  if (nof1$step.trend) {
    tmp.formula <- paste0(tmp.formula, " + nof1$step.long.matrix")
  }

  model <- glm(as.formula(tmp.formula), family = gaussian(link = nof1$model.linkfunc))

  summ.model    <- summary(model)
  coefMat.model <- coef(summ.model)
  names.coef    <- rownames(coefMat.model)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]] <- list()

    if (nof1$model.intcpt == "fixed") {
      # d
      initial.values[[n.chains.i]]$d <- NA
      for (Treat.i in 2:length(nof1$Treat.name)) {

        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.i], "$"), names.coef)
        initial.values[[n.chains.i]]$d <- c(initial.values[[n.chains.i]]$d,
                                            rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

    } else if (nof1$model.intcpt == "random") { # if (nof1$model.intcpt == "fixed") {
      # b
      initial.values[[n.chains.i]]$b <- NULL
      for (Treat.i in 1:length(nof1$Treat.name)) {
        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.i], "$"), names.coef)
        initial.values[[n.chains.i]]$b <- c(initial.values[[n.chains.i]]$b,
                                            rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
      }
    }

    # do not generate alpha and beta's and variance
    # values will become really unstable for individual participants and the model will not fit

    # adjust for trend - splines
    if (nof1$spline.trend) {
      # eta
      initial.values[[n.chains.i]]$eta <- NULL
      for (eta.i in 1:nof1$spline.df) {

        ind.coefMat <- grepl(paste0("spline.long.matrix", eta.i), names.coef)
        initial.values[[n.chains.i]]$eta <- c(initial.values[[n.chains.i]]$eta,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta > 100] <- 100
      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta < -100] <- -100
    } # if (nof1$bs.trend) {

    # adjust for trend - steps
    if (nof1$step.trend) {
      # eta
      initial.values[[n.chains.i]]$eta <- NULL
      for (eta.i in 1:nof1$step.df) {

        ind.coefMat <- grepl(paste0("step.long.matrix", eta.i), names.coef)
        initial.values[[n.chains.i]]$eta <- c(initial.values[[n.chains.i]]$eta,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta > 100] <- 100
      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta < -100] <- -100
    } # if (nof1$bs.trend) {

  } # n.chains.i

  return(initial.values)

}

nof1.nma.poisson.inits <- function(nof1, n.chains) {

  if (nof1$model.intcpt == "fixed") {
    tmp.formula <- "nof1$data.long$Y ~ nof1$data.long$Treat"
  }

  if (!is.null(nof1$cov.matrix)) {
    for (cov.i in 1:nof1$n.cov) {
      tmp.formula <- paste0(tmp.formula, " + nof1$data.long$cov", cov.i, ":nof1$data.long$Treat")
    }
  }

  if (nof1$spline.trend) {
    tmp.formula <- paste0(tmp.formula, " + nof1$spline.long.matrix")
  }

  model <- glm(as.formula(tmp.formula), family = poisson(link = nof1$model.linkfunc))

  summ.model    <- summary(model)
  coefMat.model <- coef(summ.model)
  names.coef    <- rownames(coefMat.model)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]]   <- list()

    if (nof1$model.intcpt == "fixed") {
      # d
      initial.values[[n.chains.i]]$d <- NA
      for (Treat.i in 2:length(nof1$Treat.name)) {

        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.i], "$"), names.coef)
        initial.values[[n.chains.i]]$d <- c(initial.values[[n.chains.i]]$d,
                                            rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

    }
    # else if (nof1$model.intcpt == "random") { # if (nof1$model.intcpt == "fixed") {
    #   # b
    #   initial.values[[n.chains.i]]$b <- NULL
    #   for (n.Treat.i in 1:nof1$n.Treat) {
    #     ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[n.Treat.i]), names.coef)
    #     initial.values[[n.chains.i]]$b <- c(initial.values[[n.chains.i]]$b,
    #                                         rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
    #   }
    # }

    # do not generate alpha and beta's and variance
    # values will become really unstable for individual participants and the model will not fit

    # adjust for trend - splines
    if (nof1$spline.trend) {
      # eta
      initial.values[[n.chains.i]]$eta <- NULL
      for (eta.i in 1:nof1$spline.df) {

        ind.coefMat <- grepl(paste0("spline.long.matrix", eta.i), names.coef)
        initial.values[[n.chains.i]]$eta <- c(initial.values[[n.chains.i]]$eta,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta > 100] <- 100
      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta < -100] <- -100
    } # if (nof1$bs.trend) {

  } # n.chains.i

  return(initial.values)

}

nof1.nma.binomial.inits <- function(nof1, n.chains) {

  if (nof1$model.intcpt == "fixed") {
    tmp.formula <- "nof1$data.long$Y ~ nof1$data.long$Treat"
  }
  # else if (nof1$model.intcpt == "random") {
  #   tmp.formula <- "nof1$Y.long ~ -1 + nof1$Treat"
  # }

  if (!is.null(nof1$cov.matrix)) {
    for (cov.i in 1:nof1$n.cov) {
      tmp.formula <- paste0(tmp.formula, " + nof1$data.long$cov", cov.i, ":nof1$data.long$Treat")
    }
  }

  if (nof1$spline.trend) {
      tmp.formula <- paste0(tmp.formula, " + nof1$spline.long.matrix")
  }

  model <- glm(as.formula(tmp.formula), family = binomial(link = nof1$model.linkfunc))

  summ.model    <- summary(model)
  coefMat.model <- coef(summ.model)
  names.coef    <- rownames(coefMat.model)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]]   <- list()

    if (nof1$model.intcpt == "fixed") {
      # d
      initial.values[[n.chains.i]]$d <- NA
      for (Treat.i in 2:length(nof1$Treat.name)) {

        ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.i], "$"), names.coef)
        initial.values[[n.chains.i]]$d <- c(initial.values[[n.chains.i]]$d,
                                            rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

    }
    # else if (nof1$model.intcpt == "random") { # if (nof1$model.intcpt == "fixed") {
    #   # b
    #   initial.values[[n.chains.i]]$b <- NULL
    #   for (n.Treat.i in 1:nof1$n.Treat) {
    #     ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[n.Treat.i]), names.coef)
    #     initial.values[[n.chains.i]]$b <- c(initial.values[[n.chains.i]]$b,
    #                                         rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
    #   }
    # }

    # do not generate alpha and beta's and variance
    # values will become really unstable for individual participants and the model will not fit

    # adjust for trend - splines
    if (nof1$spline.trend) {
      # eta
      initial.values[[n.chains.i]]$eta <- NULL
      for (eta.i in 1:nof1$spline.df) {

        ind.coefMat <- grepl(paste0("spline.long.matrix", eta.i), names.coef)
        initial.values[[n.chains.i]]$eta <- c(initial.values[[n.chains.i]]$eta,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))
        # initialize beta
        # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
      }

      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta > 100] <- 100
      initial.values[[n.chains.i]]$eta[initial.values[[n.chains.i]]$eta < -100] <- -100
    } # if (nof1$bs.trend) {

  } # n.chains.i

  return(initial.values)

}



