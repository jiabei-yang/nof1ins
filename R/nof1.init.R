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

    if(!is.nan(summary(model)$fstat[1])){
      for(i in 1:n.chains){

        fit <- summary(model)
        df <- fit$df[2]
        random.ISigma <- rchisq(1, df)
        resid.var <- fit$sigma^2
        sigma2 <- resid.var * df/ random.ISigma

        if(hy.prior[[1]] == "dunif"){
          if(sqrt(sigma2) > hy.prior[[3]]){
            stop("data has more variability than your prior does")
          }
        }

        if(hy.prior[[1]] == "dgamma"){
          initial.values[[i]][["prec"]] <- 1/sigma2
        } else if(hy.prior[[1]] == "dunif" || network$hy.prior[[1]] == "dhnorm"){
          initial.values[[i]][["sd"]] <- sqrt(sigma2)
        }
      }
    }
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
  }

  #if initial value generated is too big (i.e. model didn't work because of sparse data), just set inital value to be NULL
  if(any(abs(unlist(initial.values)) > 100)) {
    initial.values <- NULL
  }

  return(initial.values)
}



nof1.ma.normal.inits <- function(nof1, n.chains) {

  if (is.null(nof1$names.covariates)) {
    model <- lm(nof1$Y.long ~ -1 + nof1$ID + nof1$Treat)
  } else {
    tmp.formula <- "nof1$Y.long ~ -1 + nof1$ID + nof1$Treat"
    for (covariates.i in 1:length(nof1$names.long.covariates)) {
      tmp.formula <- paste0(tmp.formula, " + nof1$", nof1$names.long.covariates[covariates.i])
    }
    model <- lm(as.formula(tmp.formula))
  }

  summ.model    <- summary(model)
  coefMat.model <- coef(summ.model)
  names.coef    <- rownames(coefMat.model)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]]   <- list()

    # d
    initial.values[[n.chains.i]]$d <- NULL
    for (Treat.name.i in 2:nof1$n.Treat) {

      ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.name.i]), names.coef)
      initial.values[[n.chains.i]]$d <- c(initial.values[[n.chains.i]]$d,
                                          rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

      # initialize beta
      # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
    }

    # alpha, beta's
    # beta and variance of beta are not settable
    initial.values[[n.chains.i]]$alpha <- NULL
    for (id.i in 1:nof1$n.ID) {

      ind.coefMat <- grepl(nof1$uniq.ID[id.i], names.coef)
      initial.values[[n.chains.i]]$alpha <- c(initial.values[[n.chains.i]]$alpha,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

      # beta
      ind.id        <- which((nof1$ID == nof1$uniq.ID[id.i]) & (!is.na(nof1$Y.long)))
      if (length(unique(nof1$Treat[ind.id])) == 1) {
        model.id.orig <- lm(nof1$Y.long[ind.id] ~ 1)
      } else {
        model.id.orig <- lm(nof1$Y.long[ind.id] ~ nof1$Treat[ind.id])
      }
      coefMat.model.id.orig <- coef(summary(model.id.orig))

      for (Treat.name.i in 2:nof1$n.Treat) {

        tmp.ind  <- grepl(paste0("]", nof1$Treat.name[Treat.name.i]), rownames(coefMat.model.id.orig))
        # if there is no this treatment on this participant, generate grand mean
        # minus 1 because of there is a reference treatment
        if (sum(tmp.ind) == 0) {
          tmp.coef <- initial.values[[n.chains.i]]$d[Treat.name.i - 1]
        } else {
          tmp.coef <- rnorm(1, mean = coefMat.model.id.orig[tmp.ind, 1], sd = coefMat.model.id.orig[tmp.ind, 2])
        }
        # if (length(tmp.coef) == 0) {
        #   tmp.coef <- coef.model.id.orig[grepl(paste0("]", nof1$Treat.name[Treat.name.i]), names(coef.model.id.orig))]
        # }

        initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <-
          c(initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]], tmp.coef)
      }

    }

    # sigmaSq_d
    tmp.betas <- NULL
    initial.values[[n.chains.i]]$beta <- NULL
    for (Treat.name.i in 2:nof1$n.Treat) {
      tmp.betas <- c(tmp.betas, initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])

      # assign values to beta matrix and remove those treatment specific beta's
      initial.values[[n.chains.i]]$beta <- cbind(initial.values[[n.chains.i]]$beta, initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])
      initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
    }
    initial.values[[n.chains.i]]$prec_delta <- 1 / var(tmp.betas)

    # se_resid
    initial.values[[n.chains.i]]$prec_resid <- 1 / (summ.model$sigma)^2 * summ.model$df[2]/rchisq(1, df = summ.model$df[2])

  } # n.chains.i

  return(initial.values)

}



nof1.ma.poisson.inits <- function(nof1, n.chains) {

  if (is.null(nof1$names.covariates)) {
    model <- glm(nof1$Y.long ~ -1 + nof1$ID + nof1$Treat, family = poisson(link = "log"))
  } else {
    tmp.formula <- "nof1$Y.long ~ -1 + nof1$ID + nof1$Treat"
    for (covariates.i in 1:length(nof1$names.long.covariates)) {
      tmp.formula <- paste0(tmp.formula, " + nof1$", nof1$names.long.covariates[covariates.i])
    }
    model <- glm(as.formula(tmp.formula), family = poisson(link = "log"))
  }

  summ.model    <- summary(model)
  coefMat.model <- coef(summ.model)
  names.coef    <- rownames(coefMat.model)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]]   <- list()

    # d
    initial.values[[n.chains.i]]$d <- NULL
    for (Treat.name.i in 2:nof1$n.Treat) {

      ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.name.i]), names.coef)
      initial.values[[n.chains.i]]$d <- c(initial.values[[n.chains.i]]$d,
                                          rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

      # initialize beta
      # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
    }

    # alpha, beta's
    # beta and variance of beta are not settable
    initial.values[[n.chains.i]]$alpha <- NULL
    for (id.i in 1:nof1$n.ID) {

      ind.coefMat <- grepl(nof1$uniq.ID[id.i], names.coef)
      initial.values[[n.chains.i]]$alpha <- c(initial.values[[n.chains.i]]$alpha,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

      # beta
      ind.id        <- which((nof1$ID == nof1$uniq.ID[id.i]) & (!is.na(nof1$Y.long)))
      if (length(unique(nof1$Treat[ind.id])) == 1) {
        model.id.orig <- glm(nof1$Y.long[ind.id] ~ 1, family = poisson(link = "log"))
      } else {
        model.id.orig <- glm(nof1$Y.long[ind.id] ~ nof1$Treat[ind.id], family = poisson(link = "log"))
      }
      coefMat.model.id.orig <- coef(summary(model.id.orig))

      for (Treat.name.i in 2:nof1$n.Treat) {

        tmp.ind  <- grepl(paste0("]", nof1$Treat.name[Treat.name.i]), rownames(coefMat.model.id.orig))
        # if there is no this treatment on this participant, generate grand mean
        # minus 1 because of there is a reference treatment
        if (sum(tmp.ind) == 0) {
          tmp.coef <- initial.values[[n.chains.i]]$d[Treat.name.i - 1]
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
    tmp.betas <- NULL
    initial.values[[n.chains.i]]$beta <- NULL
    for (Treat.name.i in 2:nof1$n.Treat) {
      tmp.betas <- c(tmp.betas, initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])

      # assign values to beta matrix and remove those treatment specific beta's
      initial.values[[n.chains.i]]$beta <- cbind(initial.values[[n.chains.i]]$beta, initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])
      initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
    }
    initial.values[[n.chains.i]]$prec_delta <- 1 / var(tmp.betas)

  } # n.chains.i

  return(initial.values)

}



nof1.ma.binomial.inits <- function(nof1, n.chains) {

  if (is.null(nof1$names.covariates)) {
    model <- glm(nof1$Y.long ~ -1 + nof1$ID + nof1$Treat, family = binomial(link = "logit"))
  } else {
    tmp.formula <- "nof1$Y.long ~ -1 + nof1$ID + nof1$Treat"
    for (covariates.i in 1:length(nof1$names.long.covariates)) {
      tmp.formula <- paste0(tmp.formula, " + nof1$", nof1$names.long.covariates[covariates.i])
    }
    model <- glm(as.formula(tmp.formula), family = binomial(link = "logit"))
  }

  summ.model    <- summary(model)
  coefMat.model <- coef(summ.model)
  names.coef    <- rownames(coefMat.model)

  initial.values = list()
  for(n.chains.i in 1:n.chains){

    initial.values[[n.chains.i]]   <- list()

    # d
    initial.values[[n.chains.i]]$d <- NULL
    for (Treat.name.i in 2:nof1$n.Treat) {

      ind.coefMat <- grepl(paste0("Treat", nof1$Treat.name[Treat.name.i]), names.coef)
      initial.values[[n.chains.i]]$d <- c(initial.values[[n.chains.i]]$d,
                                          rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

      # initialize beta
      # initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
    }

    # do not generate alpha and beta's
    # values will become really unstable for individual participants and the model will not fit
    # alpha, beta's
    # beta and variance of beta are not settable
    initial.values[[n.chains.i]]$alpha <- NULL
    for (id.i in 1:nof1$n.ID) {

      ind.coefMat <- grepl(nof1$uniq.ID[id.i], names.coef)
      initial.values[[n.chains.i]]$alpha <- c(initial.values[[n.chains.i]]$alpha,
                                              rnorm(1, mean = coefMat.model[ind.coefMat, 1], sd = coefMat.model[ind.coefMat, 2]))

      # beta
      ind.id        <- which((nof1$ID == nof1$uniq.ID[id.i]) & (!is.na(nof1$Y.long)))
      if (length(unique(nof1$Treat[ind.id])) == 1) {
        model.id.orig <- glm(nof1$Y.long[ind.id] ~ 1, family = binomial(link = "logit"))
      } else {
        model.id.orig <- glm(nof1$Y.long[ind.id] ~ nof1$Treat[ind.id], family = binomial(link = "logit"))
      }
      coefMat.model.id.orig <- coef(summary(model.id.orig))

      for (Treat.name.i in 2:nof1$n.Treat) {

        tmp.ind  <- grepl(paste0("]", nof1$Treat.name[Treat.name.i]), rownames(coefMat.model.id.orig))
        # if there is no this treatment on this participant, generate grand mean
        # minus 1 because of there is a reference treatment
        if (sum(tmp.ind) == 0) {
          tmp.coef <- initial.values[[n.chains.i]]$d[Treat.name.i - 1]
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
    tmp.betas <- NULL
    initial.values[[n.chains.i]]$beta <- NULL
    for (Treat.name.i in 2:nof1$n.Treat) {
      tmp.betas <- c(tmp.betas, initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])

      # assign values to beta matrix and remove those treatment specific beta's
      initial.values[[n.chains.i]]$beta <- cbind(initial.values[[n.chains.i]]$beta, initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]])
      initial.values[[n.chains.i]][[paste0("beta_", nof1$Treat.name[Treat.name.i])]] <- NULL
    }
    initial.values[[n.chains.i]]$prec_delta <- 1 / var(tmp.betas)

  } # n.chains.i

  return(initial.values)

}
