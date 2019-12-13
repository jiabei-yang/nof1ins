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

  with(nof1, {

  p <- rep(NA, ncat)
  c <- rep(NA, ncat-1)

  for (i in seq(ncat)) {
    p[i] = sum(Y[!is.na(Y)]==i)/nobs
    if (!is.na(p[i]) & p[i] == 0) { p[i] = 0.05 }
  }
  if(sum(p[!is.na(p)])> 1) {
    p[max(which(p == max(p)))] = p[max(which(p == max(p)))] + 1 - sum(p)
  }

  initial.values = list()
  for(i in 1:n.chains){
    initial.values[[i]] = list()
  }

  for(i in 1:n.chains){
    p <- combinat::rmultz2(nobs,p)/nobs
    if (any(p == 0)) {
      p[which(p == 0)] <- 0.05
      p[max(which(p == max(p)))] <- p[max(which(p == max(p)))] + 1 - sum(p)
    }

    q <- 1 - cumsum(p)
    for(j in seq(ncat -1)){
      c[j] <- -log(q[j]/(1-q[j]))
    }
    dc <- c(c[1], c[-1] - c[-(ncat-1)])
    initial.values[[i]][["dc"]] <- dc
  }


  Treat.matrix <- NULL
  for(i in Treat.name){
    Treat.matrix <- cbind(Treat.matrix, nof1[[paste0("Treat_", i)]])
  }

  # if(is.na(knots)){
  #   model <- polr(as.ordered(Y) ~ Treat.matrix, Hess = TRUE)
  #   co = coef(summary(model))
  # } 
  model <- MASS::polr(as.ordered(Y) ~ Treat.matrix, Hess = TRUE)
  co = coef(summary(model))
  
  if(!is.null(model)){
    co_Treat <- co[grep('Treat.matrix', rownames(coef(summary(model)))),,drop = FALSE]
    #co_BS <- co[grep('BS', rownames(coef(summary(model)))), ]

    for(i in 1:n.chains){
      for(j in 1:length(Treat.name)){
        initial.values[[i]][[paste0("beta_", Treat.name[j])]] <- co_Treat[j,1] + rnorm(1) * co_Treat[j,2]
      }

      # if(!is.na(knots)){
      #   for(j in 1:ncol(BS)){
      #     initial.values[[i]][[paste0("gamma", j)]] <- co_BS[j, 1] + rnorm(1) * co_BS[j, 2]
      #   }
      # }
    }
  }
  
  return(initial.values)
  })
}


