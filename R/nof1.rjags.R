nof1.rjags <- function(nof1){

  response <- nof1$response
  code <- if(response == "normal"){
    nof1.normal.rjags(nof1)
  } else if(response == "ordinal"){
    nof1.ordinal.rjags(nof1)
  } else if(response == "binomial"){
    nof1.binomial.rjags(nof1)
  } else if(response == "poisson"){
    nof1.poisson.rjags(nof1)
  }
  return(code)
}

nof1.normal.rjags <- function(nof1){

  code <- paste0("model{")

  code <- paste0(code,
                 "\n\tfor (i in 1:", nof1$nobs, ") {")

  if (nof1$corr.y) {
    code <- paste0(code,
                   "\n\t\tY[i] ~ dnorm(m[i], prec*((1 - equals(i, 1)) * 1 + equals(i, 1) * (1-rho^2)))",
                   "\n\t\tm[i] <- mu[i] + rho * e[i]")
  } else {
    code <- paste0(code,
                   "\n\t\tY[i] ~ dnorm(m[i], prec)",
                   "\n\t\tm[i] <- mu[i]")
  }

  code <- paste0(code,
                 "\n\t\tmu[i] <- beta_", nof1$Treat.name[1], "*Treat_", nof1$Treat.name[1], "[i]")

  if (length(nof1$Treat.name) > 1){
    for(i in nof1$Treat.name[2:length(nof1$Treat.name)]){
      code <- paste0(code, " + beta_", i, "*Treat_", i, "[i]")
    }
  }

  if (nof1$bs.trend){
    for (i in 1:nof1$bs_df){
      code <- paste0(code, " + eta_", i, "*bs", i, "[i]")
    }
  }

  code <- paste0(code,
                 "\n\t}")

  if (nof1$corr.y) {
    code <- paste0(code,
                   "\n",
                   "\n\te[1] <- 0",
                   "\n\tfor (i in 2:", nof1$nobs, ") {",
                   "\n\t\te[i] <- (1 - equals(i, 1)) * (Y[i-1] - mu[i-1]) + equals(i, 1) * 0",
                   "\n\t}")

    code <- paste0(code,
                   "\n\trho ~ ", nof1$rho.prior[[1]], "(", nof1$rho.prior[[2]], ", ", nof1$rho.prior[[3]], ")")
  }

  for(i in nof1$Treat.name){
    code <- paste0(code,
                   "\n\tbeta_", i, " ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")")

    # Truncated normal distribution for beta's
    if (length(nof1$beta.prior) > 3) {

      if (!is.na(nof1$beta.prior[[4]])) {
        code <- paste0(code,
                       " T(", nof1$beta.prior[[4]], ",")
      } else {
        code <- paste0(code,
                       " T(,")
      }

      if(!is.na(nof1$beta.prior[[5]])) {
        code <- paste0(code,
                       nof1$beta.prior[[5]], ")")
      } else {
        code <- paste0(code,
                       ")")
      }
    } # if (length(beta.prior) > 3) {
  } #  for(i in Treat.name){

  if (nof1$bs.trend){
    for (i in 1:nof1$bs_df){
      code <- paste0(code, "\n\teta_", i, " ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")")
    }
  }

  code <- paste0(code, nof1.hy.prior.rjags(nof1$hy.prior), "\n}")
  return(code)
}

nof1.binomial.rjags <- function(nof1){

  with(nof1, {

    code <- paste0("model{")
    code <- paste0(code,
                   "\n\tfor (i in 1:", nobs, ") {",
                   "\n\t\tlogit(p[i]) <- beta_", Treat.name[1], "*Treat_", Treat.name[1], "[i]")

    if (length(Treat.name) > 1){
      for(i in Treat.name[2:length(Treat.name)]){
        code <- paste0(code, " + beta_", i, "*Treat_", i, "[i]")
      }
    }

    # if(!is.null(knots)){
    #   for(j in 1:ncol(BS)){
    #     code <- paste0(code, " + gamma", j, "* BS[i,", j, "]")
    #   }
    # }

    code <- paste0(code,
                   "\n\t\tY[i] ~ dbern(p[i])",
                   "\n\t}")
    # "\n\talpha ~ ", alpha.prior[[1]], "(", alpha.prior[[2]], ",", alpha.prior[[3]], ")")

    for(i in Treat.name){
      code <- paste0(code, "\n\tbeta_", i, " ~ ", beta.prior[[1]], "(", beta.prior[[2]], ",", beta.prior[[3]], ")")
    }

    # if(!is.null(knots)){
    #   for(j in 1:ncol(BS)){
    #     code <- paste0(code, "\n\tgamma", j, " ~ ", gamma.prior[[1]], "(", gamma.prior[[2]], ",", gamma.prior[[3]], ")")
    #   }
    # }
    code <- paste0(code, "\n}")

    return(code)
  })

}

nof1.poisson.rjags <- function(nof1){

  with(nof1, {

    code <- paste0("model{")
    code <- paste0(code,
                   "\n\tfor (i in 1:", nobs, ") {",
                   "\n\t\tlog(lambda[i]) <- beta_", Treat.name[1], "*Treat_", Treat.name[1], "[i]")

    if (length(Treat.name) > 1){
      for(i in Treat.name[2:length(Treat.name)]){
        code <- paste0(code, " + beta_", i, "*Treat_", i, "[i]")
      }
    }

    if (bs.trend){
      for (i in 1:bs_df){
        code <- paste0(code, " + eta_", i, "*bs", i, "[i]")
      }
    }

    code <- paste0(code,
                   "\n\t\tY[i] ~ dpois(lambda[i])",
                   "\n\t}")
    #               "\n\talpha ~ ", alpha.prior[[1]], "(", alpha.prior[[2]], ",", alpha.prior[[3]], ")")

    for(i in Treat.name){
      code <- paste0(code, "\n\tbeta_", i, " ~ ", beta.prior[[1]], "(", beta.prior[[2]], ",", beta.prior[[3]], ")")
    }

    if (bs.trend){
      for (i in 1:bs_df){
        code <- paste0(code, "\n\teta_", i, " ~ ", eta.prior[[1]], "(", eta.prior[[2]], ", ", eta.prior[[3]], ")")
      }
    }

    code <- paste0(code, "\n}")
    return(code)
  })
}

nof1.ordinal.rjags <- function(nof1){

  code <- paste0("model{\n")
  code <- paste0(code,
                 "\tfor (i in 1:", nof1$nobs, ") {\n")

  code <- paste0(code,
                 "\t\tfor (j in 2:", nof1$ncat, ") {\n")

  if (nof1$ord.model == "cumulative") {

    code <- paste0(code,
                   "\t\t\tlogit(q[i, j-1]) <- alpha[j-1]")
    for(Treat.name.i in 2:nof1$n.Treat){
      code <- paste0(code,
                     " + beta_", nof1$Treat.name[Treat.name.i], " * Treat_", nof1$Treat.name[Treat.name.i], "[i]")
    }
    code <- paste0(code, "\n")
    code <- paste0(code,
                   "\t\t}\n")

    code <- paste0(code,
                   "\t\tp[i, 1] <- q[i, 1]\n")
    code <- paste0(code,
                   "\t\tfor (j in 2:", nof1$ncat - 1, ") {\n")
    code <- paste0(code,
                   "\t\t\tp[i,j] <- q[i, j] - q[i, j-1]\n")
    code <- paste0(code,
                   "\t\t}\n")
    code <- paste0(code,
                   "\t\tp[i, ", nof1$ncat, "] <- 1 - q[i, ", nof1$ncat - 1, "]\n")

  } else { # (nof1$ord.model == "acat")
    code <- paste0(code,
                   "\t\t\tp[i, j] <- p[i, j-1] * c[i, j-1]\n")
    code <- paste0(code,
                   "\t\t\tlog(c[i, j-1]) <- alpha[j-1]")
    for(Treat.name.i in 2:nof1$n.Treat){
      code <- paste0(code,
                     " + beta_", nof1$Treat.name[Treat.name.i], " * Treat_", nof1$Treat.name[Treat.name.i], "[i]")
    }
    code <- paste0(code, "\n")
    code <- paste0(code,
                   "\t\t}\n")

    # at least 3 categories in the ordinal outcome because we always retain the missing levels
    # so at least 2 contrasts
    code <- paste0(code,
                   "\t\tp[i, 1] <- 1 / (1 + c[i, 1] + c[i, 1] * c[i, 2]")

    # need to test if there are more than 5 levels
    if (nof1$ncat >= 4) {
      for (i in 4:nof1$ncat) {
        code <- paste0(code,
                       " + c[i, 1]")
        for (j in 2:(i - 1)) {
          code <- paste0(code,
                         " * c[i, ", j, "]")
        }
      }
    }
    code <- paste0(code,
                   ")\n")
  }


  code <- paste0(code,
                 "\t\tY[i] ~ dcat(p[i, ])\n")
  code <- paste0(code,
                 "\t}\n\n")

  # "\n\t\tp[i,1] <- 1 - Q[i,1]",
  # "\n\t\tfor(r in 2:", ncat-1, ") {",
  # "\n\t\t\tp[i,r] <- Q[i,r-1] - Q[i,r]",
  # "\n\t\t}",
  # "\n\t\tp[i,", ncat, "] <- Q[i,", ncat-1, "]",
  # "\n\t\tfor(r in 1:", ncat-1, ") {",
  # "\n\t\t\tlogit(Q[i,r]) <- -c[r]") # + epsilon[i]")

  # if(!is.null(knots)){
  #   for(j in 1:ncol(BS)){
  #     code <- paste0(code, " + gamma", j, "* BS[i,", j, "]")
  #   }
  # }

  # priors
  if (nof1$ord.model == "cumulative") {

    # dalpha put increasing restriction on alpha
    code <- paste0(code,
                   "\talpha[1] <- dalpha[1]\n",
                   "\tdalpha[1] ~ ", nof1$alpha.prior[[1]], "(", nof1$alpha.prior[[2]], ", ", nof1$alpha.prior[[3]], ")\n")

    code <- paste0(code,
                   "\tfor (j in 2:", nof1$ncat - 1,") {\n")
    code <- paste0(code,
                   "\t\talpha[j] <- alpha[j-1] + dalpha[j]\n")
    code <- paste0(code,
                   "\t\tdalpha[j] ~ ", nof1$alpha.prior[[1]], "(", nof1$alpha.prior[[2]], ", ", nof1$alpha.prior[[3]], ") T(0, )\n")
    code <- paste0(code,
                   "\t}\n")

  } else { # nof1$ord.model == "acat"
    code <- paste0(code,
                   "\tfor (j in 2:", nof1$ncat, ") {\n")
    code <- paste0(code,
                   "\t\talpha[j-1] ~ ", nof1$alpha.prior[[1]], "(", nof1$alpha.prior[[2]], ", ", nof1$alpha.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n\n")
  }

  for(Treat.name.i in 2:nof1$n.Treat){
    code <- paste0(code,
                   "\tbeta_", nof1$Treat.name[Treat.name.i], " ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
  }

  # code <- paste0(code,
  #                "\n\t\t}",
  #                "\n\t}",
  #                "\n\tfor(i in 2:", ncat-1, ") {",
  #                "\n\t\tdc[i] ~ ", dc.prior[[1]], "(", dc.prior[[2]], ",", dc.prior[[3]], ")",
  #                "\n\t}",
  #                "\n\tc[1] <- dc[1]",
  #                "\n\tfor (i in 2:", ncat-1, ") {",
  #                "\n\t\tc[i] <- c[i-1] + dc[i]",
  #                "\n\t}",
  #                "\n\tdc[1] ~ ", c1.prior[[1]], "(", c1.prior[[2]], ",", c1.prior[[3]], ")")
  #
  # for(i in Treat.name){
  #   code <- paste0(code, "\n\tbeta_", i, " ~ ", beta.prior[[1]], "(", beta.prior[[2]], ",", beta.prior[[3]], ")")
  # }
  # if(!is.null(knots)){
  #   for(j in 1:ncol(BS)){
  #     code <- paste0(code, "\n\tgamma", j, " ~ ", gamma.prior[[1]], "(", gamma.prior[[2]], ",", gamma.prior[[3]], ")")
  #   }
  # }

  #code <- paste0(code, nof1.hy.prior.rjags(hy.prior), "\n}")
  code <- paste0(code, "}")
  return(code)

}

nof1.hy.prior.rjags <- function(hy.prior){

  code <- ""
  distr <- hy.prior[[1]]
  if (distr == "dunif") {
    code <- paste0(code,
                   "\n\tsd <- pow(prec, -0.5)",
                   "\n\tprec ~ ", hy.prior[[1]], "(", hy.prior[[2]], ", ", hy.prior[[3]], ")",
                   "\n\tlogprec <- log(prec)")
  } else if(distr == "dgamma"){
    code <- paste0(code,
                   "\n\tsd <- pow(prec, -0.5)",
                   "\n\tprec ~ dgamma(", hy.prior[[2]], ", ", hy.prior[[3]], ")",
                   "\n\tlogprec <- log(prec)")
  } else if(distr == "dhnorm"){
    code <- paste0(code,
                   "\n\tsd ~ dnorm(", hy.prior[[2]], ", ", hy.prior[[3]], ")T(0,)",
                   "\n\tprec <- pow(sd, -2)",
                   "\n\tlogprec <- log(prec)")
  }
  return(code)
}

#################################
######### Meta analysis #########
#################################
nof1.ma.rjags <- function(nof1) {

  if (nof1$response == "normal") {
    code <- nof1.ma.normal.rjags(nof1)
  } else if (nof1$response == "poisson") {
    code <- nof1.ma.poisson.rjags(nof1)
  } else if (nof1$response == "binomial") {
    code <- nof1.ma.binomial.rjags(nof1)
  } else if (nof1$response == "ordinal") {
    code <- nof1.ma.ordinal.rjags(nof1)
  }


  return(code)
}

nof1.ma.normal.rjags <- function(nof1) {

  code <- paste0("model{\n")

  code <- paste0(code,
                 "\tfor (i in 1:n.ID) {\n\n")

  if ((nof1$model.intcpt == "fixed") & (nof1$model.slp == "random")) {
    # Fixed intercepts (prior)
    code <- paste0(code,
                   "\t\talpha[i] ~ ", nof1$alpha.prior[[1]], "(", nof1$alpha.prior[[2]], ", ", nof1$alpha.prior[[3]], ")")

    # Truncated normal distribution for beta's
    if (length(nof1$alpha.prior) > 3) {

      if (!is.na(nof1$alpha.prior[[4]])) {
        code <- paste0(code,
                       " T(", nof1$alpha.prior[[4]], ",")
      } else {
        code <- paste0(code,
                       " T(,")
      }

      if(!is.na(nof1$alpha.prior[[5]])) {
        code <- paste0(code,
                       nof1$alpha.prior[[5]], ")")
      } else {
        code <- paste0(code,
                       ")")
      }
    } # if (length(nof1$alpha.prior) > 3) {

    # Random effects
    code <- paste0(code,
                   "\n\t\tbeta[i, 1:", nof1$n.Treat - 1, "] ~ dmnorm.vcov(d[1:", nof1$n.Treat - 1, "], Sigma_delta[1:",
                   nof1$n.Treat - 1, ", 1:", nof1$n.Treat - 1, "])\n")

    for (Treat.name.i in 2:nof1$n.Treat) {
      code <- paste0(code,
                     "\t\tbeta_", nof1$Treat.name[Treat.name.i], "[i] <- beta[i, ", Treat.name.i - 1, "]\n")
    }


  } else if ((nof1$model.intcpt == "random") & (nof1$model.slp == "random")) { # (nof1$model.intcpt == "random")
    code <- paste0(code,
                   "\t\tbeta[i, 1:n.Treat] ~ dmnorm.vcov(b[1:n.Treat], Sigma_beta[1:n.Treat, 1:n.Treat])\n")
    for (n.Treat.i in 1:nof1$n.Treat) {
      code <- paste0(code,
                     "\t\tbeta_", nof1$Treat.name[n.Treat.i], "[i] <- beta[i, ", n.Treat.i, "]\n")
    }
  }
  code <- paste0(code, "\n")

  # First level model
  code <- paste0(code,
                 "\t\tfor (j in 1:nobs.ID[i]) {\n")

  # if (nof1$model.linkfunc == "log") {
  #   code <- paste0(code,
  #                  "\t\t\tY[j, i] ~ dnorm(mu[j, i], prec_resid) T(0, )\n")
  # } else {
  code <- paste0(code,
                 "\t\t\tY[j, i] ~ dnorm(mu[j, i], prec_resid)\n")
  # }

  # not common intercept
  if (nof1$model.intcpt != "common") {

    # fixed intercept, identity link
    if ((nof1$model.intcpt == "fixed")) {
      code <- paste0(code,
                     "\t\t\tmu[j, i] <- alpha[i]")

      # random intercept, identity link
    } else if (nof1$model.intcpt == "random") {
      code <- paste0(code,
                     "\t\t\tmu[j, i] <- beta_", nof1$Treat.name[1], "[i] * Treat_", nof1$Treat.name[1], "[j, i]")
    }

    # random slopes
    if (nof1$n.Treat > 1) {
      for(Treat.name.i in 2:nof1$n.Treat){
        code <- paste0(code,
                       " + beta_", nof1$Treat.name[Treat.name.i], "[i] * Treat_", nof1$Treat.name[Treat.name.i], "[j, i]")
      }
    }
    # cases not working yet: fixed and common slopes; other link functions

    # common intercept
  } else { # if (nof1$model.intcpt != "common") {

    # log link
    if (nof1$model.linkfunc == "log") {
      code <- paste0(code,
                     "\t\t\tmu[j, i] <- exp(lmu[j, i])\n")
      # identity link
    } else {
      code <- paste0(code,
                     "\t\t\tmu[j, i] <- lmu[j, i]\n")
    }

    # common intercept
    code <- paste0(code,
                   "\t\t\tlmu[j, i] <- beta_", nof1$Treat.name[1], " * Treat_", nof1$Treat.name[1], "[j, i]")

    # common slopes
    if (nof1$n.Treat > 1) {
      for(Treat.name.i in 2:nof1$n.Treat){
        code <- paste0(code,
                       " + beta_", nof1$Treat.name[Treat.name.i], " * Treat_", nof1$Treat.name[Treat.name.i], "[j, i]")
      }
    }
    # cases not working yet: fixed and random slopes;
  }


  # adjust for covariates
  if (!is.null(nof1$names.covariates)) {

    for (covariates.i in 1:length(nof1$names.covariates)) {
      code <- paste0(code,
                     " + beta_", nof1$names.covariates[covariates.i], " * ", nof1$names.covariates[covariates.i], "[j, i]")
    }

  }

  # adjust for trend by basis splines
  if (nof1$spline.trend) {
    for (l in 1:nof1$spline_df){
      code <- paste0(code, " + eta[", l, "]*spline", l, "[time[j, i]]")
    }
  }

  # adjust for trend by step functions
  if (nof1$step.trend) {
    for (p in nof1$n.steps){
      code <- paste0(code, " + eta[", p-1, "]*step", p, "[j, i]")
    }
  }

  code <- paste0(code, "\n")

  code <- paste0(code, "\t\t}\n") # for (j in 1:nobs.ID[i])
  code <- paste0(code, "\t}\n\n")  # for (i in 1:n.ID)

  # Priors
  # prior for residual error
  code <- paste0(code,
                 "\tsd_resid <- pow(prec_resid, -0.5)\n")
  code <- paste0(code,
                 "\tprec_resid ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n\n")

  # prior for treatment effect
  if ((nof1$model.intcpt == "fixed") & (nof1$model.slp == "random")) {
    code <- paste0(code,
                   "\tprec_delta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n")
    code <- paste0(code,
                   "\tsigmaSq_d <- pow(prec_delta, -1)\n\n")

    for (Treat.name.i in 2:nof1$n.Treat) {
      code <- paste0(code,
                     "\td[", Treat.name.i-1, "] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
      code <- paste0(code,
                     "\tSigma_delta[", Treat.name.i-1, ", ", Treat.name.i-1, "] <- sigmaSq_d\n")

      if ((Treat.name.i-1) > 1) {
        for (Treat.name.j in 1:(Treat.name.i-1-1)) {
          code <- paste0(code,
                         "\tSigma_delta[", Treat.name.i-1, ", ", Treat.name.j, "] <- sigmaSq_d / 2\n")
        }
      }

      if (Treat.name.i <= (nof1$n.Treat - 1)) {
        for (Treat.name.j in Treat.name.i:(nof1$n.Treat-1)) {
          code <- paste0(code,
                         "\tSigma_delta[", Treat.name.i-1, ", ", Treat.name.j, "] <- sigmaSq_d / 2\n")
        }
      }
    }

  } else if ((nof1$model.intcpt == "random") & (nof1$model.slp == "random")) { # if (nof1$model.intcpt == "fixed") {

    code <- paste0(code,
                   "\tprec_beta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n")
    code <- paste0(code,
                   "\tsigmaSq_beta <- pow(prec_beta, -1)\n")
    code <- paste0(code,
                   "\trho ~ ", nof1$rho.prior[[1]], "(", nof1$rho.prior[[2]], ", ", nof1$rho.prior[[3]], ")\n\n")

    code <- paste0(code,
                   "\tfor (k in 1:n.Treat) {\n")
    code <- paste0(code,
                   "\t\tb[k] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n\n") # for (k in 1:n.Treat) {

    # row 1
    code <- paste0(code,
                   "\tSigma_beta[1, 1] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\tfor (k in 2:n.Treat) {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[1, k] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t}\n") # for (k in 2:n.Treat) {
    # row 2:(n.Treat-1)
    code <- paste0(code,
                   "\tfor (k in 2:(n.Treat - 1)) {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[k, k] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\t\tfor (l in 1:(k-1)) {\n")
    code <- paste0(code,
                   "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t\t}\n") # for (l in 1:(k-1)) {
    code <- paste0(code,
                   "\t\tfor (l in (k+1):n.Treat) {\n")
    code <- paste0(code,
                   "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t\t}\n") # for (l in (k+1):n.Treat) {
    code <- paste0(code,
                   "\t}\n") # for (k in 2:(n.Treat - 1)) {
    # row n.Treat
    code <- paste0(code,
                   "\tSigma_beta[n.Treat, n.Treat] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\tfor (k in 1:(n.Treat - 1)) {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[n.Treat, k] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t}\n") # for (k in 1:(n.Treat - 1)) {

  } else if ((nof1$model.intcpt == "common") & (nof1$model.slp == "common")) {

    for(i in nof1$Treat.name){
      code <- paste0(code,
                     "\tbeta_", i, " ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")")

      # Truncated normal distribution for beta's
      if (length(nof1$beta.prior) > 3) {

        if (!is.na(nof1$beta.prior[[4]])) {
          code <- paste0(code,
                         " T(", nof1$beta.prior[[4]], ", ")
        } else {
          code <- paste0(code,
                         " T(, ")
        }

        if(!is.na(nof1$beta.prior[[5]])) {
          code <- paste0(code,
                         nof1$beta.prior[[5]], ")")
        } else {
          code <- paste0(code,
                         ")")
        }
      } # if (length(beta.prior) > 3) {

      code <- paste0(code, "\n")
    } #  for(i in Treat.name){
  }
  code <- paste0(code, "\n")

  # prior for covariate coefficients
  if (!is.null(nof1$names.covariates)) {

    for (covariates.i in 1:length(nof1$names.covariates)) {
      code <- paste0(code,
                     "\tbeta_", nof1$names.covariates[covariates.i], " ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    }
    code <- paste0(code, "\n")
  }

  # prior for trend coefficients - splines
  if (nof1$spline.trend) {
    code <- paste0(code,
                   "\tfor (l in 1:", nof1$spline_df, ") {\n")
    code <- paste0(code, "\t\teta[l] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n")
  }

  # prior for trend coefficients - steps
  if (nof1$step.trend) {
    code <- paste0(code,
                   "\tfor (p in 1:", max(nof1$n.steps) - 1, ") {\n")
    code <- paste0(code, "\t\teta[p] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n")
  }

  code <- paste0(code, "}") # model {

  # cat(code)
  return(code)
}

nof1.ma.poisson.rjags <- function(nof1) {

  code <- paste0("model{\n")

  code <- paste0(code,
                 "\tfor (i in 1:n.ID) {\n\n")

  if (nof1$model.intcpt == "fixed") {
    # Fixed intercepts (prior)
    code <- paste0(code,
                   "\t\talpha[i] ~ ", nof1$alpha.prior[[1]], "(", nof1$alpha.prior[[2]], ", ", nof1$alpha.prior[[3]], ")\n")

    # Random effects
    code <- paste0(code,
                   "\t\tbeta[i, 1:", nof1$n.Treat - 1, "] ~ dmnorm.vcov(d[1:", nof1$n.Treat - 1, "], Sigma_delta[1:",
                   nof1$n.Treat - 1, ", 1:", nof1$n.Treat - 1, "])\n")

    for (Treat.name.i in 2:nof1$n.Treat) {
      code <- paste0(code,
                     "\t\tbeta_", nof1$Treat.name[Treat.name.i], "[i] <- beta[i, ", Treat.name.i - 1, "]\n")
    }
  } else { # nof1$model.intcpt == "random"
    code <- paste0(code,
                   "\t\tbeta[i, 1:n.Treat] ~ dmnorm.vcov(b[1:n.Treat], Sigma_beta[1:n.Treat, 1:n.Treat])\n")
    for (n.Treat.i in 1:nof1$n.Treat) {
      code <- paste0(code,
                     "\t\tbeta_", nof1$Treat.name[n.Treat.i], "[i] <- beta[i, ", n.Treat.i, "]\n")
    }
  }
  code <- paste0(code, "\n")

  # First level model
  code <- paste0(code,
                 "\t\tfor (j in 1:nobs.ID[i]) {\n")

  code <- paste0(code,
                 "\t\t\tY[j, i] ~ dpois(lambda[j, i])\n")

  if (nof1$model.intcpt == "fixed") {
    code <- paste0(code,
                   "\t\t\tlog(lambda[j, i]) <- alpha[i]")
  } else {
    code <- paste0(code,
                   "\t\t\tlog(lambda[j, i]) <- beta_", nof1$Treat.name[1], "[i] * Treat_", nof1$Treat.name[1], "[j, i]")
  }

  if (nof1$n.Treat > 1) {
    for(Treat.name.i in 2:nof1$n.Treat){
      code <- paste0(code,
                     " + beta_", nof1$Treat.name[Treat.name.i], "[i] * Treat_", nof1$Treat.name[Treat.name.i], "[j, i]")
    }
  }

  # adjust for covariates
  if (!is.null(nof1$names.covariates)) {

    for (covariates.i in 1:length(nof1$names.covariates)) {
      code <- paste0(code,
                     " + beta_", nof1$names.covariates[covariates.i], " * ", nof1$names.covariates[covariates.i], "[j, i]")
    }

  }

  # adjust for trend by basis splines
  if (nof1$spline.trend) {
    for (l in 1:nof1$spline_df){
      code <- paste0(code, " + eta[", l, "]*spline", l, "[time[j, i]]")
    }
  }
  code <- paste0(code, "\n")

  code <- paste0(code, "\t\t}\n") # for (j in 1:nobs.ID[i])
  code <- paste0(code, "\t}\n\n")  # for (i in 1:n.ID)

  # Priors
  if (nof1$model.intcpt == "fixed") {
    # prior for treatment effect
    code <- paste0(code,
                   "\tprec_delta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n")
    code <- paste0(code,
                   "\tsigmaSq_d <- pow(prec_delta, -1)\n\n")

    for (Treat.name.i in 2:nof1$n.Treat) {
      code <- paste0(code,
                     "\td[", Treat.name.i-1, "] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
      code <- paste0(code,
                     "\tSigma_delta[", Treat.name.i-1, ", ", Treat.name.i-1, "] <- sigmaSq_d\n")

      if ((Treat.name.i-1) > 1) {
        for (Treat.name.j in 1:(Treat.name.i-1-1)) {
          code <- paste0(code,
                         "\tSigma_delta[", Treat.name.i-1, ", ", Treat.name.j, "] <- sigmaSq_d / 2\n")
        }
      }

      if (Treat.name.i <= (nof1$n.Treat - 1)) {
        for (Treat.name.j in Treat.name.i:(nof1$n.Treat-1)) {
          code <- paste0(code,
                         "\tSigma_delta[", Treat.name.i-1, ", ", Treat.name.j, "] <- sigmaSq_d / 2\n")
        }
      }
    }
  } else { # if (nof1$model.intcpt == "fixed") {
    code <- paste0(code,
                   "\tprec_beta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n")
    code <- paste0(code,
                   "\tsigmaSq_beta <- pow(prec_beta, -1)\n")
    code <- paste0(code,
                   "\trho ~ ", nof1$rho.prior[[1]], "(", nof1$rho.prior[[2]], ", ", nof1$rho.prior[[3]], ")\n\n")

    code <- paste0(code,
                   "\tfor (k in 1:n.Treat) {\n")
    code <- paste0(code,
                   "\t\tb[k] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n\n") # for (k in 1:n.Treat) {

    # row 1
    code <- paste0(code,
                   "\tSigma_beta[1, 1] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\tfor (k in 2:n.Treat) {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[1, k] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t}\n") # for (k in 2:n.Treat) {
    # row 2:(n.Treat-1)
    code <- paste0(code,
                   "\tfor (k in 2:(n.Treat - 1)) {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[k, k] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\t\tfor (l in 1:(k-1)) {\n")
    code <- paste0(code,
                   "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t\t}\n") # for (l in 1:(k-1)) {
    code <- paste0(code,
                   "\t\tfor (l in (k+1):n.Treat) {\n")
    code <- paste0(code,
                   "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t\t}\n") # for (l in (k+1):n.Treat) {
    code <- paste0(code,
                   "\t}\n") # for (k in 2:(n.Treat - 1)) {
    # row n.Treat
    code <- paste0(code,
                   "\tSigma_beta[n.Treat, n.Treat] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\tfor (k in 1:(n.Treat - 1)) {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[n.Treat, k] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t}\n") # for (k in 1:(n.Treat - 1)) {

  }
  code <- paste0(code, "\n")

  # prior for covariate coefficients
  if (!is.null(nof1$names.covariates)) {

    for (covariates.i in 1:length(nof1$names.covariates)) {
      code <- paste0(code,
                     "\tbeta_", nof1$names.covariates[covariates.i], " ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    }
    code <- paste0(code, "\n")
  }

  # prior for trend coefficients - splines
  if (nof1$spline.trend) {
    code <- paste0(code,
                   "\tfor (l in 1:", nof1$spline_df, ") {\n")
    code <- paste0(code, "\t\teta[l] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n")
  }

  code <- paste0(code, "}") # model {

  # cat(code)
  return(code)
}

nof1.ma.binomial.rjags <- function(nof1) {

  code <- paste0("model{\n")

  code <- paste0(code,
                 "\tfor (i in 1:n.ID) {\n\n")

  if (nof1$model.intcpt == "fixed") {
    # Fixed intercepts (prior)
    code <- paste0(code,
                   "\t\talpha[i] ~ ", nof1$alpha.prior[[1]], "(", nof1$alpha.prior[[2]], ", ", nof1$alpha.prior[[3]], ")")

    # Truncated normal distribution for beta's
    if (length(nof1$alpha.prior) > 3) {

      if (!is.na(nof1$alpha.prior[[4]])) {
        code <- paste0(code,
                       " T(", nof1$alpha.prior[[4]], ", ")
      } else {
        code <- paste0(code,
                       " T(, ")
      }

      if(!is.na(nof1$alpha.prior[[5]])) {
        code <- paste0(code,
                       nof1$alpha.prior[[5]], ")")
      } else {
        code <- paste0(code,
                       ")")
      }
    } # if (length(nof1$alpha.prior) > 3) {
    code <- paste0(code, "\n")

    # Random effects
    code <- paste0(code,
                   "\t\tbeta[i, 1:", nof1$n.Treat - 1, "] ~ dmnorm.vcov(d[1:", nof1$n.Treat - 1, "], Sigma_delta[1:",
                   nof1$n.Treat - 1, ", 1:", nof1$n.Treat - 1, "])\n")

    for (Treat.name.i in 2:nof1$n.Treat) {
      code <- paste0(code,
                     "\t\tbeta_", nof1$Treat.name[Treat.name.i], "[i] <- beta[i, ", Treat.name.i - 1, "]\n")
    }

  } else { # nof1$model.intcpt == "random"
    code <- paste0(code,
                   "\t\tbeta[i, 1:n.Treat] ~ dmnorm.vcov(b[1:n.Treat], Sigma_beta[1:n.Treat, 1:n.Treat])\n")
    for (n.Treat.i in 1:nof1$n.Treat) {
      code <- paste0(code,
                     "\t\tbeta_", nof1$Treat.name[n.Treat.i], "[i] <- beta[i, ", n.Treat.i, "]\n")
    }
  }
  code <- paste0(code, "\n")

  # First level model
  code <- paste0(code,
                 "\t\tfor (j in 1:nobs.ID[i]) {\n")

  code <- paste0(code,
                 "\t\t\tY[j, i] ~ dbern(p[j, i])\n")

  if (nof1$model.linkfunc == "logit") {
    if (nof1$model.intcpt == "fixed") {
      code <- paste0(code,
                     "\t\t\tlogit(p[j, i]) <- alpha[i]")
    } else {
      code <- paste0(code,
                     "\t\t\tlogit(p[j, i]) <- beta_", nof1$Treat.name[1], "[i] * Treat_", nof1$Treat.name[1], "[j, i]")
    }

  } else if (nof1$model.linkfunc == "log") {
    code <- paste0(code,
                   "\t\t\tp[j, i] <- exp(rlp[j, i])\n")
    code <- paste0(code,
                   "\t\t\trlp[j, i] <- min(0, lp[j, i])\n")

    if (nof1$model.intcpt == "fixed") {
      code <- paste0(code,
                     "\t\t\tlp[j, i] <- alpha[i]")
    } else {
      code <- paste0(code,
                     "\t\t\tlp[j, i] <- beta_", nof1$Treat.name[1], "[i] * Treat_", nof1$Treat.name[1], "[j, i]")
    }
  }

  if (nof1$n.Treat > 1) {
    for(Treat.name.i in 2:nof1$n.Treat){
      code <- paste0(code,
                     " + beta_", nof1$Treat.name[Treat.name.i], "[i] * Treat_", nof1$Treat.name[Treat.name.i], "[j, i]")
    }
  }

  # adjust for covariates
  if (!is.null(nof1$names.covariates)) {

    for (covariates.i in 1:length(nof1$names.covariates)) {
      code <- paste0(code,
                     " + beta_", nof1$names.covariates[covariates.i], " * ", nof1$names.covariates[covariates.i], "[j, i]")
    }

  }

  # adjust for trend by basis splines
  if (nof1$spline.trend) {
    for (l in 1:nof1$spline_df){
      code <- paste0(code, " + eta[", l, "]*spline", l, "[time[j, i]]")
    }
  }
  code <- paste0(code, "\n")

  code <- paste0(code, "\t\t}\n") # for (j in 1:nobs.ID[i])
  code <- paste0(code, "\t}\n\n")  # for (i in 1:n.ID)

  # Priors
  if (nof1$model.intcpt == "fixed") {
    # prior for treatment effect
    code <- paste0(code,
                   "\tprec_delta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n")
    code <- paste0(code,
                   "\tsigmaSq_d <- pow(prec_delta, -1)\n\n")

    for (Treat.name.i in 2:nof1$n.Treat) {
      code <- paste0(code,
                     "\td[", Treat.name.i-1, "] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
      code <- paste0(code,
                     "\tSigma_delta[", Treat.name.i-1, ", ", Treat.name.i-1, "] <- sigmaSq_d\n")

      if ((Treat.name.i-1) > 1) {
        for (Treat.name.j in 1:(Treat.name.i-1-1)) {
          code <- paste0(code,
                         "\tSigma_delta[", Treat.name.i-1, ", ", Treat.name.j, "] <- sigmaSq_d / 2\n")
        }
      }

      if (Treat.name.i <= (nof1$n.Treat - 1)) {
        for (Treat.name.j in Treat.name.i:(nof1$n.Treat-1)) {
          code <- paste0(code,
                         "\tSigma_delta[", Treat.name.i-1, ", ", Treat.name.j, "] <- sigmaSq_d / 2\n")
        }
      }
    }

  } else { # if (nof1$model.intcpt == "fixed") {
    code <- paste0(code,
                   "\tprec_beta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n")
    code <- paste0(code,
                   "\tsigmaSq_beta <- pow(prec_beta, -1)\n")
    code <- paste0(code,
                   "\trho ~ ", nof1$rho.prior[[1]], "(", nof1$rho.prior[[2]], ", ", nof1$rho.prior[[3]], ")\n\n")

    code <- paste0(code,
                   "\tfor (k in 1:n.Treat) {\n")
    code <- paste0(code,
                   "\t\tb[k] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n\n") # for (k in 1:n.Treat) {

    # row 1
    code <- paste0(code,
                   "\tSigma_beta[1, 1] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\tfor (k in 2:n.Treat) {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[1, k] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t}\n") # for (k in 2:n.Treat) {
    # row 2:(n.Treat-1)
    code <- paste0(code,
                   "\tfor (k in 2:(n.Treat - 1)) {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[k, k] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\t\tfor (l in 1:(k-1)) {\n")
    code <- paste0(code,
                   "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t\t}\n") # for (l in 1:(k-1)) {
    code <- paste0(code,
                   "\t\tfor (l in (k+1):n.Treat) {\n")
    code <- paste0(code,
                   "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t\t}\n") # for (l in (k+1):n.Treat) {
    code <- paste0(code,
                   "\t}\n") # for (k in 2:(n.Treat - 1)) {
    # row n.Treat
    code <- paste0(code,
                   "\tSigma_beta[n.Treat, n.Treat] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\tfor (k in 1:(n.Treat - 1)) {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[n.Treat, k] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t}\n") # for (k in 1:(n.Treat - 1)) {
  }

  code <- paste0(code, "\n")

  # prior for covariate coefficients
  if (!is.null(nof1$names.covariates)) {

    for (covariates.i in 1:length(nof1$names.covariates)) {
      code <- paste0(code,
                     "\tbeta_", nof1$names.covariates[covariates.i], " ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    }
    code <- paste0(code, "\n")
  }

  # prior for trend coefficients - splines
  if (nof1$spline.trend) {
    code <- paste0(code,
                   "\tfor (l in 1:", nof1$spline_df, ") {\n")
    code <- paste0(code, "\t\teta[l] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n")
  }

  code <- paste0(code, "}") # model {

  # cat(code)
  return(code)
}

nof1.ma.ordinal.rjags <- function(nof1) {

  code <- paste0("model{\n")
  code <- paste0(code,
                 "\tfor (i in 1:n.ID) {\n")
  code <- paste0(code,
                 "\t\tfor (j in 1:nobs.ID[i]) {\n")

  code <- paste0(code,
                 "\t\t\tfor (r in 1:(ord.ncat - 1)) {\n")
  # contrasts
  code <- paste0(code,
                 "\t\t\t\tp[i, j, r+1] <- p[i, j, r] * c[i, j, r]\n")
  code <- paste0(code,
                 "\t\t\t\tlog(c[i, j, r]) <- beta_", nof1$Treat.name[1], "[i, r] * Treat_", nof1$Treat.name[1], "[j, i]")
  for(n.Treat.i in 2:nof1$n.Treat){
    code <- paste0(code,
                   " + beta_", nof1$Treat.name[n.Treat.i], "[i, r] * Treat_", nof1$Treat.name[n.Treat.i], "[j, i]")
  }
  code <- paste0(code, "\n")
  code <- paste0(code,
                 "\t\t\t}\n") # for (r in 1:(ord.ncat - 1)) {

  # at least 3 categories in the ordinal outcome because we always retain the missing levels
  # so at least 2 contrasts
  code <- paste0(code,
                 "\t\t\tp[i, j, 1] <- 1 / (1 + c[i, j, 1] + c[i, j, 1] * c[i, j, 2]")
  # need to test if there are more than 5 levels
  if (nof1$ord.ncat >= 4) {
    for (i in 4:nof1$ord.ncat) {
      code <- paste0(code,
                     " + c[i, j, 1]")
      for (j in 2:(i - 1)) {
        code <- paste0(code,
                       " * c[i, j, ", j, "]")
      }
    }
  }
  code <- paste0(code, ")\n")

  code <- paste0(code,
                 "\t\t\tY[j, i] ~ dcat(p[i, j, ])\n")
  code <- paste0(code,
                 "\t\t}\n\n") # for (j in 1:nobs.ID[i]) {

  # random effects
  if (nof1$ord.parallel) {
    code <- paste0(code,
                   "\t\tbeta[i, 1, 1:n.Treat] ~ dmnorm.vcov(b[1, 1:n.Treat], Sigma_beta[1:n.Treat, 1:n.Treat])\n")
    for(n.Treat.i in 1:nof1$n.Treat){
      code <- paste0(code,
                     "\t\tbeta_", nof1$Treat.name[n.Treat.i], "[i, 1] <- beta[i, 1, ", n.Treat.i, "]\n")
    }
    code <- paste0(code,
                   "\t\tfor (r in 2:(ord.ncat - 1)) {\n")
    code <- paste0(code,
                   "\t\t\tbeta[i, r, 1] ~ dnorm(b[r, 1], prec_beta)\n")
    code <- paste0(code,
                   "\t\t\tbeta_", nof1$Treat.name[1], "[i, r] <- beta[i, r, 1]\n")
    for(n.Treat.i in 2:nof1$n.Treat){
      code <- paste0(code,
                     "\t\t\tbeta_", nof1$Treat.name[n.Treat.i], "[i, r] <- beta_", nof1$Treat.name[1],
                     "[i, r] + (beta_", nof1$Treat.name[n.Treat.i], "[i, 1] - beta_", nof1$Treat.name[1], "[i, 1])\n")
    }
    code <- paste0(code,
                   "\t\t}\n") # for (r in 2:(ord.ncat - 1)) {
  } else {
    code <- paste0(code,
                   "\t\tfor (r in 1:(ord.ncat - 1)) {\n")
    code <- paste0(code,
                   "\t\t\tbeta[i, r, 1:n.Treat] ~ dmnorm.vcov(b[r, 1:n.Treat], Sigma_beta[1:n.Treat, 1:n.Treat])\n")
    for(n.Treat.i in 1:nof1$n.Treat){
      code <- paste0(code,
                     "\t\t\tbeta_", nof1$Treat.name[n.Treat.i], "[i, r] <- beta[i, r, ", n.Treat.i, "]\n")
    }
    code <- paste0(code,
                   "\t\t}\n") # for (r in 1:(ord.ncat - 1)) {
  }
  code <- paste0(code,
                 "\t}\n\n") # for (i in 1:n.ID) {

  # priors
  # prior for b
  if (nof1$ord.parallel) {
    code <- paste0(code,
                   "\tfor (k in 1:n.Treat) {\n")
    code <- paste0(code,
                   "\t\tb[1, k] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n") # for (k in 1:n.Treat) {
    code <- paste0(code,
                   "\tfor (r in 2:(ord.ncat - 1)) {\n")
    code <- paste0(code,
                   "\t\tb[r, 1] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n\n") # for (r in 2:(ord.ncat - 1)) {
  } else {
    code <- paste0(code,
                   "\tfor (r in 1:(ord.ncat - 1)) {\n")
    code <- paste0(code,
                   "\t\tfor (k in 1:n.Treat) {\n")
    code <- paste0(code,
                   "\t\t\tb[r, k] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t\t}\n") # for (k in 1:n.Treat) {
    code <- paste0(code,
                   "\t}\n\n") # for (r in 1:(ord.ncat - 1)) {
  }

  # prior for sigma
  # row 1
  code <- paste0(code,
                 "\tSigma_beta[1, 1] <- sigmaSq_beta\n")
  code <- paste0(code,
                 "\tfor (k in 2:n.Treat) {\n")
  code <- paste0(code,
                 "\t\tSigma_beta[1, k] <- rho * sigmaSq_beta\n")
  code <- paste0(code,
                 "\t}\n") # for (k in 2:n.Treat) {
  # row 2:(n.Treat-1)
  code <- paste0(code,
                 "\tfor (k in 2:(n.Treat - 1)) {\n")
  code <- paste0(code,
                 "\t\tSigma_beta[k, k] <- sigmaSq_beta\n")
  code <- paste0(code,
                 "\t\tfor (l in 1:(k-1)) {\n")
  code <- paste0(code,
                 "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
  code <- paste0(code,
                 "\t\t}\n") # for (l in 1:(k-1)) {
  code <- paste0(code,
                 "\t\tfor (l in (k+1):n.Treat) {\n")
  code <- paste0(code,
                 "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
  code <- paste0(code,
                 "\t\t}\n") # for (l in (k+1):n.Treat) {
  code <- paste0(code,
                 "\t}\n") # for (k in 2:(n.Treat - 1)) {
  # row n.Treat
  code <- paste0(code,
                 "\tSigma_beta[n.Treat, n.Treat] <- sigmaSq_beta\n")
  code <- paste0(code,
                 "\tfor (k in 1:(n.Treat - 1)) {\n")
  code <- paste0(code,
                 "\t\tSigma_beta[n.Treat, k] <- rho * sigmaSq_beta\n")
  code <- paste0(code,
                 "\t}\n\n") # for (k in 1:(n.Treat - 1)) {

  code <- paste0(code,
                 "\tsigmaSq_beta <- pow(prec_beta, -1)\n")
  code <- paste0(code,
                 "\tprec_beta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n")

  code <- paste0(code,
                 "\trho ~ ", nof1$rho.prior[[1]], "(", nof1$rho.prior[[2]], ", ", nof1$rho.prior[[3]], ")\n")


  code <- paste0(code,
                 "}") # model {

  # cat(code)
  return(code)
}

#########################################
######### Network meta analysis #########
#########################################
nof1.nma.rjags <- function(nof1) {

  if (nof1$response == "normal") {
    code <- nof1.nma.normal.rjags(nof1)
  } else if (nof1$response == "poisson") {
    code <- nof1.nma.poisson.rjags(nof1)
  } else if (nof1$response == "binomial") {
    code <- nof1.nma.binomial.rjags(nof1)
  }
  # else if (nof1$response == "ordinal") {
  #   code <- nof1.ma.ordinal.rjags(nof1)
  # }

  return(code)
}

nof1.nma.normal.rjags <- function(nof1) {

  code <- paste0("model{\n")

  if (nof1$summ.nID.perTreat$n_ID_perNTreat[1] != 0) {
    code <- paste0(code,
                   "\tfor (i in 1:", nof1$summ.nID.perTreat$n_ID_perNTreat[1], ") {\n")

    if (nof1$model.intcpt == "random") {
      code <- paste0(code,
                     "\t\tbeta_1[i] ~ dnorm(b[uniq.Treat.matrix[1, i]], prec_beta)\n")

      # covariate coefficients will be 0 any way for only one treatment if fixed intercept
      if (!is.null(nof1$cov.matrix)) {

        for (cov.i in 1:nrow(nof1$cov.matrix)) {
          code <- paste0(code,
                         "\t\teta_cov_i[i, 1, ", cov.i, "] <- eta_cov[uniq.Treat.matrix[1, i], ", cov.i, "]\n")
        }

      } # assign the actual eta cov
    } # if random intercept

    # First level model
    code <- paste0(code,
                   "\t\tfor (j in 1:nobs.ID[i]) {\n")

    if (!nof1$corr.y) {
      code <- paste0(code,
                     "\t\t\tY.matrix[j, i] ~ dnorm(m[j, i], prec_resid)\n")
    }
    # when correlation, Y.matrix will be the same across different # of treatment
    # will be defined differently j = 1 and j > 1 because index cannot be 0

    code <- paste0(code,
                   "\t\t\tm[j, i] <- mu[j, i]")
    if (nof1$corr.y) {
      code <- paste0(code,
                     " + e[j, i]")
    }
    code <- paste0(code,
                   "\n")

    if (nof1$model.intcpt == "fixed") {
      code <- paste0(code,
                     "\t\t\tmu[j, i] <- alpha[i]")
    } else if (nof1$model.intcpt == "random") {
      code <- paste0(code,
                     "\t\t\tmu[j, i] <- beta_1[i] * Treat.1[j, i]")
    }

    # adjust for trend by basis splines
    if (nof1$spline.trend) {
      for (l in 1:nof1$spline.df){
        code <- paste0(code, " + eta[", l, "] * spline.matrix[time.matrix[j, i], ", l,"]")
      }
    }

    # adjust for trend by step functions
    if (nof1$step.trend) {
      for (l in 1:nof1$step.df){
        code <- paste0(code, " + eta[", l, "] * step.matrix[period.matrix[j, i], ", l,"]")
      }
    }

    # adjust for level 2 covariates
    # only random intercept because the coefficients for fixed intercept when one treatment being 0 anyway
    if (!is.null(nof1$cov.matrix)) {
      if (nof1$model.intcpt == "random") {

        for (cov.i in 1:nrow(nof1$cov.matrix)) {
          code <- paste0(code,
                         " + eta_cov_i[i, 1, ", cov.i, "] * cov.matrix[", cov.i, ", i] * Treat.1[j, i]")
        }

      } # if (nof1$model.intcpt == "random") {
    } # if (!is.null(nof1$cov.matrix)) {
    code <- paste0(code, "\n")

    code <- paste0(code,
                   "\t\t}\n") # for (j in 1:nobs.ID[i])

    code <- paste0(code,
                   "\t}\n\n") # for (i in cumsum:cumsum)
  } # if (nof1$summ.nID.perTreat$n_ID_perNTreat[1] != 0) {

  # when there is only 1 treatment on all participant, do not append code
  if (nrow(nof1$summ.nID.perTreat) >= 2) {

    # summ.nID.perTreat always have 1:maximum number of treatments per participant
    for (n.Treat.ID.i in 2:nrow(nof1$summ.nID.perTreat)) {

      # only append code when there are participants has this number of treatments
      if (nof1$summ.nID.perTreat$n_ID_perNTreat[n.Treat.ID.i] != 0) {

        code <- paste0(code,
                       "\tfor (i in ", cumsum(nof1$summ.nID.perTreat$n_ID_perNTreat)[n.Treat.ID.i-1]+1, ":", cumsum(nof1$summ.nID.perTreat$n_ID_perNTreat)[n.Treat.ID.i], ") {\n")

        # random effects
        if (nof1$model.intcpt == "fixed") {
          for (Treat.i in 2:n.Treat.ID.i) {
            code <- paste0(code,
                           "\t\tbeta_", Treat.i, "[i] ~ dnorm(d[uniq.Treat.matrix[", Treat.i, ", i]] - d[uniq.Treat.matrix[1, i]], prec_beta)\n")
          }
        } else if (nof1$model.intcpt == "random") {
          code <- paste0(code,
                         "\t\tbeta[i, 1:", n.Treat.ID.i,
                         "] ~ dmnorm.vcov(b[uniq.Treat.matrix[1:", n.Treat.ID.i, ", i]], Sigma_beta[uniq.Treat.matrix[1:", n.Treat.ID.i, ", i], uniq.Treat.matrix[1:", n.Treat.ID.i, ", i]])\n")
          for (Treat.i in 1:n.Treat.ID.i) {
            code <- paste0(code,
                           "\t\tbeta_", Treat.i, "[i] <- beta[i, ", Treat.i, "]\n")
          }
        }

        # assign the actual eta cov
        if (!is.null(nof1$cov.matrix)) {

          if (nof1$model.intcpt == "fixed") {
            for (Treat.i in 2:n.Treat.ID.i) {
              for (cov.i in 1:nrow(nof1$cov.matrix)) {
                code <- paste0(code,
                               "\t\teta_cov_i[i, ", Treat.i, ", ", cov.i, "] <- eta_cov[uniq.Treat.matrix[", Treat.i, ", i], ", cov.i, "] - eta_cov[uniq.Treat.matrix[1, i], ", cov.i, "]\n")
              }
            }
          } else if (nof1$model.intcpt == "random") {
            for (Treat.i in 1:n.Treat.ID.i) {
              for (cov.i in 1:nrow(nof1$cov.matrix)) {
                code <- paste0(code,
                               "\t\teta_cov_i[i, ", Treat.i, ", ", cov.i, "] <- eta_cov[uniq.Treat.matrix[", Treat.i, ", i], ", cov.i, "]\n")
              }
            }
          } # else if (nof1$model.intcpt == "random") {

        } # assign the actual eta cov

        # First level model
        code <- paste0(code,
                       "\t\tfor (j in 1:nobs.ID[i]) {\n")
        if (!nof1$corr.y) {
          code <- paste0(code,
                         "\t\t\tY.matrix[j, i] ~ dnorm(m[j, i], prec_resid)\n")
        }
        # when correlation, Y.matrix will be the same across different # of treatment
        # will be defined differently j = 1 and j > 1 because index cannot be 0

        code <- paste0(code,
                       "\t\t\tm[j, i] <- mu[j, i]")
        # correlation
        if (nof1$corr.y) {
          code <- paste0(code,
                         " + e[j, i]")
        }
        code <- paste0(code,
                       "\n")

        if (nof1$model.intcpt == "fixed") {
          code <- paste0(code,
                         "\t\t\tmu[j, i] <- alpha[i]")
        } else if (nof1$model.intcpt == "random") {
          code <- paste0(code,
                         "\t\t\tmu[j, i] <- beta_1[i] * Treat.1[j, i]")
        }
        for (Treat.i in 2:n.Treat.ID.i){
          code <- paste0(code,
                         " + beta_", Treat.i, "[i] * Treat.", Treat.i, "[j, i]")
        }

        # adjust for trend by basis splines
        if (nof1$spline.trend) {
          for (l in 1:nof1$spline.df){
            code <- paste0(code, " + eta[", l, "] * spline.matrix[time.matrix[j, i], ", l,"]")
          }
        }

        # adjust for trend by step functions
        if (nof1$step.trend) {
          for (l in 1:nof1$step.df){
            code <- paste0(code, " + eta[", l, "] * step.matrix[period.matrix[j, i], ", l,"]")
          }
        }

        # adjust for level 2 covariates
        if (!is.null(nof1$cov.matrix)) {

          # random intercept also first intercept
          if (nof1$model.intcpt == "random") {
            for (cov.i in 1:nrow(nof1$cov.matrix)) {
              code <- paste0(code,
                             " + eta_cov_i[i, 1, ", cov.i, "] * cov.matrix[", cov.i, ", i] * Treat.1[j, i]")
            }
          }

          # fixed intercept do not need to do interaction with first treatment
          for (Treat.i in 2:n.Treat.ID.i) {
            for (cov.i in 1:nrow(nof1$cov.matrix)) {
              code <- paste0(code,
                             " + eta_cov_i[i, ", Treat.i, ", ", cov.i, "] * cov.matrix[", cov.i, ", i] * Treat.", Treat.i, "[j, i]")
            }
          }
          # eta_cov be a matrix of coefficients for the interaction between treatment and covariates
        } # if (!is.null(nof1$cov.matrix)) {
        code <- paste0(code, "\n")

        code <- paste0(code,
                       "\t\t}\n") # for (j in 1:nobs.ID[i])
        code <- paste0(code,
                       "\t}\n\n") # for (i in cumsum:cumsum)

      } # if (nof1$summ.nID.perTreat$n_ID_perNTreat[n.Treat.ID.i] != 0) {
    } # for (n.Treat.ID.i in 2:nrow(nof1$summ.nID.perTreat)) {
  }

  # e[j, i], resid[j, i] for correlation
  # take care of outcome measurements on non-consecutive days
  if (nof1$corr.y) {
    code <- paste0(code,
                   "\tfor (i in 1:", sum(nof1$summ.nID.perTreat$n_ID_perNTreat),"){\n")
    code <- paste0(code,
                   "\t\te[1, i] <- 0\n")
    code <- paste0(code,
                   "\t\tY.matrix[1, i] ~ dnorm(m[1, i], (1 - rho_resid^2) * prec_resid)\n")
    code <- paste0(code,
                   "\t\tfor (j in 2:nobs.ID[i]) {\n")
    code <- paste0(code,
                   "\t\t\te[j, i] <- rho_resid^(time.matrix[j, i] - time.matrix[j-1, i]) * (Y.matrix[j-1, i] - mu[j-1, i])\n")
    code <- paste0(code,
                   "\t\t\tY.matrix[j, i] ~ dnorm(m[j, i], prec_resid * (1 - rho_resid^2) / (1 - (rho_resid^2)^(time.matrix[j, i] - time.matrix[j-1, i])))\n")
    code <- paste0(code,
                   "\t\t}\n") # j loop
    code <- paste0(code,
                   "\t}\n\n") # i loop
  }

  # Priors
  if (nof1$model.intcpt == "fixed") {
    # alpha prior
    code <- paste0(code,
                   "\tfor (i in 1:", sum(nof1$summ.nID.perTreat$n_ID_perNTreat), ") {\n")
    code <- paste0(code,
                   "\t\talpha[i] ~ ", nof1$alpha.prior[[1]], "(", nof1$alpha.prior[[2]], ", ", nof1$alpha.prior[[3]], ")")
    # Truncated normal distribution for beta's
    if (length(nof1$alpha.prior) > 3) {

      if (!is.na(nof1$alpha.prior[[4]])) {
        code <- paste0(code,
                       " T(", nof1$alpha.prior[[4]], ", ")
      } else {
        code <- paste0(code,
                       " T(, ")
      }

      if(!is.na(nof1$alpha.prior[[5]])) {
        code <- paste0(code,
                       nof1$alpha.prior[[5]], ")")
      } else {
        code <- paste0(code,
                       ")")
      }
    } # if (length(nof1$alpha.prior) > 3) {
    code <- paste0(code, "\n")
    code <- paste0(code, "\t}\n\n") # for (i in 1:", n.ID, ") {"

    # prior for d
    code <- paste0(code,
                   "\td[1] <- 0\n")
    code <- paste0(code,
                   "\tfor (i in 2:", length(nof1$Treat.name), ") {\n")
    code <- paste0(code,
                   "\t\td[i] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    code <- paste0(code, "\t}\n\n") # for (i in 1:", n.Treat, ") {"

  } else if (nof1$model.intcpt == "random") {
    # prior for b
    code <- paste0(code,
                   "\tfor (i in 1:", length(nof1$Treat.name), ") {\n")
    code <- paste0(code,
                   "\t\tb[i] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
    code <- paste0(code, "\t}\n\n") # for (i in 1:", n.Treat, ") {"

    # prior for Sigma_beta_...
    code <- paste0(code,
                   "\tsigmaSq_beta <- pow(prec_beta, -1)\n")
    code <- paste0(code,
                   "\trho ~ ", nof1$rho.prior[[1]], "(", nof1$rho.prior[[2]], ", ", nof1$rho.prior[[3]], ")\n")

    # row 1
    code <- paste0(code,
                   "\tSigma_beta[1, 1] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\tfor (k in 2:", length(nof1$Treat.name), ") {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[1, k] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t}\n") # for (k in 2:n.Treat) {

    # row 2:(n.Treat-1), only if number of treatments greater than 2
    if (length(nof1$Treat.name) > 2) {

      code <- paste0(code,
                     "\tfor (k in 2:", length(nof1$Treat.name) - 1, ") {\n")
      code <- paste0(code,
                     "\t\tSigma_beta[k, k] <- sigmaSq_beta\n")
      code <- paste0(code,
                     "\t\tfor (l in 1:(k-1)) {\n")
      code <- paste0(code,
                     "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
      code <- paste0(code,
                     "\t\t}\n") # for (l in 1:(k-1)) {
      code <- paste0(code,
                     "\t\tfor (l in (k+1):", length(nof1$Treat.name), ") {\n")
      code <- paste0(code,
                     "\t\t\tSigma_beta[k, l] <- rho * sigmaSq_beta\n")
      code <- paste0(code,
                     "\t\t}\n") # for (l in (k+1):n.Treat) {
      code <- paste0(code,
                     "\t}\n") # for (k in 2:(n.Treat - 1)) {
    }
    # row n.Treat
    code <- paste0(code,
                   "\tSigma_beta[", length(nof1$Treat.name), ", ", length(nof1$Treat.name), "] <- sigmaSq_beta\n")
    code <- paste0(code,
                   "\tfor (k in 1:", length(nof1$Treat.name) - 1, ") {\n")
    code <- paste0(code,
                   "\t\tSigma_beta[", length(nof1$Treat.name), ", k] <- rho * sigmaSq_beta\n")
    code <- paste0(code,
                   "\t}\n") # for (k in 1:(n.Treat - 1)) {

  } # else if (nof1$model.intcpt == "random") {

  # prior for prec_beta
  code <- paste0(code,
                 "\tprec_beta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n\n")

  # prior for residual error
  code <- paste0(code,
                 "\tprec_resid ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n")
  # prior for serial correlation
  if (nof1$corr.y) {
    code <- paste0(code,
                   "\trho_resid ~ ", nof1$rho.prior[[1]], "(", nof1$rho.prior[[2]], ", ", nof1$rho.prior[[3]], ")\n")
  }
  code <- paste0(code,
                 "\n")

  # prior for trend coefficients - splines
  if (nof1$spline.trend) {
    code <- paste0(code,
                   "\tfor (l in 1:", nof1$spline.df, ") {\n")
    code <- paste0(code, "\t\teta[l] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n\n")
  }

  # prior for trend coefficients - steps
  if (nof1$step.trend) {
    code <- paste0(code,
                   "\tfor (l in 1:", nof1$step.df, ") {\n")
    code <- paste0(code, "\t\teta[l] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n\n")
  }

  # prior for coefficients for level 2 covariates
  if (!is.null(nof1$cov.matrix)) {

    code <- paste0(code,
                   "\tfor (j in 1:", nrow(nof1$cov.matrix), ") {\n")

    if (nof1$model.intcpt == "fixed") {
      code <- paste0(code,
                     "\t\teta_cov[1, j] <- 0\n")
      code <- paste0(code,
                     "\t\tfor (i in 2:", length(nof1$Treat.name), ") {\n")
    } else if (nof1$model.intcpt == "random") {
      code <- paste0(code,
                     "\t\tfor (i in 1:", length(nof1$Treat.name), ") {\n")
    }

    code <- paste0(code, "\t\t\teta_cov[i, j] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t\t}\n")
    code <- paste0(code,
                   "\t}\n\n")
  }

  code <- paste0(code,
                 "}") # model {

  # cat(code)
  return(code)
}


nof1.nma.poisson.rjags <- function(nof1) {

  code <- paste0("model{\n")

  if (nof1$summ.nID.perTreat$n_ID_perNTreat[1] != 0) {
    code <- paste0(code,
                   "\tfor (i in 1:", nof1$summ.nID.perTreat$n_ID_perNTreat[1], ") {\n")
    # First level model
    code <- paste0(code,
                   "\t\tfor (j in 1:nobs.ID[i]) {\n")
    code <- paste0(code,
                   "\t\t\tY.matrix[j, i] ~ dpois(lambda[j, i])\n")
    if (nof1$model.intcpt == "fixed") {
      code <- paste0(code,
                     "\t\t\tlog(lambda[j, i]) <- alpha[i]")
    }

    # adjust for trend by basis splines
    if (nof1$spline.trend) {
      for (l in 1:nof1$spline.df){
        code <- paste0(code, " + eta[", l, "] * spline.matrix[time.matrix[j, i], ", l,"]")
      }
    }
    code <- paste0(code, "\n")

    code <- paste0(code,
                   "\t\t}\n") # for (j in 1:nobs.ID[i])

    code <- paste0(code,
                   "\t}\n\n") # for (i in cumsum:cumsum)
  } # if (nof1$summ.nID.perTreat$n_ID_perNTreat[1] != 0) {

  # when there is only 1 treatment on all participant, do not append code
  if (nrow(nof1$summ.nID.perTreat) >= 2) {

    for (n.Treat.ID.i in 2:nrow(nof1$summ.nID.perTreat)) {

      # only append code when there are participants has this number of treatments
      if (nof1$summ.nID.perTreat$n_ID_perNTreat[n.Treat.ID.i] != 0) {

        code <- paste0(code,
                       "\tfor (i in ", cumsum(nof1$summ.nID.perTreat$n_ID_perNTreat)[n.Treat.ID.i-1]+1, ":", cumsum(nof1$summ.nID.perTreat$n_ID_perNTreat)[n.Treat.ID.i], ") {\n")

        # random effects
        for (Treat.i in 2:n.Treat.ID.i) {
          code <- paste0(code,
                         "\t\tbeta_", Treat.i, "[i] ~ dnorm(d[uniq.Treat.matrix[", Treat.i, ", i]] - d[uniq.Treat.matrix[1, i]], prec_beta)\n")
        }

        # assign the actual eta cov
        if (!is.null(nof1$cov.matrix)) {

          if (nof1$model.intcpt == "fixed") {
            for (Treat.i in 2:n.Treat.ID.i) {
              for (cov.i in 1:nrow(nof1$cov.matrix)) {
                code <- paste0(code,
                               "\t\teta_cov_i[i, ", Treat.i, ", ", cov.i, "] <- eta_cov[uniq.Treat.matrix[", Treat.i, ", i], ", cov.i, "] - eta_cov[uniq.Treat.matrix[1, i], ", cov.i, "]\n")
              }
            }
          }
          # else if (nof1$model.intcpt == "random") {
          #   for (Treat.i in 1:n.Treat.ID.i) {
          #     for (cov.i in 1:nrow(nof1$cov.matrix)) {
          #       code <- paste0(code,
          #                      "\t\teta_cov_i[i, ", Treat.i, ", ", cov.i, "] <- eta_cov[uniq.Treat.matrix[", Treat.i, ", i], ", cov.i, "]\n")
          #     }
          #   }
          # } # else if (nof1$model.intcpt == "random") {

        } # if (!is.null(nof1$cov.matrix)) {

        # First level model
        code <- paste0(code,
                       "\t\tfor (j in 1:nobs.ID[i]) {\n")
        code <- paste0(code,
                       "\t\t\tY.matrix[j, i] ~ dpois(lambda[j, i])\n")

        if (nof1$model.intcpt == "fixed") {
          code <- paste0(code,
                         "\t\t\tlog(lambda[j, i]) <- alpha[i]")
        }
        for(Treat.i in 2:n.Treat.ID.i){
          code <- paste0(code,
                         " + beta_", Treat.i, "[i] * Treat.", Treat.i, "[j, i]")
        }

        # adjust for trend by basis splines
        if (nof1$spline.trend) {
          for (l in 1:nof1$spline.df){
            code <- paste0(code, " + eta[", l, "] * spline.matrix[time.matrix[j, i], ", l,"]")
          }
        }

        # adjust for level 2 covariates
        if (!is.null(nof1$cov.matrix)) {

          # # random intercept also first intercept
          # if (nof1$model.intcpt == "random") {
          #   for (cov.i in 1:nrow(nof1$cov.matrix)) {
          #     code <- paste0(code,
          #                    " + eta_cov_i[i, 1, ", cov.i, "] * cov.matrix[", cov.i, ", i] * Treat.1[j, i]")
          #   }
          # }

          # work for fixed intercept model only
          for (Treat.i in 2:n.Treat.ID.i) {
            for (cov.i in 1:nrow(nof1$cov.matrix)) {
              code <- paste0(code,
                             " + eta_cov_i[i, ", Treat.i, ", ", cov.i, "] * cov.matrix[", cov.i, ", i] * Treat.", Treat.i, "[j, i]")
            }
          }
          # eta_cov be a matrix of coefficients for the interaction between treatment and covariates
        } # if (!is.null(nof1$cov.matrix)) {
        code <- paste0(code, "\n")

        code <- paste0(code,
                       "\t\t}\n") # for (j in 1:nobs.ID[i])
        code <- paste0(code,
                       "\t}\n\n") # for (i in cumsum:cumsum)

      } # if (nof1$summ.nID.perTreat$n_ID_perNTreat[n.Treat.ID.i] != 0) {
    } # for (n.Treat.ID.i in 2:nrow(nof1$summ.nID.perTreat)) {
  }

  # Priors
  # alpha prior
  code <- paste0(code,
                 "\tfor (i in 1:", sum(nof1$summ.nID.perTreat$n_ID_perNTreat), ") {\n")
  code <- paste0(code,
                 "\t\talpha[i] ~ ", nof1$alpha.prior[[1]], "(", nof1$alpha.prior[[2]], ", ", nof1$alpha.prior[[3]], ")")
  # Truncated normal distribution for beta's
  if (length(nof1$alpha.prior) > 3) {

    if (!is.na(nof1$alpha.prior[[4]])) {
      code <- paste0(code,
                     " T(", nof1$alpha.prior[[4]], ", ")
    } else {
      code <- paste0(code,
                     " T(, ")
    }

    if(!is.na(nof1$alpha.prior[[5]])) {
      code <- paste0(code,
                     nof1$alpha.prior[[5]], ")")
    } else {
      code <- paste0(code,
                     ")")
    }
  } # if (length(nof1$alpha.prior) > 3) {
  code <- paste0(code, "\n")
  code <- paste0(code, "\t}\n\n") # for (i in 1:", n.ID, ") {"

  # prior for d
  code <- paste0(code,
                 "\td[1] <- 0\n")
  code <- paste0(code,
                 "\tfor (i in 2:", length(nof1$Treat.name), ") {\n")
  code <- paste0(code,
                 "\t\td[i] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
  code <- paste0(code, "\t}\n\n") # for (i in 1:", n.Treat, ") {"

  # prior for prec_beta
  code <- paste0(code,
                 "\tprec_beta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n\n")

  # prior for trend coefficients - splines
  if (nof1$spline.trend) {
    code <- paste0(code,
                   "\tfor (l in 1:", nof1$spline.df, ") {\n")
    code <- paste0(code, "\t\teta[l] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n\n")
  }

  # prior for coefficients for level 2 covariates
  if (!is.null(nof1$cov.matrix)) {

    code <- paste0(code,
                   "\tfor (j in 1:", nrow(nof1$cov.matrix), ") {\n")

    if (nof1$model.intcpt == "fixed") {
      code <- paste0(code,
                     "\t\teta_cov[1, j] <- 0\n")
      code <- paste0(code,
                     "\t\tfor (i in 2:", length(nof1$Treat.name), ") {\n")
    }
    # else if (nof1$model.intcpt == "random") {
    #   code <- paste0(code,
    #                  "\t\tfor (i in 1:", length(nof1$Treat.name), ") {\n")
    # }
    code <- paste0(code, "\t\t\teta_cov[i, j] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t\t}\n")
    code <- paste0(code,
                   "\t}\n\n")
  }

  code <- paste0(code,
                 "}") # model {

  # cat(code)
  return(code)
}

nof1.nma.binomial.rjags <- function(nof1) {

  # only log link, fixed intercept work for now

  code <- paste0("model{\n")

  if (nof1$summ.nID.perTreat$n_ID_perNTreat[1] != 0) {
    code <- paste0(code,
                   "\tfor (i in 1:", nof1$summ.nID.perTreat$n_ID_perNTreat[1], ") {\n")
    # First level model
    code <- paste0(code,
                   "\t\tfor (j in 1:nobs.ID[i]) {\n")

    code <- paste0(code,
                   "\t\t\tY.matrix[j, i] ~ dbern(p[j, i])\n")

    if (nof1$model.linkfunc == "log") {
      code <- paste0(code,
                     "\t\t\tp[j, i] <- exp(rlp[j, i])\n")
      code <- paste0(code,
                     "\t\t\trlp[j, i] <- min(0, lp[j, i])\n")

      if (nof1$model.intcpt == "fixed") {
        code <- paste0(code,
                       "\t\t\tlp[j, i] <- alpha[i]")
      }

    } else if (nof1$model.linkfunc == "logit") {

      if (nof1$model.intcpt == "fixed") {
        code <- paste0(code,
                       "\t\t\tlogit(p[j, i]) <- alpha[i]")
      }

    } # else if logit link

    # adjust for trend by basis splines
    if (nof1$spline.trend) {
      for (l in 1:nof1$spline.df){
        code <- paste0(code, " + eta[", l, "] * spline.matrix[time.matrix[j, i], ", l,"]")
      }
    }
    code <- paste0(code, "\n")

    code <- paste0(code,
                   "\t\t}\n") # for (j in 1:nobs.ID[i])

    code <- paste0(code,
                   "\t}\n\n") # for (i in cumsum:cumsum)

  }

  # when there is only 1 treatment on all participant, do not append code
  if (nrow(nof1$summ.nID.perTreat) >= 2) {

    for (n.Treat.ID.i in 2:nrow(nof1$summ.nID.perTreat)) {

      # only append code when there are participants has this number of treatments
      if (nof1$summ.nID.perTreat$n_ID_perNTreat[n.Treat.ID.i] != 0) {

        code <- paste0(code,
                       "\tfor (i in ", cumsum(nof1$summ.nID.perTreat$n_ID_perNTreat)[n.Treat.ID.i-1]+1, ":", cumsum(nof1$summ.nID.perTreat$n_ID_perNTreat)[n.Treat.ID.i], ") {\n")

        # random effects
        for (Treat.i in 2:n.Treat.ID.i) {
          ## need to fix here
          code <- paste0(code,
                         "\t\tbeta_", Treat.i, "[i] ~ dnorm(d[uniq.Treat.matrix[", Treat.i, ", i]] - d[uniq.Treat.matrix[1, i]], prec_beta)\n")
        }

        # assign the actual eta cov
        if (!is.null(nof1$cov.matrix)) {

          if (nof1$model.intcpt == "fixed") {
            for (Treat.i in 2:n.Treat.ID.i) {
              for (cov.i in 1:nrow(nof1$cov.matrix)) {
                code <- paste0(code,
                               "\t\teta_cov_i[i, ", Treat.i, ", ", cov.i, "] <- eta_cov[uniq.Treat.matrix[", Treat.i, ", i], ", cov.i, "] - eta_cov[uniq.Treat.matrix[1, i], ", cov.i, "]\n")
              }
            }
          }
          # else if (nof1$model.intcpt == "random") {
          #   for (Treat.i in 1:n.Treat.ID.i) {
          #     for (cov.i in 1:nrow(nof1$cov.matrix)) {
          #       code <- paste0(code,
          #                      "\t\teta_cov_i[i, ", Treat.i, ", ", cov.i, "] <- eta_cov[uniq.Treat.matrix[", Treat.i, ", i], ", cov.i, "]\n")
          #     }
          #   }
          # } # else if (nof1$model.intcpt == "random") {

        } # if (!is.null(nof1$cov.matrix)) {

        # First level model
        code <- paste0(code,
                       "\t\tfor (j in 1:nobs.ID[i]) {\n")

        code <- paste0(code,
                       "\t\t\tY.matrix[j, i] ~ dbern(p[j, i])\n")

        if (nof1$model.linkfunc == "log") {
          code <- paste0(code,
                         "\t\t\tp[j, i] <- exp(rlp[j, i])\n")
          code <- paste0(code,
                         "\t\t\trlp[j, i] <- min(0, lp[j, i])\n")

          if (nof1$model.intcpt == "fixed") {
            code <- paste0(code,
                           "\t\t\tlp[j, i] <- alpha[i]")
          }

        } else if (nof1$model.linkfunc == "logit") {

          if (nof1$model.intcpt == "fixed") {
            code <- paste0(code,
                           "\t\t\tlogit(p[j, i]) <- alpha[i]")
          }

        } # else if logit link

        for(Treat.i in 2:n.Treat.ID.i){
          code <- paste0(code,
                         " + beta_", Treat.i, "[i] * Treat.", Treat.i, "[j, i]")
        }

        # adjust for trend by basis splines
        if (nof1$spline.trend) {
          for (l in 1:nof1$spline.df){
            code <- paste0(code, " + eta[", l, "] * spline.matrix[time.matrix[j, i], ", l,"]")
          }
        }

        # adjust for level 2 covariates
        if (!is.null(nof1$cov.matrix)) {

          # # random intercept also first intercept
          # if (nof1$model.intcpt == "random") {
          #   for (cov.i in 1:nrow(nof1$cov.matrix)) {
          #     code <- paste0(code,
          #                    " + eta_cov_i[i, 1, ", cov.i, "] * cov.matrix[", cov.i, ", i] * Treat.1[j, i]")
          #   }
          # }

          # work for fixed intercept model only
          for (Treat.i in 2:n.Treat.ID.i) {
            for (cov.i in 1:nrow(nof1$cov.matrix)) {
              code <- paste0(code,
                             " + eta_cov_i[i, ", Treat.i, ", ", cov.i, "] * cov.matrix[", cov.i, ", i] * Treat.", Treat.i, "[j, i]")
            }
          }
          # eta_cov be a matrix of coefficients for the interaction between treatment and covariates
        } # if (!is.null(nof1$cov.matrix)) {
        code <- paste0(code, "\n")

        code <- paste0(code,
                       "\t\t}\n") # for (j in 1:nobs.ID[i])
        code <- paste0(code,
                       "\t}\n\n") # for (i in cumsum:cumsum)

      } # if (nof1$summ.nID.perTreat$n_ID_perNTreat[n.Treat.ID.i] != 0) {
    } # for (n.Treat.ID.i in 2:nrow(nof1$summ.nID.perTreat)) {
  } # if (nrow(nof1$summ.nID.perTreat) >= 2) {

  # Priors
  # alpha prior
  code <- paste0(code,
                 "\tfor (i in 1:", sum(nof1$summ.nID.perTreat$n_ID_perNTreat), ") {\n")
  code <- paste0(code,
                 "\t\talpha[i] ~ ", nof1$alpha.prior[[1]], "(", nof1$alpha.prior[[2]], ", ", nof1$alpha.prior[[3]], ")")
  # Truncated normal distribution for beta's
  if (length(nof1$alpha.prior) > 3) {

    if (!is.na(nof1$alpha.prior[[4]])) {
      code <- paste0(code,
                     " T(", nof1$alpha.prior[[4]], ", ")
    } else {
      code <- paste0(code,
                     " T(, ")
    }

    if(!is.na(nof1$alpha.prior[[5]])) {
      code <- paste0(code,
                     nof1$alpha.prior[[5]], ")")
    } else {
      code <- paste0(code,
                     ")")
    }
  } # if (length(nof1$alpha.prior) > 3) {
  code <- paste0(code, "\n")
  code <- paste0(code, "\t}\n\n") # for (i in 1:", n.ID, ") {"

  # prior for d
  code <- paste0(code,
                 "\td[1] <- 0\n")
  code <- paste0(code,
                 "\tfor (i in 2:", length(nof1$Treat.name), ") {\n")
  code <- paste0(code,
                 "\t\td[i] ~ ", nof1$beta.prior[[1]], "(", nof1$beta.prior[[2]], ", ", nof1$beta.prior[[3]], ")\n")
  code <- paste0(code, "\t}\n\n") # for (i in 1:", n.Treat, ") {"

  # prior for prec_beta
  code <- paste0(code,
                 "\tprec_beta <- pow(sd_beta, -2)\n")
  code <- paste0(code,
                 "\tsd_beta ~ ", nof1$hy.prior[[1]], "(", nof1$hy.prior[[2]], ", ", nof1$hy.prior[[3]], ")\n\n")

  # prior for trend coefficients - splines
  if (nof1$spline.trend) {
    code <- paste0(code,
                   "\tfor (l in 1:", nof1$spline.df, ") {\n")
    code <- paste0(code, "\t\teta[l] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t}\n\n")
  }

  # prior for coefficients for level 2 covariates
  if (!is.null(nof1$cov.matrix)) {

    code <- paste0(code,
                   "\tfor (j in 1:", nrow(nof1$cov.matrix), ") {\n")

    if (nof1$model.intcpt == "fixed") {
      code <- paste0(code,
                     "\t\teta_cov[1, j] <- 0\n")
      code <- paste0(code,
                     "\t\tfor (i in 2:", length(nof1$Treat.name), ") {\n")
    }
    # else if (nof1$model.intcpt == "random") {
    #   code <- paste0(code,
    #                  "\t\tfor (i in 1:", length(nof1$Treat.name), ") {\n")
    # }

    code <- paste0(code, "\t\t\teta_cov[i, j] ~ ", nof1$eta.prior[[1]], "(", nof1$eta.prior[[2]], ", ", nof1$eta.prior[[3]], ")\n")
    code <- paste0(code,
                   "\t\t}\n")
    code <- paste0(code,
                   "\t}\n\n")
  }

  code <- paste0(code,
                 "}") # model {

  # cat(code)
  return(code)
}
