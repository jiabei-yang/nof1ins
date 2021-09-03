#' Make an N of 1 object containing data, priors, and a jags model file for individual analysis
#'
#' @param Y Outcome of the study. This should be a vector with \code{NA}'s included in time order.
#' @param Treat Treatment indicator vector with same length as the outcome. Can be character or numeric.
#' @param ord.baseline.Treat Used for ordinal outcome to define a reference treatment level.
#' @param ord.model Used for ordinal outcome to pick the model. Can be "cumulative" for proportional odds model
#' or "acat" for restricted adjacent category model.
#' @param response Type of outcome. Can be "normal" for continuous outcome, "binomial" for binary outcome,
#' "poisson" for count outcome, or "ordinal" for ordinal or nominal outcome.
#' @param ncat Number of categories. Used in ordinal models.
#' @param bs.trend Indicator for whether the model should adjust for trend using splines. The default
#' is \code{F}.
#' @param y.time Parameter used for modeling splines. Time when the outcome is measured.
#' @param knots.bt.block parameter used for modeling splines. Indicator for whether or not knots should be set
#' at the end of each block except for the last block. If \code{TRUE}, user should then specify \code{block.no};
#' if \code{FALSE}, user should then specify \code{bs.df}.
#' @param block.no Block number used for modeling splines for the setting where the knots are set at the end
#' of each block. Block number with the same length as the outcome.
#' @param bs.df Degrees of freedom for modeling splines when knots are not set at the end of each block.
#' @param corr.y Indicator for whether the correlation among measurements shoule be modeled. The default is
#' \code{F}.
#' @param alpha.prior Prior for the intercept of the model. Not needed now since we are using treatment-specific
#' intercept model.
#' @param beta.prior Prior for the treatment-specific intercept.
#' @param eta.prior Prior for modelling spline terms.
#' @param dc.prior Prior for the length between cutpoints. Used only for ordinal logistic models.
#' @param c1.prior Prior for the first cutpoint. Used only for ordinal logistic models.
#' @param rho.prior Prior for the correlated error model.
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for
#' normal and binomial response and wishart for multinomial response. It should be a list of length 3,
#' where first element should be the distribution (one of dunif, dgamma, dhnorm, dwish) and the next
#' two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give
#' uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter. For wishart
#' distribution, the last two parameter would be the scale matrix and the degrees of freedom.
#' @param ... Arguments to be passed to \code{bs()}.
#' @return An object of class "nof1.data" that is used to run the model using \code{\link{nof1.run}} is a list
#' containing
#' \item{Y}{Outcome}
#' \item{Treat}{Treatment}
#' \item{ncat}{Number of categories for ordinal response}
#' \item{nobs}{Total number of observations in a study}
#' \item{Treat.name}{Treatment name}
#' \item{response}{Type of outcome}
#' \item{Treat_\emph{Treat.name}}{Vector in the model matrix for \emph{Treat.name}}
#' \item{bs.trend}{Indicator for whether the model should adjust for trend using splines}
#' \item{corr.y}{Indicator for whether the correlation among measurements shoule be modeled}
#' \item{priors}{Priors that the code will be using. Default priors are used if prior was not specified}
#' \item{code}{Rjags model file code that is generated using information provided by the user. To view model file inside R, use \code{cat(nof1$code).}}
#' @examples
#' ###Blocker data example
#' laughter
#' Y <- laughter$Y
#' Treat <- laughter$Treat
#' nof1 <- nof1.data(Y, Treat, ncat = 11, response = "ordinal")
#' str(nof1)
#' cat(nof1$code)
#' @export

# @param baseline baseline Treatment name. This serves as a baseline/placebo when comparing different treatments.
# \item{baseline}{Baseline variable}
nof1.data <- function(Y, Treat, ord.baseline.Treat = NULL, ord.model = NULL, response = NULL, ncat = NULL,
                      bs.trend = F, y.time = NULL, knots.bt.block = NULL, block.no = NULL, bs.df = NULL,
                      corr.y = F,
                      alpha.prior = NULL, beta.prior = NULL, eta.prior = NULL, dc.prior = NULL, c1.prior = NULL,
                      rho.prior = NULL, hy.prior = NULL, ...){

  if(response == "ordinal"){
    if(is.null(ncat)){
      stop("ncat (number of categories) must be entered for ordinal response")
    }
  }
  nobs <- length(Y)

  # if(!baseline %in% Treat){
  #   stop("baseline treatment name is not in Treat")
  # }

  Treat <- gsub(" ", "\\_", Treat)
  # baseline <- gsub(" ", "\\_", baseline)

  # sort to make the treatment order the same every time for different participants
  Treat.name <- sort(unique(Treat))
  if (!is.null(ord.baseline.Treat)) {
    Treat.name <- c(ord.baseline.Treat, as.character(Treat.name[Treat.name != ord.baseline.Treat]))
  }
  # Treat.name <- Treat.name[Treat.name != baseline]

  nof1 = list(Y = Y, Treat = Treat, nobs = nobs, Treat.name = Treat.name, n.Treat = length(Treat.name), response = response)
  # nof1 = list(Y = Y, Treat = Treat, baseline = baseline, ncat = ncat, nobs = nobs, Treat.name = Treat.name, response = response)

  if (response == "ordinal") {
    nof1 <- c(nof1,
              ncat = ncat,
              ord.baseline.Treat = ord.baseline.Treat,
              ord.model          = ord.model)
  }

  # Treatment
  for(i in Treat.name){
    nam <- paste("Treat_", i, sep = "")
    nam <- assign(nam, as.numeric(Treat == i))
    nof1[[ paste("Treat_", i, sep = "")]] <- nam
  }

  # Indicators for whether adjusting for trend and correlation in the model
  nof1 <- c(nof1,
            bs.trend = bs.trend,
            corr.y   = corr.y)

  # for splines
  # Default knots at end of each block
  if (bs.trend){

    # center y.time
    cent.y.time <- y.time - mean(y.time, na.rm=T)

    # currently only have the option: knots at the end of each block
    if (knots.bt.block){

      knots <- cent.y.time[cumsum(rle(block.no)$lengths)]
      # remove knot at the end of last block
      knots <- knots[-length(knots)]
      bs.design.matrix <- bs(cent.y.time, knots = knots, ...)

    } else {

      bs.design.matrix <- bs(cent.y.time, df = bs.df, ...)
    }


    nof1$bs_df <- ncol(bs.design.matrix)

    # Save the columns in the bs.design.matrix in nof1 as bs*
    for (i in 1:(nof1$bs_df)){
      nof1[[paste0("bs", i)]] <- bs.design.matrix[, i]
    }

  }

  # if(!is.null(knots)){
  #   cen.knots <- (knots - mean(Time, na.rm = TRUE))/ sd(Time, na.rm = TRUE)
  #   BS <- bs(cen.Time, knots = cen.knots)
  #   nof1$BS <- BS
  #   nof1$knots <- knots
  # }

  # for correlated model, not used
  # if(!is.null(y.time)){
  #   cen.Time <- (y.time - mean(y.time, na.rm = TRUE)) / sd(y.time, na.rm = TRUE)
  #   nof1$y.time = cen.Time
  # }

  prior.param <- list(response = response, dc.prior = dc.prior, c1.prior = c1.prior, alpha.prior = alpha.prior, beta.prior = beta.prior, eta.prior = eta.prior, hy.prior = hy.prior, rho.prior = rho.prior)
  # if (response == "normal"){
  #   prior.data <- nof1.prior.default(prior.param, Y)
  # } else {
  prior.data <- nof1.prior.default(prior.param)
  # }

  nof1 <- c(nof1, prior.data)

  code <- nof1.rjags(nof1)
  nof1$code <- code
  # cat(code)

  class(nof1) <- "nof1.data"
  nof1
}

#' @export

# ID must be a complete vector with no missing values
nof1.ma.data <- function(Y, Treat, baseline.Treat, ID, response,
                         model.linkfunc = NULL, model.intcpt = "fixed", model.slp = "random",
                         ord.ncat = NULL, ord.model, ord.parallel = NULL,
                         covariates,
                         spline.trend = F, trend.type, y.time = NULL, knots = NULL, trend.df = NULL,
                         step.trend = F, y.step = NULL,
                         corr.y = F,
                         alpha.prior = NULL, beta.prior = NULL, eta.prior = NULL, dc.prior = NULL, c1.prior = NULL,
                         rho.prior = NULL, hy.prior = NULL, ...) {

  # Y: must be sorted according to ID, day
  # ID: same ID should stick together
  # model.intcpt: "random" or "fixed". Currently, only work for normal outcome.
  # covariates: lists of covariates, categorical covariates must be of factor type
  #             need to be careful about the covariate names, both treatment and covariate names are used to identify initial values
  # ord.ncat: number of levels in ordinal outcome
  # ord.parallel: logical. Whether the adjacent category model is restricted (same as parameter parallel in vglm).
  # knots: specify knots rather than knots between blocks bc people might have different block break points
  # step.trend: do step functions for each period
  # y.step: same length as y, corresponding to the periods, start from 1, integer values

  Treat.name <- sort(unique(Treat))
  if (!is.null(baseline.Treat)) {
    Treat.name <- c(baseline.Treat, as.character(Treat.name[Treat.name != baseline.Treat]))
  }
  # Relevel the treatment
  levels(Treat) <- Treat.name

  rle.ID  <- rle(ID)
  nobs.ID <- rle.ID$lengths
  uniq.ID <- rle.ID$values
  n.ID    <- length(nobs.ID)

  max.obs.ID <- max(nobs.ID)

  # uniq.ID has correspondence with nobs.ID
  nof1 <- list(Y.long     = Y,
               Treat      = Treat,
               ID         = ID,
               uniq.ID    = uniq.ID,
               nobs.ID    = nobs.ID,
               n.ID       = n.ID,
               Treat.name = Treat.name,
               n.Treat    = length(Treat.name),
               response   = response,
               model.linkfunc = model.linkfunc,
               model.intcpt   = model.intcpt,
               model.slp      = model.slp)

  if (response == "ordinal") {
    nof1 <- c(nof1,
              ord.ncat     = ord.ncat,
              ord.model    = ord.model,
              ord.parallel = ord.parallel)
  }

  # Generate outcome and treatment matrix
  y.matrix <- matrix(NA, nrow = max.obs.ID, ncol = n.ID)
  for (Treat.name.i in 1:length(Treat.name)) {
    assign(paste0("Treat_", Treat.name[Treat.name.i]), NULL)
  }

  for (ID.i in 1:n.ID) {
    y.matrix[1:nobs.ID[ID.i], ID.i] <- Y[ID == uniq.ID[ID.i]]

    for (Treat.name.i in 1:length(Treat.name)) {
      assign(paste0("Treat_", Treat.name[Treat.name.i]),
             cbind(get(paste0("Treat_", Treat.name[Treat.name.i])),
                   c(as.numeric(Treat[ID == uniq.ID[ID.i]] == Treat.name[Treat.name.i]), rep(NA, max.obs.ID - nobs.ID[ID.i]))))
    }
  }

  nof1$Y <- y.matrix
  for (Treat.name.i in 1:length(Treat.name)) {
    nof1[[paste0("Treat_", Treat.name[Treat.name.i])]] <- get(paste0("Treat_", Treat.name[Treat.name.i]))
  }

  # Covariates
  if (!is.null(covariates)) {
    names.covariates <- names(covariates)
    nof1$names.covariates      <- NULL
    nof1$names.long.covariates <- NULL

    for (covariates.i in 1:length(covariates)) {

      nof1$names.long.covariates <- c(nof1$names.long.covariates, paste0(names.covariates[covariates.i], "_long"))
      nof1[[paste0(names.covariates[covariates.i], "_long")]] <- covariates[[covariates.i]]

      # continuous covariates
      # NEED TO TEST CODE ON CONTINUOUS COVARIATES
      if (!is.factor(covariates[[covariates.i]])) {
        nof1[[names.covariates[covariates.i]]] <- NULL
        for (ID.i in 1:n.ID) {
          nof1[[names.covariates[covariates.i]]] <- cbind(nof1[[names.covariates[covariates.i]]],
                                                          c(covariates[[covariates.i]][ID == uniq.ID[ID.i]], rep(NA, max.obs.ID - nobs.ID[ID.i])))
        }
        nof1$names.covariates <- c(nof1$names.covariates, names.covariates[covariates.i])

      } else { # categorical covariates

        tmp.lvls <- levels(covariates[[covariates.i]])

        # For each level, create a covariate matrix with each column for a participant
        for (tmp.lvls.i in 2:length(tmp.lvls)) {

          tmp.cov.name <- paste0(names.covariates[covariates.i], "_", tmp.lvls.i)
          nof1[[tmp.cov.name]] <- NULL
          # add covariate name to the vector
          nof1$names.covariates <- c(nof1$names.covariates, tmp.cov.name)

          for (ID.i in 1:n.ID) {
            nof1[[tmp.cov.name]] <- cbind(nof1[[tmp.cov.name]],
                                          c(as.numeric(covariates[[covariates.i]][ID == uniq.ID[ID.i]] == tmp.lvls[tmp.lvls.i]), rep(NA, max.obs.ID - nobs.ID[ID.i])))
          }

        } # tmp.lvls.i
      } # categorical covariates
    } # covariate.i
  } # covariates exist

  # Indicators for whether adjusting for trend and correlation in the model
  nof1 <- c(nof1,
            spline.trend = spline.trend,
            step.trend = step.trend,
            corr.y     = corr.y)

  # for splines
  # Default knots at end of each block
  if (spline.trend) {

    # center y.time
    uniq.y.time <- sort(unique(y.time))
    cent.y.time <- uniq.y.time - mean(uniq.y.time, na.rm = T)

    # currently only have the option: knots at the end of each block
    if (!is.null(trend.df)) {

      spline.design.matrix <- bs(cent.y.time, df = bs.df, ...)
      spline.long.design.matrix <- spline.design.matrix[y.time, ]

    } else {

      cent.knots <- knots - mean(uniq.y.time)

      if (trend.type == "bs") {
        spline.design.matrix <- bs(cent.y.time, knots = cent.knots, ...)
      } else if (trend.type == "ns") {
        spline.design.matrix <- ns(cent.y.time, knots = cent.knots, ...)
      }

      spline.long.design.matrix <- spline.design.matrix[y.time, ]

    }

    # remove the column if there is a single value in the column on non-missing dates
    # we still have the same knots, but the coefficients for the columns with a single value will be 0
    # save in these temporary matrix because the column index will change due to removing columns
    tmp.spline.design.matrix <- NULL
    tmp.spline.long.design.matrix <- NULL

    for (i in 1:ncol(spline.design.matrix)) {

      if (length(unique(spline.long.design.matrix[!is.na(Y), i])) != 1) {
        tmp.spline.design.matrix <- cbind(tmp.spline.design.matrix, spline.design.matrix[, i])
        tmp.spline.long.design.matrix <- cbind(tmp.spline.long.design.matrix, spline.long.design.matrix[, i])
      }
    }

    spline.design.matrix      <- tmp.spline.design.matrix
    spline.long.design.matrix <- tmp.spline.long.design.matrix
    nof1$spline_df <- ncol(spline.design.matrix)

    # Save the columns in the spline.design.matrix in nof1 as bs*
    for (i in 1:(nof1$spline_df)){

      nof1[[paste0("spline", i)]] <- spline.design.matrix[, i]
      # used to find initial values
      nof1[[paste0("spline.long", i)]] <- spline.long.design.matrix[, i]

    }

    # create time index matrix incase there are multiple measurements on the same day
    time <- matrix(NA, nrow = max.obs.ID, ncol = n.ID)
    for (ID.i in 1:n.ID) {
      time[1:nobs.ID[ID.i], ID.i] <- y.time[ID == uniq.ID[ID.i]]
    }
    nof1$time <- time
  }

  # for step functions
  if (step.trend) {

    nof1$n.steps <- 2:max(y.step)
    nof1$y.step  <- y.step
    # produce steps from the second period
    for (step.i in nof1$n.steps) {
      assign(paste0("step", step.i), NULL)

      for (ID.i in 1:n.ID) {
        assign(paste0("step", step.i),
               cbind(get(paste0("step", step.i)),
                     c(as.numeric(y.step[ID == uniq.ID[ID.i]] == step.i), rep(NA, max.obs.ID - nobs.ID[ID.i]))))

      }
      nof1[[paste0("step", step.i)]] <- get(paste0("step", step.i))
    } # for (step.i in 2:max(y.step)) {

  } # if (step.trend)

  prior.param <- list(response = response, dc.prior = dc.prior, c1.prior = c1.prior, alpha.prior = alpha.prior, beta.prior = beta.prior, eta.prior = eta.prior, hy.prior = hy.prior, rho.prior = rho.prior)
  prior.data <- nof1.prior.default(prior.param)
  nof1 <- c(nof1, prior.data)

  code <- nof1.ma.rjags(nof1)
  nof1$code <- code
  # cat(code)

  class(nof1) <- "nof1.data"
  return(nof1)
}

#' Make an N of 1 object containing data, priors, and a jags model file for (network) meta-analyses
#'
#' @param Y Outcome of the study. This should be a vector with \code{NA}'s included in time order.
#' @param Treat Treatment indicator vector with the same length as the outcome. Can be character or numeric.
#' @param baseline.Treat Name of the reference treatment.
#' @param ID Participant ID vector with the same length as the outcome.
#' @param response Type of outcome. Can be "normal" for continuous outcome, "binomial" for binary outcome,
#' "poisson" for count outcome, or "ordinal" for ordinal or nominal outcome.
#' @param model.linkfunc Link function in the model.
#' @param model.intcpt Form of intercept.
#' @param model.slp Form of slope. For link function and the forms of intercept and slopes, currently implemented for
#' 1) normal response: "identity"-"fixed/random"-"random", 2) poisson response: "log"-"fixed"-"random",
#' 3) binomial response: "log/logit"-"fixed"-"random".
#' @param ord.ncat Number of categories in ordinal response. The parameters for ordinal outcomes need to be tested.
#' @param ord.model Used for ordinal outcome to pick the model. Can be "cumulative" for cumulative probability model
#' or "acat" for adjacent category model.
#' @param ord.parallel Whether or not the comparisons between categories or cumulative probabilities are parallel.
#' @param strata.cov Only applicable when random intercept model is used. Stratification covariates used during randomization.
#' Should be a data frame with the columns being the covariates and the number of rows should be equal to the length of the
#' outcome. Must be of factor type if categorical covariates.
#' @param adjust.strata.cov Only applicable when random intercept model is used. A vector with each element taking on possible
#' values from \code{c("fixed", "random")} indicating how the stratification covariate should be adjusted in the model.
#' @param lvl2.cov Participant level covariates for heterogeneous treatment effects.
#' Should be a data frame with the columns being the covariates and the number of rows should be equal to the length of the outcome.
#' For fixed-intercept model, participant level covariates will have interaction with treatment (slope) because
#' there are fixed-intercepts adjusting for each participant in the model; for random-intercept model, participant level covariates
#'  will have interaction with all treatment-specific intercepts.
#' @param spline.trend Indicator for whether the model should adjust for trend using splines. The default
#' is \code{F}.
#' @param trend.type "bs" for basis splines or "ns" for natural splines.
#' @param y.time Time when the outcome is measured. Required when adjusting for trend or correlation
#' (\code{spline.trend}/\code{step.trend}/\code{corr.y} is \code{TRUE}). Should be a vector of the same length as the outcome.
#' @param knots Knots in \code{y.time} when adjusting for trend using splines. For \code{trend.type = "bs"}, knots should be set
#' at the end of each block except for the last block; for or \code{trend.type = "ns"}, knots should be set at the end of each
#' period except for the last period.
#' @param trend.df Degrees of freedom for modeling splines when knots are not set.
#' @param step.trend Indicator for whether to adjust for trend using step functions for each period.
#' @param y.step Period number corresponding to the outcome. Should be a vector of the same length of the outcome.
#' @param corr.y Indicator for whether the correlation among measurements shoule be modeled. The default is
#' \code{F}.
#' @param alpha.prior Prior for the fixed-intercepts.
#' @param beta.prior Prior for random intercepts or slopes.
#' @param eta.prior Prior for modelling spline terms or heterogeneous treatment effects.
#' @param dc.prior Prior for the length between cutpoints. Used only for ordinal logistic models.
#' @param c1.prior Prior for the first cutpoint. Used only for ordinal logistic models.
#' @param rho.prior Prior for the correlation between random effects or the correlation among repeated measurements on each participant.
#' @param hy.prior Prior for the variance of the residual errors for normal response,
#' the standard deviation of random slopes for binary outcome, the variance of random effects for other types of outcomes.

# Supports uniform, gamma, and half normal for
# normal and binomial response and wishart for multinomial response. It should be a list of length 3,
# where first element should be the distribution (one of dunif, dgamma, dhnorm, dwish) and the next
# two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give
# uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter. For wishart
# distribution, the last two parameter would be the scale matrix and the degrees of freedom.

#' @export

# ID must be a complete vector with no missing values
nof1.nma.data <- function(Y, Treat, baseline.Treat, ID, response,
                          model.linkfunc = NULL, model.intcpt = "fixed", model.slp = "random",
                          ord.ncat = NULL, ord.model, ord.parallel = NULL,
                          strata.cov = NULL, adjust.strata.cov = NULL, lvl2.cov = NULL,
                          spline.trend = F, trend.type, y.time = NULL, knots = NULL, trend.df = NULL,
                          step.trend = F, y.step = NULL,
                          corr.y = F,
                          alpha.prior = NULL, beta.prior = NULL, eta.prior = NULL, dc.prior = NULL, c1.prior = NULL,
                          rho.prior = NULL, hy.prior = NULL,
                          na.rm = T, ...) {

  data.long <- data.frame(ID, Treat, Y)

  # Treatment levels, removed treatment with all NA's in outcome Y
  Treat.name <- as.character(sort(unique(data.long$Treat[!is.na(data.long$Y)])))
  if (!is.null(baseline.Treat)) {
    Treat.name <- c(baseline.Treat, as.character(Treat.name[Treat.name != baseline.Treat]))
  }
  # Relevel the treatment
  data.long <- data.long %>%
    mutate(Treat = factor(Treat, Treat.name))

  # summary by ID
  if (na.rm) {
    summ.ID <- data.long %>%
      group_by(ID) %>%
      summarise(nobs = sum(!is.na(Y)),
                n_Treat = length(unique(Treat[!is.na(Y)]))) %>%
      arrange(n_Treat) # order by the # of treatment
  } else {
    summ.ID <- data.long %>%
      group_by(ID) %>%
      summarise(nobs = n(),
                n_Treat = length(unique(Treat))) %>%
      arrange(n_Treat) # order by the # of treatment
  }

  max.obs.ID <- max(summ.ID$nobs)

  summ.nID.perTreat <- summ.ID %>%
    group_by(n_Treat) %>%
    summarise(n_ID_perNTreat = n()) %>%
    arrange(n_Treat)
  max.Treat.ID <- max(summ.nID.perTreat$n_Treat)

  # always have 1:maximum number of treatment per participant
  tmp.n.Treat <- data.frame(n_Treat = 1:max.Treat.ID)
  summ.nID.perTreat <- tmp.n.Treat %>%
    left_join(summ.nID.perTreat, by = "n_Treat")
  summ.nID.perTreat <- summ.nID.perTreat %>%
    mutate(n_ID_perNTreat = ifelse(is.na(n_ID_perNTreat), 0, n_ID_perNTreat))

  # summ.ID ordered by the number of treatment each participant has
  nof1 <- list(data.long  = data.long,
               summ.ID    = summ.ID,
               nobs.ID    = summ.ID$nobs,
               summ.nID.perTreat = summ.nID.perTreat,
               Treat.name = Treat.name,
               response   = response,
               model.linkfunc = model.linkfunc,
               model.intcpt   = model.intcpt,
               model.slp      = model.slp)

  # Generate outcome, treatment matrix, and treatment indicator matrix
  Y.matrix     <- matrix(NA, nrow = max.obs.ID, ncol = nrow(summ.ID))
  Treat.matrix <- matrix(NA, nrow = max.obs.ID, ncol = nrow(summ.ID))
  uniq.Treat.matrix <- matrix(NA, nrow = max.Treat.ID, ncol = nrow(summ.ID))
  for (Treat.i in 1:nrow(summ.nID.perTreat)) {
    nof1[[paste0("Treat.", Treat.i)]] <- matrix(NA, nrow = max.obs.ID, ncol = nrow(summ.ID))
  }

  if (na.rm) {
    for (ID.i in 1:nrow(summ.ID)) {
      Y.matrix[1:summ.ID$nobs[ID.i], ID.i] <- data.long$Y[(ID == summ.ID$ID[ID.i]) & (!is.na(data.long$Y))]
      Treat.matrix[1:summ.ID$nobs[ID.i], ID.i] <- as.numeric(data.long$Treat[(ID == summ.ID$ID[ID.i]) & (!is.na(data.long$Y))])

      tmp.Treat <- sort(unique(Treat.matrix[1:summ.ID$nobs[ID.i], ID.i]))
      uniq.Treat.matrix[1:summ.ID$n_Treat[ID.i], ID.i] <- tmp.Treat
      for (Treat.i in 1:summ.ID$n_Treat[ID.i]) {
        nof1[[paste0("Treat.", Treat.i)]][1:summ.ID$nobs[ID.i], ID.i] <- as.numeric(Treat.matrix[1:summ.ID$nobs[ID.i], ID.i] == tmp.Treat[Treat.i])
      }
    }
  } else { # na.rm

    for (ID.i in 1:nrow(summ.ID)) {
      Y.matrix[1:summ.ID$nobs[ID.i], ID.i] <- data.long$Y[ID == summ.ID$ID[ID.i]]
      Treat.matrix[1:summ.ID$nobs[ID.i], ID.i] <- as.numeric(data.long$Treat[ID == summ.ID$ID[ID.i]])

      tmp.Treat <- sort(unique(Treat.matrix[1:summ.ID$nobs[ID.i], ID.i]))
      uniq.Treat.matrix[1:summ.ID$n_Treat[ID.i], ID.i] <- tmp.Treat
      for (Treat.i in 1:summ.ID$n_Treat[ID.i]) {
        nof1[[paste0("Treat.", Treat.i)]][1:summ.ID$nobs[ID.i], ID.i] <- as.numeric(Treat.matrix[1:summ.ID$nobs[ID.i], ID.i] == tmp.Treat[Treat.i])
      }
    }
  }

  nof1$Y.matrix <- Y.matrix
  nof1$Treat.matrix <- Treat.matrix
  nof1$uniq.Treat.matrix <- uniq.Treat.matrix

  # Indicators for whether adjusting for trend and correlation in the model
  nof1 <- c(nof1,
            spline.trend = spline.trend,
            step.trend   = step.trend,
            corr.y       = corr.y)

  # splines
  if (spline.trend) {

    nof1$data.long$Y_time <- y.time
    # center y.time
    uniq.y.time <- sort(unique(y.time))
    cent.y.time <- uniq.y.time - mean(uniq.y.time, na.rm = T)

    # currently only have the option: knots at the end of each block
    if (!is.null(trend.df)) {

      spline.matrix <- bs(cent.y.time, df = bs.df, ...)
      spline.long.matrix <- spline.matrix[y.time, ]

    } else {

      cent.knots <- knots - mean(uniq.y.time)

      if (trend.type == "bs") {
        spline.matrix <- bs(cent.y.time, knots = cent.knots, ...)
      } else if (trend.type == "ns") {
        spline.matrix <- ns(cent.y.time, knots = cent.knots, ...)
      }

      spline.long.matrix <- spline.matrix[y.time, ]

    }

    # remove the column if there is a single value in the column on non-missing dates
    # we still have the same knots, but the coefficients for the columns with a single value will be 0
    # save in these temporary matrix because the column index will change due to removing columns
    tmp.spline.matrix <- NULL
    tmp.spline.long.matrix <- NULL

    for (i in 1:ncol(spline.matrix)) {

      if (length(unique(spline.long.matrix[!is.na(Y), i])) != 1) {
        tmp.spline.matrix <- cbind(tmp.spline.matrix, spline.matrix[, i])
        tmp.spline.long.matrix <- cbind(tmp.spline.long.matrix, spline.long.matrix[, i])
      }
    }

    spline.matrix      <- tmp.spline.matrix
    spline.long.matrix <- tmp.spline.long.matrix

    # Save spline.matrix in nof1
    nof1$spline.df          <- ncol(spline.matrix)
    nof1$spline.matrix      <- spline.matrix
    nof1$spline.long.matrix <- spline.long.matrix

    # create time index matrix incase there are multiple measurements on the same day
    time.matrix <- matrix(NA, nrow = max.obs.ID, ncol = nrow(summ.ID))
    for (ID.i in 1:nrow(summ.ID)) {
      time.matrix[1:summ.ID$nobs[ID.i], ID.i] <- nof1$data.long$Y_time[(ID == summ.ID$ID[ID.i]) & (!is.na(nof1$data.long$Y))]
    }
    nof1$time.matrix <- time.matrix
  }

  # Step function
  if (step.trend) {

    nof1$data.long$Y_step <- y.step
    # NEED TO TEST WHEN THERE ARE NA's in OUTCOME Y
    step.matrix <- diag(length(unique(y.step[!is.na(Y)])))[, -1]
    step.long.matrix <- step.matrix[y.step, ]

    # Save spline.matrix in nof1
    nof1$step.df          <- ncol(step.matrix)
    nof1$step.matrix      <- step.matrix
    nof1$step.long.matrix <- step.long.matrix

    # create period matrix
    period.matrix <- matrix(NA, nrow = max.obs.ID, ncol = nrow(summ.ID))
    for (ID.i in 1:nrow(summ.ID)) {
      period.matrix[1:summ.ID$nobs[ID.i], ID.i] <- nof1$data.long$Y_step[(ID == summ.ID$ID[ID.i]) & (!is.na(nof1$data.long$Y))]
    }
    nof1$period.matrix <- period.matrix

  }

  # serial correlation
  if (corr.y) {
    nof1$data.long$Y_time <- y.time

    # create time index matrix to account for non consecutive outcome y
    time.matrix <- matrix(NA, nrow = max.obs.ID, ncol = nrow(summ.ID))
    for (ID.i in 1:nrow(summ.ID)) {
      time.matrix[1:summ.ID$nobs[ID.i], ID.i] <- nof1$data.long$Y_time[(ID == summ.ID$ID[ID.i]) & (!is.na(nof1$data.long$Y))]
    }
    nof1$time.matrix <- time.matrix
  }

  # stratification covariates used during randomization
  if (!is.null(strata.cov)) {

    if (model.intcpt == "fixed") {
      stop("Fixed intercept model do not need to adjust for stratification covariates used in randomization!")

    } else if (model.intcpt == "random") {

      nof1$data.long <- cbind(data.long, strata.cov)
      colnames.cov <- paste0("strata.cov", 1:ncol(strata.cov))
      colnames(nof1$data.long)[(ncol(nof1$data.long) - ncol(strata.cov) + 1):ncol(nof1$data.long)] <- colnames.cov

      strata.cov.wizId <- nof1$data.long[!duplicated(nof1$data.long[, c("ID", colnames.cov)]), c("ID", colnames.cov)]

      # need to justify the covariates are participant level covariate
      if (nrow(strata.cov.wizId) != length(nof1$nobs.ID)) {
        stop("strata.cov must be participant level covariates!")
      }

      nof1$summ.ID <- nof1$summ.ID %>%
        left_join(strata.cov.wizId, by = "ID")
      nonDup.lvl2.cov <- data.frame(nof1$summ.ID[, colnames.cov])
      for (cov.i in 1:ncol(nonDup.lvl2.cov)) {

        # must be factor covariates
        # create dummy variables if fixed effects
        # create ordinal variables if random effects
        if (is.factor(nonDup.lvl2.cov[, cov.i])) {

          if (adjust.strata.cov[cov.i] == "fixed") {
            tmp.cov.lvls <- levels(nonDup.lvl2.cov[, cov.i])
            for (cov.lvls.i in 2:length(tmp.cov.lvls)) {
              nof1$fixed.strata.cov.matrix <- rbind(nof1$fixed.strata.cov.matrix,
                                                    as.numeric(nonDup.lvl2.cov[, cov.i] == tmp.cov.lvls[cov.lvls.i]))
            }
          } else if (adjust.strata.cov[cov.i] == "random") {
            nof1$random.strata.cov.matrix <- rbind(nof1$random.strata.cov.matrix,
                                                   as.numeric(nonDup.lvl2.cov[, cov.i]))
          }

        } else {
          stop("strata.cov must be factor!")
        }

      } # for (cov.i in 1:ncol(nonDup.lvl2.cov)) {

      nof1$n.strata.cov       <- ncol(strata.cov)
      nof1$adjust.strata.cov  <- adjust.strata.cov
      nof1$n.fixed.cov.model  <- nrow(nof1$fixed.strata.cov.matrix)
      nof1$n.random.cov.model <- nrow(nof1$random.strata.cov.matrix)
    } # random intercept

  } # strata.cov

  # participant level covariates
  if (!is.null(lvl2.cov)) {

    nof1$data.long <- cbind(nof1$data.long, lvl2.cov)
    colnames.cov <- paste0("lvl2.cov", 1:ncol(lvl2.cov))
    colnames(nof1$data.long)[(ncol(nof1$data.long) - ncol(lvl2.cov) + 1):ncol(nof1$data.long)] <- colnames.cov

    lvl2.cov.wizId <- nof1$data.long[!duplicated(nof1$data.long[, c("ID", colnames.cov)]), c("ID", colnames.cov)]

    # need to justify the covariates are participant level covariate
    if (nrow(lvl2.cov.wizId) != length(nof1$nobs.ID)) {
      stop("lvl2.cov must be second level covariates")
    }

    nof1$summ.ID <- nof1$summ.ID %>%
      left_join(lvl2.cov.wizId, by = "ID")
    nonDup.lvl2.cov <- data.frame(nof1$summ.ID[, colnames.cov])
    for (cov.i in 1:ncol(nonDup.lvl2.cov)) {

      # for factor covariates, create dummy variables
      if (is.factor(nonDup.lvl2.cov[, cov.i])) {

        tmp.cov.lvls <- levels(nonDup.lvl2.cov[, cov.i])
        for (cov.lvls.i in 2:length(tmp.cov.lvls)) {
          nof1$cov.matrix <- rbind(nof1$cov.matrix,
                                   as.numeric(nonDup.lvl2.cov[, cov.i] == tmp.cov.lvls[cov.lvls.i]))
        }

      } else { # continuous covariates
        nof1$cov.matrix <- rbind(nof1$cov.matrix,
                                 nonDup.lvl2.cov[, cov.i])
      }

    } # for (cov.i in 1:ncol(nonDup.lvl2.cov)) {

    nof1$n.cov <- ncol(lvl2.cov)
    nof1$n.cov.model <- nrow(nof1$cov.matrix)

    # create actual ordered treatment matrix
    # because nof1$Treat.x may correspond to different treatment, the interaction term is with the actual treatment
    # for (Treat.i in 1:length(nof1$Treat.name)) {
    #   nof1[[paste0("Treat.order.", Treat.i)]] <- matrix(NA, nrow = max.obs.ID, ncol = nrow(summ.ID))
    # }
    #
    # for (ID.i in 1:nrow(summ.ID)) {
    #   for (Treat.i in 1:length(nof1$Treat.name)) {
    #     nof1[[paste0("Treat.order.", Treat.i)]][1:summ.ID$nobs[ID.i], ID.i] <- as.numeric(Treat.matrix[1:summ.ID$nobs[ID.i], ID.i] == Treat.i)
    #   }
    # }

  } # if (!is.null(lvl2.cov)) {

  prior.param <- list(response = response, dc.prior = dc.prior, c1.prior = c1.prior, alpha.prior = alpha.prior, beta.prior = beta.prior, eta.prior = eta.prior, hy.prior = hy.prior, rho.prior = rho.prior)
  prior.data <- nof1.prior.default(prior.param)
  nof1 <- c(nof1, prior.data)

  code <- nof1.nma.rjags(nof1)
  nof1$code <- code
  # cat(code)

  class(nof1) <- "nof1.data"
  return(nof1)

}


