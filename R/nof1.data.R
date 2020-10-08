#' Make an N of 1 object containing data, priors, and a jags model file
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

    # Save the columns in the bs.design.matrix in nof1 as Treat_*
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


# ID must be a complete vector with no missing values
nof1.ma.data <- function(Y, Treat, baseline.Treat, ID, response, model.intcpt = "fixed", model.linkfunc = NULL,
                         ord.ncat = NULL, ord.model, ord.parallel = NULL,
                         covariates,
                         bs.trend = F, y.time = NULL, knots.bt.block = NULL, block.no = NULL, bs.df = NULL,
                         corr.y = F,
                         alpha.prior = NULL, beta.prior = NULL, eta.prior = NULL, dc.prior = NULL, c1.prior = NULL,
                         rho.prior = NULL, hy.prior = NULL, ...) {

  # ID: same ID should stick together
  # model.intcpt: "random" or "fixed". Currently, only work for normal outcome.
  # covariates: lists of covariates, categorical covariates must be of factor type
  #             need to be careful about the covariate names, both treatment and covariate names are used to identify initial values
  # ord.ncat: number of levels in ordinal outcome
  # ord.parallel: logical. Whether the adjacent category model is restricted (same as parameter parallel in vglm).

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
               model.intcpt   = model.intcpt,
               model.linkfunc = model.linkfunc)

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

  prior.param <- list(response = response, dc.prior = dc.prior, c1.prior = c1.prior, alpha.prior = alpha.prior, beta.prior = beta.prior, eta.prior = eta.prior, hy.prior = hy.prior, rho.prior = rho.prior)
  prior.data <- nof1.prior.default(prior.param)
  nof1 <- c(nof1, prior.data)

  code <- nof1.ma.rjags(nof1)
  nof1$code <- code
  # cat(code)

  class(nof1) <- "nof1.data"
  return(nof1)
}
