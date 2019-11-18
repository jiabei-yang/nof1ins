#' Function to present a summary of our results
#'
#' A neat function to summarize the results.
#'
#' @param result A list which contains data file created
#' by \code{nof1.data} and the result file created by \code{nof1.run}
#' @param alpha The alpha value for the confidence interval. If no value is
#' entered will give the 95\% confidence interval.
#' @return The function computes and returns a list of summary statistics of 
#' the raw data and the fitted model.
#' \item{raw.y.mean}{The raw mean of the outcome for each treatment}
#' \item{raw.y.median}{The raw median of the outcome for each treatment}
#' \item{post.coef.mean}{The posterior mean of the coefficient for each treatment}
#' \item{post.coef.median}{The posterior median of the coefficient for each treatment}
#' \item{post.y.mean}{The posterior mean of the outcome for each treatment}
#' \item{post.y.median}{The posterior median of the outcome for each treatment}
#' \item{post.coef.ci}{The credible interval of the coefficient for each treatment}
#' \item{post.y.ci}{The credible interval of the outcome for each treatment}
#' \item{comp.treat.post.coef}{The posterior quantiles of one coefficient minus the other 
#' when comparing two treatments}
#' \item{comp.treat.post.y}{The posterior quantiles of outcome minus the other when 
#' comparing two treatments}
#' \item{p.comp.coef.greater.0}{The posterior probability of one coefficient minus the other 
#' greater than 0}
#' @export

# \item{comp.treat.post.y}{The posterior quantiles of outcome minus the other when 
# comparing two treatments}
summarize_nof1 <- function(result, alpha = 0.05){
  
  treat.name <- result$nof1$Treat.name
  samples    <- do.call(rbind, result$samples)
  response   <- result$nof1$response
  
  raw.y.mean   <- NULL
  raw.y.median <- NULL
  post.coef.mean   <- NULL
  post.coef.median <- NULL
  post.y.mean   <- NULL
  post.y.median <- NULL
  post.coef.ci <- NULL
  post.y.ci    <- NULL
  
  for (i in treat.name){
    raw.y.mean   <- c(raw.y.mean, mean(result$nof1$Y[result$nof1$Treat == i], na.rm = TRUE))
    raw.y.median <- c(raw.y.median, median(result$nof1$Y[result$nof1$Treat == i], na.rm = TRUE))
    
    col.treat.name <- paste0("beta_", i)
    post.coef <- samples[, col.treat.name]
    post.coef.mean   <- c(post.coef.mean, mean(post.coef))
    post.coef.median <- c(post.coef.median, median(post.coef))
    
    post.y <- link_function(post.coef, response)
    post.y.mean   <- c(post.y.mean, mean(post.y))
    post.y.median <- c(post.y.median, median(post.y))
    
    post.coef.ci <- rbind(post.coef.ci, quantile(post.coef, c(alpha/2, 1-alpha/2)))
    post.y.ci    <- rbind(post.y.ci, quantile(post.y, c(alpha/2, 1-alpha/2)))
  }
  
  names(raw.y.mean) <- names(raw.y.median) <- names(post.coef.mean) <- names(post.coef.median) <- names(post.y.mean) <- names(post.y.median) <- treat.name
  rownames(post.coef.ci) <- rownames(post.y.ci) <- treat.name
  
  comp.treat.name <- t(utils::combn(treat.name, 2))
  comp.treat.post.coef   <- NULL
  p.comp.coef.greater.0 <- NULL
  for (i in 1:nrow(comp.treat.name)){
    
    col.treat.name.1 <- paste0("beta_", comp.treat.name[i, 1])
    col.treat.name.2 <- paste0("beta_", comp.treat.name[i, 2])
    post.coef.1 <- samples[, col.treat.name.1]
    post.coef.2 <- samples[, col.treat.name.2]
    comp.treat.post.coef <- rbind(comp.treat.post.coef, 
                                  quantile(post.coef.2 - post.coef.1, c(alpha/2, 0.5, 1-alpha/2)))
    
    p.comp.coef.greater.0 <- c(p.comp.coef.greater.0,
                                mean((post.coef.2 - post.coef.1) > 0))
    # The following does not make sense for outcome other than normal
    # post.y.1 <- link_function(post.coef.1, response)
    # post.y.2 <- link_function(post.coef.2, response)
    # comp.treat.post.y <- rbind(comp.treat.post.y, 
    #                            quantile(post.y.2 - post.y.1, c(alpha/2, 0.5, 1-alpha/2)))
  }
  
  rownames(comp.treat.post.coef) <- paste(comp.treat.name[, 2], comp.treat.name[, 1], sep = "_minus_")
    # rownames(comp.treat.post.y) <- 
  names(p.comp.coef.greater.0)  <- paste(comp.treat.name[, 2], comp.treat.name[, 1], sep = "_minus_")
  
  summ <- list(raw.y.mean   = raw.y.mean, 
               raw.y.median = raw.y.median,
               post.coef.mean   = post.coef.mean,
               post.coef.median = post.coef.median,
               post.y.mean   = post.y.mean,
               post.y.median = post.y.median,
               post.coef.ci  = post.coef.ci,
               post.y.ci     = post.y.ci, 
               comp.treat.post.coef  = comp.treat.post.coef,
               # comp.treat.post.y     = comp.treat.post.y,
               p.comp.coef.greater.0 = p.comp.coef.greater.0)
  return(summ)
  
}

#' time series plot across different interventions
#' 
#' @param nof1 nof1 object created using nof1.data
#' @param trial.start Start time of the trial specified with \code{timestamp.format}.
#' @param trial.end End time of the trial specified with \code{timestamp.format}.
#' @param timestamp.format Format of the \code{trial.start} and \code{trial.end}.
#' @examples
#' Y <- laughter$Y
#' Treat <- laughter$Treat
#' nof1 <- nof1.data(Y, Treat, ncat = 11, baseline = "Usual Routine", response = "ordinal")
#' timestamp <- seq(as.Date('2015-01-01'),as.Date('2016-01-31'), length.out = length(Y))
#' time_series_plot(nof1, timestamp = timestamp, timestamp.format = "%m-%d-%Y", Outcome.name = "Stress")
#' @export

# @param x.name used to label x-axis time variable. Default is "time".
# @param y.name used to label y-axis outcome variable
# @param normal.response.range the range of the outcome if continuous; a vector of minimum and maximum
time_series_plot <- function(result, 
                             overlay.with.model = F, plot.by.treat = T,
                             trial.start = NULL, trial.end = NULL, timestamp.format = NULL){
  
  if (!is.null(trial.start)){
    
    trial.start <- as.Date(trial.start, timestamp.format)
    trial.end   <- as.Date(trial.end, timestamp.format)
    
    time <- seq(trial.start, trial.end, length = length(result$nof1$Y))
      
    data <- data.frame(Y = as.numeric(result$nof1$Y), 
                       Treatment = gsub("\\_", " ", result$nof1$Treat), 
                       time = time)

  } else {
    data <- data.frame(Y = as.numeric(result$nof1$Y), 
                       Treatment = gsub("\\_", " ", result$nof1$Treat), 
                       time = 1:length(result$nof1$Y))
  }

  # Now only normal has been checked
  if (result$nof1$response %in% c("normal", "poisson")) {
    
    if (plot.by.treat){
      fig <- ggplot(data, aes(x = time, Y, color = Treatment)) + 
        geom_point(na.rm = T) +
        facet_wrap(. ~ Treatment) + 
        theme_bw()
    } else{
      fig <- ggplot(data, aes(x = time, Y, color = Treatment)) + 
        geom_point(na.rm = T) + 
        theme_bw()
    }
    
    if (overlay.with.model) {
      
      data$model <- 0
      treat.name <- result$nof1$Treat.name
      col.treat   <- paste0("beta_", treat.name)
      samples    <- do.call(rbind, result$samples)
      median.beta <- apply(samples[, col.treat], 2, median)
      
      for (i in 1:length(treat.name)){
        data <- data %>%
          mutate(model = model + median.beta[i] * (Treatment == treat.name[i]))
      }
      
      if (result$nof1$bs.trend) {
        samples <- do.call(rbind, result$samples)
        col.bs  <- paste0("eta_", 1:result$nof1$bs_df) 
        median.eta <- apply(samples[, col.bs], 2, median)
        
        for (i in 1:result$nof1$bs_df) {
          bs.name <- paste0("bs", i)
          data <- data %>%
            mutate(model = model + median.eta[i] * result$nof1[[bs.name]])
        }

      }
      
      trt.lth <- rle(as.vector(data$Treatment))$length
      
      fig <- fig + 
        geom_line(data = data[1:cumsum(trt.lth)[1], ],
                  aes(x = time, y = model, linetype = "Fitted model", color = Treatment))
                  # color = brewer.pal(7, "Set1")[which(levels(data$Treatment) ==  data$Treatment[1])])
      for (i in 2:length(trt.lth)){
        fig <- fig + 
          geom_line(data = data[(cumsum(trt.lth)[i-1]+1):cumsum(trt.lth)[i], ],
                    aes(x = time, y = model, linetype = "Fitted model", color = Treatment))
                    # color = brewer.pal(7, "Set1")[which(levels(data$Treatment) ==  data$Treatment[cumsum(trt.lth)[i]])]
      }
      
      fig <- fig + 
        scale_linetype(guide = FALSE)
        
    }
    fig <- fig + 
      scale_color_brewer(palette = "Set1")

  } else { # other than normal, poisson
    
    fig <- ggplot(data, aes(x=date, Y, fill = Treatment)) + 
      geom_bar(stat = "identity")  + 
      facet_grid(. ~ Treatment) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line        = element_line(colour = "black")) + 
      geom_hline(data  = data2, 
                 aes(yintercept = x, linetype = "Mean"), 
                 color ="black") + 
      scale_y_continuous(breaks = 0:nof1$ncat, 
                         oob    = rescale_none, 
                         label  = c("Low", rep("", length = nof1$ncat -1), "High")) +
      scale_fill_manual(values = c("#adc2eb", "#ffb380")) + 
      scale_linetype_manual(name   = "", 
                            values = 1, 
                            guide  = guide_legend(override.aes = list(color = c("black"))))
    
  }
  
  # if (!is.null(x.name)){
  #   fig <- fig + 
  #     xlab(x.name) 
  #     # theme(axis.title.x = element_blank())
  # } 
  # 
  # if (!is.null(y.name)) {
  #   fig <- fig +
  #     ylab(y.name) 
  # }
  # 
  # if (!is.null(title)){
  #   fig <- fig + 
  #     ggtitle(title)
  # }
  
  fig
}



#' Frequency plot for raw data
#' 
#' @param nof1 nof1 object created using nof1.data
#' @param ... parameters to pass to \code{geom_bar} or \code{geom_histogram}.
#' @examples
#' Y <- laughter$Y
#' Treat <- laughter$Treat
#' nof1 <- nof1.data(Y, Treat, ncat = 11, baseline = "Usual Routine", response = "ordinal")
#' frequency_plot(nof1)
#' @export

# @param xlab x axis label
# @param title title name
frequency_plot <- function(nof1, ...){
  
  if(nof1$response %in% c("binomial", "ordinal")){
    
    data <- aggregate(nof1$Y, list(Y = nof1$Y, Treat = nof1$Treat), length)
    ggplot(data= data, aes(x= Y, y= x, fill=Treat, color = Treat)) +  
      geom_bar(stat="identity", position="dodge", alpha = 0.5, ...) + 
      ylab("Count") + 
      # xlim(0.5, nof1$ncat +0.5) + 
      scale_x_continuous(breaks = unique(data$Y),
                         label  = as.character(unique(data$Y))) + 
      theme_bw()  
  
  } else if(nof1$response %in% c("normal", "poisson")){
    
    data <- data.frame(Y = nof1$Y, Treat = nof1$Treat)
    ggplot(data, aes(x = Y, fill = Treat, color = Treat)) + 
      geom_histogram(position = "dodge", na.rm = T, alpha = 0.5, ...) + 
      ylab("Count") + 
      theme_bw()
  }
}

#' Stacked_percent_barplot for raw data (for ordinal or binomial data)
#'
#' @param nof1 nof1 object created using nof1.data
#' @examples
#' Y <- laughter$Y
#' Treat <- laughter$Treat
#' nof1 <- nof1.data(Y, Treat, ncat = 11, baseline = "Usual Routine", response = "ordinal")
#' stacked_percent_barplot(nof1)
#' @export

# @param title title name
stacked_percent_barplot <- function(nof1){
  
  if(nof1$response %in% c("binomial", "ordinal")){
    data <- aggregate(nof1$Y, list(Y = nof1$Y, Treat = nof1$Treat), length)
    
    ggplot(data, aes(fill= factor(Y, levels = unique(data$Y)), y= x, x= Treat)) + 
      geom_bar( stat="identity", position="fill") + 
      scale_y_continuous(labels = percent_format()) +  
      theme_bw() + 
      labs(x = "Treatment", y = "Percentage", fill = "Outcomes") +
      scale_fill_manual(values = 4:(3+length(unique(data$Y))), 
                        labels = unique(data$Y), 
                        drop   = FALSE)
  } else{
    stop("only works for binomial and ordinal data")
  }
}

#' Summary data table for nof1
#'
#' Provides a summary data table for the particular outcome in a particular dataset.
#'
#' @param nof1 nof1 object created using nof1.data
#' @return Gives a comprhensive table with several statistical values.
#' Each column indicates the value given. For a normal or poisson
#' response type the following are given:
#' \item{Treat}{The treatment recieved}
#' \item{mean}{The average value of the outcome}
#' \item{sd}{The standard deviation for the outcome}
#' \item{2.5\%}{2.5\% of the data are equal to or less than this value}
#' \item{50\%}{50\% of the data are equal to or less than this value}
#' \item{97.5\%}{97.5\% of the data are equal to or less than this value}
#'
#' For a binomial or ordinal response type, returns a table where first row
#' is each treatment and the following rows are the the number of data points
#' taken at each possible value.
#' @examples
#' Y <- laughter$Y
#' Treat <- laughter$Treat
#' nof1 <- nof1.data(Y, Treat, ncat = 11, baseline = "Usual Routine", response = "ordinal")
#' raw_table(nof1)
# @export

raw_table <- function(nof1){
  
  if(nof1$response %in% c("binomial", "ordinal")){
    table(nof1$Y, nof1$Treat)
  } else if(nof1$response %in% c("normal", "poisson")){
    raw_table <- aggregate(nof1$Y, list(Treat = nof1$Treat), mean, na.rm = T)
    colnames(raw_table)[2] <- "mean"
    cbind(raw_table, 
          sd = aggregate(nof1$Y, list(Treat = nof1$Treat), sd, na.rm = T)[,-1],
          aggregate(nof1$Y, list(Treat = nof1$Treat), quantile, na.rm = T, c(0.025, 0.5, 0.975))[,-1])
  }
}

#' Kernel density estimation plot
#'
#' Creates a kernel density estimation plot for a specific outcome
#'
#' @param result An object with information about the simulation. The object is
#' derived from the output of nof1.run
#' @param ... parameters to pass to \code{geom_histogram}.
#' @export

# @param bins The number of bins the histogram will contain. Default is 30.
# @param x_min The lower limit of the x-axis. Default is to set to zero
# @param x_max The upper limit of the x-axis. Default is to set the upper limit
# to the maximum value of the data inputed
# @param title The title of the graph
kernel_plot <- function(result, comp = T, ...){
  
  samples <- do.call(rbind, result$samples)
  # beta_variable <- exp(samples[,grep("beta", colnames(samples))])
  beta_variable <- samples[, grep("beta", colnames(samples))]
  data <- as.data.frame(beta_variable)
  
  response_type = result$nof1$response
  if (response_type == "binomial") {
    xlab = "Log Odds"
  } else if (response_type == "poisson") {
    xlab = "Log Risk"
  } else if (response_type == "normal") {
    xlab = "Outcome Mean"
  } else if (response_type == "ordinal") {
    xlab = "Log Odds"
  }
  
  if (!comp){
    data <- melt(data, id.vars = NULL, variable.name = "beta", value.name = "estimate")
    data <- data %>%
      mutate(beta = gsub("beta_", "Treatment: ", beta))
    
    ggplot(data, aes(x = estimate, color = beta, fill = beta)) + 
      geom_density(na.rm = TRUE, size = 0.5, alpha = 0) +
      geom_histogram(aes(y = ..density..), alpha = 0.5, na.rm = TRUE, ...) + 
      theme_bw() + 
      facet_wrap(. ~ beta) + 
      # xlim(xlim_value[1], xlim_value[2]) + 
      labs(x = xlab, y = "Density") + 
      theme(legend.title=element_blank())
    
  } else {
    
    comp.treat.name <- t(utils::combn(colnames(data), 2))
    comp.treat.data <- NULL
    for (i in 1:nrow(comp.treat.name)) {
      tmp.colnames <- gsub("beta_", "", comp.treat.name[i, ])
      tmp.colnames <- paste0(tmp.colnames[2], " - ", tmp.colnames[1])
      
      comp.treat.data <- cbind(comp.treat.data,
                               data[, comp.treat.name[i, 2]] -  data[, comp.treat.name[i, 1]])
      comp.treat.data <- data.frame(comp.treat.data)
      colnames(comp.treat.data)[ncol(comp.treat.data)] <- tmp.colnames
    }
    
    data <- melt(comp.treat.data, id.vars = NULL, variable.name = "beta", value.name = "estimate")
    
    ggplot(data, aes(x = estimate, color = beta, fill = beta)) + 
      geom_density(na.rm = TRUE, size = 0.5, alpha = 0) +
      geom_histogram(aes(y = ..density..), alpha = 0.5, na.rm = TRUE, ...) + 
      theme_bw() + 
      facet_wrap(. ~ beta) + 
      # xlim(xlim_value[1], xlim_value[2]) + 
      labs(x = xlab, y = "Density", 
           fill = "Treatment Comparison", color = "Treatment Comparison") 
  }
  
}



trt_eff_plot <- function(result.list, level = 0.95, ...){
  
  trt.eff <- NULL
  for(i in 1:length(result.list)){
    
    result <- result.list[[i]]
    samples <- do.call(rbind, result$samples)
    # beta_variable <- exp(samples[,grep("beta", colnames(samples))])
    beta_variable <- samples[, grep("beta", colnames(samples))]
    data <- as.data.frame(beta_variable)
    
    comp.treat.name <- t(utils::combn(colnames(data), 2))
    for (j in 1:nrow(comp.treat.name)) {
      tmp.compname <- gsub("beta_", "", comp.treat.name[j, ])
      tmp.compname <- paste0(tmp.compname[2], " - ", tmp.compname[1])
      
      trt.eff <- rbind(trt.eff, 
                       quantile(data[, comp.treat.name[j, 2]] -  data[, comp.treat.name[j, 1]], 
                                probs = c((1 -level)/2, 0.5, 1 - (1 -level)/2)))
      
      trt.eff <- data.frame(trt.eff)
      rownames(trt.eff)[nrow(trt.eff)] <- paste0(names(result.list)[i], ": ", tmp.compname)
    }
  
  }
  
  trt.eff.range <- range(trt.eff, 0)
  trt.eff <- trt.eff %>%
    mutate(`Treatment Comparison` = rownames(trt.eff))
  colnames(trt.eff)[1:3] <- c("lower", "median", "upper")
  
  ggplot(trt.eff, aes(y = median, x = `Treatment Comparison`, color  = `Treatment Comparison`)) + 
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), ...) +
    # scale_y_log10(breaks = ticks, labels = ticks) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray") +
    coord_flip() +
    labs(x = "Treatment Comparison", y = "Treatment effect") +
    theme_bw() + 
    theme(legend.position = "none")
}



probability_barplot <- function(result.list){
  
  probability <- NULL
  for(i in 1:length(result.list)) {
    
    result <- result.list[[i]]
    samples <- do.call(rbind, result$samples)
    beta_variable <- samples[, grep("beta", colnames(samples))]
    data <- as.data.frame(beta_variable)
    
    comp.treat.name <- t(utils::combn(colnames(data), 2))
    for (j in 1:nrow(comp.treat.name)) {
      tmp.compname <- gsub("beta_", "", comp.treat.name[j, ])
      tmp.compname <- paste0(tmp.compname[1], " vs. ", tmp.compname[2])
      
      probability <- rbind(probability, 
                           c(mean((data[, comp.treat.name[j, 1]] -  data[, comp.treat.name[j, 2]]) >= 0), 
                             mean((data[, comp.treat.name[j, 2]] -  data[, comp.treat.name[j, 1]]) > 0)))
      
      probability <- data.frame(probability)
      rownames(probability)[nrow(probability)] <- paste0(names(result.list)[i], ": ", tmp.compname)
    }
  }
  
  probability <- probability %>%
    mutate(result.name = rownames(probability))
  probability <- probability %>%
    tidyr::separate(result.name, c("result", "trt.comp"), sep = ": ", extra = "merge") %>%
    tidyr::separate(trt.comp, c("trt.1", "trt.2"), sep = " vs. ", remove = F)
  
  probability <- probability %>%
    gather(trt, probability, X1:X2) %>%
    mutate(trt = ifelse(trt == "X1", trt.1, trt.2)) %>%
    dplyr::select(-c(trt.1, trt.2))
  
  ggplot(probability, aes(fill = factor(trt), y = probability, x = result)) + 
    geom_bar(stat="identity", position="fill") + 
    scale_y_continuous(labels = percent_format()) + 
    labs(x = "Result", y = "Percentages", fill = "Treatment") + 
    facet_wrap(. ~ trt.comp) + 
    coord_flip()+  
    theme_bw()  
}
