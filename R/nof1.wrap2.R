read_input_data2 <- function(data, metadata){ 

  Y <- unlist(data$response$afib_episode_yn)

  Treatment <- data$treatment
  length_each <- sapply(data$response$afib_episode_yn, length)
  Treat <- rep(Treatment, time = length_each)

  # Treat[Treat == "control"] = "baseline"
  # Treat[Treat == "trigger"] = "A"
  # for (i in Treatment){
  #   if (sum(Treat == i) == sum(is.na(Y[Treat == i]))){
  #     Y     <- Y[Treat != i]
  #     Treat <- Treat[Treat != i]
  #   }
  # }

  list(Treat = Treat, Y = Y)
}

# find_raw_mean2 <- function(Y, Treat, baseline, response){
# 
#   raw_mean <- c(mean(Y[Treat == baseline], na.rm = TRUE), mean(Y[Treat == "A"], na.rm = TRUE))
#   raw_mean[is.nan(raw_mean)] <- NA
#   raw_mean
# }

# check_enough_data2 <- function(Treatment, x){
#   length(table(Treatment[!is.na(x)])) == 2
# }

summarize_nof1_afib <- function(nof1, result, treat.name){

  with(c(nof1, result),{

    samples <- do.call(rbind, samples)
    # samples <- do.call(rbind, result$samples)
    # Y        <- nof1$Y
    # Treat    <- nof1$Treat
    # response <- nof1$response
    raw_mean <- find_raw_mean(Y, Treat, treat.name)
    rounded_raw_mean <- round_number(raw_mean, response)
    raw_mean <- list(control = rounded_raw_mean[1], trigger = rounded_raw_mean[2])

    # An odds ratio of 1 indicates that the condition or event under study is equally likely to occur in both groups. 
    # An odds ratio greater than 1 indicates that the condition or event is more likely to occur in the first group. 
    # And an odds ratio less than 1 indicates that the condition or event is less likely to occur in the first group.

    coef <- samples[, colnames(samples) %in% paste("beta_", treat.name, sep = "")]
    
    for (i in treat.name){
      col.treat.name <- paste("beta_", i, sep = "")
      if (col.treat.name %in% colnames(coef)){
        assign(paste("coef_", col.treat.name, sep = ""), coef[, col.treat.name, drop = F])
      }
    }
    
    # hist(coef_beta_control)
    # hist(coef_beta_trigger - coef_beta_control)
    
    # base <- inv_logit(coef_alpha)
    # trigger <- inv_logit(coef_alpha + coef_beta_A)
    
    greater_than_1 <- round(mean(coef_beta_trigger > coef_beta_control) *100)
    # greater_than_1 <- round(mean((trigger/base) > 1, na.rm = TRUE)*100)
    greater_than_1 <- change(greater_than_1)
    
    return(list(raw_mean = raw_mean, prob_afib_more_likely_with_trigger = greater_than_1))
  })
}

#' For PCORI purposes
#'
#' @export

wrap2 <- function(data, metadata){

  read_data <- tryCatch({
    read_dummy <- read_input_data2(data, metadata)
    read_dummy
  }, error = function(error){
    return(paste("input read error: ", error))
  })

  print(read_data)
  treat.name <- c("control", "trigger")

  afib <- tryCatch({
    data_afib <- read_data
    nof1_afib <- with(data_afib, {
      # Y <- data_afib$Y
      # Treat <- data_afib$Treat
      nof1.data(Y, Treat, response = "binomial")
    })
    # nof1 <- nof1_afib
    result_afib <- nof1.run(nof1_afib)
    # result <- result_afib
    summarize_nof1_afib(nof1_afib, result_afib, treat.name)
  }, error = function(error){
    return(paste("afib run error: ", error))
  })

  metadata <- list(
                   successful_input_reading = check_success(read_data),
                   successful_run_afib = check_success(afib),
                   enough_afib = check_enough_data(read_data$Treat, 
                                                   read_data$Y,
                                                   treat.name),
                   user_id = metadata$user_id,
                   trigger = metadata$trigger,
                   design = metadata$design,
                   timestamp_sammy_completed = Sys.time(),
                   sammy_version_id = 1,
                   sammy_version_date = "8/15/2017",
                   sammy_version_note = "")

  final <- list(metadata = metadata, afib = afib)
  return(final)
}
