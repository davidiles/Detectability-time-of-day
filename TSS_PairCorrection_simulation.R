library(mgcv)
library(reshape2)
library(tidyverse)

rm(list=ls())

# Simulate an effect of time-since-sunrise (TSS) that differs between each of two sexes
#  - then evaluate the ability of various approaches to correct for those effects and estimate true (or even just relative) density

#***************************************************************
#* Attempt #1 - using intercept-only removal model fit to all survey times
#***************************************************************

simulation_results <- data.frame(sim_rep = 1:1000,
                                 poptotal_true = NA,
                                 poptotal_est = NA)

for (rep in simulation_results$sim_rep){
  
  # --------------------------------
  # Simulation parameters
  # --------------------------------
  
  n_surveys <- 2500
  
  male_density <- 5 # expected number of birds per survey (assume no distance effects, same density at every location)
  female_ratio <- 1 # number of females per male
  
  survey_duration <- 5  # minutes per survey
  tint <- 1:survey_duration # 1-minute time intervals
  
  # Dataframe to hold survey data
  survey_data <- data.frame(survey_number = 1:n_surveys,
                            n_males = rpois(n_surveys,male_density),
                            n_females = rpois(n_surveys,male_density*female_ratio),
                            time_since_sunrise = runif(n_surveys,-2,6))
  
  # --------------------------------
  # Relationships between time-since-sunrise and availability
  # --------------------------------
  
  TSS_df <- data.frame(time_since_sunrise = seq(-2,10,length.out = 500))
  
  # Time curve for males
  TSS_df$phi_male <- rep(0.2,nrow(TSS_df))
  time_since_sunrise <- seq(-4,10,by = 1) 
  
  # use this to create a somewhat realistic curve and then fit a spline through it to smooth it out for simulation
  phi_male <- c(0,0.1,0.1,0.1,0.4,0.5,0.3,0.3,0.2,0.2,0.2,0.2,0.15,0.1,0.1)
  #phi_male <- c(1,1,0.8,0.5,0.2,0.1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
  #phi_male <- cumsum(rep(0.05,length(time_since_sunrise)))
  
  gam_male <- gam(phi_male~s(time_since_sunrise, k = 11))
  TSS_df$phi_male <- predict(gam_male, newdata = TSS_df)
  
  # Females assumed to have constant detection
  TSS_df$phi_female <- rep(0,nrow(TSS_df))
  
  #plot(phi_male~time_since_sunrise, data = TSS_df, type ="l", col = "blue", lwd = 2, ylab = "cue rate per minute", xlab = "time since sunrise", main = "cue rate across the morning")
  #lines(phi_female~time_since_sunrise, data = TSS_df, col = "red", lwd = 2)
  
  # --------------------------------
  # Simulate survey data; attempt to recreate the curve and estimate total density
  # --------------------------------
  
  survey_data$phi_male <- predict(gam_male, newdata = survey_data)
  survey_data$phi_female <- mean(TSS_df$phi_female)
  
  # ------------------------------------
  # Simulate bird cues, based on phi_true
  # ------------------------------------
  
  # Array to store counts within each time interval for each survey
  Y_array <- array(0,dim=c(n_surveys,length(tint)))
  
  for (i in 1:n_surveys){
    
    # ------------------------------------
    # Simulate cues produced by each bird
    # ------------------------------------
    
    n_males <- survey_data$n_males[i]
    n_females <- survey_data$n_females[i]
    
    # Simulate first cues produced by each male bird during this particular survey
    cue_males <- rexp(n_males,survey_data$phi_male[i])
    detected_males <- cue_males[cue_males<=survey_duration]
    if (sum(!is.na(detected_males))==0) detected_males <- c()
    
    # Simulate first cues produced by each female bird during this particular survey
    cue_females <- rexp(n_females,survey_data$phi_female[i])
    detected_females <- cue_females[cue_females<=survey_duration]
    if (sum(!is.na(detected_females))==0) detected_females <- c()
    
    n_birds_detected <- length(detected_males) + length(detected_females)
    if (n_birds_detected == 0) next
    
    # ------------------------------------
    # Transcribed data for this survey
    # ------------------------------------
    
    idat <- data.frame(bird_number = 1:n_birds_detected,
                       time = c(detected_males,detected_females))
    
    # Separate into distance and time bins
    idat$tint <- cut(idat$time,c(0,tint))
    
    # Store in overall array
    Y_array[i,] <- table(idat[,c("tint")])
    
  }
  
  # ------------------------------------
  # Analysis: removal model
  # ------------------------------------
  
  D_array <- matrix(rep(tint,n_surveys),n_surveys, byrow = T)
  
  #fit_p <- detect::cmulti(Y_array | D_array ~ survey_data$time_since_sunrise + I(survey_data$time_since_sunrise^2), type = "rem")
  #TSS_df$phi_pred <- exp(fit_p$coefficients[1] + fit_p$coefficients[2]*TSS_df$time_since_sunrise + fit_p$coefficients[3]*TSS_df$time_since_sunrise^2)
  #survey_data$phi_pred <- exp(fit_p$coefficients[1] + fit_p$coefficients[2]*survey_data$time_since_sunrise + fit_p$coefficients[3]*survey_data$time_since_sunrise^2)
  
  fit_p <- detect::cmulti(Y_array | D_array ~ 1, type = "rem")
  phi_pred <- exp(fit_p$coefficients[1])
  offset <- 1-exp(-survey_duration*phi_pred) # Proportion of birds that were missed
  
  TSS_df$phi_pred <- phi_pred
  TSS_df$offset <- offset
  
  survey_data$phi_pred <- phi_pred
  survey_data$offset <- offset
  
  # Check out the estimate
  testplot <- ggplot(data = TSS_df)+
    geom_line(aes(x = time_since_sunrise, y = phi_male), col = "blue", linewidth=2)+
    geom_line(aes(x = time_since_sunrise, y = phi_female), col = "red", linewidth=2)+
    geom_line(aes(x = time_since_sunrise, y = phi_pred), col = "black", linewidth=2, linetype = 2)+
    geom_vline(xintercept = mean(survey_data$time_since_sunrise))+
    theme_bw()+
    ylab("Cue rate")+
    xlab("Hours Since Sunrise")
  #print(testplot)

  # ------------------------------------
  # Attempt to estimate a temporal correction factor and estimate overall population total using QPAD and GAM effect of time of day
  # ------------------------------------
  
  # Observed counts for each survey
  survey_data$count <- rowSums(Y_array)
  
  gamfit <- gam(count ~ s(time_since_sunrise, k = -1),family = "poisson", data = survey_data)
  
  TSS_df$gam_pred <- predict(gamfit, newdata = TSS_df, type = "response")
  
  testplot2 <- ggplot(data = TSS_df)+
    geom_line(aes(x = time_since_sunrise, y = gam_pred), col = "black", linewidth=2, linetype = 2)+
    theme_bw()+
    ylab("Expected count")+
    xlab("Hours Since Sunrise")
  #print(testplot2)
  
  # Prediction of what you would have seen at mean survey timing
  survey_data$gam_pred <- predict(gamfit, 
                                  newdata = data.frame(time_since_sunrise = rep(mean(survey_data$time_since_sunrise),nrow(survey_data))), 
                                 type = "response") / offset
  
  # Estimate of total population size
  simulation_results$poptotal_est[rep] <- sum(survey_data$gam_pred)
  simulation_results$poptotal_true[rep] <- sum(survey_data$n_males)
  simulation_results$bias <- simulation_results$poptotal_est - simulation_results$poptotal_true
  simulation_results$percent_bias <- 100*(simulation_results$poptotal_est - simulation_results$poptotal_true)/simulation_results$poptotal_true
  
  print(rep)
  
  hist(simulation_results$percent_bias, breaks = seq(-100,100,1), main = "Percent bias in abundance estimates across repeated simulations", xlab = "Percent bias")
  abline(v = mean(simulation_results$percent_bias,na.rm = TRUE), col = "blue", lwd = 2)
  
}


hist(simulation_results$bias)
mean(simulation_results$bias, na.rm = TRUE)

hist(simulation_results$percent_bias)
mean(simulation_results$percent_bias, na.rm = TRUE)


#***************************************************************
#* Attempt #2 - using intercept-only removal model within 1 hour of sunrise
#* This removes the bias caused by fitting an intercept-only removal model across all times (where detectability varies a lot)
#***************************************************************

simulation_results <- data.frame(sim_rep = 1:1000,
                                 poptotal_true = NA,
                                 poptotal_est = NA)

for (rep in simulation_results$sim_rep){
  
  # --------------------------------
  # Simulation parameters
  # --------------------------------
  
  n_surveys <- 2500
  
  male_density <- 5 # expected number of birds per survey (assume no distance effects, same density at every location)
  female_ratio <- 1 # number of females per male
  
  survey_duration <- 5  # minutes per survey
  tint <- 1:survey_duration # 1-minute time intervals
  
  # Dataframe to hold survey data
  survey_data <- data.frame(survey_number = 1:n_surveys,
                            n_males = rpois(n_surveys,male_density),
                            n_females = rpois(n_surveys,male_density*female_ratio),
                            time_since_sunrise = runif(n_surveys,-2,6))
  
  # --------------------------------
  # Relationships between time-since-sunrise and availability
  # --------------------------------
  
  TSS_df <- data.frame(time_since_sunrise = seq(-2,10,length.out = 500))
  
  # Time curve for males
  TSS_df$phi_male <- rep(0.2,nrow(TSS_df))
  time_since_sunrise <- seq(-4,10,by = 1) # use this to create a somewhat realistic curve and then fit a spline through it to smooth it out for simulation
  
  # use this to create a somewhat realistic curve and then fit a spline through it to smooth it out for simulation
  #phi_male <- c(0,0.1,0.1,0.1,0.4,0.5,0.3,0.3,0.2,0.2,0.2,0.2,0.15,0.1,0.1)
  phi_male <- c(1,1,0.8,0.5,0.2,0.1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
  #phi_male <- cumsum(rep(0.05,length(time_since_sunrise)))
  
  gam_male <- gam(phi_male~s(time_since_sunrise, k = 11))
  TSS_df$phi_male <- predict(gam_male, newdata = TSS_df)
  
  # Females assumed to have constant detection
  TSS_df$phi_female <- rep(0,nrow(TSS_df))
  
  # plot(phi_male~time_since_sunrise, data = TSS_df, ylim = c(0,1), type ="l", col = "blue", lwd = 2, ylab = "cue rate per minute", xlab = "time since sunrise", main = "cue rate across the morning")
  # lines(phi_female~time_since_sunrise, data = TSS_df, col = "red", lwd = 2)
  
  # --------------------------------
  # Simulate survey data; attempt to recreate the curve and estimate total density
  # --------------------------------
  
  survey_data$phi_male <- predict(gam_male, newdata = survey_data)
  survey_data$phi_female <- mean(TSS_df$phi_female)
  
  # ------------------------------------
  # Simulate bird cues, based on phi_true
  # ------------------------------------
  
  # Array to store counts within each time interval for each survey
  Y_array <- array(0,dim=c(n_surveys,length(tint)))
  
  for (i in 1:n_surveys){
    
    # ------------------------------------
    # Simulate cues produced by each bird
    # ------------------------------------
    
    n_males <- survey_data$n_males[i]
    n_females <- survey_data$n_females[i]
    
    # Simulate first cues produced by each male bird during this particular survey
    cue_males <- rexp(n_males,survey_data$phi_male[i])
    detected_males <- cue_males[cue_males<=survey_duration]
    if (sum(!is.na(detected_males))==0) detected_males <- c()
    
    # Simulate first cues produced by each female bird during this particular survey
    cue_females <- rexp(n_females,survey_data$phi_female[i])
    detected_females <- cue_females[cue_females<=survey_duration]
    if (sum(!is.na(detected_females))==0) detected_females <- c()
    
    n_birds_detected <- length(detected_males) + length(detected_females)
    if (n_birds_detected == 0) next
    
    # ------------------------------------
    # Transcribed data for this survey
    # ------------------------------------
    
    idat <- data.frame(bird_number = 1:n_birds_detected,
                       time = c(detected_males,detected_females))
    
    # Separate into distance and time bins
    idat$tint <- cut(idat$time,c(0,tint))
    
    # Store in overall array
    Y_array[i,] <- table(idat[,c("tint")])
    
  }
  
  # ------------------------------------
  # Analysis: removal model
  # ------------------------------------
  
  surveys_for_rem_model <- which(survey_data$time_since_sunrise > -0.5 & survey_data$time_since_sunrise < 0.5)
  
  D_array <- matrix(rep(tint,n_surveys),n_surveys, byrow = T)
  
  Y_array2 <- Y_array[surveys_for_rem_model,]
  D_array2 <- D_array[surveys_for_rem_model,]
  
  fit_p <- detect::cmulti(Y_array2 | D_array2 ~ 1, type = "rem")
  phi_pred <- exp(fit_p$coefficients[1])
  offset <- 1-exp(-survey_duration*phi_pred) # Proportion of birds that were missed
  
  TSS_df$phi_pred <- phi_pred
  TSS_df$offset <- offset
  
  survey_data$phi_pred <- phi_pred
  survey_data$offset <- offset
  
  # Check out the estimate
  testplot <- ggplot(data = TSS_df)+
    geom_line(aes(x = time_since_sunrise, y = phi_male), col = "blue", linewidth=2)+
    geom_line(aes(x = time_since_sunrise, y = phi_female), col = "red", linewidth=2)+
    geom_line(aes(x = time_since_sunrise, y = phi_pred), col = "black", linewidth=2, linetype = 2)+
    theme_bw()+
    ylab("Cue rate")+
    xlab("Hours Since Sunrise")
  #print(testplot)
  
  # ------------------------------------
  # Attempt to estimate a temporal correction factor and estimate overall population total using QPAD and GAM effect of time of day
  # ------------------------------------
  
  # Observed counts for each survey
  survey_data$count <- rowSums(Y_array)
  
  gamfit <- gam(count ~ s(time_since_sunrise, k = -1),family = "poisson", data = survey_data)
  
  TSS_df$gam_pred <- predict(gamfit, newdata = TSS_df, type = "response")
  
  testplot2 <- ggplot(data = TSS_df)+
    geom_line(aes(x = time_since_sunrise, y = gam_pred), col = "black", linewidth=2, linetype = 2)+
    theme_bw()+
    ylab("Expected count")+
    xlab("Hours Since Sunrise")
  
  #print(testplot2)
  
  # Prediction of what you would have seen at sunrise
  survey_data$gam_pred <- predict(gamfit, 
                                  newdata = data.frame(time_since_sunrise = rep(0,nrow(survey_data))), 
                                  type = "response") / offset
  
  # Estimate of total population size
  simulation_results$poptotal_est[rep] <- sum(survey_data$gam_pred)
  simulation_results$poptotal_true[rep] <- sum(survey_data$n_males)
  simulation_results$bias <- simulation_results$poptotal_est - simulation_results$poptotal_true
  simulation_results$percent_bias <- 100*(simulation_results$poptotal_est - simulation_results$poptotal_true)/simulation_results$poptotal_true
  
  print(rep)
  
  hist(simulation_results$percent_bias, breaks = seq(-100,100,1), main = "Percent bias in abundance estimates across repeated simulations", xlab = "Percent bias")
  abline(v = mean(simulation_results$percent_bias,na.rm = TRUE), col = "blue", lwd = 2)
  
}


hist(simulation_results$bias)
mean(simulation_results$bias, na.rm = TRUE)

hist(simulation_results$percent_bias)
mean(simulation_results$percent_bias, na.rm = TRUE)


#***************************************************************
#* Attempt #3 - using intercept-only removal model within 1 hour of sunrise
#* This removes the bias caused by fitting an intercept-only removal model across all times (where detectability varies a lot)
#* ALSO INCORPORATE FEMALES INTO THE ANALYSIS
#***************************************************************

simulation_results <- data.frame(sim_rep = 1:1000,
                                 poptotal_true = NA,
                                 poptotal_est = NA)

for (rep in simulation_results$sim_rep){
  
  # --------------------------------
  # Simulation parameters
  # --------------------------------
  
  n_surveys <- 2500
  
  male_density <- 5 # expected number of birds per survey (assume no distance effects, same density at every location)
  female_ratio <- 1 # number of females per male
  
  survey_duration <- 5  # minutes per survey
  tint <- 1:survey_duration # 1-minute time intervals
  
  # Dataframe to hold survey data
  survey_data <- data.frame(survey_number = 1:n_surveys,
                            n_males = rpois(n_surveys,male_density),
                            n_females = rpois(n_surveys,male_density*female_ratio),
                            time_since_sunrise = runif(n_surveys,-2,6))
  
  # --------------------------------
  # Relationships between time-since-sunrise and availability
  # --------------------------------
  
  TSS_df <- data.frame(time_since_sunrise = seq(-2,10,length.out = 500))
  
  # Time curve for males
  TSS_df$phi_male <- rep(0.2,nrow(TSS_df))
  time_since_sunrise <- seq(-4,10,by = 1) # use this to create a somewhat realistic curve and then fit a spline through it to smooth it out for simulation
  
  # use this to create a somewhat realistic curve and then fit a spline through it to smooth it out for simulation
  #phi_male <- c(0,0.1,0.1,0.1,0.4,0.5,0.3,0.3,0.2,0.2,0.2,0.2,0.15,0.1,0.1)
  phi_male <- c(1,1,0.8,0.5,0.2,0.1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
  #phi_male <- cumsum(rep(0.05,length(time_since_sunrise)))
  
  gam_male <- gam(phi_male~s(time_since_sunrise, k = 11))
  TSS_df$phi_male <- predict(gam_male, newdata = TSS_df)
  
  # Females assumed to have constant detection
  TSS_df$phi_female <- rep(0.1,nrow(TSS_df))
  
  ##plot(phi_male~time_since_sunrise, data = TSS_df, ylim = c(0,1), type ="l", col = "blue", lwd = 2, ylab = "cue rate per minute", xlab = "time since sunrise", main = "cue rate across the morning")
  #lines(phi_female~time_since_sunrise, data = TSS_df, col = "red", lwd = 2)
  
  # --------------------------------
  # Simulate survey data; attempt to recreate the curve and estimate total density
  # --------------------------------
  
  survey_data$phi_male <- predict(gam_male, newdata = survey_data)
  survey_data$phi_female <- mean(TSS_df$phi_female)
  
  # ------------------------------------
  # Simulate bird cues, based on phi_true
  # ------------------------------------
  
  # Array to store counts within each time interval for each survey
  Y_array <- array(0,dim=c(n_surveys,length(tint)))
  n_sexes_detected <- array(0,dim=c(n_surveys,2))
  
  for (i in 1:n_surveys){
    
    # ------------------------------------
    # Simulate cues produced by each bird
    # ------------------------------------
    
    n_males <- survey_data$n_males[i]
    n_females <- survey_data$n_females[i]
    
    # Simulate first cues produced by each male bird during this particular survey
    cue_males <- rexp(n_males,survey_data$phi_male[i])
    detected_males <- cue_males[cue_males<=survey_duration]
    if (sum(!is.na(detected_males))==0) detected_males <- c()
    
    # Simulate first cues produced by each female bird during this particular survey
    cue_females <- rexp(n_females,survey_data$phi_female[i])
    detected_females <- cue_females[cue_females<=survey_duration]
    if (sum(!is.na(detected_females))==0) detected_females <- c()
    
    n_birds_detected <- length(detected_males) + length(detected_females)
    if (n_birds_detected == 0) next
    
    n_sexes_detected[i,1] <- length(detected_males)
    n_sexes_detected[i,2] <- length(detected_females)
    
    # ------------------------------------
    # Transcribed data for this survey
    # ------------------------------------
    
    idat <- data.frame(bird_number = 1:n_birds_detected,
                       time = c(detected_males,detected_females))
    
    # Separate into distance and time bins
    idat$tint <- cut(idat$time,c(0,tint))
    
    # Store in overall array
    Y_array[i,] <- table(idat[,c("tint")])
    
  }
  
  # ------------------------------------
  # Analysis: removal model
  # ------------------------------------
  
  surveys_for_rem_model <- which(survey_data$time_since_sunrise > -0.5 & survey_data$time_since_sunrise < 0.5)
  
  D_array <- matrix(rep(tint,n_surveys),n_surveys, byrow = T)
  
  Y_array2 <- Y_array[surveys_for_rem_model,]
  D_array2 <- D_array[surveys_for_rem_model,]
  
  fit_p <- detect::cmulti(Y_array2 | D_array2 ~ 1, type = "rem")
  phi_pred <- exp(fit_p$coefficients[1])
  offset <- 1-exp(-survey_duration*phi_pred) # Proportion of birds that were missed
  
  TSS_df$phi_pred <- phi_pred
  TSS_df$offset <- offset
  
  survey_data$phi_pred <- phi_pred
  survey_data$offset <- offset
  
  # Check out the estimate
  testplot <- ggplot(data = TSS_df)+
    geom_line(aes(x = time_since_sunrise, y = phi_male), col = "blue", linewidth=2)+
    geom_line(aes(x = time_since_sunrise, y = phi_female), col = "red", linewidth=2)+
    geom_line(aes(x = time_since_sunrise, y = phi_pred), col = "black", linewidth=2, linetype = 2)+
    theme_bw()+
    ylab("Cue rate")+
    xlab("Hours Since Sunrise")
  #print(testplot)
  
  # ------------------------------------
  # Attempt to estimate a temporal correction factor and estimate overall population total using QPAD and GAM effect of time of day
  # ------------------------------------
  
  # Observed counts for each survey
  survey_data$count <- rowSums(Y_array)
  
  gamfit <- gam(count ~ s(time_since_sunrise, k = -1),family = "poisson", data = survey_data)
  
  TSS_df$gam_pred <- predict(gamfit, newdata = TSS_df, type = "response")
  
  testplot2 <- ggplot(data = TSS_df)+
    geom_line(aes(x = time_since_sunrise, y = gam_pred), col = "black", linewidth=2, linetype = 2)+
    theme_bw()+
    ylab("Expected count")+
    xlab("Hours Since Sunrise")
  #print(testplot2)
  
  # Overlay mean estimates of number of males and females (and total) birds detected across time bins
  tseq <- seq(min(survey_data$time_since_sunrise), max(survey_data$time_since_sunrise), length.out = 21)
  
  sex_ratio_df <- data.frame()
  for (t in 2:length(tseq)){
    sdat <- subset(survey_data, time_since_sunrise > tseq[t-1] & time_since_sunrise < tseq[t])
    sex_ratio_df <- rbind(sex_ratio_df, data.frame(time_since_sunrise = mean(c(tseq[t-1],tseq[t])),
                                                   n_males = mean(n_sexes_detected[sdat$survey_number,1]),
                                                   n_females = mean(n_sexes_detected[sdat$survey_number,2])))
  }
  
  testplot2 <- testplot2 +
    geom_line(data = sex_ratio_df, aes(x = time_since_sunrise, y = n_males), col = "blue")+
    geom_line(data = sex_ratio_df, aes(x = time_since_sunrise, y = n_females), col = "red")+
    geom_line(data = sex_ratio_df, aes(x = time_since_sunrise, y = n_females+n_males), col = "gray", linewidth = 2)
  #testplot2
  
  # Assume we can estimate pair correction reliably (need to know the relative number of females detected during point counts at dawn)
  pair_correction <- subset(sex_ratio_df,time_since_sunrise > -0.5 & time_since_sunrise < 0.5)
  proportion_males_in_dawn_count <- sum(pair_correction$n_males)/sum(pair_correction$n_females + pair_correction$n_males)
  proportion_males_in_dawn_count
  
  
  # Prediction of what you would have seen at sunrise
  survey_data$gam_pred <- predict(gamfit, 
                                  newdata = data.frame(time_since_sunrise = rep(0,nrow(survey_data))), 
                                  type = "response") / offset * proportion_males_in_dawn_count
  
  # Estimate of total population size
  simulation_results$poptotal_est[rep] <- sum(survey_data$gam_pred)
  simulation_results$poptotal_true[rep] <- sum(survey_data$n_males)
  simulation_results$bias <- simulation_results$poptotal_est - simulation_results$poptotal_true
  simulation_results$percent_bias <- 100*(simulation_results$poptotal_est - simulation_results$poptotal_true)/simulation_results$poptotal_true
  
  print(rep)
  
  hist(simulation_results$percent_bias, breaks = seq(-100,100,1), main = "Percent bias in abundance estimates across repeated simulations", xlab = "Percent bias")
  abline(v = mean(simulation_results$percent_bias,na.rm = TRUE), col = "blue", lwd = 2)
  
}


hist(simulation_results$bias)
mean(simulation_results$bias, na.rm = TRUE)

hist(simulation_results$percent_bias)
mean(simulation_results$percent_bias, na.rm = TRUE)
