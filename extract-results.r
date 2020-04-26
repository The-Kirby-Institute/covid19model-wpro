## Functions for exploring and plotting results from Imperial model
#
# Richard T. Gray
# 
# This script when sourced loads various functions for extracting and
# producing results. Based on the code used in plot-3-panel.r and 
# plot-forecast.r. Attempt to make it more succinct and based on my 
# preferences.

library(ggplot2)
library(tidyr)
library(dplyr)
library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(EnvStats)
library(matrixStats)
library(scales)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)
library(bayesplot)

source("utils/geom-stepribbon.r")

# Fitting results ------------------------
StanResults <- function(countries,JOBID,out,resultsDir) {
  filename <- paste0('base-',JOBID)
  mu = (as.matrix(out$mu))
  colnames(mu) = countries
  g = (mcmc_intervals(mu,prob = .9))
  ggsave(paste0(resultsDir, filename, "-mu.png"),g,width=4,height=6)
  # ggsave(sprintf("results/%s-mu.png",filename),g,width=4,height=6)
  tmp = lapply(1:length(countries), function(i) (out$Rt_adj[,stan_data$N[i],i]))
  Rt_adj = do.call(cbind,tmp)
  colnames(Rt_adj) = countries
  g = (mcmc_intervals(Rt_adj,prob = .9))
  ggsave(paste0(resultsDir, filename, "-final-rt.png"),g,width=4,height=6)
  # ggsave(sprintf("results/%s-final-rt.png",filename),g,width=4,height=6)
  return(list("mu" = mu, "Rt_adj" = Rt_adj))
}

## Impact no NPIs -predictions if R0 stays at original value
resultsR0 <- function(stanData, stanOut) {
  
  infections <- array(rep(0, 100*stanData$N2*stanData$M), 
    c(100, stanData$N2, stanData$M)) 
  cummInfects <- infections
  y <- stanOut$prediction[,1,]
  
  for (ii in 1:stanData$M) {
    for (jj in 2:stanData$N0){
      cummInfects[,jj,ii] <- cummInfects[,jj-1,ii] + y[,ii]; 
      infections[,jj,ii] <- y[,ii]
    }
    
    for (jj in (stanData$N0+1):stanData$N2) {
      convolution0 <- rep(0, 100)
      for(kk in 1:(jj-1)) {
        convolution0 <- convolution0 + infections[,kk,ii] * 
          stanData$SI[jj-kk]; 
      }
      cummInfects[,jj,ii] = cummInfects[,jj-1,ii] + infections[,jj-1,ii]; 
      infections[,jj,ii] =  ((stanData$pop[ii]-cummInfects[,jj,ii]) / stanData$pop[ii]) * 
        stanOut$mu[,ii] * convolution0;
    }
  }
  
  return(infections)
  
}



# Sort out covariates ----------------------------------------------------
TidyCovariates <- function(ncountries, covariates) {  
  names_covariates = c('Schools + Universities','Self-isolating if ill', 
    'Public events', 'Lockdown', 'Social distancing encouraged')
  covariates <- covariates %>%
    filter((Type %in% names_covariates))
  covariates <- covariates[,c(1,2,4)]
  covariates <- spread(covariates, Type, Date.effective)
  names(covariates) <- c('Country','lockdown', 'public_events', 
    'schools_universities','self_isolating_if_ill', 
    'social_distancing_encouraged')
  covariates <- covariates[c('Country','schools_universities', 
    'self_isolating_if_ill', 'public_events', 'lockdown', 
    'social_distancing_encouraged')]
  covariates$schools_universities <- as.Date(covariates$schools_universities, 
    format = "%d.%m.%Y")
  covariates$lockdown <- as.Date(covariates$lockdown, format = "%d.%m.%Y")
  covariates$public_events <- as.Date(covariates$public_events, 
    format = "%d.%m.%Y")
  covariates$self_isolating_if_ill <- as.Date(covariates$self_isolating_if_ill, 
    format = "%d.%m.%Y")
  covariates$social_distancing_encouraged <- as.Date(covariates$social_distancing_encouraged, format = "%d.%m.%Y")
  
  return(covariates)
  
}

CountryCovariates <- function(country, covariates,rt_ui) {
  
  covariates_country <- covariates[which(covariates$Country == country), 2:6] 
  covariates_country_long <- gather(covariates_country[], key = "key", 
    value = "value")
  covariates_country_long$x <- rep(NULL, length(covariates_country_long$key))
  un_dates <- unique(covariates_country_long$value)
  
  for (k in 1:length(un_dates)){
    idxs <- which(covariates_country_long$value == un_dates[k])
    max_val <- round(max(rt_ui)) + 0.3
    for (j in idxs){
      covariates_country_long$x[j] <- max_val
      max_val <- max_val - 0.3
    }
  }
  
  covariates_country_long$value <- as_date(covariates_country_long$value) 
  covariates_country_long$country <- rep(country, 
    length(covariates_country_long$value))
  
  return(covariates_country_long)
}    

# Country estimates  ------------------------------------------------------
CountryOutputs <- function(index,country,dates,reported_cases,
  deaths_by_country,prediction,estimated.deaths,rt) {
  
  N <- length(dates)
  predicted_cases <- colMeans(prediction[,1:N,index])
  predicted_cases_li <- colQuantiles(prediction[,1:N,index], probs=.025)
  predicted_cases_ui <- colQuantiles(prediction[,1:N,index], probs=.975)
  predicted_cases_li2 <- colQuantiles(prediction[,1:N,index], probs=.25)
  predicted_cases_ui2 <- colQuantiles(prediction[,1:N,index], probs=.75)
  
  estimated_deaths <- colMeans(estimated.deaths[,1:N,index])
  estimated_deaths_li <- colQuantiles(estimated.deaths[,1:N,index], probs=.025)
  estimated_deaths_ui <- colQuantiles(estimated.deaths[,1:N,index], probs=.975)
  estimated_deaths_li2 <- colQuantiles(estimated.deaths[,1:N,index], probs=.25)
  estimated_deaths_ui2 <- colQuantiles(estimated.deaths[,1:N,index], probs=.75)
  
  rt <- colMeans(out$Rt[,1:N,index])
  rt_li <- colQuantiles(out$Rt[,1:N,index],probs=.025)
  rt_ui <- colQuantiles(out$Rt[,1:N,index],probs=.975)
  rt_li2 <- colQuantiles(out$Rt[,1:N,index],probs=.25)
  rt_ui2 <- colQuantiles(out$Rt[,1:N,index],probs=.75)
  
  data_country <- data.frame("time" = as_date(as.character(dates)),
    "country" = rep(country, length(dates)),
    "reported_cases" = reported_cases[[index]], 
    "reported_cases_c" = cumsum(reported_cases[[index]]), 
    "predicted_cases_c" = cumsum(predicted_cases),
    "predicted_min_c" = cumsum(predicted_cases_li),
    "predicted_max_c" = cumsum(predicted_cases_ui),
    "predicted_cases" = predicted_cases,
    "predicted_min" = predicted_cases_li,
    "predicted_max" = predicted_cases_ui,
    "predicted_min2" = predicted_cases_li2,
    "predicted_max2" = predicted_cases_ui2,
    "deaths" = deaths_by_country[[index]],
    "deaths_c" = cumsum(deaths_by_country[[index]]),
    "estimated_deaths_c" =  cumsum(estimated_deaths),
    "death_min_c" = cumsum(estimated_deaths_li),
    "death_max_c"= cumsum(estimated_deaths_ui),
    "estimated_deaths" = estimated_deaths,
    "death_min" = estimated_deaths_li,
    "death_max"= estimated_deaths_ui,
    "death_min2" = estimated_deaths_li2,
    "death_max2"= estimated_deaths_ui2,
    "rt" = rt,
    "rt_min" = rt_li,
    "rt_max" = rt_ui,
    "rt_min2" = rt_li2,
    "rt_max2" = rt_ui2)
  
  return(data_country)
  
}

CountryForecast <- function(index,country,dates,forecast,prediction,
  estimated.deaths) {
  N <- length(dates)
  N2 <- N + forecast
  
  predicted_cases_forecast <- colMeans(prediction[,1:N2,index])[N:N2]
  predicted_cases_li_forecast <- colQuantiles(prediction[,1:N2,index], 
    probs=.025)[N:N2]
  predicted_cases_ui_forecast <- colQuantiles(prediction[,1:N2,index], 
    probs=.975)[N:N2]
  
  estimated_deaths_forecast <- colMeans(estimated.deaths[,1:N2,index])[N:N2]
  estimated_deaths_li_forecast <- colQuantiles(estimated.deaths[,1:N2,index], 
    probs=.025)[N:N2]
  estimated_deaths_ui_forecast <- colQuantiles(estimated.deaths[,1:N2,index], 
    probs=.975)[N:N2]
  
  times <- as_date(as.character(dates))
  times_forecast <- times[length(times)] + 0:forecast
  data_country_forecast <- data.frame("time" = times_forecast,
    "country" = rep(country, forecast+1),
    "predicted_cases_forecast" = predicted_cases_forecast,
    "predicted_min_forecast" = predicted_cases_li_forecast,
    "predicted_max_forecast"= predicted_cases_ui_forecast,
    "estimated_deaths_forecast" = estimated_deaths_forecast,
    "death_min_forecast" = estimated_deaths_li_forecast,
    "death_max_forecast"= estimated_deaths_ui_forecast)
  
  return(data_country_forecast)
  
}

CountryForecastR0 <- function(index,country,dates,forecast,predictionR0) {
  
  N <- length(dates)
  N2 <- N + forecast
  
  predicted_cases <- colMeans(predictionR0[,1:N2,index])
  predicted_cases_li <- colQuantiles(predictionR0[,1:N2,index], 
    probs=.025)
  predicted_cases_ui <- colQuantiles(predictionR0[,1:N2,index], 
    probs=.975)
  predicted_cases_li2 <- colQuantiles(predictionR0[,1:N2,index], 
    probs=.25)
  predicted_cases_ui2 <- colQuantiles(predictionR0[,1:N2,index], 
    probs=.75)
  
  times <- as_date(as.character(dates))
  times_forecast <- c(times, times[length(times)] + 1:forecast)
  data_country_R0 <- data.frame("time" = times_forecast,
    "country" = rep(country, length(times_forecast)),
    "predicted_cases_c" = cumsum(predicted_cases),
    "predicted_min_c" = cumsum(predicted_cases_li),
    "predicted_max_c" = cumsum(predicted_cases_ui),
    "predicted_cases" = predicted_cases,
    "predicted_min" = predicted_cases_li,
    "predicted_max" = predicted_cases_ui,
    "predicted_min2" = predicted_cases_li2,
    "predicted_max2" = predicted_cases_ui2)
  
  return(data_country_R0)
  
}

# Plotting functions ------------------------------------------------------

make_plots <- function(data_country, covariates_country_long, 
  filename2, figuresDir, country){
  
  data_cases_95 <- data.frame(data_country$time, data_country$predicted_min, 
    data_country$predicted_max)
  names(data_cases_95) <- c("time", "cases_min", "cases_max")
  data_cases_95$key <- rep("nintyfive", length(data_cases_95$time))
  data_cases_50 <- data.frame(data_country$time, data_country$predicted_min2, 
    data_country$predicted_max2)
  names(data_cases_50) <- c("time", "cases_min", "cases_max")
  data_cases_50$key <- rep("fifty", length(data_cases_50$time))
  data_cases <- rbind(data_cases_95, data_cases_50)
  levels(data_cases$key) <- c("ninetyfive", "fifty")
  
  p1 <- ggplot(data_country) +
    geom_bar(data = data_country, aes(x = time, y = reported_cases), 
      fill = "coral4", stat='identity', alpha=0.5) + 
    geom_ribbon(data = data_cases, aes(x = time, ymin = cases_min, 
      ymax = cases_max, fill = key)) +
    geom_line(data = data_country, aes(x = time, y = predicted_cases), 
      col = "blue4") + 
    xlab("") +
    ylab("Daily number of infections") +
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) + 
    scale_fill_manual(name = "", labels = c("50%", "95%"),
      values = c(alpha("deepskyblue4", 0.55), 
        alpha("deepskyblue4", 0.45))) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
      legend.position = "None") + 
    guides(fill=guide_legend(ncol=1))
  
  data_deaths_95 <- data.frame(data_country$time, data_country$death_min, 
    data_country$death_max)
  names(data_deaths_95) <- c("time", "death_min", "death_max")
  data_deaths_95$key <- rep("nintyfive", length(data_deaths_95$time))
  data_deaths_50 <- data.frame(data_country$time, data_country$death_min2, 
    data_country$death_max2)
  names(data_deaths_50) <- c("time", "death_min", "death_max")
  data_deaths_50$key <- rep("fifty", length(data_deaths_50$time))
  data_deaths <- rbind(data_deaths_95, data_deaths_50)
  levels(data_deaths$key) <- c("ninetyfive", "fifty")
  
  
  p2 <-   ggplot(data_country, aes(x = time)) +
    geom_bar(data = data_country, aes(y = deaths, fill = "reported"),
      fill = "coral4", stat='identity', alpha=0.5) +
    geom_ribbon(
      data = data_deaths,
      aes(ymin = death_min, ymax = death_max, fill = key)) +
    geom_line(data = data_country, aes(x = time, y = estimated_deaths), 
      col = "blue4") + 
    xlab("") +
    ylab("Daily number of deaths") +
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
      values = c(alpha("deepskyblue4", 0.55), 
        alpha("deepskyblue4", 0.45))) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
      legend.position = "None") + 
    guides(fill=guide_legend(ncol=1))
  
  
  plot_labels <- c("Complete lockdown", 
    "Public events banned",
    "School closure",
    "Self isolation",
    "Social distancing")
  
  # Plotting interventions
  data_rt_95 <- data.frame(data_country$time, 
    data_country$rt_min, data_country$rt_max)
  names(data_rt_95) <- c("time", "rt_min", "rt_max")
  data_rt_95$key <- rep("nintyfive", length(data_rt_95$time))
  data_rt_50 <- data.frame(data_country$time, data_country$rt_min2, 
    data_country$rt_max2)
  names(data_rt_50) <- c("time", "rt_min", "rt_max")
  data_rt_50$key <- rep("fifty", length(data_rt_50$time))
  data_rt <- rbind(data_rt_95, data_rt_50)
  levels(data_rt$key) <- c("ninetyfive", "fifth")
  
  p3 <- ggplot(data_country) +
    geom_stepribbon(data = data_rt, aes(x = time, ymin = rt_min, ymax = rt_max, 
      group = key,
      fill = key)) +
    geom_hline(yintercept = 1, color = 'black', size = 0.1) + 
    geom_segment(data = covariates_country_long,
      aes(x = value, y = 0, xend = value, yend = max(x)), 
      linetype = "dashed", colour = "grey", alpha = 0.75) +
    geom_point(data = covariates_country_long, aes(x = value, 
      y = x, 
      group = key, 
      shape = key, 
      col = key), size = 2) +
    xlab("") +
    ylab(expression(R[t])) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
      values = c(alpha("seagreen", 0.75), alpha("seagreen", 0.5))) + 
    scale_shape_manual(name = "Interventions", labels = plot_labels,
      values = c(21, 22, 23, 24, 25, 12)) + 
    scale_colour_discrete(name = "Interventions", labels = plot_labels) + 
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b"), 
      limits = c(data_country$time[1], 
        data_country$time[length(data_country$time)])) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="right")
  
  p <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 2))
  save_plot(filename = paste0(figuresDir, country, "_three_panel_", 
    filename2, ".png"), p, base_width = 14)
}

make_single_plot <- function(data_country, data_country_forecast, filename,
  figuresDir, country, logy = TRUE, ymax = 100000){
  
  # data_deaths <- data_country %>%
  #   select(time, deaths, estimated_deaths) %>%
  #   gather("key" = key, "value" = value, -time)
  # 
  # data_deaths_forecast <- data_country_forecast %>%
  #   select(time, estimated_deaths_forecast) %>%
  #   gather("key" = key, "value" = value, -time)
  # 
  # # Force less than 1 case to zero
  # data_deaths$value[data_deaths$value < 1] <- NA
  # data_deaths_forecast$value[data_deaths_forecast$value < 1] <- NA
  # data_deaths_all <- rbind(data_deaths, data_deaths_forecast)
  
  p <- ggplot(data_country) +
    geom_bar(data = data_country, aes(x = time, y = deaths), 
      fill = "coral4", stat='identity', alpha=0.5) + 
    geom_line(data = data_country, aes(x = time, y = estimated_deaths), 
      col = "deepskyblue4") + 
    geom_line(data = data_country_forecast, 
      aes(x = time, y = estimated_deaths_forecast), 
      col = "black", alpha = 0.5) + 
    geom_ribbon(data = data_country, aes(x = time, 
      ymin = death_min, 
      ymax = death_max),
      fill="deepskyblue4", alpha=0.3) +
    geom_ribbon(data = data_country_forecast, 
      aes(x = time, 
        ymin = death_min_forecast, 
        ymax = death_max_forecast),
      fill = "black", alpha=0.35) +
    geom_vline(xintercept = data_country$time[length(data_country$time)], 
      col = "black", linetype = "dashed", alpha = 0.5) + 
    xlab("Date") +
    ylab("Daily number of deaths\n") + 
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) + 
    theme_pubr(base_family="sans") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(fill=guide_legend(ncol=1, reverse = TRUE)) + 
    annotate(geom="text", x=data_country$time[length(data_country$time)]+11,
      y=10000, label="",
      color="black")
  if (logy) {
    p <- p + scale_y_continuous(trans='log10', labels=comma) +
      coord_cartesian(ylim = c(1, ymax), expand = FALSE)
  } else {
    p <- p + scale_y_continuous(labels=comma) +
      coord_cartesian(ylim = c(0, ymax), expand = TRUE)
  }
  print(p)
  
  ggsave(file= paste0(figuresDir, country, "_forecast_", filename, ".png"), 
    p, width = 10)
}

make_comparison_plot <- function(data_country, data_country_forecast, 
  data_country_R0, filename2, figuresDir, country, 
  logy = TRUE, ymax = 10e6, include_cases = TRUE) {
  
  # Remove first zero entry
  data_country_R0[1,2:ncol(data_country_R0)] <- 
    data_country_R0[2,2:ncol(data_country_R0)]
  
  p <- ggplot(data_country) +
    geom_ribbon(data = data_country, aes(x = time, ymin = predicted_min, 
      ymax = predicted_max), fill = "deepskyblue4", alpha=0.35) +
    geom_line(data = data_country, aes(x = time, y = predicted_cases), 
      col = "deepskyblue4") + 
    geom_ribbon(data = data_country_forecast, aes(x = time, 
      ymin = predicted_min_forecast, ymax = predicted_max_forecast), 
      fill = "deepskyblue4", alpha=0.35) +
    geom_line(data = data_country_forecast, aes(x = time, 
      y = predicted_cases_forecast), col = "deepskyblue4") + 
    geom_ribbon(data = data_country_R0, aes(x = time, 
      ymin = predicted_min, ymax = predicted_max), 
      fill = "black", alpha=0.35) +
    geom_line(data = data_country_R0, aes(x = time, 
      y = predicted_cases), col = "black") +
    geom_vline(xintercept = data_country$time[length(data_country$time)], 
      col = "black", linetype = "dashed", alpha = 0.5) + 
    xlab("Date") +
    ylab("Daily number of infectons\n") + 
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) + 
    theme_pubr(base_family="sans") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # +
  # guides(fill=guide_legend(ncol=1, reverse = TRUE)) + 
  # annotate(geom="text", x=data_country$time[length(data_country$time)]+11,
  # y=10000, label="",
  # color="black")
  if (include_cases) {
    p <- p + geom_bar(data = data_country, aes(x = time, y = reported_cases), 
      fill = "coral4", stat='identity', alpha=0.5)
  }
  
  if (logy) {
    p <- p + scale_y_continuous(trans='log10', labels=comma) +
      coord_cartesian(ylim = c(1, ymax), expand = FALSE)
    fileTag <- "_R0_log_"
  } else {
    p <- p + scale_y_continuous(labels=comma) +
      coord_cartesian(ylim = c(0, ymax), expand = TRUE)
    fileTag <- "_R0_linear_"
  }
  
  print(p)
  
  ggsave(file= paste0(figuresDir, country, fileTag, filename, ".png"), 
    p, width = 10)
  
}


plot_covariate_effects <- function(resultsDir, out) {
  
  alpha = data.frame(as.matrix(out$alpha))
  plot_labels <- c("School Closure",
    "Self Isolation",
    "Public Events",
    "First Intervention",
    "Lockdown", 'Social distancing')
  colnames(alpha) = plot_labels
  first.intervention = alpha[,c(1,2,3,5,6)] + alpha[,4]
  data1 = mcmc_intervals_data(first.intervention,prob=.95,
    transformation=function(x) 1-exp(-x),point_est="mean")
  data1$type = "First Intervention"
  
  data2 = mcmc_intervals_data(alpha[,c(1,2,3,5,6)],  prob = .95,
    transformation=function(x) 1-exp(-x),point_est="mean")
  data2$type = "Later Intervention"
  # data = rbind(rbind(data2[6,],data1),data2[1:5,])
  data = rbind(data1,data2[1:5,])
  colQuantiles(data.matrix(first.intervention),probs=c(.025,.975))
  
  #data$type[1] = "First Intervention"
  
  levels(data$parameter) = gsub("t(", "", levels(data$parameter), fixed=TRUE)
  levels(data$parameter) = gsub(")", "", levels(data$parameter), fixed=TRUE)
  data$parameter = (as.character(data$parameter))
  
  no_point_est <- all(data$point_est == "none")
  x_lim <- range(c(data$ll, data$hh))
  x_range <- diff(x_lim)
  x_lim[1] <- x_lim[1] - 0.05 * x_range
  x_lim[2] <- x_lim[2] + 0.05 * x_range
  layer_vertical_line <- if (0 > x_lim[1] && 0 < x_lim[2]) {
    vline_0(color = "gray90", size = 0.5)
  } else {
    geom_ignore()
  }
  args_outer <- list(mapping = aes_(x = ~ll, xend = ~hh, y = ~parameter, 
    yend = ~parameter)) #, color = bayesplot::get_color("mid"))
  args_inner <- list(mapping = aes_(x = ~l, xend = ~h, y = ~parameter, 
    yend = ~parameter), size = 2, show.legend = FALSE)
  args_point <- list(mapping = aes_(x = ~m, y = ~parameter), 
    data = data, size = 4, shape = 21)
  
  args_point$color <- "blue" #get_color("dark_highlight")
  
  point_func <- geom_point
  layer_outer <- do.call(geom_segment, args_outer)
  layer_inner <- do.call(geom_segment, args_inner)
  layer_point <- do.call(point_func, args_point)
  
  data$parameter = factor(as.character(data$parameter),
    levels=plot_labels[order(plot_labels)[6:1]])
  # data = data[order(-data$m),]
  p = ggplot(data) +theme_pubr() +  
    geom_point(aes(x=m,y=parameter,colour=type),position = position_dodge(-.5)) + 
    geom_linerange(aes(xmin=ll,xmax=hh,y=parameter,colour=type),
      position = position_dodge(-.5)) + 
    scale_x_continuous(breaks=seq(0,1,.25),
      labels = c("0%\n(no effect on transmissibility)",
        "25%","50%","75%","100%\n(ends transmissibility)"),
      expand=c(0.005,0.005),expression(paste("Relative % reduction in  ",R[t])))  +
    scale_colour_manual(name = "", #labels = c("50%", "95%"),
      values = c(("coral4"), ("seagreen"))) + 
    
    geom_vline(xintercept=1,colour="darkgray") +
    scale_y_discrete("Governmental intervention\n") +
    #geom_vline(xintercept=0,colour="darkgray") + 
    theme(plot.margin = margin(0, 2, 0, .5, "cm"))
  
  #+ guides(fill=guide_legend(nrow=2))
  p    
  
  ggsave(filename = paste0(resultsDir , "covars-alpha-reduction.png"),
    p,height=4,width=8)
  
}
