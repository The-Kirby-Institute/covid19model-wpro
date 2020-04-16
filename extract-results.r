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
StanResults <- function(countries,JOBID,out) {
  filename <- paste0('base-',JOBID)
  plot_labels <- c("School Closure",
    "Self Isolation",
    "Public Events",
    "First Intervention",
    "Lockdown", 'Social distancing')
  alpha = (as.matrix(out$alpha))
  colnames(alpha) = plot_labels
  g = (mcmc_intervals(alpha, prob = .9))
  ggsave(sprintf("results/%s-covars-alpha-log.pdf",filename),g,width=4,
    height=6)
  g = (mcmc_intervals(alpha, prob = .9,
    transformations = function(x) exp(-x)))
  ggsave(sprintf("results/%s-covars-alpha.pdf",filename),g,width=4,height=6)
  mu = (as.matrix(out$mu))
  colnames(mu) = countries
  g = (mcmc_intervals(mu,prob = .9))
  ggsave(sprintf("results/%s-covars-mu.pdf",filename),g,width=4,height=6)
  dimensions <- dim(out$Rt)
  Rt = (as.matrix(out$Rt[,dimensions[2],]))
  colnames(Rt) = countries
  g = (mcmc_intervals(Rt,prob = .9))
  ggsave(sprintf("results/%s-covars-final-rt.pdf",filename),g,width=4,
    height=6)
  return(list("alpha" = alpha, "mu" = mu, "Rt" = Rt))
}

# Sort out covariates ----------------------------------------------------
TidyCovariates <- function(ncountries, data_interventions) {  
  covariates <- data_interventions[1:ncountries, c(1,2,3,4,5,6, 7, 8)]
  covariates[is.na(covariates)] <- "31/12/2020" 
  covariates[,2:8] <- lapply(covariates[,2:8], 
    function(x) as.Date(x, format='%d/%m/%Y'))
  
  return(covariates)
  
}

CountryCovariates <- function(country, covariates,rt_ui) {
  # delete these 2 lines
  covariates_country <- covariates[which(covariates$Country == country), 2:8]   
  
  # Remove sport
  covariates_country$sport = NULL 
  covariates_country$travel_restrictions = NULL 
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
  
  data_country <- data.frame("time" = as_date(as.character(dates[[index]])),
    "country" = rep(country, length(dates[[index]])),
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
    "estimated_deaths_forecast" = estimated_deaths_forecast,
    "death_min_forecast" = estimated_deaths_li_forecast,
    "death_max_forecast"= estimated_deaths_ui_forecast)
  
  return(data_country_forecast)
  
}

# Plotting functions ------------------------------------------------------

make_plots <- function(data_country, covariates_country_long, 
  filename2, country){
  
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
  save_plot(filename = paste0("figures/", country, "_three_panel_", 
    filename2, ".pdf"), p, base_width = 14)
}

make_single_plot <- function(data_country, data_country_forecast, filename, 
  country){
  
  data_deaths <- data_country %>%
    select(time, deaths, estimated_deaths) %>%
    gather("key" = key, "value" = value, -time)
  
  data_deaths_forecast <- data_country_forecast %>%
    select(time, estimated_deaths_forecast) %>%
    gather("key" = key, "value" = value, -time)
  
  # Force less than 1 case to zero
  data_deaths$value[data_deaths$value < 1] <- NA
  data_deaths_forecast$value[data_deaths_forecast$value < 1] <- NA
  data_deaths_all <- rbind(data_deaths, data_deaths_forecast)
  
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
    geom_vline(xintercept = data_deaths$time[length(data_deaths$time)], 
      col = "black", linetype = "dashed", alpha = 0.5) + 
    xlab("Date") +
    ylab("Daily number of deaths\n") + 
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) + 
    scale_y_continuous(trans='log10', labels=comma) + 
    coord_cartesian(ylim = c(1, 100000), expand = FALSE) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(fill=guide_legend(ncol=1, reverse = TRUE)) + 
    annotate(geom="text", x=data_country$time[length(data_country$time)]+8, 
      y=10000, label="Forecast",
      color="black")
  print(p)
  
  ggsave(file= paste0("figures/", country, "_forecast_", filename, ".pdf"), 
    p, width = 10)
}