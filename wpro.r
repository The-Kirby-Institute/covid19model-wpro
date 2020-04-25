library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats) # For gammaAlt functions
library(optparse)

## Code addapted to handle Western Pacific countries for WHO Regional 
## Office (WPRO). Additional comments added to translate methods from 
## Imperial report:
## * Seth Flaxman, Swapnil Mishra, Axel Gandy et al. Estimating the number 
## of infections and the impact of nonpharmaceutical interventions on 
## COVID-19 in 11 
## European countries. Imperial College London (30-03-2020)
## doi: https://doi.org/10.25561/77731.

# User Options ------------------------------------------------------------
countries <- c(
  "Philippines",
  "Malaysia" #,
  # "Laos"
)

options <- list("include_ncd" = TRUE,
  "npi_on" = TRUE,
  "fullRun" = FALSE,
  "debug" = FALSE)

# Commandline options and parsing
parser <- OptionParser()
parser <- add_option(parser, c("-D", "--debug"), action="store_true",
  help="Perform a debug run of the model")
parser <- add_option(parser, c("-F", "--full"), action="store_true",
  help="Perform a full run of the model")
cmdoptions <- parse_args(parser, args = commandArgs(trailingOnly = TRUE), 
  positional_arguments = TRUE)

# Default run parameters for the model
if(is.null(cmdoptions$options$debug) || options$debug) {
  DEBUG = Sys.getenv("DEBUG") == "TRUE"
} else {
  DEBUG = cmdoptions$options$debug
}

if(is.null(cmdoptions$options$full) || options$fullRun) {
  FULL = Sys.getenv("FULL") == "TRUE"
} else {
  FULL = cmdoptions$options$full
}

if(DEBUG && FULL) {
  stop("Setting both debug and full run modes at once is invalid")
}

if(length(cmdoptions$args) == 0) {
  StanModel = 'base'
} else {
  StanModel = cmdoptions$args[1]
}

print(sprintf("Running %s",StanModel))
if(DEBUG) {
  print("Running in DEBUG mode")
} else if (FULL) {
  print("Running in FULL mode")
}

## Reading all data -------------------------------------------------------
d=readRDS('data_wpro/COVID-19-up-to-date.rds')

## Ensure that output directories exist
dateResults <- max(as.Date(d1$DateRep,format='%d/%m/%Y'))
resultsDir <- paste0("results/DateRep-",dateResults, "/")
figuresDir <- paste0("figures/DateRep-",dateResults, "/")
dir.create(resultsDir, showWarnings = FALSE, recursive = TRUE)
dir.create(figuresDir, showWarnings = FALSE, recursive = TRUE)
#dir.create("web/", showWarnings = FALSE, recursive = TRUE)

## get IFR and population from same file
ifr.by.country = read.csv("data_wpro/popt_ifr.csv")
ifr.by.country$country = as.character(ifr.by.country[,2])
# ifr.by.country = read.csv("data_wpro/weighted_fatality.csv")
# ifr.by.country$country = as.character(ifr.by.country[,1])
if (options$include_ncd) {
  ifr.by.country$ifr = ifr.by.country$ifr_NCD
} else {
  ifr.by.country$ifr = ifr.by.country$ifr_noNCD
}

if (length(countries) == 1) {
  ifr.by.country = ifr.by.country[ifr.by.country$country == countries[1],]
  ifr.by.country = rbind(ifr.by.country,ifr.by.country)
}

# Get serial interval
serial.interval = read.csv("data_wpro/serial_interval.csv") 
# Not sure why we serial interval instead of just specifying mean and cv below?? Maybe because it is discretized see page 19 of report. 
# dgammaAlt(100, 6.5, cv = 0.62)

# Start sorting out NPIs which are given by dates of start.
covariates = read.csv('data_wpro/interventions.csv', 
  stringsAsFactors = FALSE)
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
covariates$self_isolating_if_ill <- as.Date(covariates$self_isolating_if_ill, format = "%d.%m.%Y")
covariates$social_distancing_encouraged <- as.Date(covariates$social_distancing_encouraged, format = "%d.%m.%Y")

# If missing put far into future
# covariates[is.na(covariates)] <- "31/12/2020" 

# Makes sure dates are right format
# covariates[,2:8] <- lapply(covariates[,2:8], 
#   function(x) as.Date(x, format='%d/%m/%Y'))

# Hack if only one country
if (length(countries) == 1) {
  covariates = covariates[covariates$Country == countries[1],]
  covariates = rbind(covariates,covariates)
}

## using covariates as dates we want
covariates$schools_universities[covariates$schools_universities > covariates$lockdown] <- covariates$lockdown[covariates$schools_universities > covariates$lockdown]
covariates$public_events[covariates$public_events > covariates$lockdown] <- covariates$lockdown[covariates$public_events > covariates$lockdown]
covariates$social_distancing_encouraged[covariates$social_distancing_encouraged > covariates$lockdown] <- covariates$lockdown[covariates$social_distancing_encouraged > covariates$lockdown]
covariates$self_isolating_if_ill[covariates$self_isolating_if_ill > covariates$lockdown] <- covariates$lockdown[covariates$self_isolating_if_ill > covariates$lockdown]

forecast = 0

N2 = 90 # increase if you need more forecast

# Hack to handle one country
if (length(countries) == 1) {
  countries <- c(countries[1], countries[1])
}

# Initialize inputs ------------------------------------------------------
dates = list()
reported_cases = list()
stan_data = list(M=length(countries),N=NULL,covariate1=NULL,
  covariate2=NULL,covariate3=NULL,covariate4=NULL,covariate5=NULL,
  covariate6=NULL,deaths=NULL,f=NULL,N0=6,cases=NULL,
  SI=serial.interval$fit[1:N2],EpidemicStart = NULL, pop = NULL) 
# N0 = 6 to make it consistent with Rayleigh
deaths_by_country = list()

# Various distributions required for modeling using gammaAlt functions 
# from EnvStats package
mean1 = 5.1; cv1 = 0.86; # infection to onset
mean2 = 18.8; cv2 = 0.45 # onset to death
x1 = rgammaAlt(1e7,mean1,cv1) # infection-to-onset distribution
x2 = rgammaAlt(1e7,mean2,cv2) # onset-to-death distribution
# Combined rate - time from infection to deaths door 
ecdf.saved = ecdf(x1+x2)

for(Country in countries) {
  IFR=ifr.by.country$ifr[ifr.by.country$country == Country]
  
  covariates1 <- covariates[covariates$Country == Country, c(2,3,4,5,6)]
  
  d1_pop = ifr.by.country[ifr.by.country$country==Country,]
  d1=d[d$Countries.and.territories==Country,c(1,5,6,7)]
  d1$date = as.Date(d1$DateRep,format='%d/%m/%Y')
  d1$t = decimal_date(d1$date) 
  d1=d1[order(d1$t),]
  
  # Sort out day of first case and day when 10 deaths are reached
  date_min <- dmy('31/12/2019') 
  if (as.Date(d1$DateRep[1], format='%d/%m/%Y') > as.Date(date_min, 
    format='%d/%m/%Y')){
    print(paste(Country,'In padding'))
    pad_days <- as.Date(d1$DateRep[1], format='%d/%m/%Y') - date_min
    pad_dates <- date_min + days(1:pad_days[[1]]-1)
    padded_data <- data.frame("Countries.and.territories" = rep(Country, 
      pad_days),
      "DateRep" = format(pad_dates, '%d/%m/%Y'),
      "t" = decimal_date(as.Date(pad_dates,format='%d/%m/%Y')),
      "date" = as.Date(pad_dates,format='%d/%m/%Y'),
      "Cases" = as.integer(rep(0, pad_days)),
      "Deaths" = as.integer(rep(0, pad_days)),
      stringsAsFactors=F)
    
    d1 <- bind_rows(padded_data, d1)
  }
  index = which(d1$Cases>0)[1]
  index1 = which(cumsum(d1$Deaths)>10)[1] # also 5
  # Assumed day of seeding of new infections. See page 20 of report. 
  index2 = index1-30 
  
  print(sprintf(paste("First non-zero cases is on day %d, and 30 days", 
    "before 10 deaths is day %d"), index, index2))
  d1=d1[index2:nrow(d1),]
  stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)
  stan_data$pop = c(stan_data$pop, d1_pop$popt)
  
  # Specify covariate being on/off using 0, 1. Assume 1 on and after date 
  # of start
  for (ii in 1:ncol(covariates1)) {
    covariate = names(covariates1)[ii]
    d1[covariate] <- (as.Date(d1$DateRep, 
      format='%d/%m/%Y') >= as.Date(covariates1[1,covariate]))*1  
    # should this be > or >=?
  }
  
  dates[[Country]] = d1$date
  
  ## Hazard function and survival -----------------------------------------
  
  # Hazard estimation for death following infection
  # Gamma distributions for time to infection from onset and time from 
  # onset to death (see parameters below). 
  # Then multiply by IFR to get death probability.
  N = length(d1$Cases)
  print(sprintf("%s has %d days of data",Country,N))
  forecast = N2 - N # number of days to forecast?
  if(forecast < 0) {
    print(sprintf("%s: %d", Country, N))
    print("ERROR!!!! increasing N2")
    N2 = N
    forecast = N2 - N
  }
  
  # IFR is the overall probability of dying given infection
  convolution = function(u) (IFR * ecdf.saved(u))
  
  f = rep(0,N2) # f is the probability of dying on day i given infection
  
  # Discretized fatality
  f[1] = (convolution(1.5) - convolution(0))
  for(i in 2:N2) {
    f[i] = (convolution(i+.5) - convolution(i-.5)) 
  }
  
  # Set-up country epi data and stan inputs -------------------------------
  # y=c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))
  reported_cases[[Country]] = as.vector(as.numeric(d1$Cases))
  deaths=c(as.vector(as.numeric(d1$Deaths)),rep(-1,forecast))
  cases=c(as.vector(as.numeric(d1$Cases)),rep(-1,forecast))
  deaths_by_country[[Country]] = as.vector(as.numeric(d1$Deaths))
  
  # Extend covaraites to forecast days
  covariates2 <- as.data.frame(d1[, colnames(covariates1)])
  # x=1:(N+forecast)
  covariates2[N:(N+forecast),] <- covariates2[N,]
  
  ## Append data to forecast days into stan_data
  stan_data$N = c(stan_data$N,N)
  # stan_data$y = c(stan_data$y,y[1]) # just the index case!
  # stan_data$x = cbind(stan_data$x,x)
  stan_data$covariate1 = cbind(stan_data$covariate1,covariates2[,1])
  stan_data$covariate2 = cbind(stan_data$covariate2,covariates2[,2])
  stan_data$covariate3 = cbind(stan_data$covariate3,covariates2[,3])
  stan_data$covariate4 = cbind(stan_data$covariate4,covariates2[,4])
  stan_data$covariate5 = cbind(stan_data$covariate5,covariates2[,4])
  stan_data$covariate6 = cbind(stan_data$covariate6,covariates2[,5])
  stan_data$f = cbind(stan_data$f,f)
  stan_data$deaths = cbind(stan_data$deaths,deaths)
  stan_data$cases = cbind(stan_data$cases,cases)
  
  stan_data$N2=N2
  stan_data$x=1:N2
  if(length(stan_data$N) == 1) {
    stan_data$N = as.array(stan_data$N)
  }
}

if (options$npi_on) {
  # create the `any intervention` covariate
  stan_data$covariate4 = 1*as.data.frame((stan_data$covariate1+
      stan_data$covariate2+
      stan_data$covariate3+
      stan_data$covariate5+
      stan_data$covariate6) >= 1)
} else {
  # Turn off NPIs
  stan_data$covariate1 = 0 * stan_data$covariate1
  stan_data$covariate2 = 0 * stan_data$covariate2
  stan_data$covariate3 = 0 * stan_data$covariate3
  stan_data$covariate4 = 0 * stan_data$covariate4
  stan_data$covariate5 = 0 * stan_data$covariate5
}

# Check NPI dates
if(DEBUG) {
  for(i in 1:length(countries)) {
    write.csv(
      data.frame(date=dates[[i]],
        `school closure`=stan_data$covariate1[1:stan_data$N[i],i],
        `self isolating if ill`=stan_data$covariate2[1:stan_data$N[i],i],
        `public events`=stan_data$covariate3[1:stan_data$N[i],i],
        `government makes any intervention`=stan_data$covariate4[1:stan_data$N[i],i],
        `lockdown`=stan_data$covariate5[1:stan_data$N[i],i],
        `social distancing encouraged`=stan_data$covariate6[1:stan_data$N[i],i]),
      file=paste0(resultsDir, countries[i], "-check-dates.csv"),row.names=F)
      # file=sprintf("results/%s-check-dates.csv",countries[i]),row.names=F)
  }
}

# Set-up and run sampling -------------------------------------------------
# stan_data$y = t(stan_data$y)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m = stan_model(paste0('stan-models/',StanModel,'.stan'))

if(DEBUG) {
  fit = sampling(m,data=stan_data,iter=40,warmup=20,chains=2)
} else if (FULL) {
  fit = sampling(m,data=stan_data,iter=4000,warmup=2000,chains=4,thin=4,
    control = list(adapt_delta = 0.95, max_treedepth = 10))
} else { 
  fit = sampling(m,data=stan_data,iter=200,warmup=100,chains=4,thin=4,
    control = list(adapt_delta = 0.95, max_treedepth = 10))
}  

# Extract outputs and save results ----------------------------------------
out = rstan::extract(fit)
prediction = out$prediction
estimated.deaths = out$E_deaths
estimated.deaths.cf = out$E_deaths0

JOBID = Sys.getenv("PBS_JOBID")
if(JOBID == "")
  JOBID = as.character(abs(round(rnorm(1) * 1000000)))
print(sprintf("Jobid = %s",JOBID))

save.image(paste0(resultsDir, StanModel,'-',JOBID,'.Rdata'))

save(fit,out,prediction,dates,reported_cases,deaths_by_country,countries,
  estimated.deaths,estimated.deaths.cf,out,covariates,
  file=paste0(resultsDir, StanModel,'-',JOBID,'-stanfit.Rdata'))

# Visualize results -------------------------------------------------------
# Run wrpo-results.r script

# library(bayesplot)
# filename <- paste0(StanModel,'-',JOBID)
# system(paste0("Rscript covariate-size-effects.r ", filename,
#   '-stanfit.Rdata'))
# mu = (as.matrix(out$mu))
# colnames(mu) = countries
# g = (mcmc_intervals(mu,prob = .9))
# ggsave(paste0(resultsDir, filename, "-mu.png"),g,width=4,height=6)
# # ggsave(sprintf("results/%s-mu.png",filename),g,width=4,height=6)
# tmp = lapply(1:length(countries), function(i) (out$Rt_adj[,stan_data$N[i],i]))
# Rt_adj = do.call(cbind,tmp)
# colnames(Rt_adj) = countries
# g = (mcmc_intervals(Rt_adj,prob = .9))
# ggsave(paste0(resultsDir, filename, "-final-rt.png"),g,width=4,height=6)
# ggsave(sprintf("results/%s-final-rt.png",filename),g,width=4,height=6)

# system(paste0("Rscript plot-3-panel.r ", filename,'-stanfit.Rdata'))
# source("plot-3-panel.r")
# summaryOutput <- make_three_panel_plot(paste0(filename, '-stanfit.Rdata'),
  # resultsDir, figuresDir)
# for (ii in 1:length(countries)) {
#   write.csv(summaryOutput[[ii]], paste0(figuresDir, 'SummaryResults-',JOBID,
#     '-',countries[[ii]],'.csv'))
# }
# system(paste0("Rscript plot-forecast.r ", paste0(resultsDir, filename),'-stanfit.Rdata'))
# system(paste0("Rscript make-table.r results/",filename,'-stanfit.Rdata'))
# verify_result <- system(paste0("Rscript web-verify-output.r ", filename,'.Rdata'),intern=FALSE)
# if(verify_result != 0){
#   stop("Verification of web output failed!")
# }



