library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats) # For gammaAlt functions
# library(optparse)

## Code addapted to handle Western Pacific countries for WHO Regional 
## Office (WPRO). Additional comments added to translate methods from 
## Imperial report:
## * Seth Flaxman, Swapnil Mishra, Axel Gandy et al. Estimating the number 
## of infections and the impact of nonpharmaceutical interventions on 
## COVID-19 in 11 
## European countries. Imperial College London (30-03-2020)
## doi: https://doi.org/10.25561/77731.

rm(list = ls()) # Don't like doing this but just in case # Restart R

source('utils/read-data.r')
source('utils/process-covariates.r')

# User Options ------------------------------------------------------------
options <- list("include_ncd" = TRUE,
  "npi_on" = TRUE,
  "fullRun" = FALSE,
  "debug" = FALSE,
  "model" = "base_wpro",
  "useJHU" = TRUE,
  "EndDate" = NULL) # EndDated default = NULL

DEBUG <- options$debug

FULL <- options$fullRun

if(DEBUG && FULL) {
  stop("Setting both debug and full run modes at once is invalid")
}

StanModel <- options$model
print(sprintf("Running %s",StanModel))
if(DEBUG) {
  print("Running in DEBUG mode")
} else if (FULL) {
  print("Running in FULL mode")
}

## Reading all data -------------------------------------------------------
# Read which countires to use
countries <- read.csv('data/regions.csv', stringsAsFactors = FALSE)
# Read deaths data for regions
d <- read_obs_data(countries)
if (options$useJHU) {
  d_jhu <- readRDS('data/COVID-19-up-to-date_JHU.rds') %>%
    select(DateRep, Cases, Deaths, Country) %>%
    filter(Country %in% countries$Regions, DateRep >= "2020-01-23") %>%
    mutate(Country = as.character(Country))
  
  d <- filter(d, DateRep <= "2020-01-22") %>%
    bind_rows(d_jhu) %>%
    arrange(Country, DateRep)
  d$Deaths[d$Deaths < 0] <- 0
  
}

if (!is.null(options$EndDate[1])) {
  d <- filter(d, DateRep <= options$EndDate)
}

# Read ifr 
ifr.by.country <- read_ifr_data()

# Read interventions
interventions <- read_interventions(countries)

## Ensure that output directories exist
dateResults <- max(as.Date(d$DateRep,format='%d/%m/%Y'))
resultsDir <- paste0("results/DateRep-",dateResults, "/")
figuresDir <- paste0("figures/DateRep-",dateResults, "/")
dir.create(resultsDir, showWarnings = FALSE, recursive = TRUE)
dir.create(figuresDir, showWarnings = FALSE, recursive = TRUE)

N2 = 100 # increase if you need more forecast

# Initialize inputs ------------------------------------------------------
processed_data <- process_covariates(countries = countries, 
  interventions = interventions, d = d , ifr.by.country = ifr.by.country, 
  N2 = N2)
stan_data = processed_data$stan_data
dates = processed_data$dates
deaths_by_country = processed_data$deaths_by_country
reported_cases = processed_data$reported_cases

# Set-up and run sampling -------------------------------------------------
# stan_data$y = t(stan_data$y)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m = stan_model(paste0('stan-models/',StanModel,'.stan'))

if(DEBUG) {
  fit = sampling(m,data=stan_data,iter=40,warmup=20,chains=2)
} else if (FULL) {
  fit = sampling(m,data=stan_data,iter=1800,warmup=1000,chains=5,thin=1,
    control = list(adapt_delta = 0.95, max_treedepth = 15))
} else { 
  fit = sampling(m,data=stan_data,iter=1000,warmup=500,chains=4,thin=1,
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

countries <- countries$Regions
save.image(paste0(resultsDir, StanModel,'-',JOBID,'.Rdata')) # don't like this
save(fit,out,prediction,dates,reported_cases,deaths_by_country,countries,
  estimated.deaths,estimated.deaths.cf,out,options,stan_data,
  file=paste0(resultsDir, StanModel,'-',JOBID,'-stanfit.Rdata'))

# Visualize results -------------------------------------------------------
print("To visualize results run code in wrpo-results.r script")
