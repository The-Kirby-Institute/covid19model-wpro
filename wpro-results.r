## Produce results for WPRO

source("extract-results.r")

##JOBIDs
dateDir <- "2020-04-24"
JOBID <- "727474"
# doublePHL <- "64383" #Full
# PHL-MYS <- "296839" # Quick
# no_NPIs <- "725240" # Quick
# PHL-MYS-noNCD <- "105368" # Quick
# PHL-MYS-noNCD-noNPIs <- "2140698" # Quick

# Load files and results
resultsDir <- paste0("results/DateRep-", dateDir, "/")
filename <- paste0('base-',JOBID,"-stanfit.Rdata")
print(sprintf("loading: %s",paste0("results/DateRep-", dateDir, "/", filename)))
load(paste0("results/DateRep-", dateDir, "/", filename))

# system(paste0("Rscript covariate-size-effects.r ", "DateRep-", dateDir, "/", filename,
  # '-stanfit.Rdata'))
results <- StanResults(countries,JOBID,out,resultsDir)

# Sort out intervention dates
covariates <- read.csv("data_wpro/interventions.csv", 
  stringsAsFactors = FALSE)

if ((length(countries) == 2) & (countries[1] == countries[2])) {
  nCountries <- 1
} else {
  nCountries <- length(countries) 
}

covariates <- TidyCovariates(nCountries, covariates)

results <- vector(mode = "list", length = nCountries)

for(ii in 1:nCountries) {
    print(ii)
    N <- length(dates[[ii]])
    country <- countries[[ii]]
    print(country)
    
    data_country <- CountryOutputs(ii,country,dates[[ii]],reported_cases,
      deaths_by_country,prediction,estimated.deaths,out$Rt)
    
    write.csv(data_country, paste0(figuresDir, 'SummaryResults-',JOBID,
      '-',countries[[ii]],'.csv'))
    
    country_covariates <- CountryCovariates(country, covariates,data_country$rt_max)
    
    make_plots(data_country = data_country, 
      covariates_country_long = country_covariates,
      filename = filename,
      figuresDir = figuresDir, 
      country = country)
    
    results[[ii]] <- data_country
}

forecast <- 10

for(ii in 1:nCountries) {
    N <- length(dates[[ii]])
    N2 <- N + forecast
    country <- countries[[ii]]
    
    data_country <- CountryOutputs(ii,country,dates[[ii]],reported_cases,
      deaths_by_country,prediction,estimated.deaths,rt)
    
    data_country_forecast <- CountryForecast(ii,country,dates[[ii]],forecast,prediction,
      estimated.deaths)
    
    make_single_plot(data_country, data_country_forecast, filename,
      country, logy = FALSE, ymax = round(max(data_country$deaths), digits = -1))

}
