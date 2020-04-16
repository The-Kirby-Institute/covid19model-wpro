## Produce results for WPRO

source("extract-results.r")

##JOBIDs
doublePHL <- "64383" #Full
PHL-MYS <- "296839" # Quick
no_NPIs <- "725240" # Quick
PHL-MYS-noNCD <- "105368" # Quick
PHL-MYS-noNCD-noNPIs <- "2140698" # Quick

# Load files and results
filename <- paste0('base-',doublePHL,".Rdata")
print(sprintf("loading: %s",paste0("results/",filename)))
load(paste0("results/", filename))

StanResults(countries,doublePHL,out)

# Sort out intervention dates
data_interventions <- read.csv("data_wpro/interventions.csv", 
  stringsAsFactors = FALSE)

if ((length(countries) == 2) & (countries[1] == countries[2])) {
  nCountries <- 1
} else {
  nCountries <- length(countries) 
}

covariates <- TidyCovariates(nCountries, data_interventions)

results <- vector(mode = "list", length = nCountries)

for(ii in 1:nCountries) {
    print(i)
    N <- length(dates[[ii]])
    country <- countries[[ii]]
    print(country)
    
    data_country <- CountryOutputs(ii,country,dates[[ii]],reported_cases,
      deaths_by_country,prediction,estimated.deaths,out$Rt)
    
    country_covariates <- CountryCovariates(country, covariates,data_country$rt_max)
    
    make_plots(data_country = data_country, 
      covariates_country_long = country_covariates,
      filename = filename,
      country = country)
    
    results[[ii]] <- data_country

}

forecast <- 7

for(ii in 1:nCountries) 
    N <- length(dates[[ii]])
    N2 <- N + forecast
    country <- countries[[ii]]
    
    data_country <- CountryOutputs(ii,country,dates[[ii]],reported_cases,
      deaths_by_country,prediction,estimated.deaths,rt)
    
    data_country_forecast <- CountryForecast(ii,country,dates[[ii]],forecast,prediction,
      estimated.deaths)
    
    make_single_plot(data_country, data_country_forecast, filename,
      country)

}