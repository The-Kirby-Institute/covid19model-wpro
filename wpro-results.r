## Produce results for WPRO

source("utils/extract-results.r")

##JOBIDs
dateDir <- "2020-05-07"
JOBID <- "171055"
StanModel <- "base_wpro"
# Load files and results
resultsDir <- paste0("results/DateRep-", dateDir, "/")
figuresDir <- paste0("figures/DateRep-", dateDir, "/")
filename <- paste0(StanModel, '-',JOBID,"-stanfit.Rdata")
# filename <- paste0('base-',JOBID,".Rdata")
print(sprintf("loading: %s",paste0("results/DateRep-", dateDir, "/", 
  filename)))
load(paste0("results/DateRep-", dateDir, "/", filename))

# system(paste0("Rscript covariate-size-effects.r ", "DateRep-", dateDir, "/", filename,
  # '-stanfit.Rdata'))
results <- StanResults(countries,JOBID,out,resultsDir)
plot_covariate_effects(resultsDir, out)

predictionR0 <- resultsR0(stan_data, out)

# Sort out intervention dates
covariates <- read.csv("data/interventions.csv", 
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
    
    country_covariates <- CountryCovariates(country, covariates,
      data_country$rt_max)
    
    write.csv(country_covariates, paste0(figuresDir, 'Interventions-',
      JOBID, '-',countries[[ii]],'.csv'))
    
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
    country <- countries[[ii]]
    
    data_country <- CountryOutputs(ii,country,dates[[ii]],reported_cases,
      deaths_by_country,prediction,estimated.deaths,rt)
    
    data_country_forecast <- CountryForecast(ii,country,dates[[ii]],
      forecast,prediction, estimated.deaths)
    
    make_single_plot(data_country, data_country_forecast, filename,
      figuresDir, country, logy = FALSE, 
      ymax = round(max(data_country$deaths), digits = -1))
}

forecastR0 <- vector(mode = "list", length = nCountries)

for(ii in 1:nCountries) {
  N <- length(dates[[ii]])
  country <- countries[[ii]]
  
  forecastR0[[ii]] <- CountryForecastR0(ii,country,dates[[ii]],
    forecast,predictionR0)
  
   write.csv(forecastR0[[ii]], paste0(figuresDir, 'R0_Results-',JOBID,
      '-',countries[[ii]],'.csv'))
  
  data_countryTest <- CountryOutputs(ii,country,dates[[ii]],reported_cases,
    deaths_by_country, prediction, estimated.deaths, out$Rt)
  
  data_country_forecastTest <- CountryForecast(ii,country,dates[[ii]],
    forecast, prediction, estimated.deaths)
  
  make_comparison_plot(data_countryTest, data_country_forecastTest, 
    forecastR0[[ii]], "infections", filename, figuresDir, country, 
    logy = TRUE, include_cases = TRUE, 
    ymax = max(forecastR0[[ii]]$predicted_max))
  
  make_comparison_plot(data_countryTest, data_country_forecastTest, 
    forecastR0[[ii]], "infections", paste0(filename, "linear"), 
    figuresDir, country, logy = FALSE, include_cases = FALSE, 
    ymax = max(forecastR0[[ii]]$predicted_max))
  
  make_comparison_plot(data_countryTest, data_country_forecastTest, 
    forecastR0[[ii]], "deaths", filename, figuresDir, country, 
    logy = TRUE, include_cases = TRUE, 
    ymax = max(forecastR0[[ii]]$death_max))
  
  make_comparison_plot(data_countryTest, data_country_forecastTest, 
    forecastR0[[ii]], "deaths", paste0(filename, "linear"), 
    figuresDir, country, logy = FALSE, include_cases = FALSE, 
    ymax = max(forecastR0[[ii]]$death_max))
}


