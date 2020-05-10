# covid19model-wpro
Application of the Imperial College covid19model to countries in the Western Pacific for the WHO Western Pacific Regional Office (WPRO) and to provide modelling support to countries in the region. 

This repository is a fork of the original code used for modelling estimated deaths and cases for COVID19 from Report 13 published by MRC Centre for Global Infectious Disease Analysis, Imperial College London: [Estimating the number of infections and the impact of nonpharmaceutical interventions on COVID-19 in 11 European countries](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-13-europe-npi-impact/) 

The code used is in the wpro_analysis branch (default) off the forked master branch of the Version 3 release. Unnecessary code for our analysis is stored in the `imperial-code/` directory. Please go to the upstream repository to see the original code and for version details and instructions in the associated README.md file. Results for the original model applied to European countries are available [here](https://mrc-ide.github.io/covid19estimates/#/) with the technical model description available [here](https://arxiv.org/abs/2004.11342) and [here](https://github.com/ImperialCollegeLondon/covid19model/blob/master/Technical_description_of_Imperial_COVID_19_Model.pdf)

## Changes to upstream code

For the analysis we made the following key changes to the model and modelling code:

* We changed the serial interval fit to have a Weibull distribution with mean = 5.5 and sd = 2.9 in line with the generation time distribution reported in [Ferretti et al](https://science.sciencemag.org/content/368/6491/eabb6936). 
* We changed the prior for $\mu$ to have a mean of 2.7 to be more in-line with meta-analyses of R0 estimates (for example see [here](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01597-8)).
* We created a separate script to fetch John Hopkins University data as an alternative to ECDC data. ECDC data seems to be lagged by a day for Western Pacific countries. 
* We created a separate results and figures generation script `wpro-results.r` rather than using the base `wpro.r` script. 

## Countries

So far we are applying this code to the following countries:
* Philippines
* Malaysia

# How to run the code

There are two ways to run the code:-
* Open the rstudio project covid19model.Rproj file in rstudio and run/source wpro.r file
* To run from the command line please enter the cloned directory and type `Rscript wpro.r base` in terminal
 
## Run mode settings 
Different run modes incorporated into the `wpro.r` script: :

* include_ncd which specifies if the weighted-IFR estimate includes the relative differences in all-cause mortality between the specific country and where the IFR estimates are from (to account for co-morbidities and environmental factors potentially affecting COVID-19 deaths)
* npi_on to turn on the effect of NPIs
* fullRun which must always be used if you want to obtain reliable results
* debug quick run and extra outputs for debugging 
* model to specify the stan model used
* useJHU set to use data from John Hopkins University rather than the ECDC. ECDC data seems to be a day late for Western Pacific countries. 
* endDate to specify the end date of data to run the model on. If set to `NULL` then runs to the last date in the fetched data. 

# Results 
To produce results we have set-up a separate script `wpro-results.r` which needs to be run/sourced after the main script has been run. Note the JOBID needs to be entered manually into the script to produce the results. We store results and figures for different runs by date of data fetching. 
* The results are stored in two folders results and figures.
* Results has the stored stan fits and data used for plotting
* Figures have the images with daily cases, daily death and Rt for all countries.
