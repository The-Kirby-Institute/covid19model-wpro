# covid19model-wpro
Application of the Imperial College covid19model to countries in the Western Pacific for the WHO Western Pacific Regional Office (WPRO) and to provide modelling support to countries in the region. 

This repository is a fork of the original code used for  modelling estimated deaths and infections for COVID-19 from ["Estimating the effects of non-pharmaceutical interventions on COVID-19 in Europe"](https://www.nature.com/articles/s41586-020-2405-7), Flaxman, Mishra, Gandy et al, Nature, 2020, the published version of the original Report 13 published by MRC Centre for Global Infectious Disease Analysis, Imperial College London: [Estimating the number of infections and the impact of nonpharmaceutical interventions on COVID-19 in 11 European countries](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-13-europe-npi-impact/) 

The code used is in the wpro_analysis branch (default) off the forked master branch of the Version 4 release. Unnecessary code for our analysis is stored in the `imperial-code/` directory. Please go to the upstream repository to see the original code and for version details and instructions in the associated README.md file. Results for the original model applied to European countries are available [here](https://mrc-ide.github.io/covid19estimates/#/) with the technical model description available [here](https://arxiv.org/abs/2004.11342) and [here](https://github.com/ImperialCollegeLondon/covid19model/blob/master/Technical_description_of_Imperial_COVID_19_Model.pdf)

## Version 7 Release [![DOI](https://zenodo.org/badge/250386901.svg)](https://zenodo.org/badge/latestdoi/250386901)

This code is the exact code that was used in Flaxman, Mishra, Gandy et al. "Estimating the effects of non-pharmaceutical interventions on COVID-19 in Europe," Nature, 2020. [https://www.nature.com/articles/s41586-020-2405-7](https://www.nature.com/articles/s41586-020-2405-7)

To run the code from the main folder in Rstudio ``source("base-nature.r")`` or from the command line ``Rscript base-nature.r``.

The code should be run in full mode to obtain results---debug mode is only to check that your environment has the required libraries; results will not be reliable as the MCMC chains will not have converged.

The repository with posterior draws of the model in Flaxman, Mishra, Gandy et al. "Estimating the effects of non-pharmaceutical interventions on COVID-19 in Europe," Nature, 2020. [https://www.nature.com/articles/s41586-020-2405-7](https://www.nature.com/articles/s41586-020-2405-7) is [here](https://github.com/ImperialCollegeLondon/covid19modelnaturedraws).

This code doesn't supersede our earlier model, it is here for everyone to have direct access to code used in Flaxman, Mishra, Gandy et al. "Estimating the effects of non-pharmaceutical interventions on COVID-19 in Europe," Nature, 2020.[https://www.nature.com/articles/s41586-020-2405-7](https://www.nature.com/articles/s41586-020-2405-7).

The instructions for European, Italy, Brazil, and USA code are the same as earlier (Look at version 3, version 4, version 5, version 6).

## Version 6 Release

This is the release related to [report 23](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-23-united-states/), where we use mobility data to estimate situation in all states of the USA. All other code is still the same.

To run this code you can directly run the base-usa.r file or from command line after seting the current directory as the repository directory run the following command `Rscript base-usa.r`

The code shold be run in full mode to obtain any results. Not running full model to estimate anything is not recommended and discouraged. Only full run should be used to get results.

The instructions for European, Italy and Brazil code are same as earlier (Look at version 3, version 4 and version 5). This release is specific to [USA report](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-23-united-states/)

## Version 5 Release

This is the release related to [report 21](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-21-brazil/), where we use mobility data to estimate situation in Brazil. All other code is still the same.

To run this code you can directly run the base-Brazil.r file or from command line after seting the current directory as the repository directory run the following command `Rscript base-Brazil.r`

The code shold be run in full mode to obtain any results. Not running full model to estimate anything is not recommended and discouraged. Only full run should be used to get results.

The instructions for European and Italy code are same as earlier (Look at version 3 and version 4). This release is specific to [Brazil report](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-21-brazil/)

## Version 4 Release

This is the release related to [report 20](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-20-italy/), where we use mobility data to estimate situation in Italy. All other code is still the same.

To run this code you can directly source the base-italy.r file in rstudio inside the project or from command line after setting the current directory as the repository directory run the following command `Rscript base-italy.r base-italy google interventions '~ -1 + residential + transit + averageMobility' '~ -1 + residential + transit + averageMobility'`

The code for scenarios runs only in full mode not in short run or debug mode. Not running full model to estimate anything is not recommended and discouraged. Only full run should be used to get results.

The instructions for European code are below. This release is specific to [Italy report](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-20-italy/)

## Changes to upstream code

For the analysis we made the following key changes to the model and modelling code:

* We changed the serial interval fit to have a Weibull distribution with mean = 5.5 and sd = 2.9 in line with the generation time distribution reported in [Ferretti et al](https://science.sciencemag.org/content/368/6491/eabb6936). 
* We changed the prior for `mu` to have a mean of 2.7 to be more in-line with meta-analyses of R0 estimates (for example see [here](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01597-8)).
* We created a separate script to fetch John Hopkins University data as an alternative to ECDC data. ECDC data seems to be lagged by a day for Western Pacific countries. 
* We created a separate results and figures generation script `wpro-results.r` rather than using the base `wpro.r` script. 

## Countries

So far we are applying this code to the following countries:
* Philippines
* Malaysia

## How to run the code

There are two ways to run the code:
* Open the rstudio project covid19model.Rproj file in rstudio and run/source wpro.r file
* To run from the command line please enter the cloned directory and type `Rscript wpro.r base` in terminal
 
## Run mode settings 
Different run modes incorporated into the `wpro.r` script: :

* include_ncd which specifies if the weighted-IFR estimate includes the relative differences in all-cause mortality between the specific country and where the IFR estimates are from (to account for co-morbidities and environmental factors potentially affecting COVID-19 deaths).
* npi_on to turn on the effect of NPIs.
* fullRun which must always be used if you want to obtain reliable results.
* debug quick run and extra outputs for debugging. 
* model to specify the stan model used.
* useJHU set to use data from John Hopkins University rather than the ECDC. ECDC data seems to be a day late for Western Pacific countries. 
* endDate to specify the end date of data to run the model on. If set to `NULL` then runs to the last date in the fetched data. 

## Results 
To produce results we have written a separate script `wpro-results.r` which needs to be run/sourced after the main `wpro.r` script has been run. Note the JOBID number needs to be entered manually into the script to produce the results. We store results and figures for different runs by last date in the fetched data.
 
* The results are stored in two folders results and figures.
* Results have the stored stan fits and data used for plotting
* Figures have the images with daily cases, daily death and Rt for all countries.
