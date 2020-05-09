library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)

covid_raw <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

deaths_raw <- read.csv("https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")

# Process cases
covid_agg_cases <- aggregate(.~ Country.Region, covid_raw,sum)
covid_agg_deaths <- aggregate(.~ Country.Region, deaths_raw,sum)

covid_deaths <- covid_agg_deaths %>%
  select(-(2:4)) %>%
  rename(Country = Country.Region) %>%
  gather("DateRep", "Deaths_c", 2:ncol(.)) %>%
  arrange(Country) %>%
  group_by(Country) %>%
  mutate(Deaths = c(0, diff(Deaths_c))) %>%
  mutate(DateRep = mdy(substring(DateRep, 2))) %>%
  select(DateRep, everything()) %>%
  ungroup() %>%
  tibble()

d_jhu <- covid_agg_cases %>%
  select(-(2:4)) %>%
  rename(Country = Country.Region) %>%
  gather("DateRep", "Cases_c", 2:ncol(.)) %>%
  arrange(Country) %>%
  group_by(Country) %>%
  mutate(Cases = c(0, diff(Cases_c))) %>%
  mutate(DateRep = mdy(substring(DateRep, 2))) %>%
  select(DateRep, everything()) %>%
  left_join(covid_deaths, by = c("DateRep", "Country")) %>%
  select(DateRep, Cases, Deaths, Country) %>%
  ungroup() %>%
  as_tibble()

# # Check numbers
# d_jhu$Cases[d_jhu$DateRep == "2020-01-22"]
# d_jhu$Deaths[d_jhu$DateRep == "2020-02-21"]
# d_jhu$Cases[d_jhu$Country == "Philippines"]
# d_jhu$Deaths[d_jhu$Country == "Philippines"]
# View(d_jhu[d_jhu$Country == "Philippines",])
# d_jhu$Cases[d_jhu$Country == "Malaysia"]
# d_jhu$Deaths[d_jhu$Country == "Malaysia"]
#Manual Entries 

# 3 cases due to imports on 22 January 2020
d_jhu$Cases[d_jhu$Country == "Philippines" & d_jhu$DateRep == "2020-01-22"] <- 3
# We get negative 2 deaths for Philippines on the 19 March so manual edit.
d_jhu$Deaths[d_jhu$Country == "Philippines" & d_jhu$DateRep == "2020-03-19"] <- 0
d_jhu$Deaths[d_jhu$Country == "Philippines" & d_jhu$DateRep == "2020-03-18"] <- 5

write.csv(d_jhu, "data/COVID-19-up-to-date_JHU.csv")
saveRDS(d_jhu, "data/COVID-19-up-to-date_JHU.rds")
