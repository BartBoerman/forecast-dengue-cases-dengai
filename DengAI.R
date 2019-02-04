####################################################################################
# Required packages                                                                #   
####################################################################################
# Data wrangling
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org") # month()
if(!require(zoo)) install.packages("zoo", repos = "http://cran.us.r-project.org") # na.locf()
if(!require(RcppRoll)) install.packages("RcppRoll", repos = "http://cran.us.r-project.org") # roll_mean
# Functions to streamline the model training process for complex regression and classification problems. 
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
# Calculate evaluation metrics
if(!require(Metrics)) install.packages("Metrics", repos = "http://cran.us.r-project.org")
####################################################################################
# Get data                                                                         #   
####################################################################################
train.df <- read.csv("dengue_features_train.csv",
                     stringsAsFactors = FALSE,
                     na.strings = "NA")
train.labels.df <- read.csv("dengue_labels_train.csv",
                            stringsAsFactors = FALSE,
                            na.strings = "NA")
train.df <- train.df %>% inner_join(train.labels.df, 
                                    by = c("year" = "year", "city" = "city", "weekofyear" = "weekofyear")) 
rm(train.labels.df)
test.df <- read.csv("dengue_features_test.csv",
                    stringsAsFactors = FALSE,
                    na.strings = "NA")
####################################################################################
# Missing values                                                                   #   
####################################################################################
missing.org <- train.df %>% 
                    select_if(function(x) any(is.na(x))) %>% 
                    summarise_all(funs(sum(is.na(.)))) %>%
                    gather(key = "column",value = "count_missing") %>%
                    arrange(desc(count_missing))
missing_by_year.org <- train.df %>%
                    group_by(city, year) %>% 
                    summarise_all(funs(sum(is.na(.))))  %>%
                    arrange(city, year)
impute.mean <- function(x) {
  replace(x, is.na(x), mean(x, na.rm = TRUE))
}
####################################################################################
# Data wrangling                                                                   #   
####################################################################################
features.org = names(train.df %>% select(matches('ndvi|precipitation|reanalysis|station', ignore.case = TRUE)))
f_data_wrangling <- function(df) {
                        df <- df %>%
                                mutate(week_start_date = as.Date(week_start_date)
                                       ,month = month(week_start_date)
                                       ,city = as.factor(city)
                                ) %>%
                                group_by(year, month) %>% 
                                mutate_at(vars(features.org),.funs = impute.mean) %>%
                                ungroup %>%
                                mutate_at(vars(features.org),.funs = na.locf) %>% # opvullen van overgebleven ontbrekende waarden
                                #na.locf() %>% 
                                arrange(city,week_start_date)
  return(df)
}
train.df <- f_data_wrangling(train.df)

colSums(is.na(train.df)) 
####################################################################################
# Correlation                                                                      #   
####################################################################################
corMatrix = round(cor(train.df %>% select(features.org,total_cases) ) ,2)

highlyCorrelated = findCorrelation(corMatrix, cutoff=0.9,names=TRUE)

# Some features are 100% correlated: 
# reanalysis_dew_point_temp_k with reanalysis_specific_humidity_g_per_kg
# reanalysis_sat_precip_amt_mm with precipitation_amt_mm
# So we can one of them.
# Reanalysis_max_air_temp_k is highly correlated (0.92) with reanalysis_tdtr_k. The latter has the highest correlation with total cases.

# find the moving average / min / max with the highest correlation

for (x in features.org) {
  cPrev <- 0
  cNew  <- 0
  nHigh <- 1
  for (n in seq(1,16)){
    df <- train.df %>% select(year,total_cases, x) %>% filter(!is.na(total_cases)) # train data
    df <- df %>% mutate(var = roll_mean(x = lag(df[[x]],1), n = n, fill = NA, align = "right")) %>% na.omit()
    cNew = round(cor(df$total_cases, df$var),2)
    if(cPrev < cNew) nHigh <- n 
    cPrev = cNew
  }
  xNew = paste(x,"mean",n, sep = "_")
  train.df[[xNew]] <- roll_mean(x = lag(train.df[[x]],1), n = n, fill = NA, align = "right")
  rm(df)
}

for (x in features.org) {
  cPrev = 0
  cNew  <- 0
  nHigh <- 1
  for (n in seq(1,16)){
    df <- train.df %>% select(year,total_cases, x) %>% filter(!is.na(total_cases)) # train data
    df <- df %>% mutate(var = roll_min(x = lag(df[[x]],1), n = n, fill = NA, align = "right")) %>% na.omit()
    cNew = round(cor(df$total_cases, df$var),2)
    if(cPrev < cNew) nHigh <- n 
    cPrev = cNew
  }
  xNew = paste(x,"min", n, sep = "_")
  train.df[[xNew]] <- roll_min(x = lag(train.df[[x]],1), n = n, fill = NA, align = "right")
  rm(df)
}

for (x in features.org) {
  cPrev = 0
  cNew  <- 0
  nHigh <- 1
  for (n in seq(1,16)){
    df <- train.df %>% select(year,total_cases, x) %>% filter(!is.na(total_cases)) # train data
    df <- df %>% mutate(var = roll_max(x = lag(df[[x]],1), n = n, fill = NA, align = "right")) %>% na.omit()
    cNew = round(cor(df$total_cases, df$var),2)
    if(cPrev < cNew) nHigh <- n 
    cPrev = cNew
  }
  xNew = paste(x,"max", n, sep = "_")
  train.df[[xNew]] <- roll_max(x = lag(train.df[[x]],1), n = n, fill = NA, align = "right")
  rm(df)
}


features.use = names(train.df %>% select(matches('ndvi|precipitation|reanalysis|station', ignore.case = TRUE))) # toegevoegde features

corMatrix = round(cor(train.df %>% na.omit %>% select(features.use,total_cases) ) ,2)

highlyCorrelated = findCorrelation(corMatrix, cutoff=0.9,names=TRUE)







