###################################################################################
# Requered packages and default theme for data visualization                      #
###################################################################################
# Data wrangling
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org") # data wrangling and visualization framework
if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org") # date functions
if(!require(imputeTS)) install.packages("imputeTS", repos = "http://cran.us.r-project.org") # impute time series
if(!require(TTR)) install.packages("TTR", repos = "http://cran.us.r-project.org") # moving averages
# Weather metrics
if(!require(weathermetrics)) install.packages("weathermetrics", repos = "http://cran.us.r-project.org")
# Anomalies
if(!require(tibbletime)) install.packages("tibbletime", repos = "http://cran.us.r-project.org") # anomalize dependency
if(!require(anomalize)) install.packages("anomalize", repos = "http://cran.us.r-project.org")
# Machine learning
if(!require(randomForest)) install.packages("randomForest", repos = "http://cran.us.r-project.org")
if(!require(prophet)) install.packages("prophet", repos = "http://cran.us.r-project.org") 
if(!require(Metrics)) install.packages("Metrics", repos = "http://cran.us.r-project.org")
# Visualisation
if(!require(GGally)) install.packages("GGally", repos = "http://cran.us.r-project.org") # 
if(!require(ggpubr)) install.packages("ggpubr", repos = "http://cran.us.r-project.org")
if(!require(ggthemes)) install.packages("ggthemes", repos = "http://cran.us.r-project.org") 
# Chose theme and update it with desired look and feel
theme_set(theme_light(base_size = 9)) 
theme_update(plot.title = element_text(hjust = 0, size=10 , face = "bold", color = "#7a7a7a" )
             ,plot.subtitle = element_text(hjust = 0, color = "#7a7a7a")   
             ,plot.caption = element_text(hjust = 0, color = "#7a7a7a")
             ,panel.border = element_blank()
             ,axis.line = element_line(colour = "black")
             ,axis.title = element_text(hjust = 0.5, color = "#7a7a7a")
             ,axis.text.y = element_text(color = "#7a7a7a")
             ,axis.text.x = element_text(color = "#7a7a7a")
             ,legend.text = element_text(color = "#7a7a7a", size = 9)
             ,legend.title = element_text(color = "#7a7a7a", size = 9)
)
####################################################################################
# Custom functions                                                                 #
####################################################################################
f_convert_temperature_to_c <- function(x){return (convert_temperature(x,old_metric = "k", new_metric = "c"))}
f_data_wrangling <- function(df) {
  df <- df %>%
    arrange(city,week_start_date) %>%
    mutate(week_start_date = as.Date(week_start_date)
           ,month = month(week_start_date)
           ,city = as.factor(city)
    ) %>%
    mutate_at(vars(features.selected),.funs = na.kalman) %>%
    mutate_at(vars(features.kelvin),.funs = f_convert_temperature_to_c) 
  return(df)
}
f_normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
f_roll_cor <- rollify(~cor(.x, .y,method = "spearman", use = "pairwise.complete.obs"), window = 53)
####################################################################################
# Start data wrangling                                                             #                                                                                                                 #
####################################################################################
# Load train data
train.df <- read.csv("dengue_features_train.csv",
                     stringsAsFactors = FALSE,
                     na.strings = "NA")
train.labels.df <- read.csv("dengue_labels_train.csv",
                            stringsAsFactors = FALSE,
                            na.strings = "NA")
# Load unseen, test, data
test.df <- read.csv("dengue_features_test.csv",
                    stringsAsFactors = FALSE,
                    na.strings = "NA")
# combine in one data set for data wrangling and data analysis 
full.df <- rbind(train.df,test.df) 
full.df <- full.df %>% left_join(train.labels.df, 
                                 by = c("year" = "year", "city" = "city", "weekofyear" = "weekofyear")) 
rm(train.labels.df)
# let's get all features on the same temperature scale, kelvin will be converted to celsius
features.kelvin <- names(full.df %>% select(matches('temp_k|tdtr_k', ignore.case = TRUE)))
# create a list of independent variables, features, available in the data set
features.selected = names(full.df %>% select(matches('ndvi|precipitation|reanalysis|station',ignore.case = TRUE)))
# base variables in each data set
names.base <- c("city", "week_start_date", "year", "month", "weekofyear","total_cases")
####################################################################################
# San Juan: Data wrangling                                                         #
####################################################################################
# select city
full.df.sj <- full.df %>% filter(city == 'sj') 
# Data wrangling
full.df.sj <- f_data_wrangling(full.df.sj) 
# select features and add class tbl_time (the latter is required by anomalize)
full.df.sj <- full.df.sj %>% select(names.base,
                                    humidity = reanalysis_relative_humidity_percent,
                                    temp = station_max_temp_c,
                                    temp_air = reanalysis_min_air_temp_k,
                                    precip = station_precip_mm,
                                    dewp = reanalysis_dew_point_temp_k,
                                    ndvi = ndvi_nw) %>%
  as_tbl_time(week_start_date)
full.df.sj.org <- full.df.sj
by_month.sj <- full.df.sj %>% select(year, month, total_cases) %>%
  group_by(month) %>% 
  summarise(min_total_cases = min(total_cases, na.rm = TRUE),
            max_total_cases = max(total_cases, na.rm = TRUE))  
####################################################################################
# San Juan feature engineering: moving average and interaction variable            #
####################################################################################
f_rolling_sum <- rollify(sum, window = 52)
full.df.sj <- full.df.sj %>% mutate(temp_sma = SMA(full.df.sj$temp, n = 25),
                                    temp_air_sma = SMA(full.df.sj$temp_air, n = 16),
                                    humidity_sma = SMA(full.df.sj$humidity, n = 18),
                                    dewp_sma = SMA(full.df.sj$dewp, n = 18),  
                                    ndvi_sma = SMA(full.df.sj$ndvi, n = 12), 
                                    precip_sum = f_rolling_sum(precip),  
                                    interaction = temp_sma * humidity_sma) 
####################################################################################
# Data wrangling: Iquitos                                                          #
####################################################################################
# select city
full.df.iq <- full.df %>% filter(city == 'iq') 
# Data wrangling
full.df.iq <- f_data_wrangling(full.df.iq) 
# select features and add class tbl_time (the latter is required by anomalize)
full.df.iq <- full.df.iq %>% select(names.base,
                                    humidity = reanalysis_relative_humidity_percent,
                                    temp = station_max_temp_c,
                                    temp_air = reanalysis_min_air_temp_k,
                                    precip = station_precip_mm,
                                    dewp = reanalysis_dew_point_temp_k) %>%
  as_tbl_time(week_start_date)
by_month.iq <- full.df.iq %>% select(year, month, total_cases) %>%
  group_by(month) %>% 
  summarise(min_total_cases = min(total_cases, na.rm = TRUE),
            max_total_cases = max(total_cases, na.rm = TRUE))  
####################################################################################
# Iquitos: feature engineering: moving average and interaction variable            #
####################################################################################
f_rolling_sum <- rollify(sum, window = 6)
full.df.iq <- full.df.iq %>% mutate(temp_sma = SMA(temp, n = 12),
                                    temp_air_sma = SMA(full.df.iq$temp_air, n = 8), 
                                    humidity_sma = SMA(humidity, n = 3), 
                                    precip_sum = f_rolling_sum(precip),
                                    dewp_sma = SMA(dewp, n = 6),
                                    interaction = temp_sma * humidity_sma) 
####################################################################################
# Machine learning: use random forest to classify years as year with outbreak      #
####################################################################################
# determine years with an outbreak of dengue, outbreaks start in second semester
outbreak_year      <- full.df.sj %>%
                                mutate(semester = semester(week_start_date)) %>%
                                group_by(year,semester) %>%
                                summarise(total_cases = max(total_cases)) %>%
                                mutate(outbreak = if_else(total_cases > 125,"Yes","No")) %>%
                                mutate(outbreak = as.factor(outbreak)) %>%
                                filter(semester == 2) %>%
                                select("year","outbreak")
# determine import data points just before outbreak period 
outbreak_data      <- full.df.sj %>%
                                mutate(quarter = quarter(week_start_date)) %>%
                                filter(quarter == 3, year > 1990) %>%
                                group_by(year) %>%
                                summarise(temp = mean(temp),
                                          precip= sum(precip)) 
# combine 
outbreak <- inner_join(outbreak_data,outbreak_year,by = c("year"="year"))
# create train and test set
train_outbreak  <- outbreak  %>% filter(year > 1990 & year < 2006)
# used added 2006 and 2007 to test if tree predicts outbreak in 2007 and not in 2006.
test_outbreak   <- outbreak %>% filter(year >= 2006) %>% select(-matches("cases")) 
# train a simple random forest model
tree.fit <- randomForest(
            outbreak ~. 
            ,data = train_outbreak
)    
# add predictions to data
outbreak$outbreak[outbreak$year >= 2006] <- predict(tree.fit, test_outbreak) 
full.df.sj <- full.df.sj %>% left_join(outbreak, by = c("year","year"))
# outbreak peak periods are august, september and oktober 
full.df.sj$outbreak[full.df.sj$month < 7 & full.df.sj$outbreak == "Yes"] <- "No"
full.df.sj$outbreak[full.df.sj$month > 10 & full.df.sj$outbreak == "Yes"] <- "No"
full.df.sj$outbreak[is.na(full.df.sj$outbreak)] <- "No"
# did not have the time create a model for Iquitos
full.df.iq$outbreak <- "Unkonwn"
####################################################################################
# Exploratory Data Analysis (EDA)                                                  #
####################################################################################
####################################################################################
# Plot anomolies per city                                                          #
####################################################################################
anomalized.sj <- full.df.sj %>% filter(year > 2001 & year <= 2008) %>% na.omit %>%
                                time_decompose(total_cases, 
                                               method = "stl",
                                               frequency = "12 months",
                                               trend = "12 months",
                                               message = FALSE) %>%
                                anomalize(remainder, alpha = 0.2, max_anoms = 0.05,  method = "gesd", verbose = TRUE) 
anomalized.iq <- full.df.iq %>% filter(year > 2001 & year <= 2008) %>% na.omit %>%
                                time_decompose(total_cases, 
                                               method = "stl",
                                               frequency = "12 months",
                                               trend = "12 months",
                                               message = FALSE) %>%
                                anomalize(remainder, alpha = 0.2, max_anoms = 0.05,  method = "gesd", verbose = TRUE) 
p1 <-           anomalized.sj$anomalized_tbl %>% 
                                time_recompose() %>% 
                                plot_anomalies(time_recompose = TRUE) +
                                scale_x_date(date_breaks = "52 weeks", labels = scales::date_format(format = "%Y"), minor_breaks = NULL) +
                                scale_y_continuous(breaks = seq(0,150, by = 50), minor_breaks = NULL) +
                                labs(title = "San Juan", x =  NULL, y = NULL) + 
                                theme(legend.position = "none")
p2 <-           anomalized.iq$anomalized_tbl %>% 
                                time_recompose() %>% 
                                plot_anomalies(time_recompose = TRUE) +
                                scale_x_date(date_breaks = "52 weeks", labels = scales::date_format(format = "%Y"), minor_breaks = NULL) +
                                scale_y_continuous(breaks = seq(0,150, by = 50), minor_breaks = NULL) +
                                labs(title = "Iquitos", x =  NULL, y = NULL) + 
                                theme(legend.position = "none")
figure_anomalies <- ggarrange(p1, p2, ncol = 1, nrow = 2) %>%
                                annotate_figure(top = text_grob("Anomalies", color = "#7a7a7a", face = "bold", size = 10),
                                                bottom = text_grob("", color = "#7a7a7a", hjust = 0.5, face = "italic", size = 9),
                                                fig.lab = NULL, fig.lab.face = "plain")
####################################################################################
# Plot seasonn per city                                                            #
####################################################################################
figure_season_sj <- ggplot(data = anomalized.sj$anomalized_tbl %>% filter (year(week_start_date) %in% c(2006)), aes(x = week_start_date, y = season)) + 
                                geom_smooth(method = "loess", formula = "y ~ x", se = FALSE, size =2) +
                                scale_x_date(date_breaks = "1 month", labels = scales::date_format(format = "%m"), minor_breaks = NULL) +
                                labs(title = "San Juan", subtitle = "season", x =  NULL, y = NULL) +
                                theme(axis.text.y=element_blank(), panel.grid = element_blank(), axis.line.y = element_blank())

figure_season_iq <- ggplot(data = anomalized.iq$anomalized_tbl %>% filter (year(week_start_date) %in% c(2006)), aes(x = week_start_date, y = season)) + 
                                geom_smooth(method = "loess", formula = "y ~ x", se = FALSE, size =2) + 
                                scale_x_date(date_breaks = "1 month", labels = scales::date_format(format = "%m"), minor_breaks = NULL) +
                                labs(title = "Iquitos", subtitle = "season", x =  NULL, y = NULL) +
                                theme(axis.text.y=element_blank(), panel.grid = element_blank(), axis.line.y = element_blank(),
                                      panel.border = element_rect(colour = "#7a7a7a", fill=NA, size=2) )
pMin.sj <- ggplot(data = by_month.sj, aes(x=month, y= min_total_cases)) +
                  geom_col(position = "dodge") +
                  scale_y_continuous(breaks = seq(0,max(by_month.sj$min_total_cases), by = 1), 
                                     minor_breaks = NULL, limits = c(0,max(by_month.sj$min_total_cases))) +
                  scale_x_continuous(breaks = seq(1,12, by = 1), minor_breaks = NULL) +
                  labs(title = "San Juan", 
                       subtitle = "minimum", 
                       caption = NULL, x =  NULL, y = NULL)
pMax.sj <- ggplot(data = by_month.sj, aes(x=month, y= max_total_cases)) +
                  geom_col(position = "dodge") +
                  scale_y_continuous(breaks = seq(0,max(by_month.sj$max_total_cases), by = 100), 
                                     minor_breaks = NULL, limits = c(0,max(by_month.sj$max_total_cases))) +
                  scale_x_continuous(breaks = seq(1,12, by = 1), minor_breaks = NULL) +
                  labs(title = "San Juan", 
                       subtitle = "maximum", 
                       caption = NULL, x =  NULL, y = NULL)
pMin.iq <- ggplot(data = by_month.iq, aes(x=month, y= min_total_cases)) +
                  geom_col(position = "dodge") +
                  scale_y_continuous(breaks = seq(0,max(by_month.iq$min_total_cases), by = 1), 
                                     minor_breaks = NULL, limits = c(0,max(by_month.iq$min_total_cases))) +
                  scale_x_continuous(breaks = seq(1,12, by = 1), minor_breaks = NULL) +
                  theme(axis.text.y=element_blank(), panel.grid = element_blank(), axis.line.y = element_blank(),
                        panel.border = element_rect(colour = "#7a7a7a", fill=NA, size=2) ) +
                  labs(title = "Iquitos", 
                       subtitle = "minimum", 
                       caption = NULL, x =  NULL, y = NULL)
pMax.iq <- ggplot(data = by_month.iq, aes(x=month, y= max_total_cases)) +
                  geom_col(position = "dodge") +
                  scale_y_continuous(breaks = seq(0,max(by_month.iq$max_total_cases), by = 25), 
                                     minor_breaks = NULL, limits = c(0,max(by_month.iq$max_total_cases))) +
                  scale_x_continuous(breaks = seq(1,12, by = 1), minor_breaks = NULL) +
                  theme(axis.text.y=element_blank(), panel.grid = element_blank(), axis.line.y = element_blank(),
                        panel.border = element_rect(colour = "#7a7a7a", fill=NA, size=2) ) +
                  labs(title = "Iquitos", 
                       subtitle = "maximum", 
                       caption = NULL, x =  NULL, y = NULL)
####################################################################################
# Plot pairwise comparison                                                         #
####################################################################################
f_ggpairs_scatter <- function(data, mapping, ...){
                              p <- ggplot(data = data, mapping = mapping) + 
                                geom_point(size = 0.1, colour = "blue")
                              return(p)
}
f_ggpairs_density <- function(data, mapping, ...){
                              p <- ggplot(data = data, mapping = mapping) + 
                                geom_density(fill = "darkblue", colour = "darkblue")
                              return(p)
}
columns <- c("total_cases","temp","temp_air","humidity","precip","dewp","ndvi")
columnLabels <- c("Cases","Temp","Air temp","Humidity","Precip","Dewp","Ndvi")
figure_comp <- ggpairs(data = full.df.sj.org %>% select(columns) %>% na.omit, 
                       upper = list(continuous = wrap("cor", method = "spearman", use = "pairwise.complete.obs",color = "grey50")),
                       diag = list(continuous = f_ggpairs_density),
                       lower = list(continuous = f_ggpairs_scatter),
                       columnLabels = columnLabels ,
                       axisLabels = "none") + 
                      labs(title = "Dengue cases and selected features") +
                      theme(legend.position = "none", 
                            axis.ticks = element_blank(), 
                            panel.grid.major = element_blank(),
                            panel.border = element_rect(linetype = "dashed", colour = "grey", fill = NA),
                            panel.spacing = grid::unit(0.5,"lines"),
                            strip.text = element_text(size = 8, colour = "darkblue"))
####################################################################################
# Plot pairwise comparison moving averages and interaction variable                #
####################################################################################
columns <- names(full.df.sj %>% select(matches("sma|sum|interaction|cases")))
columnLabels <- c("Cases","Temp.","Air temp","Humid.","Dewp.","NDVI","Precip.","Interact.")
figure_comp_sma <- ggpairs(data = full.df.sj %>% select(columns) %>% na.omit, 
                           upper = list(continuous = wrap("cor", method = "spearman", use = "pairwise.complete.obs",color = "grey50")),
                           diag = list(continuous = f_ggpairs_density),
                           lower = list(continuous = f_ggpairs_scatter),
                           columnLabels = columnLabels ,
                           axisLabels = "none") + 
                            labs(title = "Dengue cases and moving averages of selected features") +
                            theme(legend.position = "none", 
                                  axis.ticks = element_blank(), 
                                  panel.grid.major = element_blank(),
                                  panel.border = element_rect(linetype = "dashed", colour = "grey", fill = NA),
                                  panel.spacing = grid::unit(0.5,"lines"),
                                  strip.text = element_text(size = 8, colour = "darkblue"))
####################################################################################
# Plot scatter plot                                                                #
####################################################################################
pData <- full.df.sj %>% filter(!is.na(total_cases)) %>% mutate(quarter = quarter(week_start_date)) %>% na.omit 
figure_interaction_scatter <- ggplot(data = pData) + 
                                    geom_point(aes(x = total_cases, y = interaction, colour = quarter) ) +
                                    scale_y_continuous(breaks = seq(0,3000, by = 100)) +
                                    scale_x_continuous(breaks = seq(0,500, by = 100)) +
                                    labs(title = "Large outbreaks above 2500", 
                                         x =  "", 
                                         y = "Temp sma * by humidity sma",
                                         color = "Quarter") 
figure_dewp_scatter <- ggplot(data = pData) + 
                              geom_point(aes(x = total_cases, y = dewp_sma, colour = quarter) ) +
                              scale_x_continuous(breaks = seq(0,500, by = 100)) +
                              labs(title = "Large outbreaks above 21.5", 
                                   x =  "", 
                                   y = "Dewp point sma",
                                   color = "Quarter")
figure_precipitation_scatter <- ggplot(data = pData) + 
                              geom_point(aes(x = total_cases, y = precip_sum, colour = quarter) ) +
                              scale_x_continuous(breaks = seq(0,500, by = 100)) +
                              labs(title = "A plain lagging value might be better a better fit", 
                                   x =  "", 
                                   y = "Precipitation rolling sum",
                                   color = "Quarter")
figure_temp_scatter <- ggplot(data = pData) + 
                              geom_point(aes(x = total_cases, y = temp_sma, colour = quarter) ) +
                              scale_x_continuous(breaks = seq(0,500, by = 100)) +
                              labs(title = "Large outbreaks above 30.5", 
                                   x =  "", 
                                   y = "Temp sma",
                                   color = "Quarter")
figure_scatter <- ggarrange(figure_interaction_scatter, 
                            figure_dewp_scatter,
                            figure_precipitation_scatter,
                            figure_temp_scatter,
                            ncol = 2,
                            nrow = 2,
                            common.legend = TRUE)
####################################################################################
# Plot rolling correlation                                                         #
####################################################################################
full.df.sj <- full.df.sj %>% mutate(temp_rcor = f_roll_cor(temp_sma, total_cases))
LegendTitle <- "Time series:"
LineType <- c('solid', 'dashed')
LineColor <- c('blue', 'darkblue')
NameS1 <-"Dengue cases"
NameS2 <-"Rolling correlation"
pData <- full.df.sj %>% select(week_start_date, s1 = total_cases, s2 = temp_rcor) %>% na.omit
figure_cor_roll <- ggplot(data = pData) + 
                          geom_line(aes(x = week_start_date, y = f_normalize(s1), linetype = NameS1, color = NameS1)) + 
                          geom_line(aes(x = week_start_date, y = s2, linetype = NameS2, color = NameS2)) +
                          scale_color_manual(name = LegendTitle, values = LineColor) +
                          scale_linetype_manual(name = LegendTitle, values = LineType) +
                          scale_x_date(date_breaks = "2 years", labels = scales::date_format(format = "%Y")) +
                          scale_y_continuous(breaks = waiver()) +  
                          labs(title = "Rolling correlation total cases per week with temperature sma", 
                               x =  "", 
                               y = "") +
                          theme(legend.direction = 'horizontal', legend.position = 'bottom')
####################################################################################
# Machine learning: forecast total dengue cases per week with Prophet              #
####################################################################################
f_prophet_model <- function(df, cutOffDate, horizon = 5){
                                  city <- unique(df$city)
                                  cutOffDate <-as.Date(cutOffDate)
                                  cap <- ceiling(IQR(df$total_cases, na.rm = T) * 1.5) 
                                  floor <- if_else(unique(df$city) == "sj",20,4)  
                                  df$y <- df$total_cases
                                  df$ds <- df$week_start_date
                                  df$y[df$y > cap] <- NA 
                                  df$cap <- cap
                                  df$floor <- floor 
                                  train.df <- df %>% filter(week_start_date <= cutOffDate) %>% na.omit 
                                  m <- prophet(weekly.seasonality = FALSE, 
                                               daily.seasonality = FALSE, 
                                               yearly.seasonality = if_else(city == "sj",5,5),
                                               changepoint.prior.scale = if_else(city == "sj",0.001,0.001),
                                               mode = if_else(city == "sj","additive","additive"),
                                               growth = "logistic") 
                                  m <- add_regressor(m, if_else(city == "sj","temp_sma","temp_air_sma"))
                                  m <- add_regressor(m, if_else(city == "sj","humidity_sma","humidity_sma"))
                                  m <- add_regressor(m, if_else(city == "sj","precip_sum","precip_sum"))
                                  m <- add_regressor(m, if_else(city == "sj","outbreak","precip_sum"))
                                  m <- fit.prophet(m,train.df)
                                  test.df <- df %>% filter(week_start_date > cutOffDate & 
                                                             week_start_date <= (cutOffDate + years(horizon))) %>% 
                                    mutate(ds = week_start_date, y = 0)
                                  
                                  forecast <- predict(m, test.df) %>%
                                    mutate(week_start_date = as.Date(ds)) %>%
                                    select(week_start_date, yhat, yhat_upper) %>%
                                    inner_join(df,  by = "week_start_date")
                                  m$forecast <- forecast 
                                  return(m)  
}
####################################################################################
# Plot: forecast and validation                                                    #
####################################################################################
f_figure_forecast <- function(forecast){ 
                                  mae <- round(mae(forecast$total_cases[!is.na(forecast$total_cases)],
                                                   forecast$yhat[!is.na(forecast$total_cases)]),2)
                                  city <- if_else(unique(forecast$city) == "sj","San Juan","Iquitos")
                                  Title <- paste(city, ", mean absolute error on train data: ",mae,sep = "")
                                  LegendTitle = ""
                                  LineType = c("solid","dashed")
                                  LineColor = c("darkblue", "darkblue")
                                  figure_forecast <- ggplot(data = forecast) +
                                    geom_line(aes(x = ds, y = yhat,
                                                  linetype = "Forecast",
                                                  color = "Forecast")) +
                                    scale_color_manual(name = LegendTitle, values = LineColor) +
                                    scale_linetype_manual(name = LegendTitle, values = LineType) +
                                    scale_x_date(date_breaks = "1 year", 
                                                 labels = scales::date_format(format = "%Y"), minor_breaks = NULL) +
                                    labs(title = Title, x =  "", y = "...") +
                                    theme(legend.direction = 'horizontal', legend.position = 'bottom')
                                  return(figure_forecast)
}
f_figure_validate <- function(figure_forecast){ 
                                  figure_validate <- figure_forecast +
                                    geom_line(aes(x = ds, y = total_cases,
                                                  linetype = "Observed",
                                                  color = "Observed"), na.rm = TRUE)
                                  return(figure_validate)
}
forecast.sj <- f_prophet_model(full.df.sj, "2003-04-20", horizon = 5)$forecast
forecast.sj$yhat <- if_else(forecast.sj$outbreak == "Yes",forecast.sj$yhat * 3,forecast.sj$yhat)
p_sj <- f_figure_validate(f_figure_forecast(forecast.sj))
forecast.iq <- f_prophet_model(full.df.iq, "2005-06-20", horizon = 4)$forecast
p_iq <- f_figure_validate(f_figure_forecast(forecast.iq))
figure_validate <- ggarrange(p_sj, p_iq, ncol = 1, nrow = 2, common.legend = TRUE)
####################################################################################
# Create submission file with forecast                                             #
####################################################################################
cutOffDate.sj <- full.df.sj %>% filter(!is.na(total_cases)) %>% slice(which.max(week_start_date)) %>% pull(week_start_date)
forecast.sj <- f_prophet_model(full.df.sj, cutOffDate.sj, horizon = 10)$forecast %>%
                              mutate(yhat = if_else(outbreak == "Yes",yhat * 3,yhat)) %>%
                              mutate(month = month(ds)) %>% 
                              inner_join(by_month.sj, by = "month") %>% 
                              mutate(yhat = if_else(yhat < min_total_cases,min_total_cases,round(yhat))) %>% 
                              select(-c(names(by_month.sj)))
cutOffDate.iq <- full.df.iq %>% filter(!is.na(total_cases)) %>% slice(which.max(week_start_date)) %>% pull(week_start_date)
forecast.iq <- f_prophet_model(full.df.iq, cutOffDate.iq, horizon = 10)$forecast %>% 
                              mutate(month = month(ds)) %>% 
                              inner_join(by_month.iq, by = "month") %>% 
                              mutate(yhat = if_else(yhat < min_total_cases,min_total_cases,round(yhat))) %>% 
                              select(-c(names(by_month.iq)))
submission <-  rbind(select(forecast.sj, c("city","ds","yhat")) , select(forecast.iq, c("city","ds","yhat"))) %>% 
                            rename(week_start_date = ds, total_cases = yhat) %>%
                            mutate(year = lubridate::year(week_start_date), weekofyear = lubridate::isoweek(week_start_date)) %>%
                            select(-week_start_date)

submission_format <- read_csv("submission_format.csv")
submission_format$total_cases <- NULL
submission_format <- inner_join(submission_format, submission, by = c("city","year","weekofyear"))
write.csv(file = "submission.csv", x = submission_format, row.names = FALSE, quote = TRUE)





  






