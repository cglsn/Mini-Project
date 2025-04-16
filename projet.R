# Loading data
data2025<-(load(file = "/Users/camillegilson/Desktop/Risk and environmental sustainability/Mini-Projet/Sheffield.Tinsley_no2.Rdata"))
data2025<-Sheffield.Tinsley

# See structure of the object
head(data2025)

# Creation of a data-frame object
library(lubridate)
data.tmp <- data.frame("date"= data2025$date, "value"= data2025$no2,
                       "year"= year(data2025$date), "month"= month(data2025$date), "day"= day(data2025$date))
library(dplyr)
daily_maxima <-data.tmp %>%
  group_by(year, month, day) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  drop_na() #drop missing observations