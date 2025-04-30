# Loading data
data2025<-(load(file = "/Users/camillegilson/Desktop/Risk and environmental sustainability/Mini-Projet/Sheffield.Tinsley_no2.Rdata"))
data2025<-Sheffield.Tinsley

# See structure of the object
head(data2025)

# Creation of a data-frame object
library(lubridate)
data.tmp <- data.frame("date"= data2025$date, "value"= data2025$no2,
                     "year"= year(data2025$date), "month"= month(data2025$date), "week" = isoweek(data2025$date), "day"= day(data2025$date))

# compute daily maximum of hourly measurements
library(dplyr)
daily_maxima <-data.tmp %>%
  group_by(year, month, week, day) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  drop_na() #drop missing observations

# Compute monthly maximum of hourly measurements
library(dplyr)
monthly_maxima <-data.tmp %>%
  group_by(year, month) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  drop_na() #drop missing observations

# Compute weekly maximum of hourly measurements
library(dplyr)
weekly_maxima <-data.tmp %>%
  group_by(year, month, week) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  drop_na() #drop missing observations

# Compute yearly maximum of hourly measurements
library(dplyr)
yearly_maxima <-data.tmp %>%
  group_by(year) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  drop_na() #drop missing observations

# Fitting standard three-parameter GEV to annual maxima
value_year<-yearly_maxima$value
fit_year1<-fgev(value_year, prob=1/10)
fit_year2<-fgev(value_year, prob=1/100)
fit_year1
fit_year2
par(mfrow=c(1,4))
plot(fit_year1)
plot(fit_year2)
par(mfrow=c(1,3))
plot(profile(fit_year1))
plot(profile(fit_year2))


# Fitting standard three-parameter GEV to monthly maxima
value_month<-monthly_maxima$value
fit_month1<-fgev(value_month, prob=1/(10*12))
fit_month2<-fgev(value_month, prob=1/(100*12))
fit_month1
fit_month2
par(mfrow=c(1,4))
plot(fit_month1)
plot(fit_month2)
plot(profile(fit_month1))
plot(profile(fit_month2))

# Fitting standard three-parameter GEV to annual maxima
value_year<-yearly_maxima$value
fit.gev_year1<-fgev(value_year, prob=1/10)
fit.gev_year2<-fgev(value_year, prob=1/100)
fit.gev_year1
fit.gev_year2
par(mfrow=c(1,4))
plot(fit.gev_year1)
plot(fit.gev_year2)
par(mfrow=c(1,3))
plot(profile(fit.gev_year1))
plot(profile(fit.gev_year2))


# Fitting standard three-parameter GEV to monthly maxima
value_month<-monthly_maxima$value
fit.gev_month1<-fgev(value_month, prob=1/(10*12))
fit.gev_month2<-fgev(value_month, prob=1/(100*12))
fit.gev_month1
fit.gev_month2
par(mfrow=c(1,4))
plot(fit.gev_month1)
plot(fit.gev_month2)
par(mfrow=c(1,3))
plot(profile(fit.gev_month1))
plot(profile(fit.gev_month2))

# Fitting standard three-parameter POT to annual maxima
value_year<-yearly_maxima$value
fit.pot_year1<-fpot(value_year, prob=1/10)
fit.pot_year2<-fpot(value_year, prob=1/100)
fit.pot_year1
fit.pot_year2
par(mfrow=c(1,4))
plot(fit.pot_year1)
plot(fit.pot_year2)
par(mfrow=c(1,3))
plot(profile(fit.pot_year1))
plot(profile(fit.pot_year2))


# Fitting standard three-parameter POT to monthly maxima
value_month<-monthly_maxima$value
fit.pot_month1<-fpot(value_month, prob=1/(10*12))
fit.pot_month2<-fpot(value_month, prob=1/(100*12))
fit.pot_month1
fit.pot_month2
par(mfrow=c(1,4))
plot(fit.pot_month1)
plot(fit.pot_month2)
par(mfrow=c(1,3))
plot(profile(fit.pot_month1))
plot(profile(fit.pot_month2))



