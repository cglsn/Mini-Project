# Loading data
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


#J'abandonne pour l'instant le plot des résidus après des heures de recherche et de bug...


#Une fois le problème réglé avec l'assistant, faire de même pour mensuel


# Choisir une séquence de seuils à tester (par ex. percentiles 85 % à 98 %)
thresholds <- quantile(daily_maxima$value, probs = seq(0.85, 0.98, by = 0.01), na.rm = TRUE)

# Initialiser des vecteurs pour stocker les estimations
locs <- scales <- shapes <- mean_excess <- numeric(length(thresholds))

# Boucle pour ajuster fpot() à chaque seuil et stocker les paramètres
for (i in seq_along(thresholds)) {
  th <- thresholds[i]
  fit <- try(fpot(daily_maxima$value, threshold = th), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    locs[i] <- th
    scales[i] <- fit$estimate["scale"]
    shapes[i] <- fit$estimate["shape"]
    mean_excess[i] <- mean(daily_maxima$value[daily_maxima$value > th] - th)
  } else {
    locs[i] <- NA
    scales[i] <- NA
    shapes[i] <- NA
    mean_excess[i] <- NA
  }
}

par(mfrow = c(2, 2))
plot(thresholds, locs, type = "b", pch = 16, main = "Location parameter", xlab = "Threshold", ylab = "Loc")
plot(thresholds, scales, type = "b", pch = 16, main = "Scale parameter", xlab = "Threshold", ylab = "Scale")
plot(thresholds, shapes, type = "b", pch = 16, main = "Shape parameter", xlab = "Threshold", ylab = "Shape")
plot(thresholds, mean_excess, type = "b", pch = 16, main = "Mean excess function", xlab = "Threshold", ylab = "Mean excess")


# Fitting standard three-parameter POT to annual maxima
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

# Return level formula (slide 102)
pot_return_level <- function(T, u, sigma, xi, zeta_u) {
  if (abs(xi) < 1e-6) {
    return(u + sigma * log(T * zeta_u))
  } else {
    return(u + (sigma / xi) * ((T * zeta_u)^xi - 1))
  }
}

# Compute 10- and 100-year return levels
y10 <- pot_return_level(10, u_winter, sigma, xi, zeta_u)
y100 <- pot_return_level(100, u_winter, sigma, xi, zeta_u)

# Print results
cat("Winter POT model (Nov–Mar):\n")
cat("Threshold (u):", round(u_winter, 2), "\n")
cat("10-year return level:", round(y10, 2), "\n")
cat("100-year return level:", round(y100, 2), "\n")


