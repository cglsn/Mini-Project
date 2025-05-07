# Loading data
data2025<-Sheffield.Tinsley

# See structure of the object
head(data2025)

# Creation of a data-frame object
library(lubridate)
data.tmp <- data.frame("date"= data2025$date, "value"= data2025$no2,
                       "year"= year(data2025$date), "month" = month(data2025$date), "week" = isoweek(data2025$date), "day"= day(data2025$date))

# compute daily maximum of hourly measurements
library(tidyr)
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
fit_year<-fgev(value_year)
fit_year1<-fgev(value_year, prob=1/10)
fit_year2<-fgev(value_year, prob=1/100)
fit_year
fit_year1
fit_year2
par(mfrow=c(1,4))
plot(fit_year)
plot(fit_year1)
plot(fit_year2)
par(mfrow=c(1,3))
plot(profile(fit_year))
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
par(mfrow=c(1,3))
plot(profile(fit_month1))
plot(profile(fit_month2))


#Partie POT

# Choisir une séquence de seuils à tester (par ex. percentiles 85 % à 98 %)
thresholds <- quantile(data.tmp$value, probs = seq(0.85, 0.98, by = 0.01), na.rm = TRUE)

# Initialiser des vecteurs pour stocker les estimations
locs <- scales <- shapes <- mean_excess <- numeric(length(thresholds))

# Boucle pour ajuster fpot() à chaque seuil et stocker les paramètres
for (i in seq_along(thresholds)) {
  th <- thresholds[i]
  exceedances <- data.tmp$value[data.tmp$value > th]
  cat("Threshold:", th, "- Number of exceedances:", length(exceedances), "\n")
  
  if (length(exceedances) > 0) {
    fit <- try(fpot(data.tmp$value, threshold = th), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      locs[i] <- th
      scales[i] <- fit$estimate["scale"]
      shapes[i] <- fit$estimate["shape"]
      mean_excess[i] <- mean(exceedances - th, na.rm = TRUE)
    } else {
      locs[i] <- NA; scales[i] <- NA; shapes[i] <- NA; mean_excess[i] <- NA
    }
  } else {
    locs[i] <- NA; scales[i] <- NA; shapes[i] <- NA; mean_excess[i] <- NA
  }
}

par(mfrow = c(2, 2))
plot(thresholds, locs, type = "b", pch = 16, main = "Location parameter", xlab = "Threshold", ylab = "Loc")
plot(thresholds, scales, type = "b", pch = 16, main = "Scale parameter", xlab = "Threshold", ylab = "Scale")
plot(thresholds, shapes, type = "b", pch = 16, main = "Shape parameter", xlab = "Threshold", ylab = "Shape")
plot(thresholds[valid], mean_excess[valid], type = "b", pch = 16, main = "Mean excess function", xlab = "Threshold", ylab = "Mean excess")

mrlplot(data.tmp$value)
tcplot(data.tmp$value,
       tlim = c(10, 150), model = "gpd", nt = 20)



# Fitting standard three-parameter POT to annual maxima
u<-73 # Chosen based on threshold stability plots
npp <- 8766  # 24 × 365.25 = Number of observations per year (hourly data)
fit.pot<-fpot(data.tmp$value, threshold = u, npp=npp)
summary(fit.pot)
str(fit.pot)
par(mfrow=c(1,4))
plot(fit.pot)
par(mfrow=c(1,2))
plot(profile(fit.pot))

# Fit GPD with POT model for 10-year return period (slide 97)
fit_pot_10 <- fpot(data.tmp$value, threshold = u, mper = 10, npp = npp)
par(mfrow=c(1,4))
plot(fit_pot_10)

# Fit GPD with POT model for 100-year return period
fit_pot_100 <- fpot(data.tmp$value, threshold = u, mper = 100, npp = npp)
par(mfrow=c(1,4))
plot(fit_pot_100)


# Inputs
u <- 73               # Threshold
npp <- 8766            # 24 × 365.25
zeta_u <- mean(data.tmp$value > u, na.rm = TRUE)  # Proportion of exceedances
fit_pot <- fpot(data.tmp$value, threshold = u, std.err = TRUE)

# Estimated parameters from the GPD fit
xi <- fit_pot$estimate["shape"]
sigma <- fit_pot$estimate["scale"]

# Function to compute return level for a given T
pot_return_level <- function(T, u, sigma, xi, zeta_u) {
  if (abs(xi) < 1e-6) {
    # Gumbel case
    return(u + sigma * log(T * zeta_u))
  } else {
    return(u + (sigma / xi) * ((T * zeta_u)^xi - 1))
  }
}

# 10- and 100-year return levels
y10 <- pot_return_level(10, u, sigma, xi, zeta_u)
y100 <- pot_return_level(100, u, sigma, xi, zeta_u)

# Output
cat("10-year return level:", round(y10, 2), "\n")
cat("100-year return level:", round(y100, 2), "\n")

# plot residuals against time (monthly)
library(evir)







# Choose a threshold (can be adjusted depending on the data)
u <- 73 

# Fit the GPD (POT) model to the data above the threshold
fit.pot <- gpd.fit(data2025$value, threshold = u)

# Compute the CDF values under the fitted GPD model
F_vals <- pgpd(data2025$value,
               loc = u,
               scale = fit.pot$mle[1],
               shape = fit.pot$mle[2])

# Transform the CDF values to standard normal residuals
# If the model is correct, these residuals should follow a standard normal distribution
residuals <- qnorm(F_vals)

# Plot residuals over time
plot(data2025$date, residuals,
     type = 'p', col = 'darkgreen', pch = 20,
     xlab = 'Time', ylab = 'Residuals',
     main = 'Residuals vs Time')

# Add horizontal reference line at 0
abline(h = 0, col = 'red', lty = 2)




library(dplyr)

# Get top 1% of all hourly values
threshold_extreme <- quantile(data2025$no2, 0.99, na.rm = TRUE)

# Filter extreme values
extremes <- data.tmp %>% filter(value > threshold_extreme)

# Plot frequency of extremes by month
par(mfrow=c(1,1))
barplot(table(extremes$month), 
        main = "Monthly frequency of top 1% NO₂ values",
        ylab = "Count of extreme values",
        xlab = "Month")

library(lubridate)
library(evir)

# Filter observations between November and March
winter_data <- data.tmp %>% filter(month %in% c(11, 12, 1, 2, 3))

# Choose a threshold (e.g. 95th percentile of winter data)
u_winter <- quantile(winter_data$value, 0.95, na.rm = TRUE) # Souldn't we plot the same graphs as before to select the threshold ?

# Fit POT model on winter data
fit_pot_winter <- fpot(winter_data$value, threshold = u_winter, std.err = TRUE)

# Estimate exceedance probability (proportion above threshold)
zeta_u_winter <- mean(winter_data$value > u_winter, na.rm = TRUE)

# Extract parameters
xi_winter <- fit_pot_winter$estimate["shape"]
sigma_winter <- fit_pot_winter$estimate["scale"]

# Return level formula (slide 102)
pot_return_level <- function(T, u, sigma, xi, zeta_u) {
  if (abs(xi) < 1e-6) {
    return(u + sigma * log(T * zeta_u))
  } else {
    return(u + (sigma / xi) * ((T * zeta_u)^xi - 1))
  }
}

# Compute 10- and 100-year return levels
y10 <- pot_return_level(10, u_winter, sigma_winter, xi_winter, zeta_u_winter)
y100 <- pot_return_level(100, u_winter, sigma_winter, xi_winter, zeta_u_winter)

# Print results
cat("Winter POT model (Nov–Mar):\n")
cat("Threshold (u):", round(u_winter, 2), "\n")
cat("10-year return level:", round(y10, 2), "\n")
cat("100-year return level:", round(y100, 2), "\n")

