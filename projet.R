# Loading data
data2025<-(load(file = "C:/Users/micae/OneDrive/Bureau/Risk and environnemental sustainability/Mini Projet/Mini-Project/Sheffield.Tinsley_no2.Rdata"))
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

library(evd)
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
par(mfrow=c(1,3))
plot(profile(fit_year))
plot(profile(fit_year1))
plot(profile(fit_year2))

# Fitting standard three-parameter GEV to monthly maxima
value_month<-monthly_maxima$value
fit_month<-fgev(value_month)
fit_month1<-fgev(value_month, prob=1/(10*12))
fit_month2<-fgev(value_month, prob=1/(100*12))
fit_month
fit_month1
fit_month2
par(mfrow=c(1,4))
plot(fit_month)
par(mfrow=c(1,3))
plot(profile(fit_month))
plot(profile(fit_month1))
plot(profile(fit_month2))

# Plot residuals against time (monthly)

# Extract estimated parameters from GEV fit
loc_hat <- fit_month$param[1]
scale_hat <- fit_month$param[2]

# Compute standardized residuals for monthly maxima
std_resids <- (monthly_maxima$value - loc_hat) / scale_hat

# Plot residuals vs time
plot(monthly_maxima$date, std_resids,
     type = "o", pch = 16,
     xlab = "Date", ylab = "Standardized Residuals",
     main = "Residuals vs Time (Monthly Maxima)",
     col = "black")
abline(h = 0, col = "red", lty = 2)
abline(lm(std_resids ~ as.numeric(monthly_maxima$date)), col = "blue", lwd = 2)


# Plot residuals against time (yearly)

# Extract estimated parameters
loc_hat_year <- fit_year$param[1]
scale_hat_year <- fit_year$param[2]

# Compute standardized residuals
std_resids_year <- (yearly_maxima$value - loc_hat_year) / scale_hat_year

# Plot residuals vs time
plot(yearly_maxima$date, std_resids_year,
     type = "o", pch = 16,
     xlab = "Date", ylab = "Standardized Residuals",
     main = "Residuals vs Time (Annual Maxima)",
     col = "black")
abline(h = 0, col = "red", lty = 2)
abline(lm(std_resids_year ~ as.numeric(yearly_maxima$date)), col = "blue", lwd = 2)



#POT Part

# Choose a sequence of thresholds to test (e.g. 85% to 98% percentiles)
thresholds <- quantile(data.tmp$value, probs = seq(0.85, 0.98, by = 0.01), na.rm = TRUE)

# Initialize vectors to store estimates
locs <- scales <- shapes <- mean_excess <- numeric(length(thresholds))

# Loop to adjust fpot() at each threshold and store parameters
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
plot(thresholds, mean_excess, type = "b", pch = 16, main = "Mean excess function", xlab = "Threshold", ylab = "Mean excess")


#Mean Residual Life plot (mrlplot) and Threshold Choice plot (tcplot)
par(mfrow=c(1,3))
mrlplot(data.tmp$value)
tcplot(data.tmp$value,
       tlim = c(10, 150), model = "gpd", nt = 30)


# Fitting POT to annual maxima
u<-71 # Chosen based on threshold stability plots
npp <- 8766  # 24 × 365.25 = Number of observations per year (hourly data)
fit.pot<-fpot(data.tmp$value, threshold = u, npp=npp)
fit.pot
par(mfrow=c(1,4))
plot(fit.pot)
par(mfrow=c(1,2))
plot(profile(fit.pot))

# Fit GPD with POT model for 10-year return period (slide 97)
fit_pot_10 <- fpot(data.tmp$value, threshold = u, mper = 10, npp = npp)
fit_pot_10
par(mfrow=c(1,2))
plot(profile(fit_pot_10))

# Fit GPD with POT model for 100-year return period
fit_pot_100 <- fpot(data.tmp$value, threshold = u, mper = 100, npp = npp)
fit_pot_100
par(mfrow=c(1,2))
plot(profile(fit_pot_100))


# Inputs
u <- 71               # Threshold
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
    return(u + sigma * log(365.25*24*T * zeta_u))
  } else {
    return(u + (sigma / xi) * ((365.25*24*T * zeta_u)^xi - 1))
  }
}

# 10- and 100-year return levels
y10 <- pot_return_level(10, u, sigma, xi, zeta_u)
y100 <- pot_return_level(100, u, sigma, xi, zeta_u)

# Output
cat("10-year return level:", round(y10, 2), "\n")
cat("100-year return level:", round(y100, 2), "\n")

# Plot residuals against time 
u <- 71

# Extract exceedance indices and values
exceed_idx <- which(data.tmp$value > u)
exceedances <- data.tmp$value[exceed_idx]
times_exc <- data.tmp$date[exceed_idx]

# Extract scale parameter from fit.pot
scale_hat <- fit.pot$param["scale"]

# Compute standardized residuals
std_resids <- (exceedances - u) / scale_hat

# Plot residuals vs time
plot(times_exc, std_resids,
     type = "o", pch = 16,
     xlab = "Date", ylab = "Standardized Residuals",
     main = "Residuals vs Time (POT)",
     col = "black")
abline(h = 0, col = "red", lty = 2)  # horizontal line


library(dplyr)

# Get top 1% of all hourly values
threshold_extreme <- quantile(data2025$no2, 0.99, na.rm = TRUE)

# Filter extreme values
extremes <- data.tmp %>% filter(value > threshold_extreme)

# Plot frequency of extremes by month
par(mfrow=c(1,1))
barplot(table(extremes$month), 
        names.arg = month.abb,  # Jan, Feb, Mar, …
        main     = expression(paste("Monthly frequency of top 1% ", NO[2], " values")),
        xlab     = "Month",
        ylab     = "Count of extreme values")


library(lubridate)
library(evir)

# Filter observations between November and March
winter_data <- data.tmp %>% filter(month %in% c(11, 12, 1, 2, 3))

# Choose a threshold (e.g. 95th percentile of winter data)
par(mfrow=c(1,3))
mrlplot(winter_data$value)
tcplot(winter_data$value,
       tlim = c(10, 150), model = "gpd", nt = 30)
u_winter <- quantile(winter_data$value, 0.93, na.rm = TRUE) 
u_winter<- 76

# Fit POT model on winter data
fit_pot_winter <- fpot(winter_data$value, threshold = u_winter, std.err = TRUE,  npp = npp)
fit_pot_winter
par(mfrow=c(1,4))
plot(fit_pot_winter)
par(mfrow=c(1,2))
plot(profile(fit_pot_winter))

# Fit GPD with POT model for 10-year return period in winter (slide 97)
fit_pot_10_winter <- fpot(winter_data$value, threshold = u_winter, mper = 10, npp = npp)
fit_pot_10_winter
par(mfrow=c(1,2))
plot(profile(fit_pot_10_winter))

# Fit GPD with POT model for 100-year return period in winter
fit_pot_100_winter <- fpot(winter_data$value, threshold = u_winter, mper = 100, npp = npp)
fit_pot_100_winter
par(mfrow=c(1,2))
plot(profile(fit_pot_100_winter))

# Estimate exceedance probability (proportion above threshold)
zeta_u_winter <- mean(winter_data$value > u_winter, na.rm = TRUE)

# Extract parameters
xi_winter <- fit_pot_winter$estimate["shape"]
sigma_winter <- fit_pot_winter$estimate["scale"]

# Return level formula (slide 102)
pot_return_level <- function(T, u, sigma, xi, zeta_u) {
  if (abs(xi) < 1e-6) {
    return(u + sigma * log(365.25*24*T * zeta_u))
  } else {
    return(u + (sigma / xi) * ((365.25*24*T * zeta_u)^xi - 1))
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

# Plot residuals against time 

# Extract exceedance indices and values
exceed_idx <- which(data.tmp$value > u_winter)
exceedances <- data.tmp$value[exceed_idx]
times_exc <- data.tmp$date[exceed_idx]

# Extract scale parameter from fit.pot
scale_hat <- fit_pot_winter$param["scale"]

# Compute standardized residuals
std_resids <- (exceedances - u_winter) / scale_hat

# Plot residuals vs time
plot(times_exc, std_resids,
     type = "o", pch = 16,
     xlab = "Date", ylab = "Standardized Residuals",
     main = "Residuals vs Time (winter POT)",
     col = "black")
abline(h = 0, col = "red", lty = 2)  # horizontal line

# Extremal index estimation
tlim <- quantile(data.tmp$value, probs = c(0.7, 0.95), na.rm = TRUE)
par(mfrow = c(1, 4))
for (r_val in c(1, 3, 5, 10)) {
  exiplot(data.tmp$value, tlim = tlim, r = r_val, main = paste("r =", r_val))
}
# Since it is stable accross the different r, we can take r=1
par(mfrow = c(1, 1))
result <- exiplot(data.tmp$value, tlim = tlim, r = 1)
result$x  # thresholds
result$y  # extremal index estimates
# For a threshold of 71, it suggets an extremal index of 0.33

# Fitting POT and allowing for dependence and clustering:
fit_pot_dep <- fpot(data.tmp$value, threshold = u, npp=npp, r=1, cmax=TRUE)
fit_pot_dep

fit_pot_dep_10 <- fpot(data.tmp$value, threshold = u, npp=npp, r=1, mper=10, cmax=TRUE)
fit_pot_dep_10
plot(profile(fit_pot_dep_10))

fit_pot_dep_100 <- fpot(data.tmp$value, threshold = u, npp=npp, r=1, mper=100, cmax=TRUE)
fit_pot_dep_100
plot(profile(fit_pot_dep_100))

# Computation of the return levels:
sigma <- fit_pot_dep$estimate["scale"]
xi <- fit_pot_dep$estimate["shape"]
theta<-0.3261
p_u <- mean(data.tmp$value > u, na.rm = TRUE)
p_u
N_p10 <- 365.25*24*10
N_p100<- 365.25*24*100
x_10 <- u + (sigma/xi)*((p_u*theta*N_p10)^xi -1)  
x_100 <- u + (sigma/xi)*((p_u*theta*N_p100)^xi -1) 
x_10
x_100

#Part 3.3
library (ismev)

# Modelling non-stationarity
library(mgcv)
library(nlme)
library(ismev)

library(ismev)

# Define the time variable 
t <- (monthly_maxima$month - 4 + 12 * (monthly_maxima$year - 1996))  # month index starting at April 1996
trend <- t / (12 * 100)

# Create seasonal components for K = 1
cos1 <- cos(2 * pi * t / 12)
sin1 <- sin(2 * pi * t / 12)

# Create seasonal components for K = 2
cos2 <- cos(4 * pi * t / 12)
sin2 <- sin(4 * pi * t / 12)

# ydat matrices for K = 1 and K = 2
ydat_k1 <- matrix(c(trend, cos1, sin1), ncol = 3, byrow = FALSE)
ydat_k2 <- matrix(c(trend, cos1, sin1, cos2, sin2), ncol = 5, byrow = FALSE)

# M0: Stationary model
fit_M0 <- gev.fit(xdat = monthly_maxima$value)

# M1: K = 1, non-stationary location
fit_M1 <- gev.fit(xdat = monthly_maxima$value,
                  ydat = ydat_k1, mul = 1:3)

# M2: K = 1, non-stationary location + scale
fit_M2 <- gev.fit(xdat = monthly_maxima$value,
                  ydat = ydat_k1, mul = 1:3, sigl = 1:3)

# M3: K = 1, non-stationary location + scale + shape
fit_M3 <- gev.fit(xdat = monthly_maxima$value,
                  ydat = ydat_k1, mul = 1:3, sigl = 1:3, shl = 1:3)

# M4: K = 2, non-stationary location
fit_M4 <- gev.fit(xdat = monthly_maxima$value,
                  ydat = ydat_k2, mul = 1:5)

# M5: K = 2, non-stationary location + scale
fit_M5 <- gev.fit(xdat = monthly_maxima$value,
                  ydat = ydat_k2, mul = 1:5, sigl = 1:5)

# M6: K = 2, non-stationary location + scale + shape
fit_M6 <- gev.fit(xdat = monthly_maxima$value,
                  ydat = ydat_k2, mul = 1:5, sigl = 1:5, shl = 1:5)


extract_model_info <- function(fit, model_name = "Model") {
  logLik <- -fit$nllh
  dev <- 2 * fit$nllh
  n_params <- length(fit$mle)
  AIC <- dev + 2 * n_params
  
  return(data.frame(
    Model = model_name,
    LogLikelihood = round(logLik, 2),
    Deviance = round(dev, 2),
    N_Params = n_params,
    AIC = round(AIC, 2)
  ))
}


results <- rbind(
  extract_model_info(fit_M0, "M0 (stationary)"),
  extract_model_info(fit_M1, "M1 (K=1, loc)"),
  extract_model_info(fit_M2, "M2 (K=1, loc+scale)"),
  extract_model_info(fit_M3, "M3 (K=1, loc+scale+shape)"),
  extract_model_info(fit_M4, "M4 (K=2, loc)"),
  extract_model_info(fit_M5, "M5 (K=2, loc+scale)"),
  extract_model_info(fit_M6, "M6 (K=2, loc+scale+shape)")
)

print(results)


likelihood_ratio_test <- function(fit_simple, fit_complex, name_simple = "Simple", name_complex = "Complex") {
  dev_simple <- 2 * fit_simple$nllh
  dev_complex <- 2 * fit_complex$nllh
  df <- length(fit_complex$mle) - length(fit_simple$mle)
  stat <- dev_simple - dev_complex
  p_value <- pchisq(stat, df = df, lower.tail = FALSE)
  
  cat(sprintf("Likelihood Ratio Test:\n"))
  cat(sprintf("  %s deviance: %.2f\n", name_simple, dev_simple))
  cat(sprintf("  %s deviance: %.2f\n", name_complex, dev_complex))
  cat(sprintf("  Test statistic: %.2f (df = %d)\n", stat, df))
  cat(sprintf("  p-value: %.4f\n", p_value))
}

# Test if M1 (K=1, loc) is significantly better than M0 (stationary)
likelihood_ratio_test(fit_M0, fit_M1, "M0", "M1")

# Test if M2 is significantly better than M1
likelihood_ratio_test(fit_M1, fit_M2, "M1", "M2")

# Test if M3 (K=1, loc+scale+shape) improves over M2 (loc+scale)
likelihood_ratio_test(fit_M2, fit_M3, "M2", "M3")

# Test if M4  improves over M1 
likelihood_ratio_test(fit_M1, fit_M4, "M1", "M4")

# Test if M5  improves over M4
likelihood_ratio_test(fit_M4, fit_M5, "M4", "M5")


# Test if M6  improves over M5
likelihood_ratio_test(fit_M5, fit_M6, "M5", "M6")
