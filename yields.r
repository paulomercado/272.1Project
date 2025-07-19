# Time-series library
library(TSA)
library(tseries)
library(forecast)
library(astsa)
library(fGarch)
library(rugarch)

HistoricalPrices = read.csv("HistoricalPrices.csv")
data <- HistoricalPrices[[3]]                         # Extract 3rd column
data <- data[1:843]                        # Get rows 843 to end
data <- rev(data)                                     # Reverse the values
yield <- ts(data)

n = length(yield)
h = 7

train = window(yield, end = c(time(yield)[n-h]))
test  = window(yield, start = c(time(yield)[n-h+1]))

adf.test(train)
Box.test(train, type="Ljung")

acf(train,  lag.max=100)
pacf(train, lag.max=100)

arma11 <- Arima(train, order=c(1,0,1))
summary(arma11)

res11 <- residuals(arma11)

# check the residuals
acf(res11)
pacf(res11)
Box.test(res11, type="Ljung")
jarque.bera.test(res11)

# check the squared residuals
acf(res11^2)
pacf(res11^2)
Box.test(res11^2, type="Ljung")

# Fitting ARMA(1,1)-GARCH models

# GARCH student's t
spec11 <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "std"
)

arma11garch11 <- ugarchfit(spec = spec11, data = train)
show(arma11garch11)

spec21 <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(2, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "std"
)

arma11garch21 <- ugarchfit(spec = spec21, data = train)
show(arma11garch21)

spec31 <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(3, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "std"
)

arma11garch31 <- ugarchfit(spec = spec31, data = train)
show(arma11garch31)

# EGARCH student's t
e_spec11 <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "std"
)

arma11egarch11 <- ugarchfit(spec = e_spec11, data = train)
show(arma11egarch11)

e_spec21 <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(2, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "std"
)

arma11egarch21 <- ugarchfit(spec = e_spec21, data = train)
show(arma11egarch21)

e_spec31 <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(3, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "std"
)

arma11egarch31 <- ugarchfit(spec = e_spec31, data = train)
show(arma11egarch31)

# arma11egarch31
# Checking standardized residuals
fit.vol.31 = sigma(arma11egarch31)
fit.vol.sr.31 = residuals(arma11egarch31)/fit.vol.31
plot(as.numeric(fit.vol.31), type = "l", main = "Volatility", ylab = "Volatility", xlab = "")

plot(as.numeric(fit.vol.sr.31), type="l", main = "Standardized Residuals", ylab = "Standardized Residuals", xlab = "")
abline(h = 0, col = "red", lty = 2)

acf(fit.vol.sr.31)
pacf(fit.vol.sr.31)
Box.test(fit.vol.sr.31, type="Ljung", lag = 7)

# Check standardized squared residuals
plot(as.numeric(fit.vol.sr.31^2), type="l", main = "Squared Standardized Residuals", ylab = "Squared Standardized Residuals", xlab = "")

acf(fit.vol.sr.31^2)
pacf(fit.vol.sr.31^2)
Box.test(fit.vol.sr.31^2, type="Ljung", lag = 7)

# arma11egarch11
# Checking standardized residuals
fit.vol.11 = sigma(arma11egarch11)
fit.vol.sr.11 = residuals(arma11egarch31)/fit.vol.11
plot(as.numeric(fit.vol.11), type = "l", main = "Volatility", ylab = "Volatility", xlab = "")

plot(as.numeric(fit.vol.sr.11), type="l", main = "Standardized Residuals", ylab = "Standardized Residuals", xlab = "")
abline(h = 0, col = "red", lty = 2)

acf(fit.vol.sr.11)
pacf(fit.vol.sr.11)
Box.test(fit.vol.sr.11, type="Ljung", lag = 7)

# Check standardized squared residuals
plot(as.numeric(fit.vol.sr.11^2), type="l", main = "Squared Standardized Residuals", ylab = "Squared Standardized Residuals", xlab = "")

acf(fit.vol.sr.11^2)
pacf(fit.vol.sr.11^2)
Box.test(fit.vol.sr.11^2, type="Ljung", lag = 7)

# Forecast 7 steps ahead
mase <- function(actual, forecast, train) {
  n <- length(actual)
  
  # Numerator: model MAE on test
  mae_model <- mean(abs(actual - forecast))
  
  # Denominator: in-sample naive MAE (lag-1 difference)
  scale <- mean(abs(diff(train)))
  
  # Return MASE
  return(mae_model / scale)
}

f7_e_11 <- ugarchforecast(fitORspec = arma11egarch11, n.ahead = h)
f7_e_11
f7_e_11_mean <- as.numeric(fitted(f7_e_11))
accuracy(f7_e_11_mean, test)
mase(test, f7_e_11_mean, train)

f7_e_31 <- ugarchforecast(fitORspec = arma11egarch31, n.ahead = h)
f7_e_31
f7_e_31_mean <- as.numeric(fitted(f7_e_31))
accuracy(f7_e_31_mean, test)
mase(test, f7_e_31_mean, train)


yield_dates <- rev(as.Date(HistoricalPrices$Date, format = "%m/%d/%Y"))


start_date <- as.Date("2022-04-01")
start_index <- which(yield_dates == start_date)


yield_dates <- yield_dates[start_index:(start_index + length(yield) - 1)]

forecast_series <- as.numeric(fitted(f7_e_11std))
forecast_sigma <- as.numeric(sigma(f7_e_11std))

start_point <- 700  
forecast_start <- 837 

plot(yield_dates[start_point:843], yield[start_point:843], type = "l",
     main = "Log Returns of 10-Year Government Bond Yields",
     xlab = "Date", ylab = "Log Return", xaxt='n')
axis.Date(1, at = seq(yield_dates[start_point], yield_dates[843], by = "3 months"),
          format = "%b %Y")
forecast_start <- 837

lines(yield_dates[forecast_start:843], forecast_series, col = "red", lwd = 2)
lines(yield_dates[forecast_start:843], forecast_series + forecast_sigma, col = "blue", lty = 2)
lines(yield_dates[forecast_start:843], forecast_series - forecast_sigma, col = "blue", lty = 2)

abline(v = yield_dates[forecast_start - 1], col = "darkgreen", lty = 3)

