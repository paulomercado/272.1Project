# Time-series library
library(TSA)
library(tseries)
library(itsmr)
library(forecast)
library(astsa)
library(readr)
library(ggplot2)
library(scales)
# Loading the data
BudgetExp = read.csv("BudgetExp.csv")


data = BudgetExp[169:nrow(BudgetExp), 2]
exp = ts(data)

n = length(exp)
h = 12

train = window(exp, end = c(time(exp)[n-h]))
test  = window(exp, start = c(time(exp)[n-h+1]))

# perform initial tests
adf.test(train)
Box.test(train, type="Ljung")

# transformations
train.l = log(train)
test.l = log(test)
plot(train.l)
adf.test(train.l)
Box.test(train.l, type="Ljung")

train.l.d = diff(train.l)
plot(train.l.d)
adf.test(train.l.d)
Box.test(train.l.d, type="Ljung")

train.l.d.3 = diff(train.l.d, 3)
plot(train.l.d.3)
adf.test(train.l.d.3)
Box.test(train.l.d.3, type = "Ljung")

train.l.d.6 = diff(train.l.d, 6)
plot(train.l.d.6)
adf.test(train.l.d.6)
Box.test(train.l.d.6, type = "Ljung")

train.l.d.12 = diff(train.l.d, 12)
plot(train.l.d.12)
adf.test(train.l.d.12)
Box.test(train.l.d.12, type = "Ljung")

# check ACF, PACF
acf(train.l.d, lag.max = 100)
acf(train.l.d, type = "partial", lag.max = 100)

acf(train.l.d.3, lag.max = 100)
acf(train.l.d.3, type = "partial", lag.max = 100)

acf(train.l.d.6, lag.max = 100)
acf(train.l.d.6, type = "partial", lag.max = 100)

acf(train.l.d.12, lag.max = 100)
acf(train.l.d.12, type = "partial", lag.max = 100)

# fitting model using autoarima

# d = NA, no freq
auto.arima(train.l,
           d = NA,
           D = NA,
           max.p = 5,
           max.q = 5,
           max.order = 5,
           seasonal = TRUE,
           stepwise = FALSE,
           approximation = FALSE,
           trace = TRUE,
           allowdrift = FALSE,
           allowmean = TRUE
)

# d = NA, frequency = 3
train.f3 = ts(train, frequency = 3)
train.l.f3 = log(train.f3)
auto.arima(train.l.f3,
           d = NA,
           D = NA,
           max.p = 5,
           max.q = 5,
           max.order = 5,
           seasonal = TRUE,
           stepwise = FALSE,
           approximation = FALSE,
           trace = TRUE,
           allowdrift = FALSE,
           allowmean = TRUE
)

# d = NA, frequency = 6
train.f6 = ts(train, frequency = 6)
train.l.f6 = log(train.f6)
auto.arima(train.l.f6,
           d = NA,
           D = NA,
           max.p = 5,
           max.q = 5,
           max.order = 5,
           seasonal = TRUE,
           stepwise = FALSE,
           approximation = FALSE,
           trace = TRUE,
           allowdrift = FALSE,
           allowmean = TRUE
)

# d = NA, frequency = 12
train.f12 = ts(train, frequency = 12)
train.l.f12 = log(train.f12)
auto.arima(train.l.f12,
           d = NA,
           D = NA,
           max.p = 5,
           max.q = 5,
           max.order = 5,
           seasonal = TRUE,
           stepwise = FALSE,
           approximation = FALSE,
           trace = TRUE,
           allowdrift = FALSE,
           allowmean = TRUE
)

# fitting model
sarima.112_200_6 = Arima(train.l, order=c(1,1,2), seasonal=list(order=c(2,0,0), period=6))
summary(sarima.112_200_6)

# residual diagnostics
tsdiag(sarima.112_200_6)
jarque.bera.test(residuals(sarima.112_200_6)) 
hist(residuals(sarima.112_200_6), breaks = 20)
qqnorm(residuals(sarima.112_200_6)); qqline(residuals(sarima.112_200_6))

# forecast
forecast_112_200_6 = forecast(sarima.112_200_6,h=h)

#plot
start_date = as.Date("2000-01-01")
dates = seq(from = start_date, by = "month", length.out = length(data))
df = data.frame(Date = dates, Actual = as.numeric(data))

forecast_start_index = length(train) + 1
forecast_time = dates[forecast_start_index:(forecast_start_index + h - 1)]

df_forecast = data.frame(
  Time = forecast_time,
  Forecast = exp(forecast_112_200_6$mean),
  Lower = exp(forecast_112_200_6$lower[,2]),
  Upper = exp(forecast_112_200_6$upper[,2])
)


ggplot() +
  geom_line(data = df, aes(x = Date, y = Actual, color = "Actual"), linewidth = 0.4) +
  geom_line(data = df_forecast, aes(x = Time, y = Forecast, color = "Forecast"), linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(data = df_forecast, aes(x = Time, ymin = Lower, ymax = Upper, fill = "95% CI"), alpha = 0.4) +
  labs(title = "Monthly Government Expenditure Forecast",
       x = "Date", y = "Expenditure (in million PHP)",
       color = "", fill = "") +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_manual(values = c("Actual" = "black", "Forecast" = "blue")) +
  scale_fill_manual(values = c("95% CI" = "gray")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot() +
  geom_line(data = df, aes(x = Date, y = Actual, color = "Actual"), linewidth = 0.4)+
  labs(title = "Monthly Government Expenditure",
       x = "Date", y = "Expenditure (in million PHP)",
       color = "", fill = "") +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_manual(values = c("Actual" = "black")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Check Metrics
accuracy(exp(forecast_112_200_6$mean), test)

