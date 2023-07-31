library(forecast)

# 1
#
# (a)
#
# load the data in Merck.csv
mydata <- read.delim("Merck.csv", header = TRUE,  sep = ",")

mydata$Date = as.Date(mydata$Date)

sel_sample <- mydata$Date >= as.Date("2011-01-01") &
              mydata$Date <= as.Date("2012-01-31")

y <- mydata$Adj_Close[sel_sample]
y = as.matrix(y)

# (b)
#
# difference and lag operators
Dy <- diff(y)


r = diff(log(y))

# alternative:

r = log(y[-1]/y[-length(y)])
r = as.matrix(r)

# alternative 2:
r <- as.matrix(log(y[2:nrow(y)]) - log(lag(y[1:nrow(y) - 1])))



colnames(y) <- "Stock Prices of Merck (MRK)"
colnames(Dy) <- "Changes in Stock Prices"
colnames(r) <- "Log Returns" 

# (c)
#
dates = as.Date(mydata$Date[sel_sample])

# examples

plot(dates, y)
plot(dates, y, type = 'l')
plot(dates, y, type = 'l', xlab = "Time (2011)", ylab = colnames(y))


# stylized


sel2011 <- mydata$Date[sel_sample] <= as.Date("2011-12-31")

plot(dates[sel2011], y[sel2011], type = "l", xlab = "Time (2011)",
     ylab = colnames(y))
plot(dates[sel2011], Dy[sel2011], type = "l", xlab = "Time (2011)",
     ylab = colnames(Dy))
plot(dates[sel2011], r[sel2011], type = "l", xlab = "Time (2011)",
     ylab = colnames(r))

# the DGP for y is likely not stationary as its mean appears to vary over time;
# the DGP for Dy seems to have zero mean, but its variance may depend on time


# (d)
#

acf(y[sel2011], main = colnames(y))
pacf(y[sel2011], main = colnames(y))
acf(Dy[sel2011], main = colnames(Dy))
pacf(Dy[sel2011], main = colnames(Dy))

# (e)
#
# example ARMA(2,1)

fit = Arima(Dy, order = c(2, 0, 1))
fit
fit[['aic']]
fit[['bic']]

# for loop: a matrix for results


ic <- matrix( nrow = 25, ncol = 4 )
colnames(ic) <- c("p", "q", "aic", "bic")
ic

# for loop: estimations

for (p in 0:4){
  for (q in 0:4)
  {
    fit_p_q <- Arima(Dy, order = c(p, 0, q))
    ic[p * 5 + q + 1,] = c(p, q, fit_p_q[["aic"]], fit_p_q[["bic"]])
  }
}

ic

# (f)
#

print(ic)

# ordering:

order(ic[,3], decreasing = FALSE)

# create best aic and bic tables: 
ic_aic <- ic[order(ic[,3], decreasing = FALSE),]
ic_bic <- ic[order(ic[,4], decreasing = FALSE),]
print(ic_aic)
print(ic_bic)
# 10 best:

ic_aic <- ic[order(ic[,3], decreasing = FALSE),][1:10,]
ic_bic <- ic[order(ic[,4], decreasing = FALSE),][1:10,]

print(ic_aic)
print(ic_bic)

# select best ones: (no general rule)
# rule of thumbs: 
# 1 - common best models in aic and bic table
# 2 - best models in aic + bic

adq_set = list(c(1, 0, 1), c(1, 0, 2), c(2, 0, 1), c(3, 0, 0))

# (g)
#

# Example: Check residuals

adq_set[[1]]

checkresiduals(Arima(Dy[sel2011], order = adq_set[[1]]))

# Null hypothesis for LJ BOX:
# H0: `residuals are independently distributed` / `residuals are white noise`


# for loop:

for (i in 1:length(adq_set)){
  checkresiduals(Arima(Dy[sel2011], order = adq_set[[i]]))
}

# (h)
#

# forecasting horizon
hrz = sum(sel_sample) - sum(sel2011) 
# dates to be shown in plots (for visualization purpose)
xticks <- c(sum(sel_sample) - 3 * hrz + c(1, 2 * hrz, 3 * hrz))
# we use 2011 values to fit model
actual_Dy <- as.matrix(Dy[!sel2011]) 
# to save the forecast results we create an empty list
fcst_Dy <- vector(mode = "list", length(adq_set)) 
fcst_Dy

for (i in 1:length(adq_set)){
  model_p_q <- adq_set[[i]]
  fit_model = Arima(Dy[sel2011], order = model_p_q)
  fcst_Dy[[i]] <- forecast(fit_model, h = hrz, level = c(68, 95))

  title_p_q <- paste("ARMA(", as.character(model_p_q[1]), ", ",
                              as.character(model_p_q[3]), ")", sep = "")
  
  plot(fcst_Dy[[i]], include = hrz * 2, ylab = colnames(Dy),
       main = title_p_q, xaxt = "n")
  lines(sum(sel2011) + 1:hrz, actual_Dy)
  axis(1, at = xticks, labels = dates[xticks])
}


# (i)
#
# repeat the forecasting steps with the "y" variable using an ARMA(2,1)
actual_y <- as.matrix(y[!sel2011])
fcst_y_lev = forecast(Arima(y[sel2011], order = c(2, 0, 1)),
                      h = hrz, level = c(68, 95) )

plot(fcst_y_lev, include = hrz * 2, ylab = colnames(y),
     main = "ARMA(2, 1)", xaxt = "n", ylim = c(26.1, 33.4))
lines(sum(sel2011) + 1:hrz, actual_y)
axis(1, at = xticks, labels = dates[xticks])

# add first observation
y0 <- mydata$Adj_Close[sum(mydata$Date < as.Date("2011-01-01") - 1)]
y_ext = as.matrix(c(y0, y[sel2011]))

# forecasting:
fcst_y <- vector(mode = "list", length(adq_set))
for (i in 1:length(adq_set)){
  model_p_q <- adq_set[[i]]
  model_p_q[2] = 1
  fit_model = Arima(y_ext, order = model_p_q, include.constant = T)
  fcst_y[[i]] <- forecast(fit_model,
                          h = hrz, level = c(68, 95) )
  
  title_p_q <- paste("ARIMA(", as.character(model_p_q[1]), ", ",
                     as.character(model_p_q[3]), ")", sep = "")
  
  plot(fcst_y[[i]], include = hrz * 2, ylab = colnames(y),
       main = title_p_q, xaxt = "n")
  lines(1 + sum(sel2011) + 1:hrz, actual_y)
  axis(1, at = 1 + xticks, labels = dates[xticks])
}


# (j)
#

# repeat (d)-(h) with `r`, instead of `Dy`

