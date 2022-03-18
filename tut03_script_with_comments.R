#### TUTORIAL 3 SOLUTIONS ####



### Preparing 
# installing a package: install.packages('packagename')
# install package named `forecast`.
install.packages('forecast')

# load the forecast package using `library` function
library(forecast)


### 1
### (a)
# load the data in Merck.csv

mydata = read.csv('Merck.csv')

# Selecting dates from 1 January 2011 to 31 January 2012

mydata$Date = as.Date(mydata$Date)

sel_sample = (mydata$Date<=as.Date('2012-01-31')) & 
  (mydata$Date>=as.Date('2011-01-01'))

y = mydata$Adj_Close[sel_sample]

### (b)
# difference and lag operators

Dy = diff(y)
# log(a/b) = log(a) - log(b)

r = diff(log(y))


y = as.matrix(y)
Dy = as.matrix(Dy)
r = as.matrix(r)

# name the variables using `colnames()` function
colnames(y) <- "Stock Prices of Merck (MRK)"
colnames(Dy) <- "Changes in Stock Prices"
colnames(r) <- "Log Returns" 

### (c)
# we use time series plot to detect obvious trends, volatility clustering, etc
# this is typically the first thing to do when processing time series data
# Select 2011
dates = mydata$Date[sel_sample]
sel2011 = dates <= as.Date("2011-12-31")

# plot the three time series we generated.


plot(dates[sel2011], y[sel2011], type = 'l',
     xlab = 'Time (2011)', ylab = colnames(y))

plot(dates[sel2011], Dy[sel2011], type = 'l',
     xlab = 'Time (2011)', ylab = colnames(Dy))

plot(dates[sel2011], r[sel2011], type = 'l',
     xlab = 'Time (2011)', ylab = colnames(r))




### (d)
# plot ACFs and PACFs

acf(y)
pacf(y)

acf(Dy)
pacf(Dy)




# Make Comments



### (e)
# Let's fit ARMA(1,1) as an example first. 
# Example: fit = Arima(data, order = c(1, 0, 1))
# Summarize the results: summary(fit)

fit = Arima(Dy, order = c(1, 0 ,1))
summary(fit)


# How to retrive AIC and BIC?
# AIC: fit[['aic']] or fit$aic
# BIC: fit[['bic]] or fit$bic

print(fit$aic)
print(fit$bic)

# Create a matrix with size 25x4. We will 
# save results into this empty matrix:


ic = matrix( nrow = 25, ncol = 4 )

# name the columns of `ic` for simplicity
colnames(ic) = c('p','q','aic','bic')

# Now, we create a nested for loop iterating for both
# q and p from 0 to 4. We will save results to `ic` 

for (p in 0:4) {
  for (q in 0:4) {
    fit = Arima(Dy, order = c(p, 0, q))
    ic[5 *p + q+ 1,] = c(p, q, fit$aic, fit$bic)
  }
}




### (f)
# find small aic and bics

ic[order(ic[,'aic']),][1:5,] # top 5
ic[order(ic[,'bic']),][1:5,] # top 5


# next, we select the models that make top 10 in both lists;
# (not 100% objective)

adq_set = list(c(1, 0, 1), c(1, 0, 2), c(2, 0, 1), c(1, 0, 0))


# EXTRA: try to use auto.arima(data, seasonal = FALSE)

fit.auto = auto.arima(Dy, seasonal = FALSE)
summary(fit.auto)



### (g)
#
# once we have reduced the set, perform residual analysis;
# the "forecast" package simplifies this with the command "checkresiduals"

checkresiduals(Arima(Dy[sel2011], order = c(1,0,1)))

for (pq_order in adq_set) {
  checkresiduals(Arima(Dy[sel2011], order = pq_order))
}


### (h)
# Forecast changes in MRK stock prices in January, 2012.

hrz = sum(sel_sample) - sum(sel2011) # working days in january 2012
xticks = c(sum(sel_sample) - 3 * hrz + c(1, 2 * hrz, 3 * hrz))


fcast = forecast(Arima(Dy[sel2011], order = c(1, 0, 1)) , h = hrz)

plot(fcast, include = 2 * hrz)
lines(sum(sel2011) + 1:hrz, Dy[!sel2011])


for (pq_order in adq_set) {
  fcast = forecast(Arima(Dy[sel2011], order = pq_order) , h = hrz)
  plot(fcast, include = 2 * hrz)
  lines(sum(sel2011) + 1:hrz, Dy[!sel2011])
}


### (i)
# repeat the forecasting steps with the "y" variable using an ARMA(2,1)

hrz = sum(sel_sample) - sum(sel2011)
xticks = c(sum(sel_sample) - 3 * hrz + c(1, 2 * hrz, 3 * hrz))


fcast = forecast(Arima(y[sel2011], order = c(2, 0, 1)) , h = hrz)

plot(fcast, include = 2 * hrz)
lines(sum(sel2011) + 1:hrz, y[!sel2011])




fcast = forecast(Arima(y[sel2011], order = c(2, 1, 1)) , h = hrz)

plot(fcast, include = 2 * hrz)
lines(sum(sel2011) + 1:hrz, y[!sel2011])



### (j)
# Repeat parts (e)-(h) with the return variable (r)












ic = matrix( nrow = 25, ncol = 4 )

colnames(ic) = c('p','q','aic','bic')


for (p in 0:4) {
  for (q in 0:4) {
    fit = Arima(r, order = c(p, 0, q))
    ic[5 * p + q + 1,] = c(p, q, fit$aic, fit$bic)
  }
}

### (f)

ic[order(ic[,'aic']),][1:5,] # top 5
ic[order(ic[,'bic']),][1:5,] # top 5


adq_set = list(c(1, 0, 1), c(1, 0, 2), c(2, 0, 1), c(1, 0, 0))

### (g)

for (i in adq_set) {
  checkresiduals(Arima(r[sel2011], order = i))
}


### (h)

hrz = sum(sel_sample) - sum(sel2011)
xticks = c(sum(sel_sample) - 3 * hrz + c(1, 2 * hrz, 3 * hrz))

fcast = forecast(Arima(r[sel2011], order = c(1,0,1)), h = hrz, 
                 level = c(68, 95))

plot(fcast, include = 2 * hrz, ylab = colnames(r), xaxt = 'n')
lines(sum(sel2011) + 1:hrz, r[!sel2011])
axis(1, at = xticks, dates[xticks])




for (pq_order in adq_set) {
  fcast = forecast(Arima(r[sel2011], order = pq_order), h = hrz, 
                   level = c(68, 95))
  
  plot(fcast, include = 2 * hrz, ylab = colnames(r), xaxt = 'n')
  lines(sum(sel2011) + 1:hrz, r[!sel2011])
  axis(1, at = xticks, dates[xticks])
}
