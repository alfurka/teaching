library(forecast)
library(dplyr)
library(zoo)


install.packages('aTSA')
library(aTSA)

# 1
#
# (a)
#
# load the data in usdata.csv
mydata = read.csv('usdata.csv')

dates <- as.yearqtr(mydata$obs)
y <- mydata$GDP
r <- mydata$FFR

plot(dates, y, type = 'l', main = 'log real per capita GDP',
     xlab = 'Quarters', ylab = 'log(GDP)')

# (b)
#

# example 
TT <- length(y)
fit = Arima(diff(y), xreg = y[-TT],
      order = c(2, 0, 0),
      include.mean = T,
      include.drift = F)
summary(fit)

# create results matrix

ic = matrix(nrow = 30, ncol = 6)
colnames(ic) = c('i','constant','trend','p','aic','bic')
ADF_est = list()

# answer here

i = 0
for (cons in c(F, T)) {
  for (p in 0:9) {
    i = i + 1
    ADF_est[[i]] = Arima(diff(y), xreg = y[-TT],
                         order = c(p, 0, 0),
                         include.mean = cons,
                         include.drift = F)
    ic[i,] = c(i, cons, 0, p,
               AIC(ADF_est[[i]]),
               BIC(ADF_est[[i]]))
    
    if (cons) {
      i = i + 1
      ADF_est[[i]] = Arima(diff(y), xreg = y[-TT],
                           order = c(p, 0, 0),
                           include.mean = cons,
                           include.drift = T)
      ic[i,] = c(i, cons, 1, p,
                 AIC(ADF_est[[i]]),
                 BIC(ADF_est[[i]]))
    }
  }
}



# (c)
#



# (d)














# (e)
#
# We infer that GDP unlike diff(GDP) is  not empirically 
# distinguishable from an I(1) processes.

# (f)

TT = length(y)
ADF_est = list()
ic = matrix(nrow = 120, ncol = 8)
colnames(ic) = c('i','cons','trend','p','d','q','aic','bic')

