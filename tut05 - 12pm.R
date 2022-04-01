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

# github.com/alfurka/teaching

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

ic[order(ic[,'aic']),]
ic[order(ic[,'bic']),]

adq_set = ic[order(ic[,'bic']),][1:5,]

adq_set

best_models = adq_set[,'i']

for (model in best_models) {
  checkresiduals(ADF_est[[model]])
}

# (c)
#

adf.test(y)

# Type 2 and Type 3 specifications are
# empirically not distinguishable from
# I(1)


# (d)


plot(dates[-1], diff(y), type = 'l', main = 'log real per capita GDP',
     xlab = 'Quarters', ylab = 'log(GDP)')

TT <- length(diff(y))
fit = Arima(diff(diff(y)), xreg = diff(y)[-TT],
            order = c(2, 0, 0),
            include.mean = T,
            include.drift = F)
summary(fit)


ic = matrix(nrow = 30, ncol = 6)
colnames(ic) = c('i','constant','trend','p','aic','bic')
ADF_est = list()


i = 0
for (cons in c(F, T)) {
  for (p in 0:9) {
    i = i + 1
    ADF_est[[i]] = Arima(diff(diff(y)), xreg = diff(y)[-TT],
                         order = c(p, 0, 0),
                         include.mean = cons,
                         include.drift = F)
    ic[i,] = c(i, cons, 0, p,
               AIC(ADF_est[[i]]),
               BIC(ADF_est[[i]]))
    
    if (cons) {
      i = i + 1
      ADF_est[[i]] = Arima(diff(diff(y)), xreg = diff(y)[-TT],
                           order = c(p, 0, 0),
                           include.mean = cons,
                           include.drift = T)
      ic[i,] = c(i, cons, 1, p,
                 AIC(ADF_est[[i]]),
                 BIC(ADF_est[[i]]))
    }
  }
}

ic[order(ic[,'aic']),]
ic[order(ic[,'bic']),]

adq_set = ic[order(ic[,'bic']),][1:5,]

adq_set

best_models = adq_set[,'i']

for (model in best_models) {
  checkresiduals(ADF_est[[model]])
}


adf.test(diff(y))











# (e)
#
# We infer that GDP unlike diff(GDP) is  not empirically 
# distinguishable from an I(1) processes.
# we will prefer Arima(p, 1, q) models

# (f)

TT = length(y)
ADF_est = list()
ic = matrix(nrow = 200, ncol = 8)
colnames(ic) = c('i','cons','trend','p','d','q','aic','bic')

i = 0
for (cons in c(T, F)) {
  for (p in 0:3) {
    for (d in 0:1) {
      for (q in 0:3) {
        try(silent = T, expr = {
          i = i + 1
          ADF_est[[i]] = Arima(y,
                               order = c(p, d, q),
                               include.mean = cons,
                               include.drift = F)
          ic[i, ] = c(i, cons, 0, p, d, q, 
                      AIC(ADF_est[[i]]),
                      BIC(ADF_est[[i]]))
        })
        
        
        if (cons) {
          try(silent = T, expr = {
            i = i + 1
            ADF_est[[i]] = Arima(y,
                                 order = c(p, d, q),
                                 include.mean = cons,
                                 include.drift = T)
            ic[i, ] = c(i, cons, 1, p, d, q, 
                        AIC(ADF_est[[i]]),
                        BIC(ADF_est[[i]]))
          })
        }
      }
    }
  }
}

ic[order(ic[,'aic']),][1:5,]
ic[order(ic[,'bic']),][1:5,]

