library(forecast)
library(dplyr)
library(zoo)
library(aTSA)

# 1
#
# (a)
#
# load the data in usdata.csv
mydata <- read.delim("usdata.csv", header = TRUE,  sep = ",")

dates <- as.yearqtr(mydata$obs)
y <- mydata$GDP
r <- mydata$FFR

plot(dates, y, type = "l", xlab = "Time (Quarters)",
     main = "Log Real US GDP per Capita")

# (b)
#



TT <- length(y)
ADF_est <- list()
ic <- matrix( nrow = 30, ncol = 5 )
colnames(ic) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est[[i]] <- Arima(diff(y), xreg = y[-TT],
                          order = c(p, 0, 0),
                          include.mean = as.logical(const),
                          include.drift = F)
    ic[i,] <- c(const, 0, p, ADF_est[[i]]$aic,
                             ADF_est[[i]]$bic)
  }
  
  if (const)
  {
    # only add a specification with trend if there is a
    # constant (i.e., exclude no constant with trend)
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est[[i]] <- Arima(diff(y), xreg = y[-TT],
                            order = c(p, 0, 0),
                            include.mean = as.logical(const),
                            include.drift = T)
      ic[i,] <- c(const, 1, p, ADF_est[[i]]$aic,
                               ADF_est[[i]]$bic)
    }
  }
}

ic_aic <- ic[order(ic[,4]),][1:10,]
ic_bic <- ic[order(ic[,5]),][1:10,]

ic_aic
ic_bic

# AIC: constant + trend and lags 1,...,5
# BIC: constant + trend and lags 1,...,3 or
# const, no trend and lags 1,...,2;
# select top 5 BIC

adq_set <- as.matrix(arrange(as.data.frame(ic_bic[1:5,]),
                   const, trend, p))
adq_idx <- match(data.frame(t(adq_set[, 1:3])),
                 data.frame(t(ic[, 1:3])))

for (i in 1:length(adq_idx)){
  checkresiduals(ADF_est[[adq_idx[i]]])
}


# (c)
#
adq_set
adf.test(y)

# Only "Type 2" and "Type 3" specifications are in our adequate set,
# so we ignore the output related to "Type 1".
# For all specifications in our adequate set, null of unit root
# cannot be rejected.
# The process is not empirically distinguishable from a unit root
# process with a drift (and possibly nonlinear linear growth).

# (d)
#

TT <- length(diff(y))
ADF_est_diff <- list()
ic_diff <- matrix( nrow = 30, ncol = 5 )
colnames(ic_diff) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est_diff[[i]] <- Arima(diff(diff(y)),
                               xreg = diff(y)[-TT],
                               order = c(p, 0, 0),
                               include.mean = as.logical(const),
                               include.drift = F)
    ic_diff[i,] <- c(const, 0, p, ADF_est_diff[[i]]$aic,
                                  ADF_est_diff[[i]]$bic)
  }
  
  if (const)
  {
    # only add a specification with trend if there is a
    # constant (i.e., exclude no constant with trend)
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est_diff[[i]] <- Arima(diff(diff(y)),
                                 xreg = diff(y)[-TT],
                                 order = c(p, 0, 0),
                                 include.mean = as.logical(const),
                                 include.drift = T)
      ic_diff[i,] <- c(const, 1, p, ADF_est_diff[[i]]$aic,
                                    ADF_est_diff[[i]]$bic)
    }
  }
}

ic_aic_diff <- ic_diff[order(ic_diff[,4]),][1:10,]
ic_bic_diff <- ic_diff[order(ic_diff[,5]),][1:10,]


ic_aic_diff
ic_bic_diff

# AIC and BIC agree on the top 5! Great!!!
# The adequate set includes const + no trend for p = 0,...,2
# and const + trend for p = 0, 1.

adq_set_diff <- as.matrix(arrange(as.data.frame(
                      ic_bic_diff[1:5,]), const, trend, p))
adq_idx_diff <- match(data.frame(t(adq_set_diff[, 1:3])),
                      data.frame(t(ic_diff[, 1:3])))

for (i in 1:length(adq_idx_diff)){
  checkresiduals(ADF_est_diff[[adq_idx_diff[i]]])
}


adq_set_diff
adf.test(diff(y))

# unit root is rejected at very small significance level
# for all specifications

# (e)
#
# We infer that GDP is not empirically distinguishable from an
# I(1) processes.

# (f)
#
# We will consider models for p = 0,..,3; q = 0,...,3; d = 0, 1
# and either with or without constant and/or trend terms.
# There are 96 models to estimate altogether. We use the Arima function
# in a nested for loop. There are two caveats to deal with.
# 1. The Arima function with d=1 will only specify an intercept
#    when setting the include.drift = T option. Since we want
#    to include a linear growth term (i.e. "t" as a regressor)
#    in the differenced specification, we need to pass it as
#    an exogenous variable to the Arima function through
#    the "xreg" option. However, with d=1, Arima will difference
#    whatever data we pass this way, so we need to cumulatively
#    sum "t" before passing it.
# 2. Some specifications will be so bad that MLE will run into
#    numerical problems and return an error. We want to ignore
#    these specifications in the least disruptive way possible.
#    The way to do it is to embed the Arima function as an
#    argument to the "try" function with the "silent = T" opt.

# Example try():

try(silent = T, expr = {
  2+2
})

try(silent = T, expr = {
  zzz+'a'
})

zzz+'a'

# Loop

TT <- length(y)
ARIMA_est <- list()
ic_arima <- matrix( nrow = (2 * 2 + 2) * 4 ^ 2, ncol = 7 )
colnames(ic_arima) <- c("d", "cons", "trend", "p", "q", "aic", "bic")
i <- 0
for (d in 0:1)
{
  for (const in 0:1)
  {
    for (p in 0:3)
    {
      for (q in 0:3)
      {
        i <- i + 1
        d1 <- as.logical(d)
        c1 <- as.logical(const)
        
        try(silent = T, expr =
        {
        ARIMA_est[[i]] <- Arima(y, order = c(p, d, q),
                                include.constant = c1)
        
        ic_arima[i,] <- c(d, const, 0, p, q,
                          ARIMA_est[[i]]$aic,
                          ARIMA_est[[i]]$bic)
        })

        if (const)
        {
          # only add a specification with trend if there is a
          # constant (i.e., exclude no constant with trend)
          i <- i + 1
          
          if (d1)
          {
            x <- c(0,cumsum(1:(TT - 1)))
          }
          else
          {
            x <- NULL
          }

          try(silent = T, expr =
          {
          ARIMA_est[[i]] <- Arima(y, order = c(p, d, q),
                                  xreg = x,
                                  include.constant = c1,
                                  include.drift = T)
          
          ic_arima[i,] <- c(d, const, 1, p, q,
                            ARIMA_est[[i]]$aic,
                            ARIMA_est[[i]]$bic)
          })
        }
      }
    }
  }
}

ic_aic_arima <- ic_arima[order(ic_arima[,6]),][1:10,]
ic_bic_arima <- ic_arima[order(ic_arima[,7]),][1:10,]

ic_aic_arima
ic_bic_arima

# find the intersection of AIC and BIC preferred sets
ic_int_arima <- intersect(as.data.frame(ic_aic_arima),
                          as.data.frame(ic_bic_arima))

ic_int_arima


# Note that the intersection contains only specifications
# in levels (d=0); take a closer look to see if any
# specifications in differences (d=1) are worth considering.
# Show the aic and bic preferred tables and point out that
# that the AIC heavily prefers d=0 specifications, whereas
# the BIC in fact slightly prefers d=1 specifications,
# although the ranking is more balanced.
# ARIMA(1,1,0) and ARIMA(2,1,0) -- both with a constant
# only -- are in the top 3 of the BIC ranking.
# Taking into consideration also our inference that
# that GDP is not empirically distinguishable from an I(1)
# process, we will add ARIMA(1,1,0) and ARIMA(2,1,0)
# to the four in the intersecting set.

rbind(ic_int_arima, ic_bic_arima[c(1, 3),])

as.data.frame(rbind(ic_int_arima, ic_bic_arima[c(1, 3),]))


adq_set_arima <- as.matrix(arrange(as.data.frame(
                           rbind(ic_int_arima,
                                 ic_bic_arima[c(1, 3),])),
                                      d, const, trend, p))
adq_idx_arima <- match(data.frame(t(adq_set_arima[, 1:5])),
                       data.frame(t(ic_arima[, 1:5])))


adq_set_arima

# Check the residuals for specifications in the adequate set.
for (i in 1:length(adq_idx_arima)){
  checkresiduals(ARIMA_est[[adq_idx_arima[i]]])
}


# 2
#
plot(dates, r, type = "l", xlab = "Time (Quarters)",
     main = "Federal Funds Rate")

TT <- length(r)
ADF_est <- list()
ic <- matrix( nrow = 30, ncol = 5 )
colnames(ic) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est[[i]] <- Arima(diff(r), xreg = r[-TT],
                          order = c(p, 0, 0),
                          include.mean = as.logical(const),
                          include.drift = F)
    ic[i,] <- c(const, 0, p, ADF_est[[i]]$aic,
                ADF_est[[i]]$bic)
  }
  
  if (const)
  {
    # only add a specification with trend if there is a
    # constant (i.e., exclude no constant with trend)
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est[[i]] <- Arima(diff(r), xreg = r[-TT],
                            order = c(p, 0, 0),
                            include.mean = as.logical(const),
                            include.drift = T)
      ic[i,] <- c(const, 1, p, ADF_est[[i]]$aic,
                  ADF_est[[i]]$bic)
    }
  }
}

ic_aic <- ic[order(ic[,4]),][1:10,]
ic_bic <- ic[order(ic[,5]),][1:10,]

# find the intersection of AIC and BIC preferred sets
ic_int <- intersect(as.data.frame(ic_aic),
                    as.data.frame(ic_bic))

ic_int

# intersecting set has 3 models: 2 with const and p = 6, 7,
# and 1 without constant and 7 lags; no models have a trend

adq_set <- as.matrix(arrange(as.data.frame(ic_int),
                             const, trend, p))
adq_idx <- match(data.frame(t(adq_set[, 1:3])),
                 data.frame(t(ic[, 1:3])))

for (i in 1:length(adq_idx))
{
  checkresiduals(ADF_est[[adq_idx[i]]])
}

# Residuals look ok, although some ACs exceed 95% conf
# intervals very slightly (and at larger lags). For robustness,
# we might add a few larger p models into the adequate set,
# such as p = 8, 9.

# Here the default option p <= 5 is too small, so we need to
# explicitly specify longer lag lengths.

ic_int
adf.test(r, nlag = 10)

# Only "Type 1" and "Type 2" specifications are in our adequate set,
# so we ignore the output related to "Type 3".
# For all specifications in our adequate set, null of unit root
# cannot be rejected. Note however, that there are specifications
# for which the null is indeed rejected, such as constant and
# no trend with p = 5. This is the 9th ranked specification by
# the BIC and outside the top 10 in terms of AIC.
# We can summarise our inference as follows.
# The process is not empirically distinguishable from a unit root
# process with a drift, but with a small element of uncertainty.

TT <- length(diff(r))
ADF_est_diff <- list()
ic_diff <- matrix( nrow = 30, ncol = 5 )
colnames(ic_diff) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est_diff[[i]] <- Arima(diff(diff(r)),
                               xreg = diff(r)[-TT],
                               order = c(p, 0, 0),
                               include.mean = as.logical(const),
                               include.drift = F)
    ic_diff[i,] <- c(const, 0, p, ADF_est_diff[[i]]$aic,
                     ADF_est_diff[[i]]$bic)
  }
  
  if (const)
  {
    # only add a specification with trend if there is a
    # constant (i.e., exclude no constant with trend)
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est_diff[[i]] <- Arima(diff(diff(r)),
                                 xreg = diff(r)[-TT],
                                 order = c(p, 0, 0),
                                 include.mean = as.logical(const),
                                 include.drift = T)
      ic_diff[i,] <- c(const, 1, p, ADF_est_diff[[i]]$aic,
                       ADF_est_diff[[i]]$bic)
    }
  }
}

ic_aic_diff <- ic_diff[order(ic_diff[,4]),][1:10,]
ic_bic_diff <- ic_diff[order(ic_diff[,5]),][1:10,]

# find the intersection of AIC and BIC preferred sets
ic_int_diff <- intersect(as.data.frame(ic_aic_diff),
                         as.data.frame(ic_bic_diff))

ic_int_diff

# There are only two specifications in the intersecting
# set: both with no constant or trend and p = 6, 7.
# Overall, both AIC and BIC favour specifications
# without a trend, and AIC generally favours larger
# lag lengths, whereas the BIC yields a more dispersed
# ranking.

adq_set_diff <- as.matrix(arrange(
                          as.data.frame(ic_int_diff),
                                       const, trend, p))
adq_idx_diff <- match(data.frame(t(adq_set_diff[, 1:3])),
                      data.frame(t(ic_diff[, 1:3])))

for (i in 1:length(adq_idx_diff))
{
  checkresiduals(ADF_est_diff[[adq_idx_diff[i]]])
}

# As with levels, residuals look ok, but with some indication
# that longer lag lengths should be considered.

adf.test(diff(r), nlag = 10)

# unit root is rejected at very small significance level
# for all specifications.
# We infer that FFR is not empirically distinguishable from an
# I(1) process.