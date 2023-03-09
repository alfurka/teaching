rm(list = ls())
getwd()

# for this tutorial we will need the package "forecast"; show students how to
# install and be sure to load it. Note that "forecast" requires the package
# "colorspace", which is not (always) installed automatically because on CRAN,
# the source version is newer than the binary version. To install "colorspace",
# it may be necessary to begin installation then say "No" to installing from
# sources.
library(forecast)

# 1
#
# (a)
#
# load the data in Merck.csv
mydata <- read.delim("Merck.csv", header = TRUE,  sep = ",")

sel_sample <- mydata$Date >= as.Date("2011-01-01") &
              mydata$Date <= as.Date("2012-01-31")
y <- as.matrix(mydata$Adj_Close[sel_sample])

# (b)
#
# difference and lag operators
Dy <- diff(y)
r <- as.matrix(log(y[2:nrow(y)]) - log(lag(y[1:nrow(y) - 1])))

colnames(y) <- "Stock Prices of Merck (MRK)"
colnames(Dy) <- "Changes in Stock Prices"
colnames(r) <- "Log Returns" 

# (c)
#
# we use time series plot to detect obvious trends, volatility clustering, etc
# this is typically the first thing to do when processing time series data
sel2011 <- mydata$Date[sel_sample] <= as.Date("2011-12-31")
dates = as.Date(mydata$Date[sel_sample])
plot(dates[sel2011], y[sel2011], type = "l", xlab = "Time (2011)",
     ylab = colnames(y))
plot(dates[sel2011], Dy[sel2011], type = "l", xlab = "Time (2011)",
     ylab = colnames(Dy))
plot(dates[sel2011], r[sel2011], type = "l", xlab = "Time (2011)",
     ylab = colnames(r))
# the DGP for y is likely not stationary as its mean appears to vary over time;
# the DGP for Dy seems to have zero mean, but its variance may depend on time
# ignore the variance for now (we will cover volatility later) 

# (d)
#
# For prices, SACF decays but the SPACF drops to zero after one lag;
# also, the first PAC is near 1, suggesting an AR(1) model of the form
# y[t] = y[t-1] + e[t]. Point out to students the connection to Dy[t].
#
# For price changes, SPACF drops to zero after one lag but the SACF oscillates;
# also, the first AC is near 1, suggesting an MA(1) model of the form
# Dy[t] = e[t-1] + e[t]. Challenge the students think about what's going on.
acf(y[sel2011], main = colnames(y))
pacf(y[sel2011], main = colnames(y))
acf(Dy[sel2011], main = colnames(Dy))
pacf(Dy[sel2011], main = colnames(Dy))

# (e)
#
# explain students how to construct a nice table of aic and bic values for
# combinations of p and q values;
# point out that running estimation in nested for loops is easy to program
# and can be executed very quickly!
ic <- matrix( nrow = 25, ncol = 4 )
colnames(ic) <- c("p", "q", "aic", "bic")
for (p in 0:4)
{
  for (q in 0:4)
  {
    fit_p_q <- Arima(Dy, order = c(p, 0, q))
    ic[p * 5 + q + 1,] = c(p, q, fit_p_q[["aic"]], fit_p_q[["bic"]])
  }
}

# (f)
#
# a quick look through the aic/bic table constructed in (g) shows that aic and
# bic wildly disagree!
# there is no systematic approach to resolving this conflicting information;
# we will look at the top 10 specifications preferred by the aic as well as
# the top 10 preferred by the bic---this is easy to do by sorting the ic table
print(ic)
ic_aic <- ic[order(ic[,3], decreasing = FALSE),][1:10,]
ic_bic <- ic[order(ic[,4], decreasing = FALSE),][1:10,]
print(ic_aic)
print(ic_bic)

# next, we select the models that make top 10 in both lists;
# however, reinforce the point that this is just a "sensible" choice in this
# particular setting -- there may be other reasonable ways of using the info
# provided by aic/bic to reduce the set;
# here is where judgement / experience / common sense come into play!
adq_set = list(c(1, 0, 1), c(1, 0, 2), c(2, 0, 1), c(3, 0, 0))

# (g)
#
# once we have reduced the set, perform residual analysis;
# the "forecast" package simplifies this with the command "checkresiduals"
# we don't get carried away with any particular hypothesis test here --
# nothing clearly stands out to suggest a problem with correlated residuals,
# so we proceed with the four models in the set
for (i in 1:length(adq_set))
{
  checkresiduals(Arima(Dy[sel2011], order = adq_set[[i]]))
}

# (h)
#
# here, we forecast each model in the set using the "forecast" command;
# again, it is all one easily and quickly with a for loop;
# point out that in the forecast command, we pass an option "level = c(68, 95)",
# which instructs the command to construct two sets of predictive intervals:
# one with 95% coverage and another with 68% coverage;
# discuss with students why we would want to compare predictive intervals
# at different coverage levels, by connecting it to "risk" associated with
# decisions based on forecasts
hrz = sum(sel_sample) - sum(sel2011)
xticks <- c(sum(sel_sample) - 3 * hrz + c(1, 2 * hrz, 3 * hrz))
actual_Dy <- as.matrix(Dy[!sel2011])
fcst_Dy <- vector(mode = "list", length(adq_set))
for (i in 1:length(adq_set))
{
  model_p_q <- adq_set[[i]]
  fcst_Dy[[i]] <- forecast(Arima(Dy[sel2011], order = model_p_q),
                           h = hrz, level = c(68, 95))

  title_p_q <- paste("ARMA(", as.character(model_p_q[1]), ", ",
                              as.character(model_p_q[3]), ")", sep = "")
  
  plot(fcst_Dy[[i]], include = hrz * 2, ylab = colnames(Dy),
       main = title_p_q, xaxt = "n")
  lines(sum(sel2011) + 1:hrz, actual_Dy)
  axis(1, at = xticks, labels = dates[xticks])
}
# when we use the "plot" command with output from the "forecast" command, we
# get a nice depiction of how the data is extrapolated into the future,
# complete with predictive intervals to capture uncertainty;
# we also add the actual outcomes in the forecast period to help us compare
# the forecast performance of each ARMA in the adequate set;
#
# point out that predictive intervals for price changes appear to have a fixed
# width even as the forecast horizon increases (from 1 day to 20 days)
#
# comparison and inference:
# clearly, all four generate very similar forecasts for January 2012; therefore,
# the specification differences between them are not important for our purpose;
# alternatively, our forecast of price changes is ROBUST to minor differences in
# the specification of ARMA models in that we cannot clearly distinguish between
# them with our diagnostic tools.

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

# now repeat (g) using ARIMA on "y" instead of ARMA on "Dy";
# explain to students that when we use an ARIMA on "y", it is the same as
# estimating an ARMA for "Dy", but the algorithm will generate forecasts for "y"
# instead of "Dy"; in particular, it generates forecasts for "Dy", but then
# cumulatively sums them starting with an initial point "y0", which in our case
# is the last observation in the "pre-sample";
# DO NOT go into a discussion on "integration" or "unit roots" -- just let the
# students know that this is part of a broader methodology that will be
# discussed in Week 5; the important point here is that we can forecast "y"
# using either an ARIMA or an ARMA.
y0 <- mydata$Adj_Close[sum(mydata$Date < as.Date("2011-01-01") - 1)]
y_ext = as.matrix(c(y0, y[sel2011]))
fcst_y <- vector(mode = "list", length(adq_set))
for (i in 1:length(adq_set))
{
  model_p_q <- adq_set[[i]]
  model_p_q[2] = 1
  fcst_y[[i]] <- forecast(Arima(y_ext, order = model_p_q, include.constant = T),
                          h = hrz, level = c(68, 95) )
  
  title_p_q <- paste("ARIMA(", as.character(model_p_q[1]), ", ",
                     as.character(model_p_q[3]), ")", sep = "")
  
  plot(fcst_y[[i]], include = hrz * 2, ylab = colnames(y),
       main = title_p_q, xaxt = "n")
  lines(1 + sum(sel2011) + 1:hrz, actual_y)
  axis(1, at = 1 + xticks, labels = dates[xticks])
}
# the first thing to note from the results is that, again, all ARIMAs in the
# adequate set generate very similar forecasts;
# the second is to note that the predictive intervals for "y" increase as the
# horizon increases (they are narrower for forecasts in the beginning of the
# forecasting period, and wider towards the end of the forecast period);
# this is a reflection of "y" being possibly generated by a non-stationary
# process --> more uncertainty in the more distant future!
# finally, compare the ARIMA(1,1) to the ARMA(2,1) in "y"; just hint that we
# will cover later why in particular these two specifications provide an
# especially interesting comparison;
# for now, point out that the ARMA(2,1) also produces predictive intevals that
# increase as the horizon increases, but relative to the ARIMA(1,1), the
# intervals indicate that prices should fall in January 2012. This is slightly
# different from what the ARIMA(1,1) produces -- although the intervals largely
# overlap, the ARIMA(1,1) clearly puts more weight on higher prices in January;
# when comparing to the actual observations in January 2012, it is easy to see
# that the ARIMA forecasts are better (which can be confirmed by formal metrics)

# (j)
#
# only do this if there is time remaining; if there is, let students try this
# on their own and provide support as needed
acf(r[sel2011], main = colnames(y))
pacf(r[sel2011], main = colnames(y))

ic <- matrix( nrow = 25, ncol = 4 )
colnames(ic) <- c("p", "q", "aic", "bic")
for (p in 0:4)
{
  for (q in 0:4)
  {
    fit_p_q <- Arima(r, order = c(p, 0, q))
    c(p * 5 + q + 1, p, q)
    ic[p * 5 + q + 1,] = c(p, q, fit_p_q[["aic"]], fit_p_q[["bic"]])
  }
}

ic_aic <- ic[order(ic[,3], decreasing = FALSE),][1:10,]
ic_bic <- ic[order(ic[,4], decreasing = FALSE),][1:10,]

adq_set = list(c(1, 0, 1), c(1, 0, 2), c(2, 0, 1), c(3, 0, 0))

for (i in 1:length(adq_set))
{
  checkresiduals(Arima(r[sel2011], order = adq_set[[i]]))
}

hrz <- sum(sel_sample) - sum(sel2011)
xticks <- c(sum(sel_sample) - 3 * hrz + c(1, 2 * hrz, 3 * hrz))
actual_r <- as.matrix(r[!sel2011])
fcst_r <- vector(mode = "list", length(adq_set))
for (i in 1:length(adq_set))
{
  model_p_q <- adq_set[[i]]
  fcst_r[[i]] <- forecast(Arima(r[sel2011], order = model_p_q),
                          h = hrz, level = c(68, 95))
  
  title_p_q <- paste("ARMA(", as.character(model_p_q[1]), ", ",
                     as.character(model_p_q[3]), ")", sep = "")
  
  plot(fcst_r[[i]], include = hrz * 2, ylab = colnames(r),
       main = title_p_q, xaxt = "n")
  lines(sum(sel2011) + 1:hrz, actual_r)
  axis(1, at = xticks, labels = dates[xticks])
}
