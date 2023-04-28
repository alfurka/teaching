library(forecast)
library(dplyr)
library(rugarch)

# load the data in Merck.csv
mydata <- read.delim("Merck.csv", header = TRUE,  sep = ",")

date <- as.Date(mydata$Date, format = "%d/%m/%Y")
y <- mydata$y
r <- diff(log(y))

# plot the returns
plot(date[-1], r, type = "l", xlab = "", ylab = "returns")

# 1
#
# (a)
#
# same approach as with basic GARCH but now add also the threshold option;
# note: what we call "TGARCH", the "rugarch" package calls "gjrGARCH" and it
# calls "TGARCH" a different variant -- this is common, so we must always check
# the documentation to be sure!
submods <- c("sGARCH", "gjrGARCH")
ARMA_TGARCH_est <- list()
ic_arma_tgarch <- matrix( nrow = 4 ^ 2 * 2 ^ 3,
                          ncol = 7 )
colnames(ic_arma_tgarch) <- c("pm", "qm", "ph", "qh",
                              "thresh", "aic", "bic")
i <- 0; t0 <- proc.time()
for (pm in 0:3)
{
  for (qm in 0:3)
  {
    for (ph in 0:1)
    {
      for (qh in 0:1)
      {
        for (thresh in 0:1)
        {
          i <- i + 1
          ic_arma_tgarch[i, 1:5] <- c(pm, qm, ph, qh, thresh)
          
          if (ph == 0 && qh == 0)
          {
            if (!thresh) # no such thing as a homoscedastic threshold ARMA
            {
              # for models with constant variance, the ugarchspec and
              # ugarchfit functions do not work well; instead, the
              # documentation advises to use arfimaspec and arfimafit
              try(silent = T, expr =
              {
                ARMA_TGARCH_mod <- arfimaspec(
                  mean.model = list(armaOrder = c(pm, qm)))
                
                ARMA_TGARCH_est[[i]] <- arfimafit(ARMA_TGARCH_mod, r)
                
                ic_arma_tgarch[i,6:7] <- infocriteria(
                  ARMA_TGARCH_est[[i]])[1:2]
              })
            }
          }
          else
          {
            try(silent = T, expr =
            {
              ARMA_TGARCH_mod <- ugarchspec(
                mean.model = list(armaOrder = c(pm, qm)),
                variance.model = list(model = submods[1 + thresh],
                                      garchOrder = c(ph, qh)))
              
              ARMA_TGARCH_est[[i]] <- ugarchfit(ARMA_TGARCH_mod, r,
                                                solver = 'hybrid')
              
              ic_arma_tgarch[i,6:7] <- infocriteria(ARMA_TGARCH_est[[i]])[1:2]
            })
          }
          cat("\r", "Processed ARMA(",
              pm, ",", qm, ")-",
              c("", "T")[1 + thresh], "GARCH(",
              ph, ",", qh, ")",
              c(" ", "")[1 + thresh], sep = "")
          
        }
      }
    }
  }
}
cat("\n", proc.time() - t0, "\n")

ic_aic_arma_tgarch <- ic_arma_tgarch[order(ic_arma_tgarch[,6]),][1:20,]
ic_bic_arma_tgarch <- ic_arma_tgarch[order(ic_arma_tgarch[,7]),][1:20,]

ic_aic_arma_tgarch
ic_bic_arma_tgarch

# find the intersection of AIC and BIC preferred sets
ic_int_arma_tgarch <- intersect(as.data.frame(ic_aic_arma_tgarch),
                                as.data.frame(ic_bic_arma_tgarch))

ic_int_arma_tgarch

# We select the entire intersection set.
adq_set_arma_tgarch <- as.matrix(arrange(as.data.frame(
                                  ic_int_arma_tgarch), pm, qm, ph, qh, thresh))
adq_idx_arma_tgarch <- match(data.frame(t(adq_set_arma_tgarch[, 1:5])),
                                  data.frame(t(ic_arma_tgarch[, 1:5])))

# Check the residuals
nmods <- length(adq_idx_arma_tgarch)
sacf_tgarch <- matrix(nrow = nmods, ncol = 15)
colnames(sacf_tgarch) <- c("pm", "qm", "ph", "qh", "thresh", 1:10)
for (i in 1:nmods)
{
  sacf_tgarch[i,1:5] <- adq_set_arma_tgarch[i,1:5]
  sacf_tgarch[i,6:15] <-
                  acf(ARMA_TGARCH_est[[adq_idx_arma_tgarch[i]]]@fit$z,
                                       lag = 10, plot = F)$acf[2:11]
}

sacf_tgarch

# Autocorrelations appear to be relatively small for all models.

# (b)
#
# plotting estimated volatilities is the same as for basic GARCH (we still
# don't get confidence intervals)
title_tgarch <- rep(NA, times = nmods)
for (i in 1:nmods)
{
  title_tgarch[i] <- paste("ARMA(",
                     as.character(adq_set_arma_tgarch[i, 1]), ",",
                     as.character(adq_set_arma_tgarch[i, 2]),
                     ")-",
                     c("", "T")[1 + adq_set_arma_tgarch[i,5]],
                     "GARCH(",
                     as.character(adq_set_arma_tgarch[i, 3]), ",",
                     as.character(adq_set_arma_tgarch[i, 4]), ")",
                     sep = "")
  plot(date[-1], sqrt(
        ARMA_TGARCH_est[[adq_idx_arma_tgarch[i]]]@fit$var),
        type = "l", xlab = "", ylab = "volatilities",
        ylim = c(0, 0.08), main = title_tgarch[i])
}

# All volatility estimates generally look very similar.
# TGARCH vs GARCH does not make a noticeable difference
# for fixed pm, qm, ph and qh. However, larger pm and qm
# appear to be associated with larger "spikes" in volatilities
# when they do occur.

# (c)
#
# To test for leverage effects, we focus on on the threshold coefficient
# in models where the threshold is actually included. Note that rugarch calls
# the threshold parameter gamma, whereas we call it lambda.
for (i in 1:nmods)
{
  if (adq_set_arma_tgarch[i, 5] == 1)
  {
    # this is a specification with a threshold
    lambda_est <- ARMA_TGARCH_est[[
              adq_idx_arma_tgarch[i]]]@fit$coef["gamma1"]
    lambda_tvl <- ARMA_TGARCH_est[[
              adq_idx_arma_tgarch[i]]]@fit$tval["gamma1"]
    cat(paste0(title_tgarch[i], ": lambda_hat = ",
                                   round(lambda_est, 2),
                                ", t-value = ",
                                   round(lambda_tvl, 2)),
                                "\n")
  }
}
# the null of "no leverage effects" is easily rejected in favour of
# "leverage effects" for all models except ARMA(2,1)-TGARCH(1,1), where
# "reverse effects" are being confirmed at a very low significance level

# (d)
#
for (i in 1:nmods)
{
  plot(ugarchboot(ARMA_TGARCH_est[[adq_idx_arma_tgarch[i]]],
                  method = "Partial"), which = 3)
}

# There are subtle differences between the volatility forecasts generated by
# different specifications, but all agree that volatilities across the forecast
# horizon will remain well below the levels estimated around 2005 and then again
# in 2008/2009

# 2
#
# (a)
#
# same approach as with GARCH/TGARCH models; here it is somewhat simpler
# because we fix all GARCH lags to be 1.
ARMA_GARCHM_est <- list()
ic_arma_garchm <- matrix( nrow = 4 ^ 2 * 2, ncol = 5 )
colnames(ic_arma_garchm) <- c("p", "q", "m", "aic", "bic")
i <- 0; t0 <- proc.time()
for (p in 0:3)
{
  for (q in 0:3)
  {
    for (m in 0:1)
    {
      i <- i + 1
      ic_arma_garchm[i, 1:3] <- c(p, q, m)
        
      try(silent = T, expr =
      {
        ARMA_GARCHM_mod <- ugarchspec(
                            mean.model = list(armaOrder = c(p, q),
                                              archm = m),
                            variance.model = list(garchOrder = c(1, 1)))
                    
        ARMA_GARCHM_est[[i]] <- ugarchfit(ARMA_GARCHM_mod, r,
                                          solver = 'hybrid')
                  
        ic_arma_garchm[i,4:5] <- infocriteria(ARMA_GARCHM_est[[i]])[1:2]
      })
      cat("\r", "Processed ARMA(",p, ",", q, ")",
              c("-GARCH(1,1) ", "-GARCHM(1,1)")[1 + m], sep = "")
          
    }
  }
}
cat("\n", proc.time() - t0, "\n")

ic_aic_arma_garchm <- ic_arma_garchm[order(ic_arma_garchm[,4]),][1:20,]
ic_bic_arma_garchm <- ic_arma_garchm[order(ic_arma_garchm[,5]),][1:20,]

ic_aic_arma_garchm
ic_bic_arma_garchm

# find the intersection of AIC and BIC preferred sets
ic_int_arma_garchm <- intersect(as.data.frame(ic_aic_arma_garchm),
                                as.data.frame(ic_bic_arma_garchm))

# We select the entire intersection set.
adq_set_arma_garchm <- as.matrix(arrange(as.data.frame(
                                ic_int_arma_garchm), p, q, m))
adq_idx_arma_garchm <- match(data.frame(t(adq_set_arma_garchm[, 1:3])),
                                  data.frame(t(ic_arma_garchm[, 1:3])))

# Check the residuals
nmods <- length(adq_idx_arma_garchm)
sacf_garchm <- matrix(nrow = nmods, ncol = 13)
colnames(sacf_garchm) <- c("p", "q", "m", 1:10)
for (i in 1:nmods)
{
  sacf_garchm[i,1:3] <- adq_set_arma_garchm[i,1:3]
  sacf_garchm[i,4:13] <-
                  acf(ARMA_GARCHM_est[[adq_idx_arma_garchm[i]]]@fit$z,
                                       lag = 10, plot = F)$acf[2:11]
}

sacf_garchm

# Autocorrelations appear to be relatively small for all models.

# (b)
#
title_garchm <- rep(NA, times = nmods)
for (i in 1:nmods)
{
  title_garchm[i] <- paste("ARMA(",
        as.character(adq_set_arma_garchm[i, 1]), ",",
        as.character(adq_set_arma_garchm[i, 2]),
        c(")-GARCH(1,1)", ")-GARCHM(1,1)")[1 + adq_set_arma_garchm[i,3]],
        sep = "")
  plot(date[-1], sqrt(
        ARMA_GARCHM_est[[adq_idx_arma_garchm[i]]]@fit$var),
        type = "l", xlab = "", ylab = "volatilities",
        ylim = c(0, 0.08), main = title_garchm[i])
}

# All volatility estimates generally look very similar.
# GARCHM vs GARCH appears to exhibit only a minor difference
# for fixed p, q. The periods of high volatility match up to those estimated
# using TGARCH models.

# (c)
#
# To test for time-varying risk premia, we focus on on the "archm" coefficient
# in models where the GARCH-in-mean is actually included (this is "delta" in
# Engle, et al., 1987).
for (i in 1:nmods)
{
  if (adq_set_arma_garchm[i, 3] == 1)
  {
    # this is a specification with GARCH-in-mean
    archm_est <- ARMA_GARCHM_est[[
                          adq_idx_arma_garchm[i]]]@fit$coef["archm"]
    archm_tvl <- ARMA_GARCHM_est[[
                          adq_idx_arma_garchm[i]]]@fit$tval["archm"]
    cat(paste0(title_garchm[i], ": archm_hat = ",
                                   round(archm_est, 2),
                                ", t-value = ",
                                   round(archm_tvl, 2)),
                                "\n")
  }
}
# the null of "time-invariant risk premium" is easily rejected in favour of
# "time-varying risk premium" for all models.

# (d)
#
for (i in 1:nmods)
{
  plot(ugarchboot(ARMA_GARCHM_est[[adq_idx_arma_garchm[i]]],
                                   method = "Partial"), which = 3)
}

# volatility forecasts are very similar to those obtained with the TGARCH.