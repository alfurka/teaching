# Load necessary packages
library(forecast)
library(dplyr)
library(rugarch)

# load the data in cwb.csv
mydata <- read.delim("cwb.csv", header = TRUE,  sep = ",")

date <- as.Date(mydata$date, format = "%m/%d/%Y")
y <- mydata$y

# 1
# Plot the returns
plot(date, log(y), type = "l", xlab = "Time", ylab = "share prices")

# 2
# Generate the log returns of y
r <- diff(log(y))
# Plot the returns
plot(date[-1], r, type = "l", xlab = "", ylab = "returns")

# 3
# Select an adequate set of ARMA models
ARMA_est <- list()
ic_arma <- matrix(nrow = 4 * 4, ncol = 4 )
colnames(ic_arma) <- c("p", "q", "aic", "bic")
for (p in 0:3)
{
  for (q in 0:3)
  {
    i <- p * 4 + q + 1
    ARMA_est[[i]] <- Arima(r, order = c(p, 0, q))
    ic_arma[i,] <- c(p, q, ARMA_est[[i]]$aic, ARMA_est[[i]]$bic)
  }
}

ic_aic_arma <- ic_arma[order(ic_arma[,3]),][1:10,]
ic_bic_arma <- ic_arma[order(ic_arma[,4]),][1:10,]

# find the intersection of AIC and BIC preferred sets
ic_int_arma <- intersect(as.data.frame(ic_aic_arma),
                         as.data.frame(ic_bic_arma))

# select MA(1), MA(2), AR(1), AR(2), ARMA(1,1)
adq_set_arma <- as.matrix(arrange(as.data.frame(
                                 rbind(ic_int_arma[c(1:3, 6),],
                                  ic_bic_arma[2,])), p, q))
adq_idx_arma <- match(data.frame(t(adq_set_arma[, 1:2])),
                           data.frame(t(ic_arma[, 1:2])))

# Check the residuals for specifications in the adequate set.
nmods <- length(adq_idx_arma)
for (i in 1:nmods){
  checkresiduals(ARMA_est[[adq_idx_arma[i]]])
}

# 4
# Generate squared residuals and plot them.
e2_arma <- list()
for (i in 1:nmods)
{
  e2_arma[[i]] <- resid(ARMA_est[[adq_idx_arma[i]]]) ^ 2
  
  title_p_q <- paste("ARMA(",
                     as.character(adq_set_arma[i, 1]), ", ",
                     as.character(adq_set_arma[i, 2]), ")",
                     sep = "")
  plot(date[-1], e2_arma[[i]], type = "l",
                 xlab = "", ylab = "squared resid",
                 main = paste("Plot: ", title_p_q))
  
  acf(e2_arma[[i]], xlab = "", ylab = "",
                 main = paste("SACF: ", title_p_q))
}

# squared residuals are clearly autocorrelated

# 5
#
# Breusch-Pagan type (LM) test for presence of heteroscedasticity
bptest <- matrix(nrow = 10 * nmods, ncol = 5)
colnames(bptest) <- c("p", "q", "j", "LM-stat", "p-value")
for (i in 1:nmods){
  e2_i <- as.vector(e2_arma[[i]])
  f <- formula(e2_i ~ 1)
  for (j in 1:10)
  {
    # lag lengths in the auto-regression of squared residuals
    k <- 10 * (i - 1) + j
    f <- update.formula(f, paste("~ . + lag(e2_i, n = ", j, ")"));
    bp_reg_j <- lm(f)
    LM_j <- length(e2_i) * summary(bp_reg_j)$r.squared
    p_val_j <- 1 - pchisq(LM_j, df = j)
    bptest[k,] <- c(adq_set_arma[i, 1:2], j, LM_j, p_val_j)
  }
}

bptest

# 6
ARMA_GARCH_est <- list()
ic_arma_garch <- matrix( nrow = 3 ^ 4, ncol = 6 )
colnames(ic_arma_garch) <- c("pm", "qm", "ph", "qh", "aic", "bic")
i <- 0
for (pm in 0:2)
{
  for (qm in 0:2)
  {
    for (ph in 0:2)
    {
      for (qh in 0:2)
      {
        i <- i + 1
        ic_arma_garch[i, 1:4] <- c(pm, qm, ph, qh)

        if (ph == 0 && qh == 0)
        {
          # for models with constant variance, the ugarchspec and
          # ugarchfit functions do not work well; instead, the
          # documentation advises to use arfimaspec and arfimafit
          ARMA_GARCH_mod <- arfimaspec(
            mean.model = list(armaOrder = c(pm, qm)))
          
          ARMA_GARCH_est[[i]] <- arfimafit(ARMA_GARCH_mod, r)
          
          ic_arma_garch[i,5:6] <- infocriteria(ARMA_GARCH_est[[i]])[1:2]
        }
        else
        {
          try(silent = T, expr =
          {
            ARMA_GARCH_mod <- ugarchspec(
                      mean.model = list(armaOrder = c(pm, qm)),
                      variance.model = list(garchOrder = c(ph, qh)))
            
            ARMA_GARCH_est[[i]] <- ugarchfit(ARMA_GARCH_mod, r,
                                             solver = 'hybrid')
            
            ic_arma_garch[i,5:6] <- infocriteria(ARMA_GARCH_est[[i]])[1:2]
          })
        }
      }
    }
  }
}

ic_aic_arma_garch <- ic_arma_garch[order(ic_arma_garch[,5]),][1:40,]
ic_bic_arma_garch <- ic_arma_garch[order(ic_arma_garch[,6]),][1:40,]

ic_int_arma_garch <- intersect(as.data.frame(ic_aic_arma_garch),
                               as.data.frame(ic_bic_arma_garch))

adq_set_arma_garch <- as.matrix(arrange(as.data.frame(
                      ic_int_arma_garch[1:36,]), pm, qm, ph, qh))
adq_idx_arma_garch <- match(data.frame(t(adq_set_arma_garch[, 1:4])),
                            data.frame(t(ic_arma_garch[, 1:4])))

# Check the residuals 
nmods <- length(adq_idx_arma_garch)
sacf_garch <- matrix(nrow = nmods, ncol = 14)
colnames(sacf_garch) <- c("pm", "qm", "ph", "qh", 1:10)
for (i in 1:nmods)
{
  sacf_garch[i,1:4] <- adq_set_arma_garch[i,1:4]
  sacf_garch[i,5:14] <-
                  acf(ARMA_GARCH_est[[adq_idx_arma_garch[i]]]@fit$z,
                  lag = 10, plot = F)$acf[2:11]
}

sacf_garch
# Autocorrelations appear to be relatively small for all models.

# 7

#example

ARMA_GARCH_est[[adq_idx_arma_garch[1]]]@fit$var

# plot

for(i in 1:nmods){
  title_p_q <- paste("ARMA(",
                     as.character(adq_set_arma_garch[i, 1]), ", ",
                     as.character(adq_set_arma_garch[i, 2]),
                     ")-GARCH(",
                     as.character(adq_set_arma_garch[i, 3]), ", ",
                     as.character(adq_set_arma_garch[i, 4]), ")",
                     sep = "")
  plot(date[-1], ARMA_GARCH_est[[adq_idx_arma_garch[i]]]@fit$var,
                     type = "l", xlab = "", ylab = "volatilities",
                     main = title_p_q)
}
