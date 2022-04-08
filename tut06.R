library(forecast)
library(dplyr)
library(zoo)
library(aTSA)


ADF_estimate_lev <- function(y,  p_max = 9){
  TT <- length(y)
  ADF_est <- list()
  ic <- matrix(nrow = 3 * (1 + p_max), ncol = 5)
  colnames(ic) <- c("const", "trend", "p", "aic", "bic")
  i <- 0
  for (const in 0:1)
  {
    for (p in 0:p_max)
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
      for (p in 0:p_max)
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
  
  return(list(ADF_est = ADF_est, ic = ic,
              ic_aic = ic_aic, ic_bic = ic_bic))
}


ADF_estimate_diff <- function(y, p_max = 9){
  TT <- length(diff(y))
  ADF_est_diff <- list()
  ic_diff <- matrix(nrow = 3 * (1 + p_max), ncol = 5)
  colnames(ic_diff) <- c("const", "trend", "p", "aic", "bic")
  i <- 0
  for (const in 0:1)
  {
    for (p in 0:p_max)
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
      for (p in 0:p_max)
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
  
  return(list(ADF_est_diff = ADF_est_diff,
              ic_diff = ic_diff,
              ic_aic_diff = ic_aic_diff,
              ic_bic_diff = ic_bic_diff))  
}

# 1.
#
# load the data in term_structure.csv
mydata <- read.delim("term_structure.csv", header = TRUE,  sep = ",")

dates <- as.yearqtr(mydata$obs)
i3y <- mydata$I3Y
i5y <- mydata$I5Y
i90d <- mydata$I90D
i180d <- mydata$I180D

# i3y
i3y_ADF_lev <- ADF_estimate_lev(i3y, p_max = 15)
print(i3y_ADF_lev$ic_aic)
print(i3y_ADF_lev$ic_bic)

# from AIC pick lower p
# best BIC + large p

i3y_adq_set <- as.matrix(arrange(as.data.frame(
                  rbind(i3y_ADF_lev$ic_aic[c(1, 6, 9),],
                        i3y_ADF_lev$ic_bic[c(1, 3, 5:7),])),
                        const, trend, p))
i3y_adq_idx <- match(data.frame(t(i3y_adq_set[, 1:3])),
                 data.frame(t(i3y_ADF_lev$ic[, 1:3])))

for (i in 1:length(i3y_adq_idx)){
  checkresiduals(i3y_ADF_lev$ADF_est[[i3y_adq_idx[i]]])
}
# We reject white noise residuals at 5% sig level for all models with p < 10.
# Hence, we remove all models except the with p = 10 and constant +/- trend.

adf.test(i3y, nlag = 15)

# now repeat for the differenced i3y
i3y_ADF_diff <- ADF_estimate_diff(i3y, p_max = 15)
print(i3y_ADF_diff$ic_aic_diff)
print(i3y_ADF_diff$ic_bic_diff)

# take Intersection
i3y_adq_set_diff <- as.matrix(arrange(as.data.frame(
                       i3y_ADF_diff$ic_bic_diff[c(1, 3, 5),]),
                       const, trend, p))
i3y_adq_idx_diff <- match(data.frame(
                       t(i3y_adq_set_diff[, 1:3])),
                       data.frame(
                       t(i3y_ADF_diff$ic_diff[, 1:3])))

for (i in 1:length(i3y_adq_idx_diff)){
  checkresiduals(
    i3y_ADF_diff$ADF_est_diff[[i3y_adq_idx_diff[i]]])
}
# All residuals look ok.

adf.test(diff(i3y))

# i5y
i5y_ADF_lev <- ADF_estimate_lev(i5y, p_max = 15)
print(i5y_ADF_lev$ic_aic)
print(i5y_ADF_lev$ic_bic)

# best of AIC and BIC + intersection

i5y_adq_set <- as.matrix(arrange(as.data.frame(
  rbind(i5y_ADF_lev$ic_aic[1,],
        i5y_ADF_lev$ic_bic[c(1, 7:8),])),
  const, trend, p))
i5y_adq_idx <- match(data.frame(t(i5y_adq_set[, 1:3])),
                     data.frame(t(i5y_ADF_lev$ic[, 1:3])))

for (i in 1:length(i5y_adq_idx)){
  checkresiduals(i5y_ADF_lev$ADF_est[[i5y_adq_idx[i]]])
}
# We reject white noise residuals at 5% sig level for all const + trend and
# p = 2, 3. Hence we are left with only no const, no trend, p = 1 and const
# + trend, p = 10. This suggests maybe we should consider models with no const,
# no trend and small p along with const + trend and larger p. We expand the set
# to include no const, no trend, p = 0-3 along with const + trend, p = 9-12.

print(i5y_ADF_lev$ic_aic)
print(i5y_ADF_lev$ic_bic)

# large p models from AIC
# no trend no constant models from BIC 


i5y_adq_set <- as.matrix(arrange(as.data.frame(
  rbind(i5y_ADF_lev$ic_aic[c(1, 6:8),],
        i5y_ADF_lev$ic_bic[c(1, 3, 6, 10),])),
  const, trend, p))
i5y_adq_idx <- match(data.frame(t(i5y_adq_set[, 1:3])),
                     data.frame(t(i5y_ADF_lev$ic[, 1:3])))

for (i in 1:length(i5y_adq_idx)){
  checkresiduals(i5y_ADF_lev$ADF_est[[i5y_adq_idx[i]]])
}


# To accommodate all specs in the adequate set, we use nlags = 12
adf.test(i5y, nlag = 12)


# now repeat for the differenced i5y
i5y_ADF_diff <- ADF_estimate_diff(i5y, p_max = 15)
print(i5y_ADF_diff$ic_aic_diff)
print(i5y_ADF_diff$ic_bic_diff)

# select the intersection of bic and aic

i5y_adq_set_diff <- as.matrix(arrange(as.data.frame(
  i5y_ADF_diff$ic_bic_diff[c(2, 5, 7, 10),]),
  const, trend, p))
i5y_adq_idx_diff <- match(data.frame(
  t(i5y_adq_set_diff[, 1:3])),
  data.frame(
    t(i5y_ADF_diff$ic_diff[, 1:3])))

for (i in 1:length(i5y_adq_idx_diff)){
  checkresiduals(
    i5y_ADF_diff$ADF_est_diff[[i5y_adq_idx_diff[i]]])
}
# All residuals look ok.

adf.test(diff(i5y))


# i90d
egr_ADF_lev <- ADF_estimate_lev(i90d, p_max = 15)
print(egr_ADF_lev$ic_aic)
print(egr_ADF_lev$ic_bic)

# intersections and 2 prefs by BIC.

egr_adq_set <- as.matrix(arrange(as.data.frame(
  rbind(egr_ADF_lev$ic_aic[c(1, 2, 5),],
        egr_ADF_lev$ic_bic[c(1, 2),])),
  const, trend, p))
egr_adq_idx <- match(data.frame(t(egr_adq_set[, 1:3])),
                     data.frame(t(egr_ADF_lev$ic[, 1:3])))

for (i in 1:length(egr_adq_idx)){
  checkresiduals(egr_ADF_lev$ADF_est[[egr_adq_idx[i]]])
}

adf.test(i90d)

# For all specs in the adequate set we reject H0 at 5% sig level and conclude
# that i90d is empirically distinguishable from a unit root process
# Since we infer that i90d is I(0), there is no need to run the test on diffs.

# i180d
i180d_ADF_lev <- ADF_estimate_lev(i180d, p_max = 15)
print(i180d_ADF_lev$ic_aic)
print(i180d_ADF_lev$ic_bic)

# intersection + top aic/bic

i180d_adq_set <- as.matrix(arrange(as.data.frame(
  rbind(i180d_ADF_lev$ic_aic[c(1:3, 5, 7, 9),],
        i180d_ADF_lev$ic_bic[c(1, 5),])),
  const, trend, p))
i180d_adq_idx <- match(data.frame(t(i180d_adq_set[, 1:3])),
                     data.frame(t(i180d_ADF_lev$ic[, 1:3])))

for (i in 1:length(i180d_adq_idx)){
  checkresiduals(i180d_ADF_lev$ADF_est[[i180d_adq_idx[i]]])
}
# all residuals look ok

adf.test(i180d)

# Again, we reject H0 at 5% sig level for some specs in the adequate set but
# not others. The ADF results are not sufficiently accurate to ascertain the
# proximity of the DGP to a unit root process.

# now repeat for the differenced i180d
i180d_ADF_diff <- ADF_estimate_diff(i180d, p_max = 15)
print(i180d_ADF_diff$ic_aic_diff)
print(i180d_ADF_diff$ic_bic_diff)

# Both AIC and BIC agree on smaller lag lengths

adf.test(diff(i180d), nlag = 10)

# For lag lengths up to 9, H0 is universally rejected at a very small sig level.
# The process generating diff(i180d) is clearly distinguishable from a unit
# root. However, we cannot say how close i180d is to an I(1) using the ADF
# testing approach.

# 2.
#


eg_reg <- lm( i5y ~ i3y + i90d + i180d, mydata)
eg_res <- eg_reg$residuals

# same approach as in Q1 but with eg_res instead of data
egr_ADF_lev <- ADF_estimate_lev(eg_res, p_max = 15)
print(egr_ADF_lev$ic_aic)
print(egr_ADF_lev$ic_bic)

# quite similar. take bic because aic has weird lag p = 15.

egr_adq_set <- as.matrix(arrange(as.data.frame(
                                    egr_ADF_lev$ic_bic),
                                    const, trend, p))
egr_adq_idx <- match(data.frame(t(egr_adq_set[, 1:3])),
                      data.frame(t(egr_ADF_lev$ic[, 1:3])))

for (i in 1:length(egr_adq_idx)){
  checkresiduals(egr_ADF_lev$ADF_est[[egr_adq_idx[i]]])
}
# All residuals look OK.

# Now, we can use the function coint.test form the aTSA package to implement
# the E-G test for the specs in the adequate set.

# example

coint.test(i5y, cbind(i3y, i90d, i180d), nlag = 4, output = T)

# loop

eg_test <- matrix(nrow = 3, ncol = 4)
colnames(eg_test) <- rep("", 4)
rownames(eg_test) <- c("No const, no trend",
                       "Const, no trend",
                       "Const with trend")
for (l in 1:4){
  eg_l <- coint.test(i5y, cbind(i3y, i90d, i180d), nlag = l, output = F)
  eg_test[, l] <- eg_l[, 3]
  colnames(eg_test)[l] <- paste("Lag", l)
}
print(eg_test)

# The best inference we can make is that if the residual in the regression i5y
# on i3y, i90d and i180d is mean-independent, then it also does not have a
# unit root. However, if mean-independence is not enforced, then we cannot
# reject a unit root in the residuals process.

# 3.
#

# go to slides...

# 4.
#
# Dependent variable: i3y
eg_reg <- lm( i3y ~ i5y + i90d + i180d, mydata)
eg_res <- eg_reg$residuals

# same approach as in Q1 but with eg_res instead of data
egr_ADF_lev <- ADF_estimate_lev(eg_res, p_max = 15)
print(egr_ADF_lev$ic_aic)
print(egr_ADF_lev$ic_bic)

# same here

egr_adq_set <- as.matrix(arrange(as.data.frame(
                                    egr_ADF_lev$ic_bic),
                                    const, trend, p))
egr_adq_idx <- match(data.frame(t(egr_adq_set[, 1:3])),
                     data.frame(t(egr_ADF_lev$ic[, 1:3])))

for (i in 1:length(egr_adq_idx)){
  checkresiduals(egr_ADF_lev$ADF_est[[egr_adq_idx[i]]])
}
# All residuals look OK.

eg_test <- matrix(nrow = 3, ncol = 4)
colnames(eg_test) <- rep("", 4)
rownames(eg_test) <- c("No const, no trend",
                       "Const, no trend",
                       "Const with trend")
for (l in 1:4){
  eg_l <- coint.test(i5y, cbind(i3y, i90d, i180d), nlag = l, output = F)
  eg_test[, l] <- eg_l[, 3]
  colnames(eg_test)[l] <- paste("Lag", l)
}
print(eg_test)

# Dependent variable: i90d
eg_reg <- lm( i90d ~ i3y + i5y + i180d, mydata)
eg_res <- eg_reg$residuals

# same approach as in Q1 but with eg_res instead of data
egr_ADF_lev <- ADF_estimate_lev(eg_res, p_max = 15)
print(egr_ADF_lev$ic_aic)
print(egr_ADF_lev$ic_bic)

egr_adq_set <- as.matrix(arrange(as.data.frame(
  egr_ADF_lev$ic_bic),
  const, trend, p))
egr_adq_idx <- match(data.frame(t(egr_adq_set[, 1:3])),
                     data.frame(t(egr_ADF_lev$ic[, 1:3])))

for (i in 1:length(egr_adq_idx)){
  checkresiduals(egr_ADF_lev$ADF_est[[egr_adq_idx[i]]])
}
# All residuals look OK.

eg_test <- matrix(nrow = 3, ncol = 4)
colnames(eg_test) <- rep("", 4)
rownames(eg_test) <- c("No const, no trend",
                       "Const, no trend",
                       "Const with trend")
for (l in 1:4){
  eg_l <- coint.test(i5y, cbind(i3y, i90d, i180d), nlag = l, output = F)
  eg_test[, l] <- eg_l[, 3]
  colnames(eg_test)[l] <- paste("Lag", l)
}
print(eg_test)

# We get the same results as with i5y being the dependent variable.

# Dependent variable: i180d
eg_reg <- lm( i180d ~ i3y + i5y + i90d, mydata)
eg_res <- eg_reg$residuals

egr_ADF_lev <- ADF_estimate_lev(eg_res, p_max = 15)
print(egr_ADF_lev$ic_aic)
print(egr_ADF_lev$ic_bic)


egr_adq_set <- as.matrix(arrange(as.data.frame(
  egr_ADF_lev$ic_bic),
  const, trend, p))
egr_adq_idx <- match(data.frame(t(egr_adq_set[, 1:3])),
                     data.frame(t(egr_ADF_lev$ic[, 1:3])))

for (i in 1:length(egr_adq_idx)){
  checkresiduals(egr_ADF_lev$ADF_est[[egr_adq_idx[i]]])
}
# All residuals look OK.

eg_test <- matrix(nrow = 3, ncol = 6)
colnames(eg_test) <- rep("", 6)
rownames(eg_test) <- c("No const, no trend",
                       "Const, no trend",
                       "Const with trend")
for (l in 1:6){
  eg_l <- coint.test(i5y, cbind(i3y, i90d, i180d), nlag = l, output = F)
  eg_test[, l] <- eg_l[, 3]
  colnames(eg_test)[l] <- paste("Lag", l)
}
print(eg_test)

# Changing the dependent variable does not materially affect our inference about
# cointegrating relations involving these four processes.

# 5.
#
# Capital market spread: i5y - i3y
cm_ADF_lev <- ADF_estimate_lev(i5y - i3y, p_max = 15)
print(cm_ADF_lev$ic_aic)
print(cm_ADF_lev$ic_bic)

# intersection 

cm_adq_set <- as.matrix(arrange(as.data.frame(
                      cm_ADF_lev$ic_bic[c(2:3, 6:8, 10),]),
                      const, trend, p))
cm_adq_idx <- match(data.frame(t(cm_adq_set[, 1:3])),
                     data.frame(t(cm_ADF_lev$ic[, 1:3])))

for (i in 1:length(cm_adq_idx)){
  checkresiduals(cm_ADF_lev$ADF_est[[cm_adq_idx[i]]])
}

adf.test(i5y - i3y, nlag = 3)

# The test rejects a unit root for models with p=1, but fails to reject for
# models with p=2.

# Money market spread: i180d - i90d
mm_ADF_lev <- ADF_estimate_lev(i180d - i90d, p_max = 15)
print(mm_ADF_lev$ic_aic)
print(mm_ADF_lev$ic_bic)

# except 0,0

mm_adq_set <- as.matrix(arrange(as.data.frame(
                        mm_ADF_lev$ic_bic[c(1:8),]),
                        const, trend, p))
mm_adq_idx <- match(data.frame(t(mm_adq_set[, 1:3])),
                    data.frame(t(mm_ADF_lev$ic[, 1:3])))

for (i in 1:length(mm_adq_idx)){
  checkresiduals(mm_ADF_lev$ADF_est[[mm_adq_idx[i]]])
}
# All residuals are OK

adf.test(i180d - i90d)

# The test rejects a unit root for all models in the adequate set. Note that
# some specifications lead to a failure to reject, but they are not in our
# adequate set, so we can ignore them! We can confidently conclude that the
# money market spread is I(0).