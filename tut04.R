library(ARDL)


# 2
#
arma2ma <- function(a, b, h){
  # we always start here
  theta0 <- b[1] / a[1]
  
  if (h == 0)
  {
    # if the horizon is zero, then just return theta0
    # note: good practice would be to check that h is an integer
    return(theta = theta0)
  }
  
  # get the orders of a(L) and b(L); in fact, the ARMA orders are p - 1 and
  # q - 1 because we also have a0 and b0 to take into account
  p <- length(a)
  q <- length(b)
  
  # augment the AR and MA coefficients vectors to match the number of thetas we
  # are going to compute -- this makes things easier later
  if (h > p)
  {
    a = c(a, numeric(1 + h - p))
  }
  
  if (h > q)
  {
    b = c(b, numeric(1 + h - q))
  }

  # allocate space for 1 + h thetas and set theta0 = b0 / a0
  theta <- c(theta0, numeric(h))
  for (j in 1:h)
  {
    theta[1 + j] <- (b[1 + j] - sum(a[2:(1 + j)] * theta[j:1])) / a[1]
  }
  
  return(theta)
}

# 3
#
ardl_irfs <- function(ardl_est, h = 40, cumirf = T){
  # extract the lag orders and coefficient estimates from the estimated ARDL
  order <- ardl_est$order
  coefficients <- ardl_est$coefficients
  
  # extract the autoregressive coefficients and construct a(L)
  j <- 1 + order[1]
  a <- c(1, -coefficients[2:j])
  
  # get the number of exogenous variables in the ARDL: we want to get IRFs
  # to each one of these separately
  k <- length(order) - 1
  
  # allocate space for all the IRFs
  irfs <- matrix(nrow = 1 + h, ncol = k)
  colnames(irfs) <- rep("", k)
  
  # allocate space for LRMs
  lrm <- numeric(k)
  names(lrm) <- rep("", k)
  
  # now, cycle through each exogenous variable and compute IRFs/LRMs
  for (i in 1:k)
  {
    # advance the index to where the estimated coefficients are stored in the
    # variable "coefficients", then extract them and construct b(L) for the ith
    # exogenous variable in the ARDL
    j0 <- 1 + j
    j <- j0 + order[1 + i]
    b <- coefficients[j0:j]
    colnames(irfs)[i] <- names(coefficients[j0])
    names(lrm)[i] <- names(coefficients[j0])

    if (cumirf)
    {
      # compute the first "h" terms of theta(L) = b(L)/a(L) if cumulative IRFs
      # are requested, do a cumulative sum of theta coefficients
      irfs[, i] <- cumsum(arma2ma(a, b, h))
    }
    else
    {
      # compute the first "h" terms of theta(L) = b(L)/a(L) and save them
      irfs[, i] <- arma2ma(a, b, h)
    }
    lrm[i] <- sum(b) / sum(a)
  }
  
  return(list(irfs = irfs, lrm = lrm))
}

# 4
#
# (a)
#
# load the data in wealth.csv
mydata <- read.delim("wealth.csv", header = TRUE,  sep = ",")

# estimate the ardl(1, 1, 2) and compute IRFs/LRMs using our function
ardl_est <- ardl(CT ~ AT + YT, mydata, order = c(1, 2, 2))

# see results
summary(ardl_est)

# compute and plot the cumulative IRFs
irfs_lrm <- ardl_irfs(ardl_est)

# see what function is doing:
irfs_lrm


for (i in 1:ncol(irfs_lrm$irfs)){
  plot(0:40, irfs_lrm$irfs[, i], type = "l",
       ylab = "Impulse Response", xlab = "Horizon",
       main = paste("Cummulative IRFs to", colnames(irfs_lrm$irfs)[i]))
}

# report the estimated LRMs
print(irfs_lrm$lrm)

#
# (b)
#
# convert the estimated ARDL(1, 1, 2) to ECM form; case = 2 here corresponds to
# no linear trend (including "t" as a regressor) and an intercept within the
# equilbrium relation (gamma = 0, mu unrestricted in our notation)

ecm_sr <- recm(ardl_est, case = 2)
summary(ecm_sr)
# note that the "recm" function does not return LRMs, but only the aggregate
# error correction term (which it calls "ect"); to get the details of this
# ect term, we need to separately use the "multipliers" function
ecm_lrm <- multipliers(ardl_est)
ecm_lrm


# (c)
#
# include the functions provided in the file "ardl_irfs_ci.R"
source("ardl_irfs_ci.R")

# estimate the confidence intervals for cumulative IRFs to both AT and YT;
# we choose 68% as the the confidence level to get approx +/-1 std dev coverage
# (in fact, Pr(-1<z<=1) = 0.68 exactly under the normal distribution)
# however, it's always a good idea to try different levels to see how sensitive

irfs_ci <- ardl_irfs_ci(ardl_est, conf = 0.68)
# Let's see results

irfs_ci

# plot all the cumulative IRFs
for (i in 1:ncol(irfs_ci$md)){
  plot(0:40, irfs_ci$md[, i], type = "l",
       ylab = "Impulse Response", xlab = "Horizon",
       main = paste("Cumulative IRFs to", colnames(irfs_ci$md)[i]),
       ylim = c(min(irfs_ci$lb[, i]), max(irfs_ci$ub[, i])))
  lines(0:40, irfs_ci$ub[, i], type = "l", col = "blue")
  lines(0:40, irfs_ci$lb[, i], type = "l", col = "blue")
}

# (d)
#
# check that the estimated IRFs with the supplied command match the estimated
# IRFs we obtain using the ardl_irfs function created in Question 3

irfs_lrm$irfs
irfs_ci$md

print(norm(irfs_lrm$irfs - irfs_ci$md))

# (e)
#
# since we want a 68% CI, we need to obtain the
# 84th percentile (i.e., 1-(1-.68)/2 = .84)
z <- qnorm(.84)

# example: 
# see table:

ecm_lrm

# get values from table:

print('Confidence interval for AT')
0.2687276 - qnorm(.84) * 0.1213092
0.2687276 + qnorm(.84) * 0.1213092

# automate
# now, construct the CI in the usual way as estimate +/ z * std error
lrm_hat <- t(as.matrix(ecm_lrm$Estimate[2:3]))
lrm_se <- t(as.matrix(ecm_lrm$`Std. Error`[2:3]))
ones <- matrix(1, nrow = 3, ncol = 1)
zz <- as.matrix(c(-z, 0, z))
lrm_ci <-  ones %*% lrm_hat + zz %*% lrm_se
rownames(lrm_ci) <- c("lb", "md", "ub")
colnames(lrm_ci) <- c("AT", "CT")
print(lrm_ci)

# compare this to the CIs for cumulative IRFs at horizon h=40
irfs_41_ci <- rbind(irfs_ci$lb[41,], irfs_ci$md[41,], irfs_ci$ub[41,])
rownames(irfs_41_ci) <- c("lb", "md", "ub")
print(irfs_41_ci)

# (f)
#
# pre-allocate space for an information criteria comparison matrix, similar to
# what we did in tutorial 3;
# unlike tutorial 3, this time we will also store the estimated models so as to
# avoid having to re-estimate them again after we've identified the adequate set
ardl_est <- list()
ic <- matrix(nrow = 125, ncol = 5)
colnames(ic) <- c("p", "l", "s", "aic", "bic")

# cycle through all the possible variants
i <- 0;
for (p in 1:5){
  for (l in 0:4)
  {
    for (s in 0:4)
    {
      i <- i + 1
      ardl_est[[i]] <- ardl(CT ~ AT + YT, mydata, order = c(p, l, s))
      ic[i,] <- c(p, l, s, AIC(ardl_est[[i]]), BIC(ardl_est[[i]]))
    }
  }
}

# look at the preference lists in terms of both AIC and BIC
ic_aic <- ic[order(ic[,4], decreasing = FALSE),][1:10,]
ic_bic <- ic[order(ic[,5], decreasing = FALSE),][1:10,]

ic_aic
ic_bic

# the first six models preferred by the BIC also have sufficient support in
# terms of the AIC, so we take this to be the adequate set
adq_set <- ic_bic[1:6,]

adq_set

# we also need to match the orders in the adequate set to the set of all models
# that we have estimated
adq_idx <- match(data.frame(t(adq_set[, 1:3])), 
                 data.frame(t(ic[, 1:3])))
adq_idx

# example: how to access stored ECM estimations: 
adq_idx[[1]]
ardl_est[[adq_idx[[1]]]]


# do a quick scan of the residuals for obvious anomalies;
for (i in 1:length(adq_idx)){
  order <- adq_set[i,1:3]
  acf(ardl_est[[adq_idx[i]]]$residuals,
      xlim = c(1, 20), xaxp = c(1, 20, 1),
      ylim = c(-0.15, 0.15), yaxp = c(-0.15, 0.15, 2),
      main = paste("Residuals ACF for ARDL(", order[1], ", ",
                                              order[2], ", ",
                                              order[3], ")", sep = ""))
}
# nothing drastic stands out so we proceed, with the set of 6 models, but noting
# that the ARDL(1,0,2) has slightly higher autocorrelations at lags 2 and 3
# compared to the other specifications

# (g)
#
#
j <- 1 # select responses to AT
shock_name <- "AT"
y_min <- Inf
y_max <- -Inf
irfs_ci <- list()
for (i in 1:length(adq_idx)){
  irfs_ci_i <- ardl_irfs_ci(ardl_est[[adq_idx[i]]], conf = 0.68)
  irfs_ci[[i]] <- cbind(irfs_ci_i$lb[, j], irfs_ci_i$md[, j], irfs_ci_i$ub[, j])
  y_min <- min(y_min, irfs_ci_i$lb[, j])
  y_max <- max(y_max, irfs_ci_i$ub[, j])
}

for (i in 1:length(adq_idx)){
  order <- adq_set[i,1:3]
  plot(0:40, irfs_ci[[i]][, 2], type = "l", ylim = c(y_min, y_max),
       ylab = "Impulse Response", xlab = "Horizon",
       main = paste("ARDL(", order[1], ", ", order[2], ", ", order[3],
                    "): Cumulative IRFs to ", shock_name, sep = ""))
  lines(0:40, irfs_ci[[i]][, 1], type = "l", col = "blue")
  lines(0:40, irfs_ci[[i]][, 3], type = "l", col = "blue")
  lines(0:40, numeric(41), type = "l", col = "red")
}
