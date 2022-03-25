# install.packages('ARDL')
library(ARDL)


# Question 2
arma2ma <- function(a, b, h)
{
  # we always start here
  theta0 <- b[1] / a[1]
  
  if (h == 0)  {
    # if the horizon is zero, then just return theta0
    return(theta = theta0)
  }
  
  # get the orders of a(L) and b(L); in fact, the ARMA orders are p - 1 and
  # q - 1 because we also have a0 and b0 to take into account
  p <- length(a)
  q <- length(b)
  
  # augment the AR and MA coefficients vectors to match the number of thetas we
  if (h > p)  {
    a = c(a, numeric(1 + h - p))
  }
  
  if (h > q)  {
    b = c(b, numeric(1 + h - q))
  }

  # allocate space for 1 + h thetas and set theta0 = b0 / a0
  theta <- c(theta0, numeric(h))
  for (j in 1:h)  {
    theta[1 + j] <- (b[1 + j] - sum(a[2:(1 + j)] * theta[j:1])) / a[1]
  }
  
  return(theta)
}

# Question 3
#
ardl_irfs <- function(ardl_est, h = 40, cumirf = T) {
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
  for (i in 1:k)  {
    # advance the index to where the estimated coefficients are stored in the
    # variable "coefficients", then extract them and construct b(L) for the ith
    # exogenous variable in the ARDL
    j0 <- 1 + j
    j <- j0 + order[1 + i]
    b <- coefficients[j0:j]
    colnames(irfs)[i] <- names(coefficients[j0])
    names(lrm)[i] <- names(coefficients[j0])

    if (cumirf)    {
      # compute the first "h" terms of theta(L) = b(L)/a(L) if cumulative IRFs
      # are requested, do a cumulative sum of theta coefficients
      irfs[, i] <- cumsum(arma2ma(a, b, h))
    }
    else  {
      # compute the first "h" terms of theta(L) = b(L)/a(L) and save them
      irfs[, i] <- arma2ma(a, b, h)
    }
    lrm[i] <- sum(b) / sum(a)
  }
  
  return(list(irfs = irfs, lrm = lrm))
}

# Question 4
#
# (a)




# b

# Error in loadNamespace(name) : there is no package called 'msm'
# remove R and install version 3.6.3
# https://cran.r-project.org/bin/windows/base/old/3.6.3/



# c 

source('ardl_irfs_ci.R')





# d 





# e 

print(ecm_lrm)
print(ecm_lrm$estimate)
print(ecm_lrm$std.error)

z = qnorm(0.84)





# f




# g

