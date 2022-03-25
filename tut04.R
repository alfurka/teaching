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


mydata = read.csv('wealth.csv')

fit = ardl(CT ~ AT + YT, data =mydata, order = c(1, 2 ,2))

irfs_lrm = ardl_irfs(fit)

for (i in 1:ncol(irfs_lrm$irfs)) {
  h = 1:length(irfs_lrm$irfs[,i])
  
  plot(h, irfs_lrm$irfs[,i], type = 'l', 
       xlab = 'Horizon', ylab = 'IRFs', 
       main = paste0('Cumulative IRFs - ',colnames(irfs_lrm$irfs)[i]),
       col = 'red', lwd = 2)
}

print(irfs_lrm$lrm)

# b

ecm_sr = recm(fit , case = 2)
ecm_lrm = multipliers(fit)

# Error in loadNamespace(name) : there is no package called 'msm'
# remove R and install version 3.6.3
# https://cran.r-project.org/bin/windows/base/old/3.6.3/

summary(ecm_sr)
print(ecm_lrm)

# c 

source('ardl_irfs_ci.R')

irfs_ci = ardl_irfs_ci(fit, conf = 0.68)



for (i in 1:ncol(irfs_lrm$irfs)) {
  h = 1:length(irfs_lrm$irfs[,i])
  
  plot(h, irfs_ci$md[,i], type = 'l', 
       xlab = 'Horizon', ylab = 'IRFs', 
       main = paste0('Cumulative IRFs - ',colnames(irfs_ci$md)[i]),
       col = 'black', lwd = 2,
       ylim = c(min(irfs_ci$lb[,i]) , max(irfs_ci$ub[,i])))
  lines(h, irfs_ci$lb[,i], type = 'l', col = 'red')
  lines(h, irfs_ci$ub[,i], type = 'l', col = 'red')
}



# d 

norm(irfs_lrm$irfs - irfs_ci$md)



# e 

print(ecm_lrm)
print(ecm_lrm$estimate)
print(ecm_lrm$std.error)

z = qnorm(0.84)

lrm_ci = as.matrix(cbind(ecm_lrm$estimate - z *ecm_lrm$std.error,
                         ecm_lrm$estimate,
                         ecm_lrm$estimate + z *ecm_lrm$std.error))
rownames(lrm_ci) = ecm_lrm$term
colnames(lrm_ci) = c('LB','M','UB')
lrm_ci

irfs_41_ci = as.matrix(cbind(irfs_ci$lb[41,],
                             irfs_ci$md[41,],
                             irfs_ci$ub[41,]))
irfs_41_ci

colnames(irfs_41_ci) = c('LB','M','UB')
irfs_41_ci


# f


ardl_est = list()

ic = matrix(ncol = 6, nrow = 5 * 5 * 5)
colnames(ic) = c('i', 'p','l','s','aic','bic')
i = 0
for (p in 1:5) {
  for (l in 0:4) {
    for (s in 0:4) {
      i = i + 1
      ardl_est[[i]] = ardl(CT ~ AT + YT, mydata, order = c(p, l, s))
      
      ic[i, ] = c(i, p ,l ,s , AIC(ardl_est[[i]]), BIC(ardl_est[[i]]))
    }  
  }
  
}

ic[order(ic[,'aic']),][1:10,]
ic[order(ic[,'bic']),][1:10,]

adq_set = ic[order(ic[,'bic']),][1:6,]


for (i in 1:nrow(adq_set)) {
  acf_title = paste0('Residual ACF of ARDL(', adq_set[i,'p'], 
                     ', ', adq_set[i,'l'], 
                     ', ', adq_set[i,'s'], 
                     ')')
  acf(ardl_est[[adq_set[i,'i']]]$residuals, main = acf_title)
}


# g

j = 1 # select response to AT

irfs_ci = list()

irfs_ci_i = ardl_irfs_ci(ardl_est[[adq_set[1,'i']]], conf = 0.68)
irfs_ci_i$md

y_min = Inf
y_max  = -Inf

for (i in 1:nrow(adq_set)) {
  irfs_ci_i = ardl_irfs_ci(ardl_est[[adq_set[i,'i']]], conf = 0.68)
  irfs_ci[[i]] = cbind(irfs_ci_i$lb[,j], 
                       irfs_ci_i$md[,j], 
                       irfs_ci_i$ub[,j])
  
  y_min = min(y_min , irfs_ci_i$lb[,j])
  y_max = max(y_max , irfs_ci_i$ub[,j])
}

for (i in 1:nrow(adq_set)) {
  plot(0:40, irfs_ci[[i]][, 2], 
       type = 'l', 
       ylim = c(y_min, y_max),
       xlab = 'IRFs', ylab = 'Horizon',
       main = paste0('Cumulative IRFs to AT - ARDL(', adq_set[i,'p'], 
                     ', ', adq_set[i,'l'], 
                     ', ', adq_set[i,'s'], 
                     ')'))
  lines(0:40, irfs_ci[[i]][, 1], col = 'red')
  lines(0:40, irfs_ci[[i]][, 3], col = 'red')
}
