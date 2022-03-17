#### TUTORIAL 3 SOLUTIONS ####



### Preparing 
# installing a package: install.packages('packagename')
# install package named `forecast`.


# load the forecast package using `library` function


### 1
### (a)
# load the data in Merck.csv



# Selecting dates from 1 January 2011 to 31 January 2012



### (b)
# difference and lag operators



# name the variables using `colnames()` function
colnames(y) <- "Stock Prices of Merck (MRK)"
colnames(Dy) <- "Changes in Stock Prices"
colnames(r) <- "Log Returns" 

### (c)
# we use time series plot to detect obvious trends, volatility clustering, etc
# this is typically the first thing to do when processing time series data
# Select 2011


# plot the three time series we generated.



### (d)
# plot ACFs and PACFs


# Make Comments



### (e)
# Let's fit ARMA(1,1) as an example first. 
# Example: fit = Arima(data, order = c(1, 0, 1))
# Summarize the results: summary(fit)

fit = Arima(Dy, order = c(1, 0 ,1))
summary(fit)


# How to retrive AIC and BIC?
# AIC: fit[['aic']] or fit$aic
# BIC: fit[['bic]] or fit$bic

print(fit$aic)
print(fit$bic)

# Create a matrix with size 25x4. We will 
# save results into this empty matrix:


ic <- matrix( nrow = 25, ncol = 4 )

# name the columns of `ic` for simplicity
colnames(ic) = c('p','q','aic','bic')

# Now, we create a nested for loop iterating for both
# q and p from 0 to 4. We will save results to `ic` 



# EXTRA: try to use auto.arima(data, seasonal = FALSE)

fit.auto = auto.arima(Dy, seasonal = FALSE)
summary(fit.auto)

### (f)
# find small aic and bics

ic[order(ic[,'aic']),][1:5,] # top 5
ic[order(ic[,'bic']),][1:5,] # top 5


# next, we select the models that make top 10 in both lists;
# (not 100% objective)

adq_set = list(c(1, 0, 1), c(1, 0, 2), c(2, 0, 1), c(1, 0, 0))

### (g)
#
# once we have reduced the set, perform residual analysis;
# the "forecast" package simplifies this with the command "checkresiduals"




### (h)
# Forecast changes in MRK stock prices in January, 2012.




### (i)
# repeat the forecasting steps with the "y" variable using an ARMA(2,1)




### (j)
# Repeat parts (d)-(h) with the return variable (r)

