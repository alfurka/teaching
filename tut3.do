
* I-(1)
de

* I-(2)

codebook countryn
list countryn


* I-(3)
twoway (scatter lgdp prot) (lfit lgdp prot)
* help graph
* help twoway
* help scatter 
* help lfit
twoway (scatter lgdp prot) (lfitci lgdp prot)
twoway (scatter lgdp prot) (qfit lgdp prot)


twoway (scatter prot logmort) (lfit prot logmort)
twoway (scatter prot logmort) (lfitci prot logmort)


* I-(4) 
reg lgdp logmort, r

reg prot logmort, r

di (-.5697983)/(-.6213181)



* I-(5) 

ivregress 2sls lgdp (prot=logmort)


* I-(6) 
reg prot logmort, r
predict prothat, xb
* help predict
reg lgdp prothat, r



* I-(7) 

* I-(8) 



* I-(9) 
ivregress 2sls lgdp euro (prot=logmort)


****************************************************************
* Problem II: measurement error
****************************************************************

cscript

* For this question, since we generate data taking random draws from a normal distribution,
* the estimates will be slightly different each time we generate the data
* unless we fix them with the set seed command
set seed 123


* II-(1)
set obs 500
*generate x1 using a standard normal distribution 
gen x1=rnormal(0,1)
* for more details on random number generation, type in: help rnormal
*generate u using a standard normal distribution 
gen u1=rnormal(0,1)
gen y=0+1*x1+u1


* II-(2)

gen v=rnormal(0,1)

gen x=x1+v


* II-(3)

reg y x1

reg y x



* II-(4)
drop x y x1 u v
	* redo II-(1)
set obs 1000

* REDO previous steps


*II-(5)



*II-(6)



*II-(7)

