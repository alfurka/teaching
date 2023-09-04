***********************************************************
* Tutorial session 6: Stata Applications of DiD 
***********************************************************

******************************************************************************************
* Problem I: The effect of worker compensation on weekd out of work
*****************************************************************************************

cscript

* I- (1)
* Loading data
use "injury.dta", clear
* Describe data

de

* I-(2)

* Diff-in-diff estimate: (ldurat_after_highearn - ldurat_before_highearn) - (ldurat_after_lowearn - ldurat_before_lowearn)

mean ldurat if afchnge==1&highearn==1&ky==1

mean ldurat if afchnge==0&highearn==1&ky==1

mean ldurat if afchnge==1&highearn==0&ky==1

mean ldurat if afchnge==0&highearn==0&ky==1

* DID estimate:  (1.580352 -1.382094 )- (1.133273 -1.125615)
di (1.580352 -1.382094 )- (1.133273 -1.125615)

* as a regression: 


reg ldurat afchnge highearn afhigh if ky==1
reg ldurat afchnge highearn afhigh if ky==1, robust



* I-(3)

reg ldurat afchnge highearn afhigh male married i.indust i.injtype if ky==1, robust

* I-(4)

* Discuss

* I-(5)
reg ldurat afchnge highearn afhigh if mi==1, robust




***************************************************************
* Problem II: Minimun wages and employment
***************************************************************

cscript

* II-(1)
* Loading data
use "minwage.dta", clear
*Describing the data
de


* II-(2)
tab nj after, su(fte) means 

di (17.583627 - 17.301056) - (18.253846 -20.3) 
* 2.328725
* Surprisingly, employment rose in NJ relative to PA after the rise in minimum wage.


* II-(3)


reg fte nj after njafter

reg fte nj after njafter, robust



* II-(4)
xtset sheet after
xtreg fte nj after njafter, fe robust


* II-(5)

* Discuss


* II-(6)


* Discuss
