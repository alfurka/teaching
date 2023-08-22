use "fertil2.dta", clear


* I-(1)
de

* I-(2)
reg children educ age agesq, r


* I-(3)

pwcorr educ frsthalf


reg educ frsthalf age agesq,r

*******> SLIDES

* I-(4)
ivregress 2sls children age agesq (educ = frsthalf), first r


******************************************************
* Problem II: The effects of smoking on birth weight
******************************************************

cscript

use "bwght.dta", clear

* II-(1)
de

* II-(2)
reg lbwght packs faminc, r


* II-(3)

*******> SLIDES

* II-(4)

twoway (scatter packs cigprice) (lfit packs cigprice)

* II-(5)

reg lbwght cigprice faminc, r


reg packs faminc cigprice, r 
test cigprice


* II-(6)

* (5) and (3) imply that we should not use IV, but we implement IV to see what happens
ivregress 2sls lbwght faminc (packs = cigprice), r



***************************************************************
* Problem III (Conditional IV): The effect of education on wages
***************************************************************

cscript

use "card.dta", clear

* III-(1)
de
sum

* III-(2)
reg lwage educ exper expersq smsa south, r



* III-(3)


reg educ exper expersq smsa south nearc4,r
test nearc4

*******> SLIDES

* III-(4)

ivregress 2sls lwage exper expersq smsa south  (educ = nearc4)

ivregress 2sls lwage exper expersq smsa south  (educ = nearc4), r

* III-(5)
reg IQ nearc4, r

*******> SLIDES

* III-(6)


reg IQ nearc4 smsa66 reg662 reg663 reg664 reg665 reg666 reg667 reg668 reg669 , r



* III-(7)


*******> SLIDES

* III-(8)

reg educ exper expersq nearc4 smsa66 reg662 reg663 reg664 reg665 reg666 reg667 reg668 reg669,r
test nearc4

ivregress 2sls lwage exper expersq smsa66 reg662 reg663 reg664 reg665 reg666 reg667 reg668 reg669 (educ = nearc4), first r
