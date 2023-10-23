
***********************************************************
* Tutorial session 12: Limited dependent variables
***********************************************************

cscript


***********************************************************
* Problem I: racial discrimination in the mortgage market
***********************************************************

use "loanapp.dta", clear
de

* I-(1)

*SLIDES

* I-(2)
reg approve white, r


* I-(3)
 
reg approve white hrat obrat loanprc unem male married dep sch cosign chist pubrec mortlat1 mortlat2 vr, r


* I-(4)

probit approve white, r

predict p_white_probit if white==1
predict p_black_probit if white==0
sum p_*
* Difference in predicted probability
di .9083878-.7077922

* note that it is similar to LPM

* Alternative 1

di normal(_b[_cons] + _b[white]) - normal(_b[_cons]) 

* alternative 2
* Note here we use "i.white" instead of "white" to let STATA knows we are working on a categorical variable. Using STATA’s factor variable notation in the estimation command for margins is necessary to compute correct results!

probit approve i.white, r
margins, dydx(white)

* I-(5) 
probit approve white hrat obrat loanprc unem male married dep sch cosign chist pubrec mortlat1 mortlat2 vr, r


preserve
replace white = 1
* Predicted probability of approval for white
predict p_white
replace white = 0
* Predicted probability of approval for non-white
predict p_black
sum p_white p_black
di .8965419 - .7923174
restore

* The computation here seems different from what we did in (4)? Not really...
* (4) is a special case where there are no other controls, so the marginal effect is simply: "normal(b0 + b1) - normal(b0)"
* But here, we have other controls so the marginal effect is: P(Y=1|white =1, other X) - P(Y=1|white =0, other X)

* Alternative (easy)

* The same result is obtained with margins:

probit approve i.white hrat obrat loanprc unem male married dep sch cosign chist pubrec mortlat1 mortlat2 vr, r
margins, dydx(white)

* I-(6) 
logit approve white, r



predict p_white_logit if white == 1
predict p_black_logit if white == 0
sum p_white_logit p_black_logit
di .9083878 -.7077922


* And directly:
logit approve i.white, r
margins, dydx(white)




 
**********************************************
* Problem II: job training and unemployment
***********************************************

cscript
use "jtrain2.dta", clear
de 

* II-(1)
tab train, m

tab mostrn
* 24 months is the longest, ie all the way from Jan 1976 to Dec 1977.

* II-(2)
reg train unem74 unem75 age educ black hisp married, r

test unem74 unem75 age educ black hisp married

* II-(3)
probit train unem74 unem75 age educ black hisp married

test unem74 unem75 age educ black hisp married

* II-(4)

* discuss

* II-(5)
reg unem78 train, r


* II-(6)
probit unem78 train, r

probit unem78 i.train, r
margins, dydx(train)


* Or alternatively you can find the marginal effects using predictions:
predict p_train if train==1
predict p_notrain if train==0
su p*train
di .2432432-.3538462 


* II-(7)
* LPM
reg unem78 train, r
predict p_train_LPM
tab p_train_LPM
* Probit
probit unem78 train, r
predict p_train_prob
tab p_train_prob

pwcorr p_train_prob p_train_LPM

* II-(8)
* Probit
probit unem78 i.train unem74 unem75 age educ black hisp married, r
margins, dydx(train)


reg unem78 train unem74 unem75 age educ black hisp married, r



