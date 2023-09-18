**************************************************************
* Tutorial session 7: Stata applications of panel data method
**************************************************************

*****************************************************
* Problem I: Construction of panel data
*****************************************************

cscript

* I-(1)
insheet using "psidw.csv"
de

* I-(2)

browse

reshape long exper married female union educ blcak lwage, i(id) j(year)

* i(cross-section-dimension) j(time-dimension): note how reshape takes the suffixes from all the variable listed (exper, married...) and inputs that into the "year" variable.
sort id year
browse

* I-(3)

*define panel data
xtset id year
replace year=1976 if year==1
replace year=1977 if year==2
replace year=1978 if year==3
replace year=1979 if year==4
xtset id year

* I-(4)

*construct experience square 
gen exper2=exper^2

reg lwage educ union married exper exper2 blcak female, vce(cluster id)



reg lwage educ union married exper exper2 blcak female if year==1976, vce(robust)
reg lwage educ union married exper exper2 blcak female if year==1977, vce(robust)
reg lwage educ union married exper exper2 blcak female if year==1978, vce(robust)
reg lwage educ union married exper exper2 blcak female if year==1979, vce(robust)

* why do we look each year separately?

* I-(5)

xtreg lwage educ union married exper exper2 blcak female, re vce(cluster id)

xtreg lwage educ union married exper exper2 blcak female, fe vce(cluster id)


** with time fixed effects

xtreg lwage educ union married exper exper2 blcak female i.year, re vce(cluster id)

xtreg lwage educ union married exper exper2 blcak female i.year, fe vce(cluster id)


* I-(6)

xtreg lwage educ union married exper exper2 blcak female, fe 
*store the fixed effects estimator
estimates store FE


xtreg lwage educ union married exper exper2 blcak female, re 
*store the random effects estimator
estimates store RE


* Now the Hausman test
hausman FE RE, sigmamore
   

* I-(7)
* discuss



**********************************************
* Problem II: regression discontinuity design
**********************************************

* II-(1)
cscript
use "schautonomy1.dta", clear
de
su
*passrate0  is the school performance before
*passrate2 is the school performance two years later
*dpass=passrate2-passrate0, i.e. the difference in school performance after / before the introduction of GMs.

* II-(2)
reg dpass win, r


* II-(3)
*vote:vote share that expresses the percentage that voted in favour of GM schools
scatter dpass vote

* II-(4)
* Linear fit

scatter dpass vote if vote>.15&vote<.85 || lfit dpass vote if vote>.15&vote < .50 || lfit dpass vote if vote>= .50&vote<.85 ,xline(0)  legend(order(2 "vote < 50" 3 "vote >= 50"))

* let's decompose the code above:
*1 
scatter dpass vote if vote>.15&vote<.85 

*2
scatter dpass vote if vote>.15&vote<.85 || lfit dpass vote if vote>.15&vote < .50 
*3
scatter dpass vote if vote>.15&vote<.85 || lfit dpass vote if vote>.15&vote < .50 || lfit dpass vote if vote>= .50&vote<.85

*4
scatter dpass vote if vote>.15&vote<.85 || lfit dpass vote if vote>.15&vote < .50 || lfit dpass vote if vote>= .50&vote<.85 ,xline(0)  legend(order(2 "vote < 50" 3 "vote >= 50"))



* Quadratic fit
scatter dpass vote if vote>.15&vote<.85 || qfit dpass vote if vote>.15&vote < .50 || qfit dpass vote if vote>= .50&vote<.85 ,xline(0)  legend(order(2 "vote < 50" 3 "vote >= 50"))


* II-(5)

reg dpass win if vote>.15&vote<.85, r
reg dpass win vote if vote>.15&vote<.85, r

* An alternative is to control for the distance to the threshold (piecewise linearly or with quadratic terms).
* Distance to the threshold for those who lost 
g lose_vote=0 if win==1
replace lose_vote=-(vote-0.5) if win==0
* imagine a school got 0.40 votes. This school is not going to become a GM.
* 0.40 is below 0.50 so this school will not be converted to a GM school
* lose_vote=0.10 for this school (the distance from cutoff)
* Distance to the threshold for those who won
g win_vote=0 if win==0
replace win_vote=(vote-0.5) if win==1
* imagine a school got 0.60 votes. This school is going to become a GM.
* win_vote=0.10 for this school (the distance from cutoff)
reg dpass win lose_vote win_vote if vote>.15&vote<.85, r /* the coefficient on wins decreases slightly but stays positive and significant */
* Also include the 2nd order polynomials
g lose_vote_2=lose_vote^2
g win_vote_2=win_vote^2
reg dpass win lose_vote lose_vote_2 win_vote win_vote_2 if vote>.15&vote<.85, r 



* II-(6)
* SLIDES

* II-(7)
reg passrate2 win lose_vote win_vote if vote>.15&vote<.85, r

bysort win: sum passrate0 passrate2


* II-(8)
reg passrate0 win lose_vote win_vote if vote>.15&vote<.85, r

* What happens if it was significant?


****************
*****************
