***********************************************************
* Tutorial session 10: PSM
************************************************************

cscript


* I-(1) 
findit psmatch2
ssc install psmatch2
help psmatch2

* I-(2)
use "training.dta", clear
de
 

* I-(3)
tab treated /* there are 297 treated indivuduals and 425 controls */
reg age treated
* age is similar in both treatment group and control group
reg educ treated
reg nodegree treated
reg black treated
reg hisp treated
reg married treated
reg re74 treated
reg re75 treated



* Alternative ways of running the balancing test is by running ttests:
ttest age, by(treated)
* Or probit estimations (more details about this in a few weeks)
probit treated age 
* It is also useful to run an F-test with all controls simultaneously
reg treated age educ nodegree black hisp married re74 re75


* I-(4) 
reg re78 treated, r

reg re78 treated age educ nodegree black hisp married re74 re75, r

* With age square
gen age2=age^2
reg re78 treated age age2 educ black hisp married nodegree re75 age2, r

****************************************************************************************
* Propensity score matching approach
*****************************************************************************************

* I-(5) We use an alternative control group from the PSID i.e. non-experimental survey data
tab sample
	* 1 is experimental data and 3 is PSID data (non-experimental survey data)
gen treated2=1 if treated==1
replace treated2=0 if sample==3
tab treated2 sample, m 


reg age treated2
reg educ treated2
reg nodegree treated2
reg black treated2
reg hisp treated2
reg married treated2
reg re74 treated2
reg re75 treated2


* I-(6)
* OLS with controls
reg re78 treated2 age age2 educ nodegree black hisp married re74 re75, r

* I-(7) 

* See SLIDES

* I-(8) Propensity score
probit treated2 age age2 educ nodegree black hisp married re74 re75
predict ps
margins, dydx(*) /* This gives you the marginal effects which is the common way of presenting estimates from a probit or a logit */
* Propensity score histogram by treatment status
psgraph, treated(treated2) pscore(ps) bin(50) 


* I-(9) Caliper matching (matching to nearest neighbour within caliper)
psmatch2 treated2, pscore(ps) outcome(re78) neighbor(1) caliper(0.01) 
psmatch2 treated2, pscore(ps) outcome(re78) neighbor(1) caliper(0.005) 
psmatch2 treated2, pscore(ps) outcome(re78) neighbor(1) caliper(0.025) 

* I-(10) Caliper matching with strict common support
psmatch2 treated2, pscore(ps) outcome(re78) neighbor(1) caliper(0.025) common 

* I-(11) NN
psmatch2 treated2, pscore(ps) outcome(re78) 

* similar one with calipher

psmatch2 treated2, pscore(ps) outcome(re78) neighbor(1) caliper(0.04)

* I-(12) Radius (=all neighbours within the caliper)
psmatch2 treated2, radius pscore(ps) outcome(re78) neighbor(1) caliper(0.025)

* Kernel Method 
psmatch2 treated2, pscore(ps) outcome(re78) kernel k(epan) bw(0.025)
	* The kernel used all control and the ATT is close to the radius (-876)

* I-(13) OLS using PSM weights and controls
psmatch2 treated2, pscore(ps) outcome(re78) neighbor(1) caliper(0.025)
reg re78 treated2 age age2 educ nodegree black hisp married re74 re75 [fweight=_weight]

* compare 

reg re78 treated age age2 educ nodegree black hisp married re74 re75, r


* I-(14) Re-run the propensity score estimate dropping age2 and re74 from the regressors
drop ps
probit treated2 age educ nodegree black hisp married re75 
predict ps
psmatch2 treated2, pscore(ps) outcome(re78) neighbor(1) caliper(0.025)
reg re78 treated2 age age2 educ nodegree black hisp married re74 re75 [fweight=_weight]

* I-(15)


pstest age educ nodegree black hisp married re75, graph treated(treated2) mweight(_weight) both
