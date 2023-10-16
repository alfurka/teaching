***********************************************************
* Tutorial session 11: Quantile regression
***********************************************************

cscript


*****************************************************
* Problem I: OLS and quantile regression
*****************************************************

* I-(1)
use "D:\Protected\Teaching\ECON 3360\2023\Course material\Tutorials\tutorial_11\tute_11\PS1incomes.dta", clear
describe
tabstat hincome, by(year) s(mean sd p10 p25 p50 p75 p90)
* In 1996, 25% of the population earnt less than 4.36 pounds an hour (and 75% of the population earnt more).

* I-(2)
*OLS estimate:
reg hincome educ exper expersq,r
* On average over the period, individuals with 1 more year of education earned an extra income of 0.79 pounds per hour (Note that nothing suggests that the coefficients can be interpreted causally).
*Quantile estimates:
qreg hincome educ exper expersq, quantile(.1) vce(robust)
qreg hincome educ exper expersq, quantile(.25) vce(robust)
qreg hincome educ exper expersq, quantile(.50) vce(robust)
qreg hincome educ exper expersq, quantile(.75) vce(robust)
qreg hincome educ exper expersq, quantile(.9) vce(robust)
* Over the period, individuals with 1 more year of education earned an extra income of 0.42 pounds per hour at the first quartile.


* I-(3)
reg hincome educ exper expersq if year==1996,r
reg hincome educ exper expersq if year==1999,r 
reg hincome educ exper expersq if year==2002,r
reg hincome educ exper expersq if year==2005,r
* The coefficient of education increases over time, so the returns to education are not constant over time.
* 1 more year of education is associated with higher and higher earnings over time (this could be due to education having a larger impact or to OVB changing over time, for example if education selects more productive people over time).

* I-(4)
*Mean
reg hincome educ exper expersq if year==1996,r
reg hincome educ exper expersq if year==2005,r
*QR 0.1
qreg hincome educ exper expersq if year==1996, quantile(.1) vce(robust)
qreg hincome educ exper expersq if year==2005, quantile(.1) vce(robust)
*QR 0.25
qreg hincome educ exper expersq if year==1996, quantile(.25) vce(robust)
qreg hincome educ exper expersq if year==2005, quantile(.25) vce(robust)
*QR 0.5
qreg hincome educ exper expersq if year==1996, quantile(.50) vce(robust)
qreg hincome educ exper expersq if year==2005, quantile(.50) vce(robust)
*QR 0.75
qreg hincome educ exper expersq if year==1996, quantile(.75) vce(robust)
qreg hincome educ exper expersq if year==2005, quantile(.75) vce(robust)
*QR 0.9
qreg hincome educ exper expersq if year==1996, quantile(.90) vce(robust)
qreg hincome educ exper expersq if year==2005, quantile(.90) vce(robust)
* The mean returns to education increases from 0.56 to 0.90, it is multiplied by 1.6.
* No the increase in the coefficients is not similar for the different quantiles. 
* It is only multiplied by 1.27 increases for the first decile (from 0.196 to 0.249) but by 1.58 for the third quartile and the last decile (e.g. from 0.91 to 1.43 for quartile 0.9).

* I-(5)
*QR 0.1
bsqreg hincome educ exper expersq if year==1996, quantile(.1) reps(10)
bsqreg hincome educ exper expersq if year==2005, quantile(.1) reps(10)
*QR 0.25
bsqreg hincome educ exper expersq if year==1996, quantile(.25) reps(10)
bsqreg hincome educ exper expersq if year==2005, quantile(.25) reps(10)
*QR 0.5
bsqreg hincome educ exper expersq if year==1996, quantile(.50) reps(10)
bsqreg hincome educ exper expersq if year==2005, quantile(.50) reps(10)
*QR 0.75
bsqreg hincome educ exper expersq if year==1996, quantile(.75) reps(10)
bsqreg hincome educ exper expersq if year==2005, quantile(.75) reps(10)
*QR 0.9
bsqreg hincome educ exper expersq if year==1996, quantile(.90) reps(10)
bsqreg hincome educ exper expersq if year==2005, quantile(.90) reps(10)
* The coefficients are the same but the standard errors are different.

* I-(6)
*Individuals at the top of the income distribution have higher levels of education in general and each extra year of education is asociated with higher incomes at the top than at the bottom. In addition this difference increased over time. Income inequality is most likely increasing.
*To check you can draw the density of hourly income showing that in 2005 the distribution is flatter and with a fatter right-hand tail.
kdensity hincome if year==1996, lcolor(green) addplot (kdensity hincome if year==2005)


*****************************************************
* Problem II: RCT and quantile regression
*****************************************************

clear all
use "D:\Protected\Teaching\ECON 3360\2023\Course material\Tutorials\tutorial_11\tute_11\STAR_public_use.dta", clear

* II-(1)
reg signup ssp sfp sfsp  
reg used_ssp ssp sfp sfsp 
reg used_adv ssp sfp sfsp 
reg used_fsg ssp sfp sfsp  
* Being allocated to a treatment group significantly increases the probability to sign up for STAR: by 49pp for SSP, 84pp for SFP and by 71pp for the combined treatment.
* As expected only allocation to SSP or the combined treatment increases the probability to use SSP, respectively by 22pp and 38pp. Note that these coefficients are far from 100pp (=perfect compliance). In contrast if allocated to SFP the probability to use SSP is not different to that of the control group. This suggest that students do not self-transfer between programs.
* In the same way, only allocation to SSP or the combined treatment increases the probability to have interacted with an advisor and attended a FSG meeting (with a smaller effect on the latter).

* II-(2)
reg ssp female, r
*girls are not more likely to be allocated to SSP (this is what we expect in a RCT)
reg signup female, r
*girls are 5pp more likely to sign up for STAR
reg used_ssp female, r
*girls are 2.5pp more likely to receive SSP services
reg used_adv female, r
*girls are 2.8pp more likely to meet with or email an advisor
reg used_fsg female, r
*girls are NOT more likely to attend FSGs

* II-(3)
* Effect of the treatments on year 1 gpa
reg GPA_year1 ssp sfp sfsp,r 
reg GPA_year1 ssp sfp sfsp if sex=="F",r
reg GPA_year1 ssp sfp sfsp if sex=="M",r
* None of the treatment has any effect on students, boys or girls
* Effect of the treatments on year 2 gpa
reg GPA_year2 ssp sfp sfsp,r 
reg GPA_year2 ssp sfp sfsp if sex=="F",r 
reg GPA_year2 ssp sfp sfsp if sex=="M",r 
* None of the treatment has any effect on students, but the combined SFSP treatment increases girls' GPA (significant at 10%) while it decreases boys' GPA (significant at 10%)

* II-(4)
* Girls
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.10) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.25) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.50) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.75) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.90) vce(robust)
*SFP and SFSP increase girls' median GPA in year 2. It has no effect on really low or high GPAs.
* Boys
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.10) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.25) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.50) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.75) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.90) vce(robust)
*SFSP decreases boys' first decile, first quartile and median GPA in year 2. It has no effect on really high GPAs.










