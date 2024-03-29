***********************************************************
* Tutorial session 11: Quantile regression
***********************************************************

cscript


*****************************************************
* Problem I: OLS and quantile regression
*****************************************************

* I-(1)
use "PS1incomes.dta", clear
describe
tabstat hincome, by(year) s(mean sd p10 p25 p50 p75 p90)

* I-(2)
<<<<<<< HEAD
*OLS estimate:
reg hincome educ exper expersq,r

=======

reg hincome educ exper expersq,r
>>>>>>> c2c227fca7b23f3c9411145f4a90fdf21647c17e

qreg hincome educ exper expersq, quantile(.1) vce(robust)
qreg hincome educ exper expersq, quantile(.25) vce(robust)
qreg hincome educ exper expersq, quantile(.50) vce(robust)
qreg hincome educ exper expersq, quantile(.75) vce(robust)
qreg hincome educ exper expersq, quantile(.9) vce(robust)

<<<<<<< HEAD

=======
>>>>>>> c2c227fca7b23f3c9411145f4a90fdf21647c17e
* I-(3)
reg hincome educ exper expersq if year==1996,r
reg hincome educ exper expersq if year==1999,r 
reg hincome educ exper expersq if year==2002,r
reg hincome educ exper expersq if year==2005,r

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

<<<<<<< HEAD

=======
>>>>>>> c2c227fca7b23f3c9411145f4a90fdf21647c17e
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
<<<<<<< HEAD
=======
* The coefficients are the same but the standard errors are different.
>>>>>>> c2c227fca7b23f3c9411145f4a90fdf21647c17e

* I-(6)

kdensity hincome if year==1996, lcolor(green) addplot (kdensity hincome if year==2005)


*****************************************************
* Problem II: RCT and quantile regression
*****************************************************

clear all
use "STAR_public_use.dta", clear

* II-(1)
reg signup ssp sfp sfsp  
reg used_ssp ssp sfp sfsp 
reg used_adv ssp sfp sfsp 
reg used_fsg ssp sfp sfsp  

* II-(2)
reg ssp female, r
<<<<<<< HEAD
reg signup female, r
reg used_ssp female, r
reg used_adv female, r
reg used_fsg female, r

* II-(3)

=======

reg signup female, r

reg used_ssp female, r

reg used_adv female, r

reg used_fsg female, r

* II-(3)
* Effect of the treatments on year 1 gpa
>>>>>>> c2c227fca7b23f3c9411145f4a90fdf21647c17e
reg GPA_year1 ssp sfp sfsp,r 
reg GPA_year1 ssp sfp sfsp if sex=="F",r
reg GPA_year1 ssp sfp sfsp if sex=="M",r

<<<<<<< HEAD

=======
* Effect of the treatments on year 2 gpa
>>>>>>> c2c227fca7b23f3c9411145f4a90fdf21647c17e
reg GPA_year2 ssp sfp sfsp,r 
reg GPA_year2 ssp sfp sfsp if sex=="F",r 
reg GPA_year2 ssp sfp sfsp if sex=="M",r 

<<<<<<< HEAD

* II-(4)

=======
* II-(4)
* Girls
>>>>>>> c2c227fca7b23f3c9411145f4a90fdf21647c17e
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.10) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.25) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.50) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.75) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="F", quantile(0.90) vce(robust)

<<<<<<< HEAD

=======
* Boys
>>>>>>> c2c227fca7b23f3c9411145f4a90fdf21647c17e
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.10) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.25) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.50) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.75) vce(robust)
qreg GPA_year2 ssp sfp sfsp if sex=="M", quantile(0.90) vce(robust)

<<<<<<< HEAD









=======
>>>>>>> c2c227fca7b23f3c9411145f4a90fdf21647c17e
