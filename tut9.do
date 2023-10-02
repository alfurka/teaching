****************************************************************
* Tutorial session 9: Regression Discontinuity Design (fuzzy)
****************************************************************


cscript



* I-(1)
use "maimonides.dta", clear
de

* I-(2)
sum classize, detail
sum avgmath,detail
sum avgverb,detail
sum perc_disadvantaged, detail

* I-(3)

foreach var in classize avgmath avgverb perc_disadvantaged {
sum `var', detail
}


* I-(4)
sort enrollment
by enrollment: egen classize_m=mean(classize)
plot classize_m enrollment

* see the difference: 
plot classize enrollment

* alternative 2

twoway scatter classize_m enrollment

* I-(5)
su avgverb

reg avgverb classize, r
reg avgverb classize perc_disadvantaged, r
reg avgverb classize perc_disadvantaged enrollment , r


su avgmath

reg avgmath classize, r
reg avgmath classize perc_disadvantaged, r
reg avgmath classize perc_disadvantaged enrollment , r


* I-(6)
gen fsc=enrollment/(int((enrollment-1)/40)+1)


* I-(7)
gen enrollment2=(enrollment^2)/100
reg avgverb fsc perc_disadvantaged enrollment enrollment2, r
reg avgmath fsc perc_disadvantaged enrollment enrollment2, r


* I-(8)
ivregress 2sls avgverb (classize= fsc) perc_disadvantaged enrollment enrollment2, r
ivregress 2sls avgmath (classize= fsc) perc_disadvantaged enrollment enrollment2, r 


* I-(9)
* Reduced form estimation using `regression' command

reg classize fsc perc_disadvantaged enrollment enrollment2, r

* using `ivregress' command:
ivregress 2sls avgverb (classize= fsc) perc_disadvantaged enrollment enrollment2, r first

* Note the slightly different result bewteen both methods - this is dues to different computation of standard errors and slightly different samples. 
* It is usually better to use the "first" option a sthis relates to the specific 2SLS regression estimated (unless there are many outcomes in which case having many first stage estimates may become intractable).

* I-(10)
* DISCUSS
