****************************************************************
* Problem I: The deterrent effects of smoking on income
****************************************************************
cscript
use "smoke.dta", clear

* I-(1)
de

* I-(2)
reg lincome cigs educ age agesq, r


* I-(3)

*>>>> Discuss

* I-(4)

*>>>> Discuss

* I-(5) 

reg cigs educ age agesq lcigpric restaurn, r


test lcigpric restaurn

* I-(6) 

ivregress 2sls lincome educ age agesq (cigs=lcigpric restaurn), r

reg cigs educ age agesq restaurn, r
test restaurn


* I-(7) 

*>>>> Discuss



*********************************************************
* Problem II: Demand and Supply for fish
*********************************************************


cscript
* Loading data
use "fish.dta", clear

* II-(1)

*>>>> Discuss

* II-(2)

*>>>> Discuss


* II-(3)

reg lavgprc wave2 wave3, r
test wave2 wave3

* II-(4)


ivregress 2sls ltotqty (lavgprc= wave2 wave3), r


* II-(5)

*>>>> Discuss

* II-(6)


reg ltotqty mon tues wed thurs, r
test mon tues wed thurs



reg ltotqty tues wed, r
test  tues wed 
