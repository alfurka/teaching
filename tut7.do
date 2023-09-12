************************************************************************
* Practical session 6: Stata Applications of advanced panel data method
************************************************************************

****************************************************************
* Problem I: the deterrent effects of execution on murder
****************************************************************
cscript
use "murder.dta", clear

* I-(1)



* I-(2)
reg mrdrte d93 exec unem if year==90|year==93

* I-(3)


xtset id year


xtreg mrdrte d93 exec unem if year==90|year==93, fe

sum mrdrte if year==93


* I-(4) 
keep if year==90|year==93
replace year=1 if year==90
replace year=2 if year==93
xtset id year

* "D." allows to use first differenced variables
reg D.(mrdrte d93 exec unem), nocon vce(robust) 

* I-(5) 
sort state
by state: sum exec if year==2

* find the outlier: 

* I-(6) 
xtset id year
reg D.(mrdrte d93 exec unem) if state!="TX", nocon
reg D.(mrdrte d93 exec unem) if state!="TX", nocon  vce(robust)


* I-(7) 

* we deleted some observations. So, we will reload the data set to stata:

cscript

use "murder.dta", clear
xtset id year
xtreg mrdrte d90 d93 exec unem, fe vce(robust)
xtreg mrdrte d90 d93 exec unem, fe vce(cluster id)
* vce(robust) produces standard errors which are robust to heteroskedasticity
* vce(cluster cross-section dimension) produces standard errors which are robust to heteroskedasticity and serial correlation




*********************************************************************
* Problem II: A Study of the effect of small class with Project STAR
*********************************************************************

cscript
use "starwide.dta", clear


* II-(1)
*small0: if they were in a small class in kindergarden
*small1: if they were in a small class in grade 1
*small2: if they were in a small class in grade 2
*small3: if they were in a small class in grade 3
*aide0: if they were in a aide class in kindergarden
*aide1: if they were in a aide class in grade 1
*aide2: if they were in a aide class in grade 2
*aide3: if they were in a aide class in grade 3
gen small_sum0=small0
egen small_sum1=rowtotal(small0 small1) /* rowtotal treats missings as 0 */
egen small_sum2=rowtotal(small0 small1 small2)
egen small_sum3=rowtotal(small0 small1 small2 small3)
gen aide_sum0=aide0
egen aide_sum1=rowtotal(aide0 aide1)
egen aide_sum2=rowtotal(aide0 aide1 aide2)
egen aide_sum3=rowtotal(aide0 aide1 aide2 aide3)


* II-(2)

* see how data set looks like: 

bro

* how to reshape a file? 
* synthax: reshape long xVar yVar zVar ..., i(id) j(time)

reshape long small_sum aide_sum female wh_asn small regular aide frlunch inismall iniaide twhite tmaster tladder avesat schid tyears, i(stdntid) j(grade)

* See the long format: 

bro

*define variables in the new format
label variable grade "student grade 0-kinder, 1-1st grade, 2-2nd grade, 3-3rd grade"
label variable small "class type small"
label variable regular "class type regular"
label variable aide "class type regular with teacher aide(part time)"
label variable frlunch "free lunch status"
label variable inismall "initial assignment of class type small"
label variable iniaide "initial assignment of class type regular class with aide"
label variable twhite "teacher race =white"
label variable tmaster "teacher's highest degree, dummy=1 if master or above 0 otherwise"
label variable tladder "teacher career ladder level"
label variable avesat "SAT score, average over three subjects"

* see variables:

de

* II-(3)
xtset stdntid grade



* II- (4)

* I will go from simple to complex for the first plot:
*1: kdensity
kdensity avesat if small==1&grade==0
*2: lwidth
kdensity avesat if small==1&grade==0, lwidth(thick)
*3: addplot
kdensity avesat if small==1&grade==0, lwidth(thick) addplot(kdensity avesat if regular==1&grade==0, lwidth(thick))
* 4: legend, label and title
kdensity avesat if small==1&grade==0, lwidth(thick) addplot(kdensity avesat if regular==1&grade==0, lwidth(thick)) legend(label(1 "small") label(2 "regular")) title("Kindergarten")



* all are here:

kdensity avesat if small==1&grade==0, lwidth(thick) addplot(kdensity avesat if regular==1&grade==0, lwidth(thick) ) legend(label(1 "small") label(2 "regular")) title("Kindergarten")
kdensity avesat if small==1&grade==1, lwidth(thick) addplot(kdensity avesat if regular==1&grade==1, lwidth(thick) ) legend(label(1 "small") label(2 "regular")) title("1st Grade")
kdensity avesat if small==1&grade==2, lwidth(thick) addplot(kdensity avesat if regular==1&grade==2, lwidth(thick) ) legend(label(1 "small") label(2 "regular")) title("2nd Grade")
kdensity avesat if small==1&grade==3, lwidth(thick) addplot(kdensity avesat if regular==1&grade==3, lwidth(thick) ) legend(label(1 "small") label(2 "regular")) title("3rd Grade") 


