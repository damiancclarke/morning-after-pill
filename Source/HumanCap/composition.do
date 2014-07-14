* composition.do v1.00           damiancclarke             yyyy-mm-dd:2014-01-10
*---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
*
/* Takes data from the file BirthsPAE (generated in birthimport.do), and runs
regressions to determine whether the morning after pill has compositional
effects on babies or mothers.

*/

version 11
clear all
cap log close
set more off

********************************************************************************
*** (1) Set globals and locals
********************************************************************************
global PATH "~/universidades/Oxford/DPhil/Thesis/Teens"
global DATA "$PATH/Data/Nacimientos"
global OUT  "$PATH/Tables"
global LOG  "$PATH/Log"

local infile "BirthsPAE"
local outfile "HumanCapital.tex"
local stderror cluster(comcode)

log using $LOG/composition.txt, text replace

local in 1
local mreg 1
local creg 1
local write 0
********************************************************************************
*** (2) Open file, make decisions on unreliable data
********************************************************************************
if `in'==1 {
use "$DATA/`infile'"

encode comuna_names, gen(comcode)
gen trend=year-2005
keep if year>2005

gen educ_level=0 if nivel_m==5
replace educ_level=1 if nivel_m==4
replace educ_level=2 if nivel_m==3|nivel_m==2
replace educ_level=3 if nivel_m==1

gen educ_years = 0 + curso_m if nivel_m==5|nivel_m==4
replace educ_years = 8 + curso_m if nivel_m==3|nivel_m==2
replace educ_years = 12 + curso_m if nivel_m==1

gen married=est_civ_m==1
replace married=. if est_civ_m==.

gen working=1 if activ_m==1
replace working=0 if activ_m==0|activ_m==2

gen bwt=peso if peso>500&peso<5500
gen gestation=semanas if semanas>22&semanas<45
gen length=talla if talla>40&talla<55

gen group=1 if edad_m>14&edad_m<19
replace group=2 if edad_m>18&edad_m<35
replace group=3 if edad_m>34&edad_m<50

gen close15=pilldistance>0&pilldistance<15
gen close30=pilldistance>=15&pilldistance<30
}
********************************************************************************
*** (3a) Maternal Regressions
********************************************************************************
if `mreg'==1 {
foreach g of numlist 1(1)3 {
	local Mbetas`g' 
	local Mses`g'
	local MR2`g' 
	local MObs`g'
	foreach y of varlist /*educ_level*/ educ_years married working {
		xi: reg `y' pill i.trend*comcode i.comcode close* if group==`g', `stderror'
		local beta=round(1000*`=_b[pill]')/1000
		local se=round(1000*`=_se[pill]')/1000	
		local t = abs(_b[pill]/_se[pill])
		local r = round(e(r2)*1000)/1000
		local n = e(N)
		
		local Mbetas`g' `Mbetas`g''  `beta'
		local Mses`g' `Mses`g''  (`se')
		local MR2`g' `MR2`g''  `r'
		local MObs`g' `MObs`g''  `n'		
		
		if `t'>2.58 local Mbetas`g' `Mbetas`g''***
		else if `t'>1.96 local Mbetas`g' `Mbetas`g''**
		else if `t'>1.64 local Mbetas`g' `Mbetas`g''*

		dis "`Mbetas`g''"
		dis "`Mses`g''"
	}
}
}
********************************************************************************
*** (3b) Child Regressions
********************************************************************************
if `creg'==1 {
foreach g of numlist 1(1)3 {
	local Cbetas`g'
	local Cses`g'
	local CR2`g'
	local CObs`g'

	foreach y of varlist peso gestation length {
		xi: reg `y' pill i.trend*comcode i.comcode close* if group==`g', `stderror'
		local beta=round(1000*`=_b[pill]')/1000
		local se=round(1000*`=_se[pill]')/1000	
		local t = abs(_b[pill]/_se[pill])
		local r = round(e(r2)*1000)/1000
		local n = e(N)
		
		local Cbetas`g' `Cbetas`g'' `beta'
		local Cses`g' `Cses`g'' (`se')
		local CR2`g' `CR2`g'' `r'
		local CObs`g' `CObs`g'' `n'		
		
		if `t'>2.58 local Cbetas`g' `Cbetas`g''***
		else if `t'>1.96 local Cbetas`g' `Cbetas`g''**
		else if `t'>1.64 local Cbetas`g' `Cbetas`g''*

		dis "`Cbetas`g''"
		dis "`Cses`g''"
	}
}
}
********************************************************************************
*** (4) Write results out
********************************************************************************
if `write'==1 {
local title Emergency Contraception and Aggregate Human Capital

cap rm "$OUT/`outfile'"
file open  af using "$OUT/`outfile'", write
file write af "\begin{landscape}\begin{table}[htpb!]\centering" _n
file write af "\caption{`title'} \label{TEENtab:PillAgg}" _n
file write af "\begin{tabular}{@{\extracolsep{5pt}}lccccccccc} \\" _n
file write af "[-1.8ex]\hline\hline \\[-1.8ex] &"
file write af "\multicolumn{3}{c}{15-19 year olds} &"
file write af "\multicolumn{3}{c}{20-34 year olds} &"
file write af "\multicolumn{3}{c}{35-49 year olds} \\" _n
file write af "\cmidrule(r){2-4} \cmidrule(r){5-7} \cmidrule(r){8-10}" _n
file write af "\textsc{Panel A:}&(1)&(2)&(3)&(4)&(5)&(6)&(7)&(8)&(9) \\" _n
file write af "\textsc{Mother Characteristics} & Yrs Educ & Working & Married"
file write af "& Yrs Educ & Working & Married & Yrs Educ & Working & Married"
file write af "\\ \midrule" _n
file write af " & & & & & & & & & \\" _n
file write af "Morning After Pill & `Mbetas1'&`Mbetas2'&`Mbetas3' \\" _n
file write af "`Mses1'&`Mses2'&`Mses3' \\" _n
file write af " & & & & & & & & & \\" _n
file write af "Observations & `MObs1'&`MObs2'&`MObs3'\\" _n
file write af "$ R^2 $ & `MR21'&`MR22'&`MR23' \\" _n
file write af " & & & & & & & & & \\" _n
file write af " & & & & & & & & & \\ \midrule" _n	
file write af "\textsc{Panel B:}&(1)&(2)&(3)&(4)&(5)&(6)&(7)&(8)&(9) \\" _n
file write af "\textsc{Child Characteristics} & Weight & Gestation & Length"
file write af "& Weight & Gestation & Length & Weight & Gestation & Length"
file write af "\\ \midrule" _n
file write af " & & & & & & & & & \\" _n
file write af "Morning After Pill & `Cbetas1'&`Cbetas2'&`Cbetas3' \\" _n
file write af "& `Cses1'&`Cses2'&`Cses3' \\" _n
file write af " & & & & & & & & & \\" _n
file write af "Observations & `CObs1'&`CObs2'&`CObs3'\\" _n
file write af "$ R^2 $ & `CR21'&`CR22'&`CR23' \\ \hline \hline \\[-1.8ex]" _n
file write af "\multicolumn{10}{p{21cm}}{\begin{footnotesize}\textsc{Notes:}"
file write af "$^{*}$ p $<0.1$; $^{**}$ p $<0.05$; $^{***}$ p $<0.01$."
file write af "\end{footnotesize}}" _n
file write af "\normalsize\end{tabular}\end{table}\end{landscape}" _n

file close af 
}


********************************************************************************
*** (5) Ad-hoc sum stats which should be incorporated into sumStats.R file
********************************************************************************
drop if pill==.
estpost tabstat bwt educ_years working married edad_m, by(pill) ///
  statistics(mean sd) columns(statistics) listwise

esttab, main(mean) aux(sd) nostar unstack noobs nonote nomtitle nonumber tex
