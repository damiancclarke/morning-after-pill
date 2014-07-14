* birthimport.do v1.00           damiancclarke             yyyy-mm-dd:2013-11-22 
*---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
*

clear all
cap log close
set more off
version 11
set matsize 5000

********************************************************************************
*** (1) Globals and Locals
********************************************************************************
global PATH "~/universidades/Oxford/DPhil/Thesis/Teens"
global DATA "$PATH/Data/MinSal/dta"
global RESULTS "$PATH/Results/Births"
global DATAOUT "$PATH/Data/Nacimientos"
global COMUNAS "~/database/ChileRegiones/Nombres"
global PAE "$PATH/Data/PAE"
global POPULATION "$PATH/Data/Poblacion"
global LOG "$PATH/Log"

local births "Nacimientos_Chile_20002011.dta"

local birthcomunas yes
local popln no
local regs no
local regs2 no

log using $LOG/birthimport.txt, text replace
********************************************************************************
*** (2) Import births, merge in Comuna lables from official comuna file
********************************************************************************
if `"`birthcomunas'"'=="yes" {
	use "$DATA/`births'"
	keep mes_nac ano_nac dom_comuna semanas peso talla comuna hij_vivos /*
	*/ hij_total edad_m nivel_m urb_rural curso_m categ_m est_civ_m activ_m

	keep if ano_nac>2000
	
	gen mage=0014 if edad_m<15
	replace mage=1519 if edad_m>=15&edad_m<20
	replace mage=2024 if edad_m>=20&edad_m<25
	replace mage=2529 if edad_m>=25&edad_m<30
	replace mage=3034 if edad_m>=30&edad_m<35
	replace mage=3539 if edad_m>=35&edad_m<40
	replace mage=4044 if edad_m>=40&edad_m<45
	replace mage=4549 if edad_m>=45&edad_m<50
	replace mage=5000 if edad_m>=50

	drop if mage==0014|mage==5000
	
	gen bord=1 if hij_total<=1
	replace bord=2 if hij_total==2
	replace bord=3 if hij_total>2
		
 	replace comuna=5801  if ano_nac<2010 & comuna==5106	
	replace comuna=5804  if ano_nac<2010 & comuna==5108
	replace comuna=5802  if ano_nac<2010 & comuna==5505
 	replace comuna=5803  if ano_nac<2010 & comuna==5507	
 	replace comuna=15101 if ano_nac<2008 & comuna==1201
	replace comuna=15102 if ano_nac<2008 & comuna==1202
	replace comuna=15201 if ano_nac<2008 & comuna==1301
	replace comuna=15202 if ano_nac<2008 & comuna==1302
	replace comuna=1402  if ano_nac<2008 & comuna==1102	
 	replace comuna=1403  if ano_nac<2008 & comuna==1103
	replace comuna=1404  if ano_nac<2008 & comuna==1104
	replace comuna=1405  if ano_nac<2008 & comuna==1105
	replace comuna=1401  if ano_nac<2008 & comuna==1106
	replace comuna=14101 if ano_nac<2008 & comuna==10501
	replace comuna=14102 if ano_nac<2008 & comuna==10502
	replace comuna=14202 if ano_nac<2008 & comuna==10503
	replace comuna=14201 if ano_nac<2008 & comuna==10504
	replace comuna=14203 if ano_nac<2008 & comuna==10505
	replace comuna=14103 if ano_nac<2008 & comuna==10506
	replace comuna=14104 if ano_nac<2008 & comuna==10507
	replace comuna=14105 if ano_nac<2008 & comuna==10508
	replace comuna=14106 if ano_nac<2008 & comuna==10509
	replace comuna=14107 if ano_nac<2008 & comuna==10510
	replace comuna=14108 if ano_nac<2008 & comuna==10511
	replace comuna=14204 if ano_nac<2008 & comuna==10512
	
	rename comuna comunacode2010
	merge m:1 comuna using "$COMUNAS/comunacodesshort.dta", gen(_mergeComuna)
	drop _mergeComuna

	
   *****************************************************************************
   *** (3) Merge with PAE
	*** Here year is the pill year, so from now on 2010 refers to those people
	*** who had the pill in 2010, and so gave birth in 2011.
   *****************************************************************************
	gen year=ano_nac-1

	gen comuna_alt=lower(comuna)
	replace comuna_alt=subinstr(comuna_alt, "á", "a", .)
	replace comuna_alt=subinstr(comuna_alt, "é", "e", .)
	replace comuna_alt=subinstr(comuna_alt, "í", "i", .)
	replace comuna_alt=subinstr(comuna_alt, "ó", "o", .)
	replace comuna_alt=subinstr(comuna_alt, "ú", "u", .)
	replace comuna_alt=subinstr(comuna_alt, "ü", "u", .)	
	replace comuna_alt=upper(comuna_alt)
	replace comuna_alt=subinstr(comuna_alt, "ñ", "Ñ", .)
	replace comuna_alt=subinstr(comuna_alt, "Á", "A", .)
	replace comuna_alt=subinstr(comuna_alt, "É", "E", .)
	replace comuna_alt=subinstr(comuna_alt, "Í", "I", .)
	replace comuna_alt=subinstr(comuna_alt, "Ó", "O", .)
	replace comuna_alt=subinstr(comuna_alt, "Ú", "U", .)
	replace comuna_alt=subinstr(comuna_alt, "Ü", "U", .)
	replace comuna_alt="OLLAGÜE" if comuna_alt=="OLLAGUE"
	replace comuna_alt="ISLA DE PASCUA" if regexm(comuna_alt, "PASCUA")
	replace comuna_alt="PEDRO AGUIRRE CERDA" if regexm(comuna_alt, "CERDA")
	replace comuna_alt="CURAHUE" if comuna_alt=="CARAHUE"
	replace comuna_alt="PUERTO NATALES" if comuna_alt=="NATALES"
	replace comuna_alt="PUERTO SAAVEDRA" if comuna_alt=="SAAVEDRA"
	replace comuna_alt="VILLARICA" if comuna_alt=="VILLARRICA"
	replace comuna_alt="ANTARTICA" if comuna_alt=="TORRES DEL PAINE"
	gen comuna_names=comuna_alt
	
	merge m:m comuna_names year using $PAE/PillDist_noreps, gen(_merPAE)
	tab year if _merPAE==2

	replace pill=0 if year<2009
		
	save $DATAOUT/BirthsPAE, replace
}


cap log close
