********************************************************************************
README.txt for:
Assessing Plan B: The Effect of the Morning After Pill on Children and Women

contact: mailto:damian.clarke@economics.ox.ac.uk
********************************************************************************

This paper is based upon a number of data sources and, a number of scripts 
written in R, Python, Stata and LaTeX.  The public repository for this project
contains all raw data after the most crude cleaning has been carried out.  For
example, for the population data we don't include the entire meta-data, but 
rather include the region-specific sheets which are imported into R when 
generating our main database.  However, in case of interest, I am more than 
happy to share any and all completely raw data.  Contact me with queries.

The data sources are:
> Population/SalComUsuarios-REGIONTok_1_ESAcal.csv
	- Number of inhabitants per Region of Chile by year
	- where REGION should be replace with 01, 02, ..., 15
	- Full data can be found online (01/01/2013) at: 
		http://www.ine.cl/canales/chile_estadistico/familias/demograficas_vitales.php
	- Original xls files are also provided in the Population/. folder
> Births/Nacimientos_Chile_20002011.dta
	- Censal data of all births in Chile between 2000 and 2011
	- Data downloaded from the Ministry of Health's DEIS website (01/01/2013):
	http://www.deis.cl/descargar-bases-de-datos/
> Geographic/regionescomunas_short.csv
	- Official names of Chile's Regions and comunas
> Contraceptive/PillDist.csv
	- Availability of Pill by region.  Coded from Dides et al 2009, 2010, 2011
	- Distance generated as polygon centrepoint to centrepoint
> Covariates/MayorCharacteristics_years.csv
	- characterisitcs of the Mayor and political party by comuna and year
> Covariates/Comunas_datos_20062012.csv
	- characteristics such as educ spending, health spending, female workers
	- compiled from SINIM's 'ficha comunal'.  Online (01/01/2013):
		http://www.sinim.gov.cl/indicadores/busq_serie.php
> Deaths/DEFyyyy.csv
	- Censal data on fetal deaths in a year
	- where yyyy should be replaced by 2002, 2003, ..., 2011	
	- Data downloaded from the Ministry of Health's DEIS website (01/01/2013):
		http://www.deis.cl/descargar-bases-de-datos/

The source code scripts are:
> Summary/SumStats.R
	- Produce summary stats.  Generates table 1 and figures 1-4 
> Births/BirthGenerate.R
	- Generate one line per comuna, year and birth status
> Births/BirthEstimates.R
	- Estimate effects of pill on births. Generate tables 2,5,6,8, figures 5,6
> Deaths/DeathGenerate.R
	- Generate one line per comuna, year and death status
> Deaths/DeathEstimates.R
	- Estimate effects of pill on fetal deaths. Generate tables 3,5,7
> HumanCap/birthimport.do
	- Creates one line per mother (and child) with their characteristics
> HumanCap/composition.do
	- Tests whether there are large-scale composition effects on birth cohorts
> Tests/PlaceboGenerate.R
> Tests/PlaceboEstimates.R

And the final data produced is:
> Working/BirthsPAE.dta 
		- characteristics of mothers/children.  Zipped due to size
> Working/S1Data_granular_covars.csv 
		- one line per muncipality, woman, birth status
> Working/S1Data_deaths_covars.csv 
		- one line per muncipality, woman, death status
