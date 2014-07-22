# Assessing Plan B
## Bentancor and Clarke, 2014


This repository contains all source code, public data, and the full commit history for the paper "Assessing Plan B: The Effect of the Morning After Pill on Children and Women".  Below details are provided regarding replicating tables, accessing data, and running the source code to generate the entire project.   For questions or requests, contact <damian.clarke@economics.ox.ac.uk>, or fork this repository.  If you are interested in this code but are not familiar with github, source files can be downloaded individually on the [Source page](https://github.com/damiancclarke/morning-after-pill/tree/master/Source) of this site.  For further details regarding git, see <http://git-scm.com/book>.

***

## General

This paper is based upon a number of data sources and, a number of scripts written in R, Python, Stata, LaTeX and make.  The public repository for this project contains all raw data where this data isn't too large to reasonably be quickly loaded on to github.  In the case of large data files, I make links available in the Data section of this wiki, or can send upon request.

## Code

All source code for data generation and estimates is provided in the following folders:
* `./Source/Births/` : Code to generate birth data and estimate tables 2,5,7,9 and figures 5-6.
* `./Source/Deaths/` : Code to generate fetal death data and estimate tables 3,5 and 8.
* `./Source/HumanCap/` : Code to generate mother and child characteristic data and estimate table 4.
* `./Source/placebo/` : Code to generate placebo data (birth and deaths) and estimate table 6.
* `./Source/Summary/` : Code to generate summary statistics table 1, and figures 1,2,3,7 and 8.

Each file automatically exports figures and tables as .eps files (figures), and .tex files (tables) into the folders `./Figures` and `./Tables` respectively.  These can be run one-by-one in R if individual tables and images are desired.  In each case, the only change required will be replacing one line in each file which defines `proj.dir`, the location of all files.


## Data Sources
The data sources are listed below, and are available in the SData/ folder:

* Population/SalComUsuarios-REGIONTok_1_ESAcal.csv
  > Number of inhabitants per Region of Chile by year 
  > where REGION should be replace with 01, 02, ..., 15 
  > Full data can be found online (01/01/2013) at:  
    http://www.ine.cl/canales/chile_estadistico/familias/demograficas_vitales.php 
  > Original xls files are also provided in the Population/. folder
* Births/Nacimientos_Chile_20002011.dta
  > Censal data of all births in Chile between 2000 and 2011
  > Data downloaded from the Ministry of Health's DEIS website (01/01/2013):
    http://www.deis.cl/descargar-bases-de-datos/
* Geographic/regionescomunas_short.csv
  > Official names of Chile's Regions and comunas
* Contraceptive/PillDist.csv
  > Availability of Pill by region.  Coded from Dides et al 2009, 2010, 2011
  > Distance generated as polygon centrepoint to centrepoint
* Covariates/MayorCharacteristics_years.csv
  > characterisitcs of the Mayor and political party by comuna and year
* Covariates/Comunas_datos_20062012.csv
  > characteristics such as educ spending, health spending, female workers
  > compiled from SINIM's 'ficha comunal'.  Online (01/01/2013):
    http://www.sinim.gov.cl/indicadores/busq_serie.php
* Deaths/DEFyyyy.csv
  > Censal data on fetal deaths in a year
  > where yyyy should be replaced by 2002, 2003, ..., 2011	
  > Data downloaded from the Ministry of Health's DEIS website (01/01/2013):
    http://www.deis.cl/descargar-bases-de-datos/

## Source Code
The source code files are listed below, and are available in the SSource/ folder:

* Summary/SumStats.R
  > Produce summary stats.  Generates table 1 and figures 1-4 
* Births/BirthGenerate.R
  > Generate one line per comuna, year and birth status
* Births/BirthEstimates.R
  > Estimate effects of pill on births. Generate tables 2,5,6,8, figures 5,6
* Deaths/DeathGenerate.R
  > Generate one line per comuna, year and death status
* Deaths/DeathEstimates.R
  > Estimate effects of pill on fetal deaths. Generate tables 3,5,7
* HumanCap/birthimport.do
  > Creates one line per mother (and child) with their characteristics
* HumanCap/composition.do
  > Tests whether there are large-scale composition effects on birth cohorts
* Tests/PlaceboGenerate.R
* Tests/PlaceboEstimates.R

## Working Data
The final datasets produced for analysis are:

* Working/BirthsPAE.dta 
  > characteristics of mothers/children.  Zipped due to size
* Working/S1Data_granular_covars.csv 
  > one line per muncipality, woman, birth status
* Working/S1Data_deaths_covars.csv 
  > one line per muncipality, woman, death status

