# Assessing Plan B
## Bentancor and Clarke, 2014


This repository contains all source code, public data, and the full commit history for the paper "Assessing Plan B: The Effect of the Morning After Pill on Children and Women".  Below details are provided regarding replicating tables, accessing data, and running the source code to generate the entire project.   For questions or requests, contact <damian.clarke@economics.ox.ac.uk>, or fork this repository.  If you are interested in this code but are not familiar with github, source files can be downloaded individually on the [Source page](https://github.com/damiancclarke/morning-after-pill/tree/master/Source) of this site.  For further details regarding git, see <http://git-scm.com/book>.

***

## General

This paper is based upon a number of data sources and, a number of scripts written in R, Python, Stata, LaTeX and make.  The public repository for this project contains all raw data where this data isn't too large to reasonably be quickly loaded on to github.  In the case of large data files, I make links available in the Data section of this wiki, or can send upon request.

## Code
### Source Files

All source code for data generation and estimates is provided in the following folders:
* `./Source/Births/` : Code to generate birth data and estimate tables 2,5,7,9 and figures 5-6.
* `./Source/Deaths/` : Code to generate fetal death data and estimate tables 3,5 and 8.
* `./Source/HumanCap/` : Code to generate mother/child characteristic data and estimate tab 4.
* `./Source/placebo/` : Code to generate placebo data (birth and deaths) and estimate table 6.
* `./Source/Summary/` : Code to generate summary statistics table 1, and figures 1,2,3,7 and 8.

Each file automatically exports figures and tables as .eps files (figures), and .tex files (tables) into the folders `./Figures` and `./Tables` respectively.  These can be run one-by-one in R if individual tables and images are desired.  In each case, the only change required will be replacing one line in each file which defines `proj.dir`, the location of all files.

### Running All Files
Alternatively, the entire paper, all estimates, and all data creation can be run which requires only initial raw data and source files.  This can be done by using the Makefile in the `./Paper` folder.  On Unix and Mac make will be automatically installed and so this file should run without additional software.  Windows users may need to install [GnuWin](http://gnuwin32.sourceforge.net/).  To run this Makefile, simply type `make` at the command line.

## Data
### Raw Data
Raw data comes from a number of sources.  These raw data files are available in the [./Data](https://github.com/damiancclarke/morning-after-pill/tree/master/Data) folder of this page where possible, and when too large to house here can be accessed (as at July 21, 2014) at the links below.  All raw data is also available by contacting me. 

* `./Data/Population/SalComUsuarios-REGIONTok_1_ESAcal.csv`
  * Number of inhabitants per Region of Chile by year where REGION should be replace with 01, 02, ..., 15 
  * Full data can be found online [here](http://www.ine.cl/canales/chile_estadistico/familias/demograficas_vitales.php)
* `./Data/Nacimientos/Nacimientos_Chile_20002011.dta`
  * Censal data of all births in Chile between 2000 and 2011
  * Data downloaded from the Ministry of Health's [DEIS website](http://www.deis.cl/descargar-bases-de-datos/) (01/01/2013)
* `./Data/Comunas/regionescomunas_short.csv`
  * Official names of Chile's Regions and comunas
* `./Data/Comunas/Comunas_datos_20062012.csv`
  * characteristics such as educ spending, health spending, female workers
  * compiled from SINIM's 'ficha comunal'.  [Online](http://www.sinim.gov.cl/indicadores/busq_serie.php)
* `./Data/PAE/PillDist.csv`
  * Availability of Pill by region.  Coded from Dides et al 2009, 2010, 2011.
  * Distance generated as polygon centrepoint to centrepoint
* `./Data/Alcaldes/MayorCharacteristics_years.csv`
  * characterisitcs of the Mayor and political party by comuna and year
* `./Data/Deaths/DEFyyyy.csv`
  * Censal data on fetal deaths in a year, where yyyy should be replaced by 2002, 2003, ..., 2011	
  * Data downloaded from the Ministry of Health's [DEIS website](http://www.deis.cl/descargar-bases-de-datos/)

### Interim Data
Rather than compiling all data from source, interested users can also download interim (processed) datasets.  Full details of processing can be found in the processing functions BirthGenerate.R and DeathGenerate.R.  A description of interim data along with its location is provided below.

* `.Data/Nacimientos/S1_Data_granular_covars.csv`
  * Contains one line per municipality*year*pregnancy status with the number of occurrences
* `.Data/Nacimientos/S1_Data_covars_20002011.csv`
  * Same as above, but for entire period 2000-2011.  Used for placebo tests.
* `.Data/Deaths/S1Data_deaths_covars.csv`
  * Contains one line per municipality*year*pregnancy*death status with the number of occurrences
* `.Data/Deaths/S1Data_20002011.csv`
  * Same as above, but for entire period 2000-2011.  Used for placebo tests.
