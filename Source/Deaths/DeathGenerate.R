# DeathGenerate.R v1.00            KEL / DCC               yyyy-mm-dd:2013-12-24
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# Generates the file described in "filename". This file has one line per comuna,
# age group and pregnancy status.  For the total number of pregnant women, there
# is also a measure of the number of pregnancies terminating in fetal death. The
# final file also contains covariates for political controls, education and hea-
# lth controls, gender controls, and pill status (and distance to nearest comuna
# with pill if pill==0).
#
# The file accepts the arguments age_range (generally 15-49) weeks (the number
# of weeks gestation the death occurs before), and regex (a regex style pattern
# to use for searching ICD codes).  If no ICD code is desired, * should simply
# be entered as this argument.  If the full file with comuna data is desired the
# usecom argument should be set to TRUE, and the name of the file to be exported
# is listed as filename.
#
# This is a minor refactorization of managedata_deaths_KEL_20131204.R. This 
# code has been written by KEL, with minor additions by DCC to incorporate 
# additional time varying controls and pill distance data.
#
# prior to running this function the following directories must be defined:
#   > pop.dir (contains all population files from INE converted to .csv)
#   > brth.dir (directory containing data file of all births in Chile)
#   > deth.dir (directory containing data file of all fetal deaths in Chile)
#   > ma.dir   (directory where pill data is stored)
#   > pol.dir  (directory where mayor characteristics are kept)
#   > com.dir  (directory where comuna characteristics and names are kept)
#
# NOTE: FIX UP PART OF CODE WHICH CREATES EARLY/LATE DUMMIES.  THIS IS
# VERY VERBOSE - JUST MAKE INTO FUNCTION.

prep_s1_data_deaths <- function(age_range,week,regex,usecom,filename) {
  
  #*****************************************************************************
  # (1) libraries
  #*****************************************************************************
  require(foreign)
  require(reshape)
  require(gdata)
  
  #*****************************************************************************
  # (2) comuna rename function
  #*****************************************************************************  
  remove_accents <- function(string) {
    string <- trim(string)
    string <- sub("Á","A",string)
    string <- sub("É","E",string)
    string <- sub("Í","I",string)
    string <- sub("Ó","O",string)
    string <- sub("Ú","U",string)  
    return(string)  
  }
  
  #*****************************************************************************
  # (3) extract csv sheets from excel population data
  #*****************************************************************************
  if(F) {
    #system('python ~universidades/Oxford/DPhil/Thesis/Teens/Source/xls2csv.py')
  }
  
  #*****************************************************************************
  # (4) Import and format popluation data
  #*****************************************************************************
  nums <- c(paste("0",1:9,sep=""),10:15)
  f.vec <- paste("SalComUsuarios-",nums, "Tok_1_ESAcal.csv",sep="")
  f <- f.vec[1]
  
  fix_missing <- function(vec) {
    notna <- vec[!(is.na(vec))]
    notna[as.character(notna) %in% c("", " -", "-")] <- NA
    vec[!(is.na(vec))] <- notna
    return(vec)
  }
  
  clean_pop_csv <- function(f) {
    
    tmp <- read.csv(paste(pop.dir,f,sep=""))
    names(tmp) <- paste("V",1:(dim(tmp)[2]),sep="")
    
    tmp <- apply(tmp,2,fix_missing)
    
    isblankcol <- apply(tmp,2,function(vec) {sum(is.na(vec)) == length(vec)})
    if(sum(isblankcol) > 1) print("Error -- too many blank columns")
    tmp1 <- tmp[,1:(which(isblankcol)-1)]
    tmp <- tmp[,(which(isblankcol)+1):dim(tmp)[2]]
    
    tmp <- rbind(c("INICIO",rep(NA,dim(tmp)[2]-1)),tmp)
    isnewtab <- which(tmp[,1] == "INICIO")
    comunas <- (as.character(tmp[isnewtab+1,1]))
    comunas <- lapply(as.list(comunas),
                      function(c) rep(c,length((age_range[1]):(age_range[2]))))
    comunas <- do.call("c",comunas)
    year <- tmp[4,]
    tmp2 <- tmp[as.character(tmp[,1]) %in% (age_range[1]):(age_range[2]),]
    tmp2 <- tmp2[!is.na(tmp2[,1]),]
    tmp2 <- cbind(comunas,tmp2)
    rownames(tmp2) <- NULL
    colnames(tmp2) <- c("dom_comuna","age",paste("d",1990:2020,sep=""))
    tmp2 <- as.data.frame(tmp2)
    tmp2 <- tmp2[,c("dom_comuna","age",paste("d",birth_y_range,sep=""))]
    tmp3 <- melt.data.frame(tmp2,id.vars=c("dom_comuna","age"),
                            measure.vars = paste("d",birth_y_range,sep=""),
                            variable_name="year")
    tmp3$year <- as.numeric(sub("d","",tmp3$year))
    tmp3$denom <- tmp3$value
    tmp3$value <- NULL
    return(tmp3)
  }
  
  #Run over each region and stack into flat file
  tab.lst <- lapply(f.vec,clean_pop_csv) 
  tot <- do.call("rbind",tab.lst)
  tot <- as.data.frame(tot)
  
  #*****************************************************************************
  # (5) Rename comuna for merging with pill data
  #*****************************************************************************  
  tot$dom_comuna <- apply(as.matrix(tot$dom_comuna),2,remove_accents)
  
  tot$dom_comuna[tot$dom_comuna == "T AMARILLA"] <- "TIERRA AMARILLA"
  tot$dom_comuna[tot$dom_comuna == "D DE ALMAGRO"] <- "DIEGO DE ALMAGRO"
  tot$dom_comuna[tot$dom_comuna == "A DEL CARMEN"] <- "ALTO DEL CARMEN"
  tot$dom_comuna[tot$dom_comuna == "ISLA DE PASCUA"] <- "ISLA  DE PASCUA"
  tot$dom_comuna[tot$dom_comuna == "Q. TILCOCO"] <- "QUINTA DE TILCOCO"
  tot$dom_comuna[tot$dom_comuna == "MARCHIGÜE"] <- "MARCHIHUE"
  tot$dom_comuna[tot$dom_comuna == "PEDRO AGUIRRE CERDA"] <- "PEDRO AGUIRRE CERDA"
  
  #*****************************************************************************
  # (6) Birth data for denominator of probability model
  #*****************************************************************************  
  f <- paste(brth.dir,"Nacimientos_Chile_20002012.dta",sep="")
  tmp <- read.dta(f)

  tmp <- tmp[tmp$edad_m <= age_range[2],]
  tmp <- tmp[tmp$edad_m >= age_range[1],]
  tmp <- tmp[tmp$ano_nac %in% birth_y_range,]
  tmp <- tmp[,c("comuna","ano_nac","edad_m","hij_total")]
  tmp$hij_total[tmp$hij_total > 4] <- 4
  tmp$hij_total[tmp$hij_total < 1] <- 1
  
  #aggregate to comuna
  
  tmp <- as.data.frame(table(tmp$comuna,tmp$ano_nac,tmp$edad_m,tmp$hij_total))
  tmp <- tmp[tmp$Freq > 0,]
  names(tmp) <- c("comuna","ano_nac","edad_m","hij_total","n")
  birth <- tmp
  
  #*****************************************************************************
  # (7) Read in death data
  #*****************************************************************************  
  read_death_data <- function(year) {
    f <- paste(deth.dir,"DEF",year,".csv",sep="")
    tmp <- read.csv(f)

    tmp <- tmp[tmp$EDAD_M <= age_range[2],]
    tmp <- tmp[tmp$EDAD_M >= age_range[1],]
    tmp <- tmp[,c("COMUNA","EDAD_M","HIJ_TOTAL","GESTACION","DIAG1")]
    names(tmp) <- c("comuna","edad_m","hij_total","weeks","code")
    tmp$hij_total[tmp$hij_total > 4] <- 4
    tmp$hij_total[tmp$hij_total < 1] <- 1
    tmp$ano_def <- year
    return(tmp)
  }
  
  lst <- lapply(2002:2012,read_death_data)
  tmp <- do.call("rbind",lst)
  
  #aggregate to comuna
  
  tmp.eA <- tmp[tmp$weeks <= week& grepl(regex,tmp$code),]
  tmp.lA <- tmp[tmp$weeks >  week& grepl(regex,tmp$code),]
  tmp.eB <- tmp[tmp$weeks <= week& !grepl(regex,tmp$code),]
  tmp.lB <- tmp[tmp$weeks >  week& !grepl(regex,tmp$code),]
  
  tmp.all <- as.data.frame(table(tmp$comuna,tmp$ano_def,tmp$edad_m,
                                 tmp$hij_total))
  tmp.eA  <- as.data.frame(table(tmp.eA$comuna,tmp.eA$ano_def,
                                 tmp.eA$edad_m,tmp.eA$hij_total))
  tmp.eB  <- as.data.frame(table(tmp.eB$comuna,tmp.eB$ano_def,
                                 tmp.eB$edad_m,tmp.eB$hij_total))
  tmp.lA  <- as.data.frame(table(tmp.lA$comuna,tmp.lA$ano_def,
                                 tmp.lA$edad_m,tmp.lA$hij_total))
  tmp.lB  <- as.data.frame(table(tmp.lB$comuna,tmp.lB$ano_def,
                                 tmp.lB$edad_m,tmp.lB$hij_total))
  
  names(tmp.all) <-c("comuna","ano_def","edad_m","hij_total","Freq")
  names(tmp.eA)  <-c("comuna","ano_def","edad_m","hij_total","earlyP")
  names(tmp.eB)  <-c("comuna","ano_def","edad_m","hij_total","earlyQ")
  names(tmp.lA)  <-c("comuna","ano_def","edad_m","hij_total","lateP")
  names(tmp.lB)  <-c("comuna","ano_def","edad_m","hij_total","lateQ")  

  tmp<-merge(tmp.all,tmp.eA[,c("comuna","edad_m","hij_total","earlyP","ano_def")],
        by=c("comuna","edad_m","hij_total","ano_def"),all=T)
  tmp<-merge(tmp,tmp.eB[,c("comuna","edad_m","hij_total","earlyQ","ano_def")],
             by=c("comuna","edad_m","hij_total","ano_def"),all=T)
  tmp<-merge(tmp,tmp.lA[,c("comuna","edad_m","hij_total","lateP","ano_def")],
             by=c("comuna","edad_m","hij_total","ano_def"),all=T)
  tmp<-merge(tmp,tmp.lB[,c("comuna","edad_m","hij_total","lateQ","ano_def")],
             by=c("comuna","edad_m","hij_total","ano_def"),all=T)
  
  tmp <- tmp[tmp$Freq > 0,]
  tmp$earlyP[is.na(tmp$earlyP)] <- 0
  tmp$earlyQ[is.na(tmp$earlyQ)] <- 0
  tmp$lateP[is.na(tmp$lateP)] <- 0
  tmp$lateQ[is.na(tmp$lateQ)] <- 0
  
  names(tmp) <- c("comuna","edad_m","hij_total","ano_nac","death","earlyP",
                  "earlyQ","lateP","lateQ")
  
  rm(tmp.all,tmp.eA,tmp.eB,tmp.lA,tmp.lB)
  #*****************************************************************************
  # (8) Merge birth and death data
  #*****************************************************************************  
  tmp$comuna <- as.numeric(as.character(tmp$comuna))
  birth$comuna <- as.numeric(as.character(birth$comuna))
  tmp <- merge(birth,tmp,by=c("comuna","ano_nac","edad_m","hij_total"),
               all.x=T,all.y=F)
  tmp$death[is.na(tmp$death)]   <- 0
  tmp$earlyQ[is.na(tmp$earlyQ)] <- 0
  tmp$earlyP[is.na(tmp$earlyP)] <- 0
  tmp$lateQ[is.na(tmp$lateQ)]   <- 0
  tmp$lateP[is.na(tmp$lateP)]   <- 0
  
  #*****************************************************************************
  # (9) Standardise comuna names
  #*****************************************************************************
  f <- paste(com.dir, "regionescomunas_short.csv",sep="")
  map <- read.csv(f,sep=";")
  map <- map[,1:7]
  map$dom_comuna <- toupper(as.character(map$COMUNA))
  map$dom_comuna <- apply(as.matrix(map$dom_comuna),2,remove_accents)
  map$comuna <- as.numeric(as.character(map$COMUNA.CODE..2010..))
  map$frommap <- 1
  tmp$fromdat <- 1
  tmp$dom_comuna <- NULL
  tmp <- merge(tmp,map[,c("comuna","frommap","dom_comuna")],by="comuna",all.x=T,all.y=F)
  tmp$frommap[is.na(tmp$frommap)] <- 0
  tmp$fromdat[is.na(tmp$fromdat)] <- 0
  tmp$matched <- tmp$frommap & tmp$fromdat
  #Codes that match 2000-2008
  tmp08 <- tmp[!(tmp$matched),]
  map$comuna <- as.numeric(as.character(map$COMUNA.CODE..2000.2008.))
  tmp08$dom_comuna <-  NULL
  tmp08$frommap <- NULL
  tmp08$fromdat <- 1
  tmp08 <- merge(tmp08,map[,c("comuna","dom_comuna","frommap")],by="comuna",
                 all.x=T,all.y=F)
  tmp08$matched <- tmp08$frommap & tmp08$fromdat
  tmp08 <- tmp08[tmp08$matched,]
  tmp08 <- tmp08[!(is.na(tmp08$matched)),]
  
  #clean up
  tmp08$matched <- tmp08$fromdat <- tmp08$frommap <- NULL
  tmp$matched <- tmp$fromdat <- tmp$frommap <- NULL
  
  #combine
  deleteit <- c(tmp08$comuna)
  tmp <- tmp[!(tmp$comuna %in% deleteit),]
  tmp <- rbind(tmp,tmp08)
  tmp <- tmp[!(is.na(tmp$dom_comuna)),]
  tmp$comuna <- NULL
  
  rm(tmp08)
  
  #*****************************************************************************
  # (10) Merge births with total population
  #*****************************************************************************
  tmp$dom_comuna <- as.character(tmp$dom_comuna)
  tot$dom_comuna <- as.character(tot$dom_comuna)
  
  #reformat
  names(tmp) <- c("year","age","order","n","death","earlyP",
                  "earlyQ","lateP","lateQ","dom_comuna")
  tmp$age <- as.numeric(as.character(tmp$age))
  tmp$order <- as.numeric(as.character(tmp$order))
  tmp$year <- as.numeric(as.character(tmp$year))
  tot$age <- as.numeric(as.character(tot$age))
  tot$denom <- as.numeric(as.character(tot$denom))
  tot$year <- as.numeric(as.character(tot$year))
  
  #calculate not-pregnant from total-pregnant
  #aggregate tmp to comuna*year*age
  tmp2 <- aggregate(tmp$n,by=list(tmp$year,tmp$age,tmp$dom_comuna),sum)
  names(tmp2) <- c("year","age","dom_comuna","n_preg")
  tot <- merge(tot,tmp2,by=c("year","age","dom_comuna"),all=T)
  tot$n_preg[is.na(tot$n_preg)] <- 0
  tot$n <- tot$denom - tot$n_preg
  tot$order <- 1
  tot$pregnant <- 0
  tot$death <- 0
  tot$earlyP <- 0
  tot$earlyQ <- 0
  tot$lateP <- 0
  tot$lateQ <- 0
  tmp$pregnant <- 1

  keepit <- c("dom_comuna","year","age","order","n","death","earlyP",
              "earlyQ","lateP","lateQ","pregnant")
  tot <- tot[,keepit]
  tmp <- tmp[,keepit]
  tot <- tot[tot$n > 0,]
  #merge birth and pop
  fin <- rbind(tot,tmp)
  fin <- fin[!(is.na(fin$death)),]
  
  #*****************************************************************************
  # (11) Morning After Pill Data
  #*****************************************************************************  
  f <- paste(ma.dir,"PillDistTime.csv",sep="")
  ma <- read.csv(f,sep=";")
  ma <- ma[,c("comuna_names","year","disponible","pilldistance")]
  names(ma) <- c("dom_comuna","year","pill","pilldistance")
  #lag of ~9 months
  ma$year <- ma$year + 1
  
  tmp <- merge(fin,ma,by=c("dom_comuna","year"),all=T)
  
  #recode levels
  tmp$pill <- as.character(tmp$pill)
  tmp$pill[tmp$pill == "2"] <- "0" #cases of rape
  tmp$pill[!(tmp$pill %in% c("0","1"))] <- NA
  tmp$pill[tmp$year <= 2009] <- 0
  tmp$pilldistance[tmp$year <= 2009] <- 0  
  
  #subset to years with all data
  fin <- tmp[complete.cases(tmp),]
  rm(tmp)
  
  #*****************************************************************************
  # (12) Mayor data
  #*****************************************************************************  
  if(usecom) {  
    f <- paste(pol.dir,"MayorCharacteristics_years.csv",sep="")
    pol <- read.csv(f,sep=";")
    pol <- pol[,c("comuna_dom", "year", "election", "mujer", "party", "votop")]
    names(pol) <- c("dom_comuna", "year", "election", "mujer", "party", "votop")
  
    tmp <- merge(fin,pol,by=c("dom_comuna","year"),all=T)
    fin <- tmp[complete.cases(tmp),]
    rm(tmp)
  }
  #*****************************************************************************
  # (13) Comuna Covariates
  #*****************************************************************************  
  if(usecom) {  
    f <- paste(com.dir,"Comunas_datos_20062012.csv",sep="")
    com <- read.csv(f,sep=";")
    com$urbind <- 0 
    com$urbind[com$urb>70]<-1
    
    tmp <- merge(com,fin,by=c("dom_comuna","year"),all=T)  
    fin <- tmp[complete.cases(tmp),]
    
    fin <- fin[,c("dom_comuna","educretiromedia","saludtotal","saludpersonal" ,
                  "saludcapacit", "eductotal","educmunic","mujeresindigente"  ,
                  "mujeresfuncionarias","pobreza","urb","urbind","year","dens",
                  "region","condom","usingcont","pill","mujer","party","votop",
                  "election","pilldistance","roadDist","travelTime","pregnant",
                  "age","order","n","death","earlyP","earlyQ","lateP","lateQ")]
    names(fin) <- c("dom_comuna","outofschool", "healthspend", "healthstaff"  , 
                    "healthtraining", "educationspend","educationmunicip"     ,
                    "femalepoverty","femaleworkers","poverty","urban","urbBin",
                    "year","density","region","condom","usingcont","pill"     ,
                    "mujer","party","votop","election","pilldistance"         ,
                    "roadDist","travelTime","pregnant","age","order","n"      ,
                    "death","earlyP","earlyQ","lateP","lateQ")    
    rm(tmp)  
  }
  #*****************************************************************************
  # (14) Export
  #*****************************************************************************
  write.csv(fin,filename,row.names=F)  
  
  rm(fin,ma,map,tmp2,tot)
  gc()
}
