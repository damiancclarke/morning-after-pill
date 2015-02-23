# BirthGenerate.R v1.20            KEL / DCC               yyyy-mm-dd:2013-12-12
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# Generates the file named in the argument "filename".  This file has one line 
# per comuna, age group and pregnancy status.  If usecom=TRUE the file created
# will also contain covariates for political controls, education and health con- 
# trols, gender controls, and frequency of contraceptive use.The municipality's 
# pill status (and distance to nearest comuna with pill if pill==0) is always
# included.
#
# This is a refactorization of managedata_KEL_20131124.R. This code has been 
# written by KEL, with updates by DCC to incorporate additional time varying 
# controls and pill distance data.
#
# prior to running this function the following directories must be defined:
#   > pop.dir (contains all population files from INE converted to .csv)
#   > brth.dir (directory containing dta file of all births in Chile)
#   > ma.dir (directory where pill data is stored)
#   > pol.dir (directory where mayor characteristics are kept)
#   > com.dir (directory where comuna characteristics and name file are kept)
#
# Last update v1.20: Refactorise

prep_s1_data <- function(age_range,usecom,filename) {

  #=============================================================================
  # (1) libraries
  #=============================================================================
  require(foreign)
  require(reshape)
  require(gdata)
  
  #=============================================================================
  # (2) comuna rename function
  #=============================================================================
  remove_accents <- function(string) {
    string <- trim(string)
    string <- sub("Á","A",string)
    string <- sub("É","E",string)
    string <- sub("Í","I",string)
    string <- sub("Ó","O",string)
    string <- sub("Ú","U",string)  
    return(string)  
  }
  
  #=============================================================================
  # (3) extract csv sheets from excel population data
  #=============================================================================
  if(F) {
    system('python ~universidades/Oxford/DPhil/Thesis/Teens/Source/xls2csv.py')
  }
  
  #=============================================================================
  # (4) Import and format popluation data
  #=============================================================================
  regions <- c(paste("0",1:9,sep=""),10:15)
  f.vec <- paste("SalComUsuarios-",regions, "Tok_1_ESAcal.csv",sep="")
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
    
    #discard male tables
    isblankcol <- apply(tmp,2,function(vec) {sum(is.na(vec)) == length(vec)})
    if(sum(isblankcol) > 1) print("Error -- too many blank columns")
    tmp1 <- tmp[,1:(which(isblankcol)-1)]
    tmp <- tmp[,(which(isblankcol)+1):dim(tmp)[2]]
    
    #Reformat female tables
    tmp <- rbind(c("INICIO",rep(NA,dim(tmp)[2]-1)),tmp)
    isnewtab <- which(tmp[,1] == "INICIO")
    comunas <- (as.character(tmp[isnewtab+1,1]))
    comunas <- lapply(as.list(comunas),function(c) 
                rep(c,length((age_range[1]):(age_range[2]))))
    comunas <- do.call("c",comunas)

    year <- tmp[4,]

    #generate tmp2, the tmp file consisting of ages of interst 
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
  
  #=============================================================================
  # (5) Rename comuna for merging with pill data
  #=============================================================================
  tot$dom_comuna <- apply(as.matrix(tot$dom_comuna),2,remove_accents)
  
  tot$dom_comuna[tot$dom_comuna == "T AMARILLA"] <- "TIERRA AMARILLA"
  tot$dom_comuna[tot$dom_comuna == "D DE ALMAGRO"] <- "DIEGO DE ALMAGRO"
  tot$dom_comuna[tot$dom_comuna == "A DEL CARMEN"] <- "ALTO DEL CARMEN"
  tot$dom_comuna[tot$dom_comuna == "ISLA DE PASCUA"] <- "ISLA  DE PASCUA"
  tot$dom_comuna[tot$dom_comuna == "Q. TILCOCO"] <- "QUINTA DE TILCOCO"
  tot$dom_comuna[tot$dom_comuna == "MARCHIGÜE"] <- "MARCHIHUE"
  tot$dom_comuna[tot$dom_comuna == "PEDRO AGUIRRE CERDA"] <- "PEDRO AGUIRRE CERDA"
  
  #=============================================================================
  # (6) Birth data
  #=============================================================================
  f <- paste(brth.dir,"Nacimientos_Chile_20002012.dta",sep="")
  tmp <- read.dta(f)

  # keep appropriate age groups and years of birth
  tmp <- tmp[tmp$edad_m <= age_range[2],]
  tmp <- tmp[tmp$edad_m >= age_range[1],]
  tmp <- tmp[tmp$ano_nac %in% birth_y_range,]
  tmp <- tmp[,c("comuna","ano_nac","edad_m","hij_total","urb_rural")]
  tmp$hij_total[tmp$hij_total > 4] <- 4
  tmp$hij_total[tmp$hij_total < 1] <- 1
  
  # aggregate to comuna
  tmp <- as.data.frame(table(tmp$comuna,tmp$ano_nac,tmp$edad_m,tmp$hij_total))
  tmp <- tmp[tmp$Freq > 0,]
  names(tmp) <- c("comuna","ano_nac","edad_m","hij_total","n")

  #=============================================================================
  # (6a) Comuna id in birth data
  #=============================================================================
  f <- paste(com.dir, "regionescomunas_short.csv",sep="")
  map <- read.csv(f,sep=";")
  map <- map[,1:7]
  map$dom_comuna <- toupper(as.character(map$COMUNA))
  map$dom_comuna <- apply(as.matrix(map$dom_comuna),2,remove_accents)
  map$comuna <- as.numeric(as.character(map$COMUNA.CODE..2010..))
  map$frommap <- 1
  tmp$fromdat <- 1
  tmp$dom_comuna <- NULL
  tmp <- merge(tmp,map[,c("comuna","frommap","dom_comuna")],by="comuna",
               all.x=T,all.y=F)
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
  
  #Pre 2000 comuna codes
  tmp00 <- tmp[as.numeric(tmp$comuna) < 1000,]
  map$comuna <- as.numeric(as.character(map$COMUNA.CODE..PRE.1999.))
  tmp00$dom_comuna <-  NULL
  tmp00$frommap <- NULL
  tmp00$fromdat <- 1
  tmp00 <- merge(tmp00,map[,c("comuna","dom_comuna","frommap")],by="comuna",
                 all.x=T,all.y=F)
  tmp00$matched <- tmp00$frommap & tmp00$fromdat
  tmp00 <- tmp00[tmp00$matched,]
  
  #clean up
  tmp00$matched <- tmp00$fromdat <- tmp00$frommap <- NULL
  tmp08$matched <- tmp08$fromdat <- tmp08$frommap <- NULL
  tmp$matched <- tmp$fromdat <- tmp$frommap <- NULL
  
  #combine
  deleteit <- c(tmp00$comuna,tmp08$comuna)
  tmp <- tmp[!(tmp$comuna %in% deleteit),]
  tmp <- rbind(tmp,tmp00,tmp08)
  tmp <- tmp[!(is.na(tmp$dom_comuna)),]
  tmp$comuna <- NULL
  
  rm(tmp00,tmp08)
  
  
  #=============================================================================
  # (7) Combine births and population
  #=============================================================================
  tmp$dom_comuna <- as.character(tmp$dom_comuna)
  tot$dom_comuna <- as.character(tot$dom_comuna)
  
  #reformat
  names(tmp) <- c("year","age","order","n","dom_comuna")
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
  tmp$pregnant <- 1
  keepit <- c("dom_comuna","year","age","order","n","pregnant")
  tot <- tot[,keepit]
  tmp <- tmp[,keepit]
  tot <- tot[tot$n > 0,]
  #merge birth and pop
  fin <- rbind(tot,tmp)
  
  #=============================================================================
  # (8) Pill data
  #=============================================================================
  f <- paste(ma.dir,"PillDist.csv",sep="")
  ma <- read.csv(f,sep=";")
  ma <- ma[,c("comuna_names","year","disponible","pilldistance")]
  names(ma) <- c("dom_comuna","year","pill","pilldistance")
  #lag of ~9 months
  ma$year <- ma$year + 1  
  ####ma$year <- ma$year - 3  
  tmp <- merge(fin,ma,by=c("dom_comuna","year"),all=T)
  
  #recode levels
  tmp$pill <- as.character(tmp$pill)
  tmp$pill[tmp$pill == "2"] <- "0" #cases of rape
  tmp$pill[!(tmp$pill %in% c("0","1"))] <- NA
  tmp$pill[tmp$year <= 2009] <- 0
  tmp$pilldistance[tmp$year <= 2009] <- 0
  ###tmp$pill[tmp$year <= 2006] <- 0
  ###tmp$pilldistance[tmp$year <= 2006] <- 0
  
  #subset to years with all data  
  fin <- tmp[complete.cases(tmp),]
  rm(tmp)

  #=============================================================================
  # (9) Mayor data
  #=============================================================================
  if(usecom) {  
    f <- paste(pol.dir,"MayorCharacteristics_years.csv",sep="")
    pol <- read.csv(f,sep=";")
    pol <- pol[,c("comuna_dom", "year", "election", "mujer", "party", "votop")]
    names(pol) <- c("dom_comuna", "year", "election", "mujer", "party", "votop")
  
    tmp <- merge(fin,pol,by=c("dom_comuna","year"),all=T)
    fin <- tmp[complete.cases(tmp),]
    rm(tmp)
  }
  #=============================================================================
  # (10) Comuna Covariates
  #=============================================================================
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
                  "election","pilldistance","pregnant","age","order","n")]
    names(fin) <- c("dom_comuna","outofschool", "healthspend", "healthstaff"  , 
                    "healthtraining", "educationspend","educationmunicip"     ,
                    "femalepoverty","femaleworkers","poverty","urban","urbBin",
                    "year","density","region","condom","usingcont","pill"     ,
                    "mujer","party","votop","election","pilldistance"         ,
                    "pregnant","age","order","n")
    rm(tmp)
  }  
  #=============================================================================
  # (11) Export
  #=============================================================================
  write.csv(fin,filename,row.names=F)
  
  rm(fin,ma,map,tmp2,tot)
  gc()  

}
