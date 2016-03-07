# BirthsEstimates.R v2.00          KEL / DCC               yyyy-mm-dd:2013-12-29
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# Import data from S1Data_granular_covars.csv to run various models.  Initial 
# data contains pregnancies by comuna and morning after pill availability by
# comuna.
# 
# Principal model is of the form:
#   preg_{ijt} = alpha + delta*PAE_{jt} + theta_j + i.eta**g(year_t) + u_{ijt}
#
# This code was first written by KEL, with updates by DCC to incorporate additi-
# onal time varying controls and export results to TeX. When running, the switc-
# hes in section (1) determine whether or not specific sections of the code will
# be run.  Note that the user-contributed stargazer library is very slow, so the
# "full" section is best avoided, or run only once.
#
#
# aboe refers to appended back of the envelope calculation which examines
# whether effect sizes look reasonable
#
# last edit
#          v2.00: Respond to referrees -- event study + OLS
#          v1.23: Clean directory structure.
# contact: damian.clarke@economics.ox.ac.uk

rm(list=ls())

#==============================================================================
#=== (1) Parameters
#==============================================================================
create <- FALSE
preg   <- FALSE
Npreg  <- TRUE
Lpreg  <- TRUE
prTab  <- TRUE
spill  <- FALSE
aboe   <- FALSE
ranges <- FALSE
events <- FALSE
ChMund <- FALSE
invPS  <- FALSE
robust <- FALSE

birth_y_range <- 2006:2012
pill_y_range  <- birth_y_range - 1
age_range     <- c(15,49)

#==============================================================================
#=== (2) Libraries, directories
#==============================================================================
ipak <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
       install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}
packages <- c("MatchIt","xtable","rms","plyr","glmmML","sandwich","lmtest",
	      "plotrix","foreign","reshape","gdata","Hmisc")
ipak(packages)

proj.dir <- "~/universidades/Oxford/DPhil/Thesis/Teens/"

brth.dir <- paste(proj.dir, "Data/Nacimientos/",sep="")
code.dir <- paste(proj.dir, "Source/Births/"   ,sep="")
com.dir  <- paste(proj.dir, "Data/Comunas/"    ,sep="")
ma.dir   <- paste(proj.dir, "Data/PAE/"        ,sep="")
pol.dir  <- paste(proj.dir, "Data/Alcaldes/"   ,sep="")
pop.dir  <- paste(proj.dir, "Data/Population/" ,sep="")
tab.dir  <- paste(proj.dir, "Tables/"          ,sep="")
graf.dir <- paste(proj.dir, "Figures/"         ,sep="")


Names <- c("dom_comuna","trend","trend2","pill","mujer","party","votes"      ,
           "outofschool","healthspend","healthstaff","healthtraining"        , 
           "educationspend","femalepoverty","year","urban","educationmunicip",
           "condom","usingcont","femaleworkers","region","popln")

#==============================================================================
#=== (3a) Source functions
#==============================================================================
if(create){
    f <- paste(code.dir,"BirthGenerate.R",sep="")
    source(f)

    filename <- paste(brth.dir, 'S1Data_granular_covars.csv' ,sep="")
    prep_s1_data(age_range,usecom="TRUE",filename)
}

#==============================================================================
#=== (3b) Load Data
#==============================================================================
f <- paste(brth.dir, "S1Data_granular_covars.csv", sep="")
orig <- read.csv(f)

#==============================================================================
#=== (4) Main Functions
#==============================================================================
stars <- function(p,B) {
    b <- ifelse(p < 0.01,
            paste(format(round(B,3),nsmall=3),"$^{***}$",sep=""),
            ifelse(p < 0.05,
               paste(format(round(B,3),nsmall=3),"$^{**}$",sep=""),
               ifelse(p < 0.1,
                   paste(format(round(B,3),nsmall=3),"$^{*}$",sep=""),
                   format(round(B,3),nsmall=3))))
    b  <- sub('-', '$-$', b)  
    return(b)
}

pillest <- function(outresults,d,n,regex,dim) {
    pillline <- grepl(regex,rownames(summary(outresults)$coefficients))
  
    if(dim==1|dim==3|dim==10) {
        beta <- summary(outresults)$coefficients[pillline,]["Estimate"]
        se   <- outresults$coefficients2[pillline,]["Std. Error"]
        if (dim==1) {p <- outresults$coefficients2[pillline,]["Pr(>|z|)"]}
        else {p <- outresults$coefficients2[pillline,]["Pr(>|t|)"]}
    }
    else {
        beta <- summary(outresults)$coefficients[pillline,][, "Estimate"]
        se   <- outresults$coefficients2[pillline,][, "Std. Error"]
        if (dim==11) {p <- outresults$coefficients2[pillline,][,"Pr(>|t|)"]}
        else {p    <- outresults$coefficients2[pillline,][, "Pr(>|z|)"]}
    }
  
    if (dim==1|dim==4) {
        null  <- glm(cbind(successes,failures) ~ 1,
                     family="binomial",data=d,weights=WT)
        Lfull <- as.numeric(logLik(outresults))
        Lnull <- as.numeric(logLik(null))
        R2    <- 1 - Lfull/Lnull
    }
    if(dim==10|dim==11) {
        R2 <- summary(outresults)$r.squared
    }
    beta  <- stars(p,beta)
    se    <- paste("[", format(round(se,3),nsmall=3),"]", sep="")
    R2    <- format(round(R2,3),nsmall=3)
    n     <- format(n,big.mark=",",scientific=F)
    
  
    return(list("b" = beta, "s" = se, "p" = p, "r" = R2, "n" = n))
}

robust.se <- function(model, cluster) {
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- model$rank
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
    rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
    rcse.se <- coeftest(model, rcse.cov)
    return(list(rcse.cov, rcse.se))
}

closegen <- function(d1,d2,dat) {
    dat2 <- dat
		# EUCLIDEAN DISTANCE
    dat2$newvar <- NA  
    dat2$newvar[dat2$pilldistance > d1 & dat2$pilldistance <= d2 &
                !(dat2$pilldistance)==0] <- 1
    dat2$newvar[is.na(dat2$newvar)]<-0

		# ROAD DISTANCE
    dat2$newvar2 <- NA  
    dat2$newvar2[dat2$roadDist/1000 > d1 & dat2$roadDist/1000 <= d2 &
                !(dat2$roadDist)==0] <- 1
    dat2$newvar2[is.na(dat2$newvar2)]<-0

		# TRAVEL TIME
    dat2$newvar3 <- NA  
    dat2$newvar3[dat2$travelTime/60 > d1 & dat2$travelTime/60 <= d2 &
                !(dat2$travelTime)==0] <- 1
    dat2$newvar3[is.na(dat2$newvar3)]<-0

    names(dat2)<-c(names(dat),paste('close',d2,sep=""),paste('road',d2,sep=""),
                   paste('time',d2,sep=""))
    return(dat2)
}

datcollapse <- function(age_sub,order_sub,ver,dat) {

    dat <- dat[dat$age %in% age_sub,]
    dat <- dat[(dat$order %in% order_sub) | !(dat$pregnant),]
    dat$popln <- ave(dat$n,dat$dom_comuna,dat$year,FUN=sum)

    dat <- closegen(0,15,dat)
    dat <- closegen(15,30,dat)
    dat <- closegen(30,45,dat)
    
    dat$failures <- (1-dat$pregnant)*dat$n
    dat$successes <- dat$pregnant*dat$n

    if(ver==2) {
        dat$n <- 1
        full  <- aggregate.data.frame(dat$n,
                                      by=list(dat$close15,dat$close30      ,
                                          dat$close45,dat$road15,dat$road30,
                                          dat$road45,dat$time15,dat$time30 ,
                                          dat$time45,dat$dom_comuna        ,
                                          dat$year-2005,(dat$year-2005)^2  ,
                                          dat$pill,dat$mujer,dat$party     ,
                                          dat$votop,dat$outofschool        ,
                                          dat$healthspend,dat$healthstaff  ,
                                          dat$healthtraining               ,
                                          dat$educationspend               ,
                                          dat$femalepoverty,dat$urbBin     ,
                                          dat$year,dat$educationmunicip    ,
                                          dat$condom,dat$usingcont         ,
                                          dat$femaleworkers,dat$region     ,
                                          dat$popln),
                                      function(vec) {sum(na.omit(vec))})
        names(full) <- c("close15","close30","close45","road15","road30","road45",
                         "time15","time30","time45",Names,"n")
        dat   <- dat[dat$pregnant==1,]
    }

    fmod <- aggregate.data.frame(dat[,c("failures","successes")],
                                 by=list(dat$close15,dat$close30,dat$close45,
                                     dat$road15,dat$road30,dat$road45       ,
                                     dat$time15,dat$time30,dat$time45       ,
                                     dat$dom_comuna,dat$year-2005           ,
                                     (dat$year-2005)^2,dat$pill,dat$mujer   ,
                                     dat$party,dat$votop,dat$outofschool    ,
                                     dat$healthspend,dat$healthstaff        ,
                                     dat$healthtraining,dat$educationspend  ,
                                     dat$femalepoverty,dat$urbBin,dat$year  ,
                                     dat$educationmunicip,dat$condom        ,
                                     dat$usingcont,dat$femaleworkers        ,
                                     dat$region,dat$popln),
                                 function(vec) {sum(na.omit(vec))})
    names(fmod) <- c("close15","close30","close45","road15","road30","road45",
                     "time15","time30","time45",Names,"failures","successes")
    
    if(ver==2) {
        fmod<-merge(fmod,full,all.y=TRUE)
        fmod$successes[is.na(fmod$successes)]<-0
    }

    fmod$healthstaff      <- fmod$healthstaff/100000
    fmod$healthspend      <- fmod$healthspend/100000
    fmod$healthtraining   <- fmod$healthtraining/100000
    fmod$educationspend   <- fmod$educationspend/100000
    fmod$educationmunicip <- fmod$educationmunicip/100000
    fmod$meanP            <- ave(fmod$pill, group=fmod$dom_comuna)

    return(fmod)
}

inversePS <- function(age_sub,order_sub,dat) {
    #CALCULATE PROPENSITY SCORE WEIGHT FOR EACH COMUNA BASED ON
    #PILL/NON-PILL STATUS.  USE FULL SET OF CONTROLS TO PREDICT PSCORE.

    psdat <- datcollapse(age_sub,order_sub,2,dat)
    psdat$conservative[psdat$party=="UDI"|psdat$party=="RN"] <- 1
    psdat$conservative[is.na(psdat$conservative)] <- 0

    psdat <- aggregate.data.frame(psdat[,c("pill","mujer","conservative"  ,
                                           "votes","outofschool","popln"  ,
                                           "healthspend","healthstaff"    ,
                                           "healthtraining","usingcont"   ,
                                           "femalepoverty","femaleworkers",
                                           "condom","educationspend"      ,
                                           "educationmunicip")],
                                  by=list(psdat$dom_comuna), FUN=mean)
    psdat$pillind[psdat$pill>0]<-1
    psdat$pillind[is.na(psdat$pillind)]<-0
        
    PSc <- glm(pillind ~ mujer + conservative + votes + outofschool +
               popln + healthspend + healthstaff + healthtraining   +
               usingcont + femalepoverty + condom + femaleworkers   +
               educationspend+ educationmunicip,
               family=binomial, data=psdat)
    psdat$predict<-predict(PSc,type="response")
    psdat$WT[psdat$pillind==1] <- 1/psdat$predict
    psdat$WT[psdat$pillind==0] <- 1/(1-psdat$predict)
    
    colnames(psdat)[1]<-"dom_comuna"
    wts <- psdat[,(names(psdat) %in% c("dom_comuna", "WT"))] 
    return(wts)
}

#==============================================================================
#=== (4a) Run various models to test effect of PAE on pregnancy 
#==============================================================================
runmod <- function(age_sub,order_sub,num,PSwt) {
  
    formod <- datcollapse(age_sub,order_sub,1,orig)
    if(PSwt) {
        wts <- inversePS(age_sub,order_sub,orig)
        formod <- join(formod,wts,by="dom_comuna",type="left",match="all")

        xPS <- glm(cbind(successes,failures) ~ factor(year) + factor(pill)    +
                   factor(dom_comuna) + factor(dom_comuna):trend,
                   family="binomial", data=formod, weights=WT)

        clusters <- mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
        xPS$coefficients2  <- robust.se(xPS,clusters)[[2]]
        n  <- sum(formod$successes) + sum(formod$failures)
        
        s0 <- pillest(xPS,   formod, n, "pill", 1)
        return(s0)
        
    } else {
        formod$WT <- 1
    }
    
    if(num==1) {
        xnT <-  glm(cbind(successes,failures) ~ factor(year) + factor(pill)   +
                    factor(dom_comuna) + factor(region):trend,
                    family="binomial", data=formod, weights=WT)
  
        xCM <-  glm(cbind(successes,failures) ~ factor(year) + factor(pill)   +
                    meanP + factor(region):trend,
                    family="binomial", data=formod, weights=WT)

        xtr <- glm(cbind(successes,failures) ~ factor(year) + factor(pill)    +
                   factor(dom_comuna) + factor(dom_comuna):trend,
                   family="binomial", data=formod, weights=WT)

        xpol <- glm(cbind(successes,failures) ~ factor(year) + factor(pill)   +
                    factor(dom_comuna) + factor(dom_comuna):trend + votes     +
                    factor(party) + factor(mujer)                             ,
                    family="binomial", data=formod, weights=WT)


        xsh  <- glm(cbind(successes,failures) ~ factor(year) + factor(pill)     +
                    factor(dom_comuna) + factor(dom_comuna):trend + votes       +
                    factor(party) + factor(mujer) + outofschool + educationspend+
                    educationmunicip + healthspend + healthtraining + healthstaff,
                    family="binomial", data=formod, weights=WT)

        xfem <- glm(cbind(successes,failures) ~ factor(year) + factor(pill)     +
                    factor(dom_comuna) + factor(dom_comuna):trend + votes       +
                    factor(party) + factor(mujer) + outofschool + educationspend+
                    educationmunicip + healthspend + healthtraining             +
                    healthstaff + femalepoverty + femaleworkers,
                    family="binomial", data=formod, weights=WT)

        clusters <- mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
        xnT$coefficients2  <- robust.se(xnT,clusters)[[2]]
        xCM$coefficients2  <- robust.se(xCM,clusters)[[2]]
        xtr$coefficients2  <- robust.se(xtr,clusters)[[2]]
        xpol$coefficients2 <- robust.se(xpol,clusters)[[2]]
        xsh$coefficients2  <- robust.se(xsh,clusters)[[2]]
        xfem$coefficients2 <- robust.se(xfem,clusters)[[2]]
  
    }
    xcont  <- glm(cbind(successes,failures) ~ factor(year) + factor(pill)       +
                  factor(dom_comuna) + factor(dom_comuna):trend + votes         +
                  factor(party) + factor(mujer) + outofschool + educationspend  +
                  educationmunicip + healthspend + healthtraining + healthstaff +
                  femalepoverty + femaleworkers + condom,
                  family="binomial", data=formod, weights=WT)
    clusters <-mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
    xcont$coefficients2 <- robust.se(xcont,clusters)[[2]]
  
    if(num==1) {
        n  <- sum(formod$successes) + sum(formod$failures)

        s0 <- pillest(xCM,   formod, n, "pill", 1)
        st <- pillest(xnT,   formod, n, "pill", 1)
        s1 <- pillest(xtr,   formod, n, "pill", 1)
        s2 <- pillest(xpol,  formod, n, "pill", 1)
        s3 <- pillest(xsh,   formod, n, "pill", 1)
        s4 <- pillest(xfem,  formod, n, "pill", 1)
        s5 <- pillest(xcont, formod, n, "pill", 1)

        betas <- paste(s1$b, "&", s2$b, "&", s4$b, "&", s5$b, sep="")
        ses   <- paste(s1$s, "&", s2$s, "&", s4$s, "&", s5$s, sep="")
        n     <- paste(s1$n, '&', s2$n, '&', s4$n, '&' ,s5$n, sep='')
        r     <- paste(s1$r, '&', s2$r, '&', s4$r, '&', s5$r, sep='')
        return(list("b" = betas,"se" = ses, "n" = n, "r" = r, "CM" = s0, "NT" = st))
    } else {
        return(xcont)
    }
}


NumMod <- function(age_sub,order_sub,rate,logm,wt) {
  
    formod <- datcollapse(age_sub, order_sub,2,orig)
    names(formod)[32] <- "births" 
    formod$weight <- 1
    
    if(rate) {
        formod$births <- (formod$births/formod$popln)*1000
    }
    if(logm) {
        formod$births<-log(formod$births+1)
    }
    if(wt) {
        formod$weight <- formod$popln
    }
    
    xba <- lm(births ~ factor(dom_comuna) + factor(year) + factor(pill),
              data=formod, weights=weight)
    xtr <- lm(births ~ factor(dom_comuna) + factor(year) + factor(pill)    +
              factor(dom_comuna):trend, data=formod, weights=weight)
    xct <- lm(births ~ factor(dom_comuna) + factor(dom_comuna):trend       +
              factor(year) + factor(pill) + factor(party) + factor(mujer)  +
              votes + outofschool + educationspend + educationmunicip      +
              healthspend + healthtraining + healthstaff + femalepoverty   + 
              femaleworkers + condom, data=formod, weights=weight)
    xsp <- lm(births ~ factor(dom_comuna) + factor(dom_comuna):trend       +
              factor(year) + factor(pill) + factor(party) + factor(mujer)  +
              votes + outofschool + educationspend + educationmunicip      +
              healthspend + healthtraining + healthstaff + femalepoverty   +
              condom + femaleworkers + factor(close15) + factor(close30)   +
              + factor(close45), data=formod, weights=weight)
    
    clusters <-mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
    xba$coefficients2 <- robust.se(xba,clusters)[[2]]
    xtr$coefficients2 <- robust.se(xtr,clusters)[[2]]
    xct$coefficients2 <- robust.se(xct,clusters)[[2]]
    xsp$coefficients2 <- robust.se(xsp,clusters)[[2]]
  
 
    formod$WT <- 1
    n  <- nrow(formod)
    mn <- format(round(wtd.mean(formod$births,weights=formod$popln),2),nsmall=2)
    s1 <- pillest(xba, formod, n, "pill", 10)
    s2 <- pillest(xtr, formod, n, "pill", 10)
    s3 <- pillest(xct, formod, n, "pill", 10)
    s4 <- pillest(xsp, formod, n, "pill", 10)
    
    betas <- paste(s1$b, "&", s2$b, "&", s3$b, "&", s4$b, sep="")
    ses   <- paste(s1$s, "&", s2$s, "&", s3$s, "&", s4$s, sep="")
    n     <- paste(s1$n, '&', s2$n, '&', s3$n, '&', s4$n, sep='')
    r     <- paste(s1$r, '&', s2$r, '&', s3$r, '&', s4$r, sep='')
    mean  <- paste(mn  , '&', mn  , '&', mn  , '&', mn  , sep='')
    return(list("b" = betas,"se" = ses, "n" = n, "r" = r, "m"=mean))  
}


#==============================================================================
#=== (4b) Various functions to examine effect of spillover 
#==============================================================================
spillovers <- function(age_sub,order_sub,time,road) {

    formod <- datcollapse(age_sub, order_sub,1,orig)
    if (time) {
        formod$close15 <- formod$time15
        formod$close30 <- formod$time30
        formod$close45 <- formod$time45
    }
    if (road) {
        formod$close15 <- formod$road15
        formod$close30 <- formod$road30
        formod$close45 <- formod$road45
    }
    
    xspill <- glm(cbind(successes,failures) ~ factor(dom_comuna)          + 
                  factor(dom_comuna):trend + factor(year) + factor(pill)  + 
                  factor(party) + factor(mujer) + votes + outofschool     + 
                  educationspend + educationmunicip + healthspend         + 
                  healthtraining + healthstaff + femalepoverty + condom   + 
                  femaleworkers + factor(close15) + factor(close30)       +
                  factor(close45), family="binomial",data=formod)
    clusters <-mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
    xspill$coefficients2 <- robust.se(xspill,clusters)[[2]]

    formod$WT <- 1
    n  <- sum(formod$successes) + sum(formod$failures)
    s1 <- pillest(xspill,formod,n,"pill|close",4)
  
    return(list("b" = s1$b,"s" = s1$s, "n" = s1$n, "r" = s1$r))
}

countpreg <- function(age_sub,order_sub,cond) {

    count <- orig
    count <- count[count$age %in% age_sub,]
    count <- count[count$order %in% order_sub,]
    count <- count[count$year %in% 2009:2011,]
    count <- count[count$pregnant==1,]
    
    if(cond==1) {
        count <- count[count$pill==1,]
    } else if (cond==2) {
        count <- count[count$pilldistance>0&count$pilldistance<15,]
    } else {
        count <- count[count$pilldistance>=15&count$pilldistance<30,]
    }
    
    count <- format(sum(count$n),big.mark=",",scientific=F)
    return(count)
}


rangeest <- function(age_sub,order_sub,measure,title,xlabl){

    formod <- datcollapse(age_sub, order_sub,1,orig)
    drops <- c("close15","close30","close45","road15","road30","road45",
 							 "time15","time30","time45")
    formod[,!(names(formod) %in% drops)]
    
    xrange <- glm(cbind(successes,failures) ~ factor(dom_comuna)          + 
                  factor(dom_comuna):trend + factor(year) + factor(pill)  +
                  factor(party) + factor(mujer) + votes + outofschool     + 
                  educationspend + educationmunicip + healthspend         + 
                  healthtraining + healthstaff + femalepoverty + condom   +
                  femaleworkers, family="binomial",data=formod)
    pillline  <- grepl("pill",rownames(summary(xrange)$coefficients))
    closeline <- grepl("closemarg",rownames(summary(xrange)$coefficients))  
    pillbeta  <- summary(xrange)$coefficients[pillline,]["Estimate"]
    pillse    <- summary(xrange)$coefficients[pillline,]["Std. Error"]
    closebeta <- summary(xrange)$coefficients[closeline,]["Estimate"]
    closese   <- summary(xrange)$coefficients[closeline,]["Std. Error"]
    distance  <- 0
    
    for(i in seq(2.5,45,2.5)) {
        cat(i,"\n")
        dat <- orig
        dat <- dat[dat$age %in% age_sub,]
        dat <- dat[(dat$order %in% order_sub) | !(dat$pregnant),]  
		    dat$popln <- ave(dat$n,dat$dom_comuna,dat$year,FUN=sum)
        n1  <- names(dat)
        dat <- closegen(0,i,dat)
        dat <- closegen(i,i+2.5,dat)
        
        dat$failures  <- (1-dat$pregnant)*dat$n
        dat$successes <- dat$pregnant*dat$n    
        names(dat) <- c(n1,"c1","c2","r1","r2","t1","t2","failures","successes")  
        
        if(measure=="dist") { 
            formod <- aggregate.data.frame(dat[,c("failures","successes")],
               by=list(dat$c1,dat$c2,dat$dom_comuna,dat$year-2005         ,
                   (dat$year-2005)^2,dat$pill,dat$mujer,dat$party,dat$votop,
                   dat$outofschool,dat$healthspend,dat$healthstaff         ,
                   dat$healthtraining,dat$educationspend,dat$femalepoverty ,
                   dat$urbBin,dat$year,dat$educationmunicip,dat$condom     ,
                   dat$usingcont,dat$femaleworkers,dat$region,dat$popln),
                                           function(vec) {sum(na.omit(vec))})
        }
        if(measure=="road") { 
            formod <- aggregate.data.frame(dat[,c("failures","successes")],
                 by=list(dat$r1,dat$r2,dat$dom_comuna,dat$year-2005         ,
                    (dat$year-2005)^2,dat$pill,dat$mujer,dat$party,dat$votop,
                    dat$outofschool,dat$healthspend,dat$healthstaff         ,
                    dat$healthtraining,dat$educationspend,dat$femalepoverty ,
                    dat$urbBin,dat$year,dat$educationmunicip,dat$condom     ,
                    dat$usingcont,dat$femaleworkers,dat$region,dat$popln),
                                           function(vec) {sum(na.omit(vec))})
        }
        if(measure=="time") {
            formod <- aggregate.data.frame(dat[,c("failures","successes")],
               by=list(dat$t1,dat$t2,dat$dom_comuna,dat$year-2005         ,
                   (dat$year-2005)^2,dat$pill,dat$mujer,dat$party,dat$votop,
                   dat$outofschool,dat$healthspend,dat$healthstaff         ,
                   dat$healthtraining,dat$educationspend,dat$femalepoverty ,
                   dat$urbBin,dat$year,dat$educationmunicip,dat$condom     ,
                   dat$usingcont,dat$femaleworkers,dat$region,dat$popln),
                                           function(vec) {sum(na.omit(vec))})
        }
        names(formod) <- c("close1","closemarg",Names,"failures","successes")
        
        xrange <- glm(cbind(successes,failures) ~ factor(dom_comuna)          + 
                      factor(dom_comuna):trend + factor(year)  + factor(pill) + 
                      factor(party) + factor(mujer) + votes + outofschool     + 
                      educationspend + educationmunicip + healthspend         + 
                      healthtraining + healthstaff + femalepoverty + condom   + 
                      femaleworkers + factor(close1) + factor(closemarg),
                      family="binomial",data=formod)
        pline     <- grepl("pill",rownames(summary(xrange)$coefficients))
        cline     <- grepl("closemarg",rownames(summary(xrange)$coefficients))  
        pillbeta  <- c(pillbeta,summary(xrange)$coefficients[pline,]["Estimate"])
        pillse    <- c(pillse,summary(xrange)$coefficients[pline,]["Std. Error"])
        closebeta <- c(closebeta,summary(xrange)$coefficients[cline,]["Estimate"])
        closese   <- c(closese,summary(xrange)$coefficients[cline,]["Std. Error"])
        distance  <- c(distance, i)  
    }

    postscript(paste(graf.dir,title,sep=""),horizontal = FALSE, 
               onefile = FALSE, paper = "special",height=7, width=9)
    plot(distance,pillbeta, type="b",ylim=c(-0.10,-0.02),lwd=2,pch=20,
         col="darkgreen", ylab="Estimate of Effect on Treated Cluster",
         xlab=xlabl)
    points(distance,pillbeta-1.96*pillse,type="l",lty=3,pch=20)
    points(distance,pillbeta+1.96*pillse,type="l",lty=3,pch=20)
    legend("topright",legend=c("Point Estimate","95% CI"),
           text.col=c("darkgreen","black"),pch=c(20,NA),lty=c(1,3),
           col=c("darkgreen","black"))
    dev.off()
    
    
    return(data.frame(distance,pillbeta,pillse,closebeta,closese))
}

#==============================================================================
#=== (4c) Event study
#==============================================================================
event <- function(age_sub,order_sub) {
  
    formod <- datcollapse(age_sub, order_sub,1,orig)
    formod <- formod[with(formod,order(dom_comuna,trend)), ]

    formod$pillbinary <- ave(formod$pill,formod$dom_comuna,FUN=sum)
    formod$treatCom[formod$pillbinary>0]  <- 1
    formod$treatCom[formod$pillbinary==0] <- 0
    formod$pilltotal <- ave(formod$pill,formod$dom_comuna,FUN=cumsum)

    formod$nopill <- 0
    formod$nopill[formod$pilltotal==0] <- 1
    formod           <- formod[with(formod,order(dom_comuna,trend,decreasing=T)), ]
    formod$add       <- ave(formod$nopill,formod$dom_comuna,FUN=cumsum)
    
    formod$pilln5[formod$add==5 & formod$treatCom==1]   <- 1
    formod$pilln5[is.na(formod$pilln5)]                 <- 0
    formod$pilln4[formod$add==4 & formod$treatCom==1]   <- 1
    formod$pilln4[is.na(formod$pilln4)]                 <- 0
    formod$pilln3[formod$add==3 & formod$treatCom==1]   <- 1
    formod$pilln3[is.na(formod$pilln3)]                 <- 0
    formod$pilln2[formod$add==2 & formod$treatCom==1]   <- 1
    formod$pilln2[is.na(formod$pilln2)]                 <- 0
    formod$pilln1[formod$add==1 & formod$treatCom==1]   <- 1
    formod$pilln1[is.na(formod$pilln1)]                 <- 0
    formod$pillp0[formod$pill==1 & formod$pilltotal==1] <- 1
    formod$pillp0[is.na(formod$pillp0)]                 <- 0
    formod$pillp1[formod$pill==1 & formod$pilltotal==2] <- 1
    formod$pillp1[is.na(formod$pillp1)]                 <- 0
    formod$pillp2[formod$pill==1 & formod$pilltotal==3] <- 1
    formod$pillp2[is.na(formod$pillp2)]                 <- 0




    eventS  <- glm(cbind(successes,failures) ~ factor(year)                      +
                   factor(dom_comuna) + factor(dom_comuna):trend + votes         +
                   factor(party) + factor(mujer) + outofschool + educationspend  +
                   educationmunicip + healthspend + healthtraining + healthstaff +
                   femalepoverty + femaleworkers + condom + factor(pilln5)       +
                   factor(pilln4) + factor(pilln2) + factor(pilln1)              +
                   factor(pillp0) + factor(pillp1) + factor(pillp2),
                   family="binomial", data=formod)
    clusters <-mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
    eventS$coefficients2 <- robust.se(eventS,clusters)[[2]]


    pillline <- grepl("pill",rownames(summary(eventS)$coefficients))
    beta <- summary(eventS)$coefficients[pillline,][, "Estimate"]
    se   <- eventS$coefficients2[pillline,][, "Std. Error"]

    return(list("b" = beta, "s" = se, "eventyr" = c(-5,-4,-2,-1,0,1,2)))
}


#==============================================================================
#=== (5) Estimate
#==============================================================================
if(preg|ChMund){
    print("Rate of Pregnancy Logits")
    a1519 <- runmod(age_sub = 15:19, order_sub = 1:100,1,PSwt=FALSE)
    a2034 <- runmod(age_sub = 20:34, order_sub = 1:100,1,PSwt=FALSE)
    a3549 <- runmod(age_sub = 35:49, order_sub = 1:100,1,PSwt=FALSE)
    aAll  <- runmod(age_sub = 15:49, order_sub = 1:100,1,PSwt=FALSE)
    if(preg) {
        b1519 <- runmod(age_sub = 15:19, order_sub = 1,1,PSwt=FALSE)
        b2034 <- runmod(age_sub = 20:34, order_sub = 1,1,PSwt=FALSE)
        b3549 <- runmod(age_sub = 35:49, order_sub = 1,1,PSwt=FALSE)
    }
}

if(Npreg){
    print("Pregnancy Rate - weighted")
    w1519 <- NumMod(age_sub = 15:19, order_sub = 1:100,rate=T,logm=F,wt=T)
    w2034 <- NumMod(age_sub = 20:34, order_sub = 1:100,rate=T,logm=F,wt=T)
    w3549 <- NumMod(age_sub = 35:49, order_sub = 1:100,rate=T,logm=F,wt=T)
    wAll  <- NumMod(age_sub = 15:49, order_sub = 1:100,rate=T,logm=F,wt=T)

    print("Number of Pregnancy - weighted")
    x1519 <- NumMod(age_sub = 15:19, order_sub = 1:100,rate=F,logm=F,wt=T)
    x2034 <- NumMod(age_sub = 20:34, order_sub = 1:100,rate=F,logm=F,wt=T)
    x3549 <- NumMod(age_sub = 35:49, order_sub = 1:100,rate=F,logm=F,wt=T)
    xAll  <- NumMod(age_sub = 15:49, order_sub = 1:100,rate=F,logm=F,wt=T)

    print("Pregnancy Rate - unweighted")
    r1519 <- NumMod(age_sub = 15:19, order_sub = 1:100,rate=T,logm=F,wt=F)
    r2034 <- NumMod(age_sub = 20:34, order_sub = 1:100,rate=T,logm=F,wt=F)
    r3549 <- NumMod(age_sub = 35:49, order_sub = 1:100,rate=T,logm=F,wt=F)
    rAll  <- NumMod(age_sub = 15:49, order_sub = 1:100,rate=T,logm=F,wt=F)

    print("Number of Pregnancy - unweighted")
    N1519 <- NumMod(age_sub = 15:19, order_sub = 1:100,rate=F,logm=F,wt=F)
    N2034 <- NumMod(age_sub = 20:34, order_sub = 1:100,rate=F,logm=F,wt=F)
    N3549 <- NumMod(age_sub = 35:49, order_sub = 1:100,rate=F,logm=F,wt=F)
    NAll  <- NumMod(age_sub = 15:49, order_sub = 1:100,rate=F,logm=F,wt=F)
}

if(Lpreg){
    print("ln(Number of Pregnancy) OLS - unweighted")
    l1519 <- NumMod(age_sub = 15:19, order_sub = 1:100,rate=F,logm=T,wt=F)
    l2034 <- NumMod(age_sub = 20:34, order_sub = 1:100,rate=F,logm=T,wt=F)
    l3549 <- NumMod(age_sub = 35:49, order_sub = 1:100,rate=F,logm=T,wt=F)
    lAll  <- NumMod(age_sub = 15:49, order_sub = 1:100,rate=F,logm=T,wt=F)

    print("ln(Number of Pregnancy) OLS - weighted")
    m1519 <- NumMod(age_sub = 15:19, order_sub = 1:100,rate=F,logm=T,wt=T)
    m2034 <- NumMod(age_sub = 20:34, order_sub = 1:100,rate=F,logm=T,wt=T)
    m3549 <- NumMod(age_sub = 35:49, order_sub = 1:100,rate=F,logm=T,wt=T)
    mAll  <- NumMod(age_sub = 15:49, order_sub = 1:100,rate=F,logm=T,wt=T)
}

if(spill){
    print("Spillover Models")
    c1519 <- spillovers(age_sub = 15:19, order_sub = 1:100,time=FALSE,road=FALSE)
    c2034 <- spillovers(age_sub = 20:34, order_sub = 1:100,time=FALSE,road=FALSE)
    c3549 <- spillovers(age_sub = 35:49, order_sub = 1:100,time=FALSE,road=FALSE)
    cAll  <- spillovers(age_sub = 15:49, order_sub = 1:100,time=FALSE,road=FALSE)
    print("Spillover Models (Road)")
    d1519 <- spillovers(age_sub = 15:19, order_sub = 1:100,time=FALSE,road=TRUE)
    d2034 <- spillovers(age_sub = 20:34, order_sub = 1:100,time=FALSE,road=TRUE)
    d3549 <- spillovers(age_sub = 35:49, order_sub = 1:100,time=FALSE,road=TRUE)
    dAll  <- spillovers(age_sub = 15:49, order_sub = 1:100,time=FALSE,road=TRUE)
    print("Spillover Models (Time)")
    e1519 <- spillovers(age_sub = 15:19, order_sub = 1:100,time=TRUE,road=FALSE)
    e2034 <- spillovers(age_sub = 20:34, order_sub = 1:100,time=TRUE,road=FALSE)
    e3549 <- spillovers(age_sub = 35:49, order_sub = 1:100,time=TRUE,road=FALSE)
    eAll  <- spillovers(age_sub = 15:49, order_sub = 1:100,time=TRUE,road=FALSE)

}  

if(invPS){
    print("Rate of Pregnancy Models with Inverse Propensity Weights")
    ps1519 <- runmod(age_sub = 15:19, order_sub = 1:100,1,PSwt=TRUE)
    ps2034 <- runmod(age_sub = 20:34, order_sub = 1:100,1,PSwt=TRUE)
    ps3549 <- runmod(age_sub = 35:49, order_sub = 1:100,1,PSwt=TRUE)
    psAll  <- runmod(age_sub = 15:49, order_sub = 1:100,1,PSwt=TRUE)
}

if(aboe) {
    NPp18   <- countpreg(age_sub = 15:18, order_sub = 1:100,1)
    NPp19   <- countpreg(age_sub = 20:49, order_sub = 1:100,1)
    NPc1518 <- countpreg(age_sub = 15:18, order_sub = 1:100,2)
    NPc1519 <- countpreg(age_sub = 20:49, order_sub = 1:100,2)
    NPc3018 <- countpreg(age_sub = 15:18, order_sub = 1:100,3)
    NPc3019 <- countpreg(age_sub = 19:49, order_sub = 1:100,3) 
    
    boe18 <- spillovers(age_sub = 15:18, order_sub = 1:100)
    boe19 <- spillovers(age_sub = 19:49, order_sub = 1:100)
}

if(ranges) {
    raDist1519  <- rangeest(age_sub = 15:19, 1:100,"dist","Dist1519.eps",
                            "Distance From Treatment Cluster (km)")
    raRoad1519  <- rangeest(age_sub = 15:19, 1:100,"road","Dist1519road.eps",
                            "Distance From Treatment Cluster (shortest distance by road)")
    raTime1519  <- rangeest(age_sub = 15:19, 1:100,"time","Dist1519time.eps",
                            "Distance From Treatment Cluster (minutes in travel time)")
    
    raDist2034  <- rangeest(age_sub = 20:34, 1:100,"dist","Dist2034.eps",
                            "Distance From Treatment Cluster (km)")
    raRoad2034  <- rangeest(age_sub = 20:34, 1:100,"road","Dist2034road.eps",
                            "Distance From Treatment Cluster (shortest distance by road)")
    raTime2034  <- rangeest(age_sub = 20:34, 1:100,"time","Dist2034time.eps",
                            "Distance From Treatment Cluster (minutes in travel time)")
}


if(events){
    note <- "Points and confidence intervals represent etimates for a
    full event study. Each point represents treatment interacted with n years
    prior/posterior to the reform. Lag 3 is omitted as the arbitrary base category."
    note <- sub("\n","",note)
    note <- gsub("  "," ",note)
                  
    
    e1519 <- event(age_sub=15:19, order_sub=1)
    
    postscript(paste(graf.dir,"Event1519.eps",sep=""),
             horizontal = FALSE, onefile = FALSE, paper = "special",
             height=7, width=9)

    plotCI(e1519$eventyr,e1519$b,ui=e1519$b+1.96*e1519$s,li=e1519$b-1.96*e1519$s,
           ylim=c(-0.2,0.10), ylab="Estimate", xlab="Event Year",
           cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)
    
    points(e1519$eventyr,e1519$b,type="l",lwd=2,pch=20)
    abline(h =  0  , lwd=1, col="gray60", lty = 2)
    abline(v = -0.1, lwd=2, col="red")
    dev.off()
    

    e2034 <- event(age_sub=20:34, order_sub=1)
    postscript(paste(graf.dir,"Event2034.eps",sep=""),
             horizontal = FALSE, onefile = FALSE, paper = "special",
             height=7, width=9)

    plotCI(e2034$eventyr,e2034$b,ui=e2034$b+1.96*e2034$s,li=e2034$b-1.96*e2034$s,
           ylim=c(-0.2,0.10), ylab="Estimate", xlab="Event Year",
           cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)
    points(e2034$eventyr,e2034$b,type="l",lwd=2,pch=20)
    abline(h =  0  , lwd=1, col="gray60", lty = 2)
    abline(v = -0.1, lwd=2, col="red")
    dev.off()


    e3549 <- event(age_sub=35:49, order_sub=1)
    postscript(paste(graf.dir,"Event3549.eps",sep=""),
             horizontal = FALSE, onefile = FALSE, paper = "special",
             height=7, width=9)
    plotCI(e3549$eventyr,e3549$b,ui=e3549$b+1.96*e3549$s,li=e3549$b-1.96*e3549$s,
           ylab="Estimate", xlab="Event Year",
           cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)
    points(e3549$eventyr,e3549$b,type="l",lwd=2,pch=20)
    abline(h =  0  , lwd=1, col="gray60", lty = 2)
    abline(v = -0.1, lwd=2, col="red")

    dev.off()
}


#==============================================================================
#=== (6) Export results
#===            1st block: Pregnancy Results
#===            2nd block: Spillover Results
#==============================================================================
xvar <- 'Morning After Pill &'
xv2  <- 'Close $<15$ km &'
xv3  <- 'Close 15-30 km &'
xv4  <- 'Close 30-45 km &'

obs  <- 'Observations&'
R2   <- 'McFadden\'s $R^2$&'
sig  <- '$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01'

prTab2 = FALSE
if(prTab2){
to <-file(paste(tab.dir,"Births.tex", sep=""))
writeLines(c('\\begin{table}[!htbp] \\centering',
           '\\caption{The Effect of the Morning After Pill on Pregnancy}',
           '\\label{TEENtab:PillPreg}',
           '\\begin{tabular}{@{\\extracolsep{5pt}}lcccc}',
           '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
           '&Pr(Birth)&Pr(Birth)&Pr(Birth)&Pr(Birth)\\\\ ',
           '&(1)&(2)&(3)&(4)\\\\ \\hline',
           '\\multicolumn{5}{l}{\\textsc{\\noindent All Women}} \\\\',
           ' & & & & \\\\',
           paste(xvar,aAll$b,'\\\\', sep=""), 
           paste(' &', aAll$se, '\\\\',sep=""), 
           ' & & & & \\\\',
           paste(obs, aAll$n,'\\\\', sep=""), 
           paste(R2, aAll$r,'\\\\', sep=""), 
           ' & & & & \\\\',
           '\\multicolumn{5}{l}{\\textsc{\\noindent 15-19 year olds}} \\\\',
           ' & & & & \\\\',
           paste(xvar,a1519$b,'\\\\', sep=""), 
           paste(' &', a1519$se, '\\\\', sep=""), 
           ' & & & & \\\\',
           paste(obs, a1519$n, '\\\\', sep=""), 
           paste(R2, a1519$r,'\\\\', sep=""), 
           ' & & & & \\\\',
           '\\multicolumn{5}{l}{\\textsc{\\noindent 20-34 year olds}} \\\\',
           ' & & & & \\\\', 
           paste(xvar,a2034$b,'\\\\', sep=""), 
           paste(' &', a2034$se, '\\\\', sep=""), 
           ' & & & & \\\\',
           paste(obs, a2034$n,'\\\\', sep=""), 
           paste(R2, a2034$r,'\\\\', sep=""), 
           ' & & & & \\\\',
           '\\multicolumn{5}{l}{\\textsc{\\noindent 35-49 year olds}} \\\\',
           ' & & & & \\\\', 
           paste(xvar,a3549$b,'\\\\', sep=""), 
           paste(' &', a3549$se,'\\\\', sep=""), 
          ' & & & & \\\\',
           paste('Observations&',a3549$n,'\\\\', sep=""), 
           paste(R2, a3549$r,'\\\\', sep=""), 
           '\\hline \\\\[-1.8ex] ', 
           '{\\small Trends \\& FEs} & Y & Y & Y & Y \\\\',
           '{\\small Political Controls} & & Y & Y & Y \\\\', 
           '{\\small Health, Educ, Gender Controls} & & & Y & Y\\\\',
           '{\\small Condom Availability} & & & & Y\\\\', 
           '\\hline \\hline \\\\[-1.8ex]',
           '\\multicolumn{5}{p{13.2cm}}{\\begin{footnotesize}\\textsc{Notes:}',
           ' Birth is a binary variable taking 1 for women who give birth    ',
           'and 0 if the woman does not give birth. All models are estimated ',
           ' by weighted logit, where weights are the number of women who    ',
           'give birth or do not give birth (respectively). Observations     ',
           'refer to the weighted number of women. Municipality and          ',
           ' year fixed effects are included, and standard errors are        ',
           'clustered at the level of the municipality. All coefficients are ',
           'reported as log odds and in each case Pill is a binary variable  ',
           'referring to the availability of the morning after pill in the   ',
           'woman\'s municipality and (lagged) year. Political controls are  ',
           'party dummies for the mayor in power, the mayor\'s gender, and   ',
           'the vote margin of the mayor.  Health and education controls     ',
           'are the percent of girls out of highschool, education spending   ',
           'spending by both the municipality and the Ministry of Education  ',
           'and total health spending and health spending on staff and       ',
           'training.  Gender controls are the percent of female heads of    ',
           ' households living below the poverty line, and the percent of    ',
           'female workers in professional positions in the Municipality,    ',
           'and condom availability is measured from survey data at the      ',
           'level of the region.',
           paste(sig, '\\end{footnotesize}}', sep=""),
           '\\normalsize\\end{tabular}\\end{table}'),to)
close(to)

to <-file(paste(tab.dir,"BirthsBoth.tex", sep=""))
writeLines(c('\\begin{landscape}','\\begin{table}[!htbp] \\centering',
           '\\caption{The Effect of the Morning After Pill on Pregnancy}',
           '\\label{TEENtab:PillPregBoth}',
           '\\begin{tabular}{@{\\extracolsep{5pt}}lccccp{1mm}cccc}',
           '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
           '&\\multicolumn{4}{c}{All Births}&&\\multicolumn{4}{c}{First Births}',
           '\\\\ \\cmidrule(r){2-5} \\cmidrule(r){7-10}',
           '&(1)&(2)&(3)&(4)&&(5)&(6)&(7)&(8)\\\\ \\hline',
           '\\multicolumn{10}{l}{\\textsc{\\noindent 15-19 year olds}} \\\\',
           ' & & & & & & & & & \\\\',
           paste(xvar,a1519$b,'&&',b1519$b,'\\\\', sep=""), 
           paste(' &', a1519$se,'&&', b1519$se, '\\\\', sep=""), 
           ' & & & & & & & & & \\\\',
           paste(obs, a1519$n,'&&',b1519$n,'\\\\', sep=""), 
           paste(R2, a1519$r,'&&',b1519$r,'\\\\', sep=""), 
           ' & & & & & & & & & \\\\',
           '\\multicolumn{10}{l}{\\textsc{\\noindent 20-34 year olds}} \\\\',
           ' & & & & & & & & & \\\\', 
           paste(xvar,a2034$b,'&&',b2034$b,'\\\\', sep=""), 
           paste(' &', a2034$se,'&&', b2034$se, '\\\\', sep=""), 
           ' & & & & & & & & & \\\\',
           paste(obs, a2034$n,'&&',b2034$n,'\\\\', sep=""), 
           paste(R2, a2034$r,'&&',b2034$r,'\\\\', sep=""), 
           ' & & & & & & & & & \\\\',
           '\\multicolumn{10}{l}{\\textsc{\\noindent 35-49 year olds}} \\\\',
           ' & & & & & & & & & \\\\', 
           paste(xvar,a3549$b,'&&',b3549$b,'\\\\', sep=""), 
           paste(' &', a3549$se,'&&', b3549$se, '\\\\', sep=""), 
          ' & & & & & & & & & \\\\',
           paste('Observations&',a3549$n,'&&',b3549$n,'\\\\', sep=""), 
           paste(R2, a3549$r,'&&',b3549$r,'\\\\', sep=""), 
           '\\hline \\\\[-1.8ex] ', 
           '{\\small Trends \\& FEs} & Y & Y & Y & Y && Y & Y & Y & Y \\\\',
           '{\\small Political Controls} & & Y & Y & Y && & Y & Y & Y \\\\', 
           '{\\small Health, Educ, Gender Controls} & & & Y & Y && & & Y & Y \\\\',
           '{\\small Condom Availability} & & & & Y && & & & Y \\\\', 
           '\\hline \\hline \\\\[-1.8ex]',
           '\\multicolumn{10}{p{22cm}}{\\begin{footnotesize}\\textsc{Notes:}',
           'Refer to notes in table \\ref{TEENtab:PillPreg} of the paper.   ',
           'For all births, logits are estimated where women giving birth   ',
           'in the year are assigned the value 1, and women not giving birth',
           ' are assigned zero.  For first births, only women who have their',
           ' first birth are assigned one.  All women who do not give birth ',
           ' are assigned the value of zero, so results are expressed as    ',
           ' first births/all women not giving birth.                       ',
           paste(sig, '\\end{footnotesize}}', sep=""),
           '\\normalsize\\end{tabular}\\end{table}\\end{landscape}'),to)
close(to)
}

if(spill){
  to <-file(paste(tab.dir,"Spillovers_A.tex", sep=""))
  writeLines(c('\\begin{table}[!htbp] \\centering',
             '\\caption{The Morning After Pill and Treatment Spillovers}',
             '\\label{TEENtab:Spillover} \\begin{tabular}',
             '{@{\\extracolsep{5pt}}lcccc}\\\\[-1.8ex]\\hline\\hline\\\\',
             '[-1.8ex] & All & 15-19 & 20-34 & 35-49 \\\\',
             '& Women & Year olds & Year olds & Year olds \\\\ \\midrule',
             '\\multicolumn{5}{l}{\\textsc{\\noindent Panel A: Births}} \\\\',
             '& & & & \\\\',
             paste(xvar,cAll$b[1],'&',c1519$b[1],'&',c2034$b[1],'&',c3549$b[1],'\\\\',sep=""),
             paste('&',cAll$s[1],'&',c1519$s[1],'&',c2034$s[1],'&',c3549$s[1],'\\\\',sep=""),            
             paste(xv2,cAll$b[2],'&',c1519$b[2],'&',c2034$b[2],'&',c3549$b[2],'\\\\',sep=""),
             paste('&',cAll$s[2],'&',c1519$s[2],'&',c2034$s[2],'&',c3549$s[2],'\\\\',sep=""),
             paste(xv3,cAll$b[3],'&',c1519$b[3],'&',c2034$b[3],'&'           ,'\\\\',sep=""),
             paste('&',cAll$s[3],'&',c1519$s[3],'&',c2034$s[3],'&'           ,'\\\\',sep=""), 
             paste(xv4,cAll$b[4],'&',c1519$b[4],'&'           ,'&'           ,'\\\\',sep=""),
             paste('&',cAll$s[4],'&',c1519$s[4],'&'           ,'&'           ,'\\\\',sep=""), 
             '& & & & \\\\',
             paste(obs,cAll$n,'&',c1519$n,'&',c2034$n,'&',c3549$n,'\\\\',sep=""),
             paste(R2,cAll$r,'&',c1519$r,'&',c2034$r,'&',c3549$r,'\\\\ \\midrule',sep="")),
             to)
  close(to)

  to <-file(paste(tab.dir,"SpilloversROAD_A.tex", sep=""))
  writeLines(c('\\begin{table}[!htbp] \\centering',
             '\\caption{The Morning After Pill and Treatment Spillovers (Roads)}',
             '\\label{TEENtab:SpilloverRoads} \\begin{tabular}',
             '{@{\\extracolsep{5pt}}lcccc}\\\\[-1.8ex]\\hline\\hline\\\\',
             '[-1.8ex] & All & 15-19 & 20-34 & 35-49 \\\\',
             '& Women & Year olds & Year olds & Year olds \\\\ \\midrule',
             '\\multicolumn{5}{l}{\\textsc{\\noindent Panel A: Births}} \\\\',
             '& & & & \\\\',
               paste(xvar,dAll$b[1],'&',d1519$b[1],'&',
                     d2034$b[1],'&',d3549$b[1],'\\\\',sep=""),
               paste('&' ,dAll$s[1],'&',d1519$s[1],'&',
                     d2034$s[1],'&',d3549$s[1],'\\\\',sep=""),            
               paste(xv2 ,dAll$b[2],'&',d1519$b[2],'&',
                     d2034$b[2],'&',d3549$b[2],'\\\\',sep=""),
               paste('&' ,dAll$s[2],'&',d1519$s[2],'&',
                     d2034$s[2],'&',d3549$s[2],'\\\\',sep=""),
               paste(xv3, dAll$b[3],'&',d1519$b[3],'&',
                     d2034$b[3],'&',d3549$b[3],'\\\\',sep=""),
               paste('&', dAll$s[3],'&',d1519$s[3],'&',
                     d2034$s[3],'&',d3549$s[3],'\\\\',sep=""), 
             '& & & & \\\\',
               paste(obs,dAll$n,'&',d1519$n,'&',d2034$n,
                     '&',d3549$n,'\\\\',sep=""),
               paste(R2,dAll$r,'&',d1519$r,'&',d2034$r,
                     '&',d3549$r,'\\\\ \\midrule',sep="")),
             to)
  close(to)
  
  xv2  <- 'Close $<15$ mins &'
  xv3  <- 'Close 15-30 mins &'
  xv4  <- 'Close 30-45 mins &'

  to <-file(paste(tab.dir,"SpilloversTIME_A.tex", sep=""))
  writeLines(c('\\begin{table}[!htbp] \\centering',
             '\\caption{The Morning After Pill and Treatment Spillovers (Time)}',
             '\\label{TEENtab:SpilloverTime} \\begin{tabular}',
             '{@{\\extracolsep{5pt}}lcccc}\\\\[-1.8ex]\\hline\\hline\\\\',
             '[-1.8ex] & All & 15-19 & 20-34 & 35-49 \\\\',
             '& Women & Year olds & Year olds & Year olds \\\\ \\midrule',
             '\\multicolumn{5}{l}{\\textsc{\\noindent Panel A: Births}} \\\\',
             '& & & & \\\\',
               paste(xvar,eAll$b[1],'&',e1519$b[1],'&',
                     e2034$b[1],'&',e3549$b[1],'\\\\',sep=""),
               paste('&', eAll$s[1],'&',e1519$s[1],'&',
                     e2034$s[1],'&',e3549$s[1],'\\\\',sep=""),            
               paste(xv2, eAll$b[2],'&',e1519$b[2],'&',
                     e2034$b[2],'&',e3549$b[2],'\\\\',sep=""),
               paste('&', eAll$s[2],'&',e1519$s[2],'&',
                     e2034$s[2],'&',e3549$s[2],'\\\\',sep=""),
               paste(xv3, eAll$b[3],'&',e1519$b[3],'&',
                     e2034$b[3],'&'           ,'\\\\',sep=""),
               paste('&', eAll$s[3],'&',e1519$s[3],'&',
                     e2034$s[3],'&'           ,'\\\\',sep=""), 
               paste(xv4, eAll$b[4],'&',e1519$b[4],'&'
                    ,'&'           ,'\\\\',sep=""),
               paste('&', eAll$s[4],'&',e1519$s[4],'&'
                    ,'&'           ,'\\\\',sep=""), 
             '& & & & \\\\',
               paste(obs, eAll$n,'&',e1519$n,'&',
                     e2034$n,'&',e3549$n,'\\\\',sep=""),
               paste(R2,  eAll$r,'&',e1519$r,'&',
                     e2034$r,'&',e3549$r,'\\\\ \\midrule',sep="")),
             to)
  close(to)
}

if(aboe) {
to <-file(paste(tab.dir,"Consistency.tex", sep=""))
writeLines(c('\\begin{table}[!htbp] \\centering',
           '\\caption{Back of the Envelope Calculation of Effect Sizes}',
           '\\label{TEENtab:BOE}',
           '\\begin{tabular}{@{\\extracolsep{5pt}}lcc}',
           '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
           '& 18 \\& Under & 19 \\& Over\\\\ ',
           '&(1)&(2) \\\\ \\hline',
           ' & &  \\\\',
            paste(xvar,boe18$b[1],'&',boe19$b[1],'\\\\',sep=""),
            paste('&', boe18$s[1],'&',boe19$s[1],'\\\\',sep=""),            
            paste(xv2, boe18$b[2],'&',boe19$b[2],'\\\\',sep=""),
            paste('&', boe18$s[2],'&',boe19$s[2],'\\\\',sep=""),
            paste(xv3, boe18$b[3],'&',boe19$b[3],'\\\\',sep=""),
            paste('&', boe18$s[3],'&',boe19$s[3],'\\\\',sep=""),
            '& & \\\\ \\midrule',
            paste('N Preg (pill) &', NPp18,'&', NPp19 ,'\\\\',sep=""),
            paste('N Preg (close 15) &',NPc1518,'&',NPc1519,'\\\\',sep=""),             
            paste('N Preg (close 30) &',NPc3018,'&',NPc3019,'\\\\',sep=""),             
            'Pills Disbursed & 5,736 & 11,121 \\\\',
            '\\hline \\hline \\\\[-1.8ex]',
            '\\multicolumn{3}{p{8cm}}{\\begin{footnotesize}\\textsc{Notes:} ',
            'Regression coefficients and standard errors are calculated in ',
            'line with specification (\\ref{TEENeqn:spillover}). The number of ',
            'pills disbursed is calculated from administrative data described in ',
            'figure \\ref{TEENfig:Pilltime}, and number of avoided pregnancy is ',
            'based on regression estimates and total births in administrative ',
            'data. Further details are provided in appendix \\ref{TEENscn:BOE}.',
            paste(sig, '\\end{footnotesize}}', sep=""),
            '\\normalsize\\end{tabular}\\end{table}'),to)

close(to)             
}


if(prTab){
    to <-file(paste(tab.dir,"aggregateASFRunweight.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{The Effect of the EC Pill on Birth Rates (Unweighted)}',
                 '\\label{TEENtab:aggregateASFRunweight}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lcccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& Birth& Birth& Birth& Birth\\\\',
                 '& Rate & Rate & Rate & Rate \\\\',
                 '&(1)&(2)&(3)&(4) \\\\ \\hline',
                 '\\multicolumn{5}{l}{\\textbf{',
                 '\\noindent Panel A: All Women}} \\\\',
                 paste('Emergency Contraceptive Pill&',wAll$b,'\\\\',sep=""),
                 paste('            &',wAll$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',wAll$n,'\\\\',sep=""),
                 paste('Mean Birth Rate&',wAll$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 15-19 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',w1519$b,'\\\\',sep=""),
                 paste('            &',w1519$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',w1519$n,'\\\\',sep=""),
                 paste('Mean Birth Rate   &',w1519$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel C: 20-34 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',w2034$b,'\\\\',sep=""),
                 paste('            &',w2034$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',w2034$n,'\\\\',sep=""),
                 paste('Mean Birth Rate   &',w2034$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 35-49 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',w3549$b,'\\\\',sep=""),
                 paste('            &',w3549$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',w3549$n,'\\\\',sep=""),
                 paste('Mean Birth Rate&',w3549$m,'\\\\',sep=""),
                 '\\hline \\\\[-1.8ex] ', 
                 '{\\small Year \\& Comuna FEs}             &Y&Y&Y&Y \\\\',
                 '{\\small Municipal-Specific Linear Trends}& &Y&Y&Y \\\\', 
                 '{\\small Time Varying Controls}           & & &Y&Y \\\\', 
                 '{\\small Spillovers}                      & & & &Y \\\\',
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{5}{p{13.8cm}}{\\begin{footnotesize}'     ,
                 '\\textsc{Notes:} Each panel presents unweighted       ',
                 'difference-in-difference results for a regression of  ',
                 'age-specific fertility rates (ASFR) on the EC reform  ',
                 'for the age group ',
                 'in each municipality.  Specifications are identical to',
                 ' table \\ref{TEENtab:aggregateASFR}, however are not  ',
                 'weighted for population. Standard errors are clustered',
                 'at the level of the municipality.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)

    to <-file(paste(tab.dir,"aggregateASFR.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{The Effect of the EC Pill on Birth Rates}',
                 '\\label{TEENtab:aggregateASFR}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lcccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& Birth& Birth& Birth& Birth\\\\',
                 '& Rate & Rate & Rate & Rate \\\\',
                 '&(1)&(2)&(3)&(4) \\\\ \\hline',
                 '\\multicolumn{5}{l}{\\textbf{',
                 '\\noindent Panel A: All Women}} \\\\',
                 paste('Emergency Contraceptive Pill     &',rAll$b,'\\\\',sep=""),
                 paste('            &',rAll$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',rAll$n,'\\\\',sep=""),
                 paste('Mean Birth Rate&',rAll$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 15-19 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',r1519$b,'\\\\',sep=""),
                 paste('            &',r1519$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',r1519$n,'\\\\',sep=""),
                 paste('Mean Birth Rate&',r1519$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel C: 20-34 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',r2034$b,'\\\\',sep=""),
                 paste('            &',r2034$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',r2034$n,'\\\\',sep=""),
                 paste('Mean Birth Rate&',r2034$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 35-49 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',r3549$b,'\\\\',sep=""),
                 paste('            &',r3549$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',r3549$n,'\\\\',sep=""),
                 paste('Mean Birth Rate&',r3549$m,'\\\\',sep=""),
                 '\\hline \\\\[-1.8ex] ', 
                 '{\\small Year \\& Comuna FEs}             &Y&Y&Y&Y \\\\',
                 '{\\small Municipal-Specific Linear Trends}& &Y&Y&Y \\\\', 
                 '{\\small Time Varying Controls}           & & &Y&Y \\\\', 
                 '{\\small Spillovers}                      & & & &Y \\\\',
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{5}{p{12.8cm}}{\\begin{footnotesize}'       ,
                 '\\textsc{Notes:} Each panel presents population-weighted',
                 ' difference-in-difference results for a regression of '  ,
                 'age-specific fertility rates (ASFR) on the EC reform for',
                 ' the age group in',
                 ' each municipality.  ASFR is defined as the number of   ',
                 'births per 1,000 women.  In the case of all women, this ',
                 'is called the General Fertility Rate (GFR). All models  ',
                 'are estimated by OLS, and each municipalities is        ',
                 'weighted by the population of women. Time varying       ',
                 'controls included in the regression consist of party    ',
                 'dummies for the mayor in power, the mayor\'s gender, the',
                 ' vote margin of the mayor, the percent of girls out of  ',
                 'highschool, education spending spending by both the     ',
                 'municipality and the Ministry of Education, total       ',
                 'health spending and health spending on staff and        ',
                 'training, the percent of female heads of households     ',
                 'living below the poverty line, the percent of female    ',
                 'workers in professional positions in the Municipality,  ',
                 'and condom availability (measured at the level of the   ',
                 'region). Standard errors are clustered at the level of  ',
                 'the municipality.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)

    to <-file(paste(tab.dir,"aggregateBirthsunweight.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{The Effect of the EC Pill on Total Births (Unweighted)}',
                 '\\label{TEENtab:aggregateunweight}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lcccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& Number& Number & Number & Number \\\\',
                 '& Births& Births & Births & Births \\\\',
                 '&(1)&(2)&(3)&(4) \\\\ \\hline',
                 '\\multicolumn{5}{l}{\\textbf{',
                 '\\noindent Panel A: All Women}} \\\\',
                 paste('Emergency Contraceptive Pill&',xAll$b,'\\\\',sep=""),
                 paste('            &',xAll$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',xAll$n,'\\\\',sep=""),
                 paste('Mean Number of Births&',xAll$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 15-19 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',x1519$b,'\\\\',sep=""),
                 paste('            &',x1519$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',x1519$n,'\\\\',sep=""),
                 paste('Mean Number of Births&',x1519$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel C: 20-34 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',x2034$b,'\\\\',sep=""),
                 paste('            &',x2034$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',x2034$n,'\\\\',sep=""),
                 paste('Mean Number of Births&',x2034$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 35-49 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',x3549$b,'\\\\',sep=""),
                 paste('            &',x3549$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',x3549$n,'\\\\',sep=""),
                 paste('Mean Number of Births&',x3549$m,'\\\\',sep=""),
                 '\\hline \\\\[-1.8ex] ', 
                 '{\\small Year \\& Comuna FEs}             &Y&Y&Y&Y \\\\',
                 '{\\small Municipal-Specific Linear Trends}& &Y&Y&Y \\\\', 
                 '{\\small Time Varying Controls}           & & &Y&Y \\\\', 
                 '{\\small Spillovers}                      & & & &Y \\\\',
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{5}{p{14.8cm}}{\\begin{footnotesize}'     ,
                 '\\textsc{Notes:} Each panel presents unweighted       ',
                 'difference-in-difference results for a regression of  ',
                 'the total number of births for the age group in each  ',
                 ' municipality. Specifications are identical to table  ',
                 '\\ref{TEENtab:aggregate}, however unweighted          ',
                 'municipality averages are used. Standard errors are   ',
                 'clustered at the level of the municipality.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)

    to <-file(paste(tab.dir,"aggregateBirths.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{The Effect of the EC Pill on Total Births}',
                 '\\label{TEENtab:aggregate}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lcccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& Number& Number & Number & Number \\\\',
                 '& Births& Births & Births & Births \\\\',
                 '&(1)&(2)&(3)&(4) \\\\ \\hline',
                 '\\multicolumn{5}{l}{\\textbf{',
                 '\\noindent Panel A: All Women}} \\\\',
                 paste('Emergency Contraceptive Pill&',NAll$b,'\\\\',sep=""),
                 paste('            &',NAll$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',NAll$n,'\\\\',sep=""),
                 paste('Mean Number of Births&',NAll$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 15-19 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',N1519$b,'\\\\',sep=""),
                 paste('            &',N1519$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',N1519$n,'\\\\',sep=""),
                 paste('Mean Number of Births&',N1519$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel C: 20-34 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',N2034$b,'\\\\',sep=""),
                 paste('            &',N2034$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',N2034$n,'\\\\',sep=""),
                 paste('Mean Number of Births&',N2034$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 35-49 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',N3549$b,'\\\\',sep=""),
                 paste('            &',N3549$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',N3549$n,'\\\\',sep=""),
                 paste('Mean Number of Births&',N3549$m,'\\\\',sep=""),
                 '\\hline \\\\[-1.8ex] ', 
                 '{\\small Year \\& Comuna FEs}             &Y&Y&Y&Y \\\\',
                 '{\\small Municipal-Specific Linear Trends}& &Y&Y&Y \\\\', 
                 '{\\small Time Varying Controls}           & & &Y&Y \\\\', 
                 '{\\small Spillovers}                      & & & &Y \\\\',
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{5}{p{14.2cm}}{\\begin{footnotesize}'     ,
                 '\\textsc{Notes:} Each panel presents population       ',
                 'weighted difference-in-difference results for a       ',
                 'regression of the total number of births for the age  ',
                 'group in each municipality. Specifications are        ',
                 'identical to table \\ref{TEENtab:aggregateASFR},      ',
                 'however birth weights are replaced by the total number',
                 ' of births. Standard errors are clustered             ',
                 'at the level of the municipality.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)


    to <-file(paste(tab.dir,"aggregateLogsunweight.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{The Effect of the EC Pill on log Births (Unweighted)}',
                 '\\label{TEENtab:aggregateLogunweight}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lcccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& ln(Birth)&ln(Birth)&ln(Birth)&ln(Births) \\\\',
                 '&(1)&(2)&(3)&(4) \\\\ \\hline',
                 '\\multicolumn{5}{l}{\\textbf{',
                 '\\noindent Panel A: All Women}} \\\\',
                 paste('Emergency Contraceptive Pill&',mAll$b,'\\\\',sep=""),
                 paste('            &',mAll$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',mAll$n,'\\\\',sep=""),
                 paste('Mean of ln(Births+1)&',mAll$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 15-19 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',m1519$b,'\\\\',sep=""),
                 paste('            &',m1519$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',m1519$n,'\\\\',sep=""),
                 paste('Mean of ln(Births+1)&',m1519$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel C: 20-34 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',m2034$b,'\\\\',sep=""),
                 paste('            &',m2034$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',m2034$n,'\\\\',sep=""),
                 paste('Mean of ln(Births+1)&',m2034$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 35-49 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',m3549$b,'\\\\',sep=""),
                 paste('            &',m3549$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',m3549$n,'\\\\',sep=""),
                 paste('Mean of ln(Births+1)&',m3549$m,'\\\\',sep=""),
                 '\\hline \\\\[-1.8ex] ', 
                 '{\\small Year \\& Comuna FEs}             &Y&Y&Y&Y \\\\',
                 '{\\small Municipal-Specific Linear Trends}& &Y&Y&Y \\\\', 
                 '{\\small Time Varying Controls}           & & &Y&Y \\\\', 
                 '{\\small Spillovers}                      & & & &Y \\\\',
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{5}{p{13.6cm}}{\\begin{footnotesize}'     ,
                 '\\textsc{Notes:} Each panel presents unweighted       ',
                 'difference-in-difference results for a regression of  ',
                 'the total log(Births+1) for the age group in each     ',
                 ' municipality. Specifications are identical to table  ',
                 '\\ref{TEENtab:aggregateLog}, however unweighted       ',
                 'municipality averages are used. Standard errors are   ',
                 'clustered at the level of the municipality.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)

        to <-file(paste(tab.dir,"aggregateLogBirths.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{The Effect of the EC Pill on log Births}',
                 '\\label{TEENtab:aggregateLog}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lcccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& ln(Birth) & ln(Birth) & ln(Birth) & ln(Birth) \\\\',
                 '&(1)&(2)&(3)&(4) \\\\ \\hline',
                 '\\multicolumn{5}{l}{\\textbf{',
                 '\\noindent Panel A: All Women}} \\\\',
                 paste('Emergency Contraceptive Pill&',lAll$b,'\\\\',sep=""),
                 paste('            &',lAll$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',lAll$n,'\\\\',sep=""),
                 paste('Mean of ln(Births+1) &',lAll$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 15-19 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',l1519$b,'\\\\',sep=""),
                 paste('            &',l1519$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',l1519$n,'\\\\',sep=""),
                 paste('Mean of ln(Births+1) &',l1519$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel C: 20-34 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',l2034$b,'\\\\',sep=""),
                 paste('            &',l2034$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',l2034$n,'\\\\',sep=""),
                 paste('Mean of ln(Births+1) &',l2034$m,'\\\\',sep=""),
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\noindent \\textbf{',
                 'Panel B: 35-49 year olds}} \\\\',
                 paste('Emergency Contraceptive Pill&',l3549$b,'\\\\',sep=""),
                 paste('            &',l3549$s,'\\\\',sep=""),
                 ' & & & & \\\\',
                 paste('Observations&',l3549$n,'\\\\',sep=""),
                 paste('Mean of ln(Births+1) &',l3549$m,'\\\\',sep=""),
                 '\\hline \\\\[-1.8ex] ', 
                 '{\\small Year \\& Comuna FEs}             &Y&Y&Y&Y \\\\',
                 '{\\small Municipal-Specific Linear Trends}& &Y&Y&Y \\\\', 
                 '{\\small Time Varying Controls}           & & &Y&Y \\\\', 
                 '{\\small Spillovers}                      & & & &Y \\\\',
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{5}{p{13.8cm}}{\\begin{footnotesize}'     ,
                 '\\textsc{Notes:} Each panel presents population       ',
                 'weighted difference-in-difference results for a       ',
                 'regression of the log(Births+1) for the age group in  ',
                 'each  municipality. Specifications are identical to   ',
                 'table \\ref{TEENtab:aggregateASFR}, however logs are  ',
                 ' used in place of birth rates. Standard errors are    ',
                 'clustered at the level of the municipality.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)
}

if(ChMund){
    to <-file(paste(tab.dir,"ChamberlainMundlak.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{Logit Results with Chamberlian-Mundlak Device}',
                 '\\label{TEENtab:ChamberlainMundlak}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& 15--19 years & 20--34 years & 35--49 years \\\\',
                 '&(1)&(2)&(3) \\\\ \\hline',
                 ' & & & \\\\',
                 paste(xvar,a1519$CM$b,'&',a2034$CM$b,'&',a3549$CM$b,'\\\\', sep=""), 
                 paste(' &',a1519$CM$s,'&',a2034$CM$s,'&',a3549$CM$s,'\\\\', sep=""), 
                 ' & & & \\\\',
                 paste(obs, a1519$CM$n,'&',a2034$CM$n,'&',a3549$CM$n,'\\\\', sep=""), 
                 paste('$R^2$&',a1519$CM$r,'&',a2034$CM$r,'&',a3549$CM$r,'\\\\', sep=""), 
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{4}{p{9.4cm}}{\\begin{footnotesize}\\textsc{Notes:}',
                 'Regression results estimated using logit with the Chamberlain--%' ,
                 'Mundlak device.  Specification identical to logit estimates from',
                 'column 1 in table \\ref{TEENtab:PillPreg}, using Chamberlain--%' ,
                 'Mundlak rather than comuna fixed effects.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)


    to <-file(paste(tab.dir,"noComunaTrends.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{Estimates without Comuna-Specific Linear Trends}',
                 '\\label{TEENtab:NoTrend}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& 15--19 years & 20--34 years & 35--49 years \\\\',
                 '&(1)&(2)&(3) \\\\ \\hline',
                 ' & & & \\\\',
                 paste(xvar,a1519$NT$b,'&',a2034$NT$b,'&',a3549$NT$b,'\\\\', sep=""), 
                 paste(' &',a1519$NT$s,'&',a2034$NT$s,'&',a3549$NT$s,'\\\\', sep=""), 
                 ' & & & \\\\',
                 paste(obs, a1519$NT$n,'&',a2034$NT$n,'&',a3549$NT$n,'\\\\', sep=""), 
                 paste('$R^2$&',a1519$NT$r,'&',a2034$NT$r,'&',a3549$NT$r,'\\\\', sep=""), 
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{4}{p{9.4cm}}{\\begin{footnotesize}\\textsc{Notes:}' ,
                 'Regression results estimated using comuna fixed effects.        ' ,
                 'Specification identical to logit estimates from column 1 in table',
                 '\\ref{TEENtab:PillPreg}, omitting comuna specific trends.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)

}

if(invPS){
    to <-file(paste(tab.dir,"inversePS.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{Inverse Propensity Score Weighting}',
                 '\\label{TEENtab:inversePS}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& 15--19 years & 20--34 years & 35--49 years \\\\',
                 '&(1)&(2)&(3) \\\\ \\hline',
                 ' & & & \\\\',
                 paste(xvar,ps1519$b,'&',ps2034$b,'&',ps3549$b,'\\\\', sep=""), 
                 paste(' &',ps1519$s,'&',ps2034$s,'&',ps3549$s,'\\\\', sep=""), 
                 ' & & & \\\\',
                 paste(obs, ps1519$n,'&',ps2034$n,'&',ps3549$n,'\\\\', sep=""), 
                 paste('$R^2$&',ps1519$r,'&',ps2034$r,'&',ps3549$r,'\\\\', sep=""), 
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{4}{p{10.4cm}}{\\begin{footnotesize}\\textsc{Notes:}' ,
                 'Regression results estimated using inverse propensity score ',
                 'weighting based on Pr(treatment) on full observables.',
                 'Specifications and controls identical to those described in table',
                 '\\ref{TEENtab:PillPreg}.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)
}


if(robust) {
to <-file(paste(tab.dir,"BirthRobust.tex", sep=""))
writeLines(c('\\begin{landscape}','\\begin{table}[!htbp] \\centering',
             '\\caption{Alternative Specifications -- Births}',
             '\\label{TEENtab:BirthRobust}',
             '\\begin{tabular}{@{\\extracolsep{5pt}}lccccccc}',
             '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
             '&(1)&(2)&(3)&(4)&(5)&(6)&(7) \\\\',
             '&Double&No    & Inv &Chamb--& Full     & OLS   & OLS  \\\\',
             '&Diff. &Trend & PS  &Mundlak& Controls & Count & Rate  \\\\ \\midrule',
             '\\multicolumn{8}{l}{\\textsc{\\noindent All Age Groups}} \\\\',
             ' & & & & & & & \\\\',
             paste(xvar,strsplit(aAll$b,'&')[[1]][1],'&',aAll$NT$b,'&',psAll$b,
                   '&',aAll$CM$b,'&',cAll$b[1],'&',
                   strsplit(NAll$b,'&')[[1]][1],'&',strsplit(rAll$b,'&')[[1]][1],
                   '\\\\',sep=""),
             paste('&' ,strsplit(aAll$se,'&')[[1]][1],'&',aAll$NT$s,'&','&',psAll$s,
                   aAll$CM$s,'&',cAll$s[1],'&',
                   strsplit(NAll$s,'&')[[1]][1],'&',strsplit(rAll$s,'&')[[1]][1],
                   '\\\\',sep=""),
             ' & & & & & & & \\\\',
             paste('R$^2$/Pseudo-R$^2$&',strsplit(aAll$r,'&')[[1]][1],'&',aAll$NT$r,
                   '&',psAll$r,'&',aAll$CM$r,'&',cAll$r[1],'&',
                   strsplit(NAll$r,'&')[[1]][1],'&',strsplit(rAll$r,'&')[[1]][1],
                   '\\\\',sep=""),
             paste('Observations&',strsplit(aAll$n,'&')[[1]][1],'&',aAll$NT$n,'&',
                   psAll$n,'&',aAll$CM$n,'&',cAll$n[1],'&',
                   strsplit(NAll$n,'&')[[1]][1],'&',strsplit(rAll$n,'&')[[1]][1],
                   '\\\\',sep=""),
             '\\multicolumn{8}{l}{\\textsc{\\noindent 15-19 Year-Olds}} \\\\',
             ' & & & & & & & \\\\',
             paste(xvar,strsplit(a1519$b,'&')[[1]][1],'&',a1519$NT$b,'&',ps1519$b,'&',
                   a1519$CM$b,'&',c1519$b[1],'&',
                   strsplit(N1519$b,'&')[[1]][1],'&',strsplit(r1519$b,'&')[[1]][1],
                   '\\\\',sep=""),
             paste('&' ,strsplit(a1519$se,'&')[[1]][1],'&',a1519$NT$s,'&','&',ps1519$s,
                   a1519$CM$s,'&',c1519$s[1],'&',
                   strsplit(N1519$s,'&')[[1]][1],'&',strsplit(r1519$s,'&')[[1]][1],
                   '\\\\',sep=""),
             ' & & & & & & & \\\\',
             paste('R$^2$/Pseudo-R$^2$&',strsplit(a1519$r,'&')[[1]][1],'&',
                   a1519$NT$r,'&',ps1519$r,'&',a1519$CM$r,'&',c1519$r[1],'&',
                   strsplit(N1519$r,'&')[[1]][1],'&',strsplit(r1519$r,'&')[[1]][1],
                   '\\\\',sep=""),
             paste('Observations&',strsplit(a1519$n,'&')[[1]][1],'&',
                   a1519$NT$n,'&',ps1519$n,'&',a1519$CM$n,'&',c1519$n[1],'&',
                   strsplit(N1519$n,'&')[[1]][1],'&',strsplit(r1519$n,'&')[[1]][1],
                   '\\\\',sep=""),
             '\\multicolumn{8}{l}{\\textsc{\\noindent 20-34 Year-Olds}} \\\\',
             ' & & & & & & & \\\\',
             paste(xvar,strsplit(a2034$b,'&')[[1]][1],'&',a2034$NT$b,'&',ps2034$b,'&',
                   a2034$CM$b,'&',c2034$b[1],'&',
                   strsplit(N2034$b,'&')[[1]][1],'&',strsplit(r2034$b,'&')[[1]][1],
                   '\\\\',sep=""),
             paste('&' ,strsplit(a2034$se,'&')[[1]][1],'&',a2034$NT$s,'&',ps2034$s,'&',
                   a2034$CM$s,'&',c2034$s[1],'&',
                   strsplit(N2034$s,'&')[[1]][1],'&',strsplit(r2034$s,'&')[[1]][1],
                   '\\\\',sep=""),
             ' & & & & & & & \\\\',
             paste('R$^2$/Pseudo-R$^2$&',strsplit(a2034$r,'&')[[1]][1],'&',
                   a2034$NT$r,'&',ps2034$r,'&',a2034$CM$r,'&',c2034$r[1],'&',
                   strsplit(N2034$r,'&')[[1]][1],'&',strsplit(r2034$r,'&')[[1]][1],
                   '\\\\',sep=""),
             paste('Observations&',strsplit(a2034$n,'&')[[1]][1],'&',
                   a2034$NT$n,'&',ps2034$n,'&',a2034$CM$n,'&',c2034$n[1],'&',
                   strsplit(N2034$n,'&')[[1]][1],'&',strsplit(r2034$n,'&')[[1]][1],
                   '\\\\',sep=""),
             '\\multicolumn{8}{l}{\\textsc{\\noindent 35-49 Year-Olds}} \\\\',
             ' & & & & & & & \\\\',
             paste(xvar,strsplit(a3549$b,'&')[[1]][1],'&',a3549$NT$b,'&',ps3549$b,'&',
                   a3549$CM$b,'&',c3549$b[1],'&',
                   strsplit(N3549$b,'&')[[1]][1],'&',strsplit(r3549$b,'&')[[1]][1],
                   '\\\\',sep=""),
             paste('&' ,strsplit(a3549$se,'&')[[1]][1],'&',a3549$NT$s,'&',ps3549$s,'&',
                   a3549$CM$s,'&',c3549$s[1],'&',
                   strsplit(N3549$s,'&')[[1]][1],'&',strsplit(r3549$s,'&')[[1]][1],
                   '\\\\',sep=""),
             ' & & & & & & & \\\\',
             paste('R$^2$/Pseudo-R$^2$&',strsplit(a3549$r,'&')[[1]][1],'&',
                   a3549$NT$r,'&',ps3549$r,'&',a3549$CM$r,'&',c3549$r[1],'&',
                   strsplit(N3549$r,'&')[[1]][1],'&',strsplit(r3549$r,'&')[[1]][1],
                   '\\\\',sep=""),
             paste('Observations&',strsplit(a3549$n,'&')[[1]][1],'&',
                   a3549$NT$n,'&',ps3549$n,'&',a3549$CM$n,'&',c3549$n[1],'&',
                   strsplit(N3549$n,'&')[[1]][1],'&',strsplit(r3549$n,'&')[[1]][1],
                   '\\\\',sep=""),
             '\\hline \\\\[-1.8ex] ',
             '\\multicolumn{8}{p{17.2cm}}{\\begin{footnotesize}\\textsc{Notes:}  ',
             'Column 1 replicates the diff-in-diff result from table             ',
             '\\ref{TEENtab:PillPreg}. Column 2 estimates the same specification,',
             ' however without municipality fixed effects. Column 3 weights using',
             ' inverse propensity score, where propensity scores are estimated   ',
             'using a probit and full controls are described in table            ',
             '\\ref{TEENtab:PillPreg}.  Column 4 uses the Chamberlain-Mundlak    ',
             'device where municipal fixed effects are replaced with             ',
             'municipal-level mean independent variables, and column 5 presents  ',
             'regressions using full controls and correcting for spillovers (see ',
             'table \\ref{TEENtab:Spillover}). Columns 6 and 7 are estimated at  ',
             'the municipal level using continuous outcome variables: the count  ',
             'of pregnancies, or the rate of pregnancies/all women.  Full results',
             'for these specifications are included in online appendix tables.   ',
             paste(sig, '\\end{footnotesize}}', sep=""),
             '\\normalsize\\end{tabular}\\end{table}\\end{landscape}'),to)
close(to)
}
