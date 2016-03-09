# DeathEstimates.R v1.23          KEL / DCC               yyyy-mm-dd:2013-12-29
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# Import data from S1Data_granular_covars.csv to run various models.  Initial 
# data contains pregnancies by comuna and morning after pill availability by
# comuna.
# 
# Principal models are of the form:
#  P(death)_{ijt} = a + beta*PAE_{jt} + gamma*comuna_j + delta*year_t + u_{ijt}
#  deaths_{jt}   = a + beta*PAE_{jt} + gamma*comuna_j + delta*year_t + u_{jt}
#
# This code has been written by KEL, with updates by DCC to incorporate 
# additional time varying controls and export results to TeX.  When running the
# switches in section (1) determine whether or not specific sections of the
# code will be run.
#
# Last version -- 1.23: Restructure directories
# contact: damian.clarke@economics.ox.ac.uk


rm(list=ls())

#******************************************************************************
#***(1) Parameters
#******************************************************************************
create   <- FALSE
death    <- FALSE
Ndeath   <- T
Ldeath   <- FALSE
deathTab <- T
spill    <- FALSE
combine  <- FALSE
events   <- FALSE
PSweight <- FALSE

birth_y_range <- 2006:2012
pill_y_range  <- birth_y_range - 1
age_range     <- c(15,49)
week          <- 20
pat           <- "P"

#******************************************************************************
#***(2) Libraries, directories
#******************************************************************************
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if (length(new.pkg))
     	install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

pk <- c("xtable","rms","sandwich","lmtest","plyr","foreign","reshape","gdata")
ipak(pk)

proj.dir <- "~/universidades/Oxford/DPhil/Thesis/Teens/"

brth.dir <- paste(proj.dir, "Data/Nacimientos/", sep="")
code.dir <- paste(proj.dir, "Source/Deaths/"   , sep="")
com.dir  <- paste(proj.dir, "Data/Comunas/"    , sep="")
deth.dir <- paste(proj.dir, "Data/Deaths/"     , sep="")
ma.dir   <- paste(proj.dir, "Data/PAE/"        , sep="")
pol.dir  <- paste(proj.dir, "Data/Alcaldes/"   , sep="")
pop.dir  <- paste(proj.dir, "Data/Population/" , sep="")
tab.dir  <- paste(proj.dir, "Tables/"          , sep="")
graf.dir <- paste(proj.dir, "Figures/"         , sep="")

Names <- c("dom_comuna","trend","trend2","pill","mujer","party","votes"      ,
           "outofschool","healthspend","healthstaff","healthtraining"        , 
           "educationspend","femalepoverty","year","urban","educationmunicip",
           "condom","usingcont","femaleworkers","poverty","popln")


#******************************************************************************
#***(3a) Source functions
#******************************************************************************
if(create) {
    f <- paste(code.dir,"DeathGenerate.R",sep="")
    source(f)

    filename <- paste(deth.dir, 'S1Data_deaths_covars.csv', sep="")
    prep_s1_data_deaths(age_range,week,pat,TRUE,filename)
}

#******************************************************************************
#***(3b) Load Data
#******************************************************************************
f <- paste(deth.dir,"S1Data_deaths_covars.csv",sep="")
orig <- read.csv(f)

#******************************************************************************
#***(4a) Editing functions
#******************************************************************************
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
  
    if(dim==1|dim==10) {
        beta <- summary(outresults)$coefficients[pillline,]["Estimate"]
        se   <- outresults$coefficients2[pillline,]["Std. Error"]
        if (dim==1) {p <- outresults$coefficients2[pillline,]["Pr(>|z|)"]}
    else {p <- outresults$coefficients2[pillline,]["Pr(>|t|)"]}
    }
    else {
        beta <- summary(outresults)$coefficients[pillline,][, "Estimate"]
        se   <- outresults$coefficients2[pillline,][, "Std. Error"]
        p    <- outresults$coefficients2[pillline,][, "Pr(>|z|)"]    
    }
  
    if(dim==1|dim==2) {
        null  <- glm(cbind(successes,failures) ~ 1,
                     family="binomial",data=d,weights=WT)
        Lfull <- as.numeric(logLik(outresults))
        Lnull <- as.numeric(logLik(null))
        R2    <- 1 - Lfull/Lnull
    }
    if(dim==10) {
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

datcollapse <- function(age_sub,deathtype,dat) {
    dat        <- dat[dat$age %in% age_sub,]
    dat        <- dat[dat$pregnant == 1,]
    dat$popln  <- ave(dat$n,dat$dom_comuna,dat$year,FUN=sum)
    dat$Q      <- dat$earlyQ+dat$lateQ
    dat$early  <- dat$earlyQ+dat$earlyP
    dat$late   <- dat$lateQ+dat$lateP
    
    dat <- closegen(0,15,dat)
    dat <- closegen(15,30,dat)
    dat <- closegen(30,45,dat)
    dat <- dat[complete.cases(dat),]
  
    dat <- aggregate.data.frame(dat[,c("n",deathtype)],
                                by=list(dat$close15,dat$close30,dat$close45,
                                    dat$dom_comuna,dat$year-2005           ,
                                    (dat$year-2005)^2,dat$pill,dat$mujer   ,
                                    dat$party,dat$votop,dat$outofschool    ,
                                    dat$healthspend,dat$healthstaff        ,
                                    dat$healthtraining,dat$educationspend  ,
                                    dat$femalepoverty,dat$urbBin,dat$year  ,
                                    dat$educationmunicip,dat$condom        ,
                                    dat$usingcont,dat$femaleworkers        ,
                                    dat$poverty,dat$popln),
                                function(vec) {sum(na.omit(vec))})
    names(dat)    <- c("close15","close30","close45",Names,"n","death")
    dat$failures  <- dat$n 
    dat$successes <- dat$death
    
    mod <- aggregate.data.frame(dat[,c("failures","successes")],
                                by=list(dat$close15,dat$close30,dat$close45,
                                    dat$dom_comuna,dat$trend,dat$trend2    ,
                                    dat$pill,dat$mujer,dat$party,dat$votes ,
                                    dat$outofschool,dat$healthspend        ,
                                    dat$healthstaff,dat$healthtraining     ,
                                    dat$educationspend,dat$femalepoverty   ,
                                    dat$year,dat$urban,dat$educationmunicip,
                                    dat$condom,dat$usingcont,
                                    dat$femaleworkers,dat$poverty,dat$popln),
                                function(vec) {sum(na.omit(vec))})
    names(mod) <- c("close15","close30","close45",Names,"failures","successes")
    return(mod)
}

inversePS <- function(age_sub,deathtype,dat) {
    #CALCULATE PROPENSITY SCORE WEIGHT FOR EACH COMUNA BASED ON
    #PILL/NON-PILL STATUS.  USE FULL SET OF CONTROLS TO PREDICT PSCORE.
    
    psdat <- datcollapse(age_sub,deathtype,dat)
    psdat$conservative[psdat$party=="UDI"|psdat$party=="RN"] <- 1
    psdat$conservative[is.na(psdat$conservative)] <- 0
    
    psdat <- aggregate.data.frame(psdat[,c("pill","mujer","conservative"  ,
                                           "votes","outofschool"          ,
                                           "healthspend","healthstaff"    ,
                                           "healthtraining","usingcont"   ,
                                           "femalepoverty","femaleworkers",
                                           "condom","educationspend"      ,
                                           "educationmunicip","popln")],
                                  by=list(psdat$dom_comuna), FUN=mean)
    psdat$pillind[psdat$pill>0]<-1
    psdat$pillind[is.na(psdat$pillind)]<-0
    
    PSc <- glm(pillind ~ mujer + conservative + votes + outofschool +
               healthspend + healthstaff + healthtraining           +
               usingcont + femalepoverty + condom + femaleworkers   +
               educationspend+ educationmunicip + popln,
               family=binomial, data=psdat)
    psdat$predict<-predict(PSc,type="response")
    psdat$WT[psdat$pillind==1] <- 1/psdat$predict
    psdat$WT[psdat$pillind==0] <- 1/(1-psdat$predict)
    
    colnames(psdat)[1]<-"dom_comuna"
    wts <- psdat[,(names(psdat) %in% c("dom_comuna", "WT"))]
    return(wts)
}


#******************************************************************************
#***(4b) Main Functions
#******************************************************************************
death_pmod <- function(age_sub,deathtype,numb,PSwt) {

    formod <- datcollapse(age_sub,deathtype,orig)
    if(PSwt) {
        wts <- inversePS(age_sub,deathtype,orig)
        formod <- join(formod,wts,by="dom_comuna",type="left",match="all")

        xPS <- glm(cbind(successes,failures) ~ factor(year) + factor(pill)    +
                       factor(dom_comuna) + factor(dom_comuna):trend,
                   family="binomial", data=formod, weights=WT)
        
        clusters <- mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
        xPS$coefficients2  <- robust.se(xPS,clusters)[[2]]
        n  <- sum(formod$successes) + sum(formod$failures)
        
        s0 <- pillest(xPS, formod, n, "pill", 1)
        return(s0)
        
    } else {
        formod$WT <- 1
    }
    
    formod$meanP  <- ave(formod$pill, group=formod$dom_comuna) 
    
    x  <- glm(cbind(successes,failures) ~ factor(dom_comuna) + factor(year) +
              factor(dom_comuna):trend + factor(pill) + factor(party)       + 
              factor(mujer) + votes + outofschool + educationspend          + 
              educationmunicip + healthspend + healthtraining + healthstaff +
              femalepoverty + femaleworkers + condom,
              family="binomial", data=formod)
    xCM <- glm(cbind(successes,failures) ~ meanP + factor(year)              +
               factor(dom_comuna):trend + factor(pill) + factor(party)       + 
               factor(mujer) + votes + outofschool + educationspend          + 
               educationmunicip + healthspend + healthtraining + healthstaff +
               femalepoverty + femaleworkers + condom,
               family="binomial", data=formod)

    if (numb==1) {
        clusters <-mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
        x$coefficients2   <- robust.se(x,clusters)[[2]]
        xCM$coefficients2 <- robust.se(xCM,clusters)[[2]]
        
        n  <- sum(formod$successes) + sum(formod$failures)
        s1 <- pillest(x, formod, n, "pill", 1)
        s2 <- pillest(xCM, formod, n, "pill", 1)
        db <- sum(formod$successes)/sum(formod$failures)
        db <- format(round(db,3),nsmall=3)
        c  <- sum(formod$successes)
        c  <- format(c,big.mark=",",scientific=F)

        return(list("beta" = s1$b, "se" = s1$s, "R2" = s1$r, "n" = s1$n, "c" = c,
                    "db" = db, "CM" = s2))
    } else {
        return(x)
    }
}

NumDeath <- function(age_sub,deathtype,R,wt)  {
    dat <- orig
    dat <- dat[dat$age %in% age_sub,]
    dat <- dat[dat$pregnant == 1,]
    dat$popln  <- ave(dat$n,dat$dom_comuna,dat$year,FUN=sum)
    
    dat <- aggregate.data.frame(dat[,c("n",deathtype)],
                                by=list(dat$dom_comuna,dat$year-2005       ,
                                    (dat$year-2005)^2,dat$pill,dat$mujer   ,
                                    dat$party,dat$votop,dat$outofschool    ,
                                    dat$healthspend,dat$healthstaff        ,
                                    dat$healthtraining,dat$educationspend  ,
                                    dat$femalepoverty,dat$urbBin,dat$year  ,
                                    dat$educationmunicip,dat$condom        ,
                                    dat$usingcont,dat$femaleworkers        ,
                                    dat$poverty,dat$popln)                 ,
                                function(vec) {sum(na.omit(vec))})
    names(dat) <- c(Names,"n","death")
    dat$weight <- 1

    if (R)    dat$FDb    <- (dat$death/dat$n)*1000
    if (wt)   dat$weight <- dat$popln
    
    
    xFl  <- lm(FDb ~ factor(dom_comuna) + factor(year)                       +
               factor(dom_comuna):trend + factor(pill) + factor(party)       + 
               factor(mujer) + votes + outofschool + educationspend          + 
               educationmunicip + healthspend + healthtraining + healthstaff +
               femalepoverty + femaleworkers + condom,
               data=dat, weights=weight)
    clusters <-mapply(paste,"dom_comuna.",dat$dom_comuna,sep="")
    xFl$coefficients2   <- robust.se(xFl,clusters)[[2]]
    
    n  <- nrow(dat)
    s2 <- pillest(xFl, dat, n, "pill", 10)
    mn <- format(round(wtd.mean(dat$FDb,weights=dat$popln),2),nsmall=2) 
    return(list("b" = s2$b, "se" = s2$s, "R2" = s2$r, "n" = s2$n, "m"=mn))
}

#==============================================================================
#===(4c) Various functions to examine effect of spillover 
#==============================================================================
closegen <- function(d1,d2,dat) {
    dat2 <- dat
    
    dat2$newvar <- NA
    dat2$newvar[dat2$pilldistance > d1 & dat2$pilldistance <= d2 &
                !(dat2$pilldistance)==0] <- 1
    dat2$newvar[is.na(dat2$newvar)]<-0
    
    names(dat2)<-c(names(dat),paste('close',d2,sep=""))
    return(dat2)
}

spillovers <- function(age_sub,deathtype,road,time) {
    mod    <- datcollapse(age_sub,deathtype,orig)
    mod$WT <- 1
    if(road) {
        mod$close15 <- mod$road15
        mod$close30 <- mod$road30
        mod$close45 <- mod$road45
    }
    if(time) {
        mod$close15 <- mod$time15
        mod$close30 <- mod$time30
        mod$close45 <- mod$time45
    }
    
    xspill <- glm(cbind(successes,failures) ~ factor(dom_comuna) + factor(year) +
                  factor(dom_comuna):trend + factor(party) + factor(pill)     + 
                  factor(mujer) + votes + outofschool + educationspend        + 
                  educationmunicip + healthspend + healthtraining + condom    +
                  healthstaff + femalepoverty + femaleworkers                 + 
                  factor(close15) + factor(close30) + factor(close45),
                family="binomial",data=mod)

    clusters <-mapply(paste,"dom_comuna.",mod$dom_comuna,sep="")
    xspill$coefficients2 <- robust.se(xspill,clusters)[[2]]
    
    n  <- sum(mod$successes) + sum(mod$failures)
    s1 <- pillest(xspill,mod,n,"pill|close",2)
  
    return(list("b" = s1$b,"s" = s1$s, "n" = s1$n, "r" = s1$r))
}

#==============================================================================
#===(4d) Event study
#==============================================================================
event <- function(age_sub,deathtype) {
    dat <- orig
    dat <- dat[dat$age %in% age_sub,]
    dat <- dat[dat$pregnant == 1,]
  
    formod <- aggregate.data.frame(dat[,c("n",deathtype)],
                                   by=list(dat$dom_comuna,dat$year-2005        ,
                                       (dat$year-2005)^2,dat$pill,dat$mujer    ,
                                       dat$party,dat$votop,dat$outofschool     ,
                                       dat$healthspend,dat$healthstaff         ,
                                       dat$healthtraining,dat$educationspend   ,
                                       dat$femalepoverty,dat$urbBin,dat$year   ,
                                       dat$educationmunicip,dat$condom         ,
                                       dat$usingcont,dat$femaleworkers         ,
                                       dat$poverty,dat$popln)                  ,
                                   function(vec) {sum(na.omit(vec))})
    names(formod) <- c(Names,"n","death")

    formod$FDbirth <- (formod$death/formod$n)*1000

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

    eventS  <- lm(FDbirth ~ factor(year)                                         +
                   factor(dom_comuna) + factor(dom_comuna):trend + votes         +
                   factor(party) + factor(mujer) + outofschool + educationspend  +
                   educationmunicip + healthspend + healthtraining + healthstaff +
                   femalepoverty + femaleworkers + condom + factor(pilln5)       +
                   factor(pilln4) + factor(pilln2) + factor(pilln1)              +
                   factor(pillp0) + factor(pillp1) + factor(pillp2),  data=formod)
    clusters <-mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
    eventS$coefficients2 <- robust.se(eventS,clusters)[[2]]


    pillline <- grepl("pill",rownames(summary(eventS)$coefficients))
    beta <- summary(eventS)$coefficients[pillline,][, "Estimate"]
    se   <- eventS$coefficients2[pillline,][, "Std. Error"]
    
    return(list("b" = beta, "s" = se, "eventyr" = c(-5,-4,-2,-1,0,1,2)))    
}


#******************************************************************************
#***(5) Estimate
#******************************************************************************
if(death){
    print("Prob Death")
    pAll  <- death_pmod(15:49, "death",1,PSwt=FALSE)
    p1519 <- death_pmod(15:19, "death",1,PSwt=FALSE)
    p2034 <- death_pmod(20:34, "death",1,PSwt=FALSE)
    p3549 <- death_pmod(35:49, "death",1,PSwt=FALSE)
    eAll  <- death_pmod(15:49,"earlyP",1,PSwt=FALSE)
    e1519 <- death_pmod(15:19,"earlyP",1,PSwt=FALSE)
    e2034 <- death_pmod(20:34,"earlyP",1,PSwt=FALSE)
    e3549 <- death_pmod(35:49,"earlyP",1,PSwt=FALSE)
    lAll  <- death_pmod(15:49, "lateP",1,PSwt=FALSE)
    l1519 <- death_pmod(15:19, "lateP",1,PSwt=FALSE)
    l2034 <- death_pmod(20:34, "lateP",1,PSwt=FALSE)
    l3549 <- death_pmod(35:40, "lateP",1,PSwt=FALSE)
}

if(Ndeath) {
    print("Rate Death (weighted)")
    wAll   <- NumDeath(15:49, "death",T,F)
    w1519  <- NumDeath(15:19, "death",T,F)
    w2034  <- NumDeath(20:34, "death",T,F)
    w3549  <- NumDeath(35:49, "death",T,F)
    weAll  <- NumDeath(15:49,"earlyP",T,F)
    we1519 <- NumDeath(15:19,"earlyP",T,F)
    we2034 <- NumDeath(20:34,"earlyP",T,F)
    we3549 <- NumDeath(35:49,"earlyP",T,F)
    wlAll  <- NumDeath(15:49, "lateP",T,F)
    wl1519 <- NumDeath(15:19, "lateP",T,F)
    wl2034 <- NumDeath(20:34, "lateP",T,F)
    wl3549 <- NumDeath(35:40, "lateP",T,F)

    print("Rate Death (unweighted)")
    xAll   <- NumDeath(15:49, "death",T,T)
    x1519  <- NumDeath(15:19, "death",T,T)
    x2034  <- NumDeath(20:34, "death",T,T)
    x3549  <- NumDeath(35:49, "death",T,T)
    xeAll  <- NumDeath(15:49,"earlyP",T,T)
    xe1519 <- NumDeath(15:19,"earlyP",T,T)
    xe2034 <- NumDeath(20:34,"earlyP",T,T)
    xe3549 <- NumDeath(35:49,"earlyP",T,T)
    xlAll  <- NumDeath(15:49, "lateP",T,T)
    xl1519 <- NumDeath(15:19, "lateP",T,T)
    xl2034 <- NumDeath(20:34, "lateP",T,T)
    xl3549 <- NumDeath(35:40, "lateP",T,T)

    print("Num Death (weighted)")
    nAll   <- NumDeath(15:49, "death",T,F)
    n1519  <- NumDeath(15:19, "death",T,F)
    n2034  <- NumDeath(20:34, "death",T,F)
    n3549  <- NumDeath(35:49, "death",T,F)
    neAll  <- NumDeath(15:49,"earlyP",T,F)
    ne1519 <- NumDeath(15:19,"earlyP",T,F)
    ne2034 <- NumDeath(20:34,"earlyP",T,F)
    ne3549 <- NumDeath(35:49,"earlyP",T,F)
    nlAll  <- NumDeath(15:49, "lateP",T,F)
    nl1519 <- NumDeath(15:19, "lateP",T,F)
    nl2034 <- NumDeath(20:34, "lateP",T,F)
    nl3549 <- NumDeath(35:40, "lateP",T,F)

    print("Rate Death (unweighted)")
    oAll   <- NumDeath(15:49, "death",T,T)
    o1519  <- NumDeath(15:19, "death",T,T)
    o2034  <- NumDeath(20:34, "death",T,T)
    o3549  <- NumDeath(35:49, "death",T,T)
    oeAll  <- NumDeath(15:49,"earlyP",T,T)
    oe1519 <- NumDeath(15:19,"earlyP",T,T)
    oe2034 <- NumDeath(20:34,"earlyP",T,T)
    oe3549 <- NumDeath(35:49,"earlyP",T,T)
    olAll  <- NumDeath(15:49, "lateP",T,T)
    ol1519 <- NumDeath(15:19, "lateP",T,T)
    ol2034 <- NumDeath(20:34, "lateP",T,T)
    ol3549 <- NumDeath(35:40, "lateP",T,T)
}

if(spill){
    print("Spillover Models")
    sAll  <- spillovers(age_sub = 15:49,"earlyP",F,F)
    s1519 <- spillovers(age_sub = 15:19,"earlyP",F,F)
    s2034 <- spillovers(age_sub = 20:34,"earlyP",F,F)
    s3549 <- spillovers(age_sub = 35:49,"earlyP",F,F)
}

#******************************************************************************
#***(6) Export
#******************************************************************************
xvar <- 'Morning After Pill &'
xvn  <- 'Emergency Contraceptive Pill &'
xv2  <- 'Close $<15$ km &'
xv3  <- 'Close $30-60$ km &'
obs  <- 'Observations&'
R2   <- 'R-Squared'
sig  <- '$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01;'
dpb  <- 'Mean (fetal deaths/live birth)&'
s    <- '\\\\'
a    <- '&'

if(deathTab) {
    to <- file(paste(tab.dir,"DeathsPerBirth.tex", sep=""))
    writeLines(c('\\begin{table}[htpb!] \\centering'                        ,
                 '\\caption{The Effect of the EC Pill on Fetal Death Rates}',
                 '\\label{TEENtab:DeathOLS}'                                ,
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lccc}\\\\[-1.8ex]'  ,
                 '\\hline\\hline\\\\[-1.8ex]'                               ,
                 '& All    & Early     & Late      \\\\'                    ,
                 '& Deaths & Gestation & Gestation \\\\ \\midrule'          ,
                 '\\multicolumn{4}{l}{\\noindent \\textbf{'                 ,
                 'Panel A: All Women}} \\\\'                                ,
                 paste(xvn, wAll$b, a, weAll$b, a, wlAll$b, s, sep="")      ,
                 paste(a  , wAll$s, a, weAll$s, a, wlAll$s, s, sep="")      ,
                 '& & & \\\\'                                               ,
                 paste(obs, wAll$n, a, weAll$n, a, wlAll$n, s, sep="")      ,
                 paste(dpb, wAll$m, a, weAll$m, a, wlAll$m, s, sep="")      ,
                 '&&&\\\\'                                                  ,
                 '\\multicolumn{4}{l}{\\noindent \\textbf{'                 ,
                 'Panel B: 15-19 year olds}} \\\\'                          ,
                 paste(xvn, w1519$b, a, we1519$b, a, wl1519$b, s, sep="")   ,
                 paste(a  , w1519$s, a, we1519$s, a, wl1519$s, s, sep="")   ,
                 '& & & \\\\'                                               ,
                 paste(obs, w1519$n, a, we1519$n, a, wl1519$n, s, sep="")   ,
                 paste(dpb, w1519$m, a, we1519$m, a, wl1519$m, s, sep="")   ,
                 '&&&\\\\'                                                  ,
                 '\\multicolumn{4}{l}{\\noindent \\textbf{'                 ,
                 'Panel C: 20-34 year olds}} \\\\'                          ,
                 paste(xvn, w2034$b, a, we2034$b, a, wl2034$b, s, sep="")   ,
                 paste(a  , w2034$s, a, we2034$s, a, wl2034$s, s, sep="")   ,
                 '& & & \\\\'                                               ,
                 paste(obs, w2034$n, a, we2034$n, a, wl2034$n, s, sep="")   ,
                 paste(dpb, w2034$m, a, we2034$m, a, wl2034$m, s, sep="")   ,
                 '&&&\\\\'                                                  ,
                 '\\multicolumn{4}{l}{\\noindent \\textbf{'                 ,
                 'Panel C: 35-49 year olds}} \\\\'                          ,
                 paste(xvn, w3549$b, a, we3549$b, a, wl3549$b, s, sep="")   ,
                 paste(a  , w3549$s, a, we3549$s, a, wl3549$s, s, sep="")   ,
                 '& & & \\\\'                                               ,
                 paste(obs, w3549$n, a, we3549$n, a, wl3549$n, s, sep="")   ,
                 paste(dpb, w3549$m, a, we3549$m, a, wl3549$m, s, sep="")   ,
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{4}{p{10.2cm}}{\\begin{footnotesize}        ',
                 '\\textsc{Notes:} Each panel presents weighted difference-',
                 'in-difference results for a regression of the fetal      ',
                 'death rate (deaths per 1,000 live births) on the EC      ',
                 'reform for the age group in question. All models are     ',
                 'estimated by OLS, and each municipality is weighted by   ',
                 'the number of live births. All regressions include the   ',
                 'controls documented in table \\ref{TEENtab:aggregateASFR}.',
                 'Standard errors are clustered at the level of the        ',
                 'municipality.'                                            ,
                 paste(sig,'\\end{footnotesize}}',sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)

    to <- file(paste(tab.dir,"DeathsPerBirthunweight.tex", sep=""))
    writeLines(c('\\begin{table}[htpb!] \\centering'                        ,
                 '\\caption{The Effect of EC on Fetal Deaths (Unweighted) }',
                 '\\label{TEENtab:DeathOLSunweight}'                        ,
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lccc}\\\\[-1.8ex]'  ,
                 '\\hline\\hline\\\\[-1.8ex]'                               ,
                 '& All    & Early     & Late      \\\\'                    ,
                 '& Deaths & Gestation & Gestation \\\\ \\midrule'          ,
                 '\\multicolumn{4}{l}{\\noindent \\textbf{'                 ,
                 'Panel A: All Women}} \\\\'                                ,
                 paste(xvn, xAll$b, a, xeAll$b, a, xlAll$b, s, sep="")      ,
                 paste(a  , xAll$s, a, xeAll$s, a, xlAll$s, s, sep="")      ,
                 '& & & \\\\'                                               ,
                 paste(obs, xAll$n, a, xeAll$n, a, xlAll$n, s, sep="")      ,
                 paste(dpb, xAll$m, a, xeAll$m, a, xlAll$m, s, sep="")      ,
                 '&&&\\\\'                                                  ,
                 '\\multicolumn{4}{l}{\\noindent \\textbf{'                 ,
                 'Panel B: 15-19 year olds}} \\\\'                          ,
                 paste(xvn, x1519$b, a, xe1519$b, a, xl1519$b, s, sep="")   ,
                 paste(a  , x1519$s, a, xe1519$s, a, xl1519$s, s, sep="")   ,
                 '& & & \\\\'                                               ,
                 paste(obs, x1519$n, a, xe1519$n, a, xl1519$n, s, sep="")   ,
                 paste(dpb, x1519$m, a, xe1519$m, a, xl1519$m, s, sep="")   ,
                 '&&&\\\\'                                                  ,
                 '\\multicolumn{4}{l}{\\noindent \\textbf{'                 ,
                 'Panel C: 20-34 year olds}} \\\\'                          ,
                 paste(xvn, x2034$b, a, xe2034$b, a, xl2034$b, s, sep="")   ,
                 paste(a  , x2034$s, a, xe2034$s, a, xl2034$s, s, sep="")   ,
                 '& & & \\\\'                                               ,
                 paste(obs, x2034$n, a, xe2034$n, a, xl2034$n, s, sep="")   ,
                 paste(dpb, x2034$m, a, xe2034$m, a, xl2034$m, s, sep="")   ,
                 '&&&\\\\'                                                  ,
                 '\\multicolumn{4}{l}{\\noindent \\textbf{'                 ,
                 'Panel C: 35-49 year olds}} \\\\'                          ,
                 paste(xvn, x3549$b, a, xe3549$b, a, xl3549$b, s, sep="")   ,
                 paste(a  , x3549$s, a, xe3549$s, a, xl3549$s, s, sep="")   ,
                 '& & & \\\\'                                               ,
                 paste(obs, x3549$n, a, xe3549$n, a, xl3549$n, s, sep="")   ,
                 paste(dpb, x3549$m, a, xe3549$m, a, xl3549$m, s, sep="")   ,
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{4}{p{11.4cm}}{\\begin{footnotesize}        ',
                 '\\textsc{Notes:} Each panel presents unweighted          ',
                 'difference-in-difference results for a regression of the ',
                 'fetal death rate (deaths per 1,000 live births) on the EC',
                 'reform for the age group in question. Refer to table     ',
                 '\\ref{TEENtab:DeathOLS} for full notes.  Standard errors ',
                 'are clustered at the level of the municipality.          ',
                 paste(sig,'\\end{footnotesize}}',sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)
}


if(spill){  
  to <- file(paste(tab.dir,"Spillovers_B.tex", sep=""))
  writeLines(c('\\multicolumn{5}{l}{\\textsc{\\noindent Panel B: Fetal Deaths}}\\\\',
               '&&&\\\\',
               paste(xvar,sAll$b[1],a,s1519$b[1],a,s2034$b[1],a,s3549$b[1],s,sep=""),
               paste('&' ,sAll$s[1],a,s1519$s[1],a,s2034$s[1],a,s3549$s[1],s,sep=""),               
               paste(xv2 ,sAll$b[2],a,s1519$b[2],a,s2034$b[2],a,s3549$b[2],s,sep=""),
               paste('&' ,sAll$s[2],a,s1519$s[2],a,s2034$s[2],a,s3549$s[2],s,sep=""),
               '&&&\\\\',
               paste(obs,sAll$n,a,s1519$n,a,s2034$n,a,s3549$n,s,sep=""),
               paste(R2,sAll$r,a,s1519$r,a,s2034$r,a,s3549$r,s,sep=""),
               '\\hline \\hline \\\\[-1.8ex]',
               '\\multicolumn{5}{p{11.2cm}}{\\begin{footnotesize}\\textsc{Notes:}',
               'All models are estimated using logistic regressions, and',
               'coefficients are reported as log odds.  Each regression includes',
               'comuna and year fixed effects and comuna-specific trends, and',
               'the full set of time-varying controls described in table',
               '\\ref{TEENtab:PillPreg}.  \\citet{Conley1999} standard errors are',
               'reported.',
               paste(sig,'\\end{footnotesize}}',sep=""),
               '\\normalsize\\end{tabular}\\end{table}'),to)
  
  close(to)

  to <- file(paste(tab.dir,"SpilloversROAD_B.tex", sep=""))
  writeLines(c('\\multicolumn{5}{l}{\\textsc{\\noindent Panel B: Fetal Deaths}}\\\\',
               '&&&\\\\',
               paste(xvar,sAll$b[1],a,s1519$b[1],a,s2034$b[1],a,s3549$b[1],s,sep=""),
               paste('&' ,sAll$s[1],a,s1519$s[1],a,s2034$s[1],a,s3549$s[1],s,sep=""),               
               paste(xv2 ,sAll$b[2],a,s1519$b[2],a,s2034$b[2],a,s3549$b[2],s,sep=""),
               paste('&' ,sAll$s[2],a,s1519$s[2],a,s2034$s[2],a,s3549$s[2],s,sep=""),
               '&&&\\\\',
               paste(obs,sAll$n,a,s1519$n,a,s2034$n,a,s3549$n,s,sep=""),
               paste(R2,sAll$r,a,s1519$r,a,s2034$r,a,s3549$r,s,sep=""),
               '\\hline \\hline \\\\[-1.8ex]',
               '\\multicolumn{5}{p{11.8cm}}{\\begin{footnotesize}\\textsc{Notes:}',
               'All models are estimated using logistic regressions, and',
               'coefficients are reported as log odds.  Each regression includes',
               'comuna and year fixed effects and comuna-specific trends, and',
               'the full set of time-varying controls described in table',
               '\\ref{TEENtab:PillPreg}.  \\citet{Conley1999} standard errors are',
               'reported.',
               paste(sig,'\\end{footnotesize}}',sep=""),
               '\\normalsize\\end{tabular}\\end{table}'),to)
  
  close(to)

  to <- file(paste(tab.dir,"SpilloversTIME_B.tex", sep=""))
  writeLines(c('\\multicolumn{5}{l}{\\textsc{\\noindent Panel B: Fetal Deaths}}\\\\',
               '&&&\\\\',
               paste(xvar,sAll$b[1],a,s1519$b[1],a,s2034$b[1],a,s3549$b[1],s,sep=""),
               paste('&' ,sAll$s[1],a,s1519$s[1],a,s2034$s[1],a,s3549$s[1],s,sep=""),               
               paste(xv2 ,sAll$b[2],a,s1519$b[2],a,s2034$b[2],a,s3549$b[2],s,sep=""),
               paste('&' ,sAll$s[2],a,s1519$s[2],a,s2034$s[2],a,s3549$s[2],s,sep=""),
               '&&&\\\\',
               paste(obs,sAll$n,a,s1519$n,a,s2034$n,a,s3549$n,s,sep=""),
               paste(R2,sAll$r,a,s1519$r,a,s2034$r,a,s3549$r,s,sep=""),
               '\\hline \\hline \\\\[-1.8ex]',
               '\\multicolumn{5}{p{11.8cm}}{\\begin{footnotesize}\\textsc{Notes:}',
               'All models are estimated using logistic regressions, and',
               'coefficients are reported as log odds.  Each regression includes',
               'comuna and year fixed effects and comuna-specific trends, and',
               'the full set of time-varying controls described in table',
               '\\ref{TEENtab:PillPreg}.  \\citet{Conley1999} standard errors are',
               'reported.',
               paste(sig,'\\end{footnotesize}}',sep=""),
               '\\normalsize\\end{tabular}\\end{table}'),to)
  
  close(to)

}  

if(combine){
    spillA <- readLines(paste(tab.dir,"Spillovers_A.tex", sep=""))
    spillB <- readLines(paste(tab.dir,"Spillovers_B.tex", sep=""))  

    to <- file(paste(tab.dir,"Spillovers.tex", sep=""))
    writeLines(c(spillA,spillB),to)
    close(to)

    spillA <- readLines(paste(tab.dir,"SpilloversROAD_A.tex", sep=""))
    spillB <- readLines(paste(tab.dir,"Spillovers_B.tex", sep=""))  

    to <- file(paste(tab.dir,"SpilloversROAD.tex", sep=""))
    writeLines(c(spillA,spillB),to)
    close(to)

    spillA <- readLines(paste(tab.dir,"SpilloversTIME_A.tex", sep=""))
    spillB <- readLines(paste(tab.dir,"Spillovers_B.tex", sep=""))  

    to <- file(paste(tab.dir,"SpilloversTIME.tex", sep=""))
    writeLines(c(spillA,spillB),to)
    close(to)
}


if(events){
    e1519 <- event(age_sub=15:19, "earlyP")

    postscript(paste(graf.dir,"Event1519_fetaldeath.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    plot(e1519$eventyr,e1519$b, type="b",ylim=c(-10,10),
         col="darkgreen",lwd=2,pch=20, ylab="Estimate",
         xlab="Event Year")
    abline(h = 0, lwd=2, col="gray60")
    points(e1519$eventyr,e1519$b+1.96*e1519$s,type="l",lty=5,pch=20)
    points(e1519$eventyr,e1519$b-1.96*e1519$s,type="l",lty=5,pch=20)
    dev.off()

    e2034 <- event(age_sub=20:34, "earlyP")
    postscript(paste(graf.dir,"Event2034_fetaldeath.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    plot(e2034$eventyr,e2034$b, type="b",ylim=c(-10,10),
         col="darkgreen",lwd=2,pch=20, ylab="Estimate",
         xlab="Event Year")
    abline(h = 0, lwd=2, col="gray60")
    points(e2034$eventyr,e2034$b+1.96*e2034$s,type="l",lty=5,pch=20)
    points(e2034$eventyr,e2034$b-1.96*e2034$s,type="l",lty=5,pch=20)
    dev.off()

    e3549 <- event(age_sub=35:49, "earlyP")
    postscript(paste(graf.dir,"Event3549.eps_fetaldeath",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    plot(e3549$eventyr,e3549$b, type="b",ylim=c(-10,10),
         col="darkgreen",lwd=2,pch=20, ylab="Estimate",
         xlab="Event Year")
    abline(h = 0, lwd=2, col="gray60")
    points(e3549$eventyr,e3549$b+1.96*e3549$s,type="l",lty=5,pch=20)
    points(e3549$eventyr,e3549$b-1.96*e3549$s,type="l",lty=5,pch=20)
    dev.off()
}
