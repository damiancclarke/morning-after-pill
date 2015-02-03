# BirthsEstimates.R v1.23          KEL / DCC               yyyy-mm-dd:2013-12-29
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# Import data from S1Data_granular_covars.csv to run various models.  Initial 
# data contains pregnancies by comuna and morning after pill availability by
# comuna.
# 
# Principal model is of the form:
#   preg_{ijt} = alpha + delta*PAE_{jt} + theta_j + i.eta**g(year_t) + u_{ijt}
#
# This code has been written by KEL, with updates by DCC to incorporate additio-
# nal time varying controls and export results to TeX.  When running, the switc-
# hes in section (1) determine whether or not specific sections of the code will
# be run.  Note that the user-contributed stargazer library is very slow, so the
# "full" section is best avoided, or run only once.
#
# NOTE: CBPS and CBMSM should be checked out for PS match
#
# aboe refers to appended back of the envelope calculation which examines
# whether effect sizes look reasonable
#
# last edit v1.23: Clean directory structure.
# contact: damian.clarke@economics.ox.ac.uk

rm(list=ls())

#==============================================================================
#=== (1) Parameters
#==============================================================================
create <- FALSE
preg   <- FALSE
Npreg  <- TRUE
prTab  <- FALSE
spill  <- FALSE
full   <- FALSE
aboe   <- FALSE
ranges <- FALSE
events <- FALSE
ChMund <- FALSE
invPS  <- FALSE
    
birth_y_range <- 2006:2011
pill_y_range <- birth_y_range - 1
age_range <- c(15,49)

#==============================================================================
#=== (2) Libraries, directories
#==============================================================================
require("MatchIt"  )
require("xtable"   )
require("rms"      )
require("plyr"     )
require("glmmML"   )
require("sandwich" )
require("stargazer")
require("lmtest"   )

proj.dir <- "~/universidades/Oxford/DPhil/Thesis/Teens/"

brth.dir <- paste(proj.dir, "Data/Nacimientos/",sep="")
code.dir <- paste(proj.dir, "Source/Births/",sep="")
com.dir  <- paste(proj.dir, "Data/Comunas/", sep="")
ma.dir   <- paste(proj.dir, "Data/PAE/",sep="")
pol.dir  <- paste(proj.dir, "Data/Alcaldes/",sep="")
pop.dir  <- paste(proj.dir, "Data/Population/",sep="")
tab.dir  <- paste(proj.dir, "Tables/", sep="")
graf.dir <- paste(proj.dir, "Figures/", sep="")


Names <- c("dom_comuna","trend","trend2","pill","mujer","party","votes"      ,
           "outofschool","healthspend","healthstaff","healthtraining"        , 
           "educationspend","femalepoverty","year","urban","educationmunicip",
           "condom","usingcont","femaleworkers","region")

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
  
    if (dim==1) {
        null  <- glm(cbind(successes,failures) ~ 1, family="binomial",data=d)
        Lfull <- as.numeric(logLik(outresults))
        Lnull <- as.numeric(logLik(null))
        R2    <- 1 - Lfull/Lnull
    }
    if(dim==10|dim==11) {
        R2 <- summary(outresults)$r.squared
    }
    beta  <- stars(p,beta)
    se    <- paste("(", format(round(se,3),nsmall=3),")", sep="")
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

datcollapse <- function(age_sub,order_sub,ver,dat) {

    dat <- dat[dat$age %in% age_sub,]
    dat <- dat[(dat$order %in% order_sub) | !(dat$pregnant),]
    if(ver==2) {
        dat <- dat[dat$pregnant==1,]
    }
    dat <- closegen(0,15,dat)
    dat <- closegen(15,30,dat)
    dat <- closegen(30,45,dat)
    
    dat$failures <- (1-dat$pregnant)*dat$n
    dat$successes <- dat$pregnant*dat$n

    fmod <- aggregate.data.frame(dat[,c("failures","successes")],
                                 by=list(dat$close15,dat$close30,dat$close45,
                                     dat$dom_comuna,dat$year-2005           ,
                                     (dat$year-2005)^2,dat$pill,dat$mujer   ,
                                     dat$party,dat$votop,dat$outofschool    ,
                                     dat$healthspend,dat$healthstaff        ,
                                     dat$healthtraining,dat$educationspend  ,
                                     dat$femalepoverty,dat$urbBin,dat$year  ,
                                     dat$educationmunicip,dat$condom        ,
                                     dat$usingcont,dat$femaleworkers        ,
                                     dat$region),
                                 function(vec) {sum(na.omit(vec))})
    names(fmod) <- c("close15","close30","close45",Names,"failures","successes")
    fmod$healthstaff      <- fmod$healthstaff/100000
    fmod$healthspend      <- fmod$healthspend/100000
    fmod$healthtraining   <- fmod$healthtraining/100000
    fmod$educationspend   <- fmod$educationspend/100000
    fmod$educationmunicip <- fmod$educationmunicip/100000
    fmod$meanP            <- ave(fmod$pill, group=fmod$dom_comuna)

    return(fmod)
}

#==============================================================================
#=== (4a) Run various models to test effect of PAE on pregnancy 
#==============================================================================
runmod <- function(age_sub,order_sub,num,PSwt) {
  
    formod <- datcollapse(age_sub,order_sub,1,orig)
    if(PSwt) {
        preddat <- datcollapse(age_sub,order_sub,2,orig)
        
        PSc <- glm(pill ~ factor(party) + factor(mujer) + votes + outofschool + 
                   educationspend + educationmunicip + healthspend            + 
                   healthtraining + healthstaff + femalepoverty + condom      +
                   femaleworkers,
                   family=binomial, data=preddat)
        preddat$predict <- predict(PSc, type="response")
        preddat$WT      <- 1/preddat$predict
        #preddat$WT[preddat$pill==1] <- 1/preddat$predict
        #preddat$WT[preddat$pill==0] <- 1/(1-preddat$predict)

        wts    <- preddat[,(names(preddat) %in% c("dom_comuna","WT","trend"))]
        formod <- merge(formod,wts,by=c("dom_comuna","trend"))        
    } else {
        formod$WT <- 1
    }
    

    
    if(num==1) {
        xnT <-  glm(cbind(successes,failures) ~ factor(year) + factor(pill)     +
                    factor(dom_comuna) + factor(region):trend,
                    family="binomial", data=formod, weights=WT)
  
        xCM <-  glm(cbind(successes,failures) ~ factor(year) + factor(pill)     +
                    meanP + factor(region):trend,
                    family="binomial", data=formod, weights=WT)

        xtr <- glm(cbind(successes,failures) ~ factor(year) + factor(pill)      +
                   factor(dom_comuna) + factor(dom_comuna):trend,
                   family="binomial", data=formod, weights=WT)

        xpol <- glm(cbind(successes,failures) ~ factor(year) + factor(pill)     +
                    factor(dom_comuna) + factor(dom_comuna):trend + votes       +
                    factor(party) + factor(mujer)                               ,
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

NumMod <- function(age_sub,order_sub) {
  
    formod <- datcollapse(age_sub, order_sub,2,orig)
    names(formod)[25] <- "births"

    xba <- lm(births ~ factor(dom_comuna) + factor(year) + factor(pill)    +
              factor(dom_comuna):trend, data=formod)
    xct <- lm(births ~ factor(dom_comuna) + factor(dom_comuna):trend       +
              factor(year) + factor(pill) + factor(party) + factor(mujer)  +
              votes + outofschool + educationspend + educationmunicip      +
              healthspend + healthtraining + healthstaff + femalepoverty   + 
              femaleworkers + condom, data=formod  )
    xsp <- lm(births ~ factor(dom_comuna) + factor(dom_comuna):trend       +
              factor(year) + factor(pill) + factor(party) + factor(mujer)  +
              votes + outofschool + educationspend + educationmunicip      +
              healthspend + healthtraining + healthstaff + femalepoverty   +
              condom + femaleworkers + factor(close15) + factor(close30)   +
              + factor(close45),  data=formod)


    clusters <-mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
    xba$coefficients2 <- robust.se(xba,clusters)[[2]]
    xct$coefficients2 <- robust.se(xct,clusters)[[2]]
    xsp$coefficients2 <- robust.se(xsp,clusters)[[2]]
 
 
    n  <- nrow(formod)
    s1 <- pillest(xba, formod, n, "pill", 10)
    s2 <- pillest(xct, formod, n, "pill", 10)
    s3 <- pillest(xsp, formod, n, "pill", 10)
    
    betas <- paste(s1$b, "&", s2$b, "&", s3$b, sep="")
    ses   <- paste(s1$s, "&", s2$s, "&", s3$s, sep="")
    n     <- paste(s1$n, '&', s2$n, '&', s3$n, sep='')
    r     <- paste(s1$r, '&', s2$r, '&', s3$r, sep='')

    return(list("b" = betas,"se" = ses, "n" = n, "r" = r))  
}

#==============================================================================
#=== (4b) Various functions to examine effect of spillover 
#==============================================================================
closegen <- function(d1,d2,dat) {
    dat2 <- dat
    dat2$newvar <- NA  
    dat2$newvar[dat2$pilldistance > d1 & dat2$pilldistance <= d2 &
                !(dat2$pilldistance)==0] <- 1
    dat2$newvar[is.na(dat2$newvar)]<-0
    names(dat2)<-c(names(dat),paste('close',d2,sep=""))
    return(dat2)
    a}

spillovers <- function(age_sub,order_sub) {

    formod <- datcollapse(age_sub, order_sub,1,orig)
    
    xspill <- glm(cbind(successes,failures) ~ factor(dom_comuna)          + 
                  factor(dom_comuna):trend + factor(year) + factor(pill)  + 
                  factor(party) + factor(mujer) + votes + outofschool     + 
                  educationspend + educationmunicip + healthspend         + 
                  healthtraining + healthstaff + femalepoverty + condom   + 
                  femaleworkers + factor(close15) + factor(close30)       +
                  factor(close45), family="binomial",data=formod2)
    clusters <-mapply(paste,"dom_comuna.",formod2$dom_comuna,sep="")
    xspill$coefficients2 <- robust.se(xspill,clusters)[[2]]

    n  <- sum(formod2$successes) + sum(formod2$failures)
    s1 <- pillest(xspill,formod2,n,"pill|close",4)
  
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


rangeest <- function(age_sub,order_sub){

    formod <- datcollapse(age_sub, order_sub,1,orig)
    drops <- c("close15","close30","close45")
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
        n1  <- names(dat)
        dat <- dat[dat$age %in% age_sub,]
        dat <- dat[(dat$order %in% order_sub) | !(dat$pregnant),]  
        dat <- closegen(0,i,dat)
        dat <- closegen(i,i+2.5,dat)
        
        dat$failures  <- (1-dat$pregnant)*dat$n
        dat$successes <- dat$pregnant*dat$n    
        names(dat) <- c(n1,"c1","c2","failures","successes")  
        
        formod <- aggregate.data.frame(dat[,c("failures","successes")],
                                       by=list(dat$c1,dat$c2,dat$dom_comuna,
                                           dat$year-2005,(dat$year-2005)^2 ,
                                           dat$pill,dat$mujer,dat$party    ,
                                           dat$votop,dat$outofschool       ,
                                           dat$healthspend,dat$healthstaff ,
                                           dat$healthtraining              ,
                                           dat$educationspend              ,
                                           dat$femalepoverty,dat$urbBin    ,
                                           dat$year,dat$educationmunicip   ,
                                           dat$condom,dat$usingcont        ,
                                           dat$femaleworkers),
                                       function(vec) {sum(na.omit(vec))})
        
        names(formod) <- c("close1","closemarg",Names)
        
        xrange <- glm(cbind(successes,failures) ~ factor(dom_comuna)          + 
                      factor(dom_comuna):trend + factor(year)  + factor(pill) + 
                      factor(party) + factor(mujer) + votes + outofschool     + 
                      educationspend + educationmunicip + healthspend         + 
                      healthtraining + healthstaff + femalepoverty + condom   + 
                      femaleworkers + factor(close1) + factor(closemarg),
                      family="binomial",data=formod)
        pline  <- grepl("pill",rownames(summary(xrange)$coefficients))
        cline  <- grepl("closemarg",rownames(summary(xrange)$coefficients))  
        pillbeta  <- c(pillbeta,summary(xrange)$coefficients[pline,]["Estimate"])
        pillse    <- c(pillse,summary(xrange)$coefficients[pline,]["Std. Error"])
        closebeta <- c(closebeta,summary(xrange)$coefficients[cline,]["Estimate"])
        closese   <- c(closese,summary(xrange)$coefficients[cline,]["Std. Error"])
        distance  <- c(distance, i)  
    }
    
    return(data.frame(distance,pillbeta,pillse,closebeta,closese))
}

#==============================================================================
#=== (4c) Event study
#==============================================================================
event <- function(age_sub,order_sub,short) {
  
    formod <- datcollapse(age_sub, order_sub,1,orig)
    formod <- formod[with(formod,order(dom_comuna,trend)), ]

    formod$pilltotal <- ave(formod$pill,formod$dom_comuna,FUN=cumsum)
    formod <- formod[with(formod,order(dom_comuna,trend,decreasing=T)), ]
    formod$o <- 0
    formod$o[formod$pilltotal==0]<-1
    formod$add <- ave(formod$o,formod$dom_comuna,FUN=cumsum)
    formod$pilltotal <- formod$pilltotal-formod$add
    
    formod$pilln4 <- formod$pilltotal==-4
    formod$pilln3 <- formod$pilltotal==-3
    formod$pilln2 <- formod$pilltotal==-2
    formod$pilln1 <- formod$pilltotal==-1
    formod$pillp1 <- formod$pilltotal==1
    formod$pillp2 <- formod$pilltotal==2

    eventS  <- glm(cbind(successes,failures) ~ factor(year)                      +
                   factor(dom_comuna) + factor(dom_comuna):trend + votes         +
                   factor(party) + factor(mujer) + outofschool + educationspend  +
                   educationmunicip + healthspend + healthtraining + healthstaff +
                   femalepoverty + femaleworkers + condom + factor(pilln4)       +
                   factor(pilln3) + factor(pilln2) + factor(pilln1)              +
                   factor(pillp1) + factor(pillp2),
                   family="binomial", data=formod)
    clusters <-mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
    eventS$coefficients2 <- robust.se(eventS,clusters)[[2]]


    pillline <- grepl("pill",rownames(summary(eventS)$coefficients))
    beta <- summary(eventS)$coefficients[pillline,][, "Estimate"]
    se   <- summary(eventS)$coefficients[pillline,][, "Std. Error"]

    return(list("b" = beta, "s" = se, "eventyr" = c(-4,-3,-2,-1,0,1)))
}

if(events){
    e1519 <- event(age_sub=15:19, order_sub=1:100)

    plot(e1519$eventyr,e1519$b, type="b",ylim=c(-0.3,0.15),
         col="darkgreen",lwd=2,pch=20, xlab="Estimate",
         ylab="Event Year")
    points(e1519$eventyr,e1519$b+1.96*e1519$s,type="l",lty=3,pch=20)
    points(e1519$eventyr,e1519$b-1.96*e1519$s,type="l",lty=3,pch=20)

    event2034 <- event(age_sub=20:34, order_sub=1)
    plot(event2034$eventyr,event2034$b, type="b",ylim=c(-0.3,0.15),
         col="darkgreen",lwd=2,pch=20, xlab="Estimate",
         ylab="Event Year")
    points(event2034$eventyr,event2034$b+1.96*event2034$s,type="l",lty=3,pch=20)
    points(event2034$eventyr,event2034$b-1.96*event2034$s,type="l",lty=3,pch=20)
}


#==============================================================================
#=== (5) Estimate
#==============================================================================
if(preg|ChMund){
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
  N1519 <- NumMod(age_sub = 15:19, order_sub = 1:100)
  N2034 <- NumMod(age_sub = 20:34, order_sub = 1:100)
  N3549 <- NumMod(age_sub = 35:49, order_sub = 1:100)
  NAll  <- NumMod(age_sub = 15:49, order_sub = 1:100)
}

if(spill){
  c1519 <- spillovers(age_sub = 15:19, order_sub = 1:100)
  c2034 <- spillovers(age_sub = 20:34, order_sub = 1:100)
  c3549 <- spillovers(age_sub = 35:44, order_sub = 1:100)
}  

if(invPS){
    ps1519 <- runmod(age_sub = 15:19, order_sub = 1:100,1,PSwt=TRUE)
    ps2034 <- runmod(age_sub = 20:34, order_sub = 1:100,1,PSwt=TRUE)
    ps3549 <- runmod(age_sub = 35:49, order_sub = 1:100,1,PSwt=TRUE)
}

if(full) {
  full1519 <- runmod(age_sub = 15:19, order_sub = 1:100,2)
  full2034 <- runmod(age_sub = 20:34, order_sub = 1:100,2)
  full3549 <- runmod(age_sub = 35:49, order_sub = 1:100,2)
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
  ra  <- rangeest(age_sub = 15:19, order_sub = 1:100)
  
  postscript(paste(graf.dir,"Dist1519.eps",sep=""),
             horizontal = FALSE, onefile = FALSE, paper = "special",
             height=7, width=9)
  plot(ra$distance,ra$pillbeta, type="b",ylim=c(-0.10,-0.02),
       col="darkgreen",lwd=2,pch=20, xlab="Distance From Treatment Cluster (km)",
       ylab="Esimate of Effect on Treated Cluster")
  points(ra$distance,ra$pillbeta+1.96*ra$pillse,type="l",lty=3,pch=20)
  points(ra$distance,ra$pillbeta-1.96*ra$pillse,type="l",lty=3,pch=20)

  legend("topright",legend=c("Point Estimate","95% CI"),
         text.col=c("darkgreen","black"),pch=c(20,NA),lty=c(1,3),
         col=c("darkgreen","black"))
  dev.off()


  ra  <- rangeest(age_sub = 20:34, order_sub = 1:100)
  
  postscript(paste(graf.dir,"Dist2034.eps",sep=""),
             horizontal = FALSE, onefile = FALSE, paper = "special",
             height=7, width=9)
  plot(ra$distance,ra$pillbeta, type="b",ylim=c(-0.06,-0.02),
       col="darkgreen",lwd=2,pch=20, xlab="Distance From Treatment Cluster (km)",
       ylab="Esimate of Effect on Treated Cluster")
  points(ra$distance,ra$pillbeta+1.96*ra$pillse,type="l",lty=3,pch=20)
  points(ra$distance,ra$pillbeta-1.96*ra$pillse,type="l",lty=3,pch=20)
  
  legend("topright",legend=c("Point Estimate","95% CI"),
         text.col=c("darkgreen","black"),pch=c(20,NA),lty=c(1,3),
         col=c("darkgreen","black"))
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

if(prTab){
to <-file(paste(tab.dir,"Births.tex", sep=""))
writeLines(c('\\begin{landscape}','\\begin{table}[!htbp] \\centering',
           '\\caption{The Effect of the Morning After Pill on Pregnancy}',
           '\\label{TEENtab:PillPreg}',
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
           'All Births and First Births are binary variables taking the value',
           'of 1 in the case that a women gives live birth and that this',
           'occurs at any birth order, or is her first birth (respectively).',
           'All models are estimated using logistic regression and include',
           'comuna and year fixed. Standard errors are clustered at the level',
           'of the comuna.  All coefficients are',
           'reported as log odds and in each case Pill is a binary variable',
           'referring to the availability of the morning after pill in the',
           'woman\'s comuna and (lagged) year.  Political controls include',
           'party dummies for the mayor in power, the mayor\'s gender, and',
           'the vote margin of the mayor.  Health and education controls',
           'include the percent of girls out of highschool, education',
           'spending by both the municipality and the Ministry of Education',
           'and total health spending and health spending on staff and',
           'training.  Gender controls are the percent of female heads of',
           ' households living below the poverty line, and the percent of',
           'female workers in professional positions in the Municipality.',
           paste(sig, '\\end{footnotesize}}', sep=""),
           '\\normalsize\\end{tabular}\\end{table}\\end{landscape}'),to)
close(to)
}

if(spill){
  to <-file(paste(tab.dir,"Spillovers_A.tex", sep=""))
  writeLines(c('\\begin{table}[!htbp] \\centering',
             '\\caption{The Morning After Pill and Treatment Spillovers}',
             '\\label{TEENtab:Spillover} \\begin{tabular}',
             '{@{\\extracolsep{5pt}}lccc}\\\\[-1.8ex]\\hline\\hline\\\\',
             '[-1.8ex] & 15-19 & 20-34 & 35-49 \\\\',
             '& Year olds & Year olds & Year olds \\\\ \\midrule',
             '\\multicolumn{4}{l}{\\textsc{\\noindent Panel A: Births}} \\\\',
             '& & & \\\\',
             paste(xvar,c1519$b[1],'&',c2034$b[1],'&',c3549$b[1],'\\\\',sep=""),
             paste('&',c1519$s[1],'&',c2034$s[1],'&',c3549$s[1],'\\\\',sep=""),            
             paste(xv2,c1519$b[2],'&',c2034$b[2],'&',c3549$b[2],'\\\\',sep=""),
             paste('&',c1519$s[2],'&',c2034$s[2],'&',c3549$s[2],'\\\\',sep=""),
             paste(xv3,c1519$b[3],'&',c2034$b[3],'&',c3549$b[3],'\\\\',sep=""),
             paste('&',c1519$s[3],'&',c2034$s[3],'&',c3549$s[3],'\\\\',sep=""), 
             paste(xv4,c1519$b[4],'&',c2034$b[4],'&',c3549$b[4],'\\\\',sep=""),
             paste('&',c1519$s[4],'&',c2034$s[4],'&',c3549$s[4],'\\\\',sep=""), 
             '& & & \\\\',
             paste(obs,c1519$n,'&',c2034$n,'&',c3549$n,'\\\\',sep=""),
             paste(R2,c1519$r,'&',c2034$r,'&',c3549$r,'\\\\ \\midrule',sep="")),
             to)
close(to)
}

if(full) {
  stargazer(full1519, full2034,  full3549,
          title="The Morning After Pill and Pregnancy: Full Covariates",
          align=TRUE, label="TEENtabPregFull",omit.stat=c("LL","ser","f"),
          keep=c("pill","mujer","votes","outofschool","educationspend",
                 "educationmunicip","healthspend","healthtraining",
                 "healthstaff","femalepoverty","femaleworkers"), 
          column.labels=c("15-19 year olds","20-34 year olds","35-49 year olds"),
          column.separate=(c(1,1,1)),
          out=paste(tab.dir, "PregFullCovars2.tex", sep=""),
          dep.var.labels="Pregnancy",
          covariate.labels=c("Morning After Pill","Female Mayor","Mayor's Support",
                             "Out of School","Total Education Spending", 
                             "Municipal Education Spending", "Health Spending",
                             "Health Training", "Health Staff", "Female Poverty",
                             "Female Workers"),
          notes="\\begin{footnotesize} \\textsc{Notes:} Each model is identical to 
            column (4) of table \\ref{TEENtab:PillPreg}.  A description of each 
            variable is also provided in table \\ref{TEENtab:PillPreg}.  Municipality
            dummies and trends and political party dummies have been omitted for 
            clarity. ^{*}p$<$0.1; ^{**}p$<$0.05; ^{***}p$<$0.01 
            \\end{footnotesize}",
            notes.align="l", notes.append=FALSE, 
            table.placement="htpb!")

  ##Read in file and grep out certain undesired features
  replaceT <- function(filename,find,replace) {
    worfile  <- file(filename)
    editfile <- readLines(workfile)
    newfile  <- gsub(find,replace,editfile,fixed=T)
    writeLines(newfile,workfile) 
    close(workfile)  
  }
  fix <-paste(tab.dir,"PregFullCovars2.tex", sep="")
  replaceT(fix,"lD{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3}",'lccc')
  replaceT(fix,"{\\textit{Dependent variable:}}",'{Pregnancy}')
  replaceT(fix,"15-19 year olds","15-19")
  replaceT(fix,"20-34 year olds","20-34")
  replaceT(fix,"35-49 year olds","35-49")
  replaceT(fix,"35-49} \\\\","35-49} \\\\ & year olds & year olds & year olds \\\\ ")  
  replaceT(fix,"^{*}","$^{*}$")
  replaceT(fix,"^{**}","$^{**}$")
  replaceT(fix,"^{***}","$^{***}$")
  replaceT(fix,"Observations","Years $\\times$ Municipality")
  replaceT(fix,"\\textit{Note:}  & \\multicolumn{3}{l","\\multicolumn{4}{p{10.8cm}")
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
            '\\multicolumn{3}{p{6cm}}{\\begin{footnotesize}\\textsc{Notes:} ',
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
    to <-file(paste(tab.dir,"aggregateBirths.tex", sep=""))
    writeLines(c('\\begin{table}[!htbp] \\centering',
                 '\\caption{Estimates Based on Aggregate Comunal Data}',
                 '\\label{TEENtab:aggregate}',
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '& N Births & N Births & N Births \\\\',
                 '&(1)&(2)&(3) \\\\ \\hline',
                 ' & & & \\\\',
                 paste('Pill (15-19 years) &',N1519$b,'\\\\',sep=""),
                 paste('                   &',N1519$s,'\\\\',sep=""),
                 ' & & & \\\\',
                 paste('Observations       &',N1519$n,'\\\\',sep=""),
                 paste('R-squared          &',N1519$r,'\\\\',sep=""),
                 ' & & & \\\\',
                 paste('Pill (20-34 years) &',N2034$b,'\\\\',sep=""),
                 paste('                   &',N2034$s,'\\\\',sep=""),
                 ' & & & \\\\',
                 paste('Observations       &',N2034$n,'\\\\',sep=""),
                 paste('R-squared          &',N2034$r,'\\\\',sep=""),
                 ' & & & \\\\',
                 paste('Pill (35-49 years) &',N3549$b,'\\\\',sep=""),
                 paste('                   &',N3549$s,'\\\\',sep=""),
                 ' & & & \\\\',
                 paste('Observations       &',N3549$n,'\\\\',sep=""),
                 paste('R-squared          &',N3549$r,'\\\\',sep=""),
                 '\\hline \\\\[-1.8ex] ', 
                 '{\\small Year \\& Comuna FEs} & Y & Y & Y \\\\',
                 '{\\small Trend, Controls} & & Y & Y \\\\', 
                 '{\\small Spillovers} & & & Y \\\\',
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{4}{p{9.2cm}}{\\begin{footnotesize}\\textsc{Notes:}',
                 'Each panel presents difference-in-difference results for a',
                 'regression of the total count of pregnancies for the age',
                 'group in each municipality.  All models are estimated by OLS',
                 'and standard errors are clustered at the level of the comuna.',
                 'Controls are described in table \\ref{TEENtab:PillPreg}.',
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
                 '\\begin{tabular}{@{\\extracolsep{5pt}}lcccc}',
                 '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
                 '&Pr(Birth)&Pr(Birth)&Pr(Birth)&Pr(Birth)\\\\',
                 '&(1)&(2)&(3)&(4)\\\\ \\hline',
                 '\\multicolumn{5}{l}{\\textsc{\\noindent 15-19 year olds}} \\\\',
                 ' & & & & \\\\',
                 paste(xvar,ps1519$b,'\\\\', sep=""), 
                 paste(' &',ps1519$se, '\\\\', sep=""), 
                 ' & & & & \\\\',
                 '\\multicolumn{5}{l}{\\textsc{\\noindent 20-34 year olds}} \\\\',
                 ' & & & & \\\\', 
                 paste(xvar,ps2034$b,'\\\\', sep=""), 
                 paste(' &',ps2034$se,'\\\\', sep=""), 
                 ' & & & &\\\\',
                 '\\multicolumn{5}{l}{\\textsc{\\noindent 35-49 year olds}} \\\\',
                 ' & & & & \\\\', 
                 paste(xvar,ps3549$b,'\\\\', sep=""), 
                 paste(' &',ps3549$se,'\\\\', sep=""), 
                 ' & & & & \\\\',
                 '\\hline \\\\[-1.8ex] ', 
                 '{\\small Trends \\& FEs} & Y & Y & Y & Y \\\\',
                 '{\\small Political Controls} & & Y & Y & Y \\\\', 
                 '{\\small Health, Educ, Gender Controls} & & & Y & Y \\\\',
                 '{\\small Condom Availability} & & & & Y \\\\', 
                 '\\hline \\hline \\\\[-1.8ex]',
                 '\\multicolumn{5}{p{12.6cm}}{\\begin{footnotesize}\\textsc{Notes:}',
                 'Regression results estimated using inverse propensity score ',
                 'weighting based on Pr(treatment) on full observables.',
                 'Specifications and controls identical to those described in table',
                 '\\ref{TEENtab:PillPreg}.',
                 paste(sig, '\\end{footnotesize}}', sep=""),
                 '\\normalsize\\end{tabular}\\end{table}'),to)
    close(to)
}
