# BirthsEstimates-corrected.R      KEL / DCC               yyyy-mm-dd:2013-12-29
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# Corrected estimates of the effect of the morning after pill on fertility. This
# file contains the difference-in-difference estimates and the event study spec-
# ifications.
#
#
# contact: damian.clarke@usach.cl

#==============================================================================
#=== (1) Parameters
#==============================================================================
rm(list=ls())
proj.dir <- "YOUR-DIRECTORY-LOCATION-ENDING-IN-SLASH/"

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
	      "plotrix","foreign","reshape","gdata","Hmisc","multiwayvcov")
ipak(packages)


brth.dir <- paste(proj.dir, "data/"   ,sep="")
tab.dir  <- paste(proj.dir, "tables/" ,sep="")
graf.dir <- paste(proj.dir, "figures/",sep="")


Names <- c("dom_comuna","trend","trend2","pill","mujer","party","votes"      ,
           "outofschool","healthspend","healthstaff","healthtraining"        , 
           "educationspend","femalepoverty","urban","year","educationmunicip",
           "condom","usingcont","femaleworkers","region","popln")

#==============================================================================
#=== (3) Load Data
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
wild.se <- function(model,cluster) {
    set.seed(2727)
    boot <- cluster.boot(model, cluster, boot_type = "wild",
                         wild_type = function() sample(c(-1, 1), 1))
    rcse.se <- coeftest(model, boot)
    return(list(boot, rcse.se))
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
    dat2$newvar3[dat2$travelTime/60>d1*1.5 & dat2$travelTime/60<=d2*1.5 &
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

    dat <- closegen(0,10,dat)
    dat <- closegen(10,20,dat)
    dat <- closegen(20,30,dat)

    
    dat$failures <- (1-dat$pregnant)*dat$n
    dat$successes <- dat$pregnant*dat$n

    if(ver==2) {
        dat$n <- 1
        full  <- aggregate.data.frame(dat$n,
                                      by=list(dat$close10,dat$close20      ,
                                          dat$close30,dat$road10,dat$road20,
                                          dat$road30,dat$time10,dat$time20 ,
                                          dat$time30,dat$dom_comuna        ,
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
        names(full) <- c("close10","close20","close30","road10","road20","road30",
                         "time10","time20","time30",Names,"n")
        dat   <- dat[dat$pregnant==1,]
    }

    fmod <- aggregate.data.frame(dat[,c("failures","successes")],
                                 by=list(dat$close10,dat$close20,dat$close30,
                                     dat$road10,dat$road20,dat$road30       ,
                                     dat$time10,dat$time20,dat$time30       ,
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
    names(fmod) <- c("close10","close20","close30","road10","road20","road30",
                     "time10","time20","time30",Names,"failures","successes")
    
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

#==============================================================================
#=== (4a) Run various models to test effect of PAE on pregnancy 
#==============================================================================
NumMod <- function(age_sub,order_sub,rate,logm,wt,PSwt) {
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
    if(PSwt) {
        wts <- inversePS(age_sub,order_sub,orig)
        formod <- join(formod,wts,by="dom_comuna",type="left",match="all")

        xPS <- lm(births ~ factor(year) + factor(pill) +
                   factor(dom_comuna), data=formod, weights=WT)

        clusters <- mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
        clusters <- formod$dom_comuna
        xPS$coefficients2  <- robust.se(xPS,clusters)[[2]]
        n  <- nrow(formod)
        
        s0 <- pillest(xPS, formod, n, "pill", 10)
        return(s0)        
    }
    f1 <- formod
    
    xba <- lm(f1$births ~ factor(f1$dom_comuna) + factor(f1$year) + f1$pill,
              weights=f1$weight)
    xtr <- lm(f1$births ~ factor(f1$dom_comuna) + factor(f1$year) +
              factor(f1$pill) + factor(f1$dom_comuna):f1$trend,
              weights=f1$weight)

    xct <- lm(f1$births ~ factor(f1$dom_comuna) + factor(f1$year) +
              factor(f1$pill) + 
              factor(f1$party) + factor(f1$mujer) + f1$votes + f1$outofschool
              + f1$educationspend + f1$educationmunicip + f1$healthspend +
              f1$healthtraining + f1$healthstaff + f1$femalepoverty   + 
              f1$femaleworkers + f1$condom, weights=f1$weight)

    xsp <- lm(f1$births ~ factor(f1$dom_comuna) + factor(f1$year) +
              factor(f1$pill) + 
              factor(f1$party) + factor(f1$mujer) + f1$votes + f1$outofschool
              + f1$educationspend + f1$educationmunicip + f1$healthspend +
              f1$healthtraining + f1$healthstaff + f1$femalepoverty +
              f1$condom + f1$femaleworkers + factor(f1$close10) +
              factor(f1$close20) + factor(f1$close30), weights=f1$weight)
    
    clusters <- formod$dom_comuna
    print("Bootstrap 1") 
    xba$coefficients2 <- wild.se(xba,clusters)[[2]]
    print("Bootstrap 2") 
    xtr$coefficients2 <- wild.se(xtr,clusters)[[2]]
    print("Bootstrap 3") 
    xct$coefficients2 <- wild.se(xct,clusters)[[2]]
    print("Bootstrap 4") 
    xsp$coefficients2 <- wild.se(xsp,clusters)[[2]]
  
 
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
#=== (4b) Event study
#==============================================================================
event <- function(age_sub,order_sub,logm) {
  
    formod <- datcollapse(age_sub,order_sub,2,orig)
    formod <- formod[with(formod,order(dom_comuna,trend)), ]
    names(formod)[32] <- "births"

    if(logm) {
        formod$births<-log(formod$births+1)
    } else {
        formod$births <- (formod$births/formod$popln)*1000
    }
    

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
    f1 <- formod
    
    eventS <- lm(f1$births ~ factor(f1$year) + factor(f1$dom_comuna) +
                 factor(f1$dom_comuna):f1$trend + factor(f1$pilln5)  +
                 factor(f1$pilln4) + factor(f1$pilln2)               +
                 factor(f1$pilln1)  + factor(f1$pillp0)              +
                 factor(f1$pillp1) + factor(f1$pillp2),
                 weights=f1$popln)
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
r1519 <- NumMod(15:19, 1:100, T, F, F, F)
r2034 <- NumMod(20:34, 1:100, T, F, F, F)
r3549 <- NumMod(35:49, 1:100, T, F, F, F)
rAll  <- NumMod(15:49, 1:100, T, F, F, F)


eventPlot <- function(res,name,min,max) {
    postscript(paste(graf.dir,name,sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    
    plotCI(res$eventyr,res$b,ui=res$b+1.96*res$s,li=res$b-1.96*res$s,
           ylim=c(min,max), ylab="Estimate", xlab="Event Year",
           cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)
    
    points(res$eventyr,res$b,type="l",lwd=2,pch=20)
    abline(h =  0  , lwd=1, col="gray60", lty = 2)
    abline(v = -0.1, lwd=2, col="red")
    dev.off()
}
e1519 <- event(15:19,1,logm=FALSE)
e2034 <- event(20:34,1,logm=FALSE)
eAll  <- event(15:49,1,logm=FALSE)
eventPlot(e1519,"Event1519-corr.eps",-10,5)
eventPlot(e2034,"Event2034-corr.eps",-10,5)
eventPlot(eAll ,"EventAll-corr.eps" ,-5,2.5)


#==============================================================================
#=== (6) Export results
#==============================================================================
 xvar <- 'Emergency Contraceptive Pill &'
 obs  <- 'Observations&'
 R2   <- 'R-Squared&'
 sig  <- '$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01'
 
 to <-file(paste(tab.dir,"aggregateASFR-corr.tex", sep=""))
 writeLines(c('\\begin{table}[!htbp] \\centering',
              '\\caption{The Effect of the EC Pill on Birth Rates [Updated]}',
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
              '{\\small Municipal-Specific Linear Trends}& &Y& &  \\\\', 
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
              'are estimated by OLS, and each municipality is          ',
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
