# DeathEstimates-corrected.R       KEL / DCC               yyyy-mm-dd:2013-12-29
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# Corrected estimates of the effect of the morning after pill on fetal deaths.
#
# 
# contact: damian.clarke@usach.cl

#******************************************************************************
#***(1) Parameters
#******************************************************************************
rm(list=ls())
proj.dir <- "YOUR-DIRECTORY-LOCATION-ENDING-IN-SLASH/"

#******************************************************************************
#***(2) Libraries, directories
#******************************************************************************
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if (length(new.pkg))
     	install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

pk <- c("xtable","rms","sandwich","lmtest","plyr","foreign","reshape","gdata",
        "multiwayvcov")
ipak(pk)

deth.dir <- paste(proj.dir, "data/"   , sep="")
tab.dir  <- paste(proj.dir, "tables/" , sep="")
graf.dir <- paste(proj.dir, "figures/", sep="")

Names <- c("dom_comuna","trend","trend2","pill","mujer","party","votes"      ,
           "outofschool","healthspend","healthstaff","healthtraining"        , 
           "educationspend","femalepoverty","urban","year","educationmunicip",
           "condom","usingcont","femaleworkers","poverty","popln")


#******************************************************************************
#***(3) Load Data
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
wild.se <- function(model,cluster) {
    boot <- cluster.boot(model, cluster, boot_type = "wild",
                         wild_type = function() sample(c(-1, 1), 1))
    rcse.se <- coeftest(model, boot)
    return(list(boot, rcse.se))
}

closegen <- function(d1,d2,dat) {
    dat2 <- dat

    dat2$newvar <- NA
    dat2$newvar[dat2$pilldistance > d1 & dat2$pilldistance <= d2 &
                    !(dat2$pilldistance)==0] <- 1
    dat2$newvar[is.na(dat2$newvar)]<-0

    names(dat2)<-c(names(dat),paste('close',d2,sep=""))
    return(dat2)
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
    

    xFl  <- lm(dat$FDb ~ factor(dat$dom_comuna) + factor(dat$year)           +
               factor(dat$pill) + factor(dat$party) + factor(dat$mujer)      +
               dat$votes + dat$outofschool + dat$educationspend              + 
               dat$educationmunicip + dat$healthspend + dat$healthtraining   +
               dat$healthstaff + dat$femalepoverty + dat$femaleworkers       +
               dat$condom, weights=dat$weight)
    clusters <-mapply(paste,"dom_comuna.",dat$dom_comuna,sep="")
    xFl$coefficients2   <- robust.se(xFl,clusters)[[2]]
    
    n  <- nrow(dat)
    s2 <- pillest(xFl, dat, n, "pill", 10)
    mn <- format(round(wtd.mean(dat$FDb,weights=dat$popln),2),nsmall=2) 
    return(list("b" = s2$b, "se" = s2$s, "R2" = s2$r, "n" = s2$n, "m"=mn))
}


#******************************************************************************
#***(5) Estimate
#******************************************************************************
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


#******************************************************************************
#***(6) Export
#******************************************************************************
xvn  <- 'Emergency Contraceptive Pill &'
obs  <- 'Observations&'
R2   <- 'R-Squared'
sig  <- '$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01;'
dpb  <- 'Mean (fetal deaths/live birth)&'
s    <- '\\\\'
a    <- '&'

to <- file(paste(tab.dir,"DeathsPerBirth-corr.tex", sep=""))
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


