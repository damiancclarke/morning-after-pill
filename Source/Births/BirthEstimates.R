# BirthsEstimates.R v1.22          KEL / DCC               yyyy-mm-dd:2013-12-29
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# Import data from S1Data_granular_covars.csv to run various models.  Initial 
# data contains pregnancies by comuna and morning after pill availability by
# comuna.
# 
# Principal model is of the form:
#   preg_{ijt} = alpha + delta*PAE_{jt} + theta_j + i.eta**g(year_t) + u_{ijt}
#
# This code has been written by KEL, with updates by DCC to incorporate 
# additional time varying controls and export results to TeX.  When running the
# switches in section (1) determine whether or not specific sections of the
# code will be run.
#
# aboe refers to appended back of the envelope calculation which examines
# whether effect sizes look reasonable
#
# last edit v1.22: Refactorise
# contact: damian.clarke@economics.ox.ac.uk

rm(list=ls())

#==============================================================================
#=== (1) Parameters
#==============================================================================
create <- FALSE
preg   <- TRUE
spill  <- TRUE
full   <- FALSE
aboe   <- FALSE
ranges <- FALSE

birth_y_range <- 2006:2011
pill_y_range <- birth_y_range - 1
age_range <- c(15,49)

#==============================================================================
#=== (2) Libraries, directories
#==============================================================================
require("xtable")
require("rms")
require("plyr")
require("glmmML")
require("sandwich")
require("lmtest")
require("stargazer")

proj.dir <- "~/universidades/Oxford/DPhil/Thesis/Teens/"
geo.dir  <- "~/database/ChileRegiones/Nombres/"
pol.dir  <- paste(proj.dir, "Data/Alcaldes/",sep="")
ma.dir   <- paste(proj.dir, "Data/PAE/",sep="")
brth.dir <- paste(proj.dir, "Data/MinSal/dta/",sep="")
pop.dir  <- paste(proj.dir, "Data/Poblacion/proyecciones/DatCom/",sep="")
com.dir  <- paste(proj.dir, "Data/Comunas/", sep="")
work.dir <- paste(proj.dir, "Data/Nacimientos/",sep="")
code.dir <- paste(proj.dir, "Source/Births/",sep="")
tab.dir  <- paste(proj.dir, "Tables/", sep="")
graf.dir <- paste(proj.dir, "Figures/", sep="")


Names <- c("dom_comuna","trend","trend2","pill","mujer","party","votes"      ,
           "outofschool","healthspend","healthstaff","healthtraining"        , 
           "educationspend","femalepoverty","year","urban","educationmunicip",
           "condom","usingcont","femaleworkers","failures","successes")

#==============================================================================
#=== (3a) Source functions
#==============================================================================
if(create){
  f <- paste(code.dir,"BirthGenerate.R",sep="")
  source(f)

  filename <- paste(work.dir, 'S1Data_granular_covars2.csv' ,sep="")
  prep_s1_data(age_range,usecom="TRUE",filename)
}

#==============================================================================
#=== (3b) Load Data
#==============================================================================
f <- paste(work.dir, "S1Data_granular_covars2.csv", sep="")
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
  
  if(dim==1) {
    beta <- summary(outresults)$coefficients[pillline,]["Estimate"]
    se   <- outresults$coefficients2[pillline,]["Std. Error"]
    p    <- outresults$coefficients2[pillline,]["Pr(>|z|)"]    
  }
  else {
    beta <- summary(outresults)$coefficients[pillline,][, "Estimate"]
    se   <- outresults$coefficients2[pillline,][, "Std. Error"]
    p    <- outresults$coefficients2[pillline,][, "Pr(>|z|)"]    
  }
  
  null  <- glm(cbind(successes,failures) ~ 1, family="binomial",data=d)
  Lfull <- as.numeric(logLik(outresults))
  Lnull <- as.numeric(logLik(null))
  R2    <- 1 - Lfull/Lnull
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

#==============================================================================
#=== (4a) Run various models to test effect of PAE on pregnancy 
#==============================================================================
runmod <- function(age_sub,order_sub,num) {
  
  dat <- orig
  dat <- dat[dat$age %in% age_sub,]
  dat <- dat[(dat$order %in% order_sub) | !(dat$pregnant),]

  dat$failures <- (1-dat$pregnant)*dat$n
  dat$successes <- dat$pregnant*dat$n

  formod <- aggregate.data.frame(dat[,c("failures","successes")],
                                 by=list(dat$dom_comuna,dat$year-2005         ,
                                         (dat$year-2005)^2,dat$pill,dat$mujer ,
                                         dat$party,dat$votop,dat$outofschool  ,
                                         dat$healthspend,dat$healthstaff      ,
                                         dat$healthtraining,dat$educationspend,
                                         dat$femalepoverty,dat$urbBin,dat$year,
                                         dat$educationmunicip,dat$condom      ,
                                         dat$usingcont,dat$femaleworkers),
                                         function(vec) {sum(na.omit(vec))})
  names(formod)           <- Names
  formod$healthstaff      <- formod$healthstaff/100000
  formod$healthspend      <- formod$healthspend/100000
  formod$healthtraining   <- formod$healthtraining/100000
  formod$educationspend   <- formod$educationspend/100000
  formod$educationmunicip <- formod$educationmunicip/100000
  
  if(num==1) {
    xtr   <- glm(cbind(successes,failures) ~ factor(dom_comuna)        +
                 factor(dom_comuna):trend + factor(year) + factor(pill), 
                 family="binomial", data=formod)
    clusters <-mapply(paste,"dom_comuna.",formod$dom_comuna,sep="")
    xtr$coefficients2 <- robust.se(xtr,clusters)[[2]]
  
    xpol  <- glm(cbind(successes,failures) ~ factor(dom_comuna)      + 
             factor(dom_comuna):trend + factor(year) + factor(pill)  + 
             factor(party) + factor(mujer) + votes, family="binomial",
             data=formod)
    xpol$coefficients2 <- robust.se(xpol,clusters)[[2]]
  
    xsh   <- glm(cbind(successes,failures) ~ factor(dom_comuna)           + 
                  factor(dom_comuna):trend + factor(year) + factor(pill)  + 
                  factor(party) + factor(mujer) + votes + outofschool     +
                  educationspend + educationmunicip + healthspend         + 
                  healthtraining + healthstaff, family="binomial"         , 
                data=formod)
    xsh$coefficients2 <- robust.se(xsh,clusters)[[2]]
  
    xfem  <- glm(cbind(successes,failures) ~ factor(dom_comuna)          + 
                 factor(dom_comuna):trend + factor(year) + factor(pill)  + 
                 factor(party) + factor(mujer) + votes + outofschool     + 
                 educationspend + educationmunicip + healthspend         + 
                 healthtraining + healthstaff + femalepoverty            + 
                 femaleworkers, family="binomial", data=formod)
    xfem$coefficients2 <- robust.se(xfem,clusters)[[2]]
  }
  xcont  <- glm(cbind(successes,failures) ~ factor(dom_comuna)           + 
                 factor(dom_comuna):trend + factor(year) + factor(pill)  + 
                 factor(party) + factor(mujer) + votes + outofschool     + 
                 educationspend + educationmunicip + healthspend         + 
                 healthtraining + healthstaff + femalepoverty            + 
                 femaleworkers + condom, family="binomial", data=formod  )
  xcont$coefficients2 <- robust.se(xcont,clusters)[[2]]
  
  if(num==1) {
    n  <- sum(formod$successes) + sum(formod$failures)

    s1 <- pillest(xtr,   formod, n, "pill", 1)
    s2 <- pillest(xpol,  formod, n, "pill", 1)
    s3 <- pillest(xsh,   formod, n, "pill", 1)
    s4 <- pillest(xfem,  formod, n, "pill", 1)
    s5 <- pillest(xcont, formod, n, "pill", 1)
    
    betas <- paste(s1$b, "&", s2$b, "&", s4$b, "&", s5$b, sep="")
    ses   <- paste(s1$s, "&", s2$s, "&", s4$s, "&", s5$s, sep="")
    n     <- paste(s1$n, '&', s2$n, '&', s4$n, '&' ,s5$n, sep='')
    r     <- paste(s1$r, '&', s2$r, '&', s4$r, '&', s5$r, sep='')

    return(list("b" = betas,"se" = ses, "n" = n, "r" = r))  
  } else {
    return(xcont)
  }
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
}

spillovers <- function(age_sub,order_sub) {
  dat <- orig
  dat <- dat[dat$age %in% age_sub,]
  dat <- dat[(dat$order %in% order_sub) | !(dat$pregnant),]  
  
  dat <- closegen(0,15,dat)
  dat <- closegen(15,30,dat)
  dat <- closegen(30,45,dat)
  
  dat$failures  <- (1-dat$pregnant)*dat$n
  dat$successes <- dat$pregnant*dat$n
  
  formod2 <- aggregate.data.frame(dat[,c("failures","successes")],
                                  by=list(dat$close15,dat$close30,dat$close45  ,
                                          dat$dom_comuna,dat$year-2005         ,
                                          (dat$year-2005)^2,dat$pill,dat$mujer ,
                                          dat$party,dat$votop,dat$outofschool  ,
                                          dat$healthspend,dat$healthstaff      ,
                                          dat$healthtraining,dat$educationspend,
                                          dat$femalepoverty,dat$urbBin,dat$year,
                                          dat$educationmunicip,dat$condom      ,
                                          dat$usingcont,dat$femaleworkers),
                                  function(vec) {sum(na.omit(vec))})
  
  names(formod2) <- c("close15","close30","close45",Names)
    
  xspill <- glm(cbind(successes,failures) ~ factor(dom_comuna)            + 
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

  dat <- orig
  dat <- dat[dat$age %in% age_sub,]
  dat <- dat[(dat$order %in% order_sub) | !(dat$pregnant),] 
  dat$failures  <- (1-dat$pregnant)*dat$n
  dat$successes <- dat$pregnant*dat$n
  formod <- aggregate.data.frame(dat[,c("failures","successes")],
                                 by=list(dat$dom_comuna,dat$year-2005         ,
                                         (dat$year-2005)^2,dat$pill,dat$mujer ,
                                         dat$party,dat$votop,dat$outofschool  ,
                                         dat$healthspend,dat$healthstaff      ,
                                         dat$healthtraining,dat$educationspend,
                                         dat$femalepoverty,dat$urbBin,dat$year,
                                         dat$educationmunicip,dat$condom      ,
                                         dat$usingcont,dat$femaleworkers),
                                 function(vec) {sum(na.omit(vec))})
  names(formod)           <- Names
    
  xrange <- glm(cbind(successes,failures) ~ factor(dom_comuna)            + 
                  factor(dom_comuna):trend + factor(year) + factor(pill)  +
                  factor(party) + factor(mujer) + votes + outofschool     + 
                  educationspend + educationmunicip + healthspend         + 
                  healthtraining + healthstaff + femalepoverty + condom   +
                  femaleworkers, family="binomial",data=formod)
  pillline <- grepl("pill",rownames(summary(xrange)$coefficients))
  closeline <- grepl("closemarg",rownames(summary(xrange)$coefficients))  
  pillbeta <- summary(xrange)$coefficients[pillline,]["Estimate"]
  pillse <- summary(xrange)$coefficients[pillline,]["Std. Error"]
  closebeta <- summary(xrange)$coefficients[closeline,]["Estimate"]
  closese <- summary(xrange)$coefficients[closeline,]["Std. Error"]
  distance <- 0
  
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
                                    by=list(dat$c1,dat$c2,dat$dom_comuna   ,
                                            dat$year-2005,(dat$year-2005)^2,
                                            dat$pill,dat$mujer,dat$party   ,
                                            dat$votop,dat$outofschool      ,
                                            dat$healthspend,dat$healthstaff,
                                            dat$healthtraining             ,
                                            dat$educationspend,
                                            dat$femalepoverty,dat$urbBin   ,
                                            dat$year,dat$educationmunicip  ,
                                            dat$condom,dat$usingcont       ,
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
    pillline <- grepl("pill",rownames(summary(xrange)$coefficients))
    closeline <- grepl("closemarg",rownames(summary(xrange)$coefficients))  
    pillbeta <- c(pillbeta,summary(xrange)$coefficients[pillline,]["Estimate"])
    pillse <- c(pillse,summary(xrange)$coefficients[pillline,]["Std. Error"])
    closebeta <- c(closebeta,summary(xrange)$coefficients[closeline,]["Estimate"])
    closese <- c(closese,summary(xrange)$coefficients[closeline,]["Std. Error"])
    distance <- c(distance, i)  
  }

  return(data.frame(distance,pillbeta,pillse,closebeta,closese))
}

#==============================================================================
#=== (5) Estimate
#==============================================================================
if(preg){
  a1519 <- runmod(age_sub = 15:19, order_sub = 1:100,1)
  a2034 <- runmod(age_sub = 20:34, order_sub = 1:100,1)
  a3549 <- runmod(age_sub = 35:49, order_sub = 1:100,1)
  b1519 <- runmod(age_sub = 15:19, order_sub = 1,1)
  b2034 <- runmod(age_sub = 20:34, order_sub = 1,1)
  b3549 <- runmod(age_sub = 35:49, order_sub = 1,1)
}


if(spill){
  c1519 <- spillovers(age_sub = 15:19, order_sub = 1:100)
  c2034 <- spillovers(age_sub = 20:34, order_sub = 1:100)
  c3549 <- spillovers(age_sub = 35:44, order_sub = 1:100)
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

if(preg){
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
