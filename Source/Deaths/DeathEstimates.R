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
create  <- TRUE
death   <- TRUE
spill   <- TRUE
combine <- TRUE
full    <- FALSE

birth_y_range <- 2006:2011
pill_y_range  <- birth_y_range - 1
age_range     <- c(15,49)
week          <- 20
pat           <- "P"

#******************************************************************************
#***(2) Libraries, directories
#******************************************************************************
require(xtable)
require(rms)
require(sandwich)
require(lmtest)
require(stargazer)

proj.dir <- "~/universidades/Oxford/DPhil/Thesis/Teens/"

brth.dir <- paste(proj.dir, "Data/Nacimientos/",sep="")
code.dir <- paste(proj.dir, "Source/Deaths/",sep="")
com.dir  <- paste(proj.dir, "Data/Comunas/", sep="")
deth.dir <- paste(proj.dir, "Data/Deaths/",sep="")
ma.dir   <- paste(proj.dir, "Data/PAE/",sep="")
pol.dir  <- paste(proj.dir, "Data/Alcaldes/",sep="")
pop.dir  <- paste(proj.dir, "Data/Poblacion/proyecciones/DatCom/",sep="")
tab.dir  <- paste(proj.dir, "Tables/", sep="")

Names <- c("dom_comuna","trend","trend2","pill","mujer","party","votes"      ,
           "outofschool","healthspend","healthstaff","healthtraining"        , 
           "educationspend","femalepoverty","year","urban","educationmunicip",
           "condom","usingcont","femaleworkers")


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

#******************************************************************************
#***(4b) Main Functions
#******************************************************************************
death_pmod <- function(age_sub,deathtype,numb) {
 
  dat       <- orig
  dat       <- dat[dat$age %in% age_sub,]
  dat       <- dat[dat$pregnant == 1,]
  dat$Q     <- dat$earlyQ+dat$lateQ
  dat$early <- dat$earlyQ+dat$earlyP
  dat$late  <- dat$lateQ+dat$lateP
  
  dat <- aggregate.data.frame(dat[,c("n",deathtype)],
                              by=list(dat$dom_comuna,dat$year-2005         ,
                                      (dat$year-2005)^2,dat$pill,dat$mujer ,
                                      dat$party,dat$votop,dat$outofschool  ,
                                      dat$healthspend,dat$healthstaff      ,
                                      dat$healthtraining,dat$educationspend,
                                      dat$femalepoverty,dat$urbBin,dat$year,
                                      dat$educationmunicip,dat$condom      ,
                                      dat$usingcont,dat$femaleworkers),
                                 function(vec) {sum(na.omit(vec))})
  names(dat)              <- c(Names,"n","death")

  dat$healthstaff      <- dat$healthstaff/100000
  dat$healthspend      <- dat$healthspend/100000
  dat$healthtraining   <- dat$healthtraining/100000
  dat$educationspend   <- dat$educationspend/100000
  dat$educationmunicip <- dat$educationmunicip/100000
  

  dat$failures  <- dat$n 
  dat$successes <- dat$death

  formod <- aggregate.data.frame(dat[,c("failures","successes")],
                                 by=list(dat$dom_comuna,dat$trend,dat$trend2    ,
                                         dat$pill,dat$mujer,dat$party,dat$votes ,
                                         dat$outofschool,dat$healthspend        ,
                                         dat$healthstaff,dat$healthtraining     ,
                                         dat$educationspend,dat$femalepoverty   ,
                                         dat$year,dat$urban,dat$educationmunicip,
                                         dat$condom,dat$usingcont,
                                         dat$femaleworkers),
                                         function(vec) {sum(na.omit(vec))})
  names(formod) <- c(Names,"failures","successes")
  
  x  <- glm(cbind(successes,failures) ~ factor(dom_comuna) + factor(year)      +
                 factor(dom_comuna):trend + factor(pill) + factor(party)       + 
                 factor(mujer) + votes + outofschool + educationspend          + 
                 educationmunicip + healthspend + healthtraining + healthstaff +
                 femalepoverty + femaleworkers + condom,
                 family="binomial", data=formod)
  if (numb==1) {
    clusters <-mapply(paste,"dom_comuna.",dat$dom_comuna,sep="")
    x$coefficients2 <- robust.se(x,clusters)[[2]]
      
    n  <- sum(formod$successes) + sum(formod$failures)
    s1 <- pillest(x, formod, n, "pill", 1)
    db <- sum(formod$successes)/sum(formod$failures)
    db <- format(round(db,3),nsmall=3)
    c  <- sum(formod$successes)
    c  <- format(c,big.mark=",",scientific=F)

    return(list("beta" = s1$b, "se" = s1$s, "R2" = s1$r, "n" = s1$n, 
                "c" = c, "db" = db))
  } else {
    return(x)
  }
}

#******************************************************************************
#***(4c) Various functions to examine effect of spillover 
#******************************************************************************
closegen <- function(d1,d2,dat) {
  dat2 <- dat
  dat2$newvar <- NA  
  dat2$newvar[dat2$pilldistance > d1 & dat2$pilldistance <= d2 &
                !(dat2$pilldistance)==0] <- 1
  dat2$newvar[is.na(dat2$newvar)]<-0
  names(dat2)<-c(names(dat),paste('close',d2,sep=""))
  return(dat2)
}

spillovers <- function(age_sub,deathtype) {
  dat        <- orig
  dat        <- dat[dat$age %in% age_sub,]
  dat        <- dat[dat$pregnant == 1,]
  dat$Q      <- dat$earlyQ+dat$lateQ
  dat$early  <- dat$earlyQ+dat$earlyP
  dat$late   <- dat$lateQ+dat$lateP
  
  dat <- closegen(0,15,dat)
  dat <- closegen(15,30,dat)
  dat <- closegen(30,45,dat)
  dat <- dat[complete.cases(dat),]
  
  dat <- aggregate.data.frame(dat[,c("n",deathtype)],
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
  names(dat)    <- c("close15","close30","close45",Names,"n","death")
  dat$failures  <- dat$n 
  dat$successes <- dat$death
  
  mod <- aggregate.data.frame(dat[,c("failures","successes")],
                             by=list(dat$close15,dat$close30,dat$close45    ,
                                     dat$dom_comuna,dat$trend,dat$trend2    ,
                                     dat$pill,dat$mujer,dat$party,dat$votes ,
                                     dat$outofschool,dat$healthspend        ,
                                     dat$healthstaff,dat$healthtraining     ,
                                     dat$educationspend,dat$femalepoverty   ,
                                     dat$year,dat$urban,dat$educationmunicip,
                                     dat$condom,dat$usingcont,
                                     dat$femaleworkers),
                               function(vec) {sum(na.omit(vec))})
  names(mod) <- c("close15","close30","close45",Names,"failures","successes")
  
  
  xspill <- glm(cbind(successes,failures) ~ factor(dom_comuna) + factor(year) +
                  factor(dom_comuna):trend + factor(party) + factor(pill)     + 
                  factor(mujer) + votes + outofschool + educationspend        + 
                  educationmunicip + healthspend + healthtraining + condom    +
                  healthstaff + femalepoverty + femaleworkers                 + 
                  factor(close15) + factor(close30) + factor(close45),
                family="binomial",data=mod)

  clusters <-mapply(paste,"dom_comuna.",dat$dom_comuna,sep="")
  xspill$coefficients2 <- robust.se(xspill,clusters)[[2]]
  #xspill$newse<-vcovHC(xspill, type="HC")
  #xspill$coefficients2 <- coeftest(xspill,xspill$newse)
  
  n  <- sum(mod$successes) + sum(mod$failures)
  s1 <- pillest(xspill,mod,n,"pill|close",2)
  
  return(list("b" = s1$b,"s" = s1$s, "n" = s1$n, "r" = s1$r))
}

#******************************************************************************
#***(5) Estimate
#******************************************************************************
if(death){
  p1519 <- death_pmod(age_sub = 15:19, "death",1)
  p2034 <- death_pmod(age_sub = 20:34, "death",1)
  p3549 <- death_pmod(age_sub = 35:49, "death",1)
  e1519 <- death_pmod(age_sub = 15:19,"earlyP",1)
  e2034 <- death_pmod(age_sub = 20:34,"earlyP",1)
  e3549 <- death_pmod(age_sub = 35:49,"earlyP",1)
  l1519 <- death_pmod(age_sub = 15:19, "lateP",1)
  l2034 <- death_pmod(age_sub = 20:34, "lateP",1)
  l3549 <- death_pmod(age_sub = 35:40, "lateP",1)
}

if(spill){
  s1519 <- spillovers(age_sub = 15:19,"earlyP")
  s2034 <- spillovers(age_sub = 20:34,"earlyP")
  s3549 <- spillovers(age_sub = 35:40,"earlyP")
}

if(full) {
  full1519 <- death_pmod(age_sub = 15:19,"earlyP",2)
  full2034 <- death_pmod(age_sub = 20:34,"earlyP",2)
  full3549 <- death_pmod(age_sub = 35:49,"earlyP",2)
}
#******************************************************************************
#***(6) Export
#******************************************************************************
xvar <- 'Morning After Pill &'
xv2  <- 'Close $<15$ km &'
xv3  <- 'Close $30-60$ km &'
obs  <- 'Observations&'
R2   <- 'McFadden\'s $R^2$&'
sig  <- '$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01;'
dpb  <- 'Mean (deaths/live birth)&'
s    <- '\\\\'
a    <- '&'

if(death){
  to <- file(paste(tab.dir,"Deaths.tex", sep=""))
  writeLines(c('\\begin{table}[htpb!] \\centering',
               '\\caption{The Effect of the Morning After Pill on Fetal Deaths}',
               '\\label{TEENtab:PillDeath}',
               '\\begin{tabular}{@{\\extracolsep{5pt}}lccc}\\\\[-1.8ex]',
               '\\hline\\hline\\\\[-1.8ex]','& All & Early & Late \\\\',
               '& Deaths & Gestation & Gestation \\\\ \\midrule',
               '\\multicolumn{4}{l}{\\textsc{15-19 year olds}} \\\\',
               '&&&\\\\',
               paste(xvar,p1519$beta,a,e1519$beta,a,l1519$beta, s,sep=""),
               paste(a,p1519$s,a,e1519$s,a,l1519$s,s,sep=""),'&&&\\\\',
               paste(dpb,p1519$db,a,e1519$db,a,l1519$db,s,sep=""),
               paste(obs,p1519$n,a,e1519$n,a,l1519$n,s,sep=""),
               paste(R2,p1519$R2,a,e1519$R2,a,l1519$R2,s,sep=""),
               '&&&\\\\',
               '\\multicolumn{4}{l}{\\textsc{20-34 year olds}} \\\\',
               '&&&\\\\',
               paste(xvar,p2034$beta,a,e2034$beta,a,l2034$beta,s,sep=""),
               paste(a,p2034$s,a,e2034$s,a,l2034$s,s,sep=""),'&&&\\\\',
               paste(dpb,p2034$db,a,e2034$db,a,l2034$db,s,sep=""),
               paste(obs,p2034$n,a,e2034$n,a,l2034$n,s,sep=""),
               paste(R2,p2034$R2,a,e2034$R2,a,l2034$R2,s,sep=""),
               '&&&\\\\',
               '\\multicolumn{4}{l}{\\textsc{35-49 year olds}} \\\\',
               '&&&\\\\',
               paste(xvar,p3549$beta,a,e3549$beta,a,l3549$beta,s,sep=""),
               paste(a,p3549$s,a,e3549$s,a,l3549$s,s,sep=""),'&&&\\\\',
               paste(dpb,p3549$db,a,e3549$db,a,l3549$db,s,sep=""),
               paste(obs,p3549$n,a,e3549$n,a,l3549$n,s,sep=""),
               paste(R2,p3549$R2,a,e3549$R2,a,l3549$R2,s,sep=""),
               '\\hline \\hline \\\\[-1.8ex]',
               '\\multicolumn{4}{p{10cm}}{\\begin{footnotesize}\\textsc{Notes:}',
               'Total fetal deaths for each group are',
               paste(p1519$c,", ",p2034$c,", and ",p3549$c, sep=""),
               ' for 15-19, 20-34 and 35-49 year olds respectively.  All',
               'regressions include year and comuna fixed-effects, and',
               'comuna-specific trends.  Each regression also includes the',
               'full set of time varying controls described in table',
               '\\ref{TEENtab:PillPreg}.  Standard errors are clustered by',
               'comuna.',
               paste(sig,'\\end{footnotesize}}',sep=""),
               '\\normalsize\\end{tabular}\\end{table}'),to)

  close(to)
}

if(spill){  
  to <- file(paste(tab.dir,"Spillovers_B.tex", sep=""))
  writeLines(c('\\multicolumn{4}{l}{\\textsc{\\noindent Panel B: Fetal Deaths}}\\\\',
               '&&&\\\\',
               paste(xvar,s1519$b[1],a,s2034$b[1],a,s3549$b[1],s,sep=""),
               paste('&',s1519$s[1],a,s2034$s[1],a,s3549$s[1],s,sep=""),               
               paste(xv2,s1519$b[2],a,s2034$b[2],a,s3549$b[2],s,sep=""),
               paste('&',s1519$s[2],a,s2034$s[2],a,s3549$s[2],s,sep=""),
               '&&&\\\\',
               paste(obs,s1519$n,a,s2034$n,a,s3549$n,s,sep=""),
               paste(R2,s1519$r,a,s2034$r,a,s3549$r,s,sep=""),
               '\\hline \\hline \\\\[-1.8ex]',
               '\\multicolumn{4}{p{9.2cm}}{\\begin{footnotesize}\\textsc{Notes:}',
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

if(full) {
stargazer(full1519, full2034,  full3549,
          title="The Morning After Pill and Fetal Death: Full Covariates",
          align=TRUE, label="TEENtabDeathFull",omit.stat=c("LL","ser","f"),
          keep=c("pill","mujer","votes","outofschool","educationspend",
                 "educationmunicip","healthspend","healthtraining",
               "healthstaff","femalepoverty","femaleworkers"), 
          column.labels=c("15-19 year olds","20-34 year olds","35-49 year olds"),
          column.separate=(c(1,1,1)),
          out=paste(tab.dir, "DeathFullCovars.tex", sep=""),
          dep.var.labels="Fetal Death (0-20 Weeks)",
          covariate.labels=c("Morning After Pill","Female Mayor","Mayor's Support",
                             "Out of School","Total Education Spending", 
                             "Municipal Education Spending", "Health Spending",
                             "Health Training", "Health Staff", "Female Poverty",
                             "Female Workers"),
          notes="\\begin{footnotesize} \\textsc{Notes:} Each model is identical to 
          column (2) of table \\ref{TEENtab:PillDeath}.  A description of each 
          variable is also provided in table \\ref{TEENtab:PillPreg}.  Municipality
          dummies and trends and political party dummies have been omitted for 
          clarity. $^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01 
          \\end{footnotesize}", notes.align="l",
          notes.append=FALSE, table.placement="htpb!")
}

if(combine){
  spillA <- readLines(paste(tab.dir,"Spillovers_A.tex", sep=""))
  spillB <- readLines(paste(tab.dir,"Spillovers_B.tex", sep=""))  

  to <- file(paste(tab.dir,"Spillovers.tex", sep=""))
  writeLines(c(spillA,spillB),to)
  close(to)
}

