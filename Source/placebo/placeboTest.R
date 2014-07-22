# placeboTest.R v0.02            damiancclarke             yyyy-mm-dd:2014-07-10
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# Runs various placebo tests using lead births and gestational deaths as the ou-
# tcome instead of births and gestational deaths following the introduction of
# the EC pill.  Here instead of using outcomes in years following the availabil-
# ity of the pill as our outcome, we use outcomes in treatment and control muni-
# cipalities lagged by i\in{3,4,5,6} years.  If diff-in-diff assumptions are va-
# lid, we should expect that the effect of the pill on lagged births should not 
# be significantly different to zero.
#
# This code has been written by DCC.
#
# Last version 0.02: Updating directory structure.
# contact: damian.clarke@economics.ox.ac.uk

rm(list=ls())


#===============================================================================
#=== (1) Directories, libraries
#===============================================================================
proj.dir <- "~/universidades/Oxford/DPhil/Thesis/Teens/"

brth.dir <- paste(proj.dir, "Data/Nacimientos/",sep="")
codB.dir <- paste(proj.dir, "Source/Births/",sep="")
codD.dir <- paste(proj.dir, "Source/Deaths/",sep="")
com.dir  <- paste(proj.dir, "Data/Comunas/", sep="")
deth.dir <- paste(proj.dir, "Data/Deaths/",sep="")
ma.dir   <- paste(proj.dir, "Data/PAE/",sep="")
outt.dir <- paste(proj.dir, "Tables/", sep="")
pop.dir  <- paste(proj.dir, "Data/Population/",sep="")


require(sandwich)
require(lmtest)
require(foreign)
require(reshape)
require(gdata)

create <- TRUE
birthM <- TRUE
deathM <- TRUE
export <- TRUE

#===============================================================================
#=== (2) Import Raw Data with Pill comuna and all births/gestational deaths
#===============================================================================
if(create) {
  birth_y_range <- 2000:2011
  pill_y_range  <- birth_y_range - 1
  age_range     <- c(15,49)  
  week          <- 20
  pat           <- "P"
  filenameB     <- paste(brth.dir, 'S1Data_covars_20002011.csv', sep="")
  filenameD     <- paste(deth.dir, 'S1Data_20002011.csv', sep="")
  
  fb <- paste(codB.dir,"BirthGenerate.R",sep="")
  fd <- paste(codD.dir,"DeathGenerate.R",sep="")
  source(fb)
  source(fd)

  prep_s1_data(age_range,usecom="FALSE",filenameB)
  prep_s1_data_deaths(age_range,week,pat,FALSE,filenameD)
}


birth.dat <- read.csv(paste(brth.dir, "S1Data_covars_20002011.csv", sep=""))
death.dat <- read.csv(paste(deth.dir, "S1Data_20002011.csv", sep=""))

#===============================================================================
#=== (3) Functions for changing actual outcome year for pre-pill outcome
#===============================================================================
birthlag <- function(dat,lag) {
  dat$failures <- (1-dat$pregnant)*dat$n
  dat$successes <- dat$pregnant*dat$n

  dat<-aggregate.data.frame(dat[,c("failures","successes")],
                            by=list(dat$dom_comuna, dat$year, dat$age,
                                    dat$pill, dat$pilldistance      ),
                            function(vec) {sum(na.omit(vec))}       )
  names(dat) <- c("Com","year","age","pill","pilldist","failures","successes")

  births      <- dat[,c("Com","year","age","failures","successes")]
  pills       <- dat[,c("Com","year","age","pill","pilldist")]
  births$year <- births$year + lag

  out <- merge(births,pills,by=c("Com","year","age"),all.x=T,all.y=F)
  out <- out[!(is.na(out$pill)),]
  return(out)  
} 


deathlag <- function(dat,lag,deathtype) {
  dat <- dat[dat$pregnant == 1,]
  dat <- aggregate.data.frame(dat[,c("n",deathtype)],
                              by=list(dat$dom_comuna, dat$year,dat$age,
                                      dat$pill,dat$pilldistance      ),
                              sum)
  names(dat) <- c("Com","year","age","pill","pilldist","successes","failures")  

  deaths      <- dat[,c("Com","year","age","failures","successes")]
  pills       <- dat[,c("Com","year","age","pill","pilldist")]
  deaths$year <- deaths$year + lag
  
  out <- merge(deaths,pills,by=c("Com","year","age"),all.x=T,all.y=F)
  out <- out[!(is.na(out$pill)),]
  return(out)
} 

#===============================================================================
#=== (4) Auxiliary functions (cluster SEs, LaTeX stars, etc.)
#===============================================================================
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

closegen <- function(d1,d2,dat) {
  named <- names(dat)
  dat$newvar <- 0  
  dat$newvar[dat$pilldist > d1 & dat$pilldist <= d2 &
                !(dat$pilldist)==0] <- 1
  names(dat) <-c (named,paste('close',d2,sep=""))
  return(dat)
}

#===============================================================================
#=== (5a) Function to run birth DDs (pill and pilldist) but with placebos
#===============================================================================
basemodelB <- function(age_sub,dat,lag) {
  dat <- birthlag(dat,lag)
  dat <- dat[dat$age %in% age_sub,]
  dat <- dat[dat$year>=2006,]
  
  mod <- glm(cbind(successes,failures) ~ factor(Com) + factor(Com):year + 
               factor(pill) + factor(year)                              , 
               family="binomial", data=dat)
  n  <- sum(dat$successes) + sum(dat$failures)

  clusters <-mapply(paste,"Com.",dat$Com,sep="")
  mod$coefficients2 <- robust.se(mod,clusters)[[2]]
  
  results <- pillest(mod, dat, n,'pill', 1)
  return(results)
}

closemodelB <- function(age_sub,dat,lag) { 
  dat <- birthlag(dat,lag)
  dat <- dat[dat$age %in% age_sub,]
  dat <- dat[dat$year>=2006,]
  
  dat <- closegen(0,15,dat)
  dat <- closegen(15,30,dat)
  dat <- closegen(30,45,dat)
  
  mod <- glm(cbind(successes,failures) ~ factor(Com) + factor(Com):year        + 
               factor(pill) + factor(year) + factor(close15) + factor(close30) +
               factor(close45)                                                , 
               family="binomial", data=dat)
  n  <- sum(dat$successes) + sum(dat$failures)
  
  clusters <-mapply(paste,"Com.",dat$Com,sep="")
  mod$coefficients2 <- robust.se(mod,clusters)[[2]]
  
  results <- pillest(mod, dat, n,'pill|close', 4)
  return(results)
}

#===============================================================================
#=== (5b) Function to run death DDs (pill and pilldist) but with placebos
#===============================================================================
basemodelD <- function(age_sub,dat,lag,type) {
  dat <- deathlag(dat,lag,type)
  dat <- dat[dat$age %in% age_sub,]
  dat <- dat[dat$year>=2006,]
  
  mod <- glm(cbind(successes,failures) ~ factor(Com) + factor(Com):year + 
               factor(pill) + factor(year)                              , 
             family="binomial", data=dat)
  n  <- sum(dat$successes) + sum(dat$failures)
  
  clusters <-mapply(paste,"Com.",dat$Com,sep="")
  mod$coefficients2 <- robust.se(mod,clusters)[[2]]
  
  results <- pillest(mod, dat, n,'pill', 1)
  return(results)
}

closemodelD <- function(age_sub,dat,lag,type) { 
  dat <- deathlag(dat,lag,type)
  dat <- dat[dat$age %in% age_sub,]
  dat <- dat[dat$year>=2006,]
  
  dat <- closegen(0,15,dat)
  
  mod <- glm(cbind(successes,failures) ~ factor(Com) + factor(Com):year + 
               factor(pill) + factor(year) + factor(close15)            , 
             family="binomial", data=dat)
  n  <- sum(dat$successes) + sum(dat$failures)
  
  clusters <-mapply(paste,"Com.",dat$Com,sep="")
  mod$coefficients2 <- robust.se(mod,clusters)[[2]]
  
  results <- pillest(mod, dat, n,'pill|close', 4)
  return(results)
}

#===============================================================================
#=== (6) Run Models
#===============================================================================
if(birthM) {
  lagBB3t <- basemodelB(15:19,birth.dat,3)
  lagBB4t <- basemodelB(15:19,birth.dat,4)
  lagBB5t <- basemodelB(15:19,birth.dat,5)
  lagBB3w <- basemodelB(20:34,birth.dat,3)
  lagBB4w <- basemodelB(20:34,birth.dat,4)
  lagBB5w <- basemodelB(20:34,birth.dat,5)
  
  lagCB3t <- closemodelB(15:19,birth.dat,3)
  lagCB4t <- closemodelB(15:19,birth.dat,4)
  lagCB5t <- closemodelB(15:19,birth.dat,5)
  lagCB3w <- closemodelB(20:34,birth.dat,3)
  lagCB4w <- closemodelB(20:34,birth.dat,4)
  lagCB5w <- closemodelB(20:34,birth.dat,5)  
}

if(deathM) {
  lagBD3te <- basemodelD(15:19,death.dat,3,"death")
  lagBD4te <- basemodelD(15:19,death.dat,4,"death")
  lagBD5te <- basemodelD(15:19,death.dat,5,"death")
  lagBD3we <- basemodelD(20:34,death.dat,3,"death")
  lagBD4we <- basemodelD(20:34,death.dat,4,"death")
  lagBD5we <- basemodelD(20:34,death.dat,5,"death")  

  lagCD3te <- closemodelD(15:19,death.dat,3,"death")
  lagCD4te <- closemodelD(15:19,death.dat,4,"death")
  lagCD5te <- closemodelD(15:19,death.dat,5,"death")
  lagCD3we <- closemodelD(20:34,death.dat,3,"death")
  lagCD4we <- closemodelD(20:34,death.dat,4,"death")
  lagCD5we <- closemodelD(20:34,death.dat,5,"death")   
}

#==============================================================================
#=== (7) Export results
#==============================================================================
xvar <- 'Morning After Pill &'
xv2  <- 'Close $<15$ km &'
xv3  <- 'Close 15-30 km &'
xv4  <- 'Close 30-45 km &'

obs  <- 'Observations&'
R2   <- 'McFadden\'s $R^2$&'
sig  <- '$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01.'

if(export){
  to <-file(paste(outt.dir,"Placebo.tex", sep=""))
  writeLines(c('\\begin{landscape}','\\begin{table}[!htbp] \\centering',
               '\\caption{Placebo Tests}',
               '\\label{TEENtab:Placebo}',
               '\\begin{tabular}{lcccccc}',
               '\\\\[-1.8ex]\\hline \\hline \\\\[-1.8ex] ',
               '&\\multicolumn{2}{c}{Lag = 3 years}',
               '&\\multicolumn{2}{c}{Lag = 4 years}',
               '&\\multicolumn{2}{c}{Lag = 5 years}',
               '\\\\ \\cmidrule(r){2-3} \\cmidrule(r){4-5} \\cmidrule(r){6-7}',
               '&(1)&(2)&(3)&(4)&(5)&(6)\\\\ \\hline',
               '\\textsc{Panel A: 15-19 Year-Olds} &&&&&& \\\\',
               ' & & & & & & \\\\',
               paste(xvar,lagBB3t$b,'&',lagCB3t$b[1],'&',lagBB4t$b,'&',
                     lagCB4t$b[1],'&',lagBB3t$b,'&',lagCB3t$b[1],'\\\\',sep=""),
               paste('&',lagBB3t$s,'&',lagCB3t$s[1],'&',lagBB4t$s,'&',
                     lagCB4t$s[1],'&',lagBB3t$s,'&',lagCB3t$s[1],'\\\\',sep=""),
               paste(xv2,'&',lagCB3t$b[2],'&&',lagCB4t$b[2],'&&',lagCB5t$b[2],
                     '\\\\', sep=""), 
               paste('&&',lagCB3t$s[2],'&&',lagCB4t$s[2],'&&',lagCB5t$s[2],
                     '\\\\', sep=""), 
               paste(xv3,'&',lagCB3t$b[3],'&&',lagCB4t$b[3],'&&',lagCB5t$b[3],
                     '\\\\', sep=""), 
               paste('&&',lagCB3t$s[3],'&&',lagCB4t$s[3],'&&',lagCB5t$s[3],
                     '\\\\', sep=""), 
               paste(xv4,'&',lagCB3t$b[4],'&&',lagCB4t$b[4],'&&',lagCB5t$b[4],
                     '\\\\', sep=""), 
               paste('&&',lagCB3t$s[4],'&&',lagCB4t$s[4],'&&',lagCB5t$s[4],
                     '\\\\', sep=""), 
               ' & & & & & & \\\\',
               paste(obs,lagBB3t$n,'&',lagCB3t$n,'&',lagBB4t$n,'&',
                     lagCB4t$n,'&',lagBB5t$n,'&',lagCB5t$n,'\\\\', sep=""), 
               paste(R2,lagBB3t$r,'&',lagCB3t$r,'&',lagBB4t$r,'&',
                     lagCB4t$r,'&',lagBB5t$r,'&',lagCB5t$r,'\\\\ \\midrule', sep=""), 
               '\\textsc{Panel A: 20-34 Year-Olds} &&&&&& \\\\',
               ' & & & & & & \\\\',
               paste(xvar,lagBB3w$b,'&',lagCB3w$b[1],'&',lagBB4w$b,'&',
                     lagCB4w$b[1],'&',lagBB3w$b,'&',lagCB3w$b[1],'\\\\',sep=""),
               paste('&',lagBB3w$s,'&',lagCB3w$s[1],'&',lagBB4w$s,'&',
                     lagCB4w$s[1],'&',lagBB3w$s,'&',lagCB3w$s[1],'\\\\',sep=""),
               paste(xv2,'&',lagCB3w$b[2],'&&',lagCB4w$b[2],'&&',lagCB5w$b[2],
                     '\\\\', sep=""), 
               paste('&&',lagCB3w$s[2],'&&',lagCB4w$s[2],'&&',lagCB5w$s[2],
                     '\\\\', sep=""), 
               paste(xv3,'&',lagCB3w$b[3],'&&',lagCB4w$b[3],'&&',lagCB5w$b[3],
                     '\\\\', sep=""), 
               paste('&&',lagCB3w$s[3],'&&',lagCB4w$s[3],'&&',lagCB5w$s[3],
                     '\\\\', sep=""), 
               paste(xv4,'&',lagCB3w$b[4],'&&',lagCB4w$b[4],'&&',lagCB5w$b[4],
                     '\\\\', sep=""), 
               paste('&&',lagCB3w$s[4],'&&',lagCB4w$s[4],'&&',lagCB5w$s[4],
                     '\\\\', sep=""), 
               ' & & & & & & \\\\',
               paste(obs,lagBB3w$n,'&',lagCB3w$n,'&',lagBB4w$n,'&',
                     lagCB4w$n,'&',lagBB5w$n,'&',lagCB5w$n,'\\\\', sep=""), 
               paste(R2,lagBB3w$r,'&',lagCB3w$r,'&',lagBB4w$r,'&',
                     lagCB4w$r,'&',lagBB5w$r,'&',lagCB5w$r,'\\\\  \\hline \\hline', 
                     sep=""), 
               '\\multicolumn{7}{p{17.2cm}}{\\begin{footnotesize}\\textsc{Notes:}',
               'All specifications are identical to those estimated in tables ',
               '\\ref{TEENtab:PillPreg} and \\ref{TEENtab:Spillover}.  However, ',
               'instead of using births 1 year subsequent to the reform ',
               'the outcome variable in each case is births and preceeding ',
               'the reform by lag$=l\\in{3,4,5}$ years, and hence entirely',
               'unaffected in both treatment and control municipalities.',
               paste(sig, '\\end{footnotesize}}', sep=""),
               '\\normalsize\\end{tabular}\\end{table}\\end{landscape}'),to)
  close(to)
}
