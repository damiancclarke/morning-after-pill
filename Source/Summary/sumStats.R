# sumStats.R v1.00               damiancclarke             yyyy-mm-dd:2014-01-01
#---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
#
# sumStats.R creates summary stats for women and municipalities in Chile by the-
# ir emergency contraceptive status.  The summary stats are extracted from two 
# source files:
#    > S1Data_granular_covars.csv - the file recording births and pill status
#    > S1Data_deaths_covars.csv   - the file recording fetal deaths and pill
#
# It also create graphical representations of trends. This data comes from the
# two pill files:
#    > PillDist.csv
#    > Pill_MinSal.csv

rm(list=ls())

#===============================================================================
#===(1) Libraries, directories
#===============================================================================
ipak <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}
packages <- c("data.table","doBy","ggplot2","SDMTools","lattice","texreg")
ipack(packages)

proj.dir <- "~/universidades/Oxford/DPhil/Thesis/Teens/"


brth.dir <- paste(proj.dir, "Data/Nacimientos/", sep="")
com.dir  <- paste(proj.dir, "Data/Comunas/"    , sep="")
codB.dir <- paste(proj.dir, "Source/Births/"   , sep="")
codD.dir <- paste(proj.dir, "Source/Deaths/"   , sep="")
deth.dir <- paste(proj.dir, "Data/Deaths/"     , sep="")
graf.dir <- paste(proj.dir, "Figures/"         , sep="")
ma.dir   <- paste(proj.dir, "Data/PAE/"        , sep="")
outt.dir <- paste(proj.dir, "Tables/"          , sep="")
pol.dir  <- paste(proj.dir, "Data/Alcaldes/"   , sep="")
pop.dir  <- paste(proj.dir, "Data/Population/" , sep="")

create     <- FALSE
comunas    <- FALSE
kids       <- FALSE
tables     <- FALSE
pillgraph  <- FALSE
preggraph  <- FALSE
deathgraph <- FALSE
totgraph   <- FALSE
trends     <- FALSE
sumplots   <- FALSE
ptest      <- TRUE 
distplots  <- FALSE

#===============================================================================
#===(2) Load required data
#===============================================================================
births <- read.csv(paste(brth.dir,"S1Data_granular_covars.csv", sep=""))
deaths <- read.csv(paste(deth.dir,"S1Data_deaths_covars.csv", sep=""))

deaths$conserv <- 0
deaths$conserv[deaths$party=="UDI"|deaths$party=="RN"]<-1
deaths$pilldistance[deaths$pill==1]<-0
deaths$healthspend <- deaths$healthspend/1000
deaths$educationspend <- deaths$educationspend/1000

#NOTE: names don't really matter here.  We only work with death (de), age (ag)
# year (ye), pill (pill) and other death classes (eP,eQ,lP,lQ)

names(deaths) <- c("dc","os","hs","hc","ht","et","em","fp","fw","po","ur","u1",
                   "ye","dn","re","co","uc","pill","mu","pa","vo","el","pd","p",
                   "ag","or","n","de","eP","eQ","lP","lQ","cn")

if(create) {
    birth_y_range <- 2000:2012
    pill_y_range <- birth_y_range - 1
    age_range <- c(15,49)  
    week          <- 20
    pat           <- "P"
    
    fB <- paste(codB.dir,"BirthGenerate.R",sep="")
    fD <- paste(codD.dir,"DeathGenerate.R",sep="")
    source(fB)
    source(fD)  
    filenameB <- paste(brth.dir, 'S1Data_covars_20002012.csv' ,sep="")
    filenameD <- paste(deth.dir, 'S1Data_20002012.csv', sep="")
    prep_s1_data(age_range,usecom="FALSE",filenameB)
    prep_s1_data_deaths(age_range,week,pat,FALSE,filenameD)
    
}
#===============================================================================
#===(3) Main Functions
#===============================================================================
comunaSum <- function() {
    comd <- deaths
    all  <- deaths
    comd <- do.call(data.frame, aggregate(. ~ pill, comd,
                                          function(x) c(mean = mean(x), sd=sd(x))))
    all$g <- 1
    all <- do.call(data.frame, aggregate(. ~ g, all,
                                         function(x) c(mean = mean(x), sd=sd(x))))
    
    #comd <- format(round(comd,3), nsmall=3)
    #all  <- format(round(all,3), nsmall=3)
    
    comd <- format(comd, digits=3, big.mark=",", scientific=F)
    all  <- format(all,  digits=3, big.mark=",", scientific=F)
  
    all  <- all[,grepl("pill.mean|pill.sd",colnames(all))==FALSE]
  
    msdd <- cbind(t(comd[1,]),'&',t(comd[2,]),"&",t(all))
    sdd <- cbind("(",gsub("&",")&(",msdd[grep(".sd",rownames(msdd)),]),") \\\\")
    md  <- cbind(msdd[grep(".mean|pill",rownames(msdd)),],"\\\\")
    
    
    n<- with(deaths, tapply(dc, pill, FUN = function(x) length(unique(x))))
    return(list("sd" = sdd, "mean" = md, "n"=n[2]))
}

childSum <- function() {
    dat <- births
    dat$totnonpreg <- (1-births$pregnant)*births$n
    dat$totpreg <- births$pregnant*births$n
    
    birthc <- summaryBy(totnonpreg + totpreg ~ pill, FUN=sum, data=dat)
    

    nopill <- birthc$totpreg.sum[1]
    pill <- birthc$totpreg.sum[2]  
    total <- pill + nopill
    
    total <- format(total, digits=3, big.mark=",", scientific=F)
    nopill <- format(nopill, digits=3, big.mark=",", scientific=F)
    pill <- format(pill, digits=3, big.mark=",", scientific=F)
    
    Tnopill <- birthc$totpreg.sum[1]+birthc$totnonpreg.sum[1]
    Tpill   <- birthc$totpreg.sum[2]  +birthc$totnonpreg.sum[2]  
    Ttotal  <- Tpill + Tnopill
    
    deathc <- summaryBy(de ~ pill, FUN=sum, data=deaths)
    dnopill <- deathc$de.sum[1]
    dpill <- deathc$de.sum[2]  
    dtotal <- dpill + dnopill
    
    dtotal <- format(dtotal, digits=3, big.mark=",", scientific=F)
    dnopill <- format(dnopill, digits=3, big.mark=",", scientific=F)
    dpill <- format(dpill, digits=3, big.mark=",", scientific=F)
    
    
    return(list("bp"=pill, "bn"=nopill, "bt"=total, 
                "dp"=dpill, "dn"=dnopill, "dt"=dtotal,
                "tp"=Tpill, "tn"=Tnopill, "tt"=Ttotal))
}

birthtrends <- function(age_sub,dat) {
    dat <- dat[dat$age %in% age_sub,]
    
    preg <-dat[dat$pregnant==1,]
    birt <-dat[dat$pregnant==0,]
    preg <-aggregate(preg$n,by=list(preg$year,preg$pill,preg$dom_comuna),sum)
    birt <-aggregate(birt$n,by=list(birt$year,birt$pill,birt$dom_comuna),sum)
    names(preg)<-c("year","pill","comuna","Npreg")
    names(birt)<-c("year","pill","comuna","Nbirth")
    
    total <-merge(preg,birt,by=c("comuna","year","pill"))
    total$birthrate <- total$Npreg/total$Nbirth
    
    total$pill2010<-0
    total$pill2010[(total$year==2010&total$pill==1)]<-1
    
    pill2010<-aggregate(total$pill2010, by=list(total$comuna), sum)
    names(pill2010)<-c("comuna","pill2010")
    total<-total[,c("year","comuna","Npreg","Nbirth","birthrate")]
    
    final<-merge(total,pill2010,by="comuna",all=T)
    
    trends <- aggregate(final$Npreg,by=list(final$pill2010,final$year),sum)
    names(trends) <- c("nopill","year","N")
    
    return(trends)  
}

deathtrends <- function(age_sub,dat) {
    dat <- dat[dat$age %in% age_sub,]
    
    preg <-dat[dat$pregnant==1,]
    deth <-dat[dat$earlyP==1,]
    preg <-aggregate(preg$n,by=list(preg$year,preg$pill,preg$dom_comuna),sum)
    deth <-aggregate(deth$n,by=list(deth$year,deth$pill,deth$dom_comuna),sum)
    names(preg)<-c("year","pill","comuna","Npreg")
    names(deth)<-c("year","pill","comuna","Ndeath")
    
    total <-merge(preg,deth,by=c("comuna","year","pill"))
    total$deathrate <- total$Ndeath/total$Npreg
    
    total$pill2010<-0
    total$pill2010[(total$year==2010&total$pill==1)]<-1
    
    pill2010<-aggregate(total$pill2010, by=list(total$comuna), sum)
    names(pill2010)<-c("comuna","pill2010")
    total<-total[,c("year","comuna","Npreg","Ndeath","deathrate")]
    
    final<-merge(total,pill2010,by="comuna",all=T)
  
    trends <- aggregate(final$Ndeath,by=list(final$pill2010,final$year),sum)
    names(trends) <- c("nopill","year","N")
    
    return(trends)
}

pilltest <- function(dat) {
    dat <- dat[dat$age %in% 15:44,]
    dat <- dat[(dat$order %in% 1:100) | !(dat$pregnant),]

    dat$failures <- (1-dat$pregnant)*dat$n
    dat$successes <- dat$pregnant*dat$n
    
    dat <- aggregate.data.frame(dat[,c("failures","successes")],
                                by=list(dat$outofschool,dat$healthspend        ,
                                        dat$healthstaff,dat$healthtraining     ,
                                        dat$educationspend,dat$educationmunicip,
                                        dat$femalepoverty,dat$femaleworkers    ,
                                        dat$urbBin,dat$year,dat$region         ,
                                        dat$density,dat$condom,dat$usingcont   ,
                                        dat$pill,dat$mujer,dat$party,dat$votop ,
                                        dat$pilldistance,dat$dom_comuna),
                                function(vec) {sum(na.omit(vec))})

    names(dat) <- c('outofschool','healthspend','healthstaff','healthtraining',
                    'educationspend','educationmunicip','femalepoverty'       ,
                    'femaleworkers','urban','year','region','density','condom',
                    'condomUse','pill','mujer','party','votop','pilldistance' ,
                    'comuna','noPreg','Preg')
    dat$healthstaff       <- dat$healthstaff/100000
    dat$healthspend       <- dat$healthspend/100000
    #dat$healthtraining    <- dat$healthtraining/100000
    #NOTE: CHECK HOW HEALTH TRAINING DEFINED
    dat$educationspend    <- dat$educationspend/100000
    dat$educationmunicip  <- dat$educationmunicip/100000
    dat$pill              <- dat$pill*100
    dat$pillclose[dat$pilldistance<45&dat$pilldistance!=0] <- 100
    dat$pillclose[is.na(dat$pillclose)] <- 0
    dat$conservative[dat$party=="UDI"|dat$party=="RN"] <- 1
    dat$conservative[is.na(dat$conservative)] <- 0

    Pmod.noReg <- lm(pill ~ outofschool + healthspend + healthstaff         +
                     educationspend + educationmunicip                      +
                     femalepoverty + femaleworkers + urban + condomUse      +
                     condom + mujer + conservative + votop + factor(year)   ,
                     data = dat)
    Pmod.noReg <- extract(Pmod.noReg, include.adjrs = FALSE, include.rmse = FALSE)
    Pmod.Reg   <- lm(pill ~ outofschool + healthspend + healthstaff         +
                     educationspend + educationmunicip                      +
                     femalepoverty + femaleworkers +  urban + condomUse     +
                     condom + mujer + conservative  + votop + factor(year)  +
                     factor(region), data = dat)
    Pmod.Reg <- extract(Pmod.Reg, include.adjrs = FALSE, include.rmse = FALSE)
    Dmod.noReg <- lm(pillclose    ~ outofschool + healthspend + healthstaff +
                     educationspend + educationmunicip                      +
                     femalepoverty + femaleworkers + urban + condomUse      +
                     condom + mujer + conservative  + votop + factor(year)  ,
                     data = dat)
    Dmod.noReg <- extract(Dmod.noReg, include.adjrs = FALSE, include.rmse = FALSE)
    Dmod.Reg   <- lm(pillclose    ~ outofschool + healthspend + healthstaff +
                     educationspend + educationmunicip                      +
                     femalepoverty + femaleworkers + urban + condomUse      +
                     condom + mujer + conservative  + votop + factor(year)  +
                     factor(region), data = dat)
    Dmod.Reg <- extract(Dmod.Reg, include.adjrs = FALSE, include.rmse = FALSE)

    results <- texreg(list(Pmod.noReg,Pmod.Reg,Dmod.noReg,Dmod.Reg),
                      file=paste(outt.dir,"PillChoice.tex",sep=""),
                      caption="Comuna Characteristics and Pill Decisions",
                      omit.coef="(region)|(year)|(Intercept)|(Adj.)|(RMS)",
                      stars = c(0.01, 0.05,0.1),
                      custom.coef.names=c(NA,"Out of School",
                          "Health Spending","Health Staff",
                          "Education Spending","Education Level",
                          "Female Poverty","Female Workers","Urban",
                          "Condom Use","Condom Availability","Female Mayor",
                          "Conservative Mayor","Vote Margin",NA,NA,NA,NA,NA,
                          NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                      custom.model.names=c("Pill$\\times$100",
                          "Pill$\\times$100","Close$\\times$100",
                          "Close$\\times$100"),booktabs=TRUE,
                      custom.gof.names=c("R-squared","Observations"),
                      use.packages=FALSE,caption.above=TRUE,
                      label="TEENtab:pillchoice",
                      custom.note=paste("\\begin{footnotesize}\\textsc{Notes:} ",
                                        "Refer to table \\ref{TEENtab:SumStats}",
                                        " for variable definitions.            ",
                                        "$^{*}$p$<$0.1; $^{**}$p$<$0.05;       ",
                                        "$^{***}$p$<$0.01",
                                        "\\end{footnotesize}",sep=""))
    return(results)
}
if(ptest) {
    tester    <- pilltest(births)
    regfile   <- readLines(paste(outt.dir,"PillChoice.tex",sep=""),-1)
    regfile[8]  = "&(1)&(2)&(3)&(4) \\\\ \\midrule"
    regfile[38] = "Year FE&Y&Y&Y&Y\\\\ Region FE &&Y&&Y\\\\ \\bottomrule"
    writeLines(regfile,paste(outt.dir,"PillChoice.tex",sep=""))
}

#===============================================================================
#===(4) Run Summary Functions
#===============================================================================
if(comunas) {
    cs <- comunaSum()
}
if(kids) {
    ks <- childSum()
} 


pm   <-format(round(wt.mean(births$pregnant,births$n),3), nsmall=3)
psd  <-format(round(wt.sd(births$pregnant,births$n),3), nsmall=3)
pmp  <-format(round(wt.mean(births$pregnant[births$pill==1],
                            births$n[births$pill==1]),3), nsmall=3)
psdp <-format(round(wt.sd(births$pregnant[births$pill==1],
                          births$n[births$pill==1]),3), nsmall=3)
pmn  <-format(round(wt.mean(births$pregnant[births$pill==0],
                            births$n[births$pill==0]),3), nsmall=3)
psdn <-format(round(wt.sd(births$pregnant[births$pill==0],
                          births$n[births$pill==0]),3), nsmall=3)

#===============================================================================
#===(5) Export Results to Summary Table
### NOTE: Currently women and child level stats below are copy and pasted from
### estpost in the Stata file. 
#===============================================================================
if(tables) {
a  <- "&&"
v1 <- "Poverty &&"
v2 <- "Conservative &&"
v3a <- "Education Spending (Total) &&"
v3b <- "Education Spending (Municipal) &&"
v4 <- "Health Spending &&"
v5 <- "Out of School &&"
v6 <- "Female Mayor &&"
v7 <- "Female Poverty &&"
v11 <- "Condom Use &&"
v8 <- "Pill Distance &&"
v9 <- "Live Births &&"
v10 <- "Fetal Deaths &&"

sumfile <- file(paste(outt.dir, "SummaryStats.tex", sep=""))
writeLines(c('\\begin{table}[htpb!] \\centering',
             '\\caption{Summary Statistics} \\label{TEENtab:SumStats}',
             '\\begin{tabular} {@{\\extracolsep{5pt}}lp{3mm}ccc}\\\\ [-1.8ex]',
             '\\hline\\hline\\\\ [-1.8ex] &&No Pill&Pill&Total \\\\', 
             '&&Available&Available& \\\\ \\midrule ',
             '\\multicolumn{5}{l}{\\textsc{Municipality Characteristics}} \\\\',
             '&&&& \\\\',
             paste(v1,paste(t(cs$mean["po.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["po.sd",1:7]),collapse=""),sep=""),
             paste(v2,paste(t(cs$mean["cn.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["cn.sd",1:7]),collapse=""),sep=""),
             paste(v3a,paste(t(cs$mean["et.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["et.sd",1:7]),collapse=""),sep=""),
             paste(v3b,paste(t(cs$mean["em.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["em.sd",1:7]),collapse=""),sep=""),
             paste(v4,paste(t(cs$mean["hs.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["hs.sd",1:7]),collapse=""),sep=""),
             paste(v5,paste(t(cs$mean["os.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["os.sd",1:7]),collapse=""),sep=""),
             paste(v6,paste(t(cs$mean["mu.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["mu.sd",1:7]),collapse=""),sep=""),
             paste(v7,paste(t(cs$mean["fp.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["fp.sd",1:7]),collapse=""),sep=""),
             paste(v11,paste(t(cs$mean["co.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["co.sd",1:7]),collapse=""),sep=""),
             paste(v8,paste(t(cs$mean["pd.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["pd.sd",1:7]),collapse=""),sep=""),
             '&&&& \\\\',
             '\\multicolumn{5}{l}{\\textsc{Individual Characteristics}}\\\\',
             '&&&& \\\\',  
             paste(v9,pmn,'&',pmp,'&',pm,'\\\\',sep=""),
             paste(a,'(',psdn,')&(',psdp,')&(',psd,') \\\\',sep=""),
             paste(v10,paste(t(cs$mean["de.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["de.sd",1:7]),collapse=""),sep=""),
             'Birthweight &&3322.7&3334.3&3324.7\\\\',
             '&&     (540.0)&     (542.3)&     (540.4)\\\\',
             'Maternal education  &&       11.92&       12.03&       11.94\\\\',
             '&&     (2.967)&     (2.894)&     (2.955)\\\\',
             'Percent working     &&       0.295&       0.395&       0.312\\\\',
             '&&     (0.456)&     (0.489)&     (0.463)\\\\',
             'Married     &&       0.340&       0.309&       0.335\\\\',
             '&&     (0.474)&     (0.462)&     (0.472)\\\\',
             'Age at Birth      &&       27.05&       27.15&       27.07\\\\',
             '&&     (6.777)&     (6.790)&     (6.779)\\\\ \\midrule',
             paste("N Comunas && 346 &", cs$n,"& 346 \\\\",sep=""),
             paste("N Fetal Deaths &&",  ks$dn, "&", ks$dp,"&",  ks$dt, "\\\\",sep=""),
             paste("N Births &&",  ks$bn, "&", ks$bp,"&",  ks$bt, "\\\\",sep=""),
             paste("Population (sum) &&",  ks$tn, "&", ks$tp,"&",  ks$tt, "\\\\",sep=""),
             '\\hline \\hline \\\\[-1.8ex]',
             '\\multicolumn{5}{p{10cm}}{\\begin{footnotesize}\\textsc{Notes:}',
             'Group means are presented with standard deviations below in',
             'parentheses.  Poverty refers to the \\% of the municipality',
             'below the poverty line, conservative is a binary variable',
             'indicating if the mayor comes from a politically conservative',
             'party (UDI or RN), health and education spending are measured',
             ' in thousands of Chilean',
             'pesos, and pill distance measures the distance (in km) to the',
             'nearest municipality which reports prescribing emergency',
             'contraceptives.  Pregnancies are reported as \\% of all women',
             'giving live birth, while fetal deaths are reported per live',
             'birth.  All summary statistics are for the period 2006-2012.',
             '\\end{footnotesize}} \\normalsize\\end{tabular}\\end{table}'),
             sumfile)


close(sumfile)
}

#*******************************************************************************
#***(6) Graphical Results
#*******************************************************************************
if (pillgraph) {
    pillS <- read.csv(paste(ma.dir,"PillDist.csv", sep=""), sep=";")  
    pillM <- read.csv(paste(ma.dir,"Pill_MinSal.csv", sep=""), sep=";")    
    
  
    pillM$under18[is.na(pillM$under18)] <- 0
    pillM$over18[is.na(pillM$over18)] <- 0  
    pillM$total <- pillM$under18 + pillM$over18
    
    pillMcollapse <- summaryBy(total ~ Year, FUN=sum, data=pillM)  
    pillScollapse <- summaryBy(pill ~ year, FUN=sum, data=pillS)  
    
    
    preyears <- as.data.frame(t(rbind(2006:2008,0:0)))
    names(preyears)      <- c("year","Pill")
    names(pillMcollapse) <- c("year","Pill")
    names(pillScollapse) <- c("year","Pill")
    
    pillMcollapse <- rbind(preyears,pillMcollapse)
    pillScollapse <- rbind(preyears,pillScollapse)
    
    pillMcollapse <- pillMcollapse[pillMcollapse$year<=2011,]
    
    
    postscript(paste(graf.dir,"Pill.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    par(mar=c(5, 4, 4, 6) + 0.1)
    plot(pillMcollapse$year, pillMcollapse$Pill, 
         pch=16, axes=FALSE, ylim=c(0,8000), xlab="", ylab="", 
         type="b",col="black")
    axis(2, ylim=c(0,8000),col="black",las=1)  ## las=1 makes horizontal labels
    mtext("Quantity of Pills Prescribed",side=2,line=3)
    box()
    
    ## Allow a second plot on the same graph
    par(new=TRUE)
    
    plot(pillScollapse$year, pillScollapse$Pill, pch=15,  
         xlab="", ylab="", ylim=c(0,300), 
         axes=FALSE, type="b", col="darkgreen")
    ## a little farther out (line=4) to make room for labels
    mtext("Pill Municipalities",side=4,col="darkgreen",line=3) 
    axis(4, ylim=c(0,300), col="darkgreen",col.axis="darkgreen",las=1)
    
    ## Draw the time axis
    axis(1,pillScollapse$year)
    mtext("Year \n",side=1,col="black",line=2.5)  
    
    ## Add Legend
    legend("topleft",legend=c("Pills Prescribed","Pill Municipalities"),
           text.col=c("black","darkgreen"),pch=c(16,15),col=c("black","darkgreen"))
    
    note <- "Note: Prescription data is from the Ministry of Health's
    administrative data on medications and medical attention. 
    Municipality data is from an independent survey conducted by Dides et al, (2010;2011;2012)."
    note <- sub("\n","",note)
    note <- gsub("  "," ",note)  
    mtext(note, side=1, line=4, adj=0, cex=0.8)
    dev.off()
}

if (preggraph) {
    births$totnonpreg <- (1-births$pregnant)*births$n
    births$totpreg <- births$pregnant*births$n
    
    births$group[births$age>=15&births$age<=19] <- 1
    births$group[births$age>=20&births$age<=34] <- 2
    births$group[births$age>=35&births$age<=49] <- 3
    
    birthc <- summaryBy(totnonpreg + totpreg ~ year + group, FUN=sum, data=births)
    birthc$prop <- birthc$totpreg/(birthc$totnonpreg + birthc$totpreg)
    
    births1519<-birthc[birthc$group==1,]
    births2034<-birthc[birthc$group==2,]  
    
    
    postscript(paste(graf.dir,"Births.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    par(mar=c(5, 4, 4, 6) + 0.1)
    plot(births1519$year, births1519$prop, 
         pch=16, axes=FALSE, ylim=c(0.049,0.055), xlab="", ylab="", 
         type="b",col="black")
    axis(2, ylim=c(0.049,0.055),col="black",las=1)  ## las=1 makes horizontal labels
    mtext("Proportion of 15-19 year olds pregnant",side=2,line=3.2)
    box()
    
    par(new=TRUE)
    plot(births2034$year, births2034$prop, pch=15,  
         xlab="", ylab="", ylim=c(0.082,0.089), 
         axes=FALSE, type="b", col="darkgreen")
    mtext("Proportion of 20-34 year olds pregnant",side=4,col="darkgreen",line=3.5) 
    axis(4, ylim=c(0.082,0.089), col="darkgreen",col.axis="darkgreen",las=1)
    
    axis(1,births2034$year)
    mtext("Year \n",side=1,col="black",line=2.5)  
    legend("topleft",legend=c("15-19 year olds","20-34 year olds"),
           text.col=c("black","darkgreen"),pch=c(16,15),col=c("black","darkgreen"))
    
    note <- "Note: Data on pregnancies comes from the Ministry of Health's birth census"
    note <- sub("\n","",note)
    note <- gsub("  "," ",note)  
    mtext(note, side=1, line=4, adj=0, cex=0.8)
    dev.off()
}

if (deathgraph) {
    deathc       <- deaths[deaths$ag %in% 15:34,]
    deathc       <- summaryBy(p + lP + eP ~ ye, FUN=sum, data=deathc)
    deathc$late  <- deathc$lP/(deathc$p)
    deathc$early <- deathc$eP/(deathc$p)  
    
    postscript(paste(graf.dir,"Deaths.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    par(mar=c(5, 4, 4, 6) + 0.1)
    plot(deathc$ye, deathc$early, 
         pch=16, axes=FALSE, ylim=c(0.02,0.025), xlab="", ylab="", 
         type="b",col="black")
    axis(2, ylim=c(0.02,0.025),col="black",las=1)  ## las=1 makes horizontal labels
    mtext("Fetal Deaths (0-20 weeks)",side=2,line=3.2)
    box()
    
    par(new=TRUE)
    
    plot(deathc$ye, deathc$late, pch=15,  
         xlab="", ylab="", ylim=c(0.04,0.055), 
         axes=FALSE, type="b", col="darkgreen")
    mtext("Fetal Deaths (21+ weeks)",side=4,col="darkgreen",line=3.5) 
    axis(4, ylim=c(0.04,0.055), col="darkgreen",col.axis="darkgreen",las=1)
    
    axis(1,deathc$ye)
    mtext("Year \n",side=1,col="black",line=2.5)  
    legend("topleft",legend=c("Early term (0-20)","Late term (21-40)"),
           text.col=c("black","darkgreen"),pch=c(16,15),col=c("black","darkgreen"))
    
    note <- "Note: Data on pregnancies and fetal deaths comes from the Ministry of 
    Health's birth census"
    note <- sub("\n","",note)
    note <- gsub("  "," ",note)  
    mtext(note, side=1, line=4, adj=0, cex=0.8)
    dev.off()
}

if (totgraph) {
    births$births<-births$n*births$pregnant
    #total numbers cited in data section of paper
    sum(deaths$de)
    sum(births$births)
    
    birthnum<- summaryBy(births ~ year, FUN=sum, data=births)
    deathnum<- summaryBy(de ~ ye, FUN=sum, data=deaths)

    postscript(paste(graf.dir,"BirthDeath.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    par(mar=c(5, 4, 4, 6) + 0.1)
    plot(deathnum$ye, deathnum$de, 
         pch=16, axes=FALSE, ylim=c(1600,2200), xlab="", ylab="", 
         type="b",col="black")
    axis(2, ylim=c(1600,2200),col="black",las=1)
    mtext("Fetal Deaths",side=2,line=3.2)
    box()
    par(new=TRUE)
    
    plot(birthnum$year, birthnum$births, pch=15,  
         xlab="", ylab="", ylim=c(200000,250000), 
         axes=FALSE, type="b", col="darkgreen")
    mtext("Births",side=4,col="darkgreen",line=3.5) 
    axis(4, ylim=c(200000,250000), col="darkgreen",col.axis="darkgreen",las=1)
    axis(1,birthnum$year)
    mtext("Year",side=1,col="black",line=2.5)  
    
    legend("topleft",legend=c("Deaths","Births"),
           text.col=c("black","darkgreen"),pch=c(16,15),col=c("black","darkgreen"))
    
    note <- "Note: Data on pregnancies and fetal deaths comes from the Ministry of 
    Health's birth census"
    note <- sub("\n","",note)
    note <- gsub("  "," ",note)  
    mtext(note, side=1, line=4, adj=0, cex=0.8)  
    dev.off()
}

if (trends) {
    f <- paste(brth.dir, "S1Data_covars_20002012.csv", sep="")
    orig <- read.csv(f)
    
    postscript(paste(graf.dir,"Reform1519.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    trends1519 <-birthtrends(age_sub=15:19, orig)  
    nopillM <- trends1519[trends1519$nopill==1,]
    pillM   <- trends1519[trends1519$nopill==0,]
    plot(pillM$year,pillM$N, type="b",pch=16,col="darkgreen",
         ylim=c(min(trends1519$N)-1000,max(trends1519$N)+2000),
         ylab="Number of Live Births",xlab="Year")
    lines(nopillM$year,nopillM$N, type="b",pch=15)
    abline(v=2008.7)  
    legend("topleft",legend=c("Did not give EC in 2010","Gave EC in 2010"),
           text.col=c("black","darkgreen"),pch=c(15,16),col=c("black","darkgreen"))
    note <- "Note: Some municipalities which did not give the EC pill in 2010 did 
  give the EC pill in 2011 (and vice versa)."
    note <- sub("\n","",note)
    note <- gsub("  "," ",note)  
    mtext(note, side=1, line=4, adj=0, cex=0.8)  
    mtext(c("Constitutional Tribunal"),side=3,line=1, cex=0.8, at=2009)
    mtext(c("Change in 'Normas de Fertilidad'"),side=3,line=0, cex=0.8, at=2009)
    dev.off()
    
    postscript(paste(graf.dir,"Reform2034.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    trends2034 <-birthtrends(age_sub=20:34, orig)
    nopillM <- trends2034[trends2034$nopill==1,]
    pillM   <- trends2034[trends2034$nopill==0,]
    plot(pillM$year,pillM$N, type="b",pch=16,col="darkgreen",
         ylim=c(min(trends2034$N)-1000,max(trends2034$N)+8000),
         ylab="Number of Live Births",xlab="Year")
    lines(nopillM$year,nopillM$N, type="b",pch=15)
    abline(v=2008.7)
    
    legend("topleft",legend=c("Did not give EC in 2010","Gave EC in 2010"),
           text.col=c("black","darkgreen"),pch=c(15,16),col=c("black","darkgreen"))
    note <- "Note: Some municipalities which did not give the EC pill in 2010 did 
  give the EC pill in 2011 (and vice versa)."
    note <- sub("\n","",note)
    note <- gsub("  "," ",note)  
    mtext(note, side=1, line=4, adj=0, cex=0.8)  
    mtext(c("Constitutional Tribunal"),side=3,line=1, cex=0.8, at=2009)
    mtext(c("Change in 'Normas de Fertilidad'"),side=3,line=0, cex=0.8, at=2009)
    dev.off() 
    
    
    f <- paste(deth.dir, "S1Data_20002012.csv", sep="")
    orig <- read.csv(f)
    postscript(paste(graf.dir,"TrendsDeath1519.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    trends1519 <-deathtrends(age_sub=15:19, orig)  
    nopillM <- trends1519[trends1519$nopill==1,]
    pillM   <- trends1519[trends1519$nopill==0,]
    plot(pillM$year,pillM$N, type="b",pch=16,col="darkgreen",
         ylim=c(min(trends1519$N)-1000,max(trends1519$N)+2000),
         ylab="Number of Fetal Deaths",xlab="Year")
    lines(nopillM$year,nopillM$N, type="b",pch=15)
    legend("topleft",legend=c("Did not give EC in 2010","Gave EC in 2010"),
           text.col=c("black","darkgreen"),pch=c(15,16),col=c("black","darkgreen"))
    note <- "Note: Some municipalities which did not give the EC pill in 2010 did 
    give the EC pill in 2011 (and vice versa)."
    note <- sub("\n","",note)
    note <- gsub("  "," ",note)  
    mtext(note, side=1, line=4, adj=0, cex=0.8)  
    dev.off()
  
    postscript(paste(graf.dir,"TrendsDeath2034.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    trends2034 <-deathtrends(age_sub=20:34, orig)
    nopillM <- trends2034[trends2034$nopill==1,]
    pillM   <- trends2034[trends2034$nopill==0,]
    plot(pillM$year,pillM$N, type="b",pch=16,col="darkgreen",
         ylim=c(min(trends2034$N)-1000,max(trends2034$N)+8000),
         ylab="Number of Fetal Deaths",xlab="Year")
    lines(nopillM$year,nopillM$N, type="b",pch=15)
    
    legend("topleft",legend=c("Did not give EC in 2010","Gave EC in 2010"),
           text.col=c("black","darkgreen"),pch=c(15,16),col=c("black","darkgreen"))
    note <- "Note: Some municipalities which did not give the EC pill in 2010 did give the EC
    pill in 2011 (and vice versa)."
    note <- sub("\n","",note)
    note <- gsub("  "," ",note)  
    mtext(note, side=1, line=4, adj=0, cex=0.8)  
    dev.off() 
  
}

if (sumplots) {
    onlyBirths <- births[births$pregnant==1, c("age", "n")]
    histBirths <- onlyBirths[rep(row.names(onlyBirths), onlyBirths$n), ]
    b <- ggplot(histBirths, aes(x = age)) + xlab("Age") + ylab("Frequency") +
         geom_histogram(binwidth=1,fill=I("blue"),col=I("red"),alpha=I(.4))
    b + theme_bw() + theme(text = element_text(size = 15))
    ggsave(paste(graf.dir,"ageDistBirths.pdf",sep=""),width=9, height=7)

    onlyDeaths <- deaths[deaths$de!=0, c("ag", "de")]
    histDeaths <- onlyDeaths[rep(row.names(onlyDeaths), onlyDeaths$de), ]
    d <- ggplot(histDeaths, aes(x = ag)) + xlab("Age") + ylab("Frequency")  +
         geom_histogram(binwidth=1,fill=I("blue"),col=I("red"),alpha=I(.4))
    d + theme_bw() + theme(text = element_text(size = 15))
    ggsave(paste(graf.dir,"ageDistDeaths.pdf",sep=""),width=9, height=7)

    
    f <- paste(brth.dir, "S1Data_covars_20002012.csv", sep="")
    orig <- read.csv(f)

    orig$agegroup[orig$age>=15&orig$age<20] <- 1 
    orig$agegroup[orig$age>=20&orig$age<34] <- 2 
    orig$agegroup[orig$age>=35]             <- 3

    birthsC <- summaryBy(n ~ year+agegroup+pregnant, FUN=sum, data=orig)
    birthsC <- reshape(birthsC,timevar="pregnant",idvar=c("agegroup","year"),
                       direction="wide")
    names(birthsC) <- c("year","agegroup","popln","births")
    birthsC$ratio  <- birthsC$births/birthsC$popln*1000
    
    xyplot(birthsC$births~birthsC$year|birthsC$agegroup)
}

if (distplots) {
    distDat <- births[births$year>2008&births$year<2012,]

    distance <-aggregate(distDat$pilldistance,
                         by=list(distDat$year,distDat$dom_comuna),FUN="mean")

    postscript(paste(graf.dir,"EuclideanDistance.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    hist(distance$x[distance$x>0&distance$x<150], col="#CCCCFF", main="",
         xlab="Euclidean Distance (km)")
    dev.off()
    print(sum(distance$x>150))
    

    distance <-aggregate(distDat$roadDist,
                         by=list(distDat$year,distDat$dom_comuna),FUN="mean")
    distance$x <- distance$x/1000
    postscript(paste(graf.dir,"RoadDistance.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    hist(distance$x[distance$x>0&distance$x<140], col="#CCCCFF", main="",
         xlab="Distance over Roads (km)")
    print(sum(distance$x>150))
    dev.off()

    
    distance <-aggregate(distDat$travelTime,
                         by=list(distDat$year,distDat$dom_comuna),FUN="mean")
    distance$x <- distance$x/60
    postscript(paste(graf.dir,"TravelTime.eps",sep=""),
               horizontal = FALSE, onefile = FALSE, paper = "special",
               height=7, width=9)
    hist(distance$x[distance$x>0&distance$x<150], col="#CCCCFF", main="",
         xlab="Travel time (minutes)")
    print(sum(distance$x>150))
    dev.off()
    
}
