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

#*******************************************************************************
#***(1) Directories, libraries
#*******************************************************************************
proj.dir <- "~/universidades/Oxford/DPhil/Thesis/Teens/"

brth.dir <- paste(proj.dir, "Data/Nacimientos/", sep="")
com.dir  <- paste(proj.dir, "Data/Comunas/", sep="")
codB.dir <- paste(proj.dir, "Source/Births/",sep="")
codD.dir <- paste(proj.dir, "Source/Deaths/",sep="")
deth.dir <- paste(proj.dir, "Data/Deaths/", sep="")
graf.dir <- paste(proj.dir, "Figures/", sep="")
ma.dir   <- paste(proj.dir, "Data/PAE/", sep="")
outt.dir <- paste(proj.dir, "Tables/", sep="")
pol.dir  <- paste(proj.dir, "Data/Alcaldes/",sep="")
pop.dir  <- paste(proj.dir, "Data/Poblacion/proyecciones/DatCom/",sep="")

library("data.table")
library("doBy")
library("SDMTools")

create     <- TRUE
comunas    <- TRUE
kids       <- TRUE
tables     <- TRUE
pillgraph  <- TRUE
preggraph  <- TRUE
deathgraph <- TRUE
totgraph   <- TRUE
trends     <- TRUE

#*******************************************************************************
#***(2) Load required data
#*******************************************************************************
births <- read.csv(paste(brth.dir,"S1Data_granular_covars.csv", sep=""))
deaths <- read.csv(paste(deth.dir,"S1Data_deaths_covars.csv", sep=""))

deaths$conserv <- 0
deaths$conserv[deaths$party=="UDI"|deaths$party=="RN"]<-1
deaths$pilldistance[deaths$pill==1]<-0
deaths$healthspend <- deaths$healthspend/1000
deaths$educationspend <- deaths$educationspend/1000

names(deaths) <- c("dc","ye","ag","or","n","de","eP","eQ","lP","lQ","p","pill",
                   "pd","el","mu","pa","vo","os","ht","hs","hc","et","em","fp",
                   "fw","po","co")

if(create) {
  birth_y_range <- 2000:2011
  pill_y_range <- birth_y_range - 1
  age_range <- c(15,49)  
  week          <- 20
  pat           <- "P"
  
  fB <- paste(codB.dir,"BirthGenerate.R",sep="")
  fD <- paste(codD.dir,"DeathGenerate.R",sep="")
  source(fB)
  source(fD)  
  filenameB <- paste(brth.dir, 'S1Data_covars_20002011.csv' ,sep="")
  filenameD <- paste(outD.dir, 'S1Data_20002011.csv', sep="")
  prep_s1_data(age_range,usecom="FALSE",filenameB)
  prep_s1_data_deaths(age_range,week,pat,FALSE,filenameD)
  
}
#*******************************************************************************
#***(3) Main Functions
#*******************************************************************************
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
  
  
  deathc <- summaryBy(de ~ pill, FUN=sum, data=deaths)
  dnopill <- deathc$de.sum[1]
  dpill <- deathc$de.sum[2]  
  dtotal <- dpill + dnopill
  
  dtotal <- format(dtotal, digits=3, big.mark=",", scientific=F)
  dnopill <- format(dnopill, digits=3, big.mark=",", scientific=F)
  dpill <- format(dpill, digits=3, big.mark=",", scientific=F)
  
  
  return(list("bp"=pill, "bn"=nopill, "bt"=total, 
              "dp"=dpill, "dn"=dnopill, "dt"=dtotal))
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

#*******************************************************************************
#***(4) Run Summary Functions
#*******************************************************************************
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
#*******************************************************************************
#***(5) Export Results to Summary Table
### NOTE: Currently women and child level stats below are copy and pasted from
### estpost in the Stata file.  This is very ad hoc but is only as a stop-gap
### to get first draft written up (2014/01/12)
#*******************************************************************************
if(tables) {
a  <- "&&"
v1 <- "Poverty &&"
v2 <- "Conservative &&"
v3 <- "Education Spending &&"
v4 <- "Health Spending &&"
v5 <- "Out of School &&"
v6 <- "Female Mayor &&"
v7 <- "Female Poverty &&"
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
             paste(v2,paste(t(cs$mean["co.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["co.sd",1:7]),collapse=""),sep=""),
             paste(v3,paste(t(cs$mean["et.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["et.sd",1:7]),collapse=""),sep=""),
             paste(v4,paste(t(cs$mean["ht.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["ht.sd",1:7]),collapse=""),sep=""),
             paste(v5,paste(t(cs$mean["os.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["os.sd",1:7]),collapse=""),sep=""),
             paste(v6,paste(t(cs$mean["mu.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["mu.sd",1:7]),collapse=""),sep=""),
             paste(v7,paste(t(cs$mean["fp.mean",1:6]),collapse=""),sep=""),
             paste(a, paste(t(cs$sd["fp.sd",1:7]),collapse=""),sep=""),
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
             '\\hline \\hline \\\\[-1.8ex]',
             '\\multicolumn{5}{p{10cm}}{\\begin{footnotesize}\\textsc{Notes:}',
             'Group means are presented with standard deviations below in',
             'parentheses.  Poverty refers to the \\% of the municipality',
             'below the poverty line, conservative is a binary variable',
             'indicating if the mayor comes from a politically conservative',
             'party','health and education spending are measured in thousands',
             'of Chilean',
             'pesos, and pill distance measures the distance (in km) to the',
             'nearest municipality which reports prescribing emergency',
             'contraceptives.  Pregnancies are reported as \\% of all women',
             'giving live birth, while fetal deaths are reported per live',
             'birth.  All summary statistics are for the period 2006-2011.',
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
  deathc <- summaryBy(p + lP + eP ~ ye, FUN=sum, data=deathc)
  deathc$late <- deathc$lP/(deathc$p)
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
       pch=16, axes=FALSE, ylim=c(1700,2200), xlab="", ylab="", 
       type="b",col="black")
  axis(2, ylim=c(1700,2200),col="black",las=1)  ## las=1 makes horizontal labels
  mtext("Fetal Deaths",side=2,line=3.2)
  box()
  par(new=TRUE)
  
  plot(birthnum$year, birthnum$births, pch=15,  
       xlab="", ylab="", ylim=c(200000,250000), 
       axes=FALSE, type="b", col="darkgreen")
  mtext("Births",side=4,col="darkgreen",line=3.5) 
  axis(4, ylim=c(200000,250000), col="darkgreen",col.axis="darkgreen",las=1)
  axis(1,birthnum$year)
  mtext("Year \n",side=1,col="black",line=2.5)  
  
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
  f <- paste(brth.dir, "S1Data_covars_20002011.csv", sep="")
  orig <- read.csv(f)

  postscript(paste(graf.dir,"Trends1519.eps",sep=""),
             horizontal = FALSE, onefile = FALSE, paper = "special",
             height=7, width=9)
  trends1519 <-birthtrends(age_sub=15:19, orig)  
  nopillM <- trends1519[trends1519$nopill==1,]
  pillM   <- trends1519[trends1519$nopill==0,]
  plot(pillM$year,pillM$N, type="b",pch=16,col="darkgreen",
       ylim=c(min(trends1519$N)-1000,max(trends1519$N)+2000),
       ylab="Number of Live Births",xlab="Year")
  lines(nopillM$year,nopillM$N, type="b",pch=15)
  legend("topleft",legend=c("Did not give EC in 2010","Gave EC in 2010"),
         text.col=c("black","darkgreen"),pch=c(15,16),col=c("black","darkgreen"))
  note <- "Note: Some municipalities which did not give the EC pill in 2010 did give the EC
  pill in 2011 (and vice versa)."
  note <- sub("\n","",note)
  note <- gsub("  "," ",note)  
  mtext(note, side=1, line=4, adj=0, cex=0.8)  
  dev.off()

  postscript(paste(graf.dir,"Trends2034.eps",sep=""),
             horizontal = FALSE, onefile = FALSE, paper = "special",
             height=7, width=9)
  trends2034 <-birthtrends(age_sub=20:34, orig)
  nopillM <- trends2034[trends2034$nopill==1,]
  pillM   <- trends2034[trends2034$nopill==0,]
  plot(pillM$year,pillM$N, type="b",pch=16,col="darkgreen",
       ylim=c(min(trends2034$N)-1000,max(trends2034$N)+8000),
       ylab="Number of Live Births",xlab="Year")
  lines(nopillM$year,nopillM$N, type="b",pch=15)
  
  legend("topleft",legend=c("Did not give EC in 2010","Gave EC in 2010"),
         text.col=c("black","darkgreen"),pch=c(15,16),col=c("black","darkgreen"))
  note <- "Note: Some municipalities which did not give the EC pill in 2010 did give the EC
  pill in 2011 (and vice versa)."
  note <- sub("\n","",note)
  note <- gsub("  "," ",note)  
  mtext(note, side=1, line=4, adj=0, cex=0.8)  
  dev.off() 

  
  f <- paste(deth.dir, "S1Data_20002011.csv", sep="")
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
  note <- "Note: Some municipalities which did not give the EC pill in 2010 did give the EC
  pill in 2011 (and vice versa)."
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
