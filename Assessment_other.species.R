# ------ Script for running multispecies BDM and catch-MSY stock assessment on other shark species---- ###################

#note: Catch-MSY implementation of Martell & Froese 2012.
#       total catches must be used (i.e. all sources of F)

#       If catches have never been >1% carrying capacity, then unexploited status so catch series have
#       no information on productivity.

#       preliminary analysis done to determine species grouping based on catch tonnage and number of
#       years of reported catch data

#       Taiwanese pelagic gillnet (1974-1986) and longline (1989-91) catch (Stevens & Davenport 1991; Stevens 1999)
#       are accounted for by setting initial depletion at <1

#Index:
#       1. Plot catch data as reported in logbooks
#       2. Reapportion catch of hammerhead and "shark,other'
#       3. Data manipulations
#       4. Get total catch
#       5. Add reapportioned catch
#       6. Select species with enough data
#       7. Species grouping
#       8. Life history parameters of selected species  
#       9. RESILIENCE list  
#       10. Catch-MSY
#       11. Multi-species BDM



rm(list=ls(all=TRUE))
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")
source.hnld="C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/"
fn.source=function(script)source(paste(source.hnld,script,sep=""))
fn.source("fn.fig.R")
fn.source("Leslie.matrix.R") 
fn.source("Catch_MSY.R")
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
Do.jpeg="YES"
Do.tiff="NO"

library(MASS)
library(plotrix)
library(PBSmapping)
#library(MASS)
library(dplyr)
library(tidyr)
set.seed(999) 


Asses.year=2019    #enter year of assessment

hNdl=paste("C:/Matias/Analyses/Population dynamics/1.Other species/",Asses.year,sep="")
fnkr8t=function(x) if(!dir.exists(x))dir.create(x)
fnkr8t(hNdl)
fnkr8t(paste(hNdl,"Outputs",sep="/"))



#---DATA SECTION-----

#Total catch
setwd("C:/Matias/Analyses/Data_outs")

  #WA Recreational catch
source("C:/Matias/Analyses/Population dynamics/Git_Stock.assessments/Recreational.catch.recons.R")

  #Commercial catch shark fisheries TDGDLF and NSF
Data.monthly=read.csv("Data.monthly.csv",stringsAsFactors=F)
Data.monthly.north=read.csv("Data.monthly.NSF.csv",stringsAsFactors=F)
Data.monthly.north$LIVEWT.c=Data.monthly.north$LIVEWT  

  #Commercial shark catch reported in other WA fisheries
Data.request.1=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/1975_1987.csv",stringsAsFactors=F)
Data.request.2=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/1988_2002.csv",stringsAsFactors=F)
Data.request.3=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/2003_2016.csv",stringsAsFactors=F)
Data.request=rbind(Data.request.1,Data.request.2,Data.request.3)
Data.request=subset(Data.request,!METHOD%in%c("GN","LL"))
names(Data.request)[match("LIVEWT",names(Data.request))]="LIVEWT.c"

  #remove data for other gears that is already recorded in Data.monthly
TDGDLF.other.gears=Data.monthly%>%
                        filter(!METHOD%in%c("GN","LL"))%>%
                        mutate(dummy=paste(FINYEAR,MONTH,VESSEL,BLOCKX,METHOD))
NSF.other.gears=Data.monthly.north%>%
                        filter(!METHOD%in%c("GN","LL"))%>%
                        mutate(dummy=paste(FINYEAR,MONTH,VESSEL,BLOCKX,METHOD))
fn.subs=function(YEAR) substr(YEAR,start=3,stop=4)
Data.request=Data.request%>%
                  mutate(FINYEAR=ifelse(MONTH%in%1:6,paste(YEAR-1,"-",fn.subs(YEAR),sep=""),
                                        ifelse(MONTH%in%7:12,paste(YEAR,"-",fn.subs(YEAR+1),sep=""),NA)),
                         LAT=-as.numeric(substr(BLOCKX,1,2)),
                         LONG=100+as.numeric(substr(BLOCKX,3,4)),
                         zone=as.character(ifelse(LONG>=116.5 & LAT<=(-26),"Zone2",
                                           ifelse(LONG<116.5 & LAT<=(-33),"Zone1",
                                           ifelse(LAT>(-33) & LAT<=(-26) & LONG<116.5,"West",
                                           ifelse(LAT>(-26) & LONG<114,"Closed",
                                           ifelse(LAT>(-26) & LONG>=114 & LONG<123.75,"North",
                                           ifelse(LAT>(-26) & LONG>=123.75,"Joint",NA))))))),
                         dummy=paste(FINYEAR,MONTH,VESSEL,BLOCKX,METHOD))%>%
                  filter(!dummy%in%c(TDGDLF.other.gears$dummy,NSF.other.gears$dummy))
Data.request.south=Data.request%>%filter(LAT<=max(Data.monthly$LAT))
Data.request.north=Data.request%>%filter(LAT>max(Data.monthly$LAT))
rm(Data.request,Data.request.1,Data.request.2,Data.request.3)


Effort.monthly=read.csv("Annual.total.eff.days.csv",stringsAsFactors=F)
Effort.monthly.north=read.csv("Annual.total.eff_NSF.csv",stringsAsFactors=F)

#Taiwanese gillnet fishery (source Stevens 1999 and Stevens & Davenport 1991, catch in tonnes)
Taiwan.gillnet.ktch=data.frame(Year=1974:1986,
                               Total.Ktch=c(321.4,8997.6,6455.3,9970.7,5528.2,3282.1,5831.1,6694.7,
                                            5624.1,7589.9,6544.2,2929.5,2111.1),
                               Shark.Ktch=c(rep(NA,5),557.9,4234.5,4486.2,3639.4,4418.6,2731.4,2327.4,2393.9),
                               WA.Ktch=c(rep(NA,5),50,1250,720,800,790,400,10,70))  #From Figure 7 Stevens & Davenport
Mn.ktch=mean(Taiwan.gillnet.ktch$WA.Ktch/Taiwan.gillnet.ktch$Shark.Ktch,na.rm=T)
Taiwan.gillnet.ktch=Taiwan.gillnet.ktch%>%
                      mutate(Shark.Ktch=ifelse(is.na(Shark.Ktch),Total.Ktch*mean(Shark.Ktch/Total.Ktch,na.rm=T),
                                               Shark.Ktch),
                             WA.Ktch=ifelse(is.na(WA.Ktch),Shark.Ktch*mean(WA.Ktch/Shark.Ktch,na.rm=T),WA.Ktch))

Taiwan.gillnet.sp.comp=as.data.frame(matrix(c(.357,.166,.059),ncol=3))
colnames(Taiwan.gillnet.sp.comp)=c("Australian.blacktip.Shark.West.prop",
                                   "Spot.tail.Shark.West.prop","Hammerheads.West.prop")
Taiwan.gillnet.ktch=cbind(Taiwan.gillnet.ktch,Taiwan.gillnet.sp.comp) %>%
                      mutate('Australian blacktip shark'=Australian.blacktip.Shark.West.prop*WA.Ktch,
                             'Spot-tail shark'=Spot.tail.Shark.West.prop*WA.Ktch,
                             Hammerheads=Hammerheads.West.prop*WA.Ktch)%>%
                      dplyr::select(Year,'Australian blacktip shark','Spot-tail shark',Hammerheads)%>%
                      gather(Species,Ktch,-Year)


Taiwan.longline.ktch=data.frame(Year=1990:1991,Shark.Ktch=c(1700*11/(11+9),1700*9/(11+9)))%>%
                        mutate(WA.ktch=Shark.Ktch*Mn.ktch)

Taiwan.longline.sp.comp=data.frame(Species=rep(c("Spot-tail shark","Australian blacktip shark","Tiger shark",
                                              "Milk shark","Spinner shark",
                                            "Pigeye shark","Graceful shark","Hammerheads","Other"),2),
                                   Year=c(rep(1990,9),rep(1991,9)), 
                                   Percent=c(18.3,25.6,8.1,.4,6.9,17.9,1.2,2,19.5,5.6,79.7,0.7,2,4.5,1,0,0.1,6.4))
Taiwan.longline.ktch=Taiwan.longline.sp.comp%>%
                          left_join(Taiwan.longline.ktch,by="Year")%>%
                          mutate(Ktch=WA.ktch*Percent/100)%>%
                          dplyr::select(names(Taiwan.gillnet.ktch))
#detach("package:MASS", unload=TRUE)
  
#source shark bio
User="Matias"
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R")

#species codes
All.species.names=read.csv("C:/Matias/Analyses/Population dynamics/1.Other species/Species_names.csv")
b=read.csv("C:\\Matias\\Data\\Species.code.csv")

#list of life history param for demography
LH.par=read.csv("C:/Matias/Data/Life history parameters/Life_History_other_sharks.csv",stringsAsFactors=F)

#Temperature
TEMP=read.csv("C:/Matias/Data/Oceanography/SST.csv")

#Species scientific names for assessed species
Shark.species=5001:24900
School.shark= 17008
Indicator.species=c(17001,17003,18001,18003,18007)
Shar_other=22999

Scien.nm=data.frame(SPECIES=c(17008,8001,10001,13000,17006,18013,18014,18021,18022,18023,18026,
                              18029,19004,19001,19002,20000,23002),
                    Scien.nm=c("Galeorhinus galeus","Carcharias taurus","Isurus oxyrinchus","Orectolobidae",
                               "Hypogaleus hyugaensis","Carcharhinus sorrah","C. limbatus & C. tilstoni",
                               "Carcharhinus leucas","Galeocerdo cuvier","C. brevipinna","C. amboinensis","Negaprion acutidens",
                               "Sphyrna zygaena","S. lewini","S. mokarran","Squalus spp.","Pristiophorus cirratus"))

#species list for exploring spatial dist of catch
SP.list=list(Angels=24900,Bignose=18012,BlacktipReef=18036,Blacktips=18014,
             Blue=18004,Bull=18021,CommonSaw=23002,CreekWhaler=18035,Graceful=18033,
             GreyNurse=8001,GreyReef=18030,Hammerheads=19000,Lemon=18029,Milk=18006,
             Nervous=18034,OceanicWhitetip=18032,Pencil=17006,Pigeye=18026,Sawsharks=23900,
             School=17008,SevenGill=5001,ShortfinMako=10001,Silky=18008,Silvertip=18027,
             SixGill=5002,SouthernSawshark=23001,Spinner=18023,Spottail=18013,Spurdogs=20000,
             TawnyNurse=13010,Threshers=12000,Tiger=18022,White=10003,Wobbegongs=13000)

non.sharks=c(                           
  "Australian Salmon","Baldchin groper",                          
  "Buffalo Bream", "Boxfish", "Blue Groper" ,
  "Bonito","Boarfish (general)" ,"Dusky Morwong" ,"Flathead","Flounder (general)",
  "John Dory (general)", "Knife Jaw" ,
  "Leatherjacket (general)",                      
  "Mackerels" , "Moonlighter", "Mulloway","North west blowfish" ,                                                    
  "Parrotfish (general)","Pink snapper" , "Queen Snapper",                                                           
  "Rankin cod","Red-lipped Morwong" , "Red Snapper, Redfish, Bight Redfish, Nannygai",                           
  "Samson fish ","Southern Blue-fin Tuna","Sergeant Baker",                         
  "Spotted sweetlips","Stripped marlin" ,"Spanish mackerel", "Skipjack trevally",
  "WHALE",                                        
  "Unidentified","Yellow tailed kingfish", "Yellowfin tuna ","Gurnard Perch" )

non.commercial.sharks=c("Brown-banded catshark","Cobbler Wobbegong",
                        "Dwarf sawfish","Eagle ray","Fiddler ray","Freshwater sawfish",
                        "Green sawfish","Guitarfish & shovelnose ray","Narrow sawfish",
                        "Port Jackson","Spotted shovelnose","Stingrays","Tawny nurse shark",
                        "Whitespot shovelnose","Zebra shark")


#Bring in abundance data
fn.read=function(x) read.csv(paste('C:/Matias/Analyses/Data_outs',x,sep='/'),stringsAsFactors = F)
  #Naturalist abundance survey
Scal.hh.nat=fn.read('Scalloped hammerhead.Srvy.FixSt.csv')
Tiger.nat=fn.read('Tiger shark.Srvy.FixSt.csv')

  #Standardised cpue
Smuz.hh.tdgdlf_mon=fn.read('Hammerhead.annual.abundance.basecase.monthly_relative.csv') #assumed to be all smooth HH
Smuz.hh.tdgdlf_daily=fn.read('Hammerhead.annual.abundance.basecase.daily_relative.csv')
Spinr.tdgdlf_mon=fn.read('Spinner Shark.annual.abundance.basecase.monthly_relative.csv')
Spinr.tdgdlf_daily=fn.read('Spinner Shark.annual.abundance.basecase.daily_relative.csv')
Tiger.tdgdlf_mon=fn.read('Tiger Shark.annual.abundance.basecase.monthly_relative.csv')
Tiger.tdgdlf_daily=fn.read('Tiger Shark.annual.abundance.basecase.daily_relative.csv')
Wobi.tdgdlf_mon=fn.read('Wobbegongs.annual.abundance.basecase.monthly_relative.csv')
Wobi.tdgdlf_daily=fn.read('Wobbegongs.annual.abundance.basecase.daily_relative.csv')

#---PARAMETERS SECTION-----

Min.max.ktch=50  #minimum total tonnes for a species to be analysed
Min.yrs=3        #minimum years with catch records for species to be included in analysis    
Min.yr.ktch=10    #minimum tonnage per year for at least Min.yrs

#Reference points
FAO.tar.prop=1.2  #FAO target and limit proportions of Bmsy
FAO.lim.prop=0.8

#Life history parameters for selected species  
pup.sx.ratio=.5

#.. Catch-MSY arguments
Do.Ktch.MSY=F

  #simulatins
SIMS=5e4  

  #Assumed process error
ERROR=0.05   #is default. 
ERROR2=0.02

KMAX=100

  #depletion levels
STARTBIO=c(.8,.99)   #starting time series prior to any fishing
FINALBIO=c(.2,.8)   #very uncertain

  #r priors
NsimSS=1000
r.prior="USER"  #demography
r.prior2=NA    #uniform

  #Scenarios  
Nscen=3
SCENARIOS=vector('list',Nscen)
names(SCENARIOS)=c('BaseCase','S1','S2')
SCENARIOS$BaseCase=list(Error=ERROR,R.prior=r.prior,Initial.dep=STARTBIO)
SCENARIOS$S1=list(Error=ERROR2,R.prior=r.prior,Initial.dep=STARTBIO)
SCENARIOS$S2=list(Error=ERROR,R.prior=r.prior2,Initial.dep=STARTBIO)

  #Future projections
years.futures=5


#.. Surplus production arguments
#note: Only fitting species with species-specific abundance time series 
#     (e.g. wobbegongs comprise several species so not fitted)
#     Assumption, negligible exploitation at start of time series
#     Assumed some efficiency increase to account for fact that spinner
#     and hammerhead cpue increased with increase catch. Only up to 2% per
#     year assumed because standardisation partly account for this (e.g. change
#     in spatial distribution of fleet by considering block effect)

  #Initial harvest rate scenarios
HR.o.scens=c(0.01,0.02,0.05) 
HR_o.sd=0.005  #SD of HR likelihood (fixed)

  #Efficiency increase scenarios from 1995 on (done up to 1994 in cpue stand.)
Efficien.scens=c(0,.01,.02)

  #Proportional biomass (as proportion of K) at start of catch time series
B.init=0.99 #(fixed)

  #Estimate q
estim.q="YES"

  #Initial estimated par value
Init.r=list("Bull shark"=.1,"Lemon shark"=.1,"Pigeye shark"=.1,
            "Scalloped hammerhead"=.1,"Smooth hammerhead"=.2,
            "Spinner shark"=.1,"Spurdogs"=.1,"Tiger shark"=.1,
            "Wobbegongs"=.1)

MAX.CV=1.1

#---PROCEDURE SECTION-----

Average.Lat=rbind(subset(Data.monthly,select=c(SPECIES,LAT)),subset(Data.monthly.north,select=c(SPECIES,LAT)))
do.sp.table.WoE.paper="NO"
if(do.sp.table.WoE.paper=="YES")
{
  a=subset(Data.monthly,SPECIES%in%Shark.species,select=c(SPECIES,SNAME,LIVEWT.c))
  b=subset(Data.monthly.north,SPECIES%in%Shark.species,select=c(SPECIES,SNAME,LIVEWT.c))
  d=rbind(a,b)
  d=subset(d,!SPECIES==Shar_other)
  Tab=aggregate(LIVEWT.c~SPECIES,d,sum)
  Snm=d[!duplicated(d$SPECIES),-match('LIVEWT.c',names(d))]
  Tab=merge(Tab,Snm,by="SPECIES")
  Tab$prop=Tab$LIVEWT.c/sum(Tab$LIVEWT.c)
  Tab=Tab[rev(order(Tab$prop)),]
  Tab$CumSum=round(100*cumsum(Tab$prop),3)
  write.csv(Tab[,c("SNAME","CumSum")],"C:\\Matias\\Scientific manuscripts\\Population dynamics\\Weight of evidence_main commercial sharks\\Table1.csv",row.names=F)
  
}

A=table(Data.monthly$FINYEAR,Data.monthly$MONTH)
A[A>0]=1
Drop=names(which(rowSums(A)<12))
Data.monthly=subset(Data.monthly,!FINYEAR%in%Drop)
Data.monthly.north=subset(Data.monthly.north,!FINYEAR%in%Drop)
YEARS=sort(as.numeric(substr(unique(Data.monthly$FINYEAR),1,4)))  
Current=YEARS[length(YEARS)]   



#Explore spatial catch distribution to check if species reporting issues
fn.expl.sp.ktch=function(d1)
{
  if(nrow(d1)>0)
  {
    aa=aggregate(LIVEWT.c~BLOCKX+LAT+LONG,d1,sum)
    plot(aa$LONG,aa$LAT,pch=19,ylab="",ylim=c(-36,-9),xlim=c(111,129),
         col="steelblue",xlab="",cex=((aa$LIVEWT.c/max(aa$LIVEWT.c))^0.5)*3)
  }else   
  {
    plot(1,axes=F,ann=F,col="white")
  }
}
pdf(paste(hNdl,"/Outputs/reported.spatial.catch.pdf",sep=""))
for(l in 1:length(SP.list))
{
  par(mfcol=c(2,1),mar=c(1.5,2.5,1,.5),las=1,mgp=c(1,.7,0))
  fn.expl.sp.ktch(d1=subset(Data.monthly.north,SPECIES%in%SP.list[[l]]))
  mtext(names(SP.list)[l],3)
  fn.expl.sp.ktch(d1=subset(Data.monthly,SPECIES%in%SP.list[[l]]))
}
dev.off()

#Add catch from 'other methods'
Data.monthly=Data.monthly%>%
                dplyr::select(Same.return,FINYEAR,LIVEWT.c,SPECIES,SNAME,RSCommonName,
                       BLOCKX,LAT,LONG,zone)
Data.monthly.north=Data.monthly.north%>%
                dplyr::select(Same.return,FINYEAR,LIVEWT.c,SPECIES,SNAME,RSCommonName,
                       BLOCKX,LAT,LONG,zone)
RSCommonName=rbind(Data.monthly%>%distinct(SPECIES,.keep_all = T),
                   Data.monthly.north%>%distinct(SPECIES,.keep_all = T))%>%
                dplyr::select(SPECIES,RSCommonName)%>%
                  distinct(SPECIES,.keep_all = T)
Data.request.south=Data.request.south%>%
                    mutate(Same.return=paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))%>%
                dplyr::select(Same.return,FINYEAR,LIVEWT.c,SPECIES,SNAME,
                           BLOCKX,LAT,LONG,zone)%>%
                    left_join(RSCommonName,by="SPECIES")
Data.request.north=Data.request.north%>%
                    mutate(Same.return=paste(FINYEAR,MONTH,VESSEL,METHOD,BLOCKX))%>%
                dplyr::select(Same.return,FINYEAR,LIVEWT.c,SPECIES,SNAME,
                           BLOCKX,LAT,LONG,zone)%>%
                    left_join(RSCommonName,by="SPECIES")

Data.monthly=rbind(Data.monthly,Data.request.south)
Data.request.north=rbind(Data.request.north,Data.request.north)

#1. Plot catch data as reported in logbooks
ThIs=subset(Shark.species,!Shark.species%in%c(Indicator.species))
Data.monthly$Region="South"
Data.monthly.north$Region="North"
Tot.ktch=rbind(subset(Data.monthly,SPECIES%in%ThIs,select=c(FINYEAR,LIVEWT.c,SPECIES,SNAME,Region)),
               subset(Data.monthly.north,SPECIES%in%ThIs,select=c(FINYEAR,LIVEWT.c,SPECIES,SNAME,Region)))
Tot.ktch=merge(Tot.ktch,All.species.names,by="SPECIES",all.x=T)  
Tot.ktch$finyear=as.numeric(substr(Tot.ktch$FINYEAR,1,4))
Tot.ktch$LIVEWT.c=Tot.ktch$LIVEWT.c/1000   #in tonnes
Tot.ktch$Name=as.character(Tot.ktch$Name)
Tot.ktch$SNAME=as.character(Tot.ktch$SNAME)
Tot.ktch$Name=with(Tot.ktch,ifelse(SPECIES%in%c(22999,31000),"unidentified sharks",Name))
Agg=aggregate(LIVEWT.c~Name+finyear+Region,Tot.ktch,sum)
Agg.r=reshape(Agg, v.names = "LIVEWT.c", idvar = c("Name","Region"),timevar = "finyear", direction = "wide")
colnames(Agg.r)[3:ncol(Agg.r)]=substr(colnames(Agg.r)[3:ncol(Agg.r)],10,20)
PCH=rep(19,nrow(Agg.r))
COL=rep(1,nrow(Agg.r))
Agg.r=Agg.r[order(Agg.r$Name),]
Sp.fig.1=unique(Agg.r$Name)

#Commercial
fn.fig(paste(hNdl,'/Outputs/Reported_catch_all_species_commercial',sep=''),2400,2400) 
smart.par(n.plots=length(Sp.fig.1),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
for(i in 1:length(Sp.fig.1))
{
  d=subset(Agg.r,Name==Sp.fig.1[i])
  d.N=subset(d,Region=="North")
  d.S=subset(d,Region=="South")
  plot(as.numeric(names(d)[3:length(d)]),d[1,3:length(d)],pch=PCH[i],
       col='transparent',cex=.8,ann=F,ylim=c(0,max(d[,3:ncol(d)],na.rm=T)))
  if(nrow(d.N)>0) points(as.numeric(names(d.N)[3:length(d.N)]),d.N[1,3:length(d.N)],pch=PCH[i],type='o',col="grey60",cex=.8)
  if(nrow(d.S)>0) points(as.numeric(names(d.S)[3:length(d.S)]),d.S[1,3:length(d.S)],pch=PCH[i],type='o',col="grey25",cex=.8)
  mtext(paste(Sp.fig.1[i]),3,line=0.2,cex=0.8)  
}
plot(1:10,ann=F,axes=F,col='transparent')
legend('center',c("North","South"),lty=c(1,1),col=c("grey60","grey25"),lwd=2,bty='n',pch=19,cex=1.5)
mtext("Financial year",1,line=0.5,cex=1.5,outer=T)
mtext("Total catch (tonnes)",2,las=3,line=0.35,cex=1.5,outer=T)
dev.off()

#Plot recreational catch
Rec.ktch=Rec.ktch%>%mutate(Region=ifelse(zone%in%c('Gascoyne','North Coast'),'North','South'),
                           year=as.numeric(substr(FINYEAR,1,4)))%>%
                    filter(year>=min(Agg$finyear))
Rec.sp=unique(Rec.ktch$Common.Name)
fn.fig(paste(hNdl,'/Outputs/Reported_catch_all_species_recreational',sep=''),2400,2400) 
smart.par(n.plots=length(Rec.sp),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
for(i in 1:length(Rec.sp))
{
  d=subset(Rec.ktch,Common.Name==Rec.sp[i])
  d.N=d%>%
    filter(Region=="North")%>%
    group_by(year)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c/1000))
  d.S=d%>%
    filter(Region=="South")%>%
    group_by(year)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c/1000))
  
  plot(sort(unique(d$year)),sort(unique(d$year)),col='transparent',cex=.8,ann=F,ylim=c(0,max(c(d.N$LIVEWT.c,d.S$LIVEWT.c,na.rm=T))))
  if(nrow(d.N)>0) points(d.N$year,d.N$LIVEWT.c,pch=PCH[i],type='o',col="grey60",cex=.8)
  if(nrow(d.S)>0) points(d.S$year,d.S$LIVEWT.c,pch=PCH[i],type='o',col="grey25",cex=.8)
  mtext(paste(Rec.sp[i]),3,line=0.2,cex=0.8)  
}
plot(1:10,ann=F,axes=F,col='transparent')
legend('center',c("North","South"),lty=c(1,1),col=c("grey60","grey25"),lwd=2,bty='n',pch=19,cex=1.5)
mtext("Financial year",1,line=0.5,cex=1.5,outer=T)
mtext("Total catch (tonnes)",2,las=3,line=0.35,cex=1.5,outer=T)
dev.off()

#plot Taiwanese catch
Taiwan.gillnet.ktch$Method="Pelagic.gillnet"
Taiwan.longline.ktch$Method="Longline"
Taiwan=rbind(Taiwan.longline.ktch,Taiwan.gillnet.ktch)%>%
              mutate(Region="North",
                     LIVEWT.c=1000*Ktch)%>%
              rename(year=Year)%>%
              arrange(Species,year)
sp.taiwan=unique(Taiwan$Species)
fn.fig(paste(hNdl,'/Outputs/Reported_catch_all_species_Taiwan',sep=''),2400,2400) 
smart.par(n.plots=length(sp.taiwan),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
for(i in 1:length(sp.taiwan))
{
  d=subset(Taiwan,Species==sp.taiwan[i])
  d.N=d%>%
    filter(Region=="North")%>%
    group_by(year)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c/1000))
  d.S=d%>%
    filter(Region=="South")%>%
    group_by(year)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c/1000))
  
  plot(sort(unique(Taiwan$year)),sort(unique(Taiwan$year)),col='transparent',cex=.8,ann=F,ylim=c(0,max(c(d.N$LIVEWT.c,d.S$LIVEWT.c,na.rm=T))))
  if(nrow(d.N)>0) points(d.N$year,d.N$LIVEWT.c,pch=PCH[i],type='o',col="grey60",cex=.8)
  if(nrow(d.S)>0) points(d.S$year,d.S$LIVEWT.c,pch=PCH[i],type='o',col="grey25",cex=.8)
  mtext(paste(sp.taiwan[i]),3,line=0.2,cex=0.8)  
  
}
legend('center',c("North","South"),lty=c(1,1),col=c("grey60","grey25"),lwd=2,bty='n',pch=19,cex=1.5)
mtext("Calendar year",1,line=0.5,cex=1.5,outer=T)
mtext("Total catch (tonnes)",2,las=3,line=0.35,cex=1.5,outer=T)
dev.off()



#2. Reapportion catch of hammerhead, blacktip and "shark,other'

  #2.1. Split of hammerhead catch by species
Split.HH="YES"
if(Split.HH=="YES")
{
  Smooth.hh.south.p=.97   #McAuley & Simpfendorfer 2003 ratios for TDGDLF
  Scalloped.hh.south.p=.03/2
  Great.hh.south.p=.03/2
  
  Smooth.hh.north.p=.01      # from Sharks database (Naturaliste trip)
  Scalloped.hh.north.p=.67
  Great.hh.north.p=.32
}

  #2.2 blacktip sharks
#note: reported all the way to Esperance, this doesn't conform to the species distribution
#      nor with observer data. Hence set to spinner shark any blacktip record east of Cape Leuwin (Last and Stevens) 
Data.monthly=Data.monthly%>%
                  mutate(SNAME=ifelse(SPECIES==18014 & LAT<(-30) & LONG>115.75,
                                      "SHARK, SPINNER (LONG-NOSE GREY)",SNAME),
                         SPECIES=ifelse(SPECIES==18014 & LAT<(-30) & LONG>115.75,
                                        18023,SPECIES),
                         RSCommonName=ifelse(SPECIES==18014 & LAT<(-30) & LONG>115.75,
                                        "Spinner Shark",RSCommonName))


  #2.3. Split 'shark, other' based on observers data of catch composition north and south
a=subset(DATA.bio,!is.na(BLOCK))
a=subset(a,!BLOCK==0)
Res.vess=c("NAT","HOU","FLIN","HAM","RV BREAKSEA","RV GANNET","RV SNIPE 2","RV Gannet")
a=subset(a,!BOAT%in%Res.vess)
a=subset(a,!is.na(BOAT))
a=subset(a,Method%in%c("GN","LL"))

a=subset(a,!COMMON_NAME%in%non.sharks)
a=subset(a,!is.na(COMMON_NAME))
a$N=1
non.commercial.sharks=subset(b,COMMON_NAME%in%non.commercial.sharks,select=c(CAES_Code,COMMON_NAME))
b=subset(b,!is.na(CAES_Code),select=c(Species,CAES_Code))
a=merge(a,b,by.x=c("SPECIES","CAES_Code" ),by.y=c("Species","CAES_Code"),all.x=T)
a=subset(a,!CAES_Code%in%c(25000:25010,13006,26999))
Agg.n.zone=aggregate(N~CAES_Code+zone,subset(a,!is.na(CAES_Code)),sum)
Agg.n.zone1=aggregate(N~zone,subset(a,!is.na(CAES_Code)),sum)
colnames(Agg.n.zone1)[2]="Total"
Agg.n.zone=merge(Agg.n.zone,Agg.n.zone1,by="zone",all.x=T)
Agg.n.zone$Prop=Agg.n.zone$N/Agg.n.zone$Total
Agg.n.zone=subset(Agg.n.zone,select=c(zone,CAES_Code,Prop))
Agg.n.zone=subset(Agg.n.zone,!CAES_Code%in%non.commercial.sharks$CAES_Code) #remove discarded species

Shark.OtheR=subset(Data.monthly,SPECIES%in%c(22999,31000))     
Shark.OtheR.north=subset(Data.monthly.north,SPECIES%in%c(22999,31000))
Shark.OtheR.north$zone=with(Shark.OtheR.north,ifelse(zone=="Closed",'North',zone))

Shark.OtheR=aggregate(LIVEWT.c~FINYEAR+zone,Shark.OtheR,sum)
Shark.OtheR.north=aggregate(LIVEWT.c~FINYEAR+zone,Shark.OtheR.north,sum)

Shark.OtheR=merge(Shark.OtheR,Agg.n.zone,by.x=c('zone'),by.y=c('zone'))
Shark.OtheR.north=merge(Shark.OtheR.north,Agg.n.zone,by.x=c('zone'),by.y=c('zone'))

names(Shark.OtheR)[3]=names(Shark.OtheR.north)[3]="weight"
names(Shark.OtheR)[4]=names(Shark.OtheR.north)[4]="SPECIES"
Shark.OtheR$LIVEWT.c=Shark.OtheR$weight*Shark.OtheR$Prop
Shark.OtheR.north$LIVEWT.c=Shark.OtheR.north$weight*Shark.OtheR.north$Prop
Shark.OtheR=aggregate(LIVEWT.c~FINYEAR+SPECIES,Shark.OtheR,sum)
Shark.OtheR.north=aggregate(LIVEWT.c~FINYEAR+SPECIES,Shark.OtheR.north,sum)

Shark.OtheR=subset(Shark.OtheR,!SPECIES%in%Indicator.species)
Shark.OtheR.north=subset(Shark.OtheR.north,!SPECIES%in%Indicator.species)


#3. Data manipulations
  #3.1 remove indicator species and 'shark, other' (after catch reapportion)
#note: in '#5. Add reapportioned catch" catch will be reapportioned accordingly to proportions in observed catch
Shark.species=subset(Shark.species,!Shark.species%in%c(Indicator.species,Shar_other,31000))


  #3.2 Remove blacktip sharks and school sharks (as per SAFS) 
#note: 
# blacktips is a shared with NT so NT assessment is used given their much higher catches (Grubert et al 2013)
# school sharks are assumed to be a shared stock in the SESSF assessment (Thomson & Punt 2009)
blacktips=subset(Scien.nm,Scien.nm%in%c('C. limbatus & C. tilstoni','Carcharhinus sorrah'))$SPECIES
Shark.species=subset(Shark.species,!Shark.species%in%c(blacktips,School.shark))


  #3.3 keep mapping data separate
# Map=subset(Data.monthly,SPECIES%in%Shark.species,select=c(Same.return,METHOD,
#         BLOCKX,LAT,LONG,Estuary,SPECIES,SNAME,FINYEAR,LIVEWT.c))
# Map.north=subset(Data.monthly.north,SPECIES%in%Shark.species,select=c(Same.return,METHOD,
#         BLOCKX,LAT,LONG,Estuary,SPECIES,SNAME,FINYEAR,LIVEWT.c))

  #3.4 keep 'shark, other' separately
Data.monthly.other=subset(Data.monthly,SPECIES%in%Shar_other,select=c(SPECIES,SNAME,FINYEAR,LIVEWT.c,BLOCKX))
Data.monthly.north.other=subset(Data.monthly.north,SPECIES%in%Shar_other,select=c(SPECIES,SNAME,FINYEAR,LIVEWT.c,BLOCKX))

Data.monthly=subset(Data.monthly,SPECIES%in%Shark.species,select=c(SPECIES,SNAME,FINYEAR,LIVEWT.c,BLOCKX))
Data.monthly.north=subset(Data.monthly.north,SPECIES%in%Shark.species,select=c(SPECIES,SNAME,FINYEAR,LIVEWT.c,BLOCKX))

Data.monthly.other$Region="South"
Data.monthly.north.other$Region="North"

Data.monthly$Region="South"
Data.monthly.north$Region="North"


#4. Get total catch (combine north and south)
Tot.ktch=rbind(Data.monthly,Data.monthly.north)
Tot.ktch.other=rbind(Data.monthly.other,Data.monthly.north.other)


#5. Add reapportioned catch
  #5.1 Shark,other
SNAMEs=Data.monthly[!duplicated(Data.monthly$SPECIES),match(c("SPECIES","SNAME"),names(Data.monthly))]
SNAMEs.north=Data.monthly.north[!duplicated(Data.monthly.north$SPECIES),match(c("SPECIES","SNAME"),names(Data.monthly.north))]

Shark.OtheR$Region="South"
Shark.OtheR.north$Region="North"

Shark.OtheR=merge(Shark.OtheR,SNAMEs,by="SPECIES",all.x=T)
Shark.OtheR.north=merge(Shark.OtheR.north,SNAMEs.north,by="SPECIES",all.x=T)

Shark.OtheR$BLOCKX=NA
Shark.OtheR.north$BLOCKX=NA


Tot.ktch=rbind(Tot.ktch,Shark.OtheR[,match(names(Tot.ktch),names(Shark.OtheR))],
               Shark.OtheR.north[,match(names(Tot.ktch),names(Shark.OtheR.north))])


  #5.2 Hammerheads  
if(Split.HH=="YES")
{
  Dat.hh=subset(Tot.ktch,SPECIES==19000 & Region=="South")
  Dat.hh.north=subset(Tot.ktch,SPECIES==19000 & Region=="North")
  Tot.ktch=subset(Tot.ktch,!SPECIES==19000)

  #south
  Dat.hh$Lat=-(as.numeric(substr(Dat.hh$BLOCKX,1,2)))
  Dat.hh$Lon=100+as.numeric(substr(Dat.hh$BLOCKX,3,4))
  a=b=d=Dat.hh
  a$SPECIES=19004 #CSIRO CAAB code for ID
  b$SPECIES=19001
  d$SPECIES=19002
  a$SNAME="SHARK, SMOOTH HH"
  b$SNAME="SHARK, SCALLOPED HH"
  d$SNAME="SHARK, GREAT HH"
  a$LIVEWT.c=with(a,ifelse(Lat<=(-26)&Lon<116,LIVEWT.c*Smooth.hh.south.p,LIVEWT.c))
  b$LIVEWT.c=with(b,ifelse(Lat<=(-26)&Lon<116,LIVEWT.c*Scalloped.hh.south.p,0))
  d$LIVEWT.c=with(b,ifelse(Lat<=(-26)&Lon<116,LIVEWT.c*Great.hh.south.p,0))
  Dat.hh=rbind(a,b,d)
  Dat.hh=Dat.hh[,-match(c('Lat','Lon'),names(Dat.hh))]
  
  #north
  Dat.hh.north$Lat=-(as.numeric(substr(Dat.hh.north$BLOCKX,1,2)))
  Dat.hh.north$Lon=100+as.numeric(substr(Dat.hh.north$BLOCKX,3,4))
  a=b=d=Dat.hh.north
  a$SPECIES=19004 
  b$SPECIES=19001
  d$SPECIES=19002
  a$SNAME="SHARK, SMOOTH HH"
  b$SNAME="SHARK, SCALLOPED HH"
  d$SNAME="SHARK, GREAT HH"
  
  a$LIVEWT.c=with(a,ifelse(Lat<=(-26)&Lon<116,LIVEWT.c*Smooth.hh.north.p,0))
  b$LIVEWT.c=b$LIVEWT.c*Scalloped.hh.north.p
  d$LIVEWT.c=d$LIVEWT.c*Great.hh.north.p
  Dat.hh.north=rbind(a,b,d)
  Dat.hh.north=Dat.hh.north[,-match(c('Lat','Lon'),names(Dat.hh.north))]
  
  Tot.ktch=rbind(Tot.ktch,Dat.hh,Dat.hh.north)
  
  #Reapportion Taiwanese catch
  Taiwan$Species=as.character(Taiwan$Species)
  drop.HH=subset(Taiwan,Species=="Hammerheads")
  nhh=nrow(drop.HH)
  drop.HH=rbind(drop.HH,drop.HH)%>%
    mutate(Prop=c(rep(Scalloped.hh.north.p/(Scalloped.hh.north.p+Great.hh.north.p),nhh),
                  rep(Great.hh.north.p/(Scalloped.hh.north.p+Great.hh.north.p),nhh)),
           LIVEWT.c=LIVEWT.c*Prop,
           Species=c(rep('Scalloped hammerhead',nhh),rep('Great hammerhead',nhh)))%>%
    dplyr::select(names(Taiwan))
  
  Taiwan=rbind(subset(Taiwan,!Species=="Hammerheads"),drop.HH)
}


#Some manipulations  
dummy=subset(All.species.names,SPECIES%in%unique(Tot.ktch$SPECIES))
Tot.ktch=merge(Tot.ktch,dummy,by="SPECIES",all.x=T)
Remaining_shark_other=subset(Tot.ktch,SPECIES==Shar_other)
Tot.ktch=subset(Tot.ktch,!SPECIES==Shar_other)
Percent.ktc.not.reapp=100*sum(Remaining_shark_other$LIVEWT.c)/sum(Tot.ktch$LIVEWT.c)
Tot.ktch$finyear=as.numeric(substr(Tot.ktch$FINYEAR,1,4))
Tot.ktch$LIVEWT.c=Tot.ktch$LIVEWT.c


Tot.ktch$Name=as.character(Tot.ktch$Name)
Tot.ktch$SNAME=as.character(Tot.ktch$SNAME)

  #5.3. Pool sawsharks because they are reported both by speces and by family
Tot.ktch$SPECIES=ifelse(Tot.ktch$SPECIES%in%c(23002,23001,23900),23900,Tot.ktch$SPECIES)    
Tot.ktch$Name=ifelse(Tot.ktch$Name%in%c("Southern sawshark","Common sawshark","Sawsharks"),"Sawsharks",Tot.ktch$Name)
Tot.ktch$SNAME=ifelse(Tot.ktch$SNAME%in%c("SHARK, COMMON SAW","shark, common saw","shark, southern saw","SHARK, SOUTHERN SAW",
                                          "SHARK, SAW","shark, saw"),"SHARK, SAW",Tot.ktch$SNAME)


  #5.4 Plot catch of all species (in tonnes) after reapportioning
Plot.yrs=sort(unique(Tot.ktch$finyear))
fn.add=function(D)
{
  id=Plot.yrs[which(!Plot.yrs%in%D$finyear)]
  if(length(id)>0)
  {
    A=as.data.frame(matrix(nrow=length(id),ncol=ncol(D)))
    colnames(A)=colnames(D)
    A$finyear=id
    D=rbind(D,A)
    D=D[order(D$finyear),]
  }
  return(D)
}
Pt.ktch.sp=function(sp,SP,LWD)
{
  d=subset(Tot.ktch,SPECIES%in%sp)
  
  b.Tot=aggregate(LIVEWT.c/1000~finyear,d,sum)
  uno=nrow(b.Tot)
  b.Tot=fn.add(D=b.Tot)
  if(uno>2)plot(1:length(b.Tot$finyear),b.Tot$"LIVEWT.c",type='l',lwd=LWD,ylab="",xlab="",xaxt='n',
                ylim=c(0,max(b.Tot$"LIVEWT.c",na.rm=T)),cex.axis=1.25)
  if(uno<=2)plot(1:length(b.Tot$finyear),b.Tot$"LIVEWT.c",pch=19,cex=2,ylab="",xlab="",xaxt='n',
                 ylim=c(0,max(b.Tot$"LIVEWT.c",na.rm=T)),cex.axis=1.25)
  
  d1=subset(d,Region=="North")
  if(nrow(d1)>0)
  {
    b.N=aggregate(LIVEWT.c/1000~finyear,d1,sum)
    uno=nrow(b.N)
    b.N=fn.add(D=b.N)
    if(uno>2)lines(1:length(b.N$finyear),b.N$"LIVEWT.c",col="grey50",lwd=LWD)
    if(uno<=2)points(1:length(b.N$finyear),b.N$"LIVEWT.c",pch=19,cex=2,col="grey50")
  }
  
  d1=subset(d,Region=="South")
  if(nrow(d1)>0)
  {
    b.S=aggregate(LIVEWT.c/1000~finyear,d1,sum)
    uno=nrow(b.S)
    b.S=fn.add(D=b.S)
    if(uno>2)lines(1:length(b.S$finyear),b.S$"LIVEWT.c",col="Grey75",lty=4,lwd=LWD)
    if(uno<=2)points(1:length(b.S$finyear),b.S$"LIVEWT.c",pch=19,cex=2,col="Grey75")
  }
  
  mtext(SP,3,0,cex=1.5)
  legend("topleft",c("Total",expression(paste("North of 26 ",degree,"S")),
                     expression(paste("South of 26 ",degree,"S"))),bty='n',lty=c(1,1,4),
         lwd=LWD,col=c("black","grey50","Grey75"),cex=1.5,pt.lwd=4)
  mtext("Catch (tonnes)",2,3,cex=2,las=3)
  mtext("Financial year",1,2.5,cex=2)
  axis(1,1:length(b.Tot$finyear),F,tck=-0.0125)
  axis(1,seq(1,length(b.Tot$finyear),10),Plot.yrs[seq(1,length(b.Tot$finyear),10)],cex.axis=1.25,tck=-0.025)
}

test=Tot.ktch[!duplicated(Tot.ktch$SPECIES),]
Uni.sp=test$SPECIES
names(Uni.sp)=test$Name


#Add recreational catch
dumi=Tot.ktch%>%distinct(SPECIES,.keep_all = T)%>%
            dplyr::select(SPECIES,Name,SNAME)%>%
                filter(Name%in%c("Scalloped hammerhead","Smooth hammerhead","Tiger shark",
                                 "Wobbegongs","Blacktip reef shark","Spinner shark"))%>%
                mutate(Common.Name=ifelse(Name=="Wobbegongs","Wobbegong",
                                   ifelse(Name=="Tiger shark","Tiger Shark",
                                   ifelse(Name=="Blacktip reef shark","Blacktip Reef Shark",
                                   ifelse(Name=="Spinner shark","Spinner Shark",
                                   Name)))))
Rec.ktch=Rec.ktch%>%
              rename(finyear=year)%>%
              mutate(BLOCKX=NA)%>%
              filter(Common.Name%in%dumi$Common.Name)%>%
              left_join(dumi,by='Common.Name')%>%
              dplyr::select(names(Tot.ktch))

Rec.ktch$Type="Recreational"
Tot.ktch$Type="Commercial"
Tot.ktch=rbind(Tot.ktch,Rec.ktch)


#Add Taiwanese catch
dumi2=Tot.ktch%>%distinct(SPECIES,.keep_all = T)%>%
             dplyr::select(SPECIES,Name,SNAME)%>%
             rename(Common.Name=Name)

Taiwan=Taiwan%>%
          rename(finyear=year,
                 Common.Name=Species)%>%
          mutate(BLOCKX=NA,
                 Type="Taiwan",
                 FINYEAR=paste(finyear,substr(finyear+1,3,4),sep="-"))%>%
          left_join(dumi2,by='Common.Name')%>%
          rename(Name=Common.Name)%>%
          dplyr::select(names(Tot.ktch))%>%
          filter(Name%in%unique(Tot.ktch$Name))
  
Tot.ktch=rbind(Tot.ktch,Taiwan)

  #Display Catch for each Species    
HnDl=paste(hNdl,"/Outputs/Catch_all_sp/",sep="")
fnkr8t(HnDl)
for(i in 1:length(Uni.sp))
{
  fn.fig(paste(HnDl,"Used.total.ktch",names(Uni.sp)[i],sep=""),2400,2400) 
  par(las=1,mgp=c(1,.8,0),mai=c(.8,1,.3,.1))
  Pt.ktch.sp(sp=Uni.sp[i],SP=names(Uni.sp)[i],LWD=3)
  dev.off()
}


  #5.5 Remove blacktips and school shark because some were reapportioned
Tot.ktch=subset(Tot.ktch,!Name%in%c('Blacktips','Spot tail shark',"School shark" ))


  #Species together   
Agg=aggregate(LIVEWT.c~Name+finyear,Tot.ktch,sum)
Agg.r=reshape(Agg, v.names = "LIVEWT.c", idvar = "Name",timevar = "finyear", direction = "wide")
colnames(Agg.r)[2:ncol(Agg.r)]=substr(colnames(Agg.r)[2:ncol(Agg.r)],10,20)
PCH=rep(19,nrow(Agg.r))
COL=rep(1,nrow(Agg.r))

id=rowSums(Agg.r[,2:ncol(Agg.r)],na.rm=T)
names(id)=Agg.r$Name
id=rev(sort(id))
Agg.r=Agg.r[match(names(id),Agg.r$Name),]


#6. Select species with enough data  

  #6.1 at least Min.max.ktch (total ever reported)
MAX.ktch=apply(Agg.r[,2:ncol(Agg.r)],1,max,na.rm=T)/1000
names(MAX.ktch)=as.character(Agg.r$Name)
Keep.species=names(MAX.ktch[which(MAX.ktch>=Min.max.ktch)])
MAX.ktch=rev(sort(MAX.ktch))

  #6.2 with annual catches of Min.yr.ktch for at least Min.yrs
Tab.sp=subset(Agg.r,Name%in%Keep.species)
rownames(Tab.sp)=Tab.sp$Name
Tab.sp=Tab.sp[,-1]
Tab.sp[Tab.sp<Min.yr.ktch]=0
Tab.sp[Tab.sp>=Min.yr.ktch]=1
Keep.species=rowSums(Tab.sp,na.rm=T)
Keep.species=names(Keep.species[Keep.species>=Min.yrs])
Tot.ktch=subset(Tot.ktch,Name%in%Keep.species)    

#7. Species grouping   
Tot.ktch$SP.group=Tot.ktch$Name

Specs=Tot.ktch[!duplicated(Tot.ktch$SP.group),match(c("SPECIES","SNAME","Name","SP.group"),names(Tot.ktch))]
Dis.sp=unique(Specs$SP.group)
Specs=Specs[order(Specs$SP.group),]
N.sp=nrow(Specs)


#8. Life history parameters of selected species   
LH.par=LH.par[order(LH.par$SPECIES),]
LH.par=merge(Scien.nm,LH.par,by="SPECIES",all.y=T)
LH.par.sawshars=subset(LH.par,SNAME=="shark, common saw")
LH.par.sawshars$SPECIES=23900
LH.par.sawshars$Scien.nm="Pristiophorus spp"
LH.par.sawshars$SNAME="SHARK, SAW"
LH.par=rbind(LH.par,LH.par.sawshars)

LH.par.hammerhead=subset(LH.par,SPECIES%in%c(19001,19002,19004))    
LH.par.hammerhead$SPECIES=19000
LH.par.hammerhead$SNAME="shark, hammerhead"
LH.par.hammerhead$Scien.nm="Sphyrna spp"
LH.par.hammerhead=LH.par.hammerhead[1,]
LH.par=rbind(LH.par,LH.par.hammerhead)

LH.par=subset(LH.par,SPECIES%in%Specs$SPECIES)   
LH.par=merge(subset(Specs,select=c(SPECIES,Name)),LH.par,by="SPECIES")
LH.par$SP.group=LH.par$Name
LH.par=LH.par[order(LH.par$SP.group),]

  
#fill in life history pars  
SPLF=LH.par$SP.group
GROWTH.F=vector('list',N.sp)
names(GROWTH.F)=SPLF
MAX.age.F=AGE.50.mat=FECU=Repro_cycle=Aver.Lat=AVER.T=GROWTH.F
for(i in 1:N.sp)
{
  dd=subset(LH.par,SP.group==SPLF[i])
  GROWTH.F[[i]]=data.frame(k=dd$K,FL_inf=dd$FL_inf)
  MAX.age.F[[i]]=c(dd$Max_Age,round(dd$Max_Age*1.1))
  AGE.50.mat[[i]]=c(dd$Age_50_Mat_min,dd$Age_50_Mat_max)
  FECU[[i]]=c(dd$Fecu_min,dd$Fecu_max)
  Repro_cycle[[i]]=LH.par$Cycle[i]
  AVER.T[[i]]=LH.par$Temperature[i]
}


setwd(paste(hNdl,'/Outputs',sep=''))

#export proportion of unreapportioned catch
write.csv(Percent.ktc.not.reapp,"Percent.ktc.not.reapp.csv",row.names=F)

#export Life history table
#LH.par=subset(LH.par,!SP.group=="Blacktips")
TabL=LH.par[,-match(c("SNAME","SP.group","Scien.nm"),names(LH.par))]
TabL=TabL[order(TabL$Name),-match("SPECIES",names(TabL))]
TabL$K=as.numeric(as.character(TabL$K))
TabL$K=round(TabL$K,3)
TabL$FL_inf=round(TabL$FL_inf)

fn.word.table(WD=getwd(),TBL=TabL,Doc.nm="Life history pars",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")


#9. RESILIENCE list                 
RESILIENCE=vector('list',N.sp)
names(RESILIENCE)=TabL$Name
#RESILIENCE$`Very.low`="Very low"
#RESILIENCE$Low="Low"
RESILIENCE$`Wobbegongs`="Low"
#RESILIENCE$`Blacktips`="Low"
RESILIENCE$`Tiger shark`="Low"
RESILIENCE$`Spinner shark`="Low"
RESILIENCE$`Bull shark`="Very low"
RESILIENCE$`Lemon shark`="Very low"
RESILIENCE$`Pigeye shark`="Very low"
RESILIENCE$`Scalloped hammerhead`="Low"
RESILIENCE$`Smooth hammerhead`="Low"
RESILIENCE$`Spurdogs`="Very low"

#Convert catch from kg to tonnes
Tot.ktch$LIVEWT.c=Tot.ktch$LIVEWT.c/1000

#Export total catch of each species
for(s in 1: N.sp)
{
  ddd=subset(Tot.ktch,SP.group==Specs$SP.group[s])
  ddd=aggregate(LIVEWT.c~finyear,ddd,sum)  
  write.csv(ddd,paste('Catch_all_sp/Total.annual.catch_',Specs$SP.group[s],
                      '.csv',sep=""),row.names = F)
}


#Displayed catches for analysed species
all.yrs=min(Tot.ktch$finyear):max(Tot.ktch$finyear)
COLs.type=c("black","grey60","white")
names(COLs.type)=c("Commercial","Recreational","Taiwan")
#Lins.type=c(1,3)
#names(Lins.type)=c("North","South")

fn.fig(paste(hNdl,'/Outputs/Figure1_catch_analysed_species',sep=''),2400,2400) 
smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
for(s in 1: N.sp)
{
  ddd=subset(Tot.ktch,SP.group==Specs$SP.group[s])%>%
            group_by(Type,finyear)%>%
            summarise(Tot=sum(LIVEWT.c,na.rm=T))%>%
            filter(!is.na(Tot))
  plot(all.yrs,all.yrs,col="transparent",ylab="",xlab="",main=Specs$SP.group[s],
       ylim=c(0,max(ddd$Tot)),xaxt='n')
  unik.T=unique(ddd$Type)
  #unik.R=unique(ddd$Region)
  for(u in 1:length(unik.T))
  {
   # for(r in 1:length(unik.R))
  #  {
      cl=COLs.type[match(unik.T[u],names(COLs.type))]
   #   lt=Lins.type[match(unik.R[r],names(Lins.type))]
      with(subset(ddd,Type==unik.T[u]),points(finyear,Tot,type='o',cex=1.5,pch=21,bg=cl))
#      with(subset(ddd,Type==unik.T[u] &Region==unik.R[r]),points(finyear,Tot,type='o',cex=1.5,pch=21,lty=lt,bg=cl))
   # }
  }
  if(s==1)legend('topleft',names(COLs.type),pt.bg=COLs.type,
                 bty='n',pch=21,cex=1.5,pt.cex=2)
  # if(s==1)legend('topleft',c(names(Lins.type),names(COLs.type)),
  #                lty=c(Lins.type,rep(NA,3)),pt.bg=c(1,1,COLs.type),
  #                lwd=2,bty='n',pch=21,cex=1,pt.cex=2)
  axis(1,all.yrs,F,tck=-.015)
  axis(1,seq(all.yrs[1],all.yrs[length(all.yrs)],10),
       seq(all.yrs[1],all.yrs[length(all.yrs)],10),tck=-.03)

}
mtext("Financial year",1,line=0.5,cex=1.5,outer=T)
mtext("Total catch (tonnes)",2,las=3,line=0.35,cex=1.5,outer=T)
dev.off()


# Build r prior -----------------------------------------------------------------------
fun.rprior.dist=function(Nsims,K,LINF,Temp,Amax,MAT,FecunditY,Cycle)
{
  Fecu=unlist(FecunditY)
  Rprior=fun.Leslie(N.sims=Nsims,k=K,Linf=LINF,Aver.T=Temp,
                    A=Amax,first.age=0,RangeMat=MAT,Rangefec=Fecu,
                    sexratio=0.5,Reprod_cycle=Breed.cycle,Hoenig.only="NO")  
  
  #get mean and sd from lognormal distribution
  normal.pars=suppressWarnings(fitdistr(Rprior, "normal"))
  gamma.pars=suppressWarnings(fitdistr(Rprior, "gamma"))  
  shape=gamma.pars$estimate[1]        
  rate=gamma.pars$estimate[2]      
  return(list(shape=shape,rate=rate,mean=normal.pars$estimate[1],sd=normal.pars$estimate[2]))
  
  #get mean and sd from lognormal distribution
  #LogN.pars=fitdistr(Rprior, "lognormal")  
  #log_mean.r=LogN.pars$estimate[1]    #already in log space     
  #log_sd.r=LogN.pars$estimate[2]      #already in log space     
  #return(list(log_mean.r=log_mean.r,log_sd.r=log_sd.r))
}
density.fun2=function(what,MAIN)
{
  #Prob of ref point
  f=ecdf(what)
  
  P.below.target=f(B.target)
  P.below.threshold=f(B.threshold)
  P.below.limit=f(B.limit)
  
  P.above.target=1-P.below.target
  P.above.threshold=1-P.below.threshold
  P.above.limit=1-P.below.limit
  
  P.between.thre.tar=P.below.target-P.below.threshold
  P.between.lim.thre=P.below.threshold-P.below.limit
  
  SEQ=seq(0,1,0.001)
  f.range=f(SEQ)
  plot(SEQ,f.range,ylab="",xlab="",type='l',lwd=2,cex.axis=1.25,main=MAIN,cex.main=1.3)
  abline(v=B.target,lty=2,col="grey60")
  abline(v=B.threshold,lty=2,col="grey60")
  abline(v=B.limit,lty=2,col="grey60")
  
  
  #Above target
  id=which.min(abs(SEQ - 1))
  id1=which.min(abs(SEQ - B.target))
  id=(id1+1):id
  X=SEQ[id]
  Y=f.range[id]
  polygon(c(X,rev(X)),c(Y,rep(0,length(Y))),col=CL.ref.pt[1],border="white")
  text(0.8,0.6,round(P.above.target,3),cex=1.5)
  
  #Between threshold & target
  id=which.min(abs(SEQ - B.target))
  id1=which.min(abs(SEQ - B.threshold))
  id=(id1+1):id
  X=SEQ[id]
  Y=f.range[id]
  polygon(c(X,rev(X)),c(Y,rep(0,length(Y))),col=CL.ref.pt[2],border="white")
  X=mean(c(B.target,B.threshold))
  text(X,0.6,round(P.between.thre.tar,3),srt=90,cex=1.5)
  #text(X,0.6,round(P.between.thre.tar,3),,srt=35,1,adj = c(0,.5))
  #arrows(X, mean(Y), X, 0.5, length = 0.1,col=1)
  
  #Between limit & threshold
  id=which.min(abs(SEQ - B.threshold))
  id1=which.min(abs(SEQ - B.limit))
  id=(id1+1):id
  X=SEQ[id]
  Y=f.range[id]
  polygon(c(X,rev(X)),c(Y,rep(0,length(Y))),col=CL.ref.pt[3],border="white")
  X=mean(c(B.limit,B.threshold))
  text(X,0.6,round(P.between.lim.thre,3),srt=90,cex=1.5)
  #text(X,0.6,round(P.between.lim.thre,3),srt=35,1,adj = c(0,.25),pos=3)
  #arrows(X, mean(Y), X, 0.5, length = 0.1,col=1)
  
  
  #Below limit
  id=which.min(abs(SEQ - B.limit))
  X=SEQ[1:id]
  Y=f.range[1:id]
  polygon(c(X,rev(X)),c(Y,rep(0,length(Y))),col=CL.ref.pt[4],border="white")
  lines(SEQ,f.range,lwd=2)
  X=B.limit/2
  text(X,0.6,round(P.below.limit,3),srt=90,cex=1.5)
  #text(X,0.6,round(P.below.limit,3),srt=35,1,adj = c(0,.5))
  #arrows(X, mean(Y), X, 0.5, length = 0.1,col=1)
  
  
  # d.frame=data.frame(P=c("P>Tar","Thre<P<Tar","Lim<P<Thre","P<Lim"),
  #                    Value=c(round(P.above.target,3),round(P.between.thre.tar,3),round(P.between.lim.thre,3),round(P.below.limit,3)))
  # addtable2plot(0,.5,d.frame,display.colnames=F,hlines=F,vlines=F,title="",bty="n",cex=.975,text.col="white",
  #               box.col="transparent",bg=rgb(.4,.4,.4,alpha=.75))
  
}
store.species=vector('list',N.sp)
names(store.species)=Specs$SP.group
PATH=paste(getwd(),"each_species",sep='/')
if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))   
setwd(file.path(PATH))
PHat=file.path(PATH)
system.time(for(s in 1: N.sp) #get r prior    #takes 0.013 sec per iteration 
{
  #life history
  Growth.F=GROWTH.F[[s]]
  TEMP=AVER.T[[s]]
  #TEMP=19
  Max.age.F=MAX.age.F[[s]]
  Age.50.mat=AGE.50.mat[[s]]
  Fecundity=FECU[[s]]
  Breed.cycle=Repro_cycle[[s]]  #years
  #Get r prior
  r.prior.dist=fun.rprior.dist(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf,Temp=TEMP,Amax=Max.age.F,
                               MAT=unlist(Age.50.mat),FecunditY=Fecundity,Cycle=Breed.cycle)
  store.species[[s]]$r.prior=list(shape=r.prior.dist$shape,rate=r.prior.dist$rate)
  store.species[[s]]$r.prior.normal=list(mean=r.prior.dist$mean,sd=r.prior.dist$sd)
})



# Catch-MSY ---------------------------------------------------------------
#takes 0.002 sec per iteration  
if(Do.Ktch.MSY)
{
  system.time(for(s in 1: N.sp)
  {
    #total catch
    ct=subset(Tot.ktch,SP.group==Specs$SP.group[s])
    ct=aggregate(LIVEWT.c~finyear,ct,sum)  
    names(ct)[match("LIVEWT.c",names(ct))]="Total.ktch"
    dummy.yrs=sort(unique(ct$finyear))
    id=YEARS[which(!YEARS%in%dummy.yrs)]
    if(length(id)>0)
    {
      dumb=ct[1:length(id),]
      dumb$Total.ktch=0
      dumb$finyear=id
      ct=rbind(ct,dumb)
      ct=ct[order(ct$finyear),]
    }
    
    #Tested scenarios
    ktch_msy_scen=vector('list',length(SCENARIOS))
    names(ktch_msy_scen)=names(SCENARIOS)
    for(sc in 1:length(ktch_msy_scen))
    {
      if(is.na(SCENARIOS[[sc]]$R.prior)) USR="No" else
        if(SCENARIOS[[sc]]$R.prior=="USER")USR= "Yes"
        Scen.start.bio=SCENARIOS[[sc]]$Initial.dep       
        ktch_msy_scen[[sc]]=list(r.prior=SCENARIOS[[sc]]$R.prior,user=USR,k.max=KMAX,
                                 startbio=Scen.start.bio,finalbio=FINALBIO,res=RESILIENCE[[s]],
                                 niter=SIMS,sigR=SCENARIOS[[sc]]$Error)
        if(!is.na(ktch_msy_scen[[sc]]$r.prior)) ktch_msy_scen[[sc]]$r.prior=unlist(store.species[[s]]$r.prior)
    }
    
    
    #Execute Catch-MSY  
    PATH=paste(PHat,Specs$SP.group[s],sep='/')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
    setwd(file.path(PATH))
    Path.ktch_msy=getwd()
    Ktch_MSY=ktch_msy_scen
    
    Yrs=ct$finyear
    Tot.Ktch=ct$Total.ktch
    yr.future=Current+(1:years.futures)
    ct.future=rep(mean(Tot.Ktch[(length(Tot.Ktch)-4):length(Tot.Ktch)]),years.futures)
    
    for(sc in 1:length(ktch_msy_scen))
    {
      Folder=names(ktch_msy_scen)[sc]
      if(!file.exists(paste(PATH,Folder,sep="/"))) dir.create(paste(PATH,Folder,sep="/"))
      setwd(paste(PATH,Folder,sep="/"))
      
      Ktch_MSY[[sc]]=Catch_MSY(ct=Tot.Ktch,
                               yr=Yrs,
                               r.prior=ktch_msy_scen[[sc]]$r.prior,
                               user=ktch_msy_scen[[sc]]$user,
                               k.max=ktch_msy_scen[[sc]]$k.max,
                               startbio=ktch_msy_scen[[sc]]$startbio,
                               finalbio=ktch_msy_scen[[sc]]$finalbio,
                               res=ktch_msy_scen[[sc]]$res,
                               n=ktch_msy_scen[[sc]]$niter,
                               sigR=ktch_msy_scen[[sc]]$sigR,
                               ct.future=ct.future,           
                               yr.future=yr.future)
      
      #Export outputs
      Table1_ktch_MSY=with(Ktch_MSY[[sc]],data.frame(`geom. mean r`,`r +/- 1.96 SD`,`geom. mean k (tons)`,`k +/- 1.96 SD (tons)`,
                                                     `geom. mean MSY (tons)`,`MSY +/- 1.96 SD (tons)`))
      write.csv(Table1_ktch_MSY,"Table1_ktch_MSY.csv",row.names=F)
      
    }
    store.species[[s]]$KTCH.MSY=Ktch_MSY
    store.species[[s]]$Catch=ct
    store.species[[s]]$K=Ktch_MSY[[sc]]$k
    
    Tabl.scen.Ktch.MSY=vector('list',length(ktch_msy_scen))
    for(i in 1:length(ktch_msy_scen))
    {
      dummy=ktch_msy_scen[[i]]
      for(a in 1:length(ktch_msy_scen[[i]])) if(length(dummy[[a]])>1) dummy[[a]]=paste(dummy[[a]],collapse=";")
      Tabl.scen.Ktch.MSY[[i]]=unlist(dummy)
    }
    Tabl.scen.Ktch.MSY=do.call(rbind,Tabl.scen.Ktch.MSY)
    row.names(Tabl.scen.Ktch.MSY)=names(ktch_msy_scen)
    write.csv(Tabl.scen.Ktch.MSY,"Scenarios.csv")
    
  })
}


# Single-species BDM -----------------------------------------------------------------------

  #get abundance data    
cpue.list=list(
  "Bull shark"=NULL,
  "Lemon shark"=NULL,
  "Pigeye shark"=NULL,
  "Scalloped hammerhead"=list(Nat=Scal.hh.nat,
                              TDGDLF.mon=NULL,
                              TDGDLF.day=NULL),
  "Smooth hammerhead"=list(Nat=NULL,
                           TDGDLF.mon=Smuz.hh.tdgdlf_mon,
                           TDGDLF.day=Smuz.hh.tdgdlf_daily), 
  "Spinner shark"=list(Nat=NULL,
                       TDGDLF.mon=Spinr.tdgdlf_mon,
                       TDGDLF.day=Spinr.tdgdlf_daily),
  "Spurdogs"=NULL,
  "Tiger shark"=list(Nat=Tiger.nat,
                     TDGDLF.mon=Tiger.tdgdlf_mon,
                     TDGDLF.day=Tiger.tdgdlf_daily),         
  "Wobbegongs"=NULL) 

Estimable.qs=list(
  "Bull shark"=NULL,
  "Lemon shark"=NULL,
  "Pigeye shark"=NULL,
  "Scalloped hammerhead"=c(q1=.005,NA,NA),
  "Smooth hammerhead"=c(NA,q2=.005,q3=.001), 
  "Spinner shark"=c(NA,q2=.005,q3=.001),
  "Spurdogs"=NULL,
  "Tiger shark"=c(q1=.005,q2=.005,q3=.001),         
  "Wobbegongs"=NULL) 



#ACA 
  #penalty for keeping biomass positive
posfun=function(x,eps,pen)
{
  if (x>=eps) return(x) else
  {
    pen=pen+.01*(x-eps)^2
    return (list(eps/(2-x/eps),pen))
  }
}
  #population dynamics function
SPM=function(Init.propK,cpue,Qs,Ktch,theta,HR_init,HR_init.sd,r.mean,r.sd)
{
  #Population dynamics
  K=exp(theta[1])
  r=exp(theta[2])
  if(estim.q=="YES")q=exp(theta[3:length(theta)])
  Bt=rep(NA,length(Ktch)+1)
  Bpen=Bt
  Bt[1]=K*Init.propK
  pen=posfun(Bt[1],max(Ktch,na.rm=T)*Init.propK,100)
  Bt[1]=pen[[1]]
  if(length(pen)>1)Bpen[1]=pen[[2]]
  for(t in 2:length(Bt))
  {
    Bt[t]=Bt[t-1]+r*Bt[t-1]*(1-Bt[t-1]/K)-Ktch[t-1]
    pen=posfun(Bt[t],1,10)
    Bt[t]=pen[[1]]
    if(Bt[t]<0) Bpen[t]=pen[[2]]
  }
  H=Ktch/Bt[-length(Bt)]
  
  #Loglikelihoods
    #Initial harvest rate
  HR.init.negLL=-log(dnorm(H[1],HR_init,HR_init.sd))
  
    #r
  r.negLL=-log(dnorm(r,r.mean,r.sd))
  names(r.negLL)=NULL
  
    #cpue likelihood
  ln.cpue.hat=vector('list',length(cpue))
  ln.cpue=ln.cpue.hat
  cpue.negLL=rep(NA,length(cpue))
  for(ku in 1:length(Qs))
  {
    if(!is.na(Qs[ku]))
    {
      no.cpue=which(!is.na(cpue[[ku]]))
      n.cpue=length(no.cpue)
      if(estim.q=="NO")q[match(names(Qs[ku]),names(q))]=exp(sum(log(cpue[no.cpue]/Bt[no.cpue]))/n.cpue) #Haddon 2001 page 321
      cpue.hat=q[match(names(Qs[ku]),names(q))]*Bt[-length(Bt)]
      
      ln.cpue.hat[[ku]]=log(cpue.hat)
      ln.cpue[[ku]]=log(cpue[[ku]])
      ln.cpue.hat[[ku]]=ln.cpue.hat[[ku]][no.cpue]
      ln.cpue[[ku]]=ln.cpue[[ku]][no.cpue]
      
      sqres=(ln.cpue[[ku]]-ln.cpue.hat[[ku]])^2
      cpue.negLL[ku]=(n.cpue/2)*(log(2*pi)+2*log(sqrt(sum(sqres,na.rm=T)/n.cpue))+1)
      
    }
  }
  cpue.negLL=sum(cpue.negLL,na.rm=T)
  
  #Total negloglike
  negLL=HR.init.negLL+r.negLL+cpue.negLL+sum(Bpen,na.rm=T)
  
  #Calculate MSY quantities
  Bmsy=K/2
  MSY=K*r/4
  Fmsy=r/2
  
  return(list(Bmsy=Bmsy,MSY=MSY,Fmsy=Fmsy,Bt=Bt,negLL=negLL,
              ln.cpue.hat=ln.cpue.hat,ln.cpue=ln.cpue))
}
  #plotting functions
fn.plt.cpue=function(ob,pred,Yr)
{
  plot(Yr,pred,pch=19,cex=1.5,type='b',ylab="lncpue",xlab="year",
       ylim=c(min(c(ob,pred)),max(c(ob,pred))))
  points(Yr,ob,col="orange",pch=19,cex=1.5)
  legend("bottomright",c("observed","predicted"),pch=19,cex=2,col=c("orange","black"),bty='n')
}
fn.plt.bio.ktch=function(Yr,Bt,Bmsy,Ktch)
{
  plot(Yr,Bt[1:length(Yr)],type='o',pch=19,ylim=c(0,max(Bt)),
       ylab="Total biomass (tonnes)",xlab="Year")
  abline(h=Bmsy,col="orange",lwd=2)              #threshold
  abline(h=FAO.lim.prop*Bmsy,col="red",lwd=2)    #limit         
  abline(h=FAO.tar.prop*Bmsy,col="green",lwd=2)  #target
  
  par(new=T)
  plot(Yr,Ktch,type='l',col="steelblue",xlab="",ylab="",axes=F,lwd=2)
  axis(side = 4)
  mtext(side = 4, line = 2, 'Total catch (tonnes')
}
fn.fill=function(x)
{
  aa=all.iers[which(!all.iers%in%x$yr)]
  aa1=x[1:length(aa),]
  aa1[,]=NA
  aa1$yr=aa
  x=rbind(x,aa1)%>%arrange(yr)
  return(x)
}

  #Loop over all species
Store.SPM=vector('list',N.sp)
names(Store.SPM)=Specs$SP.group
Store.stuff=Store.SPM
for(s in 1: N.sp)
{
  Id=match(Specs$SP.group[s],names(cpue.list))
  
  #catch
  ct=Tot.ktch%>%filter(SP.group==Specs$SP.group[s])%>%
                group_by(finyear)%>%
                summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
  all.iers=seq(min(ct$finyear),max(ct$finyear))
  
  if(!is.null(cpue.list[[Id]]))
  {
    #Get cpues
    CPUE.1=cpue.list[[Id]]$Nat
    if(!is.null(CPUE.1))CPUE.1=CPUE.1%>%mutate(CV=CV/100)%>%
                                        filter(CV<MAX.CV)
    CPUE.2=cpue.list[[Id]]$TDGDLF.mon
    if(!is.null(CPUE.2))CPUE.2=CPUE.2%>%
                                mutate(yr=as.numeric(substr(Finyear,1,4)),
                                       MeAn=Mean)%>%
                                filter(CV<MAX.CV)
    CPUE.3=cpue.list[[Id]]$TDGDLF.day
    if(!is.null(CPUE.3))CPUE.3=CPUE.3%>%
                                  mutate(yr=as.numeric(substr(Finyear,1,4)),
                                         MeAn=Mean)%>%
                                  filter(CV<MAX.CV)
    #fill in missing cpue
    CPUE=list(CPUE.1,CPUE.2,CPUE.3)
    CPUE.CV=CPUE.yr=CPUE
    for(ci in 1:length(CPUE))
    {
      if(!is.null(CPUE[[ci]]))
      {
        x=fn.fill(CPUE[[ci]])
        CPUE.CV[[ci]]=x$CV
        CPUE[[ci]]=x$MeAn
        CPUE.yr[[ci]]=x%>%filter(!is.na(MeAn))%>%select(yr)
      }
    }
 
    #input pars
      #r prior
    r.prior=store.species[[Id]]$r.prior.normal$mean
    r.prior.sd=store.species[[Id]]$r.prior.normal$sd
    
      #estimable pars
    Mx.ktch=max(ct$LIVEWT.c,na.rm=T)
    K.init=10*Mx.ktch
    r.init=Init.r[[Id]]
    QS=Estimable.qs[[Id]]
    
    #loop over scenarios
    dummy=vector('list',length(HR.o.scens))
    names(dummy)=HR.o.scens
    for(h in 1:length(HR.o.scens))
    {
      #Initial harvest rate prior
      HR_o=HR.o.scens[h]
      
      dummy.eff=vector('list',length(Efficien.scens))
      names(dummy.eff)=Efficien.scens
      Store.CPUE.eff=dummy.eff
      Eff.yrs=ct$finyear
      id.eff.yrs=which(Eff.yrs>1994)
      for(e in 1:length(Efficien.scens))
      {
        #Apply assumed efficiency to TDGDLF
        Add.eff=data.frame(yr=Eff.yrs,Efficiency=1)
        Add.eff$Efficiency[id.eff.yrs]=Add.eff$Efficiency[id.eff.yrs]-
          cumsum(rep(Efficien.scens[e],length(id.eff.yrs)))
         CPUE.eff.scen=CPUE
        for(ss in 2:3)
        {
          if(!is.null(CPUE.eff.scen[[ss]])) CPUE.eff.scen[[ss]]=CPUE.eff.scen[[ss]]*Add.eff$Efficiency
        }
         
        #objfun to minimize 
        fn_ob=function(theta)SPM(Init.propK=B.init,cpue=CPUE.eff.scen,Qs=QS,
                                 Ktch=ct$LIVEWT.c,theta,
                                 HR_init=HR_o,HR_init.sd=HR_o.sd,
                                 r.mean=r.prior,r.sd=r.prior.sd)$negLL
        #fit model
        if(estim.q=="YES")theta= c(k=log(K.init),r=log(r.init),log(QS[which(!is.na(QS))]))
        if(estim.q=="NO")theta= c(k=log(K.init),r=log(r.init))
        Lw.bound=log(c(2*Mx.ktch,0.001,rep(1e-6,length(theta)-2)))
        Up.bound=log(c(50*Mx.ktch,0.5,rep(1,length(theta)-2)))
        dummy.eff[[e]]=optim(theta,fn_ob,method="L-BFGS-B",lower =Lw.bound,
                         upper = Up.bound,hessian=T)
        Store.CPUE.eff[[e]]=CPUE.eff.scen
      }
      dummy[[h]]=dummy.eff
     }
    Store.SPM[[s]]=dummy
    Store.stuff[[s]]=list(cpue=Store.CPUE.eff,Qs=QS,Ktch=ct$LIVEWT.c,
                          r.mean=r.prior,r.sd=r.prior.sd,yrs=all.iers,
                          cpue.yrs=CPUE.yr,n.cpues=length(CPUE))
  }
}


  #Evaluate at MLE
SPM.preds=vector('list',length(Store.SPM))
names(SPM.preds)=names(Store.SPM)
for(s in 1: N.sp)
{
  if(!is.null(Store.SPM[[s]]))
  {
    dumy.pred=vector('list',length(HR.o.scens))
    names(dumy.pred)=HR.o.scens
    for(h in 1:length(HR.o.scens))
    {
      dummy.eff=vector('list',length(Efficien.scens))
      names(dummy.eff)=Efficien.scens
      for(e in 1:length(Efficien.scens))
      {
        dummy.eff[[e]]=SPM(Init.propK=B.init,
            cpue=Store.stuff[[s]]$cpue[[e]],
            Qs=Store.stuff[[s]]$Qs,
            Ktch=Store.stuff[[s]]$Ktch,
            theta=Store.SPM[[s]][[h]][[e]]$par,
            HR_init=HR.o.scens[h],
            HR_init.sd=HR_o.sd,
            r.mean=Store.stuff[[s]]$r.mean,
            r.sd=Store.stuff[[s]]$r.sd)
      }
      dumy.pred[[h]]=dummy.eff
    }
    SPM.preds[[s]]=dumy.pred
  }
}


#Plot obs VS pred cpues   #ACA show by species, h scenario, e scenario and each cpue   
for(s in 1: N.sp)
{
  if(!is.null(SPM.preds[[s]]))
  {
    for(h in 1:length(HR.o.scens))
    {
      for(e in 1:length(Efficien.scens))
      {
        for(x in 1:Store.stuff[[s]]$n.cpues)
          if(!is.null(SPM.preds[[s]][[h]][[e]]$ln.cpue[[x]]))
          {
            fn.plt.cpue(ob=SPM.preds[[s]][[h]][[e]]$ln.cpue[[x]],
                        pred=SPM.preds[[s]][[h]][[e]]$ln.cpue.hat[[x]],
                        Yr=c(Store.stuff[[s]]$cpue.yrs[[x]]$yr))
          }
        
      }
    }
      
  }
}


#Plot biomass
fn.plt.bio.ktch(Yr=1975:2017,Bt=Preds$Bt,Bmsy=Preds$Bmsy,KK)
fn.plt.bio.ktch(Yr=1975:2017,Bt=Preds.init$Bt,Bmsy=Preds.init$Bmsy,KK)



# Multi-species BDM -----------------------------------------------------------------------



#---RESULTS SECTION------

setwd(paste(hNdl,'/Outputs',sep=''))

#Plot spatial catch
data(worldLLhigh)
xlm=c(112,130)
ylm=c(-36,-10)
Map.this=subset(Tot.ktch,Name%in%Keep.species)
Map.sp=sort(unique(Map.this$SP.group))
fn.sptial.ktch=function(d,NMs)
{
  D=d
  D.agg=aggregate(LIVEWT.c~BLOCKX,D,sum)
  D.agg$LAT.cen=-(as.numeric(substr(D.agg$BLOCKX,1,2))+.5)
  D.agg$LONG.cen=100+as.numeric(substr(D.agg$BLOCKX,3,4))+.5
  scaler=max(D.agg$LIVEWT.c)/3.5
  
  plotMap(worldLLhigh, xlim=xlm,ylim=ylm,plt = c(.001, 1, 0.075, 1),
          col="grey90",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  points(D.agg$LONG.cen,D.agg$LAT.cen,cex=D.agg$LIVEWT.c/scaler,bg="grey50",pch=21)
  axis(side = 1, at =round(xlm[1]):xlm[2], labels = F, tcl = .25)
  axis(side = 2, at = round(ylm[1]):ylm[2], labels = F,tcl = .25)
  box()
  legend('topleft',NMs,bty='n',cex=1.1,xjust=0)
  Lg=round(quantile(D.agg$LIVEWT.c,probs=c(.75,.95,1)))
  legend('right',paste(Lg),pch=21,pt.bg="grey50",bty='n',pt.cex=Lg/scaler,title="Tonnes",cex=1.1)
}
fn.fig("Figure1_Map", 1600, 2400)
smart.par(n.plots=length(Map.sp),MAR=c(3,2,.1,.1),OMA=c(2.5,2,1.5,.1),MGP=c(1,.5,0))
for(s in 1: length(Map.sp))
{
  NMs=Map.sp[s]
  if(NMs=="Low") NMs="Low resilience"
  if(NMs=="Very.low") NMs="Very low resilience"
  
  fn.sptial.ktch(d=subset(Map.this,Name==Map.sp[s]),NMs=NMs)
  if(s%in%c(7,8,9))axis(side = 1, at =seq(xlm[1],xlm[2],4), labels = seq(xlm[1],xlm[2],4), tcl = .5,las=1,cex.axis=1)
  if(s%in%c(1,4,7))axis(side = 2, at = seq(ylm[1],ylm[2],4), labels = -seq(ylm[1],ylm[2],4),tcl = .5,las=2,cex.axis=1)
}
mtext(expression(paste("Latitude ",degree,"S")),side=2,line=0.6,las=3,cex=1.2,outer=T)
mtext(expression(paste("Longitude ",degree,"E")),side=1,line=0.65,cex=1.2,outer=T)
dev.off()




#Catch-MSY
if(Do.Ktch.MSY)
{
  
  #r priors   
  fn.fig("Figure1_Prior_r", 2000, 2000)
  smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(s in 1: N.sp)
  {
    NMs=names(store.species)[s]
    if(NMs=="Low") NMs="Low resilience"
    if(NMs=="Very.low") NMs="Very low resilience"
    plot(density(rgamma(10000, shape = store.species[[s]]$r.prior$shape, rate = store.species[[s]]$r.prior$rate)),
         lwd=3,main=NMs,xlab="",ylab="",cex.lab=2,cex.axis=1.15,col=1,xlim=c(0,.4),yaxt='n')
  }
  mtext(expression(paste(plain("Intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),1,0.5,cex=1.35,outer=T)
  mtext("Density",2,0,las=3,cex=1.35,outer=T)
  dev.off()
  
  YrS=sort(unique(Tot.ktch$finyear))
  
  #Catch
  fn.fig("Figure2", 2000, 2200)
  smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(s in 1: N.sp)
  {
    NMs=names(store.species)[s]
    if(NMs=="Low") NMs="Low resilience"
    if(NMs=="Very.low") NMs="Very low resilience"
    A=subset(Tot.ktch,SP.group==Specs$SP.group[s])
    ct.region=aggregate(LIVEWT.c~finyear+Region,A,sum)
    ct=aggregate(LIVEWT.c~finyear,A,sum)
    ct=ct[order(ct$finyear),]
    mis.yr=YrS[which(!YrS%in%ct$finyear)]
    add.yr=ct[mis.yr,]
    if(nrow(add.yr)>0)
    {
      add.yr$finyear=mis.yr
      add.yr$LIVEWT.c=0
      ct=rbind(ct,add.yr)
    }
    ct=ct[order(ct$finyear),]
    
    NORTH=subset(ct.region,Region=="North")
    mis.yr=ct$finyear[which(!ct$finyear%in%NORTH$finyear)]
    add.yr=NORTH[mis.yr,]
    if(nrow(add.yr)>0)
    {
      add.yr$finyear=mis.yr
      add.yr$LIVEWT.c=0
      NORTH=rbind(NORTH,add.yr)
    }
    
    NORTH=NORTH[order(NORTH$finyear),]
    
    SOUTH=subset(ct.region,Region=="South")
    mis.yr=ct$finyear[which(!ct$finyear%in%SOUTH$finyear)]
    add.yr=SOUTH[mis.yr,]
    if(nrow(add.yr)>0)
    {
      add.yr$finyear=mis.yr
      add.yr$LIVEWT.c=0
      SOUTH=rbind(SOUTH,add.yr)
    }
    
    SOUTH=SOUTH[order(SOUTH$finyear),]
    
    plot(ct$finyear, ct$LIVEWT.c, type="l", ylim = c(0, 1.01*max(ct$LIVEWT.c)), xlab = "", ylab = "", main=NMs,
         lwd=2,cex.lab=2,cex.axis=1.15,cex.main=1.3)
    with(NORTH,lines(finyear,LIVEWT.c,col='grey40',lty=3,lwd=2))
    with(SOUTH,lines(finyear,LIVEWT.c,col='grey70',lty=5,lwd=2))
    
    if(s==6)legend("topleft",c("Total catch","North","South"),bty='n',
                   col=c("black","grey40","grey70"),lty=c(1,3,5),lwd=2,cex=1.25)
    
  }
  mtext("Financial year",1,0.25,cex=1.35,outer=T)
  mtext("Total catch (tonnes)",2,0.35,las=3,cex=1.35,outer=T)
  dev.off()
  
  #Relative biomass
  CL="grey55"
  CL.mean="transparent"
  
  Low.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (0+Nper)/100))   #get percentiles
  High.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (100-Nper)/100))
  fn.cons.po=function(low,up) c(low, tail(up, 1), rev(up), low[1])  #construct polygon
  
  do.col="NO"
  if(do.col=="NO") colfunc <- colorRampPalette(c("grey95","grey60"))
  if(do.col=="YES") colfunc <- colorRampPalette(c("aliceblue","lightblue3"))
  if(do.col=="NO") CL.ref.pt=c("black","grey30","grey75","grey90")
  if(do.col=="YES") CL.ref.pt=c('forestgreen','yellow','orange','red')
  
  COLS=colfunc(3)
  fn.plot.percentile=function(DAT,YR,ADD.prob,add.RP.txt,CEX,CX.AX)
  {
    #50% of data
    Nper=(100-50)/2
    LOW.50=Low.percentile(Nper,DAT)
    UP.50=High.percentile(Nper,DAT)
    
    #75% of data
    Nper=(100-75)/2
    LOW.75=Low.percentile(Nper,DAT)
    UP.75=High.percentile(Nper,DAT)
    
    #100% of data
    Nper=(100-100)/2
    LOW.100=Low.percentile(Nper,DAT)
    UP.100=High.percentile(Nper,DAT)
    
    #construct polygons
    Year.Vec <-  fn.cons.po(YR,YR)
    Biom.Vec.50 <- fn.cons.po(LOW.50,UP.50) 
    Biom.Vec.75 <- fn.cons.po(LOW.75,UP.75) 
    Biom.Vec.100 <-fn.cons.po(LOW.100,UP.100) 
    
    
    #plot
    plot(YR,UP.100,ylim=c(0,max(UP.100)),type="l",ylab="",xlab="",xaxt='n',col='transparent',cex.axis=CX.AX)
    
    polygon(Year.Vec, Biom.Vec.100, col = COLS[3], border = "grey20")
    polygon(Year.Vec, Biom.Vec.75, col = COLS[2], border = "grey20")
    polygon(Year.Vec, Biom.Vec.50, col = COLS[1], border = "grey20")
    
    
    #add probs
    if(ADD.prob=="YES")
    {
      add.probs(id.yr=match(Current,YR),YR,DAT,UP.100,LOW.100,SRT=0,CEX)
      
      abline(h=B.target,lwd=1.5,col='grey45',lty=3)
      abline(h=B.threshold,lwd=1.5,col='grey45',lty=3)
      abline(h=B.limit,lwd=1.5,col='grey45',lty=3)
      
      if(add.RP.txt=="YES")
      {
        text(YR[4],B.target,"Target",pos=3,cex=1.1)
        text(YR[4],B.threshold,"Threshold",pos=3,cex=1.1)
        text(YR[4],B.limit,"Limit",pos=3,cex=1.1)
      }
      
    }
    axis(1,at=YR,labels=F,tck=-0.015)
    axis(1,at=seq(YR[1],YR[length(YR)],5),labels=seq(YR[1],YR[length(YR)],5),tck=-0.03,cex.axis=CX.AX)
  }
  add.probs=function(id.yr,YR,DAT,UP.100,LOW.100,SRT,CEX)
  {
    f=ecdf(DAT[id.yr,])
    P.below.target=f(B.target)
    P.below.threshold=f(B.threshold)
    P.below.limit=f(B.limit)
    P.above.target=1-P.below.target
    P.above.threshold=1-P.below.threshold
    P.above.limit=1-P.below.limit
    P.between.thre.tar=P.below.target-P.below.threshold
    P.between.lim.thre=P.below.threshold-P.below.limit
    if(P.above.target>0)
    {
      segments(YR[id.yr],B.target,YR[id.yr],UP.100[id.yr],col=CL.ref.pt[1],lwd=8,lend="butt")  
      text(YR[id.yr],B.target*1.5,paste(round(100*P.above.target,1),"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
    if(P.between.thre.tar>0)
    {
      segments(YR[id.yr],B.target,YR[id.yr],B.threshold,col=CL.ref.pt[2],lwd=8,lend="butt")  
      text(YR[id.yr],mean(c(B.target,B.threshold))*1.025,paste(round(100*P.between.thre.tar,1),"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
    if(P.between.lim.thre>0)
    {
      segments(YR[id.yr],B.threshold,YR[id.yr],B.limit,col=CL.ref.pt[3],lwd=8,lend="butt")
      text(YR[id.yr],mean(c(B.threshold,B.limit))*1.025,paste(round(100*P.between.lim.thre,1),"%",sep=""),
           col="black",cex=CEX,srt=SRT,font=2,pos=2)
    }
    if(P.below.limit>0)
    {
      segments(YR[id.yr],B.limit,YR[id.yr],LOW.100[id.yr],col=CL.ref.pt[4],lwd=8,lend="butt")
      text(YR[id.yr],B.limit*0.8,paste(round(100*P.below.limit,1),"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
    if(P.below.limit==0)
    {
      segments(YR[id.yr],B.limit,YR[id.yr],LOW.100[id.yr],col=CL.ref.pt[4],lwd=8,lend="butt")
      text(YR[id.yr],B.limit*0.8,paste(0,"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
  }
  
  fn.fig("Figure3_Biomass_relative_Base case", 2000, 2200)
  smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(s in 1: N.sp)
  {
    Yrs=store.species[[s]]$Catch$finyear
    Ktch_MSY_Rel.bio=store.species[[s]]$KTCH.MSY$BaseCase$bt.rel
    
    #Percentile   
    fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES",add.RP.txt="NO",CEX=1.2,CX.AX=1.25)
    if(s==1) legend("bottomleft",c("50%","75%","100%"),fill=COLS,bty='n',cex=1,horiz=T)
    NMs=names(store.species)[s]
    if(NMs=="Low") NMs="Low resilience"
    if(NMs=="Very.low") NMs="Very low resilience"
    mtext(NMs,3,0)
    #Geometric mean
    # Ktch_MSY_rel_bt_mean=Ktch_MSY_rel_bt_lowSE=Ktch_MSY_rel_bt_upSE=nrow(Ktch_MSY_Rel.bio)
    # for(nr in 1:nrow(Ktch_MSY_Rel.bio))
    # {
    #   Ktch_MSY_rel_bt_mean[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])))
    #   Ktch_MSY_rel_bt_upSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) + 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
    #   Ktch_MSY_rel_bt_lowSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) - 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
    # }
    # plot(Yrs,Ktch_MSY_rel_bt_mean,cex=.95,pch=19,col=CL.mean,ylim=c(0,1),xaxt='n',xlab="",ylab="",
    #      main=names(store.species)[s],cex.axis=1.15,cex.main=1.3)
    # segments(Yrs,Ktch_MSY_rel_bt_lowSE,Yrs,Ktch_MSY_rel_bt_upSE,col=CL)
    # abline(h=B.target,lwd=1,col='black',lty=2)
    # #text(Yrs[3],B.target,"Target",pos=3,cex=1.25)
    # abline(h=B.threshold,lwd=1,col='grey30',lty=2)
    # #text(Yrs[3],B.threshold,"Threshold",pos=3,cex=1.25)
    # abline(h=B.limit,lwd=1,col='grey50',lty=2)
    # #text(Yrs[3],B.limit,"Limit",pos=3,cex=1.25)
    # axis(1,Yrs,labels=F,tck=-0.015)
    # axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
  }
  mtext("Financial year",1,0.25,cex=1.35,outer=T)
  mtext("Relative biomass",2,0.25,las=3,cex=1.35,outer=T)
  dev.off()
  
  
  
  #Current depletion of total biomass   
  fn.fig("Posterior_current_year_depletion_Base case", 2000, 2200)
  smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(s in 1: N.sp)
  {
    Yrs=store.species[[s]]$Catch$finyear
    Current=Yrs[length(Yrs)]
    Ktch_MSY_Rel.bio=store.species[[s]]$KTCH.MSY$BaseCase$bt.rel
    Ktch_MSY_current_yr=Ktch_MSY_Rel.bio[length(Yrs),]
    NMs=names(store.species)[s]
    if(NMs=="Very.low") NMs="Very low"
    density.fun2(what=Ktch_MSY_current_yr,MAIN=NMs)
  }
  mtext(paste(Current,"Relative biomass"),1,line=0.25,cex=1.35,outer=T)
  mtext("Probability",2,0.25,las=3,cex=1.35,outer=T)
  dev.off()
  
  
  #Sensitivity tests Biomass
  fn.plot.percentile.sens=function(DAT,YR,ADD.prob,add.RP.txt,CEX,CX.AX,AdYXs,AdXXs)
  {
    #50% of data
    Nper=(100-50)/2
    LOW.50=Low.percentile(Nper,DAT)
    UP.50=High.percentile(Nper,DAT)
    
    #75% of data
    Nper=(100-75)/2
    LOW.75=Low.percentile(Nper,DAT)
    UP.75=High.percentile(Nper,DAT)
    
    #100% of data
    Nper=(100-100)/2
    LOW.100=Low.percentile(Nper,DAT)
    UP.100=High.percentile(Nper,DAT)
    
    #construct polygons
    Year.Vec <-  fn.cons.po(YR,YR)
    Biom.Vec.50 <- fn.cons.po(LOW.50,UP.50) 
    Biom.Vec.75 <- fn.cons.po(LOW.75,UP.75) 
    Biom.Vec.100 <-fn.cons.po(LOW.100,UP.100) 
    
    
    #plot
    plot(YR,UP.100,ylim=c(0,max(UP.100)),type="l",ylab="",xlab="",yaxt='n',xaxt='n',col='transparent',cex.axis=CX.AX)
    
    polygon(Year.Vec, Biom.Vec.100, col = COLS[3], border = "grey20")
    polygon(Year.Vec, Biom.Vec.75, col = COLS[2], border = "grey20")
    polygon(Year.Vec, Biom.Vec.50, col = COLS[1], border = "grey20")
    
    
    #add probs
    if(ADD.prob=="YES")
    {
      add.probs(id.yr=match(Current,YR),YR,DAT,UP.100,LOW.100,SRT=0,CEX)
      
      abline(h=B.target,lwd=1.5,col='grey45',lty=3)
      abline(h=B.threshold,lwd=1.5,col='grey45',lty=3)
      abline(h=B.limit,lwd=1.5,col='grey45',lty=3)
      
      if(add.RP.txt=="YES")
      {
        text(YR[4],B.target,"Target",pos=3,cex=1.1)
        text(YR[4],B.threshold,"Threshold",pos=3,cex=1.1)
        text(YR[4],B.limit,"Limit",pos=3,cex=1.1)
      }
      
    }
    if(AdYXs=="YES")axis(2,at=seq(0,1,.2),labels=seq(0,1,.2),tck=-0.05,las=1,cex.axis=CX.AX)
    axis(1,at=YR,labels=F,tck=-0.025)
    axis(1,at=seq(YR[1],YR[length(YR)],5),labels=F,tck=-0.05,cex.axis=CX.AX)
    if(AdXXs=="YES")axis(1,at=seq(YR[1],YR[length(YR)],5),labels=seq(YR[1],YR[length(YR)],5),tck=-0.05,cex.axis=CX.AX)
  }
  fn.fig("Sensitivity_Biomass_relative", 2800,800)
  par(mfcol=c((Nscen-1),N.sp),mar=c(1,1,.15,.15),oma=c(2,3.5,.95,.25),las=1,mgp=c(1,.5,0))
  for(s in 1: N.sp)
  {
    Yrs=store.species[[s]]$Catch$finyear
    NMs=names(store.species)[s]
    if(NMs=="Scalloped hammerhead") NMs="Scalloped hh"
    if(NMs=="Smooth hammerhead") NMs="Smooth hh"
    for(sc in 2:Nscen)
    {
      Ktch_MSY_Rel.bio=store.species[[s]]$KTCH.MSY[[sc]]$bt.rel
      SCNE.nm=names(store.species[[s]]$KTCH.MSY)[sc]
      AddY="NO"
      if(s==1) AddY="YES"
      AddX='NO'
      if(sc==3) AddX="YES"
      fn.plot.percentile.sens(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES",
                              add.RP.txt="NO",CEX=.7,CX.AX=.8,AdYXs=AddY,AdXXs=AddX)
      
      
      #if(s==N.sp & sc==2) legend("bottom",c("50%","75%","100%"),fill=COLS,bty='n',cex=.65,horiz=T)
      if(s==1) mtext(SCNE.nm,2,1.5)
      if(NMs=="Low") NMs="Low resilience"
      if(NMs=="Very.low") NMs="Very low resilience"
      if(sc==2)mtext(NMs,3,0,cex=.75)
    }
  }
  mtext("Financial year",1,1,cex=1.35,outer=T)
  mtext("Relative biomass",2,2,las=3,cex=1.35,outer=T)
  dev.off()
  
  
  #Fishing mortality
  fn.fig("Fishing_mortality_Base case", 2000, 2000)
  smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(s in 1: N.sp)
  {
    Yrs=store.species[[s]]$Catch$finyear
    Ktch_MSY_Rel.bio=store.species[[s]]$KTCH.MSY$BaseCase$Fish.mort
    
    #Percentile
    fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="NO",add.RP.txt="NO",CEX=1.2,CX.AX=1.2)
    if(s==1) legend("topleft",c("50%","75%","100%"),fill=COLS,bty='n',cex=1.25)
    NMs=names(store.species)[s]
    if(NMs=="Very.low") NMs="Very low"
    mtext(NMs,3,0)
    
    
    #   #Geometric mean
    # Ktch_MSY_rel_bt_mean=Ktch_MSY_rel_bt_lowSE=Ktch_MSY_rel_bt_upSE=nrow(Ktch_MSY_Rel.bio)
    # for(nr in 1:nrow(Ktch_MSY_Rel.bio))
    # {
    #   Ktch_MSY_rel_bt_mean[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])))
    #   Ktch_MSY_rel_bt_upSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) + 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
    #   Ktch_MSY_rel_bt_lowSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) - 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
    # }
    # plot(Yrs,Ktch_MSY_rel_bt_mean,cex=1.1,pch=19,col=CL.mean,ylim=c(0,max(Ktch_MSY_rel_bt_upSE,na.rm=T)),xaxt='n',xlab="",ylab="",
    #      main=names(store.species)[s],cex.axis=1.15)
    # segments(Yrs,Ktch_MSY_rel_bt_lowSE,Yrs,Ktch_MSY_rel_bt_upSE,col=CL)
    # axis(1,Yrs,labels=F,tck=-0.015)
    # axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
  }
  mtext("Financial year",1,0.5,cex=1.35,outer=T)
  mtext(expression(paste(plain("Fishing mortality (year") ^ plain("-1"),")",sep="")),2,0,las=3,cex=1.35,outer=T)
  dev.off()
  
  
  #Catch and MSY    
  fn.fig("Figure4_CatchMSY_Plots", 2000, 2200)
  smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(s in 1: N.sp)
  {
    yr=store.species[[s]]$Catch$finyear
    ct=store.species[[s]]$Catch$Total.ktch
    mean_ln_msy=as.numeric(store.species[[s]]$KTCH.MSY$BaseCase$"geom. mean MSY (tons)")
    msy=store.species[[s]]$KTCH.MSY$BaseCase$msy
    Mean.MSY=exp(mean(log(msy)))
    Mean.MSY_UP=exp(mean(log(msy)) + 1.96 * sd(log(msy)))  
    Mean.MSY_LOW=exp(mean(log(msy)) - 1.96 * sd(log(msy)))
    NMs=names(store.species)[s]
    if(NMs=="Low") NMs="Low resilience complex"
    if(NMs=="Very.low") NMs="Very low resilience complex"
    
    plot(yr, ct, type="l", ylim = c(0, 1.01*max(c(ct,Mean.MSY_UP))), xlab = "", ylab = "", 
         lwd=2,cex.lab=2,cex.axis=1.15)
    mtext(NMs,3,0)
    all.yrs=c(yr[1]-2,yr,yr[length(yr)]+2)
    polygon(c(all.yrs,rev(all.yrs)),  c(rep(Mean.MSY_LOW,length(all.yrs)),rep(Mean.MSY_UP,length(all.yrs))),
            col=rgb(.1,.1,.1,alpha=.15),border="transparent")  
    abline(h=Mean.MSY,col="grey50", lwd=2.5,lty=3)
    if(s==9)legend("topright",c("MSY (?1.96 SE)"),bty='n',col=c("grey50"),lty=3,lwd=2.5,cex=1.25)
    
  }
  mtext("Financial year",1,0.25,cex=1.35,outer=T)
  mtext("Total catch (tonnes)",2,0.35,las=3,cex=1.35,outer=T)
  dev.off()
  
  
  #Sensitivity tests MSY
  Exprt.MSY=vector('list',length(N.sp))
  fn.fig("Sensitivity_MSY", 1000, 2400)
  par(mfrow=c(N.sp,Nscen),mar=c(2,2,.25,.35),oma=c(1.75,2.75,.95,.25),las=1,mgp=c(1,.5,0))
  for(s in 1: N.sp)
  {
    Yrs=store.species[[s]]$Catch$finyear
    NMs=names(store.species)[s]
    yr=store.species[[s]]$Catch$finyear
    ct=store.species[[s]]$Catch$Total.ktch
    dummyMSY=vector('list',length(Nscen))
    for(sc in 1:Nscen)
    {
      mean_ln_msy=as.numeric(store.species[[s]]$KTCH.MSY[[sc]]$"geom. mean MSY (tons)")
      msy=store.species[[s]]$KTCH.MSY[[sc]]$msy
      Mean.MSY=exp(mean(log(msy)))
      Mean.MSY_UP=exp(mean(log(msy)) + 1.96 * sd(log(msy)))  
      Mean.MSY_LOW=exp(mean(log(msy)) - 1.96 * sd(log(msy)))
      
      plot(yr, ct, type="l", ylim = c(0, 1.01*max(c(ct,Mean.MSY_UP))), xlab = "", ylab = "",
           lwd=2,cex.lab=2,cex.axis=1)
      all.yrs=c(yr[1]-2,yr,yr[length(yr)]+2)
      polygon(c(all.yrs,rev(all.yrs)),  c(rep(Mean.MSY_LOW,length(all.yrs)),rep(Mean.MSY_UP,length(all.yrs))),
              col=rgb(.1,.1,.1,alpha=.15),border="transparent")  
      abline(h=Mean.MSY,col="grey50", lwd=2.5,lty=3)
      #if(s==4)legend("topright",c("MSY (?1.96 SE)"),bty='n',col=c("grey50"),lty=3,lwd=2.5,cex=1.25)
      SCNE.nm=names(store.species[[s]]$KTCH.MSY)[sc]
      if(s==1) mtext(SCNE.nm,3,0)
      if(NMs=="Low") NMs="Low resilience"
      if(NMs=="Very.low") NMs="Very low resilience"
      if(NMs=="Scalloped hammerhead") NMs="Scalloped hh"
      if(NMs=="Smooth hammerhead") NMs="Smooth hh"
      
      if(sc==1)mtext(NMs,2,2,las=3,cex=.8)
      dummyMSY[[sc]]=data.frame(Species=NMs,Scenario=SCNE.nm,LOW_CI=Mean.MSY_LOW,Mean.MSY=Mean.MSY,Mean.UP_CI=Mean.MSY_UP)
    }
    
    Exprt.MSY[[s]]=do.call(rbind,dummyMSY)
  }
  mtext("Financial year",1,0.25,cex=1.35,outer=T)
  mtext("Total catch (tonnes)",2,1.25,las=3,cex=1.35,outer=T)
  dev.off()
  write.csv(do.call(rbind,Exprt.MSY),"MSY_estimates.csv",row.names=F)
  
}








