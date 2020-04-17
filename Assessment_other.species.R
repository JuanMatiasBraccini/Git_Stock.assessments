# ------ Script for running SPM and catch-MSY stock assessment on other shark species---- ###################

#note: SPM and Catch-MSY implementation of Martell & Froese 2012.
#       total catches must be used (i.e. all sources of F)

#       All sources of mortalities (commercial and recreational catches) considered thru catch reconstructions

#       If catches have never been >1% carrying capacity, then unexploited status so catch series have
#       no information on productivity.

#missing:
# Include ALL species in final risk scoring
# Need CPUE for Copper sharks (run SPM and aSPM...); review Smooth HH cpue and mako cpue...; SPM Tiger fit
# Milk shark SPM, hitting upper K boundary, no trend in cpue, crap Hessian, too uncertain....mention in text...
# aSPM: finish running for all species; issues with Tiger cpue fit...


#consider:
# Rather than standard SPM, try JABBA: Just Another Bayesian Biomass Assessment (can be 
#   run from R..see Winker et al 2018; it's what IUCN uses)


rm(list=ls(all=TRUE))
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")
source.hnld="C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/"
fn.source=function(script)source(paste(source.hnld,script,sep=""))
fn.source("fn.fig.R")
fn.source("Leslie.matrix.R") 
fn.source("Catch_MSY.R")
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
Do.jpeg="NO" 
Do.tiff="YES"
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R")

library(MASS)
library(plotrix)
library(PBSmapping)
library(tidyverse)
library(mvtnorm)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("Biobase")
library(Biobase)
library(numDeriv)
library(spatstat.utils)
library(Hmisc)
library(ggplot2)
library(ggrepel)
library(datalowSA)

Asses.year=2020    #enter year of assessment
Last.yr.ktch="2017-18"

hNdl=paste("C:/Matias/Analyses/Population dynamics/1.Other species/",Asses.year,sep="")
fnkr8t=function(x) if(!dir.exists(x))dir.create(x)
fnkr8t(hNdl)
fnkr8t(paste(hNdl,"Outputs",sep="/"))



#---DATA SECTION-----
setwd("C:/Matias/Analyses/Data_outs")

#Total effort
Effort.monthly=read.csv("Annual.total.eff.days.csv",stringsAsFactors=F)
Effort.monthly.north=read.csv("Annual.total.eff_NSF.csv",stringsAsFactors=F)

Effort.monthly_blocks=read.csv("Effort.monthly.csv",stringsAsFactors=F)
Effort.daily_blocks=read.csv("Effort.daily.csv",stringsAsFactors=F)
Effort.monthly.north_blocks=read.csv("Effort.monthly.NSF.csv",stringsAsFactors=F)
Effort.daily.north_blocks=read.csv("Effort.daily.NSF.csv",stringsAsFactors=F)


#Total catch
fn.in=function(NM)
{
  read.csv(paste('C:/Matias/Analyses/Data_outs/',NM,sep=""),stringsAsFactors = F)
}

  #1.1 Catch_WA Fisheries

#Historic
Hist.expnd=fn.in(NM='recons_Hist.expnd.csv')

#Ammended reported catch including discards
Data.monthly=fn.in(NM='recons_Data.monthly.csv')
Data.monthly.north=fn.in(NM='recons_Data.monthly.north.csv')

#TEPS
Greynurse.ktch=fn.in(NM='recons_Greynurse.ktch.csv')
TEPS_dusky=fn.in(NM='recons_TEPS_dusky.csv')

WRL.ktch=fn.in(NM='Wetline_rocklobster.csv')


  #1.2. Catch of non WA Fisheries

#Taiwanese gillnet and longline
Taiwan.gillnet.ktch=fn.in(NM='recons_Taiwan.gillnet.ktch.csv')
Taiwan.longline.ktch=fn.in(NM='recons_Taiwan.longline.ktch.csv')


#Indonesian illegal fishing in Australia waters
Indo_total.annual.ktch=fn.in(NM='recons_Indo.IUU.csv') 

#AFMA's GAB & SBT fisheries
GAB.trawl_catch=fn.in(NM='recons_GAB.trawl_catch.csv') 
WTBF_catch=fn.in(NM='recons_WTBF_catch.csv') 

#SA Marine Scalefish fishery
Whaler_SA=fn.in(NM='recons_Whaler_SA.csv') 


  #2. WA Recreational catch
Rec.ktch=fn.in(NM='recons_recreational.csv')

#species codes
All.species.names=read.csv("C:/Matias/Data/Species_names_shark.only.csv") #for catch
#b=read.csv("C:\\Matias\\Data\\Species.code.csv")


#list of life history param for demography
LH.par=read.csv("C:/Matias/Data/Life history parameters/Life_History_other_sharks.csv",stringsAsFactors=F)

#PSA scores
PSA.list=read.csv('C:/Matias/Analyses/Population dynamics/PSA/PSA_scores_other.species.csv',stringsAsFactors=F)

#Temperature
#TEMP=read.csv("C:/Matias/Data/Oceanography/SST.csv")

#Species scientific names for assessed species
Shark.species=5001:24900
School.shark= 17008
Indicator.species=c(17001,17003,18003,18007)
Shar_other=22999

Scien.nm=All.species.names[,c('SPECIES','Scien.nm')]

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

#Abundance data
fn.read=function(x) read.csv(paste('C:/Matias/Analyses/Data_outs',x,sep='/'),stringsAsFactors = F)
  #Naturalist abundance survey
Scal.hh.nat=fn.read('Scalloped hammerhead.Srvy.FixSt.csv')
Tiger.nat=fn.read('Tiger shark.Srvy.FixSt.csv')
Mil.nat=fn.read('Milk shark.Srvy.FixSt.csv')

  #Standardised TDGDLF cpue
Smuz.hh.tdgdlf_mon=fn.read('Hammerhead.annual.abundance.basecase.monthly_relative.csv') #assumed to be all smooth HH
Smuz.hh.tdgdlf_daily=fn.read('Hammerhead.annual.abundance.basecase.daily_relative.csv')
Spinr.tdgdlf_mon=fn.read('Spinner Shark.annual.abundance.basecase.monthly_relative.csv')
Spinr.tdgdlf_daily=fn.read('Spinner Shark.annual.abundance.basecase.daily_relative.csv')
Tiger.tdgdlf_mon=fn.read('Tiger Shark.annual.abundance.basecase.monthly_relative.csv')
Tiger.tdgdlf_daily=fn.read('Tiger Shark.annual.abundance.basecase.daily_relative.csv')
Greynurse.tdgdlf_mon=fn.read('Greynurse Shark.annual.abundance.basecase.monthly._relative.csv')
Pencil.tdgdlf_mon=fn.read('Pencil Shark.annual.abundance.basecase.monthly_relative.csv')
Pencil.tdgdlf_daily=fn.read('Pencil Shark.annual.abundance.basecase.daily_relative.csv')
Sawshrk.tdgdlf_daily=fn.read('Sawsharks.annual.abundance.basecase.daily_relative.csv')
Mako.tdgdlf_mon=fn.read('Mako.annual.abundance.basecase.monthly_relative.csv')
Mako.tdgdlf_daily=fn.read('Mako.annual.abundance.basecase.daily_relative.csv')

  #Standardised NSF cpue
Lemon.NSF=fn.read('Lemon shark.annual.abundance.NSF_relative.csv')
Pigeye.NSF=fn.read('Pigeye shark.annual.abundance.NSF_relative.csv')
Tiger.NSF=fn.read('Tiger shark.annual.abundance.NSF_relative.csv')

   

#Mean catch weight data

  #Standardised TDGDLF mean size
Smuz.hh.tdgdlf.size=fn.read('Smooth hammerhead.annual.mean.size_relative.csv')
Spinr.tdgdlf.size=fn.read('Spinner Shark.annual.mean.size_relative.csv')
Tiger.tdgdlf.size=fn.read('Tiger Shark.annual.mean.size_relative.csv')



#Conventional tagging data
Tag=fn.read('Tagging_conventional.data.csv')



#---PARAMETERS SECTION-----
Explor="NO"
Asses.Scalloped.HH=FALSE  #2020 scalloped HH assessment


#Criteria for selecting species for quantitative analyses
Min.yrs=5
Min.ktch=5000 #in kg


#Reference points
B.threshold=0.5  #Bmys
Tar.prop=1.3    #target and limit proportions of Bmsy. source: Haddon et al 2014. Technical
Lim.prop=0.5    #   Reviews of Formal Harvest Strategies.
B.target=Tar.prop*B.threshold
B.limit=Lim.prop*B.threshold

  #Empirical reference points
Fmsy.emp=function(M) 0.41*M     #Zhou et al 2012
SPR.thre=0.3   #Punt 2000 Extinction of marine renewable resources: a demographic analysis. 
SPR.tar=0.4                # Population Ecology 42, 19–27


#Life history parameters for selected species  
pup.sx.ratio=.5
Mn.conv.Fl.Tl=.85  #average convertion (over gummy, dusky, whiskery, sandbar) from FL to TL


#Gillnet selectivity from different nets
Estim.sel.exp='NO'  #not enough observations from different mesh sizes


#Use conventional tagging data?
use.tags=F      #too few recaptures...do not use 


#... Surplus production arguments

#note: Only fitting species with species-specific abundance time series 
#     (e.g. wobbegongs comprise several species so not fitted)
#     Assumption, negligible exploitation at start of time series

    #Initial harvest rate 
HR_o.sd=0.005  #SD of HR likelihood (fixed)

    #Efficiency increase scenarios from 1995 on (done up to 1994 in cpue stand.)
Efficien.scens=c(0)
#Efficien.scens=c(.01)

    #Proportional biomass (as proportion of K) at start of catch time series
B.init=1 #(fixed)  Starting @ virgin level

    #Estimate q
estim.q="NO"   #use Haddon's q MLE calculation

    #cpue likelihood
what.like='kernel'
#what.like='full'


    #Initial estimated par value
Init.r=list("copper shark"=.15,"great hammerhead"=.15,
            "grey nurse shark"=0.1,
            "lemon shark"=.15,"milk shark"=.2,"sawsharks"=.15,
            "scalloped hammerhead"=.15,"smooth hammerhead"=.2,
            "spinner shark"=.15,"shortfin mako"=.1,"spurdogs"=.15,
            "tiger shark"=.15,"wobbegongs"=.15)

N.monte=1000

MAX.CV=0.5    #maximum acceptable CV


    #Define which optimisation method to use
#minimizer='nlminb'
minimizer='optim'

Remove.bounds=FALSE

usePen=TRUE  

    #K bounds
Low.bound.K=10  #times the maximum catch
Up.bound.K=100 

    #K init times max ktch
k.times.mx.ktch=mean(c(Low.bound.K,Up.bound.K))

    #fix or estimate r
fix.r="NO"
r.weight=1   #weight given in the likelihood function

#... Catch-MSY arguments

    #simulatins
SIMS=5e4  

    #Assumed process error
ERROR=0   #is default. 


    #depletion level at start of catch series
STARTBIO=c(B.init*.95,B.init)   #low depletion because starting time series prior to any fishing
FINALBIO=c(.2,.9)       #very uncertain


    #r priors
NsimSS=1000
r.prior="USER"  #demography
r.prior2=NA    #uniform

    #Scenarios  
Modl.rn="standard"   #for annual assessments
#Modl.rn='first'       #for paper
if(Modl.rn=="first")
{
  Nscen=2
  SCENARIOS=vector('list',Nscen)
  names(SCENARIOS)=c('BaseCase','WorstCase')
  SCENARIOS$BaseCase=list(Error=ERROR,R.prior=r.prior,Initial.dep=STARTBIO)
  SCENARIOS$WorstCase=list(Error=ERROR,R.prior=r.prior2,Initial.dep=STARTBIO)
}else
{
  Nscen=1
  SCENARIOS=vector('list',Nscen)
  names(SCENARIOS)=c('BaseCase')
  SCENARIOS$BaseCase=list(Error=ERROR,R.prior=r.prior,Initial.dep=STARTBIO)
}  


    #Future projections
years.futures=5

    #simulation test Catch-MSY for small and large catches
Do.sim.test="NO"


#... Control which assessment methods to implement
do.length.based='NO'        #not used due to knife-edge sel. assumption 
do.mean.weight.based="NO"   #not used due to knife-edge sel. assumption 
do.length.based.SPR="NO"   #not used due to knife-edge sel. assumption 
Do.SPM="YES"
Do.Ktch.MSY="YES"
Do.aSPM="NO"
do.Gedamke_Hoenig="NO"     #not used due to knife-edge sel. assumption 
                           # and available data in weights, not length

use.size.comp=T   #Use size composition?
Min.obs=100  #minimum number of size observations to derive selectivity curve
Min.len=25  #minimum length of sharks
Min.size.sample=150  #minimum number of observation to do size assessment

#---PROCEDURE SECTION-----

#Relevant catch years
Relevant.yrs=paste(seq(as.numeric(substr(min(Hist.expnd$FINYEAR),1,4)),
                       as.numeric(substr(Last.yr.ktch,1,4))),
                   substr(seq(as.numeric(substr(min(Hist.expnd$FINYEAR),1,4))+1,
                       as.numeric(substr(Last.yr.ktch,1,4))+1),3,4),sep='-')

#Select relevant effort vars  
Effort.monthly_blocks=Effort.monthly_blocks%>%
        filter(FINYEAR%in%Relevant.yrs)%>%
        count(FINYEAR,BLOCKX)%>%
        group_by(FINYEAR,BLOCKX)%>%
        mutate(n=ifelse(n>0,1,0))
Effort.daily_blocks=Effort.daily_blocks%>%
        rename(FINYEAR=finyear,
               BLOCKX=blockx)%>%
        filter(FINYEAR%in%Relevant.yrs)%>%
        count(FINYEAR,BLOCKX)%>%
        group_by(FINYEAR,BLOCKX)%>%
        mutate(n=ifelse(n>0,1,0))
Effort_blocks=rbind(Effort.monthly_blocks,Effort.daily_blocks)%>%
        count(FINYEAR,BLOCKX)%>%
        group_by(FINYEAR,BLOCKX)%>%
        mutate(n=ifelse(n>0,1,0))%>%
        group_by(FINYEAR)%>%
        summarise(Tot=sum(n))%>%
        data.frame

#get GN equivalent of LL effort
do.GN.quiv='NO'
if(do.GN.quiv=="YES")
{
  GN.ktch.north=Data.monthly.north%>%   
    filter(METHOD=="GN")%>%
    group_by(FINYEAR,BLOCKX)%>%
    summarise(Ktch.GN=sum(LIVEWT.c,na.rm=T))%>%
    mutate(dummy=paste(FINYEAR,BLOCKX))
  LL.effort.north=Effort.monthly.north_blocks%>%    
    filter(!METHOD=="GN")%>%
    group_by(FINYEAR,BLOCKX)%>%
    summarise(Effort=sum(hook.hours,na.rm=T))%>%
    mutate(dummy=paste(FINYEAR,BLOCKX))%>%
    dplyr::filter(dummy%in%GN.ktch.north$dummy)
  LL.ktch.north=Data.monthly.north%>%   
    filter(!METHOD=="GN")%>%
    group_by(FINYEAR,BLOCKX)%>%
    summarise(Ktch.LL=sum(LIVEWT.c,na.rm=T))%>%
    mutate(dummy=paste(FINYEAR,BLOCKX))%>%
    dplyr::filter(dummy%in%LL.effort.north$dummy)
  LL.cpue.north=left_join(LL.effort.north,LL.ktch.north,by=c('FINYEAR','BLOCKX'))%>%
    mutate(cpue.LL=Ktch.LL/Effort)
  GN.equiv.LL_effort.north=left_join(GN.ktch.north,LL.cpue.north,by=c('BLOCKX'))%>%
    mutate(GN.equiv.eff=Ktch.GN/cpue.LL)
} 
Effort.monthly.north_blocks=Effort.monthly.north_blocks%>%
        filter(FINYEAR%in%Relevant.yrs)%>%
        count(FINYEAR,BLOCKX)%>%
        group_by(FINYEAR,BLOCKX)%>%
        mutate(n=ifelse(n>0,1,0))
Effort.daily.north_blocks=Effort.daily.north_blocks%>%
        rename(FINYEAR=finyear,
               BLOCKX=blockx)%>%
        filter(FINYEAR%in%Relevant.yrs)%>%
        count(FINYEAR,BLOCKX)%>%
        group_by(FINYEAR,BLOCKX)%>%
        mutate(n=ifelse(n>0,1,0))
Effort.north_blocks=rbind(Effort.monthly.north_blocks,Effort.daily.north_blocks)%>%
        count(FINYEAR,BLOCKX)%>%
        group_by(FINYEAR,BLOCKX)%>%
        mutate(n=ifelse(n>0,1,0))%>%
        group_by(FINYEAR)%>%
        summarise(Tot=sum(n))%>%
        data.frame

Effort_blocks=rbind(Effort_blocks,Effort.north_blocks)%>%
  group_by(FINYEAR)%>%
  summarise(Tot=sum(Tot))%>%
  data.frame


#Add Species names to catch data sets    
All.species.names=All.species.names%>%
                      mutate(Name=tolower(Name))%>%
                      rename(SNAME=Name)
Hist.expnd=Hist.expnd%>%left_join(All.species.names,by='SPECIES')
Data.monthly=Data.monthly%>%left_join(All.species.names,by='SPECIES')
Data.monthly.north=Data.monthly.north%>%left_join(All.species.names,by='SPECIES')
Greynurse.ktch=Greynurse.ktch%>%left_join(All.species.names,by='SPECIES')
TEPS_dusky=TEPS_dusky%>%left_join(All.species.names,by='SPECIES')
WRL.ktch=WRL.ktch%>%left_join(All.species.names,by='SPECIES')
Taiwan.gillnet.ktch=Taiwan.gillnet.ktch%>%left_join(All.species.names,by='SPECIES')
Taiwan.longline.ktch=Taiwan.longline.ktch%>%left_join(All.species.names,by='SPECIES')
Indo_total.annual.ktch=Indo_total.annual.ktch%>%left_join(All.species.names,by='SPECIES')
GAB.trawl_catch=GAB.trawl_catch%>%left_join(All.species.names,by='SPECIES')
WTBF_catch=WTBF_catch%>%left_join(All.species.names,by='SPECIES')
Whaler_SA=Whaler_SA%>%left_join(All.species.names,by='SPECIES')

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

YEARS=sort(as.numeric(substr(unique(Data.monthly$FINYEAR),1,4)))  
Current=YEARS[length(YEARS)]   


#Explore spatial catch distribution to check species reporting issues
if(Explor=="YES")    
{
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
  
  fn.prop.exp=function(d,nm)
  {
    d1=d%>%mutate(FINYEAR=as.numeric(substr(FINYEAR,1,4)))%>%
      group_by(FINYEAR,BLOCKX)%>%
      summarise(Sum=sum(LIVEWT.c))%>%
      group_by(FINYEAR)%>%
      mutate(Prop=Sum/sum(Sum))%>%
      dplyr::select(FINYEAR,BLOCKX,Prop)%>%
      spread(BLOCKX,Prop)%>%
      data.frame
    colnames(d1)[-1]= substr(colnames(d1)[-1],2,100)
    CL=rainbow(ncol(d1)-1)
    plot(d1$FINYEAR,d1[,2],type='b',ylim=c(0,1),col=CL[1],pch=19,ylab="Proportion",xlab="Year")
    for(h in 3:ncol(d1)) points(jitter(d1$FINYEAR,1),d1[,h],type='b',col=CL[h-1],pch=19)
    legend("bottomright",colnames(d1)[-1],bty='n',col=CL,lty=1,lwd=2,cex=.6)
    mtext(nm,3)
  }
  these.ones.spt=c("Hammerheads","Spinner","Spurdogs","Tiger","Wobbegongs")
  smart.par(n.plots=length(these.ones.spt),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(l in 1:length(these.ones.spt))
  {
    fn.prop.exp(d=subset(Data.monthly,SPECIES==SP.list[[match(these.ones.spt[l],
                                                              names(SP.list))]]),nm=these.ones.spt[l])
  }
  dev.off()
}


#1. Plot catch data as reported in logbooks
ThIs=subset(Shark.species,!Shark.species%in%c(Indicator.species))
Data.monthly$Region="South"
Data.monthly.north$Region="North"
Tot.ktch=rbind(subset(Data.monthly,SPECIES%in%ThIs,select=c(FishCubeCode,FINYEAR,LIVEWT.c,SPECIES,SNAME,Region)),
               subset(Data.monthly.north,SPECIES%in%ThIs,select=c(FishCubeCode,FINYEAR,LIVEWT.c,SPECIES,SNAME,Region)))%>%
        mutate(finyear=as.numeric(substr(FINYEAR,1,4)),
               LIVEWT.c=LIVEWT.c/1000,        #in tonnes
               Name=ifelse(SPECIES%in%c(22999,31000),"unidentified sharks",SNAME))%>%
        group_by(Name,FishCubeCode,finyear,Region)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c))

Agg.r=Tot.ktch%>%
        group_by(Name,finyear,Region)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c))%>%
        spread(finyear,LIVEWT.c)%>%
        data.frame%>%
        arrange(Name)
colnames(Agg.r)[3:ncol(Agg.r)]=substr(colnames(Agg.r)[3:ncol(Agg.r)],2,20)
PCH=rep(19,nrow(Agg.r))
COL=rep(1,nrow(Agg.r))
Sp.fig.1=unique(Agg.r$Name)

HnDL=paste(hNdl,"/Outputs/Catch_by.sector/",sep="")
fnkr8t(HnDL)

  #Commercial
fn.fig(paste(HnDL,'commercial',sep=''),2400,2400) 
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

  #Recreational catch
Rec.ktch=Rec.ktch%>%mutate(Region=ifelse(zone%in%c('Gascoyne','North Coast'),'North','South'),
                           year=as.numeric(substr(FINYEAR,1,4)),
                           Common.Name=tolower(Common.Name))
Rec.sp=unique(Rec.ktch$Common.Name)
fn.fig(paste(HnDL,'recreational',sep=''),2400,2400) 
smart.par(n.plots=length(Rec.sp)+1,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
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

  #Taiwanese catch              
Taiwan.gillnet.ktch$Method="Pelagic.gillnet"
Taiwan.longline.ktch$Method="Longline"
Taiwan=rbind(Taiwan.longline.ktch,Taiwan.gillnet.ktch)%>%
              mutate(Region="North",
                     LIVEWT.c=LIVEWT.c/1000)%>%
              mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
              arrange(SNAME,year)
sp.taiwan=unique(Taiwan$SNAME)
sp.taiwan=sp.taiwan[!is.na(sp.taiwan)]
fn.fig(paste(HnDL,'taiwan',sep=''),2400,2400) 
smart.par(n.plots=length(sp.taiwan),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
for(i in 1:length(sp.taiwan))
{
  d=subset(Taiwan,SNAME==sp.taiwan[i])
  d.N=d%>%
    filter(Region=="North")%>%
    group_by(year)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))
  d.S=d%>%
    filter(Region=="South")%>%
    group_by(year)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))
  
  plot(sort(unique(Taiwan$year)),sort(unique(Taiwan$year)),col='transparent',cex=.8,ann=F,ylim=c(0,max(c(d.N$LIVEWT.c,d.S$LIVEWT.c,na.rm=T))))
  if(nrow(d.N)>0) points(d.N$year,d.N$LIVEWT.c,pch=PCH[i],type='o',col="grey60",cex=.8)
  if(nrow(d.S)>0) points(d.S$year,d.S$LIVEWT.c,pch=PCH[i],type='o',col="grey25",cex=.8)
  mtext(paste(sp.taiwan[i]),3,line=0.2,cex=0.8)  
  
}
legend('center',c("North","South"),lty=c(1,1),col=c("grey60","grey25"),lwd=2,bty='n',pch=19,cex=1.5)
mtext("Calendar year",1,line=0.5,cex=1.5,outer=T)
mtext("Total catch (tonnes)",2,las=3,line=0.35,cex=1.5,outer=T)
dev.off()

  #Indonesian fishing incursions catch              
Indo_total.annual.ktch=Indo_total.annual.ktch%>%filter(!is.na(SPECIES))
sp.indo=unique(Indo_total.annual.ktch$SNAME)
fn.fig(paste(HnDL,'indo',sep=''),2400,2400) 
smart.par(n.plots=length(sp.indo),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
for(i in 1:length(sp.indo))
{
  d=subset(Indo_total.annual.ktch,SNAME==sp.indo[i])%>%
    mutate(LIVEWT.c=LIVEWT.c/1000,          #in tonnes
           year=as.numeric(substr(FINYEAR,1,4)))%>%
    arrange(year)
  plot(d$year,d$LIVEWT.c,col="grey60",cex=.8,type='o',pch=PCH[i],ann=F,ylim=c(0,max(d$LIVEWT.c,na.rm=T)))
    mtext(paste(sp.indo[i]),3,line=0.2,cex=0.8)  
  
}
legend('center',c("North","South"),lty=c(1,1),col=c("grey60","grey25"),lwd=2,bty='n',pch=19,cex=1.5)
mtext("Calendar year",1,line=0.5,cex=1.5,outer=T)
mtext("Total catch (tonnes)",2,las=3,line=0.35,cex=1.5,outer=T)
dev.off()


#2. Remove species assessed elsewhere
  #2.1 remove indicator species and 'shark, other' (after catch reapportion)
Shark.species=subset(Shark.species,!Shark.species%in%c(Indicator.species,Shar_other,31000))

  #2.2 Remove blacktip sharks (blacktips and spot-tail) and school sharks (as per SAFS) 
#note: 
# . blacktips is a shared-stock with NT so NT assessment is used given their much higher catches (Grubert et al 2013)
# . school sharks are assumed to be a shared stock in the SESSF assessment (Thomson & Punt 2009)
blacktips=subset(Scien.nm,Scien.nm%in%c('C. limbatus & C. tilstoni','Carcharhinus sorrah'))$SPECIES
Shark.species=subset(Shark.species,!Shark.species%in%c(blacktips,School.shark))
  
  #2.3 Remove indicator species, blacktips  and school sharks
Data.monthly=subset(Data.monthly,SPECIES%in%Shark.species,select=c(SPECIES,FishCubeCode,SNAME,FINYEAR,LIVEWT.c,BLOCKX,METHOD))
Data.monthly.north=subset(Data.monthly.north,SPECIES%in%Shark.species,select=c(SPECIES,FishCubeCode,SNAME,FINYEAR,LIVEWT.c,BLOCKX,METHOD))

Data.monthly$Region="South"
Data.monthly.north$Region="North"


#3. Combine north and south
Tot.ktch=rbind(Data.monthly,Data.monthly.north)


#4. Some manipulations
SNAMEs=Data.monthly[!duplicated(Data.monthly$SPECIES),match(c("SPECIES","SNAME"),names(Data.monthly))]
SNAMEs.north=Data.monthly.north[!duplicated(Data.monthly.north$SPECIES),match(c("SPECIES","SNAME"),names(Data.monthly.north))]

#pull sawsharks as reported by species and as 'sawsharks'   
Tot.ktch=Tot.ktch%>%
        mutate(finyear=as.numeric(substr(FINYEAR,1,4)),
               SPECIES=ifelse(SPECIES%in%c(23002,23001,23900),23900,SPECIES), 
               SNAME=ifelse(SNAME%in%c("southern sawshark","common sawshark","sawsharks"),"sawsharks",SNAME),
               Name=SNAME)


#5. Add recreational catch
Rec.ktch=Rec.ktch%>%
              rename(finyear=year)%>%
              mutate(BLOCKX=NA,
                     Common.Name=ifelse(Common.Name=="dogfishes","spurdogs",
                                 ifelse(Common.Name=="greynurse shark","grey nurse shark",
                                 ifelse(Common.Name=="thresher shark","thresher sharks",
                                 ifelse(Common.Name=="bronze whaler","copper shark",
                                        Common.Name)))))%>%
              filter(Common.Name%in%unique(Tot.ktch$SNAME))%>%
              left_join(All.species.names%>%dplyr::select(-Scien.nm),by=c('Common.Name'='SNAME'))%>%
              mutate(SNAME=Common.Name,
                     Name=Common.Name,
                     METHOD='Rec.line',
                     FishCubeCode="Recreational")%>%
              dplyr::select(names(Tot.ktch))
Rec.ktch$Type="Recreational"
Tot.ktch$Type="Commercial"
Tot.ktch=rbind(Tot.ktch,Rec.ktch)


#6. Add Taiwanese catch
Taiwan=Taiwan%>%
          rename(finyear=year)%>%
          mutate(BLOCKX=NA,
                 Name=SNAME,
                 Type="Taiwan",
                 METHOD=Method,
                 LIVEWT.c=LIVEWT.c*1000,   #back in kg to match other fisheries
                 FishCubeCode="Taiwan")%>%
          dplyr::select(names(Tot.ktch))%>%
          filter(SPECIES%in%unique(Tot.ktch$SPECIES))
Tot.ktch=rbind(Tot.ktch,Taiwan)


#7. Add Indonesian fishing incursions
Indo=Indo_total.annual.ktch%>%
        mutate(BLOCKX=NA,
               Region="North",
               finyear=as.numeric(substr(FINYEAR,1,4)),
               Name=SNAME,
               Type="Indonesia",
               METHOD=NA,
               FishCubeCode="Indo")%>%
        dplyr::select(names(Tot.ktch))%>%
        filter(SPECIES%in%unique(Tot.ktch$SPECIES))
Tot.ktch=rbind(Tot.ktch,Indo)


#Keep relevant years
Tot.ktch=Tot.ktch%>%
          filter(FINYEAR%in%Relevant.yrs)


#8. Display Catch for each Species (in tonnes)  
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
    if(uno>2)lines(1:length(b.N$finyear),b.N$"LIVEWT.c",col="red",lwd=LWD)
    if(uno<=2)points(1:length(b.N$finyear),b.N$"LIVEWT.c",pch=19,cex=2,col="red")
  }
  
  d1=subset(d,Region=="South")
  if(nrow(d1)>0)
  {
    b.S=aggregate(LIVEWT.c/1000~finyear,d1,sum)
    uno=nrow(b.S)
    b.S=fn.add(D=b.S)
    if(uno>2)lines(1:length(b.S$finyear),b.S$"LIVEWT.c",col="forestgreen",lty=4,lwd=LWD)
    if(uno<=2)points(1:length(b.S$finyear),b.S$"LIVEWT.c",pch=19,cex=2,col="forestgreen")
  }
  
  mtext(SP,3,0,cex=1.5)
  legend("topleft",c("Total",expression(paste("North of 26 ",degree,"S")),
                     expression(paste("South of 26 ",degree,"S"))),bty='n',lty=c(1,1,4),
         lwd=LWD,col=c("black","red","forestgreen"),cex=1.5,pt.lwd=4)
  mtext("Catch (tonnes)",2,3,cex=2,las=3)
  mtext("Financial year",1,2.5,cex=2)
  axis(1,1:length(b.Tot$finyear),F,tck=-0.0125)
  axis(1,seq(1,length(b.Tot$finyear),10),Plot.yrs[seq(1,length(b.Tot$finyear),10)],cex.axis=1.25,tck=-0.025)
}

test=Tot.ktch[!duplicated(Tot.ktch$SPECIES),]
Uni.sp=test$SPECIES
names(Uni.sp)=test$Name
HnDl=paste(hNdl,"/Outputs/Catch_all_sp/",sep="")
fnkr8t(HnDl)
for(i in 1:length(Uni.sp))
{
  fn.fig(paste(HnDl,names(Uni.sp)[i],sep=""),2400,2400) 
  par(las=1,mgp=c(1,.8,0),mai=c(.8,1,.3,.1))
  Pt.ktch.sp(sp=Uni.sp[i],SP=names(Uni.sp)[i],LWD=3)
  dev.off()
}


#9. add TEP interactions        
Greynurse.ktch=Greynurse.ktch%>%
                mutate(BLOCKX=NA,
                       Region="South",
                       finyear=as.numeric(substr(FINYEAR,1,4)),
                       Name=SNAME,
                       Type="TEP",
                       METHOD='GN',
                       FishCubeCode="TEP")%>%
                dplyr::select(names(Tot.ktch))
Tot.ktch=rbind(Tot.ktch,Greynurse.ktch)


#10. add WRL        
WRL.ktch=WRL.ktch%>%
                mutate(BLOCKX=NA,
                       Region="South",
                       finyear=as.numeric(substr(FINYEAR,1,4)),
                       Name=SNAME,
                       Type="WRL",
                       METHOD='LL',
                       FishCubeCode="WRL")%>%
                dplyr::select(names(Tot.ktch))
Tot.ktch=rbind(Tot.ktch,WRL.ktch)

#11. Add historic  
a=Hist.expnd%>%
            mutate(BLOCKX=NA,
                   Region="South",
                   finyear=as.numeric(substr(FINYEAR,1,4)),
                   Type="Commercial",
                   METHOD=NA,
                   FishCubeCode="Historic",
                   SNAME=ifelse(SNAME%in%c('southern sawshark','common sawshark'),'sawsharks',SNAME),
                   SPECIES=ifelse(SPECIES%in%c(23001,23002),23900,SPECIES),
                   Name=SNAME)%>%
            dplyr::select(names(Tot.ktch))%>%
            filter(SPECIES%in%unique(Tot.ktch$SPECIES))
Tot.ktch=rbind(Tot.ktch,a)


#12. Add GAB
GAB=GAB.trawl_catch%>%
  mutate(BLOCKX=NA,
         Region="South",
         finyear=as.numeric(substr(FINYEAR,1,4)),
         Name=SNAME,
         Type="GAB",
         METHOD="trawl",
         FishCubeCode="GAB")%>%
  dplyr::select(names(Tot.ktch))%>%
  filter(SPECIES%in%unique(Tot.ktch$SPECIES))
Tot.ktch=rbind(Tot.ktch,GAB)

#13. Add WTBF
WTB=WTBF_catch%>%
  mutate(BLOCKX=NA,
         Region="South",
         finyear=as.numeric(substr(FINYEAR,1,4)),
         Name=SNAME,
         Type="WTB",
         METHOD="line",
         FishCubeCode="WTB")%>%
  dplyr::select(names(Tot.ktch))%>%
  filter(SPECIES%in%unique(Tot.ktch$SPECIES))
Tot.ktch=rbind(Tot.ktch,WTB)


#14. Add Whaler_SA (SA Marine Scalefish Fishery)
Whaler_SA=Whaler_SA%>%
  mutate(BLOCKX=NA,
         Region="South",
         finyear=as.numeric(substr(FINYEAR,1,4)),
         Name=SNAME,
         Type="SA MSF",
         METHOD="line",
         FishCubeCode="SA MSF")%>%
  dplyr::select(names(Tot.ktch))%>%
  filter(SPECIES%in%unique(Tot.ktch$SPECIES))
Tot.ktch=rbind(Tot.ktch,Whaler_SA)


#15. Remove blacktips & school shark because they are not assessed here. Ditto white sharks
Tot.ktch=subset(Tot.ktch,!Name%in%c('Blacktips','spot tail shark',
                                    'Spot tail shark',"School shark",
                                    "white shark"))


#16. Select species with enough data  
Agg=Tot.ktch%>%
  mutate(Gear=ifelse(METHOD%in%c("BS","BH","GN","HN","Pelagic.gillnet"),"net",
              ifelse(METHOD%in%c("DL","DV","EL","GL","HL","HR","HY",
                                 "LL","Longline","Rec.line","TL"),'line',
              ifelse(METHOD%in%c("FG","TW"),'trawl',
              ifelse(METHOD%in%c("FT","PT"),'trap',
                     NA)))))
Agg.r=Agg%>%
  group_by(Name,FINYEAR)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
  spread(FINYEAR,LIVEWT.c,fill=0)%>%
  data.frame
names(Agg.r)[-(1:2)]=substr(names(Agg.r)[-(1:2)],1,4)

PCH=rep(19,length(unique(Agg.r$Name)))
COL=rep(1,length(unique(Agg.r$Name)))

id=rowSums(Agg.r[,2:ncol(Agg.r)],na.rm=T)
names(id)=Agg.r$Name
id=rev(sort(id))
Agg.r=Agg.r[match(names(id),Agg.r$Name),]


#---PSA to determine which species to assess further------------------------------------------------------  
#note: aggregates the susceptibilities of multiple fleets (Micheli et al 2014)

Agg.PSA=Agg%>%
  filter(!is.na(Gear))%>%
  group_by(Name,Gear,FINYEAR)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
  spread(FINYEAR,LIVEWT.c,fill=0)%>%
  data.frame
names(Agg.PSA)[-(1:2)]=substr(names(Agg.PSA)[-(1:2)],2,5)

Agg.sp=unique(Agg.PSA$Name)
KIP=vector('list',length(Agg.sp))
for(s in 1:length(Agg.sp))
{
  d=Agg.PSA%>%filter(Name==Agg.sp[s])
  
  d1=d[,-c(1:2)]
  d1[d1<Min.ktch]=0
  d1[d1>=Min.ktch]=1
  kip=data.frame(Gear=d$Gear)%>%
    mutate(
      Name=d$Name,
      Sum.yrs=rowSums(d1),
      Keep=ifelse(Sum.yrs>=Min.yrs,"YES","NO"))
  KIP[[s]]=kip
}
KIP=do.call(rbind,KIP)%>%
        filter(Keep=="YES")%>%
        dplyr::select(Name,Gear)

Tot.ktch=Tot.ktch%>%filter(!Name%in%c("blacktips","dusky shark"))

PSA.list=PSA.list%>%filter(Species%in%unique(Tot.ktch$Name))
PSA.fn=function(d,Low.risk=2.64,medium.risk=3.18,Exprt)  #risk thresholds from Micheli et al 2014
{
  PSA=data.frame(Species=d$Species,
                 Productivity=rep(NA,nrow(d)),
                 Susceptibility=rep(NA,nrow(d)),
                 Vulnerability=rep(NA,nrow(d)))
  for(p in 1:nrow(d))    
  {
    aa=d[p,]
    if(!aa$Species%in%KIP$Name)
    {
      k=KIP%>%filter(Name==aa$Species)
      aa=aa%>%
        mutate(Net.avail=ifelse(!'net'%in%k$Gear,1,Net.avail),
               Net.encoun=ifelse(!'net'%in%k$Gear,1,Net.encoun),
               Line.avail=ifelse(!'line'%in%k$Gear,1,Line.avail),
               Line.encoun=ifelse(!'line'%in%k$Gear,1,Line.encoun),
               Trawl.avail=ifelse(!'trawl'%in%k$Gear,1,Trawl.avail),
               Trawl.encoun=ifelse(!'trawl'%in%k$Gear,1,Trawl.encoun),
               Trap.avail=ifelse(!'trap'%in%k$Gear,1,Trap.avail),
               Trap.encoun=ifelse(!'trap'%in%k$Gear,1,Trap.encoun))
    }
      
    PSA$Productivity[p]=mean(unlist(aa[,c('Max.age','Age.mat','Fecun',
                                          'Max.size','Size.mat','Rep.strat','Troph.Lvl')]))
    S1=1+((aa$Net.avail*aa$Net.encoun*aa$Net.sel*aa$Net.PCM)-1)/40
    S2=1+((aa$Line.avail*aa$Line.encoun*aa$Line.sel*aa$Line.PCM)-1)/40
    S3=1+((aa$Trawl.avail*aa$Trawl.encoun*aa$Trawl.sel*aa$Trawl.PCM)-1)/40
    S4=1+((aa$Trap.avail*aa$Trap.encoun*aa$Trap.sel*aa$Trap.PCM)-1)/40
    Cum.susc=min(3,(1+((S1-1)^2+(S2-1)^2+(S3-1)^2+(S4-1)^2)^0.5))
    PSA$Susceptibility[p]=Cum.susc
    PSA$Vulnerability[p]=(PSA$Productivity[p]^2+Cum.susc^2)^0.5  #Euclidean distance
  }
  
  PSA=PSA%>%
    mutate(Species=Hmisc::capitalize(as.character(Species)),
           Risk=factor(ifelse(Vulnerability<=Low.risk,'low',
                              ifelse(Vulnerability>Low.risk & Vulnerability<=medium.risk,'medium',
                                     'high')),levels=c('low','medium','high')))%>%
    arrange(Vulnerability)
  cols=c(low="green",medium="yellow",high="red")
  p=ggplot(PSA,
           aes(Productivity, Susceptibility, label = Species,colour = Risk, fill = Risk)) +
    geom_point(shape = 21, size = 6,colour="black") + 
    geom_text_repel(segment.colour='black',col='black',box.padding = 0.5) + 
    scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))+ 
    xlim(0.75,3.25)+ylim(0.75,3.25)+
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  p
  
  ggsave(Exprt, width = 8,height = 8, dpi = "retina")
  
  return(as.character(PSA%>%filter(Risk=="high")%>%pull(Species)))
}
Keep.species=PSA.fn(d=PSA.list,Exprt=paste(hNdl,"/Outputs/Figure. PSA.tiff",sep=''))
Keep.species=tolower(Keep.species)

#Plot catches of all species
all.yrs=min(Tot.ktch$finyear):max(Tot.ktch$finyear)
Fishry.type=sort(unique(Tot.ktch$Type))
colfunc <- colorRampPalette(c("red","yellow","springgreen","royalblue"))
COLs.type=colfunc(length(Fishry.type))
names(COLs.type)=Fishry.type
All.N.sp=sort(unique(Tot.ktch$Name))
fn.fig(paste(hNdl,'/Outputs/Figure 1_catch_all_species',sep=''),2400,2000) 
smart.par(n.plots=length(All.N.sp),MAR=c(1,1,.85,.25),OMA=c(2,2.25,.05,.05),MGP=c(1,.5,0))
par(cex.main=1,cex.axis=.85)
for(s in 1:length(All.N.sp))
{
  ddd=subset(Tot.ktch,Name==All.N.sp[s])%>%
    group_by(Type,finyear)%>%
    summarise(Tot=sum(LIVEWT.c,na.rm=T)/1000)%>%
    filter(!is.na(Tot))
  plot(all.yrs,all.yrs,col="transparent",ylab="",xlab="",main=capitalize(All.N.sp[s]),
       ylim=c(0,max(ddd$Tot)),xaxt='n',yaxt='n')
  unik.T=unique(ddd$Type)
  for(u in 1:length(unik.T))
  {
    cl=COLs.type[match(unik.T[u],names(COLs.type))]
    with(subset(ddd,Type==unik.T[u]),points(finyear,Tot,type='o',cex=.85,pch=21,bg=cl,col="grey10"))
  }
  Yrss=seq(all.yrs[1]-1,all.yrs[length(all.yrs)],5)
  axis(1,Yrss,F,tck=-.035)
  Yrss=seq(all.yrs[1]-1,all.yrs[length(all.yrs)],10)
  axis(1,Yrss,F,tck=-.07)
  if(s%in%26:31)axis(1,Yrss,Yrss,tck=-.07)
  Yax=pretty(seq(0,max(ddd$Tot),length.out = 4))
  axis(2,Yax,Yax, cex.axis=.85,padj =0.75,tck=.07,las=3)
  
}
plot.new()
legend('topleft',names(COLs.type)[1:4],pt.bg=COLs.type[1:4],bty='n',pch=21,cex=1.15,pt.cex=2)
plot.new()
legend('topleft',names(COLs.type)[5:length(COLs.type)],pt.bg=COLs.type[5:length(COLs.type)],
       bty='n',pch=21,cex=1.15,pt.cex=2)

mtext("Financial year",1,line=0.5,cex=1.5,outer=T)
mtext("Total catch (tonnes)",2,las=3,line=0.65,cex=1.5,outer=T)
dev.off()

#Analysed selected species
Tot.ktch=subset(Tot.ktch,Name%in%Keep.species)    


#Species grouping   
Tot.ktch$SP.group=Tot.ktch$Name

Specs=Tot.ktch[!duplicated(Tot.ktch$SP.group),match(c("SPECIES","SNAME","Name","SP.group"),names(Tot.ktch))]
Dis.sp=unique(Specs$SP.group)
Specs=Specs[order(Specs$SP.group),]
N.sp=nrow(Specs)

#---Life history parameters of selected species------------------------------------------------------  
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

#Fill in life history pars  
SPLF=LH.par$SP.group
GROWTH.F=vector('list',N.sp)
names(GROWTH.F)=SPLF
MAX.age.F=AGE.50.mat=FECU=Repro_cycle=Aver.Lat=AVER.T=GROWTH.F
for(i in 1:N.sp)
{
  dd=subset(LH.par,SP.group==SPLF[i])
  GROWTH.F[[i]]=data.frame(k=dd$K,FL_inf=dd$FL_inf)
  MAX.age.F[[i]]=c(dd$Max_Age,round(dd$Max_Age*1.4))
  AGE.50.mat[[i]]=c(dd$Age_50_Mat_min,dd$Age_50_Mat_max)
  FECU[[i]]=c(dd$Fecu_min,dd$Fecu_max)
  Repro_cycle[[i]]=LH.par$Cycle[i]
  AVER.T[[i]]=LH.par$Temperature[i]
}

#Export Life history table
setwd(paste(hNdl,'/Outputs',sep=''))

TabL=LH.par[,-match(c("SNAME","SP.group","Scien.nm"),names(LH.par))]
TabL=TabL[order(TabL$Name),-match("SPECIES",names(TabL))]
TabL$K=as.numeric(as.character(TabL$K))
TabL$K=round(TabL$K,3)
TabL$FL_inf=round(TabL$FL_inf)
names(TabL)[match('FL_inf',names(TabL))]='L_inf'
fn.word.table(WD=getwd(),TBL=TabL,Doc.nm="Table 1. Life history pars",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")


#---Final manipulations------------------------------------------------------  

#Catch in tonnes
Tot.ktch$LIVEWT.c=Tot.ktch$LIVEWT.c/1000

# Check available tagging data
FL.sp=c('Bull shark','Great hammerhead','Lemon shark','Pigeye shark','Scalloped hammerhead',
        'Smooth hammerhead','Spinner shark','Spurdogs','Tiger shark',' Wobbegong (general)',
        "Grey nurse shark","Hammerheads","Shortfin mako","Pencil shark",
        "Lemon shark","Milk shark","Copper shark","Common sawshark")
if(use.tags)
{
  Tag=Tag%>%filter(COMMON_NAME%in%FL.sp & Recaptured=="Yes")
  table(Tag$COMMON_NAME)
}

# Set catch of all species starting in 1940s
TAB.dummy=with(Tot.ktch,table(FINYEAR,SP.group))
Add.yrs=paste(seq(1941,1974),substr(seq(1942,1975),3,4),sep="-")
id.dummy=match(Add.yrs[1],rownames(TAB.dummy)):match(Add.yrs[length(Add.yrs)],rownames(TAB.dummy))
TAB.dummy=TAB.dummy[id.dummy,]
Nms.dummy=names(which(colSums(TAB.dummy)<=1))
if(length(Nms.dummy)>0)
{
  Add.dummy=Tot.ktch[id.dummy,]%>%
    replace(.,,NA)%>%
    mutate(Type='Commercial',
           LIVEWT.c=0,
           FINYEAR=Add.yrs,
           finyear=as.numeric(substr(FINYEAR,1,4)))
  
  list.dummy=vector('list',length(Nms.dummy))
  for(l in 1:length(Nms.dummy))
  {
    ss=Add.dummy
    dd=Tot.ktch%>%filter(SP.group==Nms.dummy[l])%>%slice(1)
    ss=ss%>%mutate(SPECIES=dd$SPECIES,
                   SNAME=dd$SNAME,
                   Name=dd$Name,
                   SP.group=dd$SP.group)
    if(Nms.dummy[l]=="Scalloped hammerhead") ss=ss%>%filter(!FINYEAR=="1974-75")
    list.dummy[[l]]=ss
  }
  list.dummy=do.call(rbind,list.dummy)
  Tot.ktch=rbind(Tot.ktch,list.dummy)%>%arrange(SP.group,finyear)  
}

#Export total catch of each species
hn.ktch.sp=paste(hNdl,"/Outputs/Catch_assessed_sp/",sep="")
fnkr8t(hn.ktch.sp)
for(s in 1: N.sp)
{
  ddd=subset(Tot.ktch,SP.group==Specs$SP.group[s])
  ddd=aggregate(LIVEWT.c~finyear,ddd,sum)  
  write.csv(ddd,paste(hn.ktch.sp,'Total.annual.catch_',Specs$SP.group[s],
                      '.csv',sep=""),row.names = F)
}

#Displayed catches for analysed species
fn.fig(paste(hNdl,'/Outputs/Figure 2_catch_analysed_species',sep=''),2400,2400) 
smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
for(s in 1: N.sp)
{
  ddd=subset(Tot.ktch,SP.group==Specs$SP.group[s])%>%
            group_by(Type,finyear)%>%
            summarise(Tot=sum(LIVEWT.c,na.rm=T))%>%
            filter(!is.na(Tot))
  plot(all.yrs,all.yrs,col="transparent",ylab="",xlab="",main=capitalize(Specs$SP.group[s]),
       ylim=c(0,max(ddd$Tot)),xaxt='n')
  unik.T=unique(ddd$Type)
  for(u in 1:length(unik.T))
  {
       cl=COLs.type[match(unik.T[u],names(COLs.type))]
      with(subset(ddd,Type==unik.T[u]),points(finyear,Tot,type='o',cex=1.1,pch=21,bg=cl))
  }
  axis(1,all.yrs,F,tck=-.015)
  axis(1,seq(all.yrs[1],all.yrs[length(all.yrs)],10),
       seq(all.yrs[1],all.yrs[length(all.yrs)],10),tck=-.03)
}
mtext("Financial year",1,line=0.5,cex=1.5,outer=T)
mtext("Total catch (tonnes)",2,las=3,line=0.35,cex=1.5,outer=T)
plot.new()
legend('top',names(COLs.type),pt.bg=COLs.type,
               bty='n',pch=21,cex=1.1,pt.cex=2)
dev.off()


#---Collate available abundance data------------------------------------------------------  
Scal.hh.nat$CV=Scal.hh.nat$CV/100
Tiger.nat$CV=Tiger.nat$CV/100
Mil.nat$CV=Mil.nat$CV/100

cpue.list=list(
  "copper shark"=list(Nat=NULL,
                      TDGDLF.mon=NULL,
                      TDGDLF.day=NULL,
                      NSF=NULL),
  "grey nurse shark"=list(Nat=NULL,
                          TDGDLF.mon=Greynurse.tdgdlf_mon,
                          TDGDLF.day=NULL,
                          NSF=NULL),
  "lemon shark"=list(Nat=NULL,
                     TDGDLF.mon=NULL,
                     TDGDLF.day=NULL,
                     NSF=Lemon.NSF),
  "great hammerhead"=NULL,
  "milk shark"=list(Nat=Mil.nat,     #flat cpue, no signal, error in estimation
                    TDGDLF.mon=NULL,
                    TDGDLF.day=NULL,
                    NSF=NULL),
  "sawsharks"=NULL,
  "scalloped hammerhead"=list(Nat=Scal.hh.nat,
                              TDGDLF.mon=NULL,
                              TDGDLF.day=NULL,
                              NSF=NULL),
  "smooth hammerhead"=list(Nat=NULL,
                           TDGDLF.mon=Smuz.hh.tdgdlf_mon,
                           TDGDLF.day=Smuz.hh.tdgdlf_daily,
                           NSF=NULL), 
  "spinner shark"=list(Nat=NULL,
                       TDGDLF.mon=Spinr.tdgdlf_mon,
                       TDGDLF.day=Spinr.tdgdlf_daily),
  "shortfin mako"=list(Nat=NULL,
                       TDGDLF.mon=Mako.tdgdlf_mon,  
                       TDGDLF.day=Mako.tdgdlf_daily,
                       NSF=NULL),
  "spurdogs"=NULL,
  "tiger shark"=list(Nat=Tiger.nat,
                     TDGDLF.mon=Tiger.tdgdlf_mon,
                     TDGDLF.day=Tiger.tdgdlf_daily,
                     NSF=Tiger.NSF),         
  "wobbegongs"=NULL) 


#---Build r prior -----------------------------------------------------------------------
fun.rprior.dist=function(Nsims,K,LINF,Temp,Amax,MAT,FecunditY,Cycle)
{
  Fecu=unlist(FecunditY)
  Rprior=fun.Leslie(N.sims=Nsims,k=K,Linf=LINF,Aver.T=Temp,
                    A=Amax,first.age=0,RangeMat=MAT,Rangefec=Fecu,
                    sexratio=0.5,Reprod_cycle=Breed.cycle,Hoenig.only="NO")  
  
  #get mean and sd from gamma and normal distribution
  normal.pars=suppressWarnings(fitdistr(Rprior$r.prior, "normal"))
  gamma.pars=suppressWarnings(fitdistr(Rprior$r.prior, "gamma"))  
  shape=gamma.pars$estimate[1]        
  rate=gamma.pars$estimate[2]      
  return(list(shape=shape,rate=rate,
              mean=normal.pars$estimate[1],sd=normal.pars$estimate[2],
              M=Rprior$M))
  
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
store.species.M=store.species
PATH=paste(getwd(),"each_species",sep='/')
if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))   
setwd(file.path(PATH))
PHat=file.path(PATH)
system.time(for(s in 1: N.sp) #get r prior    #takes 0.013 sec per iteration 
{
  #life history
  Growth.F=GROWTH.F[[s]]
  TEMP=AVER.T[[s]]
  Max.age.F=MAX.age.F[[s]]
  Age.50.mat=AGE.50.mat[[s]]
  Fecundity=FECU[[s]]
  Breed.cycle=Repro_cycle[[s]]  #years
  print(paste("r prior ",s,"--",names(store.species)[s]))
  #Get r prior
  r.prior.dist=fun.rprior.dist(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf/.85,Temp=TEMP,Amax=Max.age.F,
                               MAT=unlist(Age.50.mat),FecunditY=Fecundity,Cycle=Breed.cycle)
  store.species[[s]]$r.prior=list(shape=r.prior.dist$shape,rate=r.prior.dist$rate)
  store.species[[s]]$r.prior.normal=list(mean=r.prior.dist$mean,sd=r.prior.dist$sd)
  store.species.M[[s]]=r.prior.dist$M
})

#---Assign Resilience -----------------------------------------------------------------------
RESILIENCE=store.species
for(r in 1:length(RESILIENCE))
{
  RESILIENCE[[r]]=with(store.species[[r]]$r.prior.normal,
                                ifelse(mean>=0.6,"High",
                                ifelse(mean<0.6 & mean>=0.2,'Medium',
                                ifelse(mean<0.2 & mean>=0.1,'Low',
                                      'Very low'))))
}

# ---Derive age frequency and selectivity from length frequency -----------------------------------------------------------------------
WD=getwd()
if(use.size.comp) 
{
  User="Matias"
  source('C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R')
  DATA=DATA%>%
      mutate(COMMON_NAME=ifelse(COMMON_NAME=="Bronze whaler","Copper shark",
                         ifelse(COMMON_NAME=="Shortfin mako ","Shortfin mako",COMMON_NAME)))
  
  Res.vess=c('FLIN','NAT',"HAM","HOU","RV BREAKSEA","RV Gannet",
             "RV GANNET","RV SNIPE 2")
  fun.check.LFQ=function(a,area)
  {
    a=a%>%
      dplyr::select(c(SPECIES,COMMON_NAME,SEX,year,BOAT,
                      Method,FL,Mid.Lat,Mid.Long,
                      MESH_SIZE))%>%
      filter(!is.na(FL))%>%
      mutate(MESH_SIZE=sub("\"","",MESH_SIZE))
    Unik.sp=unique(a$COMMON_NAME)
    
    smart.par(n.plots=length(Unik.sp),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
    b1=table(a$SEX,a$COMMON_NAME)
    for(u in 1:length(Unik.sp))
    {
      bb=b1[,match(Unik.sp[u],colnames(b1))]
      if(sum(bb)>2)
      {
        Xsq <- chisq.test(bb)
        barplot(bb,main=paste(Unik.sp[u],"(F:M=",round(bb[1]/bb[2],2),":1",", X2, p=",
                              round(Xsq$p.value,2),")"),cex.main=.9)
      }
    } 
    mtext(area,1,outer = T,cex=1.5)
    
    smart.par(n.plots=length(Unik.sp),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
    b1=table(10*(round(a$FL/10)),a$COMMON_NAME)
    for(u in 1:length(Unik.sp))
    {
      barplot(b1[,match(Unik.sp[u],colnames(b1))],
              main=paste(Unik.sp[u]," (n= ",sum(b1[,match(Unik.sp[u],colnames(b1))]),")",sep=""))
    }
    mtext(area,1,outer = T,cex=1.5)
    
    
    for(u in 1:length(Unik.sp))
    {
      b1=with(subset(a,COMMON_NAME==Unik.sp[u]),table(10*(round(FL/10)),year))
      yrs=colnames(b1)
      smart.par(n.plots=ncol(b1),MAR=c(2,2,1,1),OMA=c(1.75,2,2,.1),MGP=c(1,.5,0))
      for(s in 1:ncol(b1))
      {
        barplot(b1[,match(yrs[s],colnames(b1))],main=paste(yrs[s],"n=",
                                                           sum(b1[,match(yrs[s],colnames(b1))])))
      }
      mtext(paste(area,"-",unique(a$Method),"-",Unik.sp[u]),3,outer=T)
    }
    
    return(a)
  }
  PATH=paste(hNdl,"Outputs/Size.frequency",sep='/')
  if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
  
  pdf(paste(hNdl,"/Outputs/Size.frequency/Available_data_north.pdf",sep=""))
  LFQ.north=fun.check.LFQ(a=DATA%>%filter(Mid.Lat>(-26) & COMMON_NAME%in%FL.sp &
                                            Method%in%c('LL')  & !BOAT%in%Res.vess),
                          area="North")
  dev.off()
  
  pdf(paste(hNdl,"/Outputs/Size.frequency/Available_data_south.pdf",sep=""))  
  LFQ.south=fun.check.LFQ(a=DATA%>%filter(Mid.Lat<(-26) & COMMON_NAME%in%FL.sp &
                                            Method%in%c('GN') &!BOAT%in%Res.vess),
                          area="South")
  dev.off()
  
  #Get age from size using age-length-key constructed based on assumed variability
  source("C:/Matias/Analyses/Population dynamics/Git_Stock.assessments/Age.length.key.R")  
  Conversion=data.frame(Name=c("copper shark","great hammerhead","grey nurse shark","lemon shark",
                               "milk shark","sawsharks","scalloped hammerhead","shortfin mako",
                               "smooth hammerhead","spinner shark","spurdogs",
                               "tiger shark","wobbegongs"),
                        FL.sp=c("Copper shark","Great hammerhead","Grey nurse shark","Lemon shark",
                                "Milk shark","Common sawshark","Scalloped hammerhead","Shortfin mako",
                                "Smooth hammerhead","Spinner shark","Spurdogs",
                                "Tiger shark","Wobbegong (general)"))
  Store.age.comp=vector('list',nrow(LH.par))
  names(Store.age.comp)=LH.par$SPECIES
  for(l in 1:nrow(LH.par))
  {
    NM=as.character(Conversion%>%filter(Name==LH.par$Name[l])%>%pull(FL.sp))
    Linf=LH.par$FL_inf[l] 
    Lo=LH.par$LF_o[l] 
    k=LH.par$K[l] 
    AMAX=LH.par$Max_Age[l] 
    b_w8t=LH.par$b_w8t[l]
    a_w8t=LH.par$a_w8t[l]
    
    mn.len=Lo+(Linf-Lo)*(1-exp(-k*(0:AMAX)))
    
    GN=get.prop.at.age.from.length(
            age=0:AMAX,
            mn.len=mn.len,
            SD=mn.len*seq(.15,.025,length.out = AMAX+1),
            N=100,
            int=10,
            Obs.len=DATA%>%
              filter(Method=='GN'& 
                       COMMON_NAME==NM & 
                       !BOAT%in%Res.vess)%>%
              mutate(FL=ifelse(is.na(FL),TL*Mn.conv.Fl.Tl,FL))%>%   
              filter(!is.na(FL))%>%
              filter(FL>=Min.len)%>%
              pull(FL),
            min.obs=Min.obs)
    
    LL=get.prop.at.age.from.length(
            age=0:AMAX,
            mn.len=mn.len,
            SD=mn.len*seq(.15,.025,length.out = AMAX+1),
            N=100,
            int=10,
            Obs.len=DATA%>%
              filter(Method=='LL'& 
                       COMMON_NAME==NM & 
                       !BOAT%in%Res.vess)%>%
              mutate(FL=ifelse(is.na(FL),TL*Mn.conv.Fl.Tl,FL))%>%   
              filter(!is.na(FL))%>%
              filter(FL>=Min.len)%>%
              pull(FL),
            min.obs=Min.obs)
    
    #fill in objects for aSPM
    props=data.frame(age=floor(0:AMAX),
                     laa=mn.len,
                     waa= a_w8t*mn.len^b_w8t,  #catch in tonnes; waa in kgs 
                     maa=plogis(floor(0:AMAX),
                                floor(mean(unlist(LH.par[l,c('Age_50_Mat_min',
                                                             'Age_50_Mat_max')]))),1),
                     sela=NA,
                     feca=NA)    
    
    glb=list(maxage=AMAX,
             M=NA,
             Linf=Linf,
             K=k,
             t0=NA,
             Waa=NA,
             Wab=NA,
             M50a=floor(mean(unlist(LH.par[l,c('Age_50_Mat_min','Age_50_Mat_max')]))),
             deltaM=NA,
             steep=NA,
             R0=NA,
             sela50=NA,
             deltaS=NA,
             resilience=RESILIENCE[[l]],
             nages=length(floor(0:AMAX)),
             ages=floor(0:AMAX),
             nyrs=NA,
             spsname=names(RESILIENCE)[l])
    Store.age.comp[[l]]=list(GN=GN,LL=LL,props=props,glb=glb)
    rm(GN,LL)
  }
  
  
  fn.plt.age.comp=function(sim,obs,sel,Title)   
  {
    if(!is.null(sim))
    {
      with(sim,
           {
             plot(age,Mean.len,type='l',lwd=3,xlab="Age",ylab="Length (cm)",
                  main=Title,ylim=c(0,max(len)))
             points(age,len)
           })
      with(obs,points(age,len,pch=19,col="orange"))
      par(fig = c(0.4,0.95, 0.05, .6), new = T,mgp=c(2,.5,0),las=1)  
      aa=prop.table(table(obs$age))
      plot(aa/max(aa),cex.axis=.75,
           main=paste("Age composition (n=",nrow(obs),")",sep=''),xlab="",ylab="Proportion")
      lines(sel$x,sel$y,col="forestgreen",lwd=3)
      legend("topright","Vulnerability",bty='n',text.col="forestgreen")
    }
  }
  
  pdf(paste(hNdl,"/Outputs/Size.frequency/Inferred.age.comp.and.selectivity.pdf",sep=""))
  for(l in 1:nrow(LH.par))
  {
    par(mfcol=c(1,1))
    fn.plt.age.comp(sim=Store.age.comp[[l]]$GN$dat,
                    obs=Store.age.comp[[l]]$GN$pred.age,
                    sel=Store.age.comp[[l]]$GN$Selectivity,
                    Title=paste(LH.par$Name[l],"(gillnet)"))
    par(mfcol=c(1,1))
    fn.plt.age.comp(sim=Store.age.comp[[l]]$LL$dat,
                    obs=Store.age.comp[[l]]$LL$pred.age,
                    sel=Store.age.comp[[l]]$LL$Selectivity,
                    Title=paste(LH.par$Name[l],"(longline)"))
  }
  dev.off()
  
  #combine selectivities of gillnets and longlines
  Selectivity.at.age=vector('list',nrow(LH.par))
  names(Selectivity.at.age)=LH.par$SPECIES
  for(l in 1:length(Selectivity.at.age))
  {
    dummy.GN=dummy.LL=NULL
    L=Store.age.comp[[l]]$GN$Selectivity
    if(!is.null(L))dummy.GN=data.frame(matrix(unlist(L), ncol=length(L), byrow=F))
    L=Store.age.comp[[l]]$LL$Selectivity
    if(!is.null(L))dummy.LL=data.frame(matrix(unlist(L), ncol=length(L), byrow=F))
    dummy=rbind(dummy.GN,dummy.LL)  
    if(!is.null(dummy))
    {
      names(dummy)=c('age','relative.sel')
      dummy=dummy%>%
        group_by(age)%>%
        summarise(relative.sel=max(relative.sel))
      dummy$relative.sel=dummy$relative.sel/max(dummy$relative.sel)
      Selectivity.at.age[[l]]=dummy
    }
    print(paste("selectivity ",l,"--",names(Selectivity.at.age)[l]))
  }
  
  #Collate catch size frequency
  LFQ.south$COMMON_NAME=tolower(LFQ.south$COMMON_NAME)
  ktch.size.fq=vector('list',N.sp)
  names(ktch.size.fq)=Specs$SP.group
  for(k in 1:N.sp)
  {
    dd=LFQ.south%>%filter(COMMON_NAME==names(ktch.size.fq)[k] & MESH_SIZE%in%c("6.5","7"))
    if(nrow(dd)>Min.size.sample) ktch.size.fq[[k]]=dd
    rm(dd)
  }

}
setwd(WD)
#---Calculate steepness -----------------------------------------------------------------------
fun.steepness=function(Nsims,K,LINF,first.age,sel.age,F.mult,Temp,Amax,MAT,FecunditY,Cycle,sexratio,spawn.time)
{
  Fecu=unlist(FecunditY)
  fn.draw.samples=function(A=Amax,RangeMat=MAT,Rangefec=Fecu,Reprod_cycle=Breed.cycle)
  {
    #Max Age
    if(length(A)==1) Max.A=A
    if(length(A)>1)if(A[1]==A[2]) Max.A=A[1]
    if(length(A)>1)if(A[1]<A[2]) Max.A=ceiling(rtriangle(1,a=A[1],b=A[2],c=ceiling((A[1]+A[2])/2)))
    
    #Age vector
    age=first.age:Max.A
    
    #fecundity at age
    if(Rangefec[1]==Rangefec[2]) Meanfec.sim=rep(Rangefec[1],length(age))
    if(Rangefec[1]<Rangefec[2]) Meanfec.sim=rep(ceiling(rtriangle(1,a=Rangefec[1],b=Rangefec[2],
                                                                  c=ceiling((Rangefec[1]+Rangefec[2])/2))),length(age)) 
    
    #Age at 50% maturity
    if(RangeMat[1]==RangeMat[2]) age.mat.sim=ceiling(RangeMat[1])
    if(RangeMat[1]<RangeMat[2]) age.mat.sim=ceiling(runif(1,RangeMat[1],RangeMat[2]))  
    
    #Reproductive cycle
    if(length(Reprod_cycle)==1) Rep_cycle.sim=Reprod_cycle else
    {
      if(Reprod_cycle[1]==Reprod_cycle[2]) Rep_cycle.sim=round(Reprod_cycle[1])
      if(Reprod_cycle[1]<Reprod_cycle[2]) Rep_cycle.sim=round(runif(1,Reprod_cycle[1],Reprod_cycle[2]))    
    }
    
    return(list(Max.A=Max.A,age.mat=age.mat.sim,Meanfec=Meanfec.sim,Rep_cycle=Rep_cycle.sim))    
  }
  
  M.fun=function(Amax,age.mat,Hoenig.only="NO",Linf=LINF,k=K,Aver.T=Temp)
  {
    m.Jensen.2=1.65/age.mat
    m.Jensen.2=rep(m.Jensen.2,length(age))
    
    #Pauly (1980)  
    m.Pauly=10^(-0.0066-0.279*log10(Linf)+0.6543*log10(k)+0.4634*log10(Aver.T))
    m.Pauly=rep(m.Pauly,length(age))
    
    #Hoenig (1983), combined teleost and cetaceans    
    m.Hoenig=exp(1.44-0.982*log(Amax))      
    m.Hoenig=rep(m.Hoenig,length(age))
    
    #Then et al (2015)
    m.Then.1=4.899*Amax^(-0.916)
    m.Then.1=rep(m.Then.1,length(age))
    
    
    #STEP 2. get mean at age
    if(Hoenig.only=="NO")nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Hoenig,m.Then.1)  
    if(Hoenig.only=="YES")nat.mort=data.frame(m.Hoenig)
    
    return(rowMeans(nat.mort))
    apply(nat.mort, 1, function(x) weighted.mean(x, c(1,1,1.5,1.5)))
  }
  
  Stipns=function(max.age,M,age.mat,Meanfec,CyclE,Sel)
  {  
    #survivorship
    surv=exp(-M)
    
    #fecundity  
    fecundity=Meanfec*sexratio/CyclE
    
    #maturity
    #maturity=ifelse(age>=age.mat,1,0)   #knife edge
    maturity=plogis(age,age.mat,1)      #ogive
    
    # maximum age is plus group
    phi.o=0.0
    cum.survive=1.0
    z=0.0
    for (i in 2:(max.age)  )
    {
      z=M[i] + F.mult*Sel[i]
      z.ts=(M[i]+F.mult*Sel[i])*spawn.time
      phi.o=phi.o+cum.survive*fecundity[i]*maturity[i]*exp(-z.ts)
      cum.survive=cum.survive*exp(-z )
    }
    #plus group  
    z= M[max.age+1] + F.mult*Sel[max.age+1]
    z.ts=(M[max.age+1]+F.mult*Sel[max.age+1])*spawn.time
    phi.o=phi.o + fecundity[max.age+1]*maturity[max.age+1]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
    
    #maximum lifetime reproductive rate at low density
    alpha=phi.o*surv[1]
    
    #steepness
    h=alpha/(4+alpha)
    
    #spawning potential ratio at maximum excess recruitment (MER) (Beverton-Holt relationship)
    SPR.mer=1/alpha^0.5
    
    #optimal depletionlevel (i.e.depletion at MER, the proportional reduction from unexploited level)
    Dep.MER=((alpha^0.5)-1)/(alpha-1) 
    
    return(steepness=h)  
  }
  
  Store=rep(NA,Nsims)
  for(i in 1:Nsims)
  {
    a=fn.draw.samples()
    A.sim=a$Max.A
    age=first.age:A.sim
    Age.mat.sim=a$age.mat
    Meanfec.sim=a$Meanfec
    Reprod_cycle.sim=a$Rep_cycle
    M.sim=M.fun(Amax=A.sim,age.mat=Age.mat.sim)
    Sel=sel.age
    if(length(Sel)<length(age))
    {
      Sel=c(Sel,rep(Sel[length(Sel)],length(age)-length(Sel)))
    }
    hh=Stipns(max.age=A.sim,M=M.sim,age.mat=Age.mat.sim,
              Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim,
              Sel)
    #avoid non-sense h
    if(hh<0.21)repeat 
    {
      a=fn.draw.samples()
      A.sim=a$Max.A
      age=first.age:A.sim
      Age.mat.sim=a$age.mat
      Meanfec.sim=a$Meanfec
      Reprod_cycle.sim=a$Rep_cycle
      M.sim=M.fun(Amax=A.sim,age.mat=Age.mat.sim)
      Sel=sel.age
      if(length(Sel)<length(age))
      {
        Sel=c(Sel,rep(Sel[length(Sel)],length(age)-length(Sel)))
      }
      hh=Stipns(max.age=A.sim,M=M.sim,age.mat=Age.mat.sim,
                Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim,
                Sel)
      if(hh>0.22)break
    }
    Store[i]=hh
  }
  
  #get mean and sd from lognormal distribution
  normal.pars=suppressWarnings(fitdistr(Store, "normal"))
  gamma.pars=suppressWarnings(fitdistr(Store, "gamma"))  
  shape=gamma.pars$estimate[1]        
  rate=gamma.pars$estimate[2]  
  
  return(list(shape=shape,rate=rate,
              mean=normal.pars$estimate[1],
              sd=normal.pars$estimate[2]))
}

store.species.steepness=vector('list',N.sp)
names(store.species.steepness)=Specs$SP.group
system.time(for(s in 1: N.sp) 
{
  #life history
  Growth.F=GROWTH.F[[s]]
  TEMP=AVER.T[[s]]
  Max.age.F=MAX.age.F[[s]]
  Age.50.mat=AGE.50.mat[[s]]
  Fecundity=FECU[[s]]
  Breed.cycle=Repro_cycle[[s]]  #years
  SEL=Selectivity.at.age[[s]]$relative.sel
  
  print(paste("steepness ",s,"--",names(store.species)[s]))
  store.species.steepness[[s]]=fun.steepness(Nsims=2*NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf/.85,
                                             first.age=0,sel.age=SEL,F.mult=0,Temp=TEMP,
                                             Amax=Max.age.F,MAT=unlist(Age.50.mat),
                                             FecunditY=Fecundity,Cycle=Breed.cycle,
                                             sexratio=0.5,spawn.time = 0)
  
})  

#show steepness
COl.all.sp=colfunc(length(store.species.steepness))
names(COl.all.sp)=names(store.species.steepness)  

fn.fig(paste(hNdl,'/Outputs/Steepness',sep=''),2400,2400) 
plot(1:10,xlim=c(0,1),ylim=c(0,20),col='transparent',xlab="Steepness",yaxt='n',ylab="Density")
for(s in 1:length(store.species.steepness))
{
  dist=with(store.species.steepness[[s]],density(rnorm(100000,mean,sd),adjust = 2))
  dist=with(store.species.steepness[[s]],density(rgamma(10000,shape,rate),adjust = 2))
  lines(dist,col=COl.all.sp[s],lwd=2)
}
legend('topright',names(COl.all.sp),lty=1,lwd=2,col=COl.all.sp,bty='n')
dev.off()


#---Gear Selectivity estimation------------------------------------------------------
if(do.length.based=="YES")
{
  library(TropFishR)

  #get gillnet selectivity schedule from different experimental nets
  Estim.sel.exp='NO'  #not enough observations from different mesh sizes
  Size.sel=vector('list',length(ktch.size.fq))
  names(Size.sel)=names(ktch.size.fq)
  if(Estim.sel.exp=="YES")
  {
    fn.sel.stim=function(d,size.int)
    {
      #put data as a list
      tab=d%>%mutate(Size.class=size.int*floor(FL/size.int)+size.int/2)%>%
        group_by(MESH_SIZE,Size.class)%>%
        summarise(n=n())%>%
        spread(MESH_SIZE,n,fill=0)%>%
        data.frame
      colnames(tab)[-1]= substr(colnames(tab)[-1],2,100)
      d.list=list(midLengths=tab$Size.class,
                  type="gillnet",
                  meshSizes=2.54*as.numeric(colnames(tab)[-1]),   #inches to cm
                  CatchPerNet_mat=as.matrix(tab[,-1]))
      #apply sel estim methods
      out <- TropFishR::select(param = d.list) 
      
      #plot(out)
      return(list(out=out))
    }
    for(s in 1:length(Size.sel)) Size.sel[[s]]=fn.sel.stim(d=ktch.size.fq[[s]],size.int=10)
  }
  
  #Simple selectivity
  if(Estim.sel.exp=="NO")
  {
    fn.sel.simple=function(d)
    {
      dist=fitdistr(d$FL, "normal")
      FL=seq(min(d$FL),max(d$FL))
      Sel=dnorm(FL, dist$estimate[1], dist$estimate[2])
      
      tab=d%>%mutate(Size.class=floor(FL))%>%
        group_by(Size.class)%>%
        summarise(n=n())
      plot(tab$Size.class,tab$n,type='h',ylab="Frequency",xlab="FL (cm)",
           main=names(ktch.size.fq)[s])
      par(new=T)
      plot(FL,Sel/max(Sel),ann=F,yaxt='n',type='l',col=2,lwd=2)
      return(sel=data.frame(FL=FL,Sel=Sel/max(Sel)))
    }
    smart.par(length(Size.sel),c(3,3,4,2),c(1,1,1,1),c(2,.7,0))
    for(s in 1:length(Size.sel)) Size.sel[[s]]=fn.sel.simple(d=ktch.size.fq[[s]])
  }
  
  detach("package:TropFishR", unload=TRUE)  #remove after use because it causes conflicts with dplyr
  
}

#---Mean weight-based Mortality estimation------------------------------------------------------
if(do.mean.weight.based=="YES")
{
  library(fishmethods)
  
  # daily standardised mean weight of catch
  Mn.weit.ktch=list("Smooth hammerhead"=Smuz.hh.tdgdlf.size,
                    "Spinner shark"=Spinr.tdgdlf.size,
                    "Tiger shark"=Tiger.tdgdlf.size)
  
  # Beverton-Holt Nonequilibrium Z Estimator based on mean size of catch
  #note: based on Gedamke & Hoenig 2006. 
  #      Knife-edge selectivity. Only provides Z estimate for a period
  if(do.Gedamke_Hoenig=="YES")
  {
    #for dome-shape selectivity, select Lc as size of full vulnerability
    #   (Gedamke & Hoenig 2006 page 484)
    fn.bhnoneq=function(d,Lc,K,Linf,a,b)
    {
      d=d%>%mutate(year=as.numeric(substr(Finyear,1,4)),
                   mean=(Pred.mean/a)^(1/b))
      Z=bhnoneq(year=d$year,mlen=d$mean, ss=d$n,
                K=K,Linf=Linf,Lc=Lc,nbreaks=0,stZ=0.2,
                graph =F)
      return(Z)
    }
    store.z.bhnoneq=vector('list',length(Mn.weit.ktch))
    names(store.z.bhnoneq)=names(Mn.weit.ktch)
    for(i in 1:length(Mn.weit.ktch))
    {
      this=match(names(Mn.weit.ktch)[i],names(Size.sel))
      Lc=Size.sel[[this]]
      Lc=Lc[which.max(Lc$Sel),1]
      store.z.bhnoneq[[i]]=with(subset(LH.par,Name==names(Mn.weit.ktch)[i]),
                                fn.bhnoneq(d=Mn.weit.ktch[[i]],Lc=Lc,
                                           K=K,Linf=FL_inf*Mn.conv.Fl.Tl,  #convert TL_inf to FL_inf
                                           a=a_w8t,b=b_w8t))
    }
  }
}

#---Length-based Spawning potential ratio------------------------------------------------------ ACA
#note: based on Hordyk et al 2016. Assumptions:
#      asymptotic selectivity (overestimates F and underestimates SPR if dome-shape) so it
#      cannot be used because selectivity-at-age is normally distributed
if(do.length.based.SPR=="YES")
{
  library(LBSPR)
  LBSPR.assmnt=vector('list',length(ktch.size.fq))
  names(LBSPR.assmnt)=names(ktch.size.fq)
  apply.LBSPR.fn=function(LH,M,d,BinWidth,min.samp.size)
  {
    #get length at maturity from age at maturity
    k=LH$K
    Linf=LH$FL_inf*Mn.conv.Fl.Tl   #convert to FL as Linf is in TL 
    L50=Linf*(1-exp(-k*(LH$Age_50_Mat_min-LH$to)))
    L95=Linf*(1-exp(-k*(LH$Age_50_Mat_max-LH$to)))
    
    #populate life history object
    MyPars <- new("LB_pars")
    MyPars@Linf <- Linf
    MyPars@L50 <- L50 
    MyPars@L95 <- L95
    MyPars@MK <- M/k
    MyPars@Walpha=LH$a_w8t
    MyPars@Wbeta=LH$b_w8t
    MyPars@FecB=LH$b_w8t
    MyPars@L_units <- "cm"
    
    #populate lengths object
    tab=d%>%mutate(Size.class=BinWidth*floor(FL/BinWidth)+BinWidth/2)%>%
      group_by(year,Size.class)%>%
      summarise(n=n())%>%
      spread(year,n,fill=0)%>%
      data.frame
    LMids=tab$Size.class
    tab=tab[,-1]
    colnames(tab)=substr(colnames(tab),2,100)
    tab=tab[,colSums(tab) > min.samp.size] 
    
    MyLengths <- new("LB_lengths", LB_pars=MyPars)
    MyLengths@LMids<-LMids
    MyLengths@LData<-as.matrix(tab)
    MyLengths@L_units<- "cm"
    MyLengths@Years<-as.numeric(colnames(tab))
    MyLengths@NYears<-length(MyLengths@Years)
    
    #fit model
    Fit <- LBSPRfit(MyPars, MyLengths)
    
    return(Fit)
  }
  for(s in 1:length(LBSPR.assmnt))
  {
    M=mean(unlist(lapply(store.species.M[[match(names(Size.sel)[s],names(store.species.M))]],mean)))
    LBSPR.assmnt[[s]]=apply.LBSPR.fn(LH=LH.par%>%filter(SP.group==names(Size.sel)[s]),
                                     M=M,
                                     d=ktch.size.fq[[s]],
                                     BinWidth=10,
                                     min.samp.size=40)
  }
  
  for(s in 1:length(LBSPR.assmnt))
  {
    LBSPR.assmnt[[s]]@Ests
    plotMat(LBSPR.assmnt[[s]])
    plotEsts(LBSPR.assmnt[[s]])
  }
}

#---Single-species SPM -----------------------------------------------------------------------
if(Do.SPM=="YES")
{
  Estimable.qs=list(
    "copper shark"=NULL,
    "grey nurse shark"=c(NA,q2=.005,NA),
    "lemon shark"=NULL,
    "great hammerhead"=NULL,
    "milk shark"=c(q1=.005,NA,NA),
    "sawsharks"=NULL,
    "scalloped hammerhead"=c(q1=.005,NA,NA),
    "smooth hammerhead"=c(NA,q2=.005,q3=.001), 
    "spinner shark"=c(NA,q2=.005,q3=.001),
    "shortfin mako"=c(NA,q2=NA,q3=.001),    #CV for q2 too high for all years
    "spurdogs"=NULL,
    "tiger shark"=c(q1=.005,q2=.005,q3=.001),         
    "wobbegongs"=NULL) 
  
  
  #Loop over all species
  posfun=function(x,eps,pen)  #penalty for keeping biomass positive
  {
    if (x>=eps) return(x) else
    {
      pen=pen+.01*(x-eps)^2
      return (list(eps/(2-x/eps),pen))
    }
  }
  SPM=function(Init.propK,cpue,cpue.CV,Qs,Ktch,theta,HR_init,HR_init.sd,r.mean,r.sd,usePenalties) #population dynamics
  {
    #Population dynamics
    K=exp(theta[match('k',names(theta))])
    if(!fix.r=="YES")r=exp(theta[match('r',names(theta))])
    if(fix.r=="YES")r=r.mean
    if(estim.q=="YES")q=exp(theta[grepl("q",names(theta))])   
    if(estim.q=="NO")q=Qs
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
      if(Bt[t]<0) Bpen[t]=pen[[2]]
      Bt[t]=pen[[1]]
      
    }
    H=Ktch/Bt[-length(Bt)]
    
    #Loglikelihoods
        #Initial harvest rate
    HR.init.negLL=0
    #HR.init.negLL=-log(dnorm(H[1],HR_init,HR_init.sd))
    
        #r
    r.negLL=0
    if(!fix.r=="YES") r.negLL=-log(dnorm(r,r.mean,r.sd))*r.weight
    names(r.negLL)=NULL
    
    
        #cpue likelihood
    ln.cpue.hat=vector('list',length(cpue))
    ln.cpue.hat.full=ln.cpue=ln.cpue.hat
    cpue.negLL=rep(NA,length(cpue))
    CV.weight=cpue.negLL
    for(ku in 1:length(Qs))
    {
      if(!is.na(Qs[ku]))
      {
        no.cpue=which(!is.na(cpue[[ku]]))
        n.cpue=length(no.cpue)
        if(estim.q=="NO")q[match(names(Qs[ku]),names(q))]=exp(mean(log(cpue[[ku]][no.cpue]/Bt[no.cpue]),na.rm=T)) #Haddon 2001 page 321
        cpue.hat=q[match(names(Qs[ku]),names(q))]*Bt[-length(Bt)]
        
        ln.cpue.hat[[ku]]=log(cpue.hat)
        ln.cpue.hat.full[[ku]]=ln.cpue.hat[[ku]]
        ln.cpue[[ku]]=log(cpue[[ku]])
        ln.cpue.hat[[ku]]=ln.cpue.hat[[ku]][no.cpue]
        ln.cpue[[ku]]=ln.cpue[[ku]][no.cpue]
        cV=cpue.CV[[ku]][no.cpue]
        
        if(what.like=='full')
        {
          sqres=(ln.cpue[[ku]]-ln.cpue.hat[[ku]])^2
          cpue.negLL[ku]=(n.cpue/2)*(log(2*pi)+2*log(sqrt(sum(sqres,na.rm=T)/n.cpue))+1)
        }

        if(what.like=='kernel')
        {
          sigMa=exp(theta[match('sigma',names(theta))])
          cpue.negLL[ku]=-sum(dnorm(ln.cpue[[ku]],ln.cpue.hat[[ku]], sigMa, log = TRUE))
          
          CV.weight[ku]=sum(log(sqrt(sigMa^2+cV^2)),na.rm=T)
          
          cpue.negLL[ku]=cpue.negLL[ku]+CV.weight[ku]
        }
      }
    }
    cpue.negLL=sum(cpue.negLL,na.rm=T)
    
    #Total negloglike
    if(usePenalties)
    {
      negLL=HR.init.negLL+r.negLL+cpue.negLL+sum(Bpen,na.rm=T)
    }else
    {
      negLL=cpue.negLL
    }
      
    
    #Calculate MSY quantities
    Bmsy=K/2
    MSY=K*r/4
    Fmsy=r/2
    
    return(list(Bmsy=Bmsy,MSY=MSY,Fmsy=Fmsy,Bt=Bt,negLL=negLL,
                ln.cpue.hat=ln.cpue.hat,ln.cpue=ln.cpue,ln.cpue.hat.full=ln.cpue.hat.full))
  }
  fn.fill=function(x)   #fill in missing years function
  {
    aa=all.iers[which(!all.iers%in%x$yr)]
    aa1=x[1:length(aa),]
    aa1[,]=NA
    aa1$yr=aa
    x=rbind(x,aa1)%>%arrange(yr)
    return(x)
  }
  
  #Check max possible initial harvest rate
  Mx.init.harv=rep(NA,N.sp)
  names(Mx.init.harv)=Specs$SP.group
  for(s in 1: N.sp)
  {
    if(B.init==1)   
    {
      Mx.init.harv[s]=0
    }else
    {
      Id=match(Specs$SP.group[s],names(cpue.list))
      ct=Tot.ktch%>%filter(SP.group==Specs$SP.group[s])%>%
        group_by(finyear)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
      Mx.init.harv[s]=ct$LIVEWT.c[1]/max(ct$LIVEWT.c)
    }
  }
  
  #Estimate parameters
  Store.SPM=vector('list',N.sp)
  names(Store.SPM)=Specs$SP.group
  Store.stuff=Theta=Store.SPM
  for(s in 1: N.sp)
  {
    Id=match(Specs$SP.group[s],names(cpue.list))
    
    #catch
    ct=Tot.ktch%>%filter(SP.group==Specs$SP.group[s])%>%
      group_by(finyear)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
      data.frame
    all.iers=seq(min(ct$finyear),max(ct$finyear))
    
    if(length(which(!all.iers%in%ct$finyear))>0)
    {
      aa=all.iers[which(!all.iers%in%ct$finyear)]
      aa1=ct[1:length(aa),]
      aa1[,]=0
      aa1$finyear=aa
      ct=rbind(ct,aa1)%>%arrange(finyear)
    }
    
    
    #r prior
    Idd=match(Specs$SP.group[s],names(store.species))
    r.prior=store.species[[Idd]]$r.prior.normal$mean
    r.prior.sd=store.species[[Idd]]$r.prior.normal$sd
    
    Store.CPUE.eff.dummy=NULL
    QS.dummy=NULL
    CPUE.yr.dummy=NULL
    n.cpue.dummy=NULL
    CPUE.CV=NULL
    if(!is.null(cpue.list[[Id]]))
    {
      #Get cpues
      CPUE.1=cpue.list[[Id]]$Nat
      if(!is.null(CPUE.1))CPUE.1=CPUE.1%>%
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
          
          too.large=which(CPUE.CV[[ci]]>MAX.CV) #drop observation with too large CVs
          CPUE[[ci]][too.large]=NA
          CPUE.CV[[ci]][too.large]=NA
          
          CPUE.yr[[ci]]=x%>%filter(!is.na(MeAn))%>%
            dplyr::select(yr)
        }
      }
      
      #Initial values for estimable pars
      Mx.ktch=max(ct$LIVEWT.c,na.rm=T)
      K.init=k.times.mx.ktch*Mx.ktch
      r.init=Init.r[[match(Specs$SP.group[s],names(Init.r))]]
      QS=Estimable.qs[[match(Specs$SP.group[s],names(Estimable.qs))]]
      
      #loop over scenarios
      HR.o.scens=Mx.init.harv[match(Specs$SP.group[s],names(Mx.init.harv))]
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
          fn_ob=function(theta)SPM(Init.propK=B.init,
                                   cpue=CPUE.eff.scen,
                                   cpue.CV=CPUE.CV,
                                   Qs=QS,
                                   Ktch=ct$LIVEWT.c,
                                   theta,
                                   HR_init=HR_o,
                                   HR_init.sd=HR_o.sd,
                                   r.mean=r.prior,
                                   r.sd=r.prior.sd,
                                   usePenalties=usePen)$negLL
          #fit model
          if(estim.q=="YES")theta= c(k=log(K.init),log(QS[which(!is.na(QS))]))
          if(estim.q=="NO")theta= c(k=log(K.init))
          Lw.bound=log(c(Low.bound.K*Mx.ktch,rep(1e-6,length(theta)-1)))
          Up.bound=log(c(Up.bound.K*Mx.ktch,rep(1,length(theta)-1)))
          if(!fix.r=="YES")
          {
            theta=c(theta,r=log(r.init))
            Lw.bound=c(Lw.bound,log(0.01))
            Up.bound=c(Up.bound,log(0.75))
          }
          
          if(what.like=='kernel')
          {
            theta=c(theta,sigma=log(0.2))
            Lw.bound=c(Lw.bound,log(1e-2))
            Up.bound=c(Up.bound,log(1))
          }
          
          if(minimizer=='optim')
          {
            OptiM=optim(theta,fn_ob,method="L-BFGS-B",lower =Lw.bound,upper = Up.bound,hessian=T)
            if(Remove.bounds)
            {
              paramscale = magnitude(theta)
              OptiM=optim(jitter(OptiM$par),fn_ob,method="L-BFGS-B",hessian=T,
                          control = list(maxit = 1000, parscale = paramscale))
            }

            dummy.eff[[e]]=OptiM
          }
          if(minimizer=='nlminb')
          {
            nlmb <- nlminb(theta, fn_ob, gradient = NULL,lower =Lw.bound,upper = Up.bound)
            if(Remove.bounds) nlmb <- nlminb(jitter(nlmb$par), fn_ob, gradient = NULL)
            
            nlmb$hessian=hessian(fn_ob,nlmb$par)
            
            dummy.eff[[e]]=nlmb
          }
          Store.CPUE.eff[[e]]=CPUE.eff.scen
        }
        dummy[[h]]=dummy.eff
      }
      Store.SPM[[s]]=dummy
      Theta[[s]]=list(theta=theta,Lw.bound=Lw.bound,Up.bound=Up.bound)
      rm(HR.o.scens)
      Store.CPUE.eff.dummy=Store.CPUE.eff
      QS.dummy=QS
      CPUE.yr.dummy=CPUE.yr
      n.cpue.dummy=length(CPUE)
    }
    Store.stuff[[s]]=list(cpue=Store.CPUE.eff.dummy,CPUE.CV=CPUE.CV,Qs=QS.dummy,Ktch=ct$LIVEWT.c,
                          r.mean=r.prior,r.sd=r.prior.sd,yrs=all.iers,
                          cpue.yrs=CPUE.yr.dummy,n.cpues=n.cpue.dummy)
    
  }
  
  #Evaluate model at MLE
  SPM.preds=vector('list',length(Store.SPM))
  names(SPM.preds)=names(Store.SPM)
  for(s in 1: N.sp)
  {
    if(!is.null(Store.SPM[[s]]))
    {
      HR.o.scens=Mx.init.harv[s]
      dumy.pred=vector('list',length(HR.o.scens))
      names(dumy.pred)=paste("HR=",HR.o.scens)
      for(h in 1:length(HR.o.scens))
      {
        dummy.eff=vector('list',length(Efficien.scens))
        names(dummy.eff)=Efficien.scens
        for(e in 1:length(Efficien.scens))
        {
          dummy.eff[[e]]=SPM(Init.propK=B.init,
                             cpue=Store.stuff[[s]]$cpue[[e]],
                             cpue.CV=Store.stuff[[s]]$CPUE.CV,
                             Qs=Store.stuff[[s]]$Qs,
                             Ktch=Store.stuff[[s]]$Ktch,
                             theta=Store.SPM[[s]][[h]][[e]]$par,
                             HR_init=HR.o.scens[h],
                             HR_init.sd=HR_o.sd,
                             r.mean=Store.stuff[[s]]$r.mean,
                             r.sd=Store.stuff[[s]]$r.sd,
                             usePenalties=usePen)
        }
        dumy.pred[[h]]=dummy.eff
      }
      SPM.preds[[s]]=dumy.pred
      rm(HR.o.scens)
    }
  }
  
  #Get uncertainty  
  fn.un=function(fit,n)
  {
    SIGMA=solve(fit$hessian)	#getting the variance covariance matrix
    # std=sqrt(diag(SIGMA))		#Standard deviation of parameters
    #  R=SIGMA/(std%o%std)			#Parameter correlation, the V divided by the outer product of std
    return(rmvnorm(n,mean=fit$par,sigma=SIGMA))
  }
  SPM.preds_uncertainty=vector('list',length(Store.SPM))
  names(SPM.preds_uncertainty)=names(Store.SPM)
  Estim.par.samples=SPM.preds_uncertainty
  for(s in 1: N.sp)
  {
    if(!is.null(SPM.preds[[s]]))
    {
      #Draw random samples of estimable pars
      HR.o.scens=Mx.init.harv[s]
      dummy=vector('list',length(HR.o.scens))
      names(dummy)=HR.o.scens
      for(h in 1:length(HR.o.scens))
      {
        dummy1=vector('list',length(Efficien.scens))
        names(dummy1)=Efficien.scens
        for(e in 1:length(Efficien.scens))
        {
          dummy1[[e]]=fn.un(fit=Store.SPM[[s]][[h]][[e]],n=N.monte)
        }
        dummy[[h]]=dummy1
      }
      Estim.par.samples[[s]]=dummy
      
      #Use par sample to predict quantities of interest
      dumy.pred=vector('list',length(HR.o.scens))
      names(dumy.pred)=HR.o.scens
      for(h in 1:length(HR.o.scens))
      {
        dummy.eff=vector('list',length(Efficien.scens))
        names(dummy.eff)=Efficien.scens
        for(e in 1:length(Efficien.scens))
        {
          dum=vector('list',N.monte)
          for(x in 1:N.monte)
          {
            dum[[x]]=SPM(Init.propK=B.init,
                         cpue=Store.stuff[[s]]$cpue[[e]],
                         cpue.CV=Store.stuff[[s]]$CPUE.CV,
                         Qs=Store.stuff[[s]]$Qs,
                         Ktch=Store.stuff[[s]]$Ktch,
                         theta=Estim.par.samples[[s]][[h]][[e]][x,],
                         HR_init=HR.o.scens[h],
                         HR_init.sd=HR_o.sd,
                         r.mean=Store.stuff[[s]]$r.mean,
                         r.sd=Store.stuff[[s]]$r.sd,
                         usePenalties=usePen)
          }
          dummy.eff[[e]]=dum
        }
        dumy.pred[[h]]=dummy.eff
      }
      SPM.preds_uncertainty[[s]]=dumy.pred
      rm(HR.o.scens)
      print(paste("SPM  ",s,"--",names(store.species)[s]))
    }
  }
  
}

#---Catch-MSY ---------------------------------------------------------------
#note: uses SPM inputs (already done r prior, ct, etc)

rm(DATA.bio,DATA,DATA.ecosystems,Boat_bio,Scalefish,Boat_bio_header_sp)
if(Do.Ktch.MSY=="YES")
{
  # - Simulation testing CMSY
  if(Do.sim.test=="YES")
  {
    #1. create true population trajectories under different catch scenarios
    # operating model
    OM=function(K,r,Ktch) 
    {
      #Population dynamics
      Bt=rep(NA,N.yrs.sim.tst+1)
      Bpen=Bt
      Bt[1]=K
      for(t in 2:length(Bt)) Bt[t]=Bt[t-1]+r*Bt[t-1]*(1-Bt[t-1]/K)-Ktch[t-1]
      
      #Calculate MSY quantities
      Bmsy=K/2
      MSY=K*r/4
      Fmsy=r/2
      
      return(list(Bmsy=Bmsy,MSY=MSY,Fmsy=Fmsy,Bt=Bt))
    }
    
    # catch scenarios
    Yrs.sim.tst=1975:2018
    N.yrs.sim.tst=length(Yrs.sim.tst)
    rr=.1
    KK=1000  #in tonnes
    
    Ktchdummy=KK*.005+(1:(N.yrs.sim.tst/2))^1.2
    Ktch1=c(runif(N.yrs.sim.tst/2,.7,1.3)*Ktchdummy,
            rev(runif(N.yrs.sim.tst/2,.7,1.3)*Ktchdummy))
    Ktch1.low=Ktch1*.2
    
    # Ktch2=runif(N.yrs.sim.tst,.7,1.3)*KK*.005+(1:(N.yrs.sim.tst))^1.105
    # Ktch2.low=Ktch2*.2
    
    Ktch.sim=list(S1=Ktch1,S1.low=Ktch1.low)
    #Ktch.sim=list(S1=Ktch1,S1.low=Ktch1.low,S2=Ktch2,S2.low=Ktch2.low)
    
    OM.out=Ktch.sim
    for(l in 1:length(OM.out)) OM.out[[l]]=OM(K=KK,r=rr,Ktch=Ktch.sim[[l]])
    Mxylim=ceiling(max(unlist(Ktch.sim)))
    hndl.sim.test='C:\\Matias\\Analyses\\Population dynamics\\1.Other species\\2019\\Outputs\\Simulation_testing_catchMSY\\'
    names(OM.out)=c("High catch","Low catch")
    #names(OM.out)=c("S1","S1 low catch","S2","S2 low catch")
    fn.fig(paste(hndl.sim.test,"1.OM.outs",sep=''), 1800, 2400)
    par(mfcol=c(2,1),mar=c(2,2,1,1),oma=c(1.5,2,1,3),mgp=c(1,.5,0))
    for(l in 1:length(OM.out))
    {
      with(OM.out[[l]],
           {
             plot(Yrs.sim.tst,Bt[-length(Bt)]/KK,xlab='',main=names(OM.out)[l],
                  ylab='',type='l',lwd=2,ylim=c(0,1))
             par(new=T)
             plot(Yrs.sim.tst,Ktch.sim[[l]],ylab='',xlab='',ylim=c(0,Mxylim),
                  col=2,type='l',lwd=2,yaxt='n') 
             axis(4,seq(0,Mxylim,10),
                  seq(0,Mxylim,10))
           })
    }
    mtext('Years',1,outer=T,line=0,cex=1.5)
    mtext('Relative biomass',2,outer=T,las=3,cex=1.5)
    mtext('Catch (tonnes)',4,outer=T,line=1.5,col=2,cex=1.5)
    dev.off()
    
    
    #2. apply CMSY to simulated catches
    # Rprior=rnorm(1000,0.1,0.025)  #0.1 is mean of species studied
    # hist(Rprior)
    # fitdistr(Rprior, "gamma")
    # hist(rgamma(1000, shape=13, rate = 138))
    setwd(hndl.sim.test)
    Sim.test.out=list(Informative=Ktch.sim,Default=Ktch.sim)
    for(s in 1:length(Sim.test.out))
    {
      dummy=Ktch.sim
      for(l in 1:length(Ktch.sim))
      {
        k.lower=Low.bound.K
        k.upper=Up.bound.K
        if(names(Sim.test.out)[s]=='Informative')
        {
          Usr="Yes"
          StrtbiO=c(.98,.99)
          FinbiO=c(.2,.7)
          FinbiO.low=c(.5,.99)
        }
        if(names(Sim.test.out)[s]=='Default')
        {
          Usr="No"
          StrtbiO=c(.6,.99)
          FinbiO=c(.2,.7)
          FinbiO.low=c(.5,.99)
        }
        if(grepl(".low",names(Ktch.sim)[l]))
        {
          FinbiO=FinbiO.low
          k.lower=20
          k.upper=200
        }
        
        nm=paste(names(Sim.test.out)[s],names(Ktch.sim)[l])
        print(paste("-------------",nm))
        
        PATH=paste(hndl.sim.test,nm,sep='/')
        if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
        setwd(file.path(PATH))
        Path.ktch_msy=getwd()
        
        
        dummy[[l]]=Catch_MSY(ct=Ktch.sim[[l]],
                             yr=Yrs.sim.tst,
                             r.prior=unlist(list(shape=13,rate=138)),
                             user=Usr,
                             k.lower=k.lower,
                             k.upper=k.upper,
                             startbio=StrtbiO,
                             finalbio=FinbiO,
                             res="Very low",
                             n=SIMS,
                             sigR=ERROR,
                             ct.future=0,           
                             yr.future=1)
      }
      Sim.test.out[[s]]=dummy
    }
    
    #Plot relative biomass and MSY the different scenarios 
    fn.cons.po=function(low,up) c(low, tail(up, 1), rev(up), low[1])  #construct polygon
    Low.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (0+Nper)/100))   #get percentiles
    High.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (100-Nper)/100))
    
    
    fn.fig(paste(hndl.sim.test,"2.CatchMSY.performance",sep=''), 2400, 2400)
    par(mfcol=c(2,2),mar=c(2,2,1,1),oma=c(1.5,2,1,3),mgp=c(1,.5,0),las=1)
    names(Sim.test.out[[2]])=names(OM.out)
    for(s in 1:length(Sim.test.out))
    {
      for(l in 1:length(Ktch.sim))
      {
        Rel.bio=sweep(Sim.test.out[[s]][[l]]$bt, 2, Sim.test.out[[s]][[l]]$k, `/`)
        MED=apply(Rel.bio, 1, function(x) quantile(x, .5))
        
        
        Nper=(100-50)/2
        LOW.50=Low.percentile(Nper,Rel.bio)
        UP.50=High.percentile(Nper,Rel.bio)
        
        #100% of data
        Nper=(100-100)/2
        LOW.100=Low.percentile(Nper,Rel.bio)
        UP.100=High.percentile(Nper,Rel.bio)
        
        #construct polygons
        Year.Vec <-  fn.cons.po(Yrs.sim.tst,Yrs.sim.tst)
        Biom.Vec.50 <- fn.cons.po(LOW.50,UP.50)
        Biom.Vec.100 <- fn.cons.po(LOW.100,UP.100)
        
        plot(Yrs.sim.tst,MED,ylim=c(0,1),
             ylab="",xlab="",type='l',cex=0.7,pch=19,cex.axis=1)
        polygon(Year.Vec, Biom.Vec.100, col = rgb(.1,.1,.1,alpha=.15), border = "grey20")
        polygon(Year.Vec, Biom.Vec.50, col = rgb(.1,.1,.1,alpha=.3), border = "grey20")
        
        #true biomass
        lines(Yrs.sim.tst,OM.out[[l]]$Bt[-length(OM.out[[l]]$Bt)]/KK,lwd=2,col=2)
        
        if(s==2) mtext(names(Sim.test.out[[s]])[l],4,line=1.5,las=3)
        if(l==1) mtext(names(Sim.test.out)[s],3)
        
      }
      
    }
    legend('bottomleft',c('True trajectory','Predicted median'),bty='n',lty=1,
           col=c('red','black'),cex=1.25,lwd=2)
    
    mtext("Year",1,outer=T,cex=1.5)
    mtext("Relative biomass",2,outer=T,cex=1.5,las=3,line=0)
    dev.off()
  }
  
  # - Run CMSY for each species    takes 0.002 sec per iteration  
  system.time(for(s in 1: N.sp)
  {
 
    #Tested scenarios
    ktch_msy_scen=vector('list',length(SCENARIOS))
    names(ktch_msy_scen)=names(SCENARIOS)
    for(sc in 1:length(ktch_msy_scen))
    {
      if(is.na(SCENARIOS[[sc]]$R.prior[1]))
      {
        USR="No"
        ReS=RESILIENCE[[s]]
      }else
      {
        if(!is.na(SCENARIOS[[sc]]$R.prior[1]))USR= "Yes"
        ReS=NA
        Scen.start.bio=SCENARIOS[[sc]]$Initial.dep   
      }
    
      ktch_msy_scen[[sc]]=list(r.prior=SCENARIOS[[sc]]$R.prior,
                               user=USR,
                               k.lower=Low.bound.K,k.upper=Up.bound.K, 
                               startbio=Scen.start.bio,finalbio=FINALBIO,
                               res=ReS,
                               niter=SIMS,
                               sigR=SCENARIOS[[sc]]$Error)
        if(!is.na(ktch_msy_scen[[sc]]$r.prior[1])) ktch_msy_scen[[sc]]$r.prior=unlist(store.species[[s]]$r.prior)
      rm(USR,ReS)
    }
    
    #Execute Catch-MSY  
    print(paste("CMSY------","s=",s,names(store.species)[s]))
    
    PATH=paste(PHat,Specs$SP.group[s],sep='/')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
    setwd(file.path(PATH))
    Path.ktch_msy=getwd()
    Ktch_MSY=ktch_msy_scen
    
    Yrs=Store.stuff[[s]]$yrs   
    Tot.Ktch=Store.stuff[[s]]$Ktch
    yr.future=Current+(1:years.futures)
    ct.future=rep(mean(Tot.Ktch[(length(Tot.Ktch)-4):length(Tot.Ktch)]),years.futures)
    
    for(sc in 1:length(ktch_msy_scen))
    {
      Folder=names(ktch_msy_scen)[sc]
      print(paste("-----------------------------------","scenario=",Folder))
      if(!file.exists(paste(PATH,Folder,sep="/"))) dir.create(paste(PATH,Folder,sep="/"))
      setwd(paste(PATH,Folder,sep="/"))
      
      Ktch_MSY[[sc]]=Catch_MSY(ct=Tot.Ktch,
                               yr=Yrs,
                               r.prior=ktch_msy_scen[[sc]]$r.prior,
                               user=ktch_msy_scen[[sc]]$user,
                               k.lower=ktch_msy_scen[[sc]]$k.lower,
                               k.upper=ktch_msy_scen[[sc]]$k.upper,
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
    if(nrow(Tabl.scen.Ktch.MSY)>1)Tabl.scen.Ktch.MSY[1,'res']=""
    setwd(file.path(PATH))
    write.csv(Tabl.scen.Ktch.MSY,paste('Scenarios.csv',sep=''))
  })
}


#---Single-species age-structured aSPM (aSPM) -----------------------------------------------------------------------
#Haddon datalowSA package (https://rdrr.io/github/haddonm/datalowSA/f/vignettes/aspm.md)
if(Do.aSPM=="YES")
{
  #Tweak functions to allow 2 Qs for TDGDLF and 2 cpue series and Hessian
  dynamics= function (pars, infish, inglb, inprops) 
  {
    waa <- inprops$waa
    maa <- inprops$maa
    sela <- inprops$sela
    R0 <- exp(pars[1])
    B0 <- getB0(R0, inglb, inprops)
    if (length(pars) == 3)
    {
      dep <- doDepletion(pars[1], indepl = pars[3], inprops, 
                         inglb, inc = 0.02)
      spb <- SpB(dep$Ndepl, maa, waa)
      Rinit <- bh(spb, inglb$steep, R0, B0)
    }else
    {
      Rinit <- R0
    }
    nyrs <- length(infish[, "year"])
    nages <- inglb$nages
    maxage <- inglb$maxage
    Nt <- matrix(0, nrow = nages, ncol = (nyrs + 1), dimnames = list(0:(nages - 
                                                                          1), 0:nyrs))
    columns <- c("Year", "Catch", "PredC","SpawnB",
                 "ExploitB", "FullH","CPUE","CPUE2", 
                 "PredCE","PredCE2", "Deplete")
    fishery <- matrix(NA, nrow = (nyrs + 1), ncol = length(columns), 
                      dimnames = list(0:nyrs, columns))
    fishery[, "Year"] <- c((infish$year[1] - 1), infish$year)
    fishery[, "Catch"] <- c(NA, infish$catch)
    fishery[, "CPUE"] <- c(NA, infish$cpue)
    fishery[, "CPUE2"] <- c(NA, infish$cpue2)
    hS <- exp(-inglb$M/2)
    surv <- exp(-inglb$M)
    if (length(pars) == 3)
    {
      Nt[, 1] <- dep$Ndepl
    }else
    {
      Nt[, 1] <- Rinit
      for (age in 1:(maxage - 1)) Nt[age + 1, 1] <- Nt[age, 
                                                       1] * surv
      Nt[maxage + 1, 1] <- (Nt[maxage, 1] * surv)/(1 - surv)
    }
    for (yr in 2:(nyrs + 1)) {
      spb <- SpB(Nt[, (yr - 1)], maa, waa)
      exb <- ExB(Nt[, (yr - 1)] * hS, sela, waa)
      Nt[1, yr] <- bh(spb, inglb$steep, R0, B0)
      harvest <- min((fishery[yr, "Catch"]/exb), 0.85)
      hrate <- sela * harvest
      Ct <- (Nt[, (yr - 1)] * hS) * hrate
      Nt[2:nages, yr] <- ((Nt[1:(nages - 1), (yr - 1)] * hS) - 
                            Ct[1:(nages - 1)]) * hS
      Nt[nages, yr] <- Nt[nages, yr] + ((Nt[nages, yr - 1] * 
                                           hS) - Ct[nages]) * hS
      fishery[(yr - 1), 4:5] <- c(spb, exb)
      fishery[yr, c(3, 6)] <- c(sum(Ct * waa)/1000, hrate[nages])
    }
    spb <- SpB(Nt[, yr], maa, waa)
    exb <- ExB(Nt[, yr] * hS, sela, waa)
    fishery[yr, 4:5] <- c(spb, exb)
    fishery[, "Deplete"] <- fishery[, "SpawnB"]/B0
    ExpB <- fishery[1:nyrs, "ExploitB"]
    
    #catchability
      #fleet 1
    if(!exists(c('q1.yrs','q2.yrs')))
    {
      avq <- exp(mean(log(infish$cpue/fishery[1:nyrs, "ExploitB"]), 
                      na.rm = TRUE))
      fishery[2:(nyrs + 1), "PredCE"] <- ExpB * avq
    }
    if(exists(c('q1.yrs','q2.yrs')))
    {
      id=match(q1.yrs,infish$year)
      avq1 <- exp(mean(log(infish$cpue[id]/fishery[id+1, "ExploitB"]), 
                       na.rm = TRUE))
      fishery[id+1, "PredCE"] <- ExpB[id] * avq1
      
      id=match(q2.yrs,infish$year)
      avq2 <- exp(mean(log(infish$cpue[id]/fishery[id+1, "ExploitB"]), 
                       na.rm = TRUE))
      fishery[id+1, "PredCE"] <- ExpB[id] * avq2
    }
    
      #fleet 2
    avq.fleet2 <- exp(mean(log(infish$cpue2/fishery[1:nyrs, "ExploitB"]), 
                      na.rm = TRUE))
    fishery[2:(nyrs + 1), "PredCE2"] <- ExpB * avq.fleet2

    
    return(as.data.frame(fishery))
  }
  aspmLL= function (par, infish, inglb, inprops) 
  {
    fishery <- dynamics(par, infish, inglb, inprops)
    
    #catch
    penalty <- sum((fishery[, "Catch"] - fishery[, "PredC"])^2, 
                   na.rm = TRUE)/1000
    #cpue fleet 1
    pick <- which(fishery$CPUE > 0)
    LL <- -sum(dnorm(log(fishery[pick, "CPUE"]), 
                     log(fishery[pick,"PredCE"]), par[2], log = TRUE))
    #cpue fleet 2
    LL2=0
    pick <- which(fishery$CPUE2 > 0)
    if(length(pick)>0)LL2 <- -sum(dnorm(log(fishery[pick, "CPUE2"]), 
                     log(fishery[pick,"PredCE2"]), par[2], log = TRUE))
    
    return((LL+LL2+penalty))
  }
  fitASPM=  function (initpar, infish, inglb, inprops, callfun = aspmLL) 
  {
    paramscale = magnitude(initpar)
    bestL <- optim(initpar, callfun, method = "Nelder-Mead", 
                   infish = infish, inglb = inglb, inprops = inprops, 
                   control = list(maxit = 1000, parscale = paramscale))
    paramscale = magnitude(bestL$par)
    bestL <- optim(bestL$par, callfun, method="L-BFGS-B", 
                   infish = infish, inglb = inglb, inprops = inprops,hessian=T, 
                   control = list(maxit = 1000, parscale = paramscale))
    return(bestL)
  }
  
  
  #fill in objects required for aSPM
  for(l in 1:nrow(LH.par))   
  {
    #selectivity at age
    sel=Selectivity.at.age[[l]]$relative.sel
    sel1=Store.age.comp[[l]]$props$sela
    delta=length(sel1)-length(sel)
    if(delta>0)sel=c(sel,rep(sel[length(sel)],delta))
    if(delta<0)sel=sel[1:length(sel1)]
    Store.age.comp[[l]]$props$sela=sel
    
    #Mortality
    Store.age.comp[[l]]$glb$M=mean(unlist(lapply(store.species.M[[l]], mean)))
    
    #steepness
    Store.age.comp[[l]]$glb$steep=store.species.steepness[[l]]$mean
  }
  
  #init par values
  aSPM.init=list(
    "copper shark"=NULL,
    "grey nurse shark"=c(logR0=10,sigCE=0.3),
    "lemon shark"=NULL,
    "great hammerhead"=NULL,
    "milk shark"=NULL,
    "sawsharks"=NULL,
    "scalloped hammerhead"=c(logR0=10,sigCE=0.3),
    "smooth hammerhead"=c(logR0=10,sigCE=0.3), 
    "spinner shark"=c(logR0=10,sigCE=0.3),
    "shortfin mako"=c(logR0=10,sigCE=0.3),
    "spurdogs"=NULL,
    "tiger shark"=c(logR0=11,sigCE=0.3),         
    "wobbegongs"=NULL)    
  
  #run aspm
  Store.aSPM=vector('list',N.sp)
  names(Store.aSPM)=names(aSPM.init)
  
  for(l in 1:N.sp)    #Still issues here with l==9 ACA!!
  {
    if(!is.null(cpue.list[[l]]))
    {
      #catch
      fish=Tot.ktch%>%
        filter(SP.group==names(cpue.list)[l])%>%
        group_by(finyear)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
        data.frame%>%
        rename(year=finyear,
               catch=LIVEWT.c)
      
      #cpue
      Id=match(names(cpue.list)[l],names(cpue.list))
      CPUE=compact(cpue.list[[Id]])
      len.cpue=length(CPUE)
      if(len.cpue>1)
      {
        if("Nat"%in%names(CPUE))
        {
          a=CPUE[-match("Nat",names(CPUE))]
          
          CPUE2=CPUE[[match("Nat",names(CPUE))]]%>%    
            dplyr::select(yr,MeAn,CV)%>%
            rename(year=yr,cpue2=MeAn,se2=CV)
          
          iid=which(CPUE2$se2>MAX.CV)
          CPUE2$cpue2[iid]=NA
          CPUE2$se2[iid]=NA
          
          q1.yrs=as.numeric(substr(a$TDGDLF.mon$Finyear,1,4))
          q2.yrs=as.numeric(substr(a$TDGDLF.day$Finyear,1,4))
          CPUE=do.call(rbind,a)%>%
            mutate(yr=as.numeric(substr(Finyear,1,4)),
                   MeAn=Mean)
        }else
        {
          a=CPUE
          q1.yrs=as.numeric(substr(a$TDGDLF.mon$Finyear,1,4))
          q2.yrs=as.numeric(substr(a$TDGDLF.day$Finyear,1,4))
          CPUE=do.call(rbind,a)%>%
            mutate(yr=as.numeric(substr(Finyear,1,4)),
                   MeAn=Mean)
          iid=which(CPUE$CV>MAX.CV)
          CPUE$MeAn[iid]=NA
          CPUE$CV[iid]=NA
        }

      }
      if(len.cpue==1)
      {
        CPUE=CPUE[[1]]
        colnames(CPUE)[2]="MeAn"
        if(!'yr'%in%names(CPUE))
        {
          CPUE=CPUE%>%
            mutate(yr=as.numeric(substr(Finyear,1,4)))
        }
        iid=which(CPUE$CV>MAX.CV)
        CPUE$MeAn[iid]=NA
        CPUE$CV[iid]=NA
      }
      
      CPUE=CPUE%>%    
          dplyr::select(yr,MeAn,CV)%>%
          rename(year=yr,cpue=MeAn,se=CV)
      if(min(CPUE$se,na.rm=T)>1)  CPUE$se=CPUE$se/100
      
      fish=fish%>%left_join(CPUE,by='year')
      fish0.ktch=fish[1:20,]%>%
              mutate(year=year-20,
                     catch=0,
                     cpue=NA,
                     se=NA)
      fish=rbind(fish0.ktch,fish)
      iid=which(fish$se>MAX.CV)
      fish$cpue[iid]=NA
      fish$se[iid]=NA
      
      if(exists(c('q1.yrs','q2.yrs')))
      {
        ia=fish$year[!is.na(fish$cpue)]
        q1.yrs=q1.yrs[which(q1.yrs%in%ia)]
        q2.yrs=q2.yrs[which(q2.yrs%in%ia)]
      }
      
      if(exists('CPUE2')) fish=fish%>%left_join(CPUE2,by='year')
      if(!exists('CPUE2')) fish=fish%>%mutate(cpue2=NA,se2=NA)

      
      
      #init par values
      Id=match(Specs$SPECIES[match(names(aSPM.init)[l],Specs$Name)],names(Store.age.comp))
      
      #global pars
      glb=Store.age.comp[[Id]]$glb
      glb$nyrs=nrow(fish)
      
      #relations at age
      props=Store.age.comp[[Id]]$props
      
      #fit model
      pars=aSPM.init[[l]]
      ans <- fitASPM(initpar=pars,infish=fish,inglb=glb,inprops=props)
      fishery <- dynamics(ans$par,infish=fish,inglb=glb,inprops = props)
      
      #get uncertainty
      Par.samp=fn.un(fit=ans,n=N.monte)
      fishery.MC=vector('list',nrow(Par.samp))
      for(f in 1:nrow(Par.samp))
      {
        fishery.MC[[f]]=dynamics(Par.samp[f,],infish=fish,inglb=glb,inprops = props)
      }
      
      if(exists(c('q1.yrs','q2.yrs'))) rm(q1.yrs,q2.yrs)
      if(exists('CPUE2')) rm(CPUE2)

      Store.aSPM[[l]]=list(fit=ans,quantities=fishery,fish=fish,quantities.MC=fishery.MC)
      
      print(paste("aSPM--------- l=",l,names(Store.aSPM)[l]))
    }
  }
}


#---RESULTS SECTION------

setwd(paste(hNdl,'/Outputs',sep=''))

#Plot overall catch spatial distribution 
plot.spatial.dist=FALSE
if(plot.spatial.dist)
{
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
    legend(111,-9.75,NMs,bty='n',cex=.925,xjust=0)
    Lg=round(quantile(D.agg$LIVEWT.c,probs=c(.75,.95,1)),0)
    legend('right',paste(Lg),pch=21,pt.bg="grey50",bty='n',pt.cex=Lg/scaler,title="Tonnes",cex=1.1)
  }
  fn.fig("Figure 1_Map", 1200, 2400)
  smart.par(n.plots=length(Map.sp),MAR=c(.1,.1,.1,.1),OMA=c(2.5,2.5,1.5,.1),MGP=c(1,.5,0))
  for(s in 1: length(Map.sp))
  {
    NMs=Map.sp[s]
    if(NMs=="Low") NMs="Low resilience"
    if(NMs=="Very.low") NMs="Very low resilience"
    
    fn.sptial.ktch(d=subset(Map.this,Name==Map.sp[s]),NMs=NMs)
    if(s%in%8:10)axis(side = 1, at =seq(xlm[1],xlm[2],4), labels = seq(xlm[1],xlm[2],4), tcl = .5,las=1,cex.axis=1)
    if(s%in%seq(1,10,3))axis(side = 2, at = seq(ylm[1],ylm[2],4), labels = -seq(ylm[1],ylm[2],4),tcl = .5,las=2,cex.axis=1)
  }
  mtext(expression(paste("Latitude ",degree,"S")),side=2,line=0.85,las=3,cex=1.2,outer=T)
  mtext(expression(paste("Longitude ",degree,"E")),side=1,line=1.1,cex=1.2,outer=T)
  dev.off()

  #Catch by year
  # South.WA.lat=c(-36,-25); South.WA.long=c(112,130)
  # Long.seq=seq(South.WA.long[1]+1,South.WA.long[2]-1,by=3)
  # Lat.seq=c(-26,-28,-30,-32,-34)
  # 
  # PLATE=c(.01,.9,.075,.9)
  # numInt=20
  # Colfunc <- colorRampPalette(c("yellow","red"))
  # Couleurs=c("white",Colfunc(numInt-1))
  # 
  # fn.ctch.plot.all.yrs=function(DATA,tcl.1,tcl.2,numInt) 
  # {
  #   DATA$LAT=-as.numeric(substr(DATA$BLOCKX,1,2))
  #   DATA$LONG=100+as.numeric(substr(DATA$BLOCKX,3,4))
  #   DATA=subset(DATA,!is.na(LAT))
  #   DATA=subset(DATA,LAT>(-40))
  #   DATA$blk=substr(DATA$BLOCKX,1,4)
  #   A=aggregate(LIVEWT.c~FINYEAR+blk,DATA,sum)
  #   Ymax=max(A$LIVEWT.c)
  #   Ymin=min(A$LIVEWT.c)
  #   
  #   Breaks=quantile(A$LIVEWT.c,probs=seq(0,1,1/numInt),na.rm=T)
  #   a=South.WA.long[1]:South.WA.long[2]
  #   b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
  #   DATA$BLOCKX.c=with(DATA,paste(LAT,LONG))
  #   
  #   FINYrS=table(DATA$FINYEAR)
  #   FINYrS=FINYrS[FINYrS>20]
  #   FINYrS=sort(names(FINYrS))
  #   
  #   smart.par(length(FINYrS)+1,MAR=c(1,2.5,1.5,.1),OMA=c(2,2,.1,.1),MGP=c(.1,.7,0))
  #   for(y in 1:length(FINYrS))
  #   {
  #     A=subset(DATA,FINYEAR==FINYrS[y])
  #     MapCatch=with(A,aggregate(LIVEWT.c,list(BLOCKX.c),FUN=sum,na.rm=T))
  #     colnames(MapCatch)=c("BLOCKX.c","Total Catch")
  #     id=unique(match(MapCatch$BLOCKX.c,DATA$BLOCKX.c))
  #     MapCatch$LAT=DATA$LAT[id]
  #     MapCatch$LONG=DATA$LONG[id]
  #     msn.lat=seq(min(MapCatch$LAT),max(MapCatch$LAT))
  #     msn.lat=msn.lat[which(!msn.lat%in%MapCatch$LAT)]
  #     if(length(msn.lat)>0)
  #     {
  #       dummy=MapCatch[1:length(msn.lat),]
  #       dummy$`Total Catch`=0
  #       dummy$LAT=msn.lat
  #       dummy$BLOCKX.c=with(dummy,paste(LAT,LONG))
  #       MapCatch=rbind(MapCatch,dummy)
  #     }
  #     
  #     MapCatch$LAT.cen=MapCatch$LAT-.5
  #     MapCatch$LONG.cen=MapCatch$LONG+.5  
  #     MapCatch=MapCatch[order(MapCatch$LAT.cen,MapCatch$LONG.cen),]
  #     MapCatch=subset(MapCatch,LONG.cen<=South.WA.long[2])
  #     long=sort(unique(MapCatch$LONG.cen))
  #     lat=sort(unique(MapCatch$LAT.cen))      #latitude vector for image  
  #     MapCatch=MapCatch[,match(c("LONG.cen","LAT.cen","Total Catch"),names(MapCatch))]  
  #     Reshaped=as.matrix(reshape(MapCatch,idvar="LONG.cen",  	#transposed as matrix 	
  #                                timevar="LAT.cen",v.names="Total Catch", direction="wide"))	
  #     Reshaped=Reshaped[order(Reshaped[,1]),]
  #     Reshaped=Reshaped[,-1]	
  #     numberLab=10
  #     colLeg=(rep(c("black",rep("transparent",numberLab-1)),(numInt+1)/numberLab))
  #     
  #     plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
  #     image(long,lat,z=Reshaped,xlab="",ylab="",col =Couleurs,breaks=Breaks,axes = FALSE,add=T)			
  #     axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tcl = tcl.1)
  #     axis(side = 4, at = South.WA.lat[2]:South.WA.lat[1], labels = F,tcl =tcl.2)
  #     par(new=T)
  #     plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
  #     legend('top',FINYrS[y],bty='n',cex=1.2)
  #     axis(side = 1, at =Long.seq, labels = Long.seq, tcl = .35,las=1,cex.axis=1,padj=-.15)
  #     axis(side = 2, at = Lat.seq, labels = -Lat.seq,tcl = .35,las=2,cex.axis=1,hadj=1.1)
  #   }
  #   plot(a,b,ann=F,axes=F,col='transparent')
  #   color.legend(quantile(a,probs=.25),quantile(b,probs=.75),quantile(a,probs=.4),quantile(b,probs=.25),
  #                paste(round(Breaks,0),"kg"),rect.col=Couleurs,gradient="y",col=colLeg,cex=.5)
  # }
  # Spatial.ktch.sp=c(19000,18023,20000,18022,13000)
  # names(Spatial.ktch.sp)=c("Smooth hammerhead","Spinner shark","Spurdogs",
  #                   "Tiger shark","Wobbegongs")
  # pdf("Appendix1_Spatial catch by year.pdf")
  # for(s in 1: length(Spatial.ktch.sp))
  # {
  #   fn.ctch.plot.all.yrs(DATA=subset(Data.monthly,SPECIES==Spatial.ktch.sp[s]),tcl.1=.1,tcl.2=.1,numInt=20)
  #   legend("bottom",names(Spatial.ktch.sp)[s],bty='n',cex=.9)
  #   mtext("Longitude (E)",1,outer=T,line=1,cex=1.25)
  #   mtext("Latitude (S)",2,outer=T,line=0,cex=1.25,las=3)
  # }
  # dev.off()
}

#Plot total catch and effort
ktch.s=subset(Data.monthly,SPECIES%in%c(19000,Specs$SPECIES))%>%
  mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
  group_by(finyear)%>%
  summarise(Tot=sum(LIVEWT.c/1000,na.rm=T))%>%
  dplyr::select(finyear,Tot)
ktch.n=subset(Data.monthly.north,SPECIES%in%c(19000,Specs$SPECIES))%>%
  mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
  group_by(finyear)%>%
  summarise(Tot=sum(LIVEWT.c/1000,na.rm=T))%>%
  dplyr::select(finyear,Tot)

Effrt.s=Effort.monthly%>%mutate(finyear=as.numeric(substr(FINYEAR,1,4)))
Effrt.n=Effort.monthly.north%>%mutate(finyear=as.numeric(substr(FINYEAR,1,4)))

all.YYrs=seq(min(ktch.s$finyear),max(ktch.s$finyear))
if(length(which(!all.YYrs%in%ktch.n$finyear))>0)
{
  aa=all.YYrs[which(!all.YYrs%in%ktch.n$finyear)]
  aa1=ktch.n[1:length(aa),]
  aa1[,]=0
  aa1$finyear=aa
  ktch.n=rbind(ktch.n,aa1)%>%arrange(finyear)
}
if(length(which(!all.YYrs%in%Effrt.n$finyear))>0)
{
  aa1=Effrt.n[1:length(aa),]
  aa1[,]=0
  aa1$finyear=aa
  Effrt.n=rbind(Effrt.n,aa1)%>%arrange(finyear)
  
}

fn.fig("Figure 4. Total Catch of analysed species and effort time series", 1800, 2400)
par(mfcol=c(2,1),mar=c(1.5,2.5,1.5,.5),oma=c(2,2,.1,4),las=1,mgp=c(1,.6,0))
  #North
plot(all.YYrs,ktch.n$Tot,type='o',pch=19,col=1,cex=.75,ylab="",xlab="",main="North",
     ylim=c(0,max(c(ktch.s$Tot,ktch.n$Tot))))
par(new=T)
plot(Effrt.n$finyear,Effrt.n$Hook.days,type='l',col="grey55",xlab="",ylab="",axes=F,lwd=2.5,lty=3)
axis(side = 4)
mtext(side = 4, line = 3, 'Total effort (hook days)',las=3,cex=1.5)

legend("topleft",c("Catch","Effort"),bty='n',lty=c(1,3),col=c("black","grey55"),lwd=2.5)

  #South
plot(all.YYrs,ktch.s$Tot,type='o',pch=19,col=1,cex=.75,ylab="",xlab="",main="South",
     ylim=c(0,max(c(ktch.s$Tot,ktch.n$Tot))))
par(new=T)
plot(Effrt.s$finyear,Effrt.s$Total,type='l',col="grey55",xlab="",ylab="",axes=F,lwd=2.5,lty=3)
axis(side = 4)
mtext(side = 4, line = 3, 'Total effort (km gn days)',las=3,cex=1.5)

mtext("Year",1, line = .5,outer=T,cex=1.5)
mtext("Total catch (tonnes)",2,outer=T,las=3,cex=1.5)
dev.off()


#---Spatio-temporal catch------  
#note: bubble size is proportion of blocks fished out of maximum number of blocks fished for each species
fn.spatio.temp.catch.dist=function(d)
{
  d1=d%>% filter(SNAME%in%Keep.species)%>%
    count(FINYEAR,SPECIES,BLOCKX)%>%
    group_by(FINYEAR,SPECIES)%>%
    mutate(n=ifelse(n>0,1,0))%>%
    group_by(FINYEAR,SPECIES)%>%
    summarise(n=sum(n,na.rm=T))%>%
    spread(FINYEAR,n,fill=0)
  All.sp=d1$SPECIES 
  d1=d1%>%dplyr::select(-SPECIES)%>%as.matrix
  Mx=apply(d1,1,max)
  d1=d1/Mx
  
  yrs=as.numeric(substr(colnames(d1),1,4))
  Sp.nms=subset(All.species.names,SPECIES%in%All.sp)
  
  plot(1:nrow(d1),1:nrow(d1),col='transparent',ylab="",xlab="",yaxt='n',xlim=c(min(yrs),max(yrs)))
  for(p in 1:length(All.sp)) points(yrs,rep(p,length(yrs)),col='steelblue',cex=2*d1[p,],pch=19)
  axis(2,1:length(All.sp),capitalize(Sp.nms$SNAME),las=1)
  mtext(side = 1, line = 2, 'Financial year',cex=1.5)
  
  par(new=T)
  plot(yrs,Effort_blocks$Tot,type='l',col=rgb(0,0,1,alpha=.3),xlab="",ylab="",axes=F,lwd=5,lty=1)
  axis(side = 4,las=1)
  mtext(side = 4, line = 3, 'Number of blocks fished',las=3,cex=1.5)
  
  
}
fn.fig("Figure Spatio-temporal catch", 2400, 2400)
par(mar=c(2.5,4,.1,1),oma=c(.5,6,.1,3),mgp=c(1.5,.7,0))
fn.spatio.temp.catch.dist(d=rbind(Data.monthly%>%
                                        filter(!is.na(BLOCKX))%>%
                                        dplyr::select(SNAME,FINYEAR,SPECIES,BLOCKX),
                                  Data.monthly.north%>%
                                    filter(!is.na(BLOCKX))%>%
                                    dplyr::select(SNAME,FINYEAR,SPECIES,BLOCKX)))
dev.off()


#---Mean weight-based Mortality estimation RESULTS----
if(do.mean.weight.based=="YES")
{
  fn.plt.mn.ktch.wght=function(d)
  {
    yr=as.numeric(substr(d$Finyear,1,4))
    plot(yr,d$Pred.mean,ylab="",xlab="",pch=19,cex=1.5,cex.axis=1.25,
         ylim=c(min(d$Pred.mean-1.96*d$Pred.SE),max(d$Pred.mean+1.96*d$Pred.SE)))
    arrows(yr,d$Pred.mean,yr,d$Pred.mean-1.96*d$Pred.SE, angle=90, length=0.05)
    arrows(yr,d$Pred.mean,yr,d$Pred.mean+1.96*d$Pred.SE, angle=90, length=0.05)
  }
  tiff(file="pred Mean weight of ktch.tiff",width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(length(Mn.weit.ktch),MAR=c(1,3,3,.1),OMA=c(3,1,.1,.1),MGP=c(.1,.7,0))
  for(i in 1:length(Mn.weit.ktch))
  {
    fn.plt.mn.ktch.wght(d=Mn.weit.ktch[[i]])
    mtext(names(Mn.weit.ktch)[i],cex=1.5)
  }
  mtext("Mean catch weight (kg)",2,outer=T,cex=1.25,las=3,line=-1)
  mtext("Financial year",1,outer=T,cex=1.25,line=1.5)
  dev.off() 
}


#---SPM RESULTS------
if(Do.SPM=="YES")
{
  fn.cons.po=function(low,up) c(low, tail(up, 1), rev(up), low[1])  #construct polygon
  
  
  #Plot obs VS pred cpues  
  fn.plt.cpue=function(ob,ob.CV,pred,Convergence,pred.fit)
  {
    iid=match(pred.fit,pred)
    Yr=all.iers
    ob.CV=ob.CV[iid]
    plot(Yr,pred,pch=19,cex=1.5,ylab="",xlab="",
         ylim=c(min(c(ob-ob.CV,pred),na.rm=T),max(c(ob+ob.CV,pred),na.rm=T)))
    points(Yr[iid],ob,col="orange",pch=19,cex=1.25)
    segments(Yr[iid],ob,Yr[iid],ob+ob.CV,col="orange")
    segments(Yr[iid],ob,Yr[iid],ob-ob.CV,col="orange")
    legend("bottomright",paste("convergence=",Convergence),bty='n')
  }
  Paz=paste(hNdl,"/Outputs/SPM.fit/",sep="")
  if(!file.exists(file.path(Paz))) dir.create(file.path(Paz))
  for(s in 1: N.sp)
  {
    if(!is.null(SPM.preds[[s]]))
    {
      fn.fig(paste(Paz,names(SPM.preds)[s],sep=""),2400,2400)
      
      HR.o.scens=Mx.init.harv[s]
      nrw=length(HR.o.scens)*length(Efficien.scens)
      ncl=SPM.preds[[s]][[1]][[1]]$ln.cpue
      ncl=length(ncl[!sapply(ncl,is.null)])
      par(mfrow=c(nrw,ncl),mar=c(1.2,2,.2,.1),oma=c(1.5,1.75,1.5,1),las=1,cex.axis=.8,mgp=c(1,.42,0))
      for(h in 1:length(HR.o.scens))
      {
        for(e in 1:length(Efficien.scens))
        {
          for(x in 1:Store.stuff[[s]]$n.cpues)
            if(!is.null(SPM.preds[[s]][[h]][[e]]$ln.cpue[[x]]))
            {
              fn.plt.cpue(ob=SPM.preds[[s]][[h]][[e]]$ln.cpue[[x]],
                          ob.CV=Store.stuff[[s]]$CPUE.CV[[x]],
                          pred=SPM.preds[[s]][[h]][[e]]$ln.cpue.hat.full[[x]],
                          Convergence=Store.SPM[[s]][[h]][[e]]$convergence,
                          pred.fit=SPM.preds[[s]][[h]][[e]]$ln.cpue.hat[[x]])
              Crip=Efficien.scens[e]
              if(x==1) Crip=0
              par(font=2)
              #  legend("bottomleft",paste("HR=",HR.o.scens[h],", Creep=",Crip,
              #                            ", CPUE.",x,sep=""),bty='n',cex=.85)
            }
        }
      }
      legend("bottomleft",c("observed","predicted"),pch=19,cex=1.25,col=c("orange","black"),bty='n')
      mtext(names(SPM.preds)[s],3,outer=T)
      mtext("year",1,outer=T)
      mtext("lncpue",2,outer=T,las=3)
      rm(HR.o.scens)
      
      dev.off()
    }
  }  
  
  
  #probability of above and below reference points     
  add.probs=function(id.yr,YR,DAT,UP,LOW,SRT,CEX)
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
      segments(YR[id.yr],B.target,YR[id.yr],1,col=CL.ref.pt[1],lwd=8,lend="butt")  
      Legn=round(100*P.above.target)
      if(Legn==0)Legn="<1"
      text(YR[id.yr],mean(c(B.target,UP[id.yr])),paste(Legn,"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
    if(P.between.thre.tar>0)
    {
      Upseg=B.target
      Lwseg=B.threshold
      segments(YR[id.yr],Upseg,YR[id.yr],Lwseg,col=CL.ref.pt[2],lwd=8,lend="butt")  
      Legn=round(100*P.between.thre.tar)
      if(Legn==0)Legn="<1"
      text(YR[id.yr],mean(c(Upseg,Lwseg))*1.025,paste(Legn,"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
    if(P.between.lim.thre>0)
    {
      Upseg=B.threshold
      Lowseg=B.limit
      segments(YR[id.yr],Upseg,YR[id.yr],Lowseg,col=CL.ref.pt[3],lwd=8,lend="butt")
      Legn=round(100*P.between.lim.thre)
      if(Legn==0)Legn="<1"
      wher.txt=mean(c(Upseg,Lowseg))*1.025
      if(wher.txt>0.5) wher.txt=0.5*.9
      text(YR[id.yr],wher.txt,paste(Legn,"%",sep=""),
           col="black",cex=CEX,srt=SRT,font=2,pos=2)
    }
    if(P.below.limit>0)
    {
      segments(YR[id.yr],B.limit,YR[id.yr],0,col=CL.ref.pt[4],lwd=8,lend="butt")
      Legn=round(100*P.below.limit)
      if(Legn==0)Legn="<1"
      
      text(YR[id.yr],B.limit*0.85,paste(Legn,"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
    
    return(list(C1=P.above.target, C2=P.between.thre.tar,
                C3=P.between.lim.thre, C4=P.below.limit))
  }
  
  #Plot biomass  
  Ktch.CL=rgb(0.1,0.1,0.8,alpha=0.4)
  Col.RP=c("red","orange","forestgreen")
  fn.plt.bio.ktch=function(Yr,Bt,Bmsy,Ktch,CX.AX,CX,DAT,Add.ktch,LOW,HIGH)
  {
    Year.Vec <-  fn.cons.po(Yr,Yr)
    Biom.Vec <- fn.cons.po(Bt[match(LOW,row.names(Bt)),],
                           Bt[match(HIGH,row.names(Bt)),]) 
    plot(Yr,Bt[match("50%",row.names(Bt)),],col="black",ylim=c(0,1.05),
         type='l',lwd=1.5,xaxt='n',cex=0.3,pch=19,ylab="",xlab="",xaxs="i",yaxs="i")
    polygon(Year.Vec, Biom.Vec, col = rgb(.1,.1,.1,alpha=.1), border = "grey70")
    abline(h=B.threshold,col=Col.RP[2],lwd=1.15)               #threshold  
    abline(h=B.limit,col=Col.RP[1],lwd=1.15)         #limit         
    abline(h=B.target,col=Col.RP[3],lwd=1.15)       #target
    axis(1,at=Yr,labels=F,tck=-0.015)
    axis(1,at=seq(Yr[1],Yr[length(Yr)],5),labels=seq(Yr[1],Yr[length(Yr)],5),
         tck=-0.03,cex.axis=CX.AX)
    
    #add probs
    Store.probs=add.probs(id.yr=match(Current,Yr),YR=Yr,DAT=DAT,
                          UP=Bt[match(HIGH,row.names(Bt)),],
                          LOW=Bt[match(LOW,row.names(Bt)),],
                          SRT=0,CEX=CX)
    
    #add catch
    if(Add.ktch=="YES")
    {
      par(new=T)
      plot(Yr,Ktch,type='l',col=Ktch.CL,xlab="",ylab="",axes=F,lwd=1.5,xaxs="i")
      axis(side = 4)
    }

    
    return(Store.probs)
  }
  
  #1. Get median and percentiles
  Med.biom=vector('list',length(Store.SPM))
  names(Med.biom)=names(Store.SPM)
  for(s in 1: N.sp)
  {
    if(!is.null(SPM.preds[[s]]))
    {
      HR.o.scens=Mx.init.harv[s]
      dummy1=vector('list',length(HR.o.scens))
      names(dummy1)=HR.o.scens
      for(h in 1:length(HR.o.scens))
      {
        dummy2=vector('list',length(Efficien.scens))
        names(dummy2)=Efficien.scens
        for(e in 1:length(Efficien.scens))
        {
          dummy=subListExtract(SPM.preds_uncertainty[[s]][[h]][[e]],"Bt")
          dummy=sweep(do.call(rbind,dummy), 1, exp(Estim.par.samples[[s]][[h]][[e]][,1]), `/`)
          idd=which(round(dummy[,ncol(dummy)],1)<0.05) #remove cases where final biomass=0 
          if(length(idd)) dummy=dummy[-idd,]
          Bt.all=dummy
          Bt=apply(dummy,2,function(x) quantile(x,probs=c(0,0.5,1)))   #100% to make it comparable to CMSY
          #Bt=apply(dummy,2,function(x) quantile(x,probs=c(0.2,0.5,0.8)))   #60% as required for MSC
          Bt=Bt[,-ncol(Bt)] #remove future Bt
          dummy=subListExtract(SPM.preds_uncertainty[[s]][[h]][[e]],"Bmsy")    
          dummy=do.call(rbind,dummy)
          if(length(idd)) dummy=as.matrix(dummy[-idd,])
          Bmsy=apply(dummy,2,function(x) quantile(x,probs=c(.025,.5,.975)))
          colnames(Bmsy)="Bmsy"
          
          dummy2[[e]]=list(Bt=Bt,Bmsy=Bmsy,Bt.all=Bt.all)  
        }
        dummy1[[h]]=dummy2
      }
      Med.biom[[s]]=dummy1
      rm(HR.o.scens)
    }
  }
  
  #2. Plot
  do.col="NO"
  if(do.col=="NO") colfunc <- colorRampPalette(c("grey95","grey60"))
  if(do.col=="YES") colfunc <- colorRampPalette(c("aliceblue","lightblue3"))
  if(do.col=="NO") CL.ref.pt=c("black","grey35","grey55","grey80")
  if(do.col=="YES") CL.ref.pt=c('forestgreen','yellow','orange','red')
  
  Store.cons.Like.SPM=vector('list',N.sp)
  names(Store.cons.Like.SPM)=Specs$SP.group
  Store.cons.Like.SRM=Store.cons.Like.SPM
  
  fn.fig(paste(hNdl,"/Outputs/Figure 3_Biomass_SPM",sep=""),2400,2400)
  HR.o.scens=Mx.init.harv[1]
  smart.par(n.plots=length(compact(SPM.preds)),MAR=c(1.2,2,1.5,1.25),
            OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
  par(las=1,cex.axis=.8)
  for(s in 1: N.sp)
  {
    if(!is.null(SPM.preds[[s]]))
    {
      HR.o.scens=Mx.init.harv[s]
      for(h in 1:length(HR.o.scens))
      {
        nne=length(Efficien.scens)
        for(e in 1:nne)
        {
          if(!is.null(SPM.preds[[s]][[h]][[e]]))
          {
            DAT=t(Med.biom[[s]][[h]][[e]]$Bt.all)
            Store.cons.Like.SPM[[s]]=fn.plt.bio.ktch(Yr=Store.stuff[[s]]$yrs,
                                                     Bt=Med.biom[[s]][[h]][[e]]$Bt,
                                                     Bmsy=Med.biom[[s]][[h]][[e]]$Bmsy,
                                                     Ktch=Store.stuff[[s]]$Ktch,
                                                     CX.AX=1,
                                                     CX=1,
                                                     DAT=DAT,
                                                     Add.ktch="YES",
                                                     LOW="0%",
                                                     HIGH="100%")
            if(e==1&h==1)mtext(capitalize(names(SPM.preds)[s]),3,cex=1) 
            Crip=Efficien.scens[e]
            
          }
        }
        
      }
    }
  }
  mtext("Financial year",1,cex=1.2,line=0.75,outer=T)
  mtext("Relative biomass",2,cex=1.2,outer=T,las=3)
  mtext(side = 4, line = 0.75, 'Total catch (tonnes)',las=3,outer=T,
        col=rgb(0.1,0.1,0.8,alpha=0.6),cex=1.2)
  dev.off()
  rm(HR.o.scens)
  
  
  #Plot MSY  
  fn.fig(paste(hNdl,"/Outputs/Figure MSY_SPM",sep=""), 2400, 2400)
  HR.o.scens=Mx.init.harv[1]
  smart.par(n.plots=length(compact(SPM.preds)),MAR=c(1.2,2,1.5,1.25),
            OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
  par(las=1,cex.axis=.8)
  for(s in 1: N.sp)
  {
    if(!is.null(SPM.preds[[s]]))
    {
      HR.o.scens=Mx.init.harv[s]
      for(h in 1:length(HR.o.scens))
      {
        nne=length(Efficien.scens)
        for(e in 1:nne)
        {
          if(!is.null(SPM.preds[[s]][[h]][[e]]))
          {
            dummy=unlist(subListExtract(SPM.preds_uncertainty[[s]][[h]][[e]],"MSY"))
            plot(density(dummy,adjust = 2),main="",ylab="")
            if(e==1&h==1)mtext(capitalize(names(SPM.preds)[s]),3,cex=1) 
            Crip=Efficien.scens[e]
            legend("right",paste(round(median(dummy))," tonnes",sep=""),
                   bty='n',cex=1.1,title='Median MSY')
          }
        }
      }
      mtext("Catch (tonnes)",1,line=0.75,outer=T)
      mtext("Density",2,outer=T,las=3)
    }
  }
  dev.off()

  #Output parameter estimates
  Tab.par.estim.SPM=vector('list',length(SPM.preds))
  for(s in 1: N.sp)
  {
    if(!is.null(SPM.preds[[s]]))
    {
      for(h in 1:length(HR.o.scens))
      {
        for(e in 1:nne)
        {
          #Samp=Estim.par.samples[[s]][[h]][[e]]
          #MLE=apply(Samp,2,function(x ) round(median(exp(x)),2))
          #std=apply(Samp,2,function(x ) round(sd(exp(x)),3))
          #Nms=colnames(Samp)
          fit=Store.SPM[[s]][[h]][[e]]
          MLE=round(fit$par,2)
          std=round(sqrt(diag(solve(fit$hessian))),2)
          Nms=paste("log",names(MLE))
          
          Tab=as.data.frame(matrix(paste(MLE," (",std,")",sep=''),nrow=1))
          names(Tab)=Nms
          Tab=cbind(Species=capitalize(names(SPM.preds)[s]),Tab)
        }
      }
      Tab.par.estim.SPM[[s]]=Tab
    }
  }
  Tab.par.estim.SPM=do.call(rbind,Tab.par.estim.SPM)
  fn.word.table(WD=getwd(),TBL=Tab.par.estim.SPM,Doc.nm="Table 2. SPM estimates",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
  rm(HR.o.scens)
  
}



#---Catch-MSY RESULTS------
if(Do.Ktch.MSY=="YES")
{
  
  #r priors   
  fn.fig(paste(hNdl,"/Outputs/Figure 1_Prior_r",sep=""), 2000, 2000)
  smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(s in 1: N.sp)
  {
    NMs=capitalize(names(store.species)[s])
    if(NMs=="Low") NMs="Low resilience"
    if(NMs=="Very.low") NMs="Very low resilience"
    plot(density(rgamma(10000, shape = store.species[[s]]$r.prior$shape, rate = store.species[[s]]$r.prior$rate)),
         lwd=3,main=NMs,xlab="",ylab="",cex.lab=2,cex.axis=1.15,col=1,xlim=c(0,.4),yaxt='n')
  }
  mtext(expression(paste(plain("Intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),1,0.5,cex=1.35,outer=T)
  mtext("Density",2,0,las=3,cex=1.35,outer=T)
  dev.off()
  
  YrS=sort(unique(Tot.ktch$finyear))
  
  #Relative biomass
  CL="grey55"
  CL.mean="transparent"
  
  Low.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (0+Nper)/100))   #get percentiles
  High.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (100-Nper)/100))

  COLS=colfunc(3)
  fn.plot.percentile=function(DAT,YR,ADD.prob,add.RP.txt,CEX,CX.AX,Ktch,PERCENTIL,Add.ktch)
  {
    
    # data percentile
    #Nper=(100-60)/2 # %60%
    Nper=(100-PERCENTIL)/2 
    LOW=Low.percentile(Nper,DAT)
    UP=High.percentile(Nper,DAT)
    
    
    #construct polygons
    Year.Vec <-  fn.cons.po(YR,YR)
    Biom.Vec <-fn.cons.po(LOW,UP)
    
    #plot
    MED=apply(DAT, 1, function(x) quantile(x, .5))
    plot(YR,MED,ylim=c(0,1.05),ylab="",xlab="",xaxt='n',type='l',cex=0.3,lwd=1.5,
         pch=19,cex.axis=CX.AX,xaxs="i",yaxs="i")
    polygon(Year.Vec, Biom.Vec, col = rgb(.1,.1,.1,alpha=.1), border = "grey70")
    
    #reference points
    abline(h=B.target,lwd=1.15,col= Col.RP[3])
    abline(h=B.threshold,lwd=1.15,col=Col.RP[2])
    abline(h=B.limit,lwd=1.15,col= Col.RP[1])
    
    if(add.RP.txt=="YES")
    {
      text(YR[4],B.target,"Target",pos=3,cex=1.1)
      text(YR[4],B.threshold,"Threshold",pos=3,cex=1.1)
      text(YR[4],B.limit,"Limit",pos=3,cex=1.1)
    }
    
    #add probs
    if(ADD.prob=="YES")
    {
      store.cons.like=add.probs(id.yr=match(Current,YR),YR,DAT,UP,LOW,SRT=0,CEX)
    }
    axis(1,at=YR,labels=F,tck=-0.015)
    axis(1,at=seq(YR[1],YR[length(YR)],5),labels=seq(YR[1],YR[length(YR)],5),tck=-0.03,cex.axis=CX.AX)
    
    #add catch
    if(Add.ktch=="YES")
    {
      par(new=T)
      plot(YR,Ktch,type='l',col=Ktch.CL,xlab="",ylab="",
           axes=F,lwd=1.5)
      axis(side = 4)
    }
    
    if(ADD.prob=="YES") return(store.cons.like)
  }
  
  
  fn.fig(paste(hNdl,"/Outputs/Figure 3_Biomass_Catch_MSY",sep=""), 2400, 1800)
  smart.par(n.plots=length(compact(store.species)),MAR=c(1.2,2,1.5,1.75),
            OMA=c(2,1.75,.2,2.1),MGP=c(1,.5,0))
  par(las=1,cex.axis=1.1)
  for(s in 1: N.sp)
  {
    Yrs=Store.stuff[[s]]$yrs
    Ktch_MSY_Rel.bio=with(store.species[[s]]$KTCH.MSY$BaseCase,sweep(bt, 2, k, `/`))

    #Percentile   
    Store.cons.Like.SRM[[s]]=fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES",add.RP.txt="NO",
                                         CEX=1,CX.AX=1.1,Ktch=Store.stuff[[s]]$Ktch,
                                         PERCENTIL=100,Add.ktch="YES")
    NMs=capitalize(names(store.species)[s])
    mtext(NMs,3,0)
  }
  mtext("Financial year",1,cex=1.2,line=0.75,outer=T)
  mtext("Relative biomass",2,cex=1.2,outer=T,las=3)
  mtext(side = 4, line = 0.75, 'Total catch (tonnes)',las=3,outer=T,
        col=Ktch.CL,cex=1.2)
  dev.off()
  
  if(Modl.rn=='first')  
  {
    fn.fig("Figure 3_Biomass_Catch_MSY_WorstCase", 2400, 1800)
    smart.par(n.plots=length(compact(store.species)),MAR=c(1.2,2,1.5,1.75),
              OMA=c(2,1.75,.2,2.1),MGP=c(1,.5,0))
    par(las=1,cex.axis=1.1)
    for(s in 1: N.sp)
    {
      Yrs=Store.stuff[[s]]$yrs
      Ktch_MSY_Rel.bio=with(store.species[[s]]$KTCH.MSY$WorstCase,sweep(bt, 2, k, `/`))
      
      #Percentile   
      Store.cons.Like.SRM[[s]]=fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES",add.RP.txt="NO",
                                                  CEX=1,CX.AX=1.1,Ktch=Store.stuff[[s]]$Ktch,
                                                  PERCENTIL=100,Add.ktch="YES")
      NMs=capitalize(names(store.species)[s])
      mtext(NMs,3,0)
    }
    mtext("Financial year",1,cex=1.2,line=0.75,outer=T)
    mtext("Relative biomass",2,cex=1.2,outer=T,las=3)
    mtext(side = 4, line = 0.75, 'Total catch (tonnes)',las=3,outer=T,
          col=rgb(0.1,0.1,0.8,alpha=0.6),cex=1.2)
    dev.off()
  }
  
  fn.fig(paste(hNdl,"/Outputs/Figure MSY_Catch.MSY",sep=""), 2400, 2400)
  smart.par(n.plots=length(compact(store.species)),MAR=c(1.2,2,1.5,1.75),
            OMA=c(2,3,.2,2.1),MGP=c(1,.5,0))
  par(las=1,cex.axis=1.1)
  for(s in 1: N.sp)
  {
    dummy=store.species[[s]]$KTCH.MSY$BaseCase$msy
    plot(density(dummy,adjust = 2),main="",ylab="")
    mtext(capitalize(names(store.species)[s]),3,cex=.95)
    legend("right",paste(round(median(dummy))," tonnes",sep=""),
           bty='n',cex=1.1,title='Median MSY')
  }
  mtext("Catch (tonnes)",1,line=0.75,outer=T)
  mtext("Density",2,line=1,outer=T,las=3)
  dev.off()
}


#---RISK ------
Like.ranges=list(L1=c(0,0.0499999),L2=c(0.05,0.2),L3=c(0.20001,0.5),L4=c(0.50001,1))
Risk.tab=data.frame(Consequence=paste("C",1:4,sep=""),
                    L1=NA,L2=NA,L3=NA,L4=NA,
                    Max.Risk.Score=NA)

fn.risk=function(likelihood)
{
  consequence=names(likelihood)
  TAB=Risk.tab
  for(n in 1:length(likelihood))
  {
    id=match(names(likelihood)[n],TAB$Consequence)
    idd=which(unlist(lapply(Like.ranges,function(x) check.in.range(likelihood[n],x,fatal=F))))
    TAB[id,idd+1]="X"
    TAB$Max.Risk.Score[id]=id*idd
  }
  return(TAB)
}
Risk.SPM=Store.cons.Like.SPM
Risk.SRM=Store.cons.Like.SRM
for(s in 1: N.sp)
{
  #SPM
  if(!is.null(Store.cons.Like.SPM[[s]])) Risk.SPM[[s]]=fn.risk(likelihood=unlist(Store.cons.Like.SPM[[s]]))
  #SRM
  Risk.SRM[[s]]=fn.risk(likelihood=unlist(Store.cons.Like.SRM[[s]]))
}

fn.AD=function(add.text)
{
  X.rng=1:N.sp
  x.Vec <-  fn.cons.po(0:(N.sp+1),0:(N.sp+1))
  nn=N.sp+2
  negigible.Vec <- fn.cons.po(rep(0,nn),rep(2,nn))
  low.Vec <- fn.cons.po(rep(2,nn),rep(4,nn))
  medium.Vec <- fn.cons.po(rep(4,nn),rep(8,nn))
  high.Vec <- fn.cons.po(rep(8,nn),rep(12,nn))
  severe.Vec <- fn.cons.po(rep(12,nn),rep(16,nn))
  plot(X.rng,ylim=c(0,16),col="transparent",ylab="",yaxs="i",xlab="",xaxt='n',yaxt='n')
  
  polygon(x.Vec, negigible.Vec, col = 'cornflowerblue', border = "transparent")
  polygon(x.Vec, low.Vec, col = 'olivedrab3', border = "transparent")
  polygon(x.Vec, medium.Vec, col = 'yellow', border = "transparent")
  polygon(x.Vec, high.Vec, col = 'orange', border = "transparent")
  polygon(x.Vec, severe.Vec, col = 'red', border = "transparent")
  
  if(add.text=="YES")
  {
    text(.55,mean(negigible.Vec),"Negligible",pos=4,cex=.75,col=rgb(.1,.1,.1,.3))
    text(.55,mean(low.Vec),"Low",pos=4,cex=.75,col=rgb(.1,.1,.1,.3))
    text(.55,mean(medium.Vec),"Medium",pos=4,cex=.75,col=rgb(.1,.1,.1,.3))
    text(.55,mean(high.Vec),"High",pos=4,cex=.75,col=rgb(.1,.1,.1,.3))
    text(.55,mean(severe.Vec),"Severe",pos=4,cex=.75,col=rgb(.1,.1,.1,.3))  
  }
  box()
}

fn.fig("Figure 5_Risk", 2400, 1400)
par(mfrow=c(2,1),mar=c(.5,.5,.5,.95),oma=c(6,2,.2,.1),las=1,cex.axis=.8,mgp=c(1,.3,0))

#SPM
fn.AD(add.text="YES")
for(s in 1:N.sp)
{
  if(!is.null(Risk.SPM[[s]])) 
  {
    Mx=max(Risk.SPM[[s]]$Max.Risk.Score)[1]
    points(s,Mx,type='h',lwd = 15,lend=1)
  }
}
axis(2,seq(0,16,1),labels=F,tck=-0.015)
axis(2,seq(0,16,2),seq(0,16,2),tck=-0.03,cex.axis=.7)
legend('topright',"SPM",bty='n')

#SRM
fn.AD(add.text="NO")
for(s in 1:N.sp)
{
  Mx=max(Risk.SRM[[s]]$Max.Risk.Score)[1]
  points(s,Mx,type='h',lwd = 15,lend=1)
  
}
axis(2,seq(0,16,1),labels=F,tck=-0.015)
axis(2,seq(0,16,2),seq(0,16,2),tck=-0.03,cex.axis=.7)
axis(1,1:N.sp,labels=names(Risk.SRM),tck=-0.01,las=2,cex.axis=.7)
legend('topright',"SRM",bty='n')
mtext("Highest risk score",2,outer=T,las=3,cex=1.25,line=0.85)
dev.off()


#---National Scalloped HH Assessment -----------------------------------------------------------------------
if(Asses.Scalloped.HH)
{
  SP="scalloped hammerhead"
  B.target=.6
  B.threshold=.4
  B.limit=.2
  
  #Scenarios
  Scen.HH=data.frame(Scenario=c("Base case", "Double Indonesian catch","Double discards and rec"),
                     stringsAsFactors=F)
  
  #select scalloped HH index
  s=match(SP,names(cpue.list))
  
  #Define discarding fisheries
  Calculate.discarding=c('PFT','C019','C066','C070','CSLP','EGBS','EGP','KP','KTR','NBP',
                         'OANCGCWC','OP','SBP','SBS','SBSC','SCT','SWT')
  
  #SPM
  if(Do.SPM=="YES")
  {
    Scens=Scen.HH$Scenario
    
    #Estimate parameters
    Store.SPM.scallopedHH=vector('list',length(Scens))
    names(Store.SPM.scallopedHH)=Scens
    Store.stuff.scallopedHH=Store.SPM.scallopedHH
    
    SPM.preds.scallopedHH=vector('list',length(Store.SPM.scallopedHH))
    names(SPM.preds.scallopedHH)=names(Store.SPM.scallopedHH)
    
    SPM.preds_uncertainty.scallopedHH=vector('list',length(Store.SPM.scallopedHH))
    names(SPM.preds_uncertainty.scallopedHH)=names(Store.SPM.scallopedHH)
    Estim.par.samples.scallopedHH=SPM.preds_uncertainty.scallopedHH
    
    for(sc in 1:length(Scens))
    {
      Id=match(Specs$SP.group[s],names(cpue.list))
      
      #catch
      ct=Tot.ktch%>%filter(SP.group==Specs$SP.group[s])
      if(Scens[sc]=="Double Indonesian catch")
      {
        ct=ct%>%mutate(LIVEWT.c=ifelse(FishCubeCode=="Indo",2*LIVEWT.c,LIVEWT.c))
      }
      if(Scens[sc]=="Double discards and rec")
      {
        ct=ct%>%mutate(LIVEWT.c=ifelse(FishCubeCode=="Recreational",
                                       2*LIVEWT.c,LIVEWT.c))
        ct=ct%>%mutate(LIVEWT.c=ifelse(finyear>2006 & FishCubeCode%in%Calculate.discarding,
                                       2*LIVEWT.c,LIVEWT.c))
      }
      
      ct=ct%>%
        group_by(finyear)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
        data.frame
      
      all.iers=seq(min(ct$finyear),max(ct$finyear))
      
      if(length(which(!all.iers%in%ct$finyear))>0)
      {
        aa=all.iers[which(!all.iers%in%ct$finyear)]
        aa1=ct[1:length(aa),]
        aa1[,]=0
        aa1$finyear=aa
        ct=rbind(ct,aa1)%>%arrange(finyear)
      }
      
      #r prior
      Idd=match(Specs$SP.group[s],names(store.species))
      r.prior=store.species[[Idd]]$r.prior.normal$mean
      r.prior.sd=store.species[[Idd]]$r.prior.normal$sd
      
      Store.CPUE.eff.dummy=NULL
      QS.dummy=NULL
      CPUE.yr.dummy=NULL
      n.cpue.dummy=NULL
      CPUE.CV=NULL
      if(!is.null(cpue.list[[Id]]))
      {
        #Get cpues
        CPUE.1=cpue.list[[Id]]$Nat
        if(!is.null(CPUE.1))CPUE.1=CPUE.1%>%
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
            
            too.large=which(CPUE.CV[[ci]]>.6) #drop observation with too large CVs
            CPUE[[ci]][too.large]=NA
            CPUE.CV[[ci]][too.large]=NA
            
            CPUE.yr[[ci]]=x%>%filter(!is.na(MeAn))%>%
              dplyr::select(yr)
          }
        }
        
        #Initial values for estimable pars
        Mx.ktch=max(ct$LIVEWT.c,na.rm=T)
        K.init=k.times.mx.ktch*Mx.ktch
        r.init=Init.r[[match(Specs$SP.group[s],names(Init.r))]]
        QS=Estimable.qs[[match(Specs$SP.group[s],names(Estimable.qs))]]
        
        #loop over scenarios
        HR.o.scens=Mx.init.harv[match(Specs$SP.group[s],names(Mx.init.harv))]
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
            fn_ob=function(theta)SPM(Init.propK=B.init,
                                     cpue=CPUE.eff.scen,
                                     cpue.CV=CPUE.CV,
                                     Qs=QS,
                                     Ktch=ct$LIVEWT.c,
                                     theta,
                                     HR_init=HR_o,
                                     HR_init.sd=HR_o.sd,
                                     r.mean=r.prior,
                                     r.sd=r.prior.sd,
                                     usePenalties=usePen)$negLL

            #fit model
            if(estim.q=="YES")theta= c(k=log(K.init),log(QS[which(!is.na(QS))]))
            if(estim.q=="NO")theta= c(k=log(K.init))
            Lw.bound=log(c(Low.bound.K*Mx.ktch,rep(1e-6,length(theta)-1)))
            Up.bound=log(c(Up.bound.K*Mx.ktch,rep(1,length(theta)-1)))
            
            if(!fix.r=="YES")
            {
              theta=c(theta,r=log(r.init))
              Lw.bound=c(Lw.bound,log(0.01))
              Up.bound=c(Up.bound,log(0.75))
            }
              


            if(what.like=='kernel')
            {
              theta=c(theta,sigma=log(0.2))
              Lw.bound=c(Lw.bound,log(1e-2))
              Up.bound=c(Up.bound,log(1))
            }
            
            if(minimizer=='optim')
            {
              OptiM=optim(theta,fn_ob,method="L-BFGS-B",lower =Lw.bound,upper = Up.bound,hessian=T)
              if(Remove.bounds)
              {
                paramscale = magnitude(theta)
                OptiM=optim(jitter(OptiM$par),fn_ob,method="L-BFGS-B",hessian=T,
                            control = list(maxit = 1000, parscale = paramscale))
              }
              
              dummy.eff[[e]]=OptiM
            }
            if(minimizer=='nlminb')
            {
              nlmb <- nlminb(theta, fn_ob, gradient = NULL,lower =Lw.bound,upper = Up.bound)
              if(Remove.bounds) nlmb <- nlminb(jitter(nlmb$par), fn_ob, gradient = NULL)
              
              nlmb$hessian=hessian(fn_ob,nlmb$par)
              
              dummy.eff[[e]]=nlmb
            }
            Store.CPUE.eff[[e]]=CPUE.eff.scen
          }
          dummy[[h]]=dummy.eff
        }
        Store.SPM.scallopedHH[[sc]]=dummy
        rm(HR.o.scens)
        Store.CPUE.eff.dummy=Store.CPUE.eff
        QS.dummy=QS
        CPUE.yr.dummy=CPUE.yr
        n.cpue.dummy=length(CPUE)
      }
      Store.stuff.scallopedHH[[sc]]=list(cpue=Store.CPUE.eff.dummy,CPUE.CV=CPUE.CV,Qs=QS.dummy,Ktch=ct$LIVEWT.c,
                                         r.mean=r.prior,r.sd=r.prior.sd,yrs=all.iers,
                                         cpue.yrs=CPUE.yr.dummy,n.cpues=n.cpue.dummy)
      
      #Evaluate model at MLE
      HR.o.scens=Mx.init.harv[s]
      dumy.pred=vector('list',length(HR.o.scens))
      names(dumy.pred)=paste("HR=",HR.o.scens)
      for(h in 1:length(HR.o.scens))
      {
        dummy.eff=vector('list',length(Efficien.scens))
        names(dummy.eff)=Efficien.scens
        for(e in 1:length(Efficien.scens))
        {
          dummy.eff[[e]]=SPM(Init.propK=B.init,
                             cpue=Store.stuff.scallopedHH[[sc]]$cpue[[e]],
                             cpue.CV=Store.stuff.scallopedHH[[sc]]$CPUE.CV,
                             Qs=Store.stuff.scallopedHH[[sc]]$Qs,
                             Ktch=Store.stuff.scallopedHH[[sc]]$Ktch,
                             theta=Store.SPM.scallopedHH[[sc]][[h]][[e]]$par,
                             HR_init=HR.o.scens[h],
                             HR_init.sd=HR_o.sd,
                             r.mean=Store.stuff.scallopedHH[[sc]]$r.mean,
                             r.sd=Store.stuff.scallopedHH[[sc]]$r.sd,
                             usePenalties=usePen)
        }
        dumy.pred[[h]]=dummy.eff
      }
      SPM.preds.scallopedHH[[sc]]=dumy.pred
      
      
      #Get uncertainty  
      #Draw random samples of estimable pars
      dummy=vector('list',length(HR.o.scens))
      names(dummy)=HR.o.scens
      for(h in 1:length(HR.o.scens))
      {
        dummy1=vector('list',length(Efficien.scens))
        names(dummy1)=Efficien.scens
        for(e in 1:length(Efficien.scens))
        {
          dummy1[[e]]=fn.un(fit=Store.SPM.scallopedHH[[sc]][[h]][[e]],n=N.monte)
        }
        dummy[[h]]=dummy1
      }
      Estim.par.samples.scallopedHH[[sc]]=dummy
      
      #Use par sample to predict quantities of interest
      dumy.pred=vector('list',length(HR.o.scens))
      names(dumy.pred)=HR.o.scens
      for(h in 1:length(HR.o.scens))
      {
        dummy.eff=vector('list',length(Efficien.scens))
        names(dummy.eff)=Efficien.scens
        for(e in 1:length(Efficien.scens))
        {
          dum=vector('list',N.monte)
          for(x in 1:N.monte)
          {
            dum[[x]]=SPM(Init.propK=B.init,
                         cpue=Store.stuff.scallopedHH[[sc]]$cpue[[e]],
                         cpue.CV=Store.stuff.scallopedHH[[sc]]$CPUE.CV,
                         Qs=Store.stuff.scallopedHH[[sc]]$Qs,
                         Ktch=Store.stuff.scallopedHH[[sc]]$Ktch,
                         theta=Estim.par.samples.scallopedHH[[sc]][[h]][[e]][x,],
                         HR_init=HR.o.scens[h],
                         HR_init.sd=HR_o.sd,
                         r.mean=Store.stuff.scallopedHH[[sc]]$r.mean,
                         r.sd=Store.stuff.scallopedHH[[sc]]$r.sd,
                         usePenalties=usePen)
          }
          dummy.eff[[e]]=dum
        }
        dumy.pred[[h]]=dummy.eff
      }
      SPM.preds_uncertainty.scallopedHH[[sc]]=dumy.pred
      rm(HR.o.scens)
    }
  }
  
  #Catch-MSY
  if(Do.Ktch.MSY=="YES")
  {
    #Tested scenarios
    ktch_msy_scen=vector('list',length(Scen.HH$Scenario))
    names(ktch_msy_scen)=Scen.HH$Scenario
    for(sc in 1:length(ktch_msy_scen))
    {
      USR="Yes"
      Scen.start.bio=SCENARIOS[[1]]$Initial.dep
      FinBio=FINALBIO
      ktch_msy_scen[[sc]]=list(r.prior=SCENARIOS[[1]]$R.prior,
                               user=USR,
                               k.lower=Low.bound.K,k.upper=Up.bound.K, 
                               startbio=Scen.start.bio,finalbio=FinBio,res=RESILIENCE[[s]],
                               niter=SIMS,sigR=SCENARIOS[[1]]$Error)
      if(!is.na(ktch_msy_scen[[sc]]$r.prior[1])) ktch_msy_scen[[sc]]$r.prior=unlist(store.species[[s]]$r.prior)
    }
    
    #Execute Catch-MSY  
    PATH=paste(PHat,Specs$SP.group[s],"National assessment",sep='/')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
    setwd(file.path(PATH))
    Path.ktch_msy=getwd()
    Ktch_MSY.scallopedHH=ktch_msy_scen
    for(sc in 1:length(ktch_msy_scen))
    {
      print(paste("-------------",Scen.HH$Scenario[sc],"scenario"))
      Yrs=Store.stuff.scallopedHH[[sc]]$yrs   
      Tot.Ktch=Store.stuff.scallopedHH[[sc]]$Ktch  
      yr.future=Current+(1:years.futures)
      ct.future=rep(mean(Tot.Ktch[(length(Tot.Ktch)-4):length(Tot.Ktch)]),years.futures)
      
      
      Folder=names(ktch_msy_scen)[sc]
      if(!file.exists(paste(PATH,Folder,sep="/"))) dir.create(paste(PATH,Folder,sep="/"))
      setwd(paste(PATH,Folder,sep="/"))
      
      Ktch_MSY.scallopedHH[[sc]]=Catch_MSY(ct=Tot.Ktch,
                                           yr=Yrs,
                                           r.prior=ktch_msy_scen[[sc]]$r.prior,
                                           user=ktch_msy_scen[[sc]]$user,
                                           k.lower=ktch_msy_scen[[sc]]$k.lower,
                                           k.upper=ktch_msy_scen[[sc]]$k.upper,
                                           startbio=ktch_msy_scen[[sc]]$startbio,
                                           finalbio=ktch_msy_scen[[sc]]$finalbio,
                                           res=ktch_msy_scen[[sc]]$res,
                                           n=ktch_msy_scen[[sc]]$niter,
                                           sigR=ktch_msy_scen[[sc]]$sigR,
                                           ct.future=ct.future,           
                                           yr.future=yr.future)
      
      #Export outputs
      Table1_ktch_MSY=with(Ktch_MSY.scallopedHH[[sc]],data.frame(`geom. mean r`,`r +/- 1.96 SD`,`geom. mean k (tons)`,`k +/- 1.96 SD (tons)`,
                                                                 `geom. mean MSY (tons)`,`MSY +/- 1.96 SD (tons)`))
      write.csv(Table1_ktch_MSY,"Table1_ktch_MSY.csv",row.names=F)
      
      #output base case catch to run 90% decline scenario
      if(Folder=="Base case")
      {
        write.csv(data.frame(finyear=Yrs,Total.ktch=Tot.Ktch),
                  'C:/Matias/External collaborations/Australian Scalloped HH/Roff_90_percent_depletion/ct.csv',
                  row.names=F) 
      }
    }
  }
  
  #aSPM 
  if(Do.aSPM=="YES")
  {
    l=match("scalloped hammerhead",names(cpue.list))
    
    #catch
    fish=Tot.ktch%>%
      filter(SP.group==names(cpue.list)[l])%>%
      group_by(finyear)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
      data.frame%>%
      rename(year=finyear,
             catch=LIVEWT.c)
    
    #cpue
    Id=match(names(cpue.list)[l],names(cpue.list))
    CPUE=compact(cpue.list[[Id]])
    len.cpue=length(CPUE)
    if(len.cpue>1)
    {
      if("Nat"%in%names(CPUE))
      {
        a=CPUE[-match("Nat",names(CPUE))]
        
        CPUE2=CPUE[[match("Nat",names(CPUE))]]%>%    
          dplyr::select(yr,MeAn,CV)%>%
          rename(year=yr,cpue2=MeAn,se2=CV)
        
        iid=which(CPUE2$se2>MAX.CV)
        CPUE2$cpue2[iid]=NA
        CPUE2$se2[iid]=NA
        
        q1.yrs=as.numeric(substr(a$TDGDLF.mon$Finyear,1,4))
        q2.yrs=as.numeric(substr(a$TDGDLF.day$Finyear,1,4))
        CPUE=do.call(rbind,a)%>%
          mutate(yr=as.numeric(substr(Finyear,1,4)),
                 MeAn=Mean)
      }else
      {
        a=CPUE
        q1.yrs=as.numeric(substr(a$TDGDLF.mon$Finyear,1,4))
        q2.yrs=as.numeric(substr(a$TDGDLF.day$Finyear,1,4))
        CPUE=do.call(rbind,a)%>%
          mutate(yr=as.numeric(substr(Finyear,1,4)),
                 MeAn=Mean)
        iid=which(CPUE$CV>MAX.CV)
        CPUE$MeAn[iid]=NA
        CPUE$CV[iid]=NA
      }
      
    }
    if(len.cpue==1)
    {
      CPUE=CPUE[[1]]
      colnames(CPUE)[2]="MeAn"
      if(!'yr'%in%names(CPUE))
      {
        CPUE=CPUE%>%
          mutate(yr=as.numeric(substr(Finyear,1,4)))
      }
      iid=which(CPUE$CV>MAX.CV)
      CPUE$MeAn[iid]=NA
      CPUE$CV[iid]=NA
    }
    
    CPUE=CPUE%>%    
      dplyr::select(yr,MeAn,CV)%>%
      rename(year=yr,cpue=MeAn,se=CV)
    if(min(CPUE$se,na.rm=T)>1)  CPUE$se=CPUE$se/100
    
    fish=fish%>%left_join(CPUE,by='year')
    fish0.ktch=fish[1:20,]%>%
      mutate(year=year-20,
             catch=0,
             cpue=NA,
             se=NA)
    fish=rbind(fish0.ktch,fish)
    iid=which(fish$se>MAX.CV)
    fish$cpue[iid]=NA
    fish$se[iid]=NA
    
    if(exists(c('q1.yrs','q2.yrs')))
    {
      ia=fish$year[!is.na(fish$cpue)]
      q1.yrs=q1.yrs[which(q1.yrs%in%ia)]
      q2.yrs=q2.yrs[which(q2.yrs%in%ia)]
    }
    
    if(exists('CPUE2')) fish=fish%>%left_join(CPUE2,by='year')
    if(!exists('CPUE2')) fish=fish%>%mutate(cpue2=NA,se2=NA)
    
    
    
    #init par values
    Id=match(Specs$SPECIES[match(names(aSPM.init)[l],Specs$Name)],names(Store.age.comp))
    
    #global pars
    glb=Store.age.comp[[Id]]$glb
    glb$nyrs=nrow(fish)
    
    #relations at age
    props=Store.age.comp[[Id]]$props
    
    #fit model
    pars=aSPM.init[[l]]
    ans <- fitASPM(initpar=pars,infish=fish,inglb=glb,inprops=props)
    fishery <- dynamics(ans$par,infish=fish,inglb=glb,inprops = props)
    
    #get uncertainty
    Par.samp=fn.un(fit=ans,n=N.monte)
    fishery.MC=vector('list',nrow(Par.samp))
    for(f in 1:nrow(Par.samp))
    {
      fishery.MC[[f]]=dynamics(Par.samp[f,],infish=fish,inglb=glb,inprops = props)
    }
    if(exists(c('q1.yrs','q2.yrs'))) rm(q1.yrs,q2.yrs)
    if(exists('CPUE2')) rm(CPUE2)
    Store.aSPM.scallopedHH=list(fit=ans,quantities=fishery,fish=fish,quantities.MC=fishery.MC)
  }
  
  #Display catches
  hNdl.HH=paste(hNdl,'Outputs/each_species/scalloped hammerhead/National assessment/',sep='/')
  
  fn.fig(paste(hNdl.HH,'Figure 1_catches',sep=''),2400,2400) 
  par(mar=c(3.5,3.25,1.5,1),las=1,mgp=c(1,.75,0),cex.axis=1.25)
  
  Tot.col.bg=rgb(.1,.5,1,alpha=.75)
  Tot.col=rgb(.1,.5,1,alpha=.3)
  for(s in s)
  {
    ddd=subset(Tot.ktch,SP.group==Specs$SP.group[s])%>%
      group_by(Type,finyear)%>%
      summarise(Tot=sum(LIVEWT.c,na.rm=T))%>%
      filter(!is.na(Tot))
    
    ddd.all=ddd%>%
            group_by(finyear)%>%
            summarise(Tot=sum(Tot,na.rm=T))
    plot(all.yrs,ddd.all$Tot,col=Tot.col,ylab="",xlab="",type='o',cex=1.25,
         lwd=2,pch=21,bg=Tot.col.bg,ylim=c(0,max(ddd.all$Tot)),xaxt='n')
    unik.T=unique(ddd$Type)
    for(u in 1:length(unik.T))
    {
      cl=COLs.type[match(unik.T[u],names(COLs.type))]
      with(subset(ddd,Type==unik.T[u]),points(finyear,Tot,type='o',cex=1.75,lwd=2,pch=21,bg=cl))
    }
    legend('topleft',c(names(COLs.type)[1:4],"Total"),pt.bg=c(COLs.type[1:4],Tot.col.bg),
           bty='n',pch=21,cex=1.15,pt.cex=c(rep(2,4),1.5),lwd=2,col=c(rep('black',4),Tot.col))
    axis(1,all.yrs,F,tck=-.01)
    Yrs.lab=round(seq(all.yrs[1],all.yrs[length(all.yrs)],10)/10)*10
    axis(1,Yrs.lab,Yrs.lab,tck=-.02)
    
  }
  mtext("Financial year",1,line=-1.2,cex=1.5,outer=T)
  mtext("Total catch (tonnes)",2,las=3,line=-1.2,cex=1.5,outer=T)
  dev.off()
  
  #Display vulnerability at age
  fn.fig(paste(hNdl.HH,'Figure 1_selectivity_gillnet',sep=''),2400,2400) 
  par(mar=c(3.5,3.5,1.5,1),mgp=c(2,.75,0),cex.axis=1.25,cex.lab=1.5)
  fn.plt.age.comp(sim=Store.age.comp$`19001`$GN$dat,
                  obs=Store.age.comp$`19001`$GN$pred.age,
                  sel=Store.age.comp$`19001`$GN$Selectivity,
                  Title="Gillnet")
  dev.off()
  
  fn.fig(paste(hNdl.HH,'Figure 1_selectivity_longline',sep=''),2400,2400) 
  par(mar=c(3.5,3.5,1.5,1),mgp=c(2,.75,0),cex.axis=1.25,cex.lab=1.5)
  fn.plt.age.comp(sim=Store.age.comp$`19001`$LL$dat,
                  obs=Store.age.comp$`19001`$LL$pred.age,
                  sel=Store.age.comp$`19001`$LL$Selectivity,
                  Title="Longline")
  dev.off()

  
  #Display Steepness  
  fn.fig(paste(hNdl.HH,'Figure 1_steepness',sep=''),2400,2400) 
  dist=with(store.species.steepness$`scalloped hammerhead`,density(rgamma(10000,shape,rate),adjust = 2))
  plot(dist,xlim=c(0,1),col='red',xlab="Steepness",yaxt='n',ylab="Density",lwd=2,
       main=paste("mean=",round(mean(mean(dist$x)),2)))
  dev.off()
  
  
  #Display SPM
  if(Do.SPM=="YES")
  {
    #Plot obs VS pred cpues  
    Paz=paste(hNdl.HH,"SPM.fit/",sep="")
    if(!file.exists(file.path(Paz))) dir.create(file.path(Paz))
    fn.fig(paste(Paz,"Fit",sep=""),1600,2400)
    smart.par(n.plots=length(compact(SPM.preds.scallopedHH)),MAR=c(1.2,2,2,1.25),
              OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
    for(sc in 1:length(Scens))
    {
      if(!is.null(SPM.preds.scallopedHH[[sc]]))
      {
        HR.o.scens=Mx.init.harv[s]
        #nrw=length(HR.o.scens)*length(Efficien.scens)
        #ncl=SPM.preds.scallopedHH[[sc]][[1]][[1]]$ln.cpue
        #ncl=length(ncl[!sapply(ncl,is.null)])
        #par(mfrow=c(nrw,ncl),mar=c(1.2,2,.2,.1),oma=c(1.5,1.75,1.5,1),las=1,cex.axis=.8,mgp=c(1,.42,0))
        for(h in 1:length(HR.o.scens))
        {
          for(e in 1:length(Efficien.scens))
          {
            for(x in 1:Store.stuff.scallopedHH[[sc]]$n.cpues)
              if(!is.null(SPM.preds.scallopedHH[[sc]][[h]][[e]]$ln.cpue[[x]]))
              {
                fn.plt.cpue(ob=SPM.preds.scallopedHH[[sc]][[h]][[e]]$ln.cpue[[x]],
                            ob.CV=Store.stuff.scallopedHH[[sc]]$CPUE.CV[[x]],
                            pred=SPM.preds.scallopedHH[[sc]][[h]][[e]]$ln.cpue.hat.full[[x]],
                            Convergence=Store.SPM.scallopedHH[[sc]][[h]][[e]]$convergence,
                            pred.fit=SPM.preds.scallopedHH[[sc]][[h]][[e]]$ln.cpue.hat[[x]])
                Crip=Efficien.scens[e]
                if(x==1) Crip=0
                par(font=2)
                legend("bottomleft",Scens[sc],bty='n',cex=1.15)
              }
          }
        }
         mtext("Financial year",1,line=.6,outer=T)
        mtext("lncpue",2,outer=T,las=3)
        rm(HR.o.scens)
        
      }
    }
    legend("bottom",c("observed","predicted"),pch=19,cex=1.25,col=c("orange","black"),bty='n')
    dev.off()

    
    #Plot biomass  
    Col.RP=c("red","orange","forestgreen")
    
    #1. Get median and percentiles
    Med.biom.scallopedHH=vector('list',length(Store.SPM.scallopedHH))
    names(Med.biom.scallopedHH)=names(Store.SPM.scallopedHH)
    for(sc in 1:length(Scens))
    {
      if(!is.null(Store.SPM.scallopedHH[[sc]]))
      {
        HR.o.scens=Mx.init.harv[s]
        dummy1=vector('list',length(HR.o.scens))
        names(dummy1)=HR.o.scens
        for(h in 1:length(HR.o.scens))
        {
          dummy2=vector('list',length(Efficien.scens))
          names(dummy2)=Efficien.scens
          for(e in 1:length(Efficien.scens))
          {
            dummy=subListExtract(SPM.preds_uncertainty.scallopedHH[[sc]][[h]][[e]],"Bt")
            dummy=sweep(do.call(rbind,dummy), 1, exp(Estim.par.samples.scallopedHH[[sc]][[h]][[e]][,1]), `/`)
            idd=which(round(dummy[,ncol(dummy)],1)<0.05) #remove cases where final biomass=0 
            if(length(idd)) dummy=dummy[-idd,]
            Bt.all=dummy
            Bt=apply(dummy,2,function(x) quantile(x,probs=c(.025,0.5,.975)))   
            Bt=Bt[,-ncol(Bt)] #remove future Bt
            dummy=subListExtract(SPM.preds_uncertainty.scallopedHH[[sc]][[h]][[e]],"Bmsy")    
            dummy=do.call(rbind,dummy)
            if(length(idd)) dummy=as.matrix(dummy[-idd,])
            Bmsy=apply(dummy,2,function(x) quantile(x,probs=c(.025,.5,.975)))
            colnames(Bmsy)="Bmsy"
            
            dummy2[[e]]=list(Bt=Bt,Bmsy=Bmsy,Bt.all=Bt.all)  
          }
          dummy1[[h]]=dummy2
        }
        Med.biom.scallopedHH[[sc]]=dummy1
        rm(HR.o.scens)
      }
    }
    
    #2. Plot
    do.col="NO"
    if(do.col=="NO") colfunc <- colorRampPalette(c("grey95","grey60"))
    if(do.col=="YES") colfunc <- colorRampPalette(c("aliceblue","lightblue3"))
    if(do.col=="NO") CL.ref.pt=c("black","grey35","grey55","grey80")
    if(do.col=="YES") CL.ref.pt=c('forestgreen','yellow','orange','red')
    
    
    fn.fig(paste(hNdl.HH,"Figure 2_Biomass_SPM",sep=""),1600,2400)
    HR.o.scens=Mx.init.harv[1]
    smart.par(n.plots=length(compact(SPM.preds.scallopedHH)),MAR=c(1.2,2,2,1.25),
              OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
    par(las=1,cex.axis=1)
    for(sc in 1: length(Scens))
    {
      if(!is.null(SPM.preds.scallopedHH[[sc]]))
      {
        HR.o.scens=Mx.init.harv[s]
        for(h in 1:length(HR.o.scens))
        {
          nne=length(Efficien.scens)
          for(e in 1:nne)
          {
            if(!is.null(SPM.preds.scallopedHH[[sc]][[h]][[e]]))
            {
              DAT=t(Med.biom.scallopedHH[[sc]][[h]][[e]]$Bt.all)
              fn.plt.bio.ktch(Yr=Store.stuff.scallopedHH[[sc]]$yrs,
                              Bt=Med.biom.scallopedHH[[sc]][[h]][[e]]$Bt,
                              Bmsy=Med.biom.scallopedHH[[sc]][[h]][[e]]$Bmsy,
                              Ktch=Store.stuff.scallopedHH[[sc]]$Ktch,
                              CX.AX=1,
                              CX=1,
                              DAT=DAT,
                              Add.ktch="NO",
                              LOW="2.5%",
                              HIGH="97.5%")
              if(e==1&h==1)mtext(Scens[sc],3,bty='n',cex=1.25) 
              Crip=Efficien.scens[e]
              
            }
          }
          
        }
      }
    }
    mtext("Financial year",1,cex=1.2,line=0.75,outer=T)
    mtext("Relative biomass",2,cex=1.2,outer=T,las=3)
    #mtext(side = 4, line = 0.75, 'Total catch (tonnes)',las=3,outer=T,
    #      col=rgb(0.1,0.1,0.8,alpha=0.6),cex=1.2)
    dev.off()
    rm(HR.o.scens)
    
    
    #Plot MSY  
    fn.fig(paste(hNdl.HH,"Figure MSY_SPM",sep=""), 1600, 2400)
    HR.o.scens=Mx.init.harv[1]
    smart.par(n.plots=length(compact(SPM.preds.scallopedHH)),MAR=c(1.2,2,1.5,1.25),
              OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
    par(las=1,cex.axis=1)
    for(sc in 1: length(Scens))
    {
      if(!is.null(SPM.preds.scallopedHH[[sc]]))
      {
        HR.o.scens=Mx.init.harv[s]
        for(h in 1:length(HR.o.scens))
        {
          nne=length(Efficien.scens)
          for(e in 1:nne)
          {
            if(!is.null(SPM.preds.scallopedHH[[sc]][[h]][[e]]))
            {
              dummy=unlist(subListExtract(SPM.preds_uncertainty.scallopedHH[[sc]][[h]][[e]],"MSY"))
              plot(density(dummy,adjust = 2),main="",ylab="")
              if(e==1&h==1)mtext(Scens[sc],3,bty='n',cex=1) 
              Crip=Efficien.scens[e]
              legend("right",paste(round(median(dummy))," tonnes",sep=""),
                     bty='n',cex=1.2,title='Median MSY')
            }
          }
        }
        mtext("Catch (tonnes)",1,line=0.75,outer=T)
        mtext("Density",2,line=0.5,outer=T,las=3)
      }
    }
    dev.off()
    
    
    #Output parameter estimates
    Tab.par.estim.SPM.scallopedHH=vector('list',length(SPM.preds.scallopedHH))
    for(sc in 1: length(Scens))
    {
      for(h in 1:length(HR.o.scens))
      {
        for(e in 1:nne)
        {
          fit=Store.SPM.scallopedHH[[sc]][[h]][[e]]
          MLE=round(fit$par,2)
          std=round(sqrt(diag(solve(fit$hessian))),2)
          Nms=paste("log",names(MLE))
          
          Tab=as.data.frame(matrix(paste(MLE," (",std,")",sep=''),nrow=1))
          names(Tab)=Nms
          Tab=cbind(Species=capitalize(names(SPM.preds.scallopedHH)[sc]),Tab)
        }
      }
      Tab.par.estim.SPM.scallopedHH[[sc]]=Tab
    }
    Tab.par.estim.SPM.scallopedHH=do.call(rbind,Tab.par.estim.SPM.scallopedHH)
    setwd(hNdl.HH)
    fn.word.table(WD=getwd(),TBL=Tab.par.estim.SPM.scallopedHH,Doc.nm="estimates_SPM",caption=NA,paragph=NA,
                  HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                  Zebra='NO',Zebra.col='grey60',Grid.col='black',
                  Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
    
    
    rm(HR.o.scens)
    
  }
  
  #Display Catch-MSY
  if(Do.Ktch.MSY=="YES")
  {
    #r priors   
    fn.fig(paste(hNdl.HH,"Figure Prior_r",sep=""), 2000, 2000)
    smart.par(n.plots=1,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
    for(s in s)
    {
      NMs=capitalize(names(store.species)[s])
      if(NMs=="Low") NMs="Low resilience"
      if(NMs=="Very.low") NMs="Very low resilience"
      plot(density(rgamma(10000, shape = store.species[[s]]$r.prior$shape, rate = store.species[[s]]$r.prior$rate)),
           lwd=3,main='',xlab="",ylab="",cex.lab=2,cex.axis=1.15,col=1,xlim=c(0.05,.25),yaxt='n')
    }
    mtext(expression(paste(plain("Intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),1,0.5,cex=1.35,outer=T)
    mtext("Density",2,0,las=3,cex=1.35,outer=T)
    dev.off()
    
    YrS=sort(unique(Tot.ktch$finyear))
    

    #Relative biomass
    fn.fig(paste(hNdl.HH,"Figure 2_Biomass_Catch_MSY",sep=""), 1600, 2400)
    smart.par(n.plots=length(compact(SPM.preds.scallopedHH)),MAR=c(1.2,2,2,1.25),
              OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
    par(las=1,cex.axis=1.1)
    for(sc in 1: length(Scens))
    {
      Yrs=Store.stuff.scallopedHH[[sc]]$yrs
      Ktch_MSY_Rel.bio=with(Ktch_MSY.scallopedHH[[sc]],sweep(bt, 2, k, `/`))
      
      #Percentile   
      fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES",add.RP.txt="NO",
                         CEX=1,CX.AX=1.1,Ktch=Store.stuff.scallopedHH[[sc]]$Ktch,
                         PERCENTIL=95,Add.ktch="NO")
      mtext(Scens[sc],3,0)
    }
    mtext("Financial year",1,cex=1.2,line=0.75,outer=T)
    mtext("Relative biomass",2,cex=1.2,outer=T,las=3)
    #mtext(side = 4, line = 0.75, 'Total catch (tonnes)',las=3,outer=T,
    #      col=rgb(0.1,0.1,0.8,alpha=0.6),cex=1.2)
    dev.off()
    
    #MSY
    fn.fig(paste(hNdl.HH,"Figure MSY_Catch.MSY",sep=""), 1600, 2400)
    smart.par(n.plots=length(compact(SPM.preds.scallopedHH)),MAR=c(1.2,2,1.5,.1),
              OMA=c(2,2,.2,2.1),MGP=c(1,.62,0))
    par(las=1,cex.axis=1.1)
    for(sc in 1: length(Scens))
    {
      dummy=Ktch_MSY.scallopedHH[[sc]]$msy
      plot(density(dummy,adjust = 2),main="",ylab="")
      mtext(Scens[sc],3,cex=1)
      legend("right",paste(round(median(dummy))," tonnes",sep=""),
             bty='n',cex=1.1,title='Median MSY')
    }
    mtext("Catch (tonnes)",1,line=0.75,outer=T)
    mtext("Density",2,line=.75,outer=T,las=3)
    dev.off()
  }
  
  #Display aSPM
  if(Do.aSPM=="YES")
  {
    #Plot obs VS pred cpues  
    fn.fig(paste(Paz,"Fit_aSPM",sep=""),2400,2000)
    smart.par(n.plots=1,MAR=c(1.2,2,2,1.25),
              OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
    
    ii=match(all.yrs,Store.aSPM.scallopedHH$fish$year)
    iii=which(!is.na(Store.aSPM.scallopedHH$fish$cpue))
    fn.plt.cpue(ob=log(Store.aSPM.scallopedHH$fish$cpue[iii]),
                ob.CV=Store.aSPM.scallopedHH$fish$se[ii],
                pred=log(Store.aSPM.scallopedHH$quantities$PredCE[ii]),
                Convergence=Store.aSPM.scallopedHH$fit$convergence,
                pred.fit=log(Store.aSPM.scallopedHH$quantities$PredCE[iii]))
    
    legend("bottom",c("observed","predicted"),pch=19,cex=1.25,col=c("orange","black"),bty='n')
    dev.off()
    
    
    #Plot biomass  
    
    #1. Get median and percentiles    Store.aSPM.scallopedHH
    Med.biom.scallopedHH_aSPM=matrix(NA,nrow=NsimSS,ncol=length(all.yrs))
    for(i in 1:NsimSS)
    {
      dummy=Store.aSPM.scallopedHH$quantities.MC[[i]]$Deplete
      Med.biom.scallopedHH_aSPM[i,]=dummy[ii]
    }
    Bt=apply(Med.biom.scallopedHH_aSPM,2,function(x) quantile(x,probs=c(.025,0.5,.975)))  
    
    #2. Plot
    fn.fig(paste(hNdl.HH,"Figure 2_Biomass_aSPM",sep=""),2400,2000)
    smart.par(n.plots=1,MAR=c(1.2,2,2,1.25),OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
    par(las=1,cex.axis=1)
    DAT=t(Med.biom.scallopedHH_aSPM)
    
    fn.plt.bio.ktch(Yr=all.yrs,
                    Bt=Bt,
                    Bmsy=Med.biom.scallopedHH[[sc]][[h]][[e]]$Bmsy,
                    Ktch=Store.stuff.scallopedHH[[sc]]$Ktch,
                    CX.AX=1,
                    CX=1,
                    DAT=DAT,
                    Add.ktch="NO",
                    LOW="2.5%",
                    HIGH="97.5%")
    mtext("Financial year",1,cex=1.2,line=0.75,outer=T)
    mtext("Relative biomass",2,cex=1.2,outer=T,las=3)
    
    dev.off()
    
    
    #Output parameter estimates
    fit=Store.aSPM.scallopedHH$fit
    MLE=round(fit$par,2)
    std=round(sqrt(diag(solve(fit$hessian))),2)
    Nms=names(MLE)
    Tab=as.data.frame(matrix(paste(MLE," (",std,")",sep=''),nrow=1))
    names(Tab)=Nms
    setwd(hNdl.HH)
    fn.word.table(WD=getwd(),TBL=Tab,Doc.nm="estimates_aSPM",caption=NA,paragph=NA,
                  HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                  Zebra='NO',Zebra.col='grey60',Grid.col='black',
                  Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  }
}

#---Removed from CMSY ------
do.removed=FALSE
if(do.removed)
{
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
  
  
  #Current depletion of total biomass
  # fn.fig("Figure 3_Current.depletion_Catch_MSY", 2000, 2200)
  # smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  # for(s in 1: N.sp)
  # {
  #   Yrs=Store.stuff[[s]]$yrs
  #   Current=Yrs[length(Yrs)]
  #   Ktch_MSY_Rel.bio=with(store.species[[s]]$KTCH.MSY$BaseCase,bt,k)
  #   Ktch_MSY_current_yr=Ktch_MSY_Rel.bio[length(Yrs),]
  #   NMs=names(store.species)[s]
  #   if(NMs=="Very.low") NMs="Very low"
  #   density.fun2(what=Ktch_MSY_current_yr,MAIN=NMs)
  # }
  # mtext(paste(Current,"Relative biomass"),1,line=0.25,cex=1.35,outer=T)
  # mtext("Probability",2,0.25,las=3,cex=1.35,outer=T)
  # dev.off()
  
  
  # #Fishing mortality
  # fn.fig("Fishing_mortality_Base case", 2000, 2000)
  # smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  # for(s in 1: N.sp)
  # {
  #   YYrs=Store.stuff[[s]]$yrs
  #   Ktch_MSY_Rel.bio=store.species[[s]]$KTCH.MSY$BaseCase$Fish.mort
  #
  #   #Percentile
  #   fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="NO",add.RP.txt="NO",CEX=1.2,
  #                      CX.AX=1.2,Ktch=Store.stuff[[s]]$Ktch)
  #   #if(s==1) legend("topleft",c("50%","75%","100%"),fill=COLS,bty='n',cex=1.25)
  #   NMs=names(store.species)[s]
  #   if(NMs=="Very.low") NMs="Very low"
  #   mtext(NMs,3,0)
  #
  #
  #   #   #Geometric mean
  #   # Ktch_MSY_rel_bt_mean=Ktch_MSY_rel_bt_lowSE=Ktch_MSY_rel_bt_upSE=nrow(Ktch_MSY_Rel.bio)
  #   # for(nr in 1:nrow(Ktch_MSY_Rel.bio))
  #   # {
  #   #   Ktch_MSY_rel_bt_mean[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])))
  #   #   Ktch_MSY_rel_bt_upSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) + 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
  #   #   Ktch_MSY_rel_bt_lowSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) - 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
  #   # }
  #   # plot(Yrs,Ktch_MSY_rel_bt_mean,cex=1.1,pch=19,col=CL.mean,ylim=c(0,max(Ktch_MSY_rel_bt_upSE,na.rm=T)),xaxt='n',xlab="",ylab="",
  #   #      main=names(store.species)[s],cex.axis=1.15)
  #   # segments(Yrs,Ktch_MSY_rel_bt_lowSE,Yrs,Ktch_MSY_rel_bt_upSE,col=CL)
  #   # axis(1,Yrs,labels=F,tck=-0.015)
  #   # axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
  # }
  # mtext("Financial year",1,0.5,cex=1.35,outer=T)
  # mtext(expression(paste(plain("Fishing mortality (year") ^ plain("-1"),")",sep="")),2,0,las=3,cex=1.35,outer=T)
  # dev.off()
  
  
  #Catch and MSY
  # fn.fig("Figure4_CatchMSY_Plots", 2000, 2200)
  # smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  # for(s in 1: N.sp)
  # {
  #   yr=store.species[[s]]$Catch$finyear
  #   ct=store.species[[s]]$Catch$Total.ktch
  #   mean_ln_msy=as.numeric(store.species[[s]]$KTCH.MSY$BaseCase$"geom. mean MSY (tons)")
  #   msy=store.species[[s]]$KTCH.MSY$BaseCase$msy
  #   Mean.MSY=exp(mean(log(msy)))
  #   Mean.MSY_UP=exp(mean(log(msy)) + 1.96 * sd(log(msy)))
  #   Mean.MSY_LOW=exp(mean(log(msy)) - 1.96 * sd(log(msy)))
  #   NMs=names(store.species)[s]
  #   if(NMs=="Low") NMs="Low resilience complex"
  #   if(NMs=="Very.low") NMs="Very low resilience complex"
  #
  #   plot(yr, ct, type="l", ylim = c(0, 1.01*max(c(ct,Mean.MSY_UP))), xlab = "", ylab = "",
  #        lwd=2,cex.lab=2,cex.axis=1.15)
  #   mtext(NMs,3,0)
  #   all.yrs=c(yr[1]-2,yr,yr[length(yr)]+2)
  #   polygon(c(all.yrs,rev(all.yrs)),  c(rep(Mean.MSY_LOW,length(all.yrs)),rep(Mean.MSY_UP,length(all.yrs))),
  #           col=rgb(.1,.1,.1,alpha=.15),border="transparent")
  #   abline(h=Mean.MSY,col="grey50", lwd=2.5,lty=3)
  #   if(s==9)legend("topright",c("MSY (?1.96 SE)"),bty='n',col=c("grey50"),lty=3,lwd=2.5,cex=1.25)
  #
  # }
  # mtext("Financial year",1,0.25,cex=1.35,outer=T)
  # mtext("Total catch (tonnes)",2,0.35,las=3,cex=1.35,outer=T)
  # dev.off()
  
  
  #Sensitivity tests
  # Biomass
  # fn.plot.percentile.sens=function(DAT,YR,ADD.prob,add.RP.txt,CEX,CX.AX,AdYXs,AdXXs)
  # {
  #   #50% of data
  #   Nper=(100-50)/2
  #   LOW.50=Low.percentile(Nper,DAT)
  #   UP.50=High.percentile(Nper,DAT)
  #
  #   #75% of data
  #   Nper=(100-75)/2
  #   LOW.75=Low.percentile(Nper,DAT)
  #   UP.75=High.percentile(Nper,DAT)
  #
  #   #100% of data
  #   Nper=(100-100)/2
  #   LOW.100=Low.percentile(Nper,DAT)
  #   UP.100=High.percentile(Nper,DAT)
  #
  #   #construct polygons
  #   Year.Vec <-  fn.cons.po(YR,YR)
  #   Biom.Vec.50 <- fn.cons.po(LOW.50,UP.50)
  #   Biom.Vec.75 <- fn.cons.po(LOW.75,UP.75)
  #   Biom.Vec.100 <-fn.cons.po(LOW.100,UP.100)
  #
  #
  #   #plot
  #   plot(YR,UP.100,ylim=c(0,max(UP.100)),type="l",ylab="",xlab="",yaxt='n',xaxt='n',col='transparent',cex.axis=CX.AX)
  #
  #   polygon(Year.Vec, Biom.Vec.100, col = COLS[3], border = "grey20")
  #   polygon(Year.Vec, Biom.Vec.75, col = COLS[2], border = "grey20")
  #   polygon(Year.Vec, Biom.Vec.50, col = COLS[1], border = "grey20")
  #
  #
  #   #add probs
  #   if(ADD.prob=="YES")
  #   {
  #     add.probs(id.yr=match(Current,YR),YR,DAT,UP.100,LOW.100,SRT=0,CEX)
  #
  #     abline(h=B.target,lwd=1.5,col='grey45',lty=3)
  #     abline(h=B.threshold,lwd=1.5,col='grey45',lty=3)
  #     abline(h=B.limit,lwd=1.5,col='grey45',lty=3)
  #
  #     if(add.RP.txt=="YES")
  #     {
  #       text(YR[4],B.target,"Target",pos=3,cex=1.1)
  #       text(YR[4],B.threshold,"Threshold",pos=3,cex=1.1)
  #       text(YR[4],B.limit,"Limit",pos=3,cex=1.1)
  #     }
  #
  #   }
  #   if(AdYXs=="YES")axis(2,at=seq(0,1,.2),labels=seq(0,1,.2),tck=-0.05,las=1,cex.axis=CX.AX)
  #   axis(1,at=YR,labels=F,tck=-0.025)
  #   axis(1,at=seq(YR[1],YR[length(YR)],5),labels=F,tck=-0.05,cex.axis=CX.AX)
  #   if(AdXXs=="YES")axis(1,at=seq(YR[1],YR[length(YR)],5),labels=seq(YR[1],YR[length(YR)],5),tck=-0.05,cex.axis=CX.AX)
  # }
  # fn.fig("Sensitivity_Biomass_relative", 2800,800)
  # par(mfcol=c((Nscen-1),N.sp),mar=c(1,1,.15,.15),oma=c(2,3.5,.95,.25),las=1,mgp=c(1,.5,0))
  # for(s in 1: N.sp)
  # {
  #   Yrs=store.species[[s]]$Catch$finyear
  #   NMs=names(store.species)[s]
  #   if(NMs=="Scalloped hammerhead") NMs="Scalloped hh"
  #   if(NMs=="Smooth hammerhead") NMs="Smooth hh"
  #   for(sc in 2:Nscen)
  #   {
  #     Ktch_MSY_Rel.bio=store.species[[s]]$KTCH.MSY[[sc]]$bt.rel
  #     SCNE.nm=names(store.species[[s]]$KTCH.MSY)[sc]
  #     AddY="NO"
  #     if(s==1) AddY="YES"
  #     AddX='NO'
  #     if(sc==3) AddX="YES"
  #     fn.plot.percentile.sens(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES",
  #                             add.RP.txt="NO",CEX=.7,CX.AX=.8,AdYXs=AddY,AdXXs=AddX)
  #
  #
  #     #if(s==N.sp & sc==2) legend("bottom",c("50%","75%","100%"),fill=COLS,bty='n',cex=.65,horiz=T)
  #     if(s==1) mtext(SCNE.nm,2,1.5)
  #     if(NMs=="Low") NMs="Low resilience"
  #     if(NMs=="Very.low") NMs="Very low resilience"
  #     if(sc==2)mtext(NMs,3,0,cex=.75)
  #   }
  # }
  # mtext("Financial year",1,1,cex=1.35,outer=T)
  # mtext("Relative biomass",2,2,las=3,cex=1.35,outer=T)
  # dev.off()
  
  #MSY
  # Exprt.MSY=vector('list',length(N.sp))
  # fn.fig("Sensitivity_MSY", 1000, 2400)
  # par(mfrow=c(N.sp,Nscen),mar=c(2,2,.25,.35),oma=c(1.75,2.75,.95,.25),las=1,mgp=c(1,.5,0))
  # for(s in 1: N.sp)
  # {
  #   Yrs=store.species[[s]]$Catch$finyear
  #   NMs=names(store.species)[s]
  #   yr=store.species[[s]]$Catch$finyear
  #   ct=store.species[[s]]$Catch$Total.ktch
  #   dummyMSY=vector('list',length(Nscen))
  #   for(sc in 1:Nscen)
  #   {
  #     mean_ln_msy=as.numeric(store.species[[s]]$KTCH.MSY[[sc]]$"geom. mean MSY (tons)")
  #     msy=store.species[[s]]$KTCH.MSY[[sc]]$msy
  #     Mean.MSY=exp(mean(log(msy)))
  #     Mean.MSY_UP=exp(mean(log(msy)) + 1.96 * sd(log(msy)))
  #     Mean.MSY_LOW=exp(mean(log(msy)) - 1.96 * sd(log(msy)))
  #
  #     plot(yr, ct, type="l", ylim = c(0, 1.01*max(c(ct,Mean.MSY_UP))), xlab = "", ylab = "",
  #          lwd=2,cex.lab=2,cex.axis=1)
  #     all.yrs=c(yr[1]-2,yr,yr[length(yr)]+2)
  #     polygon(c(all.yrs,rev(all.yrs)),  c(rep(Mean.MSY_LOW,length(all.yrs)),rep(Mean.MSY_UP,length(all.yrs))),
  #             col=rgb(.1,.1,.1,alpha=.15),border="transparent")
  #     abline(h=Mean.MSY,col="grey50", lwd=2.5,lty=3)
  #     #if(s==4)legend("topright",c("MSY (?1.96 SE)"),bty='n',col=c("grey50"),lty=3,lwd=2.5,cex=1.25)
  #     SCNE.nm=names(store.species[[s]]$KTCH.MSY)[sc]
  #     if(s==1) mtext(SCNE.nm,3,0)
  #     if(NMs=="Low") NMs="Low resilience"
  #     if(NMs=="Very.low") NMs="Very low resilience"
  #     if(NMs=="Scalloped hammerhead") NMs="Scalloped hh"
  #     if(NMs=="Smooth hammerhead") NMs="Smooth hh"
  #
  #     if(sc==1)mtext(NMs,2,2,las=3,cex=.8)
  #     dummyMSY[[sc]]=data.frame(Species=NMs,Scenario=SCNE.nm,LOW_CI=Mean.MSY_LOW,Mean.MSY=Mean.MSY,Mean.UP_CI=Mean.MSY_UP)
  #   }
  #
  #   Exprt.MSY[[s]]=do.call(rbind,dummyMSY)
  # }
  # mtext("Financial year",1,0.25,cex=1.35,outer=T)
  # mtext("Total catch (tonnes)",2,1.25,las=3,cex=1.35,outer=T)
  # dev.off()
  # write.csv(do.call(rbind,Exprt.MSY),"MSY_estimates.csv",row.names=F)
  
}