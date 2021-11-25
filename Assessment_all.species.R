# ------ Script for running stock assessments for WA sharks---- ###################

#missing:  
#     Include ALL species in final risk scoring
#     review Smooth HH cpue and mako cpue...; SPM Tiger fit
#     Milk shark SPM, hitting upper K boundary, no trend in cpue, crap Hessian, too uncertain....mention in text...
#     aSPM: finish running for all species; issues with Tiger cpue fit...
#     Size-based Catch curve (for some species there's NSF size compo, not used at the moment)
#     Implement simple SS?
#     set up integrated model for dusky and sandbar
#     set up SS3 model for indicator species
#     Use .CPUE_Observer_TDGDLF.csv in likelihoods?
#     Double check that TwT calcultion uses TL and not FL, be consistent

#Notes:
#     1. The PSA filters out species from further analyses thru 'Criteria for selecting what species to assess'
#     2. A range of assessment methods are used depending on data availability: 
#               . Integrated size- and sex- structured models,
#               . SPM,
#               . Catch-MSY (Martell & Froese 2012)
#               . aSPM(Haddon SimpleSA())
#     3. Total (reconstructed) catches are used (commercial, recreational and TDGDLF discards); recons_NT_catch.csv' only includes dusky and sandbar
#     4. Assumption: If catches have never been >1% carrying capacity, then it's in unexploited 
#                    status so catch series have no information on productivity

#Steps: 
#     1. Update year of assessment in '1. DEFINE GLOBALS'  
#     2. Define arguments used in each of the shark species/species complex assessed.
#     3. Bring in updated available data and input parameters
#     4. Determine which species to assessed based on PSA
#     5. Run relevant population models according to data availability
#     6. Generate relevant outputs

# ------ Header---- 
rm(list=ls(all=TRUE))
options(dplyr.summarise.inform = FALSE)
options(stringsAsFactors = FALSE) 
options(ggrepel.max.overlaps = Inf)
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#library(ReporteRs)   #MISSING: replace by Officer from Source.script git_other
library(rlist)
library(MASS)
library(plotrix)
library(PBSmapping)
library(tidyverse)
library(dplyr)
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
library(zoo)
library(MCDA)
library(sfsmisc)   # p values from rlm
library(data.table)
library(ggpmisc)
library(ggpubr)
library(doParallel)
library(flextable)
library(officer)
library(fishmethods)
#devtools::install_github("cfree14/datalimited2")
library(datalimited2)
#devtools::install_github("jabbamodel/JABBA")
library(JABBA)

source.hnld=handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/")
fn.source=function(script)source(paste(source.hnld,script,sep=""))
fn.source("fn.fig.R")
fn.source("Catch_MSY.R")
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R"))
#source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R"))  #MISSING: replace by Officer from Source.script git_other

smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
colfunc <- colorRampPalette(c("red","yellow","springgreen","royalblue"))
fun.find.in.list=function(x,Drop)   #drop stuff from list
{
  x=x%>%discard(is.null)
  if(!is.null(Drop))x=x[-match(Drop,names(x))]
  return(x)
}

fn.get.stuff.from.list=function(lista,stuff) lapply(lista, function(x) x[[stuff]])# get stuff from list

source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/ggplot.themes.R'))  #my themes

#---1.  DEFINE GLOBALS----- 

#Year of assessment 
Year.of.assessment=2021
Asses.year=Year.of.assessment
AssessYr=Year.of.assessment


#Last complete financial year of catches
Last.yr.ktch="2019-20"

#Define is doing stand-alone sawfish assessment
do.sawfish=FALSE

#Model run
First.run="NO"
#First.run="YES"  #set to yes to create all model input data sets and data presentation
Modl.rn="standard"   #for annual assessments
#Modl.rn='first'    #for paper


#Define if calculating r, steepness and size-based catch curve
do.r.prior=do.steepness=FALSE    #it's a one off, doesn't need update
if(First.run=="YES") do.Size.based.Catch.curve=TRUE  else #needs update with each new year of data
                     do.Size.based.Catch.curve=FALSE
#Define if exporting figures as jpeg or tiff (creation of RAR requires jpeg)
Do.tiff="YES" 
Do.jpeg="NO"


#Catch in tonnes?  
KTCH.UNITS="TONNES" 
#KTCH.UNITS="KGS"    
if(KTCH.UNITS=="KGS") unitS=1000
if(KTCH.UNITS=="TONNES") unitS=1


#Criteria for selecting what species to assess
Min.yrs=5
if(KTCH.UNITS=="KGS") Min.ktch=5000 
if(KTCH.UNITS=="TONNES") Min.ktch=5
prop.disc.ER=.4  #proportion of vessels discarding eagle rays in last 5 years (from catch and effort returns)


  #PSA
PSA.min.tons=5
PSA.min.years=Min.yrs
PSA.max.ton=50
Low.risk=2.64  #risk thresholds from Micheli et al 2014
medium.risk=3.18


#Size composition
  #Initial bin size
MN.SZE=0
#MN.SZE="size.at.birth"

  #size bin
TL.bins.cm=5


  #Minimun number of samples of size composition
Min.obs=10  #at least 10 observations
Min.shts=5  #from at least 5 shots


# Global arguments for indicator species 
if(First.run=="YES") run.all.scenarios="YES" #What scenarios to run?
if(First.run=="NO") run.all.scenarios="NO"  
if(First.run=="NO") run.future="YES" #Future projection scenarios
if(First.run=="YES") run.future="NO"
Show.yrs="DATA"   #What years to show in model outputs
#Show.yrs="FUTURE"  #show data plus future projections
Present.in.log="NO"  #present cpue in log space or normal space

  #Size-based model
Do.zero.Ktch="NO" #do zero catch projections of size-based model?
Do.sim.test="NO" #do simulation testing of size-based model?
mcmc.do="YES"
mcmc.n=1e6
mcmc.thin=10
Plus.gp.size=1.25  #add 25% to max size make sure no accumulation of survivals in last size class
MIN.DAYS.LARGE=30  #minimum number of days at liberty for movement
Sim.trans.Mat="NO" #show simulated size transition matrix
add.conv.tag="YES" #define if using conventional tagging and effort in model
fn.ji=function(a) jitter(a,factor=5)  #Add randomness to pin values   
Move.mode="Individual-based"  #select type of movement modellling approach
#Move.mode="Population-based"
conv.tag.all="YES"  #Combine all size classes for tagging model
Fz.off=-22
Do.from.at.len.to.at.age="NO" #Convert from length to age
basecase="Size-based"  #select base case model
#basecase="Age-based"  
do.cols="YES"  #Do Scenario comparison in colors or greyscale
#do.cols="NO"
size.likelihood="Dirichlet" #Select type of size composition Likelihood
#size.likelihood="Multinomial"
size.sex.comb="NO" #Combined size composition?
size.comp.prop="YES" #Size compostion as proportions?
yrs.future=5 #Number of years for future projections
effective.n=300 #Effective sample size size composition
maxF=3.0 #Maximum possible F
Add_Finit_prior=0  ##Prior for initial fishing mortality
#Add_Finit_prior=1  #add prior
arguments=""  #run without arguments
#arguments=" -est -sim "  #run with these arguments
do.MSY="NO"   #Do MSY calculation using Base Case model
msy.yrs=100   #MSY years  
msy.sd.rec=0  #no recruitment deviations (in log space)
#msy.sd.rec=0.05
msy.sims=100
f.vec=seq(0,.2,by=.01)

  #SPM
r_max=0.5     #Include a prior for r in surplus production. This is the max reported value of r for sharks (blue shark)
Add.r.prior=0   #no r prior 
#Add.r.prior=1   #use r prior


# Global arguments for non-indicator species
Asses.Scalloped.HH=FALSE  #2020 scalloped HH assessment

# Minimun number of annual observations in analysis of changes in size
Min.annual.obs=300

# Control which assessment methods to implement
Do.SPM="YES"
Do.Ktch.MSY="YES"
Do.aSPM="YES"


# Reference points
#note: Historically, single reference point for the fishery (biomass at 40% unexploited conditions)
B.threshold=0.5  #Bmys
Tar.prop=1.3    #target and limit proportions of Bmsy. source: Haddon et al 2014. Technical
Lim.prop=0.5    #   Reviews of Formal Harvest Strategies.
B.target=Tar.prop*B.threshold
B.limit=Lim.prop*B.threshold

    #Empirical reference points
Fmsy.emp=function(M) 0.41*M     #Zhou et al 2012 but see Cortes & Brooks 2018
SPR.thre=0.3   #Punt 2000 Extinction of marine renewable resources: a demographic analysis. 
SPR.tar=0.4                # Population Ecology 42, 


# Define if using effort
add.effort="NO"    
What.Effort="km.gn.hours"  #What effort to display?
#What.Effort="km.gn.days" 


# Define if assessing 'other species'
do.other.ass=TRUE


# Define if doing other species paper
do.other.ass.paper=FALSE  

# Assumed PCM for reconstructed discards in TDGLDF
TDGLDF.disc.assumed.PCM="BaseCase" 
#TDGLDF.disc.assumed.PCM="100%" 


#---2.  Catch and effort data section-----   
Dat.repository=handl_OneDrive('Analyses/Data_outs/')  #location where all input data are located

#2.1. Species  
All.species.names=read.csv(handl_OneDrive("Data/Species_names_shark.only.csv"))
All.species.names=All.species.names%>%
                    mutate(Name=tolower(Name))%>%
                    filter(!(SPECIES==18014 & Name=="australian blacktip shark"))%>%
                    rename(SNAME=Name)%>%
                    arrange(SPECIES)
Shark.species=5001:24900
School.shark= 17008
Indicator.species=c(17001,17003,18003,18007)
names(Indicator.species)=c("gummy shark","whiskery shark","dusky shark","sandbar shark")
Shar_other=22999
Scien.nm=All.species.names[,c('SPECIES','Scien.nm')]

#define what species are assessed elsewhere
assessed.elsewhere=c("white shark","school shark","spot-tail shark",
                     "blacktips","unidentified sharks")  



#2.2. Import Total effort
fn.in=function(NM) fread(paste(Dat.repository,NM,sep=""),data.table=FALSE)
  #TDGDLF
Effort.monthly=fn.in(NM='Annual.total.eff.days.csv')          #all years
Effort.monthly_hours=fn.in(NM='Annual.total.eff.hours.csv')   #all years
Effort.monthly_blocks=fn.in("Effort.monthly.csv")
Effort.daily_blocks=fn.in("Effort.daily.csv")
  #NSF
Effort.monthly.north=fn.in(NM='Annual.total.eff_NSF.csv') 
Effort.monthly.north_blocks=fn.in("Effort.monthly.NSF.csv")
Effort.daily.north_blocks=fn.in("Effort.daily.NSF.csv")


#2.3. Import Total Catch
Mode <- function(x)
{
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
fn.import.catch.data=function(KTCH.UNITS)
{
  #2.1. Commercial
    #2.1.1 Catch_WA Fisheries
  
  #..Historic
  Historic.ktch=fn.in(NM='recons_Hist.expnd.csv')%>%filter(!FINYEAR=='1975-76')
  
  #..Ammended reported commercial catch including discards
  Data.monthly=fn.in(NM='recons_Data.monthly.csv')
  Data.monthly.north=fn.in(NM='recons_Data.monthly.north.csv')
  
  #..TEPS
  Sawfish_Pilbara.ktch=fn.in(NM='recons_Pilbara_sawfish.ktch.csv')
  Greynurse.ktch=fn.in(NM='recons_Greynurse.ktch.csv')
  TEPS_dusky=fn.in(NM='recons_TEPS_dusky.csv')
  Sawfish_KPTF.ktch=fn.in(NM='recons_Kimberly.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_NBPTF.ktch=fn.in(NM='recons_Nickol.Bay.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_OPTF.ktch=fn.in(NM='recons_Onslow.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_SBPTF.ktch=fn.in(NM='recons_Shark.Bay.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_EGPTF.ktch=fn.in(NM='recons_Exmouth.Gulf.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_BPTF.ktch=fn.in(NM='recons_Broome.Prawn.Trawl_sawfish.ktch.csv')
  
  #..Droplines Western Rock Lobster
  WRL.ktch=fn.in(NM='Wetline_rocklobster.csv')
  
  #..Discards from TDGDLF
  if(TDGLDF.disc.assumed.PCM=="BaseCase") discard_TDGDLF=fn.in(NM='recons_discard_TDGDLF.csv')
  if(TDGLDF.disc.assumed.PCM=="100%") discard_TDGDLF=fn.in(NM='recons_discard_TDGDLF_100.PCM')
  
  
    #2.1.2. Catch of non WA Fisheries
  
  #..Taiwanese gillnet and longline
  Taiwan.gillnet.ktch=fn.in(NM='recons_Taiwan.gillnet.ktch.csv')
  Taiwan.longline.ktch=fn.in(NM='recons_Taiwan.longline.ktch.csv')
  Taiwan.gillnet.ktch$Method="Pelagic.gillnet"
  Taiwan.longline.ktch$Method="Longline"
  Taiwan=rbind(Taiwan.longline.ktch,Taiwan.gillnet.ktch)%>%
    mutate(Region="North")%>%
    mutate(year=as.numeric(substr(FINYEAR,1,4)))
  
  #..Indonesian illegal fishing in Australia waters
  Indo_total.annual.ktch=fn.in(NM='recons_Indo.IUU.csv')
  Indo_total.annual.ktch=Indo_total.annual.ktch%>%filter(!is.na(SPECIES))
  
  #..AFMA's GAB & SBT fisheries
  GAB.trawl_catch=fn.in(NM='recons_GAB.trawl_catch.csv') 
  WTBF_catch=fn.in(NM='recons_WTBF_catch.csv')
  
  #..SA Marine Scalefish fishery
  Whaler_SA=fn.in(NM='recons_Whaler_SA.csv') 
  
  #..NT catches
  NT_catch=fn.in(NM='recons_NT_catch.csv') 
  
  
  #2.2. WA Recreational catch
  Rec.ktch=fn.in(NM='recons_recreational.csv')  
  Rec.ktch=Rec.ktch%>%mutate(Region=ifelse(zone%in%c('Gascoyne','North Coast'),'North','South'),
                             year=as.numeric(substr(FINYEAR,1,4)),
                             Common.Name=tolower(Common.Name))
  
  
  #Combine data sets
    # TDGDLF and NSF
  Data.monthly$Region="South"
  Data.monthly.north$Region="North"
  Data.monthly$Data.set="Data.monthly"
  Data.monthly.north$Data.set="Data.monthly.north"
  Tot.ktch=rbind(subset(Data.monthly,select=c(SPECIES,FishCubeCode,FINYEAR,Region,Data.set,LIVEWT.c,BLOCKX,METHOD,zone)),
                 subset(Data.monthly.north,select=c(SPECIES,FishCubeCode,FINYEAR,Region,Data.set,LIVEWT.c,BLOCKX,METHOD,zone)))%>%
    left_join(All.species.names,by='SPECIES')%>%     #add species names
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%    
    mutate(Name=SNAME)
  
    #Add Recreational
  a=Rec.ktch%>%
    rename(finyear=year)%>%
    mutate(BLOCKX=NA,
           zone=case_when(zone%in%c("North Coast","Gascoyne Coast")~'North',
                          zone=="West Coast"~'West',
                          zone=="South Coast"~'Zone1'),
           Common.Name=ifelse(Common.Name=="dogfishes","spurdogs",
                       ifelse(Common.Name=="greynurse shark","grey nurse shark",
                       ifelse(Common.Name=="thresher shark","thresher sharks",
                       ifelse(Common.Name=="bronze whaler","copper shark",
                       ifelse(Common.Name=="dusky whaler","dusky shark",
                       ifelse(Common.Name=="gummy sharks","gummy shark",
                       ifelse(Common.Name=="tawny shark","tawny nurse shark",
                       ifelse(Common.Name=="sawfishes","sawfish (general)",
                       ifelse(Common.Name=="australian blacktip shark","australian blacktip",
                       Common.Name))))))))))%>%
    left_join(All.species.names,by=c('Common.Name'='SNAME'))%>%
    mutate(SNAME=Common.Name,
           Name=Common.Name,
           METHOD='Rec.line',
           Data.set="Recreational",
           FishCubeCode="Recreational")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
    #Add Discards from TDGDLF
  a=discard_TDGDLF%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Name=SNAME,
           Type="Discards_TDGDLF",
           Data.set="Discards_TDGDLF",
           METHOD='GN',
           Region='South',
           FishCubeCode="Discards_TDGDLF",
           finyear=as.numeric(substr(FINYEAR,1,4)))%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
    #Add Taiwanese
  a=Taiwan%>%
    left_join(All.species.names,by='SPECIES')%>%
    rename(finyear=year)%>%
    mutate(BLOCKX=NA,
           Name=SNAME,
           Type="Taiwan",
           Data.set="Taiwan",
           METHOD=Method,
           FishCubeCode="Taiwan")%>%
    dplyr::select(names(Tot.ktch))%>%
    filter(SPECIES%in%unique(Tot.ktch$SPECIES))
  Tot.ktch=rbind(Tot.ktch,a)
  
    #Add Indonesian fishing incursions
  a=Indo_total.annual.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="Indonesia",
           Data.set="Indonesia",
           METHOD=NA,
           FishCubeCode="Indo")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
    #Add TEP interactions     
  a=Sawfish_Pilbara.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.Pilbara",
           Data.set="TEP_sawfish.Pilbara",
           METHOD='TW',
           FishCubeCode="PFT")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)

  a=Sawfish_KPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.KPTF",
           Data.set="TEP_sawfish.KPTF",
           METHOD='TW',
           FishCubeCode="KP")%>%   
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_NBPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.NBPTF",
           Data.set="TEP_sawfish.NBPTF",
           METHOD='TW',
           FishCubeCode="NBP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_OPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.OPTF",
           Data.set="TEP_sawfish.OPTF",
           METHOD='TW',
           FishCubeCode="OP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_SBPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.SBPTF",
           Data.set="TEP_sawfish.SBPTF",
           METHOD='TW',
           FishCubeCode="SBP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_EGPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.EGPTF",
           Data.set="TEP_sawfish.EGPTF",
           METHOD='TW',
           FishCubeCode="EGP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_BPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.BPTF",
           Data.set="TEP_sawfish.BPTF",
           METHOD='TW',
           FishCubeCode="BP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)

  a=Greynurse.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_greynurse",
           Data.set="TEP_greynurse",
           METHOD='GN',
           FishCubeCode="TEP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=TEPS_dusky%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_dusky",
           Data.set="TEP_dusky",
           METHOD='GN',
           FishCubeCode="TEP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
    #Add WRL        
  a=WRL.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="WRL",
           Data.set="WRL",
           METHOD='LL',
           FishCubeCode="WRL")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  
    #Add historic  
  a=Historic.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Type="Commercial",
           Data.set="Historic",
           METHOD=NA,
           FishCubeCode="Historic",
           SNAME=ifelse(SNAME%in%c('southern sawshark','common sawshark'),'sawsharks',SNAME),
           SPECIES=ifelse(SPECIES%in%c(23001,23002),23900,SPECIES),
           Name=SNAME)%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
    #Add GAB
  a=GAB.trawl_catch%>%
    mutate(SPECIES=ifelse(SPECIES==5000,5002,
                   ifelse(SPECIES==23000,23900,
                   ifelse(SPECIES==12901,12000,
                   ifelse(SPECIES==24000,24900,
                   ifelse(SPECIES==19000,19004,
                   SPECIES))))))%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="GAB",
           Data.set="GAB.trawl",
           METHOD="trawl",
           FishCubeCode="GAB")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
    #Add WTBF
  a=WTBF_catch%>%
    mutate(SPECIES=ifelse(SPECIES==12901,12000,
                   ifelse(SPECIES==18901,18014,
                   ifelse(SPECIES==19000,19001,
                   SPECIES))))%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="WTB",
           Data.set="WTBF",
           METHOD="line",
           FishCubeCode="WTB")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
    #Add NT_catch
  a=NT_catch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="NT",
           Data.set="NT",
           METHOD="line",
           FishCubeCode="NT")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  
    #Add Whaler_SA (SA Marine Scalefish Fishery)
  a=Whaler_SA%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="SA MSF",
           Data.set="SA MSF",
           METHOD="line",
           FishCubeCode="SA MSF")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  
  #set catch in tonnes
  if(KTCH.UNITS=="TONNES")Tot.ktch=Tot.ktch%>%mutate(LIVEWT.c=LIVEWT.c/1000) 
  
  #fix some species names
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%19003,"winghead hammerhead",SNAME),
           Name=ifelse(SPECIES%in%19003,'winghead hammerhead',Name),
           Scien.nm=ifelse(SPECIES%in%19003,'eusphyra blochii',Scien.nm))
  
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SNAME=='sawfish (general)',"sawfish",SNAME),
           Scien.nm=ifelse(Name=='sawfish (general)','Pristidae',Scien.nm),
           SPECIES=ifelse(Name=='sawfish (general)',25000,SPECIES),
           Name=ifelse( Name=='sawfish (general)','sawfish',Name))
  
  Tot.ktch=Tot.ktch%>%
    mutate(SPECIES=ifelse(Name=='rays & skates',31000,SPECIES),
           Name=ifelse( Name=='skates','rays & skates',Name))
  
  Tot.ktch=Tot.ktch%>%
    mutate(Scien.nm=ifelse(Name=="australian blacktip",'Carcharhinus tilstoni',Scien.nm))

  Tot.ktch=Tot.ktch%>%
    mutate(SPECIES=ifelse(SPECIES==20904,20906,SPECIES),
           SNAME=ifelse(SPECIES==20906,"roughskin shark",SNAME),
           Scien.nm=ifelse(SPECIES==20906,'Centroscymnus spp.',Scien.nm),
           Name=ifelse( SPECIES==20906,'roughskin shark',Name))
  
  
  #list of all species caught
  Table1=Tot.ktch%>%  
    filter(!SPECIES%in%c(31000,22999))%>%
    group_by(Name,FishCubeCode)%>%
    summarise(Catch.tons=sum(LIVEWT.c))%>%
    spread(FishCubeCode,Catch.tons,fill=0)%>%
    left_join(Tot.ktch%>%
                distinct(SPECIES, .keep_all = T)%>%
                arrange(SPECIES)%>%
                dplyr::select(Name,Scien.nm),by="Name")%>%
    mutate(Name=capitalize(Name))%>%
    relocate(where(is.numeric), .after = where(is.character))

  
  #Combine certain species with variable reporting resolution
  unik=unique(Tot.ktch$Name)
      #Catsharks
  this=unik[grep("catshark",unik)]
  this=this[-match("brown-banded catshark", this)] #not a catshark
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(Name%in%this,"catsharks",SNAME),
           Scien.nm=ifelse(Name%in%this,'Scyliorhinidae',Scien.nm),
           SPECIES=ifelse(Name%in%this,15000,SPECIES),
           Name=ifelse(Name%in%this,'catsharks',Name))
  
      #Sawsharks
  this=c(23000,23001,23002,23900)
  Tot.ktch=Tot.ktch%>%
      mutate(SNAME=ifelse(SPECIES%in%this,"sawsharks",SNAME),
             Name=ifelse(SPECIES%in%this,'sawsharks',Name),
             Scien.nm=ifelse(SPECIES%in%this,'Pristiophoridae',Scien.nm),
             SPECIES=ifelse(SPECIES%in%this,23900,SPECIES),
             Name=ifelse(Name=='sawshark','sawsharks',Name))
  
      #Threshers
  this=c(12000,12001)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'thresher sharks',SNAME),
           Name=ifelse(SPECIES%in%this,'thresher sharks',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Alopidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,12000,SPECIES))   
  
      #Wobbegongs
  this=c(13000, 13001,13003,13011,13012,13017,13022,13016,13021)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'wobbegongs',SNAME),
           Name=ifelse(SPECIES%in%this,'wobbegongs',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Orectolobidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,13000,SPECIES)) 

      #Angelsharks
  this=c(24900, 24001,24002)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'angel sharks',SNAME),
           Name=ifelse(SPECIES%in%this,'angel sharks',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Squatinidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,24900,SPECIES))  

      #Stingrays
  this=c(35001, 35000)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'smooth stingray',SNAME),
           Name=ifelse(SPECIES%in%this,'smooth stingray',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Bathytoshia brevicaudata',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,35001,SPECIES))  
  
      #Banjo rays
  this=c(26999,27001,27007,27011,27909)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'banjo rays',SNAME),
           Name=ifelse(SPECIES%in%this,'banjo rays',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Trygonorrhinidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,26999,SPECIES))     
  
      #Wedgefishes
  this=c(26000,26001,26002)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'wedgefishes',SNAME),
           Name=ifelse(SPECIES%in%this,'wedgefishes',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Rhinidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,26000,SPECIES))       
  
      #Rays and skates  
  Tot.ktch=Tot.ktch%>%
    mutate(Name=ifelse(SPECIES==31000,'rays & skates',Name),
           Scien.nm=ifelse(SPECIES==31000,'Rajidae & Arhynchobatidae',Scien.nm),
           SNAME=ifelse(SPECIES==31000,'rays & skates',SNAME))

      #Australian blacktip
  Tot.ktch=Tot.ktch%>%
    mutate(Name=ifelse(SNAME=='australian blacktip','blacktips',Name),
           Scien.nm=ifelse(SNAME=='australian blacktip','C. limbatus & C. tilstoni',Scien.nm),
           SPECIES=ifelse(SNAME=='australian blacktip',18014,SPECIES),
           SNAME=ifelse(SNAME=='australian blacktip','blacktips',SNAME))

  
      #Bull to pigeye shark (bull not likely to be taken: Heupel & McAuley 2007 page  84)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%c(18021),"pigeye shark",SNAME),
           Name=ifelse(SPECIES%in%c(18021),'pigeye shark',Name),
           Scien.nm=ifelse(SPECIES%in%c(18021),'C. amboinensis',Scien.nm),
           SPECIES=ifelse(SPECIES%in%c(18021),18026,SPECIES))

      #Reset landed 'rays and skates' in TDGDLF to eagle ray 
  #note: fishers interviews indicated that landed rays and skate are eagle rays
  #       so adjust discards estimates to only those fishers discarding them as
  #       not all fishers retain them
  Tot.ktch=Tot.ktch%>%
    mutate(LIVEWT.c=ifelse(FishCubeCode%in%c('Discards_TDGDLF')
                           & SPECIES==39001,LIVEWT.c*prop.disc.ER,LIVEWT.c),
           Name=ifelse(FishCubeCode%in%c('JASDGDL','Historic','WCDGDL')
                       & SPECIES==31000,'eagle ray',Name),
           Scien.nm=ifelse(FishCubeCode%in%c('JASDGDL','Historic','WCDGDL')
                           & SPECIES==31000,'Myliobatis tenuicaudatus',Scien.nm),
           SNAME=ifelse(FishCubeCode%in%c('JASDGDL','Historic','WCDGDL')
                        & SPECIES==31000,'eagle ray',SNAME),
           SPECIES=ifelse(FishCubeCode%in%c('JASDGDL','Historic','WCDGDL')
                          & SPECIES==31000,39001,SPECIES))
  
  
  #Disaggregate Sawfish into species  
  Prop.by.sawfish.sp=Tot.ktch%>%
    filter(SPECIES%in%25001:25020)%>%
    group_by(SPECIES,FishCubeCode)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))%>%
    group_by(FishCubeCode)%>%
    mutate(Prop=LIVEWT.c/sum(LIVEWT.c))%>%
    dplyr::select(-LIVEWT.c)%>%
    ungroup()%>%
    data.frame
  Prop.by.sawfish.sp=rbind(Prop.by.sawfish.sp,
                           Prop.by.sawfish.sp[1:2,]%>%
                             mutate(SPECIES=25001:25002,
                                    FishCubeCode=rep('Recreational',2),
                                    Prop=c(.5,.5)))
  
  sawfish=Tot.ktch%>%
    filter(SPECIES==25000)
  
  Tot.ktch=Tot.ktch%>%
    filter(!SPECIES==25000)
  
  sawfish=sawfish%>%
    mutate(SPECIES=ifelse(FishCubeCode%in%c('JANS','OANCGC'),25001,SPECIES))  #set to green sawfish based on Survey data
  
  sawfish_JANS.OANCGC=sawfish%>%
    filter(SPECIES==25001)
  
  sawfish=sawfish%>%
    filter(!SPECIES==25001)
  
  sawfish=full_join(sawfish%>%dplyr::select(-SPECIES),    #set to proportion by species
                    Prop.by.sawfish.sp%>%
                      filter(FishCubeCode%in%sawfish$FishCubeCode),
                    by='FishCubeCode')%>%
    mutate(LIVEWT.c=LIVEWT.c*Prop)%>%
    dplyr::select(-Prop)
  
  sawfish=rbind(sawfish,
                sawfish_JANS.OANCGC)%>%
    mutate(SNAME=case_when(SPECIES==25001~'green sawfish',
                           SPECIES==25002~'narrow sawfish',
                           SPECIES==25004~'dwarf sawfish'),
           Name=case_when(SPECIES==25001~'green sawfish',
                          SPECIES==25002~'narrow sawfish',
                          SPECIES==25004~'dwarf sawfish'),
           Scien.nm=case_when(SPECIES==25001~'Pristis zijsron',
                              SPECIES==25002~'Anoxypristis cuspidata',
                              SPECIES==25004~'Pristis clavata'))
  
  Tot.ktch=rbind(Tot.ktch,sawfish)
  
    #re calculate catch after reapportioning
  Tot.ktch=Tot.ktch%>%
    group_by(across(c(-LIVEWT.c)))%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))%>%
    ungroup()
  
  #Fix Port Jackson name  
  Tot.ktch=Tot.ktch%>%
    mutate(Name=ifelse(Name=="port jackson shark","port Jackson shark",Name))
 
  
  #Final touches
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(Name=='eusphyra blochii','winghead shark',Name))
  
  Tot.ktch=Tot.ktch%>%filter(!is.na(Name))
  
  
  #Aggregate by species, year and fishery
  Tot.ktch.zone=Tot.ktch%>%
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=ifelse(SPECIES%in%c(22999),"unidentified sharks",SNAME))%>%
    group_by(SPECIES,Name,FishCubeCode,Data.set,finyear,FINYEAR,Region,zone)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))
  
  Tot.ktch.method=Tot.ktch%>%
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=ifelse(SPECIES%in%c(22999),"unidentified sharks",SNAME))%>%
    group_by(SPECIES,Name,FishCubeCode,Data.set,finyear,FINYEAR,METHOD)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))%>%
    mutate(Gear=ifelse(METHOD%in%c("BS","BH","GN","HN","Pelagic.gillnet"),"net",
                ifelse(METHOD%in%c("DL","DV","EL","GL","HL","HR","HY",
                                          "LL","Longline","Rec.line","TL"),'line',
                ifelse(METHOD%in%c("FG","TW"),'trawl',
                ifelse(METHOD%in%c("FT","PT"),'trap',
                NA)))))
  
  Tot.ktch.method=Tot.ktch.method%>%
    group_by(FishCubeCode) %>%
    arrange(FishCubeCode, is.na(Gear)) %>% # in case to keep non- NA elements for a tie
    mutate(Gear = ifelse(is.na(Gear),Mode(Gear),Gear),
           Gear=ifelse(is.na(Gear) & FishCubeCode %in% c('Indo','WTB','SA MSF'),'line',
                ifelse(is.na(Gear) & FishCubeCode %in% c('GAB','PFT','SBSC'),'trawl',
                Gear)))
  
  
  Tot.ktch=Tot.ktch%>%
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
    group_by(SPECIES,Name,FishCubeCode,Data.set,finyear,FINYEAR,Region)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))
  
  return(list(Total=Tot.ktch,Zone=Tot.ktch.zone,Total.method=Tot.ktch.method,Table1=Table1))
}
Get.ktch=fn.import.catch.data(KTCH.UNITS)    
KtCh=Get.ktch$Total
KtCh.zone=Get.ktch$Zone


#2.4. Remove species assessed elsewhere
KtCh=KtCh%>%filter(!Name%in%assessed.elsewhere)
KtCh.zone=KtCh.zone%>%filter(!Name%in%assessed.elsewhere)


#2.5 Data used for analysis of changes in reported mean weight
Wei.range=read.csv(handl_OneDrive("Data/Length_Weights/Data.Ranges.csv"),stringsAsFactors = F)
Wei.range.names=read.csv(handl_OneDrive("Data/Length_Weights/Species.names.csv"),stringsAsFactors = F)
Wei.range=merge(Wei.range,Wei.range.names,by="Sname",all.x=T)

Logbook=read.csv(handl_OneDrive("Analyses/Catch and effort/Logbook.data.mean.weight.csv"))


#---3.  Life history ------------------------------------------------------  
LH.data=read.csv(handl_OneDrive('Data/Life history parameters/Life_History.csv'))

#---4.  PSA to determine which species to assess ------------------------------------------------------  
#note: run a PSA aggregating the susceptibilities of multiple fleets (Micheli et al 2014)

if(First.run=="YES")
{
  library(yarrr)
  #get catches of all species
  KtCh.method=Get.ktch$Total.method%>%filter(!Name%in%assessed.elsewhere)
  
  #biological and fishery attributes for psa
  PSA.list=read.csv(handl_OneDrive('Analyses/Population dynamics/PSA/PSA_scores_other.species.csv'))
  PSA.list=PSA.list%>%filter(!Species%in%assessed.elsewhere)
  
  
  #Show annual catches
  Exprt=handl_OneDrive("Analyses/Population dynamics/PSA")
  KtCh.method%>%
    filter(Name%in%PSA.list$Species)%>%
    mutate(Year=as.numeric(substr(FINYEAR,1,4)),
           Gear=ifelse(is.na(Gear),'unidentified',Gear),
           Gear=capitalize(Gear),
           Name=capitalize(Name))%>%
    group_by(Name,Year,Gear)%>%
    summarise(catch=sum(LIVEWT.c))%>%
    ggplot(aes(Year,catch))+
    geom_point(aes(colour = Gear),size = .8)+ylab("Catch (tonnes)")+
    facet_wrap( ~ Name, scales = "free_y")+ expand_limits(y = 0)+
    theme_PA(strx.siz=8,leg.siz=12,axs.t.siz=6.5)+
    theme(legend.position="top",
          legend.title = element_blank(),
          legend.key=element_blank(),
          title=element_text(size=12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    guides(colour = guide_legend(override.aes = list(size=5)))+
    xlab("Financial year")
    
  ggsave(paste(Exprt,'Annual_ktch_by_species.tiff',sep='/'), width = 17,height = 7.5, dpi = 300, compression = "lzw")
  
  #Export table of species catch by fishery
  write.csv(Get.ktch$Table1,paste(Exprt,'All.species.caught.by.fishery.csv',sep='/'),row.names = F)
  
  
  #which species meet PSA criteria (and hence availability and encounterability not set to 1)?
  psa.ktch=KtCh.method%>%
    filter(Name%in%PSA.list$Species)%>%
    group_by(Name,FINYEAR)%>%
    summarise(catch=sum(LIVEWT.c))%>%
    spread(Name,catch,fill=0)
  aa=as.matrix(psa.ktch[,-1])
  psa.species.max.ever=names(which(apply(aa,2,max)>=PSA.max.ton))
  aa[aa<=PSA.min.tons]=0
  aa[aa>PSA.min.tons]=1
  psa.species.series=names(which(colSums(aa)>PSA.min.years))
  species.meeting.criteria=intersect(psa.species.series,psa.species.max.ever)
  
  
  #replace missing gear info with most common value
  Agg=KtCh.method%>%
    group_by(FishCubeCode) %>%
    arrange(FishCubeCode, is.na(Gear)) %>% # in case to keep non- NA elements for a tie
    mutate(Gear = ifelse(is.na(Gear),Mode(Gear),Gear),
           Gear=ifelse(is.na(Gear) & FishCubeCode %in% c('Historic','Indo','WTB','SA MSF'),'line',
                ifelse(is.na(Gear) & FishCubeCode %in% c('GAB','PFT','SBSC'),'trawl',
                Gear)))
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
      mutate(Name=d$Name,
             Sum.yrs=rowSums(d1),
             Keep=ifelse(Sum.yrs>=Min.yrs,"YES","NO"))
    KIP[[s]]=kip
  }
  KIP=do.call(rbind,KIP)%>%
    filter(Keep=="YES")%>%
    dplyr::select(Name,Gear)
  
  #run PSA
  UniSp=unique(KtCh$Name)
  UniSp=subset(UniSp,!UniSp%in%names(Indicator.species)) 
  
  PSA.list=PSA.list%>%filter(Species%in%UniSp)  
  
  PSA.fn=function(d,line.sep,size.low,size.med,size.hig,W,H)  
  {
    set.seed(101)
    PSA=data.frame(Species=d$Species,
                   Productivity=rep(NA,nrow(d)),
                   Susceptibility=rep(NA,nrow(d)),
                   Vulnerability=rep(NA,nrow(d)))
    for(p in 1:nrow(d))    
    {
      aa=d[p,]
      if(!aa$Species%in%species.meeting.criteria) #reset availability and encounterability if not meeting criteria
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
      rename(Vulnerability.score=Vulnerability)%>%
      mutate(Species=Hmisc::capitalize(as.character(Species)),
             Vulnerability=Hmisc::capitalize(
                      ifelse(Vulnerability.score<=Low.risk,'Low',
                      ifelse(Vulnerability.score>Low.risk & Vulnerability.score<=medium.risk,'Medium',
                      'High'))),
             Vulnerability=factor(Vulnerability,levels=c('Low','Medium','High')))%>%
      arrange(Vulnerability.score)%>%
      mutate(Fnt.size=case_when(Vulnerability=='Low'~size.low,
                                Vulnerability=='Medium'~size.med,
                                Vulnerability=='High'~size.hig),
             Species=case_when(Species=='Port jackson shark'~'Port Jackson shark',
                               TRUE~Species))
    cols=c(Low="green",Medium="yellow",High="red")
    p=ggplot(PSA,
             aes(Productivity, Susceptibility, label = Species,colour = Vulnerability, fill = Vulnerability)) +
      geom_point(shape = 21, size = 5,colour="black") + 
      geom_text_repel(aes(size=Fnt.size),show.legend  = F,segment.colour=transparent("black",.75),col='black',box.padding = line.sep) + 
      scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))+ 
      xlim(1.5,3.15)+ylim(0.85,3.05)+
      theme_PA(axs.T.siz=18,axs.t.siz=12,lgT.siz=16,leg.siz=15)+
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position="top",
            legend.key=element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))+
      labs(fill = "")
    p
    ggsave(paste(Exprt,'Figure 2_PSA.tiff',sep='/'), width = W,height = H, dpi = 300, compression = "lzw")
    
    #Table PSA scores
    Table.PSA=d
    for(p in 1:nrow(Table.PSA))    
    {
      if(!Table.PSA[p,]$Species%in%species.meeting.criteria) #reset availability and encounterability if not meeting criteria
      {
        k=KIP%>%filter(Name==Table.PSA[p,]$Species)
        Table.PSA[p,]=Table.PSA[p,]%>%
          mutate(Net.avail=ifelse(!'net'%in%k$Gear,1,Net.avail),
                 Net.encoun=ifelse(!'net'%in%k$Gear,1,Net.encoun),
                 Line.avail=ifelse(!'line'%in%k$Gear,1,Line.avail),
                 Line.encoun=ifelse(!'line'%in%k$Gear,1,Line.encoun),
                 Trawl.avail=ifelse(!'trawl'%in%k$Gear,1,Trawl.avail),
                 Trawl.encoun=ifelse(!'trawl'%in%k$Gear,1,Trawl.encoun),
                 Trap.avail=ifelse(!'trap'%in%k$Gear,1,Trap.avail),
                 Trap.encoun=ifelse(!'trap'%in%k$Gear,1,Trap.encoun))
      }
    }
    Table.PSA=Table.PSA%>%
      mutate(Species=capitalize(Species))%>%
      left_join(PSA%>%
                  dplyr::select(Species,Vulnerability.score),
                by='Species')%>%
      mutate(Vulnerability.score=round(Vulnerability.score,2))
    write.csv(Table.PSA,paste(Exprt,'Table S2_PSA scores.csv',sep='/'),row.names=F)
    return(as.character(PSA%>%filter(Vulnerability=="High")%>%pull(Species)))
  }
  Keep.species=PSA.fn(d=PSA.list,line.sep=.35,size.low=2.1,size.med=2.15,size.hig=2.5,W=10,H=10)
  Keep.species=tolower(Keep.species)
  Keep.species=sort(c(Keep.species,names(Indicator.species)))
  Drop.species=UniSp[which(!UniSp%in%Keep.species)]
}
if(!First.run=="YES")
{
  Keep.species=sort(c( "angel sharks",  "copper shark",  "great hammerhead",    
                       "grey nurse shark", "lemon shark",  "pigeye shark",        
                       "sawsharks",    "scalloped hammerhead", "shortfin mako",       
                       "smooth hammerhead","spinner shark","spurdogs", "tiger shark",         
                       "wobbegongs",
                       names(Indicator.species)))
}
  
N.sp=length(Keep.species)

#---5.   Apply CMSY, DB-SRA, OCOM, JABBA-------------------------------------------------------
apply.CMSY=function(year,catch,r.range,Bo.low,Bo.hi,Int.yr,Bint.low,Bint.hi,Bf.low,Bf.hi)
{
  inputs= list(r.low=r.range[1],
               r.hi=r.range[2],
               stb.low=Bo.low,
               stb.hi=Bo.hi,
               int.yr=Int.yr,
               intb.low = Bint.low,
               intb.hi = Bint.hi,
               endb.low=Bf.low,
               endb.hi=Bf.hi)
  
  output <- cmsy2(year=year,
                  catch=catch,
                  r.low=r.range[1],
                  r.hi=r.range[2],
                  stb.low=Bo.low,
                  stb.hi=Bo.hi,
                  int.yr=Int.yr,
                  intb.low = Bint.low,
                  intb.hi = Bint.hi,
                  endb.low=Bf.low,
                  endb.hi=Bf.hi)
  
  # Extract reference points and time series from output
  output$ref_pts <- output[["ref_pts"]]
  output$ref_ts <- output[["ref_ts"]]
  
  return(list(inputs=inputs,output=output))
}

apply.DBSRA=function(year,catch,catchCV,catargs,agemat,k,b1k,btk,fmsym,bmsyk,M,graph,nsims,grout,WD)
{
  setwd(WD)  #dbsra automatically exports the biomass trajectories
  #store inputs
  inputs=list(year=year,
              catch=catch,
              catchCV=catchCV,
              catargs=catargs,
              agemat=agemat,
              k=k,
              b1k=b1k,
              btk=btk,
              fmsym=fmsym,
              bmsyk=bmsyk,
              M=M,
              graph=graph,
              nsims=nsims,
              grout=grout)
  
  #run model
  outputs <- dbsra(year=year,
                   catch=catch,
                   catchCV=catchCV,
                   catargs=catargs,
                   agemat=agemat,
                   k=k,
                   b1k=b1k,
                   btk=btk,
                   fmsym=fmsym,
                   bmsyk=bmsyk,
                   M=M,
                   graph=graph,
                   nsims=nsims,
                   grout=grout)
  
  return(list(inputs=inputs,outputs=outputs))
}

apply.OCOM=function(year,catch,M)
{
  inputs= list(m=M)
  output <- ocom(year=year,catch=catch,m=M)
  
  # Extract reference points and time series from output
  output$ref_pts <- output[["ref_pts"]]
  output$ref_ts <- output[["ref_ts"]]
  
  return(list(inputs=inputs,output=output))
}

apply.JABBA=function(Ktch,CPUE=NULL,CPUE.SE=NULL,MDL,ASS,Rdist,Rprior,Kdist,Kprior,PsiDist,Psiprior,Bprior,output.dir)
{
  # Compile JABBA JAGS model
  if(is.null(CPUE))
  {
    jbinput = build_jabba(catch=Ktch,
                          model.type = MDL,
                          assessment=ASS,
                          scenario =  "CatchOnly",
                          r.dist = Rdist,
                          r.prior = Rprior,
                          K.dist= Kdist,
                          K.prior=Kprior,
                          psi.dist=PsiDist,
                          psi.prior=Psiprior,
                          b.prior= Bprior)
  }else
  {
    jbinput = build_jabba(catch=Ktch,
                          cpue=CPUE,
                          se=CPUE.SE,
                          model.type = MDL,
                          assessment=ASS,
                          scenario =  "CatchOnly",
                          r.dist = Rdist,
                          r.prior = Rprior,
                          K.dist= Kdist,
                          K.prior=Kprior,
                          psi.dist=PsiDist,
                          psi.prior=Psiprior,
                          b.prior= Bprior)
  }

  # Fit JABBA
  return(fit_jabba(jbinput,save.csvs=TRUE,output.dir=output.dir))
}

#---6.   Bring in demography and steepness functions-------------------------------------------------------
if(do.r.prior)
{
  fn.source("Leslie.matrix.R") 
  fun.rprior.dist=function(Nsims,K,LINF,K.sd,LINF.sd,k.Linf.cor,Amax,MAT,FecunditY,
                           Cycle,BWT,AWT,LO)
  {
    Fecu=unlist(FecunditY)
    Rprior=fun.Leslie(N.sims=Nsims,k=K,Linf=LINF,k.sd=K.sd,Linf.sd=LINF.sd,k.Linf.cor=k.Linf.cor,
                      A=Amax,first.age=0,RangeMat=MAT,Rangefec=Fecu,
                      sexratio=0.5,Reprod_cycle=Cycle,
                      bwt=BWT,awt=AWT,Lo=LO,
                      Resamp=RESAMP)  
    
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
}
if(do.steepness)
{
  fn.source("Steepness.R")
}
#---7.  Sawfish stand-alone assessment-----
#notes: Run whatever assessment is appropriate see Pillans et al 2021.
# Do spatial by year presence/absence Pilbara trawl to see shrinkage
if(do.sawfish)
{
  hNdl.sawfish=handl_OneDrive(paste('Analyses/Population dynamics/Other species/Sawfishes/',Year.of.assessment,sep=''))
  if(!dir.exists(hNdl.sawfish))dir.create(hNdl.sawfish)
  
  #1. Get catch
  Sawfish.ktch=Get.ktch$Total.method%>%
    filter(SPECIES%in%25000:25020)%>%
    filter(finyear>=1960)%>%
    mutate(Name=capitalize(Name))
  
  xx=Sawfish.ktch%>%
                    data.frame%>%
                    distinct(Name,SPECIES)
  assessed.sawfish=xx%>%pull(Name)
  names(assessed.sawfish)=xx%>%pull(SPECIES)
  n.sawfish=length(assessed.sawfish)
  
  Sawfish.ktch%>%
    filter(LIVEWT.c>0)%>%
    mutate(Fishery.gear=paste(Gear,FishCubeCode,sep='-'))%>%
    ggplot(aes(finyear,LIVEWT.c,colour=Fishery.gear))+
    geom_point()+geom_line(alpha=0.25)+
    facet_wrap(~Name,nrow=3,scales='free')+
    theme_PA(strx.siz=14,leg.siz=11,axs.t.siz=13,axs.T.siz=16)+
    ylab("Catch (tonnes)")+xlab("Financial year")+
    theme(legend.position="top",
          legend.title = element_blank(),
          legend.key=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=5,linetype = 0)))
  ggsave(paste(hNdl.sawfish,'Annual_ktch_by_species.tiff',sep='/'), width = 10,height = 10, dpi = 300, compression = "lzw")
  
  #combined catches (in tonnes)
  Sawfish.ktch.combined=Sawfish.ktch%>%
                group_by(SPECIES,Name,finyear)%>%
                summarise(Tonnes=sum(LIVEWT.c,na.rm=T))
  
  #2. Remove freshwater sawfish due to extremely low catches, unaccounted customary catch, model convergence, etc
  assessed.sawfish=subset(assessed.sawfish,!assessed.sawfish=="Freshwater sawfish")
  n.sawfish=length(assessed.sawfish)
  
  #3. Get Pilbara Trawl cpue
  hNdl.sawfish.PFT.cpue=handl_OneDrive('Analyses/Data_outs')
  Sawfish.cpue.list=vector('list',n.sawfish)
  names(Sawfish.cpue.list)=assessed.sawfish
  for(n in 1:n.sawfish)
  {
    x=paste(hNdl.sawfish.PFT.cpue,names(Sawfish.cpue.list)[n],sep='/')
    if(dir.exists(x))
    {
      dumy=paste(x,'CPUE_Pilbara.trawl.csv',sep='/')
      if(file.exists(dumy))
      {
        rel.cpue=read.csv(dumy)
        rel.cpue=rel.cpue%>%
          mutate(LOW=LOW/mean(MEAN,na.rm=T),
                 UP=UP/mean(MEAN,na.rm=T),
                 MEAN=MEAN/mean(MEAN,na.rm=T))
        Sawfish.cpue.list[[n]]=rel.cpue
      }
    }
  }
  
  #4. Get Life history
  Sawfish.life.history=vector('list',n.sawfish)
  names(Sawfish.life.history)=assessed.sawfish
  for(n in 1:n.sawfish)
  {
    DD=LH.data%>%filter(SPECIES==names(assessed.sawfish)[n])
    if(is.na(DD$k.sd))DD$k.sd=DD$K*.2
    if(is.na(DD$FL_inf.sd))DD$FL_inf.sd=DD$FL_inf*.2
    if(is.na(DD$Max_Age_max))DD$Max_Age_max=round(DD$Max_Age*1.3)
    if(is.na(DD$Fecu_max))DD$Fecu_max=DD$Fecu_min
    if(is.na(DD$Cycle_max))DD$Cycle_max=DD$Cycle
    Sawfish.life.history[[n]]=DD
  }
  
  #5. Get population growth rate thru demography
  store.sawfish.r=vector('list',n.sawfish)
  names(store.sawfish.r)=assessed.sawfish
  store.sawfish.M=store.sawfish.r
  if(do.r.prior)
  {
    system.time({for(l in 1:n.sawfish)   #takes 0.013 sec per iteration per species
    {
      M.averaging<<-"min" #'min' yields rmax, Cortes pers com, but yields too high steepness for all species
      RESAMP="YES"
      linear.fec="NO"
      print(paste("r prior ","--",names(Sawfish.life.history)[l]))
      r.prior.dist=with(Sawfish.life.history[[l]],fun.rprior.dist(Nsims=1000,K=K,LINF=FL_inf,
                                                                  K.sd=k.sd,LINF.sd=FL_inf.sd,k.Linf.cor=-0.99,
                                                                  Amax=c(Max_Age,Max_Age_max),
                                                                  MAT=c(Age_50_Mat_min,Age_50_Mat_max),
                                                                  FecunditY=c(Fecu_min,Fecu_max),
                                                                  Cycle=c(Cycle,Cycle_max),
                                                                  BWT=b_w8t,AWT=a_w8t,
                                                                  LO=LF_o))
      
      #export r and M
      hndl1=paste(hNdl.sawfish,names(store.sawfish.r)[l],sep='/')
      if(!dir.exists(hndl1))dir.create(hndl1)
      out.r=data.frame(shape=r.prior.dist$shape,rate=r.prior.dist$rate,
                       mean=r.prior.dist$mean,sd=r.prior.dist$sd)
      write.csv(out.r,paste(hndl1,'r.prior.csv',sep='/'),row.names = F)
      store.sawfish.r[[l]]=r.prior.dist
      
      n.dim=max(unlist(lapply(r.prior.dist$M,length)))
      out.M=r.prior.dist$M
      for(ss in 1:length(out.M))
      {
        a=out.M[[ss]]
        delta=n.dim-length(a)
        if(delta>0) out.M[[ss]]=c(a,rep(NA,delta))
        rm(a)
      }
      out.M=do.call(rbind,out.M)
      names(out.M)=0:(n.dim-1)
      write.csv(out.M,paste(hndl1,'M.csv',sep='/'),row.names=FALSE)
      store.sawfish.M[[l]]=out.M
      rm(r.prior.dist,M.averaging,RESAMP,linear.fec)
    }})
  }else
  {
    for(l in 1:n.sawfish) 
    {
      hndl1=paste(hNdl.sawfish,names(store.sawfish.r)[l],sep='/')
      store.sawfish.r[[l]]=read.csv(paste(hndl1,"r.prior.csv",sep='/'))
      store.sawfish.M[[l]]=read.csv(paste(hndl1,"M.csv",sep='/'))
    }
  }
  
  #6. Get Maximum lifetime reproductive rate (alpha)
  fn.source("Steepness.R")
  store.sawfish.alpha=store.sawfish.r
  for(l in 1:n.sawfish) 
  {
    store.sawfish.alpha[[l]]=with(Sawfish.life.history[[l]],
                                  Alpha.Brooks(max.age=mean(c(Max_Age,Max_Age_max),na.rm=T),
                                          M=apply(store.sawfish.M[[l]],2,mean,na.rm=T),
                                          age.mat=mean(c(Age_50_Mat_min,Age_50_Mat_max),na.rm=T),
                                          Meanfec=mean(c(Fecu_min,Fecu_max),na.rm=T),
                                          CyclE=mean(c(Cycle,Cycle_max),na.rm=T)))
  }
  
  #7. Get scaler for Fmsy.to.M relationship (Fmsy= scaler x M)
  Fmsy.M.scaler=store.sawfish.alpha
  for(l in 1:n.sawfish)
  {
    Fmsy.M.scaler[[l]]=Cortes.Brooks.2018(alpha=store.sawfish.alpha[[l]])
  }
  
  
  #ACA set up scenarios for each model
  #8. Run DBSRA assessment (Dick and MAcCall (2011))
  # summary of method: http://toolbox.frdc.com.au/wp-content/uploads/sites/19/2020/07/DBSRA3.html
  DBSRA.sawfish=store.sawfish.r
  system.time({for(i in 1:length(DBSRA.sawfish))  #takes 0.3 secs per iteration
  {
    print(paste("DBSRA ","--",names(DBSRA.sawfish)[i]))
    this.wd=paste(hNdl.sawfish,names(DBSRA.sawfish)[i],'DBSRA',sep='/')
    if(!dir.exists(this.wd))dir.create(this.wd)
    ktch=Sawfish.ktch.combined%>%
      filter(Name==names(DBSRA.sawfish)[i])
    Run=apply.DBSRA(year=ktch$finyear,
                    catch=ktch$Tonnes,
                    catchCV=NULL,  #catch CV not available
                    catargs=list(dist="none",low=0,up=Inf,unit="MT"),  #not catch CV available
                    agemat=Sawfish.life.history[[i]]$Age_50_Mat_min,
                    k=list(low=max(ktch$Tonnes),up=max(c(max(ktch$Tonnes)*100,500)),
                          tol=0.01,permax=1000),
                    b1k=list(dist="unif",low=0.8,up=0.99,mean=1,sd=0.1),
                    btk=list(dist="unif",low=0.15,up=0.95,mean=1,sd=0.1,
                             refyr=max(ktch$finyear)),  #reference year
                    fmsym=list(dist="lnorm",low=0.1,up=2,
                               mean=log(Fmsy.M.scaler[[i]]),sd=0.2), # Cortes & Brooks 2018
                    bmsyk=list(dist="beta",low=0.05,up=0.95,mean=0.5,sd=0.1),
                    M=list(dist="lnorm",low=0.001,up=1,
                           mean=log(mean(apply(store.sawfish.M[[i]],2,mean,na.rm=T))),
                                          sd=mean(apply(store.sawfish.M[[i]],2,sd,na.rm=T))),
                    graph=c(1:14),
                    nsims=1e3, #MISSING, replace by 1e4
                    grout=1,
                    WD=this.wd)
    legend('topright',names(DBSRA.sawfish)[i],bty='n')
    Run$outputs$Biom.traj=read.csv(paste(this.wd,"Biotraj-dbsra.csv",sep='/'),header=FALSE)%>%
                            filter(V1==1)%>%  #select only possible runs
                            dplyr::select(-V1)

    DBSRA.sawfish[[i]]=Run
  }})
  
  #9. Run CMSY (Froese et al 2017)    #does not converge for dwarf or freshwater
  #summary of method: http://toolbox.frdc.com.au/wp-content/uploads/sites/19/2021/04/CMSY.html
  CMSY.sawfish=store.sawfish.r
  system.time({for(i in 1:length(CMSY.sawfish))  #takes 150 secs per species (1e4 iterations)
  {
    print(paste("CMSY ","--",names(CMSY.sawfish)[i]))
    
    ktch=Sawfish.ktch.combined%>%
      filter(Name==names(CMSY.sawfish)[i])
    r.range=quantile(rnorm(1e3,mean=store.sawfish.r[[i]]$mean,sd=store.sawfish.r[[i]]$sd),
                     probs=c(0.025,0.975))
    year=ktch$finyear
    catch=ktch$Tonnes
    Int.yr=round(mean(ktch$finyear))
    CMSY.sawfish[[i]]=apply.CMSY(year,catch,r.range,Bo.low=0.8,Bo.hi=0.99,
                                 Int.yr,Bint.low=0.15,Bint.hi=0.95,
                                 Bf.low=0.15,Bf.hi=0.95)
    
  }})
  
  #10. Run OCOM assessment (Zhou et al (2018))
  OCOM.sawfish=store.sawfish.r
  system.time({for(i in 1:length(OCOM.sawfish))  #takes 20 secs per species (1e4 iterations)
  {
    print(paste("OCOM ","--",names(OCOM.sawfish)[i]))
    
    ktch=Sawfish.ktch.combined%>%
      filter(Name==names(OCOM.sawfish)[i])
    OCOM.sawfish[[i]] <-apply.OCOM(year=ktch$finyear,
                                   catch=ktch$Tonnes,
                                   M=mean(apply(store.sawfish.M[[i]],2,mean,na.rm=T)))
    
  }})
  
  #11. JABBA - catch only (Winker et al 2018)
  #summary of method: https://github.com/jabbamodel/JABBA
  JABBA.sawfish=store.sawfish.r
  system.time({for(i in 1:length(JABBA.sawfish))  #takes XX secs per species
  {
    print(paste("JABBA ","--",names(JABBA.sawfish)[i]))
    this.wd=paste(hNdl.sawfish,names(JABBA.sawfish)[i],'JABBA',sep='/')
    if(!dir.exists(this.wd))dir.create(this.wd)
    
    ktch=Sawfish.ktch.combined%>%
      filter(Name==names(JABBA.sawfish)[i])%>%
      rename(Year=finyear,
             Total=Tonnes)%>%
      ungroup()%>%
      dplyr::select(Year,Total)%>%
      arrange(Year)%>%
      data.frame
    Bint=runif(1000,0.8,0.99)
    Bint.mean=mean(Bint)
    Bint.CV=sd(Bint)/Bint.mean
    
    Bfin=runif(1000,0.15,0.95)
    Bfin.mean=mean(Bfin)
    Bfin.CV=sd(Bfin)/Bfin.mean
    
    input=list(Ktch=ktch,
               MDL="Fox",
               ASS=names(JABBA.sawfish)[i],
               Rdist = "lnorm",
               Rprior = c(store.sawfish.r[[i]]$mean,store.sawfish.r[[i]]$sd),
               Kdist="range",
               Kprior=c(max(ktch$Total),max(c(max(ktch$Total)*100,500))),
               PsiDist='lnorm',
               Psiprior=c(Bint.mean,Bint.CV),
               Bprior=c(Bfin.mean,Bfin.CV,max(ktch$Year),"bk"))
    
    output=apply.JABBA(Ktch=ktch,
                       MDL="Fox",
                       ASS=names(JABBA.sawfish)[i],
                       Rdist = "lnorm",
                       Rprior = c(store.sawfish.r[[i]]$mean,store.sawfish.r[[i]]$sd),
                       Kdist="range",
                       Kprior=c(max(ktch$Total),max(c(max(ktch$Total)*100,500))),
                       PsiDist='lnorm',
                       Psiprior=c(Bint.mean,Bint.CV),
                       Bprior=c(Bfin.mean,Bfin.CV,max(ktch$Year),"bk"),
                       output.dir=this.wd)
    
    JABBA.sawfish[[i]]=list(input=unput,output=output)
  }})
  
  ##Move this to outputs evaluation
  #JABBA outputs
  jbplot_catch(output)
  jbplot_ppdist(output)
  jbplot_mcmc(output)
  jbplot_procdev(output)
  jbplot_bprior(output)
  jbplot_trj(output,type="BBmsy",add=T)
  jbplot_trj(output,type="FFmsy",add=T)
  
  # status summary
  par(mfrow=c(3,2),mar = c(3.5, 3.5, 0.5, 0.1))
  jbplot_trj(output,type="B",add=T)
  jbplot_trj(output,type="F",add=T)
  jbplot_trj(output,type="BBmsy",add=T)
  jbplot_trj(output,type="FFmsy",add=T)
  jbplot_spphase(output,add=T)
  jbplot_kobe(output,add=T)

}

#---8.  Create list of species assessed and import species-specific data-----
#note: this brings in any info on cpue, abundance, selectivity, size composition, tagging
Species.data=vector('list',length=N.sp)
names(Species.data)=Keep.species
for(s in 1:N.sp) 
{
  this.one=capitalize(Keep.species[s])
  if(this.one=="Angel sharks") this.one="Australian angelshark"
  files <- list.files(path = paste(Dat.repository,this.one,sep=""), pattern = "*.csv", full.names = T)
  file.names <- list.files(path = paste(Dat.repository,this.one,sep=""), pattern = "*.csv", full.names = F)
  removE <- c(".csv", paste(this.one,"_",sep=''), this.one)
  file.names <- gsub("^\\.","",str_remove_all(file.names, paste(removE, collapse = "|")))
  files=sapply(files, fread, data.table=FALSE)
  names(files)=file.names
  Species.data[[s]]=files
  rm(files)

}


#---9.  Import input parameters, define modeling arguments and create pin file-----
#note: For integrated model, S1 and S2 calculates pars in normal space but same order magnitude
#       Other scenarios all pars in log.
#       ln_RZERO is in 1,000 individuals so do 10 times the largest catch divided by 
#       average weight and divided by 1,000. Best units to work in are 1,000 individuals for 
#       numbers, catch in tonnes and length-weight in kg as all cancels out and predicted biomasses
#       end up being in tonnes

List.sp=vector('list',N.sp)
names(List.sp)=Keep.species
for(l in 1:N.sp)
{
  this=subset(All.species.names,SNAME==Keep.species[l])%>%
    mutate(SP=ifelse(SNAME=="sawsharks","SC",
                     ifelse(SNAME=="spurdogs","SD",
                            SP)),
           name.inputs=ifelse(SNAME=="sawsharks","common sawshark",
                              ifelse(SNAME=="spurdogs","spikey dogfish",
                                     SNAME)))
  fst.yr=ifelse(this$SPECIES%in%Indicator.species,"1975-76",min(KtCh$FINYEAR))
  List.sp[[l]]=list(Name=this$SNAME,
                    Name.inputs=this$name.inputs,
                    SP=this$SP,
                    Species=this$SPECIES,
                    First.year=fst.yr)
}

fn.mtch=function(WHAT,NMS) match(WHAT,names(NMS))
Q_phz=c("lnq","lnq2","log_Qdaily")                           
Zns.par.phz=c("lnR_prop_west","lnR_prop_zn1")
MOv.par.phz=c("log_p11","log_p22","log_p21","log_p33")
for(l in 1:N.sp)
{
  print(paste("---------",names(List.sp)[l]))
  
  LH=LH.data%>%filter(SPECIES==List.sp[[l]]$Species)
  
  #... Demography arguments
  Max_Age_max=LH$Max_Age_max
  if(is.na(Max_Age_max))Max_Age_max=round(LH$Max_Age*1.3)
  List.sp[[l]]=list.append(List.sp[[l]],
  pup.sx.ratio=0.5,
  Growth.F=data.frame(k=LH$K,FL_inf=LH$FL_inf,k.sd=LH$k.sd,FL_inf.sd=LH$FL_inf.sd),
  k.Linf.cor=-0.99,    #assumed correlation between growth parameters
  Max.age.F=c(LH$Max_Age,Max_Age_max),
  Age.50.mat=c(LH$Age_50_Mat_min,LH$Age_50_Mat_max),
  Fecundity=c(LH$Fecu_min,LH$Fecu_max),
  Breed.cycle=c(LH$Cycle,LH$Cycle_max),
  TEMP=LH$Temperature,
  BwT=LH$b_w8t,
  AwT=LH$a_w8t,
  TLmax=LH$Max.TL,
  Lzero=LH$LF_o,
  NsimSS=1000,
  r.prior="USER",                    #demography
  r.prior2=NA,                       #uniform
  a_FL.to.TL=LH$a_FL.to.TL,          # FL to TL
  b_FL.to.TL=LH$b_FL.to.TL)                

  #indicator species arguments
  if(List.sp[[l]]$Species%in%Indicator.species)
  {
    #1. Add input parameters
    List.sp[[l]]=list.append(List.sp[[l]],
                Data.yr=Last.yr.ktch,                          #last year of catch
                Frst.yr.ktch=List.sp[[l]]$First.year,          #first year of catch
                BaseCase=basecase,
                Do.cols=do.cols, 
                Max.FL.obs=LH$Max.FL.obs,                      #Maximimum observed FL
                Lo=LH$LF_o,                                    #Size at birth
                #-----Catch-MSY arguments
                SIMS=1e5,                   #simulations
                Proc.err=0,                 #sigR is PROCESS ERROR; 0 if deterministic model
                Growth.F=data.frame(k=LH$K,FL_inf=LH$FL_inf),
                TEMP=LH$Temperature,
                Max.age.F=LH$Max_Age,
                Age.50.mat=c(LH$Age_50_Mat_min,LH$Age_50_Mat_max),
                Fecundity=c(LH$Fecu_min,LH$Fecu_max),
                Breed.cycle=LH$Cycle,  #years
                years.futures=5)
                
    if(List.sp[[l]]$Name=="whiskery shark")
    {
        n.scen=13
        Drop_yr_cpue=c("1980-81","1981-82","1982-83","1983-84")  #Dropped cpue years
        Drop_yr_cpue_sens=c("1975-76","1976-77","1977-78","1978-79","1979-80",
                            "1980-81","1981-82","1982-83")
        Drop_yr_cpue.tabl=paste(substr(Drop_yr_cpue,1,4)[1],
                                substr(Drop_yr_cpue,3,4)[length(Drop_yr_cpue)],sep="-")
        Drop_yr_cpue_sens.tabl=paste(substr(Drop_yr_cpue_sens,1,4)[1],
                                     substr(Drop_yr_cpue_sens,3,4)[length(Drop_yr_cpue_sens)],sep="-")
                  
        List.sp[[l]]=list.append(List.sp[[l]],
              Prior.mean.Fo=0.01,
              Prior.SD.Log.Fo=0.5,
                  
              #Steepness
              h.M.constant=0.351,       #Braccini et al 2015 
              h.M.constant.low=0.29,     #80% percentile from Braccini et al 2015   
              h.M.constant.up=0.41, 
                  
              #Natural mortality
              M_val=0.27,          #using Hoenig set to Max age=15  (Simpfendorfer et al 2000)
              M_val.low=0.23,      #using Hoenig set to Max age=18 exp(1.46-1.01*log(18))
              M_val.high=0.35,     #using Hoenig set to Max age=12 exp(1.46-1.01*log(12))
                  
              #Initial F
              Fo=0.03,             #yields a B1975 of 90% virgin           
              Fo_Simp=0.003,       #Simpfendorfer et al 2000 estimated at 0.003
              Fo_M=0.05,            #yields a B1975 of 85% virgin 
                  
              Po_spm=0.9,  #Po for surplus production, consistent with the Fo value used in Size based model
                  
              #Data
              AREAS=c("West","Zone1","Zone2"),  #Define spatial areas; 1 is West, 2 is Zn1, 3 is Zn2.
              Yr_q_change=1982,   #last year before targeting practices changed (Simpfendorfer 2000)
              Yr_q_daily=2006,
              Do_var=0,    #How to calculate cpue variance in age-structured
              Var1=0.0296133485565278,   #fixed variances used by Simpfendorfer et al 2000
              Var2=0.023525088856464,

             #Initial value of estimable parameters
             #note: Dummy used for switching on/off phase estimation in simulation testing
             Dummy=1,   
             #----Biomass dynamics
             r_BD=fn.ji(0.2),
             k_BD=fn.ji(5000),
             Q1_BD=fn.ji(7e-4),
             Q2_BD=fn.ji(2e-4),
             daily_BD=fn.ji(4e-4),
             tau2_BD=fn.ji(0.1353),
             #----Age-structured
             RSTAR_AS=fn.ji(55),
             z_AS=fn.ji(.6),
             Q1_AS=fn.ji(.4),
             Q2_AS=fn.ji(.2),
             q_daily_AS=fn.ji(0.8),
             Fo_AS=fn.ji(0.008),
             #----Size-base
             RZERO_in_1000_individuals_SB=fn.ji(1096),
             Q1_SB=fn.ji(0.00041),
             Q2_SB=fn.ji(0.0001),
             Q_daily_SB=fn.ji(0.00024),
             Fo_SB=NA,  #no jit because it's fixed
             tau_SB=fn.ji(0.2), 
             K.F=fn.ji(0.38),
             Linf.F=fn.ji(130),
             K.M=fn.ji(0.423),
             Linf.M=fn.ji(130),
             SD.growth_SB=fn.ji(7),
             Prop.west=fn.ji(0.19),  #proportion of recruitment by zone calculated based on average prop catch by zone
             Prop.zn1=fn.ji(0.38),
             #---movement 
             p11=0.99,
             p22=0.99,
             p21=0.001,
             p33=0.99,
                  
             #Estimation phases
             Par.phases=list('Base case'=c(dummy=Fz.off,lnR_zero=1,lnR_prop_west=1,lnR_prop_zn1=1,
                                           lnq=2,lnq2=2,log_Qdaily=2,ln_Init_F=Fz.off,log_tau=5,
                                           k=3,lnLinf=3,k_M=3,lnLinf_M=3,sd_growth=4,
                                           log_p11=1,log_p22=1,log_p21=1,log_p33=1),
                              S1=c(dummy=Fz.off,r=2,ln_k=2,ln_q1=1,ln_q2=1,ln_qdaily=-1,ln_tau=3),
                              S2=c(dummy=Fz.off,Rstar=2,z=3,q1=1,q2=1,qdaily=Fz.off,Fo=4))
             )
        
        #Scenarios
          #--Integrated model
        List.sp[[l]]$N.Scens=n.scen       
        Zens=paste("S",1:(n.scen-1),sep="")
        Models=c("Base case",paste("S",1:(n.scen-1),sep=""))
        Q.scen=c(rep("three",2),"two","three",rep("two",1),rep("three",8))
        List.sp[[l]]=with(List.sp[[l]],
        {list.append(List.sp[[l]],
                    #--Integrated model 
                    Tabla.scen=data.frame(
                      Model=Models,
                      Size_comp.=c('Yes',"N/A",'No',rep('Yes',10)),
                      CPUE=rep("Stand.",13),
                      CPUE_years_dropped=c(rep(Drop_yr_cpue.tabl,2),rep("None",2),
                                           Drop_yr_cpue_sens.tabl,
                                           rep(Drop_yr_cpue.tabl,8)),
                      Age.Growth=c('Yes',"N/A",'No',rep('Yes',10)),
                      Ktch.sx.r=c('Observed','N/A','Equal',rep('Observed',2),'Equal',rep('Observed',7)),                      
                      Tagging=c('No','N/A',rep('No',10),'Yes'),                     
                      Fec.=c(rep('N/A',2),rep('constant',1),rep('N/A',10)),
                      Maturity=c('at length','N/A',rep('knife edge',1),rep("at length",10)),
                      M=c("constant","N/A",rep("constant",11)),
                      M.value=c(M_val,NA,rep(M_val,4),M_val.low,M_val.high,rep(M_val,5)),                      
                      SteepnesS=c(h.M.constant,rep("N/A",2),rep(h.M.constant,7),
                                 h.M.constant.low,h.M.constant.up,h.M.constant),
                      Q=Q.scen,   
                      Spatial_structure=c(rep('Single zone',12),'Three zones'),
                      Movement=c("No",rep("N/A",2),rep("No",9),"Yes"),
                      Fo=c(Fo,"N/A","estimated",rep(Fo,5),Fo_Simp,Fo_M,rep(Fo,3)),
                      Model_type=c('Length-based','Biomass dynamics',"Age-structured",rep("Length-based",10))
                    ),
                    #--CMSY
                    ktch_msy_scen=list(
                      'Base case'=list(r.prior="USER",user="Yes",k.max=50,startbio=c(0.7,.95),
                                       finalbio=c(0.2, 0.7),res="low",
                                       niter=SIMS,sigR=Proc.err),
                      S1=list(r.prior="USER",user="Yes",k.max=50,startbio=c(0.7,.95),
                              finalbio=c(0.2, 0.7),res="low",
                              niter=SIMS,sigR=0.02),
                      S2=list(r.prior=NA,user="No",k.max=50,startbio=c(0.7,.95),
                              finalbio=c(0.2, 0.7),res="low",
                              niter=SIMS,sigR=Proc.err))
                   )
                  })  
      }
    
    if(List.sp[[l]]$Name=="gummy shark")
    {
       n.scen=6
       List.sp[[l]]=list.append(List.sp[[l]],
             Prior.mean.Fo=0.01,
             Prior.SD.Log.Fo=0.5,
                  
             #Steepness
             h.M.constant=0.481,  #Braccini et al 2015 
             h.M.constant.low=0.461,    #80% percentile
             h.M.constant.up=0.5,
                 
             #Natural mortality
             M_val=0.283,          #Walker empirical
             M_val.low=0.22,      #using Hoenig set to Max age=19 exp(1.46-1.01*log(19))
             M_val.high=0.32,     #using Hoenig set to Max age=13 exp(1.46-1.01*log(13))
                  
             #Initial F
             Fo=0.05,               #This leaves B1975 at 95% unfished 
             Fo_Simp=0.003,              
             Fo_M=0.1,                
             Fo_AS=0.004,            #This leaves B1975 at 95% unfished 
               
             Po_spm=0.95,  #Po for surplus production, consistent with the Fo value used in Size based model
                  
             #Data
             AREAS=c("West","Zone1","Zone2"),  #Define spatial areas; 1 is West, 2 is Zn1, 3 is Zn2.
             Yr_q_change=0,   #last year before targeting practices changed (Simpfendorfer 2000)
             Yr_q_daily=2006,
             Do_var=0,     #How to calculate cpue variance in Simpfendorfer's age-structured
                   
             #Initial value of estimable parameters
             #note: Dummy is used for switching on/off phase estimation in simulation testing
             Dummy=1, 
              #---Biomass dynamics
             r_BD=fn.ji(0.3),
             k_BD=fn.ji(5000),
             Q1_BD=fn.ji(1e-7),
             Q2_BD=fn.ji(1e-7),
             daily_BD=fn.ji(1e-7),
             tau2_BD=fn.ji(0.1353),
              #---Age-structured
             RSTAR_AS=100,
             z_AS=1,
             Q1_AS=1.4,
             Q2_AS=fn.ji(0.8),
             q_daily_AS=fn.ji(0.8),
              #---Size-base
             RZERO_in_1000_individuals_SB=fn.ji(1000),
             Q1_SB=fn.ji(1e-4),
             Q2_SB=fn.ji(1e-4),
             Q_daily_SB=fn.ji(1e-4),
             Fo_SB=NA,   #no jit because it's fixed
             tau_SB=fn.ji(0.3),
             K.F=fn.ji(0.15),
             Linf.F=fn.ji(180),
             K.M=fn.ji(0.25),
             Linf.M=fn.ji(150),
             SD.growth_SB=fn.ji(10),
             Prop.west=fn.ji(0.02),  #proportion of recruitment by zone calculated based on average prop catch by zone
             Prop.zn1=fn.ji(0.08),
              #Movement 
             p11=0.999,
             p22=0.999,
             p21=0.00024,
             p33=0.999,

             #Estimation phases
             Par.phases=list('Base case'=c(dummy=Fz.off,lnR_zero=3,lnR_prop_west=-3,lnR_prop_zn1=-3,
                                           lnq=4,lnq2=Fz.off,log_Qdaily=4,ln_Init_F=Fz.off,log_tau=5,
                                            k=1,lnLinf=1,k_M=1,lnLinf_M=1,sd_growth=2,
                                            log_p11=1,log_p22=1,log_p21=1,log_p33=1),
                              S1=c(dummy=Fz.off,r=2,ln_k=2,ln_q1=1,ln_q2=Fz.off,ln_qdaily=1,ln_tau=3),
                              S2=c(dummy=Fz.off,Rstar=1,z=1,q1=1,q2=Fz.off,qdaily=Fz.off,Fo=-4))
                  
            )
       
       #Scenarios
       List.sp[[l]]$N.Scens=n.scen  
       Zens=paste("S",1:(n.scen-1),sep="")
       Models=c("Base case",paste("S",1:(n.scen-1),sep=""))
       Q.scen=c(rep("two",2),"one","N/A","two","two")
                
       List.sp[[l]]=with(List.sp[[l]],
       {list.append(List.sp[[l]],
             #--Integrated model      
          Tabla.scen=data.frame(
                   Model=Models,
                   Size_comp.=c('Yes',"N/A",'No','Yes','No','Yes'),
                   CPUE=c(rep("Stand.",3),"No","Stand.","Stand.hours"),
                   Age.Growth=c('Yes',"N/A",'No','Yes','No','Yes'),
                   Ktch.sx.r=c('Observed','N/A','Equal','Observed','Equal','Observed'),
                   Tagging=c('No','N/A',rep('No',4)),
                   Fec.=c(rep('N/A',2),'constant','N/A','constant','N/A'),
                   Maturity=c('at length','N/A','knife edge',"at length",'knife edge',"at length"),
                   M=c("constant","N/A",rep("constant",4)),
                   M.value=c(M_val,NA,rep(M_val,4)),
                   SteepnesS=c(h.M.constant,rep("N/A",2),
                   h.M.constant,"N/A",h.M.constant),
                   Q=Q.scen,   
                   Spatial_structure=rep('Single zone',6),
                   Movement=c("No","N/A",rep("No",4)),
                   Fo=c(Fo,"N/A",Fo_AS,Fo,Fo_AS,Fo),
                   Model_type=c('Length-based','Biomass dynamics',"Age-structured",
                  "Length-based","Age-structured","Length-based")),
              #--CMSY
           ktch_msy_scen=list(
                  'Base case'=list(r.prior="USER",user="Yes",k.max=50,startbio=c(0.8,.95),
                                   finalbio=c(0.2, 0.7),res="low",
                                   niter=SIMS,sigR=Proc.err),
                   S1=list(r.prior="USER",user="Yes",k.max=50,startbio=c(0.8,.95),
                           finalbio=c(0.2, 0.7),res="low",
                           niter=SIMS,sigR=0.02),
                   S2=list(r.prior=NA,user="No",k.max=50,startbio=c(0.8,.95),
                           finalbio=c(0.2, 0.7),res="low",
                           niter=SIMS,sigR=Proc.err)))
        })
    }
    
    if(List.sp[[l]]$Name=="dusky shark")   #update all input pars, currently using gummies!
    {
       n.scen=6
       List.sp[[l]]=list.append(List.sp[[l]],
               Prior.mean.Fo=0.01,
               Prior.SD.Log.Fo=0.5,
                                         
               #Steepness
               h.M.constant=0.481,  #Braccini et al 2015 
               h.M.constant.low=0.461,    #80% percentile
               h.M.constant.up=0.5,
                                           
               #Natural mortality
               M_val=0.283,          #Walker empirical
               M_val.low=0.22,      #using Hoenig set to Max age=19 exp(1.46-1.01*log(19))
               M_val.high=0.32,     #using Hoenig set to Max age=13 exp(1.46-1.01*log(13))
                                           
               #Initial F
               Fo=0.05,               #This leaves B1975 at 95% unfished 
               Fo_Simp=0.003,              
               Fo_M=0.1,                
               Fo_AS=0.004,            #This leaves B1975 at 95% unfished 
                                          
               Po_spm=0.95,  #Po for surplus production, consistent with the Fo value used in Size based model
                                          
               #Data
               Ktch.source="WA.only",  #select whether to use all catch series or only WA
               #Ktch.source="ALL",
               AREAS=c("West","Zone1","Zone2"),  #Define spatial areas; 1 is West, 2 is Zn1, 3 is Zn2.
               Yr_q_change=0,   #last year before targeting practices changed (Simpfendorfer 2000)
               Yr_q_daily=2006,
               Do_var=0,     #How to calculate cpue variance in Simpfendorfer's age-structured
                                           
               #Initial value of estimable parameters
               #note: Dummy is used for switching on/off phase estimation in simulation testing
               Dummy=1,   
                                           
               #---Biomass dynamics
              r_BD=fn.ji(0.3),
              k_BD=fn.ji(5000),
              Q1_BD=fn.ji(1e-7),
              Q2_BD=fn.ji(1e-7),
              daily_BD=fn.ji(1e-7),
              tau2_BD=fn.ji(0.1353),
                                           
               #---Age-structured
              RSTAR_AS=100,
              z_AS=1,
              Q1_AS=1.4,
              Q2_AS=fn.ji(0.8),
              q_daily_AS=fn.ji(0.8),
                               
               #---Size-base
              RZERO_in_1000_individuals_SB=fn.ji(1000),
              Q1_SB=fn.ji(1e-4),
              Q2_SB=fn.ji(1e-4),
              Q_daily_SB=fn.ji(1e-4),
              Fo_SB=NA,   #no jit because it's fixed
              tau_SB=fn.ji(0.3),
              K.F=fn.ji(0.15),
              Linf.F=fn.ji(180),
              K.M=fn.ji(0.25),
              Linf.M=fn.ji(150),
              SD.growth_SB=fn.ji(10),
              Prop.west=fn.ji(0.02),  #proportion of recruitment by zone calculated based on average prop catch by zone
              Prop.zn1=fn.ji(0.08),
                                           
               #Movement 
               p11=0.999,
               p22=0.999,
               p21=0.00024,
               p33=0.999,
                            
            #Estimation phases
            Par.phases=list('Base case'=c(dummy=Fz.off,lnR_zero=3,lnR_prop_west=-3,lnR_prop_zn1=-3,
                                          lnq=4,lnq2=Fz.off,log_Qdaily=4,ln_Init_F=Fz.off,log_tau=5,
                                           k=1,lnLinf=1,k_M=1,lnLinf_M=1,sd_growth=2,
                                           log_p11=1,log_p22=1,log_p21=1,log_p33=1),
                            S1=c(dummy=Fz.off,r=2,ln_k=2,ln_q1=1,ln_q2=Fz.off,ln_qdaily=1,ln_tau=3),
                            S2=c(dummy=Fz.off,Rstar=1,z=1,q1=1,q2=Fz.off,qdaily=Fz.off,Fo=-4))
          )
       
       #Scenarios
       List.sp[[l]]$N.Scens=n.scen  
       Zens=paste("S",1:(n.scen-1),sep="")
       Models=c("Base case",paste("S",1:(n.scen-1),sep=""))
       Q.scen=c(rep("two",2),"one","N/A","two","two")
       
       List.sp[[l]]=with(List.sp[[l]],
      {list.append(List.sp[[l]],
             #--Integrated model      
          Tabla.scen=data.frame(
          Model=Models,
          Size_comp.=c('Yes',"N/A",'No','Yes','No','Yes'),
          CPUE=c(rep("Stand.",3),"No","Stand.","Stand.hours"),
          Age.Growth=c('Yes',"N/A",'No','Yes','No','Yes'),
          Ktch.sx.r=c('Observed','N/A','Equal','Observed','Equal','Observed'),
          Tagging=c('No','N/A',rep('No',4)),
          Fec.=c(rep('N/A',2),'constant','N/A','constant','N/A'),
          Maturity=c('at length','N/A','knife edge',"at length",'knife edge',"at length"),
          M=c("constant","N/A",rep("constant",4)),
          M.value=c(M_val,NA,rep(M_val,4)),
          SteepnesS=c(h.M.constant,rep("N/A",2),h.M.constant,"N/A",h.M.constant),
          Q=Q.scen,   
          Spatial_structure=rep('Single zone',6),
          Movement=c("No","N/A",rep("No",4)),
          Fo=c(Fo,"N/A",Fo_AS,Fo,Fo_AS,Fo),
          Model_type=c('Length-based','Biomass dynamics',"Age-structured",
                       "Length-based","Age-structured","Length-based")),
            #--CMSY
          ktch_msy_scen=list(
              'Base case'=list(r.prior="USER",user="Yes",k.max=50,startbio=c(0.7,.95),
                               finalbio=c(0.2, 0.6),res="Very low",niter=SIMS,sigR=Proc.err),
              S1=list(r.prior="USER",user="Yes",k.max=50,startbio=c(0.7,.95),
                      finalbio=c(0.2, 0.6),res="Very low",niter=SIMS,sigR=0.02),
              S2=list(r.prior=NA,user="No",k.max=50,startbio=c(0.7,.95),
                      finalbio=c(0.2, 0.6),res="Very low",niter=SIMS,sigR=Proc.err)))
        })

                  
    }
    
    if(List.sp[[l]]$Name=="sandbar shark")   #update all input pars, currently using gummies!
    {
      n.scen=6
      List.sp[[l]]=list.append(List.sp[[l]],
          Prior.mean.Fo=0.01,
          Prior.SD.Log.Fo=0.5,
                               
          #Steepness
          h.M.constant=0.481,  #Braccini et al 2015 
          h.M.constant.low=0.461,    #80% percentile
          h.M.constant.up=0.5,
                               
          #Natural mortality
          M_val=0.283,          #Walker empirical
          M_val.low=0.22,      #using Hoenig set to Max age=19 exp(1.46-1.01*log(19))
          M_val.high=0.32,     #using Hoenig set to Max age=13 exp(1.46-1.01*log(13))
                               
          #Initial F
          Fo=0.05,               #This leaves B1975 at 95% unfished 
          Fo_Simp=0.003,              
          Fo_M=0.1,                
          Fo_AS=0.004,            #This leaves B1975 at 95% unfished 
                               
          Po_spm=0.95,  #Po for surplus production, consistent with the Fo value used in Size based model
                               
          #Data
          Ktch.source="WA.only",  #select whether to use all catch series or only WA
          #Ktch.source="ALL",
          AREAS=c("West","Zone1","Zone2"),  #Define spatial areas; 1 is West, 2 is Zn1, 3 is Zn2.
          Yr_q_change=0,   #last year before targeting practices changed (Simpfendorfer 2000)
          Yr_q_daily=2006,
          Do_var=0,     #How to calculate cpue variance in Simpfendorfer's age-structured
                               
          #Initial value of estimable parameters
          #note: Dummy is used for switching on/off phase estimation in simulation testing
          Dummy=1,   
                               
           #---Biomass dynamics
          r_BD=fn.ji(0.3),
          k_BD=fn.ji(5000),
          Q1_BD=fn.ji(1e-7),
          Q2_BD=fn.ji(1e-7),
          daily_BD=fn.ji(1e-7),
          tau2_BD=fn.ji(0.1353),
                               
           #---Age-structured
          RSTAR_AS=100,
          z_AS=1,
          Q1_AS=1.4,
          Q2_AS=fn.ji(0.8),
          q_daily_AS=fn.ji(0.8),
                               
           #---Size-base
          RZERO_in_1000_individuals_SB=fn.ji(1000),
          Q1_SB=fn.ji(1e-4),
          Q2_SB=fn.ji(1e-4),
          Q_daily_SB=fn.ji(1e-4),
          Fo_SB=NA,   #no jit because it's fixed
          tau_SB=fn.ji(0.3),
          K.F=fn.ji(0.15),
          Linf.F=fn.ji(180),
          K.M=fn.ji(0.25),
          Linf.M=fn.ji(150),
          SD.growth_SB=fn.ji(10),
          Prop.west=fn.ji(0.02),  #proportion of recruitment by zone calculated based on average prop catch by zone
          Prop.zn1=fn.ji(0.08),
                               
           #Movement 
          p11=0.999,
          p22=0.999,
          p21=0.00024,
          p33=0.999,
                               
          #Estimation phases
          Par.phases=list(
            'Base case'=c(dummy=Fz.off,lnR_zero=3,lnR_prop_west=-3,lnR_prop_zn1=-3,
                          lnq=4,lnq2=Fz.off,log_Qdaily=4,ln_Init_F=Fz.off,log_tau=5,
                          k=1,lnLinf=1,k_M=1,lnLinf_M=1,sd_growth=2,
                          log_p11=1,log_p22=1,log_p21=1,log_p33=1),
            S1=c(dummy=Fz.off,r=2,ln_k=2,ln_q1=1,ln_q2=Fz.off,ln_qdaily=1,ln_tau=3),
            S2=c(dummy=Fz.off,Rstar=1,z=1,q1=1,q2=Fz.off,qdaily=Fz.off,Fo=-4))
      )
      
      #Scenarios
      List.sp[[l]]$N.Scens=n.scen  
      Zens=paste("S",1:(n.scen-1),sep="")
      Models=c("Base case",paste("S",1:(n.scen-1),sep=""))
      Q.scen=c(rep("two",2),"one","N/A","two","two")
      
      List.sp[[l]]=with(List.sp[[l]],
      {list.append(List.sp[[l]],
           #--Integrated model      
         Tabla.scen=data.frame(
         Model=Models,
         Size_comp.=c('Yes',"N/A",'No','Yes','No','Yes'),
         CPUE=c(rep("Stand.",3),"No","Stand.","Stand.hours"),
         Age.Growth=c('Yes',"N/A",'No','Yes','No','Yes'),
         Ktch.sx.r=c('Observed','N/A','Equal','Observed','Equal','Observed'),
         Tagging=c('No','N/A',rep('No',4)),
         Fec.=c(rep('N/A',2),'constant','N/A','constant','N/A'),
         Maturity=c('at length','N/A','knife edge',"at length",'knife edge',"at length"),
         M=c("constant","N/A",rep("constant",4)),
         M.value=c(M_val,NA,rep(M_val,4)),
         SteepnesS=c(h.M.constant,rep("N/A",2),h.M.constant,"N/A",h.M.constant),
         Q=Q.scen,   
         Spatial_structure=rep('Single zone',6),
         Movement=c("No","N/A",rep("No",4)),
         Fo=c(Fo,"N/A",Fo_AS,Fo,Fo_AS,Fo),
         Model_type=c('Length-based','Biomass dynamics',"Age-structured",
                      "Length-based","Age-structured","Length-based")),
           #--CMSY
         ktch_msy_scen=list(
            'Base case'=list(r.prior="USER",user="Yes",k.max=50,startbio=c(0.85,.95),
                             finalbio=c(0.2, 0.6),res="Very low",niter=SIMS,sigR=Proc.err),
             S1=list(r.prior="USER",user="Yes",k.max=50,startbio=c(0.85,.95),
                     finalbio=c(0.2, 0.6),res="Very low",niter=SIMS,sigR=0.02),
             S2=list(r.prior=NA,user="No",k.max=50,startbio=c(0.85,.95),
                     finalbio=c(0.2, 0.6),res="Very low",niter=SIMS,sigR=Proc.err)))
      })
    }
    
    
    n.areas=length(List.sp[[l]]$AREAS)
    List.sp[[l]]$n.areas=n.areas
    List.sp[[l]]$Areas.zones=data.frame(area=1:n.areas,zone=List.sp[[l]]$AREAS)
    
    List.sp[[l]]$hndl=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",sep='')
    
    List.sp[[l]]$Fo_SB=List.sp[[l]]$Fo  #fixed
    
    #2. Create folders if new run
    kriat.path=paste(List.sp[[l]]$hndl,AssessYr,sep="")
    if(!dir.exists(kriat.path))dir.create(kriat.path)
    kriate.this=as.character(List.sp[[l]]$Tabla.scen$Model)
    for (i in 1:length(kriate.this)) 
    { 
      NEW=paste(kriat.path,kriate.this[i], sep="/")
      if(!dir.exists(NEW))dir.create(NEW) 
    }
    
    #drop conv tag if not using
    if(add.conv.tag=="NO") List.sp[[l]]$Tabla.scen=List.sp[[l]]$Tabla.scen[,-match("Tagging",names(List.sp[[l]]$Tabla.scen))]
    
    
    #3. Export scenarios table
    if(First.run=="YES")
    {
      inpt.pz=paste(List.sp[[l]]$hndl,AssessYr,sep="")
      setwd(inpt.pz)
      if(List.sp[[l]]$Name=="whiskery shark")
      {
        THiS=c("Model_name","Model_type","Spatial_structure"
               ,"Movement","Size_comp.","CPUE","CPUE_years_dropped","Age.Growth"         
               ,"Ktch.sx.r","Tagging","M.value","Steepness","Fo","Maturity","Q")
        HDR.span=c(2,1,1,6,4,1)
        HDR.2nd=c("Name","Type","structure",'',"Size","CPUE","CPUE years","Age &",
                  "Prop. male","Tagging","M","h","Fo","Maturity","")
        HDR.3rd=c("","","","","composition","","not used in likelihood","growth","in catch",
                  "","","","","","")
      }else
      {
        THiS=c("Model_name","Model_type","Spatial_structure"
               ,"Movement","Size_comp.","CPUE","Age.Growth"         
               ,"Ktch.sx.r","Tagging","M.value","Steepness","Fo","Maturity","Q")
        HDR.span=c(2,1,1,5,4,1)
        HDR.2nd=c("Name","Type","structure",'',"Size","CPUE","Age &",
                  "Prop. male","Tagging","M","h","Fo","Maturity","")
        HDR.3rd=c("","","","","composition","","growth","in catch",
                  "","","","","","")
      }
      
      Tabla.scen.show=List.sp[[l]]$Tabla.scen
      Tabla.scen.show=Tabla.scen.show[,-match(c("Fec.","M"),colnames(Tabla.scen.show))]
      Tabla.scen.show[is.na(Tabla.scen.show)]="N/A"
      names(Tabla.scen.show)[match(c("Model","SteepnesS"),names(Tabla.scen.show))]=c("Model_name","Steepness")
      Tabla.scen.show=Tabla.scen.show[,match(THiS,names(Tabla.scen.show))]
      Tabla.scen.show=Tabla.scen.show[match(c(Zens,"Base case" ),Tabla.scen.show$Model_name),]
      Tabla.scen.show$Q=as.character(Tabla.scen.show$Q)
      Tabla.scen.show$Q=with(Tabla.scen.show,ifelse(Q=="three",3,ifelse(Q=="two",2,ifelse(Q=="one",1,Q))))
      Tabla.scen.show$Q=with(Tabla.scen.show,ifelse(!Q=="N/A"& Q>1,paste(Q,"periods"),
                                                    ifelse(!Q=="N/A"& Q==1,paste(Q,"period"),Q))) 
      Scenarios.tbl=function(WD,Tbl,Doc.nm,caption,paragph,HdR.col,HdR.bg,Hdr.fnt.sze,Hdr.bld,
                             body.fnt.sze,Zebra,Zebra.col,Grid.col,Fnt.hdr,Fnt.body,
                             HDR.names,HDR.span,HDR.2nd,HDR.3rd)
      {
        mydoc = docx(Doc.nm)  #create r object
        mydoc = addSection( mydoc, landscape = T )   #landscape table
        # add title
        if(!is.na(caption))mydoc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
        
        # add a paragraph
        if(!is.na(paragph))mydoc = addParagraph(mydoc , paragph, stylename="Citationintense")
        
        #add table
        MyFTable=FlexTable(Tbl,header.column=F,add.rownames =F,
                           header.cell.props = cellProperties(background.color=HdR.bg), 
                           header.text.props = textProperties(color=HdR.col,font.size=Hdr.fnt.sze,
                                                              font.weight="bold",font.family =Fnt.hdr), 
                           body.text.props = textProperties(font.size=body.fnt.sze,font.family =Fnt.body))
        
        #Add header
        MyFTable = addHeaderRow(MyFTable,text.properties=textBold(),value=HDR.names,colspan=HDR.span)
        
        #Add second header
        MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value =HDR.2nd)
        
        #Add third header
        MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value =HDR.3rd)
        
        # zebra stripes - alternate colored backgrounds on table rows
        if(Zebra=="YES") MyFTable = setZebraStyle(MyFTable, odd = Zebra.col, even = "white" )
        
        # table borders
        MyFTable = setFlexTableBorders(MyFTable,
                                       inner.vertical = borderNone(),inner.horizontal = borderNone(),
                                       outer.vertical = borderNone(),
                                       outer.horizontal = borderProperties(color=Grid.col, style="solid", width=4))
        
        # set columns widths (in inches)
        #MyFTable = setFlexTableWidths( MyFTable, widths = Col.width)
        
        mydoc = addFlexTable( mydoc, MyFTable)   
        mydoc = addSection( mydoc, landscape = F ) 
        
        # write the doc 
        writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
      }
      options('ReporteRs-fontsize'= 12, 'ReporteRs-default-font'='Arial')   
      Scenarios.tbl(WD=getwd(),Tbl=Tabla.scen.show,Doc.nm="Model scenarios",
                    caption=NA,paragph=NA,HdR.col='black',HdR.bg='white',
                    Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                    Zebra='NO',Zebra.col='grey60',Grid.col='black',
                    Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
                    HDR.names=c('Model','Spatial','Movement', 'Data','Input parameters','Q'),
                    HDR.span=HDR.span,
                    HDR.2nd=HDR.2nd,
                    HDR.3rd=HDR.3rd)
    }
    
    
    #4. Create pin file
    with(List.sp[[l]],
    {
           #Population pin values
           Pin.pars=vector('list',nrow(Tabla.scen))
           names(Pin.pars)=Tabla.scen$Model
           
           #Biomass dynamics
           Pin.pars$S1=c(Dummy=Dummy,r=r_BD,log_k=log(k_BD),Q1=Q1_BD,Q2=Q2_BD,Qdaily=daily_BD,log_tau2=log(tau2_BD))
           
           #Age-structured
           Pin.pars$S2=c(Dummy=Dummy,RSTAR=RSTAR_AS,Z=z_AS,Q1=Q1_AS,Q2=Q2_AS,Qdaily=q_daily_AS,Fo=Fo_AS) 
           
           #Length-based
           Pin.pars$'Base case'=c(Dummy=Dummy,lnR_zero=log(RZERO_in_1000_individuals_SB),
                                  R_prop_west=Prop.west,R_prop_zn1=Prop.zn1,
                                  lnq=log(Q1_SB),lnq2=log(Q2_SB),lnq_daily=log(Q_daily_SB),
                                  ln_Init_F=log(Fo_SB),ln_tau=log(tau_SB),
                                  k=K.F,lnLinf=log(Linf.F),k_M=K.M,lnLinf_M=log(Linf.M),
                                  sd_growth=SD.growth_SB)
           
           #Add tagging pars       
           if(conv.tag.all=="YES") 
           {
             if(Move.mode=="Individual-based") 
             {
               Mov.pars=c(p11=p11,p22=p22,p21=p21,p33=p33)
             }
             
             if(Move.mode=="Population-based")
             {
               Mov.pars=c(mov.WC_WC=mov.WC_WC,mov.WC_ZN1=mov.WC_ZN1,
                          mov.ZN1_WC=mov.ZN1_WC,mov.ZN1_ZN1=mov.ZN1_ZN1,
                          mov.ZN2_ZN1=mov.ZN2_ZN1,mov.ZN2_ZN2=mov.ZN2_ZN2,
                          log.Q.tag.WC=log.Q.tag.WC,log.Q.tag.ZN1=log.Q.tag.ZN1,
                          log.Q.tag.ZN2=log.Q.tag.ZN2,log.tau=log.tau)
             }
             Pin.pars$'Base case'=c(Pin.pars$'Base case',Mov.pars)
           }
           
           #Set the pins from other size-based scenarios
           IDs=which(Tabla.scen$Model_type=="Length-based")
           IDs=IDs[-match("Base case",Tabla.scen$Model[IDs])]
           for(i in IDs) Pin.pars[[i]]=Pin.pars$'Base case'
           
           for(i in 1:N.Scens)     
           {
             if(Tabla.scen$Model_type[i]=="Length-based")
             {
               #higher M or Fo needs larger initial population
               if(Tabla.scen$M.value[i]==M_val.high | Tabla.scen$Fo[i]==Fo_M)    
               {
                 ss=match("lnR_zero",names(Pin.pars[[i]]))
                 Pin.pars[[i]][ss]=Pin.pars[[i]][ss]*1.15
               }
               
               #no movement outside zone
               if(Tabla.scen$Movement[i]%in%c("No","N/A"))
               {
                 Pin.pars[[i]][match(c("p11","p22","p21","p33"),names(Pin.pars[[i]]))]=c(1,1,0,1)      
               }
               
               #Change Fo if required
               if(!as.numeric(as.character(Tabla.scen$Fo[i]))==Fo) 
               {
                 Pin.pars[[i]][match("ln_Init_F",names(Pin.pars[[i]]))]=log(as.numeric(as.character(Tabla.scen$Fo[i])))
                 ss=match("lnR_zero",names(Pin.pars[[i]]))
                 #higher Fo needs larger initial population
                 if(!as.numeric(as.character(Tabla.scen$Fo[i]))>Fo)Pin.pars[[i]][ss]=Pin.pars[[i]][ss]*1.15
               }
               
             }
           }
           
           #Export as pin file
           setPath=function(Scen)setwd(paste(List.sp[[l]]$hndl,AssessYr,"/",Scen,sep=""))
           for(i in 1:nrow(Tabla.scen))
           {
             These.pars=Pin.pars[[i]]
             par.nms=names(These.pars)
             setPath(Tabla.scen[i,]$Model)
             FILE=paste(List.sp[[l]]$Name.inputs,".pin",sep="")
             write("# Input parameters",file = FILE)
             for(k in 1:length(These.pars))
             {
               Hdr=paste("#",par.nms[k])
               write(Hdr,file = FILE,append=T)
               write(These.pars[k],file = FILE,sep = "\t",append=T)
             }
           }
           List.sp[[l]]$Pin.pars=Pin.pars
           
           
           #Set estimation phases for scenarios not already assigned
           IDs=names(Pin.pars)[-match(c("Base case","S1","S2"),names(Pin.pars))]
           dummy=vector('list',length(IDs))
           names(dummy)=IDs
           for(i in 1:length(dummy)) dummy[[i]]=Par.phases$'Base case'
           Par.phases=c(Par.phases,dummy)
           
           #Turn off irrelevant pars according to scenarios
           for(i in 1:N.Scens)     
           {
             if(Tabla.scen$Model_type[i]=="Length-based")
             {
               #single zone
               if(Tabla.scen$Spatial_structure[i]=="Single zone")
               {
                 Par.phases[[i]][fn.mtch(Zns.par.phz,Par.phases[[i]])]=rep(Fz.off,length(Zns.par.phz))
                 Par.phases[[i]][fn.mtch(Q_phz,Par.phases[[i]])]=c(1,1,1)
               }
               #effective cpue
               if(Tabla.scen$CPUE[i]=="Effective")
               {
                 Par.phases[[i]][fn.mtch(Q_phz,Par.phases[[i]])]=c(1,1,Fz.off)
                 Par.phases[[i]][fn.mtch("log_tau",Par.phases[[i]])]=Fz.off
               }
               #no Qs if no cpue used
               if(Tabla.scen$CPUE[i]%in%c("N/A","No"))
               {
                 Par.phases[[i]][fn.mtch(c(Q_phz),Par.phases[[i]])]=rep(Fz.off,length(c(Q_phz)))
                 Par.phases[[i]][fn.mtch("lnR_zero",Par.phases[[i]])]=1
               }
               #no movement
               if(Tabla.scen$Movement[i]%in%c("No","N/A"))
               {
                 Par.phases[[i]][fn.mtch(MOv.par.phz,Par.phases[[i]])]=rep(Fz.off,length(MOv.par.phz))  
               }
               
               #one less q
               if(Tabla.scen$Q[i]=="two" & Tabla.scen$CPUE[i]=="Stand.")
               {
                 Par.phases[[i]][fn.mtch("lnq2",Par.phases[[i]])]=-Par.phases[[i]][fn.mtch("lnq2",Par.phases[[i]])]
               }
             }
           }
           List.sp[[l]]$Par.phases=Par.phases
         })
    
    
    #5. Define general modelling arguments
    List.sp[[l]]=list.append(List.sp[[l]],
                             
        #Select acoustic tagging model
      Acoust.format=Move.mode,
      #Acoust.format="SS3"
      
        #Select type of size composition Likelihood
      Size_like=size.likelihood,
      Dirichlet.small.value=1e-4,   #small constant for intermediate 0 observations 
    
        #Combined size composition?  
      Size.sex.comb=size.sex.comb,
    
        #Size compostion as proportions?
      Size.comp.prop=size.comp.prop,
    
        #Number of years for future projections
      Yrs.future=yrs.future,
    
        #Effective sample size size composition
      Effective.n=effective.n,
    
        #Maximum possible F
      MaxF=maxF,
      
        #Prior for initial fishing mortality
      add_Finit_prior=Add_Finit_prior,
      
        #running arguments
      Arguments=arguments,
      
      #Do MSY calculation using Base Case model   
      Do.MSY=do.MSY,
      MSY.yrs=msy.yrs,  
      MSY.sd.rec=msy.sd.rec,
      MSY.sims=msy.sims,
      F.vec=f.vec,
      
      #MCMC 
      DO.MCMC=mcmc.do,
      burning=1:(5*length(seq(1,mcmc.n,by=mcmc.thin))/100),   #5%  burning
      
      #What scenarios to run
      Run.all.Scenarios=run.all.scenarios,
      
      #Future projection scenarios
      Run.future.proj=run.future,
      Future.ktch.scen=list(mean=1,percent_50=.5,percent_150=1.5)
    )
    
    #Weight of likelihood components
    List.sp[[l]]$Rho=1  #weight of size comp 
    List.sp[[l]]$Rho2=1 #weight of age-length
    x=rep(1,length(List.sp[[l]]$Models)) #cpue likelihood
    names(x)=as.character(List.sp[[l]]$Models)
    dd=which(List.sp[[l]]$Tabla.scen$CPUE%in%c("N/A","No")) #adjust weight of cpue if cpue like not used in model
    if(length(dd)>0)x[dd]=0
    List.sp[[l]]$Rho3=x
    
    #do future projections with 0 catch?
    if(Do.zero.Ktch=="YES")
    {
      List.sp[[l]]$zero.Ktch.yrs=100
      List.sp[[l]]$zero.Ktch.this="Base case"
    }
    
    #do simulation testing?
    if(Do.sim.test=="YES")
    {
      List.sp[[l]]$N.sim.test=50
      List.sp[[l]]$CPUE.jitr=1000 
      List.sp[[l]]$Sim.Test.this="Base case"
    }

  }
  
  #'other species' arguments
  if(!List.sp[[l]]$Species%in%Indicator.species)
  {
    List.sp[[l]]=list.append(List.sp[[l]],
                             
      #... Surplus production arguments
      #note: Only fitting species with species-specific abundance time series 
      #     Assumption, negligible exploitation at start of time series
      
      #Initial harvest rate 
      HR_o.sd=0.005,  #SD of HR likelihood (fixed)
      
      #Efficiency increase scenarios from 1995 on (done up to 1994 in cpue stand.)
      Efficien.scens=c(0),
      #Efficien.scens=c(.01),
      
      #Proportional biomass (as proportion of K) at start of catch time series
      B.init=1, #(fixed)  Starting @ virgin level
      
      #Estimate q
      estim.q="NO",   #use Haddon's q MLE calculation
      
      #cpue likelihood
      what.like='kernel',
      #what.like='full',
      
      #number of MC simulations
      N.monte=1000,
      
      #Initial estimated par value
      Init.r=with(List.sp[[l]],
                   case_when(Name=="copper shark"~.05,
                             Name=="great hammerhead"~.1,
                             Name=="grey nurse shark"~0.05,
                             Name=="lemon shark"~.1,
                             Name=="milk shark"~.2,
                             Name=="pigeye shark"~0.1,
                             Name=="sawsharks"~.1,
                             Name=="scalloped hammerhead"~.1,
                             Name=="smooth hammerhead"~.1,
                             Name=="spinner shark"~.1,
                             Name=="shortfin mako"~.05,
                             Name=="spurdogs"~.05,
                             Name=="tiger shark"~.1,
                             Name=="wobbegongs"~.1,
                             TRUE~NA_real_)),
      
      #maximum acceptable CV for cpue series  
      MAX.CV=0.5,   
      
      #Define which optimisation method to use
      #minimizer='nlminb',
      minimizer='optim',
      
      Remove.bounds=FALSE,
      
      usePen=TRUE,  
      
      #K bounds
      Low.bound.K=1,  #times the maximum catch
      Up.bound.K=100, 
      
      #K init times max ktch
      k.times.mx.ktch=mean(c(1,100)),
      
      #fix or estimate r
      fix.r="NO",
      r.weight=1,   #weight given in the likelihood function
      
      #what biomass percentiles to show
      What.percentil="100%", #100% to make it comparable to CMSY  
      #What.percentil="60%" #60% as required for MSC

            
      #... Future projections
      years.futures=5
    )
    
    #... Catch-MSY arguments
    List.sp[[l]]$SIMS=5e4   #simulatins
    List.sp[[l]]$ERROR=0   #Assumed process error
    #depletion level at start of catch series
    List.sp[[l]]$STARTBIO=c(List.sp[[l]]$B.init*.95,List.sp[[l]]$B.init)   #low depletion because starting time series prior to any fishing
    List.sp[[l]]$FINALBIO=c(.2,.9)       #very uncertain
    List.sp[[l]]$Do.sim.test="NO"  #simulation test Catch-MSY for small and large catches
    
    
    
    #... CMSY Scenarios considered
    List.sp[[l]]$SCENARIOS=list(Error=List.sp[[l]]$ERROR,
                                R.prior=List.sp[[l]]$r.prior,
                                Initial.dep=List.sp[[l]]$STARTBIO)
    
  }
  
  
  
}

#Export table of life history parameters
if(First.run=="YES")
{
  #Export Life history table
  Rar.path=paste(handl_OneDrive('Reports/RARs'), AssessYr,sep="/")
  if(!dir.exists(Rar.path))dir.create(Rar.path)
  setwd(Rar.path)
  TabL=LH.data%>%
          filter(SPECIES%in%sapply(List.sp, "[[", "Species"))%>%
          dplyr::select(-SNAME)%>%
          left_join(All.species.names%>%dplyr::select(SPECIES,SNAME,Scien.nm),by="SPECIES")%>%
          mutate(K=round(as.numeric(as.character(K)),3),
                 FL_inf=round(FL_inf),
                 SNAME=capitalize(SNAME))%>%
    rename(Species=SNAME)%>%
    dplyr::select(-SPECIES,-Comment)%>%
    relocate(Species,Scien.nm)%>%
    arrange(Species)
  fn.word.table(WD=getwd(),TBL=TabL,Doc.nm="Table 1. Life history pars",caption=NA,paragph=NA,    
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  write.csv(TabL,"Table 1. Life history pars.csv",row.names = F)
}



#---10.  Export all available input data to each species assessment folder----- 
if(First.run=="YES")
{
  source(handl_OneDrive("Analyses/Population dynamics/Git_Stock.assessments/Organise data.R"))
  for(l in 1:N.sp)
  {
    print(paste("---------",names(List.sp)[l]))
    fn.input.data(Name=List.sp[[l]]$Name,
                  Name.inputs=List.sp[[l]]$Name.inputs,
                  SP=List.sp[[l]]$SP,
                  Species=List.sp[[l]]$Species,  
                  First.year=List.sp[[l]]$First.year,
                  Last.year=Last.yr.ktch,
                  Min.obs=Min.obs,
                  Min.shts=Min.shts,
                  What.Efrt=What.Effort,
                  Bin.size=TL.bins.cm,
                  Yr.assess=AssessYr,
                  Dat=Species.data[[l]],
                  LH.par=LH.data%>%filter(SPECIES==List.sp[[l]]$Species)) 
  } 
}


#---11.  Demography. r prior ----------------------------------------------------------------------- 
store.species.r=vector('list',N.sp)
names(store.species.r)=Keep.species

if(do.r.prior)
{
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
  system.time({for(l in 1:N.sp)   #takes 0.013 sec per iteration per species
  {
    print(paste("r prior ","--",List.sp[[l]]$Name))
    PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
               capitalize(List.sp[[l]]$Name),"/",AssessYr,"/demography",sep='')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))   
    setwd(PATH)
    #if no sd, replace with mean from those species with sd
    if(is.na(List.sp[[l]]$Growth.F$FL_inf.sd)) List.sp[[l]]$Growth.F$FL_inf.sd=0.038*List.sp[[l]]$Growth.F$FL_inf  
    if(is.na(List.sp[[l]]$Growth.F$k.sd)) List.sp[[l]]$Growth.F$k.sd=0.088*List.sp[[l]]$Growth.F$k     
    M.averaging<<-"min" #'min' yields rmax, Cortes pers com, but yields too high steepness for all species
    RESAMP="YES"
    
    linear.fec="NO"
    #if(names(List.sp)[l]%in%c("grey nurse shark","sandbar shark")) linear.fec="NO"
    
      
    #Get r prior
    r.prior.dist=with(List.sp[[l]],fun.rprior.dist(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf/.85,
                                   K.sd=Growth.F$k.sd,LINF.sd=Growth.F$FL_inf.sd/.85,k.Linf.cor,
                                   Amax=Max.age.F,
                                   MAT=unlist(Age.50.mat),FecunditY=Fecundity,Cycle=Breed.cycle,
                                   BWT=BwT,AWT=AwT,LO=Lzero/.85))
    #export r and M
    out.r=data.frame(shape=r.prior.dist$shape,rate=r.prior.dist$rate,
                     mean=r.prior.dist$mean,sd=r.prior.dist$sd)
    write.csv(out.r,'r.prior.csv',row.names = F)
    store.species.r[[l]]=r.prior.dist
    
    n.dim=max(unlist(lapply(r.prior.dist$M,length)))
    out.M=r.prior.dist$M
    for(ss in 1:length(out.M))
    {
      a=out.M[[ss]]
      delta=n.dim-length(a)
      if(delta>0) out.M[[ss]]=c(a,rep(NA,delta))
      rm(a)
    }
    out.M=do.call(rbind,out.M)
    names(out.M)=0:(n.dim-1)
    write.csv(out.M,"M.csv",row.names=FALSE)
    rm(r.prior.dist,M.averaging,RESAMP,linear.fec)
  }})
}else
{
  for(l in 1:N.sp) 
  {
    store.species.r[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                         AssessYr,"/demography/r.prior.csv",sep=''))
  }
}



#---12.  Assign Resilience -----------------------------------------------------------------------
RESILIENCE=vector('list',N.sp)
names(RESILIENCE)=names(List.sp)
for(r in 1:length(RESILIENCE))
{
  RESILIENCE[[r]]=with(store.species.r[[r]],ifelse(mean>=0.6,"High",
                            ifelse(mean<0.6 & mean>=0.2,'Medium',
                            ifelse(mean<0.2 & mean>=0.1,'Low',
                            'Very low'))))
}


#---13.  Extract selectivity at age and at size-----------------------------------------------------------------------
  #for species with no gillnet selectivity profile, set to closest species or family
Sel.equivalence=data.frame(
      Name=c("copper shark","great hammerhead","scalloped hammerhead",
             "grey nurse shark","wobbegongs",
             "sawsharks",
              "lemon shark","milk shark","pigeye shark","shortfin mako","spinner shark","tiger shark",
             "spurdogs"),
      Equivalence=c("Dusky shark",rep("Smooth hammerhead",2),
                    rep("Hexanchidae",2),
                    "Common sawshark",
                    rep('Carcharhinidae',6),
                    'Spikey dogfish'))
Selectivity.at.age=vector('list',N.sp)
names(Selectivity.at.age)=Keep.species
Selectivity.at.totalength=Selectivity.at.age
HandL=handl_OneDrive("Analyses/Data_outs/")
for(l in 1:N.sp)
{
  #1. Read in selectivity at age data
  if('gillnet.selectivity_len.age'%in%names(Species.data[[l]]))
  {
    GN.sel.at.age=Species.data[[l]]$gillnet.selectivity_len.age
    GN.sel.at.totalength=Species.data[[l]]$gillnet.selectivity
  }else
  {
    #allocate  selectivity from family                                  
    this.sel=Sel.equivalence%>%filter(Name==names(Species.data)[l])
    temp.wd=paste(HandL,this.sel$Equivalence,sep='')
    GN.sel.at.age=read.csv(paste(temp.wd,'/',this.sel$Equivalence,'_Gillnet.selectivity_len.age.csv',sep=''))
    GN.sel.at.totalength=read.csv(paste(temp.wd,'/',this.sel$Equivalence,'_Gillnet.selectivity.csv',sep=''))
  }

  #2. Get combined selectivity
  if(!"X16.5"%in%names(GN.sel.at.age))
  {
    GN.sel.at.age=GN.sel.at.age%>%
                    rename(X16.5='16.5',
                           X17.8='17.8')
    GN.sel.at.totalength=GN.sel.at.totalength%>%
                    rename(X16.5='16.5',
                           X17.8='17.8')
  }
  GN.sel.at.age=GN.sel.at.age%>%
                  mutate(Sum.sel=X16.5+X17.8,
                         Sel.combined=Sum.sel/max(Sum.sel),
                         Sel.combined=Sel.combined/max(Sel.combined),
                         TL=TL.mm)
  Selectivity.at.age[[l]]=GN.sel.at.age[,c('TL','Age','Sel.combined')]
  
  GN.sel.at.totalength=GN.sel.at.totalength%>%
                  mutate(Sum.sel=X16.5+X17.8,
                         Sel.combined=Sum.sel/max(Sum.sel),
                         Sel.combined=Sel.combined/max(Sel.combined),
                         TL=TL.mm)
  Selectivity.at.totalength[[l]]=GN.sel.at.totalength[,c('TL','Sel.combined')]
  
}
#display selectivities    
if(First.run=="YES")
{
  fn.fig(handl_OneDrive('Analyses\\Population dynamics\\growth.and.selectivity'),2400,2000) 
  smart.par(n.plots=N.sp,MAR=c(2,3,1,1),OMA=c(2.5,1,.05,2.5),MGP=c(1.8,.5,0))
  par(cex.lab=1.5,las=1)
  for(l in 1:N.sp)
  {
    with(Selectivity.at.age[[l]],plot(Age,Sel.combined,col=2,ylim=c(0,1),
                                      pch=19,main=capitalize(names(Selectivity.at.age)[l]),
                                      ylab='',xlab=""))
    par(new = TRUE)
    with(Selectivity.at.age[[l]],plot(Age,TL,type='l',lwd=2,
                                      xaxt = "n", yaxt = "n",
                                      ylab = "", xlab = ""))
    axis(side = 4)
    
  }
  mtext("Age", side = 1, line = 1,outer=T)
  mtext("Selectivity", side = 2, line = -.5,las=3,col=2,outer=T)
  mtext("TL (cm)", side = 4, line = 1,outer=T,las=3)
  dev.off()
  
  fn.fig(handl_OneDrive('Analyses\\Population dynamics\\growth.and.selectivity2'),2400,2000) 
  smart.par(n.plots=N.sp,MAR=c(2,3,1,1),OMA=c(2.5,1,.05,2.5),MGP=c(1.8,.5,0))
  par(cex.lab=1.5,las=1)
  for(l in 1:N.sp)
  {
    with(Selectivity.at.age[[l]],plot(TL,Sel.combined,col=2,ylim=c(0,1),
                                      pch=19,main=capitalize(names(Selectivity.at.age)[l]),
                                      ylab='',xlab=""))
    par(new = TRUE)
    with(Selectivity.at.age[[l]],plot(TL,Age,type='l',lwd=2,
                                      xaxt = "n", yaxt = "n",
                                      ylab = "", xlab = ""))
    axis(side = 4)
    
  }
  mtext("TL (cm)", side = 1, line = 1,outer=T)
  mtext("Selectivity", side = 2, line = -.5,las=3,col=2,outer=T)
  mtext("Age", side = 4, line = 1,outer=T,las=3)
  dev.off()
}

#---14. Steepness ----------------------------------------------------------------------- 
store.species.steepness=vector('list',N.sp)
names(store.species.steepness)=Keep.species

if(do.steepness)
{
  system.time(for(l in 1: N.sp) 
  {
    print(paste("steepness ","--",List.sp[[l]]$Name))
    PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
               capitalize(List.sp[[l]]$Name),"/",AssessYr,"/steepness",sep='')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
    setwd(PATH)
    
    SEL=Selectivity.at.age[[l]]$Sel.combined
    
     if(is.na(List.sp[[l]]$Growth.F$FL_inf.sd)) List.sp[[l]]$Growth.F$FL_inf.sd=0.038*List.sp[[l]]$Growth.F$FL_inf  
    if(is.na(List.sp[[l]]$Growth.F$k.sd)) List.sp[[l]]$Growth.F$k.sd=0.088*List.sp[[l]]$Growth.F$k     
    
    #Fishing mortality set at 0 so selectivity has no effect
    k.Linf.cor=List.sp[[l]]$k.Linf.cor
    
    M.averaging<<-"mean"   #'min' yields too high h values for all species
    if(names(List.sp)[l]%in%"milk shark") M.averaging<-"min" #too low h if mean, considering biology
    RESAMP="YES"
    linear.fec="NO"
    if(names(List.sp)[l]%in%c("angel sharks","grey nurse shark","spurdogs")) RESAMP="NO"   #endless loop due to life history combos <0.2


    steepNs=with(List.sp[[l]],fun.steepness(Nsims=2*NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf/.85,
                                               Linf.sd=Growth.F$FL_inf.sd/.85,k.sd=Growth.F$k.sd,
                                               first.age=0,sel.age=SEL,F.mult=0,
                                               Amax=Max.age.F,MAT=unlist(Age.50.mat),
                                               FecunditY=Fecundity,Cycle=Breed.cycle,
                                               sexratio=0.5,spawn.time = 0,
                                               AWT=AwT,BWT=BwT,LO=Lzero/.85,
                                               Resamp=RESAMP))
    
    
    
    ##
    #export h and M
    out.h=data.frame(shape=steepNs$shape,rate=steepNs$rate,
                     mean=steepNs$mean,sd=steepNs$sd)
    write.csv(out.h,'h.prior.csv',row.names = F) 
    store.species.steepness[[l]]=steepNs
    
    n.dim=max(unlist(lapply(steepNs$M,length)))
    out.M=steepNs$M
    for(ss in 1:length(out.M))
    {
      a=out.M[[ss]]
      delta=n.dim-length(a)
      if(delta>0) out.M[[ss]]=c(a,rep(NA,delta))
      rm(a)
    }
    out.M=do.call(rbind,out.M)
    names(out.M)=0:(n.dim-1)
    write.csv(out.M,"M.csv",row.names=FALSE)
    rm(steepNs,RESAMP,M.averaging,out.M,linear.fec)
  }) 
}else
{
  for(l in 1: N.sp)
  {
    store.species.steepness[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                                AssessYr,"/steepness/h.prior.csv",sep=''))
  }
}

#compare steepness and r
if(do.steepness)
{
  CompR=data.frame(Name=names(store.species.steepness),
                   h=unlist(sapply(store.species.steepness, `[`, 3)),
                   h.sd=unlist(sapply(store.species.steepness, `[`, 4)),
                   r=unlist(sapply(store.species.r, `[`, 3)),
                   r.sd=unlist(sapply(store.species.r, `[`, 4)))
  rownames(CompR)=NULL
  COL=rgb(.1,.2,.8,alpha=.45)
  CompR=CompR%>%
    arrange(r)%>%
    mutate(Name=capitalize(Name))
  my_formula=y ~ x
  Omit.these=c("Great hammerhead","Scalloped hammerhead","Smooth hammerhead","Tiger shark")
  p=CompR%>%
    ggplot(aes(r, h, label = Name)) +
    geom_point(shape = 21, size = 5,fill=COL) + 
    geom_smooth(method = "lm", data = CompR%>%filter(!Name%in%Omit.these),
                se = F, fullrange = TRUE,colour="red")+
    stat_poly_eq(data = CompR%>%filter(!Name%in%Omit.these),
                 aes(label =  paste(stat(eq.label),stat(adj.rr.label),stat(p.value.label),
                                    sep = "*\", \"*")),
                 formula = my_formula, parse = TRUE,
                 label.y = "bottom", label.x = "right", size = 4)+
    geom_text_repel(segment.colour='black',col='black',box.padding = 0.5) + 
    scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))+ 
    theme_PA(axs.T.siz=14,axs.t.siz=12)+
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ylim(0,1)
  #  geom_errorbar(aes(ymin=h-h.sd, ymax=h+h.sd),colour=COL)+
  #  geom_errorbarh(aes(xmin=r-r.sd, xmax=r+r.sd),colour=COL)
  p
  ggsave(handl_OneDrive('Analyses/Population dynamics/Steepness_vs_r.tiff'), 
         width = 8,height = 6, dpi = 300, compression = "lzw")
  
}


#Recalculate  steepness for hammerheads and tiger using linear model of r on h 
store.species.steepness.S2=fn.get.stuff.from.list(store.species.steepness,"mean")
dis.sp.h=c("great hammerhead","scalloped hammerhead","smooth hammerhead","tiger shark")

for(s in 1:length(dis.sp.h))
{
  id=match(dis.sp.h[s],names(store.species.steepness))
  store.species.steepness.S2[[id]]=store.species.r[[id]]$mean*0.856+0.197
}
  


#---15. Size-based Catch curve with specified selectivity--------------------------------------
#note: derive F from catch curve and gear selectivity
#      assume start of year in January (coincides with birht of most species)
if(do.Size.based.Catch.curve)
{
  mm.conv=10 # total length in mm #NEW 
  fn.source("Length_based.catch.curve.R")
  fn.extract.dat=function(STRING,Files) grep(paste(STRING,collapse="|"), Files, value=TRUE)
  
  
    #11.1 TDGDLF #ACA
  size.catch.curve_TDGDLF=vector('list',N.sp)
  names(size.catch.curve_TDGDLF)=Keep.species
  system.time({for(l in 1: N.sp)  
  {
    # Get F
    Outfile='TDGDLF' 
    print(paste("Size-base catch curve for --",names(Species.data)[l],"---",Outfile))
    
    this.size.comp=paste('Size_composition',c('West.6.5','West.7','Zone1.6.5','Zone1.7','Zone2.6.5','Zone2.7'),sep="_")
    outfile=paste(Outfile,'_histogram',sep='')
    iid=Species.data[[l]][fn.extract.dat(this.size.comp,names(Species.data[[l]]))]
    if(length(iid)>0)
    {
      dummy=do.call(rbind,iid)
      if(grepl("TDGDLF",outfile))
      {
        dummy=dummy%>%
          mutate(dummy=sub(".*Size_composition_", "", rownames(dummy)),
                 Mesh=word(dummy,2,sep = "\\."),
                 Mesh=ifelse(Mesh=='6','6.5',Mesh),
                 Mesh=factor(Mesh,levels=c('6.5','7')),
                 Zone=word(dummy,1,sep = "\\."))
        
      }
      dummy=dummy%>%filter(year<=as.numeric(substr(Last.yr.ktch,1,4)))
      
      N.min=dummy%>%
        group_by(year)%>%
        tally()%>%
        filter(n>=Min.annual.obs)%>%
        mutate(Keep=year)
      if(nrow(N.min)>0)
      {
        dummy=dummy%>%
          mutate(Keep=year)%>%
          filter(Keep%in%N.min$Keep)%>%
          mutate(TL=mm.conv*FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)     
        
        #1. Plot observed size frequency by year and mesh for years with minimum sample size
        if(grepl("TDGDLF",outfile))
        {
          p=dummy%>%
            ggplot( aes(x=TL/mm.conv, color=Mesh, fill=Mesh)) +
            geom_histogram(alpha=0.6, binwidth = TL.bins.cm)
          WHERE="top"
        }else
        {
          p=dummy%>%
            ggplot( aes(x=TL/mm.conv,color=year, fill=year)) +
            geom_histogram(alpha=0.6, binwidth = TL.bins.cm)
          WHERE="none"
        }
        p=p+
          facet_wrap(~year,scales='free_y')+
          xlab("Total length (cm)")+ylab("Count")+
          theme(legend.position=WHERE,
                legend.title = element_blank(),
                legend.text=element_text(size=14),
                strip.text.x = element_text(size = 12),
                axis.text=element_text(size=12),
                axis.title=element_text(size=16))
        print(p)
        ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(List.sp[[l]]$Name),"/",AssessYr,
                     "/1_Inputs/Visualise data/Size.comp.",outfile,".tiff",sep=''),
               width = 12,height = 10,compression = "lzw")
        
        #2. Calculate F
        #initial value of F
        params = log(0.2)   
        
        # age
        MaxAge = ceiling(mean(List.sp[[l]]$Max.age.F))
        
        # ages
        Ages = 0:MaxAge
        
        # Growth parameters
        Linf = mm.conv*with(List.sp[[l]],Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)  #need total length  
        vbK = List.sp[[l]]$Growth.F$k
        Lo =  mm.conv*with(List.sp[[l]],Lo*a_FL.to.TL+b_FL.to.TL)
        #tzero = 0              
        CVLenAtAge = 0.1      #assumption (MISSING... mention in Methods!!!)
        
        
        # length structure of model
        MaxLen = mm.conv*10*round(List.sp[[l]]$TLmax/10)
        LenInc = mm.conv*TL.bins.cm
        
        # length bins
        min.Lo=min(Lo-Lo*CVLenAtAge,min(dummy$TL))
        lbnd = seq(LenInc*floor((min.Lo)/LenInc),MaxLen - LenInc, LenInc)
        ubnd = lbnd + LenInc
        midpt = lbnd + (LenInc/2)
        nLenCl = length(midpt)
        
        # natural mortality
        PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(List.sp[[l]]$Name),"/",AssessYr,"/demography",sep='')
        
        NatMort =read.csv(paste(PATH,'M.csv',sep='/')) # from demography
        NatMort= mean(colMeans(NatMort,na.rm=T))
        
        # gillnet selectivity 
        
        #SelAtLength=read.csv('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/out.sel.csv')
        #SelAtLength$TL=midpt
        #SelAtLength$Sel.combined=SelAtLength$x
        #SelAtLength=SelAtLength$Sel.combined
       
         SelAtLength=Selectivity.at.totalength[[l]]%>%               
                      mutate(TL=TL*mm.conv)%>%
                      filter( TL%in%midpt)%>%
                     pull(Sel.combined)
        
        #Execute model functions           
        MeanSizeAtAge = CalcMeanSizeAtAge(Lo,Linf, vbK)
        RecLenDist = CalcSizeDistOfRecruits(GrowthCurveResults=MeanSizeAtAge, CVLenAtAge)
        
        LTM = CalcLTM(Linf, vbK, CVLenAtAge, midpt)   
        
        #Fit model for each year with length data
        nyrs=dummy%>%
          group_by(year)%>%
          tally()%>%pull(year)
        add.dummy=data.frame(bin=midpt)
        
        cl <- makeCluster(detectCores()-1)
        registerDoParallel(cl)
        clusterEvalQ(cl, .libPaths('C:/~/R/win-library/4.0'))  #added bit to point where doparallel is 
        F.at.year=foreach(q=1:length(nyrs),.packages=c('dplyr','doParallel'),.export=c('ObjectiveFunc')) %dopar%
        {
            x=dummy%>%
              filter(year==nyrs[q])%>%
              mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2)%>%
              group_by(bin)%>%
              tally()%>%
              full_join(add.dummy,by='bin')%>%
              arrange(bin)%>%
              mutate(n=ifelse(is.na(n),0,n))%>%
              filter(bin%in%midpt)
            ObsCatchFreqAtLen=x%>%pull(n)
            nlmb <- nlminb(params, ObjectiveFunc, gradient = NULL, hessian = TRUE)
            
            
            #Calculate uncertianty for parameter estimates
            #note: get variance-covariance matrix from fitted model
            hess.out = optimHess(nlmb$par, ObjectiveFunc)
            vcov.Params = solve(hess.out)
            ses = sqrt(diag(vcov.Params)) # asymptotic standard errors of parameter estimates
            EstFMort = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
            Table.check=data.frame(year=nyrs[q],
                                   objective.fun=nlmb$objective,
                                   convergence=nlmb$convergence,
                                   par.low95=EstFMort[2],
                                   par=EstFMort[1],
                                   par.up95=EstFMort[3])   
            
            #estimated
            ExpCatchAtLen = GetExpCatchAtLen(nlmb$par)
            
            return(list(Table=Table.check,
                        Obs=ObsCatchFreqAtLen/sum(ObsCatchFreqAtLen),
                        Pred=unlist(ExpCatchAtLen),
                        Size=midpt))
            rm(ExpCatchAtLen,ObsCatchFreqAtLen)
        }
        stopCluster(cl)
        names(F.at.year)=nyrs
        size.catch.curve_TDGDLF[[l]]=F.at.year
        rm(dummy,SelAtLength,NatMort,midpt,Ages,MeanSizeAtAge,RecLenDist,LTM,nyrs)
        
        
        #Export F  
        PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(List.sp[[l]]$Name),"/",AssessYr,"/Size_based.Catch.curve",sep='')
        if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))   
        out=do.call(rbind,subListExtract(size.catch.curve_TDGDLF[[l]],"Table"))%>%
          rename(Low95=par.low95,
                 Mean=par,
                 Up95=par.up95)
        write.csv(out,paste(PATH,paste(unlist(strsplit(Outfile, split='_', fixed=TRUE))[1],"F.csv",sep="."),sep='/'),row.names = F) 
        
      }  #end if  nrow(N.min)>0 statement
    }  #end if length(iid)>0 statement
  }  # end l  
  })    #takes 10 mins
  
      #obs vs pred  
  for(l in 1: N.sp)
  {
    dummy=size.catch.curve_TDGDLF[[l]]
    if(!is.null(dummy))
    {
      PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[l]]$Name),"/",AssessYr,"/Size_based.Catch.curve",sep='')
      
      fn.fig(paste(PATH,"/TDGDLF_F.fit",sep=""),2400,2400)
      smart.par(n.plots=length(dummy),MAR=c(1.2,2,2.5,1.25),OMA=c(2,2,.2,2.1),MGP=c(1,.5,0))
      for(i in 1:length(dummy))
      {
        with(dummy[[i]],{
          plot(Size/mm.conv, Obs, xlab="", ylab="",pch=19, col="orange",cex=1.25)
          lines(Size/mm.conv, Pred,lwd=1.25)
          mtext(Table$year, cex=1.25,3)
        })
      }
      legend("topright",c("Observed","Predicted"),pch=c(19,NA),lty=c(NA,1),
             col=c("orange","black"),bty='n',cex=1.25)
      mtext("Total length (cm)",1,outer=T,cex=1.25,line=.5)
      mtext("Probability",2,outer=T,las=3,cex=1.25,line=.5)
      dev.off()

    }
    rm(dummy)
  }
  
    #11.2 NSF
  #note: unknown selectivity for longline & only enough data for sandbar
  
    #11.3 Pilbara trawl
  #note: not enough observations

}else
{
  size.catch.curve_TDGDLF=vector('list',N.sp)
  names(size.catch.curve_TDGDLF)=Keep.species
  for(l in 1: N.sp)
  {
    this.file=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(List.sp[[l]]$Name),"/",AssessYr,
                    "/Size_based.Catch.curve/TDGDLF.F.csv",sep='')
    if(file.exists(this.file)) size.catch.curve_TDGDLF[[l]]=read.csv(this.file)
    rm(this.file)
  }
}





#---POPULATION DYNAMICS. 'Other' species-------------------------------------------------

if(do.other.ass)
{
  N.sp.other=N.sp-length(Indicator.species)
  Keep.species.other=subset(Keep.species,!Keep.species%in%names(Indicator.species))  
    
  # Bring in cpue  
  cpue.list=vector('list',N.sp.other)
  names(cpue.list)=Keep.species.other
  for(l in 1: N.sp.other)
  {
    id=match(Keep.species.other[l],names(List.sp))
    
    Nm=str_remove(List.sp[[id]]$Name.inputs, ' shark')
    temp.wd=paste(HandL,capitalize(Nm),'/',AssessYr,sep='')
    Files=list.files(temp.wd)
    iid=Files[fn.extract.dat(STRING="cpue.annual",Files)]
    if(length(iid)>0)
    {
      dummy=list()
      if("cpue.annual.survey.csv"%in%iid) dummy$Survey=read.csv(paste(temp.wd,'/cpue.annual.survey.csv',sep=''))
      if("cpue.annual.TDGDLF.csv"%in%iid) dummy$TDGDLF.mon=read.csv(paste(temp.wd,'/cpue.annual.TDGDLF.csv',sep=''))
      if("cpue.annual.TDGDLF.daily.csv"%in%iid) dummy$TDGDLF.day=read.csv(paste(temp.wd,'/cpue.annual.TDGDLF.daily.csv',sep=''))
      if("cpue.annual.NSF.csv"%in%iid) dummy$NSF=read.csv(paste(temp.wd,'/cpue.annual.NSF.csv',sep=''))
      cpue.list[[l]]=dummy
      rm(dummy)
    }
  }
  
  # Bring in catch  
  Tot.ktch=KtCh
  
  #---Single-species SPM -----------------------------------------------------------------------
  if(Do.SPM=="YES")                        
  {
    #Initial value of q
    Estimable.qs=vector('list',N.sp.other)
    names(Estimable.qs)=Keep.species.other
    for(l in 1: N.sp.other)
    {
      if(!is.null(cpue.list[[l]]))
      {
        q1=q2=q3=q4=NA
        if("Survey"%in%names(cpue.list[[l]])) q1=.005
        if("TDGDLF.mon"%in%names(cpue.list[[l]])) q2=.005
        if("TDGDLF.day"%in%names(cpue.list[[l]])) q3=.001
        if("NSF"%in%names(cpue.list[[l]])) q4=.005
        Estimable.qs[[l]]=c(q1=q1,q2=q2,q3=q3,q4=q4)
      }
    }
    
    #Penalty for keeping biomass positive
    posfun=function(x,eps,pen)  
    {
      if (x>=eps) return(x) else
      {
        pen=pen+.01*(x-eps)^2
        return (list(eps/(2-x/eps),pen))
      }
    }
    
    #Surplus production model
    SPM=function(Init.propK,cpue,cpue.CV,Qs,Ktch,theta,HR_init,HR_init.sd,r.mean,r.sd,usePenalties) 
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
    
    #fill in missing years
    fn.fill=function(x)   
    {
      aa=all.iers[which(!all.iers%in%x$yr)]
      aa1=x[1:length(aa),]
      aa1[,]=NA
      aa1$yr=aa
      x=rbind(x,aa1)%>%arrange(yr)
      return(x)
    }
    
    #Check max possible initial harvest rate  
    Mx.init.harv=rep(NA,N.sp.other)
    names(Mx.init.harv)=Keep.species.other
    for(s in 1: N.sp.other)
    {
      id=match(Keep.species.other[s],names(List.sp))
      if(List.sp[[id]]$B.init==1)   
      {
        Mx.init.harv[s]=0
      }else
      {
        ct=Tot.ktch%>%filter(Name==names(Mx.init.harv)[s])%>%
          group_by(finyear)%>%
          summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
        Mx.init.harv[s]=ct$LIVEWT.c[1]/max(ct$LIVEWT.c)
      }
    }
    
    #Estimate parameters                              
    Store.SPM=vector('list',N.sp.other)
    names(Store.SPM)=Keep.species.other
    Store.stuff=Theta=Store.SPM
    for(s in 1: N.sp.other)
    {
      #catch
      ct=Tot.ktch%>%filter(Name==names(Mx.init.harv)[s])%>%
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
      Idd=match(Keep.species.other[s],names(store.species.r))
      r.prior=store.species.r[[Idd]]$mean
      r.prior.sd=store.species.r[[Idd]]$sd
      
      #run if cpue available
      Store.CPUE.eff.dummy=NULL
      QS.dummy=NULL
      CPUE.yr.dummy=NULL
      n.cpue.dummy=NULL
      CPUE.CV=NULL
      n.cpues=NULL
      if(!is.null(cpue.list[[s]]))
      {
        Idd=match(Keep.species.other[s],names(List.sp))
        MAX.CV=List.sp[[Idd]]$MAX.CV
          
        #Get cpues
        CPUE.1=cpue.list[[s]]$Survey
        if(!is.null(CPUE.1))CPUE.1=CPUE.1%>%
            filter(CV<MAX.CV)
        CPUE.2=cpue.list[[s]]$TDGDLF.mon
        if(!is.null(CPUE.2))CPUE.2=CPUE.2%>%
          mutate(yr=as.numeric(substr(Finyear,1,4)),
                 MeAn=Mean)%>%
          filter(CV<MAX.CV)
        CPUE.3=cpue.list[[s]]$TDGDLF.day
        if(!is.null(CPUE.3))CPUE.3=CPUE.3%>%
          mutate(yr=as.numeric(substr(Finyear,1,4)),
                 MeAn=Mean)%>%
          filter(CV<MAX.CV)
        CPUE.4=cpue.list[[s]]$NSF                    
        if(!is.null(CPUE.4))CPUE.4=CPUE.4%>%
          mutate(yr=as.numeric(substr(FINYEAR,1,4)),
                 MeAn=Mean)%>%
          filter(CV<MAX.CV)
        
        n.cpues=length(cpue.list[[s]])
        
        #fill in missing cpue
        CPUE=list(CPUE.1,CPUE.2,CPUE.3,CPUE.4)
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
        K.init=List.sp[[Idd]]$k.times.mx.ktch*Mx.ktch
        r.init=List.sp[[Idd]]$Init.r
        QS=Estimable.qs[[s]]
        
        #loop over scenarios
        HR.o.scens=Mx.init.harv[s]
        dummy=vector('list',length(HR.o.scens))
        names(dummy)=HR.o.scens
        
        
        Efficien.scens=List.sp[[Idd]]$Efficien.scens
        
        for(h in 1:length(HR.o.scens))
        {
          #. initial harvest rate prior
          HR_o=HR.o.scens[h]  
          
          dummy.eff=vector('list',length(Efficien.scens))
          names(dummy.eff)=Efficien.scens
          Store.CPUE.eff=dummy.eff
          Eff.yrs=ct$finyear
          id.eff.yrs=which(Eff.yrs>1994)
          for(e in 1:length(Efficien.scens))
          {
            #. apply assumed efficiency to TDGDLF
            Add.eff=data.frame(yr=Eff.yrs,Efficiency=1)
            Add.eff$Efficiency[id.eff.yrs]=Add.eff$Efficiency[id.eff.yrs]-
              cumsum(rep(Efficien.scens[e],length(id.eff.yrs)))
            CPUE.eff.scen=CPUE
            for(ss in 2:3)
            {
              if(!is.null(CPUE.eff.scen[[ss]])) CPUE.eff.scen[[ss]]=CPUE.eff.scen[[ss]]*Add.eff$Efficiency
            }
            
            #. objfun to minimize
            fn_ob=function(theta)SPM(Init.propK=List.sp[[Idd]]$B.init,
                                     cpue=CPUE.eff.scen,
                                     cpue.CV=CPUE.CV,
                                     Qs=QS,
                                     Ktch=ct$LIVEWT.c,
                                     theta,
                                     HR_init=HR_o,
                                     HR_init.sd=List.sp[[Idd]]$HR_o.sd,
                                     r.mean=r.prior,
                                     r.sd=r.prior.sd,
                                     usePenalties=List.sp[[Idd]]$usePen)$negLL
            #. fit model
            estim.q=List.sp[[Idd]]$estim.q
            if(estim.q=="YES")theta= c(k=log(K.init),log(QS[which(!is.na(QS))]))
            if(estim.q=="NO")theta= c(k=log(K.init))
            Lw.bound=log(c(List.sp[[Idd]]$Low.bound.K*Mx.ktch,rep(1e-6,length(theta)-1)))
            Up.bound=log(c(List.sp[[Idd]]$Up.bound.K*Mx.ktch,rep(1,length(theta)-1)))
            if(!List.sp[[Idd]]$fix.r=="YES")
            {
              theta=c(theta,r=log(r.init))
              Lw.bound=c(Lw.bound,log(0.01))
              Up.bound=c(Up.bound,log(0.75))
            }
            
            if(List.sp[[Idd]]$what.like=='kernel')
            {
              theta=c(theta,sigma=log(0.2))
              Lw.bound=c(Lw.bound,log(1e-2))
              Up.bound=c(Up.bound,log(1))
            }
            
            #1st estimation round
            theta=nlminb(theta, fn_ob, gradient = NULL,lower =Lw.bound,upper = Up.bound)$par
            
            #2nd estimation round
            if(List.sp[[Idd]]$minimizer=='optim')
            {
              OptiM=optim(theta,fn_ob,method="L-BFGS-B",lower =Lw.bound,upper = Up.bound,hessian=T)
              if(List.sp[[Idd]]$Remove.bounds)
              {
                paramscale = magnitude(theta)
                OptiM=optim(jitter(OptiM$par),fn_ob,method="L-BFGS-B",hessian=T,
                            control = list(maxit = 1000, parscale = paramscale))
              }
              OptiM$hit.upper.boundary=names(which(round(1-exp(OptiM$par)/exp(Up.bound),1)==0))
              OptiM$hit.lower.boundary=names(which(round(1-exp(OptiM$par)/exp(Lw.bound),1)==0))
              dummy.eff[[e]]=OptiM
            }
            if(List.sp[[Idd]]$minimizer=='nlminb')
            {
              nlmb <- nlminb(theta, fn_ob, gradient = NULL,lower =Lw.bound,upper = Up.bound)
              if(List.sp[[Idd]]$Remove.bounds) nlmb <- nlminb(jitter(nlmb$par), fn_ob, gradient = NULL)
              nlmb$hit.upper.boundary=names(which(round(1-exp(nlmb$par)/exp(Up.bound),1)==0))
              nlmb$hit.lower.boundary=names(which(round(1-exp(nlmb$par)/exp(Lw.bound),1)==0))
              
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
                            cpue.yrs=CPUE.yr.dummy,n.cpues=n.cpues)
      
    }
    
    #Evaluate model at MLE
    SPM.preds=vector('list',length(Store.SPM))
    names(SPM.preds)=names(Store.SPM)
    for(s in 1: N.sp.other)
    {
      if(!is.null(Store.SPM[[s]]))
      {
        Idd=match(Keep.species.other[s],names(List.sp))
        
        HR.o.scens=Mx.init.harv[s]
        dumy.pred=vector('list',length(HR.o.scens))
        names(dumy.pred)=paste("HR=",HR.o.scens)
        for(h in 1:length(HR.o.scens))
        {
          dummy.eff=vector('list',length(Efficien.scens))
          names(dummy.eff)=Efficien.scens
          for(e in 1:length(Efficien.scens))
          {
            dummy.eff[[e]]=SPM(Init.propK=List.sp[[Idd]]$B.init,
                               cpue=Store.stuff[[s]]$cpue[[e]],
                               cpue.CV=Store.stuff[[s]]$CPUE.CV,
                               Qs=Store.stuff[[s]]$Qs,
                               Ktch=Store.stuff[[s]]$Ktch,
                               theta=Store.SPM[[s]][[h]][[e]]$par,
                               HR_init=HR.o.scens[h],
                               HR_init.sd=List.sp[[Idd]]$HR_o.sd,
                               r.mean=Store.stuff[[s]]$r.mean,
                               r.sd=Store.stuff[[s]]$r.sd,
                               usePenalties=List.sp[[Idd]]$usePen)
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
    for(s in 1: N.sp.other)
    {
      if(!is.null(SPM.preds[[s]]))
      {
        Idd=match(Keep.species.other[s],names(List.sp))
        
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
              dum[[x]]=SPM(Init.propK=List.sp[[Idd]]$B.init,
                           cpue=Store.stuff[[s]]$cpue[[e]],
                           cpue.CV=Store.stuff[[s]]$CPUE.CV,
                           Qs=Store.stuff[[s]]$Qs,
                           Ktch=Store.stuff[[s]]$Ktch,
                           theta=Estim.par.samples[[s]][[h]][[e]][x,],
                           HR_init=HR.o.scens[h],
                           HR_init.sd=List.sp[[Idd]]$HR_o.sd,
                           r.mean=Store.stuff[[s]]$r.mean,
                           r.sd=Store.stuff[[s]]$r.sd,
                           usePenalties=List.sp[[Idd]]$usePen)
            }
            dummy.eff[[e]]=dum
          }
          dumy.pred[[h]]=dummy.eff
        }
        SPM.preds_uncertainty[[s]]=dumy.pred
        rm(HR.o.scens)
        print(paste("SPM  ",s,"--",names(Store.SPM)[s]))
      }
    }
    
  }
  
  #---Catch-MSY ---------------------------------------------------------------
  #note: uses SPM inputs (already done r prior, ct, etc)
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
      hndl.sim.test=handl_OneDrive('Analyses\\Population dynamics\\1.Other species\\2019\\Outputs\\Simulation_testing_catchMSY\\')
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
    
    # - Run CMSY for each species    takes 0.03 sec per iteration per species 
    Store.CMSY=vector('list',N.sp.other)
    names(Store.CMSY)=Keep.species.other
    
    system.time(for(s in 1: N.sp.other)
    {
      Idd=match(Keep.species.other[s],names(List.sp))
      
      #Tested scenarios
      SCENARIOS=List.sp[[Idd]]$SCENARIOS
      ktch_msy_scen=vector('list',length(SCENARIOS))
      names(ktch_msy_scen)=names(SCENARIOS)
      for(sc in 1:length(ktch_msy_scen))
      {
        if(is.na(SCENARIOS[[sc]]$R.prior[1]))    #ACA, fix this, the SCENARIOS list is not right....
        {
          USR="No"
          ReS=RESILIENCE[[Idd]]
        }else
        {
          if(!is.na(SCENARIOS[[sc]]$R.prior[1]))USR= "Yes"
          ReS=NA
          Scen.start.bio=SCENARIOS[[sc]]$Initial.dep   
        }
        
        ktch_msy_scen[[sc]]=list(r.prior=SCENARIOS[[sc]]$R.prior,
                                 user=USR,
                                 k.lower=List.sp[[Idd]]$Low.bound.K,
                                 k.upper=List.sp[[Idd]]$Up.bound.K, 
                                 startbio=Scen.start.bio,
                                 finalbio=List.sp[[Idd]]$FINALBIO,
                                 res=ReS,
                                 niter=List.sp[[Idd]]$SIMS,
                                 sigR=SCENARIOS[[sc]]$Error)
        if(!is.na(ktch_msy_scen[[sc]]$r.prior[1])) ktch_msy_scen[[sc]]$r.prior=store.species.r[[Idd]]$mean
        rm(USR,ReS)
      }
      
      #Execute Catch-MSY  
      print(paste("CMSY------","s=",s,names(List.sp)[Idd]))
      
      PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[Idd]]$Name),"/",AssessYr,"/CMSY",sep='')
      if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))   
      setwd(file.path(PATH))
      PHat=file.path(PATH)
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
      Store.CMSY[[s]]$KTCH.MSY=Ktch_MSY
      Store.CMSY[[s]]$Catch=ct
      Store.CMSY[[s]]$K=Ktch_MSY[[sc]]$k
      
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
    #Tweaked functions to allow multiple CPUEs and Qs, and total biomass and Hessian calculation for uncertainty
    TotB=function (invect, WeightA) 
    {
      ans <- sum(WeightA * invect)/1000
      return(ans)
    }
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
                   "PredCE","PredCE2", "Deplete","TotalB")
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
        totb <- TotB(Nt[, (yr - 1)], waa)
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
        fishery[(yr - 1), 12] <- totb
      }
      totb <- TotB(Nt[, yr], waa)
      spb <- SpB(Nt[, yr], maa, waa)
      exb <- ExB(Nt[, yr] * hS, sela, waa)
      fishery[yr, 4:5] <- c(spb, exb)
      fishery[yr, 12] <- totb
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
        #id=match(q1.yrs,infish$year)
        id=match(infish$year[1]:2005,infish$year)
        avq1 <- exp(mean(log(infish$cpue[id]/fishery[id+1, "ExploitB"]),na.rm = TRUE))
        fishery[id+1, "PredCE"] <- ExpB[id] * avq1
        
        #id=match(q2.yrs,infish$year)
        id=match(2006:infish$year[length(infish$year)],infish$year)
        avq2 <- exp(mean(log(infish$cpue[id]/fishery[id+1, "ExploitB"]),na.rm = TRUE))
        fishery[id+1, "PredCE"] <- ExpB[id] * avq2
      }
      
      #fleet 2
      avq.fleet2 <- exp(mean(log(infish$cpue2/fishery[1:nyrs, "ExploitB"]), 
                             na.rm = TRUE))
      fishery[2:(nyrs + 1), "PredCE2"] <- ExpB * avq.fleet2
      
      #Add CVs for fit display
      fishery=cbind(fishery,rbind(c(NA,NA),infish[,c('se','se2')]))
      
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
    
    
    Store.age.comp=vector('list',N.sp.other)
    names(Store.age.comp)=Keep.species.other
    
    #fill in objects required for aSPM
    for(l in 1:length(Store.age.comp))
    {
      Idd=match(Keep.species.other[l],names(List.sp))  
      AMAX=List.sp[[Idd]]$Max.age.F[1]
      Lo=List.sp[[Idd]]$Lzero
      Linf=List.sp[[Idd]]$Growth.F$FL_inf
      k=List.sp[[Idd]]$Growth.F$k
      age=0:AMAX 

      Store.age.comp[[l]]$glb=list(
            maxage=AMAX,
            M=mean(unlist(lapply(store.species.r[[Idd]]$M, mean))),  
            Linf=Linf,
            K=k,
            t0=NA,                         
            Waa=NA,
            Wab=NA,
            M50a=floor(mean(List.sp[[Idd]]$Age.50.mat)),
            deltaM=NA,
            steep=store.species.steepness[[Idd]]$mean,
            R0=NA,
            sela50=NA,
            deltaS=NA,
            resilience=RESILIENCE[[Idd]],
            nages=length(floor(age)),
            ages=floor(age),
            nyrs=NA,
            spsname=names(RESILIENCE)[Idd])
      
      mn.len=Lo+(Linf-Lo)*(1-exp(-k*(age)))
      
      Store.age.comp[[l]]$props=data.frame(
            age=floor(age),
            laa=mn.len,
            waa= List.sp[[Idd]]$AwT*mn.len^List.sp[[Idd]]$BwT,  #catch in tonnes; waa in kgs 
            maa=plogis(floor(age),floor(mean(List.sp[[Idd]]$Age.50.mat)),1),
            sela=Selectivity.at.age[[Idd]]$Sel.combined[1:length(age)],
            feca=NA) 
      
    }

 
    #init par values
    aSPM.init=vector('list',N.sp.other)
    names(aSPM.init)=Keep.species.other
    for(l in 1: N.sp.other)
    {
      if(!is.null(cpue.list[[l]]))
      {
        logR0=10
        sigCE=0.3
        if(names(aSPM.init)[l]=="tiger shark") logR0=11
        aSPM.init[[l]]=c(logR0=logR0,sigCE=sigCE)
      }
    }

    #run aspm
    Store.aSPM=vector('list',N.sp.other)
    names(Store.aSPM)=names(aSPM.init)
    No.signal.cpue=c("milk shark","spinner shark")  #no signal in cpue
    for(l in 1:N.sp.other)    
    {
      print(paste("aSPM--------- l=",l,names(Store.aSPM)[l]))
      if(!is.null(cpue.list[[l]])&!names(Store.aSPM)[l]%in%No.signal.cpue)
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
          if("Survey"%in%names(CPUE))
          {
            a=CPUE[-match("Survey",names(CPUE))]
            
            CPUE2=CPUE[[match("Survey",names(CPUE))]]%>%    
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
        if(min(CPUE$se,na.rm=T)>1)  CPUE$se=CPUE$se/100  #reset if CV in percentage
        
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
        if(exists('CPUE')) rm(CPUE)
        if(exists('CPUE2')) rm(CPUE2)
        
        Store.aSPM[[l]]=list(fit=ans,quantities=fishery,fish=fish,quantities.MC=fishery.MC)
      }
    }
  }
  
  
  
  
}

#---POPULATION DYNAMICS. Indicator species-------------------------------------------------
#code from Assessment.R code
#note: 'Run.models.R' outputs data and parameter inputs for models,
#                     runs assessment models
#                     and displays outputs
#reset dummies
Spec=1
Pin.pars=1  #dummy to clear log
Par.phases=1
all.objects=objects()
List.objs=unique(unlist(lapply(List.sp,names)))
suppressWarnings(rm(list=all.objects[which(all.objects%in%List.objs)]))
for(l in 1:length(List.sp))
{
  attach(List.sp[[l]])
  source(handl_OneDrive("Analyses/Population dynamics/Git_Stock.assessments/Run.models.R"))
  detach(List.sp[[l]])
}




#---RESULTS.'Other' species ------------------------------------------------- 
if(do.other.ass)
{
  if(do.other.ass.paper) hNdl=handl_OneDrive('Analyses/Population dynamics/Other species/2020') else
    hNdl=paste(handl_OneDrive('Analyses/Population dynamics/Other species/'),AssessYr,sep='')
  if(!dir.exists(hNdl))dir.create(hNdl)
  
  #---RESULTS. Catches of all species by fishery ----
  Tot.ktch=KtCh %>%   
    mutate(
      Type = case_when(
        FishCubeCode=='WRL'~'WRL',
        FishCubeCode=='WTB'~'WTB',
        FishCubeCode=='TEP'~'TEP',
        FishCubeCode=='Recreational'~'Recreational',
        FishCubeCode=='SA MSF'~'SA MSF',
        FishCubeCode=='GAB'~'GAB',
        FishCubeCode=='Indo'~'Indonesia',
        FishCubeCode=='Taiwan'~'Taiwan',
        FishCubeCode=='Historic'~'Commercial',
        FishCubeCode%in%c('JASDGDL','WCDGDL','C070','OAWC')~'Commercial',
        FishCubeCode%in%c('JANS','OANCGC','WANCS')~'Commercial',   
        TRUE  ~ "Commercial"))
  all.yrs=min(Tot.ktch$finyear):max(KtCh$finyear)
  Fishry.type=sort(unique(Tot.ktch$Type))
  COLs.type=colfunc(length(Fishry.type))
  names(COLs.type)=Fishry.type
  All.N.sp=sort(unique(Tot.ktch$Name))
  All.N.sp=subset(All.N.sp,!All.N.sp%in%names(Indicator.species))
  
  fn.fig(paste(hNdl,'/Outputs/Figure 1_catch_all_species',sep=''),2400,2400) 
  smart.par(n.plots=length(All.N.sp),MAR=c(1,1,.8,.25),OMA=c(2.5,2.25,.05,.05),MGP=c(1,.5,0))
  par(cex.main=1,cex.axis=.85)
  for(s in 1:length(All.N.sp))
  {
    ddd=subset(Tot.ktch,Name==All.N.sp[s])%>%
      group_by(Type,finyear)%>%
      summarise(Tot=sum(LIVEWT.c,na.rm=T)/unitS)%>%
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
    if(s%in%31:36)axis(1,Yrss,Yrss,tck=-.07)
    Yax=pretty(seq(0,max(ddd$Tot),length.out = 4))
    axis(2,Yax,Yax, cex.axis=.85,padj =0.75,tck=.07,las=3)
    if(s==1)legend('topleft',names(COLs.type)[1:4],pt.bg=COLs.type[1:4],bty='n',pch=21,cex=1.05,pt.cex=1.5)
    if(s==7) legend('topleft',names(COLs.type)[5:length(COLs.type)],pt.bg=COLs.type[5:length(COLs.type)],
                    bty='n',pch=21,cex=1.05,pt.cex=1.5)
    
  }
  mtext("Financial year",1,line=1,cex=1.5,outer=T)
  mtext("Total catch (tonnes)",2,las=3,line=0.5,cex=1.5,outer=T)
  dev.off()
  
  
  #---RESULTS. Total catch and cpue together ----
  fn.fig(paste(hNdl,'/Outputs/Figure_Catch and cpue',sep=''),2400,1800) 
  smart.par(n.plots=N.sp,MAR=c(2,2,1,1),OMA=c(1.75,2,.5,1),MGP=c(1,.5,0))
  for(s in 1: N.sp)
  {
    ct=Tot.ktch%>%filter(SP.group==Specs$SP.group[s])%>%
      group_by(finyear)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
    plot(ct$finyear,ct$LIVEWT.c,type='o',pch=21,bg='orange',ylab="",xlab="",main=capitalize(Specs$SP.group[s]))
    
    if(!is.null(cpue.list[[s]]))
    {
      Survey=cpue.list[[s]]$Survey
      TDGDLF.mon=cpue.list[[s]]$TDGDLF.mon
      TDGDLF.day=cpue.list[[s]]$TDGDLF.day
      Survey.MAX=0
      TDGDLF.mon.MAX=0
      TDGDLF.day.MAX=0
      if(!is.null(Survey)) Survey.MAX=max(Survey$MeAn,na.rm=T)
      if(!is.null(TDGDLF.mon)) TDGDLF.mon.MAX=max(TDGDLF.mon$Mean,na.rm=T)
      if(!is.null(TDGDLF.day)) TDGDLF.day.MAX=max(TDGDLF.day$Mean,na.rm=T)
      Ylim=max(c(Survey.MAX,TDGDLF.mon.MAX,TDGDLF.day.MAX))
      
      par(new=T)
      plot(ct$finyear,ct$finyear,ylab='',xlab='',ylim=c(0,Ylim),col='transparent',yaxt='n') 
      if(!is.null(Survey)) lines(Survey$yr,Survey$MeAn,lwd=2,col='steelblue')
      if(!is.null(TDGDLF.mon)) lines(as.numeric(substr(TDGDLF.mon$Finyear,1,4)),TDGDLF.mon$Mean,lwd=2,col="red")
      if(!is.null(TDGDLF.day)) lines(as.numeric(substr(TDGDLF.day$Finyear,1,4)),TDGDLF.day$Mean,lwd=2,col="forestgreen")
      axis(4,seq(0,ceiling(Ylim),length.out=5),seq(0,ceiling(Ylim),length.out=5))
      
    }
    if(s==1)legend("topleft",c("Catch","Survey","cpue.mon","cpue.day"),pch=21,pt.bg=c("orange","steelblue","red","forestgreen"),
                   bty='n',cex=1.1)
  }
  mtext('Financial year',1,outer=T,line=0,cex=1.5)
  mtext('Catch (tonnes)',2,outer=T,las=3,cex=1.5,line=0)
  mtext('Relative cpue',4,outer=T,line=0,cex=1.5,las=3)
  dev.off() 
  
  
  setwd(paste(hNdl,'/Outputs',sep=''))
  
  
  #---RESULTS. Overall catch spatial distribution ----
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
  
  
  #---RESULTS. Total catch and effort ----
  ktch.s=subset(Data.monthly,SPECIES%in%c(19000,Specs$SPECIES))%>%
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
    group_by(finyear)%>%
    summarise(Tot=sum(LIVEWT.c/1000,na.rm=T))%>%    #MISSING, remove /1000, catch already in tonnes, see KTCH.UNITS
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
  
  fn.fig("Figure 4. Total effort time series", 2400, 2400)
  par(mfcol=c(1,1),mar=c(1.5,2.5,1.5,.5),oma=c(2,2,.1,4),las=1,mgp=c(1,.6,0))
  plot(Effrt.s$finyear,Effrt.s$Total,type='l',pch=19,col='grey65',cex=.75,ylab="",xlab="",lwd=3)
  mtext(side = 2, line = 2, 'Total effort (km gn days)',las=3,cex=1.75)
  
  par(new=T)
  plot(Effrt.n$finyear,Effrt.n$Hook.days,type='o',pch=19,col="black",xlab="",ylab="",axes=F,lwd=2.5)
  axis(side = 4)
  mtext(side = 4, line = 3, 'Total effort (hook days)',las=3,cex=1.75)
  
  legend("topleft",c("South","North"),bty='n',lty=1,cex=1.5,col=c("grey65","black"),lwd=3,pch=c(NA,19))
  mtext("Financial year",1, line = .5,outer=T,cex=1.75)
  dev.off()
  
  # fn.fig("Figure 4. Total Catch of analysed species and effort time series", 1800, 2400)
  # par(mfcol=c(2,1),mar=c(1.5,2.5,1.5,.5),oma=c(2,2,.1,4),las=1,mgp=c(1,.6,0))
  #   #North
  # plot(all.YYrs,ktch.n$Tot,type='o',pch=19,col='black',cex=.75,ylab="",xlab="",main="North",
  #      ylim=c(0,max(c(ktch.s$Tot,ktch.n$Tot))))
  # par(new=T)
  # plot(Effrt.n$finyear,Effrt.n$Hook.days,type='l',col="grey55",xlab="",ylab="",axes=F,lwd=2.5,lty=3)
  # axis(side = 4)
  # mtext(side = 4, line = 3, 'Total effort (hook days)',las=3,cex=1.5)
  # 
  # legend("topleft",c("Catch","Effort"),bty='n',lty=c(1,3),col=c("black","grey55"),lwd=2.5)
  # 
  #   #South
  # plot(all.YYrs,ktch.s$Tot,type='o',pch=19,col='black',cex=.75,ylab="",xlab="",main="South",
  #      ylim=c(0,max(c(ktch.s$Tot,ktch.n$Tot))))
  # par(new=T)
  # plot(Effrt.s$finyear,Effrt.s$Total,type='l',col="grey55",xlab="",ylab="",axes=F,lwd=2.5,lty=3)
  # axis(side = 4)
  # mtext(side = 4, line = 3, 'Total effort (km gn days)',las=3,cex=1.5)
  # 
  # mtext("Year",1, line = .5,outer=T,cex=1.5)
  # mtext("Total catch (tonnes)",2,outer=T,las=3,cex=1.5)
  # dev.off()
  
  
  #---RESULTS. Spatio-temporal catch ----
  #note: bubble size is proportion of blocks fished out of maximum number of blocks fished for each species
  CL=rgb(.5,.5,.5,alpha=.7)
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
    for(p in 1:length(All.sp)) points(yrs,rep(p,length(yrs)),col='black',cex=2*d1[p,],pch=19)
    axis(2,1:length(All.sp),capitalize(Sp.nms$SNAME),las=1)
    mtext(side = 1, line = 2, 'Financial year',cex=1.5)
    
    par(new=T)
    plot(yrs,Effort_blocks$Tot,type='l',col=CL,xlab="",ylab="",axes=F,lwd=5,lty=1)
    axis(side = 4,las=1)
    mtext(side = 4, line = 3, 'Number of blocks fished',las=3,cex=1.5,col=CL)
    
    DD=d1
    rownames(DD)=Sp.nms$SNAME
    return(DD)
  }
  fn.fig("Figure 2_Spatio-temporal catch", 2400, 2400)
  par(mar=c(2.5,4,.1,1),oma=c(.5,6,.1,3),mgp=c(1.5,.7,0))
  Store.spatial.temporal.ktch=fn.spatio.temp.catch.dist(d=rbind(Data.monthly%>%
                                                                  filter(!is.na(BLOCKX))%>%
                                                                  dplyr::select(SNAME,FINYEAR,SPECIES,BLOCKX),
                                                                Data.monthly.north%>%
                                                                  filter(!is.na(BLOCKX))%>%
                                                                  dplyr::select(SNAME,FINYEAR,SPECIES,BLOCKX)))
  dev.off()
  
  
  #---RESULTS. Changes in observed mean length for TDGDLF----
  Change.mean.length=vector('list',N.sp)
  names(Change.mean.length)=names(Species.data)
  fun.change.mean.len=function(d,NM,toMatch,min.annual.obs,XLIM)
  {
    d.list=d[grep(paste(toMatch,collapse="|"),names(d))]
    if(sum(grepl('Table',names(d.list)))>0)d.list=d.list[-grep('Table',names(d.list))]
    if(length(d.list)>0)
    {
      my_formula = y ~ x
      
      #by mesh
      d.list=do.call(rbind,d.list)
      d.list=d.list%>%
        mutate(mesh=case_when(grepl("6.5.inch",rownames(d.list))~6.5,
                              grepl("7.inch",rownames(d.list))~7))
      N.min=d.list%>%
        group_by(FINYEAR,mesh)%>%
        tally()%>%
        filter(n>=min.annual.obs)%>%
        mutate(Keep="YES")
      
      d.list=d.list%>%
        left_join(N.min,by=c('FINYEAR','mesh'))%>%
        filter(Keep=="YES")
      if(nrow(d.list)>0 & length(unique(d.list$FINYEAR))>2)
      {
        p=d.list%>%
          mutate(Finyear=as.numeric(substr(FINYEAR,1,4)),
                 Finyear.d=factor(Finyear,levels=sort(unique(Finyear))),
                 Mesh=paste(mesh,'inch'))%>%      
          ggplot(aes(x = Finyear, y = FL)) +
          geom_violin(aes(fill = Finyear.d, color = Finyear.d), alpha = 0.3) + 
          facet_wrap(~Mesh)+
          geom_smooth(color = "black", formula = my_formula, method = 'lm',se=TRUE)+ 
          stat_poly_eq(aes(label =  paste(stat(eq.label),stat(adj.rr.label),stat(p.value.label),
                                          sep = "*\", \"*")),
                       formula = my_formula, parse = TRUE,
                       label.y = "bottom", label.x = "right", size = 5) +
          xlab("")+ylab("")+
          theme_PA(axs.T.siz=18,axs.t.siz=14,str.siz=16)+
          theme(legend.position = "none",
                plot.title =element_text(size=20),
                panel.background = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                strip.background = element_rect(fill = "white"),
                panel.border = element_rect(colour = "black", fill=NA, size=1.25))+
          labs(title=NM)+xlim(XLIM)+ylim(0,max(d.list$FL,na.rm=T))  
        return(p)
      }
    }
  }
  for(s in 1:N.sp) 
  {
    print(paste("Change in mean length ","--",names(Species.data)[s]))
    dummy=print(fun.change.mean.len(d=Species.data[[s]],
                                    NM=capitalize(names(Species.data)[s]),
                                    toMatch=c("Size_composition_West","Size_composition_Zone1","Size_composition_Zone2"),
                                    min.annual.obs=Min.annual.obs,
                                    XLIM=c(1990,as.numeric(substr(Last.yr.ktch,1,4)))))
    if(!is.null(dummy)) Change.mean.length[[s]]=dummy
    rm(dummy)
  }
  figure=ggarrange(plotlist=fun.find.in.list(x=Change.mean.length, Drop=names(Indicator.species)),
            ncol=1)+
          theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
   annotate_figure(figure,
                   left = text_grob("Fork length (cm)", rot = 90,size=20,vjust=1),
                   bottom = text_grob("Financial year",size=20,vjust=-1))
   ggsave(paste(hNdl,'/Outputs/Figure 6_Changes in observed mean length_TDGDLF.tiff',sep=''),
         width = 12,height = 10,compression = "lzw")
  
  #---RESULTS. Changes in mean weight of individuals caught in the TDGDLF  ----
   #Any strong declining trend in mean weights? (Leitao 2019)
   Logbook=Logbook%>%
     filter(LatDeg<=(-26) & method=='GN')%>%
     mutate(species=ifelse(species==19000 & LongDeg>116,19004,species))%>%
     group_by(Same.return.SNo,species,finyear,mshigh)%>%
     summarise(nfish=sum(nfish),
               livewt=sum(livewt))%>%
     ungroup()%>%
     dplyr::select(finyear,livewt,nfish,species,mshigh)%>%
     mutate(Mean.wght=livewt/nfish)%>%
     left_join(All.species.names,by=c("species"="SPECIES"))%>%
     filter(SNAME%in%names(Species.data))%>%
     filter(!SNAME%in%names(Indicator.species))%>%
     left_join(Wei.range%>%dplyr::select(TW.min,TW.max,SPECIES),by=c("species"="SPECIES"))%>%
     filter(Mean.wght>=TW.min & Mean.wght<=TW.max)
   
   
   N.min=Logbook%>%
     group_by(finyear,species,mshigh)%>%
     tally()%>%
     filter(n>=Min.annual.obs)%>%
     mutate(Keep="YES")
   
   Logbook=Logbook%>%
     left_join(N.min,by=c('finyear','species','mshigh'))%>%
     filter(Keep=="YES")%>%
     mutate(Finyear=as.numeric(substr(finyear,1,4)),
            Finyear.d=factor(Finyear,levels=sort(unique(Finyear))))%>%
     filter(!is.na(Mean.wght))%>%
     mutate(Mesh=ifelse(mshigh==165,"6.5",ifelse(mshigh==178,'7',NA)))%>%
     filter(!is.na(Mesh))
   
   Logbook.sp=sort(unique(Logbook$SNAME))
   
   
   Change.mean.weight.catch=vector('list',length(Logbook.sp))
   names(Change.mean.weight.catch)=Logbook.sp
   fn.plt.mn.ktch.wght=function(d.list,NM,XLIM)
   {
     my_formula = y ~ x
     if(nrow(d.list)>0 & length(unique(d.list$finyear))>2)
     {
       p=d.list%>%
         mutate(Mesh=paste(Mesh,"inch"))%>%
         ggplot(aes(x = Finyear, y = Mean.wght)) +
         geom_violin(aes(fill = Finyear.d, color = Finyear.d), alpha = 0.3) + 
         geom_smooth(color = "black", formula = my_formula, method = 'lm',se=TRUE)+ 
         stat_poly_eq(aes(label =  paste(stat(eq.label),stat(adj.rr.label),stat(p.value.label),
                                         sep = "*\", \"*")),
                      formula = my_formula, parse = TRUE,
                      label.y = "top", label.x = "right", size = 4.5) +
         xlab("")+ylab("")+
         theme_PA(axs.T.siz=18,axs.t.siz=14,str.siz=14)+
         theme(legend.position = "none",
               plot.title =element_text(size=20),
               panel.background = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"),
               strip.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA, size=1.25))+
         labs(title=NM)+xlim(XLIM)+ylim(0,max(d.list$Mean.wght,na.rm=T))  
       return(p)
     }
     
   }
   for(s in 1:length(Logbook.sp)) 
   {
     print(paste("Change in mean weight of landed individual ","--",Logbook.sp[s]))
     dummy=print(fn.plt.mn.ktch.wght(d.list=Logbook%>%filter(SNAME==Logbook.sp[s]),
                                     NM=capitalize(Logbook.sp[s]),
                                     XLIM=c(min(Logbook$Finyear)-1,1+max(Logbook$Finyear))))
     if(!is.null(dummy)) Change.mean.weight.catch[[s]]=dummy
     rm(dummy)
   }
   figure=ggarrange(plotlist=fun.find.in.list(x=Change.mean.weight.catch, Drop=NULL))+
     theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
   annotate_figure(figure,
                   left = text_grob("Mean weight of caught individuals (kg)", rot = 90,size=25),
                   bottom = text_grob("Financial year",size=25))
   ggsave(paste(hNdl,'/Outputs/Figure 5_Changes in mean weight caught individual _TDGDLF.tiff',sep=''),
          width = 15,height = 10,compression = "lzw")
   
   
   do.model.based.mn.weight.ktch=FALSE
  if(do.model.based.mn.weight.ktch)
  {
    Mn.weit.ktch=lapply(Species.data, function(x) x[["annual.mean.size_relative"]])
    Mn.weit.ktch=fun.find.in.list(x=Mn.weit.ktch, Drop=names(Indicator.species))
    
    Change.mean.weight.catch=vector('list',length(Mn.weit.ktch))
    names(Change.mean.weight.catch)=names(Mn.weit.ktch)
    
    fn.plt.mn.ktch.wght=function(d,NM)
    {
      my_formula = y ~ x
      p=d%>%
        mutate(Finyear.f=factor(Finyear),
               Finyear=as.numeric(substr(Finyear,1,4)))%>%      
        ggplot(aes(x = Finyear, y = mean)) +
        geom_point(aes(colour=Finyear.f),size=3)+
        geom_smooth(color = "black", formula = my_formula, method = 'lm',se=F)+ 
        stat_poly_eq(aes(label =  paste(stat(eq.label),stat(adj.rr.label),stat(p.value.label),
                                        sep = "*\", \"*")),
                     formula = my_formula, parse = TRUE,
                     label.y = "bottom", label.x = "right", size = 5) +
        xlab("")+ylab("")+
        theme_PA(axs.T.siz=18,axs.t.siz=14,str.siz=16)+
        theme(legend.position = "none",
              plot.title =element_text(size=20),
              panel.background = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              strip.background = element_rect(fill = "white"),
              panel.border = element_rect(colour = "black", fill=NA, size=1.25))+
        labs(title=NM)+ylim(0,max(d$mean,na.rm=T)) 
      return(p)
    }
    for(s in 1:length(Mn.weit.ktch))
    {
      Change.mean.weight.catch[[s]]=print(fn.plt.mn.ktch.wght(d=Mn.weit.ktch[[s]],NM=capitalize(names(Mn.weit.ktch)[s])))
    }
    figure=ggarrange(plotlist=Change.mean.weight.catch)
    annotate_figure(figure,
                    left = text_grob("Relative mean weight of caught individuals", rot = 90,size=20),
                    bottom = text_grob("Financial year",size=20))
    ggsave(paste(hNdl,'/Outputs/Figure_Changes in mean weight caught individual _TDGDLF.tiff',sep=''),
           width = 10,height = 10,compression = "lzw")
    
  }
  
  #---RESULTS. SPM ------
  if(Do.SPM=="YES")
  {
    fn.cons.po=function(low,up) c(low, tail(up, 1), rev(up), low[1])  #construct polygon
    
    #Check which sp hit boundaries
    Hit.Bounds=data.frame(Name=names(SPM.preds),Hits.upper.bound=NA,Hits.lower.bound=NA)
    for(s in 1: N.sp)
    {
      HR.o.scens=Mx.init.harv[s]
      if(!is.null(SPM.preds[[s]]))
      {
        for(h in 1:length(HR.o.scens))
          for(e in 1:length(Efficien.scens))
          {
            Hit=Store.SPM[[s]][[h]][[e]]$hit.lower.boundary
            if(length(Hit)>0)Hit.Bounds$Hits.lower.bound[s]=Hit
            Hit=Store.SPM[[s]][[h]][[e]]$hit.upper.boundary
            if(length(Hit)>0)Hit.Bounds$Hits.upper.bound[s]=Hit
          }
      } 
    }
    Hit.Bounds=Hit.Bounds%>%
      filter(Name%in%names(SPM.preds[!sapply(SPM.preds, is.null)]))
    write.csv(Hit.Bounds,paste(hNdl,"/Outputs/SPM_Hitting boundaries.csv",sep=""),row.names = F)  
    
    Species.not.hitting.bounds=Hit.Bounds%>%
      filter(is.na(Hits.upper.bound) & is.na(Hits.lower.bound))%>%
      mutate(Name=as.character(Name))%>%                  
      pull(Name)
    
    #Plot obs VS pred cpues  
    fn.plt.cpue=function(ob,ob.CV,pred,Convergence,pred.fit)
    {
      iid=match(pred.fit,pred)
      Yr=all.iers
      ob.CV=ob.CV[iid]
      plot(Yr,pred,pch=19,type='l',lwd=2,ylab="",xlab="",
           ylim=c(min(c(ob-ob.CV,pred),na.rm=T),max(c(ob+ob.CV,pred),na.rm=T)))
      points(Yr[iid],ob,col="orange",pch=19,cex=1)
      segments(Yr[iid],ob,Yr[iid],ob+ob.CV,col="orange")
      segments(Yr[iid],ob,Yr[iid],ob-ob.CV,col="orange")
      #legend("bottomright",paste("convergence=",Convergence),bty='n')
    }
    #by species
    Paz=paste(hNdl,"/Outputs/SPM.fit/",sep="")
    if(!file.exists(file.path(Paz))) dir.create(file.path(Paz))
    for(s in 1: N.sp)
    {
      if(!is.null(SPM.preds[[s]]))    
      {
        fn.fig(paste(Paz,names(SPM.preds)[s],sep=""),2400,1400)
        
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
            {
              if(!is.null(SPM.preds[[s]][[h]][[e]]$ln.cpue[[x]]))
              {
                fn.plt.cpue(ob=SPM.preds[[s]][[h]][[e]]$ln.cpue[[x]],
                            ob.CV=Store.stuff[[s]]$CPUE.CV[[x]],
                            pred=SPM.preds[[s]][[h]][[e]]$ln.cpue.hat.full[[x]],
                            Convergence=Store.SPM[[s]][[h]][[e]]$convergence,
                            pred.fit=SPM.preds[[s]][[h]][[e]]$ln.cpue.hat[[x]])
                #Crip=Efficien.scens[e]
                #if(x==1) Crip=0
                #par(font=2)
                #  legend("bottomleft",paste("HR=",HR.o.scens[h],", Creep=",Crip,
                #                            ", CPUE.",x,sep=""),bty='n',cex=.85)
                legend('topleft',names(cpue.list[[s]])[x],bty='n')
              }
              
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
    
    #All species together
    nRows=length(Species.not.hitting.bounds)
    nCols=length(cpue.list$`tiger shark`)
    Shwed=which(names(SPM.preds)%in%Species.not.hitting.bounds)
    fn.fig(paste(Paz,"All.species_not.hitting.boundaries",sep=""),2400,1800)
    par(mfrow=c(nRows,nCols),mar=c(1,1,1,1),oma=c(1.8,2,.5,.1),mgp=c(1,.5,0))
    for(s in 1: N.sp)
    {
      if(names(SPM.preds[s])%in%Species.not.hitting.bounds)    
      {
        HR.o.scens=Mx.init.harv[s]
        for(h in 1:length(HR.o.scens))
        {
          for(e in 1:length(Efficien.scens))
          {
            for(x in 1:Store.stuff[[s]]$n.cpues)
            {
              if(!is.null(SPM.preds[[s]][[h]][[e]]$ln.cpue[[x]]))
              {
                fn.plt.cpue(ob=SPM.preds[[s]][[h]][[e]]$ln.cpue[[x]],
                            ob.CV=Store.stuff[[s]]$CPUE.CV[[x]],
                            pred=SPM.preds[[s]][[h]][[e]]$ln.cpue.hat.full[[x]],
                            Convergence=Store.SPM[[s]][[h]][[e]]$convergence,
                            pred.fit=SPM.preds[[s]][[h]][[e]]$ln.cpue.hat[[x]])
              }else
              {
                plot.new()
              }
              if(x==1) legend("bottomleft",capitalize(names(SPM.preds)[s]),bty='n',cex=1.5,text.font=2)
              if(s==Shwed[1])mtext(names(cpue.list[[s]])[x],3,line=0.2)
            }
          }
        }
        rm(HR.o.scens)
      }
    }  
    legend("left",c("observed","predicted"),pch=19,cex=1.25,col=c("orange","black"),bty='n')
    mtext("Financial year",1,.7,outer=T)
    mtext("lncpue",2,1,outer=T,las=3)
    dev.off()
    
    
    #Probability of above and below reference points     
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
            if(What.percentil=="100%") Bt=apply(dummy,2,function(x) quantile(x,probs=c(0,0.5,1)))   
            if(What.percentil=="60%") Bt=apply(dummy,2,function(x) quantile(x,probs=c(0.2,0.5,0.8)))   
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
    
    
    fn.fig(paste(hNdl,"/Outputs/Figure 3_Biomass_SPM_species.not.hitting.boundaries",sep=""),2000,2400)
    smart.par(n.plots=length(Species.not.hitting.bounds),MAR=c(1.2,2,1.5,1.25),
              OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
    par(las=1,cex.axis=.8)
    for(s in 1: N.sp)
    {
      if(names(SPM.preds[s])%in%Species.not.hitting.bounds)
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
    mtext("Financial year",1,cex=1.5,line=0.75,outer=T)
    mtext("Relative biomass",2,cex=1.5,outer=T,las=3)
    mtext(side = 4, line = 0.75, 'Total catch (tonnes)',las=3,outer=T,
          col=rgb(0.1,0.1,0.8,alpha=0.6),cex=1.5)
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
  #---RESULTS. Catch-MSY ------
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
           lwd=3,main=NMs,xlab="",ylab="",cex.lab=2,cex.axis=1.15,col=1,xlim=c(0,.6),yaxt='n')
    }
    mtext(expression(paste(plain("Intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),1,0.5,cex=1.35,outer=T)
    mtext("Density",2,0,las=3,cex=1.35,outer=T)
    dev.off()
    
    YrS=sort(unique(Tot.ktch$finyear))
    
    #Percentage of simulations accepted
    Per.accepted=data.frame(Species=capitalize(names(store.species)),
                            Percent.accepted=NA)
    for(s in 1: N.sp) Per.accepted$Percent.accepted[s]=100*ncol(store.species[[s]]$KTCH.MSY$BaseCase$bt)/SIMS
    write.csv(Per.accepted,paste(hNdl,"/Outputs/Per.accepted.CMSY.simulations.csv",sep=""),row.names = F)
    
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
    
    fn.fig(paste(hNdl,"/Outputs/Figure 3_Biomass_Catch_MSY",sep=""), 2400, 2000)
    smart.par(n.plots=length(compact(store.species)),MAR=c(1.2,2,1.5,2),
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
    
    #MSY
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
  #---RESULTS. aSPM ------
  if(Do.aSPM=="YES")
  {
    #MISSING: output steepness for relevant species as r prior 'Figure 1_Prior_r'
    
    #remove species with no aSPM
    no.null.Store.aSPM=Store.aSPM[!sapply(Store.aSPM, is.null)]
    
    #Plot obs VS pred cpues  
    Paz=paste(hNdl,"/Outputs/aSPM.fit/",sep="")
    if(!file.exists(file.path(Paz))) dir.create(file.path(Paz))
    
    nRows=length(no.null.Store.aSPM)
    nCols=2
    Shwed=names(Store.aSPM)
    fn.fig(paste(Paz,"All.species",sep=""),2400,1800)
    par(mfrow=c(nRows,nCols),mar=c(1,1,1,1),oma=c(1.8,2,.5,.1),mgp=c(1,.5,0))
    for(s in 1: length(Store.aSPM))
    {
      if(!is.null(Store.aSPM[[s]]))
      {
        aa=Store.aSPM[[s]]$quantities%>%filter(Year%in%all.iers)
        
        #SURVEY
        if(sum(aa$CPUE2,na.rm=T)>0)
        {
          plot(aa$Year,aa$PredCE2,type='l',lwd=2,ylab="",xlab="",xlim=c(min(aa$Year),max(aa$Year)),
               ylim=c(0,max(c(max(aa$PredCE2),aa$PredCE2+aa$se2),na.rm=T)))
          points(aa$Year,aa$CPUE2,col="orange",pch=19,cex=1)
          with(aa,segments(Year,CPUE2,Year,CPUE2+se2,col="orange"))
          with(aa,segments(Year,CPUE2,Year,CPUE2-se2,col="orange"))
        }else
          plot.new()
        legend("bottomleft",capitalize(names(Store.aSPM)[s]),bty='n',cex=1.5,text.font=2)
        if(names(Store.aSPM)[s]==names(no.null.Store.aSPM)[1]) mtext("Survey",3,line=0.2)
        
        #TDGDLF
        xx=match(2005,aa$Year)
        plot(aa$Year[1:xx],aa$PredCE[1:xx],type='l',lwd=2,ylab="",xlab="",xlim=c(min(aa$Year),max(aa$Year)),
             ylim=c(0,max(c(max(aa$PredCE),aa$PredCE+aa$se),na.rm=T)))
        lines(aa$Year[(xx+1):nrow(aa)],aa$PredCE[(xx+1):nrow(aa)],lwd=2)
        points(aa$Year,aa$CPUE,col="orange",pch=19,cex=1)
        with(aa,segments(Year,CPUE,Year,CPUE+se,col="orange"))
        with(aa,segments(Year,CPUE,Year,CPUE-se,col="orange"))
        if(names(Store.aSPM)[s]==names(no.null.Store.aSPM)[1]) mtext("TDGDLF",3,line=0.2)
      }
    }
    legend("bottomleft",c("observed","predicted"),pch=19,cex=1.25,col=c("orange","black"),bty='n')
    mtext("Financial year",1,.7,outer=T)
    mtext("cpue",2,1,outer=T,las=3)
    dev.off()
    
    
    #Plot relative spawning biomass  
    #1. Get median and percentiles    
    Med.biom_aSPM=vector('list',length(Store.aSPM))
    names(Med.biom_aSPM)=names(Store.aSPM)
    Store.cons.Like.aSPM=Store.cons.Like.aSPM.total=Tab.aSPM=Med.biom_aSPM
    dis.iers=match(all.iers,Store.aSPM$`tiger shark`$fish$year)
    for(s in 1: length(Store.aSPM))
    {
      if(!is.null(Store.aSPM[[s]]))
      {
        dummy=matrix(NA,nrow=NsimSS,ncol=length(all.yrs))
        dummy.total=dummy
        for(i in 1:NsimSS)
        {
          dummy[i,]=Store.aSPM[[s]]$quantities.MC[[i]]$Deplete[dis.iers]
          dummy.total[i,]=with(Store.aSPM[[s]]$quantities.MC[[i]],TotalB[dis.iers]/TotalB[1])
        }
        if(What.percentil=="100%")
        {
          Bt=apply(dummy,2,function(x) quantile(x,probs=c(0,0.5,1))) 
          Bt.total=apply(dummy.total,2,function(x) quantile(x,probs=c(0,0.5,1))) 
        }
        if(What.percentil=="60%")
        {
          Bt=apply(dummy,2,function(x) quantile(x,probs=c(0.2,0.5,0.8))) 
          Bt.total=apply(dummy.total,2,function(x) quantile(x,probs=c(0.2,0.5,0.8))) 
        }
        Med.biom_aSPM[[s]]=list(DAT=dummy,Bt=Bt,
                                DAT.total=dummy.total,Bt.total=Bt.total)
      }
    }
    
    #2. Plot
    fn.fig(paste(hNdl,"/Outputs/Figure 3_Biomass_aSPM_spawning",sep=""),2000,2400)
    smart.par(n.plots=nRows,MAR=c(1.2,2,1.5,1.25),
              OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
    par(las=1,cex.axis=.8)
    for(s in 1: length(Store.aSPM))
    {
      if(!is.null(Store.aSPM[[s]]))
      {
        Store.cons.Like.aSPM[[s]]=fn.plt.bio.ktch(Yr=Store.stuff[[s]]$yrs,
                                                  Bt=Med.biom_aSPM[[s]]$Bt,
                                                  Bmsy=1,
                                                  Ktch=Store.stuff[[s]]$Ktch,
                                                  CX.AX=1,
                                                  CX=1,
                                                  DAT=t(Med.biom_aSPM[[s]]$DAT),
                                                  Add.ktch="YES",
                                                  LOW="0%",
                                                  HIGH="100%")
        mtext(capitalize(names(Store.aSPM)[s]),3,cex=1) 
      }
    }
    mtext("Financial year",1,cex=1.2,line=0.75,outer=T)
    mtext("Relative spawning biomass",2,cex=1.2,outer=T,las=3)
    mtext(side = 4, line = 0.75, 'Total catch (tonnes)',las=3,outer=T,
          col=rgb(0.1,0.1,0.8,alpha=0.6),cex=1.2)
    dev.off()
    
    fn.fig(paste(hNdl,"/Outputs/Figure 3_Biomass_aSPM_total",sep=""),2000,2400)
    smart.par(n.plots=nRows,MAR=c(1.2,2,1.5,1.25),
              OMA=c(2,1.75,.2,2.1),MGP=c(1,.62,0))
    par(las=1,cex.axis=.8)
    for(s in 1: length(Store.aSPM))
    {
      if(!is.null(Store.aSPM[[s]]))
      {
        Store.cons.Like.aSPM.total[[s]]=fn.plt.bio.ktch(Yr=Store.stuff[[s]]$yrs,
                                                        Bt=Med.biom_aSPM[[s]]$Bt.total,
                                                        Bmsy=1,
                                                        Ktch=Store.stuff[[s]]$Ktch,
                                                        CX.AX=1,
                                                        CX=1,
                                                        DAT=t(Med.biom_aSPM[[s]]$DAT.total),
                                                        Add.ktch="YES",
                                                        LOW="0%",
                                                        HIGH="100%")
        mtext(capitalize(names(Store.aSPM)[s]),3,cex=1)
      }
    }
    mtext("Financial year",1,cex=1.2,line=0.75,outer=T)
    mtext("Relative total biomass",2,cex=1.2,outer=T,las=3)
    mtext(side = 4, line = 0.75, 'Total catch (tonnes)',las=3,outer=T,
          col=rgb(0.1,0.1,0.8,alpha=0.6),cex=1.2)
    dev.off()
    
    
    #Output parameter estimates
    for(s in 1: length(Store.aSPM))
    {
      if(!is.null(Store.aSPM[[s]]))
      {
        fit=Store.aSPM[[s]]$fit
        MLE=round(fit$par,2)
        std=round(sqrt(diag(solve(fit$hessian))),2)
        Nms=names(MLE)
        Tab=as.data.frame(matrix(paste(MLE," (",std,")",sep=''),nrow=1))
        names(Tab)=capitalize(Nms)
        Tab.aSPM[[s]]=Tab
      }
    }
    Tab.aSPM=do.call(rbind,Tab.aSPM)
    fn.word.table(WD=getwd(),TBL=Tab.aSPM,Doc.nm="Table 3. aSPM estimates",caption=NA,paragph=NA,
                  HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                  Zebra='NO',Zebra.col='grey60',Grid.col='black',
                  Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  }
  
  #---RESULTS. Selectivity ------
  #MISSING: display selectivity and Sel.equivalence only for those species for assessments used selectivity!!
  
  
  #---RESULTS. RISK ------
  LoE.Weights=c(psa=.2,sptemp=.2,efman=.2,spmod=.5,srmod=.5,aspmod=.5)
  
  #1. Calculate risk for each line of evidence
  #1.1. PSA (Drop.species only)
  Risk.PSA=vector('list',length(Drop.species))
  names(Risk.PSA)=Drop.species
  for(s in 1:length(Risk.PSA))
  {
    Risk.PSA[[s]]=data.frame(Max.Risk.Score=c(0,4,0,0))
  }
  
  #1.2. catch spatio-temporal      
  Risk.spatial.temporal.ktch=vector('list',N.sp)
  names(Risk.spatial.temporal.ktch)=names(Store.cons.Like.SRM)
  for(s in 1: N.sp)
  {
    sp.dat=Store.spatial.temporal.ktch[match(names(Risk.spatial.temporal.ktch)[s],
                                             rownames(Store.spatial.temporal.ktch)),]
    sp.dat=mean(sp.dat[(length(sp.dat)-4):length(sp.dat)]) #moving average
    
    
    if(sp.dat==0) dummy=data.frame(Max.Risk.Score=c(2,0,0,0))
    if(sp.dat>0 & sp.dat<=.25) dummy=data.frame(Max.Risk.Score=c(0,4,0,0))
    if(sp.dat>.25 & sp.dat<=.5) dummy=data.frame(Max.Risk.Score=c(0,0,6,0))
    if(sp.dat>.5 & sp.dat<=.75) dummy=data.frame(Max.Risk.Score=c(0,0,12,0))
    if(sp.dat>.75) dummy=data.frame(Max.Risk.Score=c(0,0,0,16))
    Risk.spatial.temporal.ktch[[s]]=dummy
  }
  
  #1.3. overall catch-effort-management               
  Risk.effort.mangmnt=vector('list',N.sp)
  names(Risk.effort.mangmnt)=names(Store.cons.Like.SRM)
  Rel.eff.n=mean(Effrt.n$Hook.hours[(nrow(Effrt.n)-4):nrow(Effrt.n)])/max(Effrt.n$Hook.hours) #moving average
  Rel.eff.s=mean(Effrt.s$Total[(nrow(Effrt.s)-4):nrow(Effrt.s)])/max(Effrt.s$Total)
  REgn=Tot.ktch%>%
    group_by(SP.group,Region)%>%
    summarise(Tot=sum(LIVEWT.c))%>%
    spread(Region,Tot)%>%
    mutate(Prop.N=North/(North+South),
           Prop.S=South/(North+South),
           Which.ef=ifelse(Prop.N>.7,'north',ifelse(Prop.S>.7,'south','north-south')))
  for(s in 1: N.sp)
  {
    Which.ef=REgn$Which.ef[match(names(Risk.effort.mangmnt)[s],REgn$SP.group)]
    if(Which.ef=='north') Which.ef=Rel.eff.n
    if(Which.ef=='south') Which.ef=Rel.eff.s
    if(Which.ef=='north-south') Which.ef=max(c(Rel.eff.s,Rel.eff.n))
    
    if(Which.ef==0) dummy=data.frame(Max.Risk.Score=c(2,0,0,0))
    if(Which.ef>0 & Which.ef<=.25) dummy=data.frame(Max.Risk.Score=c(0,4,0,0))
    if(Which.ef>.25 & Which.ef<=.5) dummy=data.frame(Max.Risk.Score=c(0,0,6,0))
    if(Which.ef>.5 & Which.ef<=.75) dummy=data.frame(Max.Risk.Score=c(0,0,12,0))
    if(Which.ef>.75) dummy=data.frame(Max.Risk.Score=c(0,0,0,16))
    
    Risk.effort.mangmnt[[s]]=dummy
  }
  
  #1.4. population dynamics models
  Like.ranges=list(L1=c(0,0.0499999),
                   L2=c(0.05,0.2),
                   L3=c(0.20001,0.5),
                   L4=c(0.50001,1))
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
  if(Do.aSPM=="YES") Risk.aSPM=Store.cons.Like.aSPM
  for(s in 1: N.sp)
  {
    #SPM
    if(names(Store.cons.Like.SPM)[s]%in%Species.not.hitting.bounds) Risk.SPM[[s]]=fn.risk(likelihood=unlist(Store.cons.Like.SPM[[s]]))
    
    #SRM
    Risk.SRM[[s]]=fn.risk(likelihood=unlist(Store.cons.Like.SRM[[s]]))
    
    #aSPM
    if(Do.aSPM=="YES" & !is.null(Store.aSPM[[s]])) Risk.aSPM[[s]]=fn.risk(likelihood=unlist(Store.cons.Like.aSPM[[s]]))
    
  }
  
  
  #2. Integrate the risk from each line of evidence
  
  #note: Use a weighted sum to aggregate the Risk Categories form the alternative lines of evidence
  #      Normalize each criterion by dividing by the highest value of each criterion. 
  #      Assign weights to each criteria 
  Order=c('Negligible','Low','Medium','High','Severe')
  Order <- factor(Order,ordered = TRUE,levels = Order)
  Integrate.LoE=function(Cons.Like.tab,criteriaMinMax,plot.data,LoE.weights,Normalised)
  {
    #Set up preference table by converting Cons-Like to Risk scores
    Preference.Table=as.data.frame(matrix(0,nrow=5,ncol=ncol(Cons.Like.tab)))
    colnames(Preference.Table)=colnames(Cons.Like.tab)
    rownames(Preference.Table)=c("Negligible","Low","Medium","High","Severe")
    for(i in 1:ncol(Preference.Table))
    {
      dd=Cons.Like.tab[,i]
      if(max(dd[1:2])<=2) Preference.Table["Negligible",i]=max(dd[1:2])
      if(max(dd)>2 & max(dd)<=4) Preference.Table["Low",i]=4
      if(max(dd[2:4])>4 & max(dd[2:4])<=8) Preference.Table["Medium",i]=max(dd[2:4])
      if(max(dd[3:4])>8 & max(dd[3:4])<=12) Preference.Table["High",i]=max(dd[3:4])
      if(dd[4]>12 & dd[4]<=16) Preference.Table["Severe",i]=dd[4]
    }
    
    #Maximise or minimise each criteria?
    criteriaMinMax=rep(criteriaMinMax,ncol(Preference.Table))
    names(criteriaMinMax) <- colnames(Preference.Table)
    
    #display data
    if(plot.data)plotRadarPerformanceTable(Preference.Table, criteriaMinMax,overlay=FALSE, bw=TRUE, lwd =5)
    
    # Normalization of the performance table
    normalizationTypes <- rep("percentageOfMax",ncol(Preference.Table))
    names(normalizationTypes) <- colnames(Preference.Table)
    if(Normalised=="YES") nPreference.Table <- normalizePerformanceTable(Preference.Table,normalizationTypes)
    if(Normalised=="NO") nPreference.Table=Preference.Table
    
    # Calculate weighted sum
    names(LoE.weights) <- colnames(nPreference.Table)
    weighted.sum<-weightedSum(nPreference.Table,LoE.weights)
    
    # Rank the scores of the alternatives
    rank.score=sort(rank(-weighted.sum))
    
    # overall risk
    risk=which(rank.score==min(rank.score))
    risk=as.character(max(Order[match(names(risk),Order)]))
    
    return(list(weighted.sum=weighted.sum,rank.score=rank.score,risk=risk))
  }
  
  All.sp=sort(c(Drop.species,Keep.species))
  Overall.risk=vector('list',length(All.sp))
  names(Overall.risk)=All.sp
  for(s in 1:length(All.sp))
  {
    dummy=list(psa=Risk.PSA[[match(All.sp[s],names(Risk.PSA))]]$Max.Risk.Score,
               sptemp=Risk.spatial.temporal.ktch[[match(All.sp[s],names(Risk.spatial.temporal.ktch))]]$Max.Risk.Score,
               efman=Risk.effort.mangmnt[[match(All.sp[s],names(Risk.effort.mangmnt))]]$Max.Risk.Score,
               spmod=Risk.SPM[[match(All.sp[s],names(Risk.SPM))]]$Max.Risk.Score,
               srmod=Risk.SRM[[match(All.sp[s],names(Risk.SRM))]]$Max.Risk.Score)
    if(Do.aSPM=="YES")dummy$aspmod=Risk.aSPM[[match(All.sp[s],names(Risk.aSPM))]]$Max.Risk.Score
    dummy=dummy[!sapply(dummy, is.null)]
    NMs=names(dummy)
    dummy=do.call(cbind,dummy)
    colnames(dummy)=NMs
    rownames(dummy)=paste("C",1:4,sep='')
    
    Overall.risk[[s]]=Integrate.LoE(Cons.Like.tab=dummy,
                                    criteriaMinMax <- "max",
                                    plot.data=FALSE,
                                    LoE.weights <- LoE.Weights[match(colnames(dummy),names(LoE.Weights))],
                                    Normalised="YES")
  }
  
  
  #3. Display overall risk for each species
  fn.each.LoE.risk=function(N.sp)
  {
    X.rng=1:N.sp
    x.Vec <-  fn.cons.po(0:(N.sp+1),0:(N.sp+1))
    nn=N.sp+2
    negigible.Vec <- fn.cons.po(rep(0,nn),rep(2,nn))
    low.Vec <- fn.cons.po(rep(2,nn),rep(4,nn))
    medium.Vec <- fn.cons.po(rep(4,nn),rep(8,nn))
    high.Vec <- fn.cons.po(rep(8,nn),rep(12,nn))
    severe.Vec <- fn.cons.po(rep(12,nn),rep(16,nn))
    plot(X.rng,xlim=c(0,16),ylim=c(0,N.sp+1),xaxs="i",yaxs="i",
         col="transparent",ylab="",xlab="",xaxt='n',yaxt='n')
    polygon(negigible.Vec, x.Vec, col = 'cornflowerblue', border = "transparent")
    polygon(low.Vec, x.Vec, col = 'olivedrab3', border = "transparent")
    polygon(medium.Vec, x.Vec, col = 'yellow', border = "transparent")
    polygon(high.Vec, x.Vec, col = 'orange', border = "transparent")
    polygon(severe.Vec, x.Vec, col = 'red', border = "transparent")
    
    axis(1,at=c(mean(negigible.Vec),mean(low.Vec),mean(medium.Vec),mean(high.Vec),mean(severe.Vec)),
         labels=c("Negl.","Low","Medium","High","Severe"))
    box()
  }
  
  fn.overall.risk=function(N,RISK,sp)
  {
    x.Vec <-  c(1,3,3,1)
    y.Vec <- c(rep((s-.5),2),rep((s+.5),2))
    if(RISK=="Negligible") CL = 'cornflowerblue'
    if(RISK=="Low") CL = 'olivedrab3'
    if(RISK=="Medium") CL ='yellow'
    if(RISK=="High") CL ='orange'
    if(RISK=="Severe") CL ='red'
    polygon(x.Vec, y.Vec, col = CL, border = "transparent")
    text(1.5,mean(y.Vec),sp,cex=1.5)
  }
  
  LoE.col=c(psa="grey85", sptemp="grey75", efman="grey55",
            spmod="grey40", srmod="grey20", aspmod="black")
  
  Sp.risk.ranking=factor(unlist(lapply(Overall.risk, '[[', 'risk')),levels=levels(Order))  
  Sp.risk.ranking=names(sort(Sp.risk.ranking))
  
  
  fn.fig("Figure 5_Risk", 2400, 2300)
  par(mar=c(.5,.5,3,1),oma=c(3,13.5,.5,.1),las=1,mgp=c(1,.5,0),cex.axis=1.5,xpd=TRUE)
  layout(matrix(c(rep(1,6),rep(2,3)),ncol=3))
  
  #Risk for each line of evidence
  fn.each.LoE.risk(N.sp=length(Sp.risk.ranking))
  for(s in 1:length(Sp.risk.ranking))
  {
    ss=Sp.risk.ranking[s]
    dummy=list(psa=Risk.PSA[[match(ss,names(Risk.PSA))]]$Max.Risk.Score,
               sptemp=Risk.spatial.temporal.ktch[[match(ss,names(Risk.spatial.temporal.ktch))]]$Max.Risk.Score,
               efman=Risk.effort.mangmnt[[match(ss,names(Risk.effort.mangmnt))]]$Max.Risk.Score,
               spmod=Risk.SPM[[match(ss,names(Risk.SPM))]]$Max.Risk.Score,
               srmod=Risk.SRM[[match(ss,names(Risk.SRM))]]$Max.Risk.Score)
    if(Do.aSPM=="YES")dummy$aspmod=Risk.aSPM[[match(ss,names(Risk.aSPM))]]$Max.Risk.Score
    dummy=dummy[!sapply(dummy, is.null)]
    NMs=names(dummy)
    dummy=do.call(cbind,dummy)
    colnames(dummy)=NMs
    dummy=apply(dummy,2,max)
    Nd=length(dummy)
    if(Nd==1)Adjst=0 
    if(Nd==3)Adjst=seq(-.3,.3,length.out = Nd)
    if(Nd==4)Adjst=seq(-.3,.3,length.out = Nd)
    if(Nd==5)Adjst=seq(-.3,.3,length.out = Nd)
    dummy=data.frame(LoE=names(dummy),
                     Start=0,
                     End=dummy,
                     y=s+Adjst)
    CLL=LoE.col[match(dummy$LoE,names(LoE.col))]
    segments(dummy$Start,dummy$y,dummy$End,dummy$y,lwd=3.75,lend=1,col=CLL)
  }
  axis(2,1:length(Sp.risk.ranking),capitalize(Sp.risk.ranking))
  mtext("Risk score",1,cex=1.25,line=2)
  legend(-1,length(Sp.risk.ranking)+3.25,c('PSA','Blocks fished','Effort management'),
         bty='n',col=LoE.col[1:3],lty=1,lwd=3,horiz = T,cex=1.5,
         text.width=c(0,1,2.25))
  if(Do.aSPM=="NO") legend(5,length(Sp.risk.ranking)+2.35,c('SPM','SRM'),bty='n',cex=1.5,
                           col=LoE.col[4:5],lty=1,lwd=3,horiz = T,text.width=c(0,.95))
  if(Do.aSPM=="YES") legend(5,length(Sp.risk.ranking)+2.35,c('SPM','SRM','aSPM'),bty='n',cex=1.5,
                            col=LoE.col[4:6],lty=1,lwd=3,horiz = T,text.width=c(0,.95,1.1))
  
  #Overall risk
  plot(0:1,ylim=c(0,length(All.sp)+1),fg='white',xaxs="i",yaxs="i",
       col="transparent",ylab="",xlab="",xaxt='n',yaxt='n')
  for(s in 1:length(Sp.risk.ranking))
  {
    ss=Sp.risk.ranking[s] 
    fn.overall.risk(N=s,RISK=Overall.risk[[match(ss,names(Overall.risk))]]$risk,sp=capitalize(Sp.risk.ranking[s]))
  }
  axis(1,1.5,"Overall risk",cex.axis=1.75,col.ticks="white",padj=.5)
  dev.off()
  
  
}


#---RESULTS. Indicator species ------------------------------------------------- 
