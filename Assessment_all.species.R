# ------ Script for running stock assessments for WA sharks---- ###################

#MISSING:
#   Include ALL species in final risk scoring
#   Double check that TwT calculation uses TL and not FL, be consistent
#   Indicator species:
#       set up integrated model for dusky and sandbar
#       set up SS3 model for indicator species
#       Use .CPUE_Observer_TDGDLF.csv in likelihoods?

#   Other species/Catch only:
#       review Smooth HH cpue and mako cpue...; SPM Tiger fit
#       Milk shark SPM, hitting upper K boundary, no trend in cpue, crap Hessian, too uncertain....mention in text...
#       aSPM: finish running for all species; issues with Tiger cpue fit...
#       Size-based Catch curve (for some species there's NSF size compo, not used at the moment)
#       Implement simple SS?
#       NT only provided catches for sandbar and dusky shark 



#Notes:
#     1. The PSA filters out species from further analyses thru 'Criteria for selecting what species to assess'
#     2. Assumption: If catches have never been >1% carrying capacity, then it's in unexploited 
#                    status so catch series have no information on productivity
#     3. Total (reconstructed) catches are used (commercial, recreational and TDGDLF discards);
#             recons_NT_catch.csv' only includes dusky and sandbar
#     4. A range of assessment methods are used depending on data availability: 
#               . Integrated size- and sex- structured models,
#               . SPM,
#               . Catch-only
#               . aSPM(Haddon SimpleSA())



#Steps: 
#     1. Update year of assessment in '1. DEFINE GLOBALS'  
#     2. Define arguments used in each of the shark species/species complex assessed.
#     3. Bring in updated available data and input parameters
#     4. Determine which species to assess based on PSA
#     5. Run relevant population models according to data availability
#     6. Generate relevant outputs

# ------ Header---- 
rm(list=ls(all=TRUE))
options(dplyr.summarise.inform = FALSE)
options(stringsAsFactors = FALSE) 
options(ggrepel.max.overlaps = Inf)
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

# install.packages("devtools")  #  unhash these if you do not have devtools or
# install.packages("Rcpp")      #  Rcpp installed
library(MASS)
library(plotrix)
library(PBSmapping)
library(tidyverse)
library(dplyr)
library(mvtnorm)
library(rlist)
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
#devtools::install_github("haddonm/datalowSA",build_vignettes=TRUE,force=TRUE)
library(scales)
library(gridExtra)
library(purrr)
library(yarrr)
library(coda)
library(ggmcmc)

source.hnld=handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/")
fn.source=function(script)source(paste(source.hnld,script,sep=""))
fn.source("fn.fig.R")
fn.source("Catch_MSY.R")
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")) 

smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
colfunc <- colorRampPalette(c("red","yellow","springgreen","royalblue"))
fun.find.in.list=function(x,Drop=NULL)   #drop stuff from list
{
  x=x%>%purrr::discard(is.null)
  if(!is.null(Drop))x=x[-match(Drop,names(x))]
  return(x)
}

fn.get.stuff.from.list=function(lista,stuff) lapply(lista, function(x) x[[stuff]])# get stuff from list

source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/ggplot.themes.R'))  #my themes

#---1.  DEFINE GLOBALS----- 

#1.1. General
  #Year of assessment 
Year.of.assessment=2022
Asses.year=Year.of.assessment
AssessYr=Year.of.assessment

  #Last complete financial year of catches
Last.yr.ktch="2019-20"

  #New assessment
New.assessment="NO"
#New.assessment="YES"   #set to 'YES' the first time new data (catch at least) are available to decide what species to assess

  #Add additional species not selected by PSA given low catch trajectories but needed for specific assessment
additional.sp=NULL  #if no additional species assessment required
if(Year.of.assessment==2022) additional.sp=c('green sawfish','narrow sawfish')   # 2022 sawfish assessment; 
                                                                  # dwarf and freshwater sawfish not assessed, 
                                                                  # recons catch likely to be incomplete (no TO catch or beach rec fishing)

  # Control which assessment methods to implement
Do.SPM="YES"
Do.Ktch.MSY="YES"
Do.aSPM="YES"

  #Model run
if(New.assessment=="YES") First.run="YES"  else #create all model input data sets and data presentation for new assessment
                          First.run="NO"
  #Paper outputs
Modl.rn="standard"   #for annual assessments
#Modl.rn='first'    #for paper

  #Define if calculating r, steepness and size-based catch curve
if(New.assessment=="YES") do.r.prior=TRUE  else 
                          do.r.prior=FALSE

do.steepness=do.r.prior  

do.Size.based.Catch.curve=FALSE  #superseded by dynamic catch-only and size model


  #Define if exporting figures as jpeg or tiff (creation of RAR requires jpeg)
Do.tiff="YES" 
Do.jpeg="NO"

  #Catch in tonnes?  
KTCH.UNITS="TONNES" 
#KTCH.UNITS="KGS"    
if(KTCH.UNITS=="KGS") unitS=1000
if(KTCH.UNITS=="TONNES") unitS=1


#1.2. Criteria for selecting what species to assess quantitatively
Min.yrs=5
if(KTCH.UNITS=="KGS") Min.ktch=5000 
if(KTCH.UNITS=="TONNES") Min.ktch=5

  #Proportion of vessels discarding eagle rays in last 5 years (extracted from catch and effort returns)
prop.disc.ER=.4  

  #PSA
PSA.min.tons=5
PSA.min.years=Min.yrs
PSA.max.ton=50
Low.risk=2.64  #risk thresholds from Micheli et al 2014
medium.risk=3.18


#1.3. Size composition
  #Initial bin size
MN.SZE=0
#MN.SZE="size.at.birth"

  #size bin
TL.bins.cm=5

  #Minimun number of samples of size composition
Min.obs=10  #at least 10 observations
Min.shts=5  #from at least 5 shots


# 1.4. Global arguments for indicator species 
if(First.run=="YES") run.all.scenarios="YES" #What scenarios to run?
if(First.run=="NO") run.all.scenarios="NO"  
if(First.run=="NO") run.future="YES" #Future projection scenarios
if(First.run=="YES") run.future="NO"
Show.yrs="DATA"   #What years to show in model outputs
#Show.yrs="FUTURE"  #show data plus future projections
Present.in.log="NO"  #present cpue in log space or normal space


  #1.4.1 Size-based model
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

  #1.4.2. SPM
r_max=0.5     #Include a prior for r in surplus production. This is the max reported value of r for sharks (blue shark)
Add.r.prior=0   #no r prior 
#Add.r.prior=1   #use r prior


#1.5. Global arguments for non-indicator species
Asses.Scalloped.HH=FALSE  #2020 scalloped HH assessment

  # Minimun number of annual observations for change in observed mean length analysis
Min.annual.obs=100

  # Minimun number of annual observations for dynamic catch-only size based analysis
Min.annual.obs_catch.curve=300

#1.6. Catch only Models
do.ensemble.simulations=FALSE
do.OCOM=FALSE    
catch.only=c('DBSRA','CMSY','JABBA')
if(do.OCOM) catch.only=c(catch.only,'OCOM')

  # Catch MSY
CMSY.method="Haddon"    # Haddon's datalowSA
#CMSY.method="Froese"  # Froese et al 2017 does not converge for dwarf or freshwater


#1.7. Reference points
#note: Historically, single reference point for the fishery (biomass at 40% unexploited conditions)
Biomass.threshold='Bmsy'
Tar.prop.bmsny=1.3    # Target and Limit proportions of Biomass.threshold. 
Lim.prop.bmsy=0.5    #    source: Haddon et al 2014. 'Technical Reviews of Formal Harvest Strategies'.

    #Empirical reference points
Fmsy.emp=function(M) 0.41*M     #Zhou et al 2012 but see Cortes & Brooks 2018
SPR.thre=0.3   #Punt 2000 Extinction of marine renewable resources: a demographic analysis. 
SPR.tar=0.4                # Population Ecology 42, 


#1.8. Define if using effort
add.effort="NO"    
What.Effort="km.gn.hours"  #What effort to display?
#What.Effort="km.gn.days" 


#1.9. Define if assessing 'other species'
do.other.ass=TRUE


#1.10. Define if doing other species paper
do.other.ass.paper=FALSE  


#1.11. Assumed PCM for reconstructed discards in TDGLDF
TDGLDF.disc.assumed.PCM="BaseCase" 
#TDGLDF.disc.assumed.PCM="100%" 


#1.12. Demography
Max.r.value=.51 # Max r. .44 for blue shark (Cortes 2016); .51 for Scyliorhinus canicula (Cortes 2002)


#---2.  Catch and effort data -----   
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
  
  #..Taiwanese gillnet, longline and trawl
  Taiwan.gillnet.ktch=fn.in(NM='recons_Taiwan.gillnet.ktch.csv')
  Taiwan.longline.ktch=fn.in(NM='recons_Taiwan.longline.ktch.csv')
  Taiwan.trawl.ktch=fn.in(NM='recons_Taiwan.trawl.ktch.csv')
  
  Taiwan.gillnet.ktch$Method="GN"   #Pelagic.gillnet
  Taiwan.longline.ktch$Method="LL"  #Longline
  Taiwan.trawl.ktch$Method="TW"     #Trawl
  
  Taiwan=rbind(Taiwan.longline.ktch,Taiwan.gillnet.ktch,Taiwan.trawl.ktch)%>%
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
  
  #..NSW Fisheries (with beach protection)
  Bronze.whaler_NSW=fn.in(NM='recons_Bronzewhaler_NSW.csv') 
  
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
           METHOD="LL",     #line
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
           METHOD="LL",   #line
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
           METHOD="LL",    #line
           FishCubeCode="SA MSF")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  
    #Add Bronze.whaler_NSW (this includes beach protection)
  a=Bronze.whaler_NSW%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="NSW fisheries",
           Data.set="NSW fisheries",
           METHOD="mixed",    #line, trawl, etc
           FishCubeCode="NSW fisheries")%>%
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

  
  #For assessments, combine certain species with unreliable reporting resolution
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
    mutate(Gear=ifelse(METHOD%in%c("BS","BH","GN","HN","Pelagic.gillnet","PS"),"net",
                ifelse(METHOD%in%c("DL","DV","EL","GL","HL","HR","HY",
                                          "LL","Longline","Rec.line","TL"),'line',
                ifelse(METHOD%in%c("FG","TW"),'trawl',
                ifelse(METHOD%in%c("FT","PT","CT"),'trap',
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

ktch.combined=KtCh%>%
        group_by(SPECIES,Name,finyear)%>%
        summarise(Tonnes=sum(LIVEWT.c,na.rm=T))

#2.4. Remove species assessed elsewhere
KtCh=KtCh%>%filter(!Name%in%assessed.elsewhere)
KtCh.zone=KtCh.zone%>%filter(!Name%in%assessed.elsewhere)
ktch.combined=ktch.combined%>%filter(!Name%in%assessed.elsewhere)

#2.5 Data used for analysis of changes in reported mean weight
Wei.range=read.csv(handl_OneDrive("Data/Length_Weights/Data.Ranges.csv"),stringsAsFactors = F)
Wei.range.names=read.csv(handl_OneDrive("Data/Length_Weights/Species.names.csv"),stringsAsFactors = F)
Wei.range=merge(Wei.range,Wei.range.names,by="Sname",all.x=T)

Logbook=read.csv(handl_OneDrive("Analyses/Catch and effort/Logbook.data.mean.weight.csv"))


#---3.  Life history data ------------------------------------------------------  
LH.data=read.csv(handl_OneDrive('Data/Life history parameters/Life_History.csv'))

#---4.  PSA (which species to assess quantitatively) -----------------------------------------------  
#note: run a PSA aggregating the susceptibilities of multiple fleets (Micheli et al 2014)

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

ggsave(paste(Exprt,'Annual_ktch_by_species.used.in.PSA.tiff',sep='/'), width = 17,height = 7.5, dpi = 300, compression = "lzw")

#Export table of species identified in the catch by fishery  (all species including indicator species)
write.csv(Get.ktch$Table1%>%
            filter(!is.na(Name)),
          paste(Exprt,'Table S1_All.species.caught.by.fishery.csv',sep='/'),row.names = F)


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
    geom_text_repel(aes(size=Fnt.size),show.legend  = F,segment.colour="grey75",col='black',box.padding = line.sep) + 
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

RAR.species=Keep.species 
if(!is.null(additional.sp)) Keep.species=c(Keep.species,additional.sp)
Keep.species=sort(Keep.species)
N.sp=length(Keep.species)
Other.species=RAR.species[-match(names(Indicator.species),RAR.species)]

#---5.  Create functions for applying CMSY, DB-SRA, OCOM, JABBA -------------------------------------------------------
apply.DBSRA=function(year,catch,catchCV,catargs,agemat,k,b1k,btk,fmsym,bmsyk,M,graph,
                     nsims=Niters,grout,WD,outfile)
{
  setwd(WD)  #dbsra automatically exports the biomass trajectories
  #store inputs
  input=list(year=year,
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
  #pdf(paste(outfile,".pdf",sep=''))
  output <- dbsra(year=year,
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
  #dev.off()
  
  
  #extract biomass trajectories
  Biom.traj=read.csv(paste(WD,"Biotraj-dbsra.csv",sep='/'),header=FALSE)
  
  output$Biom.traj=Biom.traj%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)
  
  output$Depletion.traj=Biom.traj%>%
    mutate_at(vars(- V1), ~ . / output$Values$K)%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)
  
  output$B.Bmsy=Biom.traj%>%
    mutate_at(vars(- V1), ~ . / output$Values$Bmsy)%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)  
  
  
  #extract F trajectories
  U.series=Biom.traj[,-ncol(Biom.traj)]
  U.series[,-1]=mapply('/',catch,U.series[,-1])
  
  F.series=U.series%>%
    mutate_at(vars(- V1), function(x) -log(1-x))
  
  output$F.series=F.series%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)
  #output$F.series=F.from.U((catch/apply(output$Biom.traj,2,median)[1:length(catch)]))
  
  output$F.Fmsy=F.series%>%
    mutate_at(vars(- V1), ~ . / output$Values$Fmsy)%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)
  
  #extract years
  output$Years=year
  
  
  return(list(input=input,output=output))
}
Res.fn=function(r,Def)
{
  if(Def=="Martell")
  {
    if(r<0.1) out="verylow" else if
    (r>=0.1 & r <0.6) out="low" else if
    (r>=0.6 & r <1.5 ) out="medium" else if
    (r >=1.5) out="high"
  }

  if(Def=="Haddon")
  {
    if(r<0.1) out="verylow" else if
    (r <= mean(c(0.1,0.6))) out="low" else if
    (r>mean(c(0.1,0.6)) & r <=mean(c(0.3,0.8)) ) out="medium" else if
    (r >mean(c(0.3,0.8))) out="high"
  }
  return(out)
}
apply.CMSY=function(year,catch,r.range,k.range,Bo.low,Bo.hi,Int.yr=NA,Bint.low=NA,Bint.hi=NA,
                    Bf.low,Bf.hi,outfile,CMSY=CMSY.method,nsims=Niters,Proc.error,RES)
{
  input=list(r.range=r.range,
               stb.low=Bo.low,
               stb.hi=Bo.hi,
               int.yr=Int.yr,
               intb.low = Bint.low,
               intb.hi = Bint.hi,
               endb.low=Bf.low,
               endb.hi=Bf.hi,
              RES=RES,
              Proc.error=Proc.error,
              k.range=k.range)
  
  if(CMSY=="Froese")
  {
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
    
    #Appendix
    appendix=data.frame(r=output$r_viable,k=1000*output$k_viable)
    p1=appendix%>%
      ggplot(aes(r))+
      geom_histogram()
    p2=appendix%>%
      ggplot(aes(k))+
      geom_histogram()
    p3=appendix%>%
      ggplot(aes(k,r))+
      geom_point() +
      geom_density_2d_filled(alpha = 0.5)+
      geom_density_2d(size = 0.25, colour = "black")
    ggarrange(plotlist = list(p1,p2,p3),nrow=1,widths = c(1,1,2))
    ggsave(paste(outfile,'.tiff',sep=''),width = 14,height = 8, dpi = 300, compression = "lzw")  
    
  }
  
  if(CMSY=="Haddon")
  {
    glb=list(resilience=RES,spsname="")
    indat=data.frame(year=year,catch=catch)
    output <- run_cMSY(indat=indat,
                        glob=glb,
                        n = nsims,
                        incB = 0.025,
                        sigpR = Proc.error,
                        multK = 1,
                        finaldepl = c(Bf.low,Bf.hi),
                        start_k = k.range,
                        start_r = r.range,
                        initdepl = c(Bo.low,Bo.hi),
                        maximumH = 1,
                        Hyear = NA)
    #Appendix
    ciMSY=summarycMSY(output,indat,final=TRUE)
    pdf(paste(outfile,'.pdf',sep=''))
    plotcMSY6(ciMSY,indat[,"catch"])
    dev.off()
    
    #get B.bmsy and F.fmsy
    d1=pulloutStats(output$R1,probabs=c(0.5))
    d1=d1$traject
    d1=data.frame(d1)
    Kei=d1$K
    Ar=d1$r
    Fmsy=Ar/2
    Bmsy=Kei/2
    d1=d1%>%
      dplyr::select(-c(r,K,bd))
    U.series=mapply('/',catch,d1[,1:length(catch)])
    F.series=apply(U.series,2,function(x) -log(1-x))
    
    output$B.traj=d1
    output$Depletion.traj=d1/Kei
    output$Bmsy=Bmsy
    output$B.Bmsy=d1/Bmsy
    output$F.Fmsy=F.series/Fmsy
    output$F.traj=F.series
    output$Years=year
    output$acceptance.rate=round(100*length(ciMSY$r)/nsims,2)  
  }
  return(list(input=input,output=output))
}

apply.OCOM=function(year,catch,M,outfile)
{
  input= list(m=M)
  output <- ocom(year=year,catch=catch,m=M)

  # Extract reference points and time series from output
  output$ref_pts <- output[["ref_pts"]]
  output$ref_ts <- output[["ref_ts"]]
  
  #Appendix
  pdf(paste(outfile,'.pdf',sep=''))
  plot_dlm(output)
  dev.off()
  
  appendix=data.frame(r=output$krms_draws$r,k=output$krms_draws$k)
  p1=appendix%>%
    ggplot(aes(r))+
    geom_histogram()
  p2=appendix%>%
    ggplot(aes(k))+
    geom_histogram()
  p3=appendix%>%
    ggplot(aes(k,r))+
    geom_point() +
    geom_density_2d_filled(alpha = 0.5)+
    geom_density_2d(size = 0.25, colour = "black")
  ggarrange(plotlist = list(p1,p2,p3),nrow=1,widths = c(1,1,2))
  ggsave(paste(outfile,'.tiff',sep=''),width = 14,height = 8, dpi = 300, compression = "lzw")  
  
  
  return(list(input=input,output=output))
}

apply.JABBA=function(Ktch,CPUE,CPUE.SE,MDL,Ktch.CV,ASS,Rdist,Rprior,Kdist,Kprior,
                     PsiDist,Psiprior,Bprior,BMSYK,output.dir,outfile,Sims,Proc.error.JABBA,
                     thinning = 5,nchains = 2,burn.in=5000)
{
  # Compile JABBA JAGS model
  if(is.null(CPUE))
  {
    jbinput = build_jabba(catch=Ktch,
                          model.type = MDL,
                          catch.cv=Ktch.CV,
                          assessment=ASS,
                          scenario =  "CatchOnly",
                          r.dist = Rdist,
                          r.prior = Rprior,
                          K.dist= Kdist,
                          K.prior=Kprior,
                          psi.dist=PsiDist,
                          psi.prior=Psiprior,
                          b.prior=Bprior,
                          BmsyK = BMSYK,
                          sigma.est = FALSE,
                          sigma.proc=Proc.error.JABBA)
  }else
  {
    jbinput = build_jabba(catch=Ktch,
                          cpue=CPUE,
                          se=CPUE.SE,
                          model.type = MDL,
                          catch.cv=Ktch.CV,
                          assessment=ASS,
                          scenario =  "CPUE",
                          r.dist = Rdist,
                          r.prior = Rprior,
                          K.dist= Kdist,
                          K.prior=Kprior,
                          psi.dist=PsiDist,
                          psi.prior=Psiprior,
                          b.prior= Bprior,
                          BmsyK = BMSYK,
                          sigma.est = FALSE,
                          sigma.proc=Proc.error.JABBA) # Fixed Process error (set to TRUE if estimating) 
  }

  # Fit JABBA
  output=fit_jabba(jbinput=jbinput,
                   ni = Sims,   #number of iterations    
                   nt = thinning,   #thinning interval
                   nb = burn.in,  #burn-in
                   nc = nchains,    #number of chains
                   init.values = FALSE,    #set to TRUE if specifying initial value
                   init.K = mean(Kprior),
                   init.r = Rprior[1],
                   init.q = 1e-3,
                   save.all=TRUE,
                   save.csvs=TRUE,
                   output.dir=output.dir)
  
  
  #fit diagnostics
  setwd(output.dir)
  pdf(paste(outfile,'.pdf',sep=''))
  jbplot_catch(output)
  jbplot_ppdist(output)
  jbplot_mcmc(output)
  jbplot_procdev(output)
  jbplot_bprior(output)
  jbplot_trj(output,type="BBmsy",add=T)
  jbplot_trj(output,type="FFmsy",add=T)
  
  # status summary
  par(mfrow=c(3,2),mar = c(4, 5, 0.5, 0.1))
  jbplot_trj(output,type="B",add=T)
  jbplot_trj(output,type="F",add=T)
  jbplot_trj(output,type="BBmsy",add=T)
  jbplot_trj(output,type="FFmsy",add=T)
  jbplot_spphase(output,add=T)
  jbplot_kobe(output,add=T)
  dev.off()
  
  p1=output$pars_posterior%>%
    ggplot(aes(r))+
    geom_histogram()
  p2=output$pars_posterior%>%
    ggplot(aes(K))+
    geom_histogram()
  p3=output$pars_posterior%>%
    ggplot(aes(K,r))+
    geom_point() +
    geom_density_2d_filled(alpha = 0.5)+
    geom_density_2d(size = 0.25, colour = "black")
  ggarrange(plotlist = list(p1,p2,p3),nrow=1,widths = c(1,1,2))
  ggsave(paste(outfile,'.tiff',sep=''),width = 14,height = 8, dpi = 300, compression = "lzw")  
  
  
  return(output)
}

F.from.U=function(U) -log(1-U) 

#---6.  Create function for Kobe plot  ------------------------------------------------
kobePlot <- function(f.traj,b.traj,Years,Titl,Probs=NULL,txt.col='black',YrSize=3)
{
  dta=data.frame(x=b.traj,
                 y=f.traj,
                 yr=Years)%>%
    arrange(yr)
  Mx.F=max(2,max(dta$y,na.rm=T))
  Mx.B=max(2,max(dta$x,na.rm=T))
  
  kobe <-dta%>%
    ggplot(aes(x, y))+    
    scale_x_continuous(limits=c(0,Mx.B)) +
    scale_y_continuous(limits=c(0,Mx.F))+
    geom_rect(xmin = 1, xmax = Mx.B, ymin = 0, ymax = 1, fill = 'chartreuse3', alpha = 0.05) +
    geom_rect(xmin = 0, xmax = 1, ymin = 1, ymax = Mx.F, fill = 'brown1', alpha = 0.05) +
    geom_rect(xmin = 1, xmax = Mx.B, ymin = 1, ymax = Mx.F, fill = 'orange', alpha = 0.05) +
    geom_rect(xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = 'yellow', alpha = 0.05)
  if(!is.null(Probs))
  {
    kernelF <- gplots::ci2d(Probs$x, Probs$y, nbins = 151, factor = 1.5, 
                            ci.levels = c(0.5, 0.8, 0.75, 0.9, 0.95), show = "none")
    KernelD=rbind(kernelF$contours$"0.95"%>%mutate(CI='1',col='grey30'),
                  kernelF$contours$"0.8"%>%mutate(CI='2',col='grey50'),
                  kernelF$contours$"0.5"%>%mutate(CI='3',col='grey75'))
    kernels=KernelD%>%distinct(CI,col)%>%pull(col)
    names(kernels)=KernelD%>%distinct(CI,col)%>%pull(CI)
    
    Pr.d=data.frame(
      Prob=c(sum(ifelse(Probs$x > 1 & Probs$y < 1, 1, 0))/length(Probs$x)*100,
             sum(ifelse(Probs$x < 1 & Probs$y < 1, 1, 0))/length(Probs$x)*100,
             sum(ifelse(Probs$x > 1 & Probs$y > 1, 1, 0))/length(Probs$x)*100,
             sum(ifelse(Probs$x < 1 & Probs$y > 1, 1, 0))/length(Probs$x) * 100),
      col=c("green","yellow","orange","red"),
      x=rep(-10,4),  #dummy
      y=rep(-10,4))
    pr.ds=Pr.d%>%pull(col)
    names(pr.ds)=paste(round(Pr.d%>%pull(Prob),1),'%',sep='')

    
    
    kobe <-kobe +
      geom_polygon(data=KernelD,aes(x, y,fill=CI),size=1.25,alpha=0.5)+
      scale_fill_manual(labels=c("95%","80%","50%"),values = kernels)+
      geom_point(data=Pr.d,aes(x, y,color=col),alpha = 1,size=5)+
      scale_color_manual(labels=names(pr.ds),values = pr.ds)+
      labs(CI="CI", col="Prob.")
    

    
  }
  kobe <-kobe +
          geom_path(linetype = 2, size = 0.5,color='blue')+
          geom_point(size=2,color='blue')+
          geom_point(aes(x=dta[1,'x'],y=dta[1,'y']),size=4,shape=22,fill='white',alpha=.3)+
          geom_point(aes(x=dta[nrow(dta),'x'],y=dta[nrow(dta),'y']),size=4,shape=25,fill='white',alpha=.3)+      
          geom_text_repel(data=dta[1,],aes(x=x,y=y,label=yr),size=YrSize,color=txt.col)+
          geom_text_repel(data=dta[nrow(dta),],aes(x=x,y=y,label=yr),size=YrSize,color=txt.col)+
          xlab(expression(B/~B[MSY]))+ylab(expression(F/~F[MSY]))+
          labs(title = Titl)+
          theme_bw()%+replace% 
          theme(panel.grid.minor = element_blank(),
                axis.text = element_text(size=16),
                axis.title = element_text(size=20),
                plot.title = element_text(size=20,hjust=0),
                legend.text = element_text(size=15),
                legend.title = element_text(size=17),
                legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(-10,0,-10,-10))
  return(kobe)
}

#---7.  Source demography and steepness functions -------------------------------------------------------
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
                M=Rprior$M,
                G=Rprior$G,
                Input.pars=Rprior$Input.pars))
    
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
fn.display.priors=function(d,sp,XLAB,XLIM)
{
  dummy=lapply(d[sp],function(x) rnorm(1e4,x$mean,x$sd))
  dummy=lapply(dummy,function(x)data.frame(var=x))
  dummy=do.call(rbind,dummy)%>%
    tibble::rownames_to_column(var='Species')%>%
    mutate(Species=capitalize(gsub("\\..*","",Species)))
  p=dummy%>%
    ggplot(aes(x=var))+
    geom_density(aes(color=Species),size=1.5)+
    facet_wrap(~Species,scales='free_y')+
    xlab(XLAB)+ylab("Density")+
    theme_PA(axs.T.siz=22,axs.t.siz=14,strx.siz=16)+
    theme(legend.position = "none",
          plot.title =element_text(size=17))+
    scale_x_continuous(limits = XLIM)
  
  return(p)
}

#---8.  Import species-specific data -----
#note: this brings in any available data (cpue, abundance, selectivity, size composition, tagging, etc)
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
  if(length(files)>0)
  {
    if(length(files)>1)
    {
      files=sapply(files, fread, data.table=FALSE)
    }
    else
    {
      files=list(read.csv(files))
    }
    names(files)=file.names
    Species.data[[s]]=files
  }
  rm(files)
}


#---9. Input parameters -----
  #set up list
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

  #add life history parameters
for(l in 1:N.sp)
{
  print(paste("---------Set up input parameters for --",names(List.sp)[l]))
  
  LH=LH.data%>%filter(SPECIES==List.sp[[l]]$Species)
  
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
                           NsimSS=1e4,                        #demography
                           r.prior="USER",                    #demography
                           r.prior2=NA,                       #demography uniform
                           a_FL.to.TL=LH$a_FL.to.TL,          # FL to TL
                           b_FL.to.TL=LH$b_FL.to.TL
  )              
  
}

  #Export table of life history parameters for use in RAR
Rar.path=paste(handl_OneDrive('Reports/RARs'), AssessYr,sep="/")
if(!dir.exists(Rar.path))dir.create(Rar.path)
if(First.run=="YES")
{
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
 # fn.word.table(TBL=TabL,Doc.nm="Table 1. Life history pars")
  write.csv(TabL,"Table 1. Life history pars.csv",row.names = F)
}

  #Export list of species assessed for use in RAR
if(First.run=="YES")
{
  write.csv(data.frame(Species=capitalize(RAR.species)),
            paste(Rar.path,"Assessed_species.csv",sep='/'),row.names = F)
}




#---10. Create list of species outputs ----- 
Lista.sp.outputs=list(Other.species,names(Indicator.species))
names(Lista.sp.outputs)=c('Other.sp','Indicator.sp')
if(!is.null(additional.sp))
{
  Lista.sp.outputs[[3]]=additional.sp
  names(Lista.sp.outputs)[3]='additional.sp'
}

#---11. Calculate r prior -----  
  #calculate prior
store.species.r=vector('list',N.sp)
names(store.species.r)=Keep.species
store.species.M=store.species.G=store.species.r
if(do.r.prior)
{
  system.time({for(l in 1:N.sp)   #takes 0.002 sec per iteration (NsimSS) per species
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
    r.prior.dist=with(List.sp[[l]],fun.rprior.dist(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,
                                                   K.sd=Growth.F$k.sd,LINF.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL,k.Linf.cor,
                                                   Amax=Max.age.F,
                                                   MAT=unlist(Age.50.mat),FecunditY=Fecundity,Cycle=Breed.cycle,
                                                   BWT=BwT,AWT=AwT,LO=Lzero*a_FL.to.TL+b_FL.to.TL)) #size vars as TL
    #export life history parameter distributions
    Nms=names(r.prior.dist$Input.pars[[1]])
    LH.d=matrix(unlist(list.flatten(r.prior.dist$Input.pars)),nrow=List.sp[[l]]$NsimSS,ncol=length(Nms),byrow = T)
    colnames(LH.d)=Nms
    LH.d%>%
      data.frame%>%
      gather(Variable, Value)%>%
      mutate(Variable=case_when(Variable == 'Max.age'~ 'Maximum age',
                                Variable == 'Age.mat'~ 'Age at 50% maturity',
                                Variable == 'Meanfec'~ 'Fecundity',
                                Variable == 'Reprod_cycle'~ 'Reprod cycle length',
                                Variable == 'Linf'~ 'Growth - Linf',
                                Variable == 'k'~ 'Growth - k',
                                TRUE ~ Variable))%>%
      ggplot(aes(x=Value))+
      geom_histogram(fill='#00BA38',color='black')+
      facet_wrap(~Variable,scales='free')+
      ylab('Frequency')+
      scale_x_continuous(breaks= pretty_breaks())+
      theme_PA()
    ggsave('Life_history_priors.tiff',width = 8,height = 6, dpi = 300, compression = "lzw")
    rm(LH.d)
    
    #export r
    out.r=data.frame(shape=r.prior.dist$shape,rate=r.prior.dist$rate,
                     mean=r.prior.dist$mean,sd=r.prior.dist$sd)
    write.csv(out.r,'r.prior.csv',row.names = F)
    store.species.r[[l]]=r.prior.dist
    
    #export G
    out.G=data.frame(mean=mean(r.prior.dist$G),sd=sd(r.prior.dist$G))
    write.csv(out.G,'G.prior.csv',row.names = F)
    store.species.G[[l]]=out.G
    
    #export M
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
    store.species.M[[l]]=out.M
    write.csv(out.M,"M.csv",row.names=FALSE)
    
    rm(r.prior.dist,M.averaging,RESAMP,linear.fec,out.M)
  }})
}

for(l in 1:N.sp) 
{
  hndl.dummy=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                   AssessYr,"/demography",sep='')
  store.species.r[[l]]=read.csv(paste(hndl.dummy,"/r.prior.csv",sep=''))
  store.species.G[[l]]=read.csv(paste(hndl.dummy,"/G.prior.csv",sep=''))
  store.species.M[[l]]=read.csv(paste(hndl.dummy,"/M.csv",sep=''))
  rm(hndl.dummy)
}

  #display priors   
for(l in 1:length(Lista.sp.outputs))
{
  fn.display.priors(d=store.species.r,
                    sp=Lista.sp.outputs[[l]],
                    XLAB=expression(paste(plain("Intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),
                    XLIM=c(0,NA))
  ggsave(paste(Rar.path,'/Prior_r_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = 12,height = 10,compression = "lzw")
}

  #compare M and r
if(do.r.prior)
{
  Omit.these=c("Great hammerhead","Scalloped hammerhead","Smooth hammerhead",
               "Tiger shark","Dwarf sawfish","Freshwater sawfish")
  CompR=data.frame(Name=names(store.species.r),r=NA,M=NA)
  for(l in 1:N.sp)
  {
    CompR$r[l]=store.species.r[[l]]$mean
    CompR$M[l]=mean(unlist(store.species.M[[l]]),na.rm=T)
  }
  COL=rgb(.1,.2,.8,alpha=.45)
  CompR=CompR%>%
    arrange(r)%>%
    mutate(Name=capitalize(Name))
  my_formula=y ~ x
  
  p=CompR%>%
    ggplot(aes(M, r, label = Name)) +
    geom_point(shape = 21, size = 5,fill=COL) + 
    geom_smooth(method = "lm", data = CompR%>%filter(!Name%in%Omit.these),
                se = F, fullrange = TRUE,colour="red")+
    stat_poly_eq(aes(label =  paste(stat(eq.label),stat(adj.rr.label),stat(p.value.label),
                                    sep = "*\", \"*")),
                 formula = my_formula, parse = TRUE,
                 label.y = "top", label.x = "right", size = 4)+
    geom_text_repel(segment.colour='black',col='black',box.padding = 0.5) + 
    scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))+ 
    theme_PA(axs.T.siz=14,axs.t.siz=12)+
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  p
  ggsave(handl_OneDrive('Analyses/Population dynamics/M_vs_r.tiff'), 
         width = 8,height = 6, dpi = 300, compression = "lzw")
  
}

#---12. Assign Resilience -----------------------------------------------------------------------
RESILIENCE=vector('list',N.sp)
names(RESILIENCE)=names(List.sp)
for(r in 1:length(RESILIENCE)) RESILIENCE[[r]]=Res.fn(store.species.r[[r]]$mean,Def="Haddon")


#---13. Extract selectivity at age and at size -----------------------------------------------------------------------
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
  if('gillnet.selectivity_len.age'%in%names(Species.data[[l]]) | names(Species.data)[l]%in%Sel.equivalence$Name)
  {
    #1. Read in selectivity at age data
    if('gillnet.selectivity_len.age'%in%names(Species.data[[l]]))
    {
      GN.sel.at.age=Species.data[[l]]$gillnet.selectivity_len.age%>%mutate(type='Species')
      GN.sel.at.totalength=Species.data[[l]]$gillnet.selectivity%>%mutate(type='Species')
    }else
    {
      #allocate  selectivity from family
      this.sel=Sel.equivalence%>%filter(Name==names(Species.data)[l])
      temp.wd=paste(HandL,this.sel$Equivalence,sep='')
      GN.sel.at.age=read.csv(paste(temp.wd,'/',this.sel$Equivalence,'_Gillnet.selectivity_len.age.csv',sep=''))%>%mutate(type='Family')
      GN.sel.at.totalength=read.csv(paste(temp.wd,'/',this.sel$Equivalence,'_Gillnet.selectivity.csv',sep=''))%>%mutate(type='Family')
      
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
    Selectivity.at.age[[l]]=GN.sel.at.age[,c('TL','Age','Sel.combined','X16.5','X17.8','type')]
    
    GN.sel.at.totalength=GN.sel.at.totalength%>%
      mutate(Sum.sel=X16.5+X17.8,
             Sel.combined=Sum.sel/max(Sum.sel),
             Sel.combined=Sel.combined/max(Sel.combined),
             TL=TL.mm)
    Selectivity.at.totalength[[l]]=GN.sel.at.totalength[,c('TL','Sel.combined','X16.5','X17.8','type')]
  }
}

#display selectivities    
if(First.run=="YES")
{
  N.sp.with.sel=which(sapply(Selectivity.at.age,function(x) !is.null(x)),TRUE)
  fn.fig(handl_OneDrive('Analyses\\Population dynamics\\growth.and.selectivity'),2400,2000) 
  smart.par(n.plots=length(N.sp.with.sel),MAR=c(2,3,1,1),OMA=c(2.5,1,.05,2.5),MGP=c(1.8,.5,0))
  par(cex.lab=1.5,las=1)
  for(l in N.sp.with.sel)
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
  smart.par(n.plots=length(N.sp.with.sel),MAR=c(2,3,1,1),OMA=c(2.5,1,.05,2.5),MGP=c(1.8,.5,0))
  par(cex.lab=1.5,las=1)
  for(l in N.sp.with.sel)
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

#---14. Calculate Steepness ----------------------------------------------------------------------- 
  #calculate prior
store.species.steepness=vector('list',N.sp)
names(store.species.steepness)=Keep.species
store.species.alpha=store.species.steepness
if(do.steepness)
{
  system.time(for(l in 1: N.sp)  #takes 0.003 sec per iteration (NsimSS) iteration per species
  {
    print(paste("steepness ","--",List.sp[[l]]$Name))
    PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
               capitalize(List.sp[[l]]$Name),"/",AssessYr,"/steepness",sep='')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
    setwd(PATH)
    
    SEL=Selectivity.at.age[[l]]$Sel.combined  #selectivity not used in h calculation as F.mult =0 so set to dummy if sel not available
    if(is.null(SEL)) SEL=rep(0,List.sp[[l]]$Max.age.F[2])
    
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
    
    
    #export h and alpha
    out.h=data.frame(shape=steepNs$shape,rate=steepNs$rate,
                     mean=steepNs$mean,sd=steepNs$sd)
    write.csv(out.h,'h.prior.csv',row.names = F) 
    store.species.steepness[[l]]=out.h
    
    write.csv(steepNs$Alpha,'Alpha.csv',row.names = F) 
    store.species.alpha[[l]]=steepNs$Alpha
    
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
}

for(l in 1: N.sp)
{
  store.species.steepness[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                              AssessYr,"/steepness/h.prior.csv",sep=''))
  store.species.alpha[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                          AssessYr,"/steepness/Alpha.csv",sep=''))
}

  #display priors
for(l in 1:length(Lista.sp.outputs))
{
  fn.display.priors(d=store.species.steepness,
                    sp=Lista.sp.outputs[[l]],
                    XLAB="Steepness (h)",
                    XLIM=c(0.2,NA))
  ggsave(paste(Rar.path,'/Prior_steepness_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = 12,height = 10,compression = "lzw")
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


#Recalculate steepness for hammerheads and tiger using linear model of r on h
#note: h values deemed to high (following E Cortes discussion)
store.species.steepness.S2=fn.get.stuff.from.list(store.species.steepness,"mean")
dis.sp.h=c("great hammerhead","scalloped hammerhead","smooth hammerhead",
           "tiger shark")
if("dwarf sawfish" %in% Keep.species) dis.sp.h=c(dis.sp.h,"dwarf sawfish","freshwater sawfish")

for(s in 1:length(dis.sp.h))
{
  id=match(dis.sp.h[s],names(store.species.steepness))
  store.species.steepness.S2[[id]]=store.species.r[[id]]$mean*0.643+0.229
}




#---15. Calculate scaler for Fmsy.to.M relationship (Fmsy= scaler x M) for DBSRA ----------------------------------------------------------------------- 
Cortes.Brooks.2018=function(alpha)  #source: Cortes & Brooks 2018
{
  if(alpha<=2.67)           Fmsy.to.M.scaler=0.2
  if(alpha>2.67 & alpha<=6) Fmsy.to.M.scaler=0.5
  if(alpha>6)               Fmsy.to.M.scaler=0.8
  
  return(Fmsy.to.M.scaler)
}
Fmsy.M.scaler=vector('list',N.sp)
names(Fmsy.M.scaler)=Keep.species

for(l in 1:N.sp) Fmsy.M.scaler[[l]]=Cortes.Brooks.2018(alpha=median(unlist(store.species.alpha[[l]])))

#reset dis.sp.h consistently with h resetting
for(l in 1:length(dis.sp.h))
{
  s=match(dis.sp.h[l],names(Fmsy.M.scaler))
  Fmsy.M.scaler[[s]]=0.5
}
#---16. Modelling arguments and pin files -----
#note: For integrated model, S1 and S2 calculates pars in normal space but same order magnitude
#       Other scenarios all pars in log.
#       ln_RZERO is in 1,000 individuals so do 10 times the largest catch divided by 
#       average weight and divided by 1,000. Best units to work in are 1,000 individuals for 
#       numbers, catch in tonnes and length-weight in kg as all cancels out and predicted biomasses
#       end up being in tonnes

k.fun.low=function(KTCH) max(KTCH,na.rm=T)  #K boundaries
times.max.ktch=50         #50 times: Andrade 2017
k.fun.up=function(KTCH) max(KTCH,na.rm=T)*times.max.ktch   
fn.mtch=function(WHAT,NMS) match(WHAT,names(NMS))
Q_phz=c("lnq","lnq2","log_Qdaily")                           
Zns.par.phz=c("lnR_prop_west","lnR_prop_zn1")
MOv.par.phz=c("log_p11","log_p22","log_p21","log_p33")

#species-specific proxy to Bmsy.K based on Cortes et al 2012
R=function(r,G) 0.633-0.187*log(r*G) 
BmsyK.species=data.frame(Species=names(store.species.r),r=NA,G=NA)
for(l in 1:N.sp)
{
  BmsyK.species$r[l]=store.species.r[[l]]$mean
  BmsyK.species$G[l]=store.species.G[[l]]$mean
}
BmsyK.species=BmsyK.species%>%
  mutate(BmsyK=R(r,G))
plot.bmsyK=FALSE
if(plot.bmsyK)
{
  BmsyK.species%>%
    ggplot(aes(r,BmsyK,color=G))+
    geom_point(shape=19,size=2)+
    scale_color_gradient(low = "yellow", high = "red")+
    ylim(0,1)
}

for(l in 1:N.sp)
{
  print(paste("---------Modelling arguments for --",names(List.sp)[l]))
  
  #... Surplus production arguments
  #note: Only fitting species with species-specific abundance time series 
  #Assumptions: negligible exploitation at start of time series
  
  List.sp[[l]]=list.append(List.sp[[l]],
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
                           N.monte=1e3,
                           
                           #Initial par value  "whiskery shark" "dusky shark"   "gummy shark" "sandbar shark"
                           Init.r=with(List.sp[[l]],
                                       case_when(
                                         store.species.r[[l]]$mean< .1 ~.05,
                                         store.species.r[[l]]$mean>= .1 & store.species.r[[l]]$mean<.2 ~.15,
                                         store.species.r[[l]]$mean>= .2 & store.species.r[[l]]$mean<.4 ~.3,
                                         store.species.r[[l]]$mean>= .4 & store.species.r[[l]]$mean<.6 ~.5,
                                         store.species.r[[l]]$mean>= .6  ~.7,
                                         TRUE~NA_real_)),
                           
                           #maximum acceptable CV for cpue series  
                           MAX.CV=0.5,   
                           
                           #define which optimisation method to use
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
  
  #... Catch-only arguments and sensitivity tests
  List.sp[[l]]$STARTBIO=c(List.sp[[l]]$B.init*.95,List.sp[[l]]$B.init)   #low depletion because starting time series prior to any fishing
  List.sp[[l]]$FINALBIO=c(.15,1)       #highly uncertain final depletion
  List.sp[[l]]$Do.sim.test="NO"  #simulation test Catch-MSY for small and large catches
  Proc.error=1e-05                  #default process error 
  Proc.error.JABBA=5e-02            #process error for JABBA  (Winker et al 2018 School shark)
  Proc.error.cmsy=Proc.error.JABBA   #bigger proc error for sensitivity test
  BmsyK.Cortes=BmsyK.species%>%filter(Species==names(List.sp)[l])%>%pull(BmsyK)

  #..DBSRA scenarios
  AgeMat=List.sp[[l]]$Age.50.mat[1]
  Mmean=mean(apply(store.species.M[[l]],2,mean,na.rm=T))
  Msd=mean(apply(store.species.M[[l]],2,sd,na.rm=T))
  ktch=ktch.combined%>%
    filter(Name==List.sp[[l]]$Name)  
  Klow=k.fun.low(ktch$Tonnes) 
  Kup=k.fun.up(ktch$Tonnes)  
  fmsy.m=0.41                 # Zhou et al 2012 fmsy/M=0.41 for chondrichthyans
  fmsy.m2=Fmsy.M.scaler[[l]]  #Cortes & Brooks 2018

  #..CMSY scenarios
  r.prob.min=0.01   #quantile probs for defining r range
  r.prob.max=0.99

    # some species-specific input values for CMSY Process error scenario to allow convergence
  Kup.spcfik=Kup
  if(names(List.sp)[l]%in%c("great hammerhead")) Kup.spcfik=1000
  if(names(List.sp)[l]%in%c("smooth hammerhead")) Kup.spcfik=1500
  if(names(List.sp)[l]%in%c("green sawfish")) Kup.spcfik=600
  if(names(List.sp)[l]%in%c("narrow sawfish")) Kup.spcfik=350
  if(names(List.sp)[l]%in%c("gummy shark")) Kup.spcfik=6000
  if(names(List.sp)[l]%in%c("tiger shark")) Kup.spcfik=8000
    
  Klow.spcfik=Klow
  if(names(List.sp)[l]%in%c("green sawfish")) Klow.spcfik=200
  if(names(List.sp)[l]%in%c("gummy shark")) Klow.spcfik=2500
  if(names(List.sp)[l]%in%c("tiger shark")) Klow.spcfik=2000

  
  if(names(List.sp)[l]%in%c("gummy shark","narrow sawfish")) Proc.error.cmsy=1e-3  #Proc.error.JABBA does not converge
  if(names(List.sp)[l]%in%c("tiger shark")) Proc.error.cmsy=5e-3
  if(names(List.sp)[l]%in%c("milk shark")) Proc.error.cmsy=Proc.error   #not converging otherwise
  
  
  # DBSRA not converging if set too low
  Kup.S3=Klow*(times.max.ktch/2)
  if(names(List.sp)[l]%in%c("dusky shark")) Kup.S3=Kup*.75  
  
  
  #..JABBA scenarios
  r.CV.multiplier=1
  Ktch.CV=0.05   #Winker et al 2019 set it at 0.2 for uncertain reconstructed school shark catches
  
  
  #..Put all in list  
  ensims<-3e4  #Winkner et al 2019
  #if(names(List.sp)[l]=='milk shark') ensims=1e5
  
  List.sp[[l]]$Sens.test=list(
    DBSRA=data.frame(Scenario=paste("S",1:4,sep=''),
                     Sims=rep(3e4,4), 
                     AgeMat=rep(AgeMat,4),
                     Mmean=rep(Mmean,4),
                     Msd=rep(Msd,4),
                     Klow=rep(Klow,4),
                     Kup=c(rep(Kup,2),Kup.S3,Kup),
                     fmsy.m=c(fmsy.m,fmsy.m2,rep(fmsy.m,2)),    
                     bmsyk=c(rep(0.5,3),BmsyK.Cortes)),    #Winker et al 2020 Bmsy/K=0.55 for mako shark
    
    CMSY=data.frame(Scenario=paste("S",1:3,sep=''),
                    Sims=rep(ensims,3), 
                    r.prob.min=rep(r.prob.min,3),  
                    r.prob.max=rep(r.prob.max,3),
                    Klow=c(Klow.spcfik,Klow.spcfik,Klow.spcfik),
                    Kup=c(Kup,Kup.spcfik,Klow*(times.max.ktch/2)),
                    Proc.error=c(Proc.error,Proc.error.cmsy,Proc.error)),
    
    JABBA=data.frame(Scenario=paste("S",1:4,sep=''),
                     Sims=c(2e5,rep(3e4,3)),
                     r.CV.multiplier=rep(r.CV.multiplier,4),
                     K.mean=rep(Klow*20,4),       #Winker et al 2019 mean k set at 20 times max catch
                     K.CV=rep(2,4),
                     Proc.error=c(Proc.error.JABBA,Proc.error,rep(Proc.error.JABBA,2)),
                     Ktch.CV=c(rep(Ktch.CV,2),1e-4,Ktch.CV),
                     bmsyk=c(rep(0.5,3),BmsyK.Cortes))      #Winker et al 2020 Bmsy/K=0.55 for mako shark
  )
  rm(ensims)
  
  
  
  # Arguments for indicator species only
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
      fn.word.table(TBL=Tabla.scen.show,Doc.nm="Model scenarios")
      
      
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
  
  
}

#---17. Export all available input data to each species assessment folder ----- 
if(First.run=="YES")
{
  source(handl_OneDrive("Analyses/Population dynamics/Git_Stock.assessments/Organise data.R"))
  for(l in 1:N.sp)
  {
    print(paste("----Export data inputs for ",names(List.sp)[l]))
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

#---18. Display catches by fishery ----
Tot.ktch=KtCh %>%      
  mutate(
    Type = case_when(
      FishCubeCode=='WRL'~'WRL',
      FishCubeCode=='WTB'~'WTB',
      FishCubeCode=='TEP'~'TEP',
      FishCubeCode=='Recreational'~'Recreational',
      FishCubeCode=='SA MSF'~'SA MSF',
      FishCubeCode=='NSW fisheries'~'NSW fisheries',
      FishCubeCode=='NT'~'NT',
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

for(l in 1:length(Lista.sp.outputs))
{
  SIZ=2
  WID=14
  HEI=10
  if(length(Lista.sp.outputs[[l]])>8)
  {
    SIZ=1.5
  }
  
  if(length(Lista.sp.outputs[[l]])<4)
  {
    SIZ=3
  }
  
  Tot.ktch%>%
    filter(Name%in%Lista.sp.outputs[[l]])%>%
    group_by(Type,finyear,Name)%>%
    summarise(Tot=sum(LIVEWT.c,na.rm=T)/unitS)%>%
    filter(!is.na(Tot))%>%
    filter(Tot>0)%>%
    mutate(Name=capitalize(Name))%>%
    ggplot(aes(finyear,Tot,color=Type))+
    geom_point(size=SIZ)+
    geom_line(linetype=2,alpha=0.4)+
    facet_wrap(~Name,scales='free')+
    theme_PA(strx.siz=15,leg.siz=14,axs.t.siz=14,axs.T.siz=18)+
    ylab("Total catch (tonnes)")+xlab("Financial year")+
    theme(legend.position="top",
          legend.title = element_blank(),
          legend.key=element_blank(),
          plot.margin = margin(0.1,.5,0.1,0.1, "cm"))+
    guides(colour = guide_legend(override.aes = list(size=5,linetype = 0)))
  
  ggsave(paste(Rar.path,'/Catch_all_species_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = WID,height = HEI,compression = "lzw")
  
}

#---19. Catch-only assessments. Run models --------------------------------------

fn.prior=function(N=1e4,d,MAX=NULL)
{
  if(d$dist=='unif')
  {
    out=runif(n=N, min=d$low, max=d$up) 
  }
  if(d$dist=='lnorm')
  {
    out=rlnorm(n=N, meanlog=d$mean, sdlog=d$sd)
    if(!is.null(MAX)) out=subset(out,out<=MAX)
  } 
  if(d$dist=='beta')
  {
    shape.pars=get_beta(d$mean,d$sd)
    out=rbeta(n=N, shape1=shape.pars[1], shape2=shape.pars[2])
  }
  return(out)
}
fn.show.density=function(d,NCOL)
{
  d%>%
    ggplot(aes(x=Value,fill=Distribuion))+
    geom_density(position="identity",size=1.15,alpha=0.65)+
    facet_wrap(~Parameter,scales='free',ncol=NCOL)+
    theme_PA(strx.siz=18,leg.siz=18,axs.t.siz=16,axs.T.siz=26)+
    theme(legend.position="top",
          legend.title = element_blank(),
          legend.key=element_blank(),
          title=element_text(size=12))+
    ylab("Density")+xlab("Parameter value")
}
add.probs=function(DAT,id.yr,B.threshold,plot.ranges=FALSE) #Reference points probability  
{
  DAT=DAT[,id.yr]
  B.target=Tar.prop.bmsny*B.threshold
  B.limit=Lim.prop.bmsy*B.threshold
  f=ecdf(DAT)
  P.below.target=f(B.target)
  P.below.threshold=f(B.threshold)
  P.below.limit=f(B.limit)
  P.above.target=1-P.below.target
  P.above.threshold=1-P.below.threshold
  P.above.limit=1-P.below.limit
  P.between.thre.tar=P.below.target-P.below.threshold
  P.between.lim.thre=P.below.threshold-P.below.limit
  if(plot.ranges)
  {
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
    
  }
  
  return(list(probs=data.frame(Range=c('<lim','lim–thr',
                                       'thr–tar','>tar'),
                               Probability=round(c(P.below.limit,P.between.lim.thre,
                                                   P.between.thre.tar,P.above.target),3)),
              Reference.points=data.frame(Rf.pt=c('Target','Threshold','Limit'),
                                          Value=c(B.target,B.threshold,B.limit))))
}
fn.ktch.only.get.timeseries=function(d,mods,Type,add.50=FALSE,scen,Katch)
{
  dd=d$output
  
  if(mods=='DBSRA')
  {
    Years=dd$Years
    if(Type=='Depletion')
    {
      d1=dd$Depletion.traj[1:length(Years)]
      Probs=add.probs(DAT=d1,
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=dd$Estimates[Biomass.threshold,'Median (ll=1)']/dd$Estimates['K','Median (ll=1)'])
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
    }
    if(Type=='F.series')d1=dd$F.series[1:length(Years)]
    if(Type=='B.Bmsy') d1=dd$B.Bmsy[1:length(Years)]
    if(Type=='F.Fmsy') d1=dd$F.Fmsy[1:length(Years)]
    Dat=data.frame(year=Years,
                   median=apply(d1,2,median),
                   upper.95=apply(d1,2,function(x)quantile(x,probs=0.975,na.rm=T)),
                   lower.95=apply(d1,2,function(x)quantile(x,probs=0.025,na.rm=T)))%>%
      mutate(Model=mods,
             Catch=Katch)
    if(add.50)
    {
      Dat=Dat%>%
        mutate(upper.50=apply(d1,2,function(x)quantile(x,probs=0.75,na.rm=T)),
               lower.50=apply(d1,2,function(x)quantile(x,probs=0.25,na.rm=T)))
    }
  }
  
  if(mods=='CMSY')
  {
    Years=dd$Years
    if(Type=='Depletion')
    {
      d1=dd$Depletion.traj[1:length(Years)]
      Probs=add.probs(DAT=d1,
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=median(dd$Bmsy)/dd$Statistics$output['K','50%'])
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
    }
    
    if(Type=='F.series')  d1=dd$F.traj[,1:length(Years)]
    if(Type=='B.Bmsy') d1=dd$B.Bmsy[1:length(Years)]
    if(Type=='F.Fmsy') d1=dd$F.Fmsy[,1:length(Years)]
    
    Dat=data.frame(year=as.numeric(Years),
                   median=apply(d1,2,median),
                   upper.95=apply(d1,2,function(x)quantile(x,probs=0.975,na.rm=T)),  
                   lower.95=apply(d1,2,function(x)quantile(x,probs=0.025,na.rm=T)))%>%
      mutate(Model=mods,
             Catch=Katch) 
    if(add.50)
    {
      Dat=Dat%>%
        mutate(upper.50=apply(d1,2,function(x)quantile(x,probs=0.75,na.rm=T)),
               lower.50=apply(d1,2,function(x)quantile(x,probs=0.25,na.rm=T)))
    }
  }
  
  if(mods=='JABBA')
  {
    Years=dd$yr
    
    if(Type=='Depletion') 
    {
      K=dd$pars
      K=K[match("K",rownames(K)),"Median"]
      d1=data.frame(mu=apply(dd$posteriors$P,2,median,na.rm=T),
                    lci=apply(dd$posteriors$P,2,function(x) quantile(x,probs=0.025,na.rm=T)),
                    uci=apply(dd$posteriors$P,2,function(x) quantile(x,probs=0.975,na.rm=T)))
      Probs=add.probs(DAT=sweep(dd$posteriors$SB,1,dd$posteriors$K,'/'),
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=median(dd$posteriors$SBmsy)/K)
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
    }
    if(Type=='F.series') d1=data.frame(dd$timeseries[, , "F"])
    if(Type=='B.Bmsy') d1=data.frame(dd$timeseries[, , "BBmsy"])
    if(Type=='F.Fmsy') d1=data.frame(dd$timeseries[, , "FFmsy"])
    
    Dat=data.frame(year=as.numeric(Years),
                   median=d1$mu,
                   upper.95=d1$uci,  
                   lower.95=d1$lci)%>%
      mutate(Model=mods,
             Catch=Katch) 
    
    if(add.50)
    {
      Dat=Dat%>%
        mutate(upper.50=NA,
               lower.50=NA)
    }
    
  }
  
  Dat=Dat%>%mutate(Scenario=scen)
  
  if(Type=='Depletion')
  {
    return(list(Dat=Dat,Probs=Probs))
  }else
  {
    return(list(Dat=Dat))
  }
  
}

  #19.1 Execute each COM
n.catch.only=length(catch.only)  
Catch_only=vector('list',n.catch.only)
names(Catch_only)=catch.only

for(w in 1:length(Catch_only))
{
  #11.1 DBSRA assessment (Dick and MAcCall (2011))
  # summary of method: http://toolbox.frdc.com.au/wp-content/uploads/sites/19/2020/07/DBSRA3.html
  if(names(Catch_only)[w]=="DBSRA")
  {
    dummy.store=vector('list',N.sp)     #takes 0.011 secs per iteration per species per scenario
    names(dummy.store)=Keep.species
    dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
    dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.accept.rate=
      dummy.store.ensemble=dummy.store
    for(i in 1:length(dummy.store))  
    {
      ktch=ktch.combined%>%
        filter(Name==names(dummy.store)[i])
      this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(List.sp[[i]]$Name),"/",AssessYr,"/DBSRA",sep='')
      if(!dir.exists(this.wd))dir.create(this.wd)
      
      Scens=List.sp[[i]]$Sens.test$DBSRA%>%
                mutate(Species=capitalize(names(dummy.store)[i]))
      Store.sens=vector('list',nrow(Scens))
      names(Store.sens)=Scens$Scenario
      this.wd1=this.wd
      
      Out.Scens=Scens
      Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
        Out.B.Bmsy=Out.F.Fmsy=Out.accept.rate=vector('list',length(Store.sens))
      
      for(s in 1:length(Store.sens))
      {
        print(paste("___________","DBSRA Scenario",Scens$Scenario[s],"___________",names(dummy.store)[i]))
        
        this.wd=paste(this.wd1,names(Store.sens)[s],sep='/')
        if(!dir.exists(this.wd))dir.create(this.wd)
        
        AgeMat=Scens$AgeMat[s]
        
        Mmean=Scens$Mmean[s]  
        Msd=Scens$Msd[s]
        Mcv=Msd/Mmean
        M.dist="lnorm"
        
        Klow=Scens$Klow[s]
        Kup=Scens$Kup[s]
        
        fmsy.m=Scens$fmsy.m[s]
        fmsym.dist="lnorm"
        fmsym.CV=0.1
        
        bmsyk.mean=Scens$bmsyk[s]
        bmsyk.dist="beta"
        bmsyk.CV=0.1
        
        btk.dist="unif"
        Btklow=List.sp[[i]]$FINALBIO[1]
        Btkup=List.sp[[i]]$FINALBIO[2]
        
        b1k.dist="unif"
        b1k.low=List.sp[[i]]$STARTBIO[1]
        b1k.up=List.sp[[i]]$STARTBIO[2]

        #Run model
        Store.sens[[s]]=apply.DBSRA(year=ktch$finyear,
                        catch=ktch$Tonnes,
                        catchCV=NULL,  
                        catargs=list(dist="none",low=0,up=Inf,unit="MT"),  #catch CV not available
                        agemat=AgeMat,
                        k=list(low=Klow,up=Kup,tol=0.01,permax=1000),
                        b1k=list(dist=b1k.dist,low=b1k.low,up=b1k.up,mean=1,sd=0.1),  #mean and sd not used if 'unif'
                        btk=list(dist="unif",low=Btklow,up=Btkup,mean=1,sd=0.1,refyr=max(ktch$finyear)),  #reference year
                        fmsym=list(dist="lnorm",low=0.1,up=2,mean=log(fmsy.m),sd=fmsym.CV), # Cortes & Brooks 2018. Low and up not used if 'lnorm'  
                        bmsyk=list(dist="beta",low=0.05,up=0.95,mean=bmsyk.mean,sd=bmsyk.CV),  
                        M=list(dist="lnorm",low=0.001,up=1,mean=log(Mmean),sd=Mcv),
                        graph=c(13,14),
                        nsims=Scens$Sims[s],
                        grout=1,
                        WD='C:/DummyDBSRA',  #output to dummy folder because OneDrive backup stuff things up
                        outfile="Appendix_fit")
        
        #Store acceptance rate
        Accept.tab=table(Store.sens[[s]]$output$Values$ll)
        Out.accept.rate[[s]]=data.frame(
                        Acceptance=round(100*Accept.tab[2]/sum(Accept.tab),2),
                        Scenario=Scens$Scenario[s])

        #Store scenarios
        Out.Scens$M.dist=M.dist
        Out.Scens$M.CV=Mcv
        Out.Scens$fmsym.dist=fmsym.dist
        Out.Scens$fmsym.CV=fmsym.CV
        Out.Scens$bmsyk.dist=bmsyk.dist
        Out.Scens$bmsyk.mean=Out.Scens$bmsyk
        Out.Scens$bmsyk.CV=bmsyk.CV
        Out.Scens$b1k.dist=b1k.dist
        Out.Scens$b1k.low=b1k.low
        Out.Scens$b1k.up=b1k.up
        Out.Scens$btk.dist=btk.dist
        Out.Scens$btk.low=Btklow
        Out.Scens$btk.up=Btkup
        
        #Store estimates
        d1=Store.sens[[s]]$output$Parameters
        d1=d1[,grep(paste(c('Median','2.5%','97.5%'),collapse='|'),names(d1))]
        names(d1)=c("Median","Lower.95","Upper.95")
        d2=Store.sens[[s]]$output$Estimates
        d2=d2[,grep(paste(c('Median','2.5%','97.5%'),collapse='|'),names(d2))]
        names(d2)=c("Median","Lower.95","Upper.95")
        d1=rbind(d2,d1)
        d1=d1%>%
          mutate(Parameter=rownames(d1),
                 Model='DBSRA',
                 Scenario=names(Store.sens)[s])%>%
          relocate(Model,Scenario,Parameter,Lower.95)
        Out.estimates[[s]]=d1
        
        #Store trajectories 
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                          mods=names(Catch_only)[w],
                                          Type='Depletion',
                                          scen=Scens$Scenario[s],
                                          Katch=ktch$Tonnes)
        Out.rel.biom[[s]]=dummy$Dat
        Out.probs.rel.biom[[s]]=dummy$Probs
        
        
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                    mods=names(Catch_only)[w],
                                    Type='F.series',
                                    scen=Scens$Scenario[s],
                                    Katch=ktch$Tonnes)
        Out.f.series[[s]]=dummy$Dat
        
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                    mods=names(Catch_only)[w],
                                    Type='B.Bmsy',
                                    scen=Scens$Scenario[s],
                                    Katch=ktch$Tonnes)
        Out.B.Bmsy[[s]]=dummy$Dat
        
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                    mods=names(Catch_only)[w],
                                    Type='F.Fmsy',
                                    scen=Scens$Scenario[s],
                                    Katch=ktch$Tonnes)
        Out.F.Fmsy[[s]]=dummy$Dat
        rm(dummy)
        
        #Display Priors vs Posteriors for base case scenario (S1)
        if(Scens$Scenario[s]=='S1')
        {
          par.list=c('m','fmsym','b1k','btk','bmsyk')
          dummy.prior=Store.sens[[s]]$input
          names(dummy.prior)=tolower(names(dummy.prior))
          dummy.post=Store.sens[[s]]$output$Values%>%filter(ll==1) #accepted runs
          colnames(dummy.post)=tolower(colnames(dummy.post))
          out=vector('list',length(par.list))
          for(p in 1:length(par.list))
          {
            out[[p]]=rbind(data.frame(Distribuion="Prior",
                                      Value=fn.prior(d=dummy.prior[[par.list[p]]])),
                           data.frame(Distribuion="Posterior",
                                      Value=dummy.post[[par.list[p]]]))%>%
              mutate(Parameter=par.list[p])
          }
          rm(dummy.prior,dummy.post)
          fn.show.density(d=do.call(rbind,out),NCOL=2)
          ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                       capitalize(names(dummy.store)[i]),"/",AssessYr,"/",'DBSRA',"/Prior.and.posterior.tiff",sep=''),
                 width = 12,height = 14, dpi = 300, compression = "lzw")
        }
        
        #Store quantities for ensemble model 
        if(Scens$Scenario[s]=='S1')
        {
          Years=Store.sens[[s]]$output$Years
          dummy.store.ensemble[[i]]=list(Depletion=Store.sens[[s]]$output$Depletion.traj[1:length(Years)],
                                         B.Bmsy=Store.sens[[s]]$output$B.Bmsy[1:length(Years)],
                                         F.Fmsy=Store.sens[[s]]$output$F.Fmsy[1:length(Years)])
        }

      }
      #dummy.store[[i]]=Store.sens
      Out.Scens=Out.Scens%>%
                dplyr::select(Species,Scenario,AgeMat,M.dist,Mmean,M.CV,Klow,Kup,
                              fmsym.dist,fmsy.m,fmsym.CV,bmsyk.dist,bmsyk.mean,bmsyk.CV,
                              b1k.dist,b1k.low,b1k.up,btk.dist,btk.low,btk.up)%>%
                mutate(Mmean=round(Mmean,2),
                       M.CV=round(M.CV,2),
                       Klow=round(Klow),
                       Kup=round(Kup),
                       fmsym.dist=ifelse(fmsym.dist=="lnorm","Lognormal",fmsym.dist),
                       M.dist=ifelse(M.dist=="lnorm","Lognormal",M.dist),
                       bmsyk.dist=capitalize(bmsyk.dist),
                       b1k.dist=ifelse(b1k.dist=="unif","Uniform",b1k.dist),
                       btk.dist=ifelse(btk.dist=="unif","Uniform",btk.dist))
      
      dummy.store.sens.table[[i]]=Out.Scens
      dummy.store.accept.rate[[i]]=do.call(rbind,Out.accept.rate)
      dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
      dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
      dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
      dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
      dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
      dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
    }
    
    #Catch_only[[w]]=dummy.store    #too big an object, cannot store whole model runs due to memmory size limitations
    Catch_only[[w]]$sens.table=dummy.store.sens.table
    Catch_only[[w]]$accept.rate=dummy.store.accept.rate
    Catch_only[[w]]$estimates=dummy.store.estimates
    Catch_only[[w]]$rel.biom=dummy.store.rel.biom
    Catch_only[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
    Catch_only[[w]]$f.series=dummy.store.f.series
    Catch_only[[w]]$B.Bmsy=dummy.store.B.Bmsy
    Catch_only[[w]]$F.Fmsy=dummy.store.F.Fmsy
    Catch_only[[w]]$ensemble=dummy.store.ensemble
    
    rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,
       dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,
       dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.accept.rate,
       dummy.store.ensemble)
  }
  
  #11.2. CMSY    
  #summary of method: http://toolbox.frdc.com.au/wp-content/uploads/sites/19/2021/04/CMSY.html
  #note: Milk shark need 1e5 simulations to sample enough accepted r-k combos
  if(names(Catch_only)[w]=="CMSY")
  {
    dummy.store=vector('list',N.sp)     #takes 0.008 secs per iteration per species per scenario
    names(dummy.store)=Keep.species 
    dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
      dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.accept.rate=
      dummy.store.ensemble=dummy.store
    for(i in 1:length(dummy.store))  
    {
      this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(List.sp[[i]]$Name),"/",AssessYr,"/CMSY",sep='')
      if(!dir.exists(this.wd))dir.create(this.wd)
      
      ktch=ktch.combined%>%
        filter(Name==names(dummy.store)[i])
      year=ktch$finyear
      catch=ktch$Tonnes
      
      Scens=List.sp[[i]]$Sens.test$CMSY%>%
                mutate(Species=capitalize(names(dummy.store)[i]))
      Store.sens=vector('list',nrow(Scens))
      names(Store.sens)=Scens$Scenario
      this.wd1=this.wd
      
      Out.Scens=Scens
      Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
        Out.B.Bmsy=Out.F.Fmsy=Out.accept.rate=vector('list',length(Store.sens))
      
      
      for(s in 1:length(Store.sens))
      {
        print(paste("___________","CMSY Scenario",Scens$Scenario[s],"___________",names(dummy.store)[i]))
        this.wd=paste(this.wd1,names(Store.sens)[s],sep='/')
        if(!dir.exists(this.wd))dir.create(this.wd)
        
        #Priors
        RES=RESILIENCE[[i]]
        if(Scens$r.prob.min[s]==0)
        {
          r.range=NA
          k.range=NA
          Bf.low=NA
          Bf.hi=NA
        }else
        {
          Mn=min(store.species.r[[i]]$mean,Max.r.value)  #some life history pars yield unrealistically high r
          r.range=quantile(rnorm(1e3,
                                 mean=Mn,
                                 sd=store.species.r[[i]]$sd),
                           probs=c(Scens$r.prob.min[s],Scens$r.prob.max[s]))
          
          k.range=c(Scens$Klow[s],Scens$Kup[s])
          
          Bf.low=List.sp[[i]]$FINALBIO[1]
          Bf.hi=List.sp[[i]]$FINALBIO[2]
           
        }
        Proc.error=Scens$Proc.error[s]
        
        Bo.low=List.sp[[i]]$STARTBIO[1]
        Bo.hi=List.sp[[i]]$STARTBIO[2]
        

        #Run model
        Store.sens[[s]]=apply.CMSY(year=year,
                                   catch=catch,
                                   r.range=r.range,
                                   k.range=k.range,
                                   Bo.low=Bo.low,
                                   Bo.hi=Bo.hi,
                                   Bf.low=Bf.low,
                                   Bf.hi=Bf.hi,
                                   outfile=paste(this.wd,'Appendix_fit',sep='/'),
                                   nsims=Scens$Sims[s],
                                   Proc.error=Proc.error,
                                   RES=RES)
        
        #Store acceptance rate
        Out.accept.rate[[s]]=data.frame(
          Acceptance=Store.sens[[s]]$output$acceptance.rate,
          Scenario=Scens$Scenario[s])

        #Store scenarios
        Out.Scens$Bo.low=Bo.low     
        Out.Scens$Bo.hi=Bo.hi
        Out.Scens$Bf.low=Bf.low
        Out.Scens$Bf.hi=Bf.hi
        Out.Scens$r.low=r.range[1]
        Out.Scens$r.up=r.range[2]
        Out.Scens$Klow[s]=k.range[1]
        Out.Scens$Kup[s]=k.range[2]
        
        #Store estimates
        d1=Store.sens[[s]]$output$Statistics$output  
        d1=d1[,grep(paste(c('50%','2.5%','97.5%'),collapse='|'),colnames(d1))]%>%
          data.frame
        d1=d1[,-grep('Perc',colnames(d1))]
        colnames(d1)=c("Lower.95","Median","Upper.95")
        d1=d1%>%
          mutate(Parameter=rownames(d1),
                 Model='CMSY',
                 Scenario=names(Store.sens)[s])%>%
          relocate(Model,Scenario,Parameter,Lower.95)
        Out.estimates[[s]]=d1
        
        #Store trajectories
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                          mods=names(Catch_only)[w],
                                          Type='Depletion',
                                          scen=Scens$Scenario[s],
                                          Katch=catch)
        Out.rel.biom[[s]]=dummy$Dat
        Out.probs.rel.biom[[s]]=dummy$Probs
        
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                          mods=names(Catch_only)[w],
                                          Type='F.series',
                                          scen=Scens$Scenario[s],
                                          Katch=catch)
        Out.f.series[[s]]=dummy$Dat
        
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                          mods=names(Catch_only)[w],
                                          Type='B.Bmsy',
                                          scen=Scens$Scenario[s],
                                          Katch=catch)
        Out.B.Bmsy[[s]]=dummy$Dat
        
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                          mods=names(Catch_only)[w],
                                          Type='F.Fmsy',
                                          scen=Scens$Scenario[s],
                                          Katch=catch)
        Out.F.Fmsy[[s]]=dummy$Dat
        rm(dummy)
        
        #Store quantities for ensemble model
        if(Scens$Scenario[s]=='S1')
        {
          Years=Store.sens[[s]]$output$Years
          F.Fmsy=Store.sens[[s]]$output$F.Fmsy[,1:length(Years)]%>%
            data.frame
          colnames(F.Fmsy)=Years
          dummy.store.ensemble[[i]]=list(Depletion=Store.sens[[s]]$output$Depletion.traj[1:length(Years)],
                                         B.Bmsy=Store.sens[[s]]$output$B.Bmsy[1:length(Years)],
                                         F.Fmsy=F.Fmsy) 
          
        }
        
        
        #Display Priors vs Posteriors for base case scenario (S1)
        if(Scens$Scenario[s]=='S1')
        {
          par.list=c("r.range","k.range","stb.low","stb.hi","endb.low","endb.hi")
          dummy.prior=Store.sens[[s]]$input[par.list]
          names(dummy.prior)=tolower(names(dummy.prior))
          dummy.prior$stb=c(dummy.prior$stb.low,dummy.prior$stb.hi)
          dummy.prior$endb=c(dummy.prior$endb.low,dummy.prior$endb.hi)
          dummy.prior=within(dummy.prior, rm(stb.low,stb.hi,endb.low,endb.hi))
          for(d in 1:length(dummy.prior))
          {
            dummy.prior[[d]]=list(dist='unif',
                                  low=dummy.prior[[d]][1],
                                  up=dummy.prior[[d]][2])
          }
          dummy.post=as.data.frame(Store.sens[[s]]$output$Statistics$deplet)
          dummy.post=dummy.post[,c(match(c('r','K'),names(dummy.post)),1,ncol(dummy.post)-3)]
          colnames(dummy.post)=c('r','K','Initial depletion','Final depletion')
          names(dummy.prior)=names(dummy.post)
          
          pars=colnames(dummy.post)
          out=vector('list',length(pars))
          for(p in 1:length(pars))
          {
            out[[p]]=rbind(data.frame(Distribuion="Prior",
                                      Value=fn.prior(d=dummy.prior[[pars[p]]])),
                           data.frame(Distribuion="Posterior",
                                      Value=dummy.post[[pars[p]]]))%>%
              mutate(Parameter=pars[p])
          }
          rm(dummy.prior,dummy.post)
          
          fn.show.density(d=do.call(rbind,out),NCOL=2)
          ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                       capitalize(names(dummy.store)[i]),"/",AssessYr,"/",'CMSY',"/Prior.and.posterior.tiff",sep=''),
                 width = 12,height = 14, dpi = 300, compression = "lzw")
        }
      }
      
      #dummy.store[[i]]=Store.sens
      
      Out.Scens=Out.Scens%>%
        dplyr::select(Species,Scenario,r.low,r.up,Klow,Kup,Bo.low,Bo.hi,Bf.low,Bf.hi,Proc.error)%>%
        mutate(Klow=round(Klow),
               Kup=round(Kup),
               r.low=round(r.low,3),
               r.up=round(r.up,3))
      
      dummy.store.sens.table[[i]]=Out.Scens
      dummy.store.accept.rate[[i]]=do.call(rbind,Out.accept.rate)
      dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
      dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
      dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
      dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
      dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
      dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
    }
    
    #Catch_only[[w]]=dummy.store   #too big an object, cannot store whole model runs due to memmory size limitations
    Catch_only[[w]]$sens.table=dummy.store.sens.table
    Catch_only[[w]]$accept.rate=dummy.store.accept.rate
    Catch_only[[w]]$estimates=dummy.store.estimates
    Catch_only[[w]]$rel.biom=dummy.store.rel.biom
    Catch_only[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
    Catch_only[[w]]$f.series=dummy.store.f.series
    Catch_only[[w]]$B.Bmsy=dummy.store.B.Bmsy
    Catch_only[[w]]$F.Fmsy=dummy.store.F.Fmsy
    Catch_only[[w]]$ensemble=dummy.store.ensemble
    
    rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,
       dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,
       dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.accept.rate,
       dummy.store.ensemble)
  }
  
  #11.3. OCOM assessment (Zhou et al (2018))
  #note: not used as it doesn't allow for lightly depleted stocks and catch reductions
  # due to effort reduction rather than abundance
  if(names(Catch_only)[w]=="OCOM")
  {
    dummy.store=vector('list',N.sp)  #takes 20 secs per species (1e4 iterations)    
    names(dummy.store)=Keep.species             
    system.time({for(i in 1:length(dummy.store))  
    {
      print(paste("OCOM ","--",names(dummy.store)[i]))
      this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(List.sp[[i]]$Name),"/",AssessYr,"/OCOM",sep='')
      if(!dir.exists(this.wd))dir.create(this.wd)
      ktch=ktch.combined%>%
        filter(Name==names(dummy.store)[i])
      dummy.store[[i]] <-apply.OCOM(year=ktch$finyear,
                                    catch=ktch$Tonnes,
                                    M=mean(apply(store.species.M[[i]],2,mean,na.rm=T)),
                                    outfile=paste(this.wd,'Appendix_fit',sep='/'))
    }})
    Catch_only[[w]]=dummy.store
    rm(dummy.store)
  }
  
  #11.4. JABBA - catch only (Winker et al 2018)   
  #summary of method: https://github.com/jabbamodel/JABBA
  if(names(Catch_only)[w]=="JABBA")
  {
    dummy.store=vector('list',N.sp)     #takes 0.002 secs per iteration per species per scenario
    names(dummy.store)=Keep.species
    dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
      dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.Kobe.probs=
      dummy.store.ensemble=dummy.store
    for(i in 1:length(dummy.store))  
    {
      this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(List.sp[[i]]$Name),"/",AssessYr,"/JABBA",sep='')
      if(!dir.exists(this.wd))dir.create(this.wd)
      
      #Catch
      ktch=ktch.combined%>%
        filter(Name==names(dummy.store)[i])%>%
        rename(Year=finyear,
               Total=Tonnes)%>%
        ungroup()%>%
        dplyr::select(Year,Total)%>%
        arrange(Year)%>%
        data.frame
      if(names(dummy.store)[i]=='milk shark') ktch$Total[1:10]=1 #convergence issues with first 10 years for milk shark, but catch is 0 anyways
      Scens=List.sp[[i]]$Sens.test$JABBA%>%
                mutate(Species=capitalize(names(dummy.store)[i]))
      Store.sens=vector('list',nrow(Scens))
      names(Store.sens)=Scens$Scenario
      this.wd1=this.wd
      
      Out.Scens=Scens
      Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
        Out.B.Bmsy=Out.F.Fmsy=vector('list',length(Store.sens))
      
      for(s in 1:length(Store.sens))
      {
        print(paste("___________","JABBA Scenario",Scens$Scenario[s],"___________",names(dummy.store)[i]))
        this.wd=paste(this.wd1,names(Store.sens)[s],sep='/')
        if(!dir.exists(this.wd))dir.create(this.wd)
        
        #Priors 
        bmsyk.mean=Scens$bmsyk[s]
        Proc.error=Scens$Proc.error[s]
        r.CV.multi=Scens$r.CV.multiplier[s]
        K.prior=c(Scens$K.mean[s],Scens$K.CV[s])
        
        Mn=min(store.species.r[[i]]$mean,Max.r.value)  #some life history pars yield unrealistically high r
        r.prior=c(Mn, (store.species.r[[i]]$sd/Mn)*Scens$r.CV.multiplier[s])
        Ktch.CV=Scens$Ktch.CV[s]
        Bint=runif(1000,List.sp[[i]]$STARTBIO[1],List.sp[[i]]$STARTBIO[2])
        Bint.mean=mean(Bint)
        Bint.CV=sd(Bint)/Bint.mean
        psi.prior=c(Bint.mean,Bint.CV)
        #psi.prior=c(1,0.1) #Winker et al 2018 School shark
        
        Bfin=runif(1000,List.sp[[i]]$FINALBIO[1],List.sp[[i]]$FINALBIO[2])
        Bfin.mean=mean(Bfin)          
        Bfin.CV=sd(Bfin)/Bfin.mean
        b.prior=c(Bfin.mean,Bfin.CV,max(ktch$Year),"bk")
        
        Rdist = "lnorm"
        Kdist="lnorm"  
        PsiDist='lnorm'
        
        #Put inputs together
        input=list(Ktch=ktch,
                   MDL="Schaefer",
                   Ktch.CV=Ktch.CV,
                   ASS=names(dummy.store)[i],
                   Rdist = Rdist,
                   Rprior = r.prior,
                   r.CV.multiplier=r.CV.multi,
                   Kdist= Kdist,
                   Kprior=K.prior,    
                   PsiDist= PsiDist,
                   Psiprior=psi.prior,   
                   Bprior=b.prior,    
                   BMSYK=bmsyk.mean)  
        
        #run model
        THIN=5
        CHAINS=2
        BURNIN=min(0.15*Scens$Sims[s],5000)
        output=apply.JABBA(Ktch=input$Ktch,
                           CPUE=NULL,
                           CPUE.SE=NULL,
                           MDL=input$MDL,
                           Ktch.CV=input$Ktch.CV,
                           ASS=input$ASS,
                           Rdist = input$Rdist,
                           Rprior = input$Rprior,
                           Kdist=input$Kdist,
                           Kprior=input$Kprior,
                           PsiDist=input$PsiDist,
                           Psiprior=input$Psiprior,
                           Bprior=input$Bprior,
                           BMSYK=input$BMSYK,
                           output.dir=this.wd,
                           outfile="Appendix_fit",
                           Sims=Scens$Sims[s],
                           Proc.error.JABBA=Proc.error,
                           thinning = THIN,
                           nchains = CHAINS,
                           burn.in= BURNIN)
        Store.sens[[s]]=list(input=input,output=output)
        
        #Store Scenarios
        Out.Scens$Bo.mean=Bint.mean
        Out.Scens$Bo.CV=Bint.CV
        Out.Scens$Bf.mean=Bfin.mean
        Out.Scens$Bf.CV=Bfin.CV
        Out.Scens$r.mean=store.species.r[[i]]$mean
        Out.Scens$r.cv=store.species.r[[i]]$sd/store.species.r[[i]]$mean
        Out.Scens$Rdist=Rdist
        Out.Scens$Kdist=Kdist
        Out.Scens$PsiDist=PsiDist
        Out.Scens$bmsyk.mean=Out.Scens$bmsyk
        
        #Store estimates  
        d1=Store.sens[[s]]$output$estimates
        d1=d1[,grep(paste(c('mu','lci','uci'),collapse='|'),colnames(d1))]%>%
          data.frame
        colnames(d1)=c("Median","Lower.95","Upper.95")
        d1=d1%>%
          mutate(Parameter=rownames(d1),
                 Model='JABBA',
                 Scenario=names(Store.sens)[s])%>%
          relocate(Model,Scenario,Parameter,Lower.95)%>%
          filter(Parameter%in%c("K","r","psi","Hmsy","SBmsy","MSY"))   
        Out.estimates[[s]]=d1
        
        #Store trajectories
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                          mods=names(Catch_only)[w],
                                          Type='Depletion',
                                          scen=Scens$Scenario[s],
                                          Katch=ktch$Total)
        Out.rel.biom[[s]]=dummy$Dat
        Out.probs.rel.biom[[s]]=dummy$Probs
        
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                   mods=names(Catch_only)[w],
                                                   Type='F.series',
                                                   scen=Scens$Scenario[s],
                                                   Katch=ktch$Total)
        Out.f.series[[s]]=dummy$Dat
        
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                    mods=names(Catch_only)[w],
                                                    Type='B.Bmsy',
                                                    scen=Scens$Scenario[s],
                                                    Katch=ktch$Total)
        Out.B.Bmsy[[s]]=dummy$Dat
        
        dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                    mods=names(Catch_only)[w],
                                                    Type='F.Fmsy',
                                                    scen=Scens$Scenario[s],
                                                    Katch=ktch$Total)
        Out.F.Fmsy[[s]]=dummy$Dat
        rm(dummy)
        
        if(Scens$Scenario[s]=='S1') Out.Kobe.probs=Store.sens[[s]]$output$kobe%>%dplyr::select(stock,harvest)
        
        #Evaluate posterior  
        #https://cran.r-project.org/web/packages/ggmcmc/vignettes/using_ggmcmc.html
        #MCMC trace-plots, Geweke (2003) and Gelman and Rubin (1992) diagnostics 
        if(Scens$Scenario[s]=='S1')
        {
          these.post.pars=c('K','r','psi')  #estimable pars
          NN=nrow(output$pars_posterior)/CHAINS
          Post=Store.sens[[s]]$output$pars_posterior%>%
            dplyr::select(all_of(these.post.pars))%>%
            mutate(Iteration=rep(1:NN,CHAINS),
                   Chain=rep(1:CHAINS,each=NN))%>%
            gather(Parameter,value,-Iteration,- Chain)
          tibble()
          attributes(Post)$nChains=CHAINS
          attributes(Post)$nParameters=length(these.post.pars)
          attributes(Post)$nIterations=NN
          attributes(Post)$nBurnin=BURNIN
          attributes(Post)$nThin=THIN
          
          Post.density=ggs_density(Post)+ggtitle("Density")
          Post.trace.plot=ggs_traceplot(Post)+ggtitle("Trace plot")+theme(axis.text.x = element_text(angle = 90))
          Post.running.mean=ggs_running(Post)+ggtitle('Running mean')+theme(axis.text.x = element_text(angle = 90))
          Post.autocorr=ggs_autocorrelation(Post)+ggtitle('Autocorrelation')
          Gelman.diag=ggs_Rhat(Post) + 
                        geom_point(size=3)+
                        xlab("R_hat")+
                        ggtitle('Gelman diagnostic')
          Geweke.diag= ggs_geweke(Post)+
                        geom_point(size=3)+
                        ggtitle('Geweke diagnostic')+
                        theme(legend.position = 'none')+
                        scale_color_manual(values=c("#F8766D", "#00BFC4", "#56B4E9"))
          
          pp=ggarrange(plotlist = list(Gelman.diag,Geweke.diag),ncol=1)
          ggarrange(plotlist = list(Post.density,Post.trace.plot,Post.running.mean,
                                    Post.autocorr,pp),
                    ncol=5,common.legend = TRUE)
          
          ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                       capitalize(names(dummy.store)[i]),"/",AssessYr,"/",'JABBA',
                       "/Convergence diagnostics.tiff",sep=''),
                 width = 16,height = 8, dpi = 300, compression = "lzw")
          
          
        }
        
         #Display Priors vs Posteriors for base case scenario (S1)
        if(Scens$Scenario[s]=='S1')
        {
          par.list=c("r","K","psi")
          stuff=Store.sens[[s]]$input
          dummy.prior=as.list(par.list)
          names(dummy.prior)=par.list
          dummy.prior$r=list(dist=stuff$Rdist,
                             mean=log(stuff$Rprior[1]),
                             sd=stuff$Rprior[2])
          dummy.prior$K=list(dist='lnorm',
                             mean=log(Store.sens[[s]]$output$settings$K.pr[1]),  
                             sd=Store.sens[[s]]$output$settings$K.pr[2])
          dummy.prior$psi=list(dist=stuff$PsiDist,
                               mean=log(stuff$Psiprior[1]),
                               sd=stuff$Psiprior[2])
          dummy.post=Store.sens[[s]]$output$pars_posterior%>%
            dplyr::select(r,K,psi)
          
          names(dummy.prior)=names(dummy.post)
          
          pars=colnames(dummy.post)
          out=vector('list',length(pars))
          
          
          for(p in 1:length(pars))
          {
            Value.prior=fn.prior(d=dummy.prior[[pars[p]]])
            if(pars[p]=="K") Value.prior=fn.prior(d=dummy.prior[[pars[p]]],MAX=k.fun.up(stuff$Ktch$Total))
            out[[p]]=rbind(data.frame(Distribuion="Prior",
                                      Value=Value.prior),
                           data.frame(Distribuion="Posterior",
                                      Value=dummy.post[[pars[p]]]))%>%
              mutate(Parameter=pars[p])
          }
          rm(stuff,dummy.prior,dummy.post)
          
          fn.show.density(d=do.call(rbind,out),NCOL=1)
          ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                       capitalize(names(dummy.store)[i]),"/",AssessYr,"/",'JABBA',"/Prior.and.posterior.tiff",sep=''),
                 width = 12,height = 14, dpi = 300, compression = "lzw")
        }
        
        #Store quantities for ensemble model
        if(Scens$Scenario[s]=='S1')
        {
          dummy.store.ensemble[[i]]=list(Depletion=Store.sens[[s]]$output$posteriors$P,
                                         B.Bmsy=Store.sens[[s]]$output$posteriors$BtoBmsy,
                                         F.Fmsy=Store.sens[[s]]$output$posteriors$HtoHmsy) 
        }
        
      }
      #dummy.store[[i]]=Store.sens
      
      Out.Scens=Out.Scens%>%
        dplyr::select(Species,Scenario,
                      Rdist,r.mean,r.cv,
                      Kdist,K.mean,K.CV,
                      PsiDist,Bo.mean,Bo.CV,
                      Bf.mean,Bf.CV,
                      bmsyk.mean,Proc.error,Ktch.CV)%>%
        mutate(K.mean=round(K.mean),
               r.mean=round(r.mean,3),
               bmsyk.mean=round(bmsyk.mean,3),
               r.cv=round(r.cv,3),
               Bo.mean=round(Bo.mean,3),
               Bo.CV=round(Bo.CV,3),
               Bf.mean=round(Bf.mean,3),
               Bf.CV=round(Bf.CV,3),
               Rdist=ifelse(Rdist=='lnorm','Lognormal',Rdist),
               Kdist=ifelse(Kdist=='lnorm','Lognormal',Kdist),
               PsiDist=ifelse(PsiDist=='lnorm','Lognormal',PsiDist))
      
      dummy.store.sens.table[[i]]=Out.Scens
      dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
      dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
      dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
      dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
      dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
      dummy.store.Kobe.probs[[i]]=Out.Kobe.probs  
      dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
    }
    
    #Catch_only[[w]]=dummy.store       #too big an object, cannot store whole model runs due to memmory size limitations
    Catch_only[[w]]$sens.table=dummy.store.sens.table
    Catch_only[[w]]$estimates=dummy.store.estimates
    Catch_only[[w]]$rel.biom=dummy.store.rel.biom
    Catch_only[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
    Catch_only[[w]]$f.series=dummy.store.f.series
    Catch_only[[w]]$B.Bmsy=dummy.store.B.Bmsy
    Catch_only[[w]]$F.Fmsy=dummy.store.F.Fmsy
    Catch_only[[w]]$Kobe.probs=dummy.store.Kobe.probs
    Catch_only[[w]]$ensemble=dummy.store.ensemble
    
    rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,
       dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,
       dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.ensemble)
  }
}

  #19.2. COM weighted average  
    #19.2.1. get COM weights
#r groups
All.rs=do.call(rbind,store.species.r)
r.list=All.rs$mean
names(r.list)=1:N.sp
r.groups=list(low=c(min(r.list),quantile(r.list,0.2499)),
              medium=c(quantile(r.list,0.25),quantile(r.list,0.7499)),
              high=c(quantile(r.list,0.75),0.65))
if(do.ensemble.simulations)   #579 sec per iteration    
{
  #Ks
  K.list=c(1e4)
  
  #Init depletion
  Init.depl=c(1,0.7)
  
  
  #Create exploitation histories
  #note: setting U at quantile(r.list,0.8) creates a range of exploitation histories ranging
  #       between 10 times Fmsy and Fmsy depending on life history
  YRS=1975:2020
  HarvestRate1=seq(0,quantile(r.list,0.8),length.out = (length(YRS)/2))+rnorm(length(YRS)/2,0,.01)
  HarvestRate2=seq(0,quantile(r.list,0.5),length.out = (length(YRS)/2))+rnorm(length(YRS)/2,0,.01)
  
  Inputd=data.frame(YRS=YRS)%>%
    mutate(HarvestRate1=c(HarvestRate1,rev(HarvestRate1)),
           HarvestRate2=c(HarvestRate2,rev(HarvestRate2)),
           #HarvestRate1=rep(quantile(r.list,0.5),length(YRS))+rnorm(length(YRS),0,.01),  #constant U is not representative
           #HarvestRate2=rep(quantile(r.list,0.1),length(YRS))+rnorm(length(YRS),0,.01),
           HarvestRate3=HarvestRate1, 
           HarvestRate4=HarvestRate2)
  
  n=round(length(YRS)/6)
  Inputd$HarvestRate1[(nrow(Inputd)-n):nrow(Inputd)]=Inputd$HarvestRate1[(nrow(Inputd)-n)]
  Inputd$HarvestRate2[(nrow(Inputd)-n):nrow(Inputd)]=Inputd$HarvestRate2[(nrow(Inputd)-n)]
  
  n=1:ceiling(length(YRS)*.35)
  Inputd$HarvestRate3[n]=0
  Inputd$HarvestRate4[n]=0
  n=(nrow(Inputd)-ceiling(length(YRS)*.35)):nrow(Inputd)
  Inputd$HarvestRate3[n]=0
  Inputd$HarvestRate4[n]=0
  Inputd$HarvestRate3=ifelse(Inputd$HarvestRate3<0,0,Inputd$HarvestRate3)
  Inputd$HarvestRate4=ifelse(Inputd$HarvestRate4<0,0,Inputd$HarvestRate4)
  
  
  #Factorial design
  Factorial=expand.grid(r=names(r.groups),K=K.list,p=Init.depl,ktch=paste('HarvestRate',1:(ncol(Inputd)-1),sep=''))
  Factorial$ktch=as.character(Factorial$ktch)
  

  #Operation model
  SPM=function(K,r,Init.dep,HRscen)
  {
    HarvestRate=Inputd[,HRscen]
    HarvestRate=HarvestRate+rnorm(length(YRS),0,.01)
    HarvestRate=ifelse(HarvestRate<0,0,HarvestRate)
    
    Ktch=rep(NA,length(YRS))
    Bt=rep(NA,length(Ktch))
    Bt[1]=K*Init.dep
    Ktch[1]=HarvestRate[1]*Bt[1]
    for(t in 2:length(Bt))
    {
      Bt[t]=Bt[t-1]+r*Bt[t-1]*(1-Bt[t-1]/K)-Ktch[t-1]
      Ktch[t]=HarvestRate[t]*Bt[t]
    }
    
    U=Ktch/Bt
    Ft=-log(1-U)
    Bmsy=K/2
    MSY=K*r/4
    Fmsy=r/2
    Depletion=Bt/K
    B.Bmsy=Bt/Bmsy
    F.Fmsy=Ft/Fmsy
    #ii=(length(YRS)-4):length(YRS)
    
    return(list(Ktch=Ktch, r=r, K=K, Init.dep=Init.dep,HRscen=HRscen,
                Depletion=Depletion, B.Bmsy=B.Bmsy, F.Fmsy=F.Fmsy))
  }
  
  #Test if OM can retrieve assumed pars
  this.wd='C:/Superensemble'
  setwd(this.wd)
  check.can.rekv.pars=FALSE
  if(check.can.rekv.pars)
  {
    Q=1e-3
    
    SPM.fit=function(PARS)
    {
      K=exp(PARS[1])
      r=exp(PARS[2])
      Init.dep=exp(PARS[3])
      q=exp(PARS[4])
      
      HarvestRate=Inputd[,HRscen]
      HarvestRate=HarvestRate+rnorm(length(YRS),0,.01)
      HarvestRate=ifelse(HarvestRate<0,0,HarvestRate)
      
      Ktch=rep(NA,length(YRS))
      Bt=rep(NA,length(Ktch))
      Bt[1]=K*Init.dep
      Ktch[1]=HarvestRate[1]*Bt[1]
      for(t in 2:length(Bt))
      {
        Bt[t]=Bt[t-1]+r*Bt[t-1]*(1-Bt[t-1]/K)-Ktch[t-1]
        Ktch[t]=HarvestRate[t]*Bt[t]
      }
      
      Obs=length(Bt)
      Est_CPUE = q * Bt
      sqRes = (ln_CPUE - log(Est_CPUE)) * (ln_CPUE - log(Est_CPUE))
      sumsq = sum(sqRes)
      var = sumsq / Obs
      stdev = sqrt(var)
      # calculation of negative log-likelihood value - see Haddon (2001)
      NLL = (Obs / 2) * (log(2 * pi) + (2 * log(stdev)) + 1)
      
      if(what=='fit') out=NLL
      if(what=='simulate') out= list(cpue=Est_CPUE,ktch=Ktch,Bt=Bt)
      return(out)
    }
    
    store.estimpars=store.inpupars=store.ktch=store.bt=vector('list',length=nrow(Factorial))
    Thisvec=c(1,10,22,31,43,73,85,94,106,115,127,136,148,157)
    for(s in Thisvec)
    {
      HRscen=Factorial$ktch[s]
      PARS=log(c(K=Factorial$K[s],r=Factorial$r[s],Init.dep=Factorial$p[s],q=Q))
      
      what='simulate'
      OBJ=SPM.fit(PARS)
      store.ktch[[s]]=OBJ$ktch
      store.bt[[s]]=OBJ$Bt
      Obs_CPUE=OBJ$cpue
      ln_CPUE <- log(Obs_CPUE+rnorm(length(Obs_CPUE),0,1)) 
      
      what='fit'
      init.par=PARS+rnorm(length(PARS),0,.5)
      nlmb <- nlminb(init.par, SPM.fit, gradient = NULL,hessian = TRUE)
      store.estimpars[[s]]=nlmb$par
      store.inpupars[[s]]=PARS
      
      #opt <- optim(par=init.par,fn=SPM.fit,method="Nelder-Mead",hessian=TRUE)
      #solve(opt$hessian) 
      
    }
    
    fn.fig("Recover.pars_biom&ktch.fit",2400,2400)
    par(mfcol=c(5,3),mai=c(.3,.5,.2,.5))
    for(s in Thisvec)
    {
      plot(store.ktch[[s]],cex=2,type='l',ylab='')
      par(new=T)
      plot(store.bt[[s]],col=2,type='l',axes = FALSE, xlab = '', ylab ='')
      axis(side=4, at = pretty(range(store.bt[[s]])))
      mtext(paste('K=',Factorial$K[s],'r=',round(Factorial$r[s],2),
                  'Init.dep=',Factorial$p[s],Factorial$ktch[s]),3,cex=.55)
      
    }
    mtext("Catch", side=2,-1.5,outer=T,cex=1.5)
    mtext("Biomass", side=4, line=-1,col=2,outer=T,cex=1.5)
    dev.off()
    
    fn.fig("Recover.pars_pars.fit",2400,2400)
    par(mfcol=c(5,3),mar=rep(1,4),mai=c(.3,.5,.2,.5))
    for(s in Thisvec)
    {
      plot(store.estimpars[[s]],cex=2,pch=19,ylab='Ln.par')
      points(store.inpupars[[s]],col=2,cex=2)
      mtext(paste('K=',Factorial$K[s],'r=',round(Factorial$r[s],2),
                  'Init.dep=',Factorial$p[s],Factorial$ktch[s]),3,cex=.6)
      
    }
    legend('topright',c('Estimate','Input'),pch=19,col=1:2,bty='n')
    dev.off()
  }
  
  #1. Simulate biomass trajectories
  n.reps=20
  Store.sims=vector('list',nrow(Factorial)*n.reps)
  nnn=seq(0,(length(Store.sims)-n.reps),n.reps)
  for(s in 1:nrow(Factorial))
  {
    for(i in 1:n.reps)
    {
      ar=r.groups[[match(Factorial$r[s],names(r.groups))]]
      ar=runif(1,ar[1],ar[2])
      ar=log(ar)
      x <- exp(rmvnorm(n=1, mean=c(ar,log(Factorial$K[s])), 
                       sigma=matrix(c(0.0002454185,-0.0001840755,-0.0001840755,0.0001383331),
                                    nrow=2))) #variance cov variance from fitting SPM to real data
      
      dumi=SPM(K=x[,2],r=x[,1],Init.dep=Factorial$p[s],
                    HRscen=Factorial$ktch[s])
      s.i=nnn[s]+i
      Store.sims[[s.i]]=dumi
    }
  }
  
  Inputs=Store.sims
  for(s in 1:length(Inputs))
  {
    Inputs[[s]]=data.frame(Ktch=Store.sims[[s]]$Ktch)%>%
      mutate(r=Store.sims[[s]]$r,
             K=Store.sims[[s]]$K, 
             Init.dep=Store.sims[[s]]$Init.dep,
             HRscen=Store.sims[[s]]$HRscen,
             Depletion=Store.sims[[s]]$Depletion, 
             B.Bmsy=Store.sims[[s]]$B.Bmsy,
             F.Fmsy=Store.sims[[s]]$F.Fmsy,
             Year=YRS)
  }
  
  
  # do.call(rbind,Inputs)%>%
  #   mutate(dummy1=paste(HRscen,r),
  #          dummy=paste(paste('p=',Init.dep,sep=''),paste('K=',K,sep=''),sep=', '))%>%
  #   ggplot(aes(Year,Ktch,color=dummy1))+
  #   geom_point()+
  #   geom_line(alpha=0.3)+
  #   facet_wrap(~dummy,scales='free',nrow=2)+
  #   theme(legend.position = 'none')
  
  
  #2. Fit COMs to catch trajectories as per stock assessments
  get.rs=data.frame(index=names(r.list),r=r.list)  
  system.time({for(s in 1:length(Inputs))   
  {
    Ktch=Inputs[[s]]$Ktch
    Klow=k.fun.low(Ktch) 
    Kup=k.fun.up(Ktch)
    
    Yrs=Inputs[[s]]$Year
    dis.r.range=do.call(rbind,r.groups)%>%   
      data.frame%>%
      mutate(Get=between(unique(Inputs[[s]]$r),.[[1]] , .[[2]]))%>%
      filter(Get==TRUE)
    Inputs[[s]]$r.group=row.names(dis.r.range)
    i=get.rs%>%
      filter(r>=dis.r.range$V1 & r<=dis.r.range$X25.)
    i=as.numeric(i[sample(i$index,1),'index'])
    
    Btklow=0.01
    Btkup=List.sp[[i]]$FINALBIO[2]
    b1k.low=List.sp[[i]]$STARTBIO[1]
    b1k.up=List.sp[[i]]$STARTBIO[2]
    
    
    #DBSRA
    print(paste("___________","DBSRA ","______ s=",s))
    Scens=List.sp[[i]]$Sens.test$DBSRA%>%filter(Scenario=='S1')
    
    #priors
    AgeMat=Scens$AgeMat
    Mmean=Scens$Mmean  
    Msd=Scens$Msd
    Mcv=Msd/Mmean
    M.dist="lnorm"
    fmsy.m=Scens$fmsy.m
    fmsym.dist="lnorm"
    fmsym.CV=0.1
    bmsyk.mean=Scens$bmsyk
    bmsyk.dist="beta"
    bmsyk.CV=0.1
    btk.dist="unif"
    b1k.dist="unif"
    
    #run model
    out.dummy=apply.DBSRA(year=Yrs,
                          catch=Ktch,
                          catchCV=NULL,  
                          catargs=list(dist="none",low=0,up=Inf,unit="MT"),  #catch CV not available
                          agemat=AgeMat,
                          k=list(low=Klow,up=Kup,tol=0.01,permax=1000),
                          b1k=list(dist=b1k.dist,low=b1k.low,up=b1k.up,mean=1,sd=0.1),  #mean and sd not used if 'unif'
                          btk=list(dist="unif",low=Btklow,up=Btkup,mean=1,sd=0.1,refyr=max(Yrs)),  #reference year
                          fmsym=list(dist="lnorm",low=0.1,up=2,mean=log(fmsy.m),sd=fmsym.CV), # Cortes & Brooks 2018. Low and up not used if 'lnorm'  
                          bmsyk=list(dist="beta",low=0.05,up=0.95,mean=bmsyk.mean,sd=bmsyk.CV),  
                          M=list(dist="lnorm",low=0.001,up=1,mean=log(Mmean),sd=Mcv),
                          graph=1,
                          nsims=Scens$Sims, 
                          grout=1,
                          WD=this.wd,
                          outfile="Appendix_fit")
    
    Dep=apply(out.dummy$output$Depletion.traj,2,median)
    #B.Bmsy=apply(out.dummy$output$B.Bmsy,2,median)
    #F.Fmsy=apply(out.dummy$output$F.Fmsy,2,median)
    
    Inputs[[s]]$Depletion_DBSRA=Dep[1:length(Dep)-1]
    #Inputs[[s]]$B.Bmsy_DBSRA=B.Bmsy[1:length(Dep)-1]
    #Inputs[[s]]$F.Fmsy_DBSRA=F.Fmsy[1:length(Dep)-1]
    
    #Inputs[[s]]$Depletion_DBSRA=mean(Dep[(length(Dep)-4):length(Dep)],na.rm=T)
    # Inputs[[s]]$B.Bmsy_DBSRA=mean(B.Bmsy[(length(B.Bmsy)-4):length(B.Bmsy)],na.rm=T)
    # Inputs[[s]]$F.Fmsy_DBSRA=mean(F.Fmsy[(length(F.Fmsy)-4):length(F.Fmsy)],na.rm=T)
    
    rm(Dep,B.Bmsy,F.Fmsy,out.dummy)
    
    
    #CMSY
    print(paste("___________","CMSY ","______ s=",s))
    
    #priors
    Scens=List.sp[[i]]$Sens.test$CMSY%>%filter(Scenario=='S1')
    RES=RESILIENCE[[i]]
    if(Scens$r.prob.min==0)
    {
      r.range=NA
      k.range=NA
      Bf.low=NA
      Bf.hi=NA
    }else
    {
      Mn=min(store.species.r[[i]]$mean,Max.r.value)  #some life history pars yield unrealistically high r
      r.range=quantile(rnorm(1e3,
                             mean=Mn,
                             sd=store.species.r[[i]]$sd),
                       probs=c(Scens$r.prob.min,Scens$r.prob.max))
      k.range=c(Klow,Kup)
    }
    Proc.error=Scens$Proc.error
    
    #run model
    out.dummy=apply.CMSY(year=Yrs,
                         catch=Ktch,
                         r.range=r.range,
                         k.range=k.range,
                         Bo.low=b1k.low,
                         Bo.hi=b1k.up,
                         Bf.low=Btklow,
                         Bf.hi=Btkup,
                         outfile=paste(this.wd,'Appendix_fit',sep='/'),
                         nsims=Scens$Sims,
                         Proc.error=Proc.error,
                         RES=RES)
    
    Dep=apply(out.dummy$output$Depletion.traj,2,median)
    #B.Bmsy=apply(out.dummy$output$B.Bmsy,2,median)
    #F.Fmsy=apply(out.dummy$output$F.Fmsy,2,median)
    
    Inputs[[s]]$Depletion_CMSY=Dep[1:length(Dep)-1]
    #Inputs[[s]]$B.Bmsy_CMSY=B.Bmsy[1:length(Dep)-1]
    #Inputs[[s]]$F.Fmsy_CMSY=F.Fmsy[1:length(Dep)-1]
    
    # Inputs[[s]]$Depletion_CMSY=mean(Dep[(length(Dep)-4):length(Dep)],na.rm=T)
    # Inputs[[s]]$B.Bmsy_CMSY=mean(B.Bmsy[(length(B.Bmsy)-4):length(B.Bmsy)],na.rm=T)
    # Inputs[[s]]$F.Fmsy_CMSY=mean(F.Fmsy[(length(F.Fmsy)-4):length(F.Fmsy)],na.rm=T)
    
    rm(Dep,B.Bmsy,F.Fmsy,out.dummy)
    
    
    #JABBA
    print(paste("___________","JABBA ","______ s=",s))
    
    #priors
    Scens=List.sp[[i]]$Sens.test$JABBA%>%filter(Scenario=='S1')
    bmsyk.mean=Scens$bmsyk
    Proc.error=Scens$Proc.error
    r.CV.multi=Scens$r.CV.multiplier
    K.prior=c(Klow*20,Scens$K.CV)
    
    Mn=min(store.species.r[[i]]$mean,Max.r.value)  #some life history pars yield unrealistically high r
    r.prior=c(Mn, (store.species.r[[i]]$sd/Mn)*Scens$r.CV.multiplier)
    Ktch.CV=Scens$Ktch.CV
    Bint=runif(1000,b1k.low,b1k.up)
    Bint.mean=mean(Bint)
    Bint.CV=sd(Bint)/Bint.mean
    psi.prior=c(Bint.mean,Bint.CV)
    
    Bfin=runif(1000,Btklow,Btkup)
    Bfin.mean=mean(Bfin)          
    Bfin.CV=sd(Bfin)/Bfin.mean
    b.prior=c(Bfin.mean,Bfin.CV,max(Yrs),"bk")
    
    Rdist = "lnorm"
    Kdist="lnorm"  
    PsiDist='lnorm'
    
    #Put inputs together
    Ktch=ifelse(Ktch<1,1,Ktch)   #Jabba fails if catches are negligible
    input=list(Ktch=data.frame(Year=Yrs,Total=Ktch),
               MDL="Schaefer",
               Ktch.CV=Ktch.CV,
               ASS=names(List.sp)[i],
               Rdist = Rdist,
               Rprior = r.prior,
               r.CV.multiplier=r.CV.multi,
               Kdist= Kdist,
               Kprior=K.prior,    
               PsiDist= PsiDist,
               Psiprior=psi.prior,   
               Bprior=b.prior,    
               BMSYK=bmsyk.mean)  
    
    #run model
    THIN=5
    CHAINS=2
    BURNIN=min(0.15*Scens$Sims,5000)
    out.dummy=apply.JABBA(Ktch=input$Ktch,
                          CPUE=NULL,
                          CPUE.SE=NULL,
                          MDL=input$MDL,
                          Ktch.CV=input$Ktch.CV,
                          ASS=input$ASS,
                          Rdist = input$Rdist,
                          Rprior = input$Rprior,
                          Kdist=input$Kdist,
                          Kprior=input$Kprior,
                          PsiDist=input$PsiDist,
                          Psiprior=input$Psiprior,
                          Bprior=input$Bprior,
                          BMSYK=input$BMSYK,
                          output.dir=this.wd,
                          outfile="Appendix_fit",
                          Sims=Scens$Sims,
                          Proc.error.JABBA=Proc.error,
                          thinning = THIN,
                          nchains = CHAINS,
                          burn.in= BURNIN)
    
    Dep=apply(out.dummy$posteriors$P,2,median)
    #B.Bmsy=out.dummy$timeseries[, , "BBmsy"][,1]
    #F.Fmsy=out.dummy$timeseries[, , "FFmsy"][,1]
    
    Inputs[[s]]$Depletion_JABBA=Dep
    #Inputs[[s]]$B.Bmsy_JABBA=B.Bmsy
    #Inputs[[s]]$F.Fmsy_JABBA=F.Fmsy
    
    #Inputs[[s]]$Depletion_JABBA=mean(Dep[(length(Dep)-4):length(Dep)],na.rm=T)
    # Inputs[[s]]$B.Bmsy_JABBA=mean(B.Bmsy[(length(B.Bmsy)-4):length(B.Bmsy)],na.rm=T)
    # Inputs[[s]]$F.Fmsy_JABBA=mean(F.Fmsy[(length(F.Fmsy)-4):length(F.Fmsy)],na.rm=T)
    
    rm(Dep,B.Bmsy,F.Fmsy,out.dummy)
    
    drop.files=list.files(this.wd, full.names = TRUE)
    drop.files=drop.files[-grep('Recover.pars',drop.files)]
    do.call(file.remove, list(drop.files))
  }})
  
  #Compare COMs estimates
  le.cols=c('grey40',"#F8766D", "#00BFC4", "#7CAE00")
  names(le.cols)=c('Operation model','DBSRA','CMSY','JABBA')
  Dumi=do.call(rbind,Inputs)%>%
        mutate(Grup=paste(HRscen,'&' ,capitalize(r.group),'r'),
               Grup=factor(Grup,levels=c("HarvestRate1 & Low r","HarvestRate1 & Medium r","HarvestRate1 & High r",
                                         "HarvestRate2 & Low r","HarvestRate2 & Medium r","HarvestRate2 & High r",
                                         "HarvestRate3 & Low r","HarvestRate3 & Medium r","HarvestRate3 & High r",
                                         "HarvestRate4 & Low r","HarvestRate4 & Medium r","HarvestRate4 & High r")))
  fn.out.poly=function(In.dep)
  {
    CI_lower=Dumi%>%
      filter(Init.dep==In.dep)%>%
      dplyr::select(Year,Grup,Depletion,Depletion_DBSRA,Depletion_CMSY,Depletion_JABBA)%>%
      group_by(Year,Grup)%>%
      summarise(across(everything(), function(x)quantile(x,probs=0.001)))
    CI_upper=Dumi%>%
      filter(Init.dep==In.dep)%>%
      dplyr::select(Year,Grup,Depletion,Depletion_DBSRA,Depletion_CMSY,Depletion_JABBA)%>%
      group_by(Year,Grup)%>%
      summarise(across(everything(), function(x)quantile(x,probs=0.999)))
    CI_lower%>%
      left_join(CI_upper,by=c('Year','Grup'))%>%
      arrange(Year,Grup)%>%
      ggplot(aes(Year,Depletion.x))+
      facet_wrap(~Grup,ncol=3)+
      theme(legend.position = 'none')+
      geom_ribbon(aes(ymin=Depletion_DBSRA.x,ymax=Depletion_DBSRA.y,fill="DBSRA"),alpha=0.5)+
      geom_ribbon(aes(ymin=Depletion_CMSY.x,ymax=Depletion_CMSY.y,fill="CMSY"),alpha=0.5)+
      geom_ribbon(aes(ymin=Depletion_JABBA.x,ymax=Depletion_JABBA.y,fill="JABBA"),alpha=0.5)+
      geom_ribbon(aes(ymin=Depletion.x,ymax=Depletion.y,fill='Operation model'),alpha=0.8)+
      theme_PA()+ylab('Depletion')+
      scale_fill_manual(values = le.cols)+
      theme(legend.position = 'top',
            legend.title = element_blank())
    
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),
                 paste('Compare COMs performance_Init.dep_',In.dep,'.tiff',sep=''),sep=''),
           width = 10,height = 10,compression = "lzw")
  }
  fn.out.poly(In.dep=1)
  fn.out.poly(In.dep=0.7)
  
  
  #Display harvest rates
  le.cols=c("#F8766D", "#00BFC4", "#7CAE00","#C77CFF")
  names(le.cols)=paste('HarvestRate',1:4,sep='')
  
  Inputd%>%
    mutate(HarvestRate1=HarvestRate1+rnorm(nrow(Inputd),0,.01),
           HarvestRate2=HarvestRate2+rnorm(nrow(Inputd),0,.01),
           HarvestRate3=HarvestRate3+rnorm(nrow(Inputd),0,.01),
           HarvestRate4=HarvestRate4+rnorm(nrow(Inputd),0,.01),
           HarvestRate1=ifelse(HarvestRate1<0,0,HarvestRate1),
           HarvestRate2=ifelse(HarvestRate2<0,0,HarvestRate2),
           HarvestRate3=ifelse(HarvestRate3<0,0,HarvestRate3),
           HarvestRate4=ifelse(HarvestRate4<0,0,HarvestRate4))%>%
    ggplot(aes(YRS,HarvestRate1))+
    geom_line(aes(YRS,HarvestRate1,color='HarvestRate1'),size=1.1)+
    geom_line(aes(YRS,HarvestRate2,color='HarvestRate2'),size=1.1)+
    geom_line(aes(YRS,HarvestRate3,color='HarvestRate3'),size=1.1)+
    geom_line(aes(YRS,HarvestRate4,color='HarvestRate4'),size=1.1)+
    scale_color_manual(name='Legend',values = le.cols)+
    theme_PA(axs.t.siz=14,axs.T.siz=16,leg.siz=16)+
    theme(legend.position = 'top',
          legend.title = element_blank())+
    ylab(expression('Harvest rate'~(y^-1)))+xlab('Year')
  ggsave(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'Harvest rates.tiff',sep=''),
         width = 10,height = 10,compression = "lzw")
  
  
  #3. Calculate COMs' weights
  COMs.weight=Inputs
  for(s in 1:length(Inputs))
  {
    dummy=Inputs[[s]]
    dummy=dummy[(nrow(dummy)-4):nrow(dummy),]%>%
      mutate(Depletion=mean(Depletion,na.rm=T),
             Depletion_DBSRA=mean(Depletion_DBSRA,na.rm=T),
             Depletion_CMSY=mean(Depletion_CMSY,na.rm=T),
             Depletion_JABBA=mean(Depletion_JABBA,na.rm=T))
    COMs.weight[[s]]=dummy[nrow(dummy),]%>%
                      dplyr::select(r.group,r,K,Init.dep,HRscen,Depletion,
                                    Depletion_DBSRA,Depletion_CMSY,Depletion_JABBA)
    rm(dummy)
  }
  COMs.weight=do.call(rbind,COMs.weight)
  rownames(COMs.weight)=NULL
  
  COMs.weight=COMs.weight%>%
    mutate(JABBA.error=abs((Depletion_JABBA-Depletion)/Depletion),
           CMSY.error=abs((Depletion_CMSY-Depletion)/Depletion),
           DBSRA.error=abs((Depletion_DBSRA-Depletion)/Depletion),
           JABBA.weight=1/JABBA.error,
           CMSY.weight=1/CMSY.error,
           DBSRA.weight=1/DBSRA.error,
           Sum.weight=JABBA.weight+CMSY.weight+DBSRA.weight,
           JABBA.weight=JABBA.weight/Sum.weight,
           CMSY.weight=CMSY.weight/Sum.weight,
           DBSRA.weight=DBSRA.weight/Sum.weight)
  
  #overall weight
  Overall.weight=data.frame(Model=c('JABBA','CMSY','DBSRA'),
                            Weight=c(median(COMs.weight$JABBA.weight),
                                     median(COMs.weight$CMSY.weight),
                                     median(COMs.weight$DBSRA.weight)))%>%
                            mutate(Tot.w=sum(Weight),
                                   Weight=Weight/Tot.w)%>%
                            dplyr::select(-Tot.w)
  write.csv(Overall.weight,
            paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weight.csv',sep=''),row.names = F) 

  #weight by r group
  by.r.group.weight=vector('list',length(r.groups))
  names(by.r.group.weight)=names(r.groups)
  for(x in 1:length(by.r.group.weight))
  {
    xx=names(by.r.group.weight)[x]
    dummy=COMs.weight%>%filter(r.group==xx)
    by.r.group.weight[[x]]=data.frame(Model=c('JABBA','CMSY','DBSRA'),
                                 Weight=c(median(dummy$JABBA.weight),
                                          median(dummy$CMSY.weight),
                                          median(dummy$DBSRA.weight)))%>%
                              mutate(Tot.w=sum(Weight),
                                     Weight=Weight/Tot.w)%>%
                              dplyr::select(-Tot.w)%>%
                              mutate(r.group=xx)
  }
  write.csv(do.call(rbind,by.r.group.weight),
            paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weight_by.r.group.csv',sep=''),row.names = F) 

  #weight by r group and exploitation history
  hr.group=list(c('HarvestRate1','HarvestRate2'),c('HarvestRate3','HarvestRate4'))
  by.r.group_hr.group.weight=rep(hr.group,length(r.groups))
  names(by.r.group_hr.group.weight)=rep(names(r.groups),length(hr.group))
  for(x in 1:length(by.r.group_hr.group.weight))
  {
    xx=names(by.r.group_hr.group.weight)[x]
    dummy=COMs.weight%>%filter(r.group==xx & HRscen%in%by.r.group_hr.group.weight[[x]])
    by.r.group_hr.group.weight[[x]]=data.frame(Model=c('JABBA','CMSY','DBSRA'),
                                               Weight=c(median(dummy$JABBA.weight),
                                                        median(dummy$CMSY.weight),
                                                        median(dummy$DBSRA.weight)))%>%
      mutate(Tot.w=sum(Weight),
             Weight=Weight/Tot.w)%>%
      dplyr::select(-Tot.w)%>%
      mutate(r.group=xx,
             HRscen=paste(by.r.group_hr.group.weight[[x]],collapse='_'))
  }
  write.csv(do.call(rbind,by.r.group_hr.group.weight),
            paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weight_by.r.group_hr.group.csv',sep=''),row.names = F) 
  
  
}
if(!do.ensemble.simulations)
{
  COM_weight_overall=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weight.csv',sep=''))
  COM_weight_by.r.group=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weight_by.r.group.csv',sep=''))
  COM_weight_by.r.group.h.group=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weight_by.r.group_hr.group.csv',sep=''))
  
}

    #19.2.2. sample each accepted COM proportionally to weights  #ACA. Re do weighting and redo all outputs
mod.average=function(dd,Weights)
{
  yrs=as.numeric(gsub("^.*X","",names(dd$CMSY)))
  NN=sum(sapply(dd,nrow))
  for(w in 1:length(dd))
  {
    WEI=Weights%>%filter(Model==names(dd)[w])%>%pull(Weight)
    Size=round(WEI*NN)
    dd[[w]]=dd[[w]][sample(1:nrow(dd[[w]]),Size,replace = T),]
    colnames(dd[[w]])=1:ncol(dd[[w]])
  }
  dd=do.call(rbind,dd)
  Median=apply(dd,2,median,na.rm=T)
  Lower=apply(dd,2,function(x) quantile(x,probs=0.025,na.rm=T))
  Upper=apply(dd,2,function(x) quantile(x,probs=0.975,na.rm=T)) 
  
  return(data.frame(year=yrs,median=Median,lower.95=Lower,upper.95=Upper))
}

#overall weight
Mod.AV_depletion=vector('list',length(List.sp))
names(Mod.AV_depletion)=names(List.sp)
Mod.AV_B.Bmsy=Mod.AV_F.Fmsy=Mod.AV_depletion
for(i in 1:N.sp)  
{
  dis.r.range=do.call(rbind,r.groups)%>%   
    data.frame%>%
    mutate(Get=between(r.list[i],.[[1]] , .[[2]]))%>%
    filter(Get==TRUE)
  Wei=COM_weight_by.r.group%>%
      filter(r.group==row.names(dis.r.range))%>%
    dplyr::select(-r.group)
  
  Mod.AV_depletion[[i]]=mod.average(dd=map(Catch_only, ~.x$ensemble[[i]]$Depletion),
                                    Weights=Wei)
  Mod.AV_B.Bmsy[[i]]=mod.average(dd=map(Catch_only, ~.x$ensemble[[i]]$B.Bmsy),
                                 Weights=Wei)
  Mod.AV_F.Fmsy[[i]]=mod.average(dd=map(Catch_only, ~.x$ensemble[[i]]$F.Fmsy),   
                                 Weights=Wei)
  rm(Wei)
}
Mod.AV=list(rel.biom=Mod.AV_depletion,B.Bmsy=Mod.AV_B.Bmsy,F.Fmsy=Mod.AV_F.Fmsy)
rm(Mod.AV_depletion,Mod.AV_B.Bmsy,Mod.AV_F.Fmsy)



#---20. Catch-only assessments. Generate outputs --------------------------------------

  #20.1 Table of scenarios
for(l in 1: length(Lista.sp.outputs))
{
  for(w in 1:length(Catch_only))
  {
    dummy=Catch_only[[w]]$sens.table
    dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
    
    write.csv(do.call(rbind,dummy)%>%relocate(Species),
              paste(Rar.path,paste('Table 2. Catch only scenarios_',
                                   names(Catch_only)[w],'_',names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
              row.names = F)
  }
}


  #20.2 Table of parameter estimates by species and catch-only assessment method
for(l in 1: length(Lista.sp.outputs))
{
  for(w in 1:length(Catch_only))
  {
    dummy=Catch_only[[w]]$estimates
    dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
    dummy=do.call(rbind,dummy)%>%
      rownames_to_column(var = "Species")%>%
      mutate(Species=capitalize(str_extract(Species, "[^.]+")))%>%
      relocate(Species)
    write.csv(dummy,paste(Rar.path,paste('Table 3. Catch only estimates_',
                         names(Catch_only)[w],'_',names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
              row.names = F)
  }
}


  #20.3. Time series   
#note: JABBA only outputs 95% CI so only report 95% for all COMs
make_plot <- function(da)
{
  p=da%>%
    ggplot(aes(year, median))+
    geom_line(size=1.1)  +
    geom_ribbon(aes(ymin = lower.95, ymax = upper.95), alpha = 0.3,fill='grey60') +
    facet_wrap(~Scenario,ncol=nfacets)+
    theme_PA(axs.T.siz=AXST,axs.t.siz=AXSt,strx.siz=STRs)+
    theme(plot.title =element_text(hjust = 0.5))+
    ylab(YLAB)+xlab(XLAB)+
    ylim(0,max(da$upper.95))+
    theme(panel.spacing=unit(InrMarg,"lines"),
          plot.margin = unit(c(.5, -.2, 0, 0), "cm"))
  if(any(!is.na(da$upper.50))) p=p+geom_ribbon(aes(ymin = lower.50, ymax = upper.50), alpha = 0.1)
  if(!is.null(Hline)) p=p+geom_hline(yintercept=Hline, size=1.05,alpha=0.35,
                                     color=rep(c('forestgreen','orange','red'),length(unique(da$Scenario))))
  if(addKtch)
  {
    coeff=max(da$Catch)
    da$Ktch.scaled=da$Catch/coeff
    p=p+ 
      geom_line(data=da, aes(x=year, y=Ktch.scaled),size=1.1,color='dodgerblue4',alpha=0.5,linetype ='dashed')+
      scale_y_continuous(sec.axis = sec_axis(~.*coeff, name=""))
  }
}
fn.ribbon=function(Dat,YLAB,XLAB,Titl,Hline,addKtch,nfacets=1,AXST=14,AXSt=12,STRs=14,InrMarg=.25,dropTitl=FALSE)
{
  data2 <- split(Dat, Dat$Scenario)
  p_lst <- lapply(data2, make_plot)
  figure <- ggarrange(plotlist=p_lst,ncol=nfacets,nrow=ceiling(length(p_lst)/nfacets),
                      common.legend = FALSE)
  if(!dropTitl) figure <- annotate_figure(figure,top=text_grob(Titl, size=20))
  
  return(figure) 
}
fn.plot.ktch.only.timeseries=function(d,sp,Type,YLAB,add.50=FALSE,add.sp.nm=FALSE)
{
  mods=names(d)
  store.plots=vector('list',length(mods))
  names(store.plots)=mods
  store.probs=store.plots
  id=match(sp,Keep.species)
  for(m in 1:length(mods))
  {
    Hline=NULL
    addKtch=FALSE
    if(Type=='F.series') Var=d[[m]]$f.series[[id]]
    if(Type=='B.Bmsy')   Var=d[[m]]$B.Bmsy[[id]]
    if(Type=='F.Fmsy')   Var=d[[m]]$F.Fmsy[[id]]
    
    if(Type=='Depletion')
    {
      addKtch=TRUE
      Var=d[[m]]$rel.biom[[id]]
      str.prob=compact(d[[m]]$probs.rel.biom[[id]])   
      Hline=str.prob[[1]]$Reference.points$Value
      store.probs[[m]]=do.call(rbind,sapply(str.prob,'[',1))%>%  
                          mutate(Model=names(store.probs)[m])
    }
    store.plots[[m]]=fn.ribbon(Dat=Var,
                               YLAB='',
                               XLAB="",
                               Titl=names(d)[m],
                               Hline=Hline,
                               addKtch=addKtch)
    
  }
  figure=ggarrange(plotlist=store.plots, ncol = 3,common.legend=TRUE)
  if(add.sp.nm) figure=figure+theme(plot.margin = margin(1,0,0,0, "cm"))
  figure=annotate_figure(figure,
                         bottom = text_grob('Financial year',size=26,vjust =-0.15),
                         left = text_grob(YLAB,size=26,rot = 90,vjust=0.8))
  if(add.sp.nm) figure=annotate_figure(figure,
                                       fig.lab=capitalize(sp),
                                       fig.lab.pos='top.left',
                                       fig.lab.size=28)
  if(Type=='Depletion')
  {
    figure=annotate_figure(figure,
                           right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                           color ='dodgerblue4',vjust=-.2))
  }
  
  print(figure)      
  store.probs=do.call(rbind,store.probs)
  rownames(store.probs)=NULL
  return(list(store.probs=store.probs,Ref.points=str.prob[[1]]$Reference.points))
}

    #20.3.1 Display all scenarios for each species
      #20.3.1.1 Relative biomass (i.e. Depletion)
Ref.points=vector('list',N.sp)
names(Ref.points)=Keep.species
for(i in 1:N.sp)
{
  print(paste("Relative biomass plot -----",Keep.species[i]))
  a=fn.plot.ktch.only.timeseries(d=Catch_only,
                                 sp=Keep.species[i],
                                 Type='Depletion',
                                 YLAB='Relative biomass')
  #export graph
  ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
               capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_time_series_relative_biomass.tiff",sep=''),
         width = 12,height = 10,compression = "lzw")
  
  #export current depletion probabilities
  write.csv(a$store.probs%>%
              spread(Model,Probability)%>%
              mutate(Species=Keep.species[i],
                     Range=factor(Range,levels=c("<lim","lim–thr","thr–tar",">tar")))%>%
              arrange(Range),
            paste(handl_OneDrive("Analyses/Population dynamics/1."),
                  capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_current_depletion.csv",sep=''),
            row.names = F)
  Ref.points[[i]]=a$Ref.points
}

      #20.3.1.2 Fishing mortality
do.F.series=FALSE
if(do.F.series)
{
  for(i in 1:N.sp)
  {
    print(paste("Fishing mortality plot -----",Keep.species[i]))
    a=fn.plot.ktch.only.timeseries(d=Catch_only,
                                   sp=Keep.species[i],
                                   Type='F.series',
                                   YLAB='Fishing mortality')
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_time_series_fishing_mortality.tiff",sep=''),
           width = 10,height = 10,compression = "lzw")
    
  }
}

      #20.3.1.3 B over Bmsy
do.B.over.Bmsy.series=FALSE
if(do.B.over.Bmsy.series)
{
  for(i in 1:N.sp)
  {
    print(paste("B over Bmsy plot -----",Keep.species[i]))
    a=fn.plot.ktch.only.timeseries(d=Catch_only,
                                   sp=Keep.species[i],
                                   Type='B.Bmsy',
                                   YLAB='B/Bmsy')
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_time_series_B_Bmsy.tiff",sep=''),
           width = 10,height = 10,compression = "lzw")
    
  }
}

      #20.3.1.4 F over Fmsy
do.F.over.Fmsy.series=FALSE
if(do.F.over.Fmsy.series)
{
  for(i in 1:N.sp)
  {
    print(paste("F over Fmsy plot -----",Keep.species[i]))
    a=fn.plot.ktch.only.timeseries(d=Catch_only,
                                   sp=Keep.species[i],
                                   Type='F.Fmsy',
                                   YLAB='F/Fmsy')
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_time_series_F_Fmsy.tiff",sep=''),
           width = 10,height = 10,compression = "lzw")
    
  }
}

  #20.3.2 Display Scenario 1 for combined species: each COM model as an Appendix 
fn.plot.ktch.only.timeseries_combined_Appendix=function(this.sp,d,YLAB,NM,Type,InnerMargin)
{
  mods=names(d)
  if(length(this.sp)>8)
  {
    for(m in 1:length(mods))
    {
      id=match(this.sp,names(d[[m]]$rel.biom))
      Hline=NULL
      addKtch=FALSE
      if(Type=='F.series')  Var=d[[m]]$f.series[id]
      if(Type=='B.Bmsy')    Var=d[[m]]$B.Bmsy[id]
      if(Type=='F.Fmsy')    Var=d[[m]]$F.Fmsy[id]
      if(Type=='Depletion')
      {
        Var=d[[m]]$rel.biom[id]
        addKtch=TRUE
        Hline=d[[m]]$probs.rel.biom[[id[1]]][[1]]$Reference.points$Value
      }
      figure=fn.ribbon(Dat=do.call(rbind,Var)%>%
                         rownames_to_column()%>%
                         filter(Scenario=='S1')%>%
                         mutate(Scenario=capitalize(word(rowname,1,sep = "\\."))),
                       YLAB='',
                       XLAB="",
                       Titl=mods[m],
                       Hline=Hline,
                       addKtch=addKtch,
                       nfacets=round(length(this.sp)/5),
                       AXST=15,AXSt=13,STRs=12,
                       InrMarg=InnerMargin,
                       dropTitl=TRUE)
      
      figure=annotate_figure(figure,
                             bottom = text_grob('Financial year',size=26,vjust =-0.15),
                             left = text_grob(YLAB,size=26,rot = 90,vjust=0.8))+
        theme(plot.margin = unit(c(.1,.1,0,0), "cm"))
      if(Type=='Depletion')
      {
        figure=annotate_figure(figure,
                               right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                               color ='dodgerblue4',vjust=-.2))
      }
      print(figure)
      ggsave(paste(Rar.path,'/Relative.biomass_catch.only_',NM,'_',mods[m],'_Appendix','.tiff',sep=''),
             width = 11,height = 11,compression = "lzw")
    }
    
  }else
  {
    plot.list=vector('list',length(mods))
    for(m in 1:length(mods))
    {
      id=match(this.sp,names(d[[m]]$rel.biom))
      Hline=NULL
      addKtch=FALSE
      if(Type=='F.series')  Var=d[[m]]$f.series[id]
      if(Type=='B.Bmsy')    Var=d[[m]]$B.Bmsy[id]
      if(Type=='F.Fmsy')    Var=d[[m]]$F.Fmsy[id]
      if(Type=='Depletion')
      {
        Var=d[[m]]$rel.biom[id]
        addKtch=TRUE
        Hline=d[[m]]$probs.rel.biom[[id[1]]][[1]]$Reference.points$Value
      }
      plot.list[[m]]=fn.ribbon(Dat=do.call(rbind,Var)%>%
                                     filter(Scenario=='S1')%>%
                                     rownames_to_column()%>%
                                     mutate(Scenario=capitalize(word(rowname,1,sep = "\\."))),
                               YLAB='',
                               XLAB="",
                               Titl=mods[m],
                               Hline=Hline,
                               addKtch=addKtch,
                               nfacets=1,
                               AXST=16,AXSt=14,STRs=16,
                               InrMarg=InnerMargin)
    }
    
    figure=ggarrange(plotlist=plot.list,ncol=length(mods), common.legend=TRUE)+
                      theme(plot.margin = margin(0,0,0,0, "cm"))
    figure=annotate_figure(figure,
                           bottom = text_grob('Financial year',size=26,vjust =-0.15),
                           left = text_grob(YLAB,size=26,rot = 90,vjust=0.8))
    if(Type=='Depletion')
    {
      figure=annotate_figure(figure,
                             right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                             color ='dodgerblue4',vjust=-.2))
    }
    print(figure)
    ggsave(paste(Rar.path,'/Relative.biomass_catch.only_',NM,'_Appendix','.tiff',sep=''),
           width = 15,height = 10,compression = "lzw")
    
  }
}
    # Relative biomass (i.e. Depletion) 
    #figure
for(l in 1:length(Lista.sp.outputs))
{
  if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
  fn.plot.ktch.only.timeseries_combined_Appendix(this.sp=Lista.sp.outputs[[l]],
                                        d=Catch_only,
                                        YLAB="Relative biomass",
                                        NM=names(Lista.sp.outputs)[l],
                                        Type="Depletion",
                                        InnerMargin=InMar)
}
    #table
for(l in 1:length(Lista.sp.outputs))
{
  dummy.mod=vector('list',length(Catch_only))
  for(m in 1:length(Catch_only))
  {
    str.prob=Catch_only[[m]]$probs.rel.biom
    str.prob=str.prob[match(Lista.sp.outputs[[l]],names(str.prob))]
    
    dummy=vector('list',length =length(str.prob))
    for(d in 1:length(dummy))
    {
      dummy[[d]]=str.prob[[d]][[1]]$probs%>%
                  mutate(Species=capitalize(names(str.prob)[d]))
    }
    dummy.mod[[m]]=do.call(rbind,dummy)%>%
      mutate(Model=names(Catch_only)[m])
  }
  write.csv(do.call(rbind,dummy.mod)%>%
              mutate(Range=factor(Range,levels=c("<lim","lim–thr","thr–tar",">tar")))%>%
              spread(Species,Probability)%>%
              arrange(Range),
            paste(Rar.path,'/Table 4. Current.depletion_catch.only_',names(Lista.sp.outputs)[l],'_Appendix','.csv',sep=''),
            row.names=F)
  rm(dummy.mod)
}

  #20.3.3 Display Scenario 1 for combined species: COM model-average 
fn.plot.ktch.only.timeseries_combined=function(this.sp,d,YLAB,NM,Type,InnerMargin,RefPoint,Kach)
{
  id=match(this.sp,names(d$rel.biom))
  RefPoint=RefPoint[id]
  Kach=Kach[id]
  Kach=do.call(rbind,Kach)%>%
    filter(Scenario=='S1')%>%
    rownames_to_column()%>%
    mutate(Scenario=capitalize(word(rowname,1,sep = "\\.")))%>%
    dplyr::select(year,Catch,Scenario)
  
  
  Hline=NULL
  addKtch=FALSE
  if(Type=='F.series')  Var=d$f.series[id]
  if(Type=='B.Bmsy')    Var=d$B.Bmsy[id]
  if(Type=='F.Fmsy')    Var=d$F.Fmsy[id]
  if(Type=='Depletion')
  {
    Var=d$rel.biom[id]
    addKtch=TRUE
  }
  
  if(length(this.sp)>8)  Nfast=round(length(this.sp)/5)
  if(length(this.sp)>2 & length(this.sp)<8) Nfast=2
  if(length(this.sp)<=2) Nfast=1
  
  figure=fn.ribbon(Dat=do.call(rbind,Var)%>%
                     rownames_to_column()%>%
                     mutate(Scenario=capitalize(word(rowname,1,sep = "\\.")))%>%
                     left_join(Kach,by=c('year','Scenario')),
                   YLAB='',
                   XLAB="",
                   Titl='',
                   Hline=RefPoint[[1]]$Value,
                   addKtch=addKtch,
                   nfacets=Nfast,
                   AXST=15,AXSt=13,STRs=15,
                   InrMarg=InnerMargin,
                   dropTitl=TRUE)
  
  figure=annotate_figure(figure,
                         bottom = text_grob('Financial year',size=26,vjust =-0.15),
                         left = text_grob(YLAB,size=26,rot = 90,vjust=0.8)) +
    theme(plot.margin = unit(c(.1,.1,0,0), "cm"))
  if(Type=='Depletion')
  {
    figure=annotate_figure(figure,
                           right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                           color ='dodgerblue4',vjust=0))
  }
  print(figure)
  ggsave(paste(Rar.path,'/Relative.biomass_catch.only_',NM,'.tiff',sep=''),
         width = 11,height = 11,compression = "lzw")
}
    # Relative biomass (i.e. Depletion) 
    #figure
for(l in 1:length(Lista.sp.outputs))
{
  if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
  fn.plot.ktch.only.timeseries_combined(this.sp=Lista.sp.outputs[[l]],
                                        d=Mod.AV,
                                        YLAB="Relative biomass",
                                        NM=names(Lista.sp.outputs)[l],
                                        Type="Depletion",
                                        InnerMargin=InMar,
                                        RefPoint=Ref.points,
                                        Kach=Catch_only$DBSRA$rel.biom)
}
    #table   
for(l in 1:length(Lista.sp.outputs))
{
  dummy.mod=vector('list',length(Catch_only))
  for(m in 1:length(Catch_only))
  {
    str.prob=Catch_only[[m]]$probs.rel.biom
    str.prob=str.prob[match(Lista.sp.outputs[[l]],names(str.prob))]
    
    dummy=vector('list',length =length(str.prob))
    for(d in 1:length(dummy))
    {
      dummy[[d]]=str.prob[[d]][[1]]$probs%>%
        mutate(Species=capitalize(names(str.prob)[d]))
    }
    dummy.mod[[m]]=do.call(rbind,dummy)%>%
      mutate(Model=names(Catch_only)[m])
  }
  
  dis.r.range=do.call(rbind,r.groups)%>%   
    data.frame%>%
    mutate(Get=between(r.list[i],.[[1]] , .[[2]]))%>%
    filter(Get==TRUE)
  Wei=COM_weight_by.r.group%>%
    filter(r.group==row.names(dis.r.range))%>%
    dplyr::select(-r.group)

  dd=do.call(rbind,dummy.mod)%>%
    left_join(Wei,by='Model')%>%
    group_by(Range,Scenario,Species)%>%
    summarise(Probability=weighted.mean(Probability,w=Weight))
  write.csv(dd%>%
              mutate(Range=factor(Range,levels=c("<lim","lim–thr","thr–tar",">tar")))%>%
              spread(Species,Probability)%>%
              arrange(Range),
            paste(Rar.path,'/Table 4. Current.depletion_catch.only_',names(Lista.sp.outputs)[l],'.csv',sep=''),
            row.names=F)
  rm(dummy.mod,Wei)
}


  #20.4. Kobe plots (Scenario 1)   
fn.get.Kobe.plot_appendix=function(d,sp,Scen='S1',add.sp.nm=FALSE)
{
  id=match(sp,Keep.species)
  
  #DBSRA
  dummy=d$DBSRA$B.Bmsy[[id]]%>%filter(Scenario==Scen)
  yrs=dummy%>%pull(year)
  Bmsy=dummy$median
  Fmsy=d$DBSRA$F.Fmsy[[id]]%>%filter(Scenario==Scen)%>%pull(median)
  p.DBSRA=kobePlot(f.traj=Fmsy[1:length(yrs)],
                   b.traj=Bmsy[1:length(yrs)],
                   Years=yrs,
                   Titl="DBSRA")
  rm(yrs,Fmsy,Bmsy,dummy)
  
  #CMSY
  if(CMSY.method=="Froese")
  {
    yrs=d$CMSY[[id]][[Scen]]$output$ref_ts$year
    Bmsy=d$CMSY[[id]][[Scen]]$output$ref_ts$bbmsy
    Fmsy=d$CMSY[[id]][[Scen]]$output$ref_ts$ffmsy
    
  }
  if(CMSY.method=="Haddon")
  {
    dummy=d$CMSY$B.Bmsy[[id]]%>%filter(Scenario==Scen)
    yrs=dummy%>%pull(year)
    Bmsy=dummy$median
    Fmsy=d$CMSY$F.Fmsy[[id]]%>%filter(Scenario==Scen)%>%pull(median)
  }
  p.CMSY=kobePlot(f.traj=Fmsy,
                  b.traj=Bmsy,
                  Years=yrs,
                  Titl="CMSY")
  rm(yrs,Fmsy,Bmsy,dummy)
  
  #JABBA
  dummy=d$JABBA$B.Bmsy[[id]]%>%filter(Scenario==Scen)
  yrs=dummy%>%pull(year)
  Bmsy=dummy$median
  Fmsy=d$JABBA$F.Fmsy[[id]]%>%filter(Scenario==Scen)%>%pull(median)
  p.JABBA=kobePlot(f.traj=Fmsy,
                  b.traj=Bmsy,
                  Years=yrs,
                  Titl="JABBA")
                 # Probs=data.frame(x=d$JABBA$Kobe.probs[[id]]$stock,   
                #                   y=d$JABBA$Kobe.probs[[id]]$harvest))
  
  #Combine plots
  plotlist=list(DBSRA=p.DBSRA+rremove("axis.title"),
                CMSY=p.CMSY+rremove("axis.title"),
                JABBA=p.JABBA+rremove("axis.title"))
  figure <- ggarrange(plotlist=plotlist,
                      ncol=1,nrow=3,common.legend = FALSE)
  
  if(add.sp.nm) figure=figure+theme(plot.margin = margin(1,0,0,0, "cm"))
  
  figure=annotate_figure(figure,
                         bottom = text_grob(expression(B/~B[MSY]), size=22),
                         left = text_grob(expression(F/~F[MSY]), rot = 90,size=22))
  if(add.sp.nm) figure=annotate_figure(figure,
                              fig.lab=capitalize(sp),
                              fig.lab.pos='top.left',
                              fig.lab.size=24)
    
  print(figure)
  return(plotlist)
}

    #20.4.1 by species
store.kobes=vector('list',N.sp)
names(store.kobes)=Keep.species
for(i in 1:N.sp)
{
  print(paste("Kobe plot -----",Keep.species[i]))
  store.kobes[[i]]=fn.get.Kobe.plot_appendix(d=Catch_only,
                                    sp=Keep.species[i])
  ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
               capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_Kobe_plot.tiff",sep=''),
         width = 9,height = 14, dpi = 300,compression = "lzw")
}

    #20.4.2 Display combined species: each COM model as an Appendix
do.this=FALSE
if(do.this)
{
  for(l in 1:length(Lista.sp.outputs))
  {
    if(length(Lista.sp.outputs[[l]])>8)
    {
      a=store.kobes[which(names(store.kobes)%in%Lista.sp.outputs[[l]])]
      for(m in 1:length(Catch_only))
      {
        figs=vector('list',length(a))
        for(f in 1:length(figs))
        {
          pp=a[[f]][[match(names(Catch_only)[m],names(a[[f]]))]]
          pp$labels$title=capitalize(names(a)[f])
          # figure <- ggarrange(plotlist=pp,
          #                     ncol=3,common.legend = FALSE)+
          #   theme(plot.margin = margin(1,0,0,0, "cm"))
          # figure=annotate_figure(figure,
          #                        fig.lab=capitalize(names(a)[f]),
          #                        fig.lab.pos='top.left',
          #                        fig.lab.size=22)
          figs[[f]]<-pp
          rm(pp)
        }
        figure=ggarrange(plotlist=figs)+
          theme(plot.margin = margin(0,.2,0,0, "cm"))
        annotate_figure(figure,
                        top=text_grob(names(Catch_only)[m], size=22),
                        bottom = text_grob(expression(B/~B[MSY]), size=22),
                        left = text_grob(expression(F/~F[MSY]), rot = 90,size=22))
        WID=13
        if(names(Catch_only)[m]=="JABBA") WID=16
        ggsave(paste(Rar.path,'/Kobe_plot_catch_only_',names(Lista.sp.outputs)[l],
                     '_',names(Catch_only)[m],'_Appendix.tiff',sep=''),
               width = WID,height = 12,compression = "lzw")
      }
    }else
    {
      a=store.kobes[which(names(store.kobes)%in%Lista.sp.outputs[[l]])]
      figs=vector('list',length(a))
      for(f in 1:length(figs))
      {
        figure <- ggarrange(plotlist=a[[f]],
                            ncol=3,common.legend = FALSE)+
          theme(plot.margin = margin(1,0,0,0, "cm"))
        figure=annotate_figure(figure,
                               fig.lab=capitalize(names(a)[f]),
                               fig.lab.pos='top.left',
                               fig.lab.size=22)
        figs[[f]]<-figure
      }
      figure=ggarrange(plotlist=figs,nrow=length(figs))
      annotate_figure(figure,
                      bottom = text_grob(expression(B/~B[MSY]), size=22),
                      left = text_grob(expression(F/~F[MSY]), rot = 90,size=22))
      ggsave(paste(Rar.path,'/Kobe_plot_catch_only_',names(Lista.sp.outputs)[l],'_Appendix.tiff',sep=''),
             width = 12,height = 12,compression = "lzw")
    }
  }
}


    #20.4.3 Display combined species: COM model-average 
fn.get.Kobe.plot=function(this.sp,d,NKOL,NRW)
{
  id=match(this.sp,names(d$rel.biom))
  Bmsy=d$B.Bmsy[id]
  Fmsy=d$F.Fmsy[id]
  plotlist=vector('list',length(Bmsy))
  for(x in 1:length(Bmsy))
  {
    plotlist[[x]]=kobePlot(f.traj=Fmsy[[x]]$median,
                           b.traj=Bmsy[[x]]$median,
                           Years=Bmsy[[x]]$year,
                           Titl=capitalize(names(Bmsy)[x]),
                           YrSize=6)+rremove("axis.title")
  }
  figure <- ggarrange(plotlist=plotlist,ncol=NKOL,nrow=NRW,common.legend = FALSE)
  figure=annotate_figure(figure,
                         bottom = text_grob(expression(B/~B[MSY]), size=22),
                         left = text_grob(expression(F/~F[MSY]), rot = 90,size=22))
  print(figure)
}
for(l in 1:length(Lista.sp.outputs))
{
  this.sp=Lista.sp.outputs[[l]]
  DIMS=n2mfrow(length(this.sp))
  NKOL=DIMS[2]
  NRW=DIMS[1]
  if(NKOL%in%3:4) WIZ=13
  if(NKOL==2) WIZ=11
  if(NKOL==1) WIZ=9
  fn.get.Kobe.plot(this.sp,d=Mod.AV,NKOL,NRW)
  ggsave(paste(Rar.path,'/Kobe_plot_catch_only_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = WIZ,height = 12,compression = "lzw")
}

  #20.5 Percentage acceptance rate for Scenario 1 (not applicable to JABBA)   
Per.accepted=vector('list',length(List.sp))
names(Per.accepted)=names(List.sp)
  #By species
for(i in 1: N.sp)
{
  moDs=names(Catch_only)[match(c("DBSRA","CMSY"),names(Catch_only))]  #not applicable to JABBA
  dummy=vector('list',length(moDs))
  names(dummy)=moDs
  for(m in 1:length(moDs))
  {
    dummy[[m]]=Catch_only[[m]]$accept.rate[[i]]%>%
                  filter(Scenario=='S1')%>%
                  mutate(Model=names(dummy)[m],
                         Species=Keep.species[i])
  }
  Per.accepted[[i]]=do.call(rbind,dummy)%>%
                      relocate(Species,Model,Scenario)
  rm(dummy)
}
  #By Lista.sp.outputs 
for(l in 1:length(Lista.sp.outputs))
{
  write.csv(do.call(rbind,Per.accepted[Lista.sp.outputs[[l]]])%>%
              data.frame%>%
              mutate(Species=capitalize(Species)),
            paste(Rar.path,'/Table 5. Per.acceptance.rate_catch.only_',
                  names(Lista.sp.outputs)[l],'.csv',sep=''),
            row.names = F)
}


  #20.6 MSY estimates by Lista.sp.outputs (Scenario 1)   #redundant, MSY output in #20.2
# for(l in 1: length(Lista.sp.outputs))
# {
#   dummy=vector('list',length(Catch_only))
#   for(w in 1:length(dummy))
#   {
#     a=Catch_only[[w]]$estimates[match(Lista.sp.outputs[[l]],names(Catch_only[[w]]$estimates))]
#     dummy[[w]]=do.call(rbind,a)%>%
#               filter(Scenario=='S1' & Parameter=='MSY')
#   }
#   write.csv(do.call(rbind,dummy),
#             paste(Rar.path,'/Table 6. MSY_estimates_catch.only_',
#                   names(Lista.sp.outputs)[l],'.csv',sep=''),
#             row.names = F)
# }





#---21. Spatio-temporal catch and effort. Reported TDGLF and NSF ----   
#note: bubble size is proportion of blocks fished out of maximum number of blocks fished for each species
  #get reported catch
dis.sp=All.species.names%>%filter(SNAME%in%unlist(Lista.sp.outputs))%>%pull(SPECIES)
dis.sp=c(dis.sp,19000)
Southern=fn.in(NM='Data.monthly.csv')%>%
  filter(!Shark.fishery=='non.shark.fishery' & SPECIES%in%dis.sp)%>%
  filter(!is.na(BLOCKX))%>%
  group_by(FINYEAR,BLOCKX,SPECIES)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
Northern=fn.in(NM='Data.monthly.NSF.csv')%>%
  filter(!Shark.fishery=='non.shark.fishery' & SPECIES%in%dis.sp)%>%
  filter(!is.na(BLOCKX))%>%
  group_by(FINYEAR,BLOCKX,SPECIES)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
Spatio.temp.dat=rbind(Southern,Northern)
rm(Southern,Northern)

do.reconstructed.spatial=FALSE
if(do.reconstructed.spatial)
{
  Spatio.temp.dat=rbind(fn.in(NM='recons_Data.monthly.csv')%>%
                          filter(FishCubeCode%in%c('JASDGDL','WCDGDL'))%>%
                          filter(!is.na(BLOCKX))%>%
                          dplyr::select(FINYEAR,SPECIES,BLOCKX,LIVEWT.c),
                        fn.in(NM='recons_Data.monthly.north.csv')%>%
                          filter(FishCubeCode%in%c('JANS','WANCS','OANCGC'))%>%
                          filter(!is.na(BLOCKX))%>%
                          dplyr::select(FINYEAR,SPECIES,BLOCKX,LIVEWT.c))%>%
    filter(LIVEWT.c>1) #reconstruction creates artificially large number of blocks with tiny catch, remove
  
}
  #get effort
Effort.monthly_blocks=fn.in(NM='Effort.monthly.csv')
Effort.daily_blocks=fn.in(NM='Effort.daily.csv') 
Effort.monthly.north_blocks=fn.in(NM='Effort.monthly.NSF.csv')
Effort.daily.north_blocks=fn.in(NM='Effort.daily.NSF.csv')

Relevant.yrs=sort(unique(KtCh.method$FINYEAR))
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

  #plot
fn.spatio.temp.catch.dist=function(d,Snames)
{
  this.sp=All.species.names%>%filter(SNAME%in%Snames)
  
  d1=d%>%
    filter(SPECIES%in%this.sp$SPECIES)%>%
    left_join(this.sp,by='SPECIES')%>%
    ungroup()
  Unik.yr=d1%>%
    distinct(FINYEAR)%>%
    arrange(FINYEAR)%>%
    mutate(year=as.numeric(substr(FINYEAR,1,4)))
  Unik.yr$id.x <- 1:nrow(Unik.yr) 
  
  Unik.sp=d1%>%
    distinct(SNAME)%>%
    arrange(SNAME)
  Unik.sp$id.y <- 1:nrow(Unik.sp) 
  
  Xlab=Unik.yr$year
  names(Xlab)=Unik.yr$id.x
  Xlab=Xlab[seq(1,length(Xlab),10)]
  
  Ylab=capitalize(Unik.sp$SNAME)
  names(Ylab)=Unik.sp$id.y
  
  Dummy=Effort_blocks%>%
    left_join(Unik.yr,by='FINYEAR')
  coeff=max(Dummy$Tot)

  d1=d1%>%
    count(FINYEAR,SNAME,BLOCKX)%>%
    group_by(FINYEAR,SNAME)%>%
    mutate(n=ifelse(n>0,1,0))%>%
    group_by(FINYEAR,SNAME)%>%
    summarise(n=sum(n,na.rm=T))%>%
    ungroup%>%
    group_by(SNAME)%>%
    mutate(Max.n=max(n),
           prop=n/Max.n)%>%
    ungroup%>%
    data.frame%>%
    left_join(Unik.yr,by='FINYEAR')%>%
    left_join(Unik.sp,by='SNAME')
  
  add.to.ylab=d1%>%distinct(SNAME,Max.n)%>%mutate(SNAME=capitalize(SNAME))
  add.to.ylab=add.to.ylab%>%
    left_join(data.frame(SNAME=Ylab,index=as.numeric(names(Ylab))),
                          by='SNAME')%>%
    arrange(index)%>%
    mutate(LBL=paste(SNAME,' (n=',Max.n,')',sep=''))
  
  Ylab=add.to.ylab$LBL
  names(Ylab)=add.to.ylab$index
    
  p=d1%>%
    ggplot(aes(x=id.x,y=id.y,colour=-prop)) +
    geom_point(aes(size=prop))+
    ylab('')+xlab('Financial year')+
    theme_PA(axs.t.siz=18,axs.T.siz=20)+
    theme(legend.position='none')+
    scale_x_continuous(labels=Xlab,breaks=as.numeric(names(Xlab)))+
    scale_y_continuous(labels=Ylab,breaks=as.numeric(names(Ylab)))+
    geom_line(data=Dummy, aes(y=max(Unik.sp$id.y)*Tot / coeff),size=1.5,color='black',
              alpha=0.8,linetype = 1)+
    scale_colour_gradient(low = "darkred", high = "darkgoldenrod1", na.value = NA)
  print(p)
  
  DD1=d1%>%
    dplyr::select(FINYEAR,prop,SNAME)%>%
    spread(FINYEAR,prop,fill=0)
  DD=as.matrix(DD1%>%dplyr::select(-SNAME))
  rownames(DD)=DD1$SNAME
  return(DD)
}
Store.spatial.temporal.ktch=Lista.sp.outputs[-match('additional.sp',names(Lista.sp.outputs))]
for(l in 1:length(Store.spatial.temporal.ktch))
{
  get.dis.sp=Lista.sp.outputs[[l]]
  if(names(Lista.sp.outputs)[l]=="Other.sp") get.dis.sp=c(get.dis.sp,'hammerheads')
  Store.spatial.temporal.ktch[[l]]=fn.spatio.temp.catch.dist(d=Spatio.temp.dat,
                                                             Snames=get.dis.sp)
  ggsave(paste(Rar.path,'/Spatio.temporal.catch_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = 11,height = 10,compression = "lzw")
}


#---22. Changes in observed mean length for TDGDLF----
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
        geom_violin(aes(fill = Finyear, group = Finyear.d), alpha = 0.2) +
        facet_wrap(~Mesh,scales='free')+
        stat_summary(fun = "mean",geom = "point",color = "red",size=2)+
        geom_smooth(color = "black",alpha=0.4, formula = my_formula, method = 'lm',se=TRUE)+ 
        stat_poly_eq(aes(label = paste("atop(", stat(eq.label),  ",", 
                                         paste(stat(adj.rr.label),stat(p.value.label), sep = "*\", \"*"), ")")),
                     formula = my_formula, parse = TRUE,
                     label.y = "top", label.x = "right", size = 4) +
        xlab("")+ylab("")+
        theme_PA(axs.T.siz=22,axs.t.siz=14,strx.siz=16)+
        theme(legend.position = "none",
              plot.title =element_text(size=17))+
        labs(title=NM)+
        xlim(XLIM)+
        ylim(quantile(d.list$FL,na.rm=T,probs = 0.01),quantile(d.list$FL,na.rm=T,probs = 0.975))  
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
for(l in 1:length(Store.spatial.temporal.ktch))
{
  figure=ggarrange(plotlist=fun.find.in.list(x=Change.mean.length[Lista.sp.outputs[[l]]]),
                   ncol=1)+
    theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
  annotate_figure(figure,
                  left = text_grob("Fork length (cm)", rot = 90,size=20,vjust=1),
                  bottom = text_grob("Financial year",size=20,vjust=-1))
  
  ggsave(paste(Rar.path,'/Changes in observed mean length_TDGDLF_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = 10,height = 10,compression = "lzw")
}



#---23. Changes in reported mean weight of individuals caught in the TDGDLF  ----
#Is there a strong declining trend in mean weights? (Leitao 2019)
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

Logbook=Logbook%>%filter(Finyear<=as.numeric(substr(Last.yr.ktch,1,4))) #keep only years with catch 

Change.mean.weight.catch=vector('list',length(Logbook.sp))
names(Change.mean.weight.catch)=Logbook.sp
fn.plt.mn.ktch.wght=function(d.list,NM,XLIM,show.data=FALSE)
{
  my_formula = y ~ x
  if(nrow(d.list)>0 & length(unique(d.list$finyear))>2)
  {
    p=d.list%>%
      mutate(Mesh=paste(Mesh,"inch"))%>%
      ggplot(aes(x = Finyear, y = Mean.wght)) 
   if(show.data) p=p+geom_point(position = position_jitter(seed = 1, width = 0.2),color='grey',alpha=.2)
    p=p+
      geom_violin(aes(fill = Finyear, group = Finyear.d), alpha = 0.2) + 
      facet_wrap(~Mesh,scales='free')+
      stat_summary(fun = "mean",geom = "point",color = "red",size=2)+
      geom_smooth(color = "black",alpha=0.4, formula = my_formula, method = 'lm',se=TRUE)+ 
      stat_poly_eq(aes(label = paste("atop(", stat(eq.label),  ",", 
                                     paste(stat(adj.rr.label),stat(p.value.label), sep = "*\", \"*"), ")")),
                   formula = my_formula, parse = TRUE,
                   label.y = "top", label.x = "right", size = 4.5) +
      xlab("")+ylab("")+
      ggtitle(NM)+xlim(XLIM)+
      theme_PA(axs.T.siz=22,axs.t.siz=14,str.siz=16)+
      theme(legend.position = "none",
            plot.title =element_text(size=17))+
      ylim(0,quantile(d.list$Mean.wght,na.rm=T,probs = 0.975))
       
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

for(l in 1:length(Store.spatial.temporal.ktch))
{
  figure=ggarrange(plotlist=fun.find.in.list(x=Change.mean.weight.catch[Lista.sp.outputs[[l]]]),
                   ncol=1)+
    theme(plot.margin = unit(c(0,0,0,0), "lines"))
  annotate_figure(figure,
                  left = text_grob("Mean total weight of caught individuals (kg)", rot = 90,size=20,vjust=1),
                  bottom = text_grob("Financial year",size=20,vjust=-0.5))
  
  ggsave(paste(Rar.path,'/Changes in mean weight caught individual_TDGDLF_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = 10,height = 10,compression = "lzw")
}



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


#ACA
#---24. Size-based Catch only method with dome-shaped selectivity --------------------------------------
#note: single-area, single sex, size-structured integrated model fitted to catch and size composition
#       Selectivity is assumed to be known.
#       Uncertainty derived from resampling variance-cov matrix
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/Apply Alex dynamic catch size model.R")) 
fn.extract.dat=function(STRING,Files) grep(paste(STRING,collapse="|"), Files, value=TRUE)


#24.1 TDGDLF size composition
size.catch.only_TDGDLF=vector('list',N.sp)
names(size.catch.only_TDGDLF)=Keep.species
system.time({for(l in 1: N.sp)  
{
  # Calculate F
  Outfile='TDGDLF' 
  this.size.comp=paste('Size_composition',c('West.6.5','West.7','Zone1.6.5','Zone1.7','Zone2.6.5','Zone2.7'),sep="_")
  outfile=paste(Outfile,'_histogram',sep='')
  
  #get size composition
  iid=Species.data[[l]][fn.extract.dat(this.size.comp,names(Species.data[[l]]))]
  if(length(iid)>0)
  {
    print(paste("Dynamic catch-only size based model with dome-shape selectivity for --",names(Species.data)[l],"---",Outfile))
    
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
      filter(n>=Min.annual.obs_catch.curve)%>%
      mutate(Keep=year)
    if(nrow(N.min)>0)
    {
      dummy=dummy%>%
        mutate(Keep=year)%>%
        filter(Keep%in%N.min$Keep)%>%
        mutate(TL=FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)     
      
      
      #1. Plot observed size frequency by year and mesh for years with minimum sample size
      if(grepl("TDGDLF",outfile))
      {
        p=dummy%>%
          ggplot( aes(x=TL, color=Mesh, fill=Mesh)) +
          geom_histogram(alpha=0.6, binwidth = TL.bins.cm)
        WHERE="top"
      }else
      {
        p=dummy%>%
          ggplot( aes(x=TL,color=year, fill=year)) +
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
                   "/1_Inputs/Visualise data/Size.comp_catch.curve.",outfile,".tiff",sep=''),
             width = 8,height = 8,compression = "lzw")
      
      #2. Fit model #ACA.
      
      #life history
      attach(List.sp[[l]]) 
      MaxAge = ceiling(mean(Max.age.F))
      Linf = Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL  #total length in cm
      vbK = Growth.F$k
      Lo =  Lzero*a_FL.to.TL+b_FL.to.TL   #total length in cm 
      MaxLen= 10*round(TLmax/10)
      LenInc=TL.bins.cm
      MatL50=TL.50.mat
      MatL95=TL.95.mat
      PropFemAtBirth=pup.sx.ratio
      wtlen_b=BwT
      wtlen_a=AwT
      detach(List.sp[[l]])
      NatMort=mean(colMeans(store.species.M[[l]],na.rm=T))
      NatMort_sd=sd(colMeans(store.species.M[[l]],na.rm=T))
      Steepness=store.species.steepness.S2[[l]]
      Steepness_sd=store.species.steepness[[l]]$sd
      lbnd = seq(0,MaxLen - LenInc, LenInc)
      ubnd = lbnd + LenInc
      midpt = lbnd + (LenInc/2)
      #assumed values
      UnfishRec=1
      CVLenAtAge = 0.1 
      SDGrowthRec = 20
      lnRecDev=0
      lnSigmaR=0.2 #assumed recruitment variation,assumed fairly low value for sharks
      #Missing: ask Alext about Init_F & InitRec
      
      #selectivity (combined meshes)
      SelAtLength=Selectivity.at.totalength[[l]]%>%               
        mutate(TL=TL)%>%
        filter( TL%in%midpt)%>%
        pull(Sel.combined)
      
      #total catch
      Katch=ktch.combined%>%
        filter(Name==names(size.catch.only_TDGDLF)[l])%>%
        ungroup()%>%
        dplyr::select(finyear,Tonnes)
      
      #size composition
      n.size.comp=dummy%>%
        group_by(year)%>%
        tally()
      Len_SimYr=match(n.size.comp$year,Katch$finyear)
      n_SimLen=n.size.comp$n
      
      add.dummy=data.frame(bin=midpt)
      ObsLenComp=vector('list',nrow(n.size.comp))
      for(o in 1:length(ObsLenComp))
      {
        x=dummy%>%
          filter(year==n.size.comp$year[o])%>%
          mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2)%>%
          group_by(bin)%>%
          tally()%>%
          full_join(add.dummy,by='bin')%>%
          arrange(bin)%>%
          mutate(n=ifelse(is.na(n),0,n))%>%
          filter(bin%in%midpt)
        xx=x$n
        names(xx)=x$bin
        ObsLenComp[[o]]=xx
      }
      ObsLenComp=do.call(rbind,ObsLenComp)%>%data.frame
      names(ObsLenComp)=str_remove(colnames(ObsLenComp), "[X]")
      
      size.catch.only_TDGDLF[[l]]=apply.Alex.catch.length(Init_F=0.025,InitRec=10,NatMort=NatMort,NatMort_sd=NatMort_sd,
                                                          Steepness=Steepness,Steepness_sd=Steepness_sd,lnRecDev=lnRecDev,
                                                          Lo=Lo,Linf=Linf,vbK=vbK,CVLenAtAge=CVLenAtAge,SDGrowthRec=SDGrowthRec,
                                                          MaxLen=MaxLen,LenInc=LenInc,MaxAge=MaxAge,
                                                          MatL50=MatL50,MatL95=MatL95,PropFemAtBirth=PropFemAtBirth,
                                                          wtlen_a=wtlen_a,wtlen_b=wtlen_b,
                                                          SelAtLength=SelAtLength,SelL50=MatL50, SelL95=MatL95,
                                                          Len_SimYr=Len_SimYr,n_SimLen=n_SimLen,ObsLenComp=ObsLenComp,
                                                          Katch=Katch,nsims=500,UnfishRec=UnfishRec,lnSigmaR=lnSigmaR )
      res$Table.estimates
      
      plot(res$Observed.catch)
      lines(res$Predicted.catch)
      
      for(x in 1:nrow(res$Observed.LenComp))
      {
        plot(res$midpt,res$Observed.LenComp[x,],type='h')
        lines(res$midpt,res$Predicted.LenComp[x,],col='4')
      }
      
    }
  }
}
})

######### remove this
# ObsLenComp=read.csv(handl_OneDrive("Analyses/Population dynamics/Other people's code/Alex dynamic catch & size/ObsLenComp.csv"))
# names(ObsLenComp)=str_remove(colnames(ObsLenComp), "[X]")
# ObsAnnCatch=read.csv(handl_OneDrive("Analyses/Population dynamics/Other people's code/Alex dynamic catch & size/ObsAnnCatch.csv"))
# Katch=data.frame(finyear=1:length(unlist(ObsAnnCatch)), Tonnes=unlist(ObsAnnCatch))
# 
# SelAtLength=read.csv(handl_OneDrive("Analyses/Population dynamics/Other people's code/Alex dynamic catch & size/SelAtLength.csv"))
# SelAtLength=SelAtLength$x
##
#lnSigmaR: assumed recruitment variation,assumed fairly low value for sharks



###########


if(do.Size.based.Catch.curve)
{
  #note: derive F from catch curve and gear selectivity
  #      assume start of year in January (coincides with birth of most species)  
  #      size is TL in mm
  mm.conv=10 # convert total length to mm  
  fn.source("Length_based.catch.curve.R")
 
  
    #24.1 TDGDLF 
  size.catch.curve_TDGDLF=vector('list',N.sp)
  names(size.catch.curve_TDGDLF)=Keep.species
  system.time({for(l in 1: N.sp)  
  {
    # Calculate F
    Outfile='TDGDLF' 
    this.size.comp=paste('Size_composition',c('West.6.5','West.7','Zone1.6.5','Zone1.7','Zone2.6.5','Zone2.7'),sep="_")
    outfile=paste(Outfile,'_histogram',sep='')
    
    #get size composition
    iid=Species.data[[l]][fn.extract.dat(this.size.comp,names(Species.data[[l]]))]
    if(length(iid)>0)
    {
      print(paste("Size-based catch curve with dome-shape selectivity for --",names(Species.data)[l],"---",Outfile))
      
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
        Linf = mm.conv*with(List.sp[[l]],Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)  #total length in mm 
        vbK = List.sp[[l]]$Growth.F$k
        Lo =  mm.conv*with(List.sp[[l]],Lzero*a_FL.to.TL+b_FL.to.TL)   #total length in mm   
        #tzero = 0              
        CVLenAtAge = 0.1      
        
        
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
        NatMort= mean(colMeans(store.species.M[[l]],na.rm=T))
        
        # gillnet selectivity (6.5 and 7 inch combined)
        SelAtLength=Selectivity.at.totalength[[l]]%>%               
                      mutate(TL=TL*mm.conv)%>%
                      filter( TL%in%midpt)%>%
                     pull(Sel.combined)
        #plot(midpt,SelAtLength)
        
        #Execute model functions           
        MeanSizeAtAge = CalcMeanSizeAtAge(Lo,Linf, vbK)
        RecLenDist = CalcSizeDistOfRecruits(GrowthCurveResults=MeanSizeAtAge, CVLenAtAge)
        #plot(midpt,RecLenDist)
        
        LTM = CalcLTM(Linf, vbK, CVLenAtAge, midpt)   
        #image(1:nrow(LTM),1:ncol(LTM),as.matrix(LTM))
        
        #Fit model for each year with length data using parallel processing
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
            
            
            #Calculate uncertainty for parameter estimates
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
        
      }  #end of "if  nrow(N.min)>0 "statement
      
    }  #end of "if length(iid)>0" statement
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
  
    #24.2 NSF
  #note: Not possible. Unknown selectivity for longline & only enough data for sandbar
  

}





#---25. POPULATION DYNAMICS. 'Other' species-------------------------------------------------

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
        if(is.na(SCENARIOS[[sc]]$R.prior[1]))    
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

#---26. RESULTS.'Other' species assessment ------------------------------------------------- 
  #---RESULTS. SPM    MISSING: instead of my ADMB SPM, do JABBA ------ 
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
    fn.word.table(TBL=Tab.par.estim.SPM,Doc.nm="Table 2. SPM estimates")
    rm(HR.o.scens)
    
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
    fn.word.table(TBL=Tab.aSPM,Doc.nm="Table 3. aSPM estimates")
  }
  

#---27. POPULATION DYNAMICS. Indicator species-------------------------------------------------
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




#---28. RESULTS. Indicator species assessment ------------------------------------------------- 

#---29. RESULTS. Final risk ------
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


#---30. RESULTS. Sawfish paper ------------------------------------------------- 
if(do.sawfish)
{
  hNdl.sawfish=handl_OneDrive(paste('Analyses/Population dynamics/Sawfishes/',Year.of.assessment,sep=''))
  if(!dir.exists(hNdl.sawfish))dir.create(hNdl.sawfish)
  
  #1. Plot catches
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
  
  for(i in 1:n.sawfish)  #missing
  {
    #DBSRA
    yrs=Catch_only_sawfish$DBSRA[[i]][['S1']]$output$Years
    Bmsy=apply(Catch_only_sawfish$DBSRA[[i]][['S1']]$output$B.Bmsy,2,median,na.rm=T)
    Fmsy=apply(Catch_only_sawfish$DBSRA[[i]][['S1']]$output$F.Fmsy,2,median,na.rm=T)
    p.DBSRA=kobePlot(f.traj=Fmsy[1:length(yrs)],
                     b.traj=Bmsy[1:length(yrs)],
                     Years=yrs,
                     Titl=paste("DBSRA",names(Catch_only_sawfish$DBSRA)[i],sep='-'))
    rm(yrs,Fmsy,Bmsy)
    
    #CMSY
    if(CMSY.method=="Froese")
    {
      yrs=Catch_only_sawfish$CMSY[[i]][['S1']]$output$ref_ts$year
      Bmsy=Catch_only_sawfish$CMSY[[i]][['S1']]$output$ref_ts$bbmsy
      Fmsy=Catch_only_sawfish$CMSY[[i]][['S1']]$output$ref_ts$ffmsy
      
    }
    if(CMSY.method=="Haddon")
    {
      yrs=Catch_only_sawfish$CMSY[[i]][['S1']]$output$Years
      Bmsy=apply(Catch_only_sawfish$CMSY[[i]][['S1']]$output$B.Bmsy,2,median,na.rm=T)[1:length(yrs)]
      Fmsy=apply(Catch_only_sawfish$CMSY[[i]][['S1']]$output$F.Fmsy,2,median,na.rm=T)[1:length(yrs)]
    }
    
    p.CMSY=kobePlot(f.traj=Fmsy,
                    b.traj=Bmsy,
                    Years=yrs,
                    Titl=paste("CMSY",names(Catch_only_sawfish$DBSRA)[i],sep='-'))
    rm(yrs,Fmsy,Bmsy)
    
    #JABBA
    p.JABBA=with(Catch_only_sawfish$JABBA[[i]][['S1']]$output,
                 {
                   kobePlot(f.traj=timeseries[, , "FFmsy"][,"mu"],
                            b.traj=timeseries[, , "BBmsy"][,"mu"],
                            Years=yr,
                            Titl=paste("JABBA",names(JABBA.sawfish)[i],sep='-'),
                            Probs=data.frame(x=kobe$stock,
                                             y=kobe$harvest))
                 })
    
    figure <- ggarrange(plotlist=list(p.DBSRA+rremove("axis.title"),
                                      p.CMSY+rremove("axis.title"),
                                      p.JABBA+rremove("axis.title")),
                        ncol=1,nrow=3,common.legend = FALSE)
    
    annotate_figure(figure,
                    bottom = text_grob(expression(B/~B[MSY]), size=16),
                    left = text_grob(expression(F/~F[MSY]), rot = 90,size=16))
    
    ggsave(paste(hNdl.sawfish,paste('Kobe_',names(Catch_only_sawfish$DBSRA)[i],'.tiff',sep=''),sep='/'),
           width = 8,height = 14, dpi = 300, compression = "lzw")
    
    
  }

}