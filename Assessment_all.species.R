# ------ Script for running stock assessments for WA sharks---- ###################

#MISSING:
#   Include ALL species in final risk scoring
#   Double check that TwT calculation uses TL and not FL, be consistent
#   Indicator species:
#       Run integrated models for 4 indicator species (remove irrelevant scenario folders).
#       Issues for integrated model:
#           Currently assuming 1 fleet only (TDGDLF). Ok for gummy and whiskery, not OK for 
#               dusky and sandbar (e.g. 'proportional effort for each mesh', 'selectivity')
#           Currently starting 'proportional effort for each mesh' in 1975 but new data goes back to 1940
#       set up SS3 model for indicator species
#       18. Export available data: missing dusky, gummy, sandbar integrated model
#Notes:
#     1. The PSA filters out species from further analyses thru 'Criteria for selecting what species to assess'
#     2. Assumption: If catches have never been >1% carrying capacity, then it's in unexploited 
#                    status so catch series have no information on productivity
#     3. Total (reconstructed) catches are used (commercial, recreational and TDGDLF discards);
#             recons_NT_catch.csv' only includes dusky and sandbar
#     4. A range of assessment methods are used depending on data availability: 
#               . Integrated size- and sex- structured models,
#               . State-space SPM,
#               . Catch-only

#Steps: 
#     1. For each new assessment, update 'Year.of.assessment' and 'Last.yr.ktch' in '1. DEFINE GLOBALS'  
#     2. Define arguments used in each of the shark species/species-complex assessed.
#     3. Bring in updated available data and input parameters
#     4. Determine which species to assess based on PSA
#     5. Run relevant population models according to data availability
#     6. Generate relevant outputs


# ------ Header---- 
rm(list=ls(all=TRUE))
options(dplyr.summarise.inform = FALSE)
options(stringsAsFactors = FALSE) 
options(ggrepel.max.overlaps = Inf)
fn.user=function(x1,x2)paste(x1,Sys.getenv("USERNAME"),x2,sep='/')
if(!exists('handl_OneDrive')) source(fn.user(x1='C:/Users',
                                             x2='OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R'))

# install.packages("devtools")  
# install.packages("Rcpp")      
library(MASS)
library(plotrix)
library(PBSmapping)
library(tidyverse)
library(mvtnorm)  
library(rlist)
library(readxl)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("Biobase")
library(Biobase)
library(numDeriv)
library(spatstat.utils)
library(Hmisc)
library(ggplot2)
library(ggrepel)
#devtools::install_github("haddonm/datalowSA",build_vignettes=TRUE,force=TRUE)
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
library(scales)
library(gridExtra)
library(purrr)
library(yarrr)
library(coda)
library(ggmcmc)
library("Rcpp")
library(tseries)
library(webr)
library(truncnorm)
library(TruncatedDistributions)
library(strex)
library(doSNOW)
library(tictoc)
library(ks)
library(ggwordcloud)

clear.log <- function(x, env = globalenv()) if(exists(x, envir = env))  rm(list = x, envir = env)

  #source stuff and create handy functions
source.hnld=handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/")
fn.source=function(script)source(paste(source.hnld,script,sep=""))
fn.source("fn.fig.R")
fn.source("size transition.R")
#fn.source("Catch_MSY.R")  #superseded by Haddon's package
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")) 
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/send.emails.R"))
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
colfunc <- colorRampPalette(c("red","yellow","springgreen","royalblue"))
fun.find.in.list=function(x,Drop=NULL)   #drop stuff from list
{
  x=x%>%purrr::discard(is.null)
  if(!is.null(Drop))x=x[-match(Drop,names(x))]
  return(x)
}
source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/ggplot.themes.R'))  #my themes
fn.ji=function(a) base::jitter(a,factor=5)  #Add randomness to pin values  
fn.source1=function(script)source(paste(handl_OneDrive("Analyses/Population dynamics/Git_Stock.assessments/"),script,sep=""))
#fn.source1("Create.int.mod.inputs.r")     REMOVE if not used, also remove from repository
fn.source1("auxiliary_functions.r")
fn.source1("SS3_functions.R")
fn.extract.dat=function(STRING,Files) grep(paste(STRING,collapse="|"), Files, value=TRUE)

#---1. DEFINE GLOBALS----- 

Send.email.to="matias.braccini@dpird.wa.gov.au"   #send email when model run finalised
#Send.email.to="braccinimatias@gmail.com.au"  #IT firewall doesn't allow sending to gmail :(

#1. Control which assessment method to implement
Do.Ktch.only=FALSE
Do.StateSpaceSPM=TRUE
Do.integrated=TRUE  #SS3
Do.bespoked=FALSE   #bespoked size-based
do.Andre.model=FALSE
Do.Changes.in.observed.size=FALSE
Do.Changes.in.reported.size=FALSE
do.Size.based.Catch.curve=FALSE  #superseded by dynamic catch-only and size model
do.Dynamic.catch.size.comp=FALSE #not applicable due to poor contrast in available size comp
Do.sim.Test="NO" #do simulation testing of size-based model?

do.F.series=TRUE   #plot these time series
do.B.over.Bmsy.series=TRUE
do.F.over.Fmsy.series=TRUE

Simplified.scenarios=TRUE  #only 1 base case scenario pero catch-only method


#2. Year of assessment and catches
Year.of.assessment=2022
AssessYr=Year.of.assessment
Last.yr.ktch="2021-22"
Last.yr.ktch.numeric=as.numeric(substr(Last.yr.ktch,1,4))
  #future projections
future.models=c('SS')   #c('Catch_only','State.Space.SPM','SS')
years.futures=5  #number of years to project
n.last.catch.yrs=5 #number of recent years used to calculate future catch
catches.futures="constant.last.n.yrs"
#catches.futures='upper.limit.catch.range'  #catch ranges are only available for indicator species
future.color="brown4"

#3.New assessment
New.assessment="NO"
#New.assessment="YES"   #set to 'YES' the first time a new assessment is run


#4. Model run
if(New.assessment=="YES") First.run="YES"  else #create all model input data sets and data presentation for new assessment
  First.run="NO"


#5. Add additional species of interest not selected by PSA given low catch trajectories but needed for specific assessment
additional.sp=NULL  #if no additional species assessment required
#if(Year.of.assessment==2022) additional.sp=c('green sawfish','narrow sawfish')   # 2022 sawfish assessment; 
#                                                                  # dwarf and freshwater sawfish not assessed; recons
                                                                  # catch likely to be incomplete (no TO catch or beach rec fishing)

#6. Define if calculating r & steepness
if(New.assessment=="YES") do.r.prior=TRUE  else 
                          do.r.prior=FALSE
do.steepness=do.r.prior


#7. Define if exporting figures as jpeg or tiff (creation of RAR requires jpeg)
Do.tiff="YES" 
Do.jpeg="NO"


#8. Catch units
KTCH.UNITS="TONNES" 
#KTCH.UNITS="KGS"    
if(KTCH.UNITS=="KGS") unitS=1000
if(KTCH.UNITS=="TONNES") unitS=1


#9. Criteria for selecting what species to assess quantitatively
Min.yrs=5
if(KTCH.UNITS=="KGS") Min.ktch=5000 
if(KTCH.UNITS=="TONNES") Min.ktch=5


#10. CPUEs
Min.cpue.yrs=5 #minimum number of years in abundance index
drop.large.CVs=FALSE  #drop observations with CV larger than MAX.CV or not
survey_not.representative=c("scalloped hammerhead","great hammerhead",
                            "lemon shark","pigeye shark") #Very few individuals (typically <5 per trip)  
                                                          # caught and CVs are huge.
NSF_not.representative=c("scalloped hammerhead","great hammerhead",   #NSF was dropped from assessment due to convergence issues
                          "lemon shark","pigeye shark","tiger shark",
                         "dusky shark" ,"sandbar shark")
tdgdlf_not.representative="smooth hammerhead"       #catch rates are for hammerheads
other_not.representative=c("green sawfish","narrow sawfish") #Pilbara trawl cpue, rare event, not distribution core
drop.daily.cpue='2007&2008'  #drop 2007 & 2008 from TDGDLF daily cpue (consistently higher cpues across all species due to likely effort reporting bias)


#11. Size composition
MN.SZE=0    # initial bin size
#MN.SZE="size.at.birth"
TL.bins.cm=5  # size bin
Min.obs=10  #keep records with at least 10 observations
Min.shts=5  #keep records from at least 5 shots
Min.annual.obs.ktch=50 #Minimum number of annual observations for using size composition data
Min.Nsamp=10   #Minimum number of trips for catch mean weight or size composition
fill.in.zeros=TRUE  #add missing size classes with all 0s
  

#12. Proportion of vessels discarding eagle rays in last 5 years (extracted from catch and effort returns)
prop.disc.ER=.4  


#13. PSA criteria
PSA.min.tons=5
PSA.min.years=Min.yrs
PSA.max.ton=50
Low.risk=2.64  #risk thresholds from Micheli et al 2014
medium.risk=3.18


#14. Assumed PCM for reconstructed discards in TDGLDF
TDGLDF.disc.assumed.PCM="BaseCase" 
#TDGLDF.disc.assumed.PCM="100%" 


#15. Demography
Max.Age.up.Scaler=1.3  #scaler of Max maximum age for species with only one record. 1.3 is mean across those species with a range
Max.r.value=.55 # Max r. .44 for blue shark (Cortes 2016); .51 for Scyliorhinus canicula (Cortes 2002)
Min.r.value=.025


#16. Stock recruitment
Max.h.shark=.8   #mean of h for blue shark (ICCAT 2023 assessment; Cortes 2016, Kai & Fujinami 2018).
Min.h.shark=.3  #He et al 2006, Jason Cope pers comm
Max.SR_sigmaR.shark=0.3   #maximum recruitment variability (blue shark ICCAT 2023 0.29 for North, 0.5 for South, 0.3 bigskate; 0.2 dogfish; 0.18 sandbar)
do.random.h=TRUE  #take a random sample of h and M for SS or use empirical distributions


#17. Reference points
#note: Historically, there was a single reference point (40% unexploited biomass)
Biomass.threshold='Bmsy'  #MSC sets threshold to Bmsy and limit to 0.5 Bmsy (Clinton Syers)
Tar.prop.bmsny=1.2    # Target and Limit proportions of 'Biomass.threshold' 
Lim.prop.bmsy=0.5    #    source: Haddon et al 2014. 'Technical Reviews of Formal Harvest Strategies'.
#Fmsy.emp=function(M) 0.41*M     #Zhou et al 2012 but see Cortes & Brooks 2018
#SPR.thre=0.3   #Punt 2000 Extinction of marine renewable resources: a demographic analysis. 
#SPR.tar=0.4    #   Population Ecology 42, 
Minimum.acceptable.blim=0.2  #Blim should be the max(Minimum.acceptable.blim,Lim.prop.bmsy*Biomass.threshold)


#18. Catch-only Models
do.ensemble.simulations=FALSE
do.OCOM=FALSE
do.Catch.JABBA=FALSE   #redundant, same results as CMSY
catch.only=c('DBSRA','CMSY','SSS')
if(do.OCOM) catch.only=c(catch.only,'OCOM')
if(do.Catch.JABBA) catch.only=c(catch.only,'JABBA')
CMSY.method="Haddon"    # Select which CMSY method to use. Haddon's datalowSA
#CMSY.method="Froese"  # Froese et al 2017 does not converge for dwarf or freshwater
do.parallel.SSS=TRUE   #do SSS in parallel or not (set to FALSE)
do.parallel.SS=TRUE   #do SS in parallel or not (set to FALSE)
SSS.sims=5e2   #Cope 2013 did 1e3 
ensims<-1e4   #DBSRA. 1e4 Dick & MacCall 2011
ensims.CSMY=2e4   #CMSY simulations
ensims.JABBA=3e4  # 3e4 Winkner et al 2019
Proc.Error=1e-3   #Catch-only default process error. Catch only per se yields highly uncertain estimates
Proc.Error.1=1e-2
SSS_criteria.delta.fin.dep=0.01
K_min=2000 #Min value of max K range in tonnes (based on overall catch ranges and K estimates)
r.prob.max=0.9999   #quantile probs for defining r range for CMSY
r.prob.min=1-r.prob.max
#Ensemble.weight='weighted'  #define if doing weighted or unweighted model average
Ensemble.weight='equal'
tweak.BmsyK.Cortes=FALSE #update value to increase COM acceptance rate
tweak.Final.Bio_low=FALSE #update value to increase COM acceptance rate

  #Assessed species by catch methods
Assessed.ktch.only.species='All'                #assess all species with catch only methods (level 1 assessment)
#Assessed.ktch.only.species='Only.ktch.data'   #assess species with only catch data
display.only.catch.only.sp=FALSE               #just display catch and MSY for species with only catch


#19. State space Surplus Production Models
state.space.SPM='JABBA'
do.parallel.JABBA=TRUE #do JABBA in parallel or not (set to FALSE)
use.auxiliary.effort=FALSE  #using effort as auxiliary data yielded very high RMSE and poor fits to effort
Rdist = "lnorm"
KDIST="lnorm"  
PsiDist='beta'
Whiskery.q.periods=2  # split monthly cpue into this number of periods (Simpfendorfer 2000, Braccini et al 2021)
Gummy.q.periods=1     
Obs.Err.JABBA=0.01   #JABBA uses SE2" = CPUE.se^2 + fixed.obsE^2 
increase.CV.JABBA=TRUE   #for consistency with SS3 and because not using fixed.obsE
do.MCMC.diagnostics=do.hindcasting=FALSE
if(First.run=="YES")  do.MCMC.diagnostics=do.hindcasting=TRUE
Proc.Error.cpue=1e-01   # Default process error when fitting cpue to allow enough flexibilty (0.2 showed no differences) 
Proc.Error.cpue2=5e-02   #alternative, 5e-02 process error for JABBA  (Winker et al 2018 School shark) ; Beth Babcock suggested 0.01
k.cv=2                #Carrying capacity CV (Winker et al 2018 School shark)
PEELS=5               #number of years to peel for hindcasting and retrospective   

#20. Integrated age-based model 
Integrated.age.based='SS'
SS3.run='test'  # switch to 'final' when model fitting is finalised to estimate uncertainty (Hessian, etc)
Calculate.ramp.years=FALSE  #switch to TRUE if new year of size composition available
do.Cond.age.len.SS.format=FALSE   #use age-length data to estimate growth
                                  # this is not used as age-length sandbar and dusky is for GN and LL and 
                                  # for all 4 species observations were collected over multiple years
Mean.Size.at.age.species=NULL   #  Mean.Size.at.age.species=c("gummy shark","whiskery shark" )
                                # Not implemented. Wrong SS3 format. May be applicable to gummy and 
                                # whiskery (only for these species length-@-age data collected from gillnet fishery)

Use.SEDAR.M=FALSE   #Set to TRUE if using SEDAR M @ age for dusky and sandbar
SS3_fleet.size.comp.used=c("Size_composition_West","Size_composition_Zone1","Size_composition_Zone2",
                           "Size_composition_NSF.LONGLINE","Size_composition_Survey",
                           "Size_composition_Other")
combine.scallopedHH_NSF_Survey=FALSE

  #SS model run arguments
if(SS3.run=='final') Arg=''
if(SS3.run=='test') Arg= '-nohess'   #no Hessian 
SS3.q.an.sol=TRUE   #calculate q analytically to save up pars
Arg.no.estimation='-maxfn 0 -phase 50 -nohess'  #no estimation
nMCsims=200  #number of Monte Carlo simulations for multivaritenormal
#MCMCsims=1e5; Thin=10; burning=1:(5*length(seq(1,MCMCsims,by=Thin))/100)   #5%  burning
#Arg=paste(' -mcmc',MCMCsims,' -mcsave', 100)  #MCMC
  
  #Assumed error distribution for abundance series
Abundance.error.dist='Lognormal'  #'Lognormal' if stand. cpue in normal space and CVs; 'Normal'

  #Default CV for time series with very small CVs
CV.use='loess'  #Francis 2011
#CV.use='fixed'
default.CV=0.15  #0.15 used by Punt 2009 gummy model; 0.1 Taylor Big skate; loess method 2015 ICATT blue shark 
                  # Taylor dogfish used CV as is so set Q_extraSD to 0.1 
default.Mean.weight.CV=0.2  #bit larger otherwise as it's the only signal for Southern2 selectivity

  #Drop single year size comp
Drop.single.year.size.comp=FALSE

  #Fit diagnostics
do.SS3.diagnostics=FALSE
if(First.run=="YES") do.SS3.diagnostics=TRUE    #very time consuming
Retro_start=0; Retro_end=5 #Last 5 years of observations for retrospective analysis
Number.of.jitters=5       #Number of jitters for Jitter analysis.       MISSING: bump up to 50

resample.h.greynurse=FALSE  #no need to resample h

#21. Bespoke Integrated size-based model 
Plus.gp.size=1.25  #add 25% to max size make sure no accumulation of survivals in last size class
if(Do.bespoked)
{
  if(First.run=="YES") run.all.scenarios="YES" #What scenarios to run?
  if(First.run=="NO") run.all.scenarios="NO"  
  if(First.run=="NO") run.future="YES" #Future projection scenarios
  if(First.run=="YES") run.future="NO"
  Show.yrs="DATA"   #What years to show in model outputs
  #Show.yrs="FUTURE"  #show data plus future projections
  Present.in.log="NO"  #present cpue in log space or normal space
  Do.zero.Ktch="NO" #do zero catch projections of size-based model?
  mcmc.do="YES"
  mcmc.n=1e6
  mcmc.thin=10
  MIN.DAYS.LARGE=30  #minimum number of days at liberty for movement
  Sim.trans.Mat="NO" #show simulated size transition matrix
  add.conv.tag="YES" #define if using conventional tagging and effort in model
  Move.mode="Individual-based"  #select type of movement modelling approach
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
  
}

#22. Weight of Evidence
LoE.Weights=c(psa=.25,sptemp=.25,efman=.25,COM=.5,JABBA=.75,integrated=1)  
RiskColors=c('Negligible'="cornflowerblue",
             'Low'="olivedrab3",
             'Medium'="yellow",
             'High'="orange",
             'Severe'="red2")
Choose.probability="B.over.Bmsy" #"Depletion"  #use B/Bmys or B/K probabilities
Like.ranges=list(L1=c(0,0.0499999),
                 L2=c(0.05,0.2),
                 L3=c(0.20001,0.5),
                 L4=c(0.50001,1))

#23. Average ratio L50:L95 for species with no L95 estimates
average.prop.L95_L50=mean(c(135.4/154.5,225/262,210/240,113/138,175/198,281/328,113/139,125/136,113/138))

#24. Define if using effort 
add.effort="NO"    
What.Effort="km.gn.hours"  #What effort to display?
#What.Effort="km.gn.days" 


#---2. Catch and effort data -----   

Dat.repository=handl_OneDrive('Analyses/Data_outs/')  #locations where all data are stored
Dat.repository2=handl_OneDrive('Data/Population dynamics/Data inputs for models')

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
  #NSF
Effort.monthly.north=fn.in(NM='Annual.total.eff_NSF.csv') 



#2.3. Import Total Catch
Mode <- function(x)
{
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
Get.ktch=fn.import.catch.data(KTCH.UNITS)    
KtCh=Get.ktch$Total
KtCh.zone=Get.ktch$Zone

ktch.combined=KtCh%>%
        group_by(SPECIES,Name,finyear)%>%
        summarise(Tonnes=sum(LIVEWT.c,na.rm=T))

clear.log('fn.import.catch.data')

#export catch by species for Rays Report Card
doRayRepCard=FALSE
if(doRayRepCard)
{
  hendl=handl_OneDrive('Workshops/2022_Rays Report Card')
  #Sharks
  dummy.sharks=Get.ktch$Total.method%>%
          filter(SPECIES%in%Shark.species)%>%
          group_by(Name,FishCubeCode,finyear,Gear)%>%
          summarise(Tonnes=sum(LIVEWT.c))
  uni.sp=unique(dummy.sharks$Name)
  for(u in 1:length(uni.sp))
  {
    dummy1=dummy.sharks%>%
            filter(Name==uni.sp[u])  
    dummy1%>%
    ggplot(aes(finyear,Tonnes))+
      geom_point(aes(colour = FishCubeCode),size = 3)+
      ylab("Catch (tonnes)")+xlab("Finacial year")+
      facet_wrap( ~ Gear, scales = "free_y")+ expand_limits(y = 0)+
      theme_PA(strx.siz=17,leg.siz=18,axs.t.siz=15,axs.T.siz=24)+
      theme(legend.position="top",
            legend.title = element_blank(),
            legend.key=element_blank(),
            title=element_text(size=12),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
      guides(colour = guide_legend(override.aes = list(size=5)))+
      xlab("Financial year")+ggtitle(capitalize(uni.sp[u]))
    ggsave(paste(hendl,'/Sharks/Plot_',uni.sp[u],'.tiff',sep=''), width = 12,height = 12, dpi = 300, compression = "lzw")
    
    write.csv(dummy1,paste(hendl,'/Sharks/Stats_',uni.sp[u],'.csv',sep=''),row.names = F)
    
    rm(dummy1)
  }
  rm(dummy.sharks)
  
  #Rays
  dummy.rays=Get.ktch$Total.method%>%
    filter(!SPECIES%in%Shark.species)%>%
    group_by(Name,FishCubeCode,finyear,Gear)%>%
    summarise(Tonnes=sum(LIVEWT.c))
  uni.sp=unique(dummy.rays$Name)
  for(u in 1:length(uni.sp))
  {
    dummy1=dummy.rays%>%
      filter(Name==uni.sp[u])  
    dummy1%>%
      ggplot(aes(finyear,Tonnes))+
      geom_point(aes(colour = FishCubeCode),size = 3)+
      ylab("Catch (tonnes)")+xlab("Finacial year")+
      facet_wrap( ~ Gear, scales = "free_y")+ expand_limits(y = 0)+
      theme_PA(strx.siz=17,leg.siz=18,axs.t.siz=15,axs.T.siz=24)+
      theme(legend.position="top",
            legend.title = element_blank(),
            legend.key=element_blank(),
            title=element_text(size=12),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
      guides(colour = guide_legend(override.aes = list(size=5)))+
      xlab("Financial year")+ggtitle(capitalize(uni.sp[u]))
    ggsave(paste(hendl,'/Rays/Plot_',uni.sp[u],'.tiff',sep=''), width = 12,height = 12, dpi = 300, compression = "lzw")
    
    write.csv(dummy1,paste(hendl,'/Rays/Stats_',uni.sp[u],'.csv',sep=''),row.names = F)
    
    rm(dummy1)
  }
  rm(dummy.rays)
  
}

#export catch for EPBC-listed species 
doEPBC=FALSE
if(doEPBC)
{
  EPBC_sharks_rays=tolower(c('Harrisson’s Dogfish','Southern Dogfish',
                             "whale shark","grey nurse shark",'Basking Shark',
                             'White Shark','Shortfin Mako','Longfin Mako',
                             'Porbeagle','School Shark','Silky Shark','Oceanic Whitetip Shark',
                             'Northern River Shark','Speartooth Shark',
                             "green sawfish","narrow sawfish","freshwater sawfish","dwarf sawfish",
                             'Scalloped Hammerhead','Maugean Skate','Reef Manta Ray','Giant Manta Ray',
                             'Long-horned Pygmy Devilray','Giant Devilray','Bentfin Devilray'))
  hendl=handl_OneDrive('Workshops/2022_NESP Hub Shark and Ray')
  EPBC_sharks_rays_ktch=Get.ktch$Total.method%>%   #note: missing white shark TEPS post 2006
    filter(Name%in%EPBC_sharks_rays)%>%
    group_by(Name,Data.set,finyear,Gear)%>%
    summarise(Tonnes=sum(LIVEWT.c))
  uni.sp=unique(EPBC_sharks_rays_ktch$Name)
  for(u in 1:length(uni.sp))
  {
    dummy1=EPBC_sharks_rays_ktch%>%
      filter(Name==uni.sp[u])  
    dummy1%>%
      ggplot(aes(finyear,Tonnes))+
      geom_point(aes(colour = Data.set),size = 3)+
      ylab("Catch (tonnes)")+xlab("Finacial year")+
      facet_wrap( ~ Gear, scales = "free_y")+ expand_limits(y = 0)+
      theme_PA(strx.siz=17,leg.siz=18,axs.t.siz=15,axs.T.siz=24)+
      theme(legend.position="top",
            legend.title = element_blank(),
            legend.key=element_blank(),
            title=element_text(size=12),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
      guides(colour = guide_legend(override.aes = list(size=5)))+
      xlab("Financial year")+ggtitle(capitalize(uni.sp[u]))
    ggsave(paste(hendl,'/Plot_',uni.sp[u],'.tiff',sep=''), width = 12,height = 12, dpi = 300, compression = "lzw")
    
    write.csv(dummy1,paste(hendl,'/Stats_',uni.sp[u],'.csv',sep=''),row.names = F)
    
    rm(dummy1)
  }
  rm(EPBC_sharks_rays_ktch)
}

#plot proportional catch by fishery
Shark.Fisheries=c('Joint Authority Southern Demersal Gillnet and Demersal Longline Managed Fishery'='JASDGDL',
                  'West Coast Demersal Gillnet and Demersal Longline (Interim) Managed Fishery'='WCDGDL',
                  'FBL condition 70 Power operated net hauler'='C070',
                  'Western Australian North Coast Shark Fishery'='WANCS',
                  'Open Access in the West Coast Bioregion'='OAWC',
                  'Open Access in the North Coast and Gascoyne Coast Bioregions'='OANCGC',
                  'Joint Authority Northern Shark Fishery'='JANS',
                  'Discards_TDGDLF'='Discards_TDGDLF',
                  'TEP_greynurse'='TEP_greynurse',
                  'TEP_dusky'='TEP_dusky')

Non.WA.fisheries=c('GAB.trawl','Indonesia','NSW fisheries','NT','SA MSF','Taiwan')

do.pie.donut=FALSE
if(do.pie.donut)
{
  fn.fig(handl_OneDrive("Analyses/Population dynamics/Proportional catch_WA"),2400,2400)
 KtCh.method%>%
      filter(!Data.set%in%c('Historic',Non.WA.fisheries))%>%
      mutate(Shark.Fishery=ifelse(FishCubeCode%in%Shark.Fisheries,'Shark fisheries',
                                  'Other fisheries'))%>%
      group_by(Shark.Fishery,Gear)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c))%>%
      PieDonut(aes(Shark.Fishery,Gear,count=LIVEWT.c),title='Proportional catch',
               explode = 2,r0=.2,r1=.8,maxx=1.45, titlesize=8,
               pieLabelSize=6.5,donutLabelSize=5.5,showPieName=FALSE)
  
  # ddd=KtCh.method%>%
  #   filter(!Data.set%in%c('Historic',Non.WA.fisheries))%>%
  #   mutate(Shark.Fishery=ifelse(FishCubeCode%in%Shark.Fisheries,'Shark fisheries',
  #                               'Other fisheries'),
  #          FinYr=ifelse(finyear<2006,'Pre landing ban','Post landing ban'))
  # p2=ddd%>%
  #   filter(FinYr=='Pre landing ban')%>%
  #   group_by(Shark.Fishery,Gear)%>%
  #   summarise(LIVEWT.c=sum(LIVEWT.c))%>%
  #   PieDonut(aes(Shark.Fishery,Gear,count=LIVEWT.c),title='Pre landing ban',
  #            explode = 2,r0=.2,r1=.8,maxx=1.45, titlesize=8,
  #            pieLabelSize=6.5,donutLabelSize=5.5,showPieName=FALSE)
  # 
  # p3=ddd%>%
  #   filter(FinYr=='Post landing ban')%>%
  #   group_by(Shark.Fishery,Gear)%>%
  #   summarise(LIVEWT.c=sum(LIVEWT.c))%>%
  #   PieDonut(aes(Shark.Fishery,Gear,count=LIVEWT.c),title='Post landing ban',
  #            explode = 2,r0=.2,r1=.8,maxx=1.45, titlesize=8,
  #            pieLabelSize=6.5,donutLabelSize=5.5,showPieName=FALSE)
  
  dev.off()
  
}

#export catch for NT assessments   
output.NT.ktch=FALSE
if(output.NT.ktch)
{
  Data.for.NT=KtCh%>%
    filter(SPECIES%in%c(18013,18014,18039))%>%
    group_by(SPECIES,finyear,Name,FishCubeCode)%>%
    summarise(Tons=sum(LIVEWT.c,na.rm=T))  #KTCH.UNITS are in tonnes
  write.csv(Data.for.NT,
            handl_OneDrive('Analyses/Catch and effort/Data_Resquests/Northern_Territory/Blacktips.and.spot.tail.reconstructed.catches.csv'),row.names = F)
}

#Composition of shark catch in NT (sorted by harvest)
NT.shark.compo=read.csv(handl_OneDrive('Data/Catch and Effort/NT_shark_catch_composition.csv'))


#2.4. Remove species assessed elsewhere
KtCh=KtCh%>%filter(!Name%in%assessed.elsewhere)
KtCh.zone=KtCh.zone%>%filter(!Name%in%assessed.elsewhere)
ktch.combined=ktch.combined%>%filter(!Name%in%assessed.elsewhere)


#2.5 Data used for analysis of changes in reported mean weight
Wei.range=read.csv(handl_OneDrive("Data/Length_Weights/Data.Ranges.csv"),stringsAsFactors = F)
Wei.range.names=read.csv(handl_OneDrive("Data/Length_Weights/Species.names.csv"),stringsAsFactors = F)
Wei.range=merge(Wei.range,Wei.range.names,by="Sname",all.x=T)

Logbook=read.csv(handl_OneDrive("Analyses/Catch and effort/Logbook.data.mean.weight.csv"))


#---3. Life history and selectivity data ------------------------------------------------------  
LH.data=read.csv(handl_OneDrive('Data/Life history parameters/Life_History.csv'))
SS_selectivity_init_pars=read.csv(handl_OneDrive('Analyses/Population dynamics/SS3.selectivity_pars.csv'))

#---4. PSA (which species to assess quantitatively) -----------------------------------------------  
#note: run a PSA aggregating the susceptibilities of multiple fleets (Micheli et al 2014)
#      catsharks, sawsharks, thresher sharks, wobbegongs, angel sharks, stingrays, banjo rays, 
#        wedgefishes, and rays and skates were aggregated and assessed as a group due to unreliable 
#        species reporting. hence the difference between species in the catch and in PSA assessment

#get catches of all species
KtCh.method=Get.ktch$Total.method%>%filter(!Name%in%assessed.elsewhere)

#biological and fishery attributes for PSA
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
            filter(!is.na(Name))%>%left_join(All.species.names%>%dplyr::select(SPECIES,Scien.nm),by='Scien.nm')%>%
            relocate(c(SPECIES))%>%
            arrange(SPECIES),
          paste(Exprt,'Table S1_All.species.caught.by.fishery.csv',sep='/'),row.names = F)

#Export species word cloud #ACA
Word.cloud=Get.ktch$Table1%>%
  data.frame%>%
  mutate(Tot=rowSums(across(where(is.numeric))))%>%
  dplyr::select(Name,Scien.nm,Tot)%>%
  left_join(All.species.names%>%
              distinct(Scien.nm,.keep_all=T),
            by='Scien.nm')%>%
  filter(!is.na(Scien.nm))%>%
  mutate(CITES=ifelse(CITES%in%c("") |is.na(CITES),NA,"CITES"),
         CMS=ifelse(CMS%in%c("") |is.na(CMS),NA,"CMS"),
         EPBC=ifelse(EPBC%in%c("") |is.na(EPBC),NA,"EPBC"),
         SAFS=ifelse(SAFS%in%c("") |is.na(SAFS),NA,"SAFS"),
         WA_BC_Act=ifelse(WA_BC_Act%in%c("") |is.na(WA_BC_Act),NA,"WA BC Act"),
         Indicators=ifelse(Name%in%capitalize(names(Indicator.species)),'Indicators','Non.indicators'))%>%
  unite(col='Protection',CITES,CMS,EPBC,WA_BC_Act,sep=',',na.rm = TRUE)%>%
  mutate(Listed=ifelse(!Protection=="",'Listed (CITES, CMS, EPBC, WA BC Act)',"Not listed"))

p=Word.cloud%>%
  ggplot(aes(label = Name, size = Tot,color=Indicators))+
  geom_text_wordcloud() +
  scale_size_area(max_size = 20) +
  theme_minimal()+
  scale_color_manual(values=c(Non.indicators = "black", Indicators = "steelblue"))
ggsave(paste(Exprt,'Catch_word_cloud.tiff',sep='/'), width = 10,height = 10, dpi = 300, compression = "lzw")

p1=Word.cloud%>%
  filter(Listed=='Not listed')%>%
  mutate(dummy='dummy')%>%
  ggplot(aes(label = Name, size = Tot,color=dummy))+
  geom_text_wordcloud(show.legend = FALSE, family="Purisa") +
  scale_size_area(max_size = 20) +
  theme_minimal()+
  guides(size = "none")+
  ggtitle("Not listed")+
  theme(plot.title = element_text(size=28,hjust=0.5,vjust = -8))+
  scale_color_manual(values=c(dummy = "grey50"))

p2=Word.cloud%>%
  filter(!Listed=='Not listed')%>%
  ggplot(aes(label = Name, size = Tot,color=Protection))+
  geom_text_wordcloud(show.legend = FALSE, family="Purisa") +
  scale_size_area(max_size = 20) +
  theme_minimal()+
  guides(size = "none")+
  ggtitle("Listed (CITES, CMS, EPBC, WA BC Act)")+
  theme(plot.title = element_text(size=28,hjust=0.5,vjust = -8))


ggarrange(p1,p2,nrow=2)+theme(plot.margin = unit(c(-1,-1,-1,-1), 'lines'))
ggsave(paste(Exprt,'Catch_word_cloud_protection.tiff',sep='/'),
       width = 250,height = 250,units = "mm", dpi = 300, compression = "lzw")

clear.log('Get.ktch')
clear.log('WAislands')  
clear.log('WAcoast')  

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
PSA.out=PSA.fn(d=PSA.list,line.sep=.35,size.low=2.1,size.med=2.15,size.hig=2.5,W=10,H=10)

Keep.species=tolower(as.character(PSA.out%>%filter(Vulnerability=="High")%>%pull(Species)))
Keep.species=sort(c(Keep.species,names(Indicator.species)))
Drop.species=UniSp[which(!UniSp%in%Keep.species)]

RAR.species=Keep.species 
if(!is.null(additional.sp)) Keep.species=c(Keep.species,additional.sp)
Keep.species=sort(Keep.species)
N.sp=length(Keep.species)
Other.species=RAR.species[-match(names(Indicator.species),RAR.species)]

clear.log('PSA.fn')
clear.log('Agg')
clear.log('Agg.PSA')


#---5. Source demography and steepness functions -------------------------------------------------------
if(do.r.prior)  fn.source("Leslie.matrix.R") 
if(do.steepness) fn.source("Steepness.R")

#---6. Import species-specific data -----
#note: this brings in any available data (cpue, abundance, selectivity, size composition, tagging, etc)
Species.data=vector('list',length=N.sp)
names(Species.data)=Keep.species
for(s in 1:N.sp) 
{
  print(paste('Reading in data for -----',Keep.species[s]))
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
      getfiles=sapply(files, fread, data.table=FALSE)
      if(is.matrix(getfiles))
      {
        getfiles=vector('list',length(file.names))
        for(g in 1:length(getfiles)) getfiles[[g]]=read.csv(files[g])
      }
    }else
    {
      getfiles=list(read.csv(files))
    }
    names(getfiles)=file.names
    Species.data[[s]]=getfiles
  }
  rm(files)
}

#add size composition C. brachyurus MSF (South Australia)
Species.data$`copper shark`$Size_composition_Other=read.csv(handl_OneDrive('Data/Population dynamics/SA_Cbrachyurus length data.csv'))%>%
                  mutate(DATE=as.Date(DATE,"%d/%m/%Y"),
                         Month=month(DATE),
                         year=year(DATE),
                         FINYEAR=ifelse(Month>6,paste(year,substr(year+1,3,4),sep='-'),
                                        paste(year-1,substr(year,3,4),sep='-')),
                         FL=with(LH.data%>%filter(SPECIES==18001),((TL..mm./10)-b_FL.to.TL )/a_FL.to.TL),
                         SEX=ifelse(SEX=="Male ","Male",SEX),
                         SEX=ifelse(SEX=='Female','F',ifelse(SEX=='Male','M',NA)))%>%
                         filter(SAMPLE.SOURCE=='Commercial catch sampling')%>%
                         filter(!is.na(SEX))%>%
                  dplyr::select(Month,FINYEAR,year,FL,SEX)

Species.data$`copper shark`$Size_composition_Other_Observations=data.frame(
                  FINYEAR=sort(unique(Species.data$`copper shark`$Size_composition_Other$FINYEAR)),
                  Method='LL',
                  zone="SA",
                  SPECIES=18001,
                  N.shots=30,
                  N.observations=NA)

#add size composition Spotted wobbegong and Western wobbegong to Wobbegongs
other.wobbies=c('Spotted wobbegong','Western wobbegong')
Wobbies.data=vector('list',length(other.wobbies))
names(Wobbies.data)=other.wobbies
  #1. bring in data
for(s in 1:length(other.wobbies)) 
{
  print(paste('Reading in data for -----',other.wobbies[s]))
  this.one=other.wobbies[s]
  files <- list.files(path = paste(Dat.repository,this.one,sep=""), pattern = "*.csv", full.names = T)
  file.names <- list.files(path = paste(Dat.repository,this.one,sep=""), pattern = "*.csv", full.names = F)
  removE <- c(".csv", paste(this.one,"_",sep=''), this.one)
  file.names <- gsub("^\\.","",str_remove_all(file.names, paste(removE, collapse = "|")))
  if(length(files)>0)
  {
    if(length(files)>1)
    {
      getfiles=sapply(files, fread, data.table=FALSE)
      if(is.matrix(getfiles))
      {
        getfiles=vector('list',length(file.names))
        for(g in 1:length(getfiles)) getfiles[[g]]=read.csv(files[g])
      }
    }else
    {
      getfiles=list(read.csv(files))
    }
    names(getfiles)=file.names
    Wobbies.data[[s]]=getfiles
  }
  rm(files)
}  
  #2. add to Wobbegongs
for(s in 1:length(other.wobbies)) 
{
  Species.data$wobbegongs$Size_composition_Observations=rbind(Species.data$wobbegongs$Size_composition_Observations,
                                                              Wobbies.data[[s]]$Size_composition_Observations)
  Species.data$wobbegongs$Size_composition_West.6.5.inch.raw=rbind(Species.data$wobbegongs$Size_composition_West.6.5.inch.raw,
                                                                   Wobbies.data[[s]]$Size_composition_West.6.5.inch.raw)
  Species.data$wobbegongs$Size_composition_West.7.inch.raw=rbind(Species.data$wobbegongs$Size_composition_West.7.inch.raw,
                                                                 Wobbies.data[[s]]$Size_composition_West.7.inch.raw)
  Species.data$wobbegongs$Size_composition_Zone1.6.5.inch.raw=rbind(Species.data$wobbegongs$Size_composition_Zone1.6.5.inch.raw,
                                                                    Wobbies.data[[s]]$Size_composition_Zone1.6.5.inch.raw)
  Species.data$wobbegongs$Size_composition_Zone1.7.inch.raw=rbind(Species.data$wobbegongs$Size_composition_Zone1.7.inch.raw,
                                                                  Wobbies.data[[s]]$Size_composition_Zone1.7.inch.raw)
  Species.data$wobbegongs$Size_composition_Zone2.6.5.inch.raw=rbind(Species.data$wobbegongs$Size_composition_Zone2.6.5.inch.raw,
                                                                    Wobbies.data[[s]]$Size_composition_Zone2.6.5.inch.raw)
  Species.data$wobbegongs$Size_composition_Zone2.7.inch.raw=rbind(Species.data$wobbegongs$Size_composition_Zone2.7.inch.raw,
                                                                  Wobbies.data[[s]]$Size_composition_Zone2.7.inch.raw)
}
Species.data$wobbegongs$Size_composition_Observations=Species.data$wobbegongs$Size_composition_Observations%>%
  mutate(SPECIES=13000)%>%
  group_by(FINYEAR,Method,zone,SPECIES)%>%
  summarise(N.shots=mean(N.shots),
            N.observations=sum(N.observations))%>%
  ungroup()

#remove Pilbara trawl as it's not used at all
for(i in 1:N.sp) 
{
  if(any(grepl('Size_composition_Pilbara_Trawl',names(Species.data[[i]]))))
  {
    Species.data[[i]]=Species.data[[i]][-match('Size_composition_Pilbara_Trawl',names(Species.data[[i]]))]
  }
}

#remove NA sex in length composition data
for(i in 1:N.sp) 
{
  if(any(grepl('Size_composition',names(Species.data[[i]]))))
  {
    d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[i]]))]
    if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
    if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
    if(length(d.list)>0)
    {
      for(s in 1:length(d.list))
      {
        zero.length=d.list[[s]]%>%
          filter(!is.na(SEX))
        if(nrow(zero.length)==0)
        {
          Species.data[[i]]=Species.data[[i]][-match(names(d.list)[s],names(Species.data[[i]]))]
        }
      }
    }
  }
}

#remove NA Ages in length-age data
for(s in 1:N.sp) 
{
  if('age_length'%in%names(Species.data[[s]]))
  {
    Species.data[[s]]$age_length=Species.data[[s]]$age_length%>%
      filter(!is.na(Age))
    
    if(names(Species.data)[s]=="sandbar shark") #Remove dodgy record for sandbar
    {
      Species.data[[s]]$age_length=Species.data[[s]]$age_length%>%
        mutate(Keep=ifelse(Age==0 & FL>75,'No','Yes'))%>%
        filter(Keep=='Yes')%>%
        dplyr::select(-Keep)
    }
  }
  
}

#remove nonsense size comp from milk shark
Species.data$`milk shark`$Size_composition_NSF.LONGLINE=Species.data$`milk shark`$Size_composition_NSF.LONGLINE%>%
                                                          filter(FL<=80)
Species.data$`milk shark`$Size_composition_West.7.inch.raw=Species.data$`milk shark`$Size_composition_West.7.inch.raw%>%
  filter(FL<=80)
Species.data$`milk shark`$Size_composition_Survey=Species.data$`milk shark`$Size_composition_Survey%>%
  filter(FL<=80)

#remove nonsense size comp from sandbar shark
Species.data$`sandbar shark`$Size_composition_NSF.LONGLINE=Species.data$`sandbar shark`$Size_composition_NSF.LONGLINE%>%
  filter(FL<=200)
Species.data$`sandbar shark`$Size_composition_Survey=Species.data$`sandbar shark`$Size_composition_Survey%>%
  filter(FL<=200)

#remove Observer cpue data (double dipping with standardised catch rates)
for(s in 1:N.sp)
{
  iid=grep('CPUE_Observer_TDGDLF',names(Species.data[[s]]))
  if(length(iid)>0) Species.data[[s]]=Species.data[[s]][-iid]
}


#---7. Create list of life history parameter inputs -----
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
  fst.yr=min(KtCh$FINYEAR)
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
  List.sp[[l]]=list.append(List.sp[[l]],
                           pup.sx.ratio=0.5,
                           Growth.F=data.frame(k=LH$K,FL_inf=LH$FL_inf,to=LH$to,k.sd=LH$k.sd,FL_inf.sd=LH$FL_inf.sd),
                           Growth.M=data.frame(k=LH$male_K,FL_inf=LH$male_FL_inf),
                           k.Linf.cor=-0.99,    #assumed correlation between growth parameters
                           Max.age.F=c(LH$Max_Age,LH$Max_Age_max),
                           Age.50.mat=c(LH$Age_50_Mat_min,LH$Age_50_Mat_max),
                           TL.50.mat=LH$TL.50.mat,
                           TL.95.mat=LH$TL.95.mat,
                           Fecundity=c(LH$Fecu_min,LH$Fecu_max),
                           Fecu_a=LH$Fecu_a,
                           Fecu_b=LH$Fecu_b,
                           Fecu_type_SS=LH$Fecu_type_SS,
                           Breed.cycle=c(LH$Cycle,LH$Cycle_max),
                           TEMP=LH$Temperature,
                           BwT=LH$b_w8t,
                           AwT=LH$a_w8t,
                           BwT.M=LH$male_b_w8t,
                           AwT.M=LH$male_a_w8t,
                           TLmax=LH$Max.TL,
                           Lzero=LH$LF_o,
                           NsimSS=5e3,                        #demography
                           r.prior="USER",                    #demography
                           r.prior2=NA,                       #demography uniform
                           a_FL.to.TL=LH$a_FL.to.TL,          # FL to TL
                           b_FL.to.TL=LH$b_FL.to.TL
  )              
  
}

  #To avoid inconsistencies, set Max Age to max value 
for(l in 1:N.sp)
{
  tmax=with(List.sp[[l]],((1/Growth.F$k)*log((Growth.F$FL_inf-Lzero)/( (1-0.99)*Growth.F$FL_inf)))) #theoretical lifespan (Cortes & Taylor 2023)
  Max_Age_max=List.sp[[l]]$Max.age.F[2]
  if(is.na(Max_Age_max)) List.sp[[l]]$Max.age.F[2]=round(List.sp[[l]]$Max.age.F[1]*Max.Age.up.Scaler)  
  List.sp[[l]]$tmax=tmax
  #List.sp[[l]]$Max.age.F[2]=max(tmax,List.sp[[l]]$Max.age.F[2])
  List.sp[[l]]$Max.age.F[1]=List.sp[[l]]$Max.age.F[2]
}
  
  #If no 'to' information, calculate by fitting 3 par vonB to 2 par vonB
L=function(to) LINF*(1-exp(-K*(age-to)))
sumsq <- function( x, y) {sum((x-y)^2)}
ObjFunc=function(Init.to)
{
  Observed=vonB2
  Predicted=L(to=Init.to)  
  return(sumsq(x=Observed,y=Predicted))
}
Init.to=-2
for(l in 1:N.sp)
{
  if(is.na(List.sp[[l]]$Growth.F$to)) 
  {
    age=0:List.sp[[l]]$Max.age.F[2]
    LINF=List.sp[[l]]$Growth.F$FL_inf
    K=List.sp[[l]]$Growth.F$k
    vonB2=LINF-(LINF-List.sp[[l]]$Lzero)*exp(-K*age)
    nlmb <- nlminb(Init.to, ObjFunc, gradient = NULL, hessian = TRUE)
    List.sp[[l]]$Growth.F$to=nlmb$par
    rm(age,K,LINF)
  }
}
rm(sumsq,ObjFunc,L)

  #Export table of life history parameters for use in RAR
Rar.path=paste(handl_OneDrive('Reports/RARs'), AssessYr,'Figures_and_Tables',sep="/")
if(!dir.exists(Rar.path))dir.create(Rar.path)
if(First.run=="YES")
{
  #setwd(Rar.path)
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
    arrange(Species)%>%
    mutate(Max.Age=ifelse(is.na(Max_Age_max),Max_Age,Max_Age_max),
           b_w8t=round(b_w8t,2),
           TL.50.mat=round(TL.50.mat),
           TL.95.mat=round(TL.95.mat),
           Age.50.Mat=paste('uniform(',Age_50_Mat_min,',',Age_50_Mat_max,')',sep=''),
           Mean.Fec=round(mean(c(Fecu_min,Fecu_max))),
           Litter.size=paste('triangular(',Fecu_min,',',Mean.Fec,',',Fecu_max,')',sep=''),
           Rep.cycle=paste('uniform(',Cycle,',',Cycle_max,')',sep=''))%>%
    dplyr::select(Species,K,FL_inf,Max.Age,Age.50.Mat,TL.50.mat,TL.95.mat,
                  Litter.size,Rep.cycle,a_w8t,b_w8t)
  names(TabL)[match(c("K","FL_inf","Max.Age","Age.50.Mat","TL.50.mat","TL.95.mat","Rep.cycle"),names(TabL))]=
    c("K (year-1)","FL_inf (cm)","Max.Age (year)","Age.50.Mat (year)","TL.50.mat (cm)","TL.95.mat (cm)","Rep.cycle (year)")
 # fn.word.table(TBL=TabL,Doc.nm="Table 1. Life history pars")
  write.csv(TabL,paste(Rar.path,"Table 1. Life history pars.csv",sep='/'),row.names = F)
}

  #Export table of species assessed for use in RAR
if(First.run=="YES")
{
  write.csv(data.frame(Species=capitalize(RAR.species)),
            paste(Rar.path,"Assessed_species.csv",sep='/'),row.names = F)
}

#refit indicator species growth curve  
refit.indicators.growth=FALSE
if(refit.indicators.growth) fn.source('Fit.growth.curves.R')

#Get average weight
do.this=FALSE
if(do.this)
{
  Average.weight=vector('list',N.sp)
  Keep.species
  for(p in 1:N.sp)
  {
    attach(List.sp[[p]])
    TTLL=round(Lzero*a_FL.to.TL+b_FL.to.TL):TLmax
    Kg=round(AwT*TTLL^BwT,2)
    detach(List.sp[[p]])
    Average.weight[[p]]=data.frame(Species=names(List.sp)[p],Kg=Kg,TL=TTLL)
  }
  Average.weight=do.call(rbind,Average.weight)
  Average.weight%>%
    ggplot(aes(TL,Kg))+
    geom_line()+
    facet_wrap(~Species,scales='free')
  ggsave(handl_OneDrive("Analyses/Population dynamics/Species_weights.tiff"),
         width = 10,height = 10, dpi = 300, compression = "lzw")
  Average.weight=Average.weight%>%
    group_by(Species)%>%
    summarise(AVg.weight=quantile(Kg,probs=.4))%>%
    data.frame
}

#---8. Create list of species groups for RAR outputs ----- 
Lista.sp.outputs=list(Other.species,names(Indicator.species))
names(Lista.sp.outputs)=c('Other.sp','Indicator.sp')
if(!is.null(additional.sp))
{
  Lista.sp.outputs[[3]]=additional.sp
  names(Lista.sp.outputs)[3]='additional.sp'
}


#---9. Calculate r prior -----  
  #calculate prior
store.species.r_M.min=vector('list',N.sp)
names(store.species.r_M.min)=Keep.species
store.species.M_M.min=store.species.G_M.min=store.species.r_M.min
store.species.r_M.mean=store.species.M_M.mean=store.species.G_M.mean=store.species.r_M.min
if(do.r.prior)  #0.008 sec per iteration per species
{
  set.seed(1234)
  store.species.r=vector('list',N.sp)
  names(store.species.r)=Keep.species
  tic()
  for(l in 1:N.sp)  #no parallel processing implementation as it stuffs up Mmin and Mmean
  {
    print(paste('r calculation -------------',Keep.species[l]))  
    #if no sd, replace with mean from those species with sd
    if(is.na(List.sp[[l]]$Growth.F$FL_inf.sd)) List.sp[[l]]$Growth.F$FL_inf.sd=0.038*List.sp[[l]]$Growth.F$FL_inf  
    if(is.na(List.sp[[l]]$Growth.F$k.sd)) List.sp[[l]]$Growth.F$k.sd=0.088*List.sp[[l]]$Growth.F$k     
    
    RESAMP="YES"
    
    linear.fec="NO"
    #if(names(List.sp)[l]%in%c("grey nurse shark","sandbar shark")) linear.fec="NO"
    
    #Get r prior
    M.averaging<<-"min" #'min' yields rmax, Cortes pers com
    GET.all.Ms=TRUE
    r.prior.dist_M.min=with(List.sp[[l]],fun.rprior.dist(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,
                                                         K.sd=Growth.F$k.sd,LINF.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL,k.Linf.cor,
                                                         Amax=Max.age.F,
                                                         MAT=unlist(Age.50.mat),FecunditY=Fecundity,Cycle=Breed.cycle,
                                                         BWT=BwT,AWT=AwT,LO=Lzero*a_FL.to.TL+b_FL.to.TL)) #size vars as TL
    
    M.averaging<<-"mean"
    GET.all.Ms=FALSE
    mean.donotwork=c("grey nurse shark","spurdogs") #when M.average='mean', endless loop due to 
    # low fecundity (greynurse & spurdogs) or high mortality (due to k for sandbar)
    if(!Keep.species[l]%in%mean.donotwork)  
    {
      r.prior.dist_M.mean=with(List.sp[[l]],fun.rprior.dist(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,
                                                            K.sd=Growth.F$k.sd,LINF.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL,k.Linf.cor,
                                                            Amax=Max.age.F,
                                                            MAT=unlist(Age.50.mat),FecunditY=Fecundity,Cycle=Breed.cycle,
                                                            BWT=BwT,AWT=AwT,LO=Lzero*a_FL.to.TL+b_FL.to.TL)) #size vars as TL
    }
    #r.prior.dist_M.mean=within(r.prior.dist_M.mean, rm(nat.mort.sim))
    if(Keep.species[l]%in%mean.donotwork) r.prior.dist_M.mean=r.prior.dist_M.min
    
    store.species.r[[l]]=list(r.prior.dist_M.min=r.prior.dist_M.min,r.prior.dist_M.mean=r.prior.dist_M.mean)
  }
  toc()    
  clear.log('fun.rprior.dist')
  clear.log('fun.Leslie')
  if("plyr"%in%.packages()) detach("package:plyr", unload=TRUE)
  
  #Extract quantities
  store.species.r_M.min=lapply(store.species.r, function(x) x[[1]])
  store.species.r_M.mean=lapply(store.species.r, function(x) x[[2]])
  names(store.species.r_M.min)=names(store.species.r_M.mean)=Keep.species
  for(l in 1:N.sp)  
  {
    print(paste("extract r prior ","--",List.sp[[l]]$Name))
    PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
               capitalize(List.sp[[l]]$Name),"/",AssessYr,sep='')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
    PATH=paste(PATH,"/demography",sep='')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))   
    setwd(PATH)
    
    #export life history parameter distributions
    Nms=names(store.species.r_M.min[[l]]$Input.pars[[1]])
    LH.d=matrix(unlist(list.flatten(store.species.r_M.min[[l]]$Input.pars)),nrow=List.sp[[l]]$NsimSS,ncol=length(Nms),byrow = T)
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
    write.csv(with(store.species.r_M.min[[l]],data.frame(shape=shape,rate=rate,mean=mean,sd=sd)),'r.prior_M.min.csv',row.names = F)
    write.csv(with(store.species.r_M.mean[[l]],data.frame(shape=shape,rate=rate,mean=mean,sd=sd)),'r.prior_M.mean.csv',row.names = F)
    
    #export G
    out.G_M.min=with(store.species.r_M.min[[l]],data.frame(mean=mean(G),sd=sd(G)))
    write.csv(out.G_M.min,'G.prior_M.min.csv',row.names = F)
    store.species.G_M.min[[l]]=out.G_M.min
    
    out.G_M.mean=with(store.species.r_M.mean[[l]],data.frame(mean=mean(G),sd=sd(G)))
    write.csv(out.G_M.mean,'G.prior_M.mean.csv',row.names = F)
    store.species.G_M.mean[[l]]=out.G_M.mean
    
    #export M
    out.M=store.species.r_M.min[[l]]$M
    n.dim=max(unlist(lapply(out.M,length)))
    for(ss in 1:length(out.M))
    {
      a=out.M[[ss]]
      delta=n.dim-length(a)
      if(delta>0) out.M[[ss]]=c(a,rep(NA,delta))
      rm(a)
    }
    out.M=do.call(rbind,out.M)
    names(out.M)=0:(n.dim-1)
    store.species.M_M.min[[l]]=out.M
    write.csv(out.M,"M_M.min.csv",row.names=FALSE)
    
    out.M=store.species.r_M.mean[[l]]$M
    n.dim=max(unlist(lapply(out.M,length)))
    for(ss in 1:length(out.M))
    {
      a=out.M[[ss]]
      delta=n.dim-length(a)
      if(delta>0) out.M[[ss]]=c(a,rep(NA,delta))
      rm(a)
    }
    out.M=do.call(rbind,out.M)
    names(out.M)=0:(n.dim-1)
    store.species.M_M.mean[[l]]=out.M
    write.csv(out.M,"M_M.mean.csv",row.names=FALSE)
  }

  #Plot each M estimator 
  for(l in 1:N.sp)
  {
    store.species.r_M.min[[l]]$nat.mort.sim%>%
      ggplot(aes(x=Age, y=M.mean,color=Method))+
      geom_point()+
      geom_errorbar(aes(ymin=M.mean-M.sd, ymax=M.mean+M.sd))+
      facet_wrap(~Method,scales='free_y',ncol=2)+
      theme_PA()+
      theme(legend.position="none",
            axis.text.x = element_text(size=7,angle = 90, hjust=1))+
      ylab('M (+/- SD)')
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",AssessYr,"/demography/M_all.tiff",sep=''),
           width = 10,height = 6, dpi = 300, compression = "lzw")
    
  }
  
  
  #Compare Min and Mean Natural mortality
  for(l in 1:N.sp)
  {
    sd=apply(store.species.M_M.min[[l]],2,sd,na.rm=T)
    MN=colMeans(store.species.M_M.min[[l]],na.rm=T)
    M_x=1:length(MN)
    loe.mod=loess(MN~M_x)
    MN=predict(loe.mod)
    
    sd_mean=apply(store.species.M_M.mean[[l]],2,sd,na.rm=T)
    MN_mean=colMeans(store.species.M_M.mean[[l]],na.rm=T)
    M_x=1:length(MN_mean)
    loe.mod=loess(MN_mean~M_x)
    MN_mean=predict(loe.mod)
    
    
    fn.fig(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[l]]$Name),"/",AssessYr,"/demography/M_min.vs.mean",sep=''),2400,2000) 
    
    plot(MN,ylim=c(0,max(MN_mean)*1.1),main=Keep.species[l],
         ylab='M (+/- SD)',xlab='Age',pch=19)
    lines(MN-sd)
    lines(MN+sd)
    points(MN_mean,pch=19,col=2)
    lines(MN_mean-sd_mean,col=2)
    lines(MN_mean+sd_mean,col=2)
    legend('topright',c('Min','Mean'),pch=19,col=1:2,bty='n')
    dev.off()
    
  }
  
  rm(store.species.r)
}

if(!do.r.prior)
{
  for(l in 1:N.sp) 
  {
    hndl.dummy=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                     AssessYr,"/demography",sep='')
    store.species.r_M.min[[l]]=read.csv(paste(hndl.dummy,"/r.prior_M.min.csv",sep=''))
    store.species.G_M.min[[l]]=read.csv(paste(hndl.dummy,"/G.prior_M.min.csv",sep=''))
    store.species.M_M.min[[l]]=read.csv(paste(hndl.dummy,"/M_M.min.csv",sep=''))
    
    store.species.r_M.mean[[l]]=read.csv(paste(hndl.dummy,"/r.prior_M.mean.csv",sep=''))
    store.species.G_M.mean[[l]]=read.csv(paste(hndl.dummy,"/G.prior_M.mean.csv",sep=''))
    store.species.M_M.mean[[l]]=read.csv(paste(hndl.dummy,"/M_M.mean.csv",sep=''))
    
    rm(hndl.dummy)
  }
}


  #display priors   
for(l in 1:length(Lista.sp.outputs))
{
  STXSIZ=16
  if(names(Lista.sp.outputs)[l]=="Other.sp") STXSIZ=13
  fn.display.priors(d=store.species.r_M.min,
                    sp=Lista.sp.outputs[[l]],
                    XLAB=expression(paste(plain("Maximum intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),
                    XLIM=c(0,NA),
                    Strx.siz=STXSIZ)
  ggsave(paste(Rar.path,'/Prior_r_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = 12,height = 10,compression = "lzw")
}

  #compare M and r
Omit.these=c("Great hammerhead","Scalloped hammerhead","Smooth hammerhead",
             "Tiger shark","Dwarf sawfish","Freshwater sawfish",
             "Sandbar shark")
if(do.r.prior)
{

  CompR=data.frame(Name=names(store.species.r_M.min),r=NA,M=NA)
  for(l in 1:N.sp)
  {
    CompR$r[l]=store.species.r_M.min[[l]]$mean
    CompR$M[l]=mean(unlist(store.species.M_M.min[[l]]),na.rm=T)
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

  #FishLife approach
compare.FishLife=FALSE
if(compare.FishLife)
{
  library(FishLife)
  normal.scale=c('Temperature','rho','h','r','G') #Mean_pred parameters in normal scale
  fn.convert.normal.space=function(d,normal.scale)
  {
    id=match(normal.scale,names(d[[1]]$Mean_pred))
    Median=c(exp(d[[1]]$Mean_pred[-id]),d[[1]]$Mean_pred[id])
    Mean=c(exp(d[[1]]$Mean_pred[-id]+0.5*diag(d[[1]]$Cov_pred)[-id]),d[[1]]$Mean_pred[id]) #biased-corrected mean
    return(list(Median=Median,Mean=Mean))
  }
  sp.list=data.frame(Common.name=capitalize(Keep.species),
                     Genus=c("Squatina","Carcharhinus","Carcharhinus",
                             "Sphyrna","Pristis","Carcharias",
                             "Mustelus","Negaprion","Rhizoprionodon",
                             "Anoxypristis","Carcharhinus","Carcharhinus",
                             "Pristiophorus","Sphyrna","Isurus",
                             "Sphyrna","Carcharhinus","Squalus",
                             "Galeocerdo","Furgaleus","Orectolobus"),
                     Species=c("australis","brachyurus","obscurus",
                               "mokarran","zijsron","taurus",
                               "antarcticus","acutidens","acutus",
                               "cuspidata","amboinensis","plumbeus",
                               "cirratus","lewini","oxyrinchus",
                               "zygaena","brevipinna","megalops",
                               "cuvier","macki","maculatus"))
  Store.taxa.fishlife=vector('list',nrow(sp.list))
  names(Store.taxa.fishlife)=sp.list$Common.name
  for(s in 1:nrow(sp.list))
  {
    Store.taxa.fishlife[[s]]=Plot_taxa(Search_species(Genus=sp.list$Genus[s],
                                                      Species=sp.list$Species[s],add_ancestors=FALSE)$match_taxonomy,
                                       mfrow=c(3,3))
  }
  Store.keypars.fishlife=Store.taxa.fishlife
  keypars=c('Loo','K','tmax','tm','M','Lm','ln_Fmsy_over_M','ln_Fmsy','h','r','G')
  for(s in 1:nrow(sp.list))
  {
    Median.and.mean=fn.convert.normal.space(d=Store.taxa.fishlife[[s]],normal.scale)
    dummy=rbind(Median.and.mean$Median,Median.and.mean$Mean)%>%
      data.frame%>%
      dplyr::select(all_of(keypars))%>%
      summarise_all(list(mean))%>%
      rename(Fmsy_over_M=ln_Fmsy_over_M,
             Fmsy=ln_Fmsy)%>%
      mutate(Species=sp.list$Common.name[s])%>%
      relocate(Species)%>%
      mutate(across(where(is.numeric), round, 3))
    Store.keypars.fishlife[[s]]=dummy
    rm(dummy)
  }
  a=getwd()
  setwd(handl_OneDrive('Analyses/Population dynamics/FishLife'))
  write.csv(do.call(rbind,Store.keypars.fishlife),'FishLife parameters.csv',row.names = F) 
  
  #r prior with FishLife and SPMpriors
  #devtools::install_github("henning-winker/SPMpriors")
  library(SPMpriors)
  
  #Loo optional tuning prior dist for Linf c(mu, cv)
  #K optional tuning prior dist for brody K c(mu, cv)
  #tmax optional tuning prior dist for maximum age c(mu, cv)
  #tm	optional tuning prior dist for age maturity c(mu, cv)
  #Lm	optional tuning prior dist for length maturity c(mu, cv)
  #h optional uniform tuning prior dist for steepness h c(min, max)
  Store.r.fishlife=vector('list',nrow(sp.list))
  names(Store.r.fishlife)=sp.list$Common.name
  for(s in 1:nrow(sp.list))
  {
    attach(List.sp[[s]])
    K=Growth.F$k
    LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL
    K.sd=Growth.F$k.sd
    LINF.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL
    Amax=mean(Max.age.F)
    MAT=mean(unlist(Age.50.mat))
    detach(List.sp[[s]])
    dummy=flmvn_traits(Genus=sp.list$Genus[s],Species=sp.list$Species[s],
                       Loo=c(LINF,LINF.sd/LINF),K=c(K,K.sd/K),tmax=c(Amax,0.1),tm=c(MAT,0.1),
                       Plot=T,savepng = T)
    Store.r.fishlife[[s]]=dummy$mvnstk$r
    rm(K,LINF,K.sd,LINF.sd,Amax,MAT)
  }
  setwd(a)
  
}


 
#---10. Assign Resilience -----------------------------------------------------------------------
RESILIENCE=vector('list',N.sp)
names(RESILIENCE)=names(List.sp)
for(r in 1:length(RESILIENCE)) RESILIENCE[[r]]=Res.fn(store.species.r_M.min[[r]]$mean,Def="Haddon")
clear.log('Res.fn')

#---11. Extract experimental selectivity at age and at size-----------------------------------------------------------------------
#note: not used (selectivity estimated in SS3)
#      for species with no gillnet selectivity profile, set to closest species or family
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
  #selectivity at age
  N.sp.with.sel=which(sapply(Selectivity.at.age,function(x) !is.null(x)),TRUE)
  fn.fig(handl_OneDrive('Analyses\\Population dynamics\\growth.and.combined selectivity at age'),2400,2000) 
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
  
  #selectivity at size
  fn.fig(handl_OneDrive('Analyses\\Population dynamics\\growth.and.combined selectivity at TL'),2400,2000) 
  smart.par(n.plots=length(N.sp.with.sel),MAR=c(2,3,1,1),OMA=c(2.5,1,.05,2.5),MGP=c(1.8,.5,0))
  par(cex.lab=1.5,las=1)
  for(l in N.sp.with.sel)
  {
    with(Selectivity.at.totalength[[l]],plot(TL,Sel.combined,col=2,ylim=c(0,1),xlim=c(0,200),
                                      pch=19,main=capitalize(names(Selectivity.at.totalength)[l]),
                                      ylab='',xlab=""))
  }
  mtext("TL (cm)", side = 1, line = 1,outer=T)
  mtext("Selectivity", side = 2, line = -.5,las=3,col=2,outer=T)
  dev.off()
}


#---12. Extract catch rates-----------------------------------------------------------------------
Catch.rate.series=vector('list',N.sp) 
names(Catch.rate.series)=Keep.species
for(l in 1:N.sp) 
{
  Neim=Keep.species[l]
  if(any(c(grepl('annual.abundance',names(Species.data[[l]])),
           grepl('Srvy.FixSt',names(Species.data[[l]])),
           grepl('CPUE',names(Species.data[[l]])))))
  {
    dummy=list()
    if('CPUE_Pilbara.trawl'%in%names(Species.data[[l]]))
    {
      if(Neim%in%other_not.representative)
      {
        dumi.a=Species.data[[l]][-grep('Pilbara.trawl',names(Species.data[[l]]))]
        if(length(dumi.a)==0) dumi.a='dummy'
        Species.data[[l]]=dumi.a
      }else
      {
        d=Species.data[[l]]$CPUE_Pilbara.trawl
        if(Abundance.error.dist=="Lognormal") d$siv=d$CV 
        if(Abundance.error.dist=="Normal") d$siv=d$SE
        d=d%>%
          dplyr::select(-c( CV,SE))%>%
          rename(CV=siv)%>%
          mutate(yr.f=as.numeric(substr(Finyear,1,4)))
        if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$Other=d
        rm(d)
        
      }
    }
    if('Srvy.FixSt'%in%names(Species.data[[l]]))
    {
      if(Neim%in%survey_not.representative)
      {
        dumi.a=Species.data[[l]][-grep(paste(c('Srvy.FixSt','Srvy.FixSt_relative'),collapse='|'),names(Species.data[[l]]))]
        if(length(dumi.a)==0) dumi.a=NULL
        Species.data[[l]]=dumi.a
      }else
      {
        if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$Srvy.FixSt
        if(Abundance.error.dist=="Normal") d=Species.data[[l]]$Srvy.FixSt_relative
        d=d%>%
          mutate(Finyear=paste(yr-1,substr(yr,3,4),sep='-'),
                 yr.f=as.numeric(substr(Finyear,1,4)))%>%
          rename(Mean=MeAn,
                 LOW.CI=LowCI,
                 UP.CI=UppCI)%>%
          filter(yr.f<=Last.yr.ktch.numeric)
        if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$Survey=d
        rm(d)
      }
    }
    if('annual.abundance.NSF_relative'%in%names(Species.data[[l]]))
    {
      if(Neim%in%NSF_not.representative)
      {
        dumi.a=Species.data[[l]][-grep(paste(c('annual.abundance.NSF','annual.abundance.NSF_relative'),collapse='|'),names(Species.data[[l]]))]
        if(length(dumi.a)==0) dumi.a=NULL
        Species.data[[l]]=dumi.a
      }else
      {
        if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.NSF
        if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.NSF_relative
        d=d%>%
          rename(Finyear=FINYEAR)%>%
          mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
          filter(yr.f<=Last.yr.ktch.numeric)
        if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$NSF=d
        rm(d)
      }
    }
    if('CPUE_Observer_TDGDLF'%in%names(Species.data[[l]]))
    {
      d=Species.data[[l]]$CPUE_Observer_TDGDLF%>%
        rename(yr.f=year)%>%
      filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF_observer=d
      rm(d)
    }
    if('annual.abundance.basecase.monthly_relative'%in%names(Species.data[[l]]))
    {
      if(Neim%in%tdgdlf_not.representative)
      {
        dumi.a=Species.data[[l]][-grep(paste(c('annual.abundance.basecase.monthly_relative','annual.abundance.basecase.monthly'),collapse='|'),names(Species.data[[l]]))]
        if(length(dumi.a)==0) dumi.a=NULL
        Species.data[[l]]=dumi.a
      }else
      {
        if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.monthly
        if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.monthly_relative
        d=d%>%
          mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
          filter(yr.f<=Last.yr.ktch.numeric)
        Inf.CI=which(is.infinite(d$UP.CI))
        if(length(Inf.CI)>0) d[Inf.CI,match(c('Mean','CV','LOW.CI','UP.CI'),names(d))]=NA
        if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.monthly=d
        rm(d)
      }
    }
    if('annual.abundance.basecase.monthly.West'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.monthly.West
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.monthly.West_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.monthly.West=d
      rm(d)
    }
    if('annual.abundance.basecase.monthly.Zone1'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.monthly.Zone1
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.monthly.Zone1_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.monthly.Zone1=d
      rm(d)
    }
    if('annual.abundance.basecase.monthly.Zone2'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.monthly.Zone2
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.monthly.Zone2_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.monthly.Zone2=d
      rm(d)
    }
    if('annual.abundance.basecase.daily_relative'%in%names(Species.data[[l]]))
    {
      if(Neim%in%tdgdlf_not.representative)
      {
        dumi.a=Species.data[[l]][-grep(paste(c('annual.abundance.basecase.daily_relative','annual.abundance.basecase.daily'),collapse='|'),names(Species.data[[l]]))]
        if(length(dumi.a)==0) dumi.a=NULL
        Species.data[[l]]=dumi.a
      }else
      {
        if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.daily
        if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.daily_relative
        d=d%>%
          mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
          filter(yr.f<=Last.yr.ktch.numeric)
        if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.daily=d
        rm(d)
      }
    }
    if('annual.abundance.basecase.daily.West'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.daily.West
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.daily.West_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.daily.West=d
      rm(d)
    }
    if('annual.abundance.basecase.daily.Zone1'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.daily.Zone1
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.daily.Zone1_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.daily.Zone1=d
      rm(d)
    }
    if('annual.abundance.basecase.daily.Zone2'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.daily.Zone2
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.daily.Zone2_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.daily.Zone2=d
      rm(d)
    }
    
    #reset CV if in percentage
    if(length(dummy)>0)
    {
      for(x in 1:length(dummy))
      {
        if(min(dummy[[x]]$CV,na.rm=T)>1)  dummy[[x]]$CV=dummy[[x]]$CV/100
      }
      Catch.rate.series[[l]]=dummy
    }
  }
}

# Change Naturaliste survey from numbers to weights
for(i in 1:N.sp)
{
  if("Survey"%in%names(Catch.rate.series[[i]]))
  {
    ddim=Catch.rate.series[[i]]$Survey%>%
      left_join(Species.data[[i]]$Srvy.FixSt_size.absolute%>%
                  dplyr::select(yr,MeAn)%>%
                  rename(Size=MeAn),by='yr')%>%
      mutate(Size=ifelse(is.na(Size)&!is.na(Mean),mean(Size,na.rm=T),Size))
    if(grepl(paste(c("shovelnose","zebra shark"),collapse='|'),names(Catch.rate.series)[i])) Size.var='TL' else Size.var='FL'
    if(Size.var=="FL") ddim$Size=with(List.sp[[i]],ddim$Size*a_FL.to.TL+b_FL.to.TL)
    ddim=ddim%>%
      mutate(TW=with(List.sp[[i]],AwT*Size^BwT),
             Mean=Mean*TW,
             SE=CV*Mean,
             UP.CI=UP.CI*TW,
             LOW.CI=LOW.CI*TW)%>%
      dplyr::select(names(Catch.rate.series[[i]]$Survey))
    
    fn.fig(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[i]]$Name),"/",AssessYr,'/1_Inputs/Visualise data/CPUE_Survey in weight',sep=''),
           1800,2000) 
    par(mfcol=c(2,1))
    with(Catch.rate.series[[i]]$Survey,plot(yr,Mean,main='Survey in numbers',pch=19,
                                            ylim=c(0,max(UP.CI,na.rm=T))))
    with(Catch.rate.series[[i]]$Survey,segments(x0=yr, y0=LOW.CI, x1 = yr, y1 = UP.CI))
    with(ddim,plot(yr,Mean,pch=19,main='Converted to total weight (kg)',ylim=c(0,max(ddim$UP.CI,na.rm=T))))
    with(ddim,segments(x0=yr, y0=LOW.CI, x1 = yr, y1 = UP.CI))
    dev.off()
    
    Catch.rate.series[[i]]$Survey=ddim
  }
}

# get SE if only CV available and define if using CV or SE in assessment model  
for(i in 1:N.sp)
{
  if(!is.null(Catch.rate.series[[i]]))
  {
    for(p in 1:length(Catch.rate.series[[i]]))
    {
      if(!'SE'%in%names(Catch.rate.series[[i]][[p]]))
      {
        Catch.rate.series[[i]][[p]]=Catch.rate.series[[i]][[p]]%>%
                                mutate(SE=(Mean-LOW.CI)/1.96)
        if(Abundance.error.dist=='Normal')
        {
          Catch.rate.series[[i]][[p]]=Catch.rate.series[[i]][[p]]%>%
                                dplyr::select(-CV)%>%
                                rename(CV=SE)
        }
      }
    }
  }
  
}

#consistency between daily cpue and daily mean weight
for(l in 1:N.sp)
{
  if(is.null(Catch.rate.series[[l]]) & "annual.mean.size"%in%names(Species.data[[l]]))
  {
    Species.data[[l]]=Species.data[[l]][-match("annual.mean.size",names(Species.data[[l]]))] 
  }
}

#---13. Calculate Steepness ----------------------------------------------------------------------- 
  #calculate prior
store.species.steepness_M.mean=vector('list',N.sp)
names(store.species.steepness_M.mean)=Keep.species
store.species.alpha_M.min=store.species.alpha_M.mean=store.species.steepness_M.min=store.species.steepness_M.mean
if(do.steepness)   #0.005 sec per iteration per species
{
  store.species.h=vector('list',N.sp)
  names(store.species.h)=Keep.species
  tic()
  for(l in 1:N.sp)
  {
    print(paste('Steepness calculation -------------',Keep.species[l]))  
    SEL=Selectivity.at.age[[l]]$Sel.combined  #selectivity not used in h calculation as F.mult =0 so set to dummy if sel not available
    if(is.null(SEL)) SEL=rep(0,List.sp[[l]]$Max.age.F[2])
    
    if(is.na(List.sp[[l]]$Growth.F$FL_inf.sd)) List.sp[[l]]$Growth.F$FL_inf.sd=0.038*List.sp[[l]]$Growth.F$FL_inf  
    if(is.na(List.sp[[l]]$Growth.F$k.sd)) List.sp[[l]]$Growth.F$k.sd=0.088*List.sp[[l]]$Growth.F$k     
    
    #Fishing mortality set at 0 so selectivity has no effect
    k.Linf.cor=List.sp[[l]]$k.Linf.cor
    
    RESAMP="YES"
    linear.fec="NO"
    if(names(List.sp)[l]%in%c("angel sharks","grey nurse shark","spurdogs")) RESAMP="NO"   #endless loop due to life history combos <0.2
    
    GET.all.Ms=FALSE
    
    M.averaging<<-"mean"   #'min' yields too high h values for all species
    if(names(List.sp)[l]%in%c("milk shark")) M.averaging<<-"min"  #way too high M if using average
    steepNs_M.mean=with(List.sp[[l]],fun.steepness(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,
                                                   Linf.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL,k.sd=Growth.F$k.sd,
                                                   first.age=0,sel.age=SEL,F.mult=0,
                                                   Amax=Max.age.F,MAT=unlist(Age.50.mat),
                                                   FecunditY=Fecundity,Cycle=Breed.cycle,
                                                   sexratio=0.5,spawn.time = 0,
                                                   AWT=AwT,BWT=BwT,LO=Lzero*a_FL.to.TL+b_FL.to.TL,
                                                   Resamp=RESAMP,simsout=SSS.sims))
    
    M.averaging<<-"min"
    steepNs_M.min=with(List.sp[[l]],fun.steepness(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,
                                                  Linf.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL,k.sd=Growth.F$k.sd,
                                                  first.age=0,sel.age=SEL,F.mult=0,
                                                  Amax=Max.age.F,MAT=unlist(Age.50.mat),
                                                  FecunditY=Fecundity,Cycle=Breed.cycle,
                                                  sexratio=0.5,spawn.time = 0,
                                                  AWT=AwT,BWT=BwT,LO=Lzero*a_FL.to.TL+b_FL.to.TL,
                                                  Resamp=RESAMP,simsout=SSS.sims))
    rm(RESAMP,M.averaging,linear.fec)
    store.species.h[[l]]=list(steepNs_M.mean=steepNs_M.mean,steepNs_M.min=steepNs_M.min)
  }
  toc() 

  clear.log('fun.steepness')
  clear.log('Alpha.Brooks')
  
  #Extract quantities
  steepNs_M.mean=lapply(store.species.h, function(x) x[[1]])
  steepNs_M.min=lapply(store.species.h, function(x) x[[2]])
  names(steepNs_M.mean)=names(steepNs_M.min)=Keep.species
  for(l in 1:N.sp)
  {
    print(paste("extract steepness value ","--",List.sp[[l]]$Name))
    PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
               capitalize(List.sp[[l]]$Name),"/",AssessYr,"/steepness",sep='')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
    setwd(PATH)
    
    #export h, alpha and M
    steepNs=steepNs_M.mean[[l]]
    out.h=with(steepNs,data.frame(mean=mean,sd=sd))
    write.csv(out.h,'h.prior_M.mean.csv',row.names = F) 
    store.species.steepness_M.mean[[l]]=out.h
    write.csv(steepNs$Alpha,'Alpha_M.mean.csv',row.names = F) 
    store.species.alpha_M.mean[[l]]=steepNs$Alpha
    out.M=steepNs$M
    n.dim=max(unlist(lapply(out.M,length)))
    for(ss in 1:length(out.M))
    {
      a=out.M[[ss]]
      delta=n.dim-length(a)
      if(delta>0) out.M[[ss]]=c(a,rep(NA,delta))
      rm(a)
    }
    out.M=do.call(rbind,out.M)
    names(out.M)=0:(n.dim-1)
    write.csv(out.M,"M_M.mean.csv",row.names=FALSE)
    write.csv(steepNs$Runs,"Life.history_M.mean.csv",row.names=FALSE)
    
    steepNs=steepNs_M.min[[l]]
    out.h=with(steepNs,data.frame(mean=mean,sd=sd))
    write.csv(out.h,'h.prior_M.min.csv',row.names = F) 
    store.species.steepness_M.min[[l]]=out.h
    write.csv(steepNs$Alpha,'Alpha_M.min.csv',row.names = F) 
    store.species.alpha_M.min[[l]]=steepNs$Alpha
    out.M=steepNs$M
    n.dim=max(unlist(lapply(out.M,length)))
    for(ss in 1:length(out.M))
    {
      a=out.M[[ss]]
      delta=n.dim-length(a)
      if(delta>0) out.M[[ss]]=c(a,rep(NA,delta))
      rm(a)
    }
    out.M=do.call(rbind,out.M)
    names(out.M)=0:(n.dim-1)
    write.csv(out.M,"M_M.min.csv",row.names=FALSE)
    write.csv(steepNs$Runs,"Life.history_M.min.csv",row.names=FALSE)
    
    rm(steepNs,out.M)
    
  }
  
  rm(store.species.h)
}


if(!do.steepness)
{
  for(l in 1: N.sp)
  {
    store.species.steepness_M.mean[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                                AssessYr,"/steepness/h.prior_M.mean.csv",sep=''))
    store.species.alpha_M.mean[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                            AssessYr,"/steepness/Alpha_M.mean.csv",sep=''))
    
    store.species.steepness_M.min[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                                       AssessYr,"/steepness/h.prior_M.min.csv",sep=''))
    store.species.alpha_M.min[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                                   AssessYr,"/steepness/Alpha_M.min.csv",sep=''))
  }
}


  #compare steepness and r
if(do.steepness)
{
  Omit.these.h=c("Great hammerhead","Scalloped hammerhead","Smooth hammerhead","Tiger shark",
                 "Dwarf sawfish","Freshwater sawfish","Gummy shark",
                 'Green sawfish',"Wobbegongs","Spinner shark")
  CompR=data.frame(Name=names(store.species.steepness_M.mean),
                   h=unlist(sapply(store.species.steepness_M.mean, `[`, 1)),
                   h.sd=unlist(sapply(store.species.steepness_M.mean, `[`, 2)),
                   r=unlist(sapply(store.species.r_M.min, `[`, 3)),
                   r.sd=unlist(sapply(store.species.r_M.min, `[`, 4)))
  rownames(CompR)=NULL
  COL=rgb(.1,.2,.8,alpha=.45)
  CompR=CompR%>%
    arrange(r)%>%
    mutate(Name=capitalize(Name))

  Mod=lm(h~r,data=CompR%>%filter(!Name%in%Omit.these.h))
  ResSD=sqrt(sum((CompR%>%filter(!Name%in%Omit.these.h)%>%pull(h)-predict(Mod))^2)/(length(predict(Mod))-2))
  write.csv(data.frame(intercept=coef(Mod)[1],slope=coef(Mod)[2],CV=ResSD),
            handl_OneDrive('Analyses/Population dynamics/Steepness_vs_r_coeff.csv'),row.names = F)
  
  my_formula=y ~ x
  p=CompR%>%
    ggplot(aes(r, h, label = Name)) +
    geom_point(shape = 21, size = 5,fill=COL) + 
    geom_smooth(method = "lm", data = CompR%>%filter(!Name%in%Omit.these.h),
                se = F, fullrange = TRUE,colour="red")+
    stat_poly_eq(data = CompR%>%filter(!Name%in%Omit.these.h),
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
    ylim(0,1)+
    geom_errorbar(aes(ymin=h-h.sd, ymax=h+h.sd),colour=COL)+
    geom_errorbarh(aes(xmin=r-r.sd, xmax=r+r.sd),colour=COL)
  p
  ggsave(handl_OneDrive('Analyses/Population dynamics/Steepness_vs_r.tiff'), 
         width = 8,height = 8, dpi = 300, compression = "lzw")
  
}


#Recalculate steepness for species with too high/low h using linear model of r on h
#note: some h values deemed too high (following E Cortes discussion); some deemed too low (life history mispecificaton)
#       SEDAR sandbar h=0.3; SEDAR dusky h [0.25;0.35]
Mod.Pred=read.csv(handl_OneDrive('Analyses/Population dynamics/Steepness_vs_r_coeff.csv'))

store.species.steepness.S2=fn.get.stuff.from.list(store.species.steepness_M.mean,"mean")
h_too.long.converge=c("green sawfish","wobbegongs") #SSS for greens and wobbies taking too long to converge with original h
h_too.high=c("great hammerhead","scalloped hammerhead","smooth hammerhead","tiger shark",
             h_too.long.converge) 
h_too.low=c("spurdogs","sandbar shark")  
dis.sp.h=c(h_too.high,h_too.low)
if("dwarf sawfish" %in% Keep.species) dis.sp.h=c(dis.sp.h,"dwarf sawfish","freshwater sawfish")
for(s in 1:length(dis.sp.h))
{
  set.seed(1234)
  id=match(dis.sp.h[s],names(store.species.steepness_M.mean))
  new.h=store.species.r_M.min[[id]]$mean*Mod.Pred$slope+Mod.Pred$intercept
  store.species.steepness.S2[[id]]=rnorm(1,new.h,new.h/100)
}
store.species.steepness_M.mean$spurdogs$sd=store.species.steepness_M.min$spurdogs$sd/2 #M.mean yielded NA as life history results in h<0.2
#store.species.steepness.S2$`lemon shark`=0.3 #decrease to allow SS3 to fit, too low M and too high h otherwise

#display h2 priors
fn.display.steepness=function(d,d.h2,sp,XLAB,XLIM)
{
  for(z in 1:length(d)) d[[z]]$mean=d.h2[[z]]
  dummy=lapply(d[sp],function(x) rtruncnorm(1e4,a=Min.h.shark, b=Max.h.shark,mean=x$mean,sd=max(x$sd,0.01)))
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
if(do.steepness)
{
  for(l in 1:length(Lista.sp.outputs))
  {
    fn.display.steepness(d=store.species.steepness_M.mean,
                         d.h2=store.species.steepness.S2,
                         sp=Lista.sp.outputs[[l]],
                         XLAB="Steepness (h)",
                         XLIM=c(0.2,NA))
    ggsave(paste(Rar.path,'/Prior_steepness_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = 12,height = 10,compression = "lzw")
  }
}
clear.log('fn.display.priors')
clear.log('M.fun')
clear.log('fn.display.steepness')

do.this=FALSE
if(do.this)
{
  A=vector('list',N.sp)
  for(l in 1:N.sp)A[[l]]=data.frame(Species=Keep.species[l],r=store.species.r_M.min[[l]]$mean,h=store.species.steepness.S2[[l]])
  do.call(rbind,A)%>%
    ggplot(aes(r,h, label =Species))+
    geom_point()+geom_text_repel(segment.colour='black',col='black',box.padding = 0.5)
}

#---14. Calculate scaler for Fmsy.to.M and Bmsy.K proxy ----------------------------------------------------------------------- 

  #14.1 Fmsy.to.M relationship (Fmsy= scaler x M) for DBSRA
Cortes.Brooks.2018=function(alpha)  #source: Cortes & Brooks 2018
{
  if(alpha<=2.67)           Fmsy.to.M.scaler=0.2
  if(alpha>2.67 & alpha<=6) Fmsy.to.M.scaler=0.5
  if(alpha>6)               Fmsy.to.M.scaler=0.8
  
  return(Fmsy.to.M.scaler)
}
Fmsy.M.scaler=vector('list',N.sp)
names(Fmsy.M.scaler)=Keep.species
for(l in 1:N.sp) Fmsy.M.scaler[[l]]=Cortes.Brooks.2018(alpha=median(unlist(store.species.alpha_M.mean[[l]])))
for(l in 1:length(dis.sp.h)) #reset dis.sp.h consistently with h resetting
{
  s=match(dis.sp.h[l],names(Fmsy.M.scaler))
  Fmsy.M.scaler[[s]]=0.5
}

clear.log('store.species.alpha_M.mean')
clear.log('store.species.alpha_M.min')

#1. Species-specific proxy to Bmsy.K based on Cortes et al 2012
R=function(r,G) 0.633-0.187*log(r*G) 
BmsyK.species=data.frame(Species=names(store.species.r_M.min),r=NA,G=NA)
for(l in 1:N.sp)
{
  BmsyK.species$r[l]=store.species.r_M.min[[l]]$mean
  BmsyK.species$G[l]=store.species.G_M.min[[l]]$mean
}
BmsyK.species=BmsyK.species%>%
  mutate(BmsyK=R(r,G),
         Method='Cortes et al 2012')%>%
  rename(Gen.time=G)%>%
  mutate(m=NA)
BmsyK.species=BmsyK.species%>%
                mutate(BmsyK=ifelse(Species=='milk shark',mean(BmsyK.species%>%filter(Species%in%c('gummy shark','narrow sawfish'))%>%pull(BmsyK)),
                                    BmsyK)) #Cortes method yields nonsense large value for milk shark so set

fn.minimise <- function(m)  return(abs(m^(-1/(m-1)) - R))
for(l in 1:N.sp) #get shape parameter of production models based on Cortes et al 2012 (m=2 is Schaefer)
{
  R=BmsyK.species$BmsyK[l]
  x <- optimize(fn.minimise,c(0,4))
  BmsyK.species$m[l]=x$minimum 
}

#---15. Export .pin files and define modelling arguments  -----
#note: For integrated model, all pars calculated in log space.
#       ln_RZERO is in 1,000 individuals so do 10 times the largest catch divided by 
#       average weight and divided by 1,000. Best units to work in are 1,000 individuals for 
#       numbers, catch in tonnes and length-weight in kg as all cancels out and predicted biomass
#       end up being in tonnes
fn.source1("Pin_file_and_model_arguments.r")


#---16. Export .dat (figures and assessment input files to each species' folder) and some prelim analysis----- 
library(r4ss)
if(First.run=="YES") fn.source1("Organise data.R")

#Compare empirical selectivity with observed size composition
do.this='NO'
if(do.this=="YES")
{
  fn.compare.sel.size.comp=function(Title,Sel,size,FL_TL,MX)
  {
    size=size%>%
      mutate(TL=2.5+TL.bins.cm*floor((FL*FL_TL$slope+FL_TL$inter)/TL.bins.cm))%>%
      group_by(TL,year)%>%
      count()%>%
      group_by(year)%>%
      mutate(n.rel=n/max(n),
             N=sum(n),
             Group=paste(year,' (n=',N,')',sep=''))
    
    p=size%>%
      ggplot()+
      geom_col(aes(TL,n.rel))+
      facet_wrap(~Group)+
      geom_line(data=Sel,aes(TL,Sel),col='red')+
      xlim(40,MX)+
      ylab('Selectivity / size frequency')+
      xlab('Total length (cm)')+
      ggtitle(paste(Title,' (selectivity derived for ',unique(Sel$type),')',sep=''))
    
    print(p)
  }
  for(l in 1:N.sp)
  {
    if(!is.null(Selectivity.at.totalength[[l]]))
    {
      attach(List.sp[[l]])
      if(MN.SZE=="size.at.birth") Min.size.bin=10*round((Lzero*a_FL.to.TL+b_FL.to.TL)/10)
      if(MN.SZE==0) Min.size.bin=0
      MaxLen= 10*round(TLmax*Plus.gp.size/10)
      lbnd = seq(Min.size.bin,MaxLen - TL.bins.cm, TL.bins.cm)
      ubnd = lbnd + TL.bins.cm
      midpt = lbnd + (TL.bins.cm/2)
      
      #TDGLDF
      TDGDLF.size.comp=Species.data[[l]][grep(paste(c('6.5.inch.raw','.7.inch.raw'),collapse="|"),
                                              names(Species.data[[l]]))]
      names(TDGDLF.size.comp)=str_replace_all(names(TDGDLF.size.comp), 
                                              paste(c("\\Size_composition_+",".inch.raw"),collapse='|'), "")
      if(length(TDGDLF.size.comp)>0)
      {
        pdf(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),
                  "/",AssessYr,"/1_Inputs/Visualise data","/Compare sel and obs size comp.pdf",sep=''))
        for(s in 1:length(TDGDLF.size.comp))
        {
          if(grepl('6.5',names(TDGDLF.size.comp)[s])) pull.this='X16.5'
          if(grepl('7',names(TDGDLF.size.comp)[s])) pull.this='X17.8'
          fn.compare.sel.size.comp(Title=names(TDGDLF.size.comp)[s],
                                   Sel=data.frame(type=Selectivity.at.totalength[[l]]%>%pull(type),
                                                  TL=Selectivity.at.totalength[[l]]%>%pull(TL),
                                                  Sel=Selectivity.at.totalength[[l]]%>%
                                                    pull(pull.this)),
                                   size=TDGDLF.size.comp[[s]],
                                   FL_TL=data.frame(inter=b_FL.to.TL,
                                                    slope=a_FL.to.TL),
                                   MX=10*round(TLmax/10))
        }
        dev.off()
      }
      
      #NSF (selectivity set at maturity)
      NSF.size.comp=Species.data[[l]][grep(paste(c('NSF.LONGLINE'),collapse="|"),
                                           names(Species.data[[l]]))]
      if(length(NSF.size.comp)>0)
      {
        pdf(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),
                  "/",AssessYr,"/1_Inputs/Visualise data","/Compare sel and obs size comp_NSF.pdf",sep=''))
          fn.compare.sel.size.comp(Title='NSF - longline',
                                   Sel=data.frame(type='Species',
                                                  TL=midpt,
                                                  Sel=round(1/(1+(exp(-log(19)*((midpt-TL.50.mat)/(TL.95.mat-TL.50.mat))))),3)),
                                   size=NSF.size.comp$Size_composition_NSF.LONGLINE,
                                   FL_TL=data.frame(inter=b_FL.to.TL,
                                                    slope=a_FL.to.TL),
                                   MX=10*round(TLmax/10))
        dev.off()
      }

      detach(List.sp[[l]])
      print(paste('-------------Compare selectivity and observed size comp for:',List.sp[[l]]$Name))
      
    }
  }
}

#Check if different CPUE indices show the same trend and if not biologically possible  
if(First.run=="YES")
{
  fn.source1("CPUE_correlation_and_realism.R")  
  require("ggrepel")
  require('GGally')
  for(l in 1:N.sp)
  {
    if(!is.null(Catch.rate.series[[l]]))
    {
      Neim=names(Catch.rate.series)[l]
      CPUE=compact(Catch.rate.series[[l]])
      DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
      if(length(DROP)>0)CPUE=CPUE[-DROP]
      if(Neim%in%survey_not.representative & "Survey"%in%names(CPUE)) CPUE=CPUE[-grep("Survey",names(CPUE))]
      if(Neim%in%NSF_not.representative & "NSF"%in%names(CPUE)) CPUE=CPUE[-grep("NSF",names(CPUE))]
      if(Neim%in%tdgdlf_not.representative & "TDGDLF"%in%names(CPUE)) CPUE=CPUE[-grep("TDGDLF",names(CPUE))]
      if(length(CPUE)>1)for(x in 1:length(CPUE))
      {
        dd=CPUE[[x]]%>%
          dplyr::select(Finyear,Mean)%>%
          dplyr::rename(Year=Finyear)%>%
          mutate(Year=as.numeric(substr(Year,1,4)),
                 Series=names(CPUE)[x])
        CPUE[[x]]=dd
      }
      if(length(CPUE)>1)
      {
        CPUE=do.call(rbind,CPUE)%>%
          tidyr::spread(Series,Mean)
        
        fn.cpue.corr(Dat=CPUE,
                     AREA="Western Australia",
                     WD=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(Neim),
                              "/",AssessYr,"/1_Inputs/Visualise data",sep=''),
                     R=store.species.r_M.min[[l]]$mean)
      }
    }
    print(paste("Displaying CPUE correlation for -----",names(Catch.rate.series)[l]))
  }
  detach("package:reshape", unload=TRUE)
  detach("package:FLCore", unload=TRUE)
  detach("package:plyr", unload=TRUE)
  clear.log('fn.cpue.corr')
}

#Check if observed FL is with Lo +/- CV and FLinf +/- 
#remove "Size.type"
for(i in 1:length(Species.data))
{
  nm.Dat=names(Species.data[[i]])
  iid=nm.Dat[fn.extract.dat.perl(STRING="Size_composition",nm.Dat)]
  if(length(iid)>0)
  {
    for(xx in 1:length(iid))
    {
      if("Size.type"%in%names(Species.data[[i]][[match(iid[xx],nm.Dat)]]))
      {
        Species.data[[i]][[match(iid[xx],nm.Dat)]]=Species.data[[i]][[match(iid[xx],nm.Dat)]]%>%
          dplyr::select(-Size.type)
      }
    }

  }
}
if(First.run=="YES")
{
  for(i in 1:length(Species.data))
  {
    print(paste("observed FL is with Lo +/- CV and FLinf --------",names(Species.data)[i]))
    if(any(grepl('Size_composition',names(Species.data[[i]]))))
    {
      d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                    names(Species.data[[i]]))]
      if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
      if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
      for(x in 1:length(d.list)) d.list[[x]]$Fleet=str_remove(str_remove(names(d.list)[x],'Size_composition_'),'.inch.raw')
      d.list=do.call(rbind,d.list)%>%
                mutate(Fleet=ifelse(grepl(paste(c('West','Zone'),collapse='|'),Fleet),'TDGDLF',Fleet),
                       Fleet=ifelse(Fleet=='TDGDLF' & year<=2005,'Southern.shark_1',
                             ifelse(Fleet=='TDGDLF' & year>2005,'Southern.shark_2',
                             ifelse(Fleet=='NSF.LONGLINE','Northern.shark',
                                    Fleet))))
      
      CV_young=List.sp[[i]]$Growth.CV_young
      CV_old=List.sp[[i]]$Growth.CV_old
      Lo=List.sp[[i]]$Lzero
      Lo.min=Lo-(CV_young*Lo)
      LinfF=List.sp[[i]]$Growth.F$FL_inf
      LinfM=List.sp[[i]]$Growth.M$FL_inf
      if(is.na(LinfM)) LinfM=LinfF
      if(LinfF==LinfM) LinfM=base::jitter(LinfM)
      LinfF.max=LinfF+(CV_old*LinfF)
      LinfM.max=LinfM+(CV_old*LinfM)
      nRows=length(unique(d.list$Fleet))
      nCols=1
      d.list%>%
        filter(!is.na(SEX))%>%
        ggplot(aes(FL,fill=SEX))+
        geom_histogram()+
        facet_wrap(~Fleet,nrow = nRows,scales='free_y')+
        xlim(Lo.min,max(LinfF.max,LinfM.max))+
        theme(legend.position="top",legend.title=element_blank())+
        geom_vline(xintercept=Lo,color = "black",size=1.25,alpha=.7)+
        geom_vline(xintercept=Lo.min,color = "black",alpha=.7,size=1.1,linetype = "dotted")+
        geom_vline(xintercept=LinfF,color = "#F8766D",size=1.25)+
        geom_vline(xintercept=LinfF.max,color = "#F8766D",size=1.1,linetype = "dotted")+
        geom_vline(xintercept=LinfM,color = "#00BFC4",size=1.25)+
        geom_vline(xintercept=LinfM.max,color = "#00BFC4",size=1.1,linetype = "dotted")
      
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[i]]$Name),
                   "/",AssessYr,"/1_Inputs/Visualise data","/Observed size vs Linf and Lo with CV.tiff",sep=''),
             width = 8,height = 8, dpi = 300, compression = "lzw")
    }
  }
}

# Check that meanbodywt used in SS occurs ~ at peak of selectivity
#note: SS Tiger shark model not fitting well if using meanbodywt
if(First.run=="YES")
{
  fn.source1("SS_selectivity functions.R")
  fun.check.mean.weight=function(TL,a,b,Mean.weight,SD.Mean.weight,Sel,NM)
  {
    PP=data.frame(TL=TL,TWT=a*TL^b)
    with(Sel,wrapper.fn(x=TL,a.dn=P_1,b.dn=P_2,c.dn=P_3,d.dn=P_4,e.dn=P_5,f.dn=P_6,a.log=1e5,b.log=0))
    MeanTL=PP$TL[which.min(abs(PP$TWT-Mean.weight))]
    SDMeanTL.low=PP$TL[which.min(abs(PP$TWT-(Mean.weight-SD.Mean.weight)))]
    SDMeanTL.high=PP$TL[which.min(abs(PP$TWT-(Mean.weight+SD.Mean.weight)))]
    abline(v=MeanTL,lwd=2)
    abline(v=SDMeanTL.low,lty=3,lwd=1.5,col='grey30')
    abline(v=SDMeanTL.high,lty=3,lwd=1.5,col='grey30')
    text(PP$TL[which.min(abs(PP$TWT-Mean.weight))],0,"TL at mean meanbodywt input for SS",pos=4,srt=90)
    mtext(NM)
  }
  
  for(i in 1:length(Species.data))
  {
    print(paste("meanbodywt used in SS occurs ~ at peak of selectivity --------",names(Species.data)[i]))
    if(any(grepl('annual.mean.size',names(Species.data[[i]]))))
    {
       tiff(file=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[i]]$Name),
                      "/",AssessYr,"/1_Inputs/Visualise data","/Meanbodywt VS selectivity used in SS.tiff",sep=''),
           width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
      fun.check.mean.weight(TL=with(List.sp[[i]],seq(round(Lzero*a_FL.to.TL+b_FL.to.TL),TLmax)),
                            a=List.sp[[i]]$AwT,
                            b=List.sp[[i]]$BwT,
                            Mean.weight=mean(Species.data[[i]]$annual.mean.size$mean),
                            SD.Mean.weight=sd(Species.data[[i]]$annual.mean.size$mean),
                            Sel=List.sp[[i]]$SS_selectivity%>%filter(Fleet=='Southern.shark_2'),
                            NM='')
      dev.off()
    }
  }
  clear.log('fun.check.mean.weight')
}

# Compare observed size comp and assumed SS selectivity 
if(First.run=="YES")
{
  fn.source1("SS_selectivity functions.R")
  fun.compare.sel.obs.size.comp=function(TL,Sel,size.comps,Flts)
  {
    smart.par(length(Flts),c(1,1,1,1),c(1,1,1,1),c(3, 1, 0))
    for(x in 1:length(Flts))
    {
      FLiT=Flts[x]
      if(FLiT=="Pilbara_Trawl") FLiT='Other'
      Sel=List.sp[[i]]$SS_selectivity%>%filter(Fleet==FLiT)
      Sel.ori=SS_selectivity_init_pars%>%filter(Species==names(List.sp)[i])
      a.log=1e5
      b.log=0
      if(all(is.na(Sel[,c('P_3', 'P_4', 'P_5', 'P_6')])))
      {
        a.log=Sel$P_1
        b.log=Sel$P_2
        Sel[,c('P_1', 'P_2','P_3', 'P_4', 'P_5', 'P_6')]=c(0,100,-10,-10,0,0) 

      }
      with(Sel,wrapper.fn(x=TL,a.dn=P_1,b.dn=P_2,c.dn=P_3,d.dn=P_4,e.dn=P_5,f.dn=P_6,a.log=a.log,b.log=b.log))
      
      dd=size.comps%>%
            filter(Fleet==Flts[x])
      if(nrow(dd)>2)
      {
        d <- density(dd$TL,adjust = 1.25)
        lines(d$x,d$y/max(d$y),col="forestgreen",lwd=2)
      }

      legend('bottomleft','size composition',lwd=2,col="forestgreen",bty='n')
      text(quantile(TL,.9),.9,Flts[x])
    }
  }
  for(i in 1:length(Species.data))
  {
    print(paste("Compare observed size comp and assumed SS selectivity --------",names(Species.data)[i]))
    if(any(grepl('Size_composition',names(Species.data[[i]]))))
    {
      d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                    names(Species.data[[i]]))]
      if(length(d.list)>0)
      {
        if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
        if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
        for(x in 1:length(d.list)) d.list[[x]]$Fleet=str_remove(str_remove(names(d.list)[x],'Size_composition_'),'.inch.raw')
        d.list=do.call(rbind,d.list)
        d.list=d.list%>%mutate(fleet=ifelse(grepl(paste(c('West','Zone'),collapse='|'),Fleet),'TDGDLF',Fleet),
                               Fleet=ifelse(fleet=='TDGDLF' & year<=2005,'Southern.shark_1',
                                            ifelse(fleet=='TDGDLF' & year>2005,'Southern.shark_2',
                                                   ifelse(fleet=='NSF.LONGLINE','Northern.shark',
                                                          fleet))),
                               TL=FL*List.sp[[i]]$a_FL.to.TL+List.sp[[i]]$b_FL.to.TL)
        
        tiff(file=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[i]]$Name),
                        "/",AssessYr,"/1_Inputs/Visualise data","/Observe size VS selectivity used in SS.tiff",sep=''),
             width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
        fun.compare.sel.obs.size.comp(TL=with(List.sp[[i]],seq(round(Lzero*a_FL.to.TL+b_FL.to.TL),TLmax)),
                                      Sel=List.sp[[i]]$SS_selectivity,
                                      size.comps=d.list,
                                      Flts=sort(unique(d.list$Fleet)))
        dev.off()
      }
    }
  }
  clear.log('fun.compare.sel.obs.size.comp')
}

#display cpues     
if(First.run=="YES")
{
  for(l in 1:N.sp)
  {
    if(!is.null(Catch.rate.series[[l]]))
    {
      fn.fig(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(List.sp[[l]]$Name),"/",AssessYr,
                   "/1_Inputs/Visualise data/All cpues",sep=''),2000,2000) 
      smart.par(n.plots=length(Catch.rate.series[[l]]),MAR=c(2,3,1,1),OMA=c(2.5,1,.05,2.5),MGP=c(1.8,.5,0))
      par(cex.lab=1.5,las=1)
      Mx.yr=max(sapply(Catch.rate.series[[l]], function(x) max(x$yr.f, na.rm=TRUE)))
      Min.yr=min(sapply(Catch.rate.series[[l]], function(x) min(x$yr.f, na.rm=TRUE)))
      
      for(i in 1:length(Catch.rate.series[[l]]))
      {
        with(Catch.rate.series[[l]][[i]],{
          plot(yr.f,Mean,col='orange',xlim=c(Min.yr,Mx.yr),ylim=c(0,max(UP.CI,na.rm=T)),pch=19,cex=1.15,
               main=names(Catch.rate.series[[l]])[i],ylab='',xlab="")
          segments(yr.f,LOW.CI,yr.f,UP.CI,col='orange')
        })
      }
      mtext("Financial year", side = 1, line = 1,outer=T)
      mtext("CPUE", side = 2, line = -1,las=3,outer=T)
      dev.off()
    }
  }
}

# Get sex ratio by zone
if(First.run=="YES")
{
  Dis.size.data=SS3_fleet.size.comp.used
  Sex.ratio.zone=vector('list',N.sp)
  names(Sex.ratio.zone)=Keep.species
  for(i in 1:N.sp)
  {
    print(paste("Sex ratio by zone from size comp for -----",Keep.species[i]))
    id=grep(paste(Dis.size.data,collapse="|"),names(Species.data[[i]]))
    if(length(id)>0)
    {
      Name=Keep.species[i]
      HandL=handl_OneDrive("Analyses/Population dynamics/1.")
      DiR=paste(HandL,capitalize(Name),"/",AssessYr,"/1_Inputs/Visualise data",sep='')
      Sex.ratio.zone[[i]]=fn.ktch.sex.ratio.zone(size.data=Species.data[[i]][id])
      ggsave(paste(DiR,'Sex ratio by zone from size composition.tiff',sep='/'), width = 6,height = 6, dpi = 300, compression = "lzw")
      
    }
  }
}

#Display conditional age at length
if(First.run=="YES")
{
  for(i in 1:N.sp)
  {
    print(paste('Show conditional age at length for ------',Keep.species[i]))
    if(any(grepl('age_length',names(Species.data[[i]]))))
    {
      ktch=KtCh%>%
        filter(Name==Keep.species[i])%>%
        mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                             ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                      'TEP_greynurse','TEP_dusky','Discards_TDGDLF'),'Southern.shark',
                                    'Other')))
      Flits.name=sort(unique(ktch$Fishry))  
      Flits=1:length(Flits.name)
      names(Flits)=Flits.name
      a=List.sp[[i]]$a_FL.to.TL
      b=List.sp[[i]]$b_FL.to.TL
      dd=Species.data[[i]]$age_length%>%
        mutate(TL=FL*a+b,
               LbinLo=TL.bins.cm*floor(TL/TL.bins.cm),
               LbinHi=TL.bins.cm*floor(TL/TL.bins.cm))%>%
        group_by(year,Sex,LbinLo)%>%
        mutate(Nsamps=n())%>%
        ungroup()%>%
        mutate(month=7,
               Sex=ifelse(Sex=='Male',2,ifelse(Sex=='Female',1,0)),
               part=0,
               ageErr=-1,   #- assumes no ageing error
               Fleet=Flits[match('Southern.shark',names(Flits))])%>%
        group_by(year,Sex,LbinLo,Age)%>%
        mutate(N=n())%>%
        ungroup()%>%
        distinct(year,Sex,Age,LbinLo,.keep_all = T)
      dd%>%
        ggplot(aes(Age,TL,color=year))+
        geom_point()+
        facet_wrap(~Sex,ncol=1)
      HandL=handl_OneDrive("Analyses/Population dynamics/1.")
      DiR=paste(HandL,capitalize(Keep.species[i]),"/",AssessYr,"/1_Inputs/Visualise data",sep='')
      ggsave(paste(DiR,'Conditional age at length.tiff',sep='/'), width = 6,height = 6, dpi = 300, compression = "lzw")
      
    }
  }
}

#Compare r dist and life history invariant r~ 1.5-2 M
if(First.run=="YES")
{
  r.within.m.range=vector('list',N.sp)
  for(i in 1:N.sp)
  {
    
    r.within.m.range[[i]]=data.frame(Species=capitalize(names(List.sp)[i]),
                                     r=seq(0,1,.01))%>%
      mutate(r.prob=dnorm(seq(0,1,.01),
                          mean=List.sp[[i]]$Sens.test$JABBA$r[1],
                          sd=List.sp[[i]]$Sens.test$JABBA$r.sd[1]),
             m_1.5=1.5*List.sp[[i]]$Sens.test$DBSRA$Mmean[1],
             m_2=2*List.sp[[i]]$Sens.test$DBSRA$Mmean[1])
  }
  d=do.call(rbind,r.within.m.range)
  d%>%
    ggplot(aes(r,r.prob))+
    geom_polygon()+
    geom_vline(aes(xintercept=m_1.5),color=2,size=.9)+
    geom_vline(aes(xintercept=m_2),color=2,size=.9)+
    facet_wrap(~Species,scales='free')+
    theme_PA(strx.siz=9)+ylab('')
  ggsave(handl_OneDrive('Analyses/Population dynamics/M1.5_M2_vs_r dist.tiff'), width = 9,height = 8, dpi = 300, compression = "lzw")
  
}

#---17. Display catches by fishery & display life history ----
Tot.ktch=KtCh %>%      
  mutate(
    Type = case_when(
      FishCubeCode=='WRL'~'WRL',
      FishCubeCode=='WTB'~'WTB',
      FishCubeCode%in%c('Discards_TDGDLF','TEP')~'TEP',
      FishCubeCode=='Recreational'~'Recreational',
      FishCubeCode=='SA MSF'~'SA MSF',
      FishCubeCode=='NSW fisheries'~'NSW fisheries',
      FishCubeCode=='NT'~'NT',
      FishCubeCode=='GAB'~'GAB Trawl',
      FishCubeCode=='Indo'~'Indonesia',
      FishCubeCode=='Taiwan'~'Taiwan',
      FishCubeCode=='Historic'~'Historic',
      FishCubeCode%in%c('JASDGDL','WCDGDL','C070','OAWC')~'TDGDLF',
      FishCubeCode%in%c('JANS','OANCGC','WANCS')~'NSF',   
      TRUE  ~ "Other WA Commercial"),
    Type=factor(Type,levels=c("Recreational","TDGDLF","NSF",
                              "Historic","Other WA Commercial","TEP","WRL",
                              "NT","NSW fisheries","SA MSF","GAB Trawl","WTB",
                              "Taiwan","Indonesia")))
all.yrs=min(Tot.ktch$finyear):max(KtCh$finyear)
Fishry.type=levels(Tot.ktch$Type)
COLs.type=colfunc(length(Fishry.type))
names(COLs.type)=Fishry.type
All.N.sp=sort(unique(Tot.ktch$Name))
All.N.sp=subset(All.N.sp,!All.N.sp%in%names(Indicator.species))

if(First.run=="YES")
{
  for(l in 1:length(Lista.sp.outputs))
  {
  SIZ=2
  WID=12
  HEI=10
  t.siz=15
  if(length(Lista.sp.outputs[[l]])>8)
  {
    SIZ=1.5 
    t.siz=13
    WID=14
  }
    
  if(length(Lista.sp.outputs[[l]])<4) SIZ=3
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
    theme_PA(strx.siz=15,leg.siz=14,axs.t.siz=t.siz,axs.T.siz=18)+
    ylab("Total catch (tonnes)")+xlab("Financial year")+
    theme(legend.position="top",
          legend.title = element_blank(),
          legend.key=element_blank(),
          plot.margin = margin(0.1,.5,0.1,0.1, "cm"))+
    guides(colour = guide_legend(override.aes = list(size=5,linetype = 0)))
  
  ggsave(paste(Rar.path,'/Catch_all_species_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = WID,height = HEI,compression = "lzw")
  
}
}


#Show life history  
if(First.run=="YES")
{
  Lmx<- expression("TL"['Max'])
  L50<- expression("TL"[50])
  L95<- expression("TL"[95])
  L0<- expression("TL"[0])
  A50<- expression("Age"[50])
  Amx<- expression("Age"['Max'])
  for(l in 1:N.sp)
  {
    print(paste('Displaying biology for ------', names(List.sp)[l]))
    with(List.sp[[l]],{
      eig=0:Max.age.F[2]
      lengz=Lzero:TLmax
      TLo=Lzero*a_FL.to.TL+b_FL.to.TL
      plt=vector('list',length=5)
      plt[[1]]=rbind(data.frame(Age=eig,
                                TL=TLo+(Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL-TLo)*(1-exp(-Growth.F$k*eig)),
                                Sex='F'),
                     data.frame(Age=eig,
                                TL=TLo+(Growth.M$FL_inf*a_FL.to.TL+b_FL.to.TL-TLo)*(1-exp(-Growth.M$k*eig)),
                                Sex='M'))%>%
        ggplot(aes(Age,TL,color=Sex))+
        geom_line(size=1.25)+ylab("TL (cm)")+
        theme_PA()+
        annotate("text", x = 1.5, y = TLmax*1.04,parse = T, label = as.character(Lmx))+
        geom_hline(yintercept=TLmax, linetype="dashed",alpha=0.5)+
        annotate("text", x = 1.1, y = TLo*.8,parse = T, label = as.character(L0))+
        geom_hline(yintercept=TLo, linetype="dashed",alpha=0.5)+
        annotate("text", x = mean(Max.age.F)-2, y = 0,parse = T, label = as.character(Amx),color='darkorange4')+
        geom_vline(xintercept=Max.age.F[1], linetype="dashed",alpha=0.5,color='darkorange4')+
        geom_vline(xintercept=Max.age.F[2], linetype="dashed",alpha=0.5,color='darkorange4')+
        annotate("text", x = mean(Age.50.mat), y = 0, parse = T, label = as.character(A50),col='chartreuse4')+
        geom_vline(xintercept=Age.50.mat[1], linetype="dashed",alpha=0.5,col='chartreuse4')+
        geom_vline(xintercept=Age.50.mat[2], linetype="dashed",alpha=0.5,col='chartreuse4')+
        annotate("text", y = TL.50.mat*.95, x = 1.1,parse = T, label = as.character(L50),col='brown2')+
        geom_hline(yintercept=TL.50.mat, linetype="dashed",alpha=0.5,color='brown2')+
        annotate("text", y = TL.95.mat*.95, x = 1.1,parse = T, label =as.character(L95) ,col='brown2')+
        geom_hline(yintercept=TL.95.mat, linetype="dashed",alpha=0.5,color='brown2')+
        theme(legend.title=element_blank(),
              legend.position="top")
      
      plt[[2]]=data.frame(TL=lengz,
                          Maturity=1/(1+exp(-log(19)*((lengz-TL.50.mat)/(TL.95.mat-TL.50.mat)))))%>%
        ggplot(aes(TL,Maturity))+
        geom_line(size=1.25,color="#F8766D")+xlab("TL (cm)")+ylab('Proportion mature')+
        theme_PA()+
        annotate("text", x = TL.50.mat*.75, y = .8, label = paste('Fecundity [',Fecundity[1],',',Fecundity[2],']',sep=''))+
        annotate("text", x = TL.50.mat*.75, y = .85, label = paste('Cycle [',Breed.cycle[1],',',Breed.cycle[2],']',sep=''))
      
      plt[[3]]=rbind(
        data.frame(Sex='F',
                   TL=lengz,
                   Twt=AwT*lengz^BwT),
        data.frame(Sex='M',
                   TL=lengz,
                   Twt=AwT.M*lengz^BwT.M))%>%
        ggplot(aes(TL,Twt,color=Sex))+
        geom_line(size=1.25)+xlab("TL (cm)")+ylab('Twt (kg)')+
        theme_PA() 
      
      plt[[4]]=data.frame(Age=eig,
                          M=Mmean.mean.at.age)%>%
        ggplot(aes(Age,M))+
        geom_line(size=1.25,color="black")+xlab("Age")+ylab(expression(paste(plain("Natural mortality (year") ^ plain("-1"),")",sep="")))+
        theme_PA()
      
      plt[[5]]=data.frame(TL=lengz,
                          FL=(lengz-b_FL.to.TL)/a_FL.to.TL)%>%
        ggplot(aes(FL,TL))+
        geom_line(size=1.25)+xlab("FL (cm)")+ylab('TL (cm)')+
        theme_PA()+
        annotate("text", x = mean(lengz), y = mean(lengz)*.8,
                 label = paste('TL =',a_FL.to.TL,' FL + ',b_FL.to.TL ,sep=''))

      PLT=ggarrange(plotlist = plt, common.legend = TRUE)
      #annotate_figure(PLT, top = text_grob(capitalize(names(List.sp)[l]),face = "bold", size = 14))
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(List.sp[[l]]$Name),"/",AssessYr,"/1_Inputs/Visualise data/Biology.tiff",sep=''),
             width = 12,height = 8,compression = "lzw")
      
    })
   }
}


#export catch year index for reporting
fut.yr=as.numeric(substr(Last.yr.ktch,1,4))+years.futures
write.csv(data.frame(year.current=Last.yr.ktch,
                     year.future=paste(fut.yr,substr(fut.yr+1,3,4),sep='-')),
          paste(Rar.path,'year_current_future.csv',sep='/'),row.names = FALSE)


#---18. Catch-only assessments --------------------------------------
#Use Catch-only methods for species with only catch and life history data?
Species.with.length.comp=vector('list',N.sp)
for(i in 1:N.sp)
{
  if(any(grepl('Size_composition',names(Species.data[[i]]))))
  {
    d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[i]]))]
    if(length(d.list)>0)
    {
      if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
      if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
      Life.history=List.sp[[i]]
      for(s in 1:length(d.list))
      {
        d.list[[s]]=d.list[[s]]%>%
          filter(FL>=Life.history$Lzero)%>%
          mutate(fishry=ifelse(grepl("NSF.LONGLINE",names(d.list)[s]),'NSF',
                               ifelse(grepl("Survey",names(d.list)[s]),'Survey',
                                      ifelse(grepl("Other",names(d.list)[s]),'Other',
                                             'TDGDLF'))),
                 year=as.numeric(substr(FINYEAR,1,4)),
                 TL=FL*Life.history$a_FL.to.TL+Life.history$b_FL.to.TL,
                 size.class=TL.bins.cm*floor(TL/TL.bins.cm))%>%
          filter(TL<=with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL))))%>%
          dplyr::rename(sex=SEX)%>%
          filter(!is.na(sex))%>%
          group_by(year,fishry,sex,size.class)%>%
          summarise(n=n())%>%
          ungroup()
        
      }
      d.list <- d.list[!is.na(d.list)]
      d.list=do.call(rbind,d.list)%>%
        group_by(year,fishry,sex,size.class)%>%
        summarise(n=sum(n))%>%
        ungroup()%>%
        filter(!is.na(year))
      MAXX=max(d.list$size.class)
      Tab.si.kl=table(d.list$size.class)
      if(sum(Tab.si.kl[(length(Tab.si.kl)-1):length(Tab.si.kl)])/sum(Tab.si.kl)>.05) MAXX=MAXX*1.2
      MAXX=min(MAXX,30*ceiling(with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/30))
      size.classes=seq(min(d.list$size.class),MAXX,by=TL.bins.cm)
      missing.size.classes=size.classes[which(!size.classes%in%unique(d.list$size.class))]
      if(fill.in.zeros & length(missing.size.classes)>0)
      {
        d.list=rbind(d.list,d.list[1:length(missing.size.classes),]%>%
                       mutate(size.class=missing.size.classes,
                              n=0))
      }
      Table.n=d.list%>%group_by(year,fishry,sex)%>%
        summarise(N=sum(n))%>%
        mutate(Min.accepted.N=ifelse(!fishry=='Survey',Min.annual.obs.ktch,10))%>%
        filter(N>=Min.accepted.N)%>%
        mutate(dummy=paste(year,fishry,sex))
      
      if(nrow(Table.n)>0) Species.with.length.comp[[i]]=Keep.species[i]
      rm(Life.history)
    }
  }
}
Species.with.length.comp=c(do.call(rbind,Species.with.length.comp)) 
not.fitting.SS3=c("pigeye shark")  #no Hessian 
Species.with.length.comp=subset(Species.with.length.comp,!Species.with.length.comp%in%not.fitting.SS3)
write.csv(paste('1.',capitalize(Species.with.length.comp),sep=''),paste(Rar.path,'Species.with.length.comp.csv',sep='/'),row.names = F) 
if(Assessed.ktch.only.species=='All')
{
  Catch.only.species=Keep.species
  Catch.only.species_display=Catch.only.species
  if(display.only.catch.only.sp)
  {
    Catch.only.species_display=Keep.species[which(!Keep.species%in%names(compact(Catch.rate.series)))]
    Catch.only.species_display=subset(Catch.only.species_display,!Catch.only.species_display%in%Species.with.length.comp)
  }
}

if(Assessed.ktch.only.species=='Only.ktch.data') 
{
  Catch.only.species=Keep.species[which(!Keep.species%in%names(compact(Catch.rate.series)))]
  Catch.only.species=subset(Catch.only.species,!Catch.only.species%in%Species.with.length.comp)
  Catch.only.species_display=Catch.only.species
}

Catch.only.species_only.ktch=Keep.species[which(!Keep.species%in%names(compact(Catch.rate.series)))]
Catch.only.species_only.ktch=subset(Catch.only.species_only.ktch,!Catch.only.species_only.ktch%in%Species.with.length.comp)

if(Do.Ktch.only) fn.source1("Apply_catch only assessments.R")  #~60 mins per species (DBSRA, CMSY & SSS, 1 scenario)

#Ballpark M check for SSS
check.ballpark.M=FALSE
if(check.ballpark.M)
{
  a=1.46; b=-1.01   #fish
  a=0.941; b=-0.873  #cetaceans
  M.hoenig=function(a,b,Age.max) round(exp(a+b*log(Age.max)),3)
  Expected.M=vector('list',N.sp)
  names(Expected.M)=Keep.species
  for(i in 1:N.sp) Expected.M[[i]]=M.hoenig(a=0.941,b=-0.873,Age.max=List.sp[[i]]$Max.age.F[1])
  write.csv(unlist(Expected.M),handl_OneDrive('Analyses/Population dynamics/Expected.M.for.SSS_based on Hoenig.csv'))
  
}

#Bmsy/k by species   
plot.bmsyK=FALSE
if(plot.bmsyK)
{
  # Bmsy.K based on DBSRA
  BmsyK.DBSRA=vector('list',N.sp)
  for(n in 1:N.sp)
  {
    if(Keep.species[n]%in%Catch.only.species)
    {
      x=Catch_only$DBSRA$estimates[[n]]%>%
        filter(Scenario=='S1' & Parameter=='Bmsy/K')%>%dplyr::select(Median)%>%
        mutate(Species=names(Catch_only$DBSRA$estimates)[n])
      BmsyK.DBSRA[[n]]=x
    }
  }
  BmsyK.DBSRA=do.call(rbind,BmsyK.DBSRA)%>%
                  rename(BmsyK=Median)%>%
                  mutate(Method='DBSRA')%>%
                  left_join(BmsyK.species%>%dplyr::select(Species,r,Gen.time),by='Species')%>%
                  relocate(names(BmsyK.species))
  
    # Bmsy.K based on CMSY
  BmsyK.CMSY=vector('list',N.sp)
  for(n in 1:N.sp)
  {
    if(Keep.species[n]%in%Catch.only.species)
    {
      dummy=Catch_only$CMSY$estimates[[n]]%>%
        filter(Scenario=='S1' & Parameter%in%c('K','Bmsy'))%>%dplyr::select(Median)
      dummy=data.frame(Median=dummy[rownames(dummy)=='Bmsy','Median']/dummy[rownames(dummy)=='K','Median'])%>%
        mutate(Species=names(Catch_only$CMSY$estimates)[n])
      BmsyK.CMSY[[n]]=dummy
    }
  }
  BmsyK.CMSY=do.call(rbind,BmsyK.CMSY)%>%
              rename(BmsyK=Median)%>%
              mutate(Method='CMSY')%>%
              left_join(BmsyK.species%>%dplyr::select(Species,r,Gen.time),by='Species')%>%
              relocate(names(BmsyK.species))

    
    # Bmsy.K based on SSS
  BmsyK.SSS=vector('list',N.sp)
  for(n in 1:N.sp)
  {
    if(Keep.species[n]%in%Catch.only.species)
    {
      x=Catch_only$SSS$estimates[[n]]%>%
        filter(Scenario=='S1' & Par=='bmsy/k')%>%dplyr::select(Median)%>%
        mutate(Species=Keep.species[n])
      BmsyK.SSS[[n]]=x
    }
  }
  BmsyK.SSS=do.call(rbind,BmsyK.SSS)%>%
    rename(BmsyK=Median)%>%
    mutate(Method='SSS')%>%
    left_join(BmsyK.species%>%dplyr::select(Species,r,Gen.time),by='Species')%>%
    relocate(names(BmsyK.species))

  
  #Plot
  rbind(BmsyK.species,BmsyK.DBSRA,BmsyK.CMSY,BmsyK.SSS)%>%
    ggplot(aes(r,BmsyK,color=Gen.time, label=Species))+
    geom_point(shape=19,size=2)+
    facet_wrap(~Method)+
    scale_color_gradient(low = "yellow", high = "red")+
    ylim(0,1)+ylab("Bmsy/K")+
    theme(legend.position="bottom")+
    geom_text_repel(segment.colour='grey60',col='black',size=2,box.padding = 0.5)
  ggsave(handl_OneDrive('Analyses/Population dynamics/Bmsy.over.K_vs_r.tiff'), 
         width = 8,height = 6, dpi = 300, compression = "lzw")
  
}

#Get Consequence and likelihoods
Store.cons.Like_COM=list(Depletion=NULL,B.over.Bmsy=NULL)     
DD.depletion=DD.B.over.Bmsy=vector('list',length(Lista.sp.outputs))
for(l in 1:length(Lista.sp.outputs))  
{
  ddis=names(Lista.sp.outputs)[l]
  DD.B.over.Bmsy[[l]]=read.csv(paste(Rar.path,'/Table 4. Catch only_current.B.over.Bmsy_',ddis,'.csv',sep=''))%>%
                        gather(Species,Probability,-c(Range,finyear))%>%
                        mutate(Species=gsub(".", " ", Species, fixed=TRUE))
  
  DD.depletion[[l]]=read.csv(paste(Rar.path,'/Table 4. Catch only_current.depletion_',ddis,'.csv',sep=''))%>%
                        gather(Species,Probability,-c(Range,finyear))%>%
                        mutate(Species=gsub(".", " ", Species, fixed=TRUE))
}
Store.cons.Like_COM$Depletion=do.call(rbind,DD.depletion)
Store.cons.Like_COM$B.over.Bmsy=do.call(rbind,DD.B.over.Bmsy)
  

clear.log('Store.sens')
clear.log('output')
clear.log('out')
clear.log('apply.DBSRA')
clear.log('apply.OCOM')
clear.log('apply.CMSY')
clear.log('Catch_only')
clear.log('store.kobes')
clear.log('F.Fmsy')
clear.log('fn.agg.at.level.and.exprt')
clear.log('visualize.dat')
clear.log('Mod.AV')
clear.log('mod.average')
clear.log('dummy.store.Kobe.probs')

#---19. Spatio-temporal catch and effort. Reported TDGLF and NSF ----   
fn.source1('TDGLF and NSF spatio-temporal catch and effort.r')


#---20. Changes in observed mean length ----
if(Do.Changes.in.observed.size) fn.source1("Changes_mean_length.r")

#---21. Changes in reported mean weight of individuals caught in the TDGDLF  ----
#Is there a strong declining trend in mean weights? (Leitao 2019)
if(Do.Changes.in.reported.size) fn.source1("Changes_reported_size.r")
clear.log('Logbook')
clear.log('Logbook.sp')


#---22. Dynamic catch and size composition with dome-shaped selectivity --------------------------------------
#note: single-area, two-sexes, size-structured integrated model fitted to catch and size composition
#       Selectivity is assumed to be known.
#       Uncertainty derived from resampling variance-cov matrix  
#       Only applicable to non-indicator species with representative size comp samples.
#       For indicator species, use integrated model.
#       Not used due to poor contrast in available size composition and uncertain growth parameters
if(do.Dynamic.catch.size.comp) fn.source1("Dynamic.catch.size.comp.r")


#---23. Catch curve with dome-shaped selectivity --------------------------------------
#note: superseded by Dynamic catch and size composition with dome-shaped selectivity
if(do.Size.based.Catch.curve) fn.source1("Size.based.Catch.curve.r")
clear.log('fun.get.prop.zero.plus.wght')
clear.log('fn.plt.mn.ktch.wght')
clear.log('fun.get.prop.zero.plus.length')
clear.log('fun.change.mean.len')
clear.log('recruit.cpiui')


#---24. JABBA Surplus production model -------------------------------------------------
#note: Only fitting to 'indicator and 'other species' with representative abundance time series 
#Assumptions: negligible exploitation at start of time series

  #Get m parameter value
m.par=vector('list',N.sp)
names(m.par)=Keep.species
fn2 <- function(m)
{
  pred=m^(-(1/(m-1)))
  return((bmsy.over.K-pred)^2)
}
for(l in 1:N.sp)
{
  bmsy.over.K=List.sp[[l]]$Sens.test$JABBA$bmsyk[1]
  m=optim(1.5,fn2,method ="Brent",lower=0.1,upper=5)  
  m.par[[l]]=data.frame(Species=names(m.par)[l],bmsyK=bmsy.over.K,m=m$par)
}
m.par=do.call(rbind,m.par)

if(Do.StateSpaceSPM) fn.source1("Apply_StateSpaceSPM.R")   #takes 1.3 (8 species with 8 scenarios and fit diagnostics)

Store.cons.Like_JABBA=list(Depletion=NULL,B.over.Bmsy=NULL)     
DD.depletion=DD.B.over.Bmsy=vector('list',length(Lista.sp.outputs))
for(l in 1:length(Lista.sp.outputs))  
{
  ddis=names(Lista.sp.outputs)[l]
  DD.B.over.Bmsy[[l]]=read.csv(paste(Rar.path,'/Table 8. JABBA CPUE_Current.B.over.Bmsy_',ddis,'.csv',sep=''))%>%
    dplyr::select(-c(Scenario,Model))%>%
    gather(Species,Probability,-c(Range,finyear))%>%
    mutate(Species=gsub(".", " ", Species, fixed=TRUE))
  
  DD.depletion[[l]]=read.csv(paste(Rar.path,'/Table 8. JABBA CPUE_Current.depletion_',ddis,'.csv',sep=''))%>%
    dplyr::select(-c(Scenario,Model))%>%
    gather(Species,Probability,-c(Range,finyear))%>%
    mutate(Species=gsub(".", " ", Species, fixed=TRUE))
}
Store.cons.Like_JABBA$Depletion=do.call(rbind,DD.depletion)
Store.cons.Like_JABBA$B.over.Bmsy=do.call(rbind,DD.B.over.Bmsy)  


clear.log('Store.sens')
clear.log('output')
clear.log('F.Fmsy')
clear.log('apply.JABBA')
clear.log('out')
clear.log('fn.ktch.cpue')     
clear.log('Post')
clear.log('Post.density')
clear.log('Post.trace.plot')
clear.log('Post.running.mean')
clear.log('dummy.store.Kobe.probs')
clear.log('State.Space.SPM')

#---25. Integrated Stock Synthesis (Age-based) Model-------------------------------------------------
if(Do.integrated) fn.source1("Apply_SS.R")   #Takes ~ 10 hours

clear.log('Store.sens')
clear.log('output')
clear.log('F.Fmsy')
clear.log('out')
clear.log('fn.ktch.cpue')     
clear.log('Post')
clear.log('dummy.store.Kobe.probs')

# Create SS-DL tool input files
Export.SS_DL.tool.inputs=FALSE
if(Export.SS_DL.tool.inputs)
{
  for(i in 1:N.sp)
  {
    fn.create.SS_DL_tool_inputs(Life.history=List.sp[[i]],
                                CACH=KtCh,
                                this.wd=paste(handl_OneDrive("Workshops/2023_Jason Cope_SS_DL_tool/Myinputfiles/"),
                                              capitalize(Keep.species[i]),sep=''),
                                Neim=Keep.species[i],
                                InputData=Species.data[[i]],
                                KtchRates=Catch.rate.series[[i]])
    print(paste('Export SS-DL tool inputs for ------', Keep.species[i]))
  }
}


#---26. Integrated Bespoke (Size-based) Model-------------------------------------------------
if(Do.bespoked) fn.source1("Integrated_size_based.R")


#---27. Weight of Evidence Assessment ------

#1. Calculate risk for each line of evidence  
  #1.1. PSA 
#note: only use PSA scores for species dropped by the PSA as this is the only available LoE 
Drop.species=subset(Drop.species,!Drop.species=="other sharks")
Risk.PSA=vector('list',length(Drop.species))  
names(Risk.PSA)=Drop.species
for(s in 1:length(Risk.PSA))
{
  xxx=PSA.out%>%mutate(Species=tolower(Species))%>%filter(Species==tolower(names(Risk.PSA)[s]))
  if(xxx$Vulnerability=='Low') dummy=data.frame(Species=capitalize(xxx$Species),
                                                Consequence=1:4,
                                                Likelihood=c(2,1,1,1),
                                                w=0.5)  
  if(xxx$Vulnerability=='Medium') dummy=data.frame(Species=capitalize(xxx$Species),
                                                   Consequence=1:4,
                                                   Likelihood=c(4,2,1,1),
                                                   w=0.5)  
  Risk.PSA[[s]]=dummy
  rm(xxx,dummy)
}
Risk.PSA=do.call(rbind,Risk.PSA)

  #1.2. COMS  
if(exists("Store.cons.Like_COM"))
{
  if(Choose.probability=="Depletion")   Risk.COM=fn.risk(d=Store.cons.Like_COM$Depletion,w=0.5)
  if(Choose.probability=="B.over.Bmsy") Risk.COM=fn.risk(d=Store.cons.Like_COM$B.over.Bmsy,w=0.5)
}

  #1.3. JABBA
if(exists("Store.cons.Like_JABBA"))
{
  if(Choose.probability=="Depletion")   Risk.JABBA=fn.risk(d=Store.cons.Like_JABBA$Depletion,w=0.5)
  if(Choose.probability=="B.over.Bmsy") Risk.JABBA=fn.risk(d=Store.cons.Like_JABBA$B.over.Bmsy,w=0.5)
}

  #1.4. SS3
if(exists("Store.cons.Like_Age.based"))
{
  if(Choose.probability=="Depletion")   Risk.integrated=fn.risk(d=Store.cons.Like_Age.based$Depletion,w=0.5)
  if(Choose.probability=="B.over.Bmsy") Risk.integrated=fn.risk(d=Store.cons.Like_Age.based$B.over.Bmsy,w=0.5)
}


#2. Determine overall final risk   #ACA: how to combine the risk from COM, JABBA, SS, spatial catch stuff...
Risk_Drop.species=Risk.PSA
Risk_Other.sp=Risk.COM%>%
                filter(tolower(Species)%in%Lista.sp.outputs$Other.sp)
Risk_Indicator.sp=Risk.COM%>%
                  filter(tolower(Species)%in%Lista.sp.outputs$Indicator.sp)


#3. Plot Consequence and Likelihood and extract overall risk by species   
Store.risk_Drop.species=fn.risk.figure(d=Risk_Drop.species, Risk.colors=RiskColors, out.plot=FALSE)

Store.risk_Other.sp=fn.risk.figure(d=Risk_Other.sp, Risk.colors=RiskColors, out.plot=TRUE)
ggsave(paste(Rar.path,"Risk_Other.sp.tiff",sep='/'),width = 12,height = 10,compression = "lzw")

Store.risk_Indicator.sp=fn.risk.figure(d=Risk_Indicator.sp, Risk.colors=RiskColors, out.plot=TRUE)
ggsave(paste(Rar.path,"Risk_Indicator.sp.tiff",sep='/'),width = 10,height = 10,compression = "lzw")


#4. Display final risk for all species combined   
Final.risk_Drop.species=Store.risk_Drop.species%>%
                                group_by(Score,Risk)%>%
                                tally()%>%
                                mutate(Species=paste0("PSA-only species \n", "(n=",n,")"))%>%
                                dplyr::select(Species,Score,Risk)

Final.risk_Other.sp=Store.risk_Other.sp%>%dplyr::select(Species,Score,Risk)
Final.risk_Indicator.sp=Store.risk_Indicator.sp%>%dplyr::select(Species,Score,Risk)

p=rbind(Final.risk_Drop.species,Final.risk_Other.sp,Final.risk_Indicator.sp)%>%
  mutate(Species=factor(Species,levels=p%>%arrange(Score)%>%pull(Species)))

p%>%
  ggplot(aes(Species,Score,fill=Risk))+
  geom_bar(stat="identity")+ coord_flip()+
  theme_PA(axs.t.siz=14)+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')+
  ylab('Risk score')+xlab('')+
  scale_fill_manual(values=Risk.colors)+
  scale_y_continuous(breaks = breaks_width(1))

ggsave(paste(Rar.path,"Risk_all species together.tiff",sep='/'),width = 8,height = 10,compression = "lzw")


#---Remove this stuff once risk stuff is done  ------------------------------------------------- 


  #1.2. catch spatio-temporal    #ACA
stopped.landing=c('angel sharks','sawsharks',
                  "spurdogs","grey nurse shark")  #landings ceased/decreased due to lack of market for spurdogs, 
                                                  # and marine heatwave for angel and sawsharks (industry consultation AMM 2022)
                                                  #or protection for greynurse shark
Risk.spatial.temporal.ktch=vector('list',length(Store.spatial.temporal.ktch))
names(Risk.spatial.temporal.ktch)=names(Store.spatial.temporal.ktch)
for(s in 1: length(Risk.spatial.temporal.ktch))
{
  sp.dat=Store.spatial.temporal.ktch[[s]]$prop.ktch_over.fished.bocks%>%
    mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
    filter(year<=(Last.yr.ktch.numeric))%>%
    filter(year>=(Last.yr.ktch.numeric-4))%>%  #average last 5 years  
    group_by(SNAME)%>%
    summarise(Mean=mean(prop))%>%
    ungroup()%>%
    rename(Species=SNAME)%>%
    filter(!Species%in%stopped.landing)
  dummy=vector('list',nrow(sp.dat))
  names(dummy)=sp.dat$Species
  for(x in 1:length(dummy))
  {
    dd=sp.dat%>%filter(Species==names(dummy)[x])%>%pull(Mean)
    if(dd==1) d=data.frame(Max.Risk.Score=c(2,0,0,0))
    if(dd<1 & dd>=.75) d=data.frame(Max.Risk.Score=c(0,4,0,0))
    if(dd<.75 & dd>=.5) d=data.frame(Max.Risk.Score=c(0,0,6,0))
    if(dd<.5 & dd>=.25) d=data.frame(Max.Risk.Score=c(0,0,12,0))
    if(dd<.25) d=data.frame(Max.Risk.Score=c(0,0,0,16))
    dummy[[x]]=d
  }
  Risk.spatial.temporal.ktch[[s]]=dummy
}

#1.3. overall effort-management               
Risk.effort.mangmnt=Risk.spatial.temporal.ktch
Rel.eff.s=Effort.monthly
Rel.eff.n=rbind(Effort.monthly.north%>%rename(Total='Hook hours')%>%dplyr::select(FINYEAR,Total),  #add 0 effort if no fishing
                data.frame(FINYEAR=Effort.monthly$FINYEAR[which(!Effort.monthly$FINYEAR%in%Effort.monthly.north$FINYEAR)],
                           Total=0))%>%
                arrange(FINYEAR)
Yr.max.eff.s=Rel.eff.s%>%filter(Total==max(Total))%>%pull(FINYEAR)
Yr.max.eff.n=Rel.eff.n%>%filter(Total==max(Total))%>%pull(FINYEAR)
Rel.eff.n=mean(Rel.eff.n$Total[(nrow(Rel.eff.n)-4):nrow(Rel.eff.n)])/max(Rel.eff.n$Total) #average last 5 years
Rel.eff.s=mean(Rel.eff.s$Total[(nrow(Rel.eff.s)-4):nrow(Rel.eff.s)])/max(Rel.eff.s$Total)
REgn=Tot.ktch%>% 
        filter(Name%in%names(List.sp))%>%
        group_by(Name,Region)%>%
        summarise(Tot=sum(LIVEWT.c))
for(s in 1:length(Risk.effort.mangmnt))
{
  dat=Store.spatial.temporal.ktch[[s]]$Fished.blks.Fishery
  dat.max.n=dat%>%filter(FINYEAR==Yr.max.eff.n & Fishery=='Northern')%>%pull(Tot)
  dat.max.s=dat%>%filter(FINYEAR==Yr.max.eff.s & Fishery=='Southern')%>%pull(Tot)
  dat=dat%>%  
    filter(year<=(Last.yr.ktch.numeric))%>%
    filter(year>=(Last.yr.ktch.numeric-4))%>%  #average last 5 years  
    group_by(Fishery)%>%
    summarise(Mean=mean(Tot))%>%
    ungroup
  
  Rel.blk.n=dat%>%filter(Fishery=='Northern')%>%pull(Mean)/dat.max.n
  Rel.blk.s=dat%>%filter(Fishery=='Southern')%>%pull(Mean)/dat.max.s
  
  Eff.blk.n=Rel.eff.n*Rel.blk.n
  Eff.blk.s=Rel.eff.s*Rel.blk.s
  
  dis.sp=Lista.sp.outputs[-match('additional.sp',names(Lista.sp.outputs))][[s]]
  dummy=vector('list',length(dis.sp))
  names(dummy)=dis.sp
  for(x in 1:length(dummy))
  {
    Which.ef=REgn%>%filter(Name==names(dummy)[x])%>%pull(Region)
    Which.ef=paste(Which.ef,collapse='-')
    
    if(Which.ef=='North') Which.ef=Eff.blk.n
    if(Which.ef=='South') Which.ef=Eff.blk.s
    if(Which.ef=='North-South') Which.ef=max(c(Eff.blk.s,Eff.blk.n))
    
    if(Which.ef==0) d=data.frame(Max.Risk.Score=c(2,0,0,0))
    if(Which.ef>0 & Which.ef<=.25) d=data.frame(Max.Risk.Score=c(0,4,0,0))
    if(Which.ef>.25 & Which.ef<=.5) d=data.frame(Max.Risk.Score=c(0,0,6,0))
    if(Which.ef>.5 & Which.ef<=.75) d=data.frame(Max.Risk.Score=c(0,0,12,0))
    if(Which.ef>.75) d=data.frame(Max.Risk.Score=c(0,0,0,16))
    dummy[[x]]=d
  }
  Risk.effort.mangmnt[[s]]=dummy
}

#1.4. population dynamics models

Risk.tab=data.frame(Consequence=paste("C",1:4,sep=""),
                    L1=NA,L2=NA,L3=NA,L4=NA,
                    Max.Risk.Score=NA)

Risk.COM=vector('list',N.sp)
names(Risk.COM)=names(List.sp)
Risk.JABBA=Risk.COM
if(Do.integrated) Risk.integrated=Risk.COM
for(s in 1: N.sp)
{
  #COM
  if(exists("Store.cons.Like_COM"))
  {
    Risk.COM[[s]]=fn.risk(likelihood=Store.cons.Like_COM[[s]])
  }
  
  #JABBA
  if(exists('Store.cons.Like_JABBA'))
  {
    Risk.JABBA[[s]]=fn.risk(likelihood=Store.cons.Like_JABBA[[s]]) 
  }
  
  #Integrated model
  if(exists('Store.cons.Like_Age.based'))
  {
    Risk.integrated[[s]]=fn.risk(likelihood=Store.cons.Like_Age.based[[s]])
  }
       
  
}


#2. Integrate the risk from each line of evidence    

#note: Use a weighted sum to aggregate the Risk Categories form the alternative lines of evidence
#      Normalize each criterion by dividing by the highest value of each criterion. 
#      Assign weights to each criteria 
Order=c('Negligible','Low','Medium','High','Severe')
Order <- factor(Order,ordered = TRUE,levels = Order)


All.sp=sort(c(Drop.species,Keep.species))
Overall.risk=vector('list',length(All.sp))
names(Overall.risk)=All.sp
for(s in 1:length(All.sp))
{
  nm=All.sp[s]
  nm2=nm
  if(nm=="smooth hammerhead") nm2='hammerheads'
  iid=ifelse(nm%in%Other.species,"Other.sp",ifelse(nm%in%names(Indicator.species),"Indicator.sp","Drop.species"))  
  iid=match(iid,names(Risk.spatial.temporal.ktch))
  dummy=list(psa=Risk.PSA[[match(nm,names(Risk.PSA))]]$Max.Risk.Score,
             sptemp=Risk.spatial.temporal.ktch[[iid]][[match(nm2,names(Risk.spatial.temporal.ktch[[iid]]))]]$Max.Risk.Score,
             efman=Risk.effort.mangmnt[[iid]][[match(nm,names(Risk.effort.mangmnt[[iid]]))]]$Max.Risk.Score,
             COM=Risk.COM[[match(nm,names(Risk.COM))]]$Max.Risk.Score,
             JABBA=Risk.JABBA[[match(nm,names(Risk.JABBA))]]$Max.Risk.Score)
  if(Do.integrated)dummy$integrated=Risk.integrated[[match(nm,names(Risk.integrated))]]$Max.Risk.Score
  dummy=dummy[!sapply(dummy, is.null)]
  NMs=names(dummy)
  dummy=do.call(cbind,dummy)
  colnames(dummy)=NMs
  rownames(dummy)=paste("C",1:4,sep='')
  
  Overall.risk[[s]]=Integrate.LoE(Cons.Like.tab=dummy,
                                  criteriaMinMax <- "max",
                                  plot.data=FALSE,
                                  LoE.weights = LoE.Weights[match(colnames(dummy),names(LoE.Weights))],
                                  Normalised="YES")
}


#3. Display overall risk for each species
LoE.col=c(psa="coral3", sptemp="grey50", efman="burlywood4",
          COM="grey80", JABBA="grey35", integrated="black")
LoE.names=c(psa="PSA", sptemp="Blocks fished", efman="Effort & Manag.",
            COM="COM", JABBA="BSSPM", integrated="Int. mod.")
if(!Do.integrated)
{
  LoE.col=LoE.col[-match('integrated',names(LoE.col))]
  LoE.names=LoE.names[-match('integrated',names(LoE.names))]
}
All.Sp.risk.ranking=factor(unlist(lapply(Overall.risk, '[[', 'risk')),levels=levels(Order))  
All.Sp.risk.ranking=names(sort(All.Sp.risk.ranking))
All.Sp.risk.ranking=All.Sp.risk.ranking[!duplicated(All.Sp.risk.ranking)]


#3.1. Other species
Sp.risk.ranking=subset(All.Sp.risk.ranking,!All.Sp.risk.ranking%in%names(Indicator.species))
fn.fig(paste(Rar.path,"Risk_Other.sp",sep='/'), 2000, 2400)
par(mar=c(.5,.5,3,1),oma=c(3,13.5,.5,.1),las=1,mgp=c(1,.5,0),cex.axis=1.5,xpd=TRUE)
layout(matrix(c(rep(1,6),rep(2,3)),ncol=3))

  #3.1.1. Risk for each line of evidence
fn.each.LoE.risk(N.sp=length(Sp.risk.ranking),CX.axis=0.95) #main plot
for(s in 1:length(Sp.risk.ranking))
{
  nm=Sp.risk.ranking[s]
  dummy=Overall.risk[[match(nm, names(Overall.risk))]]$Cons.Like.tab
  if(ncol(dummy)==1) Adjst=0
  if(ncol(dummy)>1) Adjst=seq(-.3,.3,length.out = ncol(dummy))

  for(nn in 1:ncol(dummy))
  {
    segments(0,s+Adjst[nn],max(dummy[,nn]),s+Adjst[nn],
             lwd=3.75,lend=1,col=LoE.col[match(colnames(dummy)[nn],names(LoE.col))])
  }

}
axis(2,1:length(Sp.risk.ranking),capitalize(Sp.risk.ranking),cex.axis=.8)
mtext("Risk score",1,cex=1.25,line=2)
legend(-1.5,length(Sp.risk.ranking)+6,LoE.names[1:3],
       bty='n',col=LoE.col[1:3],lty=1,lwd=3,horiz = T,cex=1.20,
       text.width=c(0,1,2.25))
legend(-1.5,length(Sp.risk.ranking)+4.5,LoE.names[4:length(LoE.names)],bty='n',cex=1.20,
       col=LoE.col[4:length(LoE.col)],lty=1,lwd=3,horiz = T,text.width=c(0,1.5))

  #3.1.2. Overall risk
plot(0:1,ylim=c(0,length(Sp.risk.ranking)+1),fg='white',xaxs="i",yaxs="i",
     col="transparent",ylab="",xlab="",xaxt='n',yaxt='n')
for(s in 1:length(Sp.risk.ranking))
{
  ss=Sp.risk.ranking[s] 
  fn.overall.risk(N=s,
                  RISK=Overall.risk[[match(ss,names(Overall.risk))]]$risk,
                  sp=capitalize(Sp.risk.ranking[s]),
                  CX=.9)
}
axis(1,1.5,"Overall risk",cex.axis=1.75,col.ticks="white",padj=.5)

dev.off()


#3.2.Indicator species    
Sp.risk.ranking=subset(All.Sp.risk.ranking,All.Sp.risk.ranking%in%names(Indicator.species))
fn.fig(paste(Rar.path,"Risk_Indicator.sp",sep='/'), 2200, 2400)  

par(mar=c(.5,.5,3,1),oma=c(3,9.5,.5,.1),las=1,mgp=c(1,.5,0),cex.axis=1.5,xpd=TRUE)
layout(matrix(c(rep(1,6),rep(2,3)),ncol=3))

  #3.2.1. Risk for each line of evidence
fn.each.LoE.risk(N.sp=length(Sp.risk.ranking),CX.axis=1.35) #main plot
for(s in 1:length(Sp.risk.ranking))
{
  nm=Sp.risk.ranking[s]
  dummy=Overall.risk[[match(nm, names(Overall.risk))]]$Cons.Like.tab
  if(ncol(dummy)==1) Adjst=0
  if(ncol(dummy)>1) Adjst=seq(-.3,.3,length.out = ncol(dummy))
  
  for(nn in 1:ncol(dummy))
  {
    segments(0,s+Adjst[nn],max(dummy[,nn]),s+Adjst[nn],
             lwd=3.75,lend=1,col=LoE.col[match(colnames(dummy)[nn],names(LoE.col))])
  }
  
}
axis(2,1:length(Sp.risk.ranking),capitalize(Sp.risk.ranking),cex.axis=1.5)
mtext("Risk score",1,cex=1.25,line=2)
legend(-1.5,length(Sp.risk.ranking)+1.3,LoE.names[1:3],
       bty='n',col=LoE.col[1:3],lty=1,lwd=3,horiz = T,cex=1.25,
       text.width=c(0,1,2.25))
legend(-1.5,length(Sp.risk.ranking)+1.22,LoE.names[4:length(LoE.names)],bty='n',cex=1.5,
       col=LoE.col[4:length(LoE.col)],lty=1,lwd=3,horiz = T,text.width=c(0,1.5))

  #3.2.2. Overall risk
plot(0:1,ylim=c(0,length(Sp.risk.ranking)+1),fg='white',xaxs="i",yaxs="i",
     col="transparent",ylab="",xlab="",xaxt='n',yaxt='n')
for(s in 1:length(Sp.risk.ranking))
{
  ss=Sp.risk.ranking[s] 
  fn.overall.risk(N=s,
                  RISK=Overall.risk[[match(ss,names(Overall.risk))]]$risk,
                  sp=capitalize(Sp.risk.ranking[s]),
                  CX=1.5)
}
axis(1,1.5,"Overall risk",cex.axis=1.75,col.ticks="white",padj=.5)

dev.off()



#---28. Outputs for strategic papers  ------------------------------------------------- 
# Sawfish stock assessment  paper
fn.source1("Outputs_for_strategic_papers.r")
if(do.sawfish) fn.do.sawfish()
