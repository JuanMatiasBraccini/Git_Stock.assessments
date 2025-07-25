# ------ Script for running stock assessments for WA sharks---- ###################

#Notes:
#     1. '4. Productivity Susceptibility Analyses' filters out species based on '9. Catch criteria for selecting what species to assess quantitatively'
#     2. Assumption: If catches have never been >1% carrying capacity, then it's in unexploited 
#                    status so catch series have no information on productivity
#     3. Total (reconstructed) catches are used (commercial, recreational and TDGDLF discards);
#             note that 'recons_NT_catch.csv' only includes dusky and sandbar
#     4. A range of assessment methods are used depending on data availability: 
#               . Integrated size- and sex- structured models,
#               . State-space SPM,
#               . Catch curves,
#               . Catch-only
#    5. Tagging data: 
#             . Conventional: use '_Con_tag_BLK.rec' '_Con_tag_BLK.rel' '_Con_tag_Zn.rec' '_Con_tag_Zn.rel'
#             . Acoustic: 'Acous.Tag_Rep.Recap' has recaptures (treat as conventional for F estimation!!)
#                          '_Acous.Tag_Zn.rec_Acous.Tag'  (time step of 1 week!!)

#Steps: 
#     1. For each new assessment, update 'New.assessment', 'Year.of.assessment' and 'Last.yr.ktch' in '1. DEFINE GLOBALS'  
#     2. Define arguments (inputs) used in each of the shark species/species-complex assessed.
#     3. Bring in updated available data and parameters (see 'Matias\Reports\Steps for creating RAR and SFRAR.docx')
#     4. Determine which species to assess based on PSA
#     5. Run relevant population models according to data availability
#     6. Generate relevant outputs. Outputs are stored in each species folder ('.../Analyses/Population dynamics/...')
#         and in '.../Reports/RARs...'
#     7. Run '.../Analyses/Build documents/Weight of Evidence/WoE.Rmd' with updated '.../Report/RARs/.../Risk tables' folders to
#         generate WoE report.


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
library(grid)
library("gg3D")
library(ggridges)
library(janitor)
library(units)
library(ggforce)

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

Send.email.to="matias.braccini@dpird.wa.gov.au"   #send email when model run finalised
#Send.email.to="braccinimatias@gmail.com.au"  #IT firewall doesn't allow sending to gmail :(

#---1. DEFINE GLOBALS----- 

#1.New assessment
New.assessment="NO"
#New.assessment="YES"   #set to 'YES' if a new assessment is run for the first time

#2. Control what to implement
# turn on/off assessment method as appropriate

  #2.1 Level 1
Do.Ktch.only=FALSE
Simplified.scenarios=TRUE  #simplified version of COM sensitivity tests
do.Ktch.only.pin=do.JABBA.pin=TRUE

  #2.2 Other lines of evidence
Do.Spatio.temporal.catch.effort=FALSE
Do.Changes.in.observed.size=FALSE
Do.Changes.in.reported.size=FALSE

  #2.3 Level 3
do.Size.based.Catch.curve=FALSE  
do.Dynamic.catch.size.comp=FALSE #superseded by length-only SS3 (also, not yet implemented properly)

  #2.4 Level 4
Do.StateSpaceSPM=FALSE

  #2.5 Level 5
Do.integrated=FALSE  #integrated SS3
Do.bespoke=FALSE   #bespoke integrated size-based
Do.sim.Test="NO" #do simulation testing bespoke size-based model
do.Andre.model=FALSE

  #2.6 Outputs
do.F.series=FALSE   #output Fishing mortality time series
do.B.over.Bmsy.series=TRUE
do.F.over.Fmsy.series=TRUE

#3. Assessment dimensions
Year.of.assessment=AssessYr=2022
Last.yr.ktch="2021-22"
Last.yr.ktch.numeric=as.numeric(substr(Last.yr.ktch,1,4))

  #3.1 Future projections
future.models=c('State.Space.SPM','SS')   #c('Catch_only','State.Space.SPM','SS')
years.futures=5  #number of years to project
n.last.catch.yrs=5 #number of recent years used to calculate future catch
catches.futures="constant.last.n.yrs"
#catches.futures='upper.limit.catch.range'  #catch ranges are only available for indicator species
future.color="brown4"

#4. Model run
if(New.assessment=="YES") First.run="YES"  else #create all model input data sets and data presentation for new assessment
  First.run="NO"

#5. Add additional species of interest not selected by PSA given low catch trajectories but needed for specific assessment
additional.sp=NULL  #if no additional species assessment required
#if(Year.of.assessment==2022) additional.sp=c('green sawfish','narrow sawfish')   # 2022 sawfish assessment; 
                                    # dwarf and freshwater sawfish not assessed; reconstructed
                                    # catches do not consider TO catch or beach rec fishing

#6. Define if calculating r & steepness
if(New.assessment=="YES") do.r.prior=TRUE  else 
                          do.r.prior=FALSE
do.steepness=do.r.prior

#7. Define if exporting figures as jpeg or tiff (RAR requires jpeg)
Do.tiff="YES" 
Do.jpeg="NO"

#8. Catch units
KTCH.UNITS="TONNES" 
#KTCH.UNITS="KGS"    
if(KTCH.UNITS=="KGS") unitS=1000
if(KTCH.UNITS=="TONNES") unitS=1

#9. Catch criteria for selecting what species to assess quantitatively
Min.yrs=5
if(KTCH.UNITS=="KGS") Min.ktch=5000 
if(KTCH.UNITS=="TONNES") Min.ktch=5

#10. CPUEs
Min.cpue.yrs=5 #minimum number of years in abundance index
drop.large.CVs=FALSE  #drop observations with CV larger than MAX.CV or not. Superseded by Francis CVs 

  #10.1 Define species for which cpue is not considered to be indexing abundance:
survey_not.representative=c("scalloped hammerhead","great hammerhead",
                            "lemon shark","pigeye shark","dusky shark") #few individuals caught (<5 per trip) and huge CVs 
NSF_not.representative=c("scalloped hammerhead","great hammerhead",   #NSF cpue not used as unlikely to be representative
                          "lemon shark","pigeye shark","tiger shark",
                         "dusky shark","sandbar shark")
tdgdlf_not.representative=c("smooth hammerhead","spinner shark","dusky shark")   #catch rates are for 'hammerheads' and for both species cpue tracks catch so no depletion signal
tdgdlf_monthly_not.representative=c("sandbar shark")   #increasing cpue with increasing catch and very jumpy index 
other_not.representative=c("green sawfish","narrow sawfish") #Pilbara trawl cpue, rare event & not within species distribution core
drop.daily.cpue='2007&2008'  #drop from TDGDLF daily cpue (consistently higher cpues across species due to likely effort reporting bias)

#11. Size composition
MN.SZE=0    # initial bin size
#MN.SZE="size.at.birth"
TL.bins.cm=5  # size bin
Min.obs=10  #keep records with at least 10 observations
Min.shts=5  #keep records from at least 5 shots
Min.annual.obs.ktch=150 #Minimum number of annual observations for using length composition data
Min.annual.obs.ktch_survey=30
prop.min.N.accepted_other=0.5
Min.Nsamp=10   #Minimum number of trips for catch mean weight or length composition
fill.in.zeros=TRUE  #add missing length classes with all 0s

#12. Proportion of vessels discarding eagle rays in last 5 years (extracted from catch and effort returns)
prop.disc.ER=.4  

#13. PSA criteria
PSA.min.tons=5
PSA.min.years=Min.yrs
PSA.max.ton=50
Low.risk=2.64  #risk thresholds from Hobday et al 2007 & Micheli et al 2014
medium.risk=3.18

#14. Assumed PCM for reconstructed discards in TDGLDF
TDGLDF.disc.assumed.PCM="BaseCase" 
#TDGLDF.disc.assumed.PCM="100%" 

#15. Demography
First.Age=1        #start at 1 rather than 0 as M estimators yield unrealistically high M for Age 0
Max.Age.up.Scaler=1.3  #scaler of Max maximum age for species with only one record. 1.3 is mean across those species with a range
Max.r.value=.55 # Max r = 0.44 for blue shark (Cortes 2016); 0.51 for Scyliorhinus canicula (Cortes 2002)
Min.r.value=.025
reset.max.Age=FALSE  #set max Age to mean of max.age and max age.max
fill.NA.Max.Age.Max=TRUE  #fill in NA max.Age.max with Max.Age.up.Scaler 
externally.increase.M=FALSE   #increase M in pin_file to all SS to converge
Dusky.Sedar=mean(c(0.25,0.35))  #SEDAR 21, they also obtained too high h estimates (page 30)
Sandbar.Sedar=mean(c(0.25,0.4)) #SEDAR 21
ScallopedHH.Sedar=mean(c(0.69,0.71,0.67))  #SEDAR 77
SmoothHH.Sedar=0.78 #SEDAR 77
GreatHH.Sedar=0.71 #SEDAR 77
Mako.ICCAT=0.345   #ICCAT 2019
species.too.high.M1=NULL #c("gummy shark","whiskery shark")  

  #15.1 publish demographic parameters
Demo.published.values=data.frame(Species=c("angel sharks","copper shark","grey nurse shark","gummy shark","lemon shark",
                                  "spinner shark","spurdogs","tiger shark","whiskery shark","spiny dogfish",'blue shark',
                                  "dusky shark","great hammerhead","sandbar shark","scalloped hammerhead",
                                  "shortfin mako","smooth hammerhead"),
                        r=c(0.11,0.075,0.05,0.3,0.09,0.101,0.062,0.196,0.136,0.076,0.44,
                            0.055,0.146,0.085,0.104,0.063,0.182),
                        h=c(rep(NA,11),Dusky.Sedar,GreatHH.Sedar,Sandbar.Sedar,ScallopedHH.Sedar,Mako.ICCAT,SmoothHH.Sedar),
                        Reference=c(rep('Cortes 2016',11),'SEDAR 21','SEDAR 77','SEDAR 21','SEDAR 77','ICCAT 2019','SEDAR 77'))
test.Sedar=TRUE     #consider Dusky and Sandbar h estimates used in SEDAR in the modelled scenarios
Use.SEDAR.M=FALSE   #Set to TRUE if using SEDAR M @ age for dusky and sandbar

#16. Stock recruitment
Max.h.shark=.8   #mean of h for blue shark (ICCAT 2023 assessment; Cortes 2016, Kai & Fujinami 2018).
Min.h.shark=.3  #He et al 2006, Jason Cope pers comm
Max.SR_sigmaR.shark=0.3   #maximum recruitment variability (blue shark ICCAT 2023 0.29 for North, 0.5 for South, 0.3 bigskate; 0.2 dogfish; 0.18 sandbar)
do.random.h=TRUE  #take a random sample of h and M for SS or use empirical distributions

#17. Reference points
#note: Historically, there was a single unspecified reference point (40% unexploited biomass)
#      Currently, 0.4 is used as threshold, not BMSY
Biomass.threshold='Bmsy'  #MSC sets threshold to Bmsy and limit to 0.5 Bmsy (Clinton Syers)
Biomass.threshold.min=0.4  #Andre advices against using BMSY estimates from integrated models and uses 0.4 as proxy
Tar.prop.bmsny=1.2    # Target and Limit proportions of 'Biomass.threshold' 
Lim.prop.bmsy=0.5    #    source: Haddon et al 2014. 'Technical Reviews of Formal Harvest Strategies'.
#Fmsy.emp=function(M) 0.41*M     #Zhou et al 2012 but see Cortes & Brooks 2018
#SPR.thre=0.3   #Punt 2000 Extinction of marine renewable resources: a demographic analysis. 
#SPR.tar=0.4    #   Population Ecology 42, 
Minimum.acceptable.blim=0.2  #Blim should be the max(Minimum.acceptable.blim,Lim.prop.bmsy*Biomass.threshold)
col.Target='forestgreen'
col.Threshold='orange'
col.Limit='red'
COM.limit=1  #as proportion of MSY
COM.threshold=0.8
COM.target=0.5


#18. Catch-only Models
COM_use.this.for.risk='catch' #  define if using 'catch' or 'biomass' trajectories to determine risk from catch-only methods
do.ensemble.simulations=FALSE
catch.only=c('DBSRA','CMSY','SSS')      # define model types used
do.OCOM=FALSE
do.Catch.JABBA=FALSE   #redundant, same results as CMSY
if(do.OCOM) catch.only=c(catch.only,'OCOM')
if(do.Catch.JABBA) catch.only=c(catch.only,'JABBA')
CMSY.method="Haddon"    # Select which CMSY method to use. Haddon's datalowSA
#CMSY.method="Froese"  # Froese et al 2017 does not converge for dwarf or freshwater
do.parallel.SSS=TRUE   #do SSS in parallel or not (set to FALSE)
SSS.sims=5e2   #Cope 2013 did 1e3, same results with 5e2 but faster. 
ensims<-1e4   #DBSRA. 1e4 Dick & MacCall 2011
ensims.CSMY=2e4   #CMSY simulations
ensims.JABBA=3e4  # 3e4 Winkner et al 2019
Proc.Error=1e-3   #Catch-only default process error. Catch only per se yields highly uncertain estimates
Proc.Error.1=1e-2
SSS_criteria.delta.fin.dep=0.01
K_min=2000 #Min value of max K range in tonnes (based on overall catch ranges and K estimates)
r.prob.max=0.9999   #quantile probs for defining r range for CMSY
r.prob.min=1e-3
#Ensemble.weight='weighted'  #define if doing weighted or unweighted model average
Ensemble.weight='equal'
Wei.SSS=2   #give more weight to SSS given biological realism
tweak.BmsyK.Cortes=FALSE #update value to increase COM acceptance rate
tweak.Final.Bio_low=FALSE #update value to increase COM acceptance rate
n.last.catch.yrs_MSY.catch.only=5
modify.r_k_bounds.for.convergence=TRUE   #modify the range or f and k used in CMSY for some species to allow enough combos

  #18.1 Define which species to assess using catch-only methods
Assessed.ktch.only.species='All'                #assess all species with catch only methods (level 1 assessment)
#Assessed.ktch.only.species='Only.ktch.data'   #assess species with only catch data
display.only.catch.only.sp=FALSE               #just display catch and MSY for species with only catch

#19. Catch curve and YPR
  #19.1 Catch curve inputs
Init.F.Ktch.cur=0.2
do.eye.ball=FALSE  #eye-ball CVSizeAtAge 
CVSizeAtAge = c(0.03,0.03)  #this CV is not the CV used in SS3, it's the CV of the size transition matrix so cannot be too large  
Main.zone.mesh=data.frame(Species=c("dusky shark","gummy shark",#revise when new length comp available
                                    "sandbar shark","smooth hammerhead",
                                    "whiskery shark","spinner shark"))%>%       
        mutate(Zone=case_when(Species%in%c("dusky shark","gummy shark",
                                           "sandbar shark","smooth hammerhead",
                                           "whiskery shark")~'Zone1',
                              Species%in%c("spinner shark")~'West',
                              TRUE~''),
               Mesh=case_when(Species%in%c("dusky shark","gummy shark",
                                           "sandbar shark","smooth hammerhead",
                                           "whiskery shark")~6.5,
                              Species%in%c("spinner shark")~7,
                              TRUE~NA))
TimeStep = 1 # model timestep (e.g. 1 = annual (typical long lived species), 1/12 = monthly)
LenInc = TL.bins.cm  # TL in cm   
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelParams = c(300, 50) # L50, L95-L50 for gear selectivity            NOT USED
RetenParams = c(NA, NA) # L50, L95-L50 for retention                  NOT USED
DiscMort = 0 # proportion of fish that die due to natural mortality   NOT USED
DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute 
InitDelta = 50
RefnceAges = NA
standardise.mesh.zone=TRUE  #only use records from same mesh and zone (the most representative)
if(standardise.mesh.zone) used.selectivity='Empirical for selected mesh'  #use published sel for chosen mesh
if(!standardise.mesh.zone) used.selectivity='Estimated by SS'             #use combined mesh sel estimated in SS

    #19.1.2 Species for which 'other' fleet is main fleet but no length comp available 
Other.to.NSF=c("milk shark","pigeye shark","tiger shark","scalloped hammerhead")  
Other.to.TDGDLF=c('sawsharks')

  #19.2 YPR inputs
Dummy.F.mort=1e-4   #if Catch curve estimates F.mort at 0, then reset at low value to run YPR
WLrel_Type <- 1 # 1=power, 2=log-log relationship
ReprodScale <- 1 # 1=default (standard calculations for spawning biomass), 2=hyperallometric reproductive scaling with female mass (i.e. BOFFF effects)
ReprodPattern <- 1 # 1 = gonochoristic (separate sexes), 2 = protogynous (female to male sex change), 3 = protandrous (male to female sex change)
InitRatioFem <- 0.5 # Ratio of females to males at recruitment age
FinalSex_Pmax <- NA # Logistic sex change relationship parameters (max probability of final sex)
FinalSex_L50 <- NA # Logistic sex change relationship parameters (inflection point)
FinalSex_L95 <- NA # Logistic sex change relationship parameters (95% of max probability)
SRrel_Type <- 1 # 1 = Beverton-Holt, 2=Ricker
RefPointPlotOpt <- 2 # 0=don't plot, 1=plot defaults, 2=plot BMSY ref points
nReps = 100 # 200 or 500


#20. State space Surplus Production Models
state.space.SPM='JABBA'     # define model types used
do.parallel.JABBA=TRUE      # do JABBA in parallel or not (set to FALSE)
use.auxiliary.effort=FALSE  # using effort as auxiliary data yielded very high RMSE and poor fits to effort
Rdist = "lnorm"
KDIST="lnorm"  
PsiDist='beta'
Whiskery.q.periods=2 # split monthly cpue into this number of periods (Simpfendorfer 2000, Braccini et al 2021)
drop.intermediate.yrs=FALSE  #remove whiskery inermediate Q years or not
Gummy.q.periods=1     
Obs.Err.JABBA=0.01   #JABBA uses SE2" = CPUE.se^2 + fixed.obsE^2 
increase.CV.JABBA=TRUE   #for consistency with SS3 and because not using fixed.obsE
do.MCMC.diagnostics=do.hindcasting=FALSE
if(First.run=="YES")  do.MCMC.diagnostics=do.hindcasting=TRUE
Proc.Error.cpue=1e-01   # Default process error when fitting cpue to allow enough flexibilty (0.2 showed no differences) 
Proc.Error.cpue2=5e-02   #alternative, 5e-02 process error for JABBA  (Winker et al 2018 School shark) ; Beth Babcock suggested 0.01
k.cv=2                #Carrying capacity CV (Winker et al 2018 School shark)
PEELS=5               #number of years to peel for hindcasting and retrospective   
evaluate.07.08.cpue=FALSE  #run scenario with 2007 & 08 TDGDLF cpue

#21. Integrated age-based model 
Integrated.age.based='SS'   # define model types used
do.parallel.SS=TRUE         #do SS in parallel or not
do.all.sensitivity.tests=TRUE #set to FALSE as per required
SS3.run='final' #'test'     # switch to 'final' when model fitting is finalised to estimate uncertainty (Hessian, MCMC, etc)
create.SS.inputs=FALSE       #set to FALSE once happy with SS input files and only need to run the model
run_SS_plots=FALSE          #set to TRUE once happy with model and want to plot outputs
if(SS3.run=='final') run_SS_plots=TRUE
Calculate.ramp.years=FALSE  #switch to TRUE if new year of size composition available
Run.SS=FALSE                 #switch to TRUE if want to run parameter estimation
do.Cond.age.len.SS.format=FALSE   #use age-length data to estimate growth
                                  # this is not used as age-length sandbar and dusky is for GN and LL and 
                                  # for all 4 species observations were collected over multiple years
Mean.Size.at.age.species=NULL   #  Mean.Size.at.age.species=c("gummy shark","whiskery shark" )
                                # Not implemented. Wrong SS3 format. May be applicable to gummy and 
                                # whiskery (only for these species length-@-age data collected from gillnet fishery)
SS3_fleet.size.comp.used=c("Size_composition_West","Size_composition_Zone1","Size_composition_Zone2",
                           "Size_composition_NSF.LONGLINE","Size_composition_Survey",
                           "Size_composition_Other")
combine_NSF_Survey=NULL   #combine length composition data to estimate logistic selectivity
# combine_NSF_Survey=c("dusky shark","great hammerhead","lemon shark","milk shark",
#                      "pigeye shark","sandbar shark","scalloped hammerhead","tiger shark")
combine.sexes.tdgdlf=NULL 
combine.sexes.tdgdlf.daily=NULL 
combine.sexes.survey=c("dusky shark")
combine.sexes.nsf=NULL
combine.sexes=c(combine.sexes.survey,combine.sexes.tdgdlf,combine.sexes.tdgdlf.daily,
                "angel sharks","lemon shark","milk shark","scalloped hammerhead","tiger shark")
combine.sex_type=3  #0 0 means combined male and female ; 3 3 means data from both sexes will be used and they are scaled so that they together sum to 1.0; i.e., sex ratio is preserved
#fit.to.mean.weight.Southern2=NULL
fit.to.mean.weight.Southern2=c("spinner shark","whiskery shark")  #get model to fit mean weight regardless of available length comp
drop.len.comp.like=NULL    
survey.like.weight=NULL  #"dusky shark"  
use.Gab.trawl=TRUE   #note that this has only 1 year of data
add.gummy.gab=FALSE
species.increase.terminal.age=c("gummy shark","whiskery shark") #add a few extra age classes
species.constant.fec=c("whiskery shark")   #dodgy linear fec relationship
Extract.SS.parameters=FALSE  #extract SS3 selectivity pars
rescaled.species.sel=sort(c('great hammerhead','scalloped hammerhead','grey nurse shark','milk shark',
                       'shortfin mako','sawsharks','spinner shark',
                       'tiger shark','wobbegongs'))   #this species have no species-specific empirical sel (family was used)
Plus.gp.size=1.25  #add 25% to max size make sure no accumulation of survivals in last size class

alternative.sigmaR=NULL  #Sensitivity for sigmaR (effect on rec_devs)
#alternative.sigmaR="sandbar shark"
alternative.do_recdev=NULL  #Sensitivity for do_recdev method (effect on rec_devs)
#alternative.do_recdev="sandbar shark"

alternative.forecasting=NULL  #forecasting F rather than catch
#alternative.forecasting="sandbar shark"
# F.forecasting.values=list("sandbar shark"=c('Northern.shark'=1.83343e-04,  #get from Report (EXPLOITATION report:14)
#                                             'Other'=2.11086e-02,
#                                             'Southern.shark_1'=0,
#                                             'Southern.shark_2'=1.04784e-02))  
alternative.like.weigthing=NULL  #test alternative lambdas for survey and length comps
spatial.model="sandbar shark"     #NULL, build areas-as-fleets model

  #21.1 Set WRL as a separate fleet for these species
WRL.species=c("copper shark","dusky shark","shortfin mako",
              "smooth hammerhead","spinner shark","tiger shark") 

  #21.2 Sensitivity for NSF logistic selectivity for these species
alternative.NSF.selectivity=NULL
#alternative.NSF.selectivity=c("dusky shark","great hammerhead","lemon shark","pigeye shark",
#                              "sandbar shark","scalloped hammerhead","smooth hammerhead","tiger shark")

  #21.3 No empirical Selectivity for main fleet or length comp sample size is too small
#  so cannot implement any length-based assessment (Catch curve, SS3, etc)
no.empirical.sel.main.fleet=c(GAB.main="angel sharks",
                              SA_MSF.main="copper shark",
                              Indo.main_low.n.NSF="great hammerhead",
                              NSF.main_low.n.NSF="lemon shark",
                              NSF.main_low.n.NSF="pigeye shark",
                              GAB.main="sawsharks",
                              Multispecies.no.sel="spurdogs",
                              Taiwan.main_low.n.NSF="scalloped hammerhead",
                              WRL.main="tiger shark",
                              Multispecies.no.sel="wobbegongs")
resample.h.greynurse=FALSE  #no need to resample h

  #21.4 SS model run arguments
if(SS3.run=='final') Arg=''
if(SS3.run=='test') Arg= '-nohess'   #no Hessian 
Find_Init_LnRo=FALSE   #set to TRUE first time fitting model to find Init LnRo value so that Virgin Total biomass ~ K from JABBA  
SS3.q.analit.solu=TRUE   #calculate q analytically to save up pars, set to FALSE if using block Q (time changing Q)
block.species_Q=c("whiskery shark") #"gummy shark"
Extra_Q_species=c("spinner shark","tiger shark") #needed to allow fit. Not used
do.MC.multi=FALSE #doesn't work if estimating rec devs as rec devs are not updated with random sample
nMCsims=200  #number of Monte Carlo simulations for multivaritenormal
Arg.no.estimation='-maxfn 0 -phase 50 -nohess'  #no estimation. Used for Monte Carlo simulations
#MCMCsims=1e5; Thin=10; burning=1:(5*length(seq(1,MCMCsims,by=Thin))/100)   #5%  burning
#Arg=paste(' -mcmc',MCMCsims,' -mcsave', 100)  #MCMC

  #21.5 Assumed error distribution for abundance series
Abundance.error.dist='Lognormal'  #'Lognormal' if stand. cpue in normal space and CVs; 'Normal'

  #21.6 Default CV for time series with very small CVs
CV.use='loess'  #Francis 2011
#CV.use='fixed'
default.CV=0.15  #0.15 used by Punt 2009 gummy model; 0.1 Taylor Big skate; loess method 2015 ICATT blue shark 
                  # Taylor dogfish used CV as is so set Q_extraSD to 0.1 
default.Mean.weight.CV=0.2  #bit larger otherwise as it's the only signal for Southern2 selectivity

  #21.7 Drop single year size comp
Drop.single.year.size.comp=FALSE

  #21.8 Fit diagnostics
if(SS3.run=='final') do.SS3.diagnostics=TRUE  #very time consuming. Only run once model is defined.
if(SS3.run=='test') do.SS3.diagnostics=FALSE   
Retro_start=0; Retro_end=5 #Last 5 years of observations for retrospective analysis
Number.of.jitters=50              
Number.of.likelihood.profiles=10
delta.likelihood.profiles=0.2  #margine around Ro MLE for setting range of Ro values tested in like prof.
Approach.like.prof='SE' #set to 'min.plus' for sequence between [MLE -Number.of.likelihood.profiles] and [MLE +Number.of.likelihood.profiles]
Number.of.likelihood.profiles.h=6
diag.extras=""  #set to '"-nohess" to remove hessian estimation (much faster but no uncertainty)
#like.prof.case='faster'  #faster run, no hessian estimation
#like.prof.case='standard'  #as per r4ss (not estimating Hessian by setting extras)

#22. Bespoke Integrated size-based model 
if(Do.bespoke)
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

#23. Weight of Evidence
LoE.Weights=c(Spatial=0,COM=0,JABBA=0,integrated=1)  #if no integrated, use next highest
RiskColors=c('Negligible'="cornflowerblue",
             'Low'="chartreuse3",  #olivedrab3
             'Medium'="yellow1",
             'High'="orange",
             'Severe'="brown1")   #red2
Choose.probability="Depletion" #"B.over.Bmsy"  #use B/Bmys or B/K probabilities
Like.ranges=list(L1=c(0,0.0499999),
                 L2=c(0.05,0.2),
                 L3=c(0.20001,0.5),
                 L4=c(0.50001,1))
label_colors=c(Indicator='chocolate4',Non.indicator='cadetblue4',PSA.only='black')

#24. Average ratio L50:L95 for species with no L95 estimates
average.prop.L95_L50=mean(c(135.4/154.5,225/262,210/240,113/138,175/198,281/328,113/139,125/136,113/138))

#25. Define if using effort 
add.effort="NO"    
What.Effort="km.gn.hours"  #What effort to display?
#What.Effort="km.gn.days" 


#---2. Catch and effort data -----   

FisheryCodes=read_excel(handl_OneDrive('Data/Catch and Effort/FisheryCodeTable.xlsx'), sheet = "CAEStoFISHCUBE")

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
Southern.shark.fisheries=c('Joint Authority Southern Demersal Gillnet and Demersal Longline Managed Fishery'='JASDGDL',
                           'West Coast Demersal Gillnet and Demersal Longline (Interim) Managed Fishery'='WCDGDL',
                           'FBL condition 70 Power operated net hauler'='C070',
                           'Open Access in the West Coast Bioregion'='OAWC',
                           'Discards_TDGDLF'='Discards_TDGDLF',
                           'TEP_greynurse'='TEP_greynurse',
                           'TEP_dusky'='TEP_dusky')
Northern.shark.fisheries=c('Western Australian North Coast Shark Fishery'='WANCS',
                           'Open Access in the North Coast and Gascoyne Coast Bioregions'='OANCGC',
                           'Joint Authority Northern Shark Fishery'='JANS')

Non.WA.fisheries=c('GAB.trawl','Indonesia','NSW fisheries','NT','SA MSF','Taiwan')

do.pie.donut=FALSE 
if(do.pie.donut)
{
  
 #Catch pie donut
  fn.fig(handl_OneDrive("Analyses/Population dynamics/Proportional catch_WA"),2400,2400)
  KtCh.method%>%
    filter(!Data.set%in%c('Historic',Non.WA.fisheries))%>%
    mutate(Shark.Fishery=ifelse(FishCubeCode%in%Shark.Fisheries,'Shark fisheries',
                                'Other fisheries'))%>%
    group_by(Shark.Fishery,Gear)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))%>%
    PieDonut(aes(Shark.Fishery,Gear,count=LIVEWT.c),title='Proportional catch',
             explode = 2,r0=.2,r1=.8,maxx=1.45, titlesize=9,
             pieLabelSize=6.5,donutLabelSize=4.8,showPieName=FALSE)
  dev.off()
  
  
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
  
  #WA fisheries word cloud
  fn.fig(handl_OneDrive("Analyses/Population dynamics/Proportional catch_WA_word.cloud"),2400,2400)
  KtCh.method%>%
    filter(!Data.set%in%c('Historic',Non.WA.fisheries))%>%
    mutate(FishCubeCode=ifelse(FishCubeCode=="Discards_TDGDLF","JASDGDL",FishCubeCode),
           Shark.Fishery=ifelse(FishCubeCode%in%Southern.shark.fisheries,'Southern.shark.fisheries',
                                ifelse(FishCubeCode%in%Northern.shark.fisheries,'Northern.shark.fisheries',
                                       'Other')))%>%
    group_by(Shark.Fishery,FishCubeCode)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))%>%
    ggplot(aes(label = FishCubeCode, size = LIVEWT.c^(0.5),color=Shark.Fishery))+
    geom_text_wordcloud(show.legend = FALSE, family="Purisa") +
    scale_size_area(max_size = 20)+
    guides(size = "none")+
    theme_minimal()+
    ggtitle("WA's fisheries interacting with sharks and rays")+
    theme(plot.title = element_text(size=22,hjust=0.5,vjust = -30))+
    scale_color_manual(values=c(Other = "black",
                                Southern.shark.fisheries = 'forestgreen', Northern.shark.fisheries = 'darkorange3'))
  dev.off()  
  
  #Combine the two
  require(magick)
  Prop <- image_read(handl_OneDrive("Analyses/Population dynamics/Proportional catch_WA.tiff"))%>%
    image_resize("650x")%>%image_crop("632x600+10+10")
  Cloud <- image_read(handl_OneDrive("Analyses/Population dynamics/Proportional catch_WA_word.cloud.tiff"))%>%
    image_resize("800x")%>%image_crop("700x600+50+100")
  a=image_append(c(Prop, Cloud),stack=TRUE)
  image_write(a, path = handl_OneDrive("Analyses/Population dynamics/Proportional catch_WA_combined.tiff"),
              format = "tiff")

  
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


#---3. Life history, selectivity, TDGDLF mesh proportions, recruitment and depletion  data--------------------  
LH.data=read.csv(handl_OneDrive('Data/Life history parameters/Life_History.csv'))

SS_selectivity_init_pars=read.csv(handl_OneDrive('Analyses/Population dynamics/SS3.selectivity_pars.csv'))

mesh.prop.effort=read.csv(handl_OneDrive('Analyses/Catch and effort/mesh.proportional.effort.csv'))%>%mutate(Zone='Combined')
mesh.prop.effort.West=read.csv(handl_OneDrive('Analyses/Catch and effort/mesh.proportional.effort.West.csv'))%>%mutate(Zone='West')
mesh.prop.effort.Zone1=read.csv(handl_OneDrive('Analyses/Catch and effort/mesh.proportional.effort.Zone1.csv'))%>%mutate(Zone='Zone1')
mesh.prop.effort.Zone2=read.csv(handl_OneDrive('Analyses/Catch and effort/mesh.proportional.effort.Zone2.csv'))%>%mutate(Zone='Zone2')
mesh.prop.effort=rbind(mesh.prop.effort,mesh.prop.effort.West,mesh.prop.effort.Zone1,mesh.prop.effort.Zone2)
rm(mesh.prop.effort.West,mesh.prop.effort.Zone1,mesh.prop.effort.Zone2)

SS3.Rrecruitment.inputs=read.csv(handl_OneDrive('Analyses/Population dynamics/SS3.Rrecruitment.inputs.csv'))
SS3.tune_size_comp_effective_sample=read.csv(handl_OneDrive('Analyses/Population dynamics/SS3.tune_size_comp_effective_sample.csv'))
SS3.tune_size_comp_effective_sample_spatial=read.csv(handl_OneDrive('Analyses/Population dynamics/SS3.tune_size_comp_effective_sample_spatial.csv'))

Depletion.levels<- read_excel(handl_OneDrive('Analyses/Population dynamics/K_&_depletion.levels.xlsx'),sheet = "K_&_depletion.levels")

#---4. Productivity Susceptibility Analyses -----------------------------------------------  
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
if(First.run=="YES")
{
  write.csv(KtCh.method%>%
              mutate(Year=as.numeric(substr(FINYEAR,1,4)),
                     Name=capitalize(Name))%>%
              group_by(Name,Year,Data.set)%>%
              summarise(catch=sum(LIVEWT.c)),
            paste(Exprt,'Annual_ktch_by_species.and.data.set.csv',sep='/'),row.names = F)

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
  
}

#Export table of species identified in the catch by fishery  (all species including indicator species)
if(First.run=="YES")
{
  write.csv(Get.ktch$Table1%>%
              filter(!is.na(Name))%>%left_join(All.species.names%>%dplyr::select(SPECIES,Scien.nm),by='Scien.nm')%>%
              relocate(c(SPECIES))%>%
              arrange(SPECIES),
            paste(Exprt,'Table S1_All.species.caught.by.fishery.csv',sep='/'),row.names = F)
}

#Export species word cloud 
if(First.run=="YES")
{
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
}

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
PSA.out=PSA.fn(d=PSA.list,line.sep=.45,size.low=2.1,size.med=2.15,size.hig=2.5,W=10,H=10)
if(First.run=="YES")
{
  PSA.out$p
  ggsave(paste(Exprt,'Figure_PSA.tiff',sep='/'), width = W,height = H, dpi = 300, compression = "lzw")
}
PSA.out=PSA.out$PSA
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
#      'age_length' cannot be used as conditional age-at-length data because it comes from different gears, fleets, etc
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

  #6.1 add size composition C. brachyurus MSF (South Australia, pelagic longline)  
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

  #6.2 add size composition Angel sharks (GAB Trawl)
if(use.Gab.trawl)
{
  Species.data$`angel sharks`$Size_composition_Other=read.csv(handl_OneDrive('Data/Population dynamics/GABFIS_LengthSharks.csv'))%>%
    filter(Species.Csiro==37024002)%>%
    mutate(DATE=as.Date(Start.Date,"%d/%m/%Y"),
           Month=month(DATE),
           year=year(DATE),
           FINYEAR=ifelse(Month>6,paste(year,substr(year+1,3,4),sep='-'),
                          paste(year-1,substr(year,3,4),sep='-')),
           FL=with(LH.data%>%filter(SPECIES==24900),((Species.Size)-b_FL.to.TL )/a_FL.to.TL),
           SEX=Sex,
           SEX=ifelse(SEX=="Male ","Male",SEX),
           SEX=ifelse(SEX=='Female','F',ifelse(SEX=='Male','M',NA)))%>%
    filter(!is.na(SEX))%>%
    filter(!is.na(year))%>%
    dplyr::select(Month,FINYEAR,year,FL,SEX,Species.Quantity)%>% 
    type.convert(as.is = TRUE) %>% 
    uncount(Species.Quantity)
  
  Species.data$`angel sharks`$Size_composition_Other_Observations=data.frame(
    FINYEAR=sort(unique(Species.data$`angel sharks`$Size_composition_Other$FINYEAR)),
    Method='Trawl',
    zone="GAB.trawl",
    SPECIES=24900,
    N.shots=23,
    N.observations=NA)
}

  #6.3 add size composition gummy shark (GAB Trawl)
if(use.Gab.trawl & add.gummy.gab)
{
  Species.data$`gummy shark`$Size_composition_Other=read.csv(handl_OneDrive('Data/Population dynamics/GABFIS_LengthSharks.csv'))%>%
    filter(Species.Csiro==37017001)%>%
    mutate(DATE=as.Date(Start.Date,"%d/%m/%Y"),
           Month=month(DATE),
           year=year(DATE),
           FINYEAR=ifelse(Month>6,paste(year,substr(year+1,3,4),sep='-'),
                          paste(year-1,substr(year,3,4),sep='-')),
           FL=with(LH.data%>%filter(SPECIES==17001),((Species.Size)-b_FL.to.TL )/a_FL.to.TL),
           SEX=Sex,
           SEX=ifelse(SEX=="Male ","Male",SEX),
           SEX=ifelse(SEX=='Female','F',ifelse(SEX=='Male','M',NA)))%>%
    filter(!is.na(SEX))%>%
    filter(!is.na(year))%>%
    dplyr::select(Month,FINYEAR,year,FL,SEX,Species.Quantity)%>% 
    type.convert(as.is = TRUE) %>% 
    uncount(Species.Quantity)
  
  Species.data$`gummy shark`$Size_composition_Other_Observations=data.frame(
    FINYEAR=sort(unique(Species.data$`gummy shark`$Size_composition_Other$FINYEAR)),
    Method='Trawl',
    zone="GAB.trawl",
    SPECIES=17001,
    N.shots=17,
    N.observations=NA)
}

  #6.4 add size composition Spotted wobbegong and Western wobbegong to Wobbegongs
other.wobbies=c('Spotted wobbegong','Western wobbegong')
Wobbies.data=vector('list',length(other.wobbies))
names(Wobbies.data)=other.wobbies
    #6.4.1. bring in data
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
    #6.4.2. add to Wobbegongs
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

  #6.5 remove Pilbara trawl length comp as it's not used at all
for(i in 1:N.sp) 
{
  if(any(grepl('Size_composition_Pilbara_Trawl',names(Species.data[[i]]))))
  {
    Species.data[[i]]=Species.data[[i]][-match('Size_composition_Pilbara_Trawl',names(Species.data[[i]]))]
  }
}

  #6.6 remove NA sex in length composition data
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

  #6.7 remove nonsense length comp values
nems=c('Size_composition_dropline','Size_composition_NSF.LONGLINE',
       'Size_composition_Survey','Size_composition_Other',
       'Size_composition_West.6.5.inch.raw','Size_composition_West.7.inch.raw',
       'Size_composition_Zone1.6.5.inch.raw','Size_composition_Zone1.7.inch.raw',
       'Size_composition_Zone2.6.5.inch.raw','Size_composition_Zone2.7.inch.raw')
for(i in 1:N.sp) 
{
  LH=LH.data%>%filter(SPECIES==All.species.names%>%filter(SNAME==Keep.species[i])%>%pull(SPECIES))
  Max.TL=LH$Max.TL
  Max.FL=1.15*(Max.TL-LH$b_FL.to.TL)/LH$a_FL.to.TL
  Min.FL=LH$LF_o
  for(z in 1:length(nems))
  {
    if(nems[z]%in%names(Species.data[[i]]))
    {
      Species.data[[i]][[nems[z]]]=Species.data[[i]][[nems[z]]]%>%
                                        filter(FL<=Max.FL)%>%
                                        filter(FL>=Min.FL)
    }
  }
}

  #6.8 remove NA Ages in length-age data
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

  #6.9 remove Observer cpue data (double dipping with standardised catch rates)
for(s in 1:N.sp)
{
  iid=grep('CPUE_Observer_TDGDLF',names(Species.data[[s]]))
  if(length(iid)>0) Species.data[[s]]=Species.data[[s]][-iid]
}

  #6.10 Remove F series from TDGDLF due to structural uncertainty (different growth and M estimation) and sample size
for(s in 1:N.sp)
{
  iid=grep('Fishing.mortality.TDGDLF',names(Species.data[[s]]))
  if(length(iid)>0) Species.data[[s]]=Species.data[[s]][-iid]
}

  #6.11 Look at growth Cvs
Growth.CVs=vector('list',N.sp)
names(Growth.CVs)=Keep.species
for(i in 1:N.sp)
{
  if('age_length'%in% names(Species.data[[i]]))
  {
    Di=Species.data[[i]]$age_length
    Agess=sort(unique(Di$Age))
    siVis=Agess
    for(a in 1:length(Agess))
    {
      siVis[a]=sd(Di%>%filter(Age==Agess[a])%>%pull(FL))/mean(Di%>%filter(Age==Agess[a])%>%pull(FL))
    }
    Growth.CVs[[i]]=data.frame(Age=Agess,CV=siVis)
  }
}

  #6.12 Display annual proportional effort by zone and mesh
if(First.run=="YES")
{
  mesh.prop.effort%>%
    filter(!Zone=='Combined')%>%
    gather(Mesh,Prop,-c(finyear,Zone))%>%
    mutate(Mesh=ifelse(Mesh=='X165','6.5',ifelse(Mesh=='X178','7',NA)),
           year=as.numeric(substr(finyear,1,4)))%>%
    ggplot(aes(year,Prop,color=Mesh))+
    geom_point(size=3)+
    geom_line(linewidth=1.25)+
    facet_wrap(~Zone,ncol=1)+
    theme_PA()+
    theme(legend.position = 'top')+ylab('Proportion of effort')+xlab('Financial year')
  ggsave(handl_OneDrive('Analyses/Population dynamics/mesh prop effort.tiff'),
         width = 6,height = 6, dpi = 300, compression = "lzw")
  
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
  if(names(List.sp)[l]=="sandbar shark") LH$Fecu_a=LH$Fecu_b=NA  #non significant relationship  NEW
  if(!is.na(LH$Fecu_mean)) LH$Fecu_min=LH$Fecu_max=LH$Fecu_mean
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
                           b_FL.to.TL=LH$b_FL.to.TL)              
  
  #Calculate tmax and reset Max Age to mean
  List.sp[[l]]$tmax=with(List.sp[[l]],((1/Growth.F$k)*log((Growth.F$FL_inf-Lzero)/( (1-0.99)*Growth.F$FL_inf)))) #theoretical lifespan (Cortes & Taylor 2023)
  Max_Age=List.sp[[l]]$Max.age.F[1]
  Max_Age_max=List.sp[[l]]$Max.age.F[2]
  if(is.na(Max_Age_max))
  {
    Max_Age_max=round(Max_Age*Max.Age.up.Scaler)
    if(fill.NA.Max.Age.Max) List.sp[[l]]$Max.age.F[2]=Max_Age_max
  }
    
  if(reset.max.Age) List.sp[[l]]$Max.age.F=rep(ceiling(mean(c(Max_Age,Max_Age_max))),2)
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
    age=0:max(List.sp[[l]]$Max.age.F,na.rm=T)
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
 
#Compare Max Age and tmax 
if(First.run=="YES")
{
  d1=fn.get.stuff.from.list(List.sp,"Max.age.F")
  d1=do.call(rbind,d1)%>%
    data.frame%>%
    rename(all_of(c(Amax='X1',Amax_max='X2')))%>%
    rownames_to_column(var = "Species")
  
  d2=fn.get.stuff.from.list(List.sp,"tmax")
  d2=do.call(rbind,d2)%>%data.frame%>%
    rename(all_of(c(tmax='.')))%>%
    rownames_to_column(var = "Species")
  full_join(d1,d2,by='Species')%>%
    gather(Variable,Value,-Species)%>%
    ggplot(aes(Variable,Value))+
    geom_bar(stat="identity")+
    facet_wrap(~Species,scales='free')
  ggsave(handl_OneDrive("Analyses/Population dynamics/Max age vs tmax.tiff"),width = 10,height = 8, dpi = 300, compression = "lzw")
}

#Remove A.max_max if NA
for(i in 1:N.sp)if(is.na(List.sp[[i]]$Max.age.F[2])) List.sp[[i]]$Max.age.F=List.sp[[i]]$Max.age.F[1]

# Display survey and NSF length comp vs TLmax & TLinf
if(First.run=="YES")
{
  for(i in 1:N.sp)
  {
    Name=capitalize(names(Species.data)[i])
    if('Size_composition_Survey'%in%names(Species.data[[i]]))
    {
      fn.kmpr.Linf_length.comp(d=Species.data[[i]]$Size_composition_Survey,
                               FL.TL.conv=c(a=List.sp[[i]]$a_FL.to.TL,b=List.sp[[i]]$b_FL.to.TL),
                               Linf=List.sp[[i]]$Growth.F$FL_inf,
                               Linf.m=List.sp[[i]]$Growth.M$FL_inf,
                               TLmax=mean(List.sp[[i]]$TLmax))
      ggsave(paste0(Outputs,'/Data/1.',Name,"/",AssessYr,"/1_Inputs/Visualise data/Length comps_survey_Linf.tiff"),
             width=6,height=6,compression = "lzw") 
    }
    
    if('Size_composition_NSF.LONGLINE'%in%names(Species.data[[i]]))
    {
      fn.kmpr.Linf_length.comp(d=Species.data[[i]]$Size_composition_NSF.LONGLINE,
                               FL.TL.conv=c(a=List.sp[[i]]$a_FL.to.TL,b=List.sp[[i]]$b_FL.to.TL),
                               Linf=List.sp[[i]]$Growth.F$FL_inf,
                               Linf.m=List.sp[[i]]$Growth.M$FL_inf,
                               TLmax=mean(List.sp[[i]]$TLmax))
      ggsave(paste0(Outputs,'/Data/1.',Name,"/",AssessYr,"/1_Inputs/Visualise data/Length comps_NSF_Linf.tiff"),
             width=6,height=6,compression = "lzw") 
      
    }
  }
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
store.species.r_M.age.invariant=vector('list',N.sp)
names(store.species.r_M.age.invariant)=Keep.species
store.species.M_M.age.invariant=store.species.G_M.age.invariant=store.species.r_M.age.invariant
store.species.r_M.at.age=store.species.M_M.at.age=store.species.G_M.at.age=store.species.r_M.age.invariant
if(do.r.prior)  #0.015 sec per iteration per species   
{
  set.seed(1234)
  store.species.r=vector('list',N.sp)
  names(store.species.r)=Keep.species
  tic()
  for(l in 1:N.sp)  #0.013 sec per iteration per species; no parallel processing implementation as it stuffs up Mmin and Mmean
  {
    print(paste('r calculation -------------',Keep.species[l]))  
    #if no sd, replace with mean from those species with sd
    if(is.na(List.sp[[l]]$Growth.F$FL_inf.sd)) List.sp[[l]]$Growth.F$FL_inf.sd=0.038*List.sp[[l]]$Growth.F$FL_inf  
    if(is.na(List.sp[[l]]$Growth.F$k.sd)) List.sp[[l]]$Growth.F$k.sd=0.088*List.sp[[l]]$Growth.F$k     
    
    RESAMP="YES"
    
    linear.fec="NO"
    #if(names(List.sp)[l]%in%c("grey nurse shark","sandbar shark")) linear.fec="NO"
    
    #Get r prior
      #Age invariant M
    what.M<<-'age.invariant'
    M.averaging<<-"min" #Not Applicable as only Dureiul's methods used. If using multiple methods with equal weight, 'min' yields rmax, Cortes pers com
    GET.all.Ms=TRUE
    Amx=List.sp[[l]]$Max.age.F
    if(Keep.species[l]%in%species.too.high.M1) Amx=max(List.sp[[l]]$Max.age.F)
    r.prior.dist_M.age.invariant=with(List.sp[[l]],fun.rprior.dist(Nsims=NsimSS,
                                                                   K=Growth.F$k,
                                                                   LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,
                                                                   K.sd=Growth.F$k.sd,
                                                                   LINF.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL,
                                                                   k.Linf.cor,
                                                                   Amax=Amx,
                                                                   MAT=unlist(Age.50.mat),
                                                                   FecunditY=Fecundity,
                                                                   Cycle=Breed.cycle,
                                                                   BWT=BwT,AWT=AwT,
                                                                   LO=Lzero*a_FL.to.TL+b_FL.to.TL)) #size vars as TL
      #Mortality at age
    what.M<<-'at.age'
    M.averaging<<-"mean"
    GET.all.Ms=FALSE
    mean.donotwork=NULL
    #mean.donotwork=c("grey nurse shark","spurdogs") #when M.average='mean', endless loop due to 
                                                # low fecundity (greynurse & spurdogs)
                                                # or high mortality (due to k for sandbar)
    if(!Keep.species[l]%in%mean.donotwork)  
    {
      Amx=List.sp[[l]]$Max.age.F
      if(Keep.species[l]%in%species.too.high.M1) Amx=max(List.sp[[l]]$Max.age.F)
      r.prior.dist_M.at.age=with(List.sp[[l]],fun.rprior.dist(Nsims=NsimSS,
                                                              K=Growth.F$k,
                                                              LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,
                                                              K.sd=Growth.F$k.sd,
                                                              LINF.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL,
                                                              k.Linf.cor,
                                                              Amax=Amx,
                                                              MAT=unlist(Age.50.mat),
                                                              FecunditY=Fecundity,
                                                              Cycle=Breed.cycle,
                                                              BWT=BwT,AWT=AwT,
                                                              LO=Lzero*a_FL.to.TL+b_FL.to.TL)) #size vars as TL
    }
    if(Keep.species[l]%in%mean.donotwork) r.prior.dist_M.at.age=r.prior.dist_M.age.invariant
    
    store.species.r[[l]]=list(r.prior.dist_M.age.invariant=r.prior.dist_M.age.invariant,
                              r.prior.dist_M.at.age=r.prior.dist_M.at.age)
  }
  toc()    
  clear.log('fun.rprior.dist')
  clear.log('fun.Leslie')
  if("plyr"%in%.packages()) detach("package:plyr", unload=TRUE)
  
  #Extract quantities
  store.species.r_M.age.invariant=lapply(store.species.r, function(x) x[[1]])
  store.species.r_M.at.age=lapply(store.species.r, function(x) x[[2]])
  names(store.species.r_M.age.invariant)=names(store.species.r_M.at.age)=Keep.species
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
    Nms=names(store.species.r_M.age.invariant[[l]]$Input.pars[[1]])
    LH.d=matrix(unlist(list.flatten(store.species.r_M.age.invariant[[l]]$Input.pars)),nrow=List.sp[[l]]$NsimSS,ncol=length(Nms),byrow = T)
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
    write.csv(with(store.species.r_M.age.invariant[[l]],data.frame(shape=shape,rate=rate,mean=mean,sd=sd)),'r.prior_M.age.invariant.csv',row.names = F)
    write.csv(with(store.species.r_M.at.age[[l]],data.frame(shape=shape,rate=rate,mean=mean,sd=sd)),'r.prior_M.at.age.csv',row.names = F)
    
    #export G
    out.G_M.age.invariant=with(store.species.r_M.age.invariant[[l]],data.frame(mean=mean(G),sd=sd(G)))
    write.csv(out.G_M.age.invariant,'G.prior_M.age.invariant.csv',row.names = F)
    store.species.G_M.age.invariant[[l]]=out.G_M.age.invariant
    
    out.G_M.at.age=with(store.species.r_M.at.age[[l]],data.frame(mean=mean(G),sd=sd(G)))
    write.csv(out.G_M.at.age,'G.prior_M.at.age.csv',row.names = F)
    store.species.G_M.at.age[[l]]=out.G_M.at.age
    
    #export M
    out.M=store.species.r_M.age.invariant[[l]]$M
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
    store.species.M_M.age.invariant[[l]]=out.M
    write.csv(out.M,"M_M.age.invariant.csv",row.names=FALSE)
    
    out.M=store.species.r_M.at.age[[l]]$M
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
    store.species.M_M.at.age[[l]]=out.M
    write.csv(out.M,"M_M.at.age.csv",row.names=FALSE)
  }

  #Plot each M estimator 
  for(l in 1:N.sp)
  {
    print(paste("plot M derived for each estimator","--",List.sp[[l]]$Name))
    store.species.r_M.age.invariant[[l]]$nat.mort.sim%>%
      ggplot(aes(x=Age, y=M.mean,color=Method))+
      geom_point()+
      geom_errorbar(aes(ymin=M.mean-M.sd, ymax=M.mean+M.sd))+
      facet_wrap(~Method,ncol=2)+
      theme_PA()+
      theme(legend.position="none",
            axis.text.x = element_text(size=7,angle = 90, hjust=1))+
      ylab('M (+/- SD)')+scale_y_continuous(limits = c(0, NA))
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",AssessYr,"/demography/M_all.tiff",sep=''),
           width = 8,height = 6, dpi = 300, compression = "lzw")
  }
  
  #Compare M constant and M at age Natural mortality
  for(l in 1:length(Lista.sp.outputs))
  {
    print(paste("Compare M at age VS M invariant","--",names(Lista.sp.outputs)[l]))
    d1=store.species.M_M.age.invariant[Lista.sp.outputs[[l]]] 
    d2=store.species.M_M.at.age[Lista.sp.outputs[[l]]]
    for(pp in 1:length(d1))
    {
      d1[[pp]]=d1[[pp]]%>%data.frame%>%
        mutate(iteration=row_number())%>%
        gather(Age,M,-iteration)%>%
        mutate(Age=as.numeric(gsub('X','',Age)),
               Scenario='Constant M',
               Species=capitalize(names(d1)[pp]))
      
      d2[[pp]]=d2[[pp]]%>%data.frame%>%
        mutate(iteration=row_number())%>%
        gather(Age,M,-iteration)%>%
        mutate(Age=as.numeric(gsub('X','',Age)),
               Scenario='M at age',
               Species=capitalize(names(d1)[pp]))
    }
    p=rbind(do.call(rbind,d1),do.call(rbind,d2))
    STRX=13
    WID=8
    if(names(Lista.sp.outputs)[l]=="Other.sp")
    {
      STRX=10
      WID=10
    }
      
    p%>%
      group_by(Age,Scenario,Species)%>%
      summarise(Mean=mean(M),
                SD=sd(M))%>%
      ggplot(aes(Age,Mean,color=Scenario))+
      geom_point(size=2)+geom_line(linewidth=.8,linetype=2)+ 
      geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=1.5, 
                    position=position_dodge(0.05))+
      facet_wrap(~Species,scales='free')+scale_y_continuous(limits = c(0, NA))+
      theme_PA(strx.siz=STRX)+theme(legend.position = "top",
                       legend.title=element_blank())+ylab('Mean +/- SD')
    ggsave(paste(Rar.path,'/Prior_M_',names(Lista.sp.outputs)[l],'_M.at.age_vs_Constant.M.tiff',sep=''),
           width = WID,height = 8,compression = "lzw")
    
  }
 
  rm(store.species.r)
}
if(!do.r.prior)
{
  for(l in 1:N.sp) 
  {
    hndl.dummy=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                     AssessYr,"/demography",sep='')
    store.species.r_M.age.invariant[[l]]=read.csv(paste(hndl.dummy,"/r.prior_M.age.invariant.csv",sep=''))
    store.species.G_M.age.invariant[[l]]=read.csv(paste(hndl.dummy,"/G.prior_M.age.invariant.csv",sep=''))
    store.species.M_M.age.invariant[[l]]=read.csv(paste(hndl.dummy,"/M_M.age.invariant.csv",sep=''))
    
    store.species.r_M.at.age[[l]]=read.csv(paste(hndl.dummy,"/r.prior_M.at.age.csv",sep=''))
    store.species.G_M.at.age[[l]]=read.csv(paste(hndl.dummy,"/G.prior_M.at.age.csv",sep=''))
    store.species.M_M.at.age[[l]]=read.csv(paste(hndl.dummy,"/M_M.at.age.csv",sep=''))
    
    rm(hndl.dummy)
  }
}

  #display priors
if(do.r.prior)
{
  #Base case
  for(l in 1:length(Lista.sp.outputs))
  {
    print(paste("Display r prior for","--",names(Lista.sp.outputs)[l]))
    STXSIZ=16
    Xmax=0.6
    if(names(Lista.sp.outputs)[l]=="Other.sp")
    {
      STXSIZ=12.5
      Xmax=0.96
    }
    fn.display.priors(d=store.species.r_M.age.invariant,
                      sp=Lista.sp.outputs[[l]],
                      XLAB=expression(paste(plain("Maximum intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),
                      XLIM=c(0,Xmax),
                      Strx.siz=STXSIZ)
    ggsave(paste(Rar.path,'/Prior_r_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = 12,height = 10,compression = "lzw")
  }
  
  #M.at.age vs Constant M
  for(l in 1:length(Lista.sp.outputs))
  {
    print(paste("Display M.at.age vs Constant M for","--",names(Lista.sp.outputs)[l]))
    STXSIZ=16
    Xmax=0.6
    if(names(Lista.sp.outputs)[l]=="Other.sp")
    {
      STXSIZ=12
      Xmax=0.96
    }
    out=fn.display.prior.sensitivity(d=store.species.r_M.age.invariant,
                                     d2=store.species.r_M.at.age,
                                     sp=Lista.sp.outputs[[l]],
                                     XLAB=expression(paste(plain("Maximum intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),
                                     Strx.siz=STXSIZ,
                                     Scen1='Constant M',
                                     Scen2='M.at.age')
    ggsave(paste(Rar.path,'/Prior_r_',names(Lista.sp.outputs)[l],'_M.at.age_vs_Constant.M.tiff',sep=''),
           width = 12,height = 10,compression = "lzw")
    write.csv(out,paste(Rar.path,'/Prior_r_',names(Lista.sp.outputs)[l],'_M.at.age_vs_Constant.M.csv',sep=''),row.names=F)
  }
  
  pdf(handl_OneDrive('Analyses/Population dynamics/M_constant_VS_at.age_store_M.pdf'))
  smart.par(N.sp,c(1,1.5,1,1),c(1,1,1,1),c(3, 1, 0))
  for(x in 1:N.sp)
  {
    print(paste("Display each run of M.at.age vs Constant M for","--",Keep.species[x]))
    id=match(Keep.species[x],names(store.species.M_M.at.age))
    plot(unlist(store.species.M_M.at.age[[id]][1,]),ylim=c(0,0.8),xaxt='n',main=Keep.species[x],xlab='Length',ylab="M")
    for(q in 1:nrow(store.species.M_M.at.age[[id]]))lines(unlist(store.species.M_M.at.age[[id]][q,]))
    for(q in 1:nrow(store.species.M_M.at.age[[id]]))lines(unlist(store.species.M_M.age.invariant[[id]][q,]),col=2)
  }
  dev.off()
}

  #compare M and r
if(do.r.prior)
{
  Omit.these=NULL
  CompR=data.frame(Name=names(store.species.r_M.age.invariant),r=NA,M=NA)
  for(l in 1:N.sp)
  {
    CompR$r[l]=store.species.r_M.age.invariant[[l]]$mean
    CompR$M[l]=mean(unlist(store.species.M_M.age.invariant[[l]]),na.rm=T)
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
    stat_poly_eq(aes(label =  paste(after_stat(eq.label),after_stat(adj.rr.label),after_stat(p.value.label),
                                    sep = "*\", \"*")),
                 formula = my_formula, parse = TRUE,
                 label.y = "top", label.x = "right", size = 4)+
    geom_text_repel(segment.colour='black',col='black',box.padding = 0.5) + 
    scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))+ 
    theme_PA(axs.T.siz=14,axs.t.siz=12)+
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
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
for(r in 1:length(RESILIENCE)) RESILIENCE[[r]]=Res.fn(store.species.r_M.age.invariant[[r]]$mean,Def="Haddon")
clear.log('Res.fn')

#---11. Extract experimental selectivity at age and at size-----------------------------------------------------------------------
#note: For gillnet (dome-shape) sel., fix selectivity parameters in SS3 to untangle selectivity from fishing mortality.
#      For species with no gillnet sel. profile, set to closest species or family
Sel.equivalence=data.frame(
  Name=c("great hammerhead","scalloped hammerhead",
         "grey nurse shark",
         "wobbegongs",
         "sawsharks",
         "angel sharks",
         "lemon shark","milk shark","pigeye shark","shortfin mako","spinner shark","tiger shark",
         "spurdogs"),
  Equivalence=c(rep("Smooth hammerhead",2),
                "Hexanchidae",
                "Orectolobidae",
                "Pristiophoridae",
                "Squatinidae",
                rep('Carcharhinidae',6),
                'Squalidae'))
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
    if(!"X16.5"%in%names(GN.sel.at.totalength))
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
             Sel.combined=Sel.combined/max(Sel.combined))
    Selectivity.at.age[[l]]=GN.sel.at.age[,c('TL','Age','Sel.combined','X16.5','X17.8','type')]
    
    GN.sel.at.totalength=GN.sel.at.totalength%>%
      mutate(Sum.sel=X16.5+X17.8,
             Sel.combined=Sum.sel/max(Sum.sel),
             Sel.combined=Sel.combined/max(Sel.combined)) 
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
    if('Srvy.FixSt'%in%names(Species.data[[l]]))
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
    if('annual.abundance.NSF_relative'%in%names(Species.data[[l]]))
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
    if('annual.abundance.basecase.monthlyWest'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.monthlyWest
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.monthly.West_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.monthly.West=d
      rm(d)
    }
    if('annual.abundance.basecase.monthlyZone1'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.monthlyZone1
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.monthly.Zone1_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.monthly.Zone1=d
      rm(d)
    }
    if('annual.abundance.basecase.monthlyZone2'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.monthlyZone2
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.monthly.Zone2_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.monthly.Zone2=d
      rm(d)
    }
    if('annual.abundance.basecase.daily_relative'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.daily
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.daily_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.daily=d
      rm(d)
    }
    if('annual.abundance.basecase.dailyWest'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.dailyWest
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.daily.West_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.daily.West=d
      rm(d)
    }
    if('annual.abundance.basecase.dailyZone1'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.dailyZone1
      if(Abundance.error.dist=="Normal") d=Species.data[[l]]$annual.abundance.basecase.daily.Zone1_relative
      d=d%>%
        mutate(yr.f=as.numeric(substr(Finyear,1,4)))%>%
        filter(yr.f<=Last.yr.ktch.numeric)
      if(nrow(d)>=Min.cpue.yrs & !any(apply( Filter(is.numeric, d),2,is.infinite))) dummy$TDGDLF.daily.Zone1=d
      rm(d)
    }
    if('annual.abundance.basecase.dailyZone2'%in%names(Species.data[[l]]))
    {
      if(Abundance.error.dist=="Lognormal") d=Species.data[[l]]$annual.abundance.basecase.dailyZone2
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

#display Survey and TDGDLF cpues     
if(First.run=="YES")
{
  for(l in 1:N.sp)
  {
    A=Catch.rate.series[[l]]
    IID=grep(paste(c("Survey","TDGDLF.monthly","TDGDLF.daily"),collapse='|'),names(A))
    if(length(IID)>0)
    {
      A=A[IID]
      if(!is.null(A))
      {
        print(paste("Display all cpues --------",names(Species.data)[l]))
        
        fn.fig(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(List.sp[[l]]$Name),"/",AssessYr,
                     "/1_Inputs/Visualise data/All cpues",sep=''),2000,2000) 
        smart.par(n.plots=length(A),MAR=c(2,3,1,1),OMA=c(2.5,1,.05,2.5),MGP=c(1.8,.5,0))
        par(cex.lab=1.5,las=1)
        Mx.yr=max(sapply(A, function(x) max(x$yr.f, na.rm=TRUE)))
        Min.yr=min(sapply(A, function(x) min(x$yr.f, na.rm=TRUE)))
        
        for(i in 1:length(A))
        {
          with(A[[i]],{
            plot(yr.f,Mean,col='orange',xlim=c(Min.yr,Mx.yr),ylim=c(0,max(UP.CI,na.rm=T)),pch=19,cex=1.15,
                 main=names(A)[i],ylab='',xlab="")
            segments(yr.f,LOW.CI,yr.f,UP.CI,col='orange')
          })
          if(names(A[i])=='TDGDLF.daily')
          {
            with(A[[i]]%>%filter(yr.f%in%c(2007,2008)),{
              points(yr.f,Mean,col='brown4',pch=19,cex=1.15)
              segments(yr.f,LOW.CI,yr.f,UP.CI,col='brown4')
            })
          }
        }
        mtext("Financial year", side = 1, line = 1,outer=T)
        mtext("CPUE", side = 2, line = -1,las=3,outer=T)
        dev.off()
      }
    }

    rm(A,IID)
  }
}

#Remove non-representative series  
for(l in 1:N.sp) 
{
  Neim=Keep.species[l]
  dummy=Catch.rate.series[[l]]
  drop.this=NULL
  if(Neim%in%other_not.representative)  drop.this=c(drop.this,match('Other',names(dummy)))  
  if(Neim%in%survey_not.representative) drop.this=c(drop.this,match('Survey',names(dummy)))  
  if(Neim%in%NSF_not.representative)    drop.this=c(drop.this,match('NSF',names(dummy)))  
  if(Neim%in%tdgdlf_not.representative)
  {
    drop.this=c(drop.this,match('TDGDLF.monthly',names(dummy)))
    drop.this=c(drop.this,match('TDGDLF.daily',names(dummy)))
  }
  if(!is.null(drop.this))
  {
    drop.this=subset(drop.this,!is.na(drop.this))
    dummy=dummy[-drop.this]
  }
  if(length(dummy)>0) Catch.rate.series[[l]]=dummy
  if(length(dummy)==0) Catch.rate.series[l]=list(NULL) 
  rm(dummy,drop.this)
}



#---13. Calculate Steepness ----------------------------------------------------------------------- 
  #calculate prior
store.species.steepness_M.at.age=vector('list',N.sp)
names(store.species.steepness_M.at.age)=Keep.species
store.species.alpha_M.age.invariant=store.species.alpha_M.at.age=store.species.steepness_M.age.invariant=store.species.steepness_M.at.age
if(do.steepness)   #0.011 sec per iteration per species
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
    GET.all.Ms=FALSE
    
      #1. M at age
    M.averaging<<-"mean"   #Not Applicable as only Dureiul's methods used (previously, 'min' of multiple methods yielded too high h values for all species). 
    what.M<<-'at.age'
    steepNs_M.at.age=with(List.sp[[l]],fun.steepness(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,
                                                   Linf.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL,k.sd=Growth.F$k.sd,
                                                   first.age=First.Age,sel.age=SEL,F.mult=0,
                                                   Amax=Max.age.F,MAT=unlist(Age.50.mat),
                                                   FecunditY=Fecundity,Cycle=Breed.cycle,
                                                   sexratio=0.5,spawn.time = 0,
                                                   AWT=AwT,BWT=BwT,LO=Lzero*a_FL.to.TL+b_FL.to.TL,
                                                   Resamp=RESAMP,simsout=SSS.sims))
      #2. Age invariant M
    #M.averaging<<-"min"
    what.M<<-'age.invariant'
    steepNs_M.age.invariant=with(List.sp[[l]],fun.steepness(Nsims=NsimSS,K=Growth.F$k,LINF=Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,
                                                  Linf.sd=Growth.F$FL_inf.sd*a_FL.to.TL+b_FL.to.TL,k.sd=Growth.F$k.sd,
                                                  first.age=First.Age,sel.age=SEL,F.mult=0,
                                                  Amax=Max.age.F,MAT=unlist(Age.50.mat),
                                                  FecunditY=Fecundity,Cycle=Breed.cycle,
                                                  sexratio=0.5,spawn.time = 0,
                                                  AWT=AwT,BWT=BwT,LO=Lzero*a_FL.to.TL+b_FL.to.TL,
                                                  Resamp=RESAMP,simsout=SSS.sims))
    rm(RESAMP,M.averaging,linear.fec)
    store.species.h[[l]]=list(steepNs_M.at.age=steepNs_M.at.age,steepNs_M.age.invariant=steepNs_M.age.invariant)
  }
  toc() 

  clear.log('fun.steepness')
  clear.log('Alpha.Brooks')
  
  #Extract quantities
  steepNs_M.at.age=lapply(store.species.h, function(x) x[[1]])
  steepNs_M.age.invariant=lapply(store.species.h, function(x) x[[2]])
  names(steepNs_M.at.age)=names(steepNs_M.age.invariant)=Keep.species
  for(l in 1:N.sp)
  {
    print(paste("extract steepness value ","--",List.sp[[l]]$Name))
    PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
               capitalize(List.sp[[l]]$Name),"/",AssessYr,"/steepness",sep='')
    if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))
    setwd(PATH)
    
    #export h, alpha and M
    steepNs=steepNs_M.at.age[[l]]
    out.h=with(steepNs,data.frame(mean=mean,sd=sd))
    write.csv(out.h,'h.prior_M.at.age.csv',row.names = F) 
    store.species.steepness_M.at.age[[l]]=out.h
    write.csv(steepNs$Alpha,'Alpha_M.at.age.csv',row.names = F) 
    store.species.alpha_M.at.age[[l]]=steepNs$Alpha
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
    write.csv(out.M,"M_M.at.age.csv",row.names=FALSE)
    write.csv(steepNs$Runs,"Life.history_M.at.age.csv",row.names=FALSE)
    
    steepNs=steepNs_M.age.invariant[[l]]
    out.h=with(steepNs,data.frame(mean=mean,sd=sd))
    write.csv(out.h,'h.prior_M.age.invariant.csv',row.names = F) 
    store.species.steepness_M.age.invariant[[l]]=out.h
    write.csv(steepNs$Alpha,'Alpha_M.age.invariant.csv',row.names = F) 
    store.species.alpha_M.age.invariant[[l]]=steepNs$Alpha
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
    write.csv(out.M,"M_M.age.invariant.csv",row.names=FALSE)
    write.csv(steepNs$Runs,"Life.history_M.age.invariant.csv",row.names=FALSE)
    
    rm(steepNs,out.M)
    
  }
  
  rm(store.species.h)
}
if(!do.steepness)
{
  for(l in 1: N.sp)
  {
    store.species.steepness_M.at.age[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                                AssessYr,"/steepness/h.prior_M.at.age.csv",sep=''))
    store.species.alpha_M.at.age[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                            AssessYr,"/steepness/Alpha_M.at.age.csv",sep=''))
    
    store.species.steepness_M.age.invariant[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                                       AssessYr,"/steepness/h.prior_M.age.invariant.csv",sep=''))
    store.species.alpha_M.age.invariant[[l]]=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",
                                                   AssessYr,"/steepness/Alpha_M.age.invariant.csv",sep=''))
  }
}

  #compare steepness and r   
h_too.high=c("great hammerhead","tiger shark","smooth hammerhead","whiskery shark","zebra shark")
h_too.low=c('sawsharks','green sawfish','milk shark') 
h_too.long.converge=NULL
Omit.these.h=capitalize(c(h_too.high,h_too.low))
if(do.steepness)
{
  CompR=data.frame(Name=names(store.species.steepness_M.at.age),
                   h=unlist(sapply(store.species.steepness_M.at.age, `[`, 1)),
                   h.sd=unlist(sapply(store.species.steepness_M.at.age, `[`, 2)),
                   r=unlist(sapply(store.species.r_M.age.invariant, `[`, 3)),
                   r.sd=unlist(sapply(store.species.r_M.age.invariant, `[`, 4)))
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
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    ylim(0,1)+
    geom_errorbar(aes(ymin=h-h.sd, ymax=h+h.sd),colour=COL)+
    geom_errorbarh(aes(xmin=r-r.sd, xmax=r+r.sd),colour=COL)
  p
  ggsave(handl_OneDrive('Analyses/Population dynamics/Steepness_vs_r.tiff'), 
         width = 8,height = 8, dpi = 300, compression = "lzw")
  
}


#Recalculate steepness for species with too high/low h using linear model of r on h
#note: some h values deemed too high (following E Cortes discussion); some deemed too low (life history mispecificaton)
Mod.Pred=read.csv(handl_OneDrive('Analyses/Population dynamics/Steepness_vs_r_coeff.csv'))

store.species.steepness.S2=fn.get.stuff.from.list(store.species.steepness_M.at.age,"mean")   
dis.sp.h=tolower(Omit.these.h)
if("dwarf sawfish" %in% Keep.species) dis.sp.h=c(dis.sp.h,"dwarf sawfish","freshwater sawfish")
for(s in 1:length(dis.sp.h))
{
  set.seed(1234)
  id=match(dis.sp.h[s],names(store.species.steepness_M.at.age))
  new.h=store.species.r_M.age.invariant[[id]]$mean*Mod.Pred$slope+Mod.Pred$intercept
  store.species.steepness.S2[[id]]=rnorm(1,new.h,new.h/100)
}
if(test.Sedar)
{
  store.species.steepness.S2$`dusky shark`=Dusky.Sedar
  store.species.steepness.S2$`sandbar shark`=Sandbar.Sedar
}

  #display h2 priors  
if(do.steepness)
{
  #Base case (M.at.age)
  for(l in 1:length(Lista.sp.outputs))
  {
    fn.display.steepness(d=store.species.steepness_M.at.age,
                         sp=Lista.sp.outputs[[l]],
                         XLAB="Steepness (h)",
                         XLIM=c(0.2,1))
    ggsave(paste(Rar.path,'/Prior_steepness_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = 12,height = 10,compression = "lzw")
  }
  
  #M.at.age vs Constant M vs S2 
  for(l in 1:length(Lista.sp.outputs))
  {
    STXSIZ=16
    Xmax=0.6
    if(names(Lista.sp.outputs)[l]=="Other.sp")
    {
      STXSIZ=13
      Xmax=0.96
    }
    out=fn.display.steepness.sensitivity(d=store.species.steepness.S2,
                                         d1=store.species.steepness_M.at.age,
                                         d2=store.species.steepness_M.age.invariant,
                                         sp=Lista.sp.outputs[[l]],
                                         XLAB="Steepness (h)",
                                         Strx.siz=STXSIZ,
                                         Scen1='steepness S2',
                                         Scen2='M at age',
                                         Scen3='M invariant')
    ggsave(paste(Rar.path,'/Prior_steepness_',names(Lista.sp.outputs)[l],'_M.at.age_vs_Constant.M.tiff',sep=''),
           width = 12,height = 10,compression = "lzw")
    write.csv(out,paste(Rar.path,'/Prior_steepness_',names(Lista.sp.outputs)[l],'_M.at.age_vs_Constant.M.csv',sep=''),row.names=F)
  }
}
clear.log('fn.display.priors')
clear.log('M.fun')
clear.log('fn.display.steepness')

do.this=FALSE
if(do.this)
{
  A=vector('list',N.sp)
  for(l in 1:N.sp)A[[l]]=data.frame(Species=Keep.species[l],r=store.species.r_M.age.invariant[[l]]$mean,h=store.species.steepness.S2[[l]])
  do.call(rbind,A)%>%
    ggplot(aes(r,h, label =Species))+
    geom_point()+geom_text_repel(segment.colour='black',col='black',box.padding = 0.5)
}

#Compare r and h with published estimates
if(First.run=='YES')
{
  Tab.r.and.h=data.frame(Species=Keep.species,r_Matias=NA,h_Matias=NA,h2_Matias=NA)
  for(l in 1:N.sp)
  {
    Tab.r.and.h$r_Matias[l]=round(store.species.r_M.age.invariant[[l]]$mean,2)
    Tab.r.and.h$h_Matias[l]=round(store.species.steepness_M.at.age[[l]]$mean,2)
    Tab.r.and.h$h2_Matias[l]=round(store.species.steepness.S2[[l]],2)
  }
  Tab.r.and.h=Tab.r.and.h%>%left_join(Demo.published.values,by='Species')
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
for(l in 1:N.sp) Fmsy.M.scaler[[l]]=Cortes.Brooks.2018(alpha=median(unlist(store.species.alpha_M.at.age[[l]])))
for(l in 1:length(dis.sp.h)) #reset dis.sp.h consistently with h resetting
{
  s=match(dis.sp.h[l],names(Fmsy.M.scaler))
  Fmsy.M.scaler[[s]]=0.5
}

clear.log('store.species.alpha_M.at.age')
clear.log('store.species.alpha_M.age.invariant')

#1. Species-specific proxy to Bmsy.K based on Cortes et al 2012
R=function(r,G) 0.633-0.187*log(r*G) 
BmsyK.species=data.frame(Species=names(store.species.r_M.age.invariant),r=NA,G=NA)
for(l in 1:N.sp)
{
  BmsyK.species$r[l]=store.species.r_M.age.invariant[[l]]$mean
  BmsyK.species$G[l]=store.species.G_M.age.invariant[[l]]$mean
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

#---15. Create list with modelling arguments  ----- 
#note: For integrated model, all pars calculated in log space.
#       ln_RZERO is in 1,000 individuals so do 10 times the largest catch divided by 
#       average weight and divided by 1,000. Best units to work in are 1,000 individuals for 
#       numbers, catch in tonnes and length-weight in kg as all cancels out and predicted biomass
#       end up being in tonnes
fn.source1("Pin_file_and_model_arguments.r")


#---16. Export .dat, extract SS3 selectivities and some prelim analysis----- 
#remotes::install_github("r4ss/r4ss")
library(r4ss)  
if(First.run=="YES")
{
  HandL.out=handl_OneDrive("Analyses/Population dynamics/1.")
  fn.source1("Organise data.R")
}
  

#Compare empirical selectivity with observed size composition 
if(First.run=="YES")
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
      ggtitle(paste(Title,' (',unique(Sel$type),' selectivity)',sep=''))
    
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
        #All meshes combined
        fn.compare.sel.size.comp(Title="All meshes combined",
                                 Sel=data.frame(type=Selectivity.at.totalength[[l]]%>%pull(type),
                                                TL=Selectivity.at.totalength[[l]]%>%pull(TL),
                                                Sel=Selectivity.at.totalength[[l]]%>%
                                                  pull(Sel.combined)),
                                 size=do.call(rbind,TDGDLF.size.comp),
                                 FL_TL=data.frame(inter=b_FL.to.TL,
                                                  slope=a_FL.to.TL),
                                 MX=10*round(TLmax/10))
        
        dev.off()
      }
      
      #NSF (selectivity set at maturity)
      NSF.size.comp=Species.data[[l]][grep(paste(c('NSF.LONGLINE'),collapse="|"),
                                           names(Species.data[[l]]))]
      if(length(NSF.size.comp)>0)
      {
        pdf(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),
                  "/",AssessYr,"/1_Inputs/Visualise data","/Compare sel and obs size comp_NSF.pdf",sep=''))
          fn.compare.sel.size.comp(Title='NSF - longline (sel pars set at 50% & 95% maturity)',
                                   Sel=data.frame(type='Species',
                                                  TL=midpt,
                                                  Sel=round(1/(1+(exp(-log(19)*((midpt-TL.50.mat)/(TL.95.mat-TL.50.mat))))),3)),
                                   size=NSF.size.comp$Size_composition_NSF.LONGLINE,
                                   FL_TL=data.frame(inter=b_FL.to.TL,
                                                    slope=a_FL.to.TL),
                                   MX=10*round(TLmax/10))
        dev.off()
      }
      
      #Other
      Other.size.comp=Species.data[[l]][grep('Size_composition_Other',names(Species.data[[l]]))]
      Other.size.comp=Other.size.comp[-grep("Observations",names(Other.size.comp))]
      if(length(Other.size.comp)>0)
      {
        pdf(paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),
                  "/",AssessYr,"/1_Inputs/Visualise data","/Compare sel and obs size comp_Other.pdf",sep=''))
        #All meshes combined
        fn.compare.sel.size.comp(Title="All meshes combined",
                                 Sel=data.frame(type=Selectivity.at.totalength[[l]]%>%pull(type),
                                                TL=Selectivity.at.totalength[[l]]%>%pull(TL),
                                                Sel=Selectivity.at.totalength[[l]]%>%
                                                  pull(Sel.combined)),
                                 size=do.call(rbind,Other.size.comp),
                                 FL_TL=data.frame(inter=b_FL.to.TL,
                                                  slope=a_FL.to.TL),
                                 MX=10*round(TLmax/10))
        dev.off()
      }

      detach(List.sp[[l]])
      print(paste('-------------Compare empirical selectivity & observed size comp ----',List.sp[[l]]$Name))
      
    }
  }
}

#Extract SS selectivity parameters 
#note: do this only once and update 'SS3.selectivity_pars.csv'.
if(Extract.SS.parameters) fn.source1('Re fit SS3 selectivity.R')


# Compare observed size comp and assumed SS selectivity 
if(First.run=="YES")
{
  fn.source1("SS_selectivity functions.R")
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
        for(x in 1:length(d.list))
        {
          d.list[[x]]$Fleet=str_remove(str_remove(names(d.list)[x],'Size_composition_'),'.inch.raw')
          if(!'Size.type'%in%names(d.list[[x]])) d.list[[x]]$Size.type=NA
        }
        d.list=do.call(rbind,d.list)
        d.list=d.list%>%mutate(fleet=ifelse(grepl(paste(c('West','Zone'),collapse='|'),Fleet),'TDGDLF',Fleet),
                               Fleet=ifelse(fleet=='TDGDLF' & year<=2005,'Southern.shark_1',
                                            ifelse(fleet=='TDGDLF' & year>2005,'Southern.shark_2',
                                                   ifelse(fleet=='NSF.LONGLINE','Northern.shark',
                                                          fleet))),
                               TL=FL*List.sp[[i]]$a_FL.to.TL+List.sp[[i]]$b_FL.to.TL)
        
        
        if(MN.SZE=="size.at.birth") Min.size.bin=10*round((List.sp[[i]]$Lzero*List.sp[[i]]$a_FL.to.TL+List.sp[[i]]$b_FL.to.TL)/10)
        if(MN.SZE==0) Min.size.bin=0
        MaxLen= 10*round(List.sp[[i]]$TLmax*Plus.gp.size/10)
        lbnd = seq(Min.size.bin,MaxLen - TL.bins.cm, TL.bins.cm)
        
        
        tiff(file=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[i]]$Name),
                        "/",AssessYr,"/1_Inputs/Visualise data","/Observe size VS selectivity used in SS.tiff",sep=''),
             width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
        SEl.sens_NSF=List.sp[[i]]$SS_selectivity.sensitivity
        fun.compare.sel.obs.size.comp(TL=lbnd,
                                      SEL=List.sp[[i]]$SS_selectivity,
                                      SEl.sens_NSF=List.sp[[i]]$SS_selectivity.sensitivity,
                                      size.comps=d.list,
                                      Flts=sort(unique(d.list$Fleet)))
        dev.off()
      }
    }
  }
  clear.log('fun.compare.sel.obs.size.comp')
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
                     R=store.species.r_M.age.invariant[[l]]$mean)
      }
    }
    print(paste("Displaying CPUE correlation for -----",names(Catch.rate.series)[l]))
  }
  detach("package:reshape", unload=TRUE)
  detach("package:FLCore", unload=TRUE)
  detach("package:plyr", unload=TRUE)
  clear.log('fn.cpue.corr')
}

#Check if observed FL is within Lo +/- CV and FLinf +/- 
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
    if(any(grepl('Size_composition',names(Species.data[[i]]))))
    {
      d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                    names(Species.data[[i]]))]
      if(length(d.list)>0)
      {
        print(paste("Check if observed TL is with Lo +/- CV and FLinf --------",names(Species.data)[i]))
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
        Lo=(List.sp[[i]]$Lzero)*List.sp[[i]]$a_FL.to.TL+List.sp[[i]]$b_FL.to.TL
        Lo.min=Lo-(CV_young*Lo)
        LinfF=(List.sp[[i]]$Growth.F$FL_inf)*List.sp[[i]]$a_FL.to.TL+List.sp[[i]]$b_FL.to.TL
        LinfM=(List.sp[[i]]$Growth.M$FL_inf)*List.sp[[i]]$a_FL.to.TL+List.sp[[i]]$b_FL.to.TL
        if(is.na(LinfM)) LinfM=LinfF
        if(LinfF==LinfM) LinfM=base::jitter(LinfM)
        LinfF.max=LinfF+(CV_old*LinfF)
        LinfM.max=LinfM+(CV_old*LinfM)
        nRows=length(unique(d.list$Fleet))
        nCols=1
        d.list%>%
          mutate(TL=FL*List.sp[[i]]$a_FL.to.TL+List.sp[[i]]$b_FL.to.TL)%>%
          filter(!is.na(SEX))%>%
          ggplot(aes(TL,fill=SEX))+
          geom_histogram(binwidth = TL.bins.cm)+
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
}

# Check that meanbodywt used in SS occurs ~ at peak of selectivity
if(First.run=="YES")
{
  for(i in 1:length(Species.data))
  {
    if(any(grepl('annual.mean.size',names(Species.data[[i]]))))
    {
      print(paste("check if meanbodywt used in SS occurs ~ at peak of selectivity --------",names(Species.data)[i]))
      
       tiff(file=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[i]]$Name),
                      "/",AssessYr,"/1_Inputs/Visualise data","/Meanbodywt VS selectivity used in SS.tiff",sep=''),
           width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
       TL=with(List.sp[[i]],seq(round(Lzero*a_FL.to.TL+b_FL.to.TL),TLmax))
       TL=c(seq(0,(min(TL)-1)),TL)
       fun.check.mean.weight(TL=TL,
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

# Get sex ratio by zone
if(First.run=="YES")
{
  Dis.size.data=SS3_fleet.size.comp.used
  Sex.ratio.zone=vector('list',N.sp)
  names(Sex.ratio.zone)=Keep.species
  for(i in 1:N.sp)
  {
    id=grep(paste(Dis.size.data,collapse="|"),names(Species.data[[i]]))
    if(length(id)>0)
    {
      print(paste("Display sex ratio by zone from size comp for -----",Keep.species[i]))
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
    if(any(grepl('age_length',names(Species.data[[i]]))))
    {
      print(paste('Display conditional age at length for ------',Keep.species[i]))
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
        facet_wrap(~Sex,ncol=1)+theme_PA()+ylim(0,NA)+xlim(0,NA)
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
    print(paste('Check if r within m range ------',Keep.species[i]))
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

#Define species with size composition  
Species.with.length.comp=vector('list',N.sp)
for(i in 1:N.sp)
{
  if(any(grepl('Size_composition',names(Species.data[[i]]))))
  {
    d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[i]]))]
    if(names(Species.data)[i]%in%names(Indicator.species))
    {
      Min.size=Min.annual.obs.ktch
    }else
    {
      Min.size=Min.annual.obs.ktch*prop.min.N.accepted_other
    }
    
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
        mutate(Min.accepted.N=ifelse(!fishry=='Survey',Min.size,10))%>%
        filter(N>=Min.accepted.N)%>%
        mutate(dummy=paste(year,fishry,sex))
      
      if(nrow(Table.n)>0) Species.with.length.comp[[i]]=Keep.species[i]
      rm(Life.history)
    }
  }
} #end i loop
Species.with.length.comp=c(do.call(rbind,Species.with.length.comp)) 
#not.fitting.SS3=c("pigeye shark","lemon shark")  #no Hessian 
#Species.with.length.comp=subset(Species.with.length.comp,!Species.with.length.comp%in%not.fitting.SS3)
if(First.run=="YES") write.csv(paste('1.',capitalize(Species.with.length.comp),sep=''),paste(Rar.path,'Species.with.length.comp.csv',sep='/'),row.names = F) 

#Display total catch VS cpue
if(First.run=="YES")
{
  for(x in 1:N.sp)
  {
    Neim=Keep.species[x]
    CPUE=compact(Catch.rate.series[[x]])
    if(!is.null(CPUE))
    {
      print(paste('Display Catch and cpue for ------',Keep.species[x]))
      DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
      if(length(DROP)>0)CPUE=CPUE[-DROP]
      for(p in 1:length(CPUE)) CPUE[[p]]=CPUE[[p]]%>%mutate(Mean.st=Mean/max(Mean,na.rm=T),UP.CI=UP.CI/max(Mean,na.rm=T),LOW.CI=LOW.CI/max(Mean,na.rm=T))%>%
          dplyr::select(yr.f,Mean.st,UP.CI,LOW.CI)%>%
          mutate(fleet=names(CPUE)[p])
      CPUE=do.call(rbind,CPUE)%>%rename(Year=yr.f)
      Tot=KtCh%>%
        filter(Name==Neim)%>%
        group_by(finyear)%>%
        summarise(Total=sum(LIVEWT.c))%>%
        ungroup()%>%
        mutate(Total=max(CPUE$UP.CI,na.rm=T)*(Total/max(Total)))
      
      CPUE%>%
        ggplot(aes(Year,Mean.st,color=fleet))+geom_point(size=2)+geom_errorbar(aes(x=Year,ymin=LOW.CI,ymax=UP.CI),alpha=0.5)+geom_line()+
        geom_line(data=Tot,aes(finyear,Total),linewidth=3,color='black',alpha=0.2)+
        theme_PA()+theme(legend.position = 'top')
      HandL=handl_OneDrive("Analyses/Population dynamics/1.")
      DiR=paste(HandL,capitalize(Keep.species[x]),"/",AssessYr,"/1_Inputs/Visualise data",sep='')
      ggsave(paste(DiR,'Total catch VS cpue.tiff',sep='/'), width = 8,height = 6, dpi = 300, compression = "lzw")
    }
  }
}

#---17. Display catches by fishery,life history & available time series ----
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
  #By species groups
  for(l in 1:length(Lista.sp.outputs))
  {
    print(paste('Display catch by species for ------',names(Lista.sp.outputs)[l]))
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
  
  #All assessed species
  Tot.ktch%>%
    filter(Name%in%Keep.species)%>%
    group_by(finyear,Type,Name)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm = T))%>%
    mutate(Name=capitalize(Name))%>%
    ggplot(aes(finyear,LIVEWT.c,color=Type))+
    geom_point()+geom_line()+
    facet_wrap(~Name,scales='free_y')+
    theme_PA()+
    theme(legend.position = 'top',
          legend.title=element_blank())+
    xlab('Financial year')+ylab('Catch (tonnes)')
  ggsave(paste(Rar.path,'/Catch_assessed_species.tiff',sep=''),width = 10,height = 8,compression = "lzw")

  #Display species with length comp and available selectivity by fleet  
  Export.table.sels.man.fleet=vector('list',N.sp)
  for(i in 1:N.sp)
  {
    if(Keep.species[i]%in%Species.with.length.comp)
    {
      dis.sel=SS_selectivity_init_pars%>%filter(Species==Keep.species[i])
      xx=Tot.ktch%>%
            filter(Name%in%Keep.species[i])%>%
            group_by(Type,Name)%>%
            summarise(LIVEWT.c=sum(LIVEWT.c,na.rm = T))%>%
            arrange_(~ desc(LIVEWT.c))%>%
        ungroup()%>%
        mutate(Percent=LIVEWT.c/sum(LIVEWT.c),
               CumPercent=cumsum(Percent))%>%
        filter(CumPercent<=.9)%>%
        mutate(Type=as.character(Type),
               Type=case_when(Type%in%c('NSF','Indonesia','Taiwan')~'NSF',
                              Type=='TDGDLF'~'TDGDLF',
                              Type=='WRL'~'WRL',
                              TRUE  ~ "Other"))
      
      ppp=dis.sel%>%dplyr::select(paste0("Source_",xx$Type))%>%
        gather(Type,Selectivity)%>%
        mutate(Type=sub('.*_', '',Type))
      Export.table.sels.man.fleet[[i]]=xx%>%left_join(ppp,by='Type')
    }
  }
  write.csv(do.call(rbind,compact(Export.table.sels.man.fleet)),
            paste(Rar.path,'Table_SS selectivities by fleets and prop of total catch.csv',sep='/'),row.names = F)  
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
      eig=First.Age:Max.age.F[2]
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
        annotate("text", x = mean(Max.age.F), y = 0,parse = T, label = as.character(Amx),color='darkorange4')+
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
              legend.position="top")+
        expand_limits(x = 0, y = 0)
      
      plt[[2]]=data.frame(TL=lengz,
                          Maturity=1/(1+exp(-log(19)*((lengz-TL.50.mat)/(TL.95.mat-TL.50.mat)))))%>%
        ggplot(aes(TL,Maturity))+
        geom_line(size=1.25,color="#F8766D")+xlab("TL (cm)")+ylab('Proportion mature')+
        theme_PA()+
        annotate("text", x = TL.50.mat*.75, y = .8, label = paste('Fecundity [',Fecundity[1],',',Fecundity[2],']',sep=''))+
        annotate("text", x = TL.50.mat*.75, y = .85, label = paste('Cycle [',Breed.cycle[1],',',Breed.cycle[2],']',sep=''))+
        expand_limits(x = 0, y = 0)
      
      plt[[3]]=rbind(
        data.frame(Sex='F',
                   TL=lengz,
                   Twt=AwT*lengz^BwT),
        data.frame(Sex='M',
                   TL=lengz,
                   Twt=AwT.M*lengz^BwT.M))%>%
        ggplot(aes(TL,Twt,color=Sex))+
        geom_line(size=1.25)+xlab("TL (cm)")+ylab('Twt (kg)')+
        theme_PA()+
        expand_limits(x = 0, y = 0) 
      
      plt[[4]]=data.frame(Age=eig,
                          M=Mmean.mean.at.age)%>%
        ggplot(aes(Age,M))+
        geom_line(size=1.25,color="black")+xlab("Age")+ylab(expression(paste(plain("Natural mortality (year") ^ plain("-1"),")",sep="")))+
        theme_PA()+
        expand_limits(x = 0, y = 0)
      
      plt[[5]]=data.frame(TL=lengz,
                          FL=(lengz-b_FL.to.TL)/a_FL.to.TL)%>%
        ggplot(aes(FL,TL))+
        geom_line(size=1.25)+xlab("FL (cm)")+ylab('TL (cm)')+
        theme_PA()+
        annotate("text", x = mean(lengz), y = mean(lengz)*.8,
                 label = paste('TL =',a_FL.to.TL,' FL + ',b_FL.to.TL ,sep=''))+
        expand_limits(x = 0, y = 0)

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
if(First.run=="YES")
{
  write.csv(data.frame(year.current=Last.yr.ktch,
                       year.future=paste(fut.yr,substr(fut.yr+1,3,4),sep='-')),
            paste(Rar.path,'year_current_future.csv',sep='/'),row.names = FALSE)
}

#Display available time series by species
if(First.run=='YES')
{
  Data.availability=vector('list',N.sp)
  for(i in 1:N.sp)
  {
    dd=data.frame(Species=rep(capitalize(names(Species.data)[i]),each=6),
                  Data=c('Catch','Abundance','Catch length\n composition','Catch mean\n weight','F','Tagging'),
                  Availability=NA,
                  w=1)
    Avail.ktch=1
    Avail.abun=Avail.length=Avail.ca.wei=Avail.F=Avail.Tag=0
    if(!is.null(Catch.rate.series[[i]])) Avail.abun=1
    if(Keep.species[i]%in%Species.with.length.comp) Avail.length=1
    if('annual.mean.size'%in%names(Species.data[[i]])) Avail.ca.wei=1
    if(any(grepl('Fishing.mortality',names(Species.data[[i]])))) Avail.F=1
    if(any(grepl('tag',names(Species.data[[i]])))) Avail.Tag=1
    dd$Availability=c(Avail.ktch,Avail.abun,Avail.length,Avail.ca.wei,Avail.F,Avail.Tag)
    Data.availability[[i]]=dd
  }
  Data.availability=do.call(rbind,Data.availability)%>%
    mutate(Availability=factor(Availability,levels=0:1),
           Data=factor(Data,levels=c('Catch','Abundance','Catch length\n composition',
                                     'Catch mean\n weight','F','Tagging')))
  Avail.colors=c('transparent','black')
  names(Avail.colors)=0:1
  Data.availability%>%
    ggplot(aes(Data,Species,  width = w,height=w))+
    geom_tile(aes(fill = Availability))+
    scale_fill_manual(values=Avail.colors)+
    theme_bw()%+replace% 
    theme(legend.position = "none",
          panel.background = element_rect(fill="grey95"),
          panel.border = element_rect(colour = "grey75", fill=NA, size=1.15),
          panel.grid.major = element_blank(),    
          panel.grid.minor = element_blank(),    
          axis.line = element_line(colour = "grey75"),
          strip.background = element_rect(fill = "transparent",colour = "transparent"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 13),
          axis.text.y=element_text(hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())+
    xlab('Time series')+ylab('')+
    scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0),position = "bottom")+
    geom_hline(yintercept=seq(1.5,N.sp,1),color='grey75')+
    geom_vline(xintercept=seq(1.5,5.5,1),color='grey75')
  ggsave(paste(Rar.path,"Available time series by species.tiff",sep='/'),
         width = 10,height = 8,compression = "lzw")
}

#Display available length composition and selectivity for main fleet
do.dis=FALSE
if(do.dis)
{
  Klfunc <- colorRampPalette(c("yellow", "red"))
  pdf(handl_OneDrive("Analyses/Population dynamics/Length.comps_main_fleet.pdf"))
  for(l in 1:N.sp)
  {
    this.fleet=Katch%>%filter(Name==Keep.species[l])
    this.size.comp=this.fleet%>%pull(this.size.comp)
    this.fleet=this.fleet%>%pull(FishCubeCode)
    iid=Species.data[[l]][fn.extract.dat(this.size.comp,names(Species.data[[l]]))]
    iid=iid[fn.extract.dat('Size_composition',names(iid))]
    if(any(grepl('Observations',names(iid)))) iid=iid[-grep('Observations',names(iid))]
    if(any(grepl('Table',names(iid)))) iid=iid[-grep('Table',names(iid))]
    if(length(iid)>0)
    {
      for(x in 1:length(iid))
      {
        dd=str_before_first(str_after_nth(names(iid)[x],"_",2), coll(".inch"))
        iid[[x]]=iid[[x]]%>%
          mutate(Zone=str_before_first(dd, coll(".")),  
                 Mesh=str_after_first(dd, coll(".")))
      }
      dummy=do.call(rbind,iid)%>%
        filter(year<=as.numeric(substr(Last.yr.ktch,1,4)))
      dummy=dummy%>%
        mutate(TL=FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)
      if(this.size.comp%in%c("West","Zone") & unique(Selectivity.at.totalength[[l]]$type)=='Species')
      {
        SS.flit='Southern.shark_1'
        MaxLen = 10*round(List.sp[[l]]$TLmax/10)
        min.TL=with(List.sp[[l]],Lzero*a_FL.to.TL+b_FL.to.TL)
        lbnd = seq(0,MaxLen - LenInc, LenInc)
        ubnd = lbnd + LenInc
        midpt = lbnd + (LenInc/2)
        SelectivityVec_full=with(List.sp[[l]]$SS_selectivity%>%filter(Fleet==SS.flit),doubleNorm24.fn(midpt,a=P_1,
                                                                                                      b=P_2, c=P_3, d=P_4, e=P_5, f=P_6,use_e_999=FALSE, use_f_999=FALSE))
      }else
      {
        SelectivityVec_full=NA
      }
      if(any(is.na(dummy$Zone)))
      {
        N.min=dummy%>%
          group_by(FINYEAR)%>%
          tally()
        p1=dummy%>%
          mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2)%>%
          group_by(FINYEAR,bin)%>%
          tally()%>%
          group_by(FINYEAR)%>%
          mutate(n1=n/max(n, na.rm=TRUE))%>%
          left_join(N.min%>%rename(N=n)%>%dplyr::select(FINYEAR,N),by='FINYEAR')%>%
          mutate(FINYEAR=paste0(FINYEAR," (n=",N,")"))%>%
          ggplot(aes(bin,n1))+
          geom_bar(stat="identity")+
          facet_wrap(~FINYEAR,scales='free_y')+
          ggtitle(paste(this.fleet,Keep.species[l],sep=' --- '))+
          xlab("Total length (cm)")+ylab("")
        Strip.Size=12
      }
      if(!any(is.na(dummy$Zone)))
      {
        Min.N=Min.annual.obs.ktch*prop.min.N.accepted_other
        if(Keep.species[l]%in%names(Indicator.species)) Min.N=Min.annual.obs.ktch*.5
        
        N.min=dummy%>% 
          group_by(FINYEAR,Zone,Mesh)%>%
          tally()
        p1=dummy%>%
          mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2)%>%
          group_by(FINYEAR,Zone,Mesh,bin)%>%
          tally()%>%
          group_by(FINYEAR,Zone,Mesh)%>%
          mutate(n1=n/max(n, na.rm=TRUE))%>%
          left_join(N.min%>%rename(N=n)%>%dplyr::select(FINYEAR,Zone,Mesh,N),by=c('FINYEAR','Zone','Mesh'))%>%
          mutate(FINYEAR=paste0(FINYEAR," (n=",N,")"),
                 Mesh.zone=paste0(Zone,"_",Mesh))%>%
          filter(N>=Min.N)
        N.cols=length(unique(p1$FINYEAR))
        if(nrow(p1)>0)
        {
          p1=p1%>%
            ggplot(aes(bin,n1,color=FINYEAR))+
            geom_line(linewidth=1.05)+
            facet_wrap(~Mesh.zone,scales='free_y')+
            scale_color_manual(values=Klfunc(N.cols))
          if(any(!is.na(SelectivityVec_full)))
          {
            p1=p1+
              geom_point(data=data.frame(midpt=midpt,SelectivityVec=SelectivityVec_full),
                         aes(midpt,SelectivityVec),color='blue')
            
          }
          
          p1=p1+
            ggtitle(paste(this.fleet,Keep.species[l],sep=' --- '))+
            xlab("Total length (cm)")+ylab("")+
            theme(legend.position="bottom",
                  legend.title = element_blank(),
                  legend.text=element_text(size=8),
                  strip.text.x = element_text(size = Strip.Size),
                  axis.text=element_text(size=10),
                  axis.title=element_text(size=12))
          
        }else
        {
          p1=ggplot() + theme_void()+ ggtitle(paste(this.fleet,'--',Keep.species[l],' (length comp sample size <',Min.annual.obs.ktch,')',sep=' '))
        }

        Strip.Size=12
        

        dd=dummy%>%
          mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2)%>%
          group_by(FINYEAR,Zone,Mesh,bin)%>%
          tally()%>%
          group_by(FINYEAR,Zone,Mesh)%>%
          mutate(n1=n/max(n, na.rm=TRUE))%>%
          left_join(N.min%>%rename(N=n)%>%dplyr::select(FINYEAR,Zone,Mesh,N),by=c('FINYEAR','Zone','Mesh'))%>%
          mutate(Yr.Mesh.zone=paste0(FINYEAR,Zone,Mesh),
                 Mesh.zone=paste0(Zone,"_",Mesh))%>%
          filter(N>=Min.N)
        
        p2=dummy%>%
          mutate(Yr.Mesh.zone=paste0(FINYEAR,Zone,Mesh),
                 yr.zone=paste(FINYEAR,Zone))%>%
          filter(Yr.Mesh.zone%in%unique(dd$Yr.Mesh.zone))%>%
          ggplot() +
          geom_density_ridges(aes(x = TL, y = yr.zone,fill = Mesh),scale = 1.5,jittered_points = TRUE,
                              position = position_points_jitter(width = 0.05, height = 0),
                              point_shape = '|', point_size = 3, point_alpha = 1,alpha = 0.7, rel_min_height = 0.01) +
          theme(legend.position = "top")
      }
     }
    if(length(iid)==0 | this.size.comp=="no length comp data") 
    {
      p1=ggplot() + theme_void()+ ggtitle(paste(this.fleet,Keep.species[l],sep=' --- ')) 
    }
    print(p1)
    rm(p1)
    if(exists('p2'))
    {
      print(p2)
      rm(p2)
    }
  }
  dev.off()
}

#Define length composition main mesh and zone
if(do.dis)
{
  STORE=vector('list',N.sp)
  names(STORE)=Keep.species
  for(l in 1: N.sp)
  {
    Neim=Keep.species[l]
    SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
    
    #identify fishery with most catch 
    Katch1=KtCh%>%
      filter(Name==Neim)%>%
      mutate(Fishery=case_when(FishCubeCode%in%c('Historic',Southern.shark.fisheries)~"TDGDLF",
                               FishCubeCode%in%Northern.shark.fisheries~"NSF",
                               .default="Other"))
    Katch=Katch1%>%
      group_by(Fishery)%>%
      summarise(Tonnes=sum(LIVEWT.c))
    this.size.comp=Katch%>%filter(Tonnes==max(Tonnes))%>%pull(Fishery)
    if(Neim%in%Other.to.NSF) this.size.comp='NSF'
    if(Neim%in%Other.to.TDGDLF) this.size.comp='TDGDLF'
    out.fleet=this.size.comp
    if(any(this.size.comp=='TDGDLF'))
    {
      this.size.comp=paste('Size_composition',c('West.6.5','West.7','Zone1.6.5','Zone1.7','Zone2.6.5','Zone2.7'),sep="_")
    }
    if(any(this.size.comp=='NSF'))
    {
      this.size.comp=paste('Size_composition',c('NSF.LONGLINE'),sep="_")
    }
    
    #get relevant length composition  
    iid=Species.data[[l]][fn.extract.dat(this.size.comp,names(Species.data[[l]]))]
    if(any(grepl('Observations',names(iid)))) iid=iid[-grep('Observations',names(iid))]
    if(length(iid)>0 & !Neim%in%no.empirical.sel.main.fleet)
    {
      for(x in 1:length(iid))
      {
        dd=str_before_first(str_after_nth(names(iid)[x],"_",2), coll(".inch"))
        iid[[x]]=iid[[x]]%>%
          mutate(Zone=str_before_first(dd, coll(".")),  
                 Mesh=str_after_first(dd, coll(".")))
      }
      dummy=do.call(rbind,iid)%>%
        filter(year<=as.numeric(substr(Last.yr.ktch,1,4)))
      STORE[[l]]=list(Table1=with(dummy,table(year,Zone)),Table2=with(dummy,table(year,Mesh,Zone)))
      
    }
  }
}



#---18. Catch-only assessments --------------------------------------
#Decide if using Catch-only methods for species with only catch and life history data OR all species
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

#Execute Catch only function  
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

#Get population numbers from K estimation (virgin conditions)
Get.pop.numbers=FALSE
if(Get.pop.numbers)
{
  Numbers.from.K=vector('list',length(Lista.sp.outputs))
  for(l in 1: length(Lista.sp.outputs))
  {
    dis.spicis=Lista.sp.outputs[[l]]
    dummy=vector('list',length(dis.spicis))
    K.estims_DBSRA=read.csv(paste0(Rar.path,'/Table 3. Catch only_estimates_DBSRA_',names(Lista.sp.outputs)[l],'.csv'))%>%
      filter(Scenario=='S1' & Parameter=='K')
    K.estims_CMSY=read.csv(paste0(Rar.path,'/Table 3. Catch only_estimates_CMSY_',names(Lista.sp.outputs)[l],'.csv'))%>%
      filter(Scenario=='S1' & Parameter=='K')
    for(i in 1:length(dis.spicis))
    {
      Carrying.capacity=rbind(K.estims_DBSRA%>%filter(Species==capitalize(dis.spicis[i])),
                              K.estims_CMSY%>%filter(Species==capitalize(dis.spicis[i])))
      dummy[[i]]=fn.get.Numbers.from.K(Life.history=List.sp[[match(dis.spicis[i],names(List.sp))]],K=mean(Carrying.capacity$Median)*1000)%>%
        mutate(Species=capitalize(dis.spicis[i]),
               N=sum(N.at.age))
      
      
    }
    Numbers.from.K[[l]]=do.call(rbind,dummy)
  }
  do.call(rbind,Numbers.from.K)%>%
    mutate(Species=case_when(Species=='Great hammerhead'~'Great H.',
                             Species=='Scalloped hammerhead'~'Scalloped H.',
                             Species=='Smooth hammerhead'~'Smooth H.',
                             TRUE~Species),
           Facet=paste0(Species,' (N=',format(round(N), nsmall=0, big.mark=","),')'))%>%
    ggplot(aes(Age,N.at.age))+
    geom_point()+
    facet_wrap(~Facet,scales='free')+
    ylab('Numbers')+
    theme_PA(strx.siz=9)
  ggsave(paste0(Rar.path,"/Numbers of individuals from K_Catch only.tiff"),
         width = 11.5,height = 8, dpi = 300, compression = "lzw")
  
}


#Get Consequence and likelihoods
if(COM_use.this.for.risk=='catch')
{
  if (any(grepl("Table 4. Catch.only_catch_vs_MSY",list.files(Rar.path))))
  {
    DD=vector('list',length(Lista.sp.outputs))
    for(l in 1:length(Lista.sp.outputs))  
    {
      ddis=names(Lista.sp.outputs)[l]
      DD[[l]]=read.csv(paste(Rar.path,'/Table 4. Catch.only_catch_vs_MSY_',ddis,'.csv',sep=''))%>%
        dplyr::select(Species,P.catch.below.tar,P.between.thre.tar,P.between.lim.thre,P.catch.above.limit)
    }
    Store.cons.Like_COM=do.call(rbind,DD)
  }
}
if(COM_use.this.for.risk=='biomass')
{
  if (any(grepl("Table 4. Catch only_current.B.over.Bmsy",list.files(Rar.path))))
  {
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
  }
}
  
#Clear log
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
clear.log('fn.get.Numbers.from.K')

#---19. Spatio-temporal catch and effort. Reported TDGLF and NSF ----   
if(Do.Spatio.temporal.catch.effort) fn.source1('TDGLF and NSF spatio-temporal catch and effort.r')


#---20. Changes in observed mean length ----
if(Do.Changes.in.observed.size) fn.source1("Changes_mean_length.r")

#---21. Changes in reported mean weight of individuals caught in the TDGDLF  ----
#Is there a strong declining trend in mean weights? (Leitao 2019)
if(Do.Changes.in.reported.size) fn.source1("Changes_reported_size.r")
clear.log('Logbook')
clear.log('Logbook.sp')


#---22. Catch curve with dome-shaped selectivity --------------------------------------
#note: for each species, identify main fleet and a good sampling period (i.e. fishery with highest catch and length
#       composition from exploited period to be able to determine F) and combine 2/3 years.
if(do.Size.based.Catch.curve)
{
  #Run catch curve and YPR
  fn.source1("Apply_Length.catch.curve.r")
  
  #Get Consequence and likelihoods
  if(any(grepl("Table 6. Catch.curve_YPR",list.files(Rar.path))))
  {
    Store.cons.Like_CatchCurve=list(Depletion=NULL,B.over.Bmsy=NULL)     
    DD.depletion=DD.B.over.Bmsy=vector('list',length(Lista.sp.outputs))
    for(l in 1:length(Lista.sp.outputs))  
    {
      
    }
    Store.cons.Like_CatchCurve$Depletion=do.call(rbind,DD.depletion)
    Store.cons.Like_CatchCurve$B.over.Bmsy=do.call(rbind,DD.B.over.Bmsy)  
  }
}
   

#---23. Dynamic catch and size composition with dome-shaped selectivity --------------------------------------
#note: not used due to Cpp implementation issues
if(do.Dynamic.catch.size.comp) fn.source1("Apply_Dynamic.catch.size.comp.r") 
clear.log('ObjFunc')

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

  #Run JABBA
if(Do.StateSpaceSPM) fn.source1("Apply_StateSpaceSPM.R")   #takes 0.9 hours (3e4 sims for 8 species with 6 scenarios and fit diagnostics)

  #Get Consequence and likelihoods
if (any(grepl("Table 9. JABBA CPUE_Current.B.over.Bmsy",list.files(Rar.path))))
{
  Store.cons.Like_JABBA=list(Depletion=NULL,B.over.Bmsy=NULL)     
  DD.depletion=DD.B.over.Bmsy=DD.f=vector('list',length(Lista.sp.outputs))
  for(l in 1:length(Lista.sp.outputs))  
  {
    ddis=names(Lista.sp.outputs)[l]
    DD.B.over.Bmsy[[l]]=read.csv(paste(Rar.path,'/Table 9. JABBA CPUE_Current.B.over.Bmsy_',ddis,'.csv',sep=''))%>%
      dplyr::select(-c(Scenario,Model))%>%
      gather(Species,Probability,-c(Range,finyear))%>%
      mutate(Species=gsub(".", " ", Species, fixed=TRUE))
    
    DD.depletion[[l]]=read.csv(paste(Rar.path,'/Table 9. JABBA CPUE_Current.depletion_',ddis,'.csv',sep=''))%>%
      dplyr::select(-c(Scenario,Model))%>%
      gather(Species,Probability,-c(Range,finyear))%>%
      mutate(Species=gsub(".", " ", Species, fixed=TRUE))
    
    DD.f[[l]]=read.csv(paste(Rar.path,'/Table 9. JABBA CPUE_Current.f_',ddis,'.csv',sep=''))%>%
      dplyr::select(-c(Scenario,Model))%>%
      gather(Species,Probability,-c(Range,finyear))%>%
      mutate(Species=gsub(".", " ", Species, fixed=TRUE))
  }
  Store.cons.Like_JABBA$Depletion=do.call(rbind,DD.depletion)
  Store.cons.Like_JABBA$B.over.Bmsy=do.call(rbind,DD.B.over.Bmsy)
  Store.cons.Like_JABBA$f=do.call(rbind,DD.f)
}
  
  #Clear log
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
  #Run Stock Synthesis
HandL.out=handl_OneDrive("Analyses/Population dynamics/1.")
HandL.out.RAR=Rar.path
if(Do.integrated) fn.source1("Apply_SS.R")   #Takes ~ 10 hours

  #Get Consequence and likelihoods 
if(any(grepl("Table 12. Age.based_SS_current.depletion",list.files(Rar.path))))
{
  Store.cons.Like_Age.based=list(Depletion=NULL,B.over.Bmsy=NULL)     
  DD.depletion=DD.B.over.Bmsy=DD.f=vector('list',length(Lista.sp.outputs))
  for(l in 1:length(Lista.sp.outputs))  
  {
    ddis=names(Lista.sp.outputs)[l]
    DD.B.over.Bmsy[[l]]=read.csv(paste(Rar.path,'/Table 12. Age.based_SS_current.B.over.Bmsy_',ddis,'.csv',sep=''))%>%
      dplyr::select(-c(Scenario,Model))%>%
      gather(Species,Probability,-c(Range,finyear))%>%
      mutate(Species=gsub(".", " ", Species, fixed=TRUE))
    
    DD.depletion[[l]]=read.csv(paste(Rar.path,'/Table 12. Age.based_SS_current.depletion_',ddis,'.csv',sep=''))%>%
      dplyr::select(-c(Scenario,Model))%>%
      gather(Species,Probability,-c(Range,finyear))%>%
      mutate(Species=gsub(".", " ", Species, fixed=TRUE))
    
    DD.f[[l]]=read.csv(paste(Rar.path,'/Table 12. Age.based_SS_Current.f_',ddis,'.csv',sep=''))%>%
      dplyr::select(-c(Scenario,Model))%>%
      gather(Species,Probability,-c(Range,finyear))%>%
      mutate(Species=gsub(".", " ", Species, fixed=TRUE))
  }
  Store.cons.Like_Age.based$Depletion=do.call(rbind,DD.depletion)
  Store.cons.Like_Age.based$B.over.Bmsy=do.call(rbind,DD.B.over.Bmsy)
  Store.cons.Like_Age.based$f=do.call(rbind,DD.f)
}

  #Clear log
clear.log('Store.sens')
clear.log('output')
clear.log('F.Fmsy')
clear.log('out')
clear.log('fn.ktch.cpue')     
clear.log('Post')
clear.log('dummy.store.Kobe.probs')


  #Create SS-DL tool input files
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
if(Do.bespoke) fn.source1("Integrated_size_based.R")


#---27. Weight of Evidence Assessment ------

#0. Stock assessment flow chart
Ass.flow=TRUE
if(Ass.flow)
{
  library(igraph)
  Dat1 <- tibble(from = c("Shark and ray species\n interacting with\n fisheries in WA",
                          "PSA", "PSA", "PSA",
                          "High\n vulnerability",
                          "COMs", "COMs", "COMs",
                          "Abundance\n index","Abundance\n index","Size\n composition","Size\n composition"),
                 to = c("PSA",
                        "Low\n vulnerability", "Medium\n vulnerability", "High\n vulnerability",
                        "COMs",
                        "No abundance\n or\n size composition", "Abundance\n index", "Size\n composition",
                        "JABBA","SS3","SS3","SS3.size"))
  g = graph_from_data_frame(Dat1, directed = TRUE)
  coords = layout_as_tree(g)
  colnames(coords) = c("x", "y")
  output_df = as_tibble(coords) %>%
    mutate(step = vertex_attr(g, "name"),
           label = gsub("\\d+$", "", step),
           x = x*-1,
           x=case_when(label=="Low\n vulnerability"~-1,
                       label=="Medium\n vulnerability"~0,
                       label=="Size\n composition"~2,
                       label=="No abundance\n or\n size composition"~-0.5,
                       label=="Abundance\n index"~0.85,
                       label=="JABBA"~0.75,
                       label=="SS"~1.625,
                       label=="SS3.size"~2.5,
                       TRUE~x),
           Ass.Level = case_when(grepl("SS",label)~'5',
                                 label=="JABBA"~'4',
                                 label=="COMs"~'1',
                                 TRUE~'0'),
           label=ifelse(label=="SS3.size","SS (size only)",label))   
  plot_nodes = output_df %>%
    mutate(xmin = x - 0.35,
           xmax = x + 0.35,
           ymin = y - 0.25,
           ymax = y + 0.25,
           xmin=case_when(grepl(paste(c("Shark and ray","No abundance"),collapse='|'),label)~x - 0.8,
                          grepl(paste(c("vulnerability","Abundance","composition"),collapse='|'),label)~x - 0.425,
                          TRUE~xmin),
           xmax=case_when(grepl(paste(c("Shark and ray","No abundance"),collapse='|'),label)~x + 0.8,
                          grepl(paste(c("vulnerability","Abundance","composition"),collapse='|'),label)~x + 0.425,
                          TRUE~xmax),
           ymin=case_when(grepl("Shark and ray",label)~y - 0.4,
                          label=="No abundance\n or\n size composition"~y - 0.4,
                          TRUE~ymin),
           ymax=case_when(grepl("Shark and ray",label)~y + 0.4,
                          label=="No abundance\n or\n size composition"~y + 0.4,
                          TRUE~ymax))
  
  
  plot_edges = Dat1 %>%
    mutate(id = row_number()) %>%
    pivot_longer(cols = c("from", "to"),
                 names_to = "s_e",
                 values_to = "step") %>%
    left_join(plot_nodes, by = "step") %>%
    select(-c(label, Ass.Level, y, xmin, xmax)) %>%
    mutate(y = ifelse(s_e == "from", ymin, ymax)) %>%
    select(-c(ymin, ymax))
  
  N=nrow(read.csv(paste(Exprt,'Table S1_All.species.caught.by.fishery.csv',sep='/')))
  N.COMS=N.sp
  JABBA.species=unique(Store.cons.Like_JABBA$Depletion$Species)
  N.JABBA=length(JABBA.species)
  SS.species=unique(Store.cons.Like_Age.based$Depletion$Species)
  N.SS=length(SS.species)
  SS.species.only.length=SS.species[which(!SS.species%in%JABBA.species)]
  N.SS.size.only=length(SS.species.only.length)
  N.SS=N.SS-N.SS.size.only
  
  plot_nodes=plot_nodes%>%
    mutate(label=case_when(grepl('Shark and ray',label)~paste0("Shark and ray species\n interacting with\n fisheries in WA (n= ",N,')'),
                           label=='COMs'~paste0('COMs\n (n= ',N.COMS,')'),
                           label=='JABBA'~paste0('JABBA\n (n= ',N.JABBA,')'),
                           label=='SS'~paste0('SS\n (n= ',N.SS,')'),
                           label=='SS (size only)'~paste0('SS\n (n= ',N.SS.size.only,')'),
                           TRUE~label))
  VulScore=plot_nodes%>%
    filter(grepl('vulnerability',label))%>%
    mutate(y=y+.5,
           label=ifelse(grepl('High',label),paste0('>',medium.risk),
                        ifelse(grepl('Medium',label), paste(Low.risk,'to',medium.risk),
                               ifelse(grepl('Low',label),paste0('<',Low.risk),
                                      NA))))
  ggplot() +
    geom_rect(data = plot_nodes,mapping = aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, 
                                              fill = Ass.Level),alpha = 0.5)  + 
    geom_text(data = plot_nodes, mapping = aes(x = x, y = y, label = label),color = "black") + 
    geom_path(data = plot_edges,mapping = aes(x = x, y = y, group = id),color = "grey50",
              arrow = arrow(length = unit(0, "cm"), type = "closed"))+ 
    geom_text(data = VulScore, mapping = aes(x = x, y = y, label = label),color = "black")+
    labs(title = "Stock Assessment Decision Tree",
         caption = " PSA, Productivity and Susceptibility analysis; COMs,Catch only Models\n\ JABBA, Just Another Bayesian Biomass Assessment; SS, Stock Synthesis")+
    guides(fill=guide_legend(title="Assessment level"))+
    theme_void()+
    theme(legend.position = 'bottom',
          plot.caption = element_text(hjust = 1))+
    scale_fill_manual(values = c('0'='transparent','1'='dodgerblue4','4'='forestgreen','5'='darkorange3'),limits = c('1', '4','5'))
  
  ggsave(handl_OneDrive("Analyses/Population dynamics/Stock assessment flow chart.tiff"),compression = "lzw")        
} 

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
Risk.PSA$LoE='PSA'

#1.2 Spatial distribution of catch and effort
Risk.Spatial=fun.risk.spatial.dist(d=data.frame(Species=c("angel sharks","copper shark","dusky shark","hammerheads","grey nurse shark",
                                                          "gummy shark","lemon shark","sandbar shark","sawsharks","shortfin mako",
                                                          "spinner shark","spurdogs","tiger shark","whiskery shark","wobbegongs"))%>%
                                     mutate(Market=case_when(Species%in%c("grey nurse shark")~'zero',
                                                             Species%in%c("spurdogs","sawsharks","angel sharks")~'decreasing',
                                                             Species%in%c("sandbar shark","copper shark","hammerheads",
                                                                          "spinner shark","tiger shark","wobbegongs")~'increasing',
                                                             TRUE~'stable'),
                                            Conservation=case_when(Species=='grey nurse shark'~'yes',
                                                                   TRUE~'no'),
                                            Blocks.fished=case_when(Species%in%c('grey nurse shark',"tiger shark","sandbar shark")~'decreasing',
                                                                    TRUE~'stable'),
                                            Blocks.fished.with.catch=case_when(Species%in%c('grey nurse shark',"spurdogs","sawsharks","angel sharks")~'decreasing',
                                                                               Species%in%c("sandbar shark","copper shark","hammerheads",
                                                                                            "spinner shark","tiger shark","wobbegongs")~'increasing',
                                                                               TRUE~'stable'),
                                            Effort='decreasing',
                                            Catch='decreasing'),
                                   outpath=Rar.path)
Risk.Spatial=Risk.Spatial%>%
              mutate(LoE='Spatial',
                     Species=ifelse(Species=="Hammerheads","Smooth hammerhead",Species),
                     Likelihood=ifelse(Species%in%c('Angel sharks','Copper shark') & Consequence==4,2,
                                       Likelihood))   #increase risk due to unknown trends for GAB Trawl and SA MSF



  #1.2. COMS  
if(exists("Store.cons.Like_COM"))
{
  if(COM_use.this.for.risk=='catch')
  {
    Risk.COM=fun.risk.CoMs(d=Store.cons.Like_COM,Ktch.type='dist')
  }
  if(COM_use.this.for.risk=='biomass')
  {
    if(Choose.probability=="Depletion")   Risk.COM=fn.risk(d=Store.cons.Like_COM$Depletion,w=1)
    if(Choose.probability=="B.over.Bmsy") Risk.COM=fn.risk(d=Store.cons.Like_COM$B.over.Bmsy,w=1)
  }
  Risk.COM$LoE='COM'
}

#1.3. Catch curve
if(exists("Store.cons.Like_CatchCurve"))
{
  if(Choose.probability=="Depletion")   Risk.CatchCurve=fn.risk(d=Store.cons.Like_CatchCurve$Depletion,w=1)
  if(Choose.probability=="B.over.Bmsy") Risk.CatchCurve=fn.risk(d=Store.cons.Like_CatchCurve$B.over.Bmsy,w=1)
  Risk.CatchCurve$LoE='CatchCurve'
}

  #1.4. JABBA
if(exists("Store.cons.Like_JABBA"))
{
  #Biomass based
  if(Choose.probability=="Depletion")   Risk.JABBA=fn.risk(d=Store.cons.Like_JABBA$Depletion,w=1)
  if(Choose.probability=="B.over.Bmsy") Risk.JABBA=fn.risk(d=Store.cons.Like_JABBA$B.over.Bmsy,w=1)
  Risk.JABBA$LoE='JABBA'
  
  #F based
  Risk.JABBA_f=fn.risk_f(d=Store.cons.Like_JABBA$f,w=1)  
  Risk.JABBA_f$LoE='JABBA'
}
do.this=FALSE
if(do.this) fn.compare.risk.ref.point(Depletion=fn.risk(d=Store.cons.Like_JABBA$Depletion,w=1)%>%filter(finyear==2021),
                               B.over.Bmsy=fn.risk(d=Store.cons.Like_JABBA$B.over.Bmsy,w=1)%>%filter(finyear==2021),
                               TITL='JABBA')

  #1.5. SS3
if(exists("Store.cons.Like_Age.based"))
{
  #Biomass based
  if(Choose.probability=="Depletion")   Risk.integrated=fn.risk(d=Store.cons.Like_Age.based$Depletion,w=1)
  if(Choose.probability=="B.over.Bmsy") Risk.integrated=fn.risk(d=Store.cons.Like_Age.based$B.over.Bmsy,w=1)
  Risk.integrated$LoE='integrated'
  
  #F based
  Risk.integrated_f=fn.risk_f(d=Store.cons.Like_Age.based$f,w=1)  
  Risk.integrated_f$LoE='integrated'
}


#2. Overall risk   

  #2.1 extract future risk   
Risk.JABBA_future=Risk.JABBA%>%filter(finyear==max(Risk.JABBA$finyear))
Risk.integrated_future=Risk.integrated%>%filter(finyear==max(Risk.integrated$finyear))
Risk.JABBA=Risk.JABBA%>%filter(finyear==min(Risk.JABBA$finyear))%>%dplyr::select(-finyear)
Risk.integrated=Risk.integrated%>%filter(finyear==min(Risk.integrated$finyear))%>%dplyr::select(-finyear)

  #2.2. Combine Risks from all individual LoEs 
Risk.COM=Risk.COM%>%
  relocate(names(Risk.JABBA))
LOE.risks=list(PSA=Risk.PSA,Spatial=Risk.Spatial,COM=Risk.COM,JABBA=Risk.JABBA,integrated=Risk.integrated)
if(exists("Store.cons.Like_CatchCurve")) LOE.risks$CatchCurve=Risk.CatchCurve 
Store.risks=LOE.risks
for(r in 1:length(LOE.risks))
{
  Store.risks[[r]]=fn.risk.figure(d=LOE.risks[[r]], Risk.colors=RiskColors, out.plot=FALSE)%>%
    dplyr::select(Species,LoE,Consequence,Likelihood,Risk,Score)
}
Store.risks=do.call(rbind,Store.risks)
Table.risks=Store.risks%>%
  mutate(Risk=paste0(Risk, ' (', Consequence,'x',Likelihood,')'))%>%
  dplyr::select(-c(Consequence,Likelihood,Score))%>%
  mutate(LoE=factor(LoE,levels=names(LOE.risks)))%>%
  spread(LoE,Risk)%>%
  rename(Integrated=integrated)

  #2.3. Calculate and display overall risk 
Weighted.overall.risk=do.call(rbind,LOE.risks[-match('PSA',names(LOE.risks))])%>%
  left_join(data.frame(LoE=names(LoE.Weights),
                       LoE.weight=LoE.Weights),by='LoE')%>%
  mutate(LoE.weight=ifelse(LoE=='PSA',1,LoE.weight))
  
  #if no integrated assessment, use next highest assessment
  en.spi=unique(Weighted.overall.risk$Species)
  for(x in 1:length(en.spi))
  {
    a=Weighted.overall.risk%>%filter(Species==en.spi[x])
    Weighted.overall.risk=Weighted.overall.risk%>%filter(!Species==en.spi[x])
    disLoes=unique(a$LoE)
    if(!"integrated"%in%disLoes)
    {
      if("JABBA"%in%disLoes)
      {
        a=a%>%
          mutate(LoE.weight=ifelse(LoE=="JABBA",1,LoE.weight))
      }else
      {
        if("CatchCurve"%in%disLoes)
        {
          a=a%>%
            mutate(LoE.weight=ifelse(LoE=="CatchCurve",1,LoE.weight))
        }else
        {
          a=a%>%
            mutate(LoE.weight=ifelse(LoE%in%c("Spatial","COM"),1,LoE.weight))
        }
      }
    }
    Weighted.overall.risk=rbind(Weighted.overall.risk,a)
  }
  
Weighted.overall.risk=Weighted.overall.risk%>%
                filter(!LoE.weight==0)%>%
                filter(!(LoE=='Spatial' & Consequence%in%1:3))%>%  #only assigning Risk to C = 4 for spatial
                group_by(Species,Consequence)%>%
                mutate(Weighted.Likelihood=weighted.mean(x=Likelihood,w=LoE.weight))%>%
                dplyr::select(Species,LoE,Consequence,Weighted.Likelihood,Probability,w)%>%
                rename(Likelihood=Weighted.Likelihood)%>%
                distinct(Species,Consequence,Likelihood,w)%>%
                mutate(Likelihood=round(Likelihood))

Store.risk_Drop.species=fn.risk.figure(d=LOE.risks$PSA, Risk.colors=RiskColors, out.plot=FALSE)%>%
  dplyr::select(-LoE)

Store.risk_Other.sp=fn.risk.figure(d=Weighted.overall.risk%>%filter(tolower(Species)%in%Lista.sp.outputs$Other.sp),
                                   Risk.colors=RiskColors,
                                   out.plot=TRUE)
ggsave(paste(Rar.path,"Risk_Other.sp.tiff",sep='/'),width = 10,height = 8,compression = "lzw")

Store.risk_Indicator.sp=fn.risk.figure(d=Weighted.overall.risk%>%filter(tolower(Species)%in%Lista.sp.outputs$Indicator.sp),
                                       Risk.colors=RiskColors,
                                       out.plot=TRUE)
ggsave(paste(Rar.path,"Risk_Indicator.sp.tiff",sep='/'),width = 10,height = 8,compression = "lzw")

  #2.4. Export Risk table
Out.overall.risk=rbind(Store.risk_Drop.species,Store.risk_Other.sp,Store.risk_Indicator.sp)%>%
  mutate(Risk.overall=paste0(Risk,' (',Consequence,'x',Likelihood,')'))%>%
  dplyr::select(Species,Risk.overall)
write.csv(Table.risks%>%left_join(Out.overall.risk,by='Species'),
          paste(Rar.path,'Table 13. Risk of each LoE and Overall.csv',sep='/'),row.names=F)

  #2.5. Export future Risk from JABBA and Integrated
write.csv(rbind(Risk.JABBA_future,Risk.integrated_future),
          paste(Rar.path,'Table 14. Risk_Future.csv',sep='/'),row.names=F)

  #2.6. Export current and future Risk based on F from JABBA and Integrated  
write.csv(rbind(Risk.JABBA_f,Risk.integrated_f),
          paste(Rar.path,'Table 15. Risk_based on F_Current and Future.csv',sep='/'),row.names=F)


#3. Display final risk for all species combined 
Final.risk_Drop.species=Store.risk_Drop.species%>%
                          group_by(Score,Risk)%>%
                          tally()%>%
                          mutate(Species=paste0("PSA-only species \n", "(n=",n,")"))%>%
                          dplyr::select(Species,Score,Risk)
Final.risk_Other.sp=Store.risk_Other.sp%>%dplyr::select(Species,Score,Risk)
Final.risk_Indicator.sp=Store.risk_Indicator.sp%>%dplyr::select(Species,Score,Risk)

  #3.1 Overall risk by species (combined PSA species)
fn.risk.all.sp(d=rbind(Final.risk_Drop.species,Final.risk_Other.sp,Final.risk_Indicator.sp))
ggsave(paste(Rar.path,"Risk_all species together.tiff",sep='/'),width = 8,height = 10,compression = "lzw")

  #3.2 Overall risk by species (all species)
add.non.interacting.species=FALSE
p1=rbind(Store.risk_Drop.species%>%
          dplyr::select(Species,Score,Risk),
        Final.risk_Other.sp ,
        Final.risk_Indicator.sp)
if(add.non.interacting.species)  #display species not interacting with fishing?
{
  All.WA.species= #missing
  non.interacting.species=data.frame(Species=All.WA.species,
                                     Score=2,
                                     Risk='Negligible')%>%
    filter(!Species%in%c(p$Species,capitalize(assessed.elsewhere)))
  p1=rbind(p1,non.interacting.species)
}
fn.risk.all.sp.eye(d=p1,show.all.risk.cat=TRUE)  
ggsave(paste(Rar.path,"Risk_all species together_proportion.tiff",sep='/'),width = 8,height = 8,compression = "lzw")

  #3.3 Each LoE risk and overall risk  
fn.risk.figure.all.LOE(d=Store.risks,
                       d1=Out.overall.risk,
                       lbl.cols=label_colors,
                       RiskCls=RiskColors)
ggsave(paste(Rar.path,"Risk_all LoE for each species.tiff",sep='/'),width = 10,height = 8,compression = "lzw")

fn.risk.figure.all.LOE(d=Store.risks%>%filter(!Species%in%capitalize(names(Indicator.species))),
                       d1=Out.overall.risk%>%filter(!Species%in%capitalize(names(Indicator.species))),
                       lbl.cols=label_colors,
                       RiskCls=RiskColors)
ggsave(paste(Rar.path,"Risk_all LoE for each species_non_indicators only.tiff",sep='/'),width = 10,height = 8,compression = "lzw")


#Compare MSY estimates by method  
Catch.only_MSY_Indicator.sp=read.csv(paste0(Rar.path,'/Table 4. Catch.only_catch_vs_MSY_Indicator.sp.csv'))
Catch.only_MSY_Indicator.sp.appendix=read.csv(paste0(Rar.path,'/Table 4. Catch.only_catch_vs_MSY_Indicator.sp_Appendix.csv'))
Catch.only_MSY_Other.sp=read.csv(paste0(Rar.path,'/Table 4. Catch.only_catch_vs_MSY_Other.sp.csv'))
Catch.only_MSY_Other.sp.appendix=read.csv(paste0(Rar.path,'/Table 4. Catch.only_catch_vs_MSY_Other.sp_Appendix.csv'))
JABBA_MSY_Indicator.sp=read.csv(paste0(Rar.path,'/Table 8. JABBA CPUE_estimates_Indicator.sp.csv'))
JABBA_MSY_Other.sp=read.csv(paste0(Rar.path,'/Table 8. JABBA CPUE_estimates_Other.sp.csv'))
SS_MSY_Indicator.sp=read.csv(paste0(Rar.path,'/Table 11. Age.based_SS_quantities_Indicator.sp.csv'))
SS_MSY_Other.sp=read.csv(paste0(Rar.path,'/Table 11. Age.based_SS_quantities_Other.sp.csv'))
Catch.only_MSY_Indicator.sp=Catch.only_MSY_Indicator.sp%>%
                              dplyr::select(Species,MSY_Lower.95,MSY_Median,MSY_Upper.95)%>%
                              mutate(Method='CoM_ensemble')
Catch.only_MSY_Other.sp=Catch.only_MSY_Other.sp%>%
                            dplyr::select(Species,MSY_Lower.95,MSY_Median,MSY_Upper.95)%>%
                            mutate(Method='CoM_ensemble')
Catch.only_MSY_Indicator.sp.appendix=Catch.only_MSY_Indicator.sp.appendix%>%
                        mutate(Method=paste0('CoM_',Model))%>%
                        dplyr::select(Species,MSY_Lower.95,MSY_Median,MSY_Upper.95,Method)
Catch.only_MSY_Other.sp.appendix=Catch.only_MSY_Other.sp.appendix%>%
                        mutate(Method=paste0('CoM_',Model))%>%
                        dplyr::select(Species,MSY_Lower.95,MSY_Median,MSY_Upper.95,Method)
JABBA_MSY_Indicator.sp=JABBA_MSY_Indicator.sp%>%
                        filter(Parameter=="MSY")%>%
                        mutate(Method=paste0('JABBA_',Scenario))%>%
                        rename(MSY_Lower.95=Lower.95,
                               MSY_Median=Median,
                               MSY_Upper.95=Upper.95)%>%
                    dplyr::select(Species,MSY_Lower.95,MSY_Median,MSY_Upper.95,Method)
JABBA_MSY_Other.sp=JABBA_MSY_Other.sp%>%
                  filter(Parameter=="MSY")%>%
                  mutate(Method=paste0('JABBA_',Scenario))%>%
                  rename(MSY_Lower.95=Lower.95,
                         MSY_Median=Median,
                         MSY_Upper.95=Upper.95)%>%
                  dplyr::select(Species,MSY_Lower.95,MSY_Median,MSY_Upper.95,Method)
SS_MSY_Indicator.sp=SS_MSY_Indicator.sp%>%
                    filter(Label=='MSY')%>%
                    mutate(Method=paste0('SS3_',Scenario),
                           MSY_Lower.95=Median-1.96*SE,
                           MSY_Median=Median,
                           MSY_Upper.95=Median+1.96*SE)%>%
              dplyr::select(Species,MSY_Lower.95,MSY_Median,MSY_Upper.95,Method)
SS_MSY_Other.sp=SS_MSY_Other.sp%>%
                filter(Label=='MSY')%>%
                mutate(Method=paste0('SS3_',Scenario),
                       MSY_Lower.95=Median-1.96*SE,
                       MSY_Median=Median,
                       MSY_Upper.95=Median+1.96*SE)%>%
                dplyr::select(Species,MSY_Lower.95,MSY_Median,MSY_Upper.95,Method)

MSY_combined=rbind(Catch.only_MSY_Indicator.sp,Catch.only_MSY_Other.sp,Catch.only_MSY_Indicator.sp.appendix,
                   Catch.only_MSY_Other.sp.appendix,JABBA_MSY_Indicator.sp,JABBA_MSY_Other.sp,
                   SS_MSY_Indicator.sp,SS_MSY_Other.sp)%>%
            arrange(Species,Method)
write.csv(MSY_combined,paste0(Rar.path,'/Table_Compare MSY estimates.csv'),row.names = F)

fn.compare.MSY(d=MSY_combined%>%filter(Species%in%capitalize(names(Indicator.species))),ncols=2)
ggsave(paste0(Rar.path,'/Compare MSY estimates_indicators.tiff'), width = 6,height = 6, dpi = 300, compression = "lzw")

fn.compare.MSY(d=MSY_combined%>%filter(!Species%in%capitalize(names(Indicator.species))),
               ncols=5,xlab.angle=90,xlab.size=8,Str.siz=9)
ggsave(paste0(Rar.path,'/Compare MSY estimates_other species.tiff'), width = 9,height = 10, dpi = 300, compression = "lzw")



#---28. Outputs for strategic papers  ------------------------------------------------- 
# Sawfish stock assessment  paper
fn.source1("Outputs_for_strategic_papers.r")
if(do.sawfish) fn.do.sawfish()
