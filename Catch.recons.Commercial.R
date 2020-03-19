#-- Script for reconstructing time series of commercial catch of sharks in WA and IUU

#notes: 'wetline' (as reported in Heupel & McAuley 2007) groups the byproduct of several 
#       different fisheries. This is already included in CAESS
#       Prior November 2006: 
#         all fisheries should be in Data.monthly, Data.monthly.north 
#         hence, all I need to do is to reapportion 'shark,other' using the
#         observed species composition, and if not available, the composition of
#         the reported catch that is reported in CAESS (as done in 'Assessment_other.species.R')

#       Post November 2006: 
#         All sharks from non-shark fisheries are not allowed to be retained (a few excempted fisheries),
#         hence this discarded catch is estimated for all fisheries that reported shark catches
#         prior to November 2006 as total landings X proportion shark X PCM  X BRD (trawl fisheries only)
#         If catch composition is not in weight, then convert numbers to weights

#ANNUAL UPDATES:

#1. "#Total landings time series". For fisheries listed in 'Calculate.discarding_catch' each year 
         #  download annual reported landings from FISHCUBE using the FishCube fishery code 
         #  from 'Lista.reap.FishCubeCode' or could use the link provided in 
         # Calculate.discarding.xlsx and run an update)

#2. Annual effort for "Kimberley.GBF.annual.effort"



### MISSING ### 
#   Commonwealth (waiting from AFMA's Julie Cotsell/ Ryan Murpthy): 1.2.2 Commonwealth GAB trawl and Western Tuna and Billfish Fisheries (WTBF)
#   Whaler_SA (waiting to hear from SARDI's Angelo Tsolos, Paul Rogers)
#   add MOU box catches, so far only reconstructed ILLEGAL!!
####

library(tidyverse)
library(readxl)
library(lubridate)
library(tm)

# 1 -------------------PARAMETERS SECTION------------------------------------

Asses.year=2019    #enter year of assessment
Last.yr.ktch="2017-18"  #enter year of last complete catches

Shark.protection.yr=2007   #Commercial protection in non-shark fisheries came in November 2006 (Heupel & McAuley 2007 page 74)

use.effort=FALSE  #whether to use Effort or catch to derive discards

Do.recons.paper="NO"  #outputs for reconstruction paper

Historic.yrs=1941:1975

#Hammerhead species composition North and South of 26 S
Comp.hh.south=data.frame(SPECIES=c(19004,19001,19002),   #McAuley & Simpfendorfer 2003 ratios for TDGDLF
                         Name=c('Smooth','Scalloped','Great'),
                         Prop=c(.97,.03/2,.03/2))

Comp.hh.north=data.frame(SPECIES=c(19004,19001,19002),   # from Sharks database (Naturaliste trip)
                         Name=c('Smooth','Scalloped','Great'),
                         Prop=c(.01,.67,.32))
  #convert observed ratio of numbers to ratio of weights using species maximum weigth 
Max.weight.great=(1.23e-3)*(445^3.24)/1000       #Stevens & Lyle 1989 MFR
Max.weight.scalloped=(3.99e-3)*(346^3.03)/1000   #Stevens & Lyle 1989 MFR
Max.weight.smooth=400     #Wikipedia (https://en.wikipedia.org/wiki/Smooth_hammerhead#cite_note-bester-8)
Max.weight.tot=Max.weight.great+Max.weight.scalloped+Max.weight.smooth
Comp.hh.south=Comp.hh.south%>%
              mutate(Prop=ifelse(Name=='Smooth',Prop*(Max.weight.smooth/Max.weight.tot),
                          ifelse(Name=='Scalloped',Prop*(Max.weight.scalloped/Max.weight.tot),
                          ifelse(Name=='Great',Prop*(Max.weight.great/Max.weight.tot),
                                 NA))),
                     Prop=Prop/sum(Prop))
Comp.hh.north=Comp.hh.north%>%
              mutate(Prop=ifelse(Name=='Smooth',Prop*(Max.weight.smooth/Max.weight.tot),
                          ifelse(Name=='Scalloped',Prop*(Max.weight.scalloped/Max.weight.tot),
                          ifelse(Name=='Great',Prop*(Max.weight.great/Max.weight.tot),
                                 NA))),
                     Prop=Prop/sum(Prop))

  
#Post Capture Mortality
  #no data for these species
assumed.sawfish.trawl=.3   
assumed.sawfish.gn=.1
assumed.sawshark.trawl=.3
assumed.wobbie.trawl=.3
assumed.greynurse.ll=.1
assumed.guitarfish.gn=.3

  #Assumed 30% increase from AVM to PCM
Inc.AVM=1.3  

  #Ellis et al 2016. Values are for AVM so it's bumped up by Inc.AVM
PCM=data.frame(Group=c("Sawfish","Sawsharks","Wobbegongs","Mackerel","Greynurse","Hexanchids","Whalers",
                       "Hammerheads","Triakids","Angels","Dogfish","Catsharks",
                       "Guitarfish","Numbfish","Rajids","Eagle.rays","Other.rays",
                       "Dasyatids"),
               Trawl=c(assumed.sawfish.trawl,assumed.sawshark.trawl,assumed.wobbie.trawl,rep(NA,3),
                       .519,.98,.3,.6,.21,.14,.185,.4,.24,.25,.42,.56),
               GN=c(assumed.sawfish.gn,.31,0,.66,.41,.59,.47,.94,.68,.34,.133,.28,assumed.guitarfish.gn,NA,
                    .06,.13,.5,.07),
               LL=c(NA,NA,0,.28,assumed.greynurse.ll,NA,.21,1,0,rep(NA,8),.12))%>%
      mutate(Trawl=sapply(Trawl,function(x)min(x*Inc.AVM,1)),
             GN=sapply(GN,function(x)min(x*Inc.AVM,1)),
             LL=sapply(LL,function(x)min(x*Inc.AVM,1)))



# 2 -------------------DATA SECTION------------------------------------

options(stringsAsFactors = FALSE)
fn.hndl=function(x)paste('C:/Matias/Data/Catch and Effort/',x,sep='')

## Fishery codes
FisheryCodes=read_excel(fn.hndl('FisheryCodeTable.xlsx'), sheet = "CAEStoFISHCUBE")


## Total landings time series (to be updated each year)                          
  #WA fisheries
Calculate.discarding_catch=read_excel(fn.hndl('Calculate.discarding.xlsx'), sheet = "Sheet1")%>%  #catch in kg
  data.frame%>%
  rename(KG=Weight..kg..Total,
         FINYEAR=Financial.Year,
         Group=Species.Group)

  #Non-WA fisheries
Whaler_SA=read.csv(fn.hndl("SA_marine_scalefish_whaler_ktch.csv"))  #catch in tonnes; source SARDI's Angelo Tsolos, Paul Rogers
WTBF_catch=read.csv(fn.hndl("WTBF_catch_Benseley.et.al.2010.csv"))  #catch in kg; source AFMA's Julie Cotsell/ Ryan Murpthy
GAB.trawl_catch=read.csv(fn.hndl("Gab.trawl_catch_Benseley.et.al.2010.csv"))  #catch in kg; source AFMA's Julie Cotsell/ Ryan Murpthy


## Species codes
All.species.names=read.csv("C:/Matias/Data/Species_names_shark.only.csv") #for catch


## 2.1 Catch_WA Fisheries
setwd("C:/Matias/Analyses/Data_outs")

  #-- 2.1.1 all monthly return and daily logbook fisheries from 1975-76

    #south of 26 S
Data.monthly=read.csv("Data.monthly.csv")
#note:  this has catch of sharks from monthly returns of all fisheries 
#       but from daily records of shark fisheries only 

daily.other=read.csv("Data.daily.other.fisheries.csv")                
#note:  this has shark catch from daily records from other fisheries   
  
    #north of 26 S
Data.monthly.north=read.csv("Data.monthly.NSF.csv")
Data.monthly.north$LIVEWT.c=Data.monthly.north$LIVEWT  


  #-- 2.1.2 Historic shark fisheries (pre 1975)

      #2.1.2.1 Whitley 1944 (citation in Simpfendorfer & Donohue 1998)
#note: All shark catch across southern WA (south of Abrolhos islands). Weight in pounds
Catch_1941_1943=data.frame(year=1941:1943,LIVEWT.c=c(55332,77441,109064)) 
Catch_1941_1943$LIVEWT.c=Catch_1941_1943$LIVEWT.c*0.453592/1000  #Convert to tonnes

      #2.1.2.2 Heald 1987
#note: All shark catch by port. Weight in kgs
Catch_1949=read.csv(fn.hndl("Heald(1987)_1949_total_catch_by port_allspecies.csv")) 
Catch_1949$LIVEWT.c=Catch_1949$Shark.catch.kg_live_wgt./1000   #in tonnes

      #2.1.2.3 Simpfendorfer & Donohue 1998
#note: All shark catch for southern WA (Geraldton to SA border). Weight in tonnes
Catch_1950=data.frame(year=1950,LIVEWT.c=50)   

      #2.1.2.4 Heald 1987
#note= All shark catch for southern WA (Geraldton to SA border). Weight in tonnes
Catch_1952_1975=read.csv(fn.hndl("Historical_WA_shark_catches_1952_1975.csv")) 
names(Catch_1952_1975)=c("year","LIVEWT.c")


  #-- 2.1.3 Wetline in the Western Rock lobster fishery
WRL=read.csv(fn.hndl("WRL/Number-of-vessels.csv"))
WRL.Wann=read.csv(fn.hndl("WRL/Wann_catch.csv"))

WRL.prop=0.1  #small number of operators used the gear (Taylor et al 2015)   
WRL.assumed.weight=100  #assumed weight of sharks in kg
WRL.copper.dusky.prop=0.1  #average proportion copper sharks

  #-- 2.1.4 TEPS   

    #2.1.4.1 TDGDLF
#note: data obtained from Comments in TDGDLF returns
TEPS=read.csv(paste('C:/Matias/Data/Catch and Effort',
              paste(Asses.year-2,substr(Asses.year-1,3,4),sep="_"),
              "TEPS_PROTECTEDSP.csv",sep='/'))
  #Length weigths
bwt=3.47e-06   #dusky
awt=3.10038
bwt.grey=5.4511   #greynurse Otway et al "Documentation of depth-related..."
awt.grey=3.1716
Weight=data.frame(Species=c('BT','BW','CP','GB','HZ',
                            'LG','MS','TG','WP'),    #FishBASE
                  CAES_Code=c(18014,18003,18001,18004,19004,
                              18023,10001,18022,10003),
                  awt=c(2.90,3.10038,2.9,3.24,2.86,
                        3.07,3.12,3.25,3.04),
                  bwt=c(0.0107,3.47e-06,0.01040,0.0046,0.0083,
                        0.0037,0.0054,0.0028,0.01))


  #-- 2.1.5 Pilbara trawl shark species composition (go to scientist: Corey W.) 
#source: Table 6.4 Heupel & McAuley 2007 (weight is live weight in kg)

    #shark catch composition
Pilbara.trawl.observed.comp=data.frame(
  Common.name=c('Sandbar shark','Great hammerhead','Pigeye shark',
                'Fossil shark','Blacktip sharks','Tiger shark',
                'Milk shark','Weasel shark','Sharpnose shark','Whitecheek shark',
                'Leopard shark','Tawny nurse shark','Scalloped hammerhead','Spot-tail shark',
                'Smooth hammerhead','Longtail carpet sharks','Banded catshark','Winghead',
                'Tasseled woggegong','Bignose shark','Spinner shark',
                'Sliteye','Northern wobbegong'),
  Scientific.name=c('Carcharhinus plumbeus','Sphyrna mokarran','Carcharhinus amboinensis',
                    'Hemipristis elongata','Carcharhinus tilstoni & C. limbatus','Galeocerdo cuvier',
                    'Rhizoprionodon acutus','Hemigaleus microstoma','Rhizoprionodon taylori','Carcharhinus dussumieri',
                    'Stegastoma fasciatum','Nebrius ferrugineus','Sphyna lewini','Carcharhinus sorrah',
                    'Sphyrna zygaena','Hemiscylliidae','Chiloscylium punctatum','Eusphyra blochii',
                    'Eucrossorhinus dasypogon','Carcharhinus altimus','Carcharhinus brevipinna',
                    'Loxodon macrorhinus','Orectolobus wardi'),
  SPECIES=c(18007,19002,18026,NA,18014,18022,18006,NA,NA,NA,NA,13010,19001,18013,19004,NA,NA,NA,
            NA,18012,18023,NA,NA),
  Retained=c(rep('yes',10),rep('no',13)),
  Numbers=c(90,7,3,27,80,5,238,559,114,102,101,4,115,25,1,12,43,1,4,6,1,1,2),
  Weight=c(2058.6,595.9,295.1,395.8,593.2,168.3,522.6,772.7,210.5,214.2,
           1627.8,257.1,189.5,90.7,39.4,20.0,15.6,15.1,9.5,4.9,4.2,2,.8)
  )%>%
    mutate(Prop=Weight/sum(Weight))

    #ratio observed shark catch to landed catch
Pilbara.trawl.observed.landing=12901     # in kg; Table 9 Stephenson & Chidlow 2003
Pilbara.trawl.observed.shrk= 562 + 885   # in kg; observed retained and discarded species
Pilbara.trawl_shrk.to.land=Pilbara.trawl.observed.shrk/Pilbara.trawl.observed.landing

    #BRD
BRD_pilbara.trawl_prop.shark=0.61  #Wakefield et al 2017                     
BRD_pilbara.trawl_prop.ray=0.36  
BRD_pilbara.trawl_year='2003-04'

    #Observed effort
  #Pilbara.trawl.observed.effort=100  #days (Stephenson & Chidlow 2003 fide in McAuley et al 2005 page 25)
  Pilbara.trawl.observed.effort=1601.898   #hours; derived from get.observed.Pilbara
  get.observed.Pilbara=FALSE
  if(get.observed.Pilbara)
  {
    ## Shark bio data
    User="Matias"
    source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R")
    Observed.hours.Pilbara=DATA.bio%>%
      filter(Method=='TW')%>%
      distinct(SHEET_NO)%>%
      pull(SHEET_NO)
    Observed.hours.Pilbara=DATA%>%
      filter(SHEET_NO%in%Observed.hours.Pilbara)%>%
      distinct(SHEET_NO,.keep_all = T)
    #Stephenson & Chidlow repoted 100 days observed but only 83 records in DATA
    # hence, add 17 days with mean observed hours per day
    Mean.hour=Observed.hours.Pilbara%>%
      group_by(date)%>%
      summarise(Total.hours=sum(SOAK.TIME,na.rm=T))
    Mean.hour=mean(Mean.hour$Total.hours)
    Total.observed.hours=sum(Observed.hours.Pilbara$SOAK.TIME,na.rm=T)+17*Mean.hour
  }


  #-- 2.1.6 WA scallop and prawn trawl fisheries (go to scientist: Mervi Kangas)

    #South West Trawl (Laurenson et al 1993 Table 5)
# notes: Authors grouped all sites and replicates in Table 5 
#        Authors reported a range of number of individuals so the mid point is used

    #shark catch composition
South.west.trawl.observed.comp=data.frame(
    Scientific.name=c('Squatina australis','Aulohalaelurus labiosus','Urolophus circularis','Orectolobus tentaculatus',
             'Pristiophorus cirratus','Myliobatis australis','Trygonorhina fasciata','Mustelus antarcticus',
             'Urolophus lobatus','Trygonoptera personata','Hypnos monopterygium','Heterodontus portusjacksoni',
             'Sphyrna zygaena','Dasyatis brevicaudata','Aptychotrema vincentiana','Urolophus paucimaculatus',
             'Parascyllium variolatum','Urolophus mucosus','Orectolobus sp.'),
  
    Common.name=c('Angel shark','Black-spotted catshark','Circular stingaree','Cobbler carpetshark',
            'Common sawshark','Eagle ray','Fiddler ray','Gummy shark',
            'Lobed stingaree','Masked stingaree','Numbfish','Port Jackson shark',
            'Smooth hammerhead','Smooth stingray','Southern shovelnose ray','Sparsely-spotted stingaree',
            'Varied catshark','Western stingaree','Western wobbegong'),
    SPECIES=c(24900,NA,NA,NA,23002,31000,26999,17001,NA,NA,NA,NA,19004,NA,NA,NA,
              NA,NA,13000),                              
    Commercial=c('y',rep('n',3),rep('y',4),rep('n',4),'y','n','y',rep('n',3),'y'),
    Bell.buoy=c(5,0,0,0,0,5,5,0,75,30,0,5,0,0,5,30,0,30,0),
    Cottesloe=c(5,0,0,0,0,30,5,0,75,30,0,550,0,5,5,5,0,30,0),
    Conventry.reef=c(5,0,0,0,0,5,5,5,550,30,0,5,0,5,5,75,0,5,0),
    Comet.bay=c(30,0,0,0,5,5,5,5,5,5,5,30,5,0,5,5,0,30,0),
    Preston.deep=c(5,5,5,0,5,5,5,0,75,30,5,5,0,0,5,5,0,30,0),
    Preston.shallow=c(0,5,0,0,0,5,5,0,0,5,5,5,0,0,0,0,5,5,5),
    Capel=c(5,0,0,0,0,5,5,0,30,30,0,5,0,0,5,550,0,30,0),
    Busselton=c(5,0,0,5,5,5,5,0,0,5,0,5,0,0,0,75,0,5,0),
    zoneC=c(30,0,0,0,5,5,5,0,0,5,0,30,0,0,5,75,0,75,0)
    )%>%
     mutate(N=Bell.buoy+Cottesloe+Conventry.reef+Comet.bay+Preston.deep+Preston.shallow+Capel,
            Prop=N/sum(N))

South.west.trawl.observed.PCM= 0.5 #observed for gummy sharks after 7 days of capture

    #ratio observed shark catch to landed catch
South.west.trawl.observed.landing=read.csv(fn.hndl("Observed.commercial_Laurenson et al 1993.csv"))
South.west.trawl.observed.landing=South.west.trawl.observed.landing%>%
                      mutate(N=Bell.Buoy+Cottesloe+Coventry.Reef+Comet.Bay+
                               Preston.Deep+Preston.Shallow+Capel+Busselton+Zone.C)
South.west.trawl_shrk.to.land=sum(South.west.trawl.observed.comp$N)/sum(South.west.trawl.observed.landing$N)

    #observed effort
  South.west.trawl.observed.effort=9*4*(15/60)*4  #hours, 9 sites X 4 replicates X 15 mins X 4 times in a year (Laurenson et al 1993 page 13)


    #Shark Bay, Exmouth and Onslow prawn trawl (Kangas et al 2006 Appendix 3.1 and 3.2)
# notes: Authors grouped all sites and replicates in Appendix 3.1 and 3.2 and only 
  #      report presence/absence 

    #shark catch composition
Shark.Bay.trawl.observed.comp=data.frame(
  Scientific.name=c('Chiloscyllium punctatum','Halaelurus boesemani','Mustelus sp. A',
             'Carcharhinus cautus','Carcharhinus melanopterus','Carcharhinus obscurus',
             'Carcharhinus plumbeus','Rhizoprionodon acutus Shark','Hemigaleus australiensis',
             'Rhina ancylostoma','Rhynchobatus australiae',
             'Aptychotrema vincentiana','Narcine westraliensis','Hypnos monopterygium',
             'Dasyatis kuhlii','Dasyatis leylandi','Himantura sp.','Himantura toshi',
             'Himantura uarnak','Taeniura meyeni','Gymnura australis','Trygonoptera ovalis',
             'Aetobatus narinari'),
  Common.name=c('Catshark, Brown-banded','Catshark, Speckled','Shark, Grey Gummy',
            'Shark, Nervous','Shark, Blacktip Reef','Shark, Dusky Whaler',
            'Shark, Sandbar','Shark, Milk','Shark, Weasel','Shark Ray','Shovelnose Ray, White-spotted',
            'Shovelnose Ray, Western','Numbfish, Banded','Numbfish',
            'Stingray, Blue-spotted','Stingray, Brown reticulated','Stingray, Coachwhip','Whipray, Black-spotted',
            'Whipray, Reticulate','Stingray, Black-blotched','Ray, Rat-tailed/Butterfly','Stingaree, Striped',
            'Ray, White-spotted Eagle'),
  SPECIES=c(NA,NA,NA,18034,18036,18003,18007,18006,NA,NA,26999,26999,NA,NA,rep(31000,7),NA,31000),                              
  Commercial=c('n','n',rep('y',7),'n','n','y',rep('n',10),'y'),
  Total=c(2,1,3,1,1,1,1,1,1,1,3,3,4,1,3,4,1,2,2,1,4,2,1))%>%
  mutate(Prop=Total/sum(Total))
Shark.Bay.trawl.observed.effort=26*3*(10/60)*4  #hours, 26 sites X 3 replicates X 10 mins X 4 times in a year (Kangas et al 2006 Table 2.1)
Shark.Bay.trawl.observed.number.sharks=185  #Mervi Kangas pers comm
Shark.Bay.trawl.observed.number.commercial=36735     #Mervi Kangas pers comm 
Trawl.shark.mean.weight=2                     #2 kg assumed mean weight because only small individuals are not excluded by BRDs
Trawl.prawn.scallop.mean.weight=mean(c(.025,.02))           #Mervi Kangas pers comm


Exmouth.Onslow.trawl.observed.comp=data.frame(
  Scientific.name=c('Chiloscyllium punctatum','Eucrossorhinus dasypogon','Atelomycterus sp.','Rhizoprionodon acutus',
             'Hemigaleus australiensis','Rhynchobatus australiae','Rhinobatos typus','Dasyatis kuhlii',
             'Dasyatis leylandi','Himantura toshi','Gymnura australis',
             'Aetomylaeus vespertilio','Aetomylaeus nichofii'),
  Common.name=c('Catshark, Brown-banded','Wobbegong, Tasselled','Catshark, Banded','Shark, Milk',
            'Shark, Sicklefin Weasel','Ray, White-spotted Shovelnose','Ray, Giant Shovelnose','Ray, Blue-spotted Stingray',
            'Ray, Brown Reticulated','Ray, Black-spotted Whipray','Ray, Butterfly/Rat-tailed',
            'Ray, Ornate Eagle Ray','Ray, Banded Eagle Ray'),
  SPECIES=c(NA,NA,NA,18006,NA,26999,26999,31000,31000,31000,31000,31000,31000),                              
  Commercial=c(rep('n',3),'y','y',rep('n',8)),
  Total=c(3,1,3,2,2,2,1,1,3,3,3,1,1))%>%
  mutate(Prop=Total/sum(Total))
Exmouth.Onslow.trawl.observed.effort=26*3*(10/60)*3  #hours, 26 sites X 3 replicates X 10 mins X 3 times in a year (Kangas et al 2006 Table 2.1)
Exmouth.Onslow.trawl.observed.number.sharks=92        #Mervi Kangas pers comm
Exmouth.Onslow.trawl.observed.number.commercial=15146   #Mervi Kangas pers comm  

    #ratio observed shark catch to landed catch             
Shark.Bay.trawl_shrk.to.land=(Shark.Bay.trawl.observed.number.sharks*Trawl.shark.mean.weight)/(Shark.Bay.trawl.observed.number.commercial*Trawl.prawn.scallop.mean.weight)
Exmouth.Onslow.trawl_shrk.to.land=(Exmouth.Onslow.trawl.observed.number.sharks*Trawl.shark.mean.weight)/(Exmouth.Onslow.trawl.observed.number.commercial*Trawl.prawn.scallop.mean.weight)

    #BRD
BRD_prawn.trawl_prop.shark=9/70  #Kangas & Thomson 2004 :70 sharks retained with no BRD VS 9 with BRD
BRD_prawn.trawl_prop.ray=8/65  #                         65 rays retained with no BRD VS 8 with BRD
BRD_prawn.trawl_year='2003-04'


  #-- 2.1.7 Kimberley Gillnet and Barramundi Fishery and Eighty Mile Beach Gillnet Fishery

  #source: Table 6.5 Heupel & McAuley 2007 (weights are live weight in kg)

    #shark catch composition
Kimberley.GBF.observed.comp=data.frame(
  Common.name=c('Blacktip sharks','Bull shark','Graceful shark','Hardnose shark',
                'Scalloped hammerhead','Winghead','Lemon shark','Spinner shark',
                'Milk shark','Nervous shark','Narrow sawfish','Dwarf sawfish',
                'Pigeye shark','Freshwater sawfish','Green sawfish',
                'Blacktip reef','Sliteye shark','Shovelnose','Spot-tail shark',
                'Stingray','Australian sharpnose','Tiger shark','Whitespot guitarfish'),
  Scientific.name=c('Carcharhinus tilstoni & C. limbatus','Carcharhinus leucas',
                    'Carcharhinus amblyrhynchoides','Carcharhinus macloti',
                    'Sphyrna lewini','Eusphyra blochii','Negaprion acutidens',
                    'Carcharhinus brevipinna','Rhizoprionodon acutus',
                    'Carcharhinus cautus','Anoxypristis cuspidata','Pristis clavata',
                    'Carcharhinus amboinensis','Pristis microdon','Pristis zijsron',
                    'Carcharhinus melanopterus','Loxodon macrorhinus',
                    'Rhynchobatidae/F Rhinobatidae','Carcharhinus sorrah',
                    'Dasyatididae','Rhizoprionodon taylori','Galeocerdo cuvier',
                    'Rhynchobatus australiae'),
  SPECIES=c(18014,18021,18033,NA,19001,NA,18029,18023,18006,18034,25002,25004,
            18026,25000,25001,18036,NA,26999,18013,31000,NA,18022,31000),                              
  Numbers=c(87,32,36,3,1,45,27,3,13,82,190,31,398,1,17,1,6,22,1,5,4,1,5),
  Weight=c(578.1,126.5,178.6,12.5,109.1,106.9,165.1,3.4,13.9,422.5,5239.1,370.1,
           2284.5,23.7,161.6,1.4,5.4,63.5,5.7,20.0,2.9,57.2,56.8)
)

    #observed effort
Kimberley.GBF.observed.effort=160  #days (McAuley et al 2005 page 26; 5 vessels observed)

    #total annual effort
fn.rid.efrt=function(d) paste('C:/Matias/Data/Catch and Effort/Effort_other_fisheries',d,sep='/')
Kimberley.GBF.annual.effort=read.csv(fn.rid.efrt('KGBF Annual catch and Bdays.csv'))    #days (KGBF and 80 mile beach combined)


#--Effort time series 
if(use.effort)
{
  #Pilbara trawl
  
  #Scallop and prawn trawl
  Prawn.scallop_Trawl.effort=read.csv(fn.rid.efrt('Prawn.scallop_Trawl fishing effort.csv')) #hours
  
  Exmouth.Onslow.trawl.annual.effort=Prawn.scallop_Trawl.effort[,c('Finyear','EGP','ONP')]%>%  # (note: combined Exmouth & Onslow)
    mutate(ONP=ifelse(is.na(ONP),0,ONP),
           Effort=EGP+ONP)%>%
    dplyr::select(Finyear,Effort)
  Shark.Bay.trawl.annual.effort=Prawn.scallop_Trawl.effort[,c('Finyear','SBP')]%>%
    rename(Effort=SBP)%>%
    mutate(Effort=ifelse(is.na(Effort),0,Effort))
  Shark.Bay.trawl.scallop.annual.effort=Prawn.scallop_Trawl.effort[,c('Finyear','SBS')] %>%
    rename(Effort=SBS)%>%
    mutate(Effort=ifelse(is.na(Effort),0,Effort))   
  Kimberley.trawl.annual.effort=Prawn.scallop_Trawl.effort[,c('Finyear','KP')] %>%
    rename(Effort=KP) %>%
    mutate(Effort=ifelse(is.na(Effort),0,Effort)) 
  Nickol.Bay.trawl.annual.effort=Prawn.scallop_Trawl.effort[,c('Finyear','NBP')] %>%
    rename(Effort=NBP) %>%
    mutate(Effort=ifelse(is.na(Effort),0,Effort))  
  Abrolhos.trawl.annual.effort=Prawn.scallop_Trawl.effort[,c('Finyear','AIS')]%>%
    rename(Effort=AIS)%>%
    mutate(Effort=ifelse(is.na(Effort),0,Effort))
  
  #South.west.trawl.annual.effort=Prawn.scallop_Trawl.effort[,]  
  
}


## 2.2. Catch of non WA Fisheries

  #-- 2.2.1 Taiwanese gillnet (1974-1986) and longline (1989-91) fishery (catch in tonnes)
#Sources: Stevens & Davenport 1991; Stevens 1999
Taiwan.gillnet.ktch=data.frame(           #in kg; From Figure 7 Stevens & Davenport
      Year=1974:1986,
      Total.Ktch=c(321.4,8997.6,6455.3,9970.7,5528.2,3282.1,5831.1,
                   6694.7,5624.1,7589.9,6544.2,2929.5,2111.1)*1000,
      Shark.Ktch=c(rep(NA,5),557.9,4234.5,4486.2,3639.4,4418.6,
                   2731.4,2327.4,2393.9)*1000,
      WA.Ktch=c(rep(NA,5),50,1250,720,800,790,400,10,70)*1000)  

Taiwan.gillnet.sp.comp=as.data.frame(matrix(c(.357,.166,.059),ncol=3))
colnames(Taiwan.gillnet.sp.comp)=c("Australian.blacktip.Shark.West.prop",
                                   "Spot.tail.Shark.West.prop","Hammerheads.West.prop")

Taiwan.longline.ktch=data.frame(Year=1990:1991,
                                Shark.Ktch=c(1700*11/(11+9),1700*9/(11+9))*1000)
Taiwan.longline.sp.comp=data.frame(
  Species=rep(c("Spot-tail shark","Australian blacktip shark","Tiger shark",
                "Milk shark","Spinner shark",
                "Pigeye shark","Graceful shark","Hammerheads","Other"),2),
  Year=c(rep(1990,9),rep(1991,9)), 
  Percent=c(18.3,25.6,8.1,.4,6.9,17.9,1.2,2,19.5,5.6,79.7,0.7,2,4.5,1,0,0.1,6.4))


  #-- 2.2.2 Commonwealth GAB trawl and Western Tuna and Billfish Fisheries (WTBF)   
if(!exists('WTBF_catch'))
{
  #Bensely et al 2010. Table 2 in Appendix A
  WTBF_effort=read.csv(fn.hndl("WTBF_effort_Benseley.et.al.2010.csv"))  #hook numbers
  GAB.trawl_effort=read.csv(fn.hndl("GAB.trawl_effort_Benseley.et.al.2010.csv"))  #hours trawled
  
  #Stobutski et al 2006 Table 3 
  WTBF_observed=read.csv(fn.hndl("WTBF_Stobutzki.et.al.2006.csv"))  #catch in numbers
  cpue_dusky_WTBF=37/203205  #number of individuals per 203205 hooks observed.
  cpue_sandbar_WTBF=8/203205
  dusky_WTBF.at.vessel.mortality=1-.97  # 97% discarded alive
  sandbar_WTBF.at.vessel.mortality=1-1  # 100% discarded alive
}


  #-- 2.2.3 SA Marine Scalefish fishery (source Taylor et al 2015)
#description: whaler shark catch from SA MArine Scale fishery (in tonnes). 
Whaler_SA_dusky.prop=.1  # Steer et al 2018 (page 148)
Whaler_SA_bronzie.prop=1-Whaler_SA_dusky.prop


  #-- 2.2.4 Indonesian illegal fishing in Australia waters

#source: ABC (https://www.abc.net.au/news/2019-11-12/illegal-shark-fishing-northern-territory-fishing-boat/11697036)
Indo_flesh=60 #flesh (kg)
Indo_fins=63
Indo_skins=16
Indo_prop.shark=8/36   #proportion of apprehended vessels fishing for sharks

#source: Edyvane & Penny 2017 (using all vessels, including MOU because the MOU catch is not accounted for anywhere)
Indo_sightings=data.frame(year=2000:2013,
                          Total.FFVs.inside.AEEZ=c(4867,5878,3047,9550,6638,9362,7378,4320,
                                                   6827,9117,9517,11822,13979,11455))
#source: Stacey 2007 Boats to Burn: Bajo Fishing Activity in the Australian Fishing Zone. 
Indo_apprehensions.Stacey=data.frame(year=1975:1999,
                               Apprehensions=c(3,rep(0,4),2,rep(0,4),5,0,1,46,29,43,38,15,
                                               23,111,76,97,122,rep(NA,2)))
#AFMA (sourced by Rik Buckworth)
Indo_apprehensions=data.frame(year=c(2000:2019),
                              Apprehensions=c(62,96,143,132,203,360,210,141,18,23,14,12,7,10,5,7,6,9,5,0))
Indo.prop.apprehen.with.shark=data.frame(Year=2008:2019,
                                         Total.appr=c(27,23,14,12,7,26,6,20,15,14,5,0),
                                         With.shark=c(10,14,10,7,5,3,0,6,5,1,0,0))%>%
                                mutate(prop=With.shark/Total.appr)

# ANAO 2010
Indo_jurisdiction.prop=data.frame(Jurisdiction=c('WA',"NT","QLD"),
                                  prop=c(.33,.34,.33))   #Figure 1.2 , equally split in 3 based on points scatter

#source: Marshall et al 2016 (MOU box)
Indo_shark.N.vessels=9  # 9 Type 2 vessels sampled over a period of 1 week (8th to 15th May 2015)
Indo_shark.comp=data.frame(Species=c('Silvertip shark','Bignose shark','Grey reef shark','Pigeye shark',
                                     'Spinner shark','Silky shark','Bull shark','Blacktips',
                                     'Dusky shark','Sandbar shark','Spot tail shark',
                                     'Tiger shark','Lemon shark',
                                     'Scalloped hammerhead','Great hammerhead',
                                     'Whitetip reef shark'),
                           Proportion=c(0.0022,0.0058,0.0141,0.0033,  #by weight
                                        0.0822,0.0011,0.0001,
                                        0.0259,0.0293,0.1524,
                                        0.001,0.665,0.0075,
                                        0.0029,0.0056,0.0017),
                           Proportion.by.number=c(0.0066,0.0132,0.0526,
                                                  0.0066,0.0724,0.0066,
                                                  0.0132,0.0395,0.0132,
                                                  0.4342,0.0132,0.2961,
                                                  0.0066,0.0066,0.0066,
                                                  0.0132))   
Indo_shark.weight=10486  #kg  Total weight of the shark catch for the monitored period
Indo_shark.TL.range=c(89,409)
Indo_shark.mature=.586    #Tables 5 and 6 have species-specific size and maturity ranges
#Indo_average.trip.length=23  #days for fishers fishing in MOU box


# Salini et al 2007
Indo_FFV_vessel.days.2005=sum(c(65,147,250,227,250,413,166,450,694,639,465,464)) #Table 7.2, zone 6
Indo_Prop.Shark.vesl.Salini=350/(350+75+10+30) #Figure 6-4 

#Assumptions
Indo_assumed.n.trips.per.year=10  #assumed number of trips per vessel per year 


#Vanesa Jaiteh's thesis
Indo_MOU.Vanesa=read.csv('C:/Matias/Data/Catch and Effort/Indonesia_Shark Data MasteR_updated_09_Sept_15.csv',stringsAsFactors = F)
#Indo_average.trip.length_MOU=mean(c(3*7,8*7))  #days 
#Indo_assumed.n.trips.per.year_MOU=round(28*365/639) #number of trips between March 2012 and November 2013
#Indo_average.shark.per.day_MOU=7 # 4 +/- 3 reported

#OECD
Indo_average.shark.per.trip_MOU=2600  #kg
Indo_trips.per.year_MOU=round(160/3)  #160 reported between 1992 and 1994
Indo_average.trip.length=4            #days (days in AFZ before being apprehended)
Indo_MOU.annual.trips.prop=data.frame(year=c(1975:2019),
                                      trips.prop=rep(1,length(1975:2019)))
Indo_missed.appr.rate=1.75  

# 3 -------------------PROCEDURE SECTION------------------------------------

## Fisheries code and Fisheries name
FisheryCodes=FisheryCodes%>%data.frame
Data.monthly=Data.monthly%>%left_join(FisheryCodes,by=c("FisheryCode"="SASCode"))%>%
                    mutate(DATA.type="Data.monthly")
Data.monthly.north=Data.monthly.north%>%left_join(FisheryCodes,by=c("FisheryCode"="SASCode"))%>%
                    mutate(FishCubeName=ifelse(FishCubeCode=='NCS',
                                               'Joint Authority Northern Shark Fishery',FishCubeName),
                           FishCubeCode=ifelse(FishCubeCode=='NCS','JANS',FishCubeCode))%>%
                    mutate(DATA.type="Data.monthly.north")
daily.other=daily.other%>%left_join(FisheryCodes,by=c("fishery"="SASCode"))%>%
              rename(SPECIES=species,
                     LIVEWT=livewt,
                     FisheryCode=fishery,
                     FINYEAR=finyear)%>%
              mutate(FishCubeCode=ifelse(is.na(FishCubeCode),FisheryCode,FishCubeCode),
                     DATA.type="daily.other")


## Add 'daily.other' to Data.monthly or Data.monthly.north accordingly
dummy.daily.oder=Data.monthly[1:nrow(daily.other),]
dummy.daily.oder[,]=NA
dummy.daily.oder=dummy.daily.oder%>%
                  mutate(FINYEAR=daily.other$FINYEAR,
                         MONTH=NA,
                         YEAR.c=daily.other$year,
                         SPECIES=daily.other$SPECIES,
                         LIVEWT.c=daily.other$LIVEWT,
                         FishCubeCode=daily.other$FishCubeCode,
                         FishCubeName=daily.other$FishCubeName,
                         DATA.type=daily.other$DATA.type)
Data.monthly=rbind(Data.monthly,subset(dummy.daily.oder,FishCubeCode%in%c("SWT","OASC")))
Data.monthly.north=rbind(Data.monthly.north,
                         dummy.daily.oder%>%
                            filter(FishCubeCode%in%c("NDS","PFT","OP","OANCGCWC"))%>%
                            dplyr::select(-NETLEN.c))

## Keep only sharks and rays and set NA FishCubeCode to 'unknwn'
Data.monthly=Data.monthly%>%
                  filter(SPECIES<50000)%>%
                  mutate(FishCubeCode=ifelse(is.na(FishCubeCode),"unknwn",FishCubeCode),
                         Discarded.ktch="NO")
Data.monthly.north=Data.monthly.north%>%
                  filter(SPECIES<50000)%>%
                  mutate(FishCubeCode=ifelse(is.na(FishCubeCode),"unknwn",FishCubeCode),
                         Discarded.ktch="NO")

if(Do.recons.paper=="YES")
{
  #-- Keep copy of original data with no reapportioning or discarding
  Data.monthly.original=Data.monthly
  Data.monthly.north.original=Data.monthly.north 
}


## Select fisheries for reapportioning 'shark other' & reconstructing discards post 2006
Prawn.Trawl.fisheries=c('EGP','SWT','SCT','SBP','NBP','KTR','KP','OP','SBSC','C156','AIMWT')
Scalefish.Trawl.fisheries=c('PFT')

if(Do.recons.paper=="YES")
{
  fn.hnd.out=function(x)paste('C:/Matias/Analyses/Reconstruction_catch_commercial/',x,sep='')
  source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")
  
  dit=rbind(daily.other%>%dplyr::select(SPECIES,FishCubeCode,FishCubeName,LIVEWT),
            Data.monthly%>%dplyr::select(SPECIES,FishCubeCode,FishCubeName,LIVEWT),
            Data.monthly.north%>%dplyr::select(SPECIES,FishCubeCode,FishCubeName,LIVEWT))%>%
            filter(SPECIES<50000 & !is.na(FishCubeName))%>%
            group_by(FishCubeCode,FishCubeName)%>%
            summarise(Total=round(sum(LIVEWT,na.rm=T)/1000,2))%>%    #in tonnes
            arrange(-Total)%>%
            data.frame%>%
            rename(Total.catch=Total,
                   Fishery=FishCubeName)%>%
    dplyr::select(-FishCubeCode)
  setwd('C:/Matias/Analyses/Reconstruction_catch_commercial')
  fn.word.table(WD=getwd(),TBL=dit,Doc.nm="Fisheries.reporting.shark_rays_tonnes",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
}
Calculate.discarding=c('PFT','C019','C066','C070','CSLP','EGBS','EGP','KP','KTR','NBP',
                       'OANCGCWC','OP','SBP','SBS','SBSC','SCT','SWT')
Reaportion.from.reported.ktch=c('C019','C066','C070','CSLP','EGBS','JANS','JASDGDL',
                                'KTR','NCS','OANCGCWC','OASC','SBS','WANCS','WCDGDL')
Reaportion.from.survey=c('PFT','EGP','KP','NBP','OP','SBP','SBSC','SCT','SWT')
Lista.reap.FishCubeCode=list(Pilbara.trawl='PFT',
                             Estuaries.19='C019',
                             Cockburn.Sound.fish.net='C066',
                             Power.hauler.70='C070',
                             Cockburn.Sound.line.pot='CSLP',
                             Exmouth.beach.seine.mesh='EGBS',
                             Exmouth.Gulf.prawn='EGP',
                             JANSF=c('JANS'),
                             JASDGDL='JASDGDL',
                             Kimberley.prawn='KP',
                             Kimberley.trap='KTR',
                             Nickol.Bay.prawn='NBP',
                             Open.north.gas.west='OANCGCWC',
                             Open.south='OASC',
                             Onslow.prawn='OP',
                             Shark.Bay.prawn='SBP',
                             Shark.Bay.snapper='SBS',
                             Shark.Bay.scallop='SBSC',
                             South.coast.trawl='SCT',
                             South.west.trawl='SWT',
                             WANCS='WANCS',
                             WCDGDL='WCDGDL')


## 3.1. Catch_WA Fisheries

  #3.1.1 Calculate commercial species discarding from non-shark fisheries since Shark.protection.yr
#note: this calculates the discarding of shark, then following step reapportions it to
#       species
All.data.yrs=sort(unique(Data.monthly$FINYEAR))
First.disc.yr=paste(Shark.protection.yr-1,substr(Shark.protection.yr,3,4),sep='-')
Discarding.yrs=All.data.yrs[match(First.disc.yr,All.data.yrs):length(All.data.yrs)]
Prop.reported.shark.ray=Calculate.discarding_catch%>%
  filter(Fishery.Code%in%Calculate.discarding & !FINYEAR%in%Discarding.yrs)%>%
  mutate(Group=ifelse(Group=="Sharks & Rays","Sharks_Rays",
               ifelse(Group=="Prawns & yabbies","Prawns_yabbies",Group)))%>%
  group_by(Fishery.Code,Group)%>%
  summarise(Tot=sum(KG))%>%
  spread(Group,Tot,fill=0)%>%
  data.frame%>%
  mutate(Prop.shark.ray=Sharks_Rays/(Cephalopods+Crustaceans+Echinoderms+
                                       Finfish+Molluscs+Prawns_yabbies))
for(i in 1:length(Calculate.discarding))
{
  Fishry=Calculate.discarding[[i]]
  
  #select years to calculate discarding
  Dis.yrs=Discarding.yrs
  if(Fishry%in%Prawn.Trawl.fisheries)
  {
    Dis.yrs=Calculate.discarding_catch%>%
      filter(Fishery.Code%in%Fishry)%>%
      distinct(FINYEAR)%>%pull(FINYEAR)
  }

  Landings=Calculate.discarding_catch%>%
              filter(FINYEAR%in%Dis.yrs & 
                     Fishery.Code==Fishry)%>%
              group_by(FINYEAR)%>%
              summarise(KG=sum(KG))%>%
              arrange(FINYEAR)%>%
              data.frame
  yrs.reported.ktch=unique(c(Data.monthly%>%filter(FishCubeCode%in%Fishry)%>%distinct(FINYEAR)%>%pull(FINYEAR),
                             Data.monthly.north%>%filter(FishCubeCode%in%Fishry)%>%distinct(FINYEAR)%>%pull(FINYEAR)))
  yrs.reported.ktch=yrs.reported.ktch[which(yrs.reported.ktch%in%Dis.yrs)]
  
  
  #Bycatch reduction devices
  BRD=1
  if(Fishry%in%Prawn.Trawl.fisheries) BRD=mean(c(BRD_prawn.trawl_prop.shark,BRD_prawn.trawl_prop.ray))
  if(Fishry%in%Scalefish.Trawl.fisheries) BRD=mean(c(BRD_pilbara.trawl_prop.shark,BRD_pilbara.trawl_prop.ray))
  
  #only calculate discarding for years since Shark.protection.yr 
  if(nrow(Landings)>0)  
  {
    #Calculate discards using landings ratio shark:other groups
    if(Fishry %in% Reaportion.from.reported.ktch)
    {
      Prop=Prop.reported.shark.ray%>%
                  filter(Fishery.Code==Fishry)%>%
                  pull(Prop.shark.ray)
      Shark.disc=Landings%>%mutate(Shark.ktch=KG*Prop*BRD)
      
      #remove catch if it's reported
      if(length(yrs.reported.ktch)>0)
      {
        Drop=Data.monthly.north%>%
          filter(FishCubeCode%in%Fishry & FINYEAR%in% yrs.reported.ktch)%>%
          group_by(FINYEAR)%>%
          summarise(LIVEWT.c=sum(LIVEWT.c))
        Shark.disc=Shark.disc%>%
          left_join(Drop,'FINYEAR')%>%
          mutate(LIVEWT.c=ifelse(is.na(LIVEWT.c),0,LIVEWT.c),
                 Shark.ktch=Shark.ktch-LIVEWT.c)%>%
          dplyr::select(-LIVEWT.c)
      }
    }
    
    #Calculate discards from observed composition
    if(Fishry %in% Reaportion.from.survey)
    {
      if(Fishry%in%c("PFT")) Prop=Pilbara.trawl_shrk.to.land
      if(Fishry%in%c("SWT","SCT")) Prop=South.west.trawl_shrk.to.land
      if(Fishry%in%c("SBP","SBSC")) Prop=Shark.Bay.trawl_shrk.to.land
      if(Fishry%in%c("EGP","OP","KP","NBP")) Prop=Exmouth.Onslow.trawl_shrk.to.land
        
      Shark.disc=Landings%>%mutate(Shark.ktch=KG*Prop*BRD)
      
      #remove catch if it's reported 
      if(length(yrs.reported.ktch)>0)
      {
        Drop=Data.monthly.north%>%
                filter(FishCubeCode%in%Fishry & FINYEAR%in% yrs.reported.ktch)%>%
                group_by(FINYEAR)%>%
                summarise(LIVEWT.c=sum(LIVEWT.c))
        Shark.disc=Shark.disc%>%
                    left_join(Drop,'FINYEAR')%>%
                    mutate(LIVEWT.c=ifelse(is.na(LIVEWT.c),0,LIVEWT.c),
                           Shark.ktch=Shark.ktch-LIVEWT.c)%>%
              dplyr::select(-LIVEWT.c)
      }
      
      
      #Expand to species using observed composition
      if(Fishry%in%c("PFT")) Survey=Pilbara.trawl.observed.comp
      if(Fishry%in%c("SWT","SCT")) Survey=South.west.trawl.observed.comp
      if(Fishry%in%c("SBP","SBSC")) Survey=Shark.Bay.trawl.observed.comp
      if(Fishry%in%c("EGP","OP","KP","NBP")) Survey=Exmouth.Onslow.trawl.observed.comp
      
      Survey=Survey%>%filter(!is.na(SPECIES))%>%
        dplyr::select(SPECIES,Prop)%>%
        group_by(SPECIES)%>%
        summarise(Prop=sum(Prop))%>%
        mutate(Prop=Prop/sum(Prop))

      N.dat=nrow(Shark.disc)
      NN=rep(1:N.dat,each=nrow(Survey))
      Shark.disc=Shark.disc[NN,]%>%
                  mutate(SPECIES=rep(Survey$SPECIES,N.dat))%>%
                  left_join(Survey,by='SPECIES')%>%
                  mutate(Shark.ktch=Shark.ktch*Prop)%>%
                  dplyr::select(-Prop)
    }
    
    # Add discards to Data.mmonthly or Data.monthly.north
    if(Fishry%in%c("SBP","SBS")) which.data.set="Data.monthly.north" else
    {
      if(Fishry%in%unique(Data.monthly$FishCubeCode))which.data.set="Data.monthly"
      if(Fishry%in%unique(Data.monthly.north$FishCubeCode))which.data.set="Data.monthly.north"
    }
    
    if(which.data.set=="Data.monthly")
    {
      Dumy=Data.monthly[1:nrow(Shark.disc),]
      Dumy[,]=NA
      Dumy=Dumy%>%
        mutate(LIVEWT.c=Shark.disc$Shark.ktch,
               FINYEAR=Shark.disc$FINYEAR,
               FishCubeCode=Fishry,
               Discarded.ktch="YES")
      
      if(Fishry %in% Reaportion.from.reported.ktch)
      {
        Dumy=Dumy%>%mutate(SPECIES=22999,
                           SNAME="SHARK, OTHER")
      }
      if(Fishry %in% Reaportion.from.survey)
      {
        Dumy=Dumy%>%mutate(SPECIES=Shark.disc$SPECIES)
      }
      a1=Data.monthly%>%filter(FishCubeCode==Fishry)%>%distinct(zone)
      Dumy$zone=names(rev(sort(table(a1$zone))))[1]
      Data.monthly=rbind(Data.monthly,Dumy)
    }
    if(which.data.set=="Data.monthly.north")
    {
      Dumy=Data.monthly.north[1:nrow(Shark.disc),]
      Dumy[,]=NA
      Dumy=Dumy%>%
        mutate(LIVEWT.c=Shark.disc$Shark.ktch,
               FINYEAR=Shark.disc$FINYEAR,
               FishCubeCode=Fishry,
               Discarded.ktch="YES")
      
      if(Fishry %in% Reaportion.from.reported.ktch)
      {
        Dumy=Dumy%>%mutate(SPECIES=22999,
                           SNAME="SHARK, OTHER")
      }
      if(Fishry %in% Reaportion.from.survey)
      {
        Dumy=Dumy%>%mutate(SPECIES=Shark.disc$SPECIES)
      }
      a1=Data.monthly.north%>%filter(FishCubeCode==Fishry)%>%distinct(zone)
      Dumy$zone=names(rev(sort(table(a1$zone))))[1]
      Data.monthly.north=rbind(Data.monthly.north,Dumy)
    }
  }
}


  #3.1.2 Reapportion catch of 'shark,other' for non-shark fisheries      
#note:  for TDGDLF, main 4 species already reapportioned
for(f in 1:length(Lista.reap.FishCubeCode))   
{
  Fishry=Lista.reap.FishCubeCode[[f]]
  
  #get fishery
  South=Data.monthly%>%filter(FishCubeCode==Fishry & SPECIES<50000)
  North=Data.monthly.north%>%filter(FishCubeCode==Fishry & SPECIES<50000)
  
  #Reapportion using reported catch composition
  if(Fishry %in% Reaportion.from.reported.ktch)
  {
    if(nrow(South)>0)
    {
      Shark.other=South%>%
        filter(SPECIES==22999)%>%
        mutate(SNAME='',
               RSCommonName='',
               RSSpeciesId='')
      South=South%>%filter(!SPECIES==22999)
      
      if(sum(!South$SPECIES==South$Spec.old,na.rm=T)>0) #remove indicator species if already reapportioned
      {
        Survey=South%>%filter(!SPECIES%in%c(17001,18003,17003,18007))%>%
              group_by(SPECIES)%>%
              summarise(Total=sum(LIVEWT.c))%>%
              mutate(Prop=Total/sum(Total))%>%
              dplyr::select(-Total)%>%
              data.frame
      }else
      {
        Survey=South%>%group_by(SPECIES)%>%
          summarise(Total=sum(LIVEWT.c))%>%
          mutate(Prop=Total/sum(Total))%>%
          dplyr::select(-Total)%>%
          data.frame
      }
      
      N.dat=nrow(Shark.other)
      NN=rep(1:N.dat,each=nrow(Survey))
      Shark.other=Shark.other[NN,]%>%
        mutate(SPECIES=rep(Survey$SPECIES,N.dat))%>%
        left_join(Survey,by='SPECIES')%>%
        mutate(LIVEWT=LIVEWT*Prop,
               LIVEWT.c=LIVEWT.c*Prop,
               LIVEWT.orgnl=LIVEWT.orgnl*Prop,
               LIVEWT.reap=LIVEWT.reap*Prop)%>%
        dplyr::select(-Prop)
      South=rbind(South,Shark.other)
    }
    if(nrow(North)>0)
    {
      Shark.other=North%>%
        filter(SPECIES==22999)%>%
        mutate(SNAME='',
               RSCommonName='',
               RSSpeciesId='')
      North=North%>%filter(!SPECIES==22999)
      
      Survey=North%>%group_by(SPECIES)%>%
        summarise(Total=sum(LIVEWT.c))%>%
        mutate(Prop=Total/sum(Total))%>%
        dplyr::select(-Total)%>%
        data.frame
      
      N.dat=nrow(Shark.other)
      NN=rep(1:N.dat,each=nrow(Survey))
      Shark.other=Shark.other[NN,]%>%
        mutate(SPECIES=rep(Survey$SPECIES,N.dat))%>%
        left_join(Survey,by='SPECIES')%>%
        mutate(LIVEWT=LIVEWT*Prop,
               LIVEWT.c=LIVEWT.c*Prop,
               LIVEWT.orgnl=LIVEWT.orgnl*Prop,
               LIVEWT.reap=LIVEWT.reap*Prop)%>%
        dplyr::select(-Prop)
      North=rbind(North,Shark.other)
    }
  }
  #Reapportion using observer catch composition
  if(Fishry %in% Reaportion.from.survey)
  {
    if(Fishry=="PFT") Survey=Pilbara.trawl.observed.comp
    if(Fishry%in%c("SWT","SCT")) Survey=South.west.trawl.observed.comp
    if(Fishry%in%c("SBP","SBSC")) Survey=Shark.Bay.trawl.observed.comp
    if(Fishry%in%c("EGP","OP","KP","NBP")) Survey=Exmouth.Onslow.trawl.observed.comp
    Survey=Survey%>%filter(!is.na(SPECIES))%>%
      dplyr::select(SPECIES,Prop)%>%
      group_by(SPECIES)%>%
      summarise(Prop=sum(Prop))%>%
      mutate(Prop=Prop/sum(Prop))
    if(nrow(South)>0)
    {
      Shark.other=South%>%
        filter(SPECIES==22999)%>%
        mutate(SNAME='',
               RSCommonName='',
               RSSpeciesId='')
      N.dat=nrow(Shark.other)
      NN=rep(1:N.dat,each=nrow(Survey))
      Shark.other=Shark.other[NN,]%>%
        mutate(SPECIES=rep(Survey$SPECIES,N.dat))%>%
        left_join(Survey,by='SPECIES')%>%
        mutate(LIVEWT=LIVEWT*Prop,
               LIVEWT.c=LIVEWT.c*Prop,
               LIVEWT.orgnl=LIVEWT.orgnl*Prop,
               LIVEWT.reap=LIVEWT.reap*Prop)%>%
        dplyr::select(-Prop)
      
      South=South%>%filter(!SPECIES==22999)
      South=rbind(South,Shark.other)
    }
    if(nrow(North)>0)
    {
      Shark.other=North%>%
        filter(SPECIES==22999)%>%
        mutate(SNAME='',
               RSCommonName='',
               RSSpeciesId='')
      N.dat=nrow(Shark.other)
      NN=rep(1:N.dat,each=nrow(Survey))
      Shark.other=Shark.other[NN,]%>%
        mutate(SPECIES=rep(Survey$SPECIES,N.dat))%>%
        left_join(Survey,by='SPECIES')%>%
        mutate(LIVEWT=LIVEWT*Prop,
               LIVEWT.c=LIVEWT.c*Prop,
               LIVEWT.orgnl=LIVEWT.orgnl*Prop,
               LIVEWT.reap=LIVEWT.reap*Prop)%>%
        dplyr::select(-Prop)
      
      North=North%>%filter(!SPECIES==22999)
      North=rbind(North,Shark.other)
    }
  }
  
  #Remove fishery 
  if(nrow(South)>0) Data.monthly=Data.monthly%>%filter(!FishCubeCode==Fishry & SPECIES<50000)
  if(nrow(North)>0) Data.monthly.north=Data.monthly.north%>%filter(!FishCubeCode==Fishry & SPECIES<50000)
  
  #Add fishery back with species reapportioned
  if(nrow(South)>0) Data.monthly=rbind(Data.monthly,South)
  if(nrow(North)>0) Data.monthly.north=rbind(Data.monthly.north,North)
  
  rm(North,South)
}                                     

  #final reaportion of 'shark,other' for a few fisheries not meeting criterion
fn.reap.shkr.oder=function(d)
{
  ThiS=d%>%filter(SPECIES==22999)%>%distinct(FishCubeCode)%>%pull(FishCubeCode)
  ADD=vector('list',length(ThiS))
  for(x in 1:length(ThiS))
  {
    dd=d%>%filter(FishCubeCode%in%ThiS[x])
    if(length(unique(dd$SPECIES))>1)
    {
      Shark.other=dd%>%
        filter(SPECIES==22999)%>%
        mutate(SNAME='',
               RSCommonName='',
               RSSpeciesId='')
      dd=dd%>%filter(!SPECIES==22999)
      
      Survey=dd%>%group_by(SPECIES)%>%
        summarise(Total=sum(LIVEWT.c))%>%
        mutate(Prop=Total/sum(Total))%>%
        dplyr::select(-Total)%>%
        data.frame
      
      N.dat=nrow(Shark.other)
      NN=rep(1:N.dat,each=nrow(Survey))
      Shark.other=Shark.other[NN,]%>%
        mutate(SPECIES=rep(Survey$SPECIES,N.dat))%>%
        left_join(Survey,by='SPECIES')%>%
        mutate(LIVEWT=LIVEWT*Prop,
               LIVEWT.c=LIVEWT.c*Prop,
               LIVEWT.orgnl=LIVEWT.orgnl*Prop,
               LIVEWT.reap=LIVEWT.reap*Prop)%>%
        dplyr::select(-Prop)
      dd=rbind(dd,Shark.other)
    }
    ADD[[x]]=dd
  }
  ADD=do.call(rbind,ADD)
  d=d%>%filter(!FishCubeCode%in%ThiS)
  return(rbind(d,ADD))
}
Data.monthly=fn.reap.shkr.oder(d=Data.monthly)
Data.monthly.north=fn.reap.shkr.oder(d=Data.monthly.north)

  #reapportion of hammerheads
fn.reap.hh=function(d,Survey)
{
  dd=d%>%filter(SPECIES%in%19000:19100)
  
  Shark.other=dd%>%
    filter(SPECIES==19000)%>%
    mutate(SNAME='',
           RSCommonName='',
           RSSpeciesId='')
  dd=dd%>%filter(!SPECIES==19000)
  
  N.dat=nrow(Shark.other)
  NN=rep(1:N.dat,each=nrow(Survey))
  Shark.other=Shark.other[NN,]%>%
    mutate(SPECIES=rep(Survey$SPECIES,N.dat))%>%
    left_join(Survey,by='SPECIES')%>%
    mutate(LIVEWT=LIVEWT*Prop,
           LIVEWT.c=LIVEWT.c*Prop,
           LIVEWT.orgnl=LIVEWT.orgnl*Prop,
           LIVEWT.reap=LIVEWT.reap*Prop)%>%
    dplyr::select(-c(Prop,Name))
  dd=rbind(dd,Shark.other)
  
  d=d%>%filter(!SPECIES%in%19000:19100)
  return(rbind(d,dd))
  
}
Data.monthly=fn.reap.hh(d=Data.monthly,Survey=Comp.hh.south)
Data.monthly.north=fn.reap.hh(d=Data.monthly.north,Survey=Comp.hh.north)


  #reapportion blacktip sharks reported too far south
#note: reported all the way to Esperance (doesn't conform to the species distribution / observer data) 
#      Hence set to spinner shark any blacktip record south of 32 (Harry et al 2019 as 
#      evidence of southernmost distribution, though for east coast)
Data.monthly=Data.monthly%>%
  mutate(SNAME=ifelse(SPECIES==18014 & LAT<(-31.5) & LONG >114,
                      "SHARK, SPINNER (LONG-NOSE GREY)",SNAME),
         SPECIES=ifelse(SPECIES==18014 & LAT<(-31.5) & LONG >114,
                        18023,SPECIES))


  #3.1.3 Apply fishery and species specific PCM to discarded catch
fn.x=function(x)as.numeric(x[1]):as.numeric(x[2])
PCM.sp=list(Sawfish=fn.x(c( 25000,25020)),
            Wobbegongs=fn.x(c(13000,13020)),
            Mackerel=fn.x(c(10000,10020)),
            Greynurse=fn.x(c(8000,8020)),
            Hexanchids=fn.x(c(5001,5020)),
            Whalers=fn.x(c(18000,18090)),
            Hammerheads=fn.x(c(19000,19020)),
            Triakids=fn.x(c(17000,17020)),
            Angels=fn.x(c(24900,24920)),
            Dogfish=fn.x(c(20000,20020)),
            Sawsharks=fn.x(c(23000,23920)),
            Guitarfish=fn.x(c(26990,26999)),
            Shark.other=22999,
            Ray.other=31000)
PCM.sp=unlist(PCM.sp)
names(PCM.sp)=removeNumbers(names(PCM.sp))
All.rays=c("Numbfish","Rajids","Eagle.rays","Other.rays","Dasyatids")  
All.shrks=PCM$Group[-match(c("Sawfish","Guitarfish",All.rays),PCM$Group)]
PCM.sp=data.frame(Group=names(PCM.sp),SPECIES=PCM.sp)%>%
  left_join(PCM,by="Group")%>%
  mutate(Trawl=ifelse(Group=='Shark.other',mean(PCM%>%filter(Group%in%All.shrks)%>%pull(Trawl),na.rm=T),
               ifelse(Group=='Ray.other',mean(PCM%>%filter(Group%in%All.rays)%>%pull(Trawl),na.rm=T),Trawl)),
         GN=ifelse(Group=='Shark.other',mean(PCM%>%filter(Group%in%All.shrks)%>%pull(GN),na.rm=T),
            ifelse(Group=='Ray.other',mean(PCM%>%filter(Group%in%All.rays)%>%pull(GN),na.rm=T),GN)),
         LL=ifelse(Group=='Shark.other',mean(PCM%>%filter(Group%in%All.shrks)%>%pull(LL),na.rm=T),
            ifelse(Group=='Ray.other',mean(PCM%>%filter(Group%in%All.rays)%>%pull(LL),na.rm=T),LL)))
fn.PCM=function(d)
{
  d.disc=d%>%filter(Discarded.ktch=='YES')
  d=d%>%filter(!Discarded.ktch=='YES')
  
  d.disc=d.disc%>%
      mutate(what.method=ifelse(FishCubeCode%in%c('C156','PFT','KP','SBSC','AIMWT','SBP','SCT','SWT','EGP',
                                                  'NBP','OP'),"Trawl",
                         ifelse(FishCubeCode%in%c('SWBN','MBC','CSFN','SCE','WCE','C066','C070',
                                                  'JASDGDL','OANCGCWC','KGB','SBBS','EGBS'),"GN",
                         ifelse(FishCubeCode%in%c('C019','CSC','C048','C074','C111','PT','NDS','KTR',
                                                  'OASC','SBS','CSLP','WCDGDL','WANCS'),"LL",NA))))%>%
      left_join(PCM.sp,by="SPECIES")%>%
      mutate(LIVEWT.c=ifelse(what.method=='Trawl',LIVEWT.c*Trawl,
                      ifelse(what.method=='GN',LIVEWT.c*GN,
                      ifelse(what.method=='LL',LIVEWT.c*LL,
                      LIVEWT.c))))
  d.disc=d.disc%>%dplyr::select(-c(what.method,Group,Trawl,GN,LL))
  return(rbind(d,d.disc))
}
Data.monthly=fn.PCM(d=Data.monthly)
Data.monthly.north=fn.PCM(d=Data.monthly.north)


  #3.1.4. Kimberley Gillnet and Barramundi recalculation of catch (in kg) using effort as per McAuley et al 2005
Kimberley.GBF.observed.comp=Kimberley.GBF.observed.comp%>%
          mutate(cpue=Weight/Kimberley.GBF.observed.effort)
fn.cpue.to.ktch=function(cpue,annual.effrt)
{
  spi=unique(cpue$Common.name)
  d.list=vector('list',length(spi))
  for(d in 1:length(spi))
  {
    x1=cpue[d,]
    x2=annual.effrt%>%
      mutate(total.catch=Total.Bdays*x1$cpue,
             Common.name=x1$Common.name,
             SPECIES=x1$SPECIES)
    d.list[[d]]=x2
  }
  return(do.call(rbind,d.list))
}
KGB.tot.ktch=fn.cpue.to.ktch(cpue=Kimberley.GBF.observed.comp,
                             annual.effrt=Kimberley.GBF.annual.effort%>%
                               dplyr::select(Year,Total.Bdays))
Tab.zon.KGB=table(Data.monthly.north%>%filter(FishCubeCode=='KGB')%>%pull(zone))
Data.monthly.north=Data.monthly.north%>%filter(!FishCubeCode=='KGB')
add.KGBM.to.Data.monthly.north=Data.monthly.north[1:nrow(KGB.tot.ktch),]
add.KGBM.to.Data.monthly.north[,]=NA

add.KGBM.to.Data.monthly.north=add.KGBM.to.Data.monthly.north%>%
                        mutate(LIVEWT.c=KGB.tot.ktch$total.catch,
                               SPECIES=KGB.tot.ktch$SPECIES,
                               FINYEAR=with(KGB.tot.ktch,paste(Year,substring(Year+1,3,4),sep='-')),
                               FishCubeCode='KGB',
                               METHOD="GN",
                               LAT=-16,   #mean latitude
                               Estuary="NO",
                               zone=sample(names(Tab.zon.KGB),nrow(KGB.tot.ktch),replace=T,prob=Tab.zon.KGB/sum(Tab.zon.KGB)))   #randomly allocate zone by observed records
Data.monthly.north=rbind(Data.monthly.north,add.KGBM.to.Data.monthly.north)
rm(KGB.tot.ktch,add.KGBM.to.Data.monthly.north)


  #3.1.5. Combine historic in single dataframe
Catch_1949.total=data.frame(year=1949,LIVEWT.c=sum(Catch_1949$LIVEWT.c))     
Catch_1949.total$LIVEWT.c=Catch_1950$LIVEWT.c     
Historic.ktch=rbind(Catch_1941_1943,Catch_1949.total,Catch_1950,Catch_1952_1975)
if(Do.recons.paper=="YES")
{
  tiff(file=fn.hnd.out("Figure_Reported.pre-1975.annual.shark.catch.tiff"),2400,2400,
       units="px",res=300,compression="lzw")
  par(mfcol=c(1,1),las=1,mgp=c(2.9,.6,0))
  with(subset(Historic.ktch,year<1975),
       {
         plot(year,LIVEWT.c, ylab="Total shark landings (tonnes)",xlab="Financial year",
              pch=19,col=1,cex=2,cex.lab=1.7,cex.axis=1.35)
       })
  dev.off()
}
  #add missing years through linear interpolation
Missing.yrs=Historic.yrs[which(!Historic.yrs%in%Historic.ktch$year)]
Missing.ktch=approx(Historic.ktch$year,Historic.ktch$LIVEWT.c,xout=Missing.yrs)   
Historic.ktch=rbind(Historic.ktch,data.frame(year=Missing.ktch$x,LIVEWT.c=Missing.ktch$y))
Historic.ktch=Historic.ktch[order(Historic.ktch$year),]
Historic.ktch$LIVEWT.c=Historic.ktch$LIVEWT.c*1000   #convert back to kg

  #get prop of catch for 1975-1980 to extract historic
ALL=subset(Data.monthly,SPECIES<32000 & METHOD%in%c("GN","LL") & LAT<=(-26) & Estuary=="NO"
           & FINYEAR%in%c("1975-76","1976-77","1977-78","1978-79","1979-80","1980-81"),
           select=c(FINYEAR,LIVEWT.c,SPECIES,SNAME))%>%
                  filter(!SPECIES==22999)
ALL.north=subset(Data.monthly.north,SPECIES<32000 & METHOD%in%c("GN","LL") & LAT>(-26) & 
                  FINYEAR%in%c("1975-76","1976-77","1977-78","1978-79","1979-80","1980-81"),
                  select=c(FINYEAR,LIVEWT.c,SPECIES,SNAME))%>%
                    filter(!SPECIES==22999)
Prop.sp.yr.all=rbind(ALL,ALL.north)%>%group_by(SPECIES)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c))%>%
  mutate(Proportion=LIVEWT.c/sum(LIVEWT.c))%>%
  dplyr::select(SPECIES,Proportion)%>%
  data.frame

#reapportion historic catch to species
Hist.expnd=expand.grid(year=Historic.ktch$year,SPECIES=Prop.sp.yr.all$SPECIES)
Hist.expnd=Hist.expnd%>%
  left_join(Historic.ktch,by="year")%>%
  left_join(Prop.sp.yr.all,by='SPECIES')%>%
  rename(ktch=LIVEWT.c)%>%
  mutate(LIVEWT.c=ktch*Proportion,
         FINYEAR=paste(year,substr(year+1,3,4),sep='-'))%>%
  dplyr::select(FINYEAR,SPECIES,LIVEWT.c)


  #3.1.6. TEPS    
#TDGDLF
Size.comp.Dusky.TEPS_TDGLDF=c(3.05,4,3,3.5,3.5,3.5,3.5,3,3,3,4,3)     #from Comments in TDGDLF returns
Size.comp.Greynurse.TEPS_TDGLDF=c(rep(.61,5),rep(.94,10),rep(1.22,10),
                                  rep(1.52,7),rep(1.83,7),rep(2.1,2),
                                  rep(2.4,6),rep(3.05,5),rep(1,2,2),rep(3,5),2.5,rep(1.5,3))     
TEPS=TEPS%>%
        filter(SpeciesCode%in%c(18003,8001))%>%
        left_join(PCM.sp,by=c("SpeciesCode"= "SPECIES"))%>%
        rename(FINYEAR=finyear,
               SPECIES=SpeciesCode)%>%      
        mutate(Kg=ifelse(SPECIES==18003,bwt*(100*mean(Size.comp.Dusky.TEPS_TDGLDF))^awt,
                  ifelse(SPECIES==8001,bwt.grey*mean(Size.comp.Greynurse.TEPS_TDGLDF)^awt.grey,
                         NA)),
               LIVEWT.c=ifelse(Status%in%c("D","d"),Kg*Number,
                        ifelse(Status%in%c("A","a"),Kg*Number*GN,
                               NA)))%>%
        group_by(FINYEAR,SPECIES)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
        data.frame
Greynurse.ktch=TEPS%>%filter(SPECIES==8001)
TEPS_dusky=TEPS%>%filter(SPECIES==18003)


  # 3.1.7 Wetline in the Western Rock lobster fishery
#note: extrapolate Wann's catch (for one season) to the entire fishery as a proportion of the number of boats
WRL.Wann1=WRL.Wann%>%
    mutate(TL=TL_metres*100,
           TL=ifelse(is.na(TL),TL_feet*30.48,TL),
           Species=ifelse(Species=="DW",
                         sample(c("CP","BW"),55,
                                prob=c(WRL.copper.dusky.prop,1-WRL.copper.dusky.prop), replace = TRUE),
                         Species))%>%   #convert to cm
    left_join(Weight,by="Species")%>%
    mutate(N=1,
           LIVEWT.c=(bwt*TL^awt),
           LIVEWT.c=ifelse(!Species=='BW',LIVEWT.c/1000,LIVEWT.c))%>%
    group_by(Species,CAES_Code)%>%
    summarise(N=sum(N,na.rm=T))

WRL.total.ktch=WRL[rep(1:nrow(WRL),each=nrow(WRL.Wann1)),]%>%
                mutate(Species=rep(WRL.Wann1$Species,nrow(WRL)))%>%
                left_join(WRL.Wann1,by="Species")%>%
                mutate(LIVEWT.c=N*WRL.assumed.weight*Number.of.vessels*WRL.prop)%>%
                dplyr::select(Finyear,CAES_Code,LIVEWT.c)%>%
                rename(FINYEAR=Finyear,
                       SPECIES=CAES_Code)



## 3.2. Catch of non WA Fisheries

  #-- 3.2.1 Taiwanese gillnet and longline
Mn.ktch=mean(Taiwan.gillnet.ktch$WA.Ktch/Taiwan.gillnet.ktch$Shark.Ktch,na.rm=T)
Taiwan.gillnet.ktch=Taiwan.gillnet.ktch%>%
  mutate(Shark.Ktch=ifelse(is.na(Shark.Ktch),Total.Ktch*mean(Shark.Ktch/Total.Ktch,na.rm=T),
                           Shark.Ktch),
         WA.Ktch=ifelse(is.na(WA.Ktch),Shark.Ktch*mean(WA.Ktch/Shark.Ktch,na.rm=T),
                        WA.Ktch))

Taiwan.gillnet.ktch=cbind(Taiwan.gillnet.ktch,Taiwan.gillnet.sp.comp) %>%
  mutate('Australian blacktip shark'=Australian.blacktip.Shark.West.prop*WA.Ktch,
         'Spot-tail shark'=Spot.tail.Shark.West.prop*WA.Ktch,
         Hammerheads=Hammerheads.West.prop*WA.Ktch)%>%
  dplyr::select(Year,'Australian blacktip shark','Spot-tail shark',Hammerheads)%>%
  gather(Species,Ktch,-Year)

Taiwan.longline.ktch=Taiwan.longline.ktch%>%
  mutate(WA.ktch=Shark.Ktch*Mn.ktch)
Taiwan.longline.ktch=Taiwan.longline.sp.comp%>%
  left_join(Taiwan.longline.ktch,by="Year")%>%
  mutate(Ktch=WA.ktch*Percent/100)%>%
  dplyr::select(names(Taiwan.gillnet.ktch))

    #add SPECIES
Taiwan.gillnet.ktch=Taiwan.gillnet.ktch%>%  
              mutate(SPECIES=ifelse(Species=="Australian blacktip shark",18014,
                             ifelse(Species=="Spot-tail shark",18013,
                             ifelse(Species=="Hammerheads",19000,
                                    22999))))
Taiwan.longline.ktch=Taiwan.longline.ktch%>%  
              mutate(SPECIES=ifelse(Species=="Australian blacktip shark",18014,
                             ifelse(Species=="Spot-tail shark",18013,
                             ifelse(Species=="Hammerheads",19000,
                             ifelse(Species=="Tiger shark",18022,
                             ifelse(Species=="Milk shark",18006,
                             ifelse(Species=="Spinner shark",18023,
                             ifelse(Species=="Pigeye shark",18026,
                             ifelse(Species=="Graceful shark",18033,
                                    22999)))))))))  
  
    #reapportion of hammerhead catch    
fn.reap.hh.Taiwan=function(d,Survey)
{
  d=d%>%
    rename(LIVEWT.c=Ktch)%>%
    mutate(FINYEAR=paste(Year,substr(Year+1,3,4),sep='-'))%>%
    dplyr::select(c(FINYEAR,SPECIES,LIVEWT.c))
  
  dd=d%>%filter(SPECIES%in%19000:19100)
  N.dat=nrow(dd)
  NN=rep(1:N.dat,each=nrow(Survey))
  dd=dd[NN,]%>%
    mutate(SPECIES=rep(Survey$SPECIES,N.dat))%>%
    left_join(Survey,by='SPECIES')%>%
    mutate(LIVEWT.c=LIVEWT.c*Prop)%>%
    dplyr::select(c(FINYEAR,SPECIES,LIVEWT.c))

  d=d%>%filter(!SPECIES%in%19000:19100)
  return(rbind(d,dd))
  
}
Taiwan.gillnet.ktch=fn.reap.hh.Taiwan(d=Taiwan.gillnet.ktch,Survey=Comp.hh.north)
Taiwan.longline.ktch=fn.reap.hh.Taiwan(d=Taiwan.longline.ktch,Survey=Comp.hh.north)


  #-- 3.2.2 Commonwealth GAB trawl and Western Tuna and Billfish Fisheries (WTBF)  
#WTBF_catch GAB.trawl_catch     #by species, add SPECIES, etc
WTBF_catch=WTBF_catch%>%
         rename(LIVEWT.c=Catch_kg,
               FINYEAR=Year)%>%
         mutate(SPECIES=
                 ifelse(Species=="Blacktip sharks",18014,
                 ifelse(Species=="Blue shark",18004,
                 ifelse(Species=="Hammerheads",19000,
                 ifelse(Species=="Tiger shark",18022,
                 ifelse(Species=="Bronze whaler",18001,
                 ifelse(Species=="Dusky shark",18003,
                 ifelse(Species=="Scalloped hammerhead",19001,
                 ifelse(Species=="Sandbar shark",18007,
                 ifelse(Species=="Shortfin mako",10001,
                 22999))))))))))%>%
  dplyr::select(-Species)
GAB.trawl_catch=GAB.trawl_catch%>%
          rename(LIVEWT.c=Catch_kg,
                 FINYEAR=Finyear)%>%
          mutate(SPECIES=
           ifelse(Species%in%c("Broadnose sevengill shark","Sevengilled shark"),5001,
           ifelse(Species=="Blue shark",18004,
           ifelse(Species=="Common saw shark",23002,
           ifelse(Species=="Saw sharks",23900,
           ifelse(Species=="Bronze whaler",18001,
           ifelse(Species=="Greynurse shark",8001,
           ifelse(Species=="Scalloped hammerhead",19001,
           ifelse(Species=="Gummy shark",17001,
           ifelse(Species=="Shortfin mako",10001,
           ifelse(Species=="School shark",17008,
           ifelse(Species=="Wobbegong sharks",13000,
           ifelse(Species=="Smooth hammerhead",19004,
           ifelse(Species=="Whiskery shark",17003,
           22999))))))))))))))%>%
  dplyr::select(-Species)


  #-- 3.2.3 SA Marine Scalefish fishery     
  #split dusky from bronzy                      
Whaler_SA=Whaler_SA%>%
          rename(LIVEWT.c=Live_wt.tons.)%>%
          mutate(LIVEWT.c=LIVEWT.c*1000)    #convert to kg
N.SA=nrow(Whaler_SA)
Whaler_SA=Whaler_SA[rep(1:N.SA,each=2),]%>%
  mutate(SPECIES=rep(c(18003,18001),N.SA),
         LIVEWT.c=ifelse(SPECIES==18003,
                         LIVEWT.c*Whaler_SA_dusky.prop,
                         LIVEWT.c*Whaler_SA_bronzie.prop))


#-- 3.2.4 Indonesian illegal fishing in Australia waters             

Indo_apprehensions=rbind(Indo_apprehensions.Stacey,Indo_apprehensions)
Missn.appr=which(is.na(Indo_apprehensions$Apprehensions))
if(length(Missn.appr)>0)
{
  Miss.appr.yrs=Indo_apprehensions$year[Missn.appr]
  Missing.appre=with(Indo_apprehensions,approx(year,Apprehensions,xout=Miss.appr.yrs))
  Indo_apprehensions$Apprehensions[Missn.appr]=Missing.appre$y
}
 
# Combine Vanesa's and Marshall et al data on catch compo
Indo_MOU.Vanesa=Indo_MOU.Vanesa%>%
  filter(Location%in%c('MOU Box',"MoU Box"))%>%
  rename(Species=Species.Ref_plus.Genetics)%>%
  dplyr::select(Species)%>%
  mutate(Species=paste(Species,"shark"),
         Species=ifelse(Species=="Hammerhead_unsp shark","Hammerheads",
                        ifelse(Species=="Guitarfish_unsp shark","Guitarfish",
                               ifelse(Species=="Blacktip_unsp shark","Blacktips",Species))))
Indo_MOU.Vanesa=table(Indo_MOU.Vanesa)   
Indo_MOU.Vanesa=Indo_MOU.Vanesa/sum(Indo_MOU.Vanesa)
Indo_MOU.Vanesa=data.frame(Species=names(Indo_MOU.Vanesa),
                           Proportion.by.number=c(Indo_MOU.Vanesa))%>%
  mutate(Species=as.character(Species))
rownames(Indo_MOU.Vanesa)=NULL
Indo.prop.ratio=Indo_shark.comp%>%
            mutate(ratio=Proportion/Proportion.by.number,
                   Species=as.character(Species))%>%
            dplyr::select(Species,ratio)
Indo_MOU.Vanesa=left_join(Indo_MOU.Vanesa,Indo.prop.ratio,by='Species')%>%
                  mutate(ratio=ifelse(is.na(ratio),1,ratio),
                         Proportion.Vanesa=ratio*Proportion.by.number)%>%
                  dplyr::select(Species,Proportion.Vanesa)

Indo_shark.comp=Indo_shark.comp%>%
                dplyr::select(-Proportion.by.number)%>%
                mutate(Species=as.character(Species))%>%
                full_join(Indo_MOU.Vanesa,by='Species')%>%
                mutate(Proportion=ifelse(is.na(Proportion),Proportion.Vanesa,
                                         Proportion),
                       Proportion.Vanesa=ifelse(is.na(Proportion.Vanesa),Proportion,
                                          Proportion.Vanesa))
Indo_shark.comp$Prop=rowMeans(Indo_shark.comp[,2:3])
Indo_shark.comp=Indo_shark.comp%>%
                  mutate(Prop=Prop/sum(Prop))%>%
          dplyr::select(Species,Prop)%>%
                  rename(Proportion=Prop)
#reapportion 'hammerhead'               
Indo_shark.comp=Indo_shark.comp%>%
                mutate(Proportion=ifelse(Species=='Great hammerhead',
                                         5.659801e-03+1.295742e-02+0.05597481,
                                   ifelse(Species=='Scalloped hammerhead',
                                          2.930968e-03+0.00881229,
                                      Proportion)))%>%
                filter(!Species%in%c('Great hammerhead shark','Hammerheads'))


#Indonesian. Calculate catch per vessel-day
Indo_avrg.ktch.per.vessel.day=Indo_shark.weight/(Indo_shark.N.vessels*7)
Indo_avrg.ktch.per.vessel.day=Indo_shark.comp%>%
  mutate(kg.day.vessel=Proportion*Indo_avrg.ktch.per.vessel.day)

#Indonesian. Calculate total catch for 2005
Indo_avrg.ktch.per.year=Indo_avrg.ktch.per.vessel.day%>%
  mutate(kg.year=kg.day.vessel*Indo_FFV_vessel.days.2005*Indo_Prop.Shark.vesl.Salini)

#Indonesian. Calculate Total catch by year
Indo_total.annual.ktch=Indo_avrg.ktch.per.year[rep(1:nrow(Indo_avrg.ktch.per.year),nrow(Indo_apprehensions)),]%>%
  mutate(year=rep(Indo_apprehensions$year,each=nrow(Indo_avrg.ktch.per.year)))%>%
  left_join(Indo_apprehensions,by='year')%>%
  mutate(Ap.prop.2005=Apprehensions/Indo_apprehensions$Apprehensions[which.max(Indo_apprehensions$Apprehensions)],
         LIVEWT.c=kg.year*Ap.prop.2005)  

Indo_total.annual.ktch=Indo_total.annual.ktch%>%
  mutate(Species=as.character(Species),
         FINYEAR=paste(year,substr(year+1,3,4),sep='-'))%>%
  left_join(All.species.names%>%dplyr::select(-Scien.nm),by=c('Species'='Name'))%>%
  dplyr::select(FINYEAR,SPECIES,LIVEWT.c)


  #MOU box. Calculate Total catch by year
Indo_total.annual.ktch_MOU=Indo_average.shark.per.trip_MOU*Indo_trips.per.year_MOU
Indo_total.annual.ktch_MOU=Indo_shark.comp%>%
  mutate(kg.year=Proportion*Indo_total.annual.ktch_MOU)

  #MOU box. Calculate Total catch by year
Indo_total.annual.ktch_MOU=Indo_total.annual.ktch_MOU%>%
          left_join(All.species.names%>%
          dplyr::select(-Scien.nm),by=c('Species'='Name'))
Indo_total.annual.ktch_MOU=Indo_total.annual.ktch_MOU[rep(1:nrow(Indo_total.annual.ktch_MOU),nrow(Indo_MOU.annual.trips.prop)),]%>%
          mutate(year=rep(Indo_apprehensions$year,each=nrow(Indo_total.annual.ktch_MOU)))%>%
          left_join(Indo_MOU.annual.trips.prop,by='year')%>%
          mutate(LIVEWT.c=kg.year*trips.prop,
                 Species=as.character(Species),
                 FINYEAR=paste(year,substr(year+1,3,4),sep='-'))%>%
  dplyr::select(colnames(Indo_total.annual.ktch))


  #Combine Indonesian illegal and legal Mou box
Indo_total.annual.ktch=rbind(Indo_total.annual.ktch,Indo_total.annual.ktch_MOU)%>%
              group_by(FINYEAR,SPECIES)%>%
              summarise(LIVEWT.c=sum(LIVEWT.c))%>%
              data.frame




# 4 -------------------EXPORT CATCH DATA------------------------------------
fn.out=function(d,NM)
{
  write.csv(d,paste('C:/Matias/Analyses/Data_outs/',NM,sep=""),row.names = F)
}
 
Yr.lim=c(as.numeric(substr(min(unique(Hist.expnd$FINYEAR)),1,4)),
         as.numeric(substr(Last.yr.ktch,1,4)))
This.fin.yr=paste(Yr.lim[1]:Yr.lim[2],substr((Yr.lim[1]+1):(Yr.lim[2]+1),3,4),sep='-')

#4.1 Catch_WA Fisheries

  #Historic
fn.out(d=Hist.expnd%>%filter(FINYEAR%in%This.fin.yr),NM='recons_Hist.expnd.csv')

  #Ammended reported catch including discards
fn.out(d=Data.monthly%>%filter(FINYEAR%in%This.fin.yr)%>%
         dplyr::select(FINYEAR,FishCubeCode,METHOD,BLOCKX,LAT,LONG,SPECIES,LIVEWT.c,
                       Estuary,zone),               
       NM='recons_Data.monthly.csv')
fn.out(d=Data.monthly.north%>%filter(FINYEAR%in%This.fin.yr)%>%
         dplyr::select(FINYEAR,FishCubeCode,METHOD,BLOCKX,LAT,LONG,SPECIES,LIVEWT.c,
                       Estuary,zone),                
       NM='recons_Data.monthly.north.csv')

  #TEPS
fn.out(d=Greynurse.ktch%>%filter(FINYEAR%in%This.fin.yr),NM='recons_Greynurse.ktch.csv')
fn.out(d=TEPS_dusky%>%filter(FINYEAR%in%This.fin.yr),NM='recons_TEPS_dusky.csv')

  #Wetline in the Western Rock lobster fishery
fn.out(d=WRL.total.ktch%>%filter(FINYEAR%in%This.fin.yr),NM='Wetline_rocklobster.csv')


#4.2. Catch of non WA Fisheries
 
   #Taiwanese gillnet and longline
fn.out(d=Taiwan.gillnet.ktch%>%filter(FINYEAR%in%This.fin.yr),NM='recons_Taiwan.gillnet.ktch.csv')
fn.out(d=Taiwan.longline.ktch%>%filter(FINYEAR%in%This.fin.yr),NM='recons_Taiwan.longline.ktch.csv')

  #Commonwealth GAB trawl and Western Tuna and Billfish Fisheries (WTBF) 
fn.out(d=GAB.trawl_catch%>%filter(FINYEAR%in%This.fin.yr),NM='recons_GAB.trawl_catch.csv')
fn.out(d=WTBF_catch%>%filter(FINYEAR%in%This.fin.yr),NM='recons_WTBF_catch.csv')
 
  #SA Marine Scalefish fishery
fn.out(d=Whaler_SA%>%filter(FINYEAR%in%This.fin.yr),NM='recons_Whaler_SA.csv')


  #Indonesian illegal fishing in Australia waters
fn.out(d=Indo_total.annual.ktch%>%filter(FINYEAR%in%This.fin.yr),NM='recons_Indo.IUU.csv')



# 5 -------------------REPORT SECTION------------------------------------
if(Do.recons.paper=="YES")   #for paper, report only IUU and reconstructions (not the TDGDLF or NSF)
{
  #Table 1. Species composition
  Hist.expnd$TYPE="Historic"
  Greynurse.ktch$TYPE="Protected"
  TEPS_dusky$TYPE="Protected"
  WRL.total.ktch$TYPE="WRL"
  Taiwan.gillnet.ktch$TYPE="Taiwanese"
  Taiwan.longline.ktch$TYPE="Taiwanese"
  
  Indo_total.annual.ktch$TYPE="Indonesian"
  Data.monthly.agg=Data.monthly%>%
                    group_by(FINYEAR,SPECIES)%>%
                    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
                    mutate(TYPE="WA.south")%>%
                    data.frame
  Data.monthly.north.agg=Data.monthly.north%>%
                    group_by(FINYEAR,SPECIES)%>%
                    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
                    mutate(TYPE="WA.north")%>%
                    data.frame
  
  Unified=rbind(Hist.expnd,Greynurse.ktch,TEPS_dusky,WRL.total.ktch,Indo_total.annual.ktch,
                Taiwan.gillnet.ktch,Taiwan.longline.ktch,
                Data.monthly.agg,Data.monthly.north.agg)%>%
              left_join(All.species.names,by='SPECIES')%>%
              filter(!is.na(Name) & FINYEAR%in%This.fin.yr) #report only complete years
    
  
  Table1=Unified%>%
            group_by(TYPE,SPECIES,Name)%>%
            summarise(Total=round(sum(LIVEWT.c)/1000,1))%>%
            spread(TYPE,Total,fill=0)%>%
            data.frame%>%
            mutate(Total=Historic+WA.south+WA.north+Protected+Taiwanese+Indonesian+WRL)%>%
            arrange(-Total)%>%
            left_join(All.species.names%>%dplyr::select(-SPECIES),by="Name")%>%
            rename(Common.name=Name,
                   South=WA.south,
                   North=WA.north,
                   Scientific.name=Scien.nm)
    
  
  source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")  
  setwd('C:/Matias/Analyses/Reconstruction_catch_commercial')
  fn.word.table(WD=getwd(),TBL=Table1%>%dplyr::select(Common.name,Scientific.name,Historic,South,North,Protected,Taiwanese,Indonesian),
                Doc.nm="Table1",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")

  
  #Figure 1. species annual catch by fishery
  source('C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R')
  Cum.ktch=cumsum(Table1$Total)/sum(Table1$Total)
  ID=which.min(abs(Cum.ktch-0.99))
  these.sp=Table1[1:ID,]%>%pull(Common.name)
  
  LWD=1.8
  CLs=data.frame(TYPE=c("Historic","South","North","Protected",
                        "Taiwanese","Indonesian","WRL"),
                 CL=c("black","deepskyblue2","coral3","green",
                      "dodgerblue4","darkorange","forestgreen"),
                 LT=c(1,1,1,1,3,3,1))
  fun.plt.Fig2=function(SP)
  {
    d=Unified%>%filter(Name==SP)%>%
                group_by(FINYEAR,TYPE)%>%
                summarise(Total=sum(LIVEWT.c)/1000)%>%
                spread(TYPE,Total,fill=NA)%>%
                data.frame
    yr=as.numeric(substr(d$FINYEAR,1,4))
    plot(yr,yr,xlim=Yr.lim,col='transparent',ylab="",xlab="",ylim=c(0,max(d[,-c(1)],na.rm=T)))
    with(CLs,
         {
           if("Historic"%in%names(d)) lines(yr,d$Historic,lwd=LWD,col=CL[1],lty=LT[1])
           if("WA.south"%in%names(d)) lines(yr,d$WA.south,lwd=LWD,col=CL[2],lty=LT[2])
           if("WA.north"%in%names(d)) lines(yr,d$WA.north,lwd=LWD,col=CL[3],lty=LT[3])
           if("Protected"%in%names(d)) lines(yr,d$Protected,lwd=LWD,col=CL[4],lty=LT[4])
           if("Taiwanese"%in%names(d)) lines(yr,d$Taiwanese,lwd=LWD,col=CL[5],lty=LT[5])
           if("Indonesian"%in%names(d)) lines(yr,d$Indonesian,lwd=LWD,col=CL[6],lty=LT[6])
           if("WRL"%in%names(d)) lines(yr,d$WRL,lwd=LWD,col=CL[7],lty=LT[7])
         })

    mtext(SP,3,cex=.85)
  }
  tiff(file=fn.hnd.out("Figure2_reconstructed.catch.by.species.tiff"),2400,2400,
       units="px",res=300,compression="lzw")
  smart.par(n.plots=length(these.sp),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(s in 1:length(these.sp)) fun.plt.Fig2(SP=these.sp[s])
  plot.new()
  legend("center",CLs$TYPE[1:4],col=CLs$CL[1:4],lty=CLs$LT[1:4],lwd=LWD*1.5,bty='n',cex=1.2)
  plot.new()
  legend("center",CLs$TYPE[5:7],col=CLs$CL[5:7],lty=CLs$LT[5:7],lwd=LWD*1.5,bty='n',cex=1.2)
  mtext("Financial year",1,outer=T,cex=1.15)
  mtext("Total catch (tonnes)",2,outer=T,las=3,cex=1.15)
  dev.off()


  #Difference between original and reconstructed   
  these.sp=these.sp   #combine hammerheads because original doesn't discriminate among species
  
  Data.monthly.original.agg=Data.monthly.original%>%
    group_by(FINYEAR,SPECIES)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
    mutate(TYPE="WA.south")%>%
    data.frame
  Data.monthly.north.original.agg=Data.monthly.north.original%>%
    group_by(FINYEAR,SPECIES)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
    mutate(TYPE="WA.north")%>%
    data.frame
  
  Data.monthly.original.agg=rbind(Data.monthly.original.agg,Data.monthly.north.original.agg)%>%
    group_by(FINYEAR,SPECIES)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
    data.frame%>%
    left_join(All.species.names,by='SPECIES')
  
  Unified2=Unified%>%
    mutate(Name=ifelse(Name%in%c("Smooth hammerhead","Scalloped hammerhead","Great hammerhead"),
                  "Hammerheads",Name))
  these.sp2=these.sp
  these.sp2=unique(ifelse(these.sp2%in%c("Smooth hammerhead","Scalloped hammerhead","Great hammerhead"),
                   "Hammerheads",these.sp2))

  
  fun.plt.Fig1=function(SP)
  {
    d=Unified2%>%filter(Name==SP)%>%
      group_by(FINYEAR)%>%
      summarise(Total=sum(LIVEWT.c)/1000)%>%
      data.frame
    yr=as.numeric(substr(d$FINYEAR,1,4))
    
    d1=Data.monthly.original.agg%>%filter(Name==SP)%>%
      group_by(FINYEAR)%>%
      summarise(Total=sum(LIVEWT.c)/1000)%>%
      data.frame
    yr1=as.numeric(substr(d1$FINYEAR,1,4))
    
    plot(yr,yr,xlim=Yr.lim,col='transparent',ylab="",xlab="",ylim=c(0,max(d[,-c(1)],na.rm=T)))
    lines(yr,d$Total,lwd=LWD,col='deepskyblue4',lty=1)
    lines(yr1,d1$Total,lwd=LWD,col='red',lty=3)
    mtext(SP,3,cex=.85)
  }
  tiff(file=fn.hnd.out("Figure1_reconstructed.vs.original.tiff"),2400,2400,
       units="px",res=300,compression="lzw")
  smart.par(n.plots=length(these.sp2),MAR=c(2,2,1,1),OMA=c(1.75,2,.5,.1),MGP=c(1,.5,0))
  for(s in 1:length(these.sp2)) fun.plt.Fig1(SP=these.sp2[s])
  plot.new()
  legend("right",c("Reconstructed","Original"),col=c('deepskyblue4','red'),lty=c(1,3),lwd=LWD*1.5,bty='n',cex=1.2)
  mtext("Financial year",1,outer=T,cex=1.15)
  mtext("Total catch (tonnes)",2,outer=T,las=3,cex=1.15)
  dev.off()
  
  #Indonesian illegal fishing in Western Australia waters
  # tiff(file=fn.hnd.out("Figure_WA_illegal.Indo.tiff"),2400,2400,units="px",res=300,compression="lzw")
  # par(mfcol=c(2,1),mar=c(4,4,1,.5),oma=c(2.2,7.2,1,1),las=1,mgp=c(2.5,.4,0))
  # D1=Indo_total.annual.ktch%>%
  #           group_by(FINYEAR)%>%
  #           summarise(Tot=sum(LIVEWT.c))%>%
  #           mutate(year=as.numeric(substr(FINYEAR,1,4)))
  # plot(D1$year,D1$Tot/1000,type='o',pch=19,cex=2,ylab="Total catch (tonnes)",xlab="Year")
  # D1=Indo_total.annual.ktch%>%group_by(SPECIES)%>%summarise(Tot=sum(LIVEWT.c,na.rm=T))%>%
  #   arrange(Tot)
  # barplot(D1$Tot/1000,horiz = T,xlab="Total catch (tonnes)",names.arg=D1$SPECIES,las=1)
  # dev.off()
}


