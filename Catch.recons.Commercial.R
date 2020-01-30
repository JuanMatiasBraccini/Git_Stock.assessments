#-- Script for reconstructing time series of commercial catch of sharks in WA and IUU

#notes: 'wetline' (as reported in Heupel & McAuley 2007) groups the byproduct of several 
#       different fisheries. This is already included in CAESS

#MISSING: Mervi to send species numbers from Appendix 3.1 and 3.2
#         See Table 6.3 Heupel & McAuley 2007. They report annual shark catch by fishery in the north!!
#           are all these fisheries captured in 'other fisheries' from returns????


# Annual updates for fisheries listed in: #Total landings time series
# (each year download annual reported landings from FISHCUBE using the FishCube fishery code in 'Lista.reap.FishCubeCode')

#Asses.year=2019     #delete once script is done as this is declared in "Assessment.R" & "Assessment_other.species.R"

library(tidyverse)
library(readxl)
library(lubridate)


Do.recons.paper="NO"  #outputs for reconstruction paper
Historic.yrs=1941:1975

#Catch reconstruction scenarios for ...what???            MISSING
ScenarioS=data.frame(Scenario=c('Base Case','High','Low'),
                     PCM=c(1,1.5,.5),
                     Weight=c(1,1.5,.5))


# 1 -------------------DATA SECTION------------------------------------

options(stringsAsFactors = FALSE)

#Fishery codes
FisheryCodes=read_excel('C:/Matias/Data/Catch and Effort/FisheryCodeTable.xlsx', sheet = "CAEStoFISHCUBE")


#Total landings time series                          #MISSING, have all fisheries here....
Whaler_SA=read.csv("C:/Matias/Data/Catch and Effort/SA_marine_scalefish_whaler_ktch.csv")
WTBF_catch=read.csv("C:/Matias/Data/Catch and Effort/WTBF_catch_Benseley.et.al.2010.csv")  #catch in kg
GAB.trawl_catch=read.csv("C:/Matias/Data/Catch and Effort/Gab.trawl_catch_Benseley.et.al.2010.csv")  #catch in kg


#Shark bio data
User="Matias"
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R")

#Species codes
All.species.names=read.csv("C:/Matias/Analyses/Population dynamics/1.Other species/Species_names.csv")
species.codes=read.csv("C:\\Matias\\Data\\Species.code.csv")   #previously declared as 'b'


#1.1 Catch_WA Fisheries
setwd("C:/Matias/Analyses/Data_outs")

  #-- 1.1.1 all monthly return and daily logbook fisheries from 1975-76

    #south of 26 S
Data.monthly=read.csv("Data.monthly.csv")
#note:  this has catch of sharks from monthly returns of all fisheries 
#       but from daily records of shark fisheries only 

daily.other=read.csv("Data.daily.other.fisheries.csv")                
#note:  this has shark catch from daily records from other fisheries   
  
    #north of 26 S
Data.monthly.north=read.csv("Data.monthly.NSF.csv")
Data.monthly.north$LIVEWT.c=Data.monthly.north$LIVEWT  


  #-- 1.1.2 Historic shark fisheries (pre 1975)

      #1.1.2.1 Whitley 1944 (citation in Simpfendorfer & Donohue 1998)
#note: All shark catch across southern WA (south of Abrolhos islands). Weight in pounds
Catch_1941_1943=data.frame(year=1941:1943,LIVEWT.c=c(55332,77441,109064)) 
Catch_1941_1943$LIVEWT.c=Catch_1941_1943$LIVEWT.c*0.453592/1000  #Convert to tonnes

      #1.1.2.2 Heald 1987
#note: All shark catch by port. Weight in kgs
Catch_1949=read.csv("C:/Matias/Data/Catch and Effort/Heald(1987)_1949_total_catch_by port_allspecies.csv") 
Catch_1949$LIVEWT.c=Catch_1949$Shark.catch.kg_live_wgt./1000   #in tonnes

      #1.1.2.3 Simpfendorfer & Donohue 1998
#note: All shark catch for southern WA (Geraldton to SA border). Weight in tonnes
Catch_1950=data.frame(year=1950,LIVEWT.c=50)   

      #1.1.2.4 Heald 1987
#note= All shark catch for southern WA (Geraldton to SA border). Weight in tonnes
Catch_1952_1975=read.csv("C:/Matias/Data/Catch and Effort/Historical_WA_shark_catches_1952_1975.csv") 
names(Catch_1952_1975)=c("year","LIVEWT.c")


  #-- 1.1.3 Wetline  
  #Western Rock lobster
WRL=read.csv("C:/Matias/Data/Catch and Effort/WRL/Number-of-vessels.csv")
WRL.Wann=read.csv("C:/Matias/Data/Catch and Effort/WRL/Wann_catch.csv")


  #-- 1.1.4 TEPs   
    #1.1.4.1 TDGDLF
#note: data obtained from Comments in TDGDLF returns

# Oversized dusky sharks                            
Overzd.Whler=read.csv("C:/Matias/Data/Catch and Effort/TEPS/Comments from SharkDailyTripReturns.csv")

    #1.1.4.2 Other fisheries
TEPS_other=read.csv("C:/Matias/Data/Catch and Effort/TEPS_catch.csv")

    #1.1.4.3 Greynurse catches since protection
GN.mn.wght=25
GN.pcm=.2
Under.rep.factor=2  #from white shark catch recons
grey.hndl=paste('C:/Matias/Analyses/Catch and effort/State of fisheries/',
                paste(Asses.year-2,substr(Asses.year-1,3,4),sep='-'),sep='')
Greynurse.ktch=read.csv(paste(grey.hndl,'3.Table2.TEPS.csv',sep='/'))


#-- 1.1.5 Pilbara trawl shark species composition (go to scientist: Corey W.) 
#source: Table 6.4 Heupel & McAuley 2007 (weight is live weight in kg)
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
  Retained=c(rep('yes',10),rep('no',13)),
  Numbers=c(90,7,3,27,80,5,238,559,114,102,101,4,115,25,1,12,43,1,4,6,1,1,2),
  Weight=c(2058.6,595.9,295.1,395.8,593.2,168.3,522.6,772.7,210.5,214.2,
           1627.8,257.1,189.5,90.7,39.4,20.0,15.6,15.1,9.5,4.9,4.2,2,.8)
)
#Pilbara.trawl.observed.effort=100  #days (Stephenson & Chidlow 2003 fide in McAuley et al 2005 page 25)
Pilbara.trawl.observed.effort=1601.898   #hours; derived from get.observed.Pilbara
get.observed.Pilbara=FALSE
if(get.observed.Pilbara)
{
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
#BRD_pilbara.trawl_prop.shark=xx  #source ???       MISSING!!!!!!!!!
#BRD_pilbara.trawl_prop.ray=xx  
#BRD_pilbara.trawl_year='200x-0x'



  #-- 1.1.6 WA scallop and prawn trawl fisheries (go to scientist: Mervi Kangas)

    #source Laurenson et al 1993 Table 5
# notes: Authors grouped all sites and replicates in Table 5 
#        Authors reported a range of number of individuals so the mid point is used
South.west.trawl.observed.comp=data.frame(
    Scien.nm=c('Squatina australis','Aulohalaelurus labiosus','Urolophus circularis','Orectolobus tentaculatus',
             'Pristiophorus cirratus','Myliobatis australis','Trygonorhina fasciata','Mustelus antarcticus',
             'Urolophus lobatus','Trygonoptera personata','Hypnos monopterygium','Heterodontus portusjacksoni',
             'Sphyrna zygaena','Dasyatis brevicaudata','Aptychotrema vincentiana','Urolophus paucimaculatus',
             'Parascyllium variolatum','Urolophus mucosus','Orectolobus sp.'),
  
    Species=c('Angel shark','Black-spotted catshark','Circular stingaree','Cobbler carpetshark',
            'Common sawshark','Eagle ray','Fiddler ray','Gummy shark',
            'Lobed stingaree','Masked stingaree','Numbfish','Port Jackson shark',
            'Smooth hammerhead','Smooth stingray','Southern shovelnose ray','Sparsely-spotted stingaree',
            'Varied catshark','Western stingaree','Western wobbegong'),
    Commercial=c('y',rep('n',3),rep('y',4),rep('n',4),'y','n','y',rep('n',3),'y'),
    Bell.buoy=c(5,0,0,0,0,5,5,0,75,30,0,5,0,0,5,30,0,30,0),
    Cottesloe=c(5,0,0,0,0,30,5,0,75,30,0,550,0,5,5,5,0,30,0),
    Conventry.reef=c(5,0,0,0,0,5,5,5,550,30,0,5,0,5,5,75,0,5,0),
    Comet.bay=c(30,0,0,0,5,5,5,5,5,5,5,30,5,0,5,5,0,30,0),
    Preston.deep=c(5,5,5,0,5,5,5,0,75,30,5,5,0,0,5,5,0,30,0),
    Preston.shallow=c(0,5,0,0,0,5,5,0,0,5,5,5,0,0,0,0,5,5,5),
    Capel=c(5,0,0,0,0,5,5,0,30,30,0,5,0,0,5,550,0,30,0),
    Busselton=c(5,0,0,5,5,5,5,0,0,5,0,5,0,0,0,75,0,5,0),
    zoneC=c(30,0,0,0,5,5,5,0,0,5,0,30,0,0,5,75,0,75,0))
South.west.trawl.observed.effort=9*4*(15/60)*4  #hours, 9 sites X 4 replicates X 15 mins X 4 times in a year (Laurenson et al 1993 page 13)
South.west.trawl.observed.PCM= 0.5 #observed for gummy sharks after 7 days of capture


    #Source Kangas et al 2006 Appendix 3.1 and 3.2
# notes: Authors grouped all sites and replicates in Appendix 3.1 and 3.2 and only report presence/absence
Shark.Bay.trawl.observed.comp=data.frame(
  Scien.nm=c('Chiloscyllium punctatum','Halaelurus boesemani','Mustelus sp. A',
             'Carcharhinus cautus','Carcharhinus melanopterus','Carcharhinus obscurus',
             'Carcharhinus plumbeus','Rhizoprionodon acutus Shark','Hemigaleus australiensis',
             'Rhina ancylostoma','Rhynchobatus australiae',
             'Aptychotrema vincentiana','Narcine westraliensis','Hypnos monopterygium',
             'Dasyatis kuhlii','Dasyatis leylandi','Himantura sp.','Himantura toshi',
             'Himantura uarnak','Taeniura meyeni','Gymnura australis','Trygonoptera ovalis',
             'Aetobatus narinari'),
  Species=c('Catshark, Brown-banded','Catshark, Speckled','Shark, Grey Gummy',
            'Shark, Nervous','Shark, Blacktip Reef','Shark, Dusky Whaler',
            'Shark, Sandbar','Shark, Milk','Shark, Weasel','Shark Ray','Shovelnose Ray, White-spotted',
            'Shovelnose Ray, Western','Numbfish, Banded','Numbfish',
            'Stingray, Blue-spotted','Stingray, Brown reticulated','Stingray, Coachwhip','Whipray, Black-spotted',
            'Whipray, Reticulate','Stingray, Black-blotched','Ray, Rat-tailed/Butterfly','Stingaree, Striped',
            'Ray, White-spotted Eagle'),
  Commercial=c('n','n',rep('y',7),'n','n','y',rep('n',10),'y'),
  Total=c(2,1,3,1,1,1,1,1,1,1,3,3,4,1,3,4,1,2,2,1,4,2,1))
Shark.Bay.trawl.observed.effort=26*3*(10/60)*4  #hours, 26 sites X 3 replicates X 10 mins X 4 times in a year (Kangas et al 2006 Table 2.1)


Exmouth.Onslow.trawl.observed.comp=data.frame(
  Scien.nm=c('Chiloscyllium punctatum','Eucrossorhinus dasypogon','Atelomycterus sp.','Rhizoprionodon acutus',
             'Hemigaleus australiensis','Rhynchobatus australiae','Rhinobatos typus','Dasyatis kuhlii',
             'Dasyatis leylandi','Himantura toshi','Gymnura australis',
             'Aetomylaeus vespertilio','Aetomylaeus nichofii'),
  Species=c('Catshark, Brown-banded','Wobbegong, Tasselled','Catshark, Banded','Shark, Milk',
            'Shark, Sicklefin Weasel','Ray, White-spotted Shovelnose','Ray, Giant Shovelnose','Ray, Blue-spotted Stingray',
            'Ray, Brown Reticulated','Ray, Black-spotted Whipray','Ray, Butterfly/Rat-tailed',
            'Ray, Ornate Eagle Ray','Ray, Banded Eagle Ray'),
  Commercial=c(rep('n',3),'y','y',rep('n',8)),
  Total=c(3,1,3,2,2,2,1,1,3,3,3,1,1))
Exmouth.Onslow.trawl.observed.effort=26*3*(10/60)*3  #hours, 26 sites X 3 replicates X 10 mins X 3 times in a year (Kangas et al 2006 Table 2.1)


BRD_prawn.trawl_prop.shark=9/70  #Kangas & Thomson 2004 :70 sharks retained with no BRD VS 9 with BRD
BRD_prawn.trawl_prop.ray=8/65  #                         65 rays retained with no BRD VS 8 with BRD
BRD_prawn.trawl_year='2003-04'


  #-- 1.1.7 Kimberley Gillnet and Barramundi Fishery and Eighty Mile Beach Gillnet Fishery

  #source: Table 6.5 Heupel & McAuley 2007 (weights are live weight in kg)
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
  Numbers=c(87,32,36,3,1,45,27,3,13,82,190,31,398,1,17,1,6,22,1,5,4,1,5),
  Weight=c(578.1,126.5,178.6,12.5,109.1,106.9,165.1,3.4,13.9,422.5,5239.1,370.1,
           2284.5,23.7,161.6,1.4,5.4,63.5,5.7,20.0,2.9,57.2,56.8)
)
Kimberley.GBF.observed.effort=160  #days (McAuley et al 2005 page 26; 5 vessels observed)


#Effort time series (must be updated every year)       MISSING, add all fisheries for which I use annual effort or catch
fn.rid.efrt=function(d) paste('C:/Matias/Data/Catch and Effort/Effort_other_fisheries',d,sep='/')

  #KGBF
Kimberley.GBF.annual.effort=read.csv(fn.rid.efrt('KGBF Annual catch and Bdays.csv'))    #days (KGBF and 80 mile beach combined)

  #Pilbara trawl
#Pilbara.trawl.annual.effort=read.csv()   #hours

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


#1.2. Catch of non WA Fisheries

  #-- 1.2.1 Taiwanese gillnet (1974-1986) and longline (1989-91) fishery (catch in tonnes)
#Sources: Stevens & Davenport 1991; Stevens 1999
Taiwan.gillnet.ktch=data.frame(
  Year=1974:1986,
  Total.Ktch=c(321.4,8997.6,6455.3,9970.7,5528.2,3282.1,5831.1,6694.7,
               5624.1,7589.9,6544.2,2929.5,2111.1),
  Shark.Ktch=c(rep(NA,5),557.9,4234.5,4486.2,3639.4,4418.6,2731.4,2327.4,2393.9),
  WA.Ktch=c(rep(NA,5),50,1250,720,800,790,400,10,70))  #From Figure 7 Stevens & Davenport

Taiwan.gillnet.sp.comp=as.data.frame(matrix(c(.357,.166,.059),ncol=3))
colnames(Taiwan.gillnet.sp.comp)=c("Australian.blacktip.Shark.West.prop",
                                   "Spot.tail.Shark.West.prop","Hammerheads.West.prop")

Taiwan.longline.ktch=data.frame(Year=1990:1991,
                                Shark.Ktch=c(1700*11/(11+9),1700*9/(11+9)))
Taiwan.longline.sp.comp=data.frame(
  Species=rep(c("Spot-tail shark","Australian blacktip shark","Tiger shark",
                "Milk shark","Spinner shark",
                "Pigeye shark","Graceful shark","Hammerheads","Other"),2),
  Year=c(rep(1990,9),rep(1991,9)), 
  Percent=c(18.3,25.6,8.1,.4,6.9,17.9,1.2,2,19.5,5.6,79.7,0.7,2,4.5,1,0,0.1,6.4))



  #-- 1.2.2 Commonwealth GAB trawl and Western Tuna and Billfish Fisheries (WTBF)     ACA     MISSING!!!!!!!
#source: Bensely et al 2010. Table 2 in Appendix A
WTBF_effort=read.csv("C:/Matias/Data/Catch and Effort/WTBF_effort_Benseley.et.al.2010.csv")  #hook numbers
GAB.trawl_effort=read.csv("C:/Matias/Data/Catch and Effort/GAB.trawl_effort_Benseley.et.al.2010.csv")  #hours trawled


#Stobutski et al 2006 Table 3 
WTBF_observed=read.csv("C:/Matias/Data/Catch and Effort/WTBF_Stobutzki.et.al.2006.csv")  #catch in numbers

cpue_dusky_WTBF=37/203205  #number of individuals per 203205 hooks observed.
cpue_sandbar_WTBF=8/203205   
dusky_WTBF.at.vessel.mortality=1-.97  # 97% discarded alive
sandbar_WTBF.at.vessel.mortality=1-1  # 100% discarded alive


  #-- 1.2.3 SA Marine Scalefish fishery (source Taylor et al 2015)
#description: whaler shark catch from SA MArine Scale fishery (in tonnes). 
Whaler_SA_dusky.prop=.1  # Steer et al 2018 (page 148)
Whaler_SA_bronzie.prop=1-Whaler_SA_dusky.prop


#--Length weigth
bwt=3.47e-06
awt=3.10038

#-- Add fisheries code and name
FisheryCodes=FisheryCodes%>%data.frame
Data.monthly=Data.monthly%>%left_join(FisheryCodes,by=c("FisheryCode"="SASCode"))
Data.monthly.north=Data.monthly.north%>%left_join(FisheryCodes,by=c("FisheryCode"="SASCode"))



#-- Select fisheries for reapportioning 'shark other' & reconstructing discards post 2006 if appropriate
#notes: Total of 51 fisheries have reported Shark/ray catches in WA
#       Criteria for calculating discards post 2006: > 10 tonnes (all years with reported catches combined)
#       Other fisheries reported <10 tonnes between 1975 and 2006 so discards post 2006 no calculated as 
#         catches are negligible but the reported catches from these fisheries are considered 
Calculate.discarding=c('PFT','C019','C066','C070','CSLP','EGBS','EGP','KGB','KP','KTR','NBP',
                       'OANCGCWC','OP','SBP','SBS','SBSC','SCT','SWT')
Reaportion.from.reported.ktch=c('C019','C066','C070','CSLP','EGBS','JANS','JASDGDL','KGB','KTR','NCS',
                                'OANCGCWC','OASC','SBS','WANCS','WCDGDL')
Reaportion.from.survey=c('PFT','EGP','KP','NBP','OP','SBP','SBSC','SCT','SWT')

Lista.reap.FishCubeCode=list(Pilbara.trawl='PFT',
                             Estuaries.19='C019',
                             Cockburn.Sound.fish.net='C066',
                             Power.hauler.70='C070',
                             Cockburn.Sound.line.pot='CSLP',
                             Exmouth.beach.seine.mesh='EGBS',
                             JANSF=c('JANS','NCS'),
                             JASDGDL='JASDGDL',
                             Kimberley.gillnet.barra='KGB',
                             Kimberley.prawn='KP',
                             Kimberley.trap='KTR',
                             Nickol.Bay.prawn='NBP',
                             Open.north.gas.west='OANCGCWC',
                             Open.south=' OASC',
                             Onslow.prawn='OP',
                             Shark.Bay.prawn='SBP',
                             Shark.Bay.snapper='SBS',
                             Shark.Bay.scallop='SBSC',
                             South.coast.trawl='SCT',
                             South.west.trawl='SWT',
                             WANCS='WANCS',
                             WCDGDL='WCDGDL')
  

 some.ktch.in.daily.other=c('OANCGCWC','OASC','PFT')  #note, for these, only calculate discards for years not reported in daily.other
 
 





#DEJE ACA 
# Careful when calculating total catch as 'Organise.data.R' sources this script and uses
# 'Data.monthly', etc to calculate catch by fishery and total catch. Make sure other
#  scripts like 'Assessment.other.R' source this script and Catch.recons.Recreational.R'

#  Careful if having rm(all) statements 

# Move code from 'Assessment_other.species.R' to here. Also include 'daily.other'.
#       for reapportioning, consider fisheryCode and area.....

# Oversized dusky sharks : aDD 2011-12 AND 2012-13 YEARS!!!

# Series starts in 1941 so add 0 until start of each fishery

# Prior November 2006: 
#       all fisheries should be in Data.monthly, Data.monthly.north & daily.other
#       hence, all I need to do is to reapportion 'shark,other' using the
#       observed species composition, and if not available, the composition of
#       the reported catch that is reported in CAESS (as done in 'Assessment_other.species.R')

# Post November 2006: 
#       All sharks from non-shark fisheries & Pilbara trawl, KGB and ??, are not
#       allowed to be retained, hence this discarded catch has to be estimated for
#       all fisheries that reported shark catches prior to November 2006
#       For this, if catch rate information is available, do catch rate X annual effort X PCM
#       Also factor in BRD in case of trawl fisheries.
#       If catch rates are not in weight, then need to convert numbers to weights
#       If catch rate data is not available, then do Total catch X proportion in catch reported prior 2006 X PCM

# Do TEPS



# 2 -------------------PROCEDURE SECTION------------------------------------

#2.1 Catch from WA Fisheries

  #2.1.1 combine historic in single dataframe
Catch_1949.total=data.frame(year=1949,LIVEWT.c=sum(Catch_1949$LIVEWT.c))     
Catch_1949.total$LIVEWT.c=Catch_1950$LIVEWT.c     
Historic.ktch=rbind(Catch_1941_1943,Catch_1949.total,Catch_1950,Catch_1952_1975)

#add missing years through linear interpolation
Missing.yrs=Historic.yrs[which(!Historic.yrs%in%Historic.ktch$year)]
Missing.ktch=approx(Historic.ktch$year,Historic.ktch$LIVEWT.c,xout=Missing.yrs)   

fn.fig("C:/Matias/Analyses/Catch and Effort/Pre_1975_total_shark_catch",2400, 2400) 
par(mfcol=c(1,1),las=1,mgp=c(2.9,.6,0))
with(subset(Historic.ktch,year<1975),
     {
       plot(year,LIVEWT.c, ylab="Total shark landings (tonnes)",xlab="Financial year",
            pch=19,col=1,cex=2,cex.lab=1.7,cex.axis=1.35)
     })
dev.off()

Historic.ktch=rbind(Historic.ktch,data.frame(year=Missing.ktch$x,LIVEWT.c=Missing.ktch$y))
Historic.ktch=Historic.ktch[order(Historic.ktch$year),]
Historic.ktch$LIVEWT.c=Historic.ktch$LIVEWT.c*1000   #convert back to kg

#get prop of catch for 1975-1980 to extract historic
ALL=subset(Data.monthly,SPECIES<32000 & METHOD%in%c("GN","LL") & LAT<=(-26) & Estuary=="NO"
           & FINYEAR%in%c("1975-76","1976-77","1977-78","1978-79","1979-80","1980-81"),
           select=c(FINYEAR,LIVEWT.c,SPECIES,SNAME))%>%
  filter(!SPECIES==22999)
Prop.sp.yr.all=ALL%>%group_by(SPECIES)%>%
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
  mutate(LIVEWT.c=ktch*Proportion)


#get data for other gears from Data.monthly
TDGDLF.other.gears=Data.monthly%>%
  filter(!METHOD%in%c("GN","LL"))%>%
  mutate(dummy=paste(FINYEAR,MONTH,VESSEL,BLOCKX,METHOD))
NSF.other.gears=Data.monthly.north%>%
  filter(!METHOD%in%c("GN","LL"))%>%
  mutate(dummy=paste(FINYEAR,MONTH,VESSEL,BLOCKX,METHOD))


  #2.1.2 reapportioning of 'shark,other' and 'hammerheads'
fn.subs=function(YEAR) substr(YEAR,start=3,stop=4)


  #2.1.3 Wetline catches of dusky sharks                   
  #rock lobster
#note: extrapolate Wann's catch (for one season) to the entire fishery as a proportion of the number of boats
WRL.Wann=WRL.Wann%>%mutate(TL=TL_metres*100,
                           TL=ifelse(is.na(TL),TL_feet*30.48,TL))   #convert to cm
fn.weight=function(TL,bwt,awt) bwt*TL^awt
#1.  Assumption on proportion of lobster boats setting droplines
WRL.prop=0.1  #small number of operators used the gear (Taylor et al 2015)     #is this already in CAESS??
Dusky.WRL=WRL.Wann%>%
            filter(Species=="DW")%>%
            mutate(LiveWt=fn.weight(TL,bwt,awt))   #in kg
  
if(KTCH.UNITS=="TONNES") Annual.Dusky.Ktch.WRL.Wann=sum(Dusky.WRL$LiveWt)/1000 else #in tonnes
                         Annual.Dusky.Ktch.WRL.Wann=sum(Dusky.WRL$LiveWt)
Dusky.Tot.Ktch.WRL=WRL
Dusky.Tot.Ktch.WRL$LiveWt=WRL.prop*Annual.Dusky.Ktch.WRL.Wann*Dusky.Tot.Ktch.WRL$Number.of.vessels  #in tonnes
colnames(Dusky.Tot.Ktch.WRL)[2:3]=c("FINYEAR","LIVEWT.c")
Dusky.Tot.Ktch.WRL=Dusky.Tot.Ktch.WRL[,-1]   #remove number of vessels     


#2.2. Catch from non WA Fisheries

#Taiwanese gillnet and longline
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


  #2.2.1 Commonwealth SWTBF and GAB trawl fishery               
#Option=1       #select source to use to obtain catches
Option=2
if(KTCH.UNITS=="TONNES")   #convert to tonnes
{
  WTBF_catch$Livewt=WTBF_catch$Catch_kg/1000
  GAB.trawl_catch$Livewt=GAB.trawl_catch$Catch_kg/1000
}
    

if(SP=="BW") 
{
  GAB_trawl=AFMA_catch[,c(1,2)]     #MISSING: I changed AFMA_catch for WTBF_catch GAB.trawl_catch  ditto Effort_WTBF
  Bronzey_WTBF=AFMA_catch[,c(1,5)]
  names(Bronzey_WTBF)[2]='LIVEWT.c'
}

if(SP=="GM") GAB_trawl=AFMA_catch[,c(1,3)]
if(SP=="WH") GAB_trawl=AFMA_catch[,c(1,4)]

if(!SP=="TK") names(GAB_trawl)[2]='LIVEWT.c'

if(SP%in%c("BW","TK"))
{
  #option 1 (data reported by AFMA)
  if(Option==1)    
  {
    Dusky_WTBF=AFMA_catch[,c(1,6)]  
    Sandbar_WTBF=AFMA_catch[,c(1,7)]
  }
  
  #option 2 (weigh up observed catch rate by effort in Commonwealth waters, an assumed 
  # mean weight and an assumed PCM
  Effort_WTBF=subset(Effort_WTBF,Fleet=="Australian")
  if(Option==2)    
  {
    PCM=0.3
    Mean.wt.dusky=50   #MISSING: get mean weight from size frequency composion of longlines   
    Mean.wt.sandbar=20
    Dusky_WTBF=Effort_WTBF[,1:2]
    #Dusky_WTBF=AFMA_catch[,c(1,6)]  
    Dusky_WTBF[,2]=Effort_WTBF[,2]*PCM*cpue_dusky_WTBF*Mean.wt.dusky    #MISSING, replaced cpue_dusky_WTBF with WTBF_observed
    if(KTCH.UNITS=="TONNES") Dusky_WTBF[,2]=Dusky_WTBF[,2]/1000
    Sandbar_WTBF=Effort_WTBF[,1:2]
    #Sandbar_WTBF=AFMA_catch[,c(1,7)]
    Sandbar_WTBF[,2]=Effort_WTBF[,2]*PCM*cpue_sandbar_WTBF*Mean.wt.sandbar
    if(KTCH.UNITS=="TONNES")Sandbar_WTBF[,2]=Sandbar_WTBF[,2]/1000
  }
  names(Dusky_WTBF)[2]=names(Sandbar_WTBF)[2]='LIVEWT.c'
}

#add Bronzey catch to dusky catch for trawl and WTBF?????
if(SP=="BW")
{
  ADD.bronzy_dusky="NO"        
  if(ADD.bronzy_dusky=="YES")
  {
    Dusky_GAB_trawl=Bronzey_GAB_trawl
    Dusky_WTBF[,2]=Dusky_WTBF[,2]+Bronzey_WTBF[,2] 
    names(GAB_trawl)[2]="LIVEWT.c"
    GAB_trawl=fn.exp.his(GAB_trawl,Monthly.prop)
  }
}

if(SP=="BW") Dusky_WTBF$FINYEAR=as.numeric(substr(Dusky_WTBF$FINYEAR,1,4))
# partition by month

if(SP=="BW")
{
  GAB_trawl=fn.exp.his(GAB_trawl,Monthly.prop) 
  colnames(Dusky_WTBF)[match("Calendar.Year",colnames(Dusky_WTBF))]="FINYEAR"
  Dusky_WTBF=fn.exp.his(Dusky_WTBF,Monthly.prop)
}

if(SP=="TK")
{
  colnames(Sandbar_WTBF)[match("Calendar.Year",colnames(Sandbar_WTBF))]="FINYEAR"
  Sandbar_WTBF=fn.exp.his(Sandbar_WTBF,Monthly.prop)
}


  #2.3  TEPS
#TDGDLF
if(SP=="BW")
{
  Overzd.Whler$FINYEAR=as.character(with(Overzd.Whler,
                                         ifelse(Month%in%1:6,paste(Year-1,"-",fn.subs(Year),sep=""),
                                                ifelse(Month%in%7:12,paste(Year,"-",fn.subs(Year+1),sep=""),NA))))
  
  # Also note that for most whalers, (e.g. Tigers), comments don't say if oversized or not. So 
  #assumed they are
  Lista.protected.elasmos=list( Dusky=c("dusky","Bronze","bronze","b/w","Dusky",
                                        "whalers","BW","Bronzy","b/whaler"))
  
  
  Lista.id.pro.el=Check.These.Sks=vector('list',length=length(Lista.protected.elasmos))
  names(Lista.id.pro.el)=names(Check.These.Sks)=names(Lista.protected.elasmos)
  
  fn.find=function(what,inwhat)
  {
    Find=regexpr(what, inwhat$Comments) > 0
    return(which(Find==T))
  }
  
  for (i in 1)
  {
    these=Lista.protected.elasmos[[i]]
    store=NULL
    if(length(these)>1)for(j in 1:length(these))store=c(store,fn.find(these[j],Overzd.Whler))    
    if(length(these)==1)store=fn.find(these,Overzd.Whler)
    Lista.id.pro.el[[i]]=store
  }
  
  for (i in 1:length(Lista.id.pro.el))
  {
    id=Lista.id.pro.el[[i]]
    a=Overzd.Whler[id,]
    a=a[order(a$FINYEAR),]
    Check.These.Sks[[i]]=a
  }
  
  #The only manual bit is to go thru each element of Check.These.Sks and compare comments with data
  i=1;print(Check.These.Sks[[i]][,6:7])
  
  TEPS_TDGDLF=data.frame(FINYEAR=paste(2006:2012,"-",fn.subs(2007:2013),sep=""))
  TEPS_TDGDLF$N.dusky.alive=c(5,4,77,8,16,11,2)
  TEPS_TDGDLF$N.dusky.dead=c(1,1,22,7,9,32,13)
  TEPS_TDGDLF$N.dusky.U=c(1,5,8,10,2,0,0)
  
  #Catch is number dead X average weight of 3 m sharks + PCM * (number alive + number unknown)
  Mn.wt=bwt*300^awt
  if(KTCH.UNITS=="TONNES")Mn.wt=Mn.wt/1000
  TEPS_TDGDLF$LIVEWT.c=(TEPS_TDGDLF$N.dusky.dead*Mn.wt)+ ((TEPS_TDGDLF$N.dusky.alive+TEPS_TDGDLF$N.dusky.U)*Mn.wt*PCM)       
  
  #1.8.2 Other fisheries
  these.sp=c("Shark, Bronze Whaler","Shark, Bronze Whaler (oversize)","Shark, Dusky Whaler (oversize)")
  #note: no sandbars reported
  
  TEPS_other=subset(TEPS_other,CommonName%in%these.sp & !FisheryCode%in%c('SGL','TDGDLF','WCGL')) #avoid duplication

  TEPS_dusky=TEPS_TDGDLF[,c(1,5)]
  TEPS_dusky=fn.exp.his(TEPS_dusky,Monthly.prop)
  
}

# Greynurse catches since protection      REVIEW in light of TEP stuff for dusky....
Greynurse.ktch=Greynurse.ktch%>%
  filter(Species=="SHARK, GREY NURSE")
Grey.nms=colnames(Greynurse.ktch)[-1]
Al.ded=c(NA,sapply(strsplit(Grey.nms, "[.]"), "[", 1))
Greynurse.ktch=data.frame(
  finyear=unique(as.numeric(substr(sapply(strsplit(Grey.nms, "[.]"), "[", 2),1,4))),
  N.alive=unlist(Greynurse.ktch[,which(Al.ded=="alive")]),
  N.dead=unlist(Greynurse.ktch[,which(Al.ded=="dead")]))%>%
  mutate(LIVEWT.c=Under.rep.factor*(N.dead*GN.mn.wght+N.alive*GN.mn.wght*GN.pcm))%>%
  dplyr::select(finyear,LIVEWT.c)

# 2.3 WA scallop and prawn trawl fisheries
South.west.trawl.observed.comp=South.west.trawl.observed.comp%>%
       mutate(Total = rowSums(select(., -Scien.nm,-Species,-Commercial)))%>%
      dplyr::select(Species,Scien.nm,Commercial,Total)

Shark.Bay.trawl.observed.comp

Exmouth.Onslow.trawl.observed.comp


#SA Marine Scalefish fishery        #MISSING: add these two to total catch!!!
#split dusky from bronzy                      
names(Whaler_SA)[2]="LIVEWT.c"
Whaler_SA$Dusky=Whaler_SA$LIVEWT.c*Whaler_SA_dusky.prop
Whaler_SA$Bronzie=Whaler_SA$LIVEWT.c*Whaler_SA_bronzie.prop      



# 3 -------------------REPORT SECTION------------------------------------
if(Do.recons.paper=="YES")   #for paper, report only IUU and reconstructions (not the TDGDLF or NSF)
{
  
}
