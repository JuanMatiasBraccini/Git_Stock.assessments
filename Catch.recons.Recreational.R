#-- Script for reconstructing time series of recreational catch of sharks in WA

#missing: consider "#Perth metro pilot survey (Claire Smallwood)"

#         MAKE sure all species (rec, comm, discards) are considered (e.g PJs)     VIP!!!

# Annual updates:
        # For each new I-survey, get data from Karina
        # Update Charter boat data each year
        # Update WA.population


#notes: This script uses I-Survey (boat-based) point estimates to reconstruct
#       time series considering rec fishing participation rates and WA population size
#       Shore-based fishing is added as a proportion of boat-based fishing
#       Charter boat fishing is also added


#Notes on Charter boat (from Rhonda):
#       There is no reliable charter data before this date as it was not compulsory to send in logbook sheets.
#         Compulsory logbooks were introduced in 2001.
#         We have  a little data on the sharks for 1998 and 2000
#       Bronze Whaler have been coded to Whalers general from 2014. So we have no specific lengths for Bronzies after this date.
#       Rory did not trust that the sharks were being identified correctly so they were all put under whaler general code.

#      To address this issue reconstructed catch time series are
#        calculated using only the Isurvey years (combining Isurvey, beach fishing and charter
#        boats which should be reliable) and the back calculating using WA population size
#        and participating rate


#Notes on Isurvey (from Karina):
# The data is an extract of what will be the most up-to-date iSurvey data to be made
#   available on fishcube in the near future, however, there are final checks to complete
#   before release – will let you know when final estimates are available# 
# The variable names all appear to be correct for your coding; the variable names
#   may change slightly as there may be some differences between the extract provided and fishcube 
# The variables Kept, Released and Total have always been estimates with decimals, 
#   but rounded to integers for publication / release – and should have been rounded in the extract 
options(dplyr.summarise.inform = FALSE)

library(tidyverse)
library(readxl)
library(lubridate)
library(Hmisc)

Do.recons.rec.fishn.paper="NO"


#Catch reconstruction scenarios
Scenarios=data.frame(Scenario=c('Base Case','High','Low'),
                     PCM=c(1,1.5,.5),
                     Weight=c(1,1.5,.5))

# 1 -------------------DATA SECTION------------------------------------

# I-Survey
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')
Rec.hndl=handl_OneDrive("Data/Catch and Effort/Recreational/I.Survey.")
#Rec.fish.catch.2011.12=read.csv(paste(Rec.hndl,"2011_12.csv",sep=''),stringsAsFactors=F) #Ryan et al 2013    
#Rec.fish.catch.2013.14=read.csv(paste(Rec.hndl,"2013_14.csv",sep=''),stringsAsFactors=F) #Ryan et al 2015 
#Rec.fish.catch.2015.16=read.csv(paste(Rec.hndl,"2015_16.csv",sep=''),stringsAsFactors=F) #Ryan et al 2017 
#Rec.fish.catch=rbind(Rec.fish.catch.2011.12,Rec.fish.catch.2013.14,Rec.fish.catch.2015.16)
Rec.fish.catch=read.csv(paste(Rec.hndl,"csv",sep=''),stringsAsFactors=F)
Scien.nm=Rec.fish.catch%>%
            rename(Common.Name=Lowlevelgrouping,
                   Scientific.name=ScientificName)%>%
            dplyr::select(Common.Name,Scientific.name)%>%
            distinct(Common.Name,.keep_all=T)

Rec.fish.catch=Rec.fish.catch%>%
            mutate(FinYear=paste(2000+as.numeric(substr(Year,1,2)),substr(Year,3,4),sep="-"),
                   Bioregion=paste(sapply(strsplit(Rec.fish.catch$RegionReportingGroupName, "_"), "[", 2),
                                   sapply(strsplit(Rec.fish.catch$RegionReportingGroupName, "_"), "[", 3)),
                   Common.Name=Lowlevelgrouping,
                   Scientific.Name=ScientificName,
                   Kept.Number=ceiling(Kept),
                   Kept.Number.se=se.Kept,
                   Rel.Number=ceiling(Released),
                   Rel.Number.se=se.Released,
                   RSE=rse.Total,
                   Sample.size=counts)%>%
      dplyr::select(FinYear,Bioregion,Sample.size,Common.Name,Scientific.Name,
                    Kept.Number,Kept.Number.se,Rel.Number,Rel.Number.se,RSE)
I.survey.years=unique(Rec.fish.catch$FinYear)

# Shore-based
#notes from Karina:
  #statewide estimated kept, released & total recreational catch for sharks & rays in 2000_01
  #K_SE = standard error associated with Kept
  # K_RSE = relative standard error associated with Kept
  # K_HHs = number of households that reported a Kept catch 
  # R = released, T = Total 
Shore.based=read.csv(handl_OneDrive("Data/Catch and Effort/Recreational/statewide shark 2000_01.csv"),stringsAsFactors=F)

#Perth metro pilot survey (Claire Smallwood)
Shore.based.metro.pilot=read.csv('handl_OneDrive(Data/Catch and Effort/Recreational/RawCatch_Sharks_Export.csv'),stringsAsFactors=F)


# Charter boats
Charter=read_excel("handl_OneDrive(Data\\Catch and Effort\\Charter\\Charter.xlsx"),sheet ='Data')


# WA population for rec catch recons (ABS) (Google "What is the population of Western Australia 20xx?")
#source: https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/3101.0Dec%202018?OpenDocument
WA.population=read.csv(handl_OneDrive("Data/AusBureauStatistics.csv"),stringsAsFactors=F)

#Participationg rate (Ryan et al 2017)
#Part.rate.hist=30
#Part.rate.89=26.6  
#Part.rate.00=28.5
#Part.rate=mean(c(Part.rate.hist,Part.rate.89,Part.rate.00))  #mean fishing participating rate
#MISSING: update participation rates with figures from Eva for recent years. 
#         Karina reported 26.5%
#        Eva: "The participation rate for 18/19 was about 26%. It was below 30% in 
#               the last several years. The average over the last 10 years (2009/10 
#               to 2018/19) was about 29.8%.”"
Part.rate=29.8

# 2 -------------------I-Survey------------------------------------
Rec.fish.catch=Rec.fish.catch%>%
      mutate(Bioregion=ifelse(Bioregion%in%c("Gasconye","Gascoyne","Gascoyne Coast"),
                              "Gascoyne Coast",Bioregion))%>%
      rename(FINYEAR=FinYear)

RSE.weight=Rec.fish.catch%>%
  dplyr::select(FINYEAR,Bioregion,Common.Name,RSE)  

Rec.fish.catch=Rec.fish.catch%>% 
    dplyr::select(Common.Name,Kept.Number,Rel.Number,Bioregion,FINYEAR)

#Rec catch weights and assumed post capture mortality of released individuals
# Source of data: 
  #Average weight:
#   average weight sourced from Smallwood et al 2018 FRR 278 (Appendix 4)
Mn.w.gum=4.2
Mn.w.whi=3.7
Mn.w.whalr=5.4
Mn.w.wobi=6.7
Mn.w.ray=.5

  #Reported PCM (or PCS, post capture survival):
#   scalloped HH PCS=0.4 (Gulak et al 2015)
#   Gummy PCS=0.9  (Frick et al 2010)
#   School PCS=0.9  (Rogers et al 2017)
#   Dusky PCS=0.75   (Gallagher et al 2014, though this is at vessel survival)
#   Tiger PCS=0.95   (Gallagher et al 2014, though this is at vessel survival)
#   Blue PCS= 0.8   (Gallagher et al 2014, though this is at vessel survival)
#   Oceanic white tip PCS=0.75   (Gallagher et al 2014, though this is at vessel survival)
#  Wobbegongs PCM=0   #(Ellis et al 2016, average at vessel survival for LL)
#  Whalers PCM=0.273   #(Ellis et al 2016, average at vessel survival for LL)
#  Dasyatids PCM=0.156   #(Ellis et al 2016, average at vessel survival for LL)
#  Greynurse  PCM=0.13   #(Ellis et al 2016, average at vessel survival for LL)

# Given the lack of PCM info (but expected low PCM) a precautious value is 
#  assumed (this is considered in Sensitivity Scenarios)
Asmd=0.3  

Small.shrk=c('Dogfishes',"Gummy Sharks","Pencil Shark","Nervous Shark","Port Jackson Shark",
             "Sawsharks","School Shark","Sliteye Shark","Whitetip Reef Shark")
AVG.WT=data.frame(
  Common.Name=c(
    "Australian Blacktip Shark","Bignose Shark","Blacktip Reef Shark","Blue Shark","Bronze Whaler","Bull Shark",
    "Dogfishes","Dusky Whaler","Grey Reef Shark","Greynurse Shark","Gummy Sharks","Lemon Shark","Nervous Shark",
    "Oceanic Whitetip Shark","Pencil Shark","Pigeye Shark","Port Jackson Shark","Rays & Skates","Sandbar Shark",
    "Sawfishes","Sawsharks","Scalloped hammerhead","School Shark","Silky Shark","Silvertip Shark","Sliteye Shark",
    "Smooth hammerhead","Spinner Shark","Tawny Shark","Thresher Shark","Tiger Shark","Western Shovelnose Ray",
    "Whiskery Shark","Whitetip Reef Shark","Wobbegong","Zebra Shark"),
  AVG.wt=rep(Mn.w.whalr,36), 
  PCM.rec=rep(Asmd,36))%>%
          mutate(AVG.wt=ifelse(Common.Name%in%Small.shrk,Mn.w.gum,
                        ifelse(Common.Name=="Whiskery Shark",Mn.w.whi,
                        ifelse(Common.Name=="Rays & Skates",Mn.w.ray,
                        ifelse(Common.Name=="Wobbegong",Mn.w.wobi,
                        AVG.wt)))),
                PCM.rec=ifelse(Common.Name=="Port Jackson Shark",.05,
                        ifelse(Common.Name%in%c("Scalloped hammerhead","Smooth hammerhead"),0.6,
                        ifelse(Common.Name=="Gummy Sharks",.1,
                        ifelse(Common.Name=="Rays & Skates",0.16,
                        ifelse(Common.Name=="Wobbegong",0.05,
                        ifelse(Common.Name=="Greynurse Shark",0.13,
                        ifelse(Common.Name=="School Shark",0.1,
                               PCM.rec))))))))%>%  
           arrange(Common.Name)
AVG.WT$Common.Name=as.character(AVG.WT$Common.Name)

#fix species names
Rec.fish.catch=Rec.fish.catch%>%
        mutate(Common.Name=ifelse(Bioregion%in%c("Gascoyne","Gascoyne Coast","North Coast") & 
                                         Common.Name=="Bronze Whaler","Dusky Whaler",
                           ifelse(Common.Name%in%c('Whaler Shark – other',"Whaler & Weasel Sharks"),'Whaler Sharks',
                           ifelse(Common.Name=='Hammerhead Shark','Hammerhead Sharks',
                           ifelse(Common.Name%in%c('Other Rays Skates','Unspecified Rays Skates'),'Rays & Skates',
                           ifelse(Common.Name%in%c('Other Sharks','Unspecified Shark',"Sharks"),'Other Shark',
                                  Common.Name))))))

#reapportion 'whaler sharks' among all reported whaler species in Isurvey
Reported.whalers=c("Dusky Whaler","Sandbar Shark","Bronze Whaler","Tiger Shark",
               "Blacktip Reef Shark","Lemon Shark","Whitetip Reef Shark")
Whaler.prop=Rec.fish.catch%>%
                filter(Common.Name%in%Reported.whalers)%>%
                group_by(Common.Name,Bioregion,FINYEAR)%>%
                summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
                data.frame
Whaler.prop$Kept.Number.prop=Whaler.prop$Kept.Number/sum(Whaler.prop$Kept.Number)
Whaler.prop$Rel.Number.prop=Whaler.prop$Rel.Number/sum(Whaler.prop$Rel.Number)
Whaler.prop=Whaler.prop%>%dplyr::select(-c(Kept.Number,Rel.Number))

Reap.whalers=c('Other Whaler',"Whaler Sharks","Whaler Shark - other")
Whaler.Sharks=Rec.fish.catch%>%
                filter(Common.Name%in%Reap.whalers)%>%
                summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
                data.frame
Whaler.reap=cbind(Whaler.prop,Whaler.Sharks)%>%
                mutate(Kept.Number=Kept.Number*Kept.Number.prop,
                       Rel.Number=Rel.Number*Rel.Number.prop)%>%
        dplyr::select(names(Rec.fish.catch))
     
Rec.fish.catch=Rec.fish.catch%>%
                filter(!Common.Name%in%Reap.whalers)
Rec.fish.catch=rbind(Rec.fish.catch,Whaler.reap)

#reapportion 'other sharks' among all reported shark species                            
Shark.prop=Rec.fish.catch%>%
        filter(!Common.Name%in%c("Other Shark","Western Shovelnose Ray",
                                 "Rays & Skates","Other Rays and Skates"))%>%
        group_by(Common.Name,Bioregion,FINYEAR)%>%
        summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
        data.frame
Shark.prop$Kept.Number.prop=Shark.prop$Kept.Number/sum(Shark.prop$Kept.Number)
Shark.prop$Rel.Number.prop=Shark.prop$Rel.Number/sum(Shark.prop$Rel.Number)
Shark.prop=Shark.prop%>%dplyr::select(-c(Kept.Number,Rel.Number))
Other.Sharks=Rec.fish.catch%>%
               filter(Common.Name%in%c("Other Shark"))%>%
        summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
        data.frame
Other.reap=cbind(Shark.prop,Other.Sharks)%>%
        mutate(Kept.Number=Kept.Number*Kept.Number.prop,
               Rel.Number=Rel.Number*Rel.Number.prop)%>%
        dplyr::select(names(Rec.fish.catch))
Rec.fish.catch=Rec.fish.catch%>%
        filter(!Common.Name=="Other Shark")
Rec.fish.catch=rbind(Rec.fish.catch,Other.reap)

Rec.fish.catch.alone=Rec.fish.catch   #keep copy

# 3 -------------------Shore-based------------------------------------
  #get  shore:boat ratio
Shore.to.Boat.ratio=Shore.based$Total[2]/Shore.based$Total[1]

  #add shore-based to I-survey
Rec.fish.catch=Rec.fish.catch%>%
                  mutate(Kept.Number=Kept.Number+Shore.to.Boat.ratio*Kept.Number,
                         Rel.Number=Rel.Number+Shore.to.Boat.ratio*Rel.Number)



# 4 -------------------Charter boats------------------------------------
Charter=Charter%>%data.frame
charter.shk.sp=unique(Charter$Species.Common.Name)
charter.shk.sp=charter.shk.sp[grep('Shark|Whaler|Hammerhead',charter.shk.sp)]

Charter=Charter%>%
  filter(Species.Common.Name%in%charter.shk.sp)%>%
  rename(Common.Name=Species.Common.Name,
         Bioregion=TO.Zone,
         Kept.Number=Fish.Kept.Count,
         Rel.Number=Fish.Released.Count)%>%
  mutate(Mn=month(Month),
         Yr=year(Calendar.Year),
         FINYEAR=ifelse(Mn<6,Yr-1,Yr),
         FINYEAR=paste(FINYEAR,substr(FINYEAR+1,3,4),sep="-"))%>%
  dplyr::select(names(Rec.fish.catch))%>%
         mutate(Common.Name=ifelse(Common.Name=="Gulper sharks, Sleeper Sharks & Dogfishes",
                                   "Dogfishes",Common.Name))

# reapportion "Whaler & Weasel Sharks" among all reported whaler species
Whaler.prop=Charter%>%
  filter(Common.Name%in% 
           c("Australian Blacktip Shark","Bignose Shark","Blacktip Reef Shark",                          
             "Blue Shark","Bronze Whaler","Bull Shark","Dusky Whaler","Grey Reef Shark",
             "Lemon Shark","Nervous Shark","Oceanic Whitetip Shark",
             "Pigeye Shark","Sandbar Shark","Silky Shark","Silvertip Shark",
             "Sliteye Shark","Spinner Shark","Tiger Shark","Whitetip Reef Shark"))%>%
  group_by(Common.Name,Bioregion,FINYEAR)%>%
  summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
  data.frame
Whaler.prop$Kept.Number.prop=with(Whaler.prop,Kept.Number/sum(Kept.Number))
Whaler.prop$Rel.Number.prop=with(Whaler.prop,Rel.Number/sum(Rel.Number))
Whaler.prop=Whaler.prop%>%dplyr::select(-c(Kept.Number,Rel.Number))

Whaler.Sharks=Charter%>%
  filter(Common.Name=="Whaler & Weasel Sharks")%>%
  summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
  data.frame
Whaler.reap=cbind(Whaler.prop,Whaler.Sharks)%>%
  mutate(Kept.Number=Kept.Number*Kept.Number.prop,
         Rel.Number=Rel.Number*Rel.Number.prop)%>%
  dplyr::select(names(Rec.fish.catch))

Charter=Charter%>%
  filter(!Common.Name=="Whaler & Weasel Sharks")
Charter=rbind(Charter,Whaler.reap)

# Change bull to pigeye shark
Charter=Charter%>%
          mutate(Common.Name=ifelse(Common.Name=='Bull Shark','Pigeye Shark',Common.Name))

# reapportion "Hound Sharks" among all reported hound shark species
Hound.prop=Charter%>%
  filter(Common.Name%in% 
           c("Gummy Shark","Pencil Shark","School Shark","Whiskery Shark"))%>%
  group_by(Common.Name,Bioregion,FINYEAR)%>%
  summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
  data.frame
Hound.prop$Kept.Number.prop=with(Hound.prop,Kept.Number/sum(Kept.Number))
Hound.prop$Rel.Number.prop=with(Hound.prop,Rel.Number/sum(Rel.Number))
Hound.prop=Hound.prop%>%dplyr::select(-c(Kept.Number,Rel.Number))

Hound.Sharks=Charter%>%
  filter(Common.Name=="Hound Sharks")%>%
  summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
  data.frame
Hound.reap=cbind(Hound.prop,Hound.Sharks)%>%
  mutate(Kept.Number=Kept.Number*Kept.Number.prop,
         Rel.Number=Rel.Number*Rel.Number.prop)%>%
  dplyr::select(names(Rec.fish.catch))

Charter=Charter%>%
  filter(!Common.Name=="Hound Sharks")
Charter=rbind(Charter,Hound.reap)

# add Charter to Rec.fish.catch
Rec.fish.catch=rbind(Rec.fish.catch,Charter)


# standardise names and fix mis reporting
Rec.fish.catch=Rec.fish.catch%>%
    mutate(Common.Name=ifelse(Common.Name=="Blacktip reef shark","Blacktip Reef Shark",
                       ifelse(Common.Name=="Whitetip reef shark","Whitetip Reef Shark",
                       ifelse(Common.Name=='Gulper sharks, Sleeper Sharks & Dogfishes','Dogfishes',
                       ifelse(Common.Name=='Blind, Nurse, Carpet & Zebra Sharks','Zebra Shark',
                       ifelse(Common.Name%in%c('Other Rays and Skates',"Western Shovelnose Ray"),
                              'Rays & Skates',
                       ifelse(Common.Name=='Sawshark','Sawsharks',
                       ifelse(Common.Name=='Wobbegong','Wobbegongs',
                       ifelse(Common.Name=="School Shark" & 
                                  Bioregion%in%c('North Coast'),'Gummy Sharks',
                       ifelse(Common.Name=="Blacktip Reef Shark" & 
                                  Bioregion%in%c('South Coast','West Coast'),'Spinner Shark',
                        ifelse(Common.Name=="Gummy Shark","Gummy Sharks",Common.Name)))))))))))
          
AVG.WT=AVG.WT%>%
        mutate(Common.Name=ifelse(Common.Name=='Wobbegong','Wobbegongs',Common.Name))  
    

  
# 5 -------------------Reapportion 'hammerheads' across all rec fishing sources------------------------
Rec.fish.catch=Rec.fish.catch%>%
         mutate(Common.Name=ifelse(Common.Name=='Hammerhead Sharks' & Bioregion%in%c("Gascoyne",
                                  "North Coast","Gascoyne Coast"),"Scalloped hammerhead",
                            ifelse(Common.Name=='Hammerhead Sharks' &
                                   Bioregion%in%c("West Coast","South Coast"),"Smooth hammerhead",
                            ifelse(Common.Name=='Scalloped Hammerhead',"Scalloped hammerhead",
                            ifelse(Common.Name=='Smooth Hammerhead',"Smooth hammerhead",
                            Common.Name)))))



# 6 -------------------Reconstruct time series------------------------------------

#reconstruct size of human population fishing
dummy=rbind(cbind(Year=c(1940,1950,1960),Population=c(473300,557100,722100)),WA.population)
mod=loess(Population~Year,data=dummy)
Historic.pop=predict(mod,newdata = data.frame(Year=1941:1970))
WA.population=rbind(cbind(Year=1941:1970,Population=round(Historic.pop)),
                    WA.population)

Fishing.population=(Part.rate/100)*WA.population$Population
Fishing.population=Fishing.population/Fishing.population[match(2011,WA.population$Year)] #relative to 2011-12
Fishing.population=data.frame(
    Size=Fishing.population[1:(length(Fishing.population)-1)],
    FinYear=paste(WA.population$Year[1:(length(WA.population$Year)-1)],"-",
                  substr(WA.population$Year[2:length(WA.population$Year)],start=3,stop=4),sep=""))

  #1. calculate total catch in weight (kg) for each species for Isurvey years
fn.rec=function(DAT,PCM.scen,Wght.scen)  
{
   AGG=DAT%>%
      left_join(AVG.WT,by="Common.Name")%>%
      mutate(LIVEWT.c=ceiling((Kept.Number+Rel.Number*PCM.rec*PCM.scen)*AVG.wt*Wght.scen))%>%
      group_by(FINYEAR,Common.Name,Bioregion)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c))%>%
      data.frame
  LisT=vector('list',length(unique(AGG$Common.Name)))
  names(LisT)=unique(AGG$Common.Name)
  for(i in 1:length(LisT))
  {
    LisT[[i]]=subset(AGG,Common.Name==names(LisT)[i],
                     select=c(Common.Name,LIVEWT.c,FINYEAR,Bioregion))
  }
  return(LisT)
}
Rec.ktch=vector('list',nrow(Scenarios))
names(Rec.ktch)=Scenarios$Scenario
for(s in 1:length(Rec.ktch))
{
  #select only Isurvey years
  dd=subset(Rec.fish.catch,FINYEAR%in%I.survey.years)  
  
  #add rare species not in Isurvey but in Charter (though not during Isurvey years)
  Uni=unique(Rec.fish.catch$Common.Name)
  IId=which(!Uni%in%unique(dd$Common.Name))
  if(length(IId)>0)
  {
    dd.IID=Rec.fish.catch%>%
            filter(Common.Name%in%Uni[IId])%>%
            group_by(Common.Name,Bioregion)%>%
            summarise(Kept.Number=sum(Kept.Number),
                      Rel.Number=sum(Rel.Number))%>%
            mutate(FINYEAR='2011-12')%>%
            data.frame%>%
            arrange(colnames(dd))
    dd=rbind(dd,dd.IID)
  }
  
  #run reconstruction
  Rec.ktch[[s]]=fn.rec(DAT=dd,
                       PCM.scen=Scenarios$PCM[s],
                       Wght.scen=Scenarios$Weight[s])   
}

  #2. reconstruct total catch (in kg) time series
back.fill=function(dat)
{
  Regns=unique(dat$Bioregion)
  Dummy=vector('list',length(Regns))
  for (d in 1:length(Dummy))
  {
     a=subset(dat,Bioregion==Regns[d])
     Mat=matrix(NA,nrow=length(Fishing.population$FinYear),ncol=2)
     rse.W=subset(RSE.weight,Common.Name==a$Common.Name[1] & Bioregion==a$Bioregion[1])
     if(nrow(rse.W)==0) rse.w=mean(RSE.weight$RSE)
     if(nrow(rse.W)>0)  rse.w=rse.W$RSE
     if(!length(a$LIVEWT.c)==length(rse.w))
     {
        if(nrow(rse.W)==0) rse.w=rep(rse.w,length(a$LIVEWT.c))
        if(nrow(rse.W)>0)
        {
          misn.id=which(!a$FINYEAR%in%rse.W$FINYEAR)
          misn=rse.W[1:length(misn.id),]
          misn$FINYEAR=a$FINYEAR[misn.id]
          misn$RSE=max(rse.W$RSE)
          rse.W=rbind(rse.W,misn)%>%arrange(FINYEAR)
          rse.w=rse.W$RSE
        }
     }
     #weighted mean by RSE
     VAL=weighted.mean(a$LIVEWT.c,w=1/rse.w)   
      
     Mat[,2]=VAL*Fishing.population$Size
     Mat=as.data.frame(Mat)
     names(Mat)=c("FINYEAR","LIVEWT.c")
     Mat$FINYEAR=Fishing.population$FinYear
     Mat$zone=Regns[d]
     Mat$Common.Name=dat$Common.Name[1]
     Dummy[[d]]=Mat
  }
  Dummy=do.call(rbind,Dummy)
  return(Dummy)
}
for(s in 1:length(Rec.ktch))
{
  DumY=vector('list',length(Rec.ktch[[s]]))
  names(DumY)=names(Rec.ktch[[s]])
  for(i in 1:length(DumY))DumY[[i]]=back.fill(dat=Rec.ktch[[s]][[i]])
  DumY=do.call(rbind,DumY)
  row.names(DumY)=NULL 
  Rec.ktch[[s]]=DumY
}
Rec.ktch.Upper=Rec.ktch$High
Rec.ktch.Lower=Rec.ktch$Low
Rec.ktch=Rec.ktch$`Base Case`




# 7 -------------------EXPORT CATCH DATA------------------------------------
fn.out=function(d,NM)
{
  write.csv(d,paste(handl_OneDrive('Analyses/Data_outs/'),NM,sep=""),row.names = F)
}
fn.out(d=Rec.ktch,NM='recons_recreational.csv')



# 8 -------------------Report------------------------------------
if(Do.recons.rec.fishn.paper=="YES")
{
  hndl.out=handl_OneDrive("Analyses\\Reconstruction_catch_recreational\\")
  
  #export summarised ISurvey and Charter logbooks
  Tab.2=Rec.fish.catch%>%
          group_by(Common.Name)%>%
          summarise(Kept=round(sum(Kept.Number)),
                     Rel=round(sum(Rel.Number)))%>%
          data.frame
  write.csv(Tab.2,paste(hndl.out,"Table2.Table.Kept.Rel.csv",sep=''),row.names = FALSE)
  
  #export weight and PCS table
  Tab.1=AVG.WT%>%
    filter(Common.Name%in%unique(Rec.fish.catch$Common.Name))%>%
    left_join(Scien.nm,by='Common.Name')%>%
    arrange(Common.Name)%>%
    dplyr::select(Common.Name,Scientific.name,AVG.wt,PCM.rec)
  write.csv(Tab.1,paste(hndl.out,"Table2.Table.wght.PCS.csv",sep=''),row.names = FALSE)
  
  #1. Species proportions (by number) for each bioregion based on original Isurvey and Charter data
  Rec.fish.catch.alone$source='Isurvey'
  Charter$source='Charter'
  Dat.show=rbind(Rec.fish.catch.alone,Charter)    
  Dat.show=Dat.show%>%
    mutate(Common.Name=ifelse(Common.Name=="Blacktip reef shark","Blacktip Reef Shark",
                       ifelse(Common.Name=="Whitetip reef shark","Whitetip Reef Shark",
                       ifelse(Common.Name=='Gulper sharks, Sleeper Sharks & Dogfishes','Dogfishes',
                       ifelse(Common.Name=='Blind, Nurse, Carpet & Zebra Sharks','Zebra Shark',
                       ifelse(Common.Name%in%c('Other Rays and Skates',"Western Shovelnose Ray"),
                              'Rays & Skates',
                       ifelse(Common.Name=='Sawshark','Sawsharks',
                       ifelse(Common.Name=='Wobbegong','Wobbegongs',
                       ifelse(Common.Name=="School Shark" & 
                                    Bioregion%in%c('North Coast'),'Gummy Sharks',
                       ifelse(Common.Name=="Blacktip Reef Shark" & 
                                    Bioregion%in%c('South Coast','West Coast'),'Spinner Shark',
                       ifelse(Common.Name=="Gummy Shark","Gummy Sharks",Common.Name)))))))))))
  
  Dat.show=Dat.show%>%
    mutate(Common.Name=ifelse(Common.Name=='Hammerhead Sharks' & Bioregion%in%c("Gascoyne",
                                  "North Coast","Gascoyne Coast"),"Scalloped hammerhead",
                       ifelse(Common.Name=='Hammerhead Sharks' &
                                    Bioregion%in%c("West Coast","South Coast"),"Smooth hammerhead",
                       ifelse(Common.Name=='Scalloped Hammerhead',"Scalloped hammerhead",
                       ifelse(Common.Name=='Smooth Hammerhead',"Smooth hammerhead",
                       Common.Name)))))  
  
  Ktch.by.sp.zn=Dat.show %>%
    filter(FINYEAR%in%I.survey.years)%>%
    mutate(Tot=Kept.Number+ Rel.Number)%>%
    group_by(Bioregion, Common.Name) %>%
    summarise(n = sum(Tot,na.rm=T))
  
  Prop= Ktch.by.sp.zn%>%
    mutate(freq = n / sum(n,na.rm=T))%>%
    dplyr::select(-n)%>%
    spread(Bioregion,freq,fill=0)%>%
    dplyr::select(Common.Name,'North Coast','Gascoyne Coast','West Coast','South Coast')
  
  write.csv(Prop,paste(hndl.out,"Catch.Proportion.by.weight.bioregion.csv",sep=''),row.names = FALSE)

  Ktch.by.sp.zn=Ktch.by.sp.zn%>%
                  spread(Bioregion,n,fill=0)
  Sp.tot=rowSums(Ktch.by.sp.zn[,-1])
  Ktch.by.sp.zn=Ktch.by.sp.zn%>%
    mutate(Sp.tot=Sp.tot)%>%
    mutate_at(vars(-c(Common.Name,Sp.tot)), funs(. / Sp.tot))
  
  Ktch.by.sp.zn=Ktch.by.sp.zn[order(Ktch.by.sp.zn$'North Coast',
                                    Ktch.by.sp.zn$'Gascoyne Coast',
                                    Ktch.by.sp.zn$'West Coast'),]
  Sp.tot=Ktch.by.sp.zn$Sp.tot
  Ktch.by.sp.zn=Ktch.by.sp.zn%>%dplyr::select(Common.Name,'North Coast',
                        'Gascoyne Coast','West Coast','South Coast')
  
  COls=c('grey90','grey70','grey40','grey10')
  names(COls)=c('North Coast','Gascoyne Coast','West Coast','South Coast')
  Legdns=capitalize(tolower(Ktch.by.sp.zn$Common.Name))
  Legdns=ifelse(Legdns=="Port jackson shark","Port Jackson shark",Legdns)
  
  tiff(file=paste(hndl.out,"Fig2. Proportion.tiff",sep=''),width=2400,height=2400,
       units="px",res=300,compression="lzw")
  par(mar=c(2.5,8,.5,3),oma=c(.1,.1,.1,.1),mgp=c(.1,.3,0),las=1,xpd=TRUE)
  s=barplot(t(as.matrix(Ktch.by.sp.zn[,-1])),horiz=T,col=COls,
          names.arg = Legdns,cex.names=.8,
          legend.text=T,args.legend=list(x=1,y=nrow(Ktch.by.sp.zn)+8,
                                         bty='n',horiz=T,cex=1.1))
  text(.99,s,paste("n=",round(Sp.tot)),pos=4,cex=.8)
  mtext("Proportion",1,1.5,cex=1.5)
  dev.off()
  
  #2. Multivariate
    #MDS comparing bioregions from I-survey and Charter (using same years)
  library(vegan)
  library(network)
  MDS=Dat.show %>%
    filter(FINYEAR%in%I.survey.years)%>%
    mutate(Tot=round(Kept.Number+ Rel.Number))%>%
    group_by(Bioregion, FINYEAR,source,Common.Name) %>%
    summarise(n = sum(Tot,na.rm=T))%>%
    spread(Common.Name,n,fill=0)
  Factors=MDS%>%dplyr::select(Bioregion,source,FINYEAR)
  Factors$FINYEAR=substr(Factors$FINYEAR,1,4)
  col.nms=colnames(MDS)[-(1:ncol(Factors))]
  MDS=as.matrix(MDS[,-(1:ncol(Factors))])
  MDS=sqrt(MDS) #apply sqrt transf
  #MDS=MDS^(1/4) #4th root trans
  colnames(MDS)=col.nms
  NMDS=metaMDS(MDS,k=2,trymax=100)
  
  dummi=data.frame(col=COls,Bioregion=names(COls))
  dummi$col=as.character(dummi$col)
  dummi$Bioregion=as.character(dummi$Bioregion)
  Factors=left_join(Factors,dummi,by='Bioregion')
  Factors$Shape=with(Factors,ifelse(source=='Charter',21,24))
  
  tiff(file=paste(hndl.out,"Fig3. MDS.tiff",sep=''),width=2400,height=2400,
       units="px",res=300,compression="lzw")
  par(mar=c(3,1,3,1),oma=c(.1,.1,.1,.1),mgp=c(.1,.3,0),las=1,xpd=TRUE)
  plot(NMDS$points[,1],NMDS$points[,2],bg=Factors$col,cex=2,
       pch=Factors$Shape,ann=F,xaxt='n',yaxt='n')
  #text(NMDS$points[,1],NMDS$points[,2],Factors$FINYEAR,col=as.color(Factors$Bioregion),cex=.8)
  legend('bottomright',paste("2D Stress=",round(NMDS$stress,2)),bty='n')
  legend(x=-.4,y=-.96,c("Private-boat","Charter-boat"),pch=c(24,21),pt.cex=1.5,bty='n',horiz = T,cex=1.25)
  legend(x=-1.75,y=1.4,names(COls),pch=22,pt.bg=COls,pt.cex=2,bty='n',horiz = T,cex=1.2)
  dev.off()
  
  
    #Simper to find species making the difference
  Factors.env=subset(Factors,select=c(Bioregion,source))
  Factors.env$Bioregion=factor(Factors.env$Bioregion,levels=names(COls))
  Factors.env$source=factor(Factors.env$source)
  
  SIMPER.bioregion <- with(Factors.env, simper(MDS, Bioregion))
  SIMPER.bioregion=summary(SIMPER.bioregion)
  
  SIMPER.source <- with(Factors.env, simper(MDS, source))
  SIMPER.source=summary(SIMPER.source)
  
  capture.output(SIMPER.bioregion, file = paste(hndl.out,'Appendix1.SIMPER.bioregion.txt',sep=''))
  capture.output(SIMPER.source, file = paste(hndl.out,'Appendix1.SIMPER.source.txt',sep=''))
  
  
  #3. Temporal trends in reconstructed catches      
  source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R'))
  Rec.ktch=Rec.ktch%>%
    mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
    arrange(Common.Name)
  Rec.sp=unique(Rec.ktch$Common.Name)
  LWD=3
  Legdns=capitalize(tolower(Rec.sp))
  Legdns=ifelse(Legdns=="Port jackson shark","Port Jackson shark",Legdns)
  
  tiff(file=paste(hndl.out,"Fig4. Time series.tiff",sep=''),width=2400,height=2000,
       units="px",res=300,compression="lzw") 
  smart.par(n.plots=length(Rec.sp),MAR=c(2,2,1,1.2),OMA=c(1.75,2.5,.5,.1),MGP=c(1,.5,0))
  for(i in 1:length(Rec.sp))
  {
    d=Rec.ktch%>%
      filter(Common.Name==Rec.sp[i])%>%
      mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
      group_by(year)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c/1000))
    
    d.up=Rec.ktch.Upper%>%
      filter(Common.Name==Rec.sp[i])%>%
      mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
      group_by(year)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c/1000))
    
    d.low=Rec.ktch.Lower%>%
      filter(Common.Name==Rec.sp[i])%>%
      mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
      group_by(year)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c/1000))
    
    plot(sort(unique(d$year)),sort(unique(d$year)),col='transparent',cex=.8,ann=F,
         ylim=c(0,max(d.up$LIVEWT.c,na.rm=T)))
    if(nrow(d)>0) lines(d$year,d$LIVEWT.c,col="grey55",lwd=LWD)
    if(nrow(d.up)>0) lines(d.up$year,d.up$LIVEWT.c,col="grey20",lwd=LWD)
    if(nrow(d.low)>0) lines(d.low$year,d.low$LIVEWT.c,col="grey80",lwd=LWD)
    #ySeq=seq(min(d.up$LIVEWT.c),max(d.up$LIVEWT.c),length.out = 3)
    #if(trunc(ySeq)[1]>0) ySeq.lab=round(ySeq.lab) else
    #ySeq.lab=round(ySeq,3)
    #axis(2,ySeq,ySeq.lab)
    
    nm=Legdns[i]
  #  nm=ifelse(nm=="Wobbegong","Wobbegongs",nm)
    mtext(paste(nm),3,line=0.2,cex=0.8)  
  }
  plot.new()
  legend('left',c("High","Base case","Low"),lty=1,col=c("grey20","grey55","grey80"),
         lwd=LWD,bty='n',cex=1.3)
  mtext("Financial year",1,line=0.5,cex=1.5,outer=T)
  mtext("Total harvest (tonnes)",2,las=3,line=.7,cex=1.5,outer=T)
  dev.off()
  
  
  #4. Bioregion map
  hndl.map=handl_OneDrive("Data/Mapping/Bioregions")
  library(rgdal)
  library(PBSmapping)
  data(worldLLhigh)
  source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R"))
  Bioregions=readOGR(paste(hndl.map,"Bioregions.shp",sep="/"), layer="Bioregions") 
  
  South.WA.long=c(109,129)
  South.WA.lat=c(-38,-12)
  Xlim=c(109,129)
  Ylim=South.WA.lat
  
  tiff(file=paste(hndl.out,"Fig1. Map.tiff",sep=''),width=1600,height=2400,
       units="px",res=300,compression="lzw") 
  
  par(mar=c(1,1,.5,.5),oma=c(3,3,1,.3),las=1,mgp=c(.04,.6,0))
  plot(1,xlim=Xlim,ylim=Ylim,xlab="",ylab="",axes=F,main="")
  plot(Bioregions,add=T)
  polygon(WAcoast$Longitude,WAcoast$Latitude, col="grey70")
  axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),-seq(round(Ylim[1]),round(Ylim[2]),2),cex.axis=1.25)
  axis(side = 1, seq(South.WA.long[1],South.WA.long[2],2), 
       labels =seq(South.WA.long[1],South.WA.long[2],2), tck = -.015,cex.axis=1.25)
  box()
  mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1.2,font=1,las=0,cex=1.35,outer=T)
  mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=1,font=1,las=0,cex=1.35,outer=T)
  
  text(120,-16,"North Coast",cex=1.5,srt=45)
  text(112.5,-22.8,"Gascoyne Coast",srt=75,cex=1.5)
  text(113,-31,"West",cex=1.5)
  text(113,-32.5,"Coast",cex=1.5)
  text(121,-35.8,"South Coast",cex=1.5)
  text(122,-26,"Western Australia",cex=2)
  
  # #Australia
  par(fig = c(.05, .35, .7, 1), mar=c(0,0,0,0), new=TRUE)
  OZ.lat=c(-44.5,South.WA.lat[2]);OZ.long=c(South.WA.long[1],155)
  plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
          col='black',tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  polygon(WAcoast$Longitude,WAcoast$Latitude, col="grey50")
  text(135,-25,("Australia"),col="white", cex=1.3)
  dev.off()
}
