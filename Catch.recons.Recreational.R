#-- Script for reconstructing time series of recreational catch of sharks in WA

#notes: This script uses I-Survey (boat-based) point estimates to reconstruct
#       time series considering rec fishing participation rates and WA 
#       population growth
#       Shore-based fishing is added as a proportion of boat-based fishing
#       Charter boat fishing is also added

#For each new I-survey, get data from Karina
# Also update Charter boat data each year

rm(list=ls(all=TRUE))

library(tidyverse)
library(readxl)
library(lubridate)


Do.recons.rec.fishn.paper="NO"


#Catch reconstruction scenarios
Scenarios=data.frame(Scenario=c('Base Case','Upper 95%','Lower 95%'),
                     Time.series=c('mean','upper','lower'),
                     PCM=c(1,1.5,.5))

# 1 -------------------DATA SECTION------------------------------------

# I-Survey
Rec.hndl="C:/Matias/Data/Catch and Effort/Recreational/I.Survey."
#Rec.fish.catch.2011.12=read.csv(paste(Rec.hndl,"2011_12.csv",sep=''),stringsAsFactors=F) #Ryan et al 2013    
#Rec.fish.catch.2013.14=read.csv(paste(Rec.hndl,"2013_14.csv",sep=''),stringsAsFactors=F) #Ryan et al 2015 
#Rec.fish.catch.2015.16=read.csv(paste(Rec.hndl,"2015_16.csv",sep=''),stringsAsFactors=F) #Ryan et al 2017 
#Rec.fish.catch=rbind(Rec.fish.catch.2011.12,Rec.fish.catch.2013.14,Rec.fish.catch.2015.16)
Rec.fish.catch=read.csv(paste(Rec.hndl,"csv",sep=''),stringsAsFactors=F)
Rec.fish.catch=Rec.fish.catch%>%
            mutate(FinYear=paste(2000+as.numeric(substr(Year,1,2)),substr(Year,3,4),sep="-"),
                   Bioregion=paste(sapply(strsplit(Rec.fish.catch$RegionReportingGroupName, "_"), "[", 2),
                                   sapply(strsplit(Rec.fish.catch$RegionReportingGroupName, "_"), "[", 3)),
                   Common.Name=Lowlevelgrouping,
                   Scientific.Name=ScientificName,
                   Kept.Number=Kept,
                   Kept.Number.se=se.Kept,
                   Rel.Number=Released,
                   Rel.Number.se=se.Released,
                   Caught.Number=Total,
                   Caught.Number.se=se.Total)%>%
      dplyr::select(Common.Name,Scientific.Name,Kept.Number,Kept.Number.se,Rel.Number,
                   Rel.Number.se,Caught.Number,Caught.Number.se,Bioregion,FinYear)
I.survey.years=unique(Rec.fish.catch$FinYear)

# Shore-based
#notes from Karina:
  #statewide estimated kept, released & total recreational catch for sharks & rays in 2000_01
  #K_SE = standard error associated with Kept
  # K_RSE = relative standard error associated with Kept
  # K_HHs = number of households that reported a Kept catch 
  # R = released, T = Total 
Shore.based=read.csv("C:/Matias/Data/Catch and Effort/Recreational/statewide shark 2000_01.csv",stringsAsFactors=F)

# Charter boats
Charter=read_excel("C:\\Matias\\Data\\Catch and Effort\\Charter\\Charter.xlsx",sheet ='Data')

# WA population for rec catch recons (ABS)
#source: https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/3101.0Dec%202018?OpenDocument
WA.population=read.csv("C:/Matias/Data/AusBureauStatistics.csv",stringsAsFactors=F)

#Participationg rate (Ryan et al 2012)
Part.rate.hist=30
Part.rate.89=26.6  
Part.rate.00=28.5

# 2 -------------------I-Survey------------------------------------
Rec.fish.catch=Rec.fish.catch%>%
      mutate(Bioregion=ifelse(Bioregion%in%c("Gasconye","Gascoyne","Gascoyne Coast"),"Gascoyne Coast",Bioregion),
             FINYEAR=FinYear)%>% 
    dplyr::select(Common.Name,Kept.Number,Rel.Number,Bioregion,FINYEAR)

#Rec catch weights and assumed post capture mortality of released individuals
#Sources: 
#     average weight sourced from Smallwood et al 2018 FRR 278 (Appendix 4)
#     scalloped HH PCS=0.4 (Gulak et al 2015)
#     Gummy PCS=0.9  (Frick et al 2010)
#     School PCS=0.9  (Rogers et al 2017)
#     Dusky PCS=0.75   (Gallagher et al 2014, though this is at vessel survival)
#     Tiger=0.95   (Gallagher et al 2014, though this is at vessel survival)
#     Blue=0.8   (Gallagher et al 2014, though this is at vessel survival)
#     Oceanic white tip=0.75   (Gallagher et al 2014, though this is at vessel survival)
Asmd=0.3  #Precautious PCM
AVG.WT=data.frame(Common.Name=c("Bronze Whaler","Greynurse Shark","Gummy Sharks","Scalloped hammerhead",
                                "Smooth hammerhead","Sandbar Shark","Tiger Shark","Wobbegong",             
                                "Western Shovelnose Ray","Rays & Skates","Sawsharks",             
                                "Port Jackson Shark","Whiskery Shark","School Shark","Blacktip Reef Shark",   
                                "Dusky Whaler","Lemon Shark","Whitetip Reef Shark",
                                "Australian Blacktip Shark","Bignose Shark",                            
                                "Blue Shark","Bull Shark",                               
                                "Grey Reef Shark","Dogfishes",
                                "Nervous Shark","Oceanic Whitetip Shark",                    
                                "Pencil Shark","Pigeye Shark",                             
                                "Silky Shark","Silvertip Shark",                          
                                "Sliteye Shark","Spinner Shark",                            
                                "Tawny Shark","Thresher Shark",                        
                                "Zebra Shark",
                                "Sawfishes"),
                  AVG.wt=c(6.7,6.7,4.2,6.7,6.7,6.7,6.7,6.7,6.7,5,4.2,4.2,4.2,4.2,6.7,6.7,6.7,6.7,
                           rep(6.7,5),4.2,rep(6.7,3),4.2,rep(6.7,7),10), 
                  PCM.rec=c(Asmd,Asmd,.1,.7,.7,Asmd,Asmd,Asmd,Asmd,Asmd,Asmd,.05,Asmd,Asmd,Asmd,
                            Asmd,Asmd,Asmd,
                            rep(Asmd,5),.5,rep(Asmd,2),Asmd,rep(Asmd,9)))           
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

#reapportion 'whaler sharks' among all reported whaler species
Whaler.prop=Rec.fish.catch%>%
                filter(Common.Name%in%c("Dusky Whaler","Sandbar Shark","Bronze Whaler","Tiger Shark",
                                        "Blacktip Reef Shark","Lemon Shark","Whitetip Reef Shark"))%>%
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
        filter(!Common.Name%in%c("Other Shark","Western Shovelnose Ray","Rays & Skates","Other Rays and Skates"))%>%
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

# add to Rec.fish.catch
Rec.fish.catch=rbind(Rec.fish.catch,Charter)


# standardise names and fix mis reporting
Rec.fish.catch=Rec.fish.catch%>%
    mutate(Common.Name=ifelse(Common.Name=="Blacktip reef shark","Blacktip Reef Shark",
                       ifelse(Common.Name=="Whitetip reef shark","Whitetip Reef Shark",
                       ifelse(Common.Name=='Gulper sharks, Sleeper Sharks & Dogfishes','Dogfishes',
                       ifelse(Common.Name=='Blind, Nurse, Carpet & Zebra Sharks','Zebra Shark',
                       ifelse(Common.Name=='Other Rays and Skates','Rays & Skates',
                       ifelse(Common.Name=='Sawshark','Sawsharks',
                       ifelse(Common.Name=="School Shark" & 
                                  Bioregion%in%c('North Coast'),'Gummy Sharks',
                       ifelse(Common.Name=="Blacktip Reef Shark" & 
                                  Bioregion%in%c('South Coast','West Coast'),'Spinner Shark',
                        ifelse(Common.Name=="Gummy Shark","Gummy Sharks",Common.Name))))))))))
          
  
    
   
  
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

#reconstruct size of population fishing
Part.rate=mean(c(Part.rate.hist,Part.rate.89,Part.rate.00))  #mean fishing participating rate

dummy=rbind(cbind(Year=c(1940,1950,1960),Population=c(473300,557100,722100)),WA.population)
mod=loess(Population~Year,data=dummy)
Historic.pop=predict(mod,newdata = data.frame(Year=1941:1970))
WA.population=rbind(cbind(Year=1941:1970,Population=round(Historic.pop)),
                    WA.population)

fn.rec=function(DAT,PCM.scen)
{
        AGG=DAT%>%left_join(AVG.WT,by="Common.Name")%>%
                  mutate(LIVEWT.c=ceiling((Kept.Number+Rel.Number*PCM.rec*PCM.scen)*AVG.wt))%>%
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
Rec.ktch=vector('list',length(Scenarios))
names(Rec.ktch)=Scenarios$Scenario
for(s in 1:length(Rec.ktch))
{
  Rec.ktch[[s]]=fn.rec(DAT=subset(Rec.fish.catch,FINYEAR%in%I.survey.years),
                       PCM.scen=Scenarios$PCM[s])   
}

Fishing.population=(Part.rate/100)*WA.population$Population
Fishing.population=Fishing.population/Fishing.population[match(2011,WA.population$Year)] #relative to 2011-12
Fishing.population=data.frame(Size=Fishing.population[1:(length(Fishing.population)-1)],
                             FinYear=paste(WA.population$Year[1:(length(WA.population$Year)-1)],"-",
                                        substr(WA.population$Year[2:length(WA.population$Year)],start=3,stop=4),sep=""))

back.fill=function(dat,scen)
{
        Regns=unique(dat$Bioregion)
        Dummy=vector('list',length(Regns))
        for (d in 1:length(Dummy))
        {
                a=subset(dat,Bioregion==Regns[d])
                Mat=matrix(NA,nrow=length(Fishing.population$FinYear),ncol=2)
                if(scen=='mean') VAL=mean(a$LIVEWT.c)
                if(scen=='upper') VAL=quantile(a$LIVEWT.c,.975)
                if(scen=='lower') VAL=quantile(a$LIVEWT.c,.025)
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
  for(i in 1:length(DumY))DumY[[i]]=back.fill(dat=Rec.ktch[[s]][[i]],scen=Scenarios$Time.series[s])
  DumY=do.call(rbind,DumY)
  row.names(DumY)=NULL 
  Rec.ktch[[s]]=DumY
}
Rec.ktch.Upper=Rec.ktch$`Upper 95%`
Rec.ktch.Lower=Rec.ktch$`Lower 95%`
Rec.ktch=Rec.ktch$`Base Case`



# 7 -------------------Report------------------------------------
if(Do.recons.rec.fishn.paper=="YES")
{
  hndl.out="C:\\Matias\\Analyses\\Reconstruction_catch_recreational\\"
  
  #export weight and PCS table
  write.csv(AVG.WT[order(AVG.WT$Common.Name),],paste(hndl.out,"Appendix.Table.wght.PCS.csv",sep=''),row.names = FALSE)
  
  #1. Species proportions (by number) for each bioregion based on original Isurvey and Charter data
  Rec.fish.catch.alone$source='Isurvey'
  Charter$source='Charter'
  Dat.show=rbind(Rec.fish.catch.alone,Charter)    
  Dat.show=Dat.show%>%
    mutate(Common.Name=ifelse(Common.Name=="Blacktip reef shark","Blacktip Reef Shark",
                       ifelse(Common.Name=="Whitetip reef shark","Whitetip Reef Shark",
                       ifelse(Common.Name=='Gulper sharks, Sleeper Sharks & Dogfishes','Dogfishes',
                       ifelse(Common.Name=='Blind, Nurse, Carpet & Zebra Sharks','Zebra Shark',
                       ifelse(Common.Name=='Other Rays and Skates','Rays & Skates',
                       ifelse(Common.Name=='Sawshark','Sawsharks',
                       ifelse(Common.Name=="School Shark" & 
                                    Bioregion%in%c('North Coast'),'Gummy Sharks',
                       ifelse(Common.Name=="Blacktip Reef Shark" & 
                                    Bioregion%in%c('South Coast','West Coast'),'Spinner Shark',
                       ifelse(Common.Name=="Gummy Shark","Gummy Sharks",Common.Name))))))))))
  
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
    select(-n)%>%
    spread(Bioregion,freq,fill=0)%>%
    select(Common.Name,'North Coast','Gascoyne Coast','West Coast','South Coast')
  
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
  Ktch.by.sp.zn=Ktch.by.sp.zn%>%select(Common.Name,'North Coast',
                        'Gascoyne Coast','West Coast','South Coast')
  
  COls=c('grey90','grey70','grey40','grey10')
  names(COls)=c('North Coast','Gascoyne Coast','West Coast','South Coast')
  
  tiff(file=paste(hndl.out,"Fig2. Proportion.tiff",sep=''),width=2400,height=2400,
       units="px",res=300,compression="lzw")
  par(mar=c(2.5,8,.5,3),oma=c(.1,.1,.1,.1),mgp=c(.1,.3,0),las=1,xpd=TRUE)
  s=barplot(t(as.matrix(Ktch.by.sp.zn[,-1])),horiz=T,col=COls,
          names.arg = Ktch.by.sp.zn$Common.Name,cex.names=.8,
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
  Factors=MDS%>%select(Bioregion,source,FINYEAR)
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
  legend('bottomright',paste("Stress=",round(NMDS$stress,2)),bty='n')
  legend(x=-.4,y=-.96,c("I-survey","Charter"),pch=c(24,21),pt.cex=1.5,bty='n',horiz = T,cex=1.25)
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
  
  capture.output(SIMPER.bioregion, file = paste(hndl.out,'SIMPER.bioregion.txt',sep=''))
  capture.output(SIMPER.source, file = paste(hndl.out,'SIMPER.source.txt',sep=''))
  
  
  #3. Temporal trends in reconstructed catches    
  #ACA. plot this Rec.ktch.Lower,  Rec.ktch, Rec.ktch.Upper  as done in 'Other Assessment'

}
