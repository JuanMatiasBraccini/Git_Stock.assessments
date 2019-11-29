#-- Script for reconstructing time series of recreational catch of sharks in WA

#notes: This script uses I-Survey (boat-based) point estimates to reconstruct
#       time series considering rec fishing participation rates and WA 
#       population growth
#       Shore-based fishing is added as a proportion of boat-based fishing
#       Charter boat fishing is also added

#MISSING: add 2017-18 I-Survey


library(tidyverse)
library(readxl)
library(lubridate)


Do.recons.rec.fishn.paper="NO"


# 1 -------------------DATA SECTION------------------------------------

# I-Survey
Rec.hndl="C:/Matias/Data/Catch and Effort/Recreational/I.Survey."
Rec.fish.catch.2011.12=read.csv(paste(Rec.hndl,"2011_12.csv",sep=''),stringsAsFactors=F) #Ryan et al 2013    
Rec.fish.catch.2013.14=read.csv(paste(Rec.hndl,"2013_14.csv",sep=''),stringsAsFactors=F) #Ryan et al 2015 
Rec.fish.catch.2015.16=read.csv(paste(Rec.hndl,"2015_16.csv",sep=''),stringsAsFactors=F) #Ryan et al 2015 
Rec.fish.catch=rbind(Rec.fish.catch.2011.12,Rec.fish.catch.2013.14,Rec.fish.catch.2015.16)
#Rec.fish.catch=read.csv(paste(Rec.hndl,"csv",sep=''),stringsAsFactors=F)
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
                                "Zebra Shark"),
                  AVG.wt=c(6.7,6.7,4.2,6.7,6.7,6.7,6.7,6.7,6.7,5,4.2,4.2,4.2,4.2,6.7,6.7,6.7,6.7,
                           rep(6.7,5),4.2,rep(6.7,3),4.2,rep(6.7,7)), 
                  PCM.rec=c(Asmd,Asmd,.1,.7,.7,Asmd,Asmd,Asmd,Asmd,Asmd,Asmd,.05,Asmd,Asmd,Asmd,
                            Asmd,Asmd,Asmd,
                            rep(Asmd,5),.5,rep(Asmd,2),Asmd,rep(Asmd,8)))           
AVG.WT$Common.Name=as.character(AVG.WT$Common.Name)

#reconstruct population fishing
dummy=rbind(cbind(Year=c(1940,1950,1960),Population=c(473300,557100,722100)),WA.population)
mod=loess(Population~Year,data=dummy)
Historic.pop=predict(mod,newdata = data.frame(Year=1941:1970))
WA.population=rbind(cbind(Year=1941:1970,Population=round(Historic.pop)),
                    WA.population)

#fix species names
Rec.fish.catch=Rec.fish.catch%>%
        mutate(Common.Name=ifelse(Bioregion%in%c("Gascoyne","North Coast") & 
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

Whaler.Sharks=Rec.fish.catch%>%
                filter(Common.Name%in%c("Whaler Sharks","Whaler Shark - other"))%>%
                summarise_at(vars(c(Kept.Number,Rel.Number)), sum, na.rm = TRUE)%>%
                data.frame
Whaler.reap=cbind(Whaler.prop,Whaler.Sharks)%>%
                mutate(Kept.Number=Kept.Number*Kept.Number.prop,
                       Rel.Number=Rel.Number*Rel.Number.prop)%>%
        dplyr::select(names(Rec.fish.catch))
     
Rec.fish.catch=Rec.fish.catch%>%
                filter(!Common.Name%in%c("Whaler Sharks","Whaler Shark - other"))
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
fn.rec=function(DAT)
{
        AGG=DAT%>%left_join(AVG.WT,by="Common.Name")%>%
                  mutate(LIVEWT.c=ceiling((Kept.Number+Rel.Number*PCM.rec)*AVG.wt))%>%
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
Rec.ktch=fn.rec(DAT=Rec.fish.catch)  

Part.rate=mean(c(Part.rate.hist,Part.rate.89,Part.rate.00))
Fishing.population=(Part.rate/100)*WA.population$Population
Fishing.population=Fishing.population/Fishing.population[match(2011,WA.population$Year)] #relative to 2011-12
Fishing.population=data.frame(Size=Fishing.population[1:(length(Fishing.population)-1)],
                             FinYear=paste(WA.population$Year[1:(length(WA.population$Year)-1)],"-",
                                        substr(WA.population$Year[2:length(WA.population$Year)],start=3,stop=4),sep=""))

back.fill=function(dat)
{
        id=match(unique(dat$FINYEAR),Fishing.population$FinYear)
        Regns=unique(dat$Bioregion)
        Dummy=vector('list',length(Regns))
        for (d in 1:length(Dummy))
        {
                a=subset(dat,Bioregion==Regns[d])
                Mat=matrix(NA,nrow=length(Fishing.population$FinYear),ncol=2)
                Mat[,2]=mean(a$LIVEWT.c)*Fishing.population$Size
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
for(i in 1:length(Rec.ktch)) Rec.ktch[[i]]=back.fill(dat=Rec.ktch[[i]])
Rec.ktch=do.call(rbind,Rec.ktch)
row.names(Rec.ktch)=NULL




# 7 -------------------Report------------------------------------
if(Do.recons.rec.fishn.paper=="YES")
{
  hndl.out="C:\\Matias\\Analyses\\Reconstruction_catch_recreational\\"
  
  #export weight and PCS table
  write.csv(AVG.WT[order(AVG.WT$Common.Name),],paste(hndl.out,"Appendix.Table.wght.PCS.csv",sep=''),row.names = FALSE)
  
  #Species proportions by zone based on Isurvey and Charter data
  Ktch.by.sp.zn=Rec.fish.catch %>%
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
  
  

  tiff(file=paste(hndl.out,"Fig2. Proportion.tiff",sep=''),width=2400,height=2400,
       units="px",res=300,compression="lzw")
  par(mar=c(2.5,8,1,3),oma=c(.1,.1,.1,.1),mgp=c(.1,.3,0),las=1,xpd=TRUE)
  s=barplot(t(as.matrix(Ktch.by.sp.zn[,-1])),horiz=T,
          names.arg = Ktch.by.sp.zn$Common.Name,cex.names=.8,
          legend.text=T,args.legend=list(x=1,y=nrow(Ktch.by.sp.zn)+8,
                                         bty='n',horiz=T,cex=.9))
  text(.99,s,paste("n=",round(Sp.tot)),pos=4,cex=.8)
  mtext("Proportion",1,1.5,cex=1.5)
  dev.off()
  
  
  #MDS
  library(vegan)
  library(network)
  MDS=Rec.fish.catch %>%
    filter(FINYEAR%in%I.survey.years)%>%
    mutate(Tot=Kept.Number+ Rel.Number)%>%
    group_by(Bioregion, FINYEAR,Common.Name) %>%
    summarise(n = sum(Tot,na.rm=T))%>%
    spread(Common.Name,n,fill=0)
  Factors=MDS%>%select(Bioregion,FINYEAR)
  Factors$FINYEAR=substr(Factors$FINYEAR,1,4)
  col.nms=colnames(MDS)[-(1:2)]
  MDS=as.matrix(MDS[,-(1:2)])
  MDS=sqrt(MDS) #apply sqrt transf
  #MDS=MDS^(1/4) #4th root trans
  colnames(MDS)=col.nms
  NMDS=metaMDS(MDS,k=2,trymax=100)

  plot(NMDS$points[,1],NMDS$points[,2],col='transparent')
  text(NMDS$points[,1],NMDS$points[,2],Factors$FINYEAR,col=as.color(Factors$Bioregion),cex=.8)
}
