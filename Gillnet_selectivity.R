# Script for estimating gear selectivity based on Kirkwood & Walker 1986, extended to 
#       Millar & Holst 1997, and Millar & Fryer 1999

# size classes for selectivity analysis: 5 cm bins (as used in pop dyn model)
# size type for selectivity analysis: total length (in cm)

#assumptions: different mesh sizes set at the same time in same place

#MISSING: 
#         Add Terry's 1970s to 1990s data (if not in TL, then convert to TL; check if 
#                                         southern eagle ray mean size increases with mesh)



library(tidyverse)
library(RODBC)
library(doParallel)
library(Hmisc)
library(magrittr)
library(expandFunctions)
library(stringr)
library(msm)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240) 
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R"))
source(handl_OneDrive("Analyses/Population dynamics/Git_Stock.assessments/NextGeneration.R"))
source(handl_OneDrive("Analyses/Population dynamics/Git_Stock.assessments/SelnCurveDefinitions.R")) #These can be extended by the user


Do.K_W=FALSE   #superceeded by Millar & Holst

# DATA  -------------------------------------------------------------------

#1. WA Fisheries experimental mesh selectivity studies 

  #1.1 1994-1996   Simpfendorfer & Unsworth 1998   (need to physically open this file for R to connect)
channel <- odbcConnectExcel2007("U:/Shark/ExperimentalNet.mdb")
EXP_NET<- sqlFetch(channel,"EXP_NET", colnames = F)
EXPNET_B<- sqlFetch(channel,"EXPNET_B", colnames = F)
close(channel)

  #1.2 2001-2003
channel <- odbcConnectAccess2007("U:/Shark/Sharks v20200323.mdb")  
Boat_bio=sqlFetch(channel, "Boat_bio", colnames = F) 
Boat_hdr=sqlFetch(channel, "Boat_hdr", colnames = F)   
close(channel)


#2. SESSF 2007-2008 experimental mesh selectivity and survey 
channel <- odbcConnectExcel2007(handl_OneDrive("Data/SSF_survey_07_08/SharkSurveyData_30_09_2008.xls"))
F2_Sampling<- sqlFetch(channel,"F2_Sampling", colnames = F)
F1_SamplingTwo<- sqlFetch(channel,"F1_SamplingTwo", colnames = F)
close(channel)


#3. TDGDLF observed catch composition 
#note: this is for fishers using one type of net at a time.
#      Not used as models cannot converge and not used simultaneously.
#      Mean sizes of 6.5 and 7 inch similar, or 6.5 > 7 
LFQ.south=read.csv(handl_OneDrive("Analyses/Selectivity_Gillnet/out.LFQ.south.csv"))


#4. Species names
SP.names=read.csv(handl_OneDrive('Data/Species_names_shark.only.csv'))
SP.codes=read.csv(handl_OneDrive('Data/Species.code.csv'))

#5. Life history
LH=read.csv(handl_OneDrive('Data/Life history parameters/Life_History.csv'))

#6. Rory's data published in McAuley et al 2007
Rory.d=read.csv(handl_OneDrive('Analyses/Selectivity_Gillnet/Rory.sandbar/data.csv')) #size is Fork length


# PARAMETERS  -------------------------------------------------------------------

Min.sample=30  #keep species with at least 30 observations in at least 2 mesh sizes
Min.nets=2
min.obs.per.mesh=Min.sample  #for each kept species, use nets with a minimum # of records

Size.Interval=5 #size intervals for selectivity estimation (in cm)

Min.length=0  #min and max total length considered (in cm)
Max.length=600

do.paper.figures=FALSE
Preliminary=FALSE

Published=capitalize(tolower(c("Sandbar shark","Gummy Shark","Whiskery shark","Dusky shark")))

Published.sel.pars_K.W=data.frame(Species=Published,
                                  Theta1=c(135.5,  184.3, 173.7, 130.1),
                                  Theta2=c(117695, 29739, 26415, 29237))

Add.Rory=FALSE  #add missing sandbar data published in Rory's paper but not in datbase



# Manipulate Rory's  -------------------------------------------------------------------
#not used

# Manipulate 1994-1996  -------------------------------------------------------------------
EXP_NET=EXP_NET[grep("E",EXP_NET$SHEET_NO),]%>%
           mutate(BOAT=VESSEL)%>%
          mutate(SKIPPER=NA,
                 Method='GN',
                 BOTDEPTH=NA,
                 SOAK.TIME=NA,
                 NET_LENGTH=270,   #From Simpfendorfer & Unsworth 1998
                 Lng1=END1_LNG_D+(END1_LNG_M/60),
                 Lng2=END2_LNG_D+(END2_LNG_M/60),
                 Lat1=END1_LAT_D+(END1_LAT_M/60),
                 Lat2=END2_LAT_D+(END2_LAT_M/60))
EXP_NET=EXP_NET%>%
          mutate(MID.LONG=rowMeans(select(EXP_NET,starts_with("Lng"))),
                 MID.LAT=-rowMeans(select(EXP_NET,starts_with("Lat"))))%>%
  dplyr::select(SHEET_NO,DATE,BOAT,SKIPPER,Method,START_SET,END_HAUL,BOTDEPTH,
                MID.LAT,MID.LONG,SOAK.TIME,NET_LENGTH)
    
EXPNET_B=EXPNET_B[grep("E",EXPNET_B$SHEET_NO),]%>%
          left_join(EXP_NET,by='SHEET_NO')%>%
          filter(SPP_CODE<50000)%>%
          mutate(species=SPP_CODE)%>%
          rename(SPECIES1=SPECIES)

Exp.net.94_96=EXPNET_B%>%
                left_join(SP.names,by=c("SPP_CODE" = "SPECIES"))%>%
                filter(!is.na(Name))%>%
                rename(tl=TOT_LENGTH,
                       fl=FORK_LNGTH)
colnames(Exp.net.94_96)=tolower(colnames(Exp.net.94_96))

Exp.net.94_96=Exp.net.94_96%>%
              mutate(experiment='94_96',
                     start_set1=case_when(nchar(start_set)==3~as.numeric(substr(start_set,1,1))+as.numeric(substr(start_set,2,3))/60,
                                         nchar(start_set)==4~as.numeric(substr(start_set,1,2))+as.numeric(substr(start_set,3,4))/60),
                     end_haul1=case_when(nchar(end_haul)==3~as.numeric(substr(end_haul,1,1))+as.numeric(substr(end_haul,2,3))/60,
                                          nchar(end_haul)==4~as.numeric(substr(end_haul,1,2))+as.numeric(substr(end_haul,3,4))/60),
                     soak.time=end_haul1-start_set1,
                     soak.time=ifelse(soak.time<3.9,24+end_haul1-start_set1,soak.time),
                     soak.time=ifelse(is.na(soak.time),mean(soak.time,na.rm=T),soak.time))

# Manipulate 2001-2003  -------------------------------------------------------------------
Boat_hdr=Boat_hdr[grep("E",Boat_hdr$SHEET_NO),]%>%
          dplyr::select(SHEET_NO,DATE,BOAT,SKIPPER,Method,START_SET,END_HAUL,BOTDEPTH,
                        'MID LAT','MID LONG','SOAK TIME',MESH_SIZE,MESH_DROP,NET_LENGTH)%>%
          rename(MID.LAT='MID LAT',MID.LONG='MID LONG',SOAK.TIME='SOAK TIME')%>%
          filter(Method=='GN')
Boat_bio=Boat_bio[grep("E",Boat_bio$SHEET_NO),]%>%
          dplyr::select(SHEET_NO,SPECIES,TL,FL,PL,SEX)
Exp.net.01_03=left_join(Boat_bio,Boat_hdr,by="SHEET_NO")%>%
          mutate(SEX=ifelse(SEX=="m","M",ifelse(SEX=="f","F",SEX)),
                 MESH_SIZE=ifelse(MESH_SIZE=="10\"","10",
                           ifelse(MESH_SIZE=="6\"","6",
                           ifelse(MESH_SIZE=="5\r\n5","5",
                           ifelse(MESH_SIZE=="7\"","7",
                           ifelse(MESH_SIZE=="5\"","5",
                           ifelse(MESH_SIZE=="4\"","4",
                           ifelse(MESH_SIZE=="8\"","8",
                           MESH_SIZE))))))),
                 MESH_SIZE=as.numeric(MESH_SIZE))%>%
          dplyr::select(-PL)%>%
  left_join(SP.codes,by=c("SPECIES" = "Species"))%>%
  rename(Name=COMMON_NAME)
colnames(Exp.net.01_03)=tolower(colnames(Exp.net.01_03))

Exp.net.01_03=Exp.net.01_03%>%
                mutate(experiment='01_03',
                       net_length=case_when(mesh_size==4~91,
                                            mesh_size==5~mean(c(91,113,108,105,103)),
                                            mesh_size==5.5~113,
                                            mesh_size==6~108,
                                            mesh_size==7~105,
                                            mesh_size==8~mean(c(91,113,108,105,103)),
                                            mesh_size==8.8~105,
                                            mesh_size==10~103,
                                            TRUE~NA_real_),
                       soak.time=ifelse(is.na(soak.time),mean(soak.time,na.rm=T),soak.time))

This.col=c('sheet_no','date','experiment','mid.lat','mid.long','mesh_size',
           'mesh_drop','name','tl','fl','sex','net_length','soak.time')

# Combine 1994-1996 & 2001-2003  -------------------------------------------------------------------
Exp.net.WA=rbind(Exp.net.94_96[,match(This.col,names(Exp.net.94_96))],
                 Exp.net.01_03[,match(This.col,names(Exp.net.01_03))])%>%
          mutate(name=ifelse(name=="Angel Shark (general)","Angel Shark",
                      ifelse(name=="Eagle ray","Eagle Ray",
                      ifelse(name=="Gummy shark","Gummy Shark",
                      ifelse(name%in%c("Port Jackson","Port Jackson shark"),"PortJackson shark",
                      ifelse(name=="Big eye sixgill shark","Sixgill shark",
                      ifelse(name%in%c("Wobbegong (general)",'Spotted Wobbegong','Western Wobbegong'),"Wobbegongs",name)))))))

    #export for length-length conversion study
write.csv(Exp.net.WA%>%filter(!is.na(fl) & !is.na(tl))%>%
            mutate(name=capitalize(tolower(name))),handl_OneDrive('Analyses/Length conversions/Interdorsal length conversions/Exp.net.WA.csv'),row.names = F)

TAB=table(Exp.net.WA$name,Exp.net.WA$mesh_size)
TAB[TAB<2]=0
TAB[TAB>=2]=1

Exp.net.WA=Exp.net.WA%>%
  filter(name%in%names(which(rowSums(TAB)>=2)))%>%
  mutate(tl=ifelse(tl>500,tl/10,tl))


if(Preliminary) ggplot(Exp.net.WA,aes(tl,fl,shape=name, colour=name, fill=name))+
                geom_point() + 
                geom_smooth(method = "lm", fill = NA)+
                facet_wrap(vars(name), scales = "free")


#Convert FL to TL for records with no TL
TL_FL=data.frame(name=c('Angel Shark','Dusky shark','Gummy Shark','Pencil shark',
           'PortJackson shark','Sandbar shark','Smooth hammerhead','Spurdogs',
           'Whiskery shark','Tiger shark','Sliteye shark'))%>%
          mutate(get.this=str_remove(tolower(name), " shark"))

LH=LH%>%mutate(SNAME=str_remove(tolower(SNAME), "shark, "),
              SNAME=ifelse(SNAME=="smooth hh","smooth hammerhead",
                           ifelse(SNAME=="spurdog","spurdogs",SNAME)))
TL_FL=TL_FL%>%
        left_join(LH%>%dplyr::select(SNAME,a_FL.to.TL,b_FL.to.TL),by=c('get.this'='SNAME'))%>%
        rename(intercept=b_FL.to.TL,
               slope=a_FL.to.TL)
Exp.net.WA=Exp.net.WA%>%
           left_join(TL_FL,by='name')%>%
            mutate(tl=ifelse(is.na(tl),intercept+fl*slope,tl),
                   Length=tl)%>%   #Length in cm; use Total length for selectivity estimation
            filter(!is.na(Length))%>%
            filter(!name=='Sliteye shark')#too few observations for sliteye

#Size frequency   
if(Preliminary)
{
  dummy=Exp.net.WA%>%mutate(mesh_size=factor(mesh_size))
  ggplot(dummy,aes(tl,fill=mesh_size))+
    geom_histogram(color="#e9ecef",alpha=0.6, binwidth = 5) +
    labs(fill="Mesh size")+
    xlab("Total length (cm)")+
    facet_wrap(vars(name), scales = "free") 
  ggsave(handl_OneDrive('Analyses/Selectivity_Gillnet/Size.frequency_experimental.WA.tiff'), width = 10,height = 8, dpi = 300, compression = "lzw")
}

# Manipulate SSF  -------------------------------------------------------------------
F2_Sampling=F2_Sampling%>%
                filter(Csiro>37000002 & Csiro<37110000 & !is.na(Mesh1) & !is.na(Length))%>%
                mutate(Length=Length/10,          #length in cm      
                       Length=ifelse(Length>350,Length/10,Length),
                       Mesh1=ifelse(Mesh1%in%c('C6.00', 'C6.01', 'C6.02', 'C6.03', 'C6.04'),'C6',
                                    ifelse(Mesh1%in%c('C6.51', 'C6.52', 'C6.53', 'C6.54'),'C6.5',Mesh1)),
                       Mesh.size=2.54*as.numeric(substr(Mesh1,2,10)))%>%   #inches to cm
                dplyr::select(Cruise,Station,Species,Csiro,Sex,Mesh1,Mesh.size,Length,LengthType)

if(Preliminary) ggplot(F2_Sampling, aes(x = Length)) +
                    geom_histogram(color = "grey30", fill ="salmon",binwidth=10) +
                    facet_grid(Species~Mesh.size, scales = "free")

#convert FL or DW to TL (Southern Eagle Ray left as DW as tails can be chopped)
F2_Sampling=F2_Sampling%>%
  mutate(Length=case_when(Csiro==37005002 & LengthType=='FL' ~ 1.267105*Length,
                          Csiro==37007001 & LengthType=='FL' ~ 1.105265*Length+0.4579532,
                          Csiro==37010001 & LengthType=='FL' ~ 1.07689*Length+1.7101,
                          Csiro==37017001 & LengthType=='FL' ~ 1.09528*Length+3.8499,
                          Csiro==37020008 & LengthType=='FL' ~ 1.138*Length+5.736,
                          Csiro==37043001 & LengthType=='FL' ~ 1.012*Length+13.42,
                          Csiro==37038004 & LengthType=='DW' ~ (Length-17.937)/0.589,
                          Csiro==37038007 & LengthType=='DW' ~ (Length-17.937)/0.589,
                          TRUE ~ Length))%>%
  filter(!(Csiro==37039001 & LengthType %in%c('FL','TL')))

F2_Sampling=F2_Sampling%>%  #remove too small records
   filter(!(Length<30 & LengthType%in%c('TL','tl')))%>%
  filter(!(Length<5 & LengthType%in%c('DW')))%>%
  dplyr::select(-LengthType)
  
  #add soak time
Mn_times=F1_SamplingTwo%>%
              group_by(Cruise,Station)%>%
              summarise(TimeSet_mean=mean(TimeSet,na.rm=T),
                        TimeUp_mean=mean(TimeUp,na.rm=T))%>%
              data.frame
Mn_times2=F1_SamplingTwo%>%
  dplyr::select(Cruise,Station,Mesh,TimeSet,TimeUp)%>%
  left_join(Mn_times,by=c('Cruise','Station'))%>%
  mutate(TimeSet=format(TimeSet, format = '%H:%M'),
         TimeUp=format(TimeUp, format = '%H:%M'),
         TimeSet_mean=format(TimeSet_mean, format = '%H:%M'),
         TimeUp_mean=format(TimeUp_mean, format = '%H:%M'),
         TimeSet=ifelse(is.na(TimeSet),TimeSet_mean,TimeSet),
         TimeUp=ifelse(is.na(TimeUp),TimeUp_mean,TimeUp),
         Haul=(as.numeric(substr(TimeUp,1,2))+as.numeric(substr(TimeUp,4,5))/60),
         Set=(as.numeric(substr(TimeSet,1,2))+as.numeric(substr(TimeSet,4,5))/60),
         soak.time=Haul-Set,
         soak.time=ifelse(soak.time<1,24+soak.time,soak.time))%>%
  dplyr::select(Cruise,Station,Mesh,soak.time)%>%
  mutate(Mesh=ifelse(Mesh%in%c('C6.00', 'C6.01', 'C6.02', 'C6.03', 'C6.04'),'C6',
              ifelse(Mesh%in%c('C6.51', 'C6.52', 'C6.53', 'C6.54'),'C6.5',
                     Mesh)))%>%
  distinct(Cruise,Station,Mesh,.keep_all=T)
F2_Sampling=F2_Sampling%>%
              left_join(Mn_times2,by=c('Cruise','Station','Mesh1'='Mesh'))%>%
               mutate(net_length=500,  #from Braccini et al 2009
                      Species=ifelse(Species=='Bronze Whaler',"Copper shark",Species))

#Size frequency
if(Preliminary)
{
  dummy=F2_Sampling%>%mutate(mesh_size=factor(Mesh.size))
  ggplot(dummy,aes(Length,fill=mesh_size))+
    geom_histogram(color="#e9ecef",alpha=0.6, binwidth = 5) +
    labs(fill="Mesh size")+
    xlab("Size (cm)") +
    facet_wrap(vars(Species), scales = "free") 
  ggsave(handl_OneDrive('Analyses/Selectivity_Gillnet/Size.frequency_SSF.tiff'), width = 12,height = 8, dpi = 300, compression = "lzw")
}

# Manipulate TDGLDF  -------------------------------------------------------------------
#note: not used because observations for different mesh sizes not collected @ same time
do.LFQ.south=FALSE
if(do.LFQ.south)
{
  LFQ.south=LFQ.south%>%
    mutate(Data.set="Obs.TDGDLF",
           Mesh.size=2.54*MESH_SIZE,
           Species=capitalize(COMMON_NAME),
           Species=ifelse(Species=="Sawsharks","Common sawshark",
                          ifelse(Species=="Wobbegong (general)","Wobbegong",
                                 Species)))%>%
    filter(!Species=="Wobbegong")  #remove Wobbies because size measurements are uncertain
  
  #Convert FL to TL
  TL_FL_LFQ.south=TL_FL%>%
    filter(name%in%c("Smooth hammerhead","Spurdogs","Tiger shark"))%>%
    dplyr::select(-get.this)
  
  add1=TL_FL_LFQ.south%>%filter(name=="Smooth hammerhead")
  add1=rbind(add1,add1)%>%mutate(name=c("Great hammerhead","Scalloped hammerhead"))  
  
  TL_FL_LFQ.south=rbind(TL_FL_LFQ.south,
                        add1,
                        data.frame(name="Copper shark",intercept=6.972,slope=1.214),
                        data.frame(name="Common sawshark",intercept=0,slope=1.1),
                        data.frame(name="Grey nurse shark",intercept=0,slope=1.3),
                        data.frame(name="Milk shark",intercept=5.24,slope=1.1162),
                        data.frame(name="Shortfin mako",intercept=1.71,slope=1.07689),
                        data.frame(name="Spinner shark",intercept=3.44,slope= 1.1814))  
  
  LFQ.south=LFQ.south%>%
    left_join(TL_FL_LFQ.south,by=c(Species='name'))%>%
    mutate(tl=intercept+FL*slope,
           Length=tl)%>%   #Length in cm; use Total length
    filter(!is.na(Length))%>%
    dplyr::select(Species,Mesh.size,Length,Data.set)%>%
    filter(!Species=='Spurdogs')  #mixed bag of species
  
  #Size frequency
  if(Preliminary)
  {
    dummy=LFQ.south%>%mutate(Mesh_size=as.factor(Mesh.size))
    ggplot(dummy,aes(Length,fill=Mesh_size))+
      geom_histogram(color="#e9ecef",alpha=0.6, binwidth = 5) +
      labs(fill="Mesh size")+
      xlab("Total length (cm)") +
      facet_wrap(vars(Species), scales = "free") 
    ggsave(handl_OneDrive('Analyses/Selectivity_Gillnet/Size.frequency_TDGDLF_observed.tiff'), width = 10,height = 8, dpi = 300, compression = "lzw")
  }
  
}

# Combine all data sets  -------------------------------------------------------------------
   
Exp.net.WA=Exp.net.WA%>%
  rename(Species=name)%>%
  mutate(Mesh.size=2.54*mesh_size,
         Data.set="WA")
F2_Sampling=F2_Sampling%>%
        mutate(Data.set="SSF")%>%
        filter(!Mesh1%in%c("C6.5","C6"))


Combined=rbind(F2_Sampling%>%
                 dplyr::select(Species,Mesh.size,Length,Data.set,net_length,soak.time),
               Exp.net.WA%>%
                 dplyr::select(Species,Mesh.size,Length,Data.set,net_length,soak.time))%>%
                  mutate(Species=tolower(Species),
                         Species=ifelse(Species=='portjackson shark','port jackson shark',
                                 ifelse(Species%in%c('wobbegong','banded wobbegong','spotted wobbegong'),
                                        'wobbegongs',Species)),
                         Species=capitalize(Species),
                         Species=ifelse(Species=='Port jackson shark','Port Jackson shark',Species),
                         Mesh.size=round(Mesh.size,1),
                         Species=ifelse(Species=='Angel shark','Australian angelshark',Species))%>%
  filter(!is.na(Mesh.size))%>%
  filter(!Species%in%c('Spurdogs'))%>% #remove Spurdogs from WA, could be several species
  filter(!(Species=='Eagle ray' & Data.set=='WA'))  #remove because don't know if measured TL or DW
SP.names=SP.names%>%
  mutate(Name=ifelse(Name=='PortJackson shark','Port Jackson shark',
              ifelse(Name=='Angel Shark','Australian angelshark',
              ifelse(Name=='Sevengill shark','Broadnose shark',
              ifelse(Name=='Sixgill shark','Bluntnose sixgill shark',
              ifelse(Name=='Eagle ray','Southern eagle ray',
              Name))))))

Families=SP.names%>%  
  dplyr::select(Name,Family)%>%
  rename(Species=Name)
  
Combined=Combined%>%left_join(Families,"Species")
Combined.family=Combined%>%
                rename(Species.original=Species,
                       Species=Family)  
Tab.sp.fam=table(Combined$Species,Combined$Family)
Tab.sp.fam[Tab.sp.fam>0]=1
Tab.sp.fam=colSums(Tab.sp.fam)

# Add missing Sandbar records occurring in McAuley et al 2007 but not in Data base  ----------------------------
if(Add.Rory)
{
  Rory.tot=colSums(Rory.d[,-1])
  dum=Combined%>%
    filter(Species=="Sandbar shark")
  Tab.san=dum%>%
          group_by(Mesh.size)%>%
          tally()%>%
          spread(Mesh.size,n)%>%
          data.frame
  id1=Tab.san[match(names(Rory.tot),colnames(Tab.san))]-Rory.tot
  id1=abs(id1[which(id1<0)])
  dum=dum%>%
    filter(Mesh.size%in%as.numeric(substr(names(id1),2,5)))%>%
    mutate(size.class=Size.Interval*floor(Length/Size.Interval)+Size.Interval/2,
           dummy=paste(Mesh.size,size.class))
  Rory.d1=Rory.d%>%
    gather(Mesh,n,-Size.class)%>%
    mutate(Mesh=substr(Mesh,2,5))%>%
    mutate(dummy=paste(Mesh,Size.class))%>%
    filter(n>0 & Mesh%in%as.numeric(substr(names(id1),2,5)))
  
  Rory.d1=Rory.d1[rep(seq(nrow(Rory.d1)), Rory.d1$n),]
  drops=match(dum$dummy,Rory.d1$dummy)
  drops=drops[!duplicated(drops)]
  drops=drops[!is.na(drops)]
  Rory.d1=Rory.d1[-drops,]
  Rory.d1=Rory.d1%>%
            rename(Length=Size.class)%>%
            mutate(Species='Sandbar shark',
                   Data.set="WA",
                   Mesh.size=as.numeric(Mesh),
                   Family='Carcharhinidae')%>%
    dplyr::select(names(Combined))
  Combined=rbind(Combined,Rory.d1)
}

# Tables  -------------------------------------------------------------------

#Table 1. Numbers by species and mesh
Table1=Combined%>%
  group_by(Mesh.size,Species)%>%
  tally()%>%
  spread(Mesh.size,n)%>%
  data.frame
colnames(Table1)[-1]=substr(colnames(Table1)[-1],2,10)
Table1=Table1%>%
  left_join(SP.names%>%dplyr::select(-SPECIES),by=c('Species'='Name'))%>%
  arrange(Family,Species)

# Select species  -------------------------------------------------------------------
TAB=table(Combined$Species,Combined$Mesh.size)
TAB[TAB<Min.sample]=0
TAB[TAB>=Min.sample]=1
Combined=Combined%>%
  filter(Species%in%names(which(rowSums(TAB)>=Min.nets)) & Length<=Max.length)%>%
  filter(!Species=='Wobbegongs')%>%  #unidentified Wobbegongs used in Family analysis
  filter(!Species=="Southern eagle ray")   # remove southern eagle ray because mean size decreases with mesh size

# Select families  -------------------------------------------------------------------
TAB=table(Combined.family$Species,Combined.family$Mesh.size)
TAB[TAB<20]=0
TAB[TAB>=20]=1
Combined.family=Combined.family%>%
  filter(Species%in%names(which(rowSums(TAB)>=Min.nets)) & Length<=Max.length)

Combined.family=Combined.family%>%
              filter(!Species%in%c("Scyliorhinidae",names(Tab.sp.fam[Tab.sp.fam==1])))%>%    #remove family with single species
              filter(!Species=="Squatinidae")  #remove Squatinidae because mean size decreases with mesh size


#for each selected species, remove meshes with few observations
  #species
n.sp=sort(unique(Combined$Species))
for(s in 1:length(n.sp))
{
  d=Combined%>%filter(Species==n.sp[s])
  Combined=Combined%>%filter(!Species==n.sp[s])
  
  id=table(d$Mesh.size)
  Drop=as.numeric(gsub("[^0-9.]", "",  names(which(id<min.obs.per.mesh))))
  if(length(Drop)>0)
  {
    d=d%>%filter(!Mesh.size%in%Drop)
    id=table(d$Mesh.size)
    if(length(id)==1) d=NULL
  }
  
  Combined=rbind(Combined,d)
  
}
n.sp=sort(unique(Combined$Species))

  #family
n.sp.family=sort(unique(Combined.family$Species))
for(s in 1:length(n.sp.family))
{
  d=Combined.family%>%filter(Species==n.sp.family[s])
  Combined.family=Combined.family%>%filter(!Species==n.sp.family[s])
  
  id=table(d$Mesh.size)
  Drop=as.numeric(gsub("[^0-9.]", "",  names(which(id<20))))
  if(length(Drop)>0)
  {
    d=d%>%filter(!Mesh.size%in%Drop)
    id=table(d$Mesh.size)
    if(length(id)==1) d=NULL
  }
  
  Combined.family=rbind(Combined.family,d)
  
}
n.sp.family=sort(unique(Combined.family$Species))


#Remove mesh sizes with higher mean size than next largest mesh
  #species
Tab.mean.size=Combined%>%group_by(Species,Mesh.size)%>%summarise(Mean=mean(Length))%>%spread(Mesh.size,Mean)
Tab.n=Combined%>%group_by(Species,Mesh.size)%>%tally()%>%spread(Mesh.size,n)

drop.one.mesh=data.frame(Species=c("Port Jackson shark",
                                   rep("Sandbar shark",4),
                                   "Whiskery shark",
                                   "Draughtboard shark",
                                   rep("Smooth hammerhead",2),
                                   "Spikey dogfish"),
                         Mesh.size=c(10.2,
                                     10.2,14,21.6,22.4,
                                     16.5,
                                     10.2,
                                     16.5,21.6,
                                     17.8),
                         Drop="YES")
Combined=left_join(Combined,drop.one.mesh,by=c('Species','Mesh.size'))%>%
         filter(is.na(Drop))%>%dplyr::select(-Drop)

  #family
drop.one.mesh=data.frame(Species=c(rep("Carcharhinidae",4),
                                   "Orectolobidae",
                                   "Squalidae",
                                   "Scyliorhinidae",
                                   "Triakidae"),
                         Mesh.size=c(10.2,16.5,21.6,22.4,
                                     15.2,
                                     20.3,
                                     10.2,
                                     21.6),
                         Drop="YES")
Combined.family=left_join(Combined.family,drop.one.mesh,by=c('Species','Mesh.size'))%>%
        filter(is.na(Drop))%>%dplyr::select(-Drop)

# Add scientific names to LH-------------------------------------------------------------------------
LH=LH%>%left_join(SP.names,by='SPECIES')   

# Get total length at age -------------------------------------------------------------------------
len.at.age=function(Lo,Linf,k,Age.max)
{
  VonB=data.frame(Age=0:Age.max)
  VonB$TL=(Lo+(Linf-Lo)*(1-exp(-k*VonB$Age))) 
  return(VonB)
}
LatAge=vector('list',length(n.sp))
names(LatAge)=n.sp
LatAge.family=vector('list',length(n.sp.family))
names(LatAge.family)=n.sp.family

for(s in 1:length(n.sp)) 
{
  ii=n.sp[s]
  ii=ifelse(ii=="Southern sawshark","Common sawshark",
     ifelse(ii=="Spikey dogfish","Spurdogs",ii))
  this.par=LH%>%filter(Name==ii) 
  if(!is.na(this.par$K))
  {
    if(is.na(this.par$Max_Age_max)) this.par$Max_Age_max=this.par$Max_Age*1.3
    LatAge[[s]]=with(this.par,len.at.age(Lo=LF_o*a_FL.to.TL+b_FL.to.TL,Linf=FL_inf*a_FL.to.TL+b_FL.to.TL,k=K,Age.max=Max_Age_max))
  }
  if(is.na(this.par$K) & !is.na(this.par$LF_o)) LatAge[[s]]=data.frame(Age=0,TL=this.par$LF_o)
}
for(s in 1:length(n.sp.family)) 
{
  ii=n.sp.family[s]
  this.par=LH%>%filter(Family==ii)%>%
    mutate(Max_Age_max=ifelse(is.na(Max_Age_max),Max_Age*1.3,Max_Age_max))%>%
    group_by(Family)%>%
    summarise(LF_o=mean(LF_o,na.rm=T),
              FL_inf=mean(FL_inf,na.rm=T),
              K=mean(K,na.rm=T),
              Max_Age_max=mean(Max_Age_max,na.rm=T),
              a_FL.to.TL=mean(a_FL.to.TL,na.rm=T),
              b_FL.to.TL=mean(b_FL.to.TL,na.rm=T))
  if(!is.na(this.par$K))
  {
    LatAge.family[[s]]=with(this.par,len.at.age(Lo=LF_o*a_FL.to.TL+b_FL.to.TL,Linf=FL_inf*a_FL.to.TL+b_FL.to.TL,k=K,Age.max=Max_Age_max))
  }
}  

# Remove species with no growth information -------------------------------------------------------------------------
n.sp.no.growth=names(LatAge[sapply(LatAge, is.null)])
if(length(n.sp.no.growth)>0)LatAge=LatAge[-match(n.sp.no.growth,names(LatAge))]


# Estimate selectivity parameters -------------------------------------------------------------------------
#notes: shark size (total length) in cm (in mm for Kirkwood & Walker)
#       size intervals: 5 cm
#       mesh size in cm  (in inches for Kirkwood & Walker)
#       only estimate for species with no selectivty


  #3.1 Millar & Holst 1997 
Fitfunction='gillnetfit'
#Fitfunction='NetFit'  #not used because it doesn't have gamma implemented
PlotLens=seq(Min.length+Size.Interval/2,Max.length-Size.Interval/2,by=Size.Interval)  #midpoints, 50 mm intervals  
Rtype=c("norm.loc","norm.sca","gamma","lognorm")   #consider this selection curves: normal fixed spread, 
                                                   #       normal spread proportional to mesh size,
                                                   #       gamma, lognormal
#Fit functions
Millar.Holst=function(d,size.int,length.at.age,weight.by.effort=NULL)
{
  #Create size bins
  d=d%>%mutate(Size.class=size.int*floor(Length/size.int)+size.int/2)
  
  #Tabulate observations by mid size class and mesh size
  if(is.null(weight.by.effort))
  {
    tab=d%>%
      group_by(Mesh.size,Size.class)%>%
      summarise(n=n())%>%
      spread(Mesh.size,n,fill=0)%>%
      data.frame

  }else
  {
    tab=d%>%     
      mutate(cpue=1000*1/(net_length*soak.time))%>%  #cpue in km gn hours
      group_by(Mesh.size,Size.class)%>%
      summarise(n=sum(cpue))%>%
      spread(Mesh.size,n,fill=0)%>%
      data.frame
  }

  Meshsize=as.numeric(substr(names(tab)[-1],2,10))
  
  #Fit SELECT model
  
    #1. Equal fishing power
  pwr=rep(1,length(Meshsize))
  Equal.power=vector('list',length(Rtype))
  names(Equal.power)=Rtype
  if(Fitfunction=='gillnetfit')
  {
    for(f in 1:length(Equal.power))
    {
      Equal.power[[f]]=gillnetfit(data=as.matrix(tab),
                                  meshsizes=Meshsize,
                                  type=Rtype[f],
                                  rel=pwr,
                                  plots=c(F,F),
                                  plotlens=PlotLens,
                                  plotlens_age=length.at.age,
                                  details=T)
      Equal.power[[f]]$Warnings=warnings()
      #reset.warnings()
    }
  }
  if(Fitfunction=='NetFit')
  {
    Init.par=c(mean(d$Length),sd(d$Length))
    for(f in 1:length(Equal.power))
    {
      Equal.power[[f]]=NetFit(Data=tab,
                              Meshsize=Meshsize,
                              x0=Init.par,
                              rtype=Rtype[f],
                              rel.power=pwr)
      Equal.power[[f]]$Warnings=warnings()
      reset.warnings()
    }
  }
  
    #2. fishing power proportional to mesh size
  pwr=Meshsize
  Prop.power=vector('list',length(Rtype))
  names(Prop.power)=Rtype
  if(Fitfunction=='gillnetfit')
  {
    for(f in 1:length(Prop.power))
    {
      Prop.power[[f]]=gillnetfit(data=as.matrix(tab),
                                  meshsizes=Meshsize,
                                  type=Rtype[f],
                                  rel=pwr,
                                  plots=c(F,F),
                                  plotlens=PlotLens,
                                  plotlens_age=length.at.age,
                                  details=T)
      Prop.power[[f]]$Warnings=warnings()
      #reset.warnings()
    }
  }
  if(Fitfunction=='NetFit')
  {
    Init.par=c(mean(d$Length),sd(d$Length))
    for(f in 1:length(Prop.power))
    {
      Prop.power[[f]]=NetFit(Data=tab,
                              Meshsize=Meshsize,
                              x0=Init.par,
                              rtype=Rtype[f],
                              rel.power=pwr)
      Prop.power[[f]]$Warnings=warnings()
      reset.warnings()
    }
  }

  return(list(tab=tab,Equal.power=Equal.power,Prop.power=Prop.power))
}
  #species 
Fit.M_H=vector('list',length(n.sp))
names(Fit.M_H)=n.sp
for(s in 1:length(n.sp))  
{
  if(!n.sp[s]%in%Published.sel.pars_K.W$Species)
  {
    d=Combined%>%filter(Species==n.sp[s])
    if(n.sp[s]=='Smooth hammerhead') #combined data for smooth HH yields too broad sel. profiles way beyond the observed data
    {
      d=Exp.net.WA%>%
        filter(Species=='Smooth hammerhead' & !Mesh.size%in%c(4,10))
    }
       
    Fit.M_H[[s]]=Millar.Holst(d=d,
                              size.int=Size.Interval,
                              length.at.age=LatAge[[s]]$TL,
                              weight.by.effort='Yes')
  }

}

  #family    
Fit.M_H.family=vector('list',length(n.sp.family))
names(Fit.M_H.family)=n.sp.family
for(s in 1:length(n.sp.family))  
{
  Fit.M_H.family[[s]]=Millar.Holst(d=Combined.family%>%filter(Species==n.sp.family[s]),
                                   size.int=Size.Interval,
                                   length.at.age=LatAge.family[[s]]$TL,
                                   weight.by.effort='Yes')
}

#3.2 Kirkwood & Walker
pred.Kirkwood.Walker=function(theta,pred.len,Mesh)
{
  Theta1=exp(theta[1])
  Theta2=exp(theta[2])
  if(!16.5%in%Mesh) Mesh=c(Mesh,16.5)
  if(!17.8%in%Mesh) Mesh=c(Mesh,17.8)
  Mesh=sort(Mesh)
  
  Mesh=round(Mesh*0.393701,1)  #convert cm to inches
  pred.len=pred.len*10         #length in mm
  
  d1=data.frame(Size.class=rep(pred.len,times=length(Mesh)),
                Mesh.size=rep(Mesh,each=length(pred.len)))%>%
    mutate(alpha.beta=Theta1*Mesh.size,
           beta=-0.5*((alpha.beta)-((alpha.beta*alpha.beta+4*Theta2)^0.5)),
           alpha=alpha.beta/beta,
           Rel.sel=((Size.class/(alpha*beta))^alpha)*(exp(alpha-(Size.class/beta))))%>%
    dplyr::select(-c('alpha.beta','beta','alpha'))%>%
    spread(Mesh.size,Rel.sel)
  return(d1)
}

if(Do.K_W)   #this needs TL in mm and mesh size in inches
{

  Selectivty.Kirkwood.Walker=function(d,size.int,theta)
  {
    #size in mm and mesh in inches
    size.int=size.int*10
    d=d%>%
      mutate(Mesh.size=round(Mesh.size*0.393701,1),
            Length=Length*10)
    #Create size bins
    d=d%>%mutate(Size.class=size.int*floor(Length/size.int)+size.int/2)
    
    #Tabulate observations by mid size class and mesh size
    tab=d%>%
      group_by(Mesh.size,Size.class)%>%
      summarise(n=n())%>%
      spread(Mesh.size,n,fill=0)%>%
      data.frame
    row.names(tab)=tab$Size.class
    tab=tab[,-1]
    
    
    #Calculate relative selectivity
    Theta1=exp(theta[1])
    Theta2=exp(theta[2])
    d=d%>%
      mutate(alpha.beta=Theta1*Mesh.size,
             beta=-0.5*((alpha.beta)-((alpha.beta*alpha.beta+4*Theta2)^0.5)),
             alpha=alpha.beta/beta,
             Rel.sel=((Size.class/(alpha*beta))^alpha)*(exp(alpha-(Size.class/beta))))
    
    
    #Log likelihood
    S.ij=d%>%
      distinct(Size.class,Mesh.size,.keep_all = T)%>%
      dplyr::select(Size.class,Mesh.size,Rel.sel)%>%
      spread(Mesh.size,Rel.sel,fill = 0)
    row.names(S.ij)=S.ij$Size.class
    S.ij=S.ij[,-1]
    
    mu.j=rowSums(tab)/rowSums(S.ij)
    mu.j=sapply(mu.j,function(x) max(x,0.1))
    mu.j.prop=mu.j/sum(mu.j)
    
    
    #predicted numbers
    sum.n=colSums(tab)
    NN=ncol(S.ij)
    tab.pred=(S.ij*matrix(rep(mu.j.prop,NN),ncol=NN)*matrix(rep(sum.n,each=nrow(S.ij)),ncol=NN))/
      (matrix(rep(colSums(S.ij*matrix(rep(mu.j.prop,NN))),each=nrow(S.ij)),ncol=NN))
    
    
    #Gamma log like
    negLL=min(-sum(tab*(log(mu.j*S.ij))-(mu.j*S.ij),na.rm=T),1e100)
    
    return(list(negLL=negLL,d=d,observed=tab,predicted=tab.pred))
    
  }
  theta.list=vector('list',length(n.sp))
  names(theta.list)=n.sp
  Fit.K_W=Pred.sel.K_W=theta.list
  Pred.sel.K_W_len.at.age=Pred.sel.K_W
  Pred.sel.K_W.family=vector('list',length(n.sp.family))
  names(Pred.sel.K_W.family)=n.sp.family
  Pred.sel.K_W.family_len.at.age=Fit.K_W.family=Pred.sel.K_W.family
  
  # initial parameter values
  theta.list$`Gummy shark`=c(Theta1=log(180),Theta2=log(29000))
  for(s in 1:length(theta.list)) theta.list[[s]]=jitter(theta.list$`Gummy shark`,factor=.1)
  theta.list$`Smooth hammerhead`=c(Theta1=4.5,Theta2=14.5) 
  # fit model and make predictions
  #1. Species
  for(s in 1:length(n.sp))
  {
    #. objfun to minimize
    theta=theta.list[[s]]
    D=Combined%>%filter(Species==n.sp[s])
    if(n.sp[s]=="Smooth hammerhead") D=Combined%>%filter(Species==n.sp[s])%>%filter(!Mesh.size==15.2) #not converging with this mesh size
    fn_ob=function(theta)Selectivty.Kirkwood.Walker(d=D,size.int=Size.Interval,theta)$negLL
    
    #. fit model
    Fit.K_W[[s]]=nlminb(theta, fn_ob, gradient = NULL)
    
    #. predict selectivity   
    #PlotLens
    Pred.sel.K_W[[s]]=pred.Kirkwood.Walker(theta=Fit.K_W[[s]]$par,
                                           pred.len=PlotLens,
                                           Mesh=Combined%>%
                                             filter(Species==n.sp[s])%>%
                                             distinct(Mesh.size)%>%
                                             pull(Mesh.size))
    #Lengts at age
    Pred.sel.K_W_len.at.age[[s]]=pred.Kirkwood.Walker(theta=Fit.K_W[[s]]$par,
                                                      pred.len=LatAge[[s]]$TL,
                                                      Mesh=Combined%>%
                                                        filter(Species==n.sp[s])%>%
                                                        distinct(Mesh.size)%>%
                                                        pull(Mesh.size))
    rm(D)
  }
  
  #2. Family
  for(s in 1:length(n.sp.family)) 
  {
    #. objfun to minimize
    theta=theta.list[[s]]
    fn_ob=function(theta)Selectivty.Kirkwood.Walker(d=Combined.family%>%filter(Species==n.sp.family[s]),
                                                    size.int=Size.Interval,
                                                    theta)$negLL
    
    #. fit model
    Fit.K_W.family[[s]]=nlminb(theta.list[[s]], fn_ob, gradient = NULL)
    
    #. predict selectivity   
    #PlotLens
    Pred.sel.K_W.family[[s]]=pred.Kirkwood.Walker(theta=Fit.K_W.family[[s]]$par,
                                                  pred.len=PlotLens,
                                                  Mesh=Combined%>%
                                                    filter(Species==n.sp[s])%>%
                                                    distinct(Mesh.size)%>%
                                                    pull(Mesh.size))
    #Lengts at age
    Pred.sel.K_W.family_len.at.age[[s]]=pred.Kirkwood.Walker(theta=Fit.K_W.family[[s]]$par,
                                                             pred.len=LatAge.family[[s]]$TL,
                                                             Mesh=Combined%>%
                                                               filter(Species==n.sp[s])%>%
                                                               distinct(Mesh.size)%>%
                                                               pull(Mesh.size))
  }
  
  
  # Calculate confidence intervals thru bootstrapping 
  do.boot=FALSE
  if(do.boot)
  {
    n.boot=1:1000
    cl <- makeCluster(detectCores()-1)
    registerDoParallel(cl)
    system.time({
      Fit.K_W.CI=foreach(s=1:length(n.sp),.packages=c('tidyverse','doParallel','Biobase')) %dopar%
        {
          boot=foreach(n=n.boot,.packages=c('doParallel','splitstackshape','tidyverse')) %dopar%
            {
              theta=theta.list[[s]]
              
              #bootstrapped sample
              d.samp=Combined%>%filter(Species==n.sp[s]) 
              d.samp=stratified(d.samp, "Mesh.size",size=nrow(d.samp),replace=TRUE)
              
              #. objfun to minimize
              fn_ob=function(theta)Selectivty.Kirkwood.Walker(d=d.samp,
                                                              size.int=Size.Interval,
                                                              theta)$negLL
              
              #. fit model
              return(nlminb(theta.list[[s]], fn_ob, gradient = NULL))
            }
          return(exp(do.call(rbind,subListExtract(boot,"par"))))
        }
    })    #takes 0.4 sec per iteration per species
    names(Fit.K_W.CI)=n.sp
    stopCluster(cl)
  }

  
}



  #3.3. Select best fit  

  #Species
Best.fit=vector('list',length(n.sp))
names(Best.fit)=n.sp
for(s in 1:length(n.sp))
{
  if(!names(Fit.M_H)[s]%in%Published.sel.pars_K.W$Species)
  {
    Tab=data.frame(Model=names(Fit.M_H[[s]]$Equal.power),
                   Equal_dev=NA,
                   Prop_dev=NA)
    for(f in 1:length(Rtype))
    {
      Tab$Equal_dev[f]=Fit.M_H[[s]]$Equal.power[[f]]$fit.stats['model_dev']
      Tab$Prop_dev[f]=Fit.M_H[[s]]$Prop.power[[f]]$fit.stats['model_dev']
    }
    Tab1=Tab[which.min(Tab[,2]),-3]%>%mutate(Fishing.power="Equal.power")%>%rename(Dev=Equal_dev)
    Tab2=Tab[which.min(Tab[,3]),-2]%>%mutate(Fishing.power="Prop.power")%>%rename(Dev=Prop_dev)
    Tab3=rbind(Tab1,Tab2)
    if(Do.K_W) Tab3=rbind(Tab1,Tab2,data.frame(Model="K&W",Dev=K.and.W_Dev$Deviance[s],Fishing.power='')) #removed K&W as only Angel shark selected but poor fit
    Best.fit[[s]]=Tab1   #equal fishing power only
    #Best.fit[[s]]=Tab3[which.min(Tab3[,2]),]   #both
  }
}
if(Do.K_W)
{
  K.and.W_Dev=data.frame(Species=n.sp,Model='Gamma_K&W',Deviance=NA)
  K.and.W_Residuals=vector('list',length(n.sp))
  names(K.and.W_Residuals)=n.sp
  for(s in 1:length(n.sp))
  {
    dummy=Selectivty.Kirkwood.Walker(d=Combined%>%filter(Species==n.sp[s]),
                                     size.int=Size.Interval,
                                     Fit.K_W[[s]]$par)
    K.and.W_Dev$Deviance[s]=sum((dummy$observed-dummy$predicted)^2)
    K.and.W_Residuals[[s]]=dummy$observed-dummy$predicted
  }
}

  #Family
Best.fit.family=vector('list',length(n.sp.family))
names(Best.fit.family)=n.sp.family
for(s in 1:length(n.sp.family))
{
  Tab=data.frame(Model=names(Fit.M_H.family[[s]]$Equal.power),
                 Equal_dev=NA,
                 Prop_dev=NA)
  for(f in 1:length(Rtype))
  {
    Tab$Equal_dev[f]=Fit.M_H.family[[s]]$Equal.power[[f]]$fit.stats['model_dev']
    Tab$Prop_dev[f]=Fit.M_H.family[[s]]$Prop.power[[f]]$fit.stats['model_dev']
  }
  Tab1=Tab[which.min(Tab[,2]),-3]%>%mutate(Fishing.power="Equal.power")%>%rename(Dev=Equal_dev)
  Tab2=Tab[which.min(Tab[,3]),-2]%>%mutate(Fishing.power="Prop.power")%>%rename(Dev=Prop_dev)
  Tab3=rbind(Tab1,Tab2)
  Best.fit.family[[s]]=Tab1                              #equal fishing power only
 # Best.fit.family[[s]]=Tab3[which.min(Tab3[,2]),]        #both
}

# EXPORT SELECTIVITY OGIVES  ------------------------------------------------------------------
#note: set to 'equal power' so all meshes go to 1


#function for predicting selectivity
pred.normal.fixed=function(l,k,m,sigma) exp(-((l-k*m)^2)/(2*(sigma)^2))
pred.normal.prop=function(l,m,a1,a2) exp(-(((l- a1*m)^2)/(2*a2*m^2)))
pred.gamma=function(l,m,k,alpha) ((l/((alpha-1)*k*m))^(alpha-1))*exp(alpha-1-(l/(k*m)))
pred.lognormal=function(l,m,m1,mu,sigma) (1/l)*exp(mu+(log(m/m1))-((sigma^2)/2)-((log(l)-mu-log(m/m1))^2)/(2*(sigma^2)))

#Export sel at-length and -age  
out.sel=function(d,BEST,NM,La,Fixed.equal.power)
{
  if(NM%in%Published.sel.pars_K.W$Species)
  {
    PAR=Published.sel.pars_K.W%>%
           filter(Species==NM)
    #predict and export selectivity
      #Plotlen
    dat=pred.Kirkwood.Walker(theta=log(c(PAR$Theta1,PAR$Theta2)),
                             pred.len=PlotLens,
                             Mesh=Combined%>%
                                  filter(Species==NM)%>%
                                  distinct(Mesh.size)%>%
                                  pull(Mesh.size))
    dat=dat%>%
      mutate(TL.mm=Size.class/10)%>%
      dplyr::select(-Size.class)
    colnames(dat)[-ncol(dat)]=round(as.numeric(colnames(dat)[-ncol(dat)])/0.393701,1)
    dat=dat[,c('TL.mm',colnames(dat)[-ncol(dat)])]
    
    if(NM=="Sandbar shark")    #McAuley et al 2007 recommends Lognormal over K&W
    {
      m1=12.7
      mu=4.285
      sigma=0.415 
      dat$'12.7'=pred.lognormal(l=dat$TL.mm,m=12.7,m1,mu,sigma)
      dat$'15.2'=pred.lognormal(l=dat$TL.mm,m=15.2,m1,mu,sigma)
      dat$'16.5'=pred.lognormal(l=dat$TL.mm,m=16.5,m1,mu,sigma)
      dat$'17.8'=pred.lognormal(l=dat$TL.mm,m=17.8,m1,mu,sigma)
      dat$'20.3'=pred.lognormal(l=dat$TL.mm,m=20.3,m1,mu,sigma)
    }
    
    write.csv(dat,paste(handl_OneDrive('Analyses/Data_outs/'),NM,'/',NM,"_gillnet.selectivity",".csv",sep=''),row.names = F)
    dat.out=dat
    
    #Length at age
    dat=pred.Kirkwood.Walker(theta=log(c(PAR$Theta1,PAR$Theta2)),
                             pred.len=La$TL,
                             Mesh=Combined%>%
                                  filter(Species==NM)%>%
                                  distinct(Mesh.size)%>%
                                  pull(Mesh.size))
    dat=dat%>%
      mutate(TL.mm=Size.class/10)%>%
      dplyr::select(-Size.class)%>%
      mutate(Age=La$Age)
    Id=match(c('TL.mm','Age'),colnames(dat))
    colnames(dat)[-Id]=round(as.numeric(colnames(dat)[-Id])/0.393701,1)
    dat=dat[,c('TL.mm','Age',colnames(dat)[-ncol(dat)])]
    
    if(NM=="Sandbar shark")    #McAuley et al 2007 recommends Lognormal over K&W
    {
      m1=12.7
      mu=4.285
      sigma=0.415 
      dat$'12.7'=pred.lognormal(l=dat$TL.mm,m=12.7,m1,mu,sigma)
      dat$'15.2'=pred.lognormal(l=dat$TL.mm,m=15.2,m1,mu,sigma)
      dat$'16.5'=pred.lognormal(l=dat$TL.mm,m=16.5,m1,mu,sigma)
      dat$'17.8'=pred.lognormal(l=dat$TL.mm,m=17.8,m1,mu,sigma)
      dat$'20.3'=pred.lognormal(l=dat$TL.mm,m=20.3,m1,mu,sigma)
    }
    
    
    write.csv(dat,paste(handl_OneDrive('Analyses/Data_outs/'),NM,'/',NM,"_gillnet.selectivity_len.age",".csv",sep=''),row.names = F)
  }
  
  if(!NM%in%Published.sel.pars_K.W$Species)
  {
    if(Fixed.equal.power)
    {
      DAT=d$Equal.power
    }else
    {
      if(BEST$Fishing.power=="Equal.power") DAT=d$Equal.power
      if(BEST$Fishing.power=="Prop.power")  DAT=d$Prop.power
    }
    
    id=match(BEST$Model,names(DAT))
    
    #predict and export selectivity
      #Plotlen
    dat=data.frame(TL.mm=DAT[[id]]$plotlens)
    Pars=d$Equal.power[[id]]$gear.pars
    if(BEST$Model=="norm.loc")
    {
      k=Pars[match("k",rownames(Pars)),1]
      sigma=Pars[match("sigma",rownames(Pars)),1] 
      dat$'15.2'=pred.normal.fixed(l=dat$TL.mm,k,m=15.2,sigma)
      dat$'16.5'=pred.normal.fixed(l=dat$TL.mm,k,m=16.5,sigma)
      dat$'17.8'=pred.normal.fixed(l=dat$TL.mm,k,m=17.8,sigma)
    }
    if(BEST$Model=="norm.sca")
    {
      a1=Pars[match("k1",rownames(Pars)),1]
      a2=Pars[match("k2",rownames(Pars)),1] 
      dat$'15.2'=pred.normal.prop(l=dat$TL.mm,m=15.2,a1,a2)
      dat$'16.5'=pred.normal.prop(l=dat$TL.mm,m=16.5,a1,a2)
      dat$'17.8'=pred.normal.prop(l=dat$TL.mm,m=17.8,a1,a2)
    }
    if(BEST$Model=="gamma")
    {
      k=Pars[match("k",rownames(Pars)),1]
      alpha=Pars[match("alpha",rownames(Pars)),1] 
      dat$'15.2'=pred.gamma(l=dat$TL.mm,m=15.2,k,alpha)
      dat$'16.5'=pred.gamma(l=dat$TL.mm,m=16.5,k,alpha)
      dat$'17.8'=pred.gamma(l=dat$TL.mm,m=17.8,k,alpha)
    }
    if(BEST$Model=="lognorm")
    {
      m1=min(DAT[[id]]$meshsizes)
      mu=Pars[grep('mu1',rownames(Pars)),1]
      sigma=Pars[grep('sigma',rownames(Pars)),1] 
      dat$'15.2'=pred.lognormal(l=dat$TL.mm,m=15.2,m1,mu,sigma)
      dat$'16.5'=pred.lognormal(l=dat$TL.mm,m=16.5,m1,mu,sigma)
      dat$'17.8'=pred.lognormal(l=dat$TL.mm,m=17.8,m1,mu,sigma)
    }
    write.csv(dat,paste(handl_OneDrive('Analyses/Data_outs/'),NM,'/',NM,"_gillnet.selectivity",".csv",sep=''),row.names = F)
    dat.out=dat
    
      #lenght at age
    dat=data.frame(TL.mm=La$TL,Age=La$Age)
    Pars=d$Equal.power[[id]]$gear.pars
    if(BEST$Model=="norm.loc")
    {
      k=Pars[match("k",rownames(Pars)),1]
      sigma=Pars[match("sigma",rownames(Pars)),1] 
      dat$'15.2'=pred.normal.fixed(l=dat$TL.mm,k,m=15.2,sigma)
      dat$'16.5'=pred.normal.fixed(l=dat$TL.mm,k,m=16.5,sigma)
      dat$'17.8'=pred.normal.fixed(l=dat$TL.mm,k,m=17.8,sigma)
    }
    if(BEST$Model=="norm.sca")
    {
      a1=Pars[match("k1",rownames(Pars)),1]
      a2=Pars[match("k2",rownames(Pars)),1] 
      dat$'15.2'=pred.normal.prop(l=dat$TL.mm,m=15.2,a1,a2)
      dat$'16.5'=pred.normal.prop(l=dat$TL.mm,m=16.5,a1,a2)
      dat$'17.8'=pred.normal.prop(l=dat$TL.mm,m=17.8,a1,a2)
    }
    if(BEST$Model=="gamma")
    {
      k=Pars[match("k",rownames(Pars)),1]
      alpha=Pars[match("alpha",rownames(Pars)),1] 
      dat$'15.2'=pred.gamma(l=dat$TL.mm,m=15.2,k,alpha)
      dat$'16.5'=pred.gamma(l=dat$TL.mm,m=16.5,k,alpha)
      dat$'17.8'=pred.gamma(l=dat$TL.mm,m=17.8,k,alpha)
    }
    if(BEST$Model=="lognorm")
    {
      m1=min(DAT[[id]]$meshsizes)
      mu=Pars[grep('mu1',rownames(Pars)),1]
      sigma=Pars[grep('sigma',rownames(Pars)),1] 
      dat$'15.2'=pred.lognormal(l=dat$TL.mm,m=15.2,m1,mu,sigma)
      dat$'16.5'=pred.lognormal(l=dat$TL.mm,m=16.5,m1,mu,sigma)
      dat$'17.8'=pred.lognormal(l=dat$TL.mm,m=17.8,m1,mu,sigma)
    }
    write.csv(dat,paste(handl_OneDrive('Analyses/Data_outs/'),NM,'/',NM,"_gillnet.selectivity_len.age",".csv",sep=''),row.names = F)
    
  }
  return(dat.out)
}

    #species
Pred.sels.main.mesh=vector('list',length(n.sp))
names(Pred.sels.main.mesh)=n.sp
for(s in 1:length(n.sp))
{
  NM=n.sp[s]
  if(NM=='Southern eagle ray') NM='Eagle ray'
  Pred.sels.main.mesh[[s]]=out.sel(d=Fit.M_H[[s]],
                                   BEST=Best.fit[[s]],
                                   NM=NM,
                                   La=LatAge[[s]],
                                   Fixed.equal.power=TRUE)  
}

  #family
Pred.sels.main.mesh.family=vector('list',length(n.sp.family))
names(Pred.sels.main.mesh.family)=n.sp.family
for(s in 1:length(n.sp.family))
{
  Pred.sels.main.mesh.family[[s]]=out.sel(d=Fit.M_H.family[[s]],
                                            BEST=Best.fit.family[[s]],
                                            NM=n.sp.family[s],
                                            La=LatAge.family[[s]],
                                            Fixed.equal.power=TRUE)  
}


# REPORT  ------------------------------------------------------------------
if(do.paper.figures)
{
  setwd(handl_OneDrive('Analyses/Selectivity_Gillnet'))
  
  #Export best fist
  write.csv(do.call(rbind,Best.fit),"Table best model.csv",row.names = T)

  
  colfunc <- colorRampPalette(c("white", "yellow","orange",'brown2',"darkred"))
  
  #1. Map
  do.map=FALSE
  if(do.map)
  {
    library("rnaturalearth")
    world <- ne_countries(scale = "medium", returnclass = "sf")
    Map.dat=rbind(F1_SamplingTwo%>%
                    filter(Mesh%in%c('S4','S5','S6','S7','S8'))%>%
                    distinct(Cruise,Station,.keep_all = T)%>%
                    mutate(mid.lat=-abs(StartLat/100),
                           mid.long=StartLong/100)%>%
                    dplyr::select(mid.lat,mid.long),
                  Exp.net.WA%>%
                    distinct(sheet_no,.keep_all = T)%>%
                    dplyr::select(mid.lat,mid.long)%>%
                    mutate(mid.lat=-abs(mid.lat))
    )%>%
      filter(mid.long > 111 & !is.na(mid.lat))
    Long.range=range(Map.dat$mid.long)
    Lat.range=range(Map.dat$mid.lat)
    Poly=data.frame(x=c(Long.range,rev(Long.range)),y=c(rep(Lat.range[1],2),rep(Lat.range[2],2)))
    Inset=ggplot(data = world) +
      geom_sf(color = "black", fill = "white") +
      coord_sf(xlim =c(Long.range[1],153) , ylim = c(Lat.range[1],-11), expand = T)+
      theme(axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank())+ 
      annotate(geom="text", x=134, y=-24, label="Australia",size=12)+
      geom_polygon(data=Poly,aes(x,y),fill='black',alpha=.35)
    p=ggplot(data = world) +
      geom_sf(color = "black", fill = "grey60") +
      coord_sf(xlim =Long.range , ylim = Lat.range, expand = T) +
      xlab("Longitude") + ylab("Latitude")+
      geom_point(data=Map.dat,aes(x=mid.long, y=mid.lat,size=10),shape=21,alpha=0.4,fill="brown4")+
      theme(legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.position = "none",
            legend.key=element_blank(),
            legend.background=element_blank(),
            legend.direction = "vertical",
            legend.box = "horizontal",
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            plot.margin=unit(c(.1,.1,.1,.1),"cm"))
    p+
      annotation_custom(ggplotGrob(Inset),xmin = 112, xmax = 130, ymin = -44, ymax = -36)
    ggsave("Figure 1.tiff", width = 12,height = 6,dpi = 300, compression = "lzw")
    
  }
  
  
  #2. Appendix 1
  Appendix1=TL_FL%>%
    filter(name%in%n.sp)%>%
    dplyr::select(-get.this)%>%
    mutate(intercept=round(intercept,2),
           slope=round(slope,2),
           name=capitalize(tolower(name)))%>%
    arrange(name)%>%
    distinct(name,.keep_all = T)
  Appendix1=rbind(Appendix1,LH%>%
          filter(Name%in%n.sp)%>%
          filter(!Name%in%Appendix1$name)%>%
          dplyr::select(Name,a_FL.to.TL,b_FL.to.TL)%>%
          rename(name=Name,
                 slope=a_FL.to.TL,
                 intercept=b_FL.to.TL))%>%
    filter(!name%in%Published.sel.pars_K.W$Species)%>%
    arrange(name)
  write.csv(Appendix1,'Appendix1.csv',row.names = F)
  
  
  #3. Table 1. All species observed
  write.csv(Table1,'Table1.csv',row.names = F)
  
  
  #4. Select species without selectivity published
  n.sp.pub=n.sp[-match(Published,n.sp)]
  
  
  #5. Table mean size by mesh of selected species
    #species
  Tab2.mean=Combined%>%
    group_by(Mesh.size,Species)%>%
    summarise(Mean=round(mean(Length),1))%>%
    spread(Mesh.size,Mean)%>%
    data.frame
  colnames(Tab2.mean)[-1]=substr(colnames(Tab2.mean)[-1],2,10)
  write.csv(Tab2.mean,'Table mean size by mesh selected species.csv',row.names = F)
  
  cols=colfunc(length(unique(Combined$Mesh.size)))
  
  Combined%>%
    filter(Species%in%n.sp.pub)%>%
    group_by(Mesh.size,Species)%>%
    summarise(Mean=mean(Length),
              SD=sd(Length))%>%
    mutate(Mesh=as.character(Mesh.size))%>%
    ggplot(aes(x=Mesh.size, y= Mean, colour=Mesh)) + 
    geom_errorbar(aes(ymin= Mean-SD, ymax= Mean+SD), width=.1) +
    geom_point() + ylab("Mean total length (+/-SD)") + xlab("Mesh size (cm)") +
    facet_wrap(vars(Species), scales = "free_y")+
    scale_color_manual("Mesh size (cm)",values=cols)+
    theme_dark()+ 
    theme(legend.position="top",
          strip.text.x = element_text(size = 11,face='bold',color="white"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))
  ggsave('Figure.S1.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")
  
  Combined%>%
    group_by(Mesh.size,Species)%>%
    summarise(Mean=mean(Length),
              SD=sd(Length))%>%
    mutate(Mesh=as.character(Mesh.size))%>%
    ggplot(aes(x=Mesh.size, y= Mean, colour=Mesh)) + 
    geom_errorbar(aes(ymin= Mean-SD, ymax= Mean+SD), width=.1) +
    geom_point() + ylab("Mean total length (+/-SD)") + xlab("Mesh size (cm)") +
    facet_wrap(vars(Species), scales = "free_y")+
    scale_color_manual("Mesh size (cm)",values=cols)+
    theme_dark()+ 
    theme(strip.text.x = element_text(size = 11,face='bold',color="white"))
  ggsave('Mean size by mesh.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")
  
    #family
  Tab2.mean=Combined.family%>%
    filter(Species%in%n.sp.family)%>%
    group_by(Mesh.size,Species)%>%
    summarise(Mean=round(mean(Length),1))%>%
    spread(Mesh.size,Mean)%>%
    data.frame
  colnames(Tab2.mean)[-1]=substr(colnames(Tab2.mean)[-1],2,10)
  Combined.family%>%
    filter(Species%in%n.sp.family)%>%
    group_by(Mesh.size,Species)%>%
    summarise(Mean=mean(Length),
              SD=sd(Length))%>%
    mutate(Mesh=as.character(Mesh.size))%>%
    ggplot(aes(x=Mesh.size, y= Mean, colour=Mesh)) + 
    geom_errorbar(aes(ymin= Mean-SD, ymax= Mean+SD), width=.1) +
    geom_point() + ylab("Mean total length (+/-SD)") + xlab("Mesh size (cm)") +
    facet_wrap(vars(Species), scales = "free_y")+
    scale_color_manual("Mesh size (cm)",values=colfunc(length(unique(Combined.family$Mesh.size))))+
    theme_dark()+ 
    theme(legend.position="top",
          strip.text.x = element_text(size = 11,face='bold',color="white"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))
  ggsave('Figure.S2_family.tiff', width = 10,height = 7, dpi = 300, compression = "lzw")
  

  #6. Display size frequencies 
  library(lemon)
  fig2=function(d)
  {
    cols=colfunc(length(unique(d$MESH)))
    p=d %>%
      ggplot( aes(x=TL, fill=MESH)) +
      geom_histogram(binwidth = 10, alpha=0.8,colour='grey40',size=.1)  +
      labs(fill="Mesh size (cm)")+
      facet_wrap(vars(Species), scales = "free_y")+
      xlab("Total length (cm)")+ ylab("Frequency")+scale_fill_manual(values=cols)+
      theme(legend.title = element_text(size = 10),
            legend.text = element_text(size = 9),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18))
  }

    #species in paper
  p=fig2(d=Combined%>%
         mutate(MESH=as.factor(Mesh.size),
                TL=Length)%>%
        filter(Species%in%n.sp.pub))
  p=p+theme(legend.position="top",
            legend.title=element_text(size = 14),
            legend.text = element_text(size = 12),
            strip.text.x = element_text(size = 14),
            axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 13))
  p+ guides(fill = guide_legend(nrow = 1))
  #p=p+theme(legend.position="bottom")
  #reposition_legend(p , 'center',panel='panel-4-3')   #Manually save figrue as 'Figure 2.tiff', ggsave not working with reposition_legend
  ggsave('Figure 2.tiff', width = 12,height = 10, dpi = 300, compression = "lzw")
  
      #all species
  p=fig2(d=Combined%>%
         mutate(MESH=as.factor(Mesh.size),
                TL=Length))
  p
  ggsave('Size frequency_all selected species_data set combined.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")

      #family
  p=fig2(d=Combined.family%>%
         mutate(MESH=as.factor(Mesh.size),
                TL=Length)%>%
         filter(Species%in%n.sp.family))
  p=p+theme(legend.position="top",
            legend.title=element_text(size = 14),
            legend.text = element_text(size = 12),
            strip.text.x = element_text(size = 14),
            axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 13))
  p
  #p=p+theme(legend.position="bottom")
  #reposition_legend(p , 'right',panel='panel-3-3')
  #Manually save figrue as 'Figure 3.tiff', ggsave not working with reposition_legend
  ggsave('Figure 3_Size frequency_family.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")
  
  
  
  #7. Display density
  fig.density=function(d)
  {
    cols=colfunc(length(unique(d$MESH)))
    d %>%
      ggplot( aes(x=TL, fill=MESH)) +
      geom_density(alpha=0.8)  +
      labs(fill="Mesh size (cm)")+
      facet_wrap(vars(Species), scales = "free")+
      ylab("Density")+ xlab("Total length (cm)")+
      scale_fill_manual("Mesh size (cm)", values = cols)
  }
  
      #all species
  fig.density(d=Combined%>%
         mutate(MESH=as.factor(Mesh.size),
                TL=Length))
  ggsave('Density all selected species_data set combined.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")

    #family
  fig.density(d=Combined.family%>%
                mutate(MESH=as.factor(Mesh.size),
                       TL=Length))
  ggsave('Density family_data set combined.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")

  
  #8. Output parameter estimates from all models 
  fn.rnd=function(x) sprintf(round(x,2), fmt = '%#.2f')
  
    #species
  Table.mod.fit=vector('list',length(n.sp))
  names(Table.mod.fit)=n.sp
  for(s in 1:length(n.sp))
  {
    if(!names(Fit.M_H)[s]%in%Published.sel.pars_K.W$Species)
    {
      Tab=data.frame(Model=names(Fit.M_H[[s]]$Equal.power),
                     Equal_Param1=NA,Equal_Param2=NA,Equal_Deviance=NA,
                     Prop_Param1=NA,Prop_Param2=NA,Prop_Deviance=NA)
      for(f in 1:length(Rtype))
      {
        #equal power
        Del=ifelse(sum(is.na(Fit.M_H[[s]]$Equal.power[[f]]$gear.pars[1:2,]))>0,"YES","NO")
        if(Del=="NO")
        {
          Tab$Equal_Deviance[f]=fn.rnd(Fit.M_H[[s]]$Equal.power[[f]]$fit.stats['model_dev'])
          PaR=fn.rnd(Fit.M_H[[s]]$Equal.power[[f]]$gear.pars[1:2,'estimate'])
          errOr=fn.rnd(Fit.M_H[[s]]$Equal.power[[f]]$gear.pars[1:2,'s.e.'])
          Tab$Equal_Param1[f]=paste(PaR[1]," (",errOr[1],")",sep='')
          Tab$Equal_Param2[f]=paste(PaR[2]," (",errOr[2],")",sep='')
          rm(PaR,errOr)
        }
        
        #Prop power
        Del=ifelse(sum(is.na(Fit.M_H[[s]]$Prop.power[[f]]$gear.pars[1:2,]))>0,"YES","NO")
        if(Del=="NO")
        {
          Tab$Prop_Deviance[f]=fn.rnd(Fit.M_H[[s]]$Prop.power[[f]]$fit.stats['model_dev'])
          PaR=fn.rnd(Fit.M_H[[s]]$Prop.power[[f]]$gear.pars[1:2,'estimate'])
          errOr=fn.rnd(Fit.M_H[[s]]$Prop.power[[f]]$gear.pars[1:2,'s.e.'])
          Tab$Prop_Param1[f]=paste(PaR[1]," (",errOr[1],")",sep='')
          Tab$Prop_Param2[f]=paste(PaR[2]," (",errOr[2],")",sep='')
          rm(PaR,errOr)
        }
      }
      Table.mod.fit[[s]]=Tab
    }

  }
  Table.mod.fit=do.call(rbind,Table.mod.fit)
  Table.mod.fit$Species=sub("*\\.[0-9]", "", rownames(Table.mod.fit))
  Table.mod.fit=Table.mod.fit%>%
    mutate(Model=ifelse(Model=="norm.loc","Normal (fixed spread)",
                        ifelse(Model=="norm.sca","Normal (prop. spread)",
                               ifelse(Model=="gamma","Gamma",
                                      ifelse(Model=="lognorm","Lognormal",
                                             Model)))))%>%
    dplyr::select(Species,Model,Equal_Param1,Equal_Param2,Equal_Deviance,
                  Prop_Param1,Prop_Param2,Prop_Deviance)
  if(Do.K_W) #add K&W 
  {
    Table2=data.frame(Species=n.sp,
                      Theta1=NA,Theta1.LOW95=NA,Theta1.UP95=NA,
                      Theta2=NA,Theta2.LOW95=NA,Theta2.UP95=NA)   
    for(s in 1:length(n.sp))
    {
      dummy=Fit.K_W.CI[[s]]
      Table2$Theta1[s]=round(quantile(dummy[,"Theta1"],probs=0.5),1)
      Table2$Theta1.LOW95[s]=round(quantile(dummy[,"Theta1"],probs=0.025),1)
      Table2$Theta1.UP95[s]=round(quantile(dummy[,"Theta1"],probs=0.975),1)
      Table2$Theta2[s]=round(quantile(dummy[,"Theta2"],probs=0.5),1)
      Table2$Theta2.LOW95[s]=round(quantile(dummy[,"Theta2"],probs=0.025),1)
      Table2$Theta2.UP95[s]=round(quantile(dummy[,"Theta2"],probs=0.975),1)
    }
    Table2=Table2%>%mutate(Theta1.SE=fn.rnd((Theta1.UP95-Theta1)/1.96),
                           Theta2.SE=fn.rnd((Theta2.UP95-Theta2)/1.96),
                           Theta1=fn.rnd(Theta1),
                           Theta2=fn.rnd(Theta2))
    Add.K.and.W=K.and.W_Dev%>%
      mutate(Equal_Param1=paste(Table2$Theta1," (",Table2$Theta1.SE,")",sep=''),
             Equal_Param2=paste(Table2$Theta2," (",Table2$Theta2.SE,")",sep=''),
             Prop_Param1=NA,
             Prop_Param2=NA,
             Prop_Deviance=NA,
             Deviance=fn.rnd(Deviance))%>%
      rename(Equal_Deviance=Deviance)%>%
      dplyr::select(names(Table.mod.fit))
    
    Table.mod.fit=rbind(Table.mod.fit,Add.K.and.W)%>%
      arrange(Species)
    
  }
  
  Table.mod.fit$Species[duplicated(Table.mod.fit$Species)] <- ""
  Table.mod.fit[is.na(Table.mod.fit)] <- ""
  Table.mod.fit[Table.mod.fit=="NaN"] <- ""
  write.csv(Table.mod.fit,"Table2.csv",row.names=F)
  # fn.word.table(WD=getwd(),TBL=Table.mod.fit,Doc.nm="Table2",caption=NA,paragph=NA,
  #               HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
  #               Zebra='NO',Zebra.col='grey60',Grid.col='black',
  #               Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
    #family
  Table.mod.fit.family=vector('list',length(n.sp.family))
  names(Table.mod.fit.family)=n.sp.family
  for(s in 1:length(n.sp.family))
  {
    Tab=data.frame(Model=names(Fit.M_H.family[[s]]$Equal.power),
                   Equal_Param1=NA,Equal_Param2=NA,Equal_Deviance=NA,
                   Prop_Param1=NA,Prop_Param2=NA,Prop_Deviance=NA)
    for(f in 1:length(Rtype))
    {
      #equal power
      Del=ifelse(sum(is.na(Fit.M_H.family[[s]]$Equal.power[[f]]$gear.pars[1:2,]))>0,"YES","NO")
      if(Del=="NO")
      {
        Tab$Equal_Deviance[f]=fn.rnd(Fit.M_H.family[[s]]$Equal.power[[f]]$fit.stats['model_dev'])
        PaR=fn.rnd(Fit.M_H.family[[s]]$Equal.power[[f]]$gear.pars[1:2,'estimate'])
        errOr=fn.rnd(Fit.M_H.family[[s]]$Equal.power[[f]]$gear.pars[1:2,'s.e.'])
        Tab$Equal_Param1[f]=paste(PaR[1]," (",errOr[1],")",sep='')
        Tab$Equal_Param2[f]=paste(PaR[2]," (",errOr[2],")",sep='')
        rm(PaR,errOr)
      }
      
      #Prop power
      Del=ifelse(sum(is.na(Fit.M_H.family[[s]]$Prop.power[[f]]$gear.pars[1:2,]))>0,"YES","NO")
      if(Del=="NO")
      {
        Tab$Prop_Deviance[f]=fn.rnd(Fit.M_H.family[[s]]$Prop.power[[f]]$fit.stats['model_dev'])
        PaR=fn.rnd(Fit.M_H.family[[s]]$Prop.power[[f]]$gear.pars[1:2,'estimate'])
        errOr=fn.rnd(Fit.M_H.family[[s]]$Prop.power[[f]]$gear.pars[1:2,'s.e.'])
        Tab$Prop_Param1[f]=paste(PaR[1]," (",errOr[1],")",sep='')
        Tab$Prop_Param2[f]=paste(PaR[2]," (",errOr[2],")",sep='')
        rm(PaR,errOr)
      }
    }
    Table.mod.fit.family[[s]]=Tab
  }
  Table.mod.fit.family=do.call(rbind,Table.mod.fit.family)
  Table.mod.fit.family$Species=sub("*\\.[0-9]", "", rownames(Table.mod.fit.family))
  Table.mod.fit.family=Table.mod.fit.family%>%
    mutate(Model=ifelse(Model=="norm.loc","Normal (fixed spread)",
                        ifelse(Model=="norm.sca","Normal (prop. spread)",
                               ifelse(Model=="gamma","Gamma",
                                      ifelse(Model=="lognorm","Lognormal",
                                             Model)))))%>%
    dplyr::select(Species,Model,Equal_Param1,Equal_Param2,Equal_Deviance,
                  Prop_Param1,Prop_Param2,Prop_Deviance)
  Table.mod.fit.family$Species[duplicated(Table.mod.fit.family$Species)] <- ""
  Table.mod.fit.family[is.na(Table.mod.fit.family)] <- ""
  Table.mod.fit.family[Table.mod.fit.family=="NaN"] <- ""
  write.csv(Table.mod.fit.family,"Table2.family.csv",row.names=F)
  # fn.word.table(WD=getwd(),TBL=Table.mod.fit.family,Doc.nm="Table2.family",caption=NA,paragph=NA,
  #               HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
  #               Zebra='NO',Zebra.col='grey60',Grid.col='black',
  #               Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")


  #9. Plot residuals for all model
  
  #all together
    #Published species
  nn=match(n.sp.pub,names(Fit.M_H))
  nn1=floor(length(nn)/2)
  nn2=ceiling(length(nn)/2)
  tiff(file="Figure.S3.a_Fit.tiff",width = 1800, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(nn1,length(Rtype)),mar=c(1.5,1.2,.2,.3),oma=c(1.5,2,.7,1),mgp=c(1,.5,0))
  for(s in nn[1:nn1])
  {
     for(f in 1:length(Fit.M_H[[s]]$Equal.power))
    {
      with(Fit.M_H[[s]]$Equal.power[[f]],
      {
        MAIN=""
        plot.resids(devres,meshsizes,lens,title=MAIN)
      })
       if(s==nn[1])
       {
         MAIN=with(Fit.M_H[[s]]$Equal.power[[f]],
                   ifelse(type=="norm.loc","Normal (fixed spread)",
                   ifelse(type=="norm.sca","Normal (prop. spread)",
                   ifelse(type=="gamma","Gamma",
                   ifelse(type=="lognorm","Lognormal",NA)))))
         mtext(MAIN,3,cex=.8)
       }
     }
    mtext( names(Fit.M_H)[s],4,cex=.75)
  }
  mtext("Total length (mm)",1,outer=T,line=.35,cex=1.25)
  mtext("Mesh size (cm)",2,outer=T,line=.35,cex=1.25,las=3)
  dev.off()
  tiff(file="Figure.S3.b_Fit.tiff",width = 1800, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(nn2,length(Rtype)),mar=c(1.5,1.2,.2,.3),oma=c(1.5,2,.7,1),mgp=c(1,.5,0),xpd=T)
  for(s in nn[(nn1+1):(nn1+nn2)])
  {
    for(f in 1:length(Fit.M_H[[s]]$Equal.power))
    {
      with(Fit.M_H[[s]]$Equal.power[[f]],
           {
             MAIN=""
             plot.resids(devres,meshsizes,lens,title=MAIN)
           })
      if(s==nn[(nn1+1)])
      {
        MAIN=with(Fit.M_H[[s]]$Equal.power[[f]],
                  ifelse(type=="norm.loc","Normal (fixed spread)",
                         ifelse(type=="norm.sca","Normal (prop. spread)",
                                ifelse(type=="gamma","Gamma",
                                       ifelse(type=="lognorm","Lognormal",NA)))))
        mtext(MAIN,3,cex=.8)
      }
    }
    mtext( names(Fit.M_H)[s],4,cex=.75)
  }
  mtext("Total length (mm)",1,outer=T,line=.35,cex=1.25)
  mtext("Mesh size (cm)",2,outer=T,line=.35,cex=1.25,las=3)
  dev.off()
  
  #Target species
  do.this=FALSE
  if(do.this)
  {
    nn=match(Published,names(Fit.M_H))
    tiff(file="Fit_Target.species.tiff",width = 1800, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfrow=c(length(nn),length(Rtype)),mar=c(1.5,1.2,.2,.3),oma=c(1.5,2,.7,1),mgp=c(1,.5,0))
    for(s in nn)
    {
      
      for(f in 1:length(Fit.M_H[[s]]$Equal.power))
      {
        with(Fit.M_H[[s]]$Equal.power[[f]],
             {
               MAIN=""
               plot.resids(devres,meshsizes,lens,title=MAIN)
             })
        if(s==1)
        {
          MAIN=with(Fit.M_H[[s]]$Equal.power[[f]],
                    ifelse(type=="norm.loc","Normal (fixed spread)",
                           ifelse(type=="norm.sca","Normal (prop. spread)",
                                  ifelse(type=="gamma","Gamma",
                                         ifelse(type=="lognorm","Lognormal",NA)))))
          mtext(MAIN,3,cex=.8)
        }
      }
      mtext( names(Fit.M_H)[s],4,cex=.95,line=.2)
    }
    mtext("Total length (mm)",1,outer=T,line=.35,cex=1.25)
    mtext("Mesh size (cm)",2,outer=T,line=.35,cex=1.25,las=3)
    dev.off()
    
  }
  
  
  #Each individual species
  for(s in 1:length(n.sp))
  {
    if(!names(Fit.M_H)[s]%in%Published.sel.pars_K.W$Species)
    {
      tiff(file=paste("Each species/Fit/species/Deviance residual plot_",n.sp[s],".tiff",sep=''),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
      smart.par(length(Rtype),MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.75,3,.5,.1),MGP=c(1,.5,0))
      for(f in 1:length(Fit.M_H[[s]]$Equal.power))
      {
        
        with(Fit.M_H[[s]]$Equal.power[[f]],
             {
               MAIN=ifelse(type=="norm.loc","Normal (fixed spread)",
                           ifelse(type=="norm.sca","Normal (prop. spread)",
                                  ifelse(type=="gamma","Gamma",
                                         ifelse(type=="lognorm","Lognormal",NA))))
               plot.resids(devres,meshsizes,lens,title=MAIN)
             })
      }
      mtext("Total length (mm)",1,outer=T,line=.5,cex=1.2)
      mtext("Mesh size (cm)",2,outer=T,line=1.25,cex=1.2,las=3)
      mtext( names(Fit.M_H)[s],3,outer=T,line=-.75,cex=1.2)
      dev.off() 
    }
  }

  
    #Family
      #all together
  tiff(file="Figure.S4_family.tiff",width = 1800, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(length(n.sp.family),length(Rtype)),mar=c(1.5,1.2,.2,.3),oma=c(1.5,2,.7,1.2),mgp=c(1,.5,0))
  for(s in 1:length(n.sp.family))
  {
    for(f in 1:length(Fit.M_H.family[[s]]$Equal.power))
    {
      with(Fit.M_H.family[[s]]$Equal.power[[f]],
           {
             MAIN=""
              plot.resids(devres,meshsizes,lens,title=MAIN)
           })
      if(s==1)
      {
        MAIN=with(Fit.M_H[[s]]$Equal.power[[f]],
                  ifelse(type=="norm.loc","Normal (fixed spread)",
                         ifelse(type=="norm.sca","Normal (prop. spread)",
                                ifelse(type=="gamma","Gamma",
                                       ifelse(type=="lognorm","Lognormal",NA)))))
        mtext(MAIN,3,cex=.8)
      }
    }
    mtext( names(Fit.M_H.family)[s],4,line=.2,cex=.85)
  }
  mtext("Total length (mm)",1,outer=T,line=.35,cex=1.25)
  mtext("Mesh size (cm)",2,outer=T,line=.35,cex=1.25,las=3)
  dev.off()
  
      #each individual family
  for(s in 1:length(n.sp.family))
  {
    tiff(file=paste("Each species/Fit/family/Deviance residual plot_",n.sp.family[s],".tiff",sep=''),width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
    smart.par(length(Rtype),MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.75,3,.5,.1),MGP=c(1,.5,0))
    for(f in 1:length(Fit.M_H.family[[s]]$Equal.power))
    {
      with(Fit.M_H.family[[s]]$Equal.power[[f]],
           {
             MAIN=ifelse(type=="norm.loc","Normal (fixed spread)",
                  ifelse(type=="norm.sca","Normal (prop. spread)",
                  ifelse(type=="gamma","Gamma",
                  ifelse(type=="lognorm","Lognormal",NA))))
             plot.resids(devres,meshsizes,lens,title=MAIN)
           })
    }
    mtext("Total length (mm)",1,outer=T,line=.5,cex=1.2)
    mtext("Mesh size (cm)",2,outer=T,line=1.25,cex=1.2,las=3)
    mtext( names(Fit.M_H.family)[s],3,outer=T,line=-.75,cex=1.2)
    dev.off()
  }
  
  
  #10. Observed vs predicted number at size by mesh for Kirkwood & Walker
  if(Do.K_W)
  {
    tiff(file=paste("Each species/Fit/K&W/Deviance residual plot.tiff",sep=''),width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
    smart.par(length(n.sp),MAR=c(1.5,1.2,1,.1),OMA=c(1.75,3,.1,1),MGP=c(1,.5,0))
    for(s in 1:length(n.sp))
    {
      with(Fit.M_H[[s]]$Equal.power[[1]],plot.resids(K.and.W_Residuals[[s]],meshsizes,lens,title=names(Fit.M_H)[s]))
    }
    mtext("Total length (mm)",1,outer=T,line=.35,cex=1.25)
    mtext("Mesh size (cm)",2,outer=T,line=.35,cex=1.25,las=3)
    dev.off()
    
    for(s in 1:length(n.sp))
   {
    dummy=Selectivty.Kirkwood.Walker(d=Combined%>%filter(Species==n.sp[s]),
                                     size.int=Size.Interval,
                                     Fit.K_W[[s]]$par)
    tiff(file=paste("Each species/Fit/K&W/Pre_vs_Obs_",n.sp[s],".tiff",sep=''),width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
    smart.par(n.plots=ncol(dummy$observed),MAR=c(2,1,2,1.5),OMA=c(2.5,3,.5,.1),MGP=c(1,.5,0))
    for(m in 1:ncol(dummy$observed))
    {
      DAT=data.frame(y=dummy$observed[,m],x=dummy$predicted[,m])
      mod=lm(y~x,data=DAT)
      CoF=coef(mod)
      Smry=summary(mod)
      Ndat=data.frame(x=seq(min(DAT$x),max(DAT$x),length.out = 10))
      Prd=predict(mod,newdata=Ndat,se.fit=TRUE)
      Main=paste(round(as.numeric(gsub("[^0-9.]", "",colnames(dummy$observed)[m]))/2.54,1),'inch')
      plot(dummy$predicted[,m],dummy$observed[,m],pch=19,cex=1.5,ylab='',xlab='',main=Main)
      lines(Ndat$x,Prd$fit,lwd=1.5)
      text(mean(Ndat$x),mean(Prd$fit)*.6,
           paste('y= ',round(CoF[1],2),' + ',round(CoF[2],2),'x',sep=''),pos=4,cex=1.25)
      mylabel = bquote(italic(r)^2 == .(format(round(Smry$r.squared,2), digits = 3)))
      text(mean(Ndat$x),mean(Prd$fit)*.4,mylabel,pos=4,cex=1.25)
      polygon(x=c(Ndat$x,rev(Ndat$x)),
              y=c(Prd$fit+1.96*Prd$se.fit,rev(Prd$fit-1.96*Prd$se.fit)),
              col=rgb(.1,.1,.1,alpha=.2),border='transparent')
    }
    mtext("Predicted size class (mm)",1,outer=T,line=1)
    mtext("Observed size class (mm)",2,outer=T,las=3,line=1.5)
    dev.off()  
  }
  }

  
  #11. Plot Observed size frequency VS estimated selectivity
  Cols.type=c('blue','brown','forestgreen','red')
  names(Cols.type)=Rtype
  if(Do.K_W)
  {
    Cols.type=c(Cols.type,'steelblue')
    names(Cols.type)=c(Rtype,"K&W")
  }
  
    #Published species 
  fn.freq.obs.pred1=function(d,d.KW=NULL,net.names)
  {
    #Calculated expected population frequency
    Obs=d$tab[,-match('Size.class',names(d$tab))]
    N.fish.caught=vector('list',length(d$Equal.power))
    names(N.fish.caught)=names(d$Equal.power)
    for(f in 1:length(d$Equal.power))
    {
      Lns=d$Equal.power[[f]]$lens
      iid=match(Lns,d$Equal.power[[f]]$plotlens)
      Pred.sel=d$Equal.power[[f]]$rselect[iid,]
      Total.sel=rowSums(Pred.sel)
      Total.obs=rowSums(Obs)
      Rel.num.in.pop=sapply(Total.obs/Total.sel,function(x) max(x,0.1))
      Rel.prop.in.pop=Rel.num.in.pop/sum(Rel.num.in.pop)
      N.mesh=colSums(Obs)
      Expnd.N.mesh=matrix(rep(N.mesh,each=nrow(Pred.sel)),ncol=ncol(Pred.sel))
      Expnd.Rel.prop.in.pop=matrix(rep(Rel.prop.in.pop,ncol(Pred.sel)),ncol=ncol(Pred.sel))
      Sum.prod=colSums(Pred.sel*Expnd.Rel.prop.in.pop)
      Sum.prod=matrix(rep(Sum.prod,each=nrow(Pred.sel)),ncol=ncol(Pred.sel))
      N.fish.caught[[f]]=(Expnd.N.mesh*Expnd.Rel.prop.in.pop*Pred.sel)/Sum.prod
    }
    if(!is.null(d.KW))
    {
      iid=match(Lns,d.KW$Size.class)
      
      Pred.sel=d.KW[iid,-1]
      Total.sel=rowSums(Pred.sel)
      Total.obs=rowSums(Obs)
      Rel.num.in.pop=sapply(Total.obs/Total.sel,function(x) max(x,0.1))
      Rel.prop.in.pop=Rel.num.in.pop/sum(Rel.num.in.pop)
      N.mesh=colSums(Obs)
      Expnd.N.mesh=matrix(rep(N.mesh,each=nrow(Pred.sel)),ncol=ncol(Pred.sel))
      Expnd.Rel.prop.in.pop=matrix(rep(Rel.prop.in.pop,ncol(Pred.sel)),ncol=ncol(Pred.sel))
      Sum.prod=colSums(Pred.sel*Expnd.Rel.prop.in.pop)
      Sum.prod=matrix(rep(Sum.prod,each=nrow(Pred.sel)),ncol=ncol(Pred.sel))
      N.fish.caught$'K&W'=(Expnd.N.mesh*Expnd.Rel.prop.in.pop*Pred.sel)/Sum.prod
    }
    
    #Plot by mesh
    n=ncol(d$tab)-1
    Msh=substr(names(d$tab)[-1],2,10)
    for(i in 1:length(Nets))
    {
      if(!as.character(Nets[i])%in%Msh) plot.new()
      if(as.character(Nets[i])%in%Msh)
      {
        ii=match(paste('X',Nets[i],sep=''),colnames(d$tab))
        dmax=ceiling(max(c(max(d$tab[,ii]),max(N.fish.caught[[3]][,ii-1]))))
        plot(d$tab$Size.class,d$tab[,ii],type='h',ylab='',xlab='', lwd = 4,
             col=rgb(.5,.5,.5,alpha=.5), yaxt = "n",
             ylim=c(0,dmax))
        for(m in 1:length(N.fish.caught))
        {
          lines(d$tab$Size.class,N.fish.caught[[m]][,ii-1],col=Cols.type[m],lwd=1.5)
        }
        RANGO=round(seq(0,dmax,length.out=4))
        if(nchar(RANGO[4])>=2 & RANGO[4]>30) RANGO=10*round(RANGO/10)
        axis(2, at = RANGO,las=2)
      }
      if(net.names) mtext(paste(Nets[i],'cm'),3,cex=.9)
      
    }
  }
  Nets=sort(unique(Combined$Mesh.size))
  nn=match(n.sp.pub,names(Fit.M_H))
  nn1=floor(length(nn)/2)
  nn2=ceiling(length(nn)/2)
  
  
  tiff(file="Figure.S5.a_Obs.vs.Pred.tiff",width = 2400, height = 2200,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(nn1,length(Nets)),mar=c(1.5,1.2,.2,.65),oma=c(1.5,2,.8,1.2),mgp=c(1,.5,0),las=1,xpd=T)
  for(s in nn[1:nn1])
  {
    dummy=NULL
    net.names=FALSE
    if(s==nn[1]) net.names=TRUE
    if(Do.K_W) dummy=Pred.sel.K_W[[s]]
    fn.freq.obs.pred1(d=Fit.M_H[[s]],
                      d.KW=dummy,
                      net.names=net.names)
    #if(s==nn[nn1-1])
    if(s==1)
    {
      Nms=names(Cols.type)
      Nms=ifelse(Nms=="norm.loc","Normal fixed",
                 ifelse(Nms=="norm.sca","Normal prop.",
                        ifelse(Nms=="gamma","Gamma",
                               ifelse(Nms=="lognorm","Lognormal",NA))))
      legend("left",Nms,bty='n',text.col=Cols.type,cex=1)
    }
    mtext(names(Fit.M_H)[s],4,line=.5,cex=.8,las=3)
  }
  mtext("Frequency", side = 2,outer=T, line = 0.2,las=3,cex=1.2)
  mtext("Total length (cm)",1,outer=T,line=.4,cex=1.2)
  dev.off()
  
  tiff(file="Figure.S5.b_Obs.vs.Pred.tiff",width = 2400, height = 2200,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(nn2,length(Nets)),mar=c(1.5,1.2,.2,.65),oma=c(1.5,2,.8,1.2),mgp=c(1,.5,0),las=1,xpd=T)
  for(s in nn[(nn1+1):(nn1+nn2)])
  {
    dummy=NULL
    net.names=FALSE
    if(s==nn[(nn1+1)]) net.names=TRUE
    if(Do.K_W) dummy=Pred.sel.K_W[[s]]
    fn.freq.obs.pred1(d=Fit.M_H[[s]],
                      d.KW=dummy,
                      net.names=net.names)
    #if(s==nn[length(nn)])
    if(s==8)
    {
      Nms=names(Cols.type)
      Nms=ifelse(Nms=="norm.loc","Normal fixed",
                 ifelse(Nms=="norm.sca","Normal prop.",
                        ifelse(Nms=="gamma","Gamma",
                               ifelse(Nms=="lognorm","Lognormal",NA))))
      legend("left",Nms,bty='n',text.col=Cols.type,cex=1)
    }
    
    mtext(names(Fit.M_H)[s],4,line=.5,cex=.8,las=3)
  }
  mtext("Frequency", side = 2,outer=T, line = 0.2,las=3,cex=1.2)
  mtext("Total length (cm)",1,outer=T,line=.4,cex=1.2)
  dev.off()
  
  #Target species
  if(do.this)
  {
    nn=match(Published,names(Fit.M_H))
    tiff(file="Obs.vs.Pred_Target species.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(mfrow=c(length(nn),length(Nets)),mar=c(1.5,1.2,.2,.65),oma=c(1.5,2,.8,1.2),mgp=c(1,.5,0),las=1,xpd=T)
    for(s in nn)
    {
      dummy=NULL
      net.names=FALSE
      if(s==nn[1]) net.names=TRUE
      
      if(Do.K_W) dummy=Pred.sel.K_W[[s]]
      fn.freq.obs.pred1(d=Fit.M_H[[s]],
                        d.KW=dummy,
                        net.names=net.names)
      if(s==12)
      {
        Nms=names(Cols.type)
        Nms=ifelse(Nms=="norm.loc","Normal fixed",
                   ifelse(Nms=="norm.sca","Normal prop.",
                          ifelse(Nms=="gamma","Gamma",
                                 ifelse(Nms=="lognorm","Lognormal",NA))))
        legend("left",Nms,bty='n',text.col=Cols.type,cex=1.1)
      }
      
      mtext(names(Fit.M_H)[s],4,line=.5,cex=1,las=3)
    }
    mtext("Frequency", side = 2,outer=T, line = 0.2,las=3,cex=1.2)
    mtext("Total length (cm)",1,outer=T,line=.4,cex=1.2)
    dev.off()
    
  }
 
    #Families together
  tiff(file="Figure.S6_Families.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(length(n.sp.family),length(Nets)),mar=c(1.5,1.2,.2,.65),oma=c(1.5,2,.8,1.2),mgp=c(1,.5,0),las=1,xpd=T)
  for(s in 1:length(n.sp.family))
  {
    dummy=NULL
    net.names=FALSE
    if(s==1) net.names=TRUE
    
    fn.freq.obs.pred1(d=Fit.M_H.family[[s]],
                      d.KW=dummy,
                      net.names=net.names)
    if(s==length(n.sp.family))
    {
      Nms=names(Cols.type)
      Nms=ifelse(Nms=="norm.loc","Normal fixed",
                 ifelse(Nms=="norm.sca","Normal prop.",
                        ifelse(Nms=="gamma","Gamma",
                               ifelse(Nms=="lognorm","Lognormal",NA))))
      legend("left",Nms,bty='n',text.col=Cols.type,cex=1)
    }
    
    mtext(n.sp.family[s],4,line=.5,cex=.8,las=3)
  }
  mtext("Frequency", side = 2,outer=T, line = 0.2,las=3,cex=1.2)
  mtext("Total length (cm)",1,outer=T,line=.4,cex=1.2)
  dev.off()
  
  
 
  #12. Extract mode for each mesh
  Mode.normal=function(m,k) m*k
  Mode.gamma=function(m,k,alpha) (alpha-1)*k*m
  Mode.lognormal=function(m,mu,sigma,m1) exp(mu-sigma^2)*(m/m1)
  
  fn.get.mode=function(SP,d,best.fit,meshes)
  {
    dummy=data.frame(meshes)
    dummy$Mode=NA
    if(best.fit=="lognorm")
    {
      ParS=d$Equal.power$lognorm$gear.pars[,'estimate']
      for(i in 1:nrow(dummy)) dummy$Mode[i]=Mode.lognormal(m=dummy$meshes[i],
                                                           mu=ParS[1],
                                                           sigma=ParS[2],
                                                           m1=min(meshes))
    }
    if(best.fit%in%c("norm.loc","norm.sca"))
    {
      if(best.fit=="norm.loc") ParS=d$Equal.power$norm.loc$gear.pars[,'estimate']
      if(best.fit=="norm.sca") ParS=d$Equal.power$norm.sca$gear.pars[,'estimate']
      for(i in 1:nrow(dummy)) dummy$Mode[i]=Mode.normal(m=dummy$meshes[i],k=ParS[1])
    }
    if(best.fit=="gamma")
    {
      ParS=d$Equal.power$gamma$gear.pars[,'estimate']
      for(i in 1:nrow(dummy)) dummy$Mode[i]=Mode.gamma(m=dummy$meshes[i],
                                                       k=ParS[2],
                                                       alpha=ParS[1])
    }
    return(dummy%>%
            mutate(Species=SP)%>%
            spread(meshes,Mode))
  }

  #Species
  store.mode.species=vector('list',length(n.sp))
  names(store.mode.species)=n.sp
  for(s in 1:length(n.sp))
  {
    if(!names(Fit.M_H)[s]%in%Published.sel.pars_K.W$Species)
      store.mode.species[[s]]=fn.get.mode(SP=n.sp[s],
                                          d=Fit.M_H[[s]],
                                          best.fit=Best.fit[[s]]$Model,
                                          meshes=c(10.2,12.7,15.2,16.5,17.8,20.3,21.6))
  }
  write.csv(do.call(rbind,store.mode.species),"TableS2.modes.csv",row.names = F)
  
  #Family
  store.mode.family=vector('list',length(n.sp.family))
  names(store.mode.family)=n.sp.family
  for(s in 1:length(n.sp.family)) store.mode.family[[s]]=fn.get.mode(SP=n.sp.family[s],
                                                            d=Fit.M_H.family[[s]],
                                                            best.fit=Best.fit.family[[s]]$Model,
                                                            meshes=c(10.2,12.7,15.2,16.5,17.8,20.3,21.6,25.4))
  write.csv(do.call(rbind,store.mode.family),"TableS2.modes_family.csv",row.names = F)
  
  
  #13. Display best model selectivity 
  colfunc <- colorRampPalette(c("cadetblue2", "deepskyblue4")) #Colors for displaying mesh selectivity
  unik.mesh=sort(unique(Combined$Mesh.size))
  CLS=colfunc(length(unik.mesh))
  names(CLS)=unik.mesh
  
  show.this.mesh=c('15.2','16.5','17.8')
  CLS=CLS[match(show.this.mesh,names(CLS))]
  
  plot.sel=function(d,NM,XMIN=min(d$Size.class),XMAX)  
  {
    d=d%>%rename(Size.class=TL.mm)
    plot(d$Size.class,d$'15.2',type='l',lwd=2,col=CLS[1],ylab='',xlab='',
         xlim=c(XMIN,XMAX),ylim=c(0,1))
    lines(d$Size.class,d$'16.5',lwd=2,col=CLS[2],lty=2)
    lines(d$Size.class,d$'17.8',lwd=2,col=CLS[3],lty=3)
    mtext(NM,3)
    return(d)
  }
  
    #New species  
  Store.Sels=vector('list',length(n.sp))
  names(Store.Sels)=n.sp
  nn=match(n.sp.pub,names(Fit.M_H))
  tiff(file="Figure 4. Selectivity.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=length(nn),MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
  for(s in nn)
  {
    Store.Sels[[s]]=plot.sel(d=Pred.sels.main.mesh[[s]],
                             NM=n.sp[s],
                             XMIN=min(Combined%>%
                                        filter(Species%in%n.sp)%>%
                                        pull(Length)),
                             XMAX=max(Combined%>%
                                        filter(Species%in%n.sp)%>%
                                        pull(Length)))
  }
  legend("topright",paste(round(as.numeric(names(CLS)),2)),col=CLS,bty='n',
                    lwd=2,cex=1.25,lty=1:3,title='Mesh (cm)')
  mtext("Total length (cm)",1,outer=T,line=.35,cex=1.25)
  mtext("Relative selectivity",2,outer=T,line=1,cex=1.25,las=3)
  dev.off() 
  
  #Target species (already published)
  nn=match(Published,names(Pred.sels.main.mesh))
  tiff(file="Selectivity_Target species.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=length(nn),MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
  for(s in nn)
  {
    Store.Sels[[s]]=plot.sel(d=Pred.sels.main.mesh[[s]],
                             NM=n.sp[s],
                             XMAX=1.5*max(Combined%>%
                                            filter(Species==n.sp[s])%>%
                                            pull(Length)))
    if(s==nn[length(nn)])legend("right",paste(round(as.numeric(names(CLS)),2)),col=CLS,bty='n',
                                lwd=2,cex=1.25,title='Mesh (cm)')
  }
  mtext("Total length (cm)",1,outer=T,line=.35,cex=1.25)
  mtext("Relative selectivity",2,outer=T,line=1,cex=1.25,las=3)
  dev.off() 
  
  #Each species separately
  for(s in 1:length(n.sp))
  {
    tiff(file=paste("Each species/Selectivity/",n.sp[s],".tiff",sep=''),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    par(las=1)
    dummy=plot.sel(d=Pred.sels.main.mesh[[s]],
                   NM=n.sp[s],
                   XMAX=1.5*max(Combined%>%
                                  filter(Species==n.sp[s])%>%
                                  pull(Length)))
    legend("right",paste(round(as.numeric(names(CLS)),2)),col=CLS,bty='n',
           lwd=2,cex=1.25,title='Mesh (cm)')
    mtext("Total length (cm)",1,line=2.5,cex=1.5)
    mtext("Relative selectivity",2,line=2.5,cex=1.5,las=3)
    dev.off()
  }
  
    #family
  Store.Sels.fam=vector('list',length(n.sp.family))
  names(Store.Sels.fam)=n.sp.family
  tiff(file="Figure 5.family.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=length(n.sp.family),MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.5,3,.1,.1),
            MGP=c(1,.5,0))
  par(cex.axis=1.25)
  for(s in 1:length(n.sp.family))
  {
    Store.Sels.fam[[s]]=plot.sel(d=Pred.sels.main.mesh.family[[s]],
                                 NM=n.sp.family[s],
                                 XMAX=1.5*max(Combined.family%>%
                                                filter(Species==n.sp.family[s])%>%
                                                pull(Length)))
  }
  legend("right",paste(round(as.numeric(names(CLS)),2)),col=CLS,bty='n',
         lwd=2,lty=1:3,cex=1.25,title='Mesh (cm)')
  mtext("Total length (cm)",1,outer=T,line=.3,cex=1.25)
  mtext("Relative selectivity",2,outer=T,line=1,cex=1.25,las=3)
  dev.off() 
  

  #14. Plot each species' selectivity separately for K&W
  if(Do.K_W)
  {
    #Get predicted selectivity
    fn.plt.Sel=function(Dat,theta)
    {
      Theta1=exp(theta[1])
      Theta2=exp(theta[2])
      
      sizes=seq(10*round(min(Dat$Length)*.9/10),10*(round(max(Dat$Length)*1.1/10)),by=10)
      meshes=unique(Dat$Mesh.size)
      d=data.frame(Size.class=rep(sizes,length(meshes)))%>%
        mutate(Mesh.size=rep(meshes,each=length(sizes)),
               alpha.beta=Theta1*Mesh.size,
               beta=-0.5*((alpha.beta)-((alpha.beta*alpha.beta+4*Theta2)^0.5)),
               alpha=alpha.beta/beta,
               Rel.sel=((Size.class/(alpha*beta))^alpha)*(exp(alpha-(Size.class/beta))))
      
      S.ij=d%>%
        distinct(Size.class,Mesh.size,.keep_all = T)%>%
        dplyr::select(Size.class,Mesh.size,Rel.sel)%>%
        spread(Mesh.size,Rel.sel,fill = 0)
      row.names(S.ij)=S.ij$Size.class
      S.ij=S.ij[,-1]
      
      ggplot(d, aes(Size.class,  Rel.sel)) + geom_line(aes(colour = factor(Mesh.size)),size=1.5)
      ggsave(paste("Each species/Selectivity/K&W_",n.sp[s],'.tiff',sep=''), width = 8,height = 8, dpi = 300, compression = "lzw")
      
      return(S.ij)
    }
    for(s in 1:length(n.sp)) Combined.sel[[s]]=fn.plt.Sel(Dat=Combined%>%filter(Species==n.sp[s]),theta=Fit.K_W[[s]]$par)
    
    # Show Selectivity for published species
    tiff(file="Figure 1.Selectivity_K&W.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    smart.par(n.plots=length(n.sp.pub)+1,MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
    for(s in 1:length(n.sp.pub))
    {
      d=Combined.sel[[s]]%>%
        mutate(Length=as.numeric(rownames(Combined.sel[[s]])),
               TL=Length/10)
      
      plot(d$TL,d$TL,type='l',ylab='',xlab='',
           lwd=5,col='transparent',main=n.sp.pub[s],ylim=c(0,1))
      for(l in 1:(ncol(d)))
      {
        lines(d$TL,d[,l],col=CLS[match(names(d)[l],names(CLS))],lwd=2)
      }
    }
    plot.new()
    legend("center",paste(as.numeric(names(CLS))),col=CLS,bty='n',lwd=2,cex=1.25)
    mtext("Total length (cm)",1,outer=T,line=.35,cex=1.25)
    mtext("Relative selectivity",2,outer=T,las=3,line=1.5,cex=1.25)
    dev.off() 
    
    #Combined selectivity all species together for K$W
    tiff(file="Combined selectivity all selected species_K&W.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
    smart.par(n.plots=length(n.sp)+1,MAR=c(1.5,1,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
    for(s in 1:length(n.sp))
    {
      Sum.sel=rowSums(Combined.sel[[s]])
      Combined.sel[[s]]$combined=Sum.sel/max(Sum.sel)
      if(!is.na(sum(Combined.sel[[s]]$combined)))
      {
        plot(as.numeric(rownames(Combined.sel[[s]])),Combined.sel[[s]]$combined,type='l',ylab='',xlab='',
             lwd=5,col='orange',main=n.sp[s])
        for(l in 1:(ncol(Combined.sel[[s]])-1)) lines(as.numeric(rownames(Combined.sel[[s]])),Combined.sel[[s]][,l])
      }
     }
    plot.new()
    legend('center',"combined",lwd=4,col="orange",bty='n',cex=1.25)
    mtext("Total length (mm)",1,outer=T,line=.35)
    mtext("Relative selectivity",2,outer=T,las=3,line=1.75)
    dev.off() 
  }
}


# Re fit Rory's data for sandbar to test model set up  ------------------------------------------------------------------
refit.Rory=FALSE
if(refit.Rory)
{
  Meshsize=as.numeric(substr(names(Rory.d)[-1],2,10))
  these.dist=c("lognorm","gamma","norm.sca","norm.loc")
  Rory=vector('list',length(these.dist))
  names(Rory)=these.dist
  for(d in 1:length(Rory))Rory[[d]]=gillnetfit(data=as.matrix(Rory.d),
                                               meshsizes=Meshsize,
                                               type=these.dist[d],
                                               rel=Meshsize,
                                               plots=c(F,F),
                                               plotlens=Rory.d$Size.class,
                                               plotlens_age=LatAge$`Sandbar shark`$TL,
                                               details=T)
  
  Rory.published.lognormal=pred.lognormal(Rory.d$Size.class,
                                          m=16.5,
                                          m1=min(Meshsize),
                                          mu=4.285,
                                          sigma=0.415)
  Rory.published.K.and.W=pred.Kirkwood.Walker(theta=log(c(137.3,134200)),
                                              pred.len=Rory.d$Size.class,
                                              Mesh=16.5)
  Pars=Rory$lognorm$gear.pars
  Rory.re_fitted.lognormal=pred.lognormal(Rory.d$Size.class,
                                          m=16.5,
                                          m1=min(Meshsize),
                                          mu=Pars[grep('mu1',rownames(Pars)),1],
                                          sigma=Pars[grep('sigma',rownames(Pars)),1])
  plot(Rory.d$Size.class,Rory.re_fitted.lognormal)
  lines(Rory.d$Size.class,Rory.published.lognormal)
  lines(Rory.d$Size.class,Rory.published.K.and.W$'6.5',col=2)
}

# Compare published selectivities (K&W) with those estimated here-------------------------------------------------------------------------
do.comparison=FALSE
if(do.comparison)
{
  fn.sel=function(Length,Mesh,Theta1,Theta2)
  {
    d=data.frame(Mesh.size=Mesh,
                 Size.class=Length) 
    d=d%>%
      mutate(alpha.beta=Theta1*Mesh.size,
             beta=-0.5*((alpha.beta)-((alpha.beta*alpha.beta+4*Theta2)^0.5)),
             alpha=alpha.beta/beta,
             Rel.sel=((Size.class/(alpha*beta))^alpha)*(exp(alpha-(Size.class/beta))))%>%
      dplyr::select(Mesh.size,Size.class,Rel.sel)
  }
  These.meshes=c(6,6.5,7)
  Published.sel_K.W=vector('list',length(These.meshes))
  names(Published.sel_K.W)=These.meshes
  for(m in 1:length(These.meshes))
  {
    dummy=vector('list',length(Published))
    names(dummy)=Published
    
    for(i in 1:length(Published))
    {
      dummy[[i]]=fn.sel(Length=seq(400,1600,by=10),  #in mm 
                        Mesh=These.meshes[m],   #in inches
                        Theta1=Published.sel.pars_K.W$Theta1[i],
                        Theta2=Published.sel.pars_K.W$Theta2[i])
    }
    Published.sel_K.W[[m]]=dummy
  }
  
  le.cols=c('red','forestgreen','blue4')
  names(le.cols)=c(15.2,16.5,17.8)
  Mesh.equiv=data.frame(cm=c(15.2,16.5,17.8),inches=c('6','6.5','7'))
  tiff(file="Compare my fit with published.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=length(Published),MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
  for(s in 1:length(Published))
  {
    #Plot estimated selectivity
    d=Store.Sels[[match(Published[s],names(Store.Sels))]]
    thiss=match(c(15.2,16.5,17.8),colnames(d))
    thiss=thiss[!is.na(thiss)]
    d=d[,c(1,thiss)]
    plot(d[,1],d[,2],type='l',col=le.cols[match(colnames(d)[2],names(le.cols))],ylab='',xlab='',
         lwd=2,xlim=c(min(d[,1]),200))
    for(i in 3:ncol(d)) lines(d[,1],d[,i],lwd=2,col=le.cols[match(colnames(d)[i],names(le.cols))])
    
    #Add published K&W selectivity
    id=Mesh.equiv$inches[match(colnames(d)[-1],Mesh.equiv$cm)]
    id.col=Mesh.equiv$cm[match(colnames(d)[-1],Mesh.equiv$cm)]
    for(i in 1:length(id))
    {
      d=Published.sel_K.W[[match(id[i],names(Published.sel_K.W))]][[match(Published[s],names(Published.sel_K.W[[1]]))]]
      lines(d$Size.class/10,d$Rel.sel,lwd=2,col=le.cols[match(id.col[i],names(le.cols))],lty=3)
    }
    mtext(Published[s],3,cex=1.5)
  }
  mtext("Total length (cm)",1,outer=T,line=.15,cex=1.25)
  mtext("Relative selectivity",2,outer=T,line=1,cex=1.25,las=3)
  legend('topright',names(le.cols),col=le.cols,bty='n',lwd=2,cex=1.5,title = 'Mesh (cm)')
  dev.off() 
}
