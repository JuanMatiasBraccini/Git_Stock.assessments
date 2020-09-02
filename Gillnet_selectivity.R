# Script for estimating gear selectivity based on Kirkwood & Walker 1986, extended to 
#       Millar & Holst 1997, and Millar & Fryer 1999

# size classes for selectivity analysis: 5 cm bins (as used in pop dyn model)
# size type for selectivity analysis: total length (in mm)

#assumptions: different mesh sizes set at the same time in same place
#             for greynurse, Shortfin mako, spinner and tiger shark, use commercial and mention in text the caveats
#             and future research to sort this out. Do selectivity tests.

#MISSING: 
#         For species where best is 'proportial power' model not going to 1 for all meshes... consider only 'equal power'?
#         Explore sandbar second peak for large mesh (see if fishing deeper, see Density all selected species_data set combined.tiff)
#         Squalidae is bi modal (2 different species??)
#         Compare published selectivities with those estimated here




library(tidyverse)
library(RODBC)
library(doParallel)
library(Hmisc)
library(magrittr)
library(expandFunctions)
library(stringr)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240) 
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")
source("C:/Matias/Analyses/Population dynamics/Git_Stock.assessments/NextGeneration.R")
source("C:/Matias/Analyses/Population dynamics/Git_Stock.assessments/SelnCurveDefinitions.R") #These can be extended by the user


Do.K_W=FALSE   #some issues with data in par estimation (e.g. gummy, smooth HH)

# DATA  -------------------------------------------------------------------

#1. WA Fisheries experimental mesh selectivity studies 

  #1.1 1994-1996      (need to physically open this file for R to connect)
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
channel <- odbcConnectExcel2007("C:/Matias/Data/SSF_survey_07_08/SharkSurveyData_30_09_2008.xls")
F2_Sampling<- sqlFetch(channel,"F2_Sampling", colnames = F)
F1_SamplingTwo<- sqlFetch(channel,"F1_SamplingTwo", colnames = F)
close(channel)


#3. TDGDLF observed catch composition 
#note: this is for fishers using one type of net at a time.
#      Not used as models cannot converge and not used simultaneously.
#      Mean sizes of 6.5 and 7 inch similar, or 6.5 > 7 
LFQ.south=read.csv("C:/Matias/Analyses/Selectivity_Gillnet/out.LFQ.south.csv")


#4. Species names
SP.names=read.csv('C:/Matias/Data/Species_names_shark.only.csv')
SP.codes=read.csv('C:/Matias/Data/Species.code.csv')

#5. Life history
LH=read.csv('C:/Matias/Data/Life history parameters/Life_History.csv')

# PARAMETERS  -------------------------------------------------------------------

Min.sample=25  #keep species with at least 20 observations in at least 2 mesh sizes
Min.nets=2
min.obs.per.mesh=Min.sample  #for each kept species, use nets with a minimum # of records

Size.Interval=50 #size intervals for selectivity estimation (in mm)

Min.length=250  #min and max length considered
Max.length=4000

do.paper.figures=FALSE
Preliminary=FALSE

Published=capitalize(tolower(c("Sandbar shark","Gummy Shark","Whiskery shark","Dusky shark")))


# PROCEDURE  -------------------------------------------------------------------

#--1. Data manipulation

  #1.1. Rory's

  #1.1.1 1994-1996
EXP_NET=EXP_NET[grep("E",EXP_NET$SHEET_NO),]%>%
           mutate(BOAT=VESSEL)%>%
          mutate(SKIPPER=NA,
                 Method='GN',
                 BOTDEPTH=NA,
                 SOAK.TIME=NA,
                 NET_LENGTH=NA,
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
          mutate(species=SPP_CODE)

Exp.net.94_96=EXPNET_B%>%
                left_join(SP.names,by=c("SPP_CODE" = "SPECIES"))%>%
                filter(!is.na(Name))%>%
                rename(tl=TOT_LENGTH,
                       fl=FORK_LNGTH)
colnames(Exp.net.94_96)=tolower(colnames(Exp.net.94_96))



  #1.1.2 2001-2003
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

Exp.net.94_96$experiment='94_96'
Exp.net.01_03$experiment='01_03'

This.col=c('sheet_no','date','experiment','mid.lat','mid.long','mesh_size','mesh_drop','name','tl','fl','sex')

Exp.net.WA=rbind(Exp.net.94_96[,match(This.col,names(Exp.net.94_96))],
                 Exp.net.01_03[,match(This.col,names(Exp.net.01_03))])%>%
          mutate(name=ifelse(name=="Angel Shark (general)","Angel Shark",
                      ifelse(name=="Eagle ray","Eagle Ray",
                      ifelse(name=="Gummy shark","Gummy Shark",
                      ifelse(name=="Port Jackson","PortJackson shark",
                      ifelse(name=="Big eye sixgill shark","Sixgill shark",
                      ifelse(name%in%c("Wobbegong (general)",'Spotted Wobbegong','Western Wobbegong'),"Wobbegongs",name)))))))

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


#Conversion FL to TL for records with no TL
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
                   Length=tl*10)%>%   #Length in mm; use Total length for selectivity estimation
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
  ggsave('C:/Matias/Analyses/Selectivity_Gillnet/Size.frequency_experimental.WA.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")
}


  #1.2 SSF
F2_Sampling=F2_Sampling%>%
                filter(Csiro>37000002 & Csiro<37039000 & LengthType=="TL" &
                      !is.na(Mesh1) & !is.na(Length))%>%
                mutate(Length=Length/10,              
                       Length=ifelse(Length>350,Length/10,Length),
                       Mesh1=ifelse(Mesh1%in%c('C6.00', 'C6.01', 'C6.02', 'C6.03', 'C6.04'),'C6',
                                    ifelse(Mesh1%in%c('C6.51', 'C6.52', 'C6.53', 'C6.54'),'C6.5',Mesh1)),
                       Mesh.size=2.54*as.numeric(substr(Mesh1,2,10)))%>%   #inches to cm
                filter(Length>=30)%>%
                mutate(Length=Length*10)  #length in mm

F2_Sampling=F2_Sampling%>%
              dplyr::select(Species,Csiro,Sex,Mesh1,Mesh.size,Length)

if(Preliminary) ggplot(F2_Sampling, aes(x = Length/10)) +
                    geom_histogram(color = "grey30", fill ="salmon",binwidth=10) +
                    facet_grid(Species~Mesh.size, scales = "free")

#Size frequency
if(Preliminary)
{
  dummy=F2_Sampling%>%mutate(mesh_size=factor(Mesh.size))
  ggplot(dummy,aes(Length/10,fill=mesh_size))+
    geom_histogram(color="#e9ecef",alpha=0.6, binwidth = 5) +
    labs(fill="Mesh size")+
    xlab("Total length (cm)") +
    facet_wrap(vars(Species), scales = "free") 
  ggsave('C:/Matias/Analyses/Selectivity_Gillnet/Size.frequency_SSF.tiff', width = 12,height = 8, dpi = 300, compression = "lzw")
}


  #1.3 TDGLDF
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
                filter(name%in%c("Smooth hammerhead","Spurdogs","Tiger shark"))

add1=TL_FL_LFQ.south%>%filter(name=="Smooth hammerhead")
add1=rbind(add1,add1)%>%mutate(name=c("Great hammerhead","Scalloped hammerhead"))  

TL_FL_LFQ.south=rbind(TL_FL_LFQ.south,
                      add1,
                      data.frame(name="Copper shark",intercept=5.801,slope=1.181),
                      data.frame(name="Common sawshark",intercept=0,slope=1.2),
                      data.frame(name="Grey nurse shark",intercept=0,slope=1.2),
                      data.frame(name="Milk shark",intercept=0,slope=1.2),
                      data.frame(name="Shortfin mako",intercept=0,slope=1.127),
                      data.frame(name="Spinner shark",intercept=0,slope= 1.28))  

LFQ.south=LFQ.south%>%
          left_join(TL_FL_LFQ.south,by=c(Species='name'))%>%
          mutate(tl=intercept+FL*slope,
                 Length=tl*10)%>%   #Length in mm; use Total length
          filter(!is.na(Length))%>%
    dplyr::select(Species,Mesh.size,Length,Data.set)%>%
          filter(!Species=='Spurdogs')  #mixed bag of species

#Size frequency
if(Preliminary)
{
  dummy=LFQ.south%>%mutate(Mesh_size=as.factor(Mesh.size))
  ggplot(dummy,aes(Length/10,fill=Mesh_size))+
    geom_histogram(color="#e9ecef",alpha=0.6, binwidth = 5) +
    labs(fill="Mesh size")+
    xlab("Total length (cm)") +
    facet_wrap(vars(Species), scales = "free") 
  ggsave('C:/Matias/Analyses/Selectivity_Gillnet/Size.frequency_TDGDLF_observed.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")
}

#--2. Combine all data sets    
Exp.net.WA=Exp.net.WA%>%
  rename(Species=name)%>%
  mutate(Mesh.size=2.54*mesh_size,
         Data.set="WA")
F2_Sampling=F2_Sampling%>%
        mutate(Data.set="SSF",
               Species=ifelse(Species=='Bronze Whaler',"Copper shark",Species))%>%
        filter(!Mesh1%in%c("C6.5","C6"))


Combined=rbind(F2_Sampling%>%dplyr::select(Species,Mesh.size,Length,Data.set),
               Exp.net.WA%>%dplyr::select(Species,Mesh.size,Length,Data.set))%>%
                  mutate(Species=tolower(Species),
                         Species=ifelse(Species=='portjackson shark','port jackson shark',
                                 ifelse(Species%in%c('wobbegong','banded wobbegong','spotted wobbegong'),
                                        'wobbegongs',Species)),
                         Species=capitalize(Species),
                         Species=ifelse(Species=='Port jackson shark','Port Jackson shark',Species),
                         Mesh.size=round(Mesh.size,1),
                         Species=ifelse(Species=='Angel shark','Australian angelshark',Species))%>%
  filter(!is.na(Mesh.size))%>%
  filter(!Species%in%c('Spurdogs')) #remove Spurdogs from WA, could be several species

SP.names=SP.names%>%
  mutate(Name=ifelse(Name=='PortJackson shark','Port Jackson shark',
              ifelse(Name=='Angel Shark','Australian angelshark',
              ifelse(Name=='Sevengill shark','Broadnose shark',
              ifelse(Name=='Sixgill shark','Bluntnose sixgill shark',
              ifelse(Name=='Eagle Ray','Eagle ray',
              Name))))))

Families=SP.names%>%  
  dplyr::select(Name,Family)%>%
  rename(Species=Name)
  
Combined=Combined%>%left_join(Families,"Species")
Combined.family=Combined%>%dplyr::select(-Species)%>%rename(Species=Family)  
Tab.sp.fam=table(Combined$Species,Combined$Family)
Tab.sp.fam[Tab.sp.fam>0]=1
Tab.sp.fam=colSums(Tab.sp.fam)

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

#Analyse species with at least Min.sample and Min.nets
TAB=table(Combined$Species,Combined$Mesh.size)
TAB[TAB<Min.sample]=0
TAB[TAB>=Min.sample]=1
Combined=Combined%>%
  filter(Species%in%names(which(rowSums(TAB)>=Min.nets)) & Length<=Max.length)

TAB=table(Combined.family$Species,Combined.family$Mesh.size)
TAB[TAB<Min.sample]=0
TAB[TAB>=Min.sample]=1
Combined.family=Combined.family%>%
  filter(Species%in%names(which(rowSums(TAB)>=Min.nets)) & Length<=Max.length)

#remove family with single species
Combined.family=Combined.family%>%filter(!Species%in%names(Tab.sp.fam[Tab.sp.fam==1]))

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
  Drop=as.numeric(gsub("[^0-9.]", "",  names(which(id<min.obs.per.mesh))))
  if(length(Drop)>0)
  {
    d=d%>%filter(!Mesh.size%in%Drop)
    id=table(d$Mesh.size)
    if(length(id)==1) d=NULL
  }
  
  Combined.family=rbind(Combined.family,d)
  
}
n.sp.family=sort(unique(Combined.family$Species))


#remove mesh sizes with higher mean size that next largest mesh
  #species
drop.one.mesh=data.frame(Species=c("Port Jackson shark","Sandbar shark","Whiskery shark"),
                         Mesh.size=c(10.2,10.2,16.5),
                         Drop="YES")
Combined=left_join(Combined,drop.one.mesh,by=c('Species','Mesh.size'))%>%
         filter(is.na(Drop))%>%dplyr::select(-Drop)

  #family
drop.one.mesh=data.frame(Species=c("Carcharhinidae"),
                         Mesh.size=c(10.2),
                         Drop="YES")
Combined.family=left_join(Combined.family,drop.one.mesh,by=c('Species','Mesh.size'))%>%
        filter(is.na(Drop))%>%dplyr::select(-Drop)


#Get length at age 
LatAge=vector('list',length(n.sp))
names(LatAge)=n.sp
LatAge.family=vector('list',length(n.sp.family))
names(LatAge.family)=n.sp.family

len.at.age=function(Lo,Linf,k,Age.max)
{
  VonB=data.frame(Age=0:Age.max)
  VonB$TL=10*(Lo+(Linf-Lo)*(1-exp(-k*VonB$Age))) 
  return(VonB)
}

LH=LH%>%left_join(SP.names,by='SPECIES')
  

for(s in 1:length(n.sp)) 
{
  ii=n.sp[s]
  ii=ifelse(ii=="Southern sawshark","Common sawshark",
     ifelse(ii=="Spikey dogfish","Spurdogs",ii))
  this.par=LH%>%filter(Name==ii) 
  if(is.na(this.par$Max_Age_max)) this.par$Max_Age_max=this.par$Max_Age*1.3
  LatAge[[s]]=with(this.par,len.at.age(Lo=LF_o,Linf=FL_inf/.85,k=K,Age.max=Max_Age_max))
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
              Max_Age_max=mean(Max_Age_max,na.rm=T))
  LatAge.family[[s]]=with(this.par,len.at.age(Lo=LF_o,Linf=FL_inf/.85,k=K,Age.max=Max_Age_max))
}  


#--3. Estimate selectivity parameters 
#note: shark size in mm

  #3.1 Millar & Holst 1997 
Fitfunction='gillnetfit'
#Fitfunction='NetFit'  #not used because it doesn't have gamma implemented
PlotLens=seq(Min.length+Size.Interval/2,Max.length-Size.Interval/2,by=Size.Interval)  #midpoints, 50 mm intervals  
Rtype=c("norm.loc","norm.sca","gamma","lognorm")   #consider this selection curves: normal fixed spread, 
                                                   #       normal spread proportional to mesh size,
                                                   #       gamma, lognormal
#Fit functions
Millar.Holst=function(d,size.int,length.at.age)
{
  #Create size bins
  d=d%>%mutate(Size.class=size.int*floor(Length/size.int)+size.int/2)
  
  #Tabulate observations by mid size class and mesh size   
  tab=d%>%
    group_by(Mesh.size,Size.class)%>%
    summarise(n=n())%>%
    spread(Mesh.size,n,fill=0)%>%
    data.frame
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
      reset.warnings()
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
      reset.warnings()
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
  Fit.M_H[[s]]=Millar.Holst(d=Combined%>%filter(Species==n.sp[s]),
                            size.int=Size.Interval,length.at.age=LatAge[[s]]$TL)
}

  #family
Fit.M_H.family=vector('list',length(n.sp.family))
names(Fit.M_H.family)=n.sp.family
for(s in 1:length(n.sp.family))  
{
  Fit.M_H.family[[s]]=Millar.Holst(d=Combined.family%>%filter(Species==n.sp.family[s]),
                            size.int=Size.Interval,length.at.age=LatAge.family[[s]]$TL)
}


if(Do.K_W)
{
  #3.2 Kirkwood & Walker
  Selectivty.Kirkwood.Walker=function(d,size.int,theta)
  {
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
  
  pred.Kirkwood.Walker=function(theta,pred.len,Mesh)
  {
    Theta1=exp(theta[1])
    Theta2=exp(theta[2])
    if(!16.5%in%Mesh) Mesh=c(Mesh,16.5)
    if(!17.8%in%Mesh) Mesh=c(Mesh,17.8)
    Mesh=sort(Mesh)
    
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
  
  theta.list=vector('list',length(n.sp))
  names(theta.list)=n.sp
  Fit.K_W=Pred.sel.K_W=theta.list
  Pred.sel.K_W_len.at.age=Pred.sel.K_W
  Pred.sel.K_W.family=vector('list',length(n.sp.family))
  names(Pred.sel.K_W.family)=n.sp.family
  Pred.sel.K_W.family_len.at.age=Fit.K_W.family=Pred.sel.K_W.family
  
  # initial parameter values
  theta.list$`Gummy shark`=c(Theta1=log(180),Theta2=log(29000))
  for(s in 1:length(theta.list)) theta.list[[s]]=jitter(theta.list$`Gummy shark`,factor=.25)
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
    Fit.K_W[[s]]=nlminb(theta.list[[s]], fn_ob, gradient = NULL)
    
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



  #3.3. Select best fit  

  #Species
Best.fit=vector('list',length(n.sp))
names(Best.fit)=n.sp
for(s in 1:length(n.sp))
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
  #  Tab3=rbind(Tab1,Tab2,data.frame(Model="K&W",Dev=K.and.W_Dev$Deviance[s],Fishing.power='')) #removed K&W as only Angel shark selected but poor fit
  Best.fit[[s]]=Tab3[which.min(Tab3[,2]),]
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
  Best.fit.family[[s]]=Tab3[which.min(Tab3[,2]),]
}


#ACA
# EXPORT SELECTIVITY OGIVES  ------------------------------------------------------------------

#function for predicting selectivity
pred.normal.fixed=function(l,k,m,sigma) exp(-((l-k*m)^2)/(2*(sigma)^2))
pred.normal.prop=function(l,m,a1,a2) exp(-(((l- a1*m)^2)/(2*a2*m^2)))
pred.gamma=function(l,m,k,alpha) ((l/((alpha-1)*k*m))^(alpha-1))*exp(alpha-1-(l/(k*m)))
pred.lognormal=function(l,m,m1,mu,sigma) (1/l)*exp(mu+(log(m/m1))-((sigma^2)/2)-((log(l)-mu-log(m/m1))^2)/(2*(sigma^2)))

#function for exporting
out.sel=function(d,BEST,NM,La)
{
  DAT=d$Equal.power     #set to 'equal power' so all meshes go to 1
  #if(BEST$Fishing.power=="Equal.power") DAT=d$Equal.power
  #if(BEST$Fishing.power=="Prop.power")  DAT=d$Prop.power
  
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
  write.csv(dat,paste('C:/Matias/Analyses/Data_outs/',NM,'/',NM,"_gillnet.selectivity",".csv",sep=''),row.names = F)
  
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
  write.csv(dat,paste('C:/Matias/Analyses/Data_outs/',NM,'/',NM,"_gillnet.selectivity_len.age",".csv",sep=''),row.names = F)
}

  #species
for(s in 1:length(n.sp))
{
  out.sel(d=Fit.M_H[[s]],BEST=Best.fit[[s]],
          NM=n.sp[s],La=LatAge[[s]])
}

  #family
for(s in 1:length(n.sp.family))
{
  out.sel(d=Fit.M_H.family[[s]],BEST=Best.fit.family[[s]],
          NM=n.sp.family[s],La=LatAge.family[[s]])
}
  



# REPORT  ------------------------------------------------------------------
if(do.paper.figures)
{
  setwd('C:/Matias/Analyses/Selectivity_Gillnet')
  
  colfunc <- colorRampPalette(c("white", "yellow","orange",'brown2',"darkred"))
  
  #Appendix 1
  Appendix1=TL_FL%>%
    mutate(intercept=round(intercept,2),
           slope=round(slope,2),
           name=capitalize(tolower(name)))%>%
    arrange(name)%>%
    distinct(name,.keep_all = T)
  write.csv(Appendix1,'Appendix1.csv',row.names = F)
  
  #Table 1. All species observed
  write.csv(Table1,'Table1.csv',row.names = F)
  
  #Mean size of by mesh of selected species
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
    group_by(Mesh.size,Species)%>%
    summarise(Mean=round(mean(Length),1))%>%
    spread(Mesh.size,Mean)%>%
    data.frame
  colnames(Tab2.mean)[-1]=substr(colnames(Tab2.mean)[-1],2,10)
  Combined.family%>%
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
  ggsave('Mean size by mesh_family.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")
  
  
  #Select species without selectivity published
  n.sp.pub=n.sp[-match(Published,n.sp)]
  
  
  #Display size frequencies 
  fig2=function(d)
  {
    cols=colfunc(length(unique(d$MESH)))
    d %>%
      ggplot( aes(x=TL, fill=MESH)) +
      geom_histogram(binwidth = 10, alpha=0.8,colour='grey40',size=.1)  +
      labs(fill="Mesh size (cm)")+
      facet_wrap(vars(Species), scales = "free")+
      xlab("Total length (cm)")+ ylab("Frequency")+scale_fill_manual(values=cols)
  }

    #species in paper
  fig2(d=Combined%>%
         mutate(MESH=as.factor(Mesh.size),
                TL=Length/10)%>%
        filter(Species%in%n.sp.pub))
  ggsave('Figure 2.Size frequency.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")
  
      #all species
  fig2(d=Combined%>%
         mutate(MESH=as.factor(Mesh.size),
                TL=Length/10))
  ggsave('Size frequency all selected species_data set combined.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")

      #family
  fig2(d=Combined.family%>%
         mutate(MESH=as.factor(Mesh.size),
                TL=Length/10))
  ggsave('Size frequency family_data set combined.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")
  
  
  
  #Display density
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
                TL=Length/10))
  ggsave('Density all selected species_data set combined.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")

    #family
  fig.density(d=Combined.family%>%
                mutate(MESH=as.factor(Mesh.size),
                       TL=Length/10))
  ggsave('Density family_data set combined.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")


  #Output parameter estimates from all models 
  fn.rnd=function(x) sprintf(round(x,2), fmt = '%#.2f')
  
    #species
  Table.mod.fit=vector('list',length(n.sp))
  names(Table.mod.fit)=n.sp
  for(s in 1:length(n.sp))
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
  fn.word.table(WD=getwd(),TBL=Table.mod.fit,Doc.nm="Table2",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
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
  fn.word.table(WD=getwd(),TBL=Table.mod.fit.family,Doc.nm="Table2.family",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")

  
  #Plot residuals for all model
  
    #species
      #all together
  tiff(file="Figure.3_Fit.tiff",width = 1800, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(length(n.sp),length(Rtype)),mar=c(1.5,1.2,.1,.1),oma=c(1.75,3,.1,1),mgp=c(1,.5,0))
  for(s in 1:length(n.sp))
  {
     for(f in 1:length(Fit.M_H[[s]]$Equal.power))
    {
      with(Fit.M_H[[s]]$Equal.power[[f]],
      {
        MAIN=""
        if(s==1)MAIN=ifelse(type=="norm.loc","Normal (fixed spread)",
             ifelse(type=="norm.sca","Normal (prop. spread)",
             ifelse(type=="gamma","Gamma",
             ifelse(type=="lognorm","Lognormal",NA))))
        plot.resids(devres,meshsizes,lens,title=MAIN)
      })
     }
    mtext( names(Fit.M_H)[s],4)
  }
  mtext("Total length (mm)",1,outer=T,line=.35,cex=1.25)
  mtext("Mesh size (cm)",2,outer=T,line=.35,cex=1.25,las=3)
  dev.off()
  
      #each individual species
  for(s in 1:length(n.sp))
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

    #family
      #all together
  tiff(file="Figure.3_Fit.family.tiff",width = 1800, height = 2400,units = "px", res = 300, compression = "lzw")    
  par(mfrow=c(length(n.sp.family),length(Rtype)),mar=c(1.5,1.2,.1,.1),oma=c(1.75,3,.1,1),mgp=c(1,.5,0))
  for(s in 1:length(n.sp.family))
  {
    for(f in 1:length(Fit.M_H.family[[s]]$Equal.power))
    {
      with(Fit.M_H.family[[s]]$Equal.power[[f]],
           {
             MAIN=""
             if(s==1)MAIN=ifelse(type=="norm.loc","Normal (fixed spread)",
                                 ifelse(type=="norm.sca","Normal (prop. spread)",
                                        ifelse(type=="gamma","Gamma",
                                               ifelse(type=="lognorm","Lognormal",NA))))
             plot.resids(devres,meshsizes,lens,title=MAIN)
           })
    }
    mtext( names(Fit.M_H.family)[s],4)
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
  
  #Observed vs predicted number at size by mesh for Kirkwood & Walker
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

  #Plot Observed size frequency VS estimated selectivity
  Cols.type=c('brown','red','forestgreen','orange')
  names(Cols.type)=Rtype
  if(Do.K_W)
  {
    Cols.type=c(Cols.type,'steelblue')
    names(Cols.type)=c(Rtype,"K&W")
  }
  fn.freq.obs.pred=function(d,d.KW=NULL,NME,MAR,OMA)
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
    smart.par(n,MAR,OMA,MGP=c(1,.5,0))
    Msh=substr(names(d$tab)[-1],2,10)
    for(i in 1:n)
    {
      plot(d$tab$Size.class,d$tab[,i+1],type='h',ylab='',xlab='', lwd = 4,col=rgb(.5,.5,.5,alpha=.5))
      for(m in 1:length(N.fish.caught))
      {
        lines(d$tab$Size.class,N.fish.caught[[m]][,i],col=Cols.type[m],lwd=1.5)
      }
      Where=which.max(N.fish.caught[[1]][,i])
      if(Where<10) here='topright' else here= 'topleft'
      if(NME=="Sphyrnidae") here= 'topleft'
      if(i<n)legend(here,paste(Msh[i],'cm mesh'),bty='n',cex=1.2)
      
    }
    legend(here,names(Cols.type),bty='n',lty=1,lwd=2,col=Cols.type,title=paste(Msh[i],'cm mesh'),cex=1.2)
    mtext("Frequency", side = 2,outer=T, line = 0.2,las=3,cex=1.2)
    mtext("Total length (mm)",1,outer=T,line=.45,cex=1.2)
    mtext(NME,3,outer=T,line=-.75,cex=1.2)
  }
  
    #species
  for(s in 1:length(n.sp))
  {
    tiff(file=paste("Each species/Fit/species/Observed.vs.pred_",n.sp[s],".tiff",sep=''),width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
    dummy=NULL
    if(Do.K_W) dummy=  Pred.sel.K_W[[s]]
    fn.freq.obs.pred(d=Fit.M_H[[s]],
                     d.KW=dummy,
                     NME=names(Fit.M_H)[s],
                     MAR=c(1.2,1.75,1,1),
                     OMA=c(1.75,1.75,.75,.75))
    dev.off()
  }
  
    #family
  for(s in 1:length(n.sp.family))
  {
    tiff(file=paste("Each species/Fit/family/Observed.vs.pred_",n.sp.family[s],".tiff",sep=''),width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
    fn.freq.obs.pred(d=Fit.M_H.family[[s]],
                     NME=names(Fit.M_H.family)[s],
                     MAR=c(1.2,1.85,1,1.5),
                     OMA=c(1.75,1.5,.75,1.5))
    dev.off()
  }

  #Extract model for each mesh
  Mode.normal=function(m,k) m*k
  Mode.gamma=function(m,k,alpha) (alpha-1)*k*m
  Mode.lognormal=function(m,mu,sigma,m1) exp(mu-sigma^2)*(m/m1)
  
  #Example for gummy
  #Fit.M_H$`Gummy shark`$Equal.power$norm.sca$gear.pars
  #Mode.normal(m=16.5,k=74.89336)
  #Mode.gamma(m=16.5,k=1.571537,alpha=49.218688)
  #Mode.lognormal(m=16.5,mu=6.6617214,sigma=0.1460648,m1=10.2)
  
  
  #Display best model selectivity
  
  colfunc <- colorRampPalette(c("cadetblue2", "deepskyblue4")) #Colors for displaying mesh selectivity
  unik.mesh=sort(unique(Combined$Mesh.size))
  CLS=colfunc(length(unik.mesh))
  names(CLS)=unik.mesh
  
  plot.sel=function(d,BEST,NM)
  {
    ThiS=as.numeric(substr(names(d$tab)[-1],2,10))
    cl=CLS[match(ThiS,names(CLS))]
    
    if(BEST$Fishing.power=="Equal.power") DAT=d$Equal.power
    if(BEST$Fishing.power=="Prop.power")  DAT=d$Prop.power
    id=match(BEST$Model,names(DAT))
    
    with(DAT[[id]],
         {
           plot(plotlens,rselect[,1],type='l',lwd=2,col=cl[1],ylab='',xlab='',ylim=c(0,1))
           for(m in 2:length(cl)) lines(plotlens,rselect[,m],lwd=2,col=cl[m])
         })
     mtext(NM,3)
  }
  
    #species 
  tiff(file="Figure 3. Selectivity.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=length(n.sp)+1,MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
  for(s in 1:length(n.sp))
  {
    plot.sel(d=Fit.M_H[[s]],BEST=Best.fit[[s]],NM=n.sp[s])
  }
  plot.new()
  legend("center",paste(round(as.numeric(names(CLS)),2),'cm'),col=CLS,bty='n',lwd=2,cex=1.25,title='Mesh')
  mtext("Total length (mm)",1,outer=T,line=.15,cex=1.25)
  mtext("Relative selectivity",2,outer=T,line=1,cex=1.25,las=3)
  dev.off() 
  
    #family
  tiff(file="Figure 3. Selectivity.family.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=length(n.sp.family)+1,MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
  for(s in 1:length(n.sp.family))
  {
    plot.sel(d=Fit.M_H.family[[s]],BEST=Best.fit.family[[s]],NM=n.sp.family[s])
  }
  plot.new()
  legend("center",paste(round(as.numeric(names(CLS)),2),'cm'),col=CLS,bty='n',lwd=2,cex=1.25,title='Mesh')
  mtext("Total length (mm)",1,outer=T,line=.15,cex=1.25)
  mtext("Relative selectivity",2,outer=T,line=1,cex=1.25,las=3)
  dev.off() 
  

  #Plot each species' selectivity separately for K&W
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
