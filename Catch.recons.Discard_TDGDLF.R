#                 SCRIPT FOR RECONSTRUCTING TOTAL DISCARDS IN TDGDLF    #

#notes:
# observer data must have at least Min.obs.per.block & Min.shots.per.block
# Ratio estimator method used with all years combined.  
# Model-based by species not used due to small sample sizes for most species
# Uncertainty determined using bootstrapping (as bootstrap resamples the shots, 
#       the procedure is self-weighting, i.e more weight is given to more 
#       commonly sampled blocks)
# Rarely discarded species were grouped as 'other' and then multiplied by proportion
# Longlines: Very few blocks (4 observed blocks only) so few observations and too 
#           much extrapolation...Hence, grouped commercial LL and GN and used observed GN ratio

rm(list=ls(all=TRUE))

options(stringsAsFactors = FALSE)
library(rlang)
library(tidyverse)
library(doParallel)
library(zoo)
library(abind)
library(stringr)
library(Hmisc)
library(ggrepel)
library(patchwork)

User="Matias"
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R"))




# 1. Data ---------------------------------------------------------

# Observers data
Res.ves=c("HAM","HOU","NAT","FLIN","RV BREAKSEA","RV Gannet","RV GANNET","RV SNIPE 2")
Dat_obs=DATA %>%
          filter(Mid.Lat<=(-26) & !BOAT%in%Res.ves & Taxa%in%c('Elasmobranch','Teleost') & !COMMON_NAME=='WHALE')  %>% 
          dplyr::select(c(SHEET_NO,date,Month,year,BOAT,zone,BLOCK,SOAK.TIME,MESH_SIZE,
                   MESH_DROP,NET_LENGTH,Mid.Lat,Mid.Long,Method,SPECIES,Taxa,
                   COMMON_NAME,SCIENTIFIC_NAME,CAES_Code,TL,FL,Disc.width,Number)) 
Dat_obs.LL=Dat_obs%>%
              filter(Method=="LL")
Dat_obs=Dat_obs%>%
              filter(Method=="GN")%>%
              mutate(COMMON_NAME=ifelse(SPECIES=='SM' & is.na(Disc.width),'Stingrays',COMMON_NAME),  #one smooth ray with no size info set to 'Dasyatidae'
                     SCIENTIFIC_NAME=ifelse(SPECIES=='SM' & is.na(Disc.width),'Family Dasyatidae',SCIENTIFIC_NAME),
                     CAES_Code=ifelse(SPECIES=='SM' & is.na(Disc.width),35000,CAES_Code),
                     SPECIES=ifelse(SPECIES=='SM' & is.na(Disc.width),'SR',SPECIES))   


# Commercial catch data   
Dat_total=read.csv(handl_OneDrive('Analyses\\Data_outs\\Data.monthly.csv'),stringsAsFactors = F)


# Species names
All.species.names=read.csv(handl_OneDrive("Data/Species_names_shark.only.csv")) #for catch


#List of discarded/retained species
Comm.disc.sp=read.csv(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/SPECIES+PCS+FATE.csv"),stringsAsFactors = F)


#Length weight relationships
Len.wei=read.csv(handl_OneDrive("Data/Length_Weights/length.weights.csv"),stringsAsFactors = F)

#Weight ranges
Wei.range=read.csv(handl_OneDrive("Data/Length_Weights/Data.Ranges.csv"))
Wei.range.names=read.csv(handl_OneDrive("Data/Length_Weights/Species.names.csv"))


#Post capture mortality
PCM.north=read.csv(handl_OneDrive('Analyses/Reconstruction_catch_commercial/TableS1.PCM_North.csv'))
PCM.south=read.csv(handl_OneDrive('Analyses/Reconstruction_catch_commercial/TableS1.PCM_South.csv'))


#Spatial distribution of discarded species
elasmos_dist=read.csv(handl_OneDrive('Analyses/Reconstruction_total_bycatch_TDGDLF/elasmos_dist.csv'))
teleosts_dist=read.csv(handl_OneDrive('Analyses/Reconstruction_total_bycatch_TDGDLF/teleosts_dist.csv'))


# 2. Parameter ---------------------------------------------------------

Min.obs.per.block=10  #minimum number of observations per block for use in ratio estimator
Min.shots.per.block=5  #minimum number of shots per block for use in ratio estimator

do.explrtn="NO"
STRTA.obs.disc=c("BLOCK")       #aggregating strata for observer data
STRTA.reported.ret=c("BLOCK")   #aggregating strata for total landed retained catch

n.boot=1e3

Group_rare_criteria=0.02    #criteria for grouping rare species (proportion of catch)    

Commercial.sp=subset(Comm.disc.sp,FATE=="C" & NATURE%in%c("S","R","T"))$SPECIES
Discarded.sp=subset(Comm.disc.sp,FATE=="D" & NATURE%in%c("S","R","T"))$SPECIES

do.sensitivity=TRUE
do.paper=FALSE

Stingrays=35000:40000

# Manipulate species names ---------------------------------------------------------
setwd(handl_OneDrive('Analyses/Reconstruction_total_bycatch_TDGDLF'))
All.species.names=All.species.names%>%
                    mutate(Name=capitalize(tolower(Name)),
                           Name=case_when(Name=='Port jackson shark'~'Port Jackson shark',
                                          Name=='Eagle ray'~'Southern eagle ray',
                                          TRUE~Name))

# Manipulate commercial catch ---------------------------------------------------------
Dat_total=Dat_total %>% filter(LAT<=(-26) & Estuary=="NO") %>%
      mutate(YEAR=as.numeric(substr(FINYEAR,1,4)),
             BLOCK=BLOCKX,
             Catch=LIVEWT.c,
             BLOCK=ifelse(BLOCK>=96000,paste(floor(abs(LAT)),(floor(LONG)-100)*10,sep=""),BLOCK),
             BLOCK=substring(BLOCK,1,4)) %>%
      dplyr::select(-c(ZnID,MonthlyID,ZoneID,YEAR.c,blockxFC,CONDITN,
                Factor,RSCommonName,RSSpeciesId,LANDWT,LIVEWT,
                FisheryZone,FisheryCode,Landing.Port,BDAYS,licence,
                Factor.c,LIVEWT.orgnl,LIVEWT.c,VesselID,BlockAveID,AnnualVesselAveID,
                BlockID,Reporter,Sch.or.DogS,LIVEWT.reap,Tot.shk.livewt,
                Shark.other.livewt,Reporter.old,Spec.old,Sname.old,NETLEN.c))%>% 
    left_join(All.species.names,by='SPECIES')
Dat_total.LL=Dat_total %>% filter(METHOD=="LL")
Percent.ktch.ll=100*sum(Dat_total.LL$Catch)/sum(Dat_total$Catch)
Dat_total=Dat_total %>% filter(METHOD=="GN")  #other methods discarded are catch is negligible and not expected to have discards (e.g. handline)



  #keep only elasmobranch species
DATA_total=list(GN=Dat_total%>%filter(SPECIES<50000),
                LL=Dat_total.LL%>%filter(SPECIES<50000))  

#Proportion of observed LL
prop.obs.LL=length(unique(Dat_obs.LL$SHEET_NO))/(length(unique(Dat_obs.LL$SHEET_NO))+length(unique(Dat_obs$SHEET_NO)))


# Manipulate observer data  ---------------------------------------------------------
Dat_obs=Dat_obs%>%
      mutate(COMMON_NAME=ifelse(SPECIES=='PD','Spurdogs',COMMON_NAME),
             SCIENTIFIC_NAME=ifelse(SPECIES=='PD','Squalus spp.',SCIENTIFIC_NAME),
             CAES_Code=ifelse(SPECIES=='PD',20000,CAES_Code),
             SPECIES=ifelse(SPECIES=='PD','SD',SPECIES))%>%
      left_join(All.species.names%>%rename(Species=SPECIES),by=c('SPECIES'='SP'))%>%
  mutate(Name=ifelse(is.na(Name),COMMON_NAME,Name),
         Scien.nm=ifelse(is.na(Scien.nm),SCIENTIFIC_NAME,Scien.nm),
         Name=ifelse(Name=="Leatherjacket (general)","Leatherjackets",
              ifelse(Name=="Triggerfish (general)","Triggerfishes",
              ifelse(Name=="Scorpion fish (general)","Scorpion fish",
                     ifelse(Name=="Stripped marlin","Striped marlin",Name)))))




DATA_obs=list(GN=Dat_obs) #Use only observed GN as LL has very few observations

#define discarded/retained sp   
for(i in 1:length(DATA_obs)) DATA_obs[[i]]$Discarded=with(DATA_obs[[i]],
                                          ifelse(SPECIES%in%Commercial.sp,"Retained",
                                          ifelse(SPECIES%in%Discarded.sp,"Discarded",
                                                 NA)))
#  Remove unidentified sharks
for(i in 1:length(DATA_obs)) DATA_obs[[i]]=subset(DATA_obs[[i]],!SPECIES%in%c('XX','XX.T'))


# Group Angel sharks-------------------------------------------
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]=DATA_obs[[i]]%>%
    mutate(COMMON_NAME=case_when(SPECIES%in%c("AO","AU")~"Angel Sharks (general)", TRUE~COMMON_NAME),
           SCIENTIFIC_NAME=case_when(SPECIES%in%c("AO","AU")~"Family Squatinidae", TRUE~SCIENTIFIC_NAME),
           CAES_Code=case_when(SPECIES%in%c("AO","AU")~24900, TRUE~CAES_Code),
           Species=as.double(Species),
           Species=case_when(SPECIES%in%c("AO","AU")~24900, TRUE~Species),
           Name=case_when(SPECIES%in%c("AO","AU")~"Angel sharks", TRUE~Name),
           Scien.nm=case_when(SPECIES%in%c("AO","AU")~"Squatinidae", TRUE~Scien.nm),
           SPECIES=case_when(SPECIES%in%c("AO","AU")~"AA",TRUE~SPECIES))
}


# Show Decadal species composition for Commercial catch recons paper-------------------------------------------
do.this=FALSE  
if(do.this)
{
  #Observer data
  these.common=c(23000,8001,10001,18001,18022,26999,13000,35000,24900,17006,
                 20000,39001,18023,19004,17003,17001,18007,7001,18003)
  a=DATA_obs$GN%>%
    mutate(Name=ifelse(Species%in%c(23001,23002),"Sawsharks",
                       ifelse(Species==18014,'Blacktips',
                              ifelse(Species%in%13001:13022,'Wobbegongs',Name))),
           Species=ifelse(Species%in%c(23001,23002),23000,
                          ifelse(Species==18014,18014,
                                 ifelse(Species%in%13001:13022,13000,Species))),                       
           Name1=ifelse(Taxa=='Teleost',"Teleosts",
                        ifelse(Taxa== 'Elasmobranch'& !Species%in%these.common,'Other sharks',
                               Name)),
           Decade=floor(year / 10) * 10)
  a%>%
    filter(Taxa== 'Elasmobranch')%>%
    filter(!(Decade==2020 & zone=="Zone2"))%>%
    filter(!(Decade==2020 & zone=="West"))%>%
    group_by(Name1,Decade,zone)%>%
    tally()%>%
    ggplot(aes(x=Decade, y=n, fill=Name1)) + 
    geom_bar(position="fill", stat="identity")+
    facet_wrap(~zone, scales = "free") +
    theme(legend.position="top",
          strip.text.x = element_text(size =11.5),
          axis.title=element_text(size=16)) +
    xlab("Decade") +
    ylab("Proportion")+
    labs(fill = "") +
    facet_wrap(~zone)
  ggsave(handl_OneDrive('Analyses/Reconstruction_catch_commercial/FigureS4.tiff'),width = 10,height = 6,compression = "lzw")
  
  
  #Commercial catch
  these.common=c(23000,10001,18001,18022,13000,35000,24900,17006,
                 20000,39001,18023,19004,17003,17001,18007,7001,18003)
  a=DATA_total$GN%>%
    mutate(Name=ifelse(SPECIES%in%c(23001,23002),"Sawsharks",
                       ifelse(SPECIES==18014,'Blacktips',
                              ifelse(SPECIES%in%13001:13022,'Wobbegongs',Name))),
           SPECIES=ifelse(SPECIES%in%c(23001,23002),23000,
                          ifelse(SPECIES==18014,18014,
                                 ifelse(SPECIES%in%13001:13022,13000,SPECIES))),                       
           Name1=ifelse(!SPECIES%in%these.common,'Other sharks',Name),
           Decade=floor(YEAR / 10) * 10)%>%
    filter(SPECIES<50000)
  
  a%>%
    group_by(Name1,Decade,zone)%>%
    tally()%>%
    ggplot(aes(x=Decade, y=n, fill=Name1)) + 
    geom_bar(position="fill", stat="identity")+
    facet_wrap(~zone, scales = "free") +
    theme(legend.position="top",
          strip.text.x = element_text(size =11.5),
          axis.title=element_text(size=16)) +
    xlab("Decade") +
    ylab("Proportion")+
    labs(fill = "") +
    facet_wrap(~zone)
  ggsave(handl_OneDrive('Analyses/Reconstruction_catch_commercial/FigureS5.tiff'),width = 10,height = 6,compression = "lzw")
  
}

# Show overal observed discarded and retained species-------------------------------------------
colfunc <- colorRampPalette(c("brown4","brown1","orange","yellow"))
  #Elasmos
fun.horiz.bar=function(d,gear)
{
  n=nrow(d%>%
    filter(Taxa=="Elasmobranch" & Discarded=="Discarded")%>%
    distinct(SPECIES))
  d1=d%>%
    mutate(taxa=ifelse(Discarded=="Retained","Elasmobranch",Taxa),
           dummy=ifelse(Discarded=="Retained" & Taxa=="Elasmobranch","Retained elasmobranchs",
                 ifelse(Discarded=="Retained" & Taxa=="Teleost","Retained teleosts",
                 ifelse(Discarded=="Discarded" & Taxa=="Teleost","Discarded teleosts",
                 Name))),
           dummy=ifelse(dummy=="Eagle ray","Southern eagle ray",dummy))%>%
    group_by(dummy)%>%
    tally%>%
    mutate(Percent=100*n/sum(n),
           Label=ifelse(Percent>=.5 & Percent<.1,'<1',
                 ifelse(Percent<.1,'<0.1',
                 round(Percent,1))),
           Label=paste(Label,'%',sep=''),
           Percent.original=Percent,
           Percent=ifelse(Percent<0.1,0.1,Percent))%>%
    data.frame%>%
    arrange(-Percent)%>%
    mutate(dummy=reorder(dummy, Percent),
           grp=ifelse(dummy=='Retained elasmobranchs',"#228B22",
               ifelse(dummy=='Retained teleosts',"#11BC11",
               ifelse(!dummy%in%c('Retained elasmobranchs','Retained teleosts') & Percent.original>=5,"#8B2323",
               ifelse(!dummy%in%c('Retained elasmobranchs','Retained teleosts') & Percent.original<5 & Percent.original>=1,"#D03434",
               ifelse(!dummy%in%c('Retained elasmobranchs','Retained teleosts') & Percent.original<1 & Percent.original>0.1,"#FF4A39",
                      "#FF9A06"))))))
  
  
  d1%>%
  ggplot(aes(x =  dummy,y = Percent)) + 
    geom_bar(fill = d1$grp,stat = "identity")+coord_flip()+
    geom_text(aes(label=Label), position=position_dodge(width=0.9), hjust=-.1)+
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ylim(0, 10*round(max(d1$Percent)/10)+10)
  ggsave(paste('Results/Figure2_observed.percentage_',gear,'.tiff',sep=''), 
         width = 8,height = 10,compression = "lzw")
  
  d1=d1%>%mutate(Label=ifelse(Percent>1,round(Percent,2),'<1%'))
  return(d1)

}
for(i in 1:length(DATA_obs))  Tabl.obs.GN=fun.horiz.bar(d=DATA_obs[[i]],gear=names(DATA_obs)[i])

  #Teleosts
fun.horiz.bar.scalefish=function(d,gear)
{
  n=nrow(d%>%
           filter(Taxa=="Teleost" & Discarded=="Discarded")%>%
           distinct(SPECIES))
  d1=d%>%
    mutate(taxa=ifelse(Discarded=="Retained","Elasmobranch",Taxa),
           dummy=ifelse(Discarded=="Retained" & Taxa=="Elasmobranch","Retained elasmobranchs",
                        ifelse(Discarded=="Retained" & Taxa=="Teleost","Retained teleosts",
                               ifelse(Discarded=="Discarded" & Taxa=="Elasmobranch","Discarded elasmobranchs",
                                      Name))),
           dummy=capitalize(tolower(dummy)),
           dummy=ifelse(dummy=='Sergeant baker','Sergeant Baker',dummy))%>%
    group_by(dummy)%>%
    tally%>%
    mutate(Percent=100*n/sum(n),
           Label=ifelse(Percent>=.5 & Percent<.1,'<1',
                        ifelse(Percent<.1,'<0.1',
                               round(Percent,1))),
           Label=paste(Label,'%',sep=''),
           Percent.original=Percent,
           Percent=ifelse(Percent<0.1,0.1,Percent))%>%
    data.frame%>%
    arrange(-Percent)%>%
    mutate(dummy=reorder(dummy, Percent),
           grp=ifelse(dummy=='Retained elasmobranchs',"#228B22",
                      ifelse(dummy=='Retained teleosts',"#11BC11",
                             ifelse(!dummy%in%c('Retained elasmobranchs','Retained teleosts') & Percent.original>=5,"#8B2323",
                                    ifelse(!dummy%in%c('Retained elasmobranchs','Retained teleosts') & Percent.original<5 & Percent.original>=1,"#D03434",
                                           ifelse(!dummy%in%c('Retained elasmobranchs','Retained teleosts') & Percent.original<1 & Percent.original>0.1,"#FF4A39",
                                                  "#FF9A06"))))))
  
  
  d1%>%
    ggplot(aes(x =  dummy,y = Percent)) + 
    geom_bar(fill = d1$grp,stat = "identity")+coord_flip()+
    geom_text(aes(label=Label), position=position_dodge(width=0.9), hjust=-.1)+
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ylim(0, 10*round(max(d1$Percent)/10)+10)
  ggsave(paste('Results/Recons.scalefish/Figure2_observed.percentage_',gear,'.tiff',sep=''), 
         width = 8,height = 10,compression = "lzw")
  
  d1=d1%>%mutate(Label=ifelse(Percent>1,round(Percent,2),'<1%'))
  return(d1)
  
}
for(i in 1:length(DATA_obs))  Tabl.obs.GN.scalefish=fun.horiz.bar.scalefish(d=DATA_obs[[i]],gear=names(DATA_obs)[i])

#Table of discarded and observed species
# Table.disc.obs=DATA_obs$GN%>%
#               distinct(SPECIES,.keep_all = T)%>%
#               dplyr::select(Name,Scien.nm,Discarded,Taxa)%>%
#               rename(Fate=Discarded)%>%
#               arrange(Taxa,Fate,Name)
Table.disc.obs=DATA_obs$GN%>%
              group_by(COMMON_NAME,SCIENTIFIC_NAME,Taxa,Discarded)%>%
              tally()%>%
              arrange(Taxa,Discarded,COMMON_NAME)  
  
#for teleosts, keep discarded teleosts and retained elasmos
DATA_obs_teleosts=DATA_obs
for(i in 1:length(DATA_obs_teleosts))
{
  DATA_obs_teleosts[[i]]=subset(DATA_obs_teleosts[[i]],
                                (Taxa=='Elasmobranch' & Discarded=='Retained') | 
                                (Taxa=='Teleost' & Discarded=='Discarded'))
}

#Keep only elasmobranchs from observed commercial boats
for(i in 1:length(DATA_obs)) DATA_obs[[i]]=subset(DATA_obs[[i]],Taxa=='Elasmobranch')



#Extract discarded species
Discarded.SP=DATA_obs$GN%>%
              filter(Discarded=='Discarded')%>%
              distinct(SPECIES,.keep_all = T)%>%
              pull(SPECIES)

Discarded.SP_teleosts=DATA_obs_teleosts$GN%>%
              filter(Discarded=='Discarded' & Taxa=='Teleost')%>%
              distinct(SPECIES,.keep_all = T)%>%
              pull(SPECIES)

# Set FL to disc width for stingrays for computational purposes----------------------------------------
Stingrays=35000:40000
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]=DATA_obs[[i]]%>%
    mutate(FL=ifelse(CAES_Code%in%Stingrays,Disc.width,FL),
           TL=ifelse(CAES_Code%in%Stingrays,NA,TL))
}

# Fill in missing FL info for elasmos ---------------------------------------------------------
#note: first use TL, then sample from species distribution, finally overall mean
UniK.sp=table(DATA_obs$GN$SPECIES)
UniK.sp=UniK.sp[UniK.sp>3]
UniK.sp=names(UniK.sp)
pdf("Results/Preliminary/size.frequency_observed.pdf")
for(u in 1:length(UniK.sp))
{
  a=subset(DATA_obs$GN,SPECIES==UniK.sp[u] & !is.na(FL))
  if(nrow(a)>2)hist(a$FL,main=a$Name[1],col=2, xlab="FL (cm)")
}
dev.off()

  #1. Set to NA FL or TL records less than size at birth
size.birth=30 #FL, in cm
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]$FL=with(DATA_obs[[i]],ifelse(FL<size.birth,NA,FL))
  DATA_obs[[i]]$TL=with(DATA_obs[[i]],ifelse(TL<(size.birth/.85),NA,TL))
  
  DATA_obs_teleosts[[i]]$FL=with(DATA_obs_teleosts[[i]],ifelse(FL<size.birth & Taxa=='Elasmobranch',NA,FL))
  DATA_obs_teleosts[[i]]$TL=with(DATA_obs_teleosts[[i]],ifelse(TL<(size.birth/.85) & Taxa=='Elasmobranch',NA,TL))
}
  
  #2. Derive FL as a proportion of TL
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]$FL=with(DATA_obs[[i]],ifelse(is.na(FL),TL*.875,FL))
  DATA_obs_teleosts[[i]]$FL=with(DATA_obs_teleosts[[i]],ifelse(is.na(FL) & Taxa=='Elasmobranch',TL*.875,FL))
}
  

  #3. samples from species distribution
fn.samp.dist=function(d)
{
  d=subset(d,FL>20)
  his=hist(d$FL,breaks=seq(floor(min(d$FL,na.rm=T)),ceiling(max(d$FL,na.rm=T)),by=1),plot=F)
  return(his)
}
for(i in 1:length(DATA_obs))
{
    #Elasmos
  all.sp=unique(DATA_obs[[i]]$SPECIES)
  for(s in 1:length(all.sp)) 
  {
    d=subset(DATA_obs[[i]],SPECIES==all.sp[s])

    if(sum(is.na(d$FL))>0 & sum(!is.na(d$FL))>0)
    {
      if(length(d$FL[!is.na(d$FL)])>2)Prob=fn.samp.dist(d)
      if(length(Prob$density)>2)
      {
        id=which (is.na(DATA_obs[[i]]$FL) & DATA_obs[[i]]$SPECIES==all.sp[s])
        DATA_obs[[i]]$FL[id]=sample(Prob$breaks[-1],
                                length(DATA_obs[[i]]$FL[id]),
                                replace=T,prob=Prob$density/sum(Prob$density))
      }
    }
    
    if(all.sp[s]=="SR")  #set stingrays equal to eagle ray
    {
      Prob=fn.samp.dist(d=subset(DATA_obs[[i]],SPECIES=='ER'))
      id=which (is.na(DATA_obs[[i]]$FL) & DATA_obs[[i]]$SPECIES==all.sp[s])
      DATA_obs[[i]]$FL[id]=sample(Prob$breaks[-1],
                                  length(DATA_obs[[i]]$FL[id]),
                                  replace=T,prob=Prob$density/sum(Prob$density))
      
    }
  }
  
    #Elasmos in Teleost data frame 
  all.sp=DATA_obs_teleosts[[i]]%>%filter(Taxa=='Elasmobranch')%>%distinct(SPECIES)%>%pull(SPECIES)
  for(s in 1:length(all.sp)) 
  {
    d=subset(DATA_obs_teleosts[[i]],SPECIES==all.sp[s])
    
    if(sum(is.na(d$FL))>0 & sum(!is.na(d$FL))>0)
    {
      if(length(d$FL[!is.na(d$FL)])>2)Prob=fn.samp.dist(d)
      if(length(Prob$density)>2)
      {
        id=which (is.na(DATA_obs_teleosts[[i]]$FL) & DATA_obs_teleosts[[i]]$SPECIES==all.sp[s])
        DATA_obs_teleosts[[i]]$FL[id]=sample(Prob$breaks[-1],
                                    length(DATA_obs_teleosts[[i]]$FL[id]),
                                    replace=T,prob=Prob$density/sum(Prob$density))
      }
    }
    

  }
}

  #4. overall mean (only 2 records, Crocodile shark & Grey reef shark) 
for(i in 1:length(DATA_obs))
{
  Mn.whaler=mean(DATA_obs[[i]]%>%filter(CAES_Code%in%18001:18036)%>%pull(FL),na.rm=T)
  DATA_obs[[i]]=DATA_obs[[i]]%>%
    mutate(FL=ifelse(is.na(FL) & COMMON_NAME=="Crocodile shark",Mn.whaler,
              ifelse(is.na(FL) & COMMON_NAME=="Grey reef shark",Mn.whaler,
              FL)))
  
  DATA_obs_teleosts[[i]]=DATA_obs_teleosts[[i]]%>%
    mutate(FL=ifelse(is.na(FL) & COMMON_NAME=="Crocodile shark",Mn.whaler,
                     ifelse(is.na(FL) & COMMON_NAME=="Grey reef shark",Mn.whaler,
                            FL)))
  }
  
  

pdf("Results/Preliminary/size.frequency_ammended.pdf")
for(u in 1:length(UniK.sp))
{
  a=subset(DATA_obs$GN,SPECIES==UniK.sp[u] & !is.na(FL))
  if(nrow(a)>2)hist(a$FL,main=a$Name[1],col=2, xlab="FL (cm)")
}
dev.off()


# Fix nonsense TL for teleost and fill in missing TL  ---------------------------------------------------------
for(i in 1:length(DATA_obs_teleosts))
{
  DATA_obs_teleosts[[i]]=DATA_obs_teleosts[[i]]%>%
    mutate(TL=ifelse(TL<10,NA,TL))
}

fn.samp.dist_teleost=function(d)
{
  his=hist(d$TL,breaks=seq(floor(min(d$TL,na.rm=T)),ceiling(max(d$TL,na.rm=T)),by=1),plot=F)
  return(his)
}

for(i in 1:length(DATA_obs_teleosts))
{
  for(s in 1:length(Discarded.SP_teleosts)) 
  {
    d=subset(DATA_obs_teleosts[[i]],SPECIES==Discarded.SP_teleosts[s])
    if(sum(is.na(d$TL))>0 & sum(!is.na(d$TL))>0)
    {
      if(length(d$TL[!is.na(d$TL)])>2)Prob=fn.samp.dist_teleost(d)
      if(length(Prob$density)>2)
      {
        id=which (is.na(DATA_obs_teleosts[[i]]$TL) & DATA_obs_teleosts[[i]]$SPECIES==Discarded.SP_teleosts[s])
        DATA_obs_teleosts[[i]]$TL[id]=sample(Prob$breaks[-1],
                                    length(DATA_obs_teleosts[[i]]$TL[id]),
                                    replace=T,prob=Prob$density/sum(Prob$density))
      }
    }
  }
}

for(i in 1:length(DATA_obs_teleosts))   #Marlin no TL data, only 1 record, set to tuna 90% percentile
{
  Mean.tuna=quantile(Dat_obs%>%
              filter(SPECIES=="SB.T" & TL>0)%>%
              pull(TL),probs=.9)
  DATA_obs_teleosts[[i]]=DATA_obs_teleosts[[i]]%>%
    mutate(TL=ifelse(SPECIES=='SM.T',Mean.tuna,TL))
  
}



# Convert numbers to weight  ---------------------------------------------------------
#note: calculate discard ratio using weight because total catch is in weight.
#       Use FL and a, b parameters for fork length
  #Elasmobranchs
Len.wei=Len.wei%>%
          dplyr::select(c(SPECIES,a_w8t,b_w8t,Source,Comments))
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]=DATA_obs[[i]]%>%
                 left_join(Len.wei,by="SPECIES")%>%
                 mutate(Catch=a_w8t*(FL)^b_w8t)

}

  #Teleosts  
for(i in 1:length(DATA_obs_teleosts))
{
  DATA_obs_teleosts[[i]]=DATA_obs_teleosts[[i]]%>%
                          left_join(Len.wei,by="SPECIES")%>%
                          mutate(TL=ifelse(SPECIES=="KF.T" & TL>25,25, #reset TL for some species (some too large observations)
                                    ifelse(SPECIES=="SC.T" & TL>42,42,
                                           TL)))  
  
  #Buff bream and North west blowfish length-weight is in FL, so convert TL to FL
  DATA_obs_teleosts[[i]]=DATA_obs_teleosts[[i]]%>%
                            mutate(dummy.FL=ifelse(Taxa=='Elasmobranch',FL,
                                            ifelse(Taxa=='Teleost',TL,
                                                   NA)),
                                   dummy.FL=ifelse(SPECIES%in%c("BB.T",'NW.T'),TL*.85,
                                                   dummy.FL),
                                    Catch=a_w8t*(dummy.FL)^b_w8t)%>%
                            dplyr::select(-dummy.FL)
  
}



# Remove white shark, greynurse and 'other shark' from observations ---------------------------------------------------------
#note: 2 very large white sharks acounts for a great proportion of discarded tonnage
for(i in 1:length(DATA_obs))  DATA_obs[[i]]=DATA_obs[[i]]%>%filter(!SPECIES%in%c('WP','GN','XX'))

#Reset min and maximum weights if nonsense
Wei.range=Wei.range%>%left_join(Wei.range.names,by=c("Sname"))
Wei.range=subset(Wei.range,!(is.na(SPECIES)|is.na(TW.min)|is.na(TW.max)))
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]=left_join(DATA_obs[[i]],Wei.range,by=c('Species'='SPECIES'))%>%
    mutate(Sname=case_when(SPECIES%in%c("AO","AU")~"Angel Shark (general)", TRUE~Sname))
  
  DATA_obs[[i]]=DATA_obs[[i]]%>%
                  mutate(Catch=case_when(!is.na(TW.max) & Catch> TW.max ~ TW.max,
                                         !is.na(TW.min) & Catch< TW.min ~ TW.min,
                                         TRUE~Catch))
}

pdf("Results/Recons.scalefish/Observed.Catch.frequency.pdf")
for(i in 1:length(DATA_obs_teleosts))
{
  for(s in 1:length(Discarded.SP_teleosts)) 
  {
    d=subset(DATA_obs_teleosts[[i]],SPECIES==Discarded.SP_teleosts[s])
    hist(d$Catch,xlab="Catch (kg)",main=d$Name[1],col=2)
  }
}
dev.off()


# Group rare species  ---------------------------------------------------------
  #Elasmos
for(i in 1:length(DATA_obs))
{
  Tab.sp=prop.table(with(subset(DATA_obs[[i]],Discarded=="Discarded"),table(SPECIES)))
  Rare.sp=names(Tab.sp[Tab.sp<Group_rare_criteria])
  DATA_obs[[i]]$SPECIES.ori=DATA_obs[[i]]$SPECIES
  DATA_obs[[i]]$SPECIES=with(DATA_obs[[i]],ifelse(SPECIES%in%Rare.sp,"Grouped",SPECIES))
}

  #Teleosts
for(i in 1:length(DATA_obs_teleosts))
{
  Tab.sp=prop.table(with(subset(DATA_obs_teleosts[[i]],Discarded=="Discarded" & Taxa=='Teleost'),table(SPECIES)))
  Rare.sp=names(Tab.sp[Tab.sp<Group_rare_criteria])
  DATA_obs_teleosts[[i]]$SPECIES.ori=DATA_obs_teleosts[[i]]$SPECIES
  DATA_obs_teleosts[[i]]$SPECIES=with(DATA_obs_teleosts[[i]],ifelse(SPECIES%in%Rare.sp,"Grouped",SPECIES))
}


# Define discarded_sp used to calculate discarded species proportions-------------------------------------------
  #Elasmos
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]$Discarded_sp=with(DATA_obs[[i]],ifelse(Discarded=="Retained","Retained",SPECIES))
}

  #Teleosts
for(i in 1:length(DATA_obs_teleosts))
{
  DATA_obs_teleosts[[i]]$Discarded_sp=with(DATA_obs_teleosts[[i]],ifelse(Discarded=="Retained","Retained",SPECIES))
}

# Data exploration-------------------------------------------
if(do.explrtn=="YES")
{
  A=with(Dat_obs[!duplicated(Dat_obs$SHEET_NO),],table(BLOCK))
  Exp=5*(A/max(A))
  Exp=ifelse(Exp<0.25,0.25,Exp)
  Exp=data.frame(LAT=-(as.numeric(substr(names(Exp),1,2))+.5),
                 Long=100+as.numeric(substr(names(Exp),3,4))+.5,cex=Exp)
  pdf("Results/Preliminary/shots_block_map.pdf")
  with(Exp,plot(Long,LAT,cex=cex,pch=19,col='steelblue'))
  legend('topright',paste(c(10,100,250)),pch=19,pt.cex=5*(c(10,100,250)/max(A)),
         bty='n',col='steelblue',title="# shots")
  dev.off()
  
  YRs=sort(unique(Dat_obs$year))
  pdf("Results/Preliminary/shots_block_year_map.pdf")
  par(mfcol=c(5,4),mar=c(1,1.5,1,.5),oma=c(3,1.25,.5,.1),las=1,mgp=c(1.5,.5,0),cex.axis=.85)
  for(y in 1:length(YRs))
  {
    x=subset(Dat_obs,year==YRs[y])
    A=with(x[!duplicated(x$SHEET_NO),],table(BLOCK))
    Exp=3.5*(A/max(A))
    Exp=ifelse(Exp<0.25,0.25,Exp)
    Exp=data.frame(LAT=-(as.numeric(substr(names(Exp),1,2))+.5),
                   Long=100+as.numeric(substr(names(Exp),3,4))+.5,cex=Exp)
    with(Exp,plot(Long,LAT,cex=cex,pch=19,col='steelblue',ylim=c(-36,-26),
                  xlim=c(113,129),ann=F))
    Quant=round(quantile(A,prob=c(.25,.5,.75)))
    legend('topright',paste(Quant),pch=19,pt.cex=3.5*(Quant/max(A)),
           bty='n',col='steelblue',title="# shots")
    legend('topleft',paste(YRs[y]),bty='n')
  }
  dev.off()
  
  #blocks by year
  pdf("Results/Preliminary/shots_block_year_plot.pdf")
  A=with(Dat_obs[!duplicated(Dat_obs$SHEET_NO),],table(year,BLOCK))
  plot(1,1,col="transparent",ylim=c(1,ncol(A)),xlim=c(1,nrow(A)),ylab="block",xlab="year")
  for(i in 1:nrow(A))
  {
    Exp=5*(A[i,]/max(A))
    points(rep(i,ncol(A)),1:ncol(A),cex=Exp,pch=19,col='steelblue')
  }
  legend('bottomright',paste(c(10,30,80)),pch=19,pt.cex=5*(c(10,30,80)/max(A)),
         bty='n',col='steelblue',title="# shots")
  dev.off()
  
  #Shark and ray species
  plot(sort(table(Dat_obs$SPECIES)),type='h',ylab="# individuals")
}


# Get dataframe for infographic-------------------------------------------
enough.data.yrs=Dat_obs%>%
              distinct(SHEET_NO,year)%>%
              group_by(year)%>%
              tally()%>%
              filter(n>50)
n.shots=enough.data.yrs
enough.data.yrs=enough.data.yrs%>%pull(year)
size.titl=14
size.subtitl=12
a=Dat_obs%>%
              left_join(Comm.disc.sp,by="SPECIES")%>%
              filter(!SPECIES%in%c('XX','XX.T'))%>%
              mutate(FATE=ifelse(FATE=="C","Retained",
                          ifelse(FATE=="D","Discarded",NA)),
                     FATE=factor(FATE))%>%
              filter(year%in%enough.data.yrs)%>%
              group_by(year,FATE)%>%
              summarise(count=sum(Number))%>% 
                mutate(perc = count/sum(count))%>%
              left_join(n.shots,by="year")
a$n[duplicated(a$n)] <- ""
Infographic1=a%>%
        ggplot(aes(x = year, y = perc*100, fill = FATE)) +
          geom_bar(stat="identity", width = .9) +
          labs(x = "Year",y="Percentage", y = NULL) +
        geom_text(aes(x = year,y = 100, label = n),
                  position = "identity",vjust = -0.5, 
                  fontface='bold',color="grey30",size = 2.6)+
          theme_minimal(base_size = 18) +
        theme(legend.title=element_blank(),
              legend.position="top",
              legend.direction="horizontal",
              legend.justification="right",
              legend.margin=margin(0,5,0,0),
              legend.box.margin=margin(-50,-10,-60,-10),
              panel.grid.major.y = element_line( size=1, color="grey60" ),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              legend.text = element_text(size = 12),
              axis.text.x = element_text(margin = margin(t = -10)),
              plot.title = element_text(size =size.titl,face="bold"),
              plot.subtitle = element_text(size =size.subtitl,face = "italic"),
              plot.caption = element_text(size =8,face = "italic"),
              plot.title.position = "plot")+
              labs(title ="Observed percentage of total catch (by number)",
                   subtitle = "(the annual number of observed shots is shown above each bar)",
                   caption="(Only years with more than 50 observed shots are displayed)")
  
# remove stuff-------------------------------------------
rm(DATA,Dat_obs,Dat_obs.LL,Dat_total.LL,DATA.bio,DATA.ecosystems)




# Get discard ratio by block and species (D_i_h)-------------------------------------------
STRATA_obs=function(d,Strata,Min.obs.per.block,Min.shots.per.block)
{
  #select minimum number of observations for use in analysis
  tt <- table(d$BLOCK)
  d <- d[d$BLOCK %in% names(tt[tt >= Min.obs.per.block]), ]
  
  tt <- with(d[!duplicated(d$SHEET_NO),],table(BLOCK))
  d <- d[d$BLOCK %in% names(tt[tt >= Min.shots.per.block]), ]
  
  #get discard : retain ratios
  a=d %>%
    group_by(.dots=Strata) %>%
    summarise(total = sum(Catch))%>%
    spread(key = !! parse_expr(Strata[1]), value = total)%>%
    as.data.frame
  a[is.na(a)]=0
  Vars=names(a)[-match(c(Strata[-1],"Retained"),names(a))]
  a[,Vars]=a[,Vars]/a$Retained
  names(a)[match(Vars,names(a))]=paste(names(a[,Vars]),".ratio",sep="")

  #get bycatch proportional species composition 
  prop=subset(d,Discarded=="Discarded") %>%
    group_by(.dots=c("SPECIES.ori",Strata[-1])) %>%
    summarise(total = sum(Catch))%>%
    spread(key = SPECIES.ori, value = total)%>%
    as.data.frame
  prop[is.na(prop)]=0
  Vars.prop=names(prop)[-match(c(Strata[-1]),names(prop))]
  prop[,Vars.prop]=prop[,Vars.prop]/rowSums(prop[,Vars.prop])
  
  #keep only proportion for grouped species
  id.prop=names(prop)[which(!names(prop)%in%Vars)]
  if(length(id.prop)>1)
  {
    prop=prop[,id.prop]
    
    #rescale to sum 'Grouped' to 1
    prop[,-1]=prop[,-1]/rowSums(prop[,-1])
    prop[is.na(rowSums(prop[,-1])),-1]=0
  }

  
  return(list(dat=a,prop=prop,vars=paste(Vars,".ratio",sep="")))
}
  #Elasmos
Obs_ratio.strata=STRATA_obs(d=DATA_obs$GN,
                            Strata=c("Discarded_sp",STRTA.obs.disc),
                            Min.obs.per.block=Min.obs.per.block,
                            Min.shots.per.block=Min.shots.per.block)
write.csv(Obs_ratio.strata$dat,"Results/Table1_observed.ratios.csv",row.names = F)

  #Teleosts
Obs_ratio.strata_teleosts=STRATA_obs(d=DATA_obs_teleosts$GN,
                            Strata=c("Discarded_sp",STRTA.obs.disc),
                            Min.obs.per.block=Min.obs.per.block,
                            Min.shots.per.block=Min.shots.per.block)
write.csv(Obs_ratio.strata_teleosts$dat,"Results/Recons.scalefish/Table1_observed.ratios.csv",row.names = F)


# Get stratified total reported retained catch (Y_h)-------------------------------------------
STRATA_total=function(d.gn,d.ll,Strata)
{
  if(!"YEAR"%in%Strata) Strata=c("YEAR",Strata)
  d=rbind(d.gn%>%dplyr::select(YEAR,BLOCK,Catch),d.ll%>%dplyr::select(YEAR,BLOCK,Catch))
  
  a=d %>%
    group_by(.dots=Strata) %>%
    summarise(total.retained = sum(Catch))%>%
    as.data.frame
  return(a)
}
Total_strata=STRATA_total(d.gn=DATA_total$GN,d.ll=DATA_total$LL,Strata=STRTA.reported.ret)


# Add post capture mortality-------------------------------------------
  
#Elasmos 
  #from commercial catch reconstruction
PCM=rbind(PCM.north,PCM.south)%>%
  distinct(Name,.keep_all = T)%>%
  dplyr::select(-c(Trawl,LL))%>%
  rename(PCM=GN)%>%
  left_join(All.species.names%>%dplyr::select(-c(Scien.nm,Family,SPECIES)),by='Name')%>%
  rename(SPECIES=SP)%>%
  mutate(PCM=ifelse(SPECIES=='AA',0.405,PCM))%>%
  filter(!Name=='Cobbler wobbegong')

  #collated by Agustin from literature
this.PCM=Tabl.obs.GN%>%
  filter(!dummy=='Retained')%>%
  pull(dummy)
this.PCM=DATA_obs$GN%>%
  filter(Name%in%this.PCM)%>%
  distinct(Name,.keep_all = T)%>%
  pull(SPECIES)
add.PCM=Comm.disc.sp%>%
  filter(!is.na(PCS) & NATURE%in%c("R","S") & SPECIES%in%this.PCM)%>%
  mutate(PCM=1-PCS)%>%
  dplyr::select(SPECIES,COMMON_NAME,SCIENTIFIC_NAME,PCM)%>%
  rename(Name=COMMON_NAME,
         Scien.nm=SCIENTIFIC_NAME)%>%
  dplyr::select(names(PCM))

PCM=rbind(PCM,add.PCM)%>%
  filter(SPECIES%in%Discarded.SP)%>%
  filter(!is.na(PCM))%>%
  filter(!PCM==0)%>%
  filter(!SPECIES%in%c('SE','SH'))


  #add PCM for remaining discarded species
    #Angel Shark, from Braccini et al 2012; 
    #Brown-banded catshark & catsharks, set at average of Rusty and Varied carpetshark from Braccini et al 2012
    #Crocodile shark, no information available, set = to mako shark from Braccini et al 2012 as precaution
    # Southern fiddler ray no available information, set at 10% due to hardiness
    #Sliteye shark, set at 1, always dead when brought aboard
    #Guitarfish & shovelnose ray, set equal to Whitespot shovelnose 
    #Spotted shovelnose & Western shovelnose ray, set equal to Whitespot shovelnose
    #White shark, set = to mako shark from Braccini et al 2012 
    #Whitespot shovelnose, Braodhurst & Cullis 2020 (immediate mortality)+25%
    #Cobbler wobbegong, set equal to Orectolobus parvimaculatus
    #Dwarf spotted, western and floral wobbegongs, set equal to Cobbler Wobbegong from Comm.disc.sp    
    #Spikey dogfish, set equal to Squalus spp.
    # Stingrays, set equal to Myliobatidae
PCM.remaining=data.frame(
  Name=c("Eagle ray","Brown-banded catshark","Crocodile shark",            
         "Southern fiddler ray","Sliteye shark","Guitarfish & shovelnose ray",
         "Spotted shovelnose","Dwarf spotted wobbegong","White shark",               
         "Whitespot shovelnose","Other shark",
         'Cobbler wobbegong','Western Wobbegong','Floral banded wobbegong',
         'Spikey dogfish','Stingrays',
         'Catsharks','Western shovelnose ray'),
  Scien.nm=c("Myliobatis australis","Chiloscyllium punctatum","Pseudocarcharias kamoharai",
             "Trygonorrhina dumerilii","Loxodon macrorhinus","Families Rhinobatidae & Rhynchobatidae",
             "Aptychotrema timorensis","Orectolobus parvimaculatus","Carcharodon Carcharias",
             "Rhynchobatus australiae","",
             "Sutorectus tentaculatus",'Orectolobus hutchinsi','Orectolobus floridus',
             "Squalus megalops","Dasyatidae",
             'Scyliorhinidae','Aptychotrema vincentiana'),
  PCM=c(1-0.854,1-mean(c(0.771,0.687)),1-0.242,
        1-.9,1,(0.29+0.29*.25),
        (0.29+0.29*.25),0.057,1-0.242,
        (0.29+0.29*.25),NA,
        0.0570,0.0570,0.0570,
        0.1729,0.1460,
        0.2710,0.3625),
  SPECIES=c("ER","BC","CR","FR","SE","SH","SS","WM","WP","WR","XX",
            "PD",'WW',"WF", 
            "SR","WC",
            "CA","AV"))



PCM=rbind(PCM,PCM.remaining)%>%
  filter(!Name=='Spur Dog')%>%
  mutate(PCM=ifelse(SPECIES=="PJ",2*PCM,PCM))   # bump up PCM for Port Jackson sharks

#Teleosts
PCM_teleosts=DATA_obs_teleosts$GN%>%
                filter(SPECIES.ori%in%Discarded.SP_teleosts)%>%
                distinct(SPECIES.ori,Name)%>%
                mutate(SPECIES=substr(SPECIES.ori,1,2),
                       PCM=1)


# Annual discards by species (sum(D_i_h x P_i x Y_h)) -------------------------------------------
set.seed(666) 
fn.total=function(disc.dat,tot.dat,pcm,Impute,Dist)
{

  #combine observed ratios and total reported catch
  Total.discard=merge(disc.dat$dat,tot.dat,by=STRTA.reported.ret,all.y=T)
  
  #linear interpolate missing block discard ratio
  Vars=disc.dat$vars
  for(v in 1:length(Vars)) 
  {
    id=match(Vars[v],names(Total.discard))
    
    if(Impute=='linear')
    {
      Total.discard[,id]=ifelse(is.na(Total.discard[,id]),
                                na.approx(zoo(Total.discard[,id])),
                                Total.discard[,id]) 
    }
    if(Impute=='random')
    {
      x=Total.discard[,id]
      nn=length(x[is.na(x)])
      x[is.na(x)] <- sample(x[!is.na(x)], nn, replace = TRUE)
      Total.discard[,id]=x
    }
     
  }
  
  #multiply ratio by total
  Total.discard[,Vars]=Total.discard[,Vars]*Total.discard$total.retained
  
  #rename
  i.dvar=match(Vars,names(Total.discard))
  names(Total.discard)[i.dvar]= paste("total",
      sapply(strsplit(names(Total.discard)[i.dvar], ".", fixed = TRUE), "[", 1),sep=".")
  
  #Split up 'Other'
  if(!is.na(match("total.Grouped",names(Total.discard))))
  {
    dummy=Total.discard[,c("BLOCK","YEAR","total.Grouped")]
    dummy=merge(disc.dat$prop,dummy,by="BLOCK",all.y=T)
    
    #linear interpolate missing block discard ratio
    Vars=subset(names(disc.dat$prop),!names(disc.dat$prop)=="BLOCK")
    for(v in 1:length(Vars)) 
    {
      id=match(Vars[v],names(dummy))
      
      if(Impute=='linear')
      {
        dummy[,id]=ifelse(is.na(dummy[,id]),na.approx(zoo(dummy[,id])),dummy[,id])
      }
      if(Impute=='random')
      {
        x=dummy[,id]
        nn=length(x[is.na(x)])
        x[is.na(x)] <- sample(x[!is.na(x)], nn, replace = TRUE)
        dummy[,id]=x
      }
       
    }
    
    #multiply ratio by total
    dummy[,Vars]=dummy[,Vars]*dummy$total.Grouped
    
    #rename
    i.dvar=match(Vars,names(dummy))
    names(dummy)[i.dvar]= paste("total",
        sapply(strsplit(names(dummy)[i.dvar], ".", fixed = TRUE), "[", 1),sep=".")
    
    #merge to total.discard
    Total.discard=merge(Total.discard,dummy,by=c("BLOCK","YEAR","total.Grouped"))
    Total.discard=Total.discard[,-match("total.Grouped",names(Total.discard))]
  }
  Total.discard=Total.discard[,-match("Retained",names(Total.discard))]
  
  #Add Post capture mortality
  id=colnames(Total.discard)
  id=id[-match(c("BLOCK","YEAR","total.retained"),id)]

  for(s in 1:length(id))
  {
    this=match(id[s],colnames(Total.discard))
    this.pcm=pcm%>%filter(SPECIES==str_remove(id[s], "total."))%>%pull(PCM)
    Total.discard[,this]=Total.discard[,this]*this.pcm[1]
  }
  
  #Remove blocks outside species' distribution
  Dist=Dist%>%
    mutate(SPECIES=paste("total",SPECIES,sep='.'))
  this.colmns=match(Dist$SPECIES,colnames(Total.discard))
  this.colmns=subset(this.colmns,!is.na(this.colmns))
  for(q in this.colmns)
  {
    dd=Dist%>%filter(SPECIES==colnames(Total.discard)[q])
    if(!(dd$Nor.east_block==2612 & dd$South.west_block==3330))
    {
      blks.outside.dist=which(!Total.discard$BLOCK%in%dd$Nor.east_block:(dd$South.west_block+300))
      Total.discard[blks.outside.dist,q]=0
    }
  }
  
  return(Total.discard)
}

  #Elasmos
Tot.discard.result=fn.total(disc.dat=Obs_ratio.strata,
                            tot.dat=Total_strata,
                            pcm=PCM,
                            Impute='linear',
                            Dist=elasmos_dist)

  #Teleosts
Tot.discard.result_teleosts=fn.total(disc.dat=Obs_ratio.strata_teleosts,
                                     tot.dat=Total_strata,
                                     pcm=PCM_teleosts,
                                     Impute='linear',
                                     Dist=teleosts_dist)


# Sensitivity tests ---------------------------------------------------------
if(do.sensitivity)
{
  #S1: PCM set at 100%
  PCM.S1=PCM%>%
    mutate(PCM=1)
  S1=fn.total(disc.dat=Obs_ratio.strata,
              tot.dat=Total_strata,
              pcm=PCM.S1,
              Impute='linear',
              Dist=elasmos_dist)
  
  #S2: Min obs per block = 20 & min shots per block = 10 
  Obs_ratio.strata.S2=STRATA_obs(d=DATA_obs$GN,
                                 Strata=c("Discarded_sp",STRTA.obs.disc),
                                 Min.obs.per.block=20,
                                 Min.shots.per.block=10)
  S2=fn.total(disc.dat=Obs_ratio.strata.S2,
              tot.dat=Total_strata,
              pcm=PCM,
              Impute='linear',
              Dist=elasmos_dist)
  
  
  #S3: Min obs per block = 5 & min shots per block = 2
  Obs_ratio.strata.S3=STRATA_obs(d=DATA_obs$GN,
                                 Strata=c("Discarded_sp",STRTA.obs.disc),
                                 Min.obs.per.block=5,
                                 Min.shots.per.block=2)
  S3=fn.total(disc.dat=Obs_ratio.strata.S3,
              tot.dat=Total_strata,
              pcm=PCM,
              Impute='linear',
              Dist=elasmos_dist)
  
  #S4: unobserved blocks imputed from random sample of observed blocks
  S4=fn.total(disc.dat=Obs_ratio.strata,
              tot.dat=Total_strata,
              pcm=PCM,
              Impute='random',
              Dist=elasmos_dist)
  
  
  #Display scenarios
  colfunc <- colorRampPalette(c("brown3","darkorange"))
  col.scen=colfunc(4)
  this.sp.sen=colnames(Tot.discard.result)
  this.sp.sen=this.sp.sen[-match(c("BLOCK","YEAR","total.retained"),this.sp.sen)]
  this.sp.sen=sort(this.sp.sen)
  
  names(this.sp.sen)=All.species.names%>%
                        filter(SP%in%sapply(strsplit(this.sp.sen, "total."),"[",2))%>%
                        arrange(SP)%>%
                        pull(Name)
  this.sp.sen=this.sp.sen[match(sort(names(this.sp.sen)),names(this.sp.sen))]
  
  fn.plt.sens=function(BC,s1,s2,s3,s4)
  {
    names(BC)[3]='catch'
    BC=BC%>%group_by(YEAR)%>%summarise(catch=sum(catch/1000))
    
    names(s1)[3]='catch'
    s1=s1%>%group_by(YEAR)%>%summarise(catch=sum(catch/1000))
    
    names(s2)[3]='catch'
    s2=s2%>%group_by(YEAR)%>%summarise(catch=sum(catch/1000))
    
    names(s3)[3]='catch'
    s3=s3%>%group_by(YEAR)%>%summarise(catch=sum(catch/1000))
    
    names(s4)[3]='catch'
    s4=s4%>%group_by(YEAR)%>%summarise(catch=sum(catch/1000))
    
    Ymax=max(c(max(BC$catch),max(s1$catch),max(s2$catch),max(s3$catch),max(s4$catch)))
    Ylim=c(min(log(1e-1),log(BC$catch)),log(Ymax))
    plot(BC$YEAR,log(BC$catch),ylim=Ylim,ylab='',xlab='',type='l',lwd=2,yaxt='n')
    lines(s1$YEAR,log(s1$catch),lwd=2,col=col.scen[1],lty=3)
    lines(s2$YEAR,log(s2$catch),lwd=2,col=col.scen[2],lty=3)
    lines(s3$YEAR,log(s3$catch),lwd=2,col=col.scen[3],lty=3)
    lines(s4$YEAR,log(s4$catch),lwd=2,col=col.scen[4],lty=3)
    # if(min(s1$catch)>max(BC$catch))
    # {
    #   AX=c(round(quantile(BC$catch,probs=c(0.1,.5)),1),
    #        10*round(seq(min(s1$catch),round(max(s1$catch)),length.out = 2)/10))
    #   LBL=round(AX)
    # }else
    # {
    #   AX=seq(Ymin,Ymax,length.out=4)
    #   LBL=round(AX,2)
    #   LBL=ifelse(LBL<0.01,'<0.01',LBL)
    # }
    AX=seq(Ylim[1],Ylim[2],length.out =4)
    LBL=exp(AX)
    LBL=ifelse(trunc(LBL)>=10,round(LBL),
        ifelse(trunc(LBL)<10 & trunc(LBL)>1,round(LBL,1),
        round(LBL,2)))
    LBL=ifelse(LBL<0.01,'<0.01',LBL)
    axis(2,at=AX,labels = LBL)
    #axis(2,at=log(AX),labels = LBL)
    

  }
  if(do.paper)
  {
    tiff(file="Results/Figure3_Sensitivity.tiff",width = 2400, height = 2200,
         units = "px", res = 300, compression = "lzw")    
    smart.par(n.plots=length(this.sp.sen),MAR=c(1,2,2,1),OMA=c(2,2,.1,.12),MGP=c(1,.6,0))
    par(las=1)
    for(i in 1:length(this.sp.sen))
    {
      fn.plt.sens(BC=Tot.discard.result[,c('BLOCK','YEAR',this.sp.sen[i])],
                  s1=S1[,c('BLOCK','YEAR',this.sp.sen[i])],
                  s2=S2[,c('BLOCK','YEAR',this.sp.sen[i])],
                  s3=S3[,c('BLOCK','YEAR',this.sp.sen[i])],
                  s4=S4[,c('BLOCK','YEAR',this.sp.sen[i])])
      lgn=names(this.sp.sen[i])
      mtext(lgn,3,cex=.88) 
    }
    mtext("Finacial year",1,1,outer=T,cex=1.15)
    mtext("Total discard (tonnes)",2,0.7,outer=T,las=3,cex=1.15)
    #plot.new()
    legend("bottom",c("Base case","S1","S2","S3","S4"),lty=c(1,rep(3,4)),lwd=2.5,
           bty='n',col=c("#000000",col.scen),cex=1)
    dev.off()
  }

}

# Show blocks interpolated by species ---------------------------------------------------------
fished.yrs.blocks=DATA_total$GN%>%
              dplyr::select(BLOCK,YEAR)%>%
              mutate(yr_blk=paste(YEAR,BLOCK))%>%
              distinct(yr_blk)%>%pull(yr_blk)
observed.yrs.blocks=DATA_obs$GN%>%
              dplyr::select(BLOCK,year)%>%
              mutate(yr_blk=paste(year,BLOCK))%>%
              distinct(yr_blk)%>%pull(yr_blk)
interpolated.yrs.blocks=fished.yrs.blocks[which(!fished.yrs.blocks%in%observed.yrs.blocks)]
write.csv(interpolated.yrs.blocks,"Results/Table_interpolated.year_blocks.csv",row.names=F)
write.csv(fished.yrs.blocks,"Results/Table_total.year_blocks.csv",row.names=F)

fished.blocks=unique(DATA_total$GN$BLOCK)
observed.blocks=unique(DATA_obs$GN$BLOCK)
interpolated.blocks=fished.blocks[which(!fished.blocks%in%observed.blocks)]
write.csv(interpolated.blocks,"Results/Table_interpolated.blocks.csv",row.names=F)
write.csv(fished.blocks,"Results/Table_total.blocks.csv",row.names=F)



# Uncertainty thru non-parametric bootstrap  -------------------------------------------
#note: Takes 0.4 secs per iteration
cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)
clusterCall(cl, function() {
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(zoo)
  library(stringr)
  library(tcltk)
})

  #Elasmos
Results=Tot.discard.result
system.time({
  #minimum data for analysis
  tt <- table(DATA_obs$GN$BLOCK)
  d <- DATA_obs$GN[DATA_obs$GN$BLOCK %in% names(tt[tt >= Min.obs.per.block]), ]
  
  tt <- with(d[!duplicated(d$SHEET_NO),],table(BLOCK))
  d <- d[d$BLOCK %in% names(tt[tt >= Min.shots.per.block]), ]
  
  Shots=unique(d$SHEET_NO)
  n.shots=length(Shots)
  
  store=foreach(n=1:n.boot) %dopar%
  {
    #resample observer data
    id=sample(Shots, n.shots, replace=TRUE)
    Tab=table(id)
    Tab=data.frame(SHEET_NO=names(Tab),Rep=as.numeric(Tab))
    Dat_obs.boot=d[which(d$SHEET_NO%in%id),]
    Dat_obs.boot=merge(Dat_obs.boot,Tab,by="SHEET_NO")
    Dat_obs.boot <- Dat_obs.boot[rep(row.names(Dat_obs.boot), Dat_obs.boot$Rep), ]
    
    #bootstrapped stratified observed R_h
    Obs_ratio.strata.boot=STRATA_obs(d=Dat_obs.boot,
                                     Strata=c("Discarded_sp",STRTA.obs.disc),
                                     Min.obs.per.block=Min.obs.per.block,
                                     Min.shots.per.block=Min.shots.per.block)
    
    #total discards (sum(R_h x Y_h))   
    Tot.discard.boot=fn.total(disc.dat=Obs_ratio.strata.boot,
                              tot.dat=Total_strata,
                              pcm=PCM,
                              Impute='linear',
                              Dist=elasmos_dist)
    rm(Dat_obs.boot)
    return(Tot.discard.boot)
  }
  Results=store
  rm(d)

})
rm(store)


  #Teleosts
Results_teleosts=Tot.discard.result_teleosts
system.time({
  #minimum data for analysis
  tt <- table(DATA_obs_teleosts$GN$BLOCK)
  d <- DATA_obs_teleosts$GN[DATA_obs_teleosts$GN$BLOCK %in% names(tt[tt >= Min.obs.per.block]), ]
  
  tt <- with(d[!duplicated(d$SHEET_NO),],table(BLOCK))
  d <- d[d$BLOCK %in% names(tt[tt >= Min.shots.per.block]), ]
  
  Shots=unique(d$SHEET_NO)
  n.shots=length(Shots)
  
  store=foreach(n=1:n.boot) %dopar%
    {
      #resample observer data
      id=sample(Shots, n.shots, replace=TRUE)
      Tab=table(id)
      Tab=data.frame(SHEET_NO=names(Tab),Rep=as.numeric(Tab))
      Dat_obs.boot=d[which(d$SHEET_NO%in%id),]
      Dat_obs.boot=merge(Dat_obs.boot,Tab,by="SHEET_NO")
      Dat_obs.boot <- Dat_obs.boot[rep(row.names(Dat_obs.boot), Dat_obs.boot$Rep), ]
      
      #bootstrapped stratified observed R_h 
      Obs_ratio.strata.boot=STRATA_obs(d=Dat_obs.boot,
                                       Strata=c("Discarded_sp",STRTA.obs.disc),
                                       Min.obs.per.block=Min.obs.per.block,
                                       Min.shots.per.block=Min.shots.per.block)
      
      #total discards (sum(R_h x Y_h))    
      Tot.discard.boot=fn.total(disc.dat=Obs_ratio.strata.boot,
                                tot.dat=Total_strata,
                                pcm=PCM_teleosts,
                                Impute='linear',
                                Dist=teleosts_dist)
      rm(Dat_obs.boot)
      return(Tot.discard.boot)
    }
  Results_teleosts=store
  rm(d)
  
})
rm(store)


stopCluster(cl) 





# Report ---------------------------------------------------------
tabulate.obs.eff=FALSE
if(tabulate.obs.eff)
{
  Effort_total=read.csv(handl_OneDrive('Analyses\\Data_outs\\Effort.monthly.csv'),stringsAsFactors = F)
  Effort_total.daily=read.csv(handl_OneDrive('Analyses\\Data_outs\\Effort.daily.csv'),stringsAsFactors = F)
}


#Figure 1. Map
do.map=FALSE
if(do.map)
{
  library(yarrr)
  Total.effort.year.block=Effort_total%>%
    filter(METHOD=='GN' & Shark.fishery%in%c('JASDGDL','WCDGDL') & YEAR.c>=1993)%>%
    mutate(LAT=-round(abs(LAT)),
           LAT=LAT-.5,
           LONG=LONG+.5)%>%
    group_by(Same.return,VESSEL,zone,YEAR.c,LAT,LONG)%>%
    summarise(effort=max(Km.Gillnet.Hours.c,na.rm=T))%>%
    group_by(YEAR.c,LAT,LONG)%>%
    summarise(effort=sum(effort))%>%
    rename(year=YEAR.c)%>%
    filter(year<2006)
  
  Total.effort.year.block_daily=Effort_total.daily%>%
    filter(method=='GN' & Shark.fishery%in%c('JASDGDL','WCDGDL') & !is.na(Km.Gillnet.Hours.c))%>%
    mutate(LAT=-round(abs(LAT)),
           LONG=round(LONG),
           LAT=LAT-.5,
           LONG=LONG+.5)%>%
    group_by(date,vessel,zone,year.c,LAT,LONG)%>%
    summarise(effort=max(Km.Gillnet.Hours.c,na.rm=T))%>%
    group_by(year.c,LAT,LONG)%>%
    summarise(effort=sum(effort,na.rm=T))%>%
    rename(year=year.c)
  
  Spatial.effort=rbind(Total.effort.year.block,Total.effort.year.block_daily)%>%
    group_by(LAT,LONG)%>%
    summarise(effort=sum(effort,na.rm=T))%>%
    filter(!is.na(LAT)| !is.na(LONG))%>%
    spread(LAT,effort)
  X=Spatial.effort$LONG
  Y=as.numeric(colnames(Spatial.effort)[-1])
  Z=as.matrix(Spatial.effort[,-1])
  
  
  library(PBSmapping)
  data(worldLLhigh)
  library(rgdal)
  source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R"))
  
  SDGDLL_zone1=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/SDGDLL_zone1.shp", layer="SDGDLL_zone1")) 
  SDGDLL_zone2=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/SDGDLL_zone2.shp", layer="SDGDLL_zone2")) 
  WCDGDLL=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/WCDGDLL.shp", layer="WCDGDLL")) 
  Bathymetry_120=read.table(handl_OneDrive("Data/Mapping/get_data112_120.cgi"))
  Bathymetry_138=read.table(handl_OneDrive("Data/Mapping/get_data120.05_138.cgi"))
  
  Bathymetry=rbind(Bathymetry_120,Bathymetry_138)%>%
    arrange(V1,V2)%>%
    filter(V1 >= South.WA.long[1] & V1 <=South.WA.long[2])%>%
    filter(V2 >= South.WA.lat[1] & V2 <=South.WA.lat[2])
  xbat=unique(Bathymetry$V1)
  ybat=unique(Bathymetry$V2)
  reshaped=as.matrix(Bathymetry%>%spread(V2,V3))
  
  
  South.WA.long=c(109,129)
  South.WA.lat=c(-37,-26)
  Xlim=floor(c(min(X),max(X)))
  Ylim=floor(c(min(Y),max(Y)))
  
  colfunc <- colorRampPalette(c("cyan3","deepskyblue3"))
  COLn=colfunc(3)
  
  
  tiff(file="Results/Figure1.tiff",width=2400,height=2400,units="px",res=300,compression="lzw")
  par(mfcol=c(1,1),mar=c(1,1,.5,.5),oma=c(3,3,1,.3),las=1,mgp=c(.04,.65,0))
  plot(1,xlim=Xlim,ylim=Ylim,xlab="",ylab="",axes=F,main="")
  
  #effort
  image(x=X, y=Y, z=Z,add=T)
  
  axis(side = 1, seq(112,129,1),labels =F, tck = .96,col="grey40",lty=3)
  axis(side = 2, seq(South.WA.lat[1],South.WA.lat[2],1),labels =F, tck = 1,col="grey40",lty=3)
  
  #shots
  kk=DATA_obs$GN%>%distinct(SHEET_NO,.keep_all = T)%>%dplyr::select(Mid.Lat,Mid.Long)
  points(kk$Mid.Long,kk$Mid.Lat,pch=21, cex = .9,
         bg=transparent("green2",.8),col=transparent("black",.7))
  
  #axis
  axis(side = 1, seq(South.WA.long[1],South.WA.long[2],1),labels =F, tck = -.01)
  axis(side = 2, seq(South.WA.lat[1],South.WA.lat[2],1),labels = F, tck = -.01)
  axis(side = 1, seq(112,129,2),labels =seq(112,129,2), tck = -.015)
  axis(side = 2, seq(South.WA.lat[1],South.WA.lat[2],2), 
       labels = -seq(South.WA.lat[1],South.WA.lat[2],2), tck = -.015)
  
  #bathymetry
  contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=South.WA.lat,xlim=South.WA.long, zlim=c(-1,-100),
          nlevels = 1,labcex=1,labels=c("100",""),lty = 1,col="gray20",add=T)
  polygon(WAcoast$Longitude,WAcoast$Latitude, col="grey80") 
  box()
  mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1,font=1,las=0,cex=1.35,outer=T)
  mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=1,font=1,las=0,cex=1.35,outer=T)
  
  #Inset Australia
  par(fig=c(.6,1,.475,.975), new = T,mgp=c(.1,.4,0))
  plotMap(worldLLhigh, xlim=South.WA.long,ylim=c(-39,-13),plt = c(.1, 1, 0.075, 1),
          col="Black",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
  box()
  plot(WCDGDLL,add=T,col=COLn[1])
  plot(SDGDLL_zone1,add=T,col=COLn[2])
  plot(SDGDLL_zone2,add=T,col=COLn[3])
  CL.txt='midnightblue'
  text(123,-24,"Western",col='white', cex=1.5)
  text(123,-27,"Australia",col='white', cex=1.5)
  
  text(112.5,-29,"WCDGDLF",col=CL.txt,cex=.85,srt=-70)
  text(114,-35.25,"JASDGDLF",col=CL.txt,cex=.85,srt=-60)
  text(113,-36,"(Zone 1)",col=CL.txt,cex=.85,srt=-60)
  text(123.7437,-35,"JASDGDLF",col=CL.txt,cex=.85)
  text(123.7437,-36.25,"(Zone 2)",col=CL.txt,cex=.85)
  
  dev.off()
}

#Table 2. List of retained and discarded species
if(do.paper) write.csv(Table.disc.obs,"Results/Table_2.csv",row.names=F)

#Table of PCM and weights  
if(do.paper)
{
  Table_S.2=left_join(Len.wei,
                      All.species.names[,-match('SPECIES',names(All.species.names))],by=c('SPECIES'='SP'))%>%
    filter(!is.na(Name))%>%
    full_join(PCM[,-match(c('Name','Scien.nm'),names(PCM))],by='SPECIES')%>%
    dplyr::select(Name,Scien.nm,PCM,a_w8t,b_w8t,Source)%>%
    mutate(PCM=ifelse(is.na(PCM),'',PCM))
  write.csv(Table_S.2,"Results/Table_S.2.csv",row.names=F)
  
}

# Table of effort coverage 
#note: Show porportion of effort observed by year for gillnet
if(tabulate.obs.eff)
{
  Total.effort.year=Effort_total%>%
    filter(METHOD=='GN' & Shark.fishery%in%c('JASDGDL','WCDGDL') & YEAR.c>=1993)%>%
    group_by(Same.return,VESSEL,zone,YEAR.c)%>%
    summarise(effort=max(Km.Gillnet.Hours.c,na.rm=T))%>%
    group_by(YEAR.c)%>%
    summarise(effort=sum(effort))%>%
    rename(year=YEAR.c)%>%
    filter(year<2006)
  
  Total.effort.year_daily=Effort_total.daily%>%
    filter(method=='GN' & Shark.fishery%in%c('JASDGDL','WCDGDL') & !is.na(Km.Gillnet.Hours.c))%>%
    group_by(date,vessel,zone,year.c)%>%
    summarise(effort=max(Km.Gillnet.Hours.c,na.rm=T))%>%
    group_by(year.c)%>%
    summarise(effort=sum(effort,na.rm=T))%>%
    rename(year=year.c)
  
  Total.effort.year=rbind(Total.effort.year,Total.effort.year_daily)
  
  
  Obs.effort.year=DATA_obs$GN%>%
    mutate(effort=NET_LENGTH*SOAK.TIME)%>%
    distinct(SHEET_NO,.keep_all = T)%>%
    group_by(BOAT)%>%
    mutate(mean.effort=mean(effort,na.rm=T))%>%
    mutate(effort=case_when(is.na(effort)~mean.effort,
                            TRUE~effort))%>%
    group_by(year)%>%
    summarise(effort.obs=sum(effort))
  
  
  Table_S.1=left_join(Total.effort.year,Obs.effort.year,by='year')%>%
    mutate(Percent.obs=100*effort.obs/effort)
  
  tiff(file="Results/FigureS1_observed.effort.tiff",width = 2000, height = 2000,
       units = "px", res = 300, compression = "lzw") 
  par(mar=c(2.5,2.5,1,1),oma=c(1,1,.1,.1),mgp=c(1.5,.5,0),las=1,cex.lab=1.25)
  plot(Table_S.1$year,Table_S.1$Percent.obs,ylab="Percentage observed",xlab="Year",
       pch=21,bg="brown1",cex=2)
  dev.off()
}


#Put results in right format, aggregate blocks and extract median, ci    
fn.agg.block=function(d)
{
  dummy=vector('list',n.boot)
  for(n in 1:n.boot)
  {
    Var=subset(names(d[[n]]),!names(d[[n]])%in%c("BLOCK","YEAR"))
    Var.ag=paste('cbind(',paste(Var,collapse=','),")",sep='')
    Formula=as.formula(paste(Var.ag,"YEAR",sep="~"))
    dummy[[n]]=aggregate(Formula,d[[n]],sum)
  }
  return(dummy)
}

fn.add.missing=function(d,VAR)
{
  VAR=subset(VAR,!VAR%in%c("BLOCK"))
  VAR1=subset(VAR,!VAR%in%c("YEAR","total.retained"))
  for(n in 1:n.boot)
  {
    id=which(!VAR%in%names(d[[n]]))
    if(length(id>0))
    {
      iid=which(!VAR1%in%names(d[[n]]))
      dummy=as.data.frame(d[[n]][,rep(1,length(iid))])
      names(dummy)=VAR1[iid]
      dummy[,]=0
      d[[n]]=cbind(d[[n]],dummy)
    }
    d[[n]]=d[[n]][,match(VAR,names(d[[n]]))]
  }
  return(d)
}

#aggregate blocks and add missing species     #0.018 sec per iteration
system.time({ 
  Results.show=fn.add.missing(d=fn.agg.block(d=Results),VAR=names(Tot.discard.result))
  Results.show_teleosts=fn.add.missing(d=fn.agg.block(d=Results_teleosts),
                                       VAR=names(Tot.discard.result_teleosts))
})

# Stats
PRBS=c(.025,.5,.975)
Stat.list=list(Species="Species",Total="Total")
Stats=Stats_teleosts=Stat.list
fn.stats=function(d,how)
{
  if(how=="Species")
  {
    all.matrix <- abind(d, along=3)
    Store=vector('list',ncol(all.matrix))
    names(Store)=colnames(all.matrix)
    for(nm in 2:ncol(all.matrix)) Store[[nm]]=apply(all.matrix[,nm,],1,function(x) quantile(x/1000,probs=PRBS)) #in tonnes
  }
  if(how=="Total")
  {
    id.var=which(!colnames(d[[1]])%in%c("total.retained","total.retained","YEAR"))
    for(n in 1:n.boot) d[[n]]$total.discarded=rowSums(d[[n]][,id.var])
    all.matrix <- abind(d, along=3)
    These=c("YEAR","total.discarded","total.retained")
    Store=vector('list',length(These))
    names(Store)=These
    id.these=match(These,colnames(d[[1]]))
    for(nm in 2:length(id.these)) Store[[nm]]=apply(all.matrix[,id.these[nm],],1,function(x) quantile(x/1000,probs=PRBS))
  }
  Store$YEAR=all.matrix[,,1][,1]
  
  return(Store)
}
  #Elasmos
for(l in 1:length(Stats)) Stats[[l]]=fn.stats(d=Results.show,how=Stat.list[[l]])
Plt.this=Stats$Species[-match(c("YEAR","total.retained"),names(Stats$Species))]
Plt.this=Plt.this[sort(names(Plt.this))]

  #Teleosts
for(l in 1:length(Stats_teleosts)) Stats_teleosts[[l]]=fn.stats(d=Results.show_teleosts,
                                                                how=Stat.list[[l]])
Plt.this_teleosts=Stats_teleosts$Species[-match(c("YEAR","total.retained"),names(Stats_teleosts$Species))]
Plt.this_teleosts=Plt.this_teleosts[sort(names(Plt.this_teleosts))]


# Plot statistics
yrs=sort(unique(DATA_total$GN$YEAR))
Col.totl="black"
Col.totl.retained="grey60"

  #plotting function
id.med=match("50%",rownames(Stats$Species[[2]]))
id.low=match("2.5%",rownames(Stats$Species[[2]]))
id.up=match("97.5%",rownames(Stats$Species[[2]]))

  #Elasmos
    #By species  
if(do.paper)
{
  fn.plt=function(d,CL,LWd,what)
  {
    plot(yrs,yrs,col='transparent',ylim=c(0,YMAX),ann=F)
    if(what=='polygon')
    {
      polygon(x=c(yrs,rev(yrs)),y=c(d[id.low,],rev(d[id.up,])),
              col=adjustcolor(CL, alpha.f = 0.30), border = adjustcolor(CL, alpha.f = 0.60))
      lines(yrs,d[id.med,],lwd=LWd,col=CL)
    }
    if(what=='points')
    {
      segments(yrs,d[id.low,],yrs,d[id.up,],col="grey50")
      points(yrs,d[id.med,],pch=19,col=CL,cex=.8)
      lines(yrs,d[id.med,],lty=3,col=CL)
    }
    
  }
  
  Plt.this.sorted=data.frame(Plt.this=names(Plt.this))%>%
    mutate(names=sapply(strsplit(Plt.this, "total."),"[",2))%>%
    left_join(All.species.names%>%dplyr::select(Name,SP),by=c("names"="SP"))%>%
    arrange(Name)
  Plt.this=Plt.this[Plt.this.sorted$Plt.this]
  tiff(file="Results/Figure4_by.species.tiff",width = 2400, height = 2200,
       units = "px", res = 300, compression = "lzw")    
  par(las=1)
  smart.par(n.plots=length(Plt.this),MAR=c(1,2,2,1),OMA=c(2,2,.1,.12),MGP=c(1,.6,0))
  for(i in 1:length(Plt.this))
  {
    YMAX=max(unlist(Plt.this[[i]]))
    fn.plt(d=Plt.this[[i]],CL=Col.totl,LWd=1.5,what="points")
    lgn=All.species.names%>%
      filter(SP==unlist(strsplit(names(Plt.this)[i], ".", fixed = TRUE))[2])%>%
      pull(Name)
    mtext(lgn,3,cex=.88) 
  }
  mtext("Finacial year",1,1,outer=T,cex=1.15)
  mtext("Total discard (tonnes)",2,0.5,outer=T,las=3,cex=1.15)
  dev.off()
}

    #Total
if(do.paper)
{
  fn.plt.prop=function(d,d.ret,CL,LWd)
  {
    plot(yrs,yrs,col='transparent',ylim=c(0,YMAX),ann=F)
    polygon(x=c(yrs,rev(yrs)),y=c(d[id.low,]/d.ret,rev(d[id.up,]/d.ret)),
            col=adjustcolor(CL, alpha.f = 0.30), border = adjustcolor(CL, alpha.f = 0.60))
    lines(yrs,d[id.med,]/d.ret,lwd=LWd,col=CL)
  }
  YMAX=.1
  jpeg("Results/Total.discard.estimates/Total.jpg",width=2400,height=2400,units="px",res=300)
  par(mfcol=c(1,1),mar=c(1.5,1,1,1),oma=c(1,3,1,1),mgp=c(1,.6,0),las=1)
  fn.plt.prop(d=Stats$Total$total.discarded,d.ret=Stats$Total$total.retained[id.med,],CL=Col.totl,LWd=3)
  mtext("Proportion",2,1.5,las=3,outer=T,cex=1.5)
  mtext("Year",1,0,outer=T,cex=1.5)
  dev.off()
}

  #Teleosts
    #By species 
All.species.names_teleosts=DATA_obs_teleosts$GN%>%
                              filter(SPECIES.ori%in%Discarded.SP_teleosts)%>%
                              distinct(SPECIES.ori,Name)%>%
                              mutate(SPECIES=substr(SPECIES.ori,1,2),
                                     Name=capitalize(tolower(Name)),
                                     Name=ifelse(Name=='Sergeant baker','Sergeant Baker',Name))

if(do.paper)
{
  fn.plt=function(d,CL,LWd,what)
  {
    plot(yrs,yrs,col='transparent',ylim=c(0,YMAX),ann=F)
    if(what=='polygon')
    {
      polygon(x=c(yrs,rev(yrs)),y=c(d[id.low,],rev(d[id.up,])),
              col=adjustcolor(CL, alpha.f = 0.30), border = adjustcolor(CL, alpha.f = 0.60))
      lines(yrs,d[id.med,],lwd=LWd,col=CL)
    }
    if(what=='points')
    {
      segments(yrs,d[id.low,],yrs,d[id.up,],col="grey50")
      points(yrs,d[id.med,],pch=19,col=CL,cex=.8)
      lines(yrs,d[id.med,],lty=3,col=CL)
    }
    
  }
  tiff(file="Results/Recons.scalefish/Figure4_by.species.tiff",width = 2150, height = 2000,
       units = "px", res = 300, compression = "lzw")    
  par(las=1)
  smart.par(n.plots=length(Plt.this_teleosts),MAR=c(1,1.75,1.75,1.5),OMA=c(2,2,.1,.1),MGP=c(1,.5,0))
  for(i in 1:length(Plt.this_teleosts))
  {
    YMAX=max(unlist(Plt.this_teleosts[[i]]))
    fn.plt(d=Plt.this_teleosts[[i]],CL=Col.totl,LWd=1.5,what="points")
    lgn=All.species.names_teleosts%>%
      filter(SPECIES==unlist(strsplit(names(Plt.this_teleosts)[i], ".", fixed = TRUE))[2])%>%
      pull(Name)
    mtext(lgn,3,cex=.8) 
  }
  mtext("Year",1,1,outer=T,cex=1.15)
  mtext("Total discard (tonnes)",2,0.5,outer=T,las=3,cex=1.15)
  dev.off()
}

# Plot Size frequency distributions
if(do.paper)
{
  #Elasmos
  DATA_obs$GN%>%
    filter(SPECIES.ori%in%unlist(lapply(strsplit(names(Plt.this), '.', fixed = TRUE), '[', 2)))%>%
    filter(!Name=='Stingrays')%>%  #Stringrays set to average size as not measured, no point displaying
    mutate(Name=ifelse(Name=="Eagle ray","Southern eagle ray",Name))%>%
    ggplot( aes(x=FL, color=Name, fill=Name)) +
    geom_histogram(alpha=0.6, binwidth = 5) +
    theme(legend.position="none",
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size =11.5),
          axis.title=element_text(size=16)) +
    xlab("Size (cm)") +
    ylab("Frequency") +
    facet_wrap(~Name, scales = "free")
  ggsave('Results/FigureS2_size.frequency.tiff',width = 10,height = 10,compression = "lzw")
  
  #Teleosts  
  DATA_obs_teleosts$GN%>%
    filter(SPECIES.ori%in%
             paste(unlist(lapply(strsplit(names(Plt.this_teleosts), '.', fixed = TRUE), '[', 2)),'.T',sep=''))%>%
    mutate(Name=capitalize(tolower(Name)),
           Name=ifelse(Name=="Sergeant baker","Sergeant Baker",Name))%>%
    ggplot( aes(x=TL, color=Name, fill=Name)) +
    geom_histogram(alpha=0.6, binwidth = 5) +
    theme(legend.position="none",
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size =11.5),
          axis.title=element_text(size=16)) +
    xlab("Total length (cm)") +
    ylab("Frequency") +
    facet_wrap(~Name, scales = "free")
  ggsave('Results/Recons.scalefish/FigureS2_size.frequency.tiff',width = 12,height = 10,compression = "lzw")
  
  
}

# Plot densities
if(do.paper)
{
  #Elasmos
  DATA_obs$GN%>%
    filter(SPECIES.ori%in%unlist(lapply(strsplit(names(Plt.this), '.', fixed = TRUE), '[', 2)))%>%
    filter(!Name=='Stingrays')%>%  #Stringrays set to average size as not measured, no point displaying
    mutate(Name=ifelse(Name=="Eagle ray","Southern eagle ray",Name))%>%
    ggplot(aes(x=FL,y=Name,fill=Name)) +
    ggridges::geom_density_ridges(scale = 0.95) +
    theme(legend.position="none",
          axis.text=element_text(size=14),
          axis.title=element_text(size=16)) +
    xlab("Size (cm)") +
    ylab("Density")
  ggsave('Results/Density_size.frequency.tiff',width = 10,height = 10,compression = "lzw")
  
  #Teleosts  
  DATA_obs_teleosts$GN%>%
    filter(SPECIES.ori%in%
             paste(unlist(lapply(strsplit(names(Plt.this_teleosts), '.', fixed = TRUE), '[', 2)),'.T',sep=''))%>%
    mutate(Name=capitalize(tolower(Name)),
           Name=ifelse(Name=="Sergeant baker","Sergeant Baker",Name))%>%
    ggplot(aes(x=TL,y=Name,fill=Name)) +
    ggridges::geom_density_ridges(scale = 0.95) +
    theme(legend.position="none",
          axis.text=element_text(size=14),
          axis.title=element_text(size=16)) +
    xlab("Total length (cm)") +
    ylab("Density")
  ggsave('Results/Recons.scalefish/Density_size.frequency.tiff',width = 10,height = 10,compression = "lzw")
  
}

# Export total discard estimates-----------------------------------------------------------------------
setwd(handl_OneDrive('Analyses\\Data_outs'))

  #Elasmos
    #Base Case
YR=Results.show[[1]]$YEAR
out=vector('list',length(Plt.this))
for(o in 1:length(out))
{
  dd=Plt.this[[o]]
  nm=All.species.names%>%
    filter(SP==unlist(strsplit(names(Plt.this)[o], ".", fixed = TRUE))[2])%>%
    pull(SPECIES)
  dummy=data.frame(FINYEAR=paste(YR,substr(YR+1,3,4),sep="-"),
                   SPECIES=nm,
                   LIVEWT.c=1000*dd[match('50%',rownames(dd)),],
                   zone=NA)
  out[[o]]=dummy
  rm(dd,dummy,nm)
}
Out.elasmos=do.call(rbind,out)
write.csv(Out.elasmos,"recons_discard_TDGDLF.csv",row.names = F)

  #S1
get.this=names(S1)
get.this=get.this[-match(c("BLOCK","total.retained"),get.this)]
d=S1[,get.this]%>%
  gather('SP','LIVEWT.c',-YEAR)%>%
  mutate(SP=str_extract(SP, "[^.]+$"),
         FINYEAR=paste(YEAR,substr(YEAR+1,3,4),sep="-"))%>%
  left_join(All.species.names%>%dplyr::select(SPECIES,SP),by="SP")%>%
  group_by(FINYEAR,SPECIES)%>%
  summarise(LIVEWT.c=sum(LIVEWT.c))%>%
  mutate(zone=NA)
write.csv(d,"recons_discard_TDGDLF_100.PCM.csv",row.names = F)

# for(i in 1:length(Plt.this)) 
# {
#   nm=All.species.names%>%
#     filter(SP==unlist(strsplit(names(Plt.this)[i], ".", fixed = TRUE))[2])%>%
#     pull(Name)
#   if(nm=="Australian angel shark") nm="Australian angelshark"
#   dd=as.data.frame(t(Plt.this[[i]]))
#   colnames(dd)=c("Low.CI","Median","Up.CI")
#   dd=cbind(Year=Results.show[[1]]$YEAR,dd)
#   
#   this.wd=paste(getwd(),nm,sep='/')
#   if(!dir.exists(file.path(this.wd))) dir.create(file.path(this.wd))
#  
#   write.csv(dd,paste(nm,"Total.discard.estimates_TDGDLF_tonnes.csv",sep='/'),row.names = F)
# }
  
#     #S1
# get.this=names(S1)
# get.this=get.this[-match(c("BLOCK","YEAR","total.retained"),get.this)]
# for(i in 1:length(get.this)) 
# {
#   nm=All.species.names%>%
#     filter(SP==unlist(strsplit(get.this[i], ".", fixed = TRUE))[2])%>%
#     pull(Name)
#   if(nm=="Australian angel shark") nm="Australian angelshark"
#   dd=S1[,c("YEAR",get.this[i])]
#   colnames(dd)[2]="Median"
#   dd=dd%>%
#     group_by(YEAR)%>%
#     summarise(Median=sum(Median)/1000)
#   
#   this.wd=paste(getwd(),nm,sep='/')
#   if(!dir.exists(file.path(this.wd))) dir.create(file.path(this.wd))
#   
#   write.csv(dd,paste(nm,"Total.discard.estimates_TDGDLF_100.PCM_tonnes.csv",sep='/'),row.names = F)
# }


  #Teleosts
setwd(handl_OneDrive('Analyses/Reconstruction_total_bycatch_TDGDLF/Results/Recons.scalefish'))
    #Base Case
YR=Results.show_teleosts[[1]]$YEAR
out=vector('list',length(Plt.this_teleosts))
for(o in 1:length(out))
{
  dd=Plt.this_teleosts[[o]]
  nm=All.species.names_teleosts%>%
    filter(SPECIES==unlist(strsplit(names(Plt.this_teleosts)[o], ".", fixed = TRUE))[2])%>%
    pull(SPECIES.ori)
  dummy=data.frame(FINYEAR=paste(YR,substr(YR+1,3,4),sep="-"),
                   SPECIES=nm,
                   LIVEWT.c=1000*dd[match('50%',rownames(dd)),],
                   zone=NA)
  out[[o]]=dummy
  rm(dd,dummy,nm)
}
Out.teleosts=do.call(rbind,out)
write.csv(Out.teleosts,"recons_discard_TDGDLF_teleosts.csv",row.names = F)


# Infographic-----------------------------------------------------------------------
Info3.dat=Out.elasmos%>%
        left_join(All.species.names,by='SPECIES')%>%
        dplyr::select(FINYEAR,LIVEWT.c,Name,SPECIES)%>%
        mutate(Class="elasmobranchs",
               LIVEWT.c=LIVEWT.c/1000)
Info3.dat.teleosts=Out.teleosts%>%
              left_join(All.species.names_teleosts%>%
                          mutate(SPEC=round(runif(nrow(All.species.names_teleosts),1e5,1e6))),
                        by=c('SPECIES'='SPECIES.ori'))%>%
              dplyr::select(FINYEAR,LIVEWT.c,Name,SPEC)%>%
              rename(SPECIES=SPEC)%>%
              mutate(Class="teleosts",
                     LIVEWT.c=LIVEWT.c/1000)
Info3.dat=rbind(Info3.dat,Info3.dat.teleosts)
Last.5.yrs=sort(unique(Info3.dat$FINYEAR))
Last.5.yrs=Last.5.yrs[(length(Last.5.yrs)-4):length(Last.5.yrs)]
Info3.dat=Info3.dat%>%
          group_by(Name,Class,SPECIES)%>%
          summarise(Avrg=mean(LIVEWT.c))

This.file=paste(handl_OneDrive('Analyses/Catch and effort/State of fisheries/'),
       Last.5.yrs[length(Last.5.yrs)],'/ERA_table.retained.species.csv',sep='')
if(!file.exists(This.file)) This.file=paste(handl_OneDrive('Analyses/Catch and effort/State of fisheries/'),
                                            Last.5.yrs[length(Last.5.yrs)-1],'/ERA_table.retained.species.csv',sep='')
Avrg.retained=read.csv(This.file)
                     
Avrg.dat=data.frame(Class=c("Discarded","Retained"),
                    Total=c(sum(Info3.dat$Avrg),sum(Avrg.retained$Average.catch)))                         
              
a=Avrg.dat%>% 
  mutate(perc = Total/sum(Total))%>%
  mutate(label=paste(round(100*perc/sum(perc),1),'%'),
         ymax=cumsum(perc),
         ymin = c(0, head(ymax, n=-1)),
         labelPosition =(ymax + ymin) / 2)
Infographic2=a%>%
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
  geom_rect()+
  geom_label_repel(x=3.5,aes(y = labelPosition,label = paste(round(100*perc,1),"%")),
                   size=4, color="black",fontface = 'bold',box.padding=-1.2, nudge_x = 5)+
  coord_polar(theta="y")+theme_void() + xlim(c(2, 4))+
  ggtitle(label ="Percentage of total catch (by weight)",
          subtitle = paste("(average for the last 5 years= ",round(sum(a$Total))," tonnes; discard estimates account for PCM)",sep=""))+
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.margin=unit(c(-1,-1,-1,-1), "cm"),
        plot.title = element_text(size = size.titl, face = "bold"),
        plot.subtitle = element_text(size = size.subtitl,face = "italic"),
        plot.title.position = "plot")

Info3.dat=Info3.dat%>%
      data.frame%>%
      mutate(Tot=sum(Avrg),
             perc = round(100*Avrg/Tot,2),
             Name2=ifelse(perc>2,Name,paste("Other",Class)),
             SPECIES2=ifelse(Name2=="Other elasmobranchs",6e4,
                      ifelse(Name2=="Other teleosts",1e7,
                             SPECIES)))
N.other.elas=length(Info3.dat%>%
                filter(Name2=='Other elasmobranchs')%>%
                  distinct(Name)%>%pull(Name))
N.other.teleos=length(Info3.dat%>%
                      filter(Name2=='Other teleosts')%>%
                      distinct(Name)%>%pull(Name))
Info3.dat=Info3.dat%>%
          mutate(Name2=case_when(Name2=='Other elasmobranchs'~paste(Name2," (",N.other.elas," species)",sep=""),
                                 Name2=='Other teleosts'~paste(Name2," (",N.other.teleos," species)",sep=""),
                                 TRUE~Name2))

Lvls=Info3.dat%>%
  distinct(Name2,SPECIES2)%>%
  arrange(SPECIES2)
colfunc <- colorRampPalette(c("brown3", "orange","yellow"))
MyClrs.elas=colfunc(length(Lvls%>%filter(SPECIES2<=6e4)%>%pull(SPECIES2)))
colfunc <- colorRampPalette(c("lightblue1", "dodgerblue1","cyan3"))
MyClrs.tel=colfunc(length(Lvls%>%filter(SPECIES2>6e4)%>%pull(SPECIES2)))
MyClrs=c(MyClrs.elas,MyClrs.tel)

Info3="barplot"
if(Info3=="barplot")
{
  Infographic3=Info3.dat%>%
    mutate(Name2=factor(Name2,levels=Lvls$Name2))%>%
    group_by(Name2)%>%
    summarise(perc=sum(perc))%>%
    mutate(x=1)%>%
    ggplot(aes(x=x,y = perc, fill = Name2)) +
    geom_bar(stat="identity", width = 1.2) +
    labs(x = "",y="") +
    ggtitle(label ="Percentage of discarded catch (accounting for PCM)") +
    theme_minimal(base_size = 18)  +
    theme(legend.title=element_blank(),
          legend.position = "right",
          panel.grid.major.y = element_line( size=1, color="grey60" ),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.text = element_text(size = 12),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          plot.title = element_text(size = size.titl, face = "bold"),
          plot.subtitle = element_text(size = size.subtitl,face = "italic"),
          plot.title.position = "plot") +
    scale_fill_manual(aes(name = Name2),values = MyClrs)
  
}
if(Info3=="pie")
{
  Infographic3=Info3.dat%>%
  mutate(Name2=factor(Name2,levels=Lvls$Name2))%>%
  group_by(Name2)%>%
  summarise(perc=sum(perc))%>%
  mutate(cs = rev(cumsum(rev(perc))), 
         pos = perc/2 + lead(cs, 1),
         pos = if_else(is.na(pos), perc/2, pos))%>%
  ggplot(aes(x = "" , y = perc, fill = Name2)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y", start = 0 ) +
  geom_label_repel(aes(y = pos,label = Name2), size=3.5, color="black",fontface = 'bold',
                   box.padding=.3, nudge_x = 1,
                   fill = alpha(MyClrs,0.5)) +
  guides(fill = guide_legend(title = "Status")) +
  labs(x = "",y="") +
  ggtitle(label ="Percentage of discarded catch") +
  theme_minimal(base_size = 18)  +
  theme(legend.title=element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 14),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(size = size.titl, face = "bold"),
        plot.subtitle = element_text(size = size.subtitl,face = "italic"))+
    scale_fill_manual(aes(name = Name2),values = MyClrs) 
}

#Export infographic
Infographic1 + Infographic2 / Infographic3 +plot_layout(ncol = 2, widths = unit(c(14.5, 5),c('cm','cm')))
ggsave(handl_OneDrive('Analyses/Reconstruction_total_bycatch_TDGDLF/Results/Infographic.tiff'),
       width = 13,height = 7,compression = "lzw")
