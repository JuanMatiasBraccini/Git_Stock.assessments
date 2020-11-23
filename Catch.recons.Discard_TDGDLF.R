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

User="Matias"
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R")




# 1. Data ---------------------------------------------------------

#Observers data
#keep only elasmobranchs from observed commercial boats
Res.ves=c("HAM","HOU","NAT","FLIN","RV BREAKSEA","RV Gannet","RV GANNET","RV SNIPE 2")
Dat_obs=DATA %>% filter(Mid.Lat<=(-26) & !BOAT%in%Res.ves & Taxa=='Elasmobranch' & !COMMON_NAME=='WHALE')  %>% 
                 select(c(SHEET_NO,date,Month,year,BOAT,zone,BLOCK,SOAK.TIME,MESH_SIZE,
                        MESH_DROP,NET_LENGTH,Mid.Lat,Mid.Long,Method,SPECIES,
                        COMMON_NAME,SCIENTIFIC_NAME,TL,FL,Number)) 
Dat_obs.LL=Dat_obs %>% filter(Method=="LL")
Dat_obs=Dat_obs %>% filter(Method=="GN")


#Catch and effort data   
Dat_total=read.csv('C:\\Matias\\Analyses\\Data_outs\\Data.monthly.csv',stringsAsFactors = F)


#Species names
All.species.names=read.csv("C:/Matias/Data/Species_names_shark.only.csv") #for catch


#Discarded/retained
Comm.disc.sp=read.csv("C:/Matias/Analyses/Ecosystem indices and multivariate/Shark-bycatch/SPECIES+PCS+FATE.csv",stringsAsFactors = F)


#Length weight relationships
Len.wei=read.csv("C:/Matias/Data/Length_Weights/length.weights.csv",stringsAsFactors = F)

#Weight ranges
Wei.range=read.csv("C:/Matias/Data/Length_Weights/Data.Ranges.csv")
Wei.range.names=read.csv("C:/Matias/Data/Length_Weights/Species.names.csv")

#Post capture mortality
PCM.north=read.csv('C:/Matias/Analyses/Reconstruction_catch_commercial/TableS1.PCM_North.csv')
PCM.south=read.csv('C:/Matias/Analyses/Reconstruction_catch_commercial/TableS1.PCM_South.csv')


# 2. Parameter ---------------------------------------------------------

Min.obs.per.block=10  #minimum number of observations per block for use in ratio estimator
Min.shots.per.block=5  #minimum number of shots per block for use in ratio estimator

do.explrtn="NO"
STRTA.obs.disc=c("BLOCK")       #aggregating strata for observer data
STRTA.reported.ret=c("BLOCK")   #aggregating strata for total landed retained catch

n.boot=1e3

Group_rare_criteria=0.02    #criteria for grouping rare species (proportion of catch)    

Commercial.sp=subset(Comm.disc.sp,FATE=="C" & NATURE%in%c("S","R"))$SPECIES

do.sensitivity=FALSE

# Manipulate catch ---------------------------------------------------------
setwd('C:/Matias/Analyses/Reconstruction_total_bycatch_TDGDLF')

Dat_total=Dat_total %>% filter(LAT<=(-26) & Estuary=="NO") %>%
      mutate(YEAR=YEAR.c,
             BLOCK=BLOCKX,
             Catch=LIVEWT.c,
             BLOCK=ifelse(BLOCK>=96000,paste(floor(abs(LAT)),(floor(LONG)-100)*10,sep=""),BLOCK),
             BLOCK=substring(BLOCK,1,4)) %>%
      select(-c(ZnID,MonthlyID,ZoneID,YEAR.c,blockxFC,CONDITN,
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
      left_join(All.species.names%>%rename(Species=SPECIES),by=c('SPECIES'='SP'))%>%
  mutate(Name=ifelse(is.na(Name),COMMON_NAME,Name),
         Scien.nm=ifelse(is.na(Scien.nm),SCIENTIFIC_NAME,Scien.nm))
DATA_obs=list(GN=Dat_obs) #Use only observed GN as LL has very few observations

#define discarded/retained sp
for(i in 1:length(DATA_obs)) DATA_obs[[i]]$Discarded=with(DATA_obs[[i]],
                                          ifelse(SPECIES%in%Commercial.sp,"Retained","Discarded"))

#Show overal observed discarded and retained species by year   #nice infographic: https://twitter.com/ISSF/status/1147943595942526978?s=03
fun.horiz.bar=function(d)
{
  d=d%>%
    mutate(dummy=ifelse(Discarded=="Retained","Retained",Name))%>%
    group_by(dummy)%>%
    tally%>%
    mutate(Percent=100*n/sum(n),
           Label=ifelse(Percent>1,'','<1%'))%>%
    data.frame%>%
    arrange(-Percent)
  d%>%
  ggplot(aes(x =  reorder(dummy, Percent),y = Percent)) + 
    geom_bar(colour = "black",stat = "identity")+coord_flip()+
    geom_text(aes(label=Label), position=position_dodge(width=0.9), hjust=-0.5)+
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  return(d)

}
for(i in 1:length(DATA_obs))  Tabl.obs.GN=fun.horiz.bar(d=DATA_obs[[i]])
ggsave('Results/Figure2_observed.percentage_gillnet.tiff', width = 8,height = 10,compression = "lzw")
Discarded.SP=DATA_obs$GN%>%
              filter(Discarded=='Discarded')%>%
              distinct(SPECIES,.keep_all = T)%>%
              pull(SPECIES)


# Fill in missing FL info  ---------------------------------------------------------
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
}
  
  #2. Derive FL as a proportion of TL
for(i in 1:length(DATA_obs)) DATA_obs[[i]]$FL=with(DATA_obs[[i]],ifelse(is.na(FL),TL*.875,FL))

  #3. samples from species distribution
fn.samp.dist=function(d)
{
  d=subset(d,FL>20)
  his=hist(d$FL,breaks=seq(floor(min(d$FL,na.rm=T)),ceiling(max(d$FL,na.rm=T)),by=1),plot=F)
  return(his)
}
for(i in 1:length(DATA_obs))
{
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
  }
}

  #4. overall mean (only 2 records, Crocodile shark & Grey reef shark)
for(i in 1:length(DATA_obs)) DATA_obs[[i]]$FL[is.na(DATA_obs[[i]]$FL)]=mean(DATA_obs[[i]]$FL,na.rm=T)

pdf("Results/Preliminary/size.frequency_ammended.pdf")
for(u in 1:length(UniK.sp))
{
  a=subset(DATA_obs$GN,SPECIES==UniK.sp[u] & !is.na(FL))
  if(nrow(a)>2)hist(a$FL,main=a$Name[1],col=2, xlab="FL (cm)")
}
dev.off()

# Convert numbers to weight  ---------------------------------------------------------
#note: calculate discard ratio using weight because total catch is in weight.
#       Use FL and a, b parameters for fork length
Len.wei=Len.wei%>%
          select(c(SPECIES,a_w8t,b_w8t,Source,Comments))
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]=DATA_obs[[i]]%>%
                 left_join(Len.wei,by="SPECIES")%>%
                 mutate(Catch=a_w8t*(FL)^b_w8t)

}

# Remove white shark and 'other shark' from observations ---------------------------------------------------------
#note: 2 very large white sharks acounts for a great proportion of discarded tonnage
for(i in 1:length(DATA_obs))  DATA_obs[[i]]=DATA_obs[[i]]%>%filter(!SPECIES%in%c('WP','XX'))

#Reset min and maximum weights is nonsense
Wei.range=Wei.range%>%left_join(Wei.range.names,by=c("Sname"))
Wei.range=subset(Wei.range,!(is.na(SPECIES)|is.na(TW.min)|is.na(TW.max)))
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]=left_join(DATA_obs[[i]],Wei.range,by=c('Species'='SPECIES'))
  
  DATA_obs[[i]]=DATA_obs[[i]]%>%
                  mutate(Catch=case_when(!is.na(TW.max) & Catch> TW.max ~ TW.max,
                                         !is.na(TW.min) & Catch< TW.min ~ TW.min,
                                         TRUE~Catch))
}
  



# group rare species  ---------------------------------------------------------
for(i in 1:length(DATA_obs))
{
  Tab.sp=prop.table(with(subset(DATA_obs[[i]],Discarded=="Discarded"),table(SPECIES)))
  Rare.sp=names(Tab.sp[Tab.sp<Group_rare_criteria])
  DATA_obs[[i]]$SPECIES.ori=DATA_obs[[i]]$SPECIES
  DATA_obs[[i]]$SPECIES=with(DATA_obs[[i]],ifelse(SPECIES%in%Rare.sp,"Grouped",SPECIES))
}

# Define discarded_sp used to calculate discarded species proportions-------------------------------------------
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]$Discarded_sp=with(DATA_obs[[i]],ifelse(Discarded=="Retained","Retained",SPECIES))
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

#remove stuff
rm(DATA,Dat_obs,Dat_obs.LL,Dat_total,Dat_total.LL,DATA.bio,DATA.ecosystems)


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
Obs_ratio.strata=STRATA_obs(d=DATA_obs$GN,
                            Strata=c("Discarded_sp",STRTA.obs.disc),
                            Min.obs.per.block=Min.obs.per.block,
                            Min.shots.per.block=Min.shots.per.block)
write.csv(Obs_ratio.strata$dat,"Results/Table1_observed.ratios.csv",row.names = F)


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
  #from commercial catch reconstruction
PCM=rbind(PCM.north,PCM.south)%>%
  distinct(Name,.keep_all = T)%>%
  dplyr::select(-c(Trawl,LL))%>%
  rename(PCM=GN)%>%
  left_join(All.species.names%>%dplyr::select(-c(Scien.nm,Family,SPECIES)),by='Name')%>%
  rename(SPECIES=SP)

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
  filter(SPECIES%in%Discarded.SP)


  #add PCM for remaining discarded speices
    #Angel Shark, from Braccini et al 2012; 
    #Brown-banded catshark, set at average of Rusty and Varied carpetshark from Braccini et al 2012
    #Crocodile shark, no information available, set = to mako shark from Braccini et al 2012 as precaution
    # Southern fiddler ray no available information, set at 10% due to hardiness
    #Sliteye shark, set at 1, always dead when brought aboard
    #Guitarfish & shovelnose ray, set equal to Whitespot shovelnose 
    #Spotted shovelnose, set equal to Whitespot shovelnose
    #Dwarf spotted wobbegong, set equal to Cobbler Wobbegong from Comm.disc.sp
    #White shark, set = to mako shark from Braccini et al 2012 
    #Whitespot shovelnose, Braodhurst & Cullis 2020 (immediate mortality)+25%
    #Cobbler wobbegong, set equal to Orectolobus parvimaculatus
    #Spikey dogfish, set equal to Squalus spp.
    # Stingrays, set equal to Myliobatidae
PCM.remaining=data.frame(
  Name=c("Angel Shark","Brown-banded catshark","Crocodile shark",            
         "Southern fiddler ray","Sliteye shark","Guitarfish & shovelnose ray",
         "Spotted shovelnose","Dwarf spotted wobbegong","White shark",               
         "Whitespot shovelnose","Other shark",
         'Cobbler wobbegong','Spikey dogfish','Stingrays'),
  Scien.nm=c("Squatina australis","Chiloscyllium punctatum","Pseudocarcharias kamoharai",
             "Trygonorrhina dumerilii","Loxodon macrorhinus","Families Rhinobatidae & Rhynchobatidae",
             "Aptychotrema timorensis","Orectolobus parvimaculatus","Carcharodon Carcharias",
             "Rhynchobatus australiae","",
             "Sutorectus tentaculatus","Squalus megalops","Dasyatidae"),
  PCM=c(1-0.595,1-mean(c(0.771,0.687)),1-0.242,
        1-.9,1,(0.29+0.29*.25),
        (0.29+0.29*.25),0.057,1-0.242,
        (0.29+0.29*.25),NA,
        0.0570,0.1729,0.1460),
  SPECIES=c("AU","BC","CR","FR","SE","SH","SS","WM","WP","WR","XX",
            "PD","SR","WC"))



PCM=rbind(PCM,PCM.remaining)%>%
  filter(!Name=='Spur Dog')

  

# Annual discards by species (sum(D_i_h x P_i x Y_h)) -------------------------------------------
fn.total=function(disc.dat,tot.dat,pcm)
{
  #combine observed ratios and total reported catch
  Total.discard=merge(disc.dat$dat,tot.dat,by=STRTA.reported.ret,all.y=T)
  
  #linear interpolate missing block discard ratio
  Vars=disc.dat$vars
  for(v in 1:length(Vars)) 
  {
    id=match(Vars[v],names(Total.discard))
    Total.discard[,id]=ifelse(is.na(Total.discard[,id]),
              na.approx(zoo(Total.discard[,id])),Total.discard[,id]) 
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
      dummy[,id]=ifelse(is.na(dummy[,id]),na.approx(zoo(dummy[,id])),dummy[,id]) 
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
    Total.discard[,this]=Total.discard[,this]*this.pcm
  }
  
  return(Total.discard)
}
Tot.discard.result=fn.total(disc.dat=Obs_ratio.strata,
                            tot.dat=Total_strata,
                            pcm=PCM)

# Sensitivity tests ---------------------------------------------------------
if(do.sensitivity)
{
  #S1: PCM set at 1.5 base case
  PCM.S1=PCM%>%
    mutate(PCM=1.5*PCM)
  S1=fn.total(disc.dat=Obs_ratio.strata,
              tot.dat=Total_strata,
              pcm=PCM.S1)
  
  #S2: Min obs per block = 20 & min shots per block = 10 
  Obs_ratio.strata.S2=STRATA_obs(d=DATA_obs$GN,
                                 Strata=c("Discarded_sp",STRTA.obs.disc),
                                 Min.obs.per.block=20,
                                 Min.shots.per.block=10)
  S2=fn.total(disc.dat=Obs_ratio.strata.S2,
              tot.dat=Total_strata,
              pcm=PCM)
  
  
  #S3: Min obs per block = 5 & min shots per block = 2
  Obs_ratio.strata.S3=STRATA_obs(d=DATA_obs$GN,
                                 Strata=c("Discarded_sp",STRTA.obs.disc),
                                 Min.obs.per.block=5,
                                 Min.shots.per.block=2)
  S3=fn.total(disc.dat=Obs_ratio.strata.S3,
              tot.dat=Total_strata,
              pcm=PCM)
  #ACA.
  #S4: unobserved blocks imputed from random sample of observed blocks
}

#show blocks interpolated by species
fn.number.interpolated=function(disc.dat,tot.dat)
{
  tot.dat=tot.dat[!duplicated(tot.dat$BLOCK),]
  
  #combine observed ratios and total reported catch
  Total.discard=merge(disc.dat$dat,tot.dat,by=STRTA.reported.ret,all.y=T)
  
  #count NAs by species
  return(Total.discard %>% select (-c(BLOCK, total.retained,YEAR)) %>%
           summarise_all(list(~sum(is.na(.))))%>%
           as.data.frame)
}
interpolated.blocks=fn.number.interpolated(disc.dat=Obs_ratio.strata,tot.dat=Total_strata)
write.csv(interpolated.blocks,"Results/Table2_interpolated.blocks.csv",row.names=F)


# Uncertainty thru non-parametric bootstrap  -------------------------------------------
#note: Takes 0.002 secs per iteration
cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)
clusterCall(cl, function() {
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(zoo)
})
Results=Tot.discard.result
set.seed(666)
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
    Tot.discard.boot=fn.total(disc.dat=Obs_ratio.strata.boot,tot.dat=Total_strata)
    
    rm(Dat_obs.boot)
    return(Tot.discard.boot)
  }
  Results=store
  rm(d)

})
rm(store)
stopCluster(cl) 




# Report ---------------------------------------------------------

#Table of PCM and weights
Table_S.2=left_join( Len.wei,PCM,by='SPECIES')
write.csv(Table_S.2,"Results/Table_S.2.csv",row.names=F)


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
})

# Stats
PRBS=c(.025,.5,.975)
Stat.list=list(Species="Species",Total="Total")
Stats=Stat.list
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
for(l in 1:length(Stats)) Stats[[l]]=fn.stats(d=Results.show,how=Stat.list[[l]])


#Plot statistics
yrs=sort(unique(DATA_total$GN$YEAR))
Col.totl="black"
Col.totl.retained="grey60"

#plotting function
id.med=match("50%",rownames(Stats$Species[[2]]))
id.low=match("2.5%",rownames(Stats$Species[[2]]))
id.up=match("97.5%",rownames(Stats$Species[[2]]))


#Total
fn.plt.prop=function(d,d.ret,CL,LWd)
{
  plot(yrs,yrs,col='transparent',ylim=c(0,YMAX),ann=F)
  polygon(x=c(yrs,rev(yrs)),y=c(d[id.low,]/d.ret,rev(d[id.up,]/d.ret)),
          col=adjustcolor(CL, alpha.f = 0.30), border = adjustcolor(CL, alpha.f = 0.60))
  lines(yrs,d[id.med,]/d.ret,lwd=LWd,col=CL)
}
YMAX=1
#YMAX=max(c(unlist(Stats$Total$total.retained),unlist(Stats$Total$total.discarded)))
jpeg("Results/Total.jpg",width=2400,height=2400,units="px",res=300)
par(mfcol=c(1,1),mar=c(1.5,1,1,1),oma=c(1,3,1,1),mgp=c(1,.6,0),las=1)
fn.plt.prop(d=Stats$Total$total.discarded,d.ret=Stats$Total$total.retained[id.med,],CL=Col.totl,LWd=3)
mtext("Proportion",2,1.5,las=3,outer=T,cex=1.5)
#lines(yrs,Stats$Total$total.retained[id.med,],lwd=3,col=Col.totl.retained)
#legend("topright",c("retained","discarded"),lty=1,lwd=3,
#       bty='n',col=c(Col.totl.retained,Col.totl),cex=1.5)
#mtext("Total (tonnes)",2,1.5,las=3,outer=T,cex=1.5)
mtext("Year",1,0,outer=T,cex=1.5)
dev.off()

#By species
fn.plt=function(d,CL,LWd)
{
  plot(yrs,yrs,col='transparent',ylim=c(0,YMAX),ann=F)
  polygon(x=c(yrs,rev(yrs)),y=c(d[id.low,],rev(d[id.up,])),
          col=adjustcolor(CL, alpha.f = 0.30), border = adjustcolor(CL, alpha.f = 0.60))
  lines(yrs,d[id.med,],lwd=LWd,col=CL)
}
jpeg("Results/Figure3_by.species.jpg",width=2400,height=2400,units="px",res=300)
Plt.this=Stats$Species[-match(c("YEAR","total.retained"),names(Stats$Species))]
Plt.this=Plt.this[sort(names(Plt.this))]
smart.par(n.plots=length(Plt.this),MAR=c(1,1.5,1,1),OMA=c(2,2,.1,.1),MGP=c(1,.6,0))
for(i in 1:length(Plt.this))
{
  YMAX=max(unlist(Plt.this[[i]]))
  fn.plt(d=Plt.this[[i]],CL=Col.totl,LWd=1.5)
  mtext(unlist(strsplit(names(Plt.this)[i], ".", fixed = TRUE))[2],3,-1.51,cex=1) 
}
mtext("Year",1,0.7,outer=T,cex=1.25)
mtext("Total discard (tonnes)",2,0.5,outer=T,las=3,cex=1.25)
dev.off()


# Export total discard estimates-----------------------------------------------------------------------
#MISSING: export to where comm catch recons export (i.e. each species folder!!!)
  #total
colnames(Stats$Total$total.discarded)=yrs
write.csv(Stats$Total$total.discarded,
      paste("Results/Total.discard.estimates/Total.csv",sep=""))

  #by species
for(i in 1:length(Plt.this)) 
{
  nm=unlist(strsplit(names(Plt.this)[i], ".", fixed = TRUE))[2]
  colnames(Plt.this[[i]])=yrs
  write.csv(Plt.this[[i]],paste("Results/Total.discard.estimates/",nm,".csv",sep=""))
}
  
