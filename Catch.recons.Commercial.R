#-- Script for reconstructing time series of commercial catch of sharks in WA and IUU

# Annual updates:
# Get annual effort updates for these fisheries:
#           XXX

#Sources of catch:
#     Taiwanese pelagic gillnet (1974-1986) and longline (1989-91) catch (Stevens & Davenport 1991; Stevens 1999)


library(tidyverse)
library(readxl)
library(lubridate)

Asses.year=2019    #enter year of assessment
Do.recons.com.fishn.paper="NO"
Historic.yrs=1941:1975

#Catch reconstruction scenarios
ScenarioS=data.frame(Scenario=c('Base Case','High','Low'),
                     PCM=c(1,1.5,.5),
                     Weight=c(1,1.5,.5))

# 1 -------------------DATA SECTION------------------------------------
options(stringsAsFactors = FALSE)

#Shark bio data
User="Matias"
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R")

#species codes
All.species.names=read.csv("C:/Matias/Analyses/Population dynamics/1.Other species/Species_names.csv")
b=read.csv("C:\\Matias\\Data\\Species.code.csv")


#Catch
setwd("C:/Matias/Analyses/Data_outs")

#1.1 Commercial catch shark fisheries TDGDLF and NSF
Data.monthly=read.csv("Data.monthly.csv")
Data.monthly.north=read.csv("Data.monthly.NSF.csv")
Data.monthly.north$LIVEWT.c=Data.monthly.north$LIVEWT  


#1.2 Commercial shark catch reported in other WA fisheries
Data.request.1=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/1975_1987.csv")
Data.request.2=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/1988_2002.csv")
Data.request.3=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/2003_2016.csv")
Data.request=rbind(Data.request.1,Data.request.2,Data.request.3)
Data.request=subset(Data.request,!METHOD%in%c("GN","LL"))
names(Data.request)[match("LIVEWT",names(Data.request))]="LIVEWT.c"


#1.3 Historic (pre 1975)
  #1.3.1 Whitley 1944 (citation in Simpfendorfer & Donohue 1998)
#description: All shark catch across southern WA (south of Abrolhos islands). Weight in pounds
Catch_1941_1943=data.frame(year=1941:1943,LIVEWT.c=c(55332,77441,109064)) 
Catch_1941_1943$LIVEWT.c=Catch_1941_1943$LIVEWT.c*0.453592/1000  #Convert to tonnes

  #1.3.2 Heald 1987
#description: All shark catch by port. Weight in kgs
Catch_1949=read.csv("C:/Matias/Data/Catch and Effort/Heald(1987)_1949_total_catch_by port_allspecies.csv") 
Catch_1949$LIVEWT.c=Catch_1949$Shark.catch.kg_live_wgt./1000   #in tonnes

  #1.3.3 Simpfendorfer & Donohue 1998
#description: All shark catch for southern WA (Geraldton to SA border). Weight in tonnes
Catch_1950=data.frame(year=1950,LIVEWT.c=50)   

  #1.3.4 Heald 1987
#description= All shark catch for southern WA (Geraldton to SA border). Weight in tonnes
Catch_1952_1975=read.csv("C:/Matias/Data/Catch and Effort/Historical_WA_shark_catches_1952_1975.csv") 
names(Catch_1952_1975)=c("year","LIVEWT.c")


#1.4 Taiwanese gillnet fishery (source Stevens 1999 and Stevens & Davenport 1991, catch in tonnes)
Taiwan.gillnet.ktch=data.frame(Year=1974:1986,
                               Total.Ktch=c(321.4,8997.6,6455.3,9970.7,5528.2,3282.1,5831.1,6694.7,
                                            5624.1,7589.9,6544.2,2929.5,2111.1),
                               Shark.Ktch=c(rep(NA,5),557.9,4234.5,4486.2,3639.4,4418.6,2731.4,2327.4,2393.9),
                               WA.Ktch=c(rep(NA,5),50,1250,720,800,790,400,10,70))  #From Figure 7 Stevens & Davenport
Mn.ktch=mean(Taiwan.gillnet.ktch$WA.Ktch/Taiwan.gillnet.ktch$Shark.Ktch,na.rm=T)
Taiwan.gillnet.ktch=Taiwan.gillnet.ktch%>%
  mutate(Shark.Ktch=ifelse(is.na(Shark.Ktch),Total.Ktch*mean(Shark.Ktch/Total.Ktch,na.rm=T),
                           Shark.Ktch),
         WA.Ktch=ifelse(is.na(WA.Ktch),Shark.Ktch*mean(WA.Ktch/Shark.Ktch,na.rm=T),WA.Ktch))

Taiwan.gillnet.sp.comp=as.data.frame(matrix(c(.357,.166,.059),ncol=3))
colnames(Taiwan.gillnet.sp.comp)=c("Australian.blacktip.Shark.West.prop",
                                   "Spot.tail.Shark.West.prop","Hammerheads.West.prop")
Taiwan.gillnet.ktch=cbind(Taiwan.gillnet.ktch,Taiwan.gillnet.sp.comp) %>%
  mutate('Australian blacktip shark'=Australian.blacktip.Shark.West.prop*WA.Ktch,
         'Spot-tail shark'=Spot.tail.Shark.West.prop*WA.Ktch,
         Hammerheads=Hammerheads.West.prop*WA.Ktch)%>%
  dplyr::select(Year,'Australian blacktip shark','Spot-tail shark',Hammerheads)%>%
  gather(Species,Ktch,-Year)

Taiwan.longline.ktch=data.frame(Year=1990:1991,Shark.Ktch=c(1700*11/(11+9),1700*9/(11+9)))%>%
  mutate(WA.ktch=Shark.Ktch*Mn.ktch)

Taiwan.longline.sp.comp=data.frame(Species=rep(c("Spot-tail shark","Australian blacktip shark","Tiger shark",
                                                 "Milk shark","Spinner shark",
                                                 "Pigeye shark","Graceful shark","Hammerheads","Other"),2),
                                   Year=c(rep(1990,9),rep(1991,9)), 
                                   Percent=c(18.3,25.6,8.1,.4,6.9,17.9,1.2,2,19.5,5.6,79.7,0.7,2,4.5,1,0,0.1,6.4))
Taiwan.longline.ktch=Taiwan.longline.sp.comp%>%
  left_join(Taiwan.longline.ktch,by="Year")%>%
  mutate(Ktch=WA.ktch*Percent/100)%>%
  dplyr::select(names(Taiwan.gillnet.ktch))


#1.5 Wetline  
  #Western Rock lobster
WRL=read.csv("C:/Matias/Data/Catch and Effort/WRL/Number-of-vessels.csv")
WRL.Wann=read.csv("C:/Matias/Data/Catch and Effort/WRL/Wann_catch.csv")


#1.6 Commonwealth Southern and Western Tuna and Billfish Fisheries
AFMA_catch=read.csv("Catch and Effort/AFMA_catch.csv",stringsAsFactors=F)  #source: Bensely et al 2010. Appendix A
Effort_WTBF=read.csv("Catch and Effort/AFMA full effort series.csv",stringsAsFactors=F)
#Effort_WTBF=read.csv("Catch and Effort/WTBF_effort.csv",stringsAsFactors=F)  #source: Bensely et al 2010. Appendix A. this is effort in Commonwealth waters
cpue_dusky_WTBF=37/203205  #number of individuals per 203205 hooks observed.97% discarded alive Source:Stobutski et al 2006 Table 3
cpue_sandbar_WTBF=8/203205   #100% discarded alive
#page 43 Borg & McAuley 2004: 1165 individuals of bronzey in 2001, also catches for previous years... DISCUSS!


#1.7 Pilbara trawl (Corey W.)
#44% of observed shots (n=2000) have carcharhinids reaching deck. Unknown post-capture survival and 
# unknown species composition. This could be resolved if videos reviewed...

#1.8 TEPs   
  #1.8.1 TDGDLF
#description: from Comments in TDGDLF returns
# Oversized dusky sharks   #MISING: aDD 2011-12 AND 2012-13 YEARS!!!
Overzd.Whler=read.csv("Catch and Effort/TEPS/Comments from SharkDailyTripReturns.csv",stringsAsFactors=F)

  #1.8.2 Other fisheries
TEPS_other=read.csv("Catch and Effort/TEPS_catch.csv",stringsAsFactors=F)

  #1.X Greynurse catches since protection           REVIEW in light of TEP stuff for dusky....
GN.mn.wght=25
GN.pcm=.2
Under.rep.factor=2  #from white shark catch recons
grey.hndl=paste('C:/Matias/Analyses/Catch and effort/State of fisheries/',
                paste(Asses.year-2,substr(Asses.year-1,3,4),sep='-'),sep='')
Greynurse.ktch=read.csv(paste(grey.hndl,'3.Table2.TEPS.csv',sep='/'),stringsAsFactors = F)
Greynurse.ktch=Greynurse.ktch%>%
  filter(Species=="SHARK, GREY NURSE")
Grey.nms=colnames(Greynurse.ktch)[-1]
Al.ded=c(NA,sapply(strsplit(Grey.nms, "[.]"), "[", 1))
Greynurse.ktch=data.frame(finyear=unique(as.numeric(substr(sapply(strsplit(Grey.nms, "[.]"), "[", 2),1,4))),
                          N.alive=unlist(Greynurse.ktch[,which(Al.ded=="alive")]),
                          N.dead=unlist(Greynurse.ktch[,which(Al.ded=="dead")]))%>%
  mutate(LIVEWT.c=Under.rep.factor*(N.dead*GN.mn.wght+N.alive*GN.mn.wght*GN.pcm))%>%
  dplyr::select(finyear,LIVEWT.c)



# 2 -------------------PROCEDURE SECTION------------------------------------

#2.1 #combine historic in single dataframe
Catch_1949.total=data.frame(year=1949,LIVEWT.c=sum(Catch_1949$LIVEWT.c))     
Catch_1949.total$LIVEWT.c=Catch_1950$LIVEWT.c     
Historic.ktch=rbind(Catch_1941_1943,Catch_1949.total,Catch_1950,Catch_1952_1975)

#add missing years through linear interpolation
Missing.yrs=Historic.yrs[which(!Historic.yrs%in%Historic.ktch$year)]
Missing.ktch=approx(Historic.ktch$year,Historic.ktch$LIVEWT.c,xout=Missing.yrs)   

fn.fig("C:/Matias/Analyses/Catch and Effort/Pre_1975_total_shark_catch",2400, 2400) 
par(mfcol=c(1,1),las=1,mgp=c(2.9,.6,0))
with(subset(Historic.ktch,year<1975),plot(year,LIVEWT.c, ylab="Total shark landings (tonnes)",
                                          xlab="Financial year",
                                          pch=19,col=1,cex=2,cex.lab=1.7,cex.axis=1.35))
#  points(Missing.ktch$x,Missing.ktch$y,pch=19,col=2,cex=2)
# legend("topleft",c("reported","interpolated"),pch=19,col=1:2,cex=1.5,pt.cex=2,bty='n')
dev.off()

Historic.ktch=rbind(Historic.ktch,data.frame(year=Missing.ktch$x,LIVEWT.c=Missing.ktch$y))
Historic.ktch=Historic.ktch[order(Historic.ktch$year),]
Historic.ktch$LIVEWT.c=Historic.ktch$LIVEWT.c*1000   #convert to kg

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

#allocate historic catch to species
Hist.expnd=expand.grid(year=Historic.ktch$year,SPECIES=Prop.sp.yr.all$SPECIES)
Hist.expnd=Hist.expnd%>%
  left_join(Historic.ktch,by="year")%>%
  left_join(Prop.sp.yr.all,by='SPECIES')%>%
  rename(ktch=LIVEWT.c)%>%
  mutate(LIVEWT.c=ktch*Proportion)


#2.2 remove data for other gears that is already recorded in Data.monthly
TDGDLF.other.gears=Data.monthly%>%
  filter(!METHOD%in%c("GN","LL"))%>%
  mutate(dummy=paste(FINYEAR,MONTH,VESSEL,BLOCKX,METHOD))
NSF.other.gears=Data.monthly.north%>%
  filter(!METHOD%in%c("GN","LL"))%>%
  mutate(dummy=paste(FINYEAR,MONTH,VESSEL,BLOCKX,METHOD))
fn.subs=function(YEAR) substr(YEAR,start=3,stop=4)
Data.request=Data.request%>%
  mutate(FINYEAR=ifelse(MONTH%in%1:6,paste(YEAR-1,"-",fn.subs(YEAR),sep=""),
                        ifelse(MONTH%in%7:12,paste(YEAR,"-",fn.subs(YEAR+1),sep=""),NA)),
         LAT=-as.numeric(substr(BLOCKX,1,2)),
         LONG=100+as.numeric(substr(BLOCKX,3,4)),
         zone=as.character(ifelse(LONG>=116.5 & LAT<=(-26),"Zone2",
                           ifelse(LONG<116.5 & LAT<=(-33),"Zone1",
                           ifelse(LAT>(-33) & LAT<=(-26) & LONG<116.5,"West",
                           ifelse(LAT>(-26) & LONG<114,"Closed",
                           ifelse(LAT>(-26) & LONG>=114 & LONG<123.75,"North",
                           ifelse(LAT>(-26) & LONG>=123.75,"Joint",NA))))))),
         dummy=paste(FINYEAR,MONTH,VESSEL,BLOCKX,METHOD))%>%
  filter(!dummy%in%c(TDGDLF.other.gears$dummy,NSF.other.gears$dummy))
Data.request.south=Data.request%>%filter(LAT<=max(Data.monthly$LAT))
Data.request.north=Data.request%>%filter(LAT>max(Data.monthly$LAT))
rm(Data.request,Data.request.1,Data.request.2,Data.request.3)


#ACA

#2.3 Wetline catches of dusky sharks                   
  #rock lobster
#note: extrapolate Wann's catch (for one season) to the entire fishery as a proportion of the number of boats
WRL.Wann$TL=WRL.Wann$TL_metres*100    #convert to cm
WRL.Wann$TL=with(WRL.Wann,ifelse(is.na(TL),TL_feet*30.48,TL))
Dusky.WRL=subset(WRL.Wann,Species=="DW")
Dusky.WRL$LiveWt=fn.weight(Dusky.WRL$TL,bwt,awt)   #in kg
if(KTCH.UNITS=="TONNES") Annual.Dusky.Ktch.WRL.Wann=sum(Dusky.WRL$LiveWt)/1000 else #in tonnes
                         Annual.Dusky.Ktch.WRL.Wann=sum(Dusky.WRL$LiveWt)
Dusky.Tot.Ktch.WRL=WRL
Dusky.Tot.Ktch.WRL$LiveWt=WRL.prop*Annual.Dusky.Ktch.WRL.Wann*Dusky.Tot.Ktch.WRL$Number.of.vessels  #in tonnes
colnames(Dusky.Tot.Ktch.WRL)[2:3]=c("FINYEAR","LIVEWT.c")
Dusky.Tot.Ktch.WRL=Dusky.Tot.Ktch.WRL[,-1]   #remove number of vessels     
Dusky.Tot.Ktch.WRL=fn.exp.his(Dusky.Tot.Ktch.WRL,Monthly.prop)  #partition catch by month    



#2.4 SWTBF catches of dusky and sandbar sharks               
#Option=1       #select source to use to obtain catches
Option=2

#1.6.1 Commonwealth Southern and Western Tuna and Billfish Fisheries and GAB trawl fishery
if(KTCH.UNITS=="TONNES") for(i in 2:ncol(AFMA_catch))AFMA_catch[,i]=AFMA_catch[,i]/1000   #convert to tonnes

if(SP=="BW") 
{
  GAB_trawl=AFMA_catch[,c(1,2)]
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
    Dusky_WTBF[,2]=Effort_WTBF[,2]*PCM*cpue_dusky_WTBF*Mean.wt.dusky
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


#2.5  TEPS

#1.8.1 TDGDLF
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
  
  TEPS_other=subset(TEPS_other,CommonName%in%these.sp)
  #Note: it seems all reported TEPS are for the Shark fisheries...  don't use this and only using TEPS from TDGDFL
  
  TEPS_dusky=TEPS_TDGDLF[,c(1,5)]
  TEPS_dusky=fn.exp.his(TEPS_dusky,Monthly.prop)
  
}



# 3 -------------------REPORT SECTION------------------------------------
if(Do.recons.com.fishn.paper=="YES")
{
  
}
