#------------- Script for bringing in all shark data 

library(stringr)
library(tidyr)
library(Hmisc)
library(dplyr)
options(stringsAsFactors = FALSE) 

if(!exists('fn.word.table')) source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")
if(!exists('fn.fig')) fn.fig=function(NAME,Width,Height,Do.tiff="YES",Do.jpeg="NO")
{
  if(Do.tiff=="YES") tiff(file=paste(NAME,".tiff",sep=""),width=Width,height=Height,units="px",res=300,compression="lzw")
  if(Do.jpeg=="YES") jpeg(file=paste(NAME,".jpeg",sep=""),width=Width,height=Height,units="px",res=300)
}
if(!exists('smart.par')) smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))

fn.bring.all.data=function(path)
{
  if(path%in%list.files(Dat.repository))
  {
    temp.wd=paste(Dat.repository,path,sep='')
    Files=list.files(temp.wd)
    setwd(temp.wd)
    Files1=lapply(Files,read.csv)
    names(Files1)= str_remove(Files, path)
    return(Files1)
  }
}

fn.subs=function(YEAR) substr(YEAR,start=3,stop=4)

fn.finyr=function(Finyr)paste(Finyr-1,"-",fn.subs(Finyr),sep="")

fn.month.to.fin.mn=function(dat)ifelse(dat>6,dat-6,dat+6)

fn.exp.his=function(dat,dat1)
{
  dat2=expand.grid(FINYEAR=dat$FINYEAR,MONTH=1:12)
  dat2=merge(dat2,dat1[,match(c("MONTH","prop"),names(dat1))],by="MONTH",all.x=T)
  dat2=merge(dat2,dat,by="FINYEAR",all.x=T)
  dat2=dat2[order(dat2$FINYEAR,dat2$MONTH),]
  dat2$LIVEWT.c=dat2$LIVEWT.c*dat2$prop
  return(dat2[,-match("prop",names(dat2))])
}

fn.extract.dat=function(STRING,nm.Dat) grepl(STRING, nm.Dat, perl = TRUE)

fn.rename.dat=function(x,y) str_remove_all(x, paste(y, collapse = "|"))

fn.input.data=function(Name,Name.inputs,SP,Species,First.year,Last.year,Min.obs,Min.shts,What.Efrt,Bin.size,Yr.assess)   
{
  #----DATA SECTION------ 
  
  #1. Select years of data
  Used.yrs=as.numeric(c(substr(First.year,1,4),substr(Last.year,1,4)))
  Used.yrs=seq(Used.yrs[1],Used.yrs[2])
  Used.yrs=paste(Used.yrs,substr(Used.yrs+1,3,4),sep="-")

  #2. All catches
  catch=KtCh%>%filter(SPECIES%in%Species & FINYEAR%in%Used.yrs)
  catch.zone=KtCh.zone%>%filter(SPECIES%in%Species & FINYEAR%in%Used.yrs)
  
  #3. Bring in all species-specific data files
  nm.Dat=''
  Dat=fn.bring.all.data(path=capitalize(Name.inputs))
  if(!is.null(Dat)) nm.Dat=names(Dat)
  
    #3.1 RELATIVE ABUNDANCE    
      #3.1.1  TDGDLF  
  #monthly
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*annual.abundance.basecase)(?=.*relative)(?=.*monthly)",nm.Dat)]
  if(length(iid)>0)
  {
    Ab.indx.TDGDLF= Dat[match(iid,nm.Dat)]
    Nms=fn.rename.dat(x=names(Ab.indx.TDGDLF),y=c('.annual.abundance.basecase.monthly.','relative','.csv','_'))
    Nms=ifelse(Nms=='','all',Nms)
    names(Ab.indx.TDGDLF)=Nms
    
    Ab.indx.TDGDLF.all=Ab.indx.TDGDLF$all
    Ab.indx.TDGDLF=Ab.indx.TDGDLF[names(Ab.indx.TDGDLF) != "all"] 
    if(length(Ab.indx.TDGDLF)==0) rm(Ab.indx.TDGDLF) else
    {
      for(z in 1:length(Ab.indx.TDGDLF)) Ab.indx.TDGDLF[[z]]$zone=capitalize(names(Ab.indx.TDGDLF)[z])
      Ab.indx.TDGDLF=do.call(rbind,Ab.indx.TDGDLF)
    }
  }
  
  #daily
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*annual.abundance.basecase)(?=.*relative)(?=.*daily)",nm.Dat)]
  if(length(iid)>0)
  {
    Ab.indx.TDGDLF.daily=Dat[match(iid,nm.Dat)]
    Nms=fn.rename.dat(x=names(Ab.indx.TDGDLF.daily),y=c('.annual.abundance.basecase.daily.','relative','.csv','_'))
    Nms=ifelse(Nms=='','all',Nms)
    names(Ab.indx.TDGDLF.daily)=Nms
    
    Ab.indx.TDGDLF.all.daily=Ab.indx.TDGDLF.daily$all
    Ab.indx.TDGDLF.daily=Ab.indx.TDGDLF.daily[names(Ab.indx.TDGDLF.daily) != "all"] 
    if(length( Ab.indx.TDGDLF.daily)==0) rm( Ab.indx.TDGDLF.daily) else
    {
      for(z in 1:length(Ab.indx.TDGDLF.daily)) Ab.indx.TDGDLF.daily[[z]]$zone=capitalize(names(Ab.indx.TDGDLF.daily)[z])
      Ab.indx.TDGDLF.daily=do.call(rbind,Ab.indx.TDGDLF.daily)
    }
  }
  
      #3.1.2  NSF  
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*.annual.abundance.NSF_relative)",nm.Dat)]
  if(length(iid)>0) Ab.indx.NSF= Dat[match(iid,nm.Dat)]$.annual.abundance.NSF_relative.csv
  
      #3.1.3. Naturaliste survey
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*.Srvy.FixSt)",nm.Dat)]
  if(length(iid>0))
  {
    Srvy.FixSt=Dat[match(iid,nm.Dat)]
    Ab.index.Srvy.FixSt=Srvy.FixSt$.Srvy.FixSt.csv
    Size.index.Srvy.FixSt=Srvy.FixSt$.Srvy.FixSt_size.csv
    Ab.index.Srvy.FixSt$CV=Ab.index.Srvy.FixSt$CV/100
    Ab.index.Srvy.FixSt$FINYEAR=paste(Ab.index.Srvy.FixSt$yr,"-",fn.subs(Ab.index.Srvy.FixSt$yr+1),sep="")
    Size.index.Srvy.FixSt$FINYEAR=paste(Size.index.Srvy.FixSt$yr,"-",fn.subs(Size.index.Srvy.FixSt$yr+1),sep="")
  }
  

    #3.2 CATCH SIZE COMPOSITION
      #3.2.1 Heald 1987 
  # description: Partial length in cm obtained from Perth fish marke.
  #             Unspecified fishing method, most likely gillnets.
  # Not used, very unreliable and strange size measures and unsure about how measured
  #PL_Heald_1987=read.csv("C:/Matias/Data/Size_composition/Heald(1987).csv",stringsAsFactors=F)  
  
      #3.2.2 Stevens 1990 (citation in Simpfendorfer & Donohue 1998)
  # description: TL (fish market) or FL (measured at sea) depending on period, in cm.
  #             gillnet of 6.5 and 7 inch mesh size.
  # Not used, cannot be allocated to a zone and unsure about how measured
  #TL_FL_Stevens_1990=read.csv("C:/Matias/Data/Size_composition/Stevens_1990_size_comp_6.5_7_inch.csv",stringsAsFactors=F)    
  
      #3.2.3.TDGDLF observing programs
  #description: FL (cm) composition observed as part of different research projects on commercial gillnet vessels. 
  #           6.5 and 7 inch mesh combined (also available are data by mesh size). Souce: "Shark database"
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*inch.raw)",nm.Dat)]
  if(length(iid)>0)
  {
    FL.TDGDFL= Dat[match(iid,nm.Dat)]
    Nms=fn.rename.dat(x=names(FL.TDGDFL),y=c('_Size_composition_','.inch.raw','.csv','_'))
    Nms=ifelse(Nms=='','all',Nms)
    names(FL.TDGDFL)=Nms
    
    iid=nm.Dat[fn.extract.dat(STRING="(?=.*Numb_obs_size.freq.TDGDLF)",nm.Dat)]
    TDGDFL.size.numbers= Dat[match(iid,nm.Dat)]$`_Size_composition_Numb_obs_size.freq.TDGDLF.csv`
    
    if(nrow(do.call(rbind,FL.TDGDFL))==0)
    {
      rm(FL.TDGDFL,TDGDFL.size.numbers)
    }
  }
  
  #keep years with at least 10 observations from at least 10 shots by zone
  if(exists('TDGDFL.size.numbers'))
  {
    TDGDFL.size.numbers=subset(TDGDFL.size.numbers,N.observations>=Min.obs & N.shots>=Min.shts)
    This.yr.zn=with(TDGDFL.size.numbers,paste(FINYEAR,zone))
  }
  
      #3.2.4 Pilbara trawl
  #description: FL (cm) composition observed on Pilbara trawl vessels. Source: "Shark database"
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Pilbara_Trawl)",nm.Dat)]
  if(length(iid)>0) FL_Pilbara_trawl= Dat[match(iid,nm.Dat)]$`_Size_composition_Pilbara_Trawl.csv`
  
      #3.2.5 NSF longline
  #description: FL (cm) composition observed on NSF longline vessels. Source: "Shark database"
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*NSF.LONGLINE)",nm.Dat)]
  if(length(iid)>0) FL_NSF= Dat[match(iid,nm.Dat)]$`_Size_composition_NSF.LONGLINE.csv`

      #3.2.6 TEPS_TDGLDF
  if(Name=="dusky shark") Size.comp.TEPS_TDGLDF=c(3.05,4,3,3.5,3.5,3.5,3.5,3,3,3,4,3)     #from Comments in TDGDLF returns
  

    #3.3. TAGGING
      #3.3.1 Conventional
  #note: there's also info by block but very few observations at this level....
  
  #Individual based model
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Ind.based.mod)",nm.Dat)]
  if(length(iid)>0) Rel_rec_Conv.Tag= Dat[match(iid,nm.Dat)]$`_Con_tag_Ind.based.mod.csv`

  #at age
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Zn.rel)",nm.Dat)]
  if(length(iid)>0) Zn.rel_Conv.Tag= Dat[match(iid,nm.Dat)]$`_Con_tag_Zn.rel_Conv.Tag.csv`
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Zn.rec)",nm.Dat)]
  if(length(iid)>0) Zn.rec_Conv.Tag= Dat[match(iid,nm.Dat)]$`_Con_tag_Zn.rec_Conv.Tag.csv`
  
  #at size
    #all sizes
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Zn.rel)",nm.Dat)]
  if(length(iid)>0) Zn.rel_Conv.Tag_size= Dat[match(iid,nm.Dat)]$`_Con_tag_Zn.rel_Conv.Tag_size.csv`
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Zn.rec)",nm.Dat)]
  if(length(iid)>0) Zn.rec_Conv.Tag_size= Dat[match(iid,nm.Dat)]$`_Con_tag_Zn.rec_Conv.Tag_size.csv`
  
      #adults and juvenlies
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Zn.rel.adul)",nm.Dat)]
  if(length(iid)>0) Zn.rel_Conv.Tag_size_adu= Dat[match(iid,nm.Dat)]$`_Con_tag_Zn.rel.adul_Conv.Tag_size.csv`
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Zn.rel.juv)",nm.Dat)]
  if(length(iid)>0) Zn.rel_Conv.Tag_size_juv= Dat[match(iid,nm.Dat)]$`_Con_tag_Zn.rel.juv_Conv.Tag_size.csv`
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Zn.rec.adul)",nm.Dat)]
  if(length(iid)>0) Zn.rec_Conv.Tag_size_adu= Dat[match(iid,nm.Dat)]$`_Con_tag_Zn.rec.adul_Conv.Tag_size.csv`
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Zn.rec.juv)",nm.Dat)]
  if(length(iid)>0) Zn.rec_Conv.Tag_size_juv= Dat[match(iid,nm.Dat)]$`_Con_tag_Zn.rec.juv_Conv.Tag_size.csv`
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Con_tag)(?=.*Smallest.size)",nm.Dat)]
  if(length(iid)>0) Smallest_size_tagged= Dat[match(iid,nm.Dat)]$`_Con_tag_Smallest.size_Conv.Tag_size.csv`
  

      #3.3.2 Acoustic
  #note: there's also info by block but very few observations at this level....
  
  #Taylor 2011 approach
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Acous.Tag)(?=.*Zn.rel_Acous.Tag)",nm.Dat)]
  if(length(iid)>0) Zn.rel_Acous.Tag= Dat[match(iid,nm.Dat)]$`_Acous.Tag_Zn.rel_Acous.Tag.csv`
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Acous.Tag)(?=.*Zn.rec_Acous.Tag)",nm.Dat)]
  if(length(iid)>0) Zn.rec_Acous.Tag= Dat[match(iid,nm.Dat)]$`_Acous.Tag_Zn.rec_Acous.Tag.csv`
  
  #Proportion of time approach
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Acous.Tag)(?=.*Zn.rel_Acous.Tag.prop)",nm.Dat)]
  if(length(iid)>0) Zn.rel_Acous.Tag.prop= Dat[match(iid,nm.Dat)]$`_Acous.Tag_Zn.rel_Acous.Tag.prop.csv`
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Acous.Tag)(?=.*Zn.rec_Acous.Tag.prop)",nm.Dat)]
  if(length(iid)>0) Zn.rec_Acous.Tag.prop= Dat[match(iid,nm.Dat)]$`_Acous.Tag_Zn.rec_Acous.Tag.prop.csv`

  #Individual based model
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Acous.Tag)(?=.*Ind_based)",nm.Dat)]
  if(length(iid)>0) Indiv_based_Acous.Tag= Dat[match(iid,nm.Dat)]$`_Acous.Tag_Acous.Tag.Ind_based.csv`
  
  #Reported recaptures and releases of acoustic tagging
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*Acous.Tag)(?=.*Rep.Recap)",nm.Dat)]
  if(length(iid)>0) Rep.Recap= Dat[match(iid,nm.Dat)]$`_Acous.Tag_Rep.Recap.csv`
  

    #3.4.  Age and growth data for whiskery only
  if(SP=="WH") Age.growth=read.csv("C:/Matias/Data/Age and growth/Simpfen.data.csv",stringsAsFactors=F)
  if(SP=="GM") Age.growth=read.csv("C:/Matias/Data/Age and growth/Terry/Gummy_Terry.csv",stringsAsFactors=F)
  if(SP=="BW") Age.growth=read.csv("C:/Matias/Data/Age and growth/Dusky.csv",stringsAsFactors=F)
  if(SP=="TK") Age.growth=read.csv("C:/Matias/Data/Age and growth/Sandbar.csv",stringsAsFactors=F)
  
  
    #3.5.  Standardised Mean size
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*annual.mean.size)",nm.Dat)]
  if(length(iid)>0)
  {
    Avr.wt.yr.zn=Dat[match(iid,nm.Dat)]
    Nms=fn.rename.dat(x=names(Avr.wt.yr.zn),y=c('.annual.mean.size','relative','.csv','_'))
    Nms=ifelse(Nms=='','all',Nms)
    names(Avr.wt.yr.zn)=Nms
    Avr.wt.yr=Avr.wt.yr.zn$all
    Avr.wt.yr.zn=Avr.wt.yr.zn[names(Avr.wt.yr.zn) != "all"] 
    if(length(Avr.wt.yr.zn)==0) rm(Avr.wt.yr.zn) else
    {
      for(z in 1:length(Avr.wt.yr.zn)) Avr.wt.yr.zn[[z]]$zone=capitalize(names(Avr.wt.yr.zn)[z])
      Avr.wt.yr.zn=do.call(rbind,Avr.wt.yr.zn)
    }
  }
  

    #3.6. Effort by zone (GN plus long line equivalent effort)
  #(in 1000 km gn days)
  if(What.Efrt=="km.gn.days") Eff.zn=read.csv("C:/Matias/Analyses/Data_outs/Annual.zone.eff.days.csv",stringsAsFactors=F)
  #(in 1000 km gn hours)
  if(What.Efrt=="km.gn.hours") Eff.zn=read.csv("C:/Matias/Analyses/Data_outs/Annual.zone.eff.hours.csv",stringsAsFactors=F)
  
  
    #3.7. Proportional effort by mesh size
  Mesh.prop.eff=read.csv("C:/Matias/Analyses/Catch and effort/mesh.proportional.effort.csv",stringsAsFactors=F)
  Mesh.prop.eff.West=read.csv("C:/Matias/Analyses/Catch and effort/mesh.proportional.effort.West.csv",stringsAsFactors=F)
  Mesh.prop.eff.Zn1=read.csv("C:/Matias/Analyses/Catch and effort/mesh.proportional.effort.Zone1.csv",stringsAsFactors=F)
  Mesh.prop.eff.Zn2=read.csv("C:/Matias/Analyses/Catch and effort/mesh.proportional.effort.Zone2.csv",stringsAsFactors=F)
  

    #3.8. Gillnet selectivity  
  iid=nm.Dat[fn.extract.dat(STRING="(?=.*gillnet.selectivity)",nm.Dat)]
  if(length(iid)>0)
  {
    Gillnet.selectivity=Dat[match(iid,nm.Dat)]$`_gillnet.selectivity.csv`
    Gillnet.selectivity_len.age=Dat[match(iid,nm.Dat)]$`_gillnet.selectivity_len.age.csv`
  }
    
  

  #----PARAMETERS SECTIONS ------- 
  LH.par=read.csv('C:/Matias/Data/Life history parameters/Life_History.csv')%>%filter(SPECIES==Species)
  a.TL=LH.par$a_FL.to.TL
  b.TL=LH.par$b_FL.to.TL
  Max.TL=LH.par$Max.TL


  
  #------PROCEDURE SECTION-------
  
  #1. CATCH

  #-Put all catches in list
  catch=catch %>%
    mutate(
      Fishery = case_when(
        FishCubeCode=='WRL'~'WRL',
        FishCubeCode=='WTB'~'WTBF',
        FishCubeCode=='TEP'~'TEPS',
        FishCubeCode=='Recreational'~'Rec',
        FishCubeCode=='SA MSF'~'SA_marine',
        FishCubeCode=='GAB'~'GAB_trawl',
        FishCubeCode=='Indo'~'Indonesian',
        FishCubeCode=='Taiwan'~'Taiwanese',
        FishCubeCode=='Historic'~'Historic',
        FishCubeCode%in%c('JASDGDL','WCDGDL','C070','OAWC')~'TDGDLF',
        FishCubeCode%in%c('JANS','OANCGC','WANCS')~'NSF',   
        TRUE  ~ "Other"),
      Fishery=ifelse(Fishery%in%c('NSF')& Name%in%c("gummy shark","whiskery shark"),'TDGDLF',Fishery))
  n.fishry=sort(unique(catch$Fishery))
  iid=vector('list',length(n.fishry))
  names(iid)=n.fishry
  for(v in 1:length(n.fishry)) iid[[v]]=subset(catch,Fishery==n.fishry[v])
  catch=iid
  
  catch.zone=catch.zone %>%   
    mutate(
      Fishery = case_when(
        FishCubeCode=='WRL'~'WRL',
        FishCubeCode=='WTB'~'WTBF',
        FishCubeCode=='TEP'~'TEPS',
        FishCubeCode=='Recreational'~'Rec',
        FishCubeCode=='SA MSF'~'SA_marine',
        FishCubeCode=='GAB'~'GAB_trawl',
        FishCubeCode=='Indo'~'Indonesian',
        FishCubeCode=='Taiwan'~'Taiwanese',
        FishCubeCode=='Historic'~'Historic',
        FishCubeCode%in%c('JASDGDL','WCDGDL','C070','OAWC')~'TDGDLF',
        FishCubeCode%in%c('JANS','OANCGC','WANCS')~'NSF',   
        TRUE  ~ "Other"),
      zone=ifelse(Name%in%c("gummy shark","whiskery shark") & 
                  zone%in%c('Closed','North','Joint'),'West',zone),
      Fishery=ifelse(Name%in%c("gummy shark","whiskery shark") & 
                  Fishery%in%c('NSF'),'TDGDLF',Fishery))
  n.fishry=sort(unique(catch.zone$Fishery))
  iid=vector('list',length(n.fishry))
  names(iid)=n.fishry
  for(v in 1:length(n.fishry)) iid[[v]]=subset(catch.zone,Fishery==n.fishry[v])
  catch.zone=iid
  

  #2. Put catch size composition data together (grouped by financial year and zone)
  
    #2.1. convert FL to TL 
  FL_to_TL=function(FL,a,b) TL=FL*a+b
  if(exists('FL.TDGDFL'))
  {
    for(n in 1:length(FL.TDGDFL))
    {
      FL.TDGDFL[[n]]$TL=round(FL_to_TL(FL.TDGDFL[[n]]$FL,a.TL,b.TL)) 
    }
  }
  if(exists('FL_Pilbara_trawl')) FL_Pilbara_trawl$TL=round(FL_to_TL(FL_Pilbara_trawl$FL,a.TL,b.TL))
  if(exists('FL_NSF')) FL_NSF$TL=round(FL_to_TL(FL_NSF$FL,a.TL,b.TL))
   
    #2.2. do the aggregation
  fn.add.missing.size=function(dat,Min,Max,interval)
  {
    Rango=c(Min,Max)
    SEQ=seq(Rango[1],Rango[2],interval)
    dat=dat[order(dat$year),]
    dat$FL.bin=floor(dat$FL/interval)*interval
    dat$FL.bin=factor(dat$FL.bin,levels=SEQ)
    
    dat$FINYEAR=as.character(with(dat,
           ifelse(Month%in%1:6,paste(year-1,"-",fn.subs(year),sep=""),
           ifelse(Month%in%7:12,paste(year,"-",fn.subs(year+1),sep=""),NA))))
    
    if(Conv.cal.mn.to.fin.mn=="YES") dat$Month=fn.month.to.fin.mn(dat$Month)
    
    #   dat=aggregate(Number~FL.bin+FINYEAR+Type,dat,sum)
    #   if(!is.null(dat$Sex))dat=aggregate(Number~FL.bin+FINYEAR+Type+Sex,dat,sum)
    #   
    #   dat2=subset(dat,Number>0,select=c(FL.bin,FINYEAR,Number))
    #   Table=reshape(dat2,v.names = "Number", idvar = "FINYEAR",timevar = "FL.bin", direction = "wide")
    #   X=colnames(Table)[2:ncol(Table)]
    #   colnames(Table)[2:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
    #   ID=which(!levels(dat$FL.bin)%in%colnames(Table)[2:ncol(Table)])
    #   ADD=levels(dat$FL.bin)[ID]
    #   ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
    #   names(ADD1)=ADD
    #   Table=cbind(Table,ADD1)
    #   Table=Table[,c(1,match(levels(dat$FL.bin),names(Table)))]
    #   Table[is.na(Table)]=0 
    #   
    dat=aggregate(Number~FL.bin+Month+FINYEAR+Type+Sex,dat,sum)
    
    dat2=subset(dat,Number>0,select=c(FL.bin,FINYEAR,Number,Month,Sex))
    Table=reshape(dat2,v.names = "Number", idvar = c("FINYEAR","Month","Sex"),timevar = "FL.bin", direction = "wide")
    X=colnames(Table)[4:ncol(Table)]
    colnames(Table)[4:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
    ID=which(!levels(dat$FL.bin)%in%colnames(Table)[4:ncol(Table)])
    ADD=levels(dat$FL.bin)[ID]
    ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
    names(ADD1)=ADD
    Table=cbind(Table,ADD1)
    Table=Table[,c(1:3,match(levels(dat$FL.bin),names(Table)))]
    Table[is.na(Table)]=0 
    Table=Table[order(Table$Sex,Table$FINYEAR,Table$Month),]
    
    return(Table)
    
  }
  fn.add.missing.size2=function(dat,Min,Max,interval)
  {
    Rango=c(Min,Max)
    SEQ=seq(Rango[1],Rango[2],interval)
    dat=dat[order(dat$year),]
    dat$FL.bin=floor(dat$FL/interval)*interval
    dat$FL.bin=factor(dat$FL.bin,levels=SEQ)
    dat$Sex=dat$SEX
    dat$Number=1
    dat=subset(dat,!Sex=="U")
    
    dat=aggregate(Number~FL.bin+Month+FINYEAR+Sex,dat,sum)
    
    dat2=subset(dat,Number>0,select=c(FL.bin,FINYEAR,Number,Month,Sex))
    Table=reshape(dat2,v.names = "Number", idvar = c("FINYEAR","Month","Sex"),timevar = "FL.bin", direction = "wide")
    X=colnames(Table)[4:ncol(Table)]
    colnames(Table)[4:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
    ID=which(!levels(dat$FL.bin)%in%colnames(Table)[4:ncol(Table)])
    ADD=levels(dat$FL.bin)[ID]
    ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
    names(ADD1)=ADD
    Table=cbind(Table,ADD1)
    Table=Table[,c(1:3,match(levels(dat$FL.bin),names(Table)))]
    Table[is.na(Table)]=0 
    Table=Table[order(Table$Sex,Table$FINYEAR,Table$Month),]
    return(Table)
    
  }
  fn.add.missing.size3=function(dat,Min,Max,interval)
  {
    Rango=c(Min,Max)
    SEQ=seq(Rango[1],Rango[2],interval)
    Table=dat%>%
          rename(Sex=SEX)%>%
          filter(!is.na(Sex) | !Sex=="U")%>%
          mutate(TL.bin=floor(TL/interval)*interval)%>%
          group_by(TL.bin,FINYEAR,Month,Sex)%>%
          summarise(N=n())
   add.msn=SEQ[which(!SEQ%in%Table$TL.bin)]  
   if(length(add.msn)>0)
   {
     addd=Table[1:length(add.msn),]
     addd$N=0
     addd$TL.bin=add.msn
     Table=rbind(Table,addd)
   }
   Table=Table%>%
     spread(TL.bin,N,fill=0)%>%
     data.frame%>%
     arrange(Sex,FINYEAR,Month)%>%
     filter(!is.na(FINYEAR))
   
   colnames(Table)[-match(c("FINYEAR","Month","Sex"),colnames(Table))]=
     str_remove_all(colnames(Table)[-match(c("FINYEAR","Month","Sex"),colnames(Table))], "X")
   return(Table)
   
 }
 
  if(exists('FL.TDGDFL'))   
  {
    #keep observation within logical size range and a Min.obs
    for(n in 1:length(FL.TDGDFL)) FL.TDGDFL[[n]]=FL.TDGDFL[[n]]%>%filter(TL<=Max.TL)
    
    TL_observers=vector('list',length(FL.TDGDFL))
    names(TL_observers)=names(FL.TDGDFL)
    dummy=do.call(rbind,FL.TDGDFL)
    MIN=min(dummy$TL/10,na.rm=T)
    MAX=max(dummy$TL/10,na.rm=T)
    for(n in 1:length(FL.TDGDFL))
    {
      if(nrow(FL.TDGDFL[[n]])>Min.obs)
      {
        a=fn.add.missing.size3(dat=FL.TDGDFL[[n]],
                               Min=floor(MIN)*10,
                               Max=ceiling(MAX)*10,
                               interval=Bin.size)
        ag=a[,-(2:3)]%>%
            group_by(FINYEAR) %>% 
            summarise_each(sum)
        
        ag=ag[which(rowSums(ag[,-1])>Min.obs),]
        a=a%>%filter(FINYEAR%in%ag$FINYEAR)
        TL_observers[[n]]=a
      }
    }
    TL_observers <- TL_observers[!sapply(TL_observers,is.null)]
    TL_observers <- TL_observers[sapply(TL_observers,nrow)>0]
  }
  if(exists('FL_Pilbara_trawl'))
  {
    FL_Pilbara_trawl=FL_Pilbara_trawl%>%filter(TL<=Max.TL)
    if(nrow(FL_Pilbara_trawl)>=Min.obs)
    {
      Pil.trwl_observers=fn.add.missing.size3(dat=FL_Pilbara_trawl,
                                              Min=floor(min(FL_Pilbara_trawl$TL/10,na.rm=T))*10,
                                              Max=ceiling(max(FL_Pilbara_trawl$TL/10,na.rm=T))*10,
                                              interval=Bin.size)
    }else rm(FL_Pilbara_trawl)

  }
  if(exists('FL_NSF'))
  {
    FL_NSF=FL_NSF%>%filter(TL<=Max.TL)
    if(nrow(FL_NSF)>=Min.obs)
    {
      NSF_observers=fn.add.missing.size3(dat=FL_NSF,
                                         Min=floor(min(FL_NSF$TL/10,na.rm=T))*10,
                                         Max=ceiling(max(FL_NSF$TL/10,na.rm=T))*10,
                                         interval=Bin.size)
    }else rm(FL_NSF)

  }
    
  #Wetline_WRL
  # Size.comp.Dusky.WRL=data.frame(FINYEAR="1999-00",TL=Dusky.WRL$TL)
  # Size.comp.Dusky.WRL$FINYEAR=as.character(Size.comp.Dusky.WRL$FINYEAR)
  # Size.comp.Dusky.WRL$FL=round(TL_to_FL(Dusky.WRL$TL,a.TL,b.TL))    
  # Size.comp.Dusky.WRL$Number=1
  # Size.comp.Dusky.WRL$Year=1
  # Size.comp.Dusky.WRL$Type="Fisher_measure"
  
    #2.3. Put all size composition (as TL) data in list  
  if(exists('FL.TDGDFL'))
  {
    drop.obs=function(dd,ZN) #drop non-representative observations in TDGDLF (other fisheries are ok) 
    {
      dd$dummy=paste(dd$FINYEAR,ZN)
      dd=subset(dd,dummy%in%This.yr.zn)
      dd=dd[,-match('dummy',names(dd))]
      return(dd)
    }
    
    All.size=TL_observers[grep("6.5",names(TL_observers))]
    if(length(All.size)>0)
    {
      names(All.size)=str_remove(names(All.size), ".6.5")
      for(i in 1:length(All.size)) All.size[[i]]=drop.obs(All.size[[i]],names(All.size)[i])
      names(All.size)=ifelse(names(All.size)=='West','WC',
                      ifelse(names(All.size)=='Zone1','Zn1',
                      ifelse(names(All.size)=='Zone2','Zn2',
                      NA)))
      All.size=list(TDGDLF=All.size) 
    }
    if(length(All.size)==0) rm(All.size)
    
    All.size_7=TL_observers[grep("7",names(TL_observers))]
    if(length(All.size_7)>0)
    {
      names(All.size_7)=str_remove(names(All.size_7), ".7")
      for(i in 1:length(All.size_7)) All.size_7[[i]]=drop.obs(All.size_7[[i]],names(All.size_7)[i])
      names(All.size_7)=ifelse(names(All.size_7)=='West','WC',
                        ifelse(names(All.size_7)=='Zone1','Zn1',
                        ifelse(names(All.size_7)=='Zone2','Zn2',
                        NA))) 
      All.size_7=All.size_7[sapply(All.size_7, nrow)>0]
      if(length(All.size_7)>0)
      {
        if(exists('All.size')) All.size$TDGDLF_7=All.size_7
        if(!exists('All.size')) All.size=list(TDGDLF_7=All.size_7)
      }
        
    }

  }
  
  if(exists('FL_Pilbara_trawl'))
  {
    if(exists('All.size')) All.size$Pilbara_trawl=Pil.trwl_observers
    if(!exists('All.size')) All.size=list(Pilbara_trawl=Pil.trwl_observers)
  }
  if(exists('FL_NSF'))
  {
    if(exists('All.size')) All.size$NSF=NSF_observers
    if(!exists('All.size')) All.size=list(NSF=NSF_observers)
  }
  if(exists('All.size'))
  {
    All.size=All.size[sapply(All.size, length)>0]
    if(length(All.size)==0) rm(All.size)
  }
  if(exists('All.size_7')) rm(All.size_7)

  
  #3. Put abundance data together 
  if(exists('Ab.index.Srvy.FixSt')) Naturaliste.abun=Ab.index.Srvy.FixSt
  

  #4. set working directory for outputing figures
  HandL="C:/Matias/Analyses/Population dynamics/1."
  DiR=paste(HandL,Name,"/",Yr.assess,"/1_Inputs/Visualise data",sep='')
  if(!file.exists(DiR))
  {
    mainDir=paste(HandL,Name,"/",Yr.assess,sep="")
    dir.create(mainDir)
    subDir="1_Inputs"
    dir.create(file.path(mainDir,subDir))
    subDir="/1_Inputs/Visualise data"
    dir.create(file.path(mainDir,subDir))
  }
  setwd(DiR)
  
  
  #5. Years with data for conventional tagging data
  if(exists("Zn.rel_Conv.Tag"))
  {
    a=Zn.rel_Conv.Tag
    a$Finyr=with(a,ifelse(Mn.rel>6,Yr.rel+1,Yr.rel))
    c.Tag=list(TDGDLF=data.frame(FINYEAR=paste(sort(unique(a$Finyr)),"-",fn.subs(sort(unique(a$Finyr))+1),sep="")))  
    
    #Add financial year to release and recapture files
    fn.finyr.rel=function(dat)with(dat,ifelse(Mn.rel>6,Yr.rel,Yr.rel-1))
    fn.finyr.rec=function(dat)with(dat,ifelse(Mn.rec>6,Yr.rec,Yr.rec-1))
    fn.finyr=function(Finyr)paste(Finyr,"-",fn.subs(Finyr+1),sep="")
    
    Zn.rel_Conv.Tag$FinYear.rel=fn.finyr.rel(Zn.rel_Conv.Tag)
    Zn.rec_Conv.Tag$FinYear.rec=fn.finyr.rec(Zn.rec_Conv.Tag)
    Zn.rel_Conv.Tag$FinYear.rel=fn.finyr(Zn.rel_Conv.Tag$FinYear.rel)
    Zn.rec_Conv.Tag$FinYear.rec=fn.finyr(Zn.rec_Conv.Tag$FinYear.rec)
    
    #Recode TG.zn
    fn.recode=function(dat,dat1)
    {
      dat=aggregate(Number~TG.zn+Rel.zone+Yr.rel+Mn.rel+Age+FinYear.rel,dat,sum)
      dat1=aggregate(Number~TG.zn+Rec.zone+Yr.rec+Mn.rec+FinYear.rec,dat1,sum)
      
      dat=dat[order(dat$FinYear.rel,dat$Mn.rel),]
      dat$TG=1:nrow(dat)
      dat1=merge(dat1,subset(dat,select=c(TG.zn,TG)),by="TG.zn",all.x=T)
      return(list(dat=dat,dat1=dat1))
    }
    a=fn.recode(Zn.rel_Conv.Tag,Zn.rec_Conv.Tag)
    Zn.rel_Conv.Tag=a$dat
    Zn.rec_Conv.Tag=a$dat1
    
    #Create table of proportion of recaptures by year
    fn.plot.prop.rec.yr=function(dat)
    {
      a=aggregate(Number~FinYear.rec,dat,sum)
      a$Prop=a$Number/sum(a$Number)
      Yr1=as.numeric(substr(min(dat$FinYear.rec),start=1,stop=4))
      Yr1=1993
      Yr2=as.numeric(substr(max(dat$FinYear.rec),start=1,stop=4))
      SEQ=seq(Yr1,Yr2)
      SEQ1=seq(Yr1+1,Yr2+1)
      All.yrs=paste(SEQ,"-",substr(SEQ1,3,4),sep="")
      id=All.yrs[which(!All.yrs%in%a$FinYear.rec)]
      Add=data.frame(FinYear.rec=id,Number=0,Prop=0)
      a=rbind(a,Add)
      a=a[order(a$FinYear.rec),]
      plot(1:length(a$FinYear.rec),a$Prop,xaxt='n',ylab="",xlab="",cex.axis=1.25,
           pch=19,col=2,cex=2)
      axis(1,1:length(a$FinYear.rec),F,tck=-0.03)
      axis(1,seq(1,length(a$FinYear.rec),2),a$FinYear.rec[seq(1,length(a$FinYear.rec),2)],cex.axis=.9,tck=-0.06)
    }
    fn.fig("Proportion.tag.recaptures",2400, 1800)
    par(mfcol=c(1,1),las=1,mai=c(0.3,0.55,.1,.1),oma=c(2.25,2.25,.1,.1),mgp=c(1,.5,0))
    fn.plot.prop.rec.yr(Zn.rec_Conv.Tag)
    legend("topright",Name,bty='n',cex=1.5)
    mtext("Financial year",1,line=1,cex=1.5,outer=T)
    mtext("Proportion of recaptures",2,line=1,cex=1.5,outer=T,las=3)
    dev.off()
  }
  
  
  #6. Years with data for acoustic tagging data
  if(exists("Zn.rel_Acous.Tag"))
  {
    #Add financial year
    fn.finyr.rel=function(dat)with(dat,ifelse(Month.rel>6,Year.rel,Year.rel-1))
    fn.finyr.rec=function(dat)with(dat,ifelse(Month>6,Year,Year-1))
    
    Zn.rel_Acous.Tag$FinYear.rel=fn.finyr.rel(Zn.rel_Acous.Tag)
    Zn.rec_Acous.Tag$FinYear.rec=fn.finyr.rec(Zn.rec_Acous.Tag)
    Zn.rel_Acous.Tag$FinYear.rel=fn.finyr(Zn.rel_Acous.Tag$FinYear.rel)
    Zn.rec_Acous.Tag$FinYear.rec=fn.finyr(Zn.rec_Acous.Tag$FinYear.rec)
    
    a.Tag=sort(unique(Zn.rec_Acous.Tag$FinYear.rec))
    
  }
  
  #7. Add missing years to prop.eff.mesh
  #note: back fill using temporal change in 7inch upto 2009-10
  if('TDGDLF'%in%names(catch))
  {
    dum.yr=sort(as.character(unique(catch$TDGDLF%>%filter(Data.set=='Data.monthly')%>%pull(FINYEAR))))
    id=dum.yr[which(!dum.yr%in%Mesh.prop.eff$finyear)]
    add.msh.eff=Mesh.prop.eff[(nrow(Mesh.prop.eff)+1):((nrow(Mesh.prop.eff))+length(id)),]
    add.msh.eff$finyear=id
    
    Mesh.prop.eff=rbind(add.msh.eff,Mesh.prop.eff)
    Mesh.prop.eff.West=rbind(add.msh.eff,Mesh.prop.eff.West)
    Mesh.prop.eff.Zn1=rbind(add.msh.eff,Mesh.prop.eff.Zn1)
    Mesh.prop.eff.Zn2=rbind(add.msh.eff,Mesh.prop.eff.Zn2)
    
    bck.fil.yrs=match(c("2005-06","2006-07","2007-08",
                        "2008-09","2009-10"),Mesh.prop.eff$finyear)
    fn.bck.fil=function(a,intercpt)
    {
      a$yr=1:nrow(a)
      a$OFFSET=intercpt
      DaTt=a[bck.fil.yrs,]
      #set intercept 
      Model=lm(X165~0+yr,offset=OFFSET,DaTt)
      NEWDAT=data.frame(yr=a[1:(bck.fil.yrs[1]-1),]$yr,
                        OFFSET=unique(a$OFFSET))
      P_6.5=predict(Model,newdata = NEWDAT)
      P_7=1-P_6.5
      return(data.frame(X165=P_6.5,X178=P_7))
    }
    Mesh.prop.eff[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff,intercpt=0)
    Mesh.prop.eff.West[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff.West,intercpt=0)
    Mesh.prop.eff.Zn1[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff.Zn1,intercpt=1)
    Mesh.prop.eff.Zn2[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff.Zn2,intercpt=0)
    colnames(Mesh.prop.eff)[2:3]=colnames(Mesh.prop.eff.West)[2:3]=
      colnames(Mesh.prop.eff.Zn1)[2:3]=colnames(Mesh.prop.eff.Zn2)[2:3]=c("Mesh_6.5","Mesh_7")
    
  }
  
  

  #--------------- RESULTS SECTION --------------

  #1. Visualize catch
  
  #Total color
  tot.col="tan3"
  
  #Zone colors
  Zns=c("Joint","North","Closed","West","Closed.metro","Zone1","Zone2")
  Zns.leg=c("JANSF","WANCSF","Closed","WCDGDLF","Metro closure","Zone1","Zone2")
  COL.prop=c("aquamarine2","lightseagreen",'forestgreen',"lightgreen","olivedrab4","olivedrab3","mediumseagreen")
  names(COL.prop)=names(Zns.leg)=Zns
  
  #Biological regions colors
  Bio.col=c("dodgerblue","darkorchid4","cyan4","lightpink3")
  names(Bio.col)=c("North Coast","Gascoyne","West Coast","South Coast")
 
  fn.add=function(a)
  {
    id=which(!YRS%in%a$FINYEAR)
    if(length(id)>0)
    {
      ADD=a
      ADD=ADD[1:length(id),]
      ADD$FINYEAR=YRS[id]
      xx=match("FINYEAR",colnames(ADD))
      ADD[,which(!1:ncol(ADD)%in%xx)]=NA
      a=rbind(a,ADD)
      a=a[order(a$FINYEAR),]
    }
    return(a)
  }
  fn.see.Ktch=function(DAT,DAT.zone,YRS,LWD,TCK)
  {
    for(i in 1:length(DAT))
    {
      a=DAT[[i]]
      a.zn=DAT.zone[[i]]%>%filter(!is.na(zone))
      a=a%>%
        group_by(FINYEAR)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
        arrange(FINYEAR)
      Ylim=c(0,max(a$LIVEWT.c))
      if(nrow(a)==1)
      {
        plot(1:length(a$FINYEAR),a$LIVEWT.c,xaxt='n',ylab="",xlab="",ylim=Ylim,
             main=names(DAT)[i],col=tot.col)
      }else
      {
        a=fn.add(a)
        plot(1:length(a$FINYEAR),a$LIVEWT.c,lwd=LWD*.75,xaxt='n',ylab="",xlab="",
             ylim=Ylim,type='l',main=names(DAT)[i],col=tot.col)
        if(nrow(a.zn)>0)
        {
          a.zn=aggregate(LIVEWT.c~FINYEAR+zone,a.zn,sum)
          dis.z=unique(a.zn$zone)
          for(z in 1:length(dis.z))
          {
            b=subset(a.zn,zone==dis.z[z])
            b=fn.add(b)
            lines(1:length(b$FINYEAR),b$LIVEWT.c,lwd=LWD,col=COL.prop[match(dis.z[z],names(COL.prop))])
          }
          lines(1:length(a$FINYEAR),a$LIVEWT.c,lwd=LWD*.75,col=tot.col)
        }
        LGn=c('Total')
        CLS=tot.col
        if(nrow(a.zn)>0)
        {
          Cl.z=COL.prop[match(dis.z,names(COL.prop))]
          LGn=c(LGn,Zns.leg[match(dis.z,names(Zns.leg))])
          CLS=c(CLS,Cl.z)
        }
        legend('topleft',LGn,lty=1,lwd=2,col=CLS,bty='n')
      }
      axis(1,1:length(a$FINYEAR),F,tck=TCK)
      axis(1,seq(1,length(a$FINYEAR),5),a$FINYEAR[seq(1,length(a$FINYEAR),5)],tck=2*TCK)
    }
  }
  a=do.call(rbind,catch)
  YR.span=sort(as.character(unique(a$FINYEAR)))
  rm(a)
  YRS=YR.span

  
  n.plots=length(catch)
  WIDTH=ifelse(n.plots>5,2400,2400)
  LENGTH=ifelse(n.plots>5,2000,2400)
  fn.fig("All catches",WIDTH, LENGTH)
  par(las=1)
  smart.par(n.plots=n.plots,MAR=c(1.75,1.75,1.75,1.75),OMA=c(2,2,.1,.1),MGP=c(1.5,.6,0))
  fn.see.Ktch(DAT=catch,
              DAT.zone=catch.zone,
              YRS=YR.span,
              LWD=2.5,
              TCK=-.02)
  mtext("Financial year",1,line=.75,cex=1.25,outer=T)
  mtext("Total catch (tonnes)",2,line=.4,cex=1.25,outer=T,las=3)
  dev.off()
  

  #2. Visualize size composition                             
  if(exists('All.size'))
  {
    fnx=function(x) as.numeric(substr(x,1,4))
    fn.bub=function(DAT,COL,Scale,TCK)
    {
      if(nrow(DAT)==0)
      {
        plot.new()
        legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
        return(1)
      }
         else
      {
        #aggregate data into years
        YRs=unique(DAT$FINYEAR)
        this=which(!colnames(DAT)%in%c('FINYEAR','Month','Sex'))
        SIzes=colnames(DAT[,this])
        d=matrix(ncol=length(SIzes),nrow=length(YRs))
        for(s in 1:length(YRs))
        {
          a=subset(DAT,FINYEAR==YRs[s])
          d[s,]=colSums(a[,this])
        }
        colnames(d)=SIzes
        d=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
        
        #keep years with minimum observations
        drop.yr=rowSums(d[,2:ncol(d)])
        names(drop.yr)=d$FINYEAR
        drop.yr=subset(drop.yr,drop.yr<Min.obs)
        d=subset(d,!FINYEAR%in%names(drop.yr))
        
        if(nrow(d)==0)
        {
          plot.new()
          legend('center',paste("Less than",Min.obs,"observations"),
                 cex=1.25,text.col='steelblue',bty='n')
          return(1)
        }else
        {
          DAT=d
          Ylabs=FinYrs
          ID=which(!FinYrs%in%DAT$FINYEAR)
          ADD=data.frame(FINYEAR=FinYrs[ID])
          if(nrow(ADD)>0)
          {
            ADD1=DAT[1:nrow(ADD),2:ncol(DAT)]
            ADD=cbind(ADD,ADD1)
            ADD[,2:ncol(ADD)]=0
            DAT=rbind(DAT,ADD)
          }
          DAT$FINYEAR=as.character(DAT$FINYEAR)
          DAT=DAT[order(DAT$FINYEAR),]
          DAT=DAT[,2:ncol(DAT)]
          x=as.numeric(colnames(DAT))
          y=1:nrow(DAT)
          z=(as.matrix(DAT))      
          n=length(y)      
          xo=outer(x,rep(1,length=length(y)))
          yo=t(outer(y,rep(1,length=length(x))))
          zo=z/ rowSums(z)*Scale   
          N=max(10,n/2)
          matplot(yo,xo,type="n",xlab="",ylab="",xaxt='n',yaxt='n')
          #abline(v=pretty(x),lty=3,col="black")
          for(i in 1:n) points(yo[,i],xo[,i],cex=zo[i,],pch=16,col=COL)
          axis(1,y,F,tck=TCK)
          axis(2,x,F,tck=TCK)
          
          if(length(Ylabs)>10)
          {
             Labx=y[seq(1,length(Ylabs),5)]
             Labx.lgn=Ylabs[seq(1,length(Ylabs),5)]
          }else
          {
            Labx=y[seq(1,length(Ylabs),1)]
            Labx.lgn=Ylabs[seq(1,length(Ylabs),1)]
          }

            
          axis(1,Labx,Labx.lgn,tck=TCK*2,cex.axis=1.25)

          axis(2,seq(x[1],x[length(x)],by=N),seq(x[1],x[length(x)],by=N),tck=TCK*2,cex.axis=1.25)
          legend('topright',paste('n=',sum(z)),bty='n')
          return(0)
        }
      }
    }
  
      #2.1 TDGDLF size comp of reported catch
    if(sum(c('TDGDLF','TDGDLF_7')%in%names(All.size))>0)
    {
      a=NULL
      b=NULL
      if('TDGDLF'%in%names(All.size)) a=do.call(rbind,All.size$TDGDLF)
      if('TDGDLF_7'%in%names(All.size)) b=do.call(rbind,All.size$TDGDLF_7)
      FinYrs=rbind(a,b)
      FinYrs=sort(as.character(unique(FinYrs%>%pull(FINYEAR))))
      FinYrs=paste(fnx(min(FinYrs)):(fnx(max(FinYrs))),"-",
                   substr((fnx(min(FinYrs))+1):(fnx(max(FinYrs))+1),3,4),sep="")
      
      if(sum(c('TDGDLF','TDGDLF_7')%in%names(All.size))==2)
      {
        fn.fig("Size.comp.TDGDLF",2400, 2000)
        par(mfcol=c(3,2),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.25,1,1),mgp=c(1,.85,0))
        #-6.5
        DROP_6.5=vector('list',length(All.size$TDGDLF))
        for(e in 1:length(All.size$TDGDLF))
        {
          DROP_6.5[[e]]=fn.bub(DAT=All.size$TDGDLF[[e]],COL='steelblue',Scale=7.5,TCK=-.015)
          if(e==1) mtext('6.5 inch mesh',3)
        }
        if(e<3)
        {
          misin=3-length(All.size$TDGDLF)
          for(m in 1:misin)
          {
            plot.new()
            legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
          }
        }
        
        #-7
        DROP_7=vector('list',length(All.size$TDGDLF_7))
        for(e in 1:length(All.size$TDGDLF_7))
        {
          DROP_7[[e]]=fn.bub(DAT=All.size$TDGDLF_7[[e]],COL='steelblue',Scale=7.5,TCK=-.015)
          if(e==1) mtext('7 inch mesh',3)
          mtext(names(All.size$TDGDLF_7)[e],4,las=3,line=.35)
        }
        if(e<3)
        {
          misin=3-length(All.size$TDGDLF_7)
          for(m in 1:misin)
          {
            plot.new()
            legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
          }
        }
        mtext("Financial Year",1,cex=1.35,line=0.75,outer=T)
        mtext("Total length class (cm)",2,cex=1.35,line=0.5,las=3,outer=T)
        dev.off()
        
        if(sum(c(unlist(DROP_6.5),unlist(DROP_7)))==length(All.size$TDGDLF)+length(All.size$TDGDLF_7))
        {
          rm(All.size,TDGDFL.size.numbers)
        }
        
      }
      if('TDGDLF'%in%names(All.size)& ! 'TDGDLF_7'%in%names(All.size))
      {
        fn.fig("Size.comp.TDGDLF",2400, 2000)
        par(mfcol=c(3,1),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.25,1,1),mgp=c(1,.85,0))
        #-6.5
        DROP_6.5=vector('list',length(All.size$TDGDLF))
        for(e in 1:length(All.size$TDGDLF))
        {
          DROP_6.5[[e]]=fn.bub(DAT=All.size$TDGDLF[[e]],COL='steelblue',Scale=7.5,TCK=-.015)
          if(e==1) mtext('6.5 inch mesh',3)
        }
        if(e<3)
        {
          misin=3-length(All.size$TDGDLF)
          for(m in 1:length(misin))
          {
            plot.new()
            legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
          }
        }
         mtext("Financial Year",1,cex=1.35,line=0.75,outer=T)
        mtext("Total length class (cm)",2,cex=1.35,line=0.5,las=3,outer=T)
        dev.off()
        
        if(sum(unlist(DROP_6.5))==length(All.size$TDGDLF))
        {
          rm(All.size,TDGDFL.size.numbers)
        }
        
      }
      if(!'TDGDLF'%in%names(All.size)&  'TDGDLF_7'%in%names(All.size))
      {
        fn.fig("Size.comp.TDGDLF",2400, 2000)
        par(mfcol=c(3,1),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.25,1,1),mgp=c(1,.85,0))
        #-7
        DROP_7=vector('list',length(All.size$TDGDLF_7))
        for(e in 1:length(All.size$TDGDLF_7))
        {
          DROP_7[[e]]=fn.bub(DAT=All.size$TDGDLF_7[[e]],COL='steelblue',Scale=7.5,TCK=-.015)
          if(e==1) mtext('7 inch mesh',3)
          mtext(names(All.size$TDGDLF_7)[e],4,las=3,line=.35)
        }
        if(e<3)
        {
          misin=3-length(All.size$TDGDLF)
          for(m in 1:length(misin))
          {
            plot.new()
            legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
          }
        }
        mtext("Financial Year",1,cex=1.35,line=0.75,outer=T)
        mtext("Total length class (cm)",2,cex=1.35,line=0.5,las=3,outer=T)
        mtext("Financial Year",1,cex=1.35,line=0.75,outer=T)
        mtext("Total length class (cm)",2,cex=1.35,line=0.5,las=3,outer=T)
        dev.off()
        
        if(sum(unlist(DROP_7))==length(All.size$TDGDLF_7))
        {
          rm(All.size,TDGDFL.size.numbers)
        }
      }
    }
    
    
      #2.2 NSF size comp of reported catch
    if('NSF'%in%names(All.size))
    {
      FinYrs=sort(as.character(unique(All.size$NSF$FINYEAR)))
      FinYrs=paste(fnx(min(FinYrs)):(fnx(max(FinYrs))),"-",
                   substr((fnx(min(FinYrs))+1):(fnx(max(FinYrs))+1),3,4),sep="")
      fn.fig("Size.comp.NSF",2400, 2000)
      par(mfcol=c(1,1),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.5,1,1),mgp=c(1,.65,0))
      fn.bub(DAT=All.size$NSF,COL='steelblue',Scale=7.5,TCK=-.01)
      mtext("Financial Year",1,cex=1.75,line=0.75,outer=T)
      mtext("Total length class (cm)",2,cex=1.75,line=0.95,las=3,outer=T)
      dev.off()
    }
    
      #2.3 Pilbara trawl size comp of reported catch
    if('Pilbara_trawl'%in%names(All.size))
    {
      FinYrs=sort(as.character(unique(All.size$Pilbara_trawl$FINYEAR)))
      FinYrs=paste(fnx(min(FinYrs)):(fnx(max(FinYrs))),"-",
                   substr((fnx(min(FinYrs))+1):(fnx(max(FinYrs))+1),3,4),sep="")
      fn.fig("Size.comp.Pilbara.trawl",2400, 2000)
      par(mfcol=c(1,1),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.5,1,1),mgp=c(1,.65,0))
      fn.bub(DAT=All.size$Pilbara_trawl,COL='steelblue',Scale=7.5,TCK=-.01)
      mtext("Financial Year",1,cex=1.75,line=0.75,outer=T)
      mtext("Total length class (cm)",2,cex=1.75,line=0.95,las=3,outer=T)
      dev.off()
    }
  }

    #Table of number of observations and shots
  if(exists('TDGDFL.size.numbers'))
  {
    Size.numbers=TDGDFL.size.numbers%>%dplyr::select(-Species)
    
    fn.table.shots=function(dat,FSHRY)
    {
      a=dat
      a$Number=1
      Obs=aggregate(Number~FINYEAR,a,sum)
      a$Dup=paste(a$year,a$Month,a$SHEET_NO)
      bb=a[!duplicated(a$Dup),]
      bb$Number=1
      Shots=aggregate(Number~FINYEAR,bb,sum)
      this=merge(Obs,Shots,by="FINYEAR")
      names(this)[2:3]=c("N.observations","N.shots")
      this$Species=unique(a$SPECIES)
      this$Fishery=FSHRY
      this$zone=FSHRY
      
      return(this)
    }

    if(exists('FL_NSF'))
    {
      NSF.size.numbers=fn.table.shots(dat=FL_NSF,FSHRY="NSF")
      Size.numbers=rbind(Size.numbers,NSF.size.numbers)
    }
    if(exists('FL_Pilbara_trawl'))
    {
      Pilbara_trawl.size.numbers=fn.table.shots(FL_Pilbara_trawl,FSHRY="Pilbara trawl")
      Size.numbers=rbind(Size.numbers,Pilbara_trawl.size.numbers)
    }
    
    Size.numbers=Size.numbers[,match(c("Fishery","zone","FINYEAR",
                                       "N.observations","N.shots"),names(Size.numbers))]
    Size.numbers=Size.numbers[order(Size.numbers$FINYEAR,Size.numbers$Fishery,Size.numbers$zone),]
    Numbers.SF=reshape(Size.numbers[,-match("N.shots",names(Size.numbers))],v.names = "N.observations",
                       idvar = c("Fishery","zone"),timevar = "FINYEAR", direction = "wide")
    X=colnames(Numbers.SF)[3:ncol(Numbers.SF)]
    colnames(Numbers.SF)[3:ncol(Numbers.SF)]=as.character(sapply(strsplit(X,"N.observations."), "[", 2))
    Shots.SF=reshape(Size.numbers[,-match("N.observations",names(Size.numbers))],v.names = "N.shots",
                     idvar = c("Fishery","zone"),timevar = "FINYEAR", direction = "wide")
    X=colnames(Shots.SF)[3:ncol(Shots.SF)]
    colnames(Shots.SF)[3:ncol(Shots.SF)]=as.character(sapply(strsplit(X,"N.shots."), "[", 2))
    Numbers.SF=Numbers.SF[order(Numbers.SF$Fishery,Numbers.SF$zone),]
    Shots.SF=Shots.SF[order(Shots.SF$Fishery,Shots.SF$zone),]
    
    #add total
    dum= Numbers.SF[1,]
    dum$zone="Total"
    dum[,-match(c("Fishery","zone"),names(Numbers.SF))]=colSums(Numbers.SF[,-match(c("Fishery","zone"),names(Numbers.SF))],na.rm =T)
    Numbers.SF=rbind(Numbers.SF,dum)
    dum[,-match(c("Fishery","zone"),names(Numbers.SF))]=colSums(Shots.SF[,-match(c("Fishery","zone"),names(Shots.SF))],na.rm =T)
    Shots.SF=rbind(Shots.SF,dum)
    
    #remove NAs
    Numbers.SF[is.na(Numbers.SF)]=""
    Shots.SF[is.na(Shots.SF)]=""
    
    #create nice table 
    fn.word.table(WD=getwd(),TBL=Numbers.SF,Doc.nm="Size.comp.n.observations",caption=NA,paragph=NA,
                  HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                  Zebra='NO',Zebra.col='grey60',Grid.col='black',
                  Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
    fn.word.table(WD=getwd(),TBL=Shots.SF,Doc.nm="Size.comp.n.shots",caption=NA,paragph=NA,
                  HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                  Zebra='NO',Zebra.col='grey60',Grid.col='black',
                  Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
    write.csv(Numbers.SF,"Numbers.SF.csv",row.names=F)
    write.csv(Shots.SF,"Shots.SF.csv",row.names=F)
  }

 
  #3. Visualize mean weights
 if(exists("Avr.wt.yr.zn"))
 {
   fn.see.avg.wgt=function()
   {
     b=Avr.wt.yr.zn
     zn=unique(b$zone)
     N=1:length(unique(b$Finyear))
     SD=b$mean*(b$CV)
     plot(N,ylim=c(0,max(b$mean+SD)*1.05),main="",cex.main=1.25,xaxt='n',
          ylab="",xlab="",pch=19,cex=3,cex.axis=1,col="transparent",xlim=c(0,N[length(N)]+0.5))
     
     jit=c(0,.1,.2)
     CLOS=COL.prop[match(c("West","Zone1","Zone2"),names(COL.prop))]
     
     for(x in 1:length(zn))
     {
       a=subset(b, zone==zn[x])
       N1=N+jit[x]
       SD=a$mean*(a$CV)
       
       points(N1,a$mean,ylim=c(0,max(a$mean+SD)*1.05),main=zn[x],cex.main=1.75,xaxt='n',
              ylab="",pch=19,cex=2,cex.axis=1.25,col=CLOS[x])
       segments(N1,a$mean,N1,a$mean-SD,lwd=2,col=CLOS[x])
       segments(N1,a$mean,N1,a$mean+SD,lwd=2,col=CLOS[x])
       axis(1,N,F,tck=-0.015)    
     }
     axis(1,seq(1,length(N),2),a$Finyear[seq(1,length(N),2)],tck=-0.02,cex.axis=1.25)
     legend('bottomleft',zn,pch=19,col=CLOS,pt.cex=1.5,cex=1.2,bty='n')
   }
   fn.fig("Avg.wgt.zn",2000,2000)
   par(mfcol=c(1,1),las=1,mai=c(0.45,0.35,.1,.15),oma=c(1,2.25,.1,.1),mgp=c(1,.65,0))
   fn.see.avg.wgt()
   mtext("Relative live weight",2,line=.5,cex=1.5,las=3,outer=T)
   mtext("Financial Year",1,cex=1.5,line=0,outer=T)
   dev.off()
 }

  
  #4. Select effort years
  Eff.zn=subset(Eff.zn,FINYEAR%in%YR.span)
  Eff.total=data.frame(FINYEAR=Eff.zn$FINYEAR,Total=Eff.zn$West+Eff.zn$Zone1+Eff.zn$Zone2)
  
  
  #5. Visualize data availability
  visualize.dat=function(d.list)   
  {
    cum.plots=rep(NA,length(d.list))
    for(l in 1:length(d.list)) cum.plots[l]=length(d.list[[l]])
    cum.plots=cumsum(cum.plots)
    n.rows=cum.plots[length(cum.plots)]
    LABS=vector('list',length(d.list))
    plot(as.numeric(substr(YR.span,1,4)),as.numeric(substr(YR.span,1,4)),ylim=c(1,n.rows),yaxt='n',col='transparent',ylab='',xlab='')
    All.cols=matrix(hcl.colors(length(d.list)*2,palette = 'Geyser'),ncol=2,byrow=T)
    for(l in 1:length(d.list))
    {
      dd=d.list[[l]]
      LABS[[l]]=names(d.list)[l]
      dummy=rep(NA,length(dd))
      colfunc <- function(n,col1,col2) colorRampPalette(c(col1, col2))(n)
      COL=colfunc(length(dd),col1=All.cols[l,1],col2=All.cols[l,2])
      for(ll in 1:length(dd))
      {
        if(length(dd)>1)dummy[ll]=paste(LABS[[l]]," (",names(dd)[ll],")",sep='')
        if(is.character(dd[[ll]]))x=dd[[ll]] else
        {
          colnames(dd[[ll]])=tolower(colnames(dd[[ll]]))
          x=dd[[ll]]$finyear
        }
        x=as.numeric(substr(unique(x),1,4))
        if(l==1) y=ll
        if(l>1) y=cum.plots[l-1]+ll
        points(x,rep(y,length(x)),cex=6,pch="-",col=COL[ll])
        rm(x)
      }
      if(length(dd)>1) LABS[[l]]=dummy
    }
    axis(2,1:n.rows,unlist(LABS),las=2,cex.axis=1)
    mtext("Financial year",1,cex=1.25,line=1.75)
  }
  cpue_tdgdlf.monthly=NULL
  cpue_tdgdlf.daily=NULL
  if(exists('Ab.indx.TDGDLF.all')) cpue_tdgdlf.monthly=Ab.indx.TDGDLF.all
  if(exists('Ab.indx.TDGDLF.all.daily')) cpue_tdgdlf.daily=Ab.indx.TDGDLF.all.daily
  if(!is.null(cpue_tdgdlf.monthly)|!is.null(cpue_tdgdlf.daily)) Abun=list(TDGLDF=rbind(cpue_tdgdlf.monthly,
                                                                                              cpue_tdgdlf.daily))
  if(exists('Ab.indx.NSF'))
  {
    if(exists('Abun'))Abun$NSF=Ab.indx.NSF
    if(!exists('Abun'))Abun=list(NSF=Ab.indx.NSF)
  }
    
  if(exists('Naturaliste.abun'))
  {
    if(exists('Abun')) Abun$Survey=Naturaliste.abun
    if(!exists('Abun')) Abun=list(Survey=Naturaliste.abun)
  }
  
  d.list=list(Catch=catch) 
  if(exists('All.size'))
  {
    if(sum(c('TDGDLF','TDGDLF_7')%in%names(All.size))>0)
    {
      d.list$'Catch size comp (TDGDLF)'=list(FINYEAR=unique(do.call(rbind,All.size$TDGDLF)$FINYEAR))
    }else
    {
      d.list$'Catch size comp (NSF)'=list(FINYEAR=unique(All.size$NSF$FINYEAR))
    }
  }
  if(exists('Avr.wt.yr'))d.list$'Catch mean wt (TDGDLF)'=list(Avr.wt.yr%>%distinct(Finyear))
  if(exists('Abun'))
  {
    d.list$Abundance=Abun
    if(length(d.list$Abundance)==1) names(d.list)[match('Abundance',names(d.list))]=paste("Abundance (",names(d.list$Abundance),")",sep='')
  }
  if(exists('Size.index.Srvy.FixSt'))d.list$'Survey mean size'=list(FINYEAR=Size.index.Srvy.FixSt$FINYEAR)
  if(exists('c.Tag'))d.list$'Conv tag'=c.Tag
  if(exists('a.Tag'))d.list$'Acous tag'=list(FINYEAR=a.Tag)
  
  fn.fig("avail.dat",2400,1400)
  par(mar=c(3,2,.1,.1),oma=c(.1,9,.1,1),las=1,mgp=c(1.5,.6,0))
  visualize.dat(d.list)   
  dev.off()

  
  
  #6. Table of released indivuals with conventional tags
  if(exists('Zn.rel_Conv.Tag_size_adu'))
  {
    TBL.conv.rel=merge(Zn.rel_Conv.Tag_size_adu[,-match("TG.zn",names(Zn.rel_Conv.Tag_size_adu))],
                       Zn.rel_Conv.Tag_size_juv[,-match("TG.zn",names(Zn.rel_Conv.Tag_size_juv))],
                       by=c("Rel.zone","Yr.rel"),all=T)
    names(TBL.conv.rel)[match(c("Number.x","Number.y"),
                              names(TBL.conv.rel))]=c("Adult.N","Juvenile.N")
    TBL.conv.rel=TBL.conv.rel[,-match(c("size.group.x",
                                        "size.group.y"),names(TBL.conv.rel))]
    TBL.conv.rel[is.na(TBL.conv.rel)]=0
    write.csv(TBL.conv.rel,"TBL.conv.releases.csv",row.names=F)
    fn.word.table(WD=getwd(),TBL=TBL.conv.rel,Doc.nm="TBL.conv.releases",caption=NA,paragph=NA,
                  HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                  Zebra='NO',Zebra.col='grey60',Grid.col='black',
                  Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  }
  
  
  #7. Plot acoustic tag recaptures by year
  if(exists('Rep.Recap'))
  {
    Rep.Recap=subset(Rep.Recap,!is.na(year.rel))
    i.yrs=sort(unique(Rep.Recap$year.rec))
    Y.rng=range(c(Rep.Recap$RELLATDECDEG,Rep.Recap$RECLATDECDEG),na.rm=T)
    X.rng=range(c(Rep.Recap$RELLNGDECDEG,Rep.Recap$RECLNGDECDEG),na.rm=T)
    
    fn.fig("Acoustic.tags.rel.rec.map",2400,2400)
    par(mfcol=c(ceiling(length(i.yrs)/2),2),las=1,mai=c(1,1,.1,.1),oma=c(.1,.1,.1,.1))
    dummy=Rep.Recap
    for(i.yr in 1:length(i.yrs))
    {
      recs=subset(dummy,year.rec==i.yrs[i.yr])
      na.rcs=subset(recs,is.na(RECLATDECDEG))
      if(nrow(recs)>0)
      {
        rels=subset(dummy,year.rel<=i.yrs[i.yr])
        plot(rels$RELLNGDECDEG,rels$RELLATDECDEG,ylim=Y.rng,xlim=X.rng,ylab="",xlab="",pch=19)
        points(recs$RECLNGDECDEG,recs$RECLATDECDEG,col=2,pch=4,cex=1.75)
        if(nrow(na.rcs)==0) Rc.lg=paste(nrow(recs),"recaptured")
        if(nrow(na.rcs)>0) Rc.lg=paste(nrow(recs),"recaptured (",nrow(na.rcs),"without rec. info)")
        legend('top',c(paste(nrow(rels),"with tags"),Rc.lg),bty='n',col=c(1,2),
               text.col=c(1,2),pch=c(19,4),title=paste(i.yrs[i.yr]),cex=1.2)
        dummy=subset(dummy,!ATAG_NO%in%recs$ATAG_NO)
      }
      
    }
    mtext("Longitude",1,-2,outer=T,cex=1.75)
    mtext("Latitude",2,-2,las=3,outer=T,cex=1.75)
    dev.off()
  }
 

# Export data for population dynamics modelling ---------------------------
  HandL="C:/Matias/Data/Population dynamics/Data inputs for models/"
  DiR=paste(HandL,str_remove(capitalize(Name), ' shark'),"/",Yr.assess,sep='')
  if(!file.exists(DiR)) dir.create(DiR)
  
  ff=do.call(file.remove, list(list.files(DiR, full.names = TRUE))) #remove all files
  setwd(DiR)
  
  fn.agg.at.level.and.exprt=function(DAT,Level,VAR,SOURCE)
  {
   for(i in 1:length(DAT))
   {
     if(Level=='annual')
     {
       a=aggregate(LIVEWT.c~FINYEAR,DAT[[i]],sum)
       write.csv(a,paste(VAR,".",Level,".",SOURCE[i],".csv",sep=""),row.names=F)
     }
       
     if(Level=='annual.by.zone')
     {
        if(!is.null(DAT[[i]]$zone))
        {
          DAT[[i]]$zone=with(DAT[[i]],ifelse(is.na(zone),SOURCE[i],zone))
          a=aggregate(LIVEWT.c~FINYEAR+zone,DAT[[i]],sum)
          write.csv(a,paste(VAR,".",Level,".",SOURCE[i],".csv",sep=""),row.names=F)
        }
    }
    }
 }
 
  #Catch
  fn.agg.at.level.and.exprt(DAT=catch,
                            Level='annual',
                            VAR="ktch",
                            SOURCE=names(catch))
  fn.agg.at.level.and.exprt(DAT=catch.zone,
                            Level='annual.by.zone',
                            VAR="ktch",
                            SOURCE=names(catch.zone))   
 

  #Effort
  write.csv(Eff.zn,"effort.annual.by.zone.TDGDLF.csv",row.names=F) 
  write.csv(Eff.total,"effort.annual.TDGDLF.csv",row.names=F)

   
  #avg weight 
    #TDGDLF
  if(exists('Avr.wt.yr.zn'))write.csv(Avr.wt.yr.zn,"ktch.avg.weight.annual.by.zone.csv",row.names=F) 
  if(exists('Avr.wt.yr')) write.csv(Avr.wt.yr,"ktch.avg.weight.annual.csv",row.names=F)
    #Naturaliste survey 
  if(exists('Size.index.Srvy.FixSt')) write.csv(Size.index.Srvy.FixSt,"size.annual.survey.csv",row.names=F)
  
  
  #TDGLDF cpue
    #folly
  if(exists("Ab.folly.TDGDLF.all")) write.csv(Ab.folly.TDGDLF.all,"cpue.annual.TDGDLF.folly.csv",row.names=F)
    #by zone      
  if(exists('Ab.indx.TDGDLF')) write.csv(Ab.indx.TDGDLF,"cpue.annual.by.zone.TDGDLF.csv",row.names=F)
  if(exists('Ab.indx.TDGDLF.daily')) write.csv(Ab.indx.TDGDLF.daily,"cpue.annual.by.zone.TDGDLF.daily.csv",row.names=F)
    #zones combined
  if(exists('Ab.indx.TDGDLF.all')) write.csv(Ab.indx.TDGDLF.all,"cpue.annual.TDGDLF.csv",row.names=F)
  if(exists('Ab.indx.TDGDLF.all.daily')) write.csv(Ab.indx.TDGDLF.all.daily,"cpue.annual.TDGDLF.daily.csv",row.names=F)

  #Naturaliste survey     
  if(exists('Naturaliste.abun')) write.csv(Naturaliste.abun,"cpue.annual.survey.csv",row.names=F)

  #NSF cpue     
  if(exists('Ab.indx.NSF')) write.csv(Ab.indx.NSF,"cpue.annual.NSF.csv",row.names=F)
 
  #Size composition 
  if(exists('All.size'))
  {
    fn.exp.size=function(dat,Level)
    {
      n=length(dat)
      for( i in 1:n)
      {
        a=dat[[i]]
        
        if(class(a)=="data.frame") 
        {
          YRs=sort(as.character(unique(a$FINYEAR)))
          SIzes=colnames(a[,4:ncol(a)])
          d=matrix(ncol=length(SIzes),nrow=length(YRs))
          d.f=d.m=d
          for(s in 1:length(YRs))
          {
            ddat=subset(a,FINYEAR==YRs[s])
            ddat.f=subset(ddat,Sex=="F")
            ddat.m=subset(ddat,Sex=="M")
            
            d[s,]=colSums(ddat[,4:ncol(ddat)])
            d.f[s,]=colSums(ddat.f[,4:ncol(ddat.f)])
            d.m[s,]=colSums(ddat.m[,4:ncol(ddat.m)])
          }
          #sexes combined
          colnames(d)=SIzes
          a=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
          write.csv(a,paste("size.comp.",names(dat)[i],".csv",sep=""),row.names=F)
          #by sex
          colnames(d.f)=colnames(d.m)=SIzes
          a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.f))
          write.csv(a,paste("size.comp.fem.",names(dat)[i],".csv",sep=""),row.names=F)
          
          a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.m))
          write.csv(a,paste("size.comp.mal.",names(dat)[i],".csv",sep=""),row.names=F)
        }
        
        if(class(a)=="list" & length(a)>0)
        {
          if(Level=='annual')
          {
            a=do.call(rbind,a)
            YRs=sort(as.character(unique(a$FINYEAR)))
            SIzes=colnames(a[,4:ncol(a)])
            d=matrix(ncol=length(SIzes),nrow=length(YRs))
            d.f=d.m=d
            for(s in 1:length(YRs))
            {
              ddat=subset(a,FINYEAR==YRs[s])
              ddat.f=subset(ddat,Sex=="F")
              ddat.m=subset(ddat,Sex=="M")
              
              d[s,]=colSums(ddat[,4:ncol(ddat)])
              d.f[s,]=colSums(ddat.f[,4:ncol(ddat.f)])
              d.m[s,]=colSums(ddat.m[,4:ncol(ddat.m)])
            }
            #sexes combined
            colnames(d)=SIzes
            a=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
            write.csv(a,paste("size.comp.",Level,".",names(dat)[i],".csv",sep=""),row.names=F)
            
            #by sex
            colnames(d.f)=colnames(d.m)=SIzes
            a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.f))
            write.csv(a,paste("size.comp.fem.",Level,".",names(dat)[i],".csv",sep=""),row.names=F)
            
            a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.m))
            write.csv(a,paste("size.comp.mal.",Level,".",names(dat)[i],".csv",sep=""),row.names=F)
            
          }
          
          if(Level=='annual.by.zone') 
          {
            ss=names(a)
            for(p in 1:length(ss))
            {
              b=a[[p]]
              if(nrow(b)>0)
              {
                YRs=sort(as.character(unique(b$FINYEAR)))
                SIzes=colnames(b[,4:ncol(b)])
                d=matrix(ncol=length(SIzes),nrow=length(YRs))
                d.f=d.m=d
                for(s in 1:length(YRs))
                {      
                  ddat=subset(b,FINYEAR==YRs[s])
                  ddat.f=subset(ddat,Sex=="F")
                  ddat.m=subset(ddat,Sex=="M")
                  
                  d[s,]=colSums(ddat[,4:ncol(ddat)])
                  d.f[s,]=colSums(ddat.f[,4:ncol(ddat.f)])
                  d.m[s,]=colSums(ddat.m[,4:ncol(ddat.m)])
                  
                }
                
                #sexes combined
                colnames(d)=SIzes
                b=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
                write.csv(b,paste("size.comp.annual.",ss[p],".",names(dat)[i],".csv",sep=""),row.names=F)
                
                #by sex
                colnames(d.f)=colnames(d.m)=SIzes
                b=cbind(FINYEAR=YRs,as.data.frame.matrix(d.f))
                write.csv(b,paste("size.comp.annual.fem.",ss[p],".",names(dat)[i],".csv",sep=""),row.names=F)
                
                b=cbind(FINYEAR=YRs,as.data.frame.matrix(d.m))
                write.csv(b,paste("size.comp.annual.mal.",ss[p],".",names(dat)[i],".csv",sep=""),row.names=F)
              }
              
            }
          }
        }
      }
    }
    fn.exp.size(dat=All.size,Level='annual')
    fn.exp.size(dat=All.size,Level='annual.by.zone')
  }
  
  #Conventional tagging
    #1. Individual-based model
  if(exists('Rel_rec_Conv.Tag')) write.csv(Rel_rec_Conv.Tag,"Ind_based_model.csv",row.names=F)
    #2. At age
  #note: this function puts in SS3 format
  if(exists('Zn.rel_Conv.Tag'))
  {
    tagging.fn=function(dat,THESE,THOSE,type,txt)
    {
      dat=dat[,match(THESE,names(dat))]
      if(type=='rel')
      {
        dat$gender=0
        dat=dat[,match(c("TG","Rel.zone","FinYear.rel","Mn.rel","gender","Age","Number"),names(dat))]
      }
      
      if(type=='rec')dat=dat[,match(c("TG","FinYear.rec","Mn.rec","Rec.zone","Number"),names(dat))]
      names(dat)=THOSE
      write.csv(dat,paste(txt,"csv",sep=""),row.names=F)
    }
    
    #releases
    these=c("TG","Rel.zone","FinYear.rel","Mn.rel","Age","Number")
    those=c("TG","area","yr","season","gender","Age","Nrelease")  
    tagging.fn(Zn.rel_Conv.Tag,these,those,'rel','conv.tag.rel.')
    
    #recaptures
    these=c("TG","Rec.zone","FinYear.rec","Mn.rec","Number")
    those=c("TG","year","season","area","Number")
    tagging.fn(Zn.rec_Conv.Tag,these,those,'rec',"conv.tag.reca.")
  }
    #3. At size
      #all
  if(exists('Zn.rel_Conv.Tag_size')) write.csv(Zn.rel_Conv.Tag_size,"Zn.rel_Conv.Tag_size.csv",row.names=F)
  if(exists('Zn.rec_Conv.Tag_size')) write.csv(Zn.rec_Conv.Tag_size,"Zn.rec_Conv.Tag_size.csv",row.names=F)
      #juveniles and adults
  if(exists('Zn.rel_Conv.Tag_size_adu')) write.csv(Zn.rel_Conv.Tag_size_adu,"Zn.rel_Conv.Tag_size_adu.csv",row.names=F)
  if(exists('Zn.rel_Conv.Tag_size_juv')) write.csv(Zn.rel_Conv.Tag_size_juv,"Zn.rel_Conv.Tag_size_juv.csv",row.names=F)
  if(exists('Zn.rec_Conv.Tag_size_adu')) write.csv(Zn.rec_Conv.Tag_size_adu,"Zn.rec_Conv.Tag_size_adu.csv",row.names=F)
  if(exists('Zn.rec_Conv.Tag_size_juv')) write.csv(Zn.rec_Conv.Tag_size_juv,"Zn.rec_Conv.Tag_size_juv.csv",row.names=F)
  if(exists('Smallest_size_tagged')) write.csv(Smallest_size_tagged,"Smallest_size_tagged.csv",row.names=F)

  #Acoustic tagging 
    #releases
  if(exists('Zn.rel_Acous.Tag')) write.csv(Zn.rel_Acous.Tag,"acous.tag.rel.csv",row.names=F)  
  if(exists('Zn.rel_Acous.Tag.prop')) write.csv(Zn.rel_Acous.Tag.prop,"acous.tag.rel.prop.csv",row.names=F)
    #recaptures
  if(exists('Zn.rec_Acous.Tag')) write.csv(Zn.rec_Acous.Tag,"acous.tag.reca.csv",row.names=F)
  if(exists('Zn.rec_Acous.Tag.prop')) write.csv(Zn.rec_Acous.Tag.prop,"acous.tag.reca.prop.csv",row.names=F)
    #Individual-based model
  if(exists('Indiv_based_Acous.Tag')) write.csv(Indiv_based_Acous.Tag,"acous.tag.reca.indv.based.csv",row.names=F)
  #Recapture information from acoustic tags
  if(exists('Rep.Recap')) write.csv(Rep.Recap,"Acoustic.tag.rel.rec.csv",row.names=F) 
  
  #Age and growth 
  if(exists('Age.growth'))  write.csv(Age.growth,"Age.growth.csv",row.names=F)
  
  #Mesh proportional effort
  if(exists('Mesh.prop.eff')) write.csv(Mesh.prop.eff,"Mesh.prop.eff.csv",row.names=F)
  if(exists('Mesh.prop.eff.West')) write.csv(Mesh.prop.eff.West,"Mesh.prop.eff.West.csv",row.names=F)
  if(exists('Mesh.prop.eff.Zn1')) write.csv(Mesh.prop.eff.Zn1,"Mesh.prop.eff.Zn1.csv",row.names=F)
  if(exists('Mesh.prop.eff.Zn2')) write.csv(Mesh.prop.eff.Zn2,"Mesh.prop.eff.Zn2.csv",row.names=F)
  
  #Gillnet selectivity
  if(exists('Gillnet.selectivity'))  write.csv(Gillnet.selectivity,"Gillnet.selectivity.csv",row.names=F)
  if(exists('Gillnet.selectivity_len.age'))  write.csv(Gillnet.selectivity_len.age,"Gillnet.selectivity_len.age.csv",row.names=F)

}

fn.input.data.old=function(SP,Yr.assess,Conv.cal.mn.to.fin.mn,Historic.Ktch,Bin.size,What.Efrt) #previous version
{
  source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")
  source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/fn.fig.R")
  if(!"dplyr" %in% (.packages()))library(dplyr)
  #Define species
  if(SP=="WH")
  {
    Species=17003
    SP.LABELS=c("Whiskery shark")
    SP.LABELS1=c("Whiskery")
  }

  if(SP=="GM")
  {
    Species=17001
    SP.LABELS=c("Gummy shark")
    SP.LABELS1=c("Gummy")
  }

  if(SP=="BW")
  {
    Species=c(18003)
    SP.LABELS=c("Dusky shark")
    SP.LABELS1=c("Dusky")
  }

  if(SP=="TK")
  {
    Species=18007
    SP.LABELS=c("Sandbar shark")
    SP.LABELS1=c("Sandbar")
  }



  #----DATA SECTION------
  fn.subs=function(YEAR) substr(YEAR,start=3,stop=4)
  fn.finyr=function(Finyr)paste(Finyr-1,"-",fn.subs(Finyr),sep="")
  fn.month.to.fin.mn=function(dat)ifelse(dat>6,dat-6,dat+6)
  fn.exp.his=function(dat,dat1)
  {
    dat2=expand.grid(FINYEAR=dat$FINYEAR,MONTH=1:12)
    dat2=merge(dat2,dat1[,match(c("MONTH","prop"),names(dat1))],by="MONTH",all.x=T)
    dat2=merge(dat2,dat,by="FINYEAR",all.x=T)
    dat2=dat2[order(dat2$FINYEAR,dat2$MONTH),]
    dat2$LIVEWT.c=dat2$LIVEWT.c*dat2$prop
    return(dat2[,-match("prop",names(dat2))])
  }

  #1. CATCH

  #1.1. Commercial
  fn.in=function(NM)
  {
    read.csv(paste('C:/Matias/Analyses/Data_outs/',NM,sep=""),stringsAsFactors = F)
  }

  #1.1.1 Catch_WA Fisheries

  #Historic
  Hist.expnd=fn.in(NM='recons_Hist.expnd.csv')

  #Ammended reported catch including discards
  Data.monthly=fn.in(NM='recons_Data.monthly.csv')
  Data.monthly.north=fn.in(NM='recons_Data.monthly.north.csv')

  #TEPS
  TEPS_dusky=fn.in(NM='recons_TEPS_dusky.csv')


  #1.1.2. Catch of non WA Fisheries

  #Commonwealth GAB trawl and Western Tuna and Billfish Fisheries (WTBF)
  GAB.trawl_catch=fn.in(NM='recons_GAB.trawl_catch.csv')
  WTBF_catch=fn.in(NM='recons_WTBF_catch.csv')

  #SA Marine Scalefish fishery
  Whaler_SA=fn.in(NM='recons_Whaler_SA.csv')

  #Indonesian illegal fishing in Australia waters
  Indo_total.annual.ktch=fn.in(NM='recons_Indo.IUU.csv')         


  #Proportion of gummies by region
  Gummies.prop=read.csv("C:/Matias/Analyses/Catch and effort/Gummies.prop.csv",stringsAsFactors=F)


  #Select years of data
  First.year=sort(unique(Data.monthly$FINYEAR))[1]
  Last.year=Data.yr
  Used.yrs=as.numeric(c(substr(First.year,1,4),substr(Last.year,1,4)))
  Used.yrs=seq(Used.yrs[1],Used.yrs[2])
  Used.yrs=paste(Used.yrs,substr(Used.yrs+1,3,4),sep="-")

  #Northern Shark Fisheries (WANCSF and JANSF)
  Data.monthly.north=subset(Data.monthly.north,SPECIES%in%Species)

  #Select relevant years of data
  Data.monthly=subset(Data.monthly,FINYEAR%in%Used.yrs)
  Data.monthly.north=subset(Data.monthly.north,FINYEAR%in%Used.yrs)
  if(exists('WRL'))WRL=subset(WRL,Finyear%in%Used.yrs)


  #1.2 Rec fishing
  Rec.ktch=fn.in(NM='recons_recreational.csv')


  #Convert to tonnes if required
  if(KTCH.UNITS=="TONNES")
  {
    Data.monthly$LIVEWT.c=Data.monthly$LIVEWT.c/1000
    Data.monthly.north$LIVEWT.c=Data.monthly.north$LIVEWT.c/1000
    Hist.expnd$LIVEWT.c=Hist.expnd$LIVEWT.c/1000
    TEPS_dusky$LIVEWT.c=TEPS_dusky$LIVEWT.c/1000
    GAB.trawl_catch$LIVEWT.c=GAB.trawl_catch$LIVEWT.c/1000
    WTBF_catch$LIVEWT.c=WTBF_catch$LIVEWT.c/1000
    Whaler_SA$LIVEWT.c=Whaler_SA$LIVEWT.c/1000

    Rec.ktch$LIVEWT.c=Rec.ktch$LIVEWT.c/1000
  }


  #2. ABUNDANCE
  HNDL="C:/Matias/Analyses/Data_outs/"
  if(SP=="WH")
  {
    #by zone (relative)
    #monthly
    # Ab.indx.TDGDLF.West=read.csv(paste(HNDL,"Whiskery shark.annual.abundance.basecase.monthly.West_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn1=read.csv(paste(HNDL,"Whiskery Shark.annual.abundance.basecase.monthly.Zone1_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn2=read.csv(paste(HNDL,"Whiskery Shark.annual.abundance.basecase.monthly.Zone2_relative.csv",sep=""),stringsAsFactors=F)

    #daily
    Ab.indx.TDGDLF.West.daily=read.csv(paste(HNDL,"Whiskery Shark.annual.abundance.basecase.daily.West_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn1.daily=read.csv(paste(HNDL,"Whiskery Shark.annual.abundance.basecase.daily.Zone1_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn2.daily=read.csv(paste(HNDL,"Whiskery Shark.annual.abundance.basecase.daily.Zone2_relative.csv",sep=""),stringsAsFactors=F)


    #zones combined  (relative)
    #monthly
    Ab.indx.TDGDLF.all=read.csv(paste(HNDL,"Whiskery Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
    #daily
    Ab.indx.TDGDLF.all.daily=read.csv(paste(HNDL,"Whiskery Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

  }

  if(SP=="GM")
  {
    #by zone (relative)
    #monthly
    Ab.indx.TDGDLF.Zn2=read.csv(paste(HNDL,"Gummy Shark.annual.abundance.basecase.monthly.Zone2_relative.csv",sep=""),stringsAsFactors=F)
    #daily
    Ab.indx.TDGDLF.Zn2.daily=read.csv(paste(HNDL,"Gummy Shark.annual.abundance.basecase.daily.Zone2_relative.csv",sep=""),stringsAsFactors=F)


    #zones combined (relative)
    #monthly
    Ab.indx.TDGDLF.all=read.csv(paste(HNDL,"Gummy Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
    #daily
    Ab.indx.TDGDLF.all.daily=read.csv(paste(HNDL,"Gummy Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

  }

  if(SP=="BW")
  {
    #1. TDGDLF standardised CPUE (relative)
    #by zone
    #monthly
    Ab.indx.TDGDLF.West=read.csv(paste(HNDL,"Dusky Shark.annual.abundance.basecase.monthly.West_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn1=read.csv(paste(HNDL,"Dusky Shark.annual.abundance.basecase.monthly.Zone1_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn2=read.csv(paste(HNDL,"Dusky Shark.annual.abundance.basecase.monthly.Zone2_relative.csv",sep=""),stringsAsFactors=F)

    #daily
    Ab.indx.TDGDLF.West.daily=read.csv(paste(HNDL,"Dusky shark.annual.abundance.basecase.daily.West_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn1.daily=read.csv(paste(HNDL,"Dusky shark.annual.abundance.basecase.daily.Zone1_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn2.daily=read.csv(paste(HNDL,"Dusky shark.annual.abundance.basecase.daily.Zone2_relative.csv",sep=""),stringsAsFactors=F)

    #zones combined
    #monthly
    Ab.indx.TDGDLF.all=read.csv(paste(HNDL,"Dusky Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
    #daily
    Ab.indx.TDGDLF.all.daily=read.csv(paste(HNDL,"Dusky Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

    #2. Naturaliste survey (relative)
    Ab.index.Srvy.FixSt=read.csv(paste(HNDL,"Dusky.Srvy.FixSt.csv",sep=""),stringsAsFactors=F)
    Size.index.Srvy.FixSt=read.csv(paste(HNDL,"Dusky.Srvy.FixSt_size.csv",sep=""),stringsAsFactors=F)
  }

  if(SP=="TK")
  {
    #1. TDGDLF standardised CPUE
    #by zone (relative)
    #monthly
    Ab.indx.TDGDLF.West=read.csv(paste(HNDL,"Sandbar Shark.annual.abundance.basecase.monthly.West_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn1=read.csv(paste(HNDL,"Sandbar Shark.annual.abundance.basecase.monthly.Zone1_relative.csv",sep=""),stringsAsFactors=F)

    #daily
    Ab.indx.TDGDLF.West.daily=read.csv(paste(HNDL,"Sandbar Shark.annual.abundance.basecase.daily.West_relative.csv",sep=""),stringsAsFactors=F)
    Ab.indx.TDGDLF.Zn1.daily=read.csv(paste(HNDL,"Sandbar Shark.annual.abundance.basecase.daily.Zone1_relative.csv",sep=""),stringsAsFactors=F)

    #zones combined  (relative)
    #monthly
    Ab.indx.TDGDLF.all=read.csv(paste(HNDL,"Sandbar Shark.annual.abundance.basecase.monthly_relative.csv",sep=""),stringsAsFactors=F)
    #daily
    Ab.indx.TDGDLF.all.daily=read.csv(paste(HNDL,"Sandbar Shark.annual.abundance.basecase.daily_relative.csv",sep=""),stringsAsFactors=F)

    #2. Naturaliste survey (relative)
    Ab.index.Srvy.FixSt=read.csv(paste(HNDL,"Sandbar.Srvy.FixSt.csv",sep=""),stringsAsFactors=F)
    Size.index.Srvy.FixSt=read.csv(paste(HNDL,"Sandbar.Srvy.FixSt_size.csv",sep=""),stringsAsFactors=F)
  }


  #3. CATCH SIZE COMPOSITION
  setwd("C:/Matias/Data")
  #3.1 Heald 1987
  #description: Partial length in cm obtained from Perth fish marke.
  #             Unspecified fishing method, most likely gillnets.
  PL_Heald_1987=read.csv("Size_composition/Heald(1987).csv",stringsAsFactors=F)

  #3.2 Stevens 1990 (citation in Simpfendorfer & Donohue 1998)
  #description: TL (fish market) or FL (measured at sea) depending on period, in cm.
  #             gillnet of 6.5 and 7 inch mesh size.
  TL_FL_Stevens_1990=read.csv("Size_composition/Stevens_1990_size_comp_6.5_7_inch.csv",stringsAsFactors=F)

  #3.3.TDGDLF observing programs
  #description: FL (cm) composition observed as part of different research projects on commercial gillnet vessels.
  #           6.5 and 7 inch mesh combined (also available are data by mesh size). Souce: "Shark database"
  if(SP=="WH")      #FIX. all this is each species folder in Analysis/Data outs  !!!!
  {
    FL.TDGDFL.WC=read.csv("Size_composition/Whiskery.West.6.5.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn1=read.csv("Size_composition/Whiskery.Zone1.6.5.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn2=read.csv("Size_composition/Whiskery.Zone2.6.5.inch.raw.csv",stringsAsFactors=F)

    FL.TDGDFL.WC_7=read.csv("Size_composition/Whiskery.West.7.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn1_7=read.csv("Size_composition/Whiskery.Zone1.7.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn2_7=read.csv("Size_composition/Whiskery.Zone2.7.inch.raw.csv",stringsAsFactors=F)

    TDGDFL.size.numbers=read.csv("Size_composition/Whiskery.Numb_obs_size.freq.TDGDLF.csv",stringsAsFactors=F)
  }

  if(SP=="GM")
  {
    FL.TDGDFL.WC=read.csv("Size_composition/Gummy.West.6.5.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn1=read.csv("Size_composition/Gummy.Zone1.6.5.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn2=read.csv("Size_composition/Gummy.Zone2.6.5.inch.raw.csv",stringsAsFactors=F)

    FL.TDGDFL.WC_7=read.csv("Size_composition/Gummy.West.7.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn1_7=read.csv("Size_composition/Gummy.Zone1.7.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn2_7=read.csv("Size_composition/Gummy.Zone2.7.inch.raw.csv",stringsAsFactors=F)

    TDGDFL.size.numbers=read.csv("Size_composition/Gummy.Numb_obs_size.freq.TDGDLF.csv",stringsAsFactors=F)
  }

  if(SP=="BW")
  {
    FL.TDGDFL.WC=read.csv("Size_composition/Dusky.West.6.5.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn1=read.csv("Size_composition/Dusky.Zone1.6.5.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn2=read.csv("Size_composition/Dusky.Zone2.6.5.inch.raw.csv",stringsAsFactors=F)

    FL.TDGDFL.WC_7=read.csv("Size_composition/Dusky.West.7.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn1_7=read.csv("Size_composition/Dusky.Zone1.7.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn2_7=read.csv("Size_composition/Dusky.Zone2.7.inch.raw.csv",stringsAsFactors=F)

    TDGDFL.size.numbers=read.csv("Size_composition/Dusky.Numb_obs_size.freq.TDGDLF.csv",stringsAsFactors=F)
  }

  if(SP=="TK")
  {
    FL.TDGDFL.WC=read.csv("Size_composition/Sandbar.West.6.5.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn1=read.csv("Size_composition/Sandbar.Zone1.6.5.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn2=read.csv("Size_composition/Sandbar.Zone2.6.5.inch.raw.csv",stringsAsFactors=F)

    FL.TDGDFL.WC_7=read.csv("Size_composition/Sandbar.West.7.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn1_7=read.csv("Size_composition/Sandbar.Zone1.7.inch.raw.csv",stringsAsFactors=F)
    FL.TDGDFL.Zn2_7=read.csv("Size_composition/Sandbar.Zone2.7.inch.raw.csv",stringsAsFactors=F)

    TDGDFL.size.numbers=read.csv("Size_composition/Sandbar.Numb_obs_size.freq.TDGDLF.csv",stringsAsFactors=F)
  }

  #keep years with at least 10 observations from at least 10 shots by zone
  TDGDFL.size.numbers=subset(TDGDFL.size.numbers,N.observations>=Min.obs & N.shots>=Min.shts)
  This.yr.zn=with(TDGDFL.size.numbers,paste(FINYEAR,zone))


  #3.4 Pilbara trawl
  #description: FL (cm) composition observed on Pilbara trawl vessels. Source: "Shark database"
  FL_Sandbar_Pilbara_trawl=read.csv("Size_composition/Pilbara_Trawl_sandbar.csv")

  #3.5 NSF longline
  #description: FL (cm) composition observed on NSF longline vessels. Source: "Shark database"
  FL_Sandbar_NSF=read.csv("Size_composition/NSF.LONGLINE_sandbar.csv",stringsAsFactors=F)
  FL_dusky_NSF=read.csv("Size_composition/NSF.LONGLINE_dusky.csv",stringsAsFactors=F)

  #3.6 TEPS_TDGLDF
  Size.comp.Dusky.TEPS_TDGLDF=c(3.05,4,3,3.5,3.5,3.5,3.5,3,3,3,4,3)     #from Comments in TDGDLF returns


  #4. TAGGING

  #4.1 Conventional
  #note: there's also info by block but very few observations at this level....
  hndl="C:/Matias/Analyses/Data_outs/"
  if(SP=="BW")
  {
    #Individual based model
    Rel_rec_Conv.Tag=read.csv(paste(hndl,"Dusky shark_Con_tag_Ind.based.mod_Conv.Tag.csv",sep=""),stringsAsFactors=F)

    #at age
    Zn.rel_Conv.Tag=read.csv(paste(hndl,"BW_Zn.rel_Conv.Tag.csv",sep=""),stringsAsFactors=F)   #FIX THIS
    Zn.rec_Conv.Tag=read.csv(paste(hndl,"BW_Zn.rec_Conv.Tag.csv",sep=""),stringsAsFactors=F)

    #at size
    #all sizes
    Zn.rel_Conv.Tag_size=read.csv(paste(hndl,"BW_Zn.rel_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size=read.csv(paste(hndl,"BW_Zn.rec_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)

    #adults and juvenlies
    Zn.rel_Conv.Tag_size_adu=read.csv(paste(hndl,"BW_Zn.rel.adul_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rel_Conv.Tag_size_juv=read.csv(paste(hndl,"BW_Zn.rel.juv_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size_adu=read.csv(paste(hndl,"BW_Zn.rec.adul_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size_juv=read.csv(paste(hndl,"BW_Zn.rec.juv_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Smallest_size_tagged=read.csv(paste(hndl,"BW_Smallest.size_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
  }

  if(SP=="TK")
  {
    #Individual based model
    Rel_rec_Conv.Tag=read.csv(paste(hndl,"TK_Ind.based.mod_Conv.Tag.csv",sep=""),stringsAsFactors=F)

    #at age
    Zn.rel_Conv.Tag=read.csv(paste(hndl,"TK_Zn.rel_Conv.Tag.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag=read.csv(paste(hndl,"TK_Zn.rec_Conv.Tag.csv",sep=""),stringsAsFactors=F)

    #at size
    #all sizes
    Zn.rel_Conv.Tag_size=read.csv(paste(hndl,"TK_Zn.rel_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size=read.csv(paste(hndl,"TK_Zn.rec_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)

    #adults and juvenlies
    Zn.rel_Conv.Tag_size_adu=read.csv(paste(hndl,"TK_Zn.rel.adul_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rel_Conv.Tag_size_juv=read.csv(paste(hndl,"TK_Zn.rel.juv_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size_adu=read.csv(paste(hndl,"TK_Zn.rec.adul_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size_juv=read.csv(paste(hndl,"TK_Zn.rec.juv_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Smallest_size_tagged=read.csv(paste(hndl,"TK_Smallest.size_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
  }

  if(SP=="GM")
  {
    #Individual based model
    Rel_rec_Conv.Tag=read.csv(paste(hndl,"GM_Ind.based.mod_Conv.Tag.csv",sep=""),stringsAsFactors=F)

    #at age
    Zn.rel_Conv.Tag=read.csv(paste(hndl,"GM_Zn.rel_Conv.Tag.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag=read.csv(paste(hndl,"GM_Zn.rec_Conv.Tag.csv",sep=""),stringsAsFactors=F)

    #at size
    #all sizes
    Zn.rel_Conv.Tag_size=read.csv(paste(hndl,"GM_Zn.rel_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size=read.csv(paste(hndl,"GM_Zn.rec_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)

    #adults and juvenlies
    Zn.rel_Conv.Tag_size_adu=read.csv(paste(hndl,"GM_Zn.rel.adul_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rel_Conv.Tag_size_juv=read.csv(paste(hndl,"GM_Zn.rel.juv_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size_adu=read.csv(paste(hndl,"GM_Zn.rec.adul_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size_juv=read.csv(paste(hndl,"GM_Zn.rec.juv_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Smallest_size_tagged=read.csv(paste(hndl,"GM_Smallest.size_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
  }

  if(SP=="WH")
  {
    #Individual based model
    Rel_rec_Conv.Tag=read.csv(paste(hndl,"WH_Ind.based.mod_Conv.Tag.csv",sep=""),stringsAsFactors=F)

    #at age
    Zn.rel_Conv.Tag=read.csv(paste(hndl,"WH_Zn.rel_Conv.Tag.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag=read.csv(paste(hndl,"WH_Zn.rec_Conv.Tag.csv",sep=""),stringsAsFactors=F)

    #at size
    #all sizes
    Zn.rel_Conv.Tag_size=read.csv(paste(hndl,"WH_Zn.rel_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size=read.csv(paste(hndl,"WH_Zn.rec_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)

    #adults and juvenlies
    Zn.rel_Conv.Tag_size_adu=read.csv(paste(hndl,"WH_Zn.rel.adul_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rel_Conv.Tag_size_juv=read.csv(paste(hndl,"WH_Zn.rel.juv_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size_adu=read.csv(paste(hndl,"WH_Zn.rec.adul_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Zn.rec_Conv.Tag_size_juv=read.csv(paste(hndl,"WH_Zn.rec.juv_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
    Smallest_size_tagged=read.csv(paste(hndl,"WH_Smallest.size_Conv.Tag_size.csv",sep=""),stringsAsFactors=F)
  }


  #4.2 Acoustic           FIX path, no in Analysis/Data_outs/species
  #note: there's also info by block but very few observations at this level....

  if(SP=="BW")
  {
    #Taylor 2011 approach
    Zn.rel_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Dusky_Zn.rel_Acous.Tag.csv",stringsAsFactors=F)
    Zn.rec_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Dusky_Zn.rec_Acous.Tag.csv",stringsAsFactors=F)

    #Proportion of time approach
    Zn.rel_Acous.Tag.prop=read.csv("Tagging/Pop dyn model/Acoustic/Dusky_Zn.rel_Acous.Tag.prop.csv",stringsAsFactors=F)
    Zn.rec_Acous.Tag.prop=read.csv("Tagging/Pop dyn model/Acoustic/Dusky_Zn.rec_Acous.Tag.prop.csv",stringsAsFactors=F)

    #Individual based model
    Indiv_based_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Dusky_Acous.Tag.Ind_based.csv",stringsAsFactors=F)

  }

  if(SP=="TK")
  {
    #Taylor 2011 approach
    Zn.rel_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Thickskin_Zn.rel_Acous.Tag.csv",stringsAsFactors=F)
    Zn.rec_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Thickskin_Zn.rec_Acous.Tag.csv",stringsAsFactors=F)

    #Proportion of time approach
    Zn.rel_Acous.Tag.prop=read.csv("Tagging/Pop dyn model/Acoustic/Thickskin_Zn.rel_Acous.Tag.prop.csv",stringsAsFactors=F)
    Zn.rec_Acous.Tag.prop=read.csv("Tagging/Pop dyn model/Acoustic/Thickskin_Zn.rec_Acous.Tag.prop.csv",stringsAsFactors=F)

    #Individual based model
    Indiv_based_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Thickskin_Acous.Tag.Ind_based.csv",stringsAsFactors=F)
  }

  if(SP=="GM")
  {
    #Taylor 2011 approach
    Zn.rel_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Gummy_Zn.rel_Acous.Tag.csv",stringsAsFactors=F)
    Zn.rec_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Gummy_Zn.rec_Acous.Tag.csv",stringsAsFactors=F)

    #Proportion of time approach
    Zn.rel_Acous.Tag.prop=read.csv("Tagging/Pop dyn model/Acoustic/Gummy_Zn.rel_Acous.Tag.prop.csv",stringsAsFactors=F)
    Zn.rec_Acous.Tag.prop=read.csv("Tagging/Pop dyn model/Acoustic/Gummy_Zn.rec_Acous.Tag.prop.csv",stringsAsFactors=F)

    #Individual based model
    Indiv_based_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Gummy_Acous.Tag.Ind_based.csv",stringsAsFactors=F)
  }

  if(SP=="WH")
  {
    #Taylor 2011 approach
    Zn.rel_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Whiskery_Zn.rel_Acous.Tag.csv",stringsAsFactors=F)
    Zn.rec_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Whiskery_Zn.rec_Acous.Tag.csv",stringsAsFactors=F)

    #Proportion of time approach
    Zn.rel_Acous.Tag.prop=read.csv("Tagging/Pop dyn model/Acoustic/Whiskery_Zn.rel_Acous.Tag.prop.csv",stringsAsFactors=F)
    Zn.rec_Acous.Tag.prop=read.csv("Tagging/Pop dyn model/Acoustic/Whiskery_Zn.rec_Acous.Tag.prop.csv",stringsAsFactors=F)

    #Individual based model
    Indiv_based_Acous.Tag=read.csv("Tagging/Pop dyn model/Acoustic/Whiskery_Acous.Tag.Ind_based.csv",stringsAsFactors=F)
  }

  #Reported recaptures and releases of acoustic tagging
  library(RODBC)
  library(lubridate)
  setwd("U:/Shark")  # working directory
  channel <- odbcConnectAccess2007("Sharks.mdb")
  Rep.Recap=sqlFetch(channel, "Tag data", colnames = F)
  close(channel)

  names(Rep.Recap)[match(c("Captured?","ATAG NO"),names(Rep.Recap))]=c("Recaptured","ATAG_NO")
  Rep.Recap=subset(Rep.Recap,!is.na(ATAG_NO))
  Rep.Recap$RELLATDECDEG=-with(Rep.Recap,REL_LATD+(REL_LATM/60))
  Rep.Recap$RECLATDECDEG=-with(Rep.Recap,CAP_LATD+(CAP_LATM/60))
  Rep.Recap$RELLNGDECDEG=with(Rep.Recap,REL_LNGD+(REL_LNGM/60))
  Rep.Recap$RECLNGDECDEG=with(Rep.Recap,CAP_LNGD+(CAP_LNGM/60))
  No.recap.pos=c(30926,30950,29436,30929,29463)   #Reported by J. Tindal but whithout recaptured info
  No.recap.date=c(30926,30929,30928,29499,29463)   #Reported by J. Tindal but whithout recaptured info
  #assign Zone 2 recapture position for some unknowns that were recaptured in that zone
  Rep.Recap$RECLATDECDEG=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.pos,-35,RECLATDECDEG))
  Rep.Recap$RECLNGDECDEG=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.pos,119,RECLNGDECDEG))
  Rep.Recap$DATE_CAPTR1=as.character(Rep.Recap$DATE_CAPTR)
  Rep.Recap$DATE_CAPTR1=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.date
                                              & is.na(DATE_CAPTR1),"2013-06-01",DATE_CAPTR1))
  Rep.Recap$DATE_CAPTR1=as.POSIXlt(Rep.Recap$DATE_CAPTR1)
  Rep.Recap$DATE_CAPTR=Rep.Recap$DATE_CAPTR1
  Rep.Recap=Rep.Recap[,-match("DATE_CAPTR1",names(Rep.Recap))]
  Rep.Recap$Rep.pos=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.pos,"No rec. pos","OK"))
  Rep.Recap$Rep.date=with(Rep.Recap,ifelse(ATAG_NO%in%No.recap.date,"No rec. date","OK"))

  #add recaptured gummy by J. Cooke but with no recapture information
  Add.Jeff=Rep.Recap[1,]
  Add.Jeff[,]=NA
  Add.Jeff$SPECIES="GM"
  Add.Jeff$Recaptured="YES"
  Add.Jeff$Rep.pos="No rec. pos"
  Add.Jeff$Rep.date="OK"
  Add.Jeff$DATE_CAPTR=as.POSIXlt("2015-05-15")
  Rep.Recap=rbind(Rep.Recap,Add.Jeff)

  this=match(c("SPECIES","ATAG_NO","SEX","Recaptured","Rep.pos","Rep.date","RELEASE DATE","DATE_CAPTR","RELLATDECDEG",
               "RELLNGDECDEG","RECLATDECDEG","RECLNGDECDEG"),names(Rep.Recap))
  Rep.Recap=Rep.Recap[order(Rep.Recap$SPECIES,Rep.Recap$ATAG_NO),this]

  Rep.Recap$year.rel=year(Rep.Recap$"RELEASE DATE")
  Rep.Recap=subset(Rep.Recap,SPECIES==SP)
  Rep.Recap$year.rec=year(Rep.Recap$"DATE_CAPTR")
  Rep.Recap$Recaptured=as.character(Rep.Recap$Recaptured)


  #5. Effort by zone (GN plus long line equivalent effort)
  #(in 1000 km gn days)
  if(What.Efrt=="km.gn.days") Eff.zn=read.csv("C:/Matias/Analyses/Data_outs/Annual.zone.eff.days.csv",stringsAsFactors=F)

  #(in 1000 km gn hours)
  if(What.Efrt=="km.gn.hours") Eff.zn=read.csv("C:/Matias/Analyses/Data_outs/Annual.zone.eff.hours.csv",stringsAsFactors=F)


  #6.  Age and growth data for whiskery only
  if(SP=="WH") Age.growth=read.csv("C:/Matias/Data/Age and growth/Simpfen.data.csv",stringsAsFactors=F)
  if(SP=="GM") Age.growth=read.csv("C:/Matias/Data/Age and growth/Terry/Gummy_Terry.csv",stringsAsFactors=F)
  if(SP=="BW") Age.growth=read.csv("C:/Matias/Data/Age and growth/Dusky.csv",stringsAsFactors=F)
  if(SP=="TK") Age.growth=read.csv("C:/Matias/Data/Age and growth/Sandbar.csv",stringsAsFactors=F)


  #7. Proportional effort by mesh size
  Mesh.prop.eff=read.csv("C:/Matias/Analyses/Catch and effort/mesh.proportional.effort.csv",stringsAsFactors=F)
  Mesh.prop.eff.West=read.csv("C:/Matias/Analyses/Catch and effort/mesh.proportional.effort.West.csv",stringsAsFactors=F)
  Mesh.prop.eff.Zn1=read.csv("C:/Matias/Analyses/Catch and effort/mesh.proportional.effort.Zone1.csv",stringsAsFactors=F)
  Mesh.prop.eff.Zn2=read.csv("C:/Matias/Analyses/Catch and effort/mesh.proportional.effort.Zone2.csv",stringsAsFactors=F)


  #8.  Standardised Mean size
  if(SP=="WH")
  {
    Avr.wt.yr=read.csv("C:/Matias/Analyses/Data_outs/Whiskery Shark.annual.mean.size_relative.csv",stringsAsFactors=F)

    Avr.wt.yr_west=read.csv("C:/Matias/Analyses/Data_outs/Whiskery Shark.annual.mean.size_relative_west.csv",stringsAsFactors=F)
    Avr.wt.yr_zn1=read.csv("C:/Matias/Analyses/Data_outs/Whiskery Shark.annual.mean.size_relative_zone1.csv",stringsAsFactors=F)
    Avr.wt.yr_zn2=read.csv("C:/Matias/Analyses/Data_outs/Whiskery Shark.annual.mean.size_relative_zone2.csv",stringsAsFactors=F)
    Avr.wt.yr_west$zone="West"
    Avr.wt.yr_zn1$zone="Zone1"
    Avr.wt.yr_zn2$zone="Zone2"
    Avr.wt.yr.zn=rbind(Avr.wt.yr_west,Avr.wt.yr_zn1,Avr.wt.yr_zn2)
  }
  if(SP=="GM")
  {
    Avr.wt.yr=read.csv("C:/Matias/Analyses/Data_outs/Gummy Shark.annual.mean.size_relative.csv",stringsAsFactors=F)

    Avr.wt.yr_zn2=read.csv("C:/Matias/Analyses/Data_outs/Gummy Shark.annual.mean.size_relative_zone2.csv",stringsAsFactors=F)
    Avr.wt.yr_zn2$zone="Zone2"
    Avr.wt.yr.zn=Avr.wt.yr_zn2


  }
  if(SP=="TK")
  {
    Avr.wt.yr=read.csv("C:/Matias/Analyses/Data_outs/Sandbar Shark.annual.mean.size_relative.csv",stringsAsFactors=F)

    Avr.wt.yr_west=read.csv("C:/Matias/Analyses/Data_outs/Sandbar Shark.annual.mean.size_relative_west.csv",stringsAsFactors=F)
    Avr.wt.yr_zn1=read.csv("C:/Matias/Analyses/Data_outs/Sandbar Shark.annual.mean.size_relative_zone1.csv",stringsAsFactors=F)
    Avr.wt.yr_west$zone="West"
    Avr.wt.yr_zn1$zone="Zone1"
    Avr.wt.yr.zn=rbind(Avr.wt.yr_west,Avr.wt.yr_zn1)

  }
  if(SP=="BW")
  {
    Avr.wt.yr=read.csv("C:/Matias/Analyses/Data_outs/Dusky Shark.annual.mean.size_relative.csv",stringsAsFactors=F)

    Avr.wt.yr_west=read.csv("C:/Matias/Analyses/Data_outs/Dusky Shark.annual.mean.size_relative_west.csv",stringsAsFactors=F)
    Avr.wt.yr_zn1=read.csv("C:/Matias/Analyses/Data_outs/Dusky Shark.annual.mean.size_relative_zone1.csv",stringsAsFactors=F)
    Avr.wt.yr_zn2=read.csv("C:/Matias/Analyses/Data_outs/Dusky Shark.annual.mean.size_relative_zone2.csv",stringsAsFactors=F)
    Avr.wt.yr_west$zone="West"
    Avr.wt.yr_zn1$zone="Zone1"
    Avr.wt.yr_zn2$zone="Zone2"
    Avr.wt.yr.zn=rbind(Avr.wt.yr_west,Avr.wt.yr_zn1,Avr.wt.yr_zn2)

  }


  #----PARAMETERS SECTIONS -------

  source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/Organise input parameters.R")
  ParS=fn.input.pars(SP,add.growth.cv="NO",add.Mrt.age="NO")$pars

  Max.age=ParS$Max.Age.F

  #2. LENGTHS AND WEIGHTS CONVERSION FACTORS
  #PL to TL
  if(!SP=="TK")
  {
    PL_a=ParS$PL.to.TL$a
    PL_b=ParS$PL.to.TL$b
  }
  if(SP=="TK")
  {
    PL_a=1
    PL_b=0
  }

  #TL to FL
  a.TL=ParS$TL.to.FL$a
  b.TL=ParS$TL.to.FL$b

  #FL=b+a*PL
  if(length(ParS$PL.to.FL)>0)
  {
    b.PL=ParS$PL.to.FL$b
    a.PL=ParS$PL.to.FL$a
  }

  #TL-TWT (kg)
  bwt=ParS$TL.to.TwT.F$b
  awt=ParS$TL.to.TwT.F$a

  #3. Min and Max observed FL in catch
  Min.FL=list(WH=70,GM=60,BW=60,TK=40)
  Mx.FL=list(WH=135,GM=160,BW=280,TK=160)
  id=which(names(Min.FL)==SP)
  Min.FL=Min.FL[[id]]
  Mx.FL=Mx.FL[[id]]



  #------PROCEDURE SECTION-------

  #Convert calendar month to financial month to match financial year
  if(Conv.cal.mn.to.fin.mn=="YES")
  {

    #catch
    Data.monthly$MONTH.calendar=Data.monthly$MONTH
    Data.monthly$MONTH=fn.month.to.fin.mn(Data.monthly$MONTH)

    Data.monthly.north$MONTH.calendar=Data.monthly.north$MONTH
    Data.monthly.north$MONTH=fn.month.to.fin.mn(Data.monthly.north$MONTH)

    #size
    FL.TDGDFL.WC$Month=fn.month.to.fin.mn(FL.TDGDFL.WC$Month)
    FL.TDGDFL.Zn1$Month=fn.month.to.fin.mn(FL.TDGDFL.Zn1$Month)
    FL.TDGDFL.Zn2$Month=fn.month.to.fin.mn(FL.TDGDFL.Zn2$Month)

    #avg weight
    Average.wt$month=fn.month.to.fin.mn(Average.wt$month)

  }

  fn.weight=function(TL,bwt,awt) bwt*TL^awt


  #1. CATCH

  #1.1 Historic catches in TDGDLF
  historic=Hist.expnd%>%filter(SPECIES%in%Species)


  #1.2 TDGDLF reported catch
  TDGDLF=subset(Data.monthly,SPECIES%in%Species & METHOD%in%c("GN","LL") & Estuary=="NO",
                select=c(FINYEAR,LIVEWT.c,SPECIES,zone))

  #combine TDGDLF (i.e. GN and LL) with other gears reported with TDGDLF
  other_ktch=Data.monthly%>%
    filter(SPECIES%in%Species & !METHOD%in%c("GN","LL"))%>%
    dplyr::select(c(FINYEAR,LIVEWT.c,SPECIES,zone))

  #correct gummy catches by proportion of M antarcticus per bioregion
  if(SP=="GM")
  {
    Gummy.prop=subset(Gummies.prop,Bioregion=="WC" & SPECIES=="GM" ,select=Prop)$Prop
    TDGDLF$LIVEWT.c=with(TDGDLF,ifelse(zone%in%c("West","Zone1"),LIVEWT.c*Gummy.prop,LIVEWT.c))
  }


  #1.3 Northern Shark Fisheries (WANCSF and JANSF)
  #Check nonsense catch from Data.monthly.north
  if(SP%in%c("BW","TK"))
  {
    fn.check.NSF=function(a)
    {
      unib=sort(unique(a$BLOCKX))
      this.b=unib[unib<26000]
      a=subset(a,BLOCKX%in%this.b)
      a$YEAR=as.numeric(substr(a$FINYEAR,1,4))
      b=aggregate(LIVEWT.c~YEAR,a,sum)
      dd=b$LIVEWT.c
      plot(dd,type='l',lwd=2)
    }

    fn.check.NSF(Data.monthly.north)
    NSF=subset(Data.monthly.north,SPECIES%in%Species & METHOD%in%c("LL","GN"),
               select=c(FINYEAR,LIVEWT.c,SPECIES,zone))

    other_ktch.north=subset(Data.monthly.north,SPECIES%in%Species & !METHOD%in%c("LL","GN"),
                            select=c(FINYEAR,LIVEWT.c,SPECIES,zone))
    other_ktch=rbind(other_ktch,other_ktch.north)
  }

  #construct species data frames
  Yrs.dummy=as.numeric(substr(sort(unique(Data.monthly$FINYEAR)),1,4))
  x=substr(Yrs.dummy[2]:Yrs.dummy[length(Yrs.dummy)],start=3,stop=4)
  FinYrs=paste(Yrs.dummy[1]:Yrs.dummy[length(Yrs.dummy)-1],"-",x,sep="")


  #1.4 Recreational catch
  if(SP=='WH') SPEC='Whiskery Shark'
  if(SP=='GM') SPEC='Gummy Sharks'
  if(SP=='BW') SPEC='Dusky Whaler'
  if(SP=='TK') SPEC='Sandbar Shark'
  Rec.ktch=subset(Rec.ktch,Common.Name==SPEC)
  if(SP=='GM') Rec.ktch=subset(Rec.ktch,zone%in%c("South Coast","West Coast"))


  #1.5 Non-WA Fisheries
  GAB_trawl=GAB.trawl_catch%>%filter(SPECIES%in%Species)
  WTBF=WTBF_catch%>%filter(SPECIES%in%Species)
  Whaler_SA=Whaler_SA%>%filter(SPECIES%in%Species)

  #-Put all catches in list (currently not including 'historic')
  if(SP=="WH") catch=list(TDGDLF=TDGDLF,Other=other_ktch,Rec=rec.ktch)
  if(SP=="GM") catch=list(TDGDLF=TDGDLF,Other=other_ktch,Rec=rec.ktch)
  if(SP=="BW")
  {
    if(Ktch.source=="ALL")
    {
      catch=list(TDGDLF=TDGDLF,NSF=NSF,Other=other_ktch,WRL=Dusky.Tot.Ktch.WRL,
                 WTBF=WTBF,TEPS=TEPS_dusky,SA_marine=Whaler_SA,Rec=rec.ktch)
    }
    if(Ktch.source=="WA.only")
    {
      catch=list(TDGDLF=TDGDLF,NSF=NSF,Other=other_ktch,TEPS=TEPS_dusky,Rec=rec.ktch)
    }

    if(ADD.bronzy_dusky=="YES" & Ktch.source=="ALL")
    {
      catch=list(TDGDLF=TDGDLF,NSF=NSF,Other=other_ktch,WRL=Dusky.Tot.Ktch.WRL,
                 GAB_trawl=Dusky_GAB_trawl,WTBF=WTBF,TEPS=TEPS_dusky,
                 SA_marine=Whaler_SA,Rec=rec.ktch)
    }
  }
  if(SP=="TK")
  {
    if(Ktch.source=="ALL") catch=list(TDGDLF=TDGDLF,NSF=NSF,Other=other_ktch,WTBF=WTBF,Rec=rec.ktch)
    if(Ktch.source=="WA.only") catch=list(TDGDLF=TDGDLF,NSF=NSF,Other=other_ktch,Rec=rec.ktch)
  }


  #2. Size composition
  if(SP=='BW') SP.size="Dusky"
  if(SP=='GM') SP.size="Gummy"
  if(SP=='WH') SP.size="Whiskery"


  #2.1 TDGDLF

  #2.1.1 Standardise length units
  PL_to_TL=function(PL,a,b) TL=(PL-b)/a
  TL_to_FL=function(TL,a,b) FL=(TL-b)/a
  FL_to_TL=function(FL,a,b) TL=FL*a+b

  if(!SP=="TK")
  {
    names(PL_Heald_1987)[1]="PL"
    PL_Heald_1987=subset(PL_Heald_1987,Species==SP.size)
    PL_Heald_1987$TL=PL_to_TL(PL_Heald_1987$PL,PL_a,PL_b)
    PL_Heald_1987$Month=with(PL_Heald_1987,ifelse(Season=="Summer",1,ifelse(Season=="Spring",10,NA)))
    PL_Heald_1987$FL=round(TL_to_FL(PL_Heald_1987$TL,a.TL,b.TL))
    PL_Heald_1987$Type="Fish market"

    TL_FL_Stevens_1990=subset(TL_FL_Stevens_1990,Species==SP.size)
    names(TL_FL_Stevens_1990)[1:2]=c("TL","FL")
    TL_FL_Stevens_1990$FL=round(with(TL_FL_Stevens_1990,ifelse(is.na(FL),TL_to_FL(TL,a.TL,b.TL),FL)))
  }



  #2.2. Aggregate size frequencies by financial year and month
  fn.add.missing.size=function(dat,Min,Max,interval)
  {
    Rango=c(Min,Max)
    SEQ=seq(Rango[1],Rango[2],interval)
    dat=dat[order(dat$year),]
    dat$FL.bin=floor(dat$FL/interval)*interval
    dat$FL.bin=factor(dat$FL.bin,levels=SEQ)

    dat$FINYEAR=as.character(with(dat,
                                  ifelse(Month%in%1:6,paste(year-1,"-",fn.subs(year),sep=""),
                                         ifelse(Month%in%7:12,paste(year,"-",fn.subs(year+1),sep=""),NA))))

    if(Conv.cal.mn.to.fin.mn=="YES") dat$Month=fn.month.to.fin.mn(dat$Month)

    #   dat=aggregate(Number~FL.bin+FINYEAR+Type,dat,sum)
    #   if(!is.null(dat$Sex))dat=aggregate(Number~FL.bin+FINYEAR+Type+Sex,dat,sum)
    #
    #   dat2=subset(dat,Number>0,select=c(FL.bin,FINYEAR,Number))
    #   Table=reshape(dat2,v.names = "Number", idvar = "FINYEAR",timevar = "FL.bin", direction = "wide")
    #   X=colnames(Table)[2:ncol(Table)]
    #   colnames(Table)[2:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
    #   ID=which(!levels(dat$FL.bin)%in%colnames(Table)[2:ncol(Table)])
    #   ADD=levels(dat$FL.bin)[ID]
    #   ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
    #   names(ADD1)=ADD
    #   Table=cbind(Table,ADD1)
    #   Table=Table[,c(1,match(levels(dat$FL.bin),names(Table)))]
    #   Table[is.na(Table)]=0
    #
    dat=aggregate(Number~FL.bin+Month+FINYEAR+Type+Sex,dat,sum)

    dat2=subset(dat,Number>0,select=c(FL.bin,FINYEAR,Number,Month,Sex))
    Table=reshape(dat2,v.names = "Number", idvar = c("FINYEAR","Month","Sex"),timevar = "FL.bin", direction = "wide")
    X=colnames(Table)[4:ncol(Table)]
    colnames(Table)[4:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
    ID=which(!levels(dat$FL.bin)%in%colnames(Table)[4:ncol(Table)])
    ADD=levels(dat$FL.bin)[ID]
    ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
    names(ADD1)=ADD
    Table=cbind(Table,ADD1)
    Table=Table[,c(1:3,match(levels(dat$FL.bin),names(Table)))]
    Table[is.na(Table)]=0
    Table=Table[order(Table$Sex,Table$FINYEAR,Table$Month),]

    return(Table)

  }

  fn.add.missing.size2=function(dat,Min,Max,interval)
  {
    Rango=c(Min,Max)
    SEQ=seq(Rango[1],Rango[2],interval)
    dat=dat[order(dat$year),]
    dat$FL.bin=floor(dat$FL/interval)*interval
    dat$FL.bin=factor(dat$FL.bin,levels=SEQ)
    dat$Sex=dat$SEX
    dat$Number=1
    dat=subset(dat,!Sex=="U")

    dat=aggregate(Number~FL.bin+Month+FINYEAR+Sex,dat,sum)

    dat2=subset(dat,Number>0,select=c(FL.bin,FINYEAR,Number,Month,Sex))
    Table=reshape(dat2,v.names = "Number", idvar = c("FINYEAR","Month","Sex"),timevar = "FL.bin", direction = "wide")
    X=colnames(Table)[4:ncol(Table)]
    colnames(Table)[4:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
    ID=which(!levels(dat$FL.bin)%in%colnames(Table)[4:ncol(Table)])
    ADD=levels(dat$FL.bin)[ID]
    ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
    names(ADD1)=ADD
    Table=cbind(Table,ADD1)
    Table=Table[,c(1:3,match(levels(dat$FL.bin),names(Table)))]
    Table[is.na(Table)]=0
    Table=Table[order(Table$Sex,Table$FINYEAR,Table$Month),]
    return(Table)

  }

  fn.add.missing.size3=function(dat,Min,Max,interval)
  {
    Rango=c(Min,Max)
    SEQ=seq(Rango[1],Rango[2],interval)
    dat=dat[order(dat$year),]
    dat$TL.bin=floor(dat$TL/interval)*interval
    dat$TL.bin=factor(dat$TL.bin,levels=SEQ)
    dat$Sex=dat$SEX
    dat$Number=1
    dat=subset(dat,!Sex=="U")

    dat=aggregate(Number~TL.bin+Month+FINYEAR+Sex,dat,sum)

    dat2=subset(dat,Number>0,select=c(TL.bin,FINYEAR,Number,Month,Sex))
    Table=reshape(dat2,v.names = "Number", idvar = c("FINYEAR","Month","Sex"),timevar = "TL.bin", direction = "wide")
    X=colnames(Table)[4:ncol(Table)]
    colnames(Table)[4:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
    ID=which(!levels(dat$TL.bin)%in%colnames(Table)[4:ncol(Table)])
    ADD=levels(dat$TL.bin)[ID]
    ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
    names(ADD1)=ADD
    Table=cbind(Table,ADD1)
    Table=Table[,c(1:3,match(levels(dat$TL.bin),names(Table)))]
    Table[is.na(Table)]=0
    Table=Table[order(Table$Sex,Table$FINYEAR,Table$Month),]
    return(Table)

  }

  #Heald
  if(!SP=="TK")
  {
    PL_Heald_1987$Sex="U"
    PL_Heald_1987$year=PL_Heald_1987$Year
    PL_Heald_1987$FINYEAR=as.character(with(PL_Heald_1987,
                                            ifelse(Season%in%c("Spring","Summer"),paste(Year-1,"-",fn.subs(Year),sep=""),NA)))
    PL_Heald_1987=subset(PL_Heald_1987,FL>=Min.FL)
    FL_Heald_1987=fn.add.missing.size(PL_Heald_1987,Min.FL,Mx.FL,Bin.size)
    N.Heald.1987=rowSums(FL_Heald_1987[,4:ncol(FL_Heald_1987)])
    names(N.Heald.1987)=FL_Heald_1987$Finyear

  }

  #Stevens
  if(SP%in%c("WH","BW"))
  {
    TL_FL_Stevens_1990$year=TL_FL_Stevens_1990$Year
    FL_Stevens_1990=fn.add.missing.size(TL_FL_Stevens_1990,Min.FL,Mx.FL,Bin.size)
    N.Stevens.1990=rowSums(FL_Stevens_1990[,4:ncol(FL_Stevens_1990)])
    names(N.Stevens.1990)=FL_Stevens_1990$Finyear
  }

  #convert FL to TL
  FL.TDGDFL.WC$TL=round(FL_to_TL(FL.TDGDFL.WC$FL,a.TL,b.TL))
  FL.TDGDFL.Zn1$TL=round(FL_to_TL(FL.TDGDFL.Zn1$FL,a.TL,b.TL))
  FL.TDGDFL.Zn2$TL=round(FL_to_TL(FL.TDGDFL.Zn2$FL,a.TL,b.TL))

  FL.TDGDFL.WC_7$TL=round(FL_to_TL(FL.TDGDFL.WC_7$FL,a.TL,b.TL))
  FL.TDGDFL.Zn1_7$TL=round(FL_to_TL(FL.TDGDFL.Zn1_7$FL,a.TL,b.TL))
  FL.TDGDFL.Zn2_7$TL=round(FL_to_TL(FL.TDGDFL.Zn2_7$FL,a.TL,b.TL))

  #TDGDLF observing
  # FL_observers.WC=fn.add.missing.size2(FL.TDGDFL.WC,Min.FL,Mx.FL,Bin.size)
  # FL_observers.Zn1=fn.add.missing.size2(FL.TDGDFL.Zn1,Min.FL,Mx.FL,Bin.size)
  # FL_observers.Zn2=fn.add.missing.size2(FL.TDGDFL.Zn2,Min.FL,Mx.FL,Bin.size)

  Min.TL=round(FL_to_TL(Min.FL,a.TL,b.TL)/10)*10
  Mx.TL=round(FL_to_TL(Mx.FL,a.TL,b.TL))
  TL_observers.WC=fn.add.missing.size3(FL.TDGDFL.WC,Min.TL,Mx.TL,Bin.size)
  TL_observers.Zn1=fn.add.missing.size3(FL.TDGDFL.Zn1,Min.TL,Mx.TL,Bin.size)
  TL_observers.Zn2=fn.add.missing.size3(FL.TDGDFL.Zn2,Min.TL,Mx.TL,Bin.size)

  TL_observers.WC_7=fn.add.missing.size3(FL.TDGDFL.WC_7,Min.TL,Mx.TL,Bin.size)
  TL_observers.Zn1_7=fn.add.missing.size3(FL.TDGDFL.Zn1_7,Min.TL,Mx.TL,Bin.size)
  TL_observers.Zn2_7=fn.add.missing.size3(FL.TDGDFL.Zn2_7,Min.TL,Mx.TL,Bin.size)



  if(SP=="TK")
  {
    #2.2 Pilbara observing
    FL_Sandbar_Pilbara_trawl$FINYEAR=as.character(with(FL_Sandbar_Pilbara_trawl,
                                                       ifelse(Month%in%1:6,paste(year-1,"-",fn.subs(year),sep=""),
                                                              ifelse(Month%in%7:12,paste(year,"-",fn.subs(year+1),sep=""),NA))))
    if(Conv.cal.mn.to.fin.mn=="YES") FL_Sandbar_Pilbara_trawl$Month=fn.month.to.fin.mn(FL_Sandbar_Pilbara_trawl$Month)

    FL_Sandbar_Pilbara_trawl$FL=round(FL_to_TL(FL_Sandbar_Pilbara_trawl$FL,a.TL,b.TL))
    Sandbar_FL_Pil.trwl_observers=fn.add.missing.size2(FL_Sandbar_Pilbara_trawl,Min.FL,round(FL_to_TL(Mx.FL,a.TL,b.TL)),Bin.size)

    #2.3 NSF observing
    FL_Sandbar_NSF$FINYEAR=as.character(with(FL_Sandbar_NSF,
                                             ifelse(Month%in%1:6,paste(year-1,"-",fn.subs(year),sep=""),
                                                    ifelse(Month%in%7:12,paste(year,"-",fn.subs(year+1),sep=""),NA))))
    if(Conv.cal.mn.to.fin.mn=="YES") FL_Sandbar_NSF$Month=fn.month.to.fin.mn(FL_Sandbar_NSF$Month)

    FL_Sandbar_NSF$FL=round(FL_to_TL(FL_Sandbar_NSF$FL,a.TL,b.TL))   #convert FL to TL

    FL_NSF_observers=fn.add.missing.size2(FL_Sandbar_NSF,Min.FL,round(FL_to_TL(Mx.FL,a.TL,b.TL)),Bin.size)
  }

  if(SP=="BW")
  {
    FL_dusky_NSF$FINYEAR=as.character(with(FL_dusky_NSF,
                                           ifelse(Month%in%1:6,paste(year-1,"-",fn.subs(year),sep=""),
                                                  ifelse(Month%in%7:12,paste(year,"-",fn.subs(year+1),sep=""),NA))))
    if(Conv.cal.mn.to.fin.mn=="YES") FL_dusky_NSF$Month=fn.month.to.fin.mn(FL_dusky_NSF$Month)

    FL_dusky_NSF$FL=round(FL_to_TL(FL_dusky_NSF$FL,a.TL,b.TL))   #convert FL to TL

    FL_NSF_observers=fn.add.missing.size2(FL_dusky_NSF,Min.FL,round(FL_to_TL(Mx.FL,a.TL,b.TL)),Bin.size)


    #2.4 Wetline_WRL
    Size.comp.Dusky.WRL=data.frame(FINYEAR="1999-00",TL=Dusky.WRL$TL)
    Size.comp.Dusky.WRL$FINYEAR=as.character(Size.comp.Dusky.WRL$FINYEAR)
    Size.comp.Dusky.WRL$FL=round(TL_to_FL(Dusky.WRL$TL,a.TL,b.TL))
    Size.comp.Dusky.WRL$Number=1
    Size.comp.Dusky.WRL$Year=1
    Size.comp.Dusky.WRL$Type="Fisher_measure"

  }


  #2.5 Put all size composition (as TL) data in list
  #note: don't use ..._FL_Heald_1987, very unreliable and strange size measures and unsure about how measured
  #     don't use ..._FL_Stevens_1990, as they cannot be allocated to a zone and unsure about how measured
  All.size=list(WC=TL_observers.WC,Zn1=TL_observers.Zn1,Zn2=TL_observers.Zn2)
  All.size_7=list(WC=TL_observers.WC_7,Zn1=TL_observers.Zn1_7,Zn2=TL_observers.Zn2_7)


  #drop non-representative observations in TDGDLF (other fisheries are ok)
  drop.obs=function(dd,ZN)
  {
    dd$dummy=paste(dd$FINYEAR,ZN)
    dd=subset(dd,dummy%in%This.yr.zn)
    dd=dd[,-match('dummy',names(dd))]
    return(dd)
  }
  zns=c("West","Zone1","Zone2")
  for(i in 1:length(All.size))
  {
    All.size[[i]]=drop.obs(All.size[[i]],zns[i])
    All.size_7[[i]]=drop.obs(All.size_7[[i]],zns[i])
  }


  if(SP%in%c("WH","GM")) All.size=list(TDGDLF=All.size,TDGDLF_7=All.size_7)
  if(SP=="BW")  All.size=list(TDGDLF=All.size,TDGDLF_7=All.size_7,NSF=FL_NSF_observers)
  if(SP=="TK") All.size=list(TDGDLF=All.size,TDGDLF_7=All.size_7,Pilbara_trawl=Sandbar_FL_Pil.trwl_observers,NSF=FL_NSF_observers)


  #2.6 Put abundance data together

  #TDGDLF
  if(SP %in% c("WH","BW"))
  {
    #monthly
    if(exists("Ab.indx.TDGDLF.West"))Ab.indx.TDGDLF.West$zone="West"
    Ab.indx.TDGDLF.Zn1$zone="Zone1"
    Ab.indx.TDGDLF.Zn2$zone="Zone2"
    if(exists("Ab.indx.TDGDLF.West"))Ab.indx.TDGDLF=rbind(Ab.indx.TDGDLF.West,Ab.indx.TDGDLF.Zn1,Ab.indx.TDGDLF.Zn2)
    if(!exists("Ab.indx.TDGDLF.West"))Ab.indx.TDGDLF=rbind(Ab.indx.TDGDLF.Zn1,Ab.indx.TDGDLF.Zn2)

    #daily
    Ab.indx.TDGDLF.West.daily$zone="West"
    Ab.indx.TDGDLF.Zn1.daily$zone="Zone1"
    Ab.indx.TDGDLF.Zn2.daily$zone="Zone2"
    Ab.indx.TDGDLF.daily=rbind(Ab.indx.TDGDLF.West.daily,Ab.indx.TDGDLF.Zn1.daily,Ab.indx.TDGDLF.Zn2.daily)
  }

  if(SP=="GM")
  {
    #monthly
    Ab.indx.TDGDLF.Zn2$zone="Zone2"
    Ab.indx.TDGDLF=Ab.indx.TDGDLF.Zn2

    #daily
    Ab.indx.TDGDLF.Zn2.daily$zone="Zone2"
    Ab.indx.TDGDLF.daily=Ab.indx.TDGDLF.Zn2.daily
  }

  if(SP=="TK")
  {
    #monthly
    Ab.indx.TDGDLF.West$zone="West"
    Ab.indx.TDGDLF.Zn1$zone="Zone1"
    Ab.indx.TDGDLF=rbind(Ab.indx.TDGDLF.West,Ab.indx.TDGDLF.Zn1)

    #daily
    Ab.indx.TDGDLF.West.daily$zone="West"
    Ab.indx.TDGDLF.Zn1.daily$zone="Zone1"
    Ab.indx.TDGDLF.daily=rbind(Ab.indx.TDGDLF.West.daily,Ab.indx.TDGDLF.Zn1.daily)
  }


  #Naturaliste survey
  if(exists('Ab.index.Srvy.FixSt'))
  {
    Ab.index.Srvy.FixSt$FINYEAR=paste(Ab.index.Srvy.FixSt$FINYEAR,"-",fn.subs(Ab.index.Srvy.FixSt$FINYEAR+1),sep="")
    Naturaliste.abun=Ab.index.Srvy.FixSt
  }
  if(exists('Size.index.Srvy.FixSt'))
  {
    Size.index.Srvy.FixSt$FINYEAR=paste(Size.index.Srvy.FixSt$FINYEAR,"-",fn.subs(Size.index.Srvy.FixSt$FINYEAR+1),sep="")
    Naturaliste.size=Size.index.Srvy.FixSt
  }


  #2.7 Years with data for conventional tagging data
  a=Zn.rel_Conv.Tag
  a$Finyr=with(a,ifelse(Mn.rel>6,Yr.rel+1,Yr.rel))
  c.Tag=list(TDGDLF=data.frame(FINYEAR=paste(sort(unique(a$Finyr)),"-",fn.subs(sort(unique(a$Finyr))+1),sep="")))

  #Add financial year to release and recapture files
  fn.finyr.rel=function(dat)with(dat,ifelse(Mn.rel>6,Yr.rel,Yr.rel-1))
  fn.finyr.rec=function(dat)with(dat,ifelse(Mn.rec>6,Yr.rec,Yr.rec-1))
  fn.finyr=function(Finyr)paste(Finyr,"-",fn.subs(Finyr+1),sep="")

  Zn.rel_Conv.Tag$FinYear.rel=fn.finyr.rel(Zn.rel_Conv.Tag)
  Zn.rec_Conv.Tag$FinYear.rec=fn.finyr.rec(Zn.rec_Conv.Tag)
  Zn.rel_Conv.Tag$FinYear.rel=fn.finyr(Zn.rel_Conv.Tag$FinYear.rel)
  Zn.rec_Conv.Tag$FinYear.rec=fn.finyr(Zn.rec_Conv.Tag$FinYear.rec)

  #Convert month to financial month
  if(Conv.cal.mn.to.fin.mn=="YES")
  {
    Zn.rel_Conv.Tag$Mn.rel=fn.month.to.fin.mn(Zn.rel_Conv.Tag$Mn.rel)
    Zn.rec_Conv.Tag$Mn.rec=fn.month.to.fin.mn(Zn.rec_Conv.Tag$Mn.rec)
  }

  #Recode TG.zn
  fn.recode=function(dat,dat1)
  {
    dat=aggregate(Number~TG.zn+Rel.zone+Yr.rel+Mn.rel+Age+FinYear.rel,dat,sum)
    dat1=aggregate(Number~TG.zn+Rec.zone+Yr.rec+Mn.rec+FinYear.rec,dat1,sum)

    dat=dat[order(dat$FinYear.rel,dat$Mn.rel),]
    dat$TG=1:nrow(dat)
    dat1=merge(dat1,subset(dat,select=c(TG.zn,TG)),by="TG.zn",all.x=T)
    return(list(dat=dat,dat1=dat1))
  }
  a=fn.recode(Zn.rel_Conv.Tag,Zn.rec_Conv.Tag)
  Zn.rel_Conv.Tag=a$dat
  Zn.rec_Conv.Tag=a$dat1


  #2.8 Years with data for acoustic tagging data

  #Add financial year
  fn.finyr.rel=function(dat)with(dat,ifelse(Month.rel>6,Year.rel,Year.rel-1))
  fn.finyr.rec=function(dat)with(dat,ifelse(Month>6,Year,Year-1))

  Zn.rel_Acous.Tag$FinYear.rel=fn.finyr.rel(Zn.rel_Acous.Tag)
  Zn.rec_Acous.Tag$FinYear.rec=fn.finyr.rec(Zn.rec_Acous.Tag)
  Zn.rel_Acous.Tag$FinYear.rel=fn.finyr(Zn.rel_Acous.Tag$FinYear.rel)
  Zn.rec_Acous.Tag$FinYear.rec=fn.finyr(Zn.rec_Acous.Tag$FinYear.rec)

  a.Tag=sort(unique(Zn.rec_Acous.Tag$FinYear.rec))

  #Convert month to financial month
  if(Conv.cal.mn.to.fin.mn=="YES")
  {
    Zn.rel_Acous.Tag$Month.rel=fn.month.to.fin.mn(Zn.rel_Acous.Tag$Month.rel)
    Zn.rec_Acous.Tag$Month=fn.month.to.fin.mn(Zn.rec_Acous.Tag$Month)
  }

  #Create table of proportion of recaptures by year
  fn.plot.prop.rec.yr=function(dat)
  {
    a=aggregate(Number~FinYear.rec,dat,sum)
    a$Prop=a$Number/sum(a$Number)
    Yr1=as.numeric(substr(min(dat$FinYear.rec),start=1,stop=4))
    Yr1=1993
    Yr2=as.numeric(substr(max(dat$FinYear.rec),start=1,stop=4))
    SEQ=seq(Yr1,Yr2)
    SEQ1=seq(Yr1+1,Yr2+1)
    All.yrs=paste(SEQ,"-",substr(SEQ1,3,4),sep="")
    id=All.yrs[which(!All.yrs%in%a$FinYear.rec)]
    Add=data.frame(FinYear.rec=id,Number=0,Prop=0)
    a=rbind(a,Add)
    a=a[order(a$FinYear.rec),]
    plot(1:length(a$FinYear.rec),a$Prop,xaxt='n',ylab="",xlab="",cex.axis=1.25,
         pch=19,col=2,cex=2)
    axis(1,1:length(a$FinYear.rec),F,tck=-0.03)
    axis(1,seq(1,length(a$FinYear.rec),2),a$FinYear.rec[seq(1,length(a$FinYear.rec),2)],cex.axis=.9,tck=-0.06)
  }


  #set working directory for outputing data
  HandL="C:/Matias/Analyses/Population dynamics/1."
  DiR=paste(HandL,SP.LABELS,"/",Yr.assess,"/1_Inputs/Visualise data",sep='')
  if(!file.exists(DiR))
  {
    mainDir=paste(HandL,SP.LABELS,"/",Yr.assess,sep="")
    dir.create(mainDir)
    subDir="1_Inputs"
    dir.create(file.path(mainDir,subDir))
    subDir="/1_Inputs/Visualise data"
    dir.create(file.path(mainDir,subDir))
  }
  setwd(DiR)


  #Effective sample size
  Dat.eff.n=rbind(FL.TDGDFL.WC,FL.TDGDFL.Zn1,FL.TDGDFL.Zn2)
  Eff.n=list(Female=NA,Male=NA)
  fn.get.eff=function(SeX)
  {
    aa=subset(Dat.eff.n,SEX==SeX)
    Steps=seq(2,length(aa$TL),1)
    fn.rand.samp.size=function(STEP) sample(aa$TL,STEP,replace=T)
    Store.samps.size=vector('list',length(Steps))
    for (i in 1:length(Steps)) Store.samps.size[[i]]=fn.rand.samp.size(Steps[i])
    Dummy=data.frame(mean=Steps,sd=Steps,se=Steps)
    for (i in 1:length(Steps))
    {
      Dummy$mean[i]=mean(Store.samps.size[[i]])
      Dummy$sd[i]=sd(Store.samps.size[[i]])
      Dummy$se[i]=Dummy$sd[i]/sqrt(Steps[i])
    }
    return(Dummy)
  }
  Eff.n$Female=fn.get.eff("F")
  Eff.n$Male=fn.get.eff("M")


  fn.plot.eff.n=function(Dummy,CL)
  {
    plot(1:nrow(Dummy),Dummy$mean,ylim=c(min(Dummy$mean)*1,max(Dummy$mean)*1),ylab="",xlab="",
         col=CL,pch=19,cex.axis=1.15)
    if(x==1) mtext("Mean TL (cm)",2,las=3,cex=1.3,line=2.2)
    plot(1:nrow(Dummy),Dummy$sd,ylim=c(min(Dummy$sd)*1,max(Dummy$sd)*1),ylab="",xlab="",
         col=CL,pch=19,cex.axis=1.15)
    if(x==1) mtext("SD (cm)",2,las=3,cex=1.3,line=2.2)
    #Dummy[match(300,Steps),1]/Dummy[length(Steps),1]
  }
  clss=c("pink","blue")
  fn.fig("Effective sample size",2400, 1800)
  par(mfcol=c(2,2),mai=c(.6,.6,.1,.1),las=1,mgp=c(2,.6,0))
  for(x in 1:2)fn.plot.eff.n(Eff.n[[x]],clss[x])
  mtext("Sample size",1,outer=T,cex=1.3,line=-1.5)
  dev.off()

  fn.fig("Proportion.tag.recaptures",2400, 1800)
  par(mfcol=c(1,1),las=1,mai=c(0.3,0.55,.1,.1),oma=c(2.25,2.25,.1,.1),mgp=c(1,.5,0))
  fn.plot.prop.rec.yr(Zn.rec_Conv.Tag)
  legend("topright",SP.LABELS,bty='n',cex=1.5)
  mtext("Financial year",1,line=1,cex=1.5,outer=T)
  mtext("Proportion of recaptures",2,line=1,cex=1.5,outer=T,las=3)
  dev.off()


  #2.9 Add missing years to prop.eff.mesh
  #note: back fill using temporal change in 7inch upt o 2009-10
  dum.yr=sort(as.character(unique(Data.monthly$FINYEAR)))
  id=dum.yr[which(!dum.yr%in%Mesh.prop.eff$finyear)]
  add.msh.eff=Mesh.prop.eff[(nrow(Mesh.prop.eff)+1):((nrow(Mesh.prop.eff))+length(id)),]
  add.msh.eff$finyear=id

  Mesh.prop.eff=rbind(add.msh.eff,Mesh.prop.eff)
  Mesh.prop.eff.West=rbind(add.msh.eff,Mesh.prop.eff.West)
  Mesh.prop.eff.Zn1=rbind(add.msh.eff,Mesh.prop.eff.Zn1)
  Mesh.prop.eff.Zn2=rbind(add.msh.eff,Mesh.prop.eff.Zn2)

  bck.fil.yrs=match(c("2005-06","2006-07","2007-08",
                      "2008-09","2009-10"),Mesh.prop.eff$finyear)
  fn.bck.fil=function(a,intercpt)
  {
    a$yr=1:nrow(a)
    a$OFFSET=intercpt
    DaTt=a[bck.fil.yrs,]
    #set intercept
    Model=lm(X165~0+yr,offset=OFFSET,DaTt)
    NEWDAT=data.frame(yr=a[1:(bck.fil.yrs[1]-1),]$yr,
                      OFFSET=unique(a$OFFSET))
    P_6.5=predict(Model,newdata = NEWDAT)
    P_7=1-P_6.5
    return(data.frame(X165=P_6.5,X178=P_7))
  }
  Mesh.prop.eff[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff,intercpt=0)
  Mesh.prop.eff.West[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff.West,intercpt=0)
  Mesh.prop.eff.Zn1[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff.Zn1,intercpt=1)
  Mesh.prop.eff.Zn2[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff.Zn2,intercpt=0)


  colnames(Mesh.prop.eff)[2:3]=colnames(Mesh.prop.eff.West)[2:3]=
    colnames(Mesh.prop.eff.Zn1)[2:3]=colnames(Mesh.prop.eff.Zn2)[2:3]=c("Mesh_6.5","Mesh_7")


  #--------------- RESULTS SECTION --------------

  #Zone colors
  Zns=c("Joint","North","Closed","West","Closed.metro","Zone1","Zone2")
  Zns.leg=c("JANSF","WANCSF","Ningaloo","WCDGDLF","Metro closure","Zone1","Zone2")
  COL.prop=c("aquamarine2","lightseagreen","seagreen4","lightgreen","olivedrab4","olivedrab3","mediumseagreen")
  names(COL.prop)=Zns

  #Biological regions colors
  Bio.col=c("dodgerblue","darkorchid4","cyan4","lightpink3")
  names(Bio.col)=c("North Coast","Gascoyne","West Coast","South Coast")


  #1. Visualize catch
  LTY=rep(1,4)

  fn.add=function(a)
  {
    id=which(!YRS%in%a$FINYEAR)
    if(length(id)>0)
    {
      ADD=a
      ADD=ADD[1:length(id),]
      ADD$FINYEAR=YRS[id]
      xx=match("FINYEAR",colnames(ADD))
      ADD[,which(!1:ncol(ADD)%in%xx)]=NA
      a=rbind(a,ADD)
      a=a[order(a$FINYEAR),]
    }
    return(a)
  }

  fn.see.Ktch=function(DAT,YRS,NAMES,where,cx.zn,cx.Titl,cx.Man,cx.axs,tck1,tck2,tck3)
  {
    axis1=function()axis(1,1:length(YRS),F,tck=tck1)
    if(!SP=="TK")THIS.yr=which(YRS=="1980-81")
    if(SP=="TK")THIS.yr=which(YRS=="1985-86")
    axis11=function()axis(1,seq(THIS.yr,length(YRS),10),YRS[seq(THIS.yr,length(YRS),10)],
                          tck=tck2,cex.axis=cx.axs)
    if(!SP=="TK") THIS.yr=which(YRS=="1975-76")
    axis12=function()axis(1,seq(THIS.yr,length(YRS),10),F,tck=tck3)

    Single=c("GAB_trawl","WRL","WTBF","TEPS","SA_marine")
    Byregion=c("TDGDLF","NSF","Rec","Other")
    for(i in 1:length(DAT))
    {
      a=DAT[[i]]
      if(names(DAT)[i]%in% Single) a=aggregate(LIVEWT.c~FINYEAR,a,sum)

      if(names(DAT)[i]%in%Byregion)a=aggregate(LIVEWT.c~FINYEAR+zone,a,sum)

      if(is.numeric(a$FINYEAR[1])) a$FINYEAR=paste(a$FINYEAR,"-",fn.subs(a$FINYEAR+1),sep="")
      a=a[order(a$FINYEAR),]
      if(!is.logical(a))
      {
        a$FINYEAR=factor(a$FINYEAR,levels=YRS)
        if(names(DAT)[i]%in% Single)
        {
          nm=NAMES[i]
          if(NAMES[i]=="SA_marine") nm="SA_marine (whalers catch)"
          a=fn.add(a)
          plot(1:length(YRS),a[,2],lwd=2,xaxt='n',ylab="",xlab="",type='l',main=nm)
          axis1()
          axis11()
          if(names(DAT[i])=="Historic.TDGDLF")
          {
            IDs=c(4:8,11)
            points(IDs,a[IDs,2],col="orange",cex=1.25,pch=19)
            legend("topright","interpolated",col="orange",pch=19,cex=1.25,pt.cex=2,bty='n')
          }
        }

        if(names(DAT)[i]%in%Byregion)
        {
          Tot=aggregate(LIVEWT.c~FINYEAR,a,sum)

          plot(1:length(YRS),xaxt='n',ylim=c(0,max(Tot$LIVEWT.c,na.rm=T)),
               col='transparent',ylab="",xlab="",main=NAMES[i],cex.main=cx.Man,cex.axis=cx.axs)

          zn=unique(a$zone)
          for(x in 1:length(zn))
          {
            b=subset(a,zone==zn[x])
            b=fn.add(b)
            b=b[order(as.character(b$FINYEAR)),]
            lines(1:length(YRS),b[1:length(YRS),ncol(b)],lwd=3,lty=LTY[x],col=COL[[i]][match(zn[x],names(COL[[i]]))])
          }

          #add total
          Tot$FINYEAR=as.character(Tot$FINYEAR)
          id=which(!YRS%in%Tot$FINYEAR)
          if(length(id)>0)
          {
            ddd=Tot[(nrow(Tot)+1):((nrow(Tot))+length(id)),]
            ddd$FINYEAR=YRS[id]
            Tot=rbind(Tot,ddd)
            Tot=Tot[order(Tot$FINYEAR),]
          }
          lines(1:length(YRS),Tot$LIVEWT.c,lwd=3,lty=1,col=tot.col)

          legend(where[i],c(paste(zn),"Total"),lty=LTY,col=c(COL[[i]][match(zn,names(COL[[i]]))],tot.col),bty='n',cex=cx.zn,lwd=3)
        }

        axis1()
        axis11()
        axis12()
      }

    }

  }

  YR.span=sort(as.character(unique(catch$TDGDLF$FINYEAR)))
  if(SP=="TK")YR.span=unique(sort(c("1981-82","1982-83",as.character(unique(catch$NSF$FINYEAR)),YR.span)))
  YRS=YR.span
  line.x=0.75; line.y=2

  if(SP=="WH")
  {
    PAR=function()par(mfcol=c(3,1),mai=c(.25,.1,.3,.2),oma=c(2,4,.1,1),mgp=c(1,.75,0),las=1)
    Fisheries=c("TDGDLF","Other WA fisheries","Recreational")
    File="Whiskery.catch"
    where=c('topright','topright','topleft')
    Zn.col=COL.prop[match(c("West","Zone1","Zone2"),names(COL.prop))]
    Bio.Cl=Bio.col[match(c("West Coast","South Coast"),names(Bio.col))]
    COL=list(Zn.col,Zn.col,Bio.Cl)
    CX.zn=1.75;CX.Titl=1.75;CX.Man=1.9;CX.axs=1.5;Tck1=(-0.02);Tck2=(-0.04);Tck3=(-0.03)
  }

  if(SP=="GM")
  {
    PAR=function()par(mfcol=c(3,1),mai=c(.25,.1,.3,.2),oma=c(2,4,.1,1),mgp=c(1,.75,0),las=1)
    Fisheries=c("TDGDLF","Other WA fisheries","Recreational")
    File="Gummy.catch"
    where=c('topleft','topright','topleft')
    Zn.col=COL.prop[match(c("West","Zone1","Zone2"),names(COL.prop))]
    Bio.Cl=Bio.col
    COL=list(Zn.col,Zn.col,Bio.Cl)

    CX.zn=1.75;CX.Titl=1.75;CX.Man=1.9;CX.axs=1.5;Tck1=(-0.02);Tck2=(-0.04);Tck3=(-0.03)
  }

  if(SP=="BW")
  {
    if(Ktch.source=="ALL") PAR=function()par(mfcol=c(4,2),mai=c(.1,.1,.3,.15),oma=c(4,4,.3,.75),mgp=c(2,.5,0),las=1)
    if(Ktch.source=="WA.only") PAR=function()par(mfcol=c(3,2),mai=c(.1,.1,.25,.15),oma=c(4,4,1,1),mgp=c(1,.5,0),las=1)

    Fisheries=names(catch)
    File="Dusky.catch"
    where=c('topright','topleft','topright',"top",'topleft')
    Zn.col=COL.prop[match(c("West","Zone1","Zone2"),names(COL.prop))]
    Zn.colN=c(COL.prop[match(c("Closed.ningaloo","North"),names(COL.prop))],"lightcyan3")
    Bio.Cl=Bio.col
    COL=list(Zn.col,Zn.colN,Zn.col,'black',Bio.Cl)

    CX.zn=1.2;CX.Titl=1.5;CX.Man=1.25;CX.axs=1.25;Tck1=(-0.01);Tck2=(-0.03);Tck3=(-0.02)
    line.x=1.5; line.y=2
  }

  if(SP=="TK")
  {
    if(Ktch.source=="ALL")  PAR=function()par(mfcol=c(3,2),mai=c(.1,.1,.2,.2),oma=c(4,4,1,1),mgp=c(1,.5,0),las=1)
    if(Ktch.source=="WA.only")  PAR=function()par(mfcol=c(2,2),mai=c(.1,.1,.5,.3),oma=c(3,4,1.25,1),mgp=c(1,.5,0),las=1)
    Fisheries=names(catch)
    line.x=1.75
    File="Sandbar.catch"
    where=c('topleft','topleft','topleft',"topleft")
    Zn.col=COL.prop[match(c("West","Zone1","Zone2"),names(COL.prop))]
    Zn.colN=COL.prop[match(c("Closed","North","Joint"),names(COL.prop))]
    Bio.Cl=Bio.col
    COL=list(Zn.col,Zn.colN,Zn.col,Bio.Cl)

    CX.zn=1.15;CX.Titl=1.5;CX.Man=1.25;CX.axs=1.25;Tck1=(-0.01);Tck2=(-0.03);Tck3=(-0.02)
  }

  if(SP%in%c("GM","WH")) WIDTH=LENGTH=2400
  if(SP%in%c("BW","TK"))
  {
    WIDTH=2400
    LENGTH=1800
  }
  tot.col="tan3"
  fn.fig(File,WIDTH, LENGTH)
  PAR()
  fn.see.Ktch(catch,YR.span,Fisheries,where,CX.zn,CX.Titl,CX.Man,CX.axs,Tck1,Tck2,Tck3)
  mtext("Financial year",1,line=line.x,cex=CX.Titl,outer=T)
  mtext("Total catch (tonnes)",2,line=line.y,cex=CX.Titl,outer=T,las=3)
  dev.off()

  #Compare historic and catch time series
  Ag.TDGDLF=aggregate((LIVEWT.c/1000)~FINYEAR,catch$TDGDLF,sum)
  names(Ag.TDGDLF)[2]="LIVEWT.c"
  Ag.TDGDLF$year=as.numeric(substr(Ag.TDGDLF$FINYEAR,1,4))
  fn.fig("Catch_Historic & time series",2400, 2400)
  par(las=1,cex.axis=1.25,mgp=c(2.8,.7,0),xpd=T)
  plot(Ag.TDGDLF$year,Ag.TDGDLF$LIVEWT.c,xlim=c(min(historic.prop.ktch$year),max(Ag.TDGDLF$year)),
       ylim=c(0,max(c(subset(Historic.ktch,year<1975)$LIVEWT.c,Ag.TDGDLF$LIVEWT.c))),pch=19,cex=1.5,ylab="Catch (tonnes)",xlab="Financial year",
       cex.lab=1.5)
  points(historic.prop.ktch$year,historic.prop.ktch$LIVEWT.c,pch=21,bg='orange',col='orange',cex=1.5)
  with(subset(Historic.ktch,year<1975),points(year,LIVEWT.c,col="orange",cex=1.5))
  legend("topleft",c("historic shark landings (Whitley (1944), Heald (1987), Simpfendorfer & Donohue 1998)",
                     paste("reconstructed historic",Spec,"shark landings (mean prop=",round(Mean.prop.ctch,2)," in the first",yrs.considered," years)"),
                     "model inputs"),bty="n",pch=c(21,21,21),col=c('orange','orange',"black"),
         pt.bg=c('white','orange',"black"),cex=1,inset=c(0,-0.125))
  dev.off()


  #2. Visualize size composition
  Min.yr=min(c(FL.TDGDFL.WC$year,FL.TDGDFL.Zn1$year,FL.TDGDFL.Zn2$year,
               FL.TDGDFL.WC_7$year,FL.TDGDFL.Zn1_7$year,FL.TDGDFL.Zn2_7$year))
  Mx.yr=max(c(FL.TDGDFL.WC$year,FL.TDGDFL.Zn1$year,FL.TDGDFL.Zn2$year,
              FL.TDGDFL.WC_7$year,FL.TDGDFL.Zn1_7$year,FL.TDGDFL.Zn2_7$year))
  FinYrs=paste(Min.yr:(Mx.yr-1),"-",fn.subs((Min.yr+1):Mx.yr),sep="")

  fn.bub=function(DAT,COL,Scale,bin,N)
  {
    if(nrow(DAT)==0) plot(1, type="n", axes=F, xlab="", ylab="")

    if(nrow(DAT)>0)
    {
      #aggregate data into years
      YRs=as.character(unique(DAT$FINYEAR))
      SIzes=colnames(DAT[,4:ncol(DAT)])
      d=matrix(ncol=length(SIzes),nrow=length(YRs))
      for(s in 1:length(YRs))
      {
        a=subset(DAT,FINYEAR==YRs[s])
        d[s,]=colSums(a[,4:ncol(a)])
      }
      colnames(d)=SIzes
      d=cbind(FINYEAR=YRs,as.data.frame.matrix(d))

      #keep years with minimum observations
      drop.yr=rowSums(d[,2:ncol(d)])
      names(drop.yr)=d$FINYEAR
      drop.yr=subset(drop.yr,drop.yr<Min.obs)
      d=subset(d,!FINYEAR%in%names(drop.yr))

      DAT=d

      Ylabs=FinYrs
      ID=which(!FinYrs%in%DAT$FINYEAR)
      ADD=data.frame(FINYEAR=FinYrs[ID])
      if(nrow(ADD)>0)
      {
        ADD1=DAT[1:nrow(ADD),2:ncol(DAT)]
        ADD=cbind(ADD,ADD1)
        ADD[,2:ncol(ADD)]=0
        DAT=rbind(DAT,ADD)
      }

      DAT$FINYEAR=as.character(DAT$FINYEAR)
      DAT=DAT[order(DAT$FINYEAR),]
      DAT=DAT[,2:ncol(DAT)]
      x=as.numeric(colnames(DAT))
      y=1:nrow(DAT)
      z=(as.matrix(DAT))
      n=length(y)
      xo=outer(x,rep(1,length=length(y)))
      yo=t(outer(y,rep(1,length=length(x))))
      zo=z/ rowSums(z)*Scale
      matplot(yo,xo,type="n",xlab="",ylab="",xaxt='n',yaxt='n')
      abline(v=pretty(x),lty=3,col="black")
      for(i in 1:n) points(yo[,i],xo[,i],cex=zo[i,],pch=16,col=COL)
      axis(1,y,F,tck=-0.025)
      axis(1,y[seq(1,length(Ylabs),5)],Ylabs[seq(1,length(Ylabs),5)],tck=-0.05,cex.axis=1.25)
      axis(2,x,F,tck=-0.025)
      axis(2,seq(x[1],x[length(x)],by=N),seq(x[1],x[length(x)],by=N),tck=-0.05,cex.axis=CxY)
    }
  }

  if(SP%in%c("WH","GM")) CxY=1.5
  if(SP%in%c("BW","TK")) CxY=0.9
  #COLS="steelblue"
  COLS="grey60"
  BIN=Bin.size #bin size
  Nl=length(All.size$TDGDLF)
  LEGS=c("West coast","Zone 1","Zone 2")
  if(SP%in%c("WH","GM")) WhereLEGN="bottomright"; NN=10
  if(SP%in%c("BW","TK")) WhereLEGN="topright"; NN=20

  #2.1 TDGDLF size comp of reported catch
  fn.fig("Size.comp.TDGDLF",2400, 2000)
  par(mfrow=c(3,2),las=1,mai=c(0.3,0.35,.2,.1),oma=c(2.25,2.25,.1,.1),mgp=c(1,.85,0))
  for(e in 1:Nl)
  {
    fn.bub(All.size$TDGDLF[[e]],COLS,Scale=7.5,BIN,N=NN)
    if(e==1)mtext("6.5 inch mesh",3,0,cex=1.25)
    fn.bub(All.size$TDGDLF_[[e]],COLS,Scale=7.5,BIN,N=NN)
    if(e==1)mtext("7 inch mesh",3,0,cex=1.25)
    legend(WhereLEGN,LEGS[e],bty='n',cex=1.85)
  }
  mtext("Financial Year",1,cex=1.5,line=0.75,outer=T)
  mtext("Total length class (cm)",2,cex=1.5,line=0.6,las=3,outer=T)
  dev.off()

  #2.2 NSF
  if(SP%in%c('BW',"TK"))
  {
    fn.fig("Size.comp.NSF",2000, 1800)
    par(mfcol=c(1,1),las=1,mai=c(.8,1.2,.1,.1),mgp=c(1,1.5,0))
    fn.bub(FL_NSF_observers,COLS,Scale=10,BIN,N=NN)
    mtext("Financial Year",1,cex=1.25,line=3)
    mtext("Total length class (cm)",2,cex=1.25,line=4,las=3)
    dev.off()
  }

  #2.4 Pilbara Trawl
  if(SP=="TK")
  {
    FinYrs=unique(Sandbar_FL_Pil.trwl_observers$FINYEAR)

    fn.fig("Size.comp.Pilbara_trawl",2000, 1800)
    par(las=1,mai=c(0.75,0.85,.1,.1),mgp=c(1,1.25,0))
    fn.bub(Sandbar_FL_Pil.trwl_observers,COLS,Scale=10,BIN,N=NN)
    mtext("Financial Year",1,cex=1.25,line=2.75)
    mtext("Total length class (cm)",2,cex=1.25,line=2.75,las=3)
    dev.off()
  }


  #2.5 TEPS_TDGDLF
  if(SP=="BW")
  {
    fn.fig("Size.comp.Oversized_TDGDLF",2000,1800)
    hist(Size.comp.Dusky.TEPS_TDGLDF,xlab="FL (m)",main="Oversized dusky_TDGDLF")
    dev.off()

  }

  #2.6 Table of number of observations and shots
  fn.table.shots=function(dat)
  {
    a=dat
    if(!is.na(match('CALCULATED FL',names(a))))
    {
      names(a)[match('CALCULATED FL',names(a))]='CALCULATED.FL'
      a$FL=with(a,ifelse(is.na(FL),CALCULATED.FL,FL))
    }
    a$FINYEAR=with(a,ifelse(Month>6,paste(year,"-",fn.subs(year+1),sep=""),
                            paste(year-1,"-",fn.subs(year),sep="")))
    a=subset(a,!is.na(FL))
    a$Number=1
    Obs=aggregate(Number~FINYEAR,a,sum)
    a$Dup=paste(a$year,a$Month,a$SHEET_NO)
    bb=a[!duplicated(a$Dup),]
    bb$Number=1
    Shots=aggregate(Number~FINYEAR,bb,sum)
    this=merge(Obs,Shots,by="FINYEAR")
    names(this)[2:3]=c("N.observations","N.shots")
    this$Species=unique(a$SPECIES)
    this$Fishery="NSF"
    this$zone=unique(dat$zone)

    return(this)
  }
  if(SP%in%c("BW","TK"))
  {
    Sandbar.NSF.size.numbers=fn.table.shots(FL_Sandbar_NSF)
    Dusky.NSF.size.numbers=fn.table.shots(FL_dusky_NSF)
  }
  if(SP%in%c("WH","GM")) Size.numbers=TDGDFL.size.numbers
  if(SP=="BW") Size.numbers=rbind(TDGDFL.size.numbers,Dusky.NSF.size.numbers)
  if(SP=="TK") Size.numbers=rbind(TDGDFL.size.numbers,Sandbar.NSF.size.numbers)

  Size.numbers=Size.numbers[,match(c("Fishery","zone","FINYEAR",
                                     "N.observations","N.shots"),names(Size.numbers))]
  Size.numbers=Size.numbers[order(Size.numbers$FINYEAR,Size.numbers$Fishery,Size.numbers$zone),]

  Numbers.SF=reshape(Size.numbers[,-match("N.shots",names(Size.numbers))],v.names = "N.observations",
                     idvar = c("Fishery","zone"),timevar = "FINYEAR", direction = "wide")
  X=colnames(Numbers.SF)[3:ncol(Numbers.SF)]
  colnames(Numbers.SF)[3:ncol(Numbers.SF)]=as.character(sapply(strsplit(X,"N.observations."), "[", 2))
  Shots.SF=reshape(Size.numbers[,-match("N.observations",names(Size.numbers))],v.names = "N.shots",
                   idvar = c("Fishery","zone"),timevar = "FINYEAR", direction = "wide")
  X=colnames(Shots.SF)[3:ncol(Shots.SF)]
  colnames(Shots.SF)[3:ncol(Shots.SF)]=as.character(sapply(strsplit(X,"N.shots."), "[", 2))
  Numbers.SF=Numbers.SF[order(Numbers.SF$Fishery,Numbers.SF$zone),]
  Shots.SF=Shots.SF[order(Shots.SF$Fishery,Shots.SF$zone),]

  #add total
  dum= Numbers.SF[1,]
  dum$zone="Total"
  dum[,-match(c("Fishery","zone"),names(Numbers.SF))]=colSums(Numbers.SF[,-match(c("Fishery","zone"),names(Numbers.SF))],na.rm =T)
  Numbers.SF=rbind(Numbers.SF,dum)
  dum[,-match(c("Fishery","zone"),names(Numbers.SF))]=colSums(Shots.SF[,-match(c("Fishery","zone"),names(Shots.SF))],na.rm =T)
  Shots.SF=rbind(Shots.SF,dum)

  #remove NAs
  Numbers.SF[is.na(Numbers.SF)]=""
  Shots.SF[is.na(Shots.SF)]=""

  #create nice table
  fn.word.table(WD=getwd(),TBL=Numbers.SF,Doc.nm="Size.comp.n.observations",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")

  fn.word.table(WD=getwd(),TBL=Shots.SF,Doc.nm="Size.comp.n.shots",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")

  write.csv(Numbers.SF,"Numbers.SF.csv",row.names=F)
  write.csv(Shots.SF,"Shots.SF.csv",row.names=F)


  #3. Visualize mean weights
  if(exists("Avr.wt.yr.zn"))
  {
    fn.see.avg.wgt=function()
    {
      b=Avr.wt.yr.zn
      zn=unique(b$zone)
      N=1:length(unique(b$Finyear))
      SD=b$mean*(b$CV)
      plot(N,ylim=c(0,max(b$mean+SD)*1.05),main="",cex.main=1.25,xaxt='n',
           ylab="",xlab="",pch=19,cex=3,cex.axis=1,col="transparent",xlim=c(0,N[length(N)]+0.5))

      jit=c(0,.1,.2)
      CLOS=Zn.col

      for(x in 1:length(zn))
      {
        a=subset(b, zone==zn[x])
        N1=N+jit[x]
        SD=a$mean*(a$CV)

        points(N1,a$mean,ylim=c(0,max(a$mean+SD)*1.05),main=zn[x],cex.main=1.75,xaxt='n',
               ylab="",pch=19,cex=2,cex.axis=1.25,col=CLOS[x])
        segments(N1,a$mean,N1,a$mean-SD,lwd=2,col=CLOS[x])
        segments(N1,a$mean,N1,a$mean+SD,lwd=2,col=CLOS[x])
        axis(1,N,F,tck=-0.015)
      }
      axis(1,seq(1,length(N),2),a$Finyear[seq(1,length(N),2)],tck=-0.025,cex.axis=1.25)
      legend('bottomleft',zn,pch=19,col=CLOS,pt.cex=1.5,cex=1.2,bty='n')
    }

    fn.fig("Avg.wgt.zn",2000,2000)
    par(mfcol=c(1,1),las=1,mai=c(0.45,0.35,.1,.15),oma=c(2.25,2.25,.1,.1),mgp=c(1,.5,0))
    fn.see.avg.wgt()
    mtext("Relative live weight",2,line=0,cex=1.5,las=3,outer=T)
    mtext("Financial Year",1,cex=1.5,line=0.5,outer=T)
    dev.off()
  }
  fn.see.avg.wgt=function()
  {
    a=Avr.wt.yr
    N=1:length(unique(a$Finyear))
    SD=a$mean*(a$CV)
    plot(N,a$mean,ylim=c(0,max(a$mean+SD)*1.05),main="",cex.main=1.25,xaxt='n',
         ylab="",xlab="",pch=19,cex=2,cex.axis=1,col="steelblue",xlim=c(0,N[length(N)]+0.5))
    segments(N,a$mean,N,a$mean-SD,lwd=2,col="steelblue")
    segments(N,a$mean,N,a$mean+SD,lwd=2,col="steelblue")
    axis(1,seq(1,length(N),1),labels =F)
    axis(1,seq(1,length(N),2),a$Finyear[seq(1,length(N),2)],tck=-0.02,cex.axis=1)
  }
  fn.fig("Avg.wgt",2000,2000)
  par(mfcol=c(1,1),las=1,mai=c(0.45,0.35,.1,.15),oma=c(2.25,2.25,.1,.1),mgp=c(1,.5,0))
  fn.see.avg.wgt()
  mtext("Relative live weight",2,line=0.5,cex=1.5,las=3,outer=T)
  mtext("Financial Year",1,cex=1.5,line=0.5,outer=T)
  dev.off()


  #4. Select effort years
  Eff.zn=subset(Eff.zn,FINYEAR%in%YR.span)
  Eff.total=data.frame(FINYEAR=Eff.zn$FINYEAR,Total=Eff.zn$West+Eff.zn$Zone1+Eff.zn$Zone2)


  #5. Visualize data availability
  Yrs=YR.span #Exploitation years
  fn.search.yr=function(finyear,Y)
  {
    DAT=data.frame(FINYEAR=as.character(finyear),Occur=Y)
    ID=which(!Yrs%in%finyear)
    if(length(ID)>0)
    {
      add=data.frame(FINYEAR=Yrs[ID])
      add$Occur=NA
      DAT=rbind(DAT,add)
    }
    DAT$FINYEAR=factor(DAT$FINYEAR,levels=Yrs)
    DAT=DAT[order(DAT$FINYEAR),]
    return(DAT)
  }

  visualize.dat=function(CATCH,SIZE,ABUNDANCE,C.TAGS,A.TAGS,AVG.wt,YLIM,LABS)
  {
    n.kt=length(CATCH)
    n.sz=length(SIZE)
    n.ab=1
    n.c.T=length(C.TAGS)
    n.a.T=length(A.TAGS)

    CATCH$Rec=subset(CATCH$Rec,FINYEAR%in%YRS)

    #1. Years with catch
    Catch.data=vector('list',length(CATCH))
    for (p in 1:n.kt) Catch.data[[p]]=fn.search.yr(unique(CATCH[[p]]$FINYEAR),p)

    #2. Years with size composition
    Size.data=vector('list',length(SIZE))
    for (d in 1:n.sz)
    {
      if(class(SIZE[[d]])=='list'& length(SIZE[[d]])>1)Size.data[[d]]=fn.search.yr(unique(do.call(rbind,SIZE[[d]])$FINYEAR),d+p)
      if(class(SIZE[[d]])=='data.frame')Size.data[[d]]=fn.search.yr(unique(SIZE[[d]]$FINYEAR),d+p)
    }

    #3. Years with average weight data
    # Yrs.avg.wt=fn.search.yr(as.character(AVG.wt),1+p+d)

    #4. Years with abundance index
    Abundance.data=vector('list',length(ABUNDANCE))
    for (f in 1:n.ab) Abundance.data[[f]]=fn.search.yr(unique(as.character(ABUNDANCE$Finyear)),f+p+d)

    #5. Years with conventional tagging
    ConvTag.data=vector('list',length(C.TAGS))
    for (g in 1:n.c.T) ConvTag.data[[g]]=fn.search.yr(unique(as.character(C.TAGS[[g]]$FINYEAR)),g+f+p+d)

    #6. Years with  acoustic tagging
    AcousTag.data=fn.search.yr(as.character(A.TAGS),1+g+f+p+d)

    #average weight
    Yrs.avrg.w=fn.search.yr(AVG.wt,1+g+f+p+d+1)

    #7. Years iwth effort
    Yrs.efrt=fn.search.yr(as.character(Eff.total$FINYEAR),1+g+f+p+d+2)

    #8. Plot data
    X=1:length(Yrs)
    plot(X,X,ylim=YLIM,yaxt='n',ylab="",col="transparent",xaxt='n',xlab="")

    #add catch
    for(q in 1:n.kt) points(X,Catch.data[[q]]$Occur,cex=4,pch="-",col=COL[[1]][q])

    #add size comp
    for(q in 1:n.sz)
    {
      a=Size.data[[q]]$Occur
      if(length(a)>length(X))a=a[1:(length(a)-1)]
      points(X,a,cex=4,pch="-",col=COL[[2]][q])
    }

    #add abundance
    for(q in 1:n.ab)
    {
      a=Abundance.data[[q]]$Occur      #to remove last year of data if >2012-13
      if(length(a)>length(X))a=a[1:(length(a)-1)]
      if(length(a)>length(X))a=a[1:(length(a)-1)]
      points(X,a,cex=4,pch="-",col=COL[[3]][q])
    }

    #add Conv.Tag
    if(add.conv.tag=="YES")
    {
      for(q in 1:n.c.T)
      {
        a=ConvTag.data[[q]]$Occur
        if(length(a)>length(X))a=a[1:(length(a)-1)]
        if(length(a)>length(X))a=a[1:(length(a)-1)]
        points(X,a,cex=4,pch="-",col=COL[[4]][q])
      }
    }

    #add Acous.Tag
    if(add.conv.tag=="YES") points(X,AcousTag.data$Occur[1:length(X)],cex=4,pch="-",col=COL$acoustic.t)
    if(add.conv.tag=="NO") points(X,AcousTag.data$Occur[1:length(X)]-1,cex=4,pch="-",col=COL$acoustic.t)

    #add average weight
    points(X,Yrs.avrg.w$Occur,cex=4,pch="-",col="pink")

    #add effort
    if(add.effort=="YES")points(X,Yrs.efrt$Occur,cex=4,pch="-",col="orange")

    axis(1,X,F,tck=-0.01)
    axis(1,seq(1,length(Yrs),10),Yrs[seq(1,length(Yrs),10)],tck=-0.04)
    axis(2,1:YLIM[2],LABS[1:length(1:YLIM[2])],las=2,cex.axis=1)
    mtext("Financial year",1,cex=1.75,line=2.2)
  }

  PAR=function()par(mai=c(.3,1.55,.1,.1),oma=c(2,2,.1,.1))
  COL=list(catch=c("brown","brown2","brown3","firebrick1","brown4","chocolate",
                   "coral","coral2","chocolate4","darkred","red"),
           size=c("forestgreen","darkgreen","chartreuse","darkolivegreen1","chartreuse3"),
           abundance="gold",conv.t="blue4",acoustic.t="deepskyblue")
  LWD=10


  if(SP%in%c("WH","GM"))
  {
    if(add.conv.tag=="YES")  YLIm=c(0.5,8)
    if(add.conv.tag=="NO")  YLIm=c(0.5,7)
  }
  if(SP=="BW")
  {
    if(Ktch.source=="ALL") YLIm=c(0.5,14)
    if(Ktch.source=="WA.only") YLIm=c(0.5,11)
  }
  if(SP=="TK")
  {
    if(Ktch.source=="ALL") YLIm=c(0.5,11)
    if(Ktch.source=="WA.only") YLIm=c(0.5,10)
  }

  if(add.effort=="YES") YLIm[2]=YLIm[2]+1

  if(!add.conv.tag=="NO")
  {
    if(add.effort=="YES")
    {
      lAbs=c(paste("Catch (",names(catch),")",sep=""),paste("Size comp. (",names(All.size[1]),")",sep=""),
             "Stand. cpue (TDGDLF)","Conv. tag.","Acous. tag.","Effort (TDGDLF)")
    }
    if(add.effort=="NO")
    {
      lAbs=c(paste("Catch (",names(catch),")",sep=""),paste("Size comp. (",names(All.size[1]),")",sep=""),
             "Stand. cpue (TDGDLF)","Conv. tag.","Acous. tag.","Avrg.w")
    }
  }
  if(add.conv.tag=="NO")
  {
    lAbs=c(paste("Catch (",names(catch),")",sep=""),paste("Size comp. (",names(All.size[1]),")",sep=""),
           "Stand. cpue (TDGDLF)","Acous. tag.","Effort (TDGDLF)")
  }

  fn.fig("avail.dat",2400,1400)
  PAR()
  #par(mai=c(.7,2.1,.1,.1))
  visualize.dat(CATCH=catch,SIZE=All.size[1],ABUNDANCE=rbind(Ab.indx.TDGDLF.all,Ab.indx.TDGDLF.all.daily),
                C.TAGS=c.Tag,A.TAGS=a.Tag,AVG.wt=unique(Avr.wt.yr$Finyear),YLIM=YLIm,LABS=lAbs)
  dev.off()


  #Table of released indivuals with conventional tags
  TBL.conv.rel=merge(Zn.rel_Conv.Tag_size_adu[,-match("TG.zn",names(Zn.rel_Conv.Tag_size_adu))],
                     Zn.rel_Conv.Tag_size_juv[,-match("TG.zn",names(Zn.rel_Conv.Tag_size_juv))],
                     by=c("Rel.zone","Yr.rel"),all=T)
  names(TBL.conv.rel)[match(c("Number.x","Number.y"),
                            names(TBL.conv.rel))]=c("Adult.N","Juvenile.N")
  TBL.conv.rel=TBL.conv.rel[,-match(c("size.group.x",
                                      "size.group.y"),names(TBL.conv.rel))]
  TBL.conv.rel[is.na(TBL.conv.rel)]=0
  write.csv(TBL.conv.rel,"TBL.conv.releases.csv",row.names=F)
  fn.word.table(WD=getwd(),TBL=TBL.conv.rel,Doc.nm="TBL.conv.releases",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")


  #Plot acoustic tag recaptures by year
  Rep.Recap=subset(Rep.Recap,!is.na(year.rel))
  i.yrs=sort(unique(Rep.Recap$year.rec))
  Y.rng=range(c(Rep.Recap$RELLATDECDEG,Rep.Recap$RECLATDECDEG),na.rm=T)
  X.rng=range(c(Rep.Recap$RELLNGDECDEG,Rep.Recap$RECLNGDECDEG),na.rm=T)

  fn.fig("Acoustic.tags.rel.rec.map",2400,2400)
  par(mfcol=c(ceiling(length(i.yrs)/2),2),las=1,mai=c(1,1,.1,.1),oma=c(.1,.1,.1,.1))
  dummy=Rep.Recap
  for(i.yr in 1:length(i.yrs))
  {
    recs=subset(dummy,year.rec==i.yrs[i.yr])
    na.rcs=subset(recs,is.na(RECLATDECDEG))
    if(nrow(recs)>0)
    {
      rels=subset(dummy,year.rel<=i.yrs[i.yr])
      plot(rels$RELLNGDECDEG,rels$RELLATDECDEG,ylim=Y.rng,xlim=X.rng,ylab="",xlab="",pch=19)
      points(recs$RECLNGDECDEG,recs$RECLATDECDEG,col=2,pch=4,cex=1.75)
      if(nrow(na.rcs)==0) Rc.lg=paste(nrow(recs),"recaptured")
      if(nrow(na.rcs)>0) Rc.lg=paste(nrow(recs),"recaptured (",nrow(na.rcs),"without rec. info)")
      legend('top',c(paste(nrow(rels),"with tags"),Rc.lg),bty='n',col=c(1,2),
             text.col=c(1,2),pch=c(19,4),title=paste(i.yrs[i.yr]),cex=1.2)
      dummy=subset(dummy,!ATAG_NO%in%recs$ATAG_NO)
    }

  }
  mtext("Longitude",1,-2,outer=T,cex=1.75)
  mtext("Latitude",2,-2,las=3,outer=T,cex=1.75)
  dev.off()


  #6.Export data for use in models
  HandL="C:/Matias/Data/Population dynamics/Data inputs for models/"
  DiR=paste(HandL,SP.LABELS1,"/",Yr.assess,sep='')
  if(!file.exists(DiR)) dir.create(DiR)
  setwd(DiR)

  fn.agg.at.level.and.exprt=function(DAT,Level,VAR,SOURCE)
  {
    for(i in 1:length(DAT))
    {
      if(Level=='annual')
      {
        a=aggregate(LIVEWT.c~FINYEAR,DAT,sum)
        write.csv(a,paste(VAR,".",Level,".",SOURCE,".csv",sep=""),row.names=F)
      }

      if(Level=='annual.by.zone')
      {
        if(!is.null(DAT$zone))
        {
          a=aggregate(LIVEWT.c~FINYEAR+zone,DAT,sum)
          write.csv(a,paste(VAR,".",Level,".",SOURCE,".csv",sep=""),row.names=F)
        }
      }
    }
  }

  #Catch
  for(i in 1:length(catch))
  {
    fn.agg.at.level.and.exprt(catch[[i]],'annual',"ktch",names(catch)[i])
    fn.agg.at.level.and.exprt(catch[[i]],'annual.by.zone',"ktch",names(catch)[i])
  }

  #Effort
  write.csv(Eff.zn,"effort.annual.by.zone.TDGDLF.csv",row.names=F)
  write.csv(Eff.total,"effort.annual.TDGDLF.csv",row.names=F)


  #Size composition
  fn.exp.size=function(dat,Level)
  {
    n=length(dat)
    for( i in 1:n)
    {
      a=dat[[i]]
      if(class(a)=="data.frame")
      {
        write.csv(a,paste("size.comp.",names(dat)[i],".csv",sep=""),row.names=F)
      }
      if(class(a)=="list"& length(a)>1)
      {
        if(Level=='annual')
        {
          a=do.call(rbind,a)
          YRs=as.character(unique(a$FINYEAR))
          SIzes=colnames(a[,4:ncol(a)])
          d=matrix(ncol=length(SIzes),nrow=length(YRs))
          d.f=d.m=d
          for(s in 1:length(YRs))
          {
            ddat=subset(a,FINYEAR==YRs[s])
            ddat.f=subset(ddat,Sex=="F")
            ddat.m=subset(ddat,Sex=="M")

            d[s,]=colSums(ddat[,4:ncol(ddat)])
            d.f[s,]=colSums(ddat.f[,4:ncol(ddat.f)])
            d.m[s,]=colSums(ddat.m[,4:ncol(ddat.m)])
          }
          #sexes combined
          colnames(d)=SIzes
          a=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
          write.csv(a,paste("size.comp.",Level,".",names(dat)[i],".csv",sep=""),row.names=F)

          #by sex
          colnames(d.f)=colnames(d.m)=SIzes
          a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.f))
          write.csv(a,paste("size.comp.fem.",Level,".",names(dat)[i],".csv",sep=""),row.names=F)

          a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.m))
          write.csv(a,paste("size.comp.mal.",Level,".",names(dat)[i],".csv",sep=""),row.names=F)

        }

        if(Level=='annual.by.zone')
        {
          ss=names(a)
          for(p in 1:length(ss))
          {
            b=a[[p]]
            if(nrow(b)>0)
            {
              YRs=as.character(unique(b$FINYEAR))
              SIzes=colnames(b[,4:ncol(b)])
              d=matrix(ncol=length(SIzes),nrow=length(YRs))
              d.f=d.m=d
              for(s in 1:length(YRs))
              {
                ddat=subset(b,FINYEAR==YRs[s])
                ddat.f=subset(ddat,Sex=="F")
                ddat.m=subset(ddat,Sex=="M")

                d[s,]=colSums(ddat[,4:ncol(ddat)])
                d.f[s,]=colSums(ddat.f[,4:ncol(ddat.f)])
                d.m[s,]=colSums(ddat.m[,4:ncol(ddat.m)])

              }

              #sexes combined
              colnames(d)=SIzes
              b=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
              write.csv(b,paste("size.comp.annual.",ss[p],".",names(dat)[i],".csv",sep=""),row.names=F)

              #by sex
              colnames(d.f)=colnames(d.m)=SIzes
              b=cbind(FINYEAR=YRs,as.data.frame.matrix(d.f))
              write.csv(b,paste("size.comp.annual.fem.",ss[p],".",names(dat)[i],".csv",sep=""),row.names=F)

              b=cbind(FINYEAR=YRs,as.data.frame.matrix(d.m))
              write.csv(b,paste("size.comp.annual.mal.",ss[p],".",names(dat)[i],".csv",sep=""),row.names=F)
            }

          }
        }
      }
    }
  }
  fn.exp.size(dat=All.size,Level='annual')
  fn.exp.size(dat=All.size,Level='annual.by.zone')


  #avg weight
  if(exists('Avr.wt.yr.zn'))write.csv(Avr.wt.yr.zn,"ktch.avg.weight.annual.by.zone.csv",row.names=F)
  write.csv(Avr.wt.yr,"ktch.avg.weight.annual.csv",row.names=F)


  #TDGLDF cpue
  #folly
  if(exists("Ab.folly.TDGDLF.all")) write.csv(Ab.folly.TDGDLF.all,"cpue.annual.TDGDLF.folly.csv",row.names=F)

  #by zone
  write.csv(Ab.indx.TDGDLF,"cpue.annual.by.zone.TDGDLF.csv",row.names=F)
  write.csv(Ab.indx.TDGDLF.daily,"cpue.annual.by.zone.TDGDLF.daily.csv",row.names=F)

  #zones combined
  write.csv(Ab.indx.TDGDLF.all,"cpue.annual.TDGDLF.csv",row.names=F)
  write.csv(Ab.indx.TDGDLF.all.daily,"cpue.annual.TDGDLF.daily.csv",row.names=F)


  #Naturaliste survey
  if(exists('Naturaliste.abun')) write.csv(Naturaliste.abun,"cpue.annual.survey.csv",row.names=F)
  if(exists('Naturaliste.size')) write.csv(Naturaliste.size,"size.annual.survey.csv",row.names=F)


  #Conventional tagging

  #Individual-based model
  write.csv(Rel_rec_Conv.Tag,"Ind_based_model.csv",row.names=F)


  #1. At age
  #note: this function puts in SS3 format
  tagging.fn=function(dat,THESE,THOSE,type,txt)
  {
    dat=dat[,match(THESE,names(dat))]
    if(type=='rel')
    {
      dat$gender=0
      dat=dat[,match(c("TG","Rel.zone","FinYear.rel","Mn.rel","gender","Age","Number"),names(dat))]
    }

    if(type=='rec')dat=dat[,match(c("TG","FinYear.rec","Mn.rec","Rec.zone","Number"),names(dat))]
    names(dat)=THOSE
    write.csv(dat,paste(txt,"csv",sep=""),row.names=F)
  }

  #releases
  these=c("TG","Rel.zone","FinYear.rel","Mn.rel","Age","Number")
  those=c("TG","area","yr","season","gender","Age","Nrelease")
  tagging.fn(Zn.rel_Conv.Tag,these,those,'rel','conv.tag.rel.')

  #recaptures
  these=c("TG","Rec.zone","FinYear.rec","Mn.rec","Number")
  those=c("TG","year","season","area","Number")
  tagging.fn(Zn.rec_Conv.Tag,these,those,'rec',"conv.tag.reca.")


  #2. At size
  #all
  write.csv(Zn.rel_Conv.Tag_size,"Zn.rel_Conv.Tag_size.csv",row.names=F)
  write.csv(Zn.rec_Conv.Tag_size,"Zn.rec_Conv.Tag_size.csv",row.names=F)

  #juveniles and adults
  write.csv(Zn.rel_Conv.Tag_size_adu,"Zn.rel_Conv.Tag_size_adu.csv",row.names=F)
  write.csv(Zn.rel_Conv.Tag_size_juv,"Zn.rel_Conv.Tag_size_juv.csv",row.names=F)
  write.csv(Zn.rec_Conv.Tag_size_adu,"Zn.rec_Conv.Tag_size_adu.csv",row.names=F)
  write.csv(Zn.rec_Conv.Tag_size_juv,"Zn.rec_Conv.Tag_size_juv.csv",row.names=F)
  write.csv(Smallest_size_tagged,"Smallest_size_tagged.csv",row.names=F)


  #Acoustic tagging
  #releases
  write.csv(Zn.rel_Acous.Tag,"acous.tag.rel.csv",row.names=F)
  write.csv(Zn.rel_Acous.Tag.prop,"acous.tag.rel.prop.csv",row.names=F)

  #recaptures
  write.csv(Zn.rec_Acous.Tag,"acous.tag.reca.csv",row.names=F)
  write.csv(Zn.rec_Acous.Tag.prop,"acous.tag.reca.prop.csv",row.names=F)

  #Individual-based model
  write.csv(Indiv_based_Acous.Tag,"acous.tag.reca.indv.based.csv",row.names=F)


  #Age and growth
  if(exists('Age.growth'))  write.csv(Age.growth,"Age.growth.csv",row.names=F)


  #Recapture information from acoustic tags
  write.csv(Rep.Recap,"Acoustic.tag.rel.rec.csv",row.names=F)


  #Mesh proportional effort
  write.csv(Mesh.prop.eff,"Mesh.prop.eff.csv",row.names=F)
  write.csv(Mesh.prop.eff.West,"Mesh.prop.eff.West.csv",row.names=F)
  write.csv(Mesh.prop.eff.Zn1,"Mesh.prop.eff.Zn1.csv",row.names=F)
  write.csv(Mesh.prop.eff.Zn2,"Mesh.prop.eff.Zn2.csv",row.names=F)


}