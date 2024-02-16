fn.get.stuff.from.list=function(lista,stuff) lapply(lista, function(x) x[[stuff]]) 
fn.get.and.name=function(LISTA,x)
{
  dd=fn.get.stuff.from.list(LISTA,x)
  names(dd)=names(LISTA)
  return(dd)
}
fn.extract.dat.perl=function(STRING,nm.Dat) grepl(STRING, nm.Dat, perl = TRUE)
objects.exist <- function(...)
{
  ls <- list(...)
  sapply(ls, exists)
}
fn.ktch.sex.ratio.zone=function(size.data)
{
  if(any(grepl('Observations',names(size.data)))) size.data=size.data[-grep('Observations',names(size.data))]
  if(sum(grepl('Table',names(size.data)))>0) size.data=size.data[-grep('Table',names(size.data))]
  size.data=do.call(rbind,size.data)%>%
    rownames_to_column(var = 'Zone')%>%
    mutate(N=1,
           Zone=gsub("\\..*","",str_remove(Zone,'Size_composition_')))%>%
    filter(!is.na(SEX))%>%
    group_by(Zone,year,SEX)%>%
    summarise(n=sum(N))%>%
    ungroup()%>%
    spread(SEX,n,fill=0)%>%
    mutate(Total=F+M,
           Prop.female=F/Total)%>%
    filter(Total>Min.obs)
  if(nrow(size.data)>0)
  {
    print(size.data%>%
            ggplot(aes(year,Prop.female,color=log(Total)))+
            geom_point()+
            facet_wrap(~Zone,ncol=1))+
      theme_PA()
    return(size.data%>%group_by(Zone)%>%summarise(Prop.female=mean(Prop.female)))
  }
}
# Import Total Catch------------------------------------------------------
fn.import.catch.data=function(KTCH.UNITS)
{
  #2.1. Commercial
  #2.1.1 Catch_WA Fisheries
  
  #..Historic
  Historic.ktch=fn.in(NM='recons_Hist.expnd.csv')%>%filter(!FINYEAR=='1975-76')
  
  #..Ammended reported commercial catch including discards
  Data.monthly=fn.in(NM='recons_Data.monthly.csv')
  Data.monthly.north=fn.in(NM='recons_Data.monthly.north.csv')
  
  #..TEPS
  Sawfish_Pilbara.ktch=fn.in(NM='recons_Pilbara_sawfish.ktch.csv')
  Greynurse.ktch=fn.in(NM='recons_Greynurse.ktch.csv')
  TEPS_dusky=fn.in(NM='recons_TEPS_dusky.csv')
  Sawfish_KPTF.ktch=fn.in(NM='recons_Kimberly.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_NBPTF.ktch=fn.in(NM='recons_Nickol.Bay.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_OPTF.ktch=fn.in(NM='recons_Onslow.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_SBPTF.ktch=fn.in(NM='recons_Shark.Bay.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_EGPTF.ktch=fn.in(NM='recons_Exmouth.Gulf.Prawn.Trawl_sawfish.ktch.csv')
  Sawfish_BPTF.ktch=fn.in(NM='recons_Broome.Prawn.Trawl_sawfish.ktch.csv')
  
  #..Droplines Western Rock Lobster
  WRL.ktch=fn.in(NM='Wetline_rocklobster.csv')
  
  #..Discards from TDGDLF
  if(TDGLDF.disc.assumed.PCM=="BaseCase") discard_TDGDLF=fn.in(NM='recons_discard_TDGDLF.csv')
  if(TDGLDF.disc.assumed.PCM=="100%") discard_TDGDLF=fn.in(NM='recons_discard_TDGDLF_100.PCM')
  
  
  #2.1.2. Catch of non WA Fisheries
  
  #..Taiwanese gillnet, longline and trawl
  Taiwan.gillnet.ktch=fn.in(NM='recons_Taiwan.gillnet.ktch.csv')
  Taiwan.longline.ktch=fn.in(NM='recons_Taiwan.longline.ktch.csv')
  Taiwan.trawl.ktch=fn.in(NM='recons_Taiwan.trawl.ktch.csv')
  
  Taiwan.gillnet.ktch$Method="GN"   #Pelagic.gillnet
  Taiwan.longline.ktch$Method="LL"  #Longline
  Taiwan.trawl.ktch$Method="TW"     #Trawl
  
  Taiwan=rbind(Taiwan.longline.ktch,Taiwan.gillnet.ktch,Taiwan.trawl.ktch)%>%
    mutate(Region="North")%>%
    mutate(year=as.numeric(substr(FINYEAR,1,4)))
  
  #..Indonesian illegal fishing in Australia waters
  Indo_total.annual.ktch=fn.in(NM='recons_Indo.IUU.csv')
  Indo_total.annual.ktch=Indo_total.annual.ktch%>%filter(!is.na(SPECIES))
  
  #..AFMA's GAB & SBT fisheries
  GAB.trawl_catch=fn.in(NM='recons_GAB.trawl_catch.csv') 
  WTBF_catch=fn.in(NM='recons_WTBF_catch.csv')
  
  #..SA Marine Scalefish fishery
  Whaler_SA=fn.in(NM='recons_Whaler_SA.csv') 
  
  #..NSW Fisheries (with beach protection)
  Bronze.whaler_NSW=fn.in(NM='recons_Bronzewhaler_NSW.csv') 
  
  #..NT catches
  NT_catch=fn.in(NM='recons_NT_catch.csv') 
  
  
  #2.2. WA Recreational catch
  Rec.ktch=fn.in(NM='recons_recreational.csv')  
  Rec.ktch=Rec.ktch%>%mutate(Region=ifelse(zone%in%c('Gascoyne','North Coast'),'North','South'),
                             year=as.numeric(substr(FINYEAR,1,4)),
                             Common.Name=tolower(Common.Name))
  #add dummy 0 catch for these species to get catch series starting in 1941
  if(!c("weasel shark")%in%unique(Rec.ktch$Common.Name))
  {
    Addis=Rec.ktch%>%
              filter(Common.Name=='zebra shark')%>%
              mutate(LIVEWT.c=0,
                     Common.Name="weasel shark")
    
    Rec.ktch=rbind(Rec.ktch,Addis)
  }
  if(!c("snaggletooth")%in%unique(Rec.ktch$Common.Name))
  {
    Addis=Rec.ktch%>%
      filter(Common.Name=='zebra shark')%>%
      mutate(LIVEWT.c=0,
             Common.Name="snaggletooth")
    Rec.ktch=rbind(Rec.ktch,Addis)
  }
  
  #Combine data sets
  # TDGDLF and NSF
  Data.monthly$Region="South"
  Data.monthly.north$Region="North"
  Data.monthly$Data.set="Data.monthly"
  Data.monthly.north$Data.set="Data.monthly.north"
  Tot.ktch=rbind(subset(Data.monthly,select=c(SPECIES,FishCubeCode,FINYEAR,Region,Data.set,LIVEWT.c,BLOCKX,METHOD,zone)),
                 subset(Data.monthly.north,select=c(SPECIES,FishCubeCode,FINYEAR,Region,Data.set,LIVEWT.c,BLOCKX,METHOD,zone)))%>%
    left_join(All.species.names,by='SPECIES')%>%     #add species names
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%    
    mutate(Name=SNAME)
  
  #Add Recreational
  a=Rec.ktch%>%
    rename(finyear=year)%>%
    mutate(BLOCKX=NA,
           zone=case_when(zone%in%c("North Coast","Gascoyne Coast")~'North',
                          zone=="West Coast"~'West',
                          zone=="South Coast"~'Zone1'),
           Common.Name=ifelse(Common.Name=="dogfishes","spurdogs",
                       ifelse(Common.Name=="greynurse shark","grey nurse shark",
                       ifelse(Common.Name=="thresher shark","thresher sharks",
                       ifelse(Common.Name=="bronze whaler","copper shark",
                       ifelse(Common.Name=="dusky whaler","dusky shark",
                       ifelse(Common.Name=="gummy sharks","gummy shark",
                       ifelse(Common.Name=="tawny shark","tawny nurse shark",
                       ifelse(Common.Name=="sawfishes","sawfish (general)",
                       ifelse(Common.Name=="australian blacktip shark","australian blacktip",
                       Common.Name))))))))))%>%
    left_join(All.species.names,by=c('Common.Name'='SNAME'))%>%
    mutate(SNAME=Common.Name,
           Name=Common.Name,
           METHOD='Rec.line',
           Data.set="Recreational",
           FishCubeCode="Recreational")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  #Add Discards from TDGDLF
  a=discard_TDGDLF%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Name=SNAME,
           Type="Discards_TDGDLF",
           Data.set="Discards_TDGDLF",
           METHOD='GN',
           Region='South',
           FishCubeCode="Discards_TDGDLF",
           finyear=as.numeric(substr(FINYEAR,1,4)))%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  #Add Taiwanese
  a=Taiwan%>%
    left_join(All.species.names,by='SPECIES')%>%
    rename(finyear=year)%>%
    mutate(BLOCKX=NA,
           Name=SNAME,
           Type="Taiwan",
           Data.set="Taiwan",
           METHOD=Method,
           FishCubeCode="Taiwan")%>%
    dplyr::select(names(Tot.ktch))%>%
    filter(SPECIES%in%unique(Tot.ktch$SPECIES))
  Tot.ktch=rbind(Tot.ktch,a)
  
  #Add Indonesian fishing incursions
  a=Indo_total.annual.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="Indonesia",
           Data.set="Indonesia",
           METHOD=NA,
           FishCubeCode="Indo")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  #Add TEP interactions     
  a=Sawfish_Pilbara.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.Pilbara",
           Data.set="TEP_sawfish.Pilbara",
           METHOD='TW',
           FishCubeCode="PFT")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_KPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.KPTF",
           Data.set="TEP_sawfish.KPTF",
           METHOD='TW',
           FishCubeCode="KP")%>%   
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_NBPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.NBPTF",
           Data.set="TEP_sawfish.NBPTF",
           METHOD='TW',
           FishCubeCode="NBP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_OPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.OPTF",
           Data.set="TEP_sawfish.OPTF",
           METHOD='TW',
           FishCubeCode="OP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_SBPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.SBPTF",
           Data.set="TEP_sawfish.SBPTF",
           METHOD='TW',
           FishCubeCode="SBP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_EGPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.EGPTF",
           Data.set="TEP_sawfish.EGPTF",
           METHOD='TW',
           FishCubeCode="EGP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Sawfish_BPTF.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           zone="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_sawfish.BPTF",
           Data.set="TEP_sawfish.BPTF",
           METHOD='TW',
           FishCubeCode="BP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=Greynurse.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_greynurse",
           Data.set="TEP_greynurse",
           METHOD='GN',
           FishCubeCode="TEP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  a=TEPS_dusky%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="TEP_dusky",
           Data.set="TEP_dusky",
           METHOD='GN',
           FishCubeCode="TEP")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  #Add WRL        
  a=WRL.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="WRL",
           Data.set="WRL",
           METHOD='LL',
           FishCubeCode="WRL")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  
  #Add historic  
  a=Historic.ktch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Type="Commercial",
           Data.set="Historic",
           METHOD=NA,
           FishCubeCode="Historic",
           SNAME=ifelse(SNAME%in%c('southern sawshark','common sawshark'),'sawsharks',SNAME),
           SPECIES=ifelse(SPECIES%in%c(23001,23002),23900,SPECIES),
           Name=SNAME)%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  #Add GAB
  a=GAB.trawl_catch%>%
    mutate(SPECIES=ifelse(SPECIES==5000,5002,
                   ifelse(SPECIES==23000,23900,
                   ifelse(SPECIES==12901,12000,
                   ifelse(SPECIES==24000,24900,
                   ifelse(SPECIES==19000,19004,
                   SPECIES))))))%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="GAB",
           Data.set="GAB.trawl",
           METHOD="trawl",
           FishCubeCode="GAB")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  #Add WTBF
  a=WTBF_catch%>%
    mutate(SPECIES=ifelse(SPECIES==12901,12000,
                          ifelse(SPECIES==18901,18014,
                                 ifelse(SPECIES==19000,19001,
                                        SPECIES))))%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="WTB",
           Data.set="WTBF",
           METHOD="LL",     #line
           FishCubeCode="WTB")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  #Add NT_catch
  a=NT_catch%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="North",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="NT",
           Data.set="NT",
           METHOD="LL",   #line
           FishCubeCode="NT")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  
  #Add Whaler_SA (SA Marine Scalefish Fishery)
  a=Whaler_SA%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="SA MSF",
           Data.set="SA MSF",
           METHOD="LL",    #line
           FishCubeCode="SA MSF")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  
  #Add Bronze.whaler_NSW (this includes beach protection)
  a=Bronze.whaler_NSW%>%
    left_join(All.species.names,by='SPECIES')%>%
    mutate(BLOCKX=NA,
           Region="South",
           finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=SNAME,
           Type="NSW fisheries",
           Data.set="NSW fisheries",
           METHOD="mixed",    #line, trawl, etc
           FishCubeCode="NSW fisheries")%>%
    dplyr::select(names(Tot.ktch))
  Tot.ktch=rbind(Tot.ktch,a)
  
  #set catch in tonnes
  if(KTCH.UNITS=="TONNES")Tot.ktch=Tot.ktch%>%mutate(LIVEWT.c=LIVEWT.c/1000) 
  
  #fix some species names
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%19003,"winghead hammerhead",SNAME),
           Name=ifelse(SPECIES%in%19003,'winghead hammerhead',Name),
           Scien.nm=ifelse(SPECIES%in%19003,'eusphyra blochii',Scien.nm))
  
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SNAME=='sawfish (general)',"sawfish",SNAME),
           Scien.nm=ifelse(Name=='sawfish (general)','Pristidae',Scien.nm),
           SPECIES=ifelse(Name=='sawfish (general)',25000,SPECIES),
           Name=ifelse( Name=='sawfish (general)','sawfish',Name))
  
  Tot.ktch=Tot.ktch%>%
    mutate(SPECIES=ifelse(Name=='rays & skates',31000,SPECIES),
           Name=ifelse( Name=='skates','rays & skates',Name))
  
  Tot.ktch=Tot.ktch%>%
    mutate(Scien.nm=ifelse(Name=="australian blacktip",'Carcharhinus tilstoni',Scien.nm))
  
  Tot.ktch=Tot.ktch%>%
    mutate(SPECIES=ifelse(SPECIES==20904,20906,SPECIES),
           SNAME=ifelse(SPECIES==20906,"roughskin shark",SNAME),
           Scien.nm=ifelse(SPECIES==20906,'Centroscymnus spp.',Scien.nm),
           Name=ifelse( SPECIES==20906,'roughskin shark',Name))
  
  
  #list of all species caught
  Table1=Tot.ktch%>%  
    filter(!SPECIES%in%c(31000,22999))%>%
    group_by(Name,FishCubeCode)%>%
    summarise(Catch.tons=sum(LIVEWT.c))%>%
    spread(FishCubeCode,Catch.tons,fill=0)%>%
    left_join(Tot.ktch%>%
                distinct(SPECIES, .keep_all = T)%>%
                arrange(SPECIES)%>%
                dplyr::select(Name,Scien.nm),by="Name")%>%
    mutate(Name=capitalize(Name))%>%
    relocate(where(is.numeric), .after = where(is.character))
  
  
  #For assessments, combine certain species with unreliable reporting resolution
  unik=unique(Tot.ktch$Name)
  #Catsharks
  this=unik[grep("catshark",unik)]
  this=this[-match("brown-banded catshark", this)] #not a catshark
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(Name%in%this,"catsharks",SNAME),
           Scien.nm=ifelse(Name%in%this,'Scyliorhinidae',Scien.nm),
           SPECIES=ifelse(Name%in%this,15000,SPECIES),
           Name=ifelse(Name%in%this,'catsharks',Name))
  
  #Sawsharks
  this=c(23000,23001,23002,23900)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,"sawsharks",SNAME),
           Name=ifelse(SPECIES%in%this,'sawsharks',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Pristiophoridae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,23900,SPECIES),
           Name=ifelse(Name=='sawshark','sawsharks',Name))
  
  #Threshers
  this=c(12000,12001)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'thresher sharks',SNAME),
           Name=ifelse(SPECIES%in%this,'thresher sharks',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Alopidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,12000,SPECIES))   
  
  #Wobbegongs
  this=c(13000, 13001,13003,13011,13012,13017,13022,13016,13021)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'wobbegongs',SNAME),
           Name=ifelse(SPECIES%in%this,'wobbegongs',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Orectolobidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,13000,SPECIES)) 
  
  #Angelsharks
  this=c(24900, 24001,24002)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'angel sharks',SNAME),
           Name=ifelse(SPECIES%in%this,'angel sharks',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Squatinidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,24900,SPECIES))  
  
  #Stingrays
  this=c(35001, 35000)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'smooth stingray',SNAME),
           Name=ifelse(SPECIES%in%this,'smooth stingray',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Bathytoshia brevicaudata',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,35001,SPECIES))  
  
  #Banjo rays
  this=c(26999,27001,27007,27011,27909)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'banjo rays',SNAME),
           Name=ifelse(SPECIES%in%this,'banjo rays',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Trygonorrhinidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,26999,SPECIES))     
  
  #Wedgefishes
  this=c(26000,26001,26002)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%this,'wedgefishes',SNAME),
           Name=ifelse(SPECIES%in%this,'wedgefishes',Name),
           Scien.nm=ifelse(SPECIES%in%this,'Rhinidae',Scien.nm),
           SPECIES=ifelse(SPECIES%in%this,26000,SPECIES))       
  
  #Rays and skates  
  Tot.ktch=Tot.ktch%>%
    mutate(Name=ifelse(SPECIES==31000,'rays & skates',Name),
           Scien.nm=ifelse(SPECIES==31000,'Rajidae & Arhynchobatidae',Scien.nm),
           SNAME=ifelse(SPECIES==31000,'rays & skates',SNAME))
  
  #Australian blacktip
  Tot.ktch=Tot.ktch%>%
    mutate(Name=ifelse(SNAME=='australian blacktip','blacktips',Name),
           Scien.nm=ifelse(SNAME=='australian blacktip','C. limbatus & C. tilstoni',Scien.nm),
           SPECIES=ifelse(SNAME=='australian blacktip',18014,SPECIES),
           SNAME=ifelse(SNAME=='australian blacktip','blacktips',SNAME))
  
  
  #Bull to pigeye shark (bull not likely to be taken: Heupel & McAuley 2007 page  84)
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(SPECIES%in%c(18021),"pigeye shark",SNAME),
           Name=ifelse(SPECIES%in%c(18021),'pigeye shark',Name),
           Scien.nm=ifelse(SPECIES%in%c(18021),'C. amboinensis',Scien.nm),
           SPECIES=ifelse(SPECIES%in%c(18021),18026,SPECIES))
  
  #Reset landed 'rays and skates' in TDGDLF to eagle ray 
  #note: fishers interviews indicated that landed rays and skate are eagle rays
  #       so adjust discards estimates to only those fishers discarding them as
  #       not all fishers retain them
  Tot.ktch=Tot.ktch%>%
    mutate(LIVEWT.c=ifelse(FishCubeCode%in%c('Discards_TDGDLF')
                           & SPECIES==39001,LIVEWT.c*prop.disc.ER,LIVEWT.c),
           Name=ifelse(FishCubeCode%in%c('JASDGDL','Historic','WCDGDL')
                       & SPECIES==31000,'eagle ray',Name),
           Scien.nm=ifelse(FishCubeCode%in%c('JASDGDL','Historic','WCDGDL')
                           & SPECIES==31000,'Myliobatis tenuicaudatus',Scien.nm),
           SNAME=ifelse(FishCubeCode%in%c('JASDGDL','Historic','WCDGDL')
                        & SPECIES==31000,'eagle ray',SNAME),
           SPECIES=ifelse(FishCubeCode%in%c('JASDGDL','Historic','WCDGDL')
                          & SPECIES==31000,39001,SPECIES))
  
  
  #Disaggregate Sawfish into species  
  Prop.by.sawfish.sp=Tot.ktch%>%
    filter(SPECIES%in%25001:25020)%>%
    group_by(SPECIES,FishCubeCode)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))%>%
    group_by(FishCubeCode)%>%
    mutate(Prop=LIVEWT.c/sum(LIVEWT.c))%>%
    dplyr::select(-LIVEWT.c)%>%
    ungroup()%>%
    data.frame
  Prop.by.sawfish.sp=rbind(Prop.by.sawfish.sp,
                           Prop.by.sawfish.sp[1:2,]%>%
                             mutate(SPECIES=25001:25002,
                                    FishCubeCode=rep('Recreational',2),
                                    Prop=c(.5,.5)))
  
  sawfish=Tot.ktch%>%
    filter(SPECIES==25000)
  
  Tot.ktch=Tot.ktch%>%
    filter(!SPECIES==25000)
  
  sawfish=sawfish%>%
    mutate(SPECIES=ifelse(FishCubeCode%in%c('JANS','OANCGC'),25001,SPECIES))  #set to green sawfish based on Survey data
  
  sawfish_JANS.OANCGC=sawfish%>%
    filter(SPECIES==25001)
  
  sawfish=sawfish%>%
    filter(!SPECIES==25001)
  
  sawfish=full_join(sawfish%>%dplyr::select(-SPECIES),    #set to proportion by species
                    Prop.by.sawfish.sp%>%
                      filter(FishCubeCode%in%sawfish$FishCubeCode),
                    by='FishCubeCode',relationship = "many-to-many")%>%
    mutate(LIVEWT.c=LIVEWT.c*Prop)%>%
    dplyr::select(-Prop)
  
  sawfish=rbind(sawfish,
                sawfish_JANS.OANCGC)%>%
    mutate(SNAME=case_when(SPECIES==25001~'green sawfish',
                           SPECIES==25002~'narrow sawfish',
                           SPECIES==25004~'dwarf sawfish'),
           Name=case_when(SPECIES==25001~'green sawfish',
                          SPECIES==25002~'narrow sawfish',
                          SPECIES==25004~'dwarf sawfish'),
           Scien.nm=case_when(SPECIES==25001~'Pristis zijsron',
                              SPECIES==25002~'Anoxypristis cuspidata',
                              SPECIES==25004~'Pristis clavata'))
  
  Tot.ktch=rbind(Tot.ktch,sawfish)
  
  #re calculate catch after reapportioning
  Tot.ktch=Tot.ktch%>%
    group_by(across(c(-LIVEWT.c)))%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))%>%
    ungroup()
  
  #Fix Port Jackson name  
  Tot.ktch=Tot.ktch%>%
    mutate(Name=ifelse(Name=="port jackson shark","port Jackson shark",Name))
  
  
  #Final touches
  Tot.ktch=Tot.ktch%>%
    mutate(SNAME=ifelse(Name=='eusphyra blochii','winghead shark',Name))
  
  Tot.ktch=Tot.ktch%>%filter(!is.na(Name))
  
  
  #Aggregate by species, year and fishery
  Tot.ktch.zone=Tot.ktch%>%
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=ifelse(SPECIES%in%c(22999),"unidentified sharks",SNAME))%>%
    group_by(SPECIES,Name,FishCubeCode,Data.set,finyear,FINYEAR,Region,zone)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))
  
  Tot.ktch.method=Tot.ktch%>%
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)),
           Name=ifelse(SPECIES%in%c(22999),"unidentified sharks",SNAME))%>%
    group_by(SPECIES,Name,FishCubeCode,Data.set,finyear,FINYEAR,METHOD)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))%>%
    mutate(Gear=ifelse(METHOD%in%c("BS","BH","GN","HN","Pelagic.gillnet","PS"),"net",
                       ifelse(METHOD%in%c("DL","DV","EL","GL","HL","HR","HY",
                                          "LL","Longline","Rec.line","TL"),'line',
                              ifelse(METHOD%in%c("FG","TW"),'trawl',
                                     ifelse(METHOD%in%c("FT","PT","CT"),'trap',
                                            NA)))))
  
  Tot.ktch.method=Tot.ktch.method%>%
    group_by(FishCubeCode) %>%
    arrange(FishCubeCode, is.na(Gear)) %>% # in case to keep non- NA elements for a tie
    mutate(Gear = ifelse(is.na(Gear),Mode(Gear),Gear),
           Gear=ifelse(is.na(Gear) & FishCubeCode %in% c('Indo','WTB','SA MSF'),'line',
                       ifelse(is.na(Gear) & FishCubeCode %in% c('GAB','PFT','SBSC'),'trawl',
                              Gear)))
  
  
  Tot.ktch=Tot.ktch%>%
    mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
    group_by(SPECIES,Name,FishCubeCode,Data.set,finyear,FINYEAR,Region)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c))
  
  return(list(Total=Tot.ktch,Zone=Tot.ktch.zone,Total.method=Tot.ktch.method,Table1=Table1))
}

# PSA (which species to assess quantitatively)------------------------------------------------------
PSA.fn=function(d,line.sep,size.low,size.med,size.hig,W,H)  
{
  set.seed(101)
  PSA=data.frame(Species=d$Species,
                 Productivity=rep(NA,nrow(d)),
                 Susceptibility=rep(NA,nrow(d)),
                 Vulnerability=rep(NA,nrow(d)))
  for(p in 1:nrow(d))    
  {
    aa=d[p,]
    if(!aa$Species%in%species.meeting.criteria) #reset availability and encounterability if not meeting criteria
    {
      k=KIP%>%filter(Name==aa$Species)
      if(nrow(k)==0)  #create dummy if species is not at all in KIP (i.e. does not meet criteria for any gear)
      {
        k=KIP[1,]%>%mutate(Name=aa$Species,Gear='Dummy')
      }
      aa=aa%>%
        mutate(Net.avail=ifelse(!'net'%in%k$Gear,1,Net.avail),
               Net.encoun=ifelse(!'net'%in%k$Gear,1,Net.encoun),
               Line.avail=ifelse(!'line'%in%k$Gear,1,Line.avail),
               Line.encoun=ifelse(!'line'%in%k$Gear,1,Line.encoun),
               Trawl.avail=ifelse(!'trawl'%in%k$Gear,1,Trawl.avail),
               Trawl.encoun=ifelse(!'trawl'%in%k$Gear,1,Trawl.encoun),
               Trap.avail=ifelse(!'trap'%in%k$Gear,1,Trap.avail),
               Trap.encoun=ifelse(!'trap'%in%k$Gear,1,Trap.encoun))
    }
    if(aa$Species=="longfin mako shark")  #as Productivity = 3, only option for dropping from Assessment due to no catch 
    {                                      # is to set all Susceptibility pars to 1
      aa[,-match(c('Species','Max.age','Age.mat','Fecun','Max.size','Size.mat','Rep.strat','Troph.Lvl'),names(aa))]=1
    }     
    #MSC pre accreditation geometric mean formula
    PSA$Productivity[p]=mean(unlist(aa[,c('Max.age','Age.mat','Fecun',
                                          'Max.size','Size.mat','Rep.strat','Troph.Lvl')]))
    S1=1+((aa$Net.avail*aa$Net.encoun*aa$Net.sel*aa$Net.PCM)-1)/40   
    S2=1+((aa$Line.avail*aa$Line.encoun*aa$Line.sel*aa$Line.PCM)-1)/40
    S3=1+((aa$Trawl.avail*aa$Trawl.encoun*aa$Trawl.sel*aa$Trawl.PCM)-1)/40
    S4=1+((aa$Trap.avail*aa$Trap.encoun*aa$Trap.sel*aa$Trap.PCM)-1)/40
    Cum.susc=min(3,(1+((S1-1)^2+(S2-1)^2+(S3-1)^2+(S4-1)^2)^0.5))
    PSA$Susceptibility[p]=Cum.susc
    PSA$Vulnerability[p]=(PSA$Productivity[p]^2+Cum.susc^2)^0.5  #Euclidean distance
  }
  
  PSA=PSA%>%
    rename(Vulnerability.score=Vulnerability)%>%
    mutate(Species=Hmisc::capitalize(as.character(Species)),
           Vulnerability=Hmisc::capitalize(
              ifelse(Vulnerability.score<=Low.risk,'Low',
              ifelse(Vulnerability.score>Low.risk & Vulnerability.score<=medium.risk,'Medium',
              'High'))),
           Vulnerability=factor(Vulnerability,levels=c('Low','Medium','High')))%>%
    arrange(Vulnerability.score)%>%
    mutate(Fnt.size=case_when(Vulnerability=='Low'~size.low,
                              Vulnerability=='Medium'~size.med,
                              Vulnerability=='High'~size.hig),
           Species=case_when(Species=='Port jackson shark'~'Port Jackson shark',
                             TRUE~Species))
  cols=c(Low="green",Medium="yellow",High="red")
  p=ggplot(PSA,
           aes(Productivity, Susceptibility, label = Species,colour = Vulnerability, fill = Vulnerability)) +
    geom_point(shape = 21, size = 5,colour="black") + 
    geom_text_repel(aes(size=Fnt.size),show.legend  = F,segment.colour="grey75",col='black',box.padding = line.sep) + 
    scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))+ 
    xlim(1.5,3.15)+ylim(0.85,3.05)+
    theme_PA(axs.T.siz=18,axs.t.siz=12,lgT.siz=16,leg.siz=15)+
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="top",
          legend.key=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    theme(plot.caption = element_text(hjust = 0))+
    ylab('')+xlab('')+
    coord_cartesian(clip = 'off')+
    annotation_custom(segmentsGrob(c(0.80,-0.05), c(-0.05, 0.15), 
                                   c(0.15, -0.05), c(-0.05, 0.82), gp = gpar(lwd = 3,col='dodgerblue4'),
                                   arrow = arrow(length = unit(5, 'mm'))))+
    annotation_custom(textGrob("High Productivity", gp=gpar(fontsize=12, fontface="bold",col='dodgerblue4')),
                      xmin=1.5,xmax=1.5,ymin=0.63,ymax=0.63)+
    annotation_custom(textGrob("Low Productivity", gp=gpar(fontsize=12, fontface="bold",col='dodgerblue4')),
                      xmin=3.15,xmax=3.05,ymin=0.63,ymax=0.63)+
    annotation_custom(textGrob("Low Susceptibility",rot = 90, gp=gpar(fontsize=12, fontface="bold",col='dodgerblue4')),
                      xmin=1.325,xmax=1.325,ymin=0.85,ymax=0.85)+
    annotation_custom(textGrob("High Susceptibility",rot = 90, gp=gpar(fontsize=12, fontface="bold",col='dodgerblue4')),
                      xmin=1.325,xmax=1.325,ymin=3,ymax=3)
  print(p)
  ggsave(paste(Exprt,'Figure_PSA.tiff',sep='/'), width = W,height = H, dpi = 300, compression = "lzw")
  
  #Table PSA scores
  Table.PSA=d
  for(p in 1:nrow(Table.PSA))    
  {
    if(!Table.PSA[p,]$Species%in%species.meeting.criteria) #reset availability and encounterability if not meeting criteria
    {
      k=KIP%>%filter(Name==Table.PSA[p,]$Species)
      Table.PSA[p,]=Table.PSA[p,]%>%
        mutate(Net.avail=ifelse(!'net'%in%k$Gear,1,Net.avail),
               Net.encoun=ifelse(!'net'%in%k$Gear,1,Net.encoun),
               Line.avail=ifelse(!'line'%in%k$Gear,1,Line.avail),
               Line.encoun=ifelse(!'line'%in%k$Gear,1,Line.encoun),
               Trawl.avail=ifelse(!'trawl'%in%k$Gear,1,Trawl.avail),
               Trawl.encoun=ifelse(!'trawl'%in%k$Gear,1,Trawl.encoun),
               Trap.avail=ifelse(!'trap'%in%k$Gear,1,Trap.avail),
               Trap.encoun=ifelse(!'trap'%in%k$Gear,1,Trap.encoun))
    }
  }
  Table.PSA=Table.PSA%>%   
    mutate(Species=capitalize(Species))%>%
    left_join(PSA%>%
                dplyr::select(Species,Vulnerability.score),
              by='Species')%>%
    mutate(Vulnerability.score=round(Vulnerability.score,2))
  write.csv(Table.PSA,paste(Exprt,'Table S2_PSA scores.csv',sep='/'),row.names=F)
  return(PSA)
}

# Contrast catch and cpue series   ------------------------------------------------------
fn.ktch.cpue=function(ktch,cpue)
{
  if(!is.null(cpue))
  {
    fn.fig(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(unique(ktch$Name)),"/",AssessYr,
                 "/1_Inputs/Visualise data/Total catch vs cpues",sep=''),2400,2000) 
    par(las=1,oma=c(1,1,1,2))
    plot(ktch$finyear,ktch$Tonnes,type='o',pch=19,main=capitalize(unique(ktch$Name)),
         ylab='Total catch (tonnes)',xlab='Financial year')
    DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(cpue))
    if(length(DROP)>0)cpue=cpue[-DROP]
    par(new=T)
    Mxylim=max(sapply(cpue, function(x) max(x$Mean, na.rm=TRUE)))
    clfun=colorRampPalette(c('gold',"chocolate1","firebrick4"))
    CLS=clfun(length(cpue))  
    
    dummy=cpue[[1]]%>%
      dplyr::select(yr.f,Mean)
    dummy=rbind(dummy,
                data.frame(yr.f=ktch$finyear[which(!ktch$finyear%in%dummy$yr.f)],Mean=NA))%>%
      arrange(yr.f)
    with(dummy,plot(yr.f,Mean,ylab='',xlab='',ylim=c(0,Mxylim),col='transparent',pch=15,yaxt='n'))
    for(x in 1:length(cpue))
    {
      with(cpue[[x]],points(yr.f,Mean,ylab='',xlab='',col=CLS[x],pch=15))
      with(cpue[[x]],lines(yr.f,Mean,ylab='',xlab='',col=CLS[x],lty=2))
      with(cpue[[x]],segments(yr.f,LOW.CI,yr.f,UP.CI,col=CLS[x]))
    }
    axis(4,seq(0,Mxylim,length.out=3),round(seq(0,Mxylim,length.out=3),2))
    mtext('Relative catch rate',4,3,las=3)
    legend("topleft",names(cpue),pch=15,col=CLS,bty='n')
    dev.off()
    
  }
}

# Functions for applying CMSY, DB-SRA, OCOM, JABBA ------------------------------------------------------
dbsra_tweeked=function(year = NULL, catch = NULL, catchCV = NULL, 
                       catargs = list(dist = "none",low = 0, up = Inf, unit = "MT"), agemat = NULL, maxn = 25, 
                       k = list(low = 0, up = NULL, tol = 0.01, permax = 1000), 
                       b1k = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0),
                       btk = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0, refyr = NULL),
                       fmsym = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0),
                       bmsyk = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0), 
                       M = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0),
                       nsims = 10000, catchout = 0, grout = 1, 
                       graphs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
                       grargs = list(lwd = 1, cex = 1, nclasses = 20, mains = " ", cex.main = 1, cex.axis = 1, cex.lab = 1),
                       pstats = list(ol = 1, mlty = 1, mlwd = 1.5, llty = 3, llwd = 1,ulty = 3, ulwd = 1),
                       grtif = list(zoom = 4, width = 11,height = 13, pointsize = 10),
                       Export=FALSE) 
{
  if (is.null(catch)) 
    stop("No catch data")
  if (length(year) != length(catch)) 
    stop("Length of year and catch differ")
  if (btk$refyr > c(max(year) + 1)) 
    stop("refyr is beyond the max year allowed (i.e.,max(year)+1)")
  catdef = list(dist = "none", low = 0, up = Inf, unit = "MT")
  if (any(is.na(catch))) 
    stop("There are missing values in catch")
  if (any(is.na(year))) 
    stop("There are missing values in year")
  if (any(names(catargs) == "dist")) {
    if (catdef$dist != catargs$dist) 
      catdef$dist <- catargs$dist
  }
  if (any(names(catargs) == "low")) {
    if (any(catdef$low != catargs$low)) 
      catdef$low <- catargs$low
  }
  if (any(names(catargs) == "up")) {
    if (any(catdef$up != catargs$up)) 
      catdef$up <- catargs$up
  }
  if (any(names(catargs) == "unit")) {
    if (catdef$unit != catargs$unit) 
      catdef$unit <- catargs$unit
  }
  kdef = list(low = 0, up = max(catch), tol = 1, permax = 1000)
  if (any(names(k) == "low")) {
    if (kdef$low != k$low) 
      kdef$low <- k$low
  }
  if (any(names(k) == "up")) {
    if (kdef$up != k$up) 
      kdef$up <- k$up
  }
  if (any(names(k) == "tol")) {
    if (kdef$tol != k$tol) 
      kdef$tol <- k$tol
  }
  if (any(names(k) == "permax")) {
    if (kdef$permax != k$permax) 
      kdef$permax <- k$permax
  }
  fmsymdef = list(dist = "unif", low = 0, up = 1, mean = 0, 
                  sd = 0)
  if (any(names(fmsym) == "dist")) {
    if (fmsymdef$dist != fmsym$dist) 
      fmsymdef$dist <- fmsym$dist
  }
  if (any(names(fmsym) == "low")) {
    if (fmsymdef$low != fmsym$low) 
      fmsymdef$low <- fmsym$low
  }
  if (any(names(fmsym) == "up")) {
    if (fmsymdef$up != fmsym$up) 
      fmsymdef$up <- fmsym$up
  }
  if (any(names(fmsym) == "mean")) {
    if (fmsymdef$mean != fmsym$mean) 
      fmsymdef$mean <- fmsym$mean
  }
  if (any(names(fmsym) == "sd")) {
    if (fmsymdef$sd != fmsym$sd) 
      fmsymdef$sd <- fmsym$sd
  }
  b1kdef = list(dist = "unif", low = 0, up = 1, mean = 0, 
                sd = 0)
  if (any(names(b1k) == "dist")) {
    if (b1kdef$dist != b1k$dist) 
      b1kdef$dist <- b1k$dist
  }
  if (any(names(b1k) == "low")) {
    if (b1kdef$low != b1k$low) 
      b1kdef$low <- b1k$low
  }
  if (any(names(b1k) == "up")) {
    if (b1kdef$up != b1k$up) 
      b1kdef$up <- b1k$up
  }
  if (any(names(b1k) == "mean")) {
    if (b1kdef$mean != b1k$mean) 
      b1kdef$mean <- b1k$mean
  }
  if (any(names(b1k) == "sd")) {
    if (b1kdef$sd != b1k$sd) 
      b1kdef$sd <- b1k$sd
  }
  if (any(names(b1k) == "refyr")) {
    if (b1kdef$refyr != b1k$refyr) 
      b1kdef$refyr <- b1k$refyr
  }
  btkdef = list(dist = "unif", low = 0, up = 1, mean = 0, 
                sd = 0, refyr = max(year))
  if (any(names(btk) == "dist")) {
    if (btkdef$dist != btk$dist) 
      btkdef$dist <- btk$dist
  }
  if (any(names(btk) == "low")) {
    if (btkdef$low != btk$low) 
      btkdef$low <- btk$low
  }
  if (any(names(btk) == "up")) {
    if (btkdef$up != btk$up) 
      btkdef$up <- btk$up
  }
  if (any(names(btk) == "mean")) {
    if (btkdef$mean != btk$mean) 
      btkdef$mean <- btk$mean
  }
  if (any(names(btk) == "sd")) {
    if (btkdef$sd != btk$sd) 
      btkdef$sd <- btk$sd
  }
  if (any(names(btk) == "refyr")) {
    if (btkdef$refyr != btk$refyr) 
      btkdef$refyr <- btk$refyr
  }
  bmsykdef = list(dist = "unif", low = 0.5, up = 0.5, mean = 0, 
                  sd = 0)
  if (any(names(bmsyk) == "dist")) {
    if (bmsykdef$dist != bmsyk$dist) 
      bmsykdef$dist <- bmsyk$dist
  }
  if (any(names(bmsyk) == "low")) {
    if (bmsykdef$low != bmsyk$low) 
      bmsykdef$low <- bmsyk$low
  }
  if (any(names(bmsyk) == "up")) {
    if (bmsykdef$up != bmsyk$up) 
      bmsykdef$up <- bmsyk$up
  }
  if (any(names(bmsyk) == "mean")) {
    if (bmsykdef$mean != bmsyk$mean) 
      bmsykdef$mean <- bmsyk$mean
  }
  if (any(names(bmsyk) == "sd")) {
    if (bmsykdef$sd != bmsyk$sd) 
      bmsykdef$sd <- bmsyk$sd
  }
  Mdef = list(dist = "unif", low = 0.2, up = 0.2, mean = 0, 
              sd = 0)
  if (any(names(M) == "dist")) {
    if (Mdef$dist != M$dist) 
      Mdef$dist <- M$dist
  }
  if (any(names(M) == "low")) {
    if (Mdef$low != M$low) 
      Mdef$low <- M$low
  }
  if (any(names(M) == "up")) {
    if (Mdef$up != M$up) 
      Mdef$up <- M$up
  }
  if (any(names(M) == "mean")) {
    if (Mdef$mean != M$mean) 
      Mdef$mean <- M$mean
  }
  if (any(names(M) == "sd")) {
    if (Mdef$sd != M$sd) 
      Mdef$sd <- M$sd
  }
  grdef = list(lwd = 1, cex.axis = 1, cex.lab = 1, cex = 1, 
               nclasses = 20, mains = " ", cex.main = 1)
  if (any(names(grargs) == "cex.axis")) {
    if (any(grdef$cex.axis != grargs$cex.axis)) 
      grdef$cex.axis <- grargs$cex.axis
  }
  if (any(names(grargs) == "cex.lab")) {
    if (any(grdef$cex.lab != grargs$cex.lab)) 
      grdef$cex.lab <- grargs$cex.lab
  }
  if (any(names(grargs) == "cex")) {
    if (any(grdef$cex != grargs$cex)) 
      grdef$cex <- grargs$cex
  }
  if (any(names(grargs) == "nclasses")) {
    if (any(grdef$nclasses != grargs$nclasses)) 
      grdef$nclasses <- grargs$nclasses
  }
  if (any(names(grargs) == "mains")) {
    if (any(grdef$mains != grargs$mains)) 
      grdef$mains <- grargs$mains
  }
  if (any(names(grargs) == "cex.main")) {
    if (any(grdef$cex.main != grargs$cex.main)) 
      grdef$cex.main <- grargs$cex.main
  }
  if (any(names(grargs) == "lwd")) {
    if (any(grdef$lwd != grargs$lwd)) 
      grdef$lwd <- grargs$lwd
  }
  pstdef = list(ol = 1, mlty = 1, mlwd = 1.5, llty = 3, llwd = 1, 
                ulty = 3, ulwd = 1)
  if (any(names(pstats) == "ol")) {
    if (pstdef$ol != pstats$ol) 
      pstdef$ol <- pstats$ol
  }
  if (any(names(pstats) == "mlty")) {
    if (pstdef$mlty != pstats$mlty) 
      pstdef$mlty <- pstats$mlty
  }
  if (any(names(pstats) == "mlwd")) {
    if (pstdef$mlwd != pstats$mlwd) 
      pstdef$mlwd <- pstats$mlwd
  }
  if (any(names(pstats) == "llty")) {
    if (pstdef$llty != pstats$llty) 
      pstdef$llty <- pstats$llty
  }
  if (any(names(pstats) == "llwd")) {
    if (pstdef$llwd != pstats$llwd) 
      pstdef$llwd <- pstats$llwd
  }
  if (any(names(pstats) == "ulty")) {
    if (pstdef$ulty != pstats$ulty) 
      pstdef$ulty <- pstats$ulty
  }
  if (any(names(pstats) == "ulwd")) {
    if (pstdef$ulwd != pstats$ulwd) 
      pstdef$ulwd <- pstats$ulwd
  }
  tifdef = list(zoom = 4, width = 11, height = 13, pointsize = 10)
  if (any(names(grtif) == "zoom")) {
    if (tifdef$zoom != grtif$zoom) 
      tifdef$zoom <- grtif$zoom
  }
  if (any(names(grtif) == "width")) {
    if (tifdef$width != grtif$width) 
      tifdef$width <- grtif$width
  }
  if (any(names(grtif) == "height")) {
    if (tifdef$height != grtif$height) 
      tifdef$height <- grtif$height
  }
  if (any(names(grtif) == "pointsize")) {
    if (tifdef$pointsize != grtif$pointsize) 
      tifdef$pointsize <- grtif$pointsize
  }
  if (b1kdef$low < 0 | b1kdef$low > 1) 
    stop("b1k low can only range from 0 to 1")
  if (b1kdef$up < 0 | b1kdef$up > 1) 
    stop("b1k up can only range from 0 to 1")
  if (btkdef$low < 0 | btkdef$low > 1) 
    stop("btk low can only range from 0 to 1")
  if (btkdef$up < 0 | btkdef$up > 1) 
    stop("btk up can only range from 0 to 1")
  if (any(is.na(catch))) 
    stop("Missing catch values are not allowed.")
  if (catdef$dist != "none") {
    if (catdef$dist %in% c("norm", "lnorm") & any(!is.numeric(catchCV))) 
      stop("catchCV is required for resampling")
    catdata <- as.data.frame(cbind(year, catch, catchCV))
  }
  if (catdef$dist == "none") {
    catdata <- as.data.frame(cbind(year, catch))
  }
  word.tif = function(filename = "Word_Figure_%03d.tif", zoom = tifdef$zoom, 
                      width = tifdef$width, height = tifdef$height, pointsize = tifdef$pointsize, 
                      ...) {
    if (!grepl("[.]ti[f]+$", filename, ignore.case = TRUE)) 
      filename = paste0(filename, ".tif")
    tiff(filename = filename, compression = "lzw", res = 96 * 
           zoom, width = width, height = height, units = "cm", 
         pointsize = pointsize, ...)
  }
  getrandom <- function(n, spec, a = -Inf, b = Inf, mean = NULL, 
                        sd = NULL) {
    if (!spec %in% c("unif", "beta", "lnorm", "norm", "gamma")) 
      stop("Unknown distribution name.")
    if (spec %in% c("beta", "lnorm", "norm", "gamma")) {
      if (is.null(mean) | is.null(sd)) 
        stop("mean and sd must be specified.")
      if (is.na(mean) | is.na(sd)) 
        stop("mean and sd must be specified.")
    }
    if (spec == "beta") {
      if (is.null(a) | is.null(b) | a == -Inf | b == Inf) 
        stop("lower and upper limits are required for the beta distribution.")
      if (is.na(a) | is.na(b)) 
        stop("lower and upper limits are required for the beta distribution.")
    }
    if (spec == "unif") {
      if (is.na(a) | is.na(b) | is.null(a) | is.null(b)) 
        stop("unif requires specified lower and upper values")
      if (a == -Inf | b == Inf) 
        stop("unif requires specified lower and upper values")
    }
    qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
      ttt <- p
      G <- get(paste("p", spec, sep = ""), mode = "function")
      Gin <- get(paste("q", spec, sep = ""), mode = "function")
      ttt <- Gin(G(a, ...) + p * (G(b, ...) - G(a, ...)), 
                 ...)
      return(ttt)
    }
    rtrunc <- function(n, spec, a = -Inf, b = Inf, ...) {
      x <- u <- runif(n, min = 0, max = 1)
      x <- qtrunc(u, spec, a = a, b = b, ...)
      return(x)
    }
    if (spec == "beta") {
      return(rtrunc(n, "beta", a = a, b = b, shape1 = mean * 
                      (((mean * (1 - mean))/sd^2) - 1), shape2 = (1 - 
                                                                    mean) * (((mean * (1 - mean))/sd^2) - 1)))
    }
    if (spec == "unif") 
      return(rtrunc(n, "unif", min = a, max = b))
    if (spec == "norm") 
      return(rtrunc(n, "norm", a = a, b = b, mean = mean, 
                    sd = sd))
    if (spec == "lnorm") 
      return(rtrunc(n, "lnorm", a = a, b = b, meanlog = mean, 
                    sdlog = sd))
    if (spec == "gamma") 
      return(rtrunc(n, "gamma", a = a, b = b, shape = (mean/sd)^2, 
                    scale = sd^2/mean))
  }
  if (catdef$dist == "unif") {
    if (length(catdef$low) != length(catdata[, 2]) | length(catdef$up) != 
        length(catdata[, 2])) 
      stop("The length of catargs$low and/or catargs$up should be the same length as catch")
  }
  timelen <- length(year)
  refyr <- which(btk$refyr == c(year, year[length(year)] + 
                                  1))

  fmsymX <- fmsymdef$mean
  b1kX <- b1kdef$mean
  btkX <- btkdef$mean
  bmsykX <- bmsykdef$mean
  MM <- Mdef$mean
  timeseries_B=timeseries_BK=timeseries_BBmsy=timeseries_FFmsy=vector('list',nsims)
  f <- function(d) {
    (bmsykX - (d^(1/(1 - d))))^2
  }
  findk <- function(K) {
    B <- NULL
    P <- NULL
    Fmsy <- fmsymX * MM
    Umsy <- (Fmsy/(Fmsy + MM)) * (1 - exp(-Fmsy - MM))
    MSY <- K * bmsykX * Umsy
    B[1] <- b1kX * K
    for (t in 1:timelen) {
      if (t <= agemat) 
        P[t] <- 0
      if (t > agemat) {
        if (bmsykX >= 0.5) 
          P[t] <- g * MSY * (B[t - agemat]/K) - g * 
            MSY * (B[t - agemat]/K)^n
        if (bmsykX > 0.3 & bmsykX < 0.5) {
          bjoin <- (0.75 * bmsykX - 0.075) * K
          if (B[t - agemat] < bjoin) {
            PJ <- g * MSY * (bjoin/K) - g * MSY * 
              (bjoin/K)^n
            cc <- (1 - n) * MSY * g * (bjoin^(n - 
                                                2)) * K^-n
            P[t] <- B[t - agemat] * (PJ/bjoin + cc * 
                                       (B[t - agemat] - bjoin))
          }
          if (B[t - agemat] >= bjoin) 
            P[t] <- g * MSY * (B[t - agemat]/K) - 
              g * MSY * (B[t - agemat]/K)^n
        }
        if (bmsykX <= 0.3) {
          bjoin <- (0.5 * bmsykX) * K
          if (B[t - agemat] < bjoin) {
            PJ <- g * MSY * (bjoin/K) - g * MSY * 
              (bjoin/K)^n
            cc <- (1 - n) * MSY * g * (bjoin^(n - 
                                                2)) * K^-n
            P[t] <- B[t - agemat] * (PJ/bjoin + cc * 
                                       (B[t - agemat] - bjoin))
          }
          if (B[t - agemat] >= bjoin) 
            P[t] <- g * MSY * (B[t - agemat]/K) - 
              g * MSY * (B[t - agemat]/K)^n
        }
      }
      if (P[t] < 0) 
        P[t] <- 0
      B[t + 1] <- max(0, B[t] + P[t] - dcatch[t])
    }
    (btkX - (B[refyr]/K))^2
  }
  pb <- txtProgressBar(max = nsims, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  dummy=foreach(nn= 1:nsims,.options.snow = opts) %dopar%
  {
      if (fmsymdef$dist != "none") 
        fmsymX <- getrandom(1, fmsymdef$dist, a = fmsymdef$low, 
                            b = fmsymdef$up, mean = fmsymdef$mean, sd = fmsymdef$sd)
      if (btkdef$dist != "none") 
        btkX <- getrandom(1, btkdef$dist, a = btkdef$low, 
                          b = btkdef$up, mean = btkdef$mean, sd = btkdef$sd)
      if (b1kdef$dist != "none") 
        b1kX <- getrandom(1, b1kdef$dist, a = b1kdef$low, 
                          b = b1kdef$up, mean = b1kdef$mean, sd = b1kdef$sd)
      if (bmsykdef$dist != "none") 
        bmsykX <- getrandom(1, bmsykdef$dist, a = bmsykdef$low, 
                            b = bmsykdef$up, mean = bmsykdef$mean, sd = bmsykdef$sd)
      if (Mdef$dist != "none") 
        MM <- getrandom(1, Mdef$dist, a = Mdef$low, b = Mdef$up, 
                        mean = Mdef$mean, sd = Mdef$sd)
      dcatch <- NULL
      if (catdef$dist == "none") {
        dcatch <- catdata[, 2]
        catchout <- 0
      }
      if (catdef$dist != "none") {
        for (cc in 1:c(length(catdata[, 2]))) {
          if (catdef$dist == "unif") 
            dcatch[cc] <- getrandom(1, catdef$dist, a = catdef$low[cc], 
                                    b = catdef$up[cc])
          if (catdef$dist == "norm") 
            dcatch[cc] <- getrandom(1, catdef$dist, a = catdef$low, 
                                    b = catdef$up, mean = catdata[cc, 2], sd = catdata[cc, 
                                                                                       2] * catdata[cc, 3])
          if (catdef$dist == "lnorm") 
            dcatch[cc] <- getrandom(1, catdef$dist, a = log(catdef$low), 
                                    b = log(catdef$up), mean = log(catdata[cc, 
                                                                           2]), sd = sqrt(log(catdata[cc, 3]^2 + 
                                                                                                1)))
        }
        dcatch <- ifelse(dcatch < 0 | is.infinite(dcatch), 
                         0, dcatch)
      }
      
      n <- optimize(f, c(0, 1000), tol = 1e-06)[[1]]
      g <- (n^(n/(n - 1)))/(n - 1)
      
      out <- optimize(findk, c(kdef$low, kdef$up), tol = kdef$tol)
      bigK <- out$minimum
      B <- NULL
      P <- NULL
      Fmsy <- fmsymX * MM
      Umsy <- (Fmsy/(Fmsy + MM)) * (1 - exp(-Fmsy - MM))
      MSY <- bigK * bmsykX * Umsy
      B[1] <- bigK * b1kX
      for (t in 1:timelen) {
        if (t <= agemat) 
          P[t] <- 0
        if (t > agemat) {
          if (bmsykX >= 0.5) 
            P[t] <- g * MSY * (B[t - agemat]/bigK) - g * 
              MSY * (B[t - agemat]/bigK)^n
          if (bmsykX > 0.3 & bmsykX < 0.5) {
            bjoin <- (0.75 * bmsykX - 0.075) * bigK
            if (B[t - agemat] < bjoin) {
              PJ <- g * MSY * (bjoin/bigK) - g * MSY * 
                (bjoin/bigK)^n
              cc <- (1 - n) * MSY * g * (bjoin^(n - 2)) * 
                bigK^-n
              P[t] <- B[t - agemat] * (PJ/bjoin + cc * 
                                         (B[t - agemat] - bjoin))
            }
            if (B[t - agemat] >= bjoin) 
              P[t] <- g * MSY * (B[t - agemat]/bigK) - 
                g * MSY * (B[t - agemat]/bigK)^n
          }
          if (bmsykX <= 0.3) {
            bjoin <- (0.5 * bmsykX) * bigK
            if (B[t - agemat] < bjoin) {
              PJ <- g * MSY * (bjoin/bigK) - g * MSY * 
                (bjoin/bigK)^n
              cc <- (1 - n) * MSY * g * (bjoin^(n - 2)) * 
                bigK^-n
              P[t] <- B[t - agemat] * (PJ/bjoin + cc * 
                                         (B[t - agemat] - bjoin))
            }
            if (B[t - agemat] >= bjoin) 
              P[t] <- g * MSY * (B[t - agemat]/bigK) - 
                g * MSY * (B[t - agemat]/bigK)^n
          }
        }
        if (P[t] < 0) 
          P[t] <- 0
        B[t + 1] <- max(0, B[t] + P[t] - dcatch[t])
      }
      bll <- 0
      if (min(B) > 0 && max(B) <= bigK && (out$objective <= 
                                           kdef$tol^2) && (abs((max(B) - bigK)/bigK) * 100) <= 
          kdef$permax && n <= maxn) 
        bll <- 1
      if(Export)
      {
        if (nn == 1) {
          write.table(t(c(bll, B)), file = "Biotraj-dbsra.csv", 
                      sep = ",", row.names = FALSE, col.names = FALSE, 
                      append = FALSE)
          if (catchout == 1) 
            write.table(t(c(bll, dcatch)), file = "Catchtraj-dbsra.csv", 
                        sep = ",", row.names = FALSE, col.names = FALSE, 
                        append = FALSE)
        }
        if (nn > 1) {
          write.table(t(c(bll, B)), file = "Biotraj-dbsra.csv", 
                      sep = ",", row.names = FALSE, col.names = FALSE, 
                      append = TRUE)
          if (catchout == 1) 
            write.table(t(c(bll, dcatch)), file = "Catchtraj-dbsra.csv", 
                        sep = ",", row.names = FALSE, col.names = FALSE, 
                        append = TRUE)
        }
      }

      storep=rep(NA,16)
      names(storep) <- c("ll", "FmsyM", "BtK", "BmsyK", "M", 
                         "K", "Fmsy", "Umsy", "MSY", "Bmsy", "OFLT1", "Brefyr", 
                         "BT1", "n", "g", "B1K")
      storep[1] <- bll
      storep[2] <- fmsymX
      storep[3] <- btkX
      storep[4] <- bmsykX
      storep[5] <- MM
      storep[6] <- bigK
      storep[7] <- Fmsy
      storep[8] <- Umsy
      storep[9] <- MSY
      storep[10] <- bigK * bmsykX
      storep[11] <- Umsy * B[timelen + 1]
      storep[12] <- B[refyr]
      storep[13] <- B[timelen + 1]
      storep[14] <- n
      storep[15] <- g
      storep[16] <- b1kX

      timeseries_B=B
      timeseries_BK=B/bigK
      timeseries_BBmsy=B/(bigK * bmsykX)
      FF=-log(1-(catdata$catch/B[-1]))
      timeseries_FFmsy=FF/Fmsy
      return(list(storep=storep,B=B,timeseries_BK=timeseries_BK,
                  timeseries_BBmsy=timeseries_BBmsy,FF=FF,timeseries_FFmsy=timeseries_FFmsy))
    }   #end nn loop
  close(pb)
  storep=data.frame(do.call(rbind,fn.get.stuff.from.list(dummy,'storep')))
  timeseries_B=data.frame(do.call(rbind,fn.get.stuff.from.list(dummy,'B')))
  timeseries_BK=data.frame(do.call(rbind,fn.get.stuff.from.list(dummy,'timeseries_BK')))
  timeseries_BBmsy=data.frame(do.call(rbind,fn.get.stuff.from.list(dummy,'timeseries_BBmsy')))
  timeseries_FF=data.frame(do.call(rbind,fn.get.stuff.from.list(dummy,'FF')))
  timeseries_FFmsy=data.frame(do.call(rbind,fn.get.stuff.from.list(dummy,'timeseries_FFmsy')))
  
  iiid=which(storep$ll==1)
  timeseries_B=timeseries_B[iiid,]
  timeseries_BK=timeseries_BK[iiid,]
  timeseries_BBmsy=timeseries_BBmsy[iiid,]
  timeseries_FFmsy=timeseries_FFmsy[iiid,]
  timeseries_FF=timeseries_FF[iiid,]
  
  datar <- storep[storep[, 1] == 1, ]
  if (length(datar[, 1]) > 0) {
    mMSY <- round(mean(datar$MSY), 3)
    MSY95 <- round(quantile(datar$MSY, probs = c(0.025, 
                                                 0.5, 0.975)), 4)
    mk <- round(mean(datar$K), 3)
    k95 <- round(quantile(datar$K, probs = c(0.025, 0.5, 
                                             0.975)), 4)
    mBMSY <- round(mean(datar$Bmsy), 3)
    BMSY95 <- round(quantile(datar$Bmsy, probs = c(0.025, 
                                                   0.5, 0.975)), 4)
    mFMSY <- round(mean(datar$Fmsy), 3)
    FMSY95 <- round(quantile(datar$Fmsy, probs = c(0.025, 
                                                   0.5, 0.975)), 4)
    mUMSY <- round(mean(datar$Umsy), 3)
    UMSY95 <- round(quantile(datar$Umsy, probs = c(0.025, 
                                                   0.5, 0.975)), 4)
    mOFL <- round(mean(datar$OFLT1), 3)
    OFL95 <- round(quantile(datar$OFLT1, probs = c(0.025, 
                                                   0.5, 0.975)), 4)
    mM <- round(mean(datar$M), 3)
    M95 <- round(quantile(datar$M, probs = c(0.025, 0.5, 
                                             0.975)), 4)
    mbtk <- round(mean(datar$BtK), 3)
    btk95 <- round(quantile(datar$BtK, probs = c(0.025, 
                                                 0.5, 0.975)), 4)
    mb1k <- round(mean(datar$B1K), 3)
    b1k95 <- round(quantile(datar$B1K, probs = c(0.025, 
                                                 0.5, 0.975)), 4)
    mFmsyM <- round(mean(datar$FmsyM), 3)
    FmsyM95 <- round(quantile(datar$FmsyM, probs = c(0.025, 
                                                     0.5, 0.975)), 4)
    mBmsyK <- round(mean(datar$BmsyK), 3)
    BmsyK95 <- round(quantile(datar$BmsyK, probs = c(0.025, 
                                                     0.5, 0.975)), 4)
    mBrefyr <- round(mean(datar$Brefyr), 3)
    Brefyr95 <- round(quantile(datar$Brefyr, probs = c(0.025, 
                                                       0.5, 0.975)), 4)
    if (grout > 0) {
      grunits <- data.frame(gr = c(1:15), cexa = 0, cexl = 0, 
                            cexx = 0, nclass = 0, mains = " ", cexmain = 0, 
                            lwd = 0, stringsAsFactors = FALSE)
      grunits$lwd <- ifelse(grunits$gr %in% c(1, 13), 
                            grdef$lwd, 0)
      grunits[, 2] <- grdef$cex.axis
      grunits[, 3] <- grdef$cex.lab
      grunits[, 4] <- grdef$cex
      grunits$nclass <- ifelse(grunits$gr %in% c(2, 3, 
                                                 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15), grdef$nclasses, 
                               grunits$nclass)
      grunits$cexmain <- grdef$cex.main
      grunits[graphs, 6] <- grdef$mains
      if (any(graphs == 1)) {
        plot(catch ~ year, type = "l", xlab = "Year", 
             ylab = paste("Catch (", catdef$unit, ")", 
                          sep = ""), ylim = c(0, round(max(catch, 
                                                           mMSY, MSY95[3]), 1)), cex = grunits[1, 4], 
             cex.lab = grunits[1, 3], cex.axis = grunits[1, 
                                                         2], main = grunits[1, 6], cex.main = grunits[1, 
                                                                                                      7], lwd = grunits[1, 8])
        if (pstdef$ol == 1) {
          abline(h = MSY95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(h = MSY95[1], lwd = pstdef$llwd, lty = pstdef$llty)
          abline(h = MSY95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("catch")
          plot(catch ~ year, type = "l", xlab = "Year", 
               ylab = paste("Catch (", catdef$unit, ")", 
                            sep = ""), ylim = c(0, round(max(catch, 
                                                             mMSY, MSY95[3]), 1)), cex = grunits[1, 
                                                                                                 4], cex.lab = grunits[1, 3], cex.axis = grunits[1, 
                                                                                                                                                 2], main = grunits[1, 6], cex.main = grunits[1, 
                                                                                                                                                                                              7], lwd = grunits[1, 8])
          if (pstdef$ol == 1) {
            abline(h = MSY95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(h = MSY95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(h = MSY95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 2)) {
        hist(datar$K, freq = FALSE, xlim = c(0, max(datar$K) * 
                                               1.4), xlab = paste("K (", catdef$unit, ")", 
                                                                  sep = ""), nclass = grunits[3, 5], cex.lab = grunits[3, 
                                                                                                                       3], cex.axis = grunits[3, 2], main = grunits[3, 
                                                                                                                                                                    6], cex.main = grunits[3, 7])
        if (pstdef$ol == 1) {
          abline(v = k95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(v = k95[1], lwd = pstdef$lwd, lty = pstdef$llty)
          abline(v = k95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("Kden")
          hist(datar$K, freq = FALSE, xlim = c(0, max(datar$K) * 
                                                 1.4), xlab = paste("K (", catdef$unit, ")", 
                                                                    sep = ""), nclass = grunits[3, 5], cex.lab = grunits[3, 
                                                                                                                         3], cex.axis = grunits[3, 2], main = grunits[3, 
                                                                                                                                                                      6], cex.main = grunits[3, 7])
          if (pstdef$ol == 1) {
            abline(v = k95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
            abline(v = k95[1], lwd = pstdef$lwd, lty = pstdef$llty)
            abline(v = k95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 3)) {
        hist(datar$Bmsy, freq = FALSE, xlim = c(0, max(datar$Bmsy) * 
                                                  1.4), xlab = paste("Bmsy (", catdef$unit, 
                                                                     ")", sep = ""), nclass = grunits[4, 5], cex.lab = grunits[4, 
                                                                                                                               3], cex.axis = grunits[4, 2], main = grunits[4, 
                                                                                                                                                                            6], cex.main = grunits[4, 7])
        if (pstdef$ol == 1) {
          abline(v = BMSY95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(v = BMSY95[1], lwd = pstdef$llwd, lty = pstdef$llty)
          abline(v = BMSY95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("Bmsyden")
          hist(datar$Bmsy, freq = FALSE, xlim = c(0, 
                                                  max(datar$Bmsy) * 1.4), xlab = paste("Bmsy (", 
                                                                                       catdef$unit, ")", sep = ""), nclass = grunits[4, 
                                                                                                                                     5], cex.lab = grunits[4, 3], cex.axis = grunits[4, 
                                                                                                                                                                                     2], main = grunits[4, 6], cex.main = grunits[4, 
                                                                                                                                                                                                                                  7])
          if (pstdef$ol == 1) {
            abline(v = BMSY95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = BMSY95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = BMSY95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 4)) {
        hist(datar$MSY, freq = FALSE, xlim = c(min(datar$MSY) * 
                                                 0.8, max(datar$MSY) * 1.2), xlab = paste("MSY (", 
                                                                                          catdef$unit, ")", sep = ""), nclass = grunits[5, 
                                                                                                                                        5], cex.lab = grunits[5, 3], cex.axis = grunits[5, 
                                                                                                                                                                                        2], main = grunits[5, 6], cex.main = grunits[5, 
                                                                                                                                                                                                                                     7])
        if (pstdef$ol == 1) {
          abline(v = MSY95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(v = MSY95[1], lwd = pstdef$llwd, lty = pstdef$llty)
          abline(v = MSY95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("MSYden")
          hist(datar$MSY, freq = FALSE, xlim = c(min(datar$MSY) * 
                                                   0.8, max(datar$MSY) * 1.2), xlab = paste("MSY (", 
                                                                                            catdef$unit, ")", sep = ""), nclass = grunits[5, 
                                                                                                                                          5], cex.lab = grunits[5, 3], cex.axis = grunits[5, 
                                                                                                                                                                                          2], main = grunits[5, 6], cex.main = grunits[5, 
                                                                                                                                                                                                                                       7])
          if (pstdef$ol == 1) {
            abline(v = MSY95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = MSY95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = MSY95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 5)) {
        hist(datar$Fmsy, freq = FALSE, xlim = c(0, max(datar$Fmsy) * 
                                                  1.3), xlab = "Fmsy", nclass = grunits[6, 5], 
             cex.lab = grunits[6, 3], cex.axis = grunits[6, 
                                                         2], main = grunits[6, 6], cex.main = grunits[6, 
                                                                                                      7])
        if (pstdef$ol == 1) {
          abline(v = FMSY95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(v = FMSY95[1], lwd = pstdef$llwd, lty = pstdef$llty)
          abline(v = FMSY95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("Fmsyden")
          hist(datar$Fmsy, freq = FALSE, xlim = c(0, 
                                                  max(datar$Fmsy) * 1.2), xlab = "Fmsy", nclass = grunits[6, 
                                                                                                          5], cex.lab = grunits[6, 3], cex.axis = grunits[6, 
                                                                                                                                                          2], main = grunits[6, 6], cex.main = grunits[6, 
                                                                                                                                                                                                       7])
          if (pstdef$ol == 1) {
            abline(v = FMSY95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = FMSY95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = FMSY95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 6)) {
        hist(datar$Umsy, freq = FALSE, xlim = c(0, max(datar$Umsy) * 
                                                  1.2), xlab = "Umsy", nclass = grunits[7, 5], 
             cex.lab = grunits[7, 3], cex.axis = grunits[7, 
                                                         2], main = grunits[7, 6], cex.main = grunits[7, 
                                                                                                      7])
        if (pstdef$ol == 1) {
          abline(v = UMSY95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(v = UMSY95[1], lwd = pstdef$llwd, lty = pstdef$llty)
          abline(v = UMSY95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("Umsyden")
          hist(datar$Umsy, freq = FALSE, xlim = c(0, 
                                                  max(datar$Umsy) * 1.2), xlab = "Umsy", nclass = grunits[7, 
                                                                                                          5], cex.lab = grunits[7, 3], cex.axis = grunits[7, 
                                                                                                                                                          2], main = grunits[7, 6], cex.main = grunits[7, 
                                                                                                                                                                                                       7])
          if (pstdef$ol == 1) {
            abline(v = UMSY95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = UMSY95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = UMSY95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 7)) {
        hist(datar$OFL, freq = FALSE, xlim = c(min(datar$OFL) * 
                                                 0.8, max(datar$OFL) * 1.3), xlab = paste("OFL (", 
                                                                                          catdef$unit, ")", sep = ""), nclass = grunits[8, 
                                                                                                                                        5], cex.lab = grunits[8, 3], cex.axis = grunits[8, 
                                                                                                                                                                                        2], main = grunits[8, 6], cex.main = grunits[8, 
                                                                                                                                                                                                                                     7])
        if (pstdef$ol == 1) {
          abline(v = OFL95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(v = OFL95[1], lwd = pstdef$llwd, lty = pstdef$llty)
          abline(v = OFL95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("OFLden")
          hist(datar$OFL, freq = FALSE, xlim = c(min(datar$OFL) * 
                                                   0.8, max(datar$OFL) * 1.3), xlab = paste("OFL (", 
                                                                                            catdef$unit, ")", sep = ""), nclass = grunits[8, 
                                                                                                                                          5], cex.lab = grunits[8, 3], cex.axis = grunits[8, 
                                                                                                                                                                                          2], main = grunits[8, 6], cex.main = grunits[8, 
                                                                                                                                                                                                                                       7])
          if (pstdef$ol == 1) {
            abline(v = OFL95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = OFL95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = OFL95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 8)) {
        hist(datar$M, freq = FALSE, xlim = c(0, max(datar$M) * 
                                               1.2), xlab = "M", nclass = grunits[9, 5], 
             cex.lab = grunits[9, 3], cex.axis = grunits[9, 
                                                         2], main = grunits[9, 6], cex.main = grunits[9, 
                                                                                                      7])
        if (pstdef$ol == 1) {
          abline(v = M95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(v = M95[1], lwd = pstdef$llwd, lty = pstdef$llty)
          abline(v = M95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("Mden")
          hist(datar$M, freq = FALSE, xlim = c(0, max(datar$M) * 
                                                 1.2), xlab = "M", nclass = grunits[9, 5], 
               cex.lab = grunits[9, 3], cex.axis = grunits[9, 
                                                           2], main = grunits[9, 6], cex.main = grunits[9, 
                                                                                                        7])
          if (pstdef$ol == 1) {
            abline(v = M95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
            abline(v = M95[1], lwd = pstdef$llwd, lty = pstdef$llty)
            abline(v = M95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 9)) {
        hist(datar$BtK, freq = FALSE, xlim = c(0, max(datar$BtK) * 
                                                 1.3), xlab = "Bt/K", nclass = grunits[10, 
                                                                                       5], cex.lab = grunits[10, 3], cex.axis = grunits[10, 
                                                                                                                                        2], main = grunits[10, 6], cex.main = grunits[10, 
                                                                                                                                                                                      7])
        if (pstdef$ol == 1) {
          abline(v = btk95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(v = btk95[1], lwd = pstdef$llwd, lty = pstdef$llty)
          abline(v = btk95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("BtKden")
          hist(datar$BtK, freq = FALSE, xlim = c(0, 
                                                 max(datar$BtK) * 1.3), xlab = "Bt/K", nclass = grunits[10, 
                                                                                                        5], cex.lab = grunits[10, 3], cex.axis = grunits[10, 
                                                                                                                                                         2], main = grunits[10, 6], cex.main = grunits[10, 
                                                                                                                                                                                                       7])
          if (pstdef$ol == 1) {
            abline(v = btk95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = btk95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = btk95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 15)) {
        hist(datar$B1K, freq = FALSE, xlim = c(0, max(datar$B1K) * 
                                                 1.3), xlab = "B1/K", nclass = grunits[15, 
                                                                                       5], cex.lab = grunits[15, 3], cex.axis = grunits[15, 
                                                                                                                                        2], main = grunits[15, 6], cex.main = grunits[15, 
                                                                                                                                                                                      7])
        if (pstdef$ol == 1) {
          abline(v = b1k95[2], lwd = pstdef$mlwd, lty = pstdef$mlty)
          abline(v = b1k95[1], lwd = pstdef$llwd, lty = pstdef$llty)
          abline(v = b1k95[3], lwd = pstdef$ulwd, lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("B1Kden")
          hist(datar$B1K, freq = FALSE, xlim = c(0, 
                                                 max(datar$B1K) * 1.3), xlab = "B1/K", nclass = grunits[15, 
                                                                                                        5], cex.lab = grunits[15, 3], cex.axis = grunits[15, 
                                                                                                                                                         2], main = grunits[15, 6], cex.main = grunits[15, 
                                                                                                                                                                                                       7])
          if (pstdef$ol == 1) {
            abline(v = b1k95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = b1k95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = b1k95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 10)) {
        hist(datar$FmsyM, freq = FALSE, xlim = c(0, 
                                                 max(datar$FmsyM) * 1.3), xlab = "Fmsy/M", 
             nclass = grunits[11, 5], cex.lab = grunits[11, 
                                                        3], cex.axis = grunits[11, 2], main = grunits[11, 
                                                                                                      6], cex.main = grunits[11, 7])
        if (pstdef$ol == 1) {
          abline(v = FmsyM95[2], lwd = pstdef$mlwd, 
                 lty = pstdef$mlty)
          abline(v = FmsyM95[1], lwd = pstdef$llwd, 
                 lty = pstdef$llty)
          abline(v = FmsyM95[3], lwd = pstdef$ulwd, 
                 lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("FmsyMden")
          hist(datar$FmsyM, freq = FALSE, xlim = c(0, 
                                                   max(datar$FmsyM) * 1.3), xlab = "Fmsy/M", 
               nclass = grunits[11, 5], cex.lab = grunits[11, 
                                                          3], cex.axis = grunits[11, 2], main = grunits[11, 
                                                                                                        6], cex.main = grunits[11, 7])
          if (pstdef$ol == 1) {
            abline(v = FmsyM95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = FmsyM95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = FmsyM95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 11)) {
        hist(datar$BmsyK, freq = FALSE, xlim = c(0, 
                                                 max(datar$BmsyK) * 1.3), xlab = "Bmsy/K", 
             nclass = grunits[11, 5], cex.lab = grunits[11, 
                                                        3], cex.axis = grunits[11, 2], main = grunits[11, 
                                                                                                      6], cex.main = grunits[11, 7])
        if (pstdef$ol == 1) {
          abline(v = BmsyK95[2], lwd = pstdef$mlwd, 
                 lty = pstdef$mlty)
          abline(v = BmsyK95[1], lwd = pstdef$llwd, 
                 lty = pstdef$llty)
          abline(v = BmsyK95[3], lwd = pstdef$ulwd, 
                 lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("BmsyKden")
          hist(datar$BmsyK, freq = FALSE, xlim = c(0, 
                                                   max(datar$BmsyK) * 1.3), xlab = "Bmsy/K", 
               nclass = grunits[11, 5], cex.lab = grunits[11, 
                                                          3], cex.axis = grunits[11, 2], main = grunits[11, 
                                                                                                        6], cex.main = grunits[11, 7])
          if (pstdef$ol == 1) {
            abline(v = BmsyK95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = BmsyK95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = BmsyK95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 12)) {
        hist(datar$Brefyr, freq = FALSE, xlim = c(min(datar$Brefyr) * 
                                                    0.7, max(datar$Brefyr) * 1.3), xlab = "Brefyr", 
             nclass = grunits[11, 5], cex.lab = grunits[11, 
                                                        3], cex.axis = grunits[11, 2], main = grunits[11, 
                                                                                                      6], cex.main = grunits[11, 7])
        if (pstdef$ol == 1) {
          abline(v = Brefyr95[2], lwd = pstdef$mlwd, 
                 lty = pstdef$mlty)
          abline(v = Brefyr95[1], lwd = pstdef$llwd, 
                 lty = pstdef$llty)
          abline(v = Brefyr95[3], lwd = pstdef$ulwd, 
                 lty = pstdef$ulty)
        }
        if (grout == 2) {
          word.tif("Brefyrden")
          hist(datar$Brefyr, freq = FALSE, xlim = c(min(datar$Brefyr) * 
                                                      0.7, max(datar$Brefyr) * 1.3), xlab = "Brefyr", 
               nclass = grunits[11, 5], cex.lab = grunits[11, 
                                                          3], cex.axis = grunits[11, 2], main = grunits[11, 
                                                                                                        6], cex.main = grunits[11, 7])
          if (pstdef$ol == 1) {
            abline(v = Brefyr95[2], lwd = pstdef$mlwd, 
                   lty = pstdef$mlty)
            abline(v = Brefyr95[1], lwd = pstdef$llwd, 
                   lty = pstdef$llty)
            abline(v = Brefyr95[3], lwd = pstdef$ulwd, 
                   lty = pstdef$ulty)
          }
          dev.off()
        }
      }
      if (any(graphs == 13) & Export) {
        bigB <- read.csv("Biotraj-dbsra.csv", header = FALSE)
        ar <- bigB[, 1]
        bigB <- t(bigB[, -1])
        cols <- ifelse(ar == 0, "gray85", "black")
        types <- ifelse(ar == 0, 3, 1)
        par(mfrow = c(1, 2))
        matplot(y = bigB[, c(which(ar == 1))], x = c(year, 
                                                     year[length(year)] + 1), type = "l", lty = types[c(which(ar == 
                                                                                                                1))], xlab = "Year", ylab = paste("Biomass (", 
                                                                                                                                                  catdef$unit, ")", sep = ""), ylim = c(0, max(bigB)), 
                cex = grunits[13, 4], cex.lab = grunits[13, 
                                                        3], col = cols[c(which(ar == 1))], main = "Accepted", 
                cex.axis = grunits[13, 2], cex.main = grunits[13, 
                                                              7], lwd = grunits[13, 8])
        lines(y = apply(bigB[, c(which(ar == 1))], 1, 
                        median), x = c(year, year[length(year)] + 
                                         1), lwd = 2, lty = pstdef$mlty, col = "red")
        lines(y = apply(bigB[, c(which(ar == 1))], 1, 
                        function(x) {
                          quantile(x, probs = 0.025)
                        }), x = c(year, year[length(year)] + 1), lwd = 2, 
              lty = pstdef$llty, col = "red")
        lines(y = apply(bigB[, c(which(ar == 1))], 1, 
                        function(x) {
                          quantile(x, probs = 0.975)
                        }), x = c(year, year[length(year)] + 1), lwd = 2, 
              lty = pstdef$ulty, col = "red")
        if (length(c(which(ar == 0))) > 2) {
          matplot(y = bigB[, c(which(ar == 0))], x = c(year, 
                                                       year[length(year)] + 1), type = "l", lty = types[c(which(ar == 
                                                                                                                  0))], xlab = "Year", ylab = paste("Biomass (", 
                                                                                                                                                    catdef$unit, ")", sep = ""), ylim = c(0, 
                                                                                                                                                                                          max(bigB)), cex = grunits[13, 4], cex.lab = grunits[13, 
                                                                                                                                                                                                                                              3], col = cols[c(which(ar == 0))], main = "Rejected", 
                  cex.axis = grunits[13, 2], cex.main = grunits[13, 
                                                                7], lwd = grunits[13, 8])
          lines(y = apply(bigB[, c(which(ar == 0))], 
                          1, median), x = c(year, year[length(year)] + 
                                              1), lwd = 2, lty = pstdef$mlty, col = "red")
          lines(y = apply(bigB[, c(which(ar == 0))], 
                          1, function(x) {
                            quantile(x, probs = 0.025)
                          }), x = c(year, year[length(year)] + 1), 
                lwd = 2, lty = pstdef$llty, col = "red")
          lines(y = apply(bigB[, c(which(ar == 0))], 
                          1, function(x) {
                            quantile(x, probs = 0.975)
                          }), x = c(year, year[length(year)] + 1), 
                lwd = 2, lty = pstdef$ulty, col = "red")
        }
        if (length(c(which(ar == 0))) <= 2) {
          warning("<3 runs were rejected!")
        }
        if (grout == 2) {
          word.tif("Biomasstraj")
          par(mfrow = c(1, 2))
          matplot(y = bigB[, c(which(ar == 1))], x = c(year, 
                                                       year[length(year)] + 1), type = "l", lty = types[c(which(ar == 
                                                                                                                  1))], xlab = "Year", ylab = paste("Biomass (", 
                                                                                                                                                    catdef$unit, ")", sep = ""), ylim = c(0, 
                                                                                                                                                                                          max(bigB)), cex = grunits[13, 4], cex.lab = grunits[13, 
                                                                                                                                                                                                                                              3], col = cols[c(which(ar == 1))], main = "Accepted", 
                  cex.axis = grunits[13, 2], cex.main = grunits[13, 
                                                                7], lwd = grunits[13, 8])
          lines(y = apply(bigB[, c(which(ar == 1))], 
                          1, median), x = c(year, year[length(year)] + 
                                              1), lwd = 2, lty = pstdef$mlty, col = "red")
          lines(y = apply(bigB[, c(which(ar == 1))], 
                          1, function(x) {
                            quantile(x, probs = 0.025)
                          }), x = c(year, year[length(year)] + 1), 
                lwd = 2, lty = pstdef$llty, col = "red")
          lines(y = apply(bigB[, c(which(ar == 1))], 
                          1, function(x) {
                            quantile(x, probs = 0.975)
                          }), x = c(year, year[length(year)] + 1), 
                lwd = 2, lty = pstdef$ulty, col = "red")
          if (length(c(which(ar == 0))) > 2) {
            matplot(y = bigB[, c(which(ar == 0))], x = c(year, 
                                                         year[length(year)] + 1), type = "l", lty = types[c(which(ar == 
                                                                                                                    0))], xlab = "Year", ylab = paste("Biomass (", 
                                                                                                                                                      catdef$unit, ")", sep = ""), ylim = c(0, 
                                                                                                                                                                                            max(bigB)), cex = grunits[13, 4], cex.lab = grunits[13, 
                                                                                                                                                                                                                                                3], col = cols[c(which(ar == 0))], main = "Rejected", 
                    cex.axis = grunits[13, 2], cex.main = grunits[13, 
                                                                  7], lwd = grunits[13, 8])
            lines(y = apply(bigB[, c(which(ar == 0))], 
                            1, median), x = c(year, year[length(year)] + 
                                                1), lwd = 2, lty = pstdef$mlty, col = "red")
            lines(y = apply(bigB[, c(which(ar == 0))], 
                            1, function(x) {
                              quantile(x, probs = 0.025)
                            }), x = c(year, year[length(year)] + 1), 
                  lwd = 2, lty = pstdef$llty, col = "red")
            lines(y = apply(bigB[, c(which(ar == 0))], 
                            1, function(x) {
                              quantile(x, probs = 0.975)
                            }), x = c(year, year[length(year)] + 1), 
                  lwd = 2, lty = pstdef$ulty, col = "red")
          }
          dev.off()
        }
        rm(bigB)
        par(mfrow = c(1, 1))
      }
      if (any(graphs == 14)) {
        breaks <- hist(storep[, "BmsyK"], plot = FALSE, 
                       nclass = grunits[14, 5])$breaks
        good <- table((cut(storep[storep[, 1] == 1, 
                                  "BmsyK"], breaks)))
        dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                             1)] + breaks[2:c(length(breaks))])/2
        bad <- table((cut(storep[storep[, 1] == 0, "BmsyK"], 
                          breaks)))
        dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                            1)] + breaks[2:c(length(breaks))])/2
        newd <- cbind(good, bad)
        xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                        plot = FALSE)
        barplot(t(newd), xlab = "Bmsy/K", ylab = "Frequency", 
                col = c("black", "NA"), cex.main = 0.7, main = paste("Accepted=Black", 
                                                                     "  ", "Rejected=White"))
        axis(1, at = xpos, labels = FALSE)
        breaks <- hist(storep[, "M"], plot = FALSE, 
                       nclass = grunits[14, 5])$breaks
        good <- table((cut(storep[storep[, 1] == 1, 
                                  "M"], breaks)))
        dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                             1)] + breaks[2:c(length(breaks))])/2
        bad <- table((cut(storep[storep[, 1] == 0, "M"], 
                          breaks)))
        dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                            1)] + breaks[2:c(length(breaks))])/2
        newd <- cbind(good, bad)
        xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                        plot = FALSE)
        barplot(t(newd), xlab = "M", ylab = "Frequency", 
                col = c("black", "NA"), cex.main = 0.7, main = paste("Accepted=Black", 
                                                                     "  ", "Rejected=White"))
        axis(1, at = xpos, labels = FALSE)
        breaks <- hist(storep[, "BtK"], plot = FALSE, 
                       nclass = grunits[14, 5])$breaks
        good <- table((cut(storep[storep[, 1] == 1, 
                                  "BtK"], breaks)))
        dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                             1)] + breaks[2:c(length(breaks))])/2
        bad <- table((cut(storep[storep[, 1] == 0, "BtK"], 
                          breaks)))
        dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                            1)] + breaks[2:c(length(breaks))])/2
        newd <- cbind(good, bad)
        xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                        plot = FALSE)
        barplot(t(newd), xlab = "Bt/K", ylab = "Frequency", 
                col = c("black", "NA"), cex.main = 0.7, main = paste("Accepted=Black", 
                                                                     "  ", "Rejected=White"))
        axis(1, at = xpos, labels = FALSE)
        breaks <- hist(storep[, "FmsyM"], plot = FALSE, 
                       nclass = grunits[14, 5])$breaks
        good <- table((cut(storep[storep[, 1] == 1, 
                                  "FmsyM"], breaks)))
        dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                             1)] + breaks[2:c(length(breaks))])/2
        bad <- table((cut(storep[storep[, 1] == 0, "FmsyM"], 
                          breaks)))
        dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                            1)] + breaks[2:c(length(breaks))])/2
        newd <- cbind(good, bad)
        xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                        plot = FALSE)
        barplot(t(newd), xlab = "Fmsy/M", ylab = "Frequency", 
                col = c("black", "NA"), cex.main = 0.7, main = paste("Accepted=Black", 
                                                                     "  ", "Rejected=White"))
        axis(1, at = xpos, labels = FALSE)
        breaks <- hist(storep[, "K"], plot = FALSE, 
                       nclass = grunits[14, 5])$breaks
        good <- table((cut(storep[storep[, 1] == 1, 
                                  "K"], breaks)))
        dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                             1)] + breaks[2:c(length(breaks))])/2
        bad <- table((cut(storep[storep[, 1] == 0, "K"], 
                          breaks)))
        dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                            1)] + breaks[2:c(length(breaks))])/2
        newd <- cbind(good, bad)
        xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                        plot = FALSE)
        barplot(t(newd), xlab = "K", ylab = "Frequency", 
                col = c("black", "NA"), cex.main = 0.7, main = paste("Accepted=Black", 
                                                                     "  ", "Rejected=White"))
        axis(1, at = xpos, labels = FALSE)
        breaks <- hist(storep[, "MSY"], plot = FALSE, 
                       nclass = grunits[14, 5])$breaks
        good <- table((cut(storep[storep[, 1] == 1, 
                                  "MSY"], breaks)))
        dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                             1)] + breaks[2:c(length(breaks))])/2
        bad <- table((cut(storep[storep[, 1] == 0, "MSY"], 
                          breaks)))
        dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                            1)] + breaks[2:c(length(breaks))])/2
        newd <- cbind(good, bad)
        xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                        plot = FALSE)
        barplot(t(newd), xlab = "MSY", ylab = "Frequency", 
                col = c("black", "NA"), cex.main = 0.7, main = paste("Accepted=Black", 
                                                                     "  ", "Rejected=White"))
        axis(1, at = xpos, labels = FALSE)
        breaks <- hist(storep[, "Bmsy"], plot = FALSE, 
                       nclass = grunits[14, 5])$breaks
        good <- table((cut(storep[storep[, 1] == 1, 
                                  "Bmsy"], breaks)))
        dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                             1)] + breaks[2:c(length(breaks))])/2
        bad <- table((cut(storep[storep[, 1] == 0, "Bmsy"], 
                          breaks)))
        dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                            1)] + breaks[2:c(length(breaks))])/2
        newd <- cbind(good, bad)
        xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                        plot = FALSE)
        barplot(t(newd), xlab = "Bmsy", ylab = "Frequency", 
                col = c("black", "NA"), cex.main = 0.7, main = paste("Accepted=Black", 
                                                                     "  ", "Rejected=White"))
        axis(1, at = xpos, labels = FALSE)
        breaks <- hist(storep[, "Umsy"], plot = FALSE, 
                       nclass = grunits[14, 5])$breaks
        good <- table((cut(storep[storep[, 1] == 1, 
                                  "Umsy"], breaks)))
        dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                             1)] + breaks[2:c(length(breaks))])/2
        bad <- table((cut(storep[storep[, 1] == 0, "Umsy"], 
                          breaks)))
        dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                            1)] + breaks[2:c(length(breaks))])/2
        newd <- cbind(good, bad)
        xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                        plot = FALSE)
        barplot(t(newd), xlab = "Umsy", ylab = "Frequency", 
                col = c("black", "NA"), cex.main = 0.7, main = paste("Accepted=Black", 
                                                                     "  ", "Rejected=White"))
        axis(1, at = xpos, labels = FALSE)
        if (grout == 2) {
          word.tif("BmsykAR")
          breaks <- hist(storep[, "BmsyK"], plot = FALSE, 
                         nclass = grunits[14, 5])$breaks
          good <- table((cut(storep[storep[, 1] == 1, 
                                    "BmsyK"], breaks)))
          dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                               1)] + breaks[2:c(length(breaks))])/2
          bad <- table((cut(storep[storep[, 1] == 0, 
                                   "BmsyK"], breaks)))
          dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                              1)] + breaks[2:c(length(breaks))])/2
          newd <- cbind(good, bad)
          xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                          plot = FALSE)
          barplot(t(newd), xlab = "Bmsy/K", ylab = "Frequency", 
                  col = c("black", "NA"), cex.main = 0.7, 
                  main = paste("Accepted=Black", "  ", "Rejected=White"))
          axis(1, at = xpos, labels = FALSE)
          dev.off()
          word.tif("MAR")
          breaks <- hist(storep[, "M"], plot = FALSE, 
                         nclass = grunits[14, 5])$breaks
          good <- table((cut(storep[storep[, 1] == 1, 
                                    "M"], breaks)))
          dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                               1)] + breaks[2:c(length(breaks))])/2
          bad <- table((cut(storep[storep[, 1] == 0, 
                                   "M"], breaks)))
          dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                              1)] + breaks[2:c(length(breaks))])/2
          newd <- cbind(good, bad)
          xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                          plot = FALSE)
          barplot(t(newd), xlab = "M", ylab = "Frequency", 
                  col = c("black", "NA"), cex.main = 0.7, 
                  main = paste("Accepted=Black", "  ", "Rejected=White"))
          axis(1, at = xpos, labels = FALSE)
          dev.off()
          word.tif("BtKAR")
          breaks <- hist(storep[, "BtK"], plot = FALSE, 
                         nclass = grunits[14, 5])$breaks
          good <- table((cut(storep[storep[, 1] == 1, 
                                    "BtK"], breaks)))
          dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                               1)] + breaks[2:c(length(breaks))])/2
          bad <- table((cut(storep[storep[, 1] == 0, 
                                   "BtK"], breaks)))
          dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                              1)] + breaks[2:c(length(breaks))])/2
          newd <- cbind(good, bad)
          xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                          plot = FALSE)
          barplot(t(newd), xlab = "Bt/K", ylab = "Frequency", 
                  col = c("black", "NA"), cex.main = 0.7, 
                  main = paste("Accepted=Black", "  ", "Rejected=White"))
          axis(1, at = xpos, labels = FALSE)
          dev.off()
          word.tif("B1KAR")
          breaks <- hist(storep[, "B1K"], plot = FALSE, 
                         nclass = grunits[14, 5])$breaks
          good <- table((cut(storep[storep[, 1] == 1, 
                                    "B1K"], breaks)))
          dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                               1)] + breaks[2:c(length(breaks))])/2
          bad <- table((cut(storep[storep[, 1] == 0, 
                                   "B1K"], breaks)))
          dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                              1)] + breaks[2:c(length(breaks))])/2
          newd <- cbind(good, bad)
          xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                          plot = FALSE)
          barplot(t(newd), xlab = "B1/K", ylab = "Frequency", 
                  col = c("black", "NA"), cex.main = 0.7, 
                  main = paste("Accepted=Black", "  ", "Rejected=White"))
          axis(1, at = xpos, labels = FALSE)
          dev.off()
          word.tif("FmsyMAR")
          breaks <- hist(storep[, "FmsyM"], plot = FALSE, 
                         nclass = grunits[14, 5])$breaks
          good <- table((cut(storep[storep[, 1] == 1, 
                                    "FmsyM"], breaks)))
          dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                               1)] + breaks[2:c(length(breaks))])/2
          bad <- table((cut(storep[storep[, 1] == 0, 
                                   "FmsyM"], breaks)))
          dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                              1)] + breaks[2:c(length(breaks))])/2
          newd <- cbind(good, bad)
          xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                          plot = FALSE)
          barplot(t(newd), xlab = "Fmsy/M", ylab = "Frequency", 
                  col = c("black", "NA"), cex.main = 0.7, 
                  main = paste("Accepted=Black", "  ", "Rejected=White"))
          axis(1, at = xpos, labels = FALSE)
          dev.off()
          word.tif("KAR")
          breaks <- hist(storep[, "K"], plot = FALSE, 
                         nclass = grunits[14, 5])$breaks
          good <- table((cut(storep[storep[, 1] == 1, 
                                    "K"], breaks)))
          dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                               1)] + breaks[2:c(length(breaks))])/2
          bad <- table((cut(storep[storep[, 1] == 0, 
                                   "K"], breaks)))
          dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                              1)] + breaks[2:c(length(breaks))])/2
          newd <- cbind(good, bad)
          xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                          plot = FALSE)
          barplot(t(newd), xlab = "K", ylab = "Frequency", 
                  col = c("black", "NA"), cex.main = 0.7, 
                  main = paste("Accepted=Black", "  ", "Rejected=White"))
          axis(1, at = xpos, labels = FALSE)
          dev.off()
          word.tif("MSYAR")
          breaks <- hist(storep[, "MSY"], plot = FALSE, 
                         nclass = grunits[14, 5])$breaks
          good <- table((cut(storep[storep[, 1] == 1, 
                                    "MSY"], breaks)))
          dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                               1)] + breaks[2:c(length(breaks))])/2
          bad <- table((cut(storep[storep[, 1] == 0, 
                                   "MSY"], breaks)))
          dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                              1)] + breaks[2:c(length(breaks))])/2
          newd <- cbind(good, bad)
          xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                          plot = FALSE)
          barplot(t(newd), xlab = "MSY", ylab = "Frequency", 
                  col = c("black", "NA"), cex.main = 0.7, 
                  main = paste("Accepted=Black", "  ", "Rejected=White"))
          axis(1, at = xpos, labels = FALSE)
          dev.off()
          word.tif("BmsyAR")
          breaks <- hist(storep[, "Bmsy"], plot = FALSE, 
                         nclass = grunits[14, 5])$breaks
          good <- table((cut(storep[storep[, 1] == 1, 
                                    "Bmsy"], breaks)))
          dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                               1)] + breaks[2:c(length(breaks))])/2
          bad <- table((cut(storep[storep[, 1] == 0, 
                                   "Bmsy"], breaks)))
          dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                              1)] + breaks[2:c(length(breaks))])/2
          newd <- cbind(good, bad)
          xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                          plot = FALSE)
          barplot(t(newd), xlab = "Bmsy", ylab = "Frequency", 
                  col = c("black", "NA"), cex.main = 0.7, 
                  main = paste("Accepted=Black", "  ", "Rejected=White"))
          axis(1, at = xpos, labels = FALSE)
          dev.off()
          word.tif("UmsyAR")
          breaks <- hist(storep[, "Umsy"], plot = FALSE, 
                         nclass = grunits[14, 5])$breaks
          good <- table((cut(storep[storep[, 1] == 1, 
                                    "Umsy"], breaks)))
          dimnames(good)[[1]] <- (breaks[1:c(length(breaks) - 
                                               1)] + breaks[2:c(length(breaks))])/2
          bad <- table((cut(storep[storep[, 1] == 0, 
                                   "Umsy"], breaks)))
          dimnames(bad)[[1]] <- (breaks[1:c(length(breaks) - 
                                              1)] + breaks[2:c(length(breaks))])/2
          newd <- cbind(good, bad)
          xpos <- barplot(t(newd), ylim = c(0, max(newd)), 
                          plot = FALSE)
          barplot(t(newd), xlab = "Umsy", ylab = "Frequency", 
                  col = c("black", "NA"), cex.main = 0.7, 
                  main = paste("Accepted=Black", "  ", "Rejected=White"))
          axis(1, at = xpos, labels = FALSE)
          dev.off()
        }
      }
    }
    outs <- data.frame(Distr = NA, Min = NA, Max = NA, Mean = NA, 
                       sd = NA)
    outs[1, ] <- cbind(fmsymdef$dist, fmsymdef$low, fmsymdef$up, 
                       fmsymdef$mean, fmsymdef$sd)
    outs[2, ] <- cbind(btkdef$dist, btkdef$low, btkdef$up, 
                       btkdef$mean, btkdef$sd)
    outs[3, ] <- cbind(bmsykdef$dist, bmsykdef$low, bmsykdef$up, 
                       bmsykdef$mean, bmsykdef$sd)
    outs[4, ] <- cbind(Mdef$dist, Mdef$low, Mdef$up, Mdef$mean, 
                       Mdef$sd)
    outs[5, ] <- cbind(NA, btkdef$refyr, NA, NA, NA)
    outs[6, ] <- cbind(b1kdef$dist, b1kdef$low, b1kdef$up, 
                       b1kdef$mean, b1kdef$sd)
    colnames(outs) <- c("Distr", "Lower", "Upper", "Mean ", 
                        "SD")
    rownames(outs) <- c("Fmsy/M", "Br/K", "Bmsy/K", "M", 
                        "refyr", "B1/K")
    outs1 <- data.frame(Mean = NA, Median = NA, per2_5 = NA, 
                        per97_5 = NA, min = NA, max = NA)
    outs1[1, ] <- cbind(mFmsyM, FmsyM95[[2]], FmsyM95[[1]], 
                        FmsyM95[[3]], min(datar$FmsyM), max(datar$FmsyM))
    outs1[2, ] <- cbind(mbtk, btk95[[2]], btk95[[1]], btk95[[3]], 
                        min(datar$BtK), max(datar$BtK))
    outs1[3, ] <- cbind(mBmsyK, BmsyK95[[2]], BmsyK95[[1]], 
                        BmsyK95[[3]], min(datar$BmsyK), max(datar$BmsyK))
    outs1[4, ] <- cbind(mM, M95[[2]], M95[[1]], M95[[3]], 
                        min(datar$M), max(datar$M))
    outs1[5, ] <- cbind(mb1k, b1k95[[2]], b1k95[[1]], b1k95[[3]], 
                        min(datar$B1K), max(datar$B1K))
    colnames(outs1) <- c("Mean (ll=1)", "Median (ll=1)", 
                         "2.5% (ll=1)", "97.5% (ll=1)", "min (ll=1)", "max (ll=1)")
    rownames(outs1) <- c("Fmsy/M", "Bt/K", "Bmsy/K", "M", 
                         "B1/K")
    outs2 <- data.frame(Mean = NA, Median = NA, per2_5 = NA, 
                        per97_5 = NA, min = NA, max = NA)
    outs2[1, ] <- cbind(mMSY, MSY95[[2]], MSY95[[1]], MSY95[[3]], 
                        min(datar$MSY), max(datar$MSY))
    outs2[2, ] <- cbind(mBMSY, BMSY95[[2]], BMSY95[[1]], 
                        BMSY95[[3]], min(datar$Bmsy), max(datar$Bmsy))
    outs2[3, ] <- cbind(mFMSY, FMSY95[[2]], FMSY95[[1]], 
                        FMSY95[[3]], min(datar$Fmsy), max(datar$Fmsy))
    outs2[4, ] <- cbind(mUMSY, UMSY95[[2]], UMSY95[[1]], 
                        UMSY95[[3]], min(datar$Umsy), max(datar$Umsy))
    outs2[5, ] <- cbind(mOFL, OFL95[[2]], OFL95[[1]], OFL95[[3]], 
                        min(datar$OFL), max(datar$OFL))
    outs2[6, ] <- cbind(mBrefyr, Brefyr95[[2]], Brefyr95[[1]], 
                        Brefyr95[[3]], min(datar$Brefyr), max(datar$Brefyr))
    outs2[7, ] <- cbind(mk, k95[[2]], k95[[1]], k95[[3]], 
                        min(datar$K), max(datar$K))
    colnames(outs2) <- c("Mean (ll=1)", "Median (ll=1)", 
                         "2.5% (ll=1)", "97.5% (ll=1)", "min (ll=1)", "max (ll=1)")
    rownames(outs2) <- c("MSY", "Bmsy", "Fmsy", "Umsy", 
                         "OFL", "Brefyr", "K")
    ans <- list(outs, outs1, outs2, storep, agemat, c(max(year) + 1), "dbsra",
                timeseries_B,timeseries_BK,timeseries_BBmsy,timeseries_FF,timeseries_FFmsy)
    names(ans) <- c("Initial", "Parameters", "Estimates","Values", "agemat", "end1yr", "type",
                    "Biom.traj","Depletion.traj","B.Bmsy","F.series","F.Fmsy")
  }  
  if (length(datar[, 1]) == 0) {
    outs <- data.frame(Distr = NA, Min = NA, Max = NA, Mean = NA, 
                       sd = NA)
    outs[1, ] <- cbind(fmsymdef$dist, fmsymdef$low, fmsymdef$up, 
                       fmsymdef$mean, fmsymdef$sd)
    outs[2, ] <- cbind(btkdef$dist, btkdef$low, btkdef$up, 
                       btkdef$mean, btkdef$sd)
    outs[3, ] <- cbind(bmsykdef$dist, bmsykdef$low, bmsykdef$up, 
                       bmsykdef$mean, bmsykdef$sd)
    outs[4, ] <- cbind(Mdef$dist, Mdef$low, Mdef$up, Mdef$mean, 
                       Mdef$sd)
    outs[5, ] <- cbind(NA, btkdef$refyr, NA, NA, NA)
    outs[6, ] <- cbind(b1kdef$dist, b1kdef$low, b1kdef$up, 
                       b1kdef$mean, b1kdef$sd)
    colnames(outs) <- c("Distr", "Lower", "Upper", "Mean ", 
                        "SD")
    rownames(outs) <- c("Fmsy/M", "Bt/K", "Bmsy/K", "M", 
                        "refyr", "B1/K")
    outs1 <- data.frame(Mean = NA, Median = NA, per2_5 = NA, 
                        per97_5 = NA, min = NA, max = NA)
    colnames(outs1) <- c("Mean (ll=1)", "Median (ll=1)", 
                         "2.5% (ll=1)", "97.5% (ll=1)")
    outs2 <- data.frame(Mean = NA, Median = NA, per2_5 = NA, 
                        per97_5 = NA, min = NA, max = NA)
    colnames(outs2) <- c("Mean (ll=1)", "Median (ll=1)", 
                         "2.5% (ll=1)", "97.5% (ll=1)", "min (ll=1)", "max (ll=1)")
    ans <- list(outs, outs1, outs2, storep, agemat, c(max(year) + 
                                                        1), "dbsra")
    names(ans) <- c("Initial", "Parameters", "Estimates", 
                    "Values", "agemat", "end1yr", "type")
    warning("None of the runs had a likelihood equal to 1")
  }
  return(ans)
}
dlproj.tweaked=function (dlobj = NULL, projyears = NULL, projtype = 1, projcatch = NULL, 
                         grout = 1, grargs = list(lwd = 1, unit = "MT", mains = " ", 
                                                  cex.main = 1, cex.axis = 1, cex.lab = 1), grtif = list(zoom = 4, 
                                                                                                         width = 11, height = 13, pointsize = 10)) 
{
  if (is.null(dlobj)) 
    stop("dlobj is missing")
  if (!any(names(dlobj) == "type")) 
    stop("dlobj not a catchmsy or dbsra output object file")
  if (!dlobj$type %in% c("dbsra", "catchmsy")) 
    stop("dlobj not a catchmsy or dbsra output object file")
  if (dlobj$type %in% c("dbsra", "catchmsy")) {
    if (projtype == 2 & is.null(projcatch)) 
      stop("For projtype=2, enter value for projcatch")
    if (is.null(projyears)) 
      stop("projyears is missing")
  }
  grdef = list(lwd = 1, unit = "MT", cex.axis = 1, cex.lab = 1, 
               mains = " ", cex.main = 1)
  if (any(names(grargs) == "unit")) {
    if (any(grdef$unit != grargs$unit)) 
      grdef$unit <- grargs$unit
  }
  if (any(names(grargs) == "cex.axis")) {
    if (any(grdef$cex.axis != grargs$cex.axis)) 
      grdef$cex.axis <- grargs$cex.axis
  }
  if (any(names(grargs) == "cex.lab")) {
    if (any(grdef$cex.lab != grargs$cex.lab)) 
      grdef$cex.lab <- grargs$cex.lab
  }
  if (any(names(grargs) == "mains")) {
    if (any(grdef$mains != grargs$mains)) 
      grdef$mains <- grargs$mains
  }
  if (any(names(grargs) == "cex.main")) {
    if (any(grdef$cex.main != grargs$cex.main)) 
      grdef$cex.main <- grargs$cex.main
  }
  if (any(names(grargs) == "lwd")) {
    if (any(grdef$lwd != grargs$lwd)) 
      grdef$lwd <- grargs$lwd
  }
  tifdef = list(zoom = 4, width = 11, height = 13, pointsize = 10)
  if (any(names(grtif) == "zoom")) {
    if (tifdef$zoom != grtif$zoom) 
      tifdef$zoom <- grtif$zoom
  }
  if (any(names(grtif) == "width")) {
    if (tifdef$width != grtif$width) 
      tifdef$width <- grtif$width
  }
  if (any(names(grtif) == "height")) {
    if (tifdef$height != grtif$height) 
      tifdef$height <- grtif$height
  }
  if (any(names(grtif) == "pointsize")) {
    if (tifdef$pointsize != grtif$pointsize) 
      tifdef$pointsize <- grtif$pointsize
  }
  if (dlobj$type == "dbsra") {
    agemat <- dlobj$agemat
    if (projtype < 2) {
      if (projtype == 0) 
        catch1 <- dlobj$Estimates[1, 2]
      if (projtype == 1) 
        catch1 <- dlobj$Estimates[1, 1]
      catch1 <- rep(catch1, projyears)
    }
    if (projtype == 2) {
      if (is.null(projcatch)) 
        stop("projcatch must be specified")
      if (length(projcatch) > 1 && length(projcatch) < 
          projyears) 
        stop("length of projcatch does not match length of projyears")
      if (length(projcatch) > 1) 
        catch1 <- projcatch
      if (length(projcatch) == 1) 
        catch1 <- rep(projcatch, projyears)
    }
    parms <- dlobj$Values[dlobj$Values[, 1] == 1, c(2, 4, 
                                                    5, 6, 7, 8, 9, 14, 15)]
    names(parms) <- c("FmsyM", "BmsyK", "M", 
                      "K", "Fmsy", "Umsy", "MSY", 
                      "n", "g")
    P <- NULL
    bigB <- dlobj$Biom.traj
    bigB <- t(bigB)
    pB <- as.data.frame(bigB[c(nrow(bigB) - agemat):nrow(bigB),])
    star <- nrow(pB)
    pB[c(nrow(pB) + 1):c(c(nrow(pB) + 1) + projyears - 1),] <- NA
    for (nn in 1:nrow(parms)) {
      P <- 0
      MSY <- parms$MSY[nn]
      for (t in star:c(nrow(pB) - 1)) {
        if (parms$BmsyK[nn] >= 0.5) 
          P <- parms$g[nn] * MSY * (pB[t - agemat, nn]/parms$K[nn]) - 
            parms$g[nn] * MSY * (pB[t - agemat, nn]/parms$K[nn])^parms$n[nn]
        if (parms$BmsyK[nn] > 0.3 & parms$BmsyK[nn] < 
            0.5) {
          bjoin <- (0.75 * parms$BmsyK[nn] - 0.075) * 
            parms$K[nn]
          if (pB[t - agemat, nn] < bjoin) {
            PJ <- parms$g[nn] * MSY * (bjoin/parms$K[nn]) - 
              parms$g[nn] * MSY * (bjoin/parms$K[nn])^parms$n[nn]
            cc <- (1 - parms$n[nn]) * MSY * parms$g[nn] * 
              (bjoin^(parms$n[nn] - 2)) * parms$K[nn]^-parms$n[nn]
            P <- pB[t - agemat, nn] * (PJ/bjoin + cc * 
                                         (pB[t - agemat, nn] - bjoin))
          }
          if (pB[t - agemat, nn] >= bjoin) 
            P <- parms$g[nn] * MSY * (pB[t - agemat, 
                                         nn]/parms$K[nn]) - parms$g[nn] * MSY * 
              (pB[t - agemat, nn]/parms$K[nn])^parms$n[nn]
        }
        if (parms$BmsyK[nn] <= 0.3) {
          bjoin <- (0.5 * parms$BmsyK[nn]) * parms$K[nn]
          if (pB[t - agemat, nn] < bjoin) {
            PJ <- parms$g[nn] * MSY * (bjoin/parms$K[nn]) - 
              parms$g[nn] * MSY * (bjoin/parms$K[nn])^parms$n[nn]
            cc <- (1 - parms$n[nn]) * MSY * parms$g[nn] * 
              (bjoin^(parms$n[nn] - 2)) * parms$K[nn]^-parms$n[nn]
            P <- pB[t - agemat, nn] * (PJ/bjoin + cc * 
                                         (pB[t - agemat, nn] - bjoin))
          }
          if (pB[t - agemat, nn] >= bjoin) 
            P <- parms$g[nn] * MSY * (pB[t - agemat, 
                                         nn]/parms$K[nn]) - parms$g[nn] * MSY * 
              (pB[t - agemat, nn]/parms$K[nn])^parms$n[nn]
        }
        P <- ifelse(P < 0, 0, P)
        pB[t + 1, nn] <- max(0, pB[t, nn] + P - catch1[t - 
                                                         (star - 1)])
      }
    }
    dd <- pB[star:c(length(pB[, 1])), ]
    rownames(dd) <- c(dlobj$end1yr:c(dlobj$end1yr + projyears))
    ans <- list(dlobj$type, t(dd))
    names(ans) <- c("type", "ProjBio")
    return(ans)
  }
  
}
apply.DBSRA_tweeked=function(year,catch,catchCV,catargs,agemat,k,b1k,btk,fmsym,bmsyk,M,graph,
                             nsims=Niters,grout,WD,outfile,Projections=NULL)
{
  setwd(WD)  #dbsra automatically exports the biomass trajectories
  #store inputs
  input=list(year=year,
             catch=catch,
             catchCV=catchCV,
             catargs=catargs,
             agemat=agemat,
             k=k,
             b1k=b1k,
             btk=btk,
             fmsym=fmsym,
             bmsyk=bmsyk,
             M=M,
             graph=graph,
             nsims=nsims,
             grout=grout)
  
  #run model
  #pdf(paste(outfile,".pdf",sep=''))
  output <- dbsra_tweeked(year=year,
                          catch=catch,
                          catchCV=catchCV,
                          catargs=catargs,
                          agemat=agemat,
                          k=k,
                          b1k=b1k,
                          btk=btk,
                          fmsym=fmsym,
                          bmsyk=bmsyk,
                          M=M,
                          graphs=graph,
                          nsims=nsims,
                          grout=grout)
  
  #projections
  if(!is.null(Projections))
  {
    Future.projections<-dlproj.tweaked(dlobj = output,
                                       projyears = length(Projections$finyear),
                                       projtype = 2,
                                       projcatch = Projections$Tonnes, 
                                       grout = 0)
    output$Biom.traj=cbind(output$Biom.traj,Future.projections$ProjBio)
    output$Depletion.traj=cbind(output$Depletion.traj,Future.projections$ProjBio/output$Values%>%filter(ll==1)%>%pull(K))
    adumi=output$B.Bmsy[,1:ncol(Future.projections$ProjBio)]
    adumi[,]=NA
    output$B.Bmsy=cbind(output$B.Bmsy,adumi)
    output$F.series=cbind(output$F.series,adumi)
    output$F.Fmsy=cbind(output$F.Fmsy,adumi)
    year=c(year,Projections$finyear)  
  }
  
  #extract years
  output$Years=year
  return(list(input=input,output=output))
}
apply.DBSRA=function(year,catch,catchCV,catargs,agemat,k,b1k,btk,fmsym,bmsyk,M,graph,
                     nsims=Niters,grout,WD,outfile)
{
  setwd(WD)  #dbsra automatically exports the biomass trajectories
  #store inputs
  input=list(year=year,
             catch=catch,
             catchCV=catchCV,
             catargs=catargs,
             agemat=agemat,
             k=k,
             b1k=b1k,
             btk=btk,
             fmsym=fmsym,
             bmsyk=bmsyk,
             M=M,
             graph=graph,
             nsims=nsims,
             grout=grout)
  
  #run model
  #pdf(paste(outfile,".pdf",sep=''))
  output <- dbsra(year=year,
                  catch=catch,
                  catchCV=catchCV,
                  catargs=catargs,
                  agemat=agemat,
                  k=k,
                  b1k=b1k,
                  btk=btk,
                  fmsym=fmsym,
                  bmsyk=bmsyk,
                  M=M,
                  graphs=graph,
                  nsims=nsims,
                  grout=grout)
  #dev.off()
  
  #extract biomass trajectories
  Biom.traj=read.csv(paste(WD,"Biotraj-dbsra.csv",sep='/'),header=FALSE)
  
  output$Biom.traj=Biom.traj%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)
  output$Depletion.traj=Biom.traj%>%
    mutate_at(vars(- V1), ~ . / output$Values$K)%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)
  output$B.Bmsy=Biom.traj%>%
    mutate_at(vars(- V1), ~ . / output$Values$Bmsy)%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)
  
  #extract F trajectories
  U.series=Biom.traj[,-ncol(Biom.traj)]
  U.series[,-1]=mapply('/',catch,U.series[,-1])
  
  F.series=U.series%>%
    mutate_at(vars(- V1), function(x) -log(1-x))
  
  output$F.series=F.series%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)
  #output$F.series=F.from.U((catch/apply(output$Biom.traj,2,median)[1:length(catch)]))
  
  output$F.Fmsy=F.series%>%
    mutate_at(vars(- V1), ~ . / output$Values$Fmsy)%>%
    filter(V1==1)%>%  #select only possible runs
    dplyr::select(-V1)
  
  #extract years
  output$Years=year
  
  
  return(list(input=input,output=output))
}

Res.fn=function(r,Def)
{
  if(Def=="Martell")
  {
    if(r<0.1) out="verylow" else if
    (r>=0.1 & r <0.6) out="low" else if
    (r>=0.6 & r <1.5 ) out="medium" else if
    (r >=1.5) out="high"
  }
  
  if(Def=="Haddon")
  {
    if(r<0.1) out="verylow" else if
    (r <= mean(c(0.1,0.6))) out="low" else if
    (r>mean(c(0.1,0.6)) & r <=mean(c(0.3,0.8)) ) out="medium" else if
    (r >mean(c(0.3,0.8))) out="high"
  }
  return(out)
}
apply.CMSY=function(year,catch,r.range,k.range,Bo.low,Bo.hi,Int.yr=NA,Bint.low=NA,Bint.hi=NA,
                    Bf.low,Bf.hi,outfile,CMSY=CMSY.method,nsims=Niters,Proc.error,RES,Projections=NULL)
{
  input=list(r.range=r.range,
             stb.low=Bo.low,
             stb.hi=Bo.hi,
             int.yr=Int.yr,
             intb.low = Bint.low,
             intb.hi = Bint.hi,
             endb.low=Bf.low,
             endb.hi=Bf.hi,
             RES=RES,
             Proc.error=Proc.error,
             k.range=k.range)
  
  if(CMSY=="Froese")
  {
    output <- cmsy2(year=year,
                    catch=catch,
                    r.low=r.range[1],
                    r.hi=r.range[2],
                    stb.low=Bo.low,
                    stb.hi=Bo.hi,
                    int.yr=Int.yr,
                    intb.low = Bint.low,
                    intb.hi = Bint.hi,
                    endb.low=Bf.low,
                    endb.hi=Bf.hi)
    
    # Extract reference points and time series from output
    output$ref_pts <- output[["ref_pts"]]
    output$ref_ts <- output[["ref_ts"]]
    
    #Appendix
    appendix=data.frame(r=output$r_viable,k=1000*output$k_viable)
    p1=appendix%>%
      ggplot(aes(r))+
      geom_histogram()
    p2=appendix%>%
      ggplot(aes(k))+
      geom_histogram()
    p3=appendix%>%
      ggplot(aes(k,r))+
      geom_point() +
      geom_density_2d_filled(alpha = 0.5)+
      geom_density_2d(size = 0.25, colour = "black")
    ggarrange(plotlist = list(p1,p2,p3),nrow=1,widths = c(1,1,2))
    ggsave(paste(outfile,'.tiff',sep=''),width = 14,height = 8, dpi = 300, compression = "lzw")  
    
  }
  if(CMSY=="Haddon")
  {
    glb=list(resilience=RES,spsname="")
    indat=data.frame(year=year,catch=catch)
    output <- run_cMSY(indat=indat,
                       glob=glb,
                       n = nsims,
                       incB = 0.025,
                       sigpR = Proc.error,
                       multK = 1,
                       finaldepl = c(Bf.low,Bf.hi),
                       start_k = k.range,
                       start_r = r.range,
                       initdepl = c(Bo.low,Bo.hi),
                       maximumH = 1,
                       Hyear = NA)
    #Appendix
    ciMSY=summarycMSY(output,indat,final=TRUE)
    pdf(paste(outfile,'.pdf',sep=''))
    plotcMSY6(ciMSY,indat[,"catch"])
    dev.off()
    
    #get B.bmsy and F.fmsy
    d1=pulloutStats(output$R1,probabs=c(0.5))
    d1=d1$traject
    d1=data.frame(d1)
    Kei=d1$K
    Ar=d1$r
    Fmsy=Ar/2
    Bmsy=Kei/2
    d1=d1%>%
      dplyr::select(-c(r,K,bd))
    U.series=mapply('/',catch,d1[,1:length(catch)])
    F.series=apply(U.series,2,function(x) -log(1-x))
    
    output$B.traj=d1
    output$Depletion.traj=d1/Kei
    output$Bmsy=Bmsy
    output$B.Bmsy=d1/Bmsy
    output$F.Fmsy=F.series/Fmsy
    output$F.traj=F.series
    output$Years=year
    output$acceptance.rate=round(100*length(ciMSY$r)/nsims,2)  
  }
  
  #projections
  if(!is.null(Projections))
  {
    Future.projections<-doproject(intraj=gettraject(output$R1,
                                                    projn=length(Projections$finyear)),
                                  constC=mean(Projections$Tonnes),
                                  projn=length(Projections$finyear),
                                  sigpR=Proc.error)
    Future.projections=Future.projections[,match(Projections$finyear[-1],colnames(Future.projections))]
    colnames(Future.projections)=paste("X",colnames(Future.projections),sep='')
    output$B.traj=cbind(output$B.traj,Future.projections)
    
    output$Depletion.traj=cbind(output$Depletion.traj,Future.projections/Kei)
    adumi=output$B.Bmsy[,1:ncol(Future.projections)]
    adumi[,]=NA
    colnames(adumi)=paste("X",Projections$finyear[-1],sep='')
    output$B.Bmsy=cbind(output$B.Bmsy,adumi)
    adumi=output$F.traj[,1:(ncol(Future.projections)+1)]
    adumi[,]=NA
    output$F.traj=cbind(output$F.traj,adumi)
    output$F.Fmsy=cbind(output$F.Fmsy,adumi)
    output$Years=c(output$Years,Projections$finyear)
  }
  
  return(list(input=input,output=output))
}
apply.OCOM=function(year,catch,M,outfile)
{
  input= list(m=M)
  output <- ocom(year=year,catch=catch,m=M)
  
  # Extract reference points and time series from output
  output$ref_pts <- output[["ref_pts"]]
  output$ref_ts <- output[["ref_ts"]]
  
  #Appendix
  pdf(paste(outfile,'.pdf',sep=''))
  plot_dlm(output)
  dev.off()
  
  appendix=data.frame(r=output$krms_draws$r,k=output$krms_draws$k)
  p1=appendix%>%
    ggplot(aes(r))+
    geom_histogram()
  p2=appendix%>%
    ggplot(aes(k))+
    geom_histogram()
  p3=appendix%>%
    ggplot(aes(k,r))+
    geom_point() +
    geom_density_2d_filled(alpha = 0.5)+
    geom_density_2d(size = 0.25, colour = "black")
  ggarrange(plotlist = list(p1,p2,p3),nrow=1,widths = c(1,1,2))
  ggsave(paste(outfile,'.tiff',sep=''),width = 14,height = 8, dpi = 300, compression = "lzw")  
  
  
  return(list(input=input,output=output))
}
apply.JABBA=function(Ktch,CPUE,CPUE.SE,auxil=NULL,auxil.se=NULL,auxil.type=NULL,
                     MDL,Ktch.CV,ASS,Rdist,Rprior,Kdist,Kprior,
                     PsiDist,Psiprior,Bprior,BMSYK,Shape.CV,output.dir,outfile,Sims,
                     Proc.error.JABBA,Obs.Error=NULL,
                     thinning = 5,nchains = 2,burn.in=5000,ktch.error="random",Projections=NULL)
{
  # Compile JABBA JAGS model
  if(is.null(CPUE))
  {
    jbinput = build_jabba(catch=Ktch,
                          model.type = "Schaefer",
                          catch.cv=Ktch.CV,
                          catch.error=ktch.error,
                          assessment=ASS,
                          scenario =  "CatchOnly",
                          r.dist = Rdist,
                          r.prior = Rprior,
                          K.dist= Kdist,
                          K.prior=Kprior,
                          psi.dist=PsiDist,
                          psi.prior=Psiprior,
                          b.prior=Bprior,
                          BmsyK = BMSYK,
                          shape.CV=Shape.CV,
                          sigma.est = FALSE,
                          sigma.proc=Proc.error.JABBA)
  }else
  {
    jbinput = build_jabba(catch=Ktch,
                          cpue=CPUE,
                          se=CPUE.SE,
                          auxiliary=auxil,
                          auxiliary.se=auxil.se,
                          auxiliary.type=auxil.type,
                          model.type = MDL,
                          catch.cv=Ktch.CV,
                          catch.error=ktch.error,
                          assessment=ASS,
                          scenario =  "CPUE",
                          r.dist = Rdist,
                          r.prior = Rprior,
                          K.dist= Kdist,
                          K.prior=Kprior,
                          psi.dist=PsiDist,
                          psi.prior=Psiprior,
                          b.prior= Bprior,
                          BmsyK = BMSYK,
                          shape.CV=Shape.CV,
                          sigma.est = FALSE,
                          sigma.proc=Proc.error.JABBA,  # Fixed Process error (set sigma.est to TRUE if estimating) 
                          fixed.obsE=Obs.Error) 
  }
  
  # Fit JABBA
  output=fit_jabba(jbinput=jbinput,
                   ni = Sims,   #number of iterations    
                   nt = thinning,   #thinning interval
                   nb = burn.in,  #burn-in
                   nc = nchains,    #number of chains
                   init.values = FALSE,    #set to TRUE if specifying initial value
                   init.K = Kprior[1],
                   init.r = Rprior[1],
                   init.q = 1e-3,
                   save.all=TRUE,
                   save.csvs=TRUE,
                   output.dir=output.dir)
  
  
  #fit diagnostics
  setwd(output.dir)
  pdf(paste(outfile,'.pdf',sep=''))
  jbplot_catch(output)
  jbplot_ppdist(output)
  jbplot_mcmc(output)
  jbplot_procdev(output)
  jbplot_bprior(output)
  jbplot_trj(output,type="BBmsy",add=T)
  jbplot_trj(output,type="FFmsy",add=T)
  
  # status summary
  par(mfrow=c(3,2),mar = c(4, 5, 0.5, 0.1))
  jbplot_trj(output,type="B",add=T)
  jbplot_trj(output,type="F",add=T)
  jbplot_trj(output,type="BBmsy",add=T)
  jbplot_trj(output,type="FFmsy",add=T)
  jbplot_spphase(output,add=T)
  jbplot_kobe(output,add=T)
  dev.off()
  
  p1=output$pars_posterior%>%
    ggplot(aes(r))+
    geom_histogram()
  p2=output$pars_posterior%>%
    ggplot(aes(K))+
    geom_histogram()
  p3=output$pars_posterior%>%
    ggplot(aes(K,r))+
    geom_point() +
    geom_density_2d_filled(alpha = 0.5)+
    geom_density_2d(size = 0.25, colour = "black")
  ggarrange(plotlist = list(p1,p2,p3),nrow=1,widths = c(1,1,2))
  ggsave(paste(outfile,'.tiff',sep=''),width = 14,height = 8, dpi = 300, compression = "lzw")  
  
  #projections  
  if(!is.null(Projections))
  {
    Future.projections<-fw_jabba(jabba=output,
                                 nyears=length(Projections$Year),
                                 quant="Catch",
                                 type="abs",
                                 imp.values=Projections$Total,
                                 AR1=FALSE)
    Future.projections=Future.projections%>%
      filter(year%in%Projections$Year)
      #Biomass trajectories
    adumi=Future.projections%>%
            dplyr::select(year,B,iter)%>%
            distinct()%>%
            spread(year,B)%>%
            dplyr::select(-iter)%>%
            as.matrix()
    colnames(adumi)=NULL
    Bdumi=adumi
    output$posteriors$SB=cbind(output$posteriors$SB,adumi)
      #Depletion trajectories
    adumi=Future.projections%>%
      dplyr::select(year,BB0,iter)%>%
      distinct()%>%
      spread(year,BB0)%>%
      dplyr::select(-iter)%>%
      as.matrix()
    colnames(adumi)=NULL
    output$posteriors$P=cbind(output$posteriors$P,adumi) 
      #fishing mortality trajectories
    adumi=Future.projections%>%
      dplyr::select(year,H,iter)%>%
      distinct()%>%
      spread(year,H)%>%
      dplyr::select(-iter)%>%
      as.matrix()
    colnames(adumi)=NULL
    Fdumi=adumi
    output$posteriors$H=cbind(output$posteriors$H,adumi) 
      #B.over.Bmsy
    Mat=Bdumi
    Vec=output$posteriors$SBmsy
    adumi=Mat/matrix(rep(Vec,ncol(Mat)),ncol=ncol(Mat)) 
    output$posteriors$BtoBmsy=cbind(output$posteriors$BtoBmsy,adumi)
      #F.over.Fmsy
    Mat=Fdumi
    Vec=output$posteriors$Hmsy
    adumi=Mat/matrix(rep(Vec,ncol(Mat)),ncol=ncol(Mat))
    output$posteriors$HtoHmsy=cbind(output$posteriors$HtoHmsy,adumi)

    output$yr=c(output$yr,Projections$Year)
  }
  
  return(list(fit=output,jbinput=jbinput))
}
F.from.U=function(U) -log(1-U) 
jbplot_cpuefitsNEW=function (jabba, index = NULL, output.dir = getwd(), add = FALSE, 
                             as.png = FALSE, single.plots = add, width = NULL, height = NULL, 
                             plotCIs = TRUE, verbose = TRUE) 
{
  if (as.png == TRUE) 
    add = FALSE
  if (jabba$settings$CatchOnly == FALSE | jabba$settings$Auxiliary == 
      TRUE) {
    if (verbose) 
      cat(paste0("\n", "><> jbplot_cpue() - fits to CPUE <><", 
                 "\n"))
    if (jabba$settings$Auxiliary) {
      if (jabba$settings$CatchOnly) {
        jabba$settings$nI = jabba$settings$nA
        jabba$settings$I = as.matrix(jabba$settings$A)
        jabba$settings$SE2 = as.matrix(jabba$settings$A.SE2)
        di = dim(jabba$cpue.ppd)[3]
        jabba$cpue.ppd[, , 1:(di - 1)] = jabba$cpue.ppd[, 
                                                        , 2:(di)]
        jabba$cpue.hat[, , 1:(di - 1)] = jabba$cpue.hat[, 
                                                        , 2:(di)]
      }
      else {
        jabba$settings$nI = jabba$settings$nI + jabba$settings$nA
        jabba$settings$I = cbind(jabba$settings$I, as.matrix(jabba$settings$A))
        jabba$settings$SE2 = cbind(jabba$settings$SE2, 
                                   as.matrix(jabba$settings$A.SE2))
      }
    }
    if (is.null(index)) 
      index = 1:jabba$settings$nI
    N = jabba$settings$N
    years = jabba$yr
    indices = unique(jabba$diags$name)[index]
    CPUE = jabba$settings$I
    n.indices = length(indices)
    series = 1:n.indices
    if (single.plots == TRUE) {
      if (is.null(width)) 
        width = 5
      if (is.null(height)) 
        height = 3.5
      for (i in 1:n.indices) {
        Par = list(mfrow = c(1, 1), mar = c(3.5, 3.5, 
                                            0.5, 0.1), mgp = c(2, 0.5, 0), tck = -0.02, 
                   cex = 0.8)
        if (as.png == TRUE) {
          png(file = paste0(output.dir, "/Fits", jabba$assessment, 
                            "_", jabba$scenario, "_", indices[i], ".png"), 
              width = width, height = height, res = 200, 
              units = "in")
        }
        if (!add) {
          if (as.png == TRUE | i == 1) 
            par(Par)
        }
        Yr = jabba$yr
        Yr = min(Yr):max(Yr)
        yr = Yr - min(years) + 1
        fit = t(jabba$cpue.ppd[, c(2, 1, 3), index[i]])
        fit.hat = t(jabba$cpue.hat[, c(2, 1, 3), index[i]])
        mufit = mean(fit[2, ])
        fit = fit
        fit.hat = fit.hat
        cpue.i = CPUE[is.na(CPUE[, index[i]]) == F, index[i]]
        yr.i = Yr[is.na(CPUE[, index[i]]) == F]
        se.i = sqrt(jabba$settings$SE2[is.na(CPUE[, index[i]]) == 
                                         F, index[i]])
        ylim = c(min(fit * 0.9, exp(log(cpue.i) - 1.96 * 
                                      se.i)), max(fit * 1.05, exp(log(cpue.i) + 1.96 * 
                                                                    se.i)))
        cord.x <- c(Yr, rev(Yr))
        cord.y <- c(fit[1, yr], rev(fit[3, yr]))
        cord.yhat <- c(fit.hat[1, yr], rev(fit.hat[3, 
                                                   yr]))
        plot(years, CPUE[, index[i]], ylab = "Index", 
             xlab = "Year", ylim = ylim, xlim = range(jabba$yr), 
             type = "n", xaxt = "n", yaxt = "n")
        axis(1, labels = TRUE, cex = 0.8)
        axis(2, labels = TRUE, cex = 0.8)
        polygon(cord.x, cord.y, col = grey(0.5, 0.5), 
                border = 0, lty = 2)
        polygon(cord.x, cord.yhat, col = grey(0.3, 0.5), 
                border = grey(0.3, 0.5), lty = 2)
        lines(Yr, fit[2, yr], lwd = 2, col = 1)
        if (plotCIs) {
          upper <- qlnorm(0.975, meanlog = log(cpue.i), 
                          sdlog = se.i)
          lower <- qlnorm(0.025, meanlog = log(cpue.i), 
                          sdlog = se.i)
          arrows(x0 = yr.i, y0 = lower, x1 = yr.i, y1 = upper, 
                 length = 0.02, angle = 90, code = 3, col = 1)
        }
        points(yr.i, cpue.i, pch = 21, xaxt = "n", yaxt = "n", 
               bg = "white")
        legend("top", paste(indices[i]), bty = "n", y.intersp = -0.2, 
               cex = 0.9)
        if (as.png == TRUE) 
          dev.off()
      }
    }
    else {
      if (is.null(width)) 
        width = 7
      if (is.null(height)) 
        height = ifelse(n.indices == 1, 5, ifelse(n.indices == 
                                                    2, 3, 2.5)) * round(n.indices/2 + 0.01, 0)
      Par = list(mfrow = c(round(n.indices/2 + 0.01, 0), 
                           ifelse(n.indices == 1, 1, 2)), mai = c(0.35, 
                                                                  0.15, 0, 0.15), omi = c(0.2, 0.25, 0.2, 0) + 
                   0.1, mgp = c(2, 0.5, 0), tck = -0.02, cex = 0.8)
      if (as.png == TRUE) {
        png(file = paste0(output.dir, "/Fits_", jabba$assessment, 
                          "_", jabba$scenario, ".png"), width = 7, height = ifelse(n.indices == 
                                                                                     1, 5, ifelse(n.indices == 2, 3, 2.5)) * round(n.indices/2 + 
                                                                                                                                     0.01, 0), res = 200, units = "in")
      }
      if (!add) 
        par(Par)
      for (i in 1:n.indices) {
        Yr = jabba$yr
        Yr = min(Yr):max(Yr)
        yr = Yr - min(years) + 1
        fit = t(jabba$cpue.ppd[, c(2, 1, 3), index[i]])
        fit.hat = t(jabba$cpue.hat[, c(2, 1, 3), index[i]])
        mufit = mean(fit[2, ])
        fit = fit
        fit.hat = fit.hat
        cpue.i = CPUE[is.na(CPUE[, index[i]]) == F, index[i]]
        yr.i = Yr[is.na(CPUE[, index[i]]) == F]
        se.i = sqrt(jabba$settings$SE2[is.na(CPUE[, index[i]]) == 
                                         F, index[i]])
        ylim = c(min(fit * 0.9, exp(log(cpue.i) - 1.96 * 
                                      se.i)), max(fit * 1.05, exp(log(cpue.i) + 1.96 * 
                                                                    se.i)))
        cord.x <- c(Yr, rev(Yr))
        cord.y <- c(fit[1, yr], rev(fit[3, yr]))
        cord.yhat <- c(fit.hat[1, yr], rev(fit.hat[3, 
                                                   yr]))
        plot(years, CPUE[, index[i]], ylab = "", xlab = "", 
             ylim = ylim, xlim = range(jabba$yr), type = "n", 
             xaxt = "n", yaxt = "n")
        axis(1, labels = TRUE, cex = 0.8)
        axis(2, labels = TRUE, cex = 0.8)
        polygon(cord.x, cord.y, col = grey(0.5, 0.5), 
                border = 0, lty = 2)
        polygon(cord.x, cord.yhat, col = grey(0.3, 0.5), 
                border = grey(0.3, 0.5), lty = 2)
        lines(Yr, fit[2, yr], lwd = 2, col = 1)
        if (plotCIs) {
          upper <- qlnorm(0.975, meanlog = log(cpue.i), 
                          sdlog = se.i)
          lower <- qlnorm(0.025, meanlog = log(cpue.i), 
                          sdlog = se.i)
          arrows(x0 = yr.i, y0 = lower, x1 = yr.i, y1 = upper, 
                 length = 0.02, angle = 90, code = 3, col = 1)
        }
        points(yr.i, cpue.i, pch = 21, xaxt = "n", yaxt = "n", 
               bg = "white")
        legend("top", paste(indices[i]), bty = "n", y.intersp = -0.2, 
               cex = 0.9)
      }
      mtext(paste("Year"), side = 1, outer = TRUE, at = 0.5, 
            line = 1, cex = 1)
      mtext(paste("Index"), side = 2, outer = TRUE, at = 0.5, 
            line = 1, cex = 1)
      if (as.png == TRUE) {
        dev.off()
      }
    }
  }
  else {
    if (verbose) 
      cat(paste0("\n", "><> jbplot_cpuefits() not available CatchOnly=TRUE <><", 
                 "\n"))
  }
}
# Catch-only ensemble  ------------------------------------------------------
SPM=function(K,r,Init.dep,HRscen)
{
  HarvestRate=Inputd[,HRscen]
  HarvestRate=HarvestRate+rnorm(length(YRS),0,.01)
  HarvestRate=ifelse(HarvestRate<0,0,HarvestRate)
  
  Ktch=rep(NA,length(YRS))
  Bt=rep(NA,length(Ktch))
  Bt[1]=K*Init.dep
  Ktch[1]=HarvestRate[1]*Bt[1]
  for(t in 2:length(Bt))
  {
    Bt[t]=Bt[t-1]+r*Bt[t-1]*(1-Bt[t-1]/K)-Ktch[t-1]
    Ktch[t]=HarvestRate[t]*Bt[t]
  }
  
  U=Ktch/Bt
  Ft=-log(1-U)
  Bmsy=K/2
  MSY=K*r/4
  Fmsy=r/2
  Depletion=Bt/K
  B.Bmsy=Bt/Bmsy
  F.Fmsy=Ft/Fmsy
  #ii=(length(YRS)-4):length(YRS)
  
  return(list(Ktch=Ktch, r=r, K=K, Init.dep=Init.dep,HRscen=HRscen,
              Depletion=Depletion, B.Bmsy=B.Bmsy, F.Fmsy=F.Fmsy))
}
SPM.fit=function(PARS)
{
  K=exp(PARS[1])
  r=exp(PARS[2])
  Init.dep=exp(PARS[3])
  q=exp(PARS[4])
  
  HarvestRate=Inputd[,HRscen]
  HarvestRate=HarvestRate+rnorm(length(YRS),0,.01)
  HarvestRate=ifelse(HarvestRate<0,0,HarvestRate)
  
  Ktch=rep(NA,length(YRS))
  Bt=rep(NA,length(Ktch))
  Bt[1]=K*Init.dep
  Ktch[1]=HarvestRate[1]*Bt[1]
  for(t in 2:length(Bt))
  {
    Bt[t]=Bt[t-1]+r*Bt[t-1]*(1-Bt[t-1]/K)-Ktch[t-1]
    Ktch[t]=HarvestRate[t]*Bt[t]
  }
  
  Obs=length(Bt)
  Est_CPUE = q * Bt
  sqRes = (ln_CPUE - log(Est_CPUE)) * (ln_CPUE - log(Est_CPUE))
  sumsq = sum(sqRes)
  var = sumsq / Obs
  stdev = sqrt(var)
  # calculation of negative log-likelihood value - see Haddon (2001)
  NLL = (Obs / 2) * (log(2 * pi) + (2 * log(stdev)) + 1)
  
  if(what=='fit') out=NLL
  if(what=='simulate') out= list(cpue=Est_CPUE,ktch=Ktch,Bt=Bt)
  return(out)
}
fn.out.poly=function(In.dep)
{
  CI_lower=Dumi%>%
    filter(Init.dep==In.dep)%>%
    dplyr::select(Year,Grup,Depletion,Depletion_DBSRA,Depletion_CMSY,Depletion_JABBA,Depletion_SSS)%>%
    group_by(Year,Grup)%>%
    summarise(across(everything(), function(x)quantile(x,probs=0.001)))
  CI_upper=Dumi%>%
    filter(Init.dep==In.dep)%>%
    dplyr::select(Year,Grup,Depletion,Depletion_DBSRA,Depletion_CMSY,Depletion_JABBA,Depletion_SSS)%>%
    group_by(Year,Grup)%>%
    summarise(across(everything(), function(x)quantile(x,probs=0.999)))
  CI_lower%>%
    left_join(CI_upper,by=c('Year','Grup'))%>%
    arrange(Year,Grup)%>%
    ggplot(aes(Year,Depletion.x))+
    facet_wrap(~Grup,ncol=3)+
    theme(legend.position = 'none')+
    geom_ribbon(aes(ymin=Depletion_DBSRA.x,ymax=Depletion_DBSRA.y,fill="DBSRA"),alpha=0.5)+
    geom_ribbon(aes(ymin=Depletion_CMSY.x,ymax=Depletion_CMSY.y,fill="CMSY"),alpha=0.5)+
    geom_ribbon(aes(ymin=Depletion_JABBA.x,ymax=Depletion_JABBA.y,fill="JABBA"),alpha=0.5)+
    geom_ribbon(aes(ymin=Depletion_SSS.x,ymax=Depletion_SSS.y,fill="SSS"),alpha=0.5)+
    geom_ribbon(aes(ymin=Depletion.x,ymax=Depletion.y,fill='Operation model'),alpha=0.8)+
    theme_PA()+ylab('Depletion')+
    scale_fill_manual(values = le.cols)+
    theme(legend.position = 'top',
          legend.title = element_blank())
  
  ggsave(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),
               paste('Compare COMs performance_Init.dep_',In.dep,'.tiff',sep=''),sep=''),
         width = 10,height = 10,compression = "lzw")
}
prob.wrapper=function(DAT,id.yr,B.target,B.threshold,B.limit)
{
  f=ecdf(DAT[,id.yr])
  P.below.target=f(B.target)
  P.below.threshold=f(B.threshold)
  P.below.limit=f(B.limit)
  P.above.target=1-P.below.target
  P.above.threshold=1-P.below.threshold
  P.above.limit=1-P.below.limit
  P.between.thre.tar=P.below.target-P.below.threshold
  P.between.lim.thre=P.below.threshold-P.below.limit
  return(data.frame(Range=c('<lim','lim.thr',
                           'thr.tar','>tar'),
                   Probability=round(c(P.below.limit,P.between.lim.thre,
                                       P.between.thre.tar,P.above.target),3)))
}
mod.average=function(dd,Weights,Ref.pnt=NULL)
{
  yrs=as.numeric(gsub("^.*X","",names(dd$CMSY)))
  NN=max(sapply(dd,nrow))
  for(w in 1:length(dd))
  {
    WEI=Weights%>%filter(Model==names(dd)[w])%>%pull(Weight)
    Size=round(WEI*NN)
    dd[[w]]=dd[[w]][sample(1:nrow(dd[[w]]),Size,replace = T),]
    colnames(dd[[w]])=1:ncol(dd[[w]])
  }
  dd=do.call(rbind,dd)
  Median=apply(dd,2,median,na.rm=T)
  Lower=apply(dd,2,function(x) quantile(x,probs=0.025,na.rm=T))
  Upper=apply(dd,2,function(x) quantile(x,probs=0.975,na.rm=T))
  Lower.60=apply(dd,2,function(x) quantile(x,probs=0.2,na.rm=T))
  Upper.60=apply(dd,2,function(x) quantile(x,probs=0.8,na.rm=T))
  Trajectories=data.frame(year=yrs,median=Median,
                          lower.60=Lower.60,upper.60=Upper.60,
                          lower.95=Lower,upper.95=Upper)
  
  Probs=Probs.future=NULL
  if(!is.null(Ref.pnt))
  {
    Years=Trajectories$year
    Probs=prob.wrapper(DAT=dd,
                       id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                       B.target=Ref.pnt%>%filter(Rf.pt=='Target')%>%pull(Value),
                       B.threshold=Ref.pnt%>%filter(Rf.pt=='Threshold')%>%pull(Value),
                       B.limit=Ref.pnt%>%filter(Rf.pt=='Limit')%>%pull(Value))%>%
                       mutate(Probability=Probability/sum(Probability),
                              finyear=as.numeric(substr(Last.yr.ktch,1,4)))
    
    Probs.future=prob.wrapper(DAT=dd,
                              id.yr=match(Years[length(Years)],Years),
                              B.target=Ref.pnt%>%filter(Rf.pt=='Target')%>%pull(Value),
                              B.threshold=Ref.pnt%>%filter(Rf.pt=='Threshold')%>%pull(Value),
                              B.limit=Ref.pnt%>%filter(Rf.pt=='Limit')%>%pull(Value))%>%
                      mutate(Probability=Probability/sum(Probability),
                             finyear=Years[length(Years)])
  }
  return(list(Trajectories=Trajectories,Probs=Probs,Probs.future=Probs.future))
}
mod.average.scalar=function(dd,Weights,KtcH)
{
  NN=max(sapply(dd,nrow))
  for(w in 1:length(dd))
  {
    WEI=Weights%>%filter(Model==names(dd)[w])%>%pull(Weight)
    Size=round(WEI*NN)
    dd[[w]]=dd[[w]]%>%
              filter(Scenario=='S1')%>%
              dplyr::select(Value)
    dd[[w]]=dd[[w]][sample(1:nrow(dd[[w]]),Size,replace = T),]
  }
  dd=do.call(c,dd)  
  Median=median(dd,na.rm=T)
  Lower=quantile(dd,probs=0.025,na.rm=T)
  Upper=quantile(dd,probs=0.975,na.rm=T)
  Lower.60=quantile(dd,probs=0.2,na.rm=T)
  Upper.60=quantile(dd,probs=0.8,na.rm=T)
  Threshold=Upper.60
  Limit=Upper
  Target=Median
  Catch.mean.last.n.yrs=mean(data.table::last(KtcH$Catch,n.last.catch.yrs_MSY.catch.only))
  f=ecdf(dd) 
  P.catch=f(Catch.mean.last.n.yrs)
  P.below.target=f(Target)
  P.below.threshold=f(Threshold)  
  P.below.limit=f(Limit)
  #Minor consequence (1): catch<tar; Moderate (2): tar<catch<thre; High (3): thre<catch<Limit; Major (4): Limit<catch
  P.catch.below.tar=ifelse(Catch.mean.last.n.yrs<=Target,1,0)
  P.between.thre.tar=ifelse(Target<Catch.mean.last.n.yrs & Catch.mean.last.n.yrs<=Threshold,1,0)
  P.between.lim.thre=ifelse(Threshold<Catch.mean.last.n.yrs & Catch.mean.last.n.yrs<=Limit,1,0)
  P.catch.above.limit=ifelse(Limit<Catch.mean.last.n.yrs,1,0)
  

  return(data.frame(median=Median,
                    lower.60=Lower.60,upper.60=Upper.60,
                    lower.95=Lower,upper.95=Upper,
                    Catch.mean.last.n.yrs=Catch.mean.last.n.yrs,
                    Threshold=Threshold,
                    Limit=Limit,
                    Target=Target,
                    P.catch.below.tar=P.catch.below.tar,
                    P.between.thre.tar=P.between.thre.tar,
                    P.between.lim.thre=P.between.lim.thre,
                    P.catch.above.limit=P.catch.above.limit,
                    P.catch=P.catch,
                    P.below.target=P.below.target,
                    P.below.threshold=P.below.threshold,
                    P.below.limit=P.below.limit))
}


# Display functions  ------------------------------------------------------
fn.prior=function(N=1e4,d,MAX=NULL)
{
  if(d$dist=='unif')
  {
    out=runif(n=N, min=d$low, max=d$up) 
  }
  if(d$dist=='lnorm')
  {
    out=rlnorm(n=N, meanlog=d$mean, sdlog=d$sd)
    if(!is.null(MAX)) out=subset(out,out<=MAX)
  } 
  if(d$dist=='beta')
  {
    shape.pars=get_beta(d$mean,d$sd)
    out=rbeta(n=N, shape1=shape.pars[1], shape2=shape.pars[2])
  }
  return(out)
}
fn.show.density=function(d,NCOL)
{
  d%>%
    ggplot(aes(x=Value,fill=Distribuion))+
    geom_density(position="identity",size=1.15,alpha=0.65)+
    facet_wrap(~Parameter,scales='free',ncol=NCOL)+
    theme_PA(strx.siz=18,leg.siz=18,axs.t.siz=16,axs.T.siz=26)+
    theme(legend.position="top",
          legend.title = element_blank(),
          legend.key=element_blank(),
          title=element_text(size=12))+
    ylab("Density")+xlab("Parameter value")
}
add.probs=function(DAT,id.yr,B.threshold,plot.ranges=FALSE) #Reference points probability  
{
  DAT=DAT[,id.yr]
  B.target=Tar.prop.bmsny*B.threshold
  B.limit=max(Minimum.acceptable.blim,Lim.prop.bmsy*B.threshold)
  f=ecdf(DAT)
  P.below.target=f(B.target)
  P.below.threshold=f(B.threshold)
  P.below.limit=f(B.limit)
  P.above.target=1-P.below.target
  P.above.threshold=1-P.below.threshold
  P.above.limit=1-P.below.limit
  P.between.thre.tar=P.below.target-P.below.threshold
  P.between.lim.thre=P.below.threshold-P.below.limit
  if(plot.ranges)
  {
    if(P.above.target>0)
    {
      segments(YR[id.yr],B.target,YR[id.yr],1,col=CL.ref.pt[1],lwd=8,lend="butt")
      Legn=round(100*P.above.target)
      if(Legn==0)Legn="<1"
      text(YR[id.yr],mean(c(B.target,UP[id.yr])),paste(Legn,"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
    if(P.between.thre.tar>0)
    {
      Upseg=B.target
      Lwseg=B.threshold
      segments(YR[id.yr],Upseg,YR[id.yr],Lwseg,col=CL.ref.pt[2],lwd=8,lend="butt")
      Legn=round(100*P.between.thre.tar)
      if(Legn==0)Legn="<1"
      text(YR[id.yr],mean(c(Upseg,Lwseg))*1.025,paste(Legn,"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
    if(P.between.lim.thre>0)
    {
      Upseg=B.threshold
      Lowseg=B.limit
      segments(YR[id.yr],Upseg,YR[id.yr],Lowseg,col=CL.ref.pt[3],lwd=8,lend="butt")
      Legn=round(100*P.between.lim.thre)
      if(Legn==0)Legn="<1"
      wher.txt=mean(c(Upseg,Lowseg))*1.025
      if(wher.txt>0.5) wher.txt=0.5*.9
      text(YR[id.yr],wher.txt,paste(Legn,"%",sep=""),
           col="black",cex=CEX,srt=SRT,font=2,pos=2)
    }
    if(P.below.limit>0)
    {
      segments(YR[id.yr],B.limit,YR[id.yr],0,col=CL.ref.pt[4],lwd=8,lend="butt")
      Legn=round(100*P.below.limit)
      if(Legn==0)Legn="<1"
      
      text(YR[id.yr],B.limit*0.85,paste(Legn,"%",sep=""),
           col="black",cex=CEX,srt=SRT,pos=2,font=2)
    }
    
  }
  
  return(list(probs=data.frame(Range=c('<lim','lim.thr',
                                       'thr.tar','>tar'),
                               Probability=round(c(P.below.limit,P.between.lim.thre,
                                                   P.between.thre.tar,P.above.target),3)),
              Reference.points=data.frame(Rf.pt=c('Target','Threshold','Limit'),
                                          Value=c(B.target,B.threshold,B.limit))))
}
fn.ktch.only.get.timeseries=function(d,mods,Type,add.50=FALSE,scen,Katch)
{
  if('output'%in%names(d))dd=d$output

  if(mods=='DBSRA')
  {
    Years=dd$Years
    if(Type=='Depletion')
    {
      d1=dd$Depletion.traj[1:length(Years)]
      Probs=add.probs(DAT=d1,
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=dd$Estimates[Biomass.threshold,'Median (ll=1)']/dd$Estimates['K','Median (ll=1)'])
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
      if(as.numeric(substr(Last.yr.ktch,1,4))<Years[length(Years)])  
      {
        Probs.future=add.probs(DAT=d1,
                               id.yr=match(Years[length(Years)],Years),
                               B.threshold=dd$Estimates[Biomass.threshold,'Median (ll=1)']/dd$Estimates['K','Median (ll=1)'])
        Probs$probs.future=Probs.future$probs%>%
          mutate(Scenario=scen,
                 finyear=Years[length(Years)])
      }
      
    }
    if(Type=='F.series')d1=dd$F.series[1:length(Years)]
    if(Type=='B.Bmsy')
    {
      d1=dd$B.Bmsy[1:length(Years)]
      Probs=add.probs(DAT=d1,
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=1)
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
      if(as.numeric(substr(Last.yr.ktch,1,4))<Years[length(Years)])  
      {
        Probs.future=add.probs(DAT=d1,
                               id.yr=match(Years[length(Years)],Years),
                               B.threshold=1)
        Probs$probs.future=Probs.future$probs%>%
          mutate(Scenario=scen,
                 finyear=Years[length(Years)])
      }
    }
    if(Type=='F.Fmsy') d1=dd$F.Fmsy[1:length(Years)]
    Dat=data.frame(year=Years,
                   median=apply(d1,2,median),
                   lower.60=apply(d1,2,function(x)quantile(x,probs=0.2,na.rm=T)),
                   upper.60=apply(d1,2,function(x)quantile(x,probs=0.8,na.rm=T)),
                   upper.95=apply(d1,2,function(x)quantile(x,probs=0.975,na.rm=T)),
                   lower.95=apply(d1,2,function(x)quantile(x,probs=0.025,na.rm=T)))%>%
      mutate(Model=mods,
             Catch=Katch)
    if(add.50)
    {
      Dat=Dat%>%
        mutate(upper.50=apply(d1,2,function(x)quantile(x,probs=0.75,na.rm=T)),
               lower.50=apply(d1,2,function(x)quantile(x,probs=0.25,na.rm=T)))
    }
  }
  
  if(mods=='CMSY')
  {
    Years=dd$Years
    if(Type=='Depletion')
    {
      d1=dd$Depletion.traj[1:length(Years)]
      Probs=add.probs(DAT=d1,
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=median(dd$Bmsy)/dd$Statistics$output['K','50%'])
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
      if(as.numeric(substr(Last.yr.ktch,1,4))<Years[length(Years)])  
      {
        Probs.future=add.probs(DAT=d1,
                        id.yr=match(Years[length(Years)],Years),
                        B.threshold=median(dd$Bmsy)/dd$Statistics$output['K','50%'])
        Probs$probs.future=Probs.future$probs%>%
                            mutate(Scenario=scen,
                                   finyear=Years[length(Years)])
      }
    }
    
    if(Type=='F.series')  d1=dd$F.traj[,1:length(Years)]
    if(Type=='B.Bmsy')
    {
      d1=dd$B.Bmsy[1:length(Years)]
      Probs=add.probs(DAT=d1,
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=1)
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
      if(as.numeric(substr(Last.yr.ktch,1,4))<Years[length(Years)])  
      {
        Probs.future=add.probs(DAT=d1,
                               id.yr=match(Years[length(Years)],Years),
                               B.threshold=1)
        Probs$probs.future=Probs.future$probs%>%
          mutate(Scenario=scen,
                 finyear=Years[length(Years)])
      }
    }
      
    if(Type=='F.Fmsy') d1=dd$F.Fmsy[,1:length(Years)]
    
    Dat=data.frame(year=as.numeric(Years),
                   median=apply(d1,2,median),
                   lower.60=apply(d1,2,function(x)quantile(x,probs=0.2,na.rm=T)),
                   upper.60=apply(d1,2,function(x)quantile(x,probs=0.8,na.rm=T)),
                   upper.95=apply(d1,2,function(x)quantile(x,probs=0.975,na.rm=T)),  
                   lower.95=apply(d1,2,function(x)quantile(x,probs=0.025,na.rm=T)))%>%
      mutate(Model=mods,
             Catch=Katch) 
    if(add.50)
    {
      Dat=Dat%>%
        mutate(upper.50=apply(d1,2,function(x)quantile(x,probs=0.75,na.rm=T)),
               lower.50=apply(d1,2,function(x)quantile(x,probs=0.25,na.rm=T)))
    }
  }
  
  if(mods=='JABBA')
  {
    Years=dd$yr
    
    if(Type=='Depletion') 
    {
      K=dd$pars
      K=K[match("K",rownames(K)),"Median"]
      d1=data.frame(mu=apply(dd$posteriors$P,2,median,na.rm=T),
                    lci.60=apply(dd$posteriors$P,2,function(x) quantile(x,probs=0.2,na.rm=T)),
                    uci.60=apply(dd$posteriors$P,2,function(x) quantile(x,probs=0.8,na.rm=T)),
                    lci=apply(dd$posteriors$P,2,function(x) quantile(x,probs=0.025,na.rm=T)),
                    uci=apply(dd$posteriors$P,2,function(x) quantile(x,probs=0.975,na.rm=T)))
      Probs=add.probs(DAT=sweep(dd$posteriors$SB,1,dd$posteriors$K,'/'),
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=median(dd$posteriors$SBmsy)/K)
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
      
      if(as.numeric(substr(Last.yr.ktch,1,4))<Years[length(Years)])  
      {
        Probs.future=add.probs(DAT=sweep(dd$posteriors$SB,1,dd$posteriors$K,'/'),
                               id.yr=match(Years[length(Years)],Years),
                               B.threshold=median(dd$posteriors$SBmsy)/K)
        Probs$probs.future=Probs.future$probs%>%
          mutate(Scenario=scen,
                 finyear=Years[length(Years)])
      }
    }
    if(Type=='F.series')
    {
      d1=data.frame(mu=apply(dd$posteriors$H,2,median,na.rm=T),
                    lci.60=apply(dd$posteriors$H,2,function(x) quantile(x,probs=0.2,na.rm=T)),
                    uci.60=apply(dd$posteriors$H,2,function(x) quantile(x,probs=0.8,na.rm=T)),
                    lci=apply(dd$posteriors$H,2,function(x) quantile(x,probs=0.025,na.rm=T)),
                    uci=apply(dd$posteriors$H,2,function(x) quantile(x,probs=0.975,na.rm=T)))
    }
    if(Type=='B.Bmsy')
    {
      d1=data.frame(mu=apply(dd$posteriors$BtoBmsy,2,median,na.rm=T),
                    lci.60=apply(dd$posteriors$BtoBmsy,2,function(x) quantile(x,probs=0.2,na.rm=T)),
                    uci.60=apply(dd$posteriors$BtoBmsy,2,function(x) quantile(x,probs=0.8,na.rm=T)),
                    lci=apply(dd$posteriors$BtoBmsy,2,function(x) quantile(x,probs=0.025,na.rm=T)),
                    uci=apply(dd$posteriors$BtoBmsy,2,function(x) quantile(x,probs=0.975,na.rm=T)))
      Probs=add.probs(DAT=sweep(dd$posteriors$SB,1,dd$posteriors$SBmsy,'/'),
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=1)
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
      if(as.numeric(substr(Last.yr.ktch,1,4))<Years[length(Years)])  
      {
        Probs.future=add.probs(DAT=sweep(dd$posteriors$SB,1,dd$posteriors$SBmsy,'/'),
                               id.yr=match(Years[length(Years)],Years),
                               B.threshold=1)
        Probs$probs.future=Probs.future$probs%>%
          mutate(Scenario=scen,
                 finyear=Years[length(Years)])
      }
    }
    if(Type=='F.Fmsy')
    {
      d1=data.frame(mu=apply(dd$posteriors$HtoHmsy,2,median,na.rm=T),
                    lci.60=apply(dd$posteriors$HtoHmsy,2,function(x) quantile(x,probs=0.2,na.rm=T)),
                    uci.60=apply(dd$posteriors$HtoHmsy,2,function(x) quantile(x,probs=0.8,na.rm=T)),
                    lci=apply(dd$posteriors$HtoHmsy,2,function(x) quantile(x,probs=0.025,na.rm=T)),
                    uci=apply(dd$posteriors$HtoHmsy,2,function(x) quantile(x,probs=0.975,na.rm=T)))
    }
    
    
    Dat=data.frame(year=as.numeric(Years),
                   median=d1$mu[1:length(Years)],
                   lower.60=d1$lci.60[1:length(Years)],
                   upper.60=d1$uci.60[1:length(Years)],
                   upper.95=d1$uci[1:length(Years)],  
                   lower.95=d1$lci[1:length(Years)])%>%
      mutate(Model=mods,
             Catch=Katch) 
    
    if(add.50)
    {
      Dat=Dat%>%
        mutate(upper.50=NA,
               lower.50=NA)
    }
    
  }
  
  if(mods%in%c('SSS','SS3'))
  {
    if(mods=='SSS') biom.type='total'
    if(mods=='SS3') biom.type='spawning'
    
    Years=with(d[[1]],startyr:endyr)
    if(d[[1]]$N_forecast_yrs>1) Years=with(d[[1]],startyr:(endyr+d[[1]]$N_forecast_yrs))
    
    dum=do.call(rbind,fn.get.stuff.from.list(d,"derived_quants"))
    
    SSB_Virgin=dum[grep("SSB_Virgin",dum$Label),]
    SSB_Virgin=SSB_Virgin%>%
      mutate(Simulation=readr::parse_number(paste(rownames(SSB_Virgin),'0',sep='')))%>%
      rename(Value.virgin=Value)%>%
      dplyr::select(Value.virgin,Simulation)
    SSB_MSY=dum[grep("SSB_MSY",dum$Label),]
    SSB_MSY=SSB_MSY%>%
      mutate(Simulation=readr::parse_number(paste(rownames(SSB_MSY),'0',sep='')))%>%
      rename(Value.MSY=Value)%>%
      dplyr::select(Value.MSY,Simulation)%>%
      left_join(SSB_Virgin,by='Simulation')%>%
      mutate(Depletion=Value.MSY/Value.virgin)
    
    if(Type=='Depletion') 
    {
      if(biom.type=='spawning')
      {
        SSB=dum[grep(paste(paste("SSB",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
        SSB=SSB%>%
          mutate(year=readr::parse_number(Label),
                 Simulation=readr::parse_number(paste(rownames(SSB),'0',sep='')),
                 Simulation=as.numeric(substr(Simulation,5,10)))%>%
          dplyr::select(-Label)%>%
          left_join(SSB_Virgin,by='Simulation')%>%
          mutate(Depletion=Value/Value.virgin)%>%
          filter(year%in%Years)
        d1=SSB%>%
          group_by(year)%>%
          summarise(mu=median(Depletion,na.rm=T),
                    lci.60=quantile(Depletion,probs=0.2,na.rm=T),
                    uci.60=quantile(Depletion,probs=0.8,na.rm=T),
                    lci=quantile(Depletion,probs=0.025,na.rm=T),
                    uci=quantile(Depletion,probs=0.975,na.rm=T))
        Probs=add.probs(DAT=SSB%>%
                          dplyr::select(Simulation,year,Depletion)%>%
                          spread(year,Depletion)%>%
                          dplyr::select(-Simulation),
                        id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                        B.threshold=median(SSB_MSY$Depletion))
        Probs$probs=Probs$probs%>%mutate(Scenario=scen)
        
        if(as.numeric(substr(Last.yr.ktch,1,4))<Years[length(Years)])  
        {
          Probs.future=add.probs(DAT=SSB%>%
                                   dplyr::select(Simulation,year,Depletion)%>%
                                   spread(year,Depletion)%>%
                                   dplyr::select(-Simulation),
                                 id.yr=match(Years[length(Years)],Years),
                                 B.threshold=median(SSB_MSY$Depletion))
          Probs$probs.future=Probs.future$probs%>%
            mutate(Scenario=scen,
                   finyear=Years[length(Years)])
        }
        
      }
      if(biom.type=='total')
      {
        TotB=do.call(rbind,fn.get.stuff.from.list(d,"timeseries"))
        Virgin.B=TotB%>%
          filter(Era=='VIRG')%>%
          dplyr::select(Bio_all)%>%
          mutate(Simulation=row_number())%>%
          rename(Virgin.B=Bio_all)
        nnsim=length(unique(TotB$Yr))
        TotB=TotB%>%
          dplyr::select(Yr,Bio_all)%>%
          rename(year=Yr)%>%
          mutate(Simulation=rep(Virgin.B$Simulation,each=nnsim))%>%
          left_join(Virgin.B,by='Simulation')%>%
          mutate(Depletion=Bio_all/Virgin.B)%>%
          filter(year%in%Years)
        d1=TotB%>%
          group_by(year)%>%
          summarise(mu=median(Depletion,na.rm=T),
                    lci.60=quantile(Depletion,probs=0.2,na.rm=T),
                    uci.60=quantile(Depletion,probs=0.8,na.rm=T),
                    lci=quantile(Depletion,probs=0.025,na.rm=T),
                    uci=quantile(Depletion,probs=0.975,na.rm=T))
        Probs=add.probs(DAT=TotB%>%
                          dplyr::select(Simulation,year,Depletion)%>%
                          spread(year,Depletion)%>%
                          dplyr::select(-Simulation),
                        id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                        B.threshold=median(SSB_MSY$Depletion))
        Probs$probs=Probs$probs%>%mutate(Scenario=scen)
        
        if(as.numeric(substr(Last.yr.ktch,1,4))<Years[length(Years)])  
        {
          Probs.future=add.probs(DAT=TotB%>%
                                   dplyr::select(Simulation,year,Depletion)%>%
                                   spread(year,Depletion)%>%
                                   dplyr::select(-Simulation),
                                 id.yr=match(Years[length(Years)],Years),
                                 B.threshold=median(SSB_MSY$Depletion))
          Probs$probs.future=Probs.future$probs%>%
            mutate(Scenario=scen,
                   finyear=Years[length(Years)])
        }
      }
    }
    
    if(Type=='F.series')
    {
      EF=dum[grep(paste(paste("F",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
      d1=EF%>%
        mutate(year=readr::parse_number(Label))%>%
        group_by(year)%>%
        summarise(mu=median(Value,na.rm=T),
                  lci.60=quantile(Value,probs=0.2,na.rm=T),
                  uci.60=quantile(Value,probs=0.8,na.rm=T),
                  lci=quantile(Value,probs=0.025,na.rm=T),
                  uci=quantile(Value,probs=0.975,na.rm=T))
    }
    
    if(Type=='B.Bmsy') 
    {
      SSB=dum[grep(paste(paste("SSB",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
      SSB=SSB%>%
        mutate(year=readr::parse_number(Label),
               Simulation=readr::parse_number(paste(rownames(SSB),'0',sep='')),
               Simulation=as.numeric(substr(Simulation,5,10)))%>%
        dplyr::select(-Label)%>%
        left_join(SSB_MSY%>%dplyr::select(Value.MSY,Simulation),by='Simulation')%>%
        mutate(B_Bmys=Value/Value.MSY)%>%
        filter(year%in%Years)
        
        d1=SSB%>%
            group_by(year)%>%
            summarise(mu=median( B_Bmys,na.rm=T),
                      lci.60=quantile( B_Bmys,probs=0.2,na.rm=T),
                      uci.60=quantile( B_Bmys,probs=0.8,na.rm=T),
                      lci=quantile( B_Bmys,probs=0.025,na.rm=T),
                      uci=quantile( B_Bmys,probs=0.975,na.rm=T))
      
      Probs=add.probs(DAT=SSB%>%
                          dplyr::select(Simulation,year,B_Bmys)%>%
                          spread(year,B_Bmys)%>%
                          dplyr::select(-Simulation),
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=1)
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
      if(as.numeric(substr(Last.yr.ktch,1,4))<Years[length(Years)])  
      {
        Probs.future=add.probs(DAT=SSB%>%
                                 dplyr::select(Simulation,year,B_Bmys)%>%
                                 spread(year,B_Bmys)%>%
                                 dplyr::select(-Simulation),
                               id.yr=match(Years[length(Years)],Years),
                               B.threshold=1)
        Probs$probs.future=Probs.future$probs%>%
          mutate(Scenario=scen,
                 finyear=Years[length(Years)])
      }
    }
    
    if(Type=='F.Fmsy')
    {
      EF=dum[grep(paste(paste("F",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
      
      annF_MSY=dum[grep("annF_MSY",dum$Label),]
      annF_MSY=annF_MSY%>%
        mutate(Simulation=readr::parse_number(paste(rownames(annF_MSY),'0',sep='')))%>%
        rename(Value.F_MSY=Value)%>%
        dplyr::select(Value.F_MSY,Simulation)
      
      d1=EF%>%
        mutate(year=readr::parse_number(Label),
               Simulation=readr::parse_number(paste(rownames(EF),'0',sep='')),
               Simulation=as.numeric(substr(Simulation,5,10)))%>%
        dplyr::select(-Label)%>%
        left_join(annF_MSY,by='Simulation')%>%
        mutate(F_Fmsy=Value/Value.F_MSY)%>%
        group_by(year)%>%
        summarise(mu=median(F_Fmsy,na.rm=T),
                  lci.60=quantile(F_Fmsy,probs=0.2,na.rm=T),
                  uci.60=quantile(F_Fmsy,probs=0.8,na.rm=T),
                  lci=quantile(F_Fmsy,probs=0.025,na.rm=T),
                  uci=quantile(F_Fmsy,probs=0.975,na.rm=T))
      
    }
    
    Dat=data.frame(year=as.numeric(Years),
                   median=d1$mu[1:length(Years)],
                   lower.60=d1$lci.60[1:length(Years)],
                   upper.60=d1$uci.60[1:length(Years)],
                   upper.95=d1$uci[1:length(Years)],  
                   lower.95=d1$lci[1:length(Years)])%>%
      mutate(Model=mods,
             Catch=Katch) 
    
    if(add.50)
    {
      Dat=Dat%>%
        mutate(upper.50=NA,
               lower.50=NA)
    }
  }
  
  Dat=Dat%>%mutate(Scenario=scen)
  
  if(Type%in%c('Depletion','B.Bmsy'))
  {
    return(list(Dat=Dat,Probs=Probs))
  }else
  {
    return(list(Dat=Dat))
  }
  
}

add.probs.integrated.asympt.error=function(DAT,B.threshold) #Reference points probability based on asymptotic error 
{
  B.target=Tar.prop.bmsny*B.threshold
  B.limit=max(Minimum.acceptable.blim,Lim.prop.bmsy*B.threshold)
  f=ecdf(DAT)
  P.below.target=f(B.target)
  P.below.threshold=f(B.threshold)
  P.below.limit=f(B.limit)
  P.above.target=1-P.below.target
  P.above.threshold=1-P.below.threshold
  P.above.limit=1-P.below.limit
  P.between.thre.tar=P.below.target-P.below.threshold
  P.between.lim.thre=P.below.threshold-P.below.limit
  return(list(probs=data.frame(Range=c('<lim','lim.thr',
                                       'thr.tar','>tar'),
                               Probability=round(c(P.below.limit,P.between.lim.thre,
                                                   P.between.thre.tar,P.above.target),3)),
              Reference.points=data.frame(Rf.pt=c('Target','Threshold','Limit'),
                                          Value=c(B.target,B.threshold,B.limit))))
}
fn.integrated.mod.get.timeseries=function(d,mods,Type,add.50=FALSE,scen)
{
  if(mods=='SS3')
  {
    Years=with(d,startyr:endyr)  
    
    dum=d[["derived_quants"]]
    
    if(Type=='Depletion') 
    { 
      d1=dum[grep(paste(paste("Bratio",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value','StdDev')]%>%
        mutate(year=readr::parse_number(Label))%>%
        filter(year%in%Years)%>%
        rename(mu=Value)%>%
        mutate(lci=mu-1.96*StdDev,
               uci=mu+1.96*StdDev)%>%
        relocate(year,mu,lci,uci)%>%
        `rownames<-`( NULL )%>%
        rename(lower.95=lci,
               upper.95=uci)
      
      Probs=add.probs.integrated.asympt.error(DAT=d1[nrow(d1),]%>%dplyr::select(mu,StdDev),
                                              B.threshold=dum[grep('B_MSY/SSB_unfished',dum$Label),'Value'])     
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
      
      d1=d1%>%dplyr::select(-c(Label,StdDev))
    }
    
    if(Type=='F.series')
    {
      EF=dum[grep(paste(paste("F",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value','StdDev')]
      d1=EF%>%
          mutate(year=readr::parse_number(Label))%>%
          filter(year%in%Years)%>%
          rename(mu=Value)%>%
          mutate(lci=mu-1.96*StdDev,
                 uci=mu+1.96*StdDev)%>%
          relocate(year,mu,lci,uci)%>%
        `rownames<-`( NULL )%>%
        dplyr::select(-c(Label,StdDev))%>%
        rename(lower.95=lci,
               upper.95=uci)
      
    }
    
    if(Type=='B.Bmsy') 
    {
      SSB=dum[grep(paste(paste("SSB",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
      SSB_MSY=dum[dum$Label=="SSB_MSY",c('Label','Value')]
      d1=SSB%>%
        mutate(year=readr::parse_number(Label),
               mu=Value/SSB_MSY$Value)%>%
        relocate(year,mu)%>%
        `rownames<-`( NULL )%>%
        dplyr::select(-c(Label,Value))
      
    }
    
    if(Type=='F.Fmsy')
    {
      EF=dum[grep(paste(paste("F",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
      annF_MSY=dum[dum$Label=="annF_MSY",c('Label','Value')]
      d1=EF%>%
        mutate(year=readr::parse_number(Label),
               mu=Value/annF_MSY$Value)%>%
        relocate(year,mu)%>%
        `rownames<-`( NULL )%>%
        dplyr::select(-c(Label,Value))
    }
    
    Dat=d1%>%
      mutate(Model=mods) 
    
    if(add.50)
    {
      Dat=Dat%>%
        mutate(upper.50=NA,
               lower.50=NA)
    }
  }
  
  Dat=Dat%>%
          mutate(Scenario=scen)%>%
          rename(median=mu)
  
  if(Type=='Depletion')
  {
    return(list(Dat=Dat,Probs=Probs))
  }else
  {
    return(list(Dat=Dat))
  }
  
}

fn.display.priors=function(d,sp,XLAB,XLIM,Strx.siz=16)
{
  dummy=lapply(d[sp],function(x) rnorm(1e4,x$mean,x$sd))
  dummy=lapply(dummy,function(x)data.frame(var=x))
  dummy=do.call(rbind,dummy)%>%
    tibble::rownames_to_column(var='Species')%>%
    mutate(Species=capitalize(gsub("\\..*","",Species)))
  p=dummy%>%
    ggplot(aes(x=var))+
    geom_density(aes(color=Species),linewidth=1.5)+
    facet_wrap(~Species,scales='free_y')+
    xlab(XLAB)+ylab("Density")+
    theme_PA(axs.T.siz=22,axs.t.siz=14,strx.siz=Strx.siz)+
    theme(legend.position = "none",
          plot.title =element_text(size=17))+
    scale_x_continuous(limits = XLIM)
  
  return(p)
}
fn.display.prior.sensitivity=function(d,d2,sp,XLAB,Strx.siz=16,Scen1,Scen2)
{
  dummy=lapply(d[sp],function(x) rnorm(1e4,x$mean,x$sd))
  dummy=lapply(dummy,function(x)data.frame(var=x))
  
  dummy1=lapply(d2[sp],function(x) rnorm(1e4,x$mean,x$sd))
  dummy1=lapply(dummy1,function(x)data.frame(var=x))
  
  dummy=do.call(rbind,dummy)%>%
    tibble::rownames_to_column(var='Species')%>%
    mutate(Species=capitalize(gsub("\\..*","",Species)),
           Scenario=Scen1)
  
  dummy1=do.call(rbind,dummy1)%>%
    tibble::rownames_to_column(var='Species')%>%
    mutate(Species=capitalize(gsub("\\..*","",Species)),
           Scenario=Scen2)
  ddd=rbind(dummy,dummy1)
  p=ddd%>%
    ggplot(aes(x=var))+
    geom_density(aes(color=Scenario),linewidth=1.5)+
    facet_wrap(~Species,scales='free')+
    xlab(XLAB)+ylab("Density")+
    theme_PA(axs.T.siz=22,axs.t.siz=14,strx.siz=Strx.siz)+
    theme(legend.position = "top",
          legend.title=element_blank(),
          plot.title =element_text(size=17))
  print(p)
  return(ddd%>%group_by(Species,Scenario)%>%summarise(Mean=mean(var))%>%spread(Scenario,Mean))
}
fn.display.steepness=function(d,sp,XLAB,XLIM)
{
  MinH=Min.h.shark
  MinH=0.2
  dummy=lapply(d[sp],function(x) rtruncnorm(1e4,a=MinH, b=Max.h.shark,mean=x$mean,sd=max(x$sd,0.01)))
  dummy=lapply(dummy,function(x)data.frame(var=x))
  dummy=do.call(rbind,dummy)%>%
    tibble::rownames_to_column(var='Species')%>%
    mutate(Species=capitalize(gsub("\\..*","",Species)))
  p=dummy%>%
    ggplot(aes(x=var))+
    geom_density(aes(color=Species),size=1.5)+
    facet_wrap(~Species,scales='free_y')+
    xlab(XLAB)+ylab("Density")+
    theme_PA(axs.T.siz=22,axs.t.siz=14,strx.siz=16)+
    theme(legend.position = "none",
          plot.title =element_text(size=17))+
    scale_x_continuous(limits = XLIM)
  
  return(p)
}
fn.display.steepness.sensitivity=function(d,d1,d2,sp,XLAB,Strx.siz=16,Scen1,Scen2,Scen3)
{
  for(z in 1:length(d))
  {
    d[[z]]$mean=d[[z]]
    d[[z]]$sd=d1[[z]]$sd
  }
  MinH=Min.h.shark
  MinH=0.2
  dummy=lapply(d[sp],function(x) rtruncnorm(1e4,a=MinH, b=Max.h.shark,mean=x$mean,sd=max(x$sd,0.01)))
  dummy=lapply(dummy,function(x)data.frame(var=x))
  
  dummy1=lapply(d1[sp],function(x) rtruncnorm(1e4,a=MinH, b=Max.h.shark,mean=x$mean,sd=max(x$sd,0.01)))
  dummy1=lapply(dummy1,function(x)data.frame(var=x))
  
  dummy2=lapply(d2[sp],function(x) rtruncnorm(1e4,a=MinH, b=Max.h.shark,mean=x$mean,sd=max(x$sd,0.01)))
  dummy2=lapply(dummy2,function(x)data.frame(var=x))
  
  
  dummy=do.call(rbind,dummy)%>%
    tibble::rownames_to_column(var='Species')%>%
    mutate(Species=capitalize(gsub("\\..*","",Species)),
           Scenario=Scen1)
  
  dummy1=do.call(rbind,dummy1)%>%
    tibble::rownames_to_column(var='Species')%>%
    mutate(Species=capitalize(gsub("\\..*","",Species)),
           Scenario=Scen2)
  
  dummy2=do.call(rbind,dummy2)%>%
    tibble::rownames_to_column(var='Species')%>%
    mutate(Species=capitalize(gsub("\\..*","",Species)),
           Scenario=Scen3)
  
  ddd=rbind(dummy,dummy1,dummy2)
  dA=rbind(data.frame(V=unlist(lapply(d[sp],function(x) x$mean)))%>%
             rownames_to_column(var = "Species")%>%
             mutate(Scenario="steepness S2"),
           data.frame(V=unlist(lapply(d1[sp],function(x) x$mean)))%>%
             rownames_to_column(var = "Species")%>%
             mutate(Scenario="M at age"),
           data.frame(V=unlist(lapply(d2[sp],function(x) x$mean)))%>%
             rownames_to_column(var = "Species")%>%
             mutate(Scenario="M invariant"))%>%
    mutate(Species=capitalize(Species))
  
  p=ddd%>%
    ggplot(aes(x=var))+
    geom_density(aes(color=Scenario),linewidth=1.5)+
    facet_wrap(~Species,scales='free')+
    xlab(XLAB)+ylab("Density")+
    theme_PA(axs.T.siz=22,axs.t.siz=14,strx.siz=Strx.siz)+
    theme(legend.position = "top",
          legend.title=element_blank(),
          plot.title =element_text(size=17))+
    geom_vline(data=dA,aes(xintercept = V,color=Scenario))+
    geom_text(data = dA,mapping = aes(x = V, y = -Inf, label = round(V,2)),hjust=-0.1,vjust= 0,angle = 90)
  print(p)
  return(ddd%>%group_by(Species,Scenario)%>%summarise(Mean=mean(var))%>%spread(Scenario,Mean))
}


fn.compare.dist=function(MEAN,CV,XLIM)
{
  Beta.pars=get_beta(MEAN,CV)
  Gama.pars=get_gamma(MEAN,CV)
  Log.normal.pars=c(log(MEAN),CV)
  plot(density(rnorm(1e4,MEAN,CV)),xlim=XLIM,main=paste0('Mean=',MEAN,' CV=',CV))
  lines(density(rbeta(1e4,Beta.pars[1], Beta.pars[2])),col=2)
  lines(density(rgamma(1e4,Gama.pars[1], Gama.pars[2])),col=3)
  lines(density(rlnorm(1e4,Log.normal.pars[1], Log.normal.pars[2])),col=4)
  legend('topleft',c('normal','beta','gamma','log.normal'),lty=1,col=1:4,bty='n')
}
#fn.compare.dist(MEAN=0.317,CV=0.5,XLIM=c(0,1))
make_plot <- function(da,nfacets,AXST,AXSt,STRs,InrMarg,dropTitl,addKtch,YLAB='',HLine)
{
  p=da%>%
    ggplot(aes(year, median))+
    geom_ribbon(aes(ymin = lower.95, ymax = upper.95), alpha = 0.2,fill='grey60')+
    geom_ribbon(aes(ymin = lower.60, ymax = upper.60), alpha = 0.2,fill='grey20')+
    geom_line()  +
    facet_wrap(~Scenario,ncol=nfacets)+
    theme_PA(axs.T.siz=AXST,axs.t.siz=AXSt,strx.siz=STRs)+
    theme(plot.title =element_text(hjust = 0.5))+
    ylab(YLAB)+xlab(XLAB)+   
    #ylim(0,max(da$upper.95))+
    theme(panel.spacing=unit(InrMarg,"lines"),
          plot.margin = unit(c(.5, -.2, 0, 0), "cm"))
  if(any(!is.na(da$upper.50))) p=p+geom_ribbon(aes(ymin = lower.50, ymax = upper.50), alpha = 0.1)
  if(!is.null(HLine))
  {
    if('Target'%in%names(da)) p=p+geom_hline(data=da,aes(yintercept=unique(Target)),color='forestgreen', size=1.05)
    p=p+geom_hline(data=da,aes(yintercept=unique(Threshold)),color='orange', size=1.05)
    if('Limit'%in%names(da)) p=p+geom_hline(data=da,aes(yintercept=unique(Limit)),color='red', size=1.05)
  }
  if(!addKtch) p=p+scale_y_continuous(limits = c(0, NA))
  if(addKtch)
  {
    sklr=1
    if(unique(da$Threshold)==1) sklr=2
    coeff=max(da$Catch)
    da$Ktch.scaled=da$Catch/(coeff/sklr)
    lbs=seq(0,10*ceiling(coeff/10),length.out=4)
    lbs=pretty(lbs,n=(length(lbs)-1))
    bks=seq(0,10*ceiling(sklr*coeff/10),length.out=length(lbs))
    bks=pretty(bks,n=(length(bks)-1))
    min.lb=min(length(lbs),length(bks))
    lbs=lbs[1:min.lb]
    bks=bks[1:min.lb]
    p=p+ 
      geom_line(data=da, aes(x=year, y=Ktch.scaled),size=1.1,color='dodgerblue4',alpha=0.5,linetype ='dashed')+
      scale_y_continuous(limits = c(0, NA),sec.axis = sec_axis(~.*coeff, name="",labels = lbs,breaks = bks))
  }
  p=p+geom_line(size=.7)
  future.ktch=da%>%filter(year>Last.yr.ktch.numeric)
  if(nrow(future.ktch)>0)
  {
    p=p+
      geom_ribbon(data=future.ktch,aes(ymin = lower.95, ymax = upper.95), alpha = 0.2,fill=future.color) +
      geom_line(data=future.ktch,aes(year,median),color=future.color,size=1.1)
    if(addKtch)p=p+geom_line(data=future.ktch, aes(x=year, y=Ktch.scaled),size=1.1,color=future.color)
    
  }
  
  return(p)
}
fn.ribbon=function(Dat,YLAB,XLAB,Titl,HLIne,addKtch,nfacets=1,AXST=14,AXSt=12,STRs=14,InrMarg=.25,dropTitl=FALSE)
{
  data2 <- split(Dat, Dat$Scenario)
  p_lst <- lapply(data2, make_plot,nfacets,AXST,AXSt,STRs,InrMarg,dropTitl,addKtch,HLine=HLIne)
  for(ii in 1:length(p_lst))
  {
    p_lst[[ii]]=p_lst[[ii]]+theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
  }
  figure <- ggarrange(plotlist=p_lst,ncol=nfacets,nrow=ceiling(length(p_lst)/nfacets),
                      common.legend = FALSE)
  if(!dropTitl) figure <- annotate_figure(figure,top=text_grob(Titl, size=20))
  
  return(figure) 
}
fn.plot.timeseries=function(d,sp,Type,YLAB,add.50=FALSE,add.sp.nm=FALSE)
{
  mods=names(d)
  store.plots=vector('list',length(mods))
  names(store.plots)=mods
  store.probs=str.ref.point=store.plots
  id=match(sp,Keep.species)
  for(m in 1:length(mods))
  {
    if(!is.null(d[[m]]$estimates[[id]]))
    {
      Hline=NULL 
      addKtch=FALSE
      if(Type=='F.series')
      {
        Hline=TRUE
        Var=d[[m]]$f.series[[id]]
        didi=d[[m]]$estimates[[id]]%>%
          dplyr::select(Parameter,Scenario,Median)%>%
          filter(Parameter%in%c('F.target','F.threshold','F.limit'))
        dis=data.frame(Target=didi%>%filter(Parameter=='F.target')%>%pull(Median),
                       Threshold=didi%>%filter(Parameter=='F.threshold')%>%pull(Median),
                       Limit=didi%>%filter(Parameter=='F.limit')%>%pull(Median),
                       Scenario=unique(didi$Scenario))
        Var=Var%>%left_join(dis,by='Scenario')
      }
      if(Type=='B.Bmsy')
      {
        Hline=TRUE  
        Var=d[[m]]$B.Bmsy[[id]]
        dis=data.frame(Limit=Lim.prop.bmsy,Target=Tar.prop.bmsny,Threshold=1,Scenario=unique(Var$Scenario))
        Var=Var%>%left_join(dis,by='Scenario')
        YLAB=expression("B/B"[MSY])
      }
      if(Type=='F.Fmsy')
      {
        Hline=TRUE  
        Var=d[[m]]$F.Fmsy[[id]]
        dis=data.frame(Target=1-(Tar.prop.bmsny-1),
                       Threshold=1,
                       Limit=1+Lim.prop.bmsy,
                       Scenario=unique(Var$Scenario))
        Var=Var%>%left_join(dis,by='Scenario')
        YLAB=expression("F/F"[MSY])
      }
      if(Type=='Depletion')
      {
        addKtch=TRUE
        Var=d[[m]]$rel.biom[[id]]
        str.prob=compact(d[[m]]$probs.rel.biom[[id]])   
        #dis=str.prob[[1]]$Reference.points   
        dis=fn.get.stuff.from.list(str.prob,'Reference.points')
        for(y in 1:length(dis)) dis[[y]]$Scenario=paste('S',y,sep='')
        Var=Var%>%left_join(do.call(rbind,dis)%>%spread(Rf.pt,Value),by='Scenario')
        Hline='YES'
        str.ref.point[[m]]=dis
        #Hline=str.ref.point[[m]]$Value
        store.probs[[m]]=do.call(rbind,sapply(str.prob,'[',1))%>%  
          mutate(Model=names(store.probs)[m])
      }
      Nfacets=1
      if(length(unique(Var$Scenario))>6) Nfacets=2
      store.plots[[m]]=fn.ribbon(Dat=Var,  
                                 YLAB='',
                                 XLAB="",
                                 Titl=names(d)[m],
                                 HLIne=Hline,
                                 addKtch=addKtch,
                                 nfacets=Nfacets)
      
    }
  }
  store.plots=compact(store.plots)
  if(length(store.plots)>0)
  {
    figure=ggarrange(plotlist=store.plots, nrow = length(mods),common.legend=TRUE)
    if(add.sp.nm) figure=figure+theme(plot.margin = margin(1,0,0,0, "cm"))

    figure=annotate_figure(figure,
                           bottom = text_grob('Financial year',size=26,vjust =-0.15),
                           left = text_grob(YLAB,size=26,rot = 90,vjust=0.8))
    if(add.sp.nm) figure=annotate_figure(figure,
                                         fig.lab=capitalize(sp),
                                         fig.lab.pos='top.left',
                                         fig.lab.size=28)
    if(addKtch)
    {
      figure=annotate_figure(figure,
                             right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                             color ='dodgerblue4',vjust=0))
    }
    
    print(figure)      
    store.probs=do.call(rbind,store.probs)
    rownames(store.probs)=NULL
    return(list(store.probs=store.probs,Ref.points=str.ref.point))
    
  }
}
fn.plot.timeseries_combined=function(this.sp,d,YLAB,Type,InnerMargin,RefPoint,Kach,SCEn='S1')
{
  id=match(this.sp,names(d$rel.biom))
  RefPoint=compact(RefPoint[id])
  Kach=compact(Kach[id])
  
  if(length(Kach)>0)
  {
    Kach=do.call(rbind,Kach)%>%
      filter(Scenario=='S1')%>%
      rownames_to_column()%>%
      mutate(Scenario=capitalize(word(rowname,1,sep = "\\.")))%>%
      dplyr::select(year,Catch,Scenario)
    
    Hline=NULL
    addKtch=FALSE
    if(Type=='F.series')
    {
      Hline='YES'
      Var=compact(d$f.series[id])
      Est=compact(d$estimates[id])
      for(q in 1:length(Var))
      {
        didi=Est[[q]]%>%
          dplyr::select(Parameter,Scenario,Median)%>%
          filter(Scenario=='S1' & Parameter%in%c('F.target','F.threshold','F.limit'))
        Var[[q]]=Var[[q]]%>%
          mutate(Target=didi%>%filter(Parameter=='F.target')%>%pull(Median),
                 Threshold=didi%>%filter(Parameter=='F.threshold')%>%pull(Median),
                 Limit=didi%>%filter(Parameter=='F.limit')%>%pull(Median))
      }
    }
    if(Type=='B.Bmsy')
    {
      Hline='YES'
      Var=d$B.Bmsy[id]
      Var=compact(Var)
      for(q in 1:length(Var))
      {
        Var[[q]]=Var[[q]]%>%
          mutate(Target=Tar.prop.bmsny,
                 Threshold=1,
                 Limit=Lim.prop.bmsy)
      }
      YLAB=expression("B/B"[MSY])
    }
    if(Type=='F.Fmsy')
    {
      Hline='YES'
      Var=d$F.Fmsy[id]
      Var=compact(Var)
      for(q in 1:length(Var))
      {
        Var[[q]]=Var[[q]]%>%
          mutate(Target=1-(Tar.prop.bmsny-1),
                 Threshold=1,
                 Limit=1+Lim.prop.bmsy)
      }
      YLAB=expression("F/F"[MSY])
    }
    if(Type=='Depletion')
    {
      Var=compact(d$rel.biom[id])
      addKtch=TRUE
      for(q in 1:length(Var))
      {
        RP=RefPoint[[q]]%>%data.frame
        if(!is.data.frame(RP))
        {
          if(length(RP)>1)
          {
            RP=do.call(rbind,RP)%>%
              filter(Scenario==SCEn)
          }
        }
        Var[[q]]=Var[[q]]%>%
          mutate(Target=RP%>%filter(Rf.pt=='Target')%>%pull( Value),
                 Threshold=RP%>%filter(Rf.pt=='Threshold')%>%pull( Value),
                 Limit=RP%>%filter(Rf.pt=='Limit')%>%pull( Value))
      }
      Hline='YES'
    }
    Var=compact(Var)
    
    if(length(Var)>=8)  Nfast=round(length(Var)/5)
    if(length(Var)>3 & length(Var)<8) Nfast=2
    if(length(Var)<=3) Nfast=1
    
    Dat=do.call(rbind,Var)
    if('Catch'%in%names(Dat))Dat=Dat%>%dplyr::select(-Catch)
    if('Scenario'%in%names(Dat))Dat=Dat%>%filter(Scenario==SCEn)
    STRs=16
    aXSt=13
    if(length(Var)>8)
    {
      STRs=13
      aXSt=10
    }
      
    figure=fn.ribbon(Dat=Dat%>%
                       rownames_to_column()%>%
                       mutate(Scenario=capitalize(word(rowname,1,sep = "\\.")))%>%
                       left_join(Kach,by=c('year','Scenario')),
                     YLAB='',
                     XLAB="",
                     Titl='',
                     HLIne=Hline,
                     addKtch=addKtch,
                     nfacets=Nfast,
                     AXST=15,AXSt=aXSt,STRs=STRs,
                     InrMarg=InnerMargin,
                     dropTitl=TRUE)
    figure=annotate_figure(figure,
                           bottom = text_grob('Financial year',size=26,vjust =-0.15),
                           left = text_grob(YLAB,size=26,rot = 90,vjust=0.8)) +
      theme(plot.margin = unit(c(.1,.1,0,0), "cm"))
    if(addKtch)
    {
      figure=annotate_figure(figure,
                             right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                             color ='dodgerblue4',vjust=0))
    }
    return(figure)
  }
}
fn.plot.timeseries_combined_Appendix=function(this.sp,d,YLAB,NM,Type,InnerMargin)
{
  mods=names(d)
  Out.name='Relative.biomass'
  if(Type=="F.series") Out.name='F.series'
  if(Type=="B.Bmsy") Out.name='B.over.Bmsy'
  if(Type=="F.Fmsy") Out.name='F.over.Fmsy'
  if(length(this.sp)>8)
  {
    for(m in 1:length(mods))
    {
      id=match(this.sp,names(d[[m]]$rel.biom))
      Hline=NULL
      addKtch=FALSE
      if(Type=='F.series')  Var=d[[m]]$f.series[id]
      if(Type=='B.Bmsy')
      {
        Hline=TRUE  
        addKtch=FALSE
        Var=d[[m]]$B.Bmsy[id]
        dis=data.frame(Limit=Lim.prop.bmsy,Target=Tar.prop.bmsny,Threshold=1,Scenario=unique(Var[[1]]$Scenario))
        for(q in 1:length(Var)) Var[[q]]=Var[[q]]%>%left_join(dis,by='Scenario')
        YLAB=expression("B/B"[MSY])
      }
      if(Type=='F.Fmsy')
      {
        Hline=TRUE  
        addKtch=FALSE
        Var=d[[m]]$F.Fmsy[id]
        dis=data.frame(Threshold=1,Scenario=unique(Var[[1]]$Scenario))
        for(q in 1:length(Var)) Var[[q]]=Var[[q]]%>%left_join(dis,by='Scenario')
        YLAB=expression("F/F"[MSY])
      }
      if(Type=='Depletion')
      {
        addKtch=TRUE
        Var=d[[m]]$rel.biom[id]
        str.prob=compact(d[[m]]$probs.rel.biom[id]) 
        for(q in 1:length(str.prob))
        {
          dis=fn.get.stuff.from.list(str.prob[[q]],'Reference.points')
          for(y in 1:length(dis)) dis[[y]]$Scenario=paste('S',y,sep='')
          Var[[q]]=Var[[q]]%>%left_join(do.call(rbind,dis)%>%spread(Rf.pt,Value),by='Scenario')
        }
        Hline='YES' #Hline=d[[m]]$probs.rel.biom[[id[1]]][[1]]$Reference.points$Value 
      }
      figure=fn.ribbon(Dat=do.call(rbind,Var)%>%
                         rownames_to_column()%>%
                         filter(Scenario=='S1')%>%
                         mutate(Scenario=capitalize(word(rowname,1,sep = "\\."))),
                       YLAB='',
                       XLAB="",
                       Titl=mods[m],
                       HLIne=Hline,
                       addKtch=addKtch,
                       nfacets=round(length(this.sp)/5),
                       AXST=15,AXSt=12,STRs=12,
                       InrMarg=InnerMargin,
                       dropTitl=TRUE)

      figure=annotate_figure(figure,
                             bottom = text_grob('Financial year',size=26,vjust =-0.15),
                             left = text_grob(YLAB,size=26,rot = 90,vjust=0.8))+
        theme(plot.margin = unit(c(.1,.1,0,0), "cm"))
      if(addKtch)
      {
        figure=annotate_figure(figure,
                               right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                               color ='dodgerblue4',vjust=-.2))
      }
      print(figure)
      ggsave(paste(Rar.path,'/',Out.name,'_catch.only_',NM,'_',mods[m],'_Appendix_S1','.tiff',sep=''),
             width = 14,height = 11,compression = "lzw")
    }
    
  }else
  {
    plot.list=vector('list',length(mods))
    for(m in 1:length(mods))
    {
      id=match(this.sp,names(d[[m]]$rel.biom))
      Hline=NULL
      addKtch=FALSE
      if(Type=='F.series')  Var=d[[m]]$f.series[id] 
      if(Type=='B.Bmsy')
      {
        Hline=TRUE  
        addKtch=FALSE
        Var=d[[m]]$B.Bmsy[id]
        dis=data.frame(Limit=Lim.prop.bmsy,Target=Tar.prop.bmsny,Threshold=1,Scenario=unique(Var[[1]]$Scenario))
        for(q in 1:length(Var)) Var[[q]]=Var[[q]]%>%left_join(dis,by='Scenario')
      }
      if(Type=='F.Fmsy')
      {
        Hline=TRUE  
        addKtch=FALSE
        Var=d[[m]]$F.Fmsy[id]
        dis=data.frame(Threshold=1,Scenario=unique(Var[[1]]$Scenario))
        for(q in 1:length(Var)) Var[[q]]=Var[[q]]%>%left_join(dis,by='Scenario')
      }
      if(Type=='Depletion')
      {
        addKtch=TRUE
        Var=d[[m]]$rel.biom[id]
        str.prob=compact(d[[m]]$probs.rel.biom[id]) 
        for(q in 1:length(str.prob))
        {
          dis=fn.get.stuff.from.list(str.prob[[q]],'Reference.points')
          for(y in 1:length(dis)) dis[[y]]$Scenario=paste('S',y,sep='')
          Var[[q]]=Var[[q]]%>%left_join(do.call(rbind,dis)%>%spread(Rf.pt,Value),by='Scenario')
        }
        Hline='YES'
        #Hline=d[[m]]$probs.rel.biom[[id[1]]][[1]]$Reference.points$Value 
      }
      plot.list[[m]]=fn.ribbon(Dat=do.call(rbind,Var)%>%
                                 filter(Scenario=='S1')%>%
                                 rownames_to_column()%>%
                                 mutate(Scenario=capitalize(word(rowname,1,sep = "\\."))),
                               YLAB='',
                               XLAB="",
                               Titl=mods[m],
                               HLIne=Hline,
                               addKtch=addKtch,
                               nfacets=1,
                               AXST=16,AXSt=14,STRs=16,
                               InrMarg=InnerMargin)
    }
    
    figure=ggarrange(plotlist=plot.list,ncol=length(mods), common.legend=TRUE)+
      theme(plot.margin = margin(0,0,0,0, "cm"))
    figure=annotate_figure(figure,
                           bottom = text_grob('Financial year',size=26,vjust =-0.15),
                           left = text_grob(YLAB,size=26,rot = 90,vjust=0.8))
    if(addKtch)
    {
      figure=annotate_figure(figure,
                             right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                             color ='dodgerblue4',vjust=-.2))
    }
    print(figure)
    ggsave(paste(Rar.path,'/',Out.name,'_catch.only_',NM,'_Appendix_S1','.tiff',sep=''),
           width = 14,height = 10,compression = "lzw")
    
  }
}
fn.plot.timeseries_combined_sensitivity=function(this.sp,d,InnerMargin,RefPoint,Kach,SCEn='S1')
{
  id=match(this.sp,names(d$rel.biom))
  RefPoint=compact(RefPoint[id])
  Kach=compact(Kach[id])
  if(length(Kach)>0)
  {
    ach=do.call(rbind,Kach)%>%
      rownames_to_column()%>%
      mutate(Species=capitalize(sub("\\..*", "", rowname)))
    if(length(RefPoint[[1]])>1)
    {
      for(ref in 1:length(RefPoint)) RefPoint[[ref]]=do.call(rbind,RefPoint[[ref]])
    }
    Hline=do.call(rbind,RefPoint)%>%  
      rownames_to_column()%>%
      mutate(Species=capitalize(sub("\\..*", "", rowname)))%>%
      dplyr::select(-rowname)%>%
      filter(Scenario==SCEn)%>%
      spread(Rf.pt,Value)%>%
      dplyr::select(-Scenario)
    ach=ach%>%left_join(Hline,by='Species')
    figure=ach%>%
      ggplot(aes(year,median,color=Scenario))+
      facet_wrap(~Species)+
      geom_ribbon(aes(ymin = lower.95, ymax = upper.95,fill=Scenario), alpha = 0.1)+
      geom_line(size=2)+
      guides(color = guide_legend(nrow = 1),fill= guide_legend(nrow = 1))+
      ylim(0,max(ach$upper.95))+
      theme_PA(strx.siz=14)+
      theme(legend.position = 'top',
            legend.title = element_blank())+
      ylab('')+xlab('')+
      geom_hline(aes(yintercept=Limit), size=1.05,alpha=0.35,color='red')+
      geom_hline(aes(yintercept=Threshold), size=1.05,alpha=0.35,color='orange')+
      geom_hline(aes(yintercept=Target), size=1.05,alpha=0.35,color='forestgreen')
    figure=annotate_figure(figure,
                           bottom = text_grob('Financial year',size=26,vjust =-0.15),
                           left = text_grob("Relative biomass",size=26,rot = 90,vjust=0.8)) +
      theme(plot.margin = unit(c(.1,.1,0,0), "cm"))
    return(figure)
  }
}
fn.plot.catch.vs.MSY=function(d,sp,YLAB)
{
  mods=names(d)
  store.plots=vector('list',length(mods))
  names(store.plots)=mods
  id=match(sp,Keep.species)
  for(m in 1:length(mods))
  {
    if(!is.null(d[[m]]$estimates[[id]]))
    {
      MSY.estim=d[[m]]$estimates[[id]]
      if(!'Parameter'%in%names(MSY.estim)) MSY.estim$Parameter=MSY.estim$Par
      MSY.estimS2=MSY.estim[grep('MSY',MSY.estim$Parameter),]%>%filter(Scenario=='S2')
      MSY.estim=MSY.estim[grep('MSY',MSY.estim$Parameter),]%>%filter(Scenario=='S1')
      
      MSY.estim.60.per=d[[m]]$ensemble[[id]]$MSY%>%
                        group_by(Scenario)%>%
                        summarise(Upper.60=quantile(Value,probs=0.8,na.rm=T),
                                  Lower.60=quantile(Value,probs=0.2,na.rm=T))   
      MSY.estimS2.60.per=MSY.estim.60.per%>%filter(Scenario=='S2')
      MSY.estim.60.per=MSY.estim.60.per%>%filter(Scenario=='S1')
      if(nrow(MSY.estimS2)==0)
      {
        MSY.estimS2=MSY.estim
        MSY.estimS2.60.per=MSY.estim.60.per
      }
        
      Ktch=d[[m]]$rel.biom[[id]]%>%filter(Scenario=='S1')
      
      suppressWarnings({store.plots[[m]]=data.frame(Year=Ktch$year,
                                                    Catch=Ktch$Catch,
                                                    MSY_Lower.95=MSY.estim$Lower.95,
                                                    MSY_Median=MSY.estim$Median,
                                                    MSY_Upper.95=MSY.estim$Upper.95,
                                                    MSY_Lower.95_S2=MSY.estimS2$Lower.95,
                                                    MSY_Median_S2=MSY.estimS2$Median,
                                                    MSY_Upper.95_S2=MSY.estimS2$Upper.95,
                                                    MSY_Upper.60=MSY.estim.60.per$Upper.60,
                                                    MSY_Lower.60=MSY.estim.60.per$Lower.60,
                                                    MSY_Upper.60_S2=MSY.estimS2.60.per$Upper.60,
                                                    MSY_Lower.60_S2=MSY.estimS2.60.per$Lower.60,
                                                    Model=mods[m])}) 
    }
  }
  store.plots=compact(store.plots)
  if(length(store.plots)>0)
  {
    CoLS1="black"
    ColS2="darkcyan"
    if(!"S2"%in%unique(d$SSS$estimates[[id]]$Scenario))ColS2="transparent"  
    store.plots=do.call(rbind,store.plots)
    
    figure=store.plots%>%
      ggplot(aes(Year,Catch))+
      geom_line(aes(Year,MSY_Median),color=CoLS1,size=.65)+
      geom_ribbon(aes(ymin = MSY_Lower.95, ymax = MSY_Upper.95), alpha=.1,fill =CoLS1)+
      geom_ribbon(aes(ymin = MSY_Lower.60, ymax = MSY_Upper.60), alpha=.3,fill =CoLS1)+
      geom_line(aes(Year,MSY_Median_S2),color=ColS2,size=.65)+
      geom_ribbon(aes(ymin = MSY_Lower.95_S2, ymax = MSY_Upper.95_S2), alpha=.1,fill =ColS2)+
      geom_ribbon(aes(ymin = MSY_Lower.60_S2, ymax = MSY_Upper.60_S2), alpha=.3,fill =ColS2)+
      geom_line(size=.75,color='dodgerblue4')+
      facet_wrap(~Model,ncol=1)+
      theme_PA()+ylab(YLAB)+xlab('Financial year')
    future.ktch=store.plots%>%filter(Year>Last.yr.ktch.numeric)
    if(nrow(future.ktch)>0)
    {
      figure=figure+
              geom_line(data=future.ktch,aes(Year,Catch),size=.75,color=future.color)
    }
    suppressWarnings({print(figure)})     
    return(store.plots%>%mutate(Species=sp))
  }
}
fn.plot.catch.vs.MSY_combined_Appendix=function(all.this.sp,this.sp,d,NM,YLAB)
{
  #table
  d1=d%>%
    filter(Species%in%all.this.sp)%>%
    filter(Year<=Last.yr.ktch.numeric)%>%
    mutate(Species=capitalize(Species))%>%
    group_by(Model,Species)%>%
    summarise(Catch.last.yr=mean(data.table::last(Catch,1)),
              Catch.mean.last.3.yrs=mean(data.table::last(Catch,3)),
              MSY_Lower.95=mean(MSY_Lower.95),
              MSY_Lower.60=mean(MSY_Lower.60),
              MSY_Median=mean(MSY_Median),
              MSY_Upper.60=mean(MSY_Upper.60),
              MSY_Upper.95=mean(MSY_Upper.95))
  write.csv(d1,
            paste(Rar.path,'/Table 4. Catch.only_catch_vs_MSY_',NM,'_Appendix','.csv',sep=''),
            row.names=F)
  
  #figure
  if(length(this.sp)>0)
  {
    d=d%>%
      filter(Species%in%this.sp)%>%
      mutate(Species=capitalize(Species))
    mods=unique(d$Model)
    nnn=length(unique(d$Species))
    Nkl=ifelse(nnn<4,1,ifelse(nnn>=4 & nnn<=8,2,3))
    WID=ifelse(Nkl<3,7,7)
    aXS.si=ifelse(Nkl<3,12,9)
    STR.siz=ifelse(Nkl<3,15,10.5)

    for(m in 1:length(mods))
    {
      figure=d%>%
        filter(Model==mods[m])%>%
        ggplot(aes(Year,Catch))+
        geom_line(aes(Year,MSY_Median),color="black",size=.65)+
        geom_ribbon(aes(ymin = MSY_Lower.95, ymax = MSY_Upper.95), alpha=.1,fill ="black")+
        geom_ribbon(aes(ymin = MSY_Lower.60, ymax = MSY_Upper.60), alpha=.3,fill ="black")+
        geom_line(size=.75,color='dodgerblue4')+
        facet_wrap(~Species,scales='free_y',ncol=Nkl)+
        theme_PA(axs.t.siz=aXS.si,strx.siz=STR.siz)+
        ylab(YLAB)+xlab("Financial year")
      future.ktch=d%>%filter(Year>Last.yr.ktch.numeric)
      if(nrow(future.ktch)>0)
      {
        figure=figure+
          geom_line(data=future.ktch,aes(Year,Catch),size=.75,color=future.color)
      }
      print(figure)
      ggsave(paste(Rar.path,'/Catch_vs_MSY_catch.only_',NM,'_',mods[m],'_Appendix_S1','.tiff',sep=''),
             width = WID,height = 7,compression = "lzw")
    }
  }

}
fn.plot.catch.vs.MSY_combined=function(all.this.sp,this.sp,d,NM,YLAB,MSY)
{
  for(k in 1:length(MSY)) MSY[[k]]$Species=names(MSY)[k]
  MSY=do.call(rbind,MSY)%>%
    rename(ensembled.median=median,
           ensembled.lower.60=lower.60,
           ensembled.upper.60=upper.60,
           ensembled.lower.95=lower.95,
           ensembled.upper.95=upper.95)
  
  #table
  d1=d%>%
    #dplyr::select(-c(MSY_Lower.95,MSY_Median,MSY_Upper.95))%>%
    dplyr::select(Year,Catch, Species)%>%
    filter(Species%in%all.this.sp)%>%
    filter(Year<=Last.yr.ktch.numeric)%>%
    left_join(MSY,by='Species')%>%
    mutate(Species=capitalize(Species))%>%
    group_by(Species)%>%
    summarise(Catch.last.yr=mean(data.table::last(Catch,1)),
              MSY_Lower.95=mean(ensembled.lower.95),
              MSY_Lower.60=mean(ensembled.lower.60),
              MSY_Median=mean(ensembled.median),
              MSY_Upper.60=mean(ensembled.upper.60),
              MSY_Upper.95=mean(ensembled.upper.95))
  d1=d1%>%left_join(MSY%>%dplyr::select(-c(ensembled.median,ensembled.lower.60,
                                           ensembled.upper.60,ensembled.lower.95,ensembled.upper.95))%>%
                      mutate(Species=capitalize(Species)),by='Species')
  write.csv(d1,
            paste(Rar.path,'/Table 4. Catch.only_catch_vs_MSY_',NM,'.csv',sep=''),
            row.names=F)
  
  #figure
  if(length(this.sp)>0)
  {
    MSY=MSY%>%filter(Species%in%this.sp)
    d=d%>%
      filter(Species%in%this.sp)
    mods=unique(d$Model)
    nnn=length(unique(d$Species))
    Nkl=ifelse(nnn<4,1,ifelse(nnn>=4 & nnn<=8,2,4))
    WID=ifelse(Nkl<3,7,9)
    aXS.si=ifelse(Nkl<3,12,9)
    STR.siz=ifelse(Nkl<3,15,10.5)
    d=d%>%
      filter(Model==mods[1])%>%
      left_join(MSY,by='Species')%>%
      mutate(Species=capitalize(Species))
    
    figure=d%>%
      ggplot(aes(Year,Catch))+
      geom_line(aes(Year,ensembled.median),color="black",size=.65)+
      geom_ribbon(aes(ymin = ensembled.lower.95, ymax = ensembled.upper.95), alpha=.1,fill ="black")+
      geom_ribbon(aes(ymin = ensembled.lower.60, ymax = ensembled.upper.60), alpha=.3,fill ="black")+
      geom_line(size=.75,color='dodgerblue4')+
      facet_wrap(~Species,scales='free',ncol=Nkl)+
      theme_PA(axs.t.siz=aXS.si,strx.siz=STR.siz)+
      ylab(YLAB)+xlab("Financial year")
    future.ktch=d%>%filter(Year>Last.yr.ktch.numeric)
    if(nrow(future.ktch)>0)
    {
      figure=figure+
        geom_line(data=future.ktch,aes(Year,Catch),size=.75,color=future.color)
    }
    print(figure)
    ggsave(paste(Rar.path,'/Catch_vs_MSY_catch.only_',NM,'.tiff',sep=''),
           width = WID,height = 7,compression = "lzw")
    
  }
}
kobePlot <- function(f.traj,b.traj,Years,Titl,Probs=NULL,txt.col='black',YrSize=4)
{
  dta=data.frame(x=b.traj,
                 y=f.traj,
                 yr=Years)%>%
    arrange(yr)
  Mx.F=max(2,max(dta$y,na.rm=T))
  Mx.B=max(2,max(dta$x,na.rm=T))
  
  kobe <-dta%>%
    ggplot(aes(x, y))+    
    scale_x_continuous(limits=c(0,Mx.B),expand = c(0, 0)) +
    scale_y_continuous(limits=c(0,Mx.F),expand = c(0, 0))+
    geom_rect(xmin = 1, xmax = Mx.B, ymin = 0, ymax = 1, fill = RiskColors['Low'], alpha = 0.05) +
    geom_rect(xmin = 0, xmax = 1, ymin = 1, ymax = Mx.F, fill = RiskColors['Severe'], alpha = 0.05) +
    geom_rect(xmin = 1, xmax = Mx.B, ymin = 1, ymax = Mx.F, fill = RiskColors['High'], alpha = 0.05) +
    geom_rect(xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = RiskColors['Medium'], alpha = 0.05)
  if(!is.null(Probs))
  {
    kernelF <- gplots::ci2d(Probs$x, Probs$y, nbins = 151, factor = 1.5, 
                            ci.levels = c(0.5, 0.8, 0.75, 0.9, 0.95), show = "none")
    KernelD=rbind(kernelF$contours$"0.95"%>%mutate(CI='1',col='grey30'),
                  kernelF$contours$"0.8"%>%mutate(CI='2',col='grey60'),
                  kernelF$contours$"0.5"%>%mutate(CI='3',col='grey85'))
    kernels=KernelD%>%distinct(CI,col)%>%pull(col)
    names(kernels)=KernelD%>%distinct(CI,col)%>%pull(CI)
    
    Pr.d=data.frame(
                Prob=c(sum(ifelse(Probs$x >= 1 & Probs$y <= 1, 1, 0))/length(Probs$x)*100,
                       sum(ifelse(Probs$x < 1 & Probs$y <= 1, 1, 0))/length(Probs$x)*100,
                       sum(ifelse(Probs$x >= 1 & Probs$y > 1, 1, 0))/length(Probs$x)*100,
                       sum(ifelse(Probs$x < 1 & Probs$y > 1, 1, 0))/length(Probs$x) * 100),
                col=RiskColors[c('Low','Medium','High','Severe')],
                x=rep(-10,4),  #dummy
                y=rep(-10,4))
    pr.ds=Pr.d%>%pull(col)
    Nms.probs=round(Pr.d%>%pull(Prob),1)
    Nms.probs=ifelse(Nms.probs>0 & Nms.probs<0.1,'<0.1',Nms.probs)
    names(pr.ds)=paste(Nms.probs,'%',sep='')
    #pr.ds=pr.ds[which(!names(pr.ds)=='0%')]
    kobe <-kobe +
      geom_polygon(data=KernelD%>%dplyr::filter(y<=Mx.F & x<=Mx.B),aes(x, y,fill=CI),size=1.25,alpha=0.5)+
      scale_fill_manual(labels=c("95%","80%","50%"),values = kernels)+
      geom_point(data=Pr.d,aes(x, y,color=col),alpha = 1,size=5)+
      scale_color_manual(labels=names(pr.ds),values =Pr.d$col)+
      labs(CI="", col=dta[nrow(dta),'yr'])
    
  }
  kobe <-kobe +
    geom_path(linetype = 2, size = 0.5,color='deepskyblue4')+
    geom_point(size=2,color='deepskyblue4')+
    geom_point(aes(x=dta[1,'x'],y=dta[1,'y']),size=4,shape=22,fill='white',alpha=.3)+
    geom_point(aes(x=dta[nrow(dta),'x'],y=dta[nrow(dta),'y']),size=4,shape=25,fill='white',alpha=.3)+      
    geom_text_repel(data=dta[1,],aes(x=x,y=y,label=yr),size=YrSize,color=txt.col)+
    geom_text_repel(data=dta[nrow(dta),],aes(x=x,y=y,label=yr),size=YrSize,color=txt.col)+
    xlab(expression(B/~B[MSY]))+ylab(expression(F/~F[MSY]))+
    labs(title = Titl)+
    theme_bw()%+replace% 
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size=16),
          axis.title = element_text(size=20),
          plot.title = element_text(size=20,hjust=0),
          legend.text = element_text(size=15),
          legend.title = element_text(size=17),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-10))
  if(!is.null(Probs)) kobe=kobe+theme(legend.title = element_blank(),legend.text.align = 1)
  return(kobe)
}
fn.get.Kobe.plot_appendix=function(d,sp,Scen='S1',add.sp.nm=FALSE,do.probs=FALSE,txt.size=6)
{
  id=match(sp,Keep.species)
  
  plotlist=vector('list',length(d))
  names(plotlist)=names(d)
  for(pp in 1:length(d))
  {
    #DBSRA
    if('DBSRA'%in%names(d)[[pp]])
    {
      dummy=d$DBSRA$B.Bmsy[[id]]%>%
                  filter(Scenario==Scen & year<=Last.yr.ktch.numeric)
      yrs=dummy%>%pull(year)
      Bmsy=dummy$median
      Fmsy=d$DBSRA$F.Fmsy[[id]]%>%
              filter(Scenario==Scen& year<=Last.yr.ktch.numeric)%>%
              pull(median)
      p.plot=kobePlot(f.traj=Fmsy[1:length(yrs)],
                       b.traj=Bmsy[1:length(yrs)],
                       Years=yrs,
                       Titl="DBSRA",
                       YrSize=txt.size)
      rm(yrs,Fmsy,Bmsy,dummy)
    }
    
    #CMSY
    if('CMSY'%in%names(d)[[pp]])
    {
      if(CMSY.method=="Froese")
      {
        yrs=d$CMSY[[id]][[Scen]]$output$ref_ts$year
        Bmsy=d$CMSY[[id]][[Scen]]$output$ref_ts$bbmsy
        Fmsy=d$CMSY[[id]][[Scen]]$output$ref_ts$ffmsy
        
      }
      if(CMSY.method=="Haddon")
      {
        dummy=d$CMSY$B.Bmsy[[id]]%>%
                  filter(Scenario==Scen & year<=Last.yr.ktch.numeric)
        yrs=dummy%>%pull(year)
        Bmsy=dummy$median
        Fmsy=d$CMSY$F.Fmsy[[id]]%>%
              filter(Scenario==Scen & year<=Last.yr.ktch.numeric)%>%
              pull(median)
      }
      p.plot=kobePlot(f.traj=Fmsy,
                      b.traj=Bmsy,
                      Years=yrs,
                      Titl="CMSY",
                      YrSize=txt.size)
      rm(yrs,Fmsy,Bmsy,dummy)
    }
    
    #JABBA
    if('JABBA'%in%names(d)[[pp]])
    {
      dummy=d$JABBA$B.Bmsy[[id]]%>%
                filter(Scenario==Scen & year<=Last.yr.ktch.numeric)
      yrs=dummy%>%pull(year)
      Bmsy=dummy$median
      Fmsy=d$JABBA$F.Fmsy[[id]]%>%
            filter(Scenario==Scen & year<=Last.yr.ktch.numeric)%>%
            pull(median)
      if(!do.probs)
      {
        p.plot=kobePlot(f.traj=Fmsy,
                         b.traj=Bmsy,
                         Years=yrs,
                         Titl="JABBA",
                        YrSize=txt.size)
      }
      if(do.probs)
      {
        p.plot=kobePlot(f.traj=Fmsy,
                         b.traj=Bmsy,
                         Years=yrs,
                         Titl="JABBA",
                         Probs=data.frame(x=d$JABBA$Kobe.probs[[id]]$stock,  
                                          y=d$JABBA$Kobe.probs[[id]]$harvest),
                        YrSize=txt.size)
      }
    }
    
    #SSS
    if('SSS'%in%names(d)[[pp]])
    {
      dummy=d$SSS$B.Bmsy[[id]]%>%
              filter(Scenario==Scen & year<=Last.yr.ktch.numeric)
      yrs=dummy%>%pull(year)
      Bmsy=dummy$median
      Fmsy=d$SSS$F.Fmsy[[id]]%>%
              filter(Scenario==Scen & year<=Last.yr.ktch.numeric)%>%
              pull(median)
      p.plot=kobePlot(f.traj=Fmsy[1:length(yrs)],
                     b.traj=Bmsy[1:length(yrs)],
                     Years=yrs,
                     Titl="SSS",
                     YrSize=txt.size)
      rm(yrs,Fmsy,Bmsy,dummy)
    }
    
    #SS
    if('SS'%in%names(d)[[pp]])
    {
      dummy=d$SS$B.Bmsy[[id]]%>%
              filter(Scenario==Scen & year<=Last.yr.ktch.numeric)
      yrs=dummy%>%pull(year)
      Bmsy=dummy$median
      Fmsy=d$SS$F.Fmsy[[id]]%>%
              filter(Scenario==Scen & year<=Last.yr.ktch.numeric)%>%
              pull(median)
      p.plot=kobePlot(f.traj=Fmsy[1:length(yrs)],
                    b.traj=Bmsy[1:length(yrs)],
                    Years=yrs,
                    Titl="SS",
                    YrSize=txt.size)
      rm(yrs,Fmsy,Bmsy,dummy)
    }
    
    plotlist[[pp]]=p.plot+rremove("axis.title")
  }
  #Combine plots 
  # if('DBSRA'%in%names(d) & 'CMSY'%in%names(d)& 'SSS'%in%names(d))
  # {
  #   dis.pis=c('p.DBSRA','p.CMSY','p.JABBA','p.SSS')
  #   dis.pis=dis.pis[objects.exist('p.DBSRA','p.CMSY','p.JABBA','p.SSS')]
  #   plotlist=vector('list',length(dis.pis))
  #   names(plotlist)=dis.pis
  #   for(pp in 1:length(dis.pis)) plotlist[[pp]]=dis.pis[pp] 
  #   plotlist=list(DBSRA=p.DBSRA+rremove("axis.title"),
  #                 CMSY=p.CMSY+rremove("axis.title"),
  #                 JABBA=p.JABBA+rremove("axis.title"),
  #                 SSS=p.SSS+rremove("axis.title"))
  # }else
  # {
  #   if('JABBA'%in%names(d)) plotlist=list(JABBA=p.JABBA+rremove("axis.title"))
  #   if('SS'%in%names(d))  plotlist=list(SS=p.SS+rremove("axis.title"))
  # }
  if(length(plotlist)<4)
  {
    nKl=1
    nRw=length(plotlist)
  }
  if(length(plotlist)>=4)
  {
    nKl=2
    nRw=length(plotlist)/2
  }
  
  figure <- ggarrange(plotlist=plotlist,
                      ncol=nKl,nrow=nRw,common.legend = FALSE)
  if(add.sp.nm) figure=figure+theme(plot.margin = margin(1,0,0,0, "cm"))
  figure=annotate_figure(figure,
                         bottom = text_grob(expression(B/~B[MSY]), size=22),
                         left = text_grob(expression(F/~F[MSY]), rot = 90,size=22))
  if(add.sp.nm) figure=annotate_figure(figure,
                                       fig.lab=capitalize(sp),
                                       fig.lab.pos='top.left',
                                       fig.lab.size=24)
  
  print(figure)
  return(plotlist)
}
fn.get.Kobe.plot=function(this.sp,d,NKOL,NRW,do.probs=FALSE)
{
  id=match(this.sp,names(d$rel.biom))
  Bmsy=d$B.Bmsy[id]
  Fmsy=d$F.Fmsy[id]
  plotlist=vector('list',length(Bmsy))
  for(x in 1:length(Bmsy))
  {
    dis.fmsy=Fmsy[[x]]%>%filter(year<=Last.yr.ktch.numeric)
    dis.bmsy=Bmsy[[x]]%>%filter(year<=Last.yr.ktch.numeric)
    
    if('Scenario'%in%names(dis.fmsy))
    {
      dis.fmsy=dis.fmsy%>%filter(Scenario=='S1')
      dis.bmsy=dis.bmsy%>%filter(Scenario=='S1')
    }
    
    if(!do.probs)
    {
      
      dd=kobePlot(f.traj=dis.fmsy$median,
                  b.traj=dis.bmsy$median,
                  Years=dis.bmsy$year,
                  Titl=capitalize(names(Bmsy)[x]),
                  YrSize=5)+rremove("axis.title")
    }
    
    if(do.probs)
    {
      dd=kobePlot(f.traj=dis.fmsy$median,
                  b.traj=dis.bmsy$median,
                  Years=dis.bmsy$year,
                  Titl=capitalize(names(Bmsy)[x]),
                  YrSize=6,
                  Probs=data.frame(x=d$Kobe.probs[[id[x]]]$stock,
                                   y=d$Kobe.probs[[id[x]]]$harvest))+rremove("axis.title")
    }
    plotlist[[x]]=dd
  }
  figure <- ggarrange(plotlist=plotlist,ncol=NKOL,nrow=NRW,common.legend = FALSE)
  figure=annotate_figure(figure,
                         bottom = text_grob(expression(B/~B[MSY]), size=22),
                         left = text_grob(expression(F/~F[MSY]), rot = 90,size=22))
  print(figure)
}
kobePlot_WA.Fisheries=function(f.traj,b.traj,Years,Titl,txt.col='black',YrSize=4,
                               B.lim,F.lim,B.thr,F.thr,B.tar,F.tar)
{
  dta=data.frame(x=b.traj,
                 y=f.traj,
                 yr=Years)%>%
    arrange(yr)
  Mx.F=max(F.lim,max(dta$y,na.rm=T))
  Mx.B=max(1.1,max(dta$x,na.rm=T))
  
  kobe <-dta%>%
    ggplot(aes(x, y))+    
    scale_x_continuous(limits=c(0,Mx.B),expand = c(0, 0)) +
    scale_y_continuous(limits=c(0,Mx.F),expand = c(0, 0))+
    geom_rect(xmin = 0, xmax = Mx.B, ymin = 0, ymax = Mx.F, fill = RiskColors['Severe'], alpha = 0.05)+
    geom_rect(xmin = B.lim, xmax = Mx.B, ymin = 0, ymax = F.lim, fill = RiskColors['High'], alpha = 0.05)+
    geom_rect(xmin = B.thr, xmax = Mx.B, ymin = 0, ymax = F.thr, fill = RiskColors['Medium'], alpha = 0.05)+
    geom_rect(xmin = B.tar, xmax = Mx.B, ymin = 0, ymax = F.tar, fill = RiskColors['Low'], alpha = 0.05)
  
  kobe <-kobe +
    geom_path(linetype = 2, size = 0.5,color='deepskyblue4')+
    geom_point(size=2,color='deepskyblue4')+
    geom_point(aes(x=dta[1,'x'],y=dta[1,'y']),size=4,shape=22,fill='white',alpha=.3)+
    geom_point(aes(x=dta[nrow(dta),'x'],y=dta[nrow(dta),'y']),size=4,shape=25,fill='white',alpha=.3)+      
    geom_text_repel(data=dta[1,],aes(x=x,y=y,label=yr),size=YrSize,color=txt.col)+
    geom_text_repel(data=dta[nrow(dta),],aes(x=x,y=y,label=yr),size=YrSize,color=txt.col)+
    xlab(expression(B/~B[0]))+ylab(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")))+
    labs(title = Titl)+
    theme_bw()%+replace% 
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size=16),
          axis.title = element_text(size=20),
          plot.title = element_text(size=20,hjust=0),
          legend.text = element_text(size=15),
          legend.title = element_text(size=17),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-10))
  return(kobe)
}
fn.get.Kobe.plot_appendix_WA.Fisheries=function(d,sp,Scen='S1',add.sp.nm=FALSE,txt.size=6,RF=Ref.points)
{
  id=match(sp,Keep.species)
  plotlist=vector('list',length(d))
  names(plotlist)=names(d)
  for(pp in 1:length(d))
  {
    B.ref.points=do.call(rbind,RF[[id]])%>%filter(Scenario==Scen)
    
    #JABBA
    if('JABBA'%in%names(d)[[pp]])
    {
      dummy=d$JABBA$rel.biom[[id]]%>%
        filter(Scenario==Scen & year<=Last.yr.ktch.numeric)
      yrs=dummy%>%pull(year)
      B=dummy$median
      FF=d$JABBA$f.series[[id]]%>%
        filter(Scenario==Scen & year<=Last.yr.ktch.numeric)%>%
        pull(median)
      F.ref.points=d$JABBA$estimates[[id]]%>%
        filter(Scenario==Scen & Parameter%in%c('F.target','F.threshold','F.limit'))
      
      p.plot=kobePlot_WA.Fisheries(f.traj=FF,
                                   b.traj=B,
                                   Years=yrs,
                                   Titl="JABBA",
                                   YrSize=txt.size,
                                   B.lim=B.ref.points%>%filter(Rf.pt=='Limit')%>%pull(Value),
                                   F.lim=F.ref.points%>%filter(Parameter=='F.limit')%>%pull(Median),
                                   B.thr=B.ref.points%>%filter(Rf.pt=='Threshold')%>%pull(Value),
                                   F.thr=F.ref.points%>%filter(Parameter=='F.threshold')%>%pull(Median),
                                   B.tar=B.ref.points%>%filter(Rf.pt=='Target')%>%pull(Value),
                                   F.tar=F.ref.points%>%filter(Parameter=='F.target')%>%pull(Median))
      
      
    }
    
    #SS
    if('SS'%in%names(d)[[pp]])
    {
      dummy=d$SS$rel.biom[[id]]%>%
        filter(Scenario==Scen & year<=Last.yr.ktch.numeric)
      yrs=dummy%>%pull(year)
      B=dummy$median
      FF=d$SS$F.Fmsy[[id]]%>%
        filter(Scenario==Scen & year<=Last.yr.ktch.numeric)%>%
        pull(median)
      F.ref.points=d$SS$estimates[[id]]%>%
        filter(Scenario==Scen & Parameter%in%c('F.target','F.threshold','F.limit'))
      
      p.plot=kobePlot_WA.Fisheries(f.traj=FF[1:length(yrs)],
                                   b.traj=B[1:length(yrs)],
                                   Years=yrs,
                                   Titl="SS",
                                   YrSize=txt.size,
                                   B.lim=B.ref.points%>%filter(Rf.pt=='Limit')%>%pull(Value),
                                   F.lim=F.ref.points%>%filter(Parameter=='F.limit')%>%pull(Median),
                                   B.thr=B.ref.points%>%filter(Rf.pt=='Threshold')%>%pull(Value),
                                   F.thr=F.ref.points%>%filter(Parameter=='F.threshold')%>%pull(Median),
                                   B.tar=B.ref.points%>%filter(Rf.pt=='Target')%>%pull(Value),
                                   F.tar=F.ref.points%>%filter(Parameter=='F.target')%>%pull(Median))
      rm(yrs,Fmsy,Bmsy,dummy)
    }
    
    if(length(id)>1)
    {
      plotlist[[pp]]=p.plot+rremove("axis.title")
    }else
    {
      plotlist[[pp]]=p.plot
    }
    
  }
  
  if(length(plotlist)<4)
  {
    nKl=1
    nRw=length(plotlist)
  }
  if(length(plotlist)>=4)
  {
    nKl=2
    nRw=length(plotlist)/2
  }
  
  figure <- ggarrange(plotlist=plotlist,
                      ncol=nKl,nrow=nRw,common.legend = FALSE)
  if(add.sp.nm) figure=figure+theme(plot.margin = margin(1,0,0,0, "cm"))
  # figure=annotate_figure(figure,
  #                        bottom = text_grob(expression(B/~B[0]), size=22),
  #                        left = text_grob('Fishing mortality', rot = 90,size=22))
  if(add.sp.nm) figure=annotate_figure(figure,
                                       fig.lab=capitalize(sp),
                                       fig.lab.pos='top.left',
                                       fig.lab.size=24)
  
  print(figure)
  return(plotlist)
}
fn.get.Kobe.plot_WA.Fisheries=function(this.sp,d,NKOL,NRW,RF=Ref.points)
{
  id=match(this.sp,names(d$rel.biom))
  Bmsy=d$rel.biom[id]
  Fmsy=d$f.series[id]
  plotlist=vector('list',length(Bmsy))
  for(x in 1:length(Bmsy))
  {
    dis.fmsy=Fmsy[[x]]%>%filter(year<=Last.yr.ktch.numeric)
    dis.bmsy=Bmsy[[x]]%>%filter(year<=Last.yr.ktch.numeric)
    if('Scenario'%in%names(dis.fmsy))
    {
      dis.fmsy=dis.fmsy%>%filter(Scenario=='S1')
      dis.bmsy=dis.bmsy%>%filter(Scenario=='S1')
    }
    B.ref.points=do.call(rbind,RF[id][[x]])%>%filter(Scenario=='S1')
    F.ref.points=d$estimates[id][[x]]%>%
      filter(Scenario=='S1' & Parameter%in%c('F.target','F.threshold','F.limit'))
    
    dd=kobePlot_WA.Fisheries(f.traj=dis.fmsy$median,
                             b.traj=dis.bmsy$median,
                             Years=dis.bmsy$year,
                             Titl=capitalize(names(Bmsy)[x]),
                             YrSize=5,
                             B.lim=B.ref.points%>%filter(Rf.pt=='Limit')%>%pull(Value),
                             F.lim=F.ref.points%>%filter(Parameter=='F.limit')%>%pull(Median),
                             B.thr=B.ref.points%>%filter(Rf.pt=='Threshold')%>%pull(Value),
                             F.thr=F.ref.points%>%filter(Parameter=='F.threshold')%>%pull(Median),
                             B.tar=B.ref.points%>%filter(Rf.pt=='Target')%>%pull(Value),
                             F.tar=F.ref.points%>%filter(Parameter=='F.target')%>%pull(Median))
    plotlist[[x]]=dd+rremove("axis.title")
  }
  figure <- ggarrange(plotlist=plotlist,ncol=NKOL,nrow=NRW,common.legend = FALSE)
  figure=annotate_figure(figure,
                         bottom = text_grob(expression(B/~B[0]), size=22),
                         left = text_grob(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")), rot = 90,size=22))
  print(figure)
}
kobePlot_SAFS=function(f.traj,b.traj,Years,Titl,txt.col='black',YrSize=4,
                       B.lim,F.lim)
{
  dta=data.frame(x=b.traj,
                 y=f.traj,
                 yr=Years)%>%
    arrange(yr)
  Mx.F=max(F.lim,max(dta$y,na.rm=T))*1.2
  Mx.B=max(1.01,max(dta$x,na.rm=T))
  
  kobe <-dta%>%
    ggplot(aes(x, y))+    
    scale_x_continuous(limits=c(0,Mx.B),expand = c(0, 0)) +
    scale_y_continuous(limits=c(0,Mx.F),expand = c(0, 0))+
    geom_rect(xmin = 0, xmax = B.lim, ymin = F.lim, ymax = Mx.F, fill = RiskColors['Severe'], alpha = 0.05)+
    geom_rect(xmin = 0, xmax = B.lim, ymin = 0, ymax = F.lim, fill = RiskColors['High'], alpha = 0.05)+
    geom_rect(xmin = B.lim, xmax = Mx.B, ymin = 0, ymax = F.lim, fill = RiskColors['Low'], alpha = 0.05)+
    geom_rect(xmin = B.lim, xmax = Mx.B, ymin = F.lim, ymax = Mx.F, fill = RiskColors['Medium'], alpha = 0.05)+
    annotate(geom = "text",x=mean(c(B.lim,Mx.B)),y=mean(c(F.lim,Mx.F)),label="Depleting",size=4,alpha=0.5)+
    annotate(geom = "text",x=mean(c(B.lim,Mx.B)),y=mean(c(0,F.lim)),label="Sustainable",size=4,alpha=0.5)+
    annotate(geom = "text", x = mean(c(0,B.lim)), y = mean(c(F.lim,Mx.F)), label = "Depleted",size=4,angle = 90,alpha=0.5)+
    annotate(geom = "text", x = mean(c(0,B.lim)), y = mean(c(0,F.lim)), label = "Recovering",size=4,angle = 90,alpha=0.5)
  
  kobe <-kobe +
    geom_path(linetype = 2, size = 0.5,color='deepskyblue4')+
    geom_point(size=2,color='deepskyblue4')+
    geom_point(aes(x=dta[1,'x'],y=dta[1,'y']),size=4,shape=22,fill='white',alpha=.3)+
    geom_point(aes(x=dta[nrow(dta),'x'],y=dta[nrow(dta),'y']),size=4,shape=25,fill='white',alpha=.3)+      
    geom_text_repel(data=dta[1,],aes(x=x,y=y,label=yr),size=YrSize,color=txt.col)+
    geom_text_repel(data=dta[nrow(dta),],aes(x=x,y=y,label=yr),size=YrSize,color=txt.col)+
    xlab(expression(B/~B[0]))+ylab(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")))+
    labs(title = Titl)+
    theme_bw()%+replace% 
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size=16),
          axis.title = element_text(size=20),
          plot.title = element_text(size=20,hjust=0),
          legend.text = element_text(size=15),
          legend.title = element_text(size=17),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-10))
  return(kobe)
}
fn.get.Kobe.plot_appendix_SAFS=function(d,sp,Scen='S1',add.sp.nm=FALSE,txt.size=6,RF=Ref.points)
{
  id=match(sp,Keep.species)
  plotlist=vector('list',length(d))
  names(plotlist)=names(d)
  for(pp in 1:length(d))
  {
    B.ref.points=do.call(rbind,RF[[id]])%>%filter(Scenario==Scen)
    
    #JABBA
    if('JABBA'%in%names(d)[[pp]])
    {
      dummy=d$JABBA$rel.biom[[id]]%>%
        filter(Scenario==Scen & year<=Last.yr.ktch.numeric)
      yrs=dummy%>%pull(year)
      B=dummy$median
      FF=d$JABBA$f.series[[id]]%>%
        filter(Scenario==Scen & year<=Last.yr.ktch.numeric)%>%
        pull(median)
      F.ref.points=d$JABBA$estimates[[id]]%>%
        filter(Scenario==Scen & Parameter%in%c('F.target','F.threshold','F.limit'))
      
      p.plot=kobePlot_SAFS(f.traj=FF,
                           b.traj=B,
                           Years=yrs,
                           Titl="JABBA",
                           YrSize=txt.size,
                           B.lim=B.ref.points%>%filter(Rf.pt=='Limit')%>%pull(Value),
                           F.lim=F.ref.points%>%filter(Parameter=='F.limit')%>%pull(Median))
      
      
    }
    
    #SS
    if('SS'%in%names(d)[[pp]])
    {
      dummy=d$SS$rel.biom[[id]]%>%
        filter(Scenario==Scen & year<=Last.yr.ktch.numeric)
      yrs=dummy%>%pull(year)
      B=dummy$median
      FF=d$SS$F.Fmsy[[id]]%>%
        filter(Scenario==Scen & year<=Last.yr.ktch.numeric)%>%
        pull(median)
      F.ref.points=d$SS$estimates[[id]]%>%
        filter(Scenario==Scen & Parameter%in%c('F.target','F.threshold','F.limit'))
      
      p.plot=kobePlot_SAFS(f.traj=FF[1:length(yrs)],
                           b.traj=B[1:length(yrs)],
                           Years=yrs,
                           Titl="SS",
                           YrSize=txt.size,
                           B.lim=B.ref.points%>%filter(Rf.pt=='Limit')%>%pull(Value),
                           F.lim=F.ref.points%>%filter(Parameter=='F.limit')%>%pull(Median))
      rm(yrs,Fmsy,Bmsy,dummy)
    }
    
    if(length(id)>1)
    {
      plotlist[[pp]]=p.plot+rremove("axis.title")
    }else
    {
      plotlist[[pp]]=p.plot
    }
    
  }
  
  if(length(plotlist)<4)
  {
    nKl=1
    nRw=length(plotlist)
  }
  if(length(plotlist)>=4)
  {
    nKl=2
    nRw=length(plotlist)/2
  }
  
  figure <- ggarrange(plotlist=plotlist,
                      ncol=nKl,nrow=nRw,common.legend = FALSE)
  if(add.sp.nm) figure=figure+theme(plot.margin = margin(1,0,0,0, "cm"))
  # figure=annotate_figure(figure,
  #                        bottom = text_grob(expression(B/~B[0]), size=22),
  #                        left = text_grob('Fishing mortality', rot = 90,size=22))
  if(add.sp.nm) figure=annotate_figure(figure,
                                       fig.lab=capitalize(sp),
                                       fig.lab.pos='top.left',
                                       fig.lab.size=24)
  
  print(figure)
  return(plotlist)
}
fn.get.f.ref.points=function(Report,propTar=Tar.prop.bmsny,propLim=Lim.prop.bmsy)
{
  SSB_MSY=Report$derived_quants%>%filter(Label=='SSB_MSY')%>%pull(Value)
  
  B.threshold=SSB_MSY
  B.target=SSB_MSY*propTar
  B.limit=SSB_MSY*propLim
  
  Equil.yield=Report$equil_yield%>%
    distinct(SSB,.keep_all=T)%>%
    dplyr::select(SPRloop,Iter,SSB,Tot_Catch,Fmult,F_report)%>%
    arrange(SSB)
  id.threshold=which.min(abs(Equil.yield$SSB - B.threshold))
  id.target=which.min(abs(Equil.yield$SSB - B.target))
  id.limit=which.min(abs(Equil.yield$SSB - B.limit))
  F.threshold=Equil.yield$F_report[id.threshold]
  #Report$derived_quants%>%filter(Label=='annF_MSY')%>%pull(Value)
  F.target=Equil.yield$F_report[id.target]
  F.limit=Equil.yield$F_report[id.limit]
  
  dd=Equil.yield[c(id.target,id.threshold,id.limit),]%>%
    mutate(Labl=paste0(c('F[Tar]','F[Thr]','F[Lim]'),paste0(' (',round(F_report,3),')')))
  
  Kls=c('forestgreen','orange','red')
  names(Kls)=dd$Labl
  
  Equil.yield%>%
    ggplot(aes(SSB,Tot_Catch))+
    geom_line()+
    ylab('Catch (tonnes)')+xlab('Spawning stock biomass (tonnes)')+
    theme_PA()+theme(legend.position="none")+
    geom_text_repel(data=dd,aes(SSB,Tot_Catch,label = Labl,color=Labl),
                    box.padding=3,size=5, parse = TRUE)+
    scale_colour_manual(values=Kls,labels=names(Kls))
  ggsave(paste(this.wd1,'F reference points.tiff',sep='/'),width = 6,height = 6, dpi = 300, compression = "lzw")
  return(list(F.target=F.target,F.threshold=F.threshold,F.limit=F.limit)) 
}
# Store consequence-likelihood  ------------------------------------------------------
fn.get.cons.like=function(lista,Mod.Avrg.weight=NULL)
{
  #get likelihood for each model
  out=vector('list',length(lista))
  names(out)=names(lista)
  out.bmsy=out
  for(l in 1:length(lista))
  {
    #Depletion
    d=lista[[l]]$probs.rel.biom
    d=compact(d)
    out1=vector('list',length(d))
    names(out1)=names(d)
    for(x in 1:length(d))
    {
      dumx=d[[x]][[1]]$probs
      if(!is.null(dumx))
      {
        out1[[x]]=dumx%>%filter(Scenario=='S1')%>%mutate(Model=names(out)[l])
      }
    }
    out[[l]]=out1
    
    #B over Bmsy
    d=lista[[l]]$probs.B.Bmsy
    d=compact(d)
    out1=vector('list',length(d))
    names(out1)=names(d)
    for(x in 1:length(d))
    {
      dumx=d[[x]][[1]]$probs
      if(!is.null(dumx))
      {
        out1[[x]]=dumx%>%filter(Scenario=='S1')%>%mutate(Model=names(out)[l])
      }
    }
    out.bmsy[[l]]=out1
  }
  
  #get weighted likelihood   
  out2=vector('list',length(out[[1]]))
  names(out2)=names(out[[1]])
  out2.bmsy=vector('list',length(out.bmsy[[1]]))
  names(out2.bmsy)=names(out.bmsy[[1]])
  if(!is.null(Mod.Avrg.weight))
  {
    for(i in 1:length(out[[1]]))  
    {
      dis.r.range=do.call(rbind,r.groups)%>%   
        data.frame%>%
        mutate(Get=between(r.list[i],.[[1]] , .[[2]]))%>%
        filter(Get==TRUE)
      Wei=Mod.Avrg.weight%>%
        filter(r.group==row.names(dis.r.range))%>%
        dplyr::select(-r.group)
      if('SSS'%in%names(Catch_only) & !'SSS'%in%Wei$Model)
      {
        Wei.SSS=Wei%>%filter(Model=='DBSRA')%>%mutate(Model='SSS')
        Wei=rbind(Wei,Wei.SSS)
      }
      #Depletion
      dd=map(out, ~.x[[i]])
      dd=do.call(rbind,dd)%>%
        left_join(Wei,by=c('Model'))%>%
        group_by(Range)%>%
        summarise(Probability=weighted.mean(Probability,Weight))%>%
        ungroup()%>%
        mutate(Probability=Probability/sum(Probability))%>%
        mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
        arrange(Range)
      out2[[i]]=dd
      
      #B over Bmsy
      dd=map(out.bmsy, ~.x[[i]])
      dd=do.call(rbind,dd)%>%
        left_join(Wei,by=c('Model'))%>%
        group_by(Range)%>%
        summarise(Probability=weighted.mean(Probability,Weight))%>%
        ungroup()%>%
        mutate(Probability=Probability/sum(Probability))%>%
        mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
        arrange(Range)
      out2.bmsy[[i]]=dd
    }
  }
  if(is.null(Mod.Avrg.weight))
  {
    this.sp=vector('list',length(out))
    for(x in 1:length(out)) this.sp[[x]]=names(out[[x]])
    this.sp=unlist(this.sp)
    #Depletion
    for(i in 1:length(out[[1]]))  
    {
      ii=match(names(out2)[i],this.sp)
      if(!is.na(ii))
      {
        out2[[i]]=do.call(rbind,map(out, ~.x[[ii]]))%>%
                    mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
                    arrange(Range)%>%
          dplyr::select(Range,Probability)
      }
    }
    out2=compact(out2)
    
    #B over Bmsy
    for(i in 1:length(out.bmsy[[1]]))  
    {
      ii=match(names(out2.bmsy)[i],this.sp)
      if(!is.na(ii))
      {
        out2.bmsy[[i]]=do.call(rbind,map(out.bmsy, ~.x[[ii]]))%>%
          mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
          arrange(Range)%>%
          dplyr::select(Range,Probability)
      }
    }
    out2.bmsy=compact(out2.bmsy)
  }
  return(list(Depletion=out2,B.over.Bmsy=out2.bmsy))
}
# Weight of Evidence assessment   ------------------------------------------------------
fun.risk.spatial.dist=function(d,Optns=c('increasing','stable','decreasing'),Mult1=.9,Mult2=.75,outpath)
{
  DD=expand.grid(Effort=Optns,
                 Blocks.fished=Optns,
                 Blocks.fished.with.catch=Optns)%>%
    mutate(Effort.score=case_when(Effort=='increasing' ~ 3,
                                  Effort=='stable' ~ 2,
                                  Effort=='decreasing' ~ 1),
           Blocks.fished.score=case_when(Blocks.fished=='increasing' ~ 3,
                                         Blocks.fished=='stable' ~ 2,
                                         Blocks.fished=='decreasing' ~ 1),
           Blocks.fished.with.catch.score=case_when(Blocks.fished.with.catch=='increasing' ~ 1,
                                                    Blocks.fished.with.catch=='stable' ~ 2,
                                                    Blocks.fished.with.catch=='decreasing' ~ 3),
           Score=Effort.score*Blocks.fished.score*Blocks.fished.with.catch.score,
           Score=case_when(Score==18 & Blocks.fished.with.catch=='decreasing' & Effort=='increasing'~27*Mult2,
                           Score==18 & Blocks.fished.with.catch=='decreasing' & Blocks.fished=='increasing'~27*Mult1,
                           Score==12 & Blocks.fished.with.catch=='decreasing' ~18*Mult1,
                           Score==12 & Blocks.fished.with.catch=='stable' & Blocks.fished=='increasing' ~18*Mult2,
                           Score==9 & Blocks.fished.with.catch=='decreasing' & Blocks.fished=='increasing' ~12*0.95,
                           Score==9 & Blocks.fished.with.catch=='decreasing' & Effort=='increasing'~12*0.85,
                           TRUE~Score),
           Consequence=round(rescale((4*Score)/max(Score),to=c(1,4))),
           Likelihood=4)
  
  out.table=DD%>%
    arrange(Consequence)%>%
    mutate(across(where(is.factor), as.character))%>%
    dplyr::select(Effort,Blocks.fished,Blocks.fished.with.catch,Consequence,Likelihood)%>%
    mutate(Effort=capitalize(Effort),
           Blocks.fished=capitalize(Blocks.fished),
           Blocks.fished.with.catch=capitalize(Blocks.fished.with.catch))
  write.csv(out.table,paste(outpath,'Spatial_LoE_Table.csv',sep='/'),row.names = F)
  
  out.table%>%
    flextable()%>%
      set_header_labels(Blocks.fished="Blocks\nfished",
                        Blocks.fished.with.catch="Blocks fished\nwith catch")%>%
    width(j='Blocks.fished.with.catch',width=3,unit='cm')%>%
    align(align="left")%>% save_as_docx( path = paste(outpath,'Spatial_LoE_Table.docx',sep='/'))
  
  
  risk.tab=d%>%
    mutate(Blocks.fished.with.catch=ifelse(Market=='decreasing' | Conservation=='yes','stable',Blocks.fished.with.catch))%>%
    dplyr::select(Species,Effort,Blocks.fished,Blocks.fished.with.catch)%>%
    left_join(DD%>%dplyr::select(Effort,Blocks.fished,Blocks.fished.with.catch,Consequence,Likelihood),
              by=c('Effort','Blocks.fished','Blocks.fished.with.catch'))%>%
    dplyr::select(Species,Consequence,Likelihood)%>%
    mutate(Probability=NA,
           w=1,
           Species=capitalize(Species),
           dummy=paste(Species,Consequence,Likelihood))
  
  add.dummy=data.frame(Species=rep(risk.tab$Species,each=4))%>%
    mutate(Consequence=rep(seq(1,4),nrow(risk.tab)),
           Likelihood=rep(c(4,2,1,1),nrow(risk.tab)),
           Probability=NA,
           w=1,
           dummy=paste(Species,Consequence,Likelihood))%>%
    filter(!dummy%in%risk.tab$dummy)
  risk.tab=rbind(risk.tab,add.dummy)%>%
    dplyr::select(-dummy)%>%
    arrange(Species)
  
  
  OLD.way=FALSE
  if(OLD.way)
  {
    risk.tab=d%>%
      mutate(Like.C1=case_when(Conservation=='yes' & Catch=='zero' & Blocks.with.catch%in%c('zero')~4,
                               Conservation=='yes' & Catch=='decreasing' & Blocks.with.catch%in%c('decreasing')~3,
                               Market%in%c('zero','decreasing') & Catch=='decreasing' & Management_driven.effort=='decreasing'~4,
                               Market=='increasing' & Catch=='decreasing' & Blocks.with.catch=='increasing' & Management_driven.effort=='decreasing'~4,
                               Market%in%c('stable') & Catch=='decreasing' & Blocks.with.catch=='stable' & Management_driven.effort=='decreasing'~4,
                               TRUE~NA_real_),
             Like.C2=case_when(Conservation=='yes' & Catch=='zero' & Blocks.with.catch%in%c('zero')~2,
                               Conservation=='yes' & Catch=='decreasing' & Blocks.with.catch%in%c('decreasing')~2,
                               Market%in%c('zero','decreasing') & Catch=='decreasing' & Management_driven.effort=='decreasing'~2,
                               Market=='increasing' & Catch=='decreasing' & Blocks.with.catch=='increasing' & Management_driven.effort=='decreasing'~2,
                               Market%in%c('stable') & Catch=='decreasing' & Blocks.with.catch=='stable' & Management_driven.effort=='decreasing'~2,
                               TRUE~NA_real_),
             Like.C3=case_when(Conservation=='yes' & Catch=='zero' & Blocks.with.catch%in%c('zero')~1,
                               Conservation=='yes' & Catch=='decreasing' & Blocks.with.catch%in%c('decreasing')~1,
                               Market%in%c('zero','decreasing') & Catch=='decreasing' & Management_driven.effort=='decreasing'~1,
                               Market=='increasing' & Catch=='decreasing' & Blocks.with.catch=='increasing' & Management_driven.effort=='decreasing'~1,
                               Market%in%c('stable') & Catch=='decreasing' & Blocks.with.catch=='stable' & Management_driven.effort=='decreasing'~1,
                               TRUE~NA_real_),
             Like.C4=case_when(Conservation=='yes' & Catch=='zero' & Blocks.with.catch%in%c('zero')~1,
                               Conservation=='yes' & Catch=='decreasing' & Blocks.with.catch%in%c('decreasing')~1,
                               Market%in%c('zero','decreasing') & Catch=='decreasing' & Management_driven.effort=='decreasing'~1,
                               Market=='increasing' & Catch=='decreasing' & Blocks.with.catch=='increasing' & Management_driven.effort=='decreasing'~1,
                               Market%in%c('stable') & Catch=='decreasing' & Blocks.with.catch=='stable' & Management_driven.effort=='decreasing'~1,
                               TRUE~NA_real_))%>%
      dplyr::select(Species,Like.C1,Like.C2,Like.C3,Like.C4)%>%
      gather(Consequence,Likelihood,-Species)%>%
      mutate(Consequence=as.integer(gsub('[Like.C]','',Consequence)),
             Probability=NA,
             w=1,
             Species=capitalize(Species))%>%
      arrange(Species)
    
  }
  
  return(risk.tab)
  
}
fun.risk.CoMs=function(d)
{
  dd=vector('list',nrow(d))
  for(q in 1:length(dd))
  {
    d[q,]
    dd[[q]]=data.frame(Species=rep(d[q,]$Species,4),
                       Consequence=1:4,
                       Likelihood=1,
                       Probability=0)%>%
      mutate(Likelihood=case_when(d[q,]$P.catch.below.tar==1 & Consequence==1~4,
                                  d[q,]$P.between.thre.tar==1 & Consequence==2~4,
                                  d[q,]$P.between.lim.thre==1 & Consequence==3~4,
                                  d[q,]$P.catch.above.limit==1 & Consequence==4~4,
                                  TRUE~Likelihood),
             Probability=case_when(d[q,]$P.catch.below.tar==1 & Consequence==1~1,
                                   d[q,]$P.between.thre.tar==1 & Consequence==2~1,
                                   d[q,]$P.between.lim.thre==1 & Consequence==3~1,
                                   d[q,]$P.catch.above.limit==1 & Consequence==4~1,
                                   TRUE~Probability),
             w=1)
  }
  return(do.call(rbind,dd))
}
fn.risk=function(d,w)
{
  d=d%>%
    mutate(w=w,
           Consequence=case_when(Range=='>tar'~1,
                                 Range=='thr.tar'~2,
                                 Range=='lim.thr'~3,
                                 Range=='<lim'~4))
  for(n in 1:nrow(d))
  {
    d$Likelihood[n]=which(unlist(lapply(Like.ranges,function(x) check.in.range(d$Probability[n],x,fatal=F))))
  }
  return(d%>%dplyr::select(Species,Consequence,Likelihood,Probability,w))
}
Integrate.LoE=function(Cons.Like.tab,criteriaMinMax,plot.data,LoE.weights,Normalised)
{
  #Set up preference table by converting Cons-Like to Risk scores
  Preference.Table=as.data.frame(matrix(0,nrow=5,ncol=ncol(Cons.Like.tab)))
  colnames(Preference.Table)=colnames(Cons.Like.tab)
  rownames(Preference.Table)=c("Negligible","Low","Medium","High","Severe")
  for(i in 1:ncol(Preference.Table))
  {
    dd=Cons.Like.tab[,i]
    if(max(dd[1:2])<=2) Preference.Table["Negligible",i]=max(dd[1:2])
    if(max(dd)>2 & max(dd)<=4) Preference.Table["Low",i]=4
    if(max(dd[2:4])>4 & max(dd[2:4])<=8) Preference.Table["Medium",i]=max(dd[2:4])
    if(max(dd[3:4])>8 & max(dd[3:4])<=12) Preference.Table["High",i]=max(dd[3:4])
    if(dd[4]>12 & dd[4]<=16) Preference.Table["Severe",i]=dd[4]
  }
  
  #Maximise or minimise each criteria?
  criteriaMinMax=rep(criteriaMinMax,ncol(Preference.Table))
  names(criteriaMinMax) <- colnames(Preference.Table)
  
  #display data
  if(plot.data)plotRadarPerformanceTable(Preference.Table, criteriaMinMax,overlay=FALSE, bw=TRUE, lwd =5)
  
  # Normalization of the performance table
  normalizationTypes <- rep("percentageOfMax",ncol(Preference.Table))
  names(normalizationTypes) <- colnames(Preference.Table)
  if(Normalised=="YES") nPreference.Table <- normalizePerformanceTable(Preference.Table,normalizationTypes)
  if(Normalised=="NO") nPreference.Table=Preference.Table
  
  # Calculate weighted sum
  names(LoE.weights) <- colnames(nPreference.Table)
  weighted.sum<-weightedSum(nPreference.Table,LoE.weights)
  
  # Rank the scores of the alternatives
  rank.score=sort(rank(-weighted.sum))
  
  # overall risk
  risk=which(rank.score==min(rank.score))
  risk=as.character(max(Order[match(names(risk),Order)]))
  
  return(list(weighted.sum=weighted.sum,rank.score=rank.score,risk=risk,Cons.Like.tab=Cons.Like.tab))
}
fn.cons.po=function(low,up) c(low, tail(up, 1), rev(up), low[1])  
fn.each.LoE.risk=function(N.sp,CX.axis)
{
  X.rng=1:N.sp
  x.Vec <-  fn.cons.po(0:(N.sp+1),0:(N.sp+1))
  nn=N.sp+2
  negigible.Vec <- fn.cons.po(rep(0,nn),rep(2,nn))
  low.Vec <- fn.cons.po(rep(2,nn),rep(4,nn))
  medium.Vec <- fn.cons.po(rep(4,nn),rep(8,nn))
  high.Vec <- fn.cons.po(rep(8,nn),rep(12,nn))
  severe.Vec <- fn.cons.po(rep(12,nn),rep(16,nn))
  plot(X.rng,xlim=c(0,16),ylim=c(0,N.sp+1),xaxs="i",yaxs="i",
       col="transparent",ylab="",xlab="",xaxt='n',yaxt='n')
  polygon(negigible.Vec, x.Vec, col = RiskColors['Negligible'], border = "transparent")
  polygon(low.Vec, x.Vec, col = RiskColors['Low'], border = "transparent")
  polygon(medium.Vec, x.Vec, col = RiskColors['Medium'], border = "transparent")
  polygon(high.Vec, x.Vec, col = RiskColors['High'], border = "transparent")
  polygon(severe.Vec, x.Vec, col = RiskColors['Severe'], border = "transparent")
  
  axis(1,at=c(mean(negigible.Vec),mean(low.Vec),mean(medium.Vec),mean(high.Vec),mean(severe.Vec)),
       labels=c("Negl.","Low","Med.","High","Severe"),cex.axis=CX.axis)
  box()
}
fn.overall.risk=function(N,RISK,sp,CX=1.5)
{
  x.Vec <-  c(1,3,3,1)
  y.Vec <- c(rep((s-.5),2),rep((s+.5),2))
  if(RISK=="Negligible") CL = RiskColors['Negligible']
  if(RISK=="Low") CL = RiskColors['Low']
  if(RISK=="Medium") CL =RiskColors['Medium']
  if(RISK=="High") CL =RiskColors['High']
  if(RISK=="Severe") CL =RiskColors['Severe']
  polygon(x.Vec, y.Vec, col = CL, border = "transparent")
  text(1.5,mean(y.Vec),sp,cex=CX)
}
fn.risk.figure=function(d,Risk.colors,out.plot)   
{
  d=d%>%
    mutate(Score=round(Consequence*Likelihood),
           Cons=case_when(Consequence==1~'Minor',
                          Consequence==2~'Moderate',
                          Consequence==3~'High',
                          Consequence==4~'Major'),
           Like=case_when(Likelihood==1~'Remote',
                          Likelihood==2~'Unlikely',
                          Likelihood==3~'Possible',
                          Likelihood==4~'Likely'),
           Cons.lbl=paste(paste(Cons,'\n',sep=''),' (',Consequence,')',sep=''),
           Cons.lbl=factor(Cons.lbl,levels=rev(c("Minor\n (1)","Moderate\n (2)","High\n (3)","Major\n (4)"))),
           Like.lbl=paste(paste(Like,'\n',sep=''),' (',Likelihood,')',sep=''),
           Like.lbl=factor(Like.lbl,c("Remote\n (1)","Unlikely\n (2)","Possible\n (3)","Likely\n (4)")),
           Risk=case_when(Score<=2~'Negligible',
                          Score>2 & Score<=4~'Low',
                          Score>4 & Score<=8~'Medium',
                          Score>8 & Score<12~'High',
                          Cons=='High' & Like=='Likely'~'High',
                          Cons=='Major' & Like=='Possible'~'Severe',
                          Score>12~'Severe'),
           Risk=factor(Risk,levels=c('Negligible','Low','Medium','High','Severe')))%>%
    group_by(Species)%>%
    mutate(Max.val=max(Score))%>%
    ungroup()%>%
    mutate(Max.val=ifelse(Max.val==Score,Max.val,''),
           id = row_number())
  iid.dups=which(!d$Max.val=='')
  dumi=d[iid.dups,]
  tab1=table(dumi$Species)
  if(any(tab1>1))
  {
    drop.dis=dumi%>%
      filter(Species%in%names(tab1)[tab1>1])%>%
      group_by(Species)%>%
      mutate(Max.cons=max(Consequence))%>%
      filter(!Consequence==Max.cons)%>%
      pull(id)
    iid.dups=iid.dups[-match(drop.dis,iid.dups)]
    d=d%>%
      mutate(Max.val=ifelse(id%in%drop.dis,'',Max.val))
  }
  
  if(out.plot)
  {
    NN=n2mfrow(length(unique(d$Species)))
    if(NN[2]<=3)
    {
      HJUST=1
      Strip.size=16
      Axis.titl=19
      Axis.txt=14
    }
      
    if(NN[2]>3)
    {
      HJUST=1
      Strip.size=14
      Axis.titl=17
      Axis.txt=10
    }
      
    p=d%>%
      ggplot(aes(Like.lbl,Cons.lbl,  width = w,height=w)) +
      geom_tile(aes(fill = Risk))+
      facet_wrap(~Species,nrow=NN[1],ncol=NN[2],strip.position="bottom")+
      scale_fill_manual(values=Risk.colors)+
      theme_bw()%+replace% 
      theme(legend.position = "none",
            panel.background = element_rect(fill="grey95"),
            panel.border = element_rect(colour = "grey75", fill=NA, size=1.15),
            panel.grid.major = element_blank(),    
            panel.grid.minor = element_blank(),    
            axis.line = element_line(colour = "grey75"),
            strip.background = element_rect(fill = "transparent",colour = "transparent"),
            strip.text = element_text(size = Strip.size),
            axis.title = element_text(size = Axis.titl),
            axis.text = element_text(size = Axis.txt),
            axis.text.y=element_text(hjust=0.5),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            plot.caption = element_text(hjust = HJUST,size=10,vjust=-3))+
      xlab('Likelihood')+ylab('Consequence')+
      scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0),position = "top")+
      geom_hline(yintercept=c(1.5,2.5,3.5),color='grey75')+
      geom_vline(xintercept=c(1.5,2.5,3.5),color='grey75')+
      geom_text(aes(Like.lbl,Cons.lbl,label=Max.val),size=6,fontface='bold')+
       labs(caption=expression(atop('Cons.: Minor, B'[Cur]*'>B'[Tar]*'; Moderate, B'[Tar]*'>B'[Cur]*'>B'[Thr]*
                                      '; High, B'[Thr]*'>B'[Cur]*'>B'[Lim]*'; Major, B'[Lim]*'>B'[Cur]*
                                    '   Like.: Remote, <5%, Unlikely, 5-20%, Possible, 20-50%, Likely, >50%')))
      # labs(caption=expression(atop('Consequence: Minor, B'[Cur]*'>B'[Tar]*'; Moderate, B'[Tar]*'>B'[Cur]*'>B'[Thr]*
      #                                '; High, B'[Thr]*'>B'[Cur]*'>B'[Lim]*'; Major, B'[Lim]*'>B'[Cur],
      #                              'Likelihood: Remote, <5%, Unlikely, 5-20%, Possible, 20-50%, Likely, >50%                                  ')))
    print(p)
  }
  return(d[iid.dups,]%>%dplyr::select(-c(w,Cons.lbl,Like.lbl,Max.val,id)))
}


#---Old WoE risk integration  ------------------------------------------------- 

do.old.WOE=FALSE
if(do.old.WOE)
{
  #1.2. catch spatio-temporal    
  stopped.landing=c('angel sharks','sawsharks',
                    "spurdogs","grey nurse shark")  #landings ceased/decreased due to lack of market for spurdogs, 
  # and marine heatwave for angel and sawsharks (industry consultation AMM 2022)
  #or protection for greynurse shark
  Risk.spatial.temporal.ktch=vector('list',length(Store.spatial.temporal.ktch))
  names(Risk.spatial.temporal.ktch)=names(Store.spatial.temporal.ktch)
  for(s in 1: length(Risk.spatial.temporal.ktch))
  {
    sp.dat=Store.spatial.temporal.ktch[[s]]$prop.ktch_over.fished.bocks%>%
      mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
      filter(year<=(Last.yr.ktch.numeric))%>%
      filter(year>=(Last.yr.ktch.numeric-4))%>%  #average last 5 years  
      group_by(SNAME)%>%
      summarise(Mean=mean(prop))%>%
      ungroup()%>%
      rename(Species=SNAME)%>%
      filter(!Species%in%stopped.landing)
    dummy=vector('list',nrow(sp.dat))
    names(dummy)=sp.dat$Species
    for(x in 1:length(dummy))
    {
      dd=sp.dat%>%filter(Species==names(dummy)[x])%>%pull(Mean)
      if(dd==1) d=data.frame(Max.Risk.Score=c(2,0,0,0))
      if(dd<1 & dd>=.75) d=data.frame(Max.Risk.Score=c(0,4,0,0))
      if(dd<.75 & dd>=.5) d=data.frame(Max.Risk.Score=c(0,0,6,0))
      if(dd<.5 & dd>=.25) d=data.frame(Max.Risk.Score=c(0,0,12,0))
      if(dd<.25) d=data.frame(Max.Risk.Score=c(0,0,0,16))
      dummy[[x]]=d
    }
    Risk.spatial.temporal.ktch[[s]]=dummy
  }
  
  #1.3. overall effort-management               
  Risk.effort.mangmnt=Risk.spatial.temporal.ktch
  Rel.eff.s=Effort.monthly
  Rel.eff.n=rbind(Effort.monthly.north%>%rename(Total='Hook hours')%>%dplyr::select(FINYEAR,Total),  #add 0 effort if no fishing
                  data.frame(FINYEAR=Effort.monthly$FINYEAR[which(!Effort.monthly$FINYEAR%in%Effort.monthly.north$FINYEAR)],
                             Total=0))%>%
    arrange(FINYEAR)
  Yr.max.eff.s=Rel.eff.s%>%filter(Total==max(Total))%>%pull(FINYEAR)
  Yr.max.eff.n=Rel.eff.n%>%filter(Total==max(Total))%>%pull(FINYEAR)
  Rel.eff.n=mean(Rel.eff.n$Total[(nrow(Rel.eff.n)-4):nrow(Rel.eff.n)])/max(Rel.eff.n$Total) #average last 5 years
  Rel.eff.s=mean(Rel.eff.s$Total[(nrow(Rel.eff.s)-4):nrow(Rel.eff.s)])/max(Rel.eff.s$Total)
  REgn=Tot.ktch%>% 
    filter(Name%in%names(List.sp))%>%
    group_by(Name,Region)%>%
    summarise(Tot=sum(LIVEWT.c))
  for(s in 1:length(Risk.effort.mangmnt))
  {
    dat=Store.spatial.temporal.ktch[[s]]$Fished.blks.Fishery
    dat.max.n=dat%>%filter(FINYEAR==Yr.max.eff.n & Fishery=='Northern')%>%pull(Tot)
    dat.max.s=dat%>%filter(FINYEAR==Yr.max.eff.s & Fishery=='Southern')%>%pull(Tot)
    dat=dat%>%  
      filter(year<=(Last.yr.ktch.numeric))%>%
      filter(year>=(Last.yr.ktch.numeric-4))%>%  #average last 5 years  
      group_by(Fishery)%>%
      summarise(Mean=mean(Tot))%>%
      ungroup
    
    Rel.blk.n=dat%>%filter(Fishery=='Northern')%>%pull(Mean)/dat.max.n
    Rel.blk.s=dat%>%filter(Fishery=='Southern')%>%pull(Mean)/dat.max.s
    
    Eff.blk.n=Rel.eff.n*Rel.blk.n
    Eff.blk.s=Rel.eff.s*Rel.blk.s
    
    dis.sp=Lista.sp.outputs[-match('additional.sp',names(Lista.sp.outputs))][[s]]
    dummy=vector('list',length(dis.sp))
    names(dummy)=dis.sp
    for(x in 1:length(dummy))
    {
      Which.ef=REgn%>%filter(Name==names(dummy)[x])%>%pull(Region)
      Which.ef=paste(Which.ef,collapse='-')
      
      if(Which.ef=='North') Which.ef=Eff.blk.n
      if(Which.ef=='South') Which.ef=Eff.blk.s
      if(Which.ef=='North-South') Which.ef=max(c(Eff.blk.s,Eff.blk.n))
      
      if(Which.ef==0) d=data.frame(Max.Risk.Score=c(2,0,0,0))
      if(Which.ef>0 & Which.ef<=.25) d=data.frame(Max.Risk.Score=c(0,4,0,0))
      if(Which.ef>.25 & Which.ef<=.5) d=data.frame(Max.Risk.Score=c(0,0,6,0))
      if(Which.ef>.5 & Which.ef<=.75) d=data.frame(Max.Risk.Score=c(0,0,12,0))
      if(Which.ef>.75) d=data.frame(Max.Risk.Score=c(0,0,0,16))
      dummy[[x]]=d
    }
    Risk.effort.mangmnt[[s]]=dummy
  }
  
  #1.4. population dynamics models
  
  Risk.tab=data.frame(Consequence=paste("C",1:4,sep=""),
                      L1=NA,L2=NA,L3=NA,L4=NA,
                      Max.Risk.Score=NA)
  
  Risk.COM=vector('list',N.sp)
  names(Risk.COM)=names(List.sp)
  Risk.JABBA=Risk.COM
  if(Do.integrated) Risk.integrated=Risk.COM
  for(s in 1: N.sp)
  {
    #COM
    if(exists("Store.cons.Like_COM"))
    {
      Risk.COM[[s]]=fn.risk(likelihood=Store.cons.Like_COM[[s]])
    }
    
    #JABBA
    if(exists('Store.cons.Like_JABBA'))
    {
      Risk.JABBA[[s]]=fn.risk(likelihood=Store.cons.Like_JABBA[[s]]) 
    }
    
    #Integrated model
    if(exists('Store.cons.Like_Age.based'))
    {
      Risk.integrated[[s]]=fn.risk(likelihood=Store.cons.Like_Age.based[[s]])
    }
    
    
  }
  
  
  #2. Integrate the risk from each line of evidence    
  
  #note: Use a weighted sum to aggregate the Risk Categories form the alternative lines of evidence
  #      Normalize each criterion by dividing by the highest value of each criterion. 
  #      Assign weights to each criteria 
  Order=c('Negligible','Low','Medium','High','Severe')
  Order <- factor(Order,ordered = TRUE,levels = Order)
  
  
  All.sp=sort(c(Drop.species,Keep.species))
  Overall.risk=vector('list',length(All.sp))
  names(Overall.risk)=All.sp
  for(s in 1:length(All.sp))
  {
    nm=All.sp[s]
    nm2=nm
    if(nm=="smooth hammerhead") nm2='hammerheads'
    iid=ifelse(nm%in%Other.species,"Other.sp",ifelse(nm%in%names(Indicator.species),"Indicator.sp","Drop.species"))  
    iid=match(iid,names(Risk.spatial.temporal.ktch))
    dummy=list(psa=Risk.PSA[[match(nm,names(Risk.PSA))]]$Max.Risk.Score,
               sptemp=Risk.spatial.temporal.ktch[[iid]][[match(nm2,names(Risk.spatial.temporal.ktch[[iid]]))]]$Max.Risk.Score,
               efman=Risk.effort.mangmnt[[iid]][[match(nm,names(Risk.effort.mangmnt[[iid]]))]]$Max.Risk.Score,
               COM=Risk.COM[[match(nm,names(Risk.COM))]]$Max.Risk.Score,
               JABBA=Risk.JABBA[[match(nm,names(Risk.JABBA))]]$Max.Risk.Score)
    if(Do.integrated)dummy$integrated=Risk.integrated[[match(nm,names(Risk.integrated))]]$Max.Risk.Score
    dummy=dummy[!sapply(dummy, is.null)]
    NMs=names(dummy)
    dummy=do.call(cbind,dummy)
    colnames(dummy)=NMs
    rownames(dummy)=paste("C",1:4,sep='')
    
    Overall.risk[[s]]=Integrate.LoE(Cons.Like.tab=dummy,
                                    criteriaMinMax <- "max",
                                    plot.data=FALSE,
                                    LoE.weights = LoE.Weights[match(colnames(dummy),names(LoE.Weights))],
                                    Normalised="YES")
  }
  
  
  #3. Display overall risk for each species
  LoE.col=c(psa="coral3", sptemp="grey50", efman="burlywood4",
            COM="grey80", JABBA="grey35", integrated="black")
  LoE.names=c(psa="PSA", sptemp="Blocks fished", efman="Effort & Manag.",
              COM="COM", JABBA="BSSPM", integrated="Int. mod.")
  if(!Do.integrated)
  {
    LoE.col=LoE.col[-match('integrated',names(LoE.col))]
    LoE.names=LoE.names[-match('integrated',names(LoE.names))]
  }
  All.Sp.risk.ranking=factor(unlist(lapply(Overall.risk, '[[', 'risk')),levels=levels(Order))  
  All.Sp.risk.ranking=names(sort(All.Sp.risk.ranking))
  All.Sp.risk.ranking=All.Sp.risk.ranking[!duplicated(All.Sp.risk.ranking)]
  
  
  #3.1. Other species
  Sp.risk.ranking=subset(All.Sp.risk.ranking,!All.Sp.risk.ranking%in%names(Indicator.species))
  fn.fig(paste(Rar.path,"Risk_Other.sp",sep='/'), 2000, 2400)
  par(mar=c(.5,.5,3,1),oma=c(3,13.5,.5,.1),las=1,mgp=c(1,.5,0),cex.axis=1.5,xpd=TRUE)
  layout(matrix(c(rep(1,6),rep(2,3)),ncol=3))
  
  #3.1.1. Risk for each line of evidence
  fn.each.LoE.risk(N.sp=length(Sp.risk.ranking),CX.axis=0.95) #main plot
  for(s in 1:length(Sp.risk.ranking))
  {
    nm=Sp.risk.ranking[s]
    dummy=Overall.risk[[match(nm, names(Overall.risk))]]$Cons.Like.tab
    if(ncol(dummy)==1) Adjst=0
    if(ncol(dummy)>1) Adjst=seq(-.3,.3,length.out = ncol(dummy))
    
    for(nn in 1:ncol(dummy))
    {
      segments(0,s+Adjst[nn],max(dummy[,nn]),s+Adjst[nn],
               lwd=3.75,lend=1,col=LoE.col[match(colnames(dummy)[nn],names(LoE.col))])
    }
    
  }
  axis(2,1:length(Sp.risk.ranking),capitalize(Sp.risk.ranking),cex.axis=.8)
  mtext("Risk score",1,cex=1.25,line=2)
  legend(-1.5,length(Sp.risk.ranking)+6,LoE.names[1:3],
         bty='n',col=LoE.col[1:3],lty=1,lwd=3,horiz = T,cex=1.20,
         text.width=c(0,1,2.25))
  legend(-1.5,length(Sp.risk.ranking)+4.5,LoE.names[4:length(LoE.names)],bty='n',cex=1.20,
         col=LoE.col[4:length(LoE.col)],lty=1,lwd=3,horiz = T,text.width=c(0,1.5))
  
  #3.1.2. Overall risk
  plot(0:1,ylim=c(0,length(Sp.risk.ranking)+1),fg='white',xaxs="i",yaxs="i",
       col="transparent",ylab="",xlab="",xaxt='n',yaxt='n')
  for(s in 1:length(Sp.risk.ranking))
  {
    ss=Sp.risk.ranking[s] 
    fn.overall.risk(N=s,
                    RISK=Overall.risk[[match(ss,names(Overall.risk))]]$risk,
                    sp=capitalize(Sp.risk.ranking[s]),
                    CX=.9)
  }
  axis(1,1.5,"Overall risk",cex.axis=1.75,col.ticks="white",padj=.5)
  
  dev.off()
  
  
  #3.2.Indicator species    
  Sp.risk.ranking=subset(All.Sp.risk.ranking,All.Sp.risk.ranking%in%names(Indicator.species))
  fn.fig(paste(Rar.path,"Risk_Indicator.sp",sep='/'), 2200, 2400)  
  
  par(mar=c(.5,.5,3,1),oma=c(3,9.5,.5,.1),las=1,mgp=c(1,.5,0),cex.axis=1.5,xpd=TRUE)
  layout(matrix(c(rep(1,6),rep(2,3)),ncol=3))
  
  #3.2.1. Risk for each line of evidence
  fn.each.LoE.risk(N.sp=length(Sp.risk.ranking),CX.axis=1.35) #main plot
  for(s in 1:length(Sp.risk.ranking))
  {
    nm=Sp.risk.ranking[s]
    dummy=Overall.risk[[match(nm, names(Overall.risk))]]$Cons.Like.tab
    if(ncol(dummy)==1) Adjst=0
    if(ncol(dummy)>1) Adjst=seq(-.3,.3,length.out = ncol(dummy))
    
    for(nn in 1:ncol(dummy))
    {
      segments(0,s+Adjst[nn],max(dummy[,nn]),s+Adjst[nn],
               lwd=3.75,lend=1,col=LoE.col[match(colnames(dummy)[nn],names(LoE.col))])
    }
    
  }
  axis(2,1:length(Sp.risk.ranking),capitalize(Sp.risk.ranking),cex.axis=1.5)
  mtext("Risk score",1,cex=1.25,line=2)
  legend(-1.5,length(Sp.risk.ranking)+1.3,LoE.names[1:3],
         bty='n',col=LoE.col[1:3],lty=1,lwd=3,horiz = T,cex=1.25,
         text.width=c(0,1,2.25))
  legend(-1.5,length(Sp.risk.ranking)+1.22,LoE.names[4:length(LoE.names)],bty='n',cex=1.5,
         col=LoE.col[4:length(LoE.col)],lty=1,lwd=3,horiz = T,text.width=c(0,1.5))
  
  #3.2.2. Overall risk
  plot(0:1,ylim=c(0,length(Sp.risk.ranking)+1),fg='white',xaxs="i",yaxs="i",
       col="transparent",ylab="",xlab="",xaxt='n',yaxt='n')
  for(s in 1:length(Sp.risk.ranking))
  {
    ss=Sp.risk.ranking[s] 
    fn.overall.risk(N=s,
                    RISK=Overall.risk[[match(ss,names(Overall.risk))]]$risk,
                    sp=capitalize(Sp.risk.ranking[s]),
                    CX=1.5)
  }
  axis(1,1.5,"Overall risk",cex.axis=1.75,col.ticks="white",padj=.5)
  
  dev.off()
  
  
  
}
