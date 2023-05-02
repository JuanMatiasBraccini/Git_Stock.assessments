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
                    by='FishCubeCode')%>%
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
    labs(fill = "")
  p
  ggsave(paste(Exprt,'Figure 2_PSA.tiff',sep='/'), width = W,height = H, dpi = 300, compression = "lzw")
  
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
                 "/1_Inputs/Visualise data/Total catch vs cpues.tiff",sep=''),2400,2000) 
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
                  graph=graph,
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
                    Bf.low,Bf.hi,outfile,CMSY=CMSY.method,nsims=Niters,Proc.error,RES)
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
                     PsiDist,Psiprior,Bprior,BMSYK,output.dir,outfile,Sims,
                     Proc.error.JABBA,Obs.Error=NULL,
                     thinning = 5,nchains = 2,burn.in=5000,ktch.error="random")
{
  # Compile JABBA JAGS model
  if(is.null(CPUE))
  {
    jbinput = build_jabba(catch=Ktch,
                          model.type = MDL,
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
  
  
  return(output)
}
F.from.U=function(U) -log(1-U) 

# Functions for SS3  ------------------------------------------------------
# Create SS input files
fn.get.in.betwee=function(x,PATRN="_") str_before_nth(str_after_nth(x, PATRN, 1), PATRN, 2)
fn.set.up.SS=function(Templates,new.path,Scenario,Catch,life.history,depletion.yr,
                      fleets=NULL,fleetinfo=NULL,abundance=NULL,size.comp=NULL,meanbodywt=NULL,
                      F.tagging=NULL,cond.age.len=NULL,MeanSize.at.Age.obs=NULL,Lamdas=NULL)
{
  # 1.Copy templates
  copy_SS_inputs(dir.old = Templates, dir.new = new.path,overwrite = TRUE)
  
  
  # 2.Read in templates 
  start <- r4ss::SS_readstarter(file = file.path(new.path, "starter.ss"), verbose = FALSE)
  dat <- r4ss::SS_readdat(file = file.path(new.path, start$datfile), verbose = FALSE)
  ctl <- r4ss::SS_readctl(file = file.path(new.path, start$ctlfile), verbose = FALSE, use_datlist = TRUE, datlist = dat)
  fore <- r4ss::SS_readforecast(file = file.path(new.path, "forecast.ss"),  verbose = FALSE)
  
  
  # 3.Update template with species-specific information
  
  #3.1. dat file
  dat$Comments[3]=paste('#C file write time:',Sys.time())
  dat$spawn_month=1.0  #integer is month (1-12), decimal is fraction of days in month (e.g. 15 March is 3.5)
  
  #general age info
  dat$Nages=max(life.history$Max.age.F)
  ageError=as.data.frame(matrix(nrow=2,ncol=dat$Nages+1))
  ageError[1,]=-1.00
  ageError[2,]=0.001
  names(ageError)=seq(0,dat$Nages)
  dat$ageerror=ageError
  
  #population size classes
  dat$binwidth=TL.bins.cm
  dat$minimum_size=10*floor(with(life.history,Lzero*a_FL.to.TL+b_FL.to.TL)/10)
  dat$maximum_size=30*ceiling(with(life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/30)
  
  
  styr=min(Catch$finyear)
  endyr=max(Catch$finyear)
  dat$styr=styr
  dat$endyr=endyr
  
  if(Scenario$Model=='SS')
  {
    #fleets & catch
    dis.flits=fleetinfo$fleetname
    get.fleet.ktch=vector('list',length(fleets))
    for(f in 1:length(get.fleet.ktch))
    {
      xx=Catch[,c('finyear',fleets[f])]
      names(xx)[2]='catch'
      xx=xx%>%
           rename(year=finyear)%>%
           mutate(seas=1,
                  fleet=fleets[f],
                  catch_se=0.01)
      if(f==1)
      {
        dummy=rbind(xx[1,]%>%mutate(year=-999,catch=0),xx)
      }else
      {
        dummy=xx
      }
      get.fleet.ktch[[f]]=dummy
      rm(dummy)
    }
    get.fleet.ktch=do.call(rbind,get.fleet.ktch)%>%
                      relocate(year,seas,fleet,catch,catch_se)%>%
                      data.frame
    dat$catch=get.fleet.ktch
    dat$Nfleets=nrow(fleetinfo)
    dat$fleetinfo=fleetinfo  

    #cpue
    if(!is.null(abundance))   
    {
      ddumy=dat$CPUEinfo%>%
          rownames_to_column('fleetname')%>%
          mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                           ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                           fleetname)))%>%
          filter(fleetname%in%dis.flits)%>%
          mutate(Fleet=row_number())
      rownames(ddumy)=ddumy$fleetname
      dat$CPUEinfo=ddumy%>%dplyr::select(-fleetname)
      dat$CPUE=abundance%>%mutate(Mean=ifelse(Mean<1e-4,1e-4,Mean))
    } 
    
    #Fishing mortality from tagging  
    if(!is.null(F.tagging))
    {
      #notes: to use an F series need to add an additional fleet, cpue series, q and mirror selectivity
      #add fleet
      F.fleet=paste('F.series_',dis.flits[unique(F.tagging$fleet)],sep='')
      dat$Nfleets=dat$Nfleets+length(F.fleet)
      F.fleet.number=dat$Nfleets
      F.catch=dat$catch[1:nrow(F.tagging),]%>%
                mutate(year=F.tagging$year,
                       catch=0,
                       fleet=F.fleet.number,
                       catch_se=0.2)
      dat$catch=rbind(dat$catch%>%filter(!year==-9999),F.catch)
      F.fleetinfo=dat$fleetinfo[1:length(F.fleet.number),]%>%
                    mutate(type=1,
                           fleetname=F.fleet)
      dat$fleetinfo=rbind(dat$fleetinfo,F.fleetinfo)  
      
      #add F series as a cpue series
      Finfo=dat$CPUEinfo[1:length(unique(F.tagging$fleet)),]%>%
                mutate(Fleet=F.fleet.number,
                       Units=2,
                       Errtype=-1)
      rownames(Finfo)=F.fleet
      dat$CPUEinfo=rbind(dat$CPUEinfo,Finfo)
      
      Fcpue=F.tagging%>%
              rename(Year=year,
                     seas=month,
                     index=fleet)%>%
              mutate(index=dat$Nfleets)
      rownames(Fcpue)=paste('F.series',1:nrow(Fcpue),sep='')
      dat$CPUE=rbind(dat$CPUE,Fcpue)
    }

    #meanbodywt
    if(is.null(meanbodywt))
    {
      dat$use_meanbodywt=0
      dat$DF_for_meanbodywt=''
      dat=within(dat, rm(meanbodywt)) 
    }
    if(!is.null(meanbodywt))
    {
      dat$use_meanbodywt=1
      dat$DF_for_meanbodywt=nrow(meanbodywt)-1
      dat$meanbodywt=meanbodywt
    }
    
    #size composition
    if(is.null(size.comp))
    {
      dat$use_lencomp=0
      dat=within(dat, rm(len_info,N_lbins,lbin_vector,lencomp)) 
    }
    if(!is.null(size.comp))   
    {
      ddumy=dat$len_info%>%
        rownames_to_column('fleetname')%>%
        mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                                ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                                       fleetname)))%>%
        filter(fleetname%in%dis.flits)
      rownames(ddumy)=ddumy$fleetname
      dat$len_info=ddumy%>%dplyr::select(-fleetname)
      lbin_vector=sort(as.numeric(gsub('f', '', names(size.comp)[grep("f",names(size.comp))])))
      dat$lbin_vector=lbin_vector
      dat$N_lbins=length(dat$lbin_vector)
      dat$lencomp=size.comp%>%arrange(Sex,Fleet,year)
    }
    if(!is.null(F.tagging))
    {
      addfbit=dat$len_info[length(unique(F.tagging$fleet)),]
      rownames(addfbit)=rownames(Finfo)
      dat$len_info=rbind(dat$len_info,addfbit)
    }
    
    #conditional age at length  
    ddumy=dat$age_info%>%
            rownames_to_column('fleetname')%>%
            mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                                    ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                                           fleetname)))%>%
            filter(fleetname%in%dis.flits)
    rownames(ddumy)=ddumy$fleetname
    dat$age_info=ddumy%>%dplyr::select(-fleetname)
    
    if(!is.null(F.tagging))
    {
      addfbit=dat$age_info[length(unique(F.tagging$fleet)),]
      rownames(addfbit)=rownames(Finfo)
      dat$age_info=rbind(dat$age_info,addfbit)
    }
    if(is.null(cond.age.len))
    {
      dat$agebin_vector=seq(0,dat$Nages-1)
    }
    if(!is.null(cond.age.len))  
    {
      dat$agecomp=cond.age.len
      agebin_vector=sort(as.numeric(gsub('f', '', names(cond.age.len)[grep("f",names(cond.age.len))])))
      dat$agebin_vector=agebin_vector
    }
    
    
    #MeanSize at Age obs   
    if(!is.null(MeanSize.at.Age.obs))
    {
      dat$use_MeanSize_at_Age_obs=1
      dat$MeanSize_at_Age_obs=MeanSize.at.Age.obs
      agebin_vector=suppressWarnings(sort(as.numeric(gsub('f', '', names(MeanSize.at.Age.obs)[grep("f",names(MeanSize.at.Age.obs))]))))
      dat$agebin_vector=agebin_vector
    }
    
    dat$N_agebins=length(dat$agebin_vector)
    
  }
  if(Scenario$Model=='SSS')
  {
    #catch
    dat$catch=data.frame(year=c(-999,Catch$finyear),
                         seas=1,
                         fleet=1,
                         catch=c(0,Catch$Tonnes),
                         catch_se=0.01)
    #dummy cpue
    dat$CPUE=dat$CPUE%>%
      mutate(year=c(styr,endyr),
             obs=c(Scenario$Initial.dpl,Scenario$Final.dpl))  #specify depletion. See page 50 SS manual
    if(!depletion.yr==endyr) dat$CPUE$year[2]=depletion.yr
  }

  
  #3.2. ctl file   
  
  ctl$Comments[2]=dat$Comments[3]
  
  #maturity & fecundity pars
  ctl$First_Mature_Age=0   # Set to 0 and leave Mat ogive take control
  #ctl$First_Mature_Age=life.history$First_Mature_Age 
  fec.option=4   #options: (1)eggs=Wt*(a+b*Wt); (2)eggs=a*L^b; (3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
  if(life.history$Fecu_type_SS=="2_exponential")  fec.option=2
  ctl$fecundity_option=fec.option
  fec.alpha=life.history$Fecu_a    #intercept
  fec.beta=life.history$Fecu_b     #slope
  if(is.na(fec.alpha)| is.na(fec.beta))
  {
    fec.alpha=mean(life.history$Fecundity)  #set fec-@-age to mean fec if no relationship available
    fec.beta=0    
  }
    #divide fecundity by cycle (Spiny dogfish assessment page 46)
  fec.alpha=fec.alpha/mean(life.history$Breed.cycle)  
  fec.beta=fec.beta/mean(life.history$Breed.cycle)
  
  #growth pars
  ctl$Growth_Age_for_L1=0
    #females
  ctl$MG_parms["NatM_p_1_Fem_GP_1", c("INIT","PRIOR")]=rep(Scenario$Mmean,2)
  ctl$MG_parms["L_at_Amin_Fem_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Lzero*a_FL.to.TL+b_FL.to.TL),2)  #size for specified _Age(post-settlement)_for_L1
  ctl$MG_parms["L_at_Amin_Fem_GP_1", "HI"]=max(c(ctl$MG_parms["L_at_Amin_Fem_GP_1", "HI"],ctl$MG_parms["L_at_Amin_Fem_GP_1", "INIT"]*1.5))
  ctl$MG_parms["L_at_Amax_Fem_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL),2) #life.history$TLmax  #set at Linf as _Growth_Age_for_L2 was set at 999
  ctl$MG_parms["L_at_Amax_Fem_GP_1", "HI"]=max(c(ctl$MG_parms["L_at_Amax_Fem_GP_1", "HI"],ctl$MG_parms["L_at_Amax_Fem_GP_1", "INIT"]*1.5))
  ctl$MG_parms["Wtlen_1_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$AwT,2)
  ctl$MG_parms["Wtlen_2_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$BwT,2)
  ctl$MG_parms["VonBert_K_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.F$k,2)
  ctl$MG_parms["VonBert_K_Fem_GP_1", "HI"]=max(c(ctl$MG_parms["VonBert_K_Fem_GP_1", "HI"],life.history$Growth.F$k*1.3))
  ctl$MG_parms["CV_young_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_young,2)
  ctl$MG_parms["CV_old_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_old,2)
  ctl$MG_parms["Mat50%_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$TL.mat.inf.slope[2],2) #life.history$TL.50.mat
  ctl$MG_parms["Mat_slope_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$TL.mat.inf.slope[1],2)
  ctl$MG_parms["Eggs_alpha_Fem_GP_1", c("LO","INIT","HI","PRIOR")]=c(0,fec.alpha,100,fec.alpha)   
  ctl$MG_parms["Eggs_beta_Fem_GP_1", c("INIT","PRIOR")]=rep(fec.beta,2)
  ctl$MG_parms["FracFemale_GP_1", c("INIT","PRIOR")]=rep(life.history$pup.sx.ratio,2)
  
    #males
  ctl$MG_parms["NatM_p_1_Mal_GP_1", c("INIT","PRIOR")]=rep(Scenario$Mmean,2)
  ctl$MG_parms["L_at_Amin_Mal_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Lzero*a_FL.to.TL+b_FL.to.TL),2)
  ctl$MG_parms["L_at_Amin_Mal_GP_1", "HI"]=ctl$MG_parms["L_at_Amin_Fem_GP_1", "HI"]
  ctl$MG_parms["L_at_Amax_Mal_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Growth.M$FL_inf*a_FL.to.TL+b_FL.to.TL),2)
  ctl$MG_parms["L_at_Amax_Mal_GP_1", "HI"]=max(c(ctl$MG_parms["L_at_Amax_Mal_GP_1", "HI"],ctl$MG_parms["L_at_Amax_Mal_GP_1", "INIT"]*1.5))
  ctl$MG_parms["Wtlen_1_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$AwT.M,2)
  ctl$MG_parms["Wtlen_2_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$BwT.M,2)
  ctl$MG_parms["VonBert_K_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.M$k,2)
  ctl$MG_parms["VonBert_K_Mal_GP_1", "HI"]=max(c(ctl$MG_parms["VonBert_K_Mal_GP_1", "HI"],life.history$Growth.M$k*1.3))
  ctl$MG_parms["CV_young_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_young,2)
  ctl$MG_parms["CV_old_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_old,2)
  
  #recruitment pars
  ctl$SR_function=3 # 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
  ctl$SR_parms["SR_LN(R0)", c('LO','INIT','HI')]=with(Scenario,c(Ln_R0_min,Ln_R0_init,Ln_R0_max))
  ctl$SR_parms["SR_BH_steep", "INIT"]=Scenario$Steepness
  ctl$SR_parms["SR_sigmaR", c('LO','HI','INIT')]=c(.01,1,.2) #Spiny dogfish SS assessment
  
  if(Scenario$Model=='SSS') ctl$do_recdev=0  #do_recdev:  0=none; 1=devvector; 2=simple deviations
  if(Scenario$Model=='SS')
  {
    ctl$do_recdev=2
    ctl$recdev_phase=-3
    ctl$recdev_early_start=0
    ctl$max_bias_adj=0.8
    ctl$min_rec_dev=-2
    ctl$max_rec_dev=abs(ctl$min_rec_dev)
    ctl$MainRdevYrFirst=styr
    ctl$MainRdevYrLast=endyr
    ctl$last_early_yr_nobias_adj=styr
    ctl$first_yr_fullbias_adj=min(abundance$Year)+5
    ctl$last_yr_fullbias_adj=endyr-2
    ctl$first_recent_yr_nobias_adj=endyr
  }
  
  #fishing mortality pars
  if(Scenario$use_F_ballpark)
  {
    ctl$F_ballpark=Scenario$F_ballpark
    ctl$F_ballpark_year=styr
  }
  
  #Q pars  
  if(Scenario$Model=='SS')
  {
    iis=sort(unique(dat$CPUE$index))  
    if(!'Northern.shark'%in%fleetinfo$fleetname)
    {
      ctl$Q_options$fleet=ctl$Q_options$fleet-1
    }
    ctl$Q_options=ctl$Q_options%>%filter(fleet%in%iis) 
    addis=iis[which(!iis%in%ctl$Q_options$fleet)]
    if(length(addis)>0)
    {
      addis=dat$CPUEinfo[addis,]
      addis=addis[!rownames(addis)=='Other',]
     
      ctl.q_opt.add=ctl$Q_options[1:nrow(addis),]%>%
                      mutate(fleet=addis$Fleet)
      rownames(ctl.q_opt.add)=rownames(addis)
      ctl$Q_options=rbind(ctl$Q_options,ctl.q_opt.add)%>%
                      arrange(fleet)
      
      ctl.q_par.add=ctl$Q_parms[1:nrow(addis),]
      rownames(ctl.q_par.add)=paste('LnQ_base_',rownames(addis),'(',addis$Fleet,')',sep='')
      ctl$Q_parms=rbind(ctl$Q_parms,ctl.q_par.add)
    }
    Nms=rownames(ctl$Q_parms)
    Nms=gsub(r"{\s*\([^\)]+\)}","",gsub("^.*?base_","",Nms))
    rownames(ctl$Q_parms)=Nms
    ctl$Q_parms=ctl$Q_parms[rownames(ctl$Q_parms)%in%rownames(ctl$Q_options),]
    these.qs=life.history$Q.inits%>%arrange(Fleet.n)%>%pull(Fleet)
    these.qs=subset(these.qs,these.qs%in%rownames(ctl$Q_options))
    ctl$Q_parms=ctl$Q_parms[match(these.qs,row.names(ctl$Q_parms)),]
    ctl$Q_parms=ctl$Q_parms[which(rownames(ctl$Q_parms)%in%rownames(ctl$Q_options)),]
     
    Q.inits=left_join(data.frame(Fleet=rownames(ctl$Q_parms),Order=1:nrow(ctl$Q_parms)),
                      life.history$Q.inits,by='Fleet')%>%
              arrange(Fleet.n)
      
    ctl$Q_parms[,"INIT"]=Q.inits%>%pull(Q)
    ctl$Q_parms[,"PHASE"]=2
    
    #Add extra SD to Q
    #note: Andre suggested leaving original CVs and estimating extraSD if more than one index available
    #      If only 1 index available, then do not estimate, just increase CV before fitting model
    n.indices=nrow(ctl$Q_options)
    if(n.indices>1)
    {
      ctl$Q_options$extra_se=1
      Row.order=rownames(ctl$Q_parms)
      names(Row.order)=seq(1,by=2,length.out=length(Row.order))
      
      
      Q_parms_estraSD=ctl$Q_parms%>%
        mutate(LO=0,
               HI=1,
               INIT=0.3,
               PRIOR=0.3,
               PHASE=3)
      
      rownames(Q_parms_estraSD)=paste(rownames(Q_parms_estraSD),"_Q_extraSD",sep='')
      Row.order_Q_parms=rownames(Q_parms_estraSD)
      names(Row.order_Q_parms)=seq(2,by=2,length.out=length(Row.order_Q_parms))
      
      ctl$Q_parms=rbind(ctl$Q_parms,Q_parms_estraSD)
      Row.order=c(Row.order,Row.order_Q_parms)
      Row.order=Row.order[order(as.numeric(names(Row.order)))]
      ctl$Q_parms=ctl$Q_parms[match(Row.order,rownames(ctl$Q_parms)),]
    }

  }
  
  #selectivity pars  
  if(Scenario$Model=='SS')
  {
    #1. size_selex   
    ddumy=ctl$size_selex_types%>%
      rownames_to_column('fleetname')%>%
      mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                              ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                                     fleetname)))%>%
      filter(fleetname%in%dis.flits)%>%
      mutate(Fleet=row_number())
    rownames(ddumy)=ddumy$fleetname
    if(!'Northern.shark'%in%rownames(ddumy))
    {
      ddumy=ddumy%>%
              mutate(Pattern=ifelse(fleetname=='Other',24,Pattern),
                     Special=ifelse(fleetname=='Other',0,ifelse(fleetname=='Southern.shark_2',2,Special)))
    }
    ctl$size_selex_types=ddumy%>%dplyr::select(-fleetname,-Fleet)

    Sel.ptrn=ctl$size_selex_types$Pattern
    names(Sel.ptrn)=paste('Fishery',1:length(Sel.ptrn),sep='')
    
    row_nm_size_selex_parms=gsub("\\s*\\([^\\)]+\\)","",str_after_nth(rownames(ctl$size_selex_parms), "_", 3)) 
    row_nm_size_selex_parms=ifelse(row_nm_size_selex_parms=='Southern.shark_monthly','Southern.shark_1',
                                   ifelse(row_nm_size_selex_parms=='Southern.shark_daily','Southern.shark_2',
                                          row_nm_size_selex_parms))
    
    chosen.sel.patrns=ctl$size_selex_parms[which(row_nm_size_selex_parms%in%dis.flits),]
    if(!'Northern.shark'%in%dis.flits)
    {
      chosen.sel.patrns=ctl$size_selex_parms[which(row_nm_size_selex_parms%in%c('Northern.shark',dis.flits)),]
      rownames(chosen.sel.patrns)[grep('Northern.shark',rownames(chosen.sel.patrns))]=
        str_replace(rownames(chosen.sel.patrns)[grep('Northern.shark',rownames(chosen.sel.patrns))], "Northern.shark", "Other")
      
      row_nm_size_selex_parms[row_nm_size_selex_parms=="Northern.shark"]="Other"
    }
    ctl$size_selex_parms=chosen.sel.patrns  
    added.bit=str_before_nth(rownames(ctl$size_selex_parms), "_", 3)
    row_nm_size_selex_parms=subset(row_nm_size_selex_parms,row_nm_size_selex_parms%in%dis.flits)

    #turn off irrelevant sel pars
    dummy.sel.pat=ctl$size_selex_types$Pattern
    names(dummy.sel.pat)=rownames(ctl$size_selex_types)
    names(dummy.sel.pat)=ifelse(names(dummy.sel.pat)=='Southern.shark_1','Southern.shark_monthly',
                         ifelse(names(dummy.sel.pat)=='Southern.shark_2','Southern.shark_daily',
                         names(dummy.sel.pat)))
    for(y in 1:length(dummy.sel.pat)) 
    {
      if(dummy.sel.pat[y]==24)ctl$size_selex_parms[grep('P_5',row.names(ctl$size_selex_parms)),"PHASE"]=-2
    }
    
    rownames(ctl$size_selex_parms)=paste(added.bit,row_nm_size_selex_parms,sep="_")
    
    #turn on Southern.shark_2 if size compo data  
    if(!is.null(size.comp))
    {
      Tab.size.comp.dat=with(dat$lencomp%>%filter(Sex==1),table(Fleet))
      names(Tab.size.comp.dat)=fleetinfo$fleetname[as.numeric(names(Tab.size.comp.dat))]
      nn.Southern.shark_2=subset(Tab.size.comp.dat,names(Tab.size.comp.dat)=="Southern.shark_2")
      if(length(nn.Southern.shark_2)>0)
      {
        if(nn.Southern.shark_2>1)
        {
          ctl$size_selex_types[rownames(ctl$size_selex_types)=="Southern.shark_2",]=ctl$size_selex_types[rownames(ctl$size_selex_types)=="Southern.shark_1",]
          
          add.Southern.shark_2.pars=ctl$size_selex_parms[grepl('Southern.shark_1',rownames(ctl$size_selex_parms)),]
          rownames(add.Southern.shark_2.pars)=str_replace(rownames(add.Southern.shark_2.pars), "k_1", "k_2")
          ctl$size_selex_parms=rbind(ctl$size_selex_parms[!grepl('Survey',rownames(ctl$size_selex_parms)),],
                                     add.Southern.shark_2.pars,
                                     ctl$size_selex_parms[grepl('Survey',rownames(ctl$size_selex_parms)),])
        } 
      }
    }
    
    #allocated species specific values to sel pars
    Mirrored.sels=rownames(ctl$size_selex_types%>%filter(Pattern==15))
    life.history$SS_selectivity=life.history$SS_selectivity%>%
                                  filter(Fleet%in%dis.flits)
    if(length(Mirrored.sels)>0) life.history$SS_selectivity=life.history$SS_selectivity%>%filter(!Fleet%in%Mirrored.sels)
    id.fleets=fn.get.in.betwee(x=rownames(ctl$size_selex_parms))
    pis=unique(id.fleets)
    for(px in 1:length(pis))
    {
      iid=which(id.fleets==pis[px])
      this.par=life.history$SS_selectivity[,match(pis[px],colnames(life.history$SS_selectivity))]
      ctl$size_selex_parms[iid,"INIT"]=this.par  
      
      multiplr=rep(0.1,length(this.par))
      multiplr=ifelse(this.par<0,2,multiplr)
      low.bound=multiplr*ctl$size_selex_parms[iid,"INIT"]
      if(pis[px]=="P_1")
      {
        low.bound=sapply(low.bound, function(x) max(dat$minimum_size*1.25,x))
      }
      ctl$size_selex_parms[iid,"LO"]=low.bound
        
      multiplr=rep(2,length(this.par))
      multiplr=ifelse(this.par<0,-2,multiplr)
      up.bound=multiplr*ctl$size_selex_parms[iid,"INIT"]
      if(pis[px]=="P_1")
      {
        up.bound=sapply(up.bound, function(x) min(dat$maximum_size*.975,x))
      }
      
      ctl$size_selex_parms[iid,"HI"]=up.bound
    }
    
    #set phases for estimable selectivities
    if(is.null(size.comp))ctl$size_selex_parms[,"PHASE"]=-2
    if(!is.null(size.comp))
    {
      flit.size.comp.obs=sort(unique(dat$lencomp$Fleet))
      flit.no.size.comp.obs=fleetinfo%>%filter(!fleetname%in%rownames(dat$len_info)[flit.size.comp.obs])%>%pull(fleetname)
      if(length(Mirrored.sels)>0) flit.no.size.comp.obs=subset(flit.no.size.comp.obs,!flit.no.size.comp.obs%in%Mirrored.sels)
      if(length(flit.no.size.comp.obs)>0)
      {
        for(px in 1:length(flit.no.size.comp.obs))
        {
          iid=grep(flit.no.size.comp.obs[px],rownames(ctl$size_selex_parms))
          ctl$size_selex_parms[iid,]$PHASE=-2
        }
      }
    }

    
    #2. age_selex
    ddumy=ctl$age_selex_types%>%
      rownames_to_column('fleetname')%>%
      mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                       ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                       fleetname)))%>%
      filter(fleetname%in%dis.flits)%>%
      mutate(Fleet=row_number())
    rownames(ddumy)=ddumy$fleetname
    ctl$age_selex_types=ddumy%>%dplyr::select(-fleetname,-Fleet)
    
    
    #3. Fishing mortality from tagging  
    if(!is.null(F.tagging))
    {
      F.size.sel.pat=ctl$size_selex_types[match('Southern.shark_2',rownames(ctl$size_selex_types)),]%>%
                        mutate(Pattern=15,Special=3)
      rownames(F.size.sel.pat)=F.fleet
      ctl$size_selex_types=rbind(ctl$size_selex_types,F.size.sel.pat)
      
      F.age.sel.pat=ctl$age_selex_types[match('Southern.shark_2',rownames(ctl$age_selex_types)),]
      rownames(F.age.sel.pat)=F.fleet
      ctl$age_selex_types=rbind(ctl$age_selex_types,F.age.sel.pat)
    }  
  }
  if(Scenario$Model=='SSS')  #SSS assumes selectivity = maturity
  {
   #  ctl$size_selex_types['Fishery','Pattern']=1
     ctl$size_selex_types['Depl','Pattern']=0 
     ctl$age_selex_types['Depl','Pattern']=10
     
     ctl$size_selex_parms['SizeSel_P_1_Fishery(1)',c('INIT','PRIOR')]=life.history$Logistic.selectivity[1]
     ctl$size_selex_parms['SizeSel_P_2_Fishery(1)',c('INIT','PRIOR')]=life.history$Logistic.selectivity[2]
   }
  
  #set prior to init value 
  ctl$SR_parms[,"PRIOR"]=ctl$SR_parms[,"INIT"]
  ctl$MG_parms[,"PRIOR"]=ctl$MG_parms[,"INIT"]
  ctl$Q_parms[,"PRIOR"]=ctl$Q_parms[,"INIT"]
  ctl$size_selex_parms[,"PRIOR"]=ctl$size_selex_parms[,"INIT"]
  
  # Likelihood components (lambdas)
  # Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
  # 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
  if(Scenario$Model=='SS')
  {
    Avail.dat=dat.code=dis.dat=fliit=NULL
    if(!is.null(abundance))
    {
      nn=unique(dat$CPUE$index)
      fliit=nn
      Avail.dat=rep("CPUE",length(nn))
      dat.code=rep(1,length(nn))
      dis.dat=paste('CPUE_',rownames(ctl$Q_options%>%filter(fleet%in%nn)),sep='')
    }
    if(!is.null(meanbodywt))
    {
      nn=unique(dat$meanbodywt$fleet)
      fliit=c(fliit,nn)
      Avail.dat=c(Avail.dat,rep("meanbodywt",length(nn)))
      dat.code=c(dat.code,rep(3,length(nn)))
      dis.dat=c(dis.dat,paste('meanbodywt_',fleetinfo$fleetname[nn],sep='')  )
    }
    if(!is.null(size.comp))
    {
      nn=unique(dat$lencomp$Fleet)
      fliit=c(fliit,nn)
      Avail.dat=c(Avail.dat,rep("size.comp",length(nn)))
      dat.code=c(dat.code,rep(4,length(nn)))
      dis.dat=c(dis.dat,paste('size.comp_',fleetinfo$fleetname[nn],sep='')  )
    }
    if(!is.null(MeanSize.at.Age.obs))
    {
      nn=unique(MeanSize.at.Age.obs$Fleet)
      fliit=c(fliit,nn)
      Avail.dat=c(Avail.dat,rep("meanSize.at.Age",length(nn)))
      dat.code=c(dat.code,rep(7,length(nn)))
      dis.dat=c(dis.dat,paste('meanSize.at.Age_',fleetinfo$fleetname[nn],sep=''))
    }
    Like_comp=ctl$lambdas[1:length(Avail.dat),]%>%
      mutate(like_comp=dat.code,
             fleet=fliit,
             phase=1,
             value=1,
             sizefreq_method=1)
    rownames(Like_comp)=dis.dat
    if(!is.null(Lamdas))  
    {
      Like_comp=Like_comp%>%
                left_join(Lamdas%>%
                            rename(new.value=value),
                          by=c('like_comp','fleet'))%>%
                mutate(value=ifelse(!is.na(new.value),new.value,value))%>%
                dplyr::select(-new.value)
    }
    ctl$lambdas=Like_comp
    ctl$N_lambdas=nrow(ctl$lambdas)
  }

  
  #3.3. starter file
  start$datfile='data.dat'
  start$ctlfile='control.ctl'
  
  #3.4. forecast file
  fore$Fcast_selex=0
  #fore$    #consider updating reference points, years, etc   
  
  
  # 4.Export updated templates
  r4ss::SS_writestarter(start, dir = new.path, overwrite = TRUE,verbose = FALSE)
  r4ss::SS_writedat(dat, outfile = file.path(new.path, start$datfile), overwrite = TRUE, verbose = FALSE)
  r4ss::SS_writectl(ctl, outfile = file.path(new.path, start$ctlfile), overwrite = TRUE, verbose = FALSE)
  r4ss::SS_writeforecast(fore, dir = new.path, file = "forecast.ss", overwrite = TRUE, verbose = FALSE)
}

#Run SS
fn.run.SS=function(where.inputs,where.exe,args=FALSE)
{
  wd_orig=getwd()
  setwd(where.inputs)
  if(!isFALSE(args)) system(paste(shQuote(where.exe),args))else
  {
    system(paste(shQuote(where.exe)))
  }
  setwd(wd_orig)
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
mod.average=function(dd,Weights)
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
  
  return(data.frame(year=yrs,median=Median,lower.95=Lower,upper.95=Upper))
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
  B.limit=Lim.prop.bmsy*B.threshold
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
  
  return(list(probs=data.frame(Range=c('<lim','lim–thr',
                                       'thr–tar','>tar'),
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
    }
    if(Type=='F.series')d1=dd$F.series[1:length(Years)]
    if(Type=='B.Bmsy') d1=dd$B.Bmsy[1:length(Years)]
    if(Type=='F.Fmsy') d1=dd$F.Fmsy[1:length(Years)]
    Dat=data.frame(year=Years,
                   median=apply(d1,2,median),
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
    }
    
    if(Type=='F.series')  d1=dd$F.traj[,1:length(Years)]
    if(Type=='B.Bmsy') d1=dd$B.Bmsy[1:length(Years)]
    if(Type=='F.Fmsy') d1=dd$F.Fmsy[,1:length(Years)]
    
    Dat=data.frame(year=as.numeric(Years),
                   median=apply(d1,2,median),
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
                    lci=apply(dd$posteriors$P,2,function(x) quantile(x,probs=0.025,na.rm=T)),
                    uci=apply(dd$posteriors$P,2,function(x) quantile(x,probs=0.975,na.rm=T)))
      Probs=add.probs(DAT=sweep(dd$posteriors$SB,1,dd$posteriors$K,'/'),
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=median(dd$posteriors$SBmsy)/K)
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
    }
    if(Type=='F.series')
    {
      d1=data.frame(mu=apply(dd$posteriors$H,2,median,na.rm=T),
                    lci=apply(dd$posteriors$H,2,function(x) quantile(x,probs=0.025,na.rm=T)),
                    uci=apply(dd$posteriors$H,2,function(x) quantile(x,probs=0.975,na.rm=T)))
    }
    if(Type=='B.Bmsy')
    {
      d1=data.frame(mu=apply(dd$posteriors$BtoBmsy,2,median,na.rm=T),
                    lci=apply(dd$posteriors$BtoBmsy,2,function(x) quantile(x,probs=0.025,na.rm=T)),
                    uci=apply(dd$posteriors$BtoBmsy,2,function(x) quantile(x,probs=0.975,na.rm=T)))
    }
    if(Type=='F.Fmsy')
    {
      d1=data.frame(mu=apply(dd$posteriors$HtoHmsy,2,median,na.rm=T),
                    lci=apply(dd$posteriors$HtoHmsy,2,function(x) quantile(x,probs=0.025,na.rm=T)),
                    uci=apply(dd$posteriors$HtoHmsy,2,function(x) quantile(x,probs=0.975,na.rm=T)))
    }
    
    
    Dat=data.frame(year=as.numeric(Years),
                   median=d1$mu[1:length(Years)],
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
  
  if(mods=='SSS')
  {
    Years=with(d[[1]],startyr:endyr)  
    
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
      SSB=dum[grep(paste(paste("SSB",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
      SSB=SSB%>%
        mutate(year=readr::parse_number(Label),
               Simulation=readr::parse_number(paste(rownames(SSB),'0',sep='')),
               Simulation=as.numeric(substr(Simulation,5,10)))%>%
        dplyr::select(-Label)%>%
        left_join(SSB_Virgin,by='Simulation')%>%
        mutate(Depletion=Value/Value.virgin)
      
      d1=SSB%>%
        group_by(year)%>%
        summarise(mu=median(Depletion,na.rm=T),
                  lci=quantile(Depletion,probs=0.025,na.rm=T),
                  uci=quantile(Depletion,probs=0.975,na.rm=T))
      
      Probs=add.probs(DAT=SSB%>%dplyr::select(Simulation,year,Depletion)%>%spread(year,Depletion)%>%dplyr::select(-Simulation),
                      id.yr=match(as.numeric(substr(Last.yr.ktch,1,4)),Years),
                      B.threshold=median(SSB_MSY$Depletion))
      Probs$probs=Probs$probs%>%mutate(Scenario=scen)
    }
    
    if(Type=='F.series')
    {
      EF=dum[grep(paste(paste("F",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
      d1=EF%>%
        mutate(year=readr::parse_number(Label))%>%
        group_by(year)%>%
        summarise(mu=median(Value,na.rm=T),
                  lci=quantile(Value,probs=0.025,na.rm=T),
                  uci=quantile(Value,probs=0.975,na.rm=T))
    }
    
    if(Type=='B.Bmsy') 
    {
      SSB=dum[grep(paste(paste("SSB",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
      d1=SSB%>%
        mutate(year=readr::parse_number(Label),
               Simulation=readr::parse_number(paste(rownames(SSB),'0',sep='')),
               Simulation=as.numeric(substr(Simulation,5,10)))%>%
        dplyr::select(-Label)%>%
        left_join(SSB_MSY%>%dplyr::select(Value.MSY,Simulation),by='Simulation')%>%
        mutate(B_Bmys=Value/Value.MSY)%>%
        group_by(year)%>%
        summarise(mu=median(B_Bmys,na.rm=T),
                  lci=quantile(B_Bmys,probs=0.025,na.rm=T),
                  uci=quantile(B_Bmys,probs=0.975,na.rm=T))
      
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
                  lci=quantile(F_Fmsy,probs=0.025,na.rm=T),
                  uci=quantile(F_Fmsy,probs=0.975,na.rm=T))
      
    }
    
    Dat=data.frame(year=as.numeric(Years),
                   median=d1$mu[1:length(Years)],
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
  
  if(Type=='Depletion')
  {
    return(list(Dat=Dat,Probs=Probs))
  }else
  {
    return(list(Dat=Dat))
  }
  
}
fn.display.priors=function(d,sp,XLAB,XLIM)
{
  dummy=lapply(d[sp],function(x) rnorm(1e4,x$mean,x$sd))
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
make_plot <- function(da,nfacets,AXST,AXSt,STRs,InrMarg,dropTitl,Hline,addKtch)
{
  p=da%>%
    ggplot(aes(year, median))+
    geom_line(size=1.1)  +
    geom_ribbon(aes(ymin = lower.95, ymax = upper.95), alpha = 0.3,fill='grey60') +
    facet_wrap(~Scenario,ncol=nfacets)+
    theme_PA(axs.T.siz=AXST,axs.t.siz=AXSt,strx.siz=STRs)+
    theme(plot.title =element_text(hjust = 0.5))+
    ylab("")+xlab(XLAB)+   #ylab(YLAB)
    ylim(0,max(da$upper.95))+
    theme(panel.spacing=unit(InrMarg,"lines"),
          plot.margin = unit(c(.5, -.2, 0, 0), "cm"))
  if(any(!is.na(da$upper.50))) p=p+geom_ribbon(aes(ymin = lower.50, ymax = upper.50), alpha = 0.1)
  if(!is.null(Hline)) p=p+geom_hline(yintercept=Hline, size=1.05,alpha=0.35,
                                     color=rep(c('forestgreen','orange','red'),length(unique(da$Scenario))))
  if(addKtch)
  {
    coeff=max(da$Catch)
    da$Ktch.scaled=da$Catch/coeff
    p=p+ 
      geom_line(data=da, aes(x=year, y=Ktch.scaled),size=1.1,color='dodgerblue4',alpha=0.5,linetype ='dashed')+
      scale_y_continuous(sec.axis = sec_axis(~.*coeff, name=""))
  }
}
fn.ribbon=function(Dat,YLAB,XLAB,Titl,Hline,addKtch,nfacets=1,AXST=14,AXSt=12,STRs=14,InrMarg=.25,dropTitl=FALSE)
{
  data2 <- split(Dat, Dat$Scenario)
  p_lst <- lapply(data2, make_plot,nfacets,AXST,AXSt,STRs,InrMarg,dropTitl,Hline,addKtch)
  figure <- ggarrange(plotlist=p_lst,ncol=nfacets,nrow=ceiling(length(p_lst)/nfacets),
                      common.legend = FALSE)
  if(!dropTitl) figure <- annotate_figure(figure,top=text_grob(Titl, size=20))
  
  return(figure) 
}
fn.plot.ktch.only.timeseries=function(d,sp,Type,YLAB,add.50=FALSE,add.sp.nm=FALSE)
{
  mods=names(d)
  store.plots=vector('list',length(mods))
  names(store.plots)=mods
  store.probs=store.plots
  id=match(sp,Keep.species)
  for(m in 1:length(mods))
  {
    if(!is.null(d[[m]]$estimates[[id]]))
    {
      Hline=NULL
      addKtch=FALSE
      if(Type=='F.series') Var=d[[m]]$f.series[[id]]
      if(Type=='B.Bmsy')   Var=d[[m]]$B.Bmsy[[id]]
      if(Type=='F.Fmsy')   Var=d[[m]]$F.Fmsy[[id]]
      
      if(Type=='Depletion')
      {
        addKtch=TRUE
        Var=d[[m]]$rel.biom[[id]]
        str.prob=compact(d[[m]]$probs.rel.biom[[id]])   
        Hline=str.prob[[1]]$Reference.points$Value
        store.probs[[m]]=do.call(rbind,sapply(str.prob,'[',1))%>%  
          mutate(Model=names(store.probs)[m])
      }
      store.plots[[m]]=fn.ribbon(Dat=Var,
                                 YLAB='',
                                 XLAB="",
                                 Titl=names(d)[m],
                                 Hline=Hline,
                                 addKtch=addKtch)
      
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
    if(Type=='Depletion')
    {
      figure=annotate_figure(figure,
                             right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                             color ='dodgerblue4',vjust=0))
    }
    
    print(figure)      
    store.probs=do.call(rbind,store.probs)
    rownames(store.probs)=NULL
    return(list(store.probs=store.probs,Ref.points=str.prob[[1]]$Reference.points))
    
  }
}
fn.plot.ktch.only.timeseries_combined_Appendix=function(this.sp,d,YLAB,NM,Type,InnerMargin)
{
  mods=names(d)
  if(length(this.sp)>8)
  {
    for(m in 1:length(mods))
    {
      id=match(this.sp,names(d[[m]]$rel.biom))
      Hline=NULL
      addKtch=FALSE
      if(Type=='F.series')  Var=d[[m]]$f.series[id]
      if(Type=='B.Bmsy')    Var=d[[m]]$B.Bmsy[id]
      if(Type=='F.Fmsy')    Var=d[[m]]$F.Fmsy[id]
      if(Type=='Depletion')
      {
        Var=d[[m]]$rel.biom[id]
        addKtch=TRUE
        Hline=d[[m]]$probs.rel.biom[[id[1]]][[1]]$Reference.points$Value
      }
      figure=fn.ribbon(Dat=do.call(rbind,Var)%>%
                         rownames_to_column()%>%
                         filter(Scenario=='S1')%>%
                         mutate(Scenario=capitalize(word(rowname,1,sep = "\\."))),
                       YLAB='',
                       XLAB="",
                       Titl=mods[m],
                       Hline=Hline,
                       addKtch=addKtch,
                       nfacets=round(length(this.sp)/5),
                       AXST=15,AXSt=13,STRs=12,
                       InrMarg=InnerMargin,
                       dropTitl=TRUE)
      
      figure=annotate_figure(figure,
                             bottom = text_grob('Financial year',size=26,vjust =-0.15),
                             left = text_grob(YLAB,size=26,rot = 90,vjust=0.8))+
        theme(plot.margin = unit(c(.1,.1,0,0), "cm"))
      if(Type=='Depletion')
      {
        figure=annotate_figure(figure,
                               right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                               color ='dodgerblue4',vjust=-.2))
      }
      print(figure)
      ggsave(paste(Rar.path,'/Relative.biomass_catch.only_',NM,'_',mods[m],'_Appendix','.tiff',sep=''),
             width = 11,height = 11,compression = "lzw")
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
      if(Type=='B.Bmsy')    Var=d[[m]]$B.Bmsy[id]
      if(Type=='F.Fmsy')    Var=d[[m]]$F.Fmsy[id]
      if(Type=='Depletion')
      {
        Var=d[[m]]$rel.biom[id]
        addKtch=TRUE
        Hline=d[[m]]$probs.rel.biom[[id[1]]][[1]]$Reference.points$Value
      }
      plot.list[[m]]=fn.ribbon(Dat=do.call(rbind,Var)%>%
                                 filter(Scenario=='S1')%>%
                                 rownames_to_column()%>%
                                 mutate(Scenario=capitalize(word(rowname,1,sep = "\\."))),
                               YLAB='',
                               XLAB="",
                               Titl=mods[m],
                               Hline=Hline,
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
    if(Type=='Depletion')
    {
      figure=annotate_figure(figure,
                             right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                             color ='dodgerblue4',vjust=-.2))
    }
    print(figure)
    ggsave(paste(Rar.path,'/Relative.biomass_catch.only_',NM,'_Appendix','.tiff',sep=''),
           width = 15,height = 10,compression = "lzw")
    
  }
}
fn.plot.ktch.only.timeseries_combined=function(this.sp,d,YLAB,Type,InnerMargin,RefPoint,Kach)
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
    if(Type=='F.series')  Var=d$f.series[id]
    if(Type=='B.Bmsy')    Var=d$B.Bmsy[id]
    if(Type=='F.Fmsy')    Var=d$F.Fmsy[id]
    if(Type=='Depletion')
    {
      Var=d$rel.biom[id]
      addKtch=TRUE
    }
    Var=compact(Var)
    
    if(length(Var)>8)  Nfast=round(length(Var)/5)
    if(length(Var)>3 & length(Var)<8) Nfast=2
    if(length(Var)<=3) Nfast=1
    
    Dat=do.call(rbind,Var)
    if('Catch'%in%names(Dat))Dat=Dat%>%dplyr::select(-Catch)
    
    figure=fn.ribbon(Dat=Dat%>%
                       rownames_to_column()%>%
                       mutate(Scenario=capitalize(word(rowname,1,sep = "\\.")))%>%
                       left_join(Kach,by=c('year','Scenario')),
                     YLAB='',
                     XLAB="",
                     Titl='',
                     Hline=RefPoint[[1]]$Value,
                     addKtch=addKtch,
                     nfacets=Nfast,
                     AXST=15,AXSt=13,STRs=15,
                     InrMarg=InnerMargin,
                     dropTitl=TRUE)
    
    figure=annotate_figure(figure,
                           bottom = text_grob('Financial year',size=26,vjust =-0.15),
                           left = text_grob(YLAB,size=26,rot = 90,vjust=0.8)) +
      theme(plot.margin = unit(c(.1,.1,0,0), "cm"))
    if(Type=='Depletion')
    {
      figure=annotate_figure(figure,
                             right=text_grob('Total catch (tonnes)',size=26,rot = 90,
                                             color ='dodgerblue4',vjust=0))
    }
    return(figure)
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
    scale_x_continuous(limits=c(0,Mx.B)) +
    scale_y_continuous(limits=c(0,Mx.F))+
    geom_rect(xmin = 1, xmax = Mx.B, ymin = 0, ymax = 1, fill = 'chartreuse3', alpha = 0.05) +
    geom_rect(xmin = 0, xmax = 1, ymin = 1, ymax = Mx.F, fill = 'brown1', alpha = 0.05) +
    geom_rect(xmin = 1, xmax = Mx.B, ymin = 1, ymax = Mx.F, fill = 'orange', alpha = 0.05) +
    geom_rect(xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = 'yellow', alpha = 0.05)
  if(!is.null(Probs))
  {
    kernelF <- gplots::ci2d(Probs$x, Probs$y, nbins = 151, factor = 1.5, 
                            ci.levels = c(0.5, 0.8, 0.75, 0.9, 0.95), show = "none")
    KernelD=rbind(kernelF$contours$"0.95"%>%mutate(CI='1',col='grey30'),
                  kernelF$contours$"0.8"%>%mutate(CI='2',col='grey50'),
                  kernelF$contours$"0.5"%>%mutate(CI='3',col='grey75'))
    kernels=KernelD%>%distinct(CI,col)%>%pull(col)
    names(kernels)=KernelD%>%distinct(CI,col)%>%pull(CI)
    
    Pr.d=data.frame(
      Prob=c(sum(ifelse(Probs$x > 1 & Probs$y < 1, 1, 0))/length(Probs$x)*100,
             sum(ifelse(Probs$x < 1 & Probs$y < 1, 1, 0))/length(Probs$x)*100,
             sum(ifelse(Probs$x > 1 & Probs$y > 1, 1, 0))/length(Probs$x)*100,
             sum(ifelse(Probs$x < 1 & Probs$y > 1, 1, 0))/length(Probs$x) * 100),
      col=c("green","yellow","orange","red"),
      x=rep(-10,4),  #dummy
      y=rep(-10,4))
    pr.ds=Pr.d%>%pull(col)
    names(pr.ds)=paste(round(Pr.d%>%pull(Prob),1),'%',sep='')
    
    
    
    kobe <-kobe +
      geom_polygon(data=KernelD,aes(x, y,fill=CI),size=1.25,alpha=0.5)+
      scale_fill_manual(labels=c("95%","80%","50%"),values = kernels)+
      geom_point(data=Pr.d,aes(x, y,color=col),alpha = 1,size=5)+
      scale_color_manual(labels=names(pr.ds),values = pr.ds)+
      labs(CI="CI", col="Prob.")
    
    
    
  }
  kobe <-kobe +
    geom_path(linetype = 2, size = 0.5,color='steelblue')+
    geom_point(size=2,color='steelblue')+
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
  return(kobe)
}
fn.get.Kobe.plot_appendix=function(d,sp,Scen='S1',add.sp.nm=FALSE,do.probs=FALSE)
{
  id=match(sp,Keep.species)
  
  #DBSRA
  if('DBSRA'%in%names(d))
  {
    dummy=d$DBSRA$B.Bmsy[[id]]%>%filter(Scenario==Scen)
    yrs=dummy%>%pull(year)
    Bmsy=dummy$median
    Fmsy=d$DBSRA$F.Fmsy[[id]]%>%filter(Scenario==Scen)%>%pull(median)
    p.DBSRA=kobePlot(f.traj=Fmsy[1:length(yrs)],
                     b.traj=Bmsy[1:length(yrs)],
                     Years=yrs,
                     Titl="DBSRA")
    rm(yrs,Fmsy,Bmsy,dummy)
  }
  
  #CMSY
  if('CMSY'%in%names(d))
  {
    if(CMSY.method=="Froese")
    {
      yrs=d$CMSY[[id]][[Scen]]$output$ref_ts$year
      Bmsy=d$CMSY[[id]][[Scen]]$output$ref_ts$bbmsy
      Fmsy=d$CMSY[[id]][[Scen]]$output$ref_ts$ffmsy
      
    }
    if(CMSY.method=="Haddon")
    {
      dummy=d$CMSY$B.Bmsy[[id]]%>%filter(Scenario==Scen)
      yrs=dummy%>%pull(year)
      Bmsy=dummy$median
      Fmsy=d$CMSY$F.Fmsy[[id]]%>%filter(Scenario==Scen)%>%pull(median)
    }
    p.CMSY=kobePlot(f.traj=Fmsy,
                    b.traj=Bmsy,
                    Years=yrs,
                    Titl="CMSY")
    rm(yrs,Fmsy,Bmsy,dummy)
  }
  
  #JABBA
  dummy=d$JABBA$B.Bmsy[[id]]%>%filter(Scenario==Scen)
  yrs=dummy%>%pull(year)
  Bmsy=dummy$median
  Fmsy=d$JABBA$F.Fmsy[[id]]%>%filter(Scenario==Scen)%>%pull(median)
  if(!do.probs)
  {
    p.JABBA=kobePlot(f.traj=Fmsy,
                     b.traj=Bmsy,
                     Years=yrs,
                     Titl="JABBA")
  }
  if(do.probs)
  {
    p.JABBA=kobePlot(f.traj=Fmsy,
                     b.traj=Bmsy,
                     Years=yrs,
                     Titl="JABBA",
                     Probs=data.frame(x=d$JABBA$Kobe.probs[[id]]$stock,
                                      y=d$JABBA$Kobe.probs[[id]]$harvest))
  }
  
  #SSS
  if('SSS'%in%names(d))
  {
    dummy=d$SSS$B.Bmsy[[id]]%>%filter(Scenario==Scen)
    yrs=dummy%>%pull(year)
    Bmsy=dummy$median
    Fmsy=d$SSS$F.Fmsy[[id]]%>%filter(Scenario==Scen)%>%pull(median)
    p.SSS=kobePlot(f.traj=Fmsy[1:length(yrs)],
                     b.traj=Bmsy[1:length(yrs)],
                     Years=yrs,
                     Titl="SSS")
    rm(yrs,Fmsy,Bmsy,dummy)
  }
  
  #Combine plots
  if('DBSRA'%in%names(d) & 'CMSY'%in%names(d)& 'SSS'%in%names(d))
  {
    plotlist=list(DBSRA=p.DBSRA+rremove("axis.title"),
                  CMSY=p.CMSY+rremove("axis.title"),
                  JABBA=p.JABBA+rremove("axis.title"),
                  SSS=p.SSS+rremove("axis.title"))
  }else
  {
    plotlist=list(JABBA=p.JABBA+rremove("axis.title"))
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
    if(!do.probs)
    {
      dd=kobePlot(f.traj=Fmsy[[x]]$median,
                  b.traj=Bmsy[[x]]$median,
                  Years=Bmsy[[x]]$year,
                  Titl=capitalize(names(Bmsy)[x]),
                  YrSize=6)+rremove("axis.title")
    }
    
    if(do.probs)
    {
      dd=kobePlot(f.traj=Fmsy[[x]]$median,
                  b.traj=Bmsy[[x]]$median,
                  Years=Bmsy[[x]]$year,
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

# Store consequence-likelihood  ------------------------------------------------------
fn.get.cons.like=function(lista,Mod.Avrg.weight=NULL)
{
  #get likelihood for each model
  out=vector('list',length(lista))
  names(out)=names(lista)
  for(l in 1:length(lista))
  {
    d=lista[[l]]$probs.rel.biom
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
  }
  
  #get weighted likelihood
  out2=vector('list',N.sp)
  names(out2)=names(List.sp)
  if(!is.null(Mod.Avrg.weight))
  {
    for(i in 1:N.sp)  
    {
      dis.r.range=do.call(rbind,r.groups)%>%   
        data.frame%>%
        mutate(Get=between(r.list[i],.[[1]] , .[[2]]))%>%
        filter(Get==TRUE)
      Wei=Mod.Avrg.weight%>%
        filter(r.group==row.names(dis.r.range))%>%
        dplyr::select(-r.group)
      Wei.SSS=Wei%>%filter(Model=='DBSRA')%>%mutate(Model='SSS')
      if('SSS'%in%names(Catch_only)) Wei=rbind(Wei,Wei.SSS)
      dd=map(out, ~.x[[i]])
      dd=do.call(rbind,dd)%>%
        left_join(Wei,by=c('Model'))%>%
        group_by(Range)%>%
        summarise(Probability=weighted.mean(Probability,Weight))%>%
        ungroup()%>%
        mutate(Probability=Probability/sum(Probability))%>%
        mutate(Range=factor(Range,levels=c('>tar','thr–tar','lim–thr','<lim')))%>%
        arrange(Range)
      out2[[i]]=dd
    }
  }
  if(is.null(Mod.Avrg.weight))
  {
    for(i in 1:N.sp)  
    {
      dd=map(out, ~.x[[i]])
      dd=do.call(rbind,dd)
      if(!is.null(dd))out2[[i]]=dd
    }
  }
  return(out2)
}
# Weight of Evidence assessment   ------------------------------------------------------
fn.risk=function(likelihood)
{
  likelihood=likelihood%>%
    mutate(consequence=case_when(Range=='>tar'~'C1',
                                 Range=='thr–tar'~'C2',
                                 Range=='lim–thr'~'C3',
                                 Range=='<lim'~'C4'))
  TAB=Risk.tab
  for(n in 1:nrow(likelihood))
  {
    id=match(likelihood$consequence[n],TAB$Consequence)
    idd=which(unlist(lapply(Like.ranges,function(x) check.in.range(likelihood$Probability[n],x,fatal=F))))
    TAB[id,idd+1]="X"
    TAB$Max.Risk.Score[id]=id*idd
  }
  return(TAB)
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
  polygon(negigible.Vec, x.Vec, col = 'cornflowerblue', border = "transparent")
  polygon(low.Vec, x.Vec, col = 'olivedrab3', border = "transparent")
  polygon(medium.Vec, x.Vec, col = 'yellow', border = "transparent")
  polygon(high.Vec, x.Vec, col = 'orange', border = "transparent")
  polygon(severe.Vec, x.Vec, col = 'red', border = "transparent")
  
  axis(1,at=c(mean(negigible.Vec),mean(low.Vec),mean(medium.Vec),mean(high.Vec),mean(severe.Vec)),
       labels=c("Negl.","Low","Med.","High","Severe"),cex.axis=CX.axis)
  box()
}
fn.overall.risk=function(N,RISK,sp,CX=1.5)
{
  x.Vec <-  c(1,3,3,1)
  y.Vec <- c(rep((s-.5),2),rep((s+.5),2))
  if(RISK=="Negligible") CL = 'cornflowerblue'
  if(RISK=="Low") CL = 'olivedrab3'
  if(RISK=="Medium") CL ='yellow'
  if(RISK=="High") CL ='orange'
  if(RISK=="Severe") CL ='red'
  polygon(x.Vec, y.Vec, col = CL, border = "transparent")
  text(1.5,mean(y.Vec),sp,cex=CX)
}