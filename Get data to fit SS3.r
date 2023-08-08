# Output files for SS-DL-tool-------------------------------------------------------------------------
THIS=match(c("copper shark","dusky shark","gummy shark","milk shark",
             "sandbar shark","smooth hammerhead","spinner shark",
             "tiger shark","whiskery shark"),Keep.species)
compress.tail=FALSE
for(i in THIS)
{
  Neim=Keep.species[i]
  print(Neim)
  WDD=paste(handl_OneDrive('SS3/SS-DL-tool-master/My data files/'),capitalize(Neim),sep='')
  
  
  #1. Catch
  ktch=KtCh%>%
    filter(Name==Neim)%>%
    mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                         ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                  'TEP_greynurse','TEP_dusky','Discards_TDGDLF'),'Southern.shark',
                                'Other')))%>%
    group_by(SPECIES,Name,finyear,Fishry)%>%
    summarise(Tonnes=sum(LIVEWT.c,na.rm=T))%>%
    mutate(Fishry=case_when(Fishry=="Southern.shark" & finyear<2006 ~'Southern.shark_1',
                            Fishry=="Southern.shark" & finyear>=2006~'Southern.shark_2',
                            TRUE~Fishry))
  Combined.ktch=ktch.combined%>%
    filter(Name==Neim)%>%
    rename(Year=finyear,
           Total=Tonnes)%>%
    ungroup()%>%
    dplyr::select(Year,Total)%>%
    arrange(Year)%>%
    data.frame

  Flits.name=sort(unique(ktch$Fishry))
  Flits=1:length(Flits.name)
  names(Flits)=Flits.name
  Flits.and.survey=data.frame(Fleet.number=c(Flits,1+length(Flits)),
                              Fleet.name=c(names(Flits),"Survey"))
  if(!"Survey"%in% names(Catch.rate.series[[i]])) Flits.and.survey=Flits.and.survey%>%filter(!Fleet.name=="Survey")
  if("Survey"%in% names(Catch.rate.series[[i]]))
  {
    Flits=c(Flits,length(Flits)+1)
    names(Flits)[length(Flits)]="Survey"
  }
  ktch=ktch%>%
    spread(Fishry,Tonnes,fill=0)%>%
    ungroup()%>%
    rename(Year=finyear)%>%
    dplyr::select(- c(SPECIES,Name))
  write.csv(ktch,paste(WDD,'Catches.csv',sep='/'),row.names = F)
  
  #2. Size composition   
  #note: commercial catch and survey. Nsamp set at number of shots
  Size.compo.SS.format=NULL
  if(any(grepl('Size_composition',names(Species.data[[i]]))))
  {
    d.list.n.shots=Species.data[[i]][grep(paste(c("Size_composition_Survey_Observations","Size_composition_Observations",
                                                  "Size_composition_Other_Observations"),collapse="|"),
                                          names(Species.data[[i]]))]
    d.list=Species.data[[i]][grep(paste(c("Size_composition_West","Size_composition_Zone1","Size_composition_Zone2",
                                          "Size_composition_NSF.LONGLINE","Size_composition_Survey",
                                          "Size_composition_Other"),collapse="|"),
                                  names(Species.data[[i]]))]
    if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
    if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
    
    for(s in 1:length(d.list))
    {
      d.list[[s]]=d.list[[s]]%>%
        filter(FL>=List.sp[[i]]$Lzero)%>%
        mutate(fishry=ifelse(grepl("NSF.LONGLINE",names(d.list)[s]),'NSF',
                             ifelse(grepl("Survey",names(d.list)[s]),'Survey',
                                    ifelse(grepl("Other",names(d.list)[s]),'Other',
                                           'TDGDLF'))),
               year=as.numeric(substr(FINYEAR,1,4)),
               TL=FL*List.sp[[i]]$a_FL.to.TL+List.sp[[i]]$b_FL.to.TL,
               size.class=TL.bins.cm*floor(TL/TL.bins.cm))%>%
        filter(TL<=with(List.sp[[i]],max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL))))%>%
        dplyr::rename(sex=SEX)%>%
        filter(!is.na(sex))%>%
        group_by(year,fishry,sex,size.class)%>%
        summarise(n=n())%>%
        ungroup()
    }
    d.list <- d.list[!is.na(d.list)]
    d.list=do.call(rbind,d.list)%>%
      group_by(year,fishry,sex,size.class)%>%
      summarise(n=sum(n))%>%
      ungroup()
    
    size.classes=seq(min(d.list$size.class),max(d.list$size.class),by=TL.bins.cm)
    missing.size.classes=size.classes[which(!size.classes%in%unique(d.list$size.class))]
    if(fill.in.zeros & length(missing.size.classes)>0)
    {
      d.list=rbind(d.list,d.list[1:length(missing.size.classes),]%>%
                     mutate(size.class=missing.size.classes,
                            n=0))
    }
    
    Table.n=d.list%>%group_by(year,fishry,sex)%>%
      summarise(N=sum(n))%>%
      mutate(Min.accepted.N=ifelse(!fishry=='Survey',Min.annual.obs.ktch,10))%>%
      filter(N>=Min.accepted.N)%>%
      mutate(dummy=paste(year,fishry,sex))
    
    if(nrow(Table.n)>0)
    {
      vars.head=c('year','Seas','Fleet','Sex','Part','Nsamp')
      
      if(compress.tail)
      {
        compressed=d.list%>%
          spread(size.class,n,fill=0)
        KolSam=colSums(compressed%>%dplyr::select(-c(year,fishry,sex)))
        
        compressed=d.list%>%
          dplyr::select(size.class,n)%>%
          group_by(size.class)%>%
          summarise(n=sum(n))%>%
          ungroup()%>%
          mutate(Plus.group=size.class[which.min(abs(0.999-cumsum(n)/sum(n)))],
                 Plus.group=ifelse(size.class<Plus.group,size.class,Plus.group))
        
        d.list=d.list%>%
          left_join(compressed%>%dplyr::select(-n),
                    by='size.class')%>%
          group_by(year,fishry,sex,Plus.group)%>%
          summarise(n=sum(n))%>%
          rename(size.class=Plus.group)%>%
          ungroup()
      }
      
      d.list=d.list%>%
        spread(size.class,n,fill=0)%>%
        mutate(month=1)
      
      for(pai in 1:length(d.list.n.shots))
      {
        dd=d.list.n.shots[[pai]]
        iii=length(grep('Survey',names(d.list.n.shots)[pai]))
        if(iii>0)dd$Method='LL'
        if(!'zone'%in%names(dd)) dd$zone=NA
        dd=dd%>%
          mutate(fishry=ifelse(Method=='GN','TDGDLF',
                               ifelse(Method=='LL' & zone!='SA','NSF',
                                      ifelse(zone=='SA','Other',
                                             NA))),
                 year=as.numeric(substr(FINYEAR,1,4)))%>%
          group_by(year,fishry)%>%
          summarise(Nsamp=sum(N.shots))
        if(iii>0) dd$fishry='Survey'
        d.list.n.shots[[pai]]=dd
        rm(dd)
      }
      d.list.n.shots=do.call(rbind,d.list.n.shots)
      
      d.list=left_join(d.list,d.list.n.shots,by=c('year','fishry'))%>%
        mutate(Part=0)%>%
        dplyr::rename(Fleet=fishry,
                      Sex=sex,
                      Seas=month)%>%
        relocate(all_of(vars.head))
      
      d.list=d.list%>%
        mutate(dummy2=paste(year,Fleet,Sex))%>%
        filter(dummy2%in%unique(Table.n$dummy))%>%
        dplyr::select(-dummy2)%>%
        mutate(Sex=ifelse(Sex=='F',1,
                          ifelse(Sex=='M',2,
                                 0)))
      
      d.list=d.list%>%
        mutate(dummy.fleet=case_when(Fleet=="NSF"~'Northern.shark',
                                     Fleet=="Other"~'Other',
                                     Fleet=="TDGDLF" & year<2006~'Southern.shark_1',
                                     Fleet=="TDGDLF" & year>=2006~'Southern.shark_2',
                                     Fleet=="Survey"~'Survey'))%>%
        left_join(Flits.and.survey,by=c('dummy.fleet'='Fleet.name'))%>%
        mutate(Fleet=Fleet.number)%>%
        dplyr::select(-c(dummy.fleet,Fleet.number))%>%
        arrange(Sex,Fleet,year)
      
      d.list.f=d.list%>%filter(Sex==1)%>%arrange(year)  
      d.list.m=d.list%>%filter(Sex==2)%>%arrange(year)
      
      id.var.nms=which(names(d.list.f)%in%vars.head)
      id.var.nms.f=which(!names(d.list.f)%in%vars.head)
      id.var.nms.m=length(id.var.nms.f)+which(!names(d.list.m)%in%vars.head)
      
      
      dummy.Size.compo.SS.format=rbind(d.list.f,d.list.m)
      
      dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
        filter(year<=max(ktch$Year) & Nsamp>=Min.Nsamp)  
      if(nrow(dummy.Size.compo.SS.format)>0) Size.compo.SS.format=dummy.Size.compo.SS.format
      
      names(Size.compo.SS.format)[c(1:2,6)]=c('Year','Month','Nsamps')
      Size.compo.SS.format=Size.compo.SS.format%>%dplyr::select(-Part)
      Size.compo.SS.format=Size.compo.SS.format%>%filter(!is.na(Fleet))
      
      write.csv(Size.compo.SS.format,paste(WDD,'Lengths.csv',sep='/'),row.names = F)
      
    }
  }
  
  
  #4. Abundance series
  Abundance.SS.format=NULL
  CPUE=compact(Catch.rate.series[[i]])
  if(!is.null(CPUE))
  {
    DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
    if(length(DROP)>0)CPUE=CPUE[-DROP]
    if(Neim%in%survey_not.representative & "Survey"%in%names(CPUE)) CPUE=CPUE[-grep("Survey",names(CPUE))]
    if(Neim%in%NSF_not.representative & "NSF"%in%names(CPUE)) CPUE=CPUE[-grep("NSF",names(CPUE))]
    if(Neim%in%tdgdlf_not.representative & "TDGDLF"%in%names(CPUE)) CPUE=CPUE[-grep("TDGDLF",names(CPUE))]

    #reset very low CVs
    #note: Andre suggested leaving original CVs and estimating extraSD if more than one index available
    #      If only 1 index available, then do not estimate, just increase CV before fitting model
    if(length(CPUE)==1)
    {
      for(j in 1:length(CPUE))
      {
        dd=CPUE[[j]]%>%filter(!is.na(Mean))
        reset.CVs=FALSE
        if(reset.CVs)
        {
          loes.mod=loess(log(Mean)~yr.f,data=dd)  
          loes.pred=predict(loes.mod)
          if(CV.use=='loess') use.this.CV=sqrt(sum((log(dd$Mean)-loes.pred)^2)/(length(loes.pred)-2))
          if(CV.use=='fixed') use.this.CV=default.CV
          if(use.this.CV<default.CV) use.this.CV=default.CV
          CPUE[[j]]=CPUE[[j]]%>%
            mutate(CV=ifelse(CV<default.CV,use.this.CV,CV))
        }
        CPUE[[j]]=CPUE[[j]]%>%
          mutate(Label=names(CPUE)[j])
        
        rm(dd)
      }
    }
    
  }
  MAX.CV=List.sp[[i]]$MAX.CV
  for(x in 1:length(CPUE))    
  {
    nm=names(CPUE)[x]
    if(nm=="NSF") nm="Northern.shark"
    if(nm=="TDGDLF.monthly") nm="Southern.shark_1"
    if(nm=="TDGDLF.daily") nm="Southern.shark_2"
    
    dd=CPUE[[x]][,grep(paste(c('yr.f','Mean','MeAn','CV','Label'),collapse="|"),names(CPUE[[x]]))]%>%
      relocate(yr.f)
    if(drop.large.CVs)
    {
      iid=which(dd$CV>MAX.CV)
      dd$Mean[iid]=NA
      dd$CV[iid]=NA 
    }
    dd=dd%>%
      dplyr::rename(Year=yr.f)%>%
      mutate(seas=1,
             index.dummy=nm)%>%
      left_join(Flits.and.survey,by=c('index.dummy'='Fleet.name'))%>%
      mutate(index=Fleet.number)%>%
      dplyr::select(-c(index.dummy,Fleet.number))
    

    CPUE[[x]]=dd%>%filter(!is.na(Mean))
  }
  if(!is.null(CPUE))
  {
    Abundance.SS.format=do.call(rbind,CPUE)%>%
      relocate(Year,seas,index,Mean,CV)%>%
      arrange(index,Year)
    Abundance.SS.format=Abundance.SS.format%>%
      rename(Month=seas,
             Fleet=index,
             Index=Mean)
    write.csv(Abundance.SS.format,paste(WDD,'Index.csv',sep='/'),row.names = F)
  }
  
  
  
  #Life history
  Life.history=List.sp[[i]]
  Life.history$Fecundity=ceiling(mean(Life.history$Fecundity))
  Life.history$Max.age.F=ceiling(mean(Life.history$Max.age.F))
  Life.history$Breed.cycle=mean(Life.history$Breed.cycle)
  LH=data.frame(Parameter=c('M','h','Linf','k','TL.50','TL.95','AwT','BwT','Max.age'),
                Value=with(Life.history,c(Sens.test$SS$Mmean[1],Sens.test$SS$Steepness[1],
                                          Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,Growth.F$k,
                                          TL.50.mat,TL.95.mat,AwT,BwT,Max.age.F)))
  
  
  
  write.csv(LH,paste(WDD,'Life.history.csv',sep='/'),row.names = F)
  
  
  
}

# Run SS by species to improve fit-------------------------------------------------------------------------

tweak.CV=TRUE
Keep.species[THIS]
for(i in THIS)
{
  Neim=Keep.species[i]
  print(Neim)
  
  this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
  if(!dir.exists(this.wd))dir.create(this.wd)
  
  Life.history=List.sp[[i]]
  
  #1. Catch
  ktch=KtCh%>%
    filter(Name==Neim)%>%
    mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                         ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                  'TEP_greynurse','TEP_dusky','Discards_TDGDLF'),'Southern.shark',
                                'Other')))%>%
    group_by(SPECIES,Name,finyear,Fishry)%>%
    summarise(Tonnes=sum(LIVEWT.c,na.rm=T))%>%
    mutate(Fishry=case_when(Fishry=="Southern.shark" & finyear<2006 ~'Southern.shark_1',
                            Fishry=="Southern.shark" & finyear>=2006~'Southern.shark_2',
                            TRUE~Fishry))
  Combined.ktch=ktch.combined%>%
    filter(Name==Neim)%>%
    rename(Year=finyear,
           Total=Tonnes)%>%
    ungroup()%>%
    dplyr::select(Year,Total)%>%
    arrange(Year)%>%
    data.frame
  
  Flits.name=sort(unique(ktch$Fishry))  
  Flits=1:length(Flits.name)
  names(Flits)=Flits.name
  Flits.and.survey=data.frame(Fleet.number=c(Flits,1+length(Flits)),
                              Fleet.name=c(names(Flits),"Survey"))
  if(!"Survey"%in% names(Catch.rate.series[[i]])) Flits.and.survey=Flits.and.survey%>%filter(!Fleet.name=="Survey")
  if("Survey"%in% names(Catch.rate.series[[i]]))
  {
    Flits=c(Flits,length(Flits)+1)
    names(Flits)[length(Flits)]="Survey"
  }
  ktch=ktch%>%
    spread(Fishry,Tonnes,fill=0)
  
  
  #2. Size composition   
  #note: commercial catch and survey. Nsamp set at number of shots
  Size.compo.SS.format=NULL
  if(any(grepl('Size_composition',names(Species.data[[i]]))))
  {
    d.list.n.shots=Species.data[[i]][grep(paste(c("Size_composition_Survey_Observations","Size_composition_Observations",
                                                  "Size_composition_Other_Observations"),collapse="|"),
                                          names(Species.data[[i]]))]
    d.list=Species.data[[i]][grep(paste(c("Size_composition_West","Size_composition_Zone1","Size_composition_Zone2",
                                          "Size_composition_NSF.LONGLINE","Size_composition_Survey",
                                          "Size_composition_Other"),collapse="|"),
                                  names(Species.data[[i]]))]
    if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
    if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
    
    for(s in 1:length(d.list))
    {
      d.list[[s]]=d.list[[s]]%>%
        filter(FL>=List.sp[[i]]$Lzero)%>%
        mutate(fishry=ifelse(grepl("NSF.LONGLINE",names(d.list)[s]),'NSF',
                             ifelse(grepl("Survey",names(d.list)[s]),'Survey',
                                    ifelse(grepl("Other",names(d.list)[s]),'Other',
                                           'TDGDLF'))),
               year=as.numeric(substr(FINYEAR,1,4)),
               TL=FL*List.sp[[i]]$a_FL.to.TL+List.sp[[i]]$b_FL.to.TL,
               size.class=TL.bins.cm*floor(TL/TL.bins.cm))%>%
        filter(TL<=with(List.sp[[i]],max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL))))%>%
        dplyr::rename(sex=SEX)%>%
        filter(!is.na(sex))%>%
        group_by(year,fishry,sex,size.class)%>%
        summarise(n=n())%>%
        ungroup()
    }
    d.list <- d.list[!is.na(d.list)]
    d.list=do.call(rbind,d.list)%>%
      group_by(year,fishry,sex,size.class)%>%
      summarise(n=sum(n))%>%
      ungroup()
    
    MAXX=max(d.list$size.class)
    Tab.si.kl=table(d.list$size.class)
    if(sum(Tab.si.kl[(length(Tab.si.kl)-1):length(Tab.si.kl)])/sum(Tab.si.kl)>.05) MAXX=MAXX*1.2
    MAXX=min(MAXX,30*ceiling(with(List.sp[[i]],max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/30))
    size.classes=seq(min(d.list$size.class),MAXX,by=TL.bins.cm)

    missing.size.classes=size.classes[which(!size.classes%in%unique(d.list$size.class))]
    if(fill.in.zeros & length(missing.size.classes)>0)
    {
      d.list=rbind(d.list,d.list[1:length(missing.size.classes),]%>%
                     mutate(size.class=missing.size.classes,
                            n=0))
    }
    
    Table.n=d.list%>%group_by(year,fishry,sex)%>%
      summarise(N=sum(n))%>%
      mutate(Min.accepted.N=ifelse(!fishry=='Survey',Min.annual.obs.ktch,10))%>%
      filter(N>=Min.accepted.N)%>%
      mutate(dummy=paste(year,fishry,sex))
    
    if(nrow(Table.n)>0)
    {
      vars.head=c('year','Seas','Fleet','Sex','Part','Nsamp')
      
      
      if(Life.history$compress.tail)
      {
        compressed=d.list%>%
          spread(size.class,n,fill=0)
        KolSam=colSums(compressed%>%dplyr::select(-c(year,fishry,sex)))
        
        compressed=d.list%>%
          dplyr::select(size.class,n)%>%
          group_by(size.class)%>%
          summarise(n=sum(n))%>%
          ungroup()%>%
          mutate(Plus.group=size.class[which.min(abs(0.99-cumsum(n)/sum(n)))],
                 Plus.group=ifelse(size.class<Plus.group,size.class,Plus.group))
        
        d.list=d.list%>%
          left_join(compressed%>%dplyr::select(-n),
                    by='size.class')%>%
          group_by(year,fishry,sex,Plus.group)%>%
          summarise(n=sum(n))%>%
          rename(size.class=Plus.group)%>%
          ungroup()
      }
      
      d.list=d.list%>%
        spread(size.class,n,fill=0)%>%
        mutate(month=1)
      
      for(pai in 1:length(d.list.n.shots))
      {
        dd=d.list.n.shots[[pai]]
        iii=length(grep('Survey',names(d.list.n.shots)[pai]))
        if(iii>0)dd$Method='LL'
        if(!'zone'%in%names(dd)) dd$zone=NA
        dd=dd%>%
          mutate(fishry=ifelse(Method=='GN','TDGDLF',
                               ifelse(Method=='LL' & zone!='SA','NSF',
                                      ifelse(zone=='SA','Other',
                                             NA))),
                 year=as.numeric(substr(FINYEAR,1,4)))%>%
          group_by(year,fishry)%>%
          summarise(Nsamp=sum(N.shots))
        if(iii>0) dd$fishry='Survey'
        d.list.n.shots[[pai]]=dd
        rm(dd)
      }
      d.list.n.shots=do.call(rbind,d.list.n.shots)
      
      d.list=left_join(d.list,d.list.n.shots,by=c('year','fishry'))%>%
        mutate(Part=0)%>%
        dplyr::rename(Fleet=fishry,
                      Sex=sex,
                      Seas=month)%>%
        relocate(all_of(vars.head))
      
      d.list=d.list%>%
        mutate(dummy2=paste(year,Fleet,Sex))%>%
        filter(dummy2%in%unique(Table.n$dummy))%>%
        dplyr::select(-dummy2)%>%
        mutate(Sex=ifelse(Sex=='F',1,
                          ifelse(Sex=='M',2,
                                 0)))
      
      d.list=d.list%>%
        mutate(dummy.fleet=case_when(Fleet=="NSF"~'Northern.shark',
                                     Fleet=="Other"~'Other',
                                     Fleet=="TDGDLF" & year<2006~'Southern.shark_1',
                                     Fleet=="TDGDLF" & year>=2006~'Southern.shark_2',
                                     Fleet=="Survey"~'Survey'))%>%
        left_join(Flits.and.survey,by=c('dummy.fleet'='Fleet.name'))%>%
        mutate(Fleet=Fleet.number)%>%
        dplyr::select(-c(dummy.fleet,Fleet.number))%>%
        arrange(Sex,Fleet,year)
      
      d.list.f=d.list%>%filter(Sex==1)%>%arrange(year)  
      d.list.m=d.list%>%filter(Sex==2)%>%arrange(year)
      
      id.var.nms=which(names(d.list.f)%in%vars.head)
      id.var.nms.f=which(!names(d.list.f)%in%vars.head)
      id.var.nms.m=length(id.var.nms.f)+which(!names(d.list.m)%in%vars.head)
      
      dummy.zeros.f=d.list.f[,-match(vars.head,names(d.list.f))]
      dummy.zeros.f[,]=0
      dummy.zeros.m=d.list.m[,-match(vars.head,names(d.list.m))]
      dummy.zeros.m[,]=0
      dummy.Size.compo.SS.format=rbind(cbind(d.list.f,dummy.zeros.f),
                                       cbind(d.list.m[,match(vars.head,names(d.list.m))],
                                             dummy.zeros.m,
                                             d.list.m[,-match(vars.head,names(d.list.m))]))
      names(dummy.Size.compo.SS.format)[id.var.nms.f]=paste('f',names(dummy.Size.compo.SS.format)[id.var.nms.f],sep='')
      names(dummy.Size.compo.SS.format)[id.var.nms.m]=paste('m',names(dummy.Size.compo.SS.format)[id.var.nms.m],sep='')
      
      dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
        filter(year<=max(ktch$finyear) & Nsamp>=Min.Nsamp)  
      if(nrow(dummy.Size.compo.SS.format)>0) Size.compo.SS.format=dummy.Size.compo.SS.format
      
    }
  }
  #keep only observations for fleets with more than 1 year of data
  if(Drop.single.year.size.comp)
  {
    if(!is.null(Size.compo.SS.format))
    {
      Fleet.more.one.year.obs=Size.compo.SS.format%>%
        distinct(Fleet,year)%>%
        group_by(Fleet)%>%
        tally()%>%
        filter(n>1)%>%
        pull(Fleet)
      Size.compo.SS.format=Size.compo.SS.format%>%
        filter(Fleet%in%Fleet.more.one.year.obs)
    }
  }
  
  
  #3. meanbodywt
  meanbodywt.SS.format=NULL
  if(any(grepl('annual.mean.size',names(Species.data[[i]]))))
  {
    meanbodywt.SS.format=Species.data[[i]]$annual.mean.size%>%
      mutate(year=as.numeric(substr(Finyear,1,4)),
             month=1,
             Fleet='Southern.shark_2',
             part=0,   #0, combined; 1: discard only; 2: retained only
             type=2)%>%
      filter(year<=max(ktch$finyear))%>%
      dplyr::select(-Finyear)%>%
      left_join(Flits.and.survey,by=c('Fleet'='Fleet.name'))%>%
      mutate(fleet=Fleet.number)%>%
      dplyr::select(-c(Fleet.number,Fleet))%>%
      relocate(year,month,fleet,part,type,mean,CV)
    
    #CVs can be very small so apply loess SE to free up SS fit (Andre Punt advice)
    dd=meanbodywt.SS.format%>%filter(!is.na(mean))
    loes.mod=loess(log(mean)~year,data=dd)  
    loes.pred=predict(loes.mod)
    if(CV.use=='loess') use.this.CV=sqrt(sum((log(dd$mean)-loes.pred)^2)/(length(loes.pred)-2))
    if(CV.use=='fixed') use.this.CV=default.Mean.weight.CV
    if(use.this.CV<default.Mean.weight.CV) use.this.CV=default.Mean.weight.CV
    #if(use.this.CV>0.6) use.this.CV=0.6
    if(tweak.CV)
      {      meanbodywt.SS.format=meanbodywt.SS.format%>%
                  mutate(CV=ifelse(CV<default.Mean.weight.CV,use.this.CV,CV))
      }
    rm(dd)                   
  }
  
  
  #4. Abundance series
  Abundance.SS.format=NULL
  CPUE=compact(Catch.rate.series[[i]])
  if(!is.null(CPUE))
  {
    DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
    if(length(DROP)>0)CPUE=CPUE[-DROP]
    if(Neim%in%survey_not.representative) CPUE=CPUE[-grep("Survey",names(CPUE))]
    if(Neim%in%NSF_not.representative) CPUE=CPUE[-grep("NSF",names(CPUE))]
    if(Neim%in%tdgdlf_not.representative) CPUE=CPUE[-grep("TDGDLF",names(CPUE))]
    
    #reset very low CVs
    #note: Andre suggested leaving original CVs and estimating extraSD if more than one index available
    #      If only 1 index available, then do not estimate, just increase CV before fitting model
    if(length(CPUE)==1)
    {
      for(j in 1:length(CPUE))
      {
        dd=CPUE[[j]]%>%filter(!is.na(Mean))
        loes.mod=loess(log(Mean)~yr.f,data=dd)  
        loes.pred=predict(loes.mod)
        if(CV.use=='loess') use.this.CV=sqrt(sum((log(dd$Mean)-loes.pred)^2)/(length(loes.pred)-2))
        if(CV.use=='fixed') use.this.CV=default.CV
        if(use.this.CV<default.CV) use.this.CV=default.CV
        if(tweak.CV)
        {
          CPUE[[j]]=CPUE[[j]]%>%
            mutate(CV=ifelse(CV<default.CV,use.this.CV,CV))
        }

        rm(dd)
      }
    }
    
  }
  
  
  #5. F from tagging studies on TDGDLF (1994-95 and 2001-03)
  F.SS.format=NULL  
  if(any(grepl('Fishing.mortality.TDGDLF',names(Species.data[[i]]))))
  {
    F.SS.format=Species.data[[i]][[grep('Fishing.mortality.TDGDLF',names(Species.data[[i]]))]]%>%
      mutate(year=as.numeric(substr(Finyear,1,4)),
             month=1,
             Fleet='Southern.shark_1')%>%
      filter(year<=max(ktch$finyear))%>%
      dplyr::select(-Finyear)%>%
      left_join(Flits.and.survey,by=c('Fleet'='Fleet.name'))%>%
      mutate(fleet=Fleet.number)%>%
      dplyr::select(-c(Fleet.number,Fleet))%>%
      relocate(year,month,fleet,Mean,CV)%>%
      arrange(year)
    
  }
  
  
  #6. Conditional age at length
  #note: this is not used as age-length sandbar and dusky is for GN and LL and for all 4 species
  #      observations were collected over multiple year
  Cond.age.len.SS.format=NULL
  do.Cond.age.len.SS.format=FALSE
  if(do.Cond.age.len.SS.format)
  {
    if(any(grepl('age_length',names(Species.data[[i]]))))
    {
      a=List.sp[[i]]$a_FL.to.TL
      b=List.sp[[i]]$b_FL.to.TL
      Cond.age.len.SS.format=Species.data[[i]]$age_length%>%
        mutate(TL=FL*a+b,
               LbinLo=TL.bins.cm*floor(TL/TL.bins.cm),
               LbinHi=TL.bins.cm*floor(TL/TL.bins.cm))%>%
        group_by(year,Sex,LbinLo)%>%
        mutate(Nsamps=n())%>%
        ungroup()%>%
        mutate(month=7,
               Sex=ifelse(Sex=='Male',2,ifelse(Sex=='Female',1,0)),
               part=0,
               ageErr=1,
               Fleet=Flits[match('Southern.shark_1',names(Flits))])%>%
        group_by(year,Sex,LbinLo,Age)%>%
        mutate(N=n())%>%
        ungroup()%>%
        distinct(year,Sex,Age,LbinLo,.keep_all = T)%>%
        dplyr::select(year,month,Fleet,Sex,part,ageErr,LbinLo,LbinHi,Nsamps,Age,N)%>%
        spread(Age,N,fill=0)
      join.vars=c("year","month","Fleet","Sex","part","ageErr","LbinLo","LbinHi","Nsamps")
      fem.dat=Cond.age.len.SS.format%>%filter(Sex==1)
      diz=match(join.vars,names(fem.dat))
      dumii=fem.dat[,-diz]
      dumii[,]=0
      names(dumii)=paste('m',names(dumii),sep='')
      names(fem.dat)[-diz]=paste('f',names(fem.dat)[-diz],sep='')
      fem.dat=cbind(fem.dat,dumii)
      mal.dat=Cond.age.len.SS.format%>%filter(Sex==2)
      dumii=mal.dat[,-diz]
      dumii[,]=0
      names(dumii)=paste('f',names(dumii),sep='')
      names(mal.dat)[-diz]=paste('m',names(mal.dat)[-diz],sep='')
      mal.dat=cbind(mal.dat,dumii)
      Cond.age.len.SS.format=rbind(fem.dat,mal.dat[,names(fem.dat)])%>%
        arrange(year,month,Fleet,LbinLo,Sex)
    }
  }
  
  
  #7. MeanSize at Age obs
  # Not implemented. Wrong SS3 format. May be applicable to gummy and whiskery (only for these species length-@-age data collected from gillnet fishery)
  MeanSize.at.Age.obs.SS.format=NULL
  Mean.Size.at.age.species=NULL #  Mean.Size.at.age.species=c("gummy shark","whiskery shark" )
  if(any(grepl('age_length',names(Species.data[[i]]))) & names(Species.data)[i]%in%Mean.Size.at.age.species)
  {
    Lvls=c(paste('nf',sort(unique(Species.data[[i]]$age_length%>%filter(Sex=='Female')%>%pull(Age))),sep=''),
           paste('nm',sort(unique(Species.data[[i]]$age_length%>%filter(Sex=='Male')%>%pull(Age))),sep=''))
    Numbs=Species.data[[i]]$age_length%>%
      mutate(Sex=ifelse(Sex=='Male','nm',ifelse(Sex=='Female','nf',NA)),
             Sex.age=paste(Sex,Age,sep=''))
    Numbs$Sex.age=factor(Numbs$Sex.age,levels=Lvls)
    Numbs=Numbs%>%
      group_by(Sex.age)%>%
      tally()%>%
      spread(Sex.age,n)
    
    a=List.sp[[i]]$a_FL.to.TL
    b=List.sp[[i]]$b_FL.to.TL
    MeanSize.at.Age.obs.SS.format=Species.data[[i]]$age_length%>%
      mutate(TL=FL*a+b)%>%
      mutate(Sex=ifelse(Sex=='Male','m',ifelse(Sex=='Female','f',NA)),
             Sex.age=paste(Sex,Age,sep=''))
    MeanSize.at.Age.obs.SS.format$Sex.age=factor(MeanSize.at.Age.obs.SS.format$Sex.age,levels=substr(Lvls,2,5))
    MeanSize.at.Age.obs.SS.format=MeanSize.at.Age.obs.SS.format%>%
      group_by(year,Sex.age)%>%
      summarise(Mean=mean(TL,na.rm=T))%>%
      ungroup()%>%
      mutate(month=7,
             Sex=3,
             part=0,
             ageErr=1,
             ignore=99,
             Fleet=Flits[match('Southern.shark_1',names(Flits))])%>%
      dplyr::select(year,month,Fleet,Sex,part,ageErr,ignore,Sex.age,Mean)%>%
      spread(Sex.age,Mean,fill=99.9)
    
    add.misn.age=substr(Lvls,2,5)
    unk.fem.age=substr(add.misn.age[grepl('f',add.misn.age)],2,4)
    unk.mal.age=substr(add.misn.age[grepl('m',add.misn.age)],2,4)
    add.misn.age=unk.fem.age[which(!unk.fem.age%in%unk.mal.age)]
    if(length(add.misn.age)>0)
    {
      add.misn.age=paste('m',add.misn.age,sep='')
      ddmi=data.frame(t(matrix(add.misn.age)))
      colnames(ddmi)=add.misn.age
      ddmi[,]=99.99
      MeanSize.at.Age.obs.SS.format=cbind(MeanSize.at.Age.obs.SS.format,ddmi)
      ddmi[,]=0
      colnames(ddmi)=paste('n',colnames(ddmi),sep='')
      Numbs=cbind(Numbs,ddmi)
    }
    MeanSize.at.Age.obs.SS.format=cbind(MeanSize.at.Age.obs.SS.format,Numbs)
  }
  
  
  #8. Fleet info
  flitinfo=data.frame(fleet=Flits)%>%
    mutate(type=1,
           surveytiming=-1,
           area=1,
           units=1,
           need_catch_mult=0,
           fleetname=names(Flits))%>%
    dplyr::select(-fleet)%>%
    mutate(type=ifelse(fleetname=="Survey",3,type),
           surveytiming=ifelse(fleetname=="Survey",1,surveytiming))
  rownames(flitinfo)=NULL
  
  
  #9. Run scenarios if available abundance index 
  len.cpue=length(CPUE)
  MAX.CV=List.sp[[i]]$MAX.CV
  for(x in 1:len.cpue)    
  {
    nm=names(CPUE)[x]
    if(nm=="NSF") nm="Northern.shark"
    if(nm=="TDGDLF.monthly") nm="Southern.shark_1"
    if(nm=="TDGDLF.daily") nm="Southern.shark_2"
    
    dd=CPUE[[x]][,grep(paste(c('yr.f','Mean','MeAn','CV'),collapse="|"),names(CPUE[[x]]))]%>%
      relocate(yr.f)
    if(drop.large.CVs)
    {
      iid=which(dd$CV>MAX.CV)
      dd$Mean[iid]=NA
      dd$CV[iid]=NA 
    }
    dd=dd%>%
      dplyr::rename(Year=yr.f)%>%
      mutate(seas=1,
             index.dummy=nm)%>%
      left_join(Flits.and.survey,by=c('index.dummy'='Fleet.name'))%>%
      mutate(index=Fleet.number)%>%
      dplyr::select(-c(index.dummy,Fleet.number))
    

    CPUE[[x]]=dd%>%filter(!is.na(Mean))
  }
  if(!is.null(CPUE))
  {
    Abundance.SS.format=do.call(rbind,CPUE)%>%
      relocate(Year,seas,index,Mean,CV)%>%
      arrange(index,Year)
  }
  
  #Scenarios
  Scens=List.sp[[i]]$Sens.test$SS%>%
    mutate(Species=capitalize(Neim))
  Store.sens=vector('list',nrow(Scens))
  names(Store.sens)=Scens$Scenario
  Out.Scens=Scens
  Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=Out.B.Bmsy=
    Out.F.Fmsy=Out.Kobe.probs=store.warnings=store.convergence=vector('list',length(Store.sens))
  
  #Life history
  Life.history$Fecundity=ceiling(mean(Life.history$Fecundity))
  Life.history$Max.age.F=ceiling(mean(Life.history$Max.age.F))
  Life.history$Breed.cycle=mean(Life.history$Breed.cycle)
  
  #Likelihood lambdas
  Lamdas.SS.lambdas=Life.history$SS_lambdas
  
  #rename fleets following SS nomenclature
  names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))]=match(names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))],names(Flits))
  
  s=1
  this.wd1=paste(this.wd,names(Store.sens)[s],sep='/')
  if(!dir.exists(this.wd1))dir.create(this.wd1)
  
  LH=Life.history
  # LH$Growth.F$k=LH$Growth.F$k*.8
  # LH$Growth.F$FL_inf=LH$Growth.F$FL_inf/.9
  
  SCns=Scens[s,]%>%
    mutate(Model='SS')
  #SCns$Mmean=SCns$Mmean*.8
  
  #cpues
  #remove daily years  
  if(!is.na(SCns$Daily.cpues) & 'TDGDLF.daily'%in%names(CPUE))
  {
    rid.of=as.numeric(unlist(str_split(SCns$Daily.cpues, "&")))
    drop.dis=Abundance.SS.format%>%
      mutate(this=grepl('TDGDLF.daily',rownames(Abundance.SS.format)) & Year%in%rid.of)
    Abundance.SS.format=Abundance.SS.format[-which(drop.dis$this),]
  }
  #a. Create input files
  fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SS'),   
               new.path=this.wd1,
               Scenario=SCns,
               Catch=ktch,
               life.history=LH,
               depletion.yr=NULL,
               fleets=names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))],
               fleetinfo=flitinfo,
               abundance=Abundance.SS.format,   
               size.comp=Size.compo.SS.format,
               meanbodywt=meanbodywt.SS.format,
               F.tagging=F.SS.format,
               cond.age.len=Cond.age.len.SS.format,
               MeanSize.at.Age.obs=MeanSize.at.Age.obs.SS.format,
               Lamdas=Lamdas.SS.lambdas)
  
  #b. Run SS3   
  if(Calculate.ramp.years)
  {
    fn.run.SS(where.inputs=this.wd1,
              where.exe=handl_OneDrive('SS3/ss_win.exe'),
              args='')
    Report=SS_output(this.wd1)
    ramp_years=SS_fitbiasramp(Report)
    out=ramp_years$df
    out=rbind(out,data.frame(value=unique(Report$sigma_R_info$alternative_sigma_R),label='Alternative_sigma_R'))
    write.csv(out,paste(this.wd,'Ramp_years.csv',sep='/'),row.names = F)
  }
  #Arg=''
  Arg='-nohess'
  fn.run.SS(where.inputs=this.wd1,
            where.exe=handl_OneDrive('SS3/ss_win.exe'),
            args=Arg)  
  
  #c. Bring in outputs
  if(Arg=='-nohess')Report=SS_output(this.wd1,covar=F,forecast=F,readwt=F,checkcor=F) else
    Report=SS_output(this.wd1)
  #these.plots=c(1,2,3,6,10,11,12,16,26)
  these.plots=1:26
  SS_plots(Report, plot=these.plots, png=T)
  
  #Check if sigmaR in ball park
  Report$sigma_R_in
  Report$sigma_R_info$alternative_sigma_R
  
  #Posterior vs prior
  if(Report$parameters%>%filter(Label=='SR_BH_steep')%>%pull(Phase)>0)
  {
    pdf(paste(this.wd,paste('Steepness prior vs posterior_',names(Store.sens)[s],'.pdf',sep=''),sep='/'))
    fn.compare.prior.post(d=Report$parameters%>%filter(Label=='SR_BH_steep'),
                          Par='Steepness',
                          prior_type='beta')
    dev.off()
  }
  

}

