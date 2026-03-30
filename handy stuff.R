#-----------  manually set up fn.set.up.SS()-------------------------------------------------------------------------
Templates=handl_OneDrive('SS3/Examples/SS')   
new.path=this.wd1
Scenario=Scens[s,]%>%mutate(Model='SS')
Catch=KAtch
Catch.ret.disc=KAtch.ret.disc
life.history=Life.history
depletion.yr=NULL
fleets=names(KAtch)[which(!names(KAtch)%in%c("SPECIES","Name","finyear"))]
fleetinfo=FLitinFO
abundance=Abund   
size.comp=Size.com
meanbodywt=meanbody
Tags=Tags.SS.format.zone
F.tagging=F.SS.format
cond.age.len=Cond.age.len.SS.format
MeanSize.at.Age.obs=MeanSize.at.Age.obs.SS.format
Lamdas=Lamdas.SS.lambdas
Var.adjust.factor=Var.ad
Future.project=add.future

RecDev_Phase=-3
SR_sigmaR=0.2
first.age=0
age.comp=NULL


#-----------  set up Scens-------------------------------------------------------------------------
for(i in 1:N.sp)
{
  Neim=Keep.species[i]
  
  this.wd=paste(HandL.out,capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
  if(!dir.exists(this.wd))dir.create(this.wd)
  
  Life.history=List.sp[[i]]
  
  if(Neim%in%drop.min.pop.bin.size) Life.history$Lzero=Life.history$Lzero/1.046 
  
  #1. Catch
  #1.1. zones together
  ktch=KtCh%>%
    filter(Name==Neim)
  retained.discarded.ktch=NULL
  discard.specs=c('TEP','Discards_TDGDLF')
  if(Neim%in%retained.discarded.sp)
  {
    if(Neim=="dusky shark") this.discard.sp='TEP'
    if(exists('this.discard.sp'))
    {
      retained.discarded.ktch=ktch%>%
                  filter(FishCubeCode==this.discard.sp)%>%
                  mutate(Fishry='Southern.shark_2')%>%
                  group_by(SPECIES,Name,finyear,Fishry)%>%
                  summarise(Tonnes=sum(LIVEWT.c,na.rm=T)) 
      discard.specs=subset(discard.specs,!discard.specs==this.discard.sp)
      ktch=ktch%>%
        filter(!FishCubeCode==this.discard.sp)
      if(retained.discarded.units=='numbers' & Neim=="dusky shark")  #NEW set to units 3 if retained.discarded.units=='numbers'
      {
        retained.discarded.ktch=retained.discarded.ktch%>%
                  left_join(TEPS_dusky_n.discards%>%
                                  mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
                                  group_by(finyear)%>%
                                  summarise(Discards.n_1000s=sum(Discards.n_1000s))%>%
                              ungroup(),
                            by='finyear')%>%
          mutate(Discards.n_1000s=ifelse(is.na(Discards.n_1000s),0,Discards.n_1000s),
                 Tonnes=0,
                 Tonnes=Discards.n_1000s)%>%
          dplyr::select(-Discards.n_1000s)
          
      }
    }
  }
  ktch=ktch%>%
    mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                         ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                  discard.specs),'Southern.shark',
                                ifelse(FishCubeCode%in%c('WRL') & Neim%in%WRL.species,'WRL',
                                       'Other'))))%>%
    group_by(SPECIES,Name,finyear,Fishry)%>%
    summarise(Tonnes=sum(LIVEWT.c,na.rm=T))%>%
    mutate(Fishry=case_when(Fishry=="Southern.shark" & finyear<2006 ~'Southern.shark_1',
                            Fishry=="Southern.shark" & finyear>=2006~'Southern.shark_2',
                            TRUE~Fishry))
  #1.2. by zone
  ktch.zone=KtCh.zone%>%
    ungroup()%>%
    filter(Name==Neim)
  
  #allocated historic Southern.shark to zones
  if('Historic'%in%unique(ktch.zone$FishCubeCode))
  {
    historic=ktch.zone%>%
      filter(FishCubeCode=='Historic')%>%
      dplyr::select(-zone)
    ktch.zone=ktch.zone%>%filter(!FishCubeCode=='Historic')
    First.year.k=ktch.zone%>%filter(FishCubeCode%in%c("JASDGDL","WCDGDL"))
    Dis.Yr=min(First.year.k$finyear)
    if(Neim=='sandbar shark') Dis.Yr=1985
    prop.k=First.year.k%>%
      filter(finyear==Dis.Yr)%>%
      group_by(zone)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c))%>%
      mutate(Prop=LIVEWT.c/sum(LIVEWT.c))%>%
      dplyr::select(-LIVEWT.c)
    dum=expand.grid(finyear=historic%>%distinct(finyear)%>%pull(finyear),
                    zone=unique(prop.k$zone))
    prop.k=full_join(prop.k,dum,by='zone')%>%arrange(finyear)
    historic=full_join(historic,prop.k,by='finyear')%>%
      mutate(LIVEWT.c=LIVEWT.c*Prop)%>%
      dplyr::select(-Prop)%>%
      relocate(names(ktch.zone))
    ktch.zone=rbind(ktch.zone,historic)
    
  }
  retained.discarded.ktch.zone=NULL  
  if(Neim%in%retained.discarded.sp)
  {
    if(Neim=="dusky shark") this.discard.sp='TEP'
    if(exists('this.discard.sp'))
    {
      retained.discarded.ktch.zone=ktch.zone%>%
                        filter(FishCubeCode==this.discard.sp)%>%
                        mutate(Fishry=paste('Southern.shark_2',zone,sep='_'))%>%
                        group_by(SPECIES,Name,finyear,Fishry)%>%
                        summarise(Tonnes=sum(LIVEWT.c,na.rm=T)) 
      ktch.zone=ktch.zone%>%
        filter(!FishCubeCode==this.discard.sp)
      
      if(retained.discarded.units=='numbers' & Neim=="dusky shark")  
      {
        retained.discarded.ktch.zone=retained.discarded.ktch.zone%>%
                              left_join(TEPS_dusky_n.discards%>%
                                          mutate(Fishry=paste('Southern.shark_2',zone,sep='_'),
                                                 finyear=as.numeric(substr(FINYEAR,1,4)))%>%
                                          group_by(finyear,Fishry)%>%
                                          summarise(Discards.n_1000s=sum(Discards.n_1000s))%>%
                                          ungroup(),
                                        by=c('finyear','Fishry'))%>%
                              mutate(Discards.n_1000s=ifelse(is.na(Discards.n_1000s),0,Discards.n_1000s),
                                     Tonnes=0,
                                     Tonnes=Discards.n_1000s)%>%
                              dplyr::select(-Discards.n_1000s)
        
      }
    }
  }
  ktch.zone=ktch.zone%>%
    mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                         ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                  discard.specs),'Southern.shark',
                                ifelse(FishCubeCode%in%c('WRL') & Neim%in%WRL.species,'WRL',
                                       'Other'))),
           Fishry=case_when(Fishry=="Southern.shark" & finyear<2006 ~'Southern.shark_1',
                            Fishry=="Southern.shark" & finyear>=2006~'Southern.shark_2',
                            TRUE~Fishry),
           Fishry=ifelse(grepl('Southern.shark',Fishry),paste(Fishry,zone,sep='_'),Fishry))%>%
    group_by(SPECIES,Name,finyear,Fishry)%>%
    summarise(Tonnes=sum(LIVEWT.c,na.rm=T))
  
  #1.3 combined catch for displaying only
  Combined.ktch=ktch.combined%>%
    filter(Name==Neim)%>%
    rename(Year=finyear,
           Total=Tonnes)%>%
    ungroup()%>%
    dplyr::select(Year,Total)%>%
    arrange(Year)%>%
    data.frame
  
  #fleets
  #zones together
  Flits.name=sort(unique(ktch$Fishry))  
  Flits=1:length(Flits.name)
  names(Flits)=Flits.name
  Flits.and.survey=data.frame(Fleet.number=c(Flits,1+length(Flits)),
                              Fleet.name=c(names(Flits),"Survey"))
  if(!"Survey"%in% names(Catch.rate.series[[i]]) & !"Size_composition_Survey"%in%names(Species.data[[i]]))
  {
    Flits.and.survey=Flits.and.survey%>%filter(!Fleet.name=="Survey")
  }
  if("Survey"%in% names(Catch.rate.series[[i]]) | "Size_composition_Survey"%in%names(Species.data[[i]]))
  {
    Flits=c(Flits,length(Flits)+1)
    names(Flits)[length(Flits)]="Survey"
  }
  ktch=ktch%>%
    spread(Fishry,Tonnes,fill=0)
  
  #by zone
  Flits.name.zone=sort(unique(ktch.zone$Fishry))  
  Flits.zone=1:length(Flits.name.zone)
  names(Flits.zone)=Flits.name.zone
  Flits.and.survey.zone=data.frame(Fleet.number=c(Flits.zone,1+length(Flits.zone)),
                                   Fleet.name=c(names(Flits.zone),"Survey"))
  if(!"Survey"%in% names(Catch.rate.series[[i]]) & !"Size_composition_Survey"%in%names(Species.data[[i]]))
  {
    Flits.and.survey.zone=Flits.and.survey.zone%>%filter(!Fleet.name=="Survey")
  }
  if("Survey"%in% names(Catch.rate.series[[i]]) | "Size_composition_Survey"%in%names(Species.data[[i]]))
  {
    Flits.zone=c(Flits.zone,length(Flits.zone)+1)
    names(Flits.zone)[length(Flits.zone)]="Survey"
  }
  ktch.zone=ktch.zone%>%
    spread(Fishry,Tonnes,fill=0)
  
  
  #2. Size composition   
  #note: commercial catch and survey. Nsamp set at number of shots
  #zones together
  Size.compo.SS.format=NULL
  Life.history$Max.population.TL=Life.history$Min.population.TL=NULL
  if(any(grepl('Size_composition',names(Species.data[[i]]))))
  {
    if(Neim=="sandbar shark")
    {
      Species.data[[i]]$Size_composition_Survey=Species.data[[i]]$Size_composition_Survey%>%
        mutate(SEX=ifelse(is.na(SEX) & FINYEAR=='2005-06' & FL==200,'F',SEX))
    }
    d.list.n.shots=Species.data[[i]][grep(paste(c("Size_composition_Survey_Observations","Size_composition_Observations",
                                                  "Size_composition_Other_Observations"),collapse="|"),
                                          names(Species.data[[i]]))]
    d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[i]]))]
    if(length(d.list)>0)
    {
      if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
      if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
      
      Max.obs.FL=max(unlist(lapply(d.list,function(x) max(x$FL,na.rm=T))))
      Max.obs.TL=with(Life.history,Max.obs.FL*a_FL.to.TL+b_FL.to.TL)
      Min.obs.FL=min(unlist(lapply(d.list,function(x) min(x$FL,na.rm=T))))
      Min.obs.TL=with(Life.history,Min.obs.FL*a_FL.to.TL+b_FL.to.TL)
      
      Min.population.TL=min(Min.obs.TL,with(Life.history,Lzero*a_FL.to.TL+b_FL.to.TL)) 
      Max.population.TL=max(Max.obs.TL,with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL))))
      Max.population.TL=TL.bins.cm*ceiling(Max.population.TL/TL.bins.cm)
      Min.population.TL=TL.bins.cm*floor(Min.population.TL/TL.bins.cm)
      
      Life.history$Max.population.TL=Max.population.TL
      Life.history$Min.population.TL=Min.population.TL
      
      for(s in 1:length(d.list))
      {
        d.list[[s]]=d.list[[s]]%>%
          mutate(fishry=ifelse(grepl("NSF.LONGLINE",names(d.list)[s]),'NSF',
                        ifelse(grepl("Survey",names(d.list)[s]),'Survey',
                        ifelse(grepl("Other",names(d.list)[s]),'Other',
                                    'TDGDLF'))),
                 year=as.numeric(substr(FINYEAR,1,4)),
                 TL=FL*Life.history$a_FL.to.TL+Life.history$b_FL.to.TL,
                 size.class=TL.bins.cm*floor(TL/TL.bins.cm))%>%
          filter(TL>=Min.population.TL & TL<=Max.population.TL)%>%
          dplyr::rename(sex=SEX)%>%
          filter(!is.na(sex))%>%
          group_by(year,fishry,sex,size.class)%>%
          summarise(n=n())%>%
          ungroup()
        
        #add extra bins for smooth fit to size comps  
        Maximum_size=ceiling(max(Max.population.TL,with(Life.history,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)*1.06))
        Mx.size=max(d.list[[s]]$size.class)
        extra.bins=NA
        if(Maximum_size>(Mx.size))
        {
          extra.bins=seq((Mx.size+TL.bins.cm),Maximum_size,by=TL.bins.cm)
          #extra.bins=seq(Mx.size+TL.bins.cm,10*round(Maximum_size/10),by=TL.bins.cm) 
        }
          
        if(any(!is.na(extra.bins)))
        {
          add.dumi.size=d.list[[s]][1:length(extra.bins),]%>%
            mutate(size.class=extra.bins,
                   n=0)
          d.list[[s]]=rbind(d.list[[s]],add.dumi.size)
        }
      }
      d.list <- d.list[!is.na(d.list)]
      
      if(Neim%in%names(Indicator.species))
      {
        Min.size=Min.annual.obs.ktch
      }else
      {
        Min.size=Min.annual.obs.ktch*prop.min.N.accepted_other
      }
      Min.size.NSF=Min.annual.obs.ktch_NSF
      if(Neim%in%c("dusky shark")) Min.size.NSF=20
      
      # Display sex ratio by zone used in SS  
      if(First.run=="YES")
      {
        HandL=handl_OneDrive("Analyses/Population dynamics/1.")
        DiR=paste(HandL,capitalize(Neim),"/",AssessYr,"/1_Inputs/Visualise data",sep='')
        add.n.samps=Species.data[[i]]$Size_composition_Observations%>%
                          filter(Method=='GN')%>%
                          mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
                          rename(Zone=zone)%>%
                          group_by(Zone,year)%>%
                          summarise(N.shots=sum(N.shots))%>%ungroup()
        fn.ktch.sex.ratio.zone_SS(size.data=d.list,Min.size=Min.size,N_sampleS=add.n.samps)
        ggsave(paste(DiR,'Sex ratio by zone_SS size comps data_single area.tiff',sep='/'), width = 5,height = 6, dpi = 300, compression = "lzw")
      }
      
      d.list=do.call(rbind,d.list)   
      if(Neim%in%combine_NSF_Survey) 
      {
        d.list=d.list%>%
          mutate(fishry=ifelse(fishry=='NSF',"Survey",fishry))
      }
      if(Neim%in%combine.sexes)
      {
        if(Neim%in%combine.sexes.survey)
        {
          d.list$sex=ifelse(d.list$fishry=="Survey",combine.sex_type,d.list$sex)
        }
        if(Neim%in%combine.sexes.nsf)
        {
          d.list$sex=ifelse(d.list$fishry=="NSF",combine.sex_type,d.list$sex)
        }
        if(Neim%in%combine.sexes.tdgdlf)
        {
          d.list=d.list%>%mutate(sex=ifelse(fishry=="TDGDLF",combine.sex_type,sex))
        }
        if(Neim%in%combine.sexes.tdgdlf.daily)
        {
          d.list=d.list%>%mutate(sex=ifelse(fishry=="TDGDLF" & year>2005,combine.sex_type,sex))
        }
        if(!Neim%in%c(combine.sexes.survey,combine.sexes.tdgdlf,combine.sexes.tdgdlf.daily))
        {
          d.list$sex=combine.sex_type 
        }
      }
      
      #combine sexes if number of obs per year >Min.size but per sex <Min.size
      Table.n=d.list%>%  
              group_by(year,fishry)%>%
              summarise(N=sum(n))%>%
              mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey,
                                              fishry=='NSF'~Min.size.NSF,
                                              TRUE~Min.size))%>%   
              filter(N>=Min.accepted.N)%>%
              mutate(dummy=paste(year,fishry)) 
      
      Table.n.sex=d.list%>%  
                group_by(year,fishry,sex)%>%
                summarise(N=sum(n))%>%
                mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey,
                                                fishry=='NSF'~Min.size.NSF,
                                                TRUE~Min.size))%>%   
                filter(N<Min.accepted.N)%>%
                mutate(dummy=paste(year,fishry))%>%
                distinct(dummy)
              
      d.list=d.list%>%
        mutate(dummy=paste(year,fishry),
               sex=ifelse(dummy%in%Table.n.sex$dummy,combine.sex_type,sex))%>%
        dplyr::select(-dummy)
      
      d.list=d.list%>%
        group_by(year,fishry,sex,size.class)%>%
        summarise(n=sum(n))%>%
        ungroup()%>%
        filter(!is.na(year))%>%
        filter(!is.na(fishry))
      
      MAXX=max(d.list$size.class)
      Tab.si.kl=table(d.list$size.class)
      if(sum(Tab.si.kl[(length(Tab.si.kl)-1):length(Tab.si.kl)])/sum(Tab.si.kl)>.05) MAXX=MAXX*1.2
      MAXX=min(MAXX,30*ceiling(with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/30))
      size.classes=seq(min(d.list$size.class),MAXX,by=TL.bins.cm)
      missing.size.classes=size.classes[which(!size.classes%in%unique(d.list$size.class))]
      if(fill.in.zeros & length(missing.size.classes)>0)
      {
        d.list=rbind(d.list,d.list[1:length(missing.size.classes),]%>%
                       mutate(size.class=missing.size.classes,
                              n=0))
      }

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
                                        ifelse(zone%in%c('SA','GAB.trawl'),'Other',
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
          mutate(Part=SS.part_length.comps)%>%
          dplyr::rename(Fleet=fishry,
                        Sex=sex,
                        Seas=month)%>%
          relocate(all_of(vars.head))
        #keep years with minimum number of observations
        d.list=d.list%>%
          mutate(dummy2=paste(year,Fleet))%>%
          filter(dummy2%in%unique(Table.n$dummy))%>%
          dplyr::select(-dummy2)%>%
          mutate(Sex=ifelse(Sex=='F',1,
                     ifelse(Sex=='M',2,
                            Sex)))
        d.list=d.list%>%
          mutate(dumi.n=rowSums(d.list[,-match(c('year','Seas','Fleet','Sex','Part','Nsamp'),names(d.list))]),
                 Nsamp=ifelse(Nsamp>dumi.n,dumi.n,Nsamp))%>%
          dplyr::select(-dumi.n)
        size.flits=Flits.and.survey
        if(!"Survey"%in% names(Catch.rate.series[[i]]) & "Survey"%in%unique(d.list$Fleet))
        {
          ddummis=size.flits[1,]%>%mutate(Fleet.number=1+size.flits$Fleet.number[nrow(size.flits)],
                                          Fleet.name="Survey")
          rownames(ddummis)="Survey"
          if(!'Survey'%in%size.flits$Fleet.name) size.flits=rbind(size.flits,ddummis)
        }
        d.list=d.list%>%
          mutate(dummy.fleet=case_when(Fleet=="NSF"~'Northern.shark',
                                       Fleet=="Other"~'Other',
                                       Fleet=="TDGDLF" & year<2006~'Southern.shark_1',
                                       Fleet=="TDGDLF" & year>=2006~'Southern.shark_2',
                                       Fleet=="Survey"~'Survey'))%>%
          left_join(size.flits,by=c('dummy.fleet'='Fleet.name'))%>%
          mutate(Fleet=Fleet.number)%>%
          dplyr::select(-c(dummy.fleet,Fleet.number))%>%
          arrange(Sex,Fleet,year)
        
        d.list.0=d.list%>%filter(Sex==combine.sex_type)%>%arrange(year)
        d.list.f=d.list%>%filter(Sex==1)%>%arrange(year)  
        d.list.m=d.list%>%filter(Sex==2)%>%arrange(year)
        
        id.var.nms=which(names(d.list.f)%in%vars.head)
        id.var.nms.f=which(!names(d.list.f)%in%vars.head)
        id.var.nms.m=length(id.var.nms.f)+which(!names(d.list.m)%in%vars.head)
        
        dummy.zeros.0=d.list.0[,-match(vars.head,names(d.list.0))]
        dummy.zeros.0[,]=0
        dummy.zeros.f=d.list.f[,-match(vars.head,names(d.list.f))]
        dummy.zeros.f[,]=0
        dummy.zeros.m=d.list.m[,-match(vars.head,names(d.list.m))]
        dummy.zeros.m[,]=0
        
        if(nrow(d.list.0)>0)  
        {
          dummy.Size.compo.SS.format_Sex0=cbind(d.list.0,dummy.zeros.0)
          names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.f]=paste('f',names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.f],sep='')
          names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.m]=paste('m',names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.m],sep='')
        }
        if(nrow(d.list.f)>0)
        {
          dummy.Size.compo.SS.format_Sex=rbind(cbind(d.list.f,dummy.zeros.f),
                                               cbind(d.list.m[,match(vars.head,names(d.list.m))],
                                                     dummy.zeros.m,
                                                     d.list.m[,-match(vars.head,names(d.list.m))]))
          names(dummy.Size.compo.SS.format_Sex)[id.var.nms.f]=paste('f',names(dummy.Size.compo.SS.format_Sex)[id.var.nms.f],sep='')
          names(dummy.Size.compo.SS.format_Sex)[id.var.nms.m]=paste('m',names(dummy.Size.compo.SS.format_Sex)[id.var.nms.m],sep='')
        }  
        if(!exists('dummy.Size.compo.SS.format_Sex0') & exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format_Sex
        }
        if(exists('dummy.Size.compo.SS.format_Sex0') & !exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format_Sex0
        }
        if(exists('dummy.Size.compo.SS.format_Sex0') & exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=rbind(dummy.Size.compo.SS.format_Sex0,dummy.Size.compo.SS.format_Sex)
        }
        clear.log('dummy.Size.compo.SS.format_Sex0')
        clear.log('dummy.Size.compo.SS.format_Sex')
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          arrange(Fleet,year,Sex)
        
        #select min sample size (shots) 
        min.nsamp=Min.Nsamp
        if(!Neim%in%names(Indicator.species)) min.nsamp=ceiling(min.nsamp/2)
        size.flits.min.samp=size.flits%>%
          mutate(Min.nsamp=case_when(Fleet.name=='Northern.shark'~Min.Nsamp.NSF,
                                     Fleet.name=='Survey'~Min.Nsamp.Survey,
                                     TRUE~min.nsamp))%>%
          dplyr::select(-Fleet.name)%>%
          rename(Fleet=Fleet.number)%>%
          filter(Fleet%in%unique(dummy.Size.compo.SS.format$Fleet))
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          left_join(size.flits.min.samp,by='Fleet')
        
        dummy.Size.compo.SS.format.all=dummy.Size.compo.SS.format%>%
                                        dplyr::select(-Min.nsamp)
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%  
          filter(year<=max(ktch$finyear) & Nsamp>=Min.nsamp)%>%
          dplyr::select(-Min.nsamp)
        
        if(nrow(dummy.Size.compo.SS.format)>0) Size.compo.SS.format=dummy.Size.compo.SS.format
        
      }
    }
  }
  # by zones 
  Size.compo.SS.format.zone=NULL
  if(any(grepl('Size_composition',names(Species.data[[i]]))))
  {
    d.list.n.shots=Species.data[[i]][grep(paste(c("Size_composition_Survey_Observations","Size_composition_Observations",
                                                  "Size_composition_Other_Observations"),collapse="|"),
                                          names(Species.data[[i]]))]
    d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[i]]))]
    if(Neim%in%names(drop.dodgy.len.comp))
    {
      diszone=drop.dodgy.len.comp[[match(Neim,names(drop.dodgy.len.comp))]]
      dis.yrs=unique(word(diszone, 1, sep = "-"))
      diszone=unique(str_remove(diszone, ".*-"))
      for(xx in 1:length(diszone))
      {
        for(yy in 1:length(dis.yrs))
        {
          for(qq in 1:length(d.list))
          {
            if(grepl(diszone[xx],names(d.list)[qq]))
            {
              d.list[[qq]]=d.list[[qq]]%>%
                            filter(!FINYEAR%in%paste(as.numeric(dis.yrs[yy]),
                                                     substr(as.numeric(dis.yrs[yy])+1,3,4),sep='-'))
            }
          }
        }
      }
    }
      
    if(length(d.list)>0)
    {
      if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
      if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
      
      for(s in 1:length(d.list))
      {
        NM=names(d.list)[s]
        d.list[[s]]=d.list[[s]]%>%
          filter(FL>=Life.history$Lzero)%>%
          mutate(fishry=ifelse(grepl("NSF.LONGLINE",NM),'NSF',
                               ifelse(grepl("Survey",NM),'Survey',
                                      ifelse(grepl("Other",NM),'Other',
                                             'TDGDLF'))),
                 year=as.numeric(substr(FINYEAR,1,4)),
                 TL=FL*Life.history$a_FL.to.TL+Life.history$b_FL.to.TL,
                 size.class=TL.bins.cm*floor(TL/TL.bins.cm))%>%
          filter(TL<=with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL))))%>%
          dplyr::rename(sex=SEX)%>%
          filter(!is.na(sex))%>%
          group_by(year,fishry,sex,size.class)%>%
          summarise(n=n())%>%
          ungroup()
        if(grepl(paste(c('West','Zone1','Zone2'),collapse='|'),NM))
        {
          ZnE=sub("\\..*", "", str_remove(NM,"Size_composition_"))
          d.list[[s]]=d.list[[s]]%>%
            mutate(fishry=paste(fishry,ZnE,sep='_'))
        }
        
        
        #add extra bins for smooth fit to size comps
        Maximum_size=ceiling(max(Max.population.TL,with(Life.history,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)*1.06))
        Mx.size=max(d.list[[s]]$size.class)
        extra.bins=NA
        if(Maximum_size>(Mx.size))
        {
          extra.bins=seq((Mx.size+TL.bins.cm),Maximum_size,by=TL.bins.cm)
          #extra.bins=seq(Mx.size+TL.bins.cm,10*round(Maximum_size/10),by=TL.bins.cm)
        }
           
        if(any(!is.na(extra.bins)))
        {
          add.dumi.size=d.list[[s]][1:length(extra.bins),]%>%
            mutate(size.class=extra.bins,
                   n=0)
          d.list[[s]]=rbind(d.list[[s]],add.dumi.size)
        }
      }
      d.list <- d.list[!is.na(d.list)]
      if(Neim%in%names(Indicator.species))
      {
        Min.size=Min.annual.obs.ktch.zone
      }else
      {
        Min.size=Min.annual.obs.ktch.zone*prop.min.N.accepted_other
      }
      Min.size.NSF=Min.annual.obs.ktch_NSF
      if(Neim%in%c("dusky shark")) Min.size.NSF=20
      
      
      # Display sex ratio by zone used in SS  
      if(First.run=="YES")
      {
        HandL=handl_OneDrive("Analyses/Population dynamics/1.")
        DiR=paste(HandL,capitalize(Neim),"/",AssessYr,"/1_Inputs/Visualise data",sep='')
        add.n.samps=Species.data[[i]]$Size_composition_Observations%>%
                        filter(Method=='GN')%>%
                        mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
                        rename(Zone=zone)%>%
                        group_by(Zone,year)%>%
                        summarise(N.shots=sum(N.shots))%>%ungroup()
        fn.ktch.sex.ratio.zone_SS(size.data=d.list,Min.size=Min.size,N_sampleS=add.n.samps)
        ggsave(paste(DiR,'Sex ratio by zone_SS size comps data_areas as fleets.tiff',sep='/'), width = 5,height = 6, dpi = 300, compression = "lzw")
      }
      d.list=do.call(rbind,d.list)   
      if(Neim%in%combine_NSF_Survey) 
      {
        d.list=d.list%>%
          mutate(fishry=ifelse(fishry=='NSF',"Survey",fishry))
      }
      if(Neim%in%combine.sexes)
      {
        if(Neim%in%combine.sexes.survey)
        {
          d.list$sex=ifelse(d.list$fishry=="Survey",combine.sex_type,d.list$sex)
        }
        if(Neim%in%combine.sexes.nsf)
        {
          d.list$sex=ifelse(d.list$fishry=="NSF",combine.sex_type,d.list$sex)
        }
        if(Neim%in%combine.sexes.tdgdlf)
        {
          d.list=d.list%>%mutate(sex=ifelse(grepl("TDGDLF",fishry),combine.sex_type,sex))
        }
        if(Neim%in%combine.sexes.tdgdlf.daily)
        {
          d.list=d.list%>%mutate(sex=ifelse(grepl("TDGDLF",fishry) & year>2005,combine.sex_type,sex))
        }
        if(!Neim%in%c(combine.sexes.survey,combine.sexes.tdgdlf,combine.sexes.tdgdlf.daily))
        {
          d.list$sex=combine.sex_type 
        }
      }
      #combine sexes if number of obs per year >Min.size but per sex <Min.size
      Table.n=d.list%>%
                group_by(year,fishry)%>%
                summarise(N=sum(n))%>%
                mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey,
                                                fishry=='NSF'~Min.size.NSF,
                                                TRUE~Min.size))%>%
                filter(N>=Min.accepted.N)%>%
                mutate(dummy=paste(year,fishry))
      Table.n.sex=d.list%>%  
        group_by(year,fishry,sex)%>%
        summarise(N=sum(n))%>%
        mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey,
                                        fishry=='NSF'~Min.size.NSF,
                                        TRUE~Min.size))%>%   
        filter(N<Min.accepted.N)%>%
        mutate(dummy=paste(year,fishry))%>%
        distinct(dummy)
      d.list=d.list%>%
        mutate(dummy=paste(year,fishry),
               sex=ifelse(dummy%in%Table.n.sex$dummy,combine.sex_type,sex))%>%
        dplyr::select(-dummy)
      
      d.list=d.list%>%
        group_by(year,fishry,sex,size.class)%>%
        summarise(n=sum(n))%>%
        ungroup()%>%
        filter(!is.na(year))%>%
        filter(!is.na(fishry))
      
      MAXX=max(d.list$size.class)
      Tab.si.kl=table(d.list$size.class)
      if(sum(Tab.si.kl[(length(Tab.si.kl)-1):length(Tab.si.kl)])/sum(Tab.si.kl)>.05) MAXX=MAXX*1.2
      MAXX=min(MAXX,30*ceiling(with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/30))
      size.classes=seq(min(d.list$size.class),MAXX,by=TL.bins.cm)
      missing.size.classes=size.classes[which(!size.classes%in%unique(d.list$size.class))]
      if(fill.in.zeros & length(missing.size.classes)>0)
      {
        d.list=rbind(d.list,d.list[1:length(missing.size.classes),]%>%
                       mutate(size.class=missing.size.classes,
                              n=0))
      }
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
            mutate(fishry=ifelse(Method=='GN',paste('TDGDLF',zone,sep='_'),
                                 ifelse(Method=='LL' & zone!='SA','NSF',
                                        ifelse(zone%in%c('SA','GAB.trawl'),'Other',
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
          mutate(Part=SS.part_length.comps)%>%
          dplyr::rename(Fleet=fishry,
                        Sex=sex,
                        Seas=month)%>%
          relocate(all_of(vars.head))
        #keep years with minimum number of observations
        d.list=d.list%>%
          mutate(dummy2=paste(year,Fleet))%>%
          filter(dummy2%in%unique(Table.n$dummy))%>%
          dplyr::select(-dummy2)%>%
          mutate(Sex=ifelse(Sex=='F',1,
                     ifelse(Sex=='M',2,
                            Sex)))
        d.list=d.list%>%
          mutate(dumi.n=rowSums(d.list[,-match(c('year','Seas','Fleet','Sex','Part','Nsamp'),names(d.list))]),
                 Nsamp=ifelse(Nsamp>dumi.n,dumi.n,Nsamp))%>%
          dplyr::select(-dumi.n)
        size.flits.zone=Flits.and.survey.zone
        if(!"Survey"%in% names(Catch.rate.series[[i]]) & "Survey"%in%unique(d.list$Fleet))
        {
          ddummis=size.flits.zone[1,]%>%mutate(Fleet.number=1+size.flits.zone$Fleet.number[nrow(size.flits.zone)],
                                               Fleet.name="Survey")
          rownames(ddummis)="Survey"
          if(!'Survey'%in%size.flits.zone$Fleet.name) size.flits.zone=rbind(size.flits.zone,ddummis)
        }
        d.list=d.list%>%
          mutate(dummy.fleet=case_when(Fleet=="NSF"~'Northern.shark',
                                       Fleet=="Other"~'Other',
                                       Fleet=="TDGDLF_West" & year<2006~'Southern.shark_1_West',
                                       Fleet=="TDGDLF_West" & year>=2006~'Southern.shark_2_West',
                                       Fleet=="TDGDLF_Zone1" & year<2006~'Southern.shark_1_Zone1',
                                       Fleet=="TDGDLF_Zone1" & year>=2006~'Southern.shark_2_Zone1',
                                       Fleet=="TDGDLF_Zone2" & year<2006~'Southern.shark_1_Zone2',
                                       Fleet=="TDGDLF_Zone2" & year>=2006~'Southern.shark_2_Zone2',
                                       Fleet=="Survey"~'Survey'))%>%
          left_join(size.flits.zone,by=c('dummy.fleet'='Fleet.name'))%>%
          mutate(Fleet=Fleet.number)%>%
          dplyr::select(-c(dummy.fleet,Fleet.number))%>%
          arrange(Sex,Fleet,year)
        
        d.list.0=d.list%>%filter(Sex==combine.sex_type)%>%arrange(year)
        d.list.f=d.list%>%filter(Sex==1)%>%arrange(year)  
        d.list.m=d.list%>%filter(Sex==2)%>%arrange(year)
        
        id.var.nms=which(names(d.list.f)%in%vars.head)
        id.var.nms.f=which(!names(d.list.f)%in%vars.head)
        id.var.nms.m=length(id.var.nms.f)+which(!names(d.list.m)%in%vars.head)
        
        dummy.zeros.0=d.list.0[,-match(vars.head,names(d.list.0))]
        dummy.zeros.0[,]=0
        dummy.zeros.f=d.list.f[,-match(vars.head,names(d.list.f))]
        dummy.zeros.f[,]=0
        dummy.zeros.m=d.list.m[,-match(vars.head,names(d.list.m))]
        dummy.zeros.m[,]=0
        
        if(nrow(d.list.0)>0)  
        {
          dummy.Size.compo.SS.format_Sex0=cbind(d.list.0,dummy.zeros.0)
          names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.f]=paste('f',names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.f],sep='')
          names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.m]=paste('m',names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.m],sep='')
        }
        if(nrow(d.list.f)>0)
        {
          dummy.Size.compo.SS.format_Sex=rbind(cbind(d.list.f,dummy.zeros.f),
                                               cbind(d.list.m[,match(vars.head,names(d.list.m))],
                                                     dummy.zeros.m,
                                                     d.list.m[,-match(vars.head,names(d.list.m))]))
          names(dummy.Size.compo.SS.format_Sex)[id.var.nms.f]=paste('f',names(dummy.Size.compo.SS.format_Sex)[id.var.nms.f],sep='')
          names(dummy.Size.compo.SS.format_Sex)[id.var.nms.m]=paste('m',names(dummy.Size.compo.SS.format_Sex)[id.var.nms.m],sep='')
        }  
        if(!exists('dummy.Size.compo.SS.format_Sex0') & exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format_Sex
        }
        if(exists('dummy.Size.compo.SS.format_Sex0') & !exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format_Sex0
        }
        if(exists('dummy.Size.compo.SS.format_Sex0') & exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=rbind(dummy.Size.compo.SS.format_Sex0,dummy.Size.compo.SS.format_Sex)
        }
        clear.log('dummy.Size.compo.SS.format_Sex0')
        clear.log('dummy.Size.compo.SS.format_Sex')
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          arrange(Fleet,year,Sex)
        
        #select min sample size (shots)
        min.nsamp=Min.Nsamp.zone
        if(!Neim%in%names(Indicator.species)) min.nsamp=ceiling(min.nsamp/2)
        size.flits.min.samp=size.flits.zone%>%
          mutate(Min.nsamp=case_when(Fleet.name=='Northern.shark'~Min.Nsamp.NSF,
                                     Fleet.name=='Survey'~Min.Nsamp.Survey,
                                     TRUE~min.nsamp))%>%
          dplyr::select(-Fleet.name)%>%
          rename(Fleet=Fleet.number)%>%
          filter(Fleet%in%unique(dummy.Size.compo.SS.format$Fleet))
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          left_join(size.flits.min.samp,by='Fleet')
        
        dummy.Size.compo.SS.format.all_zone=dummy.Size.compo.SS.format%>%
                                              dplyr::select(-Min.nsamp)
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          filter(year<=max(ktch$finyear) & Nsamp>=Min.nsamp)%>%
          dplyr::select(-Min.nsamp)
        
        
        if(nrow(dummy.Size.compo.SS.format)>0) Size.compo.SS.format.zone=dummy.Size.compo.SS.format
        
      }
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
    if(!is.null(Size.compo.SS.format.zone))
    {
      Fleet.more.one.year.obs=Size.compo.SS.format.zone%>%
        distinct(Fleet,year)%>%
        group_by(Fleet)%>%
        tally()%>%
        filter(n>1)%>%
        pull(Fleet)
      Size.compo.SS.format.zone=Size.compo.SS.format.zone%>%
        filter(Fleet%in%Fleet.more.one.year.obs)
    }
  }
  
  #Reset sex type if required  
  if(SS.sex.length.type==3)
  {
    #zones combined  
    Kombo=Size.compo.SS.format%>%distinct(year,Seas,Fleet)
    Kombo.list=vector('list',nrow(Kombo))  
    for(kk in 1:length(Kombo.list))
    {
      dumi.kk=Size.compo.SS.format%>%
        filter(year==Kombo$year[kk] & Seas==Kombo$Seas[kk] & Fleet==Kombo$Fleet[kk])
      if(nrow(dumi.kk)==1 & SS.sex.3_use.missing.sex) 
      {
        if(!dumi.kk$Sex==0 & !SS.sex.1_2.dont.change.to_3) 
        {
          msin.sx=which(!1:2%in%dumi.kk$Sex)
          dumi.kk.missing=dumi.kk%>%
            mutate(Sex=msin.sx,
                   across(-all_of(c('year','Seas','Fleet','Sex','Part','Nsamp')), ~0))
          dumi.kk=rbind(dumi.kk,dumi.kk.missing)%>%
            arrange(Sex)
        }
      }
      if(!0%in%dumi.kk$Sex)
      {
        if(!SS.sex.1_2.dont.change.to_3) if(nrow(dumi.kk)<2) dumi.kk=NA   
        if(any(!is.na(dumi.kk)))
        {
          if(nrow(dumi.kk)==2)
          {
            dumi.kk.1=dumi.kk%>%filter(Sex==1)
            dumi.kk.2=dumi.kk%>%filter(Sex==2)
            dumi.kk.3=dumi.kk.1
            dumi.kk.3[,grepl('m',colnames(dumi.kk.3))]=dumi.kk.2[,grepl('m',colnames(dumi.kk.2))]
            dumi.kk.3$Sex=SS.sex.length.type
            dumi.kk=dumi.kk.3
          }
        }
      }
      Kombo.list[[kk]]=dumi.kk
    }
    Kombo.list=Kombo.list[!is.na(Kombo.list)]
    Size.compo.SS.format=do.call(rbind,Kombo.list)
    
    #by zone
    Kombo=Size.compo.SS.format.zone%>%distinct(year,Seas,Fleet)
    Kombo.list=vector('list',nrow(Kombo))  
    for(kk in 1:length(Kombo.list))
    {
      dumi.kk=Size.compo.SS.format.zone%>%
        filter(year==Kombo$year[kk] & Seas==Kombo$Seas[kk] & Fleet==Kombo$Fleet[kk])
      if(nrow(dumi.kk)==1 & SS.sex.3_use.missing.sex.zone)
      {
        if(!dumi.kk$Sex==0 & !SS.sex.1_2.dont.change.to_3)
        {
          msin.sx=which(!1:2%in%dumi.kk$Sex)
          dumi.kk.missing=dumi.kk%>%
            mutate(Sex=msin.sx,
                   across(-all_of(c('year','Seas','Fleet','Sex','Part','Nsamp')), ~0))
          dumi.kk=rbind(dumi.kk,dumi.kk.missing)%>%
            arrange(Sex) 
        }
      }
      if(!0%in%dumi.kk$Sex)
      {
        if(!SS.sex.1_2.dont.change.to_3) if(nrow(dumi.kk)<2) dumi.kk=NA
        if(any(!is.na(dumi.kk)))
        {
          if(nrow(dumi.kk)==2)
          {
            dumi.kk.1=dumi.kk%>%filter(Sex==1)
            dumi.kk.2=dumi.kk%>%filter(Sex==2)
            dumi.kk.3=dumi.kk.1
            dumi.kk.3[,grepl('m',colnames(dumi.kk.3))]=dumi.kk.2[,grepl('m',colnames(dumi.kk.2))]
            dumi.kk.3$Sex=SS.sex.length.type
            dumi.kk=dumi.kk.3
          }
        }
      }
      Kombo.list[[kk]]=dumi.kk
    }
    Kombo.list=Kombo.list[!is.na(Kombo.list)]
    Size.compo.SS.format.zone=do.call(rbind,Kombo.list)
  }
  
  
  #3. meanbodywt
  meanbody.part=SS.part_meanbodywt
  if(!Neim%in%retained.discarded.sp) meanbody.part=0   
  #zones together
  meanbodywt.SS.format=NULL
  if(any(grepl('annual.mean.size',names(Species.data[[i]]))))
  {
    meanbodywt.SS.format=Species.data[[i]]$annual.mean.size%>%
      mutate(year=as.numeric(substr(Finyear,1,4)),
             month=1,
             Fleet='Southern.shark_2',
             part=meanbody.part,   
             type=2)%>%
      filter(year<=max(ktch$finyear))%>%
      dplyr::select(-Finyear)%>%
      left_join(Flits.and.survey,by=c('Fleet'='Fleet.name'))%>%
      mutate(fleet=Fleet.number)%>%
      dplyr::select(-c(Fleet.number,Fleet))%>%
      relocate(year,month,fleet,part,type,mean,CV)
    
    #reset very low CVs
    NewCVs=Francis.function(cipiuis=meanbodywt.SS.format%>%
                              dplyr::select(year,mean,fleet)%>%
                              spread(fleet,mean),
                            cvs=meanbodywt.SS.format%>%
                              dplyr::select(year,CV,fleet)%>%
                              spread(fleet,CV),
                            mininum.mean.CV=default.Mean.weight.CV)
    if(mean(meanbodywt.SS.format$CV,na.rm=T)<default.Mean.weight.CV)
    {
      newcv=NewCVs$CV.Adjusted%>%
        filter(year%in%meanbodywt.SS.format$year)
      meanbodywt.SS.format=meanbodywt.SS.format%>%mutate(CV=newcv[,2])
    }
    #CV variance adjustment if small CVs
    if(mean(meanbodywt.SS.format$CV,na.rm=T)>=default.Mean.weight.CV) NewCVs$CV.var.adj=0
    Var.ad.factr_meanbodywt=data.frame(Factor=3,
                                       Fleet=unique(meanbodywt.SS.format$fleet),
                                       Value=NewCVs$CV.var.adj)
  }
  # by zones 
  meanbodywt.SS.format.zone=NULL
  if(any(grepl('annual.mean.size',names(Species.data[[i]]))))
  {
    xxx=Species.data[[i]][grep('annual.mean.size',names(Species.data[[i]]))]
    xxx=xxx[-match("annual.mean.size",names(xxx))]
    if(Neim%in%names(tdgdlf_daily_not.representative))
    {
      id.drop=paste0('annual.mean.size_',tolower(tdgdlf_daily_not.representative[[Neim]]))
      xxx=xxx[-grep(paste(id.drop,collapse='|'),names(xxx))] 
    }
    if(length(xxx)>0)
    {
      for(y in 1:length(xxx))
      {
        ZnE=capitalize(str_remove(names(xxx)[y],"annual.mean.size_"))
        xxx[[y]]=xxx[[y]]%>%
          mutate(year=as.numeric(substr(Finyear,1,4)),
                 month=1,
                 Fleet=paste('Southern.shark_2',ZnE,sep='_'),
                 part=meanbody.part,      
                 type=2)%>%
          filter(year<=max(ktch.zone$finyear))%>%
          dplyr::select(-Finyear)
      }
      
      
      meanbodywt.SS.format.zone=do.call(rbind,xxx)%>%
        left_join(Flits.and.survey.zone,by=c('Fleet'='Fleet.name'))%>%
        mutate(fleet=Fleet.number)%>%
        dplyr::select(-c(Fleet.number,Fleet,zone))%>%
        relocate(year,month,fleet,part,type,mean,CV)
      
      #reset very low CVs
      NewCVs=Francis.function(cipiuis=meanbodywt.SS.format.zone%>%
                                dplyr::select(year,mean,fleet)%>%
                                spread(fleet,mean),
                              cvs=meanbodywt.SS.format.zone%>%
                                dplyr::select(year,CV,fleet)%>%
                                spread(fleet,CV),
                              mininum.mean.CV=default.Mean.weight.CV)
      if(mean(meanbodywt.SS.format.zone$CV,na.rm=T)<default.Mean.weight.CV)
      {
        newcv=NewCVs$CV.Adjusted%>%
          filter(year%in%meanbodywt.SS.format.zone$year)%>%
          gather(fleet,CV,-year)%>%mutate(fleet=as.numeric(fleet))
        meanbodywt.SS.format.zone=meanbodywt.SS.format.zone%>%
          dplyr::select(-CV)%>%
          left_join(newcv,by=c('year','fleet'))%>%
          relocate(names(meanbodywt.SS.format))
      }
      #CV variance adjustment if small CVs
      if(mean(meanbodywt.SS.format.zone$CV,na.rm=T)>=default.Mean.weight.CV) NewCVs$CV.var.adj=0
      Var.ad.factr_meanbodywt.zone=data.frame(Factor=3,
                                              Fleet=unique(meanbodywt.SS.format.zone$fleet),
                                              Value=NewCVs$CV.var.adj)
      
    }
  }
  
  
  #4. Abundance series           
  Abundance.SS.format=Abundance.SS.format.zone=NULL
  Var.ad.factr=Var.ad.factr.zone=NULL
  CPUE=compact(Catch.rate.series[[i]])
  if(!is.null(CPUE))
  {
    if(Neim%in%NSF_not.representative & any(grepl("NSF",names(CPUE)))) CPUE=CPUE[-grep("NSF",names(CPUE))]
    if(Neim%in%tdgdlf_monthly_not.representative & "TDGDLF.monthly"%in%names(CPUE)) CPUE=CPUE[-grep("TDGDLF.monthly",names(CPUE))]
    if(Neim%in%names(drop.dodgy.cpue))
    {
      diszone.yr=drop.dodgy.cpue[[match(Neim,names(drop.dodgy.cpue))]]
      dis.data.set=unique(word(diszone.yr, 1, sep = "-"))
      dis.yrs=unique(word(diszone.yr, 2, sep = "-"))
      for(xx in 1:length(dis.data.set))
      {
        for(yy in 1:length(dis.yrs))
        {
          CPUE[[paste0('TDGDLF.',dis.data.set[xx])]]=CPUE[[paste0('TDGDLF.',dis.data.set[xx])]]%>%
                                                      filter(!yr.f==as.numeric(dis.yrs))

        }
      }
    }
    if(Neim%in%names(tdgdlf_daily_not.representative))
    {
      id.drop=paste0('daily.',tdgdlf_daily_not.representative[[Neim]])
      CPUE=CPUE[-grep(paste(id.drop,collapse='|'),names(CPUE))] 
    }
    CPUE.zone=CPUE
    DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
    if(length(DROP)>0)CPUE=CPUE[-DROP]
    DROP.zone=match(c('observer',"TDGDLF.daily","TDGDLF.monthly"),names(CPUE.zone))
    DROP.zone=subset(DROP.zone,!is.na(DROP.zone))
    if(length(DROP.zone)>0)CPUE.zone=CPUE.zone[-DROP.zone]
    
    #reset very low CVs
    #note: Andre suggested leaving original CVs and estimating extraSD if more than one index available
    #      If only 1 index available, then do not estimate, just increase CV before fitting model
    #      ICCAT and SEDAR leave CVs as is and don't estimate extraSD but add variance adjustments factors to the control file
    #zones together
    if(length(CPUE)>0)
    {
      dumi.cpue=CPUE
      for(j in 1:length(CPUE))
      {
        dumi.cpue[[j]]=dumi.cpue[[j]]%>%
          mutate(Fleet=names(dumi.cpue)[j])%>%
          dplyr::select(yr.f,Mean,Fleet,CV)
      }
      NewCVs=Francis.function(cipiuis=do.call(rbind,dumi.cpue)%>%
                                dplyr::select(yr.f,Mean,Fleet)%>%
                                spread(Fleet,Mean),
                              cvs=do.call(rbind,dumi.cpue)%>%
                                dplyr::select(yr.f,CV,Fleet)%>%
                                spread(Fleet,CV),
                              mininum.mean.CV=default.CV) 
      for(j in 1:length(CPUE))
      {
        #CPUE[[j]]=CPUE[[j]]%>%mutate(CV.Original=CV)
        if(mean(CPUE[[j]]$CV,na.rm=T)<default.CV)
        {
          newcv=NewCVs$CV.Adjusted[,match(c("yr.f",names(CPUE)[j]),names(NewCVs$CV.Adjusted))]%>%
            filter(yr.f%in%CPUE[[j]]$yr.f)
          CPUE[[j]]=CPUE[[j]]%>%mutate(CV=newcv[,2])
        }
      }
      #CV variance adjustment if small CVs
      for(j in 1:length(CPUE))
      {
        if(mean(CPUE[[j]]$CV,na.rm=T)>=default.CV)  NewCVs$CV.var.adj[j]=0
      }
      Var.ad.factr_cpue=data.frame(Factor=1,Fleet=names(NewCVs$CV.var.adj),Value=NewCVs$CV.var.adj)%>%
        mutate(Fleet=case_when(Fleet=="NSF"~"Northern.shark",
                               Fleet=="TDGDLF.monthly"~"Southern.shark_1",
                               Fleet=="TDGDLF.daily"~"Southern.shark_2",
                               TRUE~Fleet))%>%
        left_join(Flits.and.survey,by=c('Fleet'='Fleet.name'))%>%
        mutate(Fleet=Fleet.number)%>%
        dplyr::select(-Fleet.number)
      Var.ad.factr=Var.ad.factr_cpue           
      if(exists("Var.ad.factr_meanbodywt")) Var.ad.factr=rbind(Var.ad.factr,Var.ad.factr_meanbodywt)
      clear.log("Var.ad.factr_meanbodywt")  
      clear.log("Var.ad.factr_cpue")
      if(drop.intermediate.yrs)if(Whiskery.q.periods==2 & Neim=='whiskery shark') 
      {
        CPUE$TDGDLF.monthly=CPUE$TDGDLF.monthly%>%
          filter(!yr.f%in%Life.history$Yr_q_change_transition)
      }
    }
    # by zones 
    if(length(CPUE.zone)>0)
    {
      dumi.cpue=CPUE.zone
      for(j in 1:length(CPUE.zone))
      {
        dumi.cpue[[j]]=dumi.cpue[[j]]%>%
          mutate(Fleet=names(dumi.cpue)[j])%>%
          dplyr::select(yr.f,Mean,Fleet,CV)
      }
      NewCVs=Francis.function(cipiuis=do.call(rbind,dumi.cpue)%>%
                                dplyr::select(yr.f,Mean,Fleet)%>%
                                spread(Fleet,Mean),
                              cvs=do.call(rbind,dumi.cpue)%>%
                                dplyr::select(yr.f,CV,Fleet)%>%
                                spread(Fleet,CV),
                              mininum.mean.CV=default.CV) 
      for(j in 1:length(CPUE.zone))
      {
        #CPUE.zone[[j]]=CPUE.zone[[j]]%>%mutate(CV.Original=CV)
        if(mean(CPUE.zone[[j]]$CV,na.rm=T)<default.CV)
        {
          newcv=NewCVs$CV.Adjusted[,match(c("yr.f",names(CPUE.zone)[j]),names(NewCVs$CV.Adjusted))]%>%
            filter(yr.f%in%CPUE.zone[[j]]$yr.f)
          CPUE.zone[[j]]=CPUE.zone[[j]]%>%mutate(CV=newcv[,2])
        }
      }
      #CV variance adjustment if small CVs
      for(j in 1:length(CPUE.zone))
      {
        if(mean(CPUE.zone[[j]]$CV,na.rm=T)>=default.CV)  NewCVs$CV.var.adj[j]=0
      }
      Var.ad.factr_cpue.zone=data.frame(Factor=1,
                                        Fleet=names(NewCVs$CV.var.adj),
                                        Value=NewCVs$CV.var.adj)%>%
        mutate(Fleet=case_when(Fleet=="NSF"~"Northern.shark",
                               Fleet=="TDGDLF.monthly.West"~"Southern.shark_1_West",
                               Fleet=="TDGDLF.monthly.Zone1"~"Southern.shark_1_Zone1",
                               Fleet=="TDGDLF.monthly.Zone2"~"Southern.shark_1_Zone2",
                               Fleet=="TDGDLF.daily.West"~"Southern.shark_2_West",
                               Fleet=="TDGDLF.daily.Zone1"~"Southern.shark_2_Zone1",
                               Fleet=="TDGDLF.daily.Zone2"~"Southern.shark_2_Zone2",
                               TRUE~Fleet))%>%
        left_join(Flits.and.survey.zone,by=c('Fleet'='Fleet.name'))%>%
        mutate(Fleet=Fleet.number)%>%
        dplyr::select(-Fleet.number)
      Var.ad.factr.zone=Var.ad.factr_cpue.zone           
      if(exists("Var.ad.factr_meanbodywt.zone")) Var.ad.factr.zone=rbind(Var.ad.factr.zone,Var.ad.factr_meanbodywt.zone)
      clear.log("Var.ad.factr_meanbodywt.zone")
      clear.log("Var.ad.factr_cpue.zone")
      if(drop.intermediate.yrs)if(Whiskery.q.periods==2 & Neim=='whiskery shark') 
      {
        CPUE.zone$TDGDLF.monthly=CPUE.zone$TDGDLF.monthly%>%
          filter(!yr.f%in%Life.history$Yr_q_change_transition)
      }
    }
  }
  
  
  #5. Add size comp effective sample size bias adjustment    
  #note: a Value of 0 means no effect
  #zones together
  if(!is.null(Var.ad.factr))
  {
    tuned.siz.comp=Life.history$tuned_size_comp
    if(!is.null(tuned.siz.comp))Var.ad.factr=rbind(Var.ad.factr,tuned.siz.comp)
  }
  if(is.null(Var.ad.factr) &!is.null(Size.compo.SS.format))
  {
    tuned.siz.comp=Life.history$tuned_size_comp
    Var.ad.factr=tuned.siz.comp
  }
  #by zones  
  if(!is.null(Var.ad.factr.zone))
  {
    tuned.siz.comp=Life.history$tuned_size_comp.zone
    if(!is.null(tuned.siz.comp))Var.ad.factr.zone=rbind(Var.ad.factr.zone,tuned.siz.comp)
  }
  if(is.null(Var.ad.factr.zone) &!is.null(Size.compo.SS.format.zone))
  {
    tuned.siz.comp=Life.history$tuned_size_comp.zone
    Var.ad.factr.zone=tuned.siz.comp
  }
  
  
  #6. F from tagging studies on TDGDLF (1994-95 and 2001-03)
  F.SS.format=NULL  
  if(any(grepl('Fishing.mortality.TDGDLF',names(Species.data[[i]]))) & length(CPUE)>0)
  {
    F.SS.format=Species.data[[i]][[grep('Fishing.mortality.TDGDLF',names(Species.data[[i]]))]]%>%
      mutate(year=as.numeric(substr(Finyear,1,4)),
             month=1,
             Fleet='Southern.shark_1')%>%
      filter(year<=max(ktch$finyear))%>%
      dplyr::select(-Finyear)%>%
      left_join(Flits.and.survey,by=c('Fleet'='Fleet.name'))%>%
      mutate(fleet=Fleet.number,
             CV=0.05)%>%
      dplyr::select(-c(Fleet.number,Fleet))%>%
      relocate(year,month,fleet,Mean,CV)%>%
      arrange(year)
    
  }
  
  
  #7. Tagging data
  # zones together
  Tags.SS.format=NULL
  if(names(Species.data)[i]%in%use.tag.data)
  {
    #Extract data
    releases=Species.data[[i]]$Con_tag_SS.format_releases%>%
      filter(Yr.rel<=Last.yr.ktch.numeric)
    recaptures=Species.data[[i]]$Con_tag_SS.format_recaptures%>%
      filter(Yr.rec<=Last.yr.ktch.numeric)
    
    #Keep relevant finyears and zones  
    dis.yrs.tag=Use.these.tag.year_zones[[match(names(Species.data)[i],names(Use.these.tag.year_zones))]]
    releases=releases%>%
              mutate(dummy=paste(Yr.rel,Rel.zone))%>%
              filter(dummy%in%dis.yrs.tag)%>%
              dplyr::select(-dummy)
    recaptures=recaptures%>%filter(Tag.group%in%unique(releases$Tag.group))
    if(use.tag.rec.yrs.90percent.rec)
    {
      Table.yr.releases=table(releases$Yr.rel)
      vec=cumsum(Table.yr.releases)/sum(Table.yr.releases)
      Last.yr.rel=names(vec[which.min(abs(vec - 0.9))])
      
      Table.yr.recaptures=table(recaptures$Yr.rec)
      vec=cumsum(Table.yr.recaptures)/sum(Table.yr.recaptures)
      Last.yr.rec=names(vec[which.min(abs(vec - 0.9))])
      recaptures=recaptures%>%filter(Yr.rec<=as.numeric(Last.yr.rec))
    }

    #Allocate West, Zone1 and Zone 2 to Southern 1 or 2
    releases=releases%>%
      mutate(Rel.zone=case_when(Yr.rel<=2005 ~'Southern.shark_1',
                                Yr.rel>2005 ~'Southern.shark_2'))
    recaptures=recaptures%>%
      mutate(Rec.zone=case_when(Yr.rec<=2005 ~'Southern.shark_1',
                                Yr.rec>2005 ~'Southern.shark_2'))
    
    #Recalculate TagGroup
    releases=releases%>%
      arrange(Rel.zone,Yr.rel,Sex,Age)%>%
      mutate(rowid = row_number())
    TagGroup=releases%>%distinct(Tag.group,rowid)
    recaptures=recaptures%>%left_join(TagGroup,by='Tag.group')
    releases=releases%>%
      mutate(Tag.group=rowid)%>%
      dplyr::select(-rowid)
    recaptures=recaptures%>%
      mutate(Tag.group=rowid)%>%
      dplyr::select(-rowid)
    
    #group sex  
    if(taggroup.sex.combined)   
    {
      releases=releases%>%
        mutate(Sex=0)%>%
        group_by(Rel.zone,Yr.rel,season,t.fill,Sex,Age)%>%
        mutate(N.release=sum(N.release))%>%
        ungroup()%>%
        group_by(Rel.zone,Yr.rel,season,t.fill,Sex,Age) %>%
        mutate(Tag.groups = paste(Tag.group, collapse = ", "))%>%
        ungroup()
      recaptures=recaptures%>%   
        left_join(releases%>%
                    dplyr::select(Tag.group,Tag.groups),
                  by='Tag.group')%>%
        group_by(Yr.rec,season,Rec.zone,Tag.groups)%>%
        summarise(N.recapture=sum(N.recapture))%>%
        ungroup()
      releases=releases%>%
        distinct(Tag.groups,Rel.zone,Yr.rel,season,t.fill,Sex,Age,N.release)%>%
        arrange(Rel.zone,Yr.rel,season,t.fill,Sex,Age)%>%
        mutate(Tag.group = row_number())%>%
        relocate(Tag.group)
      recaptures=recaptures%>%
        left_join(releases%>%distinct(Tag.group,Tag.groups),
                  by='Tag.groups')%>%
        dplyr::select(-Tag.groups)%>%
        relocate(Tag.group)
      releases=releases%>%dplyr::select(-Tag.groups)
    }
    
    releases=releases%>%
      rename(Area=Rel.zone)%>%
      mutate(Area=1)
    get.fleet=recaptures%>%
      distinct(Yr.rec,Rec.zone)
    Rec.ZonEs=unique(recaptures$Rec.zone)
    a1=ktch%>%
      ungroup()%>%
      dplyr::select(-c(SPECIES,Name))%>%
      gather(Fleet,Ktch,-finyear)%>%
      mutate(zone=case_when(Fleet=="Northern.shark"~'North',
                            grepl("Southern.shark_1",Fleet)~"Southern.shark_1",
                            grepl("Southern.shark_2",Fleet)~"Southern.shark_2",
                            TRUE~''))%>%
      filter(Ktch>0)%>%
      filter(zone%in%unique(get.fleet$Rec.zone))%>%
      filter(finyear%in%unique(get.fleet$Yr.rec))%>%
      distinct(finyear,Fleet,zone)%>%
      left_join(data.frame(Fleet.ID=Flits,Fleet=names(Flits)),
                by='Fleet')
    get.fleet=get.fleet%>%
      left_join(a1,by=c('Rec.zone'='zone','Yr.rec'='finyear'))
    recaptures=recaptures%>%
      left_join(get.fleet%>%dplyr::select(-Fleet)%>%rename(Fleet=Fleet.ID),
                by=c('Rec.zone','Yr.rec'))%>%
      dplyr::select(-Rec.zone)%>%
      relocate(Tag.group,Yr.rec,season,Fleet,N.recapture)
    
    Initial.tag.loss=1e-8  #tag-induced mortality immediately after tagging 
    Chronic.tag.loss=Species.data[[i]]$Con_tag_shedding_from_F.estimation.R_$x  #annual rate of tag loss; McAuley et al 2007 tag shedding
    
    
    if(Reporting.rate.type[[Neim]]=='published') Initial.reporting.rate=Species.data[[i]]$Con_tag_non_reporting_from_F.estimation.R_   #NEW
    if(Reporting.rate.type[[Neim]]=='calculated') Initial.reporting.rate=Species.data[[i]]$Con_tag_non_reporting_from_F.estimation.R_calculated
    Initial.reporting.rate=Initial.reporting.rate%>%
      dplyr::select(-Species)%>%
      gather(Zone,Non.reporting,-Finyear)%>%
      mutate(Reporting=1-Non.reporting,
             Zone=case_when(Zone%in%c('South','South.west','West') & Finyear<=2005 ~'Southern.shark_1',
                            Zone%in%c('South','South.west','West') & Finyear>2005 ~'Southern.shark_2',
                            Zone=='North'~'North'))%>%
      filter(Zone%in%Rec.ZonEs)%>%
      group_by(Finyear,Zone)%>%
      summarise(Non.reporting=mean(Non.reporting,na.rm=T),
                Reporting=mean(Reporting,na.rm=T))%>%
      ungroup()%>%
      mutate(Reporting.logit=fn.inv.logit(Reporting))%>%
      dplyr::select(-Non.reporting)
    
    get.fleet1=get.fleet%>%
      dplyr::select(-Fleet)%>%
      rename(Fleet=Fleet.ID)%>%
      mutate(dummy=paste(Yr.rec,Rec.zone))
    not.in.init.rep=paste(Initial.reporting.rate$Finyear,Initial.reporting.rate$Zone)
    not.in.init.rep=not.in.init.rep[which(!not.in.init.rep%in%get.fleet1$dummy)]
    if(length(not.in.init.rep)>0)
    {
      ad.get.flit=get.fleet1[1:length(not.in.init.rep),]%>%
        mutate(dummy=not.in.init.rep,
               Yr.rec=word(dummy, 1),
               Rec.zone=word(dummy, 2),
               Fleet=NA)
      get.fleet1=rbind(get.fleet1,ad.get.flit)%>%
        arrange(Rec.zone,Yr.rec)%>%
        fill(Fleet, .direction = "down") #NEW
      
      get.fleet1=get.fleet1%>%  
        left_join(a1%>%
                    mutate(finyear=as.character(finyear))%>%
                    distinct(finyear,zone,Fleet.ID),
                  by=c('Yr.rec'='finyear','Rec.zone'='zone'))%>%
        mutate(Fleet=case_when(!is.na(Fleet.ID) & !Fleet==Fleet.ID~Fleet.ID,
                               !is.na(Fleet.ID) & is.na(Fleet)~Fleet.ID,
                               TRUE~Fleet))%>%
        dplyr::select(-Fleet.ID)
    }
    get.fleet1=get.fleet1%>%dplyr::select(-dummy)%>%mutate(Yr.rec=as.numeric(Yr.rec))
    Initial.reporting.rate=Initial.reporting.rate%>%
      left_join(get.fleet1,
                by=c('Zone'='Rec.zone','Finyear'='Yr.rec'))%>%
                filter(!is.na(Reporting))
    if(estimate.tag.report.decay)
    {
      rep.dec.flit=sort(unique(Initial.reporting.rate$Fleet))
      Rep.decay=vector('list',length(rep.dec.flit))
      Rep.decay_p=Rep.decay
      for(re in 1:length(Rep.decay))
      {
        d.init.rep=Initial.reporting.rate%>%
          filter(Fleet==rep.dec.flit[re])%>%
          arrange(Finyear)%>%
          mutate(time=Finyear-Finyear[1])
        diKay=0  #NEW
        if(nrow(d.init.rep)>1)#NEW
        {
          Init.rep=d.init.rep$Reporting[1]
          fit_nls <- nls(Reporting ~ Init.rep * exp(-k * time), 
                         data = d.init.rep, 
                         start = list(k = 0.01))
          diKay=round(coef(fit_nls),4)              #NEW
          if(!allow.increase.tag.rep.rate) diKay=max(0,diKay) #NEW
        }
        
        Rep.decay[[re]]=data.frame(Fleet=rep.dec.flit[re],decay=diKay)
        Rep.decay_p[[re]]=ggplot(data=d.init.rep,aes(Finyear,Reporting))+
          geom_point(size=3)+
          ylim(0,1)+
          geom_line(data=data.frame(Finyear=d.init.rep$Finyear,
                                    Reporting=fn.SS3.tag.reporting.rate(init.rep.rate=d.init.rep$Reporting[1],exp.decay.rate=diKay,time=d.init.rep$time)),
                    aes(Finyear,Reporting),color=2,linewidth = 2)+
          theme_PA()+
          ggtitle(paste0('Fleet ',rep.dec.flit[re],' (',unique(d.init.rep$Zone),' decay=',diKay,')'))
        
      }
      Reporting.rate.decay=do.call(rbind,Rep.decay)
      if(First.run=='YES')
      {
        ggarrange(plotlist=Rep.decay_p,ncol=1,nrow=length(Rep.decay_p))
        ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(List.sp[[i]]$Name),"/",AssessYr,"/1_Inputs/Visualise data/Tagging_report rate decay_single.area.tiff",sep=''),
               width = 6,height = 8,compression = "lzw")
      }
      if(pass.rep.rate.decay.negative)
      {
        Reporting.rate.decay=Reporting.rate.decay%>%
                              mutate(decay=ifelse(decay>0,-decay,abs(decay)))
      }
    }
    if(!estimate.tag.report.decay) Reporting.rate.decay=0 #Andre's Gummy and (Spatial SS3 workshop) Lecture D '4 areas' models
    
    if(logit.transform.tag.pars)   
    {
      Initial.tag.loss=fn.inv.logit(Initial.tag.loss)
      Chronic.tag.loss=fn.inv.logit(Chronic.tag.loss)
      Initial.reporting.rate=Initial.reporting.rate%>%
                              dplyr::select(-Reporting)%>%
                              rename(Reporting=Reporting.logit)
    }
    
    Tags.SS.format=list(
      releases=releases%>%data.frame,
      recaptures=recaptures%>%data.frame,
      Initial.tag.loss=Initial.tag.loss,   
      Chronic.tag.loss=Chronic.tag.loss,   
      Initial.reporting.rate=Initial.reporting.rate%>%
                              filter(Finyear==Initial.reporting.rate$Finyear[which.min(abs(Initial.reporting.rate$Finyear - min(releases$Yr.rel)))]),   
      Reporting.rate.decay=Reporting.rate.decay,
      overdispersion=SS_overdispersion,       # Andre's Gummy model
      mixing_latency_period=SS_mixing_latency_period, 
      max_periods=ceiling((max(recaptures$Yr.rec)-min(releases$Yr.rel))*Extend.mx.period))  # 30 Andre's Gummy model  
    
    rm(releases,recaptures,Chronic.tag.loss,Initial.reporting.rate,Reporting.rate.decay)
  }
  
  # by zones
  Tags.SS.format.zone=NULL
  if(names(Species.data)[i]%in%use.tag.data)
  {
    #Extract data
    releases=Species.data[[i]]$Con_tag_SS.format_releases%>%
      filter(Yr.rel<=Last.yr.ktch.numeric)
    recaptures=Species.data[[i]]$Con_tag_SS.format_recaptures%>%
      filter(Yr.rec<=Last.yr.ktch.numeric)
    
    #Keep relevant finyear zones   
    dis.yrs.tag=Use.these.tag.year_zones[[match(names(Species.data)[i],names(Use.these.tag.year_zones))]]
    releases=releases%>%
              mutate(dummy=paste(Yr.rel,Rel.zone))%>%
              filter(dummy%in%dis.yrs.tag)%>%
              dplyr::select(-dummy)
    recaptures=recaptures%>%filter(Tag.group%in%unique(releases$Tag.group))
    if(use.tag.rec.yrs.90percent.rec)
    {
      Table.yr.releases=table(releases$Yr.rel)
      vec=cumsum(Table.yr.releases)/sum(Table.yr.releases)
      Last.yr.rel=names(vec[which.min(abs(vec - 0.9))])
      
      Table.yr.recaptures=table(recaptures$Yr.rec)
      vec=cumsum(Table.yr.recaptures)/sum(Table.yr.recaptures)
      Last.yr.rec=names(vec[which.min(abs(vec - 0.9))])
      recaptures=recaptures%>%filter(Yr.rec<=as.numeric(Last.yr.rec))
    }

    #Recalculate TagGroup
    releases=releases%>%
      arrange(Rel.zone,Yr.rel,Sex,Age)%>%
      mutate(rowid = row_number())
    TagGroup=releases%>%distinct(Tag.group,rowid)
    recaptures=recaptures%>%left_join(TagGroup,by='Tag.group')
    releases=releases%>%
      mutate(Tag.group=rowid)%>%
      dplyr::select(-rowid)
    recaptures=recaptures%>%
      mutate(Tag.group=rowid)%>%
      dplyr::select(-rowid)
    
    #group sex  
    if(taggroup.sex.combined)   
    {
      releases=releases%>%
        mutate(Sex=0)%>%
        group_by(Rel.zone,Yr.rel,season,t.fill,Sex,Age)%>%
        mutate(N.release=sum(N.release))%>%
        ungroup()%>%
        group_by(Rel.zone,Yr.rel,season,t.fill,Sex,Age) %>%
        mutate(Tag.groups = paste(Tag.group, collapse = ", "))%>%
        ungroup()
      recaptures=recaptures%>%   
        left_join(releases%>%
                    dplyr::select(Tag.group,Tag.groups),
                  by='Tag.group')%>%
        group_by(Yr.rec,season,Rec.zone,Tag.groups)%>%
        summarise(N.recapture=sum(N.recapture))%>%
        ungroup()
      releases=releases%>%
        distinct(Tag.groups,Rel.zone,Yr.rel,season,t.fill,Sex,Age,N.release)%>%
        arrange(Rel.zone,Yr.rel,season,t.fill,Sex,Age)%>%
        mutate(Tag.group = row_number())%>%
        relocate(Tag.group)
      recaptures=recaptures%>%
        left_join(releases%>%distinct(Tag.group,Tag.groups),
                  by='Tag.groups')%>%
        dplyr::select(-Tag.groups)%>%
        relocate(Tag.group)
      releases=releases%>%dplyr::select(-Tag.groups)
    }
    
    releases=releases%>%
      rename(Area=Rel.zone)%>%
      mutate(Area=1)
    get.fleet=recaptures%>%
      distinct(Yr.rec,Rec.zone)
    Rec.ZonEs=unique(recaptures$Rec.zone)
    a1=ktch.zone%>%
      ungroup()%>%
      dplyr::select(-c(SPECIES,Name))%>%
      gather(Fleet,Ktch,-finyear)%>%
      mutate(zone=case_when(Fleet=="Northern.shark"~'North',
                            grepl("Zone1",Fleet)~"Zone1",
                            grepl("West",Fleet)~"West",
                            grepl("Zone2",Fleet)~"Zone2",
                            TRUE~''))%>%
      filter(Ktch>0)%>%
      filter(zone%in%unique(get.fleet$Rec.zone))%>%
      filter(finyear%in%unique(get.fleet$Yr.rec))%>%
      distinct(finyear,Fleet,zone)%>%
      left_join(data.frame(Fleet.ID=Flits.zone,Fleet=names(Flits.zone)),
                by='Fleet')
    get.fleet=get.fleet%>%
      left_join(a1,by=c('Rec.zone'='zone','Yr.rec'='finyear'))
    recaptures=recaptures%>%
      left_join(get.fleet%>%dplyr::select(-Fleet)%>%rename(Fleet=Fleet.ID),
                by=c('Rec.zone','Yr.rec'))%>%
      dplyr::select(-Rec.zone)%>%
      relocate(Tag.group,Yr.rec,season,Fleet,N.recapture)
    
    Initial.tag.loss=1e-4  #tag-induced mortality immediately after tagging 
    Chronic.tag.loss=Species.data[[i]]$Con_tag_shedding_from_F.estimation.R_$x  #annual rate of tag loss; McAuley et al 2007 tag shedding
    
    if(Reporting.rate.type[[Neim]]=='published') Initial.reporting.rate=Species.data[[i]]$Con_tag_non_reporting_from_F.estimation.R_   #NEW
    if(Reporting.rate.type[[Neim]]=='calculated') Initial.reporting.rate=Species.data[[i]]$Con_tag_non_reporting_from_F.estimation.R_calculated
    Initial.reporting.rate=Initial.reporting.rate%>%
      dplyr::select(-Species)%>%
      gather(Zone,Non.reporting,-Finyear)%>%
      mutate(Reporting=1-Non.reporting,
             Zone=case_when(Zone=='South'~'Zone2',
                            Zone=='South.west'~'Zone1',
                            Zone=='West'~'West',
                            Zone=='North'~'North'))%>%
      filter(Zone%in%Rec.ZonEs)%>%
      mutate(Reporting.logit=fn.inv.logit(Reporting))%>%
      dplyr::select(-Non.reporting)
    get.fleet1=get.fleet%>%
      dplyr::select(-Fleet)%>%
      rename(Fleet=Fleet.ID)%>%
      mutate(dummy=paste(Yr.rec,Rec.zone))
    not.in.init.rep=paste(Initial.reporting.rate$Finyear,Initial.reporting.rate$Zone)
    not.in.init.rep=not.in.init.rep[which(!not.in.init.rep%in%get.fleet1$dummy)]
    if(length(not.in.init.rep)>0)
    {
      ad.get.flit=get.fleet1[1:length(not.in.init.rep),]%>%
        mutate(dummy=not.in.init.rep,
               Yr.rec=word(dummy, 1),
               Rec.zone=word(dummy, 2),
               Fleet=NA)
      get.fleet1=rbind(get.fleet1,ad.get.flit)%>%
        arrange(Rec.zone,Yr.rec)%>%
        fill(Fleet, .direction = "down") 
      
      get.fleet1=get.fleet1%>%  
        left_join(a1%>%
                    mutate(finyear=as.character(finyear))%>%
                    distinct(finyear,zone,Fleet.ID),
                  by=c('Yr.rec'='finyear','Rec.zone'='zone'))%>%
        mutate(Fleet=case_when(!is.na(Fleet.ID) & !Fleet==Fleet.ID~Fleet.ID,
                               !is.na(Fleet.ID) & is.na(Fleet)~Fleet.ID,
                               TRUE~Fleet))%>%#NEW
        dplyr::select(-Fleet.ID)
    }
    get.fleet1=get.fleet1%>%dplyr::select(-dummy)%>%mutate(Yr.rec=as.numeric(Yr.rec))
    Initial.reporting.rate=Initial.reporting.rate%>%
      left_join(get.fleet1,
                by=c('Zone'='Rec.zone','Finyear'='Yr.rec'))%>%
      filter(!is.na(Reporting))
    
    if(estimate.tag.report.decay)
    {
      rep.dec.flit=sort(unique(Initial.reporting.rate$Fleet))
      Rep.decay=vector('list',length(rep.dec.flit))
      Rep.decay_p=Rep.decay
      for(re in 1:length(Rep.decay))
      {
        d.init.rep=Initial.reporting.rate%>%
          filter(Fleet==rep.dec.flit[re])%>%
          arrange(Finyear)%>%
          mutate(time=Finyear-Finyear[1])
        diKay=0  #NEW
        if(nrow(d.init.rep)>1)#NEW
        {
          Init.rep=d.init.rep$Reporting[1]
          fit_nls <- nls(Reporting ~ Init.rep * exp(-k * time), 
                         data = d.init.rep, 
                         start = list(k = 0.01))
          diKay=round(coef(fit_nls),4)              #NEW
          if(!allow.increase.tag.rep.rate) diKay=max(0,diKay) #NEW
        }
        
        Rep.decay[[re]]=data.frame(Fleet=rep.dec.flit[re],decay=diKay)
        Rep.decay_p[[re]]=ggplot(data=d.init.rep,aes(Finyear,Reporting))+
          geom_point(size=3)+
          ylim(0,1)+
          geom_line(data=data.frame(Finyear=d.init.rep$Finyear,
                                    Reporting=fn.SS3.tag.reporting.rate(init.rep.rate=d.init.rep$Reporting[1],exp.decay.rate=diKay,time=d.init.rep$time)),
                    aes(Finyear,Reporting),color=2,linewidth = 2)+
          theme_PA()+
          ggtitle(paste0('Fleet ',rep.dec.flit[re],' (',unique(d.init.rep$Zone),' decay=',diKay,')'))
        
      }
      Reporting.rate.decay=do.call(rbind,Rep.decay)
      if(First.run=='YES')
      {
        ggarrange(plotlist=Rep.decay_p,ncol=1,nrow=length(Rep.decay_p))
        ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(List.sp[[i]]$Name),"/",AssessYr,"/1_Inputs/Visualise data/Tagging_report rate decay_area.as.fleets.tiff",sep=''),
               width = 6,height = 8,compression = "lzw")
      }
      if(pass.rep.rate.decay.negative)
      {
        Reporting.rate.decay=Reporting.rate.decay%>%
          mutate(decay=ifelse(decay>0,-decay,abs(decay)))
      }
    }
    if(!estimate.tag.report.decay) Reporting.rate.decay=0 #Andre's Gummy and (Spatial SS3 workshop) Lecture D '4 areas' models
    if(logit.transform.tag.pars)   
    {
      Initial.tag.loss=fn.inv.logit(Initial.tag.loss)
      Chronic.tag.loss=fn.inv.logit(Chronic.tag.loss)
      Initial.reporting.rate=Initial.reporting.rate%>%
        dplyr::select(-Reporting)%>%
        rename(Reporting=Reporting.logit)
    }
    Tags.SS.format.zone=list(
                        releases=releases%>%data.frame,
                        recaptures=recaptures%>%data.frame,
                        Initial.tag.loss=Initial.tag.loss,   
                        Chronic.tag.loss=Chronic.tag.loss,   
                        Initial.reporting.rate=Initial.reporting.rate%>%
                          filter(Finyear==Initial.reporting.rate$Finyear[which.min(abs(Initial.reporting.rate$Finyear - min(releases$Yr.rel)))]),
                        Reporting.rate.decay=Reporting.rate.decay,
                        overdispersion=SS_overdispersion,       # Andre's Gummy model
                        mixing_latency_period=SS_mixing_latency_period,  
                        max_periods=ceiling((max(recaptures$Yr.rec)-min(releases$Yr.rel))*Extend.mx.period))  # 30 Andre's Gummy model  
    
    rm(releases,recaptures,Chronic.tag.loss,Initial.reporting.rate,Reporting.rate.decay)
  }
  
  
  #8. Conditional age at length
  Cond.age.len.SS.format=NULL
  if(do.Cond.age.len.SS.format)
  {
    if(any(grepl('conditional_age_length',names(Species.data[[i]]))))
    {
      a=Life.history$a_FL.to.TL
      b=Life.history$b_FL.to.TL
      Cond.age.len.SS.format=Species.data[[i]]$conditional_age_length%>%
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
  
  
  #9. MeanSize at Age obs
  MeanSize.at.Age.obs.SS.format=NULL
  if(any(grepl('conditional_age_length',names(Species.data[[i]]))) & names(Species.data)[i]%in%Mean.Size.at.age.species)
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
    
    a=Life.history$a_FL.to.TL
    b=Life.history$b_FL.to.TL
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
  
  
  #10. Fleet info
  #zones together
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
  #by zone
  flitinfo.zone=data.frame(fleet=Flits.zone)%>%
    mutate(type=1,
           surveytiming=-1,
           area=1,
           units=1,
           need_catch_mult=0,
           fleetname=names(Flits.zone))%>%
    dplyr::select(-fleet)%>%
    mutate(type=ifelse(fleetname=="Survey",3,type),
           surveytiming=ifelse(fleetname=="Survey",1,surveytiming))
  rownames(flitinfo.zone)=NULL
  
  
  #11. Run scenarios if available abundance index or size comps
  len.cpue=length(CPUE)
  len.cpue.zone=length(CPUE.zone)
  len.size.comp=length(Size.compo.SS.format) 
  if(len.cpue>0 | len.cpue.zone>0 | len.size.comp>0)
  {
    #Put CPUE in SS format
    if(len.cpue>0 | len.cpue.zone>0)
    {
      MAX.CV=Life.history$MAX.CV
      
      #zones together
      if(len.cpue>0)
      {
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
        Abundance.SS.format=do.call(rbind,CPUE)%>%
          relocate(Year,seas,index,Mean,CV)%>%
          arrange(index,Year)
        
      }
       
      #by zone
      if(len.cpue.zone>0)
      {
        for(x in 1:len.cpue.zone)    
        {
          nm=names(CPUE.zone)[x]
          if(nm=="NSF") nm="Northern.shark"
          if(nm=="TDGDLF.monthly.West") nm="Southern.shark_1_West"
          if(nm=="TDGDLF.monthly.Zone1") nm="Southern.shark_1_Zone1"
          if(nm=="TDGDLF.monthly.Zone2") nm="Southern.shark_1_Zone2"
          if(nm=="TDGDLF.daily.West") nm="Southern.shark_2_West"
          if(nm=="TDGDLF.daily.Zone1") nm="Southern.shark_2_Zone1"
          if(nm=="TDGDLF.daily.Zone2") nm="Southern.shark_2_Zone2"
          dd=CPUE.zone[[x]][,grep(paste(c('yr.f','Mean','MeAn','CV'),collapse="|"),names(CPUE.zone[[x]]))]%>%
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
            left_join(Flits.and.survey.zone,by=c('index.dummy'='Fleet.name'))%>%
            mutate(index=Fleet.number)%>%
            dplyr::select(-c(index.dummy,Fleet.number))
          CPUE.zone[[x]]=dd%>%filter(!is.na(Mean))
        }
        Abundance.SS.format.zone=do.call(rbind,CPUE.zone)%>%
          relocate(Year,seas,index,Mean,CV)%>%
          arrange(index,Year)
      }

    }
    
    #Scenarios
    Scens=Life.history$Sens.test$SS
    if(!do.all.sensitivity.tests) Scens=Scens%>%filter(Scenario=='S1')
    Scens=Scens%>%
      mutate(Species=capitalize(Neim))
    Store.sens=vector('list',nrow(Scens))
    names(Store.sens)=Scens$Scenario
    Out.Scens=Scens
    Out.estimates=Out.likelihoods=Out.quantities=Out.rel.biom=Out.probs.rel.biom=Out.f.series=Out.probs.f.series=Out.B.Bmsy=
      Out.F.Fmsy=Out.probs.B.Bmsy=Out.Kobe.probs=store.warnings=store.convergence=vector('list',length(Store.sens))
    
    #Life history
    Life.history$Fecundity=ceiling(mean(Life.history$Fecundity))
    if(Neim%in%species.constant.fec)
    {
      Life.history$Fecu_a=NA
      Life.history$Fecu_b=NA
    }
    if(Neim%in%species.increase.terminal.age) Life.history$Max.age.F=max(Life.history$Max.age.F)
    Life.history$Max.age.F=ceiling(mean(Life.history$Max.age.F))
    Life.history$Breed.cycle=mean(Life.history$Breed.cycle)
    
    #Likelihood lambdas
    Lamdas.SS.lambdas=Life.history$SS_lambdas
    if(!is.null(Lamdas.SS.lambdas))
    {
      Lamdas.SS.lambdas$fleet=match(Lamdas.SS.lambdas$fleet,flitinfo$fleetname)
    }
    
    #Rename fleets following SS nomenclature
    names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))]=match(names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))],names(Flits))
    names(ktch.zone)[which(!names(ktch.zone)%in%c("SPECIES","Name","finyear"))]=match(names(ktch.zone)[which(!names(ktch.zone)%in%c("SPECIES","Name","finyear"))],names(Flits.zone))
    
    #Future catches/F  
    if("SS"%in%future.models)
    {
      #zones together
      NN=nrow(ktch)
      add.ct.future=ktch[1:years.futures,]
      add.ct.future$finyear=ktch$finyear[NN]+(1:years.futures)
      if("Survey"%in%names(Flits))
      {
        lef.flits=match(Flits[-match("Survey",names(Flits))],colnames(add.ct.future))
      }else
      {
        lef.flits=match(Flits,colnames(add.ct.future))
      }
      for(lf in 1:length(lef.flits))
      {
        add.ct.future[,lef.flits[lf]]=mean(unlist(ktch[(NN-years.futures+1):NN,lef.flits[lf]]),na.rm=T)
      }
      
      #by zone
      NN=nrow(ktch.zone)
      add.ct.future.zone=ktch.zone[1:years.futures,]
      add.ct.future.zone$finyear=ktch.zone$finyear[NN]+(1:years.futures)
      if("Survey"%in%names(Flits.zone))
      {
        lef.flits.zone=match(Flits.zone[-match("Survey",names(Flits.zone))],colnames(add.ct.future.zone))
      }else
      {
        lef.flits.zone=match(Flits.zone,colnames(add.ct.future.zone))
      }
      for(lf in 1:length(lef.flits.zone))
      {
        add.ct.future.zone[,lef.flits.zone[lf]]=mean(unlist(ktch.zone[(NN-years.futures+1):NN,lef.flits.zone[lf]]),na.rm=T)
      }
    }
    
    #Execute SS 
    for(s in 1:length(Store.sens))
    {
      this.wd1=paste(this.wd,names(Store.sens)[s],sep='/')
      if(!dir.exists(this.wd1))dir.create(this.wd1)
      
      #cpues
      Abund.single.area=Abundance.SS.format
      Abund.areas.as.fleets=Abundance.SS.format.zone
      
      #remove monthly years  
      if('drop.monthly.cpue.min'%in%names(Scens))
      {
        if(!is.na(Scens$drop.monthly.cpue.min[s]))
        {
          drop.yrS=Scens$drop.monthly.cpue.min[s]:Scens$drop.monthly.cpue.max[s]
          #zones combined
          if(!is.null(Abund.single.area))
          {
            Abund.single.area=Abund.single.area%>%
              mutate(dummy = rownames(Abund.single.area))%>%
              filter(!(grepl('TDGDLF.monthly',dummy) & Year%in%drop.yrS))%>%
              dplyr::select(-dummy)
          }
          #by zone
          if(!is.null(Abund.areas.as.fleets))
          {
            Abund.areas.as.fleets=Abund.areas.as.fleets%>%
              mutate(dummy = rownames(Abund.areas.as.fleets))%>%
              filter(!(grepl('TDGDLF.monthly',dummy) & Year%in%drop.yrS))%>%
              dplyr::select(-dummy)
          }
        }
      }
      
      
      #remove daily years  
      if(!is.na(Scens$Daily.cpues[s]))
      {
        if('TDGDLF.daily'%in%names(CPUE) | any(grep('TDGDLF.daily',names(CPUE.zone))))
        {
          rid.of=as.numeric(unlist(str_split(Scens$Daily.cpues[s], "&")))
          
          #zones combined
          if(!is.null(Abund.single.area))
          {
            drop.dis=Abund.single.area%>%
              mutate(this=grepl('TDGDLF.daily',rownames(Abund.single.area)) & Year%in%rid.of)
            if(any(drop.dis$this)) Abund.single.area=Abund.single.area[-which(drop.dis$this),]
          }
          
          #by zone
          if(!is.null(Abund.areas.as.fleets))
          {
            drop.dis=Abund.areas.as.fleets%>%
              mutate(this=grepl('TDGDLF.daily',rownames(Abund.areas.as.fleets)) & Year%in%rid.of)
            if(any(drop.dis$this)) Abund.areas.as.fleets=Abund.areas.as.fleets[-which(drop.dis$this),]
          }
        }
      }
      
      #remove TDGDLF CPUE 
      if(Scens$CPUE[s]=='No') 
      {
        Abund.single.area=Abund.areas.as.fleets=NULL
      }
      if('test.use.cpue'%in%names(Scens[s,]))
      {
        if(!(Scens$test.use.cpue[s]=='Yes'))
        {
          if(any(grepl("TDGDLF",names(CPUE))))
          {
            drop.dis.rows=grep("TDGDLF",row.names(Abund.single.area))
            if(length(drop.dis.rows)>0) Abund.single.area=Abund.single.area[-drop.dis.rows,]
          }
          if(any(grepl("TDGDLF",names(CPUE.zone))))
          {
            drop.dis.rows=grep("TDGDLF",row.names(Abund.areas.as.fleets))
            if(length(drop.dis.rows)>0) Abund.areas.as.fleets=Abund.areas.as.fleets[-drop.dis.rows,] 
          }
        }
        if(!is.null(Abund.single.area)) if(nrow(Abund.single.area)==0) Abund.single.area=NULL
        if(!is.null(Abund.areas.as.fleets)) if(nrow(Abund.areas.as.fleets)==0) Abund.areas.as.fleets=NULL
        
      }
      
      #remove Length.comps
      Size.comp.single.area=Size.compo.SS.format
      Size.comp.areas.as.fleets=Size.compo.SS.format.zone
      if(Scens$Length.comps[s]=='No') 
      {
        Size.comp.single.area=Size.comp.areas.as.fleets=NULL
      }
      
      #remove Mean.body
      Meanbodywt.single.area=meanbodywt.SS.format
      Meanbodywt.areas.as.fleets=meanbodywt.SS.format.zone
      if(Scens$Mean.body[s]=='No') 
      {
        Meanbodywt.single.area=Meanbodywt.areas.as.fleets=NULL
      }
      
      #use cpue or length comp in likelihood?
        #note: superseded
      if(Life.history$drop.length.comp)
      {
        Size.compo.SS.format=NULL
        Size.compo.SS.format.zone=NULL
      }
      if(Life.history$drop.cpue)
      {
        Abundance.SS.format=NULL
        Abundance.SS.format.zone=NULL
      }
      
      #Set future F instead of catch if required by scenario  
      if(exists('add.ct.future'))   
      {
        #zones combined
        add.ct.or.F_future=add.ct.future
        if(Scens[s,'Forecasting']=="F")
        {
          Nms=names(add.ct.or.F_future)
          mutate.these.fleets=subset(Nms,!Nms%in%c("SPECIES","Name","finyear"))
          id.mutate.these.fleets=which(!Nms%in%c("SPECIES","Name","finyear"))
          for(mu in 1:length(mutate.these.fleets))
          {
            id.mu=match(names(Life.history$F.forecasting.value)[mu],flitinfo$fleetname)
            F.val=Life.history$F.forecasting.value[id.mu]
            mu.id=match(id.mu,Nms)
            add.ct.or.F_future[,mu.id]=ifelse(add.ct.or.F_future[,mu.id]>0,F.val,0)
          }
        }
        
        #by zone 
        add.ct.or.F_future.zone=add.ct.future.zone
        if(Scens[s,'Forecasting']=="F")
        {
          Nms=names(add.ct.or.F_future.zone)
          mutate.these.fleets=subset(Nms,!Nms%in%c("SPECIES","Name","finyear"))
          id.mutate.these.fleets=which(!Nms%in%c("SPECIES","Name","finyear"))
          for(mu in 1:length(mutate.these.fleets))
          {
            id.mu=match(names(Life.history$F.forecasting.value)[mu],flitinfo$fleetname)
            F.val=Life.history$F.forecasting.value[id.mu]
            mu.id=match(id.mu,Nms)
            add.ct.or.F_future.zone[,mu.id]=ifelse(add.ct.or.F_future.zone[,mu.id]>0,F.val,0)
          }
        }
      }
      
      #a. Create SS input files   
      if(create.SS.inputs)
      {
        #a.1 Specify ktch,flitinfo,abundance,comps,mean.weight,var.adjs & future based on scenario   
        if(Scens$Spatial[s]=='single area')  
        {
          KAtch=ktch
          FLitinFO=flitinfo
          Abund=Abund.single.area
          Size.com=Size.comp.single.area
          Size.com_all=dummy.Size.compo.SS.format.all
          meanbody=Meanbodywt.single.area
          tags=Tags.SS.format 
          Var.ad=Var.ad.factr
          add.future=add.ct.or.F_future
          KAtch.ret.disc=NULL
          if(!is.null(retained.discarded.ktch))   
          {
            KAtch.ret.disc=retained.discarded.ktch%>%
                            ungroup()%>%
                            left_join(FLitinFO%>%dplyr::select(fleetname)%>%mutate(fleet=1:nrow(FLitinFO)),
                                      by=c('Fishry'='fleetname'))%>%
                            dplyr::select(-c('SPECIES','Name','Fishry'))%>%
                            mutate(seas=1,stderr=CV.discards)%>%
                            rename(yr=finyear,
                                   obs=Tonnes)%>%
                            relocate(yr,seas,fleet,obs,stderr)%>%
                            arrange(fleet,yr)%>%
                            data.frame
            Length.limit=Life.history$SS_retention$P_5
            discard.flits=unique(KAtch.ret.disc$fleet)
            
            if(!is.null(Size.com))
            {
              id.rel.cols=names(Size.com)[grep(paste(c('f','m'),collapse = '|'),names(Size.com))]
              id.rel.cols=unique(as.numeric(gsub("[a-zA-Z]", "", subset(id.rel.cols,!id.rel.cols%in%c("Nsamp")))))
              from.to=seq(id.rel.cols[which.min(abs(id.rel.cols - Length.limit))],max(id.rel.cols),by=TL.bins.cm)
              id.rel.cols=grep(paste(from.to,collapse = '|'),names(Size.com))
              id.irrel.cols=grep(paste(names(Size.com)[-c(which(c("year","Seas","Fleet","Sex","Part","Nsamp")%in%
                                                                  names(Size.com)),id.rel.cols)],collapse = '|'),
                                 names(Size.com))
              Size.com.discards=Size.com%>%filter(Fleet%in%discard.flits)
              Size.com=Size.com%>%filter(!Fleet%in%discard.flits)
              Size.com.discards_retained=Size.com.discards
              Size.com.discards_discarded=Size.com.discards
              Size.com.discards_retained[,id.rel.cols]=0
              Size.com.discards_retained$Part=2
              Size.com.discards_discarded[,id.irrel.cols]=0
              Size.com.discards_discarded$Part=1
              Size.com.discards=Size.com.discards_retained
              drop.yrs=rowSums(Size.com.discards_discarded[,id.rel.cols])
              drop.yrs=which(drop.yrs<Min.annual.obs.ktch)
              if(length(drop.yrs)>0)
              {
                Size.com.discards_discarded=Size.com.discards_discarded[-drop.yrs,]
              }
              if(nrow(Size.com.discards_discarded)>0)
              {
                Size.com.discards=rbind(Size.com.discards_retained,Size.com.discards_discarded)
              }
              Size.com=rbind(Size.com,Size.com.discards)%>%
                arrange(Fleet,year)
            }
              
            if(!is.null(meanbody))
            {
              meanbody=meanbody%>%
                mutate(part=ifelse(fleet%in%discard.flits,2,part)) 
            }
          }
        }
        if(Scens$Spatial[s]=='areas-as-fleets')   
        {
          KAtch=ktch.zone
          FLitinFO=flitinfo.zone
          Abund=Abund.areas.as.fleets
          Size.com=Size.comp.areas.as.fleets
          Size.com_all=dummy.Size.compo.SS.format.all_zone
          meanbody=Meanbodywt.areas.as.fleets
          tags=Tags.SS.format.zone  
          Var.ad=Var.ad.factr.zone
          add.future=add.ct.or.F_future.zone
          KAtch.ret.disc=NULL
          if(!is.null(retained.discarded.ktch.zone))
          {
            KAtch.ret.disc=retained.discarded.ktch.zone%>%
              ungroup()%>%
              left_join(FLitinFO%>%dplyr::select(fleetname)%>%mutate(fleet=1:nrow(FLitinFO)),
                        by=c('Fishry'='fleetname'))%>%
              dplyr::select(-c('SPECIES','Name','Fishry'))%>%
              mutate(seas=1,stderr=CV.discards)%>% 
              rename(yr=finyear,
                     obs=Tonnes)%>%
              relocate(yr,seas,fleet,obs,stderr)%>%
              arrange(fleet,yr)%>%
              data.frame

            Length.limit=Life.history$SS_retention$P_5
            discard.flits=unique(KAtch.ret.disc$fleet)
            if(!is.null(Size.com))
            {
              id.rel.cols=names(Size.com)[grep(paste(c('f','m'),collapse = '|'),names(Size.com))]
              id.rel.cols=unique(as.numeric(gsub("[a-zA-Z]", "", subset(id.rel.cols,!id.rel.cols%in%c("Nsamp")))))
              from.to=seq(id.rel.cols[which.min(abs(id.rel.cols - Length.limit))],max(id.rel.cols),by=TL.bins.cm)
              id.rel.cols=grep(paste(from.to,collapse = '|'),names(Size.com))
              id.irrel.cols=grep(paste(names(Size.com)[-c(which(c("year","Seas","Fleet","Sex","Part","Nsamp")%in%
                                                                  names(Size.com)),id.rel.cols)],collapse = '|'),
                                 names(Size.com))
              Size.com.discards=Size.com%>%filter(Fleet%in%discard.flits)
              Size.com=Size.com%>%filter(!Fleet%in%discard.flits)
              Size.com.discards_retained=Size.com.discards
              Size.com.discards_discarded=Size.com.discards
              Size.com.discards_retained[,id.rel.cols]=0
              Size.com.discards_retained$Part=2
              Size.com.discards_discarded[,id.irrel.cols]=0
              Size.com.discards_discarded$Part=1
              Size.com.discards=Size.com.discards_retained
              drop.yrs=rowSums(Size.com.discards_discarded[,id.rel.cols])
              drop.yrs=which(drop.yrs<Min.annual.obs.ktch.zone)
              if(length(drop.yrs)>0)
              {
                Size.com.discards_discarded=Size.com.discards_discarded[-drop.yrs,]
              }
              if(nrow(Size.com.discards_discarded)>0)
              {
                Size.com.discards=rbind(Size.com.discards_retained,Size.com.discards_discarded)
              }
              Size.com=rbind(Size.com,Size.com.discards)%>%
                arrange(Fleet,year)
            }
            
            if(!is.null(meanbody))
            {
              meanbody=meanbody%>%
                mutate(part=ifelse(fleet%in%discard.flits,2,part))
            }

          }
          
          #Change species specific sel pars for spatial model
          if(Neim=="sandbar shark" & !is.null(Life.history$SS_offset_selectivity))   
          {
            Life.history$SS_offset_selectivity=Life.history$SS_offset_selectivity%>%
              mutate(P_1=ifelse(Fleet=='Survey',6.5,P_1),
                     P_3=ifelse(Fleet=='Survey',0.22,P_3),
                     P_4=ifelse(Fleet=='Survey',0.09,P_4))
          }
        }
        
        #a.2 Indo IUU - F estimation
        if(Scens$Estim.Indo.IUU[s]=="Yes")
        {
          #get Indo catch
          Indo.ktch=KtCh%>%
            filter(Name==Neim & Data.set=='Indonesia')
          Indo.ktch.mean.future=mean(Indo.ktch$LIVEWT.c[(nrow(Indo.ktch)-years.futures):nrow(Indo.ktch)])
          Indo.ktch.years=Indo.ktch%>%filter(finyear%in% Indo.years.sel)   #select some years
          
          #remove Indo catch from 'other'
          id.flit.other=match('Other',FLitinFO$fleetname)
          id.yrs.indo=match(Indo.ktch$finyear,KAtch$finyear)
          
          KAtch[id.yrs.indo,match(id.flit.other,colnames(KAtch))]=KAtch[id.yrs.indo,match(id.flit.other,colnames(KAtch))]-Indo.ktch$LIVEWT.c
          add.future[,match(id.flit.other,colnames(add.future))]=add.future[,match(id.flit.other,colnames(add.future))]-Indo.ktch.mean.future
          
          #id original fleets with data 
          if(!is.null(Abund))
          {
            id.abun.flit=FLitinFO[sort(unique(Abund$index)),]%>%
              mutate(old.fleet=sort(unique(Abund$index)),new.fleet=NA)
          }
          if(!is.null(Size.com))
          {
            id.Size.com.flit=FLitinFO[sort(unique(Size.com$Fleet)),]%>%
              mutate(old.fleet=sort(unique(Size.com$Fleet)),new.fleet=NA)
          }
          if(!is.null(meanbody))
          {
            id.meanbody.flit=FLitinFO[sort(unique(meanbody$fleet)),]%>%
              mutate(old.fleet=sort(unique(meanbody$fleet)),new.fleet=NA)
            
          }
          if(!is.null(tags))
          {
            id.tag.flit=FLitinFO[sort(unique(tags$recaptures$Fleet)),]%>%
              mutate(old.fleet=sort(unique(tags$recaptures$Fleet)),new.fleet=NA)
          }
          
          #add Indo as separate fleet
          indo.fleet=FLitinFO%>%filter(fleetname=='Other')%>%mutate(fleetname='Indo.IUU')
          FLitinFO=add_row(FLitinFO, indo.fleet, .after = id.flit.other)
          
          #catch history
          Indo.dummy=KAtch%>%
            dplyr::select(finyear,as.character(id.flit.other))%>%
            rename(Indo=as.character(id.flit.other))%>%
            mutate(Indo=0)
          Indo.dummy$Indo[match(Indo.ktch$finyear,Indo.dummy$finyear)]=Indo.ktch$LIVEWT.c
          
          #Set  unknown catches to ball park      
          if(set.indo.catches.for.unknown.years)
          {
            id.unknown.yrs=Indo.dummy$finyear[which(Indo.dummy$finyear%in%indo.unknown.catch.years)]
            find.comparable.years=Indo_apprehensions%>%filter(year%in%id.unknown.yrs)
            all_values=Indo_apprehensions%>%filter(!year%in%indo.unknown.catch.years)
            for(fi in 1:nrow(find.comparable.years))
            {
              iii=find.comparable.years$year[fi]
              ref_value=Indo_apprehensions%>%filter(year==iii)%>%pull(Apprehensions)
              dumi.ktch=Indo.dummy%>%
                filter(finyear==all_values$year[which.min(abs(all_values$Apprehensions - ref_value))])%>%
                pull(Indo)
              Indo.dummy$Indo[match(iii,Indo.dummy$finyear)]=dumi.ktch
            }
          }
          
          #Set catches to 'unknown' except for some years
          if(set.indo.catches.to.unknown) Indo.dummy=Indo.dummy%>%mutate(Indo=-999)
          
          #Set catches to very low
          if(set.indo.catches.to.very.low)  Indo.dummy=Indo.dummy%>%mutate(Indo=0.001) 
          
          #keep some years catch
          if(keep.some.Indo.yrs)
          {
            Indo.dummy=Indo.dummy%>%
              mutate(Indo=replace(Indo,match(Indo.ktch.years$finyear,KAtch$finyear),Indo.ktch.years$LIVEWT.c))
          }
          
          KAtch=KAtch%>%mutate(Indo=Indo.dummy$Indo,.after = as.character(id.flit.other))
          id.Kls=match('Indo',names(KAtch)):ncol(KAtch)
          names(KAtch)[id.Kls]=seq(from = match('Indo.IUU',FLitinFO$fleetname), length.out = length(id.Kls))
          
          
          #future catch
          Indo.dummy=add.future%>%
            dplyr::select(as.character(id.flit.other))%>%
            rename(Indo=as.character(id.flit.other))%>%
            mutate(Indo=Indo.ktch.mean.future)
          add.future=add.future%>%mutate(Indo=Indo.dummy$Indo,.after = as.character(id.flit.other))
          id.Kls=match('Indo',names(add.future)):ncol(add.future)
          names(add.future)[id.Kls]=seq(from = match('Indo.IUU',FLitinFO$fleetname), length.out = length(id.Kls))
          
          #Abund  
          if(!is.null(Abund))
          {
            id.abun.flit$new.fleet=FLitinFO%>%
              mutate(row_id = row_number())%>%
              filter(fleetname%in%id.abun.flit$fleetname)%>%
              pull(row_id)
            id.abun.flit=id.abun.flit%>%dplyr::select(old.fleet,new.fleet)
            Abund=Abund%>%
              rownames_to_column(var = "row_names")%>%
              left_join(id.abun.flit,by=c('index'='old.fleet'))%>%
              mutate(index=new.fleet)%>%
              dplyr::select(-new.fleet)%>%
              column_to_rownames(var = "row_names")
          }
          #Size.com
          if(!is.null(Size.com))
          {
            id.Size.com.flit$new.fleet=FLitinFO%>%
              mutate(row_id = row_number())%>%
              filter(fleetname%in%id.Size.com.flit$fleetname)%>%
              pull(row_id)
            id.Size.com.flit=id.Size.com.flit%>%dplyr::select(old.fleet,new.fleet)
            Size.com=Size.com%>%
              left_join(id.Size.com.flit,by=c('Fleet'='old.fleet'))%>%
              mutate(Fleet=new.fleet)%>%
              dplyr::select(-new.fleet)
            
          }
          #meanbody
          if(!is.null(meanbody))
          {
            id.meanbody.flit$new.fleet=FLitinFO%>%
              mutate(row_id = row_number())%>%
              filter(fleetname%in%id.meanbody.flit$fleetname)%>%
              pull(row_id)
            id.meanbody.flit=id.meanbody.flit%>%dplyr::select(old.fleet,new.fleet)
            meanbody=meanbody%>%
              left_join(id.meanbody.flit,by=c('fleet'='old.fleet'))%>%
              mutate(fleet=new.fleet)%>%
              dplyr::select(-new.fleet)
          }
          #Tags
          if(!is.null(tags))
          {
            id.tag.flit$new.fleet=FLitinFO%>%
              mutate(row_id = row_number())%>%
              filter(fleetname%in%id.tag.flit$fleetname)%>%
              pull(row_id)
            id.tag.flit=id.tag.flit%>%
              dplyr::select(old.fleet,new.fleet)
            
            tags$Initial.reporting.rate=tags$Initial.reporting.rate%>%
              left_join(id.tag.flit,by=c('Fleet'='old.fleet'))%>%
              dplyr::select(-Fleet)%>%
              rename(Fleet=new.fleet)
            tags$recaptures=tags$recaptures%>%
              left_join(id.tag.flit,by=c('Fleet'='old.fleet'))%>%
              mutate(Fleet=new.fleet)%>%
              dplyr::select(-new.fleet)
          }
          #Variance adjustment
          if(any(!is.null(Abund),!is.null(Size.com),!is.null(meanbody)))
          {
            New.var.adj=rbind(id.abun.flit,id.Size.com.flit,id.meanbody.flit)%>%
              distinct(old.fleet,.keep_all = T)
            Var.ad=Var.ad%>%
              left_join(New.var.adj,by=c('Fleet'='old.fleet'))%>%
              mutate(Fleet=new.fleet)%>%
              dplyr::select(-new.fleet)
          }
          
          #Add Apprehensions as index of abundance
          if(!is.null(Abund))
          {
            Indo.abund=Indo_apprehensions%>%
              rename(Year=year,
                     Mean=Apprehensions)%>%
              filter(Year<=max(Abund$Year))%>%
              filter(Year%in%Indo.years.cpue)%>%
              mutate(seas=1,
                     CV=CV_apprehensions,
                     index=match('Indo.IUU',FLitinFO$fleetname))%>%
              dplyr::select(names(Abund))
            if(scale.Indo.appre) Indo.abund$Mean=Indo.abund$Mean/mean(Indo.abund$Mean,na.rm=T)
            rownames(Indo.abund)=paste('Indo',1:nrow(Indo.abund),sep='.')
            Abund=rbind(Abund,Indo.abund)%>%
              arrange(index,Year)
            
          }
          
          #Add Q
          Life.history$Q.inits=Life.history$Q.inits%>%
            filter(Fleet%in%FLitinFO$fleetname)%>%
            mutate(Fleet.n.new=match(Fleet,FLitinFO$fleetname))
          Indo.Q=Life.history$Q.inits%>%
            filter(Fleet=='Other')%>%
            mutate(Fleet='Indo.IUU',
                   Fleet.n.new=match('Indo.IUU',FLitinFO$fleetname))
          Life.history$Q.inits=rbind(Life.history$Q.inits,Indo.Q)%>%
            mutate(Fleet.n=Fleet.n.new)%>%
            arrange(Fleet.n)%>%
            dplyr::select(-Fleet.n.new)%>%
            distinct(Fleet,.keep_all = TRUE)
        }
        
        #a.3 Indo IUU- Test effect of using only Apprehensions for catch reconstruction
        if(Scens$Test.Indo.IUU.catch[s]=='Yes')  
        {
          Indo.IUU=Indo.IUU.apprehensions%>%
            filter(SPECIES==List.sp[[i]]$Species)
          Indo.IUU.app=Indo.IUU.apprehensions.only%>%
            filter(SPECIES==List.sp[[i]]$Species)
          des.flt=match(match('Other',FLitinFO$fleetname),colnames(KAtch))
          des.yrs=match(as.numeric(substr(Indo.IUU$FINYEAR,1,4)),KAtch$finyear)
          new.ktch=unlist(KAtch[des.yrs,des.flt])-(Indo.IUU$LIVEWT.c/1000)+(Indo.IUU.app$LIVEWT.c/1000)
          
          KAtch[des.yrs,c(3,des.flt)]%>%
            data.frame()%>%
            mutate(Apprehensions.only=new.ktch)%>%
            rename(max.Apprehensions.Forfeitures=X2)%>%
            gather(Method,Tons,-finyear)%>%
            ggplot(aes(finyear,Tons,color=Method))+
            geom_point()+geom_line()+ylim(0,NA)+
            theme_PA()+
            theme(legend.position = 'top',
                  legend.title=element_blank())
          ggsave(paste0(this.wd,"/Catch Other fleet_Indo catch recons_Apprehensions or Forfeitures.tiff"),width=7,height=6,compression = "lzw")
          
          KAtch[des.yrs,des.flt]=new.ktch
        }
        
        #a.4 set MainRdevYrFirst   
        #note: align with data-rich years (cpue, comps, meanbody, etc)
        Abund1=Abund
        if(!is.null(Abund1)) Abund1=Abund1%>%rename_with(tolower)
        Max.yr.obs=max(unlist(lapply(list(Abund1,Size.com,meanbody),function(x) if(!is.null(x))max(x$year))))
        Life.history$MainRdevYrLast=min(Max.yr.obs,max(KAtch$finyear,na.rm=T))
        
        Min.yr.obs=min(unlist(lapply(list(Abund1,Size.com,meanbody),function(x) if(!is.null(x))min(x$year))))
        if(Life.history$First.yr.main.rec.dev=='min.obs')
        {
          Min.yR=Min.yr.obs
          if(Life.history$First.yr.main.rec.dev_buffer)  Min.yR=Min.yr.obs-round(min(Life.history$Age.50.mat))
          MainRdevYrFirst=Min.yR
        }
         if(Life.history$First.yr.main.rec.dev=='min.ktch') MainRdevYrFirst=min(ktch$finyear)
        Life.history$MainRdevYrFirst=MainRdevYrFirst
        
        #a.5 need to reset rec pars for tuning
        #if(Scens$Scenario[s]=='S1' & Tune.SS.model)  #set to NULL in exploratory phase to fully see effect of data
        {
          Life.history$recdev_early_start=0
          Life.history$SR_sigmaR=0.2
          Life.history$RecDev_Phase=3
          
          #Ramp:
          # The model linearly interpolates the adjustment fraction between these four year-markers
          
          #The last year of the early recruitment period where no bias adjustment (0%) is applied. 
          # Typically used for very early years with no data
          Life.history$last_early_yr_nobias_adj_in_MPD=MainRdevYrFirst-1
          
          #The year when the model transitions to full bias adjustment (100%). 
          # This should align with the start of informative composition data
          Life.history$first_yr_fullbias_adj_in_MPD=MainRdevYrFirst 
          
          #The last year where full bias adjustment is applied
          # This usually marks the point where recent data (like small fish in length comps) still strongly informs recruitment
          Life.history$last_yr_fullbias_adj_in_MPD=Max.yr.obs-2
          
          Life.history$first_recent_yr_nobias_adj_in_MPD=Max.yr.obs-1
          
          #The maximum fraction of the bias adjustment to apply (typically set to 1.0 for full correction). 
          # If set lower, the model never applies the full theoretical correction
          Life.history$max_bias_adj_in_MPD=1
          
          #Comps variance adjustment
          Var.ad=NULL
        }
        
        #a.6 remove 2001 and 2002 from survey as too many juveniles, shots done all over the place 
        if(!is.null(Size.com))
        {
          id.survey=match('Survey',FLitinFO$fleetname)
          if(!is.na(id.survey)) Size.com=Size.com%>%filter(!(Fleet==id.survey & year%in%c(2001,2002)))
        }

        #a.7 plot Length comps used and all 
        do.this=FALSE
        if(do.this)
        {
          #Used on model  
          p=see.SS3.length.comp.matrix(dd=Size.com%>%
                                         dplyr::select(year,Fleet,Nsamp,names(Size.com)[grep('f',names(Size.com))]),
                                       LVLS=sort(unique(Size.com_all$Fleet)),
                                       yr.LVLS=sort(unique(Size.com_all$year)))
          base::print(p)
          ggsave(paste0(this.wd1,"/Length comps used in SS_female.tiff"),width=7,height=6,compression = "lzw")
          p=see.SS3.length.comp.matrix(dd=Size.com%>%
                                         dplyr::select(year,Fleet,Nsamp,names(Size.com)[grep('m',names(Size.com))][-1]),
                                       LVLS=sort(unique(Size.com_all$Fleet)),
                                       yr.LVLS=sort(unique(Size.com_all$year)))
          base::print(p)
          ggsave(paste0(this.wd1,"/Length comps used in SS_male.tiff"),width=7,height=6,compression = "lzw")
          
          #All available
          p=see.SS3.length.comp.matrix(dd=Size.com_all%>%
                                         dplyr::select(year,Fleet,Nsamp,names(Size.com)[grep('f',names(Size.com))]),
                                       LVLS=sort(unique(Size.com_all$Fleet)),
                                       yr.LVLS=sort(unique(Size.com_all$year)))
          base::print(p)
          ggsave(paste0(this.wd1,"/Length comps used in SS_female_all years.tiff"),width=7,height=6,compression = "lzw")
          p=see.SS3.length.comp.matrix(dd=Size.com_all%>%
                                         dplyr::select(year,Fleet,Nsamp,names(Size.com)[grep('m',names(Size.com))][-1]),
                                       LVLS=sort(unique(Size.com_all$Fleet)),
                                       yr.LVLS=sort(unique(Size.com_all$year)))
          base::print(p)
          ggsave(paste0(this.wd1,"/Length comps used in SS_male_all years.tiff"),width=7,height=6,compression = "lzw")
          
        }
            
        
        
        #a.8 remove Tag data is tested in scenario
        if(Scens$Tagging[s]=='No')
        {
          tags=NULL
        }
        #a.9 remove selectivity offsets 
        Life.history1=Life.history
        if(Scens$Use.male.sel.offset[s]=='No' & any(grepl('offset_selectivity',names(Life.history1))))
        {
          Life.history1=Life.history1[-grep('offset_selectivity',names(Life.history1))]
        }
        
        #a.10 remove growth estimation    
        if(isTRUE(Life.history1$SS3.estim.growth.pars))
        {
          if(Scens$estim.growth[s]=='No')  Life.history1$SS3.estim.growth.pars=FALSE
        }
        
        #a.11 create file  
        fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SS'),   
                     new.path=this.wd1,
                     Scenario=Scens[s,]%>%mutate(Model='SS'),
                     Catch=KAtch,
                     Catch.ret.disc=KAtch.ret.disc,
                     life.history=Life.history1,
                     depletion.yr=NULL,
                     fleets=names(KAtch)[which(!names(KAtch)%in%c("SPECIES","Name","finyear"))],
                     fleetinfo=FLitinFO,
                     abundance=Abund,   
                     size.comp=Size.com,
                     meanbodywt=meanbody,
                     Tags=tags,
                     F.tagging=F.SS.format,
                     cond.age.len=Cond.age.len.SS.format,
                     MeanSize.at.Age.obs=MeanSize.at.Age.obs.SS.format,
                     Lamdas=Lamdas.SS.lambdas,
                     Var.adjust.factor=Var.ad,    
                     Future.project=add.future)    
      }      
      print(paste("Creating input pars for---------",Neim,"----------scenario",s))
    } # end s scenario
  
  }
  
}# end i species loop


#-----------  Run all scenarios and species-------------------------------------------------------------------------
Arg='-nohess'
Where.exe=handl_OneDrive('SS3/ss_win_exe_v3.30.24.1/ss3.exe')
tic("timer")
for(i in 1:N.sp)
{
  Neim=Keep.species[i]
  this.wd=paste(HandL.out,capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
  Scens=List.sp[[i]]$Sens.test$SS
  
  for(s in 1:nrow(Scens))
  {
    this.wd1=paste(this.wd,Scens$Scenario[s],sep='/')
    fn.run.SS(where.inputs=this.wd1,  where.exe=Where.exe, args=Arg)
    COVAR=FORECAST=FALSE
    if(Arg=="") COVAR=TRUE
    if("SS"%in%future.models) FORECAST=TRUE
    Report=SS_output(this.wd1,covar=COVAR,forecast=FORECAST,readwt=F)
    this.plot=1:26 # this.plot=this.plot[-21]
    SS_plots(Report,plot=this.plot,  png=T)
  }

}
toc(log = TRUE, quiet = TRUE)
computation.time <- tic.log(format = TRUE)
tic.clearlog()
computation.time

#-----------  Run specific models  (copy input files from 2026 folder into the New folder)-------------------------------------------------------------------------
i=4
Neim=Keep.species[i]

#this.wd='C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/test scenarios/tunning_Whiskery'
this.wd='C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/New folder'

Scens=list.files(this.wd)
for(s in 1:length(Scens))
{
  this.wd1=paste(this.wd,Scens[s],sep='/')
  fn.run.SS(where.inputs=this.wd1,  where.exe=Where.exe, args=Arg)
  COVAR=FORECAST=FALSE
  if(Arg=="") COVAR=TRUE
  if("SS"%in%future.models) FORECAST=TRUE
  Report=SS_output(this.wd1,covar=COVAR,forecast=FORECAST,readwt=F)
  this.plot=1:26 
  #this.plot=this.plot[-21]
  SS_plots(Report,plot=this.plot,  png=T)
}

#-----------  tune model and calculate RAMP years-------------------------------------------------------------------------

this.wd='C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/test scenarios/tunning_Whiskery'
#this.wd1=this.wd2
#tune ramp years (blue and red lines should match)
Scens=list.files(this.wd)
tic("timer")
for(s in 1:length(Scens))
{
  this.wd1=paste(this.wd,Scens[s],sep='/')
  print(paste('Tunning ------------------',Scens[s],'-------------------------'))
  
  fn.run.SS(where.inputs=this.wd1,
            where.exe=Where.exe,
            args='')
  Report=SS_output(this.wd1)
  these.plots=1:26 
  SS_plots(Report, plot=these.plots, png=T,printfolder = "1_plots_not tuned")
  tiff(file=paste(this.wd1,'Ramp_years.tiff',sep='/'),
       width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
  ramp_years=SS_fitbiasramp(Report) 
  dev.off()
  out=ramp_years$df
  out=rbind(out,data.frame(value=unique(Report$sigma_R_info$alternative_sigma_R),label='Alternative_sigma_R'))
  write.csv(out,paste(this.wd1,'Ramp_years_first round.csv',sep='/'),row.names = F)
  
  Likelihoods.not.tuned=Report$likelihoods_used%>%mutate(type='not tuned')
  
  #tune composition data
  tune_info <- tune_comps(option = "Francis",
                          niters_tuning = 1,
                          dir = this.wd1,
                          exe=Where.exe,
                          allow_up_tuning = TRUE,
                          verbose = FALSE)
  Tuned.var.adjust=tune_info$weights[[1]]%>%mutate(Method='Francis')
  write.csv(Tuned.var.adjust,paste(this.wd1,'Tuned_size_comp.csv',sep='/'),row.names = F)
  
  #3rd. Re tune ramp years (blue and red lines should match) with updtated tuned sample sizes
  #replace var adj factor with tuned values
  start <- r4ss::SS_readstarter(file = file.path(this.wd1, "starter.ss"), verbose = FALSE)
  dat <- r4ss::SS_readdat(file = file.path(this.wd1, start$datfile), verbose = FALSE)
  ctl <- r4ss::SS_readctl(file = file.path(this.wd1, start$ctlfile), verbose = FALSE, use_datlist = TRUE, datlist = dat)
  ctl$Variance_adjustment_list=with(Tuned.var.adjust,data.frame(factor=factor,fleet=fleet,value=value)) 
  
  #replace ramp years
  ctl$last_early_yr_nobias_adj= out%>%filter(grepl('last_early_yr_nobias_adj',label))%>%pull(value)%>%as.numeric()
  ctl$first_yr_fullbias_adj= out%>%filter(grepl('first_yr_fullbias_adj',label))%>%pull(value)%>%as.numeric() 
  ctl$last_yr_fullbias_adj= out%>%filter(grepl('last_yr_fullbias_adj',label))%>%pull(value)%>%as.numeric() 
  ctl$first_recent_yr_nobias_adj= out%>%filter(grepl('first_recent_yr_nobias_adj',label))%>%pull(value) %>%as.numeric()
  ctl$max_bias_adj= out%>%filter(grepl('max_bias_adj',label))%>%pull(value)%>%as.numeric()
  ctl$recdev_early_start=0
  
  r4ss::SS_writectl(ctl, outfile = file.path(this.wd1, start$ctlfile), overwrite = TRUE, verbose = FALSE)
  
  #re run ramp
  fn.run.SS(where.inputs=this.wd1,
            where.exe=Where.exe,
            args='')
  Report=SS_output(this.wd1)
  tiff(file=paste(this.wd1,'Ramp_years_second round.tiff',sep='/'),
       width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
  ramp_years=SS_fitbiasramp(Report) 
  dev.off()
  out=ramp_years$df
  out=rbind(out,data.frame(value=unique(Report$sigma_R_info$alternative_sigma_R),label='Alternative_sigma_R'))
  write.csv(out,paste(this.wd1,'Ramp_years_second round.csv',sep='/'),row.names = F)
  SS_plots(Report, plot=these.plots, png=T,printfolder = "2_plots_tuned")
  
  Likelihoods.tuned=Report$likelihoods_used%>%mutate(type='tuned')
  
  #compare tuned and not tuned likelihoods  
  rbind(Likelihoods.not.tuned%>%
          rownames_to_column(var = "Component"),
        Likelihoods.tuned%>%
          rownames_to_column(var = "Component"))%>%
    ggplot(aes(type,values))+
    geom_bar(stat = "identity")+
    facet_wrap(~Component,scales='free')+
    theme_PA()+xlab('')
  ggsave(paste(this.wd1,'Likelihoods before and after tuning.tiff',sep='/'), width = 8,height = 6,compression = "lzw")
  
  #flush
  rm(ramp_years,out,tune_info)
  
}
toc(log = TRUE, quiet = TRUE)
computation.time <- tic.log(format = TRUE)
tic.clearlog()
computation.time


#-----------  Get Reports from different models-----------------------------------------------------
CHECK.these.mods=list(S1="C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/Population dynamics/1.Dusky shark/2026/SS3 integrated/S1",
                      S1_original="C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/Population dynamics/1.Dusky shark/2026/SS3 integrated/S1_original",
                      S1_no.NDS.length="C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/Population dynamics/1.Dusky shark/2026/SS3 integrated/S1_no NDS length",
                      Tuned_desktop="C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/test scenarios/tunning_Dusky",
                      '2022'='C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/Population dynamics/1.Dusky shark/2022/SS3 integrated/S1')


for( yy in 1:length(CHECK.these.mods))
{
  CHECK.these.mods[[yy]]=SS_output(CHECK.these.mods[[yy]],covar=COVAR,forecast=FORECAST,readwt=F)
}


#-----------  Check Estimated biomass (report)-----------------------------------------------------
Kmpr.bio=vector('list',length(CHECK.these.mods))
names(Kmpr.bio)=names(CHECK.these.mods)
for( yy in 1:length(Kmpr.bio))
{
  Rep=CHECK.these.mods[[yy]]
  Kmpr.bio[[yy]]=Rep$timeseries%>%
    dplyr::select(Yr,SpawnBio)%>%
    mutate(Rel.SpawnBio=SpawnBio/SpawnBio[1],
           Mod=names(CHECK.these.mods)[yy])
  rm(Rep)
  
}
do.call(rbind,Kmpr.bio)%>%
  ggplot(aes(Yr,Rel.SpawnBio,color=Mod))+geom_point()+
  geom_line()+theme_PA()+
  theme(legend.position = 'top')+ylim(0,1)

#-----------  Check Estimated selectivities (report)-----------------------------------------------------
Kmpr.sels=vector('list',length(CHECK.these.mods))
names(Kmpr.sels)=names(CHECK.these.mods)
for( yy in 1:length(Kmpr.sels))
{
  Rep=CHECK.these.mods[[yy]]
  FLITID=data.frame(Fleet=as.character(Rep$fleet_ID),Name=Rep$FleetNames)
  Kmpr.sels[[yy]]=Rep$sizeselex%>%
                            filter(Factor=='Lsel' & Yr==max(Rep$sizeselex$Yr))%>%
                            dplyr::select(-c(Factor,Label))%>%
                            gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
                            mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))%>%
                            mutate(Mod=names(CHECK.these.mods)[yy])%>%
                            left_join(FLITID,by='Fleet')
  rm(Rep)
  
}
do.call(rbind,Kmpr.sels)%>%
  ggplot(aes(TL,Sel,color=Mod,linetype=Sex))+geom_point(aes(shape=Sex))+
  geom_line()+theme_PA()+
  facet_wrap(~Name,ncol=3)+theme(legend.position = 'top')


#-----------  See SS selectivities inputs-----------------------------------------------------
if(!exists('doubleNorm24.fn')) fn.source1("SS_selectivity functions.R")
x=seq(30,200,5)
Mod1=data.frame(p1=113.6,p2=-1.0896,p3=4.907,p4=5.6988,Name='Southern 1 West')  
Mod2=data.frame(p1=117,p2=-12.82,p3=4.79,p4=6.31,Name='2022')   
Mod3=data.frame(p1=110.2,p2=-1.11,p3=4.65,p4=6.93,Name='Southern 1 Zone1')   
Sel.mod1=with(Mod1,doubleNorm24.fn(x,a=p1,b=p2, c=p3, d=p4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))
Sel.mod2=with(Mod2,doubleNorm24.fn(x,a=p1,b=p2, c=p3, d=p4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))
Sel.mod3=with(Mod3,doubleNorm24.fn(x,a=p1,b=p2, c=p3, d=p4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))

plot(x,Sel.mod1,type='l',lwd=2)
lines(x,Sel.mod2,col=2,lwd=2,lty=2)
lines(x,Sel.mod3,col=3,lwd=2,lty=3)
legend('topright',unique(c(Mod1$Name,Mod2$Name,Mod3$Name)),col=1:3,lty=c(1,2,3),bty='n')


Sel.log=logistic1.fn(len=x,a=265, b=30)
plot(x,Sel.log)
lines(x,logistic1.fn(len=x,a=265, b=100))


#-----------  Check retentions --------------------------------------------------------

#Retention
i=1
dd=SS_selectivity_init_pars%>%filter(Species==Keep.species[i] & !is.na(Ret_p1))
L.vec=round(with(List.sp[[i]],Lzero*a_FL.to.TL+b_FL.to.TL)*.9):
  round(with(List.sp[[i]],Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL))
A=fn.SS3.retention(p1=87,   #dd$Ret_p1
                   p2=0.1,  #dd$Ret_p2
                   p3=999,  #dd$Ret_p3
                   p4=0,   #dd$Ret_p4
                   p5=193,  #dd$Ret_p5
                   p6=1,   #dd$Ret_p6
                   p7=0,   #dd$Ret_p7
                   len.vec=L.vec)
if(is.na(dd$Ret_p5))
{
  dat1=data.frame(TL=A$len.vec,Retention=A$logistic)
}
if(!is.na(dd$Ret_p5))
{
  dat1=data.frame(TL=A$len.vec,Retention=A$dome.shaped)
}
p1=dat1%>%
  ggplot(aes(TL,Retention))+
  geom_line(color=2)+ylim(0,1)+theme_PA()+xlab('')

#Discard Mortality
A=fn.SS3.discard.mort(p1=dd$Disc_Fleet_L1, 
                      p2=dd$Disc_Fleet_L2,
                      p3=dd$Disc_Fleet_L3,
                      p4=dd$Disc_Fleet_L4,
                      len.vec=L.vec)
p2=data.frame(TL=A$len.vec,Discard.mort=A$discard.mort)%>%
  ggplot(aes(TL,Discard.mort))+
  geom_line(color=2)+ylim(0,1)+theme_PA()+
  xlab('Total length (cm)')+ylab('Discard mortality')

ggarrange(p1,p2,nrow=2)  


#-----------  Check likelihood values -------------------------------------------------
this.wd1="C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/SS3/Species/Gummy 2009 tagging"
Report_gummy=SS_output(this.wd1,covar=COVAR,forecast=FORECAST,readwt=F)
Report_gummy$likelihoods_used%>%arrange(-values)

this.wd1="C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/test tagging/1_S1 drop 94_95"
Report_1_S1=SS_output(this.wd1,covar=COVAR,forecast=FORECAST,readwt=F)
Report_1_S1$likelihoods_used%>%arrange(-values)



this.wd1="C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/test tagging/1_S1 drop 94_95"
Report_tweek=SS_output(this.wd1,covar=COVAR,forecast=FORECAST,readwt=F)
Report_tweek$likelihoods_used%>%arrange(-values)

Report_1_S1$tagreportrates
Report_tweek$tagreportrates

fn.logit(-30)
fn.inv.logit(0.5)

#-----------  Check reporting rate -------------------------------------------------
time.vec=0:5
plot(time.vec,sapply(time.vec,function(x)fn.SS3.tag.reporting.rate(init.rep.rate=1,exp.decay.rate=0.1335,time=x)),
     ylim=c(0,1))

fn.logit(-0.031)
fn.logit(0.571);fn.logit(-0.248);fn.logit(-0.603);fn.logit(-20)
fn.inv.logit(0.5)





#-----------  Create and run tweaked models -------------------------------------------
#Change this values
Min.annual.obs.ktch=150  #50
Min.annual.obs.ktch.zone=100  #100
Min.Nsamp=5    #10
Min.Nsamp.zone=5   #10
drop.dodgy.len.comp=NULL
SS.sex.length.type=3
#1. Create SS files
i=4
Neim=Keep.species[i]

#this.wd='C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/test scenarios/Gummy explore/S1'
this.wd='C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/New folder'
if(!dir.exists(this.wd))dir.create(this.wd)

if(Neim=='gummy shark')
{
  Life.history=List.sp[[i]]
  
  if(Neim%in%drop.min.pop.bin.size) Life.history$Lzero=Life.history$Lzero/1.046 
  
  #1. Catch
  #1.1. zones together
  ktch=KtCh%>%
    filter(Name==Neim)
  retained.discarded.ktch=NULL
  discard.specs=c('TEP','Discards_TDGDLF')
  if(Neim%in%retained.discarded.sp)
  {
    if(Neim=="dusky shark") this.discard.sp='TEP'
    if(exists('this.discard.sp'))
    {
      retained.discarded.ktch=ktch%>%
        filter(FishCubeCode==this.discard.sp)%>%
        mutate(Fishry='Southern.shark_2')%>%
        group_by(SPECIES,Name,finyear,Fishry)%>%
        summarise(Tonnes=sum(LIVEWT.c,na.rm=T)) 
      discard.specs=subset(discard.specs,!discard.specs==this.discard.sp)
      ktch=ktch%>%
        filter(!FishCubeCode==this.discard.sp)
      if(retained.discarded.units=='numbers' & Neim=="dusky shark")  #NEW set to units 3 if retained.discarded.units=='numbers'
      {
        retained.discarded.ktch=retained.discarded.ktch%>%
          left_join(TEPS_dusky_n.discards%>%
                      mutate(finyear=as.numeric(substr(FINYEAR,1,4)))%>%
                      group_by(finyear)%>%
                      summarise(Discards.n_1000s=sum(Discards.n_1000s))%>%
                      ungroup(),
                    by='finyear')%>%
          mutate(Discards.n_1000s=ifelse(is.na(Discards.n_1000s),0,Discards.n_1000s),
                 Tonnes=0,
                 Tonnes=Discards.n_1000s)%>%
          dplyr::select(-Discards.n_1000s)
        
      }
    }
  }
  ktch=ktch%>%
    mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                         ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                  discard.specs),'Southern.shark',
                                ifelse(FishCubeCode%in%c('WRL') & Neim%in%WRL.species,'WRL',
                                       'Other'))))%>%
    group_by(SPECIES,Name,finyear,Fishry)%>%
    summarise(Tonnes=sum(LIVEWT.c,na.rm=T))%>%
    mutate(Fishry=case_when(Fishry=="Southern.shark" & finyear<2006 ~'Southern.shark_1',
                            Fishry=="Southern.shark" & finyear>=2006~'Southern.shark_2',
                            TRUE~Fishry))
  #1.2. by zone
  ktch.zone=KtCh.zone%>%
    ungroup()%>%
    filter(Name==Neim)
  
  #allocated historic Southern.shark to zones
  if('Historic'%in%unique(ktch.zone$FishCubeCode))
  {
    historic=ktch.zone%>%
      filter(FishCubeCode=='Historic')%>%
      dplyr::select(-zone)
    ktch.zone=ktch.zone%>%filter(!FishCubeCode=='Historic')
    First.year.k=ktch.zone%>%filter(FishCubeCode%in%c("JASDGDL","WCDGDL"))
    Dis.Yr=min(First.year.k$finyear)
    if(Neim=='sandbar shark') Dis.Yr=1985
    prop.k=First.year.k%>%
      filter(finyear==Dis.Yr)%>%
      group_by(zone)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c))%>%
      mutate(Prop=LIVEWT.c/sum(LIVEWT.c))%>%
      dplyr::select(-LIVEWT.c)
    dum=expand.grid(finyear=historic%>%distinct(finyear)%>%pull(finyear),
                    zone=unique(prop.k$zone))
    prop.k=full_join(prop.k,dum,by='zone')%>%arrange(finyear)
    historic=full_join(historic,prop.k,by='finyear')%>%
      mutate(LIVEWT.c=LIVEWT.c*Prop)%>%
      dplyr::select(-Prop)%>%
      relocate(names(ktch.zone))
    ktch.zone=rbind(ktch.zone,historic)
    
  }
  retained.discarded.ktch.zone=NULL  
  if(Neim%in%retained.discarded.sp)
  {
    if(Neim=="dusky shark") this.discard.sp='TEP'
    if(exists('this.discard.sp'))
    {
      retained.discarded.ktch.zone=ktch.zone%>%
        filter(FishCubeCode==this.discard.sp)%>%
        mutate(Fishry=paste('Southern.shark_2',zone,sep='_'))%>%
        group_by(SPECIES,Name,finyear,Fishry)%>%
        summarise(Tonnes=sum(LIVEWT.c,na.rm=T)) 
      ktch.zone=ktch.zone%>%
        filter(!FishCubeCode==this.discard.sp)
      
      if(retained.discarded.units=='numbers' & Neim=="dusky shark")  
      {
        retained.discarded.ktch.zone=retained.discarded.ktch.zone%>%
          left_join(TEPS_dusky_n.discards%>%
                      mutate(Fishry=paste('Southern.shark_2',zone,sep='_'),
                             finyear=as.numeric(substr(FINYEAR,1,4)))%>%
                      group_by(finyear,Fishry)%>%
                      summarise(Discards.n_1000s=sum(Discards.n_1000s))%>%
                      ungroup(),
                    by=c('finyear','Fishry'))%>%
          mutate(Discards.n_1000s=ifelse(is.na(Discards.n_1000s),0,Discards.n_1000s),
                 Tonnes=0,
                 Tonnes=Discards.n_1000s)%>%
          dplyr::select(-Discards.n_1000s)
        
      }
    }
  }
  ktch.zone=ktch.zone%>%
    mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                         ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                  discard.specs),'Southern.shark',
                                ifelse(FishCubeCode%in%c('WRL') & Neim%in%WRL.species,'WRL',
                                       'Other'))),
           Fishry=case_when(Fishry=="Southern.shark" & finyear<2006 ~'Southern.shark_1',
                            Fishry=="Southern.shark" & finyear>=2006~'Southern.shark_2',
                            TRUE~Fishry),
           Fishry=ifelse(grepl('Southern.shark',Fishry),paste(Fishry,zone,sep='_'),Fishry))%>%
    group_by(SPECIES,Name,finyear,Fishry)%>%
    summarise(Tonnes=sum(LIVEWT.c,na.rm=T))
  
  #1.3 combined catch for displaying only
  Combined.ktch=ktch.combined%>%
    filter(Name==Neim)%>%
    rename(Year=finyear,
           Total=Tonnes)%>%
    ungroup()%>%
    dplyr::select(Year,Total)%>%
    arrange(Year)%>%
    data.frame
  
  #fleets
  #zones together
  Flits.name=sort(unique(ktch$Fishry))  
  Flits=1:length(Flits.name)
  names(Flits)=Flits.name
  Flits.and.survey=data.frame(Fleet.number=c(Flits,1+length(Flits)),
                              Fleet.name=c(names(Flits),"Survey"))
  if(!"Survey"%in% names(Catch.rate.series[[i]]) & !"Size_composition_Survey"%in%names(Species.data[[i]]))
  {
    Flits.and.survey=Flits.and.survey%>%filter(!Fleet.name=="Survey")
  }
  if("Survey"%in% names(Catch.rate.series[[i]]) | "Size_composition_Survey"%in%names(Species.data[[i]]))
  {
    Flits=c(Flits,length(Flits)+1)
    names(Flits)[length(Flits)]="Survey"
  }
  ktch=ktch%>%
    spread(Fishry,Tonnes,fill=0)
  
  #by zone
  Flits.name.zone=sort(unique(ktch.zone$Fishry))  
  Flits.zone=1:length(Flits.name.zone)
  names(Flits.zone)=Flits.name.zone
  Flits.and.survey.zone=data.frame(Fleet.number=c(Flits.zone,1+length(Flits.zone)),
                                   Fleet.name=c(names(Flits.zone),"Survey"))
  if(!"Survey"%in% names(Catch.rate.series[[i]]) & !"Size_composition_Survey"%in%names(Species.data[[i]]))
  {
    Flits.and.survey.zone=Flits.and.survey.zone%>%filter(!Fleet.name=="Survey")
  }
  if("Survey"%in% names(Catch.rate.series[[i]]) | "Size_composition_Survey"%in%names(Species.data[[i]]))
  {
    Flits.zone=c(Flits.zone,length(Flits.zone)+1)
    names(Flits.zone)[length(Flits.zone)]="Survey"
  }
  ktch.zone=ktch.zone%>%
    spread(Fishry,Tonnes,fill=0)
  
  
  #2. Size composition   
  #note: commercial catch and survey. Nsamp set at number of shots
  #zones together
  Size.compo.SS.format=NULL
  Life.history$Max.population.TL=Life.history$Min.population.TL=NULL
  if(any(grepl('Size_composition',names(Species.data[[i]]))))
  {
    if(Neim=="sandbar shark")
    {
      Species.data[[i]]$Size_composition_Survey=Species.data[[i]]$Size_composition_Survey%>%
        mutate(SEX=ifelse(is.na(SEX) & FINYEAR=='2005-06' & FL==200,'F',SEX))
    }
    d.list.n.shots=Species.data[[i]][grep(paste(c("Size_composition_Survey_Observations","Size_composition_Observations",
                                                  "Size_composition_Other_Observations"),collapse="|"),
                                          names(Species.data[[i]]))]
    d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[i]]))]
    if(length(d.list)>0)
    {
      if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
      if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
      
      Max.obs.FL=max(unlist(lapply(d.list,function(x) max(x$FL,na.rm=T))))
      Max.obs.TL=with(Life.history,Max.obs.FL*a_FL.to.TL+b_FL.to.TL)
      Min.obs.FL=min(unlist(lapply(d.list,function(x) min(x$FL,na.rm=T))))
      Min.obs.TL=with(Life.history,Min.obs.FL*a_FL.to.TL+b_FL.to.TL)
      
      Min.population.TL=min(Min.obs.TL,with(Life.history,Lzero*a_FL.to.TL+b_FL.to.TL)) 
      Max.population.TL=max(Max.obs.TL,with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL))))
      Max.population.TL=TL.bins.cm*ceiling(Max.population.TL/TL.bins.cm)
      Min.population.TL=TL.bins.cm*floor(Min.population.TL/TL.bins.cm)
      
      Life.history$Max.population.TL=Max.population.TL
      Life.history$Min.population.TL=Min.population.TL
      
      for(s in 1:length(d.list))
      {
        d.list[[s]]=d.list[[s]]%>%
          mutate(fishry=ifelse(grepl("NSF.LONGLINE",names(d.list)[s]),'NSF',
                               ifelse(grepl("Survey",names(d.list)[s]),'Survey',
                                      ifelse(grepl("Other",names(d.list)[s]),'Other',
                                             'TDGDLF'))),
                 year=as.numeric(substr(FINYEAR,1,4)),
                 TL=FL*Life.history$a_FL.to.TL+Life.history$b_FL.to.TL,
                 size.class=TL.bins.cm*floor(TL/TL.bins.cm))%>%
          filter(TL>=Min.population.TL & TL<=Max.population.TL)%>%
          dplyr::rename(sex=SEX)%>%
          filter(!is.na(sex))%>%
          group_by(year,fishry,sex,size.class)%>%
          summarise(n=n())%>%
          ungroup()
        
        #add extra bins for smooth fit to size comps  
        Maximum_size=ceiling(max(Max.population.TL,with(Life.history,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)*1.06))
        Mx.size=max(d.list[[s]]$size.class)
        extra.bins=NA
        if(Maximum_size>(Mx.size))
        {
          extra.bins=seq((Mx.size+TL.bins.cm),Maximum_size,by=TL.bins.cm)
          #extra.bins=seq(Mx.size+TL.bins.cm,10*round(Maximum_size/10),by=TL.bins.cm) 
        }
        
        if(any(!is.na(extra.bins)))
        {
          add.dumi.size=d.list[[s]][1:length(extra.bins),]%>%
            mutate(size.class=extra.bins,
                   n=0)
          d.list[[s]]=rbind(d.list[[s]],add.dumi.size)
        }
      }
      d.list <- d.list[!is.na(d.list)]
      
      if(Neim%in%names(Indicator.species))
      {
        Min.size=Min.annual.obs.ktch
      }else
      {
        Min.size=Min.annual.obs.ktch*prop.min.N.accepted_other
      }
      Min.size.NSF=Min.annual.obs.ktch_NSF
      if(Neim%in%c("dusky shark")) Min.size.NSF=20
      
      # Display sex ratio by zone used in SS  
      if(First.run=="YES")
      {
        HandL=handl_OneDrive("Analyses/Population dynamics/1.")
        DiR=paste(HandL,capitalize(Neim),"/",AssessYr,"/1_Inputs/Visualise data",sep='')
        add.n.samps=Species.data[[i]]$Size_composition_Observations%>%
          filter(Method=='GN')%>%
          mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
          rename(Zone=zone)%>%
          group_by(Zone,year)%>%
          summarise(N.shots=sum(N.shots))%>%ungroup()
        fn.ktch.sex.ratio.zone_SS(size.data=d.list,Min.size=Min.size,N_sampleS=add.n.samps)
        ggsave(paste(DiR,'Sex ratio by zone_SS size comps data_single area.tiff',sep='/'), width = 5,height = 6, dpi = 300, compression = "lzw")
      }
      
      d.list=do.call(rbind,d.list)   
      if(Neim%in%combine_NSF_Survey) 
      {
        d.list=d.list%>%
          mutate(fishry=ifelse(fishry=='NSF',"Survey",fishry))
      }
      if(Neim%in%combine.sexes)
      {
        if(Neim%in%combine.sexes.survey)
        {
          d.list$sex=ifelse(d.list$fishry=="Survey",combine.sex_type,d.list$sex)
        }
        if(Neim%in%combine.sexes.nsf)
        {
          d.list$sex=ifelse(d.list$fishry=="NSF",combine.sex_type,d.list$sex)
        }
        if(Neim%in%combine.sexes.tdgdlf)
        {
          d.list=d.list%>%mutate(sex=ifelse(fishry=="TDGDLF",combine.sex_type,sex))
        }
        if(Neim%in%combine.sexes.tdgdlf.daily)
        {
          d.list=d.list%>%mutate(sex=ifelse(fishry=="TDGDLF" & year>2005,combine.sex_type,sex))
        }
        if(!Neim%in%c(combine.sexes.survey,combine.sexes.tdgdlf,combine.sexes.tdgdlf.daily))
        {
          d.list$sex=combine.sex_type 
        }
      }
      
      #combine sexes if number of obs per year >Min.size but per sex <Min.size
      Table.n=d.list%>%  
        group_by(year,fishry)%>%
        summarise(N=sum(n))%>%
        mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey,
                                        fishry=='NSF'~Min.size.NSF,
                                        TRUE~Min.size))%>%   
        filter(N>=Min.accepted.N)%>%
        mutate(dummy=paste(year,fishry)) 
      
      Table.n.sex=d.list%>%  
        group_by(year,fishry,sex)%>%
        summarise(N=sum(n))%>%
        mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey,
                                        fishry=='NSF'~Min.size.NSF,
                                        TRUE~Min.size))%>%   
        filter(N<Min.accepted.N)%>%
        mutate(dummy=paste(year,fishry))%>%
        distinct(dummy)
      
      d.list=d.list%>%
        mutate(dummy=paste(year,fishry),
               sex=ifelse(dummy%in%Table.n.sex$dummy,combine.sex_type,sex))%>%
        dplyr::select(-dummy)
      
      d.list=d.list%>%
        group_by(year,fishry,sex,size.class)%>%
        summarise(n=sum(n))%>%
        ungroup()%>%
        filter(!is.na(year))%>%
        filter(!is.na(fishry))
      
      MAXX=max(d.list$size.class)
      Tab.si.kl=table(d.list$size.class)
      if(sum(Tab.si.kl[(length(Tab.si.kl)-1):length(Tab.si.kl)])/sum(Tab.si.kl)>.05) MAXX=MAXX*1.2
      MAXX=min(MAXX,30*ceiling(with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/30))
      size.classes=seq(min(d.list$size.class),MAXX,by=TL.bins.cm)
      missing.size.classes=size.classes[which(!size.classes%in%unique(d.list$size.class))]
      if(fill.in.zeros & length(missing.size.classes)>0)
      {
        d.list=rbind(d.list,d.list[1:length(missing.size.classes),]%>%
                       mutate(size.class=missing.size.classes,
                              n=0))
      }
      
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
                                        ifelse(zone%in%c('SA','GAB.trawl'),'Other',
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
          mutate(Part=SS.part_length.comps)%>%
          dplyr::rename(Fleet=fishry,
                        Sex=sex,
                        Seas=month)%>%
          relocate(all_of(vars.head))
        #keep years with minimum number of observations
        d.list=d.list%>%
          mutate(dummy2=paste(year,Fleet))%>%
          filter(dummy2%in%unique(Table.n$dummy))%>%
          dplyr::select(-dummy2)%>%
          mutate(Sex=ifelse(Sex=='F',1,
                            ifelse(Sex=='M',2,
                                   Sex)))
        d.list=d.list%>%
          mutate(dumi.n=rowSums(d.list[,-match(c('year','Seas','Fleet','Sex','Part','Nsamp'),names(d.list))]),
                 Nsamp=ifelse(Nsamp>dumi.n,dumi.n,Nsamp))%>%
          dplyr::select(-dumi.n)
        size.flits=Flits.and.survey
        if(!"Survey"%in% names(Catch.rate.series[[i]]) & "Survey"%in%unique(d.list$Fleet))
        {
          ddummis=size.flits[1,]%>%mutate(Fleet.number=1+size.flits$Fleet.number[nrow(size.flits)],
                                          Fleet.name="Survey")
          rownames(ddummis)="Survey"
          if(!'Survey'%in%size.flits$Fleet.name) size.flits=rbind(size.flits,ddummis)
        }
        d.list=d.list%>%
          mutate(dummy.fleet=case_when(Fleet=="NSF"~'Northern.shark',
                                       Fleet=="Other"~'Other',
                                       Fleet=="TDGDLF" & year<2006~'Southern.shark_1',
                                       Fleet=="TDGDLF" & year>=2006~'Southern.shark_2',
                                       Fleet=="Survey"~'Survey'))%>%
          left_join(size.flits,by=c('dummy.fleet'='Fleet.name'))%>%
          mutate(Fleet=Fleet.number)%>%
          dplyr::select(-c(dummy.fleet,Fleet.number))%>%
          arrange(Sex,Fleet,year)
        
        d.list.0=d.list%>%filter(Sex==combine.sex_type)%>%arrange(year)
        d.list.f=d.list%>%filter(Sex==1)%>%arrange(year)  
        d.list.m=d.list%>%filter(Sex==2)%>%arrange(year)
        
        id.var.nms=which(names(d.list.f)%in%vars.head)
        id.var.nms.f=which(!names(d.list.f)%in%vars.head)
        id.var.nms.m=length(id.var.nms.f)+which(!names(d.list.m)%in%vars.head)
        
        dummy.zeros.0=d.list.0[,-match(vars.head,names(d.list.0))]
        dummy.zeros.0[,]=0
        dummy.zeros.f=d.list.f[,-match(vars.head,names(d.list.f))]
        dummy.zeros.f[,]=0
        dummy.zeros.m=d.list.m[,-match(vars.head,names(d.list.m))]
        dummy.zeros.m[,]=0
        
        if(nrow(d.list.0)>0)  
        {
          dummy.Size.compo.SS.format_Sex0=cbind(d.list.0,dummy.zeros.0)
          names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.f]=paste('f',names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.f],sep='')
          names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.m]=paste('m',names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.m],sep='')
        }
        if(nrow(d.list.f)>0)
        {
          dummy.Size.compo.SS.format_Sex=rbind(cbind(d.list.f,dummy.zeros.f),
                                               cbind(d.list.m[,match(vars.head,names(d.list.m))],
                                                     dummy.zeros.m,
                                                     d.list.m[,-match(vars.head,names(d.list.m))]))
          names(dummy.Size.compo.SS.format_Sex)[id.var.nms.f]=paste('f',names(dummy.Size.compo.SS.format_Sex)[id.var.nms.f],sep='')
          names(dummy.Size.compo.SS.format_Sex)[id.var.nms.m]=paste('m',names(dummy.Size.compo.SS.format_Sex)[id.var.nms.m],sep='')
        }  
        if(!exists('dummy.Size.compo.SS.format_Sex0') & exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format_Sex
        }
        if(exists('dummy.Size.compo.SS.format_Sex0') & !exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format_Sex0
        }
        if(exists('dummy.Size.compo.SS.format_Sex0') & exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=rbind(dummy.Size.compo.SS.format_Sex0,dummy.Size.compo.SS.format_Sex)
        }
        clear.log('dummy.Size.compo.SS.format_Sex0')
        clear.log('dummy.Size.compo.SS.format_Sex')
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          arrange(Fleet,year,Sex)
        
        #select min sample size (shots) 
        min.nsamp=Min.Nsamp
        if(!Neim%in%names(Indicator.species)) min.nsamp=ceiling(min.nsamp/2)
        size.flits.min.samp=size.flits%>%
          mutate(Min.nsamp=case_when(Fleet.name=='Northern.shark'~Min.Nsamp.NSF,
                                     Fleet.name=='Survey'~Min.Nsamp.Survey,
                                     TRUE~min.nsamp))%>%
          dplyr::select(-Fleet.name)%>%
          rename(Fleet=Fleet.number)%>%
          filter(Fleet%in%unique(dummy.Size.compo.SS.format$Fleet))
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          left_join(size.flits.min.samp,by='Fleet')
        
        dummy.Size.compo.SS.format.all=dummy.Size.compo.SS.format%>%
          dplyr::select(-Min.nsamp)
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%  
          filter(year<=max(ktch$finyear) & Nsamp>=Min.nsamp)%>%
          dplyr::select(-Min.nsamp)
        
        if(nrow(dummy.Size.compo.SS.format)>0) Size.compo.SS.format=dummy.Size.compo.SS.format
        
      }
    }
  }
  # by zones 
  Size.compo.SS.format.zone=NULL
  if(any(grepl('Size_composition',names(Species.data[[i]]))))
  {
    d.list.n.shots=Species.data[[i]][grep(paste(c("Size_composition_Survey_Observations","Size_composition_Observations",
                                                  "Size_composition_Other_Observations"),collapse="|"),
                                          names(Species.data[[i]]))]
    d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[i]]))]
    if(Neim%in%names(drop.dodgy.len.comp))
    {
      diszone=drop.dodgy.len.comp[[match(Neim,names(drop.dodgy.len.comp))]]
      dis.yrs=unique(word(diszone, 1, sep = "-"))
      diszone=unique(str_remove(diszone, ".*-"))
      for(xx in 1:length(diszone))
      {
        for(yy in 1:length(dis.yrs))
        {
          for(qq in 1:length(d.list))
          {
            if(grepl(diszone[xx],names(d.list)[qq]))
            {
              d.list[[qq]]=d.list[[qq]]%>%
                filter(!FINYEAR%in%paste(as.numeric(dis.yrs[yy]),
                                         substr(as.numeric(dis.yrs[yy])+1,3,4),sep='-'))
            }
          }
        }
      }
    }
    
    if(length(d.list)>0)
    {
      if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
      if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
      
      for(s in 1:length(d.list))
      {
        NM=names(d.list)[s]
        d.list[[s]]=d.list[[s]]%>%
          filter(FL>=Life.history$Lzero)%>%
          mutate(fishry=ifelse(grepl("NSF.LONGLINE",NM),'NSF',
                               ifelse(grepl("Survey",NM),'Survey',
                                      ifelse(grepl("Other",NM),'Other',
                                             'TDGDLF'))),
                 year=as.numeric(substr(FINYEAR,1,4)),
                 TL=FL*Life.history$a_FL.to.TL+Life.history$b_FL.to.TL,
                 size.class=TL.bins.cm*floor(TL/TL.bins.cm))%>%
          filter(TL<=with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL))))%>%
          dplyr::rename(sex=SEX)%>%
          filter(!is.na(sex))%>%
          group_by(year,fishry,sex,size.class)%>%
          summarise(n=n())%>%
          ungroup()
        if(grepl(paste(c('West','Zone1','Zone2'),collapse='|'),NM))
        {
          ZnE=sub("\\..*", "", str_remove(NM,"Size_composition_"))
          d.list[[s]]=d.list[[s]]%>%
            mutate(fishry=paste(fishry,ZnE,sep='_'))
        }
        
        
        #add extra bins for smooth fit to size comps
        Maximum_size=ceiling(max(Max.population.TL,with(Life.history,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)*1.06))
        Mx.size=max(d.list[[s]]$size.class)
        extra.bins=NA
        if(Maximum_size>(Mx.size))
        {
          extra.bins=seq((Mx.size+TL.bins.cm),Maximum_size,by=TL.bins.cm)
          #extra.bins=seq(Mx.size+TL.bins.cm,10*round(Maximum_size/10),by=TL.bins.cm)
        }
        
        if(any(!is.na(extra.bins)))
        {
          add.dumi.size=d.list[[s]][1:length(extra.bins),]%>%
            mutate(size.class=extra.bins,
                   n=0)
          d.list[[s]]=rbind(d.list[[s]],add.dumi.size)
        }
      }
      d.list <- d.list[!is.na(d.list)]
      if(Neim%in%names(Indicator.species))
      {
        Min.size=Min.annual.obs.ktch.zone
      }else
      {
        Min.size=Min.annual.obs.ktch.zone*prop.min.N.accepted_other
      }
      Min.size.NSF=Min.annual.obs.ktch_NSF
      if(Neim%in%c("dusky shark")) Min.size.NSF=20
      
      
      # Display sex ratio by zone used in SS  
      if(First.run=="YES")
      {
        HandL=handl_OneDrive("Analyses/Population dynamics/1.")
        DiR=paste(HandL,capitalize(Neim),"/",AssessYr,"/1_Inputs/Visualise data",sep='')
        add.n.samps=Species.data[[i]]$Size_composition_Observations%>%
          filter(Method=='GN')%>%
          mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
          rename(Zone=zone)%>%
          group_by(Zone,year)%>%
          summarise(N.shots=sum(N.shots))%>%ungroup()
        fn.ktch.sex.ratio.zone_SS(size.data=d.list,Min.size=Min.size,N_sampleS=add.n.samps)
        ggsave(paste(DiR,'Sex ratio by zone_SS size comps data_areas as fleets.tiff',sep='/'), width = 5,height = 6, dpi = 300, compression = "lzw")
      }
      d.list=do.call(rbind,d.list)   
      if(Neim%in%combine_NSF_Survey) 
      {
        d.list=d.list%>%
          mutate(fishry=ifelse(fishry=='NSF',"Survey",fishry))
      }
      if(Neim%in%combine.sexes)
      {
        if(Neim%in%combine.sexes.survey)
        {
          d.list$sex=ifelse(d.list$fishry=="Survey",combine.sex_type,d.list$sex)
        }
        if(Neim%in%combine.sexes.nsf)
        {
          d.list$sex=ifelse(d.list$fishry=="NSF",combine.sex_type,d.list$sex)
        }
        if(Neim%in%combine.sexes.tdgdlf)
        {
          d.list=d.list%>%mutate(sex=ifelse(grepl("TDGDLF",fishry),combine.sex_type,sex))
        }
        if(Neim%in%combine.sexes.tdgdlf.daily)
        {
          d.list=d.list%>%mutate(sex=ifelse(grepl("TDGDLF",fishry) & year>2005,combine.sex_type,sex))
        }
        if(!Neim%in%c(combine.sexes.survey,combine.sexes.tdgdlf,combine.sexes.tdgdlf.daily))
        {
          d.list$sex=combine.sex_type 
        }
      }
      #combine sexes if number of obs per year >Min.size but per sex <Min.size
      Table.n=d.list%>%
        group_by(year,fishry)%>%
        summarise(N=sum(n))%>%
        mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey,
                                        fishry=='NSF'~Min.size.NSF,
                                        TRUE~Min.size))%>%
        filter(N>=Min.accepted.N)%>%
        mutate(dummy=paste(year,fishry))
      Table.n.sex=d.list%>%  
        group_by(year,fishry,sex)%>%
        summarise(N=sum(n))%>%
        mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey,
                                        fishry=='NSF'~Min.size.NSF,
                                        TRUE~Min.size))%>%   
        filter(N<Min.accepted.N)%>%
        mutate(dummy=paste(year,fishry))%>%
        distinct(dummy)
      d.list=d.list%>%
        mutate(dummy=paste(year,fishry),
               sex=ifelse(dummy%in%Table.n.sex$dummy,combine.sex_type,sex))%>%
        dplyr::select(-dummy)
      
      d.list=d.list%>%
        group_by(year,fishry,sex,size.class)%>%
        summarise(n=sum(n))%>%
        ungroup()%>%
        filter(!is.na(year))%>%
        filter(!is.na(fishry))
      
      MAXX=max(d.list$size.class)
      Tab.si.kl=table(d.list$size.class)
      if(sum(Tab.si.kl[(length(Tab.si.kl)-1):length(Tab.si.kl)])/sum(Tab.si.kl)>.05) MAXX=MAXX*1.2
      MAXX=min(MAXX,30*ceiling(with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/30))
      size.classes=seq(min(d.list$size.class),MAXX,by=TL.bins.cm)
      missing.size.classes=size.classes[which(!size.classes%in%unique(d.list$size.class))]
      if(fill.in.zeros & length(missing.size.classes)>0)
      {
        d.list=rbind(d.list,d.list[1:length(missing.size.classes),]%>%
                       mutate(size.class=missing.size.classes,
                              n=0))
      }
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
            mutate(fishry=ifelse(Method=='GN',paste('TDGDLF',zone,sep='_'),
                                 ifelse(Method=='LL' & zone!='SA','NSF',
                                        ifelse(zone%in%c('SA','GAB.trawl'),'Other',
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
          mutate(Part=SS.part_length.comps)%>%
          dplyr::rename(Fleet=fishry,
                        Sex=sex,
                        Seas=month)%>%
          relocate(all_of(vars.head))
        #keep years with minimum number of observations
        d.list=d.list%>%
          mutate(dummy2=paste(year,Fleet))%>%
          filter(dummy2%in%unique(Table.n$dummy))%>%
          dplyr::select(-dummy2)%>%
          mutate(Sex=ifelse(Sex=='F',1,
                            ifelse(Sex=='M',2,
                                   Sex)))
        d.list=d.list%>%
          mutate(dumi.n=rowSums(d.list[,-match(c('year','Seas','Fleet','Sex','Part','Nsamp'),names(d.list))]),
                 Nsamp=ifelse(Nsamp>dumi.n,dumi.n,Nsamp))%>%
          dplyr::select(-dumi.n)
        size.flits.zone=Flits.and.survey.zone
        if(!"Survey"%in% names(Catch.rate.series[[i]]) & "Survey"%in%unique(d.list$Fleet))
        {
          ddummis=size.flits.zone[1,]%>%mutate(Fleet.number=1+size.flits.zone$Fleet.number[nrow(size.flits.zone)],
                                               Fleet.name="Survey")
          rownames(ddummis)="Survey"
          if(!'Survey'%in%size.flits.zone$Fleet.name) size.flits.zone=rbind(size.flits.zone,ddummis)
        }
        d.list=d.list%>%
          mutate(dummy.fleet=case_when(Fleet=="NSF"~'Northern.shark',
                                       Fleet=="Other"~'Other',
                                       Fleet=="TDGDLF_West" & year<2006~'Southern.shark_1_West',
                                       Fleet=="TDGDLF_West" & year>=2006~'Southern.shark_2_West',
                                       Fleet=="TDGDLF_Zone1" & year<2006~'Southern.shark_1_Zone1',
                                       Fleet=="TDGDLF_Zone1" & year>=2006~'Southern.shark_2_Zone1',
                                       Fleet=="TDGDLF_Zone2" & year<2006~'Southern.shark_1_Zone2',
                                       Fleet=="TDGDLF_Zone2" & year>=2006~'Southern.shark_2_Zone2',
                                       Fleet=="Survey"~'Survey'))%>%
          left_join(size.flits.zone,by=c('dummy.fleet'='Fleet.name'))%>%
          mutate(Fleet=Fleet.number)%>%
          dplyr::select(-c(dummy.fleet,Fleet.number))%>%
          arrange(Sex,Fleet,year)
        
        d.list.0=d.list%>%filter(Sex==combine.sex_type)%>%arrange(year)
        d.list.f=d.list%>%filter(Sex==1)%>%arrange(year)  
        d.list.m=d.list%>%filter(Sex==2)%>%arrange(year)
        
        id.var.nms=which(names(d.list.f)%in%vars.head)
        id.var.nms.f=which(!names(d.list.f)%in%vars.head)
        id.var.nms.m=length(id.var.nms.f)+which(!names(d.list.m)%in%vars.head)
        
        dummy.zeros.0=d.list.0[,-match(vars.head,names(d.list.0))]
        dummy.zeros.0[,]=0
        dummy.zeros.f=d.list.f[,-match(vars.head,names(d.list.f))]
        dummy.zeros.f[,]=0
        dummy.zeros.m=d.list.m[,-match(vars.head,names(d.list.m))]
        dummy.zeros.m[,]=0
        
        if(nrow(d.list.0)>0)  
        {
          dummy.Size.compo.SS.format_Sex0=cbind(d.list.0,dummy.zeros.0)
          names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.f]=paste('f',names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.f],sep='')
          names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.m]=paste('m',names(dummy.Size.compo.SS.format_Sex0)[id.var.nms.m],sep='')
        }
        if(nrow(d.list.f)>0)
        {
          dummy.Size.compo.SS.format_Sex=rbind(cbind(d.list.f,dummy.zeros.f),
                                               cbind(d.list.m[,match(vars.head,names(d.list.m))],
                                                     dummy.zeros.m,
                                                     d.list.m[,-match(vars.head,names(d.list.m))]))
          names(dummy.Size.compo.SS.format_Sex)[id.var.nms.f]=paste('f',names(dummy.Size.compo.SS.format_Sex)[id.var.nms.f],sep='')
          names(dummy.Size.compo.SS.format_Sex)[id.var.nms.m]=paste('m',names(dummy.Size.compo.SS.format_Sex)[id.var.nms.m],sep='')
        }  
        if(!exists('dummy.Size.compo.SS.format_Sex0') & exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format_Sex
        }
        if(exists('dummy.Size.compo.SS.format_Sex0') & !exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format_Sex0
        }
        if(exists('dummy.Size.compo.SS.format_Sex0') & exists('dummy.Size.compo.SS.format_Sex'))
        {
          dummy.Size.compo.SS.format=rbind(dummy.Size.compo.SS.format_Sex0,dummy.Size.compo.SS.format_Sex)
        }
        clear.log('dummy.Size.compo.SS.format_Sex0')
        clear.log('dummy.Size.compo.SS.format_Sex')
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          arrange(Fleet,year,Sex)
        
        #select min sample size (shots)
        min.nsamp=Min.Nsamp.zone
        if(!Neim%in%names(Indicator.species)) min.nsamp=ceiling(min.nsamp/2)
        size.flits.min.samp=size.flits.zone%>%
          mutate(Min.nsamp=case_when(Fleet.name=='Northern.shark'~Min.Nsamp.NSF,
                                     Fleet.name=='Survey'~Min.Nsamp.Survey,
                                     TRUE~min.nsamp))%>%
          dplyr::select(-Fleet.name)%>%
          rename(Fleet=Fleet.number)%>%
          filter(Fleet%in%unique(dummy.Size.compo.SS.format$Fleet))
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          left_join(size.flits.min.samp,by='Fleet')
        
        dummy.Size.compo.SS.format.all_zone=dummy.Size.compo.SS.format%>%
          dplyr::select(-Min.nsamp)
        dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
          filter(year<=max(ktch$finyear) & Nsamp>=Min.nsamp)%>%
          dplyr::select(-Min.nsamp)
        
        
        if(nrow(dummy.Size.compo.SS.format)>0) Size.compo.SS.format.zone=dummy.Size.compo.SS.format
        
      }
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
    if(!is.null(Size.compo.SS.format.zone))
    {
      Fleet.more.one.year.obs=Size.compo.SS.format.zone%>%
        distinct(Fleet,year)%>%
        group_by(Fleet)%>%
        tally()%>%
        filter(n>1)%>%
        pull(Fleet)
      Size.compo.SS.format.zone=Size.compo.SS.format.zone%>%
        filter(Fleet%in%Fleet.more.one.year.obs)
    }
  }
  
  #Reset sex type if required  
  if(SS.sex.length.type==3)
  {
    #zones combined  
    Kombo=Size.compo.SS.format%>%distinct(year,Seas,Fleet)
    Kombo.list=vector('list',nrow(Kombo))  
    for(kk in 1:length(Kombo.list))
    {
      dumi.kk=Size.compo.SS.format%>%
        filter(year==Kombo$year[kk] & Seas==Kombo$Seas[kk] & Fleet==Kombo$Fleet[kk])
      if(nrow(dumi.kk)==1 & SS.sex.3_use.missing.sex) 
      {
        if(!dumi.kk$Sex==0 & !SS.sex.1_2.dont.change.to_3) 
        {
          msin.sx=which(!1:2%in%dumi.kk$Sex)
          dumi.kk.missing=dumi.kk%>%
            mutate(Sex=msin.sx,
                   across(-all_of(c('year','Seas','Fleet','Sex','Part','Nsamp')), ~0))
          dumi.kk=rbind(dumi.kk,dumi.kk.missing)%>%
            arrange(Sex)
        }
      }
      if(!0%in%dumi.kk$Sex)
      {
        if(!SS.sex.1_2.dont.change.to_3) if(nrow(dumi.kk)<2) dumi.kk=NA   
        if(any(!is.na(dumi.kk)))
        {
          if(nrow(dumi.kk)==2)
          {
            dumi.kk.1=dumi.kk%>%filter(Sex==1)
            dumi.kk.2=dumi.kk%>%filter(Sex==2)
            dumi.kk.3=dumi.kk.1
            dumi.kk.3[,grepl('m',colnames(dumi.kk.3))]=dumi.kk.2[,grepl('m',colnames(dumi.kk.2))]
            dumi.kk.3$Sex=SS.sex.length.type
            dumi.kk=dumi.kk.3
          }
        }
      }
      Kombo.list[[kk]]=dumi.kk
    }
    Kombo.list=Kombo.list[!is.na(Kombo.list)]
    Size.compo.SS.format=do.call(rbind,Kombo.list)
    
    #by zone
    Kombo=Size.compo.SS.format.zone%>%distinct(year,Seas,Fleet)
    Kombo.list=vector('list',nrow(Kombo))  
    for(kk in 1:length(Kombo.list))
    {
      dumi.kk=Size.compo.SS.format.zone%>%
        filter(year==Kombo$year[kk] & Seas==Kombo$Seas[kk] & Fleet==Kombo$Fleet[kk])
      if(nrow(dumi.kk)==1 & SS.sex.3_use.missing.sex.zone)
      {
        if(!dumi.kk$Sex==0 & !SS.sex.1_2.dont.change.to_3)
        {
          msin.sx=which(!1:2%in%dumi.kk$Sex)
          dumi.kk.missing=dumi.kk%>%
            mutate(Sex=msin.sx,
                   across(-all_of(c('year','Seas','Fleet','Sex','Part','Nsamp')), ~0))
          dumi.kk=rbind(dumi.kk,dumi.kk.missing)%>%
            arrange(Sex) 
        }
      }
      if(!0%in%dumi.kk$Sex)
      {
        if(!SS.sex.1_2.dont.change.to_3) if(nrow(dumi.kk)<2) dumi.kk=NA
        if(any(!is.na(dumi.kk)))
        {
          if(nrow(dumi.kk)==2)
          {
            dumi.kk.1=dumi.kk%>%filter(Sex==1)
            dumi.kk.2=dumi.kk%>%filter(Sex==2)
            dumi.kk.3=dumi.kk.1
            dumi.kk.3[,grepl('m',colnames(dumi.kk.3))]=dumi.kk.2[,grepl('m',colnames(dumi.kk.2))]
            dumi.kk.3$Sex=SS.sex.length.type
            dumi.kk=dumi.kk.3
          }
        }
      }
      Kombo.list[[kk]]=dumi.kk
    }
    Kombo.list=Kombo.list[!is.na(Kombo.list)]
    Size.compo.SS.format.zone=do.call(rbind,Kombo.list)
  }
  
  
  #3. meanbodywt
  meanbody.part=SS.part_meanbodywt
  if(!Neim%in%retained.discarded.sp) meanbody.part=0   
  #zones together
  meanbodywt.SS.format=NULL
  if(any(grepl('annual.mean.size',names(Species.data[[i]]))))
  {
    meanbodywt.SS.format=Species.data[[i]]$annual.mean.size%>%
      mutate(year=as.numeric(substr(Finyear,1,4)),
             month=1,
             Fleet='Southern.shark_2',
             part=meanbody.part,   
             type=2)%>%
      filter(year<=max(ktch$finyear))%>%
      dplyr::select(-Finyear)%>%
      left_join(Flits.and.survey,by=c('Fleet'='Fleet.name'))%>%
      mutate(fleet=Fleet.number)%>%
      dplyr::select(-c(Fleet.number,Fleet))%>%
      relocate(year,month,fleet,part,type,mean,CV)
    
    #reset very low CVs
    NewCVs=Francis.function(cipiuis=meanbodywt.SS.format%>%
                              dplyr::select(year,mean,fleet)%>%
                              spread(fleet,mean),
                            cvs=meanbodywt.SS.format%>%
                              dplyr::select(year,CV,fleet)%>%
                              spread(fleet,CV),
                            mininum.mean.CV=default.Mean.weight.CV)
    if(mean(meanbodywt.SS.format$CV,na.rm=T)<default.Mean.weight.CV)
    {
      newcv=NewCVs$CV.Adjusted%>%
        filter(year%in%meanbodywt.SS.format$year)
      meanbodywt.SS.format=meanbodywt.SS.format%>%mutate(CV=newcv[,2])
    }
    #CV variance adjustment if small CVs
    if(mean(meanbodywt.SS.format$CV,na.rm=T)>=default.Mean.weight.CV) NewCVs$CV.var.adj=0
    Var.ad.factr_meanbodywt=data.frame(Factor=3,
                                       Fleet=unique(meanbodywt.SS.format$fleet),
                                       Value=NewCVs$CV.var.adj)
  }
  # by zones 
  meanbodywt.SS.format.zone=NULL
  if(any(grepl('annual.mean.size',names(Species.data[[i]]))))
  {
    xxx=Species.data[[i]][grep('annual.mean.size',names(Species.data[[i]]))]
    xxx=xxx[-match("annual.mean.size",names(xxx))]
    if(length(xxx)>0)
    {
      for(y in 1:length(xxx))
      {
        ZnE=capitalize(str_remove(names(xxx)[y],"annual.mean.size_"))
        xxx[[y]]=xxx[[y]]%>%
          mutate(year=as.numeric(substr(Finyear,1,4)),
                 month=1,
                 Fleet=paste('Southern.shark_2',ZnE,sep='_'),
                 part=meanbody.part,      
                 type=2)%>%
          filter(year<=max(ktch.zone$finyear))%>%
          dplyr::select(-Finyear)
      }
      
      
      meanbodywt.SS.format.zone=do.call(rbind,xxx)%>%
        left_join(Flits.and.survey.zone,by=c('Fleet'='Fleet.name'))%>%
        mutate(fleet=Fleet.number)%>%
        dplyr::select(-c(Fleet.number,Fleet,zone))%>%
        relocate(year,month,fleet,part,type,mean,CV)
      
      #reset very low CVs
      NewCVs=Francis.function(cipiuis=meanbodywt.SS.format.zone%>%
                                dplyr::select(year,mean,fleet)%>%
                                spread(fleet,mean),
                              cvs=meanbodywt.SS.format.zone%>%
                                dplyr::select(year,CV,fleet)%>%
                                spread(fleet,CV),
                              mininum.mean.CV=default.Mean.weight.CV)
      if(mean(meanbodywt.SS.format.zone$CV,na.rm=T)<default.Mean.weight.CV)
      {
        newcv=NewCVs$CV.Adjusted%>%
          filter(year%in%meanbodywt.SS.format.zone$year)%>%
          gather(fleet,CV,-year)%>%mutate(fleet=as.numeric(fleet))
        meanbodywt.SS.format.zone=meanbodywt.SS.format.zone%>%
          dplyr::select(-CV)%>%
          left_join(newcv,by=c('year','fleet'))%>%
          relocate(names(meanbodywt.SS.format))
      }
      #CV variance adjustment if small CVs
      if(mean(meanbodywt.SS.format.zone$CV,na.rm=T)>=default.Mean.weight.CV) NewCVs$CV.var.adj=0
      Var.ad.factr_meanbodywt.zone=data.frame(Factor=3,
                                              Fleet=unique(meanbodywt.SS.format.zone$fleet),
                                              Value=NewCVs$CV.var.adj)
      
    }
  }
  
  
  #4. Abundance series           
  Abundance.SS.format=Abundance.SS.format.zone=NULL
  Var.ad.factr=Var.ad.factr.zone=NULL
  CPUE=compact(Catch.rate.series[[i]])
  if(!is.null(CPUE))
  {
    if(Neim%in%NSF_not.representative & any(grepl("NSF",names(CPUE)))) CPUE=CPUE[-grep("NSF",names(CPUE))]
    if(Neim%in%tdgdlf_monthly_not.representative & "TDGDLF.monthly"%in%names(CPUE)) CPUE=CPUE[-grep("TDGDLF.monthly",names(CPUE))]
    CPUE.zone=CPUE
    DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
    if(length(DROP)>0)CPUE=CPUE[-DROP]
    DROP.zone=match(c('observer',"TDGDLF.daily","TDGDLF.monthly"),names(CPUE.zone))
    DROP.zone=subset(DROP.zone,!is.na(DROP.zone))
    if(length(DROP.zone)>0)CPUE.zone=CPUE.zone[-DROP.zone]
    
    #reset very low CVs
    #note: Andre suggested leaving original CVs and estimating extraSD if more than one index available
    #      If only 1 index available, then do not estimate, just increase CV before fitting model
    #      ICCAT and SEDAR leave CVs as is and don't estimate extraSD but add variance adjustments factors to the control file
    #zones together
    if(length(CPUE)>0)
    {
      dumi.cpue=CPUE
      for(j in 1:length(CPUE))
      {
        dumi.cpue[[j]]=dumi.cpue[[j]]%>%
          mutate(Fleet=names(dumi.cpue)[j])%>%
          dplyr::select(yr.f,Mean,Fleet,CV)
      }
      NewCVs=Francis.function(cipiuis=do.call(rbind,dumi.cpue)%>%
                                dplyr::select(yr.f,Mean,Fleet)%>%
                                spread(Fleet,Mean),
                              cvs=do.call(rbind,dumi.cpue)%>%
                                dplyr::select(yr.f,CV,Fleet)%>%
                                spread(Fleet,CV),
                              mininum.mean.CV=default.CV) 
      for(j in 1:length(CPUE))
      {
        #CPUE[[j]]=CPUE[[j]]%>%mutate(CV.Original=CV)
        if(mean(CPUE[[j]]$CV,na.rm=T)<default.CV)
        {
          newcv=NewCVs$CV.Adjusted[,match(c("yr.f",names(CPUE)[j]),names(NewCVs$CV.Adjusted))]%>%
            filter(yr.f%in%CPUE[[j]]$yr.f)
          CPUE[[j]]=CPUE[[j]]%>%mutate(CV=newcv[,2])
        }
      }
      #CV variance adjustment if small CVs
      for(j in 1:length(CPUE))
      {
        if(mean(CPUE[[j]]$CV,na.rm=T)>=default.CV)  NewCVs$CV.var.adj[j]=0
      }
      Var.ad.factr_cpue=data.frame(Factor=1,Fleet=names(NewCVs$CV.var.adj),Value=NewCVs$CV.var.adj)%>%
        mutate(Fleet=case_when(Fleet=="NSF"~"Northern.shark",
                               Fleet=="TDGDLF.monthly"~"Southern.shark_1",
                               Fleet=="TDGDLF.daily"~"Southern.shark_2",
                               TRUE~Fleet))%>%
        left_join(Flits.and.survey,by=c('Fleet'='Fleet.name'))%>%
        mutate(Fleet=Fleet.number)%>%
        dplyr::select(-Fleet.number)
      Var.ad.factr=Var.ad.factr_cpue           
      if(exists("Var.ad.factr_meanbodywt")) Var.ad.factr=rbind(Var.ad.factr,Var.ad.factr_meanbodywt)
      clear.log("Var.ad.factr_meanbodywt")  
      clear.log("Var.ad.factr_cpue")
      if(drop.intermediate.yrs)if(Whiskery.q.periods==2 & Neim=='whiskery shark') 
      {
        CPUE$TDGDLF.monthly=CPUE$TDGDLF.monthly%>%
          filter(!yr.f%in%Life.history$Yr_q_change_transition)
      }
    }
    # by zones 
    if(length(CPUE.zone)>0)
    {
      dumi.cpue=CPUE.zone
      for(j in 1:length(CPUE.zone))
      {
        dumi.cpue[[j]]=dumi.cpue[[j]]%>%
          mutate(Fleet=names(dumi.cpue)[j])%>%
          dplyr::select(yr.f,Mean,Fleet,CV)
      }
      NewCVs=Francis.function(cipiuis=do.call(rbind,dumi.cpue)%>%
                                dplyr::select(yr.f,Mean,Fleet)%>%
                                spread(Fleet,Mean),
                              cvs=do.call(rbind,dumi.cpue)%>%
                                dplyr::select(yr.f,CV,Fleet)%>%
                                spread(Fleet,CV),
                              mininum.mean.CV=default.CV) 
      for(j in 1:length(CPUE.zone))
      {
        #CPUE.zone[[j]]=CPUE.zone[[j]]%>%mutate(CV.Original=CV)
        if(mean(CPUE.zone[[j]]$CV,na.rm=T)<default.CV)
        {
          newcv=NewCVs$CV.Adjusted[,match(c("yr.f",names(CPUE.zone)[j]),names(NewCVs$CV.Adjusted))]%>%
            filter(yr.f%in%CPUE.zone[[j]]$yr.f)
          CPUE.zone[[j]]=CPUE.zone[[j]]%>%mutate(CV=newcv[,2])
        }
      }
      #CV variance adjustment if small CVs
      for(j in 1:length(CPUE.zone))
      {
        if(mean(CPUE.zone[[j]]$CV,na.rm=T)>=default.CV)  NewCVs$CV.var.adj[j]=0
      }
      Var.ad.factr_cpue.zone=data.frame(Factor=1,
                                        Fleet=names(NewCVs$CV.var.adj),
                                        Value=NewCVs$CV.var.adj)%>%
        mutate(Fleet=case_when(Fleet=="NSF"~"Northern.shark",
                               Fleet=="TDGDLF.monthly.West"~"Southern.shark_1_West",
                               Fleet=="TDGDLF.monthly.Zone1"~"Southern.shark_1_Zone1",
                               Fleet=="TDGDLF.monthly.Zone2"~"Southern.shark_1_Zone2",
                               Fleet=="TDGDLF.daily.West"~"Southern.shark_2_West",
                               Fleet=="TDGDLF.daily.Zone1"~"Southern.shark_2_Zone1",
                               Fleet=="TDGDLF.daily.Zone2"~"Southern.shark_2_Zone2",
                               TRUE~Fleet))%>%
        left_join(Flits.and.survey.zone,by=c('Fleet'='Fleet.name'))%>%
        mutate(Fleet=Fleet.number)%>%
        dplyr::select(-Fleet.number)
      Var.ad.factr.zone=Var.ad.factr_cpue.zone           
      if(exists("Var.ad.factr_meanbodywt.zone")) Var.ad.factr.zone=rbind(Var.ad.factr.zone,Var.ad.factr_meanbodywt.zone)
      clear.log("Var.ad.factr_meanbodywt.zone")
      clear.log("Var.ad.factr_cpue.zone")
      if(drop.intermediate.yrs)if(Whiskery.q.periods==2 & Neim=='whiskery shark') 
      {
        CPUE.zone$TDGDLF.monthly=CPUE.zone$TDGDLF.monthly%>%
          filter(!yr.f%in%Life.history$Yr_q_change_transition)
      }
    }
  }
  
  
  #5. Add size comp effective sample size bias adjustment    
  #note: a Value of 0 means no effect
  #zones together
  if(!is.null(Var.ad.factr))
  {
    tuned.siz.comp=Life.history$tuned_size_comp
    if(!is.null(tuned.siz.comp))Var.ad.factr=rbind(Var.ad.factr,tuned.siz.comp)
  }
  if(is.null(Var.ad.factr) &!is.null(Size.compo.SS.format))
  {
    tuned.siz.comp=Life.history$tuned_size_comp
    Var.ad.factr=tuned.siz.comp
  }
  #by zones  
  if(!is.null(Var.ad.factr.zone))
  {
    tuned.siz.comp=Life.history$tuned_size_comp.zone
    if(!is.null(tuned.siz.comp))Var.ad.factr.zone=rbind(Var.ad.factr.zone,tuned.siz.comp)
  }
  if(is.null(Var.ad.factr.zone) &!is.null(Size.compo.SS.format.zone))
  {
    tuned.siz.comp=Life.history$tuned_size_comp.zone
    Var.ad.factr.zone=tuned.siz.comp
  }
  
  
  #6. F from tagging studies on TDGDLF (1994-95 and 2001-03)
  F.SS.format=NULL  
  if(any(grepl('Fishing.mortality.TDGDLF',names(Species.data[[i]]))) & length(CPUE)>0)
  {
    F.SS.format=Species.data[[i]][[grep('Fishing.mortality.TDGDLF',names(Species.data[[i]]))]]%>%
      mutate(year=as.numeric(substr(Finyear,1,4)),
             month=1,
             Fleet='Southern.shark_1')%>%
      filter(year<=max(ktch$finyear))%>%
      dplyr::select(-Finyear)%>%
      left_join(Flits.and.survey,by=c('Fleet'='Fleet.name'))%>%
      mutate(fleet=Fleet.number,
             CV=0.05)%>%
      dplyr::select(-c(Fleet.number,Fleet))%>%
      relocate(year,month,fleet,Mean,CV)%>%
      arrange(year)
    
  }
  
  
  #7. Tagging data
  # zones together
  Tags.SS.format=NULL
  if(names(Species.data)[i]%in%use.tag.data)
  {
    #Extract data
    releases=Species.data[[i]]$Con_tag_SS.format_releases%>%
      filter(Yr.rel<=Last.yr.ktch.numeric)
    recaptures=Species.data[[i]]$Con_tag_SS.format_recaptures%>%
      filter(Yr.rec<=Last.yr.ktch.numeric)
    
    #Keep relevant finyears and zones  
    dis.yrs.tag=Use.these.tag.year_zones[[match(names(Species.data)[i],names(Use.these.tag.year_zones))]]
    releases=releases%>%
      mutate(dummy=paste(Yr.rel,Rel.zone))%>%
      filter(dummy%in%dis.yrs.tag)%>%
      dplyr::select(-dummy)
    recaptures=recaptures%>%filter(Tag.group%in%unique(releases$Tag.group))
    if(use.tag.rec.yrs.90percent.rec)
    {
      Table.yr.releases=table(releases$Yr.rel)
      vec=cumsum(Table.yr.releases)/sum(Table.yr.releases)
      Last.yr.rel=names(vec[which.min(abs(vec - 0.9))])
      
      Table.yr.recaptures=table(recaptures$Yr.rec)
      vec=cumsum(Table.yr.recaptures)/sum(Table.yr.recaptures)
      Last.yr.rec=names(vec[which.min(abs(vec - 0.9))])
      recaptures=recaptures%>%filter(Yr.rec<=as.numeric(Last.yr.rec))
    }

    
    #Allocate West, Zone1 and Zone 2 to Southern 1 or 2
    releases=releases%>%
      mutate(Rel.zone=case_when(Yr.rel<=2005 ~'Southern.shark_1',
                                Yr.rel>2005 ~'Southern.shark_2'))
    recaptures=recaptures%>%
      mutate(Rec.zone=case_when(Yr.rec<=2005 ~'Southern.shark_1',
                                Yr.rec>2005 ~'Southern.shark_2'))
    
    #Recalculate TagGroup
    releases=releases%>%
      arrange(Rel.zone,Yr.rel,Sex,Age)%>%
      mutate(rowid = row_number())
    TagGroup=releases%>%distinct(Tag.group,rowid)
    recaptures=recaptures%>%left_join(TagGroup,by='Tag.group')
    releases=releases%>%
      mutate(Tag.group=rowid)%>%
      dplyr::select(-rowid)
    recaptures=recaptures%>%
      mutate(Tag.group=rowid)%>%
      dplyr::select(-rowid)
    
    #group sex  
    if(taggroup.sex.combined)   
    {
      releases=releases%>%
        mutate(Sex=0)%>%
        group_by(Rel.zone,Yr.rel,season,t.fill,Sex,Age)%>%
        mutate(N.release=sum(N.release))%>%
        ungroup()%>%
        group_by(Rel.zone,Yr.rel,season,t.fill,Sex,Age) %>%
        mutate(Tag.groups = paste(Tag.group, collapse = ", "))%>%
        ungroup()
      recaptures=recaptures%>%   
        left_join(releases%>%
                    dplyr::select(Tag.group,Tag.groups),
                  by='Tag.group')%>%
        group_by(Yr.rec,season,Rec.zone,Tag.groups)%>%
        summarise(N.recapture=sum(N.recapture))%>%
        ungroup()
      releases=releases%>%
        distinct(Tag.groups,Rel.zone,Yr.rel,season,t.fill,Sex,Age,N.release)%>%
        arrange(Rel.zone,Yr.rel,season,t.fill,Sex,Age)%>%
        mutate(Tag.group = row_number())%>%
        relocate(Tag.group)
      recaptures=recaptures%>%
        left_join(releases%>%distinct(Tag.group,Tag.groups),
                  by='Tag.groups')%>%
        dplyr::select(-Tag.groups)%>%
        relocate(Tag.group)
      releases=releases%>%dplyr::select(-Tag.groups)
    }
    
    releases=releases%>%
      rename(Area=Rel.zone)%>%
      mutate(Area=1)
    get.fleet=recaptures%>%
      distinct(Yr.rec,Rec.zone)
    Rec.ZonEs=unique(recaptures$Rec.zone)
    a1=ktch%>%
      ungroup()%>%
      dplyr::select(-c(SPECIES,Name))%>%
      gather(Fleet,Ktch,-finyear)%>%
      mutate(zone=case_when(Fleet=="Northern.shark"~'North',
                            grepl("Southern.shark_1",Fleet)~"Southern.shark_1",
                            grepl("Southern.shark_2",Fleet)~"Southern.shark_2",
                            TRUE~''))%>%
      filter(Ktch>0)%>%
      filter(zone%in%unique(get.fleet$Rec.zone))%>%
      filter(finyear%in%unique(get.fleet$Yr.rec))%>%
      distinct(finyear,Fleet,zone)%>%
      left_join(data.frame(Fleet.ID=Flits,Fleet=names(Flits)),
                by='Fleet')
    get.fleet=get.fleet%>%
      left_join(a1,by=c('Rec.zone'='zone','Yr.rec'='finyear'))
    recaptures=recaptures%>%
      left_join(get.fleet%>%dplyr::select(-Fleet)%>%rename(Fleet=Fleet.ID),
                by=c('Rec.zone','Yr.rec'))%>%
      dplyr::select(-Rec.zone)%>%
      relocate(Tag.group,Yr.rec,season,Fleet,N.recapture)
    
    Initial.tag.loss=1e-8  #tag-induced mortality immediately after tagging 
    Chronic.tag.loss=Species.data[[i]]$Con_tag_shedding_from_F.estimation.R_$x  #annual rate of tag loss; McAuley et al 2007 tag shedding
    
    if(Reporting.rate.type[[Neim]]=='published') Initial.reporting.rate=Species.data[[i]]$Con_tag_non_reporting_from_F.estimation.R_   #NEW
    if(Reporting.rate.type[[Neim]]=='calculated') Initial.reporting.rate=Species.data[[i]]$Con_tag_non_reporting_from_F.estimation.R_calculated
    Initial.reporting.rate=Initial.reporting.rate%>%
      dplyr::select(-Species)%>%
      gather(Zone,Non.reporting,-Finyear)%>%
      mutate(Reporting=1-Non.reporting,
             Zone=case_when(Zone%in%c('South','South.west','West') & Finyear<=2005 ~'Southern.shark_1',
                            Zone%in%c('South','South.west','West') & Finyear>2005 ~'Southern.shark_2',
                            Zone=='North'~'North'))%>%
      filter(Zone%in%Rec.ZonEs)%>%
      group_by(Finyear,Zone)%>%
      summarise(Non.reporting=mean(Non.reporting,na.rm=T),
                Reporting=mean(Reporting,na.rm=T))%>%
      ungroup()%>%
      mutate(Reporting.logit=fn.inv.logit(Reporting))%>%
      dplyr::select(-Non.reporting)
    get.fleet1=get.fleet%>%
      dplyr::select(-Fleet)%>%
      rename(Fleet=Fleet.ID)%>%
      mutate(dummy=paste(Yr.rec,Rec.zone))
    not.in.init.rep=paste(Initial.reporting.rate$Finyear,Initial.reporting.rate$Zone)
    not.in.init.rep=not.in.init.rep[which(!not.in.init.rep%in%get.fleet1$dummy)]
    if(length(not.in.init.rep)>0)
    {
      ad.get.flit=get.fleet1[1:length(not.in.init.rep),]%>%
        mutate(dummy=not.in.init.rep,
               Yr.rec=word(dummy, 1),
               Rec.zone=word(dummy, 2),
               Fleet=NA)
      get.fleet1=rbind(get.fleet1,ad.get.flit)%>%
        arrange(Rec.zone,Yr.rec)%>%
        fill(Fleet, .direction = "down") #NEW
      
      get.fleet1=get.fleet1%>%  
        left_join(a1%>%
                    mutate(finyear=as.character(finyear))%>%
                    distinct(finyear,zone,Fleet.ID),
                  by=c('Yr.rec'='finyear','Rec.zone'='zone'))%>%
        mutate(Fleet=case_when(!is.na(Fleet.ID) & !Fleet==Fleet.ID~Fleet.ID,
                               !is.na(Fleet.ID) & is.na(Fleet)~Fleet.ID,
                               TRUE~Fleet))%>%
        dplyr::select(-Fleet.ID)
    }
    get.fleet1=get.fleet1%>%dplyr::select(-dummy)%>%mutate(Yr.rec=as.numeric(Yr.rec))
    Initial.reporting.rate=Initial.reporting.rate%>%
      left_join(get.fleet1,
                by=c('Zone'='Rec.zone','Finyear'='Yr.rec'))%>%
      filter(!is.na(Reporting))
    if(estimate.tag.report.decay)
    {
      rep.dec.flit=sort(unique(Initial.reporting.rate$Fleet))
      Rep.decay=vector('list',length(rep.dec.flit))
      Rep.decay_p=Rep.decay
      for(re in 1:length(Rep.decay))
      {
        d.init.rep=Initial.reporting.rate%>%
          filter(Fleet==rep.dec.flit[re])%>%
          arrange(Finyear)%>%
          mutate(time=Finyear-Finyear[1])
        diKay=0  #NEW
        if(nrow(d.init.rep)>1)#NEW
        {
          Init.rep=d.init.rep$Reporting[1]
          fit_nls <- nls(Reporting ~ Init.rep * exp(-k * time), 
                         data = d.init.rep, 
                         start = list(k = 0.01))
          diKay=round(coef(fit_nls),4)              #NEW
          if(!allow.increase.tag.rep.rate) diKay=max(0,diKay) #NEW
        }
        
        Rep.decay[[re]]=data.frame(Fleet=rep.dec.flit[re],decay=diKay)
        Rep.decay_p[[re]]=ggplot(data=d.init.rep,aes(Finyear,Reporting))+
          geom_point(size=3)+
          ylim(0,1)+
          geom_line(data=data.frame(Finyear=d.init.rep$Finyear,
                                    Reporting=fn.SS3.tag.reporting.rate(init.rep.rate=d.init.rep$Reporting[1],exp.decay.rate=diKay,time=d.init.rep$time)),
                    aes(Finyear,Reporting),color=2,linewidth = 2)+
          theme_PA()+
          ggtitle(paste0('Fleet ',rep.dec.flit[re],' (',unique(d.init.rep$Zone),' decay=',diKay,')'))
        
      }
      Reporting.rate.decay=do.call(rbind,Rep.decay)
      if(First.run=='YES')
      {
        ggarrange(plotlist=Rep.decay_p,ncol=1,nrow=length(Rep.decay_p))
        ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(List.sp[[i]]$Name),"/",AssessYr,"/1_Inputs/Visualise data/Tagging_report rate decay_single.area.tiff",sep=''),
               width = 6,height = 8,compression = "lzw")
      }
      if(pass.rep.rate.decay.negative)
      {
        Reporting.rate.decay=Reporting.rate.decay%>%
          mutate(decay=ifelse(decay>0,-decay,abs(decay)))
      }
    }
    if(!estimate.tag.report.decay) Reporting.rate.decay=0 #Andre's Gummy and (Spatial SS3 workshop) Lecture D '4 areas' models
    
    if(logit.transform.tag.pars)   
    {
      Initial.tag.loss=fn.inv.logit(Initial.tag.loss)
      Chronic.tag.loss=fn.inv.logit(Chronic.tag.loss)
      Initial.reporting.rate=Initial.reporting.rate%>%
        dplyr::select(-Reporting)%>%
        rename(Reporting=Reporting.logit)
    }
    
    Tags.SS.format=list(
      releases=releases%>%data.frame,
      recaptures=recaptures%>%data.frame,
      Initial.tag.loss=Initial.tag.loss,   
      Chronic.tag.loss=Chronic.tag.loss,   
      Initial.reporting.rate=Initial.reporting.rate%>%
        filter(Finyear==Initial.reporting.rate$Finyear[which.min(abs(Initial.reporting.rate$Finyear - min(releases$Yr.rel)))]),   
      Reporting.rate.decay=Reporting.rate.decay,
      overdispersion=SS_overdispersion,       # Andre's Gummy model
      mixing_latency_period=SS_mixing_latency_period, 
      max_periods=ceiling((max(recaptures$Yr.rec)-min(releases$Yr.rel))*Extend.mx.period))  # 30 Andre's Gummy model  
    
    rm(releases,recaptures,Chronic.tag.loss,Initial.reporting.rate,Reporting.rate.decay)
  }
  
  # by zones
  Tags.SS.format.zone=NULL
  if(names(Species.data)[i]%in%use.tag.data)
  {
    #Extract data
    releases=Species.data[[i]]$Con_tag_SS.format_releases%>%
      filter(Yr.rel<=Last.yr.ktch.numeric)
    recaptures=Species.data[[i]]$Con_tag_SS.format_recaptures%>%
      filter(Yr.rec<=Last.yr.ktch.numeric)
    
    #Keep relevant finyear zones   
    dis.yrs.tag=Use.these.tag.year_zones[[match(names(Species.data)[i],names(Use.these.tag.year_zones))]]
    releases=releases%>%
      mutate(dummy=paste(Yr.rel,Rel.zone))%>%
      filter(dummy%in%dis.yrs.tag)%>%
      dplyr::select(-dummy)
    recaptures=recaptures%>%filter(Tag.group%in%unique(releases$Tag.group))
    if(use.tag.rec.yrs.90percent.rec)
    {
      Table.yr.releases=table(releases$Yr.rel)
      vec=cumsum(Table.yr.releases)/sum(Table.yr.releases)
      Last.yr.rel=names(vec[which.min(abs(vec - 0.9))])
      
      Table.yr.recaptures=table(recaptures$Yr.rec)
      vec=cumsum(Table.yr.recaptures)/sum(Table.yr.recaptures)
      Last.yr.rec=names(vec[which.min(abs(vec - 0.9))])
      recaptures=recaptures%>%filter(Yr.rec<=as.numeric(Last.yr.rec))
    }

    
    #Recalculate TagGroup
    releases=releases%>%
      arrange(Rel.zone,Yr.rel,Sex,Age)%>%
      mutate(rowid = row_number())
    TagGroup=releases%>%distinct(Tag.group,rowid)
    recaptures=recaptures%>%left_join(TagGroup,by='Tag.group')
    releases=releases%>%
      mutate(Tag.group=rowid)%>%
      dplyr::select(-rowid)
    recaptures=recaptures%>%
      mutate(Tag.group=rowid)%>%
      dplyr::select(-rowid)
    
    #group sex  
    if(taggroup.sex.combined)   
    {
      releases=releases%>%
        mutate(Sex=0)%>%
        group_by(Rel.zone,Yr.rel,season,t.fill,Sex,Age)%>%
        mutate(N.release=sum(N.release))%>%
        ungroup()%>%
        group_by(Rel.zone,Yr.rel,season,t.fill,Sex,Age) %>%
        mutate(Tag.groups = paste(Tag.group, collapse = ", "))%>%
        ungroup()
      recaptures=recaptures%>%   
        left_join(releases%>%
                    dplyr::select(Tag.group,Tag.groups),
                  by='Tag.group')%>%
        group_by(Yr.rec,season,Rec.zone,Tag.groups)%>%
        summarise(N.recapture=sum(N.recapture))%>%
        ungroup()
      releases=releases%>%
        distinct(Tag.groups,Rel.zone,Yr.rel,season,t.fill,Sex,Age,N.release)%>%
        arrange(Rel.zone,Yr.rel,season,t.fill,Sex,Age)%>%
        mutate(Tag.group = row_number())%>%
        relocate(Tag.group)
      recaptures=recaptures%>%
        left_join(releases%>%distinct(Tag.group,Tag.groups),
                  by='Tag.groups')%>%
        dplyr::select(-Tag.groups)%>%
        relocate(Tag.group)
      releases=releases%>%dplyr::select(-Tag.groups)
    }
    
    releases=releases%>%
      rename(Area=Rel.zone)%>%
      mutate(Area=1)
    get.fleet=recaptures%>%
      distinct(Yr.rec,Rec.zone)
    Rec.ZonEs=unique(recaptures$Rec.zone)
    a1=ktch.zone%>%
      ungroup()%>%
      dplyr::select(-c(SPECIES,Name))%>%
      gather(Fleet,Ktch,-finyear)%>%
      mutate(zone=case_when(Fleet=="Northern.shark"~'North',
                            grepl("Zone1",Fleet)~"Zone1",
                            grepl("West",Fleet)~"West",
                            grepl("Zone2",Fleet)~"Zone2",
                            TRUE~''))%>%
      filter(Ktch>0)%>%
      filter(zone%in%unique(get.fleet$Rec.zone))%>%
      filter(finyear%in%unique(get.fleet$Yr.rec))%>%
      distinct(finyear,Fleet,zone)%>%
      left_join(data.frame(Fleet.ID=Flits.zone,Fleet=names(Flits.zone)),
                by='Fleet')
    get.fleet=get.fleet%>%
      left_join(a1,by=c('Rec.zone'='zone','Yr.rec'='finyear'))
    recaptures=recaptures%>%
      left_join(get.fleet%>%dplyr::select(-Fleet)%>%rename(Fleet=Fleet.ID),
                by=c('Rec.zone','Yr.rec'))%>%
      dplyr::select(-Rec.zone)%>%
      relocate(Tag.group,Yr.rec,season,Fleet,N.recapture)
    
    Initial.tag.loss=1e-4  #tag-induced mortality immediately after tagging 
    Chronic.tag.loss=Species.data[[i]]$Con_tag_shedding_from_F.estimation.R_$x  #annual rate of tag loss; McAuley et al 2007 tag shedding
    
    if(Reporting.rate.type[[Neim]]=='published') Initial.reporting.rate=Species.data[[i]]$Con_tag_non_reporting_from_F.estimation.R_   #NEW
    if(Reporting.rate.type[[Neim]]=='calculated') Initial.reporting.rate=Species.data[[i]]$Con_tag_non_reporting_from_F.estimation.R_calculated
    Initial.reporting.rate=Initial.reporting.rate%>%
      dplyr::select(-Species)%>%
      gather(Zone,Non.reporting,-Finyear)%>%
      mutate(Reporting=1-Non.reporting,
             Zone=case_when(Zone=='South'~'Zone2',
                            Zone=='South.west'~'Zone1',
                            Zone=='West'~'West',
                            Zone=='North'~'North'))%>%
      filter(Zone%in%Rec.ZonEs)%>%
      mutate(Reporting.logit=fn.inv.logit(Reporting))%>%
      dplyr::select(-Non.reporting)
    get.fleet1=get.fleet%>%
      dplyr::select(-Fleet)%>%
      rename(Fleet=Fleet.ID)%>%
      mutate(dummy=paste(Yr.rec,Rec.zone))
    not.in.init.rep=paste(Initial.reporting.rate$Finyear,Initial.reporting.rate$Zone)
    not.in.init.rep=not.in.init.rep[which(!not.in.init.rep%in%get.fleet1$dummy)]
    if(length(not.in.init.rep)>0)
    {
      ad.get.flit=get.fleet1[1:length(not.in.init.rep),]%>%
        mutate(dummy=not.in.init.rep,
               Yr.rec=word(dummy, 1),
               Rec.zone=word(dummy, 2),
               Fleet=NA)
      get.fleet1=rbind(get.fleet1,ad.get.flit)%>%
        arrange(Rec.zone,Yr.rec)%>%
        fill(Fleet, .direction = "down") 
      
      get.fleet1=get.fleet1%>%  
        left_join(a1%>%
                    mutate(finyear=as.character(finyear))%>%
                    distinct(finyear,zone,Fleet.ID),
                  by=c('Yr.rec'='finyear','Rec.zone'='zone'))%>%
        mutate(Fleet=case_when(!is.na(Fleet.ID) & !Fleet==Fleet.ID~Fleet.ID,
                               !is.na(Fleet.ID) & is.na(Fleet)~Fleet.ID,
                               TRUE~Fleet))%>%#NEW
        dplyr::select(-Fleet.ID)
    }
    get.fleet1=get.fleet1%>%dplyr::select(-dummy)%>%mutate(Yr.rec=as.numeric(Yr.rec))
    Initial.reporting.rate=Initial.reporting.rate%>%
      left_join(get.fleet1,
                by=c('Zone'='Rec.zone','Finyear'='Yr.rec'))%>%
      filter(!is.na(Reporting))
    
    if(estimate.tag.report.decay)
    {
      rep.dec.flit=sort(unique(Initial.reporting.rate$Fleet))
      Rep.decay=vector('list',length(rep.dec.flit))
      Rep.decay_p=Rep.decay
      for(re in 1:length(Rep.decay))
      {
        d.init.rep=Initial.reporting.rate%>%
          filter(Fleet==rep.dec.flit[re])%>%
          arrange(Finyear)%>%
          mutate(time=Finyear-Finyear[1])
        diKay=0  #NEW
        if(nrow(d.init.rep)>1)#NEW
        {
          Init.rep=d.init.rep$Reporting[1]
          fit_nls <- nls(Reporting ~ Init.rep * exp(-k * time), 
                         data = d.init.rep, 
                         start = list(k = 0.01))
          diKay=round(coef(fit_nls),4)              #NEW
          if(!allow.increase.tag.rep.rate) diKay=max(0,diKay) #NEW
        }
        
        Rep.decay[[re]]=data.frame(Fleet=rep.dec.flit[re],decay=diKay)
        Rep.decay_p[[re]]=ggplot(data=d.init.rep,aes(Finyear,Reporting))+
          geom_point(size=3)+
          ylim(0,1)+
          geom_line(data=data.frame(Finyear=d.init.rep$Finyear,
                                    Reporting=fn.SS3.tag.reporting.rate(init.rep.rate=d.init.rep$Reporting[1],exp.decay.rate=diKay,time=d.init.rep$time)),
                    aes(Finyear,Reporting),color=2,linewidth = 2)+
          theme_PA()+
          ggtitle(paste0('Fleet ',rep.dec.flit[re],' (',unique(d.init.rep$Zone),' decay=',diKay,')'))
        
      }
      Reporting.rate.decay=do.call(rbind,Rep.decay)
      if(First.run=='YES')
      {
        ggarrange(plotlist=Rep.decay_p,ncol=1,nrow=length(Rep.decay_p))
        ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(List.sp[[i]]$Name),"/",AssessYr,"/1_Inputs/Visualise data/Tagging_report rate decay_area.as.fleets.tiff",sep=''),
               width = 6,height = 8,compression = "lzw")
      }
      if(pass.rep.rate.decay.negative)
      {
        Reporting.rate.decay=Reporting.rate.decay%>%
          mutate(decay=ifelse(decay>0,-decay,abs(decay)))
      }
    }
    if(!estimate.tag.report.decay) Reporting.rate.decay=0 #Andre's Gummy and (Spatial SS3 workshop) Lecture D '4 areas' models
    if(logit.transform.tag.pars)   
    {
      Initial.tag.loss=fn.inv.logit(Initial.tag.loss)
      Chronic.tag.loss=fn.inv.logit(Chronic.tag.loss)
      Initial.reporting.rate=Initial.reporting.rate%>%
        dplyr::select(-Reporting)%>%
        rename(Reporting=Reporting.logit)
    }
    Tags.SS.format.zone=list(
      releases=releases%>%data.frame,
      recaptures=recaptures%>%data.frame,
      Initial.tag.loss=Initial.tag.loss,   
      Chronic.tag.loss=Chronic.tag.loss,   
      Initial.reporting.rate=Initial.reporting.rate%>%
        filter(Finyear==Initial.reporting.rate$Finyear[which.min(abs(Initial.reporting.rate$Finyear - min(releases$Yr.rel)))]),
      Reporting.rate.decay=Reporting.rate.decay,
      overdispersion=SS_overdispersion,       # Andre's Gummy model
      mixing_latency_period=SS_mixing_latency_period,  
      max_periods=ceiling((max(recaptures$Yr.rec)-min(releases$Yr.rel))*Extend.mx.period))  # 30 Andre's Gummy model  
    
    rm(releases,recaptures,Chronic.tag.loss,Initial.reporting.rate,Reporting.rate.decay)
  }
  
  
  #8. Conditional age at length
  Cond.age.len.SS.format=NULL
  if(do.Cond.age.len.SS.format)
  {
    if(any(grepl('conditional_age_length',names(Species.data[[i]]))))
    {
      a=Life.history$a_FL.to.TL
      b=Life.history$b_FL.to.TL
      Cond.age.len.SS.format=Species.data[[i]]$conditional_age_length%>%
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
  
  
  #9. MeanSize at Age obs
  MeanSize.at.Age.obs.SS.format=NULL
  if(any(grepl('conditional_age_length',names(Species.data[[i]]))) & names(Species.data)[i]%in%Mean.Size.at.age.species)
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
    
    a=Life.history$a_FL.to.TL
    b=Life.history$b_FL.to.TL
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
  
  
  #10. Fleet info
  #zones together
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
  #by zone
  flitinfo.zone=data.frame(fleet=Flits.zone)%>%
    mutate(type=1,
           surveytiming=-1,
           area=1,
           units=1,
           need_catch_mult=0,
           fleetname=names(Flits.zone))%>%
    dplyr::select(-fleet)%>%
    mutate(type=ifelse(fleetname=="Survey",3,type),
           surveytiming=ifelse(fleetname=="Survey",1,surveytiming))
  rownames(flitinfo.zone)=NULL
  
  
  #11. Run scenarios if available abundance index or size comps
  len.cpue=length(CPUE)
  len.cpue.zone=length(CPUE.zone)
  len.size.comp=length(Size.compo.SS.format) 
  if(len.cpue>0 | len.cpue.zone>0 | len.size.comp>0)
  {
    #Put CPUE in SS format
    if(len.cpue>0 | len.cpue.zone>0)
    {
      MAX.CV=Life.history$MAX.CV
      
      #zones together
      if(len.cpue>0)
      {
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
        Abundance.SS.format=do.call(rbind,CPUE)%>%
          relocate(Year,seas,index,Mean,CV)%>%
          arrange(index,Year)
        
      }
      
      #by zone
      if(len.cpue.zone>0)
      {
        for(x in 1:len.cpue.zone)    
        {
          nm=names(CPUE.zone)[x]
          if(nm=="NSF") nm="Northern.shark"
          if(nm=="TDGDLF.monthly.West") nm="Southern.shark_1_West"
          if(nm=="TDGDLF.monthly.Zone1") nm="Southern.shark_1_Zone1"
          if(nm=="TDGDLF.monthly.Zone2") nm="Southern.shark_1_Zone2"
          if(nm=="TDGDLF.daily.West") nm="Southern.shark_2_West"
          if(nm=="TDGDLF.daily.Zone1") nm="Southern.shark_2_Zone1"
          if(nm=="TDGDLF.daily.Zone2") nm="Southern.shark_2_Zone2"
          dd=CPUE.zone[[x]][,grep(paste(c('yr.f','Mean','MeAn','CV'),collapse="|"),names(CPUE.zone[[x]]))]%>%
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
            left_join(Flits.and.survey.zone,by=c('index.dummy'='Fleet.name'))%>%
            mutate(index=Fleet.number)%>%
            dplyr::select(-c(index.dummy,Fleet.number))
          CPUE.zone[[x]]=dd%>%filter(!is.na(Mean))
        }
        Abundance.SS.format.zone=do.call(rbind,CPUE.zone)%>%
          relocate(Year,seas,index,Mean,CV)%>%
          arrange(index,Year)
      }
      
    }
    
    #Scenarios
    Scens=Life.history$Sens.test$SS
    if(!do.all.sensitivity.tests) Scens=Scens%>%filter(Scenario=='S1')
    Scens=Scens%>%
      mutate(Species=capitalize(Neim))
    Store.sens=vector('list',nrow(Scens))
    names(Store.sens)=Scens$Scenario
    Out.Scens=Scens
    Out.estimates=Out.likelihoods=Out.quantities=Out.rel.biom=Out.probs.rel.biom=Out.f.series=Out.probs.f.series=Out.B.Bmsy=
      Out.F.Fmsy=Out.probs.B.Bmsy=Out.Kobe.probs=store.warnings=store.convergence=vector('list',length(Store.sens))
    
    #Life history
    Life.history$Fecundity=ceiling(mean(Life.history$Fecundity))
    if(Neim%in%species.constant.fec)
    {
      Life.history$Fecu_a=NA
      Life.history$Fecu_b=NA
    }
    if(Neim%in%species.increase.terminal.age) Life.history$Max.age.F=max(Life.history$Max.age.F)
    Life.history$Max.age.F=ceiling(mean(Life.history$Max.age.F))
    Life.history$Breed.cycle=mean(Life.history$Breed.cycle)
    
    #Likelihood lambdas
    Lamdas.SS.lambdas=Life.history$SS_lambdas
    if(!is.null(Lamdas.SS.lambdas))
    {
      Lamdas.SS.lambdas$fleet=match(Lamdas.SS.lambdas$fleet,flitinfo$fleetname)
    }
    
    #Rename fleets following SS nomenclature
    names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))]=match(names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))],names(Flits))
    names(ktch.zone)[which(!names(ktch.zone)%in%c("SPECIES","Name","finyear"))]=match(names(ktch.zone)[which(!names(ktch.zone)%in%c("SPECIES","Name","finyear"))],names(Flits.zone))
    
    #Future catches/F  
    if("SS"%in%future.models)
    {
      #zones together
      NN=nrow(ktch)
      add.ct.future=ktch[1:years.futures,]
      add.ct.future$finyear=ktch$finyear[NN]+(1:years.futures)
      if("Survey"%in%names(Flits))
      {
        lef.flits=match(Flits[-match("Survey",names(Flits))],colnames(add.ct.future))
      }else
      {
        lef.flits=match(Flits,colnames(add.ct.future))
      }
      for(lf in 1:length(lef.flits))
      {
        add.ct.future[,lef.flits[lf]]=mean(unlist(ktch[(NN-years.futures+1):NN,lef.flits[lf]]),na.rm=T)
      }
      
      #by zone
      NN=nrow(ktch.zone)
      add.ct.future.zone=ktch.zone[1:years.futures,]
      add.ct.future.zone$finyear=ktch.zone$finyear[NN]+(1:years.futures)
      if("Survey"%in%names(Flits.zone))
      {
        lef.flits.zone=match(Flits.zone[-match("Survey",names(Flits.zone))],colnames(add.ct.future.zone))
      }else
      {
        lef.flits.zone=match(Flits.zone,colnames(add.ct.future.zone))
      }
      for(lf in 1:length(lef.flits.zone))
      {
        add.ct.future.zone[,lef.flits.zone[lf]]=mean(unlist(ktch.zone[(NN-years.futures+1):NN,lef.flits.zone[lf]]),na.rm=T)
      }
    }
    
    #Execute SS 
    for(s in 1:1)
    {
      
      #cpues
      Abund.single.area=Abundance.SS.format
      Abund.areas.as.fleets=Abundance.SS.format.zone
      
      #remove monthly years  
      if('drop.monthly.cpue.min'%in%names(Scens))
      {
        if(!is.na(Scens$drop.monthly.cpue.min[s]))
        {
          drop.yrS=Scens$drop.monthly.cpue.min[s]:Scens$drop.monthly.cpue.max[s]
          #zones combined
          if(!is.null(Abund.single.area))
          {
            Abund.single.area=Abund.single.area%>%
              mutate(dummy = rownames(Abund.single.area))%>%
              filter(!(grepl('TDGDLF.monthly',dummy) & Year%in%drop.yrS))%>%
              dplyr::select(-dummy)
          }
          #by zone
          if(!is.null(Abund.areas.as.fleets))
          {
            Abund.areas.as.fleets=Abund.areas.as.fleets%>%
              mutate(dummy = rownames(Abund.areas.as.fleets))%>%
              filter(!(grepl('TDGDLF.monthly',dummy) & Year%in%drop.yrS))%>%
              dplyr::select(-dummy)
          }
        }
      }
      
      #remove daily years  
      if(!is.na(Scens$Daily.cpues[s]))
      {
        if('TDGDLF.daily'%in%names(CPUE) | any(grep('TDGDLF.daily',names(CPUE.zone))))
        {
          rid.of=as.numeric(unlist(str_split(Scens$Daily.cpues[s], "&")))
          
          #zones combined
          if(!is.null(Abund.single.area))
          {
            drop.dis=Abund.single.area%>%
              mutate(this=grepl('TDGDLF.daily',rownames(Abund.single.area)) & Year%in%rid.of)
            if(any(drop.dis$this)) Abund.single.area=Abund.single.area[-which(drop.dis$this),]
          }
          
          #by zone
          if(!is.null(Abund.areas.as.fleets))
          {
            drop.dis=Abund.areas.as.fleets%>%
              mutate(this=grepl('TDGDLF.daily',rownames(Abund.areas.as.fleets)) & Year%in%rid.of)
            if(any(drop.dis$this)) Abund.areas.as.fleets=Abund.areas.as.fleets[-which(drop.dis$this),]
          }
        }
      }
      
      #remove TDGDLF CPUE 
      if(Scens$CPUE[s]=='No') 
      {
        Abund.single.area=Abund.areas.as.fleets=NULL
      }
      if('test.use.cpue'%in%names(Scens[s,]))
      {
        if(!(Scens$test.use.cpue[s]=='Yes'))
        {
          if(any(grepl("TDGDLF",names(CPUE))))
          {
            drop.dis.rows=grep("TDGDLF",row.names(Abund.single.area))
            if(length(drop.dis.rows)>0) Abund.single.area=Abund.single.area[-drop.dis.rows,]
          }
          if(any(grepl("TDGDLF",names(CPUE.zone))))
          {
            drop.dis.rows=grep("TDGDLF",row.names(Abund.areas.as.fleets))
            if(length(drop.dis.rows)>0) Abund.areas.as.fleets=Abund.areas.as.fleets[-drop.dis.rows,] 
          }
        }
        if(!is.null(Abund.single.area)) if(nrow(Abund.single.area)==0) Abund.single.area=NULL
        if(!is.null(Abund.areas.as.fleets)) if(nrow(Abund.areas.as.fleets)==0) Abund.areas.as.fleets=NULL
        
      }
      
      #remove Length.comps
      Size.comp.single.area=Size.compo.SS.format
      Size.comp.areas.as.fleets=Size.compo.SS.format.zone
      if(Scens$Length.comps[s]=='No') 
      {
        Size.comp.single.area=Size.comp.areas.as.fleets=NULL
      }
      
      #remove Mean.body
      Meanbodywt.single.area=meanbodywt.SS.format
      Meanbodywt.areas.as.fleets=meanbodywt.SS.format.zone
      if(Scens$Mean.body[s]=='No') 
      {
        Meanbodywt.single.area=Meanbodywt.areas.as.fleets=NULL
      }
      
      #use cpue or length comp in likelihood?
      #note: superseded
      if(Life.history$drop.length.comp)
      {
        Size.compo.SS.format=NULL
        Size.compo.SS.format.zone=NULL
      }
      if(Life.history$drop.cpue)
      {
        Abundance.SS.format=NULL
        Abundance.SS.format.zone=NULL
      }
      
      #Set future F instead of catch if required by scenario  
      if(exists('add.ct.future'))   
      {
        #zones combined
        add.ct.or.F_future=add.ct.future
        if(Scens[s,'Forecasting']=="F")
        {
          Nms=names(add.ct.or.F_future)
          mutate.these.fleets=subset(Nms,!Nms%in%c("SPECIES","Name","finyear"))
          id.mutate.these.fleets=which(!Nms%in%c("SPECIES","Name","finyear"))
          for(mu in 1:length(mutate.these.fleets))
          {
            id.mu=match(names(Life.history$F.forecasting.value)[mu],flitinfo$fleetname)
            F.val=Life.history$F.forecasting.value[id.mu]
            mu.id=match(id.mu,Nms)
            add.ct.or.F_future[,mu.id]=ifelse(add.ct.or.F_future[,mu.id]>0,F.val,0)
          }
        }
        
        #by zone 
        add.ct.or.F_future.zone=add.ct.future.zone
        if(Scens[s,'Forecasting']=="F")
        {
          Nms=names(add.ct.or.F_future.zone)
          mutate.these.fleets=subset(Nms,!Nms%in%c("SPECIES","Name","finyear"))
          id.mutate.these.fleets=which(!Nms%in%c("SPECIES","Name","finyear"))
          for(mu in 1:length(mutate.these.fleets))
          {
            id.mu=match(names(Life.history$F.forecasting.value)[mu],flitinfo$fleetname)
            F.val=Life.history$F.forecasting.value[id.mu]
            mu.id=match(id.mu,Nms)
            add.ct.or.F_future.zone[,mu.id]=ifelse(add.ct.or.F_future.zone[,mu.id]>0,F.val,0)
          }
        }
      }
      
      #a. Create SS input files   
      if(create.SS.inputs)
      {
        #a.1 Specify ktch,flitinfo,abundance,comps,mean.weight,var.adjs & future based on scenario   
        if(Scens$Spatial[s]=='single area')  
        {
          KAtch=ktch
          FLitinFO=flitinfo
          Abund=Abund.single.area
          Size.com=Size.comp.single.area
          Size.com_all=dummy.Size.compo.SS.format.all
          meanbody=Meanbodywt.single.area
          tags=Tags.SS.format 
          Var.ad=Var.ad.factr
          add.future=add.ct.or.F_future
          KAtch.ret.disc=NULL
          if(!is.null(retained.discarded.ktch))   
          {
            KAtch.ret.disc=retained.discarded.ktch%>%
              ungroup()%>%
              left_join(FLitinFO%>%dplyr::select(fleetname)%>%mutate(fleet=1:nrow(FLitinFO)),
                        by=c('Fishry'='fleetname'))%>%
              dplyr::select(-c('SPECIES','Name','Fishry'))%>%
              mutate(seas=1,stderr=CV.discards)%>%
              rename(yr=finyear,
                     obs=Tonnes)%>%
              relocate(yr,seas,fleet,obs,stderr)%>%
              arrange(fleet,yr)%>%
              data.frame
            Length.limit=Life.history$SS_retention$P_5
            discard.flits=unique(KAtch.ret.disc$fleet)
            
            if(!is.null(Size.com))
            {
              id.rel.cols=names(Size.com)[grep(paste(c('f','m'),collapse = '|'),names(Size.com))]
              id.rel.cols=unique(as.numeric(gsub("[a-zA-Z]", "", subset(id.rel.cols,!id.rel.cols%in%c("Nsamp")))))
              from.to=seq(id.rel.cols[which.min(abs(id.rel.cols - Length.limit))],max(id.rel.cols),by=TL.bins.cm)
              id.rel.cols=grep(paste(from.to,collapse = '|'),names(Size.com))
              id.irrel.cols=grep(paste(names(Size.com)[-c(which(c("year","Seas","Fleet","Sex","Part","Nsamp")%in%
                                                                  names(Size.com)),id.rel.cols)],collapse = '|'),
                                 names(Size.com))
              Size.com.discards=Size.com%>%filter(Fleet%in%discard.flits)
              Size.com=Size.com%>%filter(!Fleet%in%discard.flits)
              Size.com.discards_retained=Size.com.discards
              Size.com.discards_discarded=Size.com.discards
              Size.com.discards_retained[,id.rel.cols]=0
              Size.com.discards_retained$Part=2
              Size.com.discards_discarded[,id.irrel.cols]=0
              Size.com.discards_discarded$Part=1
              Size.com.discards=Size.com.discards_retained
              drop.yrs=rowSums(Size.com.discards_discarded[,id.rel.cols])
              drop.yrs=which(drop.yrs<Min.annual.obs.ktch)
              if(length(drop.yrs)>0)
              {
                Size.com.discards_discarded=Size.com.discards_discarded[-drop.yrs,]
              }
              if(nrow(Size.com.discards_discarded)>0)
              {
                Size.com.discards=rbind(Size.com.discards_retained,Size.com.discards_discarded)
              }
              Size.com=rbind(Size.com,Size.com.discards)%>%
                arrange(Fleet,year)
            }
            
            if(!is.null(meanbody))
            {
              meanbody=meanbody%>%
                mutate(part=ifelse(fleet%in%discard.flits,2,part)) 
            }
          }
        }
        if(Scens$Spatial[s]=='areas-as-fleets')   
        {
          KAtch=ktch.zone
          FLitinFO=flitinfo.zone
          Abund=Abund.areas.as.fleets
          Size.com=Size.comp.areas.as.fleets
          Size.com_all=dummy.Size.compo.SS.format.all_zone
          meanbody=Meanbodywt.areas.as.fleets
          tags=Tags.SS.format.zone  
          Var.ad=Var.ad.factr.zone
          add.future=add.ct.or.F_future.zone
          KAtch.ret.disc=NULL
          if(!is.null(retained.discarded.ktch.zone))
          {
            KAtch.ret.disc=retained.discarded.ktch.zone%>%
              ungroup()%>%
              left_join(FLitinFO%>%dplyr::select(fleetname)%>%mutate(fleet=1:nrow(FLitinFO)),
                        by=c('Fishry'='fleetname'))%>%
              dplyr::select(-c('SPECIES','Name','Fishry'))%>%
              mutate(seas=1,stderr=CV.discards)%>% 
              rename(yr=finyear,
                     obs=Tonnes)%>%
              relocate(yr,seas,fleet,obs,stderr)%>%
              arrange(fleet,yr)%>%
              data.frame
            
            Length.limit=Life.history$SS_retention$P_5
            discard.flits=unique(KAtch.ret.disc$fleet)
            if(!is.null(Size.com))
            {
              id.rel.cols=names(Size.com)[grep(paste(c('f','m'),collapse = '|'),names(Size.com))]
              id.rel.cols=unique(as.numeric(gsub("[a-zA-Z]", "", subset(id.rel.cols,!id.rel.cols%in%c("Nsamp")))))
              from.to=seq(id.rel.cols[which.min(abs(id.rel.cols - Length.limit))],max(id.rel.cols),by=TL.bins.cm)
              id.rel.cols=grep(paste(from.to,collapse = '|'),names(Size.com))
              id.irrel.cols=grep(paste(names(Size.com)[-c(which(c("year","Seas","Fleet","Sex","Part","Nsamp")%in%
                                                                  names(Size.com)),id.rel.cols)],collapse = '|'),
                                 names(Size.com))
              Size.com.discards=Size.com%>%filter(Fleet%in%discard.flits)
              Size.com=Size.com%>%filter(!Fleet%in%discard.flits)
              Size.com.discards_retained=Size.com.discards
              Size.com.discards_discarded=Size.com.discards
              Size.com.discards_retained[,id.rel.cols]=0
              Size.com.discards_retained$Part=2
              Size.com.discards_discarded[,id.irrel.cols]=0
              Size.com.discards_discarded$Part=1
              Size.com.discards=Size.com.discards_retained
              drop.yrs=rowSums(Size.com.discards_discarded[,id.rel.cols])
              drop.yrs=which(drop.yrs<Min.annual.obs.ktch.zone)
              if(length(drop.yrs)>0)
              {
                Size.com.discards_discarded=Size.com.discards_discarded[-drop.yrs,]
              }
              if(nrow(Size.com.discards_discarded)>0)
              {
                Size.com.discards=rbind(Size.com.discards_retained,Size.com.discards_discarded)
              }
              Size.com=rbind(Size.com,Size.com.discards)%>%
                arrange(Fleet,year)
            }
            
            if(!is.null(meanbody))
            {
              meanbody=meanbody%>%
                mutate(part=ifelse(fleet%in%discard.flits,2,part))
            }
            
          }
          
          #Change species specific sel pars for spatial model
          if(Neim=="sandbar shark" & !is.null(Life.history$SS_offset_selectivity))   
          {
            Life.history$SS_offset_selectivity=Life.history$SS_offset_selectivity%>%
              mutate(P_1=ifelse(Fleet=='Survey',6.5,P_1),
                     P_3=ifelse(Fleet=='Survey',0.22,P_3),
                     P_4=ifelse(Fleet=='Survey',0.09,P_4))
          }
        }
        
        #a.2 Indo IUU - F estimation
        if(Scens$Estim.Indo.IUU[s]=="Yes")
        {
          #get Indo catch
          Indo.ktch=KtCh%>%
            filter(Name==Neim & Data.set=='Indonesia')
          Indo.ktch.mean.future=mean(Indo.ktch$LIVEWT.c[(nrow(Indo.ktch)-years.futures):nrow(Indo.ktch)])
          Indo.ktch.years=Indo.ktch%>%filter(finyear%in% Indo.years.sel)   #select some years
          
          #remove Indo catch from 'other'
          id.flit.other=match('Other',FLitinFO$fleetname)
          id.yrs.indo=match(Indo.ktch$finyear,KAtch$finyear)
          
          KAtch[id.yrs.indo,match(id.flit.other,colnames(KAtch))]=KAtch[id.yrs.indo,match(id.flit.other,colnames(KAtch))]-Indo.ktch$LIVEWT.c
          add.future[,match(id.flit.other,colnames(add.future))]=add.future[,match(id.flit.other,colnames(add.future))]-Indo.ktch.mean.future
          
          #id original fleets with data 
          if(!is.null(Abund))
          {
            id.abun.flit=FLitinFO[sort(unique(Abund$index)),]%>%
              mutate(old.fleet=sort(unique(Abund$index)),new.fleet=NA)
          }
          if(!is.null(Size.com))
          {
            id.Size.com.flit=FLitinFO[sort(unique(Size.com$Fleet)),]%>%
              mutate(old.fleet=sort(unique(Size.com$Fleet)),new.fleet=NA)
          }
          if(!is.null(meanbody))
          {
            id.meanbody.flit=FLitinFO[sort(unique(meanbody$fleet)),]%>%
              mutate(old.fleet=sort(unique(meanbody$fleet)),new.fleet=NA)
            
          }
          if(!is.null(tags))
          {
            id.tag.flit=FLitinFO[sort(unique(tags$recaptures$Fleet)),]%>%
              mutate(old.fleet=sort(unique(tags$recaptures$Fleet)),new.fleet=NA)
          }
          
          #add Indo as separate fleet
          indo.fleet=FLitinFO%>%filter(fleetname=='Other')%>%mutate(fleetname='Indo.IUU')
          FLitinFO=add_row(FLitinFO, indo.fleet, .after = id.flit.other)
          
          #catch history
          Indo.dummy=KAtch%>%
            dplyr::select(finyear,as.character(id.flit.other))%>%
            rename(Indo=as.character(id.flit.other))%>%
            mutate(Indo=0)
          Indo.dummy$Indo[match(Indo.ktch$finyear,Indo.dummy$finyear)]=Indo.ktch$LIVEWT.c
          
          #Set  unknown catches to ball park      
          if(set.indo.catches.for.unknown.years)
          {
            id.unknown.yrs=Indo.dummy$finyear[which(Indo.dummy$finyear%in%indo.unknown.catch.years)]
            find.comparable.years=Indo_apprehensions%>%filter(year%in%id.unknown.yrs)
            all_values=Indo_apprehensions%>%filter(!year%in%indo.unknown.catch.years)
            for(fi in 1:nrow(find.comparable.years))
            {
              iii=find.comparable.years$year[fi]
              ref_value=Indo_apprehensions%>%filter(year==iii)%>%pull(Apprehensions)
              dumi.ktch=Indo.dummy%>%
                filter(finyear==all_values$year[which.min(abs(all_values$Apprehensions - ref_value))])%>%
                pull(Indo)
              Indo.dummy$Indo[match(iii,Indo.dummy$finyear)]=dumi.ktch
            }
          }
          
          #Set catches to 'unknown' except for some years
          if(set.indo.catches.to.unknown) Indo.dummy=Indo.dummy%>%mutate(Indo=-999)
          
          #Set catches to very low
          if(set.indo.catches.to.very.low)  Indo.dummy=Indo.dummy%>%mutate(Indo=0.001) 
          
          #keep some years catch
          if(keep.some.Indo.yrs)
          {
            Indo.dummy=Indo.dummy%>%
              mutate(Indo=replace(Indo,match(Indo.ktch.years$finyear,KAtch$finyear),Indo.ktch.years$LIVEWT.c))
          }
          
          KAtch=KAtch%>%mutate(Indo=Indo.dummy$Indo,.after = as.character(id.flit.other))
          id.Kls=match('Indo',names(KAtch)):ncol(KAtch)
          names(KAtch)[id.Kls]=seq(from = match('Indo.IUU',FLitinFO$fleetname), length.out = length(id.Kls))
          
          
          #future catch
          Indo.dummy=add.future%>%
            dplyr::select(as.character(id.flit.other))%>%
            rename(Indo=as.character(id.flit.other))%>%
            mutate(Indo=Indo.ktch.mean.future)
          add.future=add.future%>%mutate(Indo=Indo.dummy$Indo,.after = as.character(id.flit.other))
          id.Kls=match('Indo',names(add.future)):ncol(add.future)
          names(add.future)[id.Kls]=seq(from = match('Indo.IUU',FLitinFO$fleetname), length.out = length(id.Kls))
          
          #Abund  
          if(!is.null(Abund))
          {
            id.abun.flit$new.fleet=FLitinFO%>%
              mutate(row_id = row_number())%>%
              filter(fleetname%in%id.abun.flit$fleetname)%>%
              pull(row_id)
            id.abun.flit=id.abun.flit%>%dplyr::select(old.fleet,new.fleet)
            Abund=Abund%>%
              rownames_to_column(var = "row_names")%>%
              left_join(id.abun.flit,by=c('index'='old.fleet'))%>%
              mutate(index=new.fleet)%>%
              dplyr::select(-new.fleet)%>%
              column_to_rownames(var = "row_names")
          }
          #Size.com
          if(!is.null(Size.com))
          {
            id.Size.com.flit$new.fleet=FLitinFO%>%
              mutate(row_id = row_number())%>%
              filter(fleetname%in%id.Size.com.flit$fleetname)%>%
              pull(row_id)
            id.Size.com.flit=id.Size.com.flit%>%dplyr::select(old.fleet,new.fleet)
            Size.com=Size.com%>%
              left_join(id.Size.com.flit,by=c('Fleet'='old.fleet'))%>%
              mutate(Fleet=new.fleet)%>%
              dplyr::select(-new.fleet)
            
          }
          #meanbody
          if(!is.null(meanbody))
          {
            id.meanbody.flit$new.fleet=FLitinFO%>%
              mutate(row_id = row_number())%>%
              filter(fleetname%in%id.meanbody.flit$fleetname)%>%
              pull(row_id)
            id.meanbody.flit=id.meanbody.flit%>%dplyr::select(old.fleet,new.fleet)
            meanbody=meanbody%>%
              left_join(id.meanbody.flit,by=c('fleet'='old.fleet'))%>%
              mutate(fleet=new.fleet)%>%
              dplyr::select(-new.fleet)
          }
          #Tags
          if(!is.null(tags))
          {
            id.tag.flit$new.fleet=FLitinFO%>%
              mutate(row_id = row_number())%>%
              filter(fleetname%in%id.tag.flit$fleetname)%>%
              pull(row_id)
            id.tag.flit=id.tag.flit%>%
              dplyr::select(old.fleet,new.fleet)
            
            tags$Initial.reporting.rate=tags$Initial.reporting.rate%>%
              left_join(id.tag.flit,by=c('Fleet'='old.fleet'))%>%
              dplyr::select(-Fleet)%>%
              rename(Fleet=new.fleet)
            tags$recaptures=tags$recaptures%>%
              left_join(id.tag.flit,by=c('Fleet'='old.fleet'))%>%
              mutate(Fleet=new.fleet)%>%
              dplyr::select(-new.fleet)
          }
          #Variance adjustment
          if(any(!is.null(Abund),!is.null(Size.com),!is.null(meanbody)))
          {
            New.var.adj=rbind(id.abun.flit,id.Size.com.flit,id.meanbody.flit)%>%
              distinct(old.fleet,.keep_all = T)
            Var.ad=Var.ad%>%
              left_join(New.var.adj,by=c('Fleet'='old.fleet'))%>%
              mutate(Fleet=new.fleet)%>%
              dplyr::select(-new.fleet)
          }
          
          #Add Apprehensions as index of abundance
          if(!is.null(Abund))
          {
            Indo.abund=Indo_apprehensions%>%
              rename(Year=year,
                     Mean=Apprehensions)%>%
              filter(Year<=max(Abund$Year))%>%
              filter(Year%in%Indo.years.cpue)%>%
              mutate(seas=1,
                     CV=CV_apprehensions,
                     index=match('Indo.IUU',FLitinFO$fleetname))%>%
              dplyr::select(names(Abund))
            if(scale.Indo.appre) Indo.abund$Mean=Indo.abund$Mean/mean(Indo.abund$Mean,na.rm=T)
            rownames(Indo.abund)=paste('Indo',1:nrow(Indo.abund),sep='.')
            Abund=rbind(Abund,Indo.abund)%>%
              arrange(index,Year)
            
          }
          
          #Add Q
          Life.history$Q.inits=Life.history$Q.inits%>%
            filter(Fleet%in%FLitinFO$fleetname)%>%
            mutate(Fleet.n.new=match(Fleet,FLitinFO$fleetname))
          Indo.Q=Life.history$Q.inits%>%
            filter(Fleet=='Other')%>%
            mutate(Fleet='Indo.IUU',
                   Fleet.n.new=match('Indo.IUU',FLitinFO$fleetname))
          Life.history$Q.inits=rbind(Life.history$Q.inits,Indo.Q)%>%
            mutate(Fleet.n=Fleet.n.new)%>%
            arrange(Fleet.n)%>%
            dplyr::select(-Fleet.n.new)%>%
            distinct(Fleet,.keep_all = TRUE)
        }
        
        #a.3 Indo IUU- Test effect of using only Apprehensions for catch reconstruction
        if(Scens$Test.Indo.IUU.catch[s]=='Yes')  
        {
          Indo.IUU=Indo.IUU.apprehensions%>%
            filter(SPECIES==List.sp[[i]]$Species)
          Indo.IUU.app=Indo.IUU.apprehensions.only%>%
            filter(SPECIES==List.sp[[i]]$Species)
          des.flt=match(match('Other',FLitinFO$fleetname),colnames(KAtch))
          des.yrs=match(as.numeric(substr(Indo.IUU$FINYEAR,1,4)),KAtch$finyear)
          new.ktch=unlist(KAtch[des.yrs,des.flt])-(Indo.IUU$LIVEWT.c/1000)+(Indo.IUU.app$LIVEWT.c/1000)
          
          KAtch[des.yrs,c(3,des.flt)]%>%
            data.frame()%>%
            mutate(Apprehensions.only=new.ktch)%>%
            rename(max.Apprehensions.Forfeitures=X2)%>%
            gather(Method,Tons,-finyear)%>%
            ggplot(aes(finyear,Tons,color=Method))+
            geom_point()+geom_line()+ylim(0,NA)+
            theme_PA()+
            theme(legend.position = 'top',
                  legend.title=element_blank())
          ggsave(paste0(this.wd,"/Catch Other fleet_Indo catch recons_Apprehensions or Forfeitures.tiff"),width=7,height=6,compression = "lzw")
          
          KAtch[des.yrs,des.flt]=new.ktch
        }
        
        #a.4 set MainRdevYrFirst   
        #note: align with data-rich years (cpue, comps, meanbody, etc)
        Abund1=Abund
        if(!is.null(Abund1)) Abund1=Abund1%>%rename_with(tolower)
        Max.yr.obs=max(unlist(lapply(list(Abund1,Size.com,meanbody),function(x) if(!is.null(x))max(x$year))))
        Life.history$MainRdevYrLast=min(Max.yr.obs,max(KAtch$finyear,na.rm=T))
        
        Min.yr.obs=min(unlist(lapply(list(Abund1,Size.com,meanbody),function(x) if(!is.null(x))min(x$year))))
        if(Life.history$First.yr.main.rec.dev=='min.obs')
        {
          Min.yR=Min.yr.obs
          if(Life.history$First.yr.main.rec.dev_buffer)  Min.yR=Min.yr.obs-round(min(Life.history$Age.50.mat))
          MainRdevYrFirst=Min.yR
        }
        if(Life.history$First.yr.main.rec.dev=='min.ktch') MainRdevYrFirst=min(ktch$finyear)
        Life.history$MainRdevYrFirst=MainRdevYrFirst
        
        #a.5 need to reset rec pars for tuning
        #if(Scens$Scenario[s]=='S1' & Tune.SS.model)  #set to NULL in exploratory phase to fully see effect of data
        {
          Life.history$recdev_early_start=0
          Life.history$SR_sigmaR=0.2
          Life.history$RecDev_Phase=3
          
          #Ramp:
          # The model linearly interpolates the adjustment fraction between these four year-markers
          
          #The last year of the early recruitment period where no bias adjustment (0%) is applied. 
          # Typically used for very early years with no data
          Life.history$last_early_yr_nobias_adj_in_MPD=MainRdevYrFirst-1
          
          #The year when the model transitions to full bias adjustment (100%). 
          # This should align with the start of informative composition data
          Life.history$first_yr_fullbias_adj_in_MPD=MainRdevYrFirst 
          
          #The last year where full bias adjustment is applied
          # This usually marks the point where recent data (like small fish in length comps) still strongly informs recruitment
          Life.history$last_yr_fullbias_adj_in_MPD=Max.yr.obs-2
          
          Life.history$first_recent_yr_nobias_adj_in_MPD=Max.yr.obs-1
          
          #The maximum fraction of the bias adjustment to apply (typically set to 1.0 for full correction). 
          # If set lower, the model never applies the full theoretical correction
          Life.history$max_bias_adj_in_MPD=1
          
          #Comps variance adjustment
          Var.ad=NULL
        }
        
        #a.6 remove 2001 and 2002 from survey as too many juveniles, shots done all over the place 
        if(!is.null(Size.com))
        {
          id.survey=match('Survey',FLitinFO$fleetname)
          if(!is.na(id.survey)) Size.com=Size.com%>%filter(!(Fleet==id.survey & year%in%c(2001,2002)))
        }
        
        
        
        #a.8 remove Tag data is tested in scenario
        if(Scens$Tagging[s]=='No')
        {
          tags=NULL
        }
        #a.9 remove selectivity offsets 
        Life.history1=Life.history
        if(Scens$Use.male.sel.offset[s]=='No' & any(grepl('offset_selectivity',names(Life.history1))))
        {
          Life.history1=Life.history1[-grep('offset_selectivity',names(Life.history1))]
        }
        
        #a.10 remove growth estimation    
        if(isTRUE(Life.history1$SS3.estim.growth.pars))
        {
          if(Scens$estim.growth[s]=='No')  Life.history1$SS3.estim.growth.pars=FALSE
        }
      }
    }  #end s
  }
  
}# end i species loop

#Tweak inputs
Life.history2=Life.history1
#Life.history2=Life.history2[-match('SS_offset_selectivity',names(Life.history2))]
Abund2=Abund #%>%filter(index%in%c(4,7))
#Abund2=NULL
Size.com2=Size.com
#Size.com2=Size.com%>%filter(Fleet%in%c(5,8)) 
meanbody2=meanbody
#meanbody2=NULL
tags2=tags
tags2=NULL

this.wd2=this.wd1
#this.wd2=paste(this.wd,'S2',sep='/')
#this.wd2='C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/New folder/3. Sel body Cpue Zn1'


fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SS'),   
             new.path=this.wd2,
             Scenario=Scens[s,]%>%mutate(Model='SS'),
             Catch=KAtch,
             Catch.ret.disc=KAtch.ret.disc,
             life.history=Life.history2,
             depletion.yr=NULL,
             fleets=names(KAtch)[which(!names(KAtch)%in%c("SPECIES","Name","finyear"))],
             fleetinfo=FLitinFO,
             abundance=Abund2,   
             size.comp=Size.com2,
             meanbodywt=meanbody2,
             Tags=tags2,
             F.tagging=F.SS.format,
             cond.age.len=Cond.age.len.SS.format,
             MeanSize.at.Age.obs=MeanSize.at.Age.obs.SS.format,
             Lamdas=Lamdas.SS.lambdas,
             Var.adjust.factor=Var.ad,    
             Future.project=add.future) 

#2. Change phases of some pars
do.dis=FALSE
if(do.dis)
{
  start <- r4ss::SS_readstarter(file = file.path(this.wd2, "starter.ss"), verbose = FALSE)
  dat <- r4ss::SS_readdat(file = file.path(this.wd2, start$datfile), verbose = FALSE)
  ctl <- r4ss::SS_readctl(file = file.path(this.wd2, start$ctlfile), verbose = FALSE, use_datlist = TRUE, datlist = dat)
  
  only.esti.sel.pars=FALSE
  if(only.esti.sel.pars)
  {
    ctl$SR_parms[match('SR_LN(R0)',rownames(ctl$SR_parms)),c('INIT','PHASE')]=c(7.4,-1)
    ctl$size_selex_parms[grep('P_1_Southern.shark_1_Zone2',rownames(ctl$size_selex_parms)),c('PHASE')]=1
    
  }
  
  #ctl$size_selex_parms[grepl('PMalOff',rownames(ctl$size_selex_parms)),c('PHASE')]=-5
  #ctl$size_selex_parms[grepl('PMalOff_1',rownames(ctl$size_selex_parms)),c('PHASE')]=5
  #ctl$size_selex_parms[grepl('PMalOff_5',rownames(ctl$size_selex_parms)),c('PHASE')]=5
  
  #ctl$size_selex_parms[grep('P_1_Southern.shark_1_West',rownames(ctl$size_selex_parms)),c('PHASE')]=2  #why all phases <0?
  
  r4ss::SS_writectl(ctl, outfile = file.path(this.wd2, start$ctlfile), overwrite = TRUE, verbose = FALSE)
  
}


#----------- run and plot -----------------
#this.wd2=paste(this.wd,'S4',sep='/')

fn.run.SS(where.inputs=this.wd2,  where.exe=Where.exe, args=Arg)
COVAR=FORECAST=FALSE
if(Arg=="") COVAR=TRUE
if("SS"%in%future.models) FORECAST=TRUE
Report=SS_output(this.wd2,covar=COVAR,forecast=FORECAST,readwt=F)
this.plot=1:26 
SS_plots(Report,plot=this.plot,  png=T)



#----------- check selectivities -----------------
if(!exists('doubleNorm24.fn')) fn.source1("SS_selectivity functions.R")
x=seq(30,200,5)

Mod1=data.frame(p1=108.22,p2=-20.0,p3=5.322,p4=5.60,Name='Kirkwood & Walker')  
Mod2=data.frame(p1=110.2,p2=-11.11,p3=4.65,p4=6.93,Name='Southern 1 West')   
Mod3=data.frame(p1=106.9,p2=-19.99,p3=4.38,p4=6.445,Name='Southern 1 Zone1')
Mod4=data.frame(p1=118.14,p2=-19.99,p3=5.026,p4=5.8759,Name='Southern 1 Zone2')
Mod5=data.frame(p1=109.46,p2=-11,p3=4.78,p4=6.5,Name='Southern 2 Zone1')
Mod6=data.frame(p1=112.17,p2=-11,p3=4.41,p4=6.5,Name='Southern 2 Zone2')
Sel.mod1=with(Mod1,doubleNorm24.fn(x,a=p1,b=p2, c=p3, d=p4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))
Sel.mod2=with(Mod2,doubleNorm24.fn(x,a=p1,b=p2, c=p3, d=p4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))
Sel.mod3=with(Mod3,doubleNorm24.fn(x,a=p1,b=p2, c=p3, d=p4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))
Sel.mod4=with(Mod4,doubleNorm24.fn(x,a=p1,b=p2, c=p3, d=p4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))
Sel.mod5=with(Mod5,doubleNorm24.fn(x,a=p1,b=p2, c=p3, d=p4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))
Sel.mod6=with(Mod6,doubleNorm24.fn(x,a=p1,b=p2, c=p3, d=p4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))

male.offset=Report$sizeselex%>%
  filter(Factor=='Lsel' & Yr=='2023')%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))%>%
  filter(Fleet==5 & Sex==2)

plot(x,Sel.mod1,type='l',lwd=2)
lines(x,Sel.mod2,col=2,lwd=2,lty=2)
lines(x,Sel.mod3,col=3,lwd=2,lty=2)
lines(x,Sel.mod4,col=4,lwd=2,lty=2)
lines(x,Sel.mod5,col=5,lwd=2,lty=3)
lines(x,Sel.mod6,col=6,lwd=2,lty=3)
lines(male.offset$TL,male.offset$Sel,col=7)
legend('left',unique(c(Mod1$Name,Mod2$Name,Mod3$Name,
                       Mod4$Name,Mod5$Name,Mod6$Name,'Southern1.zone2_male')),
       col=1:7,lty=c(1,2,2,2,3,3),bty='n')



a=Report$sizeselex%>%
  filter(Factor=='Lsel' & Yr=='2023')%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Fleet2=Fleet,
         Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))

# Report2022=SS_output(paste(HandL.out,capitalize(Neim),"/",2022,"/SS3 integrated/S1",sep=''),
#                       covar=COVAR,forecast=FORECAST,readwt=F)
b=Report2022$sizeselex%>%
  filter(Factor=='Lsel' & Yr==max(Yr))%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))
b.zn1=b.zn2=b%>%filter(Fleet==3)
b.zn1$Fleet=4
b.zn2$Fleet=5
b.2_zn1=b.2_zn2=b%>%filter(Fleet==4)
b.2_zn1$Fleet=7
b.2_zn2$Fleet=8
b=b%>%mutate(Fleet=ifelse(Fleet==4,6,Fleet))
b=rbind(b,b.zn1,b.zn2,b.2_zn1,b.2_zn2)%>%
  mutate(Fleet2=Fleet,
         Fleet=paste0(Fleet,'_2022'))

rbind(a,
      #data.frame(Fleet='K&W',Yr=2023,Sex=1,TL=x,Sel=Sel.mod1)%>%mutate(Fleet2=Fleet),
      b)%>%
  filter(!Fleet2%in%c('1','2'))%>%
  ggplot(aes(TL,Sel,color=Fleet,linetype=Sex))+geom_point(aes(shape=Sex))+
  geom_line()+theme_PA()+facet_wrap(~Fleet2)



#----------- compare cpues ------------------
CPUE=compact(Catch.rate.series[[i]])
rbind(CPUE$TDGDLF.monthly%>%mutate(type='Monthly'),
      CPUE$TDGDLF.monthly.Zone1%>%mutate(type='Zone1'),
      CPUE$TDGDLF.monthly.Zone2%>%mutate(type='Zone2'))%>%
  ggplot(aes( yr.f,Mean,color=type))+
  geom_line()+geom_errorbar(aes(ymin = LOW.CI, ymax = UP.CI))+geom_point()+
  theme_PA()+ylim(0,NA)

rbind(CPUE$TDGDLF.daily%>%mutate(type='Daily'),
      CPUE$TDGDLF.daily.Zone1%>%mutate(type='Zone1'),
      CPUE$TDGDLF.daily.Zone2%>%mutate(type='Zone2'))%>%
  ggplot(aes( yr.f,Mean,color=type))+
  geom_line()+geom_errorbar(aes(ymin = LOW.CI, ymax = UP.CI))+geom_point()+
  theme_PA()+ylim(0,NA)


this.wd2022=paste(HandL.out,capitalize(Neim),"/",2022,"/SS3 integrated/S1",sep='')
start2022 <- r4ss::SS_readstarter(file = file.path(this.wd2022, "starter.ss"), verbose = FALSE)
Dat2022 <- r4ss::SS_readdat(file = file.path(this.wd2022, start2022$datfile), verbose = FALSE)

rbind(Dat2022$CPUE%>%
        rename(yr.f=year,Mean=obs,CV=se_log)%>%
        dplyr::select(yr.f,Mean,CV)%>%mutate(type='2022'),
      CPUE$TDGDLF.monthly%>%
        dplyr::select(yr.f,Mean,CV)%>%mutate(type='Monthly'),
      CPUE$TDGDLF.daily%>%
        dplyr::select(yr.f,Mean,CV)%>%mutate(type='Daily'))%>%
  ggplot(aes( yr.f,Mean,color=type))+
  geom_line()+geom_errorbar(aes(ymin = Mean-CV, ymax = Mean-CV))+geom_point()+
  theme_PA()+ylim(0,NA)



start2026 <- r4ss::SS_readstarter(file = file.path(this.wd2, "starter.ss"), verbose = FALSE)
Dat2026 <- r4ss::SS_readdat(file = file.path(this.wd2, start2026$datfile), verbose = FALSE)
Dat2026$CPUE%>%
  mutate(index=as.character(index))%>%
  ggplot(aes(year,obs,color=index))+
  geom_point()+
  geom_errorbar(aes(ymin = obs -   se_log, ymax = obs +se_log))

#-----------  Check Estimated selectivities (report)-----------------------------------------------------
xx="C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/New folder/S1"
Report_S1=SS_output(xx,covar=COVAR,forecast=FORECAST,readwt=F)
a=Report_S1$sizeselex%>%
  filter(Factor=='Lsel' & Yr=='2023')%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))

xx="C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/New folder/S4"
Report_S4=SS_output(xx,covar=COVAR,forecast=FORECAST,readwt=F)
a_S4=Report_S4$sizeselex%>%
  filter(Factor=='Lsel' & Yr=='2023')%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))

xx='C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/Population dynamics/1.Whiskery shark/2022/SS3 integrated/S1'
Report2022=SS_output(xx,covar=COVAR,forecast=FORECAST,readwt=F)
a2022=Report2022$sizeselex%>%
  filter(Factor=='Lsel' & Yr=='2026')%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))

Expnded.S4=a_S4%>%mutate(Fleet=case_when(Fleet==4~'6',TRUE~Fleet))
Expnded.S4=rbind(Expnded.S4,
                 Expnded.S4%>%filter(Fleet==3)%>%mutate(Fleet=4),
                 Expnded.S4%>%filter(Fleet==3)%>%mutate(Fleet=5),
                 Expnded.S4%>%filter(Fleet==6)%>%mutate(Fleet=7),
                 Expnded.S4%>%filter(Fleet==6)%>%mutate(Fleet=8))
Expnded.2022=a2022%>%mutate(Fleet=case_when(Fleet==4~'6',TRUE~Fleet))
Expnded.2022=rbind(Expnded.2022,
                   Expnded.2022%>%filter(Fleet==3)%>%mutate(Fleet=4),
                   Expnded.2022%>%filter(Fleet==3)%>%mutate(Fleet=5),
                   Expnded.2022%>%filter(Fleet==6)%>%mutate(Fleet=7),
                   Expnded.2022%>%filter(Fleet==6)%>%mutate(Fleet=8))


rbind(a%>%mutate(Ass='S1'),
      Expnded.S4%>%mutate(Ass='S4'),
      Expnded.2022%>%mutate(Ass='2022'))%>%
  ggplot(aes(TL,Sel,color=Ass,shape=Sex))+
  geom_point()+geom_line()+
  facet_wrap(~Fleet)+
  theme_PA()

#-----------  Check fixed Ro estim Sel--------
Report_fixedSel=SS_output("C:\\Users\\myb\\OneDrive - Department of Primary Industries And Regional Development\\Desktop\\test scenarios\\tunning_Whiskery\\S1_fixed Sel pars",
                          covar=COVAR,forecast=FORECAST,readwt=F)
a_fixedSel=Report_fixedSel$sizeselex%>%
  filter(Factor=='Lsel' & Yr=='2023')%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))



Report_estim.Len_Meanbody=SS_output("C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/New folder/2.Sel and meanbody",
                                    covar=COVAR,forecast=FORECAST,readwt=F)
a_estim.Len_Meanbody=Report_estim.Len_Meanbody$sizeselex%>%
  filter(Factor=='Lsel' & Yr=='2023')%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))

Report_estim.Len.only=SS_output("C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/New folder/1.Sel only",
                                covar=COVAR,forecast=FORECAST,readwt=F)
a_estim.Len.only=Report_estim.Len.only$sizeselex%>%
  filter(Factor=='Lsel' & Yr=='2023')%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))


Report_estim.Len.body.Zn1.index=SS_output("C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/New folder/3. Sel body Cpue Zn1",
                                          covar=COVAR,forecast=FORECAST,readwt=F)
a_estim.Len.body.Zn1.index=Report_estim.Len.body.Zn1.index$sizeselex%>%
  filter(Factor=='Lsel' & Yr=='2023')%>%
  dplyr::select(-c(Factor,Label))%>%
  gather(TL,Sel,-c( Fleet,Yr,Sex))%>%
  mutate(Sel=as.numeric(Sel),TL=as.numeric(TL),Fleet=as.character(Fleet),Sex=as.character(Sex))


rbind(Expnded.2022%>%mutate(Ass='2022'),
      a%>%mutate(Ass='S1'),
      a_estim.Len.only%>%mutate(Ass='Len.only'),
      a_estim.Len_Meanbody%>%mutate(Ass='Len & body'),
      a_estim.Len.body.Zn1.index%>%mutate(Ass='Len & body & Zn1 cpue'),
      a_fixedSel%>%mutate(Ass='fixed sel'))%>%
  filter(Fleet%in%3:8)%>%
  ggplot(aes(TL,Sel,color=Ass,shape=Sex))+
  geom_point()+geom_line()+
  facet_wrap(~Fleet,ncol=2)+
  theme_PA()+theme(legend.position = 'bottom')

#-----------  Check Spawning biomass (report)-----------------------------------------------------
rbind(Report_S1$timeseries%>%dplyr::select(Yr,SpawnBio)%>%mutate(SpawnBio=SpawnBio/SpawnBio[1],Ass='S1'),
      Report_S4$timeseries%>%dplyr::select(Yr,SpawnBio)%>%mutate(SpawnBio=SpawnBio/SpawnBio[1],Ass='S4'),
      Report2022$timeseries%>%dplyr::select(Yr,SpawnBio)%>%mutate(SpawnBio=SpawnBio/SpawnBio[1],Ass='2022'))%>%
  ggplot(aes(Yr,SpawnBio,color=Ass))+
  geom_point()+geom_line()+
  theme_PA()+geom_hline(yintercept = c(0.2,0.4))+
  xlim(range(Report_S1$catch$Yr))


#-----------  Check CPUE(report)-----------------------------------------------------
Expnded.S4.cpue=Report_S4$cpue%>%mutate(Fleet=case_when(Fleet==4~6,TRUE~Fleet))%>%dplyr::select(Fleet,Yr,Obs)
Expnded.S4.cpue=rbind(Expnded.S4.cpue,
                      Expnded.S4.cpue%>%filter(Fleet==3)%>%mutate(Fleet=4),
                      Expnded.S4.cpue%>%filter(Fleet==3)%>%mutate(Fleet=5),
                      Expnded.S4.cpue%>%filter(Fleet==6)%>%mutate(Fleet=7),
                      Expnded.S4.cpue%>%filter(Fleet==6)%>%mutate(Fleet=8))
Expnded.2022.cpue=Report2022$cpue%>%mutate(Fleet=case_when(Fleet==4~6,TRUE~Fleet))%>%dplyr::select(Fleet,Yr,Obs)
Expnded.2022.cpue=rbind(Expnded.2022.cpue,
                        Expnded.2022.cpue%>%filter(Fleet==3)%>%mutate(Fleet=4),
                        Expnded.2022.cpue%>%filter(Fleet==3)%>%mutate(Fleet=5),
                        Expnded.2022.cpue%>%filter(Fleet==6)%>%mutate(Fleet=7),
                        Expnded.2022.cpue%>%filter(Fleet==6)%>%mutate(Fleet=8))


rbind(Report_S1$cpue%>%dplyr::select(Fleet,Yr,Obs)%>%mutate(Ass='S1'),
      Expnded.S4.cpue%>%mutate(Ass='S4'),
      Expnded.2022.cpue%>%mutate(Ass='2022'))%>%
  ggplot(aes(Yr,Obs,color=Ass))+
  geom_point()+geom_line()+
  theme_PA()+facet_wrap(~Fleet,scales='free')

#-----------  Get probabilities and Risk-----------------------------------------------------
handl_desktop=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries And Regional Development/Desktop/test scenarios',x,sep='/')
Model.list.location=list(
  S1=handl_desktop('tunning_Whiskery/S1'),
  S1_zone1.only=handl_desktop('tunning_Whiskery/S1_daily cpue only Zn1'),
  S1_zone1.weight2=handl_desktop('tunning_Whiskery/S1_Zone1_w2'),
  S1_old.tagging=handl_desktop('tunning_Whiskery/S1_tagging'),
  S1_tagging_new_old=handl_desktop('tunning_Whiskery/S1_tagging new'),
  S4=handl_desktop('tunning_Whiskery/S4'),
  '2022'=handl_OneDrive('Analyses/Population dynamics/1.Whiskery shark/2022/SS3 integrated/S1')
)
store.prob.like.risk=Model.list.location
fn.get.prob.risk=function(report.location, SP, ASS)
{
  #get report file
  Report=SS_output(report.location,covar=TRUE,forecast=TRUE,readwt=F)
  
  #calculate probabilities
  dummy=fn.integrated.mod.get.timeseries(d=Report,mods="SS3",Type='Depletion',scen="S1") 
  
  #calculate cons and likes
  Cons.Like=rbind(dummy$Probs$probs,dummy$Probs$probs.future)%>%
    mutate(Species=SP,
           Range=factor(Range,levels=c('<lim','lim.thr','thr.tar','>tar')))%>%
    dplyr::select(Range,finyear,Species,Probability)%>%
    arrange(Range)
  #Calculate risk
  Risk=fn.risk(d=Cons.Like,w=1)%>%
    mutate(Ass=ASS,
           when=ifelse(finyear==min(finyear),'now','future'))
  Risk_now=fn.risk.figure(d=Risk%>%filter(finyear==min(finyear)), Risk.colors=RiskColors, out.plot=FALSE)%>%
    dplyr::select(Species,Consequence,Likelihood,Risk,Score)%>%
    mutate(Ass=ASS,when='now')
  Risk_future=fn.risk.figure(d=Risk%>%filter(finyear==max(finyear)), Risk.colors=RiskColors, out.plot=FALSE)%>%
    dplyr::select(Species,Consequence,Likelihood,Risk,Score)%>%
    mutate(Ass=ASS,when='future')
  
  return(list(Risk=Risk,Overall.risk=rbind(Risk_now,Risk_future)))
}
for(q in 1:length(Model.list.location))
{
  store.prob.like.risk[[q]]=fn.get.prob.risk(report.location=Model.list.location[[q]],
                                             SP='Whiskery shark',
                                             ASS=names(Model.list.location)[q])
}

do.call(rbind,lapply(store.prob.like.risk, function(x) x$Risk))%>%
  mutate(when=factor(when,levels=c('now','future')))%>%
  ggplot(aes(Consequence,Probability,color=Ass))+
  geom_point(size=4)+facet_wrap(~when,ncol=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

do.call(rbind,lapply(store.prob.like.risk, function(x) x$Overall.risk))%>%
  mutate(when=factor(when,levels=c('now','future')))%>%
  ggplot(aes( Ass,Score,color=Risk))+
  geom_point(size=6)+facet_wrap(~when,ncol=1)+
  scale_color_manual(values = RiskColors)+
  geom_text_repel(aes(label=paste0('C=',Consequence, ', L=',Likelihood)),color='black')+
  ylim(0,NA)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



