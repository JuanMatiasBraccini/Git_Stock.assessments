dummy.store=vector('list',N.sp)
names(dummy.store)=Keep.species
dummy.store.estimates=dummy.store.likelihoods=dummy.store.quantities=dummy.store.rel.biom=
  dummy.store.probs.rel.biom=dummy.store.probs.f.series=dummy.store.probs.B.Bmsy=dummy.store.f.series=
  dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.Kobe.probs=
  dummy.store.sens.table=dummy.store.ensemble=dummy.store

for(i in 1:length(dummy.store))
{
  Neim=names(dummy.store)[i]
  
  if((!is.null(Catch.rate.series[[i]]) | Neim%in%Species.with.length.comp) & !(Neim%in%no.empirical.sel.main.fleet))
  {
    this.wd=paste(HandL.out,capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
    if(!dir.exists(this.wd))dir.create(this.wd)
    
    Life.history=List.sp[[i]]
    
    #1. Catch
    #1.1. zones together
    ktch=KtCh%>%
      filter(Name==Neim)%>%
      mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                           ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                    'TEP_greynurse','TEP_dusky','Discards_TDGDLF'),'Southern.shark',
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
    ktch.zone=ktch.zone%>%
      mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                           ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                    'TEP_greynurse','TEP_dusky','Discards_TDGDLF'),'Southern.shark',
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
    if(any(grepl('Size_composition',names(Species.data[[i]]))))
    {
      d.list.n.shots=Species.data[[i]][grep(paste(c("Size_composition_Survey_Observations","Size_composition_Observations",
                                                    "Size_composition_Other_Observations"),collapse="|"),
                                            names(Species.data[[i]]))]
      d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                    names(Species.data[[i]]))]
      
      if(length(d.list)>0)
      {
        if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
        if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
        
        for(s in 1:length(d.list))
        {
          d.list[[s]]=d.list[[s]]%>%
            filter(FL>=Life.history$Lzero*1.05)%>%
            mutate(fishry=ifelse(grepl("NSF.LONGLINE",names(d.list)[s]),'NSF',
                                 ifelse(grepl("Survey",names(d.list)[s]),'Survey',
                                        ifelse(grepl("Other",names(d.list)[s]),'Other',
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
          
          #add extra bins for smooth fit to size comps
          Maximum_size=10*ceiling(with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/10)
          Mx.size=max(d.list[[s]]$size.class)
          extra.bins=seq(Mx.size+TL.bins.cm,10*round(Maximum_size/10),by=TL.bins.cm) 
          if(length(extra.bins))
          {
            add.dumi.size=d.list[[s]][1:length(extra.bins),]%>%
              mutate(size.class=extra.bins,
                     n=0)
            d.list[[s]]=rbind(d.list[[s]],add.dumi.size)
          }
        }
        d.list <- d.list[!is.na(d.list)]
        d.list=do.call(rbind,d.list)   
        if(Neim%in%combine_NSF_Survey) #assume same sel NSF and survey so combined to increase sample size
        {
          d.list=d.list%>%
            mutate(fishry=ifelse(fishry=='NSF',"Survey",fishry))
        }
        if(Neim%in%combine.sexes)
        {
          if(Neim%in%combine.sexes.survey)
          {
            d.list$sex=ifelse(d.list$fishry=="Survey",0,d.list$sex)
          }
          if(Neim%in%combine.sexes.tdgdlf)
          {
            d.list=d.list%>%mutate(sex=ifelse(fishry=="TDGDLF",0,sex))
          }
          if(Neim%in%combine.sexes.tdgdlf.daily)
          {
            d.list=d.list%>%mutate(sex=ifelse(fishry=="TDGDLF" & year>2005,0,sex))
          }
          if(!Neim%in%c(combine.sexes.survey,combine.sexes.tdgdlf,combine.sexes.tdgdlf.daily))
          {
            d.list$sex=0 
          }
        }
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
        
        if(Neim%in%names(Indicator.species))
        {
          Min.size=Min.annual.obs.ktch
        }else
        {
          Min.size=Min.annual.obs.ktch*prop.min.N.accepted_other
        }
        Table.n=d.list%>%group_by(year,fishry,sex)%>%
          summarise(N=sum(n))%>%
          mutate(Min.accepted.N=ifelse(!fishry=='Survey',Min.size,Min.annual.obs.ktch_survey))%>%
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
          
          d.list.0=d.list%>%filter(Sex==0)%>%arrange(year)
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
          
          min.nsamp=Min.Nsamp
          if(!Neim%in%Indicator.species) min.nsamp=min.nsamp/2
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
            filter(year<=max(ktch$finyear) & Nsamp>=min.nsamp)  
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
      if(length(d.list)>0)
      {
        if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
        if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
        
        for(s in 1:length(d.list))
        {
          NM=names(d.list)[s]
          d.list[[s]]=d.list[[s]]%>%
            filter(FL>=Life.history$Lzero*1.05)%>%
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
          Maximum_size=10*ceiling(with(Life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/10)
          Mx.size=max(d.list[[s]]$size.class)
          extra.bins=seq(Mx.size+TL.bins.cm,10*round(Maximum_size/10),by=TL.bins.cm) 
          if(length(extra.bins))
          {
            add.dumi.size=d.list[[s]][1:length(extra.bins),]%>%
              mutate(size.class=extra.bins,
                     n=0)
            d.list[[s]]=rbind(d.list[[s]],add.dumi.size)
          }
        }
        d.list <- d.list[!is.na(d.list)]
        d.list=do.call(rbind,d.list)   
        if(Neim%in%combine_NSF_Survey) #assume same sel NSF and survey so combined to increase sample size
        {
          d.list=d.list%>%
            mutate(fishry=ifelse(fishry=='NSF',"Survey",fishry))
        }
        if(Neim%in%combine.sexes)
        {
          if(Neim%in%combine.sexes.survey)
          {
            d.list$sex=ifelse(d.list$fishry=="Survey",0,d.list$sex)
          }
          if(Neim%in%combine.sexes.tdgdlf)
          {
            d.list=d.list%>%mutate(sex=ifelse(grepl("TDGDLF",fishry),0,sex))
          }
          if(Neim%in%combine.sexes.tdgdlf.daily)
          {
            d.list=d.list%>%mutate(sex=ifelse(grepl("TDGDLF",fishry) & year>2005,0,sex))
          }
          if(!Neim%in%c(combine.sexes.survey,combine.sexes.tdgdlf,combine.sexes.tdgdlf.daily))
          {
            d.list$sex=0 
          }
        }
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
        if(Neim%in%names(Indicator.species))
        {
          Min.size=Min.annual.obs.ktch
        }else
        {
          Min.size=Min.annual.obs.ktch*prop.min.N.accepted_other
        }
        Table.n=d.list%>%group_by(year,fishry,sex)%>%
          summarise(N=sum(n))%>%
          mutate(Min.accepted.N=ifelse(!fishry=='Survey',Min.size,Min.annual.obs.ktch_survey))%>%
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
          
          d.list.0=d.list%>%filter(Sex==0)%>%arrange(year)
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
          
          min.nsamp=Min.Nsamp
          if(!Neim%in%Indicator.species) min.nsamp=min.nsamp/2
          dummy.Size.compo.SS.format=dummy.Size.compo.SS.format%>%
            filter(year<=max(ktch$finyear) & Nsamp>=min.nsamp)  
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
    
    #3. meanbodywt
    #zones together
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
      for(y in 1:length(xxx))
      {
        ZnE=capitalize(str_remove(names(xxx)[y],"annual.mean.size_"))
        xxx[[y]]=xxx[[y]]%>%
          mutate(year=as.numeric(substr(Finyear,1,4)),
                 month=1,
                 Fleet=paste('Southern.shark_2',ZnE,sep='_'),
                 part=0,   #0, combined; 1: discard only; 2: retained only
                 type=2)%>%
          filter(year<=max(ktch.zone$finyear))%>%
          dplyr::select(-Finyear)
      }
      
      
      meanbodywt.SS.format.zone=do.call(rbind,xxx)%>%
        left_join(Flits.and.survey.zone,by=c('Fleet'='Fleet.name'))%>%
        mutate(fleet=Fleet.number)%>%
        dplyr::select(-c(Fleet.number,Fleet))%>%
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
    
    
    #4. Abundance series
    Abundance.SS.format=NULL
    Var.ad.factr=NULL
    CPUE=compact(Catch.rate.series[[i]])
    if(!is.null(CPUE))
    {
      if(Neim%in%survey_not.representative & any(grepl("Survey",names(CPUE)))) CPUE=CPUE[-grep("Survey",names(CPUE))]
      if(Neim%in%NSF_not.representative & any(grepl("NSF",names(CPUE)))) CPUE=CPUE[-grep("NSF",names(CPUE))]
      if(Neim%in%tdgdlf_not.representative & any(grepl("TDGDLF",names(CPUE)))) CPUE=CPUE[-grep("TDGDLF",names(CPUE))]
      if(Neim%in%tdgdlf_monthly_not.representative & "TDGDLF.monthly"%in%names(CPUE)) CPUE=CPUE[-grep("TDGDLF.monthly",names(CPUE))]
      if(!is.null(Life.history$drop.monthly.cpue)) CPUE$TDGDLF.monthly=CPUE$TDGDLF.monthly%>%filter(!yr.f%in%Life.history$drop.monthly.cpue)
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
    
    
    #7. Conditional age at length
    Cond.age.len.SS.format=NULL
    if(do.Cond.age.len.SS.format)
    {
      if(any(grepl('age_length',names(Species.data[[i]]))))
      {
        a=Life.history$a_FL.to.TL
        b=Life.history$b_FL.to.TL
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
    
    
    #8. MeanSize at Age obs
    MeanSize.at.Age.obs.SS.format=NULL
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
    
    
    #9. Fleet info
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
    
    
    #10. Run scenarios if available abundance index or size comps 
    len.cpue=length(CPUE)
    len.size.comp=length(Size.compo.SS.format)
    if(len.cpue>0|len.size.comp>0)
    {
      #Put CPUE in SS format
      if(len.cpue>0)
      {
        MAX.CV=Life.history$MAX.CV
        
        #zones combined
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
        
        #by zone
        len.cpue.zone=length(CPUE.zone)
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
      
      #Scenarios
      if(!do.all.sensitivity.tests) Life.history$Sens.test$SS=Life.history$Sens.test$SS%>%filter(Scenario=='S1')
      Scens=Life.history$Sens.test$SS%>%
        mutate(Species=capitalize(Neim))
      Store.sens=vector('list',nrow(Scens))
      names(Store.sens)=Scens$Scenario
      Out.Scens=Scens
      Out.estimates=Out.likelihoods=Out.quantities=Out.rel.biom=Out.probs.rel.biom=Out.f.series=Out.B.Bmsy=
        Out.probs.f.series=Out.F.Fmsy=Out.probs.B.Bmsy=Out.Kobe.probs=store.warnings=store.convergence=vector('list',length(Store.sens))
      
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
      
      #rename fleets following SS nomenclature
      names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))]=match(names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))],names(Flits))
      names(ktch.zone)[which(!names(ktch.zone)%in%c("SPECIES","Name","finyear"))]=match(names(ktch.zone)[which(!names(ktch.zone)%in%c("SPECIES","Name","finyear"))],names(Flits.zone))
      
      #future catches
      if("SS"%in%future.models)
      {
        #zones combined
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
        #remove daily years  
        if(!is.na(Scens$Daily.cpues[s]) & 'TDGDLF.daily'%in%names(CPUE))
        {
          rid.of=as.numeric(unlist(str_split(Scens$Daily.cpues[s], "&")))
          #zones combined
          drop.dis=Abundance.SS.format%>%
            mutate(this=grepl('TDGDLF.daily',rownames(Abundance.SS.format)) & Year%in%rid.of)
          if(any(drop.dis$this)) Abundance.SS.format=Abundance.SS.format[-which(drop.dis$this),]
          #by zone
          drop.dis=Abundance.SS.format.zone%>%
            mutate(this=grepl('TDGDLF.daily',rownames(Abundance.SS.format.zone)) & Year%in%rid.of)
          if(any(drop.dis$this)) Abundance.SS.format.zone=Abundance.SS.format.zone[-which(drop.dis$this),]
        }
        
        #use cpue or length comp in likelihood?
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
        
        #a. Create input files
        if(create.SS.inputs)
        {
          #Specify ktch,flitinfo,abundance,comps,mean.weight,var.adjs & future based on scenario   
          if(is.na(Scens$Spatial[s]))
          {
            KAtch=ktch
            FLitinFO=flitinfo
            Abund=Abundance.SS.format
            Size.com=Size.compo.SS.format
            meanbody=meanbodywt.SS.format
            Var.ad=Var.ad.factr
            add.future=add.ct.or.F_future
          }
          if(Scens$Spatial[s]=='areas-as-fleets')
          {
            KAtch=ktch.zone
            FLitinFO=flitinfo.zone
            Abund=Abundance.SS.format.zone
            Size.com=Size.compo.SS.format.zone
            meanbody=meanbodywt.SS.format.zone
            Var.ad=Var.ad.factr.zone
            add.future=add.ct.or.F_future.zone
          }
          
          fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SS'),   
                       new.path=this.wd1,
                       Scenario=Scens[s,]%>%
                         mutate(Model='SS'),
                       Catch=KAtch,  
                       life.history=Life.history,
                       depletion.yr=NULL,
                       fleets=names(KAtch)[which(!names(KAtch)%in%c("SPECIES","Name","finyear"))],
                       fleetinfo=FLitinFO,
                       abundance=Abund,   
                       size.comp=Size.com,
                       meanbodywt=meanbody,
                       F.tagging=F.SS.format,
                       cond.age.len=Cond.age.len.SS.format,
                       MeanSize.at.Age.obs=MeanSize.at.Age.obs.SS.format,
                       Lamdas=Lamdas.SS.lambdas,
                       Var.adjust.factor=Var.ad,
                       Future.project=add.future) 
        }
        
        #b. Run SS3
        #run this first time fitting model to define LnRo init value
        if(Find_Init_LnRo)
        {
          start <- r4ss::SS_readstarter(file = file.path(this.wd1, "starter.ss"), verbose = FALSE)
          start$last_estimation_phase=0
          r4ss::SS_writestarter(start, dir = this.wd1, overwrite = TRUE,verbose = FALSE)
          fn.run.SS(where.inputs=this.wd1,
                    where.exe=handl_OneDrive('SS3/ss_win.exe'),
                    args="-nohess") 
          Report=SS_output(this.wd1,covar=F,forecast=F,readwt=F)  
          Report$timeseries%>%filter(Era=='VIRG')%>%pull(Bio_all) #JABBA K= 6800 tonnes
          rm(Report)
        }
        #run this to tune model and calculate RAMP years
        if(Scens$Scenario[s]=='S1' & Calculate.ramp.years)
        {
          #tune ramp years
          fn.run.SS(where.inputs=this.wd1,
                    where.exe=handl_OneDrive('SS3/ss_win.exe'),
                    args='')
          Report=SS_output(this.wd1)
          tiff(file=paste(this.wd,'Ramp_years.tiff',sep='/'),
               width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
          ramp_years=SS_fitbiasramp(Report) 
          dev.off()
          out=ramp_years$df
          out=rbind(out,data.frame(value=unique(Report$sigma_R_info$alternative_sigma_R),label='Alternative_sigma_R'))
          write.csv(out,paste(this.wd,'Ramp_years.csv',sep='/'),row.names = F)
          
          these.plots=c(1:7,10,11,16,26)  #biol, selectivity, timeseries,rec devs,S-R,catch,mean weight, indices, size comp
          SS_plots(Report, plot=these.plots, png=T)
          
          #tune composition data
          tune_info <- tune_comps(option = "Francis",
                                  niters_tuning = 1,
                                  dir = this.wd1,
                                  exe=handl_OneDrive('SS3/ss_win.exe'),
                                  allow_up_tuning = TRUE,
                                  verbose = FALSE)
          write.csv(tune_info$weights[[1]]%>%mutate(Method='Francis'),paste(this.wd,'Tuned_size_comp.csv',sep='/'),row.names = F)
          
          rm(ramp_years,out,tune_info)
        }
        #run this to estimate parameters
        if(Run.SS)
        {
          fn.run.SS(where.inputs=this.wd1,
                    where.exe=handl_OneDrive('SS3/ss_win.exe'),
                    args=Arg) 
        }
        
        
        #c. Bring in outputs
        COVAR=FORECAST=FALSE
        if(Arg=="") COVAR=TRUE
        if("SS"%in%future.models) FORECAST=TRUE
        Report=SS_output(this.wd1,covar=COVAR,forecast=FORECAST,readwt=F)
        
        if(run_SS_plots) SS_plots(Report,  png=T)
        
        #d. Store estimates
        Estims=Report[["estimated_non_dev_parameters"]]
        Out.estimates[[s]]=Estims%>%
          mutate(Par=rownames(Estims),
                 Scenario=Scens$Scenario[s])%>%
          relocate(Scenario,Par)%>%
          `rownames<-`( NULL )
        F.ref.points=fn.get.f.ref.points(Report)                 
        dummi.F=Out.estimates[[s]][1:length(F.ref.points),]
        dummi.F[,]=NA
        dummi.F=dummi.F%>%
          mutate(Scenario=unique(Out.estimates[[s]]$Scenario),
                 Par=names(F.ref.points),
                 Value=unlist(F.ref.points))
        Out.estimates[[s]]=rbind(Out.estimates[[s]],dummi.F)
        Likelihoods=Report[["likelihoods_used"]]
        Out.likelihoods[[s]]=Likelihoods%>%
          mutate(Likelihood=rownames(Likelihoods),
                 Scenario=Scens$Scenario[s])%>%
          relocate(Scenario,Likelihood)%>%
          `rownames<-`( NULL )
        Bratio_current=paste0('Bratio_',max(ktch$finyear))
        Out.quantities[[s]]=Report[["derived_quants"]]%>%
          filter(Label%in%c(Bratio_current,"Dead_Catch_MSY"))%>%
          mutate(Label=ifelse(Label==Bratio_current,'Current depletion',
                              ifelse(Label=="Dead_Catch_MSY",'MSY',
                                     Label)),
                 Species=capitalize(Neim),
                 Scenario=Scens$Scenario[s])%>%
          rename(Median=Value,
                 SE=StdDev)%>%
          dplyr::select(Species,Label,Median,SE,Scenario)%>%
          `rownames<-`( NULL )
        
        #e. Store trajectories
        
        #e.1 relative biomass
        dummy=fn.integrated.mod.get.timeseries(d=Report,
                                               mods="SS3",
                                               Type='Depletion',
                                               scen=Scens$Scenario[s])
        Out.rel.biom[[s]]=dummy$Dat
        Out.probs.rel.biom[[s]]=dummy$Probs   #based on asymptotic error
        
        #e.2 F
        dummy=fn.integrated.mod.get.timeseries(d=Report,
                                               mods="SS3",
                                               Type='F.series',
                                               scen=Scens$Scenario[s])
        Out.f.series[[s]]=dummy$Dat
        Out.probs.f.series[[s]]=dummy$Probs
        
        #e.3 B/Bmsy
        dummy=fn.integrated.mod.get.timeseries(d=Report,
                                               mods="SS3",
                                               Type='B.Bmsy',
                                               scen=Scens$Scenario[s],
                                               get.uncertainty='CV')
        Out.B.Bmsy[[s]]=dummy$Dat
        Out.probs.B.Bmsy[[s]]=dummy$Probs
        Kobe.stock=dummy$Out.Kobe
        
        #e.4 F/Fmsy
        dummy=fn.integrated.mod.get.timeseries(d=Report,
                                               mods="SS3",
                                               Type='F.Fmsy',
                                               scen=Scens$Scenario[s],
                                               get.uncertainty='ratio.first')
        Out.F.Fmsy[[s]]=dummy$Dat
        Kobe.harvest=dummy$Out.Kobe
        
        #Calculate posterior depletion   #takes 4.5 secs per nMC simulation
        if(Scens$Scenario[s]=='S1' & SS3.run=='final' & do.MC.multi)
        {
          dd=fn.MC.sims(this.wd1,
                        nMC=nMCsims,
                        arg=Arg.no.estimation,
                        B.th=Report$derived_quants%>%filter(Label=="SSB_MSY")%>%pull(Value)/Report$derived_quants%>%filter(Label=="SSB_Virgin")%>%pull(Value),
                        scen=Scens$Scenario[s])
          
          Out.rel.biom[[s]]=dd$Depletion
          Out.probs.rel.biom[[s]]=dd[c("probs","Reference.points")]
          Out.f.series[[s]]=dd$F.series
          Out.B.Bmsy[[s]]=dd$B.Bmsy
          Out.F.Fmsy[[s]]=dd$F.Fmsy
          
          rm(dd)
        }
        
        #add catch for display purposes
        Add.katch=Combined.ktch
        if(any(add.ct.future$finyear%in%unique(Out.rel.biom[[s]]$year)))
        {
          TC.future=add.ct.future%>%ungroup()%>%dplyr::select(-c(SPECIES,Name,finyear))%>%rowSums()
          Add.katch=rbind(Add.katch,data.frame(Year=add.ct.future$finyear,Total=TC.future))
        }
        Out.rel.biom[[s]]=Out.rel.biom[[s]]%>%
          left_join(Add.katch,by=c('year'='Year'))%>%
          rename(Catch=Total)
        
        #Kobe
        if(Scens$Scenario[s]=='S1') Out.Kobe.probs[[s]]=data.frame(stock=Kobe.stock,
                                                                   harvest=Kobe.harvest)%>%
          filter(stock>=0 & harvest>=0)
        
        
        print(paste("___________","SS3 Scenario",Scens$Scenario[s],"___________",Neim))
        
        
        # Evaluate fit diagnostics
        GoodnessFit=function.goodness.fit_SS(Rep=Report)  
        write.csv(GoodnessFit,paste(this.wd1,"/GoodnessFit_",Neim,".csv",sep=''))
        
      } #end s loop
      Out.Scens=Out.Scens%>%
        rename(h.mean=Steepness,
               h.sd=Steepness.sd)%>%
        mutate(Mmean=round(Mmean,3),
               h.mean=round(h.mean,3),
               h.sd=round(h.sd,3))%>%
        dplyr::select(-c(Model,use_F_ballpark))%>%
        relocate(Species,Scenario,Mmean,F_ballpark,h.mean,h.sd)
      
      #add extra scenario info
      fec.alpha=Life.history$Fecu_a    #intercept
      fec.beta=Life.history$Fecu_b     #slope
      if(is.na(fec.alpha)| is.na(fec.beta))
      {
        fec.alpha=mean(Life.history$Fecundity)  #set fec-@-age to mean fec if no relationship available
        fec.beta=0    
      }
      Out.Scens=Out.Scens%>%
        mutate(TLo=round(with(Life.history,Lzero*a_FL.to.TL+b_FL.to.TL)),
               Fem.TL.inf=round(with(Life.history,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)),
               Fem.K=Life.history$Growth.F$k,
               Fem.awt=Life.history$AwT,
               Fem.bwt=Life.history$BwT,
               Fem.TL50=round(Life.history$TL.mat.inf.slope[2]),
               Fem.fec.a=fec.alpha/mean(Life.history$Breed.cycle),
               Fem.fec.b=fec.beta/mean(Life.history$Breed.cycle),
               Mal.TL.inf=round(with(Life.history,Growth.M$FL_inf*a_FL.to.TL+b_FL.to.TL)),
               Mal.K=Life.history$Growth.M$k,
               Mal.awt=Life.history$AwT.M,
               Mal.bwt=Life.history$BwT.M,
               a_FL.to.TL=Life.history$a_FL.to.TL,
               b_FL.to.TL=Life.history$b_FL.to.TL,
               Fleets=paste(flitinfo$fleetname,collapse=', '))
      
      
      #5. Store quantities
      dummy.store.sens.table[[i]]=Out.Scens
      dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
      dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
      dummy.store.probs.B.Bmsy[[i]]=Out.probs.B.Bmsy
      dummy.store.probs.f.series[[i]]=Out.probs.f.series
      dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
      dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
      dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
      dummy.store.Kobe.probs[[i]]=Out.Kobe.probs  
      dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
      dummy.store.likelihoods[[i]]=do.call(rbind,Out.likelihoods)
      dummy.store.quantities[[i]]=do.call(rbind,Out.quantities)
      
      rm(Out.Scens,Out.rel.biom,Out.probs.rel.biom,Out.f.series,
         Out.B.Bmsy,Out.F.Fmsy,Out.estimates,Out.Kobe.probs,Out.likelihoods,Out.quantities)
      
      #plot overall figure with estimates and scenarios
      yrS=dummy.store.rel.biom[[i]]%>%filter(Scenario=='S1')%>%pull(year)
      if(length(yrS)<50) delta=10 else
        delta=17
      xmin=min(yrS)+delta
      xmax=xmin+delta
      p=dummy.store.rel.biom[[i]]%>%
        ggplot(aes(year,median,color=Scenario))+
        geom_line(size=2)+
        geom_line(aes(year,upper.95),linetype=2)+
        geom_line(aes(year,lower.95),linetype=2)+
        ggtitle(Keep.species[i])+ylim(0,1)+
        theme_PA()+theme(legend.position = 'bottom')+
        annotation_custom(tableGrob(dummy.store.sens.table[[i]]%>%dplyr::select(Scenario,Mmean,h.mean)),
                          xmin=xmin+5, xmax=xmax+5, ymin=0, ymax=0.3)+
        annotation_custom(tableGrob(dummy.store.estimates[[i]]%>%filter(Par=='SR_LN(R0)')),xmin=xmin, xmax=xmax, ymin=0.35, ymax=0.55)
      print(p)
      ggsave(paste(this.wd,"/Rel.biomass&estimates.tiff",sep=''),compression = "lzw")
    }
    
    clear.log("Var.ad.factr")
    clear.log("Var.ad.factr.zone")
  }
} #end i loop

Age.based[[w]]$sens.table=dummy.store.sens.table
Age.based[[w]]$estimates=dummy.store.estimates
Age.based[[w]]$rel.biom=dummy.store.rel.biom
Age.based[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
Age.based[[w]]$probs.B.Bmsy=dummy.store.probs.B.Bmsy
Age.based[[w]]$probs.f.series=dummy.store.probs.f.series
Age.based[[w]]$f.series=dummy.store.f.series
Age.based[[w]]$B.Bmsy=dummy.store.B.Bmsy
Age.based[[w]]$F.Fmsy=dummy.store.F.Fmsy
Age.based[[w]]$Kobe.probs=dummy.store.Kobe.probs
Age.based[[w]]$likelihoods=dummy.store.likelihoods
Age.based[[w]]$quantities=dummy.store.quantities


rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,dummy.store.probs.f.series,
   dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,dummy.store.probs.B.Bmsy,
   dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.Kobe.probs,dummy.store.likelihoods,
   dummy.store.quantities)
