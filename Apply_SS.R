#---Create input files -------------------------------------------------
n.SS=length(Integrated.age.based)  
Age.based=vector('list',n.SS)
names(Age.based)=Integrated.age.based
if(create.SS.inputs)
{
  for(w in 1:n.SS)
  {
    # SS3
    if(names(Age.based)[w]=="SS") 
    {
      for(i in 1:N.sp)
      {
        Neim=Keep.species[i]
        
        if((!is.null(Catch.rate.series[[i]]) | Neim%in%Species.with.length.comp) & !(Neim%in%no.empirical.sel.main.fleet))
        {
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
              if(retained.discarded.units=='numbers' & Neim=="dusky shark")  
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
                if(Maximum_size>Mx.size)
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
              #combine sexes if number of obs per year >Min.size but per sex <Min.size/2
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
                mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey/2,
                                                fishry=='NSF'~Min.size.NSF/2,
                                                TRUE~Min.size/2))
              Len.comps.sex.ratios=Table.n.sex
              Len.comps.below.Min.accepted.N=Len.comps.sex.ratios%>%
                filter(N<Min.accepted.N)
              Len.comps.sex.ratios=Len.comps.sex.ratios%>%
                dplyr::select(-Min.accepted.N)%>%
                mutate(sex=ifelse(sex=='F','Female','Male'))%>%
                spread(sex,N,fill=0)%>%
                mutate(Total=Female+Male,
                       Fem.prop=Female/Total)
              Table.n.sex=Table.n.sex%>%   
                filter(N<Min.accepted.N)%>%
                mutate(dummy=paste(year,fishry,sex))%>%
                distinct(dummy)
              d.list=d.list%>%  
                mutate(dummy=paste(year,fishry,sex))%>%
                filter(!dummy%in% with(Len.comps.below.Min.accepted.N,paste(year,fishry,sex)))%>%
                mutate(sex=ifelse(dummy%in%Table.n.sex$dummy,combine.sex_type,sex))%>%
                dplyr::select(-dummy)
              if(Neim%in%combine.sexes.tdgdlf.single.area) d.list$sex=combine.sex_type
              Apical.prop.female=Len.comps.sex.ratios%>%
                filter(year<2015)%>%
                group_by(fishry) %>%
                summarise(w_mean = weighted.mean(Fem.prop, Total),
                          w_sd = sqrt(wtd.var(Fem.prop, Total, na.rm = TRUE)))
              d.list=d.list%>%
                group_by(year,fishry,sex,size.class)%>%
                summarise(n=sum(n))%>%
                ungroup()%>%
                filter(!is.na(year))%>%
                filter(!is.na(fishry))%>%
                data.frame()
              
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
                
                if(nrow(d.list.0)>0)  
                {
                  dummy.zeros.0=d.list.0[,-match(vars.head,names(d.list.0))]
                  dummy.zeros.0[,]=0
                  names(d.list.0)[id.var.nms.f]=paste('f',names(d.list.0)[id.var.nms.f],sep='')  
                  names(dummy.zeros.0)=paste('m',names(dummy.zeros.0),sep='')
                  dummy.Size.compo.SS.format_Sex0=bind_cols(d.list.0,dummy.zeros.0)%>%
                    mutate(Sex=as.numeric(Sex))
                  
                }
                if(nrow(d.list.f)>0)
                {
                  dummy.zeros.f=d.list.f[,-match(vars.head,names(d.list.f))]
                  dummy.zeros.f[,]=0
                  dummy.zeros.m=d.list.m[,-match(vars.head,names(d.list.m))]
                  dummy.zeros.m[,]=0
                  
                  Fem.Da=d.list.f%>%data.frame(check.names = FALSE)       
                  names(Fem.Da)[id.var.nms.f]=paste('f',names(Fem.Da)[id.var.nms.f],sep='')
                  Fem.Da.0=dummy.zeros.f%>%data.frame(check.names = FALSE)
                  names(Fem.Da.0)=paste('m',names(Fem.Da.0),sep='')
                  Fem.Da=bind_cols(Fem.Da,Fem.Da.0)
                  
                  Mal.hdr=d.list.m[,match(vars.head,names(d.list.m))]%>%data.frame(check.names = FALSE)
                  Mal.Da.0=dummy.zeros.m%>%data.frame(check.names = FALSE)
                  names(Mal.Da.0)=paste('f',names(Mal.Da.0),sep='')
                  Mal.Da=d.list.m[,-match(vars.head,names(d.list.m))]%>%data.frame(check.names = FALSE)
                  names(Mal.Da)=paste('m',names(Mal.Da),sep='')
                  Mal.Da=bind_cols(Mal.hdr,Mal.Da.0,Mal.Da)
                  
                  dummy.Size.compo.SS.format_Sex=bind_rows(Fem.Da,Mal.Da)%>%
                    mutate(Sex=as.numeric(Sex))
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
                  dummy.Size.compo.SS.format=bind_rows(dummy.Size.compo.SS.format_Sex0,dummy.Size.compo.SS.format_Sex)
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
            d.list=Species.data[[i]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),names(Species.data[[i]]))]
            
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
              
              #explore length comps
              if(First.run=='YES')
              {
                p=do.call(rbind,d.list)%>%
                  rownames_to_column(var = "fishry")%>%
                  filter(!is.na(SEX))%>%
                  mutate(fishry=str_remove_all(fishry,"^Size_composition_|\\.inch.*"),
                         fishry=sub("\\..*", "", fishry))%>%
                  mutate(fleet=ifelse(year<2006,paste('1.Monthly',fishry),paste('2.Daily',fishry)),
                         year=factor(year),
                         TL=FL*Life.history$a_FL.to.TL+Life.history$b_FL.to.TL)
                p=p%>%
                  ggplot(aes(TL,fill=year))+
                  geom_density(alpha=.5)+
                  facet_grid(SEX~fleet)+theme_PA()+
                  scale_fill_manual(values = colorRampPalette(c("red",'orange', "steelblue"))(length(unique(p$year))))
                print(p)
                ggsave(paste(DiR,'Length comps considered for SS_areas as fleets.tiff',sep='/'), 
                       width = 10,height = 5, dpi = 300, compression = "lzw")
                
              }
              #calculate TL
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
                  ungroup()%>%
                  data.frame()
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
                if(Maximum_size>Mx.size)
                {
                  extra.bins=seq((Mx.size+TL.bins.cm),Maximum_size,by=TL.bins.cm)
                  #extra.bins=seq(Mx.size+TL.bins.cm,10*round(Maximum_size/10),by=TL.bins.cm) 
                }
                
                if(any(!is.na(extra.bins)))
                {
                  add.dumi.size=d.list[[s]][1:length(extra.bins),]%>%
                    mutate(size.class=extra.bins,
                           n=0)
                  d.list[[s]]=bind_rows(d.list[[s]],add.dumi.size)
                }
              }
              d.list <- d.list[!is.na(d.list)]
              # Display sex ratio by zone used in SS  
              if(First.run=="YES")
              {
                add.n.samps=Species.data[[i]]$Size_composition_Observations%>%
                  filter(Method=='GN')%>%
                  mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
                  rename(Zone=zone)%>%
                  group_by(Zone,year)%>%
                  summarise(N.shots=sum(N.shots))%>%ungroup()
                fn.ktch.sex.ratio.zone_SS(size.data=d.list,Min.size=Min.size,N_sampleS=add.n.samps)
                ggsave(paste(DiR,'Sex ratio by zone_SS size comps data_areas as fleets.tiff',sep='/'), width = 5,height = 6, dpi = 300, compression = "lzw")
              }
              # Get min sample size
              if(Neim%in%names(Indicator.species))
              {
                Min.size=Min.annual.obs.ktch.zone
              }else
              {
                Min.size=Min.annual.obs.ktch.zone*prop.min.N.accepted_other
              }
              Min.size.NSF=Min.annual.obs.ktch_NSF
              if(Neim%in%c("dusky shark")) Min.size.NSF=20
              
              
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
              #combine sexes if number of obs per year >Min.size but per sex <Min.size/2
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
                mutate(Min.accepted.N=case_when(fishry=='Survey'~Min.annual.obs.ktch_survey/2,
                                                fishry=='NSF'~Min.size.NSF/2,
                                                TRUE~Min.size/2))
              Len.comps.sex.ratios=Table.n.sex
              Len.comps.below.Min.accepted.N=Len.comps.sex.ratios%>%
                filter(N<Min.accepted.N)
              Len.comps.sex.ratios=Len.comps.sex.ratios%>%
                dplyr::select(-Min.accepted.N)%>%
                mutate(sex=ifelse(sex=='F','Female','Male'))%>%
                spread(sex,N,fill=0)%>%
                mutate(Total=Female+Male,
                       Fem.prop=Female/Total)
              Table.n.sex=Table.n.sex%>%   
                filter(N<Min.accepted.N)%>%
                mutate(dummy=paste(year,fishry,sex))%>%
                distinct(dummy)
              d.list=d.list%>%  
                mutate(dummy=paste(year,fishry,sex))%>%
                filter(!dummy%in% with(Len.comps.below.Min.accepted.N,paste(year,fishry,sex)))%>%
                mutate(sex=ifelse(dummy%in%Table.n.sex$dummy,combine.sex_type,sex))%>%
                dplyr::select(-dummy)
              Apical.prop.female.zone=Len.comps.sex.ratios%>%
                filter(year<2015)%>%
                group_by(fishry) %>%
                summarise(w_mean = weighted.mean(Fem.prop, Total),
                          w_sd = sqrt(wtd.var(Fem.prop, Total, na.rm = TRUE)))
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
                
                if(nrow(d.list.0)>0)  
                {
                  dummy.zeros.0=d.list.0[,-match(vars.head,names(d.list.0))]
                  dummy.zeros.0[,]=0
                  names(d.list.0)[id.var.nms.f]=paste('f',names(d.list.0)[id.var.nms.f],sep='')  
                  names(dummy.zeros.0)=paste('m',names(dummy.zeros.0),sep='')  
                  dummy.Size.compo.SS.format_Sex0=bind_cols(d.list.0,dummy.zeros.0)%>%
                    mutate(Sex=as.numeric(Sex))
                  
                  
                }
                if(nrow(d.list.f)>0)
                {
                  dummy.zeros.f=d.list.f[,-match(vars.head,names(d.list.f))]
                  dummy.zeros.f[,]=0
                  dummy.zeros.m=d.list.m[,-match(vars.head,names(d.list.m))]
                  dummy.zeros.m[,]=0
                  
                  Fem.Da=d.list.f%>%data.frame(check.names = FALSE)       
                  names(Fem.Da)[id.var.nms.f]=paste('f',names(Fem.Da)[id.var.nms.f],sep='')
                  Fem.Da.0=dummy.zeros.f%>%data.frame(check.names = FALSE)
                  names(Fem.Da.0)=paste('m',names(Fem.Da.0),sep='')
                  Fem.Da=bind_cols(Fem.Da,Fem.Da.0)
                  
                  Mal.hdr=d.list.m[,match(vars.head,names(d.list.m))]%>%data.frame(check.names = FALSE)
                  Mal.Da.0=dummy.zeros.m%>%data.frame(check.names = FALSE)
                  names(Mal.Da.0)=paste('f',names(Mal.Da.0),sep='')
                  Mal.Da=d.list.m[,-match(vars.head,names(d.list.m))]%>%data.frame(check.names = FALSE)
                  names(Mal.Da)=paste('m',names(Mal.Da),sep='')
                  Mal.Da=bind_cols(Mal.hdr,Mal.Da.0,Mal.Da)
                  
                  dummy.Size.compo.SS.format_Sex=bind_rows(Fem.Da,Mal.Da)%>%
                    mutate(Sex=as.numeric(Sex)) 
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
                  dummy.Size.compo.SS.format=bind_rows(dummy.Size.compo.SS.format_Sex0,dummy.Size.compo.SS.format_Sex)
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
          if(SS.sex.length.type==3 & !Neim%in%do.not.set.SS.sex.length.type.to.3)
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
          CPUE.zone=NULL
          if(!is.null(CPUE))
          {
            #Create cpue by zone list
            CPUE.zone=CPUE
            
            #Keep relevant abundance objects
            DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
            if(length(DROP)>0) CPUE=CPUE[-DROP]
            DROP.zone=match(c('observer',"TDGDLF.daily","TDGDLF.monthly"),names(CPUE.zone))
            DROP.zone=subset(DROP.zone,!is.na(DROP.zone))
            if(length(DROP.zone)>0) CPUE.zone=CPUE.zone[-DROP.zone]
            
            #Reset very low CVs 
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
            
            #Keep relevant indices and sensitivity tests
            if(Neim%in%tdgdlf_not.representative & any(grepl("TDGDLF",names(CPUE))))  
            {
              CPUE.test.tdgdlf=CPUE   
              CPUE.zone.test.tdgdlf=CPUE.zone
              CPUE=CPUE[-grep("TDGDLF.monthly",names(CPUE))]
              CPUE=CPUE[-grep("TDGDLF.daily",names(CPUE))]
              CPUE.zone=CPUE.zone[-grep("TDGDLF.monthly",names(CPUE.zone))]
              CPUE.zone=CPUE.zone[-grep("TDGDLF.daily",names(CPUE.zone))]
            }
            if(Neim%in%tdgdlf_monthly_not.representative & "TDGDLF.monthly"%in%names(CPUE))
            {
              CPUE.test.tdgdlf=CPUE   
              CPUE.zone.test.tdgdlf=CPUE.zone
              CPUE=CPUE[-grep("TDGDLF.monthly",names(CPUE))]
              CPUE.zone=CPUE.zone[-grep("TDGDLF.monthly",names(CPUE.zone))]
            }
            if(Neim%in%names(drop.dodgy.cpue))
            {
              diszone.yr=drop.dodgy.cpue[[match(Neim,names(drop.dodgy.cpue))]]
              dis.data.set=unique(word(diszone.yr, 1, sep = "-"))
              dis.yrs=unique(word(diszone.yr, 2, sep = "-"))
              for(xx in 1:length(dis.data.set))
              {
                for(yy in 1:length(dis.yrs))
                {
                  CPUE.zone[[paste0('TDGDLF.',dis.data.set[xx])]]=CPUE.zone[[paste0('TDGDLF.',dis.data.set[xx])]]%>%
                    filter(!yr.f==as.numeric(dis.yrs))
                }
              }
            }
            if(Neim%in%names(tdgdlf_daily_not.representative))
            {
              id.drop=paste0('daily.',tdgdlf_daily_not.representative[[Neim]])
              CPUE.zone=CPUE.zone[-grep(paste(id.drop,collapse='|'),names(CPUE.zone))] 
            }
            if(Neim%in%names(tdgdlf_zone_not.representative))   
            {
              id.drop=tdgdlf_zone_not.representative[[Neim]]
              CPUE.zone=CPUE.zone[-grep(paste(id.drop,collapse='|'),names(CPUE.zone))] 
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
            #1. Extract data
            releases=Species.data[[i]]$Con_tag_SS.format_releases%>%
              filter(Yr.rel<=Last.yr.ktch.numeric)
            recaptures=Species.data[[i]]$Con_tag_SS.format_recaptures%>%
              filter(Yr.rec<=Last.yr.ktch.numeric)
            
            #2. Keep relevant finyear zones
            dis.yrs.tag=Use.these.tag.year_zones[[match(names(Species.data)[i],names(Use.these.tag.year_zones))]]
            releases=releases%>%
              mutate(dummy=paste(Yr.rel,Rel.zone))%>%
              filter(dummy%in%dis.yrs.tag)%>%
              dplyr::select(-dummy)
            recaptures=recaptures%>%filter(Tag.group%in%unique(releases$Tag.group))
            
            #3. Keep tag groups from years accounting for X% of recaptures or not
            if(!is.null(use.tag.rec.yrs.percent.rec))
            {
              Table.yr.releases=table(releases$Yr.rel)
              vec=cumsum(Table.yr.releases)/sum(Table.yr.releases)
              Last.yr.rel=names(vec[which.min(abs(vec - use.tag.rec.yrs.percent.rec))])
              
              Table.yr.recaptures=table(recaptures$Yr.rec)
              vec=cumsum(Table.yr.recaptures)/sum(Table.yr.recaptures)
              Last.yr.rec=names(vec[which.min(abs(vec - use.tag.rec.yrs.percent.rec))])
              recaptures=recaptures%>%filter(Yr.rec<=as.numeric(Last.yr.rec))
            }
            
            #4. Use minimum number of observations per release tag group  
            releases=releases%>%filter(N.release>=Min.annual.Tag.group)
            recaptures=recaptures%>%filter(Tag.group%in%unique(releases$Tag.group))
            
            #5. Remove dodgy releases and recaptures beyond gillnet selectivity
            if(Neim%in%names(remove.dodgy.tag.group_age))
            {
              drop.tag.group=releases%>%filter(Age>=remove.dodgy.tag.group_age[[Neim]])
              recaptures=recaptures%>%filter(!Tag.group%in%drop.tag.group$Tag.group)
              releases=releases%>%filter(!Tag.group%in%drop.tag.group$Tag.group)
            }
            
            #5. Allocate West, Zone1 and Zone 2 to Southern 1 or 2
            releases=releases%>%
              mutate(Rel.zone=case_when(Yr.rel<=2005 ~'Southern.shark_1',
                                        Yr.rel>2005 ~'Southern.shark_2'))
            recaptures=recaptures%>%
              mutate(Rec.zone=case_when(Yr.rec<=2005 ~'Southern.shark_1',
                                        Yr.rec>2005 ~'Southern.shark_2'))
            
            #7. Recalculate TagGroup
            if(taggroup.recalculate)
            {
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
            }
            
            #8. Group sex  
            if(Neim%in%taggroup.sex.combined)   
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
            get.fleet=recaptures%>%distinct(Yr.rec,Rec.zone)
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
            
            if(Reporting.rate.type[[Neim]]=='published') Initial.reporting.rate=Species.data[[i]]$Con_tag_non_reporting_from_F.estimation.R_   
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
                fill(Fleet, .direction = "down") 
              
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
                diKay=0  
                if(nrow(d.init.rep)>1)
                {
                  Init.rep=d.init.rep$Reporting[1]
                  fit_nls <- nls(Reporting ~ Init.rep * exp(-k * time), 
                                 data = d.init.rep, 
                                 start = list(k = 0.01))
                  diKay=round(coef(fit_nls),4)              
                  if(!allow.increase.tag.rep.rate) diKay=max(0,diKay)
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
                filter(Finyear==Initial.reporting.rate$Finyear[which.min(abs(Initial.reporting.rate$Finyear - min(releases$Yr.rel)))]),   #NEW
              Reporting.rate.decay=Reporting.rate.decay,
              overdispersion=Life.history$SS_overdispersion,       
              mixing_latency_period=SS_mixing_latency_period, 
              max_periods=ceiling((max(recaptures$Yr.rec)-min(releases$Yr.rel))*Extend.mx.period))    
            
            rm(releases,recaptures,Chronic.tag.loss,Initial.reporting.rate,Reporting.rate.decay)
          }
          # by zones
          Tags.SS.format.zone=NULL
          if(names(Species.data)[i]%in%use.tag.data)
          {
            #1. Extract data
            releases=Species.data[[i]]$Con_tag_SS.format_releases%>%
              filter(Yr.rel<=Last.yr.ktch.numeric)
            recaptures=Species.data[[i]]$Con_tag_SS.format_recaptures%>%
              filter(Yr.rec<=Last.yr.ktch.numeric)
            
            get.fleet=recaptures%>%distinct(Yr.rec,Rec.zone)
            
            #2. Keep relevant finyear zones
            dis.yrs.tag=Use.these.tag.year_zones[[match(names(Species.data)[i],names(Use.these.tag.year_zones))]]
            releases=releases%>%
              mutate(dummy=paste(Yr.rel,Rel.zone))%>%
              filter(dummy%in%dis.yrs.tag)%>%
              dplyr::select(-dummy)
            recaptures=recaptures%>%filter(Tag.group%in%unique(releases$Tag.group))
            
            #3. Keep tag groups from years accounting for X% of recaptures or not
            if(!is.null(use.tag.rec.yrs.percent.rec))
            {
              Table.yr.releases=table(releases$Yr.rel)
              vec=cumsum(Table.yr.releases)/sum(Table.yr.releases)
              Last.yr.rel=names(vec[which.min(abs(vec - use.tag.rec.yrs.percent.rec))])
              
              Table.yr.recaptures=table(recaptures$Yr.rec)
              vec=cumsum(Table.yr.recaptures)/sum(Table.yr.recaptures)
              Last.yr.rec=names(vec[which.min(abs(vec - use.tag.rec.yrs.percent.rec))])
              recaptures=recaptures%>%filter(Yr.rec<=as.numeric(Last.yr.rec))
            }
            
            #4. Use minimum number of observations per release tag group  
            releases=releases%>%filter(N.release>=Min.annual.Tag.group)
            recaptures=recaptures%>%filter(Tag.group%in%unique(releases$Tag.group))
            
            #5. Remove dodgy releases and recaptures beyond gillnet selectivity
            if(Neim%in%names(remove.dodgy.tag.group_age))
            {
              drop.tag.group=releases%>%filter(Age>=remove.dodgy.tag.group_age[[Neim]])
              recaptures=recaptures%>%filter(!Tag.group%in%drop.tag.group$Tag.group)
              releases=releases%>%filter(!Tag.group%in%drop.tag.group$Tag.group)
            }
            
            #6. Group sex  
            if(Neim%in%taggroup.sex.combined)   
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
            
            #explore tag groups consistency with gear selectivity
            if(First.run=='YES')
            {
              #add mean length
              TLzero=with(Life.history,Lzero*a_FL.to.TL+b_FL.to.TL)
              K=Life.history$Growth.F$k
              K.m=Life.history$Growth.M$k
              TLINF=with(Life.history,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)
              TLINF.m=with(Life.history,Growth.M$FL_inf*a_FL.to.TL+b_FL.to.TL)
              dd=releases%>%
                mutate(Mean.length=ifelse(Sex==1,TLINF-(TLINF-TLzero)*exp(-K*Age),
                                          TLINF.m-(TLINF.m-TLzero)*exp(-K.m*Age)))
              
              TL.vec=seq(TLzero,Life.history$TLmax)
              p=dd%>%
                mutate(LBL=paste0('TG ',Tag.group,' (n=',N.release,')'),
                       Yr.rel=as.character(Yr.rel))%>%
                ggplot(aes(Mean.length,N.release/max(N.release),color=Yr.rel))+
                geom_point()+
                geom_text_repel(aes(label = LBL))+
                facet_wrap(~Rel.zone,ncol=1)+xlim(0,max(TL.vec))+
                theme_PA()+ylab('Relative number of releases (points) or selectivity (line)')
              
              #add mean SS selectivity  
              if(!exists('doubleNorm24.fn')) fn.source1("SS_selectivity functions.R")
              West.sel=Life.history$SS_selectivity%>%filter(Fleet=='Southern.shark_1_West')
              Zn1.sel=Life.history$SS_selectivity%>%filter(Fleet=='Southern.shark_1_Zone1')
              Zn2.sel=Life.history$SS_selectivity%>%filter(Fleet=='Southern.shark_1_Zone2')
              
              d.sel=data.frame(Mean.length=rep(TL.vec,each=3),
                               Rel.zone=rep(c("West","Zone1","Zone2"),length(TL.vec)))%>%
                mutate(Selectivity=case_when(Rel.zone=="West"~doubleNorm24.fn(Mean.length,a=West.sel$P_1,b=West.sel$P_2, c=West.sel$P_3, d=West.sel$P_4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE),
                                             Rel.zone=="Zone1"~doubleNorm24.fn(Mean.length,a=Zn1.sel$P_1,b=Zn1.sel$P_2, c=Zn1.sel$P_3, d=Zn1.sel$P_4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE),
                                             Rel.zone=="Zone2"~doubleNorm24.fn(Mean.length,a=Zn2.sel$P_1,b=Zn2.sel$P_2, c=Zn2.sel$P_3, d=Zn2.sel$P_4, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE),
                                             TRUE~0),
                       Yr.rel=as.character(max(dd$Yr.rel)))
              p=p+
                geom_line(data=d.sel,aes(Mean.length,Selectivity,group = 1))
              base::print(p)
              ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                           capitalize(Life.history$Name),"/",AssessYr,"/1_Inputs/Visualise data/Tagging_Tag groups and selectivity_area.as.fleets.tiff",sep=''),
                     width = 6,height = 8,compression = "lzw")
            }
            
            #7. Allocate West, Zone1 and Zone 2 to Southern 1 or 2
            releases=releases%>%
              mutate(Rel.zone=case_when(Yr.rel<=2005 ~'Southern.shark_1',
                                        Yr.rel>2005 ~'Southern.shark_2'))%>%
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
            
            #8. Recalculate TagGroup
            if(taggroup.recalculate)
            {
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
            }
            
            #9. Final calcs
            releases=releases%>%
              rename(Area=Rel.zone)%>%
              mutate(Area=1)
            #get.fleet=recaptures%>%distinct(Yr.rec,Rec.zone)
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
                diKay=0  
                if(nrow(d.init.rep)>1)
                {
                  Init.rep=d.init.rep$Reporting[1]
                  fit_nls <- nls(Reporting ~ Init.rep * exp(-k * time), 
                                 data = d.init.rep, 
                                 start = list(k = 0.01))
                  diKay=round(coef(fit_nls),4)              
                  if(!allow.increase.tag.rep.rate) diKay=max(0,diKay)
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
                             capitalize(Life.history$Name),"/",AssessYr,"/1_Inputs/Visualise data/Tagging_report rate decay_area.as.fleets.tiff",sep=''),
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
              overdispersion=Life.history$SS_overdispersion,       
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
            if(exists('CPUE.test.tdgdlf'))
            {
              len.cpue.test=length(CPUE.test.tdgdlf)
              len.cpue.zone.test=length(CPUE.zone.test.tdgdlf)
              Abundance.SS.format.test.tdgdlf=Abundance.SS.format
              Abundance.SS.format.zone.test.tdgdlf=Abundance.SS.format.zone
              if(len.cpue.test>0 | len.cpue.zone.test>0)
              {
                MAX.CV=Life.history$MAX.CV
                #zones together
                if(len.cpue.test>0)
                {
                  for(x in 1:len.cpue.test)    
                  {
                    nm=names(CPUE.test.tdgdlf)[x]
                    if(nm=="NSF") nm="Northern.shark"
                    if(nm=="TDGDLF.monthly") nm="Southern.shark_1"
                    if(nm=="TDGDLF.daily") nm="Southern.shark_2"
                    
                    dd=CPUE.test.tdgdlf[[x]][,grep(paste(c('yr.f','Mean','MeAn','CV'),collapse="|"),names(CPUE.test.tdgdlf[[x]]))]%>%
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
                    CPUE.test.tdgdlf[[x]]=dd%>%filter(!is.na(Mean))
                  }
                  Abundance.SS.format.test.tdgdlf=do.call(rbind,CPUE.test.tdgdlf)%>%
                    relocate(Year,seas,index,Mean,CV)%>%
                    arrange(index,Year)
                  
                }
                #by zone
                if(len.cpue.zone.test>0)
                {
                  for(x in 1:len.cpue.zone.test)    
                  {
                    nm=names(CPUE.zone.test.tdgdlf)[x]
                    if(nm=="NSF") nm="Northern.shark"
                    if(nm=="TDGDLF.monthly.West") nm="Southern.shark_1_West"
                    if(nm=="TDGDLF.monthly.Zone1") nm="Southern.shark_1_Zone1"
                    if(nm=="TDGDLF.monthly.Zone2") nm="Southern.shark_1_Zone2"
                    if(nm=="TDGDLF.daily.West") nm="Southern.shark_2_West"
                    if(nm=="TDGDLF.daily.Zone1") nm="Southern.shark_2_Zone1"
                    if(nm=="TDGDLF.daily.Zone2") nm="Southern.shark_2_Zone2"
                    dd=CPUE.zone.test.tdgdlf[[x]][,grep(paste(c('yr.f','Mean','MeAn','CV'),collapse="|"),names(CPUE.zone.test.tdgdlf[[x]]))]%>%
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
                    CPUE.zone.test.tdgdlf[[x]]=dd%>%filter(!is.na(Mean))
                  }
                  Abundance.SS.format.zone.test.tdgdlf=do.call(rbind,CPUE.zone.test.tdgdlf)%>%
                    relocate(Year,seas,index,Mean,CV)%>%
                    arrange(index,Year)
                }
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
            
            #drop irrelevant lambda scenarios  
            if('like_comp.w'%in%names(Scens))
            {
              if(is.null(Abundance.SS.format) & is.null(Abundance.SS.format.zone))
              {
                DROP.ABU=which(Scens$like_comp.w==1)
                if(length(DROP.ABU)>0)
                {
                  Scens=Scens[-DROP.ABU,]
                  List.sp[[i]]$Sens.test$SS=List.sp[[i]]$Sens.test$SS[-DROP.ABU,]
                }
                
              }
              if(is.null(Size.compo.SS.format) & is.null(Size.compo.SS.format.zone))
              {
                DROP.SIZ=which(Scens$like_comp.w==4)
                if(length(DROP.SIZ)>0)
                {
                  Scens=Scens[-DROP.SIZ,]
                  List.sp[[i]]$Sens.test$SS=List.sp[[i]]$Sens.test$SS[-DROP.SIZ,]
                }
                
              }
              Scens$Scenario=paste0('S',1:nrow(Scens))
              List.sp[[i]]$Sens.test$SS$Scenario=paste0('S',1:nrow(List.sp[[i]]$Sens.test$SS))
              Store.sens=Store.sens[Scens$Scenario]
            }
            
            #Run over each scenario
            for(s in 1:length(Store.sens))
            {
              print(paste0('Creating input files for model: ',names(Age.based)[w],' --- Species: ' ,Neim,' --- Scenario: ',names(Store.sens)[s]))
              
              this.wd1=paste(this.wd,names(Store.sens)[s],sep='/')
              if(!dir.exists(this.wd1))dir.create(this.wd1)
              
              #cpues
              Abund.single.area=Abundance.SS.format
              Abund.areas.as.fleets=Abundance.SS.format.zone
              if(Scens$Use.dropped.cpue[s]=='Yes') 
              {
                Abund.single.area=Abundance.SS.format.test.tdgdlf
                Abund.areas.as.fleets=Abundance.SS.format.zone.test.tdgdlf
              }
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
                sens.test.cpue=sens.test.cpue.zone=NULL
                if(exists('CPUE.test.tdgdlf')) sens.test.cpue=CPUE.test.tdgdlf
                if(exists('CPUE.zone.test.tdgdlf')) sens.test.cpue.zone=CPUE.zone.test.tdgdlf
                if('TDGDLF.daily'%in%names(CPUE) | any(grep('TDGDLF.daily',names(CPUE.zone))) |
                   'TDGDLF.daily'%in%names(sens.test.cpue) | any(grep('TDGDLF.daily',names(sens.test.cpue.zone))))
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
              
              #Completely remove a data type
              #CPUE  
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
              #Length.comps
              Size.comp.single.area=Size.compo.SS.format
              Size.comp.areas.as.fleets=Size.compo.SS.format.zone
              if(Scens$Length.comps[s]=='No') 
              {
                Size.comp.single.area=Size.comp.areas.as.fleets=NULL
                Scens$Use.male.sel.offset[s]="No"
              }
              #Mean.body
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
                Life.history$Apical.prop.male=data.frame(Zone='',Prop=1-Apical.prop.female$w_mean)
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
                Apical.prop.male=Apical.prop.female.zone%>%
                  mutate(Zone=sub(".*_", "", fishry),
                         Prop=round(1-w_mean,2))%>%
                  dplyr::select(Zone,Prop)
                if(Neim%in%c("gummy shark","whiskery shark"))
                {
                  Apical.prop.male[match('West',Apical.prop.male$Zone),'Prop']=Apical.prop.male[match('Zone1',Apical.prop.male$Zone),'Prop']
                }
                Life.history$Apical.prop.male=Apical.prop.male
                
                
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
                new.ktch=as.numeric(unlist(KAtch[des.yrs,des.flt])-(Indo.IUU$LIVEWT.c/1000)+(Indo.IUU.app$LIVEWT.c/1000))
                
                p=KAtch[des.yrs,c(3,des.flt)]%>%
                  data.frame()%>%
                  mutate(Apprehensions.only=new.ktch)%>%
                  rename(max.Apprehensions.Forfeitures=X2)%>%
                  gather(Method,Tons,-finyear)
                p%>%
                  ggplot(aes(finyear,Tons,color=Method))+
                  geom_point()+geom_line()+ylim(0,NA)+
                  theme_PA()+
                  theme(legend.position = 'top',
                        legend.title=element_blank())
                ggsave(paste0(this.wd,"/Catch Other fleet_Indo catch recons_Apprehensions or Forfeitures.tiff"),width=7,height=6,compression = "lzw")
                
                KAtch[des.yrs,des.flt]=new.ktch
              }
              
              #a.4 set recruitment deviations parameters
              #note: Align with data-rich years , mostly comps that have info on recruitment 
              #      These pars determine how the model applies Bias Adjustment to ensure the mean recruitment
              #         isn't skewed when shifting between "data-poor" and "data-rich" periods.
              #      As data becomes weaker (at the very beginning or very end of the time series), the model can't estimate the 
              #         full variance, so you "ramp" the adjustment down
              #      recdev_early_start: first year the model starts estimating deviations before the "Main" period begins.
              #                           This allows the initial age structure (at the model start year) to deviate from equilibrium,
              #                           making the starting population more realistic
              #       First.yr.main.rec.dev: first year of (informative) data-rich years
              #       last.year.of.main.recr_devs: last year of (informative) data-rich years to reliably estimate rec devs
              #       last_yr_nobias_adj_in_MPD: the point in the early period where the model begins to transition 
              #                             from zero adjustment to full adjustment
              #       first_yr_fullbias_adj_in_MPD: the year where the model starts applying the maximum bias adjustment (set to a few years 
              #                             into data-rich period), usually when your composition data becomes consistent
              #       last_yr_fullbias_adj_in_MPD: the last year of "strong" data. After this, the model starts ramping the adjustment down  
              #                             because the latest recruits haven't been seen by the survey/fishery enough times to be well-estimated
              #       end_yr_for_ramp_in_MPD: The year the ramp hits zero again (often the end of the forecast or the last year of the model)
              Abund1=Abund
              if(!is.null(Abund1)) Abund1=Abund1%>%rename_with(tolower)
              
              Obs.list.rec.dev=list(Abund1=Abund1,Size.com=Size.com,meanbody=meanbody)
              rec_dev_data.types1=rec_dev_data.types
              if(is.null(Size.com) & is.null(meanbody) & !is.null(Abund1)) rec_dev_data.types1='Abund1'
              Obs.list.rec.dev=Obs.list.rec.dev[rec_dev_data.types1]
              Min.yr.obs=min(unlist(lapply(Obs.list.rec.dev,function(x) if(!is.null(x))min(x$year)))) 
              Max.yr.obs=max(unlist(lapply(Obs.list.rec.dev,function(x) if(!is.null(x))max(x$year))))
              
              
              Life.history$MainRdevYrLast=min(Max.yr.obs,max(KAtch$finyear,na.rm=T))
              if(Life.history$First.yr.main.rec.dev_buffer)
              {
                if(is.numeric(Life.history$recdev_early_start))
                {
                  if(Life.history$recdev_early_start>0) Life.history$recdev_early_start=Min.yr.obs-Life.history$recdev_early_start
                }
                if(Life.history$recdev_early_start=='mat')
                {
                  Min.mat=max(Early_rec_dev_start.min.yrs,round(min(Life.history$Age.50.mat)))
                  Life.history$recdev_early_start=Min.yr.obs-Min.mat
                }
              }
              if(Life.history$First.yr.main.rec.dev=='min.obs') MainRdevYrFirst=Min.yr.obs
              if(Life.history$First.yr.main.rec.dev=='min.ktch') MainRdevYrFirst=min(ktch$finyear)
              Life.history$MainRdevYrFirst=MainRdevYrFirst
              
              #a.5 Reset rec pars for tuning
              if(Scens$Scenario[s]=='S1' & Tune.SS.model)
              {
                try.dis.sigmaR=tuning_sigmaR
                if(Neim%in%sp.low.productivity) try.dis.sigmaR=tuning_sigmaR_low.prod 
                Life.history$SR_sigmaR=Scens[s,]$SR_sigmaR=try.dis.sigmaR
                Life.history$RecDev_Phase=3
                
                #Ramp:
                # The model linearly interpolates the adjustment fraction between these four year-markers
                
                #The last year of the early recruitment period where no bias adjustment (0%) is applied. 
                # Typically used for very early years with no data
                Life.history$last_early_yr_nobias_adj_in_MPD=MainRdevYrFirst-1
                
                #The year when the model transitions to full bias adjustment (100%). 
                # This should align with the start of informative composition data
                Life.history$first_yr_fullbias_adj_in_MPD=MainRdevYrFirst+1 
                
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
              if(First.run=='YES')  
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
              
              #a.8 Completely remove Tag data
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
              
              #a.11 reset CV's if required by scenario
              if(!is.null(Abund))
              {
                if(!Scens[s,'extra.CV']=='NO')
                {
                  ID.cvS=grep(Scens[s,'extra.CV'],rownames(Abund))
                  Abund[ID.cvS,'CV']=Abund[ID.cvS,'CV']*extra.CV.factor
                }
              }
              
              #a.12 remove all 0s length comps if artificially created
              if(!is.null(Size.com))
              {
                Size.com=Size.com%>%data.frame()
                id.all.zeros=rowSums(Size.com[,-match(c('year','Seas','Fleet','Sex','Part','Nsamp'),names(Size.com))])
                names(id.all.zeros)=1:length(id.all.zeros)
                id.all.zeros=names(which(id.all.zeros==0))
                if(length(id.all.zeros)>0) Size.com=Size.com[-id.all.zeros,]
              }
              
              #a.13 create file  
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
              clear.log("Abund")
              clear.log("Size.com")
              clear.log("meanbody")
              clear.log("tags")
              clear.log("Var.ad")
            }  #end s loop
            
          } # end 'if(len.cpue>0 | etc'
          
          clear.log('CPUE')
          clear.log('CPUE.zone')
          clear.log('CPUE.test.tdgdlf')
          clear.log('CPUE.zone.test.tdgdlf')
          clear.log("Var.ad.factr")
          clear.log("Var.ad.factr.zone")
          clear.log("Combined.ktch")
          clear.log("ktch.zone")
          clear.log("ktch")
          clear.log("add.ct.future")
          clear.log("add.ct.future.zone")
        } # end '!is.null(Catch.rate.series[[i]]) | etc'
      } #end i loop
    }
  } #end w loop
}

#---Run models and get outputs -------------------------------------------------
# first remove irrelevant scenarios
for(i in 1:N.sp)
{
  this.wd=paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated",sep='')
  run.scens=list.dirs(path = this.wd, full.names = FALSE, recursive = FALSE)
  run.scens <- run.scens[run.scens != ""]
  List.sp[[i]]$Sens.test$SS=List.sp[[i]]$Sens.test$SS%>%filter(Scenario%in%run.scens)
}
# then run models
tic("timer")
for(w in 1:n.SS)
{
  # SS3
  if(names(Age.based)[w]=="SS") #takes 1000 secs per species per scenario (with Hessian and MC simulation)
  {
    #Parallel processing over combination of species and scenarios
    if(do.parallel.SS) 
    {
      #1. Create combination of species and scenarios matrix
      Species.scenarios<- vector('list',length(List.sp))
      names(Species.scenarios) = names(List.sp)
      for(v in 1:length(Species.scenarios))  
      {
        Species.scenarios[[v]]=data.frame(Species=names(Species.scenarios)[v],
                                          Scenario=List.sp[[names(Species.scenarios)[v]]]$Sens.test$SS$Scenario)
      }
      Species.scenarios=do.call(rbind,Species.scenarios)
      if(!do.all.sensitivity.tests) Species.scenarios=Species.scenarios%>%filter(Scenario=='S1')
      N.sp.scen=nrow(Species.scenarios)
      Specs.scens=paste(Species.scenarios$Species, Species.scenarios$Scenario,sep="_")
      
      #2. Estimate parameters and generate outputs
      set.seed(1234)
      prog.labl=paste(Species.scenarios$Species,"----------------Scenario:",Species.scenarios$Scenario)
      progress <- function(n) cat(sprintf(": SS3 fit complete for %d-------------", n),prog.labl[n],"\n")
      opts <- list(progress = progress)
      cl <- makeCluster(detectCores()-1)
      registerDoSNOW(cl)
      out.species=foreach(i=1:N.sp.scen,   
                          .options.snow = opts,
                          .errorhandling = 'pass',
                          .packages=c('gridExtra','Hmisc','JABBA','strex',
                                      'TruncatedDistributions','tidyverse','r4ss',
                                      'mvtnorm','ggrepel','ss3diags','ggpubr')) %dopar%
        {
          Neim=Species.scenarios$Species[i]
          Scen=Species.scenarios$Scenario[i]
          
          #a set path and get life history data
          this.wd=paste(HandL.out,capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
          this.wd1=paste(this.wd,Scen,sep='/')
          Life.history=List.sp[[Neim]]
          
          #b. Run SS3
          if(file.exists(paste0(this.wd1,'/control.ctl')))
          {
            #b.1 find LnRo init value
            if(Find_Init_LnRo)  
            {
              start <- r4ss::SS_readstarter(file = file.path(this.wd1, "starter.ss"), verbose = FALSE)
              start$last_estimation_phase=0
              r4ss::SS_writestarter(start, dir = this.wd1, overwrite = TRUE,verbose = FALSE)
              fn.run.SS(where.inputs=this.wd1,
                        where.exe=Where.exe,
                        args="-nohess") 
              Report=SS_output(this.wd1,covar=F,forecast=F,readwt=F)  
              LnRo_init_value=Report$timeseries%>%filter(Era=='VIRG')%>%pull(Bio_all) #JABBA K= 6800 tonnes 
              write.csv(LnRo_init_value,paste(this.wd1,'LnRo_init_value.csv',sep='/'),row.names = F)
              rm(Report)
            }
            
            #b.2 find tag overdispersion
            #note: update 'SS3.Tagging_pars.csv' with the Phi value
            if(Scen=='S1' & Find_overdispersion)
            {
              start <- r4ss::SS_readstarter(file = file.path(this.wd1, "starter.ss"), verbose = FALSE)
              dat <- r4ss::SS_readdat(file = file.path(this.wd1, start$datfile), verbose = FALSE)
              if(!is.null(dat$tag_releases))
              {
                ctl <- r4ss::SS_readctl(file = file.path(this.wd1, start$ctlfile), verbose = FALSE, use_datlist = TRUE, datlist = dat)
                ctl$TG_overdispersion=ctl$TG_overdispersion%>%
                  mutate(LO=1,
                         INIT=1.001,
                         PRIOR=1)
                r4ss::SS_writectl(ctl, outfile = file.path(this.wd1, start$ctlfile), overwrite = TRUE, verbose = FALSE)
                fn.run.SS(where.inputs=this.wd1,
                          where.exe=Where.exe,
                          args="-nohess") 
                Report=SS_output(this.wd1,covar=F,forecast=F,readwt=F)
                Phi=fn.SS3_calc.tag.overdispersion(Report=Report)
                write.csv(Phi,paste(this.wd1,'Tag_overdispersion.csv',sep='/'),row.names = F)
                rm(Report) 
              }
            }
            
            #b.3 tune model and calculate RAMP years (takes 3444 secs for indicator species)
            #note: var adjust and ramp already reset in '#a.5 Reset rec pars for tuning'
            #      update ramp years in 'SS3.Recruitment.inputs.csv' using 'Ramp_years_first round'  
            #      and sample sizes in 'SS3.tune_size_comp_effective_sample.csv' if single area model or
            #         'SS3.tune_size_comp_effective_sample_spatial.csv' if areas as fleets or spatial model
            #          using 'Tuned_size_comp.csv'.
            if(Scen=='S1' & Tune.SS.model)
            {
              tune.folder=paste(this.wd,'tuning',sep='/')
              if(!file.exists(file.path(tune.folder))) dir.create(file.path(tune.folder))
              
              #1st. Tune ramp years (blue and red lines should match)
              fn.run.SS(where.inputs=this.wd1,
                        where.exe=Where.exe,
                        args='')
              Report=SS_output(this.wd1)
              tiff(file=paste(tune.folder,'Ramp_years_first round.tiff',sep='/'),
                   width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
              ramp_years=SS_fitbiasramp(Report) 
              dev.off()
              out=ramp_years$df
              out=rbind(out,data.frame(value=unique(Report$sigma_R_info$alternative_sigma_R),label='Alternative_sigma_R'))
              write.csv(out,paste(tune.folder,'Ramp_years_first round.csv',sep='/'),row.names = F)
              these.plots=c(1:7,10,11,16,26)  #biol, selectivity, timeseries,rec devs,S-R,catch,mean weight, indices, size comp
              SS_plots(Report, plot=these.plots, png=T,dir=tune.folder,printfolder = "1_plots_not tuned")
              Likelihoods.not.tuned=Report$likelihoods_used%>%mutate(type='not tuned')
              
              #2nd. Tune composition data
              tune_info <- tune_comps(option = "Francis",
                                      niters_tuning = 1,
                                      dir = this.wd1,
                                      exe=Where.exe,
                                      allow_up_tuning = TRUE,
                                      verbose = FALSE)
              Tuned.var.adjust=tune_info$weights[[1]]%>%mutate(Method='Francis')
              write.csv(Tuned.var.adjust,paste(tune.folder,'Tuned_size_comp.csv',sep='/'),row.names = F)
              
              #3rd. Re run model with updated tuned pars
              #3.1 bring in SS input files
              start <- r4ss::SS_readstarter(file = file.path(this.wd1, "starter.ss"), verbose = FALSE)
              dat <- r4ss::SS_readdat(file = file.path(this.wd1, start$datfile), verbose = FALSE)
              ctl <- r4ss::SS_readctl(file = file.path(this.wd1, start$ctlfile), verbose = FALSE, use_datlist = TRUE, datlist = dat)
              
              #3.2 replace var adj factor with tuned values
              ctl$Variance_adjustment_list=with(Tuned.var.adjust,data.frame(factor=factor,fleet=fleet,value=value)) 
              
              #3.3 replace ramp years
              ctl$last_early_yr_nobias_adj= out%>%filter(grepl('last_early_yr_nobias_adj',label))%>%pull(value)%>%as.numeric()
              ctl$first_yr_fullbias_adj= out%>%filter(grepl('first_yr_fullbias_adj',label))%>%pull(value)%>%as.numeric() 
              ctl$last_yr_fullbias_adj= out%>%filter(grepl('last_yr_fullbias_adj',label))%>%pull(value)%>%as.numeric() 
              ctl$first_recent_yr_nobias_adj= out%>%filter(grepl('first_recent_yr_nobias_adj',label))%>%pull(value) %>%as.numeric()
              ctl$max_bias_adj= out%>%filter(grepl('max_bias_adj',label))%>%pull(value)%>%as.numeric()
              ctl$recdev_early_start=0
              
              #3.4 replace sigmaR 
              ctl$SR_parms$INIT[match('SR_sigmaR',rownames(ctl$SR_parms))]=out%>%filter(grepl('Alternative_sigma_R',label))%>%pull(value)%>%as.numeric()
              
              #3.5 export new control
              r4ss::SS_writectl(ctl, outfile = file.path(this.wd1, start$ctlfile), overwrite = TRUE, verbose = FALSE)
              
              #3.6 re run model with updated parametes
              fn.run.SS(where.inputs=this.wd1,
                        where.exe=Where.exe,
                        args='')
              Report=SS_output(this.wd1)
              tiff(file=paste(tune.folder,'Ramp_years_second round.tiff',sep='/'),
                   width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
              ramp_years=SS_fitbiasramp(Report) 
              dev.off()
              out=ramp_years$df
              out=rbind(out,data.frame(value=unique(Report$sigma_R_info$alternative_sigma_R),label='Alternative_sigma_R'))
              write.csv(out,paste(tune.folder,'Ramp_years_first round.csv',sep='/'),row.names = F)
              SS_plots(Report, plot=these.plots, png=T,dir=tune.folder,printfolder = "2_plots_tuned")
              Likelihoods.tuned=Report$likelihoods_used%>%mutate(type='tuned')
              
              #3rd compare tuned and not tuned likelihoods
              #compare tuned and not tuned likelihoods  
              rbind(Likelihoods.not.tuned%>%rownames_to_column(var = "Component"),
                    Likelihoods.tuned%>%rownames_to_column(var = "Component"))%>%
                ggplot(aes(type,values))+
                geom_bar(stat = "identity")+
                facet_wrap(~Component)+
                theme_PA()+xlab('')
              ggsave(paste(tune.folder,'Likelihoods before and after tuning.tiff',sep='/'), width = 8,height = 6,compression = "lzw")
              
              
              #4th flush
              rm(ramp_years,out,tune_info)
            }
            
            #b.4 estimate population parameters on tuned model
            Arg.sens=Arg
            if(!Scen=='S1' & !do.all.sensitivity.tests.uncertainty) Arg.sens='-nohess'
            if(Run.SS & !Tune.SS.model)
            {
              fn.run.SS(where.inputs=this.wd1,
                        where.exe=Where.exe,
                        args=Arg.sens) 
            }
            
            #c. Bring in outputs, plot Report and store estimates
            if(file.exists(paste0(this.wd1,'/Report.sso')))
            {
              #c.1 Bring report
              COVAR=FORECAST=FALSE
              if(Arg.sens=="") COVAR=TRUE
              if("SS"%in%future.models) FORECAST=TRUE
              Report=SS_output(this.wd1,covar=COVAR,forecast=FORECAST,readwt=F)
              
              #c.2. plot
              if(run_SS_plots)
              {
                this.plot=1:26
                if(!is.null(Report[['tagdbase1']]))
                {
                  dumii=data.frame(Exp=Report$tagdbase1$Exp,
                                   Expected=rep(Report$tagdbase2$Exp,each = max(Report$tagdbase1$Fleet)))%>%
                    mutate(Numbers_Exp = round(Exp * Expected))
                  if(isTRUE(unique(dumii$Numbers_Exp)==0)) this.plot=this.plot[-21]  #remove Tag plot if this condition
                }
                SS_plots(Report,plot=this.plot,  png=T)
              }
              
              #c.3. Store estimates and likelihoods  
              Estims=Report[["estimated_non_dev_parameters"]]
              Out.estimates=Estims%>%
                mutate(Par=rownames(Estims),
                       Scenario=Scen)%>%
                relocate(Scenario,Par)%>%
                `rownames<-`( NULL )
              F.ref.points=fn.get.f.ref.points(Report)                 
              dummi.F=Out.estimates[1:length(F.ref.points),]
              dummi.F[,]=NA
              dummi.F=dummi.F%>%
                mutate(Scenario=unique(Out.estimates$Scenario),
                       Par=names(F.ref.points),
                       Value=unlist(F.ref.points))
              Out.estimates=rbind(Out.estimates,dummi.F)
              Likelihoods=Report[["likelihoods_used"]]
              Out.likelihoods=Likelihoods%>%
                mutate(Likelihood=rownames(Likelihoods),
                       Scenario=Scen)%>%
                relocate(Scenario,Likelihood)%>%
                `rownames<-`( NULL )
              
              Bratio_current=paste0('Bratio_',max(Report$catch$Yr))  
              Out.quantities=Report[["derived_quants"]]%>%
                filter(Label%in%c(Bratio_current,"Dead_Catch_MSY"))%>%
                mutate(Label=ifelse(Label==Bratio_current,'Current depletion',
                                    ifelse(Label=="Dead_Catch_MSY",'MSY',
                                           Label)),
                       Species=capitalize(Neim),
                       Scenario=Scen)%>%
                rename(Median=Value,
                       SE=StdDev)%>%
                dplyr::select(Species,Label,Median,SE,Scenario)%>%
                `rownames<-`( NULL )
              
              
              #c.4. Store trajectories
              #note: uncertainty is based on asymptotic error
              
              #c.4.1 relative biomass  
              dummy=fn.integrated.mod.get.timeseries(d=Report,
                                                     mods="SS3",
                                                     Type='Depletion',
                                                     scen=Scen)
              Out.rel.biom=dummy$Dat
              Out.probs.rel.biom=dummy$Probs 
              
              #c.4.2 F 
              dummy=fn.integrated.mod.get.timeseries(d=Report,
                                                     mods="SS3",
                                                     Type='F.series',
                                                     scen=Scen)
              Out.f.series=dummy$Dat
              Out.probs.f.series=dummy$Probs
              
              #c.4.3 B/Bmsy
              dummy=fn.integrated.mod.get.timeseries(d=Report,
                                                     mods="SS3",
                                                     Type='B.Bmsy',
                                                     scen=Scen,
                                                     get.uncertainty='CV') 
              Out.B.Bmsy=dummy$Dat
              Out.probs.B.Bmsy=dummy$Probs
              Kobe.stock=dummy$Out.Kobe
              
              #c.4.4 F/Fmsy
              dummy=fn.integrated.mod.get.timeseries(d=Report,
                                                     mods="SS3",
                                                     Type='F.Fmsy',
                                                     scen=Scen,
                                                     get.uncertainty='ratio.first')
              Out.F.Fmsy=dummy$Dat
              Kobe.harvest=dummy$Out.Kobe
              
              #c.5. Calculate posterior depletion   #takes 4.5 secs per nMC simulation
              if(Scen=='S1' & SS3.run=='final' & do.MC.multi)
              {
                dd=fn.MC.sims(this.wd1,
                              nMC=nMCsims,
                              arg=Arg.no.estimation,
                              B.th=Report$derived_quants%>%filter(Label=="SSB_MSY")%>%pull(Value)/Report$derived_quants%>%filter(Label=="SSB_Virgin")%>%pull(Value),
                              scen=Scen)
                Out.rel.biom=dd$Depletion
                Out.probs.rel.biom=dd[c("probs","Reference.points")]
                Out.f.series=dd$F.series
                Out.B.Bmsy=dd$B.Bmsy
                Out.F.Fmsy=dd$F.Fmsy
                rm(dd)
              }
              
              #c.6. Add catch for display purposes 
              Add.katch=Report$catch%>%
                filter(Yr>min(Report$catch$Yr))%>%
                rename(Year=Yr)%>%
                group_by(Year)%>%
                summarise(Total=sum(Obs,na.rm=T))%>%
                ungroup()%>%
                data.frame()
              TC.future=Report$derived_quants%>%
                filter(grepl('ForeCatch',Label))%>%
                filter(!grepl('ForeCatchret',Label))
              Add.katch=bind_rows(Add.katch,
                                  data.frame(Year=as.numeric(str_remove(TC.future$Label, "ForeCatch_")),
                                             Total=TC.future%>%pull(Value)))
              Out.rel.biom=Out.rel.biom%>%
                left_join(Add.katch,by=c('year'='Year'))%>%
                rename(Catch=Total)
              
              #c.7. Kobe
              Out.Kobe.probs=NULL
              if(Scen=='S1')
              {
                Out.Kobe.probs=data.frame(stock=Kobe.stock,harvest=Kobe.harvest)%>%
                  filter(stock>=0 & harvest>=0)
              }
              
              #c.8. Plot Posterior vs Prior 
              #c.8.1 Steepness
              if(Report$parameters%>%filter(Label=='SR_BH_steep')%>%pull(Phase)>0)
              {
                pdf(paste(this.wd1,paste('Prior vs Posterior_Steepness.pdf',sep=''),sep='/'))
                fn.compare.prior.post(d=Report$parameters%>%filter(Label=='SR_BH_steep'),
                                      Par='Steepness',
                                      prior_type='beta')
                dev.off()
              }
              #c.8.2 Growth pars
              Growth.prior.post=Report$parameters%>%
                filter(grepl(paste(c('L_at_Amax','VonBert_K'),collapse='|'),Label))%>%
                filter(Phase>0)
              if(nrow(Growth.prior.post)>0)
              {
                Growth.prior.post=Growth.prior.post%>%
                  dplyr::select(Label,Value,Parm_StDev,Pr_type,Prior,Pr_SD)%>%
                  rename(Parameter=Label)
                
                df.prior=Growth.prior.post%>%
                  dplyr::select(Parameter,Prior,Pr_SD)%>%
                  rename(mu=Prior,
                         sigma=Pr_SD)
                df.post=Growth.prior.post%>%
                  dplyr::select(Parameter,Value,Parm_StDev)%>%
                  rename(mu=Value,
                         sigma=Parm_StDev)
                
                df_long.prior <- df.prior %>%
                  rowwise() %>%
                  reframe(
                    x = seq(mu - 4*sigma, mu + 4*sigma, length.out = 100),
                    y = dnorm(x, mu, sigma),
                    id = Parameter)%>%
                  mutate(type='Prior')
                
                df_long.post <- df.post %>%
                  rowwise() %>%
                  reframe(
                    x = seq(mu - 4*sigma, mu + 4*sigma, length.out = 100),
                    y = dnorm(x, mu, sigma),
                    id = Parameter)%>%
                  mutate(type='Posterior')
                
                p=rbind(df_long.prior,df_long.post)%>%
                  ggplot(aes(x, y, color = type)) +
                  geom_line(linewidth = 1.15)+
                  facet_wrap(~id,scales='free') +
                  theme_PA()+ylab('Density')+xlab('')+
                  theme(legend.position = 'top',
                        legend.title = element_blank())
                if(Neim%in%unique(List.published.growth$NAME)) 
                {
                  add.publsh=List.published.growth%>%
                    filter(NAME==Neim)%>%
                    dplyr::select(-NAME)%>%
                    rename(L_at_Amax_Fem_GP_1=TLinf.f,
                           VonBert_K_Fem_GP_1=k.f,
                           L_at_Amax_Mal_GP_1=TLinf.m,
                           VonBert_K_Mal_GP_1=k.m)%>%
                    gather(id,x,-Reference)
                  add.publsh=rbind(add.publsh,
                                   data.frame(Reference='TLmax',
                                              id='L_at_Amax_Fem_GP_1',
                                              x=Life.history$TLmax))
                  p=p+
                    geom_vline(data = add.publsh, aes(xintercept = x,color=Reference),linetype = "dotted")
                  
                  
                }
                print(p)
                ggsave(paste(this.wd1,paste('Prior vs Posterior_growth.tiff',sep=''),sep='/'),width=7,height=6,compression = "lzw")
              }
              
              #c.9. Evaluate fit diagnostics
              GoodnessFit=function.goodness.fit_SS(Rep=Report)  
              write.csv(GoodnessFit,paste(this.wd1,"/GoodnessFit_",Neim,".csv",sep=''))
              
              #c.10. Check Life history invariants (M/k=1.5; L50/Linf=0.66; MxA50=1.65)
              if(Scen=='S1')
              {
                Mat=Report$Growth_Parameters%>%filter(Sex==1)%>%pull(Mat1)
                Linf=Report$Growth_Parameters%>%filter(Sex==1)%>%pull(Linf)
                Ratio_L50_Linf=Mat/Linf
                K=Report$Growth_Parameters%>%filter(Sex==1)%>%pull(K)
                AgE.kls=names(Report$Natural_Mortality)
                AgE.kls=AgE.kls[grepl("^[0-9]+$", AgE.kls)]
                Nat.M=Report$Natural_Mortality%>%filter(Sex==1)
                Nat.M=Nat.M[,match(AgE.kls,names(Nat.M))]
                Nat.M=unlist(Nat.M)
                Ratio_M_k=mean(Nat.M)/K
                Prod_M_A50=mean(Nat.M)*Life.history$Age.50.mat[1]
                #Prod_M_A50=exp(1.44-0.982*log(Life.history$Max.age.F))*Life.history$Age.50.mat[1]
                
                write.csv(data.frame(Type=c('L50/Linf','M/k','MxA50'),
                                     Value=c(Ratio_L50_Linf,Ratio_M_k,Prod_M_A50)),
                          paste0(this.wd1,"/Life.history.invariants.csv"),row.names = F)
                
              }
              
              #c.11. Check if pups per female is not biologically excessively large 
              if(Scen=='S1')
              {
                Report$recruit%>%
                  mutate(Pups.per.female=exp_recr/mature_num)%>%
                  ggplot(aes(Yr,Pups.per.female))+
                  geom_point()+
                  geom_line()+theme_PA()+
                  labs(title = paste('Biological annual pups per female=',Life.history$Fecundity/Life.history$Breed.cycle),
                       caption = 'Pups per female should be <= average annual fecundity (SS adds extra put mortality)')
                ggsave(paste0(this.wd1,"/Pups per female.tiff"),width=6,height=6,compression = "lzw") 
              }
              
              #c.12. Check cryptic biomass (The ratio should be less than 25%)
              if(Scen=='S1')
              {
                Cryptic=fn.cryptic(yr=as.numeric(substr(Last.yr.ktch,1,4)))
                ggarrange(plotlist = Cryptic$p.criptic,ncol=1)
                ggsave(paste0(this.wd1,"/Cryptic biomass.tiff"),width=5,height=6,compression = "lzw")
                ggarrange(plotlist = Cryptic$p.criptic1,ncol=1)
                ggsave(paste0(this.wd1,"/Cryptic biomass1.tiff"),width=5,height=6,compression = "lzw")
              }
              
              #c.13. Check what sels by fleet SS is applying
              if(Scen=='S1')
              {
                pp=fn.check.SS.sel.used(d=Report$sizeselex%>%filter(Factor=='Lsel')%>%dplyr::select(-c(Factor,Label)),
                                        check.fleet=3)
                pp$p1
                ggsave(paste0(this.wd1,"/SS_selectivity_used.tiff"),width=7,height=6,compression = "lzw")
                
                pp$p2
                ggsave(paste0(this.wd1,"/SS_selectivity_used_fleet 3.tiff"),width=7,height=6,compression = "lzw")
              }
              
              #c.14 get fleets
              Scen.flits=Report$FleetNames
              
              #c.15 return stuff
              return(list(dummy.store.rel.biom=Out.rel.biom,
                          dummy.store.probs.rel.biom=Out.probs.rel.biom,
                          dummy.store.probs.B.Bmsy=Out.probs.B.Bmsy,
                          dummy.store.probs.f.series=Out.probs.f.series,
                          dummy.store.f.series=Out.f.series,
                          dummy.store.B.Bmsy=Out.B.Bmsy,
                          dummy.store.F.Fmsy=Out.F.Fmsy,
                          dummy.store.Kobe.probs=Out.Kobe.probs,  
                          dummy.store.estimates=Out.estimates,
                          dummy.store.likelihoods=Out.likelihoods,
                          dummy.store.quantities=Out.quantities,
                          dummy.store.flits=Scen.flits))
            } # end 'file.exists(paste0(this.wd1,'/Report.sso'))'
            
          } # end 'file.exists(paste0(this.wd1,'/control.ctl'))'
          
        } #end i loop
      stopCluster(cl)
      names(out.species)=Specs.scens
      
      #3. Extract species and scenario outputs
      #here, pull out species from out.species...which should be a list of species
      out.species1<- vector('list',length(List.sp))
      names(out.species1) = names(List.sp)
      for(v in 1:length(out.species1))
      {
        Neim=names(out.species1)[v]
        Life.history=List.sp[[Neim]]
        ID.sp=grep(Neim,names(out.species))
        
        #Get scenarios info
        Dis.SceNs=Species.scenarios%>%filter(Species==Neim)%>%pull(Scenario)
        Out.Scens=List.sp[[Neim]]$Sens.test$SS%>%
          filter(Scenario%in%Dis.SceNs)%>%
          mutate(Species=capitalize(Neim))%>%  
          rename(h.mean=Steepness,
                 h.sd=Steepness.sd)%>%
          mutate(Mmean=round(Mmean,3),
                 h.mean=round(h.mean,3),
                 h.sd=round(h.sd,3))%>%
          dplyr::select(-c(Model,use_F_ballpark))%>%
          relocate(Species,Scenario,Mmean,F_ballpark,h.mean,h.sd)
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
                 b_FL.to.TL=Life.history$b_FL.to.TL)
        
        #Get relevant model outputs
        out.species1[[v]]$dummy.store.rel.biom=do.call(rbind,fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.rel.biom'))
        out.species1[[v]]$dummy.store.probs.rel.biom=fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.probs.rel.biom')
        out.species1[[v]]$dummy.store.probs.B.Bmsy=fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.probs.B.Bmsy')
        out.species1[[v]]$dummy.store.probs.f.series=fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.probs.f.series')
        out.species1[[v]]$dummy.store.f.series=do.call(rbind,fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.f.series'))
        out.species1[[v]]$dummy.store.B.Bmsy=do.call(rbind,fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.B.Bmsy'))
        out.species1[[v]]$dummy.store.F.Fmsy=do.call(rbind,fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.F.Fmsy'))
        out.species1[[v]]$dummy.store.Kobe.probs=fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.Kobe.probs')
        out.species1[[v]]$dummy.store.estimates=do.call(rbind,fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.estimates'))
        out.species1[[v]]$dummy.store.likelihoods=do.call(rbind,fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.likelihoods'))
        out.species1[[v]]$dummy.store.quantities=do.call(rbind,fn.get.and.name(LISTA=out.species[ID.sp],x='dummy.store.quantities'))
        out.species1[[v]]$dummy.store.sens.table=Out.Scens%>%
          mutate(Fleets=sapply(fn.get.and.name(LISTA=out.species[ID.sp],
                                               x='dummy.store.flits'),paste, collapse = ", "))
      }
      out.species=out.species1
      rm(out.species1)
      
      #4. Store relevant quantities
      Age.based[[w]]$sens.table=fn.get.and.name(LISTA=out.species,x='dummy.store.sens.table')
      Age.based[[w]]$estimates=fn.get.and.name(LISTA=out.species,x='dummy.store.estimates')
      Age.based[[w]]$rel.biom=fn.get.and.name(LISTA=out.species,x='dummy.store.rel.biom')
      Age.based[[w]]$probs.rel.biom=fn.get.and.name(LISTA=out.species,x='dummy.store.probs.rel.biom')
      Age.based[[w]]$probs.B.Bmsy=fn.get.and.name(LISTA=out.species,x='dummy.store.probs.B.Bmsy')
      Age.based[[w]]$probs.f.series=fn.get.and.name(LISTA=out.species,x='dummy.store.probs.f.series') 
      Age.based[[w]]$f.series=fn.get.and.name(LISTA=out.species,x='dummy.store.f.series')
      Age.based[[w]]$B.Bmsy=fn.get.and.name(LISTA=out.species,x='dummy.store.B.Bmsy')
      Age.based[[w]]$F.Fmsy=fn.get.and.name(LISTA=out.species,x='dummy.store.F.Fmsy')
      Age.based[[w]]$Kobe.probs=fn.get.and.name(LISTA=out.species,x='dummy.store.Kobe.probs')
      Age.based[[w]]$likelihoods=fn.get.and.name(LISTA=out.species,x='dummy.store.likelihoods')
      Age.based[[w]]$quantities=fn.get.and.name(LISTA=out.species,x='dummy.store.quantities')
      
      #5. Plot overall figure with estimates and scenarios  
      if(do.all.sensitivity.tests)
      {
        for(i in 1:N.sp)
        {
          print(paste("SS3 rel biom figure and estimates ____",Keep.species[i]))
          this.wd=paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated",sep='')
          if(!is.null(Age.based[[w]]$rel.biom[[i]]))
          {
            yrS=Age.based[[w]]$rel.biom[[i]]%>%filter(Scenario=='S1')%>%pull(year)
            if(length(yrS)<50) delta=10 else
              delta=5
            xmin=min(yrS)+delta
            xmax=max(yrS)-delta
            TAB=Age.based[[w]]$estimates[[i]]%>%
              filter(Par=='SR_LN(R0)')%>%
              dplyr::select(Par,Value,Init,Parm_StDev,Status,Gradient)%>%
              `rownames<-`( NULL )
            TAB.scen=Age.based[[w]]$sens.table[[i]]
            this.col=rep('no',ncol(TAB.scen))
            for(u in 1:ncol(TAB.scen)) if(length(unique(TAB.scen[,u]))>1) this.col[u]='Yes'
            TAB.scen=TAB.scen[,which(this.col=='Yes')]
            Tab.combined=cbind(TAB.scen,TAB)
            Fnt.size=10
            yMIN=0; yMAX=0.3
            if(ncol(Tab.combined)>10) Fnt.size=6
            if(nrow(Tab.combined)>6)
            {
              yMIN=0.1; yMAX=0.4
            }
            Levls=(unique(Age.based[[w]]$rel.biom[[i]]$Scenario))
            cols=colfunc1(length(Levls))
            names(cols)=Levls
            p=Age.based[[w]]$rel.biom[[i]]%>% 
              mutate(Scenario=factor(Scenario,levels=Levls))%>%
              ggplot(aes(year,median,color=Scenario))+
              annotation_custom(tableGrob(Tab.combined,rows=NULL, theme = ttheme_default(base_size = Fnt.size)),
                                xmin=xmin, xmax=xmax, ymin=yMIN, ymax=yMAX)+
              geom_line(linewidth=2,alpha=0.6,na.rm = TRUE)+
              geom_line(aes(year,upper.95),linetype=2,alpha=0.6,na.rm = TRUE)+
              geom_line(aes(year,lower.95),linetype=2,alpha=0.6,na.rm = TRUE)+
              ggtitle(Keep.species[i])+ylim(0,1)+
              theme_PA(leg.siz=9)+
              theme(legend.position = 'bottom',
                    legend.title = element_blank())+
              ylab('Relative biomass')+xlab('Year')+
              guides(colour = guide_legend(nrow = 1))+
              scale_color_manual(values = cols)
            print(p)
            ggsave(paste(this.wd,"/Rel.biomass&estimates.tiff",sep=''),width=10,height=7,compression = "lzw")  
            
            if(nrow(Tab.combined)>6)
            {
              yMIN=0.2; yMAX=1
              if(nrow(Tab.combined)>10) yMAX=2
            }
            p=Age.based[[w]]$B.Bmsy[[i]]%>%
              mutate(Scenario=factor(Scenario,levels=Levls))%>%
              ggplot(aes(year,median,color=Scenario))+
              annotation_custom(tableGrob(Tab.combined,rows=NULL, theme = ttheme_default(base_size = Fnt.size)),
                                xmin=xmin+2, xmax=xmax, ymin=yMIN, ymax=yMAX)+
              geom_line(linewidth=2,alpha=0.6,na.rm = TRUE)+
              geom_line(aes(year,upper.95),linetype=2,alpha=0.6,na.rm = TRUE)+
              geom_line(aes(year,lower.95),linetype=2,alpha=0.6,na.rm = TRUE)+
              ggtitle(Keep.species[i])+ylim(0,NA)+
              theme_PA(leg.siz=9)+
              theme(legend.position = 'bottom',
                    legend.title = element_blank())+
              ylab('B/Bmsy')+xlab('Year')+
              guides(colour = guide_legend(nrow = 1))+
              scale_color_manual(values = cols)
            print(p)
            ggsave(paste(this.wd,"/B.over.BMSY&estimates.tiff",sep=''),width=10,height=7,compression = "lzw")  
            
            rm(TAB)
            
          }
        }
      }
      
      #6. Compare likelihoods from different scenarios
      if(do.all.sensitivity.tests & Compare.sensitivity.test.likelihoods)
      {
        for(i in 1:N.sp)
        {
          this.wd=paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated",sep='')
          Scens=List.sp[[i]]$Sens.test$SS  
          dis.files=paste(this.wd,Scens$Scenario,sep='/')
          Report.list=vector('list',length(dis.files))
          names(Report.list)=Scens$Scenario
          for(k in 1:length(Report.list)) Report.list[[k]]=SS_output(dis.files[[k]],covar=FALSE,forecast=TRUE,readwt=F)
          p=fn.SS3.compare.likelihoods(Report.list=Report.list)
          print(p)
          ggsave(paste(this.wd,"/Compare likelihoods.tiff",sep=''),width=10,height=7,compression = "lzw")
        }
      }
      
      rm(out.species)
    }
    
    #Sequential processing. Superseded
    if(!do.parallel.SS)
    {
      fn.source1('Run.SS.sequentially.R')   
    }
  }
}

toc(log = TRUE, quiet = TRUE)
computation.time <- tic.log(format = TRUE)
tic.clearlog()
Arg.sens=Arg
if(!do.all.sensitivity.tests.uncertainty)
{
  Arg.sens='-nohess'
  how.many.scens='only S1 scenario,'
}
if(do.all.sensitivity.tests.uncertainty) how.many.scens='all scenarios,'

run.with.hess=ifelse(Arg.sens=='','with Hessian estimation',
                     ifelse(Arg.sens=="-nohess","without Hessian estimation",
                            NA))
send.email(TO=Send.email.to,
           CC='',
           Subject=paste("SS3 models for ",N.sp,'species,',how.many.scens,run.with.hess,"finished running at",Sys.time()),
           Body= paste("Computation",computation.time),  
           Attachment=NULL) 

#---Do SS3 diagnostics -------------------------------------------------
#notes: 
        # runs test only displays for series with >1 year
        # Hindcasting cross-validation only available for species with abundance series
        # takes ~30 mins per species (in parallel)
if(do.SS3.diagnostics)
{
  options(warn=-1)
  tic("timer")
  set.seed(1234)
  pb <- txtProgressBar(min = 0,max =N.sp, style = 3,width = 50)
  SCEN="S1"
  for(l in 1:N.sp)
  {
    Neim=Keep.species[l]
    setTxtProgressBar(pb, l)
    this.wd=paste(HandL.out,capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
    this.wd1=paste(this.wd,SCEN,sep='/')
    if(file.exists(this.wd1))
    {
      MLE=read.admbFit(paste(this.wd1,'ss3',sep='/'))
      #MLE=read.admbFit(paste(this.wd1,'ss',sep='/')) #if old SS version
      dat_temp <- SS_readdat(file.path(this.wd1, "data.dat"), verbose = FALSE)
      ctl_temp <- SS_readctl(file = file.path(this.wd1, "control.ctl"), datlist = dat_temp, verbose = FALSE)
      replist <- SS_output(dir = this.wd1,verbose = FALSE,printstats = FALSE,covar = TRUE)
      derived_quants <- replist$derived_quants
      
      #get sequence for profile likelihoods 
      Estim.LnRo=MLE$est[grep("SRparm",MLE$names)]
        #Ro
      std.LnRo=MLE$std[grep("SRparm",MLE$names)]
      R0.range=fn.like.range(Par.mle=Estim.LnRo,
                             min.par=2,
                             Par.SE=std.LnRo,
                             up=delta.likelihood.profiles,
                             low=delta.likelihood.profiles,
                             ln.out=Number.of.likelihood.profiles,
                             seq.approach=Approach.like.prof)
        #h
      Input.h=List.sp[[Neim]]$Sens.test$SS%>%filter(Scenario==SCEN)%>%pull(Steepness)
      h.range=fn.like.range(Par.mle=Input.h,
                            min.par=Min.h.shark,
                            Par.SE=0.06,
                            up=0.4,
                            low=0.4,
                            ln.out=Number.of.likelihood.profiles.h,
                            seq.approach='min.plus')
        #Ma
      M.range=fn.like.range(Par.mle=unlist(ctl_temp$natM[1,]),
                            min.par=0.01,
                            Par.SE=3*List.sp[[Neim]]$Sens.test$SS%>%filter(Scenario==SCEN)%>%pull(Msd),
                            up=0.03,
                            low=0.03,
                            ln.out=Number.of.likelihood.profiles.h,
                            seq.approach='SE')
  
        #Current depletion
      Depl.range=fn.like.range(Par.mle=derived_quants[paste0("Bratio_",replist$endyr),'Value'],
                               min.par=0.1,
                               Par.SE=derived_quants[paste0("Bratio_",replist$endyr),'StdDev'],
                               up=0.4,
                               low=0.4,
                               ln.out=Number.of.likelihood.profiles,
                               seq.approach='SE')
        
        #Current spawning biomass
      CurSB.range=fn.like.range(Par.mle=derived_quants[paste0("SSB_",replist$endyr),'Value'],
                                min.par=0.1,
                                Par.SE=derived_quants[paste0("SSB_",replist$endyr),'StdDev'],
                                up=0.4,
                                low=0.4,
                                ln.out=Number.of.likelihood.profiles,
                                seq.approach='SE')
      Par_VAR_PROF=c("R0","h","M","Depl","CurSB")
      
        #Linf
      Linf.VEC.F=Linf.VEC.M=NA
      if(any(grepl("MGparm",MLE$names)))  
      {
        Linf.est=MLE$est[grep("MGparm",MLE$names)]
        Linf.std=MLE$std[grep("MGparm",MLE$names)]
        
        Linf.VEC.F=fn.like.range(Par.mle=Linf.est[1],
                                 min.par=0,
                                 Par.SE=5*Linf.std[1],
                                 up=0,
                                 low=0,
                                 ln.out=Number.of.likelihood.profiles,
                                 seq.approach='SE')
        Linf.VEC.M=fn.like.range(Par.mle=Linf.est[2],
                                 min.par=0,
                                 Par.SE=5*Linf.std[2],
                                 up=0,
                                 low=0,
                                 ln.out=Number.of.likelihood.profiles,
                                 seq.approach='SE')
        Par_VAR_PROF=c(Par_VAR_PROF,'L_at_Amax_Fem','L_at_Amax_Mal')
      }
      
      
        #Execute wrapping function to generate diagnostics
      fn.fit.diag_SS3(WD=this.wd1,
                      do.like.prof=TRUE,
                      disfiles=c("control.ss_new", "data.dat","forecast.ss","starter.ss","Report.sso"),
                      R0.vec=R0.range,
                      h.vec=h.range,
                      M.vec=M.range,
                      depl.vec=Depl.range,
                      curSB.vec=CurSB.range,
                      Linf.vec.F=Linf.VEC.F,
                      Linf.vec.M=Linf.VEC.M,
                      exe_path=Where.exe,
                      start.retro=Retro_start,
                      end.retro=Retro_end,
                      do.retros=TRUE,
                      do.jitter=TRUE,
                      numjitter=Number.of.jitters,
                      fracjitter=jitters.frction,
                      outLength.Cross.Val=FALSE,
                      run.in.parallel=TRUE,
                      flush.files=TRUE,
                      COVAR=TRUE,
                      h.input=Input.h,
                      Par_var_profile=Par_VAR_PROF)
      rm(MLE,this.wd1,R0.range,Estim.LnRo)
    }
    cat(paste(" completed --- Species:", Neim))
  }
  close(pb)
  toc(log = TRUE, quiet = TRUE)
  computation.time <- tic.log(format = TRUE)
  tic.clearlog()
  send.email(TO=Send.email.to,
             CC='',
             Subject=paste("SS3 model diagnostics for",N.sp,
                           "species and Scenario 1 finished running at",Sys.time()),
             Body= paste("Computation time was",computation.time),  
             Attachment=NULL)
  options(warn=0)
}


#---Compare different scenarios for each species -------------------------------------------------
if(SS3.run=='final')
{
  for(i in 1:N.sp)
  {
    Neim=Keep.species[i]
    this.wd=paste(HandL.out,capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
    if(file.exists(this.wd))
    {
      print(paste("Compare scenarios for ------------",Neim))
      Scens=List.sp[[i]]$Sens.test$SS
      LIST=vector('list',nrow(Scens))
      names(LIST)=Scens$Scenario
      for(s in 1:nrow(Scens))
      {
        this.wd1=paste(this.wd,Scens$Scenario[s],sep='/')
        LIST[[s]]=SS_output(this.wd1)
      }
      # create list summarizing model results
      mod.sum <- SSsummarize(LIST,verbose = FALSE)
      # plot comparisons
      SUPLT=1:16
      SUPLT=subset(SUPLT,!SUPLT%in%c(1,3,5,7,9))
      this.wd.compare=paste(this.wd,'Scenario comparison',sep='/')
      if(!dir.exists(this.wd.compare))dir.create(this.wd.compare) 
      SSplotComparisons(mod.sum, legendlabels = names(LIST),pheight=4.5,plot = FALSE,png=TRUE,
                        plotdir=this.wd.compare,subplots=SUPLT,legendloc='bottomleft')
    }
  }
}


#---Generate outputs -------------------------------------------------
#25.1 Table of scenarios
for(l in 1: length(Lista.sp.outputs))    
{
  for(w in 1:length(Age.based))
  {
    dummy=Age.based[[w]]$sens.table
    dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
    dummy=compact(dummy)
    d.tbl=bind_rows(dummy)%>%relocate(Species)
    if(!evaluate.07.08.cpue) d.tbl=d.tbl%>%dplyr::select(-Daily.cpues)
    write.csv(d.tbl,
              paste(HandL.out.RAR,paste('Table 10. Age.based_',names(Age.based)[w],'_scenarios_',
                                   names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
              row.names = F)
  }
}

#25.2 Table of parameter estimates, likelihoods and quantities of interest
for(l in 1: length(Lista.sp.outputs))
{
  for(w in 1:length(Age.based))
  {
    #Parameter estimates
    dummy=Age.based[[w]]$estimates
    dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
    dummy=compact(dummy)
    if(length(dummy)>0)
    {
      dummy=bind_rows(dummy)%>%
        rownames_to_column(var = "Species")%>%
        mutate(Species=capitalize(str_extract(Species, "[^.]+")))%>%
        relocate(Species)
      write.csv(dummy,paste(HandL.out.RAR,paste('Table 11. Age.based_',
                                           names(Age.based)[w],'_estimates_',
                                           names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
                row.names = F)
    }
    
    #Likelihoods
    dummy=Age.based[[w]]$likelihoods
    dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
    dummy=compact(dummy)
    if(length(dummy)>0)
    {
      dummy=bind_rows(dummy)%>%
        rownames_to_column(var = "Species")%>%
        mutate(Species=capitalize(str_extract(Species, "[^.]+")))%>%
        relocate(Species)
      write.csv(dummy,paste(HandL.out.RAR,paste('Table 11. Age.based_',
                                           names(Age.based)[w],'_likelihoods_',
                                           names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
                row.names = F)
    }
    
    #Quantities
    dummy=Age.based[[w]]$quantities
    dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
    dummy=compact(dummy)
    if(length(dummy)>0)
    {
      write.csv(bind_rows(dummy),paste(HandL.out.RAR,paste('Table 11. Age.based_',
                                           names(Age.based)[w],'_quantities_',
                                           names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
                row.names = F)
    }

  }
}


#25.3. Time series   
YLAB=XLAB=''

#25.3.1 Display all scenarios for each species
  #25.3.1.1 Relative biomass (i.e. Depletion)
Ref.points=vector('list',N.sp)
names(Ref.points)=Keep.species
for(i in 1:N.sp)
{
  if(!is.null(Age.based$SS$rel.biom[[i]]))
  {
    print(paste("SS3 --- Relative biomass plot -----",Keep.species[i])) 
    a=fn.plot.timeseries(d=Age.based,    
                         sp=Keep.species[i],
                         Type='Depletion',
                         YLAB='Relative biomass')
    if(!is.null(a))
    {
      #export graph
      ggsave(paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_time_series_relative_biomass.tiff",sep=''),
             width = 6,height = 10,compression = "lzw")
      
      #export current depletion probabilities
      write.csv(a$store.probs%>%
                  spread(Model,Probability)%>%
                  mutate(Species=Keep.species[i],
                         Scenario=factor(Scenario,levels=unique(a$store.probs$Scenario)),
                         Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                  arrange(Scenario,Range),
                paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_current_depletion_probabilities.csv",sep=''),
                row.names = F)
      
      #export current depletion 
      xx=Age.based$SS$rel.biom[[i]]
      write.csv(xx,
                paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_depletion.csv",sep=''),
                row.names = F)
      
      Ref.points[[i]]=a$Ref.points$SS
    }
  }

}

  #25.3.1.2 Fishing mortality
if(do.F.series)
{
  for(i in 1:N.sp)
  {
    if(!is.null(Age.based$SS$f.series[[i]]))
    {
      print(paste("SS3 --- Fishing mortality plot -----",Keep.species[i]))
      a=fn.plot.timeseries(d=Age.based,
                           sp=Keep.species[i],
                           Type='F.series',
                           YLAB=expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")))
      if(!is.null(a))
      {
        ggsave(paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_time_series_fishing_mortality.tiff",sep=''),
               width = 8,height = 10,compression = "lzw")
      }
    }
  }
}

  #25.3.1.3 B over Bmsy
if(do.B.over.Bmsy.series)
{
  for(i in 1:N.sp)
  {
    if(!is.null(Age.based$SS$B.Bmsy[[i]]))
    {
      print(paste("SS3 --- B over Bmsy plot -----",Keep.species[i]))
      a=fn.plot.timeseries(d=Age.based,
                           sp=Keep.species[i],
                           Type='B.Bmsy',
                           YLAB='B/Bmsy')
      if(!is.null(a))
      {
        ggsave(paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_time_series_B_Bmsy.tiff",sep=''),
               width = 8,height = 10,compression = "lzw")
      }
    }
  }
}

  #25.3.1.4 F over Fmsy
if(do.F.over.Fmsy.series)
{
  for(i in 1:N.sp)
  {
    if(!is.null(Age.based$SS$F.Fmsy[[i]]))
    {
      print(paste("SS3 --- F over Fmsy plot -----",Keep.species[i]))
      a=fn.plot.timeseries(d=Age.based,
                           sp=Keep.species[i],
                           Type='F.Fmsy',
                           YLAB='F/Fmsy')
      if(!is.null(a))
      {
        ggsave(paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_time_series_F_Fmsy.tiff",sep=''),
               width = 8,height = 10,compression = "lzw")
      }
    }
  }
}


#25.3.2 Display Scenario 1 for combined species 
  #25.3.2.1 Relative biomass (i.e. Depletion) 
    #figure
for(l in 1:length(Lista.sp.outputs))
{
  print(paste("RAR --- Depletion_SS3 plot S1 -----",names(Lista.sp.outputs)[l],"----- single plot"))
  
  if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
  a=fn.plot.timeseries_combined(this.sp=Lista.sp.outputs[[l]],
                                d=Age.based$SS,
                                YLAB="Relative biomass",
                                Type="Depletion",
                                InnerMargin=InMar,
                                RefPoint=Ref.points,
                                Kach=Age.based$SS$rel.biom)
  WIDt=10
  if(length(compact(Age.based$SS$sens.table))<=3) WIDt=7
  if(!is.null(a))ggsave(paste(HandL.out.RAR,'/Relative.biomass_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
                        width = WIDt,height = 10,compression = "lzw")
}
    #table 
for(l in 1:length(Lista.sp.outputs))
{
  dummy.mod=vector('list',length(Age.based))
  for(m in 1:length(Age.based))
  {
    str.prob=Age.based[[m]]$probs.rel.biom
    str.prob=str.prob[match(Lista.sp.outputs[[l]],names(str.prob))]
    str.prob=compact(str.prob)
    if(length(str.prob)>0)
    {
      dummy=vector('list',length =length(str.prob))
      for(d in 1:length(dummy))
      {
        dummy[[d]]=str.prob[[d]][[1]]$probs%>%
          mutate(Species=capitalize(names(str.prob)[d]))
        
        if('probs.future'%in%names(str.prob[[d]][[1]]))
        {
          dummy[[d]]=rbind(dummy[[d]],
          str.prob[[d]][[1]]$probs.future%>%
            mutate(Species=capitalize(names(str.prob)[d])))
        }
      }
      dummy.mod[[m]]=bind_rows(dummy)%>%
        mutate(Model=names(Age.based)[m])
    }
  }
  dummy.mod=compact(dummy.mod)
  if(length(dummy.mod)>0)
  {
    write.csv(bind_rows(dummy.mod)%>%
                mutate(Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                spread(Species,Probability)%>%
                arrange(Range),
              paste(HandL.out.RAR,'/Table 12. Age.based_SS_current.depletion_',names(Lista.sp.outputs)[l],'.csv',sep=''),
              row.names=F)
    rm(dummy.mod)
    
  }
}
  #25.3.2.2 B over Bmsy  
if(do.B.over.Bmsy.series)
{
  #figure
  for(l in 1:length(Lista.sp.outputs))
  {
    print(paste("RAR --- B.over.Bmsy_SS3 plot S1 -----",names(Lista.sp.outputs)[l],"----- single plot")) 
    if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
    a=fn.plot.timeseries_combined(this.sp=Lista.sp.outputs[[l]],
                                  d=Age.based$SS,
                                  YLAB="B/Bmsy",
                                  Type="B.Bmsy",
                                  InnerMargin=InMar,
                                  RefPoint=NULL,
                                  Kach=Age.based$SS$rel.biom)
    WIDt=10
    if(length(compact(Age.based$SS$sens.table))<=3) WIDt=7
    if(!is.null(a))ggsave(paste(HandL.out.RAR,'/B.over.Bmsy_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
                          width = WIDt,height = 10,compression = "lzw")
  }
  #table 
  for(l in 1:length(Lista.sp.outputs))
  {
    dummy.mod=vector('list',length(Age.based))
    for(m in 1:length(Age.based))
    {
      str.prob=Age.based[[m]]$probs.B.Bmsy
      str.prob=str.prob[match(Lista.sp.outputs[[l]],names(str.prob))]
      str.prob=compact(str.prob)
      if(length(str.prob)>0)
      {
        dummy=vector('list',length =length(str.prob))
        for(d in 1:length(dummy))
        {
          dummy[[d]]=str.prob[[d]][[1]]$probs%>%
            mutate(Species=capitalize(names(str.prob)[d]))
          
          if('probs.future'%in%names(str.prob[[d]][[1]]))
          {
            dummy[[d]]=rbind(dummy[[d]],
                             str.prob[[d]][[1]]$probs.future%>%
                               mutate(Species=capitalize(names(str.prob)[d])))
          }
        }
        dummy.mod[[m]]=bind_rows(dummy)%>%
          mutate(Model=names(Age.based)[m])
      }
    }
    dummy.mod=compact(dummy.mod)
    if(length(dummy.mod)>0)
    {
      write.csv(bind_rows(dummy.mod)%>%
                  mutate(Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                  spread(Species,Probability)%>%
                  arrange(Range),
                paste(HandL.out.RAR,'/Table 12. Age.based_SS_Current.B.over.Bmsy_',names(Lista.sp.outputs)[l],'.csv',sep=''),
                row.names=F)
      rm(dummy.mod)
      
    }
  }
  
}
  #25.3.2.3 F over Fmsy  
if(do.F.over.Fmsy.series)
{
  #figure
  for(l in 1:length(Lista.sp.outputs))
  {
    print(paste("RAR --- F.over.Fmsy_SS3 plot S1 -----",names(Lista.sp.outputs)[l],"----- single plot")) 
    if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
    a=fn.plot.timeseries_combined(this.sp=Lista.sp.outputs[[l]],
                                  d=Age.based$SS,
                                  YLAB="F/Fmsy",
                                  Type="F.Fmsy",
                                  InnerMargin=InMar,
                                  RefPoint=NULL,
                                  Kach=Age.based$SS$rel.biom)
    WIDt=10
    if(length(compact(Age.based$SS$sens.table))<=3) WIDt=7
    if(!is.null(a))ggsave(paste(HandL.out.RAR,'/F.over.Fmsy_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
                          width = WIDt,height = 10,compression = "lzw")
  }
}
  #25.3.2.4 F series  
if(do.F.series)
{
  #figure
  for(l in 1:length(Lista.sp.outputs))
  {
    print(paste("RAR --- F_SS3 plot S1 -----",names(Lista.sp.outputs)[l],"----- single plot")) 
    if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
    a=fn.plot.timeseries_combined(this.sp=Lista.sp.outputs[[l]],
                                  d=Age.based$SS,
                                  YLAB=expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),
                                  Type="F.series",
                                  InnerMargin=InMar,
                                  RefPoint=NULL,
                                  Kach=Age.based$SS$rel.biom)
    WIDt=10
    if(length(compact(Age.based$SS$sens.table))<=3) WIDt=7
    if(!is.null(a))ggsave(paste(HandL.out.RAR,'/F.series_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
                          width = WIDt,height = 10,compression = "lzw")
  }
  
  #table 
  for(l in 1:length(Lista.sp.outputs))
  {
    dummy.mod=vector('list',length(Age.based))
    for(m in 1:length(Age.based))
    {
      str.prob=Age.based[[m]]$probs.f.series
      str.prob=str.prob[match(Lista.sp.outputs[[l]],names(str.prob))]
      str.prob=compact(str.prob)
      if(length(str.prob)>0)
      {
        dummy=vector('list',length =length(str.prob))
        for(d in 1:length(dummy))
        {
          dummy[[d]]=str.prob[[d]][[1]]$probs%>%
            mutate(Species=capitalize(names(str.prob)[d]))
          
          if('probs.future'%in%names(str.prob[[d]][[1]]))
          {
            dummy[[d]]=rbind(dummy[[d]],
                             str.prob[[d]][[1]]$probs.future%>%
                               mutate(Species=capitalize(names(str.prob)[d])))
          }
        }
        dummy.mod[[m]]=bind_rows(dummy)%>%
          mutate(Model=names(Age.based)[m])
      }
    }
    dummy.mod=compact(dummy.mod)
    if(length(dummy.mod)>0)
    {
      write.csv(bind_rows(dummy.mod)%>%
                  mutate(Range=factor(Range,levels=c(">lim","lim.thr","thr.tar","<tar")))%>%
                  spread(Species,Probability)%>%
                  arrange(Range),
                paste(HandL.out.RAR,'/Table 12. Age.based_SS_Current.f_',names(Lista.sp.outputs)[l],'.csv',sep=''),
                row.names=F)
      rm(dummy.mod)
      
    }
  }
}

#25.3.3 Display sensitivity tests for combined species
disspisis=names(Age.based$SS$rel.biom)[!sapply(Age.based$SS$rel.biom,is.null)] 
for(l in 1:length(Lista.sp.outputs))
{
  print(paste("SS3 --- Relative biomass plot by Scenario -----",names(Lista.sp.outputs)[l],"----- single plot"))
  if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
  a=fn.plot.timeseries_combined_sensitivity(this.sp=Lista.sp.outputs[[l]], 
                                            d=Age.based$SS,
                                            InnerMargin=InMar,
                                            RefPoint=Ref.points,
                                            Kach=Age.based$SS$rel.biom)
  HEIT=8
  WIDt=8
  if(length(which(Lista.sp.outputs[[l]]%in%disspisis))<=3) WIDt=6
  if(!is.null(a))ggsave(paste(HandL.out.RAR,'/Relative.biomass_SS3 integrated_',names(Lista.sp.outputs)[l],'_sensitivity.tiff',sep=''),
                        width = WIDt,height = HEIT,compression = "lzw")
}

#25.4. Compare goodness of fit for different Scenarios  
#note: compare DIC, RMSE & SDNR 
for(l in 1:length(Lista.sp.outputs))
{
  Nms=names(compact(Age.based$SS$estimates))
  this.sp=Lista.sp.outputs[[l]]
  this.sp=this.sp[which(this.sp%in%Nms)]
  print(paste("SS3 --- Compare Scenarios goodness of fit -----",names(Lista.sp.outputs)[l])) 
  if(length(this.sp)>0)
  {
    stor.plt=vector('list',length(this.sp))
    for(q in 1:length(this.sp))
    {
      i=match(this.sp[q],Keep.species)
      this.wd=paste(HandL.out,capitalize(List.sp[[i]]$Name),"/",AssessYr,"/SS3 integrated",sep='')
      if(file.exists(this.wd))
      {
        Scens=List.sp[[i]]$Sens.test$SS$Scenario 
        a=vector('list',length(Scens))
        for(s in 1:length(Scens))
        {
          a[[s]]=read.csv(paste(paste(this.wd,"/",Scens[s],sep=''),
                                paste0('GoodnessFit_',List.sp[[i]]$Name,'.csv'),sep='/'))%>%
            dplyr::select(-X)%>%mutate(Scenario=Scens[s])
        }
        stor.plt[[q]]=bind_rows(a)%>%
          filter(Stastistic%in%c('RMSE','AIC'))%>%
          ggplot(aes(Scenario,Value))+
          geom_col(fill='brown')+
          facet_wrap(~Stastistic,scales='free',ncol=2)+
          theme_PA()+ylab('')+xlab('')+
          ggtitle(capitalize(List.sp[[i]]$Name))
      }
    }
    figure <- ggarrange(plotlist=stor.plt,ncol=1,nrow=length(this.sp),common.legend = TRUE)
    figure=annotate_figure(figure,bottom = text_grob('Scenario', size=22),left = text_grob('Value',rot = 90, size=22))
    ggsave(paste(HandL.out.RAR,'/GoodnessFit per scenario_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = 5,height = 8,compression = "lzw")
  }
}


#25.5. Kobe plots (Scenario 1) 
  #25.5.1 by species 
store.kobes=vector('list',N.sp)
names(store.kobes)=Keep.species
for(i in 1:N.sp)
{
  if(!is.null(Age.based$SS$estimates[[i]]))
  {
    print(paste("SS3 --- Kobe plot -----",Keep.species[i]))
    store.kobes[[i]]=fn.get.Kobe.plot_appendix(d=Age.based,
                                               sp=Keep.species[i],
                                               do.probs=TRUE)
    ggsave(paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_Kobe_plot.tiff",sep=''),
           width = 10,height = 10, dpi = 300,compression = "lzw")
  }
  
}
store.kobes=compact(store.kobes)

  #25.5.2 Display combined species     
for(l in 1:length(Lista.sp.outputs))
{
  print(paste("SS3 --- Kobe plot -----",names(Lista.sp.outputs)[l]))
  Nms=names(compact(Age.based$SS$estimates))
  this.sp=Lista.sp.outputs[[l]]
  this.sp=this.sp[which(this.sp%in%Nms)]
  if(length(this.sp)>0)
  {
    DIMS=n2mfrow(length(this.sp))
    NKOL=DIMS[2]
    NRW=DIMS[1]
    if(NKOL%in%3:4) WIZ=14
    if(NKOL==2) WIZ=11
    if(NKOL==1) WIZ=9
    fn.get.Kobe.plot(this.sp,
                     d=Age.based$SS,
                     NKOL,
                     NRW,
                     do.probs=TRUE)
    ggsave(paste(HandL.out.RAR,'/Kobe_plot_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = WIZ,height = 12,compression = "lzw")
    
  }
}

#25.6. Kobe plots WA.Fisheries style (Scenario 1)
  #25.6.1 by species  
for(i in 1:N.sp)
{
  if(!is.null(Age.based$SS$estimates[[i]]))
  {
    print(paste("SS3 --- Kobe plot WA Fisheries-----",Keep.species[i]))
    fn.get.Kobe.plot_appendix_WA.Fisheries(d=Age.based,
                                           sp=Keep.species[i],
                                           RF=Ref.points)
    ggsave(paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_Kobe_plot_WA_Fisheries.tiff",sep=''),
           width = 10,height = 10, dpi = 300,compression = "lzw")
  }
}
  #25.6.2 Display combined species  
for(l in 1:length(Lista.sp.outputs))
{
  Nms=names(compact(Age.based$SS$estimates))
  this.sp=Lista.sp.outputs[[l]]
  this.sp=this.sp[which(this.sp%in%Nms)]
  if(length(this.sp)>0)
  {
    print(paste("SS3 --- Kobe plot WA Fisheries-----",names(Lista.sp.outputs)[l]))
    DIMS=n2mfrow(length(this.sp))
    NKOL=DIMS[2]
    NRW=DIMS[1]
    if(NKOL%in%3:4) WIZ=14
    if(NKOL==2) WIZ=11
    if(NKOL==1) WIZ=9
    fn.get.Kobe.plot_WA.Fisheries(this.sp,
                                  d=Age.based$SS,
                                  NKOL,
                                  NRW,
                                  RF=Ref.points)
    ggsave(paste(HandL.out.RAR,'/Kobe_plot_SS3 integrated_WA_Fisheries_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = WIZ,height = 12,compression = "lzw")
  }
}

#25.7 Kobe plots SAFS style (Scenario 1)
  #25.7.1 by species  
for(i in 1:N.sp)
{
  if(!is.null(Age.based$SS$estimates[[i]]))
  {
    print(paste("SS3 --- Kobe plot SAFS-----",Keep.species[i]))
    fn.get.Kobe.plot_appendix_SAFS(d=Age.based,
                                   sp=Keep.species[i],
                                   RF=Ref.points)
    ggsave(paste(HandL.out,capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_Kobe_plot_SAFS.tiff",sep=''),
           width = 10,height = 10, dpi = 300,compression = "lzw")
  }
}

  #25.7.2 Display combined species  
for(l in 1:length(Lista.sp.outputs))
{
  Nms=names(compact(Age.based$SS$estimates))
  this.sp=Lista.sp.outputs[[l]]
  this.sp=this.sp[which(this.sp%in%Nms)]
  if(length(this.sp)>0)
  {
    print(paste("SS3 --- Kobe plot SAFS-----",names(Lista.sp.outputs)[l]))
    DIMS=n2mfrow(length(this.sp))
    NKOL=DIMS[2]
    NRW=DIMS[1]
    if(NKOL%in%3:4) WIZ=14
    if(NKOL==2) WIZ=11
    if(NKOL==1) WIZ=9
    fn.get.Kobe.plot_SAFS(this.sp,
                          d=Age.based$SS,
                          NKOL,
                          NRW,
                          RF=Ref.points,
                          Scen='S1')
    ggsave(paste(HandL.out.RAR,'/Kobe_plot_SS3 integrated_SAFS_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = WIZ,height = 12,compression = "lzw")
  }
}

#25.8 Store Consequence and likelihood for WoE 
get.cons.like.SS=FALSE   #Superseded. Now this is extracted from exported tables
if(get.cons.like.SS) Store.cons.Like_Age.based=fn.get.cons.like(lista=Age.based) 



