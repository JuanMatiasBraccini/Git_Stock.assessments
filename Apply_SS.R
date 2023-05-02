#---Run models -------------------------------------------------
n.SS=length(Integrated.age.based)  
Age.based=vector('list',n.SS)
names(Age.based)=Integrated.age.based

for(w in 1:n.SS)
{
  # SS3
  if(names(Age.based)[w]=="SS") #takes XX secs per iteration per species per scenario
  {
    dummy.store=vector('list',N.sp)
    names(dummy.store)=Keep.species
    dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
      dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.Kobe.probs=
      dummy.store.sens.table=dummy.store.ensemble=dummy.store
    
    for(i in 1:length(dummy.store))
    {
      Neim=names(dummy.store)[i]
      
      if(!is.null(Catch.rate.series[[i]]))
      {
        this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                      capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
        if(!dir.exists(this.wd))dir.create(this.wd)
        
        
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
          if(CV.use=='fixed') use.this.CV=default.CV
          if(use.this.CV<default.CV) use.this.CV=default.CV
          #if(use.this.CV>0.6) use.this.CV=0.6
          meanbodywt.SS.format=meanbodywt.SS.format%>%
                                mutate(CV=ifelse(CV<default.CV,use.this.CV,CV))
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
              CPUE[[j]]=CPUE[[j]]%>%
                mutate(CV=ifelse(CV<default.CV,use.this.CV,CV))
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
        # only applicable to gummy and whiskery (only for these species length-@-age data collected from gillnet fishery)
        MeanSize.at.Age.obs.SS.format=NULL
        if(any(grepl('age_length',names(Species.data[[i]]))) & names(Species.data)[i]%in%c("gummy shark","whiskery shark" ))
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
        if(len.cpue>0)
        {
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
            
            #remove daily 2007  
            if(drop2007.daily & nm=='Southern.shark_2')
            {
              dd=dd%>%
                mutate(Mean=ifelse(Year%in%2007,NA,Mean))
            }
            CPUE[[x]]=dd%>%filter(!is.na(Mean))
           }
          if(!is.null(CPUE))
          {
            Abundance.SS.format=do.call(rbind,CPUE)%>%
                                    relocate(Year,seas,index,Mean,CV)%>%
                                    arrange(index,Year)
          }
          
          
          #Scenarios
          Scens=List.sp[[i]]$Sens.test$SS3%>%
            mutate(Species=capitalize(Neim))
          Store.sens=vector('list',nrow(Scens))
          names(Store.sens)=Scens$Scenario
          Out.Scens=Scens
          Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
            Out.B.Bmsy=Out.F.Fmsy=store.warnings=store.convergence=vector('list',length(Store.sens))
          
          #Life history
          Life.history=List.sp[[i]]
          Life.history$Fecundity=ceiling(mean(Life.history$Fecundity))
          Life.history$Max.age.F=ceiling(mean(Life.history$Max.age.F))
          Life.history$Breed.cycle=mean(Life.history$Breed.cycle)
          
          #Likelihood lambdas
          Lamdas.SS.lambdas=Life.history$SS_lambdas
          
          #Execute SS
          for(s in 1:length(Store.sens))
          {
            this.wd1=paste(this.wd,names(Store.sens)[s],sep='/')
            if(!dir.exists(this.wd1))dir.create(this.wd1)
            
            names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))]=match(names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))],names(Flits))
            
            #a. Create input files
            fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SS'),   
                         new.path=this.wd1,
                         Scenario=Scens[s,]%>%
                                   mutate(Model='SS',
                                          Ln_R0_init=runif(1,Ln_R0_max*.3,Ln_R0_max*.6)),
                         Catch=ktch,
                         life.history=Life.history,
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
            #Arg='-nohess'
            Arg=''
            fn.run.SS(where.inputs=this.wd1,
                      where.exe=handl_OneDrive('SS3/ss_win.exe'),
                      args=Arg)  
            
            #c. Bring in outputs
            Report=SS_output(this.wd1,covar=F,forecast=F,readwt=F,checkcor=F)
            
            #d. Store estimates
            Estims=do.call(rbind,fn.get.stuff.from.list(Report,"estimated_non_dev_parameters"))
            Out.estimates[[s]]=Estims%>%
              mutate(Par=gsub('[0-9]+', '', rownames(Estims)))%>%
              group_by(Par)%>%
              summarise(Median=median(Value),
                        Lower.95=quantile(Value,probs=0.025),
                        Upper.95=quantile(Value,probs=0.975))%>%
              dplyr::select(Par,Median,Lower.95,Upper.95)
            
            #e. Store trajectories
            dummy=fn.ktch.only.get.timeseries(d=Report,
                                              mods=names(Catch_only)[w],
                                              Type='Depletion',
                                              scen=Scens$Scenario[s],
                                              Katch=ktch$Tonnes)
            Out.rel.biom[[s]]=dummy$Dat
            Out.probs.rel.biom[[s]]=dummy$Probs
            
            dummy=fn.ktch.only.get.timeseries(d=Report,
                                              mods=names(Catch_only)[w],
                                              Type='F.series',
                                              scen=Scens$Scenario[s],
                                              Katch=ktch$Tonnes)
            Out.f.series[[s]]=dummy$Dat
            
            dummy=fn.ktch.only.get.timeseries(d=Report,
                                              mods=names(Catch_only)[w],
                                              Type='B.Bmsy',
                                              scen=Scens$Scenario[s],
                                              Katch=ktch$Tonnes)
            Out.B.Bmsy[[s]]=dummy$Dat
            
            dummy=fn.ktch.only.get.timeseries(d=Report,
                                              mods=names(Catch_only)[w],
                                              Type='F.Fmsy',
                                              scen=Scens$Scenario[s],
                                              Katch=ktch$Tonnes)
            Out.F.Fmsy[[s]]=dummy$Dat
            
            
            
            # if(Scens$Scenario[s]=='S1') Out.Kobe.probs=Store.sens[[s]]$output$kobe%>%dplyr::select(stock,harvest)  MISSING: where are Kobe plots?
            
            
            print(paste("___________","SS3 Scenario",Scens$Scenario[s],"___________",Neim))
            
            
          }
          Out.Scens=Out.Scens%>%
            rename(h.mean=Steepness,
                   h.sd=Steepness.sd)%>%
            mutate(Mmean=round(Mmean,3),
                   h.mean=round(h.mean,3),
                   h.sd=round(h.sd,3))%>%
            dplyr::select(-c(Model,use_F_ballpark, Sims,Final.dpl))%>%
            relocate(Species,Scenario,Mmean,F_ballpark,h.mean,h.sd)
          
          
          #5. Store quantities
          dummy.store.sens.table[[i]]=Out.Scens
          dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
          dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
          dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
          dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
          dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
          #dummy.store.Kobe.probs[[i]]=Out.Kobe.probs  
          dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
          
          rm(Out.Scens,Out.rel.biom,Out.probs.rel.biom,Out.f.series,
             Out.B.Bmsy,Out.F.Fmsy,Out.estimates)
          
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
            annotation_custom(tableGrob(dummy.store.sens.table[[i]]%>%dplyr::select(Scenario,Mmean,h.mean,Final.bio.low,Final.bio.up)),
                              xmin=xmin+5, xmax=xmax+5, ymin=0, ymax=0.3)+
            annotation_custom(tableGrob(dummy.store.estimates[[i]]),xmin=xmin, xmax=xmax, ymin=0.35, ymax=0.55)
          print(p)
          ggsave(paste(this.wd,"/Rel.biomass&estimates.tiff",sep=''))
        }
      }
    }
    
    Age.based[[w]]$sens.table=dummy.store.sens.table
    Age.based[[w]]$estimates=dummy.store.estimates
    Age.based[[w]]$rel.biom=dummy.store.rel.biom
    Age.based[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
    Age.based[[w]]$f.series=dummy.store.f.series
    Age.based[[w]]$B.Bmsy=dummy.store.B.Bmsy
    Age.based[[w]]$F.Fmsy=dummy.store.F.Fmsy
    Age.based[[w]]$Kobe.probs=dummy.store.Kobe.probs
    
    
    rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,
       dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,
       dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.Kobe.probs)
    
  }
}




#---Generate outputs -------------------------------------------------