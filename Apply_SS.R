#---Run models -------------------------------------------------
n.SS=length(Integrated.age.based)  
Age.based=vector('list',n.SS)
names(Age.based)=Integrated.age.based
tic("timer")
for(w in 1:n.SS)
{
  # SS3
  if(names(Age.based)[w]=="SS") #takes 1000 secs per species per scenario (with Hessian and MC simulation)
  {
    if(do.parallel.SS) 
    {
      set.seed(1234)
      progress <- function(n) cat(sprintf(": SS3 fit complete for -----------", n),Keep.species[n],"\n")
      opts <- list(progress = progress)
      cl <- makeCluster(detectCores()-1)
      registerDoSNOW(cl)
      out.species=foreach(i= 1:N.sp,.options.snow = opts,.packages=c('gridExtra','Hmisc','JABBA','strex',
                                                                     'TruncatedDistributions','tidyverse','r4ss',
                                                                     'mvtnorm','ggrepel','ss3diags')) %dopar%
      {
        Neim=Keep.species[i]
        
        if((!is.null(Catch.rate.series[[i]]) | Neim%in%Species.with.length.comp) & !(Neim%in%no.empirical.sel.main.fleet))
        {
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
                          ifelse(FishCubeCode%in%c('WRL') & Neim%in%WRL.species,'WRL',
                          'Other'))))%>%
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
          
          
          #2. Size composition   
          #note: commercial catch and survey. Nsamp set at number of shots
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
          
          
          #4. Abundance series
          Abundance.SS.format=NULL
          Var.ad.factr=NULL
          CPUE=compact(Catch.rate.series[[i]])
          if(!is.null(CPUE))
          {
            DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
            if(length(DROP)>0)CPUE=CPUE[-DROP]
            if(Neim%in%survey_not.representative & any(grepl("Survey",names(CPUE)))) CPUE=CPUE[-grep("Survey",names(CPUE))]
            if(Neim%in%NSF_not.representative & any(grepl("NSF",names(CPUE)))) CPUE=CPUE[-grep("NSF",names(CPUE))]
            if(Neim%in%tdgdlf_not.representative & any(grepl("TDGDLF",names(CPUE)))) CPUE=CPUE[-grep("TDGDLF",names(CPUE))]
            if(Neim%in%tdgdlf_monthly_not.representative & "TDGDLF.monthly"%in%names(CPUE)) CPUE=CPUE[-match("TDGDLF.monthly",names(CPUE))]
            if(!is.null(Life.history$drop.monthly.cpue)) CPUE$TDGDLF.monthly=CPUE$TDGDLF.monthly%>%filter(!yr.f%in%Life.history$drop.monthly.cpue)
            
            if(length(CPUE)>0)
            {
              #reset very low CVs
              #note: Andre suggested leaving original CVs and estimating extraSD if more than one index available
              #      If only 1 index available, then do not estimate, just increase CV before fitting model
              #      ICCAT and SEDAR leave CVs as is and don't estimate extraSD but add variance adjustments factors to the control file
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
           }
          
          #Add size comp effective sample size bias adjustment    
          #note: a Value of 0 means no effect
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
          
          
          #5. F from tagging studies on TDGDLF (1994-95 and 2001-03)
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
          
          
          #6. Conditional age at length
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
          
          
          #7. MeanSize at Age obs
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
          len.size.comp=length(Size.compo.SS.format) 
          if(len.cpue>0|len.size.comp>0)
          {
            if(len.cpue>0)
            {
              MAX.CV=Life.history$MAX.CV
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
            
            #Scenarios
            Scens=Life.history$Sens.test$SS%>%
              mutate(Species=capitalize(Neim))
            Store.sens=vector('list',nrow(Scens))
            names(Store.sens)=Scens$Scenario
            Out.Scens=Scens
            Out.estimates=Out.quantities=Out.likelihoods=Out.rel.biom=Out.probs.rel.biom=Out.f.series=Out.probs.f.series=Out.B.Bmsy=
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
            
            #rename fleets following SS nomenclature
            names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))]=match(names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))],names(Flits))
            
            #future catches
            if("SS"%in%future.models)
            {
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
                drop.dis=Abundance.SS.format%>%
                          mutate(this=grepl('TDGDLF.daily',rownames(Abundance.SS.format)) & Year%in%rid.of)
                if(any(drop.dis$this)) Abundance.SS.format=Abundance.SS.format[-which(drop.dis$this),]
              }
              
              #use cpue or length comp in likelihood?
              if(Life.history$drop.length.comp) Size.compo.SS.format=NULL
              if(Life.history$drop.cpue) Abundance.SS.format=NULL
              
              #a. Create input files
              if(create.SS.inputs)
              {
                fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SS'),   
                             new.path=this.wd1,
                             Scenario=Scens[s,]%>%
                               mutate(Model='SS'),
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
                             Lamdas=Lamdas.SS.lambdas,
                             Var.adjust.factor=Var.ad.factr,
                             Future.project=add.ct.future)
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
              if(!Calculate.ramp.years)
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
              
              #d. Store estimates, likelihoods and management quantities
              Out.quantities[[s]]=Report[["derived_quants"]]%>%
                filter(Label%in%c("Dead_Catch_MSY",paste0("Bratio_",Last.yr.ktch.numeric)))%>%
                dplyr::select(Label,Value,StdDev)%>%
                mutate(Label=case_when(grepl('Bratio',Label)~'Current depletion',
                                       Label=='Dead_Catch_MSY'~'MSY'),
                       Scenario=Scens$Scenario[s])
              
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
                                    
              #e. Store trajectories
              #note: uncertainty is based on asymptotic error
              
                #e.1 relative biomass
              dummy=fn.integrated.mod.get.timeseries(d=Report,
                                                     mods="SS3",
                                                     Type='Depletion',
                                                     scen=Scens$Scenario[s])
              Out.rel.biom[[s]]=dummy$Dat
              Out.probs.rel.biom[[s]]=dummy$Probs 
              
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
              
              #Posterior vs prior
              if(Report$parameters%>%filter(Label=='SR_BH_steep')%>%pull(Phase)>0)
              {
                pdf(paste(this.wd,paste('Steepness prior vs posterior_',names(Store.sens)[s],'.pdf',sep=''),sep='/'))
                fn.compare.prior.post(d=Report$parameters%>%filter(Label=='SR_BH_steep'),
                                      Par='Steepness',
                                      prior_type='beta')
                dev.off()
              }
              
              # Evaluate fit diagnostics
              GoodnessFit=function.goodness.fit_SS(Rep=Report)  
              write.csv(GoodnessFit,paste(this.wd1,"/GoodnessFit_",Neim,".csv",sep=''))
            }  #end s loop
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
                     SR_sigmaR=Life.history$SR_sigmaR,
                     a_FL.to.TL=Life.history$a_FL.to.TL,
                     b_FL.to.TL=Life.history$b_FL.to.TL,
                     Fleets=paste(flitinfo$fleetname,collapse=', '))
            
            
            #5. Store quantities
            return(list(dummy.store.sens.table=Out.Scens,
                        dummy.store.rel.biom=do.call(rbind,Out.rel.biom),
                        dummy.store.probs.rel.biom=Out.probs.rel.biom,
                        dummy.store.probs.B.Bmsy=Out.probs.B.Bmsy,
                        dummy.store.probs.f.series=Out.probs.f.series,
                        dummy.store.f.series=do.call(rbind,Out.f.series),
                        dummy.store.B.Bmsy=do.call(rbind,Out.B.Bmsy),
                        dummy.store.F.Fmsy=do.call(rbind,Out.F.Fmsy),
                        dummy.store.Kobe.probs=Out.Kobe.probs,  
                        dummy.store.estimates=do.call(rbind,Out.estimates),
                        dummy.store.quantities=do.call(rbind,Out.quantities),
                        dummy.store.likelihoods=do.call(rbind,Out.likelihoods)  
                        ))
            
            rm(Out.Scens,Out.rel.biom,Out.probs.rel.biom,Out.f.series,Out.probs.f.series,
               Out.B.Bmsy,Out.F.Fmsy,Out.estimates,Out.Kobe.probs,Out.likelihoods) 
            
          }
          clear.log("Var.ad.factr")
        }
      }
      stopCluster(cl)
      names(out.species)=Keep.species
      
      Age.based[[w]]$sens.table=fn.get.and.name(LISTA=out.species,x='dummy.store.sens.table')
      Age.based[[w]]$estimates=fn.get.and.name(LISTA=out.species,x='dummy.store.estimates')
      Age.based[[w]]$quantities=fn.get.and.name(LISTA=out.species,x='dummy.store.quantities')
      Age.based[[w]]$rel.biom=fn.get.and.name(LISTA=out.species,x='dummy.store.rel.biom')
      Age.based[[w]]$probs.rel.biom=fn.get.and.name(LISTA=out.species,x='dummy.store.probs.rel.biom')
      Age.based[[w]]$probs.B.Bmsy=fn.get.and.name(LISTA=out.species,x='dummy.store.probs.B.Bmsy') 
      Age.based[[w]]$f.series=fn.get.and.name(LISTA=out.species,x='dummy.store.f.series')
      Age.based[[w]]$B.Bmsy=fn.get.and.name(LISTA=out.species,x='dummy.store.B.Bmsy')
      Age.based[[w]]$F.Fmsy=fn.get.and.name(LISTA=out.species,x='dummy.store.F.Fmsy')
      Age.based[[w]]$Kobe.probs=fn.get.and.name(LISTA=out.species,x='dummy.store.Kobe.probs')
      Age.based[[w]]$likelihoods=fn.get.and.name(LISTA=out.species,x='dummy.store.likelihoods')
      
      #plot overall figure with estimates and scenarios
      for(i in 1:N.sp)
      {
        print(paste("SS3 rel biom figure and estimates ____",Keep.species[i]))
        if(!is.null(Age.based[[w]]$rel.biom[[i]]))
        {
          yrS=Age.based[[w]]$rel.biom[[i]]%>%filter(Scenario=='S1')%>%pull(year)
          if(length(yrS)<50) delta=10 else
            delta=17
          xmin=min(yrS)+delta
          xmax=xmin+delta
          TAB=Age.based[[w]]$estimates[[i]]%>%
            filter(Par=='SR_LN(R0)')%>%
            dplyr::select(Par,Value,Min,Max,Init,Status,Gradient)%>%
            `rownames<-`( NULL )
          p=Age.based[[w]]$rel.biom[[i]]%>%
            ggplot(aes(year,median,color=Scenario))+
            annotation_custom(tableGrob(Age.based[[w]]$sens.table[[i]]%>%
                                          dplyr::select(Scenario,Mmean,h.mean)),
                              xmin=xmin+5, xmax=xmax+5, ymin=0.45, ymax=0.65)+
            annotation_custom(tableGrob(TAB),xmin=xmin+2, xmax=xmax, ymin=0 , ymax=0.3)+
            geom_line(size=2)+
            geom_line(aes(year,upper.95),linetype=2)+
            geom_line(aes(year,lower.95),linetype=2)+
            ggtitle(Keep.species[i])+ylim(0,1)+
            theme_PA()+theme(legend.position = 'bottom')
          print(p)
          ggsave(paste(this.wd,"/Rel.biomass&estimates.tiff",sep=''),compression = "lzw")  
          rm(TAB)
          
        }
       }

      rm(out.species)
      
    } 
    if(!do.parallel.SS)
    {
      dummy.store=vector('list',N.sp)
      names(dummy.store)=Keep.species
      dummy.store.estimates=dummy.store.likelihoods=dummy.store.rel.biom=dummy.store.probs.rel.biom=dummy.store.probs.f.series=
        dummy.store.probs.B.Bmsy=dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.Kobe.probs=
        dummy.store.sens.table=dummy.store.ensemble=dummy.store.quantities=dummy.store
      
      for(i in 1:length(dummy.store))
      {
        Neim=names(dummy.store)[i]
        
        if((!is.null(Catch.rate.series[[i]]) | Neim%in%Species.with.length.comp) & !(Neim%in%no.empirical.sel.main.fleet))
        {
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
                          ifelse(FishCubeCode%in%c('WRL') & Neim%in%WRL.species,'WRL',
                          'Other'))))%>%
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
          
          
          #2. Size composition   
          #note: commercial catch and survey. Nsamp set at number of shots
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
          
          
          #4. Abundance series
          Abundance.SS.format=NULL
          Var.ad.factr=NULL
          CPUE=compact(Catch.rate.series[[i]])
          if(!is.null(CPUE))
          {
            DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
            if(length(DROP)>0)CPUE=CPUE[-DROP]
            if(Neim%in%survey_not.representative & any(grepl("Survey",names(CPUE)))) CPUE=CPUE[-grep("Survey",names(CPUE))]
            if(Neim%in%NSF_not.representative & any(grepl("NSF",names(CPUE)))) CPUE=CPUE[-grep("NSF",names(CPUE))]
            if(Neim%in%tdgdlf_not.representative & any(grepl("TDGDLF",names(CPUE)))) CPUE=CPUE[-grep("TDGDLF",names(CPUE))]
            if(Neim%in%tdgdlf_monthly_not.representative & "TDGDLF.monthly"%in%names(CPUE)) CPUE=CPUE[-match("TDGDLF.monthly",names(CPUE))]
            if(!is.null(Life.history$drop.monthly.cpue)) CPUE$TDGDLF.monthly=CPUE$TDGDLF.monthly%>%filter(!yr.f%in%Life.history$drop.monthly.cpue)
            
            if(length(CPUE)>0)
            {
              #reset very low CVs
              #note: Andre suggested leaving original CVs and estimating extraSD if more than one index available
              #      If only 1 index available, then do not estimate, just increase CV before fitting model
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
           }

          
          #Add size comp effective sample size bias adjustment    
          #note: a Value of 0 means no effect
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
          
          
          #5. F from tagging studies on TDGDLF (1994-95 and 2001-03)
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
          
          
          #6. Conditional age at length
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
          
          
          #7. MeanSize at Age obs
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
          len.size.comp=length(Size.compo.SS.format)
          if(len.cpue>0|len.size.comp>0)
          {
            if(len.cpue>0)
            {
              MAX.CV=Life.history$MAX.CV
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
            
            #Scenarios
            Scens=Life.history$Sens.test$SS%>%
              mutate(Species=capitalize(Neim))
            Store.sens=vector('list',nrow(Scens))
            names(Store.sens)=Scens$Scenario
            Out.Scens=Scens
            Out.estimates=Out.quantities=Out.likelihoods=Out.rel.biom=Out.probs.rel.biom=Out.f.series=Out.B.Bmsy=Out.probs.f.series=
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
            
            #rename fleets following SS nomenclature
            names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))]=match(names(ktch)[which(!names(ktch)%in%c("SPECIES","Name","finyear"))],names(Flits))
            
            #future catches
            if("SS"%in%future.models)
            {
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
                drop.dis=Abundance.SS.format%>%
                  mutate(this=grepl('TDGDLF.daily',rownames(Abundance.SS.format)) & Year%in%rid.of)
                if(any(drop.dis$this)) Abundance.SS.format=Abundance.SS.format[-which(drop.dis$this),]
              }
              
              #use cpue or length comp in likelihood?
              if(Life.history$drop.length.comp) Size.compo.SS.format=NULL
              if(Life.history$drop.cpue) Abundance.SS.format=NULL
              
              #a. Create input files
              if(create.SS.inputs)
              {
                fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SS'),   
                             new.path=this.wd1,
                             Scenario=Scens[s,]%>%
                               mutate(Model='SS'),
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
                             Lamdas=Lamdas.SS.lambdas,
                             Var.adjust.factor=Var.ad.factr,
                             Future.project=add.ct.future)
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
              if(!Calculate.ramp.years)
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
              
              #d. Store estimates, likelihoods and management quantities
              Out.quantities[[s]]=Report[["derived_quants"]]%>%
                filter(Label%in%c("Dead_Catch_MSY",paste0("Bratio_",Last.yr.ktch.numeric)))%>%
                dplyr::select(Label,Value,StdDev)%>%
                mutate(Label=case_when(grepl('Bratio',Label)~'Current depletion',
                                       Label=='Dead_Catch_MSY'~'MSY'),
                       Scenario=Scens$Scenario[s])
              
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
                     SR_sigmaR=Life.history$SR_sigmaR,
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
            dummy.store.quantities[[i]]=do.call(rbind,Out.quantities)
            dummy.store.likelihoods[[i]]=do.call(rbind,Out.likelihoods)
            
            rm(Out.Scens,Out.rel.biom,Out.probs.rel.biom,Out.f.series,
               Out.B.Bmsy,Out.F.Fmsy,Out.estimates,Out.Kobe.probs,Out.likelihoods)
            
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
        }
      } #end i loop
      
      Age.based[[w]]$sens.table=dummy.store.sens.table
      Age.based[[w]]$estimates=dummy.store.estimates
      Age.based[[w]]$quantities=dummy.store.quantities
      Age.based[[w]]$rel.biom=dummy.store.rel.biom
      Age.based[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
      Age.based[[w]]$probs.B.Bmsy=dummy.store.probs.B.Bmsy
      Age.based[[w]]$probs.f.series=dummy.store.probs.f.series
      Age.based[[w]]$f.series=dummy.store.f.series
      Age.based[[w]]$B.Bmsy=dummy.store.B.Bmsy
      Age.based[[w]]$F.Fmsy=dummy.store.F.Fmsy
      Age.based[[w]]$Kobe.probs=dummy.store.Kobe.probs
      Age.based[[w]]$likelihoods=dummy.store.likelihoods
      
      
      rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,dummy.store.probs.f.series,
         dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,dummy.store.probs.B.Bmsy,
         dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.Kobe.probs,dummy.store.likelihoods)
      
    }
  }
}
toc(log = TRUE, quiet = TRUE)
computation.time <- tic.log(format = TRUE)
tic.clearlog()
send.email(TO=Send.email.to,
           CC='',
           Subject=paste("SS3 models finished running at",Sys.time()),
           Body= paste("Computation",computation.time),  
           Attachment=NULL) 


#---Do SS3 diagnostics -------------------------------------------------
#notes: 
        # runs test only displays for series with >1 year
        # Hindcasting cross-validation only available for species with abundance series
if(do.SS3.diagnostics)
{
  tic("timer")
  set.seed(1234)
  progress <- function(n) cat(sprintf(": SS3 fit diagnostics complete for -----------", n),Keep.species[n],"\n")
  #progress <- function(n) setTxtProgressBar(txtProgressBar(max = N.sp, style = 3), n)
  opts <- list(progress = progress)
  cl <- makeCluster(detectCores()-1)
  registerDoSNOW(cl)
  out=foreach(l= 1:N.sp,.options.snow = opts,.packages=c('Hmisc','tidyverse','r4ss','ss3diags')) %dopar%
    {
      Neim=Keep.species[l]
      this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
      this.wd1=paste(this.wd,"S1",sep='/')
      if(file.exists(this.wd1))
      {
          MLE=read.admbFit(paste(this.wd1,'ss',sep='/'))
          Estim.LnRo=MLE$est[grep("SR_parm",MLE$names)]
          R0.range=seq(Estim.LnRo*(1-delta.likelihood.profiles),Estim.LnRo*(1.3+delta.likelihood.profiles),length.out=Number.of.likelihood.profiles)
          fn.fit.diag_SS3(WD=this.wd1,
                          do.like.prof=TRUE,
                          disfiles=c("control.ss_new", "data.dat","forecast.ss","starter.ss","Report.sso"),
                          R0.vec=R0.range,
                          exe_path=handl_OneDrive('SS3/ss_win.exe'),
                          start.retro=Retro_start,
                          end.retro=Retro_end,
                          do.retros=TRUE,
                          do.jitter=TRUE,
                          numjitter=Number.of.jitters)
          rm(MLE,this.wd1,R0.range,Estim.LnRo)
        }
    }
  stopCluster(cl)
  toc(log = TRUE, quiet = TRUE)
  computation.time <- tic.log(format = TRUE)
  tic.clearlog()
  send.email(TO=Send.email.to,
             CC='',
             Subject=paste("SS3 model diagnostics finished running at",Sys.time()),
             Body= paste("Computation",computation.time),  
             Attachment=NULL) 
}


#---Compare different scenarios for each species -------------------------------------------------
if(SS3.run=='final')
{
  for(i in 1:N.sp)
  {
    Neim=Keep.species[i]
    this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                  capitalize(Neim),"/",AssessYr,"/SS3 integrated",sep='')
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
    d.tbl=do.call(rbind,dummy)%>%relocate(Species)
    if(!evaluate.07.08.cpue) d.tbl=d.tbl%>%dplyr::select(-Daily.cpues)
    write.csv(d.tbl,
              paste(Rar.path,paste('Table 10. Age.based_',names(Age.based)[w],'_scenarios_',
                                   names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
              row.names = F)
  }
}

#25.2 Table of parameter estimates, likelihoods and quantites of interest by species, scenario and method
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
      dummy=do.call(rbind,dummy)%>%
        rownames_to_column(var = "Species")%>%
        mutate(Species=capitalize(str_extract(Species, "[^.]+")))%>%
        relocate(Species)
      write.csv(dummy,paste(Rar.path,paste('Table 11. Age.based_',
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
      dummy=do.call(rbind,dummy)%>%
        rownames_to_column(var = "Species")%>%
        mutate(Species=capitalize(str_extract(Species, "[^.]+")))%>%
        relocate(Species)
      write.csv(dummy,paste(Rar.path,paste('Table 11. Age.based_',
                                           names(Age.based)[w],'_likelihoods_',
                                           names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
                row.names = F)
    }
    
    #quantities of interest 
    dummy=Age.based[[w]]$quantities
    dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
    dummy=compact(dummy)
    if(length(dummy)>0)
    {
      dummy=do.call(rbind,dummy)%>%
        rownames_to_column(var = "Species")%>%
        mutate(Species=capitalize(str_extract(Species, "[^.]+")))%>%
        relocate(Species)%>%
        rename(Median=Value,
               SE=StdDev)
      write.csv(dummy,paste(Rar.path,paste('Table 11. Age.based_',
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
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_time_series_relative_biomass.tiff",sep=''),
             width = 6,height = 10,compression = "lzw")
      
      #export current depletion probabilities
      write.csv(a$store.probs%>%
                  spread(Model,Probability)%>%
                  mutate(Species=Keep.species[i],
                         Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                  arrange(Scenario,Range),
                paste(handl_OneDrive("Analyses/Population dynamics/1."),
                      capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_current_depletion_probabilities.csv",sep=''),
                row.names = F)
      
      #export current depletion 
      xx=Age.based$SS$rel.biom[[i]]
      write.csv(xx,
                paste(handl_OneDrive("Analyses/Population dynamics/1."),
                      capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_depletion.csv",sep=''),
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
        ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_time_series_fishing_mortality.tiff",sep=''),
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
        ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_time_series_B_Bmsy.tiff",sep=''),
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
        ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_time_series_F_Fmsy.tiff",sep=''),
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
  if(!is.null(a))ggsave(paste(Rar.path,'/Relative.biomass_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
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
      dummy.mod[[m]]=do.call(rbind,dummy)%>%
        mutate(Model=names(Age.based)[m])
    }
  }
  dummy.mod=compact(dummy.mod)
  if(length(dummy.mod)>0)
  {
    write.csv(do.call(rbind,dummy.mod)%>%
                mutate(Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                spread(Species,Probability)%>%
                arrange(Range),
              paste(Rar.path,'/Table 12. Age.based_SS_current.depletion_',names(Lista.sp.outputs)[l],'.csv',sep=''),
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
    if(!is.null(a))ggsave(paste(Rar.path,'/B.over.Bmsy_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
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
        dummy.mod[[m]]=do.call(rbind,dummy)%>%
          mutate(Model=names(Age.based)[m])
      }
    }
    dummy.mod=compact(dummy.mod)
    if(length(dummy.mod)>0)
    {
      write.csv(do.call(rbind,dummy.mod)%>%
                  mutate(Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                  spread(Species,Probability)%>%
                  arrange(Range),
                paste(Rar.path,'/Table 12. Age.based_SS_Current.B.over.Bmsy_',names(Lista.sp.outputs)[l],'.csv',sep=''),
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
    if(!is.null(a))ggsave(paste(Rar.path,'/F.over.Fmsy_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
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
    if(!is.null(a))ggsave(paste(Rar.path,'/F.series_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
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
        dummy.mod[[m]]=do.call(rbind,dummy)%>%
          mutate(Model=names(Age.based)[m])
      }
    }
    dummy.mod=compact(dummy.mod)
    if(length(dummy.mod)>0)
    {
      write.csv(do.call(rbind,dummy.mod)%>%
                  mutate(Range=factor(Range,levels=c(">lim","lim.thr","thr.tar","<tar")))%>%
                  spread(Species,Probability)%>%
                  arrange(Range),
                paste(Rar.path,'/Table 12. Age.based_SS_Current.f_',names(Lista.sp.outputs)[l],'.csv',sep=''),
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
  if(!is.null(a))ggsave(paste(Rar.path,'/Relative.biomass_SS3 integrated_',names(Lista.sp.outputs)[l],'_sensitivity.tiff',sep=''),
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
      this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(List.sp[[i]]$Name),"/",AssessYr,"/SS3 integrated",sep='')
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
        stor.plt[[q]]=do.call(rbind,a)%>%
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
    ggsave(paste(Rar.path,'/GoodnessFit per scenario_SS3 integrated_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
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
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_Kobe_plot.tiff",sep=''),
           width = 9,height = 14, dpi = 300,compression = "lzw")
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
    ggsave(paste(Rar.path,'/Kobe_plot_Age.based_SS_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
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
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_Kobe_plot_WA_Fisheries.tiff",sep=''),
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
    ggsave(paste(Rar.path,'/Kobe_plot_Age.based_SS_WA_Fisheries_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
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
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/SS3 integrated/SS3_integrated_Kobe_plot_SAFS.tiff",sep=''),
           width = 10,height = 10, dpi = 300,compression = "lzw")
  }
}

#25.8 Store Consequence and likelihood for WoE 
get.cons.like.SS=FALSE   #Superseded. Now this is extracted from exported tables
if(get.cons.like.SS) Store.cons.Like_Age.based=fn.get.cons.like(lista=Age.based) 



