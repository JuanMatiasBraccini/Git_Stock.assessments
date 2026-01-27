#remotes::install_github("JABBAmodel/ss3diags")
library(ss3diags)  
# Create SS-DL tool input files ------------------------------------------------------
fn.create.SS_DL_tool_inputs=function(Life.history,CACH,this.wd,Neim,InputData,KtchRates)
{
  # iid=match(c("Size_composition_Observations","Size_composition_Pilbara_Trawl"),names(InputData))
  # iid=iid[!is.na(iid)]
  # if(length(iid)>0)InputData=InputData[-iid]
  if(!is.null(KtchRates) | any(grepl('Size_composition',names(InputData))))
  {
    
    #1. Catch
    ktch=CACH%>%
      filter(Name==Neim)%>%
      mutate(Fishry=ifelse(FishCubeCode%in%c('OANCGC','JANS','WANCS'),'Northern.shark',
                           ifelse(FishCubeCode%in%c('Historic','JASDGDL','WCDGDL','C070','OAWC',
                                                    'TEP_greynurse','TEP_dusky','Discards_TDGDLF'),'Southern.shark',
                                  'Other')))%>%
      group_by(SPECIES,Name,finyear,Fishry)%>%
      summarise(Tonnes=sum(LIVEWT.c,na.rm=T))%>%
      mutate(Fishry=case_when(Fishry=="Southern.shark" & finyear<2006 ~'Southern.shark_1',
                              Fishry=="Southern.shark" & finyear>=2006~'Southern.shark_2',
                              TRUE~Fishry),
             Year=finyear)%>%
      ungroup()
     
    
    #2. Fleets
    Flits.name=sort(unique(ktch$Fishry))  
    Flits=1:length(Flits.name)
    names(Flits)=Flits.name
    Flits.and.survey=data.frame(Fleet.number=c(Flits,1+length(Flits)),
                                Fleet.name=c(names(Flits),"Survey"))
    if(!"Survey"%in% names(KtchRates)) Flits.and.survey=Flits.and.survey%>%filter(!Fleet.name=="Survey")
    if("Survey"%in% names(KtchRates))
    {
      Flits=c(Flits,length(Flits)+1)
      names(Flits)[length(Flits)]="Survey"
    }
    
    
    #3. Size composition   
    #note: commercial catch and survey. Nsamp set at number of shots
    Size.compo.SS.format=NULL
    if(any(grepl('Size_composition',names(InputData))))
    {
      d.list.n.shots=InputData[grep(paste(c("Size_composition_Survey_Observations","Size_composition_Observations",
                                            "Size_composition_Other_Observations"),collapse="|"),
                                    names(InputData))]
      d.list=InputData[grep(paste(c("Size_composition_West","Size_composition_Zone1","Size_composition_Zone2",
                                    "Size_composition_NSF.LONGLINE","Size_composition_Survey",
                                    "Size_composition_Other"),collapse="|"),
                            names(InputData))]
      if(length(d.list)>0)
      {
        if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
        if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
        
        for(s in 1:length(d.list))
        {
          d.list[[s]]=d.list[[s]]%>%
            filter(FL>=Life.history$Lzero)%>%
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
        }
        d.list <- d.list[!is.na(d.list)]
        d.list=do.call(rbind,d.list)%>%
          group_by(year,fishry,sex,size.class)%>%
          summarise(n=sum(n))%>%
          ungroup()
        
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
                  mutate(dumi.n=rowSums(d.list[,-match(c('year','Seas','Fleet','Sex','Part','Nsamp'),names(d.list))]),
                         Nsamp=ifelse(Nsamp>dumi.n,dumi.n,Nsamp))%>%
                  dplyr::select(-dumi.n)
          size.flits=Flits.and.survey
          if(!"Survey"%in% names(KtchRates) & "Survey"%in%unique(d.list$Fleet))
          {
            ddummis=size.flits[1,]%>%mutate(Fleet.number=1+size.flits$Fleet.number[nrow(size.flits)],
                                            Fleet.name="Survey")
            rownames(ddummis)="Survey"
            size.flits=rbind(size.flits,ddummis)
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
          dummy.Size.compo.SS.format=d.list
          if(nrow(dummy.Size.compo.SS.format)>0) Size.compo.SS.format=dummy.Size.compo.SS.format
          
          Size.compo.SS.format=Size.compo.SS.format%>%
            rename(Year=year,
                   Month=Seas,
                   Nsamps=Nsamp)%>%
            dplyr::select(-Part)
          iidd=colSums(Size.compo.SS.format)
          iidd=iidd[which(iidd==0)]
          iidd=match(names(iidd),names(Size.compo.SS.format))
          if(length(iidd)>0) Size.compo.SS.format=Size.compo.SS.format[,-iidd]
          Size.compo.SS.format=Size.compo.SS.format%>%filter(!is.na(Fleet))
        }
      }
      

    }
    
    #4. Abundance series
    Abundance.SS.format=NULL
    CPUE=compact(KtchRates)
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
      
      MAX.CV=Life.history$MAX.CV
      len.cpue=length(CPUE)
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
      rid.of=c(2007,2008)
      drop.dis=Abundance.SS.format%>%
        mutate(this=grepl('TDGDLF.daily',rownames(Abundance.SS.format)) & Year%in%rid.of)
      ithis=which(drop.dis$this)
      if(length(ithis)>0)Abundance.SS.format=Abundance.SS.format[-ithis,]
      
      Abundance.SS.format=Abundance.SS.format%>%
        rename(Month=seas,
               Fleet=index,
               Index=Mean)
      Abundance.SS.format=Abundance.SS.format%>%
        left_join(Flits.and.survey%>%dplyr::rename(Fleet=Fleet.number),by="Fleet")%>%
        dplyr::rename(Label=Fleet.name)
      
    }
    
    #5. Life history
    out.LH=with(Life.history,
         {
           data.frame(Parameter=c('M','Linf','k','to','CV',
                                         'L50%','L95%','W.L_a','W.L_b','Fec_a','Fec_b',
                                         'M','Linf','k','to','CV','W.L_a','W.L_b','h'),
                             Sex=c(rep('F',11),rep('M',7),''),
                             Mean=c(Sens.test$SS3$Mmean[1],Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,Growth.F$k,Growth.F$to,0.1,
                                    TL.50.mat,TL.95.mat,AwT,BwT,mean(Life.history$Fecundity),0,
                                    Sens.test$SS3$Mmean[1],Growth.M$FL_inf*a_FL.to.TL+b_FL.to.TL,Growth.M$k,Growth.F$to,0.1,
                                    AwT.M,BwT.M,Sens.test$SS3$Steepness[1]),
                             SD=c(round(Sens.test$SS3$Msd[1],3),rep(0,4),rep('',6),
                                  round(Sens.test$SS3$Msd[1],3),rep(0,4),rep('',2),Sens.test$SS3$Steepness.sd[1]),
                             Prior.type=c('Lognormal',rep('',10),'Lognormal',rep('Normal',3),rep('',4)),
                             Phase=c(2,rep(3,3),-1,rep('',6),
                                     2,rep(3,3),-1,rep('',3)))
         })

    #export stuff
    if(!is.null(Abundance.SS.format) | !is.null(Size.compo.SS.format))
    {
      if(!dir.exists(this.wd))dir.create(this.wd)
      
      #export catch
      write.csv(ktch%>%
                  dplyr::select(Year,Fishry,Tonnes)%>%
                  spread(Fishry,Tonnes,fill=0),paste(this.wd,'Catches.csv',sep='/'),row.names = FALSE)
      
      #export lengths
      if(!is.null(Size.compo.SS.format)) write.csv(Size.compo.SS.format,paste(this.wd,'Lengths.csv',sep='/'),row.names = FALSE)
      
      #export abundance
      if(!is.null(Abundance.SS.format)) write.csv(Abundance.SS.format,paste(this.wd,'Index.csv',sep='/'),row.names = FALSE)
        
      #export life history
      write.csv(out.LH,paste(this.wd,'Life history.csv',sep='/'),row.names = FALSE)
    }
   }
}
# Create SS input files ------------------------------------------------------
fn.get.in.betwee=function(x,PATRN="_") str_before_nth(str_after_nth(x, PATRN, 1), PATRN, 2)
fn.set.up.SS=function(Templates,new.path,Scenario,Catch,life.history,depletion.yr,fleets=NULL,
                      fleetinfo=NULL,abundance=NULL,size.comp=NULL,age.comp=NULL,meanbodywt=NULL,
                      Tags=NULL,F.tagging=NULL,cond.age.len=NULL,MeanSize.at.Age.obs=NULL,Lamdas=NULL,
                      RecDev_Phase=-3,SR_sigmaR=0.2,Var.adjust.factor=NULL,Future.project=NULL,
                      first.age=0) #SS starts at age 0, if using 1, then Nages must be changed
{
  # 1.Copy templates
  copy_SS_inputs(dir.old = Templates, dir.new = new.path,overwrite = TRUE,verbose=FALSE)
  
  
  # 2.Read in templates 
  start <- r4ss::SS_readstarter(file = file.path(new.path, "starter.ss"), verbose = FALSE)
  dat <- r4ss::SS_readdat(file = file.path(new.path, start$datfile), verbose = FALSE)
  ctl <- r4ss::SS_readctl(file = file.path(new.path, start$ctlfile), verbose = FALSE, use_datlist = TRUE, datlist = dat)
  fore <- r4ss::SS_readforecast(file = file.path(new.path, "forecast.ss"),  verbose = FALSE)
  
  #dat.new <- r4ss::SS_readdat(file = file.path(new.path, "data_echo.ss_new"), verbose = FALSE)
  #ctl.new <- r4ss::SS_readctl(file = file.path(new.path, "control.ss_new"), verbose = FALSE, use_datlist = TRUE, datlist = dat.new)
  
  
  # 3.Update template with species-specific information
  
  #3.1. dat file
  dat$Comments[3]=paste('#C file write time:',Sys.time())
  dat$spawn_month=1.0  #integer is month (1-12), decimal is fraction of days in month (e.g. 15 March is 3.5)
  
  #general age info
  dat$Nages=max(life.history$Max.age.F)
  nages=dat$Nages+1  #including plus group
  ageError=as.data.frame(matrix(nrow=2,ncol=nages))   
  ageError[1,]=-1.00
  ageError[2,]=0.001
  names(ageError)=seq(first.age,dat$Nages)  
  dat$ageerror=ageError
  
  #population size classes
  dat$binwidth=TL.bins.cm
  dat$minimum_size=TL.bins.cm*floor(1*with(life.history,Lzero*a_FL.to.TL+b_FL.to.TL)/TL.bins.cm)
  dat$maximum_size=TL.bins.cm*ceiling(1*with(life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/TL.bins.cm)
  if(!is.null(life.history$Min.population.TL)) dat$minimum_size=life.history$Min.population.TL
  if(!is.null(life.history$Max.population.TL)) dat$maximum_size=life.history$Max.population.TL
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
    if(is.null(abundance)) get.fleet.ktch$catch_se=1e-3
    dat$catch=get.fleet.ktch
    dat$Nfleets=nrow(fleetinfo)
    dat$fleetinfo=fleetinfo  
    
    #cpue
    names(dat$CPUEinfo)=capitalize(names(dat$CPUEinfo))
    ddumy=dat$CPUEinfo[rep(1,length(dis.flits)),]
    row.names(ddumy)=dis.flits
    ddumy=ddumy%>%
            rownames_to_column('fleetname')%>%
            mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                             ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                                    fleetname)))%>%
            mutate(Fleet=row_number())
    if(Abundance.error.dist=='Normal') ddumy$Errtype=-1  
    if(Abundance.error.dist=='Lognormal') ddumy$Errtype=0
    rownames(ddumy)=ddumy$fleetname
    dat$CPUEinfo=ddumy%>%dplyr::select(-fleetname)
    if(!is.null(abundance)) dat$CPUE=abundance%>%mutate(Mean=ifelse(Mean<1e-6,1e-6,Mean))
    if(is.null(abundance))  dat$CPUE=NULL
    
    #Add WRL selectivity if appropriate 
    WRL.fleet=dat$fleetinfo%>%filter(fleetname=='WRL')
    if(nrow(WRL.fleet)>0)  
    {
      WRLinfo=dat$CPUEinfo[1,]%>%
        mutate(Fleet=match('WRL',dat$fleetinfo$fleetname))
      rownames(WRLinfo)='WRL'

      if("Survey"%in%dat$fleetinfo$fleetname)
      {
        dat$CPUEinfo=rbind(dat$CPUEinfo[-match('Survey',rownames(dat$CPUEinfo)),],WRLinfo,
                           dat$CPUEinfo[match('Survey',rownames(dat$CPUEinfo)),]%>%mutate(Fleet=match('Survey',dat$fleetinfo$fleetname)))
      }else
      {
        dat$CPUEinfo=rbind(dat$CPUEinfo,WRLinfo)
      }

    }
    
    #Fishing mortality from tagging
    #notes: to use an F series need to add an additional fleet, cpue series, q and mirror selectivity
    if(!is.null(F.tagging))
    {
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
      if(Abundance.error.dist=='Normal') Finfo$Errtype=-1  
      if(Abundance.error.dist=='Lognormal') Finfo$Errtype=0
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
      ddumy=dat$len_info[rep(1,length(dis.flits)),]
      row.names(ddumy)=dis.flits
      ddumy=ddumy%>%
        rownames_to_column('fleetname')%>%
        mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                         ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                         fleetname)))
      rownames(ddumy)=ddumy$fleetname
      dat$len_info=ddumy%>%dplyr::select(-fleetname)
      lbin_vector=sort(as.numeric(gsub('f', '', names(size.comp)[grep("f",names(size.comp))])))
      dat$lbin_vector=lbin_vector
      dat$N_lbins=length(dat$lbin_vector)
      dat$lencomp=size.comp%>%arrange(Fleet,Sex,year)
      dat$maximum_size=max(dat$maximum_size,max(lbin_vector))
    }
    if(!is.null(F.tagging) & !is.null(size.comp))
    {
      addfbit=dat$len_info[length(unique(F.tagging$fleet)),]
      rownames(addfbit)=rownames(Finfo)
      dat$len_info=rbind(dat$len_info,addfbit)
    }
    
    #conditional age at length
    ddumy=dat$age_info[rep(1,length(dis.flits)),]
    row.names(ddumy)=dis.flits
    ddumy=ddumy%>%
              rownames_to_column('fleetname')%>%
              mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                               ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                               fleetname)))
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
      dat$agebin_vector=seq(first.age,dat$Nages)
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
    
    #fleetinfo1
    addis=which(!dat$fleetinfo$fleetname%in%colnames(dat$fleetinfo1))
    if(length(addis)>0)
    {
      ddumy=dat$fleetinfo1[,which(colnames(dat$fleetinfo1)%in%dat$fleetinfo$fleetname)]
      ddumy1=dat$fleetinfo1[rep(1,length(addis))]
      colnames(ddumy1)=dat$fleetinfo$fleetname[addis]
      ddumy=cbind(ddumy,ddumy1)%>%
              relocate(dat$fleetinfo$fleetname)
      dat$fleetinfo1=ddumy
    }
    
    #fleetinfo2
    addis=which(!dat$fleetinfo$fleetname%in%colnames(dat$fleetinfo2))
    if(length(addis)>0)
    {
      ddumy=dat$fleetinfo2[,which(colnames(dat$fleetinfo2)%in%dat$fleetinfo$fleetname)]
      ddumy1=dat$fleetinfo2[rep(1,length(addis))]
      colnames(ddumy1)=dat$fleetinfo$fleetname[addis]
      ddumy=cbind(ddumy,ddumy1)%>%
              relocate(dat$fleetinfo$fleetname)
      dat$fleetinfo2=ddumy
    }
    
    #all these has 5 elements   
    if(nrow(WRL.fleet)>0)
    {
      dumifleetinfo=data.frame(dat$fleetinfo1[,1])
      colnames(dumifleetinfo)='WRL'
      dat$fleetinfo1=cbind(dat$fleetinfo1,dumifleetinfo)
      dumifleetinfo=data.frame(dat$fleetinfo2[,1])
      colnames(dumifleetinfo)='WRL'
      dat$fleetinfo2=cbind(dat$fleetinfo2,dumifleetinfo)
    }
    if(!is.null(F.tagging))
    {
      Nme=dat$fleetinfo$fleetname[grep('F.series',dat$fleetinfo$fleetname)]   
      dumifleetinfo=data.frame(dat$fleetinfo1[,1])
      colnames(dumifleetinfo)=Nme
      dat$fleetinfo1=cbind(dat$fleetinfo1,dumifleetinfo)
      dumifleetinfo=data.frame(dat$fleetinfo2[,1])
      colnames(dumifleetinfo)=Nme
      dat$fleetinfo2=cbind(dat$fleetinfo2,dumifleetinfo)
    }
    dat$fleetinfo1=dat$fleetinfo1[,match(dat$fleetinfo$fleetname,colnames(dat$fleetinfo1))]
    dat$fleetinfo2=dat$fleetinfo2[,match(dat$fleetinfo$fleetname,colnames(dat$fleetinfo2))]
    
    dat$surveytiming=as.numeric(dat$fleetinfo1[match('surveytiming',rownames(dat$fleetinfo1)),])
    dat$units_of_catch=as.numeric(dat$fleetinfo2[match('units',rownames(dat$fleetinfo2)),])
    dat$areas=as.numeric(dat$fleetinfo1[match('areas',rownames(dat$fleetinfo1)),])

    dat$comp_tail_compression=rep(-1,ncol(dat$fleetinfo1))
    dat$add_to_comp=rep(0.001,ncol(dat$fleetinfo1))
    dat$max_combined_lbin=rep(0,ncol(dat$fleetinfo1))
    
    
    #Tagging   
    if(Scenario$Tagging=='Yes' & !is.null(Tags))
    {
      dat$do_tags=1
      dat$N_tag_groups=nrow(Tags$releases)
      dat$N_recap_events=nrow(Tags$recaptures)
      dat$mixing_latency_period=Tags$mixing_latency_period    
      dat$max_periods=Tags$max_periods             
      dat$tag_releases=Tags$releases
      dat$tag_recaps=Tags$recaptures%>%arrange(Tag.group,Yr.rec)
    }
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
  
  #block patterns    
  #if(!life.history$Nblock_Patterns==0 & !is.null(abundance)) 
  if(!life.history$Nblock_Patterns==0)
  {
    ctl$N_Block_Designs=life.history$Nblock_Patterns
    ctl$blocks_per_pattern=life.history$blocks_per_pattern
    ctl$Block_Design=life.history$block_pattern_begin_end
    autgen.default=ctl$time_vary_auto_generation
    autgen.input=life.history$autogen
    names(autgen.input)=names(autgen.default)
    ctl$time_vary_auto_generation=autgen.input
  }
  
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
  Linf.multi=Scenario$alternative.Linf
  if(is.null(Linf.multi)) Linf.multi=1
  #females
  if("Eggs/kg_inter_Fem_GP_1"%in%rownames(ctl$MG_parms) & !"Eggs_alpha_Fem_GP_1"%in%rownames(ctl$MG_parms))
  {
    rownames(ctl$MG_parms)[match("Eggs/kg_inter_Fem_GP_1",rownames(ctl$MG_parms))]="Eggs_alpha_Fem_GP_1" 
    rownames(ctl$MG_parms)[match("Eggs/kg_slope_wt_Fem_GP_1",rownames(ctl$MG_parms))]="Eggs_beta_Fem_GP_1"
  }
  ctl$MG_parms["NatM_p_1_Fem_GP_1", c("INIT","PRIOR")]=rep(Scenario$Mmean,2)
  ctl$MG_parms["NatM_p_1_Fem_GP_1", c("LO","HI")]=c(ctl$MG_parms["NatM_p_1_Fem_GP_1", "INIT"]*.1,ctl$MG_parms["NatM_p_1_Fem_GP_1", "INIT"]*4)
  ctl$MG_parms["L_at_Amin_Fem_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Lzero*a_FL.to.TL+b_FL.to.TL),2)  #size for specified _Age(post-settlement)_for_L1
  ctl$MG_parms["L_at_Amin_Fem_GP_1", c("LO","HI")]=round(c(ctl$MG_parms["L_at_Amin_Fem_GP_1", "INIT"]*.25,ctl$MG_parms["L_at_Amin_Fem_GP_1", "INIT"]*2))
  ctl$MG_parms["L_at_Amax_Fem_GP_1", c("INIT","PRIOR")]=Linf.multi*rep(with(life.history,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL),2) #life.history$TLmax  #set at Linf as _Growth_Age_for_L2 was set at 999
  ctl$MG_parms["L_at_Amax_Fem_GP_1", c("LO","HI")]=round(c(ctl$MG_parms["L_at_Amax_Fem_GP_1", "INIT"]*.8,ctl$MG_parms["L_at_Amax_Fem_GP_1", "INIT"]*1.4))
  ctl$MG_parms["Wtlen_1_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$AwT,2)
  ctl$MG_parms["Wtlen_2_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$BwT,2)
  ctl$MG_parms["VonBert_K_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.F$k,2)
  ctl$MG_parms["VonBert_K_Fem_GP_1", c("LO","HI")]=c(life.history$Growth.F$k*.5,life.history$Growth.F$k*2)
  ctl$MG_parms["CV_young_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_young,2)
  ctl$MG_parms["CV_old_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_old,2)
  ctl$MG_parms["Mat50%_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$TL.mat.inf.slope[2],2) #life.history$TL.50.mat
  ctl$MG_parms["Mat_slope_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$TL.mat.inf.slope[1],2)
  min_Eggs_alpha=0
  if(fec.alpha<0) min_Eggs_alpha=fec.alpha*1.5
  ctl$MG_parms["Eggs_alpha_Fem_GP_1", c("LO","INIT","HI","PRIOR")]=c(min_Eggs_alpha,fec.alpha,100,fec.alpha)   
  ctl$MG_parms["Eggs_beta_Fem_GP_1", c("HI","INIT","PRIOR")]=c(fec.beta*1.5,rep(fec.beta,2))
  ctl$MG_parms["FracFemale_GP_1", c("INIT","PRIOR")]=rep(life.history$pup.sx.ratio,2)
  
  #males
  ctl$MG_parms["NatM_p_1_Mal_GP_1", c("INIT","PRIOR")]=rep(Scenario$Mmean,2)
  ctl$MG_parms["NatM_p_1_Mal_GP_1", c("LO","HI")]=c(ctl$MG_parms["NatM_p_1_Mal_GP_1", "INIT"]*.1,ctl$MG_parms["NatM_p_1_Mal_GP_1", "INIT"]*4)
  ctl$MG_parms["L_at_Amin_Mal_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Lzero*a_FL.to.TL+b_FL.to.TL),2)
  ctl$MG_parms["L_at_Amin_Mal_GP_1", c("LO","HI")]=ctl$MG_parms["L_at_Amin_Fem_GP_1", c("LO","HI")]
  ctl$MG_parms["L_at_Amax_Mal_GP_1", c("INIT","PRIOR")]=Linf.multi*rep(with(life.history,Growth.M$FL_inf*a_FL.to.TL+b_FL.to.TL),2)
  ctl$MG_parms["L_at_Amax_Mal_GP_1", c("LO","HI")]=round(c(ctl$MG_parms["L_at_Amax_Mal_GP_1", "INIT"]*.8,ctl$MG_parms["L_at_Amax_Mal_GP_1", "INIT"]*1.4))
  ctl$MG_parms["Wtlen_1_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$AwT.M,2)
  ctl$MG_parms["Wtlen_2_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$BwT.M,2)
  ctl$MG_parms["VonBert_K_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.M$k,2)
  ctl$MG_parms["VonBert_K_Mal_GP_1", c("LO","HI")]=c(life.history$Growth.M$k*.5,life.history$Growth.M$k*2)
  ctl$MG_parms["CV_young_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_young,2)
  ctl$MG_parms["CV_old_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_old,2)
  
  #estimate growth params
  if(Scenario$Model=='SS' & life.history$SS3.estim.growth.pars)  
  {
    if('Linf'%in%names(life.history$Growth.F.prior))
    {
      ctl$MG_parms["L_at_Amax_Fem_GP_1", "PHASE"]=2
      ctl$MG_parms["L_at_Amax_Fem_GP_1", c("PRIOR")]=life.history$Growth.F.prior$Linf
      ctl$MG_parms["L_at_Amax_Fem_GP_1", c("LO","HI")]=round(c(life.history$Growth.F.prior$Linf*.7,life.history$Growth.F.prior$Linf*1.3))
      ctl$MG_parms["L_at_Amax_Fem_GP_1", "PR_SD"]=life.history$Growth.F.prior$Linf.se  
      ctl$MG_parms["L_at_Amax_Fem_GP_1", "PR_type"]=life.history$Growth.prior.type$Linf
    }
    if('Linf'%in%names(life.history$Growth.M.prior))
    {
      ctl$MG_parms["L_at_Amax_Mal_GP_1", "PHASE"]=2
      ctl$MG_parms["L_at_Amax_Mal_GP_1", c("PRIOR")]=life.history$Growth.M.prior$Linf 
      ctl$MG_parms["L_at_Amax_Mal_GP_1", c("LO","HI")]=round(c(life.history$Growth.M.prior$Linf*.7,life.history$Growth.M.prior$Linf*1.3))
      ctl$MG_parms["L_at_Amax_Mal_GP_1", "PR_SD"]=life.history$Growth.M.prior$Linf.se   
      ctl$MG_parms["L_at_Amax_Mal_GP_1", "PR_type"]=life.history$Growth.prior.type$Linf 
    }
    
    if('k'%in%names(life.history$Growth.F.prior))
    {
      ctl$MG_parms["VonBert_K_Fem_GP_1", "PHASE"]=3
      ctl$MG_parms["VonBert_K_Fem_GP_1", c("PRIOR")]=life.history$Growth.F.prior$k
      ctl$MG_parms["VonBert_K_Fem_GP_1", c("LO","HI")]=round(c(life.history$Growth.F.prior$k*.60,life.history$Growth.F.prior$k*1.60),3)
      ctl$MG_parms["VonBert_K_Fem_GP_1", "PR_SD"]=life.history$Growth.F.prior$k.se 
      ctl$MG_parms["VonBert_K_Fem_GP_1", "PR_type"]=life.history$Growth.prior.type$k 
    }
    if('k'%in%names(life.history$Growth.M.prior))
    {
      ctl$MG_parms["VonBert_K_Mal_GP_1", "PHASE"]=3
      ctl$MG_parms["VonBert_K_Mal_GP_1", c("PRIOR")]=life.history$Growth.M.prior$k
      ctl$MG_parms["VonBert_K_Mal_GP_1", c("LO","HI")]=round(c(life.history$Growth.M.prior$k*.60,life.history$Growth.M.prior$k*1.60),3)
      ctl$MG_parms["VonBert_K_Mal_GP_1", "PR_SD"]=life.history$Growth.M.prior$k.se  
      ctl$MG_parms["VonBert_K_Mal_GP_1", "PR_type"]=life.history$Growth.prior.type$k
    }
  }
  
  #M at age
  if(Scenario$Model=='SS' & 'M.at.age'%in%names(Scenario))
  {
    ctl$natM_type=3 #0=1Parm; 1=N_breakpoints; 2=Lorenzen; 3=agespecific; 4=agespec_withseasinterpolate
    if(Scenario$M.at.age=="Mmean.mean.at.age")
    {
      M.age.fem=life.history$Mmean.mean.at.age[1:length(ageError)]
      M.age.male=life.history$Mmean.mean.at.age[1:length(ageError)]
    }
    if(Scenario$M.at.age=="constant") 
    {
      M.age.fem=rep(Scenario$Mmean,length(1:length(ageError))) 
      M.age.male=rep(Scenario$Mmean,length(1:length(ageError)))
    }
    if(is.na(M.age.fem[length(M.age.fem)])) M.age.fem[length(M.age.fem)]=M.age.fem[length(M.age.fem)-1]
    if(is.na(M.age.male[length(M.age.male)])) M.age.male[length(M.age.male)]=M.age.male[length(M.age.male)-1]
    natM=data.frame(matrix(c(M.age.fem,M.age.male),nrow=2,byrow = T))
    colnames(natM)=paste('Age',first.age:max(life.history$Max.age.F),sep='_')
    ctl$natM=natM
    ctl$MG_parms=ctl$MG_parms[-grep('NatM',rownames(ctl$MG_parms)),]
  }
  
  #recruitment pars  
  ctl$SR_function=Scenario$SR_type # 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
  ctl$SR_parms["SR_LN(R0)", c('LO','INIT','HI')]=with(Scenario,c(Ln_R0_min,Ln_R0_init,Ln_R0_max))
  ctl$SR_parms["SR_BH_steep", "INIT"]=Scenario$Steepness
  ctl$SR_parms["SR_BH_steep", "LO"]=Min.h.shark  #0.25
  if(Scenario$Steepness<Min.h.shark) ctl$SR_parms["SR_BH_steep", "LO"]=0.2
  if(Scenario$Model=='SS')  #estimate h with strong prior (Punt 2023)?
  {
    ctl$SR_parms["SR_BH_steep", "PHASE"]=life.history$Steepness_Phase
    ctl$SR_parms["SR_BH_steep", "PR_SD"]=Scenario$Steepness.sd  
    ctl$SR_parms["SR_BH_steep", "PR_type"]=1 #0, no prior; 1, symmetric beta; 2, beta; 3, lognormal; 4, lognormal with bias correction; 5, gamma; 6, normal
    if(ctl$SR_parms["SR_BH_steep", "PR_type"]==1) ctl$SR_parms["SR_BH_steep", "PR_SD"]=0.5  #increase CV to avoid bounds
    if(is.null(abundance)) ctl$SR_parms["SR_BH_steep", "PHASE"]=-4
  }
  if(Scenario$Model=='SS') SR_sigmaR=Scenario$SR_sigmaR  
  ctl$SR_parms["SR_sigmaR", c('LO','HI','INIT')]=c(.01,1,SR_sigmaR) #Spiny dogfish SS assessment
  if(ctl$SR_function==7)
  {
    id.h=which(rownames(ctl$SR_parms)=='SR_BH_steep')
    add.dumi.sr=ctl$SR_parms[rep(id.h,2),]
    rownames(add.dumi.sr)=c('SR_surv_zfrac','SR_surv_Beta')
    add.dumi.sr[rownames(add.dumi.sr)=='SR_surv_zfrac',c('LO','HI','INIT')]=c(0,1,Scenario$SR_surv_zfrac)
    add.dumi.sr[rownames(add.dumi.sr)=='SR_surv_Beta',c('LO','HI','INIT')]=c(0.2,ceiling(Scenario$SR_surv_Beta*1.5),Scenario$SR_surv_Beta)
    add.dumi.sr=add.dumi.sr%>%
                  mutate(PRIOR=INIT,
                         PR_type=0)
    ctl$SR_parms=ctl$SR_parms[-id.h,]
    ctl$SR_parms=rbind(ctl$SR_parms,add.dumi.sr)
    ctl$SR_parms=ctl$SR_parms[match(c("SR_LN(R0)","SR_surv_zfrac","SR_surv_Beta",
                                      "SR_sigmaR","SR_regime","SR_autocorr"),rownames(ctl$SR_parms)),]

  }
    #rec_devs
  if(Scenario$Model=='SSS') ctl$do_recdev=0  #do_recdev:  0=none; 1=devvector; 2=simple deviations
  if(Scenario$Model=='SS')
  {
    ctl$do_recdev=Scenario$do_recdev   # 0=none; 1=devvector (R=F(SSB)+dev). Sum to 0; 2=deviations (R=F(SSB)+dev). Does not sum to 0; 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
    RecDev_Phase=life.history$RecDev_Phase
    ctl$recdev_phase=RecDev_Phase
     
    if(is.null(abundance))
    {
      ctl$recdev_phase=1
      if(is.null(size.comp))
      {
           ctl$do_recdev=0
           ctl$recdev_phase=-1
      }
      if(!is.null(size.comp))   #don't calculate rec devs if only a few years of size comp and no abundance
      {
        if(nrow(size.comp)<3)
        {
          ctl$do_recdev=0
          ctl$recdev_phase=-1
        }
      }
    }
 
    if(RecDev_Phase<0)
    {
      ctl$recdev_early_start=0
      ctl$recdev_early_phase=-1
    }
      
    if(RecDev_Phase>0) 
    {
      ctl$recdev_early_start=life.history$recdev_early_start  
      ctl$recdev_early_phase=life.history$recdev_early_phase
    }
    ctl$max_bias_adj=0.8
    ctl$min_rec_dev=-1
    ctl$max_rec_dev=abs(ctl$min_rec_dev)
    ctl$MainRdevYrFirst=life.history$MainRdevYrFirst  
    if(is.null(size.comp) & !is.null(abundance))
    {
      ctl$MainRdevYrLast=max(abundance$Year)
      ctl$last_early_yr_nobias_adj=min(abundance$Year)-1
      ctl$first_yr_fullbias_adj=min(abundance$Year)+5 #first year with full bias adj. should be a few years into the data-rich period
    }
    if(!is.null(size.comp))
    {
      ctl$MainRdevYrLast=max(size.comp$year)
     if(!is.null(abundance) & life.history$Name=="smooth hammerhead") ctl$MainRdevYrLast=max(abundance$Year)#to allow fitting Southern2 cpue
      ctl$last_early_yr_nobias_adj=min(size.comp$year)-1
      ctl$first_yr_fullbias_adj=min(ctl$MainRdevYrLast-1,min(size.comp$year)+5) 
    }
    ctl$last_yr_fullbias_adj=endyr-2
    ctl$first_recent_yr_nobias_adj=endyr   #end_yr_for_ramp_in_MPD
    if(!is.na(life.history$last_early_yr_nobias_adj_in_MPD))
    {
      ctl$last_early_yr_nobias_adj=life.history$last_early_yr_nobias_adj_in_MPD
      ctl$first_yr_fullbias_adj=life.history$first_yr_fullbias_adj_in_MPD
      ctl$last_yr_fullbias_adj=life.history$last_yr_fullbias_adj_in_MPD
      ctl$first_recent_yr_nobias_adj=life.history$first_recent_yr_nobias_adj_in_MPD
      if(!is.null(abundance) & life.history$Name=="smooth hammerhead") ctl$first_recent_yr_nobias_adj=max(abundance$Year)#to allow fitting Southern2 cpue
      ctl$max_bias_adj=life.history$max_bias_adj_in_MPD
    }
  }
  
  #fishing mortality pars
  if(Scenario$use_F_ballpark)
  {
    ctl$F_ballpark=Scenario$F_ballpark
    ctl$F_ballpark_year=styr
  }
  
  # Get flit number
  if(Scenario$Model=='SS')
  {
    flit.numb=fleetinfo%>%
      dplyr::select(fleetname)%>%
      mutate(Fleet.n=row_number())
  }
  
  #Q pars
  if(Scenario$Model=='SS')
  {
    if(!is.null(abundance))
    {
      life.history$Q.inits=life.history$Q.inits%>%
                            filter(Fleet%in%fleetinfo$fleetname)%>%
                            dplyr::select(-Fleet.n)%>%
                            left_join(flit.numb,by=c('Fleet'='fleetname'))%>%
                            arrange(Fleet.n)%>%
                            relocate(Q,.after=Fleet.n)
      
      ddumy=ctl$Q_options[rep(1,length(dis.flits)),]%>%
                mutate(fleet=row_number())
      row.names(ddumy)=dis.flits
      ctl$Q_options=ddumy[-match(c('Northern.shark','Other'),rownames(ddumy)),]%>%
                      filter(fleet%in%unique(abundance$index))
      ctl$Q_options=ctl$Q_options[which(rownames(ctl$Q_options)%in%dat$fleetinfo$fleetname),]
      
      flits.with.cpue=flit.numb%>%filter(Fleet.n%in%unique(abundance$index))
      Nms=rownames(ctl$Q_parms)
      Nms=gsub(r"{\s*\([^\)]+\)}","",gsub("^.*?base_","",Nms))
      rownames(ctl$Q_parms)=Nms
      ddumy=ctl$Q_parms[rep(1,nrow(flits.with.cpue)),]
      rownames(ddumy)=flits.with.cpue$fleetname
      ddumy[grep('Southern.shark_2',rownames(ddumy)),'INIT']=-8
      ctl$Q_parms=ddumy
      #ctl$Q_parms=ctl$Q_parms[grep(paste(rownames(ctl$Q_options),collapse='|'),rownames(ctl$Q_parms)),]
       
      if(!'Northern.shark'%in%fleetinfo$fleetname)
      {
        ctl$Q_options$fleet=ctl$Q_options$fleet-1
      }
      if('Northern.shark'%in%fleetinfo$fleetname & !'Northern.shark'%in%rownames(ctl$Q_options))
      {
        addNSFQ=ctl$Q_options[1,]%>%
          mutate(fleet=match('Northern.shark',fleetinfo$fleetname))
        rownames(addNSFQ)="Northern.shark"
        ctl$Q_options=rbind(addNSFQ,ctl$Q_options)
      }
      if('Northern.shark'%in%fleetinfo$fleetname & !any(grepl('Northern.shark',rownames(ctl$Q_parms))))
      {
        addNSFQ=ctl$Q_parms[1,]%>%mutate(INIT=-7)
        rownames(addNSFQ)="Northern.shark"
        #rownames(addNSFQ)=paste('LnQ_base_',"Northern.shark(",match('Northern.shark',fleetinfo$fleetname),')',sep='')
        ctl$Q_parms=rbind(addNSFQ,ctl$Q_parms)
      }
      if('Other'%in%fleetinfo$fleetname & !'Other'%in%rownames(ctl$Q_options) & fleets[match('Other',fleetinfo$fleetname)]%in%unique(abundance$index))
      {
        addOther=ctl$Q_options[1,]%>%
          mutate(fleet=match('Other',fleetinfo$fleetname))
        rownames(addOther)="Other"
        ctl$Q_options=rbind(addOther,ctl$Q_options)
      }
      if('Other'%in%fleetinfo$fleetname & !any(grepl('Other',rownames(ctl$Q_parms))) & fleets[match('Other',fleetinfo$fleetname)]%in%unique(abundance$index))
      {
        addOther=ctl$Q_parms[1,]
        rownames(addOther)=paste('LnQ_base_',"Other(",match('Other',fleetinfo$fleetname),')',sep='')
        ctl$Q_parms=rbind(addOther,ctl$Q_parms)
      }
      if(!is.null(F.tagging))
      {
        if(length(grep('F.series',rownames(ctl$Q_options)))==0)
        {
          id.F.series=grep('F.series',dat$fleetinfo$fleetname)
          add.F.series.dummy=ctl$Q_options[1,]%>%mutate(fleet=id.F.series)
          rownames(add.F.series.dummy)=dat$fleetinfo$fleetname[id.F.series]
          ctl$Q_options=rbind(ctl$Q_options,add.F.series.dummy)
        }
        if(!any(grepl('F.series',rownames(ctl$Q_parms))))
        {
          add.F.series.dummy=ctl$Q_parms[grepl('Southern.shark_1',rownames( ctl$Q_parms)),]
          rownames(add.F.series.dummy)=paste('LnQ_base_',dat$fleetinfo$fleetname[id.F.series],"(",id.F.series,')',sep='')
          ctl$Q_parms=rbind(ctl$Q_parms,add.F.series.dummy)
        }
      }
      
      iis=sort(unique(dat$CPUE$index)) 
      Qdummi=left_join(ctl$Q_options%>%
                          tibble::rownames_to_column(),
                          fleetinfo%>%mutate(fleet1=row_number()),
                       by=c("rowname"="fleetname"))%>%
              mutate(fleet=fleet1)%>%
              tibble::column_to_rownames('rowname')%>%
              dplyr::select(names(ctl$Q_options))%>%
              filter(fleet%in%iis)
      if(nrow(Qdummi)>0)
      {
        ctl$Q_options=Qdummi
       # rownames(ctl$Q_options)=dat$fleetinfo$fleetname[iis]
      }
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
      
      #Add extra SD to Q if CV too small  
      Difalcivi=default.CV
      SS3.q.an.sol=life.history$SS3.q.an.sol
      if(Scenario$extra.SD.Q=='YES') Difalcivi=1
      n.indices=nrow(ctl$Q_options)
      Indx.small.CV=dat$CPUE%>%group_by(index)%>%summarise(Mean.CV=mean(CV))  
      if(life.history$Name%in%Extra_Q_species & !SS3.q.an.sol) Indx.small.CV$Mean.CV=Difalcivi*.9 #need extra_Q for to fit cpue
      Indx.small.CV=Indx.small.CV%>%filter(Mean.CV<Difalcivi)%>%pull(index)
      if(length(Indx.small.CV)>0)
      {
        if(is.data.frame(Indx.small.CV)) Indx.small.CV=Indx.small.CV$index
        ID.fleets.extraSD=match(Indx.small.CV,ctl$Q_options$fleet)
        ctl$Q_options$extra_se[ID.fleets.extraSD]=1
        rownames(ctl$Q_parms)=paste('fleet_',match(rownames(ctl$Q_parms),dat$fleetinfo$fleetname),sep='')
        Q_parms_estraSD=ctl$Q_parms[ID.fleets.extraSD,]%>%
          mutate(LO=0,
                 HI=1,
                 INIT=0.3,
                 PRIOR=0.3,
                 PHASE=3)
        rownames(Q_parms_estraSD)=paste(rownames(Q_parms_estraSD),"_Q_extraSD",sep='')
        ctl$Q_parms=rbind(ctl$Q_parms,Q_parms_estraSD)
        ctl$Q_parms=ctl$Q_parms[order(rownames(ctl$Q_parms)),]
      }
      
      
      #Block patterns
      if(!life.history$Nblock_Patterns==0 & ctl$time_vary_auto_generation[3]==0)
      {
        Q_param_Block=life.history$Q_param_Block
        Q_param_Blk_Fxn=life.history$Q_param_Blk_Fxn
        id_qs=grep('Southern.shark_1',rownames(ctl$Q_parms))
        if(length(id_qs)>1)
        {
          Q_param_Block=rep(Q_param_Block,length(id_qs))
          Q_param_Blk_Fxn=rep(Q_param_Blk_Fxn,length(id_qs))
        }
        ctl$Q_parms[id_qs,'Block']=Q_param_Block
        ctl$Q_parms[id_qs,'Block_Fxn']=Q_param_Blk_Fxn
      }
      
      #Analytical solution
      if(SS3.q.an.sol)
      {
        ctl$Q_options$float=1
      }
    }
    if(is.null(abundance))
    {
      ctl$Q_options=NULL
      ctl$Q_parms=NULL
    }
  }
  
  #selectivity pars  
  if(Scenario$Model=='SS') 
  {
    #1. size_selex   
    ddumy=ctl$size_selex_types[rep(1,length(dis.flits)),]%>%
              mutate(Fleet=row_number())%>%
              left_join(flit.numb,by=c('Fleet'='Fleet.n'))%>%
              mutate(Pattern=case_when(fleetname=='Other'~15,
                                       grepl('Southern.shark_2',fleetname)~15,
                                       TRUE~Pattern),
                     Special=case_when(fleetname=='Other'~match('Northern.shark',fleetname),
                                       fleetname=='Southern.shark_2'~match('Southern.shark_1',fleetname),
                                       fleetname=='Southern.shark_2_West'~match('Southern.shark_1_West',fleetname),
                                       fleetname=='Southern.shark_2_Zone1'~match('Southern.shark_1_Zone1',fleetname),
                                       fleetname=='Southern.shark_2_Zone2'~match('Southern.shark_1_Zone2',fleetname),
                                       TRUE~Special))%>%
              dplyr::select(-fleetname)
    row.names(ddumy)=dis.flits
    
    if(!'Northern.shark'%in%rownames(ddumy))   
    {
      ddumy=ddumy%>%
        mutate(Pattern=ifelse(fleetname=='Other',24,Pattern),
               Special=case_when(fleetname=='Other'~0,
                                 fleetname=='Southern.shark_2'~2,
                                 TRUE~Special))
    }
    if(!is.null(F.tagging))
    {
      add.F.series.dummy=ddumy[grep('Southern.shark_2',ddumy$fleetname),]%>%
        mutate(Fleet=id.F.series,
               fleetname=dat$fleetinfo$fleetname[id.F.series])
      rownames(add.F.series.dummy)=dat$fleetinfo$fleetname[id.F.series]
      ddumy=rbind(ddumy,add.F.series.dummy)
    }
    
      #replace sels that mimic other fleet
    if(!is.null(life.history$SS_selectivity_mimic))
    {
      id.FleeT=match(life.history$SS_selectivity_mimic$Fleet,rownames(ddumy))
      id.FleeT=id.FleeT[!is.na(id.FleeT)]
      if(length(id.FleeT)>0)
      {
        id.FleeT.mimic=match(life.history$SS_selectivity_mimic$Fleet.mimic,rownames(ddumy))
        ddumy[id.FleeT,'Pattern']=15
        ddumy[id.FleeT,'Special']=ddumy$Fleet[id.FleeT.mimic]
      }
    }
    
    #select relevant selectivity fleets
    life.history$SS_selectivity=life.history$SS_selectivity%>%
                                      filter(Fleet%in%fleetinfo$fleetname)%>%
                                      left_join(flit.numb,by=c('Fleet'='fleetname'))%>%
                                      arrange(Fleet.n)%>%
                                      dplyr::select(-Fleet.n)
     life.history$SS_selectivity_phase=life.history$SS_selectivity_phase%>%
                                            filter(Fleet%in%fleetinfo$fleetname)%>%
                                            left_join(flit.numb,by=c('Fleet'='fleetname'))%>%
                                            arrange(Fleet.n)%>%
                                            dplyr::select(-Fleet.n)
     if(!is.null(life.history$SS_selectivity.sensitivity))
     {
       life.history$SS_selectivity.sensitivity=life.history$SS_selectivity.sensitivity%>%
                                                 filter(Fleet%in%fleetinfo$fleetname)%>%
                                                 left_join(flit.numb,by=c('Fleet'='fleetname'))%>%
                                                 arrange(Fleet.n)%>%
                                                 dplyr::select(-Fleet.n)
       
       life.history$SS_selectivity.sensitivity_phase=life.history$SS_selectivity.sensitivity_phase%>%
                                                       filter(Fleet%in%fleetinfo$fleetname)%>%
                                                       left_join(flit.numb,by=c('Fleet'='fleetname'))%>%
                                                       arrange(Fleet.n)%>%
                                                       dplyr::select(-Fleet.n)
       
     }
     
     #select if using sensitivity selectivity   
    if(!is.na(Scenario$NSF.selectivity))
    {
      life.history$SS_selectivity=life.history$SS_selectivity.sensitivity
      life.history$SS_selectivity_phase=life.history$SS_selectivity.sensitivity_phase
    }
    
    #change selectivity pattern to logistic if specified in SS_selectivity 
    Logis.sel=life.history$SS_selectivity%>%filter(Fleet=='Northern.shark')
    if('Northern.shark'%in%flitinfo$fleetname & all(is.na(Logis.sel[,c('P_3','P_4','P_5','P_6')])))
    {
      ddumy[rownames(ddumy)=='Northern.shark','Pattern']=1
    }
    Logis.sel=life.history$SS_selectivity%>%filter(Fleet=='Other')
    if(!'Northern.shark'%in%dis.flits & all(is.na(Logis.sel[,c('P_3','P_4','P_5','P_6')])))
    {
      ddumy[rownames(ddumy)=='Other','Pattern']=1
    }
    Logis.sel=life.history$SS_selectivity%>%filter(Fleet=='Survey')
    if(all(is.na(Logis.sel[,c('P_3','P_4','P_5','P_6')])))
    {
      ddumy[rownames(ddumy)=='Survey','Pattern']=1
    }
    
    #Add WRL selectivity if appropriate 
    WRL.sel=life.history$SS_selectivity%>%filter(Fleet=='WRL')
    if(nrow(WRL.sel)>0)  
    {
      ddumy.sel=ddumy[rownames(ddumy)=='Survey',]%>%
                mutate(Fleet=dat$CPUEinfo[match('Survey',rownames(dat$CPUEinfo)),'Fleet'])
      ddumy=ddumy[!rownames(ddumy)=='Survey',]
      if(all(!is.na(WRL.sel[,c('P_1','P_2','P_3','P_4','P_5','P_6')]))) Ptrn=24
      if(all(is.na(WRL.sel[,c('P_3','P_4','P_5','P_6')]))) Ptrn=1
      dumi.lbster=data.frame(fleetname='WRL', Pattern=Ptrn, Discard=0, Male=0, Special=0,
                             Fleet=dat$CPUEinfo[match('WRL',rownames(dat$CPUEinfo)),'Fleet'])
      row.names(dumi.lbster)='WRL'
      ddumy=rbind(ddumy,dumi.lbster,ddumy.sel)%>%
            arrange(Fleet)  
      
      if(!is.null(dat$len_info))
      {
        ddummy=dat$len_info[1,] 
        rownames(ddummy)='WRL'
        if("Survey"%in%rownames(dat$len_info))
        {
          dat$len_info=rbind(dat$len_info[-match("Survey",rownames(dat$len_info)),],ddummy,dat$len_info[match("Survey",rownames(dat$len_info)),])
        }else
        {
          dat$len_info=rbind(dat$len_info,ddummy)
        }
      }

      ddummy=dat$age_info[1,]
      rownames(ddummy)='WRL'
      if("Survey"%in%rownames(dat$age_info))
      {
        dat$age_info=rbind(dat$age_info[-match("Survey",rownames(dat$age_info)),],ddummy,dat$age_info[match("Survey",rownames(dat$age_info)),]) 
      }else
      {
        dat$age_info=rbind(dat$age_info,ddummy)
      }
      
      if(!is.null(dat$len_info)) dat$len_info=dat$len_info[match(rownames(dat$CPUEinfo),rownames(dat$len_info)),]
      dat$age_info=dat$age_info[match(rownames(dat$CPUEinfo),rownames(dat$age_info)),]
    }
    
    #sawsharks other fleet
    if(life.history$Name=="sawsharks")
    {
      ddumy=ddumy%>%
              mutate(Pattern=ifelse(fleetname=='Other',24,Pattern),
                     Special=ifelse(fleetname=='Other',0,Special))
    }
    
    #turn on Southern2 if mean.weight of length comp data  
    if(any(grepl('Southern.shark_2',rownames(ddumy))))
    {
      id.Southern.shark2.fleet=grep('Southern.shark_2',fleetinfo$fleetname)
      id.mean.body.w=unique(dat$meanbodywt$fleet)
      id.length.comp=unique(dat$lencomp$Fleet)
      turn.on=id.Southern.shark2.fleet[which(id.Southern.shark2.fleet%in%unique(c(id.mean.body.w,id.length.comp)))]
      if(length(turn.on)>0)
      {
        ddumy=ddumy%>%
              mutate(Pattern=ifelse(Fleet%in%turn.on & Pattern==15,24,Pattern),
                     Special=ifelse(Fleet%in%turn.on & Special>0,0,Special))
      }
    }

    
    #replace size_selex_types
    ctl$size_selex_types=ddumy%>%dplyr::select(-Fleet)  
    
    Sel.ptrn=ctl$size_selex_types$Pattern
    names(Sel.ptrn)=paste('Fishery',1:length(Sel.ptrn),sep='')
    
    #replace size_selex_parms 
    ddami=ctl$size_selex_parms%>%
            mutate(Fleet=gsub("\\s*\\([^\\)]+\\)","",str_after_nth(rownames(ctl$size_selex_parms), "_", 2)))
    dis.vec=paste0(1:6,'_Southern.shark_1')
    Zones.condition=length(grep(paste(c("West","Zone1","Zone2"),collapse='|'),fleetinfo$fleetname))>0
    Souther2.condition=FALSE
    if(!Zones.condition)
    {
      Souther2.condition=ddumy[grep('Southern.shark_2',rownames(ddumy)),'Pattern']
      Souther2.condition=!Souther2.condition==15
    }

    if(Zones.condition | Souther2.condition)
    {
      if(Zones.condition)
      {
        y.vec=c('Southern.shark_1_West','Southern.shark_1_Zone1','Southern.shark_1_Zone2',
                'Southern.shark_2_West','Southern.shark_2_Zone1','Southern.shark_2_Zone2')
      }
      if(Souther2.condition)
      {
        y.vec='Southern.shark_2'
      }
      ddami_extra=vector('list',length(dis.vec))
      for(v in 1:length(dis.vec))
      {
        aaa=fn.add.fleet.zone.sel(d=ddami,
                                  x=dis.vec[v],
                                  y=y.vec)
        rownames(aaa)=paste0('SizeSel_P_',v,"_",aaa$Fleet)  
        ddami_extra[[v]]=aaa%>%mutate(Fleet=paste0(v,'_',Fleet))
      }
      ddami_extra=do.call(rbind,ddami_extra)%>%arrange(Fleet)
      ddami=rbind(ddami,ddami_extra)%>%
              mutate(Fleet=sub(".*?_","",Fleet))%>%
              filter(Fleet%in%fleetinfo$fleetname)%>%
              arrange(Fleet)%>%
              dplyr::select(-Fleet)
    }
    ctl$size_selex_parms=ddami 
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
    if("Other"%in%dis.flits & life.history$Name=="sawsharks")
    {
      other.sel.patrn=chosen.sel.patrns[1:6,]
      rownames(other.sel.patrn)=str_replace(rownames(other.sel.patrn)[grep('Northern.shark',rownames(other.sel.patrn))], "Northern.shark", "Other")
      chosen.sel.patrns=rbind(chosen.sel.patrns,other.sel.patrn)
      row_nm_size_selex_parms=c(row_nm_size_selex_parms,rep("Other",6))
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
    
    #turn on Southern.shark_2 if size compo data OR if meanbodywt and specified in life.history
    supersid.this=FALSE
    if(length(grep(paste(c("West","Zone1","Zone2"),collapse='|'),fleetinfo$fleetname))>0) supersid.this=TRUE
    if(!supersid.this)
    {
      if(!is.null(size.comp) | isTRUE(life.history$fit.Southern.shark_2.to.meanbodywt))
      {
        if(!is.null(size.comp))
        {
          Tab.size.comp.dat=with(dat$lencomp,table(Fleet))
          names(Tab.size.comp.dat)=fleetinfo$fleetname[as.numeric(names(Tab.size.comp.dat))]
          dis.vec=c("Southern.shark_2","Southern.shark_2_West","Southern.shark_2_Zone1","Southern.shark_2_Zone2")
          names(dis.vec)=c("Southern.shark_1","Southern.shark_1_West","Southern.shark_1_Zone1","Southern.shark_1_Zone2")
          for(v in 1:length(dis.vec))
          {
            input=dis.vec[v]
            output=names(dis.vec[v])  
            nn.dis=subset(Tab.size.comp.dat,names(Tab.size.comp.dat)==input)
            if(length(nn.dis)>0)
            {
              ctl$size_selex_types[rownames(ctl$size_selex_types)==input,]=ctl$size_selex_types[rownames(ctl$size_selex_types)==output,]
              
              add.Southern.shark_2.pars=ctl$size_selex_parms[grepl(output,rownames(ctl$size_selex_parms)),]
              rownames(add.Southern.shark_2.pars)=str_replace(rownames(add.Southern.shark_2.pars), "k_1", "k_2")
              
              ctl$size_selex_parms=rbind(ctl$size_selex_parms,add.Southern.shark_2.pars)
              ctl$size_selex_parms$fleet=sub(".*?_","",str_remove_all(rownames(ctl$size_selex_parms), "SizeSel_P_"))
              ctl$size_selex_parms$order=gsub("\\_.*", "", str_remove_all(rownames(ctl$size_selex_parms), "SizeSel_P_"))
              ctl$size_selex_parms=ctl$size_selex_parms%>%
                left_join(data.frame(fleet=rownames(ctl$size_selex_types%>%filter(Special==0)))%>%
                            mutate(Fleet.order=row_number()),
                          by='fleet')%>%
                arrange(Fleet.order,order)
              rownames(ctl$size_selex_parms)=paste0('SizeSel_P_',ctl$size_selex_parms$order,'_',ctl$size_selex_parms$fleet)
              ctl$size_selex_parms=ctl$size_selex_parms%>%
                dplyr::select(-order,-fleet,-Fleet.order)
            }
          }
          
        }
        if( !any(grepl('Southern.shark_2',rownames(ctl$size_selex_parms))) &  (!is.null(meanbodywt) & isTRUE(life.history$fit.Southern.shark_2.to.meanbodywt)))
        {
          ctl$size_selex_types[rownames(ctl$size_selex_types)=="Southern.shark_2",]=ctl$size_selex_types[rownames(ctl$size_selex_types)=="Southern.shark_1",]
          add.Southern.shark_2.pars=ctl$size_selex_parms[grepl('Southern.shark_1',rownames(ctl$size_selex_parms)),]
          rownames(add.Southern.shark_2.pars)=str_replace(rownames(add.Southern.shark_2.pars), "k_1", "k_2")
          ctl$size_selex_parms=rbind(ctl$size_selex_parms,add.Southern.shark_2.pars)
          ctl$size_selex_parms$fleet=gsub("^\\.","",str_remove_all(rownames(ctl$size_selex_parms), paste(c("SizeSel_P_", paste0(1:6,"_")), collapse = "|")))
          ctl$size_selex_parms$order=gsub("^\\.","",str_remove_all(rownames(ctl$size_selex_parms), paste(c("SizeSel_P_",paste0("_",ctl$size_selex_parms$fleet)), collapse = "|")))
          ctl$size_selex_parms=ctl$size_selex_parms%>%
            left_join(data.frame(fleet=rownames(ctl$size_selex_types%>%filter(Special==0)))%>%mutate(Fleet.order=row_number()),
                      by='fleet')%>%
            arrange(Fleet.order,order)
          rownames(ctl$size_selex_parms)=paste0('SizeSel_P_',ctl$size_selex_parms$order,'_',ctl$size_selex_parms$fleet)
          ctl$size_selex_parms=ctl$size_selex_parms%>%
            dplyr::select(-order,-fleet,-Fleet.order)
        }
      } 
    }

    #add WRL  
    if(nrow(WRL.sel)>0)
    {
      add.this=ctl$size_selex_parms[which(row_nm_size_selex_parms%in%dis.flits),]
      add.WRL=add.this[which(row_nm_size_selex_parms==dis.flits[1]),]
      rownames(add.WRL)=str_replace(rownames(add.WRL), dis.flits[1], "WRL")
      if(Ptrn==1) add.WRL=add.WRL[1:2,]
      if(any(!rownames(add.WRL)%in%rownames(ctl$size_selex_parms)))
      {
        ctl$size_selex_parms=rbind(ctl$size_selex_parms,add.WRL)
        ctl$size_selex_parms$fleet=gsub("^\\.","",str_remove_all(rownames(ctl$size_selex_parms), paste(c("SizeSel_P_", paste0(1:6,"_")), collapse = "|")))
        ctl$size_selex_parms$order=gsub("^\\.","",str_remove_all(rownames(ctl$size_selex_parms), paste(c("SizeSel_P_",paste0("_",ctl$size_selex_parms$fleet)), collapse = "|")))
        ctl$size_selex_parms=ctl$size_selex_parms%>%
          left_join(data.frame(fleet=rownames(ctl$size_selex_types%>%filter(Special==0)))%>%mutate(Fleet.order=row_number()),
                    by='fleet')%>%
          arrange(Fleet.order,order)
        rownames(ctl$size_selex_parms)=paste0('SizeSel_P_',ctl$size_selex_parms$order,'_',ctl$size_selex_parms$fleet)
        ctl$size_selex_parms=ctl$size_selex_parms%>%
          dplyr::select(-order,-fleet,-Fleet.order)
      }
    }
    
    #allocated species specific values to sel pars 
    Mirrored.sels=rownames(ctl$size_selex_types%>%filter(Pattern==15))
    if(length(Mirrored.sels)>0)
    {
      life.history$SS_selectivity=life.history$SS_selectivity%>%filter(!Fleet%in%Mirrored.sels)
      id.drop=grep(paste(Mirrored.sels,collapse='|'),rownames(ctl$size_selex_parms))
      if(length(id.drop)>0) ctl$size_selex_parms=ctl$size_selex_parms[-id.drop,]
    }
    id.fleets=fn.get.in.betwee(x=rownames(ctl$size_selex_parms))
    pis=unique(id.fleets)  
    for(px in 1:length(pis))
    {
      iid=which(id.fleets==pis[px])
      this.par=life.history$SS_selectivity[,match(pis[px],colnames(life.history$SS_selectivity))]  
      names(this.par)=life.history$SS_selectivity$Fleet
      this.par=subset(this.par,!is.na(this.par))
      if(length(this.par)>0)
      {
        ctl$size_selex_parms[iid,"INIT"]=this.par[match( str_remove_all(rownames(ctl$size_selex_parms[iid,]),paste0("SizeSel_", pis[px], '_')),names(this.par))]
        enen=length(ctl$size_selex_parms[iid,"INIT"])
        multiplr=rep(0.1,enen)
        multiplr[which(!is.na(ctl$size_selex_parms[iid,"INIT"]))]=ifelse(this.par<0,2,multiplr)
        low.bound=multiplr*ctl$size_selex_parms[iid,"INIT"]
        if(pis[px]=="P_1")
        {
          #bump=1.05 #1.25
          #low.bound=sapply(low.bound, function(x) max(dat$minimum_size*bump,x))
          bump=seq(dat$minimum_size,dat$minimum_size+dat$binwidth,by=dat$binwidth)
          bump=bump[length(bump)]+(dat$binwidth/2)
          low.bound=sapply(low.bound, function(x) max(bump,x))
          if(!is.null(size.comp))
          {
            low.bound=sapply(low.bound, function(x) max(min(dat$lbin_vector),x))  
            #if(life.history$Name%in%c("whiskery shark","tiger shark")) low.bound=min(dat$lbin_vector) 
          }
          #low.bound=min(low.bound,mean(seq(dat$minimum_size,(dat$minimum_size+1*TL.bins.cm),by=TL.bins.cm)))
        }
        ctl$size_selex_parms[iid,"LO"]=low.bound 
        
        multiplr=rep(1.3,enen)
        multiplr[which(!is.na(ctl$size_selex_parms[iid,"INIT"]))]=ifelse(this.par<0,-2,multiplr)
        up.bound=multiplr*ctl$size_selex_parms[iid,"INIT"]
        if(pis[px]=="P_1")
        {
          up.bound=sapply(up.bound, function(x) min(dat$maximum_size*.975,x))
        }
        ctl$size_selex_parms[iid,"HI"]=up.bound
      }else
      {
        ctl$size_selex_parms[iid,]=NA
      }
    }
    ctl$size_selex_parms=ctl$size_selex_parms%>%filter(!is.na(INIT))
    
    #set phases for estimable selectivities  
    if(!is.null(size.comp)|!is.null(meanbodywt))  
    {
      flit.size.comp.obs=sort(unique(c(meanbodywt$fleet,dat$lencomp$Fleet)))
      flit.no.size.comp.obs=fleetinfo%>%filter(!fleetname%in%rownames(dat$len_info)[flit.size.comp.obs])%>%pull(fleetname)
      # life.history$SS_selectivity_phase=life.history$SS_selectivity_phase%>%
      #                                     filter(Fleet%in%life.history$SS_selectivity$Fleet)
      no.length.comp=which(life.history$SS_selectivity_phase$Fleet%in%flit.no.size.comp.obs)
      if(length(no.length.comp)>0)
      {
        for(n in 1:length(no.length.comp))
        {
          life.history$SS_selectivity_phase[no.length.comp[n],-1]=-2
        }
      }


      Sel_phase=life.history$SS_selectivity_phase%>%
                  gather(P,PHASE1,-Fleet)%>%
                  mutate(dumy=paste('SizeSel',P,Fleet,sep="_"))
      
      xx=ctl$size_selex_parms%>%
                tibble::rownames_to_column(var = "dumy")%>%
                left_join(Sel_phase%>%dplyr::select(PHASE1,dumy),by='dumy')%>%
                mutate(PHASE=PHASE1)%>%
                dplyr::select(-PHASE1)%>%
                `row.names<-`(., NULL)%>%
                column_to_rownames(var = "dumy")
      ctl$size_selex_parms=xx
    }
    if(is.null(size.comp))ctl$size_selex_parms[,"PHASE"]=-2
    if(!is.null(size.comp))   
    {
      if(length(Mirrored.sels)>0) flit.no.size.comp.obs=subset(flit.no.size.comp.obs,!flit.no.size.comp.obs%in%Mirrored.sels)
      if("Southern.shark_2"%in%flit.no.size.comp.obs & isTRUE(life.history$fit.Southern.shark_2.to.meanbodywt))
      {
        flit.no.size.comp.obs=subset(flit.no.size.comp.obs,!flit.no.size.comp.obs=="Southern.shark_2") 
      }
      if(length(flit.no.size.comp.obs)>0)
      {
        for(px in 1:length(flit.no.size.comp.obs))
        {
          iid=grep(flit.no.size.comp.obs[px],rownames(ctl$size_selex_parms))
          ctl$size_selex_parms[iid,]$PHASE=-2
        }
      }
    }
    # if(isTRUE(life.history$fit.Southern.shark_2.to.meanbodywt) & !is.null(meanbodywt))
    # {
    #   id.South2=match(paste('SizeSel',colnames(life.history$SS_selectivity_phase)[-1],'Southern.shark_2',sep='_'),rownames(ctl$size_selex_parms))
    #   ctl$size_selex_parms[id.South2,'PHASE']=unlist(life.history$SS_selectivity_phase%>%filter(Fleet=='Southern.shark_2')%>%dplyr::select(-Fleet))
    #   
    # }
    # if(any(is.null(abundance) | life.history$Name=="spinner shark"))   #fixed most sel pars if no abundance (SS-CL)  
    # {
    #   if(life.history$Name=="spinner shark")
    #   {
    #     ctl$size_selex_parms$PHASE=-abs(ctl$size_selex_parms$PHASE)
    #     ctl$size_selex_parms$PHASE[grep('P_3_Southern.shark_1',rownames(ctl$size_selex_parms))]=2
    #   }
    # }
    # if(life.history$Name=="whiskery shark")
    # {
    #   ctl$size_selex_parms$PHASE[grep('P_2_Southern.shark',rownames(ctl$size_selex_parms))]=-4
    # }
    
      #turn on Southern.shark_1 if available meanbodywt (because Southern.shark_2 mirrors Southern.shark_1) 
    if(!is.null(meanbodywt) & !is.null(size.comp))
    {
      
      xx=size.comp%>%filter(Fleet%in%grep('Southern.shark',fleetinfo$fleetname))
      if(nrow(xx)==0)   #only need to turn off if there is no size comp for Southern.shark_1 or Southern.shark_2
      {
        iiD=grepl('P_1_Southern.shark_1',rownames(ctl$size_selex_parms))
        ctl$size_selex_parms[iiD,"PHASE"]=abs(ctl$size_selex_parms[iiD,"PHASE"])
        ctl$size_selex_parms[iiD,'PRIOR']=ctl$size_selex_parms[iiD,c('INIT')]
        ctl$size_selex_parms[iiD,'PR_SD']=ctl$size_selex_parms[iiD,c('INIT')]*.2
        ctl$size_selex_parms[iiD,'PR_type']=0 #0, no prior; 1, symmetric beta; 2, beta; 3, lognormal; 4, lognormal with bias correction; 5, gamma; 6, normal
        if(life.history$Name=="tiger shark")
        {
          ctl$size_selex_parms[iiD,'PR_type']=0
          iiD=grepl(paste(c('P_4_Southern.shark_1'),collapse='|'),rownames(ctl$size_selex_parms))
          ctl$size_selex_parms[iiD,"PHASE"]=abs(ctl$size_selex_parms[iiD,"PHASE"])
        }
      }
    }
    ctl$size_selex_parms[grep('Southern.shark_2',rownames(ctl$size_selex_parms)),'PHASE']=-2
    if(!is.null(size.comp) | !is.null(meanbodywt))
    {
      dis.daily.size.comp.flit=NA
      if(!is.null(size.comp))
      {
        daily.fleets=grep('Southern.shark_2',fleetinfo$fleetname)
        daily.size.comp=size.comp%>%filter(year>2005)%>%filter(Fleet%in%daily.fleets)
        if(nrow(daily.size.comp)>0) dis.daily.size.comp.flit=unique(daily.size.comp$Fleet)
      }
      dis.mean.w.flit=NA
      if(!is.null(meanbodywt)) dis.mean.w.flit=fleetinfo$fleetname[unique(meanbodywt$fleet)]
      
      id.fleet.to.estim=unique(c(dis.daily.size.comp.flit,dis.mean.w.flit))  
      id.fleet.to.estim=id.fleet.to.estim[!is.na(id.fleet.to.estim)]
      if(length(id.fleet.to.estim)>0)
      {
        id.fleet.to.estim=paste0(c('SizeSel_P_1_','SizeSel_P_3_'),id.fleet.to.estim)
        id.fleet.to.estim=paste(id.fleet.to.estim,collapse='|')
        ctl$size_selex_parms[grep(id.fleet.to.estim,rownames(ctl$size_selex_parms)),'PHASE']=3
      }
    }
    
    #set NSF and Survey to logistic if specified in scenario 
    if(!is.na(Scenario$NSF.selectivity))
    {
      if(Scenario$NSF.selectivity=='Logistic')
      {
        idflits=match(c('Northern.shark','Survey'),rownames(ctl$size_selex_types))
        idflits=subset(idflits,!is.na(idflits))
        ctl$size_selex_types[idflits,'Pattern']=1
        idselpars=match(c('SizeSel_P_3_Northern.shark','SizeSel_P_4_Northern.shark',
                 'SizeSel_P_5_Northern.shark','SizeSel_P_6_Northern.shark',
                 'SizeSel_P_3_Survey','SizeSel_P_4_Survey',
                 'SizeSel_P_5_Survey','SizeSel_P_6_Survey'),rownames(ctl$size_selex_parms))
        idselpars=subset(idselpars,!is.na(idselpars))
        if(length(idselpars)>0)ctl$size_selex_parms=ctl$size_selex_parms[-idselpars,]
      }   
    }
    
    #Retention and discard mortality  
    if('SS_retention'%in%names(life.history))
    {
        #retention
      xx=life.history$SS_retention%>%
        select_if(~ !any(is.na(.)))
      n.disc.pars=length(grep('P',colnames(xx)))
      if(n.disc.pars<4) Discard_option=1 #0=none; 1= retention and all discards Dead;
      if(n.disc.pars==4) Discard_option=2 #2= 4 retention parameters (i.e. logistic) and 4 discard mortality parameters
      if(n.disc.pars==7) Discard_option=4 #4= 7 retention parameters (i.e. dome=shaped) and 4 discard mortality parameters
      id_disc=grep(xx$Fleet,rownames(ctl$size_selex_types))
      if(length(id_disc)>1)
      {
        xx=xx %>%
          slice(rep(1:n(), each = length(id_disc)))%>%
          mutate(Fleet=rownames(ctl$size_selex_types)[id_disc])
      }
      colnames.xx=names(xx)[-1]
      ctl$size_selex_types$Discard[id_disc]=Discard_option 
      
      retention_params=ctl$size_selex_parms[grep(paste(xx$Fleet,collapse='|'),rownames(ctl$size_selex_parms)),]
      retention_params=retention_params%>%
                          mutate(rowname=rownames(retention_params),
                                 Fleet=str_after_nth(rowname,'_',3),
                                 P=str_before_nth(str_after_nth(rowname, "_", 1),"_",2))
     adiss=colnames.xx[which(!colnames.xx%in%retention_params$P)]
     if(length(adiss)>0)
     {
       adiss=rep(adiss,nrow(xx))
       ddmI=retention_params[1:length(adiss),]%>%
         mutate(Fleet=xx$Fleet,   
                P=adiss)
       rownames(ddmI)=paste0('SizeSel_',adiss,'_',xx$Fleet)
       retention_params=rbind(retention_params,ddmI)
     }
     retention_params=retention_params%>%
                        arrange(Fleet,P)
     retention_params=retention_params%>%dplyr::select(-c(rowname,Fleet,P))
     retention_params$dumi=paste(rep(paste(grep(paste(xx$Fleet,collapse='|'),fleetinfo$fleetname),2,sep='_'),each=ncol(xx)-1),
                                 rep(colnames.xx,nrow(xx)),sep='_')
      retention_params=retention_params%>%    
                                  mutate(INIT=c(t(xx[,2:(ncol(xx))])),
                                         PRIOR=INIT)
      this.par=retention_params$INIT
      multiplr=rep(0.1,length(this.par))
      multiplr=ifelse(this.par<0,2,multiplr)
      low.bound=multiplr*retention_params$INIT
      id0=which(this.par==0)
      id999=which(this.par==999)
      if(length(id0)>0) low.bound[id0]=-1
      if(length(id999)>0) low.bound[id999]=-10
      multiplr=rep(1.5,length(this.par))
      multiplr=ifelse(this.par<0,-2,multiplr)
      up.bound=multiplr*retention_params$INIT
      if(length(id0)>0) up.bound[id0]=1
      if(length(id999)>0) up.bound[id999]=1e3
      retention_params=retention_params%>%
                        mutate(LO=low.bound,
                               HI=up.bound,
                               PHASE=-2,
                               PR_SD=99)%>%
                          replace(is.na(.), 0)
      rownames(retention_params)=str_replace(rownames(retention_params),'SizeSel','Retention')
      
      
      #discard mortality 
      xx=life.history$SS_discard_mortality
      if(length(id_disc)>1)
      {
        xx=xx %>%
          slice(rep(1:n(), each = length(id_disc)))%>%
          mutate(Fleet=rownames(ctl$size_selex_types)[id_disc])
      }
      colnames.xx=names(xx)[-1]
      discard_mortality=ctl$size_selex_parms[grep(paste(xx$Fleet,collapse='|'),rownames(ctl$size_selex_parms)),]
      discard_mortality=discard_mortality%>%
                          mutate(rowname=rownames(discard_mortality),
                                 P=str_before_nth(str_after_nth(rowname, "_", 1),"_",2))%>%
                          filter(P%in%colnames.xx)%>%
                          dplyr::select(-rowname,-P)%>%  
                          mutate(INIT=c(t(xx[,2:(ncol(xx))])),
                                 PRIOR=INIT,
                                 dumi=paste(rep(paste(grep(paste(xx$Fleet,collapse='|'),fleetinfo$fleetname),3,sep='_'),each=ncol(xx)-1),
                                            rep(colnames.xx,nrow(xx)),sep='_'))
      this.par=discard_mortality$INIT
      multiplr=rep(0.1,length(this.par))
      multiplr=ifelse(this.par<0,2,multiplr)
      low.bound=multiplr*discard_mortality$INIT
      id0=which(this.par==0)
      id999=which(this.par==999)
      if(length(id0)>0) low.bound[id0]=-1
      if(length(id999)>0) low.bound[id999]=-10
      low.bound=sapply(low.bound,function(x) max(0,x))
      multiplr=rep(1.5,length(this.par))
      multiplr=ifelse(this.par<0,-2,multiplr)
      up.bound=multiplr*discard_mortality$INIT
      if(length(id0)>0) up.bound[id0]=1
      if(length(id999)>0) up.bound[id999]=10
      discard_mortality=discard_mortality%>%
        mutate(LO=low.bound,
               HI=up.bound,
               PHASE=-2,
               PR_SD=99)%>%
              replace(is.na(.), 0)
      rownames(discard_mortality)=str_replace(rownames(discard_mortality),'SizeSel','DiscMort')
      
      
      #combine with sel pars
      dumiSel=ctl$size_selex_parms
      kkk=strsplit(rownames(dumiSel), '_') 
      dumi.f=data.frame(A=sapply(kkk, `[`, 4),
                        B=sapply(kkk, `[`, 5),
                        Z=sapply(kkk, `[`, 6))%>%
                        mutate(B=ifelse(is.na(B),'',B),
                               Z=ifelse(is.na(Z),'',Z),
                               C=paste(A,B,sep='_'),
                               CZ=paste(A,B,Z,sep='_'),
                               D=ifelse(B=="",A,ifelse(!Z=='',CZ,C)),
                               E=match(D,fleetinfo$fleetname))
      dumiSel$dumi=paste(paste(paste(dumi.f$E,1,sep='_'),sapply(kkk, `[`, 2),sep='_'),sapply(kkk, `[`, 3),sep='_')
      
      ctl$size_selex_parms=rbind(dumiSel,retention_params,discard_mortality)%>%
                              arrange(dumi)%>%
                              dplyr::select(-dumi)
    }
    
    #Male offset 
    if('SS_offset_selectivity'%in%names(life.history) & is.na(Scenario$NSF.selectivity))  
    {
      xx=life.history$SS_offset_selectivity%>%select_if(~ !any(is.na(.)))
      xx.phase=life.history$SS_offset_selectivity_phase%>%select_if(~ !any(is.na(.)))
      id.male.off=grep(paste(xx$Fleet,collapse='|'),rownames(ctl$size_selex_types))
      if(length(id.male.off)>1 & length(unique(xx$Fleet))==1)
      {
        xx=xx %>%
              slice(rep(1:n(), each = length(id.male.off)))%>%
              mutate(Fleet=rownames(ctl$size_selex_types)[id.male.off])
        xx.phase= xx.phase %>%
              slice(rep(1:n(), each = length(id.male.off)))%>%
          mutate(Fleet=rownames(ctl$size_selex_types)[id.male.off])
      }
      ctl$size_selex_types$Male[id.male.off]=3
      ctl$size_selex_types=ctl$size_selex_types%>%
                              mutate(Male=ifelse(Pattern==15,0,Male))
      
      xx=xx%>%
            filter(grepl(paste(rownames(ctl$size_selex_types[ctl$size_selex_types$Male>0,]),collapse='|'),Fleet))
      xx.phase=xx.phase%>%
            filter(grepl(paste(rownames(ctl$size_selex_types[ctl$size_selex_types$Male>0,]),collapse='|'),Fleet))
      if(is.null(size.comp) & is.null(meanbodywt)) xx.phase[,-1]=-1
      
      xx.min=xx%>%mutate(P_1=-50,P_3=-5,P_4=-5,P_5=-8,P_6=0)
      xx.max=xx%>%mutate(P_1=50,P_3=5,P_4=5,P_5=5,P_6=1.5)

      id.offset.patrn=grep(pattern = paste(unique(xx$Fleet),collapse="|"), x = rownames(ctl$size_selex_parms))
      offset_params=ctl$size_selex_parms[id.offset.patrn[grep(pattern = paste(colnames(xx)[-1],collapse="|"),
                                                              x = rownames(ctl$size_selex_parms[id.offset.patrn,]))],]
      rownames(offset_params)=paste(rownames(offset_params),"offset.male",sep='_')
      
      offset_params=offset_params%>%
                      mutate(INIT=c(t(xx[,2:(ncol(xx))])),
                             PRIOR=INIT,
                             PHASE=c(t(xx.phase[,2:(ncol(xx.phase))])),
                             LO=c(t(xx.min[,2:(ncol(xx.min))])),
                             HI=c(t(xx.max[,2:(ncol(xx.max))])))%>%
                      replace(is.na(.), 0)
      offset_params$dumi=sort(paste(paste(match(xx$Fleet,fleetinfo$fleetname),2,sep='_'),rep(names(xx)[-1],nrow(xx)),sep='_'))
      
      #combine with sel pars
      dumiSel=ctl$size_selex_parms
      kkk=strsplit(rownames(dumiSel), '_') 
      dumi.f=data.frame(A=sapply(kkk, `[`, 4),
                        B=sapply(kkk, `[`, 5),
                        Z=sapply(kkk, `[`, 6))%>%
                          mutate(B=ifelse(is.na(B),'',B),
                                 Z=ifelse(is.na(Z),'',Z),
                                 C=paste(A,B,sep='_'),
                                 CZ=paste(A,B,Z,sep='_'),
                                 D=ifelse(B=="",A,ifelse(!Z=='',CZ,C)),
                                 E=match(D,fleetinfo$fleetname))
      dumiSel$dumi=paste(paste(paste(dumi.f$E,1,sep='_'),sapply(kkk, `[`, 2),sep='_'),sapply(kkk, `[`, 3),sep='_')
      ctl$size_selex_parms=rbind(dumiSel,offset_params)%>%
                            arrange(dumi)%>%
                            dplyr::select(-dumi)
    }
    
    #order
    if(life.history$Name=="sawsharks")
    {
      ctl$size_selex_parms$fleet=gsub("^\\.","",str_remove_all(rownames(ctl$size_selex_parms), paste(c("SizeSel_P_", paste0(1:6,"_")), collapse = "|")))
      ctl$size_selex_parms$order=gsub("^\\.","",str_remove_all(rownames(ctl$size_selex_parms), paste(c("SizeSel_P_",paste0("_",ctl$size_selex_parms$fleet)), collapse = "|")))
      ctl$size_selex_parms=ctl$size_selex_parms%>%
        left_join(data.frame(fleet=rownames(ctl$size_selex_types%>%filter(Special==0)))%>%mutate(Fleet.order=row_number()),
                  by='fleet')%>%
        arrange(Fleet.order,order)
      rownames(ctl$size_selex_parms)=paste0('SizeSel_P_',ctl$size_selex_parms$order,'_',ctl$size_selex_parms$fleet)
      ctl$size_selex_parms=ctl$size_selex_parms%>%
        dplyr::select(-order,-fleet,-Fleet.order)
    }
    
    #Add timevary selex parameters values if not estimated
    if('timevary_selex_parameters'%in%names(life.history))
    {
      ctl$size_selex_parms_tv=life.history$timevary_selex_parameters
    }
    
    #add priors if required
    if(!is.null(life.history$Sel.prior.sd_type))
    {
      #P1 Logistic
      if(ctl$size_selex_types[grep("Northern.shark",rownames(ctl$size_selex_types)),'Pattern']==1)
      {
        aidII=grep("SizeSel_P_1_Northern.shark",rownames(ctl$size_selex_parms))
        ctl$size_selex_parms[aidII,'PR_SD']=life.history$Sel.prior.sd_type$P1.sd
        ctl$size_selex_parms[aidII,'PR_type']=life.history$Sel.prior.sd_type$P1.type
        
        aidII=grep("SizeSel_P_2_Northern.shark",rownames(ctl$size_selex_parms))
        ctl$size_selex_parms[aidII,'PR_SD']=life.history$Sel.prior.sd_type$P2.sd
        ctl$size_selex_parms[aidII,'PR_type']=life.history$Sel.prior.sd_type$P2.type
      }
      
      if(ctl$size_selex_types[grep("Survey",rownames(ctl$size_selex_types)),'Pattern']==1)
      {
        aidII=grep("SizeSel_P_1_Survey",rownames(ctl$size_selex_parms))
        ctl$size_selex_parms[aidII,'PR_SD']=life.history$Sel.prior.sd_type$P1.sd
        ctl$size_selex_parms[aidII,'PR_type']=life.history$Sel.prior.sd_type$P1.type
        
        aidII=grep("SizeSel_P_2_Survey",rownames(ctl$size_selex_parms))
        ctl$size_selex_parms[aidII,'PR_SD']=life.history$Sel.prior.sd_type$P2.sd
        ctl$size_selex_parms[aidII,'PR_type']=life.history$Sel.prior.sd_type$P2.type
      }
      
    }
    
    if('Fleet'%in%colnames(ctl$size_selex_parms))ctl$size_selex_parms=ctl$size_selex_parms%>%dplyr::select(-Fleet)
    
    #2. age_selex
    ctl$age_selex_types=ctl$size_selex_types%>%
                          mutate(Pattern=life.history$age_selex_pattern,
                                 Discard=0,
                                 Male=0,
                                 Special=0)
    # ddumy=ctl$age_selex_types%>%
    #   rownames_to_column('fleetname')%>%
    #   mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
    #                           ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
    #                                  fleetname)))%>%
    #   filter(fleetname%in%dis.flits)%>%
    #   mutate(Fleet=row_number())
    # rownames(ddumy)=ddumy$fleetname
    # ddumy$Pattern=life.history$age_selex_pattern
    # idd=rownames(ctl$size_selex_types)[which(!rownames(ctl$size_selex_types)%in%rownames(ddumy))]
    # if(length(idd>0))
    # {
    #   add=ddumy[1:length(idd),]%>%mutate(Fleet=match(idd,rownames(ctl$size_selex_types)))
    #   rownames(add)=idd
    #   ddumy=rbind(ddumy,add)%>%arrange(Fleet)
    # }
    # ctl$age_selex_types=ddumy%>%dplyr::select(-fleetname,-Fleet)
    # if(!is.null(F.tagging))
    # {
    #   F.age.sel.pat=ctl$age_selex_types[match('Southern.shark_2',rownames(ctl$age_selex_types)),]
    #   rownames(F.age.sel.pat)=F.fleet
    #   ctl$age_selex_types=rbind(ctl$age_selex_types,F.age.sel.pat)
    # }  
    
    #Block pattern - time changing selectivity (must estimate the par)
    if((!life.history$Nblock_Patterns==0 & ctl$time_vary_auto_generation[5]==0)|
       ('timevary_selex_parameters'%in%names(life.history)))
    {
      dis1=fleetinfo$fleetname[grep(life.history$Sel.Block.fleet,fleetinfo$fleetname)]
      if(Scenario$Spatial=="areas-as-fleets") dis1=life.history$areas.as.fleet.zone.block
      dis2=life.history$SizeSelex_Block
      id.dis.flit=grep(paste(dis1,collapse='|'),rownames(ctl$size_selex_parms))
      active.pars=paste0('P_',str_extract(row.names(ctl$size_selex_parms[id.dis.flit,]%>%filter(PHASE>0)), "\\d+(?=_+.+$)"))
      dis2[]=0
      id.active.pars=match(life.history$areas.as.fleet.zone.block_pars,names(dis2))
      dis2[id.active.pars]=1
      Block_Fxn=dis2
      Block_Fxn[id.active.pars]=life.history$Sel_param_Blk_Fxn
      dis3=paste(names(dis2),rep(dis1,each=length(dis2)),sep='_')
      dis.blok=match(paste('SizeSel',dis3,sep='_'),rownames(ctl$size_selex_parms))
      
      ctl$size_selex_parms$Block[dis.blok]=rep(dis2,length(dis1))  
      ctl$size_selex_parms$Block_Fxn[dis.blok]=rep(Block_Fxn,length(dis1))
      
      if(ctl$size_selex_parms[id.dis.flit,]$PHASE[id.active.pars]<0 &
         !'timevary_selex_parameters'%in%names(life.history))
      {
        ctl$size_selex_parms[id.dis.flit,]$PHASE[id.active.pars]=4
      }
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
  
  #Tagging   
  if(Scenario$Model=='SS' & Scenario$Tagging=='Yes' & !is.null(Tags))
  {
    ctl$TG_custom=1
    
    dummy.tg.matrx=ctl$size_selex_parms[1,]
    dummy.tg.matrx[,]=NA
    rownames(dummy.tg.matrx)=NULL
    
    seg.TG=seq(1,dat$N_tag_groups)
    
    ctl$TG_Loss_init=dummy.tg.matrx[seg.TG,]%>%
                        mutate(INIT=round(Tags$Initial.tag.loss,3),
                               LO=-15,
                               HI=10,
                               PRIOR=INIT,
                               PR_SD=99,
                               PR_type=6,
                               PHASE=-4,
                               'env_var&link'=0, dev_link=0, dev_minyr=0, dev_maxyr=0, dev_PH=0, Block=0, Block_Fxn=0)
    rownames(ctl$TG_Loss_init)=paste0('TG_Loss_init_',seg.TG)
    
    ctl$TG_Loss_chronic=dummy.tg.matrx[seg.TG,]%>%
                      mutate(INIT=round(Tags$Chronic.tag.loss,3),
                             LO=-15,
                             HI=10,
                             PRIOR=INIT,
                             PR_SD=99,
                             PR_type=6,
                             PHASE=-4,
                             'env_var&link'=0, dev_link=0, dev_minyr=0, dev_maxyr=0, dev_PH=0, Block=0, Block_Fxn=0)
    rownames(ctl$TG_Loss_chronic)=paste0('TG_Loss_chronic_',seg.TG)
    
    ctl$TG_overdispersion=dummy.tg.matrx[seg.TG,]%>%
                      mutate(INIT=round(Tags$overdispersion,3),
                             LO=1,
                             HI=100,
                             PRIOR=INIT,
                             PR_SD=99,
                             PR_type=6,
                             PHASE=-4,
                             'env_var&link'=0, dev_link=0, dev_minyr=0, dev_maxyr=0, dev_PH=0, Block=0, Block_Fxn=0)
    rownames(ctl$TG_overdispersion)=paste0('TG_overdispersion_',seg.TG)
    
    N.flits.tag=nrow(dat$fleetinfo%>%filter(type==1)) 
    ctl$TG_Report_fleet=dummy.tg.matrx[seq(1,N.flits.tag),]%>%
                          mutate(LO=-20,
                                 HI=50,
                                 INIT=-10,   #set to low all fleets other than relevant ones
                                 PRIOR=INIT,
                                 PR_SD=99,
                                 PR_type=6,
                                 PHASE=-4,
                                 'env_var&link'=0, dev_link=0, dev_minyr=0, dev_maxyr=0, dev_PH=0, Block=0, Block_Fxn=0)
    rownames(ctl$TG_Report_fleet)=paste0('TG_report_fleet_par',seq(1,N.flits.tag))
    ctl$TG_Report_fleet$INIT[Tags$Initial.reporting.rate$Fleet]=round(Tags$Initial.reporting.rate$Reporting,3)
    
    ctl$TG_Report_fleet_decay=dummy.tg.matrx[seq(1,N.flits.tag),]%>%
                                  mutate(INIT=round(Tags$Reporting.rate.decay,3),  
                                         LO=-10,
                                         HI=10,
                                         PRIOR=INIT,
                                         PR_SD=99,
                                         PR_type=6,
                                         PHASE=-4,
                                         'env_var&link'=0, dev_link=0, dev_minyr=0, dev_maxyr=0, dev_PH=0, Block=0, Block_Fxn=0)
    rownames(ctl$TG_Report_fleet_decay)=paste0('TG_report_decay_par',seq(1,N.flits.tag))
  }
  
  #Set prior to init value 
  ctl$SR_parms[,"PRIOR"]=ctl$SR_parms[,"INIT"]
  ctl$MG_parms=ctl$MG_parms%>%mutate(PRIOR=ifelse(PHASE<0,INIT,PRIOR))
  ctl$Q_parms[,"PRIOR"]=ctl$Q_parms[,"INIT"]
  ctl$size_selex_parms[,"PRIOR"]=ctl$size_selex_parms[,"INIT"]
  
  #Input variance adjustments factors 
  #Factors: 1=add_to_survey_CV; 2=add_to_discard_stddev; 3=add_to_bodywt_CV; 4=mult_by_lencomp_N
  #         5=mult_by_agecomp_N; 6=mult_by_size-at-age_N; 7=mult_by_generalized_sizecomp
  if(Scenario$Model=='SS')
  {
    if(!is.null(Var.adjust.factor))
    {
      ctl$DoVar_adjust=1
      ctl$Variance_adjustment_list=Var.adjust.factor
    }
  }
  
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
    if(!is.na(Scenario$like_comp.w))
    {
      Like_comp[which( Like_comp$like_comp ==Scenario$like_comp.w & Like_comp$fleet==Scenario$like_comp_fleet.w),'value']=Scenario$like_comp.w.val
    }
    ctl$lambdas=Like_comp
    ctl$N_lambdas=nrow(ctl$lambdas)
  }
  
  
  #3.3. starter file
  start$datfile='data.dat'
  start$ctlfile='control.ctl'
  
  #3.4. forecast file 
  index.F=grep('F.series',dat$fleetinfo$fleetname)
  if(length(index.F)>0)
  {
    add.future.f.series=Future.project[ncol(Future.project)]
    add.future.f.series[,]=0
    colnames(add.future.f.series)=index.F
    Future.project=cbind(Future.project,add.future.f.series)
  }
  
  fore$benchmarks=1 # 0 = skip; 1 = calc F_spr,F_btgt & F_msy
  fore$MSY=2 #1 = set to F(SPR); 2= calc F(MSY); 3 = set to F(Btgt); 4 = set to F(endyr) 
  fore$SPRtarget=0.6
  fore$Btarget=0.4
  fore$Bmark_relF_Basis=1
  fore$Forecast=2 #0=none; 1=F(SPR); 2=F(MSY) 3=F(Btgt); 4=Ave F (uses first-last relF yrs); 5=input annual F scalar
  if(!is.null(Future.project))fore$Nforecastyrs=length(Future.project$finyear)
  fore$Fcast_selex=0 #0=mean from year range, 1=annual time varying
  fore$ControlRuleMethod=1 #0=none, 1=catch as function of SSB, 2=F as function of SSB
  fore$BforconstantF=0.4
  fore$BfornoF=0.1
  fore$Flimitfraction=0.913
  fore$N_forecast_loops=3 #1=OFL only; 2=ABC; 3=get F from forecast ABC catch with allocations applied
  fore$First_forecast_loop_with_stochastic_recruitment=3
  fore$fcast_rec_option=0
  fore$fcast_rec_val=1 #0 = spawner recruit curve,1 = value*(spawner recruit curve), 2 = value*(virgin recruitment), 3 = recent mean from year range above
  fore$Forecast_loop_control_5=0
  if(!is.null(Future.project)) fore$FirstYear_for_caps_and_allocations=Future.project$finyear[nrow(Future.project)]+1
  fore$stddev_of_log_catch_ratio=0
  fore$Do_West_Coast_gfish_rebuilder_output=0
  fore$Ydecl=-1
  fore$Yinit=-1
  fore$fleet_relative_F=1
  
  fore$basis_for_fcast_catch_tuning=2 #2=total catch biomass; 3=retained catch biomass; 5=total catch numbers; 6=retained total numbers
  fore$N_allocation_groups=0
  if(Scenario$Model=='SS')
  {
    if(Scenario$Forecasting=='catch') fore$InputBasis=2 #-1 = Read basis with each observation, 2 = Dead catch (retained + discarded),3 = Retained catch, 99 = Input apical F
    if(Scenario$Forecasting=='F') fore$InputBasis=99
  }

  if(!is.null(Future.project))  # future catch
  {
    if(Scenario$Model=='SSS')
    {
      Future.catch=data.frame(year=Future.project$finyear,
                              seas=1,
                              fleet=1,
                              catch=Future.project$Tonnes)
    }
    if(Scenario$Model=='SS')
    {
      dis.kls=(ncol(Future.project)-nrow(dat$fleetinfo)+1):ncol(Future.project)
      Future.catch=Future.project[,-match(c('SPECIES','Name'),names(Future.project))]%>%
                      gather(fleet,catch,-finyear)%>%
                      mutate(seas=1)%>%
                      arrange(finyear)%>%
                      relocate(seas,.before=fleet)%>%
                      rename(Year=finyear)
    }
    fore$ForeCatch=Future.catch 
  }

  #age composition 
  if(is.null(age.comp))
  {
    dat$N_agebins=0
    if(Scenario$Model=='SSS') dat=dat[-match(c("ageerror"),names(dat))]
    if(Scenario$Model=='SS') dat=dat[-match(c("agebin_vector","N_ageerror_definitions","ageerror"),names(dat))]
  }
  
  
  # 4.Export updated templates
  r4ss::SS_writestarter(start, dir = new.path, overwrite = TRUE,verbose = FALSE)
  r4ss::SS_writedat(dat, outfile = file.path(new.path, start$datfile), overwrite = TRUE, verbose = FALSE)
  r4ss::SS_writectl(ctl, outfile = file.path(new.path, start$ctlfile), overwrite = TRUE, verbose = FALSE)
  r4ss::SS_writeforecast(fore, dir = new.path, file = "forecast.ss", overwrite = TRUE, verbose = FALSE)
}

#Run SS ------------------------------------------------------
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

#MC simulations for posterior probabilities ------------------------------------------------------  
#note: this approach does not work if estimating rec devs 
fn.med.ci.SS=function(d)
{
  dframe=data.frame(year=as.numeric(names(d)),
                    median=apply(d,2,function(x) median(x)),
                    lower.95=apply(d,2,function(x) quantile(x,probs=0.025)),
                    upper.95=apply(d,2,function(x) quantile(x,probs=0.975)))
  
  return(dframe)
}
fn.MC.sims=function(this.wd1,nMC=nMCsims,arg=Arg.no.estimation,B.th,scen,
                    dis.pars=c('SR_parm','selparm','Q_parm'),turn.off.recdevs=TRUE)   
{
  #Get var-covar matrix for relevant parameters
  MLE=read.admbFit(paste(this.wd1,'ss',sep='/'))
  #n.mle=grep(paste(dis.pars,collapse='|'),MLE$names)
  n.mle=1:MLE$nopar
  Nms=MLE$names[n.mle]

  #Monte Carlo sampling
  Depletion=vector('list',nMC)
  F.series=B.Bmsy=F.Fmsy=Depletion
  New.dir=paste(this.wd1,'MonteCarlo',sep='/')
  if(!dir.exists(New.dir))dir.create(New.dir)

  for(n in 1:nMC)   
  {
    unlink(paste(New.dir,'*',sep='/'), recursive = TRUE, force = TRUE)
    copy_SS_inputs(dir.old = this.wd1, dir.new = New.dir,overwrite = TRUE)

    #1. Sample from multivariate distribution
    Rand.par=c(mvtnorm::rmvnorm(1,mean=MLE$est[n.mle],sigma=MLE$cov[n.mle,n.mle]))
    names(Rand.par)=Nms
    
    #2. Create new SS ctl file
    ctl <- r4ss::SS_readctl(file = file.path(this.wd1, "control.ctl")) 
    id=which(ctl$SR_parms$PHASE>0)
    if(length(id)>0)
    {
      dispars=grep('SR_parm',Nms)
      if(length(dispars)>0) ctl$SR_parms$INIT[id]=Rand.par[dispars]
    }
    id=which(ctl$Q_parms$PHASE>0)
    id.Q.float=which(ctl$Q_options$float==0)
    if(length(id)>0 & length(id.Q.float)>0)
    {
      dispars=grep('Q_parm',Nms)[1:length(id)]
      if(length(dispars)>0) ctl$Q_parms$INIT[id]=Rand.par[dispars]
    }
    id=which(ctl$size_selex_parms$PHASE>0)
    if(length(id)>0)
    {
      dispars=grep('selparm',Nms)
      if(length(dispars)>0) ctl$size_selex_parms$INIT[id]=Rand.par[dispars]
    }
    
    #turn off recruitment devs
    if(turn.off.recdevs)
    {
      ctl$recdev_phase=-abs(ctl$recdev_phase)
      ctl$recdev_early_start=0
      ctl$recdev_early_phase=-1
    }

    
    #export ctl
    r4ss::SS_writectl(ctl, outfile = file.path(New.dir, "control.ctl"), overwrite = TRUE, verbose = FALSE)
    
    #3. Run SS without estimation
    fn.run.SS(where.inputs=New.dir,
              where.exe=handl_OneDrive('SS3/ss_win.exe'),
              args=arg)  
    
    #4. Store quantities of interest (depletion)
    Report=SS_output(New.dir,covar=F,forecast=F,readwt=F)
    Years=with(Report,startyr:endyr)
    dum=Report[["derived_quants"]]
    
    Depletion[[n]]=dum[grep(paste(paste("Bratio",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]%>%
      mutate(year=readr::parse_number(Label))%>%
      filter(year%in%Years)%>%
      `rownames<-`( NULL )%>%
      relocate(year)%>%
      dplyr::select(-Label)%>%
      spread(year,Value)
    
    F.series[[n]]=dum[grep(paste(paste("F",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]%>%
      mutate(year=readr::parse_number(Label))%>%
      filter(year%in%Years)%>%
      `rownames<-`( NULL )%>%
      relocate(year)%>%
      dplyr::select(-Label)%>%
      spread(year,Value)
    
    SSB_MSY=dum[dum$Label=="SSB_MSY",'Value']
    B.Bmsy[[n]]=dum[grep(paste(paste("SSB",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]%>%
      mutate(year=readr::parse_number(Label),
             Value=Value/SSB_MSY)%>%
      filter(year%in%Years)%>%
      `rownames<-`( NULL )%>%
      relocate(year)%>%
      dplyr::select(-Label)%>%
      spread(year,Value)
    
    annF_MSY=dum[dum$Label=="annF_MSY",'Value']
    F.Fmsy[[n]]=dum[grep(paste(paste("F",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]%>%
      mutate(year=readr::parse_number(Label),
             Value=Value/annF_MSY)%>%
      filter(year%in%Years)%>%
      `rownames<-`( NULL )%>%
      relocate(year)%>%
      dplyr::select(-Label)%>%
      spread(year,Value)
  }
  Depletion=do.call(rbind,Depletion)
  F.series=do.call(rbind,F.series)
  B.Bmsy=do.call(rbind,B.Bmsy)
  F.Fmsy=do.call(rbind,F.Fmsy)
  
  #Calculate reference point probs
  Probs=add.probs(DAT=Depletion,
                  id.yr=ncol(Depletion),
                  B.threshold=B.th,
                  plot.ranges=FALSE)
  Probs$probs=Probs$probs%>%mutate(Scenario=scen)
  
  #Get quantities of interest  
  Probs$Depletion=fn.med.ci.SS(Depletion)%>%
    mutate(Model='SS3',
           Scenario=scen)
  Probs$F.series=fn.med.ci.SS(F.series)%>%
    mutate(Model='SS3',
           Scenario=scen)
  Probs$B.Bmsy=fn.med.ci.SS(B.Bmsy)%>%
    mutate(Model='SS3',
           Scenario=scen)
  Probs$F.Fmsy=fn.med.ci.SS(F.Fmsy)%>%
    mutate(Model='SS3',
           Scenario=scen)
  
  
  return(Probs)
}

#Compare prior and posterior ------------------------------------------------------
fn.compare.prior.post=function(d,Par,prior_type)
{
  if(prior_type=='beta')
  {
    Post.beta=get_beta(d$Value,d$Parm_StDev)
    Posterior=density(rbeta(1e4,Post.beta[1],Post.beta[2]))
    
    Post.beta=get_beta(d$Prior,d$Pr_SD)
    Prior=density(rbeta(1e4,Post.beta[1],Post.beta[2]))
  }
  if(prior_type=='normal')
  {
    Posterior=density(rnorm(1e4,d$Value,d$Parm_StDev))
    Prior=density(rnorm(1e4,d$Prior,d$Pr_SD))
  }
  
  plot(Posterior$x,Posterior$y,type='l',col=2,xlab=Par,ylab='Density',lwd=2,
       xlim=c(min(c(Posterior$x,Prior$x)),max(c(Posterior$x,Prior$x))),
       ylim=c(min(c(Posterior$y,Prior$y)),max(c(Posterior$y,Prior$y))))
  lines(Prior$x,Prior$y,col=3,lwd=2)
  legend("topleft",c('Posterior','Prior'),col=2:3,lty=1,lwd=2,bty='n')
}

#Fit diagnostic functions ------------------------------------------------------
#note: this analysis follows Carvalho et al 2021
fn.like.range=function(Par.mle,min.par,Par.SE,up,low,ln.out,seq.approach='SE')
{
  if(length(Par.mle)==1)
  {
    if(seq.approach=='min.plus')
    {
      Rango=seq(max(min.par,Par.mle)*(1-up),Par.mle*(1+up),length.out=ln.out-1)
    }
    if(seq.approach=='SE')
    {
      Rango=seq(max(min.par,Par.mle)-1.96*Par.SE,Par.mle+1.96*Par.SE,length.out=ln.out-1)
    }
    return(sort(unique(c(Rango,round(Par.mle,2)))))
  }else
  {
    if(seq.approach=='min.plus') offsets <- seq(-up,low, length.out=ln.out-1)
    if(seq.approach=='SE') offsets <- seq(-Par.SE, Par.SE, length.out=ln.out-1)
    if(0 %in%offsets)
    {
      offsets=subset(offsets,!offsets==0)
      offsets=c(offsets,offsets[length(offsets)]*1.1)
    }
    Rango <- offsets %>%
      purrr::map(~ Par.mle + .x)
    Rango =do.call(rbind,Rango)
    Rango=rbind(Rango,Par.mle)
    row.names(Rango)=NULL
    Rango=Rango[order(Rango[,1], decreasing = TRUE),]

    return(Rango)
  }

}
fn.fit.diag_SS3=function(WD,disfiles,R0.vec,h.vec,M.vec,depl.vec,curSB.vec,Linf.vec.F=NA,Linf.vec.M=NA,
                         exe_path,start.retro=0,end.retro=5,
                         do.like.prof=FALSE,do.retros=FALSE,do.jitter=FALSE,numjitter,
                         outLength.Cross.Val=FALSE,run.in.parallel=TRUE,flush.files=TRUE,
                         COVAR=FALSE,h.input=NULL,drop_LP_CurSB=TRUE,
                         Par_var_profile=c("R0","h","M","Depl","CurSB"))
{
  Report=SS_output(dir=WD,covar=COVAR,verbose=FALSE,printstats=FALSE)
  
  dirname.diagnostics <- paste(WD,"Diagnostics",sep='/')
  if(!dir.exists(dirname.diagnostics)) dir.create(path=dirname.diagnostics, showWarnings = TRUE, recursive = TRUE)
  
  cpue.series=length(unique(Report$cpue$Fleet))
  length.series=length(unique(unique(Report$len_comp_fit_table$Fleet)))
  age.series=length(unique(unique(Report$agedbase$Fleet)))
  
  #1. Goodness-of-fit diagnostic
  #1.1. residuals with smoothing function showing trends
  NrowS=1
  if(cpue.series>0 & length.series>0) NrowS=2
  tiff(file.path(dirname.diagnostics,"jabbaresidual.tiff"),
       width = 1500, height = 2000,units = "px", res = 300, compression = "lzw")
  sspar(mfrow=c(NrowS,1),labs=T,plot.cex=0.9)
  if(cpue.series>0) ss3diags::SSplotJABBAres(ss3rep=Report,subplots = "cpue",add=TRUE,verbose=FALSE)
  if(length.series>0) ss3diags::SSplotJABBAres(ss3rep=Report,subplots = "len",add=TRUE,verbose=FALSE)
  dev.off()
  
  #1.2. runs test
  dis.dat=NULL
  if(cpue.series>0) dis.dat=c(dis.dat,"cpue")
  if(length.series>0) dis.dat=c(dis.dat,"len")
  if(age.series>0) dis.dat=c(dis.dat,"age")
  if(!is.null(dis.dat))
  {
    for(pp in 1:length(dis.dat))
    {
      if(dis.dat[pp]=="cpue")
      {
        nRws=cpue.series
        TAb=table(Report$cpue$Fleet)
        if(any(TAb<=1)) nRws=nRws-1
      }
      if(dis.dat[pp]=="len")
      {
        nRws=length.series
        TAb=table(Report$len_comp_fit_table$Fleet)
        if(any(TAb<=1)) nRws=nRws-1
      }
      if(dis.dat[pp]=="age") nRws=age.series   
      tiff(file.path(dirname.diagnostics,paste0("runs_tests_",dis.dat[pp],".tiff")),
           width = 2000, height = 2000,units = "px", res = 300, compression = "lzw")
      sspar(mfrow=n2mfrow(nRws),labs=T,plot.cex=0.9)
      ss3diags::SSplotRunstest(ss3rep=Report,subplots = dis.dat[[pp]],add=TRUE,verbose =FALSE)
      dev.off()
    }
    runs.test.value=vector('list',length(dis.dat))
    for(pp in 1:length(dis.dat))  runs.test.value[[pp]]=SSrunstest(Report,quants =  dis.dat[[pp]],verbose = FALSE)
    write.csv(do.call(rbind,runs.test.value),paste(dirname.diagnostics,"runs_tests.csv",sep='/'),row.names = FALSE)
  }
  
  
  #2. Model consistency  
    #2.1. Likelihood profile on selected parameters or quantities
  if(do.like.prof) # 30 secs per species-parameter-iteration if in parallel and with hessian; 
  {
    Par_var.vec_profile=list(R0.vec,
                             h.vec,
                             data.frame(M.vec),
                             depl.vec,
                             curSB.vec,
                             Linf.vec.F,Linf.vec.M)
    Par_var_string_profile=list("SR_LN(R0)",
                                "SR_BH_steep",
                                c("NatM_uniform_Fem_GP_1", "NatM_uniform_Mal_GP_1"),
                                "Depl",
                                "CurSB",
                                "L_at_Amax_Fem","L_at_Amax_Mal")
    Par_prof_string=list("R0",'steep',"M","Depl","CurSB","L_at_Amax_Fem","L_at_Amax_Mal")
    Par_prof_label=list(expression(log(italic(R)[0])),'h',expression(italic(M)[a]),"Current depletion",
                        "Current SB","L_at_Amax_Fem","L_at_Amax_Mal")
    
    Par_var_profile.full=c("R0","h","M","Depl","CurSB","L_at_Amax_Fem","L_at_Amax_Mal")
    id.par=match(Par_var_profile,Par_var_profile.full)
    Par_var.vec_profile=Par_var.vec_profile[id.par]
    Par_var_string_profile=Par_var_string_profile[id.par]
    Par_prof_string=Par_prof_string[id.par]
    Par_prof_label=Par_prof_label[id.par]
    names(Par_var.vec_profile)=names(Par_var_string_profile)=names(Par_prof_string)=names(Par_prof_label)=Par_var_profile
    
    if(drop_LP_CurSB)
    {
      drop.this=match("CurSB",Par_var_profile)
      if(!is.na(drop.this))
      {
        Par_var_profile=Par_var_profile[-drop.this]
        Par_var.vec_profile=Par_var.vec_profile[-drop.this]
        Par_var_string_profile=Par_var_string_profile[-drop.this]
        Par_prof_string=Par_prof_string[-drop.this]
        Par_prof_label=Par_prof_label[-drop.this]
      }
    }
    saveoutput=TRUE
    overwrite=TRUE
    use_par_file=TRUE
    for(pp in 1:length(Par_var_profile))
    {
      Par_var=Par_var_profile[pp]
      Par_var.vec=Par_var.vec_profile[[pp]]
      Par_var_string=Par_var_string_profile[[pp]]
      prof_string=Par_prof_string[[pp]]
      prof_label=Par_prof_label[[pp]]
      
      baseval=NULL
      if(Par_var%in%c("R0")) baseval <- round(Report$parameters$Value[grep(Par_var_profile[pp],Report$parameters$Label)],2)
      if(Par_var=="h") baseval <- round(h.input,2) 
      if(Par_var=="M") baseval <- Report$Natural_Mortality[1,match(0,names(Report$Natural_Mortality))]
      if(Par_var=="Depl") baseval<- Report$derived_quants[paste0("Bratio_",Report$endyr),'Value']
      if(Par_var%in%c("L_at_Amax_Fem","L_at_Amax_Mal"))
      {
        baseval <- round(Report$parameters$Value[grep(Par_var_profile[pp],Report$parameters$Label)],2)
      }
        
      
      # Step 1. Identify a directory for the profile likelihood model run(s)
      dirname.base <- paste(dirname.diagnostics,paste0("Profile_",Par_var),sep='/')
      
      # Step 2. Identify a directory where the completed base model run is located
      dirname.completed.model.run<-WD
      
      # Step 3. Create a "x_profile" subdirectory and set as the working directory
      dirname.Par_var.profile<- dirname.base
      if(!dir.exists(dirname.Par_var.profile)) dir.create(path=dirname.Par_var.profile, showWarnings = TRUE, recursive = TRUE)
      mydir <- dirname.Par_var.profile
      setwd(dirname.Par_var.profile)
      
      # Step 4. Create a "Figures_Tables" subdirectory
      plotdir=paste0(dirname.Par_var.profile, "/Figures & Tables")
      if(!dir.exists(plotdir)) dir.create(path=plotdir, showWarnings = TRUE, recursive = TRUE)
      
      
      # Step 5. Create a "Reference_run" subdirectory and copy completed base model output to this directory
      #reference.dir <- paste0(mydir,'/Reference_run') 
      #if(!dir.exists(reference.dir)) dir.create(path=reference.dir, showWarnings = TRUE, recursive = TRUE)
      #file.copy(from=Sys.glob(paste(dirname.completed.model.run, "*.*", sep="/"), dirmark = FALSE),
      #          to=reference.dir)
      #for(nn in disfiles){file.copy(paste(dirname.completed.model.run,"/", nn, sep='')  ,     to=reference.dir)}
      
      
      # Step 6. Copy necessary files from the "Reference_run" subdirectory to the "x_profile" working directory 
      copylst<-disfiles[-match("Report.sso",disfiles)]
      for(nn in copylst){file.copy(paste(WD,"/", nn, sep=''),file.path(dirname.Par_var.profile))}
      
      # Step 7. Edit "control.ss" in the working directory to estimate at least one parameter in each phase
      # E.g., 
      control.file <- readLines(paste(dirname.Par_var.profile, "/control.ss_new", sep=""))
      linen <- NULL
      linen <- grep("#_recdev phase", control.file)
      control.file[linen] <- paste0("1 #_recdev phase")
      write(control.file, paste(dirname.Par_var.profile, "/control.ss_new", sep=""))
      
      # Step 8. Edit "starter.ss" in the working directory to read from init values from control_modified.ss
      starter.file <- readLines(paste(dirname.Par_var.profile, "/starter.ss", sep=""))
      linen <- NULL
      linen <- grep(paste(c("# 0=use init values in control file; 1=use ss.par","#_init_values_src"),collapse='|'), starter.file)
      starter.file[linen] <- paste0("0 # 0=use init values in control file; 1=use ss.par")
      #if(like.prof.case=='faster') starter.file[grep("#_converge_criterion", starter.file)] <- paste0("0.001 #_converge_criterion")
      write(starter.file, paste(dirname.Par_var.profile, "/starter.ss", sep=""))
      
      
      # Step 9. Begin Likelihood profile_x_example.R
      
      ###Working directory
      setwd(dirname.Par_var.profile)
      
      
      # Step 10. Parameter profile
      
      # vector of values to profile over
      if(is.vector(Par_var.vec))
      {
        Nprof <- length(Par_var.vec)
        profilevec=Par_var.vec
      }else
      {
        Nprof <- nrow(Par_var.vec)
        profilevec=Par_var.vec[,1]
      }
      names(profilevec)=paste0('replist',1:Nprof)
      
      #Define the starter file
      starter <- SS_readstarter(file.path(mydir, "starter.ss"), verbose = FALSE)
      
      #Change control file name in the starter file
      starter$ctlfile <- "control_modified.ss" 
      
      # Make sure the prior likelihood is calculated for non-estimated quantities
      starter$prior_like <- 1                           
      SS_writestarter(starter, dir=mydir, overwrite=TRUE,verbose = FALSE)
      
      #Run SS_profile command
      if(run.in.parallel)
      {
        ncores <- parallelly::availableCores(omit = 1)
        future::plan(future::multisession, workers = ncores)
      }
      if(Par_var%in%c("R0","h"))
      {
        Like.profile <- r4ss::profile(dir=mydir, 
                                      oldctlfile="control.ss_new",
                                      newctlfile="control_modified.ss",
                                      string=Par_var_string,
                                      profilevec=Par_var.vec,
                                      exe=exe_path,
                                      saveoutput = saveoutput,
                                      overwrite = overwrite,
                                      verbose = FALSE,
                                      extras = diag.extras)
      }
      if(Par_var%in%c("M","Depl","CurSB","L_at_Amax_Fem","L_at_Amax_Mal"))
      {
        whichruns=1:Nprof
        
        #Create all iteration folder
        for(n in whichruns)
        {
          run_dir=paste(mydir,paste0('profile',n),sep='/')
          dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
          FILEs=list.files(run_dir)
          if("ss.par"%in%FILEs)
          {
            for(f in 1:length(FILEs)) unlink(paste(run_dir, FILEs[f], sep="/"), recursive = TRUE, force = TRUE) 
          }
          for(nn in copylst){file.copy(paste(file.path(dirname.Par_var.profile),"/", nn, sep=''),run_dir)}
          if(Par_var%in%c("Depl","CurSB")) file.copy(paste(file.path(WD),"/", "ss.par", sep=''),run_dir)
          
          
          #Modify value of profiled quantity
          dat_temp <- SS_readdat(file.path(run_dir, "data.dat"), verbose = FALSE)
          ctl_temp <- SS_readctl(file = file.path(run_dir, "control.ss_new"),datlist = dat_temp,verbose = FALSE)
          
          if(Par_var=="M")  ctl_temp$natM["natM1", ]=ctl_temp$natM["natM2", ] = Par_var.vec[n,]
          
          if(Par_var%in%c("L_at_Amax_Fem","L_at_Amax_Mal"))
          {
            id.L_at_amax=grep(Par_var,rownames(ctl_temp$MG_parms))
            ctl_temp$MG_parms[id.L_at_amax,"INIT"]= Par_var.vec[n]
            ctl_temp$MG_parms[id.L_at_amax,'LO']=min(ctl_temp$MG_parms[id.L_at_amax,'LO'],0.9*Par_var.vec[n])
            ctl_temp$MG_parms[id.L_at_amax,'HI']=max(ctl_temp$MG_parms[id.L_at_amax,'HI'],1.1*Par_var.vec[n])
            ctl_temp$MG_parms[id.L_at_amax,"PHASE"]=-abs(ctl_temp$MG_parms[id.L_at_amax,"PHASE"])
          }
            
          
          if(Par_var%in%c("Depl","CurSB"))
          {
            # --- Modify data file ---
            dat_temp$Nfleets <- dat_temp$Nfleets + 1
            new_fleet_num <- dat_temp$Nfleets
            
            # Add fleet info (duplicating a survey is a safe way)
            dat_temp$fleetinfo <- rbind(dat_temp$fleetinfo, dat_temp$fleetinfo[1, ])
            dat_temp$fleetinfo[new_fleet_num, "fleetname"] <- "Depletion_Survey"
            dat_temp$fleetinfo[new_fleet_num, "type"] <- 3
            
            # Make sure its fishery_timing is 1
            dat_temp$fleetinfo[new_fleet_num, "surveytiming"] <- 1
            
            # determine which unit for indices, 34 = depletion, 30= CurSB
            indices_units <- ifelse(Par_var == "Depl",34,30)
            
            # Add CPUE info
            dat_temp$CPUEinfo <- rbind(as.data.frame(dat_temp$CPUEinfo), 
                                       c(new_fleet_num, indices_units, 0, 0)) 
            
            # Add settings row for the new fleet to lencomp and agecomp info
            new_comp_info_row <- data.frame(
              mintailcomp = -1, addtocomp = 0.001, combine_M_F = 0,
              CompressBins = 0, CompError = 0, ParmSelect = 0, minsamplesize = 0.001)
            row.names(new_comp_info_row) <- "Depletion_Survey"
            dat_temp$len_info <- rbind(dat_temp$len_info, new_comp_info_row)
            dat_temp$age_info <- rbind(dat_temp$age_info, new_comp_info_row)
            
            # Add the actual index data lines, using the value from the profile vector `vec` 
            if (Par_var %in% c("Depl"))
            {
              new_indices <- data.frame(
                year = c(dat_temp$styr - 1, dat_temp$endyr), month = 1,
                index = new_fleet_num, obs = c(1.0, Par_var.vec[n]), se_log = 0.0001)
            }
            if (Par_var %in% c("CurSB"))
            {
              new_indices <- data.frame(
                year = dat_temp$endyr, month = 1,
                index = new_fleet_num, obs = Par_var.vec[n], se_log = 0.0001)
            }
            dat_temp$CPUE <- rbind(dat_temp$CPUE, new_indices)
            SS_writedat(dat_temp, file.path(run_dir, "data.dat"), overwrite = TRUE, verbose = FALSE)
            
            # --- Modify control file ---
            # ctl_temp$Q_options <- rbind(ctl_temp$Q_options, c(new_fleet_num, 1, 0, 0, 0, 0))
            
            # Create the vector with the correct column names
            new_q_row <- c(new_fleet_num, 1, 0, 0, 0, 0)
            names(new_q_row) <- names(ctl_temp$Q_options)
            
            # Add the row to the data frame and assign the row name in one step
            ctl_temp$Q_options <- rbind(ctl_temp$Q_options, Depletion_Survey = new_q_row)
            new_q_parm <- c(-15, 15, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0) # Phase -1 makes it non-estimated
            ctl_temp$Q_parms <- rbind(ctl_temp$Q_parms, Depletion_Survey = new_q_parm)
            ctl_temp$size_selex_types <- rbind(ctl_temp$size_selex_types, Depletion_Survey = c(0, 0, 0, 0)) # Non-selective
            ctl_temp$age_selex_types <- rbind(ctl_temp$age_selex_types, Depletion_Survey = c(0, 0, 0, 0))   # Non-selective
            
            # --- Modify par file --- 
            starter <- SS_readstarter(file.path(run_dir, "starter.ss"), verbose = FALSE)
            if(use_par_file)
            {
              starter[["init_values_src"]] <- 1 # use par file as initial values instead of ctl file
              SS_writestarter(starter, dir = run_dir, overwrite = TRUE, verbose = FALSE)
              
              # read the par file
              # 2. Read all lines from the original file
              lines <- readLines(file.path(run_dir, "ss.par"))
              
              # 3. Find the indices of all lines containing a Q_parm entry
              q_parm_indices <- grep("# Q_parm\\[", lines)
              
              # Check if any Q_parm lines were found
              if (length(q_parm_indices) > 0)
              {
                # 4. Get the index and content of the *last* Q_parm line
                last_q_parm_label_index <- tail(q_parm_indices, 1)
                last_q_parm_label_line <- lines[last_q_parm_label_index]
                
                # 5. Extract the number from the last Q_parm line and increment it
                # This regular expression extracts the digits from inside "[ ]"
                last_q_number <- as.numeric(gsub(".*Q_parm\\[(\\d+)\\].*", "\\1", last_q_parm_label_line))
                new_q_number <- last_q_number + 1
                
                # 6. Create the new lines to be added
                new_content <- c(paste0("# Q_parm[", new_q_number, "]:"),"0")
                
                # 7. Insert the new content right after the value of the last Q_parm
                # The insertion point is after the last Q_parm's value line
                insertion_point <- last_q_parm_label_index + 1
                modified_lines <- append(lines, new_content, after = insertion_point)
                
                # 8. Write the modified lines to a new file
                writeLines(modified_lines, file.path(run_dir, "ss.par"))
                
                #cat("Successfully added '# Q_parm[", new_q_number, "]' to the file '", file.path(run_dir, "ss3.par"), "'.\n", sep = "")
                
              }else
              {
                cat("No '# Q_parm' entries were found in the file. No changes were made.\n")
              }
            }
          }
          
          SS_writectl(ctl_temp, file.path(run_dir, "control_modified.ss"), overwrite = TRUE, verbose = FALSE)
          
        }
        
        #Run estimation in parallel
        res <- furrr::future_map(whichruns, function(i)
        {
          run_dir=paste(mydir,paste0('profile',i),sep='/')
          fn.run.SS(where.inputs=run_dir,where.exe=handl_OneDrive('SS3/ss_win.exe'),args=diag.extras)
          #run(dir = profile_dir, verbose = FALSE, exe = handl_OneDrive('SS3/ss_win.exe'),...)
          repfile_loc <- file.path(run_dir, "Report.sso")
          if (file.exists(repfile_loc) & file.info(repfile_loc)[["size"]] > 0)
          {
            goodrep <- TRUE
            Rep <- readLines(repfile_loc, n = 400)
            convergence_line <- grep("Convergence_Level",Rep)
            max_grad <- as.numeric(stringr::str_extract(Rep[convergence_line], 
                                                        "[[:digit:]|\\.|\\-|e]{2,}"))
            converged <- max_grad <= 1e-4
            skip <- grep("LIKELIHOOD", Rep)[2]
            nrows <- grep("Crash_Pen", Rep) - skip - 1
            like <- read.table(repfile_loc, skip = skip, 
                               nrows = nrows, header = TRUE, fill = TRUE)
            likevec <- as.numeric(like[["logL.Lambda"]])
            names(likevec) <- like[["Component"]]
          } else
          {
            goodrep <- FALSE
            converged <- FALSE
            max_grad <- NA
            likevec <- rep(NA, 10)
          }
          return(list(goodrep = goodrep, converged = converged, max_grad = max_grad, likevec = likevec))
          
        })
        
        if (saveoutput) {
          purrr::walk(whichruns, function(i) {
            profile_dir <- file.path(mydir, paste0("profile", i))
            if (file.exists(file.path(profile_dir, "Report.sso")) & 
                file.info(file.path(profile_dir, "Report.sso"))[["size"]] > 0)
            {
              file.copy(file.path(profile_dir, "Report.sso"), 
                        file.path(mydir, paste0("Report", i, ".sso")), 
                        overwrite = overwrite)
              file.copy(file.path(profile_dir, "CompReport.sso"), 
                        file.path(mydir, paste0("CompReport", i, ".sso")), 
                        overwrite = overwrite)
              file.copy(file.path(profile_dir, "covar.sso"), 
                        file.path(mydir, paste0("covar", i, ".sso")), 
                        overwrite = overwrite)
              file.copy(file.path(profile_dir, "warning.sso"), 
                        file.path(mydir, paste0("warning", i, ".sso")), 
                        overwrite = overwrite)
              file.copy(file.path(profile_dir, "admodel.hes"), 
                        file.path(mydir, paste0("admodel", i, ".hes")), 
                        overwrite = overwrite)
              # file.copy(file.path(profile_dir, parfile), 
              #            file.path(mydir,paste0(parfile, "_", i, ".sso")), overwrite = overwrite)
            }
          })
        }
        if(flush.files) purrr::walk(whichruns, ~unlink(file.path(mydir, paste0("profile",.x)), recursive = TRUE))
        
        res_keep <- which(!sapply(res, is.null))
        res_clean <- res[res_keep]
        goodrep <- sapply(res_clean, function(x) x[["goodrep"]])
        if (!any(goodrep)) 
          stop("Error: no good Report.sso files created in profile")
        liketable <- as.data.frame(t(sapply(res_clean, function(x) x[["likevec"]])))
        Like.profile <- cbind(Value = profilevec[whichruns[res_keep]], 
                              converged = sapply(res_clean, function(x) x[["converged"]]), 
                              liketable, max_grad = sapply(res_clean, function(x) x[["max_grad"]]))
      }
      
      future::plan(future::sequential)
      
      # export like prof stats
      write.csv(Like.profile,
                file.path(dirname.diagnostics,paste0("profile_summary_table_",Par_var,".csv")),
                row.names = TRUE)
      
      # read the output files (with names like Report1.sso, Report2.sso, etc.)
      prof.models <- SSgetoutput(dirvec=mydir, keyvec=c('',1:Nprof), getcovar = FALSE, verbose=FALSE) 
      if(any(is.na(prof.models)))prof.models=prof.models[-which(is.na(prof.models))]  #remove empty list elements (no convergence)
      
      
      # Step 11.  summarize output
      prof.summary <- SSsummarize(prof.models,verbose = FALSE)
      
      # Likelihood components 
      mainlike_components <- c('TOTAL',"Survey", "Catch", 'Length_comp',
                               "Age_comp","Mean_body_wt",'Recruitment') 
      
      mainlike_components_labels  <- c('Total likelihood','Index likelihood',"Catch",'Length likelihood',
                                       "Age likelihood","Mean body weight",'Recruitment likelihood') 
      
      
      # Plot profile using summary created above
      PR_st=prof_string
      PR_la=prof_label
      VEC.prof=NULL
      if(Par_var%in%c("M","Depl","CurSB"))
      {
        PR_st="R0"
        VEC.prof=profilevec
      }
      
      tiff(file.path(dirname.diagnostics,paste0("profile_plot_",Par_var,".tiff")),
           width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
      par(mar=c(5,4,1,1))
      SSplotProfile1(summaryoutput=prof.summary,
                     profile.string = PR_st,
                     profile.label=PR_la,
                     Xvec=VEC.prof,
                     minfraction = 0.001,
                     pheight=4.5,
                     print=FALSE,
                     plotdir=plotdir,
                     components = mainlike_components,
                     component.labels = mainlike_components_labels,
                     add_cutoff = TRUE,
                     cutoff_prob = 0.95,
                     verbose = FALSE)
      if(!is.null(baseval)) abline(v = baseval, lty=2,col='orange',lwd=2)
      legend('right','Base value',lty = 2,col='orange',lwd=2,bty='n')
      dev.off()
      
      # make timeseries plots comparing models in profile
      labs <- paste0(Par_var_string,"= ",Par_var.vec)
      if(!is.null(baseval))labs[which(round(Par_var.vec,2)==baseval)] <- paste0(Par_var_string,"= ",baseval,"(Base model)")
      SSplotComparisons(prof.summary,legendlabels=labs,pheight=4.5,plot = FALSE,png=TRUE,
                        plotdir=plotdir,legendloc='bottomleft',verbose = FALSE)
      
      if(Par_var%in%c("R0"))
      {
        #Piner plot
        #Size comp
        if(length.series>0)
        {
          tiff(file.path(dirname.diagnostics,paste0("profile_plot_",Par_var,"_Length_like.tiff")),
               width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
          par(mar=c(5,4,1,1))
          PinerPlot(prof.summary, 
                    profile.string = PR_st, 
                    component = "Length_like",
                    main = "Changes in length-composition likelihoods by fleet",
                    add_cutoff = TRUE,
                    cutoff_prob = 0.95,
                    verbose = FALSE)
          if(!is.null(baseval))abline(v = baseval, lty=2,col='orange',lwd=2)
          legend('right','Base value',lty = 2,col='orange',lwd=2,bty='n')
          dev.off()
        }
        #Survey
        if(cpue.series>0)
        {
          tiff(file.path(dirname.diagnostics,paste0("profile_plot_",Par_var,"_Survey_like.tiff")),
               width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
          par(mar=c(5,4,1,1))
          PinerPlot(prof.summary,
                    profile.string = PR_st,
                    component = "Surv_like",
                    main = "Changes in Index likelihoods by fleet",
                    add_cutoff = TRUE,
                    cutoff_prob = 0.95,
                    legendloc="topleft",
                    verbose = FALSE)
          if(!is.null(baseval))abline(v = baseval, lty=2,col='orange',lwd=2)
          legend('right','Base value',lty = 2,col='orange',lwd=2,bty='n')
          dev.off()
        }
      }
      
      
      # Step 12.  Remove prof likelihood files
      if(flush.files) 
      {
        dropfiles=list.files(dirname.Par_var.profile)
        dropfiles=subset(dropfiles,!dropfiles=='Figures & Tables')
        for(f in 1:length(dropfiles)) unlink(paste(dirname.Par_var.profile, dropfiles[f], sep="/"), recursive = TRUE, force = TRUE) 
      }
      
    }

    
    #display M.vec
    if("M"%in%Par_var_profile)  
    {
      pi=M.vec%>%
        data.frame()%>%
        mutate(Profile=as.character(1:nrow(M.vec)))%>%
        gather(Age,M,-Profile)%>%
        mutate(Age=as.numeric(substr(Age,5,7)))
      pi=rbind(pi,
               Report$Natural_Mortality[1,match(unique(pi$Age),names(Report$Natural_Mortality))]%>%
                 mutate(Profile='Base value')%>%
                 gather(Age,M,-Profile)%>%
                 mutate(Age=as.numeric(Age)))
      pi%>%
        ggplot(aes(Age,M,color=Profile))+
        geom_line()+
        geom_point(size=2.5)+
        theme_PA()+theme(legend.position = 'top')
      ggsave(file.path(dirname.diagnostics,"profile_M_vec.tiff"),width = 8,height = 8,compression = "lzw")
    }
    
  }
  
    #2.2. Retrospective analysis and Predicting skills
  if(do.retros)      
  {
    dirname.Retrospective <- paste(dirname.diagnostics,"Retrospective",sep='/')
    if(!dir.exists(dirname.Retrospective)) dir.create(path=dirname.Retrospective, showWarnings = TRUE, recursive = TRUE)
    plots.Retrospective <- paste(dirname.Retrospective,"Plots",sep='/')
    if(!dir.exists(plots.Retrospective)) dir.create(path=plots.Retrospective, showWarnings = TRUE, recursive = TRUE)
    file.copy(Sys.glob(paste(WD, "*.*", sep="/"), dirmark = FALSE),dirname.Retrospective)
    if(run.in.parallel)
    {
      ncores <- parallelly::availableCores(omit = 1)
      future::plan(future::multisession, workers = ncores)
    }
    retro(dir = dirname.Retrospective,
          years = start.retro:-end.retro,
          exe=exe_path,
          verbose = FALSE,
          extras= diag.extras)
    future::plan(future::sequential)
    retroModels <- SSgetoutput(dirvec = file.path(dirname.Retrospective, "retrospectives", paste("retro", start.retro:-end.retro, sep = "")))
    if(any(is.na(retroModels)))retroModels=retroModels[-which(is.na(retroModels))]
    if(Neim=='milk shark') retroModels=retroModels[-1]   
    retroSummary <- r4ss::SSsummarize(retroModels,verbose = FALSE)
    if("len"%in%dis.dat & outLength.Cross.Val) retroComp= SSretroComps(retroModels)
    
    # make timeseries plots comparing models in profile
    endyrvec <- retroSummary[["endyrs"]] + start.retro:-end.retro
    SSplotComparisons(summaryoutput=retroSummary,
                      subplots= c(2,4,6,12,14),
                      endyrvec = endyrvec,
                      plot = FALSE,
                      png=TRUE,
                      plotdir=plots.Retrospective,
                      legendlabels = paste("Data", start.retro:-end.retro, "years"),
                      verbose = FALSE)
    tiff(file.path(dirname.diagnostics,"retro_Mohns_Rho.tiff"),
         width = 2000, height = 1800,units = "px", res = 300, compression = "lzw")
    sspar(mfrow=c(2,2),plot.cex=0.8)
    #full series
    rb.full = SSplotRetro(retroSummary,add=T,forecast = F,legend = F,verbose=F)
    rf.full = SSplotRetro(retroSummary,add=T,subplots="F", forecast = F,legendloc="topleft",legendcex = 0.8,verbose=F)
    #last 1- years
    rb = SSplotRetro(retroSummary,add=T,forecast = T,legend = F,verbose=F,xmin=min(endyrvec)-10)
    rf = SSplotRetro(retroSummary,add=T,subplots="F", forecast = T,legendloc="topleft",legendcex = 0.8,verbose=F,xmin=min(endyrvec)-10)
    dev.off()
    
    #Get MASE as metric of prediction skill  
    if(cpue.series>0)
    {
      hcI = SSmase(retroSummary,verbose = FALSE)
      out.MASE=hcI
      write.csv(out.MASE,paste(dirname.diagnostics,"retro_hcxval_MASE.csv",sep='/'),row.names = FALSE)
      write.csv(SShcbias(retroSummary,verbose=FALSE),paste(dirname.diagnostics,"retro_Mohns_Rho.csv",sep='/'),row.names = FALSE)
      if(exists('retroComp'))
      {
        hcL = SSmase(retroSummary=retroComp,quants = "len",verbose = FALSE)
        out.MASE=rbind(out.MASE,hcL)
        write.csv(out.MASE,paste(dirname.diagnostics,"retro_hcxval_MASE.csv",sep='/'),row.names = FALSE)
        write.csv(SShcbias(retroSummary,verbose = FALSE),paste(dirname.diagnostics,"retro_Mohns_Rho.csv",sep='/'),row.names = FALSE)
      }
      
      #Hindcasting cross validation  
      tiff(file.path(dirname.diagnostics,"Hindcasting cross-validation_Survey.tiff"),
           width = 2000, height = 1800,units = "px", res = 300, compression = "lzw")
      nRws=length(unique(Report$cpue%>%filter(Yr%in%endyrvec)%>%pull(Fleet)))
      sspar(mfrow=n2mfrow(nRws),plot.cex=0.8)
      hci = SSplotHCxval(retroSummary,add=T,verbose=F,ylimAdj = 1.2,legendcex = 0.7)
      dev.off()
      
      if(exists('retroComp'))
      {
        nRws=length(unique(Report$len_comp_fit_table%>%filter(Yr%in%endyrvec)%>%pull(Fleet)))
        if(nRws>0)
        {
          tiff(file.path(dirname.diagnostics,"Hindcasting cross-validation_Length.tiff"),
               width = 2000, height = 1800,units = "px", res = 300, compression = "lzw")
          sspar(mfrow=n2mfrow(nRws),plot.cex=0.8)
          hci = SSplotHCxval(retroComp,subplots="len",add=T,verbose=F,ylimAdj = 1.2,legendcex = 0.7)
          dev.off()
        }
      }
    }
    
    #remove retro files
    if(flush.files)
    {
      dropfiles=list.files(dirname.Retrospective)
      dropfiles=subset(dropfiles,!dropfiles%in%c('Plots'))
      if(length(dropfiles)>0)for(f in 1:length(dropfiles)) unlink(paste(dirname.Retrospective, dropfiles[f], sep="/"), recursive = TRUE, force = TRUE) 
    }
   }
  
  
  #3. Convergence
    #3.1. Jittering 
  if(do.jitter)
  {
    #Create directories
    dirname.Jitter <- paste(dirname.diagnostics,"Jitter",sep='/')
    if(!dir.exists(dirname.Jitter)) dir.create(path=dirname.Jitter, showWarnings = TRUE, recursive = TRUE)
    file.copy(Sys.glob(paste(WD, "*.*", sep="/"), dirmark = FALSE),dirname.Jitter)
    
    #Run jitter in parallel
    if(run.in.parallel)
    {
      ncores <- parallelly::availableCores(omit = 1)
      future::plan(future::multisession, workers = ncores)
    }
    jit.likes <- r4ss::jitter(dir = dirname.Jitter,
                              Njitter = numjitter,
                              jitter_fraction =0.1 ,  #0.05
                              init_values_src = 0,
                              exe=exe_path,
                              verbose = FALSE,
                              extras=diag.extras) # 5 times faster with "-nohess" but no uncertainty estim
    future::plan(future::sequential)
    
    #Read in results using other r4ss functions
    keyvec_1=0 #0 is basecase
    profilemodels <- SSgetoutput(dirvec = dirname.Jitter, keyvec = keyvec_1:numjitter, getcovar = FALSE,verbose=FALSE) 
    if(any(is.na(profilemodels)))profilemodels=profilemodels[-which(is.na(profilemodels))]
    profilesummary <- SSsummarize(profilemodels,verbose = FALSE)
    Total.likelihoods=profilesummary[["likelihoods"]][1, -match('Label',names(profilesummary[["likelihoods"]]))]
    Params=profilesummary[["pars"]]
    Base.model.like=Report[["likelihoods_used"]][match("TOTAL",rownames(Report[["likelihoods_used"]])),'values']
    
    #Plot
    if(!'replist0'%in%names(profilesummary[["SpawnBioLower"]])) profilesummary[["SpawnBioLower"]]$replist0=profilesummary[["SpawnBioLower"]]$replist1
    if(!'replist0'%in%names(profilesummary[["SpawnBioUpper"]])) profilesummary[["SpawnBioUpper"]]$replist0=profilesummary[["SpawnBioUpper"]]$replist1

      #Biomass
    YLIM=c(0,max(profilesummary[["SpawnBioUpper"]]$replist0))
    dummy.dat=profilesummary[["SpawnBio"]][,grep('replist',names(profilesummary[["SpawnBio"]]))]%>%
                mutate(Yr=profilesummary[["SpawnBio"]]$Yr)%>%
                gather(replist,value,-Yr)
    p1=profilesummary[["SpawnBioLower"]]%>%
      ggplot(aes(Yr,replist0))+
      geom_line(color='transparent')+
      ylim(YLIM)+ylab('SSB (t)')+xlab('Year')+
      theme_PA()+
      geom_polygon(data=data.frame(Yr=c(profilesummary[["SpawnBioLower"]]$Yr,rev(profilesummary[["SpawnBioLower"]]$Yr)),
                                   replist0=c(profilesummary[["SpawnBioLower"]]$replist0,rev(profilesummary[["SpawnBioUpper"]]$replist0))),
                   colour='grey60',fill='grey60',alpha=.6)+
      geom_line(data=dummy.dat,aes(Yr,value,color=replist),linewidth=0.9)+
      theme(legend.position = 'none')
    
      #Jitters
    YLIM=c(min(Total.likelihoods)*.9,min(Total.likelihoods)*1.1)
    First.jit=-1
    dumydat=data.frame(Jitter.run=1:length(Total.likelihoods[First.jit]),
                       Tot.like=unlist(Total.likelihoods[First.jit]))%>%
      mutate( Col=ifelse(Tot.like==Base.model.like,'Base model',
                  ifelse(Tot.like<Base.model.like,'Smaller',
                  'Higher')))
    p2=dumydat%>%
      ggplot(aes(Jitter.run,Tot.like,color=Col))+
      geom_hline(yintercept=Base.model.like, linetype='longdash')+
      geom_point(size=4)+
      theme_PA()+
      theme(legend.position = 'top')+
      ylim(YLIM)+
      labs(colour = NULL)+
      xlab('Jitter runs at a converged solution')+ylab('Total likelihood')+
      geom_text_repel(aes(label=round(Tot.like,2)),box.padding=1)
    
    ggarrange(plotlist = list(p1,p2),ncol=1)
    ggsave(file.path(dirname.diagnostics,"Jitter.tiff"),width = 6,height = 7,compression = "lzw")
    
    
    #remove jitter files
    if(flush.files) unlink(paste(WD, "Diagnostics/Jitter", sep="/"), recursive = TRUE, force = TRUE)
      
  }
  
}
SSplotProfile1=function (summaryoutput, plot = TRUE, print = FALSE, models = "all", 
                         profile.string = "steep", profile.label = NULL, Xvec=NULL, exact = FALSE, 
                         ylab = "Change in -log-likelihood", components = c("TOTAL", 
                                                                            "Catch", "Equil_catch", "Survey", "Discard", "Mean_body_wt", 
                                                                            "Length_comp", "Age_comp", "Size_at_age", "SizeFreq", 
                                                                            "Morphcomp", "Tag_comp", "Tag_negbin", "Recruitment", 
                                                                            "InitEQ_Regime", "Forecast_Recruitment", "Parm_priors", 
                                                                            "Parm_softbounds", "Parm_devs", "F_Ballpark", "Crash_Pen"), 
                         component.labels = c("Total", "Catch", "Equilibrium catch", 
                                              "Index data", "Discard", "Mean body weight", "Length data", 
                                              "Age data", "Size-at-age data", "Generalized size data", 
                                              "Morph composition data", "Tag recapture distribution", 
                                              "Tag recapture total", "Recruitment", "Initital equilibrium recruitment", 
                                              "Forecast recruitment", "Priors", "Soft bounds", "Parameter deviations", 
                                              "F Ballpark", "Crash penalty"), minfraction = 0.01, sort.by.max.change = TRUE, 
                         col = "default", pch = "default", lty = 1, lty.total = 1, 
                         lwd = 2, lwd.total = 3, cex = 1, cex.total = 1.5, xlim = "default", 
                         ymax = "default", xaxs = "r", yaxs = "r", type = "o", legend = TRUE, 
                         legendloc = "topright", pwidth = 6.5, pheight = 5, punits = "in", 
                         res = 300, ptsize = 10, cex.main = 1, plotdir = NULL, add_cutoff = FALSE, 
                         cutoff_prob = 0.95, add_no_prior_line = TRUE, verbose = TRUE, 
                         ...) 
{
  if (print) {
    if (is.null(plotdir)) {
      stop("to print PNG files, you must supply a directory as 'plotdir'")
    }
    if (!file.exists(plotdir)) {
      if (verbose) {
        message("creating directory:", plotdir)
      }
      dir.create(plotdir, recursive = TRUE)
    }
  }
  if (length(components) != length(component.labels))
  {
    stop("Inputs 'components' and 'component.labels' should have equal length")
  }
  n <- summaryoutput[["n"]]
  likelihoods <- summaryoutput[["likelihoods"]]
  if (is.null(likelihoods)) {
    stop("Input 'summaryoutput' needs to be a list output from SSsummarize\n", 
         "and have an element named 'likelihoods'.")
  }
  pars <- summaryoutput[["pars"]]
  par_prior_likes <- summaryoutput[["par_prior_likes"]]
  if (models[1] == "all")
  {
    models <- 1:n
  }else
  {
    if (!all(models %in% 1:n)) {
      stop("Input 'models' should be a vector of values from 1 to n=", 
           n, " (for your inputs).\n")
    }
  }
  if (exact)
  {
    parnumber <- match(profile.string, pars[["Label"]])
  }else {
    parnumber <- grep(profile.string, pars[["Label"]])
  }
  if (length(parnumber) <= 0) {
    stop("No parameters matching profile.string='", profile.string, 
         "'", sep = "")
  }
  parlabel <- pars[["Label"]][parnumber]
  if (length(parlabel) > 1) {
    stop("Multiple parameters matching profile.string='", 
         profile.string, "':\n", paste(parlabel, collapse = ", "), 
         "\nYou may need to use 'exact=TRUE'.", sep = "")
  }
  parvec <- as.numeric(pars[pars[["Label"]] == parlabel, models])
  names(parvec)=paste0('replist',1:length(parvec))
  if (verbose) {
    message("Parameter matching profile.string=", profile.string, 
            ": ", parlabel)
    message("Parameter values (after subsetting based on input 'models'): ", 
            paste0(parvec, collapse = ", "))
  }
  par_prior_like_vec <- as.numeric(par_prior_likes[par_prior_likes[["Label"]] == 
                                                     parlabel, models])
  if (all(is.na(par_prior_like_vec))) {
    add_no_prior_line <- FALSE
  }
  par_prior_like_vec[is.na(par_prior_like_vec)] <- 0
  if (all(par_prior_like_vec == 0)) {
    add_no_prior_line <- FALSE
  }
  if (verbose & add_no_prior_line) {
    message("Parameter prior likelihoods: ", paste0(par_prior_like_vec, 
                                                    collapse = ", "))
  }
  if (xlim[1] == "default") 
    xlim <- range(parvec)
  prof.table <- data.frame(t(likelihoods[likelihoods[["Label"]] %in% 
                                           components, models]))
  names(prof.table) <- likelihoods[likelihoods[["Label"]] %in% 
                                     components, ncol(likelihoods)]
  component.labels.good <- rep("", ncol(prof.table))
  for (icol in 1:ncol(prof.table)) {
    ilabel <- which(components == names(prof.table)[icol])
    component.labels.good[icol] <- component.labels[ilabel]
  }
  TOTAL_no_prior <- prof.table[["TOTAL"]] - par_prior_like_vec
  subset <- parvec >= xlim[1] & parvec <= xlim[2]
  for (icol in 1:ncol(prof.table))
  {
    prof.table[, icol] <- prof.table[, icol] - min(prof.table[subset,icol])
  }
  TOTAL_no_prior <- TOTAL_no_prior - min(TOTAL_no_prior)
  if (ymax == "default") {
    ymax <- 1.1 * max(prof.table[subset, ])
  }
  ylim <- c(0, ymax)
  column.max <- apply(prof.table[subset, ], 2, max)
  change.fraction <- column.max/column.max[1]
  include <- change.fraction >= minfraction
  nlines <- sum(include)
  if (verbose) {
    message("Likelihood components showing max change as fraction of total change.\n", 
            "To change which components are included, change input 'minfraction'.\n")
    print(data.frame(frac_change = round(change.fraction, 
                                         4), include = include, label = component.labels.good))
  }
  if (nlines == 0) {
    stop("No components included, 'minfraction' should be smaller.")
  }
  component.labels.used <- component.labels.good[include]
  prof.table <- prof.table[order(parvec), include]
  TOTAL_no_prior <- TOTAL_no_prior[order(parvec)]
  parvec <- parvec[order(parvec)]
  change.fraction <- change.fraction[include]
  if (nlines > 1) {
    if (sort.by.max.change) {
      neworder <- c(1, 1 + order(change.fraction[-1], decreasing = TRUE))
      prof.table <- prof.table[, neworder]
      component.labels.used <- component.labels.used[neworder]
    }
  }
  prof.table <- data.frame(prof.table, TOTAL_no_prior)
  if (add_no_prior_line) {
    component.labels.used <- c(component.labels.used, "Total without prior")
  }
  if (col[1] == "default") {
    col <- rich.colors.short(nlines)
  }
  if (pch[1] == "default") {
    pch <- 1:nlines
  }
  if (add_no_prior_line) {
    col <- c(col, col[1])
    pch <- c(pch, NA)
  }
  lwd <- c(lwd.total, rep(lwd, nlines - 1), switch(add_no_prior_line + 
                                                     1, NULL, lwd))
  cex <- c(cex.total, rep(cex, nlines - 1), switch(add_no_prior_line + 
                                                     1, NULL, cex.total))
  lty <- c(lty.total, rep(lty, nlines - 1), switch(add_no_prior_line + 
                                                     1, NULL, 2))
  if (is.null(profile.label)) {
    if (grepl("steep", parlabel)) {
      profile.label <- "Spawner-recruit steepness (h)"
    }
    if (grepl("R0", parlabel)) {
      profile.label <- paste0("Log of unfished equilibrium recruitment, ", 
                              expression(log(R[0])))
    }
    if (grepl("NatM", parlabel) && grepl("Fem", parlabel)) {
      profile.label <- "Female natural mortality (M)"
    }
    if (grepl("NatM", parlabel) && grepl("Mal", parlabel)) {
      profile.label <- "Male natural mortality (M)"
    }
    if (grepl("LnQ", parlabel)) {
      profile.label <- paste0("Log of catchability, ", 
                              expression(log(q)))
    }
    if (grepl("sigmaR", parlabel)) {
      profile.label <- "SigmaR"
    }
    if (grepl("L_at_Amax", parlabel) && grepl("Fem", parlabel)) {
      profile.label <- "Female length at Amax"
    }
    if (grepl("L_at_Amax", parlabel) && grepl("Mal", parlabel)) {
      profile.label <- "Male length at Amax"
    }
    if (is.null(profile.label)) {
      profile.label <- parlabel
      message("The input profile.label = NULL and the parameter label doesn't ", 
              "correspond to an automatically generated label. ", 
              "Setting profile.label equal to the parameter label.")
    }
  }
  if(!is.null(Xvec))
  {
    Xvec=Xvec[match(names(Xvec),names(parvec))]
    xlim=range(Xvec)
  }
  if(is.null(Xvec)) Xvec=parvec
  plotprofile <- function() {
    plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = profile.label, 
         ylab = ylab, yaxs = yaxs, xaxs = xaxs, ...)
    abline(h = 0, col = "grey")
    if (add_cutoff) {
      abline(h = 0.5 * qchisq(p = cutoff_prob, df = 1), 
             lty = 2)
    }
    matplot(x = Xvec, y = prof.table, type = type, pch = pch, 
            col = col, cex = cex, lty = lty, lwd = lwd, add = T)
    if (legend) {
      legend(legendloc, bty = "n", legend = component.labels.used, 
             lwd = lwd, pt.cex = cex, lty = lty, pch = pch, 
             col = col)
    }
    box()
  }
  if (plot) 
  {
    plotprofile()
  }
  
  if (print)
  {
    save_png(plotinfo = NULL, file = "profile_plot_likelihood.png", 
             plotdir = plotdir, pwidth = pwidth, pheight = pheight, 
             punits = punits, res = res, ptsize = ptsize)
    plotprofile()
    dev.off()
  }
  out <- data.frame(parvec = parvec, prof.table)
  names(out)[1] <- parlabel
  return(invisible(out))
}
function.goodness.fit_SS=function(Rep)
{
  Res.cpue=NA  
  if('cpue'%in%names(Rep))
  {
    if(nrow(Rep$cpue)>0)
    {
      Res.cpue = Rep$cpue%>%
        mutate(residuals = ifelse(is.na(Obs), NA, log(Obs) - log(Exp)))%>%
        pull(residuals)
    }
  }

  Res.length.comps=NA
  if(nrow(Rep$lendbase)>0)
  {
    Res.length.comps = SScompsTA1.8(Rep, fleet = NULL, type = "len", plotit = FALSE)$runs_dat%>%
      mutate(residuals = ifelse(is.na(Obs), NA, log(Obs) -log(Exp)))%>%pull(residuals)
  }
  Res.age.comps=NA
  if(nrow(Rep$agedbase)>0)
  {
    Res.age.comps = SScompsTA1.8(Rep, fleet = NULL, type = "age", plotit = FALSE)$runs_dat%>%
      mutate(residuals = ifelse(is.na(Obs), NA, log(Obs) -log(Exp)))%>%pull(residuals)
  }
  Res.mnwgt=NA  
  if(!is.null(nrow(Rep$mnwgt)))
  {
    if(nrow(Rep$mnwgt)>0)
    {
      Res.mnwgt = Rep$mnwgt%>%
      mutate(residuals = ifelse(is.na(Obs), NA, log(Obs) - log(Exp)))%>%
      pull(residuals)
    }
  }
  Resids=c(Res.cpue,Res.length.comps,Res.age.comps,Res.mnwgt)
  Resids=Resids[!is.na(Resids)]
  
  npar = nrow(Rep$estimated_non_dev_parameters)
  Nobs = length(Resids)
  DF = Nobs - npar
  RMSE = round(100 * sqrt(sum(Resids^2, na.rm = TRUE)/DF),1)
  NegLogLike = Rep$likelihoods_used[match("TOTAL",rownames(Rep$likelihoods_used)),]$values

  AIC= 2*npar+2*NegLogLike
  AICc =AIC+2*npar*(npar + 1) / (Nobs - npar - 1)
  GoodnessFit= data.frame(Stastistic = c("N", "p","DF", "RMSE","AIC", "AICc"), 
                          Value = c(Nobs, npar, DF, RMSE,AIC, AICc))
  return(GoodnessFit)
}
#Francis re weighting of cpue CVs ------------------------------------------------------
Francis.function=function(cipiuis,cvs,mininum.mean.CV=NULL)  
{
  #Get RMSE
  rmseloess=rep(NA,ncol(cvs)-1)
  for(i in 1:length(rmseloess))
  {
    ii=i+1
    trend = cipiuis[,ii]/mean(cipiuis[,ii], na.rm=T)
    sm = loess( log(trend) ~ c(1:nrow(cipiuis)), na.action=na.omit )
    rmseloess[i] <- sqrt( sum(sm$residuals^2)/sm$n)
  }
  cpue_included=names(cipiuis)[-1]
  names(rmseloess)=cpue_included
  if(!is.null(mininum.mean.CV)) for(i in 1:length(rmseloess))rmseloess[i]=max(mininum.mean.CV,rmseloess[i])
  
  #Get adjusted CVs
  sdlogI <- data.frame(cvs[,1],apply(as.data.frame(cvs[,-1]), 2, function(x){ sqrt(log(1+x^2)) } ))
  names(sdlogI) <- names(cvs)
  
  if(ncol(sdlogI)>2){sdlogImean <- apply(sdlogI[,-1],2, function(x){mean(x, na.rm=T)})}else
                    {sdlogImean <- mean(sdlogI[,2],na.rm=T)}
  auxi <- cpue_included
  for(j in 1:length(sdlogImean)) sdlogImean[j] <- max(sdlogImean[j], rmseloess[auxi==names(sdlogI)[j+1]])/sdlogImean[j] 
  sdlogI_not.amended=sdlogI
  
  #these are the ammended CVs used in JABBA
  if(ncol(sdlogI)>2){sdlogI[,-1] <- t(t(sdlogI[,-1])*sdlogImean)}else
                    {sdlogI[,2] <- sdlogI[,2]*sdlogImean} 
  
  #Variance adjustment for SS3
  CV.var.adj=rmseloess
  CV.var.adj[]=NA
  for(v in 1:length(CV.var.adj))
  {
    Mn.CV=mean(cvs[,v+1],na.rm=T)
    if(Mn.CV<=rmseloess[v])
    {
      CV.var.adj[v]=rmseloess[v]-Mn.CV
    }else
    {
      CV.var.adj[v]=0
    }
  }

  return(list(CV.Adjusted=sdlogI,CV.Original=sdlogI_not.amended,CV.var.adj=CV.var.adj))
}
  


# Cryptic mortality -------------------------------------------------------
fn.cryptic=function(yr)
{
  Kls=c('black','steelblue','blue4','brown4','chocolate1','darkolivegreen','chartreuse3')
  names(Kls)=c("Total","Selected","Not sel","Mature","Immature","Mature sel","Mature not sel")
  #Numbers of female at length by year
  Num.fem=Report$natlen%>%
    rename(Beg_Mid='Beg/Mid')%>%
    filter(!Beg_Mid=='M')%>%
    dplyr::select(-c(Area,Bio_Pattern,BirthSeas,Settlement,Platoon,Morph,Seas,Time,Beg_Mid,Era))%>%
    gather(Length,N,-c(Sex,Yr))%>%
    arrange(Yr,Length,Sex)%>%
    mutate(Length=as.numeric(Length))%>%
    filter(Sex==1)
  Yr.LVLs=sort(unique(Num.fem$Yr))
  Yr.kls=colfunc1(length(Yr.LVLs))
  names(Yr.kls)=Yr.LVLs
  
  p_Numbers.at.length.year=Num.fem%>%
    mutate(Yr=factor(Yr,levels=Yr.LVLs))%>%
    ggplot(aes(Length,N,color=Yr))+
    geom_line()+
    theme_PA(leg.siz=6)+ylab('Number of females')+
    scale_color_manual(values = Yr.kls)+
    guides(color = guide_legend(ncol = 2))
  
  Num.fem=Num.fem%>%
    filter(Yr==yr)
  
  #Maturity at length
  Fem.mat=Report$biology%>%
    rename(Length=Len_lo)
  p_maturity.ogive=Fem.mat%>%
    ggplot(aes(Length,Mat))+
    geom_line()+
    theme_PA()
  
  #Selectivity at length
  Fem.sel=Report$sizeselex%>%
    filter(Factor=='Lsel' & Sex==1)%>%
    dplyr::select(-c(Factor,Label,Sex))%>%
    gather(Length,N,-c(Fleet,Yr))%>%
    arrange(Yr,Length,Fleet)%>%
    mutate(Length=as.numeric(Length),
           Length=Length-2.5,
           Fleet=as.character(Fleet))
  Yr.LVLs=sort(unique(Fem.sel$Yr))
  Yr.kls=colfunc1(length(Yr.LVLs))
  names(Yr.kls)=Yr.LVLs
  Fem.sel=Fem.sel%>%
    filter(Yr==as.numeric(substr(Last.yr.ktch,1,4)))
  
  p_sel.by.fleet=Fem.sel%>%
    mutate(Yr=factor(Yr,levels=Yr.LVLs))%>%
    ggplot(aes(Length,N,color=Fleet))+
    geom_line(linewidth=1.25)+
    theme_PA()+ylab('Selectivity')+
    theme(legend.position = 'top')
  
  
  
  #plot cryptic mature biomass by fleet
  N.pop=Num.fem%>%arrange(Length)
  Mat=Fem.mat%>%dplyr::select(Length,Mat)%>%arrange(Length)
  
  dis.fleets=match(c("Northern.shark","Southern.shark_1","Survey"),Report$FleetNames)
  names(dis.fleets)=c("Northern.shark","Southern.shark_1","Survey")
  dis.fleets=dis.fleets[!is.na(dis.fleets)]
  p_list=vector('list',length(dis.fleets))
  p1_list=p_list
  for(kek in 1:length(dis.fleets))
  {
    Sel=Fem.sel%>%filter(Fleet==dis.fleets[kek])%>%arrange(Length)
    
    #Numbers at length in population
    DAT=N.pop%>%
      mutate(N=N/max(N))%>%
      dplyr::select(Length,N)%>%
      mutate(Type='Population size')
    DAT=rbind(DAT,Sel%>%
                dplyr::select(Length,N)%>%
                mutate(Type='Selectivity'))
    DAT=rbind(DAT,Mat%>%
                rename(N=Mat)%>%
                dplyr::select(Length,N)%>%
                mutate(Type='Maturity'))%>%
      mutate(Type=factor(Type,levels=c('Population size','Maturity','Selectivity')))
    
    p=DAT%>%
      ggplot(aes(Length,N,color=Type))+
      geom_point(size=3)+geom_line(linewidth=1.05)+theme_PA(Ttl.siz=14)+
      theme(legend.title = element_blank(),
            legend.position = 'bottom')
    
    #Calculate mature cryptic
    N.sel.mat=full_join(N.pop%>%
                          dplyr::select(Length,N),
                        Sel%>%
                          dplyr::select(Length,N)%>%
                          rename(Select=N),
                        by=c('Length'))%>%
      full_join(Mat,by=c('Length'))%>%
      mutate(N.selected=N*Select,
             N.not.selected=N-N.selected,
             N.mature=N*Mat,
             N.immature=N-N.mature,
             N.mature.selected=N.selected*Mat,
             N.mature.not.selected=N.not.selected*Mat)%>%
      arrange(Length)
    id.sel.peak=match(max(N.sel.mat$Select),N.sel.mat$Select)
    Mature.cryptic=round(100*sum(N.sel.mat$N.mature.not.selected[id.sel.peak:nrow(N.sel.mat)])/sum(N.sel.mat$N.mature),2)
    
    p=p+
      labs(title=paste0(names(dis.fleets)[kek],'      ',Mature.cryptic,'%',' cryptic mature female'))
    
    p_list[[kek]]=p
    
    p1=N.sel.mat%>%
      dplyr::select(-c(Select,Mat))%>%
      gather(Type,Value,-Length)%>%
      mutate(Type=case_when(
        Type=="N"~"Total",
        Type=="N.selected"~"Selected",
        Type=="N.not.selected"~"Not sel",
        Type=="N.mature"~"Mature",
        Type=="N.immature"~"Immature",
        Type=="N.mature.selected"~"Mature sel",
        Type=="N.mature.not.selected"~"Mature not sel"),
        Type=factor(Type,levels=names(Kls)))%>%
      ggplot(aes(Length,Value,color=Type))+
      geom_point(size=3)+
      geom_line(linewidth=1.05)+
      theme_PA(Ttl.siz=14,leg.siz=10)+
      theme(legend.title = element_blank())+
      scale_color_manual(values = Kls)+
      geom_line(data=Sel%>%rename(Value=N)%>%mutate(Type='',Value=Value*max(N.sel.mat$N)),
                aes(Length,Value,color=Type),linewidth=1.25)+
      geom_line(data=Mat%>%rename(Value=Mat)%>%mutate(Type='',Value=Value*max(N.sel.mat$N)),
                aes(Length,Value,color=Type),linetype='longdash',linewidth=1.25)+
      labs(title=paste0(names(dis.fleets)[kek],'      ',Mature.cryptic,'%',' cryptic mature female'),
           caption='solid=Selectivity\n dashed=Maturity')
    p1_list[[kek]]=p1
    
  }
  
  return(list(p.criptic=p_list,p.criptic1=p1_list,p_Numbers.at.length.year=p_Numbers.at.length.year,
              p_maturity.ogive=p_maturity.ogive,p_sel.by.fleet=p_sel.by.fleet))
}


# Some random functions-------------------------------------------------------------------------
fn.check.SS.sel.used=function(d,check.fleet=NULL)
{
  p=d%>%
    gather(TL,Sel,-c(Fleet,Yr,Sex))%>%
    mutate(TL=as.numeric(TL),
           Sel=as.numeric(Sel),
           Sex=ifelse(Sex==1,'F','M'),
           Year.sex=paste(Yr,Sex))
  
  p1=p%>%
    ggplot(aes(TL,Sel,color=Year.sex))+
    geom_line()+
    geom_point(aes(shape=Sex))+
    facet_wrap(~Fleet)+theme_PA(leg.siz=7)+
    theme(legend.position = 'top',legend.title = element_blank())+
    guides(color=guide_legend(nrow=3,byrow=TRUE))
  
  p2=NULL
  if(!is.null(check.fleet))
  {
    p2=p%>%
      filter(Fleet%in%check.fleet)%>%
      mutate(Yr=as.character(Yr))%>%
      ggplot(aes(TL,Sel,color=Sex))+
      geom_line()+
      facet_wrap(~Yr)+theme_PA(leg.siz=10)+
      theme(legend.position = 'top',legend.title = element_blank())+
      guides(color=guide_legend(nrow=1,byrow=TRUE))+ggtitle(paste('Fleet',check.fleet))
  }
  p.list=list(p1=p1)
  if(!is.null(p2)) p.list=list(p1=p1,p2=p2)
  return(p.list)
}
see.SS3.length.comp.matrix=function(dd)
{
  return(dd%>%
           gather(Size,n,-c(year ,Fleet))%>%
           mutate(Fleet=as.character(Fleet),
                  Size=as.numeric(substr(Size,2,10)))%>%
           ggplot(aes(Size,n))+
           geom_line(aes(color=Fleet))+
           facet_wrap(~year)+
           theme_PA()+theme(legend.position = 'top')+xlab('TL (cm)'))
}