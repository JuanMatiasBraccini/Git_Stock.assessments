#---Run models --------------------------------------
#18.1 Execute each COM
n.catch.only=length(catch.only)  
Catch_only=vector('list',n.catch.only)
names(Catch_only)=catch.only
tic("timer")
for(w in 1:length(Catch_only))
{
  #1. DBSRA assessment (Dick and MAcCall (2011))
  # summary of method: http://toolbox.frdc.com.au/wp-content/uploads/sites/19/2020/07/DBSRA3.html
  #note: this was parallelised to improve computation time
  if(names(Catch_only)[w]=="DBSRA")     #parallelised: 0.004 secs per iteration-species-scenario (otherwise 0.025 secs)
  {
    dummy.store=vector('list',N.sp)     
    names(dummy.store)=Keep.species
    dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
      dummy.store.probs.B.Bmsy=dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=
      dummy.store.accept.rate=dummy.store.ensemble=dummy.store
    for(i in 1:length(dummy.store))  
    {
      if(names(dummy.store)[i]%in%Catch.only.species)    
      {
        #Catch
        ktch=ktch.combined%>%
          filter(Name==names(dummy.store)[i])
        all.years=seq(min(ktch$finyear),max(ktch$finyear))
        misn.yr=all.years[which(!all.years%in%ktch$finyear)]
        if(length(misn.yr)>0)
        {
          ktch=rbind(ktch,ktch[length(misn.yr),]%>%mutate(finyear=misn.yr,Tonnes=0))%>%arrange(finyear)
        }
        
        #future catches
        outktch=ktch$Tonnes
        add.ct.future=NULL
        if('Catch_only'%in%future.models)
        {
          if(catches.futures=="constant.last.n.yrs")
          {
            NN=nrow(ktch)
            add.ct.future=ktch[1:years.futures,]%>%
              mutate(finyear=ktch$finyear[NN]+(1:years.futures),
                     Tonnes=rep(mean(ktch$Tonnes[(NN-(n.last.catch.yrs-1)):NN]),years.futures))
            outktch=c(outktch,add.ct.future$Tonnes)
            #ktch=rbind(ktch,add.ct.future)
          } 
        }

        
        #working directory
        this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                      capitalize(List.sp[[i]]$Name),"/",AssessYr,"/DBSRA",sep='')
        if(!dir.exists(this.wd))dir.create(this.wd)
        
        #Scenarios
        Scens=List.sp[[i]]$Sens.test$DBSRA%>%
          mutate(Species=capitalize(names(dummy.store)[i]))
        Store.sens=vector('list',nrow(Scens))
        names(Store.sens)=Scens$Scenario
        this.wd1=this.wd
        
        Out.Scens=Scens%>%
          mutate(M.dist=NA,M.CV=NA,fmsym.dist=NA,fmsym.CV=NA,bmsyk.dist=NA,bmsyk.mean=NA,bmsyk.CV=NA,
                 b1k.dist=NA,b1k.low=NA,b1k.up=NA,btk.dist=NA,btk.low=NA,btk.up=NA)
        Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.probs.B.Bmsy=Out.f.series=
          Out.B.Bmsy=Out.F.Fmsy=Out.accept.rate=Out.ensemble=vector('list',length(Store.sens))
        names(Out.ensemble)=names(Store.sens)
        
        cl <- makeCluster(detectCores()-1)
        registerDoSNOW(cl)
        for(s in 1:length(Store.sens))
        {
          print(paste("___________","DBSRA Scenario",Scens$Scenario[s],"___________",names(dummy.store)[i]))
          
          this.wd=paste(this.wd1,names(Store.sens)[s],sep='/')
          if(!dir.exists(this.wd))dir.create(this.wd)
          
          AgeMat=Scens$AgeMat[s]
          
          Mmean=Scens$Mmean[s]  
          Msd=Scens$Msd[s]
          Mcv=Msd/Mmean
          M.dist="lnorm"
          
          Klow=Scens$Klow[s]
          Kup=Scens$Kup[s]
          
          fmsy.m=Scens$fmsy.m[s]
          fmsym.dist="lnorm"
          fmsym.CV=0.1
          
          bmsyk.mean=Scens$bmsyk[s]
          bmsyk.dist="beta"
          bmsyk.CV=0.1
          
          btk.dist="unif"
          Btklow=List.sp[[i]]$FINALBIO[1]
          Btkup=List.sp[[i]]$FINALBIO[2]
          
          b1k.dist="unif"
          b1k.low=List.sp[[i]]$STARTBIO[1]
          b1k.up=List.sp[[i]]$STARTBIO[2]
          
          Depletion.year=max(ktch$finyear)
          # if(names(dummy.store)[i]=="milk shark")
          # {
          #   Depletion.year=1991  #set to end of catch period to allow convergence
          #   Btklow=0.2
          # }
          
          #Run model
          Store.sens[[s]]=apply.DBSRA_tweeked(year=ktch$finyear,
                                              catch=ktch$Tonnes,
                                              catchCV=NULL,  
                                              catargs=list(dist="none",low=0,up=Inf,unit="MT"),  #catch CV not available
                                              agemat=AgeMat,
                                              k=list(low=Klow,up=Kup,tol=0.01,permax=1000),
                                              b1k=list(dist=b1k.dist,low=b1k.low,up=b1k.up,mean=1,sd=0.1),  #mean and sd not used if 'unif'
                                              btk=list(dist="unif",low=Btklow,up=Btkup,mean=1,sd=0.1,refyr=Depletion.year),  #reference year
                                              fmsym=list(dist="lnorm",low=0.1,up=2,mean=log(fmsy.m),sd=fmsym.CV), # Cortes & Brooks 2018. Low and up not used if 'lnorm'  
                                              bmsyk=list(dist="beta",low=0.05,up=0.95,mean=bmsyk.mean,sd=bmsyk.CV),  
                                              M=list(dist="lnorm",low=0.001,up=1,mean=log(Mmean),sd=Mcv),
                                              graph=c(),
                                              nsims=Scens$Sims[s],
                                              grout=1,
                                              WD='C:/DummyDBSRA',  #output to dummy folder because OneDrive backup stuff things up
                                              outfile="Appendix_fit",
                                              Projections=add.ct.future)
          
          #Store acceptance rate
          Accept.tab=table(Store.sens[[s]]$output$Values$ll)
          Out.accept.rate[[s]]=data.frame(
            Acceptance=round(100*Accept.tab[2]/sum(Accept.tab),2),
            Scenario=Scens$Scenario[s])
          
          #Store scenarios
          Out.Scens$M.dist[s]=M.dist
          Out.Scens$M.CV[s]=Mcv
          Out.Scens$fmsym.dist[s]=fmsym.dist
          Out.Scens$fmsym.CV[s]=fmsym.CV
          Out.Scens$bmsyk.dist[s]=bmsyk.dist
          Out.Scens$bmsyk.mean[s]=bmsyk.mean
          Out.Scens$bmsyk.CV[s]=bmsyk.CV
          Out.Scens$b1k.dist[s]=b1k.dist
          Out.Scens$b1k.low[s]=b1k.low
          Out.Scens$b1k.up[s]=b1k.up
          Out.Scens$btk.dist[s]=btk.dist
          Out.Scens$btk.low[s]=Btklow
          Out.Scens$btk.up[s]=Btkup
          
          #Store estimates
          d1=Store.sens[[s]]$output$Parameters
          d1=d1[,grep(paste(c('Median','2.5%','97.5%'),collapse='|'),names(d1))]
          names(d1)=c("Median","Lower.95","Upper.95")
          d2=Store.sens[[s]]$output$Estimates
          d2=d2[,grep(paste(c('Median','2.5%','97.5%'),collapse='|'),names(d2))]
          names(d2)=c("Median","Lower.95","Upper.95")
          d1=rbind(d2,d1)
          d1=d1%>%
            mutate(Parameter=rownames(d1),
                   Model='DBSRA',
                   Scenario=names(Store.sens)[s])%>%
            relocate(Model,Scenario,Parameter,Lower.95)
          Out.estimates[[s]]=d1
          
          #Store trajectories 
          dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                            mods=names(Catch_only)[w],
                                            Type='Depletion',
                                            scen=Scens$Scenario[s],
                                            Katch=outktch)
          Out.rel.biom[[s]]=dummy$Dat
          Out.probs.rel.biom[[s]]=dummy$Probs
          
          
          dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                            mods=names(Catch_only)[w],
                                            Type='F.series',
                                            scen=Scens$Scenario[s],
                                            Katch=outktch)
          Out.f.series[[s]]=dummy$Dat
          
          dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                            mods=names(Catch_only)[w],
                                            Type='B.Bmsy',
                                            scen=Scens$Scenario[s],
                                            Katch=outktch)
          Out.B.Bmsy[[s]]=dummy$Dat
          Out.probs.B.Bmsy[[s]]=dummy$Probs
          
          dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                            mods=names(Catch_only)[w],
                                            Type='F.Fmsy',
                                            scen=Scens$Scenario[s],
                                            Katch=outktch)
          Out.F.Fmsy[[s]]=dummy$Dat
          rm(dummy)
          
          #Display Priors vs Posteriors for base case scenario (S1) 
          if(Scens$Scenario[s]=='S1')
          {
            par.list=c('m','fmsym','b1k','btk','bmsyk','k')
            dummy.prior=Store.sens[[s]]$input
            names(dummy.prior)=tolower(names(dummy.prior))
            dummy.post=Store.sens[[s]]$output$Values%>%filter(ll==1) #accepted runs
            colnames(dummy.post)=tolower(colnames(dummy.post))
            out=vector('list',length(par.list))
            for(p in 1:length(par.list))
            {
              prior.dis=dummy.prior[[par.list[p]]] 
              if(par.list[p]=='k')
              {
                prior.dis$dist="unif" 
              }
              out[[p]]=rbind(data.frame(Distribuion="Prior",
                                        Value=fn.prior(d=prior.dis)),   
                             data.frame(Distribuion="Posterior",
                                        Value=dummy.post[[par.list[p]]]))%>%
                mutate(Parameter=par.list[p])
            }
            rm(dummy.prior,dummy.post)
            fn.show.density(d=do.call(rbind,out),NCOL=2)
            ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                         capitalize(names(dummy.store)[i]),"/",AssessYr,"/",'DBSRA',"/Prior.and.posterior.tiff",sep=''),
                   width = 12,height = 14, dpi = 300, compression = "lzw")
          }
          
          #Store quantities for ensemble model  
          Years=Store.sens[[s]]$output$Years
          Out.ensemble[[s]]=list(Depletion=Store.sens[[s]]$output$Depletion.traj[1:length(Years)],
                                 B.Bmsy=Store.sens[[s]]$output$B.Bmsy[1:length(Years)],
                                 F.Fmsy=Store.sens[[s]]$output$F.Fmsy[1:length(Years)],
                                 MSY=data.frame(Par='MSY',
                                                Value=Store.sens[[s]]$output$Values$MSY[which(Store.sens[[s]]$output$Values$ll==1)],
                                                Scenario=Scens$Scenario[s]))
        }
        stopCluster(cl)
        #output quantities for ensemble
        Out.ens=list(Depletion=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'Depletion')),
                     B.Bmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'B.Bmsy')),
                     F.Fmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'F.Fmsy')),
                     MSY=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'MSY')))
        if(nrow(Scens)>1)
        {
          for(f in 1:length(Out.ens))   
          {
            ddmi=Out.ens[[f]]
            TAB=table(sub("\\..*", "", rownames(ddmi)))
            Out.ens[[f]]=ddmi[sample(1:nrow(ddmi),size=max(TAB),replace =FALSE),]
          }
        }
        
        #output scenarios
        #dummy.store[[i]]=Store.sens
        Out.Scens=Out.Scens%>%
          dplyr::select(Species,Scenario,AgeMat,M.dist,Mmean,M.CV,Klow,Kup,
                        fmsym.dist,fmsy.m,fmsym.CV,bmsyk.dist,bmsyk.mean,bmsyk.CV,
                        b1k.dist,b1k.low,b1k.up,btk.dist,btk.low,btk.up)%>%
          mutate(Mmean=round(Mmean,2),
                 M.CV=round(M.CV,2),
                 Klow=round(Klow),
                 Kup=round(Kup),
                 fmsym.dist=ifelse(fmsym.dist=="lnorm","Lognormal",fmsym.dist),
                 M.dist=ifelse(M.dist=="lnorm","Lognormal",M.dist),
                 bmsyk.dist=capitalize(bmsyk.dist),
                 b1k.dist=ifelse(b1k.dist=="unif","Uniform",b1k.dist),
                 btk.dist=ifelse(btk.dist=="unif","Uniform",btk.dist))
        
        dummy.store.sens.table[[i]]=Out.Scens
        dummy.store.accept.rate[[i]]=do.call(rbind,Out.accept.rate)
        dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
        dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
        dummy.store.probs.B.Bmsy[[i]]=Out.probs.B.Bmsy
        dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
        dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
        dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
        dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
        dummy.store.ensemble[[i]]=Out.ens
      }
     } #end i
    
    #Catch_only[[w]]=dummy.store    #too big an object, cannot store whole model runs due to memory size limitations
    Catch_only[[w]]$sens.table=dummy.store.sens.table
    Catch_only[[w]]$accept.rate=dummy.store.accept.rate
    Catch_only[[w]]$estimates=dummy.store.estimates
    Catch_only[[w]]$rel.biom=dummy.store.rel.biom
    Catch_only[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
    Catch_only[[w]]$probs.B.Bmsy=dummy.store.probs.B.Bmsy
    Catch_only[[w]]$f.series=dummy.store.f.series
    Catch_only[[w]]$B.Bmsy=dummy.store.B.Bmsy
    Catch_only[[w]]$F.Fmsy=dummy.store.F.Fmsy
    Catch_only[[w]]$ensemble=dummy.store.ensemble
    
    rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,
       dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,
       dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.accept.rate,
       dummy.store.ensemble,dummy.store.probs.B.Bmsy)
  }
  
  #2. CMSY    
  #summary of method: http://toolbox.frdc.com.au/wp-content/uploads/sites/19/2021/04/CMSY.html
  #note: no need of parallelisation
  if(names(Catch_only)[w]=="CMSY")     #takes 0.008 secs per iteration-species-scenario
  {
    dummy.store=vector('list',N.sp)     
    names(dummy.store)=Keep.species 
    dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
      dummy.store.probs.B.Bmsy=dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=
      dummy.store.accept.rate=dummy.store.ensemble=dummy.store
    for(i in 1:length(dummy.store))  
    {
      if(names(dummy.store)[i]%in%Catch.only.species) 
      {
        this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                      capitalize(List.sp[[i]]$Name),"/",AssessYr,"/CMSY",sep='')
        if(!dir.exists(this.wd))dir.create(this.wd)
        
        #Catch
        ktch=ktch.combined%>%
          filter(Name==names(dummy.store)[i])
        all.years=seq(min(ktch$finyear),max(ktch$finyear))
        misn.yr=all.years[which(!all.years%in%ktch$finyear)]
        if(length(misn.yr)>0)
        {
          ktch=rbind(ktch,ktch[length(misn.yr),]%>%mutate(finyear=misn.yr,Tonnes=0))%>%arrange(finyear)
        }
        
        #future catches
        outktch=ktch$Tonnes
        add.ct.future=NULL
        if('Catch_only'%in%future.models)
        {
          if(catches.futures=="constant.last.n.yrs")
          {
            NN=nrow(ktch)
            add.ct.future=ktch[1:years.futures,]%>%
              mutate(finyear=ktch$finyear[NN]+(1:years.futures),
                     Tonnes=rep(mean(ktch$Tonnes[(NN-(n.last.catch.yrs-1)):NN]),years.futures))
            outktch=c(outktch,add.ct.future$Tonnes)
            #ktch=rbind(ktch,add.ct.future)
          }
        }

        
        year=ktch$finyear
        catch=ktch$Tonnes
        
        #Scenarios
        Scens=List.sp[[i]]$Sens.test$CMSY%>%
          mutate(Species=capitalize(names(dummy.store)[i]))
        Store.sens=vector('list',nrow(Scens))
        names(Store.sens)=Scens$Scenario
        this.wd1=this.wd
        
        Out.Scens=Scens%>%   
          mutate(Bo.low=NA,Bo.hi=NA,Bf.low=NA,Bf.hi=NA,r.low=NA,r.up=NA)
        Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.probs.B.Bmsy=Out.f.series=
          Out.B.Bmsy=Out.F.Fmsy=Out.accept.rate=Out.ensemble=vector('list',length(Store.sens))
        names(Out.ensemble)=names(Store.sens)
        for(s in 1:length(Store.sens))
        {
          print(paste("___________","CMSY Scenario",Scens$Scenario[s],"___________",names(dummy.store)[i]))
          this.wd=paste(this.wd1,names(Store.sens)[s],sep='/')
          if(!dir.exists(this.wd))dir.create(this.wd)
          
          #Priors
          RES=RESILIENCE[[i]]
          if(Scens$r.prob.min[s]==0)
          {
            r.range=NA
            k.range=NA
            Bf.low=NA
            Bf.hi=NA
          }else
          {
            Mn=max(min(Scens$r[s],Max.r.value),Min.r.value)  #some life history pars yield unrealistically high r
            r.range=quantile(rnorm(1e3,
                                   mean=Mn,
                                   sd=Scens$r.sd[s]),
                             probs=c(Scens$r.prob.min[s],Scens$r.prob.max[s]))

            k.range=c(Scens$Klow[s],Scens$Kup[s])
            
            Bf.low=List.sp[[i]]$FINALBIO[1]
            Bf.hi=List.sp[[i]]$FINALBIO[2]
            
            #need to modify bound a bit to allow enough combos (i.e. convergence)
            if(names(dummy.store)[i]%in%c("milk shark","narrow sawfish","snaggletooth","zebra shark")) r.range[1]=0.1
            if(names(dummy.store)[i]%in%c("pigeye shark","smooth hammerhead")) k.range[2]=k.range[1]*3
            if(names(dummy.store)[i]%in%c("pigeye shark","sandbar shark")) r.range[2]=r.range[2]*1.2
            if(names(dummy.store)[i]%in%c("zebra shark")) r.range[2]=min(r.range[2],Scens$r[s]*1.1)
            
            r.range[1]=max(r.range[1],Min.r.value)  #minimum level of recruitment
            r.range[2]=min(r.range[2],Max.r.value)  #maximum level biological possible
            
          }
          Proc.error=Scens$Proc.error[s]
          
          Bo.low=List.sp[[i]]$STARTBIO[1]
          Bo.hi=List.sp[[i]]$STARTBIO[2]
          
          
          #Run model
          Store.sens[[s]]=apply.CMSY(year=year,
                                     catch=catch,
                                     r.range=r.range,
                                     k.range=k.range,
                                     Bo.low=Bo.low,
                                     Bo.hi=Bo.hi,
                                     Bf.low=Bf.low,
                                     Bf.hi=Bf.hi,
                                     outfile=paste(this.wd,'Appendix_fit',sep='/'),
                                     nsims=Scens$Sims[s],
                                     Proc.error=Proc.error,
                                     RES=RES,
                                     Projections=add.ct.future)  
          
          #Store acceptance rate
          Out.accept.rate[[s]]=data.frame(
            Acceptance=Store.sens[[s]]$output$acceptance.rate,
            Scenario=Scens$Scenario[s])
          
          #Store scenarios
          Out.Scens$Bo.low[s]=Bo.low     
          Out.Scens$Bo.hi[s]=Bo.hi
          Out.Scens$Bf.low[s]=Bf.low
          Out.Scens$Bf.hi[s]=Bf.hi
          Out.Scens$r.low[s]=r.range[1]
          Out.Scens$r.up[s]=r.range[2]
          Out.Scens$Klow[s]=k.range[1]
          Out.Scens$Kup[s]=k.range[2]
          
          #Store estimates
          d1=Store.sens[[s]]$output$Statistics$output  
          d1=d1[,grep(paste(c('50%','2.5%','97.5%'),collapse='|'),colnames(d1))]%>%
            data.frame
          d1=d1[,-grep('Perc',colnames(d1))]
          colnames(d1)=c("Lower.95","Median","Upper.95")
          d1=d1%>%
            mutate(Parameter=rownames(d1),
                   Model='CMSY',
                   Scenario=names(Store.sens)[s])%>%
            relocate(Model,Scenario,Parameter,Lower.95)
          
          add.Bmsy=d1[1,]%>%
            mutate(Parameter='Bmsy',
                   Lower.95=quantile(Store.sens[[s]]$output$Bmsy,prob=0.025),
                   Median=median(Store.sens[[s]]$output$Bmsy),
                   Upper.95=quantile(Store.sens[[s]]$output$Bmsy,prob=0.975))
          rownames(add.Bmsy)='Bmsy'
          
          Out.estimates[[s]]=rbind(d1,add.Bmsy)
          
          #Store trajectories
          dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                            mods=names(Catch_only)[w],
                                            Type='Depletion',
                                            scen=Scens$Scenario[s],
                                            Katch=outktch)
          Out.rel.biom[[s]]=dummy$Dat
          Out.probs.rel.biom[[s]]=dummy$Probs
          
          dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                            mods=names(Catch_only)[w],
                                            Type='F.series',
                                            scen=Scens$Scenario[s],
                                            Katch=outktch)
          Out.f.series[[s]]=dummy$Dat
          
          dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                            mods=names(Catch_only)[w],
                                            Type='B.Bmsy',
                                            scen=Scens$Scenario[s],
                                            Katch=outktch)
          Out.B.Bmsy[[s]]=dummy$Dat
          Out.probs.B.Bmsy[[s]]=dummy$Probs
          
          dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                            mods=names(Catch_only)[w],
                                            Type='F.Fmsy',
                                            scen=Scens$Scenario[s],
                                            Katch=outktch)
          Out.F.Fmsy[[s]]=dummy$Dat
          rm(dummy)
          
          #Store quantities for ensemble model
          Years=Store.sens[[s]]$output$Years
          F.Fmsy=Store.sens[[s]]$output$F.Fmsy[,1:length(Years)]%>%
            data.frame
          colnames(F.Fmsy)=Years
          Out.ensemble[[s]]=list(Depletion=Store.sens[[s]]$output$Depletion.traj[1:length(Years)], 
                                 B.Bmsy=Store.sens[[s]]$output$B.Bmsy[1:length(Years)],
                                 F.Fmsy=F.Fmsy,
                                 MSY=data.frame(Par='MSY',
                                                Value=Store.sens[[s]]$output$Statistics$traject[,'r']*Store.sens[[s]]$output$Statistics$traject[,'K']/4,
                                                Scenario=Scens$Scenario[s]))
          
          
          #Display Priors vs Posteriors for base case scenario (S1)
          if(Scens$Scenario[s]=='S1')
          {
            par.list=c("r.range","k.range","stb.low","stb.hi","endb.low","endb.hi")
            dummy.prior=Store.sens[[s]]$input[par.list]
            names(dummy.prior)=tolower(names(dummy.prior))
            dummy.prior$stb=c(dummy.prior$stb.low,dummy.prior$stb.hi)
            dummy.prior$endb=c(dummy.prior$endb.low,dummy.prior$endb.hi)
            dummy.prior=within(dummy.prior, rm(stb.low,stb.hi,endb.low,endb.hi))
            for(d in 1:length(dummy.prior))
            {
              dummy.prior[[d]]=list(dist='unif',
                                    low=dummy.prior[[d]][1],
                                    up=dummy.prior[[d]][2])
            }
            dummy.post=as.data.frame(Store.sens[[s]]$output$Statistics$deplet)
            dummy.post=dummy.post[,c(match(c('r','K'),names(dummy.post)),1,ncol(dummy.post)-3)]
            colnames(dummy.post)=c('r','K','Initial depletion','Final depletion')
            names(dummy.prior)=names(dummy.post)
            
            pars=colnames(dummy.post)
            out=vector('list',length(pars))
            for(p in 1:length(pars))
            {
              out[[p]]=rbind(data.frame(Distribuion="Prior",
                                        Value=fn.prior(d=dummy.prior[[pars[p]]])),
                             data.frame(Distribuion="Posterior",
                                        Value=dummy.post[[pars[p]]]))%>%
                mutate(Parameter=pars[p])
            }
            rm(dummy.prior,dummy.post)
            
            fn.show.density(d=do.call(rbind,out),NCOL=2)
            ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                         capitalize(names(dummy.store)[i]),"/",AssessYr,"/",'CMSY',"/Prior.and.posterior.tiff",sep=''),
                   width = 12,height = 14, dpi = 300, compression = "lzw")
          }
        }
        #output quantities for ensemble
        Out.ens=list(Depletion=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'Depletion')),
                     B.Bmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'B.Bmsy')),
                     F.Fmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'F.Fmsy')),
                     MSY=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'MSY')))
        if(nrow(Scens)>1)
        {
          for(f in 1:length(Out.ens))   
          {
            ddmi=Out.ens[[f]]
            TAB=table(sub("\\..*", "", rownames(ddmi)))
            Out.ens[[f]]=ddmi[sample(1:nrow(ddmi),size=max(TAB),replace =FALSE),]
          }
        }
        
        #output scenarios
        Out.Scens=Out.Scens%>%
          dplyr::select(Species,Scenario,r.low,r.up,Klow,Kup,Bo.low,Bo.hi,Bf.low,Bf.hi,Proc.error)%>%
          mutate(Klow=round(Klow),
                 Kup=round(Kup),
                 r.low=round(r.low,3),
                 r.up=round(r.up,3))
        
        dummy.store.sens.table[[i]]=Out.Scens
        dummy.store.accept.rate[[i]]=do.call(rbind,Out.accept.rate)
        dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
        dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
        dummy.store.probs.B.Bmsy[[i]]=Out.probs.B.Bmsy
        dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
        dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
        dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
        dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
        dummy.store.ensemble[[i]]=Out.ens
      }
     }#end i
    
    Catch_only[[w]]$sens.table=dummy.store.sens.table
    Catch_only[[w]]$accept.rate=dummy.store.accept.rate
    Catch_only[[w]]$estimates=dummy.store.estimates
    Catch_only[[w]]$rel.biom=dummy.store.rel.biom
    Catch_only[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
    Catch_only[[w]]$probs.B.Bmsy=dummy.store.probs.B.Bmsy
    Catch_only[[w]]$f.series=dummy.store.f.series
    Catch_only[[w]]$B.Bmsy=dummy.store.B.Bmsy
    Catch_only[[w]]$F.Fmsy=dummy.store.F.Fmsy
    Catch_only[[w]]$ensemble=dummy.store.ensemble
    
    rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,
       dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,
       dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.accept.rate,
       dummy.store.ensemble,dummy.store.probs.B.Bmsy)
  }
  
  #3. OCOM assessment (Zhou et al (2018))
  #note: not used as it doesn't allow for lightly depleted stocks and catch reductions
  # due to effort reduction rather than abundance. Also cannot define depletion ranges
  if(names(Catch_only)[w]=="OCOM")     #takes 20 secs per species (1e4 iterations) 
  {
    dummy.store=vector('list',N.sp)     
    names(dummy.store)=Keep.species             
    system.time({for(i in 1:length(dummy.store))  
    {
      print(paste("OCOM ","--",names(dummy.store)[i]))
      this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(List.sp[[i]]$Name),"/",AssessYr,"/OCOM",sep='')
      if(!dir.exists(this.wd))dir.create(this.wd)
      ktch=ktch.combined%>%
        filter(Name==names(dummy.store)[i])
      dummy.store[[i]] <-apply.OCOM(year=ktch$finyear,
                                    catch=ktch$Tonnes,
                                    M=mean(apply(store.species.M[[i]],2,mean,na.rm=T)),
                                    outfile=paste(this.wd,'Appendix_fit',sep='/'))
    }})
    Catch_only[[w]]=dummy.store
    rm(dummy.store)
  }
  
  #4. JABBA - catch only (Winker et al 2018)   
  #summary of method: https://github.com/jabbamodel/JABBA
  #                   https://www.youtube.com/watch?v=9icz7S4VpPU
  #      no need of parallelisation
  if(names(Catch_only)[w]=="JABBA")     #takes 0.002 secs per iteration-species-scenario
  {
    dummy.store=vector('list',N.sp)     
    names(dummy.store)=Keep.species
    dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
      dummy.store.probs.B.Bmsy=dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=
      dummy.store.Kobe.probs=dummy.store.ensemble=dummy.store
    for(i in 1:length(dummy.store))  
    {
      if(names(dummy.store)[i]%in%Catch.only.species) 
      {
        this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                      capitalize(List.sp[[i]]$Name),"/",AssessYr,"/JABBA",sep='')
        if(!dir.exists(this.wd))dir.create(this.wd)
        
        #Catch
        ktch=ktch.combined%>%
          filter(Name==names(dummy.store)[i])%>%
          rename(Year=finyear,
                 Total=Tonnes)%>%
          ungroup()%>%
          dplyr::select(Year,Total)%>%
          arrange(Year)%>%
          data.frame
        all.years=seq(min(ktch$Year),max(ktch$Year))
        misn.yr=all.years[which(!all.years%in%ktch$Year)]
        if(length(misn.yr)>0)
        {
          ktch=rbind(ktch,ktch[length(misn.yr),]%>%mutate(Year=misn.yr,Total=0))%>%arrange(Year)
        }
        if(names(dummy.store)[i]=='milk shark') ktch$Total[1:10]=1 #convergence issues with first 10 years for milk shark, but catch is 0 anyways
        
        #future catches
        outktch=ktch$Total
        add.ct.future=NULL
        if('Catch_only'%in%future.models)
        {
          if(catches.futures=="constant.last.n.yrs")
          {
            NN=nrow(ktch)
            add.ct.future=ktch[1:years.futures,]%>%
              mutate(Year=ktch$Year[NN]+(1:years.futures),
                     Total=rep(mean(ktch$Total[(NN-(n.last.catch.yrs-1)):NN]),years.futures))
            outktch=c(outktch,add.ct.future$Total)
            #ktch=rbind(ktch,add.ct.future)
          }
        }

        
        #Scenarios
        Scens=List.sp[[i]]$Sens.test$JABBA_catch.only%>%
          mutate(Species=capitalize(names(dummy.store)[i]))
        Store.sens=vector('list',nrow(Scens))
        names(Store.sens)=Scens$Scenario
        this.wd1=this.wd
        
        Out.Scens=Scens%>%
          mutate(Bo.mean=NA,Bo.CV=NA,Bf.mean=NA,Bf.CV=NA,r.mean=NA,r.cv=NA,Rdist=NA,
                 Kdist=NA,PsiDist=NA,bmsyk.mean=NA)
        Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.probs.B.Bmsy=Out.f.series=
          Out.B.Bmsy=Out.F.Fmsy=Out.ensemble=vector('list',length(Store.sens))
        names(Out.ensemble)=names(Store.sens)
        
        CMSY.K.estim=mean(Catch_only$CMSY$estimates[[i]]%>%
                            filter(Parameter=='K')%>%
                            pull(Median)) 
        
        for(s in 1:length(Store.sens))
        {
          print(paste("___________","JABBA Scenario",Scens$Scenario[s],"___________",names(dummy.store)[i]))
          this.wd=paste(this.wd1,names(Store.sens)[s],sep='/')
          if(!dir.exists(this.wd))dir.create(this.wd)
          
          #Priors 
          bmsyk.mean=Scens$bmsyk[s]
          Proc.error=Scens$Proc.error[s]
          r.CV.multi=Scens$r.CV.multiplier[s]
          K.prior=c(Scens$K.mean[s],log(Scens$K.CV[s]))
          K.prior[1]=CMSY.K.estim*exp(rnorm(1,0,.1))
          if(names(dummy.store)[i]=="smooth hammerhead")  #bump up K a bit as not enough biomass if using CMSY estimate
          {
            K.prior[1]=mean(c(K.prior[1],CMSY.K.estim))
          }
          Mn=min(Scens$r[s],Max.r.value)  #some life history pars yield unrealistically high r
          Mn=max(Mn,Min.r.value)          #allow minimum recruitment
          r.prior=c(Mn, (Scens$r.sd[s]/Mn)*Scens$r.CV.multiplier[s])
          Ktch.CV=Scens$Ktch.CV[s]
          Bint=runif(1000,List.sp[[i]]$STARTBIO[1],List.sp[[i]]$STARTBIO[2])
          Bint.mean=mean(Bint)
          Bint.CV=sd(Bint)/Bint.mean
          psi.prior=c(Bint.mean,Bint.CV)
          #psi.prior=c(1,0.1) #Winker et al 2018 School shark
          
          Bfin=runif(1000,List.sp[[i]]$FINALBIO[1],List.sp[[i]]$FINALBIO[2])
          Bfin.mean=mean(Bfin)          
          Bfin.CV=sd(Bfin)/Bfin.mean
          b.prior=c(Bfin.mean,Bfin.CV,max(ktch$Year),"bk")
          
          Rdist = "lnorm"
          Kdist="lnorm"  
          PsiDist='beta'
          
          #Put inputs together
          input=list(Ktch=ktch,
                     MDL="Schaefer",
                     Ktch.CV=Ktch.CV,
                     ASS=names(dummy.store)[i],
                     Rdist = Rdist,
                     Rprior = r.prior,
                     r.CV.multiplier=r.CV.multi,
                     Kdist= Kdist,
                     Kprior=K.prior,    
                     PsiDist= PsiDist,
                     Psiprior=psi.prior,   
                     Bprior=b.prior,    
                     BMSYK=bmsyk.mean)  
          
          #run model
          THIN=5
          CHAINS=2
          BURNIN=min(0.15*Scens$Sims[s],5000)
          JABBA.run=apply.JABBA(Ktch=input$Ktch,
                                CPUE=NULL,
                                CPUE.SE=NULL,
                                MDL=input$MDL,
                                Ktch.CV=input$Ktch.CV,
                                ASS=input$ASS,
                                Rdist = input$Rdist,
                                Rprior = input$Rprior,
                                Kdist=input$Kdist,
                                Kprior=input$Kprior,
                                PsiDist=input$PsiDist,
                                Psiprior=input$Psiprior,
                                Bprior=input$Bprior,
                                BMSYK=input$BMSYK,
                                output.dir=this.wd,
                                outfile="Appendix_fit",
                                Sims=Scens$Sims[s],
                                Proc.error.JABBA=Proc.error,
                                thinning = THIN,
                                nchains = CHAINS,
                                burn.in= BURNIN,
                                Projections= add.ct.future)
          output=JABBA.run$fit
          output$posteriors$P=output$posteriors$P[,1:length(outktch)]
          output$posteriors$BtoBmsy=output$posteriors$BtoBmsy[,1:length(outktch)]
          output$posteriors$SB=output$posteriors$SB[,1:length(outktch)]
          
          
          #keep results if final depletion within specified priors
          keep.runs=which(output$posteriors$P[,nrow(ktch)]>=List.sp[[i]]$FINALBIO[1])
          if(length(keep.runs)>0)
          {
            #only select runs where depletion > prior's lower bound
            for(o in 1:length(output$posteriors))
            {
              if(length(dim(output$posteriors[[o]]))==2) output$posteriors[[o]]=output$posteriors[[o]][keep.runs,]
              if(length(dim(output$posteriors[[o]]))==3) output$posteriors[[o]]=output$posteriors[[o]][keep.runs,,]
            }
            output$kobe=output$kobe[keep.runs,]
            output$pars_posterior=output$pars_posterior[keep.runs,]
            
            Store.sens[[s]]=list(input=input,output=output)
            
            #Store Scenarios
            Out.Scens$Bo.mean[s]=Bint.mean
            Out.Scens$Bo.CV[s]=Bint.CV
            Out.Scens$Bf.mean[s]=Bfin.mean
            Out.Scens$Bf.CV[s]=Bfin.CV
            Out.Scens$r.mean[s]=Scens$r[s]
            Out.Scens$r.cv[s]=Scens$r.sd[s]/Scens$r[s]
            Out.Scens$Rdist[s]=Rdist
            Out.Scens$Kdist[s]=Kdist
            Out.Scens$PsiDist[s]=PsiDist
            Out.Scens$bmsyk.mean[s]=bmsyk.mean
            
            #Store estimates  
            d1=with(Store.sens[[s]]$output$posteriors,data.frame(K=K,r=r,psi=psi,
                                                                 Hmsy=Hmsy,SBmsy=SBmsy,MSY=MSY))%>%
              summarise_all(funs(list(quantile(., probs = c(0.25, 0.5, 0.75))))) %>%
              unnest(c(K, r, psi, Hmsy, SBmsy, MSY)) %>%
              transpose %>%
              setNames(., c("Lower.95","Median","Upper.95"))%>%
              mutate(Parameter=c('K', 'r', 'psi', 'Hmsy', 'SBmsy', 'MSY'),
                     Model='JABBA',
                     Scenario=names(Store.sens)[s])%>%
              relocate(Model,Scenario,Parameter,Lower.95) 
            Out.estimates[[s]]=d1
            
            #Store trajectories  
            dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                              mods=names(Catch_only)[w],
                                              Type='Depletion',
                                              scen=Scens$Scenario[s],
                                              Katch=outktch)
            Out.rel.biom[[s]]=dummy$Dat
            Out.probs.rel.biom[[s]]=dummy$Probs
            
            dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                              mods=names(Catch_only)[w],
                                              Type='F.series',
                                              scen=Scens$Scenario[s],
                                              Katch=outktch)
            Out.f.series[[s]]=dummy$Dat
            
            dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                              mods=names(Catch_only)[w],
                                              Type='B.Bmsy',
                                              scen=Scens$Scenario[s],
                                              Katch=outktch)
            Out.B.Bmsy[[s]]=dummy$Dat
            Out.probs.B.Bmsy[[s]]=dummy$Probs
            
            dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                              mods=names(Catch_only)[w],
                                              Type='F.Fmsy',
                                              scen=Scens$Scenario[s],
                                              Katch=outktch)
            Out.F.Fmsy[[s]]=dummy$Dat
            rm(dummy)
            
            if(Scens$Scenario[s]=='S1') Out.Kobe.probs=Store.sens[[s]]$output$kobe%>%dplyr::select(stock,harvest)
            
            #Evaluate posterior  
            #https://cran.r-project.org/web/packages/ggmcmc/vignettes/using_ggmcmc.html
            #MCMC trace-plots, Geweke (2003) and Gelman and Rubin (1992) diagnostics 
            if(Scens$Scenario[s]=='S1')
            {
              these.post.pars=c('K','r','psi')  #estimable pars
              NN=nrow(output$pars_posterior)/CHAINS
              Post=Store.sens[[s]]$output$pars_posterior%>%
                dplyr::select(all_of(these.post.pars))%>%
                mutate(Iteration=rep(1:NN,CHAINS),
                       Chain=rep(1:CHAINS,each=NN))%>%
                gather(Parameter,value,-Iteration,- Chain)%>%
                tibble()
              attributes(Post)$nChains=CHAINS
              attributes(Post)$nParameters=length(these.post.pars)
              attributes(Post)$nIterations=NN
              attributes(Post)$nBurnin=BURNIN
              attributes(Post)$nThin=THIN
              
              Post.density=ggs_density(Post)+ggtitle("Density")
              Post.trace.plot=ggs_traceplot(Post)+ggtitle("Trace plot")+theme(axis.text.x = element_text(angle = 90))
              Post.running.mean=ggs_running(Post)+ggtitle('Running mean')+theme(axis.text.x = element_text(angle = 90))
              Post.autocorr=ggs_autocorrelation(Post)+ggtitle('Autocorrelation')
              Gelman.diag=ggs_Rhat(Post) + 
                geom_point(size=3)+
                xlab("R_hat")+
                ggtitle('Gelman diagnostic')
              Geweke.diag= ggs_geweke(Post)+
                geom_point(size=3)+
                ggtitle('Geweke diagnostic')+
                theme(legend.position = 'none')+
                scale_color_manual(values=c("#F8766D", "#00BFC4", "#56B4E9"))
              
              pp=ggarrange(plotlist = list(Gelman.diag,Geweke.diag),ncol=1)
              ggarrange(plotlist = list(Post.density,Post.trace.plot,Post.running.mean,
                                        Post.autocorr,pp),
                        ncol=5,common.legend = TRUE)
              
              ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                           capitalize(names(dummy.store)[i]),"/",AssessYr,"/",'JABBA',
                           "/Convergence diagnostics.tiff",sep=''),
                     width = 16,height = 8, dpi = 300, compression = "lzw")
              
              rm(Post,Post.density,Post.trace.plot,Post.running.mean,
                 Post.autocorr,Gelman.diag,Geweke.diag)
            }
            
            #Display Priors vs Posteriors for base case scenario (S1)
            if(Scens$Scenario[s]=='S1')
            {
              par.list=c("r","K","psi")
              stuff=Store.sens[[s]]$input
              dummy.prior=as.list(par.list)
              names(dummy.prior)=par.list
              dummy.prior$r=list(dist=stuff$Rdist,
                                 mean=log(stuff$Rprior[1]),
                                 sd=stuff$Rprior[2])
              dummy.prior$K=list(dist='lnorm',
                                 mean=log(Store.sens[[s]]$output$settings$K.pr[1]),  
                                 sd=Store.sens[[s]]$output$settings$K.pr[2])
              dummy.prior$psi=list(dist=stuff$PsiDist,
                                   mean=stuff$Psiprior[1],
                                   sd=stuff$Psiprior[2])
              dummy.post=Store.sens[[s]]$output$pars_posterior%>%
                dplyr::select(r,K,psi)
              
              names(dummy.prior)=names(dummy.post)
              
              pars=colnames(dummy.post)
              out=vector('list',length(pars))
              
              
              for(p in 1:length(pars))
              {
                Value.prior=fn.prior(d=dummy.prior[[pars[p]]])
                if(pars[p]=="K") Value.prior=fn.prior(d=dummy.prior[[pars[p]]],MAX=k.fun(stuff$Ktch$Total,25))
                out[[p]]=rbind(data.frame(Distribuion="Prior",
                                          Value=Value.prior),
                               data.frame(Distribuion="Posterior",
                                          Value=dummy.post[[pars[p]]]))%>%
                  mutate(Parameter=pars[p])
              }
              rm(stuff,dummy.prior,dummy.post)
              
              fn.show.density(d=do.call(rbind,out),NCOL=1)
              ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                           capitalize(names(dummy.store)[i]),"/",AssessYr,"/",'JABBA',"/Prior.and.posterior.tiff",sep=''),
                     width = 12,height = 14, dpi = 300, compression = "lzw")
            }
            
            #Store quantities for ensemble model
            Out.ensemble[[s]]=list(Depletion=Store.sens[[s]]$output$posteriors$P,
                                   B.Bmsy=Store.sens[[s]]$output$posteriors$BtoBmsy,
                                   F.Fmsy=Store.sens[[s]]$output$posteriors$HtoHmsy,
                                   MSY=data.frame(Par='MSY',
                                                  Value=Store.sens[[s]]$output$posteriors$MSY,
                                                  Scenario=Scens$Scenario[s]))  
          }
          
        }
        #output quantities for ensemble
        Out.ens=list(Depletion=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'Depletion')),
                     B.Bmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'B.Bmsy')),
                     F.Fmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'F.Fmsy')),
                     MSY=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'MSY')))
        if(nrow(Scens)>1)
        {
          for(f in 1:length(Out.ens))   
          {
            ddmi=Out.ens[[f]]
            TAB=table(sub("\\..*", "", rownames(ddmi)))
            Out.ens[[f]]=ddmi[sample(1:nrow(ddmi),size=max(TAB),replace =FALSE),]
          }
        }
        
        #output scenarios
        #dummy.store[[i]]=Store.sens
        Out.Scens=Out.Scens%>%
          dplyr::select(Species,Scenario,
                        Rdist,r,r.sd,
                        Kdist,K.mean,K.CV,
                        PsiDist,Bo.mean,Bo.CV,
                        Bf.mean,Bf.CV,
                        bmsyk.mean,Proc.error,Ktch.CV)%>%
          mutate(K.mean=round(K.mean),
                 r=round(r,3),
                 bmsyk.mean=round(bmsyk.mean,3),
                 r.sd=round(r.sd,3),
                 Bo.mean=round(Bo.mean,3),
                 Bo.CV=round(Bo.CV,3),
                 Bf.mean=round(Bf.mean,3),
                 Bf.CV=round(Bf.CV,3),
                 Rdist=ifelse(Rdist=='lnorm','Lognormal',Rdist),
                 Kdist=ifelse(Kdist=='lnorm','Lognormal',Kdist),
                 PsiDist=ifelse(PsiDist=='lnorm','Lognormal',PsiDist))
        
        dummy.store.sens.table[[i]]=Out.Scens
        dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
        dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
        dummy.store.probs.B.Bmsy[[i]]=Out.probs.B.Bmsy
        dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
        dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
        dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
        dummy.store.Kobe.probs[[i]]=Out.Kobe.probs  
        dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
        dummy.store.ensemble[[i]]=Out.ens
      }
    } #end i
    
    #Catch_only[[w]]=dummy.store       #too big an object, cannot store whole model runs due to memory size limitations
    Catch_only[[w]]$sens.table=dummy.store.sens.table
    Catch_only[[w]]$estimates=dummy.store.estimates
    Catch_only[[w]]$rel.biom=dummy.store.rel.biom
    Catch_only[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
    Catch_only[[w]]$probs.B.Bmsy=dummy.store.probs.B.Bmsy
    Catch_only[[w]]$f.series=dummy.store.f.series
    Catch_only[[w]]$B.Bmsy=dummy.store.B.Bmsy
    Catch_only[[w]]$F.Fmsy=dummy.store.F.Fmsy
    Catch_only[[w]]$Kobe.probs=dummy.store.Kobe.probs
    Catch_only[[w]]$ensemble=dummy.store.ensemble
    
    rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,
       dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,
       dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.ensemble,
       dummy.store.probs.B.Bmsy)
  }
  
  #5. SSS-MC assessment (Cope (2013))   
  # summary of method: https://github.com/shcaba/SSS
  #note: this was parallelised to improve computation time
  if(names(Catch_only)[w]=="SSS")     #parallelised: 3.2 secs per iteration-species-scenario (otherwise 12 secs)
  {
    if(do.parallel.SSS)
    {
      set.seed(1234)
      progress <- function(n) cat(Keep.species[n],sprintf(": SSS fit complete-----------\n", n))
      opts <- list(progress = progress)
      cl <- makeCluster(detectCores()-1)
      registerDoSNOW(cl)
      out.species=foreach(i= 1:N.sp,.options.snow = opts,.errorhandling="pass",
                          .packages=c('gridExtra','Hmisc','JABBA','TruncatedDistributions','tidyverse','r4ss','zoo')) %dopar%
      {
        Neim=names(List.sp)[i]
        if(Neim%in%Catch.only.species) 
        {
          this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                        capitalize(Neim),"/",AssessYr,"/SS3",sep='')
          if(!dir.exists(this.wd))dir.create(this.wd)
          
          #Criteria for accepting run
          Criteria.delta.fin.dep=SSS_criteria.delta.fin.dep  
          if(Neim%in%c('grey nurse shark',"pigeye shark","shortfin mako"))Criteria.delta.fin.dep=0.1
          if(Neim=='milk shark')Criteria.delta.fin.dep=0.2
          
          #Catch
          ktch=ktch.combined%>%
            filter(Name==Neim)
          all.years=seq(min(ktch$finyear),max(ktch$finyear))
          misn.yr=all.years[which(!all.years%in%ktch$finyear)]
          if(length(misn.yr)>0)
          {
            ktch=rbind(ktch,ktch[length(misn.yr),]%>%mutate(finyear=misn.yr,Tonnes=0))%>%arrange(finyear)
          }
          
          #future catches
          outktch=ktch$Tonnes
          add.ct.future=NULL
          if('Catch_only'%in%future.models)
          {
            if(catches.futures=="constant.last.n.yrs")
            {
              NN=nrow(ktch)
              add.ct.future=ktch[1:years.futures,]%>%
                mutate(finyear=ktch$finyear[NN]+(1:years.futures),
                       Tonnes=rep(mean(ktch$Tonnes[(NN-(n.last.catch.yrs-1)):NN]),years.futures))
              outktch=c(outktch,add.ct.future$Tonnes)
              #ktch=rbind(ktch,add.ct.future)
            }
          }

          
          #random LH samples
          LH.sim=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                                capitalize(List.sp[[i]]$Name),"/",AssessYr,"/steepness/Life.history_M.mean.csv",sep=''))
          
          #Scenario
          Scens=List.sp[[i]]$Sens.test$SS3%>%
            mutate(Species=capitalize(Neim))
          Out.Scens=Scens%>%
            mutate(h.mean=NA,h.sd=NA,Initial.bio.low=NA,Initial.bio.up=NA,Final.bio.low=NA,Final.bio.up=NA)
          Store.sens=vector('list',nrow(Scens))
          names(Store.sens)=Scens$Scenario
          Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.probs.B.Bmsy=Out.f.series=Out.B.Bmsy=
            Out.F.Fmsy=store.warnings=store.convergence=Out.ensemble=
            Out.accept.rate=Out.prior.post=vector('list',length(Store.sens))
          names(Out.ensemble)=names(Store.sens)
          
          for(s in 1:length(Store.sens))
          {
            this.wd1=paste(this.wd,Scens$Scenario[s],sep='/')
            if(!dir.exists(this.wd1))dir.create(this.wd1)
            
            n.sims_SS=Scens$Sims[s]
            Report=vector('list',n.sims_SS)
            h.sim=dep.f=m.sim=rep(NA,n.sims_SS)
            
            #get h beta distribution pars
            # h.trunc.normal=rtruncnorm(1e3,a=Min.h.shark,b=Max.h.shark,
            #                           mean=Out.Scens$Steepness[s],sd=max(Out.Scens$Steepness.sd[s],0.01))
            # h.beta=get_beta(mean(h.trunc.normal),sd(h.trunc.normal)/h.trunc.normal)
            h.beta=get_beta(Out.Scens$Steepness[s],max(Out.Scens$Steepness.sd[s],0.01))
            
            #get M lognormal distribution pars
            m.mean=Scens$Mmean[s]
            m.sd=Scens$Msd[s]
            
            #loop over simulations only accepting good runs
            n=trials=1
            while(n <=n.sims_SS)  
            {
              
              #1. Get random values of interest
              if(do.random.h)
              {
                h.sim[n]=rtbeta(1,h.beta[1],h.beta[2],a=Min.h.shark,b=Max.h.shark)
                m.sim[n]=rlnorm(1,mean=log(m.mean),sd=m.sd/m.mean)
              }else
              {
                h.sim[n]=LH.sim$h[n]
                if(Keep.species[i]=="grey nurse shark" & resample.h.greynurse)  #too low h for greynurse
                {
                  h.sim[n]=rtruncnorm(1,a=Min.h.shark,b=Max.h.shark, 
                                      mean=Out.Scens$Steepness[s],
                                      sd=max(Out.Scens$Steepness.sd[s],0.01))
                }
                m.sim[n]=LH.sim$M[n] 
              }
              
              Depletion.year=ktch$finyear[nrow(ktch)]
              Btklow=List.sp[[i]]$FINALBIO[1]
              # if(Neim=="milk shark")
              # {
              #   Depletion.year=1991  #set to end of catch period to allow convergence
              #   Btklow=0.2
              # }
              
              dep.f[n]=runif(1,Btklow,List.sp[[i]]$FINALBIO[2])
              Scenario=Scens[s,]%>%
                mutate(Mmean=m.sim[n],
                       Steepness=h.sim[n],
                       Initial.dpl=runif(1,List.sp[[i]]$STARTBIO[1],List.sp[[i]]$STARTBIO[2]),
                       Final.dpl=dep.f[n],
                       Ln_R0_init=runif(1,Ln_R0_min,Ln_R0_max))
              Life.history=List.sp[[i]]
              if(do.random.h)
              {
                Life.history$Fecundity=mean(Life.history$Fecundity)
                Life.history$Max.age.F=Life.history$Max.age.F[1]
                Life.history$Breed.cycle=mean(Life.history$Breed.cycle)
              }else
              {
                if(!Keep.species[i]=="grey nurse shark")    #no convergence with resampled pars, misspecification in life history
                {
                  Life.history$Fecundity=LH.sim$Fecundity[n]
                  Life.history$Max.age.F=LH.sim$Max.age[n]
                  Life.history$Breed.cycle=LH.sim$Rep.cycle[n]
                  Life.history$Growth.F[1:2]=c(LH.sim$k[n],LH.sim$Linf[n])
                }
              }
              
              
              #2. Create input files
              fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SSS'),
                           new.path=this.wd1,
                           Scenario=Scenario,
                           Catch=ktch,
                           life.history=Life.history,
                           depletion.yr=Depletion.year,
                           Future.project=add.ct.future)
              
              #3. Run SS3
              fn.run.SS(where.inputs=this.wd1,
                        where.exe=handl_OneDrive('SS3/ss_win.exe'),
                        args='-nohess')  #no Hessian estimation as uncertainty generated thru MC procedure
              
              #4. Bring in outputs
              Report[[n]]=SS_output(this.wd1,covar=F,forecast=F,readwt=F,checkcor=F)
              #ss.std=read.table(paste(this.wd1,'ss.std',sep='/')) https://vlab.noaa.gov/web/stock-synthesis/public-forums/-/message_boards/view_message/11664630
              Report[[n]]$Final.dpl=data.frame(Distribuion=c('Prior','Posterior'),
                                               Value=c(Scenario$Final.dpl,Report[[n]]$current_depletion),
                                               Parameter=rep('Final depletion',2))
              Report[[n]]$Initial.dpl=data.frame(Distribuion=c('Prior','Posterior'),
                                                 Value=c(Scenario$Initial.dpl,Report[[n]]$derived_quants[grep('Bratio',Report[[n]]$derived_quants$Label)[1],'Value']),
                                                 Parameter=rep('Initial depletion',2))
              Report[[n]]$Ln.R0=data.frame(Distribuion=c('Prior','Posterior'),
                                           Value=with(Report[[n]]$estimated_non_dev_parameters,c(Init,Value)),
                                           Parameter=rep('Ln.R0',2))
              
              #5. Keep only runs that converged
              Converged=data.frame(Max.grad=Report[[n]]$maximum_gradient_component,
                                   dep.f_estim=Report[[n]]$current_depletion,
                                   dep.f_obs=dep.f[[n]],
                                   LnRo=Report[[n]]$estimated_non_dev_parameters$Value,
                                   LnRo.min=Report[[n]]$estimated_non_dev_parameters$Min,
                                   LnRo.max=Report[[n]]$estimated_non_dev_parameters$Max)%>%
                mutate(Delta.conv.criteria=abs(Max.grad - 1e-4),
                       Delta.fin.dep=abs(dep.f_obs-dep.f_estim),
                       Delta.lower.LnRo=LnRo/LnRo.min,
                       Delta.upper.LnRo=LnRo/LnRo.max,
                       Keep=ifelse(Delta.fin.dep>Criteria.delta.fin.dep |
                                     Delta.upper.LnRo>0.9 | Delta.lower.LnRo<5,'No','Yes'))
              #     Keep=ifelse(Delta.conv.criteria>0.001 | Delta.fin.dep>Criteria.delta.fin.dep |
              #          Delta.upper.LnRo<0.5 | Delta.lower.LnRo>0.5,'No','Yes'))
              trials=trials+1
              if(Converged$Keep=='Yes')
              {
                Add.bmsy=Report[[n]]$derived_quants%>%filter(Label%in%c('SSB_Virgin','SSB_MSY'))
                Add.bmsy.val=Add.bmsy[Add.bmsy$Label=='SSB_MSY','Value']/Add.bmsy[Add.bmsy$Label=='SSB_Virgin','Value']
                Add.bmsy=Report[[n]]$estimated_non_dev_parameters[1,]
                rownames(Add.bmsy)='bmsy/k'
                Add.bmsy[,]=NA
                Add.bmsy=Add.bmsy%>%
                  mutate(Value=Add.bmsy.val)
                Report[[n]]$estimated_non_dev_parameters=rbind(Report[[n]]$estimated_non_dev_parameters,Add.bmsy)
                
                n=n+1
              }
            }
            
            #Store estimates
            if(length(Report)>0)
            {
              #Store acceptance rate
              Out.accept.rate[[s]]=data.frame(Acceptance=100*n.sims_SS/trials, 
                                              Scenario=Scens$Scenario[s])
              
              #Store estimates (recalculate with retained runs)
              Estims=do.call(rbind,fn.get.stuff.from.list(Report,"estimated_non_dev_parameters"))
              Estims=Estims%>%
                mutate(Par=gsub('[0-9]+', '', rownames(Estims)),
                       Scenario=Scens$Scenario[s])
              Out.estimates[[s]]=Estims%>%
                group_by(Par,Scenario)%>%
                summarise(Median=median(Value),
                          Lower.95=quantile(Value,probs=0.025),
                          Upper.95=quantile(Value,probs=0.975))%>%
                dplyr::select(Par,Median,Lower.95,Upper.95,Scenario)%>%
                ungroup()
              add.MSY=do.call(rbind,fn.get.stuff.from.list(Report,"derived_quants"))%>%
                filter(grepl('Dead_Catch_MSY',Label))%>%
                rename(Par=Label)%>%
                mutate(Scenario=Scens$Scenario[s])
              MSY.ensemble=add.MSY%>%dplyr::select(Par,Value,Scenario)%>%mutate(Par='MSY')
              add.MSY=add.MSY%>%
                group_by(Par,Scenario)%>%
                summarise(Median=median(Value),
                          Lower.95=quantile(Value,probs=0.025),
                          Upper.95=quantile(Value,probs=0.975))%>%
                dplyr::select(Par,Median,Lower.95,Upper.95,Scenario)%>%
                ungroup()
              Out.estimates[[s]]=rbind(Out.estimates[[s]],add.MSY)
              
              #Store Scenarios
              Out.Scens$h.mean[s]=Out.Scens$Steepness[s]
              Out.Scens$h.sd[s]=max(Out.Scens$Steepness.sd[s],0.01)
              Out.Scens$Initial.bio.low[s]=List.sp[[i]]$STARTBIO[1]
              Out.Scens$Initial.bio.up[s]=List.sp[[i]]$STARTBIO[2]
              Out.Scens$Final.bio.low[s]=List.sp[[i]]$FINALBIO[1]
              Out.Scens$Final.bio.up[s]=List.sp[[i]]$FINALBIO[2]
              
              #Store priors and posteriors
              out.pri.post=list(Final.dpl=do.call(rbind,fn.get.stuff.from.list(Report,'Final.dpl')),
                                Initial.dpl=do.call(rbind,fn.get.stuff.from.list(Report,'Initial.dpl')),
                                Ln.R0=do.call(rbind,fn.get.stuff.from.list(Report,'Ln.R0'))%>%
                                      mutate(Parameter='Ln.R0')) 
              out.pri.post$Final.dpl=out.pri.post$Final.dpl%>%
                            mutate(Parameter='Final.dpl',
                                   Value=ifelse(Distribuion=='Prior',runif(nrow(out.pri.post$Final.dpl),List.sp[[i]]$FINALBIO[1],List.sp[[i]]$FINALBIO[2]),
                                                Value))
              out.pri.post$Initial.dpl=out.pri.post$Initial.dpl%>%
                        mutate(Parameter='Initial.dpl',
                               Value=ifelse(Distribuion=='Prior',runif(nrow(out.pri.post$Initial.dpl),List.sp[[i]]$STARTBIO[1],List.sp[[i]]$STARTBIO[2]),
                                            Value))
              Out.prior.post[[s]]=do.call(rbind,out.pri.post)   
                
              #Store trajectories
              dummy=fn.ktch.only.get.timeseries(d=Report,
                                                mods=names(Catch_only)[w],
                                                Type='Depletion',
                                                scen=Scens$Scenario[s],
                                                Katch=outktch)
              Out.rel.biom[[s]]=dummy$Dat
              Out.probs.rel.biom[[s]]=dummy$Probs
              
              dummy=fn.ktch.only.get.timeseries(d=Report,
                                                mods=names(Catch_only)[w],
                                                Type='F.series',
                                                scen=Scens$Scenario[s],
                                                Katch=outktch)
              Out.f.series[[s]]=dummy$Dat
              
              dummy=fn.ktch.only.get.timeseries(d=Report,
                                                mods=names(Catch_only)[w],
                                                Type='B.Bmsy',
                                                scen=Scens$Scenario[s],
                                                Katch=outktch)
              Out.B.Bmsy[[s]]=dummy$Dat
              Out.probs.B.Bmsy[[s]]=dummy$Probs
              
              dummy=fn.ktch.only.get.timeseries(d=Report,
                                                mods=names(Catch_only)[w],
                                                Type='F.Fmsy',
                                                scen=Scens$Scenario[s],
                                                Katch=outktch)
              Out.F.Fmsy[[s]]=dummy$Dat
              
              
              #Store quantities for ensemble model
              Years=c(ktch$finyear,add.ct.future$finyear)
              
              dum=do.call(rbind,fn.get.stuff.from.list(Report,"derived_quants"))
              
              SSB_Virgin=dum[grep("SSB_Virgin",dum$Label),]
              SSB_Virgin=SSB_Virgin%>%
                mutate(Simulation=readr::parse_number(paste(rownames(SSB_Virgin),'0',sep='')))%>%
                rename(Value.virgin=Value)%>%
                dplyr::select(Value.virgin,Simulation)
              
              SSB_MSY=dum[grep("SSB_MSY",dum$Label),]
              SSB_MSY=SSB_MSY%>%
                mutate(Simulation=readr::parse_number(paste(rownames(SSB_MSY),'0',sep='')))%>%
                rename(Value.MSY=Value)%>%
                dplyr::select(Value.MSY,Simulation)
              
              annF_MSY=dum[grep("annF_MSY",dum$Label),]
              annF_MSY=annF_MSY%>%
                mutate(Simulation=readr::parse_number(paste(rownames(annF_MSY),'0',sep='')))%>%
                rename(Value.F_MSY=Value)%>%
                dplyr::select(Value.F_MSY,Simulation)  
              
              
              SSB=dum[grep(paste(paste("SSB",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
              SSB=SSB%>%
                mutate(year=readr::parse_number(Label),
                       Simulation=readr::parse_number(paste(rownames(SSB),'0',sep='')),
                       Simulation=as.numeric(substr(Simulation,5,10)))%>%
                dplyr::select(-Label)
              
              EF=dum[grep(paste(paste("F",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
              EF=EF%>%
                mutate(year=readr::parse_number(Label),
                       Simulation=readr::parse_number(paste(rownames(EF),'0',sep='')),
                       Simulation=as.numeric(substr(Simulation,5,10)))%>%
                dplyr::select(-Label)
              Misn.yr=which(!paste(SSB$year,SSB$Simulation)%in%paste(EF$year,EF$Simulation)) 
              All.misn.yrs=SSB$year[Misn.yr]
              Misn.yr_no.fishing=All.misn.yrs[which(All.misn.yrs<min(EF$year))]
              if(length(Misn.yr_no.fishing)>0)
              {
                add.dis=EF[1:length(Misn.yr_no.fishing),]%>%
                  mutate(Value=0,
                         Simulation=rep(sort(unique(SSB$Simulation)),each=length(unique(Misn.yr_no.fishing))),
                         year=Misn.yr_no.fishing)
                kkk=str_sub(add.dis$Simulation,1,-2)
                rownames(add.dis)=paste('F_',add.dis$year,kkk,sep='')
                EF=rbind(EF,add.dis)%>%
                  arrange(Simulation,year)
                
                Misn.yr=which(!paste(SSB$year,SSB$Simulation)%in%paste(EF$year,EF$Simulation))
                if(length(Misn.yr)>0)
                {
                  add.dis=EF[Misn.yr,]%>%
                    mutate(Value=NA,
                           Simulation=sort(unique(SSB$Simulation)),
                           year=SSB$year[Misn.yr])
                  rownames(add.dis)=paste('F_',SSB$year[Misn.yr],c('',1:(nrow(add.dis)-1)),sep='')
                  EF=rbind(EF,add.dis)%>%
                    arrange(Simulation,year)%>%
                    mutate(Value = na.approx(Value))
                }
              }
              
              Depletion=SSB%>%
                left_join(SSB_Virgin%>%dplyr::select(Value.virgin,Simulation),by='Simulation')%>%
                mutate(Depletion=Value/Value.virgin)%>%
                dplyr::select(Simulation,year,Depletion)%>%
                spread(year,Depletion)%>%
                dplyr::select(-Simulation)
              
              B.Bmsy=SSB%>%
                left_join(SSB_MSY%>%dplyr::select(Value.MSY,Simulation),by='Simulation')%>%
                mutate(B_Bmys=Value/Value.MSY)%>%
                dplyr::select(Simulation,year,B_Bmys)%>%
                spread(year,B_Bmys)%>%
                dplyr::select(-Simulation)
              
              F.Fmsy=EF%>%
                left_join(annF_MSY%>%dplyr::select(Value.F_MSY,Simulation),by='Simulation')%>%
                mutate(F_Fmys=Value/Value.F_MSY)%>%
                dplyr::select(Simulation,year,F_Fmys)%>%
                spread(year,F_Fmys)%>%
                dplyr::select(-Simulation)
              Out.ensemble[[s]]=list(Depletion=Depletion,  
                                     B.Bmsy=B.Bmsy,
                                     F.Fmsy=F.Fmsy,
                                     MSY=MSY.ensemble) 
              rm(Depletion,B.Bmsy,F.Fmsy,dum,MSY.ensemble)
              
              
              #Flag warnings
              n.sims_SS=length(Report)
              warnings=fn.get.stuff.from.list(Report,"warnings")
              dummy.warnings=vector('list',n.sims_SS)
              for(n in 1:n.sims_SS)
              {
                dummy.warnings[[n]]=data.frame(warning=warnings[[n]])%>%
                  mutate(warning.type=case_when(grepl('bad Beta prior',warning)~'bad Beta prior',
                                                grepl('Final gradient',warning)~'Final gradient',
                                                grepl('midsize bin',warning)~'midsize bin',
                                                grepl('Fmsy/mey is <',warning)~'Fmsy/mey',
                                                grepl('Maximum pop size bin',warning)~'Maximum pop size bin'))%>%
                  distinct(warning.type,.keep_all = T)%>%
                  filter(!is.na(warning.type))%>%
                  mutate(Steepness=h.sim[n],
                         Final.depl=dep.f[n])
                
                
              }
              store.warnings[[s]]=do.call(rbind,dummy.warnings)
              
              TAB=data.frame(Bad.beta.prior=rep(NA,n.sims_SS))%>%
                mutate(Final.grad=Bad.beta.prior,
                       midsize.bin=Bad.beta.prior)
              for(n in 1:n.sims_SS)
              {
                TAB$Bad.beta.prior[n]=any(grepl('bad Beta prior',Report[[n]]$warnings))
                TAB$Final.grad[n]=any(grepl('Final gradient',Report[[n]]$warnings))
                TAB$midsize.bin[n]=any(grepl('midsize bin',Report[[n]]$warnings))
              }
              
            }
            
            rm(this.wd1,Report)
          } # end s loop
          #output quantities for ensemble
          Out.ens=list(Depletion=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'Depletion')),
                       B.Bmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'B.Bmsy')),
                       F.Fmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'F.Fmsy')),
                       MSY=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'MSY')))
          if(nrow(Scens)>1)
          {
            for(f in 1:length(Out.ens))   
            {
              ddmi=Out.ens[[f]]
              TAB=table(sub("\\..*", "", rownames(ddmi)))
              Out.ens[[f]]=ddmi[sample(1:nrow(ddmi),size=max(TAB),replace =FALSE),]
            }
          }
          
          #output scenarios
          Out.Scens=Out.Scens%>%
            mutate(Mmean=round(Mmean,3),
                   h.mean=round(h.mean,3),
                   h.sd=round(h.sd,3))%>%
            dplyr::select(-c(Model,Final.dpl, Steepness,Steepness.sd,use_F_ballpark, Sims))%>%
            relocate(Species,Scenario,Mmean,F_ballpark,h.mean,h.sd,
                     Initial.bio.low,Initial.bio.up,Final.bio.low,Final.bio.up)
          
          
          return(list(dummy.store.sens.table=Out.Scens,dummy.store.estimates=do.call(rbind,Out.estimates),
                      dummy.store.rel.biom=do.call(rbind,Out.rel.biom), dummy.store.probs.rel.biom=Out.probs.rel.biom,
                      dummy.store.probs.B.Bmsy=Out.probs.B.Bmsy, 
                      dummy.store.f.series=do.call(rbind,Out.f.series),dummy.store.B.Bmsy=do.call(rbind,Out.B.Bmsy),
                      dummy.store.F.Fmsy=do.call(rbind,Out.F.Fmsy),dummy.store.Kobe.probs=NULL,
                      dummy.store.ensemble=Out.ens,dummy.store.warnings=do.call(rbind,store.warnings),
                      dummy.store.convergence=do.call(rbind,store.convergence),
                      dummy.store.accept.rate=do.call(rbind,Out.accept.rate),
                      dummy.store.out.prior.post=Out.prior.post))   
        }
      }  # end i parallel loop
      stopCluster(cl)
      names(out.species)=Keep.species
      
      Catch_only[[w]]$sens.table=fn.get.and.name(LISTA=out.species,x="dummy.store.sens.table")
      Catch_only[[w]]$estimates=fn.get.and.name(LISTA=out.species,x="dummy.store.estimates")
      Catch_only[[w]]$rel.biom=fn.get.and.name(LISTA=out.species,x="dummy.store.rel.biom")
      Catch_only[[w]]$probs.rel.biom=fn.get.and.name(LISTA=out.species,x="dummy.store.probs.rel.biom")
      Catch_only[[w]]$probs.B.Bmsy=fn.get.and.name(LISTA=out.species,x="dummy.store.probs.B.Bmsy")
      Catch_only[[w]]$f.series=fn.get.and.name(LISTA=out.species,x="dummy.store.f.series")
      Catch_only[[w]]$B.Bmsy=fn.get.and.name(LISTA=out.species,x="dummy.store.B.Bmsy")
      Catch_only[[w]]$F.Fmsy=fn.get.and.name(LISTA=out.species,x="dummy.store.F.Fmsy")
      Catch_only[[w]]$Kobe.probs=fn.get.and.name(LISTA=out.species,x="dummy.store.Kobe.probs")
      Catch_only[[w]]$ensemble=fn.get.and.name(LISTA=out.species,x="dummy.store.ensemble")
      Catch_only[[w]]$warnings=fn.get.and.name(LISTA=out.species,x="dummy.store.warnings")
      Catch_only[[w]]$convergence=fn.get.and.name(LISTA=out.species,x="dummy.store.convergence")
      Catch_only[[w]]$accept.rate=fn.get.and.name(LISTA=out.species,x="dummy.store.accept.rate")
      Catch_only[[w]]$out.prior.post=fn.get.and.name(LISTA=out.species,x="dummy.store.out.prior.post")
      
      #plot overall figure with estimates and scenarios
      for(i in 1:N.sp)
      {
        Neim=names(List.sp)[i]
        if(Neim%in%Catch.only.species) 
        {
          print(paste("SSS rel biom figure and estimates ____",Keep.species[i]))
          yrS=Catch_only[[w]]$rel.biom[[i]]%>%filter(Scenario=='S1')%>%pull(year)
          if(length(yrS)<50) delta=10 else
            delta=17
          xmin=min(yrS)+delta
          xmax=xmin+delta
          p=Catch_only[[w]]$rel.biom[[i]]%>%
            ggplot(aes(year,median,color=Scenario))+
            annotation_custom(tableGrob(Catch_only[[w]]$sens.table[[i]]%>%
                                          dplyr::select(Scenario,Mmean,h.mean,Final.bio.low,Final.bio.up)),
                              xmin=xmin+5, xmax=xmax+5, ymin=0, ymax=0.3)+
            annotation_custom(tableGrob(Catch_only[[w]]$estimates[[i]]),xmin=xmin, xmax=xmax, ymin=0.45, ymax=0.65)+
            geom_line(size=2)+
            geom_line(aes(year,upper.95),linetype=2)+
            geom_line(aes(year,lower.95),linetype=2)+
            ggtitle(Keep.species[i])+ylim(0,1)+
            theme_PA()+theme(legend.position = 'bottom')
          print(p)
          ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                       capitalize(List.sp[[i]]$Name),"/",AssessYr,"/SS3/Rel.biomass&estimates.tiff",sep=''),compression = "lzw")
          
        }
      }
      #Plot priors and posterior
      for(i in 1:N.sp)
      {
        Neim=names(List.sp)[i]
        if(Neim%in%Catch.only.species) 
        {
          print(paste("SSS priors and posteriors ____",Keep.species[i]))
          fn.show.density(d=do.call(rbind,Catch_only[[w]]$out.prior.post[[i]]),NCOL=2) 
          ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[i]]$Name),"/",AssessYr,"/SS3/Prior.and.posterior_",Scens$Scenario[s],".tiff",sep=''),compression = "lzw")
          
        }
      }

      rm(out.species)
    }
    if(!do.parallel.SSS)
    {
      dummy.store=vector('list',N.sp)
      names(dummy.store)=Keep.species
      dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
        dummy.store.probs.B.Bmsy=dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.Kobe.probs=
        dummy.store.sens.table=dummy.store.ensemble=dummy.store.warnings=dummy.store.accept.rate=
        dummy.store.convergence=dummy.store
      for(i in 1:N.sp)
      {
        Neim=names(List.sp)[i]
        if(Neim%in%Catch.only.species)
        {
          this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                        capitalize(Neim),"/",AssessYr,"/SS3",sep='')
          if(!dir.exists(this.wd))dir.create(this.wd)
          
          #Criteria for accepting run
          Criteria.delta.fin.dep=SSS_criteria.delta.fin.dep  
          if(Neim%in%c('grey nurse shark',"pigeye shark","shortfin mako"))Criteria.delta.fin.dep=0.1
          if(Neim=='milk shark')Criteria.delta.fin.dep=0.2
          
          #Catch
          ktch=ktch.combined%>%
            filter(Name==Neim)
          all.years=seq(min(ktch$finyear),max(ktch$finyear))
          misn.yr=all.years[which(!all.years%in%ktch$finyear)]
          if(length(misn.yr)>0)
          {
            ktch=rbind(ktch,ktch[length(misn.yr),]%>%mutate(finyear=misn.yr,Tonnes=0))%>%arrange(finyear)
          }
          #future catches
          outktch=ktch$Tonnes
          add.ct.future=NULL
          if('Catch_only'%in%future.models)
          {
            if(catches.futures=="constant.last.n.yrs")
            {
              NN=nrow(ktch)
              add.ct.future=ktch[1:years.futures,]%>%
                mutate(finyear=ktch$finyear[NN]+(1:years.futures),
                       Tonnes=rep(mean(ktch$Tonnes[(NN-(n.last.catch.yrs-1)):NN]),years.futures))
              outktch=c(outktch,add.ct.future$Tonnes)
              #ktch=rbind(ktch,add.ct.future)
            }
          }

          
          
          #random LH samples
          LH.sim=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                                capitalize(List.sp[[i]]$Name),"/",AssessYr,"/steepness/Life.history_M.mean.csv",sep=''))
          
          #Scenario
          Scens=List.sp[[i]]$Sens.test$SS3%>%
            mutate(Species=capitalize(Neim))
          Out.Scens=Scens%>%
            mutate(h.mean=NA,h.sd=NA,Initial.bio.low=NA,Initial.bio.up=NA,Final.bio.low=NA,Final.bio.up=NA)
          Store.sens=vector('list',nrow(Scens))
          names(Store.sens)=Scens$Scenario
          Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.probs.B.Bmsy=Out.f.series=Out.B.Bmsy=
            Out.F.Fmsy=store.warnings=store.convergence=Out.ensemble=Out.accept.rate=vector('list',length(Store.sens))
          names(Out.ensemble)=names(Store.sens)
          for(s in 1:length(Store.sens))
          {
            this.wd1=paste(this.wd,Scens$Scenario[s],sep='/')
            if(!dir.exists(this.wd1))dir.create(this.wd1)
            
            n.sims_SS=Scens$Sims[s]
            Report=vector('list',n.sims_SS)
            h.sim=dep.f=m.sim=rep(NA,n.sims_SS)
            
            #get h beta distribution pars
            # h.trunc.normal=rtruncnorm(1e3,a=Min.h.shark,b=Max.h.shark,
            #                           mean=Out.Scens$Steepness[s],sd=max(Out.Scens$Steepness.sd[s],0.01))
            # h.beta=get_beta(mean(h.trunc.normal),sd(h.trunc.normal)/h.trunc.normal)
            h.beta=get_beta(Out.Scens$Steepness[s],max(Out.Scens$Steepness.sd[s],0.01))
            
            #get M lognormal distribution pars
            m.mean=Scens$Mmean[s]
            m.sd=Scens$Msd[s]
            
            #loop over simulations only accepting good runs
            n=trials=1
            while(n <=n.sims_SS)  
            {
              print(paste("___________","SS3 Scenario",Scens$Scenario[s],"___________",Neim," --  Simulation =",n ))
              
              #1. Get random values of interest
              if(do.random.h)
              {
                h.sim[n]=rtbeta(1,h.beta[1],h.beta[2],a=Min.h.shark,b=Max.h.shark)
                m.sim[n]=rlnorm(1,mean=log(m.mean),sd=m.sd/m.mean)
              }else
              {
                h.sim[n]=LH.sim$h[n]
                if(Keep.species[i]=="grey nurse shark" & resample.h.greynurse)  #too low h for greynurse
                {
                  h.sim[n]=rtruncnorm(1,a=Min.h.shark,b=Max.h.shark, 
                                      mean=Out.Scens$Steepness[s],
                                      sd=max(Out.Scens$Steepness.sd[s],0.01))
                }
                m.sim[n]=LH.sim$M[n] 
              }
              
              Depletion.year=ktch$finyear[nrow(ktch)]
              Btklow=List.sp[[i]]$FINALBIO[1]
              # if(Neim=="milk shark")
              # {
              #   Depletion.year=1991  #set to end of catch period to allow convergence
              #   Btklow=0.2
              # }
              
              dep.f[n]=runif(1,Btklow,List.sp[[i]]$FINALBIO[2])
              Scenario=Scens[s,]%>%
                mutate(Mmean=m.sim[n],
                       Steepness=h.sim[n],
                       Initial.dpl=runif(1,List.sp[[i]]$STARTBIO[1],List.sp[[i]]$STARTBIO[2]),
                       Final.dpl=dep.f[n],
                       Ln_R0_init=runif(1,Ln_R0_min,Ln_R0_max))
              Life.history=List.sp[[i]]
              if(do.random.h)
              {
                Life.history$Fecundity=mean(Life.history$Fecundity)
                Life.history$Max.age.F=Life.history$Max.age.F[1]
                Life.history$Breed.cycle=mean(Life.history$Breed.cycle)
              }else
              {
                if(!Keep.species[i]=="grey nurse shark")    #no convergence with re-sampled pars, misspecification in life history
                {
                  Life.history$Fecundity=LH.sim$Fecundity[n]
                  Life.history$Max.age.F=LH.sim$Max.age[n]
                  Life.history$Breed.cycle=LH.sim$Rep.cycle[n]
                  Life.history$Growth.F[1:2]=c(LH.sim$k[n],LH.sim$Linf[n])
                }
              }
              
              
              #2. Create input files
              fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SSS'),
                           new.path=this.wd1,
                           Scenario=Scenario,
                           Catch=ktch,
                           life.history=Life.history,
                           depletion.yr=Depletion.year,
                           Future.project=add.ct.future) 
              
              #3. Run SS3
              fn.run.SS(where.inputs=this.wd1,
                        where.exe=handl_OneDrive('SS3/ss_win.exe'),
                        args='-nohess')  #no Hessian estimation as uncertainty generated thru MC procedure
              
              #4. Bring in outputs
              Report[[n]]=SS_output(this.wd1,covar=F,forecast=F,readwt=F,checkcor=F)
              #ss.std=read.table(paste(this.wd1,'ss.std',sep='/')) https://vlab.noaa.gov/web/stock-synthesis/public-forums/-/message_boards/view_message/11664630
              Report[[n]]$Final.dpl=data.frame(Distribuion=c('Prior','Posterior'),
                                               Value=c(Scenario$Final.dpl,Report[[n]]$current_depletion),
                                               Parameter=rep('Final depletion',2))
              Report[[n]]$Initial.dpl=data.frame(Distribuion=c('Prior','Posterior'),
                                                 Value=c(Scenario$Initial.dpl,Report[[n]]$derived_quants[grep('Bratio',Report[[n]]$derived_quants$Label)[1],'Value']),
                                                 Parameter=rep('Initial depletion',2))
              Report[[n]]$Ln.R0=data.frame(Distribuion=c('Prior','Posterior'),
                                           Value=with(Report[[n]]$estimated_non_dev_parameters,c(Init,Value)),
                                           Parameter=rep('Ln.R0',2))
              
              #5. Keep only runs that converged
              Converged=data.frame(Max.grad=Report[[n]]$maximum_gradient_component,
                                   dep.f_estim=Report[[n]]$current_depletion,
                                   dep.f_obs=dep.f[[n]],
                                   LnRo=Report[[n]]$estimated_non_dev_parameters$Value,
                                   LnRo.min=Report[[n]]$estimated_non_dev_parameters$Min,
                                   LnRo.max=Report[[n]]$estimated_non_dev_parameters$Max)%>%
                mutate(Delta.conv.criteria=abs(Max.grad - 1e-4),
                       Delta.fin.dep=abs(dep.f_obs-dep.f_estim),
                       Delta.lower.LnRo=LnRo/LnRo.min,
                       Delta.upper.LnRo=LnRo/LnRo.max,
                       Keep=ifelse(Delta.fin.dep>Criteria.delta.fin.dep |
                                     Delta.upper.LnRo>0.9 | Delta.lower.LnRo<5,'No','Yes'))
              #     Keep=ifelse(Delta.conv.criteria>0.001 | Delta.fin.dep>Criteria.delta.fin.dep |
              #          Delta.upper.LnRo<0.5 | Delta.lower.LnRo>0.5,'No','Yes'))
              trials=trials+1
              if(Converged$Keep=='Yes')
              {
                Add.bmsy=Report[[n]]$derived_quants%>%filter(Label%in%c('SSB_Virgin','SSB_MSY'))
                Add.bmsy.val=Add.bmsy[Add.bmsy$Label=='SSB_MSY','Value']/Add.bmsy[Add.bmsy$Label=='SSB_Virgin','Value']
                Add.bmsy=Report[[n]]$estimated_non_dev_parameters[1,]
                rownames(Add.bmsy)='bmsy/k'
                Add.bmsy[,]=NA
                Add.bmsy=Add.bmsy%>%
                  mutate(Value=Add.bmsy.val)
                Report[[n]]$estimated_non_dev_parameters=rbind(Report[[n]]$estimated_non_dev_parameters,Add.bmsy)
                
                n=n+1
              }
            }
            
            #Store estimates
            Estims=do.call(rbind,fn.get.stuff.from.list(Report,"estimated_non_dev_parameters"))
              
            Redundant=TRUE
            if(!Redundant)
            {
              store.convergence[[s]]=data.frame(m=m.sim)%>%
                mutate(h=h.sim,
                       Max.grad=unlist(fn.get.stuff.from.list(Report,'maximum_gradient_component')),
                       Delta.conv.criteria=abs(Max.grad-1e-4),
                       dep.f=dep.f,
                       Current.depletion=unlist(fn.get.stuff.from.list(Report,'current_depletion')),
                       Delta.fin.dep=abs(dep.f-Current.depletion),
                       LnRo=Estims$Value,
                       Delta.lower.LnRo=abs(LnRo-1),
                       Delta.upper.LnRo=abs(LnRo-20),
                       Keep=ifelse(Delta.fin.dep>Criteria.delta.fin.dep,'No','Yes'))
              #Keep=ifelse(Delta.conv.criteria>0.001 | Delta.fin.dep>Criteria.delta.fin.dep |
              #              Delta.upper.LnRo<1 | Delta.lower.LnRo<0.5,'No','Yes'))
              
              #show retained and discarded runs
              p=store.convergence[[s]]%>%
                mutate(Keep2=case_when(Keep =='No'  & Delta.conv.criteria>0.001 ~ "No - convergence",
                                       Keep =='No'  & Delta.fin.dep>Criteria.delta.fin.dep ~ "No - final depletion",
                                       Keep =='No'  & Delta.upper.LnRo<1 ~ "No - LnRo upper bound",
                                       Keep =='No'  & Delta.lower.LnRo<0.5 ~ "No - LnRo lower bound",
                                       TRUE ~ Keep))%>%
                ggplot(aes(h,dep.f,color=LnRo))+
                geom_point()+
                facet_wrap(~Keep2,ncol=1)+ggtitle(Keep.species[i])+
                theme(legend.position="top")+ylim(0,1)+geom_hline(yintercept = 0.2)
              pdf(paste(this.wd1,'Retained_discarded_runs.pdf',sep='/'))
              print(p)
              dev.off()
              clear.log('p')
              
              Report=Report[which(store.convergence[[s]]$Keep=='Yes')]
              
            }
            
            if(length(Report)>0)
            {
              #Store acceptance rate
              Out.accept.rate[[s]]=data.frame(Acceptance=100*n.sims_SS/trials, 
                                              Scenario=Scens$Scenario[s])
              
              #Store estimates (recalculate with retained runs)
              Estims=do.call(rbind,fn.get.stuff.from.list(Report,"estimated_non_dev_parameters"))
              Out.estimates[[s]]=Estims%>%
                mutate(Par=gsub('[0-9]+', '', rownames(Estims)))%>%
                group_by(Par)%>%
                summarise(Median=median(Value),
                          Lower.95=quantile(Value,probs=0.025),
                          Upper.95=quantile(Value,probs=0.975))%>%
                dplyr::select(Par,Median,Lower.95,Upper.95)
              
              add.MSY=do.call(rbind,fn.get.stuff.from.list(Report,"derived_quants"))%>%
                filter(grepl('Dead_Catch_MSY',Label))%>%
                rename(Par=Label)%>%
                mutate(Scenario=Scens$Scenario[s])
              MSY.ensemble=add.MSY%>%dplyr::select(Par,Value,Scenario)%>%mutate(Par='MSY')
              add.MSY=add.MSY%>%
                group_by(Par)%>%
                summarise(Median=median(Value),
                          Lower.95=quantile(Value,probs=0.025),
                          Upper.95=quantile(Value,probs=0.975))%>%
                dplyr::select(Par,Median,Lower.95,Upper.95)   
              Out.estimates[[s]]=rbind(Out.estimates[[s]],add.MSY)
              
              #Store Scenarios
              Out.Scens$h.mean[s]=Out.Scens$Steepness[s]
              Out.Scens$h.sd[s]=max(Out.Scens$Steepness.sd[s],0.01)
              Out.Scens$Initial.bio.low[s]=List.sp[[i]]$STARTBIO[1]
              Out.Scens$Initial.bio.up[s]=List.sp[[i]]$STARTBIO[2]
              Out.Scens$Final.bio.low[s]=List.sp[[i]]$FINALBIO[1]
              Out.Scens$Final.bio.up[s]=List.sp[[i]]$FINALBIO[2]
              
              
              #Store trajectories  
              dummy=fn.ktch.only.get.timeseries(d=Report,
                                                mods=names(Catch_only)[w],
                                                Type='Depletion',
                                                scen=Scens$Scenario[s],
                                                Katch=outktch)
              Out.rel.biom[[s]]=dummy$Dat
              Out.probs.rel.biom[[s]]=dummy$Probs
              
              dummy=fn.ktch.only.get.timeseries(d=Report,
                                                mods=names(Catch_only)[w],
                                                Type='F.series',
                                                scen=Scens$Scenario[s],
                                                Katch=outktch)
              Out.f.series[[s]]=dummy$Dat
              
              dummy=fn.ktch.only.get.timeseries(d=Report,
                                                mods=names(Catch_only)[w],
                                                Type='B.Bmsy',
                                                scen=Scens$Scenario[s],
                                                Katch=outktch)
              Out.B.Bmsy[[s]]=dummy$Dat
              Out.probs.B.Bmsy[[s]]=dummy$Probs
              
              dummy=fn.ktch.only.get.timeseries(d=Report,
                                                mods=names(Catch_only)[w],
                                                Type='F.Fmsy',
                                                scen=Scens$Scenario[s],
                                                Katch=outktch)
              Out.F.Fmsy[[s]]=dummy$Dat
              
              
              #Store quantities for ensemble model
              Years=c(ktch$finyear,add.ct.future$finyear)  
              dum=do.call(rbind,fn.get.stuff.from.list(Report,"derived_quants"))
              
              SSB_Virgin=dum[grep("SSB_Virgin",dum$Label),]
              SSB_Virgin=SSB_Virgin%>%
                mutate(Simulation=readr::parse_number(paste(rownames(SSB_Virgin),'0',sep='')))%>%
                rename(Value.virgin=Value)%>%
                dplyr::select(Value.virgin,Simulation)
              
              SSB_MSY=dum[grep("SSB_MSY",dum$Label),]
              SSB_MSY=SSB_MSY%>%
                mutate(Simulation=readr::parse_number(paste(rownames(SSB_MSY),'0',sep='')))%>%
                rename(Value.MSY=Value)%>%
                dplyr::select(Value.MSY,Simulation)
              
              annF_MSY=dum[grep("annF_MSY",dum$Label),]
              annF_MSY=annF_MSY%>%
                mutate(Simulation=readr::parse_number(paste(rownames(annF_MSY),'0',sep='')))%>%
                rename(Value.F_MSY=Value)%>%
                dplyr::select(Value.F_MSY,Simulation)  
              
              
              SSB=dum[grep(paste(paste("SSB",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
              SSB=SSB%>%
                mutate(year=readr::parse_number(Label),
                       Simulation=readr::parse_number(paste(rownames(SSB),'0',sep='')),
                       Simulation=as.numeric(substr(Simulation,5,10)))%>%
                dplyr::select(-Label)
              
              EF=dum[grep(paste(paste("F",Years,sep='_'),collapse="|"),dum$Label),c('Label','Value')]
              EF=EF%>%
                mutate(year=readr::parse_number(Label),
                       Simulation=readr::parse_number(paste(rownames(EF),'0',sep='')),
                       Simulation=as.numeric(substr(Simulation,5,10)))%>%
                dplyr::select(-Label)
              Misn.yr=which(!paste(SSB$year,SSB$Simulation)%in%paste(EF$year,EF$Simulation)) 
              All.misn.yrs=SSB$year[Misn.yr]
              Misn.yr_no.fishing=All.misn.yrs[which(All.misn.yrs<min(EF$year))]
              if(length(Misn.yr_no.fishing)>0)
              {
                add.dis=EF[1:length(Misn.yr_no.fishing),]%>%
                  mutate(Value=0,
                         Simulation=rep(sort(unique(SSB$Simulation)),each=length(unique(Misn.yr_no.fishing))),
                         year=Misn.yr_no.fishing)
                kkk=str_sub(add.dis$Simulation,1,-2)
                rownames(add.dis)=paste('F_',add.dis$year,kkk,sep='')
                EF=rbind(EF,add.dis)%>%
                  arrange(Simulation,year)
                
                Misn.yr=which(!paste(SSB$year,SSB$Simulation)%in%paste(EF$year,EF$Simulation))
                if(length(Misn.yr)>0)
                {
                  add.dis=EF[Misn.yr,]%>%
                    mutate(Value=NA,
                           Simulation=sort(unique(SSB$Simulation)),
                           year=SSB$year[Misn.yr])
                  rownames(add.dis)=paste('F_',SSB$year[Misn.yr],c('',1:(nrow(add.dis)-1)),sep='')
                  EF=rbind(EF,add.dis)%>%
                    arrange(Simulation,year)%>%
                    mutate(Value = na.approx(Value))
                }
              }

              Depletion=SSB%>%
                left_join(SSB_Virgin%>%dplyr::select(Value.virgin,Simulation),by='Simulation')%>%
                mutate(Depletion=Value/Value.virgin)%>%
                dplyr::select(Simulation,year,Depletion)%>%
                spread(year,Depletion)%>%
                dplyr::select(-Simulation)
              
              B.Bmsy=SSB%>%
                left_join(SSB_MSY%>%dplyr::select(Value.MSY,Simulation),by='Simulation')%>%
                mutate(B_Bmys=Value/Value.MSY)%>%
                dplyr::select(Simulation,year,B_Bmys)%>%
                spread(year,B_Bmys)%>%
                dplyr::select(-Simulation)
              
              F.Fmsy=EF%>%
                left_join(annF_MSY%>%dplyr::select(Value.F_MSY,Simulation),by='Simulation')%>%
                mutate(F_Fmys=Value/Value.F_MSY)%>%
                dplyr::select(Simulation,year,F_Fmys)%>%
                spread(year,F_Fmys)%>%
                dplyr::select(-Simulation)
              Out.ensemble[[s]]=list(Depletion=Depletion,
                                     B.Bmsy=B.Bmsy,
                                     F.Fmsy=F.Fmsy,
                                     MSY=MSY.ensemble) 
              rm(Depletion,B.Bmsy,F.Fmsy,dum,MSY.ensemble)
              
              
              #Flag warnings
              n.sims_SS=length(Report)
              warnings=fn.get.stuff.from.list(Report,"warnings")
              dummy.warnings=vector('list',n.sims_SS)
              for(n in 1:n.sims_SS)
              {
                dummy.warnings[[n]]=data.frame(warning=warnings[[n]])%>%
                  mutate(warning.type=case_when(grepl('bad Beta prior',warning)~'bad Beta prior',
                                                grepl('Final gradient',warning)~'Final gradient',
                                                grepl('midsize bin',warning)~'midsize bin',
                                                grepl('Fmsy/mey is <',warning)~'Fmsy/mey',
                                                grepl('Maximum pop size bin',warning)~'Maximum pop size bin'))%>%
                  distinct(warning.type,.keep_all = T)%>%
                  filter(!is.na(warning.type))%>%
                  mutate(Steepness=h.sim[n],
                         Final.depl=dep.f[n])
                
                
              }
              store.warnings[[s]]=do.call(rbind,dummy.warnings)
              
              TAB=data.frame(Bad.beta.prior=rep(NA,n.sims_SS))%>%
                mutate(Final.grad=Bad.beta.prior,
                       midsize.bin=Bad.beta.prior)
              for(n in 1:n.sims_SS)
              {
                TAB$Bad.beta.prior[n]=any(grepl('bad Beta prior',Report[[n]]$warnings))
                TAB$Final.grad[n]=any(grepl('Final gradient',Report[[n]]$warnings))
                TAB$midsize.bin[n]=any(grepl('midsize bin',Report[[n]]$warnings))
              }
              
            }
            
            out=list(Final.dpl=do.call(rbind,fn.get.stuff.from.list(Report,'Final.dpl')),
                     Initial.dpl=do.call(rbind,fn.get.stuff.from.list(Report,'Initial.dpl')),
                     Ln.R0=do.call(rbind,fn.get.stuff.from.list(Report,'Ln.R0'))) 
            out$Final.dpl=out$Final.dpl%>%
              mutate(Value=ifelse(Distribuion=='Prior',runif(nrow(out$Final.dpl),List.sp[[i]]$FINALBIO[1],List.sp[[i]]$FINALBIO[2]),
                                  Value))
            out$Initial.dpl=out$Initial.dpl%>%
              mutate(Value=ifelse(Distribuion=='Prior',runif(nrow(out$Initial.dpl),List.sp[[i]]$STARTBIO[1],List.sp[[i]]$STARTBIO[2]),
                                  Value))
            fn.show.density(d=do.call(rbind,out),NCOL=2) 
            ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                         capitalize(names(dummy.store)[i]),"/",AssessYr,"/",'SS3',"/Prior.and.posterior_",Scens$Scenario[s],".tiff",sep=''),
                   width = 12,height = 14, dpi = 300, compression = "lzw")
            
            rm(this.wd1,Report)
          } # end s loop
          
          #output quantities for ensemble
          Out.ens=list(Depletion=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'Depletion')),
                       B.Bmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'B.Bmsy')),
                       F.Fmsy=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'F.Fmsy')),
                       MSY=do.call(rbind,fn.get.stuff.from.list(Out.ensemble,'MSY')))
          if(nrow(Scens)>1)
          {
            for(f in 1:length(Out.ens))   
            {
              ddmi=Out.ens[[f]]
              TAB=table(sub("\\..*", "", rownames(ddmi)))
              Out.ens[[f]]=ddmi[sample(1:nrow(ddmi),size=max(TAB),replace =FALSE),]
            }
          }
          
          #output scenarios
          Out.Scens=Out.Scens%>%
            mutate(Mmean=round(Mmean,3),
                   h.mean=round(h.mean,3),
                   h.sd=round(h.sd,3))%>%
            dplyr::select(-c(Model,Final.dpl, Steepness,Steepness.sd,use_F_ballpark, Sims))%>%
            relocate(Species,Scenario,Mmean,F_ballpark,h.mean,h.sd,
                     Initial.bio.low,Initial.bio.up,Final.bio.low,Final.bio.up)
          dummy.store.sens.table[[i]]=Out.Scens
          dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
          dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
          dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
          dummy.store.probs.B.Bmsy[[i]]=Out.probs.B.Bmsy
          dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
          dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
          dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
          dummy.store.warnings[[i]]=do.call(rbind,store.warnings)
          dummy.store.convergence[[i]]=do.call(rbind,store.convergence)
          dummy.store.ensemble[[i]]=Out.ens
          dummy.store.accept.rate[[i]]=do.call(rbind,Out.accept.rate)
          
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
          ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                       capitalize(List.sp[[i]]$Name),"/",AssessYr,"/SS3/Rel.biomass&estimates.tiff",sep=''),compression = "lzw")
        }
      } #end i
      
      Catch_only[[w]]$sens.table=dummy.store.sens.table
      Catch_only[[w]]$estimates=dummy.store.estimates
      Catch_only[[w]]$rel.biom=dummy.store.rel.biom
      Catch_only[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
      Catch_only[[w]]$probs.B.Bmsy=dummy.store.probs.B.Bmsy
      Catch_only[[w]]$f.series=dummy.store.f.series
      Catch_only[[w]]$B.Bmsy=dummy.store.B.Bmsy
      Catch_only[[w]]$F.Fmsy=dummy.store.F.Fmsy
      Catch_only[[w]]$Kobe.probs=dummy.store.Kobe.probs
      Catch_only[[w]]$ensemble=dummy.store.ensemble
      Catch_only[[w]]$warnings=dummy.store.warnings
      Catch_only[[w]]$convergence=dummy.store.convergence
      Catch_only[[w]]$accept.rate=dummy.store.accept.rate
      
      rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,
         dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,
         dummy.store.B.Bmsy,dummy.store.F.Fmsy,dummy.store.ensemble,
         dummy.store.warnings,dummy.store.convergence,dummy.store.accept.rate,
         dummy.store.probs.B.Bmsy)
      
    }
  }
    
}
toc(log = TRUE, quiet = TRUE)
computation.time <- tic.log(format = TRUE)
tic.clearlog()
send.email(TO=Send.email.to,
           CC='',
           Subject=paste("Catch Only models finished running at",Sys.time()),
           Body= paste("Computation",computation.time),  
           Attachment=NULL) 


#18.2. COM weighted average  
  #18.2.1. get COM weights
#r groups
All.rs=do.call(rbind,store.species.r_M.min)
r.list=All.rs$mean
names(r.list)=1:N.sp
r.groups=list(low=c(min(r.list),quantile(r.list,0.2499)),
              medium=c(quantile(r.list,0.25),quantile(r.list,0.7499)),
              high=c(quantile(r.list,0.75),0.65))
if(do.ensemble.simulations)   #56 sec per iteration under do.full.sims=FALSE (iterations= r.groups(n=3) * HarvestRates (n=4) * n.reps (20))   
{
  #Ks
  K.list=c(1e4)
  
  #Init depletion
  Init.depl=1 #Init.depl=c(1,0.7)
  
  #do.full.sims=TRUE
  do.full.sims=FALSE
  
  #Create exploitation histories
  #note: setting U at quantile(r.list,0.8) creates a range of exploitation histories ranging
  #       between 10 times Fmsy and Fmsy depending on life history
  YRS=1974:2021
  HarvestRate1=seq(0,quantile(r.list,0.8),length.out = (length(YRS)/2))+rnorm(length(YRS)/2,0,.01)
  HarvestRate2=seq(0,quantile(r.list,0.5),length.out = (length(YRS)/2))+rnorm(length(YRS)/2,0,.01)
  
  Inputd=data.frame(YRS=YRS)%>%
    mutate(HarvestRate1=c(HarvestRate1,rev(HarvestRate1)),
           HarvestRate2=c(HarvestRate2,rev(HarvestRate2)),
           #HarvestRate1=rep(quantile(r.list,0.5),length(YRS))+rnorm(length(YRS),0,.01),  #constant U is not representative
           #HarvestRate2=rep(quantile(r.list,0.1),length(YRS))+rnorm(length(YRS),0,.01),
           HarvestRate3=HarvestRate1, 
           HarvestRate4=HarvestRate2)
  
  n=round(length(YRS)/6)
  Inputd$HarvestRate1[(nrow(Inputd)-n):nrow(Inputd)]=Inputd$HarvestRate1[(nrow(Inputd)-n)]
  Inputd$HarvestRate2[(nrow(Inputd)-n):nrow(Inputd)]=Inputd$HarvestRate2[(nrow(Inputd)-n)]
  
  n=1:ceiling(length(YRS)*.35)
  Inputd$HarvestRate3[n]=0
  Inputd$HarvestRate4[n]=0
  n=(nrow(Inputd)-ceiling(length(YRS)*.35)):nrow(Inputd)
  Inputd$HarvestRate3[n]=0
  Inputd$HarvestRate4[n]=0
  Inputd$HarvestRate3=ifelse(Inputd$HarvestRate3<0,0,Inputd$HarvestRate3)
  Inputd$HarvestRate4=ifelse(Inputd$HarvestRate4<0,0,Inputd$HarvestRate4)
  
  Inputd=Inputd%>%
    mutate(HarvestRate1=HarvestRate1+rnorm(nrow(Inputd),0,.01),
           HarvestRate2=HarvestRate2+rnorm(nrow(Inputd),0,.01),
           HarvestRate3=HarvestRate3+rnorm(nrow(Inputd),0,.01),
           HarvestRate4=HarvestRate4+rnorm(nrow(Inputd),0,.01),
           HarvestRate1=ifelse(HarvestRate1<0,0,HarvestRate1),
           HarvestRate2=ifelse(HarvestRate2<0,0,HarvestRate2),
           HarvestRate3=ifelse(HarvestRate3<0,0,HarvestRate3),
           HarvestRate4=ifelse(HarvestRate4<0,0,HarvestRate4))
  
  #Factorial design
  Factorial=expand.grid(r=names(r.groups),K=K.list,p=Init.depl,ktch=paste('HarvestRate',1:(ncol(Inputd)-1),sep=''))
  Factorial$ktch=as.character(Factorial$ktch)
  
  #Test if OM can retrieve assumed pars
  this.wd='C:/Superensemble'
  setwd(this.wd)
  check.can.rekv.pars=FALSE
  if(check.can.rekv.pars)
  {
    Q=1e-3
    store.estimpars=store.inpupars=store.ktch=store.bt=vector('list',length=nrow(Factorial))
    Thisvec=c(1,10,22,31,43,73,85,94,106,115,127,136,148,157)
    for(s in Thisvec)
    {
      HRscen=Factorial$ktch[s]
      PARS=log(c(K=Factorial$K[s],r=Factorial$r[s],Init.dep=Factorial$p[s],q=Q))
      
      what='simulate'
      OBJ=SPM.fit(PARS)
      store.ktch[[s]]=OBJ$ktch
      store.bt[[s]]=OBJ$Bt
      Obs_CPUE=OBJ$cpue
      ln_CPUE <- log(Obs_CPUE+rnorm(length(Obs_CPUE),0,1)) 
      
      what='fit'
      init.par=PARS+rnorm(length(PARS),0,.5)
      nlmb <- nlminb(init.par, SPM.fit, gradient = NULL,hessian = TRUE)
      store.estimpars[[s]]=nlmb$par
      store.inpupars[[s]]=PARS
      
      #opt <- optim(par=init.par,fn=SPM.fit,method="Nelder-Mead",hessian=TRUE)
      #solve(opt$hessian) 
    }
    
    fn.fig("Recover.pars_biom&ktch.fit",2400,2400)
    par(mfcol=c(5,3),mai=c(.3,.5,.2,.5))
    for(s in Thisvec)
    {
      plot(store.ktch[[s]],cex=2,type='l',ylab='')
      par(new=T)
      plot(store.bt[[s]],col=2,type='l',axes = FALSE, xlab = '', ylab ='')
      axis(side=4, at = pretty(range(store.bt[[s]])))
      mtext(paste('K=',Factorial$K[s],'r=',round(Factorial$r[s],2),
                  'Init.dep=',Factorial$p[s],Factorial$ktch[s]),3,cex=.55)
      
    }
    mtext("Catch", side=2,-1.5,outer=T,cex=1.5)
    mtext("Biomass", side=4, line=-1,col=2,outer=T,cex=1.5)
    dev.off()
    
    fn.fig("Recover.pars_pars.fit",2400,2400)
    par(mfcol=c(5,3),mar=rep(1,4),mai=c(.3,.5,.2,.5))
    for(s in Thisvec)
    {
      plot(store.estimpars[[s]],cex=2,pch=19,ylab='Ln.par')
      points(store.inpupars[[s]],col=2,cex=2)
      mtext(paste('K=',Factorial$K[s],'r=',round(Factorial$r[s],2),
                  'Init.dep=',Factorial$p[s],Factorial$ktch[s]),3,cex=.6)
      
    }
    legend('topright',c('Estimate','Input'),pch=19,col=1:2,bty='n')
    dev.off()
  }
  
  #1. Simulate biomass trajectories
  n.reps=10
  Store.sims=vector('list',nrow(Factorial)*n.reps)
  nnn=seq(0,(length(Store.sims)-n.reps),n.reps)
  for(s in 1:nrow(Factorial))
  {
    for(i in 1:n.reps)
    {
      ar=r.groups[[match(Factorial$r[s],names(r.groups))]]
      ar=runif(1,ar[1],ar[2])
      ar=log(ar)
      x <- exp(rmvnorm(n=1, mean=c(ar,log(Factorial$K[s])), 
                       sigma=matrix(c(0.0002454185,-0.0001840755,-0.0001840755,0.0001383331),
                                    nrow=2))) #variance cov variance from fitting SPM to real data
      
      x[,1]=min(max(min(r.groups$low),x[,1]),r.groups$high[2])
      dumi=SPM(K=x[,2],r=x[,1],Init.dep=Factorial$p[s],HRscen=Factorial$ktch[s])
      s.i=nnn[s]+i
      Store.sims[[s.i]]=dumi
    }
  }
  Inputs=Store.sims
  for(s in 1:length(Inputs))
  {
    Inputs[[s]]=data.frame(Ktch=Store.sims[[s]]$Ktch)%>%
                          mutate(r=Store.sims[[s]]$r,
                                 K=Store.sims[[s]]$K, 
                                 Init.dep=Store.sims[[s]]$Init.dep,
                                 HRscen=Store.sims[[s]]$HRscen,
                                 Depletion=Store.sims[[s]]$Depletion, 
                                 B.Bmsy=Store.sims[[s]]$B.Bmsy,
                                 F.Fmsy=Store.sims[[s]]$F.Fmsy,
                                 Year=YRS)
  }
  
  #2. Fit COMs to catch trajectories as per stock assessments 83 sec per simulation
  get.rs=data.frame(index=names(r.list),r=r.list)  
  system.time({for(s in 1:length(Inputs))   
  {
    Ktch=Inputs[[s]]$Ktch
    Klow=k.fun(Ktch,5) 
    Kup=k.fun(Ktch,25)
    
    Yrs=Inputs[[s]]$Year
    dis.r.range=do.call(rbind,r.groups)%>%   
      data.frame%>%
      mutate(Get=between(unique(Inputs[[s]]$r),.[[1]] , .[[2]]))%>%
      filter(Get==TRUE)
    Inputs[[s]]$r.group=row.names(dis.r.range)
    i=get.rs%>%
      filter(r>=dis.r.range[,1] & r<=dis.r.range[,2])
    i=as.numeric(i[sample(i$index,1),'index'])
    
    Btklow=0  #must allow highly depleted stock 
    Btkup=1
    #Btklow=List.sp[[i]]$FINALBIO[1]
    #Btkup=List.sp[[i]]$FINALBIO[2]
    b1k.low=List.sp[[i]]$STARTBIO[1]
    b1k.up=List.sp[[i]]$STARTBIO[2]
    
    
    #1. DBSRA
    print(paste("___________","DBSRA ","______ s=",s))
    Scens=List.sp[[i]]$Sens.test$DBSRA%>%filter(Scenario=='S1')
    
    #priors
    AgeMat=Scens$AgeMat
    Mmean=Scens$Mmean  
    Msd=Scens$Msd
    Mcv=Msd/Mmean
    M.dist="lnorm"
    fmsy.m=Scens$fmsy.m
    fmsym.dist="lnorm"
    fmsym.CV=0.1
    bmsyk.mean=Scens$bmsyk
    bmsyk.dist="beta"
    bmsyk.CV=0.1
    btk.dist="unif"
    b1k.dist="unif"
    
    if(do.full.sims) SIMS=Scens$Sims else
    SIMS=1e3
    
    #run model
    out.dummy=apply.DBSRA(year=Yrs,
                          catch=Ktch,
                          catchCV=NULL,  
                          catargs=list(dist="none",low=0,up=Inf,unit="MT"),  #catch CV not available
                          agemat=AgeMat,
                          k=list(low=Klow,up=Kup,tol=0.01,permax=1000),
                          b1k=list(dist=b1k.dist,low=b1k.low,up=b1k.up,mean=1,sd=0.1),  #mean and sd not used if 'unif'
                          btk=list(dist="unif",low=Btklow,up=Btkup,mean=1,sd=0.1,refyr=max(Yrs)),  #reference year
                          fmsym=list(dist="lnorm",low=0.1,up=2,mean=log(fmsy.m),sd=fmsym.CV), # Cortes & Brooks 2018. Low and up not used if 'lnorm'  
                          bmsyk=list(dist="beta",low=0.05,up=0.95,mean=bmsyk.mean,sd=bmsyk.CV),  
                          M=list(dist="lnorm",low=0.001,up=1,mean=log(Mmean),sd=Mcv),
                          graph=1,
                          nsims=SIMS, 
                          grout=1,
                          WD=this.wd,
                          outfile="Appendix_fit")
    
    Dep=apply(out.dummy$output$Depletion.traj,2,median)
    #B.Bmsy=apply(out.dummy$output$B.Bmsy,2,median)
    #F.Fmsy=apply(out.dummy$output$F.Fmsy,2,median)
    
    Inputs[[s]]$Depletion_DBSRA=Dep[1:length(Dep)-1]
    #Inputs[[s]]$B.Bmsy_DBSRA=B.Bmsy[1:length(Dep)-1]
    #Inputs[[s]]$F.Fmsy_DBSRA=F.Fmsy[1:length(Dep)-1]
    
    #Inputs[[s]]$Depletion_DBSRA=mean(Dep[(length(Dep)-4):length(Dep)],na.rm=T)
    # Inputs[[s]]$B.Bmsy_DBSRA=mean(B.Bmsy[(length(B.Bmsy)-4):length(B.Bmsy)],na.rm=T)
    # Inputs[[s]]$F.Fmsy_DBSRA=mean(F.Fmsy[(length(F.Fmsy)-4):length(F.Fmsy)],na.rm=T)
    
    rm(Dep,out.dummy)
    
    
    #2. CMSY
    print(paste("___________","CMSY ","______ s=",s))
    
    #priors
    Scens=List.sp[[i]]$Sens.test$CMSY%>%filter(Scenario=='S1')
    RES=RESILIENCE[[i]]
    if(Scens$r.prob.min==0)
    {
      r.range=NA
      k.range=NA
      Bf.low=NA
      Bf.hi=NA
    }else
    {
      Mn=min(Scens$r,Max.r.value)  #some life history pars yield unrealistically high r
      r.range=quantile(rnorm(1e3,
                             mean=Mn,
                             sd=Scens$r.sd),
                       probs=c(Scens$r.prob.min,Scens$r.prob.max))
      r.range[1]=max(r.range[1],Min.r.value)
      k.range=c(Klow,Kup)
    }
    Proc.error=Scens$Proc.error
    if(do.full.sims) SIMS=Scens$Sims else
      SIMS=1e3
    
    #run model
    out.dummy=apply.CMSY(year=Yrs,
                         catch=Ktch,
                         r.range=r.range,
                         k.range=k.range,
                         Bo.low=b1k.low,
                         Bo.hi=b1k.up,
                         Bf.low=Btklow,
                         Bf.hi=Btkup,
                         outfile=paste(this.wd,'Appendix_fit',sep='/'),
                         nsims=SIMS,
                         Proc.error=Proc.error,
                         RES=RES)
    
    Dep=apply(out.dummy$output$Depletion.traj,2,median)
    #B.Bmsy=apply(out.dummy$output$B.Bmsy,2,median)
    #F.Fmsy=apply(out.dummy$output$F.Fmsy,2,median)
    
    Inputs[[s]]$Depletion_CMSY=Dep[1:length(Dep)-1]
    #Inputs[[s]]$B.Bmsy_CMSY=B.Bmsy[1:length(Dep)-1]
    #Inputs[[s]]$F.Fmsy_CMSY=F.Fmsy[1:length(Dep)-1]
    
    # Inputs[[s]]$Depletion_CMSY=mean(Dep[(length(Dep)-4):length(Dep)],na.rm=T)
    # Inputs[[s]]$B.Bmsy_CMSY=mean(B.Bmsy[(length(B.Bmsy)-4):length(B.Bmsy)],na.rm=T)
    # Inputs[[s]]$F.Fmsy_CMSY=mean(F.Fmsy[(length(F.Fmsy)-4):length(F.Fmsy)],na.rm=T)
    
    rm(Dep,out.dummy)
    
    
    #3. JABBA
    print(paste("___________","JABBA ","______ s=",s))
    
    #priors
    Scens=List.sp[[i]]$Sens.test$JABBA_catch.only%>%filter(Scenario=='S1')
    bmsyk.mean=Scens$bmsyk
    Proc.error=Scens$Proc.error
    r.CV.multi=Scens$r.CV.multiplier
    K.prior=c(Klow*20,Scens$K.CV)
    
    Mn=min(Scens$r,Max.r.value)  #some life history pars yield unrealistically high r
    r.prior=c(Mn, (Scens$r.sd/Mn)*Scens$r.CV.multiplier)
    Ktch.CV=Scens$Ktch.CV
    Bint=runif(1000,b1k.low,b1k.up)
    Bint.mean=mean(Bint)
    Bint.CV=sd(Bint)/Bint.mean
    psi.prior=c(Bint.mean,Bint.CV)
    
    Bfin=runif(1000,Btklow,Btkup)
    Bfin.mean=mean(Bfin)          
    Bfin.CV=sd(Bfin)/Bfin.mean
    b.prior=c(Bfin.mean,Bfin.CV,max(Yrs),"bk")
    
    Rdist = "lnorm"
    Kdist="lnorm"  
    PsiDist='beta'
    
    #Put inputs together
    KtchJ=ifelse(Ktch<1,1,Ktch)   #Jabba fails if catches are negligible
    input=list(Ktch=data.frame(Year=Yrs,Total=KtchJ),
               MDL="Schaefer",
               Ktch.CV=Ktch.CV,
               ASS=names(List.sp)[i],
               Rdist = Rdist,
               Rprior = r.prior,
               r.CV.multiplier=r.CV.multi,
               Kdist= Kdist,
               Kprior=K.prior,    
               PsiDist= PsiDist,
               Psiprior=psi.prior,   
               Bprior=b.prior,    
               BMSYK=bmsyk.mean)  
    
    if(do.full.sims) SIMS=Scens$Sims else
      SIMS=5e3
    
    #run model
    THIN=5
    CHAINS=2
    BURNIN=min(0.15*Scens$Sims,5000)
    JABBA.run=apply.JABBA(Ktch=input$Ktch,
                          CPUE=NULL,
                          CPUE.SE=NULL,
                          MDL=input$MDL,
                          Ktch.CV=input$Ktch.CV,
                          ASS=input$ASS,
                          Rdist = input$Rdist,
                          Rprior = input$Rprior,
                          Kdist=input$Kdist,
                          Kprior=input$Kprior,
                          PsiDist=input$PsiDist,
                          Psiprior=input$Psiprior,
                          Bprior=input$Bprior,
                          BMSYK=input$BMSYK,
                          output.dir=this.wd,
                          outfile="Appendix_fit",
                          Sims=SIMS,
                          Proc.error.JABBA=Proc.error,
                          thinning = THIN,
                          nchains = CHAINS,
                          burn.in= BURNIN)
    out.dummy=JABBA.run$fit
    Dep=apply(out.dummy$posteriors$P,2,median)
    #B.Bmsy=out.dummy$timeseries[, , "BBmsy"][,1]
    #F.Fmsy=out.dummy$timeseries[, , "FFmsy"][,1]
    
    Inputs[[s]]$Depletion_JABBA=Dep[-length(Dep)]
    #Inputs[[s]]$B.Bmsy_JABBA=B.Bmsy
    #Inputs[[s]]$F.Fmsy_JABBA=F.Fmsy
    
    #Inputs[[s]]$Depletion_JABBA=mean(Dep[(length(Dep)-4):length(Dep)],na.rm=T)
    # Inputs[[s]]$B.Bmsy_JABBA=mean(B.Bmsy[(length(B.Bmsy)-4):length(B.Bmsy)],na.rm=T)
    # Inputs[[s]]$F.Fmsy_JABBA=mean(F.Fmsy[(length(F.Fmsy)-4):length(F.Fmsy)],na.rm=T)
    
    rm(Dep,out.dummy)
    

    #4. SSS
    print(paste("___________","SSS ","______ s=",s))
    Scens=List.sp[[i]]$Sens.test$SS3%>%filter(Scenario=='S1')
    #random LH samples
    LH.sim=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                          capitalize(List.sp[[i]]$Name),"/",AssessYr,"/steepness/Life.history.csv",sep=''))
    this.wd1=this.wd
    
    if(do.full.sims) n.sims_SS=100 else
    n.sims_SS=40

    Report=vector('list',n.sims_SS)
    h.sim=dep.f=m.sim=rep(NA,n.sims_SS)
    for(n in 1:n.sims_SS)  
    {
      if(do.random.h)
      {
        h.mean=Out.Scens$Steepness[1]
        h.sim[n]=rtruncnorm(1,a=Min.h.shark, b=Max.h.shark,
                            mean=h.mean,
                            sd=max(Out.Scens$Steepness.sd[1],0.01))
        m.mean=List.sp[[i]]$Sens.test$DBSRA$Mmean[1]
        m.sim[n]=rtruncnorm(1,a=0.02, b=5,
                            mean=m.mean,
                            sd=List.sp[[i]]$Sens.test$DBSRA$Msd[1])
      }else
      {
        h.sim[n]=LH.sim$h[n]
        if(Keep.species[i]=="grey nurse shark" & resample.h.greynurse)  #too low h for greynurse
        {
          h.sim[n]=rtruncnorm(1,a=Min.h.shark,b=Max.h.shark, 
                              mean=Scens$Steepness,
                              sd=max(Scens$Steepness.sd,0.01))
        }
        m.sim[n]=LH.sim$M[n] 
      }
      
      dep.f[n]=runif(1,Btklow,Btkup)
      Scenario=Scens%>%
        mutate(Mmean=m.sim[n],
               Steepness=h.sim[n],
               Initial.dpl=runif(1,b1k.low,b1k.up),
               Final.dpl=dep.f[n],
               Ln_R0_init=rtruncnorm(1,a=Ln_R0_min, b=10,mean=4,sd=1))
      
      Life.history=List.sp[[i]]
      if(!Keep.species[i]=="grey nurse shark")    #no convergence with resampled pars, misspecification in life history
      {
        Life.history$Fecundity=LH.sim$Fecundity[n]
        Life.history$Max.age.F=LH.sim$Max.age[n]
        Life.history$Breed.cycle=LH.sim$Rep.cycle[n]
        Life.history$Growth.F[1:2]=c(LH.sim$k[n],LH.sim$Linf[n])
      }
      
      
      #1. Create input files
      fn.set.up.SS(Templates=handl_OneDrive('SS3/Examples/SSS'),  
                   new.path=this.wd1,
                   Scenario=Scenario,
                   Catch=data.frame(SPECIES=1, Name=Keep.species[i],finyear=Yrs,Tonnes=Ktch),     
                   life.history=Life.history,
                   depletion.yr=Yrs[length(Yrs)])
      
      
      #2. Run SS3
      fn.run.SS(where.inputs=this.wd1,
                where.exe=handl_OneDrive('SS3/ss_win.exe'),
                args='-nohess')  #no Hessian estimation as uncertainty generated thru MC procedure
      
      
      #3. Bring in outputs
      Report[[n]]=SS_output(this.wd1,covar=F,forecast=F,readwt=F,checkcor=F)
      #ss.std=read.table(paste(this.wd1,'ss.std',sep='/')) https://vlab.noaa.gov/web/stock-synthesis/public-forums/-/message_boards/view_message/11664630
     }
    
    #Store estimates
    Estims=do.call(rbind,fn.get.stuff.from.list(Report,"estimated_non_dev_parameters"))
    
    #Check convergence and keep only runs that converge
    Criteria.delta.fin.dep=0.01  
    if(Keep.species[i]%in%c('grey nurse shark',"pigeye shark"))Criteria.delta.fin.dep=0.1
    if(Keep.species[i]=='milk shark')Criteria.delta.fin.dep=0.2
    store.convergence=data.frame(m=m.sim)%>%
                        mutate(h=h.sim,
                               Max.grad=unlist(fn.get.stuff.from.list(Report,'maximum_gradient_component')),
                               Delta.conv.criteria=abs(Max.grad-1e-4),
                               dep.f=dep.f,
                               Current.depletion=unlist(fn.get.stuff.from.list(Report,'current_depletion')),
                               Delta.fin.dep=abs(dep.f-Current.depletion),
                               LnRo=Estims$Value,
                               Delta.lower.LnRo=abs(LnRo-1),
                               Delta.upper.LnRo=abs(LnRo-20),
                               Keep=ifelse(Delta.conv.criteria>0.001 | Delta.fin.dep>Criteria.delta.fin.dep |
                                             Delta.upper.LnRo<1 | Delta.lower.LnRo<0.5,'No','Yes'))
    Report=Report[which(store.convergence$Keep=='Yes')]
    dummy=fn.ktch.only.get.timeseries(d=Report,
                                      mods="SSS",
                                      Type='Depletion',
                                      scen=Scenario$Scenario,
                                      Katch=Ktch)  
    Inputs[[s]]$Depletion_SSS=dummy$Dat$median
    rm(Report)
    

    
    drop.files=list.files(this.wd, full.names = TRUE)
    drop.files=drop.files[-grep('Recover.pars',drop.files)]
    do.call(file.remove, list(drop.files))
  }})
  
  #Compare COMs estimates
  le.cols=c('grey40',"#F8766D", "#00BFC4", "#7CAE00", "#C77CFF")
  names(le.cols)=c('Operation model','DBSRA','CMSY','JABBA','SSS')
  Dumi=do.call(rbind,Inputs)%>%
    mutate(Grup=paste(HRscen,'&' ,capitalize(r.group),'r'),
           Grup=factor(Grup,levels=c("HarvestRate1 & Low r","HarvestRate1 & Medium r","HarvestRate1 & High r",
                                     "HarvestRate2 & Low r","HarvestRate2 & Medium r","HarvestRate2 & High r",
                                     "HarvestRate3 & Low r","HarvestRate3 & Medium r","HarvestRate3 & High r",
                                     "HarvestRate4 & Low r","HarvestRate4 & Medium r","HarvestRate4 & High r")))
  

  for(y in 1:length(Init.depl)) fn.out.poly(In.dep=Init.depl[y])  
  
  
  #Display harvest rates
  le.cols=c("#F8766D", "#00BFC4", "#7CAE00","#C77CFF")
  names(le.cols)=paste('HarvestRate',1:4,sep='')
  
  Inputd%>%
    ggplot(aes(YRS,HarvestRate1))+
    geom_line(aes(YRS,HarvestRate1,color='HarvestRate1'),size=1.1)+
    geom_line(aes(YRS,HarvestRate2,color='HarvestRate2'),size=1.1)+
    geom_line(aes(YRS,HarvestRate3,color='HarvestRate3'),size=1.1)+
    geom_line(aes(YRS,HarvestRate4,color='HarvestRate4'),size=1.1)+
    scale_color_manual(name='Legend',values = le.cols)+
    theme_PA(axs.t.siz=14,axs.T.siz=16,leg.siz=16)+
    theme(legend.position = 'top',
          legend.title = element_blank())+
    ylab(expression('Harvest rate'~(y^-1)))+xlab('Year')
  ggsave(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'Harvest rates.tiff',sep=''),
         width = 10,height = 10,compression = "lzw")
  
  
  #3. Calculate COMs' weights
  COMs.weight=Inputs
  for(s in 1:length(COMs.weight))
  {
    dummy=Inputs[[s]]
    dummy=dummy[(nrow(dummy)-4):nrow(dummy),]%>%
      mutate(Depletion=mean(Depletion,na.rm=T),
             Depletion_DBSRA=mean(Depletion_DBSRA,na.rm=T),
             Depletion_CMSY=mean(Depletion_CMSY,na.rm=T),
             Depletion_JABBA=mean(Depletion_JABBA,na.rm=T),
             Depletion_SSS=mean(Depletion_SSS,na.rm=T))
    COMs.weight[[s]]=dummy[nrow(dummy),]%>%
      dplyr::select(r.group,r,K,Init.dep,HRscen,Depletion,
                    Depletion_DBSRA,Depletion_CMSY,Depletion_JABBA,Depletion_SSS)
    rm(dummy)
  }
  COMs.weight=do.call(rbind,COMs.weight)
  rownames(COMs.weight)=NULL
  
  COMs.weight=COMs.weight%>%
    mutate(SSS.error=abs((Depletion_SSS-Depletion)/Depletion),
           JABBA.error=abs((Depletion_JABBA-Depletion)/Depletion),
           CMSY.error=abs((Depletion_CMSY-Depletion)/Depletion),
           DBSRA.error=abs((Depletion_DBSRA-Depletion)/Depletion),
           SSS.weight=1/SSS.error,
           JABBA.weight=1/JABBA.error,
           CMSY.weight=1/CMSY.error,
           DBSRA.weight=1/DBSRA.error,
           Sum.weight=SSS.weight+JABBA.weight+CMSY.weight+DBSRA.weight,
           SSS.weight=SSS.weight/Sum.weight,
           JABBA.weight=JABBA.weight/Sum.weight,
           CMSY.weight=CMSY.weight/Sum.weight,
           DBSRA.weight=DBSRA.weight/Sum.weight)
  
  #overall weight
  Overall.weight=data.frame(Model=c('SSS','JABBA','CMSY','DBSRA'),
                            Weight=c(median(COMs.weight$SSS.weight),
                                     median(COMs.weight$JABBA.weight),
                                     median(COMs.weight$CMSY.weight),
                                     median(COMs.weight$DBSRA.weight)))%>%
    mutate(Tot.w=sum(Weight),
           Weight=Weight/Tot.w)%>%
    dplyr::select(-Tot.w)
  write.csv(Overall.weight,
            paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weight.csv',sep=''),row.names = F) 
  
  #weight by r group
  by.r.group.weight=vector('list',length(r.groups))
  names(by.r.group.weight)=names(r.groups)
  for(x in 1:length(by.r.group.weight))
  {
    xx=names(by.r.group.weight)[x]
    dummy=COMs.weight%>%filter(r.group==xx)
    by.r.group.weight[[x]]=data.frame(Model=c('SSS','JABBA','CMSY','DBSRA'),
                                      Weight=c(median(dummy$SSS.weight),
                                               median(dummy$JABBA.weight),
                                               median(dummy$CMSY.weight),
                                               median(dummy$DBSRA.weight)))%>%
      mutate(Tot.w=sum(Weight),
             Weight=Weight/Tot.w)%>%
      dplyr::select(-Tot.w)%>%
      mutate(r.group=xx)
  }
  write.csv(do.call(rbind,by.r.group.weight),
            paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weight_by.r.group.csv',sep=''),row.names = F) 
  
  #weight by r group and exploitation history
  hr.group=list(c('HarvestRate1','HarvestRate2'),c('HarvestRate3','HarvestRate4'))
  by.r.group_hr.group.weight=rep(hr.group,length(r.groups))
  names(by.r.group_hr.group.weight)=rep(names(r.groups),length(hr.group))
  for(x in 1:length(by.r.group_hr.group.weight))
  {
    xx=names(by.r.group_hr.group.weight)[x]
    dummy=COMs.weight%>%filter(r.group==xx & HRscen%in%by.r.group_hr.group.weight[[x]])
    by.r.group_hr.group.weight[[x]]=data.frame(Model=c('SSS','JABBA','CMSY','DBSRA'),
                                               Weight=c(median(dummy$SSS.weight),
                                                        median(dummy$JABBA.weight),
                                                        median(dummy$CMSY.weight),
                                                        median(dummy$DBSRA.weight)))%>%
      mutate(Tot.w=sum(Weight),
             Weight=Weight/Tot.w)%>%
      dplyr::select(-Tot.w)%>%
      mutate(r.group=xx,
             HRscen=paste(by.r.group_hr.group.weight[[x]],collapse='_'))
  }
  write.csv(do.call(rbind,by.r.group_hr.group.weight),
            paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weight_by.r.group_hr.group.csv',sep=''),row.names = F) 
  
  
}
if(!do.ensemble.simulations)
{
  if(Ensemble.weight=='weighted')
  {
    COM_weight_overall=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weightedl_weight.csv',sep=''))
    COM_weight_by.r.group=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weighted_weight_by.r.group.csv',sep=''))
    COM_weight_by.r.group.h.group=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_weighted_weight_by.r.group_hr.group.csv',sep=''))
    
  }
  if(Ensemble.weight=='equal')
  {
    COM_weight_overall=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_equal_weight.csv',sep=''))
    COM_weight_by.r.group=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_equal_weight_by.r.group.csv',sep=''))
    COM_weight_by.r.group.h.group=read.csv(paste(handl_OneDrive("Analyses/Population dynamics/Ensemble/"),'COMs_equal_weight_by.r.group_hr.group.csv',sep=''))
  }
  
}


#---Generate outputs --------------------------------------
#18.3 Table of scenarios
for(l in 1: length(Lista.sp.outputs))  
{
  dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species)]
  if(length(dis.spicis)>0)
  {
    for(w in 1:length(Catch_only))
    {
      dummy=Catch_only[[w]]$sens.table
      dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
      
      write.csv(do.call(rbind,dummy)%>%relocate(Species),
                paste(Rar.path,paste('Table 2. Catch only_scenarios_',
                                     names(Catch_only)[w],'_',names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
                row.names = F)
    }
  }  
}


#18.4 Table of parameter estimates by species and catch-only assessment method
for(l in 1: length(Lista.sp.outputs))
{
  dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species)]
  if(length(dis.spicis)>0)
  {
    for(w in 1:length(Catch_only))
    {
      dummy=Catch_only[[w]]$estimates
      dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
      dummy=do.call(rbind,dummy)%>%
        rownames_to_column(var = "Species")%>%
        mutate(Species=capitalize(str_extract(Species, "[^.]+")))%>%
        relocate(Species)
      write.csv(dummy,paste(Rar.path,paste('Table 3. Catch only_estimates_',
                                           names(Catch_only)[w],'_',names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
                row.names = F)
    }
  }
}


#18.5. Time series & catch vs MSY by species  
#note: JABBA only outputs 95% CI so only report 95% for all COMs

  #18.5.1 Relative biomass (i.e. Depletion)  
YLAB=XLAB=''
Ref.points=vector('list',N.sp)
names(Ref.points)=Keep.species
for(i in 1:N.sp)
{
  if(Keep.species[i]%in%Catch.only.species)
  {
    print(paste("Catch-only --- Relative biomass plot -----",Keep.species[i]))
    a=fn.plot.timeseries(d=Catch_only,
                         sp=Keep.species[i],
                         Type='Depletion',
                         YLAB='Relative biomass')
    #export graph
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_time_series_relative_biomass.tiff",sep=''),
           width = 6,height = 10,compression = "lzw")
    
    #export current depletion probabilities
    write.csv(a$store.probs%>%
                spread(Model,Probability)%>%
                mutate(Species=Keep.species[i],
                       Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                arrange(Scenario,Range),
              paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_current_depletion.csv",sep=''),
              row.names = F)
    Ref.points[[i]]=a$Ref.points
    
  }
}

  #18.5.2 Fishing mortality
if(do.F.series)
{
  for(i in 1:N.sp)
  {
    if(Keep.species[i]%in%Catch.only.species)
    {
      print(paste("Catch-only --- Fishing mortality plot -----",Keep.species[i]))
      a=fn.plot.timeseries(d=Catch_only,
                           sp=Keep.species[i],
                           Type='F.series',
                           YLAB='Fishing mortality')
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_time_series_fishing_mortality.tiff",sep=''),
             width = 5,height = 10,compression = "lzw")
      
    }
  }
}

  #18.5.3 B over Bmsy
if(do.B.over.Bmsy.series)
{
  for(i in 1:N.sp)
  {
    if(Keep.species[i]%in%Catch.only.species)
    {
      print(paste("Catch-only --- B over Bmsy plot -----",Keep.species[i]))
      a=fn.plot.timeseries(d=Catch_only,
                           sp=Keep.species[i],
                           Type='B.Bmsy',
                           YLAB='B/Bmsy')
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_time_series_B_Bmsy.tiff",sep=''),
             width = 6,height = 10,compression = "lzw")
    }
  }
}

  #18.5.4 F over Fmsy
if(do.F.over.Fmsy.series)
{
  for(i in 1:N.sp)
  {
    if(Keep.species[i]%in%Catch.only.species)
    {
      print(paste("Catch-only --- F over Fmsy plot -----",Keep.species[i]))
      a=fn.plot.timeseries(d=Catch_only,
                           sp=Keep.species[i],
                           Type='F.Fmsy',
                           YLAB='F/Fmsy')
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_time_series_F_Fmsy.tiff",sep=''),
             width = 5,height = 10,compression = "lzw")
    }
  }
}

  #18.5.5 Catch vs MSY
Store.catch.MSY=vector('list',N.sp)
for(i in 1:N.sp)
{
  if(Keep.species[i]%in%Catch.only.species)
  {
    print(paste("Catch-only --- Catch vs MSY -----",Keep.species[i]))  
    Store.catch.MSY[[i]]=fn.plot.catch.vs.MSY(d=Catch_only,
                                               sp=Keep.species[i],
                                               YLAB='Catch (tonnes)')
    
    #export graph
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_catch_vs_MSY.tiff",sep=''),
           width = 6,height = 10,compression = "lzw")
    
  }
}

#18.6. Get ensemble. Sample each accepted COM proportionally to weights  
#note: first sample equally within each model type (done in'#output quantities for ensemble')
#      second sample models proportionally to weights (done here)
Mod.AV_depletion=vector('list',length(List.sp))
names(Mod.AV_depletion)=names(List.sp)
Mod.AV_B.Bmsy=Mod.AV_F.Fmsy=Mod.AV_MSY=Mod.AV_Probs=Mod.AV_Probs_B.Bmsy=Ref.points.Mod.AV=Mod.AV_depletion
for(i in 1:N.sp)  
{
  if(Keep.species[i]%in%Catch.only.species)
  {
    print(paste("Model average for--------",Keep.species[i]))
    dis.r.range=do.call(rbind,r.groups)%>%   
      data.frame%>%
      mutate(Get=between(r.list[i],.[[1]] , .[[2]]))%>%
      filter(Get==TRUE)
    Wei=COM_weight_by.r.group%>%
      filter(r.group==row.names(dis.r.range))%>%
      dplyr::select(-r.group)
    if(!'SSS'%in%Wei$Model)
    {
      Wei.SSS=Wei%>%filter(Model=='DBSRA')%>%mutate(Model='SSS')
      if('SSS'%in%names(Catch_only)) Wei=rbind(Wei,Wei.SSS)
    }
    Wei=Wei%>%filter(Model%in%names(Catch_only)) 
    
    #Reference points  
    Avg.Ref.point=Ref.points[[i]]
    Avg.Ref.point=unlist(Avg.Ref.point,recursive=FALSE)
    for(av in 1:length(Avg.Ref.point))
    {
      Avg.Ref.point[[av]]=Avg.Ref.point[[av]]%>%
                            mutate(Model=names(Avg.Ref.point)[av],
                                   Model=gsub('[0-9]+', '', Model))
    }
      
    Avg.Ref.point=do.call(rbind,Avg.Ref.point)%>%
                    left_join(Wei,by='Model')%>%
                    group_by(Rf.pt)%>%
                    summarise(Value=weighted.mean(Value,Weight,na.rm=T))%>%
                    mutate(Rf.pt=factor(Rf.pt,levels=c('Target','Threshold','Limit')))%>%
                    arrange(Rf.pt)
    Ref.points.Mod.AV[[i]]=Avg.Ref.point
    
    #depletion
    dumi=mod.average(dd=map(Catch_only, ~.x$ensemble[[i]]$Depletion),
                     Weights=Wei,
                     Ref.pnt=Avg.Ref.point)
    Mod.AV_depletion[[i]]=dumi$Trajectories
    dumi$Probs$Species=capitalize(Keep.species[i])
    dumi$Probs.future$Species=capitalize(Keep.species[i])
    Mod.AV_Probs[[i]]=dumi[c("Probs" ,"Probs.future")] 
    
    #B.Bmsy
    dumi=mod.average(dd=map(Catch_only, ~.x$ensemble[[i]]$B.Bmsy),
                     Weights=Wei,
                     Ref.pnt=data.frame(Rf.pt=c('Target','Threshold','Limit'),
                                        Value=c(Tar.prop.bmsny,1,Lim.prop.bmsy)))
    Mod.AV_B.Bmsy[[i]]=dumi$Trajectories
    dumi$Probs$Species=capitalize(Keep.species[i])
    dumi$Probs.future$Species=capitalize(Keep.species[i])
    Mod.AV_Probs_B.Bmsy[[i]]=dumi[c("Probs" ,"Probs.future")]
    
    #F.Fmsy
    dumi=mod.average(dd=map(Catch_only, ~.x$ensemble[[i]]$F.Fmsy),   
                     Weights=Wei)
    Mod.AV_F.Fmsy[[i]]=dumi$Trajectories
    
    #MSY
    Mod.AV_MSY[[i]]=mod.average.scalar(dd=map(Catch_only, ~.x$ensemble[[i]]$MSY),   
                                       Weights=Wei)  

    
    rm(Wei)
    
  }
}
Mod.AV=list(rel.biom=Mod.AV_depletion,Probs=Mod.AV_Probs,
            B.Bmsy=Mod.AV_B.Bmsy,Probs_B.Bmsy=Mod.AV_Probs_B.Bmsy,
            F.Fmsy=Mod.AV_F.Fmsy,MSY=Mod.AV_MSY)
rm(Mod.AV_depletion,Mod.AV_Probs,Mod.AV_B.Bmsy,Mod.AV_F.Fmsy,Mod.AV_MSY,Mod.AV_Probs_B.Bmsy)


#18.7 Display Scenarios for combined species: each COM model as an Appendix 

  #18.7.1 Relative biomass (i.e. Depletion) 
    #figure  (Scenario 1 only)
for(l in 1:length(Lista.sp.outputs))
{
  if(length(Lista.sp.outputs[[l]])>8) InMar=.65 else InMar=.5
  dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
  if(length(dis.spicis)>0)
  {
    print(paste("RAR ---- relative.biomass_catch.only by model-----",names(Lista.sp.outputs)[l]))
    fn.plot.timeseries_combined_Appendix(this.sp=dis.spicis,   
                                         d=Catch_only,
                                         YLAB="Relative biomass",
                                         NM=names(Lista.sp.outputs)[l],
                                         Type="Depletion",
                                         InnerMargin=InMar)
  }
}
    #table (all scenarios)
for(l in 1:length(Lista.sp.outputs))
{
  dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
  if(length(dis.spicis)>0)
  {
    dummy.mod=vector('list',length(Catch_only))
    for(m in 1:length(Catch_only))
    {
      str.prob=Catch_only[[m]]$probs.rel.biom
      str.prob=str.prob[match(Lista.sp.outputs[[l]],names(str.prob))]
      str.prob=compact(str.prob)
      dummy=vector('list',length =length(str.prob))
      for(d in 1:length(dummy))
      {
        ddmmi=fn.get.stuff.from.list(str.prob[[d]],'probs')
        dummy[[d]]=do.call(rbind,ddmmi)%>%
          mutate(Species=capitalize(names(str.prob)[d]))
      }
      dummy.mod[[m]]=do.call(rbind,dummy)%>%
        mutate(Model=names(Catch_only)[m])
    }
    write.csv(do.call(rbind,dummy.mod)%>%
                mutate(Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                spread(Species,Probability)%>%
                arrange(Range),
              paste(Rar.path,'/Table 4. Catch only_current.depletion_',names(Lista.sp.outputs)[l],'_Appendix','.csv',sep=''),
              row.names=F)
    rm(dummy.mod)
  }
}

  #18.7.2 B over Bmsy        
if(do.B.over.Bmsy.series)
{
  #figure  (Scenario 1 only)
  for(l in 1:length(Lista.sp.outputs))
  {
    if(length(Lista.sp.outputs[[l]])>8) InMar=.65 else InMar=.5
    dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
    if(length(dis.spicis)>0)
    {
      print(paste("RAR ---- B.over.Bmsy_catch.only by model-----",names(Lista.sp.outputs)[l]))
      fn.plot.timeseries_combined_Appendix(this.sp=dis.spicis,   
                                           d=Catch_only,
                                           YLAB="B/Bmsy",
                                           NM=names(Lista.sp.outputs)[l],
                                           Type="B.Bmsy",
                                           InnerMargin=InMar)
    }
  }
  
  #table (all scenarios)
  for(l in 1:length(Lista.sp.outputs))
  {
    dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
    if(length(dis.spicis)>0)
    {
      dummy.mod=vector('list',length(Catch_only))
      for(m in 1:length(Catch_only))
      {
        str.prob=Catch_only[[m]]$probs.B.Bmsy
        str.prob=str.prob[match(Lista.sp.outputs[[l]],names(str.prob))]
        str.prob=compact(str.prob)
        dummy=vector('list',length =length(str.prob))
        for(d in 1:length(dummy))
        {
          ddmmi=fn.get.stuff.from.list(str.prob[[d]],'probs')
          dummy[[d]]=do.call(rbind,ddmmi)%>%
            mutate(Species=capitalize(names(str.prob)[d]))
        }
        dummy.mod[[m]]=do.call(rbind,dummy)%>%
          mutate(Model=names(Catch_only)[m])
      }
      write.csv(do.call(rbind,dummy.mod)%>%
                  mutate(Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                  spread(Species,Probability)%>%
                  arrange(Range),
                paste(Rar.path,'/Table 4. Catch only_current.B.over.Bmsy_',names(Lista.sp.outputs)[l],'_Appendix','.csv',sep=''),
                row.names=F)
      rm(dummy.mod)
    }
  }
}

  #18.7.3 Fishing mortality 
#figure (Scenario 1 only)
if(do.F.series)
{
  for(l in 1:length(Lista.sp.outputs))
  {
    if(length(Lista.sp.outputs[[l]])>8) InMar=.65 else InMar=.5
    dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
    if(length(dis.spicis)>0)
    {
      print(paste("RAR ---- Fishing mortality_catch.only by model-----",names(Lista.sp.outputs)[l]))
      fn.plot.timeseries_combined_Appendix(this.sp=dis.spicis,   
                                           d=Catch_only,
                                           YLAB="Fishing mortality (y-1)",
                                           NM=names(Lista.sp.outputs)[l],
                                           Type="F.series",
                                           InnerMargin=InMar)
    }
  }
}

  #18.7.4 F over Fmsy   
  #figure  (Scenario 1 only)
if(do.F.over.Fmsy.series)
{
  for(l in 1:length(Lista.sp.outputs))
  {
    if(length(Lista.sp.outputs[[l]])>8) InMar=.65 else InMar=.5
    dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
    if(length(dis.spicis)>0)
    {
      print(paste("RAR ---- F.over.Fmsy_catch.only by model-----",names(Lista.sp.outputs)[l]))
      fn.plot.timeseries_combined_Appendix(this.sp=dis.spicis,   
                                           d=Catch_only,
                                           YLAB="F/Fmsy",
                                           NM=names(Lista.sp.outputs)[l],
                                           Type="F.Fmsy",
                                           InnerMargin=InMar)
    }
  }
}

  #18.7.5 Catch vs MSY       
All.catch.MSY=do.call(rbind,compact(Store.catch.MSY))
All.catch.MSY.species=unique(All.catch.MSY$Species)
#figure & table (Scenario 1 only)
for(l in 1:length(Lista.sp.outputs))
{
  if(length(Lista.sp.outputs[[l]])>8) InMar=.65 else InMar=.5
  all.dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%All.catch.MSY.species)]
  dis.spicis=subset(all.dis.spicis,all.dis.spicis%in%Catch.only.species_display) 
  print(paste("RAR ---- Catch vs MSY_catch.only by model-----",names(Lista.sp.outputs)[l]))
  fn.plot.catch.vs.MSY_combined_Appendix(all.this.sp=all.dis.spicis,
                                         this.sp=dis.spicis,
                                         d=All.catch.MSY,
                                         NM=names(Lista.sp.outputs)[l],
                                         YLAB='Catch (tonnes)')
}


#18.8 Display COM model-average for combined species   

# 18.8.1 Relative biomass (i.e. Depletion) 
    #figure
for(l in 1:length(Lista.sp.outputs))
{
  if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
  dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
  if(length(dis.spicis)>0)
  {
    print(paste("RAR ---- relative.biomass_catch.only Ensemble-----",names(Lista.sp.outputs)[l]))
    fn.plot.timeseries_combined(this.sp=dis.spicis,
                                d=Mod.AV,
                                YLAB="Relative biomass",
                                Type="Depletion",
                                InnerMargin=InMar,
                                RefPoint=Ref.points.Mod.AV,
                                Kach=Catch_only$DBSRA$rel.biom)  
    ggsave(paste(Rar.path,'/Relative.biomass_catch.only_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = 13,height = 11,compression = "lzw")
  }
  
  #Display species with only catch separately
  if(names(Lista.sp.outputs)[l]=="Other.sp")
  {
    if(length(Catch.only.species_only.ktch)>8) InMar=1.25 else InMar=.5
    print(paste("RAR ---- relative.biomass_catch.only Ensemble-----","Other.sp_species_only.ktch"))
    fn.plot.timeseries_combined(this.sp=Catch.only.species_only.ktch,
                                d=Mod.AV,
                                YLAB="Relative biomass",
                                Type="Depletion",
                                InnerMargin=InMar,
                                RefPoint=Ref.points.Mod.AV,
                                Kach=Catch_only$DBSRA$rel.biom)
    ggsave(paste(Rar.path,'/Relative.biomass_catch.only_Other.sp_species_only.ktch.tiff',sep=''),
           width = 8,height = 10,compression = "lzw")
  }
}
    #table
for(l in 1:length(Lista.sp.outputs))
{
  dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
  if(length(dis.spicis)>0)
  {
    Probs=do.call(rbind,fn.get.stuff.from.list(Mod.AV$Probs[dis.spicis],'Probs'))%>%
      mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
      spread(Species,Probability)%>%
      arrange(Range)
    Probs.future=do.call(rbind,fn.get.stuff.from.list(Mod.AV$Probs[dis.spicis],'Probs.future'))%>%
      mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
      spread(Species,Probability)%>%
      arrange(Range)
    dd=rbind(Probs,Probs.future)%>%
      mutate(dummy=paste(Range, finyear))%>%
      distinct(dummy,.keep_all = T)%>%
      dplyr::select(-dummy)
    write.csv(dd,
              paste(Rar.path,'/Table 4. Catch only_current.depletion_',names(Lista.sp.outputs)[l],'.csv',sep=''),
              row.names=F)
  }
}

  #18.8.2 B over Bmsy        
if(do.B.over.Bmsy.series) 
{
  #figure
  for(l in 1:length(Lista.sp.outputs))
  {
    if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
    dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
    if(length(dis.spicis)>0)
    {
      print(paste("RAR ---- B.over.Bmsy_catch.only Ensemble-----",names(Lista.sp.outputs)[l]))
      fn.plot.timeseries_combined(this.sp=dis.spicis,
                                  d=Mod.AV,
                                  YLAB="B/Bmsy",
                                  Type="B.Bmsy",
                                  InnerMargin=InMar,
                                  RefPoint=NULL,
                                  Kach=Catch_only$DBSRA$rel.biom)  
      ggsave(paste(Rar.path,'/B.over.Bmsy_catch.only_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
             width = 13,height = 11,compression = "lzw")
    }
    
    #Display species with only catch separately
    if(names(Lista.sp.outputs)[l]=="Other.sp")
    {
      if(length(Catch.only.species_only.ktch)>8) InMar=1.25 else InMar=.5
      print(paste("RAR ---- B.over.Bmsy_catch.only Ensemble-----","Other.sp_species_only.ktch"))
      fn.plot.timeseries_combined(this.sp=Catch.only.species_only.ktch,
                                  d=Mod.AV,
                                  YLAB="B/Bmsy",
                                  Type="B.Bmsy",
                                  InnerMargin=InMar,
                                  RefPoint=NULL,
                                  Kach=Catch_only$DBSRA$rel.biom)
      ggsave(paste(Rar.path,'/B.over.Bmsy_catch.only_Other.sp_species_only.ktch.tiff',sep=''),
             width = 8,height = 10,compression = "lzw")
    }
  }
  #table
  for(l in 1:length(Lista.sp.outputs))
  {
    dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
    if(length(dis.spicis)>0)
    {
      Probs=do.call(rbind,fn.get.stuff.from.list(Mod.AV$Probs_B.Bmsy[dis.spicis],'Probs'))%>%
        mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
        spread(Species,Probability)%>%
        arrange(Range)
      Probs.future=do.call(rbind,fn.get.stuff.from.list(Mod.AV$Probs_B.Bmsy[dis.spicis],'Probs.future'))%>%
        mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
        spread(Species,Probability)%>%
        arrange(Range)
      dd=rbind(Probs,Probs.future)%>%
        mutate(dummy=paste(Range, finyear))%>%
        distinct(dummy,.keep_all = T)%>%
        dplyr::select(-dummy)
      write.csv(dd,
                paste(Rar.path,'/Table 4. Catch only_current.B.over.Bmsy_',names(Lista.sp.outputs)[l],'.csv',sep=''),
                row.names=F)
    }
  }
}


  #18.8.3 F over Fmsy   
if(do.F.over.Fmsy.series)
{
  for(l in 1:length(Lista.sp.outputs))
  {
    if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
    dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%Catch.only.species_display)]
    if(length(dis.spicis)>0)
    {
      print(paste("RAR ---- F.over.Fmsy_catch.only Ensemble-----",names(Lista.sp.outputs)[l]))
      fn.plot.timeseries_combined(this.sp=dis.spicis,
                                  d=Mod.AV,
                                  YLAB="F/Fmsy",
                                  Type="F.Fmsy",
                                  InnerMargin=InMar,
                                  RefPoint=NULL,
                                  Kach=Catch_only$DBSRA$rel.biom)  
      ggsave(paste(Rar.path,'/F.over.Fmsy_catch.only_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
             width = 13,height = 11,compression = "lzw")
    }
    
    
  }
}

# 18.8.4 Catch vs MSY  
  #figure & table
for(l in 1:length(Lista.sp.outputs))
{
  all.dis.spicis=Lista.sp.outputs[[l]][which(Lista.sp.outputs[[l]]%in%All.catch.MSY.species)]
  dis.spicis=subset(all.dis.spicis,all.dis.spicis%in%Catch.only.species_display)
  print(paste("RAR ---- Catch vs MSY_catch.only Ensemble-----",names(Lista.sp.outputs)[l]))
  fn.plot.catch.vs.MSY_combined(all.this.sp=all.dis.spicis,
                                this.sp=dis.spicis,
                                d=All.catch.MSY,
                                NM=names(Lista.sp.outputs)[l],
                                YLAB='Catch (tonnes)',
                                MSY=compact(Mod.AV$MSY))
}


#18.9. Kobe plots (Scenario 1)   

  #18.9.1 by species
store.kobes=vector('list',N.sp)
names(store.kobes)=Keep.species
for(i in 1:N.sp)
{
  if(Keep.species[i]%in%Catch.only.species)
  {
    print(paste("Catch-only --- Kobe plot S1-----",Keep.species[i]))
    store.kobes[[i]]=fn.get.Kobe.plot_appendix(d=Catch_only,
                                               sp=Keep.species[i])
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/Catch_only_Kobe_plot.tiff",sep=''),
           width = 9,height = 14, dpi = 300,compression = "lzw")
  }

}

  #18.9.2 Display combined species: each COM model as an Appendix
do.this=FALSE
if(do.this)
{
  for(l in 1:length(Lista.sp.outputs))
  {
    if(length(Lista.sp.outputs[[l]])>8)
    {
      a=store.kobes[which(names(store.kobes)%in%Lista.sp.outputs[[l]])]
      for(m in 1:length(Catch_only))
      {
        figs=vector('list',length(a))
        for(f in 1:length(figs))
        {
          pp=a[[f]][[match(names(Catch_only)[m],names(a[[f]]))]]
          pp$labels$title=capitalize(names(a)[f])
          # figure <- ggarrange(plotlist=pp,
          #                     ncol=3,common.legend = FALSE)+
          #   theme(plot.margin = margin(1,0,0,0, "cm"))
          # figure=annotate_figure(figure,
          #                        fig.lab=capitalize(names(a)[f]),
          #                        fig.lab.pos='top.left',
          #                        fig.lab.size=22)
          figs[[f]]<-pp
          rm(pp)
        }
        figure=ggarrange(plotlist=figs)+
          theme(plot.margin = margin(0,.2,0,0, "cm"))
        annotate_figure(figure,
                        top=text_grob(names(Catch_only)[m], size=22),
                        bottom = text_grob(expression(B/~B[MSY]), size=22),
                        left = text_grob(expression(F/~F[MSY]), rot = 90,size=22))
        WID=13
        if(names(Catch_only)[m]=="JABBA") WID=16
        ggsave(paste(Rar.path,'/Kobe_plot_catch_only_',names(Lista.sp.outputs)[l],
                     '_',names(Catch_only)[m],'_Appendix.tiff',sep=''),
               width = WID,height = 12,compression = "lzw")
      }
    }else
    {
      a=store.kobes[which(names(store.kobes)%in%Lista.sp.outputs[[l]])]
      figs=vector('list',length(a))
      for(f in 1:length(figs))
      {
        figure <- ggarrange(plotlist=a[[f]],
                            ncol=3,common.legend = FALSE)+
          theme(plot.margin = margin(1,0,0,0, "cm"))
        figure=annotate_figure(figure,
                               fig.lab=capitalize(names(a)[f]),
                               fig.lab.pos='top.left',
                               fig.lab.size=22)
        figs[[f]]<-figure
      }
      figure=ggarrange(plotlist=figs,nrow=length(figs))
      annotate_figure(figure,
                      bottom = text_grob(expression(B/~B[MSY]), size=22),
                      left = text_grob(expression(F/~F[MSY]), rot = 90,size=22))
      ggsave(paste(Rar.path,'/Kobe_plot_catch_only_',names(Lista.sp.outputs)[l],'_Appendix.tiff',sep=''),
             width = 12,height = 12,compression = "lzw")
    }
  }
}

  #18.9.3 Display COM model-average combined species   
for(l in 1:length(Lista.sp.outputs))
{
  this.sp=Lista.sp.outputs[[l]]
  dis.spicis=this.sp[which(this.sp%in%Catch.only.species_display)]
  if(length(dis.spicis)>0)
  {
    print(paste("RAR ---- Catch-only --- Kobe plot Ensemble-----",names(Lista.sp.outputs)[l]))
    DIMS=n2mfrow(length(dis.spicis))
    #NKOL=DIMS[2]
    #NRW=DIMS[1]
    NKOL=round(length(dis.spicis)/5)
    if(names(Lista.sp.outputs)[l]=="Indicator.sp") NKOL=2
    NRW=length(dis.spicis)/NKOL
    if(NKOL%in%3:4) WIZ=14
    if(NKOL==2) WIZ=11
    if(NKOL==1) WIZ=6
    fn.get.Kobe.plot(this.sp=dis.spicis,
                     d=Mod.AV,
                     NKOL,
                     NRW)
    ggsave(paste(Rar.path,'/Kobe_plot_catch_only_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = WIZ,height = 12,compression = "lzw")
  }
}


#18.10 Percentage acceptance rate for Scenario 1 (not applicable to JABBA)   
Per.accepted=vector('list',length(List.sp))
names(Per.accepted)=names(List.sp)
  #18.10.1. Get for each species
for(i in 1: N.sp)
{
  if(Keep.species[i]%in%Catch.only.species)
  {
    print(paste("Catch-only --- Get acceptance rate S1-----",Keep.species[i]))
    moDs=names(Catch_only)[match(c("DBSRA","CMSY","SSS"),names(Catch_only))]  
    dummy=vector('list',length(moDs))
    names(dummy)=moDs
    for(m in 1:length(moDs))
    {
      dummy[[m]]=Catch_only[[m]]$accept.rate[[i]]%>%
        filter(Scenario=='S1')%>%
        mutate(Model=names(dummy)[m],
               Species=Keep.species[i])
    }
    Per.accepted[[i]]=do.call(rbind,dummy)%>%
      relocate(Species,Model,Scenario)
    rm(dummy)
  }

}

  #18.10.2. Export by Lista.sp.outputs 
for(l in 1:length(Lista.sp.outputs))
{
  this.sp=Lista.sp.outputs[[l]]
  dis.spicis=this.sp[which(this.sp%in%Catch.only.species)]
  if(length(dis.spicis)>0)
  {
    write.csv(do.call(rbind,Per.accepted[dis.spicis])%>%
                data.frame%>%
                mutate(Species=capitalize(Species))%>%
                spread(Model,Acceptance),
              paste(Rar.path,'/Table 5. Catch.only_per.acceptance.rate_',
                    names(Lista.sp.outputs)[l],'.csv',sep=''),
              row.names = F)
  }
}


#18.11 MSY estimates by Lista.sp.outputs (Scenario 1)   
get.MSY.estimates=FALSE  #redundant, MSY output in #18.4
if(get.MSY.estimates)
{
  for(l in 1: length(Lista.sp.outputs))
  {
    dummy=vector('list',length(Catch_only))
    for(w in 1:length(dummy))
    {
      a=Catch_only[[w]]$estimates[match(Lista.sp.outputs[[l]],names(Catch_only[[w]]$estimates))]
      dummy[[w]]=do.call(rbind,a)%>%
        filter(Scenario=='S1' & Parameter=='MSY')
    }
    write.csv(do.call(rbind,dummy),
              paste(Rar.path,'/Table 6. MSY_estimates_catch.only_',
                    names(Lista.sp.outputs)[l],'.csv',sep=''),
              row.names = F)
  }
}


#18.12 Store Consequence and likelihood for WoE 
  #Ensemble of models for S1 only
Probs=do.call(rbind,fn.get.stuff.from.list(Mod.AV$Probs,'Probs'))%>%
  mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
  arrange(Range)
Probs.future=do.call(rbind,fn.get.stuff.from.list(Mod.AV$Probs,'Probs.future'))%>%
  mutate(Range=factor(Range,levels=c('>tar','thr.tar','lim.thr','<lim')))%>%
  arrange(Range)
Store.cons.Like_COM=vector('list',N.sp)
names(Store.cons.Like_COM)=Keep.species
for(i in 1:N.sp)
{
  Store.cons.Like_COM[[i]]=rbind(Probs%>%
                                   filter(Species==capitalize(Keep.species[i]))%>%
                                   dplyr::select(-Species),
                                  Probs.future%>%
                                    filter(Species==capitalize(Keep.species[i]))%>%
                                    dplyr::select(-Species))
}
rm(Probs,Probs.future)