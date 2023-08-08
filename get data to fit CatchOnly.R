# Output files for SS-DL-tool-------------------------------------------------------------------------
THIS=match(Catch.only.species,Keep.species)
for(i in THIS)
{
  Neim=Keep.species[i]
  print(Neim)
  WDD=paste(handl_OneDrive('SS3/SS-DL-tool-master/My data files/'),capitalize(Neim),sep='')
  if (!dir.exists(WDD)) dir.create(WDD)
  
  #1. Catch
  ktch=ktch.combined%>%
    filter(Name==Neim)
  all.years=seq(min(ktch$finyear),max(ktch$finyear))
  misn.yr=all.years[which(!all.years%in%ktch$finyear)]
  if(length(misn.yr)>0)
  {
    ktch=rbind(ktch,ktch[length(misn.yr),]%>%mutate(finyear=misn.yr,Tonnes=0))%>%arrange(finyear)
  }
  ktch=ktch%>%
    rename(Year=finyear,
           Total=Tonnes)%>%
    ungroup()%>%
    dplyr::select(Year,Total)
  
   write.csv(ktch,paste(WDD,'Catches.csv',sep='/'),row.names = F)
  
  #Life history
  Life.history=List.sp[[i]]
  Life.history$Fecundity=ceiling(mean(Life.history$Fecundity))
  Life.history$Max.age.F=ceiling(mean(Life.history$Max.age.F))
  Life.history$Breed.cycle=mean(Life.history$Breed.cycle)
  LH=data.frame(Parameter=c('M','h','Linf','k','TL.50','TL.95','AwT','BwT','Max.age','Final.dpl'),
                Value=with(Life.history,c(Sens.test$SS$Mmean[1],Sens.test$SS$Steepness[1],
                                          Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL,Growth.F$k,
                                          TL.50.mat,TL.95.mat,AwT,BwT,Max.age.F,Sens.test$SS3$Final.dpl[1])),
                CV=with(Life.history,c(Sens.test$DBSRA$Msd[1],Sens.test$SS$Steepness.sd[1],
                                       Growth.F$FL_inf.sd,Growth.F$k.sd,NA,NA,NA,NA,NA,NA)))
  
  write.csv(LH,paste(WDD,'Life.history.csv',sep='/'),row.names = F)
  
  
  
}

# Fit catch only-------------------------------------------------------------------------

tic()
for(g in 1:length(Catch.only.species))
{
  i=match(Catch.only.species[g],Keep.species)
  
  out.here=paste("C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/Population dynamics/Exploration_ctch only/",Keep.species[i],"_",sep='')  
  
  #1. DBSRA assessment (Dick and MAcCall (2011))
  dothis=FALSE
  if(dothis)     #parallelised: 0.004 secs per iteration-species-scenario (otherwise 0.025 secs)
  {
    dummy.store=vector('list',1)     
    names(dummy.store)=Keep.species[i]
    dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
      dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.accept.rate=
      dummy.store.ensemble=dummy.store
    if(names(dummy.store)%in%Catch.only.species)    
    {
      #Catch
      ktch=ktch.combined%>%
        filter(Name==names(dummy.store))
      all.years=seq(min(ktch$finyear),max(ktch$finyear))
      misn.yr=all.years[which(!all.years%in%ktch$finyear)]
      if(length(misn.yr)>0)
      {
        ktch=rbind(ktch,ktch[length(misn.yr),]%>%mutate(finyear=misn.yr,Tonnes=0))%>%arrange(finyear)
      }
      
      this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(List.sp[[i]]$Name),"/",AssessYr,"/DBSRA",sep='')
      if(!dir.exists(this.wd))dir.create(this.wd)
      
      ktch%>%
        ggplot(aes(finyear,Tonnes))+
        geom_line()+geom_point()+ggtitle(Keep.species[i])
      ggsave(paste(out.here,"total catch.tiff",sep=''),width = 6,height = 6, dpi = 300, compression = "lzw")
      
      
      #Scenarios
      Scens=List.sp[[i]]$Sens.test$DBSRA%>%
        mutate(Species=capitalize(names(dummy.store)))
      Store.sens=vector('list',nrow(Scens))
      names(Store.sens)=Scens$Scenario
      this.wd1=this.wd
      
      Out.Scens=Scens%>%
        mutate(M.dist=NA,M.CV=NA,fmsym.dist=NA,fmsym.CV=NA,bmsyk.dist=NA,bmsyk.mean=NA,bmsyk.CV=NA,
               b1k.dist=NA,b1k.low=NA,b1k.up=NA,btk.dist=NA,btk.low=NA,btk.up=NA)
      Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
        Out.B.Bmsy=Out.F.Fmsy=Out.accept.rate=Out.ensemble=vector('list',length(Store.sens))
      names(Out.ensemble)=names(Store.sens)
      
      cl <- makeCluster(detectCores()-1)
      registerDoSNOW(cl)
      for(s in 1:length(Store.sens))
      {
        print(paste("___________","DBSRA Scenario",Scens$Scenario[s],"___________",names(dummy.store)))
        
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
        if(names(dummy.store)=="milk shark")
        {
          Depletion.year=1991  #set to end of catch period to allow convergence
          Btklow=0.2
        }
        
        #Run model
        OUT=apply.DBSRA_tweeked(year=ktch$finyear,
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
                                outfile="Appendix_fit")
        
        #Acceptance rate
        Accept.tab=table(OUT$output$Values$ll)
        Accept.rate_DBSRA=data.frame(
          Acceptance=round(100*Accept.tab[2]/sum(Accept.tab),2),
          Scenario=Scens$Scenario[s])
        
        #Store trajectories 
        dummy=fn.ktch.only.get.timeseries(d=OUT,
                                          mods="DBSRA",
                                          Type='Depletion',
                                          scen=Scens$Scenario[s],
                                          Katch=ktch$Tonnes)
        Rel.biomas_DBSRA=dummy$Dat
        
        #Display Priors vs Posteriors for base case scenario (S1)
        par.list=c('m','fmsym','b1k','btk','bmsyk')
        dummy.prior=OUT$input
        names(dummy.prior)=tolower(names(dummy.prior))
        dummy.post=OUT$output$Values%>%filter(ll==1) #accepted runs
        colnames(dummy.post)=tolower(colnames(dummy.post))
        out=vector('list',length(par.list))
        for(p in 1:length(par.list))
        {
          out[[p]]=rbind(data.frame(Distribuion="Prior",
                                    Value=fn.prior(d=dummy.prior[[par.list[p]]])),
                         data.frame(Distribuion="Posterior",
                                    Value=dummy.post[[par.list[p]]]))%>%
            mutate(Parameter=par.list[p])
        }
        rm(dummy.prior,dummy.post)
      }
      stopCluster(cl)
      
    }
    fn.show.density(d=do.call(rbind,out),NCOL=2)
    ggsave(paste(out.here,"Prior.and.posterior_DBSRA.tiff",sep=''),width = 12,height = 14, dpi = 300, compression = "lzw")
  }
  
  
  #2. CMSY 
  dothis=FALSE
  if(dothis)     #takes 0.008 secs per iteration-species-scenario
  {
    dummy.store=vector('list',1)     
    names(dummy.store)=Keep.species[i] 
    dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
      dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.accept.rate=
      dummy.store.ensemble=dummy.store
    for(i in i)  
    {
      if(names(dummy.store)%in%Catch.only.species) 
      {
        this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                      capitalize(List.sp[[i]]$Name),"/",AssessYr,"/CMSY",sep='')
        if(!dir.exists(this.wd))dir.create(this.wd)
        
        #Catch
        ktch=ktch.combined%>%
          filter(Name==names(dummy.store))
        all.years=seq(min(ktch$finyear),max(ktch$finyear))
        misn.yr=all.years[which(!all.years%in%ktch$finyear)]
        if(length(misn.yr)>0)
        {
          ktch=rbind(ktch,ktch[length(misn.yr),]%>%mutate(finyear=misn.yr,Tonnes=0))%>%arrange(finyear)
        }
        year=ktch$finyear
        catch=ktch$Tonnes
        
        #Scenarios
        Scens=List.sp[[i]]$Sens.test$CMSY%>%
          mutate(Species=capitalize(names(dummy.store)))
        Store.sens=vector('list',nrow(Scens))
        names(Store.sens)=Scens$Scenario
        this.wd1=this.wd
        
        Out.Scens=Scens%>%   
          mutate(Bo.low=NA,Bo.hi=NA,Bf.low=NA,Bf.hi=NA,r.low=NA,r.up=NA)
        Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
          Out.B.Bmsy=Out.F.Fmsy=Out.accept.rate=Out.ensemble=vector('list',length(Store.sens))
        names(Out.ensemble)=names(Store.sens)
        for(s in 1:length(Store.sens))
        {
          print(paste("___________","CMSY Scenario",Scens$Scenario[s],"___________",names(dummy.store)))
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
            if(names(dummy.store)%in%c("milk shark","great hammerhead")) r.range[1]=0.1
            if(names(dummy.store)=="zebra shark") r.range=c(0.07,Scens$r[s]) 
            if(names(dummy.store)=="green sawfish") r.range[1]=0.09
            if(names(dummy.store)=="narrow sawfish") r.range=c(0.2,Scens$r[s])  
            if(names(dummy.store)=="snaggletooth")  r.range=c(0.2,Scens$r[s])  
            if(names(dummy.store)=="weasel shark")  r.range=c(0.2,Scens$r[s])
            if(names(dummy.store)=="pigeye shark") r.range[2]=0.15
            
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
                                     RES=RES)
          #acceptance rate
          Accept.rate_CMSY=data.frame(
            Acceptance=Store.sens[[s]]$output$acceptance.rate,
            Scenario=Scens$Scenario[s])
          
          #Store trajectories
          dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                            mods="CMSY",
                                            Type='Depletion',
                                            scen=Scens$Scenario[s],
                                            Katch=catch)
          Rel.biomas_CMSY=dummy$Dat
          
          #Display Priors vs Posteriors for base case scenario (S1)
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
          
          
          
        }
      }
    }#end i
    fn.show.density(d=do.call(rbind,out),NCOL=2)
    
    CMSY.K.estim=mean(out[[2]]%>%filter(Parameter=='K' &Distribuion=='Posterior')%>%pull(Value))    
    
    ggsave(paste(out.here,"Prior.and.posterior_CMSY.tiff",sep=''),width = 12,height = 14, dpi = 300, compression = "lzw")
  }
  
  
  #4. JABBA - catch only (Winker et al 2018)   
  dothis=FALSE
  if(dothis)     #takes 0.002 secs per iteration-species-scenario
  {
    dummy.store=vector('list',1)     
    names(dummy.store)=Keep.species[i]
    dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
      dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.Kobe.probs=
      dummy.store.ensemble=dummy.store
    for(i in i)  
    {
      if(names(dummy.store)%in%Catch.only.species) 
      {
        this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                      capitalize(List.sp[[i]]$Name),"/",AssessYr,"/JABBA",sep='')
        if(!dir.exists(this.wd))dir.create(this.wd)
        
        #Catch
        ktch=ktch.combined%>%
          filter(Name==names(dummy.store))%>%
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
        if(names(dummy.store)=='milk shark') ktch$Total[1:10]=1 #convergence issues with first 10 years for milk shark, but catch is 0 anyways
        
        #Scenarios
        Scens=List.sp[[i]]$Sens.test$JABBA_catch.only%>%
          mutate(Species=capitalize(names(dummy.store)))
        Store.sens=vector('list',nrow(Scens))
        names(Store.sens)=Scens$Scenario
        this.wd1=this.wd
        
        Out.Scens=Scens%>%
          mutate(Bo.mean=NA,Bo.CV=NA,Bf.mean=NA,Bf.CV=NA,r.mean=NA,r.cv=NA,Rdist=NA,
                 Kdist=NA,PsiDist=NA,bmsyk.mean=NA)
        Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
          Out.B.Bmsy=Out.F.Fmsy=Out.ensemble=vector('list',length(Store.sens))
        names(Out.ensemble)=names(Store.sens)
        
        for(s in 1:length(Store.sens))
        {
          print(paste("___________","JABBA Scenario",Scens$Scenario[s],"___________",names(dummy.store)))
          this.wd=paste(this.wd1,names(Store.sens)[s],sep='/')
          if(!dir.exists(this.wd))dir.create(this.wd)
          
          #Priors 
          bmsyk.mean=Scens$bmsyk[s]
          Proc.error=Scens$Proc.error[s]
          r.CV.multi=Scens$r.CV.multiplier[s]
          K.prior=c(Scens$K.mean[s],log(Scens$K.CV[s]))
          K.prior[1]=CMSY.K.estim*exp(rnorm(1,0,.1))
          if(names(dummy.store)=="smooth hammerhead")  #bump up K a bit as not enough biomass if using CMSY estimate
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
                     ASS=names(dummy.store),
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
          output=apply.JABBA(Ktch=input$Ktch,
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
                             burn.in= BURNIN)
          
          output$posteriors$P=output$posteriors$P[,1:nrow(ktch)]
          output$posteriors$BtoBmsy=output$posteriors$BtoBmsy[,1:nrow(ktch)]
          output$posteriors$SB=output$posteriors$SB[,1:nrow(ktch)]
          
          
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
            
            
            #Store trajectories  
            dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                              mods='JABBA',
                                              Type='Depletion',
                                              scen=Scens$Scenario[s],
                                              Katch=ktch$Total)
            Rel.biomas_JABBA=dummy$Dat
            
            #Display Priors vs Posteriors for base case scenario (S1)
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
            
            
            
          }
          
        }
      }
    } #end i
    fn.show.density(d=do.call(rbind,out),NCOL=1)
    ggsave(paste(out.here,"Prior.and.posterior_JABBA.tiff",sep=''),width = 12,height = 14, dpi = 300, compression = "lzw")
    
  }
  
  
  #5. SSS-MC assessment (Cope (2013)) 
  dothis=TRUE
  if(dothis)     #parallelised: 1.4 secs per iteration-species-scenario (otherwise 6 secs)
  {
    do.parallel.SS3=FALSE
    if(!do.parallel.SS3)
    {
      dummy.store=vector('list',1)
      names(dummy.store)=Keep.species[i]
      dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
        dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.Kobe.probs=
        dummy.store.sens.table=dummy.store.ensemble=dummy.store.warnings=
        dummy.store.convergence=dummy.store
      for(i in i)
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
          Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=Out.B.Bmsy=
            Out.F.Fmsy=store.warnings=store.convergence=Out.ensemble=vector('list',length(Store.sens))
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
            m.sd=List.sp[[i]]$Sens.test$DBSRA$Msd[s]
            
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
                if(Keep.species[i]=="grey nurse shark")  #too low h for greynurse
                {
                  h.sim[n]=rtruncnorm(1,a=Min.h.shark,b=Max.h.shark, 
                                      mean=Out.Scens$Steepness[s],
                                      sd=max(Out.Scens$Steepness.sd[s],0.01))
                }
                m.sim[n]=LH.sim$M[n] 
              }
              
              Depletion.year=ktch$finyear[nrow(ktch)]
              Btklow=List.sp[[i]]$FINALBIO[1]
              if(Neim=="milk shark")
              {
                Depletion.year=1991  #set to end of catch period to allow convergence
                Btklow=0.2
              }
              
              dep.f[n]=runif(1,Btklow,List.sp[[i]]$FINALBIO[2])
              Scenario=Scens[s,]%>%
                mutate(Mmean=m.sim[n],
                       Steepness=h.sim[n],
                       Initial.dpl=runif(1,List.sp[[i]]$STARTBIO[1],List.sp[[i]]$STARTBIO[2]),
                       Final.dpl=dep.f[n],
                       Ln_R0_init=runif(1,Ln_R0_max*.2,Ln_R0_max*.6))
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
                           depletion.yr=Depletion.year)
              
              #3. Run SS3
              fn.run.SS(where.inputs=this.wd1,
                        where.exe=handl_OneDrive('SS3/ss_win.exe'),
                        args='-nohess')  #no Hessian estimation as uncertainty generated thru MC procedure
              
              #4. Bring in outputs
              Report[[n]]=SS_output(this.wd1,covar=F,forecast=F,readwt=F,checkcor=F)
              Report[[n]]$Final.dpl=data.frame(Distribuion=c('Prior','Posterior'),
                                               Value=c(Scenario$Final.dpl,Report[[n]]$current_depletion),
                                               Parameter=rep('Final depletion',2))
              Report[[n]]$Initial.dpl=data.frame(Distribuion=c('Prior','Posterior'),
                                                 Value=c(Scenario$Initial.dpl,Report[[n]]$derived_quants[grep('Bratio',Report[[n]]$derived_quants$Label)[1],'Value']),
                                                 Parameter=rep('Initial depletion',2))
              Report[[n]]$Ln.R0=data.frame(Distribuion=c('Prior','Posterior'),
                                           Value=with(Report[[n]]$estimated_non_dev_parameters,c(Init,Value)),
                                           Parameter=rep('Ln.R0',2))
              
              
              #ss.std=read.table(paste(this.wd1,'ss.std',sep='/')) https://vlab.noaa.gov/web/stock-synthesis/public-forums/-/message_boards/view_message/11664630
              
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
              if(Converged$Keep=='Yes') n=n+1
            }
            Accept.rate_SSS=data.frame(
              Acceptance=100*n.sims_SS/trials,
              Scenario=Scens$Scenario[s])
            
            if(length(Report)>0)
            {
              
              
              #Store trajectories
              dummy=fn.ktch.only.get.timeseries(d=Report,
                                                mods="SSS",
                                                Type='Depletion',
                                                scen=Scens$Scenario[s],
                                                Katch=ktch$Tonnes)
              Rel.biomas_SSS=dummy$Dat
            }
            
            
          } # end s loop
        }
      } #end i
      
      
      
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
      ggsave(paste(out.here,"Prior.and.posterior_SS3.tiff",sep=''),width = 12,height = 14, dpi = 300, compression = "lzw")
      rm(this.wd1,Report)
      
    }
  }
  
  
  if(exists('Rel.biomas_DBSRA'))
  {
    pipi=rbind(Rel.biomas_DBSRA%>%mutate(Model=paste(Model," (acc. rate=",round(Accept.rate_DBSRA$Acceptance),'%)',sep='')),
               Rel.biomas_CMSY%>%mutate(Model=paste(Model," (acc. rate=",round(Accept.rate_CMSY$Acceptance),'%)',sep='')))
    
  }
  if(exists('Rel.biomas_JABBA')) pipi=rbind(pipi,Rel.biomas_JABBA)
  if(exists('Rel.biomas_DBSRA') & exists('Rel.biomas_SSS')) pipi=rbind(pipi,Rel.biomas_SSS)
  if(!exists('Rel.biomas_DBSRA') & exists('Rel.biomas_SSS')) pipi=Rel.biomas_SSS%>%mutate(Model=paste(Model," (acc. rate=",round(Accept.rate_SSS$Acceptance),'%)',sep=''))
  pipi%>%
  ggplot(aes(year,median,color=Model,fill=Model))+facet_wrap(~Model,ncol=1)+
    geom_line()+geom_ribbon(aes(ymin = lower.95, ymax = upper.95), alpha = 0.3)+ylim(0,1)+
    theme(legend.position = "none")
  ggsave(paste(out.here,"Biomass_SSS.tiff",sep=''),width = 8,height = 8, dpi = 300, compression = "lzw")
  
}
toc()