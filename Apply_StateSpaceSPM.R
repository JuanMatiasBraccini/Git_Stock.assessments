#---Run models -------------------------------------------------
#note: r prior tail should be less than 1.5 (Alex Hesp)
n.State.Space.SPM=length(state.space.SPM)  
State.Space.SPM=vector('list',n.State.Space.SPM)
names(State.Space.SPM)=state.space.SPM


#Auxiliary data - effort trends
#note: allocated fishery effort according to species distribution
Fish.power=data.frame(FINYEAR=Effort.monthly$FINYEAR)%>%
              arrange(FINYEAR)%>%
              mutate(Fish.pow.inc=as.numeric(substr(FINYEAR,1,4)),
                     Fish.pow.inc=1+c(seq(0,0.4,by=.02),rep(0.4,length(1996:Fish.pow.inc[length(Fish.pow.inc)]))))
auxiliary.effort.tdgdlf=Effort.monthly%>%
              left_join(Fish.power,by='FINYEAR')%>% #add the effort creep assumed for cpue standardisation
              mutate(Total=Total*Fish.pow.inc,
                     Rel.effort=Total/max(Total),
                     Year=as.numeric(substr(FINYEAR,1,4)))%>%
              filter(Year%in%unique(ktch.combined$finyear))%>%
              dplyr::select(Year,Rel.effort)
auxiliary.effort.nsf=Effort.monthly.north%>%
              rename(Total='Hook days')%>%
              mutate(Rel.effort=Total/max(Total),
                     Year=as.numeric(substr(FINYEAR,1,4)))%>%
              filter(Year%in%unique(ktch.combined$finyear))%>%
              dplyr::select(Year,Rel.effort)
auxiliary.effort.nsf.0=auxiliary.effort.nsf[1:length(2009:Last.yr.ktch.numeric),]%>%    #add 0 effort post 2008-09
              mutate(Year=2009:Last.yr.ktch.numeric,
                     Rel.effort=0)
auxiliary.effort.nsf=rbind(auxiliary.effort.nsf,auxiliary.effort.nsf.0)

Species.with.cpue=names(which(!sapply(Catch.rate.series,is.null)))
Species.non.representative.cpue=c("lemon shark","pigeye shark","scalloped hammerhead",
                                  "narrow sawfish","green sawfish")
Species.with.cpue=subset(Species.with.cpue,!Species.with.cpue%in%Species.non.representative.cpue)

Auxiliary.eff.dist=vector('list',length(Species.with.cpue))
for(o in 1:length(Species.with.cpue))
{
  dis.fish=unique(gsub("\\..*","",names(Catch.rate.series[[Species.with.cpue[o]]])))
  Fishries=dis.fish[grep(paste(c('TDGDLF','NSF'),collapse='|'),dis.fish)]
  if(length(Fishries)==0 & Species.with.cpue[o]%in%c("milk shark")) Fishries='NSF'
  if(length(Fishries)>1) Fishries=paste(Fishries,collapse=' & ')
  Auxiliary.eff.dist[[o]]=data.frame(Species=Species.with.cpue[o],Fishery=Fishries)
}
Auxiliary.eff.dist=do.call(rbind,Auxiliary.eff.dist)

tic("timer")
for(w in 1:length(State.Space.SPM))
{
  # JABBA (Winker et al 2018)   
  #summary of method: https://github.com/jabbamodel/JABBA
  if(names(State.Space.SPM)[w]=="JABBA")  #113 sec per species-scenario
  {
    if(do.parallel.JABBA)
    {
      set.seed(1234)
      progress <- function(n) cat(Keep.species[n],sprintf(": JABBA-CPUE fit complete-----------\n", n))
      opts <- list(progress = progress)
      cl <- makeCluster(detectCores()-1)
      registerDoSNOW(cl)
      out.species=foreach(i= 1:N.sp,.options.snow = opts,.errorhandling="stop",
                          .packages=c('ggmcmc','tseries','webr','ggpubr','Hmisc','JABBA','tidyverse')) %dopar%
      {
        Neim=Keep.species[i]
        if(!is.null(Catch.rate.series[[i]]))
        {
          this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                        capitalize(List.sp[[i]]$Name),"/",AssessYr,"/JABBA CPUE",sep='')
          if(!dir.exists(this.wd))dir.create(this.wd)
          
          #Catch
          ktch=ktch.combined%>%
            filter(Name==Neim)%>%
            rename(Year=finyear,
                   Total=Tonnes)%>%
            ungroup()%>%
            dplyr::select(Year,Total)%>%
            arrange(Year)%>%
            data.frame
          if(Neim=='milk shark') ktch$Total[1:10]=1 #convergence issues with first 10 years for milk shark, but catch is 0 anyways
          
          #CPUE series
          CPUE=compact(Catch.rate.series[[i]])
          DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
          if(length(DROP)>0)CPUE=CPUE[-DROP]
          if(Neim%in%survey_not.representative & "Survey"%in%names(CPUE)) CPUE=CPUE[-grep("Survey",names(CPUE))]
          if(Neim%in%NSF_not.representative & "NSF"%in%names(CPUE)) CPUE=CPUE[-grep("NSF",names(CPUE))]
          if(Neim%in%tdgdlf_not.representative & "TDGDLF"%in%names(CPUE)) CPUE=CPUE[-grep("TDGDLF",names(CPUE))]
           
          #reset very low CVs       
          #note:very low CVs in stand cpues, hence use Francis CVs 
          if(increase.CV.JABBA)
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
              CPUE[[j]]=CPUE[[j]]%>%mutate(CV.Original=CV)
              if(mean(CPUE[[j]]$CV,na.rm=T)>=default.CV)
              {
                NewCVs$CV.var.adj[j]=0
              }
              if(mean(CPUE[[j]]$CV,na.rm=T)<default.CV)
              {
                newcv=NewCVs$CV.Adjusted[,match(c("yr.f",names(CPUE)[j]),names(NewCVs$CV.Adjusted))]%>%
                  filter(yr.f%in%CPUE[[j]]$yr.f)
                CPUE[[j]]=CPUE[[j]]%>%mutate(CV=newcv[,2])
              }
            }
            
            tiff(file=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[i]]$Name),
                            "/",AssessYr,"/Adjusted cpue CVs.tiff",sep=''),
                 width = 1600, height = 2400,units = "px", res = 300, compression = "lzw")
            par(mfrow=n2mfrow(length(CPUE)))
            for(x in 1:length(CPUE))
            {
              plot(CPUE[[x]]$yr.f,CPUE[[x]]$CV,pch=19,
                   main=paste(names(CPUE)[x]," (Var. Adj. factor=",round(NewCVs$CV.var.adj[x],3),")",sep=''),
                   ylim=c(0,max(CPUE[[x]]$CV,na.rm=TRUE)),ylab='CV',xlab='Year')
              points(CPUE[[x]]$yr.f,CPUE[[x]]$CV.Original)
              abline(h=default.CV,col=2)
              if(x==1) legend("bottomleft",c('Adjusted CV','Original CV',"Default mean CV"),bty='n',
                              pch=c(19,1,NA),text.col=c(1,1,2))
              CPUE[[x]]=CPUE[[x]]%>%dplyr::select(-CV.Original)
            }
            dev.off() 
          }
          

          len.cpue=length(CPUE)
          MAX.CV=List.sp[[i]]$MAX.CV
          if(len.cpue>0)
          {
            Cpues=data.frame(Year=ktch$Year)
            Cpues.SE=Cpues
            for(x in 1:length(CPUE))    
            {
              dd=CPUE[[x]][,grep(paste(c('yr.f','Mean','MeAn','CV'),collapse="|"),names(CPUE[[x]]))]%>%
                relocate(yr.f)
              if(drop.large.CVs)
              {
                iid=which(dd$CV>MAX.CV)
                dd$Mean[iid]=NA
                dd$CV[iid]=NA 
              }
              add.yrs=sort(unique(ktch$Year))
              AdD=data.frame(x=add.yrs[which(!add.yrs%in%dd$yr.f)],y=NA,z=NA)
              colnames(AdD)=colnames(dd)
              dd=rbind(dd,AdD)%>%
                arrange(yr.f)%>%
                rename(Year=yr.f)
              colnames(dd)[2:3]=c(names(CPUE)[x],paste(names(CPUE)[x],'CV',sep='.'))
              Cpues=left_join(Cpues,dd%>%dplyr::select(Year,names(CPUE)[x]),by='Year')
              Cpues.SE=left_join(Cpues.SE,dd%>%dplyr::select(Year,paste(names(CPUE)[x],'CV',sep='.')),by='Year')
            }
            #block qs
            if(Neim=="whiskery shark")  
            {
              if(Whiskery.q.periods==2) #2 catchability periods
              {
                Cpues=Cpues%>%
                  mutate(TDGDLF.monthly=ifelse(Year%in%List.sp[[i]]$Yr_q_change_transition,NA,TDGDLF.monthly),   
                         TDGDLF.monthly2=ifelse(Year>List.sp[[i]]$Yr_q_change,TDGDLF.monthly,NA),   
                         TDGDLF.monthly=ifelse(Year<=List.sp[[i]]$Yr_q_change,TDGDLF.monthly,NA))
                
                Cpues.SE=Cpues.SE%>%
                  mutate(TDGDLF.monthly.CV=ifelse(Year%in%List.sp[[i]]$Yr_q_change_transition,NA,TDGDLF.monthly.CV),
                         TDGDLF.monthly2.CV=ifelse(Year>List.sp[[i]]$Yr_q_change,TDGDLF.monthly.CV,NA),
                         TDGDLF.monthly.CV=ifelse(Year<=List.sp[[i]]$Yr_q_change,TDGDLF.monthly.CV,NA))
              }
              if(Whiskery.q.periods==3) #3 catchability periods
              {
                Cpues=Cpues%>%
                  mutate(TDGDLF.monthly2=ifelse(Year%in%List.sp[[i]]$Yr_q_change_transition,TDGDLF.monthly,NA), 
                         TDGDLF.monthly3=ifelse(Year>max(List.sp[[i]]$Yr_q_change_transition),TDGDLF.monthly,NA),
                         TDGDLF.monthly=ifelse(Year<min(List.sp[[i]]$Yr_q_change_transition),TDGDLF.monthly,NA))
                
                Cpues.SE=Cpues.SE%>%
                  mutate(TDGDLF.monthly2.CV=ifelse(Year%in%List.sp[[i]]$Yr_q_change_transition,TDGDLF.monthly.CV,NA), 
                         TDGDLF.monthly3.CV=ifelse(Year>max(List.sp[[i]]$Yr_q_change_transition),TDGDLF.monthly.CV,NA),
                         TDGDLF.monthly.CV=ifelse(Year<min(List.sp[[i]]$Yr_q_change_transition),TDGDLF.monthly.CV,NA))
              }
            }
            if(Neim=='gummy shark' & Gummy.q.periods>1)     
            {
              Cpues=Cpues%>%
                mutate(TDGDLF.monthly2=ifelse(Year%in%List.sp[[i]]$Yr_second_q,TDGDLF.monthly,NA),   #second monthly q due to increase in cpue and catch unaccounted by cpue stand.
                       TDGDLF.monthly=ifelse(Year%in%List.sp[[i]]$Yr_second_q,NA,TDGDLF.monthly))
              Cpues.SE=Cpues.SE%>%
                mutate(TDGDLF.monthly2.CV=ifelse(Year%in%List.sp[[i]]$Yr_second_q,TDGDLF.monthly.CV,NA),
                       TDGDLF.monthly.CV=ifelse(Year%in%List.sp[[i]]$Yr_second_q,NA,TDGDLF.monthly.CV))
            }
            
            #Scenarios
            Scens=List.sp[[i]]$Sens.test$JABBA%>%
              mutate(Species=capitalize(Neim))
            Store.sens=vector('list',nrow(Scens))
            names(Store.sens)=Scens$Scenario
            this.wd1=this.wd
            Out.Scens=Scens%>%
              mutate(Bo.mean=NA,Bo.CV=NA,r.mean=NA,r.cv=NA,Rdist=NA,Kdist=NA,PsiDist=NA,bmsyk.mean=NA)
            Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
              Out.B.Bmsy=Out.F.Fmsy=vector('list',length(Store.sens))
            
            for(s in 1:length(Store.sens))
            {
              print(paste("___________","JABBA-CPUE Scenario",Scens$Scenario[s],"___________",Neim))
              this.wd=paste(this.wd1,names(Store.sens)[s],sep='/')
              if(!dir.exists(this.wd))dir.create(this.wd)
              
              #cpues
              #remove daily years 
              Cpues1=Cpues
              Cpues1.SE=Cpues.SE
              if(!is.na(Scens$Daily.cpues[s]) & 'TDGDLF.daily'%in%names(Cpues))
              {
                rid.of=as.numeric(unlist(str_split(Scens$Daily.cpues[s], "&")))
                Cpues1=Cpues%>%
                  mutate(TDGDLF.daily=ifelse(Year%in%rid.of,NA,TDGDLF.daily))
                Cpues1.SE=Cpues.SE%>%
                  mutate(TDGDLF.daily.CV=ifelse(Year%in%rid.of,NA,TDGDLF.daily.CV))
              }
              
              #Priors 
              bmsyk.mean=Scens$bmsyk[s]
              Proc.error=Scens$Proc.error[s]
              r.CV.multi=Scens$r.CV.multiplier[s]
              K.prior=c(Scens$K.mean[s],log(Scens$K.CV[s]))  
              Mn=min(Scens$r[s],Max.r.value)  #some life history pars yield unrealistically high r
              Mn=max(Mn,Min.r.value)          #allow minimum recruitment
              r.prior=c(Mn, (Scens$r.sd[s]/Mn)*Scens$r.CV.multiplier[s])
              Ktch.CV=Scens$Ktch.CV[s]
              Bint=runif(1000,List.sp[[i]]$STARTBIO[1],List.sp[[i]]$STARTBIO[2])
              Bint.mean=mean(Bint)
              Bint.CV=sd(Bint)/Bint.mean
              psi.prior=c(Bint.mean,Bint.CV)
              b.prior=FALSE  #N.A., only use b.prior for catch-only approach options: FALSE, 0.3, NA, "bk"
              
              #allocate relevant fishing effort 
              if(use.auxiliary.effort)
              {
                applicable.fishery=Auxiliary.eff.dist%>%filter(Species==Neim)%>%pull(Fishery)
                if(applicable.fishery=='TDGDLF')
                {
                  auxiliary.effort=auxiliary.effort.tdgdlf
                }
                if(applicable.fishery=='NSF')
                {
                  auxiliary.effort=auxiliary.effort.nsf
                }
                if(applicable.fishery=="TDGDLF & NSF")
                {
                  auxiliary.effort=full_join(auxiliary.effort.tdgdlf,auxiliary.effort.nsf,by='Year')%>%
                    rowwise()%>% 
                    mutate(Rel.effort = sum(c_across(Rel.effort.x:Rel.effort.y), na.rm = T))%>%
                    ungroup()%>%
                    data.frame%>%
                    mutate(Rel.effort=Rel.effort/max(Rel.effort))%>%
                    dplyr::select(Year,Rel.effort)
                }
                add.yrs=sort(unique(ktch$Year))
                AdD=data.frame(x=add.yrs[which(!add.yrs%in%auxiliary.effort$Year)],Rel.effort=NA)
                colnames(AdD)=colnames(auxiliary.effort)
                auxiliary.effort=rbind(auxiliary.effort,AdD)%>%
                  arrange(Year)
                auxiliary.effort.se=auxiliary.effort%>% 
                  mutate(Rel.effort=ifelse(Rel.effort==0,NA,Rel.effort))%>%
                  rename(se=Rel.effort)%>%
                  mutate(se=ifelse(!is.na(se),Ktch.CV,se))   #assume same effort error as for catch error
                auxiliary.type="effort"
              }else
              {
                auxiliary.effort=auxiliary.effort.se=auxiliary.type=NULL
              }
              
              #Put inputs together 
              input=list(Ktch=ktch,
                         MDL="Pella_m",  #Pella_m if estimating m; other options: "Schaefer", "Fox", "Pella"
                         Ktch.CV=Ktch.CV,
                         ASS=Neim,
                         Rdist = Rdist,
                         Rprior = r.prior,
                         r.CV.multiplier=r.CV.multi,
                         Kdist= Kdist,
                         Kprior=K.prior,    
                         PsiDist= PsiDist,
                         Psiprior=psi.prior,   
                         Bprior=b.prior,    
                         BMSYK=bmsyk.mean,
                         shape.CV=r.prior[2])  
              
              #Run model    
              THIN=5
              CHAINS=2
              BURNIN=min(0.15*Scens$Sims[s],5000)
              JABBA.run=apply.JABBA(Ktch=input$Ktch,
                                     CPUE=Cpues1,
                                     CPUE.SE=Cpues1.SE,
                                     auxil=auxiliary.effort,
                                     auxil.se=auxiliary.effort.se,
                                     auxil.type=auxiliary.type,
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
                                     Shape.CV=input$shape.CV,
                                     output.dir=this.wd,
                                     outfile="fit_Appendix",
                                     Sims=Scens$Sims[s],
                                     Proc.error.JABBA=Proc.error,
                                     Obs.Error=Obs.Err.JABBA,
                                     thinning = THIN,
                                     nchains = CHAINS,
                                     burn.in= BURNIN)
              output=JABBA.run$fit
                
              #Display model fits
              fn.fig(paste(this.wd,"Fit_stdresiduals",sep='/'), 2400, 2400)
              jbplot_stdresiduals(output)
              dev.off()
              fn.fig(paste(this.wd,"Fit_dresiduals",sep='/'), 2400, 2400)
              jbplot_residuals(output)
              dev.off()
              fn.fig(paste(this.wd,"Fit_cpues",sep='/'), 2400, 2400)
              jbplot_cpuefits(output)
              dev.off()
              fn.fig(paste(this.wd,"Fit_cpues_log",sep='/'), 2400, 2400)
              jbplot_logfits(output) 
              dev.off()
              fn.fig(paste(this.wd,"Fit_catcherror",sep='/'), 2400, 2400)
              jbplot_catcherror(output)
              dev.off()
              fn.fig(paste(this.wd,"Fit_runstest",sep='/'), 2400, 2400)
              jbplot_runstest(output)
              dev.off()
              fn.fig(paste(this.wd,"Summary",sep='/'), 2400, 2400)
              jbplot_summary(output)
              dev.off()
              
              #Get posteriors
              output$posteriors$P=output$posteriors$P[,1:nrow(ktch)]
              output$posteriors$BtoBmsy=output$posteriors$BtoBmsy[,1:nrow(ktch)]
              output$posteriors$SB=output$posteriors$SB[,1:nrow(ktch)]
              Store.sens[[s]]=list(input=input,output=output)
              
              #find pattern in residuals
              RES=as.data.frame(t(output$residuals))%>%
                mutate(across(everything(),~ifelse(.x<0,-1,1)),
                       across(everything(),~as.factor(.x)))
              store.runs.test=vector('list',ncol(RES))
              names(store.runs.test)=names(RES)
              for(x in 1:ncol(RES)) store.runs.test[[x]]=runs.test(RES[!is.na(RES[,x]),x])
              
              #Store Scenarios
              Out.Scens$Bo.mean[s]=Bint.mean
              Out.Scens$Bo.CV[s]=Bint.CV
              Out.Scens$r.mean[s]=Mn
              Out.Scens$r.cv[s]=Scens$r.sd[s]/Mn
              Out.Scens$Rdist[s]=Rdist
              Out.Scens$Kdist[s]=Kdist
              Out.Scens$PsiDist[s]=PsiDist
              Out.Scens$bmsyk.mean[s]=bmsyk.mean
              
              #Store estimates  
              d1=Store.sens[[s]]$output$estimates
              d1=d1[,grep(paste(c('mu','lci','uci'),collapse='|'),colnames(d1))]%>%
                data.frame
              colnames(d1)=c("Median","Lower.95","Upper.95")
              d1=d1%>%
                mutate(Parameter=rownames(d1),
                       Model='JABBA CPUE',
                       Scenario=names(Store.sens)[s])%>%
                relocate(Model,Scenario,Parameter,Lower.95)%>%
                filter(Parameter%in%c("K","r","psi","Hmsy","SBmsy","MSY"))   
              Out.estimates[[s]]=d1
              
              #Store trajectories  
              dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                mods=names(State.Space.SPM)[w],  
                                                Type='Depletion',
                                                scen=Scens$Scenario[s],
                                                Katch=ktch$Total)
              Out.rel.biom[[s]]=dummy$Dat
              Out.probs.rel.biom[[s]]=dummy$Probs
              
              dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                mods=names(State.Space.SPM)[w],
                                                Type='F.series',
                                                scen=Scens$Scenario[s],
                                                Katch=ktch$Total)
              Out.f.series[[s]]=dummy$Dat
              
              dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                mods=names(State.Space.SPM)[w],
                                                Type='B.Bmsy',
                                                scen=Scens$Scenario[s],
                                                Katch=ktch$Total)
              Out.B.Bmsy[[s]]=dummy$Dat
              
              dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                mods=names(State.Space.SPM)[w],
                                                Type='F.Fmsy',
                                                scen=Scens$Scenario[s],
                                                Katch=ktch$Total)
              Out.F.Fmsy[[s]]=dummy$Dat
              rm(dummy)
              
              if(Scens$Scenario[s]=='S1') Out.Kobe.probs=Store.sens[[s]]$output$kobe%>%dplyr::select(stock,harvest)
              
              #Evaluate posteriors and fit diagnostics
              if(Scens$Scenario[s]=='S1')
              {
                outfilename=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                                  capitalize(Neim),"/",AssessYr,"/",'JABBA CPUE/',Scens$Scenario[s],sep='')
                
                #MCMC trace-plots, Geweke (2003) and Gelman and Rubin (1992) diagnostics
                #https://cran.r-project.org/web/packages/ggmcmc/vignettes/using_ggmcmc.html
                if(do.MCMC.diagnostics)
                {
                  these.post.pars=c('K','r','psi')  #estimable pars
                  NN=nrow(output$pars_posterior)/CHAINS
                  Post=Store.sens[[s]]$output$pars_posterior%>%
                    dplyr::select(all_of(these.post.pars))%>%
                    mutate(Iteration=rep(1:NN,CHAINS),
                           Chain=rep(1:CHAINS,each=NN))%>%
                    gather(Parameter,value,-Iteration,- Chain)
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
                  
                  ggsave(paste(outfilename,"/Fit_Bayesian convergence diagnostics.tiff",sep=''),
                         width = 16,height = 8, dpi = 300, compression = "lzw")
                  
                  
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
                  if(dummy.prior$psi$dist=='lnorm') dummy.prior$psi$mean=log(dummy.prior$psi$mean)
                  dummy.post=Store.sens[[s]]$output$pars_posterior%>%
                    dplyr::select(r,K,psi)
                  names(dummy.prior)=names(dummy.post)
                  pars=colnames(dummy.post)
                  out=vector('list',length(pars))
                  for(p in 1:length(pars))
                  {
                    Value.prior=fn.prior(d=dummy.prior[[pars[p]]])
                    if(pars[p]=="K") Value.prior=fn.prior(d=dummy.prior[[pars[p]]],MAX=max(stuff$Ktch$Total,na.rm=T)*50)
                    out[[p]]=rbind(data.frame(Distribuion="Prior",
                                              Value=Value.prior),
                                   data.frame(Distribuion="Posterior",
                                              Value=dummy.post[[pars[p]]]))%>%
                      mutate(Parameter=pars[p])
                  }
                  rm(stuff,dummy.prior,dummy.post)
                  
                  fn.show.density(d=do.call(rbind,out),NCOL=1)
                  ggsave(paste(outfilename,"/Prior.and.posterior.tiff",sep=''),
                         width = 12,height = 14, dpi = 300, compression = "lzw")
                  
                }
                
                #Hindcasting and retrospective analysis  
                if(do.hindcasting)
                {
                  #retrospective
                  hc = hindcast_jabba(jbinput=JABBA.run$jbinput,fit=JABBA.run$fit,peels=1:5)
                  
                  tiff(paste(outfilename,"/Fit_retrospective.tiff",sep=''),
                       width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
                  mohns=jbplot_retro(hc) 
                  dev.off()
                  write.csv(mohns,paste(outfilename,"/Fit_retrospective_mohns.csv",sep=''))
                  
                  #hindcasting  
                  tiff(paste(outfilename,"/Fit_hindcasting.tiff",sep=''),
                       width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
                  jbplot_hcxval(hc)
                  dev.off()
                  
                }
                
                #Plot observed (black) and model-predicted (red) abundance series
                nn.s=dim(output$cpue.hat)[3]
                x.dummy=vector('list',nn.s)
                names(x.dummy)=colnames(Cpues)[-1]
                for(pp in 1:nn.s)
                {
                  x.dummy[[pp]]=as.data.frame(output$cpue.hat[,,pp])%>%
                    mutate(Data=names(x.dummy)[pp])
                }
                x.dummy=do.call(rbind,x.dummy)%>%
                  tibble::rownames_to_column(var='Year')%>%
                  mutate(Year=as.numeric(str_extract(Year, "[[:digit:]]+")))%>%
                  rename(mu.se=se)%>%
                  dplyr::select(Year,Data,mu,mu.se)
                
                Cpues%>%
                  gather(Data,cpue,-Year)%>%
                  left_join(Cpues.SE%>%
                              gather(Data,cpue.cv,-Year)%>%
                              mutate(Data=str_replace(Data,".CV", "")),
                            by=c('Year','Data'))%>%
                  filter(!is.na(cpue))%>%
                  left_join(x.dummy,by=c('Year','Data'))%>%
                  ggplot(aes(Year,cpue))+
                  geom_point()+
                  geom_errorbar(aes(ymin = cpue-cpue.cv, ymax = cpue+cpue.cv), width = 0.2)+
                  geom_line(aes(Year+.25,mu),color='red', alpha=.5)+
                  geom_errorbar(aes(x=Year+.25,ymin = mu-mu.se, ymax = mu+mu.se), alpha=.5,width = 0.2,color='red')+
                  facet_wrap(~Data,scales='free')
                ggsave(paste(outfilename,"/Observed_predicted_cpue.tiff",sep=''),
                       width = 6,height = 6, dpi = 300, compression = "lzw")
                
              }
              
            }# end s loop
            
            Out.Scens=Out.Scens%>%
              dplyr::select(Species,Scenario,
                            Rdist,r.mean,r.cv,
                            Kdist,K.mean,K.CV,
                            PsiDist,Bo.mean,Bo.CV,
                            bmsyk.mean,Proc.error,Ktch.CV,
                            Daily.cpues)%>%
              mutate(K.mean=round(K.mean),
                     r.mean=round(r.mean,3),
                     bmsyk.mean=round(bmsyk.mean,3),
                     r.cv=round(r.cv,2),
                     Bo.mean=round(Bo.mean,2),
                     Bo.CV=round(Bo.CV,2),
                     Rdist=ifelse(Rdist=='lnorm','Lognormal',Rdist),
                     Kdist=ifelse(Kdist=='lnorm','Lognormal',Kdist),
                     PsiDist=ifelse(PsiDist=='lnorm','Lognormal',PsiDist))

            return(list(dummy.store.sens.table=Out.Scens,
                        dummy.store.rel.biom=do.call(rbind,Out.rel.biom),
                        dummy.store.probs.rel.biom=Out.probs.rel.biom,
                        dummy.store.f.series=do.call(rbind,Out.f.series),
                        dummy.store.B.Bmsy=do.call(rbind,Out.B.Bmsy),
                        dummy.store.F.Fmsy=do.call(rbind,Out.F.Fmsy),
                        dummy.store.Kobe.probs=Out.Kobe.probs, 
                        dummy.store.estimates=do.call(rbind,Out.estimates)))
                        
            rm(Out.Scens,Out.rel.biom,Out.probs.rel.biom,Out.f.series,
               Out.B.Bmsy,Out.F.Fmsy,Out.Kobe.probs,Out.estimates)
          }
        } 
      }
      stopCluster(cl)
      names(out.species)=Keep.species

      State.Space.SPM[[w]]$sens.table=fn.get.and.name(LISTA=out.species,x="dummy.store.sens.table")
      State.Space.SPM[[w]]$estimates=fn.get.and.name(LISTA=out.species,x="dummy.store.estimates")
      State.Space.SPM[[w]]$rel.biom=fn.get.and.name(LISTA=out.species,x="dummy.store.rel.biom")
      State.Space.SPM[[w]]$probs.rel.biom=fn.get.and.name(LISTA=out.species,x="dummy.store.probs.rel.biom")
      State.Space.SPM[[w]]$f.series=fn.get.and.name(LISTA=out.species,x="dummy.store.f.series")
      State.Space.SPM[[w]]$B.Bmsy=fn.get.and.name(LISTA=out.species,x="dummy.store.B.Bmsy")
      State.Space.SPM[[w]]$F.Fmsy=fn.get.and.name(LISTA=out.species,x="dummy.store.F.Fmsy")
      State.Space.SPM[[w]]$Kobe.probs=fn.get.and.name(LISTA=out.species,x="dummy.store.Kobe.probs")
      
      rm(out.species)
      
    }
    
    if(!do.parallel.JABBA)
    {
      dummy.store=vector('list',N.sp)     
      names(dummy.store)=Keep.species
      dummy.store.sens.table=dummy.store.estimates=dummy.store.rel.biom=dummy.store.probs.rel.biom=
        dummy.store.f.series=dummy.store.B.Bmsy=dummy.store.F.Fmsy=dummy.store.Kobe.probs=dummy.store
      
      for(i in 1:length(dummy.store))
      {
        Neim=names(dummy.store)[i]
        
        if(!is.null(Catch.rate.series[[i]]))
        {
          this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                        capitalize(List.sp[[i]]$Name),"/",AssessYr,"/JABBA CPUE",sep='')
          if(!dir.exists(this.wd))dir.create(this.wd)
          
          #Catch
          ktch=ktch.combined%>%
            filter(Name==Neim)%>%
            rename(Year=finyear,
                   Total=Tonnes)%>%
            ungroup()%>%
            dplyr::select(Year,Total)%>%
            arrange(Year)%>%
            data.frame
          if(Neim=='milk shark') ktch$Total[1:10]=1 #convergence issues with first 10 years for milk shark, but catch is 0 anyways
          
          #CPUE series
          CPUE=compact(Catch.rate.series[[i]])
          DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
          if(length(DROP)>0)CPUE=CPUE[-DROP]
          if(Neim%in%survey_not.representative & "Survey"%in%names(CPUE)) CPUE=CPUE[-grep("Survey",names(CPUE))]
          if(Neim%in%NSF_not.representative & "NSF"%in%names(CPUE)) CPUE=CPUE[-grep("NSF",names(CPUE))]
          if(Neim%in%tdgdlf_not.representative & "TDGDLF"%in%names(CPUE)) CPUE=CPUE[-grep("TDGDLF",names(CPUE))]
          
          #reset very low CVs       
          #note:very low CVs in stand cpues, hence increase obs error a bit 
          if(increase.CV.JABBA)
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
              CPUE[[j]]=CPUE[[j]]%>%mutate(CV.Original=CV)
              if(mean(CPUE[[j]]$CV,na.rm=T)>=default.CV)
              {
                NewCVs$CV.var.adj[j]=0
              }
              if(mean(CPUE[[j]]$CV,na.rm=T)<default.CV)
              {
                newcv=NewCVs$CV.Adjusted[,match(c("yr.f",names(CPUE)[j]),names(NewCVs$CV.Adjusted))]%>%
                  filter(yr.f%in%CPUE[[j]]$yr.f)
                CPUE[[j]]=CPUE[[j]]%>%mutate(CV=newcv[,2])
              }
            }
            
            tiff(file=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[i]]$Name),
                            "/",AssessYr,"/Adjusted cpue CVs.tiff",sep=''),
                 width = 1600, height = 2400,units = "px", res = 300, compression = "lzw")
            par(mfrow=n2mfrow(length(CPUE)))
            for(x in 1:length(CPUE))
            {
              plot(CPUE[[x]]$yr.f,CPUE[[x]]$CV,pch=19,
                   main=paste(names(CPUE)[x]," (Var. Adj. factor=",round(NewCVs$CV.var.adj[x],3),")",sep=''),
                   ylim=c(0,max(CPUE[[x]]$CV,na.rm=TRUE)),ylab='CV',xlab='Year')
              points(CPUE[[x]]$yr.f,CPUE[[x]]$CV.Original)
              abline(h=default.CV,col=2)
              if(x==1) legend("bottomleft",c('Adjusted CV','Original CV',"Default mean CV"),bty='n',
                              pch=c(19,1,NA),text.col=c(1,1,2))
              CPUE[[x]]=CPUE[[x]]%>%dplyr::select(-CV.Original)
            }
            dev.off() 
          }
          
          len.cpue=length(CPUE)
          MAX.CV=List.sp[[i]]$MAX.CV
          if(len.cpue>0)
          {
            Cpues=data.frame(Year=ktch$Year)
            Cpues.SE=Cpues
            for(x in 1:length(CPUE))    
            {
              dd=CPUE[[x]][,grep(paste(c('yr.f','Mean','MeAn','CV'),collapse="|"),names(CPUE[[x]]))]%>%
                relocate(yr.f)
              if(drop.large.CVs)
              {
                iid=which(dd$CV>MAX.CV)
                dd$Mean[iid]=NA
                dd$CV[iid]=NA 
              }
              add.yrs=sort(unique(ktch$Year))
              AdD=data.frame(x=add.yrs[which(!add.yrs%in%dd$yr.f)],y=NA,z=NA)
              colnames(AdD)=colnames(dd)
              dd=rbind(dd,AdD)%>%
                arrange(yr.f)%>%
                rename(Year=yr.f)
              colnames(dd)[2:3]=c(names(CPUE)[x],paste(names(CPUE)[x],'CV',sep='.'))
              Cpues=left_join(Cpues,dd%>%dplyr::select(Year,names(CPUE)[x]),by='Year')
              Cpues.SE=left_join(Cpues.SE,dd%>%dplyr::select(Year,paste(names(CPUE)[x],'CV',sep='.')),by='Year')
            }
            #block qs
            if(Neim=="whiskery shark")  
            {
              if(Whiskery.q.periods==2) #2 catchability periods
              {
                Cpues=Cpues%>%
                  mutate(TDGDLF.monthly=ifelse(Year%in%List.sp[[i]]$Yr_q_change_transition,NA,TDGDLF.monthly),   
                         TDGDLF.monthly2=ifelse(Year>List.sp[[i]]$Yr_q_change,TDGDLF.monthly,NA),   
                         TDGDLF.monthly=ifelse(Year<=List.sp[[i]]$Yr_q_change,TDGDLF.monthly,NA))
                
                Cpues.SE=Cpues.SE%>%
                  mutate(TDGDLF.monthly.CV=ifelse(Year%in%List.sp[[i]]$Yr_q_change_transition,NA,TDGDLF.monthly.CV),
                         TDGDLF.monthly2.CV=ifelse(Year>List.sp[[i]]$Yr_q_change,TDGDLF.monthly.CV,NA),
                         TDGDLF.monthly.CV=ifelse(Year<=List.sp[[i]]$Yr_q_change,TDGDLF.monthly.CV,NA))
              }
              if(Whiskery.q.periods==3) #3 catchability periods
              {
                Cpues=Cpues%>%
                  mutate(TDGDLF.monthly2=ifelse(Year%in%List.sp[[i]]$Yr_q_change_transition,TDGDLF.monthly,NA), 
                         TDGDLF.monthly3=ifelse(Year>max(List.sp[[i]]$Yr_q_change_transition),TDGDLF.monthly,NA),
                         TDGDLF.monthly=ifelse(Year<min(List.sp[[i]]$Yr_q_change_transition),TDGDLF.monthly,NA))
                
                Cpues.SE=Cpues.SE%>%
                  mutate(TDGDLF.monthly2.CV=ifelse(Year%in%List.sp[[i]]$Yr_q_change_transition,TDGDLF.monthly.CV,NA), 
                         TDGDLF.monthly3.CV=ifelse(Year>max(List.sp[[i]]$Yr_q_change_transition),TDGDLF.monthly.CV,NA),
                         TDGDLF.monthly.CV=ifelse(Year<min(List.sp[[i]]$Yr_q_change_transition),TDGDLF.monthly.CV,NA))
              }
            }
            if(Neim=='gummy shark' & Gummy.q.periods>1)     
            {
              Cpues=Cpues%>%
                mutate(TDGDLF.monthly2=ifelse(Year%in%List.sp[[i]]$Yr_second_q,TDGDLF.monthly,NA),   #second monthly q due to increase in cpue and catch unaccounted by cpue stand.
                       TDGDLF.monthly=ifelse(Year%in%List.sp[[i]]$Yr_second_q,NA,TDGDLF.monthly))
              Cpues.SE=Cpues.SE%>%
                mutate(TDGDLF.monthly2.CV=ifelse(Year%in%List.sp[[i]]$Yr_second_q,TDGDLF.monthly.CV,NA),
                       TDGDLF.monthly.CV=ifelse(Year%in%List.sp[[i]]$Yr_second_q,NA,TDGDLF.monthly.CV))
            }
            
            #Scenarios
            Scens=List.sp[[i]]$Sens.test$JABBA%>%
              mutate(Species=capitalize(Neim))
            Store.sens=vector('list',nrow(Scens))
            names(Store.sens)=Scens$Scenario
            this.wd1=this.wd
            Out.Scens=Scens%>%
              mutate(Bo.mean=NA,Bo.CV=NA,r.mean=NA,r.cv=NA,Rdist=NA,Kdist=NA,PsiDist=NA,bmsyk.mean=NA)
            Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
              Out.B.Bmsy=Out.F.Fmsy=vector('list',length(Store.sens))
            
            for(s in 1:length(Store.sens))
            {
              print(paste("___________","JABBA-CPUE Scenario",Scens$Scenario[s],"___________",Neim))
              this.wd=paste(this.wd1,names(Store.sens)[s],sep='/')
              if(!dir.exists(this.wd))dir.create(this.wd)
              
              #cpues
              #remove daily years
              Cpues1=Cpues
              Cpues1.SE=Cpues.SE
              if(!is.na(Scens$Daily.cpues[s]) & 'TDGDLF.daily'%in%names(Cpues))
              {
                rid.of=as.numeric(unlist(str_split(Scens$Daily.cpues[s], "&")))
                Cpues1=Cpues%>%
                  mutate(TDGDLF.daily=ifelse(Year%in%rid.of,NA,TDGDLF.daily))
                Cpues1.SE=Cpues.SE%>%
                  mutate(TDGDLF.daily.CV=ifelse(Year%in%rid.of,NA,TDGDLF.daily.CV))
              }
              
              #Priors 
              bmsyk.mean=Scens$bmsyk[s]
              Proc.error=Scens$Proc.error[s]
              r.CV.multi=Scens$r.CV.multiplier[s]
              K.prior=c(Scens$K.mean[s],log(Scens$K.CV[s]))  
              Mn=min(Scens$r[s],Max.r.value)  #some life history pars yield unrealistically high r
              Mn=max(Mn,Min.r.value)          #allow minimum recruitment
              r.prior=c(Mn, (Scens$r.sd[s]/Mn)*Scens$r.CV.multiplier[s])
              Ktch.CV=Scens$Ktch.CV[s]
              Bint=runif(1000,List.sp[[i]]$STARTBIO[1],List.sp[[i]]$STARTBIO[2])
              Bint.mean=mean(Bint)
              Bint.CV=sd(Bint)/Bint.mean
              psi.prior=c(Bint.mean,Bint.CV)
              b.prior=FALSE  #N.A. only use b.prior for catch-only approach optons: c(FALSE, 0.3, NA, "bk")
              
              #allocate relevant fishing effort 
              if(use.auxiliary.effort)
              {
                applicable.fishery=Auxiliary.eff.dist%>%filter(Species==Neim)%>%pull(Fishery)
                if(applicable.fishery=='TDGDLF')
                {
                  auxiliary.effort=auxiliary.effort.tdgdlf
                }
                if(applicable.fishery=='NSF')
                {
                  auxiliary.effort=auxiliary.effort.nsf
                }
                if(applicable.fishery=="TDGDLF & NSF")
                {
                  auxiliary.effort=full_join(auxiliary.effort.tdgdlf,auxiliary.effort.nsf,by='Year')%>%
                    rowwise()%>% 
                    mutate(Rel.effort = sum(c_across(Rel.effort.x:Rel.effort.y), na.rm = T))%>%
                    ungroup()%>%
                    data.frame%>%
                    mutate(Rel.effort=Rel.effort/max(Rel.effort))%>%
                    dplyr::select(Year,Rel.effort)
                }
                add.yrs=sort(unique(ktch$Year))
                AdD=data.frame(x=add.yrs[which(!add.yrs%in%auxiliary.effort$Year)],Rel.effort=NA)
                colnames(AdD)=colnames(auxiliary.effort)
                auxiliary.effort=rbind(auxiliary.effort,AdD)%>%
                  arrange(Year)
                auxiliary.effort.se=auxiliary.effort%>% 
                  mutate(Rel.effort=ifelse(Rel.effort==0,NA,Rel.effort))%>%
                  rename(se=Rel.effort)%>%
                  mutate(se=ifelse(!is.na(se),Ktch.CV,se))   #assume same effort error as for catch error
                auxiliary.type="effort"
              }else
              {
                auxiliary.effort=auxiliary.effort.se=auxiliary.type=NULL
              }
              
              #Put inputs together 
              input=list(Ktch=ktch,
                         MDL="Pella_m",  #Pella_m if estimating m; other options: "Schaefer", "Fox", "Pella"
                         Ktch.CV=Ktch.CV,
                         ASS=Neim,
                         Rdist = Rdist,
                         Rprior = r.prior,
                         r.CV.multiplier=r.CV.multi,
                         Kdist= Kdist,
                         Kprior=K.prior,    
                         PsiDist= PsiDist,
                         Psiprior=psi.prior,   
                         Bprior=b.prior,    
                         BMSYK=bmsyk.mean,
                         shape.CV=r.prior[2])  #use same CV from r as shape parameter comes from there  
              
              #Run model    
              THIN=5
              CHAINS=2
              BURNIN=min(0.15*Scens$Sims[s],5000)
              JABBA.run=apply.JABBA(Ktch=input$Ktch,
                                 CPUE=Cpues1,
                                 CPUE.SE=Cpues1.SE,
                                 auxil=auxiliary.effort,
                                 auxil.se=auxiliary.effort.se,
                                 auxil.type=auxiliary.type,
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
                                 Shape.CV=input$shape.CV,
                                 output.dir=this.wd,
                                 outfile="fit_Appendix",
                                 Sims=Scens$Sims[s],
                                 Proc.error.JABBA=Proc.error,
                                 Obs.Error=Obs.Err.JABBA,
                                 thinning = THIN,
                                 nchains = CHAINS,
                                 burn.in= BURNIN)
              output=JABBA.run$fit
              
              #Display model fits
              fn.fig(paste(this.wd,"Fit_stdresiduals",sep='/'), 2400, 2400)
              jbplot_stdresiduals(output)
              dev.off()
              fn.fig(paste(this.wd,"Fit_dresiduals",sep='/'), 2400, 2400)
              jbplot_residuals(output)
              dev.off()
              fn.fig(paste(this.wd,"Fit_cpues",sep='/'), 2400, 2400)
              jbplot_cpuefits(output)
              dev.off()
              fn.fig(paste(this.wd,"Fit_cpues_log",sep='/'), 2400, 2400)
              jbplot_logfits(output) 
              dev.off()
              fn.fig(paste(this.wd,"Fit_catcherror",sep='/'), 2400, 2400)
              jbplot_catcherror(output)
              dev.off()
              fn.fig(paste(this.wd,"Summary",sep='/'), 2400, 2400)
              jbplot_summary(output)
              dev.off()
              
              #Get posteriors
              output$posteriors$P=output$posteriors$P[,1:nrow(ktch)]
              output$posteriors$BtoBmsy=output$posteriors$BtoBmsy[,1:nrow(ktch)]
              output$posteriors$SB=output$posteriors$SB[,1:nrow(ktch)]
              Store.sens[[s]]=list(input=input,output=output)
              
              #find pattern in residuals
              RES=as.data.frame(t(output$residuals))%>%
                mutate(across(everything(),~ifelse(.x<0,-1,1)),
                       across(everything(),~as.factor(.x)))
              store.runs.test=vector('list',ncol(RES))
              names(store.runs.test)=names(RES)
              for(x in 1:ncol(RES)) store.runs.test[[x]]=runs.test(RES[!is.na(RES[,x]),x])
              
              
              
              #Store Scenarios
              Out.Scens$Bo.mean[s]=Bint.mean
              Out.Scens$Bo.CV[s]=Bint.CV
              Out.Scens$r.mean[s]=Mn
              Out.Scens$r.cv[s]=Scens$r.sd[s]/Mn
              Out.Scens$Rdist[s]=Rdist
              Out.Scens$Kdist[s]=Kdist
              Out.Scens$PsiDist[s]=PsiDist
              Out.Scens$bmsyk.mean[s]=bmsyk.mean
              
              #Store estimates  
              d1=Store.sens[[s]]$output$estimates
              d1=d1[,grep(paste(c('mu','lci','uci'),collapse='|'),colnames(d1))]%>%
                data.frame
              colnames(d1)=c("Median","Lower.95","Upper.95")
              d1=d1%>%
                mutate(Parameter=rownames(d1),
                       Model='JABBA CPUE',
                       Scenario=names(Store.sens)[s])%>%
                relocate(Model,Scenario,Parameter,Lower.95)%>%
                filter(Parameter%in%c("K","r","psi","Hmsy","SBmsy","MSY"))   
              Out.estimates[[s]]=d1
              
              #Store trajectories  
              dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                mods=names(State.Space.SPM)[w],  
                                                Type='Depletion',
                                                scen=Scens$Scenario[s],
                                                Katch=ktch$Total)
              Out.rel.biom[[s]]=dummy$Dat
              Out.probs.rel.biom[[s]]=dummy$Probs
              
              dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                mods=names(State.Space.SPM)[w],
                                                Type='F.series',
                                                scen=Scens$Scenario[s],
                                                Katch=ktch$Total)
              Out.f.series[[s]]=dummy$Dat
              
              dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                mods=names(State.Space.SPM)[w],
                                                Type='B.Bmsy',
                                                scen=Scens$Scenario[s],
                                                Katch=ktch$Total)
              Out.B.Bmsy[[s]]=dummy$Dat
              
              dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
                                                mods=names(State.Space.SPM)[w],
                                                Type='F.Fmsy',
                                                scen=Scens$Scenario[s],
                                                Katch=ktch$Total)
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
                  gather(Parameter,value,-Iteration,- Chain)
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
                             capitalize(Neim),"/",AssessYr,"/",'JABBA CPUE',
                             "/Convergence diagnostics.tiff",sep=''),
                       width = 16,height = 8, dpi = 300, compression = "lzw")
                
                
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
                if(dummy.prior$psi$dist=='lnorm') dummy.prior$psi$mean=log(dummy.prior$psi$mean)
                dummy.post=Store.sens[[s]]$output$pars_posterior%>%
                  dplyr::select(r,K,psi)
                
                names(dummy.prior)=names(dummy.post)
                
                pars=colnames(dummy.post)
                out=vector('list',length(pars))
                
                
                for(p in 1:length(pars))
                {
                  Value.prior=fn.prior(d=dummy.prior[[pars[p]]])
                  if(pars[p]=="K") Value.prior=fn.prior(d=dummy.prior[[pars[p]]],MAX=max(stuff$Ktch$Total,na.rm=T)*50)
                  out[[p]]=rbind(data.frame(Distribuion="Prior",
                                            Value=Value.prior),
                                 data.frame(Distribuion="Posterior",
                                            Value=dummy.post[[pars[p]]]))%>%
                    mutate(Parameter=pars[p])
                }
                rm(stuff,dummy.prior,dummy.post)
                
                fn.show.density(d=do.call(rbind,out),NCOL=1)
                ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                             capitalize(Neim),"/",AssessYr,"/",'JABBA CPUE/',Scens$Scenario[s],"/Prior.and.posterior.tiff",sep=''),
                       width = 12,height = 14, dpi = 300, compression = "lzw")
              }
              
            }
            
            Out.Scens=Out.Scens%>%
              dplyr::select(Species,Scenario,
                            Rdist,r.mean,r.cv,
                            Kdist,K.mean,K.CV,
                            PsiDist,Bo.mean,Bo.CV,
                            bmsyk.mean,Proc.error,Ktch.CV,
                            Daily.cpues)%>%
              mutate(K.mean=round(K.mean),
                     r.mean=round(r.mean,3),
                     bmsyk.mean=round(bmsyk.mean,3),
                     r.cv=round(r.cv,2),
                     Bo.mean=round(Bo.mean,2),
                     Bo.CV=round(Bo.CV,2),
                     Rdist=ifelse(Rdist=='lnorm','Lognormal',Rdist),
                     Kdist=ifelse(Kdist=='lnorm','Lognormal',Kdist),
                     PsiDist=ifelse(PsiDist=='lnorm','Lognormal',PsiDist))

            dummy.store.sens.table[[i]]=Out.Scens
            dummy.store.rel.biom[[i]]=do.call(rbind,Out.rel.biom)
            dummy.store.probs.rel.biom[[i]]=Out.probs.rel.biom
            dummy.store.f.series[[i]]=do.call(rbind,Out.f.series)
            dummy.store.B.Bmsy[[i]]=do.call(rbind,Out.B.Bmsy)
            dummy.store.F.Fmsy[[i]]=do.call(rbind,Out.F.Fmsy)
            dummy.store.Kobe.probs[[i]]=Out.Kobe.probs  
            dummy.store.estimates[[i]]=do.call(rbind,Out.estimates)
            
            rm(Out.Scens,Out.rel.biom,Out.probs.rel.biom,Out.f.series,
               Out.B.Bmsy,Out.F.Fmsy,Out.Kobe.probs,Out.estimates)
          }
        }
      }
      
      State.Space.SPM[[w]]$sens.table=dummy.store.sens.table
      State.Space.SPM[[w]]$estimates=dummy.store.estimates
      State.Space.SPM[[w]]$rel.biom=dummy.store.rel.biom
      State.Space.SPM[[w]]$probs.rel.biom=dummy.store.probs.rel.biom
      State.Space.SPM[[w]]$f.series=dummy.store.f.series
      State.Space.SPM[[w]]$B.Bmsy=dummy.store.B.Bmsy
      State.Space.SPM[[w]]$F.Fmsy=dummy.store.F.Fmsy
      State.Space.SPM[[w]]$Kobe.probs=dummy.store.Kobe.probs
      
      rm(dummy.store,dummy.store.sens.table,dummy.store.estimates,
         dummy.store.rel.biom,dummy.store.probs.rel.biom,dummy.store.f.series,
         dummy.store.B.Bmsy,dummy.store.F.Fmsy)
      
    }
  }
}
toc(log = TRUE, quiet = TRUE)
computation.time <- tic.log(format = TRUE)
tic.clearlog()
send.email(TO=Send.email.to,
           CC='',
           Subject=paste("JABBA models finished running at",Sys.time()),
           Body= paste("Computation",computation.time),  
           Attachment=NULL) 


#---Generate outputs -------------------------------------------------

#24.1 Table of scenarios
for(l in 1: length(Lista.sp.outputs))  
{
  for(w in 1:length(State.Space.SPM))
  {
    dummy=State.Space.SPM[[w]]$sens.table
    dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
    
    write.csv(do.call(rbind,dummy)%>%relocate(Species),
              paste(Rar.path,paste('Table 6. ',names(State.Space.SPM)[w],' CPUE_scenarios_',
                                   names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
              row.names = F)
  }
}


#24.2 Table of parameter estimates by species and SPM method
for(l in 1: length(Lista.sp.outputs))
{
  for(w in 1:length(State.Space.SPM))
  {
    dummy=State.Space.SPM[[w]]$estimates
    dummy=dummy[match(Lista.sp.outputs[[l]],names(dummy))]
    dummy=compact(dummy)
    if(length(dummy)>0)
    {
      dummy=do.call(rbind,dummy)%>%
        rownames_to_column(var = "Species")%>%
        mutate(Species=capitalize(str_extract(Species, "[^.]+")))%>%
        relocate(Species)
      write.csv(dummy,paste(Rar.path,paste('Table 7. ',names(State.Space.SPM)[w],' CPUE_estimates_',
                                           names(Lista.sp.outputs)[l],'.csv',sep=''),sep='/'),
                row.names = F)
    }
    
  }
}


#24.3. Time series  
YLAB=XLAB=''

#24.3.1 Display all scenarios for each species
  #24.3.1.1 Relative biomass (i.e. Depletion)
Ref.points=vector('list',N.sp)
names(Ref.points)=Keep.species
for(i in 1:N.sp)
{
  print(paste("JABBA --- Relative biomass plot -----",Keep.species[i])) 
  a=fn.plot.timeseries(d=State.Space.SPM,
                       sp=Keep.species[i],
                       Type='Depletion',
                       YLAB='Relative biomass')
  if(!is.null(a))
  {
    #export graph
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/JABBA CPUE/JABBA_CPUE_time_series_relative_biomass.tiff",sep=''),
           width = 6,height = 10,compression = "lzw")
    
    #export current depletion probabilities
    write.csv(a$store.probs%>%
                spread(Model,Probability)%>%
                mutate(Species=Keep.species[i],
                       Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                arrange(Scenario,Range),
              paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(Keep.species[i]),"/",AssessYr,"/JABBA CPUE/JABBA_CPUE_current_depletion.csv",sep=''),
              row.names = F)
    Ref.points[[i]]=a$Ref.points$JABBA
  }
}
  #24.3.1.2 Fishing mortality
if(do.F.series)
{
  for(i in 1:N.sp)
  {
    print(paste("JABBA --- Fishing mortality plot -----",Keep.species[i]))
    a=fn.plot.timeseries(d=State.Space.SPM,
                                   sp=Keep.species[i],
                                   Type='F.series',
                                   YLAB='Fishing mortality')
    if(!is.null(a))
    {
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(Keep.species[i]),"/",AssessYr,"/JABBA CPUE/time_series_fishing_mortality.tiff",sep=''),
             width = 6,height = 10,compression = "lzw")
    }
  }
}
  #24.3.1.3 B over Bmsy
if(do.B.over.Bmsy.series)
{
  for(i in 1:N.sp)
  {
    print(paste("JABBA --- B over Bmsy plot -----",Keep.species[i]))
    a=fn.plot.timeseries(d=State.Space.SPM,
                                   sp=Keep.species[i],
                                   Type='B.Bmsy',
                                   YLAB='B/Bmsy')
    if(!is.null(a))
    {
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(Keep.species[i]),"/",AssessYr,"/JABBA CPUE/time_series_B_Bmsy.tiff",sep=''),
             width = 6,height = 10,compression = "lzw")
    }
  }
}
  #24.3.1.4 F over Fmsy
if(do.F.over.Fmsy.series)
{
  for(i in 1:N.sp)
  {
    print(paste("JABBA --- F over Fmsy plot -----",Keep.species[i]))
    a=fn.plot.timeseries(d=State.Space.SPM,
                                   sp=Keep.species[i],
                                   Type='F.Fmsy',
                                   YLAB='F/Fmsy')
    if(!is.null(a))
    {
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(Keep.species[i]),"/",AssessYr,"/JABBA CPUE/time_series_F_Fmsy.tiff",sep=''),
             width = 6,height = 10,compression = "lzw")
    }
  }
}

  #24.3.2 Display Scenario 1 for combined species

#Relative biomass (i.e. Depletion) 
  #figure
for(l in 1:length(Lista.sp.outputs))
{
  if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
  a=fn.plot.timeseries_combined(this.sp=Lista.sp.outputs[[l]],
                                d=State.Space.SPM$JABBA,
                                YLAB="Relative biomass",
                                Type="Depletion",
                                InnerMargin=InMar,
                                RefPoint=Ref.points,
                                Kach=State.Space.SPM$JABBA$rel.biom)
  WIDt=10
  if(length(compact(State.Space.SPM$JABBA$sens.table))<=3) WIDt=7
  if(!is.null(a))ggsave(paste(Rar.path,'/Relative.biomass_JABBA CPUE_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
                        width = WIDt,height = 10,compression = "lzw")
}
  #table 
for(l in 1:length(Lista.sp.outputs))
{
  dummy.mod=vector('list',length(State.Space.SPM))
  for(m in 1:length(State.Space.SPM))
  {
    str.prob=State.Space.SPM[[m]]$probs.rel.biom
    str.prob=str.prob[match(Lista.sp.outputs[[l]],names(str.prob))]
    str.prob=compact(str.prob)
    if(length(str.prob)>0)
    {
      dummy=vector('list',length =length(str.prob))
      for(d in 1:length(dummy))
      {
        dummy[[d]]=str.prob[[d]][[1]]$probs%>%
          mutate(Species=capitalize(names(str.prob)[d]))
      }
      dummy.mod[[m]]=do.call(rbind,dummy)%>%
        mutate(Model=names(State.Space.SPM)[m])
    }
  }
  dummy.mod=compact(dummy.mod)
  if(length(dummy.mod)>0)
  {
    write.csv(do.call(rbind,dummy.mod)%>%
                mutate(Range=factor(Range,levels=c("<lim","lim.thr","thr.tar",">tar")))%>%
                spread(Species,Probability)%>%
                arrange(Range),
              paste(Rar.path,'/Table 8. JABBA CPUE_Current.depletion_',names(Lista.sp.outputs)[l],'.csv',sep=''),
              row.names=F)
    rm(dummy.mod)
    
  }
}

  #24.3.3 Display sensitivity tests for combined species
for(l in 1:length(Lista.sp.outputs))
{
  if(length(Lista.sp.outputs[[l]])>8) InMar=1.25 else InMar=.5
  a=fn.plot.timeseries_combined_sensitivity(this.sp=Lista.sp.outputs[[l]], 
                                            d=State.Space.SPM$JABBA,
                                            InnerMargin=InMar,
                                            RefPoint=Ref.points,
                                            Kach=State.Space.SPM$JABBA$rel.biom)
  WIDt=10
  if(length(compact(State.Space.SPM$JABBA$sens.table))<=3) WIDt=7
  if(!is.null(a))ggsave(paste(Rar.path,'/Relative.biomass_JABBA CPUE_',names(Lista.sp.outputs)[l],'_sensitivity.tiff',sep=''),
                        width = WIDt,height = 10,compression = "lzw")
}

#24.4. Kobe plots (Scenario 1)  
  #24.4.1 by species 
store.kobes=vector('list',N.sp)
names(store.kobes)=Keep.species
for(i in 1:N.sp)
{
  if(!is.null(State.Space.SPM$JABBA$estimates[[i]]))
  {
    print(paste("JABBA --- Kobe plot -----",Keep.species[i]))
    store.kobes[[i]]=fn.get.Kobe.plot_appendix(d=State.Space.SPM,
                                               sp=Keep.species[i],
                                               do.probs=TRUE)
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(Keep.species[i]),"/",AssessYr,"/JABBA CPUE/Kobe_plot.tiff",sep=''),
           width = 10,height = 10, dpi = 300,compression = "lzw")
  }
  
}
store.kobes=compact(store.kobes)
  #24.4.2 Display combined species  
for(l in 1:length(Lista.sp.outputs))
{
  Nms=names(compact(State.Space.SPM$JABBA$estimates))
  this.sp=Lista.sp.outputs[[l]]
  this.sp=this.sp[which(this.sp%in%Nms)]
  if(length(this.sp)>0)
  {
    print(paste("JABBA --- Kobe plot -----",names(Lista.sp.outputs)[l]))
    DIMS=n2mfrow(length(this.sp))
    NKOL=DIMS[2]
    NRW=DIMS[1]
    if(NKOL%in%3:4) WIZ=14
    if(NKOL==2) WIZ=11
    if(NKOL==1) WIZ=9
    fn.get.Kobe.plot(this.sp,
                     d=State.Space.SPM$JABBA,
                     NKOL,
                     NRW,
                     do.probs=TRUE)
    ggsave(paste(Rar.path,'/Kobe_plot_JABBA CPUE_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
           width = WIZ,height = 12,compression = "lzw")
    
  }
}

#24.5 store Consequence and likelihood for WoE
Store.cons.Like_JABBA=fn.get.cons.like(lista=State.Space.SPM) 



