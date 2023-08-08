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

tic()
for(w in 1:length(State.Space.SPM))
{
  # JABBA (Winker et al 2018)   
  #summary of method: https://github.com/jabbamodel/JABBA
  if(names(State.Space.SPM)[w]=="JABBA")  #113 sec per species-scenario
  {
    do.parallel.JABBA=FALSE
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
          
          len.cpue=length(CPUE)
          MAX.CV=List.sp[[i]]$MAX.CV
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
            if(Neim=='gummy shark')     
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
            Out.Scens=Scens
            Out.estimates=Out.rel.biom=Out.probs.rel.biom=Out.f.series=
              Out.B.Bmsy=Out.F.Fmsy=vector('list',length(Store.sens))
            
            s=1
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
              b.prior=c(FALSE, 0.3, NA, "bk")  #only use b.prior for catch-only approach
              

              
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
                         MDL="Schaefer",
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
                         BMSYK=bmsyk.mean)  
              
              #Run model    
              THIN=5
              CHAINS=2
              BURNIN=min(0.15*Scens$Sims[s],5000)
              output=apply.JABBA(Ktch=input$Ktch,
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
                                 output.dir=this.wd,
                                 outfile="Appendix_fit",
                                 Sims=Scens$Sims[s],
                                 Proc.error.JABBA=Proc.error,
                                 Obs.Error=Obs.Err.JABBA,
                                 thinning = THIN,
                                 nchains = CHAINS,
                                 burn.in= BURNIN)
              
              #Display model fits
              fn.fig(paste(this.wd,"Fit_dresiduals",sep='/'), 2400, 2400)
              jbplot_residuals(output)
              dev.off()
              fn.fig(paste(this.wd,"Fit_cpues",sep='/'), 2400, 2400)
              jbplot_cpuefits(output)
              dev.off()
              fn.fig(paste(this.wd,"Fit_cpues_log",sep='/'), 2400, 2400)
              jbplot_logfits(output) 
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
              # RES=as.data.frame(t(output$residuals))%>%
              #   mutate(across(everything(),~ifelse(.x<0,-1,1)),
              #          across(everything(),~as.factor(.x)))
              # store.runs.test=vector('list',ncol(RES))
              # names(store.runs.test)=names(RES)
              # for(x in 1:ncol(RES)) store.runs.test[[x]]=runs.test(RES[!is.na(RES[,x]),x])
              # 
              
              
              #Store Scenarios
              # Out.Scens$Bo.mean=Bint.mean
              # Out.Scens$Bo.CV=Bint.CV
              # Out.Scens$r.mean=Mn
              # Out.Scens$r.cv=Scens$r.sd[s]/Mn
              # Out.Scens$Rdist=Rdist
              # Out.Scens$Kdist=Kdist
              # Out.Scens$PsiDist=PsiDist
              # Out.Scens$bmsyk.mean=Out.Scens$bmsyk
              
              #Store estimates  
              # d1=Store.sens[[s]]$output$estimates
              # d1=d1[,grep(paste(c('mu','lci','uci'),collapse='|'),colnames(d1))]%>%
              #   data.frame
              # colnames(d1)=c("Median","Lower.95","Upper.95")
              # d1=d1%>%
              #   mutate(Parameter=rownames(d1),
              #          Model='JABBA CPUE',
              #          Scenario=names(Store.sens)[s])%>%
              #   relocate(Model,Scenario,Parameter,Lower.95)%>%
              #   filter(Parameter%in%c("K","r","psi","Hmsy","SBmsy","MSY"))   
              # Out.estimates[[s]]=d1
              
              #Store trajectories  
              # dummy=fn.ktch.only.get.timeseries(d=Store.sens[[s]],
              #                                   mods=names(State.Space.SPM)[w],  
              #                                   Type='Depletion',
              #                                   scen=Scens$Scenario[s],
              #                                   Katch=ktch$Total)
              # plot(dummy$Dat$year,dummy$Dat$median,ylim=c(0,1))
              # lines(dummy$Dat$year,dummy$Dat$upper.95)
              # lines(dummy$Dat$year,dummy$Dat$lower.95)
              # 
  

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
                  if(pars[p]=="K") Value.prior=fn.prior(d=dummy.prior[[pars[p]]],MAX= max(stuff$Ktch$Total,na.rm=T)*50)
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
                ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                             capitalize(Neim),"/",AssessYr,"/",'JABBA CPUE/',Scens$Scenario[s],"/Observed_predicted_cpue.tiff",sep=''),
                       width = 6,height = 6, dpi = 300, compression = "lzw")
                
                
            rm(output)

        }
      }
      
    }
  }
}
toc()

