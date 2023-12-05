#---15. Export .pin files and define modelling arguments  -----
#note: For integrated model, all pars calculated in log space.
#       ln_RZERO is in 1,000 individuals so do 10 times the largest catch divided by 
#       average weight and divided by 1,000. Best units to work in are 1,000 individuals for 
#       numbers, catch in tonnes and length-weight in kg as all cancels out and predicted biomass
#       end up being in tonnes

k.fun=function(KTCH,times.ktch) max(KTCH,na.rm=T)*times.ktch  #use species-specific min and max times k based on exploitation history
                                                               #upper bound: 50 times (Andrade 2017)
fn.mtch=function(WHAT,NMS) match(WHAT,names(NMS))
if(Do.bespoked)
{
  Q_phz=c("lnq","lnq2","log_Qdaily")                           
  Zns.par.phz=c("lnR_prop_west","lnR_prop_zn1")
  MOv.par.phz=c("log_p11","log_p22","log_p21","log_p33")
}

#Recruitment inputs
SS3.Rrecruitment.inputs=read.csv(handl_OneDrive('Analyses/Population dynamics/SS3.Rrecruitment.inputs.csv'))
SS3.tune_size_comp_effective_sample=read.csv(handl_OneDrive('Analyses/Population dynamics/SS3.tune_size_comp_effective_sample.csv'))

#Final depletion
#note: modify depletion levels for species with negligible catches in recent years (based on generation time and r)
Depletion.levels<- read_excel(handl_OneDrive('Analyses/Population dynamics/K_&_depletion.levels.xlsx'),sheet = "K_&_depletion.levels")
for(i in 1:N.sp)
{
  G=round(store.species.G_M.min[[i]]$mean) #generation time
  R=store.species.r_M.min[[i]]$mean        #r
  ktch=ktch.combined%>%
    filter(Name==Keep.species[i])%>%
    mutate(Rel.ktch=Tonnes/max(Tonnes),
           yr.index=row_number(),
           rel.ktch.id=ifelse(Rel.ktch>0.15,0,1))
  id.yr=which(ktch$Rel.ktch==1)
  id.G=(nrow(ktch)-G+1):nrow(ktch)
  ktch=ktch%>%mutate(G.yr=ifelse(yr.index%in%id.G,1,0))
  
  ktch=ktch[id.yr:nrow(ktch),]%>%
    filter(rel.ktch.id==1)%>%
    mutate(Depletion=0.2)
  if(nrow(ktch)>0)
  {
    alpha=R*.85
    if(nrow(ktch)>1)for(q in 2:nrow(ktch)) ktch$Depletion[q]=(ktch$Depletion[q-1]+ktch$Depletion[q-1]*R)*exp(-alpha*ktch$Depletion[q-1])
    ktch=ktch%>%filter(G.yr==1)
    
    if(tweak.Final.Bio_low)
    {
      if(ktch$Depletion[nrow(ktch)]>=0.8) updated.FINALBIO1=0.6
      if(ktch$Depletion[nrow(ktch)]>=0.6 & ktch$Depletion[nrow(ktch)]<0.8) updated.FINALBIO1=0.4
      if(ktch$Depletion[nrow(ktch)]<0.6)  updated.FINALBIO1=Depletion.levels%>%filter(Species==Keep.species[i])%>%pull(FINALBIO1)
      if(Depletion.levels$Species[i]=="narrow sawfish" & updated.FINALBIO1==0.4) updated.FINALBIO1=0.3
      if(Depletion.levels$Species[i]=="scalloped hammerhead" & updated.FINALBIO1==0.4) updated.FINALBIO1=0.2
      if(Depletion.levels$Species[i]%in%c("green sawfish","snaggletooth","weasel shark","tiger shark","zebra shark") & updated.FINALBIO1==0.6) updated.FINALBIO1=0.3
      Depletion.levels[i,]=Depletion.levels[i,]%>%
        mutate(FINALBIO1=ifelse(Species==Keep.species[i],updated.FINALBIO1,
                                FINALBIO1))
    }
  }
}

#Create pin files
for(l in 1:N.sp)    
{
  NeiM=names(List.sp)[l]
  
  print(paste("---------Modelling arguments for --",NeiM))
  
  LH=LH.data%>%filter(SPECIES==List.sp[[l]]$Species)
  
  # 1...General arguments
  List.sp[[l]]$B.init=1         #initial depletion
  List.sp[[l]]$MAX.CV=0.5       #maximum acceptable CV for cpue series
  List.sp[[l]]$Lzero_SD=5
  
  
    #steepness
  h.min=round(min(store.species.steepness.S2[[l]],store.species.steepness_M.mean[[l]]$mean),3)
  h.M_mean2=round(max(Min.h.shark,store.species.steepness.S2[[l]]),3)
  h.M_mean=round(store.species.steepness_M.mean[[l]]$mean,3)
  if(NeiM%in%h_too.long.converge) h.M_mean=h.M_mean2
  h.M.mean_low=round(max(Min.h.shark,h.M_mean*.9),3)
  h.M_mean=round(max(Min.h.shark,h.M_mean),3)
  
    #Effective sample size (Francis method as default)
  List.sp[[l]]$tuned_size_comp=NULL
  tuned_size_comp=SS3.tune_size_comp_effective_sample%>%
                            filter(Species==NeiM)%>%
                            dplyr::select(-Species)%>%
                            gather(Fleet,Value,-Data_type)%>%
                            filter(!is.na(Value))%>%
                            mutate(Fleet=as.numeric(str_extract(Fleet, "[[:digit:]]+")))%>%
                            rename(Factor=Data_type)
  if(nrow(tuned_size_comp)>0) List.sp[[l]]$tuned_size_comp=tuned_size_comp
  
    #recruitment
  ramp.yrs=SS3.Rrecruitment.inputs%>%filter(Species==NeiM)
  Min.logR0=ramp.yrs$Ln_R0_min
  Upper.bound.LnR0=ramp.yrs$Ln_R0_max
  if(is.na(Upper.bound.LnR0))
  {
    if(NeiM%in%names(Indicator.species)) Upper.bound.LnR0=Upper.bound.LnR0*1.3
    
    #tweak to improve SSS acceptance rate
    if(NeiM%in%c("great hammerhead","shortfin mako","weasel shark","angel sharks",
                 "sawsharks","wobbegongs","spurdogs","zebra shark")) Upper.bound.LnR0=7 
    if(NeiM%in%c("green sawfish","narrow sawfish","pigeye shark","scalloped hammerhead")) Upper.bound.LnR0=6
    if(NeiM%in%c("grey nurse shark","lemon shark","snaggletooth")) Upper.bound.LnR0=5
  }

  
    #M
  #age-invariant
  Mmean=mean(unlist(store.species.M_M.min[[l]]),na.rm=T)
  Mmean.mean=mean(unlist(store.species.M_M.mean[[l]]),na.rm=T)
  
      #increase M (based on Hoenig's maximum age) a bit to allow SS3 to fit, too low M and too high h otherwise
  if(NeiM=="lemon shark")
  {
    Mmean=0.154   
    Mmean.mean=Mmean.mean*1.5
  }
  if(NeiM=="grey nurse shark") Mmean.mean=0.102
  if(NeiM=="spurdogs") Mmean.mean=0.149
  if(NeiM=="sawsharks") Mmean.mean=0.18
  if(NeiM%in%c("grey nurse shark","shortfin mako"))
  {
    Mmean=Mmean*2
    Mmean.mean=Mmean.mean*2
  }
  Msd=max(sd(unlist(store.species.M_M.min[[l]]),na.rm=T),.01)
  
  #age-variant
  M_y=apply(store.species.M_M.min[[l]],2,mean,na.rm=T)
  M_x=1:length(M_y)
  loe.mod=loess(M_y~M_x)
  List.sp[[l]]$Mmean.min.at.age=predict(loe.mod)
  
  M_y=apply(store.species.M_M.mean[[l]],2,mean,na.rm=T)
  M_x=1:length(M_y)
  loe.mod=loess(M_y~M_x)
  List.sp[[l]]$Mmean.mean.at.age=predict(loe.mod)
  if(Use.SEDAR.M)
  {
    if(NeiM=="sandbar shark")
    {
      #Age-dependent M methods yield very high M given new k of 0.12 so set to the one use in SEDAR 54
      List.sp[[l]]$Mmean.mean.at.age=c(rep(0.1604,6),0.1578,rep(0.1168,length(8:length(List.sp[[l]]$Mmean.mean.at.age))))
    }
    
    if(NeiM=="dusky shark")
    {
      List.sp[[l]]$Mmean.mean.at.age=c(rep(0.104,5),0.098,0.092,0.088,0.084,0.08,0.077,0.074,
                                       0.072,0.07,0.068,0.066,0.064,0.063,0.061,0.06,0.059,
                                       0.058,0.057,0.056,0.055,0.054,0.053,0.052,0.052,0.051,
                                       rep(0.048,length(31:length(List.sp[[l]]$Mmean.mean.at.age))))
    }
    
  }

  
  # 2...Catch-only arguments and sensitivity tests
  Depl.lvl=Depletion.levels%>%filter(Species==NeiM)
  
  List.sp[[l]]$STARTBIO=c(List.sp[[l]]$B.init*.95,List.sp[[l]]$B.init)   #low initial depletion because starting time series prior to any fishing
  FINALBIO=c(Depl.lvl$FINALBIO1,Depl.lvl$FINALBIO2)       #uncertain though ranges informed by fitting SSS and excluding B.final cases with poor fit
  List.sp[[l]]$FINALBIO=FINALBIO
  List.sp[[l]]$Do.sim.test=Do.sim.Test  #simulation test Catch-MSY for small and large catches
  bmsyK=0.5   #Schaefer
  if(NeiM=="lemon shark") bmsyK=0.45    #to improve convergence of DBSRA
  BmsyK.Cortes=BmsyK.species%>%filter(Species==names(List.sp)[l])%>%pull(BmsyK)
  if(tweak.BmsyK.Cortes)
  {
    if(NeiM%in%c("angel sharks","lemon shark")) BmsyK.Cortes=0.4   
    if(NeiM=="grey nurse shark") BmsyK.Cortes=0.49
    if(NeiM=="milk shark") BmsyK.Cortes=0.375
    if(NeiM=="pigeye shark") BmsyK.Cortes=0.45
    if(NeiM=="sawsharks") BmsyK.Cortes=0.475
    if(NeiM=="sandbar shark") BmsyK.Cortes=0.49
  }
  
  #2.1 DBSRA scenarios
  AgeMat=List.sp[[l]]$Age.50.mat[1]

  ktch=ktch.combined%>%
    filter(Name==List.sp[[l]]$Name)  
  Klow=k.fun(ktch$Tonnes,Depl.lvl%>%pull(K_times_max_ktch_min)) 
  Kup=max(K_min,k.fun(ktch$Tonnes,Depl.lvl%>%pull(K_times_max_ktch_max)))  
  fmsy.m=0.41                 # Zhou et al 2012 fmsy/M=0.41 for chondrichthyans
  fmsy.m2=Fmsy.M.scaler[[l]]  #Cortes & Brooks 2018
  
  
  #2.2 CMSY scenarios
  r.mean=store.species.r_M.min[[l]]$mean
  r.mean2=store.species.r_M.mean[[l]]$mean
  if(NeiM%in%c("dusky shark","sandbar shark","grey nurse shark","spurdogs"))
  {
    r.mean2=r.mean*0.6  #r.mean2 couldn't be estimated due to life history mis-specification so 
                        # set at fraction of r.mean based on other species ratio
  }
  r.mean=round(r.mean,3)
  r.mean2=round(r.mean2,3)
  r.mean.sd=store.species.r_M.min[[l]]$sd
  r.mean.sd2=store.species.r_M.mean[[l]]$sd
  r.mean_lower=round(r.mean*.9,3)
  r.mean_upper=round(r.mean*1.2,3)
  ensims.csmy=ensims.CSMY
  Kup.spcfik=Kup
  Klow.spcfik=Klow
  Proc.error.cmsy=Proc.Error.1


  #2.3 JABBA - catch only scenarios
  r.CV.multiplier=1
  Ktch.CV=0.1   #Winker et al 2019 set it at 0.2 for uncertain reconstructed school shark catches
  
  #2.4 Scenarios 
  #combine all input pars and assumptions in list 
  if(Simplified.scenarios)
  {
    List.sp[[l]]$Sens.test=list(
      DBSRA=data.frame(Scenario=paste("S",1,sep=''),
                       Sims=ensims, 
                       AgeMat=AgeMat,
                       Mmean=Mmean,
                       Msd=Msd,
                       Klow=Klow,
                       Kup=Kup.spcfik,
                       fmsy.m=fmsy.m,    
                       bmsyk=BmsyK.Cortes),    
      CMSY=data.frame(Scenario=paste("S",1,sep=''),
                      Sims=ensims.csmy, 
                      r=r.mean,
                      r.sd=r.mean.sd,
                      r.prob.min=r.prob.min,  
                      r.prob.max=r.prob.max,
                      Klow=Klow.spcfik,
                      Kup=Kup.spcfik,
                      Proc.error=Proc.Error),
      JABBA_catch.only=data.frame(Scenario=paste("S",1:1,sep=''),
                                  Sims=ensims.JABBA,
                                  r=r.mean,
                                  r.sd=r.mean.sd,
                                  r.CV.multiplier=r.CV.multiplier,
                                  K.mean=mean(c(Klow.spcfik,Kup.spcfik)), 
                                  K.CV=k.cv,
                                  Proc.error=Proc.Error,
                                  Ktch.CV=Ktch.CV,
                                  bmsyk=BmsyK.Cortes),      
      SS3=data.frame(Scenario=paste('S',1:2,sep=''),
                     Model='SSS',
                     Mmean=Mmean.mean,  #using Mmean.mean as base case because h.M_mean2 is used as basecase (both using Mmean.mean)    
                     Msd=Msd,
                     Final.dpl=mean(FINALBIO),
                     Steepness=c(h.M_mean,h.M_mean2),
                     Steepness.sd=store.species.steepness_M.min[[l]]$sd,
                     F_ballpark=1e-3,
                     use_F_ballpark=FALSE,
                     Sims=SSS.sims,
                     Ln_R0_min=Min.logR0,
                     Ln_R0_max=Upper.bound.LnR0))
    if(List.sp[[l]]$Sens.test$SS3$Steepness[1]==List.sp[[l]]$Sens.test$SS3$Steepness[2])
    {
      List.sp[[l]]$Sens.test$SS3=List.sp[[l]]$Sens.test$SS3[1,]
    }
  }else
  {
    List.sp[[l]]$Sens.test=list(
      DBSRA=data.frame(Scenario=paste("S",1:4,sep=''),
                       Sims=rep(ensims,4), 
                       AgeMat=rep(AgeMat,4),
                       Mmean=c(rep(Mmean,3),Mmean.mean),
                       Msd=rep(Msd,4),
                       Klow=rep(Klow,4),
                       Kup=rep(Kup.spcfik,4),
                       fmsy.m=c(fmsy.m,fmsy.m2,rep(fmsy.m,2)),    
                       bmsyk=c(rep(0.5,2),BmsyK.Cortes,0.5)),    #Winker et al 2020 Bmsy/K=0.55 for mako shark
      CMSY=data.frame(Scenario=paste("S",1:3,sep=''),
                      Sims=rep(ensims.csmy,3), 
                      r=c(rep(r.mean,2),r.mean2),
                      r.sd=c(rep(r.mean.sd,2),r.mean.sd2),
                      r.prob.min=rep(r.prob.min,3),  
                      r.prob.max=rep(r.prob.max,3),
                      Klow=rep(Klow.spcfik,3),
                      Kup=rep(Kup.spcfik,3),
                      Proc.error=c(Proc.Error,Proc.error.cmsy,Proc.Error)),
      JABBA_catch.only=data.frame(Scenario=paste("S",1:4,sep=''),
                                  Sims=rep(ensims.JABBA,4),
                                  r=c(r.mean,r.mean2,rep(r.mean,2)),
                                  r.sd=c(r.mean.sd,r.mean.sd2,rep(r.mean.sd,2)),
                                  r.CV.multiplier=rep(r.CV.multiplier,4),
                                  K.mean=rep(mean(c(Klow.spcfik,Kup.spcfik)),4), #Winker et al 2019 mean k set at 20 times max catch
                                  K.CV=rep(2,4),
                                  Proc.error=c(rep(Proc.Error,2),Proc.error.cmsy,Proc.Error),
                                  Ktch.CV=rep(Ktch.CV,4),
                                  bmsyk=c(rep(0.5,3),BmsyK.Cortes)),      #Winker et al 2020 Bmsy/K=0.55 for mako shark
      SS3=data.frame(Scenario=paste('S',1:4,sep=''),
                     Model=rep('SSS',4),
                     Mmean=c(rep(Mmean.mean,3),Mmean),      #using Mmean.mean as base case because h.M_mean2 is used as basecase (both using Mmean.mean)
                     Final.dpl=c(mean(FINALBIO)),
                     Steepness=c(h.M_mean,h.M_mean2,h.M.mean_low,h.M_mean),
                     Steepness.sd=rep(store.species.steepness_M.min[[l]]$sd,4),
                     F_ballpark=rep(1e-3,4),
                     use_F_ballpark=rep(FALSE,4),
                     Sims=rep(SSS.sims,4),
                     Ln_R0_min=rep(Min.logR0,4),
                     Ln_R0_max=rep(Upper.bound.LnR0,4))
    ) 
  }
  
  #2.5 JABBA - cpue scenarios
  List.sp[[l]]$Sens.test$JABBA=List.sp[[l]]$Sens.test$JABBA_catch.only[1,]%>%
                                        mutate(Proc.error=Proc.Error.cpue)   
  tested.r=c(r.mean,r.mean_lower,r.mean_upper,rep(r.mean,5))
  #tested.r=rep(r.mean,4)
  List.sp[[l]]$Sens.test$JABBA=do.call("rbind", replicate(length(tested.r), List.sp[[l]]$Sens.test$JABBA, simplify = FALSE))%>%
                                  mutate(r=tested.r,
                                         Scenario=paste('S',row_number(),sep=''),
                                         Klow=NA,
                                         Kup=NA)
  List.sp[[l]]$Sens.test$JABBA$Proc.error[4]=Proc.Error.cpue2
  List.sp[[l]]$Sens.test$JABBA$Kdist=KDIST
  List.sp[[l]]$Sens.test$JABBA$Kdist[5]="range"
  List.sp[[l]]$Sens.test$JABBA$K.mean[5]=NA  
  List.sp[[l]]$Sens.test$JABBA$K.CV[5]=NA
  List.sp[[l]]$Sens.test$JABBA$Klow[5]=Klow.spcfik
  List.sp[[l]]$Sens.test$JABBA$Kup[5]=Kup.spcfik
  List.sp[[l]]$Sens.test$JABBA$Daily.cpues=c(rep(drop.daily.cpue,length(tested.r)-1),NA)
  if(NeiM%in%c("gummy shark","whiskery shark")) List.sp[[l]]$Sens.test$JABBA$K.mean=6000  #improve K mean as very wide range used for catch only
  if(NeiM=="tiger shark") List.sp[[l]]$Sens.test$JABBA$K.mean=9000
  if(NeiM=="smooth hammerhead") List.sp[[l]]$Sens.test$JABBA$K.mean=1500
  List.sp[[l]]$Sens.test$JABBA$use.these.abundances=NA
  List.sp[[l]]$Sens.test$JABBA$model.type=c(rep("Pella_m",nrow(List.sp[[l]]$Sens.test$JABBA)-3),"Schaefer", "Fox","Pella_m")
  if(NeiM%in%c("dusky shark","sandbar shark"))
  {
    List.sp[[l]]$Sens.test$JABBA$use.these.abundances=c('Survey') # Due to selectivity, TDGDLF catches small juveniles not spawning stock "Survey_TDGDLF.monthly"
    List.sp[[l]]$Sens.test$JABBA=List.sp[[l]]$Sens.test$JABBA%>%   #redundant scenario
                                    filter(!is.na(Daily.cpues))
    List.sp[[l]]$Sens.test$JABBA$model.type=c(rep("Pella_m",nrow(List.sp[[l]]$Sens.test$JABBA)-2),"Schaefer", "Fox") 
  }
  if(!evaluate.07.08.cpue) List.sp[[l]]$Sens.test$JABBA=List.sp[[l]]$Sens.test$JABBA%>%filter(!is.na(Daily.cpues))
  
  
  # 3... Dynamic catch and size model
  InRec=NULL
  if(NeiM=="smooth hammerhead") InRec= 15000
  if(NeiM=="spinner shark") InRec= 15000
  if(NeiM=="dusky shark") InRec= 300000
  if(NeiM=="sandbar shark") InRec= 50000
  if(NeiM=="gummy shark") InRec= 200000
  if(NeiM=="whiskery shark") InRec= 300000
  List.sp[[l]]$InitRec_dynKtchSize =  InRec  
  
  
  # 4... Integrated models
  
  #-- 4.1 SS
  
  #4.1.1 Scenarios
  List.sp[[l]]$Sens.test$SS=List.sp[[l]]$Sens.test$SS3[1,]%>%
                              mutate(NSF.selectivity=NA)
  tested.h=c(List.sp[[l]]$Sens.test$SS3$Steepness,h.min,List.sp[[l]]$Sens.test$SS$Steepness[1])
  if(NeiM%in%h_too.high & !NeiM%in%h_too.long.converge)  tested.h=c(List.sp[[l]]$Sens.test$SS3$Steepness,List.sp[[l]]$Sens.test$SS$Steepness[1]) 
  if(NeiM%in%h_too.low)  tested.h=c(h.M_mean2,store.species.steepness_M.mean[[l]]$mean,h.M_mean2) 
  tested.h=unique(tested.h)
  List.sp[[l]]$Sens.test$SS=do.call("rbind", replicate(length(tested.h), List.sp[[l]]$Sens.test$SS, simplify = FALSE))%>%
                              mutate(Steepness=tested.h,
                                     Scenario=paste('S',row_number(),sep=''))
    #LnR0 init
  List.sp[[l]]$Sens.test$SS=List.sp[[l]]$Sens.test$SS%>%
                              mutate(Ln_R0_init=runif(1,Ln_R0_max*.3,Ln_R0_max*.6))
    #Ramp
  if(!is.na(ramp.yrs$LnRo))
  {
    List.sp[[l]]$Sens.test$SS=List.sp[[l]]$Sens.test$SS%>%
                                mutate(Ln_R0_init=ramp.yrs$LnRo,
                                       Ln_R0_max=ramp.yrs$Ln_R0_max)
  }
    #M at age
    List.sp[[l]]$Sens.test$SS$M.at.age=rep('Mmean.mean.at.age',nrow(List.sp[[l]]$Sens.test$SS)) #use Mean M for consistency with h

    #Cpues
    if(!evaluate.07.08.cpue) List.sp[[l]]$Sens.test$SS$Daily.cpues=rep(drop.daily.cpue,length(tested.h))
    if(evaluate.07.08.cpue) List.sp[[l]]$Sens.test$SS$Daily.cpues=c(rep(drop.daily.cpue,length(tested.h)-1),NA)
  
    #Always calculate extra SD for Q or only when CV is small
  List.sp[[l]]$Sens.test$SS$extra.SD.Q='always'
  
    #Remove SSS inputs
  List.sp[[l]]$Sens.test$SS=List.sp[[l]]$Sens.test$SS%>%
                                  mutate(Model='SS')%>%
                                  dplyr::select(-c(Final.dpl,Sims))
  
  if(!evaluate.07.08.cpue) List.sp[[l]]$Sens.test$SS=List.sp[[l]]$Sens.test$SS%>%filter(!is.na(Daily.cpues))
  
    #Alternative selectivity for NSF and survey
  if(NeiM%in%alternative.NSF.selectivity)
  {
    add.dumi=List.sp[[l]]$Sens.test$SS[1,]%>%
              mutate(NSF.selectivity='Logistic',
                     Scenario=paste0('S',nrow(List.sp[[l]]$Sens.test$SS)+1))
    List.sp[[l]]$Sens.test$SS=rbind(List.sp[[l]]$Sens.test$SS,add.dumi)
  }
  
    #4.1.2 Growth
  List.sp[[l]]$Growth.CV_young=0.1  #constant CV 8.5%-10% Tremblay-Boyer et al 2019
  List.sp[[l]]$Growth.CV_old=0.1
  if(NeiM%in%c("milk shark","wobbegongs")) List.sp[[l]]$Growth.CV_old=0.15  #widen the CV to allow observed max sizes be within Linf + CV
  List.sp[[l]]$SS3.estim.growth.pars=FALSE
  #if(NeiM%in%c("milk shark")) List.sp[[l]]$SS3.estim.growth.pars=TRUE
  List.sp[[l]]$compress.tail=FALSE  #compress the size distribution tail into plus group. 
    
    #Qs (in log space)
  List.sp[[l]]$Q.inits=data.frame(Fleet=c('Northern.shark','Other','Southern.shark_1','Southern.shark_2','Survey',
                                          'F.series_Southern.shark_1'),
                                  Fleet.n=1:6,
                                  Q=log(c(fn.ji(1e-2),fn.ji(1e-2),fn.ji(1e-3),fn.ji(4e-4),fn.ji(1e-2),fn.ji(1e-3))))

    #4.1.3 selectivity
  #SizeSelex
  #note: unknown for NSF and other fisheries so set to logistic using available size compo as a default unless otherwise
  #      double normal: p1. peak, beginning size for plateau
  #                     p2. top, width of plateau
  #                     p3. ascending width
  #                     p4. descending width
  #                     p5. selectivity at first bin
  #                     p6. selectivity at last bin (set to 1 for logistic)
  #       logistic: p1. inflection
  #                 p2. 95% width (>0)
  
  SS.sel.init.pars=SS_selectivity_init_pars%>%filter(Species==NeiM)
    #init values
  attach(SS.sel.init.pars)
      #NSF, Survey and Others
    p1.sel=NSF_p1       
    p2.sel=NSF_p2 
    p3.sel=NSF_p3
    p4.sel=NSF_p4 
    p5.sel=NSF_p5 
    p6.sel=NSF_p6 
    if(NSF_Alternative_sel_type=='Logistic')
    {
      p1.sel=Logistic_p1       
      p2.sel=Logistic_p2 
      p3.sel=NA
      p4.sel=NA 
      p5.sel=NA 
      p6.sel=NA 
    }
    if(NeiM%in%alternative.NSF.selectivity)
    {
      p1.sel_sens=Logistic_p1       
      p2.sel_sens=Logistic_p2 
      p3.sel_sens=NA
      p4.sel_sens=NA 
      p5.sel_sens=NA 
      p6.sel_sens=NA 
    }
    p1.sel_NSF=p1.sel_Other=p1.sel_Survey=p1.sel       
    p2.sel_NSF=p2.sel_Other=p2.sel_Survey=p2.sel 
    p3.sel_NSF=p3.sel_Other=p3.sel_Survey=p3.sel
    p4.sel_NSF=p4.sel_Other=p4.sel_Survey=p4.sel 
    p5.sel_NSF=p5.sel_Other=p5.sel_Survey=p5.sel 
    p6.sel_NSF=p6.sel_Other=p6.sel_Survey=p6.sel 

      #TDGLDF
    p1.sel_TDGDLF=TDGDLF_p1
    p2.sel_TDGDLF=TDGDLF_p2
    p3.sel_TDGDLF=TDGDLF_p3
    p4.sel_TDGDLF=TDGDLF_p4
    p5.sel_TDGDLF=TDGDLF_p5
    p6.sel_TDGDLF=TDGDLF_p6
  detach(SS.sel.init.pars)
  
  if(NeiM=="tiger shark")
  {
    p1.sel_Other=220   #Lobster fishery targeting large sharks but not jumbos!
    p2.sel_Other=-2.0
    p3.sel_Other=7
    p4.sel_Other=8.0
    p5.sel_Other=-999
    p6.sel_Other=-999
  }
    
    List.sp[[l]]$SS_selectivity=data.frame(Fleet=c("Northern.shark","Other","Southern.shark_1","Southern.shark_2","Survey"),
                                           P_1=c(p1.sel_NSF,p1.sel_Other,p1.sel_TDGDLF,p1.sel_TDGDLF,p1.sel_Survey),
                                           P_2=c(p2.sel_NSF,p2.sel_Other,p2.sel_TDGDLF,p2.sel_TDGDLF,p2.sel_Survey),
                                           P_3=c(p3.sel_NSF,p3.sel_Other,p3.sel_TDGDLF,p3.sel_TDGDLF,p3.sel_Survey),
                                           P_4=c(p4.sel_NSF,p4.sel_Other,p4.sel_TDGDLF,p4.sel_TDGDLF,p4.sel_Survey),
                                           P_5=c(p5.sel_NSF,p5.sel_Other,p5.sel_TDGDLF,p5.sel_TDGDLF,p5.sel_Survey),
                                           P_6=c(p6.sel_NSF,p6.sel_Other,p6.sel_TDGDLF,p6.sel_TDGDLF,p6.sel_Survey))
    if(NeiM%in%alternative.NSF.selectivity)
    {
      p1.sel_NSF=p1.sel_Other=p1.sel_Survey=p1.sel_sens       
      p2.sel_NSF=p2.sel_Other=p2.sel_Survey=p2.sel_sens 
      p3.sel_NSF=p3.sel_Other=p3.sel_Survey=p3.sel_sens
      p4.sel_NSF=p4.sel_Other=p4.sel_Survey=p4.sel_sens 
      p5.sel_NSF=p5.sel_Other=p5.sel_Survey=p5.sel_sens
      p6.sel_NSF=p6.sel_Other=p6.sel_Survey=p6.sel_sens
      List.sp[[l]]$SS_selectivity.sensitivity=data.frame(Fleet=c("Northern.shark","Other","Southern.shark_1","Southern.shark_2","Survey"),
                                             P_1=c(p1.sel_NSF,p1.sel_Other,p1.sel_TDGDLF,p1.sel_TDGDLF,p1.sel_Survey),
                                             P_2=c(p2.sel_NSF,p2.sel_Other,p2.sel_TDGDLF,p2.sel_TDGDLF,p2.sel_Survey),
                                             P_3=c(p3.sel_NSF,p3.sel_Other,p3.sel_TDGDLF,p3.sel_TDGDLF,p3.sel_Survey),
                                             P_4=c(p4.sel_NSF,p4.sel_Other,p4.sel_TDGDLF,p4.sel_TDGDLF,p4.sel_Survey),
                                             P_5=c(p5.sel_NSF,p5.sel_Other,p5.sel_TDGDLF,p5.sel_TDGDLF,p5.sel_Survey),
                                             P_6=c(p6.sel_NSF,p6.sel_Other,p6.sel_TDGDLF,p6.sel_TDGDLF,p6.sel_Survey))
      
    }
    #retention & discard mortality
    if(any(!is.na(SS.sel.init.pars[,grepl('Ret_',names(SS.sel.init.pars))])))
    {
      xx=SS.sel.init.pars[,grepl('Ret_',names(SS.sel.init.pars))]
      if(!is.na(xx$Ret_p1))
      {
        List.sp[[l]]$SS_retention=data.frame(Fleet=xx$Ret_Fleet,
                                             P_1=xx$Ret_p1,
                                             P_2=xx$Ret_p2,
                                             P_3=xx$Ret_p3,
                                             P_4=xx$Ret_p4,
                                             P_5=xx$Ret_p5,
                                             P_6=xx$Ret_p6,
                                             P_7=xx$Ret_p7)
      }
        
      yy=SS.sel.init.pars[,grepl('Disc_Fleet',names(SS.sel.init.pars))]
      if(!is.na(yy$Disc_Fleet_L1))
      {
        List.sp[[l]]$SS_discard_mortality=data.frame(Fleet=xx$Ret_Fleet,
                                                     P_1=yy$Disc_Fleet_L1,
                                                     P_2=yy$Disc_Fleet_L2,
                                                     P_3=yy$Disc_Fleet_L3,
                                                     P_4=yy$Disc_Fleet_L4)
      }
      
    }
    
    #phases
    attach(SS.sel.init.pars)
      #NSF, Survey and Others
    p1.sel=Phase_NSF_p1       
    p2.sel=Phase_NSF_p2 
    p3.sel=Phase_NSF_p3
    p4.sel=Phase_NSF_p4 
    p5.sel=Phase_NSF_p5 
    p6.sel=Phase_NSF_p6 
    if(NSF_Alternative_sel_type=='Logistic')
    {
      p1.sel=Phase_Logistic_p1       
      p2.sel=Phase_Logistic_p2 
    }
    if(NeiM%in%alternative.NSF.selectivity)
    {
      p1.sel_sens=Phase_Logistic_p1       
      p2.sel_sens=Phase_Logistic_p2
      p3.sel_sens=p4.sel_sens=p5.sel_sens=p6.sel_sens=NA
    }
    p1.sel_NSF=p1.sel_Other=p1.sel_Survey=p1.sel       
    p2.sel_NSF=p2.sel_Other=p2.sel_Survey=p2.sel 
    p3.sel_NSF=p3.sel_Other=p3.sel_Survey=p3.sel
    p4.sel_NSF=p4.sel_Other=p4.sel_Survey=p4.sel 
    p5.sel_NSF=p5.sel_Other=p5.sel_Survey=p5.sel 
    p6.sel_NSF=p6.sel_Other=p6.sel_Survey=p6.sel 
      #TDGLDF
    p1.sel_TDGDLF=Phase_TDGDLF_p1
    p2.sel_TDGDLF=Phase_TDGDLF_p2
    p3.sel_TDGDLF=Phase_TDGDLF_p3
    p4.sel_TDGDLF=Phase_TDGDLF_p4
    p5.sel_TDGDLF=Phase_TDGDLF_p5
    p6.sel_TDGDLF=Phase_TDGDLF_p6
    detach(SS.sel.init.pars)
    List.sp[[l]]$SS_selectivity_phase=data.frame(Fleet=c("Northern.shark","Other","Southern.shark_1","Southern.shark_2","Survey"),
                                                 P_1=c(p1.sel_NSF,p1.sel_Other,p1.sel_TDGDLF,p1.sel_TDGDLF,p1.sel_Survey),
                                                 P_2=c(p2.sel_NSF,p2.sel_Other,p2.sel_TDGDLF,p2.sel_TDGDLF,p2.sel_Survey),
                                                 P_3=c(p3.sel_NSF,p3.sel_Other,p3.sel_TDGDLF,p3.sel_TDGDLF,p3.sel_Survey),
                                                 P_4=c(p4.sel_NSF,p4.sel_Other,p4.sel_TDGDLF,p4.sel_TDGDLF,p4.sel_Survey),
                                                 P_5=c(p5.sel_NSF,p5.sel_Other,p5.sel_TDGDLF,p5.sel_TDGDLF,p5.sel_Survey),
                                                 P_6=c(p6.sel_NSF,p6.sel_Other,p6.sel_TDGDLF,p6.sel_TDGDLF,p6.sel_Survey))%>%
                                          replace(is.na(.), -2)
    if(NeiM%in%alternative.NSF.selectivity)
    {
      p1.sel_NSF=p1.sel_Other=p1.sel_Survey=p1.sel_sens       
      p2.sel_NSF=p2.sel_Other=p2.sel_Survey=p2.sel_sens
      p3.sel_NSF=p3.sel_Other=p2.sel_Survey=p3.sel_sens
      p4.sel_NSF=p4.sel_Other=p2.sel_Survey=p4.sel_sens
      p5.sel_NSF=p5.sel_Other=p2.sel_Survey=p5.sel_sens
      p6.sel_NSF=p6.sel_Other=p2.sel_Survey=p6.sel_sens
      List.sp[[l]]$SS_selectivity.sensitivity_phase=data.frame(Fleet=c("Northern.shark","Other","Southern.shark_1","Southern.shark_2","Survey"),
                                                   P_1=c(p1.sel_NSF,p1.sel_Other,p1.sel_TDGDLF,p1.sel_TDGDLF,p1.sel_Survey),
                                                   P_2=c(p2.sel_NSF,p2.sel_Other,p2.sel_TDGDLF,p2.sel_TDGDLF,p2.sel_Survey),
                                                   P_3=c(p3.sel_NSF,p3.sel_Other,p3.sel_TDGDLF,p3.sel_TDGDLF,p3.sel_Survey),
                                                   P_4=c(p4.sel_NSF,p4.sel_Other,p4.sel_TDGDLF,p4.sel_TDGDLF,p4.sel_Survey),
                                                   P_5=c(p5.sel_NSF,p5.sel_Other,p5.sel_TDGDLF,p5.sel_TDGDLF,p5.sel_Survey),
                                                   P_6=c(p6.sel_NSF,p6.sel_Other,p6.sel_TDGDLF,p6.sel_TDGDLF,p6.sel_Survey))%>%
        replace(is.na(.), -2)
    }
    
    #AgeSelex 
    List.sp[[l]]$age_selex_pattern=0  #0,Selectivity = 1.0 for ages 0+; 10, Selectivity = 1.0 for all ages beginning at age 1
    if(NeiM=="spinner shark") List.sp[[l]]$age_selex_pattern=10  #to allow convergence
    
    #4.1.4 lambdas (i.e. likelihood weighting)
  #note: only specify which lambdas are different to 1
  # Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
  # 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
  
  SS.lambdas=NULL
  if(NeiM%in%c("sandbar shark"))
  {
    SS.lambdas=data.frame(like_comp=1,
                          fleet=5,
                          value=2)
  }
  List.sp[[l]]$SS_lambdas=SS.lambdas
  
  #4.1.5 Recruitment pars
  #note: Rec Devs needed for model convergence for smooth hh, spinners; and for other species needed for improved cpue fit
  if(is.na(ramp.yrs$RecDev_Phase)) List.sp[[l]]$RecDev_Phase=-3 
  if(!is.na(ramp.yrs$RecDev_Phase)) List.sp[[l]]$RecDev_Phase=ramp.yrs$RecDev_Phase
  
  List.sp[[l]]$Steepness_Phase=-4    
  if(!is.na(ramp.yrs$Steepness_Phase)) List.sp[[l]]$Steepness_Phase=ramp.yrs$Steepness_Phase 
  
  List.sp[[l]]$SR_sigmaR=0.2    #Default (Spiny dogfish SS assessment)
  if(!is.na(ramp.yrs$SR_sigmaR)) List.sp[[l]]$SR_sigmaR=min(ramp.yrs$SR_sigmaR,Max.SR_sigmaR.shark) #from Report$sigma_R_info$alternative_sigma_R
  List.sp[[l]]$last_early_yr_nobias_adj_in_MPD=ramp.yrs$last_early_yr_nobias_adj_in_MPD
  List.sp[[l]]$first_yr_fullbias_adj_in_MPD=ramp.yrs$first_yr_fullbias_adj_in_MPD
  List.sp[[l]]$last_yr_fullbias_adj_in_MPD=ramp.yrs$last_yr_fullbias_adj_in_MPD
  List.sp[[l]]$first_recent_yr_nobias_adj_in_MPD=ramp.yrs$first_recent_yr_nobias_adj_in_MPD
  List.sp[[l]]$max_bias_adj_in_MPD=ramp.yrs$max_bias_adj_in_MPD
  
  #4.1.6 Catchabilities
  if(NeiM=="whiskery shark")
  {
    List.sp[[l]]$Nblock_Patterns=1
    if(Whiskery.q.periods==2)
    {
      List.sp[[l]]$blocks_per_pattern=2
      List.sp[[l]]$block_pattern_begin_end=list(c(1975,1980),c(1984,2005)) #targeting, not targeting periods
    }
    if(Whiskery.q.periods==3)
    {
      List.sp[[l]]$blocks_per_pattern=3
      List.sp[[l]]$block_pattern_begin_end=list(c(1975,1980),c(1981,1983),c(1984,2005)) #targeting, transitioning, not targeting periods
    }
    List.sp[[l]]$autogen=c(1,1,0,1,1) #1st biology, 2nd SR, 3rd Q, 4th reserved, 5th selex
    List.sp[[l]]$Q_param_Block=1
    List.sp[[l]]$Q_param_Blk_Fxn=2 #as per BigSkate assessment. 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
    
    List.sp[[l]]$Yr_q_change=1982   #last year before targeting practices changed (Simpfendorfer 2000)
    List.sp[[l]]$Yr_q_change_transition=1981:1983  #Disregard 80-83 as transitional years (Braccini et al 2021)
    List.sp[[l]]$Yr_q_daily=2006
  }else
  {
    List.sp[[l]]$Nblock_Patterns=0
    List.sp[[l]]$autogen=rep(1,5)
  }
  
  if(NeiM=="gummy shark")
  {
    if(Gummy.q.periods==2)
    {
      List.sp[[l]]$blocks_per_pattern=2
      List.sp[[l]]$block_pattern_begin_end=list(c(1975,2000),c(2001,2005)) #targeting, not targeting periods?
      List.sp[[l]]$Yr_second_q=2001:2005    #second monthly q due to increase in cpue and catch unaccounted by cpue stand.
    }
    
    
  }

  
  #-- 4.2 Bespoked Size-based integrated model
  if(Do.bespoked)
  {
    if(List.sp[[l]]$Species%in%Indicator.species)
    {
      # General input pars
      n.scen=2
      List.sp[[l]]=list.append(List.sp[[l]],
                               Data.yr=Last.yr.ktch,                          #last year of catch
                               Frst.yr.ktch=List.sp[[l]]$First.year,          #first year of catch
                               BaseCase=basecase,
                               Do.cols=do.cols, 
                               Max.FL.obs=LH$Max.FL.obs,                      #Maximum observed FL
                               Lo=LH$LF_o                                    #Size at birth
      ) 
      
      #Observed proportion of male sharks in TDGDLF catch
      List.sp[[l]]$Prop.males.ktch.zone=read.csv(paste(handl_OneDrive("Data/Population dynamics/Prop.males.in.catch/prop.males."),
                                                       List.sp[[l]]$SP,".csv",sep=""))
      List.sp[[l]]$Prop.males.ktch=read.csv(paste(handl_OneDrive("Data/Population dynamics/Prop.males.in.catch/prop.males"),
                                                  '.All.',List.sp[[l]]$SP,".csv",sep=""))%>%pull(x)
      
      
      h.M.constant.low=quantile(rlnorm(1000,log(store.species.steepness.S2[[l]]),store.species.steepness_M.mean[[l]]$sd),probs = 0.2)
      h.M.constant.up=quantile(rlnorm(1000,log(store.species.steepness.S2[[l]]),store.species.steepness_M.mean[[l]]$sd),probs = 0.8)
      
      
      # Species-specific input pars
      if(NeiM=="whiskery shark")
      {
        Drop_yr_cpue=c("1980-81","1981-82","1982-83","1983-84")  #Dropped cpue years (Braccini et al 2021)
        Drop_yr_cpue_sens=c("1975-76","1976-77","1977-78","1978-79","1979-80",
                            "1980-81","1981-82","1982-83")
        Drop_yr_cpue.tabl=paste(substr(Drop_yr_cpue,1,4)[1],
                                substr(Drop_yr_cpue,1,4)[length(Drop_yr_cpue)],sep="-")
        Drop_yr_cpue_sens.tabl=paste(substr(Drop_yr_cpue_sens,1,4)[1],
                                     substr(Drop_yr_cpue_sens,3,4)[length(Drop_yr_cpue_sens)],sep="-")
        
        List.sp[[l]]=list.append(List.sp[[l]],
                                 Prior.mean.Fo=0.01,
                                 Prior.SD.Log.Fo=0.5,
                                 
                                 #Steepness
                                 h.M.constant=h.M_mean2,       #0.351 Braccini et al 2015
                                 h.M.constant.low=h.M.constant.low,    #80% percentile
                                 h.M.constant.up=h.M.constant.up,
                                 
                                 #Natural mortality
                                 M_val=Mmean,          #0.27 using Hoenig set to Max age=15  (Simpfendorfer et al 2000)
                                 M_val.low=min(apply(store.species.M_M.min[[l]],2,mean,na.rm=T)),   
                                 M_val.high=max(apply(store.species.M_M.min[[l]],2,mean,na.rm=T)),  
                                 
                                 #Initial F in 1974
                                 Fo=0.03,             #yields a B1975 of 90% virgin           
                                 Fo_Simp=0.003,       #Simpfendorfer et al 2000 estimated at 0.003
                                 Fo_M=0.05,            #yields a B1975 of 85% virgin 
                                 
                                 #      Po_spm=0.9,  #Po for surplus production, consistent with the Fo value used in Size based model
                                 
                                 #Data
                                 AREAS=c("West","Zone1","Zone2"),  #Define spatial areas; 1 is West, 2 is Zn1, 3 is Zn2.
                                 Yr_q_change=1982,   #last year before targeting practices changed (Simpfendorfer 2000)
                                 Yr_q_change_transition=1981:1983,  #Disregard 80-83 as transitional years (Braccini et al 2021)
                                 Yr_q_daily=2006,
                                 Do_var=0,    #How to calculate cpue variance in age-structured
                                 Var1=0.0296133485565278,   #fixed variances used by Simpfendorfer et al 2000
                                 Var2=0.023525088856464,
                                 
                                 #Initial value of estimable parameters
                                 #note: Dummy used for switching on/off phase estimation in simulation testing
                                 Dummy=1,   
                                 #----Size-base
                                 RZERO_in_1000_individuals_SB=fn.ji(1096),
                                 Q1_SB=fn.ji(0.00041),
                                 Q2_SB=fn.ji(0.0001),
                                 Q_daily_SB=fn.ji(0.00024),
                                 Fo_SB=NA,  #no jit because it's fixed
                                 tau_SB=fn.ji(0.2),
                                 K.F=fn.ji(List.sp[[l]]$Growth.F$k),  
                                 Linf.F=fn.ji(List.sp[[l]]$Growth.F$FL_inf),
                                 K.M=fn.ji(List.sp[[l]]$Growth.M$k),
                                 Linf.M=fn.ji(List.sp[[l]]$Growth.M$FL_inf),
                                 SD.growth_SB=fn.ji(7),
                                 Prop.west=fn.ji(0.19),  #proportion of recruitment by zone calculated based on average prop catch by zone
                                 Prop.zn1=fn.ji(0.38),
                                 #---movement 
                                 p11=0.99,
                                 p22=0.99,
                                 p21=0.001,
                                 p33=0.99,
                                 
                                 #Estimation phases
                                 Par.phases=list('Base case'=c(dummy=Fz.off,lnR_zero=1,lnR_prop_west=1,lnR_prop_zn1=1,
                                                               lnq=2,lnq2=2,log_Qdaily=2,ln_Init_F=Fz.off,log_tau=5,
                                                               k=3,lnLinf=3,k_M=3,lnLinf_M=3,sd_growth=4,
                                                               log_p11=1,log_p22=1,log_p21=1,log_p33=1))
        )
        
      }
      
      if(NeiM=="gummy shark")
      {
        List.sp[[l]]=list.append(List.sp[[l]],
                                 Prior.mean.Fo=0.01,
                                 Prior.SD.Log.Fo=0.5,
                                 
                                 #Steepness
                                 h.M.constant=h.M_mean2,  #0.481 Braccini et al 2015 
                                 h.M.constant.low=h.M.constant.low,    #80% percentile
                                 h.M.constant.up=h.M.constant.up,
                                 
                                 #Natural mortality
                                 M_val=Mmean,          
                                 M_val.low=min(apply(store.species.M_M.min[[l]],2,mean,na.rm=T)),   
                                 M_val.high=max(0.283,max(apply(store.species.M_M.min[[l]],2,mean,na.rm=T))),  #0.283 Walker empirical
                                 
                                 #Initial F in 1974
                                 Fo=0.05,               #This leaves B1975 at 95% unfished 
                                 Fo_Simp=0.003,              
                                 Fo_M=0.1,                
                                 Fo_AS=0.004,            #This leaves B1975 at 95% unfished 
                                 
                                 #       Po_spm=0.95,  #Po for surplus production, consistent with the Fo value used in Size based model
                                 
                                 #Data
                                 AREAS=c("West","Zone1","Zone2"),  #Define spatial areas; 1 is West, 2 is Zn1, 3 is Zn2.
                                 Yr_q_change=0,   #last year before targeting practices changed (Simpfendorfer 2000)
                                 Yr_q_daily=2006,
                                 Do_var=0,     #How to calculate cpue variance in Simpfendorfer's age-structured
                                 
                                 #Initial value of estimable parameters
                                 #note: Dummy is used for switching on/off phase estimation in simulation testing
                                 Dummy=1, 
                                 
                                 #---Size-base
                                 RZERO_in_1000_individuals_SB=fn.ji(1000),
                                 Q1_SB=fn.ji(1e-4),
                                 Q2_SB=fn.ji(1e-4),
                                 Q_daily_SB=fn.ji(1e-4),
                                 Fo_SB=NA,   #no jit because it's fixed
                                 tau_SB=fn.ji(0.3),
                                 K.F=fn.ji(List.sp[[l]]$Growth.F$k),  
                                 Linf.F=fn.ji(List.sp[[l]]$Growth.F$FL_inf),
                                 K.M=fn.ji(List.sp[[l]]$Growth.M$k),
                                 Linf.M=fn.ji(List.sp[[l]]$Growth.M$FL_inf),
                                 SD.growth_SB=fn.ji(10),
                                 Prop.west=fn.ji(0.02),  #proportion of recruitment by zone calculated based on average prop catch by zone
                                 Prop.zn1=fn.ji(0.08),
                                 #Movement 
                                 p11=0.999,
                                 p22=0.999,
                                 p21=0.00024,
                                 p33=0.999,
                                 
                                 #Estimation phases
                                 Par.phases=list('Base case'=c(dummy=Fz.off,lnR_zero=3,lnR_prop_west=-3,lnR_prop_zn1=-3,
                                                               lnq=4,lnq2=Fz.off,log_Qdaily=4,ln_Init_F=Fz.off,log_tau=5,
                                                               k=1,lnLinf=1,k_M=1,lnLinf_M=1,sd_growth=2,
                                                               log_p11=1,log_p22=1,log_p21=1,log_p33=1))
                                 
        )
      }
      
      if(NeiM=="dusky shark")    #update all input pars, currently using gummies!  
      {
        List.sp[[l]]=list.append(List.sp[[l]],
                                 Prior.mean.Fo=0.01,
                                 Prior.SD.Log.Fo=0.5,
                                 
                                 #Steepness
                                 h.M.constant=h.M_mean2,  #0.234 Braccini et al 2015 
                                 h.M.constant.low=h.M.constant.low,    #80% percentile
                                 h.M.constant.up=h.M.constant.up,
                                 
                                 #Natural mortality
                                 M_val=Mmean,          
                                 M_val.low=min(apply(store.species.M_M.min[[l]],2,mean,na.rm=T)),   
                                 M_val.high=max(apply(store.species.M_M.min[[l]],2,mean,na.rm=T)),     
                                 
                                 #Initial F in 1974
                                 Fo=0.05,               #This leaves B1975 at 95% unfished 
                                 Fo_Simp=0.003,              
                                 Fo_M=0.1,                
                                 Fo_AS=0.004,            #This leaves B1975 at 95% unfished 
                                 
                                 #    Po_spm=0.95,  #Po for surplus production, consistent with the Fo value used in Size based model
                                 
                                 #Data
                                 Ktch.source="WA.only",  #select whether to use all catch series or only WA
                                 #Ktch.source="ALL",
                                 AREAS=c("West","Zone1","Zone2"),  #Define spatial areas; 1 is West, 2 is Zn1, 3 is Zn2.
                                 Yr_q_change=0,   #last year before targeting practices changed (Simpfendorfer 2000)
                                 Yr_q_daily=2006,
                                 Do_var=0,     #How to calculate cpue variance in Simpfendorfer's age-structured
                                 
                                 #Initial value of estimable parameters
                                 #note: Dummy is used for switching on/off phase estimation in simulation testing
                                 Dummy=1,   
                                 
                                 #---Size-base
                                 RZERO_in_1000_individuals_SB=fn.ji(1000),
                                 Q1_SB=fn.ji(1e-4),
                                 Q2_SB=fn.ji(1e-4),
                                 Q_daily_SB=fn.ji(1e-4),
                                 Fo_SB=NA,   #no jit because it's fixed
                                 tau_SB=fn.ji(0.3),
                                 K.F=fn.ji(List.sp[[l]]$Growth.F$k),  
                                 Linf.F=fn.ji(List.sp[[l]]$Growth.F$FL_inf),
                                 K.M=fn.ji(List.sp[[l]]$Growth.M$k),
                                 Linf.M=fn.ji(List.sp[[l]]$Growth.M$FL_inf),
                                 SD.growth_SB=fn.ji(10),
                                 Prop.west=fn.ji(0.02),  #proportion of recruitment by zone calculated based on average prop catch by zone
                                 Prop.zn1=fn.ji(0.08),
                                 
                                 #Movement 
                                 p11=0.999,
                                 p22=0.999,
                                 p21=0.00024,
                                 p33=0.999,
                                 
                                 #Estimation phases
                                 Par.phases=list('Base case'=c(dummy=Fz.off,lnR_zero=3,lnR_prop_west=-3,lnR_prop_zn1=-3,
                                                               lnq=4,lnq2=Fz.off,log_Qdaily=4,ln_Init_F=Fz.off,log_tau=5,
                                                               k=1,lnLinf=1,k_M=1,lnLinf_M=1,sd_growth=2,
                                                               log_p11=1,log_p22=1,log_p21=1,log_p33=1))
        )
        
      }
      
      if(NeiM=="sandbar shark")   #update all input pars, currently using gummies!
      {
        List.sp[[l]]=list.append(List.sp[[l]],
                                 Prior.mean.Fo=0.01,
                                 Prior.SD.Log.Fo=0.5,
                                 
                                 #Steepness
                                 h.M.constant=h.M_mean2,   #0.216 Braccini et al 2015
                                 h.M.constant.low=h.M.constant.low,    #80% percentile
                                 h.M.constant.up=h.M.constant.up,
                                 
                                 #Natural mortality
                                 M_val=Mmean,          
                                 M_val.low=min(apply(store.species.M_M.min[[l]],2,mean,na.rm=T)),   
                                 M_val.high=max(apply(store.species.M_M.min[[l]],2,mean,na.rm=T)),    
                                 
                                 #Initial F in 1974
                                 Fo=0.05,               #This leaves B1975 at 95% unfished 
                                 Fo_Simp=0.003,              
                                 Fo_M=0.1,                
                                 Fo_AS=0.004,            #This leaves B1975 at 95% unfished 
                                 
                                 # Po_spm=0.95,  #Po for surplus production, consistent with the Fo value used in Size based model
                                 
                                 #Data
                                 Ktch.source="WA.only",  #select whether to use all catch series or only WA
                                 #Ktch.source="ALL",
                                 AREAS=c("West","Zone1","Zone2"),  #Define spatial areas; 1 is West, 2 is Zn1, 3 is Zn2.
                                 Yr_q_change=0,   #last year before targeting practices changed (Simpfendorfer 2000)
                                 Yr_q_daily=2006,
                                 Do_var=0,     #How to calculate cpue variance in Simpfendorfer's age-structured
                                 
                                 #Initial value of estimable parameters
                                 #note: Dummy is used for switching on/off phase estimation in simulation testing
                                 Dummy=1,   
                                 
                                 #---Size-base
                                 RZERO_in_1000_individuals_SB=fn.ji(1000),
                                 Q1_SB=fn.ji(1e-4),
                                 Q2_SB=fn.ji(1e-4),
                                 Q_daily_SB=fn.ji(1e-4),
                                 Fo_SB=NA,   #no jit because it's fixed
                                 tau_SB=fn.ji(0.3),
                                 K.F=fn.ji(List.sp[[l]]$Growth.F$k),  
                                 Linf.F=fn.ji(List.sp[[l]]$Growth.F$FL_inf),
                                 K.M=fn.ji(List.sp[[l]]$Growth.M$k),
                                 Linf.M=fn.ji(List.sp[[l]]$Growth.M$FL_inf),
                                 SD.growth_SB=fn.ji(10),
                                 Prop.west=fn.ji(0.02),  #proportion of recruitment by zone calculated based on average prop catch by zone
                                 Prop.zn1=fn.ji(0.08),
                                 
                                 #Movement 
                                 p11=0.999,
                                 p22=0.999,
                                 p21=0.00024,
                                 p33=0.999,
                                 
                                 #Estimation phases
                                 Par.phases=list(
                                   'Base case'=c(dummy=Fz.off,lnR_zero=3,lnR_prop_west=-3,lnR_prop_zn1=-3,
                                                 lnq=4,lnq2=Fz.off,log_Qdaily=4,ln_Init_F=Fz.off,log_tau=5,
                                                 k=1,lnLinf=1,k_M=1,lnLinf_M=1,sd_growth=2,
                                                 log_p11=1,log_p22=1,log_p21=1,log_p33=1))
        )
      }
      
      # Model scenarios
      List.sp[[l]]$N.Scens=n.scen       
      Zens=paste("S",1:(n.scen-1),sep="")
      Models=c("Base case",paste("S",1:(n.scen-1),sep=""))
      Q.scen=c(rep("three",2))
      
      Tabla.scen=data.frame(
        Model=Models,
        Size_comp.=c(rep('Yes',n.scen)),
        CPUE=rep("Stand.",n.scen),
        Age.Growth=rep('Yes',n.scen),
        Ktch.sx.r=rep('Observed',n.scen),                      
        Tagging=rep('No',n.scen),                     
        Fec.=rep('N/A',n.scen),
        Maturity=rep("at length",n.scen),
        M=rep("at length",n.scen),
        SteepnesS=c(h.M_mean2,h.M.constant.low),
        Q=Q.scen,   
        Spatial_structure=rep('1 zone',n.scen),
        Movement=rep("No",n.scen),
        Fo=rep(List.sp[[l]]$Fo,n.scen),
        Model_type=rep("Length-based",n.scen)) 
      
      if(NeiM=="whiskery shark") Tabla.scen$CPUE_years_dropped=rep(Drop_yr_cpue.tabl,n.scen)
      List.sp[[l]]$Tabla.scen=Tabla.scen  
      
      n.areas=length(List.sp[[l]]$AREAS)
      List.sp[[l]]$n.areas=n.areas
      List.sp[[l]]$Areas.zones=data.frame(area=1:n.areas,zone=List.sp[[l]]$AREAS)
      List.sp[[l]]$hndl=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(List.sp[[l]]$Name),"/",sep='')
      List.sp[[l]]$Fo_SB=List.sp[[l]]$Fo  #fixed
      
      rm(n.areas,Tabla.scen)
      
      # Create folders if new run
      if(First.run=='YES')
      {
        kriat.path=paste(List.sp[[l]]$hndl,AssessYr,'/Integrated sized-structured',sep="")
        if(!dir.exists(kriat.path))dir.create(kriat.path)
        kriate.this=as.character(List.sp[[l]]$Tabla.scen$Model)
        for (i in 1:length(kriate.this)) 
        { 
          NEW=paste(kriat.path,kriate.this[i], sep="/")
          if(!dir.exists(NEW))dir.create(NEW) 
        }
      }
      
      #drop conv tag if not using
      if(add.conv.tag=="NO") List.sp[[l]]$Tabla.scen=List.sp[[l]]$Tabla.scen[,-match("Tagging",names(List.sp[[l]]$Tabla.scen))]
      
      # Export scenarios table
      inpt.pz=paste(List.sp[[l]]$hndl,AssessYr,'/Integrated sized-structured/',sep="")
      setwd(inpt.pz)
      if(NeiM=="whiskery shark")
      {
        THiS=c("Model_name","Model_type","Spatial_structure"
               ,"Movement","Size_comp.","CPUE","CPUE_years_dropped","Age.Growth"         
               ,"Ktch.sx.r","Tagging","Steepness","Fo","Maturity","Q")
        HDR.span=c(2,1,1,6,4,1)
        HDR.2nd=c("Name","Type","structure",'',"Size","CPUE","CPUE years","Age &",
                  "Prop. male","Tagging","M","h","Fo","Maturity","")
        HDR.3rd=c("","","","","composition","","not used in likelihood","growth","in catch",
                  "","","","","","")
      }else
      {
        THiS=c("Model_name","Model_type","Spatial_structure"
               ,"Movement","Size_comp.","CPUE","Age.Growth"         
               ,"Ktch.sx.r","Tagging","Steepness","Fo","Maturity","Q")
        HDR.span=c(2,1,1,5,4,1)
        HDR.2nd=c("Name","Type","structure",'',"Size","CPUE","Age &",
                  "Prop. male","Tagging","M","h","Fo","Maturity","")
        HDR.3rd=c("","","","","composition","","growth","in catch",
                  "","","","","","")
      }
      if(First.run=='YES')
      {
        Tabla.scen.show=List.sp[[l]]$Tabla.scen
        Tabla.scen.show=Tabla.scen.show[,-match(c("Fec.","M"),colnames(Tabla.scen.show))]
        Tabla.scen.show[is.na(Tabla.scen.show)]="N/A"
        names(Tabla.scen.show)[match(c("Model","SteepnesS"),names(Tabla.scen.show))]=c("Model_name","Steepness")
        Tabla.scen.show=Tabla.scen.show[,match(THiS,names(Tabla.scen.show))]
        Tabla.scen.show=Tabla.scen.show[match(c(Zens,"Base case" ),Tabla.scen.show$Model_name),]
        Tabla.scen.show$Q=as.character(Tabla.scen.show$Q)
        Tabla.scen.show$Q=with(Tabla.scen.show,ifelse(Q=="three",3,ifelse(Q=="two",2,ifelse(Q=="one",1,Q))))
        Tabla.scen.show$Q=with(Tabla.scen.show,ifelse(!Q=="N/A"& Q>1,paste(Q,"periods"),
                                                      ifelse(!Q=="N/A"& Q==1,paste(Q,"period"),Q)))
        Tabla.scen.show=Tabla.scen.show%>%
          mutate(Steepness=ifelse(!Steepness=='N/A',as.character(round(as.numeric(Steepness),3)),Steepness))
        fn.word.table(TBL=Tabla.scen.show,Doc.nm="Integrated Model_scenarios")
        
      }
      
      # Create pin file
      if(First.run=='YES')
      {
        with(List.sp[[l]],
             {
               #Population pin values
               Pin.pars=vector('list',nrow(Tabla.scen))
               names(Pin.pars)=Tabla.scen$Model
               
               #Length-based
               Pin.pars$'Base case'=c(Dummy=Dummy,lnR_zero=log(RZERO_in_1000_individuals_SB),
                                      R_prop_west=Prop.west,R_prop_zn1=Prop.zn1,
                                      lnq=log(Q1_SB),lnq2=log(Q2_SB),lnq_daily=log(Q_daily_SB),
                                      ln_Init_F=log(Fo_SB),ln_tau=log(tau_SB),
                                      k=K.F,lnLinf=log(Linf.F),k_M=K.M,lnLinf_M=log(Linf.M),
                                      sd_growth=SD.growth_SB)
               
               #Add tagging pars       
               if(conv.tag.all=="YES") 
               {
                 if(Move.mode=="Individual-based") 
                 {
                   Mov.pars=c(p11=p11,p22=p22,p21=p21,p33=p33)
                 }
                 
                 if(Move.mode=="Population-based")
                 {
                   Mov.pars=c(mov.WC_WC=mov.WC_WC,mov.WC_ZN1=mov.WC_ZN1,
                              mov.ZN1_WC=mov.ZN1_WC,mov.ZN1_ZN1=mov.ZN1_ZN1,
                              mov.ZN2_ZN1=mov.ZN2_ZN1,mov.ZN2_ZN2=mov.ZN2_ZN2,
                              log.Q.tag.WC=log.Q.tag.WC,log.Q.tag.ZN1=log.Q.tag.ZN1,
                              log.Q.tag.ZN2=log.Q.tag.ZN2,log.tau=log.tau)
                 }
                 Pin.pars$'Base case'=c(Pin.pars$'Base case',Mov.pars)
               }
               
               #Set the pins from other size-based scenarios
               IDs=which(Tabla.scen$Model_type=="Length-based")
               IDs=IDs[-match("Base case",Tabla.scen$Model[IDs])]
               for(i in IDs) Pin.pars[[i]]=Pin.pars$'Base case'
               
               for(i in 1:N.Scens)     
               {
                 if(Tabla.scen$Model_type[i]=="Length-based")
                 {
                   #higher M or Fo needs larger initial population
                   if(any(Tabla.scen$M.value[i]==M_val.high,as.numeric(Tabla.scen$Fo[i])==Fo_M))    
                   {
                     ss=match("lnR_zero",names(Pin.pars[[i]]))
                     Pin.pars[[i]][ss]=Pin.pars[[i]][ss]*1.15
                   }
                   
                   #no movement outside zone
                   if(Tabla.scen$Movement[i]%in%c("No","N/A"))
                   {
                     Pin.pars[[i]][match(c("p11","p22","p21","p33"),names(Pin.pars[[i]]))]=c(1,1,0,1)      
                   }
                   
                   #Change Fo if required
                   if(!as.numeric(as.character(Tabla.scen$Fo[i]))==Fo) 
                   {
                     Pin.pars[[i]][match("ln_Init_F",names(Pin.pars[[i]]))]=log(as.numeric(as.character(Tabla.scen$Fo[i])))
                     ss=match("lnR_zero",names(Pin.pars[[i]]))
                     #higher Fo needs larger initial population
                     if(!as.numeric(as.character(Tabla.scen$Fo[i]))>Fo)Pin.pars[[i]][ss]=Pin.pars[[i]][ss]*1.15
                   }
                   
                 }
               }
               
               #Export as pin file
               setPath=function(Scen)setwd(paste(List.sp[[l]]$hndl,AssessYr,"/Integrated sized-structured/",Scen,sep=""))
               for(i in 1:nrow(Tabla.scen))
               {
                 These.pars=Pin.pars[[i]]
                 par.nms=names(These.pars)
                 setPath(Tabla.scen[i,]$Model)
                 FILE=paste(List.sp[[l]]$Name.inputs,".pin",sep="")
                 write("# Input parameters",file = FILE)
                 for(k in 1:length(These.pars))
                 {
                   Hdr=paste("#",par.nms[k])
                   write(Hdr,file = FILE,append=T)
                   write(These.pars[k],file = FILE,sep = "\t",append=T)
                 }
               }
               List.sp[[l]]$Pin.pars=Pin.pars
               
               
               #Set estimation phases for scenarios not already assigned
               IDs=names(Pin.pars)[-match(c("Base case"),names(Pin.pars))]
               dummy=vector('list',length(IDs))
               names(dummy)=IDs
               for(i in 1:length(dummy)) dummy[[i]]=Par.phases$'Base case'
               Par.phases=c(Par.phases,dummy)
               
               #Turn off irrelevant pars according to scenarios 
               for(i in 1:N.Scens)     
               {
                 if(Tabla.scen$Model_type[i]=="Length-based")
                 {
                   #single zone
                   if(Tabla.scen$Spatial_structure[i]=="1 zone")
                   {
                     Par.phases[[i]][fn.mtch(Zns.par.phz,Par.phases[[i]])]=rep(Fz.off,length(Zns.par.phz))
                     Par.phases[[i]][fn.mtch(Q_phz,Par.phases[[i]])]=c(1,1,1)
                   }
                   #effective cpue
                   if(Tabla.scen$CPUE[i]=="Effective")
                   {
                     Par.phases[[i]][fn.mtch(Q_phz,Par.phases[[i]])]=c(1,1,Fz.off)
                     Par.phases[[i]][fn.mtch("log_tau",Par.phases[[i]])]=Fz.off
                   }
                   #no Qs if no cpue used
                   if(Tabla.scen$CPUE[i]%in%c("N/A","No"))
                   {
                     Par.phases[[i]][fn.mtch(c(Q_phz),Par.phases[[i]])]=rep(Fz.off,length(c(Q_phz)))
                     Par.phases[[i]][fn.mtch("lnR_zero",Par.phases[[i]])]=1
                   }
                   #no movement
                   if(Tabla.scen$Movement[i]%in%c("No","N/A"))
                   {
                     Par.phases[[i]][fn.mtch(MOv.par.phz,Par.phases[[i]])]=rep(Fz.off,length(MOv.par.phz))  
                   }
                   
                   #one less q
                   if(Tabla.scen$Q[i]=="two" & Tabla.scen$CPUE[i]=="Stand.")
                   {
                     Par.phases[[i]][fn.mtch("lnq2",Par.phases[[i]])]=-Par.phases[[i]][fn.mtch("lnq2",Par.phases[[i]])]
                   }
                 }
               }
               List.sp[[l]]$Par.phases=Par.phases
             })
      }
      
      #Define general arguments
      List.sp[[l]]=list.append(List.sp[[l]],
                               
                               #Select acoustic tagging model
                               Acoust.format=Move.mode,
                               #Acoust.format="SS3"
                               
                               #Select type of size composition Likelihood
                               Size_like=size.likelihood,
                               Dirichlet.small.value=1e-4,   #small constant for intermediate 0 observations 
                               
                               #Combined size composition?  
                               Size.sex.comb=size.sex.comb,
                               
                               #Size composition as proportions?
                               Size.comp.prop=size.comp.prop,
                               
                               #Effective sample size size composition
                               Effective.n=Min.annual.obs,
                               
                               #Maximum possible F
                               MaxF=maxF,
                               
                               #Prior for initial fishing mortality
                               add_Finit_prior=Add_Finit_prior,
                               
                               #running arguments
                               Arguments=arguments,
                               
                               #Do MSY calculation using Base Case model   
                               Do.MSY=do.MSY,
                               MSY.yrs=msy.yrs,  
                               MSY.sd.rec=msy.sd.rec,
                               MSY.sims=msy.sims,
                               F.vec=f.vec,
                               
                               #MCMC 
                               DO.MCMC=mcmc.do,
                               burning=1:(5*length(seq(1,mcmc.n,by=mcmc.thin))/100),   #5%  burning
                               
                               #What scenarios to run
                               Run.all.Scenarios=run.all.scenarios,
                               
                               #Number of years for future projections (use generation time)
                               Yrs.future=round(store.species.G_M.min[[l]]$mean),  
                               
                               #Future projection scenarios
                               Run.future.proj=run.future,
                               Future.ktch.scen=list(mean=1,percent_50=.5,percent_150=1.5)
      )
      
      #Weight of likelihood components
      List.sp[[l]]$Rho=1  #weight of size comp 
      List.sp[[l]]$Rho2=1 #weight of age-length
      List.sp[[l]]$Rho3=1  #weight of cpue
      
      #Do future projections with 0 catch?
      if(Do.zero.Ktch=="YES")
      {
        List.sp[[l]]$zero.Ktch.yrs=100
        List.sp[[l]]$zero.Ktch.this="Base case"
      }
      
      #Do simulation testing?
      if(Do.sim.Test=="YES")
      {
        List.sp[[l]]$N.sim.test=50
        List.sp[[l]]$CPUE.jitr=1000 
        List.sp[[l]]$Sim.Test.this="Base case"
      }
    } # end Indicator.species
  } #end Bespoked

  rm(NeiM)
}

#remove duplicated objects from multiple runs of List.sp loops
for(l in 1:N.sp) List.sp[[l]]=List.sp[[l]][!duplicated(names(List.sp[[l]]))]


#Get age at first maturity
fn.age.mat.per=function(d,mat.per=0.1)
{
  eig=0:d$Max.age.F[2]
  TLo=d$Lzero*d$a_FL.to.TL+d$b_FL.to.TL
  TL=TLo+(d$Growth.F$FL_inf*d$a_FL.to.TL+d$b_FL.to.TL-TLo)*(1-exp(-d$Growth.F$k*eig))
  Maturity=1/(1+exp(-log(19)*((TL-d$TL.50.mat)/(d$TL.95.mat-d$TL.50.mat))))
  mature.percent=TL[which(abs(Maturity - mat.per) == min(abs(Maturity - mat.per)))]
  age.percent.mature=log((1-((mature.percent-TLo)/(d$Growth.F$FL_inf*d$a_FL.to.TL+d$b_FL.to.TL-TLo))))/(-d$Growth.F$k)
  return(age.percent.mature)
}
for(l in 1:N.sp) List.sp[[l]]$First_Mature_Age=fn.age.mat.per(d=List.sp[[l]])

#Get maturity slope and inflection for SS3
fn.predlog=function(pars) 1/(1+exp(pars[1]*(dat$x-pars[2])))
fn.fit.log.infl=function(pars)
{
  pred= fn.predlog(pars)
  return(sum((dat$y-pred)^2))
}
for(l in 1:N.sp)
{
  dat=data.frame(x=with(List.sp[[l]],(Lzero*a_FL.to.TL+b_FL.to.TL):List.sp[[l]]$TLmax))%>%
    mutate(y=with(List.sp[[l]],1/(1+exp(-log(19)*((x-TL.50.mat)/(TL.95.mat-TL.50.mat))))),
           x=base::jitter(x,50),
           y=base::jitter(y,50))
  fit=nlminb(start=c(-1,List.sp[[l]]$TL.50.mat), objective=fn.fit.log.infl)   
  dd=fit$par
  names(dd)=c('slope','inflection')
  List.sp[[l]]$TL.mat.inf.slope=dd
  #dat=dat%>%arrange(x)
  #plot(dat$x,dat$y,main=names(List.sp)[l])
  #lines(dat$x,fn.predlog(fit$par),col=2,lwd=2)
  rm(dat,dd)
}

#Get selectivity p1 and p2 for SS3
fn.predlog=function(pars) 1/(1+exp(-log(19)*((dat$x-pars[1])/pars[2])))
for(l in 1:N.sp)
{
  dat=data.frame(x=with(List.sp[[l]],(Lzero*a_FL.to.TL+b_FL.to.TL):List.sp[[l]]$TLmax))%>%
    mutate(y=with(List.sp[[l]],1/(1+exp(TL.mat.inf.slope[1]*(x-TL.mat.inf.slope[2])))))
  fit=nlminb(start=c(List.sp[[l]]$TL.50.mat,List.sp[[l]]$TL.50.mat*1.1), objective=fn.fit.log.infl)   
  dd=fit$par
  names(dd)=c('p1','p2')
  List.sp[[l]]$Logistic.selectivity=dd
  rm(dat,dd)
}

