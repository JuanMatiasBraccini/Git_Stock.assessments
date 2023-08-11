#remotes::install_github("JABBAmodel/ss3diags")
library(ss3diags)  
# Create SS input files ------------------------------------------------------
fn.get.in.betwee=function(x,PATRN="_") str_before_nth(str_after_nth(x, PATRN, 1), PATRN, 2)
fn.set.up.SS=function(Templates,new.path,Scenario,Catch,life.history,depletion.yr,
                      fleets=NULL,fleetinfo=NULL,abundance=NULL,size.comp=NULL,meanbodywt=NULL,
                      F.tagging=NULL,cond.age.len=NULL,MeanSize.at.Age.obs=NULL,Lamdas=NULL,
                      RecDev_Phase=-3,SR_sigmaR=0.2,Var.adjust.factor=NULL)
{
  # 1.Copy templates
  copy_SS_inputs(dir.old = Templates, dir.new = new.path,overwrite = TRUE)
  
  
  # 2.Read in templates 
  start <- r4ss::SS_readstarter(file = file.path(new.path, "starter.ss"), verbose = FALSE)
  dat <- r4ss::SS_readdat(file = file.path(new.path, start$datfile), verbose = FALSE)
  ctl <- r4ss::SS_readctl(file = file.path(new.path, start$ctlfile), verbose = FALSE, use_datlist = TRUE, datlist = dat)
  fore <- r4ss::SS_readforecast(file = file.path(new.path, "forecast.ss"),  verbose = FALSE)
  
  
  # 3.Update template with species-specific information
  
  #3.1. dat file
  dat$Comments[3]=paste('#C file write time:',Sys.time())
  dat$spawn_month=1.0  #integer is month (1-12), decimal is fraction of days in month (e.g. 15 March is 3.5)
  
  #general age info
  dat$Nages=max(life.history$Max.age.F)
  ageError=as.data.frame(matrix(nrow=2,ncol=dat$Nages+1))
  ageError[1,]=-1.00
  ageError[2,]=0.001
  names(ageError)=seq(0,dat$Nages)
  dat$ageerror=ageError
  
  #population size classes
  dat$binwidth=TL.bins.cm
  dat$minimum_size=10*floor(0.95*with(life.history,Lzero*a_FL.to.TL+b_FL.to.TL)/10)
  dat$maximum_size=10*ceiling(1.01*with(life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/10)
  #dat$maximum_size=30*ceiling(with(life.history,max(c(TLmax,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)))/30)
  
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
    ddumy=dat$CPUEinfo%>%
      rownames_to_column('fleetname')%>%
      mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                              ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                                     fleetname)))%>%
      filter(fleetname%in%dis.flits)%>%
      mutate(Fleet=row_number())
    if(Abundance.error.dist=='Normal') ddumy$Errtype=-1  
    if(Abundance.error.dist=='Lognormal') ddumy$Errtype=0
    rownames(ddumy)=ddumy$fleetname
    dat$CPUEinfo=ddumy%>%dplyr::select(-fleetname)
    if(!is.null(abundance)) dat$CPUE=abundance%>%mutate(Mean=ifelse(Mean<1e-6,1e-6,Mean))
    if(is.null(abundance))  dat$CPUE=NULL
    
    #Fishing mortality from tagging  
    if(!is.null(F.tagging))
    {
      #notes: to use an F series need to add an additional fleet, cpue series, q and mirror selectivity
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
      ddumy=dat$len_info%>%
        rownames_to_column('fleetname')%>%
        mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                                ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                                       fleetname)))%>%
        filter(fleetname%in%dis.flits)
      rownames(ddumy)=ddumy$fleetname
      dat$len_info=ddumy%>%dplyr::select(-fleetname)
      lbin_vector=sort(as.numeric(gsub('f', '', names(size.comp)[grep("f",names(size.comp))])))
      dat$lbin_vector=lbin_vector
      dat$N_lbins=length(dat$lbin_vector)
      dat$lencomp=size.comp%>%arrange(Sex,Fleet,year)
    }
    if(!is.null(F.tagging))
    {
      addfbit=dat$len_info[length(unique(F.tagging$fleet)),]
      rownames(addfbit)=rownames(Finfo)
      dat$len_info=rbind(dat$len_info,addfbit)
    }
    
    #conditional age at length  
    ddumy=dat$age_info%>%
      rownames_to_column('fleetname')%>%
      mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                              ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                                     fleetname)))%>%
      filter(fleetname%in%dis.flits)
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
      dat$agebin_vector=seq(0,dat$Nages-1)
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
  if(!life.history$Nblock_Patterns==0 & !is.null(abundance))
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
  #females
  ctl$MG_parms["NatM_p_1_Fem_GP_1", c("INIT","PRIOR")]=rep(Scenario$Mmean,2)
  ctl$MG_parms["NatM_p_1_Fem_GP_1", c("LO","HI")]=c(ctl$MG_parms["NatM_p_1_Fem_GP_1", "INIT"]*.1,ctl$MG_parms["NatM_p_1_Fem_GP_1", "INIT"]*4)
  ctl$MG_parms["L_at_Amin_Fem_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Lzero*a_FL.to.TL+b_FL.to.TL),2)  #size for specified _Age(post-settlement)_for_L1
  ctl$MG_parms["L_at_Amin_Fem_GP_1", c("LO","HI")]=round(c(ctl$MG_parms["L_at_Amin_Fem_GP_1", "INIT"]*.25,ctl$MG_parms["L_at_Amin_Fem_GP_1", "INIT"]*2))
  ctl$MG_parms["L_at_Amax_Fem_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL),2) #life.history$TLmax  #set at Linf as _Growth_Age_for_L2 was set at 999
  ctl$MG_parms["L_at_Amax_Fem_GP_1", c("LO","HI")]=round(c(ctl$MG_parms["L_at_Amax_Fem_GP_1", "INIT"]*.25,ctl$MG_parms["L_at_Amax_Fem_GP_1", "INIT"]*2))
  ctl$MG_parms["Wtlen_1_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$AwT,2)
  ctl$MG_parms["Wtlen_2_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$BwT,2)
  ctl$MG_parms["VonBert_K_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.F$k,2)
  ctl$MG_parms["VonBert_K_Fem_GP_1", c("LO","HI")]=c(life.history$Growth.F$k*.25,life.history$Growth.F$k*3)
  ctl$MG_parms["CV_young_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_young,2)
  ctl$MG_parms["CV_old_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_old,2)
  ctl$MG_parms["Mat50%_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$TL.mat.inf.slope[2],2) #life.history$TL.50.mat
  ctl$MG_parms["Mat_slope_Fem_GP_1", c("INIT","PRIOR")]=rep(life.history$TL.mat.inf.slope[1],2)
  ctl$MG_parms["Eggs_alpha_Fem_GP_1", c("LO","INIT","HI","PRIOR")]=c(0,fec.alpha,100,fec.alpha)   
  ctl$MG_parms["Eggs_beta_Fem_GP_1", c("INIT","PRIOR")]=rep(fec.beta,2)
  ctl$MG_parms["FracFemale_GP_1", c("INIT","PRIOR")]=rep(life.history$pup.sx.ratio,2)
  
  #males
  ctl$MG_parms["NatM_p_1_Mal_GP_1", c("INIT","PRIOR")]=rep(Scenario$Mmean,2)
  ctl$MG_parms["NatM_p_1_Mal_GP_1", c("LO","HI")]=c(ctl$MG_parms["NatM_p_1_Mal_GP_1", "INIT"]*.1,ctl$MG_parms["NatM_p_1_Mal_GP_1", "INIT"]*4)
  ctl$MG_parms["L_at_Amin_Mal_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Lzero*a_FL.to.TL+b_FL.to.TL),2)
  ctl$MG_parms["L_at_Amin_Mal_GP_1", c("LO","HI")]=ctl$MG_parms["L_at_Amin_Fem_GP_1", c("LO","HI")]
  ctl$MG_parms["L_at_Amax_Mal_GP_1", c("INIT","PRIOR")]=rep(with(life.history,Growth.M$FL_inf*a_FL.to.TL+b_FL.to.TL),2)
  ctl$MG_parms["L_at_Amax_Mal_GP_1", c("LO","HI")]=round(c(ctl$MG_parms["L_at_Amax_Mal_GP_1", "INIT"]*.25,ctl$MG_parms["L_at_Amax_Mal_GP_1", "INIT"]*2))
  ctl$MG_parms["Wtlen_1_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$AwT.M,2)
  ctl$MG_parms["Wtlen_2_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$BwT.M,2)
  ctl$MG_parms["VonBert_K_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.M$k,2)
  ctl$MG_parms["VonBert_K_Mal_GP_1", c("LO","HI")]=c(life.history$Growth.M$k*.25,life.history$Growth.M$k*3)
  ctl$MG_parms["CV_young_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_young,2)
  ctl$MG_parms["CV_old_Mal_GP_1", c("INIT","PRIOR")]=rep(life.history$Growth.CV_old,2)
  
  #estimate growth params
  if(Scenario$Model=='SS' & life.history$SS3.estim.growth.pars)  
  {
    ctl$MG_parms["L_at_Amax_Fem_GP_1", "PHASE"]=3
    ctl$MG_parms["L_at_Amax_Fem_GP_1", "PR_SD"]=ctl$MG_parms["L_at_Amax_Fem_GP_1", "INIT"]*.15  
    ctl$MG_parms["L_at_Amax_Fem_GP_1", "PR_type"]=6 #0, no prior; 1, symmetric beta; 2, beta; 3, lognormal; 4, lognormal with bias correction; 5, gamma; 6, normal
    
    ctl$MG_parms["L_at_Amax_Mal_GP_1", "PHASE"]=3
    ctl$MG_parms["L_at_Amax_Mal_GP_1", "PR_SD"]=ctl$MG_parms["L_at_Amax_Mal_GP_1", "INIT"]*.15  
    ctl$MG_parms["L_at_Amax_Mal_GP_1", "PR_type"]=6 
    
    ctl$MG_parms["VonBert_K_Fem_GP_1", "PHASE"]=3
    ctl$MG_parms["VonBert_K_Fem_GP_1", "PR_SD"]=ctl$MG_parms["VonBert_K_Fem_GP_1", "INIT"]*.15 
    ctl$MG_parms["VonBert_K_Fem_GP_1", "PR_type"]=6 
    
    ctl$MG_parms["VonBert_K_Mal_GP_1", "PHASE"]=3
    ctl$MG_parms["VonBert_K_Mal_GP_1", "PR_SD"]=ctl$MG_parms["VonBert_K_Mal_GP_1", "INIT"]*.15  
    ctl$MG_parms["VonBert_K_Mal_GP_1", "PR_type"]=6 
    
  }
  
  #M at age
  if(Scenario$Model=='SS' & 'M.at.age'%in%names(Scenario))
  {
    ctl$natM_type=3 #0=1Parm; 1=N_breakpoints; 2=Lorenzen; 3=agespecific; 4=agespec_withseasinterpolate
    if(Scenario$M.at.age=="Mmean.mean.at.age")
    {
      M.age.fem=life.history$Mmean.mean.at.age
      M.age.male=life.history$Mmean.mean.at.age
    }
    if(Scenario$M.at.age=="Mmean.min.at.age")
    {
      M.age.fem=life.history$Mmean.min.at.age
      M.age.male=life.history$Mmean.min.at.age
    }
    natM=data.frame(matrix(c(M.age.fem,M.age.male),nrow=2,byrow = T))
    colnames(natM)=paste('Age',0:unique(life.history$Max.age.F),sep='_')
    ctl$natM=natM
    ctl$MG_parms=ctl$MG_parms[-grep('NatM',rownames(ctl$MG_parms)),]
  }
  #recruitment pars
  ctl$SR_function=3 # 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
  ctl$SR_parms["SR_LN(R0)", c('LO','INIT','HI')]=with(Scenario,c(Ln_R0_min,Ln_R0_init,Ln_R0_max))
  ctl$SR_parms["SR_BH_steep", "INIT"]=Scenario$Steepness
  ctl$SR_parms["SR_BH_steep", "LO"]=0.25
  if(Scenario$Model=='SS')  #estimate h with strong prior (Punt 2023)?
  {
    ctl$SR_parms["SR_BH_steep", "PHASE"]=life.history$Steepness_Phase
    ctl$SR_parms["SR_BH_steep", "PR_SD"]=Scenario$Steepness.sd  
    ctl$SR_parms["SR_BH_steep", "PR_type"]=1 #0, no prior; 1, symmetric beta; 2, beta; 3, lognormal; 4, lognormal with bias correction; 5, gamma; 6, normal
    if(is.null(abundance)) ctl$SR_parms["SR_BH_steep", "PHASE"]=-4
  }
  if(Scenario$Model=='SS') SR_sigmaR=life.history$SR_sigmaR
  ctl$SR_parms["SR_sigmaR", c('LO','HI','INIT')]=c(.01,1,SR_sigmaR) #Spiny dogfish SS assessment
  
  if(Scenario$Model=='SSS') ctl$do_recdev=0  #do_recdev:  0=none; 1=devvector; 2=simple deviations
  if(Scenario$Model=='SS')
  {
    ctl$do_recdev=1 #0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
    RecDev_Phase=life.history$RecDev_Phase
    ctl$recdev_phase=RecDev_Phase
    if(is.null(abundance)) ctl$recdev_phase=1
    ctl$recdev_early_start=0
    ctl$max_bias_adj=0.8
    ctl$min_rec_dev=-2
    ctl$max_rec_dev=abs(ctl$min_rec_dev)
    if(is.null(size.comp) & !is.null(abundance))
    {
      ctl$MainRdevYrFirst=min(abundance$Year)
      ctl$MainRdevYrLast=max(abundance$Year)
      ctl$last_early_yr_nobias_adj=min(abundance$Year)-1
      ctl$first_yr_fullbias_adj=min(abundance$Year)+5 #first year with full bias adj. should be a few years into the data-rich period
    }
    if(!is.null(size.comp))
    {
      ctl$MainRdevYrFirst=min(size.comp$year) #-dat$Nages
      ctl$MainRdevYrLast=max(size.comp$year)
      if(life.history$Name=="smooth hammerhead") ctl$MainRdevYrLast=max(abundance$Year)#to allow fitting Southern2 cpue
      ctl$last_early_yr_nobias_adj=ctl$MainRdevYrFirst-1
      ctl$first_yr_fullbias_adj=min(ctl$MainRdevYrLast-1,ctl$MainRdevYrFirst+5) 
    }
    ctl$last_yr_fullbias_adj=endyr-2
    ctl$first_recent_yr_nobias_adj=endyr   #end_yr_for_ramp_in_MPD
    if(!is.na(life.history$last_early_yr_nobias_adj_in_MPD))
    {
      ctl$last_early_yr_nobias_adj=life.history$last_early_yr_nobias_adj_in_MPD
      ctl$first_yr_fullbias_adj=life.history$first_yr_fullbias_adj_in_MPD
      ctl$last_yr_fullbias_adj=life.history$last_yr_fullbias_adj_in_MPD
      ctl$first_recent_yr_nobias_adj=life.history$first_recent_yr_nobias_adj_in_MPD
      if(life.history$Name=="smooth hammerhead") ctl$first_recent_yr_nobias_adj=max(abundance$Year)#to allow fitting Southern2 cpue
      ctl$max_bias_adj=life.history$max_bias_adj_in_MPD
    }
  }
  
  #fishing mortality pars
  if(Scenario$use_F_ballpark)
  {
    ctl$F_ballpark=Scenario$F_ballpark
    ctl$F_ballpark_year=styr
  }
  
  #Q pars  
  if(Scenario$Model=='SS')
  {
    if(!is.null(abundance))
    {
      ctl$Q_options=ctl$Q_options[which(rownames(ctl$Q_options)%in%dat$fleetinfo$fleetname),]
      ctl$Q_parms=ctl$Q_parms[grep(paste(rownames(ctl$Q_options),collapse='|'),rownames(ctl$Q_parms)),]
       
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
        addNSFQ=ctl$Q_parms[1,]
        rownames(addNSFQ)=paste('LnQ_base_',"Northern.shark(",match('Northern.shark',fleetinfo$fleetname),')',sep='')
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

      Nms=rownames(ctl$Q_parms)
      Nms=gsub(r"{\s*\([^\)]+\)}","",gsub("^.*?base_","",Nms))
      rownames(ctl$Q_parms)=Nms
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
      n.indices=nrow(ctl$Q_options)
      Indx.small.CV=dat$CPUE%>%group_by(index)%>%summarise(Mean.CV=mean(CV))  
      if(life.history$Name%in%c("spinner shark","tiger shark")) Indx.small.CV$Mean.CV=default.CV*.9 #need extra_Q for both indices to fit spinner cpue
      if(!Scenario$extra.SD.Q=='always' & n.indices>1) Indx.small.CV=Indx.small.CV%>%filter(Mean.CV<default.CV)%>%pull(index)
      if(Scenario$extra.SD.Q=='always') Indx.small.CV=Indx.small.CV$index
      if(length(Indx.small.CV)>0)
      {
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
        iiq_block=match(paste('fleet_',ctl$Q_options$fleet[match('Southern.shark_1',rownames(ctl$Q_options))],sep=''),rownames(ctl$Q_parms))
        if(is.na(iiq_block)) iiq_block=match('Southern.shark_1',rownames(ctl$Q_parms))
        if(!is.na(iiq_block)) ctl$Q_parms[iiq_block,c('Block','Block_Fxn')]=c(life.history$Q_param_Block,life.history$Q_param_Blk_Fxn)
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
    ddumy=ctl$size_selex_types%>%
      rownames_to_column('fleetname')%>%
      mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                              ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                                     fleetname)))%>%
      filter(fleetname%in%dis.flits)%>%
      mutate(Fleet=row_number())
    rownames(ddumy)=ddumy$fleetname
    if(!'Northern.shark'%in%rownames(ddumy))
    {
      ddumy=ddumy%>%
        mutate(Pattern=ifelse(fleetname=='Other',24,Pattern),
               Special=ifelse(fleetname=='Other',0,ifelse(fleetname=='Southern.shark_2',2,Special)))
    }
    if(!is.null(F.tagging))
    {
      add.F.series.dummy=ddumy[grep('Southern.shark_2',ddumy$fleetname),]%>%
        mutate(Fleet=id.F.series,
               fleetname=dat$fleetinfo$fleetname[id.F.series])
      rownames(add.F.series.dummy)=dat$fleetinfo$fleetname[id.F.series]
      ddumy=rbind(ddumy,add.F.series.dummy)
    }
    
    #change to logistic if required
    Logis.sel=life.history$SS_selectivity%>%filter(Fleet=='Northern.shark')
    if(all(is.na(Logis.sel[,c('P_3','P_4','P_5','P_6')])))
    {
      ddumy[ddumy$fleetname=='Northern.shark','Pattern']=1
    }
    Logis.sel=life.history$SS_selectivity%>%filter(Fleet=='Survey')
    if(all(is.na(Logis.sel[,c('P_3','P_4','P_5','P_6')])))
    {
      ddumy[ddumy$fleetname=='Survey','Pattern']=1
    }
    #change tiger shark 'other' as unlikely to be landing very large sharks
    if(life.history$Name=="tiger shark")
    {
      ddumy[ddumy$fleetname=="Other",c('Pattern','Special')]=c(24,0)
    }
    
    ctl$size_selex_types=ddumy%>%dplyr::select(-fleetname,-Fleet)
    
    Sel.ptrn=ctl$size_selex_types$Pattern
    names(Sel.ptrn)=paste('Fishery',1:length(Sel.ptrn),sep='')
    
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
    
    #turn on Southern.shark_2 if size compo data  
    if(!is.null(size.comp))
    {
      Tab.size.comp.dat=with(dat$lencomp,table(Fleet))
      #Tab.size.comp.dat=with(dat$lencomp%>%filter(Sex==1),table(Fleet))
      names(Tab.size.comp.dat)=fleetinfo$fleetname[as.numeric(names(Tab.size.comp.dat))]
      nn.Southern.shark_2=subset(Tab.size.comp.dat,names(Tab.size.comp.dat)=="Southern.shark_2")
      if(length(nn.Southern.shark_2)>0)
      {
        #if(nn.Southern.shark_2>1)
        #{
        ctl$size_selex_types[rownames(ctl$size_selex_types)=="Southern.shark_2",]=ctl$size_selex_types[rownames(ctl$size_selex_types)=="Southern.shark_1",]
        add.Southern.shark_2.pars=ctl$size_selex_parms[grepl('Southern.shark_1',rownames(ctl$size_selex_parms)),]
        rownames(add.Southern.shark_2.pars)=str_replace(rownames(add.Southern.shark_2.pars), "k_1", "k_2")
        ctl$size_selex_parms=rbind(ctl$size_selex_parms[!grepl('Survey',rownames(ctl$size_selex_parms)),],
                                   add.Southern.shark_2.pars,
                                   ctl$size_selex_parms[grepl('Survey',rownames(ctl$size_selex_parms)),])
        #} 
      }
    }
    #change tiger shark 'other'
    if(life.history$Name=="tiger shark")
    {
      add.this=ctl$size_selex_parms[which(row_nm_size_selex_parms%in%dis.flits),]
      add.Other=add.this[which(row_nm_size_selex_parms=='Northern.shark'),]
      rownames(add.Other)=str_replace(rownames(add.Other), "Northern.shark", "Other")
      ctl$size_selex_parms=rbind(add.this[grep('Northern.shark',rownames(add.this)),],add.Other,
                                 add.this[-grep('Northern.shark',rownames(add.this)),])
      
    }
    
    #allocated species specific values to sel pars
    Mirrored.sels=rownames(ctl$size_selex_types%>%filter(Pattern==15))
    life.history$SS_selectivity=life.history$SS_selectivity%>%
      filter(Fleet%in%dis.flits)
    if(length(Mirrored.sels)>0) life.history$SS_selectivity=life.history$SS_selectivity%>%filter(!Fleet%in%Mirrored.sels)
    id.fleets=fn.get.in.betwee(x=rownames(ctl$size_selex_parms))
    pis=unique(id.fleets)
    for(px in 1:length(pis))
    {
      iid=which(id.fleets==pis[px])
      this.par=life.history$SS_selectivity[,match(pis[px],colnames(life.history$SS_selectivity))]
      ctl$size_selex_parms[iid,"INIT"]=this.par  
      
      multiplr=rep(0.1,length(this.par))
      multiplr=ifelse(this.par<0,2,multiplr)
      low.bound=multiplr*ctl$size_selex_parms[iid,"INIT"]
      if(pis[px]=="P_1")
      {
        bump=1.1 #1.25
        low.bound=sapply(low.bound, function(x) max(dat$minimum_size*bump,x))
        if(!is.null(size.comp))
        {
          low.bound=sapply(low.bound, function(x) max(min(dat$lbin_vector),x))
          if(life.history$Name=="tiger shark") low.bound=min(dat$lbin_vector) 
        }
        
      }
      ctl$size_selex_parms[iid,"LO"]=low.bound
      
      multiplr=rep(2,length(this.par))
      multiplr=ifelse(this.par<0,-2,multiplr)
      up.bound=multiplr*ctl$size_selex_parms[iid,"INIT"]
      if(pis[px]=="P_1")
      {
        up.bound=sapply(up.bound, function(x) min(dat$maximum_size*.975,x))
      }
      
      ctl$size_selex_parms[iid,"HI"]=up.bound
    }
    ctl$size_selex_parms=ctl$size_selex_parms%>%filter(!is.na(INIT))
    
    #set phases for estimable selectivities
    if(is.null(size.comp))ctl$size_selex_parms[,"PHASE"]=-2
    if(!is.null(size.comp))
    {
      flit.size.comp.obs=sort(unique(dat$lencomp$Fleet))
      flit.no.size.comp.obs=fleetinfo%>%filter(!fleetname%in%rownames(dat$len_info)[flit.size.comp.obs])%>%pull(fleetname)
      if(length(Mirrored.sels)>0) flit.no.size.comp.obs=subset(flit.no.size.comp.obs,!flit.no.size.comp.obs%in%Mirrored.sels)
      if(length(flit.no.size.comp.obs)>0)
      {
        for(px in 1:length(flit.no.size.comp.obs))
        {
          iid=grep(flit.no.size.comp.obs[px],rownames(ctl$size_selex_parms))
          ctl$size_selex_parms[iid,]$PHASE=-2
        }
      }
    }
    if(life.history$Name=="tiger shark")
    {
      ctl$size_selex_parms[grep(paste(c('P_6_Northern.shark','SizeSel_P_6_Survey'),collapse='|'),rownames(ctl$size_selex_parms)),'PHASE']=4
    }
    if(any(is.null(abundance) | life.history$Name=="spinner shark"))   #fixed most sel pars if no abundance (SS-CL)  
    {
      if(life.history$Name=="spinner shark")
      {
        ctl$size_selex_parms$PHASE=-abs(ctl$size_selex_parms$PHASE)
        ctl$size_selex_parms$PHASE[grep('P_3_Southern.shark_1',rownames(ctl$size_selex_parms))]=2
      }
    }
    if(life.history$Name=="whiskery shark")
    {
      ctl$size_selex_parms$PHASE[grep('P_2_Southern.shark',rownames(ctl$size_selex_parms))]=-4
    }
    
    #turn on Southern.shark_1 if available meanbodywt (because Southern.shark_2 mirrors Southern.shark_1) 
    if(!is.null(meanbodywt) & !is.null(size.comp))
    {
      xx=size.comp%>%filter(Fleet%in%c(3:4))
      if(nrow(xx)==0)   #only need to turn of if there is no size comp for Southern.shark_1 or Southern.shark_2
      {
        iiD=grepl('P_1_Southern.shark_1',rownames(ctl$size_selex_parms))
        ctl$size_selex_parms[iiD,"PHASE"]=abs(ctl$size_selex_parms[iiD,"PHASE"])
        ctl$size_selex_parms[iiD,'PRIOR']=ctl$size_selex_parms[iiD,c('INIT')]
        ctl$size_selex_parms[iiD,'PR_SD']=ctl$size_selex_parms[iiD,c('INIT')]*.2
        ctl$size_selex_parms[iiD,'PR_type']=6 #0, no prior; 1, symmetric beta; 2, beta; 3, lognormal; 4, lognormal with bias correction; 5, gamma; 6, normal
        if(life.history$Name=="tiger shark")
        {
          ctl$size_selex_parms[iiD,'PR_type']=0
          iiD=grepl(paste(c('P_4_Southern.shark_1'),collapse='|'),rownames(ctl$size_selex_parms))
          ctl$size_selex_parms[iiD,"PHASE"]=abs(ctl$size_selex_parms[iiD,"PHASE"])
        }
      }
    }
    
    
    #2. age_selex
    ddumy=ctl$age_selex_types%>%
      rownames_to_column('fleetname')%>%
      mutate(fleetname=ifelse(fleetname=='Southern.shark_monthly','Southern.shark_1',
                              ifelse(fleetname=='Southern.shark_daily','Southern.shark_2',
                                     fleetname)))%>%
      filter(fleetname%in%dis.flits)%>%
      mutate(Fleet=row_number())
    rownames(ddumy)=ddumy$fleetname
    ddumy$Pattern=life.history$age_selex_pattern
    ctl$age_selex_types=ddumy%>%dplyr::select(-fleetname,-Fleet)
    
    if(!is.null(F.tagging))
    {
      F.age.sel.pat=ctl$age_selex_types[match('Southern.shark_2',rownames(ctl$age_selex_types)),]
      rownames(F.age.sel.pat)=F.fleet
      ctl$age_selex_types=rbind(ctl$age_selex_types,F.age.sel.pat)
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
  
  #set prior to init value 
  ctl$SR_parms[,"PRIOR"]=ctl$SR_parms[,"INIT"]
  ctl$MG_parms[,"PRIOR"]=ctl$MG_parms[,"INIT"]
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
    ctl$lambdas=Like_comp
    ctl$N_lambdas=nrow(ctl$lambdas)
  }
  
  
  #3.3. starter file
  start$datfile='data.dat'
  start$ctlfile='control.ctl'
  
  #3.4. forecast file
  fore$Fcast_selex=0
  #fore$    #consider updating reference points, years, etc   
  
  
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
fn.med.ci.SS=function(d)
{
  dframe=data.frame(year=as.numeric(names(d)),
                    median=apply(d,2,function(x) median(x)),
                    lower.95=apply(d,2,function(x) quantile(x,probs=0.025)),
                    upper.95=apply(d,2,function(x) quantile(x,probs=0.975)))
  
  return(dframe)
}
fn.MC.sims=function(this.wd1,nMC=nMCsims,arg=Arg.no.estimation,B.th,scen)   
{
  #Get var-covar matrix
  MLE=read.admbFit(paste(this.wd1,'ss',sep='/'))
  n.mle=1:MLE$nopar
  Nms=MLE$names[n.mle]
  
  #Monte Carlo sampling
  Depletion=vector('list',nMC)
  F.series=B.Bmsy=F.Fmsy=Depletion
  New.dir=paste(this.wd1,'MonteCarlo',sep='/')
  if(!dir.exists(New.dir))dir.create(New.dir)
  for(n in 1:nMC)   
  {
    #1. Sample from multivariate distribution
    Rand.par=c(rmvnorm(1,mean=MLE$est[n.mle],sigma=MLE$cov[n.mle,n.mle]))
    names(Rand.par)=Nms
    
    #2. Create new SS input files
    copy_SS_inputs(dir.old = this.wd1, dir.new = New.dir,overwrite = TRUE)
    ctl <- r4ss::SS_readctl(file = file.path(this.wd1, "control.ctl"))
    id=which(ctl$SR_parms$PHASE>0)
    if(length(id)>0)
    {
      dispars=grep('SR_parm',Nms)
      if(length(dispars)>0) ctl$SR_parms$INIT[id]=Rand.par[dispars]
    }
    id=which(ctl$Q_parms$PHASE>0)
    if(length(id)>0)
    {
      dispars=grep('Q_parm',Nms)
      if(length(dispars)>0) ctl$Q_parms$INIT[id]=Rand.par[dispars]
    }
    id=which(ctl$size_selex_parms$PHASE>0)
    if(length(id)>0)
    {
      dispars=grep('selparm',Nms)
      if(length(dispars)>0) ctl$size_selex_parms$INIT[id]=Rand.par[dispars]
    }
    r4ss::SS_writectl(ctl, outfile = file.path(New.dir, "control.ctl"), overwrite = TRUE, verbose = FALSE)
    
    #3. Run SS without estimation
    fn.run.SS(where.inputs=New.dir,
              where.exe=handl_OneDrive('SS3/ss_win.exe'),
              args=arg)  
    
    #4. Store quantities of interest (depletion)
    Report=SS_output(New.dir,covar=F,forecast=F,readwt=F,checkcor=F)
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
#note: this follows Carvalho et al 2021
fn.fit.diag_SS3=function(WD,disfiles,R0.vec,exe_path,start.retro=0,end.retro=5,
                         do.like.prof=FALSE,do.retros=FALSE,do.jitter=FALSE,numjitter)
{
  Report=SS_output(dir=WD,covar=T)
  
  dirname.diagnostics <- paste(WD,"Diagnostics",sep='/')
  if(!dir.exists(dirname.diagnostics)) dir.create(path=dirname.diagnostics, showWarnings = TRUE, recursive = TRUE)
  
  cpue.series=length(unique(Report$cpue$Fleet))
  length.series=length(unique(unique(Report$len_comp_fit_table$Fleet)))
  
  #1. Goodness-of-fit diagnostic
  #1.1. residuals with smoothing function showing trends
  tiff(file.path(dirname.diagnostics,"jabbaresidual.tiff"),
       width = 1000, height = 2000,units = "px", res = 300, compression = "lzw")
  sspar(mfrow=c(2,1),labs=T,plot.cex=0.9)
  ss3diags::SSplotJABBAres(ss3rep=Report,subplots = "cpue",add=TRUE)
  ss3diags::SSplotJABBAres(ss3rep=Report,subplots = "len",add=TRUE)
  dev.off()
  
  #1.2. runs test
  dis.dat=c("cpue")
  if(any(unique(Report$lendbase$Yr)%in%Report$catch$Yr[nrow(Report$catch)-end.retro+1]:
         Report$catch$Yr[nrow(Report$catch)])) dis.dat=c(dis.dat,"len")      
  if(any(unique(Report$agedbase$Yr)%in%Report$catch$Yr[nrow(Report$catch)-end.retro+1]:
         Report$catch$Yr[nrow(Report$catch)])) dis.dat=c(dis.dat,"age")  
  tiff(file.path(dirname.diagnostics,"runs_tests.tiff"),
       width = 2000, height = 2000,units = "px", res = 300, compression = "lzw")
  sspar(mfrow=n2mfrow(length(dis.dat)),labs=T,plot.cex=0.9)
  for(pp in 1:length(dis.dat)) ss3diags::SSplotRunstest(ss3rep=Report,subplots = dis.dat[[pp]],add=TRUE)
  dev.off()
  
  runs.test.value=vector('list',length(dis.dat))
  for(pp in 1:length(dis.dat))  runs.test.value[[pp]]=SSrunstest(Report,quants =  dis.dat[[pp]])
  write.csv(do.call(rbind,runs.test.value),paste(dirname.diagnostics,"runs_tests.csv",sep='/'),row.names = F)
  
  
  #2. Model consistency  
  #2.1. Likelihood profile on Ro
  if(do.like.prof) # very time consuming, takes 2.5 minutes per R0 value
  {
    library(doParallel)
    registerDoParallel(detectCores()-1)
    
    # Step 1. Identify a directory for the profile likelihood model run(s)
    dirname.base <- paste(dirname.diagnostics,"R0_profile",sep='/')
    
    # Step 2. Identify a directory where the completed base model run is located
    dirname.completed.model.run<-WD
    
    # Step 3. Create a "R0_profile" subdirectory and set as the working directory
    dirname.R0.profile <- dirname.base
    if(!dir.exists(dirname.R0.profile)) dir.create(path=dirname.R0.profile, showWarnings = TRUE, recursive = TRUE)
    
    mydir <- dirname.R0.profile
    setwd(dirname.R0.profile)
    
    
    # Step 4. Create a "Figures_Tables" subdirectory
    plotdir=paste0(dirname.R0.profile, "/Figures & Tables")
    if(!dir.exists(plotdir)) dir.create(path=plotdir, showWarnings = TRUE, recursive = TRUE)
    
    
    # Step 5. Create a "Reference_run" subdirectory and copy completed base model output to this directory
    #reference.dir <- paste0(mydir,'/Reference_run') 
    #if(!dir.exists(reference.dir)) dir.create(path=reference.dir, showWarnings = TRUE, recursive = TRUE)
    #file.copy(from=Sys.glob(paste(dirname.completed.model.run, "*.*", sep="/"), dirmark = FALSE),
    #          to=reference.dir)
    #for(nn in disfiles){file.copy(paste(dirname.completed.model.run,"/", nn, sep='')  ,     to=reference.dir)}
    
    
    # Step 6. Copy necessary files from the "Reference_run" subdirectory to the "R0_profile" working directory 
    copylst<-disfiles[-match("Report.sso",disfiles)]
    #copylst <-  c("control.ss_new", "data.ss",  "forecast.ss",  "ss.exe", "starter.ss")
    #for(nn in copylst){file.copy(  paste(reference.dir,"/", nn, sep='')  ,     file.path(dirname.R0.profile))}
    for(nn in copylst){file.copy(  paste(WD,"/", nn, sep='')  ,     file.path(dirname.R0.profile))}
    
    # Step 7. Edit "control.ss" in the "R0_profile" working directory to estimate at least one parameter in each phase
    # E.g., 
    control.file <- readLines(paste(dirname.R0.profile, "/control.ss_new", sep=""))
    linen <- NULL
    linen <- grep("#_recdev phase", control.file)
    control.file[linen] <- paste0("1 #_recdev phase")
    write(control.file, paste(dirname.R0.profile, "/control.ss_new", sep=""))
    
    # Step 8. Edit "starter.ss" in the "R0_profile" working directory to read from init values from control_modified.ss
    starter.file <- readLines(paste(dirname.R0.profile, "/starter.ss", sep=""))
    linen <- NULL
    linen <- grep("# 0=use init values in control file; 1=use ss.par", starter.file)
    starter.file[linen] <- paste0("0 # 0=use init values in control file; 1=use ss.par")
    write(starter.file, paste(dirname.R0.profile, "/starter.ss", sep=""))
    ###############
    
    # Step 9. Begin Likelihood profile_R0_example.R
    
    ###Working directory
    setwd(dirname.R0.profile)
    
    ####Set the plotting directory
    plotdir=paste0(dirname.R0.profile, "/Figures & Tables")
    
    
    #########################################################
    ### R0 or any other parameter profile
    #########################################################
    
    # vector of values to profile over
    Nprof.R0 <- length(R0.vec)
    #Define directory
    #mydir <- mydir
    
    #Define the starter file
    starter <- SS_readstarter(file.path(mydir, "starter.ss"))
    
    #Change control file name in the starter file
    starter$ctlfile <- "control_modified.ss" 
    
    # Make sure the prior likelihood is calculated for non-estimated quantities
    starter$prior_like <- 1                                 
    
    SS_writestarter(starter, dir=mydir, overwrite=TRUE)
    
    #Run SS_profile command   
    profile <- profile(dir=mydir, # directory
                       oldctlfile="control.ss_new",
                       newctlfile="control_modified.ss",
                       string="SR_LN(R0)",
                       profilevec=R0.vec,
                       exe=exe_path)
    
    # read the output files (with names like Report1.sso, Report2.sso, etc.)
    prof.R0.models <- SSgetoutput(dirvec=mydir, keyvec=1:Nprof.R0, getcovar = FALSE) # 
    
    # Step 10.  summarize output
    prof.R0.summary <- SSsummarize(prof.R0.models)
    
    # Likelihood components 
    mainlike_components         <- c('TOTAL',"Survey", "Catch", 'Length_comp',
                                     "Age_comp","Mean_body_wt",'Recruitment') 
    
    mainlike_components_labels  <- c('Total likelihood','Index likelihood',"Catch",'Length likelihood',
                                     "Age likelihood","Mean body weight",'Recruitment likelihood') 
    
    # END OPTIONAL COMMANDS
    
    # plot profile using summary created above
    tiff(file.path(dirname.diagnostics,"R0_profile_plot.tiff"),
         width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
    par(mar=c(5,4,1,1))
    
    SSplotProfile(prof.R0.summary,           # summary object
                  profile.string = "R0",     # substring of profile parameter
                  profile.label=expression(log(italic(R)[0])), ymax=150,minfraction = 0.001,
                  pheight=4.5, 
                  print=FALSE, 
                  plotdir=plotdir, 
                  components = mainlike_components, 
                  component.labels = mainlike_components_labels,
                  add_cutoff = TRUE,
                  cutoff_prob = 0.95)
    Baseval <- round(Report$parameters$Value[grep("R0",Report$parameters$Label)],2)
    abline(v = Baseval, lty=2,col='orange',lwd=2)
    legend('bottomright','Base value',lty = 2,col='orange',lwd=2)
    dev.off()
    
    # make timeseries plots comparing models in profile
    labs <- paste("SR_Ln(R0) = ",R0.vec)
    labs[which(round(R0.vec,2)==Baseval)] <- paste("SR_Ln(R0) = ",Baseval,"(Base model)")
    
    SSplotComparisons(prof.R0.summary,legendlabels=labs,
                      pheight=4.5,png=TRUE,plotdir=plotdir,legendloc='bottomleft')
    
    dev.off()
    
    ###Piner plot
    #Size comp
    tiff(file.path(dirname.diagnostics,"R0_profile_plot_Length_like.tiff"),
         width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
    par(mar=c(5,4,1,1))
    PinerPlot(prof.R0.summary, 
              profile.string = "R0", 
              component = "Length_like",
              main = "Changes in length-composition likelihoods by fleet",
              add_cutoff = TRUE,
              cutoff_prob = 0.95)
    abline(v = Baseval, lty=2,col='orange',lwd=2)
    legend('bottomright','Base value',lty = 2,col='orange',lwd=2)
    dev.off()
    
    #Survey
    tiff(file.path(dirname.diagnostics,"R0_profile_plot_Survey_like.tiff"),
         width = 2100, height = 2400,units = "px", res = 300, compression = "lzw")
    par(mar=c(5,4,1,1))
    PinerPlot(prof.R0.summary, profile.string = "R0", component = "Surv_like",main = "Changes in Index likelihoods by fleet",
              add_cutoff = TRUE,
              cutoff_prob = 0.95, legendloc="topleft")
    abline(v = Baseval, lty=2,col='orange',lwd=2)
    legend('bottomright','Base value',lty = 2,col='orange',lwd=2)
    dev.off()
  }
  
  #2.2. Retrospective analysis and Predicting skills
  if(do.retros)      #took 9 mins for 5 years
  {
    dirname.Retrospective <- paste(dirname.diagnostics,"Retrospective",sep='/')
    if(!dir.exists(dirname.Retrospective)) dir.create(path=dirname.Retrospective, showWarnings = TRUE, recursive = TRUE)
    plots.Retrospective <- paste(dirname.Retrospective,"Plots",sep='/')
    if(!dir.exists(plots.Retrospective)) dir.create(path=plots.Retrospective, showWarnings = TRUE, recursive = TRUE)
    file.copy(Sys.glob(paste(WD, "*.*", sep="/"), dirmark = FALSE),dirname.Retrospective)
    retro(dir = dirname.Retrospective,
          years = start.retro:-end.retro,
          exe=exe_path)
    retroModels <- SSgetoutput(dirvec = file.path(dirname.Retrospective, "retrospectives", paste("retro", start.retro:-end.retro, sep = "")))
    retroSummary <- SSsummarize(retroModels,verbose = FALSE)
    if("len"%in%dis.dat) retroComp= SSretroComps(retroModels)
    #endyrvec <- retroSummary[["endyrs"]] + start.retro:-end.retro
    #SSplotComparisons(summaryoutput=retroSummary,
    #                  subplots= c(2,4,6,12,14),
    #                  endyrvec = endyrvec,
    #                  png=TRUE,
    #                  plotdir=plots.Retrospective,
    #                  legendlabels = paste("Data", start.retro:-end.retro, "years"))
    
    tiff(file.path(dirname.diagnostics,"retro_Mohns_Rho.tiff"),
         width = 2000, height = 1800,units = "px", res = 300, compression = "lzw")
    sspar(mfrow=c(2,2),plot.cex=0.8)
    rb.full = SSplotRetro(retroSummary,add=T,forecast = F,legend = F,verbose=F)
    rf.full = SSplotRetro(retroSummary,add=T,subplots="F", ylim=c(0,0.4),
                          forecast = F,legendloc="topleft",legendcex = 0.8,verbose=F)
    
    rb = SSplotRetro(retroSummary,add=T,forecast = T,legend = F,verbose=F,xmin=2000)
    rf = SSplotRetro(retroSummary,add=T,subplots="F", ylim=c(0,0.4),
                     forecast = T,legendloc="topleft",legendcex = 0.8,verbose=F,xmin=2000)
    
    dev.off()
    
    #Get MASE as metric of prediction skill  
    hcI = SSmase(retroSummary)
    out.MASE=hcI
    if(exists('retroComp'))
    {
      hcL = SSmase(retroComp,quants = "len")
      out.MASE=rbind(out.MASE,hcL)
    }
    write.csv(out.MASE,paste(dirname.diagnostics,"retro_hcxval_MASE.csv",sep='/'),row.names = F)
    write.csv(SShcbias(retroSummary),paste(dirname.diagnostics,"retro_Mohns_Rho.csv",sep='/'),row.names = F)
    
    
    #Hindcasting cross validation
    tiff(file.path(dirname.diagnostics,"Hindcasting cross-validation.tiff"),
         width = 2000, height = 1800,units = "px", res = 300, compression = "lzw")
    sspar(mfrow=n2mfrow(length(dis.dat)),plot.cex=0.8)
    hci = SSplotHCxval(retroSummary,add=T,verbose=F,ylimAdj = 1.3,legendcex = 0.7)
    if(exists('retroComp'))hci = SSplotHCxval(retroComp,subplots="len",add=T,verbose=F,ylimAdj = 1.3,legendcex = 0.7)
    dev.off()
  }
  
  
  #3. Convergence
  #3.1. Jittering 
  if(do.jitter)
  {
    #Create directories
    dirname.Jitter <- paste(dirname.diagnostics,"Jitter",sep='/')
    if(!dir.exists(dirname.Jitter)) dir.create(path=dirname.Jitter, showWarnings = TRUE, recursive = TRUE)
    file.copy(Sys.glob(paste(WD, "*.*", sep="/"), dirmark = FALSE),dirname.Jitter)
    
    #Run jitter   
    jit.likes <- r4ss::jitter(dir = dirname.Jitter,
                              Njitter = numjitter,
                              jitter_fraction = 0.1,
                              exe=exe_path)
    
    #Read in results using other r4ss functions
    # (note that un-jittered model can be read using keyvec=0:numjitter)
    profilemodels <- SSgetoutput(dirvec = dirname.Jitter, keyvec = 0:numjitter, getcovar = FALSE) #0 is basecase
    profilesummary <- SSsummarize(profilemodels)
    
    Total.likelihoods=profilesummary[["likelihoods"]][1, -match('Label',names(profilesummary[["likelihoods"]]))]
    Params=profilesummary[["pars"]]
    
    tiff(file.path(dirname.diagnostics,"Jitter.tiff"), 
         width = 2100, height = 2100,units = "px", res = 300, compression = "lzw")
    
    sspar(mfrow=c(2,1),labs=T,plot.cex=0.9)
    
    dis.cols=grep('replist',colnames(profilesummary[["SpawnBio"]]))
    YLIM=c(0,max(profilesummary[["SpawnBioUpper"]]$replist0))
    with(profilesummary[["SpawnBioLower"]],plot(Yr,replist0,ylim=YLIM,col='transparent',ylab='SSB (t)',xlab='Year'))
    polygon(x=c(profilesummary[["SpawnBioLower"]]$Yr,rev(profilesummary[["SpawnBioLower"]]$Yr)),
            y=c(profilesummary[["SpawnBioLower"]]$replist0,rev(profilesummary[["SpawnBioUpper"]]$replist0)),
            col='grey60')
    CL=rainbow(length(dis.cols))
    for(x in dis.cols)
    {
      lines(profilesummary[["SpawnBio"]]$Yr,profilesummary[["SpawnBio"]][,x],col=CL[x])
    }
    
    plot(1:length(Total.likelihoods[-1]),Total.likelihoods[-1],
         xlab='Jitter runs at a converged solution',ylab='Total likelihood')
    abline(h=Total.likelihoods[1],lty=2,col=2)
    points(1:length(Total.likelihoods[-1]),Total.likelihoods[-1],cex=2,pch=19)
    
    dev.off()
    
    
    
  }
  
  
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
  
