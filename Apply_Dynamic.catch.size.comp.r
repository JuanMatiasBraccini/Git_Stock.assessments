#notes: single-area, two-sexes, size-structured integrated model fitted to catch and size composition
#       Selectivity is assumed to be known.
#       Uncertainty derived from resampling variance-cov matrix  
#       Only applicable to species with representative size comp samples and a history of exploitation.


sourceCpp(handl_OneDrive("/Analyses/Population dynamics/Other people's code/Alex dynamic catch & size/DynamicCatchLenModfunctions.cpp"))

# Function for nlmb (which can only have one input, i.e. parameter list)
ObjFunc <- function(params)
{
  res=ObjectiveFunc_cpp(params, nYears, nLen_SimYrs, Len_SimYr, ObsLenCompVec1, ObsLenCompVec2, ObsLenCompVec3, 
                        ObsAnnCatch, MaxAge, nLenCl, midpt, LTM_Fem, LTM_Mal, WtAtLen, FemMatAtLen, RecLenDist, Init_F, 
                        PropFemAtBirth, SelAtLength, NatMort_mean, NatMort_sd, Steepness_mean, Steepness_sd,
                        lnSigmaR)
  return(res)
}

#1. TDGDLF size composition
Dynamic_catch_size.comp_TDGDLF=vector('list',N.sp)
names(Dynamic_catch_size.comp_TDGDLF)=Keep.species
system.time({for(l in 1: N.sp)  
{
  # Calculate F
  Outfile='TDGDLF' 
  this.size.comp=paste('Size_composition',c('West.6.5','West.7','Zone1.6.5','Zone1.7','Zone2.6.5','Zone2.7'),sep="_")
  outfile=paste(Outfile,'_histogram',sep='')
  
  #get size composition  
  iid=Species.data[[l]][fn.extract.dat(this.size.comp,names(Species.data[[l]]))]
  if(length(iid)>0 & names(Species.data)[l]%in%Lista.sp.outputs$Other.sp)
  {
    #get size composition data
    dummy=do.call(rbind,iid)
    if(grepl("TDGDLF",outfile))
    {
      dummy=dummy%>%
        mutate(dummy=sub(".*Size_composition_", "", rownames(dummy)),
               Mesh=word(dummy,2,sep = "\\."),
               Mesh=ifelse(Mesh=='6','6.5',Mesh),
               Mesh=factor(Mesh,levels=c('6.5','7')),
               Zone=word(dummy,1,sep = "\\."))
      
    }
    dummy=dummy%>%filter(year<=as.numeric(substr(Last.yr.ktch,1,4)))
    
    #keep years with minimum sample size
    N.min=dummy%>%
      group_by(FINYEAR)%>%
      tally()%>%
      filter(n>=Min.annual.obs.ktch)%>%
      mutate(Keep=FINYEAR)
    if(nrow(N.min)>0)
    {
      print(paste("Dynamic catch and size comp model with dome-shape selectivity for --",names(Species.data)[l],"---",Outfile))
      
      dummy=dummy%>%
        mutate(Keep=FINYEAR,
               year.f=as.numeric(substr(FINYEAR,1,4)))%>%
        filter(Keep%in%N.min$Keep)%>%
        mutate(TL=FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)     #calculate total length
      
      
      #1. Plot observed size frequency by year and mesh for years with minimum sample size
      if(grepl("TDGDLF",outfile))
      {
        p=dummy%>%
          ggplot( aes(x=TL, color=Mesh, fill=Mesh)) +
          geom_histogram(alpha=0.6, binwidth = TL.bins.cm)
        WHERE="top"
      }else
      {
        p=dummy%>%
          ggplot( aes(x=TL,color=year, fill=year)) +
          geom_histogram(alpha=0.6, binwidth = TL.bins.cm)
        WHERE="none"
      }
      p=p+
        facet_wrap(~FINYEAR,scales='free_y')+
        xlab("Total length (cm)")+ylab("Count")+
        theme(legend.position=WHERE,
              legend.title = element_blank(),
              legend.text=element_text(size=14),
              strip.text.x = element_text(size = 12),
              axis.text=element_text(size=12),
              axis.title=element_text(size=16))
      print(p)
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(List.sp[[l]]$Name),"/",AssessYr,
                   "/1_Inputs/Visualise data/Size.comp_Dynamic catch and size comp.",outfile,".tiff",sep=''),
             width = 8,height = 8,compression = "lzw")
      
      
      #2. Fit model and derive quantities of interest
      
      #Fixed input pars
      lnSigmaR = 0.1 # specifying fairly low value here (as might be assumed for sharks)
      Init_F = 0.001 # Currently specified parameter
      SDGrowthRec = c(20, 20) # sd for initial size distns. for female and male recruits
      CVLenAtAge = c(0.1, 0.1)
      lnRecDev=0  # mean recruitment deviation across years
      
      
      #Life history   
      MaxAge = ceiling(mean(List.sp[[l]]$Max.age.F))
      Linf = c(List.sp[[l]]$Growth.F$FL_inf*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL,
               List.sp[[l]]$Growth.M$FL_inf*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL) #total length in cm for females and males   
      vbK = c(List.sp[[l]]$Growth.F$k,List.sp[[l]]$Growth.M$k)          # k for females and males
      Lo =  c(List.sp[[l]]$Lzero*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL,
              List.sp[[l]]$Lzero*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)    #total length in cm for females and males
      MaxLen= 10*round(List.sp[[l]]$TLmax/10)
      LenInc=TL.bins.cm
      MatL50=List.sp[[l]]$TL.50.mat
      MatL95=List.sp[[l]]$TL.95.mat
      PropFemAtBirth=List.sp[[l]]$pup.sx.ratio
      wtlen_b=c(List.sp[[l]]$BwT,List.sp[[l]]$BwT.M)   #length-weight pars for females and males
      wtlen_a=c(List.sp[[l]]$AwT,List.sp[[l]]$AwT.M)
      Ages = 1:MaxAge
      nAges = length(Ages)
      NatMort_mean=List.sp[[l]]$Sens.test$SS3$Mmean[1]
      NatMort_sd=sd(colMeans(store.species.M_M.mean[[l]],na.rm=T))
      Steepness_mean=List.sp[[l]]$Sens.test$SS3$Steepness[1]
      Steepness_sd=List.sp[[l]]$Sens.test$SS3$Steepness.sd[1]
      lbnd = seq(0,MaxLen - LenInc, LenInc)
      ubnd = lbnd + LenInc
      midpt = lbnd + (LenInc/2)
      nLenCl = length(midpt)
      
      #Selectivity (combined meshes)
      SelectivityOption = 3 # 1 = gillnet selectivity (Kirkwood and Walker, 1986), 
      # 2 = logistic selectivity, 3 = input vector
      # gillnet selectivity inputs for Kirkwood and Walker (1986) model 
      nMeshes = NA # 2
      MeshSize_mm = c(NA, NA) # c(115, 127)
      theta1 = c(NA, NA) # c(20, 20)
      theta2 = c(NA, NA) # c(200000, 200000)
      # logistic selectivity (F, M)
      SelL50 = c(NA, NA)
      SelL95 = c(NA, NA)
      
      # SelAtLength=Selectivity.at.totalength[[l]]%>%               
      #   mutate(TL=TL)%>%
      #   filter( TL%in%midpt)%>%
      #   pull(Sel.combined)
      SelAtLength=with(List.sp[[l]]$SS_selectivity%>%filter(Fleet=='Southern.shark_1'),doubleNorm24.fn(midpt,a=P_1,
                                      b=P_2, c=P_3, d=P_4, e=P_5, f=P_6,use_e_999=FALSE, use_f_999=FALSE))
      Selectivity <- t(data.frame(FemSelectivity=SelAtLength,
                                  MalSelectivity=SelAtLength)) # selectivity inputted as matrix (F, M)
      
      
      #Catch
      used.ktch='fishery'
        #total
      if(used.ktch=='Total')
      {
        Katch=ktch.combined%>%
          filter(Name==names(Dynamic_catch_size.comp_TDGDLF)[l])%>%
          ungroup()%>%
          dplyr::select(finyear,Tonnes)
      }
        #TDGLDF only
      if(used.ktch=='fishery')
      {
        Katch=KtCh%>%
          filter(FishCubeCode%in%c("Discards_TDGDLF","JASDGDL","WCDGDL","Historic"))%>%
          filter(Name==names(Dynamic_catch_size.comp_TDGDLF)[l])%>%
          group_by(finyear)%>%
          summarise(Tonnes=sum(LIVEWT.c))%>%
          ungroup()
      }
      ObsAnnCatch=Katch$Tonnes
      nYears=length(ObsAnnCatch)
      
      
      #Size composition
      n.size.comp=dummy%>%
        group_by(year.f)%>%
        tally()
      if(nrow(n.size.comp)>3)     #select max of 3 years of size compo. Model currently set to 3 years, ask Alex to relax this
      {
        n.size.comp=n.size.comp[sample(1:nrow(n.size.comp),3),]%>%
          arrange(year.f)
      }
      Len_SimYr=match(n.size.comp$year.f,Katch$finyear) # year index for which length sample is taken
      n_SimYrs =length(Len_SimYr)
      nLen_SimYrs = n_SimYrs
      add.dummy=data.frame(bin=midpt)
      ObsCatchFreqAtLen=vector('list',nrow(n.size.comp))
      for(o in 1:length(ObsCatchFreqAtLen))
      {
        x=dummy%>%
          filter(year.f==n.size.comp$year.f[o])%>%
          mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2)%>%
          group_by(bin)%>%
          tally()%>%
          full_join(add.dummy,by='bin')%>%
          arrange(bin)%>%
          mutate(n=ifelse(is.na(n),0,n))%>%
          filter(bin%in%midpt)
        xx=x$n
        names(xx)=x$bin
        ObsCatchFreqAtLen[[o]]=xx
      }
      ObsCatchFreqAtLen=do.call(rbind,ObsCatchFreqAtLen)%>%data.frame
      names(ObsCatchFreqAtLen)=str_remove(colnames(ObsCatchFreqAtLen), "[X]")
      
      ObsLenCompVec1 = rep(0,nLenCl)
      ObsLenCompVec2 = rep(0,nLenCl)
      ObsLenCompVec3 = rep(0,nLenCl)
      ObsLenCompVec1 = as.vector(unlist(ObsCatchFreqAtLen[1,]))
      if (n_SimYrs > 1) ObsLenCompVec2 = as.vector(unlist(ObsCatchFreqAtLen[2,]))
      if (n_SimYrs > 2) ObsLenCompVec3 = as.vector(unlist(ObsCatchFreqAtLen[3,])) 
      
      plot(midpt,SelAtLength)
      lines(midpt,ObsLenCompVec1/max(ObsLenCompVec1))
      lines(midpt,ObsLenCompVec2/max(ObsLenCompVec2))
      
      #Get model inputs into right format
      res=GetBiologyAndFisheryParams_cpp(Ages, nAges, nLenCl, midpt, ubnd, lbnd, 
                                         Lo, Linf, vbK, CVLenAtAge, wtlen_a, wtlen_b, 
                                         MatL50, MatL95, SelectivityOption, Selectivity, theta1, theta2, 
                                         MeshSize_mm, nMeshes, SelL50, SelL95)
      RecLenDist=res$RecLenDist
      LTM_Fem=res$LTM_Fem
      LTM_Mal=res$LTM_Mal
      WtAtLen=res$WtAtLen
      FemMatAtLen=res$FemMatAtLen
      SelAtLength=res$SelAtLength
      
      # Initial values of estimated parameters
      InitRec = List.sp[[l]]$InitRec_dynKtchSize   
      NatMort = NatMort_mean*exp(rnorm(1, 0,NatMort_sd))
      Steepness = min(max(Steepness_mean*exp(rnorm(1, 0,Steepness_sd)),0.2),1)
      params = c(log(InitRec),log(NatMort),log(Steepness),lnRecDev) # Estimated 
      names(params)=c("ln_R0", "ln_M","ln_h","lnRecDev")
      #ObjFunc(params)  
      
      # Fit model (twice)
      InitTime = Sys.time()
      nlmb <- nlminb(params, ObjFunc, gradient = NULL, hessian = TRUE)
      params = nlmb$par
      nlmb <- nlminb(params, ObjFunc, gradient = NULL, hessian = TRUE)
      Duration = Sys.time() - InitTime
      cat("Duration: model fit",Duration,'\n')
      #nlmb$objective
      if(!nlmb$convergence==0)
      {
        print("Model did not converge")
        break
      }
      
      
      # Calculation estimation uncertainty
      hess.out = optimHess(params, ObjFunc)
      vcov.Params = solve(hess.out)
      ses = sqrt(diag(vcov.Params))
      nsims=500
      sims = data.frame(mvrnorm(n = nsims, params, vcov.Params))
      names(sims) = names(params)
      Estimates=data.frame(Parameter=names(sims),
                           Lower.95=apply(sims[,1:4], MARGIN=2, function(x) quantile(x, 0.025,na.rm=T)),
                           Median=apply(sims[,1:4], MARGIN=2, function(x) quantile(x, 0.5,na.rm=T)),
                           Upper.95=apply(sims[,1:4], MARGIN=2, function(x) quantile(x, .975,na.rm=T)))
      
      # Fit diagnostics
      Dyn.ktch.length.path=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                                 capitalize(Keep.species[l]),"/",AssessYr,'/Dyn.ktch.length',sep='')
      if(!dir.exists(Dyn.ktch.length.path))dir.create(Dyn.ktch.length.path)
      
      # 1. Observed vs Expected 
      res=AssessmentModel_cpp(params, nYears, ObsAnnCatch, MaxAge, nLenCl, midpt, 
                              LTM_Fem, LTM_Mal, WtAtLen, FemMatAtLen, RecLenDist, 
                              Init_F, PropFemAtBirth, SelAtLength, lnSigmaR)
      
      pdf(paste(Dyn.ktch.length.path,"/Fit_diagnostics.pdf",sep='')) 
      # catches
      par(mfrow=c(1+nrow(ObsCatchFreqAtLen),1))
      plot(1:nYears,ObsAnnCatch, "b")
      lines(1:nYears,res$Catch_Biom, col=2)
      legend('topleft',c('observed','expected'),bty='n',col=c(1,2),lty=1)
      
      # length composition
      for (i in 1:length(Len_SimYr))
      {
        ExpCatch = res$CatchN[Len_SimYr[i],] / sum(res$CatchN[Len_SimYr[i],])
        plot(midpt,ObsCatchFreqAtLen[i,]/sum(ObsCatchFreqAtLen[i,]), "b",
             main=n.size.comp$year.f[i],ylab='Frequency',xlab='Size class (cm)')
        lines(midpt,ExpCatch, col=2)
      }
      
      # 2. Distribution for estimated parameter 
      par(mfrow=c(2,2))
      sims$R0 = exp(sims$ln_R0)
      sims$M = exp(sims$ln_M)
      sims$h = exp(sims$ln_h)
      hist(sims$R0,main=paste("input R0=",round(InitRec,3)))
      hist(sims$M,main=paste("input M=",round(NatMort,3)))
      hist(sims$h,main=paste("input h=",round(Steepness,3)))
      hist(sims$lnRecDev,main=paste("input lnRecDev=",lnRecDev))
      dev.off()
      
      
      # Calculate quantities of interest   
      BiomEst = data.frame(matrix(nrow=nsims,ncol=nYears))  #female breeding biomass
      colnames(BiomEst) = 1:nYears
      FMortEst = BiomEst       #fishing mortality
      RelBiomEst = BiomEst     #relative female breeding biomass
      for (i in 1:nsims)
      {
        res=AssessmentModel_cpp(as.vector(unlist(sims[i,1:4])),
                                nYears, ObsAnnCatch, MaxAge, nLenCl, midpt, 
                                LTM_Fem, LTM_Mal, WtAtLen, FemMatAtLen, RecLenDist, 
                                Init_F, PropFemAtBirth, SelAtLength, lnSigmaR)
        
        BiomEst[i,] = res$FemSpBiom
        FMortEst[i,] = res$Ann_FMort
        RelBiomEst[i,] = res$FemSpBiom / (res$Unfish_FemSpBiomPerRec * exp(params[1]))
      }
      
      # relative female spawning biomass
      relBiomEst.median = apply(RelBiomEst[,], MARGIN=2, function(x) quantile(x, 0.5,na.rm=T))
      relBiomEst.lowCL = apply(RelBiomEst[,], MARGIN=2, function(x) quantile(x, 0.025,na.rm=T))
      relBiomEst.uppCL = apply(RelBiomEst[,], MARGIN=2, function(x) quantile(x, 0.975,na.rm=T))
      rel.biom=data.frame(year=Katch$finyear,
                          median=relBiomEst.median,
                          lower.95=relBiomEst.lowCL,
                          upper.95=relBiomEst.uppCL,
                          Catch=Katch$Tonnes)  
      
      # fishing mortality
      FMortEst.median = apply(FMortEst[,], MARGIN=2, function(x) quantile(x, 0.5,na.rm=T))
      FMortEst.lowCL = apply(FMortEst[,], MARGIN=2, function(x) quantile(x, 0.025,na.rm=T))
      FMortEst.uppCL = apply(FMortEst[,], MARGIN=2, function(x) quantile(x, 0.975,na.rm=T))
      FMort=data.frame(year=Katch$finyear,
                       median=FMortEst.median,
                       lower.95=FMortEst.lowCL,
                       upper.95=FMortEst.uppCL,
                       Catch=Katch$Tonnes) 
      
      
      #Store quantities
      Dynamic_catch_size.comp_TDGDLF[[l]]=list(Estimates=Estimates,
                                       rel.biom=rel.biom,
                                       FMort=FMort)
      
      plot(rel.biom$year,rel.biom$median,ylim=c(0,max(rel.biom$upper.95)))
      lines(rel.biom$year,rel.biom$lower.95)
      lines(rel.biom$year,rel.biom$upper.95)
      
      rm(dummy,RecLenDist,LTM_Fem,LTM_Mal,WtAtLen,FemMatAtLen,SelAtLength,
         params,sims,nlmb,res,BiomEst,FMortEst,RelBiomEst,
         NatMort,NatMort_sd,Steepness,Steepness_sd,lnRecDev,
         Lo,Linf,vbK,CVLenAtAge,SDGrowthRec,MaxLen,LenInc,MaxAge,
         MatL50,MatL95,PropFemAtBirth,wtlen_a,wtlen_b,
         Len_SimYr,ObsLenCompVec1,
         Katch,lnSigmaR,ObsAnnCatch,
         Estimates,rel.biom,FMort)
    }
  }
}
})