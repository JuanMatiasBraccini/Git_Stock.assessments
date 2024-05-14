#to install latest version of L3Assess
# library(devtools)
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)
# install.packages("C:/~/L3Assess_0.1.0.tar.gz", source = TRUE, repos=NULL) #Install directly from .gz file

library(L3Assess) 
if(!exists('doubleNorm24.fn')) fn.source1("SS_selectivity functions.R")
Alex.Logis.sel=function(pars) 1/(1 + exp(-log(19) * (midpt - pars[1])/(pars[2])))

# Estimated SS selectivity parameters -----------------------------------------------------------------------
List.estimated.SS.selectivity=vector('list',N.sp)
names(List.estimated.SS.selectivity)=Keep.species
for(l in 1:N.sp)
{
  x=match(names(List.estimated.SS.selectivity)[l],names(List.sp))
  List.estimated.SS.selectivity[[l]]=List.sp[[x]]$SS_selectivity
}
#update if re-estimating selectivity with new assessment 
List.estimated.SS.selectivity$`dusky shark`=List.estimated.SS.selectivity$`dusky shark`%>%
                                        mutate(P_1=case_when(Fleet=='Southern.shark_1'~88.646,
                                                             Fleet=='Southern.shark_2'~88.958,
                                                             TRUE~P_1),
                                               P_2=case_when(Fleet=='Southern.shark_1'~-17.942,
                                                             Fleet=='Southern.shark_2'~-18.15,
                                                             TRUE~P_2),
                                               P_4=case_when(Fleet=='Southern.shark_1'~6.524,
                                                             Fleet=='Southern.shark_2'~6.327,
                                                             TRUE~P_4))
List.estimated.SS.selectivity$`sandbar shark`=List.estimated.SS.selectivity$`sandbar shark`%>%
                                      mutate(P_1=case_when(Fleet=='Southern.shark_1'~100.62,
                                                           Fleet=='Southern.shark_2'~100.62,
                                                           TRUE~P_1),
                                             P_2=case_when(Fleet=='Southern.shark_1'~-6.0618,
                                                           Fleet=='Southern.shark_2'~-6.0618,
                                                           TRUE~P_2),
                                             P_3=case_when(Fleet=='Southern.shark_1'~4.166,
                                                           Fleet=='Southern.shark_2'~4.166,
                                                           TRUE~P_3))
List.estimated.SS.selectivity$`whiskery shark`=List.estimated.SS.selectivity$`whiskery shark`%>%
                                            mutate(P_1=case_when(Fleet=='Southern.shark_2'~123.385,
                                                                 TRUE~P_1),
                                                   P_2=case_when(Fleet=='Southern.shark_2'~-3.093,
                                                                 TRUE~P_2),
                                                   P_3=case_when(Fleet=='Southern.shark_2'~4.776,
                                                                 TRUE~P_3),
                                                   P_4=case_when(Fleet=='Southern.shark_2'~5.22,
                                                                 TRUE~P_4))
List.estimated.SS.selectivity$`milk shark`=List.estimated.SS.selectivity$`milk shark`%>%
                                            mutate(P_1=case_when(Fleet=='Northern.shark'~81.123,
                                                                 Fleet=='Survey'~81.063,
                                                                 TRUE~P_1),
                                                   P_2=case_when(Fleet=='Northern.shark'~9.948,
                                                                 Fleet=='Survey'~5.753,
                                                                 TRUE~P_2))
List.estimated.SS.selectivity$`spinner shark`=List.estimated.SS.selectivity$`spinner shark`%>%
                                            mutate(P_1=case_when(Fleet=='Southern.shark_2'~78.72,
                                                                 TRUE~P_1),
                                                   P_2=case_when(Fleet=='Southern.shark_2'~-1.54,
                                                                 TRUE~P_2))
List.estimated.SS.selectivity$`smooth hammerhead`=List.estimated.SS.selectivity$`smooth hammerhead`%>%
                                            mutate(P_1=case_when(Fleet=='Southern.shark_1'~115.79,
                                                                 TRUE~P_1))



# Useful functions -----------------------------------------------------------------------
#getAnywhere()  #see hidden function
#SimLenAndAgeFreqData() #Simulate length and age data, given specified growth, mortality and selectivity
#VisualiseGrowthApplyingLTM() 

# Catch Curve -----------------------------------------------------------------------
Length.catch.curve=vector('list',N.sp)
names(Length.catch.curve)=Keep.species
tic()
for(l in 1: N.sp)
{
  Neim=Keep.species[l]
  SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
  
  #identify fishery with most catch 
  Katch1=KtCh%>%
        filter(Name==Neim)%>%
        mutate(Fishery=case_when(FishCubeCode%in%c('Historic',Southern.shark.fisheries)~"TDGDLF",
                                 FishCubeCode%in%Northern.shark.fisheries~"NSF",
                                 .default="Other"))
  Katch=Katch1%>%
        group_by(Fishery)%>%
        summarise(Tonnes=sum(LIVEWT.c))
  this.size.comp=Katch%>%filter(Tonnes==max(Tonnes))%>%pull(Fishery)
  if(Neim%in%Other.to.NSF) this.size.comp='NSF'
  if(Neim%in%Other.to.TDGDLF) this.size.comp='TDGDLF'
  out.fleet=this.size.comp
  if(any(this.size.comp=='TDGDLF'))
  {
    this.size.comp=paste('Size_composition',c('West.6.5','West.7','Zone1.6.5','Zone1.7','Zone2.6.5','Zone2.7'),sep="_")
  }
  if(any(this.size.comp=='NSF'))
  {
    this.size.comp=paste('Size_composition',c('NSF.LONGLINE'),sep="_")
  }
  
  #get relevant length composition  
  iid=Species.data[[l]][fn.extract.dat(this.size.comp,names(Species.data[[l]]))]
  if(any(grepl('Observations',names(iid)))) iid=iid[-grep('Observations',names(iid))]
  if(length(iid)>0 & !Neim%in%no.empirical.sel.main.fleet)
  {
    for(x in 1:length(iid))
    {
      dd=str_before_first(str_after_nth(names(iid)[x],"_",2), coll(".inch"))
      iid[[x]]=iid[[x]]%>%
        mutate(Zone=str_before_first(dd, coll(".")),  
               Mesh=str_after_first(dd, coll(".")))
    }
    dummy=do.call(rbind,iid)%>%
              filter(year<=as.numeric(substr(Last.yr.ktch,1,4)))
    
    Min.sample.size=Min.annual.obs.ktch*prop.min.N.accepted_other
    if(Neim%in%names(Indicator.species)) Min.sample.size=Min.annual.obs.ktch
    
    #Standardise to same mesh and zone
    if(standardise.mesh.zone)
    {
      #export observations by sex, mesh and zone
      if(any(!is.na(dummy$Mesh)))
      {
        dummy%>%
          mutate(year=as.numeric(substr(FINYEAR,1,4)),
                 year=ifelse(Zone=='West',year-0.15,ifelse(Zone=='Zone2',year+0.15,year)))%>%
          group_by(year,SEX,Zone,Mesh)%>%
          tally()%>%
          ungroup()%>%
          ggplot(aes(year,n))+
          geom_point(aes(color=Zone,shape=Mesh),size=3)+
          facet_wrap(~SEX,ncol=1)+
          theme(legend.position = 'top') 
        ggsave(paste(this.wd,"/Number of length comps observations by sex, mesh and zone.tiff",sep=''),
               width = 8,height = 8,compression = "lzw")
        
      }
        
      dd=Main.zone.mesh%>%filter(Species==Neim)
      if(!is.na(dd$Mesh))
      {
        dummy=dummy%>%filter(Zone==dd$Zone & Mesh==dd$Mesh)
      }
      Min.sample.size=50
    }
    
    #Keep records with min sample size
    N.min=dummy%>%
      group_by(FINYEAR)%>%
      tally()%>%
      filter(n>=Min.sample.size)%>%
      mutate(Keep=FINYEAR)
    if(nrow(N.min)>0)
    {
      this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                    capitalize(Neim),"/",AssessYr,"/Catch curve and per recruit",sep='')
      if(!dir.exists(this.wd))dir.create(this.wd)

      print(paste("Length-based catch curve with input selectivity for --",names(Species.data)[l],'--',out.fleet,'fishery'))
      
      #keep years with min sample size and calculate total length
      dummy=dummy%>%
        mutate(Keep=FINYEAR,
               year.f=as.numeric(substr(FINYEAR,1,4)))%>%
        filter(Keep%in%N.min$Keep)%>%
        mutate(TL=FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)     
      
      #Model inputs
      MaxAge = ceiling(mean(List.sp[[l]]$Max.age.F)) 
      NatMort = List.sp[[l]]$Sens.test$SS3$Mmean[1]   #Hoenig: exp(1.44-0.982*log(MaxAge)); Dureuil: exp(1.551-1.066*log(MaxAge));4.22/MaxAge
      MaxLen = 10*round(List.sp[[l]]$TLmax/10)
      min.TL=with(List.sp[[l]],Lzero*a_FL.to.TL+b_FL.to.TL)
      lbnd = seq(0,MaxLen - LenInc, LenInc)
      ubnd = lbnd + LenInc
      midpt = lbnd + (LenInc/2)
      if(out.fleet=="TDGDLF") SS.flit='Southern.shark_1'
      if(out.fleet=="NSF") SS.flit='Northern.shark'
      if(out.fleet=="Other") SS.flit='Other'
      
      if(used.selectivity=='Estimated by SS')
      {
        dis.Sel=List.estimated.SS.selectivity[[l]]%>%filter(Fleet==SS.flit)
        if(!is.na(dis.Sel$P_3))
        {
          SelectivityVec_full=with(dis.Sel,doubleNorm24.fn(midpt,a=P_1,b=P_2, c=P_3, d=P_4, e=P_5, f=P_6,use_e_999=FALSE, use_f_999=FALSE))
        }
        if(is.na(dis.Sel$P_3))
        {
          SelectivityVec_full=with(dis.Sel,logistic1.fn(midpt,a=P_1,b=P_2))
        }
      }
      if(used.selectivity=='Empirical for selected mesh')  
      {
        dd=Main.zone.mesh%>%filter(Species==Neim)
        if(!is.na(dd$Mesh) & Neim%in%names(Indicator.species))
        {
          dis.mesh=ifelse(dd$Mesh==6.5,'X16.5',ifelse(dd$Mesh==7,'X17.8',NA))
          SelectivityVec_full=Selectivity.at.totalength[[Neim]]%>%
            filter(TL%in%midpt)%>%select(matches(dis.mesh))%>%pull(dis.mesh)
        }
        if(is.na(dd$Mesh) | !Neim%in%names(Indicator.species))
        {
          dis.Sel=List.estimated.SS.selectivity[[l]]%>%filter(Fleet==SS.flit)
          if(!is.na(dis.Sel$P_3))
          {
            SelectivityVec_full=with(dis.Sel,doubleNorm24.fn(midpt,a=P_1,b=P_2, c=P_3, d=P_4, e=P_5, f=P_6,use_e_999=FALSE, use_f_999=FALSE))
          }
          if(is.na(dis.Sel$P_3))
          {
            SelectivityVec_full=with(dis.Sel,logistic1.fn(midpt,a=P_1,b=P_2))
          }
        }
      }
      
      SelectivityVec=SelectivityVec_full
      SelectivityVec[which(lbnd<min.TL)]=1e-20
      
      Linf = c(List.sp[[l]]$Growth.F$FL_inf*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL,
               List.sp[[l]]$Growth.M$FL_inf*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL) #total length in cm for females and males   
      vbK = c(List.sp[[l]]$Growth.F$k,List.sp[[l]]$Growth.M$k)          # k for females and males
      tzero = List.sp[[l]]$Growth.F$to
      GrowthParams = data.frame(Linf=Linf, vbK=vbK)
      InitL50 = List.sp[[l]]$TL.50.mat
      
      
      #Plot length comp,selectivity and catch
      p1=dummy%>%
        mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2)%>%
        group_by(FINYEAR,bin)%>%
        tally()%>%
        group_by(FINYEAR)%>%
        mutate(n1=n/max(n, na.rm=TRUE))%>%
        left_join(N.min%>%rename(N=n)%>%dplyr::select(FINYEAR,N),by='FINYEAR')%>%
        mutate(FINYEAR=paste0(FINYEAR," (n=",N,")"))%>%
        ggplot(aes(bin,n1))+
        geom_bar(stat="identity")+
        facet_wrap(~FINYEAR,scales='free_y')+
        geom_point(data=data.frame(midpt=midpt,SelectivityVec=SelectivityVec),
                  aes(midpt,SelectivityVec),color='red')+
        xlab("Total length (cm)")+ylab("")+
        theme(legend.position="none",
              legend.title = element_blank(),
              legend.text=element_text(size=14),
              strip.text.x = element_text(size = 12),
              axis.text=element_text(size=12),
              axis.title=element_text(size=16))
      
      p2=Katch1%>%
        filter(Fishery==out.fleet)%>%
        group_by(finyear)%>%
        summarise(Tonnes=sum(LIVEWT.c))%>%
        ggplot(aes(finyear,Tonnes))+
        geom_line(linewidth=1.25)+
        xlab("Financial year")+ylab(paste(out.fleet,"catch (tonnes)"))+
        theme(legend.position="none",
              legend.title = element_blank(),
              legend.text=element_text(size=14),
              strip.text.x = element_text(size = 12),
              axis.text=element_text(size=12),
              axis.title=element_text(size=16))+
        geom_text_repel(data=Katch1%>%
                          filter(Fishery==out.fleet & FINYEAR%in%unique(dummy$FINYEAR))%>%
                          group_by(finyear,FINYEAR)%>%
                          summarise(Tonnes=sum(LIVEWT.c)),
                        aes(finyear, Tonnes, label = FINYEAR, colour = factor(FINYEAR)),box.padding = 1,size=5)
      
      ggarrange(p1, p2, ncol=1,nrow=2, heights = c(1,0.5))
      
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(List.sp[[l]]$Name),"/",AssessYr,
                   "/1_Inputs/Visualise data/Size.comp_Catch curve_",out.fleet,".tiff",sep=''),
             width = 8,height = 8,compression = "lzw")
      
      #Plot length comp by zone and mesh selectivity and proportion of mesh
      if(out.fleet=="TDGDLF")
      {
        p1=dummy%>%
          mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2,
                 Zone.Mesh=paste(Zone, Mesh,sep='-'))%>%
          group_by(FINYEAR,bin,Zone.Mesh)%>%
          tally()%>%
          group_by(FINYEAR)%>%
          mutate(n1=n/max(n, na.rm=TRUE))%>%
          left_join(N.min%>%rename(N=n)%>%dplyr::select(FINYEAR,N),by='FINYEAR')%>%
          mutate(FINYEAR=paste0(FINYEAR," (n=",N,")"))%>%
          ggplot(aes(bin,n1,color=Zone.Mesh))+
          geom_line()+
          facet_wrap(~FINYEAR,scales='free_y')+
          geom_point(data=data.frame(midpt=midpt,SelectivityVec=SelectivityVec),
                     aes(midpt,SelectivityVec),color='red',size=.9,alpha=.5)+
          xlab("Total length (cm)")+ylab("")+
          theme(legend.position="top",
                legend.title = element_blank(),
                legend.text=element_text(size=12),
                strip.text.x = element_text(size = 12),
                axis.text=element_text(size=12),
                axis.title=element_text(size=16))+ guides(colour = guide_legend(nrow = 1))
        
        p2=mesh.prop.effort%>%
          filter(Zone=='Combined')%>%
          mutate(Finyr=as.numeric(substr(finyear,1,4)))%>%
          gather(Mesh,Prop,-c(Zone,Finyr,finyear))%>%
          mutate(Mesh=ifelse(Mesh=='X165','6.5',ifelse(Mesh=='X178','7',NA)))%>%
          ggplot(aes(Finyr,Prop,fill=Mesh))+
          geom_bar(stat="identity")+
          xlab("Financial year")+ylab(paste(out.fleet,"Proportional effort out of the 2 meshes"))+
          theme(legend.position="top",
                legend.title = element_blank(),
                legend.text=element_text(size=14),
                strip.text.x = element_text(size = 12),
                axis.text=element_text(size=12),
                axis.title=element_text(size=16)) 
        p3=mesh.prop.effort%>%
          filter(!Zone=='Combined')%>%
          mutate(Finyr=as.numeric(substr(finyear,1,4)))%>%
          gather(Mesh,Prop,-c(Zone,Finyr,finyear))%>%
          mutate(Mesh=ifelse(Mesh=='X165','6.5',ifelse(Mesh=='X178','7',NA)))%>%
          ggplot(aes(Finyr,Prop,fill=Mesh))+
          geom_bar(stat="identity")+
          facet_wrap(~Zone,ncol=1)+
          xlab("Financial year")+ylab("Proportional effort out of the 2 meshes")+
          theme(legend.position="top",
                legend.title = element_blank(),
                legend.text=element_text(size=14),
                strip.text.x = element_text(size = 12),
                axis.text=element_text(size=12),
                axis.title=element_text(size=13)) 
        ggarrange(p1, p3, ncol=1,nrow=2, heights = c(1,0.7))
        
        ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                     capitalize(List.sp[[l]]$Name),"/",AssessYr,
                     "/1_Inputs/Visualise data/Size.comp_Catch curve_by zone & mesh_",out.fleet,".tiff",sep=''),
               width = 8,height = 10,compression = "lzw")
        
      }

      #Get observed length composition for selected period 
      Catch.curve.years=List.sp[[l]]$Catch.curve.years 
      F.estims=vector('list',length(Catch.curve.years))
      names(F.estims)=sapply(Catch.curve.years, paste, collapse="_")
      for(g in 1:length(Catch.curve.years))
      {
        dummy1=dummy%>%filter(FINYEAR%in%Catch.curve.years[[g]]) 
        ObsCatchFreqAtLen=vector('list',1)
        add.dummy=data.frame(bin=midpt)
        x=dummy1%>%
          mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2)%>%
          group_by(bin)%>%
          tally()%>%
          full_join(add.dummy,by='bin')%>%
          arrange(bin)%>%
          mutate(n=ifelse(is.na(n),0,n))%>%
          filter(bin%in%midpt)
        xx=x$n
        names(xx)=x$bin
        ObsCatchFreqAtLen[[1]]=xx
        
        ObsCatchFreqAtLen=do.call(rbind,ObsCatchFreqAtLen)%>%data.frame
        names(ObsCatchFreqAtLen)=str_remove(colnames(ObsCatchFreqAtLen), "[X]")
        
        tiff(file=paste0(this.wd,"/Length comp used in catch curve_",Catch.curve.years[[g]],".tiff"),
             width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
        plot(midpt,ObsCatchFreqAtLen/max(ObsCatchFreqAtLen),type='h',main=paste0(Catch.curve.years[[g]],collapse='; '))
        points(midpt,SelectivityVec,pch=19,col=2)
        dev.off()
        
        #Get F for selected years
        ObsRetCatchFreqAtLen=ObsCatchFreqAtLen[1,]
        ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
        PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
        InitFishMort = Init.F.Ktch.cur # specify starting parameters
        InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
        params = c(InitFishMort_logit)
        if(List.sp[[l]]$estim.logis.sel.catch.curve)
        {
          params=c(InitFishMort_logit,log(quantile(midpt,.7)),log(10))
          SelectivityType=2
          SelectivityVec=NA
        }
        
        FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, 
                                                  SelectivityType, ObsRetCatchFreqAtLen,lbnd, ubnd, midpt, 
                                                  SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, 
                                                  DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
        FittedRes$Finyear=Catch.curve.years[[g]]
        FittedRes$out.fleet=out.fleet
        FittedRes$MaxAge=MaxAge
        FittedRes$MaxLen=MaxLen
        FittedRes$min.TL=min.TL
        FittedRes$Linf=Linf
        FittedRes$vbK=vbK
        FittedRes$tzero=tzero
        FittedRes$CVSizeAtAge=CVSizeAtAge
        FittedRes$NatMort=NatMort
        FittedRes$SelectivityVec=SelectivityVec
        
        F.estims[[g]]=FittedRes
        if(FittedRes$ParamEst[1,1]>0 & !FittedRes$ParamEst[1,2]==0)
        {
          pdf(paste0(this.wd,"/CatchCurve_RetCatch_Mortality_",paste(Catch.curve.years[[g]],collapse='_'),".pdf"))
          PlotLengthBasedCatchCurve_RetCatch(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt,
                                             SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams,
                                             RefnceAges, MaxAge, NatMort, TimeStep, MainLabel=NA,
                                             xaxis_lab=NA, yaxis_lab=NA, xmax=300, xint=50,
                                             ymax=0.15, yint=0.05, PlotCLs=TRUE, FittedRes, nReps=200)
          
          PlotLengthBasedCatchCurve_Mortality(params, DistnType, MLL, SelectivityType, ObsRetCatchFreqAtLen, lbnd, ubnd, midpt, SelectivityVec,
                                              PropReleased, ObsDiscCatchFreqAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges,
                                              MaxAge, NatMort, TimeStep, xmax=NA, xint=NA, ymax=NA, yint=NA, FittedRes)
          dev.off()
        }
        rm(dummy1,ObsCatchFreqAtLen,FittedRes)
      }
      Length.catch.curve[[l]]=F.estims
      
      rm(Neim,dummy,iid,this.size.comp,SelectivityVec,params)
    }
  }
}
toc()


# Yield per Recruit & Spawning Potential Ratio -----------------------------------------------------------------------
#1. Calculate YPR, SPR and Brel
YPR=vector('list',N.sp)
names(YPR)=Keep.species
tic()
for(l in 1:N.sp)  #takes 3.7 secs per nReps per species
{
  Neim=Keep.species[l]
  if(!is.null(Length.catch.curve[[l]]))
  {
    this.wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                  capitalize(Neim),"/",AssessYr,"/Catch curve and per recruit",sep='')
    
    #identify fishery with most catch 
    out.fleet=Length.catch.curve[[l]][[1]]$out.fleet
    
    print(paste("Yield per Recruit for -------------------------------",Neim,'--',out.fleet,'fishery'))
    
    #life history parameters
    MaxModelAge=Length.catch.curve[[l]][[1]]$MaxAge
    MaxLen = Length.catch.curve[[l]][[1]]$MaxLen
    min.TL=Length.catch.curve[[l]][[1]]$min.TL
    lbnd = seq(0,MaxLen - LenInc, LenInc)
    ubnd = lbnd + LenInc
    midpt = lbnd + (LenInc/2)
    nLenCl = length(midpt)
    
    Linf = Length.catch.curve[[l]][[1]]$Linf
    vbK = Length.catch.curve[[l]][[1]]$vbK
    tzero = Length.catch.curve[[l]][[1]]$tzero
    GrowthParams = data.frame(Linf=Linf, vbK=vbK, tzero=tzero)
    CVSizeAtAge = Length.catch.curve[[l]][[1]]$CVSizeAtAge
    lenwt_a <- List.sp[[l]]$AwT/0.001   #from cm-kg to cm-grams
    ln_lenwt_a <- NA # for log-log relationship
    lenwt_b <- List.sp[[l]]$BwT
    EstWtAtLen <- data.frame(EstFemWtAtLen=NA,
                             EstMalWtAtLen=NA) # weight at length, inputted as values in data frame
    
    mat_L50 <- rep(List.sp[[l]]$TL.50.mat,2)
    mat_L95 <- rep(List.sp[[l]]$TL.95.mat,2)
    EstMatAtLen <- data.frame(EstFemMatAtLen=NA,
                              EstMalMatAtLen=NA) # maturity at length, inputted as values in data frame
    sel_L50 <- NA # females, males - Logistic length selectivity relationship parameters
    sel_L95 <- NA # females, males - Logistic length selectivity relationship parameters
    
    NatMort = Length.catch.curve[[l]][[1]]$NatMort
    NatMort_sd <- List.sp[[l]]$Sens.test$SS3$Msd[1]

    Steepness <- List.sp[[l]]$Sens.test$SS3$Steepness[1]
    Steepness_sd <- List.sp[[l]]$Sens.test$SS3$Steepness.sd[1]
    
    #selectivity
    SelectivityVec=Length.catch.curve[[l]][[1]]$SelectivityVec
    if(any(is.na(SelectivityVec))) SelectivityVec= with(Length.catch.curve[[l]][[1]],Alex.Logis.sel(ParamEst[2:3,1]))
    EstGearSelAtLen <- data.frame(EstFemGearSelAtLen=SelectivityVec,
                                  EstMalGearSelAtLen=SelectivityVec)
    ret_Pmax <- NA # maximum retention, values lower than 1 imply discarding of fish above MLL
    ret_L50 <- NA # females, males - Logistic fish retention at length parameters
    ret_L95 <- NA # females, males - Logistic fish retention at length parameters
    EstRetenAtLen <- data.frame(EstFemRetenAtLen=rep(1,nLenCl),
                                EstMalRetenAtLen=rep(1,nLenCl))
    #Current F
    dumi=vector('list',length(Length.catch.curve[[l]]))
    names(dumi)=names(Length.catch.curve[[l]])
    for(g in 1:length(Length.catch.curve[[l]]))
    {
      Current_F = Length.catch.curve[[l]][[g]]$ParamEst[1,1] # estimate of fishing mortality, e.g. from catch curve analysis
      Current_F_sd <- (Length.catch.curve[[l]][[g]]$ParamEst[1,1]-Length.catch.curve[[l]][[g]]$ParamEst[1,2])/1.96
      FMort <- Current_F 
      if(FMort==0) FMort=Dummy.F.mort
      
      #Run Yield per Recruit
      if(FMort>0)
      {
        Res=CalcYPRAndSPRForFMort_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                     RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                     EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                     FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax, 
                                     ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, FMort)
        
        Output_Opt = 1 # 1=standard output, 2=with added length and weight outputs (slower)
        Res=GetPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                    RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                    EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                    FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                    ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, Output_Opt)
        
        Res$F_lim=Res$FishMort[which.min(abs(Res$Eq_FemRelSpBiomResults - Res$BMSY_Lim))]
        Res$F_tar=Res$FishMort[which.min(abs(Res$Eq_FemRelSpBiomResults - Res$BMSY_Targ))]
        
        # explore various per recruit outputs
        pdf(paste0(this.wd,"/CatchCurve_RetCatch_Mortality_",names(Length.catch.curve[[l]])[g],".pdf"))
        
        plot(Res$FishMort, Res$Eq_CatchResults)
        abline(v=Res$F_MSY)
        abline(v=Res$FishMort[which.min(abs(Res$Eq_FemRelSpBiomResults - Res$BMSY_Thresh))],col=2) #another way of searching for Fmsy
        
        PlotOpt <- 0 # 0=all plots, 1=len at-age, 2=wt at length, 3=fem mat/sel/ret at length, 4=mal mat/sel/ret at length,
        # 5=fem F at length, 6=mal F at length, 7=fem rel surv, 8=mal rel surv, 9=fem biom at age, 10=fem biom at age,
        # 11=ypr/eq catch, 12=fem SPR/Brel, 13=mal SPR/Brel, 14=eq recruit
        
        
        PlotPerRecruitResults_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                 RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                 EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                 FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax, 
                                 ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, Current_F, PlotOpt)
        
        PlotOpt <- 1 # 1=females, 2=males, 3=combined sex
        par(mfrow=c(1,1))
        PlotPerRecruit_Biom_no_err_LB(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                      RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                      EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                      FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                      ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, SRrel_Type, NatMort, PlotOpt, RefPointPlotOpt, Current_F)
        
        PlotPerRecruit_Param_Err_Distns(NatMort, NatMort_sd, Current_F, Current_F_sd, Steepness, Steepness_sd)
        

        FittedRes=GetPerRecruitResults_LB_with_err(MaxModelAge, TimeStep, lbnd, ubnd, midpt, nLenCl, GrowthCurveType, GrowthParams,
                                                   RefnceAges, CVSizeAtAge, lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type,
                                                   EstWtAtLen, ReprodScale, ReprodPattern, InitRatioFem, FinalSex_Pmax, FinalSex_L50,
                                                   FinalSex_L95, mat_L50, mat_L95, EstMatAtLen, sel_L50, sel_L95, EstGearSelAtLen, ret_Pmax,
                                                   ret_L50, ret_L95, EstRetenAtLen, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                                   Current_F, Current_F_sd, nReps)
        
        PlotPerRecruit_Biom_with_err_LB(MaxModelAge, TimeStep, Linf, vbK, tzero, EstLenAtAge,
                                        lenwt_a, ln_lenwt_a, lenwt_b, WLrel_Type, EstWtAtAge, ReprodScale, ReprodPattern,
                                        InitRatioFem, FinalSex_Pmax, FinalSex_A50, FinalSex_A95, mat_A50, mat_A95,
                                        EstMatAtAge, Gear_sel_A50, Gear_sel_A95, EstGearSelAtAge, ret_Pmax, ret_A50, ret_A95,
                                        EstRetenAtAge, DiscMort, Steepness, Steepness_sd, SRrel_Type, NatMort, NatMort_sd,
                                        Current_F, Current_F_sd, PlotOpt, RefPointPlotOpt, FittedRes, nReps, MainLabel=NA,
                                        xaxis_lab=NA, yaxis_lab=NA, xmax=NA, xint=NA, ymax=NA, yint=NA)
        
        PlotPerRecruit_RiskMatrix_RelBiom(RefPointOpt=0, FittedRes, nReps=10)
        
        dev.off()
        
        dumi[[g]]=list(Res=Res,FittedRes=FittedRes)
      }
      rm(FMort)
    }
    YPR[[l]]=dumi
  }
}
toc()

#---Generate outputs --------------------------------------
#1. Export table of results
Store.per.recruit.results=vector('list',N.sp)
names(Store.per.recruit.results)=Keep.species
for(l in 1:N.sp)
{
  Neim=Keep.species[l]
  if(!is.null(Length.catch.curve[[l]]))
  {
    dd=vector('list',length(Length.catch.curve[[l]]))
    for(g in 1:length(dd))
    {
      Fs=unlist(Length.catch.curve[[l]][[g]]$ParamEst[1,])
      names(Fs)=c('F','F_Low','F_Upp')
      dd[[g]]=data.frame(value=c(Fs,
                                 unlist(YPR[[l]][[g]]$Res[c('BMSY_Targ','BMSY_Thresh','BMSY_Lim','F_tar','F_MSY','F_lim','Eq_FemRelSpBiom')]), 
                                 unlist(YPR[[l]][[g]]$FittedRes$'ResSummary_with_err')))%>%
        tibble::rownames_to_column(var='Parameter')%>%
        mutate(Year=names(Length.catch.curve[[l]])[g])%>%
        spread(Parameter,value)%>%
        mutate(Species=names(Length.catch.curve)[l])
    }
    Store.per.recruit.results[[l]]=do.call(rbind,dd)
  }
}
Out.table=do.call(rbind,Store.per.recruit.results)%>%
  relocate(Species)%>%
  rename(Period=Year)%>%
  adorn_rounding(digits = 3, rounding = "half up")
write.csv(Out.table,paste0(Rar.path,'/Table 6. Catch.curve_YPR.csv'),row.names = F)

#2. By species group
for(l in 1:length(Lista.sp.outputs))
{
  d=Out.table%>%
    filter(Species%in%Lista.sp.outputs[[l]])%>%
    mutate(Species=capitalize(Species))%>%
    dplyr::select(Species,Period,Fem_LowEquilSB,Fem_EquilSB,Fem_UppEquilSB,
                  BMSY_Lim,BMSY_Thresh,BMSY_Targ,F_Low,F,F_Upp,F_lim,F_MSY,F_tar)%>%
    mutate(Period=paste(sub("\\-.*", "", Period),'to',sub('.*-', '', Period)))
  NKoL=2
  if(length(unique(d$Species))<4) NKoL=1
  d%>%
    ggplot(aes(Fem_EquilSB,F))+
    geom_hline(aes(yintercept=F_tar),color=col.Target, size=.5)+
    geom_hline(aes(yintercept=F_MSY),color=col.Threshold, size=.5)+
    geom_hline(aes(yintercept=F_lim),color=col.Limit, size=.5)+
    geom_vline(aes(xintercept=BMSY_Targ),color=col.Target, size=.5)+
    geom_vline(aes(xintercept=BMSY_Thresh),color=col.Threshold, size=.5)+
    geom_vline(aes(xintercept=BMSY_Lim),color=col.Limit, size=.5)+
    geom_point(size=3)+
    geom_errorbarh(aes(xmin = Fem_LowEquilSB,xmax = Fem_UppEquilSB),linewidth=.1,position=position_dodge(.9))+
    geom_errorbar(aes(ymin=F_Low, ymax=F_Upp),width=.1,position=position_dodge(.9))+
    facet_wrap(~Species,ncol = NKoL)+ ylim(0, NA)+xlim(0,1)+
    theme_PA(strx.siz=12,axs.T.siz=14)+
    xlab(expression(B/~B[0]))+  #Relative female equilibrium breeding biomass
    ylab(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")))+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
    geom_text_repel(aes(label = Period),size = 4,col='steelblue')
  if(NKoL==2) DIMS=6
  if(NKoL==1) DIMS=3
  ggsave(paste(Rar.path,'/YPR_catch.curve_',names(Lista.sp.outputs)[l],'.tiff',sep=''),
         width = DIMS,height = 6,compression = "lzw")
}

# Simulate growth to eye-ball CVs -----------------------------------------------------------------------
if(do.eye.ball)
{
  #note: tweak CVSizeAtAge and FishMort to compare to observation dispersion in paper
  CV.scens=list(c(0.025,0.025),c(0.05,0.05),c(0.1,0.1),c(0.15,0.15))
  F.scens = c(0.01,0.05,0.1,0.2) # 0.05 low level of mortality
  TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
  LenInc = 5
  MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
  SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
  SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
  RetenParams = c(NA, NA) # L50, L95-L50 for retention
  DiscMort = 0 # proportion of fish that die due to natural mortality
  # single sex, von Bertalanffy
  GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
  RefnceAges = NA
  for(l in 1:N.sp)
  {
    print(paste("Eye ball CVs for --",names(Species.data)[l]))
    if(names(Species.data)[l]%in%c("dusky shark","gummy shark","sandbar shark","whiskery shark",
                                   "milk shark","smooth hammerhead","spinner shark"))
    {
      pdf(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                capitalize(List.sp[[l]]$Name),"/",AssessYr,
                "/Catch curve and per recruit/Effect of CVSizeAtAge.pdf",sep=''))
      for(f in 1:length(F.scens))
      {
        FishMort = F.scens[f]
        for(v in 1:length(CV.scens))
        {
          CVSizeAtAge = CV.scens[[v]]
          MaxAge = ceiling(List.sp[[l]]$Max.age.F[1])
          NatMort = List.sp[[l]]$Sens.test$SS3$Mmean[1]
          MaxLen = 10*round(List.sp[[l]]$TLmax/10)
          min.TL=with(List.sp[[l]],Lzero*a_FL.to.TL+b_FL.to.TL)
          lbnd = seq(0,MaxLen - LenInc, LenInc)
          ubnd = lbnd + LenInc
          midpt = lbnd + (LenInc/2)
          if(out.fleet=="TDGDLF") SS.flit='Southern.shark_1'
          if(out.fleet=="NSF") SS.flit='Northern.shark'
          if(out.fleet=="Other") SS.flit='Other'
          SelectivityVec_full=with(List.sp[[l]]$SS_selectivity%>%filter(Fleet==SS.flit),doubleNorm24.fn(midpt,a=P_1,
                                                                                                        b=P_2, c=P_3, d=P_4, e=P_5, f=P_6,use_e_999=FALSE, use_f_999=FALSE))
          SelectivityVec=SelectivityVec_full
          SelectivityVec[which(lbnd<min.TL)]=1e-6
          Linf = c(List.sp[[l]]$Growth.F$FL_inf*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL,
                   List.sp[[l]]$Growth.M$FL_inf*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL) #total length in cm for females and males   
          vbK = c(List.sp[[l]]$Growth.F$k,List.sp[[l]]$Growth.M$k)          # k for females and males
          GrowthParams = data.frame(Linf=Linf, vbK=vbK)
          
          
          #1.2 Simulated data
          SampleSize=5000
          set.seed(123)
          Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                                   SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
          
          #1.3 plot data
          PlotOpt=1 # 0=all plots, 1=retained lengths at age, 2=retained plus discarded lengths at age, 3=length frequency, 6
          # 7=female length frequency, 8=male length frequency, 9=female age frequency, 10=male age frequency,
          # 11=selectivity/retention, 12=F-at-age reten + disc, 13=F-at-age reten, 14=F-at-age disc
          PlotSimLenAndAgeFreqData(MaxAge, MaxLen, Res, PlotOpt)
          legend('topright',paste0(c('F= ','CV= '),c(FishMort,unique(CVSizeAtAge))),bty='n')
        }
      }
      dev.off()
    }
  }
}