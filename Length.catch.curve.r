#to install latest version of L3Assess
# library(devtools)
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)
# install.packages("C:/~/L3Assess_0.1.0.tar.gz", source = TRUE, repos=NULL) #Install directly from .gz file

library(L3Assess) 
#note: identify relevant fishery per species and a good sampling period (i.e. fishery with highest catch and length
#       composition from exploited period to be able to determine F) and combine 2/3 years.
if(!exists('doubleNorm24.fn')) fn.source1("SS_selectivity functions.R")

Length.catch.curve=vector('list',N.sp)
names(Length.catch.curve)=Keep.species
tic()
for(l in 1: N.sp)
{
  #identify fishery with most catch 
  Katch1=KtCh%>%
        filter(Name==names(Length.catch.curve)[l])%>%
        mutate(Fishery=case_when(FishCubeCode%in%c('Historic',Southern.shark.fisheries)~"TDGDLF",
                                 FishCubeCode%in%Northern.shark.fisheries~"NSF",
                                 .default="Other"))
  Katch=Katch1%>%
        group_by(Fishery)%>%
        summarise(Tonnes=sum(LIVEWT.c))
  this.size.comp=Katch%>%filter(Tonnes==max(Tonnes))%>%pull(Fishery)
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
  if(length(iid)>0)
  {
    dummy=do.call(rbind,iid)%>%
      filter(year<=as.numeric(substr(Last.yr.ktch,1,4)))
    
    Min.sample.size=Min.annual.obs.ktch
    if(names(Length.catch.curve)[l]%in%names(Indicator.species)) Min.sample.size=300
    
    N.min=dummy%>%
      group_by(FINYEAR)%>%
      tally()%>%
      filter(n>=Min.sample.size)%>%
      mutate(Keep=FINYEAR)
    if(nrow(N.min)>0)
    {
      print(paste("Length-based catch curve with input selectivity for --",names(Species.data)[l],'--',out.fleet))
      
      #keep years with min sample size and calculate total length
      dummy=dummy%>%
        mutate(Keep=FINYEAR,
               year.f=as.numeric(substr(FINYEAR,1,4)))%>%
        filter(Keep%in%N.min$Keep)%>%
        mutate(TL=FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)     
      
      #Model inputs
      MaxAge = ceiling(mean(List.sp[[l]]$Max.age.F))
      TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
      NatMort = List.sp[[l]]$Sens.test$SS3$Mmean[1]
      MaxLen = 10*round(List.sp[[l]]$TLmax/10)
      LenInc = TL.bins.cm
      MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
      SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
      lbnd = seq(0,MaxLen - LenInc, LenInc)
      ubnd = lbnd + LenInc
      midpt = lbnd + (LenInc/2)
      if(out.fleet=="TDGDLF") SS.flit='Southern.shark_1'
      if(out.fleet=="NSF") SS.flit='Northern.shark'
      if(out.fleet=="Other") SS.flit='Other'
      SelectivityVec=with(List.sp[[l]]$SS_selectivity%>%filter(Fleet==SS.flit),doubleNorm24.fn(midpt,a=P_1,
                                      b=P_2, c=P_3, d=P_4, e=P_5, f=P_6,use_e_999=FALSE, use_f_999=FALSE))
      SelParams = c(300, 50) # L50, L95-L50 for gear selectivity  NOT USED
      RetenParams = c(NA, NA) # L50, L95-L50 for retention
      DiscMort = 0 # proportion of fish that die due to natural mortality
      DistnType = 1 # 1 = Multinomial, 2 = Dirichlet multinomial
      GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute 
      Linf = c(List.sp[[l]]$Growth.F$FL_inf*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL,
               List.sp[[l]]$Growth.M$FL_inf*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL) #total length in cm for females and males   
      vbK = c(List.sp[[l]]$Growth.F$k,List.sp[[l]]$Growth.M$k)          # k for females and males
      CVSizeAtAge = c(0.1,0.1)
      GrowthParams = data.frame(Linf=Linf, vbK=vbK)
      RefnceAges = NA
      InitL50 = List.sp[[l]]$TL.50.mat
      InitDelta = 50
      
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
        geom_line(data=data.frame(midpt=midpt,SelectivityVec=SelectivityVec),
                  aes(midpt,SelectivityVec),color='red',linewidth=1.25)+
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
              axis.title=element_text(size=16))
      
      ggarrange(p1, p2, ncol=1,nrow=2, heights = c(1,0.5))
      
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(List.sp[[l]]$Name),"/",AssessYr,
                   "/1_Inputs/Visualise data/Size.comp_Catch curve_",out.fleet,".tiff",sep=''),
             width = 8,height = 8,compression = "lzw")

      #Get observed length composition for selected period 
      Catch.curve.years=List.sp[[l]]$Catch.curve.years 
      dummy=dummy%>%filter(FINYEAR%in%Catch.curve.years) 
      ObsCatchFreqAtLen=vector('list',1)
      add.dummy=data.frame(bin=midpt)
      x=dummy%>%
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
      
      plot(midpt,ObsCatchFreqAtLen/max(ObsCatchFreqAtLen));lines(midpt,SelectivityVec)
      
      #Get F for selected years
      ObsRetCatchFreqAtLen=ObsCatchFreqAtLen[1,]
      ObsDiscCatchFreqAtLen = NA # (or set to Res$ObsDiscCatchFreqAtLen)
      PropReleased = NA # proportion of fish released, vector including mean and sd (option probably now obselete)
      InitFishMort = 0.4 # specify starting parameters
      InitFishMort_logit = log(InitFishMort/(1-InitFishMort)) # logit transform
      params = c(InitFishMort_logit)
      FittedRes=GetLengthBasedCatchCurveResults(params, DistnType, GrowthCurveType, GrowthParams, RefnceAges, MLL, 
                                                SelectivityType, ObsRetCatchFreqAtLen,lbnd, ubnd, midpt, 
                                                SelectivityVec, PropReleased, ObsDiscCatchFreqAtLen, 
                                                DiscMort, CVSizeAtAge, MaxAge, NatMort, TimeStep)
      FittedRes$Finyear=Catch.curve.years
      Length.catch.curve[[l]]=FittedRes
      rm(dummy,iid,this.size.comp,ObsCatchFreqAtLen,SelectivityVec,FittedRes,params)
    }
  }
}
toc()