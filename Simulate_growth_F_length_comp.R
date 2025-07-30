# For Matias
# suggest updating both packages (as have made many changes recently). 
# install package from github
# library(devtools)
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)
# devtools::install_github("SAlexHesp/WAFishBiologyRPackage", build_vignettes=TRUE, force=TRUE)


# library(WAFishBiology)
library(L3Assess)

# Simulate data
Alex_sim.data=function(SampleSize=5000, # sample size for retained catches (and same number for released fish, if an MLL is specified)
                       MaxAge = 31,
                       TimeStep = 1, # model timestep (e.g. 1 = annual, 1/12 = monthly)
                       FishMort = 0.02,
                       MaxLen = 2500,
                       LenInc = 20,
                       MLL=NA, # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
                       SelectivityType=2, # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
                       SelectivityAtLen = NA, # selectivity vector
                       SelParams = c(1500, 200), # L50, L95-L50 for gear selectivity
                       RetenParams = c(NA, NA), # L50, L95-L50 for retention
                       DiscMort = 0, # proportion of fish that die due to natural mortality
                       GrowthCurveType = 1, # 1 = von Bertalanffy, 2 = Schnute
                       Linf = c(1800,1650),
                       vbK = c(0.42,0.45),
                       CVSizeAtAge = c(0.04,0.04))
{
  set.seed(123)
  NatMort = 4.22/MaxAge
  
  GrowthParams = data.frame(Linf=Linf, vbK=vbK)
  RefnceAges = NA
  Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                                 SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
  
  # scatterplot data
  p.growth=data.frame(x=Res$ObsDecAgeRetCatch_Fem,
                      y=Res$ObsRandLenRetCatch_Fem)%>%
              ggplot(aes(x,y))+
              geom_point(color='pink')+
              xlim(0,MaxAge)+ylim(0,MaxLen)+
              theme_PA(Ttl.siz=10)+
              geom_point(data=data.frame(x=Res$ObsDecAgeRetCatch_Mal+.2,
                                         y=Res$ObsRandLenRetCatch_Mal),
                         aes(x,y),color='blue')+
              geom_hline(yintercept =SelParams[1])+
              geom_hline(yintercept =SelParams[1]+SelParams[2],linetype="dotted")+
              ylab('Length')+xlab('Age')+
    ggtitle(paste(paste0('F=',FishMort),
                  paste0('MaxLen=',MaxLen),
                  paste0('Linf [',paste(Linf,collapse=","),']'),
                  paste0('vbK [',paste(vbK,collapse=","),']'),
                  paste0('CV [',paste(CVSizeAtAge,collapse=","),']'),
                  paste0('SelPars [',paste(SelParams,collapse=","),']'),
                  sep='; '))
    

  # length frequency dist
  p.len.fre=data.frame(x=Res$midpt,
                       y=Res$ObsRetCatchFreqAtLen_Fem)%>%
                      ggplot(aes(x,y))+
                      geom_line(color='pink')+
                      xlim(0,MaxLen)+ylim(0,max(c(Res$ObsRetCatchFreqAtLen_Mal,Res$ObsRetCatchFreqAtLen_Fem)))+
                      theme_PA()+
                      geom_line(data=data.frame(x=Res$midpt,
                                                 y=Res$ObsRetCatchFreqAtLen_Mal),
                                 aes(x,y),color='blue')+
    geom_vline(xintercept = Linf[1],color="pink",linewidth=2, linetype="dashed")+
    geom_vline(xintercept = Linf[2],color="blue",linewidth=2, linetype="dashed")+
    ylab('Frequency')+xlab('Length class (mm)')
  

  #Selectivity
  p.sel=data.frame(x=Res$midpt,
                   y=Res$ModelDiag$SelAtLength)%>%
                ggplot(aes(x,y))+
                geom_line()+theme_PA()+
    ylab('Selectivity')+xlab('Length class (mm)')
  
  p=ggarrange(p.growth,p.sel,p.len.fre,ncol=1)
  
  return(p)
}


# # Example with selectivity specified as a vector
# # Simulate data
# SampleSize=1000
# set.seed(123)
# MaxAge = 30
# TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
# NatMort = 4.22/MaxAge
# FishMort = 0.2
# MaxLen = 1200
# LenInc = 20
# MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
# SelectivityType=1 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
# lbnd = seq(0,MaxLen - LenInc, LenInc)
# midpt = lbnd + (LenInc/2)
# # # single sex selectivity input
# # SelectivityAtLen = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
# # two sex selectivity input
# FemSelAtLen = 1 / (1 + exp(-log(19)*(midpt-400)/(500-400)))
# MalSelAtLen = 1 / (1 + exp(-log(19)*(midpt-450)/(550-450)))
# SelectivityAtLen = t(data.frame(FemSelAtLen=FemSelAtLen,MalSelAtLen=MalSelAtLen))
# colnames(SelectivityAtLen) = midpt
# SelParams = c(NA, NA) # L50, L95-L50 for gear selectivity
# RetenParams = c(NA, NA) # L50, L95-L50 for retention
# DiscMort = 0 # proportion of fish that die due to natural mortality
# # 2 sexes, von Bertalanffy
# GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
# Linf = c(700,850)
# vbK = c(0.3,0.2)
# CVSizeAtAge = c(0.06,0.06)
# GrowthParams = data.frame(Linf=Linf, vbK=vbK)
# RefnceAges = NA
# SimRes=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#                                   SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
