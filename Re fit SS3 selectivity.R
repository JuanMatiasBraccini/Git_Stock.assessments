# Get SS3 selectivity parameters from refitting empirical selectivities ------------------------------------
fn.source1("SS_selectivity functions.R")
sumsq <- function( x, y) {sum((x-y)^2)}


#TDGDLF  
SS.double.normal.pars_empirical=vector('list',N.sp)
names(SS.double.normal.pars_empirical)=Keep.species
ObjFunc=function(par)
{
  Predicted=doubleNorm24.fn(TL,par[1],par[2],par[3],par[4],e=-999, f=-999, use_e_999=TRUE, use_f_999=TRUE)  
  return(sumsq(x=Obs,y=Predicted))
}
for(l in 1:N.sp)
{
  if(!is.null(Selectivity.at.totalength[[l]]) & any(grepl(paste(c('West','Zone'),collapse='|'),names(Species.data[[l]]))))
  {
    SP=Keep.species[l]
    print(paste('Extractiong SS_double normal sel pars for -----',SP))
    TLmax=List.sp[[l]]$TLmax
    if(Keep.species[l]%in%c('milk shark','sawsharks','spurdogs',
                            'gummy shark','whiskery shark')) TLmax=TLmax*1.5
    AA=Selectivity.at.totalength[[l]]%>%
      filter(TL<TLmax)
    TL=AA$TL-(TL.bins.cm/2)
    Obs=AA$Sel.combined
    par=c(90,-4.5,7,9)
    if(SP=='angel sharks') par=c(65,-4.5,5.5,7)
    if(SP=='copper shark') par=c(110,-2.5,6,8)
    if(SP=='dusky shark') par=c(90,-4.5,5.2,5.2)
    if(SP=='gummy shark') par=c(130,-4.5,6,7)
    if(grepl("hammerhead",SP)) par=c(100,-4,9,10.2)
    if(SP=='sawsharks') par=c(100,-4.5,6,7)
    if(SP=='whiskery shark') par=c(130,-4.5,6,7)
    
    Init.val=doubleNorm24.fn(TL,par[1],par[2],par[3],par[4],e=1e-5, f=1e-5, use_e_999=TRUE, use_f_999=TRUE)
    
    d.list=Species.data[[l]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[l]]))]
    if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
    if(any(grepl('Other',names(d.list)))) d.list=d.list[-grep('Other',names(d.list))]
    if(any(grepl('NSF',names(d.list)))) d.list=d.list[-grep('NSF',names(d.list))]
    if(any(grepl('Survey',names(d.list)))) d.list=d.list[-grep('Survey',names(d.list))]
    if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
    for(x in 1:length(d.list)) d.list[[x]]$Fleet=str_remove(str_remove(names(d.list)[x],'Size_composition_'),'.inch.raw')
    d.list=do.call(rbind,d.list)
    d.list=d.list%>%
      mutate(Region=Fleet,
             fleet=ifelse(grepl(paste(c('West','Zone'),collapse='|'),Fleet),'TDGDLF',Fleet),
             Fleet=ifelse(fleet=='TDGDLF' & year<=2005,'Southern.shark_1',
                          ifelse(fleet=='TDGDLF' & year>2005,'Southern.shark_2',
                                 ifelse(fleet=='NSF.LONGLINE','Northern.shark',
                                        fleet))),
             TL=FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)
    dd=d.list%>%
      filter(Fleet%in%c('Southern.shark_1','Southern.shark_2'))
    dd1=dd
    
    #display size and selectivity
    fn.fig(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[l]]$Name),"/",AssessYr,
                 "/1_Inputs/Visualise data/Estimated SS3 selectivity",sep=''),2000,2000)  
    plot(TL,Obs,main=SP,pch=19,ylab='Relative selectivity',xlab='TL (cm)')
    lines(TL,Init.val,col='red',lwd=1)
    if(SP=='copper shark')dd=dd%>%filter(Region=='Zone2.7')
    if(SP=='scalloped hammerhead') dd=dd%>%filter(Region%in%c('West.6.5','West.7'))
    if(nrow(dd)>2)
    {
      dd$bin=TL.bins.cm*floor(dd$TL/TL.bins.cm)
      y=table(dd$bin)
      y=y/max(y)
      points(as.numeric(names(y)),y,type='h',col="green")
      d <- density(dd$TL)
      #  lines(d$x,d$y/max(d$y),col="forestgreen",lwd=2)
    }
    nlmb <- nlminb(par, ObjFunc, gradient = NULL, hessian = TRUE)
    if(SP=='gummy shark')nlmb$par[1]=112
    if(SP=='spurdogs')nlmb$par[1]=56
    TL1=TL
    Predicted=with(nlmb,doubleNorm24.fn(TL1,par[1],par[2],par[3],par[4],e=1e-5, f=1e-5, use_e_999=TRUE, use_f_999=TRUE))
    lines(TL1,Predicted,col='black',lwd=3)
    Obs=d$y/max(d$y)
    TL=d$x
    if(SP=='great hammerhead') par=c(300,-1,10,10)
    
    #rescale peak
    nlmb_fitted_to_obs <- nlminb(par, ObjFunc, gradient = NULL, hessian = TRUE) 
    if(SP=='great hammerhead')nlmb_fitted_to_obs$par[3:4]=9 
    if(SP=='milk shark')nlmb_fitted_to_obs$par[2:4]=c(-9,4.5,4.5) 
    Predicted_fitted_to_obs=with(nlmb_fitted_to_obs,doubleNorm24.fn(TL1,par[1],par[2],par[3],par[4],e=1e-5, f=1e-5, use_e_999=TRUE, use_f_999=TRUE))
    lines(TL1,Predicted_fitted_to_obs,col='forestgreen',lwd=3)
    legend('topright',c('empirical sel','init values','predicted SS3',
                        'predicted SS3_rescaled','observed lengths'),
           pch=c(19,NA,NA,NA,NA),lty=c(NA,1,1,1,1),lwd=3,
           col=c('black','red','black','forestgreen',"green"),bty='n')
    dev.off()
    
    #display size comp by region
    pp=dd1%>%
      mutate(bin=TL.bins.cm*floor(TL/TL.bins.cm),
             Mesh=sub("^[^.]+.","",Region),
             Region=sub("\\..*", "", Region))%>%
      group_by(Region,bin,Mesh)%>%
      tally()%>%
      ggplot(aes(bin,n,color=Region))+
      geom_line()+
      facet_wrap(~Mesh,ncol=1)+xlab('TL (cm)')
    print(pp)  
    ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[l]]$Name),"/",AssessYr,
                 "/1_Inputs/Visualise data/Estimated SS3 selectivity_length comps region.tiff",sep=''),
           width = 8,height = 8, dpi = 300, compression = "lzw")
    dev.off()
    
    
    #store SS3 sel pars for use in assessment
    usedis=nlmb
    SS3.par.calculation='SS3 pars fitted to empirical selectivity'
    rescaled.species=c('great hammerhead','scalloped hammerhead','copper shark',
                       'grey nurse shark','milk shark','shortfin mako','sawsharks',
                       'spinner shark','tiger shark','wobbegongs')
    if(SP%in%rescaled.species)
    {
      usedis=nlmb_fitted_to_obs
      SS3.par.calculation='SS3 pars fitted to empirical selectivity rescaled to length comp'
    }
    SS.double.normal.pars_empirical[[l]]=with(usedis,data.frame(Species=SP,
                                                                p1=par[1],
                                                                p2=par[2],
                                                                p3=par[3],
                                                                p4=par[4],
                                                                p5=-999,
                                                                p6=-999,
                                                                Source=unique(AA$type),
                                                                SS3.par.calculation=SS3.par.calculation))
    rm(AA,nlmb,nlmb_fitted_to_obs,dd,dd1,usedis)
  }
}
write.csv(do.call(rbind,SS.double.normal.pars_empirical),
          handl_OneDrive("Analyses/Population dynamics/Refit SS3 selectivity/fitted to empirical selectivity_TDGDLF.csv"),
          row.names = F)

#NSF
SS.logistic.pars_empirical=vector('list',N.sp)  
names(SS.logistic.pars_empirical)=Keep.species
ObjFunc=function(par)
{
  Predicted=logistic1.fn(midpt, par[1],par[2])
  return(sumsq(x=Obs,y=Predicted))
}
for(l in 1:N.sp)
{
  if(!is.null(Selectivity.at.totalength[[l]]))
  {
    print(paste('Extractiong SS_logistic sel pars for -----',Keep.species[l]))
    attach(List.sp[[l]])
    if(MN.SZE=="size.at.birth") Min.size.bin=10*round((Lzero*a_FL.to.TL+b_FL.to.TL)/10)
    if(MN.SZE==0) Min.size.bin=0
    MaxLen= 10*ceiling(TLmax/10)
    lbnd = seq(Min.size.bin,MaxLen - TL.bins.cm, TL.bins.cm)
    ubnd = lbnd + TL.bins.cm
    midpt = lbnd + (TL.bins.cm/2)
    Obs=round(1/(1+(exp(-log(19)*((midpt-TL.50.mat)/(TL.95.mat-TL.50.mat))))),3)
    par2=20
    par=c(TL.50.mat,par2)
    detach(List.sp[[l]])
    
    nlmb <- nlminb(par, ObjFunc, gradient = NULL, hessian = TRUE)
    Predicted=logistic1.fn(midpt,nlmb$par[1],nlmb$par[2])
    Init.val=logistic1.fn(midpt,par[1],par[2])
    
    fn.fig(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[l]]$Name),"/",AssessYr,
                 "/1_Inputs/Visualise data/Estimated SS3 selectivity_NSF",sep=''),2000,2000)  
    plot(midpt,Obs,main=Keep.species[l],pch=19,ylab='Relative selectivity',xlab='TL (cm)')
    lines(midpt,Init.val,col='grey90',lwd=1)
    lines(midpt,Predicted,col='red',lwd=2)
    legend('topright',c('observed','init values','predicted'),pch=c(19,NA,NA),lty=c(NA,1,1),lwd=3,
           col=c('black','grey90','red'),bty='n')
    dev.off()
    
    SS.logistic.pars_empirical[[l]]=data.frame(Species=Keep.species[l],
                                               p1=nlmb$par[1],
                                               p2=nlmb$par[2],
                                               Source='Set at maturity ogive')
    rm(nlmb)
  }
}
write.csv(do.call(rbind,SS.logistic.pars_empirical),
          handl_OneDrive("Analyses/Population dynamics/Refit SS3 selectivity/fitted to empirical selectivity_NSF_logistic.csv"),
          row.names = F)

#Other   
SS.double.normal.pars_empirical_Other=vector('list',N.sp)
names(SS.double.normal.pars_empirical_Other)=Keep.species
ObjFunc=function(par)
{
  Predicted=doubleNorm24.fn(TL,par[1],par[2],par[3],par[4],e=-999, f=-999, use_e_999=TRUE, use_f_999=TRUE)  
  return(sumsq(x=Obs,y=Predicted))
}
for(l in 1:N.sp)
{
  if(!is.null(Selectivity.at.totalength[[l]]) & any(grepl(paste(c('Other'),collapse='|'),names(Species.data[[l]]))))
  {
    SP=Keep.species[l]
    print(paste('Extractiong SS_double normal sel pars Other for -----',SP))
    TLmax=List.sp[[l]]$TLmax
    if(Keep.species[l]%in%c('milk shark','sawsharks','spurdogs',
                            'gummy shark','whiskery shark')) TLmax=TLmax*1.5
    AA=Selectivity.at.totalength[[l]]%>%
      filter(TL<TLmax)
    TL=AA$TL-(TL.bins.cm/2)
    Obs=AA$Sel.combined
    par=c(90,-4.5,7,9)
    if(SP=='angel sharks') par=c(65,-4.5,5.5,7)
    if(SP=='copper shark') par=c(110,-2.5,6,8)
    
    
    fn.fig(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[l]]$Name),"/",AssessYr,
                 "/1_Inputs/Visualise data/Estimated SS3 selectivity_Other",sep=''),2000,2000)  
    
    Init.val=doubleNorm24.fn(TL,par[1],par[2],par[3],par[4],e=1e-5, f=1e-5, use_e_999=TRUE, use_f_999=TRUE)
    
    plot(TL,Obs,main=SP,pch=19,ylab='Relative selectivity',xlab='TL (cm)')
    lines(TL,Init.val,col='red',lwd=1)
    
    #size comp
    d.list=Species.data[[l]][grep("Size_composition_Other",names(Species.data[[l]]))]
    if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
    for(x in 1:length(d.list)) d.list[[x]]$Fleet=str_remove(str_remove(names(d.list)[x],'Size_composition_'),'.inch.raw')
    d.list=do.call(rbind,d.list)
    dd=d.list%>%mutate(fleet=ifelse(grepl(paste(c('West','Zone'),collapse='|'),Fleet),'TDGDLF',Fleet),
                       Fleet=ifelse(fleet=='TDGDLF' & year<=2005,'Southern.shark_1',
                                    ifelse(fleet=='TDGDLF' & year>2005,'Southern.shark_2',
                                           ifelse(fleet=='NSF.LONGLINE','Northern.shark',
                                                  fleet))),
                       TL=FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)
    
    if(nrow(dd)>2)
    {
      dd$bin=TL.bins.cm*floor(dd$TL/TL.bins.cm)
      y=table(dd$bin)
      y=y/max(y)
      points(as.numeric(names(y)),y,type='h',col="green")
      d <- density(dd$TL)
      #  lines(d$x,d$y/max(d$y),col="forestgreen",lwd=2)
    }
    
    nlmb <- nlminb(par, ObjFunc, gradient = NULL, hessian = TRUE)
    TL1=TL
    Predicted=with(nlmb,doubleNorm24.fn(TL1,par[1],par[2],par[3],par[4],e=1e-5, f=1e-5, use_e_999=TRUE, use_f_999=TRUE))
    lines(TL1,Predicted,col='black',lwd=3)
    
    Obs=d$y/max(d$y)
    TL=d$x
    nlmb_fitted_to_obs <- nlminb(par, ObjFunc, gradient = NULL, hessian = TRUE) #rescale peak
    if(SP=='copper shark') nlmb_fitted_to_obs$par[c(2,4)]=c(-25,8.5)
    Predicted_fitted_to_obs=with(nlmb_fitted_to_obs,doubleNorm24.fn(TL1,par[1],par[2],par[3],par[4],e=1e-5, f=1e-5, use_e_999=TRUE, use_f_999=TRUE))
    lines(TL1,Predicted_fitted_to_obs,col='forestgreen',lwd=3)
    
    
    legend('topright',c('empirical sel','init values','predicted SS3',
                        'predicted SS3_rescaled','length comps'),
           pch=c(19,NA,NA,NA,NA),lty=c(NA,1,1,1,1),lwd=3,
           col=c('black','red','black','forestgreen',"green"),bty='n')
    dev.off()
    
    usedis=nlmb
    SS3.par.calculation='SS3 pars fitted to empirical selectivity'
    if(grepl(paste(c('copper','angel'),collapse='|'),SP))
    {
      usedis=nlmb_fitted_to_obs
      SS3.par.calculation='SS3 pars fitted to empirical selectivity rescaled to length comp'
    }
    
    SS.double.normal.pars_empirical_Other[[l]]=with(usedis,data.frame(Species=SP,
                                                                      p1=par[1],
                                                                      p2=par[2],
                                                                      p3=par[3],
                                                                      p4=par[4],
                                                                      p5=-999,
                                                                      p6=-999,
                                                                      Source=unique(AA$type),
                                                                      SS3.par.calculation=SS3.par.calculation))
    rm(AA,nlmb,nlmb_fitted_to_obs)
  }
  if(Keep.species[l]=='sawsharks')
  {
    SP=Keep.species[l]
    print(paste('Extractiong SS_double normal sel pars Other for -----',SP))
    
    fn.fig(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[l]]$Name),"/",AssessYr,
                 "/1_Inputs/Visualise data/Estimated SS3 selectivity_Other",sep=''),2000,2000)  
    
    dummi=data.frame(TL=sort(runif(1e3,103*.6,129*1.4)),
                     Obs=rep(1,1e3))%>%
      mutate(Obs=case_when(TL >= 103*.9 & TL < 103~ .9,
                           TL >= 103*.8 & TL < 103*.9~ .75,
                           TL >= 103*.7 & TL < 103*.8~ .5,
                           TL < 103*.7 ~ .01,
                           TL > 129*1.3 ~ .01,
                           TL >= 129*1.2 & TL <= 129*1.3~ .5,
                           TL >= 129*1.1 & TL < 129*1.2 ~ .75,
                           TL >= 129 & TL < 129*1.1~ .9,
                           TRUE~Obs))
    
    TL=dummi$TL
    Obs=dummi$Obs
    par=c(115,-2,8,8)
    Init.val=doubleNorm24.fn(TL,par[1],par[2],par[3],par[4],e=1e-5, f=1e-5, use_e_999=TRUE, use_f_999=TRUE)
    
    plot(TL,Obs,main=SP,pch=19,ylab='Relative selectivity',xlab='TL (cm)')
    lines(TL,Init.val,col='red',lwd=1)
    nlmb <- nlminb(par, ObjFunc, gradient = NULL, hessian = TRUE)
    TL1=TL
    Predicted=with(nlmb,doubleNorm24.fn(TL1,par[1],par[2],par[3],par[4],e=1e-5, f=1e-5, use_e_999=TRUE, use_f_999=TRUE))
    lines(TL1,Predicted,col='black',lwd=3)
    dev.off()
    
    SS.double.normal.pars_empirical_Other[[l]]=with(nlmb,data.frame(Species=SP,
                                                                    p1=par[1],
                                                                    p2=par[2],
                                                                    p3=par[3],
                                                                    p4=par[4],
                                                                    p5=-999,
                                                                    p6=-999,
                                                                    Source="Koopman",
                                                                    SS3.par.calculation='SS3 pars fitted to size range'))
    
  }
}
write.csv(do.call(rbind,SS.double.normal.pars_empirical_Other),
          handl_OneDrive("Analyses/Population dynamics/Refit SS3 selectivity/fitted to empirical selectivity_Other.csv"),
          row.names = F)



# Fit SS3 selectivity to length composition data ------------------------------------

#NSF - logistic 
SS.logistic_length.comps=vector('list',N.sp)  
names(SS.logistic_length.comps)=Keep.species
ObjFunc=function(par)
{
  Predicted=logistic1.fn(len=TL, par[1],par[2])
  return(sumsq(x=Obs,y=Predicted))
}
pdf(handl_OneDrive("Analyses/Population dynamics/Refit SS3 selectivity/fitted to length comps_NSF_logistic.pdf"))
for(l in 1:N.sp)
{
  if(any(grepl(paste(c('NSF.LONGLINE','Survey'),collapse='|'),names(Species.data[[l]]))))
  {
    SP=Keep.species[l]
    print(paste('Extractiong logistic sel pars for -----',SP))

    
    d.list=Species.data[[l]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[l]]))]
    if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
    if(any(grepl('Other',names(d.list)))) d.list=d.list[-grep('Other',names(d.list))]
    if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
    for(x in 1:length(d.list)) d.list[[x]]$Fleet=str_remove(str_remove(names(d.list)[x],'Size_composition_'),'.inch.raw')
    d.list=do.call(rbind,d.list)
    d.list=d.list%>%
      mutate(TL=FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)
    dd=d.list%>%
      filter(Fleet%in%c('NSF.LONGLINE','Survey'))

    
    #display size and selectivity
    if(nrow(dd)>10)
    {
      dd$bin=TL.bins.cm*floor(dd$TL/TL.bins.cm)
      y=table(dd$bin)
      Obs=y/max(y)
      TL=as.numeric(names(y))
      id=match(max(Obs),Obs)
      Obs1=Obs
      Obs[(id+1):length(Obs)]=1
      par=c(List.sp[[l]]$TL.50.mat,20)
      
      nlmb <- nlminb(par, ObjFunc, gradient = NULL, hessian = TRUE)
      Predicted=logistic1.fn(TL,nlmb$par[1],nlmb$par[2])
      
      
      #fn.fig(handl_OneDrive(paste("Analyses/Population dynamics/Refit SS3 selectivity/fitted to length comps_NSF_logistic",SP,sep='_')),2000,2000)  
      plot(TL,Obs,type='h',col="green",main=SP)
      points(TL,Obs1,type='h',col="red")
      lines(TL,Predicted,col='black',lwd=3)
      legend('topleft',c('obs','obs dummy to fit logis'),lty=1,col=c('red','green'),bty='n')
      
      
      #store SS3 sel pars for use in assessment
      usedis=nlmb
      SS.logistic_length.comps[[l]]=with(usedis,data.frame(Species=SP,
                                                                  p1=par[1],
                                                                  p2=par[2],
                                                                  Source='length composition'))
      rm(nlmb,dd,usedis)
    }
  }
}
dev.off()
write.csv(do.call(rbind,SS.logistic_length.comps),
          handl_OneDrive("Analyses/Population dynamics/Refit SS3 selectivity/fitted to length comps_NSF_logistic.csv"),
          row.names = F)


#NSF - double normal 
SS.double.normal_length.comps=vector('list',N.sp)
names(SS.double.normal_length.comps)=Keep.species
ObjFunc=function(par)
{
  Predicted=doubleNorm24.fn(TL,par[1],par[2],par[3],par[4],e=-999, f=-999, use_e_999=TRUE, use_f_999=TRUE)  
  return(sumsq(x=Obs,y=Predicted))
}
pdf(handl_OneDrive("Analyses/Population dynamics/Refit SS3 selectivity/fitted to length comps_NSF_double_normal.pdf"))
for(l in 1:N.sp)
{
  if(any(grepl(paste(c('NSF.LONGLINE','Survey'),collapse='|'),names(Species.data[[l]]))))
  {
    SP=Keep.species[l]
    print(paste('Extracting double normal sel pars for -----',SP))
    
    
    d.list=Species.data[[l]][grep(paste(SS3_fleet.size.comp.used,collapse="|"),
                                  names(Species.data[[l]]))]
    if(any(grepl('Observations',names(d.list)))) d.list=d.list[-grep('Observations',names(d.list))]
    if(any(grepl('Other',names(d.list)))) d.list=d.list[-grep('Other',names(d.list))]
    if(sum(grepl('Table',names(d.list)))>0) d.list=d.list[-grep('Table',names(d.list))]
    for(x in 1:length(d.list)) d.list[[x]]$Fleet=str_remove(str_remove(names(d.list)[x],'Size_composition_'),'.inch.raw')
    d.list=do.call(rbind,d.list)
    d.list=d.list%>%
      mutate(TL=FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)
    dd=d.list%>%
      filter(Fleet%in%c('NSF.LONGLINE','Survey'))
    
    
    #display size and selectivity
    if(nrow(dd)>20)
    {
      dd$bin=TL.bins.cm*floor(dd$TL/TL.bins.cm)
      y=table(dd$bin)
      if(SP=="great hammerhead") y[17]=7
      Obs=y/max(y)
      TL=as.numeric(names(y))
      id=match(max(Obs),Obs)
      par=c(List.sp[[l]]$TL.50.mat,-8,8,4)
      if(SP=="great hammerhead") par=c(260,-20,7.5,6.5)
      if(SP=="pigeye shark") par=c(215,-6,6.5,5.5)
      if(SP=="sandbar shark") par=c(155,-6,6,5.25)
      if(SP=="scalloped hammerhead") par=c(180,-20,9,8)
      if(SP=="spinner shark") par=c(180,-2,10,10)
      if(SP=="tiger shark") par=c(110,-10,5,10)
      
      nlmb <- nlminb(par, ObjFunc, gradient = NULL, hessian = TRUE)
      Predicted=with(nlmb,doubleNorm24.fn(TL,par[1],par[2],par[3],par[4],e=-999, f=-999, use_e_999=TRUE, use_f_999=TRUE))
      
      plot(TL,Obs,type='h',col="green",main=SP,ylim=c(0,1.1))
      lines(TL,Predicted,col='green',lwd=3)
      lines(TL,doubleNorm24.fn(TL,par[1],par[2],par[3],par[4],e=-999, f=-999, use_e_999=TRUE, use_f_999=TRUE),col='red',lwd=3)
      legend('topleft',c(paste(c('obs (',paste(paste0('p',1:4),round(nlmb$par,2),sep='='),')'),collapse=''),
                         paste(c('init par (',paste(paste0('p',1:4),round(par,2),sep='='),')'),collapse='')),
             lty=1,col=c('green','red'),bty='n',cex=0.8)
      
      
      #store SS3 sel pars for use in assessment
      usedis=nlmb
      SS.double.normal_length.comps[[l]]=with(usedis,data.frame(Species=SP,
                                                                 p1=par[1],
                                                                 p2=par[2],
                                                                 p3=par[3],
                                                                 p4=par[4],
                                                                 Source='length composition'))
      rm(nlmb,dd,usedis)
    }
  }
}
dev.off()
write.csv(do.call(rbind,SS.double.normal_length.comps),
          handl_OneDrive("Analyses/Population dynamics/Refit SS3 selectivity/fitted to length comps_NSF_double_normal.csv"),
          row.names = F)
