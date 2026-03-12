library(AquaticLifeHistory)
library(BayesGrowth)
library(tidybayes)
Geraghty.dusky=read.csv(handl_OneDrive('Data/Age and growth/Geraghty et al 2013/Dusky.csv'))
Geraghty.sandbar=read.csv(handl_OneDrive('Data/Age and growth/Geraghty et al 2013/Sandbar.csv'))
plot.Bayes.growth=function(gc,DD)
{
  p=gc%>%
    ggplot(aes(Age, LAA))+
    geom_point(data = DD, aes(Age, Length), alpha = .3)+
    geom_lineribbon(aes( ymin = .lower, ymax = .upper, fill = factor(.width)), size = .8) +
    labs(y = "Fork length (cm)", x = "Age (yrs)")+
    scale_fill_brewer(palette="BuPu", direction=-1,name = "Credibility interval")+
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,max(DD$Age),1))+
    theme_bw()+
    theme(text = element_text(size = 14),
          legend.position = c(0.8,0.325),
          legend.background = element_rect(colour = "transparent"))+
    ylim(0,NA)
  return(p)
}
add.Frequentist.plot=function(p,Estimates)
{
  Estimates=Estimates%>%
    mutate(LAA=AVG)
  p <- p+
    geom_line(data=Estimates, aes(Age, AVG, col = Model),size = 1) + 
    geom_line(data=Estimates, aes(Age, low, col = Model),size = 1,alpha=0.4,linetype = "dotted") +
    geom_line(data=Estimates, aes(Age, upp, col = Model),size = 1,alpha=0.4,linetype = "dotted")
#    geom_ribbon(data=Estimates,aes(ymin = low, ymax = upp, fill = Model), alpha = 0.4)  
  return(p)

}
fit.growth.curve=function(SP,dat,LH,N.sims=1e4,K.max=1)
{
  ii=match(names(SP),names(dat))
  dd=dat[[ii]]$age_length%>%
    rename(Length=FL)%>%
    mutate(Sex=ifelse(Sex=='Male','M',ifelse(Sex=='Female','F',NA)))%>%
    dplyr::select(Age,Length,Sex)
  ii=match(names(SP),names(LH))
  Lh=LH[[ii]]
  birth_size=Lh$Lzero
  birth_size_se=1  
  max_size_se=5    
  wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(names(SP)),"/",AssessYr,sep='')
  if(names(SP)%in%c("dusky shark","sandbar shark"))
  {
    if(names(SP)=="dusky shark") add.geraghty=Geraghty.dusky
    if(names(SP)=="sandbar shark") add.geraghty=Sandbar.dusky
    add.geraghty=add.geraghty%>%
      filter(!Readability..amended.==1)%>%
      dplyr::select(AgeAgree,FL,Sex)%>%
      rename(Age=AgeAgree,
             Length=FL)
    
    rbind(add.geraghty%>%mutate(Data.set='Geraghty'),dd%>%mutate(Data.set='WA'))%>%
      mutate(Age=as.numeric(Age))%>%
      ggplot(aes(Age,Length,color=Data.set))+
      geom_point()+
      facet_wrap(~Sex,ncol = 1)
    ggsave(paste(wd,"/1_Inputs/Visualise data/refit_growth_compare_WA_Geraghty.tiff",sep=''), 
           width = 8,height = 8, dpi = 300, compression = "lzw")
    
    ns=table(dd$Sex)
    n.fem=ns[match('F',names(ns))]
    n.mal=ns[match('M',names(ns))]
    
    dd.gera.f=add.geraghty%>%filter(Sex=='F')
    dd.gera.f=dd.gera.f[sample(1:nrow(dd.gera.f),n.fem,replace = T),]
    dd.gera.m=add.geraghty%>%filter(Sex=='M')
    dd.gera.m=dd.gera.m[sample(1:nrow(dd.gera.m),n.mal,replace = T),]
    
    dd=rbind(dd,dd.gera.f,dd.gera.m)%>%mutate(Age=as.numeric(Age))
    
  }
  
  #Females
  SeXes=unique(dd$Sex)
  SeXes=subset(SeXes,!SeXes%in%c('','U'))
  for(x in 1:length(SeXes))
  {
    Nme=ifelse(SeXes[x]=='F','female',ifelse(SeXes[x]=='M','male',''))
    dd1=dd%>%filter(Sex==SeXes[x])
    if(SeXes[x]=='F')
    {
      max_size=max(c(max(dd1$Length),
                     with(Lh,(TLmax_obs-b_FL.to.TL)/a_FL.to.TL),
                     with(Lh,(TLmax-b_FL.to.TL)/a_FL.to.TL)))
    }
    if(SeXes[x]=='M')
    {
      max_size=max(c(max(dd1$Length),
                          with(Lh,(TLmax_obs_male-b_FL.to.TL)/a_FL.to.TL)))
    }
    
    #frequentist
    Growth=Estimate_Growth(dd1,Birth.Len=birth_size,plots = FALSE)
    #bayesian
    Growth.B <- Estimate_MCMC_Growth(data = dd1%>%dplyr::select(Length,Age), 
                                         Model = "VB" ,iter = 5e3,Linf = max_size,Linf.se = max_size_se,
                                         L0 = birth_size,L0.se = birth_size_se,k.max = K.max,sigma.max = 100)
    #export plots
    fit.summary=summary(Growth.B)
    fit.summary=fit.summary$summary%>%
      data.frame%>%
      dplyr::select(mean,sd)%>%
      rename( Parameter=mean,SE=sd)%>%
      mutate(Model='Bayesian.VonB')
    
    p=plot.Bayes.growth(gc= Calculate_MCMC_growth_curve(Growth.B, Model = "VB",max.age = max(dd1$Age), probs = c(.5,.95)),DD=dd1)
    p1=add.Frequentist.plot(p=p,Estimates=Growth$Estimates)
    
    p2=data.frame(Linf.Prior=rnorm(N.sims,max_size,max_size_se),
                  L0.Prior=rnorm(N.sims,birth_size,birth_size_se),
                  K.Prior=runif(N.sims,min=0,max=1),
                  Linf.Posterior=rnorm(N.sims,fit.summary$Parameter[1],fit.summary$SE[1]),
                  L0.Posterior=rnorm(N.sims,fit.summary$Parameter[3],fit.summary$SE[3]),
                  K.Posterior=rnorm(N.sims,fit.summary$Parameter[2],fit.summary$SE[2]))%>%
      pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")%>%
      mutate(type=str_extract(Parameter, "(?<=\\.).*"),
             Parameter = str_extract(Parameter , "^[^.]+"))%>%
      ggplot(aes(x = Value, fill = type)) +
      geom_density(alpha = 0.5) +
      facet_wrap(~Parameter, scales = "free") + 
      theme_PA()+
      theme(legend.position = 'top',
            legend.title = element_blank())
    
    ggarrange(plotlist = list(p1,p2),ncol=1,nrow=2)  
    ggsave(paste(wd,paste0("/1_Inputs/Visualise data/refit_growth_",Nme,".tiff"),sep=''), 
           width = 8,height = 8, dpi = 300, compression = "lzw")
    
    #Export table of par estimates
    out=rbind(Growth$VonB%>%data.frame%>%mutate(Model='VonB'),
              Growth$Logistic%>%data.frame %>%mutate(Model='Logistic'),
              Growth$Gompertz%>%data.frame%>%mutate(Model='Gompertz'))
    write.csv(rbind(out,fit.summary),paste(wd,"/1_Inputs/Visualise data/refit_growth_female.csv",sep=''))
    
    
  }
  
  

}
