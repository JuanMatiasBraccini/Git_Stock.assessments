# SCRIPT FOR ORGANISING INPUT PARAMETERS AND RELATIONSHIPS
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

fn.input.pars=function(SP,add.growth.cv,add.Mrt.age)
{
  #Read in parameters from centralised parameter data base
  Pars=read.csv(handl_OneDrive("Data/Life history parameters/Shark input parameters.csv"),stringsAsFactors=F)
  Pars=subset(Pars,Population%in%c("WA","Australia","southern Australia"))
  Pars$SP=with(Pars,ifelse(Species=='Sandbar',"TK",ifelse(Species=='Dusky',"BW",
          ifelse(Species=='Gummy',"GM",ifelse(Species=='Whiskery',"WH",NA)))))
  
  Pars=subset(Pars,SP==SP)
  Pars.WA=subset(Pars,Population=="WA")
  Pars.southern.Oz=subset(Pars,Population=="southern Australia")
  Pars.Oz=subset(Pars,Population=="Australia")
  
  hndl=handl_OneDrive("Data/Population dynamics/Prop.males.in.catch/prop.males.")
  Prop.males.in.ktch=read.csv(paste(hndl,SP,".csv",sep=""))
  Prop.males.in.ktch.all=read.csv(paste(hndl,"All.",SP,".csv",sep=""))
  names(Prop.males.in.ktch.all)="All"
  Prop.males.in.ktch=cbind(Prop.males.in.ktch,Prop.males.in.ktch.all)
  
  
  #functions
  fn.get.par=function(sp,pars,sex)
  {
    a=subset(Pars.WA,SP==sp & Parameter%in%pars & Sex%in%sex)
    if(nrow(a)==0) a=subset(Pars.southern.Oz,SP==sp & Parameter%in%pars & Sex%in%sex)
    if(nrow(a)==0) a=subset(Pars.Oz,SP==sp & Parameter%in%pars & Sex%in%sex)
    return(list(mean.par=a$Mean,sd.par=a$SD,min.par=a$Min,max.par=a$Max,
                Low_95_CI.par=a$Low_95_CI,Up_95_CI.par=a$Up_95_CI,Source=a$Source))
  }
  fn.names=function(d,nam)
  {
    for(s in 1:length(d)) 
    {
      if(length(d[[s]])>0)
      {
        names(d[[s]])=nam
        d[[s]]=data.frame(t(d[[s]]))
      }
    }
    
    return(d)
  }
  
  
  #---PARAMETERS----
  N=length(SP)
  PARAMETERS=unique(Pars$Parameter)
  
  #1. Max age
  Max.Age.F=vector('list',length(SP))
  names(Max.Age.F)=SP
  Max.Age.F.source=Max.Age.M=Max.Age.M.source=Max.Age.F
  what='min.par'
  for (i in 1:N) 
  {
    dummy=fn.get.par(SP[i],"maximum age",'F')
    id=which(names(dummy)==what)
    Max.Age.F[[i]]=dummy[id]
    Max.Age.F.source[[i]]=dummy$Source
    
    dummy=fn.get.par(SP[i],"maximum age",'M')
    Max.Age.M[[i]]=dummy[id]
    Max.Age.M.source[[i]]=dummy$Source
  }
  
  
  #2. Natural mortality (yr-1)
  #2.1 Age independent
  M=M.source=Max.Age.F
  what='mean.par'
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],"Natural mortality",'Combined')
    id=which(names(dummy)==what)
    M[[i]]=dummy[id]
    M.source[[i]]=dummy$Source
  }
  
  #2.2. Age dependent
  setwd(handl_OneDrive("Data/Population dynamics/Parameter inputs for models"))
  if(add.Mrt.age=="YES")
  {
    #Median values from reference point paper. Calculated in "1.MSY proxy reference points.R"
    M.w.age=read.csv("whiskery.M_at_age.csv")
    M.g.age=read.csv("gummy.M_at_age.csv")
    M.d.age=read.csv("dusky.M_at_age.csv")
    M.s.age=read.csv("sandbar.M_at_age.csv")
    M.age=list(WH=M.w.age$Median.M,GM=M.g.age$Median.M,
               BW=M.d.age$Median.M,TK=M.s.age$Median.M)
    id=which(names(M.age)==SP)
    M.age=M.age[id]
    M.age.source=M.age
    for (i in 1:N) M.age.source[[i]]="Braccini et al 2015 (ref points)"
    
  }
  
  #3. Allometry
  #3.1. PL to TL (Heald 1987; TL=(PL-b_PL)/a_PL)
  PL.to.TL=PL.to.TL.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],c("a_PL","b_PL"),'Combined')
    PL.to.TL[[i]]=dummy$mean.par
    PL.to.TL.source[[i]]=dummy$Source
  }
  PL.to.TL=fn.names(PL.to.TL,c("a","b"))
  
  #3.2. TL to FL (McAuley unpublished; FL=(TL-b)/a or TL= a FL +b)
  TL.to.FL=TL.to.FL.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],c("a_TL","b_TL"),'Combined')
    TL.to.FL[[i]]=dummy$mean.par
    TL.to.FL.source[[i]]=dummy$Source
  }
  TL.to.FL=fn.names(TL.to.FL,c("a","b"))
  
  #3.3. PL to FL (Simpfendorfer et al 1999; FL=b+a*PL) 
  PL.to.FL=PL.to.FL.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],c("a_PL_FL","b_PL_FL"),'Combined')
    PL.to.FL[[i]]=dummy$mean.par
    PL.to.FL.source[[i]]=dummy$Source
  }  
  PL.to.FL=fn.names(PL.to.FL,c("a","b"))
  
  #3.4. TL-TWT (kg)     (TWT=b*TL^a)
  #Females
  TL.to.TwT.F=TL.to.TwT.F.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],c("a_w","b_w"),'Combined')
    TL.to.TwT.F[[i]]=dummy$mean.par
    TL.to.TwT.F.source[[i]]=dummy$Source
  }  
  TL.to.TwT.F=fn.names(TL.to.TwT.F,c("a","b"))
  if(SP=="WH") TL.to.TwT.F[["WH"]]=TL.to.TwT.F[["WH"]][1:2]  #use updated values provided by Rory
  if(SP=="WH") TL.to.TwT.F.source[["WH"]]=TL.to.TwT.F.source[["WH"]][1:2]  
  
  #Males
  TL.to.TwT.M=TL.to.TwT.F  #assume same length-weight for WH, BW, TK
  if(SP=="GM") TL.to.TwT.M[['GM']]=fn.get.par('GM',c("a_w","b_w"),'M')$mean.par #update male pars GM
  
  #3.5. Min and Max FL in population
  Min.Max.FL=Min.Max.FL.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],c("FL at birth"),c('Combined','F'))
    dummy1=fn.get.par(SP[i],c("maximum FL"),c('Combined','F'))
    Min.Max.FL[[i]]=c(dummy$min.par,dummy1$mean.par)
    #Min.Max.FL[[i]]=c(dummy$mean.par,dummy1$mean.par)
    Min.Max.FL.source[[i]]=c(dummy$Source,dummy1$Source)
  }  
  Min.Max.FL=fn.names(Min.Max.FL,c("Min","Max"))
  
  #3.6. Size at birth (FL,cm)
  Size.birth=Min.Max.FL
  Size.birth.source=Min.Max.FL.source
  for (i in 1:N)
  {
    Size.birth[[i]]=Size.birth[[i]][1]
    Size.birth.source[[i]]=Size.birth.source[[i]][1]
  }
  
  
  #4. Growth
  #note: 3 par von B:   TL=TL.Linf*(1-exp(-K*(age-to)))
  #     two par von B:  TL=size.birth+(TL.Linf-size.birth)*(1-exp(-k*age))
  
  #Females
  Growth.F=Growth.source=Growth.M=Max.Age.F
  Grw.par=list(WH=c("k_growth","FL_infinity","to","std.dev_growth"),GM=c("k_growth","TL_infinity","to","std.dev_growth"),
               BW=c("k_growth_2par","TL_infinity_2par"),TK=c("k_growth_2par","TL_infinity_2par"))
  id=which(names(Grw.par)==SP)
  Grw.par=Grw.par[id]
  
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],Grw.par[[i]],'F')
    Growth.F[[i]]=dummy$mean.par
    Growth.source[[i]]=dummy$Source
    dummy=fn.get.par(SP[i],Grw.par[[i]],'M')
    Growth.M[[i]]=dummy$mean.par  
  }  
  nam=list(WH=c("k","FL_inf",'to','SD'),GM=c("k","TL_inf","to","SD"),
           BW=c("k_2par","TL_inf_2par"),TK=c("k_2par","TL_inf_2par"))
  id=which(names(nam)==SP)
  nam=nam[id]
  
  for(s in 1:length(Growth.F)) 
  {
    names(Growth.F[[s]])=nam[[s]]
    Growth.F[[s]]=data.frame(t(Growth.F[[s]]))
    names(Growth.M[[s]])=nam[[s]]
    Growth.M[[s]]=data.frame(t(Growth.M[[s]]))
  }
  
  
  #5. CV of size at birth and Linf (SS3 inputs)
  if(add.growth.cv=="YES")
  {
    cv.growth.w=c(L1.f=.1,L2.f=.3,L1.m=.1,L2.m=.3)
    cv.growth.g=c(L1.f=.1,L2.f=.2,L1.m=.1,L2.m=.2)
    cv.growth.d=c(L1.f=.1,L2.f=.2,L1.m=.1,L2.m=.2)
    cv.growth.s=c(L1.f=.1,L2.f=.2,L1.m=.1,L2.m=.2)
    Growth.CV=list(WH=cv.growth.w,GM=cv.growth.g,BW=cv.growth.d,TK=cv.growth.s)
    Growth.CV.source=list(WH=NA,GM=NA,BW=NA,TK=NA)
    id=which(names(Growth.CV)==SP)
    Growth.CV=Growth.CV[id]
    id=which(names(Growth.CV.source)==SP)
    Growth.CV.source=Growth.CV.source[id]  
  }
  
  
  #6. Reproduction
  #6.1. Breeding frequency (1/cycle)
  Breed.freq=Breed.freq.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],"rep_cycl",'F')
    Breed.freq[[i]]=1/c(dummy$max.par,dummy$min.par)
    Breed.freq.source[[i]]=dummy$Source
  }  
  Breed.freq=fn.names(Breed.freq,c("Min","Max"))
  
  #6.2. Age at 50% maturity
  Age.50.mat=Age.50.mat.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],"age at 50% maturity",'F')
    Age.50.mat[[i]]=c(dummy$min.par,dummy$max.par)
    Age.50.mat.source[[i]]=dummy$Source
  }  
  Age.50.mat=fn.names(Age.50.mat,c("Min","Max"))
  
  
  #6.3. litte size 
  #Constant
  Litter.sz=Litter.sz.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],"litter size",'F')
    Litter.sz[[i]]=c(dummy$min.par,dummy$max.par)
    Litter.sz.source[[i]]=dummy$Source
  }  
  Litter.sz=fn.names(Litter.sz,c("Min","Max"))
  
  #At size
  #exponential relation : pups=exp(b+a*TL) (gummy)
  #linear relation:       pups=b+a*FL      (whiskery and sandbar)
  Litter.sz.at.size=Litter.sz.at.size.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],"a_pup",'F')
    dummy1=fn.get.par(SP[i],"b_pup",'F')
    Litter.sz.at.size[[i]]=c(dummy$mean.par,dummy1$mean.par)
    Litter.sz.at.size.source[[i]]=dummy$Source
  }  
  nam=list(WH=c("a","b"),GM=c("a","b"),BW=NA,TK=c("a","b"))
  id=which(names(nam)==SP)
  nam=nam[id]
  
  for(s in 1:length(Litter.sz.at.size)) 
  {
    if(length(Litter.sz.at.size[[s]])>0)names(Litter.sz.at.size[[s]])=nam[[s]]
  }
  
  #6.5. embryo sex ratio
  Sex.ratio=Sex.ratio.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],"pup sex ratio",'Combined')
    Sex.ratio[[i]]=dummy$mean.par
    Sex.ratio.source[[i]]=dummy$Source
  }  
  
  
  #6.6 Female size at maturity
  #note: Prop.Mat= 1/(1+exp(-(log(19)*((L-L50)/(L95-L50))))) (L= FL whiskery, TL gummy)
  Mat.50.95=Mat.50.95.source=Max.Age.F
  Mat.par.lst=list(WH=c("TL at 50% maturity","TL at 95% maturity"),
                   GM=c("TL at 50% maturity","TL at 95% maturity"),
                   BW=c("FL at 50% maturity",NA),TK=c("FL at 50% maturity",NA))
  id=which(names(Mat.par.lst)==SP)
  Mat.par.lst=Mat.par.lst[id]
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],Mat.par.lst[[i]][1],'F')
    dummy1=fn.get.par(SP[i],Mat.par.lst[[i]][2],'F')
    Mat.50.95[[i]]=c(dummy$mean.par,dummy1$mean.par)
    Mat.50.95.source[[i]]=dummy$Source
  }  
  nam=list(WH=c("L50","L95"),GM=c("L50","L95"),BW=c("L50"),TK=c("L50"))
  id=which(names(nam)==SP)
  nam=nam[id]
  for(s in 1:length(Mat.50.95)) 
  {
    if(length(Mat.50.95[[s]])>0)names(Mat.50.95[[s]])=nam[[s]]
  }
  
  
  #7. Steepness
  # mean and SD values from reference point paper. Calculated in "1.MSY proxy reference points.R"
  Steep.w=read.csv("whiskery.Steepness.csv")
  Steep.g=read.csv("gummy.Steepness.csv")
  Steep.d=read.csv("dusky.Steepness.csv")
  Steep.s=read.csv("sandbar.Steepness.csv")
  STEEP=list(WH=Steep.w,GM=Steep.g,BW=Steep.d,TK=Steep.s)
  STEEP.source=list(WH="Braccini 2015 (ref points)",GM="Braccini 2015 (ref points)",
                    BW="Braccini 2015 (ref points)",TK="Braccini 2015 (ref points)")
  id=which(names(STEEP)==SP)
  STEEP=STEEP[id]
  id=which(names(STEEP.source)==SP)
  STEEP.source=STEEP.source[id]
  
  #8. Selectivity of 6.5 inch mesh (16.5 cm) and 7 inch mesh (17.8 cm)
  #note: (sel=((FL*10/alphabeta)^alpha)*(exp(alpha-(FL*10/beta))) ; Fork length in mm)
  Selectivity=Selectivity.source=Max.Age.F
  Selectivity_7=Selectivity.source_7=Max.Age.F
  
  Sel.par=c("gillnet_selectivity_alpha","gillnet_selectivity_beta")
  Sel.par_7=c("gillnet_selectivity_alpha_7","gillnet_selectivity_beta_7")
  
  for (i in 1:N)
  {
    #6.5 inch
    dummy=fn.get.par(SP[i],Sel.par[1],'Combined')
    dummy1=fn.get.par(SP[i],Sel.par[2],'Combined')
    Selectivity[[i]]=c(dummy$mean.par,dummy1$mean.par)
    Selectivity.source[[i]]=dummy$Source
    
    #7 inch
    dummy=fn.get.par(SP[i],Sel.par_7[1],'Combined')
    dummy1=fn.get.par(SP[i],Sel.par_7[2],'Combined')
    Selectivity_7[[i]]=c(dummy$mean.par,dummy1$mean.par)
    Selectivity.source_7[[i]]=dummy$Source
    
  }  
  Selectivity=fn.names(Selectivity,c("alpha","beta"))
  Selectivity_7=fn.names(Selectivity_7,c("alpha","beta"))
  
  
  #SS3
  #Missing: estimate the 5 parameters for double logistic based on empirical curve
  Cal.double.logistic="NO"
  if(Cal.double.logistic=="YES")
  {
    library(r4ss)
    sel.pars.double.log=Selectivity
    sel.parS=list(WH=c(117.5, -5, 6, 5.9, -10, -3),GM=c(120, -3, 6.5, 6.5, -10, -4),
                  BW=NA,TK=NA)
    names(sel.parS$WH)=names(sel.parS$GM)=c("peak","top","ascwidth","descwidth","init","final")
    TL=list(WH=60:180,GM=60:180,BW=NA,TK=NA)
    id=which(names(TL)==SP)
    TL=TL[id]
    for (i in 1:2)
    {
      alpha=Selectivity[[i]]$alpha
      beta=Selectivity[[i]]$beta
      alphabeta=alpha*beta
      Tl=TL[[i]]
      
      sel.fem=((Tl*10/alphabeta)^alpha)*(exp(alpha-(Tl*10/beta)))
      plot(Tl,sel.fem,type='l',lwd=2)
      
      y=sel.line(model = 'Double_Normal', min.dist = 60, max.dist = 150,sp = sel.parS[[i]])
      sel.pars.double.log[[i]]=sel.parS[[i]]
    }
  }
  
  
  #8. Movement/Tagging
  
  #8.1 Tag shedding
  Shedding=Shedding.source=Max.Age.F
  for (i in 1:N)
  {
    dummy=fn.get.par(SP[i],"TG_shd",'Combined')
    Shedding[[i]]=dummy$mean.par
    Shedding.source[[i]]=dummy$Source
  }  
  
  #8.2 Tag reporting
  #note: these vary by zone and time. Need Rory's expertise on what values to use
  
  #non-reporting rate, Simpfendorfer et al 1999 (page 91)
  Rep.rate=read.csv(handl_OneDrive("Data/Reporting_rates_Simpfendorfer1999.csv")) 
  Rep.rate$Zn1=1-Rep.rate$Zn1
  Rep.rate$Zn2=1-Rep.rate$Zn2
  Rep.rate$WC=1-Rep.rate$WC
  Rep.rate.ag=aggregate(cbind(Zn1,Zn2,WC)~FinYear,Rep.rate,mean)
  
  #non-reporting rate, McAuley et al 2005 (page 79)
  Re.rate2=data.frame(FinYear=c("2001-02","2002-03","2003-04"),
                      WC=c(.23,.13,.09),Zn1_2=c(.5,.47,.67))   
  Re.rate2$WC=1-Re.rate2$WC
  Re.rate2$Zn1_2=1-Re.rate2$Zn1_2
  Re.rate2$Zn1=Re.rate2$Zn1_2
  Re.rate2$Zn2=Re.rate2$Zn1_2
  
  REP.RATE=rbind(Rep.rate.ag,Re.rate2[,c(1,4:5,2)])
  
  Reporting=list(WH=REP.RATE,GM=REP.RATE,BW=REP.RATE,TK=REP.RATE)              
  id=which(names(Reporting)==SP)
  Reporting=Reporting[id]
  Reporting.source=list(WH=c("Simpfendorfer et al 1999","McAuley et al 2005"),
                        GM=c("Simpfendorfer et al 1999","McAuley et al 2005"),
                        BW=c("Simpfendorfer et al 1999","McAuley et al 2005"),
                        TK=c("Simpfendorfer et al 1999","McAuley et al 2005"))
  #rep.g=0.7       #Walker 2010 (page 102), though this is for SA-Vic-Tas mostly
  id=which(names(Reporting.source)==SP)
  Reporting.source=Reporting.source[id]
  
  #9. Proportion of males in catch
  Prop.males.in.ktch=list(Prop.males.in.ktch)
  names(Prop.males.in.ktch)=SP
  Prop.males.in.ktch.source=list("DoF unpublished")
  names(Prop.males.in.ktch.source)=SP
  
  #10. List all pars
  All.pars=list(Max.Age.F=Max.Age.F,Max.Age.M=Max.Age.M,
                M=M,
                PL.to.TL=PL.to.TL,TL.to.FL=TL.to.FL,PL.to.FL=PL.to.FL,
                TL.to.TwT.F=TL.to.TwT.F,TL.to.TwT.M=TL.to.TwT.M,Min.Max.FL=Min.Max.FL,
                Size.birth=Size.birth,Growth.F=Growth.F,Growth.M=Growth.M,
                Breed.freq=Breed.freq,Age.50.mat=Age.50.mat,Litter.sz=Litter.sz,
                Litter.sz.at.size=Litter.sz.at.size,Sex.ratio=Sex.ratio,Mat.50.95=Mat.50.95,
                STEEP=STEEP,Selectivity=Selectivity,Selectivity_7=Selectivity_7,
                Shedding=Shedding,Reporting=Reporting,Prop.males.in.ktch=Prop.males.in.ktch)
  
  All.pars.source=list(Max.Age.F=Max.Age.F.source,Max.Age.M=Max.Age.M.source,
                       M=M.source,PL.to.TL=PL.to.TL.source,TL.to.FL=TL.to.FL.source,
                       PL.to.FL=PL.to.FL.source,TL.to.TwT.F=TL.to.TwT.F.source,
                       TL.to.TwT.M=TL.to.TwT.F.source,Min.Max.FL=Min.Max.FL.source,
                       Size.birth=Size.birth.source,Growth.F=Growth.source,Growth.M=Growth.source,
                       Breed.freq=Breed.freq.source,
                       Age.50.mat=Age.50.mat.source,Litter.sz=Litter.sz.source,
                       Litter.sz.at.size=Litter.sz.at.size.source,Sex.ratio=Sex.ratio.source,
                       Mat.50.95=Mat.50.95.source,STEEP=STEEP.source,
                       Selectivity=Selectivity.source,Selectivity_7=Selectivity.source_7,
                       Shedding=Shedding.source,Reporting=Reporting.source,
                       Prop.males.in.ktch.source=Prop.males.in.ktch.source)

  if(add.growth.cv=="YES" & add.Mrt.age=="YES")
  {
    All.pars=list(Max.Age.F=Max.Age.F,Max.Age.M=Max.Age.M,
                  M=M,M.age=M.age,
                  PL.to.TL=PL.to.TL,TL.to.FL=TL.to.FL,PL.to.FL=PL.to.FL,
                  TL.to.TwT.F=TL.to.TwT.F,TL.to.TwT.M=TL.to.TwT.M,Min.Max.FL=Min.Max.FL,
                  Size.birth=Size.birth,Growth.F=Growth.F,Growth.M=Growth.M,Growth.CV=Growth.CV,
                  Breed.freq=Breed.freq,Age.50.mat=Age.50.mat,Litter.sz=Litter.sz,
                  Litter.sz.at.size=Litter.sz.at.size,Sex.ratio=Sex.ratio,Mat.50.95=Mat.50.95,
                  STEEP=STEEP,Selectivity=Selectivity,Selectivity_7=Selectivity_7,
                  Shedding=Shedding,Reporting=Reporting,Prop.males.in.ktch=Prop.males.in.ktch)
    
    All.pars.source=list(Max.Age.F=Max.Age.F.source,Max.Age.M=Max.Age.M.source,
                         M=M.source,M.age=M.age.source,PL.to.TL=PL.to.TL.source,TL.to.FL=TL.to.FL.source,
                         PL.to.FL=PL.to.FL.source,TL.to.TwT.F=TL.to.TwT.F.source,
                         TL.to.TwT.M=TL.to.TwT.F.source,Min.Max.FL=Min.Max.FL.source,
                         Size.birth=Size.birth.source,Growth.F=Growth.source,Growth.M=Growth.source,
                         Growth.CV=Growth.CV.source,Breed.freq=Breed.freq.source,
                         Age.50.mat=Age.50.mat.source,Litter.sz=Litter.sz.source,
                         Litter.sz.at.size=Litter.sz.at.size.source,Sex.ratio=Sex.ratio.source,
                         Mat.50.95=Mat.50.95.source,STEEP=STEEP.source,
                         Selectivity=Selectivity.source,Selectivity_7=Selectivity.source_7,
                         Shedding=Shedding.source,Reporting=Reporting.source,
                         Prop.males.in.ktch.source=Prop.males.in.ktch.source)
   }
   
  
  #Create separate list by species
  #parameters
  PARS=vector('list',length(All.pars))
  names(PARS)=names(All.pars)
  PARS.s=PARS
  for(i in 1:length(All.pars))
  {
    PARS[[i]]=unlist(All.pars[[i]][[1]])
    if(is.data.frame(All.pars[[i]][[1]])) PARS[[i]]=All.pars[[i]][[1]]
  }
  
  #source
  for(i in 1:length(All.pars.source)) PARS.s[[i]]=unlist(All.pars.source[[i]][[1]])
  
  return(list(pars=PARS,par.source=PARS.s))
}
