# ------ Script for running shark stock assessments ---- ###################

#note:  This script firstdefine arguments used in each of the shark species/species complex assessed.
#       It then run the relevant population models according to data availability

rm(list=ls(all=TRUE))
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")
library(rlist)
library(dplyr)


# 1. Define input data and input parameters -------------------------------

  #Is this the first time the model is run?
First.run="NO"
#First.run="YES"  #set to yes to create all model input data sets and data presentation

  #Define if exporting figures as jpeg or tiff (creation of RAR requires jpeg)
Do.jpeg="YES"
Do.tiff="NO"

  #Catch in tonnes?  
KTCH.UNITS="KGS"    #used in Organise.data.R; converted to tonnes for size-based models in Run.models.R
#KTCH.UNITS="TONNES"

  #Initial bin size
MN.SZE=0
#MN.SZE="size.at.birth"

  #Minimun number of samples of size composition
Min.obs=10
Min.shts=10

  #add 25% to max size make sure no accumulation of survivals in last size class
Plus.gp.size=1.25 

  #Convert from length to age
Do.from.at.len.to.at.age="NO"

  #Reference points
B.target=.4   #single reference point for the fishery (biomass at 40% unexploited conditions)
B.threshold=.3  #default one
B.limit=.2      #default one

  #Minimum number of days at liberty for movement
MIN.DAYS.LARGE=30

  #Include a prior for r in surplus production
r_max=0.5     #max reported value of r for sharks (blue shark)
Add.r.prior=0   #no r prior 
#Add.r.prior=1   # r prior

#Define if using conventional tagging and effort in model
add.conv.tag="YES"
add.effort="NO"

  #Catch in fisheries other than NSF and TDGDLF
Data.request.1=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/1975_1987.csv",stringsAsFactors=F)
Data.request.2=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/1988_2002.csv",stringsAsFactors=F)
Data.request.3=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/2003_2016.csv",stringsAsFactors=F)
Data.request=rbind(Data.request.1,Data.request.2,Data.request.3)
Data.request=subset(Data.request,!METHOD%in%c("GN","LL"))
names(Data.request)[match("LIVEWT",names(Data.request))]="LIVEWT.c"

  #WA population for rec catch recons
WA.population=read.csv("C:/Matias/Data/AusBureauStatistics.csv",stringsAsFactors=F)

  #Jitter for estimable parameters
N=5  #how much jitter
fn.ji=function(a) jitter(a,factor=N)  #Add randomness to pin values     

  #select type of movement modellling approach
Move.mode="Individual-based"
#Move.mode="Population-based"

#Combine all size classes for tagging model
conv.tag.all="YES"    

  #Define species and year of assessment
Spec=c("Whiskery","Gummy","Dusky","Sandbar","Other")

List.sp=vector('list',length(Spec))
names(List.sp)=Spec

  #Fill in Whiskery
N.Scens=14   #Scenarios tested
List.sp$Whiskery=list(
  
  #Biology
  Spec="Whiskery",
  SPEC=17003,
  Sp2='whiskery',
  species="WH",
  TL.bins.cm=5,
  Max.FL.obs=150,          #Maximimum observed FL
  Lo=25,                   #Size at birth
  Prior.mean.Fo=0.01,
  Prior.SD.Log.Fo=0.5,
  
  #Steepness
  h.M.constant=0.351,        #Braccini et al 2015 
  h.M.constant.low=0.29,     #80% percentile from Braccini et al 2015   
  h.M.constant.up=0.41, 
  
  #Natural mortality
  M_val=0.27,          #using Hoenig set to Max age=15  (Simpfendorfer et al 2000)
  M_val.low=0.23,      #using Hoenig set to Max age=18 exp(1.46-1.01*log(18))
  M_val.high=0.35,     #using Hoenig set to Max age=12 exp(1.46-1.01*log(12))
  
  #Initial F
  Fo=0.03,             #yields a B1975 of 90% virgin           
  Fo_Simp=0.003,       #Simpfendorfer et al 2000 estimated at 0.003
  Fo_M=0.05,            #yields a B1975 of 85% virgin 
  
  Po_spm=0.9,  #Po for surplus production, consistent with the Fo value used in Size based model
  
  
  #Data
  AssessYr=2019,
  Data.yr="2017-18",         #last year of catch
  Frst.yr.ktch="1975-76",    #first year of catch)
  AREAS=c("West","Zone1","Zone2"),  #Define spatial areas; 1 is West, 2 is Zn1, 3 is Zn2.
  Yr_q_change=1982,   #last year before targeting practices changed (Simpfendorfer 2000)
  Yr_q_daily=2006,
  Drop_yr_cpue=c("1980-81","1981-82","1982-83","1983-84"),  #Dropped cpue years
  Drop_yr_cpue_sens=c("1975-76","1976-77","1977-78","1978-79","1979-80",
                                          "1980-81","1981-82","1982-83"),
  Do_var=0,    #How to calculate cpue variance in age-structured
  Var1=0.0296133485565278,   #fixed variances used by Simpfendorfer et al 2000
  Var2=0.023525088856464,
  
  #Scenarios
  N.Scens=N.Scens,       
  Zens=paste("S",1:(N.Scens-1),sep=""),
  Models=c("Base case",paste("S",1:(N.Scens-1),sep="")),
  Q.scen=c(rep("three",2),"two","three",rep("two",2),rep("three",8)),
  
  #Initial value of estimable parameters
    #Dummy for switching on/off phase estimation in simulation testing
  Dummy=1,   
  
    #Biomass dynamics
  r_BD=fn.ji(0.2),
  k_BD=fn.ji(5000),
  Q1_BD=fn.ji(7e-4),
  Q2_BD=fn.ji(2e-4),
  daily_BD=fn.ji(4e-4),
  tau2_BD=fn.ji(0.1353),
  
    #Age-structured
  RSTAR_AS=fn.ji(55),
  z_AS=fn.ji(.6),
  Q1_AS=fn.ji(.4),
  Q2_AS=fn.ji(.2),
  q_daily_AS=fn.ji(0.8),
  Fo_AS=fn.ji(0.008),
  
    #Size-base
  RZERO_in_1000_individuals_SB=fn.ji(1096),
  Q1_SB=fn.ji(0.00041),
  Q2_SB=fn.ji(0.0001),
  Q_daily_SB=fn.ji(0.00024),
  Fo_SB=NA,  #no jit because it's fixed
  tau_SB=fn.ji(0.2) ,
  K.F=fn.ji(0.38),
  Linf.F=fn.ji(130),
  K.M=fn.ji(0.423),
  Linf.M=fn.ji(130),
  SD.growth_SB=fn.ji(7),
  Prop.west=fn.ji(0.19),  #proportion of recruitment by zone calculated based on average prop catch by zone
  Prop.zn1=fn.ji(0.38),
  
    #movement 
  p11=0.99,
  p22=0.99,
  p21=0.001,
  p33=0.99,
  #mov.WC_WC=.9,
  #mov.WC_ZN1=.1,
  #mov.ZN1_WC=.1,
  #mov.ZN1_ZN1=.9,
  #mov.ZN2_ZN1=.1,
  #mov.ZN2_ZN2=.9,
  #log.Q.tag.WC=-5,
  #log.Q.tag.ZN1=-5,
  #log.Q.tag.ZN2=-5,
  #log.tau=1
  
    #catch-MSY arguments
  Growth.F=data.frame(k=0.369,FL_inf=121),
  TEMP=18,
  Max.age.F=c(15,19),
  Age.50.mat=c(6,7),
  Fecundity=c(4,28),
  Breed.cycle=2,  #years
  years.futures=5
  
  
  )


  #fill in Gummy
  
  #fill in Dusky
  
  #fill in Sandbar
  
  #fill in Other shark species



# 2. Define Scenarios for sensitivity tests -------------------------------

BaseCase="Size-based"
#BaseCase="Age-based"

  #Do Scenario comparison in colors or greyscale
Do.cols="YES"  
#Do.cols="NO"

  #Define scenarios
SIMS=1e5
Proc.err=0.05 #sigR is PROCESS ERROR; 0 if deterministic model

      #fill in Whiskery
List.sp$Whiskery=with(List.sp$Whiskery,
     {
       list.append(List.sp$Whiskery,
          Drop_yr_cpue.tabl=paste(substr(Drop_yr_cpue,1,4)[1],substr(Drop_yr_cpue,3,4)[length(Drop_yr_cpue)],sep="-"),
        Drop_yr_cpue_sens.tabl=paste(substr(Drop_yr_cpue_sens,1,4)[1],substr(Drop_yr_cpue_sens,3,4)[length(Drop_yr_cpue_sens)],sep="-"))
      })
List.sp$Whiskery=with(List.sp$Whiskery,
    {
      list.append(List.sp$Whiskery,
              Tabla.scen=data.frame(
                 Model=Models,
                 Size_comp.=c('Yes',"N/A",'No',rep('Yes',11)),
                 CPUE=c(rep("Stand.",2),"Effective",rep('Stand.',2),"Effective",rep('Stand.',8)),
                 CPUE_years_dropped=c(rep(Drop_yr_cpue.tabl,2),rep("None",2),Drop_yr_cpue_sens.tabl,"None",rep(Drop_yr_cpue.tabl,8)),
                 Age.Growth=c('Yes',"N/A",'No',rep('Yes',11)),
                 Ktch.sx.r=c('Observed','N/A','Equal',rep('Observed',3),'Equal',rep('Observed',7)),                      
                 Tagging=c('No','N/A',rep('No',11),'Yes'),                     
                 Fec.=c(rep('N/A',2),rep('constant',1),rep('N/A',11)),
                 Maturity=c('at length','N/A',rep('knife edge',1),rep("at length",11)),
                 M=c("constant","N/A",rep("constant",12)),
                 M.value=c(M_val,NA,rep(M_val,5),M_val.low,M_val.high,rep(M_val,5)),                      
                 SteepnesS=c(h.M.constant,rep("N/A",2),rep(h.M.constant,8),h.M.constant.low,h.M.constant.up,h.M.constant),
                 Q=Q.scen,   
                 Spatial_structure=c(rep('Single zone',13),'Three zones'),
                 Movement=c("No",rep("N/A",2),rep("No",10),"Yes"),
                 Fo=c(Fo,"N/A","estimated",rep(Fo,6),Fo_Simp,Fo_M,rep(Fo,3)),
                 Model_type=c('Length-based','Biomass dynamics',"Age-structured",rep("Length-based",11))
             )
        )
    })
ktch_msy_scen=list()
ktch_msy_scen$'Base case'=list(
  r.prior="USER",    
  user="Yes",
  k.max=50,
  startbio=c(0.7,.95),
  finalbio=c(0.2, 0.7),
  res="low",
  niter=SIMS,
  sigR=Proc.err)
ktch_msy_scen$S1=list(
  r.prior="USER",    
  user="Yes",
  k.max=50,
  startbio=c(0.7,.95),
  finalbio=c(0.2, 0.7),
  res="low",
  niter=SIMS,
  sigR=0.02)
ktch_msy_scen$S2=list(
  r.prior=NA,  
  user="No",
  k.max=50,
  startbio=c(0.7,.95),
  finalbio=c(0.2, 0.7),
  res="low",     
  niter=SIMS,
  sigR=Proc.err)
List.sp$Whiskery=list.append(List.sp$Whiskery,ktch_msy_scen=ktch_msy_scen)

      #fill in Gummy

      #fill in Dusky
      
      #fill in Sandbar
      
      #fill in Other shark species



  #Add number of areas and handle, create path/folders, add fixed Fo
for(l in 1:length(List.sp))
{
  n.areas=length(List.sp[[l]]$AREAS)
  List.sp[[l]]$n.areas=n.areas
  List.sp[[l]]$Areas.zones=data.frame(area=1:n.areas,zone=List.sp[[l]]$AREAS)
  
  List.sp[[l]]$hndl=paste("C:/Matias/Analyses/Population dynamics/",Spec[l]," shark/",sep='')
  
  List.sp[[l]]$Fo_SB=List.sp[[l]]$Fo  #fixed
  
  with(List.sp[[l]],
  {
    #create folders if new run
     kriat.path=paste(hndl,AssessYr,sep="")
     if(!dir.exists(kriat.path))dir.create(kriat.path)
     kriate.this=as.character(Tabla.scen$Model)
     for (i in 1:length(kriate.this)) 
     { 
         NEW=paste(kriat.path,kriate.this[i], sep="/")
         if(!dir.exists(NEW))dir.create(NEW) 
     }
     
     #drop conv tag if not using
     if(add.conv.tag=="NO") List.sp[[l]]$Tabla.scen=Tabla.scen[,-match("Tagging",names(Tabla.scen))]
     
     
     #Edit scenarios table
     inpt.pz=paste(hndl,AssessYr,"/1_Inputs",sep="")
     if(!dir.exists(inpt.pz))dir.create(inpt.pz)
     setwd(inpt.pz)
  })
  
  #output scenarios table
  if(First.run=="YES" )
  {
    
    Tabla.scen.show=List.sp[[l]]$Tabla.scen
    Tabla.scen.show=Tabla.scen.show[,-match(c("Fec.","M"),colnames(Tabla.scen.show))]
    Tabla.scen.show[is.na(Tabla.scen.show)]="N/A"
    names(Tabla.scen.show)[match(c("Model","SteepnesS"),names(Tabla.scen.show))]=c("Model_name","Steepness")
    Tabla.scen.show=Tabla.scen.show[,match(c("Model_name","Model_type","Spatial_structure"
                                             ,"Movement","Size_comp.","CPUE","CPUE_years_dropped","Age.Growth"         
                                             ,"Ktch.sx.r","Tagging","M.value","Steepness","Fo","Maturity","Q"),names(Tabla.scen.show))]
    Tabla.scen.show=Tabla.scen.show[match(c(Zens,"Base case" ),Tabla.scen.show$Model_name),]
    Tabla.scen.show$Q=as.character(Tabla.scen.show$Q)
    Tabla.scen.show$Q=with(Tabla.scen.show,ifelse(Q=="three",3,ifelse(Q=="two",2,ifelse(Q=="one",1,Q))))
    Tabla.scen.show$Q=with(Tabla.scen.show,ifelse(!Q=="N/A"& Q>1,paste(Q,"periods"),
                                                  ifelse(!Q=="N/A"& Q==1,paste(Q,"period"),Q))) 
    
    
    Scenarios.tbl=function(WD,Tbl,Doc.nm,caption,paragph,HdR.col,HdR.bg,Hdr.fnt.sze,Hdr.bld,
                           body.fnt.sze,Zebra,Zebra.col,Grid.col,Fnt.hdr,Fnt.body,
                           HDR.names,HDR.span,HDR.2nd,HDR.3rd)
    {
      mydoc = docx(Doc.nm)  #create r object
      mydoc = addSection( mydoc, landscape = T )   #landscape table
      # add title
      if(!is.na(caption))mydoc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
      
      # add a paragraph
      if(!is.na(paragph))mydoc = addParagraph(mydoc , paragph, stylename="Citationintense")
      
      #add table
      MyFTable=FlexTable(Tbl,header.column=F,add.rownames =F,
                         header.cell.props = cellProperties(background.color=HdR.bg), 
                         header.text.props = textProperties(color=HdR.col,font.size=Hdr.fnt.sze,
                                                            font.weight="bold",font.family =Fnt.hdr), 
                         body.text.props = textProperties(font.size=body.fnt.sze,font.family =Fnt.body))
      
      #Add header
      MyFTable = addHeaderRow(MyFTable,text.properties=textBold(),value=HDR.names,colspan=HDR.span)
      
      #Add second header
      MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value =HDR.2nd)
      
      #Add third header
      MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value =HDR.3rd)
      
      # zebra stripes - alternate colored backgrounds on table rows
      if(Zebra=="YES") MyFTable = setZebraStyle(MyFTable, odd = Zebra.col, even = "white" )
      
      # table borders
      MyFTable = setFlexTableBorders(MyFTable,
                                     inner.vertical = borderNone(),inner.horizontal = borderNone(),
                                     outer.vertical = borderNone(),
                                     outer.horizontal = borderProperties(color=Grid.col, style="solid", width=4))
      
      # set columns widths (in inches)
      #MyFTable = setFlexTableWidths( MyFTable, widths = Col.width)
      
      mydoc = addFlexTable( mydoc, MyFTable)   
      mydoc = addSection( mydoc, landscape = F ) 
      
      # write the doc 
      writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
    }
    options('ReporteRs-fontsize'= 12, 'ReporteRs-default-font'='Arial')   
    Scenarios.tbl(WD=getwd(),Tbl=Tabla.scen.show,Doc.nm="Model scenarios",
                  caption=NA,paragph=NA,HdR.col='black',HdR.bg='white',
                  Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                  Zebra='NO',Zebra.col='grey60',Grid.col='black',
                  Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
                  HDR.names=c('Model','Spatial','Movement', 'Data','Input parameters','Q'),
                  HDR.span=c(2,1,1,6,4,1),
                  HDR.2nd=c("Name","Type","structure",'',"Size","CPUE","CPUE years","Age &",
                            "Prop. male","Tagging","M","h","Fo","Maturity",""),
                  HDR.3rd=c("","","","","composition","","not used in likelihood","growth","in catch",
                            "","","","","",""))
  }
}
  

#Define groups for comparing model outputs 
if(First.run=="YES" )
{
  MatCh=function(what,where) as.character(Tabla.scen$Model[match(what,where)])
  Compr.grup=list(
    c("Base case",MatCh(c("Biomass dynamics","Age-structured"),Tabla.scen$Model_type),
      MatCh("Three zones",Tabla.scen$Spatial_structure)),
    c("Base case","S3","S4"),
    c("Base case","S5"),
    c("Base case","S6"),
    c("Base case",MatCh(c(M_val.low,M_val.high),Tabla.scen$M.value)),
    c("Base case",MatCh(c(Fo_Simp,Fo_M),Tabla.scen$Fo)),
    c("Base case",MatCh(c(h.M.constant.low,h.M.constant.up),Tabla.scen$SteepnesS))
  )
  names(Compr.grup)=c("Model_type","CPUE_years_dropped","CPUE","Ktch.sx.r","M.value","Fo","SteepnesS")
  Title.Compr.grup=c("Model structure","CPUE years not used in likelihood","Effective vs standardised CPUE",
                     "Proportion of males in catch","M","Fo","Steepness")
}



# 3. Export pin file and determine estimation phases-------------------------------

#note: S1 and S2 calculates pars in normal space but same order magnitude
#       Other scenarios all pars in log.
#       ln_RZERO is in 1,000 individuals so do 10 times the largest catch divided by 
#       average weight and divided by 1,000. Best units to work in are 1,000 individuals for 
#       numbers, catch in tonnes and length-weight in kg as all cancels out and predicted biomasses
#       end up being in tonnes

Fz.off=-22
fn.mtch=function(WHAT,NMS) match(WHAT,names(NMS))
Q_phz=c("lnq","lnq2","log_Qdaily")                           
Zns.par.phz=c("lnR_prop_west","lnR_prop_zn1")
MOv.par.phz=c("log_p11","log_p22","log_p21","log_p33")
for(l in 1:length(List.sp))
{
  attach(List.sp[[l]])
  
  #Population pin values
  Pin.pars=vector('list',nrow(Tabla.scen))
  names(Pin.pars)=Tabla.scen$Model
  
    #Biomass dynamics
  Pin.pars$S1=c(Dummy=Dummy,r=r_BD,log_k=log(k_BD),Q1=Q1_BD,Q2=Q2_BD,Qdaily=daily_BD,log_tau2=log(tau2_BD))
  
    #Age-structured
  Pin.pars$S2=c(Dummy=Dummy,RSTAR=RSTAR_AS,Z=z_AS,Q1=Q1_AS,Q2=Q2_AS,Qdaily=q_daily_AS,Fo=Fo_AS) 
  
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
  
    #set the pins from other size-based scenarios
  IDs=which(Tabla.scen$Model_type=="Length-based")
  IDs=IDs[-match("Base case",Tabla.scen$Model[IDs])]
  for(i in IDs) Pin.pars[[i]]=Pin.pars$'Base case'
  
  for(i in 1:N.Scens)     
  {
    if(Tabla.scen$Model_type[i]=="Length-based")
    {
      #higher M or Fo needs larger initial population
      if(Tabla.scen$M.value[i]==M_val.high | Tabla.scen$Fo[i]==Fo_M)    
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
  setPath=function(Scen)setwd(paste(List.sp[[l]]$hndl,List.sp[[l]]$AssessYr,"/",Scen,sep=""))
  for(i in 1:nrow(Tabla.scen))
  {
    These.pars=Pin.pars[[i]]
    par.nms=names(These.pars)
    setPath(Tabla.scen[i,]$Model)
    FILE=paste(List.sp[[l]]$Spec,".pin",sep="")
    write("# Input parameters",file = FILE)
    for(k in 1:length(These.pars))
    {
      Hdr=paste("#",par.nms[k])
      write(Hdr,file = FILE,append=T)
      write(These.pars[k],file = FILE,sep = "\t",append=T)
    }
  }
  List.sp[[l]]$Pin.pars=Pin.pars
  
  
  #Determine estimation phases
  Par.phases=Pin.pars
  Par.phases$S1=c(dummy=Fz.off,r=2,ln_k=2,ln_q1=1,ln_q2=1,ln_qdaily=-1,ln_tau=3)
  Par.phases$S2=c(dummy=Fz.off,Rstar=2,z=3,q1=1,q2=1,qdaily=Fz.off,Fo=4)
  Par.phases$'Base case'=c(dummy=-1,
                           lnR_zero=1,
                           lnR_prop_west=1,
                           lnR_prop_zn1=1,
                           lnq=2,                           
                           lnq2=2,
                           log_Qdaily=2,
                           ln_Init_F=Fz.off,      #cannot estimate Fo
                           log_tau=5,
                           k=3,
                           lnLinf=3,
                           k_M=3,
                           lnLinf_M=3,
                           sd_growth=4,
                           log_p11=1,
                           log_p22=1,
                           log_p21=1,
                           log_p33=1)
  
  IDs=(1:length(Par.phases))[-match(c("Base case","S1","S2"),names(Par.phases))]
  for(i in IDs) Par.phases[[i]]=Par.phases$'Base case'
  
  
  #Turn off irrelevant pars according to scenarios
  for(i in 1:N.Scens)     
  {
    if(Tabla.scen$Model_type[i]=="Length-based")
    {
      #single zone
      if(Tabla.scen$Spatial_structure[i]=="Single zone")
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
    
    
  detach(List.sp[[l]])
  
}

Pin.pars=1
Par.phases=1

# 4. Define modelling arguments -------------------------------

  #What effort to use?
#What.Effort="km.gn.days" 
What.Effort="km.gn.hours" 


  #Select acoustic tagging model
Acoust.format=Move.mode
#Acoust.format="SS3"

  #Select type of size composition Likelihood
Size_like="Dirichlet"
#Size_like="Multinomial"
Dirichlet.small.value=1e-4   #small constant for intermediate 0 observations 

  #Combined size composition?
Size.sex.comb="NO"

  #Size compostion as proportions?
Size.comp.prop="YES"

  #Number of years for future projections
Yrs.future=5

  #Effective sample size size composition
Effective.n=300

  #Maximum possible F
MaxF=3.0

  #Weight of likelihood components
for(l in 1:length(List.sp))
{
  #weight of size comp 
  List.sp[[l]]$Rho=1
  
  #weight of age-length 
  List.sp[[l]]$Rho2=1
  
  #cpue likelihood
  x=rep(1,length(List.sp[[l]]$Models))
  names(x)=as.character(List.sp[[l]]$Models)
  dd=which(List.sp[[l]]$Tabla.scen$CPUE%in%c("N/A","No")) #adjust weight of cpue if cpue like not used in model
  if(length(dd)>0)x[dd]=0
  List.sp[[l]]$Rho3=x
}

  #Prior for initial fishing mortality
add_Finit_prior=0  #no prior
#add_Finit_prior=1  #add prior


  #What scenarios to run
if(First.run=="YES") Run.all.Scenarios="YES"
if(First.run=="NO") Run.all.Scenarios="NO"  #base case only


Arguments=""  #run without arguments
#Arguments=" -est -sim "  #run with these arguments

Do.zero.Ktch="NO" #do future projections with 0 catch?
if(Do.zero.Ktch=="YES")
{
  zero.Ktch.yrs=100
  zero.Ktch.this="Base case"
}

Do.sim.test="NO" #do simulation testing?
if(Do.sim.test=="YES")
{
  N.sim.test=50
  CPUE.jitr=1000 
  Sim.Test.this="Base case"
}

  #Control if doing MSY   
Do.MSY="NO"
MSY.yrs=100  
MSY.sd.rec=0    #no recruitment deviations (in log space)
#MSY.sd.rec=0.05
MSY.sims=100
F.vec=seq(0,.2,by=.01)

  #MCMC 
DO.MCMC="YES"
nSims=1e6
Thin=10   
burning=1:(5*length(seq(1,nSims,by=Thin))/100)   #5%  burning

  #Future projection scenarios
if(First.run=="NO") Run.future.proj="YES"
Future.ktch.scen=list(mean=1,percent_50=.5,percent_150=1.5)




# 5. Define how to display model outputs -------------------------------
  #What years to show in model outputs
Show.yrs="DATA"  
#Show.yrs="FUTURE"  #show data plus future projections

  #present cpue in log space or normal space
Present.in.log="NO"   



# 6. Execute population models -------------------------------
#note: 'Run.models.R' outputs data and parameter inputs for models,
#                     runs assessment models
#                     and displays outputs
for(l in 1:length(List.sp))
{
  attach(List.sp[[l]])
  source("C:/Matias/Analyses/SOURCE_SCRIPTS/Population dynamics/Run.models.R")
  detach(List.sp[[l]])
}



