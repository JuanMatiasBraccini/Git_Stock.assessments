# ------ Script for running gummy shark stock assessment ---- ###################

#Index:
    #A. Define arguments used in 'Run.models.R' (source script)
        #A.1. Data inputs and input parameter arguments
        #A.2. Define Scenarios for sensitivity tests
        #A.3. Define initial values of estimable parameters
        #A.4. Modelling arguments
        #A.5. Output display
        
    #B. Execute 'Run.models.R'



rm(list=ls(all=TRUE))
source("C:/Matias/Analyses/SOURCE_SCRIPTS/MS.Office.outputs.R")


#-- A. Define arguments used in 'Run.models.R' (source script) 


#A.1. Data inputs and input parameter arguments

#Is this the first time the model is run?
First.run="NO"
#First.run="YES"

#Select species and year of assessment
SPEC=17001     
Spec="Gummy"
Sp2='gummy'
species="GM"
AssessYr=2017
TL.bins.cm=5
Data.yr="2015-16"         #last year of catch
Frst.yr.ktch="1975-76"    #first year of catch

#Define if exporting figures as jpeg or tiff (creation of RAR requires jpeg)
Do.jpeg="YES"
Do.tiff="NO"

#Define base case model
BaseCase="Size-based"
#BaseCase="Age-based"

#Define if using conventional tagging and effort in model
add.conv.tag="YES"

#Catch in tonnes?  
KTCH.UNITS="KGS"
#KTCH.UNITS="TONNES"

#Use effort or not?
add.effort="NO"

#what effort to use?
What.Effort="km.gn.days" 
#What.Effort="km.gn.hours" 

#Combine all size classes for tagging model
conv.tag.all="YES"    

#select type of movement modellling approach
Move.mode="Individual-based"
#Move.mode="Population-based"

  #Select acoustic tagging model
Acoust.format=Move.mode
#Acoust.format="SS3"

  #Initial bin size
MN.SZE=0
#MN.SZE="size.at.birth"


  #Minimun number of samples of size composition
Min.obs=10
Min.shts=10

  #Select type of size composition Likelihood
Size_like="Dirichlet"
#Size_like="Multinomial"
Dirichlet.small.value=1e-4   #small constant for intermediate 0 observations 

  #Maximimum observed FL
Max.FL.obs=175 

  #Size at birth
Lo=33

#add 25% to max size make sure no accumulation of survivals in last size class
Plus.gp.size=1.25           

#Convert from length to age
Do.from.at.len.to.at.age="NO"

#combined size composition?
Size.sex.comb="NO"

#size compostion as proportions?
Size.comp.prop="YES"

#Define spatial areas
AREAS=c("West","Zone1","Zone2") #note: 1 is West, 2 is Zn1, 3 is Zn2.


#Number of years for future projections
Yrs.future=5

#Effective sample size size composition
Effective.n=300

#Scenarios tested
N.Scens=6
Zens=paste("S",1:(N.Scens-1),sep="")
Models=c("Base case",Zens)

#weight of size comp likelihood
Rho=1 

#weight of age-length likelihood
Rho2=1 

  #weight of cpue likelihood              
Rho3=rep(1,length(Models))
names(Rho3)=as.character(Models)


#maximum possible F
MaxF=3.0

#Prior for initial fishing mortality
add_Finit_prior=0  #no prior
#add_Finit_prior=1  #add prior
Prior.mean.Fo=0.01
Prior.SD.Log.Fo=0.5

#Reference points
B.target=.4   #single reference point for the fishery (biomass at 40% unexploited conditions)
B.threshold=.3  #default one
B.limit=.2      #default one

#Minimum number of days at liberty for movement
MIN.DAYS.LARGE=30

Yr_q_change=0   #dummy
Yr_q_daily=2006

#Po for surplus production
Po_spm=0.95  #consistent with the Fo value used in Size based model

#Include a prior for r in surplus production
r_max=0.5     #max reported value of r for sharks (blue shark)
Add.r.prior=0   #no r prior 
#Add.r.prior=1   # r prior


#How to calculate cpue variance in age-structured
Do_var=0
Var1=0.029
Var2=Var1


#A.2. Define Scenarios for sensitivity tests

#Select number of areas
n.areas=length(AREAS)  
Areas.zones=data.frame(area=1:n.areas,zone=AREAS)


#Catch in fisheries other than NSF and TDGDLF
Data.request.1=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/1975_1987.csv",stringsAsFactors=F)
Data.request.2=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/1988_2002.csv",stringsAsFactors=F)
Data.request.3=read.csv("C:/Matias/Data/Catch and Effort/Data_request_12_2014/2003_2016.csv",stringsAsFactors=F)
Data.request=rbind(Data.request.1,Data.request.2,Data.request.3)
Data.request=subset(Data.request,!METHOD%in%c("GN","LL"))
names(Data.request)[match("LIVEWT",names(Data.request))]="LIVEWT.c"


#WA population for rec catch recons
WA.population=read.csv("C:/Matias/Data/AusBureauStatistics.csv",stringsAsFactors=F)

#Model scenarios                    
hndl=paste("C:/Matias/Analyses/Population dynamics/",Spec," shark/",sep='')

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


h.M.constant=0.481  #Braccini et al 2015 
h.M.constant.low=0.461    #80% percentile
h.M.constant.up=0.5
M_val=0.283          #Walker empirical
M_val.low=0.22      #using Hoenig set to Max age=19 exp(1.46-1.01*log(19))
M_val.high=0.32     #using Hoenig set to Max age=13 exp(1.46-1.01*log(13))


Fo=0.05                #This leaves B1975 at 95% unfished 
Fo_Simp=0.003              
Fo_M=0.1                
Fo_AS=0.004             #This leaves B1975 at 95% unfished 


Q.scen=c(rep("two",2),"one","N/A","two","two")

Tabla.scen=data.frame(
  Model=Models,
  Size_comp.=c('Yes',"N/A",'No','Yes','No','Yes'),
  CPUE=c(rep("Stand.",2),"Effective","No","Stand.","Stand.hours"),
  Age.Growth=c('Yes',"N/A",'No','Yes','No','Yes'),
  Ktch.sx.r=c('Observed','N/A','Equal','Observed','Equal','Observed'),
  Tagging=c('No','N/A',rep('No',4)),
  Fec.=c(rep('N/A',2),'constant','N/A','constant','N/A'),
  Maturity=c('at length','N/A','knife edge',"at length",'knife edge',"at length"),
  M=c("constant","N/A",rep("constant",4)),
  M.value=c(M_val,NA,rep(M_val,4)),
  SteepnesS=c(h.M.constant,rep("N/A",2),h.M.constant,"N/A",h.M.constant),
  Q=Q.scen,   
  Spatial_structure=rep('Single zone',6),
  Movement=c("No","N/A",rep("No",4)),
  Fo=c(Fo,"N/A",Fo_AS,Fo,Fo_AS,Fo),
  Model_type=c('Length-based','Biomass dynamics',"Age-structured","Length-based","Age-structured","Length-based")
)

#create folders if new run
kriat.path=paste(hndl,AssessYr,sep="")
kriate.this=as.character(Tabla.scen$Model)
for (i in 1:length(kriate.this)) 
{ 
  NEW=paste(kriat.path,kriate.this[i], sep="/")
  if(!dir.exists(NEW))dir.create(NEW) 
}

if(add.conv.tag=="NO") Tabla.scen=Tabla.scen[,-match("Tagging",names(Tabla.scen))]
setwd(paste(hndl,AssessYr,"/1_Inputs",sep=""))
Tabla.scen.show=Tabla.scen
Tabla.scen.show=Tabla.scen.show[,-match(c("Fec.","M"),colnames(Tabla.scen.show))]
Tabla.scen.show[is.na(Tabla.scen.show)]="N/A"
names(Tabla.scen.show)[match(c("Model","SteepnesS"),names(Tabla.scen.show))]=c("Model_name","Steepness")
Tabla.scen.show=Tabla.scen.show[,match(c("Model_name","Model_type","Spatial_structure"
             ,"Movement","Size_comp.","CPUE","Age.Growth"         
             ,"Ktch.sx.r","Tagging","M.value","Steepness","Fo","Maturity","Q"),names(Tabla.scen.show))]
Tabla.scen.show=Tabla.scen.show[match(c(Zens,"Base case" ),Tabla.scen.show$Model_name),]
Tabla.scen.show$Q=as.character(Tabla.scen.show$Q)
Tabla.scen.show$Q=with(Tabla.scen.show,ifelse(Q=="three",3,ifelse(Q=="two",2,ifelse(Q=="one",1,Q))))
Tabla.scen.show$Q=with(Tabla.scen.show,ifelse(!Q=="N/A"& Q>1,paste(Q,"periods"),
                                              ifelse(!Q=="N/A"& Q==1,paste(Q,"period"),Q))) 

options('ReporteRs-fontsize'= 12, 'ReporteRs-default-font'='Arial')   
Scenarios.tbl(WD=getwd(),Tbl=Tabla.scen.show,Doc.nm="Model scenarios",
              caption=NA,paragph=NA,HdR.col='black',HdR.bg='white',
              Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
              HDR.names=c('Model','Spatial','Movement', 'Data','Input parameters','Q'),
              HDR.span=c(2,1,1,5,4,1),
              HDR.2nd=c("Name","Type","structure",'',"Size","CPUE","Age &",
                        "Prop. male","Tagging","M","h","Fo","Maturity",""),
              HDR.3rd=c("","","","","composition","","growth","in catch",
                        "","","","","",""))

#adjust weight of cpue if cpue like not used in model
dd=which(Tabla.scen$CPUE%in%c("N/A","No"))
if(length(dd)>0)Rho3[dd]=0


#Define groups for comparing model outputs     
MatCh=function(what,where) as.character(Tabla.scen$Model[match(what,where)])
Compr.grup=list(
  c("Base case","S1","S2","S4"),
  c("Base case","S3","S5")
)
names(Compr.grup)=c("Model_type","CPUE")
Title.Compr.grup=c("Model structure","CPUE")

Do.cols="YES"  #Do Scenario comparison in colors or greyscale
#Do.cols="NO"


#A.3. Define initial values of estimable parameters
N=5  #how much jitter
fn.ji=function(a) jitter(a,factor=N)  #Add randomness to pin values 


  #Dummy    #used for switching on/off phase estimation in simulation testing   
Dummy=1

  #Biomass dynamics
r_BD=fn.ji(0.3)
k_BD=fn.ji(5000)
Q1_BD=fn.ji(1e-7)
Q2_BD=fn.ji(1e-7)
daily_BD=fn.ji(1e-7)
tau2_BD=fn.ji(0.1353)


  #Age-structured
RSTAR_AS=100
z_AS=1
Q1_AS=1.4
Q2_AS=fn.ji(0.8)
q_daily_AS=fn.ji(0.8)
      

  #Size-base
RZERO_in_1000_individuals_SB=fn.ji(1000)
Q1_SB=fn.ji(1e-4)
Q2_SB=fn.ji(1e-4)
Q_daily_SB=fn.ji(1e-4)
Fo_SB=Fo   #no jit because it's fixed
tau_SB=fn.ji(0.3) 

K.F=fn.ji(0.15)
Linf.F=fn.ji(180)
K.M=fn.ji(0.25)
Linf.M=fn.ji(150)
SD.growth_SB=fn.ji(10)
Prop.west=fn.ji(0.02)   #proportion of recruitment by zone calculated based on average prop catch by zone
Prop.zn1=fn.ji(0.08)

#add tagging pars 
if(Move.mode=="Individual-based")
{
  p11=0.999
  p22=0.999
  p21=0.00024
  p33=0.999
}
if(Move.mode=="Population-based")
{
  mov.WC_WC=.9;mov.WC_ZN1=.1
  mov.ZN1_WC=.1;mov.ZN1_ZN1=.9
  mov.ZN2_ZN1=.1;mov.ZN2_ZN2=.9
  log.Q.tag.WC=-5;log.Q.tag.ZN1=-5;log.Q.tag.ZN2=-5;log.tau=1  
}



#Put pin pars in list                    
#note: S1 and S2 calculates pars in normal space but same order magnitude
#       Other scenarios all pars in log.
#       ln_RZERO is in 1,000 individuals so do 10 times the largest catch divided by 
#       average weight and divided by 1,000. Best units to work in are 1,000 individuals for 
#       numbers, catch in tonnes and length-weight in kg as all cancels out and predicted biomasses
#       end up being in tonnes

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
                       #lnq_west=log(Q1_SB),lnq_zn1=log(Q1_SB),lnq_zn2=log(Q1_SB),
                       #lnq2_west=log(Q2_SB),lnq2_zn1=log(Q2_SB),lnq2_zn2=log(Q2_SB),
                       #lnq_daily_west=log(Q_daily_SB),lnq_daily_zn1=log(Q_daily_SB),lnq_daily_zn2=log(Q_daily_SB),
                       ln_Init_F=log(Fo_SB),log_tau=log(tau_SB),
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

#set the pins from other age-strucutred scenarios
IDs=which(Tabla.scen$Model_type=="Age-structured")
IDs=IDs[-match("S2",Tabla.scen$Model[IDs])]
for(i in IDs) Pin.pars[[i]]=Pin.pars$'S2'


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
    if(as.numeric(as.character(Tabla.scen$Fo[i]))>Fo) 
    {
      Pin.pars[[i]][match("ln_Init_F",names(Pin.pars[[i]]))]=log(as.numeric(as.character(Tabla.scen$Fo[i])))
      ss=match("lnR_zero",names(Pin.pars[[i]]))
      Pin.pars[[i]][ss]=Pin.pars[[i]][ss]*1.15#higher Fo needs larger initial population
    }
    
  }
}


#Export as pin file
setPath=function(Scen)setwd(paste(hndl,AssessYr,"/",Scen,sep=""))
for(i in 1:nrow(Tabla.scen))
{
  These.pars=Pin.pars[[i]]
  par.nms=names(These.pars)
  setPath(Tabla.scen[i,]$Model)
  FILE=paste(Spec,".pin",sep="")
  write("# Input parameters",file = FILE)
  for(k in 1:length(These.pars))
  {
    Hdr=paste("#",par.nms[k])
    write(Hdr,file = FILE,append=T)
    write(These.pars[k],file = FILE,sep = "\t",append=T)
  }
}


#Estimable parameter phases
Fz.off=-22
Par.phases=Pin.pars
Par.phases$S1=c(dummy=Fz.off,r=2,ln_k=2,ln_q1=1,ln_q2=Fz.off,ln_qdaily=1,ln_tau=3)
Par.phases$S2=c(dummy=Fz.off,Rstar=1,z=1,q1=1,q2=Fz.off,qdaily=Fz.off,Fo=-4) 

Par.phases$'Base case'=c(dummy=Fz.off,
                         lnR_zero=3,
                         lnR_prop_west=-3,
                         lnR_prop_zn1=-3,
                         lnq=4,
                         lnq2=Fz.off,
                         # lnq_west=5,
                         # lnq_zn1=5,
                         # lnq_zn2=5,
                         # lnq2_west=-5,
                         # lnq2_zn1=-5,
                         # lnq2_zn2=-5,
                         log_Qdaily=4,
                         # log_Qdaily_west=5,
                         # log_Qdaily_zn1=5,
                         # log_Qdaily_zn2=5,
                         ln_Init_F=Fz.off,              #cannot estimate Fo
                         log_tau=5,
                         k=1,
                         lnLinf=1,
                         k_M=1,
                         lnLinf_M=1,
                         sd_growth=2,
                         log_p11=-1,
                         log_p22=-1,
                         log_p21=-1,
                         log_p33=-1)

IDs=(1:length(Par.phases))[-match(c("Base case","S1","S2"),names(Par.phases))]
for(i in IDs) Par.phases[[i]]=Par.phases$'Base case'


#Turn off irrelevant pars according to scenarios
fn.mtch=function(WHAT,NMS) match(WHAT,names(NMS))
Q_phz=c("lnq","lnq2","log_Qdaily")                          
# Q_zn_phz=c("lnq_west","lnq_zn1","lnq_zn2",
#            #"lnq2_west","lnq2_zn1","lnq2_zn2",
#            "log_Qdaily_west","log_Qdaily_zn1","log_Qdaily_zn2")
#Zns.par.phz=c("lnR_prop_west","lnR_prop_zn1",Q_zn_phz)
Zns.par.phz=c("lnR_prop_west","lnR_prop_zn1")
MOv.par.phz=c("log_p11","log_p22","log_p21","log_p33")


for(i in 1:N.Scens)     
{
  if(Tabla.scen$Model_type[i]=="Length-based")
  {
    #single zone
    if(Tabla.scen$Spatial_structure[i]=="Single zone")
    {
      Par.phases[[i]][fn.mtch(Zns.par.phz,Par.phases[[i]])]=rep(Fz.off,length(Zns.par.phz))
      #Par.phases[[i]][fn.mtch(Q_phz,Par.phases[[i]])]=c(1,Fz.off,1)
    }
    #effective cpue
    if(Tabla.scen$CPUE[i]=="Effective")
    {
      Par.phases[[i]][fn.mtch(Q_phz,Par.phases[[i]])]=c(1,Fz.off,Fz.off)
      ##Par.phases[[i]][fn.mtch(Q_zn_phz,Par.phases[[i]])]=rep(-5,length(Q_zn_phz))      
      Par.phases[[i]][fn.mtch("log_tau",Par.phases[[i]])]=Fz.off
    }
    #no Qs if no cpue used
    if(Tabla.scen$CPUE[i]%in%c("N/A","No"))
    {
      #Par.phases[[i]][fn.mtch(c(Q_phz,Q_zn_phz),Par.phases[[i]])]=rep(-1,length(c(Q_phz,Q_zn_phz)))    
      Par.phases[[i]][fn.mtch(c(Q_phz),Par.phases[[i]])]=rep(Fz.off,length(c(Q_phz)))
      Par.phases[[i]][fn.mtch("lnR_zero",Par.phases[[i]])]=1
    }
    #no movement
    if(Tabla.scen$Movement[i]%in%c("No","N/A"))
    {
      Par.phases[[i]][fn.mtch(MOv.par.phz,Par.phases[[i]])]=rep(Fz.off,length(MOv.par.phz))  
    }
  }
}
if(Tabla.scen$Model_type[match("S4",Tabla.scen$Model)]=="Age-structured" & Tabla.scen$CPUE[match("S4",Tabla.scen$Model)]=="Stand.")
{
  Par.phases$S4=Par.phases$S2
  Par.phases$S4[match("qdaily",names(Par.phases$S4))]=1
}



#A.4. Modelling arguments

  #What scenarios to run
if(First.run=="YES") Run.all.Scenarios="YES"
if(First.run=="NO") Run.all.Scenarios="NO"  #base case only

Arguments=""  #run without arguments
#Arguments=" -est -sim "  #run with these arguments

Do.zero.Ktch="NO" #do future projections with 0 catch?
if(Do.zero.Ktch=="YES")
{
  zero.Ktch.yrs=50
  zero.Ktch.this=c("S2","S4")
}

Do.sim.test="NO" #do simulation testing?
if(Do.sim.test=="YES")
{
  N.sim.test=50
  CPUE.jitr=1000  
  Sim.Test.this=c("S2","S4")
}

  #MCMC 
DO.MCMC="NO"
nSims=1e6
Thin=10   
burning=1:(5*length(seq(1,nSims,by=Thin))/100)   #5%  burning


#Future projection scenarios
if(First.run=="NO") Run.future.proj="YES"
Future.ktch.scen=list(mean=1,no_ktch=0,percent_50=.5,percent_150=1.5)

#A.5. Output display

#What years to show in model outputs
Show.yrs="DATA"  
#Show.yrs="FUTURE"  #show data plus future projections

#present cpue in log space or normal space
Present.in.log="NO"   

#Catch-MSY arguments   
Growth.F=data.frame(k=0.123,FL_inf=202)
TEMP=18
#Max.age.F=c(20,20)
Max.age.F=c(16,20)
Age.50.mat=c(4,6)
Fecundity=c(3,31)
Breed.cycle=1  #years
Proc.err=0.05 #sigR is PROCESS ERROR; 0 if deterministic model; 0.2 is too high for a shark
SIMS=1e5
years.futures=5

ktch_msy_scen=list()
ktch_msy_scen$'Base case'=list(
  r.prior="USER",  
  user="Yes",
  k.max=50,
  startbio=c(0.8,.95),
  finalbio=c(0.2, 0.7),
  res="low",
  niter=SIMS,
  sigR=Proc.err
)

ktch_msy_scen$S1=list(
  r.prior="USER",  
  user="Yes",
  k.max=50,
  startbio=c(0.8,.95),
  finalbio=c(0.2, 0.7),
  res="low",
  niter=SIMS,
  sigR=0.02 
)

ktch_msy_scen$S2=list(
  r.prior=NA,  
  user="No",
  k.max=50,
  startbio=c(0.8,.95),
  finalbio=c(0.2, 0.7),
  res="low",
  niter=SIMS,
  sigR=Proc.err 
)

#-- B. Execute 'Run.models.R'
#note: this source code brings data and parameter inputs, runs assessment models and displays outputs

source("C:/Matias/Analyses/SOURCE_SCRIPTS/Population dynamics/Run.models.R")





#Previous stuff
# Q.scen=c(rep("two",2),rep("one",2),rep("two",8))
# 
# Tabla.scen=data.frame(
#   Model=Models,
#   Size_comp.=c('Yes',"N/A",'No',rep('Yes',9)),
#   CPUE=c(rep("Stand.",2),"Effective","Effective",rep('Stand.',8)),
#   Age.Growth=c('Yes',"N/A",'No',rep('Yes',9)),
#   Ktch.sx.r=c('Observed','N/A','Equal','Observed','Equal',rep('Observed',7)),                      
#   Tagging=c('No','N/A',rep('No',9),'Yes'),                     
#   Fec.=c(rep('N/A',2),rep('constant',1),rep('N/A',9)),
#   Maturity=c('at length','N/A',rep('knife edge',1),rep("at length",9)),
#   M=c("constant","N/A",rep("constant",10)),
#   M.value=c(M_val,NA,rep(M_val,3),M_val.low,M_val.high,rep(M_val,5)),                      
#   SteepnesS=c(h.M.constant,rep("N/A",2),rep(h.M.constant,6),h.M.constant.low,h.M.constant.up,h.M.constant),
#   Q=Q.scen,   
#   Spatial_structure=c(rep('Single zone',11),'Three zones'),
#   Movement=c("No",rep("N/A",2),rep("No",8),"Yes"),
#   Fo=c(Fo,rep("N/A",2),rep(Fo,4),Fo_Simp,Fo_M,rep(Fo,3)),
#   Model_type=c('Length-based','Biomass dynamics',"Age-structured",rep("Length-based",9))
# )