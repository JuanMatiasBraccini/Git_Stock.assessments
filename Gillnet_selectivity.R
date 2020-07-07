# Script for estimating gear selectivity based on Kirkwood & Walker 1986, extended to 
#       Millar & Holst 1997, and Millar & Fryer 1999

# size classes for selectivity analysis: 5 cm bins (as used in pop dyn model)
# size type for selectivity analysis: total length (in mm)

#assumptions: different mesh sizes set at the same time in same place

#note: are different mesh sizes  used equally?

#MISSING: revise input data (numbers at size class, species names, etc). 
#         Remove commercial data as different mesh sizes not fished at same time??? check with Alex
#         Something out of wack between Figure 1 and Combined selectivity figures (e.g. PortJackson)
#         Compare published selectivities with those estimated here
#         Issues with some species rel. sel. not going to 1 (grey nurse, spinner, etc)
#         Tiger size frequency seems wider than selectivites

#         Implemente Millar method. 

library(tidyverse)
library(RODBC)
library(doParallel)
library(Hmisc)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240) 
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R")
source("C:/Matias/Analyses/Population dynamics/Git_Stock.assessments/NextGeneration.R")
source("C:/Matias/Analyses/Population dynamics/Git_Stock.assessments/SelnCurveDefinitions.R") #These can be extended by the user

Min.sample=25
Min.nets=2

# DATA  -------------------------------------------------------------------

#1. WA Fisheries experimental mesh selectivity studies 

  #1.1 1994-1996
channel <- odbcConnectExcel2007("U:/Shark/ExperimentalNet.mdb")
EXP_NET<- sqlFetch(channel,"EXP_NET", colnames = F)
EXPNET_B<- sqlFetch(channel,"EXPNET_B", colnames = F)
close(channel)

  #1.2 2001-2003
channel <- odbcConnectAccess2007("U:/Shark/Sharks v20200323.mdb")  
Boat_bio=sqlFetch(channel, "Boat_bio", colnames = F) 
Boat_hdr=sqlFetch(channel, "Boat_hdr", colnames = F)   
close(channel)


#2. SESSF 2007-2008 experimental mesh selectivity and survey 
channel <- odbcConnectExcel2007("C:/Matias/Data/SSF_survey_07_08/SharkSurveyData_30_09_2008.xls")
F2_Sampling<- sqlFetch(channel,"F2_Sampling", colnames = F)
F1_SamplingTwo<- sqlFetch(channel,"F1_SamplingTwo", colnames = F)
close(channel)


#3. TDGDLF observed catch composition
LFQ.south=read.csv("C:/Matias/Analyses/Selectivity/out.LFQ.south.csv")



#Species names
SP.names=read.csv('C:/Matias/Data/Species_names_shark.only.csv')
SP.codes=read.csv('C:/Matias/Data/Species.code.csv')


# PROCEDURE  -------------------------------------------------------------------

#1. Data manipulation

  #1.1. Rory's

  #1.1.1 1994-1996
EXP_NET=EXP_NET[grep("E",EXP_NET$SHEET_NO),]%>%
           mutate(BOAT=VESSEL)%>%
          mutate(SKIPPER=NA,
                 Method='GN',
                 BOTDEPTH=NA,
                 SOAK.TIME=NA,
                 NET_LENGTH=NA,
                 Lng1=END1_LNG_D+(END1_LNG_M/60),
                 Lng2=END2_LNG_D+(END2_LNG_M/60),
                 Lat1=END1_LAT_D+(END1_LAT_M/60),
                 Lat2=END2_LAT_D+(END2_LAT_M/60))
EXP_NET=EXP_NET%>%
          mutate(MID.LONG=rowMeans(select(EXP_NET,starts_with("Lng"))),
                 MID.LAT=-rowMeans(select(EXP_NET,starts_with("Lat"))))%>%
  dplyr::select(SHEET_NO,DATE,BOAT,SKIPPER,Method,START_SET,END_HAUL,BOTDEPTH,
                MID.LAT,MID.LONG,SOAK.TIME,NET_LENGTH)
    
EXPNET_B=EXPNET_B[grep("E",EXPNET_B$SHEET_NO),]%>%
          left_join(EXP_NET,by='SHEET_NO')%>%
          filter(SPP_CODE<50000)%>%
          mutate(species=SPP_CODE)

Exp.net.94_96=EXPNET_B%>%
                left_join(SP.names,by=c("SPP_CODE" = "SPECIES"))%>%
                filter(!is.na(Name))%>%
                rename(tl=TOT_LENGTH,
                       fl=FORK_LNGTH)
colnames(Exp.net.94_96)=tolower(colnames(Exp.net.94_96))



  #1.1.2 2001-2003
Boat_hdr=Boat_hdr[grep("E",Boat_hdr$SHEET_NO),]%>%
          dplyr::select(SHEET_NO,DATE,BOAT,SKIPPER,Method,START_SET,END_HAUL,BOTDEPTH,
                        'MID LAT','MID LONG','SOAK TIME',MESH_SIZE,MESH_DROP,NET_LENGTH)%>%
          rename(MID.LAT='MID LAT',MID.LONG='MID LONG',SOAK.TIME='SOAK TIME')%>%
          filter(Method=='GN')
Boat_bio=Boat_bio[grep("E",Boat_bio$SHEET_NO),]%>%
          dplyr::select(SHEET_NO,SPECIES,TL,FL,PL,SEX)
Exp.net.01_03=left_join(Boat_bio,Boat_hdr,by="SHEET_NO")%>%
          mutate(SEX=ifelse(SEX=="m","M",ifelse(SEX=="f","F",SEX)),
                 MESH_SIZE=ifelse(MESH_SIZE=="10\"","10",
                           ifelse(MESH_SIZE=="6\"","6",
                           ifelse(MESH_SIZE=="5\r\n5","5",
                           ifelse(MESH_SIZE=="7\"","7",
                           ifelse(MESH_SIZE=="5\"","5",
                           ifelse(MESH_SIZE=="4\"","4",
                           ifelse(MESH_SIZE=="8\"","8",
                           MESH_SIZE))))))),
                 MESH_SIZE=as.numeric(MESH_SIZE))%>%
          dplyr::select(-PL)%>%
  left_join(SP.codes,by=c("SPECIES" = "Species"))%>%
  rename(Name=COMMON_NAME)
colnames(Exp.net.01_03)=tolower(colnames(Exp.net.01_03))

Exp.net.94_96$experiment='94_96'
Exp.net.01_03$experiment='01_03'

This.col=c('sheet_no','date','experiment','mid.lat','mid.long','mesh_size','mesh_drop','name','tl','fl','sex')

Exp.net.WA=rbind(Exp.net.94_96[,match(This.col,names(Exp.net.94_96))],
                 Exp.net.01_03[,match(This.col,names(Exp.net.01_03))])%>%
          mutate(name=ifelse(name=="Angel Shark (general)","Angel Shark",
                      ifelse(name=="Eagle ray","Eagle Ray",
                      ifelse(name=="Gummy shark","Gummy Shark",
                      ifelse(name=="Port Jackson","PortJackson shark",
                      ifelse(name=="Big eye sixgill shark","Sixgill shark",
                      ifelse(name%in%c("Wobbegong (general)",'Spotted Wobbegong','Western Wobbegong'),"Wobbegongs",name)))))))

Published=capitalize(tolower(c("Sandbar shark","Gummy Shark","Whiskery shark","Dusky shark")))

TAB=table(Exp.net.WA$name,Exp.net.WA$mesh_size)
TAB[TAB<2]=0
TAB[TAB>=2]=1

Exp.net.WA=Exp.net.WA%>%
  filter(name%in%names(which(rowSums(TAB)>=2)))%>%
  mutate(tl=ifelse(tl>500,tl/10,tl))

Preliminary=FALSE
if(Preliminary) ggplot(Exp.net.WA,aes(tl,fl,shape=name, colour=name, fill=name))+
                geom_point() + 
                geom_smooth(method = "lm", fill = NA)+
                facet_wrap(vars(name), scales = "free")

#Conversion TL-FL
TL_FL=data.frame(name=c('Angel Shark','Dusky shark','Gummy Shark','Pencil shark',
           'PortJackson shark','Sandbar shark','Smooth hammerhead','Spurdogs',
           'Whiskery shark','Tiger shark'),intercept=NA,slope=NA)
for(l in 1:nrow(TL_FL))
{
  a=Exp.net.WA%>%filter(name==TL_FL$name[l])
  mod=lm(tl~fl,a)
  COEF=coef(mod)
  TL_FL$intercept[l]=COEF[1]
  TL_FL$slope[l]=COEF[2]
}
TL_FL=rbind(TL_FL,
            data.frame(name='Sliteye shark',intercept=7.0195,slope=1.134)) # add sliteye conversion manually (Gutridge et al 2011)

#Derive TL from FL if no TL observation
Exp.net.WA=Exp.net.WA%>%
           left_join(TL_FL,by='name')%>%
            mutate(tl=ifelse(is.na(tl),intercept+fl*slope,tl),
                   Length=tl*10)%>%   #Length in mm; use Total length
            filter(!is.na(Length))%>%
            filter(!name=='Sliteye shark')#too few observations for sliteye

#Size frequency
ggplot(Exp.net.WA,aes(tl,colour=name, fill=name))+
                geom_histogram(alpha=0.6, binwidth = 5, show.legend = FALSE) +
                  xlab("Total length (cm)")+
                  facet_wrap(vars(name), scales = "free") 
ggsave('C:/Matias/Analyses/Selectivity/Size.frequency_experimental.WA.tiff', width = 8,height = 8, dpi = 300, compression = "lzw")


  #1.2 SSF
F2_Sampling=F2_Sampling%>%
                filter(Csiro>37000002 & Csiro<37039000 & LengthType=="TL" &
                      !is.na(Mesh1) & !is.na(Length))%>%
                mutate(Length=Length/10,              
                       Length=ifelse(Length>350,Length/10,Length),
                       Mesh1=ifelse(Mesh1%in%c('C6.00', 'C6.01', 'C6.02', 'C6.03', 'C6.04'),'C6',
                                    ifelse(Mesh1%in%c('C6.51', 'C6.52', 'C6.53', 'C6.54'),'C6.5',Mesh1)),
                       Mesh.size=2.54*as.numeric(substr(Mesh1,2,10)))%>%   #inches to cm
                filter(Length>=30)%>%
                mutate(Length=Length*10)  #length in mm

F2_Sampling=F2_Sampling%>%
              dplyr::select(Species,Csiro,Sex,Mesh1,Mesh.size,Length)

if(Preliminary) ggplot(F2_Sampling, aes(x = Length/10)) +
                    geom_histogram(color = "grey30", fill ="salmon",binwidth=10) +
                    facet_grid(Species~Mesh.size, scales = "free")

#Size frequency
ggplot(F2_Sampling,aes(Length/10,colour=Species, fill=Species))+
  geom_histogram(alpha=0.6, binwidth = 5, show.legend = FALSE) +
  xlab("Total length (cm)") +
  facet_wrap(vars(Species), scales = "free") 
ggsave('C:/Matias/Analyses/Selectivity/Size.frequency_SSF.tiff', width = 12,height = 8, dpi = 300, compression = "lzw")



  #1.3 TDGLDF
LFQ.south=LFQ.south%>%
          mutate(Data.set="Obs.TDGDLF",
                 Mesh.size=2.54*MESH_SIZE,
                 Species=capitalize(COMMON_NAME),
                 Species=ifelse(Species=="Sawsharks","Common sawshark",
                         ifelse(Species=="Wobbegong (general)","Wobbegong",
                         Species)))%>%
          filter(!Species=="Wobbegong")  #remove Wobbies because size measurements are uncertain

#Convert FL to TL
TL_FL_LFQ.south=TL_FL%>%
                filter(name%in%c("Smooth hammerhead","Spurdogs","Tiger shark"))

add1=TL_FL_LFQ.south%>%filter(name=="Smooth hammerhead")
add1=rbind(add1,add1)%>%mutate(name=c("Great hammerhead","Scalloped hammerhead"))  

TL_FL_LFQ.south=rbind(TL_FL_LFQ.south,
                      add1,
                      data.frame(name="Copper shark",intercept=5.801,slope=1.181),
                      data.frame(name="Common sawshark",intercept=0,slope=1.2),
                      data.frame(name="Grey nurse shark",intercept=0,slope=1.2),
                      data.frame(name="Milk shark",intercept=0,slope=1.2),
                      data.frame(name="Shortfin mako",intercept=0,slope=1.127),
                      data.frame(name="Spinner shark",intercept=0,slope= 1.28))  

LFQ.south=LFQ.south%>%
          left_join(TL_FL_LFQ.south,by=c(Species='name'))%>%
          mutate(tl=intercept+FL*slope,
                 Length=tl*10)%>%   #Length in mm; use Total length
          filter(!is.na(Length))%>%
    dplyr::select(Species,Mesh.size,Length,Data.set)

#Size frequency
ggplot(LFQ.south,aes(Length/10,colour=Species, fill=Species))+
  geom_histogram(alpha=0.6, binwidth = 5, show.legend = FALSE) +
  xlab("Total length (cm)") +
  facet_wrap(vars(Species), scales = "free") 
ggsave('C:/Matias/Analyses/Selectivity/Size.frequency_TDGDLF_observed.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")



#2. Combine all data sets
Exp.net.WA=Exp.net.WA%>%
  rename(Species=name)%>%
  mutate(Mesh.size=2.54*mesh_size,
         Data.set="WA")
F2_Sampling=F2_Sampling%>%
        mutate(Data.set="SSF",
               Species=ifelse(Species=='Bronze Whaler',"Copper shark",Species))

####remove this when U drive works
Combined=F2_Sampling%>%dplyr::select(Species,Mesh.size,Length,Data.set)%>%
  mutate(Species=tolower(Species),
         Species=ifelse(Species=='portjackson shark','port jackson shark',
                        ifelse(Species=='wobbegong','wobbegongs',
                               Species)),
         Species=capitalize(Species))

#####

###switch back when U drive works
# Combined=rbind(F2_Sampling%>%dplyr::select(Species,Mesh.size,Length,Data.set),
#                Exp.net.WA%>%dplyr::select(Species,Mesh.size,Length,Data.set))%>%
#                   mutate(Species=tolower(Species),
#                          Species=ifelse(Species=='portjackson shark','port jackson shark',
#                                  ifelse(Species=='wobbegong','wobbegongs',
#                                   Species)),
#                          Species=capitalize(Species))
############
Combined=rbind(Combined,LFQ.south)%>%
  mutate(Species=ifelse(Species=='Angel shark','Australian angelshark',Species))%>%
  filter(!is.na(Mesh.size))


#Analyse species with at least Min.sample and Min.nets
TAB=table(Combined$Species,Combined$Mesh.size)
TAB[TAB<Min.sample]=0
TAB[TAB>=Min.sample]=1
Combined=Combined%>%
  filter(Species%in%names(which(rowSums(TAB)>=Min.nets)) & Length<=3500 & Mesh.size<22)
  

#for each selected species, remove meshes with few observations
min.obs.per.mesh=Min.sample
n.sp=sort(unique(Combined$Species))
for(s in 1:length(n.sp))
{
  d=Combined%>%filter(Species==n.sp[s])
  Combined=Combined%>%filter(!Species==n.sp[s])
  
  id=table(d$Mesh.size)
  Drop=as.numeric(gsub("[^0-9.]", "",  names(which(id<min.obs.per.mesh))))
  if(length(Drop)>0)
  {
    d=d%>%filter(!Mesh.size%in%Drop)
    id=table(d$Mesh.size)
    if(length(id)==1) d=NULL
  }
  
  Combined=rbind(Combined,d)
  
}
n.sp=sort(unique(Combined$Species))


#3. Estimate selectivity parameters 


  #3.1 Millar & Holst 1997 
Fitfunction='gillnetfit'
#Fitfunction='NetFit'
  
Rtype=c("norm.loc","norm.sca","gamma","lognorm")

Millar.Holst=function(d,size.int)
{
  #Create size bins
  d=d%>%mutate(Size.class=size.int*floor(Length/size.int)+size.int/2)
  
  #Tabulate observations by mid size class and mesh size   
  tab=d%>%
    group_by(Mesh.size,Size.class)%>%
    summarise(n=n())%>%
    spread(Mesh.size,n,fill=0)%>%
    data.frame
  Meshsize=as.numeric(substr(names(tab)[-1],2,10))
  
  #Fit SELECT model  ACA. Gamma not implemented in NetFit, can I use gamma estim pars from gillnetfit?
  
  #Equal fishing power
  pwr=rep(1,length(Meshsize))
  Equal.power=vector('list',length(Rtype))
  names(Equal.power)=Rtype
  
  #gillnetfit approah
  if(Fitfunction=='gillnetfit')
  {
    for(f in 1:length(Equal.power))Equal.power[[f]]=gillnetfit(data=as.matrix(tab),
                                                               meshsizes=Meshsize,
                                                               type=Rtype[f],
                                                               plots=c(T,T),
                                                               plotlens=NULL,
                                                               details=F)
  }
  
  #NetFit approach (gamma not implemented)
  if(Fitfunction=='NetFit')
  {
    Init.par=c(mean(d$Length),sd(d$Length))
    for(f in 1:length(Equal.power))Equal.power[[f]]=NetFit(Data=tab,Meshsize=Meshsize,
                                                               x0=Init.par,rtype=Rtype[f],
                                                               rel.power=pwr)
  }

  
}
for(s in 1:length(n.sp))
{
  Millar.Holst(d=Combined%>%filter(Species==n.sp[s]),size.int=50)
}


  #3.2 Kirkwood & Walker
Selectivty.Kirkwood.Walker=function(d,size.int,theta)
{
  #Create size bins
  d=d%>%mutate(Size.class=size.int*floor(Length/size.int)+size.int/2)
  
  #Tabulate observations by mid size class and mesh size
  tab=d%>%
    group_by(Mesh.size,Size.class)%>%
    summarise(n=n())%>%
    spread(Mesh.size,n,fill=0)%>%
    data.frame
  row.names(tab)=tab$Size.class
  tab=tab[,-1]
  
  
  #Calculate relative selectivity
  Theta1=exp(theta[1])
  Theta2=exp(theta[2])
  d=d%>%
    mutate(alpha.beta=Theta1*Mesh.size,
           beta=-0.5*((alpha.beta)-((alpha.beta*alpha.beta+4*Theta2)^0.5)),
           alpha=alpha.beta/beta,
           Rel.sel=((Size.class/(alpha*beta))^alpha)*(exp(alpha-(Size.class/beta))))
  
  
  #Log likelihood
  S.ij=d%>%
    distinct(Size.class,Mesh.size,.keep_all = T)%>%
    dplyr::select(Size.class,Mesh.size,Rel.sel)%>%
    spread(Mesh.size,Rel.sel,fill = 0)
  row.names(S.ij)=S.ij$Size.class
  S.ij=S.ij[,-1]
  
  mu.j=rowSums(tab)/rowSums(S.ij)
  mu.j=sapply(mu.j,function(x) max(x,0.1))
  mu.j.prop=mu.j/sum(mu.j)
  
  
  #predicted numbers
  sum.n=colSums(tab)
  NN=ncol(S.ij)
  tab.pred=(S.ij*matrix(rep(mu.j.prop,NN),ncol=NN)*matrix(rep(sum.n,each=nrow(S.ij)),ncol=NN))/
    (matrix(rep(colSums(S.ij*matrix(rep(mu.j.prop,NN))),each=nrow(S.ij)),ncol=NN))
  
  
  #Gamma log like
  negLL=min(-sum(tab*(log(mu.j*S.ij))-(mu.j*S.ij),na.rm=T),1e100)
  
  return(list(negLL=negLL,d=d,observed=tab,predicted=tab.pred))
  
}
theta.list=vector('list',length(n.sp))
names(theta.list)=n.sp
Fit=Combined.sel=theta.list

# initial parameter values
theta.list$`Gummy shark`=c(Theta1=log(80),Theta2=log(29000))
for(s in 1:length(theta.list)) theta.list[[s]]=jitter(theta.list$`Gummy shark`,factor=.1)

# fit model
for(s in 1:length(n.sp))
{
  theta=theta.list[[s]]
  
  #. objfun to minimize
  fn_ob=function(theta)Selectivty.Kirkwood.Walker(d=Combined%>%filter(Species==n.sp[s]),
                                                  size.int=50,
                                                  theta)$negLL
  
  #. fit model
  Fit[[s]]=nlminb(theta.list[[s]], fn_ob, gradient = NULL)
}

# Calculate confidence intervals thru bootstrapping 
n.boot=1:1000
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
system.time({
  Fit.CI=foreach(s=1:length(n.sp),.packages=c('tidyverse','doParallel','Biobase')) %dopar%
    {
      boot=foreach(n=n.boot,.packages=c('doParallel','splitstackshape','tidyverse')) %dopar%
      {
          theta=theta.list[[s]]
          
          #bootstrapped sample
          d.samp=Combined%>%filter(Species==n.sp[s]) 
          d.samp=stratified(d.samp, "Mesh.size",size=nrow(d.samp),replace=TRUE)
          
          #. objfun to minimize
          fn_ob=function(theta)Selectivty.Kirkwood.Walker(d=d.samp,
                                                          size.int=100,
                                                          theta)$negLL
          
          #. fit model
          return(nlminb(theta.list[[s]], fn_ob, gradient = NULL))
      }
      return(exp(do.call(rbind,subListExtract(boot,"par"))))
    }
})    #takes 0.5 sec per iteration per species
names(Fit.CI)=n.sp
stopCluster(cl)



# REPORT  -------------------------------------------------------------------
setwd('C:/Matias/Analyses/Selectivity')

#Plot each species' selectivity
fn.plt.Sel=function(Dat,theta)
{
  Theta1=exp(theta[1])
  Theta2=exp(theta[2])
  
  sizes=seq(10*round(min(Dat$Length)*.9/10),10*(round(max(Dat$Length)*1.1/10)),by=10)
  meshes=unique(Dat$Mesh.size)
  d=data.frame(Size.class=rep(sizes,length(meshes)))%>%
    mutate(Mesh.size=rep(meshes,each=length(sizes)),
           alpha.beta=Theta1*Mesh.size,
           beta=-0.5*((alpha.beta)-((alpha.beta*alpha.beta+4*Theta2)^0.5)),
           alpha=alpha.beta/beta,
           Rel.sel=((Size.class/(alpha*beta))^alpha)*(exp(alpha-(Size.class/beta))))
  
  S.ij=d%>%
    distinct(Size.class,Mesh.size,.keep_all = T)%>%
    dplyr::select(Size.class,Mesh.size,Rel.sel)%>%
    spread(Mesh.size,Rel.sel,fill = 0)
  row.names(S.ij)=S.ij$Size.class
  S.ij=S.ij[,-1]
  
  ggplot(d, aes(Size.class,  Rel.sel)) + geom_line(aes(colour = factor(Mesh.size)),size=1.5)
  ggsave(paste("Each species/Selectivity/",n.sp[s],'.tiff',sep=''), width = 8,height = 8, dpi = 300, compression = "lzw")
  
  return(S.ij)
}
for(s in 1:length(n.sp)) Combined.sel[[s]]=fn.plt.Sel(Dat=Combined%>%filter(Species==n.sp[s]),theta=Fit[[s]]$par)

#Output parameter estimates            
Table1=data.frame(Species=n.sp,
                  Theta1=NA,Theta1.LOW95=NA,Theta1.UP95=NA,
                  Theta2=NA,Theta2.LOW95=NA,Theta2.UP95=NA)
for(s in 1:length(n.sp))
{
  dummy=Fit.CI[[s]]
  Table1$Theta1[s]=round(quantile(dummy[,"Theta1"],probs=0.5),1)
  Table1$Theta1.LOW95[s]=round(quantile(dummy[,"Theta1"],probs=0.025),1)
  Table1$Theta1.UP95[s]=round(quantile(dummy[,"Theta1"],probs=0.975),1)
  Table1$Theta2[s]=round(quantile(dummy[,"Theta2"],probs=0.5),1)
  Table1$Theta2.LOW95[s]=round(quantile(dummy[,"Theta2"],probs=0.025),1)
  Table1$Theta2.UP95[s]=round(quantile(dummy[,"Theta2"],probs=0.975),1)
}
write.csv(Table1,'Table1.csv',row.names = F) 

Table11=Table1%>%
        mutate(Theta1=paste(Theta1," (",paste(Theta1.LOW95,Theta1.UP95,sep="-"),")",sep=''),
               Theta2=paste(Theta2," (",paste(Theta2.LOW95,Theta2.UP95,sep="-"),")",sep=''))%>%
  dplyr::select(Species,Theta1,Theta2)
#colnames(Table11)[2]=intToUtf8(0x03B8)
#colnames(Table11)[2]=intToUtf8(U+03B8)
fn.word.table(WD=getwd(),TBL=Table11,Doc.nm="Table1",caption=NA,paragph=NA,
              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
              Zebra='NO',Zebra.col='grey60',Grid.col='black',
              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")



#Combined selectivity
tiff(file="Combined selectivity all selected species.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
smart.par(n.plots=length(n.sp),MAR=c(1.5,1,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
for(s in 1:length(n.sp))
{
  Sum.sel=rowSums(Combined.sel[[s]])
  Combined.sel[[s]]$combined=Sum.sel/max(Sum.sel)
  plot(as.numeric(rownames(Combined.sel[[s]])),Combined.sel[[s]]$combined,type='l',ylab='',xlab='',
       lwd=5,col='orange',main=n.sp[s])
  for(l in 1:(ncol(Combined.sel[[s]])-1)) lines(as.numeric(rownames(Combined.sel[[s]])),Combined.sel[[s]][,l])
}
plot.new()
legend('center',"combined",lwd=4,col="orange",bty='n',cex=1.25)
mtext("Total length (mm)",1,outer=T,line=.35)
mtext("Relative selectivity",2,outer=T,las=3,line=1.75)
dev.off()  


#Display size frequencies
d=Combined%>%mutate(MESH=as.factor(Mesh.size/2.54),
                    TL=Length/10)
d %>%
  ggplot( aes(x=TL, fill=MESH)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 10)  +
  labs(fill="Mesh size")+
  facet_wrap(vars(Species), scales = "free_y")+
  xlab("Total length (cm)")+ ylab("Frequency")
ggsave('Size frequency all selected species_data set combined.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")

#Plot fit (observed vs predicted number at size by mesh)
for(s in 1:length(n.sp))
{
  dummy=Selectivty.Kirkwood.Walker(d=Combined%>%filter(Species==n.sp[s]),
                                   size.int=100,
                                   Fit[[s]]$par)
  tiff(file=paste("Each species/Fit/Pre_vs_Obs_",n.sp[s],".tiff",sep=''),width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")    
  smart.par(n.plots=ncol(dummy$observed),MAR=c(2,1,2,1.5),OMA=c(2.5,3,.5,.1),MGP=c(1,.5,0))
  for(m in 1:ncol(dummy$observed))
  {
    DAT=data.frame(y=dummy$observed[,m],x=dummy$predicted[,m])
    mod=lm(y~x,data=DAT)
    CoF=coef(mod)
    Smry=summary(mod)
    Ndat=data.frame(x=seq(min(DAT$x),max(DAT$x),length.out = 10))
    Prd=predict(mod,newdata=Ndat,se.fit=TRUE)
    Main=paste(as.numeric(gsub("[^0-9.]", "",colnames(dummy$observed)[m]))/2.54,'inch')
    plot(dummy$predicted[,m],dummy$observed[,m],pch=19,cex=1.5,ylab='',xlab='',main=Main)
    lines(Ndat$x,Prd$fit,lwd=1.5)
    text(mean(Ndat$x),mean(Prd$fit)*.6,
         paste('y= ',round(CoF[1],2),' + ',round(CoF[2],2),'x',sep=''),pos=4,cex=1.25)
    mylabel = bquote(italic(r)^2 == .(format(round(Smry$r.squared,2), digits = 3)))
    text(mean(Ndat$x),mean(Prd$fit)*.4,mylabel,pos=4,cex=1.25)
    polygon(x=c(Ndat$x,rev(Ndat$x)),
            y=c(Prd$fit+1.96*Prd$se.fit,rev(Prd$fit-1.96*Prd$se.fit)),
            col=rgb(.1,.1,.1,alpha=.2),border='transparent')
  }
  mtext("Predicted size class (mm)",1,outer=T,line=1)
  mtext("Observed size class (mm)",2,outer=T,las=3,line=1.5)
  dev.off()  
}

#Select species without selectivity published
n.sp=n.sp[-match(Published,n.sp)]


#Selectivity
colfunc <- colorRampPalette(c("cadetblue2", "deepskyblue4"))
unik.mesh=sort(unique(Combined$Mesh.size))
CLS=colfunc(length(unik.mesh))
names(CLS)=unik.mesh
tiff(file="Figure 1.Selectivity.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
smart.par(n.plots=length(n.sp),MAR=c(1.5,1.2,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
for(s in 1:length(n.sp))
{
  d=Combined.sel[[s]]%>%
    mutate(Length=as.numeric(rownames(Combined.sel[[s]])),
           TL=Length/10)
  
  plot(d$TL,d$TL,type='l',ylab='',xlab='',
       lwd=5,col='transparent',main=n.sp[s],ylim=c(0,1))
  for(l in 1:(ncol(d)))
  {
    lines(d$TL,d[,l],col=CLS[match(names(d)[l],names(CLS))],lwd=2)
  }
}
plot.new()
legend("center",paste(as.numeric(names(CLS))/2.54,'inch'),col=CLS,bty='n',lwd=2,cex=1.25)
mtext("Total length (cm)",1,outer=T,line=.35,cex=1.25)
mtext("Relative selectivity",2,outer=T,las=3,line=1.5,cex=1.25)
dev.off() 


#Display size frequencies
d=Combined%>%mutate(MESH=as.factor(Mesh.size/2.54),
                    TL=Length/10)%>%
  filter(Species%in%n.sp)
d %>%
  ggplot( aes(x=TL, fill=MESH)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', binwidth = 10)  +
  labs(fill="Mesh size")+
  facet_wrap(vars(Species), scales = "free_y")+
  xlab("Total length (cm)")+ ylab("Frequency")
ggsave('Figure 2.Size frequency.tiff', width = 10,height = 8, dpi = 300, compression = "lzw")

