# Script for estimating gear selectivity based on Kirkwood & Walker 1986

#assumptions: different mesh sizes set at the same time in same place

#note: are different mesh sizes  used equally?
#      fish size in mm



library(tidyverse)
library(RODBC)
library(doParallel)
library(Hmisc)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240) 
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))

Min.sample=10

# DATA  -------------------------------------------------------------------

#1. Rory's gillnet selectivity study 

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


#2. SSF Shark survey 2007-2008
channel <- odbcConnectExcel2007("H:/Backups/Matias_4_13/Data/Shark survey/2007-2008/SharkSurveyData_30_09_2008.xls")
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

# add sliteye conversion manually (Gutridge et al 2011)
TL_FL=rbind(TL_FL,
            data.frame(name='Sliteye shark',intercept=7.0195,slope=1.134))

Exp.net.WA=Exp.net.WA%>%
           left_join(TL_FL,by='name')%>%
            mutate(tl=ifelse(is.na(tl),intercept+fl*slope,tl),
                   Length=tl*10)%>%   #Length in mm; use Total length
            filter(!is.na(Length))%>%
            filter(!name=='Sliteye shark')#too few observations for sliteye


if(Preliminary) ggplot(Exp.net.WA,aes(tl,colour=name, fill=name))+
                geom_histogram(alpha=0.6, binwidth = 5) +
                  facet_wrap(vars(name), scales = "free")


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


  #1.3 TDGLDF
LFQ.south=LFQ.south%>%
          mutate(Data.set="Obs.TDGDLF",
                 Mesh.size=2.54*MESH_SIZE,
                 Species=capitalize(COMMON_NAME),
                 Species=ifelse(Species=="Sawsharks","Common sawshark",
                         ifelse(Species=="Wobbegong (general)","Wobbegong",
                         Species)))%>%
          filter(!Species=="Wobbegong")  #remove Wobbies because measurement ucertain
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



#2. Combine SSF, Rory's experimental and TDGDLF observed
Exp.net.WA=Exp.net.WA%>%
  rename(Species=name)%>%
  mutate(Mesh.size=2.54*mesh_size,
         Data.set="WA")
F2_Sampling=F2_Sampling%>%
        mutate(Data.set="SSF",
               Species=ifelse(Species=='Bronze Whaler',"Copper shark",Species))

Combined=rbind(F2_Sampling%>%dplyr::select(Species,Mesh.size,Length,Data.set),
               Exp.net.WA%>%dplyr::select(Species,Mesh.size,Length,Data.set))%>%
                  mutate(Species=tolower(Species),
                         Species=ifelse(Species=='portjackson shark','port jackson shark',
                                 ifelse(Species=='wobbegong','wobbegongs',
                                  Species)),
                         Species=capitalize(Species))

Combined=rbind(Combined,LFQ.south)%>%
  mutate(Species=ifelse(Species=='Angel shark','Australian angelshark',Species))

TAB=table(Combined$Species,Combined$Mesh.size)
TAB[TAB<Min.sample]=0
TAB[TAB>=Min.sample]=1
Combined=Combined%>%
  filter(Species%in%names(which(rowSums(TAB)>=2)))
  


#3. Estimate selectivity parameters 
min.obs.per.mesh=20
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

n.sp=unique(Combined$Species)

#remove meshes with few observations
Chek=vector('list',length(n.sp))
names(Chek)=n.sp
size.int=100
for(s in 1:length(n.sp))
{
  #Create size bins
  d=Combined%>%filter(Species==n.sp[s])%>%mutate(Size.class=size.int*floor(Length/size.int)+size.int/2)
  
  #Tabulate observations by mid size class and mesh size
  tab=d%>%
    group_by(Mesh.size,Size.class)%>%
    summarise(n=n())%>%
    spread(Mesh.size,n,fill=0)%>%
    data.frame
  row.names(tab)=tab$Size.class
  tab=tab[,-1]
  
  
  id=colSums(tab)
  Drop=as.numeric(gsub("[^0-9.]", "",  names(which(id<min.obs.per.mesh))))
  if(length(Drop)>0)
  {
    d=d%>%filter(!Mesh.size%in%Drop)
    tab=d%>%
      group_by(Mesh.size,Size.class)%>%
      summarise(n=n())%>%
      spread(Mesh.size,n,fill=0)%>%
      data.frame
    row.names(tab)=tab$Size.class
    tab=tab[,-1]
  } 
  if(! is.numeric(tab))
  {
    if(nrow(tab)>0) AA=n.sp[s]
    if(nrow(tab)==0) AA=NA
  }
  if(is.numeric(tab))AA=NA
  n.sp[s]=AA
  rm(AA)
}

n.sp=n.sp[!is.na(n.sp)]

theta.list=vector('list',length(n.sp))
names(theta.list)=n.sp
Fit=Combined.sel=theta.list

#initial parameter values
theta.list$`Gummy shark`=c(Theta1=log(80),Theta2=log(29000))
for(s in 1:length(theta.list)) theta.list[[s]]=jitter(theta.list$`Gummy shark`,factor=.1)
#theta=c(Theta1=log(80),Theta2=log(29000))

# fit model
for(s in 1:length(n.sp))
{
  theta=theta.list[[s]]
  
  #. objfun to minimize
  fn_ob=function(theta)Selectivty.Kirkwood.Walker(d=Combined%>%filter(Species==n.sp[s]),
                                                  size.int=100,
                                                  theta)$negLL
  
  #. fit model
  Fit[[s]]=nlminb(theta.list[[s]], fn_ob, gradient = NULL)
}


# 4. Calculate confidence intervals thru bootstrapping 
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

#Plot selectivity
fn.plt.Sel=function(Dat,theta)
{
  Theta1=exp(theta[1])
  Theta2=exp(theta[2])
  
  sizes=seq(round(min(Dat$Length)*.9),round(max(Dat$Length)*1.2))
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
  ggsave(paste("Selectivity_",n.sp[s],'.tiff',sep=''), width = 8,height = 8, dpi = 300, compression = "lzw")
  
  return(S.ij)
}
for(s in 1:length(n.sp)) Combined.sel[[s]]=fn.plt.Sel(Dat=Combined%>%filter(Species==n.sp[s]),theta=Fit[[s]]$par)


#Combined selectivity
tiff(file="Combined selectivity.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
smart.par(n.plots=length(n.sp),MAR=c(1.5,.5,1.5,1.5),OMA=c(1.5,3,.1,.1),MGP=c(1,.5,0))
for(s in 1:length(n.sp))
{
  Sum.sel=rowSums(Combined.sel[[s]])
  Combined.sel[[s]]$combined=Sum.sel/max(Sum.sel)
  plot(as.numeric(rownames(Combined.sel[[s]])),Combined.sel[[s]]$combined,type='l',ylab='',xlab='',
       lwd=5,col='orange',main=n.sp[s])
  for(l in 1:(ncol(Combined.sel[[s]])-1)) lines(as.numeric(rownames(Combined.sel[[s]])),Combined.sel[[s]][,l])
}
mtext("Total length (mm)",1,outer=T,line=.5)
mtext("Relative selectivity",2,outer=T,las=3,line=1.75)
dev.off()  


#Plot predicted numbers
for(s in 1:length(n.sp))
{
  dummy=Selectivty.Kirkwood.Walker(d=Combined%>%filter(Species==n.sp[s]),
                                   size.int=100,
                                   Fit[[s]]$par)
  tiff(file=paste("Pre_vs_Obs_",n.sp[s],".tiff",sep=''),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  
  smart.par(n.plots=ncol(dummy$observed),MAR=c(.1,.1,.1,.1),OMA=c(2.5,2.5,1.5,.1),MGP=c(1,.5,0))
  for(m in 1:ncol(dummy$observed))
  {
    plot(dummy$observed[,m],dummy$predicted[,m],pch=19)
    lines(dummy$observed[,m],dummy$observed[,m],col=2)
  }
  dev.off()  
}
