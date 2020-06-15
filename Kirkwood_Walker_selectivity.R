# Script for estimating gear selectivity based on Kirkwood & Walker 1986

#assumptions: different mesh sizes set at the same time in same place

#note: it is important to know if different mesh sizes are used equally
#      fish size in mm

#Missing: Rory's gillnet selectivity data


library(tidyverse)
library(RODBC)
library(doParallel)

options(stringsAsFactors = FALSE,"max.print"=50000,"width"=240) 
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))

# DATA  -------------------------------------------------------------------

#1. Rory's gillnet selectivity study (need to enter data)           #MISSING


#2. SSF Shark survey 2007-2008
channel <- odbcConnectExcel2007("H:/Backups/Matias_4_13/Data/Shark survey/2007-2008/SharkSurveyData_30_09_2008.xls")
F2_Sampling<- sqlFetch(channel,"F2_Sampling", colnames = F)
F1_SamplingTwo<- sqlFetch(channel,"F1_SamplingTwo", colnames = F)
close(channel)


#3. TDGDLF catch composition
Do.TDGDLF=FALSE
if(Do.TDGDLF)
{
  LFQ.south=LFQ.south%>%mutate(COMMON_NAME=ifelse(COMMON_NAME=="Common sawshark","Sawsharks",
                                                  ifelse(COMMON_NAME=="Wobbegong (general)","Wobbegongs",
                                                         COMMON_NAME)))
  for(k in 1:N.sp)
  {
    dd=LFQ.south%>%
      mutate(Name=tolower(COMMON_NAME))%>%
      filter(Name==Specs$SP.group[k])
  }
}


# PROCEDURE  -------------------------------------------------------------------

#1. Data manipulation

  #1.1. Rory's
  
  #1.2. SSF
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

TAB=table(F2_Sampling$Species,F2_Sampling$Mesh.size)
TAB[TAB<50]=0
TAB[TAB>=50]=1
F2_Sampling=F2_Sampling%>%
              filter(Species%in%names(which(rowSums(TAB)>=2)))%>%
              dplyr::select(Species,Csiro,Sex,Mesh1,Mesh.size,Length)
ggplot(F2_Sampling, aes(x = Length/10)) +
  geom_histogram(color = "grey30", fill ="salmon",binwidth=10) +
  facet_grid(Species~Mesh.size, scales = "free")




#2. Selectivity parameter estimation
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
  tab.pred=(S.ij*matrix(rep(mu.j.prop,3),ncol=3)*matrix(rep(sum.n,each=nrow(S.ij)),ncol=3))/
    (matrix(rep(colSums(S.ij*matrix(rep(mu.j.prop,3))),each=nrow(S.ij)),ncol=3))
  
  
  
  #Gamma log like
  negLL=min(-sum(tab*(log(mu.j*S.ij))-(mu.j*S.ij),na.rm=T),1e100)
  
  return(list(negLL=negLL,d=d,observed=tab,predicted=tab.pred))
}

  #2.1

  #2.2 SSF

SSF.sp=unique(F2_Sampling$Species)

theta.list=vector('list',length(SSF.sp))
names(theta.list)=SSF.sp
SSF.fit=Combined.sel=theta.list

#initial parameter values
theta.list$`Gummy shark`=c(Theta1=log(80),Theta2=log(29000))
for(s in 1:length(theta.list)) theta.list[[s]]=jitter(theta.list$`Gummy shark`,factor=.1)

#theta=c(Theta1=log(80),Theta2=log(29000))



# fit model
for(s in 1:length(SSF.sp))
{
  theta=theta.list[[s]]
  
  #. objfun to minimize
  fn_ob=function(theta)Selectivty.Kirkwood.Walker(d=F2_Sampling%>%filter(Species==SSF.sp[s]),
                                                  size.int=100,
                                                  theta)$negLL
  
  #. fit model
  SSF.fit[[s]]=nlminb(theta.list[[s]], fn_ob, gradient = NULL)
}



# Calculate confidence intervals thru bootstrapping 
n.boot=1:1000
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
system.time({
  SSF.fit.CI=foreach(s=1:length(SSF.sp),.packages=c('tidyverse','doParallel','Biobase')) %dopar%
    {
      boot=foreach(n=n.boot,.packages=c('doParallel','splitstackshape','tidyverse')) %dopar%
      {
          theta=theta.list[[s]]
          
          #bootstrapped sample
          d.samp=F2_Sampling%>%filter(Species==SSF.sp[s]) 
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
names(SSF.fit.CI)=SSF.sp
stopCluster(cl)




# REPORT  -------------------------------------------------------------------

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
  
  return(S.ij)
}

#2. SSF
for(s in 1:length(SSF.sp)) Combined.sel[[s]]=fn.plt.Sel(Dat=F2_Sampling%>%filter(Species==SSF.sp[s]),theta=SSF.fit[[s]]$par)


#Combined selectivity

#2. SSF
for(s in 1:length(SSF.sp))
{
  Sum.sel=rowSums(Combined.sel[[s]])
  Combined.sel[[s]]$combined=Sum.sel/max(Sum.sel)
  plot(as.numeric(rownames(Combined.sel[[s]])),Combined.sel[[s]]$combined,type='l',lwd=5,col='orange')
  for(l in 1:(ncol(Combined.sel[[s]])-1)) lines(as.numeric(rownames(Combined.sel[[s]])),Combined.sel[[s]][,l])
}
  


#Plot predicted numbers

#2. SSF
for(s in 1:length(SSF.sp))
{
  dummy=Selectivty.Kirkwood.Walker(d=F2_Sampling%>%filter(Species==SSF.sp[s]),
                                   size.int=100,
                                   SSF.fit[[s]]$par)
  smart.par(n.plots=ncol(dummy$observed),MAR=c(.1,.1,.1,.1),OMA=c(2.5,2.5,1.5,.1),MGP=c(1,.5,0))
  for(m in 1:ncol(dummy$observed))
  {
    plot(dummy$observed[,m],dummy$predicted[,m],pch=19)
    lines(dummy$observed[,m],dummy$observed[,m],col=2)
  }
    
}
