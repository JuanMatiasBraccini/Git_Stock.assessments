#                 SCRIPT FOR RECONSTRUCTING TOTAL DISCARDS IN TDGDLF    #

#notes:
# observer data must have at least Min.obs.per.block & Min.shots.per.block
# Ratio estimator method used with all years combined.  
# Model-based by species not used due to small sample sizes for most species
# Uncertainty determined using bootstrapping (as bootstrap resamples the shots, 
#       the procedure is self-weighting, i.e more weight is given to more 
#       commonly sampled blocks)
# Rarely discarded species were grouped as 'other' and then multiplied by proportion
# Longlines: Very few blocks (4 observed blocks only) so few observations and too 
#           much extrapolated!! Hence, grouped commercial LL and GN and use observed GN


#To do:
#  Add post capture mortality to total discard calculation
#  Run Sensitivity tests on key assumptions
# Update Discarding csv file with Parks Australia findings
# Tabulate effort coverage from the observer data
# combine GN and LL total catch by year block
# Show porportion of effort observed by year/gear. 


rm(list=ls(all=TRUE))

options(stringsAsFactors = FALSE)
library(rlang)
#library(dplyr)
#library(tidyr)
library(tidyverse)
library(doParallel)
library(zoo)
library(abind)

User="Matias"
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R")




# 1. Data ---------------------------------------------------------

#Observers data
#keep only elasmobranchs from observed commercial boats
Res.ves=c("HAM","HOU","NAT","FLIN","RV BREAKSEA","RV Gannet","RV GANNET","RV SNIPE 2")
Dat_obs=DATA %>% filter(Mid.Lat<=(-26) & !BOAT%in%Res.ves & Taxa=='Elasmobranch' & !COMMON_NAME=='WHALE')  %>% 
                 select(c(SHEET_NO,date,Month,year,BOAT,zone,BLOCK,SOAK.TIME,MESH_SIZE,
                        MESH_DROP,NET_LENGTH,Mid.Lat,Mid.Long,Method,SPECIES,
                        COMMON_NAME,SCIENTIFIC_NAME,TL,FL,Number)) 
Dat_obs.LL=Dat_obs %>% filter(Method=="LL")
Dat_obs=Dat_obs %>% filter(Method=="GN")


#Catch and effort data   
Dat_total=read.csv('C:\\Matias\\Analyses\\Data_outs\\Data.monthly.csv',stringsAsFactors = F)
Dat_total=Dat_total %>% filter(LAT<=(-26) & Estuary=="NO") %>%
                        mutate(YEAR=YEAR.c,
                               BLOCK=BLOCKX,
                               Catch=LIVEWT.c,
                               BLOCK=ifelse(BLOCK>=96000,paste(floor(abs(LAT)),(floor(LONG)-100)*10,sep=""),BLOCK),
                               BLOCK=substring(BLOCK,1,4)) %>%
                        select(-c(ZnID,MonthlyID,ZoneID,YEAR.c,blockxFC,CONDITN,
                                  Factor,RSCommonName,RSSpeciesId,LANDWT,LIVEWT,
                                  FisheryZone,FisheryCode,Landing.Port,BDAYS,licence,
                                  Factor.c,LIVEWT.orgnl,LIVEWT.c,VesselID,BlockAveID,AnnualVesselAveID,
                                  BlockID,Reporter,Sch.or.DogS,LIVEWT.reap,Tot.shk.livewt,
                                  Shark.other.livewt,Reporter.old,Spec.old,Sname.old,NETLEN.c)) 
Dat_total.LL=Dat_total %>% filter(METHOD=="LL")
Dat_total=Dat_total %>% filter(METHOD=="GN")

#Discarded/retained
Comm.disc.sp=read.csv("C:/Matias/Analyses/Ecosystem indices and multivariate/Shark-bycatch/SPECIES+PCS+FATE.csv",stringsAsFactors = F)

#Length weight relationships
setwd('C:/Matias/Analyses/Reconstruction_total_bycatch_TDGDLF')

Len.wei=read.csv("length.weights.csv",stringsAsFactors = F)



# 2. Parameter ---------------------------------------------------------

Min.obs.per.block=10  #minimum number of observations per block for use in ratio estimator
Min.shots.per.block=5  #minimum number of shots per block for use in ratio estimator

do.explrtn="NO"
STRTA.obs.disc=c("BLOCK")       #aggregating strata for observer data
STRTA.reported.ret=c("BLOCK")   #aggregating strata for total reported retained catch

n.boot=1e3

Group_rare_criteria=0.02    #criteria for grouping rare species (proportion of catch)    

Commercial.sp=subset(Comm.disc.sp,FATE=="C" & NATURE=="S")$SPECIES



#Use only observed GN. Combine the total catch of GN and LL 
DATA_obs=list(GN=Dat_obs,LL=Dat_obs.LL)
DATA_total=list(GN=Dat_total,LL=Dat_total.LL)  



# 3. Procedure ---------------------------------------------------------

#define discarded/retained sp
for(i in 1:length(DATA_obs)) DATA_obs[[i]]$Discarded=with(DATA_obs[[i]],
                                          ifelse(SPECIES%in%Commercial.sp,"Retained","Discarded"))

#Show overal observed discarded and retained species by year   #Replace by this infographic: https://twitter.com/ISSF/status/1147943595942526978?s=03
fun.horiz.bar=function(d)
{
  d$dummy=with(d,ifelse(Discarded=="Retained","Retained",SPECIES))
  TAB=table(d$dummy)
  TAB=TAB/sum(TAB)
  TAB=sort(TAB)
  a=barplot(TAB,horiz=T,xlim=c(0,1),las=1)
  box()
}
axlbl=c("Gillnet","Longline")
jpeg("Results/Figure1_observed.proportions.jpg",width=2400,height=2400,units="px",res=300)
par(mfcol=c(1,2),mar=c(1,2.5,1,1),oma=c(1,2.75,1,1),mgp=c(1.5,.4,0))
for(i in 1:length(DATA_obs))
{
  fun.horiz.bar(d=DATA_obs[[i]])
  mtext(axlbl[i],3)
}
mtext("Proportion",1,outer=T)
dev.off()


#fill in missing length info 
#note: first use TL, then sample from species distribution, finally overall mean
UniK.sp=table(DATA_obs$GN$SPECIES)
UniK.sp=UniK.sp[UniK.sp>3]
UniK.sp=names(UniK.sp)
pdf("Results/Preliminary/lengths.pdf")
for(u in 1:length(UniK.sp))
{
  a=subset(DATA_obs$GN,SPECIES==UniK.sp[u] & !is.na(FL))
  if(nrow(a)>2)hist(a$FL,main=UniK.sp[u],col=2)
}
dev.off()

  #1. proportion of TL
for(i in 1:length(DATA_obs)) DATA_obs[[i]]$FL=with(DATA_obs[[i]],ifelse(is.na(FL),TL*.875,FL))

  #2. samples from species distribution
fn.samp.dist=function(d)
{
  d=subset(d,FL>20)
  his=hist(d$FL,breaks=seq(floor(min(d$FL,na.rm=T)),ceiling(max(d$FL,na.rm=T)),by=1),plot=F)
  return(his)
}
for(i in 1:length(DATA_obs))
{
  all.sp=unique(DATA_obs[[i]]$SPECIES)
  for(s in 1:length(all.sp)) 
  {
    d=subset(DATA_obs[[i]],SPECIES==all.sp[s])
    if(sum(is.na(d$FL))>0 & sum(!is.na(d$FL))>0)
    {
      if(length(d$FL[!is.na(d$FL)])>2)Prob=fn.samp.dist(d)
      if(length(Prob$density)>2)
      {
        id=which (is.na(DATA_obs$FL) & DATA_obs$SPECIES==all.sp[s])
        DATA_obs[[i]]$FL[id]=sample(Prob$breaks[-1],
                                length(DATA_obs[[i]]$FL[id]),
                                replace=T,prob=Prob$density/sum(Prob$density))
      }
    }
  }
}

  #3. overall mean
for(i in 1:length(DATA_obs)) DATA_obs[[i]]$FL[is.na(DATA_obs[[i]]$FL)]=mean(DATA_obs[[i]]$FL,na.rm=T)

pdf("Results/Preliminary/lengths_ammended.pdf")
for(u in 1:length(UniK.sp))
{
  a=subset(DATA_obs$GN,SPECIES==UniK.sp[u] & !is.na(FL))
  if(nrow(a)>2)hist(a$FL,main=UniK.sp[u],col=2)
}
dev.off()

#Convert number to weight               
#note: calculate discard ratio using weight as total catch is in weight
Len.wei=Len.wei%>%select(c(SPECIES,a,b))
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]=merge(DATA_obs[[i]],Len.wei,by="SPECIES",all.x=T)
  DATA_obs[[i]]$Catch=with(DATA_obs[[i]],b*(FL)^a)
}


#group rare species
  #tabulate rare species
for(i in 1:length(DATA_obs))
{
  Tab.sp=prop.table(with(subset(DATA_obs[[i]],Discarded=="Discarded"),table(SPECIES)))
  Rare.sp=names(Tab.sp[Tab.sp<Group_rare_criteria])
  DATA_obs[[i]]$SPECIES.ori=DATA_obs[[i]]$SPECIES
  DATA_obs[[i]]$SPECIES=with(DATA_obs[[i]],ifelse(SPECIES%in%Rare.sp,"Grouped",SPECIES))
}


  #Discarded_sp used to calculate discarded species proportions
for(i in 1:length(DATA_obs))
{
  DATA_obs[[i]]$Discarded_sp=with(DATA_obs[[i]],ifelse(Discarded=="Retained","Retained",SPECIES))
}


#3.1 Data exploration
if(do.explrtn=="YES")
{
  A=with(Dat_obs[!duplicated(Dat_obs$SHEET_NO),],table(BLOCK))
  Exp=5*(A/max(A))
  Exp=ifelse(Exp<0.25,0.25,Exp)
  Exp=data.frame(LAT=-(as.numeric(substr(names(Exp),1,2))+.5),
                 Long=100+as.numeric(substr(names(Exp),3,4))+.5,cex=Exp)
  pdf("Results/Preliminary/shots_block_map.pdf")
  with(Exp,plot(Long,LAT,cex=cex,pch=19,col='steelblue'))
  legend('topright',paste(c(10,100,250)),pch=19,pt.cex=5*(c(10,100,250)/max(A)),
         bty='n',col='steelblue',title="# shots")
  dev.off()
  
  YRs=sort(unique(Dat_obs$year))
  pdf("Results/Preliminary/shots_block_year_map.pdf")
  par(mfcol=c(5,4),mar=c(1,1.5,1,.5),oma=c(3,1.25,.5,.1),las=1,mgp=c(1.5,.5,0),cex.axis=.85)
  for(y in 1:length(YRs))
  {
    x=subset(Dat_obs,year==YRs[y])
    A=with(x[!duplicated(x$SHEET_NO),],table(BLOCK))
    Exp=3.5*(A/max(A))
    Exp=ifelse(Exp<0.25,0.25,Exp)
    Exp=data.frame(LAT=-(as.numeric(substr(names(Exp),1,2))+.5),
                   Long=100+as.numeric(substr(names(Exp),3,4))+.5,cex=Exp)
    with(Exp,plot(Long,LAT,cex=cex,pch=19,col='steelblue',ylim=c(-36,-26),
                  xlim=c(113,129),ann=F))
    Quant=round(quantile(A,prob=c(.25,.5,.75)))
    legend('topright',paste(Quant),pch=19,pt.cex=3.5*(Quant/max(A)),
           bty='n',col='steelblue',title="# shots")
    legend('topleft',paste(YRs[y]),bty='n')
  }
  dev.off()
  
  #blocks by year
  pdf("Results/Preliminary/shots_block_year_plot.pdf")
  A=with(Dat_obs[!duplicated(Dat_obs$SHEET_NO),],table(year,BLOCK))
  plot(1,1,col="transparent",ylim=c(1,ncol(A)),xlim=c(1,nrow(A)),ylab="block",xlab="year")
  for(i in 1:nrow(A))
  {
    Exp=5*(A[i,]/max(A))
    points(rep(i,ncol(A)),1:ncol(A),cex=Exp,pch=19,col='steelblue')
  }
  legend('bottomright',paste(c(10,30,80)),pch=19,pt.cex=5*(c(10,30,80)/max(A)),
         bty='n',col='steelblue',title="# shots")
  dev.off()
  
  #Shark and ray species
  plot(sort(table(Dat_obs$SPECIES)),type='h',ylab="# individuals")
}

#remove stuff
rm(DATA,Dat_obs,Dat_obs.LL,Dat_total,Dat_total.LL,DATA.bio,DATA.ecosystems)


#3.2 Get stratified observed discarded and retained catch (R_h)
STRATA_obs=function(d,Strata)
{
  #select minimum number of observations for use in analysis
  tt <- table(d$BLOCK)
  d <- d[d$BLOCK %in% names(tt[tt >= Min.obs.per.block]), ]
  
  tt <- with(d[!duplicated(d$SHEET_NO),],table(BLOCK))
  d <- d[d$BLOCK %in% names(tt[tt >= Min.shots.per.block]), ]
  
  #get discard : retain ratios
  a=d %>%
    group_by(.dots=Strata) %>%
    summarise(total = sum(Catch))%>%
    spread(key = !! parse_expr(Strata[1]), value = total)%>%
    as.data.frame
  a[is.na(a)]=0
  Vars=names(a)[-match(c(Strata[-1],"Retained"),names(a))]
  a[,Vars]=a[,Vars]/a$Retained
  names(a)[match(Vars,names(a))]=paste(names(a[,Vars]),".ratio",sep="")

  #get bycatch proportional species composition 
  prop=subset(d,Discarded=="Discarded") %>%
    group_by(.dots=c("SPECIES.ori",Strata[-1])) %>%
    summarise(total = sum(Catch))%>%
    spread(key = SPECIES.ori, value = total)%>%
    as.data.frame
  prop[is.na(prop)]=0
  Vars.prop=names(prop)[-match(c(Strata[-1]),names(prop))]
  prop[,Vars.prop]=prop[,Vars.prop]/rowSums(prop[,Vars.prop])
  
  #keep only proportion for grouped species
  id.prop=names(prop)[which(!names(prop)%in%Vars)]
  if(length(id.prop)>1)
  {
    prop=prop[,id.prop]
    
    #rescale to sum 'Grouped' to 1
    prop[,-1]=prop[,-1]/rowSums(prop[,-1])
    prop[is.na(rowSums(prop[,-1])),-1]=0
  }

  
  return(list(dat=a,prop=prop,vars=paste(Vars,".ratio",sep="")))
}
Obs_ratio.strata=STRATA_obs(d=DATA_obs$GN,Strata=c("Discarded_sp",STRTA.obs.disc))
write.csv(Obs_ratio.strata$dat,"Results/Table1_observed.ratios.csv",row.names = F)



#3.3 Get stratified total reported retained catch (Y_h)
STRATA_total=function(d,Strata)
{
  if(!"YEAR"%in%Strata) Strata=c("YEAR",Strata)
  a=d %>%
    group_by(.dots=Strata) %>%
    summarise(total.retained = sum(Catch))%>%
    as.data.frame
  return(a)
}
Total_strata=STRATA_total(d=DATA_total$GN,Strata=STRTA.reported.ret)


#3.4 Annual discards by species (sum(R_h x Y_h))       
fn.total=function(disc.dat,tot.dat)
{
  #combine observed ratios and total reported catch
  Total.discard=merge(disc.dat$dat,tot.dat,by=STRTA.reported.ret,all.y=T)
  
  #linear interpolate missing block discard ratio
  Vars=disc.dat$vars
  for(v in 1:length(Vars)) 
  {
    id=match(Vars[v],names(Total.discard))
    Total.discard[,id]=ifelse(is.na(Total.discard[,id]),
              na.approx(zoo(Total.discard[,id])),Total.discard[,id]) 
  }
  
  #multiply ratio by total
  Total.discard[,Vars]=Total.discard[,Vars]*Total.discard$total.retained
  
  #rename
  i.dvar=match(Vars,names(Total.discard))
  names(Total.discard)[i.dvar]= paste("total",
      sapply(strsplit(names(Total.discard)[i.dvar], ".", fixed = TRUE), "[", 1),sep=".")
  
  #Split up 'Other'
  if(!is.na(match("total.Grouped",names(Total.discard))))
  {
    dummy=Total.discard[,c("BLOCK","YEAR","total.Grouped")]
    dummy=merge(disc.dat$prop,dummy,by="BLOCK",all.y=T)
    
    #linear interpolate missing block discard ratio
    Vars=subset(names(disc.dat$prop),!names(disc.dat$prop)=="BLOCK")
    for(v in 1:length(Vars)) 
    {
      id=match(Vars[v],names(dummy))
      dummy[,id]=ifelse(is.na(dummy[,id]),na.approx(zoo(dummy[,id])),dummy[,id]) 
    }
    
    #multiply ratio by total
    dummy[,Vars]=dummy[,Vars]*dummy$total.Grouped
    
    #rename
    i.dvar=match(Vars,names(dummy))
    names(dummy)[i.dvar]= paste("total",
        sapply(strsplit(names(dummy)[i.dvar], ".", fixed = TRUE), "[", 1),sep=".")
    
    #merge to total.discard
    Total.discard=merge(Total.discard,dummy,by=c("BLOCK","YEAR","total.Grouped"))
    Total.discard=Total.discard[,-match("total.Grouped",names(Total.discard))]
  }
  Total.discard=Total.discard[,-match("Retained",names(Total.discard))]
  return(Total.discard)
}
Tot.discard.result=fn.total(disc.dat=Obs_ratio.strata,tot.dat=Total_strata)

#show blocks interpolated by species
fn.number.interpolated=function(disc.dat,tot.dat)
{
  tot.dat=tot.dat[!duplicated(tot.dat$BLOCK),]
  
  #combine observed ratios and total reported catch
  Total.discard=merge(disc.dat$dat,tot.dat,by=STRTA.reported.ret,all.y=T)
  
  #count NAs by species
  return(Total.discard %>% select (-c(BLOCK, total.retained,YEAR)) %>%
           summarise_all(list(~sum(is.na(.))))%>%
           as.data.frame)
}
interpolated.blocks=fn.number.interpolated(disc.dat=Obs_ratio.strata,tot.dat=Total_strata)
write.csv(interpolated.blocks,"Results/Table2_interpolated.blocks.csv",row.names=F)


#3.5 Uncertainty thru non-parametric bootstrap    Takes 0.002 secs per iteration
cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)
clusterCall(cl, function() {
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(zoo)
})
Results=Tot.discard.result
set.seed(666)
system.time({
  #minimum data for analysis
  tt <- table(DATA_obs$GN$BLOCK)
  d <- DATA_obs$GN[DATA_obs$GN$BLOCK %in% names(tt[tt >= Min.obs.per.block]), ]
  
  tt <- with(d[!duplicated(d$SHEET_NO),],table(BLOCK))
  d <- d[d$BLOCK %in% names(tt[tt >= Min.shots.per.block]), ]
  
  Shots=unique(d$SHEET_NO)
  n.shots=length(Shots)
  
  store=foreach(n=1:n.boot) %dopar%
  {
    #resample observer data
    id=sample(Shots, n.shots, replace=TRUE)
    Tab=table(id)
    Tab=data.frame(SHEET_NO=names(Tab),Rep=as.numeric(Tab))
    Dat_obs.boot=d[which(d$SHEET_NO%in%id),]
    Dat_obs.boot=merge(Dat_obs.boot,Tab,by="SHEET_NO")
    Dat_obs.boot <- Dat_obs.boot[rep(row.names(Dat_obs.boot), Dat_obs.boot$Rep), ]
    
    #bootstrapped stratified observed R_h
    Obs_ratio.strata.boot=STRATA_obs(d=Dat_obs.boot,Strata=c("Discarded_sp",STRTA.obs.disc))
    
    #total discards (sum(R_h x Y_h))   
    Tot.discard.boot=fn.total(disc.dat=Obs_ratio.strata.boot,tot.dat=Total_strata)
    
    rm(Dat_obs.boot)
    return(Tot.discard.boot)
  }
  Results=store
  rm(d)

})
rm(store)
stopCluster(cl) 


# 4. Report ---------------------------------------------------------

#Put results in right format, aggregate blocks and extract median, ci    
fn.agg.block=function(d)
{
  dummy=vector('list',n.boot)
  for(n in 1:n.boot)
  {
    Var=subset(names(d[[n]]),!names(d[[n]])%in%c("BLOCK","YEAR"))
    Var.ag=paste('cbind(',paste(Var,collapse=','),")",sep='')
    Formula=as.formula(paste(Var.ag,"YEAR",sep="~"))
    dummy[[n]]=aggregate(Formula,d[[n]],sum)
  }
  return(dummy)
}

fn.add.missing=function(d,VAR)
{
  VAR=subset(VAR,!VAR%in%c("BLOCK"))
  VAR1=subset(VAR,!VAR%in%c("YEAR","total.retained"))
  for(n in 1:n.boot)
  {
    id=which(!VAR%in%names(d[[n]]))
    if(length(id>0))
    {
      iid=which(!VAR1%in%names(d[[n]]))
      dummy=as.data.frame(d[[n]][,rep(1,length(iid))])
      names(dummy)=VAR1[iid]
      dummy[,]=0
      d[[n]]=cbind(d[[n]],dummy)
    }
    d[[n]]=d[[n]][,match(VAR,names(d[[n]]))]
  }
  return(d)
}

#aggregate blocks and add missing species     #0.018 sec per iteration
system.time({ 
  Results.show=fn.add.missing(d=fn.agg.block(d=Results),VAR=names(Tot.discard.result))
})

# Stats
PRBS=c(.025,.5,.975)
Stat.list=list(Species="Species",Total="Total")
Stats=Stat.list
fn.stats=function(d,how)
{
  if(how=="Species")
  {
    all.matrix <- abind(d, along=3)
    Store=vector('list',ncol(all.matrix))
    names(Store)=colnames(all.matrix)
    for(nm in 2:ncol(all.matrix)) Store[[nm]]=apply(all.matrix[,nm,],1,function(x) quantile(x/1000,probs=PRBS)) #in tonnes
  }
  if(how=="Total")
  {
    id.var=which(!colnames(d[[1]])%in%c("total.retained","total.retained","YEAR"))
    for(n in 1:n.boot) d[[n]]$total.discarded=rowSums(d[[n]][,id.var])
    all.matrix <- abind(d, along=3)
    These=c("YEAR","total.discarded","total.retained")
    Store=vector('list',length(These))
    names(Store)=These
    id.these=match(These,colnames(d[[1]]))
    for(nm in 2:length(id.these)) Store[[nm]]=apply(all.matrix[,id.these[nm],],1,function(x) quantile(x/1000,probs=PRBS))
  }
  Store$YEAR=all.matrix[,,1][,1]

   return(Store)
}
for(l in 1:length(Stats)) Stats[[l]]=fn.stats(d=Results.show,how=Stat.list[[l]])


#Plot statistics
yrs=sort(unique(DATA_total$GN$YEAR))
Col.totl="black"
Col.totl.retained="grey60"

#plotting function
id.med=match("50%",rownames(Stats$Species[[2]]))
id.low=match("2.5%",rownames(Stats$Species[[2]]))
id.up=match("97.5%",rownames(Stats$Species[[2]]))


#Total
fn.plt.prop=function(d,d.ret,CL,LWd)
{
  plot(yrs,yrs,col='transparent',ylim=c(0,YMAX),ann=F)
  polygon(x=c(yrs,rev(yrs)),y=c(d[id.low,]/d.ret,rev(d[id.up,]/d.ret)),
          col=adjustcolor(CL, alpha.f = 0.30), border = adjustcolor(CL, alpha.f = 0.60))
  lines(yrs,d[id.med,]/d.ret,lwd=LWd,col=CL)
}
YMAX=1
#YMAX=max(c(unlist(Stats$Total$total.retained),unlist(Stats$Total$total.discarded)))
jpeg("Results/Figure2_total.jpg",width=2400,height=2400,units="px",res=300)
par(mfcol=c(1,1),mar=c(1.5,1,1,1),oma=c(1,3,1,1),mgp=c(1,.6,0),las=1)
fn.plt.prop(d=Stats$Total$total.discarded,d.ret=Stats$Total$total.retained[id.med,],CL=Col.totl,LWd=3)
mtext("Proportion",2,1.5,las=3,outer=T,cex=1.5)
#lines(yrs,Stats$Total$total.retained[id.med,],lwd=3,col=Col.totl.retained)
#legend("topright",c("retained","discarded"),lty=1,lwd=3,
#       bty='n',col=c(Col.totl.retained,Col.totl),cex=1.5)
#mtext("Total (tonnes)",2,1.5,las=3,outer=T,cex=1.5)
mtext("Year",1,0,outer=T,cex=1.5)
dev.off()

#By species
fn.plt=function(d,CL,LWd)
{
  plot(yrs,yrs,col='transparent',ylim=c(0,YMAX),ann=F)
  polygon(x=c(yrs,rev(yrs)),y=c(d[id.low,],rev(d[id.up,])),
          col=adjustcolor(CL, alpha.f = 0.30), border = adjustcolor(CL, alpha.f = 0.60))
  lines(yrs,d[id.med,],lwd=LWd,col=CL)
}
jpeg("Results/Figure3_by.species.jpg",width=2400,height=2400,units="px",res=300)
Plt.this=Stats$Species[-match(c("YEAR","total.retained"),names(Stats$Species))]
Plt.this=Plt.this[sort(names(Plt.this))]
smart.par(n.plots=length(Plt.this),MAR=c(1,1.5,1,1),OMA=c(2,2,.1,.1),MGP=c(1,.6,0))
for(i in 1:length(Plt.this))
{
  YMAX=max(unlist(Plt.this[[i]]))
  fn.plt(d=Plt.this[[i]],CL=Col.totl,LWd=1.5)
  mtext(unlist(strsplit(names(Plt.this)[i], ".", fixed = TRUE))[2],3,-1.51,cex=1) 
}
mtext("Year",1,0.7,outer=T,cex=1.25)
mtext("Total discard (tonnes)",2,0.5,outer=T,las=3,cex=1.25)
dev.off()

#Output total discard estimates
  #total
colnames(Stats$Total$total.discarded)=yrs
write.csv(Stats$Total$total.discarded,
      paste("Results/Total.discard.estimates/Total.csv",sep=""))

  #by species
for(i in 1:length(Plt.this)) 
{
  nm=unlist(strsplit(names(Plt.this)[i], ".", fixed = TRUE))[2]
  colnames(Plt.this[[i]])=yrs
  write.csv(Plt.this[[i]],paste("Results/Total.discard.estimates/",nm,".csv",sep=""))
}
  
