#------------- Script for bringing in all shark data  ----------- 

library(stringr)
library(tidyr)
library(Hmisc)
library(dplyr)
options(stringsAsFactors = FALSE) 

#Total catch color
tot.col="tan3"

#Zone catch colors
Zns=c("Joint","North","Closed","West","Closed.metro","Zone1","Zone2")
Zns.leg=c("JANSF","WANCSF","Closed","WCDGDLF","Metro closure","Zone1","Zone2")
COL.prop=c("aquamarine2","lightseagreen",'forestgreen',"lightgreen","olivedrab4","olivedrab3","mediumseagreen")
names(COL.prop)=names(Zns.leg)=Zns

#Biological catch regions colors
Bio.col=c("dodgerblue","darkorchid4","cyan4","lightpink3")
names(Bio.col)=c("North Coast","Gascoyne","West Coast","South Coast")


if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

if(!exists('fn.word.table')) source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/MS.Office.outputs.R"))

#------------- Functions  ----------- 
if(!exists('fn.fig')) fn.fig=function(NAME,Width,Height,Do.tiff="YES",Do.jpeg="NO")
{
  if(Do.tiff=="YES") tiff(file=paste(NAME,".tiff",sep=""),width=Width,height=Height,units="px",res=300,compression="lzw")
  if(Do.jpeg=="YES") jpeg(file=paste(NAME,".jpeg",sep=""),width=Width,height=Height,units="px",res=300)
}
if(!exists('smart.par')) smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))

fn.bring.all.data=function(path)
{
  if(path%in%list.files(Dat.repository))
  {
    temp.wd=paste(Dat.repository,path,sep='')
    Files=list.files(temp.wd)
    setwd(temp.wd)
    Files1=lapply(Files,read.csv)
    names(Files1)= str_remove(Files, path)
    return(Files1)
  }
}

fn.subs=function(YEAR) substr(YEAR,start=3,stop=4)

fn.finyr=function(Finyr)paste(Finyr-1,"-",fn.subs(Finyr),sep="")

fn.month.to.fin.mn=function(dat)ifelse(dat>6,dat-6,dat+6)

fn.exp.his=function(dat,dat1)
{
  dat2=expand.grid(FINYEAR=dat$FINYEAR,MONTH=1:12)
  dat2=merge(dat2,dat1[,match(c("MONTH","prop"),names(dat1))],by="MONTH",all.x=T)
  dat2=merge(dat2,dat,by="FINYEAR",all.x=T)
  dat2=dat2[order(dat2$FINYEAR,dat2$MONTH),]
  dat2$LIVEWT.c=dat2$LIVEWT.c*dat2$prop
  return(dat2[,-match("prop",names(dat2))])
}

fn.rename.dat=function(x,y) str_remove_all(x, paste(y, collapse = "|"))

#Functions for adding missing size
fn.add.missing.size=function(dat,Min,Max,interval)
{
  Rango=c(Min,Max)
  SEQ=seq(Rango[1],Rango[2],interval)
  dat=dat[order(dat$year),]
  dat$FL.bin=floor(dat$FL/interval)*interval
  dat$FL.bin=factor(dat$FL.bin,levels=SEQ)
  
  dat$FINYEAR=as.character(with(dat,
                                ifelse(Month%in%1:6,paste(year-1,"-",fn.subs(year),sep=""),
                                       ifelse(Month%in%7:12,paste(year,"-",fn.subs(year+1),sep=""),NA))))
  
  if(Conv.cal.mn.to.fin.mn=="YES") dat$Month=fn.month.to.fin.mn(dat$Month)
  
  #   dat=aggregate(Number~FL.bin+FINYEAR+Type,dat,sum)
  #   if(!is.null(dat$Sex))dat=aggregate(Number~FL.bin+FINYEAR+Type+Sex,dat,sum)
  #   
  #   dat2=subset(dat,Number>0,select=c(FL.bin,FINYEAR,Number))
  #   Table=reshape(dat2,v.names = "Number", idvar = "FINYEAR",timevar = "FL.bin", direction = "wide")
  #   X=colnames(Table)[2:ncol(Table)]
  #   colnames(Table)[2:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
  #   ID=which(!levels(dat$FL.bin)%in%colnames(Table)[2:ncol(Table)])
  #   ADD=levels(dat$FL.bin)[ID]
  #   ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
  #   names(ADD1)=ADD
  #   Table=cbind(Table,ADD1)
  #   Table=Table[,c(1,match(levels(dat$FL.bin),names(Table)))]
  #   Table[is.na(Table)]=0 
  #   
  dat=aggregate(Number~FL.bin+Month+FINYEAR+Type+Sex,dat,sum)
  
  dat2=subset(dat,Number>0,select=c(FL.bin,FINYEAR,Number,Month,Sex))
  Table=reshape(dat2,v.names = "Number", idvar = c("FINYEAR","Month","Sex"),timevar = "FL.bin", direction = "wide")
  X=colnames(Table)[4:ncol(Table)]
  colnames(Table)[4:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
  ID=which(!levels(dat$FL.bin)%in%colnames(Table)[4:ncol(Table)])
  ADD=levels(dat$FL.bin)[ID]
  ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
  names(ADD1)=ADD
  Table=cbind(Table,ADD1)
  Table=Table[,c(1:3,match(levels(dat$FL.bin),names(Table)))]
  Table[is.na(Table)]=0 
  Table=Table[order(Table$Sex,Table$FINYEAR,Table$Month),]
  
  return(Table)
  
}
fn.add.missing.size2=function(dat,Min,Max,interval)
{
  Rango=c(Min,Max)
  SEQ=seq(Rango[1],Rango[2],interval)
  dat=dat[order(dat$year),]
  dat$FL.bin=floor(dat$FL/interval)*interval
  dat$FL.bin=factor(dat$FL.bin,levels=SEQ)
  dat$Sex=dat$SEX
  dat$Number=1
  dat=subset(dat,!Sex=="U")
  
  dat=aggregate(Number~FL.bin+Month+FINYEAR+Sex,dat,sum)
  
  dat2=subset(dat,Number>0,select=c(FL.bin,FINYEAR,Number,Month,Sex))
  Table=reshape(dat2,v.names = "Number", idvar = c("FINYEAR","Month","Sex"),timevar = "FL.bin", direction = "wide")
  X=colnames(Table)[4:ncol(Table)]
  colnames(Table)[4:ncol(Table)]=as.character(sapply(strsplit(X,"Number."), "[", 2))
  ID=which(!levels(dat$FL.bin)%in%colnames(Table)[4:ncol(Table)])
  ADD=levels(dat$FL.bin)[ID]
  ADD1=as.data.frame(matrix(nrow=nrow(Table),ncol=length(ADD)))
  names(ADD1)=ADD
  Table=cbind(Table,ADD1)
  Table=Table[,c(1:3,match(levels(dat$FL.bin),names(Table)))]
  Table[is.na(Table)]=0 
  Table=Table[order(Table$Sex,Table$FINYEAR,Table$Month),]
  return(Table)
  
}
fn.add.missing.size3=function(dat,Min,Max,interval)
{
  Rango=c(Min,Max)
  SEQ=seq(Rango[1],Rango[2],interval)
  Table=dat%>%
    rename(Sex=SEX)%>%
    filter(!is.na(Sex) | !Sex=="U")%>%
    mutate(TL.bin=floor(TL/interval)*interval)%>%
    group_by(TL.bin,FINYEAR,Month,Sex)%>%
    summarise(N=n())
  add.msn=SEQ[which(!SEQ%in%Table$TL.bin)]  
  if(length(add.msn)>0)
  {
    addd=Table[1:length(add.msn),]
    addd$N=0
    addd$TL.bin=add.msn
    Table=rbind(Table,addd)
  }
  Table=Table%>%
    spread(TL.bin,N,fill=0)%>%
    data.frame%>%
    arrange(Sex,FINYEAR,Month)%>%
    filter(!is.na(FINYEAR))
  
  colnames(Table)[-match(c("FINYEAR","Month","Sex"),colnames(Table))]=
    str_remove_all(colnames(Table)[-match(c("FINYEAR","Month","Sex"),colnames(Table))], "X")
  return(Table)
  
}

#catch visualization functions
fn.add=function(a,YRS)
{
  id=which(!YRS%in%a$FINYEAR)
  if(length(id)>0)
  {
    ADD=a
    ADD=ADD[1:length(id),]
    ADD$FINYEAR=YRS[id]
    xx=match("FINYEAR",colnames(ADD))
    ADD[,which(!1:ncol(ADD)%in%xx)]=NA
    a=rbind(a,ADD)
    a=a[order(a$FINYEAR),]
  }
  return(a)
}
fn.see.Ktch=function(DAT,DAT.zone,YR.span,LWD,TCK)
{
  DAT$Total=do.call(rbind,DAT)%>%
    mutate(FishCubeCode='Total',
           Data.set='Total',
           Region='Total',
           Fishery='Total')%>%
    group_by(SPECIES,Name,FishCubeCode,Data.set,finyear,FINYEAR,Region,Fishery)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
  DAT.zone$Total=do.call(rbind,DAT.zone)%>%
    mutate(FishCubeCode='Total',
           Data.set='Total',
           Region='Total',
           Fishery='Total',
           zone=NA)%>%
    group_by(SPECIES,Name,FishCubeCode,Data.set,finyear,FINYEAR,Region,zone,Fishery)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
  for(i in 1:length(DAT))
  {
    if(names(DAT)[i]=='Historic') names(DAT)[i]='Historical'
    a=DAT[[i]]
    a.zn=DAT.zone[[i]]%>%filter(!is.na(zone))
    a=a%>%
      group_by(FINYEAR)%>%
      summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))%>%
      arrange(FINYEAR)
    Ylim=c(0,max(a$LIVEWT.c))
    if(nrow(a)==1)
    {
      plot(1:length(a$FINYEAR),a$LIVEWT.c,xaxt='n',ylab="",xlab="",ylim=Ylim,
           main=names(DAT)[i],col=tot.col)
    }else
    {
      disKol=tot.col
      disLWD=LWD*.75
      if(names(DAT)[i]=="Total")
      {
        disKol='black'
        disLWD=LWD*1.25
      }
        
      a=fn.add(a,YR.span)
      plot(1:length(a$FINYEAR),a$LIVEWT.c,lwd=disLWD,xaxt='n',ylab="",xlab="",
           ylim=Ylim,type='l',main=names(DAT)[i],col=disKol)
      if(names(DAT)[i]=="TDGDLF")
      {
      if(nrow(a.zn)>0)
      {
        a.zn=aggregate(LIVEWT.c~FINYEAR+zone,a.zn,sum)
        dis.z=unique(a.zn$zone)
        for(z in 1:length(dis.z))
        {
          b=subset(a.zn,zone==dis.z[z])
          b=fn.add(b,YR.span)
          lines(1:length(b$FINYEAR),b$LIVEWT.c,lwd=LWD,col=COL.prop[match(dis.z[z],names(COL.prop))])
        }
        lines(1:length(a$FINYEAR),a$LIVEWT.c,lwd=LWD*.75,col=tot.col)
      }
      LGn=c('Total')
      CLS=tot.col
      if(nrow(a.zn)>0)
      {
        Cl.z=COL.prop[match(dis.z,names(COL.prop))]
        LGn=c(LGn,Zns.leg[match(dis.z,names(Zns.leg))])
        CLS=c(CLS,Cl.z)
      }
      LGN.location='topright'
      if(which.max(a$LIVEWT.c)>20)LGN.location='topleft'
      legend(LGN.location,LGn,lty=1,lwd=2,col=CLS,bty='n')
      }
    }
    axis(1,1:length(a$FINYEAR),F,tck=TCK)
    axis(1,seq(1,length(a$FINYEAR),5),a$FINYEAR[seq(1,length(a$FINYEAR),5)],tck=2*TCK)
  }
}

#function for visualizing available data
visualize.dat=function(d.list,YR.span)   
{
  cum.plots=rep(NA,length(d.list))
  for(l in 1:length(d.list)) cum.plots[l]=length(d.list[[l]])
  cum.plots=cumsum(cum.plots)
  n.rows=cum.plots[length(cum.plots)]
  LABS=vector('list',length(d.list))
  plot(as.numeric(substr(YR.span,1,4)),as.numeric(substr(YR.span,1,4)),ylim=c(1,n.rows),yaxt='n',col='transparent',ylab='',xlab='')
  All.cols=matrix(hcl.colors(length(d.list)*2,palette = 'Geyser'),ncol=2,byrow=T)
  aidi=which(names(d.list$Catch)=="Historic")
  if(length(aidi)>0) names(d.list$Catch)[aidi]="Historical"
  for(l in 1:length(d.list))
  {
    dd=d.list[[l]]
    LABS[[l]]=names(d.list)[l]
    dummy=rep(NA,length(dd))
    colfunc <- function(n,col1,col2) colorRampPalette(c(col1, col2))(n)
    COL=colfunc(length(dd),col1=All.cols[l,1],col2=All.cols[l,2])
    for(ll in 1:length(dd))
    {
      if(length(dd)>1)dummy[ll]=paste(LABS[[l]]," (",names(dd)[ll],")",sep='')
      if(is.character(dd[[ll]]))x=dd[[ll]] else
      {
        colnames(dd[[ll]])=tolower(colnames(dd[[ll]]))
        x=dd[[ll]]$finyear
      }
      x=as.numeric(substr(unique(x),1,4))
      if(l==1) y=ll
      if(l>1) y=cum.plots[l-1]+ll
      points(x,rep(y,length(x)),cex=6,pch="-",col=COL[ll])
      rm(x)
    }
    if(length(dd)>1) LABS[[l]]=dummy
  }
  axis(2,1:n.rows,unlist(LABS),las=2,cex.axis=1)
  mtext("Financial year",1,cex=1.25,line=1.75)
}

#function for exporting used data
fn.agg.at.level.and.exprt=function(DAT,Level,VAR,SOURCE)
{
  for(i in 1:length(DAT))
  {
    if(Level=='annual')
    {
      a=aggregate(LIVEWT.c~FINYEAR,DAT[[i]],sum)
      write.csv(a,paste(VAR,".",Level,".",SOURCE[i],".csv",sep=""),row.names=F)
    }
    
    if(Level=='annual.by.zone')
    {
      if(!is.null(DAT[[i]]$zone))
      {
        DAT[[i]]$zone=with(DAT[[i]],ifelse(is.na(zone),SOURCE[i],zone))
        a=aggregate(LIVEWT.c~FINYEAR+zone,DAT[[i]],sum)
        write.csv(a,paste(VAR,".",Level,".",SOURCE[i],".csv",sep=""),row.names=F)
      }
    }
  }
}

#function for exporting size composition
fn.exp.size=function(dat,Level)
{
  n=length(dat)
  for( i in 1:n)
  {
    a=dat[[i]]
    
    if(class(a)=="data.frame") 
    {
      YRs=sort(as.character(unique(a$FINYEAR)))
      SIzes=colnames(a[,4:ncol(a)])
      d=matrix(ncol=length(SIzes),nrow=length(YRs))
      d.f=d.m=d
      for(s in 1:length(YRs))
      {
        ddat=subset(a,FINYEAR==YRs[s])
        ddat.f=subset(ddat,Sex=="F")
        ddat.m=subset(ddat,Sex=="M")
        
        d[s,]=colSums(ddat[,4:ncol(ddat)])
        d.f[s,]=colSums(ddat.f[,4:ncol(ddat.f)])
        d.m[s,]=colSums(ddat.m[,4:ncol(ddat.m)])
      }
      #sexes combined
      colnames(d)=SIzes
      a=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
      write.csv(a,paste("size.comp.",names(dat)[i],".csv",sep=""),row.names=F)
      #by sex
      colnames(d.f)=colnames(d.m)=SIzes
      a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.f))
      write.csv(a,paste("size.comp.fem.",names(dat)[i],".csv",sep=""),row.names=F)
      
      a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.m))
      write.csv(a,paste("size.comp.mal.",names(dat)[i],".csv",sep=""),row.names=F)
    }
    
    if(class(a)=="list" & length(a)>0)
    {
      if(Level=='annual')
      {
        a=do.call(rbind,a)
        YRs=sort(as.character(unique(a$FINYEAR)))
        SIzes=colnames(a[,4:ncol(a)])
        d=matrix(ncol=length(SIzes),nrow=length(YRs))
        d.f=d.m=d
        for(s in 1:length(YRs))
        {
          ddat=subset(a,FINYEAR==YRs[s])
          ddat.f=subset(ddat,Sex=="F")
          ddat.m=subset(ddat,Sex=="M")
          
          d[s,]=colSums(ddat[,4:ncol(ddat)])
          d.f[s,]=colSums(ddat.f[,4:ncol(ddat.f)])
          d.m[s,]=colSums(ddat.m[,4:ncol(ddat.m)])
        }
        #sexes combined
        colnames(d)=SIzes
        a=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
        write.csv(a,paste("size.comp.",Level,".",names(dat)[i],".csv",sep=""),row.names=F)
        
        #by sex
        colnames(d.f)=colnames(d.m)=SIzes
        a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.f))
        write.csv(a,paste("size.comp.fem.",Level,".",names(dat)[i],".csv",sep=""),row.names=F)
        
        a=cbind(FINYEAR=YRs,as.data.frame.matrix(d.m))
        write.csv(a,paste("size.comp.mal.",Level,".",names(dat)[i],".csv",sep=""),row.names=F)
        
      }
      
      if(Level=='annual.by.zone') 
      {
        ss=names(a)
        for(p in 1:length(ss))
        {
          b=a[[p]]
          if(nrow(b)>0)
          {
            YRs=sort(as.character(unique(b$FINYEAR)))
            SIzes=colnames(b[,4:ncol(b)])
            d=matrix(ncol=length(SIzes),nrow=length(YRs))
            d.f=d.m=d
            for(s in 1:length(YRs))
            {      
              ddat=subset(b,FINYEAR==YRs[s])
              ddat.f=subset(ddat,Sex=="F")
              ddat.m=subset(ddat,Sex=="M")
              
              d[s,]=colSums(ddat[,4:ncol(ddat)])
              d.f[s,]=colSums(ddat.f[,4:ncol(ddat.f)])
              d.m[s,]=colSums(ddat.m[,4:ncol(ddat.m)])
              
            }
            
            #sexes combined
            colnames(d)=SIzes
            b=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
            write.csv(b,paste("size.comp.annual.",ss[p],".",names(dat)[i],".csv",sep=""),row.names=F)
            
            #by sex
            colnames(d.f)=colnames(d.m)=SIzes
            b=cbind(FINYEAR=YRs,as.data.frame.matrix(d.f))
            write.csv(b,paste("size.comp.annual.fem.",ss[p],".",names(dat)[i],".csv",sep=""),row.names=F)
            
            b=cbind(FINYEAR=YRs,as.data.frame.matrix(d.m))
            write.csv(b,paste("size.comp.annual.mal.",ss[p],".",names(dat)[i],".csv",sep=""),row.names=F)
          }
          
        }
      }
    }
  }
}

#function for exporting tagging data
tagging.fn=function(dat,THESE,THOSE,type,txt)
{
  dat=dat[,match(THESE,names(dat))]
  if(type=='rel')
  {
    dat$gender=0
    dat=dat[,match(c("TG","Rel.zone","FinYear.rel","Mn.rel","gender","Age","Number"),names(dat))]
  }
  
  if(type=='rec')dat=dat[,match(c("TG","FinYear.rec","Mn.rec","Rec.zone","Number"),names(dat))]
  names(dat)=THOSE
  write.csv(dat,paste(txt,"csv",sep=""),row.names=F)
}


#Visualize average weight
fn.see.avg.wgt=function(b=Avr.wt.yr)
{
  N=1:length(unique(b$Finyear))
  SD=b$mean*(b$CV)
  plot(N,b$mean,main="",cex.main=1.25,xaxt='n',ylim=c(0,max(b$mean+SD)*1.05),
       ylab="",xlab="",pch=19,cex=3,cex.axis=1,col=tot.col,xlim=c(0,N[length(N)]+0.5))
  segments(N,b$mean,N,b$mean-SD,lwd=2,col=tot.col)
  segments(N,b$mean,N,b$mean+SD,lwd=2,col=tot.col)
  axis(1,N,F,tck=-0.015)
  axis(1,seq(1,length(N),2),b$Finyear[seq(1,length(N),2)],tck=-0.02,cex.axis=1.25)
}

fn.see.avg.wgt.zn=function(b=Avr.wt.yr.zn)
{
  zn=unique(b$zone)
  N=1:length(unique(b$Finyear))
  SD=b$mean*(b$CV)
  plot(N,ylim=c(0,max(b$mean+SD)*1.05),main="",cex.main=1.25,xaxt='n',
       ylab="",xlab="",pch=19,cex=3,cex.axis=1,col="transparent",xlim=c(0,N[length(N)]+0.5))
  
  jit=c(0,.1,.2)
  CLOS=COL.prop[match(c("West","Zone1","Zone2"),names(COL.prop))]
  
  for(x in 1:length(zn))
  {
    a=subset(b, zone==zn[x])
    N1=N+jit[x]
    SD=a$mean*(a$CV)
    
    points(N1,a$mean,ylim=c(0,max(a$mean+SD)*1.05),main=zn[x],cex.main=1.75,xaxt='n',
           ylab="",pch=19,cex=2,cex.axis=1.25,col=CLOS[x])
    segments(N1,a$mean,N1,a$mean-SD,lwd=2,col=CLOS[x])
    segments(N1,a$mean,N1,a$mean+SD,lwd=2,col=CLOS[x])
    axis(1,N,F,tck=-0.015)    
  }
  axis(1,seq(1,length(N),2),a$Finyear[seq(1,length(N),2)],tck=-0.02,cex.axis=1.25)
  legend('bottomleft',zn,pch=19,col=CLOS,pt.cex=1.5,cex=1.2,bty='n')
}

#Recode TG.zn
fn.recode=function(dat,dat1)
{
  dat=aggregate(Number~TG.zn+Rel.zone+Yr.rel+Mn.rel+Age+FinYear.rel,dat,sum)
  dat1=aggregate(Number~TG.zn+Rec.zone+Yr.rec+Mn.rec+FinYear.rec,dat1,sum)
  
  dat=dat[order(dat$FinYear.rel,dat$Mn.rel),]
  dat$TG=1:nrow(dat)
  dat1=merge(dat1,subset(dat,select=c(TG.zn,TG)),by="TG.zn",all.x=T)
  return(list(dat=dat,dat1=dat1))
}

#Create table of proportion of recaptures by year
fn.plot.prop.rec.yr=function(dat)
{
  a=aggregate(Number~FinYear.rec,dat,sum)
  a$Prop=a$Number/sum(a$Number)
  Yr1=as.numeric(substr(min(dat$FinYear.rec),start=1,stop=4))
  Yr1=1993
  Yr2=as.numeric(substr(max(dat$FinYear.rec),start=1,stop=4))
  SEQ=seq(Yr1,Yr2)
  SEQ1=seq(Yr1+1,Yr2+1)
  All.yrs=paste(SEQ,"-",substr(SEQ1,3,4),sep="")
  id=All.yrs[which(!All.yrs%in%a$FinYear.rec)]
  Add=data.frame(FinYear.rec=id,Number=0,Prop=0)
  a=rbind(a,Add)
  a=a[order(a$FinYear.rec),]
  plot(1:length(a$FinYear.rec),a$Prop,xaxt='n',ylab="",xlab="",cex.axis=1.25,
       pch=19,col=2,cex=2)
  axis(1,1:length(a$FinYear.rec),F,tck=-0.03)
  axis(1,seq(1,length(a$FinYear.rec),2),a$FinYear.rec[seq(1,length(a$FinYear.rec),2)],cex.axis=.9,tck=-0.06)
}

#Convert FL to TL 
FL_to_TL=function(FL,a,b) TL=FL*a+b

#Drop non-representative observations in TDGDLF (other fisheries are ok) 
drop.obs=function(dd,ZN,This.yr.zn) 
{
  dd$dummy=paste(dd$FINYEAR,ZN)
  dd=subset(dd,dummy%in%This.yr.zn)
  dd=dd[,-match('dummy',names(dd))]
  return(dd)
}

#back fill mesh using temporal change in 7inch upto 2009-10
fn.bck.fil=function(a,intercpt,bck.fil.yrs)
{
  a$yr=1:nrow(a)
  a$OFFSET=intercpt
  DaTt=a[bck.fil.yrs,]
  #set intercept 
  Model=lm(X165~0+yr,offset=OFFSET,DaTt)
  NEWDAT=data.frame(yr=a[1:(bck.fil.yrs[1]-1),]$yr,
                    OFFSET=unique(a$OFFSET))
  P_6.5=predict(Model,newdata = NEWDAT)
  P_7=1-P_6.5
  return(data.frame(X165=P_6.5,X178=P_7))
}

#Visualize size composition
fnx=function(x) as.numeric(substr(x,1,4))
fn.bub=function(DAT,COL,Scale,TCK,FinYrs)
{
  if(nrow(DAT)==0)
  {
    plot.new()
    legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
    return(1)
  }
  else
  {
    #aggregate data into years
    YRs=unique(DAT$FINYEAR)
    this=which(!colnames(DAT)%in%c('FINYEAR','Month','Sex'))
    SIzes=colnames(DAT[,this])
    d=matrix(ncol=length(SIzes),nrow=length(YRs))
    for(s in 1:length(YRs))
    {
      a=subset(DAT,FINYEAR==YRs[s])
      d[s,]=colSums(a[,this])
    }
    colnames(d)=SIzes
    d=cbind(FINYEAR=YRs,as.data.frame.matrix(d))
    
    #keep years with minimum observations
    drop.yr=rowSums(d[,2:ncol(d)])
    names(drop.yr)=d$FINYEAR
    drop.yr=subset(drop.yr,drop.yr<Min.obs)
    d=subset(d,!FINYEAR%in%names(drop.yr))
    
    if(nrow(d)==0)
    {
      plot.new()
      legend('center',paste("Less than",Min.obs,"observations"),
             cex=1.25,text.col='steelblue',bty='n')
      return(1)
    }else
    {
      DAT=d
      Ylabs=FinYrs
      ID=which(!FinYrs%in%DAT$FINYEAR)
      ADD=data.frame(FINYEAR=FinYrs[ID])
      if(nrow(ADD)>0)
      {
        ADD1=DAT[1:nrow(ADD),2:ncol(DAT)]
        ADD=cbind(ADD,ADD1)
        ADD[,2:ncol(ADD)]=0
        DAT=rbind(DAT,ADD)
      }
      DAT$FINYEAR=as.character(DAT$FINYEAR)
      DAT=DAT[order(DAT$FINYEAR),]
      DAT=DAT[,2:ncol(DAT)]
      x=as.numeric(colnames(DAT))
      y=1:nrow(DAT)
      z=(as.matrix(DAT))      
      n=length(y)      
      xo=outer(x,rep(1,length=length(y)))
      yo=t(outer(y,rep(1,length=length(x))))
      zo=z/ rowSums(z)*Scale   
      N=max(10,n/2)
      matplot(yo,xo,type="n",xlab="",ylab="",xaxt='n',yaxt='n')
      #abline(v=pretty(x),lty=3,col="black")
      for(i in 1:n) points(yo[,i],xo[,i],cex=zo[i,],pch=16,col=COL)
      axis(1,y,F,tck=TCK)
      axis(2,x,F,tck=TCK)
      
      if(length(Ylabs)>10)
      {
        Labx=y[seq(1,length(Ylabs),5)]
        Labx.lgn=Ylabs[seq(1,length(Ylabs),5)]
      }else
      {
        Labx=y[seq(1,length(Ylabs),1)]
        Labx.lgn=Ylabs[seq(1,length(Ylabs),1)]
      }
      
      
      axis(1,Labx,Labx.lgn,tck=TCK*2,cex.axis=1.25)
      
      axis(2,seq(x[1],x[length(x)],by=N),seq(x[1],x[length(x)],by=N),tck=TCK*2,cex.axis=1.25)
      legend('topright',paste('n=',sum(z)),bty='n')
      return(0)
    }
  }
}

#Table of number of observations and shots
fn.table.shots=function(dat,FSHRY)
{
  a=dat
  a$Number=1
  Obs=aggregate(Number~FINYEAR,a,sum)
  a$Dup=paste(a$year,a$Month,a$SHEET_NO)
  bb=a[!duplicated(a$Dup),]
  bb$Number=1
  Shots=aggregate(Number~FINYEAR,bb,sum)
  this=merge(Obs,Shots,by="FINYEAR")
  names(this)[2:3]=c("N.observations","N.shots")
  this$Species=unique(a$SPECIES)
  this$Fishery=FSHRY
  this$zone=FSHRY
  
  return(this)
}

fn.input.data=function(Name,Name.inputs,SP,Species,First.year,Last.year,Min.obs,Min.shts,
                       What.Efrt,Bin.size,Yr.assess,Dat,LH.par,remove.old.figures=FALSE,
                       HandL)   
{
  #---DATA SECTION--- 
  
  #1. Select years of data
  Used.yrs=as.numeric(c(substr(First.year,1,4),substr(Last.year,1,4)))
  Used.yrs=seq(Used.yrs[1],Used.yrs[2])
  Used.yrs=paste(Used.yrs,substr(Used.yrs+1,3,4),sep="-")

  #2. All catches
  catch=KtCh%>%filter(SPECIES%in%Species & FINYEAR%in%Used.yrs)
  catch.zone=KtCh.zone%>%filter(SPECIES%in%Species & FINYEAR%in%Used.yrs)
  
  #3. Bring in all species-specific data files
  #nm.Dat=''
  #Dat=fn.bring.all.data(path=capitalize(Name.inputs)) superseded
  if(!is.null(Dat))
  {
    nm.Dat=names(Dat)
    
    #3.1 RELATIVE ABUNDANCE    
      #3.1.1  TDGDLF  
        #monthly
    iid=nm.Dat[fn.extract.dat.perl(STRING="(?=.*annual.abundance.basecase)(?=.*relative)(?=.*monthly)",nm.Dat)]
    if(length(iid)>0)
    {
      Ab.indx.TDGDLF= Dat[match(iid,nm.Dat)]
      Nms=fn.rename.dat(x=names(Ab.indx.TDGDLF),y=c('.annual.abundance.basecase.monthly.','relative','.csv','_'))
      Nms=ifelse(Nms=='annual.abundance.basecase.monthly','all',Nms)
      names(Ab.indx.TDGDLF)=Nms
      
      Ab.indx.TDGDLF.all=Ab.indx.TDGDLF$all
      Ab.indx.TDGDLF=Ab.indx.TDGDLF[names(Ab.indx.TDGDLF) != "all"] 
      if(length(Ab.indx.TDGDLF)==0) rm(Ab.indx.TDGDLF) else
      {
        for(z in 1:length(Ab.indx.TDGDLF)) Ab.indx.TDGDLF[[z]]$zone=capitalize(names(Ab.indx.TDGDLF)[z])
        Ab.indx.TDGDLF=do.call(rbind,Ab.indx.TDGDLF)
      }
    }
        #daily
    iid=nm.Dat[fn.extract.dat.perl(STRING="(?=.*annual.abundance.basecase)(?=.*relative)(?=.*daily)",nm.Dat)]
    if(length(iid)>0)
    {
      Ab.indx.TDGDLF.daily=Dat[match(iid,nm.Dat)]
      Nms=fn.rename.dat(x=names(Ab.indx.TDGDLF.daily),y=c('.annual.abundance.basecase.daily.','relative','.csv','_'))
      Nms=ifelse(Nms=='annual.abundance.basecase.daily','all',Nms)
      names(Ab.indx.TDGDLF.daily)=Nms
      
      Ab.indx.TDGDLF.all.daily=Ab.indx.TDGDLF.daily$all
      Ab.indx.TDGDLF.daily=Ab.indx.TDGDLF.daily[names(Ab.indx.TDGDLF.daily) != "all"] 
      if(length( Ab.indx.TDGDLF.daily)==0) rm( Ab.indx.TDGDLF.daily) else
      {
        for(z in 1:length(Ab.indx.TDGDLF.daily)) Ab.indx.TDGDLF.daily[[z]]$zone=capitalize(names(Ab.indx.TDGDLF.daily)[z])
        Ab.indx.TDGDLF.daily=do.call(rbind,Ab.indx.TDGDLF.daily)
      }
    }
    
      #3.1.2  NSF 
    if('annual.abundance.NSF_relative'%in%nm.Dat) Ab.indx.NSF=Dat$annual.abundance.NSF_relative
    
      #3.1.3. Naturaliste survey
    if('Srvy.FixSt'%in%nm.Dat)
    {
      Ab.index.Srvy.FixSt=Dat$Srvy.FixSt
      Size.index.Srvy.FixSt=Dat$Srvy.FixSt_size
      Ab.index.Srvy.FixSt$CV=Ab.index.Srvy.FixSt$CV/100
      Ab.index.Srvy.FixSt$FINYEAR=paste(Ab.index.Srvy.FixSt$yr,"-",fn.subs(Ab.index.Srvy.FixSt$yr+1),sep="")
      Size.index.Srvy.FixSt$FINYEAR=paste(Size.index.Srvy.FixSt$yr,"-",fn.subs(Size.index.Srvy.FixSt$yr+1),sep="")
    }
    
      #3.1.4  Observer TDGDLF 
    if('CPUE_Observer_TDGDLF'%in%nm.Dat) Ab.indx.observer.TDGDLF=Dat$CPUE_Observer_TDGDLF
    
      #3.1.5 Pilbara trawil 
    if('CPUE_Pilbara.trawl'%in%nm.Dat) Ab.indx.PFT=Dat$CPUE_Pilbara.trawl
    
    
    #3.2 SIZE COMPOSITION
      #3.2.1 Heald 1987 
    # description: Partial length in cm obtained from Perth fish marke.
    #             Unspecified fishing method, most likely gillnets.
    # Not used, very unreliable and strange size measures and unsure about how measured
    #PL_Heald_1987=read.csv("C:/Matias/Data/Size_composition/Heald(1987).csv",stringsAsFactors=F)  
    
      #3.2.2 Stevens 1990 (citation in Simpfendorfer & Donohue 1998)
    # description: TL (fish market) or FL (measured at sea) depending on period, in cm.
    #             gillnet of 6.5 and 7 inch mesh size.
    # Not used, cannot be allocated to a zone and unsure about how measured
    #TL_FL_Stevens_1990=read.csv("C:/Matias/Data/Size_composition/Stevens_1990_size_comp_6.5_7_inch.csv",stringsAsFactors=F)    
    
      #3.2.3.TDGDLF observing programs
    #description: FL (cm); composition observed as part of different research projects on commercial gillnet vessels. 
    #           6.5 and 7 inch mesh combined (also available are data by mesh size). Souce: "Shark database"
    iid=nm.Dat[fn.extract.dat.perl(STRING="(?=.*inch.raw)",nm.Dat)]
    if(length(iid)>0)
    {
      FL.TDGDFL= Dat[match(iid,nm.Dat)]
      Nms=fn.rename.dat(x=names(FL.TDGDFL),y=c('Size_composition_','.inch.raw'))
      Nms=ifelse(Nms=='','all',Nms)
      names(FL.TDGDFL)=Nms
      if('Size_composition_Observations'%in%nm.Dat)
      {
        TDGDFL.size.numbers=Dat$Size_composition_Observations%>%filter(Method=='GN')%>%mutate(Fishery='TDGDLF')
        NSF.size.numbers=Dat$Size_composition_Observations%>%filter(Method=='LL')%>%mutate(Fishery='NSF')
        if(nrow(NSF.size.numbers)==0) rm(NSF.size.numbers)
      }
      
      if(nrow(do.call(rbind,FL.TDGDFL))==0)
      {
        rm(FL.TDGDFL,TDGDFL.size.numbers)
      }
    }
    
    #keep years with at least 10 observations from at least 10 shots by zone
    if(exists('TDGDFL.size.numbers'))
    {
      TDGDFL.size.numbers=subset(TDGDFL.size.numbers,N.observations>=Min.obs & N.shots>=Min.shts)
      This.yr.zn=with(TDGDFL.size.numbers,paste(FINYEAR,zone))
    }
    
      #3.2.4 Pilbara trawl
    #description: FL (cm); composition observed on Pilbara trawl vessels. Source: "Shark database"
    if('Size_composition_Pilbara_Trawl'%in%nm.Dat) FL_Pilbara_trawl=Dat$Size_composition_Pilbara_Trawl
    
      #3.2.5 NSF longline
    #description: FL (cm); composition observed on NSF longline vessels. Source: "Shark database"
    if('Size_composition_NSF.LONGLINE'%in%nm.Dat) FL_NSF=Dat$Size_composition_NSF.LONGLINE
    
      #3.2.6 Naturaliste survey
    #description: FL (cm); annual longline survey. Source: "Shark database"
    if('Size_composition_Survey'%in%nm.Dat) FL_Survey=Dat$Size_composition_Survey
    if('Size_composition_Survey_Observations'%in%nm.Dat)
    {
      Survey.size.numbers=Dat$Size_composition_Survey_Observations%>%
                            mutate(Fishery='Survey',zone='Survey')%>%
                            dplyr::select(FINYEAR,zone,N.shots,N.observations,Fishery)
    }
    
    #3.2.7 South Australia
    #description: FL (cm)
    if('Size_composition_Other'%in%nm.Dat) FL_Other=Dat$Size_composition_Other
    if('Size_composition_Other_Observations'%in%nm.Dat)
    {
      Other.size.numbers=Dat$Size_composition_Other_Observations%>%
        mutate(Fishery='SA',zone='SA')%>%
        dplyr::select(FINYEAR,zone,N.shots,N.observations,Fishery)
    }
    
      #3.2.8 TEPS_TDGLDF
    if(Name=="dusky shark") Size.comp.TEPS_TDGLDF=c(3.05,4,3,3.5,3.5,3.5,3.5,3,3,3,4,3)     #raw data from Comments in TDGDLF returns
    
    
    #3.3. TAGGING
      #3.3.1 Conventional
    #note: there's also info by block but very few observations at this level....
    
    #Individual based model
    if('Con_tag_Ind.based.mod'%in%nm.Dat) Rel_rec_Conv.Tag=Dat$Con_tag_Ind.based.mod
    
    #at age
    if('Con_tag_Zn.rel_Conv.Tag'%in%nm.Dat) Zn.rel_Conv.Tag=Dat$Con_tag_Zn.rel_Conv.Tag
    if('Con_tag_Zn.rec_Conv.Tag'%in%nm.Dat) Zn.rec_Conv.Tag=Dat$Con_tag_Zn.rec_Conv.Tag
    
    #at size
    #all sizes
    if('Con_tag_Zn.rel_Conv.Tag_size'%in%nm.Dat) Zn.rel_Conv.Tag_size=Dat$Con_tag_Zn.rel_Conv.Tag_size
    if('Con_tag_Zn.rec_Conv.Tag_size'%in%nm.Dat) Zn.rec_Conv.Tag_size=Dat$Con_tag_Zn.rec_Conv.Tag_size
    
    #adults and juvenlies
    if('Con_tag_Zn.rel.adul_Conv.Tag_size'%in%nm.Dat) Zn.rel_Conv.Tag_size_adu=Dat$Con_tag_Zn.rel.adul_Conv.Tag_size
    if('Con_tag_Zn.rel.juv_Conv.Tag_size'%in%nm.Dat) Zn.rel_Conv.Tag_size_juv=Dat$Con_tag_Zn.rel.juv_Conv.Tag_size
    if('Con_tag_Zn.rec.adul_Conv.Tag_size'%in%nm.Dat) Zn.rec_Conv.Tag_size_adu=Dat$Con_tag_Zn.rec.adul_Conv.Tag_size
    if('Con_tag_Zn.rec.juv_Conv.Tag_size'%in%nm.Dat) Zn.rec_Conv.Tag_size_juv=Dat$Con_tag_Zn.rec.juv_Conv.Tag_size
    if('Con_tag_Smallest.size_Conv.Tag_size'%in%nm.Dat) Smallest_size_tagged=Dat$Con_tag_Smallest.size_Conv.Tag_size
    
    
      #3.3.2 Acoustic
    #note: there's also info by block but very few observations at this level....
    
    #Taylor 2011 approach
    if('Acous.Tag_Zn.rel_Acous.Tag'%in%nm.Dat) Zn.rel_Acous.Tag=Dat$Acous.Tag_Zn.rel_Acous.Tag
    if('Acous.Tag_Zn.rec_Acous.Tag'%in%nm.Dat) Zn.rec_Acous.Tag=Dat$Acous.Tag_Zn.rec_Acous.Tag
    
    #Proportion of time approach
    if('Acous.Tag_Zn.rel_Acous.Tag.prop'%in%nm.Dat) Zn.rel_Acous.Tag.prop=Dat$Acous.Tag_Zn.rel_Acous.Tag.prop
    if('Acous.Tag_Zn.rec_Acous.Tag.prop'%in%nm.Dat) Zn.rec_Acous.Tag.prop=Dat$Acous.Tag_Zn.rec_Acous.Tag.prop
    
    #Individual based model
    iid=nm.Dat[fn.extract.dat.perl(STRING="(?=.*Acous.Tag)(?=.*Ind_based)",nm.Dat)]
    if(length(iid)>0) Indiv_based_Acous.Tag= Dat[match(iid,nm.Dat)]$`_Acous.Tag_Acous.Tag.Ind_based.csv`
    
    #Reported recaptures and releases of acoustic tagging
    if('Acous.Tag_Rep.Recap'%in%nm.Dat) Rep.Recap=Dat$Acous.Tag_Rep.Recap
    
    
    #3.4.  Age and growth data for whiskery only
    if(SP=="WH") Age.growth=read.csv(handl_OneDrive("Data/Age and growth/Simpfen.data.csv"))
    if(SP=="GM") Age.growth=read.csv(handl_OneDrive("Data/Age and growth/Terry/Gummy_Terry.csv"))
    if(SP=="BW") Age.growth=read.csv(handl_OneDrive("Data/Age and growth/Dusky.csv"))
    if(SP=="TK") Age.growth=read.csv(handl_OneDrive("Data/Age and growth/Sandbar.csv"))
    
    
    #3.5.  Standardised Mean size
    iid=nm.Dat[fn.extract.dat.perl(STRING="(?=.*annual.mean.size)",nm.Dat)]
    if(length(iid)>0)
    {
      Avr.wt.yr.zn=Dat[match(iid,nm.Dat)]
      Nms=fn.rename.dat(x=names(Avr.wt.yr.zn),y=c('annual.mean.size','relative','_'))
      Nms=ifelse(Nms=='','all',Nms)
      names(Avr.wt.yr.zn)=Nms
      Avr.wt.yr=Avr.wt.yr.zn$all
      Avr.wt.yr.zn=Avr.wt.yr.zn[names(Avr.wt.yr.zn) != "all"] 
      if(length(Avr.wt.yr.zn)==0) rm(Avr.wt.yr.zn) else
      {
        for(z in 1:length(Avr.wt.yr.zn)) Avr.wt.yr.zn[[z]]$zone=capitalize(names(Avr.wt.yr.zn)[z])
        Avr.wt.yr.zn=do.call(rbind,Avr.wt.yr.zn)
      }
    }
    
    #3.6. F from tagging
    iid=nm.Dat[fn.extract.dat.perl(STRING="Fishing.mortality.TDGDLF",nm.Dat)]
    if(length(iid)>0)
    {
      F.from.tagging=Dat[[match(iid,nm.Dat)]]
    }
    
    
    #3.7. Effort by zone (GN plus long line equivalent effort)
    #(in 1000 km gn days)
    if(What.Efrt=="km.gn.days") Eff.zn=read.csv(handl_OneDrive("Analyses/Data_outs/Annual.zone.eff.days.csv"))
    #(in 1000 km gn hours)
    if(What.Efrt=="km.gn.hours") Eff.zn=read.csv(handl_OneDrive("Analyses/Data_outs/Annual.zone.eff.hours.csv"))
    
    
    #3.8. Proportional effort by mesh size
    Mesh.prop.eff=read.csv(handl_OneDrive("Analyses/Catch and effort/mesh.proportional.effort.csv"))
    Mesh.prop.eff.West=read.csv(handl_OneDrive("Analyses/Catch and effort/mesh.proportional.effort.West.csv"))
    Mesh.prop.eff.Zn1=read.csv(handl_OneDrive("Analyses/Catch and effort/mesh.proportional.effort.Zone1.csv"))
    Mesh.prop.eff.Zn2=read.csv(handl_OneDrive("Analyses/Catch and effort/mesh.proportional.effort.Zone2.csv"))
    
    
    #3.9. Gillnet selectivity 
    if('gillnet.selectivity'%in%nm.Dat) Gillnet.selectivity=Dat$gillnet.selectivity
    if('gillnet.selectivity_len.age'%in%nm.Dat) Gillnet.selectivity_len.age=Dat$gillnet.selectivity_len.age
    
  }
    

  #----PARAMETERS SECTIONS -- 
  a.TL=LH.par$a_FL.to.TL
  b.TL=LH.par$b_FL.to.TL
  Max.TL=LH.par$Max.TL


  
  #------PROCEDURE SECTION---
  
  #1. CATCH

  #-Put all catches in list   
  catch=catch %>%
    mutate(
      Fishery = case_when(
        FishCubeCode=='WRL'~'WRL',
        FishCubeCode=='WTB'~'WTBF',
        FishCubeCode=='TEP'~'TEPS',
        FishCubeCode=='Recreational'~'Recreational',
        FishCubeCode=='SA MSF'~'SA_marine',
        FishCubeCode=='NSW fisheries'~'NSW fisheries',
        FishCubeCode=='GAB'~'GAB_trawl',
        FishCubeCode=='Indo'~'Indonesian',
        FishCubeCode=='Taiwan'~'Taiwanese',
        FishCubeCode=='Historic'~'Historic',
        FishCubeCode%in%c('JASDGDL','WCDGDL','C070','OAWC')~'TDGDLF',
        FishCubeCode%in%c('Discards_TDGDLF')~'TDGDLF (discards)',
        FishCubeCode%in%c('JANS','OANCGC','WANCS')~'NSF',   
        TRUE  ~ "Other"),
      Fishery=ifelse(Fishery%in%c('NSF')& Name%in%c("gummy shark","whiskery shark"),'TDGDLF',Fishery))
  n.fishry=sort(unique(catch$Fishery))
  iid=vector('list',length(n.fishry))
  names(iid)=n.fishry
  for(v in 1:length(n.fishry)) iid[[v]]=subset(catch,Fishery==n.fishry[v])
  catch=iid
  
  catch.zone=catch.zone %>%   
    mutate(
      Fishery = case_when(
        FishCubeCode=='WRL'~'WRL',
        FishCubeCode=='WTB'~'WTBF',
        FishCubeCode=='TEP'~'TEPS',
        FishCubeCode=='Recreational'~'Rec',
        FishCubeCode=='SA MSF'~'SA_marine',
        FishCubeCode=='NSW fisheries'~'NSW fisheries',
        FishCubeCode=='GAB'~'GAB_trawl',
        FishCubeCode=='Indo'~'Indonesian',
        FishCubeCode=='Taiwan'~'Taiwanese',
        FishCubeCode=='Historic'~'Historic',
        FishCubeCode%in%c('JASDGDL','WCDGDL','C070','OAWC')~'TDGDLF',
        FishCubeCode%in%c('Discards_TDGDLF')~'TDGDLF (discards)',
        FishCubeCode%in%c('JANS','OANCGC','WANCS')~'NSF',   
        TRUE  ~ "Other"),
      zone=ifelse(Name%in%c("gummy shark","whiskery shark") & 
                  zone%in%c('Closed','North','Joint'),'West',zone),
      Fishery=ifelse(Name%in%c("gummy shark","whiskery shark") & 
                  Fishery%in%c('NSF'),'TDGDLF',Fishery))
  n.fishry=sort(unique(catch.zone$Fishery))
  iid=vector('list',length(n.fishry))
  names(iid)=n.fishry
  for(v in 1:length(n.fishry)) iid[[v]]=subset(catch.zone,Fishery==n.fishry[v])
  catch.zone=iid
  

  #2. Put catch size composition data together (grouped by financial year and zone)
  
    #2.1. convert FL to TL 
  if(exists('FL.TDGDFL'))
  {
    for(n in 1:length(FL.TDGDFL))
    {
      FL.TDGDFL[[n]]$TL=round(FL_to_TL(FL.TDGDFL[[n]]$FL,a.TL,b.TL)) 
    }
  }
  if(exists('FL_Pilbara_trawl')) FL_Pilbara_trawl$TL=round(FL_to_TL(FL_Pilbara_trawl$FL,a.TL,b.TL))
  if(exists('FL_NSF')) FL_NSF$TL=round(FL_to_TL(FL_NSF$FL,a.TL,b.TL))
  if(exists('FL_Survey')) FL_Survey$TL=round(FL_to_TL(FL_Survey$FL,a.TL,b.TL))
  if(exists('FL_Other')) FL_Other$TL=round(FL_to_TL(FL_Other$FL,a.TL,b.TL)) 
  
    #2.2. do the aggregation
  if(exists('FL.TDGDFL'))   
  {
    #keep observation within logical size range and a Min.obs
    for(n in 1:length(FL.TDGDFL)) FL.TDGDFL[[n]]=FL.TDGDFL[[n]]%>%filter(TL<=Max.TL)
    
    TL_observers=vector('list',length(FL.TDGDFL))
    names(TL_observers)=names(FL.TDGDFL)
    dummy=do.call(rbind,FL.TDGDFL)
    if(nrow(dummy)==0) rm(FL.TDGDFL)
    if(nrow(dummy)>0)
    {
      MIN=min(dummy$TL/10,na.rm=T)
      MAX=max(dummy$TL/10,na.rm=T)
      for(n in 1:length(FL.TDGDFL))
      {
        if(nrow(FL.TDGDFL[[n]])>Min.obs)
        {
          a=fn.add.missing.size3(dat=FL.TDGDFL[[n]],
                                 Min=floor(MIN)*10,
                                 Max=ceiling(MAX)*10,
                                 interval=Bin.size)
          ag=a[,-(2:3)]%>%
            group_by(FINYEAR) %>% 
            summarise(across(where(is.numeric),sum))  
          
          ag=ag[which(rowSums(ag[,-1])>Min.obs),]
          a=a%>%filter(FINYEAR%in%ag$FINYEAR)
          TL_observers[[n]]=a
        }
      }
      TL_observers <- TL_observers[!sapply(TL_observers,is.null)]
      TL_observers <- TL_observers[sapply(TL_observers,nrow)>0]
    }

  }
  if(exists('FL_Pilbara_trawl'))
  {
    FL_Pilbara_trawl=FL_Pilbara_trawl%>%filter(TL<=Max.TL)
    if(nrow(FL_Pilbara_trawl)>=Min.obs)
    {
      Pil.trwl_observers=fn.add.missing.size3(dat=FL_Pilbara_trawl,
                                              Min=floor(min(FL_Pilbara_trawl$TL/10,na.rm=T))*10,
                                              Max=ceiling(max(FL_Pilbara_trawl$TL/10,na.rm=T))*10,
                                              interval=Bin.size)
    }else rm(FL_Pilbara_trawl)

  }
  if(exists('FL_NSF'))
  {
    FL_NSF=FL_NSF%>%filter(TL<=Max.TL)
    if(nrow(FL_NSF)>=Min.obs)
    {
      NSF_observers=fn.add.missing.size3(dat=FL_NSF,
                                         Min=floor(min(FL_NSF$TL/10,na.rm=T))*10,
                                         Max=ceiling(max(FL_NSF$TL/10,na.rm=T))*10,
                                         interval=Bin.size)
    }else rm(FL_NSF)

  }
  if(exists('FL_Survey'))
  {
    FL_Survey=FL_Survey%>%filter(TL<=Max.TL)
    if(nrow(FL_Survey)>=Min.obs)
    {
      Survey_observers=fn.add.missing.size3(dat=FL_Survey,
                                         Min=floor(min(FL_Survey$TL/10,na.rm=T))*10,
                                         Max=ceiling(max(FL_Survey$TL/10,na.rm=T))*10,
                                         interval=Bin.size)
    }else rm(FL_Survey)
    
  }
  if(exists('FL_Other'))
  {
    FL_Other=FL_Other%>%filter(TL<=Max.TL)
    if(nrow(FL_Other)>=Min.obs)
    {
      Other_observers=fn.add.missing.size3(dat=FL_Other,
                                            Min=floor(min(FL_Other$TL/10,na.rm=T))*10,
                                            Max=ceiling(max(FL_Other$TL/10,na.rm=T))*10,
                                            interval=Bin.size)
    }else rm(FL_Other)
    
  }
    
  #Wetline_WRL
  # Size.comp.Dusky.WRL=data.frame(FINYEAR="1999-00",TL=Dusky.WRL$TL)
  # Size.comp.Dusky.WRL$FINYEAR=as.character(Size.comp.Dusky.WRL$FINYEAR)
  # Size.comp.Dusky.WRL$FL=round(TL_to_FL(Dusky.WRL$TL,a.TL,b.TL))    
  # Size.comp.Dusky.WRL$Number=1
  # Size.comp.Dusky.WRL$Year=1
  # Size.comp.Dusky.WRL$Type="Fisher_measure"
  
    #2.3. Put all size composition (as TL) data in list  
  if(exists('FL.TDGDFL'))
  {
    All.size=TL_observers[grep("6.5",names(TL_observers))]
    if(length(All.size)>0)
    {
      names(All.size)=str_remove(names(All.size), ".6.5")
      for(i in 1:length(All.size)) All.size[[i]]=drop.obs(All.size[[i]],names(All.size)[i],This.yr.zn)
      names(All.size)=ifelse(names(All.size)=='West','WC',
                      ifelse(names(All.size)=='Zone1','Zn1',
                      ifelse(names(All.size)=='Zone2','Zn2',
                      NA)))
      All.size=list(TDGDLF=All.size) 
    }
    if(length(All.size)==0) rm(All.size)
    
    All.size_7=TL_observers[grep("7",names(TL_observers))]
    if(length(All.size_7)>0)
    {
      names(All.size_7)=str_remove(names(All.size_7), ".7")
      for(i in 1:length(All.size_7)) All.size_7[[i]]=drop.obs(All.size_7[[i]],names(All.size_7)[i],This.yr.zn)
      names(All.size_7)=ifelse(names(All.size_7)=='West','WC',
                        ifelse(names(All.size_7)=='Zone1','Zn1',
                        ifelse(names(All.size_7)=='Zone2','Zn2',
                        NA))) 
      All.size_7=All.size_7[sapply(All.size_7, nrow)>0]
      if(length(All.size_7)>0)
      {
        if(exists('All.size')) All.size$TDGDLF_7=All.size_7
        if(!exists('All.size')) All.size=list(TDGDLF_7=All.size_7)
      }
        
    }

  }
  if(exists('FL_Pilbara_trawl'))
  {
    if(exists('All.size')) All.size$Pilbara_trawl=Pil.trwl_observers
    if(!exists('All.size')) All.size=list(Pilbara_trawl=Pil.trwl_observers)
  }
  if(exists('FL_NSF'))
  {
    if(exists('All.size')) All.size$NSF=NSF_observers
    if(!exists('All.size')) All.size=list(NSF=NSF_observers)
  }
  if(exists('FL_Survey'))
  {
    if(exists('All.size')) All.size$Survey=Survey_observers
    if(!exists('All.size')) All.size=list(Survey=Survey_observers)
  }
  if(exists('FL_Other'))
  {
    if(exists('All.size')) All.size$Other=Other_observers
    if(!exists('All.size')) All.size=list(Other=Other_observers)
  }
  if(exists('All.size'))
  {
    All.size=All.size[sapply(All.size, length)>0]
    if(length(All.size)==0) rm(All.size)
  }
  if(exists('All.size_7')) rm(All.size_7)

  
  #3. Put abundance data together 
  if(exists('Ab.index.Srvy.FixSt')) Naturaliste.abun=Ab.index.Srvy.FixSt
  

  #4. Set working directory for outputting figures
  DiR=paste(HandL,capitalize(Name),"/",Yr.assess,"/1_Inputs/Visualise data",sep='')
  if(!file.exists(DiR))
  {
    mainDir=paste(HandL,capitalize(Name),"/",Yr.assess,sep="")
    dir.create(mainDir)
    subDir="1_Inputs"
    dir.create(file.path(mainDir,subDir))
    subDir="/1_Inputs/Visualise data"
    dir.create(file.path(mainDir,subDir))
  }
  setwd(DiR)
  
  #4.1 remove old figures from path
  if(remove.old.figures) ff=do.call(file.remove, list(list.files(DiR, full.names = TRUE)))
  
  
  #5. Years with data for conventional tagging data
  if(exists("Zn.rel_Conv.Tag"))
  {
    a=Zn.rel_Conv.Tag
    a$Finyr=with(a,ifelse(Mn.rel>6,Yr.rel+1,Yr.rel))
    c.Tag=list(TDGDLF=data.frame(FINYEAR=paste(sort(unique(a$Finyr)),"-",fn.subs(sort(unique(a$Finyr))+1),sep="")))  
    
    #Add financial year to release and recapture files
    fn.finyr.rel=function(dat)with(dat,ifelse(Mn.rel>6,Yr.rel,Yr.rel-1))
    fn.finyr.rec=function(dat)with(dat,ifelse(Mn.rec>6,Yr.rec,Yr.rec-1))
    fn.finyr=function(Finyr)paste(Finyr,"-",fn.subs(Finyr+1),sep="")
    
    Zn.rel_Conv.Tag$FinYear.rel=fn.finyr.rel(Zn.rel_Conv.Tag)
    Zn.rec_Conv.Tag$FinYear.rec=fn.finyr.rec(Zn.rec_Conv.Tag)
    Zn.rel_Conv.Tag$FinYear.rel=fn.finyr(Zn.rel_Conv.Tag$FinYear.rel)
    Zn.rec_Conv.Tag$FinYear.rec=fn.finyr(Zn.rec_Conv.Tag$FinYear.rec)
    
    #Recode TG.zn
    a=fn.recode(Zn.rel_Conv.Tag,Zn.rec_Conv.Tag)
    Zn.rel_Conv.Tag=a$dat
    Zn.rec_Conv.Tag=a$dat1
    
    #Create table of proportion of recaptures by year
    fn.fig("Proportion.tag.recaptures",2400, 1800)
    par(mfcol=c(1,1),las=1,mai=c(0.3,0.55,.1,.1),oma=c(2.25,2.25,.1,.1),mgp=c(1,.5,0))
    fn.plot.prop.rec.yr(Zn.rec_Conv.Tag)
    legend("topright",Name,bty='n',cex=1.5)
    mtext("Financial year",1,line=1,cex=1.5,outer=T)
    mtext("Proportion of recaptures",2,line=1,cex=1.5,outer=T,las=3)
    dev.off()
  }
  
  #6. Years with data for acoustic tagging data
  if(exists("Zn.rel_Acous.Tag"))
  {
    #Add financial year
    fn.finyr.rel=function(dat)with(dat,ifelse(Month.rel>6,Year.rel,Year.rel-1))
    fn.finyr.rec=function(dat)with(dat,ifelse(Month>6,Year,Year-1))
    
    Zn.rel_Acous.Tag$FinYear.rel=fn.finyr.rel(Zn.rel_Acous.Tag)
    Zn.rec_Acous.Tag$FinYear.rec=fn.finyr.rec(Zn.rec_Acous.Tag)
    Zn.rel_Acous.Tag$FinYear.rel=fn.finyr(Zn.rel_Acous.Tag$FinYear.rel)
    Zn.rec_Acous.Tag$FinYear.rec=fn.finyr(Zn.rec_Acous.Tag$FinYear.rec)
    
    a.Tag=sort(unique(Zn.rec_Acous.Tag$FinYear.rec))
    
  }
  
  #7. Add missing years to prop.eff.mesh
  #note: back fill using temporal change in 7inch upto 2009-10
  if('TDGDLF'%in%names(catch) & exists('Mesh.prop.eff'))
  {
    dum.yr=sort(as.character(unique(catch$TDGDLF%>%filter(Data.set=='Data.monthly')%>%pull(FINYEAR))))
    id=dum.yr[which(!dum.yr%in%Mesh.prop.eff$finyear)]
    add.msh.eff=Mesh.prop.eff[(nrow(Mesh.prop.eff)+1):((nrow(Mesh.prop.eff))+length(id)),]
    add.msh.eff$finyear=id
    
    Mesh.prop.eff=rbind(add.msh.eff,Mesh.prop.eff)
    Mesh.prop.eff.West=rbind(add.msh.eff,Mesh.prop.eff.West)
    Mesh.prop.eff.Zn1=rbind(add.msh.eff,Mesh.prop.eff.Zn1)
    Mesh.prop.eff.Zn2=rbind(add.msh.eff,Mesh.prop.eff.Zn2)
    
    bck.fil.yrs=match(c("2005-06","2006-07","2007-08",
                        "2008-09","2009-10"),Mesh.prop.eff$finyear)

    Mesh.prop.eff[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff,intercpt=0,bck.fil.yrs)
    Mesh.prop.eff.West[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff.West,intercpt=0,bck.fil.yrs)
    Mesh.prop.eff.Zn1[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff.Zn1,intercpt=1,bck.fil.yrs)
    Mesh.prop.eff.Zn2[1:length(id),2:3]=fn.bck.fil(Mesh.prop.eff.Zn2,intercpt=0,bck.fil.yrs)
    colnames(Mesh.prop.eff)[2:3]=colnames(Mesh.prop.eff.West)[2:3]=
      colnames(Mesh.prop.eff.Zn1)[2:3]=colnames(Mesh.prop.eff.Zn2)[2:3]=c("Mesh_6.5","Mesh_7")
    
  }
  
 
  #------ RESULTS SECTION ---

  #1. Visualize catch 
  a=do.call(rbind,catch)
  YR.span=sort(as.character(unique(a$FINYEAR)))
  rm(a)
  #YRS=YR.span
  n.plots=length(catch)+1
  WIDTH=ifelse(n.plots>5,2400,2400)
  LENGTH=ifelse(n.plots>5,2000,2400)
  fn.fig("All recons catches",WIDTH, LENGTH)
  par(las=1)
  smart.par(n.plots=n.plots,MAR=c(1.75,1.75,1.75,1.75),OMA=c(2,2,.1,.1),MGP=c(1.5,.6,0))
  fn.see.Ktch(DAT=catch,
              DAT.zone=catch.zone,
              YR.span=YR.span,
              LWD=2.5,
              TCK=-.02)
  mtext("Financial year",1,line=.75,cex=1.25,outer=T)
  mtext("Total catch (tonnes)",2,line=.4,cex=1.25,outer=T,las=3)
  dev.off()
  

  #2. Visualize size composition                           
  if(exists('All.size'))
  {
    #Size frequency
    All.size.gathered=All.size
    for(ka in 1:length(All.size.gathered))
    {
      dummy=All.size.gathered[[ka]]
      if(class(dummy)=="data.frame")
      {
        if('SEX'%in%names(dummy)) All.size.gathered[[ka]]=All.size.gathered[[ka]]%>%rename(Sex=SEX)
        All.size.gathered[[ka]]=All.size.gathered[[ka]]%>%
                                  filter(!is.na(Sex))%>%
                                  gather(Size.class,Numbers,-c(FINYEAR,Month,Sex))%>%
                                  mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
                                  group_by(year,Sex,Size.class)%>%
                                  summarise(Numbers=sum(Numbers))%>%
                                  ungroup()%>%
                                  mutate(Zone=names(All.size.gathered)[ka],
                                         Fishery=names(All.size.gathered)[ka])
      }
      if(class(dummy)=="list")
      {
        for(pp in 1:length(All.size.gathered[[ka]]))
        {
          All.size.gathered[[ka]][[pp]]=All.size.gathered[[ka]][[pp]]%>%
            gather(Size.class,Numbers,-c(FINYEAR,Month,Sex))%>%
            mutate(year=as.numeric(substr(FINYEAR,1,4)))%>%
            group_by(year,Sex,Size.class)%>%
            summarise(Numbers=sum(Numbers))%>%
            ungroup()%>%
            mutate(Zone=names(All.size.gathered[[ka]])[[pp]])
        }
        All.size.gathered[[ka]]=do.call(rbind,All.size.gathered[[ka]])%>%
          mutate(Fishery=names(All.size.gathered)[ka])
      }
    }
    All.size.gathered=do.call(rbind,All.size.gathered)%>%
      mutate(Size.class=as.numeric(Size.class))
    
    colfunc <- colorRampPalette(c("red","orange","forestgreen","dodgerblue4"))
    Fem.linf.TL=with(LH.par,FL_inf*a_FL.to.TL+b_FL.to.TL)
    Male.linf.TL=with(LH.par,male_FL_inf*a_FL.to.TL+b_FL.to.TL)
    TL.ratio=paste('TL_50%mat : TL_inf ratio =',round(100*LH.par$TL.50.mat/Fem.linf.TL),"%",sep='')
    
    tdgdlf.size=All.size.gathered%>%filter(grepl('TDGDLF',Fishery))%>%ungroup()
    if(nrow(tdgdlf.size)>0)
    {
      nKl=length(unique(tdgdlf.size$year))
      tdgdlf.size%>%
        mutate(Fishery=case_when(Fishery=='TDGDLF'~'6.5 inch mesh',
                                 Fishery=='TDGDLF_7'~'7 inch mesh'),
               year=factor(year,levels=sort(unique(year))),
               Zone.Fishery=paste(Zone,Fishery,sep=' & '))%>%
        ggplot(aes(Size.class,Numbers,color=year))+
        geom_line(aes(linetype=Sex),size=1.05)+
        facet_wrap(~Zone.Fishery,scales='free_y')+
        ylab("Frequency")+xlab("TL size class (cm)")+
        theme_PA(leg.siz=14,axs.t.siz=12,axs.T.siz=20,strx.siz=16)+
        scale_color_manual(values=colfunc(nKl))+
        geom_vline(xintercept=LH.par$TL.50.mat,color='deeppink')+
        annotate("text",LH.par$TL.50.mat,0,label="50% maturity",angle='90',size=5,vjust=-0.1,hjust=-0.1,color='deeppink')+
        geom_vline(xintercept=Fem.linf.TL)+
        annotate("text",Fem.linf.TL,0,label="Linf",angle='90',size=5,vjust=-0.1,hjust=-0.1)+
        geom_vline(xintercept=Male.linf.TL,linetype="dashed")+
        xlim(min(tdgdlf.size$Size.class),max(tdgdlf.size$Size.class,Fem.linf.TL))+
        labs(caption = TL.ratio)+theme(plot.caption = element_text(size = 15))
 
      ggsave("Size.comp.TDGDLF_dist.tiff",width = 14,height = 10,compression = "lzw")  
    }
    
    nsf.size=All.size.gathered%>%filter(grepl('NSF',Fishery))%>%ungroup()
    if(nrow(nsf.size)>0)
    {
      nKl=length(unique(nsf.size$year))
      nsf.size%>%
        mutate(year=factor(year,levels=sort(unique(year))))%>%
        ggplot(aes(Size.class,Numbers,color=year))+
        geom_line(aes(linetype=Sex),size=1.05)+
        ylab("Frequency")+xlab("TL size class (cm)")+
        theme_PA(leg.siz=14,axs.t.siz=12,axs.T.siz=20,strx.siz=18)+
        scale_color_manual(values=colfunc(nKl))+
        geom_vline(xintercept=LH.par$TL.50.mat,color='deeppink')+
        annotate("text",LH.par$TL.50.mat,0,label="50% maturity",angle='90',size=5,vjust=-0.1,hjust=-0.1,color='deeppink')+
        geom_vline(xintercept=Fem.linf.TL)+
        annotate("text",Fem.linf.TL,0,label="Linf",angle='90',size=5,vjust=-0.1,hjust=-0.1)+
        geom_vline(xintercept=Male.linf.TL,linetype="dashed")+
        xlim(min(nsf.size$Size.class),max(nsf.size$Size.class,Fem.linf.TL))+
        labs(caption = TL.ratio)+theme(plot.caption = element_text(size = 15))
      
      ggsave("Size.comp.NSF_dist.tiff",width = 10,height = 10,compression = "lzw")
    }
    
    Survey.size=All.size.gathered%>%filter(grepl('Survey',Fishery))%>%ungroup()
    if(nrow(Survey.size)>0)
    {
      nKl=length(unique(Survey.size$year))
      Survey.size%>%
        mutate(year=factor(year,levels=sort(unique(year))))%>%
        ggplot(aes(Size.class,Numbers,color=year))+
        geom_line(aes(linetype=Sex),size=1.05)+
        ylab("Frequency")+xlab("TL size class (cm)")+
        theme_PA(leg.siz=14,axs.t.siz=12,axs.T.siz=20,strx.siz=18)+
        scale_color_manual(values=colfunc(nKl))+
        geom_vline(xintercept=LH.par$TL.50.mat,color='deeppink')+
        annotate("text",LH.par$TL.50.mat,0,label="50% maturity",angle='90',size=5,vjust=-0.1,hjust=-0.1,color='deeppink')+
        geom_vline(xintercept=Fem.linf.TL)+
        annotate("text",Fem.linf.TL,0,label="Linf",angle='90',size=5,vjust=-0.1,hjust=-0.1)+
        geom_vline(xintercept=Male.linf.TL,linetype="dashed")+
        xlim(min(Survey.size$Size.class),max(Survey.size$Size.class,Fem.linf.TL))+
        labs(caption = TL.ratio)+theme(plot.caption = element_text(size = 15))
      
      ggsave("Size.comp.Survey_dist.tiff",width = 10,height = 10,compression = "lzw")
    }
    
    pilbara.size=All.size.gathered%>%filter(grepl('Pilbara_trawl',Fishery))%>%ungroup()
    if(nrow(pilbara.size)>0)
    {
      nKl=length(unique(pilbara.size$year))
      pilbara.size%>%
        mutate(year=factor(year,levels=sort(unique(year))))%>%
        ggplot(aes(Size.class,Numbers,color=year))+
        geom_line(aes(linetype=Sex),size=1.05)+
        ylab("Frequency")+xlab("TL size class (cm)")+
        theme_PA(leg.siz=14,axs.t.siz=12,axs.T.siz=20,strx.siz=18)+
        scale_color_manual(values=colfunc(nKl))+
        geom_vline(xintercept=LH.par$TL.50.mat,color='deeppink')+
        annotate("text",LH.par$TL.50.mat,0,label="50% maturity",angle='90',size=5,vjust=-0.1,hjust=-0.1,color='deeppink')+
        geom_vline(xintercept=Fem.linf.TL)+
        annotate("text",Fem.linf.TL,0,label="Linf",angle='90',size=5,vjust=-0.1,hjust=-0.1)+
        geom_vline(xintercept=Male.linf.TL,linetype="dashed")+
        xlim(min(pilbara.size$Size.class),max(pilbara.size$Size.class,Fem.linf.TL))+
        labs(caption = TL.ratio)+theme(plot.caption = element_text(size = 15))
      
      ggsave("Size.comp.Pilbara.trawl_dist.tiff",width = 10,height = 10,compression = "lzw")
    }
    
    other.size=All.size.gathered%>%filter(grepl('Other',Fishery))%>%ungroup()
    if(nrow(other.size)>0)
    {
      nKl=length(unique(other.size$year))
      other.size%>%
        mutate(year=factor(year,levels=sort(unique(year))))%>%
        ggplot(aes(Size.class,Numbers,color=year))+
        geom_line(aes(linetype=Sex),size=1.05)+
        ylab("Frequency")+xlab("TL size class (cm)")+
        theme_PA(leg.siz=14,axs.t.siz=12,axs.T.siz=20,strx.siz=18)+
        scale_color_manual(values=colfunc(nKl))+
        geom_vline(xintercept=LH.par$TL.50.mat)+
        annotate("text",LH.par$TL.50.mat,0,label="50% maturity",size=5,vjust=1.5)
      ggsave("Size.comp.Other_dist.tiff",width = 10,height = 10,compression = "lzw")
    }
    
    
    #bubble plots
      #2.1 TDGDLF size comp of reported catch
    if(sum(c('TDGDLF','TDGDLF_7')%in%names(All.size))>0)
    {
      a=NULL
      b=NULL
      if('TDGDLF'%in%names(All.size)) a=do.call(rbind,All.size$TDGDLF)
      if('TDGDLF_7'%in%names(All.size)) b=do.call(rbind,All.size$TDGDLF_7)
      FinYrs=rbind(a,b)
      FinYrs=sort(as.character(unique(FinYrs%>%pull(FINYEAR))))
      FinYrs=paste(fnx(min(FinYrs)):(fnx(max(FinYrs))),"-",
                   substr((fnx(min(FinYrs))+1):(fnx(max(FinYrs))+1),3,4),sep="")
      
      if(sum(c('TDGDLF','TDGDLF_7')%in%names(All.size))==2)
      {
        fn.fig("Size.comp.TDGDLF",2400, 2000)
        par(mfcol=c(3,2),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.25,1,1),mgp=c(1,.85,0))
        #-6.5
        DROP_6.5=vector('list',length(All.size$TDGDLF))
        for(e in 1:length(All.size$TDGDLF))
        {
          DROP_6.5[[e]]=fn.bub(DAT=All.size$TDGDLF[[e]],COL='steelblue',Scale=7.5,TCK=-.015,FinYrs)
          if(e==1) mtext('6.5 inch mesh',3)
        }
        if(e<3)
        {
          misin=3-length(All.size$TDGDLF)
          for(m in 1:misin)
          {
            plot.new()
            legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
          }
        }
        
        #-7
        DROP_7=vector('list',length(All.size$TDGDLF_7))
        for(e in 1:length(All.size$TDGDLF_7))
        {
          DROP_7[[e]]=fn.bub(DAT=All.size$TDGDLF_7[[e]],COL='steelblue',Scale=7.5,TCK=-.015,FinYrs)
          if(e==1) mtext('7 inch mesh',3)
          mtext(names(All.size$TDGDLF_7)[e],4,las=3,line=.35)
        }
        if(e<3)
        {
          misin=3-length(All.size$TDGDLF_7)
          for(m in 1:misin)
          {
            plot.new()
            legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
          }
        }
        mtext("Financial Year",1,cex=1.35,line=0.75,outer=T)
        mtext("Total length class (cm)",2,cex=1.35,line=0.5,las=3,outer=T)
        dev.off()
        
        if(sum(c(unlist(DROP_6.5),unlist(DROP_7)))==length(All.size$TDGDLF)+length(All.size$TDGDLF_7))
        {
          rm(All.size,TDGDFL.size.numbers)
        }
        
      }
      if('TDGDLF'%in%names(All.size)& ! 'TDGDLF_7'%in%names(All.size))
      {
        fn.fig("Size.comp.TDGDLF",2400, 2000)
        par(mfcol=c(3,1),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.25,1,1),mgp=c(1,.85,0))
        #-6.5
        DROP_6.5=vector('list',length(All.size$TDGDLF))
        for(e in 1:length(All.size$TDGDLF))
        {
          DROP_6.5[[e]]=fn.bub(DAT=All.size$TDGDLF[[e]],COL='steelblue',Scale=7.5,TCK=-.015,FinYrs)
          if(e==1) mtext('6.5 inch mesh',3)
        }
        if(e<3)
        {
          misin=3-length(All.size$TDGDLF)
          for(m in 1:length(misin))
          {
            plot.new()
            legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
          }
        }
         mtext("Financial Year",1,cex=1.35,line=0.75,outer=T)
        mtext("Total length class (cm)",2,cex=1.35,line=0.5,las=3,outer=T)
        dev.off()
        
        if(sum(unlist(DROP_6.5))==length(All.size$TDGDLF))
        {
          rm(All.size,TDGDFL.size.numbers)
        }
        
      }
      if(!'TDGDLF'%in%names(All.size)&  'TDGDLF_7'%in%names(All.size))
      {
        fn.fig("Size.comp.TDGDLF",2400, 2000)
        par(mfcol=c(3,1),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.25,1,1),mgp=c(1,.85,0))
        #-7
        DROP_7=vector('list',length(All.size$TDGDLF_7))
        for(e in 1:length(All.size$TDGDLF_7))
        {
          DROP_7[[e]]=fn.bub(DAT=All.size$TDGDLF_7[[e]],COL='steelblue',Scale=7.5,TCK=-.015,FinYrs)
          if(e==1) mtext('7 inch mesh',3)
          mtext(names(All.size$TDGDLF_7)[e],4,las=3,line=.35)
        }
        if(e<3)
        {
          misin=3-length(All.size$TDGDLF)
          for(m in 1:length(misin))
          {
            plot.new()
            legend('center',paste("No observations"),cex=1.25,text.col='steelblue',bty='n')
          }
        }
        mtext("Financial Year",1,cex=1.35,line=0.75,outer=T)
        mtext("Total length class (cm)",2,cex=1.35,line=0.5,las=3,outer=T)
        mtext("Financial Year",1,cex=1.35,line=0.75,outer=T)
        mtext("Total length class (cm)",2,cex=1.35,line=0.5,las=3,outer=T)
        dev.off()
        
        if(sum(unlist(DROP_7))==length(All.size$TDGDLF_7))
        {
          rm(All.size,TDGDFL.size.numbers)
        }
      } 
    }
    
      #2.2 NSF size comp of reported catch
    if('NSF'%in%names(All.size))
    {
      FinYrs=sort(as.character(unique(All.size$NSF$FINYEAR)))
      FinYrs=paste(fnx(min(FinYrs)):(fnx(max(FinYrs))),"-",
                   substr((fnx(min(FinYrs))+1):(fnx(max(FinYrs))+1),3,4),sep="")
      fn.fig("Size.comp.NSF",2400, 2000)
      par(mfcol=c(1,1),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.5,1,1),mgp=c(1,.65,0))
      fn.bub(DAT=All.size$NSF,COL='steelblue',Scale=7.5,TCK=-.01,FinYrs)
      mtext("Financial Year",1,cex=1.75,line=0.75,outer=T)
      mtext("Total length class (cm)",2,cex=1.75,line=0.95,las=3,outer=T)
      dev.off()
    }
    
      #2.3 Pilbara trawl size comp of reported catch
    if('Pilbara_trawl'%in%names(All.size))
    {
      FinYrs=sort(as.character(unique(All.size$Pilbara_trawl$FINYEAR)))
      FinYrs=paste(fnx(min(FinYrs)):(fnx(max(FinYrs))),"-",
                   substr((fnx(min(FinYrs))+1):(fnx(max(FinYrs))+1),3,4),sep="")
      fn.fig("Size.comp.Pilbara.trawl",2400, 2000)
      par(mfcol=c(1,1),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.5,1,1),mgp=c(1,.65,0))
      fn.bub(DAT=All.size$Pilbara_trawl,COL='steelblue',Scale=7.5,TCK=-.01,FinYrs)
      mtext("Financial Year",1,cex=1.75,line=0.75,outer=T)
      mtext("Total length class (cm)",2,cex=1.75,line=0.95,las=3,outer=T)
      dev.off()
    }
    
    #2.4 Other size comp of reported catch
    if('Other'%in%names(All.size))
    {
      FinYrs=sort(as.character(unique(All.size$Other$FINYEAR)))
      FinYrs=paste(fnx(min(FinYrs)):(fnx(max(FinYrs))),"-",
                   substr((fnx(min(FinYrs))+1):(fnx(max(FinYrs))+1),3,4),sep="")
      fn.fig("Size.comp.Other",2400, 2000)
      par(mfcol=c(1,1),las=1,mai=c(0.3,0.35,.1,.1),oma=c(2.25,2.5,1,1),mgp=c(1,.65,0))
      fn.bub(DAT=All.size$Other,COL='steelblue',Scale=7.5,TCK=-.01,FinYrs)
      mtext("Financial Year",1,cex=1.75,line=0.75,outer=T)
      mtext("Total length class (cm)",2,cex=1.75,line=0.95,las=3,outer=T)
      dev.off()
    }
  }

    #Table of number of observations and shots
  if(exists('TDGDFL.size.numbers'))
  {
    Size.numbers=TDGDFL.size.numbers%>%dplyr::select(Fishery,FINYEAR,zone,N.observations)
    Size.shots=TDGDFL.size.numbers%>%dplyr::select(Fishery,FINYEAR,zone,N.shots)
    
    if(exists('NSF.size.numbers'))
    {
      Size.numbers.NSF=NSF.size.numbers%>%dplyr::select(Fishery,FINYEAR,zone,N.observations)
      Size.shots.NSF=NSF.size.numbers%>%dplyr::select(Fishery,FINYEAR,zone,N.shots)
      
      Size.numbers=rbind(Size.numbers,Size.numbers.NSF)
      Size.shots=rbind(Size.shots,Size.shots.NSF)
    }
    if(exists('FL_Pilbara_trawl'))  
    {
      Pilbara_trawl.obs=fn.table.shots(FL_Pilbara_trawl,FSHRY="Pilbara trawl")
      Pilbara_trawl.size.numbers=Pilbara_trawl.obs%>%dplyr::select(Fishery,FINYEAR,zone,N.observations)
      Pilbara_trawl.size.shots=Pilbara_trawl.obs%>%dplyr::select(Fishery,FINYEAR,zone,N.shots)
      
      Size.numbers=rbind(Size.numbers,Pilbara_trawl.size.numbers)
      Size.shots=rbind(Size.shots,Pilbara_trawl.size.shots)
    }
    if(exists('Survey.size.numbers'))
    {
      Size.numbers.Survey=Survey.size.numbers%>%dplyr::select(Fishery,FINYEAR,zone,N.observations)
      Size.shots.Survey=Survey.size.numbers%>%dplyr::select(Fishery,FINYEAR,zone,N.shots)
      
      Size.numbers=rbind(Size.numbers,Size.numbers.Survey)
      Size.shots=rbind(Size.shots,Size.shots.Survey)
    }
    if(exists('Other.size.numbers'))
    {
      Size.numbers.Other=Other.size.numbers%>%dplyr::select(Fishery,FINYEAR,zone,N.observations)
      Size.shots.Other=Other.size.numbers%>%dplyr::select(Fishery,FINYEAR,zone,N.shots)
      
      Size.numbers=rbind(Size.numbers,Size.numbers.Other)
      Size.shots=rbind(Size.shots,Size.shots.Other)
    }
    
    Numbers.SF=Size.numbers%>%
                    mutate(Fishery.zone=paste(Fishery,zone,sep='-'))%>%
                    dplyr::select(-c(Fishery,zone))%>%
                    spread(Fishery.zone,N.observations,fill='')
    Shots.SF=Size.shots%>%
      mutate(Fishery.zone=paste(Fishery,zone,sep='-'))%>%
      dplyr::select(-c(Fishery,zone))%>%
      spread(Fishery.zone,N.shots,fill='')
 

    #create nice table 
    fn.word.table(TBL=Numbers.SF,Doc.nm="Size.comp.n.observations")
    fn.word.table(TBL=Shots.SF,Doc.nm="Size.comp.n.shots")
  }

 
  #3. Visualize mean weights 
  if(exists("Avr.wt.yr"))
  {

    fn.fig("Avg.wgt",2000,2000)
    par(mfcol=c(1,1),las=1,mai=c(0.45,0.35,.1,.15),oma=c(1,2.25,.1,.1),mgp=c(1,.65,0))
    fn.see.avg.wgt(b=Avr.wt.yr)
    mtext("Relative live weight",2,line=.5,cex=1.5,las=3,outer=T)
    mtext("Financial Year",1,cex=1.5,line=0,outer=T)
    dev.off()
  }
  if(exists("Avr.wt.yr.zn"))
  {
   fn.fig("Avg.wgt.zn",2000,2000)
   par(mfcol=c(1,1),las=1,mai=c(0.45,0.35,.1,.15),oma=c(1,2.25,.1,.1),mgp=c(1,.65,0))
   fn.see.avg.wgt.zn(b=Avr.wt.yr.zn)
   mtext("Relative live weight",2,line=.5,cex=1.5,las=3,outer=T)
   mtext("Financial Year",1,cex=1.5,line=0,outer=T)
   dev.off()
 }

  
  #4. Visualize F from tagging 
  if(exists("F.from.tagging"))
  {
    
    fn.fig("F from tagging_TDGDLF",2000,2000)
    par(mfcol=c(1,1),las=1,mai=c(0.45,0.35,.1,.15),oma=c(1,2.25,.1,.1),mgp=c(1,.65,0))
    SD=F.from.tagging$Mean*F.from.tagging$CV
    with(F.from.tagging,plot(Finyear,Mean,main="",cex.main=1.25,xaxt='n',ylim=c(0,max(Mean+SD)*1.05),
         ylab="",xlab="",pch=19,cex=3,cex.axis=1,col=tot.col))
    with(F.from.tagging,segments(Finyear,Mean+SD,Finyear,Mean-SD,lwd=2,col=tot.col))
    LABs=with(F.from.tagging,paste(Finyear,substr(Finyear+1,3,4),sep='-'))
    axis(1,F.from.tagging$Finyear,LABs,tck=-0.015)
    mtext("F (yr-1)",2,line=.5,cex=1.5,las=3,outer=T)
    mtext("Financial Year",1,cex=1.5,line=0,outer=T)
    dev.off()
  }

  
  #5. Select effort years
  if(exists("Eff.zn")) Eff.zn=subset(Eff.zn,FINYEAR%in%YR.span)
  if(exists("Eff.total")) Eff.total=data.frame(FINYEAR=Eff.zn$FINYEAR,Total=Eff.zn$West+Eff.zn$Zone1+Eff.zn$Zone2)
  
  
  #6. Visualize data availability
  cpue_tdgdlf.monthly=NULL
  cpue_tdgdlf.daily=NULL
  if(exists('Ab.indx.TDGDLF.all')) cpue_tdgdlf.monthly=Ab.indx.TDGDLF.all
  if(exists('Ab.indx.TDGDLF.all.daily')) cpue_tdgdlf.daily=Ab.indx.TDGDLF.all.daily
  if(!is.null(cpue_tdgdlf.monthly)|!is.null(cpue_tdgdlf.daily))
  {
    Abun=list(TDGLDF=rbind(cpue_tdgdlf.monthly,cpue_tdgdlf.daily)) 
  }
  if(exists('Ab.indx.NSF'))
  {
    if(exists('Abun'))Abun$NSF=Ab.indx.NSF
    if(!exists('Abun'))Abun=list(NSF=Ab.indx.NSF)
  }
  if(exists('Naturaliste.abun'))
  {
    if(exists('Abun')) Abun$Survey=Naturaliste.abun
    if(!exists('Abun')) Abun=list(Survey=Naturaliste.abun)
  }
  if(exists('Ab.indx.observer.TDGDLF'))
  {
    if(exists('Abun'))Abun$observer.TDGDLF=Ab.indx.observer.TDGDLF
    if(!exists('Abun'))Abun=list(observer.TDGDLF=Ab.indx.observer.TDGDLF)
  }
  if(exists('Ab.indx.PFT'))
  {
    if(exists('Abun'))Abun$PFT=Ab.indx.PFT
    if(!exists('Abun'))Abun=list(PFT=Ab.indx.PFT)
  }
  d.list=list(Catch=catch) 
  if(exists('All.size'))
  {
    if(sum(c('TDGDLF','TDGDLF_7')%in%names(All.size))>0)
    {
      d.list$'Catch size comp (TDGDLF)'=list(FINYEAR=unique(do.call(rbind,All.size$TDGDLF)$FINYEAR))
    }
    if(sum(c('NSF')%in%names(All.size))>0)
    {
      d.list$'Catch size comp (NSF)'=list(FINYEAR=unique(All.size$NSF$FINYEAR))
    }
    if(sum(c('Pilbara_trawl')%in%names(All.size))>0)
    {
      d.list$'Catch size comp (PFT)'=list(FINYEAR=unique(All.size$Pilbara_trawl$FINYEAR))
    }
    if(sum(c('Other')%in%names(All.size))>0)
    {
      d.list$'Catch size comp (Other)'=list(FINYEAR=unique(All.size$Other$FINYEAR))
    }
    
  }
  if(exists('Avr.wt.yr'))d.list$'Catch mean wt (TDGDLF)'=list(Avr.wt.yr%>%distinct(Finyear))
  if(exists('F.from.tagging'))d.list$'Fishing mortality (TDGDLF)'=list(F.from.tagging%>%
                                                                         mutate(Finyear=paste(Finyear,substr(Finyear+1,3,4),sep='-'))%>%
                                                                         distinct(Finyear))
  if(exists('Abun'))
  {
    d.list$Abundance=Abun
    if(length(d.list$Abundance)==1) names(d.list)[match('Abundance',names(d.list))]=paste("Abundance (",names(d.list$Abundance),")",sep='')
  }
  if(exists('Size.index.Srvy.FixSt'))d.list$'Survey mean size'=list(FINYEAR=Size.index.Srvy.FixSt$FINYEAR)
  if(exists('c.Tag'))d.list$'Conv tag'=c.Tag
  if(exists('a.Tag'))d.list$'Acous tag'=list(FINYEAR=a.Tag)
  Left.mar=9
  if(exists('Ab.indx.observer.TDGDLF')) Left.mar=11
  fn.fig("avail.dat",2400,1400)
  par(mar=c(3,2,.1,.1),oma=c(.1,Left.mar,.1,1),las=1,mgp=c(1.5,.6,0))
  visualize.dat(d.list,YR.span)   
  dev.off()
  
  
  #7. Table of released individuals with conventional tags 
  if(exists('Zn.rel_Conv.Tag_size_adu'))
  {
    TBL.conv.rel=merge(Zn.rel_Conv.Tag_size_adu[,-match("TG.zn",names(Zn.rel_Conv.Tag_size_adu))],
                       Zn.rel_Conv.Tag_size_juv[,-match("TG.zn",names(Zn.rel_Conv.Tag_size_juv))],
                       by=c("Rel.zone","Yr.rel"),all=T)
    names(TBL.conv.rel)[match(c("Number.x","Number.y"),
                              names(TBL.conv.rel))]=c("Adult.N","Juvenile.N")
    TBL.conv.rel=TBL.conv.rel[,-match(c("size.group.x",
                                        "size.group.y"),names(TBL.conv.rel))]
    TBL.conv.rel[is.na(TBL.conv.rel)]=0
    write.csv(TBL.conv.rel,"TBL.conv.releases.csv",row.names=F)
    fn.word.table(TBL=TBL.conv.rel,Doc.nm="TBL.conv.releases")
  }
  
  
  #8. Plot acoustic tag recaptures by year
  if(exists('Rep.Recap'))
  {
    Rep.Recap=subset(Rep.Recap,!is.na(year.rel))
    i.yrs=sort(unique(Rep.Recap$year.rec))
    Y.rng=range(c(Rep.Recap$RELLATDECDEG,Rep.Recap$RECLATDECDEG),na.rm=T)
    X.rng=range(c(Rep.Recap$RELLNGDECDEG,Rep.Recap$RECLNGDECDEG),na.rm=T)
    
    fn.fig("Acoustic.tags.rel.rec.map",2400,2400)
    par(mfcol=c(ceiling(length(i.yrs)/2),2),las=1,mai=c(1,1,.1,.1),oma=c(.1,.1,.1,.1))
    dummy=Rep.Recap
    for(i.yr in 1:length(i.yrs))
    {
      recs=subset(dummy,year.rec==i.yrs[i.yr])
      na.rcs=subset(recs,is.na(RECLATDECDEG))
      if(nrow(recs)>0)
      {
        rels=subset(dummy,year.rel<=i.yrs[i.yr])
        plot(rels$RELLNGDECDEG,rels$RELLATDECDEG,ylim=Y.rng,xlim=X.rng,ylab="",xlab="",pch=19)
        points(recs$RECLNGDECDEG,recs$RECLATDECDEG,col=2,pch=4,cex=1.75)
        if(nrow(na.rcs)==0) Rc.lg=paste(nrow(recs),"recaptured")
        if(nrow(na.rcs)>0) Rc.lg=paste(nrow(recs),"recaptured (",nrow(na.rcs),"without rec. info)")
        legend('top',c(paste(nrow(rels),"with tags"),Rc.lg),bty='n',col=c(1,2),
               text.col=c(1,2),pch=c(19,4),title=paste(i.yrs[i.yr]),cex=1.2)
        dummy=subset(dummy,!ATAG_NO%in%recs$ATAG_NO)
      }
      
    }
    mtext("Longitude",1,-2,outer=T,cex=1.75)
    mtext("Latitude",2,-2,las=3,outer=T,cex=1.75)
    dev.off()
  }
 

  #9. Export data for population dynamics modelling ---
  HandL=handl_OneDrive("Data/Population dynamics/Data inputs for models/")
  DiR=paste(HandL,str_remove(capitalize(Name), ' shark'),"/",Yr.assess,sep='')
  if(!file.exists(DiR)) dir.create(DiR)
  
  ff=do.call(file.remove, list(list.files(DiR, full.names = TRUE))) #remove all old files first
  setwd(DiR)

  #Catch
  fn.agg.at.level.and.exprt(DAT=catch,
                            Level='annual',
                            VAR="ktch",
                            SOURCE=names(catch))
  fn.agg.at.level.and.exprt(DAT=catch.zone,
                            Level='annual.by.zone',
                            VAR="ktch",
                            SOURCE=names(catch.zone))   

  #Effort
  if(exists("Eff.zn")) write.csv(Eff.zn,"effort.annual.by.zone.TDGDLF.csv",row.names=F) 
  if(exists("Eff.total")) write.csv(Eff.total,"effort.annual.TDGDLF.csv",row.names=F)

  #F from tagging
  #TDGDLF
  if(exists('F.from.tagging'))write.csv(F.from.tagging,"F.from.tagging_TDGDLF.csv",row.names=F) 
  
  #avg weight 
    #TDGDLF
  if(exists('Avr.wt.yr.zn'))write.csv(Avr.wt.yr.zn,"ktch.avg.weight.annual.by.zone.csv",row.names=F) 
  if(exists('Avr.wt.yr')) write.csv(Avr.wt.yr,"ktch.avg.weight.annual.csv",row.names=F)
    #Naturaliste survey 
  if(exists('Size.index.Srvy.FixSt')) write.csv(Size.index.Srvy.FixSt,"size.annual.survey.csv",row.names=F)
  
  #TDGLDF cpue
    #folly
  if(exists("Ab.folly.TDGDLF.all")) write.csv(Ab.folly.TDGDLF.all,"cpue.annual.TDGDLF.folly.csv",row.names=F)
    #by zone      
  if(exists('Ab.indx.TDGDLF')) write.csv(Ab.indx.TDGDLF,"cpue.annual.by.zone.TDGDLF.csv",row.names=F)
  if(exists('Ab.indx.TDGDLF.daily')) write.csv(Ab.indx.TDGDLF.daily,"cpue.annual.by.zone.TDGDLF.daily.csv",row.names=F)
    #zones combined
  if(exists('Ab.indx.TDGDLF.all')) write.csv(Ab.indx.TDGDLF.all,"cpue.annual.TDGDLF.csv",row.names=F)
  if(exists('Ab.indx.TDGDLF.all.daily')) write.csv(Ab.indx.TDGDLF.all.daily,"cpue.annual.TDGDLF.daily.csv",row.names=F)

  #Naturaliste survey     
  if(exists('Naturaliste.abun')) write.csv(Naturaliste.abun,"cpue.annual.survey.csv",row.names=F)

  #NSF cpue     
  if(exists('Ab.indx.NSF')) write.csv(Ab.indx.NSF,"cpue.annual.NSF.csv",row.names=F)
 
  #Observer TDGDLF cpue     
  if(exists('Ab.indx.observer.TDGDLF')) write.csv(Ab.indx.observer.TDGDLF,"cpue.annual.observer.TDGDLF.csv",row.names=F)
  
  #Pilbara fish trawl cpue     
  if(exists('Ab.indx.PFT')) write.csv(Ab.indx.PFT,"cpue.annual.PFT.csv",row.names=F)
  
  #Size composition 
  if(exists('All.size'))
  {
    fn.exp.size(dat=All.size,Level='annual')
    fn.exp.size(dat=All.size,Level='annual.by.zone')
    
    if(exists('Numbers.SF')) write.csv(Numbers.SF,"Numbers_size_frequency.csv",row.names=F)
    if(exists('Shots.SF')) write.csv(Shots.SF,"Shots_size_frequency.csv",row.names=F)
  }
  
  #Conventional tagging
    #1. Individual-based model
  if(exists('Rel_rec_Conv.Tag')) write.csv(Rel_rec_Conv.Tag,"Ind_based_model.csv",row.names=F)
    #2. At age
  #note: this function puts in SS3 format
  if(exists('Zn.rel_Conv.Tag'))
  {
    #releases
    these=c("TG","Rel.zone","FinYear.rel","Mn.rel","Age","Number")
    those=c("TG","area","yr","season","gender","Age","Nrelease")  
    tagging.fn(Zn.rel_Conv.Tag,these,those,'rel','conv.tag.rel.')
    
    #recaptures
    these=c("TG","Rec.zone","FinYear.rec","Mn.rec","Number")
    those=c("TG","year","season","area","Number")
    tagging.fn(Zn.rec_Conv.Tag,these,those,'rec',"conv.tag.reca.")
  }
    #3. At size
      #all
  if(exists('Zn.rel_Conv.Tag_size')) write.csv(Zn.rel_Conv.Tag_size,"Zn.rel_Conv.Tag_size.csv",row.names=F)
  if(exists('Zn.rec_Conv.Tag_size')) write.csv(Zn.rec_Conv.Tag_size,"Zn.rec_Conv.Tag_size.csv",row.names=F)
      #juveniles and adults
  if(exists('Zn.rel_Conv.Tag_size_adu')) write.csv(Zn.rel_Conv.Tag_size_adu,"Zn.rel_Conv.Tag_size_adu.csv",row.names=F)
  if(exists('Zn.rel_Conv.Tag_size_juv')) write.csv(Zn.rel_Conv.Tag_size_juv,"Zn.rel_Conv.Tag_size_juv.csv",row.names=F)
  if(exists('Zn.rec_Conv.Tag_size_adu')) write.csv(Zn.rec_Conv.Tag_size_adu,"Zn.rec_Conv.Tag_size_adu.csv",row.names=F)
  if(exists('Zn.rec_Conv.Tag_size_juv')) write.csv(Zn.rec_Conv.Tag_size_juv,"Zn.rec_Conv.Tag_size_juv.csv",row.names=F)
  if(exists('Smallest_size_tagged')) write.csv(Smallest_size_tagged,"Smallest_size_tagged.csv",row.names=F)

  #Acoustic tagging 
    #releases
  if(exists('Zn.rel_Acous.Tag')) write.csv(Zn.rel_Acous.Tag,"acous.tag.rel.csv",row.names=F)  
  if(exists('Zn.rel_Acous.Tag.prop')) write.csv(Zn.rel_Acous.Tag.prop,"acous.tag.rel.prop.csv",row.names=F)
    #recaptures
  if(exists('Zn.rec_Acous.Tag')) write.csv(Zn.rec_Acous.Tag,"acous.tag.reca.csv",row.names=F)
  if(exists('Zn.rec_Acous.Tag.prop')) write.csv(Zn.rec_Acous.Tag.prop,"acous.tag.reca.prop.csv",row.names=F)
    #Individual-based model
  if(exists('Indiv_based_Acous.Tag')) write.csv(Indiv_based_Acous.Tag,"acous.tag.reca.indv.based.csv",row.names=F)
  #Recapture information from acoustic tags
  if(exists('Rep.Recap')) write.csv(Rep.Recap,"Acoustic.tag.rel.rec.csv",row.names=F) 
  
  #Age and growth 
  if(exists('Age.growth'))  write.csv(Age.growth,"Age.growth.csv",row.names=F)
  
  #Mesh proportional effort
  if(exists('Mesh.prop.eff')) write.csv(Mesh.prop.eff,"Mesh.prop.eff.csv",row.names=F)
  if(exists('Mesh.prop.eff.West')) write.csv(Mesh.prop.eff.West,"Mesh.prop.eff.West.csv",row.names=F)
  if(exists('Mesh.prop.eff.Zn1')) write.csv(Mesh.prop.eff.Zn1,"Mesh.prop.eff.Zn1.csv",row.names=F)
  if(exists('Mesh.prop.eff.Zn2')) write.csv(Mesh.prop.eff.Zn2,"Mesh.prop.eff.Zn2.csv",row.names=F)
  
  #Gillnet selectivity
  if(exists('Gillnet.selectivity'))  write.csv(Gillnet.selectivity,"Gillnet.selectivity.csv",row.names=F)
  if(exists('Gillnet.selectivity_len.age'))  write.csv(Gillnet.selectivity_len.age,"Gillnet.selectivity_len.age.csv",row.names=F)

  
  #remove stuff
  Rem.stuff=c('FL.TDGDFL','FL_Pilbara_trawl','FL_NSF','FL_Survey','All.size','All.size_7',
    'Ab.index.Srvy.FixSt','Zn.rel_Conv.Tag','Zn.rel_Acous.Tag','Mesh.prop.eff',
    'TDGDFL.size.numbers','Avr.wt.yr','Avr.wt.yr.zn','F.from.tagging','Eff.zn',
    'Eff.total','Ab.indx.TDGDLF.all','Ab.indx.TDGDLF.all.daily','Ab.indx.NSF',
    'Naturaliste.abun','Ab.indx.observer.TDGDLF','Ab.indx.PFT','Abun',
    'Size.index.Srvy.FixSt','Survey mean size','c.Tag','a.Tag',
    'Zn.rel_Conv.Tag_size_adu','Rep.Recap')
  pp=sapply(Rem.stuff,clear.log)
}

fn.add.future.sel=function(d)
{
  Add.future=d[rep(nrow(d),List.sp[[l]]$Yrs.future),]
  lst.yr=d$finyear[nrow(d)]
  y1=as.numeric(substr(lst.yr,1,4))+1
  y2=as.numeric(substr(lst.yr,6,7))+1
  Add.future$finyear=paste(y1:(y1+List.sp[[l]]$Yrs.future-1),y2:(y2+List.sp[[l]]$Yrs.future-1),sep="-")
  return(rbind(d,Add.future))
}
fn.varying.sel=function(Sel_6.5,Sel_7,TIME)
{
  yrs=TIME$finyear
  Str=matrix(rep(rep(NA,length(Sel_6.5)),length(yrs)),ncol=length(yrs),byrow=T)
  Str=as.data.frame(Str)
  names(Str)=yrs
  for(yyy in 1:length(yrs))
  {
    prop.eff=subset(TIME,finyear==yrs[yyy],select=c(Mesh_6.5,Mesh_7))
    Sel.dummy=data.frame(Mesh_6.5=Sel_6.5,Mesh_7=Sel_7)
    Sel.dummy$Mesh_6.5=Sel.dummy$Mesh_6.5*prop.eff$Mesh_6.5
    Sel.dummy$Mesh_7=Sel.dummy$Mesh_7*prop.eff$Mesh_7
    Sel.sum=rowSums(Sel.dummy)
    Str[,yyy]=round(Sel.sum/max(Sel.sum),4)
  }
  return(Str)
}

#------------- Execute functions  ----------- 
for(l in 1:N.sp)
{
  Neim=names(List.sp)[l]
  print(paste("---------Export data inputs for -- ",Neim))
  fn.input.data(Name=Neim,
                Name.inputs=List.sp[[l]]$Name.inputs,
                SP=List.sp[[l]]$SP,
                Species=List.sp[[l]]$Species,  
                First.year=List.sp[[l]]$First.year,
                Last.year=Last.yr.ktch,
                Min.obs=Min.obs,
                Min.shts=Min.shts,
                What.Efrt=What.Effort,
                Bin.size=TL.bins.cm,
                Yr.assess=AssessYr,
                Dat=Species.data[[l]],
                LH.par=LH.data%>%filter(SPECIES==List.sp[[l]]$Species),
                HandL=HandL.out) 
  
  
  #Output total catch and abundance indices to see contrast 
  fn.ktch.cpue(ktch=ktch.combined%>%
                      filter(Name==names(Catch.rate.series)[l]),
               cpue=Catch.rate.series[[l]],
               HandL=HandL.out)
  
  
  #Create .dat & .ctl for bespoke integrated assessments
  #Fix this but should use  Auxiliary_functons.R Superseded, not doing bespoke model, only SS3
  do.this=FALSE
  if(do.this)if(List.sp[[l]]$Species%in%Indicator.species)
  {
    all.objects=objects()
    List.objs=unique(names(List.sp[[l]]))
    suppressWarnings(rm(list=all.objects[which(all.objects%in%List.objs)]))
    with(List.sp[[l]],{
      #Current catch
      ktch=ktch.combined%>%
        filter(Name==Neim)%>%
        rename(Year=finyear,
               Total=Tonnes)%>%
        ungroup()%>%
        dplyr::select(Year,Total)%>%
        arrange(Year)%>%
        data.frame
      
      yr.start=min(ktch$Year)
      yr.end=max(ktch$Year)
      
      #Future catch
      Future.ktch=rep(mean(ktch$Total[(nrow(ktch)-List.sp[[l]]$Yrs.future+1):nrow(ktch)]),List.sp[[l]]$Yrs.future)
      
      #Population Length bins (cm)
      if(MN.SZE=="size.at.birth") Min.size.bin=10*round((Lzero*a_FL.to.TL+b_FL.to.TL)/10)
      if(MN.SZE==0) Min.size.bin=0
      MaxLen= 10*round(TLmax*Plus.gp.size/10)
      lbnd = seq(Min.size.bin,MaxLen - TL.bins.cm, TL.bins.cm)
      ubnd = lbnd + TL.bins.cm
      midpt = lbnd + (TL.bins.cm/2)
      
      #Maturity at total length
      Mat.fem=round(1/(1+(exp(-log(19)*((midpt-TL.50.mat)/(TL.95.mat-TL.50.mat))))),3)
      
      #total weight (kg) at length
      TwT.M=AwT.M*midpt^BwT.M
      TwT.F=AwT*midpt^BwT
      
      #Length at age and growth pars
      TL.zero=Lzero*a_FL.to.TL+b_FL.to.TL
      Growth.fem=Growth.F%>%mutate(TL_inf=FL_inf*a_FL.to.TL+b_FL.to.TL)
      Growth.mal=Growth.M%>%mutate(TL_inf=FL_inf*a_FL.to.TL+b_FL.to.TL)
      TL.F=TL.zero+(Growth.fem$TL_inf-TL.zero)*(1-exp(-Growth.fem$k*(1:Max.age.F[2])))
      
      #Natural mortality at length
      Nat.mort=rbind(data.frame(TL=TL.F,
                                M=apply(store.species.M_M.min[[l]],2,mean,na.rm=T)),
                     data.frame(TL=midpt,M=NA))%>%  #interpolate for size classes
        arrange(TL)%>%
        mutate(M=na.approx(M,maxgap=4))%>%
        fill(M, .direction = 'up')%>%
        fill(M, .direction = 'down')%>%
        filter(TL%in%midpt)%>%
        pull(M)
      
      #TDGDLF proportional effort for each mesh
      hndl.Rep2=paste(Dat.repository2,capitalize(str_remove(Neim, " shark")),"/",AssessYr,sep='')
      Mesh.prop.eff=read.csv(paste(hndl.Rep2,"Mesh.prop.eff.csv",sep='/'))
      iid=match(Last.yr.ktch,Mesh.prop.eff$finyear)
      Mesh.prop.eff=Mesh.prop.eff[1:iid,]
      Mesh.prop.eff.West=read.csv(paste(hndl.Rep2,"Mesh.prop.eff.West.csv",sep='/'))[1:iid,]
      Mesh.prop.eff.Zn1=read.csv(paste(hndl.Rep2,"Mesh.prop.eff.Zn1.csv",sep='/'))[1:iid,]
      Mesh.prop.eff.Zn2=read.csv(paste(hndl.Rep2,"Mesh.prop.eff.Zn2.csv",sep='/'))[1:iid,]
      
      Mesh.prop.eff=fn.add.future.sel(Mesh.prop.eff)  #add future Mesh props for projections
      Mesh.prop.eff.West=fn.add.future.sel(Mesh.prop.eff.West)
      Mesh.prop.eff.Zn1=fn.add.future.sel(Mesh.prop.eff.Zn1)
      Mesh.prop.eff.Zn2=fn.add.future.sel(Mesh.prop.eff.Zn2)
      
      # TDGDLF selectivity by year-region-mesh   
      Sel_6.5=round(Selectivity.at.totalength[[l]]%>%filter(TL%in%midpt)%>%pull(X16.5),4)
      Sel_7=round(Selectivity.at.totalength[[l]]%>%filter(TL%in%midpt)%>%pull(X17.8),4)
      Total.sel.all=fn.varying.sel(Sel_6.5, Sel_7, Mesh.prop.eff)
      Total.sel.West=fn.varying.sel(Sel_6.5, Sel_7, Mesh.prop.eff.West)    
      Total.sel.Zn1=fn.varying.sel(Sel_6.5, Sel_7, Mesh.prop.eff.Zn1)
      Total.sel.Zn2=fn.varying.sel(Sel_6.5, Sel_7, Mesh.prop.eff.Zn2)
      
      #Size composition (original data is FL)
      #TDGDLF
      Size.comp=Species.data[[l]][grep(paste(c('Size_composition_West','Size_composition_Zone1',
                                               'Size_composition_Zone2'),collapse="|"),
                                       names(Species.data[[l]]))]
      Size.comp=Size.comp[!grepl('Table',names(Size.comp))]
      for(p in 1:length(Size.comp))
      {
        Size.comp[[p]]=Size.comp[[p]]%>%
          mutate(dummy=gsub("\\.inch.raw*","",sub(" *Size_composition_", "", names(Size.comp)[p])),
                 Zone=gsub("\\..*","",dummy),
                 Mesh=sub(paste(unique(Zone),'.',sep=''),'',dummy))%>%
          dplyr::select(-dummy)
      }
      Size.comp=do.call(rbind,Size.comp)%>%
        mutate(Fishery='TDGDLF')
      
      #NSF
      if(length(grep('Size_composition_NSF',names(Species.data[[l]]))>0))
      {
        NSF.size.comp=Species.data[[l]]$Size_composition_NSF.LONGLINE%>%
          mutate(Fishery='NSF',
                 Zone=NA,
                 Mesh=NA)
        Size.comp=rbind(Size.comp,NSF.size.comp)
      }
      Size.comp=Size.comp%>%
        filter(FL<=quantile(Size.comp$FL,probs=.9995))%>%
        mutate(Year=as.numeric(substr(FINYEAR,1,4)),
               SEX=ifelse(SEX=='M',1,2),
               TL=FL*a_FL.to.TL+b_FL.to.TL,
               bin.mid=TL.bins.cm*floor(TL/TL.bins.cm)+TL.bins.cm/2,
               bin.lower=TL.bins.cm*floor(TL/TL.bins.cm))%>%
        filter(!is.na(SEX))
      
      #Abundance indices 
      CPUE=compact(Catch.rate.series[[l]])
      DROP=grep(paste(c('observer','West','Zone'),collapse="|"),names(CPUE))   
      if(length(DROP)>0)CPUE=CPUE[-DROP]
      if(Neim%in%survey_not.representative) CPUE=CPUE[-grep("Survey",names(CPUE))]
      if(Neim%in%NSF_not.representative) CPUE=CPUE[-grep("NSF",names(CPUE))]
      if(Neim%in%tdgdlf_not.representative) CPUE=CPUE[-grep("TDGDLF",names(CPUE))]
      for(x in 1:length(CPUE)) 
      {
        if(drop.large.CVs)
        {
          iid=which(CPUE[[x]]$CV>List.sp[[l]]$MAX.CV)
          CPUE[[x]]$Mean[iid]=NA
          CPUE[[x]]$CV[iid]=NA 
        }
      }
      if("TDGDLF.daily"%in%names(CPUE))
      {
        CPUE$TDGDLF.daily=rbind(data.frame(Finyear='2006-07', Mean=-100, CV=-100, LOW.CI=-100, UP.CI=-100, yr.f=2006),
                                CPUE$TDGDLF.daily)
      }
      
      #Conventional tagging
      if(Move.mode=="Individual-based")
      {
        Conv.tag=Species.data[[l]]$Con_tag_Ind.based.mod%>%
          filter(DaysAtLarge>=MIN.DAYS.LARGE)%>%
          dplyr::select(DaysAtLarge,Rel.zone,Rec.zone)%>%
          mutate(Rel.zone=case_when(Rel.zone=='West'~1,
                                    Rel.zone=='Zone1'~2,
                                    Rel.zone=='Zone2'~3),
                 Rec.zone=case_when(Rec.zone=='West'~1,
                                    Rec.zone=='Zone1'~2,
                                    Rec.zone=='Zone2'~3))
      }
      if(Move.mode=="Population-based")
      {
        Conv.tag.rel=Species.data[[l]]$Con_tag_Zn.rel_Conv.Tag_size
        Conv.tag.rec=Species.data[[l]]$Con_tag_Zn.rec_Conv.Tag_size
      }
      
      #acoustic tagging
      if(Acoust.format=="Individual-based")
      {
        Acous.tg.rec=Species.data[[l]]$Acous.Tag_Acous.Tag.Ind_based
        Acous.nzones=3
        Acous.detections=nrow(Acous.tg.rec)
        Acous.ntags=length(unique(Acous.tg.rec$TagID))
        Acous.TagID_observations=c(table(Acous.tg.rec$TagID))
        names(Acous.TagID_observations)=NULL
        Acous.tg.rec$Index=1:nrow(Acous.tg.rec)
        Acous.TagID_start=do.call("rbind",by(Acous.tg.rec,Acous.tg.rec$TagID,head,1))$Index
        Acous.TagID_end=do.call("rbind",by(Acous.tg.rec,Acous.tg.rec$TagID,tail,1))$Index
        Acous.tg.rec=Acous.tg.rec%>%
          mutate(TagID=as.numeric(substr(TagID,4,6)))%>%
          dplyr::select(DaysAtLarge,Rel.zn,Rec.zn,TagID)%>%
          mutate(Rel.zn=case_when(Rel.zn=='West'~1,
                                  Rel.zn=='Zone1'~2,
                                  Rel.zn=='Zone2'~3),
                 Rec.zn=case_when(Rec.zn=='West'~1,
                                  Rec.zn=='Zone1'~2,
                                  Rec.zn=='Zone2'~3))
      }
      
      #Recruitment error for future projections
      Rec.error=rlnorm(nrow(ktch)+length(Future.ktch), meanlog = log(1), sdlog = Sens.test$CMSY$Proc.error[2])  
      
      #Age growth
      Age.growth=Species.data[[l]]$age_length%>%
        mutate(TL=FL*a_FL.to.TL+b_FL.to.TL)
      Age.growth.F=Age.growth%>%filter(Sex=='Female')%>%dplyr::select(Age,TL)
      Age.growth.M=Age.growth%>%filter(Sex=='Male')%>%dplyr::select(Age,TL)
      
      if(do.Andre.model)
      {
        #Fleets
        Andre.fleets=data.frame(Fleet.name=c("Survey","NSF","TDGDLF.monthly","TDGDLF.daily"),
                                Fleet=0:3)
        
        #Size composition (original data is FL)
        Andre.SBIM.size=Size.comp%>%
          group_by(Year,bin.lower,SEX,Fishery)%>%
          tally()%>%
          ungroup()
        Misn.siz=lbnd[which(!lbnd%in%Andre.SBIM.size$bin.lower)] #add missing size classes
        if(length(Misn.siz)>0)
        {
          Misn.siz=data.frame(Year=NA,
                              bin.lower=Misn.siz,
                              SEX=NA,
                              Fishery=NA,
                              n=0)
          Andre.SBIM.size=rbind(Andre.SBIM.size,Misn.siz)
        }
        Andre.SBIM.size=Andre.SBIM.size%>%
          spread(bin.lower,n,fill=0)%>%
          arrange(Fishery,SEX,Year)%>%
          rename(SEASON=Year,
                 Sex=SEX)%>%
          mutate(Ind=rowSums(select(., -c(Fishery,SEASON,Sex))),
                 Fleet=case_when(SEASON<=2005 & Fishery=='TDGDLF' ~ Andre.fleets%>%filter(Fleet.name=='TDGDLF.monthly')%>%pull(Fleet),
                                 SEASON>2005 & Fishery=='TDGDLF' ~ Andre.fleets%>%filter(Fleet.name=='TDGDLF.daily')%>%pull(Fleet),
                                 Fishery=='NSF' ~ Andre.fleets%>%filter(Fleet.name=='NSF')%>%pull(Fleet)),
                 tstep=1)%>%
          dplyr::select(-"Fishery")%>%
          relocate(Fleet,Sex,SEASON,tstep,Ind)
        Andre.SBIM.size$''=0  #add extra column of 0s for Andre Model
        
        #Weight (kg) at total length
        Andre.SBIM.kg.at.TL=t(data.frame(male=round(TwT.M,3),
                                         female=round(TwT.F,3)))
        Andre.SBIM.mat.at.TL=t(data.frame(female=Mat.fem))  
        
        #Fecundity at total length  (maturity X # of pups)      MISSING: use fecundity relationship?; Cycle length not considered???
        Andre.SBIM.fecun.at.TL=t(data.frame(female=mean(Fecundity)*Mat.fem))
        
        #Catch (kg) 
        Andre.SBIM.ktch=ktch.combined%>%
          filter(Name==Neim)%>%
          ungroup()%>%
          mutate(Year=finyear,
                 step=1,
                 fleet=1,	   #all fleets combined
                 catch=round(Tonnes*1000))%>%
          dplyr::select(Year,step,fleet,catch)
        
        #Discard mortality  
        Andre.SBIM.disc.mort=cbind(data.frame(Age=0,
                                              Fleet=Andre.fleets$Fleet,
                                              tstep=0),
                                   as.data.frame(matrix(0,
                                                        nrow=length(Andre.fleets$Fleet),
                                                        ncol=length(Andre.SBIM.ktch$Year))))
        colnames(Andre.SBIM.disc.mort)=c('Age','Fleet','tstep',Andre.SBIM.ktch$Year)
        
        #Fleet-Area
        Fleet.Area=data.frame(Fleet=Andre.fleets$Fleet,
                              Area=0)
        
        #Recruitment_deviations
        startseason=min(Andre.SBIM.ktch$Year)
        endseason=max(Andre.SBIM.ktch$Year)
        Rec.dev.phase=1
        Rec.dev.phase_spatial=-1
        
        #Prespecify_rec_devs
        Prespecify_rec_devs=data.frame(col1=0,
                                       Numeral.in.between='#',
                                       Years=c(Andre.SBIM.ktch$Year,seq(endseason+1,endseason+List.sp[[l]]$Yrs.future)))
        #Prespecify_spatial_rec_devs
        Prespecify_spatial_rec_devs=data.frame(col1=0,
                                               Numeral.in.between='#',
                                               Years=c(Andre.SBIM.ktch$Year,seq(endseason+1,endseason+List.sp[[l]]$Yrs.future)))
        Prespecify_spatial_rec_devs=rbind(Prespecify_spatial_rec_devs%>%mutate(Col2=2),
                                          Prespecify_spatial_rec_devs%>%mutate(Col2=1))
        
        #Weights on the data
        Weight_cpue=1
        Weight_ktch=1
        Weight_lenfreq=0.001
        Weight_larval=1
        Weight_tag1=0
        Weight_tag2=1
        
        #Weights by fleet
        
        #Basic parameters
        Basic.pars=data.frame(lower=c(.01, 0, 0, 1, 0, .001),
                              upper=c(2, 1, 10, 5, 1, 1),
                              estimate=c(.05, .12, 1, 3, 1, .5),
                              phase=c(2, -1, -1, -1, -1, --1),
                              Pars=paste('#',c('Rbar', "M area 1","Multiplier for age 0",
                                               "WhiteScaleM", "RedScaleQ","SigmaR")))
        
        #Initial_dev_option
        Initial_dev_option=0 
        Initial_dev_option.val=1 
        
        #Number of movement patterns
        N.move.zones=2
        Move.pattern.mat=data.frame(Pattern=0:1,
                                    Type=0:1,
                                    Dest=rep(0,2),
                                    Extra=rep(0,2))
        Movement.specifications=cbind(data.frame(Age=0,Area=0,TStep=0),
                                      as.data.frame(matrix(0,ncol=length(startseason:endseason))))
        Movement.parameters=data.frame(lower=0.1,
                                       upper=1,
                                       estimate=0.5,
                                       phase=-1)
        
        #Number of growth patterns
        N.growth.patrn=2
        Growth.pattern.mat=data.frame(Pattern=0:1,
                                      Type=rep(1,2),
                                      Sex=c(0,1),
                                      Extra=rep(0,2),
                                      Pointer=c(0,1),
                                      Mpower=rep(1,2),
                                      Months=c('# 12','# 12'),
                                      growthareas=rep(1,2))
        Spec.for.growth.mat=cbind(data.frame(Sex=c(0,1),
                                             Age=rep(0,2),
                                             Area=rep(0,2),
                                             Step=rep(0,2)),
                                  as.data.frame(matrix(0,nrow=2,ncol=length(startseason:endseason))))
        
        #Size transition matrix
        Num.size.trans=2
        Andre.SBIM.STM.F=Sadovy.size.trans.mat(Linf=Growth.fem$TL_inf,
                                               K=Growth.fem$k,
                                               SD=SD.growth_SB,  
                                               bin.low=lbnd,
                                               bin.up=ubnd,
                                               Truncate="YES")
        Andre.SBIM.STM.M=Sadovy.size.trans.mat(Linf=Growth.mal$TL_inf,
                                               K=Growth.mal$k,
                                               SD=SD.growth_SB,
                                               bin.low=lbnd,
                                               bin.up=ubnd,
                                               Truncate="YES")
        for(qq in 1:ncol(Andre.SBIM.STM.F))
        {
          Andre.SBIM.STM.F[,qq]=round(Andre.SBIM.STM.F[,qq],7)
          Andre.SBIM.STM.F[,qq]=Andre.SBIM.STM.F[,qq]/sum(Andre.SBIM.STM.F[,qq])
          
          Andre.SBIM.STM.M[,qq]=round(Andre.SBIM.STM.M[,qq],7)
          Andre.SBIM.STM.M[,qq]=Andre.SBIM.STM.M[,qq]/sum(Andre.SBIM.STM.M[,qq])
        }
        
        #Selectivity patterns
        Selectivity.patterns=2  #TDGDLF gillnet & NSF longline
        Andre.SBIM.sel.GN=Selectivity.at.totalength[[l]]%>%               
          mutate(TL=TL)%>%
          filter( TL%in%midpt)%>%
          pull(Sel.combined)
        Andre.SBIM.sel.LL=Mat.fem   #set NSF longline selectivity to maturity given species distribution
        
        #Abundance indices
        Andre.SBIM.cpue=CPUE
        for(x in 1:length(Andre.SBIM.cpue))    
        {
          dd=Andre.SBIM.cpue[[x]][,grep(paste(c('yr.f','Mean','MeAn','CV'),collapse="|"),names(Andre.SBIM.cpue[[x]]))]%>%
            relocate(yr.f)%>%
            mutate(CpueInd=x-1,
                   Fleet=Andre.fleets%>%filter(Fleet.name==names(Andre.SBIM.cpue)[x])%>%pull(Fleet),
                   Sex=0,
                   Step=1,
                   Year=yr.f,
                   CV=round(CV,2))
          names(dd)[2]='Index'
          dd=dd%>%
            filter(!is.na(Index))%>%
            mutate(Index=round(Index,4))%>%
            dplyr::select(CpueInd,Fleet,Sex,Year,Step,Index,CV)
          
          #remove daily 2007  
          if(drop2007.daily & 'TDGDLF.daily'%in%names(Andre.SBIM.cpue)[x]) dd=dd%>%filter(!Year%in%2007)
          
          Andre.SBIM.cpue[[x]]=dd
        }
        Andre.SBIM.cpue=do.call(rbind,Andre.SBIM.cpue)
        Type.of.index=ifelse(names(Andre.SBIM.cpue)=='Survey',2,1)
        
        
        #Create .dat file
        create.ass.inputs(WD=handl_OneDrive('Workshops/2022_Andre Punt_LengthBased/Sandbar'),
                          Outfile=capitalize( gsub(" ", "", List.sp[[l]]$Name, fixed = TRUE)),
                          Input_list=list(
                            '# First year of the assessment\n',  min(KtCh$finyear),
                            '\n# Last year of the assessment\n',  Last.yr.ktch.numeric,
                            '\n# Maximum projection years\n',  List.sp[[l]]$Yrs.future,
                            '\n# Burn-in\n',15,
                            '\n# Time steps per year\n', 1,
                            '\n# Number of areas of data included in the file\n', 1,
                            '\n# Number of sexes\n', 2,
                            '\n# Number of ages\n', 1,
                            '\n# Number of fleets\n', length(Type.of.index),
                            '\n# Number of size-classes (males then females)\n', rep(length(lbnd)-1,2),
                            '\n# The Time steps\n', 1,
                            '\n# Loop counter for initial conditions\n', 2,
                            '\n# Years over which to tune\n', 15,
                            c('\n# Lower Length Bins (one more than number of size-classes)\n',
                              paste(lbnd, collapse = "\t"),
                              "\n",
                              paste(lbnd, collapse = "\t")),
                            '\n# Catch data (kg) - Number of observations\n',  nrow(Andre.SBIM.ktch),
                            '\n#Year	step	fleet	catch\n', Andre.SBIM.ktch,
                            '\n# Index data',
                            '\n#Number of cpue datasets\n', length(unique(Andre.SBIM.cpue$Fleet)),
                            '\n# Type of index (1=weight;2=numbers)\n', paste(Type.of.index,collapse='\t'),
                            '\n# Treatment of sigma (unique value represents a unique SS for the series)\n', paste(seq(0,length(Type.of.index)-1),collapse='\t'),
                            '\n# Treatment of q (a value represents a unique q for that series)\n', paste(seq(0,length(Type.of.index)-1),collapse='\t'),
                            '\n# Environmental Index / Fishing Efficiency (value points to index, 0 = no index)\n', paste(rep(0,length(Type.of.index)-1),collapse='\t'),
                            '\n# Minimum sigma\n', 0.05,
                            '\n# The cpue data\n',  nrow(Andre.SBIM.cpue),
                            '\n#CpueInd Fleet	Sex	Year	Step	Index	CV\n', Andre.SBIM.cpue,
                            '\n# Numbers data (recreational fishery only)',
                            '\n# Number of catch data sets\n', 0,
                            '\n# Treatment of catch in numbers series\n', '#0',
                            '\n# Minimum sigma\n', 0.05,
                            '\n# The numbers data\n', 0,
                            '\n#Group  Fleet  Year  Step  Catch  CV',
                            '\n',
                            '\n# Length composition\n', nrow(Andre.SBIM.size),
                            fn.col.nms.c(colnames(Andre.SBIM.size)),
                            Andre.SBIM.size%>%
                              data.frame%>%
                              filter(!is.na(Fleet))%>%
                              dplyr::select(-Fleet,Sex,SEASON,tstep,Ind),
                            '\n# Larval index (puerulus)',
                            '\n# Likelihood for larval (puerulus) data (0=lognormal, else normal\n', 0,
                            '\n# Delay from puerulus to entering the model (years)\n', 0,
                            '\n# Number of data points (puerulus samples)\n', 0,
                            '\n# Area	Year	Index	SD\n', 0,
                            '\n# Environmental Data - ln(Efficiency creep)',
                            '\n# Number of environmental series (commercial efficiency creep for each area)\n', 0,
                            '\n# Years of data per series\n', 0,
                            '\n# Year	Tstep	ln(creep)	Series\n', 0,
                            '\n',
                            '\n# Final check\n', 123456))
        
        #Create .CTL file
        create.ass.inputs(WD=handl_OneDrive('Workshops/2022_Andre Punt_LengthBased/Sandbar'),
                          Outfile=capitalize( gsub(" ", "", List.sp[[l]]$Name, fixed = TRUE)),
                          Input_list=list(
                            '\n# weight-at-length (W=aL^b) (kg)\n', Andre.SBIM.kg.at.TL,
                            '\n# Maturity (by area)\n', Andre.SBIM.mat.at.TL,
                            '\n# Egg production (by area)\n', Andre.SBIM.fecun.at.TL,
                            '\n# Egg time step (Assume this is when to determine egg production?? time step 6 [6-1])\n', 1,
                            '\n# Fleet specification',
                            fn.col.nms.c(colnames(Fleet.Area)), Fleet.Area,
                            '\n# Number of Zones\n', 1,
                            '\n# Areas in each Zone\n', 1,
                            '\n# The Zones\n', 0,
                            '\n# Discard mortality',
                            fn.col.nms.c(colnames(Andre.SBIM.disc.mort)),Andre.SBIM.disc.mort%>%
                              data.frame%>%
                              dplyr::select(-Age,Fleet,tstep),
                            fn.paste.c(startseason, "first year with estimated recruitment deviations\n"),
                            fn.paste.c(endseason, "last year with estimated recruitment deviations\n"),
                            fn.paste.c(Rec.dev.phase,"phase for recruitment deviations\n"),
                            '# Spatial_deviations_in_recruitment\n',
                            fn.paste.c(startseason, "first year with estimates spatial recruitment deviations\n"),
                            fn.paste.c(endseason, "last year with estimates spatial recruitment deviations\n"),
                            fn.paste.c(Rec.dev.phase_spatial,"phase for spatial recruitment deviations\n\n"),
                            c("# Prespecify_rec_devs :  dev # Year\n",1,"\t\t\t# 1 = rec_devs are to be pre-specified\n"),
                            Prespecify_rec_devs,
                            c("\n# Prespecify_spatial_rec_devs : dev # Year Area\n",0,"\t\t\t# 1 = Spatial_rec_devs are to be pre-specified\n"),
                            Prespecify_spatial_rec_devs,
                            '\n# Weights on the data (simple)\n',
                            fn.paste.c(Weight_cpue, "Weight on CPUE data\n"),
                            fn.paste.c(Weight_ktch, "Weight on catch-numbers data\n"),
                            fn.paste.c(Weight_lenfreq, "Weight on Length-frequency data\n"),
                            fn.paste.c(Weight_larval, "Weight on larval data\n"),
                            fn.paste.c(Weight_tag1, "Weight on Tag1 data\n"),
                            fn.paste.c(Weight_tag2, "Weight on Tag2 data\n"),
                            '# Weights by fleet (Type: 1=Cpue,2=Numbers,3=Length;4=Larva;) - Can tweak by Fleet and index.  Set values to -1 for them to encompass all options.  e.g. time-step set to -1 covers all timesteps. our fleets',
                            '\n# Type	Fleet	Time-step	sex\n', fn.paste.c(0,"set to >0 for yes"),
                            '\n# Basic parameters (lower, upper, estimate, phase)\n',
                            Basic.pars,
                            '\n# Initial_dev_option\n',
                            fn.paste.c(Initial_dev_option, "0=convetional; 1=alternative"),
                            fn.paste.c(Initial_dev_option.val, "Initial value options (0=default; 1=same for all)"),
                            '\n# Initial size parameters    - 0=convetional; 1=alternative\n', paste(c(-100,100,0,1),collapse = "\t"),
                            '\n# Q parameters    - Multiplier for Efficiency creep - Keep to 1\n', paste(c(-100,100,1,1),collapse = "\t"),
                            '\n# variance specification parameters (1=SSB;2=Recruitment;3legal biomass by area;4:exploitation rate by zone)',
                            '\n# Number of variance specifications\n', 4,
                            '\n# variance components\n', paste(1:4,collapse = "\t"),
                            '\n',
                            fn.paste.c(1,"use the pin file for specifying parameters (ADMB)\n"),
                            fn.paste.c(1,"last function call (ADMB)\n"),
                            '\n# Final check\n',123456),
                          extension=".CTL")
        
        #Create MOVE file
        create.ass.inputs(WD=handl_OneDrive('Workshops/2022_Andre Punt_LengthBased/Sandbar'),
                          Outfile=paste(capitalize( gsub(" ", "", List.sp[[l]]$Name, fixed = TRUE)),'MOVE',sep=''),
                          Input_list=list(
                            '# Number of movement patterns\n',N.move.zones,
                            c('\n#',paste(colnames(Move.pattern.mat),collapse = '\t'),'(Type 0: none; 1 constant [prespecified or estimated; 1 parameter]; 2 knife-eded-specific [pre-specified or estimated; 2 parameters])\n'),
                            Move.pattern.mat,
                            '\n# Movement specifications',
                            c(paste("\n#Age\tArea\tTStep\t", paste(startseason:endseason,collapse = "\t")),'\n'),
                            Movement.specifications,
                            '# Movement parameters',
                            c('\n#',paste(colnames(Movement.parameters),collapse='\t'),'\tSource to Dest\n'),
                            Movement.parameters,
                            '\n# Final check\n',123456))
        
        #Create GROWTH file
        create.ass.inputs(WD=handl_OneDrive('Workshops/2022_Andre Punt_LengthBased/Sandbar'),
                          Outfile=paste(capitalize( gsub(" ", "", List.sp[[l]]$Name, fixed = TRUE)),'GROWTH',sep=''),
                          Input_list=list(
                            '# Growth specification (Months represent the number of times a monthly-STM has been compounded)',
                            '\n# Number of growth Patterns (Mpower is legacy and needed for power function on growth)\n',N.growth.patrn,
                            fn.col.nms.c(colnames(Growth.pattern.mat)),Growth.pattern.mat,
                            '\n# Specifications for growth',
                            fn.col.nms.c(colnames(Spec.for.growth.mat)),Spec.for.growth.mat,
                            '\n# Number prespecified size-transition\n',Num.size.trans,
                            '\n# Sex of matrices\n',matrix(c(0,1),nrow=1),
                            '\n# Prespecified size-transition',
                            '\n# Matrix 0 Mstm1 months_12\n', Andre.SBIM.STM.M,
                            '\n# Matrix 1 Fstm1 months_12\n', Andre.SBIM.STM.F,
                            '\n# Growth parameters',
                            '\n#LB	UP	Estimate	Phase',
                            '\n# Final check\n',123456))
        
      }
      
      if(Do.bespoked)
      {
        SIZE.comp.nobs=rbind(Size.comp%>%
                               dplyr::select(Year,Mesh,SEX,Fishery),  #1,male; 2, female
                             expand.grid(Year=ktch$Year,
                                         Mesh=unique(Size.comp$Mesh),
                                         SEX=unique(Size.comp$SEX),
                                         Fishery=unique(Size.comp$Fishery)))%>%
          group_by(Year,Mesh,SEX,Fishery)%>%
          tally()%>%
          ungroup()%>%
          mutate(n=ifelse(n==1,0,n-1))%>%
          filter(Year%in%ktch$Year)
        First.yr.size.obs=match(min(SIZE.comp.nobs%>%filter(n>0)%>%pull(Year)),ktch$Year)
        Last.yr.size.obs=match(max(SIZE.comp.nobs%>%filter(n>0)%>%pull(Year)),ktch$Year)
        
        SIZE.comp.prop=rbind(Size.comp%>%
                               dplyr::select(Year,bin.mid,Mesh,SEX,Fishery),
                             expand.grid(Year=ktch$Year,
                                         bin.mid=midpt,
                                         Mesh=unique(Size.comp$Mesh),
                                         SEX=unique(Size.comp$SEX),
                                         Fishery=unique(Size.comp$Fishery)))%>%
          group_by(Year,bin.mid,Mesh,SEX,Fishery)%>%
          tally()%>%
          ungroup()%>%
          mutate(n=ifelse(n==1,0,n-1))%>%
          filter(Year%in%ktch$Year)%>%
          group_by(Year,Mesh,SEX,Fishery)%>%
          mutate(Sum=sum(n))%>%
          ungroup()%>%
          mutate(Prop=ifelse(Sum>0,n/Sum,0))
        
        First.size.obs=match(min(Size.comp%>%filter(Mesh==6.5)%>%pull(bin.mid)),midpt)
        Last.size.obs=match(max(Size.comp%>%filter(Mesh==6.5)%>%pull(bin.mid)),midpt)
        First.size.obs_7=match(min(Size.comp%>%filter(Mesh==7)%>%pull(bin.mid)),midpt)
        Last.size.obs_7=match(max(Size.comp%>%filter(Mesh==7)%>%pull(bin.mid)),midpt)  
        
        if(Size_like=="Multinomial") Type_of_Size_like=0
        if(Size_like=="Dirichlet") Type_of_Size_like=1
        
        CPUE.int=do.call(rbind,CPUE)
        
        #populate scenarios
        for(s in 1:N.Scens)
        {
          nzones=str_extract(Tabla.scen$Spatial_structure[s], "[[:digit:]]+")
          syr_cpue=match(min(CPUE.int$yr.f),ktch$Year) 
          CPUE_eff=rep(-100,nrow(CPUE.int))
          if('CPUE_years_dropped'%in%colnames(Tabla.scen))
          {
            drop.this=as.numeric(unlist(strsplit(Tabla.scen$CPUE_years_dropped[s],'-')))
            CPUE.int=CPUE.int%>%
              mutate(Mean=ifelse(between(yr.f,drop.this[1],drop.this[2]),-100,Mean),
                     CV=ifelse(between(yr.f,drop.this[1],drop.this[2]),-100,CV))
            
          }
          CPUE.int_mean=CPUE.int$Mean
          CPUE.int_CV=CPUE.int$CV
          if(Tabla.scen$CPUE[s]=="Effective") CPUE_eff=CPUE.int$Mean  #dummy, not used  
          
          #create .dat file  
          create.ass.inputs(WD=paste(hndl,AssessYr,"/Integrated sized-structured/",
                                     Tabla.scen$Model[s],sep=""),
                            Outfile=capitalize(Name),
                            Input_list=list(
                              '# Basic model dimensions (yr.start, yr.end, nzones)\n',fn.paste.vec(c(yr.start,yr.end,nzones)),
                              '\n# yrs_ktch\n',nrow(ktch)+List.sp[[l]]$Yrs.future,
                              '\n# nTL\n',length(midpt),
                              '\n# size_mat\n',fn.paste.vec(Mat.fem),
                              '\n# sx_ratio\n',pup.sx.ratio,
                              '\n# mean.size.birth\n',round(TL.zero),
                              '\n# SD.mean.size.birth\n',Lzero_SD,
                              '\n# TL.bin.mid\n',fn.paste.vec(midpt),
                              '\n# TL.bin.low\n',fn.paste.vec(lbnd),
                              '\n# TL.bin.up\n',fn.paste.vec(ubnd),
                              '\n# Max.age.size\n',2*Max.age.F[1],       #max age for plus group
                              '\n# TwT.bin.F\n',fn.paste.vec(TwT.F),
                              '\n# TwT.bin.M\n',fn.paste.vec(TwT.M),
                              '\n# M\n',fn.paste.vec(Nat.mort),
                              fn.paste.vec(c('\n# Selectivity',substr(colnames(Total.sel.all),1,4),'\n')),
                              Total.sel.all,
                              fn.paste.vec(c('# Selectivity.west',substr(colnames(Total.sel.all),1,4),'\n')),
                              Total.sel.West,
                              fn.paste.vec(c('# Selectivity.zn1',substr(colnames(Total.sel.all),1,4),'\n')),
                              Total.sel.Zn1,
                              fn.paste.vec(c('# Selectivity.zn2',substr(colnames(Total.sel.all),1,4),'\n')),
                              Total.sel.Zn2,
                              '# SEL_6.5\n',fn.paste.vec(Sel_6.5),
                              '\n# SEL_7\n',fn.paste.vec(Sel_7),
                              '\n# STEEPns\n',Tabla.scen$SteepnesS[s],
                              '\n# iyr\n',fn.paste.vec(ktch$Year),
                              '\n# ct_F\n',fn.paste.vec(round(c(ktch$Total,Future.ktch)*(1-Prop.males.ktch),4)),
                              '\n# ct_M\n',fn.paste.vec(round(c(ktch$Total,Future.ktch)*(Prop.males.ktch),4)),    #assuming total catch has the same sex comp of TDGDLF
                              '\n# n.SIZE.comp.yrs.F\n',fn.paste.vec(SIZE.comp.nobs%>%
                                                                       filter(Mesh==6.5 & SEX==2 & Fishery=='TDGDLF')%>%
                                                                       pull(n)),
                              '\n# n.SIZE.comp.yrs.M\n',fn.paste.vec(SIZE.comp.nobs%>%
                                                                       filter(Mesh==6.5 & SEX==1 & Fishery=='TDGDLF')%>%
                                                                       pull(n)),
                              '\n# n.SIZE.comp.yrs.F_7\n',fn.paste.vec(SIZE.comp.nobs%>%
                                                                         filter(Mesh==7 & SEX==2 & Fishery=='TDGDLF')%>%
                                                                         pull(n)),
                              '\n# n.SIZE.comp.yrs.M_7\n',fn.paste.vec(SIZE.comp.nobs%>%
                                                                         filter(Mesh==7 & SEX==1 & Fishery=='TDGDLF')%>%
                                                                         pull(n)),
                              '\n# First.size.obs\n',First.size.obs,
                              '\n# Last.size.obs\n',Last.size.obs,
                              '\n# First.size.obs_7\n',First.size.obs_7,
                              '\n# Last.size.obs_7\n',Last.size.obs_7,
                              '\n# First.yr.size.obs\n',First.yr.size.obs,
                              '\n# Last.yr.size.obs\n',Last.yr.size.obs,
                              fn.paste.vec(c('\n# SIZE.comp.F',midpt,'\n')),
                              SIZE.comp.prop%>%
                                filter(Mesh==6.5  & SEX==2 & Fishery=='TDGDLF')%>%
                                dplyr::select(Year,bin.mid,Prop)%>%
                                spread(bin.mid,Prop)%>%
                                dplyr::select(-Year),
                              '# SIZE.comp.M\n',
                              SIZE.comp.prop%>%
                                filter(Mesh==6.5  & SEX==1 & Fishery=='TDGDLF')%>%
                                dplyr::select(Year,bin.mid,Prop)%>%
                                spread(bin.mid,Prop)%>%
                                dplyr::select(-Year),
                              '# SIZE.comp.F_7\n',
                              SIZE.comp.prop%>%
                                filter(Mesh==7  & SEX==2 & Fishery=='TDGDLF')%>%
                                dplyr::select(Year,bin.mid,Prop)%>%
                                spread(bin.mid,Prop)%>%
                                dplyr::select(-Year),
                              '# SIZE.comp.M_7\n',
                              SIZE.comp.prop%>%
                                filter(Mesh==7  & SEX==1 & Fishery=='TDGDLF')%>%
                                dplyr::select(Year,bin.mid,Prop)%>%
                                spread(bin.mid,Prop)%>%
                                dplyr::select(-Year),
                              '# rho\n',Rho,
                              '\n# Type_of_Size_like\n',Type_of_Size_like,
                              '\n# do_cpue\n',ifelse(Tabla.scen$CPUE[s]%in%c("N/A","No"),0,1),
                              '\n# effective_cpue\n',ifelse(Tabla.scen$CPUE[s]=="Effective",1,0),
                              '\n# syr_cpue\n',syr_cpue,
                              '\n# CPUE_eff\n',fn.paste.vec(CPUE_eff),
                              '\n# CPUE\n',fn.paste.vec(CPUE.int_mean),
                              '\n# CPUE.CV\n',fn.paste.vec(CPUE.int_CV),
                              '\n# Conv.tg.nzones\n',gsub(".*?([0-9]+).*", "\\1", Tabla.scen$Spatial_structure[s]),
                              '\n# Conv.tg.recaptures (DaysAtLarge Rel.zn Rec.zn)\n',nrow(Conv.tag),
                              '\n# Conv.tg.rec (DaysAtLarge Rel.zn Rec.zn)\n',Conv.tag,
                              '# Smallest_move\n',match(round((Species.data[[l]]$Con_tag_Smallest.size_Conv.Tag_size$x*a_FL.to.TL+b_FL.to.TL)/10)*10,round(midpt/10)*10),
                              '\n# Acous.tg.detections\n',Acous.detections,
                              '\n# Acous.tg.tags\n',Acous.ntags,
                              '\n# Acous.tg.obs.per.tags\n',fn.paste.vec(Acous.TagID_observations),
                              '\n# Acous.tg.start\n',fn.paste.vec(Acous.TagID_start),
                              '\n# Acous.tg.end\n',fn.paste.vec(Acous.TagID_end),
                              '\n# Acous.tg.data (DaysAtLarge Rel.zn Rec.zn TagID)\n',Acous.tg.rec,
                              '# Rec.error\n',fn.paste.vec(Rec.error),
                              '\n# q_change\n',ifelse(Yr_q_change>0,match(Yr_q_change,CPUE.int$yr.f),-100),
                              '\n# q_daily\n',match(Yr_q_daily,CPUE.int$yr.f),
                              '\n# Nobs.Age.Grw.F\n',nrow(Age.growth.F),
                              '\n# Nobs.Age.Grw.M\n',nrow(Age.growth.M),
                              '\n# Age.Grw.F (Counts TL)\n',Age.growth.F,
                              '# Age.Grw.M (Counts TL)\n',Age.growth.M,
                              '# Rho2\n',Rho2,
                              '\n# Rho3\n',Rho3,
                              '\n# MaxF\n',MaxF,
                              '\n# add_Finit_prior\n',add_Finit_prior,
                              '\n# mu_F_init\n',Prior.mean.Fo,
                              '\n# Ln_SD_F_init\n',Prior.SD.Log.Fo,
                              '\n# Phases\n',fn.paste.vec(Par.phases$`Base case`),
                              '\n# Do.move\n',ifelse(Tabla.scen$Movement[s]=="No",0,1),
                              '\n# eof\n',999))
        }
      }
      
    })
  }
} 


#------------- Remove functions  ----------- 
drop.funs=c('fn.exp.his', 'fn.add.missing.size', 'fn.add.missing.size2', 'fn.add.missing.size3', 'fn.add', 'fn.see.Ktch',
'visualize.dat', 'fn.agg.at.level.and.exprt', 'fn.exp.size', 'tagging.fn', 'fn.see.avg.wgt', 'fn.see.avg.wgt.zn',
'fn.recode', 'fn.plot.prop.rec.yr', 'drop.obs', 'fn.bck.fil', 'fn.bub', 'fn.table.shots','fn.input.data',
'fn.add.future.sel','fn.varying.sel')
for(d in 1:length(drop.funs))clear.log(drop.funs[d])