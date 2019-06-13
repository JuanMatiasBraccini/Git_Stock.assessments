library(plotrix)  #add table to plot

#Define species
#SPEC="Whiskery"
SPEC="Gummy"

#Define scenario
Scen="S1"
#Scen="S2"
Last.yr=2015

SCEN=paste("C:/Matias/Analyses/Population dynamics/",SPEC," shark/2017/",sep="")


#Bring in ADMB outputs
fn.get=function(x) round(LST[[match(c(x),names(LST))]],2)
reptoRlist = function(fn)
{
  ifile=scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
  idx=sapply(as.double(ifile),is.na)
  vnam=ifile[idx]  #list names
  nv=length(vnam)	#number of objects
  A=list()
  ir=0
  for(i in 1:nv)
  {
    ir=match(vnam[i],ifile)
    if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
    dum=NA
    if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=T,what=""))
    if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
    
    #Logical test to ensure dealing with numbers
    if(is.numeric(dum)) A[[vnam[i]]]=dum
  }
  return(A)
}
fn.which=function(D,n,x) which(substr(D[[1]]$std$name,1,n)==x)

if(Scen=="S1")
{
  if(SPEC=="Whiskery")Mods=c("S1","S1_pin_sensitivity","S1_no_fishing_efficiency","S1_po_0.7","S1_po_0.99")
  if(SPEC=="Gummy")Mods=c("S1","S1_pin_sensitivity","S1_po_0.7","S1_po_0.99")
  LIST.mod=vector('list',length(Mods))
  names(LIST.mod)=Mods
}
  
for(i in 1:length(LIST.mod))
{
  mod=names(LIST.mod)[i]
  setwd(paste(SCEN,mod,sep="/"))
  LIST.mod[[i]]=list(rep=reptoRlist(paste(SPEC,".rep",sep="")),
                     std=read.table(paste(SPEC,".std",sep=""), header =TRUE))
  
}

fn.plt=function(d,Compare.what,CEx.tab,x.tab,y.tab)
{
  setwd(paste(SCEN,"2_Outputs/Compare_scenarios",sep="/"))
  
  tiff(file=paste(SPEC,"_compare_plots_",Compare.what,".tiff",sep=''),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(3,1),mai=c(.35,.5,.05,.2),oma=c(.1,.2,1,1),mgp=c(1.5,.5,0),las=1)
  
  #Biomass
  Bio.1=d[[1]]$std$value[fn.which(D=d,n=10,x="Total_biom")]/1000
  Bio.1.SE=d[[1]]$std$std.dev[fn.which(D=d,n=10,x="Total_biom")]/1000
  Year=1975:(1975+length(Bio.1)-1)
   YMAX=max(Bio.1+1.96*Bio.1.SE)
   YMAX=rep(NA,length(d))
   for(i in 1:length(d))YMAX[i]=max(c(1.96*d[[i]]$std$std.dev+d[[i]]$std$value)/1000)
   YMAX=max(YMAX)
  plot(Year,Bio.1,ylab="",ylim=c(0,YMAX),xlim=c(1975,Last.yr),cex.lab=1.5,col="transparent")
  CL=c("black","red","blue","green","orange","gold")
  Nms=1:length(d)
  for(n in 1:length(d))
  {
    Bio=d[[n]]$std$value[fn.which(D=d,n=10,x="Total_biom")]/1000
    Bio.SE=d[[n]]$std$std.dev[fn.which(D=d,n=10,x="Total_biom")]/1000
    
    points(Year,Bio,col=CL[n],pch=19)
    segments(Year,Bio,Year,Bio+1.96*Bio.SE,col=CL[n])
    segments(Year,Bio,Year,Bio-1.96*Bio.SE,col=CL[n])
    Nms[n]=paste(names(d)[n]," (depletion=",round(d[[n]]$rep$depletion,2),")",sep="")
  }
  legend("topright",Nms,bty="n",pch=19,col=CL[1:length(d)],cex=1.25)
  mtext("Total biomass (tonnes)",2,las=3,line=2.5,cex=1.25)
  
  #Table of par estimates
  TABL=vector('list',length(d))
  names(TABL)=names(d)
  for(n in 1:length(TABL))
  {
      idx=fn.which(D=d,n=10,x="Total_biom")
      EST=d[[n]]$std[1:(idx[1]-1),match(c("name","value","std.dev"),names(d[[n]]$std))]
      EST=cbind(Scenario=names(TABL)[n],EST)
      TABL[[n]]=EST
    }
  TABL=do.call(rbind,TABL)
  TABL$MLE=paste(TABL$value," ±",TABL$std.dev,sep="")
  TABL=TABL[,-match(c("value","std.dev"),names(TABL))]
  #TABL <- reshape(TABL, v.names = "MLE", idvar = "Scenario",
   #               timevar = "name", direction = "wide")
  plot(1,ann=F,axes=F,col="transparent")
  addtable2plot(x.tab,y.tab,TABL,bty="o",display.rownames=F,hlines=F,vlines=TRUE,title="",cex=CEx.tab)
  
  #cpue fit
  CPUE.1=d[[1]]$rep$"CPUE(syr_cpue,nyr)"
  Est_CPUE.1=d[[1]]$rep$"Est_CPUE(syr_cpue,nyr)"
  Year=1975:(1975+length(CPUE.1)-1)
  
  plot(Year,CPUE.1,ylab="",ylim=c(0,max(c(CPUE.1,Est_CPUE.1))),xlim=c(1975,Last.yr),cex.lab=1.5,col="transparent")
  CL=c("black","red","blue","green","orange","gold")
  for(n in 1:length(d))
  {
    CPUE=d[[n]]$rep$"CPUE(syr_cpue,nyr)"
    Est_CPUE=d[[n]]$rep$"Est_CPUE(syr_cpue,nyr)"
    points(Year,CPUE,col=CL[n],pch=19)
    lines(Year,Est_CPUE,col=CL[n])
  }
  mtext("CPUE",2,las=3,line=2.5,cex=1.25)
  
  dev.off()
  
}

#S1 vs S1 pin sensitivity
fn.plt(d=LIST.mod[match(c("S1","S1_pin_sensitivity"),names(LIST.mod))],
       Compare.what="S1 vs S1_pin_sensitivity",CEx.tab=1.3,x.tab=.75,y.tab=0.5)


#S1 vs S1 no fishing efficiency
if(SPEC=="Whiskery")
{
  fn.plt(d=LIST.mod[match(c("S1","S1_no_fishing_efficiency"),names(LIST.mod))],
         Compare.what="S1 vs S1_no_fishing_efficiency",CEx.tab=1.3,x.tab=.75,y.tab=0.5)
}


#S1 vs S1_po_0.7 vs S1_po_0.99
fn.plt(d=LIST.mod[match(c("S1","S1_po_0.7","S1_po_0.99"),names(LIST.mod))],
       Compare.what="S1 vs S1_Po_sensitivity",CEx.tab=1,x.tab=.6,y.tab=.5)



