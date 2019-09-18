#script for creating input files and running stock assessment models of WA shark resources

#Steps: Section A: Bring in data (this brings in all the data series from "Organise data.R")
#       Section B: Bring in parameters (this brings in all parameters from "Organise parameters.R")
#       Section C: Put data together for scenarios
#       Section D: Run models
#       Section E: Display model outputs 
#       Section F: Run base case MCMC and display outputs
  #     Section G: RUN CATCH-msy method and display outputs
#       Section H: Dropped code


#source Handy function for plotting
source.hnld="C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/"
fn.source=function(script)source(paste(source.hnld,script,sep=""))
fn.source("fn.fig.R")
library(expm)
library(MASS)

if(First.run=="YES")
{
  set.seed(999)  #for reproducibility
  
# Section A: BRING IN INPUT DATA -------------------------------------------
  fn.source("Organise data.R")
  fn.input.data(SP=species,Yr.assess=AssessYr,Conv.cal.mn.to.fin.mn="NO",           
                Historic.Ktch="NO",Bin.size=TL.bins.cm,What.Efrt=What.Effort)  
  
  #Read in all input data files created by fn.input.data()
  if(file.exists(".Rhistory"))file.remove(".Rhistory") #remove the created '.Rhistory' from the folder to read in all .csv files
  All.input.data=list.files()
  data.list <- lapply(All.input.data, read.csv,stringsAsFactors=F)
  names(data.list)=unlist(lapply(All.input.data, function(x) unlist(strsplit(x, split='.csv', fixed=TRUE))))
  
  fn.add.mis.yr=function(what)
  {
    missing.yr=TC.TDGDLF$FINYEAR[which(!TC.TDGDLF$FINYEAR%in%what$FINYEAR)]
    add=what[1:length(missing.yr),]
    add$FINYEAR=missing.yr
    add[,names(add)[-match("FINYEAR",names(add))]]=0
    return(rbind(what,add))
  }
  fn.two.cols=function(what)
  {
    zn=unique(what$zone)
    what=what[,-match('zone',names(what))]
    names(what)[2]=zn
    return(what)
  }
  
  #A.1. Effort                        
    #A.1.1. TDGDLF
  Eff.tdgdlf=data.list$effort.annual.TDGDLF 
  Ef.zn=data.list$effort.annual.by.zone.TDGDLF
  Eff.wst=Ef.zn[match(c("FINYEAR","West"),names(Ef.zn))]
  Eff.zn1=Ef.zn[match(c("FINYEAR","Zone1"),names(Ef.zn))]
  Eff.zn2=Ef.zn[match(c("FINYEAR","Zone2"),names(Ef.zn))]
  
  #A.2. Catch                           
    #A.2.1 TDGDLF
  TC.TDGDLF=data.list$ktch.annual.TDGDLF
  kt=data.list$ktch.annual.by.zone.TDGDLF
  TC.TDGDLF.zn1=subset(kt,zone=="Zone1")
  TC.TDGDLF.zn2=subset(kt,zone=="Zone2")
  TC.TDGDLF.wst=subset(kt,zone=="West")
  
  TC.TDGDLF.zn1=fn.two.cols(TC.TDGDLF.zn1) #select only two columns
  TC.TDGDLF.zn2=fn.two.cols(TC.TDGDLF.zn2)
  TC.TDGDLF.wst=fn.two.cols(TC.TDGDLF.wst)
  
    #A.2.2. Other
  TC.other=data.list$ktch.annual.Other
  kt=data.list$ktch.annual.by.zone.Other
  TC.other=fn.add.mis.yr(TC.other)  #add missing years
  
  TC.other.zn1=subset(kt,zone=="Zone1")
  TC.other.zn2=subset(kt,zone=="Zone2")
  TC.other.wst=subset(kt,zone=="West")
  
  TC.other.zn1=fn.two.cols(TC.other.zn1) #select only two columns
  TC.other.zn2=fn.two.cols(TC.other.zn2)
  TC.other.wst=fn.two.cols(TC.other.wst)
  
  TC.other.zn1=fn.add.mis.yr(TC.other.zn1)  #add missing years
  TC.other.zn2=fn.add.mis.yr(TC.other.zn2)
  TC.other.wst=fn.add.mis.yr(TC.other.wst)
  
    #A.2.3. Rec fishing 
  #note: West Coast bioregion split in half West, half Zone 1
  TC.rec=data.list$ktch.annual.Rec
  kt=data.list$ktch.annual.by.zone.Rec
  
  TC.rec=subset(TC.rec,FINYEAR%in%TC.TDGDLF$FINYEAR)
  kt=subset(kt,FINYEAR%in%TC.TDGDLF$FINYEAR)
  
  TC.rec.zn2=subset(kt,zone=="South Coast")
  TC.rec.zn2$zone="Zone2"
  
  TC.rec.WC.bio=subset(kt,zone=="West Coast")
  n1=nrow(TC.rec.WC.bio)
  TC.rec.WC.bio=rbind(TC.rec.WC.bio,TC.rec.WC.bio)
  TC.rec.WC.bio$zone=c(rep("West",n1),rep("Zone1",n1))
  TC.rec.zn1=subset(TC.rec.WC.bio,zone=="Zone1")
  TC.rec.wst=subset(TC.rec.WC.bio,zone=="West")
  rm(kt)
  TC.rec.zn1=fn.two.cols(TC.rec.zn1) #select only two columns
  TC.rec.zn2=fn.two.cols(TC.rec.zn2)
  TC.rec.wst=fn.two.cols(TC.rec.wst)
  
    #A.2.4 combine  catches from all fisheries          
  fn.combo.ktch=function(LISTA)
  {
    dummy=LISTA[[1]]
    dat=do.call(cbind,LISTA)
    dat=dat[,-which(names(dat)=="FINYEAR")]
    dummy[,2]=rowSums(dat)
    return(dummy)
  }
  Ktch.All.1975=fn.combo.ktch(list(TC.TDGDLF,TC.other,TC.rec))
  Ktch.All.West.1975=fn.combo.ktch(list(TC.TDGDLF.wst,TC.other.wst,TC.rec.wst))
  Ktch.All.zn1.1975=fn.combo.ktch(list(TC.TDGDLF.zn1,TC.other.zn1,TC.rec.zn1))
  Ktch.All.zn2.1975=fn.combo.ktch(list(TC.TDGDLF.zn2,TC.other.zn2,TC.rec.zn2))
  
    #A.3. TDGDLF Standardised CPUE
  if('cpue.annual.TDGDLF.folly'%in%names(data.list))Cpue.folly=data.list$cpue.annual.TDGDLF.folly
  
      #Monthly
  Cpue.all=data.list$cpue.annual.TDGDLF
  Cpue.all.hours=data.list$cpue.annual.TDGDLF.hours
  Cpue.zns=data.list$cpue.annual.by.zone.TDGDLF
  Cpue.West=subset(Cpue.zns,zone=="West")
  Cpue.zn1=subset(Cpue.zns,zone=="Zone1")
  Cpue.zn2=subset(Cpue.zns,zone=="Zone2")
  rm(Cpue.zns)
  
    #Daily
  Cpue.all.daily=data.list$cpue.annual.TDGDLF.daily
  Cpue.all.daily.hours=data.list$cpue.annual.TDGDLF.daily.hours
  Cpue.zns.daily=data.list$cpue.annual.by.zone.TDGDLF.daily
  Cpue.West.daily=subset(Cpue.zns.daily,zone=="West")
  Cpue.zn1.daily=subset(Cpue.zns.daily,zone=="Zone1")
  Cpue.zn2.daily=subset(Cpue.zns.daily,zone=="Zone2")
  rm(Cpue.zns.daily)

    #A.4. TDGDLF Size composition
  #combined sexes
    #6.5 inch
  size.all=data.list$size.comp.annual.TDGDLF
  size.wst=data.list$size.comp.annual.WC.TDGDLF
  size.zn1=data.list$size.comp.annual.Zn1.TDGDLF
  size.zn2=data.list$size.comp.annual.Zn2.TDGDLF
  Nms=colnames(size.all)
  Nms[2:length(Nms)]=substr(Nms[2:length(Nms)],2,10)
  colnames(size.all)=colnames(size.wst)=colnames(size.zn1)=colnames(size.zn2)=Nms
  
    #7 inch
  size.all_7=data.list$size.comp.annual.TDGDLF_7
  size.wst_7=data.list$size.comp.annual.WC.TDGDLF_7
  size.zn1_7=data.list$size.comp.annual.Zn1.TDGDLF_7
  size.zn2_7=data.list$size.comp.annual.Zn2.TDGDLF_7
  colnames(size.all_7)=colnames(size.wst_7)=colnames(size.zn1_7)=colnames(size.zn2_7)=Nms
  
  #By sex
    #6.5 inch
      #females
  size.all.fem=data.list$size.comp.fem.annual.TDGDLF
  size.wst.fem=data.list$size.comp.annual.fem.WC.TDGDLF
  size.zn1.fem=data.list$size.comp.annual.fem.Zn1.TDGDLF
  size.zn2.fem=data.list$size.comp.annual.fem.Zn2.TDGDLF
  colnames(size.all.fem)=colnames(size.wst.fem)=colnames(size.zn1.fem)=colnames(size.zn2.fem)=Nms
  
      #males
  size.all.mal=data.list$size.comp.mal.annual.TDGDLF
  size.wst.mal=data.list$size.comp.annual.mal.WC.TDGDLF
  size.zn1.mal=data.list$size.comp.annual.mal.Zn1.TDGDLF
  size.zn2.mal=data.list$size.comp.annual.mal.Zn2.TDGDLF
  colnames(size.all.mal)=colnames(size.wst.mal)=colnames(size.zn1.mal)=colnames(size.zn2.mal)=Nms
  
    #7 inch
      #females
  size.all.fem_7=data.list$size.comp.fem.annual.TDGDLF_7
  size.wst.fem_7=data.list$size.comp.annual.fem.WC.TDGDLF_7
  size.zn1.fem_7=data.list$size.comp.annual.fem.Zn1.TDGDLF_7
  size.zn2.fem_7=data.list$size.comp.annual.fem.Zn2.TDGDLF_7
  colnames(size.all.fem_7)=colnames(size.wst.fem_7)=colnames(size.zn1.fem_7)=colnames(size.zn2.fem_7)=Nms
  
      #males
  size.all.mal_7=data.list$size.comp.mal.annual.TDGDLF_7
  size.wst.mal_7=data.list$size.comp.annual.mal.WC.TDGDLF_7
  size.zn1.mal_7=data.list$size.comp.annual.mal.Zn1.TDGDLF_7
  size.zn2.mal_7=data.list$size.comp.annual.mal.Zn2.TDGDLF_7
  colnames(size.all.mal_7)=colnames(size.wst.mal_7)=colnames(size.zn1.mal_7)=colnames(size.zn2.mal_7)=Nms
  
  
    #A.5 TDGDLF catch average weight 
  avg.wt=data.list$ktch.avg.weight.annual
  avg.wt.zn=data.list$ktch.avg.weight.annual.by.zone
  avg.wt.wst=subset(avg.wt.zn,zone=="West")
  avg.wt.zn1=subset(avg.wt.zn,zone=="Zone1")
  avg.wt.zn2=subset(avg.wt.zn,zone=="Zone2")
  rm(avg.wt.zn)
  
    #A.6. Conventional tagging  
  
  #individual based model inputs      
  
  #all 
  Conv.tg.rel_size=data.list$Zn.rel_Conv.Tag_size
  Conv.tg.rec_size=data.list$Zn.rec_Conv.Tag_size
  if(Move.mode=="Individual-based")
  {
    Conv.tg.rec.exp=data.list$Ind_based_model
    Conv.tg.rec.exp=subset(Conv.tg.rec.exp,DaysAtLarge>=MIN.DAYS.LARGE)
    Conv.tg.rec.exp=Conv.tg.rec.exp[,-match("Tag.no",names(Conv.tg.rec.exp))]
    names(Conv.tg.rec.exp)[match(c("Rel.zone","Rec.zone"),names(Conv.tg.rec.exp))]=c("Rel.zn","Rec.zn")
    Conv.tg.nzones=3
    Conv.tg.recs=nrow(Conv.tg.rec.exp)
  }
  if(Move.mode=="Population-based")
  {
    Minimum.number.released=10
    
    
    dropped=subset(Conv.tg.rel_size,Number<Minimum.number.released)
    dropped=dropped$TG.zn
    Conv.tg.rel_size=subset(Conv.tg.rel_size,Number>=Minimum.number.released)
    Conv.tg.rec_size=subset(Conv.tg.rec_size,!TG.zn%in%dropped)
    
    Conv.tg.rel_plot=Conv.tg.rel_size
    Conv.tg.rec_plot=Conv.tg.rec_size
    
    Tag.groups=unique(Conv.tg.rel_size$TG.zn)
    
    Conv.tg.StartYear=min(Conv.tg.rel_size$Yr.rel)
    Conv.tg.EndYear=max(Conv.tg.rec_size$Yr.rec)
    Conv.tg.nyrs=length(Conv.tg.StartYear:Conv.tg.EndYear)
    Conv.tg.Areas=length(unique(Conv.tg.rel_size$Rel.zone))
    Conv.tg.numTagGp=nrow(Conv.tg.rel_size)
    
    Conv.tg.rel_size$Rel.zone=factor(Conv.tg.rel_size$Rel.zone)
    Conv.tg.rel_size=Conv.tg.rel_size[order(Conv.tg.rel_size$Yr.rel,Conv.tg.rel_size$Rel.zone),]   
    Conv.tg.rel_size=reshape(Conv.tg.rel_size,v.names = "Number", idvar = c("TG.zn","Yr.rel"),
                             timevar = "Rel.zone", direction = "wide")
    
    Conv.tg.Releases_yr=as.matrix(Conv.tg.rel_size$Yr.rel)
    Conv.tg.rel_size=Conv.tg.rel_size[,-match(c("TG.zn","Yr.rel"),names(Conv.tg.rel_size))]
    Conv.tg.rel_size[is.na(Conv.tg.rel_size)]=0
    
    Conv.tg.rec_size$Rec.zone=factor(Conv.tg.rec_size$Rec.zone)
    Conv.tg.rec_size=Conv.tg.rec_size[order(Conv.tg.rec_size$TG.zn,Conv.tg.rec_size$Yr.rec,Conv.tg.rec_size$Rec.zone),]   
    Conv.tg.rec_size=reshape(Conv.tg.rec_size,v.names = "Number", idvar = c("TG.zn","Yr.rec"),
                             timevar = "Rec.zone", direction = "wide")
    
    #fill in no observations to have square matrices
    dummy=vector('list',Conv.tg.numTagGp)
    dummy1=dummy
    Conv.tg.all.yrs=Conv.tg.StartYear:Conv.tg.EndYear
    for(i in 1:Conv.tg.numTagGp)
    {
      a=subset(Conv.tg.rec_size,TG.zn==Tag.groups[i])
      id=Conv.tg.all.yrs[which(!Conv.tg.all.yrs%in%a$Yr.rec)]
      nn=nrow(a)
      a=rbind(a,a[(nn+1):(nn+length(id)),])
      a$TG.zn=a$TG.zn[1]
      a$Yr.rec[(nn+1):nrow(a)]=id
      a=a[order(a$Yr.rec),]
      a[is.na(a)]=0
      a=a[,-match(c("TG.zn","Yr.rec"),names(a))]
      dummy[[i]]=a
      
      #get index
      a[a>0]=1
      dummy1[[i]]=a
      
    }
    
    Conv.tg.rec_size=do.call(rbind,dummy)
    Conv.tg.rec_size_index=do.call(rbind,dummy1)
    Conv.tg.mov.pars=Conv.tg.Areas^2-Conv.tg.Areas
    Conv.tg.Like.type=3   #3: negative binomial, 2: Poisson, 1:normal
    Conv.tg.constant=0.00001
    
  }
  
  #juveniles and adults
  Conv.tg.rel_size_adul=data.list$Zn.rel_Conv.Tag_size_adu
  Conv.tg.rel_size_juv=data.list$Zn.rel_Conv.Tag_size_juv
  Conv.tg.rec_size_adul=data.list$Zn.rec_Conv.Tag_size_adu
  Conv.tg.rec_size_juv=data.list$Zn.rec_Conv.Tag_size_juv
  Smallest_size_tagged=data.list$Smallest_size_tagged[1,1]
  
    #A.7. Acoustic tagging  
  if(Acoust.format=="Individual-based")
  {
    Acous.tg.rec=data.list$acous.tag.reca.indv.based
    Acous.nzones=3
    Acous.detections=nrow(Acous.tg.rec)
    Acous.ntags=length(unique(Acous.tg.rec$TagID))
    Acous.TagID_observations=c(table(Acous.tg.rec$TagID))
    names(Acous.TagID_observations)=NULL
    Acous.tg.rec$Index=1:nrow(Acous.tg.rec)
    Acous.TagID_start=do.call("rbind",by(Acous.tg.rec,Acous.tg.rec$TagID,head,1))$Index
    Acous.TagID_end=do.call("rbind",by(Acous.tg.rec,Acous.tg.rec$TagID,tail,1))$Index
    Acous.tg.rec$TagID=as.numeric(substr(Acous.tg.rec$TagID,4,6))
    Acous.tg.rec=subset(Acous.tg.rec,select=c(DaysAtLarge,Rel.zn,Rec.zn,TagID))
  }
  if(Acoust.format=="SS3")
  {
    Type.acous.dat='number'
    #Type.acous.dat='proportion'
    
    if(Type.acous.dat=='number') 
    {
      Acous.tg.rel=data.list$acous.tag.rel
      Acous.tg.rec=data.list$acous.tag.reca
    }
    
    if(Type.acous.dat=='proportion') 
    {
      Acous.tg.rel=data.list$acous.tag.rel.prop
      Acous.tg.rec=data.list$acous.tag.reca.prop
    }
  }
  

    #A.8. Age and growth data 
  Age.growth=data.list$Age.growth

  
    #A.9. Mesh proportional effort
  Mesh.prop.eff=data.list$Mesh.prop.eff
  Mesh.prop.eff.West=data.list$Mesh.prop.eff.West
  Mesh.prop.eff.Zn1=data.list$Mesh.prop.eff.Zn1
  Mesh.prop.eff.Zn2=data.list$Mesh.prop.eff.Zn2
  
  
  Mesh.prop.eff=subset(Mesh.prop.eff,finyear%in%unique(Ktch.All.1975$FINYEAR))
  Mesh.prop.eff.West=subset(Mesh.prop.eff.West,finyear%in%unique(Ktch.All.1975$FINYEAR))
  Mesh.prop.eff.Zn1=subset(Mesh.prop.eff.Zn1,finyear%in%unique(Ktch.All.1975$FINYEAR))
  Mesh.prop.eff.Zn2=subset(Mesh.prop.eff.Zn2,finyear%in%unique(Ktch.All.1975$FINYEAR))
  
  Yrs.sel=Mesh.prop.eff$finyear
  
  
    #add future selectivity for future projections
  fn.add.future.sel=function(d)
  {
    Add.future=d[rep(nrow(d),Yrs.future),]
    lst.yr=d$finyear[nrow(d)]
    y1=as.numeric(substr(lst.yr,1,4))+1
    y2=as.numeric(substr(lst.yr,6,7))+1
    Add.future$finyear=paste(y1:(y1+Yrs.future-1),y2:(y2+Yrs.future-1),sep="-")
    return(rbind(d,Add.future))
  }
  Mesh.prop.eff=fn.add.future.sel(Mesh.prop.eff)
  Mesh.prop.eff.West=fn.add.future.sel(Mesh.prop.eff.West)
  Mesh.prop.eff.Zn1=fn.add.future.sel(Mesh.prop.eff.Zn1)
  Mesh.prop.eff.Zn2=fn.add.future.sel(Mesh.prop.eff.Zn2)
  
    #rm list after extracting data
  rm(data.list)
    
  

# Section B: BRING IN PARAMETERS ------------------------------------------
  hndl=paste("C:/Matias/Analyses/Population dynamics/1.",Spec," shark/",sep='')
  
    #B.1 Source all input parameters
  fn.source("Organise input parameters.R")
  ParS=fn.input.pars(SP=species,add.growth.cv="NO",add.Mrt.age="NO")$pars
  
  #Maximum age
  Max.age.F=ParS$Max.Age.F
  Max.age.M=ParS$Max.Age.M
  Max.age=Max.age.F
  
  #Natural mortality
  M=ParS$M
  
  
  #TL-FL
  TL.FL=ParS$TL.to.FL
  
  #Min and max size in population
  Min.Max.FL=ParS$Min.Max.FL
  Min.size.FL=Min.Max.FL$Min 
  Max.size.FL=Min.Max.FL$Max
  
  Min.size.TL=TL.FL[1]*Min.size.FL+TL.FL[2]
  Max.size.TL=TL.FL[1]*Max.size.FL+TL.FL[2]
  Min.size.TL=floor(Min.size.TL/10)*10
  Max.size.TL=floor(Max.size.TL/10)*10
  
  #Growth           
  Growth.F=ParS$Growth.F
  Growth.M=ParS$Growth.M
  
  if(is.na(match("FL_inf",names(Growth.F))))
  {
    Growth.F$FL_inf=(Growth.F$TL_inf-TL.FL$b)/TL.FL$a
    Growth.M$FL_inf=(Growth.M$TL_inf-TL.FL$b)/TL.FL$a
  }
  if(!"TL_inf"%in%names(Growth.M)) 
  {
    Growth.M$TL_inf=TL.FL$a*Growth.M$FL_inf+TL.FL$b
    Growth.F$TL_inf=TL.FL$a*Growth.F$FL_inf+TL.FL$b
  }
  
  
  #TL-TWT
  Tl_Twlt.fem=ParS$TL.to.TwT.F
  Tl_Twlt.mal=ParS$TL.to.TwT.M
  if(!is.data.frame(Tl_Twlt.mal)) Tl_Twlt.mal=data.frame(a=Tl_Twlt.mal[1],b=Tl_Twlt.mal[2])
  
  #Fecundity
  Fecundity=ParS$Litter.sz
  Fecundity.min=Fecundity$Min
  Fecundity.max=Fecundity$Max
  Fecundity.at.size=ParS$Litter.sz.at.size
  
  #Breeding frequency
  Breed.freq=ParS$Breed.freq
  Breed.freq.min=Breed.freq$Min
  Breed.freq.max=Breed.freq$Max
  
  #Maturity
  Age.50.mat=ParS$Age.50.mat
  Age.50.mat.min=Age.50.mat$Min
  Age.50.mat.max=Age.50.mat$Max
  Maturity.at.size=ParS$Mat.50.95
  
  #Pups sex ratio
  pup.sx.ratio=ParS$Sex.ratio
  
  #Proportion of males in catch
  Prop.males.ktch=ParS$Prop.males.in.ktch
  Prop.males.ktch.wst=Prop.males.ktch$WC  #these are used for by zone model
  Prop.males.ktch.zn1=Prop.males.ktch$Zn1
  Prop.males.ktch.zn2=Prop.males.ktch$Zn2
  Prop.males.ktch=Prop.males.ktch$All     #this is used for single zone model
  if(Spec=="Gummy") Prop.males.ktch.zn1=mean(c(Prop.males.ktch.wst,Prop.males.ktch.zn2))    #Zn2 fishers input from Esperance meeting May 2017
  
  #Size at birth (FL)
  Size.birth=ParS$Size.birth
  
  #Gillnet selectivity
  Selectivity=ParS$Selectivity
  Selectivity_7=ParS$Selectivity_7
  
  #Tagg shedding
  Shedding=ParS$Shedding
  
  #Tag reporting
  Reporting=ParS$Reporting 
  if(Move.mode=="Population-based")
  {
    Reporting$YR=as.numeric(substr(Reporting$FinYear,1,4))+1
    a=lm(Zn1~YR,Reporting)
    Rep.predZn1=predict(a,type='response',newdata=data.frame(YR=Conv.tg.all.yrs))
    a=lm(Zn2~YR,Reporting)
    Rep.predZn2=predict(a,type='response',newdata=data.frame(YR=Conv.tg.all.yrs))
    a=lm(WC~YR,Reporting)
    Rep.predWC=predict(a,type='response',newdata=data.frame(YR=Conv.tg.all.yrs))
    Reporting.pred=data.frame(YR=Conv.tg.all.yrs,WC.pred=Rep.predWC,ZN1.pred=Rep.predZn1,ZN2.pred=Rep.predZn2)
    
    Reporting=merge(Reporting.pred,Reporting,by="YR",all.x=T)
    Reporting$Zn1=with(Reporting,ifelse(is.na(Zn1),ZN1.pred,Zn1))
    Reporting$Zn2=with(Reporting,ifelse(is.na(Zn2),ZN2.pred,Zn2))
    Reporting$WC=with(Reporting,ifelse(is.na(WC),WC.pred,WC))
    Reporting_plot=Reporting[,-match("FinYear",names(Reporting))]
    Reporting=Reporting[,match(c("WC","Zn1","Zn2"),names(Reporting))]
    Reporting[Reporting<0.1]=0.1
    Reporting_plot[Reporting_plot<0.1]=0.1
    # plot(Reporting$YR,Reporting$Zn2,pch=19,xlim=c(1994,2014),ylim=c(0,1))
    #  lines(Reporting.pred$YR,Reporting.pred$ZN2.pred)  
  }
  
  #Steeepnes
  Steepness=ParS$STEEP
  
  Hndl="C:/Matias/Data/Population dynamics/Parameter inputs for models/"
  Mortality.from.Ref.Point=read.csv(paste(Hndl,Sp2,".M_at_age.csv",sep=""))
  
  
   #B.2. Create input pars table
  fn.source("Table.input.pars.R")
  TBl=fn.par.ref.tbl(species,add.growth.cv="NO",add.Mrt.age="NO")
  
  #select parameters of interest
  Size.birth_SD=TBl[1,]
  rownames(Size.birth_SD)=Size.birth_SD$Parameter="Size.birth_SD"
  Size.birth_SD$Value=5
  Size.birth_SD$Comment=""
  Size.birth_SD$Source='assumed'
  TBl=rbind(TBl,Size.birth_SD)
  
  Smal_siz_tag=TBl[1,]
  rownames(Smal_siz_tag)=Smal_siz_tag$Parameter="Smallest_size_tagged"
  Smal_siz_tag$Value=Smallest_size_tagged
  Smal_siz_tag$Comment=""
  Smal_siz_tag$Source='conv tagging observations'
  TBl=rbind(TBl,Smal_siz_tag)
  
  Mean.size.birth.TL=Min.size.TL
  if(species=="WH") Mean.size.birth.TL=25
  if(species=="GM") Mean.size.birth.TL=33
  
  Basic.pars=c("TL.to.FL.1","TL.to.FL.2","TL.to.TwT.F.2","TL.to.TwT.F.1","TL.to.TwT.M.2","TL.to.TwT.M.1","Min.Max.FL.2",
               "Max.Age.M","Max.Age.F","Growth.F.1","Growth.F.2","Growth.F.3","Growth.F.4","Growth.M.1","Growth.M.2","Growth.M.3",
               "Breed.freq.1","Size.birth","Size.birth_SD","Mat.50.95.1","Mat.50.95.2","Age.50.mat.1","Litter.sz.1","Litter.sz.2",
               "Litter.sz.at.size.1","Litter.sz.at.size.2","Sex.ratio","Selectivity.1","Selectivity.2",
               "Selectivity_7.1","Selectivity_7.2",
               "Prop.males.in.ktch.3","Prop.males.in.ktch.1","Prop.males.in.ktch.2","Prop.males.in.ktch.4",
               "M","STEEP.1")
  if(Move.mode=="Population-based") Conv.T.pars="Shedding"
  if(Move.mode=="Individual-based") Conv.T.pars="Smallest_size_tagged"  
  These.rows=Basic.pars
  if(add.conv.tag=="YES")  These.rows=c(Basic.pars,Conv.T.pars)
  
  TBl=TBl[match(These.rows,rownames(TBl)),]
  TBl$Value[match("Size.birth",TBl$Parameter)]=Mean.size.birth.TL
  TBl$Value[match("Min.Max.FL.Max",TBl$Parameter)]=Max.size.TL
  
  Stipnes=as.numeric(TBl$Value[[which(TBl$Parameter=="STEEP.mean")]])
  #TBl=TBl[-which(TBl$Parameter=="STEEP.mean"),]
  
  # if(Spec=="Whiskery") TBl=TBl[-match(c("Growth.F.k","Growth.F.FL_inf","Growth.F.to",
  #       "Growth.F.SD","Growth.M.k","Growth.M.FL_inf","Growth.M.to"),TBl$Parameter),] #remove growth input pars as they are estimated
  
  setwd(paste(hndl,AssessYr,"/1_Inputs",sep=""))
  fn.word.table(WD=getwd(),TBL=TBl,Doc.nm="Input pars",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
  rm(ParS)
  


# Section C: PUT DATA TOGETHER FOR SCENARIOS ------------------------------
  
  #C.1 Construct relationships at TL and at age
  fn.plt.rel=function(one,two,add,three,four,CL,CL1)
  {
    YLIM=c(0,max(two)*1.15)
    if(add=="YES") YLIM=c(0,max(c(two,four)))
    plot(one,two,xlab="",ylab="",type='l',lwd=2,col=CL,ylim=YLIM)
    if(add=="YES") lines(three,four,lwd=2,col=CL1)
  }
  fn.litter=function(L,a,b,MaxLitter)
  {
    if(species=="WH") litter=L*a+b
    if(species=="GM")  litter=exp(b+a*L)
    litter=ifelse(litter>MaxLitter,MaxLitter,litter)   #cap to within observations
    return(litter)
  }
  fn.mature=function(L,L50,L95) 1/(1+exp(-log(19)*(L-L50)/(L95-L50)))
  fn.sel=function(L,alpha,beta) ((L/(alpha*beta))^alpha)*(exp(alpha-(L/beta)))  
  Age.from.TL=function(TL,Linf,k,to)
  {
    TL=ifelse(TL>Linf,Linf*.99,TL)
    Age=-(log(1-(TL/Linf))/k)+to
    Age=ifelse(Age<=0,0,Age)
    return(Age)
  }
  TL.from.age=function(age,Linf,k,to) Linf*(1-exp(-k*(age-to)))  
  
    #C.1.1 At length relationships
      #C.1.1.1. set length vector
  if(MN.SZE=="size.at.birth") Min.size.bin=floor(Min.size.TL$a)
  if(MN.SZE==0) Min.size.bin=0
  TL.bin.low=seq(Min.size.bin,Max.size.TL$a*Plus.gp.size,by=5)        
  TL.bin.mid=TL.bin.low+(TL.bins.cm/2)     
  TL.bin.up=TL.bin.low+TL.bins.cm
    
      #C.1.1.2. at TL relationships
  #TL and FL (cm) (mid year)
  FL.bin.mid=(TL.bin.mid-TL.FL$b)/TL.FL$a
  
  #size at birth in TL
  Size.birth.TL=TL.FL[1]*Size.birth+TL.FL[2]
  
  #Twt (kg)
  TwT.bin.mal=Tl_Twlt.mal$b*TL.bin.mid^Tl_Twlt.mal$a
  TwT.bin.fem=Tl_Twlt.fem$b*TL.bin.mid^Tl_Twlt.fem$a
  
  #M
  #note: obtained from Reference point paper using a spline for interpolating ages derived
  #       from size bins
  library(splines)
  age=Age.from.TL(TL.bin.mid,Growth.F$TL_inf,Growth.F$k,Growth.F$to)
  age=ifelse(age>max(Mortality.from.Ref.Point[,2]),max(Mortality.from.Ref.Point[,2]),age)
  Mortality.from.Ref.Point$dif=c(-1,diff(Mortality.from.Ref.Point$Median.M))
  Mortality.from.Ref.Point=subset(Mortality.from.Ref.Point,dif<0)
  Spline.fn=splinefun(Mortality.from.Ref.Point$Age,Mortality.from.Ref.Point$Median.M, method = "natural")
  Mort.bin=Spline.fn(age)
  
  #Maturity
  Mat.bin=fn.mature(TL.bin.mid,Maturity.at.size[1],Maturity.at.size[2])
    
  #Litter size
  if(species=="WH")Litter.size.bin=fn.litter(FL.bin.mid,Fecundity.at.size[1],Fecundity.at.size[2],Fecundity.max)
  if(species=="GM")Litter.size.bin=fn.litter(TL.bin.mid,Fecundity.at.size[1],Fecundity.at.size[2],Fecundity.max)
  Litter.size.bin=Litter.size.bin*Mat.bin
  Litter.size.bin[Litter.size.bin<0]=0
  
  #Selectivity (selectivity was calculated in mm)
  Select.bin=fn.sel(FL.bin.mid*10,Selectivity$alpha,Selectivity$beta)  
  Select.bin=ifelse(is.na(Select.bin),0,Select.bin)
  
  Select.bin_7=fn.sel(FL.bin.mid*10,Selectivity_7$alpha,Selectivity_7$beta)  
  Select.bin_7=ifelse(is.na(Select.bin_7),0,Select.bin_7)
  
  #Selectivity by year-region-mesh
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
      Str[,yyy]=Sel.sum/max(Sel.sum)
    }
    return(Str)
  }  
  Total.sel.all=fn.varying.sel(Sel_6.5=Select.bin,Sel_7=Select.bin_7,TIME=Mesh.prop.eff)    
  
  Total.sel.West=fn.varying.sel(Sel_6.5=Select.bin,Sel_7=Select.bin_7,TIME=Mesh.prop.eff.West)    
  Total.sel.Zn1=fn.varying.sel(Sel_6.5=Select.bin,Sel_7=Select.bin_7,TIME=Mesh.prop.eff.Zn1)
  Total.sel.Zn2=fn.varying.sel(Sel_6.5=Select.bin,Sel_7=Select.bin_7,TIME=Mesh.prop.eff.Zn2)
  

  
  #Age and growth
  if(Spec=="Whiskery") Age.growth$TL=TL.FL$a*Age.growth$FL+TL.FL$b #calculate TL for whiskery
  Age.growth.F=subset(Age.growth,Sex=="Female",select=c(Counts,TL))
  Age.growth.M=subset(Age.growth,Sex=="Male",select=c(Counts,TL))
  
  #Size transition matrix
  fn.source("size transition.R")
  STM.F=Sadovy.size.trans.mat(Linf=Growth.F$TL_inf,K=Growth.F$k,SD=Growth.F$SD,bin.low=TL.bin.low,
                              bin.up=TL.bin.up,Truncate="YES")
  STM.M=Sadovy.size.trans.mat(Linf=Growth.M$TL_inf,K=Growth.M$k,SD=Growth.M$SD,bin.low=TL.bin.low,
                              bin.up=TL.bin.up,Truncate="YES")
  if(MN.SZE==0) int=2
  if(MN.SZE=="size.at.birth") int=3 
  
  #Test size transition matrix
  CLS=rainbow(ncol(STM.F))
  fn.fig(paste(getwd(),"/Visualise data/Input.STM.F",sep=""),2400,2400) 
  plot(TL.bin.low,STM.F[,1],ylim=c(0,.75),type='l')
  for(i in 2:ncol(STM.F)) lines(TL.bin.low,STM.F[,i],col=CLS[i])
  dev.off()

  #show growth data       
  fn.fig(paste(getwd(),"/Visualise data/Growth_data",sep=""),2400,2400)  
  par(xpd=T,las=1,mgp=c(2.5,.7,0))
  with(subset(Age.growth,Sex=="Male"),plot(Counts, TL,ylab="TL (cm)",xlab="Age",pch=19,col="blue",
        cex=1.5,xlim=c(0,max(Age.growth$Counts)*1.05),ylim=c(0,max(Age.growth$TL)*1.05),
        cex.axis=1.5,cex.lab=2,yaxs="i",xaxs="i"))
  with(subset(Age.growth,Sex=="Female"),points(Counts+.05, TL,pch=19,col="pink",cex=1.5))
  points(0,Mean.size.birth.TL,pch=19,cex=1.5)
  text(3.75,Mean.size.birth.TL,"Size at birth",offset=1,cex=1.75)
  arrows(2,Mean.size.birth.TL,0.25,Mean.size.birth.TL,length=0.15)
  dev.off()
  

    #C.1.2 At age relationships
  
      #C.1.2.1. set age vector
  age.M=0:Max.age.M
  age.F=0:Max.age.F
  
      #C.1.2.2. constant relationships
  #M
  Mort.M.k=rep(M,length(age.M))
  Mort.F.k=rep(M,length(age.F))
  
  #Maturity
  Mat.k=ifelse(age.F<Age.50.mat.min,0,1)
  
  #Litter size
  Litter.size.k= mean(c(Fecundity.min,Fecundity.max)) 
  if(species=="WH") Litter.size.k=19
  FecU=Litter.size.k
  Litter.size.k=Litter.size.k*Mat.k
  
  
      #C.1.2.3. at age relationships
  #TL and FL (cm) (mid year)
  mid.TL.M=Growth.M$TL_inf*(1-exp(-Growth.M$k*(age.F+0.5-Growth.M$to)))  #use same number of ages
  mid.TL.F=Growth.F$TL_inf*(1-exp(-Growth.F$k*(age.F+0.5-Growth.F$to)))
  mid.FL.F=(mid.TL.F-TL.FL$b)/TL.FL$a
  mid.FL.M=(mid.TL.M-TL.FL$b)/TL.FL$a
  
  #Twt (kg)
  mid.TwT.M=Tl_Twlt.mal$b*mid.TL.M^Tl_Twlt.mal$a
  mid.TwT.F=Tl_Twlt.fem$b*mid.TL.F^Tl_Twlt.fem$a
  
  #M  
  Mort.M.age=Spline.fn(age.M)
  Mort.F.age=Spline.fn(age.F)
  
  #Maturity
  Mat.age=fn.mature(mid.TL.F,Maturity.at.size[1],Maturity.at.size[2])
  Mat.age=Mat.age/max(Mat.age)
  
  #Litter size
  Whis.FL.fec=mid.FL.F    #correct artefact of low Linf to allow for observed litter size 
  Whis.FL.fec=Whis.FL.fec/max(Whis.FL.fec)*Max.FL.obs
  if(species=="WH")Litter.size.age=fn.litter(Whis.FL.fec,Fecundity.at.size[1],Fecundity.at.size[2],Fecundity.max)
  if(species=="GM")Litter.size.age=fn.litter(mid.TL.F,Fecundity.at.size[1],Fecundity.at.size[2],Fecundity.max)
  FecU_age=Litter.size.age
  
    #from at length to at age
  if(Do.from.at.len.to.at.age=="YES")
  {
    fn.from.at.length.to.at.age=function(VAR,Age.grow.F,Age.grow.M)
    {
      fit.fn=function(pars,Obs,AGE)
      {
        Par=exp(pars)
        if(VAR=="litter") Pred=Par[1]*(1-exp(-Par[2]*(AGE-Par[3])))    
        if(VAR=="weight") Pred=Par[1]*(1-exp(-Par[2]*AGE))^Par[3] 
        
        if(length(Obs)==length(Pred))
        {
          epsilon = Obs-Pred  
          nloglike=-1.0*sum(dnorm(epsilon,0,Par[4]))
          return(list(nloglike=nloglike,Pred=Pred,epsilon=epsilon))  
        }else  return(list(Pred=Pred))      
      }
      
      #length-litter size to age-litter size
      VAR="litter"
      pars=c(a=log(25),k=log(0.4),b=log(1.5),SD=log(5))
      
      if(species=="WH")Litter.age.F.obs=fn.litter(Age.grow.F$FL,Fecundity.at.size[1],Fecundity.at.size[2],Fecundity.max)
      if(species=="GM")Litter.age.F.obs=fn.litter(Age.grow.F$TL,Fecundity.at.size[1],Fecundity.at.size[2],Fecundity.max)
      
      fn=function(pars) fit.fn(pars,Litter.age.F.obs,Age.grow.F$Counts)$nloglike
      fit = optim(pars,fn,method="BFGS",hessian=T)
      Pred.litter.age.F=fit.fn(fit$par,Litter.age.F.obs,age.F)$Pred
      plot(Age.grow.F$Counts,Litter.age.F.obs,pch=19)
      lines(age.F,Pred.litter.age.F,col=2,lwd=3)
      
      #length-weight to age-weight
      VAR="weight"
      pars=c(a=log(9),k=log(0.7),b=log(5),SD=log(5))
      
      #females
      Wt.age.F.obs=Tl_Twlt.fem$b*Age.grow.F$TL^Tl_Twlt.fem$a
      fn=function(pars) fit.fn(pars,Wt.age.F.obs,Age.grow.F$Counts)$nloglike
      fit = optim(pars,fn,method="BFGS",hessian=T)
      Pred.wei.age.F=fit.fn(fit$par,Wt.age.F.obs,age.F)$Pred
      
      plot(Age.grow.F$Counts,Wt.age.F.obs,pch=19)
      lines(age.F,Pred.wei.age.F,col=2,lwd=3)
      
      
      return(list(Litter.age=Pred.litter.age.F,Wei.Age=Pred.wei.age.F))
    }
    
    STRE=fn.from.at.length.to.at.age(Age.grow.F=subset(Age.growth,Sex=="Female"),
                                     Age.grow.M=subset(Age.growth,Sex=="Male"))
    mid.TwT.F=STRE$Wei.Age
    Litter.size.age=STRE$Litter.age
  }
  
  Litter.size.age=Litter.size.age*Mat.age
  Litter.size.age[Litter.size.age<0]=0
   
  #Selectivity (selectivity was calculated in mm)
  Select.age.M=fn.sel(mid.FL.M*10,Selectivity$alpha,Selectivity$beta)  
  Select.age.F=fn.sel(mid.FL.F*10,Selectivity$alpha,Selectivity$beta)
  
  
  #D2. Display relationships  
  fn.mtxt=function(TxT) mtext(TxT,2,2,las=3,cex=0.85)
  if(BaseCase=="Size-based")
  {
    fn.fig("Visualise data/Input relations",2000, 2600)  
    
    par(mfcol=c(6,2),mai=c(.15,.125,.1,.195),oma=c(3,2.5,.25,.1),las=1,mgp=c(2,.6,0))
    
    #1. At length relationships
    
    #Size transition matrix
    model.estimated="YES"
    if(model.estimated=="NO")
    {
      f <- function(m) t(m)[,nrow(m):1]  #function for rotating matrix
      STM.t=f(STM.F)
      N.int=50
      colfunc <- colorRampPalette(c("navy", "cadetblue","white"))
      couleurs=rev(colfunc(N.int))
      BREAKS=seq(0,1,length.out=N.int+1)
      
      library(plotrix)
      xx=1:nrow(STM.t)
      yy=1:ncol(STM.t)
      image(xx,yy,STM.t,ylab="",xlab="",xaxt='n',
            yaxt='n',col =couleurs,breaks=BREAKS)
      axis(3,xx,F,tck=0.025)
      axis(2,yy,F,tck=0.025)
      axis(3,seq(1,nrow(STM.t),5),colnames(STM.F)[seq(1,nrow(STM.t),5)],tck=0.05,padj=0.65)
      SS=1+nrow(STM.t)-seq(1,nrow(STM.t),5)
      
      axis(2,SS,rev(colnames(STM.F)[SS-int]),tck=0.05,hadj=.8)
      fn.mtxt("To TL bin class (cm)")
      mtext("From TL bin class (cm)",3,1.25,cex=0.85)
      box()
      SQ=rev(seq(BREAKS[1],BREAKS[length(BREAKS)],.25))
      color.legend(xx[length(xx)*.9],yy[length(yy)*.35],xx[length(xx)*.99],yy[length(yy)*.99],
                   SQ,rect.col=rev(couleurs),gradient="y",col=1,cex=.7)    
    }
    
    #Growth
    Growth.size.based=NA
    
    plot(1,1,xaxt='n',yaxt='n',ann=F,col="transparent")
    text(1,1,"Model estimated",cex=2)
    
    fn.mtxt("TL (cm)")
    
    #TL-TwT
    fn.plt.rel(TL.bin.mid,TwT.bin.fem,add="YES",TL.bin.mid,TwT.bin.mal,"pink",4)
    fn.mtxt("Mid Twt (kg)")
    
    #TL- female maturity
    ADD="NO"
    if(length(which(Tabla.scen$Model_type=="Length-based" & Tabla.scen$Maturity=="knife edge"))>0) ADD="YES"
    fn.plt.rel(TL.bin.mid,Mat.bin,add=ADD,TL.bin.mid,ifelse(TL.bin.mid>=Maturity.at.size[1],1,0),"grey20","grey70")
    fn.mtxt("Prop. mat.")
    legend("topleft","Female",bty='n',cex=1.25)
    if(ADD=="YES")legend("left",c("at length","knife-edge"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)
    
    
    #TL-litter size * proportion mature
    Liter.len=NA
    if(is.na(Liter.len))
    {
      plot(1,1,xaxt='n',yaxt='n',ann=F,col="transparent")
      text(1,1,"Not applicable",cex=2)
    }
    if(!is.na(Liter.len))
    {
      fn.plt.rel(TL.bin.mid,Litter.size.bin,add="YES",
                 TL.bin.mid,ifelse(TL.bin.mid>=Maturity.at.size[1],Litter.size.k,0),"grey20","grey70")
    }
    fn.mtxt("Litter size")
    
    #TL-mortality
    ADD="NO"
    if(length(which(Tabla.scen$Model_type=="Length-based" & Tabla.scen$M=="at length"))>0) ADD="YES"
    fn.plt.rel(TL.bin.mid,rep(unique(Mort.F.k),length(TL.bin.mid)),add=ADD,TL.bin.mid,Mort.bin,"grey20","grey70")
    fn.mtxt(expression("M  " (year^-1)))
    if(ADD=="YES")legend("bottomleft",c("constant","at length"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)
    
    #TL-selectivity 
    fn.plt.rel(TL.bin.mid,Select.bin,add="NO",Dum,Dum,1,4)  
    fn.mtxt("Rel. selectivity")
    mtext("mid TL (cm)",1,2)
    
    
    #2.Relationships at age
    
    #Age-TL
    fn.plt.rel(age.F,mid.TL.F,add="YES",age.F,mid.TL.M,"pink",4)
    legend("bottomright",c("F","M"),lty=1,lwd=2,col=c("pink","blue"),bty='n',cex=1.25)
        
    #Age-TwT
    fn.plt.rel(age.F,mid.TwT.F,add="YES",age.F,mid.TwT.M,"pink",4)
        
    #Age- female maturity
    ADD="NO"
    if(length(which(Tabla.scen$Model_type=="Age-structured" & Tabla.scen$Maturity=="at age"))>0) ADD="YES"
     fn.plt.rel(age.F,Mat.k,add=ADD,age.F,Mat.age,"grey20","grey70")
    if(ADD=="YES")legend("bottomright",c("knife-edge","at age"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)
    legend("topleft","Female",bty='n',cex=1.25)
    
    
    #Litter size
    ADD="NO"
    if(length(which(Tabla.scen$Model_type=="Age-structured" & Tabla.scen$Fec.=="at age"))>0) ADD="YES"
    fn.plt.rel(age.F,FecU*Mat.k,add=ADD,age.F,FecU_age,"grey20","grey70")
    if(ADD=="YES")legend("bottomright",c("knife-edge","at age"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)
    
    
    #Age-female mortality
    ADD="NO"
    if(length(which(Tabla.scen$Model_type=="Age-structured" & Tabla.scen$M=="at age"))>0) ADD="YES"
    fn.plt.rel(age.F,Mort.F.k,add=ADD,age.F,Mort.F.age,"grey20","grey70")   
    if(ADD=="YES")legend("bottomright",c("constant","at age"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)  
    
    #Age-selectivity 
    fn.plt.rel(age.F,Select.age.F,add="YES",age.F,Select.age.M,"pink",4)
    mtext("Age",1,2) 
    
    dev.off()
  }
  
  if(BaseCase=="Age-based")
  {
    fn.fig("Visualise data/Input relations",1000, 2400)  
    par(mfcol=c(6,1),mai=c(.15,.125,.1,.15),oma=c(3,2.5,.25,.1),las=1,mgp=c(2,.6,0))  
    
    #2.Relationships at age
    
    #Age-TL
    fn.plt.rel(age.F,mid.TL.F,add="YES",age.F,mid.TL.M,"pink",4)
    legend("bottomright",c("F","M"),lty=1,lwd=2,col=c("pink","blue"),bty='n',cex=1.25)
    fn.mtxt("TL (cm)")
    
    #Age-TwT
    fn.plt.rel(age.F,mid.TwT.F,add="YES",age.F,mid.TwT.M,"pink",4)
    fn.mtxt("Twt (kg)")
    
    #Age- female maturity
    fn.plt.rel(age.F,Mat.age,add="YES",age.F,Mat.k,"grey20","grey70")
    legend("bottomright",c("at age","knife-edge"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)
    legend("topleft","Female",bty='n',cex=1.25)
    fn.mtxt("Prop. mat.")
    
    #Age-litter size * proportion mature
    fn.plt.rel(age.F,Litter.size.age,add="YES",age.F,Litter.size.k,"grey20","grey70")
    legend("bottomright",c("at age","knife-edge"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)
    fn.mtxt("Litter size x prop. mat.")
    
    #Age-female mortality
    fn.plt.rel(age.F,Mort.F.age,add="YES",age.F,Mort.F.k,"white",3)
    fn.mtxt(expression("M  " (year^-1)))
    
    
    #Age-selectivity 
    fn.plt.rel(age.F,Select.age.F,add="YES",age.F,Select.age.M,"pink",4)
    mtext("Age",1,2)
    fn.mtxt("Rel. selectivity")
    
    dev.off()
    
  }
  
  
  #Display relations separately
  if(BaseCase=="Size-based")
  {
    fn.plt.rel1=function(one,two,add,three,four,CL,CL1)
    {
      YLIM=c(0,max(two)*1.15)
      if(add=="YES") YLIM=c(0,max(c(two,four)))
      plot(one,two,xlab="",ylab="",type='l',lwd=2,col=CL,ylim=YLIM,cex.axis=1.3)
      if(add=="YES") lines(three,four,lwd=2,col=CL1)
    }
    txt.off=2.65
    
    #1. At length relationships
    MaxLeN=round(max(Age.growth$TL))
    ID.MaxLeN=1:which.min(abs(TL.bin.mid-MaxLeN))
    
    fn.fig("Visualise data/Input relations_at_length",1800, 2400)  
    par(mfcol=c(3,1),mai=c(.3,.4,.1,.4),oma=c(3,2.5,.25,.1),las=1,mgp=c(2,.6,0))
    #par(mfcol=c(2,2),mai=c(.3,.4,.1,.4),oma=c(3,2.5,.25,.1),las=1,mgp=c(2,.6,0))
    
    #TL-TwT 
    fn.plt.rel1(TL.bin.mid[ID.MaxLeN],TwT.bin.fem[ID.MaxLeN],add="YES",TL.bin.mid[ID.MaxLeN],TwT.bin.mal[ID.MaxLeN],"pink",4)
    mtext("Total weight (kg)",2,txt.off,las=3,cex=1.5)
    legend("bottomright",c("F","M"),lty=1,lwd=2,col=c("pink","blue"),bty='n',cex=1.75)
    
    #TL- female maturity
    ADD="NO"
    if(length(which(Tabla.scen$Model_type=="Length-based" & Tabla.scen$Maturity=="knife edge"))>0) ADD="YES"
    fn.plt.rel1(TL.bin.mid[ID.MaxLeN],Mat.bin[ID.MaxLeN],add=ADD,TL.bin.mid,ifelse(TL.bin.mid>=Maturity.at.size[1],1,0),"grey20","grey70")
    mtext("Proportion mature",2,txt.off,las=3,cex=1.5)
    legend("topleft","Female",bty='n',cex=1.75)
    if(ADD=="YES")legend("left",c("at length","knife-edge"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)
    #mtext("Total length (cm)",1,2,cex=1.5)
    
    #TL-mortality
    ADD="NO"
    if(length(which(Tabla.scen$Model_type=="Length-based" & Tabla.scen$M=="at length"))>0) ADD="YES"
    fn.plt.rel1(TL.bin.mid[ID.MaxLeN],rep(unique(Mort.F.k),length(TL.bin.mid[ID.MaxLeN])),add=ADD,TL.bin.mid[ID.MaxLeN],Mort.bin[ID.MaxLeN],"grey20","grey70")
    mtext(expression("Natural mortality " (year^-1)),2,txt.off,las=3,cex=1.5)
    if(ADD=="YES")legend("bottomleft",c("constant","at length"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)
    
    #TL-selectivity 
    # fn.plt.rel1(TL.bin.mid,Select.bin,add="NO",Dum,Dum,1,4)  
    # mtext("Relative selectivity",2,txt.off,las=3,cex=1.5)
     mtext("Total length (cm)",1,3,cex=1.5)
     
    dev.off()
    

    #Varying selectivities
    #ClS=grey.colors(ncol(Total.sel.all),start=0.1,end=0.80)
    colfunc <- colorRampPalette(c("red","yellow","springgreen","royalblue"))
    Id.Sel=match(Yrs.sel,colnames(Total.sel.all))
    ClS=colfunc(ncol(Total.sel.all[,Id.Sel])) 
    
    #All
    fn.fig("Visualise data/Input relations_at_length_varying_selectivity_all",2400, 2400) 
    par(mai=c(.8,.8,.1,.01),las=1,xpd=T,mgp=c(.65,.8,0))
    plot(TL.bin.mid,Total.sel.all[,1],col=1,pch=19,ylab="",xlab="",cex.axis=1.25)
    for(i in 2:ncol(Total.sel.all[,Id.Sel])) lines(TL.bin.mid, Total.sel.all[,i],col=ClS[i],lwd=2)
    mtext("Relative selectivity",2,las=3,cex=1.75,line=2.6)
    mtext("Total length (cm)",1,cex=1.75,line=-1.5,outer=T)
    legend("topleft",colnames(Total.sel.all[,Id.Sel]),bty='n',lwd=2,lty=1,col=ClS,cex=.8)
    dev.off()
    
    #By zone
    rnd=round(length(ClS)/3)
    SPlt=list(1:rnd,(rnd+1):(rnd+rnd),(rnd+rnd+1):length(ClS))
    
    fn.fig("Visualise data/Input relations_at_length_varying_selectivity_zones",2000, 2400)
    par(mfcol=c(3,1),mai=c(.1,.5,.15,.01),oma=c(4,2,1.5,.1),las=1,xpd=T,mgp=c(.65,0.8,0))
    plot(TL.bin.mid,Total.sel.West[,1],col=1,pch=19,ylab="",xlab="",cex.axis=1.25)
    for(i in 2:ncol(Total.sel.all[,Id.Sel])) lines(TL.bin.mid, Total.sel.West[,i],col=ClS[i],lwd=2)
    legend("topright",legend ="West" ,bty='n',cex=1.75)
    legend("topleft",colnames(Total.sel.all[,Id.Sel])[SPlt[[1]]],bty='n',cex=1.1,lwd=2,lty=1,col=ClS[SPlt[[1]]])
    
    plot(TL.bin.mid,Total.sel.Zn1[,1],col=1,pch=19,ylab="",xlab="",cex.axis=1.25)
    for(i in 2:ncol(Total.sel.all[,Id.Sel])) lines(TL.bin.mid, Total.sel.Zn1[,i],col=ClS[i],lwd=2)
    legend("topright",legend ="Zone1" ,bty='n',cex=1.75)
    legend("topleft",colnames(Total.sel.all[,Id.Sel])[SPlt[[2]]],bty='n',cex=1.1,lwd=2,lty=1,col=ClS[SPlt[[2]]])
    
    plot(TL.bin.mid,Total.sel.Zn2[,1],col=1,pch=19,ylab="",xlab="",cex.axis=1.25)
    for(i in 2:ncol(Total.sel.all[,Id.Sel])) lines(TL.bin.mid, Total.sel.Zn2[,i],col=ClS[i],lwd=2)
    legend("topright",legend ="Zone2" ,bty='n',cex=1.75)
    legend("topleft",colnames(Total.sel.all[,Id.Sel])[SPlt[[3]]],bty='n',cex=1.1,lwd=2,lty=1,col=ClS[SPlt[[3]]])
    
    mtext("Relative selectivity",2,las=3,cex=1.75,line=-.85,outer=T)
    mtext("Total length (cm)",1,cex=1.75,line=2,outer=T)
    
    dev.off()
    
    
    
    #2.Relationships at age
    fn.fig("Visualise data/Input relations_at_age",2400, 2400)  
    par(mfcol=c(3,2),mai=c(.3,.4,.1,.4),oma=c(3,2.5,.25,.1),las=1,mgp=c(2.2,.6,0))
    
    #Age-TL
    fn.plt.rel1(age.F,mid.TL.F,add="YES",age.F,mid.TL.M,"pink",4)
    mtext("Total length (cm)",2,txt.off,las=3,cex=1.5)
    legend("bottomright",c("F","M"),lty=1,lwd=2,col=c("pink","blue"),bty='n',cex=1.35)
    
    #Age-TwT
    fn.plt.rel1(age.F,mid.TwT.F,add="YES",age.F,mid.TwT.M,"pink",4)
    mtext("Total weight (kg)",2,txt.off,las=3,cex=1.5)
    
    #Age- female maturity  
    plot(age.F,Mat.k,pch=19,col="grey20",ylim=c(0,max(Mat.k)*1.15),cex.axis=1.3,ylab="",xlab="")
    legend("topleft","Female",bty='n',cex=1.35)
    mtext("Proportion mature",2,txt.off,las=3,cex=1.5)
    mtext("Age",1,2,cex=1.5)  
    
    #Litter size
    plot(age.F,FecU*Mat.k,pch=19,col="grey20",ylim=c(0,max(FecU*Mat.k)*1.15),cex.axis=1.3,ylab="",xlab="")
    legend("topleft","Female",bty='n',cex=1.35)
    mtext("Litter size",2,txt.off,las=3,cex=1.5) 
    
    #Age-female mortality
    ADD="NO"
    if(length(which(Tabla.scen$Model_type=="Age-structured" & Tabla.scen$M=="at age"))>0) ADD="YES"
    fn.plt.rel1(age.F,Mort.F.k,add=ADD,age.F,Mort.F.age,"grey20","grey70")   
    if(ADD=="YES")legend("bottomright",c("constant","at age"),lty=1,lwd=2,col=c("grey20","grey70"),bty='n',cex=1.25)  
    mtext(expression("Natural mortality " (year^-1)),2,txt.off,las=3,cex=1.5)
    
    #Age-selectivity 
    fn.plt.rel1(age.F,Select.age.F,add="YES",age.F,Select.age.M,"pink",4)
    mtext("Relative selectivity",2,txt.off,las=3,cex=1.5)
    mtext("Age",1,2,cex=1.5)  
    
    dev.off()
  }
  
  Zns=c("Joint","North","Closed","West","Closed.metro","Zone1","Zone2")
  Zns.leg=c("JANSF","WANCSF","Ningaloo","WCDGDLF","Metro closure","Zone1","Zone2")
  COL.prop=c("aquamarine2","lightseagreen","seagreen4","lightgreen","olivedrab4","olivedrab3","mediumseagreen")
  names(COL.prop)=Zns
  
  #Conventional tagging data 
  if(Move.mode=="Population-based")      
  {
    n.zns=unique(Conv.tg.rel_plot$Rel.zone)
    
    col.zn=COL.prop[match(Zns,names(COL.prop))]
    
    fn.fig("Visualise data/Conventional tagging",2400, 2400)  
    par(mfcol=c(3,1),mai=c(.4,.5,.05,.1),oma=c(2,2,.01,.1),las=1)  
    #releases
    plot(1:10,ylim=c(0,max(Conv.tg.rel_plot$Number)+1),xlim=c(Conv.tg.StartYear,Conv.tg.EndYear),ylab="",xlab="",cex.axis=1.5)
    for(z in 1:length(n.zns))
    {
      a=subset(Conv.tg.rel_plot,Rel.zone==n.zns[[z]])
      points(a$Yr.rel,a$Number,pch=19,col=col.zn[z],cex=1.75)
    }
    legend("topright",c("WC", "ZN1", "ZN2"),bty='n',pch=c(19,19,19),cex=2.5,col=col.zn)
    legend("top","Releases",bty='n',cex=2.5)
    
    #recaptures
    plot(1:10,ylim=c(0,max(Conv.tg.rec_plot$Number)+1),xlim=c(Conv.tg.StartYear,Conv.tg.EndYear),ylab="",xlab="",cex.axis=1.5)
    for(z in 1:length(n.zns))
    {
      a=subset(Conv.tg.rec_plot,Rec.zone==n.zns[[z]])
      points(a$Yr.rec,a$Number,pch=19,col=col.zn[z],cex=1.75)
    }
    legend("top","Recaptures",bty='n',cex=2.5)
    
    #reporting                           
    plot(1:10,ylim=c(0,1),xlim=c(Conv.tg.StartYear,Conv.tg.EndYear),ylab="",xlab="",cex.axis=1.5)
    Reporting_plot$Pre.ZN1=with(Reporting_plot,ifelse(Zn1==ZN1.pred,"YES","NO"))
    Reporting_plot$Pre.ZN2=with(Reporting_plot,ifelse(Zn2==ZN2.pred,"YES","NO"))
    Reporting_plot$Pre.WC=with(Reporting_plot,ifelse(WC==WC.pred,"YES","NO"))
    
    with(Reporting_plot,points(YR,WC,col=col.zn[1],cex=1.5))
    with(Reporting_plot,points(YR,Zn1,col=col.zn[2],cex=1.5))
    with(Reporting_plot,points(YR,Zn2,col=col.zn[3],cex=1.5))
    
    Reporting_plot.ob=subset(Reporting_plot,Pre.ZN1=="NO")
    with(Reporting_plot.ob,points(YR,WC,col=col.zn[1],pch=19,cex=1.75))
    with(Reporting_plot.ob,points(YR,jitter(Zn1,5),col=col.zn[2],pch=19,cex=1.75))
    with(Reporting_plot.ob,points(YR,Zn2,col=col.zn[3],pch=19,cex=1.75))
    legend("bottomleft","Reporting rate",bty='n',cex=2.5)
    
    mtext("Year",1,0.5,outer=T,cex=2.75)
    mtext("                        Numbers",2,-0.75,outer=T,las=3,cex=2.75)
    mtext("Proportion",2,3,las=3,cex=2.75)
    
    dev.off()
    
  }
  
  if(Move.mode=="Individual-based")
  {
    Zns=unique(Conv.tg.rec.exp$Rec.zn)
    col.zn=COL.prop[match(sort(Zns),names(COL.prop))]
    
    names(col.zn)=AREAS
    jit.zn=c(.9,1,1.1)
    names(jit.zn)=names(col.zn)
    
    fn.fig("Visualise data/Conventional tagging_recaptures only",2000, 2400)  
     par(mfcol=c(3,1),mai=c(.4,.5,.3,.2),oma=c(2,2,.01,.1),mgp=c(2,.65,0),las=1)  
    
    #recaptures
     dummy=Conv.tg.rec.exp
     dummy$N=1
    Ymax=max(aggregate(N~Rec.zn+DaysAtLarge+Rel.zn,dummy,sum)$N)
    fn.rec.ind.base=function(where)
    {
      dat=subset(Conv.tg.rec.exp,Rel.zn==where)
      dat$N=1
      dat=aggregate(N~Rec.zn+DaysAtLarge,dat,sum)
      dummy=data.frame(CL=col.zn,Rec.zn=names(col.zn),jitr=jit.zn)
      dat=merge(dat,dummy,by="Rec.zn")
      dat$CL=as.character(dat$CL)
      with(dat,plot(DaysAtLarge,N*jitr,ylim=c(0,Ymax),xlim=c(0,max(Conv.tg.rec.exp$DaysAtLarge)),
                    ylab="",cex=2,xlab="",cex.axis=2,bg=CL,pch=21))
      mtext(paste("released in",where,sep=" "),3,0.5,cex=1.5)
    }
    fn.rec.ind.base("West")
    legend("right",c("recaptured in West", "recaptured in Zone1", "recaptured in Zone2"),
           bty='n',pch=c(19,19,19),cex=2,col=col.zn)
    fn.rec.ind.base("Zone1")
    fn.rec.ind.base("Zone2")  
    mtext("Days at liberty",1,0.5,outer=T,cex=2.25)
    mtext("Numbers",2,-0.75,outer=T,las=3,cex=2.25)  
    
    dev.off()
    
  }
  
  #Acoustic tagging data    
  if(Move.mode=="Individual-based")
  {
    Zns=unique(Acous.tg.rec$Rec.zn)
    col.zn=COL.prop[match(sort(Zns),names(COL.prop))]
    
    names(col.zn)=AREAS
    jit.zn=c(.9,1,1.1)
    names(jit.zn)=names(col.zn)
    
    
    fn.fig("Visualise data/Acoustic tagging_recaptures only",2000, 2400)  
     par(mfcol=c(3,1),mai=c(.4,.5,.3,.2),oma=c(2,2,.01,.1),mgp=c(2,.65,0),las=1)  
    
    #recaptures
     dummy=Acous.tg.rec
     dummy$N=1
     Ymax=max(aggregate(N~Rec.zn+DaysAtLarge+Rel.zn,dummy,sum)$N)
    fn.rec.ind.base=function(where)
    {
      dat=subset(Acous.tg.rec,Rel.zn==where)
      YRS=1:3
      if(nrow(dat)==0)
      {
        plot(YRS[1]:YRS[length(YRS)],YRS[1]:YRS[length(YRS)],ylim=c(0,1),
             ylab="",xlab="",cex.axis=2,col='transparent',xaxt='n',yaxt='n')        
      }
      if(nrow(dat)>0)
      {
        dat$N=1
        dat=aggregate(N~Rec.zn+DaysAtLarge,dat,sum)
        dummy=data.frame(CL=col.zn,Rec.zn=names(col.zn),jitr=jit.zn)
        dat=merge(dat,dummy,by="Rec.zn")
        dat$CL=as.character(dat$CL)
        with(dat,plot(DaysAtLarge,N*jitr,ylim=c(0,Ymax),xlim=c(0,max(Acous.tg.rec$DaysAtLarge)),
                ylab="",xlab="",cex=2,cex.axis=2,bg=CL,pch=21))
       }
      mtext(paste("released in",where,sep=" "),3,0.5,cex=1.25)
    }
    
    fn.rec.ind.base("West")
    legend("right",c("detected in West", "detected in Zone1", "detected in Zone2"),
           bty='n',pch=c(19,19,19),cex=2,col=col.zn)
    fn.rec.ind.base("Zone1")
    fn.rec.ind.base("Zone2")
    mtext("Days at liberty",1,0.5,outer=T,cex=2.25)
    mtext("Number of events",2,-0.75,outer=T,las=3,cex=2.25)
    dev.off()
    
    #Multiple recaptures by tagID
    # fn.fig("Acoustic tagging_recaptures only_multirecap",2000, 2400)  
    #  par(mfcol=c(1,1),mai=c(.4,.5,.3,.2),oma=c(2,2,.01,.1),las=1,mgp=c(2,0.7,0))  
    # barplot(table(Acous.tg.rec$TagID))
    # box()
    # mtext("Tag ID",1,0.75,outer=T,cex=2.25)
    # mtext("Number of detections",2,-0.75,outer=T,las=3,cex=2.25)
    # dev.off()
  }
  
  
  #C.3 Create input files (data and parameters) for population models as required by ADMB
  
    #C.3.1. Input data and parameters
  setPath=function(Scen)setwd(paste(hndl,AssessYr,"/",Scen,sep="")) #set paths
  
  Scenarios=Tabla.scen
  Scenarios <- data.frame(lapply(Scenarios, as.character), stringsAsFactors=FALSE) #all strings
  
  Inputs=vector('list',nrow(Scenarios))
  names(Inputs)=Scenarios$Model
  
  get.yr=function(dat) as.numeric(substr(dat,start=1,stop=4)) 
  Yr=sort(get.yr(Ktch.All.1975$FINYEAR))
  yr.start=Yr[1]                      #first year
  yr.end=Yr[length(Yr)]        	    #last year  
  
  
  #make sure data ordered by year
  Ktch.All.1975=Ktch.All.1975[order(Ktch.All.1975$FINYEAR),]
  Ktch.All.West.1975=Ktch.All.West.1975[order(Ktch.All.West.1975$FINYEAR),]
  Ktch.All.zn1.1975=Ktch.All.zn1.1975[order(Ktch.All.zn1.1975$FINYEAR),]
  Ktch.All.zn2.1975 =Ktch.All.zn2.1975[order(Ktch.All.zn2.1975$FINYEAR),]
  
    #Monthly
  Cpue.all=Cpue.all[order(Cpue.all$Finyear),]
  if(exists('Cpue.folly'))Cpue.folly=Cpue.folly[order(Cpue.folly$FINYEAR),]
  if(nrow(Cpue.West)>0) Cpue.West=Cpue.West[order(Cpue.West$Finyear),]
  if(nrow(Cpue.zn1)>0)  Cpue.zn1=Cpue.zn1[order(Cpue.zn1$Finyear),]
  if(nrow(Cpue.zn2)>0)  Cpue.zn2=Cpue.zn2[order(Cpue.zn2$Finyear),]

    #Daily
  Cpue.all.daily=Cpue.all.daily[order(Cpue.all.daily$Finyear),]
  if(nrow(Cpue.West.daily)>0) Cpue.West.daily=Cpue.West.daily[order(Cpue.West.daily$Finyear),]
  if(nrow(Cpue.zn1.daily)>0)  Cpue.zn1.daily=Cpue.zn1.daily[order(Cpue.zn1.daily$Finyear),]
  if(nrow(Cpue.zn2.daily)>0)  Cpue.zn2.daily=Cpue.zn2.daily[order(Cpue.zn2.daily$Finyear),]
  
  Eff.tdgdlf=Eff.tdgdlf[order(Eff.tdgdlf$FINYEAR),]
  Eff.wst=Eff.wst[order(Eff.wst$FINYEAR),]
  Eff.zn1=Eff.zn1[order(Eff.zn1$FINYEAR),]
  Eff.zn2=Eff.zn2[order(Eff.zn2$FINYEAR),]
  
  
  size.all=size.all[order(size.all$FINYEAR),]
  size.wst=size.wst[order(size.wst$FINYEAR),]
  size.zn1=size.zn1[order(size.zn1$FINYEAR),]
  size.zn2=size.zn2[order(size.zn2$FINYEAR),]
  
  size.all_7=size.all_7[order(size.all_7$FINYEAR),]
  size.wst_7=size.wst_7[order(size.wst_7$FINYEAR),]
  size.zn1_7=size.zn1_7[order(size.zn1_7$FINYEAR),]
  size.zn2_7=size.zn2_7[order(size.zn2_7$FINYEAR),]
  
  avg.wt=avg.wt%>%rename(finyear=Finyear)%>%
                  arrange(finyear)
  
  if(nrow(avg.wt.wst)>0) avg.wt.wst=avg.wt.wst%>%rename(finyear=Finyear)%>%
                                                   arrange(finyear)
  if(nrow(avg.wt.zn1)>0) avg.wt.zn1=avg.wt.zn1%>%rename(finyear=Finyear)%>%
                                                   arrange(finyear)
  if(nrow(avg.wt.zn2)>0) avg.wt.zn2=avg.wt.zn2%>%rename(finyear=Finyear)%>%
                                                   arrange(finyear)


  if(Acoust.format=="SS3")
  {
    Acous.tg.rel=Acous.tg.rel[order(Acous.tg.rel$FinYear.rel),]
    Acous.tg.rec=Acous.tg.rec[order(Acous.tg.rec$FinYear.rec),]  
  }
  
  #Add future projections
  #note: future catch is set at the average catch of last 5 years
  Ktch.future=function(KTCH) rep(mean(KTCH[(length(KTCH)-4):length(KTCH)]),Yrs.future)
  
  #Prior for F 
  fn.fig("Visualise data/Prior_F.init",2000, 2000)    
  par(las=1,mai=c(1,1.15,.1,.15),mgp=c(3.5,.75,0))
  plot(density(rlnorm(10000, meanlog = log(Prior.mean.Fo), sdlog = Prior.SD.Log.Fo)), lwd=3,
       main="",xlab=expression('F'['init']),cex.lab=2,cex.axis=1.25,col=1)
  dev.off()
  
  #Prior for r (biomass dynamics) using Leslie matrix
  Fecu=unlist(Fecundity)
  Fecu[1]=mean(Fecu)
  fn.source("Leslie.matrix.R")
  Rprior=fun.Leslie(N.sims=1000,k=Growth.F$k,Linf=Growth.F$FL_inf,Aver.T=18,A=Max.age.F,first.age=0,
                    RangeMat=unlist(Age.50.mat),Rangefec=Fecu,sexratio=pup.sx.ratio,
                    Reprod_cycle=1/Breed.freq.min,Hoenig.only="YES")   
  
  #get mean and SD from lognormal distribution
  LogN.pars=fitdistr(Rprior, "lognormal")  
  log_mean.r=LogN.pars$estimate[1]    #already in log space     
  log_sd.r=LogN.pars$estimate[2]      #already in log space     
  
  write.csv(data.frame(log_mean.r=log_mean.r,log_sd.r=log_sd.r),"log_r_priors.csv",row.names=F)
  
  fn.fig("Visualise data/Prior_r", 2000, 2000)
   par(las=1,mai=c(1,1.15,.1,.15),mgp=c(3.5,.75,0))
   plot(density(rlnorm(10000, meanlog = log_mean.r, sdlog = log_sd.r)),
       lwd=3,main="",xlab='Intrinsic rate of increase',cex.lab=2,cex.axis=1.25,col=1)
  dev.off()
  
  
  #Movement transition matrix from acoustic tagging data
  if(Sim.trans.Mat=="NO")
  {
    if(species=="WH")MOV.TRANS.MAT=matrix(c(0.95,0.05,0,0.01,.95,0.04,0,0.28,0.72),nrow=3,byrow=T)
    if(species=="GM")MOV.TRANS.MAT=matrix(c(0.44,0.56,0,0.17,0.6,0.23,0,0.04,0.96),nrow=3,byrow=T)  
    rownames(MOV.TRANS.MAT)=paste("from",c("WC","ZN1","ZN2"))
    colnames(MOV.TRANS.MAT)=paste("to",c("WC","ZN1","ZN2"))  
    write.csv(MOV.TRANS.MAT,"MOV.TRANS.MAT.csv",row.names=T)  
    MOV.TRANS.MAT.no.movement=matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,byrow=T)
  }
  
  
  #Combine monthly and daily cpue into single file
  Cpue.all.monthly=Cpue.all
  Cpue.all.monthly.hours=Cpue.all.hours
  Cpue.West.monthly=Cpue.West
  Cpue.zn1.monthly=Cpue.zn1
  Cpue.zn2.monthly=Cpue.zn2

  Cpue.all=rbind(Cpue.all.monthly,Cpue.all.daily)
  Cpue.all.hours=rbind(Cpue.all.monthly.hours,Cpue.all.daily.hours)
  Cpue.West=rbind(Cpue.West.monthly,Cpue.West.daily)
  Cpue.zn1=rbind(Cpue.zn1.monthly,Cpue.zn1.daily)
  Cpue.zn2=rbind(Cpue.zn2.monthly,Cpue.zn2.daily)

  
  #Show cpues used 
  id.Yrs=which(Yr%in%as.numeric(substr(Cpue.all$Finyear,1,4)))
  fn.fig("Visualise data/Input_CPUE",2000, 2000)   
  par(mfcol=c(2,1),mai=c(.75,.75,.05,.175),las=1,mgp=c(2,.6,0),cex.axis=1.25,cex.lab=1.75)
  
  plot(Yr[id.Yrs],Cpue.all[,match('Mean',names(Cpue.all))],ylab="",xlab="")
  points(Yr[which(Yr%in%as.numeric(substr(Cpue.all.monthly$Finyear,1,4)))],Cpue.all.monthly$Mean,pch=19)
  legend("topright",c("Standardised"),text.col=1:2,cex=1.5,bty='n')
  
  plot(Yr[id.Yrs],Cpue.all[,match('Mean',names(Cpue.all))],ylab="",xlab="",col="transparent",ylim=c(0,max(c(Cpue.zn1$Mean,Cpue.zn2$Mean,Cpue.West$Mean))))
  if(nrow(Cpue.West)>0)   
  {
    clrS=COL.prop[match("West",names(COL.prop))]
    IIDi=which(Yr%in%as.numeric(substr(Cpue.West$Finyear,1,4)))
    if(nrow(Cpue.West)>0)points(Yr[IIDi],Cpue.West[,match('Mean',names(Cpue.West))],ylab="",xlab="",xlim=c(min(Yr),max(Yr)),pch=21,bg=clrS)
    if(nrow(Cpue.West.monthly)>0)points(Yr[IIDi][1:nrow(Cpue.West.monthly)],Cpue.West.monthly$Mean,pch=21,bg=clrS)    
  }
  if(nrow(Cpue.zn1)>0)
  {
    clrS=COL.prop[match("Zone1",names(COL.prop))]
    IIDi=which(Yr%in%as.numeric(substr(Cpue.zn1$Finyear,1,4)))
    if(nrow(Cpue.zn1)>0)points(Yr[IIDi],Cpue.zn1[,match('Mean',names(Cpue.zn1))],ylab="",xlab="",pch=21,bg=clrS)
    if(nrow(Cpue.West.monthly)>0)points(Yr[IIDi][1:nrow(Cpue.zn1.monthly)],Cpue.zn1.monthly$Mean,pch=21,bg=clrS)    
  }
  if(nrow(Cpue.zn2)>0)
  {
    clrS=COL.prop[match("Zone2",names(COL.prop))]
    IIDi=which(Yr%in%as.numeric(substr(Cpue.zn2$Finyear,1,4)))
    if(nrow(Cpue.zn2)>0)points(Yr[IIDi],Cpue.zn2[,match('Mean',names(Cpue.zn2))],ylab="",xlab="",pch=21,bg=clrS)
    if(nrow(Cpue.West.monthly)>0)points(Yr[IIDi][1:nrow(Cpue.zn2.monthly)],Cpue.zn2.monthly$Mean,pch=21,bg=clrS)    
  }
  legend("topright",c("Standard. WC","Standard. ZN1","Standard. ZN2"),cex=1.5,bty='n',
         text.col=COL.prop[match(c("West","Zone1","Zone2"),names(COL.prop))])
  mtext("Financial year",1,outer=T,line=-1.5,cex=2)
  mtext("Relative cpue",2,outer=T,line=-2.0,las=3,cex=2) 
  dev.off()
  
  #convert character zone to number #ACA
  if(Move.mode=="Individual-based")
  {
    Conv.tg.rec.exp$Rel.zn=with(Conv.tg.rec.exp,
            ifelse(Rel.zn=="West",1,ifelse(Rel.zn=="Zone1",2,
            ifelse(Rel.zn=="Zone2",3,NA))))
    Conv.tg.rec.exp$Rec.zn=with(Conv.tg.rec.exp,
            ifelse(Rec.zn=="West",1,ifelse(Rec.zn=="Zone1",2,
            ifelse(Rec.zn=="Zone2",3,NA))))
    
    Acous.tg.rec$Rel.zn=with(Acous.tg.rec,
            ifelse(Rel.zn=="West",1,ifelse(Rel.zn=="Zone1",2,
            ifelse(Rel.zn=="Zone2",3,NA))))
    Acous.tg.rec$Rec.zn=with(Acous.tg.rec,
            ifelse(Rec.zn=="West",1,ifelse(Rec.zn=="Zone1",2,
            ifelse(Rec.zn=="Zone2",3,NA))))
  }
  
  
  #Get bin of minimum size moving     
  Smallest_size_tagged=TL.bin.up[match(round(Smallest_size_tagged/10)*10,round(TL.bin.mid/10)*10)]
  Smallest_size_tagged.index=match(Smallest_size_tagged, TL.bin.low)
  
  #Export conventional and acoustic tagging data used for movement rates
  a=getwd()
      #data for stand-alone movement rate estimation
  if(Acoust.format=="Individual-based")
  {
    Movement.rate.data=list()
    Movement.rate.data$Conv.tg.nzones=Conv.tg.nzones
    Movement.rate.data$Conv.tg.recaptures=Conv.tg.recs
    Movement.rate.data$Conv.tg.rec=Conv.tg.rec.exp 
    Movement.rate.data$Acous.tg.detections=Acous.detections
    Movement.rate.data$Acous.tg.tags=Acous.ntags
    Movement.rate.data$Acous.tg.obs.per.tags=Acous.TagID_observations
    Movement.rate.data$Acous.tg.start=Acous.TagID_start
    Movement.rate.data$Acous.tg.end=Acous.TagID_end
    Movement.rate.data$Acous.tg.data=Acous.tg.rec
    
    #Create .dat file
    n=Movement.rate.data
    setwd(paste("C:/Matias/Analyses/Movement rate estimation/Joint.estim_ind.base.mod/",Spec,sep=""))
    FILE="model.dat"
    Hdr="#Model inputs"
    write(Hdr,file = FILE)
    for(k in 1:length(n))
    {
      nn=n[[k]]
      if(is.data.frame(nn)|is.matrix(nn))
      {
        Hdr=paste("#",paste(c(names(n)[k],"(",names(nn),")"),collapse=' '))
        write(Hdr,file = FILE,append=T)      
        write.table(nn,file = FILE,row.names=F,col.names=F,append=T)
      }else
      {
        Hdr=paste("#",names(n)[k],sep='')
        write(Hdr,file = FILE,append=T)
        write(n[[k]],file = FILE,sep = "\t",append=T)
      }
    }
    write("#eof",file = FILE,append=T)
    write(999,file = FILE,sep = "\t",append=T)
    rm(n)
  }
  
      #data for displaying inputs
  hndl.mov.rte="C:/Matias/Analyses/Movement rate estimation/Joint.estim_ind.base.mod/Show Gummy and whiskery outputs/"
  write.csv(Conv.tg.rec.exp,paste(hndl.mov.rte,"Conv.",Spec,".csv",sep=""),row.names=F)
  write.csv(Acous.tg.rec,paste(hndl.mov.rte,"Acous.",Spec,".csv",sep=""),row.names=F)
  setwd(a)
  
  #Get last years of monthly and daily cpue
  Yr_q_daily=as.numeric(substr(Cpue.all.daily$Finyear[1],1,4))

  #create all .dat files        
    #fn for adding missing size classes and years
  fn.add.missing.size.year=function(dummy)
  {
    drop.yr=rowSums(dummy[,2:ncol(dummy)])
    names(drop.yr)=dummy$FINYEAR
    drop.yr=subset(drop.yr,drop.yr<Min.obs)
    dummy=subset(dummy,!FINYEAR%in%names(drop.yr))
    Yrs.dat=dummy$FINYEAR
    All.yrs=Ktch.All.1975$FINYEAR
    misin.yr=sort(All.yrs[which(!All.yrs%in%Yrs.dat)])
    
    #add missing size
    dummy=dummy[,-match('FINYEAR',names(dummy))]
    id=TL.bin.low[which(!TL.bin.low%in%as.numeric(colnames(dummy)))]
    AdD=matrix(0,nrow=nrow(dummy),ncol=length(id))
    colnames(AdD)=id
    dummy=cbind(dummy,AdD)
    dummy=dummy[,match(TL.bin.low,colnames(dummy))]
    colnames(dummy)=TL.bin.mid
    rownames(dummy)=Yrs.dat
    
    #add missing years
    Msin.yrs=matrix(0,nrow=length(misin.yr),ncol=ncol(dummy))
    colnames(Msin.yrs)=TL.bin.mid
    rownames(Msin.yrs)=misin.yr
    dummy=rbind(dummy,Msin.yrs)
    dummy=dummy[order(rownames(dummy)),]
    return(dummy)
  }
  
    #fn for adding missing year
  fn.add.missing.year=function(dummy)
  {
    Yrs.dat=dummy$finyear
    All.yrs=Ktch.All.1975$FINYEAR
    misin.yr=sort(All.yrs[which(!All.yrs%in%Yrs.dat)])
    dummyCV=dummy$CV
    dummy=dummy$Mean.wgt
    names(dummy)=names(dummyCV)=Yrs.dat
    Msin.yrs=rep(0,length(misin.yr))
    names(Msin.yrs)=misin.yr
    dummy=c(dummy,Msin.yrs)
    dummy=dummy[order(names(dummy))]
    dummyCV=c(dummyCV,Msin.yrs)
    dummyCV=dummyCV[order(names(dummyCV))]
    return(list(Mean=dummy,CV=dummyCV))
  }
  
    #fn for calculating proportions
  fn.size.prop=function(dat)
  {
    dat=dat/rowSums(dat)
    dat[is.na(dat)]=0
    
    #add small value for 0 observations within range of data if Dirichlet is used
    if(Size_like=="Dirichlet")
    {
      drop.yr=rowSums(dat)
      indx=match(names(which(drop.yr>0)),names(drop.yr))
      dat[indx,][dat[indx,]==0]=Dirichlet.small.value
      dat=dat/rowSums(dat)  #rescale
      dat[is.na(dat)]=0
    }
    return(dat)
  }
  
    #fn for filling in missing cpue years
  fn.add.missing.cpue.year=function(dummy)
  {
    Yrs.dat=dummy$Finyear
    All.yrs=Ktch.All.1975$FINYEAR
    misin.yr=sort(All.yrs[which(!All.yrs%in%Yrs.dat)])
    if(length(misin.yr)>0)
    {
      ADD=dummy[1:length(misin.yr),]
      ADD[,]=-100
      ADD$Finyear=misin.yr
      dummy=rbind(dummy,ADD)
      dummy=dummy[order(dummy$Finyear),]
    }
    return(dummy)
  }
  Cpue.all=fn.add.missing.cpue.year(Cpue.all)
  Cpue.all.hours=fn.add.missing.cpue.year(Cpue.all.hours)
  if(nrow(Cpue.West)>0)Cpue.West=fn.add.missing.cpue.year(Cpue.West)
  if(nrow(Cpue.zn1)>0)Cpue.zn1=fn.add.missing.cpue.year(Cpue.zn1)
  if(nrow(Cpue.zn2)>0)Cpue.zn2=fn.add.missing.cpue.year(Cpue.zn2)  
  
  fn.add.missing.cpue.year.monthly=function(dummy)
  {
    Yrs.dat=dummy$Finyear
    All.yrs=Ktch.All.1975$FINYEAR[1:31]
    misin.yr=sort(All.yrs[which(!All.yrs%in%Yrs.dat)])
    if(length(misin.yr)>0)
    {
      ADD=dummy[1:length(misin.yr),]
      ADD[,]=-100
      ADD$Finyear=misin.yr
      dummy=rbind(dummy,ADD)
      dummy=dummy[order(dummy$Finyear),]
    }
    return(dummy)
  }
  Cpue.all.monthly=fn.add.missing.cpue.year.monthly(Cpue.all.monthly)
  
   
  
    #Loop over all scenarios  
  for(i in 1:nrow(Scenarios))
  {
    d=Scenarios[i,]
    La.lista=list(yrs_ktch=NA,nage_F=NA,nage_M=NA,nTL=NA,
                  fec=NA,size_mat=NA,sx_ratio=NA,breed=NA,age_mat=NA,
                  STM.F=NA,STM.M=NA,mean.size.birth=NA,SD.mean.size.birth=NA,
                  TL.bin.mid=NA,TL.bin.low=NA,TL.bin.up=NA,
                  Max.age.size=NA,TwT.bin.F=NA,TwT.bin.M=NA,
                  M=NA,Selectivity=NA,
                  Selectivity.west=NA,Selectivity.zn1=NA,Selectivity.zn2=NA,
                  SEL_F=NA,SEL_M=NA,SEL_6.5=NA,SEL_7=NA,
                  STEEPns=NA,age_F=NA,age_M=NA,TWT_F=NA,TWT_M=NA,
                  iyr=NA,CATCH=NA,ct_F=NA,ct_M=NA,
                  n.size.yrs.F=NA,n.size.yrs.M=NA,
                  n.size.yrs.F_7=NA,n.size.yrs.M_7=NA,
                  SIZE.comp.yrs=NA,
                  n.SIZE.comp.yrs.F=NA,n.SIZE.comp.yrs.M=NA,
                  n.SIZE.comp.yrs.F_7=NA,n.SIZE.comp.yrs.M_7=NA,
                  First.size.obs=NA,Last.size.obs=NA,First.size.obs_7=NA,Last.size.obs_7=NA,
                  First.yr.size.obs=NA,Last.yr.size.obs=NA,
                  SIZE.comp=NA,SIZE.comp.F=NA,SIZE.comp.M=NA,
                  SIZE.comp_7=NA,SIZE.comp.F_7=NA,SIZE.comp.M_7=NA,
                  size.wst.yrs=NA,size.zn1.yrs=NA,size.zn2.yrs=NA,
                  size.wst=NA,size.zn1=NA,size.zn2=NA,
                  size.wst.f=NA,size.zn1.f=NA,size.zn2.f=NA,
                  size.wst.m=NA,size.zn1.m=NA,size.zn2.m=NA,
                  size.wst_7=NA,size.zn1_7=NA,size.zn2_7=NA,
                  size.wst.f_7=NA,size.zn1.f_7=NA,size.zn2.f_7=NA,
                  size.wst.m_7=NA,size.zn1.m_7=NA,size.zn2.m_7=NA,
                  rho=NA,Type_of_Size_like=NA,
                  n.KTCH.wght.yrs=NA,KTCH.wght.yrs=NA,KTCH.wght=NA,KTCH.wght.CV=NA,
                  n.KTCH.wght.yrs.WC=NA,n.KTCH.wght.yrs.zn1=NA,n.KTCH.wght.yrs.zn2=NA,
                  KTCH.wght.yrs.WC=NA,KTCH.wght.yrs.zn1=NA,KTCH.wght.yrs.zn2=NA,
                  KTCH.wght.WC=NA,KTCH.wght.zn1=NA,KTCH.wght.zn2=NA,
                  KTCH.wght.CV.WC=NA,KTCH.wght.CV.zn1=NA,KTCH.wght.CV.zn2=NA,
                  EFFORT=NA,
                  do_cpue=NA,effective_cpue=NA,syr_cpue=NA,CPUE_eff=NA,CPUE=NA,CPUE.CV=NA,
                  Conv.tg.StartYear=NA,Conv.tg.EndYear=NA,Conv.tg.nyrs=NA,Conv.tg.Areas=NA,
                  Conv.tg.numTagGp=NA,Conv.tg.rel=NA,Conv.tg.Releases_yr=NA,
                  Conv.tg.nzones=NA,Conv.tg.recaptures=NA,Conv.tg.rec=NA,
                  Conv.tg.rec_size_index=NA,Reporting=NA,                
                  Shedding=NA,Conv.tg.mov.pars=NA,Conv.tg.Like.type=NA,Conv.tg.constant=NA,                   
                  Conv.tg.rel_adul=NA,Conv.tg.rel_juv=NA,
                  Conv.tg.rec_adul=NA,Conv.tg.rec_juv=NA,Smallest_move=NA, 
                  Acous.tg.detections=NA,Acous.tg.tags=NA,
                  Acous.tg.obs.per.tags=NA,Acous.tg.start=NA,Acous.tg.end=NA,Acous.tg.data=NA,
                  Acous.tg.rel=NA,Acous.tg.rec=NA,
                  Rec.error=NA,
                  q_change=NA,q_daily=NA,
                  Nobs.Age.Grw.F=NA,Nobs.Age.Grw.M=NA,Age.Grw.F=NA,Age.Grw.M=NA,Rho2=NA,Rho3=NA,
                  MaxF=NA,add_Finit_prior=NA,
                  mu_F_init=NA,Ln_SD_F_init=NA,Po=NA,r_max=NA,add_r_prior=NA,mu_r=NA,SD_r=NA,Mov_trans_mat=NA,Move_age=NA,
                  Phases=NA,n_par=NA,Do_var=NA,Var1=NA,Var2=NA,Do.move=NA,
                  Calc_MSY=0,yrs_Fishing.mort=1,Fishing.mort=0,Rec.error_MSY=0)
    
    #1. Add values
    #1.1. Number of years of catch
    La.lista$iyr=as.numeric(substr(Ktch.All.1975[,match('FINYEAR',names(Ktch.All.1975))],1,4))
    
    #use cpue or not in Length-based model
    if(d$Model_type=="Length-based")La.lista$do_cpue=ifelse(d$CPUE%in%c("N/A","No"),0,1)
    
    #1.2 Use effective or stand cpue
    La.lista$effective_cpue=ifelse(d$CPUE=="Effective",1,0)
    
    #years of cpue data
    if(d$Model_type%in%c("Biomass dynamics","Age-structured"))
    {
      La.lista$syr_cpue=yr.start
     # if(La.lista$effective_cpue==0 & !Yr_q_change==0) La.lista$syr_cpue=Yr_q_change
    }
      
    if(d$Model_type%in%c("Length-based"))
    {
      La.lista$syr_cpue=1
     # if(La.lista$effective_cpue==0 & !Yr_q_change==0) La.lista$syr_cpue=match(Yr_q_change,La.lista$iyr)
    }
      
    
       
    #1.3 Phases for estimable pars
    La.lista$Phases=Par.phases[[i]]
    
    #number of estimable pars for Age-structured
    if(d$Model_type=="Age-structured")
    {
      La.lista$n_par=length(La.lista$Phases[La.lista$Phases>0])
      La.lista$Do_var=Do_var
      La.lista$Var1=Var1
      La.lista$Var2=Var2
    }
      
    
    #1.4 Select years of cpue used
    NN=1:nrow(Cpue.all)
    if(d$CPUE=='Stand.1988')
    {      
      ID=match("1988-89",Cpue.all$Finyear)   
      NN=ID:nrow(Cpue.all)  
      La.lista$syr_cpue=1988
    }
    
    #1.5 Growth data
    if(d$Age.Growth=="Yes")  
    {
      if(exists("Age.growth.F"))
      {
        La.lista$Nobs.Age.Grw.F=nrow(Age.growth.F)  
        La.lista$Age.Grw.F=Age.growth.F
      }
      
      if(exists("Age.growth.M"))
      {
        La.lista$Nobs.Age.Grw.M=nrow(Age.growth.M)    
        La.lista$Age.Grw.M=Age.growth.M
      }
      
      La.lista$Rho2=Rho2
      La.lista$Rho3=Rho3[i]
    }
    
    #1.6 Size composition data
    if(d$Size_comp.=="Yes")
    {
      #number of size bins
      La.lista$nTL=length(TL.bin.mid)
      
      #weight of size composition likelihood
      La.lista$rho=Rho
      
      #define what type of LL to use
      if(Size_like=="Multinomial") Type_of_Size_like=0
      if(Size_like=="Dirichlet") Type_of_Size_like=1
      La.lista$Type_of_Size_like=Type_of_Size_like
      
      #define sizes with data
      First.size.obs=colSums(fn.add.missing.size.year(size.all))
      dd=names(which(First.size.obs>0))
      La.lista$First.size.obs=match(dd[1],names(First.size.obs))
      La.lista$Last.size.obs=match(dd[length(dd)],names(First.size.obs))
        
      First.size.obs_7=colSums(fn.add.missing.size.year(size.all_7))
      dd=names(which(First.size.obs_7>0))
      La.lista$First.size.obs_7=match(dd[1],names(First.size.obs_7))
      La.lista$Last.size.obs_7=match(dd[length(dd)],names(First.size.obs_7))
 
      #define years with data   
      First.yr.obs=c(match(size.all$FINYEAR,Ktch.All.1975$FINYEAR),
                     match(size.all_7$FINYEAR,Ktch.All.1975$FINYEAR))
      
      La.lista$First.yr.size.obs=min(First.yr.obs)
      La.lista$Last.yr.size.obs=max(First.yr.obs)
      
  
        
      #Combined sexes     
      if(Size.sex.comb=="YES") 
      {
        if(d$Spatial_structure=="Single zone")
        {
          La.lista$SIZE.comp=fn.add.missing.size.year(size.all)
          La.lista$SIZE.comp_7=fn.add.missing.size.year(size.all_7)
          if(Size.comp.prop=="YES")
          {
            La.lista$SIZE.comp=fn.size.prop(La.lista$SIZE.comp)
            La.lista$SIZE.comp_7=fn.size.prop(La.lista$SIZE.comp_7)
          }
          
        }
        if(d$Spatial_structure=="Three zones")
        {
          La.lista$size.wst=fn.add.missing.size.year(size.wst)
          La.lista$size.zn1=fn.add.missing.size.year(size.zn1)
          La.lista$size.zn2=fn.add.missing.size.year(size.zn2)
          
          La.lista$size.wst_7=fn.add.missing.size.year(size.wst_7)
          La.lista$size.zn1_7=fn.add.missing.size.year(size.zn1_7)
          La.lista$size.zn2_7=fn.add.missing.size.year(size.zn2_7)
          
          
          if(Size.comp.prop=="YES")
          {
            La.lista$size.wst=fn.size.prop(La.lista$size.wst)
            La.lista$size.zn1=fn.size.prop(La.lista$size.zn1)
            La.lista$size.zn2=fn.size.prop(La.lista$size.zn2)
            
            La.lista$size.wst_7=fn.size.prop(La.lista$size.wst_7)
            La.lista$size.zn1_7=fn.size.prop(La.lista$size.zn1_7)
            La.lista$size.zn2_7=fn.size.prop(La.lista$size.zn2_7)
            
          }
          
        } 
      }

      #By sex  
      if(Size.sex.comb=="NO")    
      {
          #Combined zones
        if(d$Spatial_structure=="Single zone")
        {
          La.lista$SIZE.comp.F=fn.add.missing.size.year(size.all.fem)
          La.lista$SIZE.comp.M=fn.add.missing.size.year(size.all.mal)
          
          La.lista$SIZE.comp.F_7=fn.add.missing.size.year(size.all.fem_7)
          La.lista$SIZE.comp.M_7=fn.add.missing.size.year(size.all.mal_7)
          
          N.yrs.F=rowSums(La.lista$SIZE.comp.F)
          N.yrs.M=rowSums(La.lista$SIZE.comp.M)
          La.lista$n.SIZE.comp.yrs.F=mapply(function(x) min(x,Effective.n),N.yrs.F)
          La.lista$n.SIZE.comp.yrs.M=mapply(function(x) min(x,Effective.n),N.yrs.M)
          
          N.yrs.F=rowSums(La.lista$SIZE.comp.F_7)
          N.yrs.M=rowSums(La.lista$SIZE.comp.M_7)
          La.lista$n.SIZE.comp.yrs.F_7=mapply(function(x) min(x,Effective.n),N.yrs.F)
          La.lista$n.SIZE.comp.yrs.M_7=mapply(function(x) min(x,Effective.n),N.yrs.M)
          
          
          if(Size.comp.prop=="YES")
          {
            La.lista$SIZE.comp.F=fn.size.prop(La.lista$SIZE.comp.F)
            La.lista$SIZE.comp.M=fn.size.prop(La.lista$SIZE.comp.M)   
            La.lista$SIZE.comp.F_7=fn.size.prop(La.lista$SIZE.comp.F_7)
            La.lista$SIZE.comp.M_7=fn.size.prop(La.lista$SIZE.comp.M_7)            
            
          }   
        }
        
          #By zone
        if(d$Spatial_structure=="Three zones")
        {
          La.lista$size.wst.f=fn.add.missing.size.year(size.wst.fem)
          La.lista$size.zn1.f=fn.add.missing.size.year(size.zn1.fem)
          La.lista$size.zn2.f=fn.add.missing.size.year(size.zn2.fem)
          
          La.lista$size.wst.m=fn.add.missing.size.year(size.wst.mal)
          La.lista$size.zn1.m=fn.add.missing.size.year(size.zn1.mal)
          La.lista$size.zn2.m=fn.add.missing.size.year(size.zn2.mal)
          
          La.lista$size.wst.f_7=fn.add.missing.size.year(size.wst.fem_7)
          La.lista$size.zn1.f_7=fn.add.missing.size.year(size.zn1.fem_7)
          La.lista$size.zn2.f_7=fn.add.missing.size.year(size.zn2.fem_7)
          
          La.lista$size.wst.m_7=fn.add.missing.size.year(size.wst.mal_7)
          La.lista$size.zn1.m_7=fn.add.missing.size.year(size.zn1.mal_7)
          La.lista$size.zn2.m_7=fn.add.missing.size.year(size.zn2.mal_7)
          
          
          N.yrs.wst.F=rowSums(La.lista$size.wst.f) 
          N.yrs.zn1.F=rowSums(La.lista$size.zn1.f)
          N.yrs.zn2.F=rowSums(La.lista$size.zn2.f)
          
          N.yrs.wst.M=rowSums(La.lista$size.wst.m) 
          N.yrs.zn1.M=rowSums(La.lista$size.zn1.m)
          N.yrs.zn2.M=rowSums(La.lista$size.zn2.m)
          
          n.size.wst.yrs.F=mapply(function(x) min(x,Effective.n),N.yrs.wst.F)
          n.size.zn1.yrs.F=mapply(function(x) min(x,Effective.n),N.yrs.zn1.F)
          n.size.zn2.yrs.F=mapply(function(x) min(x,Effective.n),N.yrs.zn2.F)    
          
          n.size.wst.yrs.M=mapply(function(x) min(x,Effective.n),N.yrs.wst.M) 
          n.size.zn1.yrs.M=mapply(function(x) min(x,Effective.n),N.yrs.zn1.M)
          n.size.zn2.yrs.M=mapply(function(x) min(x,Effective.n),N.yrs.zn2.M)
          
          La.lista$n.size.yrs.F=cbind(n.size.wst.yrs.F,n.size.zn1.yrs.F,n.size.zn2.yrs.F)
          La.lista$n.size.yrs.M=cbind(n.size.wst.yrs.M,n.size.zn1.yrs.M,n.size.zn2.yrs.M)
          
          
          N.yrs.wst.F=rowSums(La.lista$size.wst.f_7) 
          N.yrs.zn1.F=rowSums(La.lista$size.zn1.f_7)
          N.yrs.zn2.F=rowSums(La.lista$size.zn2.f_7)
          
          N.yrs.wst.M=rowSums(La.lista$size.wst.m_7) 
          N.yrs.zn1.M=rowSums(La.lista$size.zn1.m_7)
          N.yrs.zn2.M=rowSums(La.lista$size.zn2.m_7)
          
          n.size.wst.yrs.F=mapply(function(x) min(x,Effective.n),N.yrs.wst.F)
          n.size.zn1.yrs.F=mapply(function(x) min(x,Effective.n),N.yrs.zn1.F)
          n.size.zn2.yrs.F=mapply(function(x) min(x,Effective.n),N.yrs.zn2.F)   
          
          n.size.wst.yrs.M=mapply(function(x) min(x,Effective.n),N.yrs.wst.M) 
          n.size.zn1.yrs.M=mapply(function(x) min(x,Effective.n),N.yrs.zn1.M)
          n.size.zn2.yrs.M=mapply(function(x) min(x,Effective.n),N.yrs.zn2.M)
          
          La.lista$n.size.yrs.F_7=cbind(n.size.wst.yrs.F,n.size.zn1.yrs.F,n.size.zn2.yrs.F)
          La.lista$n.size.yrs.M_7=cbind(n.size.wst.yrs.M,n.size.zn1.yrs.M,n.size.zn2.yrs.M)
          
          
          if(Size.comp.prop=="YES")
          {
            La.lista$size.wst.f=fn.size.prop(La.lista$size.wst.f)
            La.lista$size.zn1.f=fn.size.prop(La.lista$size.zn1.f)
            La.lista$size.zn2.f=fn.size.prop(La.lista$size.zn2.f)
            
            La.lista$size.wst.m=fn.size.prop(La.lista$size.wst.m)
            La.lista$size.zn1.m=fn.size.prop(La.lista$size.zn1.m)
            La.lista$size.zn2.m=fn.size.prop(La.lista$size.zn2.m)
            
            La.lista$size.wst.f_7=fn.size.prop(La.lista$size.wst.f_7)
            La.lista$size.zn1.f_7=fn.size.prop(La.lista$size.zn1.f_7)
            La.lista$size.zn2.f_7=fn.size.prop(La.lista$size.zn2.f_7)
            
            La.lista$size.wst.m_7=fn.size.prop(La.lista$size.wst.m_7)
            La.lista$size.zn1.m_7=fn.size.prop(La.lista$size.zn1.m_7)
            La.lista$size.zn2.m_7=fn.size.prop(La.lista$size.zn2.m_7)
            
          }
        }
      }

        
        #La.lista$size.wst.yrs=get.yr(size.wst[,match('FINYEAR',names(size.wst))])
        #La.lista$size.zn1.yrs=get.yr(size.zn1[,match('FINYEAR',names(size.zn1))])
        #La.lista$size.zn2.yrs=get.yr(size.zn2[,match('FINYEAR',names(size.zn2))])
      
      # La.lista$SIZE.comp.yrs=get.yr(size.all[,match('FINYEAR',names(size.all))])
      # if(d$Catch_ave._weight=="Yes")
      #  {
      #   aaa=fn.add.missing.year(avg.wt)
      #   La.lista$KTCH.wght=aaa$Mean  
      #   La.lista$KTCH.wght.CV=aaa$CV
      #   #La.lista$KTCH.wght.yrs=get.yr(avg.wt[,match('finyear',names(avg.wt))])
      #   #La.lista$n.KTCH.wght.yrs=length(La.lista$KTCH.wght.yrs)
      #  }
      #if(d$Tagging=="Yes") La.lista$EFFORT=Ef.zn[,-match('FINYEAR',names(Ef.zn))]
    }
    
    
    #1.7 Steepness
    if(!d$SteepnesS=="N/A") La.lista$STEEPns=as.numeric(d$SteepnesS)
    
    #1.8 By model inputs
    if(d$Model_type=="Biomass dynamics")
    {
      #Q
      La.lista$q_change=Yr_q_change  
      La.lista$q_daily=Yr_q_daily
      
      #Catch in tonnes
      La.lista$CATCH=Ktch.All.1975[,match('LIVEWT.c',names(Ktch.All.1975))]/1000   
      
      #add future years 
      La.lista$CATCH=c(La.lista$CATCH,Ktch.future(La.lista$CATCH))
      La.lista$yrs_ktch=(yr.start:(yr.end+Yrs.future))[length(La.lista$CATCH)]
      
      #initial F for start of catch series
      La.lista$Po=Po_spm
      
      #r priors
      La.lista$r_max=r_max
      La.lista$add_r_prior=Add.r.prior
      La.lista$mu_r=log_mean.r 
      La.lista$SD_r=log_sd.r
      
      #CPUE
      if(d$CPUE=='Stand.')
      {
        La.lista$CPUE=Cpue.all[NN,match('Mean',names(Cpue.all))]     
        La.lista$CPUE.CV=Cpue.all[NN,match('CV',names(Cpue.all))]
        
      }
    }
    if(d$Model_type=="Age-structured")
    {
      # Q
      La.lista$q_change=Yr_q_change
      La.lista$q_daily=Yr_q_daily

      # Life history parameters
      La.lista$age_F=age.F
      La.lista$age_M=age.M
      La.lista$TWT_F=mid.TwT.F
      La.lista$TWT_M=mid.TwT.M
      La.lista$nage_F=length(age.F)
      La.lista$nage_M=length(age.M)
      if(d$Fec.=="constant") La.lista$fec=Litter.size.k 
      if(d$Fec.=="at age") La.lista$fec=Litter.size.age 
      if(d$Maturity=="knife edge")La.lista$age_mat=Age.50.mat.min
      if(d$Maturity=="at age")La.lista$age_mat=Mat.age 
      if(d$M=="constant") La.lista$M=rep(unique(d$M.value),length(Mort.F.k))
      if(d$M=="at age") La.lista$M=Mort.F.age
      La.lista$breed=Breed.freq.min
      La.lista$sx_ratio=pup.sx.ratio
      
      # Selectivity   
      La.lista$SEL_F=Select.age.F
      La.lista$SEL_M=Select.age.M
      
      # Catches by sex (in kg)
        # single zone
      if(d$Spatial_structure=="Single zone")
      {
        if(d$Ktch.sx.r=="Observed")
        {
          La.lista$ct_F=Ktch.All.1975[,match('LIVEWT.c',names(Ktch.All.1975))]*(1-Prop.males.ktch)
          La.lista$ct_M=Ktch.All.1975[,match('LIVEWT.c',names(Ktch.All.1975))]*Prop.males.ktch
        }
        if(d$Ktch.sx.r=="Equal")
        {
          La.lista$ct_F=Ktch.All.1975[,match('LIVEWT.c',names(Ktch.All.1975))]*0.5
          La.lista$ct_M=Ktch.All.1975[,match('LIVEWT.c',names(Ktch.All.1975))]*0.5
        }
        
        # add future years of catch by sex
        La.lista$ct_F=c(La.lista$ct_F,Ktch.future(La.lista$ct_F))
        La.lista$ct_M=c(La.lista$ct_M,Ktch.future(La.lista$ct_M)) 
        
        La.lista$yrs_ktch=(yr.start:(yr.end+Yrs.future))[length(La.lista$ct_F)]
      }
        # Three zones
      if(d$Spatial_structure=="Three zones")
      {
        if(d$Ktch.sx.r=="Observed")
        {
          Ktch.F.WC=(1-Prop.males.ktch.wst)*Ktch.All.West.1975[,match('West',names(Ktch.All.West.1975))]
          Ktch.F.zn1=(1-Prop.males.ktch.zn1)*Ktch.All.zn1.1975[,match('Zone1',names(Ktch.All.zn1.1975))]
          Ktch.F.zn2=(1-Prop.males.ktch.zn2)*Ktch.All.zn2.1975[,match('Zone2',names(Ktch.All.zn2.1975))]
          
          Ktch.M.WC=Prop.males.ktch.wst*Ktch.All.West.1975[,match('West',names(Ktch.All.West.1975))]
          Ktch.M.zn1=Prop.males.ktch.zn1*Ktch.All.zn1.1975[,match('Zone1',names(Ktch.All.zn1.1975))]
          Ktch.M.zn2=Prop.males.ktch.zn2*Ktch.All.zn2.1975[,match('Zone2',names(Ktch.All.zn2.1975))]
        }
        if(d$Ktch.sx.r=="Equal")
        {
          Ktch.F.WC=0.5*Ktch.All.West.1975[,match('West',names(Ktch.All.West.1975))]
          Ktch.F.zn1=0.5*Ktch.All.zn1.1975[,match('Zone1',names(Ktch.All.zn1.1975))]
          Ktch.F.zn2=0.5*Ktch.All.zn2.1975[,match('Zone2',names(Ktch.All.zn2.1975))]
          
          Ktch.M.WC=0.5*Ktch.All.West.1975[,match('West',names(Ktch.All.West.1975))]
          Ktch.M.zn1=0.5*Ktch.All.zn1.1975[,match('Zone1',names(Ktch.All.zn1.1975))]
          Ktch.M.zn2=0.5*Ktch.All.zn2.1975[,match('Zone2',names(Ktch.All.zn2.1975))]
        }
        #add future years to base case
        Ktch.F.WC=c(Ktch.F.WC,Ktch.future(Ktch.F.WC))
        Ktch.F.zn1=c(Ktch.F.zn1,Ktch.future(Ktch.F.zn1))
        Ktch.F.zn2=c(Ktch.F.zn2,Ktch.future(Ktch.F.zn2))   
        
        Ktch.M.WC=c(Ktch.M.WC,Ktch.future(Ktch.M.WC))
        Ktch.M.zn1=c(Ktch.M.zn1,Ktch.future(Ktch.M.zn1))
        Ktch.M.zn2=c(Ktch.M.zn2,Ktch.future(Ktch.M.zn2))      
        
        
        La.lista$ct_F=data.frame(WC=Ktch.F.WC,zn1=Ktch.F.zn1,zn2=Ktch.F.zn2)
        La.lista$ct_M=data.frame(WC=Ktch.M.WC,zn1=Ktch.M.zn1,zn2=Ktch.M.zn2)
        
        La.lista$yrs_ktch=(yr.start:(yr.end+Yrs.future))[nrow(La.lista$ct_F)]
      }
      
      # CPUE
        # Single zone
      if(d$Spatial_structure=="Single zone")
      {
        if(d$CPUE=="Effective") La.lista$CPUE=Cpue.folly[,match("cpue",names(Cpue.folly))] 
        if(d$CPUE=="Stand.") La.lista$CPUE=Cpue.all$Mean
        if(d$CPUE=="Stand.hours") La.lista$CPUE=Cpue.all.hours$Mean
      }
      
      # Recruitment error for future projections
      La.lista$Rec.error=rlnorm(1+La.lista$yrs_ktch-yr.start, meanlog = log(1), sdlog = MSY.sd.rec)
      
    }
    if(d$Model_type=="Length-based")
    {
      # Fint
      La.lista$MaxF=MaxF
      La.lista$add_Finit_prior=add_Finit_prior
      La.lista$mu_F_init=Prior.mean.Fo
      La.lista$Ln_SD_F_init=Prior.SD.Log.Fo
      
      # add the empirical selectivity of 6.5 and 7 inch mesh
      if(d$Size_comp.=="Yes")  
      {
        La.lista$SEL_6.5=Select.bin
        La.lista$SEL_7=Select.bin_7
      }
      
      # apply movement or not
      if(d$Movement%in%c("N/A","No")) La.lista$Do.move=0
      if(d$Movement%in%c("Yes")) La.lista$Do.move=1
      
      # Q
      QQQ=match(Yr_q_change,La.lista$iyr)
      if(is.na(QQQ)) QQQ=Yr_q_change
      La.lista$q_change=QQQ
      La.lista$q_daily=match(Yr_q_daily,La.lista$iyr)
      
        # Life history parameters
      La.lista$mean.size.birth=Mean.size.birth.TL 
      La.lista$SD.mean.size.birth=Size.birth_SD$Value
      La.lista$TL.bin.mid=TL.bin.mid
      La.lista$TL.bin.low=TL.bin.low
      La.lista$TL.bin.up=TL.bin.up
      La.lista$Max.age.size=2*Max.age.F       #max age for plus group    
      La.lista$TwT.bin.F=TwT.bin.fem
      La.lista$TwT.bin.M=TwT.bin.mal
      if(d$Maturity=="at length")  La.lista$size_mat=Mat.bin      
      if(d$Fec.=="at length") La.lista$fec=Litter.size.bin 
      if(d$M=="at length")  La.lista$M=Mort.bin   
      if(d$Maturity=="knife edge") La.lista$size_mat=ifelse(TL.bin.mid>=Maturity.at.size[1],1,0)
      if(d$Fec.=="constant") La.lista$fec=ifelse(TL.bin.mid>=Maturity.at.size[1],Litter.size.k,0)
      if(d$M=="constant") La.lista$M=rep(unique(d$M.value),length(TL.bin.mid))
      La.lista$sx_ratio=pup.sx.ratio
      
      # Selectivity
      La.lista$Selectivity=Total.sel.all
      La.lista$Selectivity.west=Total.sel.West
      La.lista$Selectivity.zn1=Total.sel.Zn1
      La.lista$Selectivity.zn2=Total.sel.Zn2
      
      # Conventional tagging data
      La.lista$Smallest_move=Smallest_size_tagged.index 
      if(conv.tag.all=="YES" & Move.mode=="Individual-based")
      {
        La.lista$Conv.tg.nzones=Conv.tg.nzones
        La.lista$Conv.tg.recaptures=Conv.tg.recs
        La.lista$Conv.tg.rec=Conv.tg.rec.exp 
      }
      if(conv.tag.all=="YES" & Move.mode=="Population-based")                                         
      {
        La.lista$Reporting=Reporting
        La.lista$Shedding=Shedding
        La.lista$Conv.tg.rel=Conv.tg.rel_size
        La.lista$Conv.tg.rec=Conv.tg.rec_size
        
        La.lista$Conv.tg.StartYear=Conv.tg.StartYear
        La.lista$Conv.tg.EndYear=Conv.tg.EndYear
        La.lista$Conv.tg.nyrs=Conv.tg.nyrs
        La.lista$Conv.tg.Areas=Conv.tg.Areas
        La.lista$Conv.tg.numTagGp=Conv.tg.numTagGp
        La.lista$Conv.tg.Releases_yr=Conv.tg.Releases_yr
        La.lista$Conv.tg.rec_size_index= Conv.tg.rec_size_index    
        La.lista$Conv.tg.mov.pars=Conv.tg.mov.pars
        La.lista$Conv.tg.Like.type=Conv.tg.Like.type
        La.lista$Conv.tg.constant =Conv.tg.constant
      }
      if(conv.tag.all=="NO")
      {
        La.lista$Conv.tg.rel_adul=Conv.tg.rel_size_adul
        La.lista$Conv.tg.rel_juv=Conv.tg.rel_size_juv
        La.lista$Conv.tg.rec_adul=Conv.tg.rec_size_adul
        La.lista$Conv.tg.rec_juv=Conv.tg.rec_size_juv
      }
      
      # Acoustic tagging
      if(Acoust.format=="Individual-based")
      {
        La.lista$Acous.tg.detections=Acous.detections
        La.lista$Acous.tg.tags=Acous.ntags
        La.lista$Acous.tg.obs.per.tags=Acous.TagID_observations
        La.lista$Acous.tg.start=Acous.TagID_start
        La.lista$Acous.tg.end=Acous.TagID_end
        La.lista$Acous.tg.data=Acous.tg.rec
      }
      if(Acoust.format=="SS3")
      {
        La.lista$Acous.tg.rel=Acous.tg.rel             
        La.lista$Acous.tg.rec=Acous.tg.rec
      }
       
      # Catches by sex (in tonnes)
        # single zone
      if(d$Spatial_structure=="Single zone")
      {
        if(d$Ktch.sx.r=="Observed")
        {
          La.lista$ct_F=(Ktch.All.1975[,match('LIVEWT.c',names(Ktch.All.1975))]*(1-Prop.males.ktch))/1000
          La.lista$ct_M=(Ktch.All.1975[,match('LIVEWT.c',names(Ktch.All.1975))]*Prop.males.ktch)/1000
        }
        if(d$Ktch.sx.r=="Equal")
        {
          La.lista$ct_F=(Ktch.All.1975[,match('LIVEWT.c',names(Ktch.All.1975))]*0.5)/1000
          La.lista$ct_M=(Ktch.All.1975[,match('LIVEWT.c',names(Ktch.All.1975))]*0.5)/1000
        }
        
        # add future years of catch by sex
        La.lista$ct_F=c(La.lista$ct_F,Ktch.future(La.lista$ct_F))
        La.lista$ct_M=c(La.lista$ct_M,Ktch.future(La.lista$ct_M))      
        
        La.lista$yrs_ktch=length(La.lista$ct_F)
        
      }
        # Three zones
      if(d$Spatial_structure=="Three zones")
      {
        if(d$Ktch.sx.r=="Observed")
        {
          Ktch.F.WC=(1-Prop.males.ktch.wst)*Ktch.All.West.1975[,match('West',names(Ktch.All.West.1975))]/1000
          Ktch.F.zn1=(1-Prop.males.ktch.zn1)*Ktch.All.zn1.1975[,match('Zone1',names(Ktch.All.zn1.1975))]/1000
          Ktch.F.zn2=(1-Prop.males.ktch.zn2)*Ktch.All.zn2.1975[,match('Zone2',names(Ktch.All.zn2.1975))]/1000
          
          Ktch.M.WC=Prop.males.ktch.wst*Ktch.All.West.1975[,match('West',names(Ktch.All.West.1975))]/1000
          Ktch.M.zn1=Prop.males.ktch.zn1*Ktch.All.zn1.1975[,match('Zone1',names(Ktch.All.zn1.1975))]/1000
          Ktch.M.zn2=Prop.males.ktch.zn2*Ktch.All.zn2.1975[,match('Zone2',names(Ktch.All.zn2.1975))]/1000
        }
        if(d$Ktch.sx.r=="Equal")
        {
          Ktch.F.WC=0.5*Ktch.All.West.1975[,match('West',names(Ktch.All.West.1975))]/1000
          Ktch.F.zn1=0.5*Ktch.All.zn1.1975[,match('Zone1',names(Ktch.All.zn1.1975))]/1000
          Ktch.F.zn2=0.5*Ktch.All.zn2.1975[,match('Zone2',names(Ktch.All.zn2.1975))]/1000
          
          Ktch.M.WC=0.5*Ktch.All.West.1975[,match('West',names(Ktch.All.West.1975))]/1000
          Ktch.M.zn1=0.5*Ktch.All.zn1.1975[,match('Zone1',names(Ktch.All.zn1.1975))]/1000
          Ktch.M.zn2=0.5*Ktch.All.zn2.1975[,match('Zone2',names(Ktch.All.zn2.1975))]/1000
        }
        #add future years to base case
        Ktch.F.WC=c(Ktch.F.WC,Ktch.future(Ktch.F.WC))
        Ktch.F.zn1=c(Ktch.F.zn1,Ktch.future(Ktch.F.zn1))
        Ktch.F.zn2=c(Ktch.F.zn2,Ktch.future(Ktch.F.zn2))   
        
        Ktch.M.WC=c(Ktch.M.WC,Ktch.future(Ktch.M.WC))
        Ktch.M.zn1=c(Ktch.M.zn1,Ktch.future(Ktch.M.zn1))
        Ktch.M.zn2=c(Ktch.M.zn2,Ktch.future(Ktch.M.zn2))      
        
        
        La.lista$ct_F=data.frame(WC=Ktch.F.WC,zn1=Ktch.F.zn1,zn2=Ktch.F.zn2)
        La.lista$ct_M=data.frame(WC=Ktch.M.WC,zn1=Ktch.M.zn1,zn2=Ktch.M.zn2)
        
        La.lista$yrs_ktch=nrow(La.lista$ct_F)
      }
      
      # CPUE
         #Single zone
      if(d$Spatial_structure=="Single zone")
      {
        La.lista$CPUE=Cpue.all[NN,match('Mean',names(Cpue.all))]     
        La.lista$CPUE.CV=Cpue.all[NN,match('CV',names(Cpue.all))]  
        if(exists('Cpue.folly')) La.lista$CPUE_eff=Cpue.folly[,match("cpue",names(Cpue.folly))] else
          La.lista$CPUE_eff=rep(1,length(La.lista$CPUE))
        if(d$CPUE=="Stand.hours")
        {
          La.lista$CPUE=Cpue.all.hours[NN,match('Mean',names(Cpue.all.hours))]
          La.lista$CPUE.CV=Cpue.all.hours[NN,match('CV',names(Cpue.all.hours))]
        }
           
      }
      
        # Three zones 
      if(d$Spatial_structure=="Three zones")   
      {
        list_cpue=list(Cpue.West,Cpue.zn1,Cpue.zn2)
        names(list_cpue)=c("WC","Zn1","Zn2")
        ThiS=sapply(list_cpue, function(x) nrow(x)>0)
        copy.this=which(ThiS==TRUE)[1]
        create.this=which(!ThiS==TRUE)
        if(length(create.this)>0)
        {
          for(i in 1:length(create.this))
          {
            a=list_cpue[[copy.this]]
            a[,2:ncol(a)]=-100
            list_cpue[[create.this[i]]]=a
          }  
        }
        All.YrS=data.frame(Finyear=unique(Cpue.all$Finyear))
        D=merge(All.YrS,list_cpue$WC[,match(c('Finyear','Mean'),names(list_cpue$WC))],by="Finyear",all.x=T)
        D=merge(D,list_cpue$Zn1[,match(c('Finyear','Mean'),names(list_cpue$Zn1))],by="Finyear",all.x=T)
        D=merge(D,list_cpue$Zn2[,match(c('Finyear','Mean'),names(list_cpue$Zn2))],by="Finyear",all.x=T)
        colnames(D)=c("Finyear","WC","zn1","zn2")
        La.lista$CPUE=D[,-match("Finyear",names(D))]
        La.lista$CPUE[is.na(La.lista$CPUE)]=-100
        
        D=merge(All.YrS,list_cpue$WC[,match(c('Finyear','CV'),names(list_cpue$WC))],by="Finyear",all.x=T)
        D=merge(D,list_cpue$Zn1[,match(c('Finyear','CV'),names(list_cpue$Zn1))],by="Finyear",all.x=T)
        D=merge(D,list_cpue$Zn2[,match(c('Finyear','CV'),names(list_cpue$Zn2))],by="Finyear",all.x=T)
        colnames(D)=c("Finyear","WC","zn1","zn2")
        La.lista$CPUE.CV=D[,-match("Finyear",names(D))] 
        La.lista$CPUE.CV[is.na(La.lista$CPUE.CV)]=-100


        if(exists('Cpue.folly')) La.lista$CPUE_eff=Cpue.folly[,match("cpue",names(Cpue.folly))] else
          La.lista$CPUE_eff=rep(1,length(La.lista$CPUE))
      }
      
      # Recruitment error for future projections
      La.lista$Rec.error=rlnorm(La.lista$yrs_ktch, meanlog = log(1), sdlog = MSY.sd.rec)
      
    }
    
    #Remove cpue years if applicable   
    match("CPUE_years_dropped",names(d))
    if(!is.na(match("CPUE_years_dropped",names(d)))) if(!d$CPUE_years_dropped=="None")
    {
      ID.dropped.yr=c(substr(d$CPUE_years_dropped,1,4),paste("19",substr(d$CPUE_years_dropped,6,7),sep=""))
      Chck.yr=substr(Cpue.all$Finyear,1,4)
      DRop=match(ID.dropped.yr,Chck.yr)
      if(is.numeric(La.lista$CPUE))
      {
        La.lista$CPUE[DRop[1]:DRop[2]]=-100
        if(!is.na(La.lista$CPUE.CV[1]))La.lista$CPUE.CV[DRop[1]:DRop[2]]=-100
      }else
      {
        La.lista$CPUE[DRop[1]:DRop[2],]=-100
        if(d$CPUE=="Stand.")La.lista$CPUE.CV[DRop[1]:DRop[2],]=-100
      }

      
    }
    
    
    #Catch average weight   ACA use if(nrow(avg.wt.wst)>0)
    #     if(d$Catch_ave._weight=="Yes") 
    #     {
    #       aaa=fn.add.missing.year(avg.wt.wst)
    #       La.lista$KTCH.wght.WC=aaa$Mean    
    #       La.lista$KTCH.wght.CV.WC=aaa$CV    
    #       
    #       aaa=fn.add.missing.year(avg.wt.zn1)
    #       La.lista$KTCH.wght.zn1=aaa$Mean    
    #       La.lista$KTCH.wght.CV.zn1=aaa$CV   
    #       
    #       aaa=fn.add.missing.year(avg.wt.zn2)
    #       La.lista$KTCH.wght.zn2=aaa$Mean    
    #       La.lista$KTCH.wght.CV.zn2=aaa$CV  
    #       
    # #       La.lista$KTCH.wght.yrs.WC=get.yr(avg.wt.wst[,match('finyear',names(avg.wt.wst))])
    # #       La.lista$KTCH.wght.yrs.zn1=get.yr(avg.wt.zn1[,match('finyear',names(avg.wt.zn1))])
    # #       La.lista$KTCH.wght.yrs.zn2=get.yr(avg.wt.zn2[,match('finyear',names(avg.wt.zn2))])
    #       
    # #       La.lista$n.KTCH.wght.yrs.WC=length(La.lista$KTCH.wght.yrs.WC)
    # #       La.lista$n.KTCH.wght.yrs.zn1=length(La.lista$KTCH.wght.yrs.zn1)
    # #       La.lista$n.KTCH.wght.yrs.zn2=length(La.lista$KTCH.wght.yrs.zn2)
    #       
    #     }
    
    #2. Remove NA data
    La.lista=subset(La.lista,!is.na(La.lista))
    
    #3. Store inputs
    Inputs[[i]]=La.lista
    
    #4. Create .dat file
    n=Inputs[[i]]
    
    setPath(d$Model)
    FILE=paste(Spec,".dat",sep="")
    if( d$Spatial_structure=="Single zone") nzones=1
    if( d$Spatial_structure=="Three zones") nzones=3
    if(is.data.frame(n$CATCH))nfisheries=ncol(n$CATCH)	
    ModDims=unlist(c(yr.start,yr.end,nzones))
    Hdr="#Basic model dimensions (yr.start, yr.end, nzones)"
    write(Hdr,file = FILE)
    write(ModDims,file = FILE,sep = "\t",append=T)
    for(k in 1:length(n))
    {
      nn=n[[k]]
      if(is.data.frame(nn)|is.matrix(nn))
      {
        Hdr=paste("#",paste(c(names(n)[k],"(",names(nn),")"),collapse=' '))
        write(Hdr,file = FILE,append=T)      
        write.table(nn,file = FILE,row.names=F,col.names=F,append=T)
      }else
      {
        Hdr=paste("#",names(n)[k],sep='')
        write(Hdr,file = FILE,append=T)
        write(n[[k]],file = FILE,sep = "\t",append=T)
      }
    }
    write("#eof",file = FILE,append=T)
    write(999,file = FILE,sep = "\t",append=T)
    
    #create data for scenarios on future projections
    if(d$Model=="Base case")
    {
      for(future in 1:length(Future.ktch.scen))
      {
        n=Inputs[[i]]
        Folder=paste("Future.projections",names(Future.ktch.scen)[future],sep="_")
        if(!file.exists(file.path(getwd(), Folder))) dir.create(file.path(getwd(), Folder)) 
        setwd(file.path(getwd(), Folder))
        
        #set future catches
        inDx=(length(n$iyr)+1):(length(n$iyr)+Yrs.future)
        n$ct_F[inDx]=n$ct_F[inDx]*Future.ktch.scen[[future]]
        n$ct_M[inDx]=n$ct_M[inDx]*Future.ktch.scen[[future]]
        
        
        Hdr="#Basic model dimensions (yr.start, yr.end, nzones)"
        write(Hdr,file = FILE)
        write(ModDims,file = FILE,sep = "\t",append=T)
        for(k in 1:length(n))
        {
          nn=n[[k]]
          if(is.data.frame(nn)|is.matrix(nn))
          {
            Hdr=paste("#",paste(c(names(n)[k],"(",names(nn),")"),collapse=' '))
            write(Hdr,file = FILE,append=T)      
            write.table(nn,file = FILE,row.names=F,col.names=F,append=T)
          }else
          {
            Hdr=paste("#",names(n)[k],sep='')
            write(Hdr,file = FILE,append=T)
            write(n[[k]],file = FILE,sep = "\t",append=T)
          }
        }
        write("#eof",file = FILE,append=T)
        write(999,file = FILE,sep = "\t",append=T)
        setPath(d$Model)
      }
      
    }
    rm(n)
    Inputs[[i]]=c(syr=yr.start,nyr=yr.end,nzone=nzones,Inputs[[i]],eof=999)
  }

  # Create all .tpl files
  fn.source("Copy_and_paste_files.R")
  where.tpls.are="C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Repository_of_tpl"
  FIND.NM=list.files(where.tpls.are)
  for(i in 1:nrow(Scenarios))
  {
    
    d=Scenarios[i,]
    OLD.name=FIND.NM[which(FIND.NM==paste(d$Model_type,".tpl",sep=""))]
    New.name=paste(Spec,".tpl",sep="")
    fn.copy.tpl(dir.from=where.tpls.are,
                dir.to=paste(hndl,AssessYr,"/",d$Model,sep=""),
                orgnl.nm=OLD.name,new.nm=New.name)
    
  }
  
  
  # Show catch Vs cpue 
  fn.fig(paste(hndl,AssessYr,"/1_Inputs","/Catch_Vs_cpue",sep=""),2400,2400)  
  par(xpd=T,las=1,mgp=c(2.5,1,0),mai=c(1,1,.1,1.25))
  yrs.ktch.cpue=1:length(Ktch.All.1975$FINYEAR)
  XX=Ktch.All.1975$LIVEWT.c
  if(nchar(max(round(Ktch.All.1975$LIVEWT.c)))>=5) XX=XX/1000
  plot(yrs.ktch.cpue,XX,type='l',ylab="",xlab="",xaxt='n',lwd=2,cex.axis=1.25,ylim=c(0,max(XX)))
  par(new=T)
  plot(yrs.ktch.cpue,Cpue.all$Mean, axes=F, xlab=NA, ylab=NA,pch=21,bg=2,cex=2,cex.axis=1.25,ylim=c(0,max(Cpue.all$Mean)))
  
  #add polygons corresponding to different q periods
  UP=max(c(Cpue.all.monthly$Mean,Cpue.all.daily$Mean))
  LOW=0
  y.Vec <- c(UP, UP, LOW, LOW)
  
  #monthly
  YR=match(Cpue.all.monthly$Finyear,Ktch.All.1975$FINYEAR)
  x.Vec <- c(YR[1], tail(YR, 1), tail(YR, 1), YR[1])
  polygon(x.Vec, y.Vec, col = rgb(red=0.1, green=0.3, blue=0.1, alpha=.2), border = "grey80")
  
  if(Scenarios[1,]$Q=="three")
  {
    YR=match(Yr_q_change,as.numeric(substr(Ktch.All.1975$FINYEAR,1,4)))
    x.Vec <- c(1, YR, YR, 1)
    polygon(x.Vec, y.Vec, col = rgb(red=0.1, green=0.1, blue=0.3, alpha=.2), border = "grey80")
    
  }
  
  #daily
  YR=match(Cpue.all.daily$Finyear,Ktch.All.1975$FINYEAR)
  x.Vec <- c(YR[1], tail(YR, 1), tail(YR, 1), YR[1])
  polygon(x.Vec, y.Vec, col = rgb(red=0.3, green=0.1, blue=0.1, alpha=.2), border = "grey80")
  
  mtext("Total catch (tonnes)",2,3,las=3,cex=2.5)
  mtext("Relative cpue",4,4,las=3,cex=2.5)
  mtext("Financial year",1,3,las=1,cex=2.5)
  
  axis(4,round(seq(0,max(Cpue.all$Mean),length.out = 4),2),cex.axis=1.25)
  axis(1,yrs.ktch.cpue,F,tck=-0.02)
  axis(1,seq(yrs.ktch.cpue[1],length(yrs.ktch.cpue),10),
       Ktch.All.1975$FINYEAR[seq(yrs.ktch.cpue[1],length(yrs.ktch.cpue),10)],tck=-0.04,cex.axis=1.25)
  dev.off()
  
  
  # Create table of proportion male in catch 
  Prop.male.tl=subset(TBl,Parameter=="Prop.males.in.ktch",select=c("Comment","Value"))
  names(Prop.male.tl)=c("Zone","Proportion of males")
  row.names(Prop.male.tl)=NULL
  Prop.male.tl$Zone=with(Prop.male.tl,ifelse(Zone=="WC","West Coast",
                                             ifelse(Zone=="Zn1","Zone 1",ifelse(Zone=="Zn2","Zone 2","Zones combined"))))
  
  setwd(paste(hndl,AssessYr,"/1_Inputs",sep=""))
  write.csv(as.matrix(Prop.male.tl),"Prop_male_catch.csv",row.names=F)
  
  # Create table of model scenarios
  Tab.scen=Tabla.scen.show
  write.csv(as.matrix(Tab.scen),"Tab.scen.csv",row.names=F)
  
}


# Section D: RUN MODELS ---------------------------------------------------
if(BaseCase=="Age-based") ID.base.Model="Base case_age"
if(BaseCase=="Size-based") ID.base.Model="Base case"

library(R2admb)
library(beepr)
fn.source("reptoRlist.R")
fn.source("ADMB_read.fit.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/send.emails.R")

#1. Run model scenarios
if(!First.run=="YES") Scenarios=Tabla.scen[match(ID.base.Model,Tabla.scen$Model),]   
Store.Models=vector('list',nrow(Scenarios))
names(Store.Models)=Scenarios$Model
Store.Reports=Store.std=Model.Out.sumry=Store.Models

Start.time=Sys.time()

if(Run.all.Scenarios=="NO")      
{
  #run .tpl
  i=match(ID.base.Model,Scenarios$Model)
  d=Scenarios[i,]
  setPath(d$Model)  
  args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep=""),Arguments,sep=""), sep="")
  clean_admb(Spec)
  compile_admb(Spec,verbose=T)
  system(args)
  
  #Bring in .tpl outputs
  Store.Models[[i]]=read.fit(paste(Spec))    # param estimates
  Store.Reports[[i]]=suppressWarnings(reptoRlist(paste(Spec,".rep",sep="")))  # Report file  
  Store.std[[i]]=read.table(paste(Spec,".std",sep=""), header =TRUE)
  m1=Store.Models[[i]]
  Coefs=m1$est 
  names(Coefs)=m1$names
  Correlation=m1$cor 
  NLL=m1$nlogl
  Max.gradient=m1$maxgrad
  Model.Out.sumry[[i]]=list(Coefs=Coefs,Correlation=Correlation,NLL=NLL,Max.gradient=Max.gradient)
}
if(Run.all.Scenarios=="YES")   #Size structured with movement takes 29 mins
{
  #run .tpl
  for(i in 1:nrow(Scenarios))
  {
      d=Scenarios[i,]
      setPath(d$Model)  
      args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep=""),Arguments,sep=""), sep="")
      clean_admb(Spec)
      compile_admb(Spec,verbose=T)
      system(args)
     print(paste("just finished Model",Scenarios$Model[i]))
        #Store
      Store.Models[[i]]=read.fit(paste(Spec))    # param estimates
      Store.Reports[[i]]=suppressWarnings(reptoRlist(paste(Spec,".rep",sep="")))  # Report file  
      Store.std[[i]]=read.table(paste(Spec,".std",sep=""), header =TRUE)
  }  
  
  #Bring in .tpl outputs
  for( i in 1:nrow(Scenarios))
  {
    m1=Store.Models[[i]]
    Coefs=m1$est 
    names(Coefs)=m1$names
    Correlation=m1$cor 
    NLL=m1$nlogl
    Max.gradient=m1$maxgrad
    Model.Out.sumry[[i]]=list(Coefs=Coefs,Correlation=Correlation,NLL=NLL,Max.gradient=Max.gradient)
  }
}

End.time=Sys.time()
Tot.time=round(difftime(End.time,Start.time,units="mins"),2)
  #email running time
function.send.email(
  to ="Matias.Braccini@fish.wa.gov.au",
  subject ="Model run",
    body =paste("Model run finished. All these models: ",paste(Scenarios$Model, collapse=", "),
              ", took",Tot.time,"minutes to run"),                     
  Attachment=NULL)

#2. Project population into the future with no catch
  #note: is the population rebuilding to equilibrium unfished conditions?
if(First.run=="YES")
{
  if(Do.zero.Ktch=="YES")
  {
    fn.check.model=function(MODEL,yrs.projections)
    {
      setPath(Scenarios[match(MODEL,Scenarios$Model),]$Model)
      A=getwd()
      
      if(!file.exists(file.path(getwd(), "Zero.catch"))) dir.create(file.path(getwd(), "Zero.catch"))  
      setwd(file.path(getwd(), "Zero.catch"))
      
      #create new data and pin files
      #.pin file
      These.pars=Pin.pars[[match(MODEL,names(Pin.pars))]]
      par.nms=names(These.pars)
      FILE=paste(Spec,".pin",sep="")
      write("# Input parameters",file = FILE)
      for(k in 1:length(These.pars))
      {
        Hdr=paste("#",par.nms[k])
        write(Hdr,file = FILE,append=T)
        write(These.pars[k],file = FILE,sep = "\t",append=T)
      }
      
      #.dat file  
      n=Inputs[[match(MODEL,names(Inputs))]]
      
      #add future catch and rec deviations
      Extra=rep(0,yrs.projections)  
      if(is.data.frame(n$ct_F)|is.matrix(n$ct_F))
      {
        Rep.ktch=as.data.frame(matrix(rep(Extra,ncol(n$ct_F)),ncol=ncol(n$ct_F)))
        names(Rep.ktch)=names(n$ct_F)
        n$ct_F=rbind(n$ct_F,Rep.ktch)
        n$ct_M=rbind(n$ct_M,Rep.ktch)
      }else
      {
        n$ct_F=c(n$ct_F,Extra)
        n$ct_M=c(n$ct_M,Extra)    
      }
      #expand selectivities if applicable
      x=c("Selectivity","Selectivity.west","Selectivity.zn1","Selectivity.zn2")
      Exist=match(x,names(n))
      if(!is.na(sum(Exist)))
      {
        for(ex in 1:length(Exist))
        {
          ZZ=n[[Exist[ex]]]
          ADs=ZZ[,1]
          ADs=as.data.frame(matrix(ADs,nrow=length(ADs),ncol=yrs.projections,byrow=F))
          a=as.numeric(substr(names(ZZ)[ncol(ZZ)],1,4))
          a=seq(a+1,a+yrs.projections)
          names(ADs)=paste(a,substr(a+1,3,4),sep="-")
          n[[Exist[ex]]]=cbind(ZZ,ADs)
        }
      }
      n$yrs_ktch=n$yrs_ktch+yrs.projections
      n$Rec.error=c(n$Rec.error,rlnorm(yrs.projections, meanlog = log(1), sdlog = MSY.sd.rec))   
      
      
      #Don't estimate
      #n$Phases[1]=-n$Phases[1]
      #n$Phases[2:length(n$Phases)]=-abs(n$Phases[2:length(n$Phases)])
      
      FILE=paste(Spec,".dat",sep="")
      nzones=n$nzone
      ModDims=unlist(c(yr.start,yr.end,nzones))
      Hdr="#Basic model dimensions (yr.start, yr.end, nzones)"
      write(Hdr,file = FILE)
      write(ModDims,file = FILE,sep = "\t",append=T)
      
      for(k in (length(ModDims)+1):length(n))
      {
        nn=n[[k]]
        if(is.data.frame(nn)|is.matrix(nn))
        {
          Hdr=paste("#",paste(c(names(n)[k],"(",names(nn),")"),collapse=' '))
          write(Hdr,file = FILE,append=T)      
          write.table(nn,file = FILE,row.names=F,col.names=F,append=T)
        }else
        {
          Hdr=paste("#",names(n)[k],sep='')
          write(Hdr,file = FILE,append=T)
          write(n[[k]],file = FILE,sep = "\t",append=T)
        }
      }
      #write("#eof",file = FILE,append=T)
      #write(999,file = FILE,sep = "\t",append=T)
      
      #copy .tpl
      file.copy(paste(A,paste("/",Spec,".tpl",sep=""),sep=""), paste(Spec,".tpl",sep=""),
                overwrite =T)
      
      #run .tpl
      args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep=""),sep=""), sep="")
      #args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep="")," -est",sep=""), sep="")
      clean_admb(Spec)
      compile_admb(Spec,verbose=T)
      system(args)
      NoCatch.future=reptoRlist(paste(Spec,".rep",sep=""))  
      NoCatch.future.std=read.table(paste(Spec,".std",sep=""), header =TRUE)
      
      
      #export figures
      fn.showX=function(Y,YLB,REL,Y.rel)
      {
        y=NoCatch.future[[match(Y,names(NoCatch.future))]]
        if(is.matrix(y)) y=rowSums(y)
        if(REL=="YES")y=y/NoCatch.future.std$value[match(Y.rel,NoCatch.future.std$name)]
        x<<-1:length(y)
        plot(x,y,type='l',ylim=c(0,max(y)),cex.lab=1.5,xlab="",ylab=YLB,lwd=2,xaxt='n',cex.axis=1.25)
        axis(1,x,F,tck=-0.02)
        axis(1,seq(x[1],x[length(x)],10),F,tck=-0.04)
      }
      
      fn.fig("Zer.catch.projections",2400, 2400)
      par(mfcol=c(2,1),mai=c(.4,.6,.05,.1),oma=c(2,2,.01,.1),las=1,mgp=c(1.6,.6,0))
      
      #Total catch
      fn.showX(Y="TC_out",YLB="",REL="NO",Y.rel="")
      mtext("Total catch (tonnes)",2,line=3.5,las=3,cex=2.1)
      
      #Spawning biomass
      fn.showX(Y="Fem_spawn_biom",YLB="",REL="YES",Y.rel="Virgin_Spawn_Biom")
      #legend("right","Total biomass",cex=1.75,bty="n")
      Yrss=as.numeric(substr(Frst.yr.ktch,1,4))
      Yrss=seq(Yrss,(Yrss+length(x)))
      Yrss=paste(Yrss,substr((Yrss+1),3,4),sep="-")
      axis(1,seq(x[1],x[length(x)],10),Yrss[seq(x[1],x[length(x)],10)],tck=-0.04,cex.axis=1.25)
      rm(x)
      mtext("Relative biomass",2,line=3.5,las=3,cex=2.1)
      mtext("Financial year",1,line=0.75,outer=T,cex=2.1)
      dev.off()
      
      exprt.pdf="NO"
      if(exprt.pdf=="YES")
      {
        #create report
        pdf("REPORT.0.catch.pdf") 
        
        #Spawning biomass per recruit
        SPR=NoCatch.future$Fem_spawn_biom/NoCatch.future$Annual_rec
        plot(SPR,type='l',lwd=2)
        points(1,NoCatch.future$Spawn_biom_per_recruit[2],pch=19,col=2,cex=2)
        mtext("Spawning biomass per recruit",3,cex=1.25)
        
        #CPUE
        plot(NoCatch.future$CPUE,pch=19,cex=2,ylab="CPUE",xlab="Year")
        lines(NoCatch.future$Est_CPUE,col=3,lwd=2)
        
        #Numbers at age from per recruit
        plot(NoCatch.future$PerRec_Exp_Age_Comp[1,],type='l',lwd=2,col=2,xlab="Age",
             main=("Per recruit numbers at age"),ylab="Numbers")
        #lines(NoCatch.future$PerRec_Exp_Age_Comp[2,],lwd=2,col='blue')
        
        
        #Size composition at equilibrium and initial
        #   plot(NoCatch.future$N_L[1,],type="l",lwd=2,xlab="Size clase",ylab="Frequency",
        #        main=("Per recruit size distribution"))
        #   lines(NoCatch.future$N_L[2,],col=2)
        #   legend("topleft",c("N_L(0)","N_L(1)"),bty="n",lty=1,col=1:2)
        
        #Time series of predictions
        fn.X=function(Y)
        {
          y=NoCatch.future[[match(Y,names(NoCatch.future))]]
          plot(1:length(y),y,type='l',ylab=Y,ylim=c(0,max(y)),cex.lab=1.5,xlab="",lwd=2)
        }
        
        par(mfcol=c(3,2),mai=c(.35,.4,.15,.1),mgp=c(1.6,.5,0))
        
        plot(NoCatch.future$TC_out,ylab="total catch",pch=19,cex.lab=1.5,xlab="")
        lines(NoCatch.future$Est_TC_out,col=3,lwd=2)
        mtext("total catch",3,cex=1.25)
        
        #   a=NoCatch.future$"Annual_Z(1)"
        #   YRSS=1:nrow(a)
        #   colfunc <- colorRampPalette(c("yellow", "green","blue","red"))
        #   CLSS=colfunc(length(YRSS))  
        #   plot(a[1,],type='l',ylab="Z",ylim=c(.2,.5),xlab="",
        #        cex.axis=1.25,cex.lab=2,col=CLSS[1],xaxt="n")
        #   for(n in 2:nrow(a)) lines(jitter(a[n,],1),col=CLSS[n])
        #   axis(1,1:length(TL.bin.mid),F,tck=-0.01)
        #   axis(1,seq(1,length(TL.bin.mid),2),TL.bin.mid[seq(1,length(TL.bin.mid),2)],las=2,tck=-0.02)
        #   mtext("Female total mortality at size",3,cex=1.25)   
        #   nn=round(length(YRSS)/5)
        #   nn2=nn+nn
        #   nn3=nn+nn+nn
        #   nn4=nn+nn+nn+nn
        #   NN=length(YRSS)
        #   legend("left",paste(YRSS[1:nn]),bty='n',lty=1,col=CLSS[1:nn],cex=0.4,xjust=0)
        #   legend("left",paste(YRSS[(nn+1):nn2]),bty='n',lty=1,col=CLSS[(nn+1):nn2],cex=0.4,inset=0.1)
        #   legend("left",paste(YRSS[(nn2+1):nn3]),bty='n',lty=1,col=CLSS[(nn2+1):nn3],cex=0.4,inset=0.2)
        #   legend("left",paste(YRSS[(nn3+1):nn4]),bty='n',lty=1,col=CLSS[(nn3+1):nn4],cex=0.4,inset=0.3)
        #   legend("left",paste(YRSS[(nn4+1):NN]),bty='n',lty=1,col=CLSS[(nn4+1):NN],cex=0.4,inset=0.4)
        #   
        #   fn.X("Annual_F")
        #   mtext("fishing mortality",3,cex=1.25)
        
        plot(rowSums(NoCatch.future$"Num_in_len_class(1)"),type='l',col="pink",xlab="",ylab="Numbers",cex.lab=1.5,lwd=2)
        lines(rowSums(NoCatch.future$"Num_in_len_class(2)"),col="blue",lwd=2)
        mtext("Population numbers (females and males)",3,cex=1)
        
        
        fn.X("Annual_rec")
        mtext("Annual recruitment",3,cex=1.25)
        
        fn.X("Total_biom")
        mtext("Total biomass",3,cex=1.25)
        
        fn.X("Fem_spawn_biom")
        mtext("Spawning biomass",3,cex=1.25)
        legend("bottomright",paste("virgin_spawning_B=",NoCatch.future$Virgin_Spawn_Biom),bty='n')
        
        with(NoCatch.future,plot(Fem_spawn_biom,Annual_rec,ylab="Recruits",pch=19,cex.lab=1.5,
                                 xlab="Spawning stock",xlim=c(0,max(Fem_spawn_biom)),ylim=c(0,max(Annual_rec))))
        mtext("Stock-recruitment",3,cex=1.25)
        dev.off()
        
      }
    }
    for(tt in 1:length(zero.Ktch.this)) fn.check.model(MODEL=zero.Ktch.this[tt],yrs.projections=zero.Ktch.yrs)  
  }
}
                   
#3 Simulation testing                      
  #note: is the model returning the initial parameter values?
if(First.run=="YES")
{
  if(Do.sim.test=="YES")
  {
    fn.source("Simulation.testing.R")
    for(tt in 1:length(Sim.Test.this))
    {
      MODEL=Sim.Test.this[tt]
      
      #run simulations
      Sim.parS=Simul.Eval(MODEL,Spec,PARS=names(Pin.pars[[match(MODEL,names(Pin.pars))]]),
                          n.sims=N.sim.test,Jitr=CPUE.jitr)
      
      #Plot parameters
      #ID=match("Dummy",colnames(Sim.parS))
      ID=which(Par.phases[[match(MODEL,names(Par.phases))]]<0)
      PINpar=Pin.pars[[match(MODEL,names(Pin.pars))]]
      if(length(ID)>0)
      {
        Sim.parS=Sim.parS[,-ID] #remove dummy par
        PINpar=PINpar[-ID]
      }
      fn.fig("boxplot",2400, 1200)
      par(las=1,mai=c(.5,.6,.1,.1),mgp=c(2,.7,0))
      a=Sim.parS-matrix(rep(PINpar,N.sim.test),nrow=N.sim.test,byrow=T)
      boxplot(a,ylim=c(quantile(a,prob=.016), quantile(a,prob=.9925)),cex.axis=.8,
              ylab="Estimated value - initial value")
      dev.off()
    }
    
  }
}

#4. Derive MSY quantities
if(Do.MSY=="YES")
{
  library(r4ss)
  library(mvtnorm)      #for multivariate normal pdf
  
  fn.MSY=function(MODEL,F.mort,yrs.projections,sdlog.rec,n.SIM)   
  {
    #set up files and directories
    setPath(Scenarios[match(MODEL,Scenarios$Model),]$Model)  
    A=getwd()
    if(!file.exists(file.path(getwd(), "MSY"))) dir.create(file.path(getwd(), "MSY")) 
    setwd(file.path(getwd(), "MSY"))
    
    #copy .tpl
    if(!file.exists(paste(Spec,"tpl",sep="")))file.copy(paste(A,paste("/",Spec,".tpl",sep=""),sep=""), 
                                paste(Spec,".tpl",sep=""), overwrite =T)
    
    ouT=vector('list',n.SIM)
    for(m in 1:n.SIM)
    {
      #1. create new pin file (use multivariate sample from MLE)
      MLE=read.admbFit(paste(A,paste("/",Spec,sep=""),sep=""))
      n.mle=1:MLE$nopar
      Nms=MLE$names[n.mle]
      Rand.par=c(rmvnorm(1,mean=MLE$est[n.mle],sigma=MLE$cov[n.mle,n.mle]))
      names(Rand.par)=Nms
      Pin.pars=scan(paste(A,paste("/",Spec,".pin",sep=""),sep=""),what = list("", "", ""))
      Names.pin.pars=Pin.pars[[2]][-1]
      Val.pin.pars=as.numeric(Pin.pars[[3]][-1])
      names(Val.pin.pars)=Names.pin.pars
      id.estim.pars=match(Nms,Names.pin.pars)  
      Val.pin.pars[id.estim.pars]=Rand.par
      par.nms=names(Val.pin.pars)
      FILE=paste(Spec,".pin",sep="")
      write("# Input parameters",file = FILE)
      for(k in 1:length(Val.pin.pars))
      {
        Hdr=paste("#",par.nms[k])
        write(Hdr,file = FILE,append=T)
        write(Val.pin.pars[k],file = FILE,sep = "\t",append=T)
      }
      
      #create dat file and run model under each assumed F value
      dummy1=vector('list',length(F.mort))
      for(f in 1:length(F.mort))
      {
        #2. create new .dat file  
        n=Inputs[[match(MODEL,names(Inputs))]]
        
        #add future F and rec deviations
        n$Fishing.mort=rep(F.mort[f],yrs.projections)
        n$yrs_Fishing.mort=yrs.projections
        n$Rec.error_MSY=rlnorm(yrs.projections, meanlog = log(1), sdlog = sdlog.rec) 
        
        #don't estimate
        n$Phases=-abs(n$Phases)
        n$Phases[1]=1
        
        #turn on switch for calculating MSY
        n$Calc_MSY=1
        
        #export new .dat file
        FILE=paste(Spec,".dat",sep="")
        nzones=n$nzone
        ModDims=unlist(c(yr.start,yr.end,nzones))
        Hdr="#Basic model dimensions (yr.start, yr.end, nzones)"
        write(Hdr,file = FILE)
        write(ModDims,file = FILE,sep = "\t",append=T)
        for(k in (length(ModDims)+1):length(n))
        {
          nn=n[[k]]
          if(is.data.frame(nn)|is.matrix(nn))
          {
            Hdr=paste("#",paste(c(names(n)[k],"(",names(nn),")"),collapse=' '))
            write(Hdr,file = FILE,append=T)      
            write.table(nn,file = FILE,row.names=F,col.names=F,append=T)
          }else
          {
            Hdr=paste("#",names(n)[k],sep='')
            write(Hdr,file = FILE,append=T)
            write(n[[k]],file = FILE,sep = "\t",append=T)
          }
        }
        
        #4. run .tpl
        args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep="")," -est",sep=""), sep="")
        if(!file.exists(paste(Spec,".exe",sep="")))compile_admb(Spec,verbose=T)
        suppressWarnings({Report=system(args, intern = TRUE)})
        
        
        #5. Extract equilibrium Fmsy, Bmsy, MSY
        B_virgin=as.numeric(Report[match("Virgin_Total_biom",Report)+1])
        B_virgin_spawning=as.numeric(Report[match("Virgin_Spawn_Biom",Report)+1])
        MSY=as.numeric(Report[match("Total_ktch_MSY",Report)+yrs.projections])
        Bmsy_total=as.numeric(Report[match("Total_biom_MSY",Report)+yrs.projections])
        Bmsy_spawning=as.numeric(Report[match("Female_Spawning_biom_MSY",Report)+yrs.projections])
        dummy1[[f]]=data.frame(Fishing=F.mort[f],MSY=MSY,
                               B_virgin=B_virgin,Bmsy_total=Bmsy_total,
                               B_virgin_spawning=B_virgin_spawning,Bmsy_spawning=Bmsy_spawning)
      }
      
      res=do.call(rbind,dummy1)
      
      #coarse F sequence
      # fit=smooth.spline(res$Fishing,res$MSY)
      # pred.prime <- data.frame(predict(fit, deriv=1))
      # id.fmsy=which.min(abs(pred.prime$y))
      # MSYs=data.frame(Fmsy = pred.prime[id.fmsy, c('x')],
      #                 BMSY_total=res$Bmsy_total[id.fmsy],
      #                 B_virgin=res$B_virgin[id.fmsy],
      #                 BMSY_spawning=res$Bmsy_spawning[id.fmsy],
      #                 B_virgin_spawning=res$B_virgin_spawning[id.fmsy],
      #                 MSY=res$MSY[id.fmsy])
      
      #finer F sequence
      F.fine.scale=seq(0,1,.01)
      fit=smooth.spline(res$Fishing,res$MSY)
      pred.prime <- data.frame(predict(fit,x=F.fine.scale, deriv=1))
      id.fmsy=which.min(abs(pred.prime$y))
      fit.msy <- data.frame(predict(fit,x=F.fine.scale))
      fit.bmsy=smooth.spline(res$Fishing,res$Bmsy_total)
      fit.bmsy <- data.frame(predict(fit.bmsy,x=F.fine.scale))
      fit.B.virgin=smooth.spline(res$Fishing,res$B_virgin)
      fit.B.virgin <- data.frame(predict(fit.B.virgin,x=F.fine.scale))
      fit.bmsy.spawn=smooth.spline(res$Fishing,res$Bmsy_spawning)
      fit.bmsy.spawn <- data.frame(predict(fit.bmsy.spawn,x=F.fine.scale))
      fit.B.virgin.spawn=smooth.spline(res$Fishing,res$B_virgin_spawning)
      fit.B.virgin.spawn <- data.frame(predict(fit.B.virgin.spawn,x=F.fine.scale))
      MSYs=data.frame(Fmsy = pred.prime[id.fmsy, c('x')],
                      BMSY_total=fit.bmsy[id.fmsy, c('y')],
                      B_virgin=fit.B.virgin[id.fmsy, c('y')],
                      BMSY_spawning=fit.bmsy.spawn[id.fmsy, c('y')],
                      B_virgin_spawning=fit.B.virgin.spawn[id.fmsy, c('y')],
                      MSY=fit.msy[id.fmsy, c('y')])
      
      ouT[[m]]=list(res=res,MSYs=MSYs)
    }
    return(ouT)
  }
    Store.MSY=fn.MSY(MODEL="Base case",F.mort=F.vec,yrs.projections=MSY.yrs,
                     sdlog.rec=MSY.sd.rec,n.SIM=MSY.sims)

  
  #Export results
  library(purrr)
  library(reshape2)
  library(plotrix)
  
  RES=map(Store.MSY, 1)
  d=array(as.numeric(unlist(RES)), dim=c(nrow(RES[[1]]),ncol(RES[[1]]), length(RES)))
  RES_stats <- apply(d, c(1,2), quantile, probs=c(0.025,.5,.975), na.rm=TRUE)
  RES_stats=dcast(melt(RES_stats), Var1+Var2~Var3) 
  colnames(RES_stats)=c('Percentile','index',names(RES[[1]]))
  
  
  MSYs=map(Store.MSY, 2)
  d=array(as.numeric(unlist(MSYs)), dim=c(nrow(MSYs[[1]]),ncol(MSYs[[1]]), length(MSYs)))
  MSY_stats <- apply(d, c(1,2), quantile, probs=c(0.025,.5,.975), na.rm=TRUE)
  MSY_stats=dcast(melt(MSY_stats), Var1~Var3) 
  colnames(MSY_stats)=c('Percentile',names(MSYs[[1]]))
  idd=match(c('Percentile','Fmsy'),names(MSY_stats))
  MSY_stats[,-idd]=round( MSY_stats[,-idd])
  
  inx=length(F.vec)
  fn.fig("MSY_quantities",2400, 2400)
  par(mfcol=c(1,1),mar=c(3,3,.1,.1),oma=c(1,1,.01,.1),las=1,mgp=c(1.95,.6,0),cex.lab=1.25)
  
  plot(RES_stats$Fishing[1:inx],RES_stats$MSY[(inx+1):(2*inx)],ylab="Equilibrium catch (tonnes)",
       xlab='Equilibrium fishing mortality',type='l',lwd=2,ylim=c(0,max(RES_stats$MSY)))
  polygon(x=c(RES_stats$Fishing[1:inx],rev(RES_stats$Fishing[1:inx])),
          y=c(RES_stats$MSY[1:inx],rev(RES_stats$MSY[(2*inx+1):(3*inx)])),
          col=rgb(.1,.4,.2,alpha=.3))
  polygon(x=c(MSY_stats$Fmsy[1],MSY_stats$Fmsy[3],MSY_stats$Fmsy[3],MSY_stats$Fmsy[1]),
          y=c(-20,-20,max(RES_stats$MSY),max(RES_stats$MSY)),col=rgb(.3,.2,.1,alpha=.3))
  
  addtable2plot(0 ,1,MSY_stats,bty="o",display.rownames=F,hlines=F,
                vlines=TRUE,title="",cex=1,xjust=-0.05)
  
  dev.off()
  
  fn.fig("Relative spawning biomass",2400, 2400)
  x=RES_stats$Fishing[1:inx]
  y=RES_stats$Bmsy_spawning[(inx+1):(2*inx)]/RES_stats$B_virgin_spawning[(inx+1):(2*inx)]  
  plot(x,y,type='l',lwd=2,col="steelblue",ylab="Relative spawning biomass",
       xlab='Equilibrium fishing mortality',ylim=c(0,1.05),yaxs = "i",xaxs = "i")
  #Rel.b=MSY_stats$BMSY_spawning[2]/MSY_stats$B_virgin_spawning[2]
  #iidd=which.min(abs(y-Rel.b))
  #lines(c(MSY_stats$Fmsy[2],MSY_stats$Fmsy[2]),c(0,y[iidd]),lwd=2,col='firebrick')
  #lines(c(MSY_stats$Fmsy[2],0),c(y[iidd],y[iidd]),lwd=2,col='firebrick')
  dev.off()
}

#Profile likelikhood for Fo (run with args -lprof)
# note: ADMB cannot calculate the profile for S3
#dd=read.table("lp_ln_In.plt", skip=2, nrow=70)
#par(mfcol=c(1, 1), cex.lab=1, cex.axis=1)
#plot(dd, xlab="Estimate", ylab="Density", type="l") 


# Section E: DISPLAY MODEL OUTPUTS ----------------------------------------
fn.source("Pearson.Residuals.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R")
library(plotrix)
library(gridExtra)
library(grid)

yr.start=as.numeric(substr(Frst.yr.ktch,1,4))
yr.end=as.numeric(substr(Data.yr,1,4))

CL="forestgreen"
CLzone=2:4

#2. Some useful functions
fn.lay.multi=function()
{
  par(mfcol=c(4,3),mai=c(.55,.6,.05,.05),las=1)
}

fn.lay.single=function(MAI) par(mfcol=c(1,1),mai=MAI,las=1,mgp=c(2.5, .7, 0))

fn.get.CI=function(Coef,SE,LOGGED)
{
  if(LOGGED=="YES")
    {
     LOW=exp(Coef-(1.96*SE))
     UP=exp(Coef+(1.96*SE))
     MLE=exp(Coef)
   }else
   {
     LOW=(Coef-(1.96*SE))
     UP=(Coef+(1.96*SE))
     MLE=Coef
   }
  return(list(low=LOW,MLE=MLE,Up=UP))
}

fn.comp.ob.pred=function(what,what1,where.leg) 
{
  OBS=MOD[[match(what,names(MOD))]]
  PRED=MOD[[match(what1,names(MOD))]]
  id=match("CPUE_SD_out",names(MOD))
  YLIM=c(0,max(c(OBS,PRED)))
  if(what=="CPUE_out" & !is.na(id)) YLIM=c(0,max(c(PRED,OBS+MOD$CPUE_SD_out)))
  Id=which(OBS<0)
  OBS[Id]=PRED[Id]=NA
  plot(Yrs,OBS[1:length(Yrs)],ylab="",xlab="",pch=19,cex.axis=1.5,xaxt='n',cex=1.75,ylim=YLIM)
  if(what=="CPUE_out")
  {
    ii=match(Yr_q_change,Yrs)
    iii=match(Yr_q_daily,Yrs)
    if(Yr_q_change>0)
    {
      lines(Yrs[1:ii],PRED[1:ii],col=CL,lwd=3,lty=LTYp[1])
      lines(Yrs[(ii+1):(iii-1)],PRED[(ii+1):(iii-1)],col=CL,lwd=3,lty=LTYp[2])
    }
    if(Yr_q_change==0)lines(Yrs[1:(iii-1)],PRED[1:(iii-1)],col=CL,lwd=3,lty=LTYp[1])
    lines(Yrs[iii:length(Yrs)],PRED[iii:length(Yrs)],col=CL,lwd=3,lty=LTYp[3])
  }else lines(Yrs,PRED[1:length(Yrs)],col=CL,lwd=3)
  if(what%in%c("CPUE_out","CPUE_eff"))
  {
    if(!is.na(match("CPUE_SD_out",names(MOD))))legend(where.leg,c("observed (? CV)","predicted"),bty='n',pch=c(19,NA),cex=2,lty=c(NA,1),col=c("black",CL),lwd=3)else
      legend(where.leg,c("observed","predicted"),bty='n',pch=c(19,NA),cex=2,lty=c(NA,1),col=c("black",CL),lwd=3)
  } else
  {
    legend(where.leg,c("observed","predicted"),bty='n',pch=c(19,NA),cex=2,lty=c(NA,1),col=c("black",CL),lwd=3)
  }
  axis(1,Yrs,F,tck=-0.015)
  axis(1,seq(Yrs[1],Yrs[length(Yrs)],5),seq(Yrs[1],Yrs[length(Yrs)],5),tck=-0.03,cex.axis=1.25)
  if(what=="CPUE_out"  & !is.na(id))
  {
    segments(Yrs,OBS,Yrs,(OBS+MOD$CPUE_SD_out),lwd=2)
    segments(Yrs,OBS,Yrs,(OBS-MOD$CPUE_SD_out),lwd=2)
  }
}

fn.poly=function(XX,YY,CL) polygon(XX,YY,col=CL,border="transparent")

comp.fn.comp.ob.pred=function(what,what1,LW)
{
  OBS=MOD[[match(what,names(MOD))]]
  PRED=MOD[[match(what1,names(MOD))]]
  id=match("CPUE_SD_out",names(MOD))
  YLIM=c(0,max(c(OBS,PRED)))
  YY=c(YLIM[1],YLIM[1],YLIM[2]*1.058,YLIM[2]*1.058)   
  if(what=="CPUE_out" & !is.na(id)) YLIM=c(0,max(c(PRED,OBS+MOD$CPUE_SD_out)))
  Id=which(OBS<0)
  OBS[Id]=PRED[Id]=NA
  plot(Yrs,OBS[1:length(Yrs)],ylab="",xlab="",pch=21,bg='grey80',col="grey40",     
       cex.axis=1.5,xaxt='n',cex=1.25,ylim=YLIM)
  id.cv=match("CPUE_SD_out",names(MOD))
  if(what=="CPUE_out"  & !is.na(id.cv))
  {
    segments(Yrs,OBS,Yrs,(OBS+MOD$CPUE_SD_out),lwd=LW,col="grey40")
    segments(Yrs,OBS,Yrs,(OBS-MOD$CPUE_SD_out),lwd=LW,col="grey40")
    points(Yrs,OBS[1:length(Yrs)],pch=21,bg='grey80',col="grey40",cex=1.25)
  }
  if(Scenarios$Q[i]=="one")
  {
    XX=Yrs     
    lines(XX,PRED,col=CL,lwd=3,lty=LTYp[1])
    XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
    fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
    
  }else
  {
    if(what=="CPUE_out"){
      ii=match(Yr_q_change,Yrs)
      iii=match(Yr_q_daily,Yrs)
      if(Yr_q_change>0)
      {
        Nms.par.fz=names(Fze)
        id.q2=c("lnq2","ln_q2","q2")
        id.q2=id.q2[which(id.q2%in%Nms.par.fz)]
        if(Fze[match(id.q2,names(Fze))]>0)
        {
          XX=Yrs[1:ii]    
          lines(XX,PRED[1:ii],col=CL,lwd=LW,lty=LTYp[1])
          XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
          fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
          
          
          id.qdaily=c("log_Qdaily","ln_qdaily","qdaily")
          id.qdaily=id.qdaily[which(id.qdaily%in%Nms.par.fz)]
          
          if(Fze[match(id.qdaily,names(Fze))]>0)
          {
            XX=Yrs[(ii+1):(iii-1)]
            lines(XX,PRED[(ii+1):(iii-1)],col="firebrick4",lwd=LW,lty=LTYp[2])
            XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
            fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
            
            XX=Yrs[iii:length(Yrs)]
            lines(XX,PRED[iii:length(Yrs)],col="royalblue4",lwd=LW,lty=LTYp[3])
            XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)  
            fn.poly(XX=XX,YY=YY,CL=rgb(.1,.1,.4,alpha=.15))
            
            
          }else
          {
            XX=Yrs[(ii+1):length(Yrs)]
            lines(XX,PRED[(ii+1):length(Yrs)],col="firebrick4",lwd=LW,lty=LTYp[2])
            XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
            fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
          }
        }else
        {
          XX=Yrs[1:(iii-1)]     
          lines(XX,PRED[1:(iii-1)],col=CL,lwd=LW,lty=LTYp[1]) 
          XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
          fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
          
          XX=Yrs[iii:length(Yrs)]     
          lines(XX,PRED[iii:length(Yrs)],col="firebrick4",lwd=LW,lty=LTYp[2]) 
          XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)  
          fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
          
        }
        
        
      }
      if(Yr_q_change==0)
      {
        XX=Yrs[1:(iii-1)]     
        lines(XX,PRED[1:(iii-1)],col=CL,lwd=LW,lty=LTYp[1])
        XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
        fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
        
        XX=Yrs[iii:length(Yrs)]     
        lines(XX,PRED[iii:length(Yrs)],col="firebrick4",lwd=LW,lty=LTYp[2])
        XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
        fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
        
      }
      
    }else if(what=="CPUE_eff"){
      ii=match(Yr_q_change,Yrs)
      if(Yr_q_change>0)
      {
        XX=Yrs[1:ii]     
        lines(XX,PRED[1:ii],col=CL,lwd=3,lty=LTYp[1])
        XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
        fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
        
        XX=Yrs[(ii+1):length(Yrs)]     
        lines(XX,PRED[(ii+1):length(Yrs)],col="firebrick4",lwd=LW,lty=LTYp[2])
        XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
        fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
        
      }
    }else
    {
      XX=Yrs     
      lines(XX,PRED[1:length(Yrs)],col=CL,lwd=LW)
      XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
      fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
    }
  }
  axis(1,Yrs,F,tck=-0.015)
  axis(1,seq(Yrs[1],Yrs[length(Yrs)],5),seq(Yrs[1],Yrs[length(Yrs)],5),tck=-0.03,cex.axis=1.5)
}
comp.fn.comp.ob.pred_logcpue=function(what,what1,LW)
{
  OBS=MOD[[match(what,names(MOD))]]
  OBS[OBS<0]=NA
  OBS=log(OBS)
  PRED=subset(STD,name==what1)
  id=match("CPUE_SD_out",names(MOD))
  YLIM=c(min(c(OBS,PRED$value-1.96*PRED$std.dev),na.rm=T),max(c(OBS,PRED$value+1.96*PRED$std.dev),na.rm=T))
  YY=c(YLIM[1],YLIM[1],YLIM[2]*1.058,YLIM[2]*1.058)   
  if(what=="CPUE_out" & !is.na(id))
  {
    OBS_CV=MOD$CPUE_SD_out
    YLIM=c(min(c(OBS-OBS_CV,PRED$value-1.96*PRED$std.dev),na.rm=T),max(c(OBS+OBS_CV,PRED$value+1.96*PRED$std.dev),na.rm=T))
    YY=c(YLIM[1],YLIM[1],YLIM[2]*1.058,YLIM[2]*1.058)
  }
  
  Id=which(is.na(OBS))
  PRED[Id,]=NA
  plot(Yrs,OBS[1:length(Yrs)],ylab="",xlab="",pch=21,bg='grey80',col="grey40",     
       cex.axis=1.5,xaxt='n',cex=1.25,ylim=YLIM)
  id.cv=match("CPUE_SD_out",names(MOD))
  if(what=="CPUE_out"  & !is.na(id.cv))
  {
    segments(Yrs,(OBS-OBS_CV),Yrs,(OBS+OBS_CV),lwd=LW,col="grey40")
    points(Yrs,OBS[1:length(Yrs)],pch=21,bg='grey80',col="grey40",cex=1.25)
  }
  if(Scenarios$Q[i]=="one")
  {
    XX=Yrs     
    lines(XX,PRED,col=CL,lwd=3,lty=LTYp[1])
    XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
    fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
    
  }else
  {
    if(what=="CPUE_out"){
      ii=match(Yr_q_change,Yrs)
      iii=match(Yr_q_daily,Yrs)
      if(Yr_q_change>0)
      {
        Nms.par.fz=names(Fze)
        id.q2=c("lnq2","ln_q2","q2")
        id.q2=id.q2[which(id.q2%in%Nms.par.fz)]
        if(Fze[match(id.q2,names(Fze))]>0)
        {
          XX=Yrs[1:ii]    
          points(XX,PRED$value[1:ii] ,col=CL,pch=19)
          segments(XX,PRED$value[1:ii]-1.96*PRED$std.dev[1:ii],
                   XX,PRED$value[1:ii]+1.96*PRED$std.dev[1:ii],lwd=LW,col=CL)
          XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
          fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
          id.qdaily=c("log_Qdaily","ln_qdaily","qdaily")
          id.qdaily=id.qdaily[which(id.qdaily%in%Nms.par.fz)]
          
          if(Fze[match(id.qdaily,names(Fze))]>0)
          {
            XX=Yrs[(ii+1):(iii-1)]
            points(XX,PRED$value[(ii+1):(iii-1)] ,col="firebrick4",pch=19)
            segments(XX,PRED$value[(ii+1):(iii-1)]-1.96*PRED$std.dev[(ii+1):(iii-1)],
                     XX,PRED$value[(ii+1):(iii-1)]+1.96*PRED$std.dev[(ii+1):(iii-1)],lwd=LW,col="firebrick4")
            XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
            fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
            
            XX=Yrs[iii:length(Yrs)]
            points(XX,PRED$value[iii:length(Yrs)] ,col="royalblue4",pch=19)
            segments(XX,PRED$value[iii:length(Yrs)]-1.96*PRED$std.dev[iii:length(Yrs)],
                     XX,PRED$value[iii:length(Yrs)]+1.96*PRED$std.dev[iii:length(Yrs)],lwd=LW,col="royalblue4")
            XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)  
            fn.poly(XX=XX,YY=YY,CL=rgb(.1,.1,.4,alpha=.15))
            
            
          }else
          {
            XX=Yrs[(ii+1):length(Yrs)]
            points(XX,PRED$value[(ii+1):length(Yrs)] ,col="firebrick4",pch=19)
            segments(XX,PRED$value[(ii+1):length(Yrs)]-1.96*PRED$std.dev[(ii+1):length(Yrs)],
                     XX,PRED$value[(ii+1):length(Yrs)]+1.96*PRED$std.dev[(ii+1):length(Yrs)],lwd=LW,col="firebrick4")
            
            XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
            fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
          }
        }else
        {
          XX=Yrs[1:(iii-1)]  
          points(XX,PRED$value[1:(iii-1)] ,col=CL,pch=19)
          segments(XX,PRED$value[1:(iii-1)]-1.96*PRED$std.dev[1:(iii-1)],
                   XX,PRED$value[1:(iii-1)]+1.96*PRED$std.dev[1:(iii-1)],lwd=LW,col=CL)
          XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
          fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
          
          XX=Yrs[iii:length(Yrs)]   
          points(XX,PRED$value[iii:length(Yrs)] ,col="firebrick4",pch=19)
          segments(XX,PRED$value[iii:length(Yrs)]-1.96*PRED$std.dev[iii:length(Yrs)],
                   XX,PRED$value[iii:length(Yrs)]+1.96*PRED$std.dev[iii:length(Yrs)],lwd=LW,col="firebrick4")
          XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)  
          fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
        }
      }
      if(Yr_q_change==0)
      {
        XX=Yrs[1:(iii-1)] 
        points(XX,PRED$value[1:(iii-1)] ,col=CL,pch=19)
        segments(XX,PRED$value[1:(iii-1)]-1.96*PRED$std.dev[1:(iii-1)],
                 XX,PRED$value[1:(iii-1)]+1.96*PRED$std.dev[1:(iii-1)],lwd=LW,col=CL)
        XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
        fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
        
        XX=Yrs[iii:length(Yrs)] 
        points(XX,PRED$value[iii:length(Yrs)] ,col="firebrick4",pch=19)
        segments(XX,PRED$value[iii:length(Yrs)]-1.96*PRED$std.dev[iii:length(Yrs)],
                 XX,PRED$value[iii:length(Yrs)]+1.96*PRED$std.dev[iii:length(Yrs)],lwd=LW,col="firebrick4")
        XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
        fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
        
      }
      
    }else if(what=="CPUE_eff"){
      ii=match(Yr_q_change,Yrs)
      if(Yr_q_change>0)
      {
        XX=Yrs[1:ii]  
        points(XX,PRED$value[1:ii] ,col=CL,pch=19)
        segments(XX,PRED$value[1:ii]-1.96*PRED$std.dev[1:ii],
                 XX,PRED$value[1:ii]+1.96*PRED$std.dev[1:ii],lwd=LW,col=CL)
        XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
        fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
        
        XX=Yrs[(ii+1):length(Yrs)]  
        points(XX,PRED$value[(ii+1):length(Yrs)] ,col="firebrick4",pch=19)
        segments(XX,PRED$value[(ii+1):length(Yrs)]-1.96*PRED$std.dev[(ii+1):length(Yrs)],
                 XX,PRED$value[(ii+1):length(Yrs)]+1.96*PRED$std.dev[(ii+1):length(Yrs)],lwd=LW,col="firebrick4")
        XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
        fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
        
      }
    }else
    {
      XX=Yrs   
      points(XX,PRED$value[1:length(Yrs)] ,col=CL,pch=19)
      segments(XX,PRED$value[1:length(Yrs)]-1.96*PRED$std.dev[1:length(Yrs)],
               XX,PRED$value[1:length(Yrs)]+1.96*PRED$std.dev[1:length(Yrs)],lwd=LW,col=CL)
      XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
      fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
    }
  }
  axis(1,Yrs,F,tck=-0.015)
  axis(1,seq(Yrs[1],Yrs[length(Yrs)],5),seq(Yrs[1],Yrs[length(Yrs)],5),tck=-0.03,cex.axis=1.5)
}

fn.comp.ob.pred.spatial=function(what,what1,where.leg,Plot.dims)   
{
  OBS=MOD[[match(what,names(MOD))]]
  PRED=MOD[[match(what1,names(MOD))]]
  peep=which(OBS<0)
  PRED[peep]=NA
  OBS[OBS<0]=NA 
  id=match("CPUE_SD_out",names(MOD))
  if(Plot.dims=="YES") par(mfcol=c(ncol(OBS),1),mai=c(.75,.75,.1,.1),las=1)
  for(t in 1:ncol(OBS))
  {
    YLIM=c(0,max(c(OBS[1:length(Yrs),t],PRED[1:length(Yrs),t]),na.rm=T))
    if(what=="CPUE_out" & !is.na(id)) YLIM=c(0,max(c(OBS[1:length(Yrs),t]+MOD$CPUE_SD_out[1:length(Yrs),t],PRED[1:length(Yrs),t]),na.rm=T))
    YY=c(YLIM[1],YLIM[1],YLIM[2]*1.058,YLIM[2]*1.058)   
    Id=which(OBS<0)
    OBS[Id,t]=PRED[Id,t]=NA
    
    plot(Yrs,OBS[1:length(Yrs),t],ylab="",xlab="",pch=21,bg='grey80',col="grey40",cex.axis=1.5,xaxt='n',cex=1.75,ylim=YLIM)
    if(what=="CPUE_out"  & !is.na(id))
    {
      segments(Yrs,OBS[1:length(Yrs),t],Yrs,(OBS[1:length(Yrs),t]+MOD$CPUE_SD_out[,t]),lwd=LW,col="grey40")
      segments(Yrs,OBS[1:length(Yrs),t],Yrs,(OBS[1:length(Yrs),t]-MOD$CPUE_SD_out[,t]),lwd=LW,col="grey40")
      points(Yrs,OBS[1:length(Yrs),t],pch=21,bg='grey80',col="grey40",cex=1.75)
    }
    if(what=="CPUE_out")
    {
      ii=match(Yr_q_change,Yrs)
      iii=match(Yr_q_daily,Yrs)
      if(Yr_q_change>0)
      {
        XX=Yrs[1:ii]     
        lines(XX,PRED[1:ii,t],col=CL,lwd=LW,lty=LTYp[1])
        XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
        if(any(!is.na(PRED[1:ii,t])))fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
        
        XX=Yrs[(ii+1):(iii-1)]
        lines(XX,PRED[(ii+1):(iii-1),t],col="firebrick4",lwd=LW,lty=LTYp[2])
        XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
        if(any(!is.na(PRED[1:ii,t])))fn.poly(XX=XX,YY=YY,CL=rgb(.4,.1,.1,alpha=.15))
      }
      
      if(Yr_q_change==0)
      {
        XX=Yrs[1:(iii-1)]
        lines(XX,PRED[1:(iii-1),t],col=CL,lwd=LW,lty=LTYp[1])
        XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
        if(any(!is.na(PRED[1:ii,t])))fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
      }
      
      XX=Yrs[iii:length(Yrs)]  
      lines(XX,PRED[iii:length(Yrs),t],col="royalblue4",lwd=LW,lty=LTYp[3])
      XX=c(XX[1]-.25,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-.25)   
      if(any(!is.na(PRED[iii:length(Yrs),t])))fn.poly(XX=XX,YY=YY,CL=rgb(.1,.1,.4,alpha=.15))
    }else
    {
      XX=Yrs
      lines(XX,PRED[1:length(Yrs),t],col=CL,lwd=LW)
      XX=c(XX[1]-1,XX[length(XX)]+.25,XX[length(XX)]+.25,XX[1]-1)
      if(any(!is.na(PRED[1:ii,t])))fn.poly(XX=XX,YY=YY,CL=rgb(.1,.4,.1,alpha=.15))
    }
      
    axis(1,Yrs,F,tck=-0.02)
    axis(1,seq(Yrs[1],Yrs[length(Yrs)],5),F,tck=-0.04,cex.axis=1.25)
    legend("top",as.character(Areas.zones$zone[t]),bty='n',cex=1.5)
    if(Plot.dims=="NO")axis(1,seq(Yrs[1],Yrs[length(Yrs)],5),seq(Yrs[1],Yrs[length(Yrs)],5),tck=-0.04,cex.axis=1.5)
  }
  if(Plot.dims=="YES")
  {
    legend(where.leg,c("observed (? CV)","predicted"),bty='n',pch=c(21,NA),cex=1.75,lty=c(NA,1),
           col=c("grey40","black"),pt.bg=c("grey80",""),lwd=3)
    axis(1,seq(Yrs[1],Yrs[length(Yrs)],5),seq(Yrs[1],Yrs[length(Yrs)],5),tck=-0.04,cex.axis=1.5)
  }
}

fn.comp.ob.pred.size.comp=function(what,what1,BINS,where.leg,N,COL)
{
  OBS=MOD[[match(what,names(MOD))]]
  PRED=MOD[[match(what1,names(MOD))]]
  bins=MOD[[match(BINS,names(MOD))]]
  Effective.n=MOD[[match(N,names(MOD))]][,COL]
  rownames(OBS)=rownames(PRED)=Yrs
  KEEP=rowSums(OBS)
  KEEP=subset(KEEP,KEEP>0)
  id=match(names(KEEP),rownames(OBS))
  OBS=OBS[id,]
  Effective.n=Effective.n[id]
  PRED=PRED[id,]
  YLIM=c(0,max(c(OBS,PRED)))  
  smart.par(n.plots=length(Effective.n)+1,MAR=c(4,4,1,1),OMA=rep(1,4),MGP=c(2.5,.7,0))
  for(n in 1:nrow(OBS))
  {
    pp <- barplot(OBS[n,],ylim=YLIM,col="grey60",border="grey60",xaxt='n',cex.axis=1.25)
    lines(x = pp, y = PRED[n,], col=CL,lwd=3)
    box()
    legend("topleft",paste("Observations=",Effective.n[n]),bty='n',cex=1.25,title=rownames(OBS)[n]) 
    axis(1,at=pp,labels=NA,line=0,tck=-0.025)
    axis(1,at=pp[seq(1,length(pp),5)],labels=NA,tck=-0.05)
    axis(1,at=pp[seq(1,length(pp),5)],labels=bins[seq(1,length(pp),5)],
         cex.axis=1.25,tck=-0.05)
  }
  plot(1:10,axes=F,col="transparent",ylab='',xlab='')
  legend(where.leg,c("observed","predicted"),bty='n',pch=c(15,NA),cex=1.5,lty=c(NA,1),
         col=c("grey60",CL),lwd=3,pt.cex=c(3,NA))
}

fn.obs.age.len=function(X,what,what1,cl)     
{
  OBS=MOD[[match(what,names(MOD))]]
  PRED=MOD[[match(what1,names(MOD))]]
  PRED=data.frame(X=X,PRED=PRED)
  PRED=PRED[order(PRED$X),]
  YLIM=c(0,1.1*max(c(OBS,PRED$PRED)))
  plot(X,OBS,ylab="",xlab="",pch=19,cex.axis=1.5,cex=2,ylim=YLIM,xlim=c(0,max(X)*1.1),col=cl, yaxs="i",xaxs="i")
  with(PRED,lines(X,PRED,col=CL,lwd=3))
  arrows(1,Lo,0,Lo,lwd=2,length = 0.1)
}

fn.see.pred=function(X,what,lista.columns)
{
  OBS=MOD[[match(what,names(MOD))]]
  YLIM=c(0,max(OBS))
  if(!class(OBS)=="matrix")  plot(X,OBS[1:length(X)],ylab="",xlab="",type='l',cex.axis=1.5,xaxt='n',cex=1.5,lwd=3,ylim=YLIM,col=CL)
  if(class(OBS)=="matrix")
  {
    if(ncol(OBS)==1)
    {
      if(lista.columns>1)
      {
        OBS=matrix(OBS,ncol=lista.columns)
        plot(X,OBS[1:length(X),1],ylab="",xlab="",type='l',cex.axis=1.5,xaxt='n',cex=1.5,lwd=3,ylim=YLIM,col=CL)
        lines(X,OBS[1:length(X),2],lty=3,lwd=2,col=CL)
      }else plot(X,OBS[1:length(X),1],ylab="",xlab="",type='l',cex.axis=1.5,xaxt='n',cex=1.5,lwd=3,ylim=YLIM,col=CL)
    }else
    {
      plot(X,OBS[1,1:length(X)],ylab="",xlab="",type='l',cex.axis=1.5,xaxt='n',cex=1.5,lwd=3,ylim=YLIM,col="pink")
      lines(X,OBS[2,1:length(X)],lwd=3,col="blue")
      
    }
  }
  
  axis(1,X,F,tck=-0.0075)
  axis(1,seq(X[1],X[length(X)],5),seq(X[1],X[length(X)],5),tck=-0.015,cex.axis=1.25)
}

fn.see.F=function(X,what,what1,deplet,Spatial)
{
  OBS=subset(STD,name==what,select=c(value, std.dev))
  OBS1=subset(STD,name==what1,select=c(value, std.dev))
  YLIM=c(0,max(c(OBS$value+1.96*OBS$std.dev,OBS1$value+1.96*OBS$std.dev)))
  if(Spatial=="Single zone")
  {
    OBS=OBS[1:length(X),]
    OBS1=OBS1[1:length(X),]
    Now=length(OBS$value)
    plot(X,OBS$value,ylab="",xlab="",cex.axis=1.5,xaxt='n',cex=1.5,pch=19,ylim=YLIM,col=CL)
    with(OBS,segments(X,value-1.96*std.dev,X,value+1.96*std.dev,col=CL))
    
    with(OBS1,segments(X,value-1.96*std.dev,X,value+1.96*std.dev,col=CL))
    points(X,OBS1$value,cex=1.5,pch=21,bg="white",col=CL)
    axis(1,X,F,tck=-0.01)
  }else
  {
    Zns=as.character(Areas.zones$zone)
    zz=length(Zns)
    par_zn()
    for(z in 1:zz)
    {
      DD=OBS[seq(z,nrow(OBS),by=zz),]
      DD=DD[1:length(X),]
      Now=length(DD$value)
      plot(X,DD$value,ylab="",xlab="",cex.axis=1.5,xaxt='n',cex=1.5,pch=19,ylim=YLIM,col=CL)
      with(DD,segments(X,value-1.96*std.dev,X,value+1.96*std.dev,col=CL))
      with(OBS1,segments(X,value-1.96*std.dev,X,value+1.96*std.dev,col=CL))
      points(X,OBS1$value,cex=1.5,pch=21,bg="white",col=CL)
      axis(1,X,F,tck=-0.01)
      legend('topright',Zns[z],bty='n',cex=1.75)
      axis(1,seq(X[1],X[length(X)],5),F,tck=-0.02)
    }
  }
  legend("topright",c("females","males"),bty='n',pch=c(19,21),col=CL,bg=c(NA,"white"),cex=1.5)
  axis(1,seq(X[1],X[length(X)],5),seq(X[1],X[length(X)],5),tck=-0.02,cex.axis=1.5)
}

fn.see.pred_F=function(X,what,lista.columns)
{
  OBS=MOD[[match(what,names(MOD))]]
  YLIM=c(0,max(OBS))
  if(!class(OBS)=="matrix")  plot(X,OBS[1:length(X)],ylab="",xlab="",type='l',cex.axis=1.5,xaxt='n',cex=1.5,lwd=3,ylim=YLIM,col=CL)
  if(class(OBS)=="matrix")
  {
    if(ncol(OBS)==1)
    {
      if(lista.columns>1)
      {
        OBS=matrix(OBS,ncol=lista.columns)
        plot(X,OBS[1:length(X),1],ylab="",xlab="",type='l',cex.axis=1.5,xaxt='n',cex=1.5,lwd=3,ylim=YLIM,col=CL)
        lines(X,OBS[1:length(X),2],lty=3,lwd=2,col=CL)
      }else
      {
        Nrw=nrow(OBS)
        Fem=OBS[1:(Nrw/2),]
        Male=OBS[(Nrw/2+1):nrow(OBS),]
        plot(X,Fem[1:length(X)],ylab="",xlab="",type='l',cex.axis=1.5,xaxt='n',cex=1.5,lwd=3,ylim=YLIM,col=CL)
        lines(X,Male[1:length(X)],lty=3,lwd=2,col=CL)
      }
      axis(1,X,F,tck=-0.01)
    }else
    {
      Nrw=nrow(OBS)
      Fem=OBS[1:(Nrw/2),]
      Male=OBS[(Nrw/2+1):nrow(OBS),]
      
      Zns=as.character(Areas.zones$zone)
      zz=length(Zns)
      par_zn()
      for(z in 1:zz)
      {
        plot(X,Fem[1:length(X),z],ylab="",xlab="",type='l',cex.axis=1.5,xaxt='n',cex=1.5,lwd=3,ylim=YLIM,col=CL)
        lines(X,Male[1:length(X),z],lty=3,lwd=2,col=CL)
        legend('topright',Zns[z],bty='n',cex=1.5)
        axis(1,X,F,tck=-0.01)
        axis(1,seq(X[1],X[length(X)],5),F,tck=-0.02)
      }
    }
  }
  axis(1,seq(X[1],X[length(X)],5),seq(X[1],X[length(X)],5),tck=-0.02,cex.axis=1.5)
}

fn.see.pred.spatial=function(X,what,where.leg,yLIm)
{
  OBS=MOD[[match(what,names(MOD))]]
  if(is.na(yLIm[1]))YLIM=c(0,max(OBS))
  if(!is.na(yLIm[1]))YLIM=yLIm
  plot(X,OBS[1:length(X),1],ylab="",xlab="",type='l',cex.axis=1.5,xaxt='n',cex=1.5,lwd=3,ylim=YLIM,col=CLzone[1])
  for(q in 2:3)lines(X,OBS[1:length(X),q],lwd=3,col=CLzone[q])  
  axis(1,X,F,tck=-0.015)
  axis(1,seq(X[1],X[length(X)],5),seq(X[1],X[length(X)],5),tck=-0.03,cex.axis=1.25)
  legend(where.leg,as.character(Areas.zones$zone),bty='n',lty=1,col=CLzone,lwd=3,cex=1.5)  
}

fn.see.Biom=function(X,what,Spatial)
{
  OBS=subset(STD,name==what,select=c(value, std.dev))
  YLIM=c(0,max(OBS$value+1.96*OBS$std.dev))
  if(Spatial=="Single zone")
  {
    OBS=OBS[1:length(X),]
    Now=length(OBS$value)
    plot(X,OBS$value,ylab="",xlab="",cex.axis=1.5,xaxt='n',cex=1.25,pch=19,ylim=YLIM,col=CL)
    with(OBS,segments(X,value-1.96*std.dev,X,value+1.96*std.dev,col=CL))
    axis(1,X,F,tck=-0.015)
  }else
  {
    Zns=as.character(Areas.zones$zone)
    zz=length(Zns)
    par_zn()
    for(z in 1:zz)
    {
      DD=OBS[seq(z,nrow(OBS),by=zz),]
      DD=DD[1:length(X),]
      Now=length(DD$value)
      plot(X,DD$value,ylab="",xlab="",cex.axis=1.5,xaxt='n',cex=1.5,pch=19,ylim=YLIM,col=CL)
      with(DD,segments(X,value-1.96*std.dev,X,value+1.96*std.dev,col=CL))
      axis(1,X,F,tck=-0.015)
      legend('topright',Zns[z],bty='n',cex=1.75)
      axis(1,seq(X[1],X[length(X)],5),F,tck=-0.03)
    }
  }
  axis(1,seq(X[1],X[length(X)],5),seq(X[1],X[length(X)],5),tck=-0.03,cex.axis=1.5)
}

fn.plt.STM=function(Size.Tran.Mat,Sex,TL.bin.mid)
{
  N.int=50
  if(MN.SZE==0) int=2
  if(MN.SZE=="size.at.birth") int=3  
  colfunc <- colorRampPalette(c("navy", "cadetblue","white"))
  couleurs=rev(colfunc(N.int))
  BREAKS=seq(0,1,length.out=N.int+1)
  f <- function(m) t(m)[,nrow(m):1]  #function for rotating matrix
  STM.t=f(Size.Tran.Mat)
  SS=1+nrow(STM.t)-seq(1,nrow(STM.t),5)
  xx=1:nrow(STM.t)
  yy=1:ncol(STM.t)
  image(xx,yy,STM.t,ylab="",xlab="",xaxt='n',
        yaxt='n',col =couleurs,breaks=BREAKS)
  axis(2,yy,F,tck=0.025)
  axis(2,seq(1,nrow(STM.t),5),rev(TL.bin.mid[seq(1,nrow(STM.t),5)]),tck=0.05,padj=0.65,cex.axis=1.25)
  SS=1+nrow(STM.t)-seq(1,nrow(STM.t),5)
  axis(3,xx,F,tck=0.025)
  axis(3,seq(1,nrow(STM.t),5),F,tck=0.05,padj=0.65) 
  SQ=rev(seq(BREAKS[1],BREAKS[length(BREAKS)],.25))              
  if(Sex=="Females")
  {
    axis(3,seq(1,nrow(STM.t),5),TL.bin.mid[seq(1,nrow(STM.t),5)],tck=0.05,padj=0.65,cex.axis=1.25)    
    color.legend(xx[length(xx)*.9],yy[length(yy)*.35],xx[length(xx)*.99],yy[length(yy)*.99],
                 SQ,rect.col=rev(couleurs),gradient="y",col=1,cex=1)
  }
  legend("bottomleft",Sex,bty='n',cex=1.5)
  box()
  
}

Test.STM.fn=function(D,LEG,Yr_bin,TL.bin.mid)
{
  par(mfcol=c(2,1),mai=c(.4,.2,.175,.1),oma=c(1,4,2,1),las=1,mgp=c(1,0.7,0))
  SEQ=seq(1,length(TL.bin.mid),Yr_bin)
  CLS=rainbow(length(SEQ))
  fn.get.50=function(x) TL.bin.mid[findInterval(0.5,cumsum(x))]
  Ln.nm=rep(NA,length(SEQ))
  for(l in 1:length(D))
  {
    plot(TL.bin.mid,D[[1]][,1],ylim=c(0,.75),type='l',ylab="",xlab="",col=CLS[1])
    Ln.nm[1]=fn.get.50(MOD$Prob_at_len_rec)
    for(u in 2:length(SEQ)) 
    {
      id=SEQ[u]
      lines(TL.bin.mid,D[[l]][,id],col=CLS[u])
      Ln.nm[u]=fn.get.50(D[[1]][,id-1])
    }
    if(l==1) legend("topleft",as.character(Ln.nm),lty=1,col=CLS,bty='n')
    legend("topright",LEG[l],bty='n',cex=1.35)
  }
  mtext("Probability of growing from a specified initial length",3,line=0,cex=1.5,outer=T)
  mtext("Total length (cm)",1,line=-0.2,cex=1.65,outer=T)
  mtext("Probability",2,line=2,cex=1.65,las=3,outer=T)
}

fn.cpue.simple.res=function(what,what1)
{
  OBS=MOD[[match(what,names(MOD))]]
  PRED=MOD[[match(what1,names(MOD))]]  
  if(!is.matrix(OBS))
  {
    Dropd=which(OBS[]<0)
    OBS[Dropd]=PRED[Dropd]=NA
    Res=OBS-PRED   
    plot(PRED,Res,ylab="",xlab="",pch=19,cex.axis=1.5,cex=2,col=CL)
    abline(h=0,lty=2,lwd=3)
  }  
  if(is.matrix(OBS))
  {
    if(ncol(OBS)==1)
    {
      Dropd=which(OBS[,1]<0)
      OBS[Dropd,]=PRED[Dropd,]=NA
      Res=OBS-PRED   
      plot(PRED,Res,ylab="",xlab="",pch=19,cex.axis=1.5,cex=2,col=CL)
      abline(h=0,lty=2,lwd=3)
    } else
    {
      
      for(t in 1:ncol(OBS))
      {
        Ob=OBS[,t]
        Ob[Ob<0]=NA
        Prd=PRED[,t]
        Res=Ob-Prd    
        plot(Prd,Res,ylab="",xlab="",pch=19,cex.axis=1.5,cex=2,col=CL,
             ylim=c(min(c(0,Res),na.rm=T),max(c(0,Res),na.rm=T)))
        abline(h=0,lty=2,lwd=3)
        legend("top",as.character(Areas.zones$zone[t]),bty='n',cex=1.75)
      }
    }
  }  
}

fn.cpue.simple.res.time=function(what,what1)
{
  OBS=MOD[[match(what,names(MOD))]]
  PRED=MOD[[match(what1,names(MOD))]]  
  if(!is.matrix(OBS))
  {
    Dropd=which(OBS[]<0)
    OBS[Dropd]=PRED[Dropd]=NA
    Res=OBS-PRED   
    plot(1:length(Res),Res,ylab="",xlab="",pch=19,cex.axis=1.5,cex=2,col=CL,xaxt='n')
    abline(h=0,lty=2,lwd=3)
  }  
  if(is.matrix(OBS))
  {
    if(ncol(OBS)==1)
    {
      Dropd=which(OBS[,1]<0)
      OBS[Dropd,]=PRED[Dropd,]=NA
      Res=OBS-PRED   
      plot(1:length(Res),Res,ylab="",xlab="",pch=19,cex.axis=1.5,cex=2,col=CL,xaxt='n')
      abline(h=0,lty=2,lwd=3)
    } else
    {
      for(t in 1:ncol(OBS))
      {
        Ob=OBS[,t]
        Ob[Ob<0]=NA
        Prd=PRED[,t]
        Res=Ob-Prd    
        plot(1:length(Res),Res,ylab="",xlab="",pch=19,cex.axis=1.5,cex=2,col=CL,xaxt='n',
             ylim=c(min(c(0,Res),na.rm=T),max(c(0,Res),na.rm=T)))
        abline(h=0,lty=2,lwd=3)
        legend("top",as.character(Areas.zones$zone[t]),bty='n',cex=1.75)
        axis(1,1:length(Res),F,tck=-0.02)
        axis(1,seq(1,length(Res),5),F,tck=-0.04)
      }
    }
  }  
  axis(1,seq(1,length(Res),5),seq(Yrs[1],Yrs[length(Yrs)],5),tck=-0.04,cex.axis=1.5)
}

fn.normalised.cpue.res=function(what,what1,SDs)
{
  OBS=MOD[[match(what,names(MOD))]]
  PRED=MOD[[match(what1,names(MOD))]]
  Sigma=MOD[[match(SDs,names(MOD))]]
  if(is.matrix(OBS))
  {
    if(ncol(OBS)>1) par(mfrow=c(3,2),mai=c(.75,.75,.05,.05),las=1,mgp=c(2,.6,0),cex.axis=1.25)else
      par(mfrow=c(2,1),mai=c(.75,.75,.05,.05),las=1,mgp=c(2,.6,0),cex.axis=1.25)
    
    for(t in 1:ncol(OBS))
    {
      Ob=OBS[,t]
      Prd=PRED[,t]
      Sig=Sigma[,t]
      
      Id=which(Ob<0)
      if(length(Id)>0)
      {
        Ob=Ob[-Id]
        Prd=Prd[-Id]
        Sig=Sig[-Id]
      }
      
      Norm.res=log(Ob/Prd)/Sig
      hist(Norm.res,mai="",xlab="",ylab="",col=CL)
      legend("topright",paste("SDSR=",round(sd(Norm.res,na.rm=T),2)),bty='n',cex=1.5)
      box()
      if(ncol(OBS)==1)
      {
        mtext("Frequency",2,las=3,cex=1.35,line=2)
        mtext("Standardised residuals",1,cex=1.35,line=2)
      }
      if(t==2) mtext("Frequency",2,las=3,cex=1.35,line=2)
      if(t==3) mtext("Standardised residuals",1,cex=1.35,line=2)
      
      plot(Prd,Norm.res,pch=19,col=CL,cex=1.5,xlab="",ylab="")
      abline(h=0,lty=2,lwd=2,col="grey70")
      if(ncol(OBS)>1)legend("top",as.character(Areas.zones$zone[t]),bty='n',cex=1.5)
      if(ncol(OBS)==1)
      {
        mtext("Standardised residuals",2,cex=1.35,line=2,las=3)
        mtext("Predicted values",1,cex=1.35,line=2)
      }
      
      if(t==2) mtext("Standardised residuals",2,cex=1.35,line=2,las=3)
      if(t==3) mtext("Predicted values",1,cex=1.35,line=2)
    }
  }else
  {
    Id=which(OBS<0)
    if(length(Id)>0)
    {
      OBS=OBS[-Id]
      PRED=PRED[-Id]
      Sigma=Sigma[-Id]
    }
    
    Norm.res=log(OBS/PRED)/Sigma
    par(mfcol=c(2,1),mai=c(.75,.75,.05,.05),las=1,mgp=c(2,.6,0),cex.axis=1.25,cex.lab=1.5)
    hist(Norm.res,mai="",xlab="Standardised residuals",col=CL)
    legend("topright",paste("SDSR=",round(sd(Norm.res),2)),bty='n',cex=1.5)
    box()
    plot(PRED,Norm.res,pch=19,col=CL,cex=1.5,xlab="Predicted values",ylab="Standardised residuals")
    abline(h=0,lty=2,lwd=2,col="grey70")
  }
}

fn.plt.number.size.cls=function(X,bins.mid,a,CLSS,CEx)
{
  YLIM=c(0,max(a))
  if(Spatial.Str=="Single zone")
  {
    par(mfcol=c(1,1),mai=c(.85,1.,.15,.05),las=1,mgp=c(3,.6,0))
    plot(bins.mid,a[1,],type='l',ylab="Numbers (1000s of individuals)",ylim=YLIM,
         xlab="Total length class mid point (cm)",cex.axis=1.25,cex.lab=2)
    for(n in 2:nrow(a)) lines(bins.mid,a[n,],col=CLSS[n])
    legend("topright",c("Initial",Yrs),lty=1,col=c("black",CLSS),bty='n',cex=CEx)
  }else
  {
    NZn=nrow(Areas.zones)
    smart.par(n.plots=NZn,MAR=c(4,4,1,1),OMA=rep(1,4),MGP=c(2.5,.7,0))
    Nrw=nrow(a)
    NMBRS=c("Initial",Yrs)
    smrt.indx=1:length(NMBRS)
    smrt.indx=ceiling(seq(1,length(c("Initial",Yrs)),length.out=NZn+1))
    CLReS=c("black",CLSS)
    LISTA=vector('list',NZn)
    LISTA[[1]]=seq(smrt.indx[1],(smrt.indx[2]))
    for(l in 2:NZn) LISTA[[l]]=seq((smrt.indx[l]+1),(smrt.indx[l+1]))
    for(z in 1:NZn)
    {
      DD=a[seq(z,Nrw,by=NZn),]
      DD=DD[1:(length(X)+1),]
      plot(bins.mid,DD[1,],type='l',ylab="",ylim=YLIM,xlab="",cex.axis=1.25,cex.lab=2)
      for(n in 2:nrow(DD)) lines(bins.mid,DD[n,],col=CLSS[n])
      legend("top",paste(as.character(Areas.zones$zone[z])),bty='n',cex=1.5)
      legend("topright",NMBRS[LISTA[[z]]],lty=1,col=CLReS[LISTA[[z]]],bty='n',cex=1)
      
    }
    mtext("Numbers (1000s of individuals)",2,line=-1,las=3,outer=T,cex=1.5)
    mtext("Total length class mid point (cm)",1,line=-1,outer=T,cex=1.5)
    
  }
  
}

fn.gt=function(x)as.numeric(substr(x, nchar(x), nchar(x)))


fn.show.MTM=function(MATRIZ,N.int,TITLE,add.legend)
{
  colfunc <- colorRampPalette(c("grey10", "grey60","white"))
  couleurs=rev(colfunc(N.int))
  BREAKS=seq(0,1,length.out=N.int+1)
  xx=1:nrow(MATRIZ)
  yy=1:ncol(MATRIZ)
  image(xx,yy,t(MATRIZ),ylab="",xlab="",xaxt='n',yaxt='n',col =couleurs,breaks=BREAKS)
  axis(1,xx,F,tck=0.025)
  axis(2,yy,F,tck=0.025)
  mtext(TITLE,3,.25,cex=1.5)
  box()
  SQ=rev(seq(BREAKS[1],BREAKS[length(BREAKS)],.25))
  if(add.legend=="YES")color.legend(3,1,3.35,2,SQ,rect.col=rev(couleurs),gradient="y",col="black",cex=1)
  axis(1,1:nrow(MATRIZ),colnames(MATRIZ)[1:nrow(MATRIZ)],las=2,cex.axis=1.25,hadj=0.75,tck=0.025)
  axis(2,1:nrow(MATRIZ),colnames(MATRIZ)[1:nrow(MATRIZ)],cex.axis=1.25,tck=0.025,hadj=.75,las=1) 
}

fn.model.outputs=function(SCENARIO)
{
  Pth=paste("2_Outputs/Model_outputs/",SCENARIO,sep='')    
  #if(!file.exists(file.path(setPath("2_Outputs/Model_outputs"), SCENARIO))) dir.create(file.path(setPath("2_Outputs/Model_outputs"), SCENARIO))  
  setPath(Pth)
  
  
  par_zn<<-function() par(mfcol=c(3,1),mai=c(.3,.65,.1,.1),oma=c(2,2,.1,.1),las=1,mgp=c(2,.65,0),cex.axis=1.25)
  
  #Observed vs Predicted cpue
  if(!Scenarios$CPUE[i]%in%c("N/A","No"))
  {
    fn.fig("CPUE",2400,2400) 
    
    id=match("CPUE_SD_out",names(MOD))
    if(!is.na(id))
    {
      if(Spatial.Str=="Single zone")
      {
        fn.lay.single(c(.8,1,.25,.05))
        comp.fn.comp.ob.pred(what="CPUE_out",what1="Est_CPUE_out",LW=LW)    
        legend('topright',c("observed (? CV)","predicted"),bty='n',pch=c(21,NA),cex=1.5,lty=c(NA,1),
               pt.bg=c("grey80",NA),col=c("grey40","black"),lwd=2)
      }
      if(Spatial.Str=="Three zones")
      {   
        fn.comp.ob.pred.spatial("CPUE_out","Est_CPUE_out","topright",Plot.dims="YES")   
      }
    }
    if(!is.na(match("CPUE_eff",names(MOD))))
    {
      fn.lay.single(c(.8,.8,.5,.05))
      comp.fn.comp.ob.pred(what="CPUE_eff",what1="Est_CPUE_eff",LW=LW)
      legend('topright',c("observed","predicted"),bty='n',pch=c(21,NA),cex=1.25,lty=c(NA,1),
             pt.bg=c("grey80",NA),col=c("grey40","black"),lwd=2)
    }
    if(is.na(match("CPUE_eff",names(MOD))) & is.na(id))
    {
      fn.lay.single(c(.8,.8,.5,.05))
      comp.fn.comp.ob.pred(what="CPUE_out",what1="Est_CPUE_out",LW=LW)    
      if(!is.na(id))legend('topright',c("observed (? CV)","predicted"),bty='n',pch=c(21,NA),cex=1.25,lty=c(NA,1),
             pt.bg=c("grey80",NA),col=c("grey40","black"),lwd=2)else
               legend('topright',c("observed","predicted"),bty='n',pch=c(21,NA),cex=1.25,lty=c(NA,1),
                      pt.bg=c("grey80",NA),col=c("grey40","black"),lwd=2)
    }
    if(Spatial.Str=="Three zones")
    { 
      mtext("Financial year",1,outer=T,line=-2,cex=2)
      mtext("Relative cpue",2,outer=T,line=-2.25,las=3,cex=2)    
    }else
    {
      mtext("Financial year",1,outer=T,line=-1.25,cex=2)
      mtext("Relative cpue",2,outer=T,line=-1.85,las=3,cex=2)  
    }
    
    dev.off()
    
    Yrs=Rest.Yrs()  
    
    #Simple residuals
    fn.fig("CPUE.simple.residuals",2400,2400)
    if(Spatial.Str=="Single zone")
    {
      par(mfcol=c(1,1),mai=c(.4,.65,.1,.1),oma=c(2,2,.1,.1),las=1)
      if(!is.na(match("CPUE_eff",names(MOD))))  fn.cpue.simple.res(what="CPUE_eff",what1="Est_CPUE_eff")else
        fn.cpue.simple.res(what="CPUE_out",what1="Est_CPUE_out")     
      mtext("Expected value",1,outer=T,line=0.75,cex=2) 
      mtext("Residual",2,outer=T,line=0,las=3,cex=2)  
      mtext("Cpue residuals",3,line=.5,cex=2)    
    }
    if(Spatial.Str=="Three zones")
    {
      par_zn()
      fn.cpue.simple.res(what="CPUE_out",what1="Est_CPUE_out")
      mtext("Expected value",1,outer=T,line=0.75,cex=2) 
      mtext("Residual",2,outer=T,line=-0.5,las=3,cex=2)  
    }
    dev.off()
    
    #Simple residuals thru time  
    fn.fig("CPUE.simple.residuals.time",2400,2400)
    if(Spatial.Str=="Single zone")
    {
      par(mfcol=c(1,1),mai=c(.4,.65,.1,.1),oma=c(2,2,.1,.1),las=1)
      if(!is.na(match("CPUE_eff",names(MOD))))  fn.cpue.simple.res.time(what="CPUE_eff",what1="Est_CPUE_eff")else
        fn.cpue.simple.res.time(what="CPUE_out",what1="Est_CPUE_out") 
      mtext("Financial year",1,outer=T,line=0.75,cex=2)
      mtext("Residual",2,outer=T,line=0,las=3,cex=2)  
      mtext("Cpue residuals",3,line=.5,cex=2)    
    }
    if(Spatial.Str=="Three zones")
    {
      par_zn()
      fn.cpue.simple.res.time(what="CPUE_out",what1="Est_CPUE_out")
      mtext("Year",1,outer=T,line=0.75,cex=2) 
      mtext("Residual",2,outer=T,line=-.5,las=3,cex=2)  
    }
    dev.off()
    
    #Normalised cpue residuals (Francis 2011 Appendix B)
    if(!is.na(id))
    {
      fn.fig("CPUE.normalised.residuals",2400,2400)  
      fn.normalised.cpue.res(what="CPUE_out",what1="Est_CPUE_out",SDs="CPUE_SD_out")
      dev.off()
    }
  }
  #comp.fn.comp.ob.pred_logcpue(what="CPUE_out",what1="log_Est_CPUE_out",LW=LW)
  
  #Observed vs Predicted catch
  id=match("Est_TC_out",names(MOD))
  if(!is.na(id))
  {
    fn.fig("Catch",2400,2400)    
    fn.lay.single(c(.8,1.1,.05,.05))
    par(mgp=c(2.5, .9, 0))
    fn.comp.ob.pred("TC_out","Est_TC_out","topright")
    mtext("Financial year",1,outer=T,line=-1.5,cex=2)
    mtext("Catch (tonnes)",2,outer=T,line=-2.0,las=3,cex=2)
    dev.off() 
  }
  
  
  #Total biomass   
    #absolute  
  id=match("Total_biom",names(MOD))
  if(!is.na(id))
  {
    fn.fig("Biomass_total",2400,2400)
    fn.lay.single(c(.8,1.3,.05,.05))
    fn.see.Biom(X=Yrs,what="Total_biom",Spatial=Spatial.Str)
    if(Spatial.Str=="Three zones")
    { 
      mtext("Financial year",1,outer=T,line=.5,cex=2)
      mtext("Total biomass (tonnes)",2,outer=T,line=-.5,las=3,cex=2)
    }else
    {
      mtext("Financial year",1,outer=T,line=-1.5,cex=2)
      mtext("Total biomass (tonnes)",2,outer=T,line=-2,las=3,cex=2)
    }
    dev.off()
  }
    #relative
  if(!is.na(id))
  {
    fn.fig("Biomass_total_rel",2400,2400)
    fn.lay.single(c(.8,1.3,.05,.05))
    fn.see.Biom(X=Yrs,what="Total_biom_rel",Spatial=Spatial.Str)
    if(Spatial.Str=="Three zones")
    { 
      mtext("Financial year",1,outer=T,line=.5,cex=2)
      mtext("Relative biomass",2,outer=T,line=-.5,las=3,cex=2)
    }else
    {
      mtext("Financial year",1,outer=T,line=-1.5,cex=2)
      mtext("Relative biomass",2,outer=T,line=-2,las=3,cex=2)
    }
    dev.off()
  }
  
  #Female mature biomass
    #absolute
  id=match("Fem_spawn_biom",names(MOD))
  if(!is.na(id))
  {
    fn.fig("Biomass_female_mature",2400,2400) 
    fn.lay.single(c(.8,1.3,.05,.05))
    fn.see.Biom(X=Yrs,what="Fem_spawn_biom",Spatial=Spatial.Str)
    if(Spatial.Str=="Three zones")
    { 
      mtext("Financial year",1,outer=T,line=.5,cex=2)
      mtext("Female mature biomass (tonnes)",2,outer=T,line=-.5,las=3,cex=2)
    }else
    {
      mtext("Financial year",1,outer=T,line=-1.5,cex=2)
      mtext("Female mature biomass (tonnes)",2,outer=T,line=-2,las=3,cex=2)
    }
    dev.off()
  }
  #relative
  if(!is.na(id))
  {
    fn.fig("Biomass_female_mature_rel",2400,2400) 
    fn.lay.single(c(.8,1.3,.05,.05))
    fn.see.Biom(X=Yrs,what="Fem_spawn_biom_rel",Spatial=Spatial.Str)
    if(Spatial.Str=="Three zones")
    { 
      mtext("Financial year",1,outer=T,line=.5,cex=2)
      mtext("Relative biomass",2,outer=T,line=-.5,las=3,cex=2)
    }else
    {
      mtext("Financial year",1,outer=T,line=-1.5,cex=2)
      mtext("Relative biomass",2,outer=T,line=-2,las=3,cex=2)
    }
    dev.off()
  }
  
  
  #Annual fishing mortality
  id=match("Annual_F_sex",names(MOD))
  if(!is.na(id))
  {
    fn.fig("Annual_fishing_mortality",2400,2400)
    fn.lay.single(c(.8,1.1,.05,.05))
    if(SCENARIO=="Base case")fn.see.F(X=Yrs,what="Annual_F_female",what1="Annual_F_male",Spatial=Spatial.Str)else
        fn.see.pred_F(X=Yrs,what="Annual_F_sex",lista.columns=ncol(MOD$Annual_F_sex))
    if(Spatial.Str=="Three zones")
    { 
      mtext("Financial year",1,outer=T,line=.65,cex=2)
      mtext(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),
            2,outer=T,line=-1,las=3,cex=2)
    }else
    {
      mtext("Financial year",1,outer=T,line=-1.5,cex=2)
      mtext(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),
            2,outer=T,line=-2.15,las=3,cex=2)
    }
    dev.off()
  }
  id=match("Annual_F",names(MOD))
  if(!is.na(id))
  {
    fn.fig("Annual_fishing_mortality",2400,2400) 
    fn.lay.single(c(.8,1.05,.5,.05))
    fn.see.pred(Yrs,"Annual_F","topright")
    mtext("Financial year",1,outer=T,line=-1.5,cex=2)
    mtext(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),2,outer=T,line=-2.0,las=3,cex=2)
    dev.off()
  }
  id=match("Annual_F_F",names(MOD))
  if(!is.na(id))
  {
    fn.fig("Annual_fishing_mortality",2400,2400) 
    par(mfcol=c(2,1),mai=c(.8,1,.05,.1),las=1,mgp=c(1,.6,0))
    YY=max(c(max(MOD[[match("Annual_F_F",names(MOD))]]),max(MOD[[match("Annual_F_M",names(MOD))]])))
    fn.see.pred.spatial(Yrs,"Annual_F_F","topright",c(0,YY))
    legend("topleft","Females",bty='n',cex=1.5)
    
    fn.see.pred.spatial(Yrs,"Annual_F_M","topright",c(0,YY))
    legend("topleft","Males",bty='n',cex=1.5)
    
    mtext("Financial year",1,outer=T,line=-1.5,cex=2)
    mtext(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),2,
          outer=T,line=-2.25,las=3,cex=2)
    dev.off()
  }
  
  
  #Annual number of recruits
  id=match("Annual_rec",names(MOD))
  if(!is.na(id))
  {
    fn.fig("Recruitment_annual",2400,2400)
    fn.lay.single(c(.7,1,.05,.05))
    if(Scenarios[i,]$Model_type=="Length-based")fn.see.Biom(X=Yrs,what="Annual_rec",Spatial=Spatial.Str)else
    {
      if(Spatial.Str=="Single zone") fn.see.pred(Yrs,"Annual_rec",lista.columns=1)    
      if(Spatial.Str=="Three zones") fn.see.pred.spatial(Yrs,"Annual_rec","bottomleft",NA) 
    }
    mtext("Financial year",1,outer=T,line=-1.25,cex=2)
    mtext("Number of recruits (1000s of individuals)",2,outer=T,line=-1.8,las=3,cex=2)
    dev.off()
    
    #Stock-recruitment
    REC=MOD[[match("Annual_rec",names(MOD))]]
    FMB=MOD[[match("Fem_spawn_biom",names(MOD))]]
    fn.fig("Recruitment_Stock",2400,2400)
    if(Spatial.Str=="Single zone")
    {
      fn.lay.single(c(.7,1,.05,.05))
      if(class(REC)=='numeric')
      {
        REC=REC[1:length(Yrs)]
        FMB=FMB[1:length(Yrs)]
      }
      if(class(REC)=='matrix')
      {
        REC=REC[1:length(Yrs),]
        FMB=FMB[1:length(Yrs),]
      }
      SRL=data.frame(REC=REC,FMB=FMB,Year=Yrs)
      SRL=SRL[order(SRL$Year),]
      with(SRL,plot(FMB,REC,pch=19,cex=1.5,ylab="",xlab="",col=CL,
                    xlim=c(0,max(FMB,na.rm=T)),ylim=c(0,max(REC,na.rm=T))))
      #with(SRL,text(FMB,REC,Year,adj=1.25,col=1,cex=0.75,srt=-15))
      with(SRL,points(FMB[nrow(SRL)],REC[nrow(SRL)],col=2,cex=1.5,pch=19))
    }
    if(Spatial.Str=="Three zones")
    {
      par(mfcol=c(3,1),mai=c(.3,.4,.1,.2),oma=c(2,2,.1,.1),las=1,mgp=c(2,.65,0),cex.axis=1.35)
      for(o in 1:3)
      {
        SRL=data.frame(REC=REC[1:length(Yrs),o],FMB=FMB[1:length(Yrs),o],Year=Yrs)
        SRL=SRL[order(SRL$Year),]
        with(SRL,plot(FMB,REC,pch=19,cex=1.5,ylab="",xlab="",col=CL,
                      xlim=c(0,max(FMB,na.rm=T)),ylim=c(0,max(REC,na.rm=T))))
        with(SRL,points(FMB[nrow(SRL)],REC[nrow(SRL)],col=2,cex=1.5,pch=19))
        legend("bottomright",as.character(Areas.zones$zone[o]),bty="n",cex=2)
      }
    }
    legend("topleft",paste(SRL$Year[nrow(SRL)]),pch=19,col=2,text.col=2,cex=1.5,bty='n')
    
    mtext("Female mature biomass (tonnes)",1,outer=T,line=-1.25,cex=1.75)
    mtext("Number of recruits (1000s of individuals)",2,outer=T,line=-2,las=3,cex=1.75)
    dev.off()
  }
  
  
  #Size stuff
  if(!is.na(match("siz_comp_F",substr(names(MOD),1,10))))
  {
    CLSS=rainbow(length(Yrs))
    mesh_6.5.ob=c("siz_comp_F(z)","siz_comp_M(z)")
    mesh_7.ob=c("siz_comp_F_7(z)","siz_comp_M_7(z)")
    mesh_6.5.ob.N=c("N_siz_comp" ,"N_siz_comp_M")
    mesh_7.ob.N=c("N_siz_comp_7" ,"N_siz_comp_M_7")
    mesh_6.5.pred=c("Est_Prop_at_len_6_5_out(1,z)","Est_Prop_at_len_6_5_out(2,z)")
    mesh_7.pred=c("Est_Prop_at_len_7_out(1,z)","Est_Prop_at_len_7_out(2,z)")
    
    if(Spatial.Str=="Single zone")COLumn=1 
    
    if(!Spatial.Str=="Single zone")
    {
      mesh_6.5.ob=c(sapply(mesh_6.5.ob, paste, Areas.zones$area, sep="_"))
      mesh_7.ob=c(sapply(mesh_7.ob, paste, Areas.zones$area, sep="_"))
      mesh_6.5.pred=c(sapply(mesh_6.5.pred, paste, Areas.zones$area, sep="_"))
      mesh_7.pred=c(sapply(mesh_7.pred, paste, Areas.zones$area, sep="_"))
    }
    
    nm.lst.obs=list(mesh_6.5=mesh_6.5.ob,mesh_7=mesh_7.ob)
    nm.lst.obs.sam.size=list(mesh_6.5=mesh_6.5.ob.N,mesh_7=mesh_7.ob.N)
    nm.lst.pred=list(mesh_6.5=mesh_6.5.pred,mesh_7=mesh_7.pred)
    
    
    #Observed vs Predicted size composition
    for(m in 1:length(nm.lst.obs))
    {
      for(sx in 1:length(nm.lst.obs[[m]]))
      {
        if(!Spatial.Str=="Single zone")
        {
          COLumn=fn.gt(nm.lst.obs[[m]][sx])
          ZnE=as.character(Areas.zones$zone[COLumn])
        }
        SEcx=substr(nm.lst.obs[[m]][sx],10,10)
        if(SEcx=="F") S=1 else S=2
        SEcX=ifelse(SEcx=="F","female","male")
        if(Spatial.Str=="Single zone") PLT.nm=paste(names(nm.lst.obs)[m],SEcX,sep="_")else
          PLT.nm=paste(names(nm.lst.obs)[m],SEcX,ZnE,sep="_")
        
        fn.fig(paste("Size_composition",PLT.nm,sep="_"),2400,2400) 
        fn.comp.ob.pred.size.comp(what=nm.lst.obs[[m]][sx],what1=nm.lst.pred[[m]][sx],
                                  BINS="Len_bin_mdpt",where.leg="topright",N=nm.lst.obs.sam.size[[m]][S],COL=COLumn)
        mtext("Total length class mid point (cm)",1,outer=T,line=-.5,cex=1.8)
        mtext("Proportion",2,outer=T,line=-1,las=3,cex=1.8)
        dev.off()
      }
    }
    
    #Pearson residuals
    for(m in 1:length(nm.lst.obs))
    {
      for(sx in 1:length(nm.lst.obs[[m]]))
      {
        if(!Spatial.Str=="Single zone")
        {
          COLumn=fn.gt(nm.lst.obs[[m]][sx])
          ZnE=as.character(Areas.zones$zone[COLumn])
        }
        SEcx=substr(nm.lst.obs[[m]][sx],10,10)
        if(SEcx=="F") S=1 else S=2
        SEcX=ifelse(SEcx=="F","female","male")
        if(Spatial.Str=="Single zone") PLT.nm=paste(names(nm.lst.obs)[m],SEcX,sep="_")else
          PLT.nm=paste(names(nm.lst.obs)[m],SEcX,ZnE,sep="_")
        
        fn.fig(paste("Size_Pearson",PLT.nm,sep="_"),2400,2400) 
        Ns=MOD[[match(nm.lst.obs.sam.size[[m]][S],names(MOD))]][,COLumn]  
        Pearson.res(OBS=MOD[[match(nm.lst.obs[[m]][sx],names(MOD))]],
                    PRED=MOD[[match(nm.lst.pred[[m]][sx],names(MOD))]],
                    Ns,BINS=MOD[[match("Len_bin_mdpt",names(MOD))]],
                    YRS=Yrs,YLAB="Size class (cm)",
                    scaler=.5,CL.pos=CL,CL.neg="white",SeqMin=1, SeqMax=7)
        dev.off()
        
      }
    }
    
    #Numbers in size class per year
    
    #Females
    fn.fig("Size_Numbers_size_class_female",2400,2400) 
    fn.plt.number.size.cls(X=Yrs,bins.mid=MOD[[match("Len_bin_mdpt",names(MOD))]],
                           a=MOD[[match("Num_in_len_class(1)",names(MOD))]],CLSS,CEx=.8)
    dev.off()
    
    #Males
    fn.fig("Size_Numbers_size_class_male",2400,2400) 
    fn.plt.number.size.cls(X=Yrs,bins.mid=MOD[[match("Len_bin_mdpt",names(MOD))]],
                           a=MOD[[match("Num_in_len_class(2)",names(MOD))]],CLSS,CEx=.8)
    dev.off()
    
    
    #Length distribution of recruits
    fn.fig("Recruits_size_distribution",2400,2400)
    bins.mid=MOD[[match("Len_bin_mdpt",names(MOD))]]
    fn.lay.single(c(.8,1.3,.05,.05))
    fn.see.pred(bins.mid,"Prob_at_len_rec","topright")
    mtext("Total length class mid point (cm)",1,outer=T,line=-1.5,cex=2)
    mtext("Probability",2,outer=T,line=-2.0,las=3,cex=2)
    dev.off()
    
    
    #Obs vs Pred age-length 
    fn.fig("Age_length",2400,2400)
    par(mfcol=c(2,1),mai=c(.25,.2,.1,.1),oma=c(3,4,1,1),las=1,mgp=c(1,0.7,0))
    fn.obs.age.len(MOD$AgE,"L","Lpred","black")
    legend("right","Females",bty='n',cex=1.5)
    
    fn.obs.age.len(MOD$AgE_M,"L_M","Lpred_M","black")
    mtext("Age",1,outer=T,line=1,cex=2)
    mtext("TL (cm)",2,outer=T,line=2,las=3,cex=2)  
    legend("right","Males",bty='n',cex=1.5)
    dev.off()
    
    #Plot size transition matrix
    STM_F=MOD$"STM(1)"
    STM_M=MOD$"STM(2)"
    fn.fig("Size_transition_matrix",2400,2400)
    par(mfcol=c(2,1),mai=c(.2,.2,.175,.1),oma=c(1,4,3,1),las=1,mgp=c(1,0.7,0))
    fn.plt.STM(STM_F,"Females",TL.bin.mid=MOD$Len_bin_mdpt)
    fn.plt.STM(STM_M,"Males",TL.bin.mid=MOD$Len_bin_mdpt)
    mtext("To TL bin class (cm)",2,2.4,cex=1.5,outer=T,las=3)
    mtext("From TL bin class (cm)",3,1.25,cex=1.5,outer=T)
    dev.off()
    
    #Test size transition matrix
    fn.fig("Size_transition_matrix_test",2400,2400)
    Test.STM.fn(D=list(STM_F,STM_M),LEG=c("Females","Males"),Yr_bin=5,TL.bin.mid=MOD$Len_bin_mdpt)
    dev.off()
  }
  
  #Movement  
  if(Spatial.Str=="Three zones")     
  {
    MOve_mat=MOD[[match("Mov_mat",names(MOD))]]
    #Move_after_year=MOve_mat %^% 365
    Move_after_year=MOD[[match("Mov_mat_annual",names(MOD))]]
    rownames(Move_after_year)=colnames(Move_after_year)=c("West","Zone 1","Zone 2")
    fn.fig("Annual_movement",2400,2400)
    fn.lay.single(c(1,1,.05,.05))
    par(mgp=c(2,1,0))
    fn.show.MTM(MATRIZ=Move_after_year,N.int=50,"",'YES')  
    mtext("To",1,-2,cex=1.85,outer=T)
    mtext("From ",2,-1.6,cex=1.85,outer=T,las=3) 
    dev.off()
    write.csv(Move_after_year,"Annual_movement.csv")
  }
}

fn.pdf.model.outputs=function(SCENARIO)
{
  Pth=paste("2_Outputs/Model_outputs/",SCENARIO,sep='')    
  #if(!file.exists(file.path(setPath("2_Outputs/Model_outputs"), SCENARIO))) dir.create(file.path(setPath("2_Outputs/Model_outputs"), SCENARIO))  
  setPath(Pth)
  
  pdf(paste("REPORT.",SCENARIO,".pdf",sep="")) 
  
  #Parameters, gradient and negloglike
  Coefss=STD[match(estim.pars,STD$name),match(c("name","value","std.dev"),names(STD))]
  Coefss=subset(Coefss,!is.na(value))
  Loged_coef=which(substr(Coefss$name,1,2)%in%c("ln","lo"))   
  CoefsS=vector('list',nrow(Coefss))
  for(nn in 1:nrow(Coefss))
  {
    if(nn%in%Loged_coef) LGeD="YES" else LGeD="NO" 
    dummy=fn.get.CI(Coefss[nn,]$value,Coefss[nn,]$std.dev,LOGGED=LGeD)
    NME=as.character(Coefss[nn,]$name)
    if(substr(NME,1,3)%in%c("ln_")) NME=substr(NME,4,50)  
    if(substr(NME,1,2)%in%c("ln")) NME=substr(NME,3,50)
    if(substr(NME,1,4)%in%c("log_")) NME=substr(NME,5,50) 
    if(substr(NME,1,3)%in%c("ln_","log")) NME=substr(NME,4,50)
    if(NME=="tau2") NME="tau"
    CoefsS[[nn]]=data.frame(name=NME,Lower.95=dummy$low,MLE=dummy$MLE,Upper.95=dummy$Up)
  }
  CoefsS=do.call(rbind,CoefsS)
  nId=match(c("Lower.95","MLE","Upper.95"),names(CoefsS))
  CoefsS[,nId]=apply(CoefsS[,nId],2,function(x) as.character(format(round(x,8), scientific =T)))
  Max.Gradient=PARAMS$Max.gradient
  NegLogLike=PARAMS$NLL 
  Dummy=data.frame(name="",Lower.95="",MLE="",Upper.95="")
  Dummy1=data.frame(name="",Lower.95="value",MLE="",Upper.95="")
  Grad.Like=data.frame(name=c("Max.Gradient","NegLogLike"),Lower.95=c(Max.Gradient,NegLogLike),MLE="",Upper.95="")
  Par.table=as.data.frame.matrix(rbind(CoefsS,Dummy,Dummy1,Grad.Like))
  row.names(Par.table)=NULL
  grid.table(Par.table)
  
  par_zn<<-function() par(mfcol=c(3,1),mai=c(.4,.65,.2,.1),oma=c(2,2,.1,.1),las=1)
  
  #Observed vs Predicted cpue
  if(!Scenarios$CPUE[i]%in%c("N/A","No"))
  {
    id=match("CPUE_SD_out",names(MOD))
    if(!is.na(id))
    {
      if(Spatial.Str=="Single zone")
      {
        fn.lay.single(c(.8,1,.5,.05))
        fn.comp.ob.pred("CPUE_out","Est_CPUE_out","topright")  
        mtext("Observed vs Predicted cpue",3,line=.5,cex=2)
      }
      if(Spatial.Str=="Three zones")
      {   
        fn.comp.ob.pred.spatial("CPUE_out","Est_CPUE_out","topright",Plot.dims="YES")
        mtext("Observed vs Predicted cpue",3,line=.5,cex=2,outer=T)
      }
    }
    if(!is.na(match("CPUE_eff",names(MOD))))
    {
      fn.lay.single(c(.8,.8,.5,.05))
      fn.comp.ob.pred("CPUE_eff","Est_CPUE_eff","topright")  
      mtext("Observed vs Predicted cpue",3,line=.5,cex=2)
    }
    if(is.na(match("CPUE_eff",names(MOD))) & is.na(id))
    {
      fn.lay.single(c(.8,.8,.5,.05))
      fn.comp.ob.pred("CPUE_out","Est_CPUE_out","topright")  
      mtext("Observed vs Predicted cpue",3,line=.5,cex=2)
    }
    mtext("Financial year",1,outer=T,line=-1.25,cex=1.75)
    mtext("Relative cpue",2,outer=T,line=-1.85,las=3,cex=1.75)    
    
    #Reset Yrs
    Yrs=Rest.Yrs() 
    
    
    #Simple residuals
    if(Spatial.Str=="Single zone")
    {
      par(mfcol=c(1,1),mai=c(.4,.65,.1,.1),oma=c(2,2,.1,.1),las=1)
      if(!is.na(match("CPUE_eff",names(MOD))))  fn.cpue.simple.res(what="CPUE_eff",what1="Est_CPUE_eff")else
        fn.cpue.simple.res(what="CPUE_out",what1="Est_CPUE_out")  
      mtext("Expected value",1,outer=T,line=0.75,cex=2) 
      mtext("Residual",2,outer=T,line=0,las=3,cex=2)  
      mtext("Cpue residuals",3,line=.5,cex=2)    
    }
    if(Spatial.Str=="Three zones")
    {
      par_zn()
      fn.cpue.simple.res(what="CPUE_out",what1="Est_CPUE_out")
      mtext("Expected value",1,outer=T,line=0.75,cex=2) 
      mtext("Residual",2,outer=T,line=-.5,las=3,cex=2)  
      mtext("Cpue residuals",3,line=.5,cex=2,outer=T)
    }
    
    
    #Simple residuals thru time  
    if(Spatial.Str=="Single zone")
    {
      par(mfcol=c(1,1),mai=c(.4,.65,.1,.1),oma=c(2,2,.1,.1),las=1)
      if(!is.na(match("CPUE_eff",names(MOD))))  fn.cpue.simple.res.time(what="CPUE_eff",what1="Est_CPUE_eff")else
        fn.cpue.simple.res.time(what="CPUE_out",what1="Est_CPUE_out") 
      mtext("Financial Year",1,outer=T,line=0.75,cex=2) 
      mtext("Residual",2,outer=T,line=0,las=3,cex=2)  
      mtext("Cpue residuals",3,line=.5,cex=2)    
    }
    if(Spatial.Str=="Three zones")
    {
      par_zn()
      fn.cpue.simple.res.time(what="CPUE_out",what1="Est_CPUE_out")
      mtext("Financial Year",1,outer=T,line=0.75,cex=2) 
      mtext("Residual",2,outer=T,line=-.5,las=3,cex=2)  
      mtext("Cpue residuals",3,line=.5,cex=2,outer=T)
    }
    
    
    #Normalised cpue residuals (Francis 2011 Appendix B)
    if(!is.na(id)) fn.normalised.cpue.res(what="CPUE_out",what1="Est_CPUE_out",SDs="CPUE_SD_out")
  }
  
  #Observed vs Predicted catch
  id=match("Est_TC_out",names(MOD))
  if(!is.na(id))
  {
    fn.lay.single(c(.8,1.1,.5,.05))
    fn.comp.ob.pred("TC_out","Est_TC_out","topright")
    mtext("Financial year",1,outer=T,line=-1.5,cex=2)
    mtext("Catch (tonnes)",2,outer=T,line=-2.0,las=3,cex=2)
    mtext("Observed vs Predicted catch",3,line=.5,cex=2)
  }
  
  #Total biomass
  id=match("Total_biom",names(MOD))
  if(!is.na(id))
  {
    fn.lay.single(c(.8,1.3,2,.05))
    fn.see.Biom(X=Yrs,what="Total_biom",Spatial=Spatial.Str)
    mtext("Financial year",1,outer=T,line=-0.5,cex=2)
    mtext("Tonnes",2,outer=T,line=-1.0,las=3,cex=2)
    mtext("Total biomass",3,outer=T,line=-1.75,cex=1.5)
  }
  
  #Female mature biomass
  id=match("Fem_spawn_biom",names(MOD))
  if(!is.na(id))
  {
    fn.lay.single(c(.8,1.3,1.25,.05))
    fn.see.Biom(X=Yrs,what="Fem_spawn_biom",Spatial=Spatial.Str)
    mtext("Financial year",1,outer=T,line=-0.5,cex=2)
    mtext("Tonnes",2,outer=T,line=-1.0,las=3,cex=2)
    mtext("Female mature biomass",3,outer=T,line=-1.75,cex=1.5)
  }
  

  #Annual fishing mortality    
  id=match("Annual_F_sex",names(MOD))
  if(!is.na(id))
  {
    fn.lay.single(c(.8,1.15,1.25,.05))
    fn.see.pred_F(X=Yrs,what="Annual_F_sex",lista.columns=ncol(MOD$Annual_F_sex))
    mtext("Financial year",1,outer=T,line=0,cex=2)
    mtext(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),2,outer=T,line=-1,las=3,cex=2)
    mtext("Fishing mortality",3,line=-1.75,cex=1.5,outer=T)  
  }
  id=match("Annual_F",names(MOD))
  if(!is.na(id))
  {
    fn.lay.single(c(.8,1.15,.5,.05))
    fn.see.pred(Yrs,"Annual_F","topright")
    mtext("Financial year",1,outer=T,line=-1.5,cex=2)
    mtext(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),2,outer=T,line=-2.0,las=3,cex=2)
    mtext("Fishing mortality",3,line=.5,cex=2)    
  }
  id=match("Annual_F_F",names(MOD))
  if(!is.na(id))
  {
    par(mfcol=c(2,1),mai=c(.75,.75,.1,.1),las=1)
    YY=max(c(max(MOD[[match("Annual_F_F",names(MOD))]]),max(MOD[[match("Annual_F_M",names(MOD))]])))
    fn.see.pred.spatial(Yrs,"Annual_F_F","topright",c(0,YY))
    legend("topleft","Females",bty='n',cex=1.5)
    
    fn.see.pred.spatial(Yrs,"Annual_F_M","topright",c(0,YY))
    legend("topleft","Males",bty='n',cex=1.5)
    
    mtext("Financial year",1,outer=T,line=-1.5,cex=2)
    mtext(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),2,
          outer=T,line=0,las=3,cex=2)
    
  }
  
  #Annual number of recruits    
  id=match("Annual_rec",names(MOD))
  if(!is.na(id))
  {
    par(mfcol=c(1,1),mai=c(.6,.7,.5,.05),las=1,mgp=c(1,.6,0))
    if(Scenarios[i,]$Model_type=="Length-based")fn.see.Biom(X=Yrs,what="Annual_rec",Spatial=Spatial.Str)else
    {
      if(Spatial.Str=="Single zone") fn.see.pred(Yrs,"Annual_rec",lista.columns=1)    
      if(Spatial.Str=="Three zones") fn.see.pred.spatial(Yrs,"Annual_rec","bottomleft",NA)   
    }
    mtext("Financial year",1,outer=T,line=-1,cex=2)
    mtext("Number of recruits (1000s of individuals)",2,outer=T,line=-0.5,las=3,cex=2)
    mtext("Annual number of recruits",3,line=.5,cex=2)
    
    #Stock-recruitment
    REC=MOD[[match("Annual_rec",names(MOD))]]
    FMB=MOD[[match("Fem_spawn_biom",names(MOD))]]
    if(Spatial.Str=="Single zone")
    {
      par(mfcol=c(1,1),mai=c(.6,.7,.5,.1),las=1,mgp=c(1,.6,0))
      if(class(REC)=='numeric')
      {
        REC=REC[1:length(Yrs)]
        FMB=FMB[1:length(Yrs)]
      }
      if(class(REC)=='matrix')
      {
        REC=REC[1:length(Yrs),]
        FMB=FMB[1:length(Yrs),]
      }
      SRL=data.frame(REC=REC,FMB=FMB,Year=Yrs)
      SRL=SRL[order(SRL$Year),]
      with(SRL,plot(FMB,REC,pch=19,cex=1.5,ylab="",xlab="",col=CL,cex.axis=1.25,
                    xlim=c(0,max(FMB,na.rm=T)),ylim=c(0,max(REC,na.rm=T))))
      #with(SRL,text(FMB,REC,Year,adj=1.25,col=1,cex=0.75,srt=-15))
      with(SRL,points(FMB[nrow(SRL)],REC[nrow(SRL)],col=2,cex=1.5,pch=19))
      with(SRL,text(FMB[nrow(SRL)],REC[nrow(SRL)],Year[nrow(SRL)],adj=1.25,col=2,cex=0.75,srt=-15))
    }
    if(Spatial.Str=="Three zones")
    {
      par_zn()
      for(o in 1:3)
      {
        SRL=data.frame(REC=REC[1:length(Yrs),o],FMB=FMB[1:length(Yrs),o],Year=Yrs)
        SRL=SRL[order(SRL$Year),]
        with(SRL,plot(FMB,REC,pch=19,cex=1.5,ylab="",xlab="",col=CL,cex.axis=1.25,
                      xlim=c(0,max(FMB,na.rm=T)),ylim=c(0,max(REC,na.rm=T))))
        with(SRL,points(FMB[nrow(SRL)],REC[nrow(SRL)],col=2,cex=1.5,pch=19))
        with(SRL,text(FMB[nrow(SRL)],REC[nrow(SRL)],Year[nrow(SRL)],adj=1.25,col=2,cex=0.75,srt=-15))
        legend("bottomright",as.character(Areas.zones$zone[o]),bty="n",cex=2)
      }
    }
    mtext("Female mature biomass (tonnes)",1,outer=T,line=-1,cex=1.75)
    mtext("Number of recruits (1000s of individuals)",2,outer=T,line=-0.5,las=3,cex=1.75)
    mtext("Stock-recruitment",3,line=.5,cex=2,outer=T)
  }
  
  
  if(!is.na(match("siz_comp_F",substr(names(MOD),1,10))))
  {
    CLSS=rainbow(length(Yrs))
    
    mesh_6.5.ob=c("siz_comp_F(z)","siz_comp_M(z)")
    mesh_7.ob=c("siz_comp_F_7(z)","siz_comp_M_7(z)")
    mesh_6.5.ob.N=c("N_siz_comp" ,"N_siz_comp_M")
    mesh_7.ob.N=c("N_siz_comp_7" ,"N_siz_comp_M_7")
    mesh_6.5.pred=c("Est_Prop_at_len_6_5_out(1,z)","Est_Prop_at_len_6_5_out(2,z)")
    mesh_7.pred=c("Est_Prop_at_len_7_out(1,z)","Est_Prop_at_len_7_out(2,z)")
    
    if(Spatial.Str=="Single zone")COLumn=1 
    
    if(!Spatial.Str=="Single zone")
    {
      mesh_6.5.ob=c(sapply(mesh_6.5.ob, paste, Areas.zones$area, sep="_"))
      mesh_7.ob=c(sapply(mesh_7.ob, paste, Areas.zones$area, sep="_"))
      mesh_6.5.pred=c(sapply(mesh_6.5.pred, paste, Areas.zones$area, sep="_"))
      mesh_7.pred=c(sapply(mesh_7.pred, paste, Areas.zones$area, sep="_"))
    }
    
    nm.lst.obs=list(mesh_6.5=mesh_6.5.ob,mesh_7=mesh_7.ob)
    nm.lst.obs.sam.size=list(mesh_6.5=mesh_6.5.ob.N,mesh_7=mesh_7.ob.N)
    nm.lst.pred=list(mesh_6.5=mesh_6.5.pred,mesh_7=mesh_7.pred)
    
    
    #Observed vs Predicted size composition
    for(m in 1:length(nm.lst.obs))
    {
      for(sx in 1:length(nm.lst.obs[[m]]))
      {
        if(!Spatial.Str=="Single zone")
        {
          COLumn=fn.gt(nm.lst.obs[[m]][sx])
          ZnE=as.character(Areas.zones$zone[COLumn])
        }
        SEcx=substr(nm.lst.obs[[m]][sx],10,10)
        if(SEcx=="F") S=1 else S=2
        SEcX=ifelse(SEcx=="F","female","male")
        
        fn.comp.ob.pred.size.comp(what=nm.lst.obs[[m]][sx],what1=nm.lst.pred[[m]][sx],
                                  BINS="Len_bin_mdpt",where.leg="topright",N=nm.lst.obs.sam.size[[m]][S],COL=COLumn)
        mtext("Total length class mid point (cm)",1,outer=T,line=-.5,cex=1.8)
        mtext("Proportion",2,outer=T,line=-1,las=3,cex=1.8)
        if(Spatial.Str=="Single zone")
        {
          mtext(paste(names(nm.lst.obs)[m]," Observed vs Predicted size composition (",
                      SEcX,")",sep=""),3,line=-.65,cex=1.25,outer=T)
        }else
          mtext(paste(names(nm.lst.obs)[m]," Observed vs Predicted size composition (",
                      SEcX,", ",ZnE,")",sep=""),3,line=-.65,cex=1.25,outer=T)
        
      }
    }
    
    #Pearson residuals
    for(m in 1:length(nm.lst.obs))
    {
      for(sx in 1:length(nm.lst.obs[[m]]))
      {
        if(!Spatial.Str=="Single zone")
        {
          COLumn=fn.gt(nm.lst.obs[[m]][sx])
          ZnE=as.character(Areas.zones$zone[COLumn])
        }
        SEcx=substr(nm.lst.obs[[m]][sx],10,10)
        if(SEcx=="F") S=1 else S=2
        SEcX=ifelse(SEcx=="F","female","male")
        
        Ns=MOD[[match(nm.lst.obs.sam.size[[m]][S],names(MOD))]][,COLumn]  
        Pearson.res(OBS=MOD[[match(nm.lst.obs[[m]][sx],names(MOD))]],
                    PRED=MOD[[match(nm.lst.pred[[m]][sx],names(MOD))]],
                    BINS=MOD[[match("Len_bin_mdpt",names(MOD))]],Ns,
                    YRS=Yrs,YLAB="Size class (cm)",
                    scaler=.5,CL.pos=CL,CL.neg="white",SeqMin=1, SeqMax=7)
        
        if(Spatial.Str=="Single zone")
        {
          mtext(paste(names(nm.lst.obs)[m]," Pearson residuals (",
                      SEcX,")",sep=""),3,line=-.65,cex=1.5,outer=T)
        }else
          mtext(paste(names(nm.lst.obs)[m]," Pearson residuals (",
                      SEcX,", ",ZnE,")",sep=""),3,line=-.65,cex=1.5,outer=T)
      }
    }
    
    #Numbers in size class per year
    
      #Females
    fn.plt.number.size.cls(X=Yrs,bins.mid=MOD[[match("Len_bin_mdpt",names(MOD))]],
                           a=MOD[[match("Num_in_len_class(1)",names(MOD))]],CLSS,CEx=.6)
    mtext(paste("Numbers in size class per year (females)",sep=""),3,line=-.5,cex=1.35,outer=T)
    
      #Males
    fn.plt.number.size.cls(X=Yrs,bins.mid=MOD[[match("Len_bin_mdpt",names(MOD))]],
                           a=MOD[[match("Num_in_len_class(2)",names(MOD))]],CLSS,CEx=.6)
    mtext(paste("Numbers in size class per year (males)",sep=""),3,line=-.5,cex=1.35,outer=T)
    
    
    #Length distribution of recruits
    bins.mid=MOD[[match("Len_bin_mdpt",names(MOD))]]
    fn.lay.single(c(.8,1.3,.5,.05))
    fn.see.pred(bins.mid,"Prob_at_len_rec","topright")
    mtext("Total length class mid point (cm)",1,outer=T,line=-1.5,cex=2)
    mtext("Probability",2,outer=T,line=-2.0,las=3,cex=2)
    mtext("Length distribution of recruits",3,line=.5,cex=2)
    
    
    #Obs vs Pred age-length    
    par(mfcol=c(2,1),mai=c(.25,.2,.3,.1),oma=c(3,4,2,1),las=1,mgp=c(1,0.7,0))
    fn.obs.age.len(MOD$AgE,"L","Lpred","black")
    #legend("bottomleft",c("observed","predicted"),bty='n',pch=c(19,NA),cex=1.5,lty=c(NA,1),col=c("black",CL),lwd=3)
    legend("bottomright",paste("SD=",MOD$sd_growth),bty='n',cex=1.25)
    legend("right","Females",bty='n',cex=1.5)
    
    fn.obs.age.len(MOD$AgE_M,"L_M","Lpred_M","black")
    legend("bottomright",paste("SD=",MOD$sd_growth),bty='n',cex=1.25)
    mtext("Age",1,outer=T,line=1,cex=2)
    mtext("TL (cm)",2,outer=T,line=2,las=3,cex=2)  
    legend("right","Males",bty='n',cex=1.5)
    mtext("Obs vs Pred age-length",3,line=.5,cex=2,outer=T)
    
    #Plot size transition matrix
    STM_F=MOD$"STM(1)"
    STM_M=MOD$"STM(2)"
    par(mfcol=c(2,1),mai=c(.2,.2,.175,.1),oma=c(1,4,5,1),las=1,mgp=c(1,0.7,0))
    fn.plt.STM(STM_F,"Females",TL.bin.mid=MOD$Len_bin_mdpt)
    fn.plt.STM(STM_M,"Males",TL.bin.mid=MOD$Len_bin_mdpt)
    mtext("To TL bin class (cm)",2,2.4,cex=1.5,outer=T,las=3)
    mtext("From TL bin class (cm)",3,1.25,cex=1.5,outer=T)
    mtext("Size transition matrix",3,line=3,cex=2,outer=T)
    
    #Test size transition matrix
    Test.STM.fn(D=list(STM_F,STM_M),LEG=c("Females","Males"),Yr_bin=5,TL.bin.mid=MOD$Len_bin_mdpt)
    
  }
  
  dev.off()
}

fn.par.table=function(Scens,Rshp)
{
  Par.tbl=NULL
  #Parameters MLE and SE
  for(ii in 1:nrow(Scens))
  {
    estim.pars=names(Pin.pars[[match(Scens$Model[ii],names(Pin.pars))]])
    STD=Store.std[[match(Scens$Model[ii],names(Store.std))]]
    Coefss=STD[match(estim.pars,STD$name),match(c("name","value","std.dev"),names(STD))]
    Coefss=subset(Coefss,!is.na(value))
    Coefss$Scenario=Scens$Model[ii]
    names(Coefss)=c("Parameter","MLE","SD","Scenario")
    Coefss=Coefss[,match(c("Scenario","Parameter","MLE","SD"),names(Coefss))]
    Coefss[,match(c("MLE","SD"),names(Coefss))]=Coefss[,match(c("MLE","SD"),names(Coefss))]
    Par.tbl=rbind(Par.tbl,Coefss)
  }  
  Loged_coef=which(substr(Par.tbl$Parameter,1,2)%in%c("ln","lo")) 
  CoefsS=vector('list',nrow(Coefss))
  for(nn in 1:nrow(Par.tbl))
  {
    if(nn%in%Loged_coef) LGeD="YES" else LGeD="NO" 
    dummy=fn.get.CI(Par.tbl[nn,]$MLE,Par.tbl[nn,]$SD,LOGGED=LGeD)
    NME=as.character(Par.tbl[nn,]$Parameter)
    if(substr(NME,1,3)%in%c("ln_")) NME=substr(NME,4,50)  
    if(substr(NME,1,2)%in%c("ln")) NME=substr(NME,3,50)
    if(substr(NME,1,4)%in%c("log_")) NME=substr(NME,5,50) 
    if(substr(NME,1,3)%in%c("ln_","log")) NME=substr(NME,4,50)
    if(NME=="tau2") NME="tau"
    CoefsS[[nn]]=data.frame(Scenario=Par.tbl[nn,]$Scenario,Parameter=NME,Lower.95=dummy$low,MLE=dummy$MLE,Upper.95=dummy$Up)
  }
  CoefsS=do.call(rbind,CoefsS)
  if(Rshp=="YES")
  {
    Par.tbl$MLE_SD=with(Par.tbl,paste(MLE,SD,sep=" ?"))
    wide <- reshape(subset(Par.tbl,select=c(Scenario,Parameter,MLE_SD)), 
                    v.names = "MLE_SD", idvar = "Scenario",timevar = "Parameter", direction = "wide")
    wide[is.na(wide)]="N/A"
  }
  fn.word.table(WD=getwd(),TBL=CoefsS,Doc.nm="MLE parameter estimates",caption=NA,paragph=NA,
                HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
                Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
}

Rest.Yrs=function()
{
  if(Show.yrs=="DATA")Yrs=yr.start:yr.end  
  if(Show.yrs=="FUTURE") Yrs=yr.start:(yr.end+Yrs.future)  
  return(Yrs)
}

density.fun=function(what)
{
  D=density(what,adjust = 2,from=0, to=1)
  plot(D$x,D$y,main='',ylab="",xlab="",type='l',lwd=2,cex.axis=1.5,yaxt='n')
  
  #Prob above ref point
  id=which(D$x<B.target)
  X=D$x[id]
  Y=D$y[id]
  polygon(c( X, B.target),  c(Y,0 ), col=CL)
  
  #Prob above ref point
  f=ecdf(what)
  Prob.below=f(B.target)
  Prob.above=round(1-Prob.below,2)
  Prob.below=round(Prob.below,2)
  
  abline(v=B.target,lty=2)
  text(B.target*.95,mean(D$y),paste("P(<",B.target,")=",Prob.below,sep=""),pos=2,cex=1.75,col="brown4")
  text(B.target*1.05,mean(D$y),paste("P(>",B.target,")=",Prob.above,sep=""),pos=4,cex=1.75,col="brown4")
  
  
}
fn.exp.less=function(x,symb,y) bquote("P"["<"][.(x)]~.(paste0(symb,y)))
fn.exp.more=function(x,symb,y) bquote("P"[">"][.(x)]~.(paste0(symb,y)))
density.fun2=function(what,B.ref,CEX)
{
  #Prob above ref point
  f=ecdf(what)
  Prob.below=f(B.ref)
  SEQ=seq(0,1,0.001)
  f.range=f(SEQ)
  id=which.min(abs(SEQ - B.ref))
  X=SEQ[1:id]
  Y=f.range[1:id]
  plot(SEQ,f.range,main='',ylab="",xlab="",type='l',lwd=2,cex.axis=1)
  polygon(c(X,rev(X)),c(Y,rep(0,length(Y))),col=CL)
  Prob.above=round(1-Prob.below,3)
  Prob.below=round(Prob.below,3)
 # if(Prob.below<=0.4 & Prob.below>0.01)SRT=65
 # if(Prob.below<=0.01) SRT=70 else SRT=0
  if(Prob.below<0.4) SRT=70 else SRT=0
  
  abline(v=B.ref,lty=2,lwd=2,col="grey50")
  if(Prob.below>=0.01|Prob.below==0)text(B.ref*.99,mean(f.range),fn.exp.less(B.ref,"=",Prob.below),pos=2,cex=CEX,col="brown4",srt=SRT)else
    text(B.ref*.99,mean(f.range),fn.exp.less(B.ref,"<","0.01"),pos=2,cex=CEX,col="brown4",srt=SRT)
  
  if(Prob.above>=0.01|Prob.above==0)text(B.ref*1.05,mean(f.range),fn.exp.more(B.ref,"=",Prob.above),pos=4,cex=CEX,col="brown4",srt=0)else
    text(B.ref*1.05,mean(f.range),fn.exp.more(B.ref,"<","0.01"),pos=4,cex=CEX,col="brown4",srt=0)
}
Probs.ref.point=function(what)
{
  f=ecdf(what)
  
  P.below.target=f(B.target)
  P.below.threshold=f(B.threshold)
  P.below.limit=f(B.limit)
  
  P.above.target=1-P.below.target
  P.above.threshold=1-P.below.threshold
  P.above.limit=1-P.below.limit
  
  P.between.thre.tar=P.below.target-P.below.threshold
  P.between.lim.thre=P.below.threshold-P.below.limit
  
  return(data.frame(P.above.target=P.above.target,
                    P.between.thre.tar=P.between.thre.tar,
                    P.between.lim.thre=P.between.lim.thre,
                    P.below.limit=P.below.limit))
}

#3. Create output pdf for each model
#LTYp=c(1,3,5)
LTYp=c(1,1,1)
Lcol=2:4
LW=2
for( i in 1:nrow(Scenarios))
{
    MOD=Store.Reports[[i]]
    PARAMS=Model.Out.sumry[[i]]
    estim.pars=names(Pin.pars[[match(Scenarios$Model[i],names(Pin.pars))]])
    STD=Store.std[[i]]
    Spatial.Str=as.character(Scenarios$Spatial_structure[i])
    
    #convert quantities in kg to tonnes for Age structured model
    if(Scenarios$Model_type[i]=="Age-structured")
    {
      kg.to.tons=c("Est_TC_out","TC_out","Total_biom","Fem_spawn_biom","Vul_biom",
                   "Virgin_Total_biom","Virgin_Vul_biom","Virgin_Spawn_Biom")
      for(kk in 1:length(kg.to.tons))
      {
        ID.kk=match(kg.to.tons[kk],names(MOD))
        if(!is.na(ID.kk))MOD[[ID.kk]]=MOD[[ID.kk]]/1000
      }
      STD$name=as.character(STD$name)
      STD$value=with(STD,ifelse(name%in%kg.to.tons,value/1000,value))
      STD$std.dev=with(STD,ifelse(name%in%kg.to.tons,std.dev/1000,std.dev))
    }

    #set yrs
    if(Show.yrs=="DATA")
    {
      Yrs=yr.start:yr.end  
      if(Scenarios$CPUE[i]=="Stand.1988")Yrs=Yrs[match("1988-89",Cpue.all$Finyear):nrow(Cpue.all)]     
    }
    
    if(Show.yrs=="FUTURE") 
    {
      Yrs=yr.start:(yr.end+Yrs.future)   
      if(Scenarios$CPUE[i]=="Stand.1988")Yrs=Yrs[match("1988-89",Cpue.all$Finyear):nrow(Cpue.all)]
    }
    
    fn.pdf.model.outputs(SCENARIO=as.character(Scenarios$Model[i])) 
  }

#reset years   
Yrs=Rest.Yrs()


#4 Compare Parameter estimates and biomass among scenarios   
setPath("2_Outputs/Compare_scenarios")
if(Run.all.Scenarios=="YES")
{
  #Table of parameter estimates
  fn.par.table(Scens=Scenarios,Rshp="NO")
  
  #Biomass comparison
  if(Do.cols=="YES")colfunc <- colorRampPalette(c("black", "darkolivegreen3"))
  if(Do.cols=="NO")colfunc <- colorRampPalette(c("black", "grey80"))

  fn.comp.bios=function(VAR,scaler,VAR1,WHERE,add.lgn,CeX)
  {
    ThisZen=match(Zen[[1]],names(Store.Reports))
    DaTT=Store.Reports[ThisZen]
    DaTT.std=Store.std[ThisZen]
    VAR1=VAR1[ThisZen]
    miniL= sapply(DaTT, '[[', VAR, simplify = FALSE, USE.NAMES = FALSE)
    Init=rep(NA,length(VAR1))
    for( i in 1:length(Init)) 
    {
      d=subset(DaTT.std[[i]],name==VAR1[i])
      if(nrow(d)==1) d=d$value else
        d=sum(d$value)
      Init[i]=d
      nn=ncol(miniL[[i]])
      if(is.matrix(miniL[[i]])) if(nn>1) miniL[[i]]=rowSums(miniL[[i]])
    }
    for( i in 1:length(miniL)) miniL[[i]]=miniL[[i]]/scaler
    
    COLs.SCEN=colfunc(length(miniL))
    LTy.SCEN=1:length(miniL)
    
    fn.plt=function(Y) plot(Yrs,Y,ylim=c(0,MaxY),type='l',lwd=2.5,xlab='',xaxt='n',cex.axis=1.25,ylab='',col=COLs.SCEN[1])
    MaxY=1
    
    fn.plt(miniL[[1]][1:length(Yrs)]/Init[[1]]) #relative terms
    for( i in 2:length(miniL))
    {
      aa=miniL[[i]][1:length(Yrs)]/Init[[i]]  #relative terms
      aa[aa>1]=1
      lines(Yrs,aa,lwd=2.5,col=COLs.SCEN[i],lty=LTy.SCEN[i])
    }
    if(add.lgn=="YES")
    {
      MTCH=subset(Scenarios,Model%in%names(miniL))
      This.var=names(Zen)
      if(This.var=="Model_type") This.var=c(This.var,"Spatial_structure")
      LGn=MTCH[,match(This.var,names(MTCH))]
      if(is.data.frame(LGn))
      {
        LGn$Spatial_structure=with(LGn,ifelse(Spatial_structure=="Single zone","1 zone",
                            ifelse(Spatial_structure=="Three zones","3 zones",NA)))
        LGn$Model_type=with(LGn,ifelse(Model_type=="Length-based","Len. based",
                            ifelse(Model_type=="Biomass dynamics","Bio. dyn.",
                            ifelse(Model_type=="Age-structured","Age str.",NA))))
        LGn=paste("(",LGn[,1],", ",LGn[,2],")",sep="")
        
      }else
      {LGn=paste("(",LGn,")",sep="")}
      LGn=paste(names(miniL),LGn,sep=" ")
      No.scen=which(Init==0)
      if(length(No.scen)>0)
      {
        LGn=LGn[-No.scen]
        legend(WHERE,LGn,bty='n',cex=CeX,lwd=3,col=COLs.SCEN[-No.scen],lty=LTy.SCEN[-No.scen])
      }else
      {
        legend(WHERE,LGn,bty='n',cex=CeX,lwd=3,col=COLs.SCEN,lty=LTy.SCEN)
      }  
    }
      
    axis(1,Yrs,F,tck=-0.01)
    axis(1,seq(Yrs[1],Yrs[length(Yrs)],5),F,tck=-0.02,cex.axis=1.25)
  }
  fn.mTxt=function()mtext("Financial year",1,1,outer=T,cex=2) 
  fn.mTxtY=function()mtext("Relative biomass",2,0.45,outer=T,cex=2,las=3)
  fn.Axs=function()axis(1,seq(Yrs[1],Yrs[length(Yrs)],5),seq(Yrs[1],Yrs[length(Yrs)],5),tck=-0.02,cex.axis=1.25)
  
  
    #Total biomass  
      #figures combined
  fn.fig("Compare_biomass_total",2400,2200)
  smart.par(n.plots=length(Compr.grup),MAR=c(2.2,2,1,1),OMA=c(2.5,2.5,.5,.5),MGP=c(1,.6,0))
  Bo=Scenarios$Model_type
  names(Bo)=names(Store.Reports)
  Bo=ifelse(Bo=="Biomass dynamics","Virgin_Total_biom",
            ifelse(Bo=="Age-structured","Virgin_Total_biom",
            ifelse(Bo=="Length-based","Virgin_Total_biom",
            NA)))
  for(n in 1:length(Compr.grup))
  {
    Zen=Compr.grup[n]
    fn.comp.bios(VAR="Total_biom",scaler=1e0,VAR1=Bo,WHERE="bottomleft",add.lgn="YES",CeX=1.1)
    fn.mTxt();fn.mTxtY();fn.Axs()
    mtext(Title.Compr.grup[n],3,cex=1)
  }
  dev.off() 
  
  # by figure  
  NEW=paste(getwd(),"Biomasses", sep="/")
  if(!dir.exists(NEW))dir.create(NEW)
  for(n in 1:length(Compr.grup))
  {
    fn.fig(paste("Biomasses/Compare_biomass_total",Title.Compr.grup[n],sep="_"),2400,2200)
    par(las=1,mgp=c(2,.7,0))
    Zen=Compr.grup[n]
    fn.comp.bios(VAR="Total_biom",scaler=1e0,VAR1=Bo,WHERE="bottomleft",add.lgn="YES",CeX=1.25)
    mtext("Financial year",1,3,cex=2)
    mtext("Relative biomass",2,2.5,cex=2,las=3)
    fn.Axs()
    mtext(Title.Compr.grup[n],3,1,cex=2)
    dev.off()
  }

  
    #Mature biomass
      #figures combined
  fn.fig("Compare_biomass_mature_female",2400,2200)
  smart.par(n.plots=length(Compr.grup),MAR=c(2.2,2,1,1),OMA=c(2.5,2.5,.5,.5),MGP=c(1,.6,0))
  Bo=Scenarios$Model_type
  names(Bo)=names(Store.Reports)
  Bo=ifelse(Bo=="Biomass dynamics",NA,
            ifelse(Bo=="Age-structured","Virgin_Spawn_Biom",
                   ifelse(Bo=="Length-based","Virgin_Spawn_Biom",
                          NA)))
  for(n in 1:length(Compr.grup))
  {
    Zen=Compr.grup[n]
    fn.comp.bios(VAR="Fem_spawn_biom",scaler=1e0,VAR1=Bo,WHERE="bottomleft",add.lgn="YES",CeX=1.1)
    fn.mTxt();fn.mTxtY();fn.Axs()
    mtext(Title.Compr.grup[n],3,cex=1)
  }
  dev.off() 
  
      # by figure     
  for(n in 1:length(Compr.grup))
  {
    fn.fig(paste("Biomasses/Compare_biomass_mature_female",Title.Compr.grup[n],sep="_"),2400,2200)
    par(las=1,mgp=c(2,.7,0))
    Zen=Compr.grup[n]
    fn.comp.bios(VAR="Fem_spawn_biom",scaler=1e0,VAR1=Bo,WHERE="bottomleft",add.lgn="YES",CeX=1.25)
    mtext("Financial year",1,3,cex=2)
    mtext("Relative biomass",2,2.5,cex=2,las=3)
    fn.Axs()
    mtext(Title.Compr.grup[n],3,1,cex=2)
    dev.off()
  }
  

   #4.2 Compare Observed vs Predicted cpue among scenarios 
  THeSe=which(!Scenarios$CPUE%in%c("No","N/A"))
  fn.fig("CPUEs",2400,2400) 
  par(mfrow=n2mfrow(length(THeSe)),mai=c(.3,.4,.1,.1),oma=c(2,2,.5,.5),las=1,mgp=c(2,.85,0))
  for( i in THeSe)
  {
      MOD=Store.Reports[[i]]
      Fze=Par.phases[[i]]
      Spatial.Str=as.character(Scenarios$Spatial_structure[i])
      id=match("CPUE_SD_out",names(MOD))
      
      if(Spatial.Str=="Three zones")
      {
        fn.comp.ob.pred.spatial("CPUE_out","Est_CPUE_out","topright",Plot.dims="NO")
      }else
      {
        if(!is.na(match("CPUE_eff",names(MOD))))
        {
          what="CPUE_eff"
          what1="Est_CPUE_eff"
        }else
        {
          what="CPUE_out"
          what1="Est_CPUE_out"
        }
        comp.fn.comp.ob.pred(what=what,what1=what1,LW=LW)
      }
      legend("topright",paste(Scenarios$Model[i]),bty='n',cex=1.5,inset=c(-.01,-0.05),text.col="black")
   }
  mtext("Financial year",1,outer=T,line=0.75,cex=1.75)
  mtext("Relative cpue",2,outer=T,line=0.1,las=3,cex=1.75)  
  dev.off()

}

#reset yrs
Yrs=Rest.Yrs()  


#5. Extract base case MLE parameter estimates and quantities
  #Base case only
if(Run.all.Scenarios=="NO")
{
  i=match(ID.base.Model,Scenarios$Model)
  MOD=Store.Reports[[i]]
  Fze=Par.phases[[i]]
  PARAMS=Model.Out.sumry[[i]]
  estim.pars=names(Pin.pars[[match(Scenarios$Model[i],names(Pin.pars))]])
  STD=Store.std[[i]]
  Spatial.Str=as.character(Scenarios$Spatial_structure[i])
  
  #set yrs
  if(Show.yrs=="DATA")
  {
    Yrs=yr.start:yr.end  
    if(Scenarios$CPUE[i]=="Stand.1988")Yrs=Yrs[match("1988-89",Cpue.all$Finyear):nrow(Cpue.all)]     
  }
  if(Show.yrs=="FUTURE") 
  {
    Yrs=yr.start:(yr.end+Yrs.future)   
    if(Scenarios$CPUE[i]=="Stand.1988")Yrs=Yrs[match("1988-89",Cpue.all$Finyear):nrow(Cpue.all)]
  }
  
  setPath(paste("2_Outputs/Model_outputs/",Scenarios[i,]$Model,sep=""))
  fn.par.table(Scens=Scenarios[i,],Rshp="NO") 
  fn.model.outputs(Scenarios[i,]$Model)
}

  #All scenarios separately    
if(Run.all.Scenarios=="YES")
{
  for( i in 1:nrow(Scenarios))
  {
    MOD=Store.Reports[[i]]
    Fze=Par.phases[[i]]
    PARAMS=Model.Out.sumry[[i]]
    estim.pars=names(Pin.pars[[match(Scenarios$Model[i],names(Pin.pars))]])
    STD=Store.std[[i]]
    Spatial.Str=as.character(Scenarios$Spatial_structure[i])
    
    #convert quantities in kg to tonnes for Age structured model
    if(Scenarios$Model_type[i]=="Age-structured")
    {
      kg.to.tons=c("Est_TC_out","TC_out","Total_biom","Fem_spawn_biom","Vul_biom",
                   "Virgin_Total_biom","Virgin_Vul_biom","Virgin_Spawn_Biom")
      for(kk in 1:length(kg.to.tons))
      {
        ID.kk=match(kg.to.tons[kk],names(MOD))
        if(!is.na(ID.kk))MOD[[ID.kk]]=MOD[[ID.kk]]/1000
      }
      STD$name=as.character(STD$name)
      STD$value=with(STD,ifelse(name%in%kg.to.tons,value/1000,value))
      STD$std.dev=with(STD,ifelse(name%in%kg.to.tons,std.dev/1000,std.dev))
    }
    
    #set yrs
    if(Show.yrs=="DATA")
    {
      Yrs=yr.start:yr.end  
      if(Scenarios$CPUE[i]=="Stand.1988")Yrs=Yrs[match("1988-89",Cpue.all$Finyear):nrow(Cpue.all)]     
    }
    if(Show.yrs=="FUTURE") 
    {
      Yrs=yr.start:(yr.end+Yrs.future)   
      if(Scenarios$CPUE[i]=="Stand.1988")Yrs=Yrs[match("1988-89",Cpue.all$Finyear):nrow(Cpue.all)]
    }
    setPath(paste("2_Outputs/Model_outputs/",Scenarios[i,]$Model,sep=""))
    fn.par.table(Scens=Scenarios[i,],Rshp="NO") 

    fn.model.outputs(SCENARIO=Scenarios[i,]$Model)
  }
}

#set yrs
Yrs=Rest.Yrs() 

#6. Simulate data from size transition matrix
if(Sim.trans.Mat=="YES")
{
  fn.simul.test.STM=function(MAT,n.sim,MinAge,MaxAge,bin.low,bin.up,GROWTH.pars,Lo,plt.what,Dist.Birth)
  {
    #1. get culumative distribution for each size class
    MAT.cum=apply(MAT, 2, cumsum)
    colnames(MAT.cum)=bin.low
    
    #2. generate random ages
    Age=round(runif(n.sim,MinAge+1,MaxAge))
    
    #3. generate random size at birth
    if(Dist.Birth=='unif') Birth.s=round(runif(n.sim,bin.low[1],bin.up[3]))    
    if(Dist.Birth=='norm')
    {
      Birth.s=round(rnorm(n.sim,Lo*1.25,Lo*0.35))
      #Birth.s=ifelse(Birth.s<Lo,Lo,Birth.s)
    }
    
    
    #4. find initial size bin
    Init.s=cbind(Birth.s,findInterval(Birth.s, c(bin.low, bin.up[length(bin.up)])))
    Init.s=Init.s[Init.s[,1]>0,]
    n.sim=nrow(Init.s)
    #5. define function for finding size class
    find.size=function(SIZE)
    {
      r=runif(1,0.1,1)
      get.cum.trans=MAT.cum[,match(paste(SIZE),colnames(MAT.cum))]  
      names(get.cum.trans)=colnames(MAT.cum)
      new.Size=as.numeric(names(get.cum.trans[findInterval(r, get.cum.trans)]))
      return(new.Size)
    }
    
    #6. loop thru timesteps (ages as per vonB)
    Store=vector('list',n.sim)
    for(i in 1:n.sim)
    {
      n=Age[i]+1
      Size.dummy=rep(NA,n)
      Size.dummy[1]=bin.low[Init.s[i,2]]
      if(Age[i]>1)for(a in 2:n) Size.dummy[a]=find.size(Size.dummy[a-1])
      Store[[i]]=data.frame(Age=0:Age[i],Len=Size.dummy)
    }
    AGE=0:MaxAge
    von.B=Lo+(GROWTH.pars$TL_inf-Lo)*(1-exp(-(GROWTH.pars$k*AGE)))
    CLS=rainbow(n.sim)
    plot(AGE,von.B,lwd=2,xlab="",ylab="",type='l',ylim=c(0,max(von.B)*1.5)) 
    if (plt.what=="all")for(i in 1:n.sim) points(Store[[i]]$Age,Store[[i]]$Len,pch=19,col=2,cex=1.25)
    if (plt.what=="last")for(i in 1:n.sim) points(Store[[i]]$Age[Age[i]],Store[[i]]$Len[Age[i]],pch=19,col=2,cex=1.25)
  }
  
  i=match(ID.base.Model,Scenarios$Model)
  MOD=Store.Reports[[i]]
  STD=Store.std[[i]]
  setPath(Scen="1_Inputs/Visualise data")
  
  fn.fig(paste(getwd(),"/preds.from.Input.STM",sep=""),2400,2400) 
  
  par(mfcol=c(2,1),mai=c(.4,.4,.1,.1),oma=c(2,2,2,1),las=1,mgp=c(1.5,.6,0))
  
  #females  
  STM.F=MOD$"STM(1)"
  Growth.F=as.data.frame(t(subset(STD,name%in%c("k","lnLinf","sd_growth"), select=c(value))))
  colnames(Growth.F)=c("k","TL_inf","SD")
  Growth.F$TL_inf=exp(Growth.F$TL_inf)
  
  fn.simul.test.STM(MAT=STM.F,n.sim=5000,MinAge=0,MaxAge=Max.age.F,bin.low=TL.bin.low,bin.up=TL.bin.up,
                    GROWTH.pars=Growth.F,Lo=Min.size.TL$a,plt.what="all",Dist.Birth='norm')
  legend("bottomright","Females",bty='n',cex=2)
  
  
  #males
  STM.M=MOD$"STM(2)"
  Growth.M=as.data.frame(t(subset(STD,name%in%c("k_M","lnLinf_M","sd_growth"), select=c(value))))
  colnames(Growth.M)=c("k","TL_inf","SD")
  Growth.M$TL_inf=exp(Growth.M$TL_inf)
  
  fn.simul.test.STM(MAT=STM.M,n.sim=5000,MinAge=0,MaxAge=Max.age.M,bin.low=TL.bin.low,bin.up=TL.bin.up,
                    GROWTH.pars=Growth.M,Lo=Min.size.TL$a,plt.what="all",Dist.Birth='norm')
  legend("bottomright","Males",bty='n',cex=2)
  
  mtext("TL (cm)",2,outer=T,las=3,line=.5,cex=1.5)
  mtext("Age",1,outer=T,cex=1.5)
  dev.off()
  
}

#7. Run scenarios of future projections 
if(Run.future.proj=="YES")
{
  Store.std.future=vector('list',length(Future.ktch.scen))
  names(Store.std.future)=names(Future.ktch.scen)
  Store.Reports.future=Store.std.future
  
  for(future in 1:length(Future.ktch.scen))
  {
    setPath(Scenarios[match("Base case",Scenarios$Model),]$Model)
    A=getwd()
    Folder=paste("Future.projections",names(Future.ktch.scen)[future],sep="_")
    setwd(file.path(getwd(), Folder))
    
    
    #copy .pin and .tpl to folder
    file.copy(paste(A,paste("/",Spec,".tpl",sep=""),sep=""), paste(Spec,".tpl",sep=""),overwrite =T)
    file.copy(paste(A,paste("/",Spec,".pin",sep=""),sep=""), paste(Spec,".pin",sep=""),overwrite =T)
    
    #run .tpl
    args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep=""),Arguments,sep=""), sep="")
    clean_admb(Spec)
    compile_admb(Spec,verbose=T)
    system(args)
    
    #Bring in outputs
    Store.std.future[[future]]=read.table(paste(Spec,".std",sep=""), header =TRUE)
    Store.Reports.future[[future]]=reptoRlist(paste(Spec,".rep",sep=""))  
  }
  
  setPath(paste("2_Outputs/Model_outputs/","Base case",sep=""))
    #total biomass
  fn.fig("Future_catch_scenarios_total_biomass",2000,2400) 
  par(mfcol=c(length(Future.ktch.scen),1),mai=c(.4,.4,.1,.1),oma=c(2,2,2,1),las=1,mgp=c(1.5,.6,0))
  for(future in 1:length(Future.ktch.scen))
  {
    STD=Store.std.future[[future]]
    XX=c(Yrs,(Yrs[length(Yrs)]+1):(Yrs[length(Yrs)]+Yrs.future))
    fn.see.Biom(X=XX,what="Total_biom_rel",Spatial="Single zone")
    abline(h=B.target)
    legend("bottomleft",paste("future catches= ",round(Store.Reports.future[[future]]$TC_out[length(XX)]),
              " (tonnes)",sep=""),bty='n',cex=2)
  }
  mtext("Financial year",1,outer=T,line=.5,cex=2)
  mtext("Relative biomass",2,outer=T,line=-.3,las=3,cex=2)
  dev.off()
  
    #mature female biomass
  fn.fig("Future_catch_scenarios_mature_female_biomass",2000,2400) 
  par(mfcol=c(length(Future.ktch.scen),1),mai=c(.4,.4,.1,.1),oma=c(2,2,2,1),las=1,mgp=c(1.5,.6,0))
  for(future in 1:length(Future.ktch.scen))
  {
    STD=Store.std.future[[future]]
    XX=c(Yrs,(Yrs[length(Yrs)]+1):(Yrs[length(Yrs)]+Yrs.future))
    fn.see.Biom(X=XX,what="Fem_spawn_biom_rel",Spatial="Single zone")
    abline(h=B.target)
    legend("bottomleft",paste("future catches= ",round(Store.Reports.future[[future]]$TC_out[length(XX)]),
                              " (tonnes)",sep=""),bty='n',cex=2)
  }
  mtext("Financial year",1,outer=T,line=.5,cex=2)
  mtext("Relative biomass",2,outer=T,line=-.3,las=3,cex=2)
  dev.off()
  

}


# Section F: RUN BASE CASE MCMC AND DISPLAY OUTPUTS -----------------------
#notes:    First compile the ADMB code once to create the .exe file    
#         mceval stops abruptly, execute manually and then import chains
if(DO.MCMC=="YES")
{
  #clean up
  if(exists("Store.Models"))rm(Store.Models,Store.Reports,Store.std,Model.Out.sumry)
  
  
  library(coda)  
  library(lattice)
  
  #1. Run MCMC
  fn.source("ADMB_MCMC.estimation.R")
  
  #define chain
  #note: it takes ~0.047 secs per iteration
  
  #Choose scenario to run
  MDL=Scenarios[match(ID.base.Model,Scenarios$Model),]$Model
  setPath(MDL)
  Start.time=Sys.time()
  system.time(fn.run.MCMC(MODEL=Spec,DATA=paste(Spec,".dat",sep=""),nSims,Thin))  
  End.time=Sys.time()
  
  #email running time
  Tot.time=round(difftime(End.time,Start.time,units="mins"),2)
  function.send.email(
    to ="Matias.Braccini@fish.wa.gov.au",
    subject ="Model run_MCMC",
    body =paste("MCMC finished. Base Case model with",nSims,"simulations took",Tot.time,"minutes to run"),                     
    Attachment=NULL)
  
  
  #2. Analyse chain and show posteriors
  fn.source("ADMB_MCMC_vect.post.format.R")
  
  #2.1 Read in MCMC posteriors of params and quantities of interest 
  #parameter estimates 
  Output=read.delim("PARS.mcmc")    #The PARS.mcmc object is created within the .tpl during mceval phase
  Output=Output[-burning,]        #removing burning period
  mcmcObject=as.mcmc(Output)
  
  #get quantities                                  
  Post.names <- Sys.glob("*.mcmc")
  Post.names=Post.names[-match("PARS.mcmc",Post.names)]
  Posteriors=vector('list',length(Post.names))
  names(Posteriors)=unlist(lapply(strsplit(Post.names,"[.]"), function(x) x[1]))
  
  #remove burning period and put in nice format   
  for(p in 1:length(Posteriors))
  {
    dummy=read.delim(Post.names[[p]]) 
    dummy=vect.post.format(dummy,burning)
    Posteriors[[p]]=dummy
    rm(dummy)
  }
  
  
  #2.2 Stats
  setPath(paste("2_Outputs/Model_outputs/",ID.base.Model,sep=""))
  if(!file.exists(file.path(getwd(), "/MCMC"))) dir.create(file.path(getwd(), "/MCMC"))   #create folder
  setwd(file.path(getwd(), "/MCMC"))
  
  Post.stats=summary(mcmcObject)
  Table.param.posterior.stats=cbind(Post.stats$statistics,Post.stats$quantiles)
  Table.param.posterior.stats=Table.param.posterior.stats[,match(c("Mean","SD","Time-series SE",
                    "2.5%","50%","97.5%"),colnames(Table.param.posterior.stats))]
  Table.param.posterior.stats=as.data.frame(Table.param.posterior.stats)
  Table.param.posterior.stats=cbind(data.frame(Parameter=rownames(Table.param.posterior.stats)),
                                    round(Table.param.posterior.stats,4))
  rownames(Table.param.posterior.stats)=NULL
  names(Table.param.posterior.stats)=c("Parameter","Mean","SD","SE","2.5% Cred. I.","Median","97.5% Cred. I.")
  
  fn.word.table(WD=getwd(),TBL=Table.param.posterior.stats,Doc.nm="Posterior estimated pars",
                caption=NA,paragph=NA,HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',
                body.fnt.sze=10,Zebra='NO',Zebra.col='grey60',Grid.col='black',
                Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman")
  
  
  
  #2.3 Diagnostics
  Plt.dim=function(MFROW) par(mfrow=MFROW,las=1,mar=c(4,2,1,1),omi=c(.5,.5,0.1,0.1),cex.lab=1.5)
  
  #plot chain for estimable pars
  fn.fig("Parameter_chains",2000,2600)
  xyplot(mcmcObject,lwd=1.5)
  dev.off()
  
  #plot posterior for estimable pars
  fn.fig("Parameter_densities",2000,2600)
  densityplot(mcmcObject,lwd=2)
  dev.off()
  
  
  #correlatin between parameters
  fn.fig("Parameter_correlation",2400,2400)
  plot(Output)
  dev.off()
  

  #Chain convergence diagnostics
    # Raftery
  #note: estimates the number of iterations needed for a given level of precision in 
  #     posterior samples, as well as estimating burn-in,. 
  #     Values of I (the 'dependence factor') >5 indicate strong autocorrelation which may be 
  #     due to a poor choice of starting value. If the number of iterations in data is 
  #     too small, an error message is printed indicating the minimum length of pilot run. 
  #     The minimum length is the required sample size for a chain with no correlation 
  #     between consecutive samples. 
  capture.output(raftery.diag(mcmcObject), file = "Chain_diag_Raftery.diag_I.less.5.no.autocorrelation.txt")
  
    # Geweke
  #note: geweke.diag returns Z scores for the equality of (by default) the first 10%
  #     and the last 50% of the chain. Do two-tailed Z-test to test for Sig difference between the
  #     two parts of the chain ( p value 0 in this case)
  Gewek=geweke.diag(mcmcObject)
  Gewek=Gewek[[1]]
  Gewek.out=sapply(Gewek,function(x) pnorm(abs(x),lower.tail=FALSE)*2)
  capture.output(Gewek.out, file = "Chain_diag_Gewek.out_Prob.more.0.05.ok.txt")
  
    #Effective size
  #note: This value should be at least 100, and probably greater than 200, for reasonable
  #          estimation of confidence intervals
  capture.output(effectiveSize(mcmcObject), file = "Chain_diag_Effective.size_should.be.more.than.100.txt")
  
  # Highest posterior density (i.e. Bayesian credible) intervals:
  capture.output(HPDinterval(mcmcObject), file = "Chain_diag_HPDinterval.txt")
  

  #Check for chain autocorrelation  
  Post.par=as.character(Table.param.posterior.stats$Parameter)
  fn.fig("Chain_diag_autocorrelation",2400,2400)
  par(mfrow=n2mfrow(length(Post.par)),mar=c(2,2.5,.5,.5),oma=c(2,2,1,1),las=1,mgp=c(2,.8,0),cex.axis=1.25,col.lab="white")
  for(x in 1:length(Post.par))
  {
    autocorr.plot(mcmcObject[ , x], auto.layout = FALSE,lwd = 4, col = "red", lag.max = 100)
    mtext(Post.par[x], side = 3,cex = 1.25, line = -1.5)
    abline(h=0)     
  }
  mtext("Autocorrelation",2,line=0,las=3,outer=T,cex=1.5)
  dev.off()



  #- Posterior time series of quantities of interest                            
  Display.plot.index=function(Post,Post.init,AXS,Add.Re.Pnt,Relative,Do.polys)
  {
    Post.rel=Post[1:length(All.yrs),]
    if(!is.na(Post.init[1]) & Relative=="YES") for(nnn in 1:ncol(Post)) Post.rel[,nnn]=Post.rel[,nnn]/Post.init[nnn]
    if(!is.na(Post.init[1]) & Relative=='NO')  btarget=median(B.target*Post.init)
    if(Relative=='YES') btarget=B.target
    Plot.Index(DAT=Post.rel,FINYEAR=All.yrs,COL,Colr,Colr2,AXS,Add.Re.Pnt,Do.polys,btarget,do.rel="NO")  
  }
  
  Display.plot.index.rel=function(Post,AXS,Add.Re.Pnt,Do.polys)
  {
    Post.rel=Post[1:length(All.yrs),]
    btarget=B.target
    bthreshold=B.threshold
    blimit=B.limit
    Plot.Index(DAT=Post.rel,FINYEAR=All.yrs,COL,Colr,Colr2,AXS,Add.Re.Pnt,Do.polys,btarget,bthreshold,blimit,do.rel="YES")  
  }
  
  Plot.Index=function(DAT,FINYEAR,COL,Colr,Colr2,Do.axis,Add.Re.Pnt,do.polygons,btarget,bthreshold,blimit,do.rel)
  {
    MEAN=rowMeans(DAT,na.rm=T)
    MEDIAN=apply(DAT,1,median,na.rm=T)
    SD=apply(DAT,1,sd,na.rm=T)
    CV=SD/MEAN
    LOW95=apply(DAT,1, function(x) quantile(x, 0.025,na.rm=T))
    UP95=apply(DAT,1, function(x) quantile(x, 0.975,na.rm=T))
    LOW50=apply(DAT,1, function(x) quantile(x, 0.25,na.rm=T))
    UP50=apply(DAT,1, function(x) quantile(x, 0.75,na.rm=T))
    
    MAX.y=max(UP95)
    if(do.rel=="YES") MAX.y=1
    N.YRs=length(FINYEAR)
    YR=1:N.YRs
    
    if(do.polygons=="YES")
    {
      #polygons
      Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1])
      Biom.Vec <- c(LOW50, tail(UP50, 1), rev(UP50), LOW50[1])
      plot(YR,MEDIAN,type='l',lty=1,col=COL,lwd=3,ylim=c(0,MAX.y),xlim=c(YR[1]-1,YR[length(YR)]+1),
           xaxt='n',xlab="",ylab="",cex.axis=1.25,xaxs="i")
      polygon(Year.Vec, Biom.Vec, col = Colr, border = Colr2,lwd=1)
      lines(YR,MEDIAN,lty=1,col=COL,lwd=3)
      lines(YR,LOW95,lty=2,col=Colr2,lwd=2)
      lines(YR,UP95,lty=2,col=Colr2,lwd=2)
    }else
    {
      #points  
      plot(YR,MEDIAN,cex=1.5,pch=19,col=CL,lwd=3,ylim=c(0,MAX.y),xlim=c(YR[1]-1,YR[length(YR)]+1),
           xaxt='n',xlab="",ylab="",cex.axis=1.25,xaxs="i")
      segments(YR,LOW95,YR,UP95,col=CL)
      axis(1,YR,labels=F,tck=-0.015)
      
      #future projections
      Yr.future=(yr.end+1):(yr.end+Yrs.future)
      ftre.indx=match(Yr.future,FINYEAR)
      points(YR[ftre.indx],MEDIAN[ftre.indx],pch=19,col="brown4",cex=1.7)
      segments(YR[ftre.indx],LOW95[ftre.indx],YR[ftre.indx],UP95[ftre.indx],col="brown4")
    }
    if(Do.axis=="NO")axis(1,seq(1,N.YRs,5),labels=F,tck=-0.030,cex.axis=1.25)
    if(Do.axis=="YES")axis(1,seq(1,N.YRs,5),labels=FINYEAR[seq(1,N.YRs,5)],tck=-0.030,cex.axis=1.25)
    
    #add reference point
    if(Add.Re.Pnt=="YES")
    {
      if(do.polygons=="YES")
      {
        B.target.Vec <- c(rep(0,N.YRs), btarget, rep(btarget,N.YRs), 0)
        polygon(Year.Vec, B.target.Vec, col = rgb(0.01, 0, 0.4,alpha =0.15), border = "transparent")
      }else
      {
        abline(h=btarget,lwd=2,col='grey30',lty=2)
        text(YR[3],btarget,"Target",pos=3,cex=1.25)
        abline(h=bthreshold,lwd=2,col='grey30',lty=2)
        text(YR[5],bthreshold,"Threshold",pos=3,cex=1.25)  
        abline(h=blimit,lwd=2,col='grey30',lty=2)
        text(YR[3],blimit,"Limit",pos=3,cex=1.25)
      }
        
      
      #Prob above ref point
      # Now=nrow(DAT)
      # f=ecdf(DAT[Now,])
      # Prob.below=f(btarget)
      # Prob.above=round(1-Prob.below,2)
      # arrows(YR[Now]*.9,btarget*.65,YR[Now],btarget,length=0.1,col="brown4",lwd=1.5)
      # text(YR[Now]*.8,btarget*.65,paste("P(>",B.target,")=",Prob.above,sep=""),pos=1,cex=1.5,col="brown4")
    }
    
  }
  
  hist.fun=function(what,MAXY)
  {
    hist(what,col=Colr2,ylab="",xlab="",main="",cex.axis=1.75,xlim=c(0,MAXY))
    box()
  }
  

  All.yrs=yr.start:(yr.end+Yrs.future)       
  COL="dodgerblue4"
  Colr="lightsteelblue1"
  Colr2="lightsteelblue4"
  
  #biomasses
    #a. absolute
  fn.fig("Biomass",2000,2400)
  par(mfcol=c(2,1),mai=c(.1,.1,.1,.1),oma=c(3,4.5,2,.1),las=1,mgp=c(1,.6,0))
  
  # total biomass     
  Display.plot.index(Post=Posteriors[[match("Total_biom",names(Posteriors))]],
                     Post.init=Posteriors[[match("Virgin_Total_biom",names(Posteriors))]],
                     AXS="NO",Add.Re.Pnt="NO",Relative="NO",Do.polys="NO")
  mtext("Total  (tonnes)",2,cex=1.75,las=3,line=3.5)
  
  # Mature biomass
  Display.plot.index(Post=Posteriors[[match("Fem_spawn_biom",names(Posteriors))]],
              Post.init=Posteriors[[match("Virgin_Spawn_Biom",names(Posteriors))]],
              AXS="YES",Add.Re.Pnt="NO",Relative="NO",Do.polys="NO")
  mtext("Mature female (tonnes)",2,cex=1.75,las=3,line=3.5)
  
  mtext("Biomass",3,line=0,cex=1.5,outer=T)
  mtext("Financial year",1,line=2,cex=1.75)
  
  dev.off()
  
  
  #b. relative
  fn.fig("Biomass_relative",2000,2400)
  par(mfcol=c(2,1),mai=c(.1,.1,.1,.1),oma=c(3,4,.1,.1),las=1,mgp=c(1,.6,0))
  
  # total biomass     
  Display.plot.index.rel(Post=Posteriors[[match("Total_biom_rel",names(Posteriors))]],
                     AXS="NO",Add.Re.Pnt="YES",Do.polys="NO")
  mtext("Total",2,cex=1.75,las=3,line=2.5)
  
  # Mature biomass
  Display.plot.index.rel(Post=Posteriors[[match("Fem_spawn_biom_rel",names(Posteriors))]],
                     AXS="YES",Add.Re.Pnt="YES",Do.polys="NO")
  mtext("Mature female",2,cex=1.75,las=3,line=2.5)
  
  #mtext("Relative biomass",3,line=0,cex=1.5,outer=T)
  mtext("Financial year",1,line=2,cex=1.75)
  dev.off()
  
  
  #Annual recruitment
  fn.fig("Annual recruitment",2000,2000)
  par(mfcol=c(1,1),mai=c(.1,.1,.1,.1),oma=c(3,4,.1,.1),las=1,mgp=c(1,.6,0))
  Display.plot.index(Post=Posteriors[[match("Annual_rec",names(Posteriors))]],
                     Post.init=NA,AXS="YES",Add.Re.Pnt="NO",Relative="NO",Do.polys="NO")
  mtext("Number of recruits (1000s of individuals)",2,cex=1.5,las=3,line=2.5)
  mtext("Financial year",1,line=2,cex=1.5)
  dev.off()

    
  #Fishing mortality
  Display.plot.F=function(DAT,FINYEAR,Do.axis,MAX.y)
  {
    DAT=DAT[1:length(All.yrs),]
    MEAN=rowMeans(DAT,na.rm=T)
    MEDIAN=apply(DAT,1,median,na.rm=T)
    SD=apply(DAT,1,sd,na.rm=T)
    CV=SD/MEAN
    LOW95=apply(DAT,1, function(x) quantile(x, 0.025,na.rm=T))
    UP95=apply(DAT,1, function(x) quantile(x, 0.975,na.rm=T))
    LOW50=apply(DAT,1, function(x) quantile(x, 0.25,na.rm=T))
    UP50=apply(DAT,1, function(x) quantile(x, 0.75,na.rm=T))
    
    N.YRs=length(FINYEAR)
    YR=1:N.YRs
    
    plot(YR,MEDIAN,cex=1.5,pch=19,col=CL,lwd=3,ylim=c(0,MAX.y),xlim=c(YR[1]-1,YR[length(YR)]+1),
         xaxt='n',xlab="",ylab="",cex.axis=1.25,xaxs="i")
    segments(YR,LOW95,YR,UP95,col=CL)
    axis(1,YR,labels=F,tck=-0.015)
    
    if(Do.axis=="NO")axis(1,seq(1,N.YRs,5),labels=F,tck=-0.030,cex.axis=1.25)
    if(Do.axis=="YES")axis(1,seq(1,N.YRs,5),labels=FINYEAR[seq(1,N.YRs,5)],tck=-0.030,cex.axis=1.25)
  }
  
  Y.F=Posteriors[[match("Annual_F_F",names(Posteriors))]]
  Y.F=Y.F[1:length(All.yrs),]
  Y.M=Posteriors[[match("Annual_F_M",names(Posteriors))]]
  Y.M=Y.M[1:length(All.yrs),]
  

  fn.fig("Fishing mortality",2000,2000)
  par(mfcol=c(2,1),mai=c(.1,.1,.1,.1),oma=c(3,4,.1,.1),las=1,mgp=c(1,.6,0))
  
    #females
  Display.plot.F(DAT=Posteriors[[match("Annual_F_F",names(Posteriors))]],
                 FINYEAR=All.yrs,Do.axis="NO",MAX.y=quantile(c(Y.F,Y.M),.9999))
  
  legend('topright',"Females",bty='n',cex=1.5)

    #males
  Display.plot.F(DAT=Posteriors[[match("Annual_F_M",names(Posteriors))]],
                 FINYEAR=All.yrs,Do.axis="YES",MAX.y=quantile(c(Y.F,Y.M),.9999))
  legend('topright',"Males",bty='n',cex=1.5)
  
  mtext(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),2,las=3,
        line=2,cex=1.75,outer=T)   
  mtext("Financial year",1,line=1.75,cex=1.75,outer=T)
  dev.off()
  
  
  
  #- Posterior distributions current year  
  Current=Yrs[length(Yrs)]
  id.current=match(Current,All.yrs)  

    #Depletion    
  fn.fig("Posterior_current_year_depletion",1400,2400)
  par(mfcol=c(3,2),mar=c(3.5,2,.1,.1),oma=c(.5,1.5,1.5,.5),mgp=c(1,.6,0),las=1)
  
  # Total biomass
  density.fun2(what=Posteriors$Total_biom_rel[length(Yrs),],B.ref=B.target,CEX=1.2)
  mtext("Total",3,cex=1.25,line=0)
  density.fun2(what=Posteriors$Total_biom_rel[length(Yrs),],B.ref=B.threshold,CEX=1.2) 
  density.fun2(what=Posteriors$Total_biom_rel[length(Yrs),],B.ref=B.limit,CEX=1.2) 
  
  # Mature biomass
  density.fun2(what=Posteriors$Fem_spawn_biom_rel[length(Yrs),],B.ref=B.target,CEX=1.2) 
  legend('bottomright',"Target",bty='n',cex=2)
  mtext("Mature female",3,cex=1.25,line=0)
  density.fun2(what=Posteriors$Fem_spawn_biom_rel[length(Yrs),],B.ref=B.threshold,CEX=1.2) 
  legend('bottomright',"Threshold",bty='n',cex=2)
  density.fun2(what=Posteriors$Fem_spawn_biom_rel[length(Yrs),],B.ref=B.limit,CEX=1.2) 
  legend('bottomright',"Limit",bty='n',cex=2)
  
  mtext("Probability",2,cex=1.5,las=3,line=-0.2,outer=T)
  mtext(paste(Current,"relative biomass"),1,line=-1,cex=1.5,outer=T)
  
  dev.off()
  
    #export probabilities for Weight of Evidence
  write.csv(Probs.ref.point(what=Posteriors$Total_biom_rel[length(Yrs),]),"Consequence_likelihood_total.csv",row.names=F)
  write.csv(Probs.ref.point(what=Posteriors$Fem_spawn_biom_rel[length(Yrs),]),"Consequence_likelihood_mature.female.csv",row.names=F)
  
    # total biomass
  # fn.fig("Posterior_current_year_depletion",2000,2400)
  # par(mfcol=c(2,1),mai=c(.2,.1,.2,.1),oma=c(3,3,.1,.1),las=1,mgp=c(1,.75,0))
  # density.fun(what=Posteriors[[match("depletion",names(Posteriors))]])
  # legend("topright","Total",cex=2,bty='n')
  # 
  #   # Mature biomass
  #  density.fun(what=Posteriors[[match("depletion_spawn",names(Posteriors))]])
  # legend("topright","Female mature",cex=2,bty='n')
  # mtext("Density",2,cex=2.5,line=1,outer=T,las=3)
  # mtext(paste(Current,"relative biomass"),1,line=2.5,cex=2.5)
  # dev.off()
  
  

  #- Posterior distributions future year  
  Future=yr.end+Yrs.future
  id.future=match(Future,All.yrs)  
  
    #Depletion    
  fn.fig("Posterior_future_year_depletion",1400,2400)
  par(mfcol=c(3,2),mar=c(3.5,2,.1,.1),oma=c(.5,1.5,1.5,.5),mgp=c(1,.6,0),las=1)
  
      # Total biomass
  density.fun2(what=Posteriors$Total_biom_rel[id.future,],B.ref=B.target,CEX=1.2)
  mtext("Total",3,cex=1.25,line=0)
  density.fun2(what=Posteriors$Total_biom_rel[id.future,],B.ref=B.threshold,CEX=1.2) 
  density.fun2(what=Posteriors$Total_biom_rel[id.future,],B.ref=B.limit,CEX=1.2) 
  
      # Mature biomass
  density.fun2(what=Posteriors$Fem_spawn_biom_rel[id.future,],B.ref=B.target,CEX=1.2) 
  legend('bottomright',"Target",bty='n',cex=2)
  mtext("Mature female",3,cex=1.25,line=0)
  density.fun2(what=Posteriors$Fem_spawn_biom_rel[id.future,],B.ref=B.threshold,CEX=1.2) 
  legend('bottomright',"Threshold",bty='n',cex=2)
  density.fun2(what=Posteriors$Fem_spawn_biom_rel[id.future,],B.ref=B.limit,CEX=1.2) 
  legend('bottomright',"Limit",bty='n',cex=2)
  
  mtext("Probability",2,cex=1.5,las=3,line=-0.2,outer=T)
  mtext(paste(Future,"relative biomass"),1,line=-1,cex=1.5,outer=T)
  
  dev.off()
  
    #export probabilities for Weight of Evidence
  write.csv(Probs.ref.point(what=Posteriors$Total_biom_rel[id.future,]),"Consequence_likelihood_total_future.csv",row.names=F)
  write.csv(Probs.ref.point(what=Posteriors$Fem_spawn_biom_rel[id.future,]),"Consequence_likelihood_mature.female_future.csv",row.names=F)
  
  
  #Output mean annual realtive biomass
  Get.Index.stats=function(DAT,FINYEAR)
  {
    MEAN=rowMeans(DAT,na.rm=T)
    MEDIAN=apply(DAT,1,median,na.rm=T)
    LOW95=apply(DAT,1, function(x) quantile(x, 0.025,na.rm=T))
    UP95=apply(DAT,1, function(x) quantile(x, 0.975,na.rm=T))
    return(cbind(Yrs=FINYEAR,MEDIAN=MEDIAN,MEAN=MEAN,LOW95=LOW95,UP95=UP95))
    
  }
      #total 
  write.csv(Get.Index.stats(Posteriors$Total_biom_rel[1:length(All.yrs),],All.yrs),"rel.biom.base.case.csv",row.names=F)
  
    #Female mature
  write.csv(Get.Index.stats(Posteriors$Fem_spawn_biom_rel[1:length(All.yrs),],All.yrs),"rel.biom.base.case_female_mature.csv",row.names=F)
}


# Section G: RUN CATCH-MSY METHOD AND DISPLAY OUTPUTS ---------------------
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/Catch_MSY.R")
if(DO.CatchMSY=="YES")
{
  #r prior
  fn.source("Leslie.matrix.R") 
  library(MASS)
  Rprior=fun.Leslie(N.sims=10000,k=Growth.F$k,Linf=Growth.F$FL_inf,Aver.T=TEMP,
                    A=Max.age.F,first.age=0,RangeMat=Age.50.mat,Rangefec=Fecundity,
                    sexratio=0.5,Reprod_cycle=Breed.cycle,Hoenig.only="NO")  
  
  #get mean and SD from gamma distribution
  gamma.pars=suppressWarnings(fitdistr(Rprior, "gamma"))  
  shape=gamma.pars$estimate[1]        
  rate=gamma.pars$estimate[2]      
  r.pars=c(shape,rate)
  
  #get mean and SD from lognormal distribution
  #LogN.pars=fitdistr(Rprior, "lognormal")  
  #log_mean.r=LogN.pars$estimate[1]    #already in log space     
  #log_sd.r=LogN.pars$estimate[2]      #already in log space     
  #r.pars=c(log_mean.r,log_sd.r)
  
  for(sc in 1:length(ktch_msy_scen)) if(!is.na(ktch_msy_scen[[sc]]$r.prior)[1]) ktch_msy_scen[[sc]]$r.prior=r.pars
  
  
  #run catch_msy function                   #takes 0.0001 seconds per iteration per scenario
  setPath(paste("2_Outputs/Model_outputs/",ID.base.Model,sep=""))
  if(!file.exists(file.path(getwd(), "/Catch_MSY"))) dir.create(file.path(getwd(), "/Catch_MSY"))   
  setwd(file.path(getwd(), "/Catch_MSY"))
  Path.ktch_msy=getwd()
  
  Ktch_MSY=ktch_msy_scen
  system.time(for(sc in 1:length(ktch_msy_scen))
  {
    Folder=names(Ktch_MSY)[sc]
    if(!file.exists(paste(Path.ktch_msy,Folder,sep="/"))) dir.create(paste(Path.ktch_msy,Folder,sep="/"))   
    setwd(paste(Path.ktch_msy,Folder,sep="/"))
    
    ct=Store.Reports$`Base case`$TC_out[1:length(Yrs)]
    Current=Yrs[length(Yrs)]
    
    #forward projections
    ct.future=rep(mean(ct[(length(ct)-4):length(ct)]),years.futures)
    yr.future=Current+(1:years.futures)
    
    
    Ktch_MSY[[sc]]=Catch_MSY(ct=ct,
                             yr=Yrs,
                             r.prior=ktch_msy_scen[[sc]]$r.prior,
                             user=ktch_msy_scen[[sc]]$user,
                             k.max=ktch_msy_scen[[sc]]$k.max,
                             startbio=ktch_msy_scen[[sc]]$startbio,
                             finalbio=ktch_msy_scen[[sc]]$finalbio,
                             res=ktch_msy_scen[[sc]]$res,
                             n=ktch_msy_scen[[sc]]$niter,
                             sigR=ktch_msy_scen[[sc]]$sigR,
                             ct.future=ct.future,           
                             yr.future=yr.future)
    
    #Export outputs
    Table1_ktch_MSY=with(Ktch_MSY[[sc]],data.frame(`geom. mean r`,`r +/- 1.96 SD`,`geom. mean k (tons)`,`k +/- 1.96 SD (tons)`,
                                                   `geom. mean MSY (tons)`,`MSY +/- 1.96 SD (tons)`))
    write.csv(Table1_ktch_MSY,"Table1_ktch_MSY.csv",row.names=F)
    write.csv(cbind(Yrs=c(Yrs,yr.future),ct=c(ct,ct.future)),"ct.future.csv",row.names=F)
    
  })
  
  
  #Outputs
  setwd(Path.ktch_msy)
  
  
  #Output table of scenarios
  Tabl.scen.Ktch.MSY=vector('list',length(ktch_msy_scen))
  for(i in 1:length(ktch_msy_scen))
  {
    dummy=ktch_msy_scen[[i]]
    for(a in 1:length(ktch_msy_scen[[i]])) if(length(dummy[[a]])>1) dummy[[a]]=paste(dummy[[a]],collapse=";")
    Tabl.scen.Ktch.MSY[[i]]=unlist(dummy)
  }
  Tabl.scen.Ktch.MSY=do.call(rbind,Tabl.scen.Ktch.MSY)
  row.names(Tabl.scen.Ktch.MSY)=names(ktch_msy_scen)
  write.csv(Tabl.scen.Ktch.MSY,"Scenarios.csv")
  
  
  
  #Plot r prior dist
  fn.fig("Prior_r", 2000, 2000)
  par(las=1,mai=c(1,1.15,.1,.15),mgp=c(3.5,.75,0))
  plot(density(rgamma(10000, shape = r.pars[1], rate = r.pars[2])),
       lwd=3,main="",xlab=expression(paste(plain("Intrinsic rate of increase (year") ^ plain("-1"),")",sep="")),
       cex.lab=2,cex.axis=1.25,col=1)
  # plot(density(rlnorm(10000, meanlog = r.pars[1], sdlog = r.pars[2])),
  #      lwd=3,main="",xlab=expression(paste(plain("Intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),
  #      cex.lab=2,cex.axis=1.25,col=1)
  dev.off()
  
  
  
  
  
  #Relative biomass trend
  CL.mean="transparent"
  Yrs=c(Yrs,yr.future)
  indx.ftur=(length(Yrs)-years.futures+1):length(Yrs)
  
  Low.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (0+Nper)/100))   #get percentiles
  High.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (100-Nper)/100))
  fn.cons.po=function(low,up) c(low, tail(up, 1), rev(up), low[1])  #construct polygon
  
  #colfunc <- colorRampPalette(c("grey90","grey50"))
  colfunc <- colorRampPalette(c("aliceblue","lightblue3"))
  COLS=colfunc(3)
  colfunc.f <- colorRampPalette(c("white","burlywood3"))
  COLS.f=colfunc.f(3)
  fn.plot.percentile=function(DAT,YR,ADD.prob)
  {
    #50% of data
    Nper=(100-50)/2
    LOW.50=Low.percentile(Nper,DAT)
    UP.50=High.percentile(Nper,DAT)
    
    #75% of data
    Nper=(100-75)/2
    LOW.75=Low.percentile(Nper,DAT)
    UP.75=High.percentile(Nper,DAT)
    
    #100% of data
    Nper=(100-100)/2
    LOW.100=Low.percentile(Nper,DAT)
    UP.100=High.percentile(Nper,DAT)
    
    #construct polygons
    Year.Vec <-  fn.cons.po(YR[-indx.ftur],YR[-indx.ftur])
    Biom.Vec.50 <- fn.cons.po(LOW.50[-indx.ftur],UP.50[-indx.ftur]) 
    Biom.Vec.75 <- fn.cons.po(LOW.75[-indx.ftur],UP.75[-indx.ftur]) 
    Biom.Vec.100 <-fn.cons.po(LOW.100[-indx.ftur],UP.100[-indx.ftur]) 
    
    id.futr=c((indx.ftur[1]-1),indx.ftur)
    Year.Vec.f <-  fn.cons.po(YR[id.futr],YR[id.futr])
    Biom.Vec.50.f <- fn.cons.po(LOW.50[id.futr],UP.50[id.futr]) 
    Biom.Vec.75.f <- fn.cons.po(LOW.75[id.futr],UP.75[id.futr]) 
    Biom.Vec.100.f <-fn.cons.po(LOW.100[id.futr],UP.100[id.futr]) 
    
    
    #plot
    plot(YR,UP.100,ylim=c(0,max(UP.100)),type="l",ylab="",xlab="",xaxt='n',col='transparent',cex.axis=1.25)
    
    polygon(Year.Vec, Biom.Vec.100, col = COLS[3], border = "grey20")
    polygon(Year.Vec, Biom.Vec.75, col = COLS[2], border = "grey20")
    polygon(Year.Vec, Biom.Vec.50, col = COLS[1], border = "grey20")
    
    polygon(Year.Vec.f, Biom.Vec.100.f, col = COLS.f[3], border = "grey20")
    polygon(Year.Vec.f, Biom.Vec.75.f, col = COLS.f[2], border = "grey20")
    polygon(Year.Vec.f, Biom.Vec.50.f, col = COLS.f[1], border = "grey20")
    
    
    #add probs
    if(ADD.prob=="YES")
    {
      add.probs(id.yr=match(Current,YR),YR,DAT,UP.100,LOW.100)
      add.probs(id.yr=length(YR),YR,DAT,UP.100,LOW.100)
      
      abline(h=B.target,lwd=2,col='grey30',lty=2)
      text(YR[4],B.target,"Target",pos=3,cex=1.1)
      abline(h=B.threshold,lwd=2,col='grey30',lty=2)
      text(YR[4],B.threshold,"Threshold",pos=3,cex=1.1)
      abline(h=B.limit,lwd=2,col='grey30',lty=2)
      text(YR[4],B.limit,"Limit",pos=3,cex=1.1)
    }
    axis(1,at=YR,labels=F,tck=-0.01)
    axis(1,at=seq(YR[1],YR[length(YR)],5),labels=seq(YR[1],YR[length(YR)],5),tck=-0.02,cex.axis=1.25)
  }
  add.probs=function(id.yr,YR,DAT,UP.100,LOW.100)
  {
    f=ecdf(DAT[id.yr,])
    P.below.target=f(B.target)
    P.below.threshold=f(B.threshold)
    P.below.limit=f(B.limit)
    P.above.target=1-P.below.target
    P.above.threshold=1-P.below.threshold
    P.above.limit=1-P.below.limit
    P.between.thre.tar=P.below.target-P.below.threshold
    P.between.lim.thre=P.below.threshold-P.below.limit
    if(P.above.target>0)
    {
      segments(YR[id.yr],B.target,YR[id.yr],UP.100[id.yr],col="forestgreen",lwd=8,lend="butt")  
      text(YR[id.yr],B.target*1.5,paste(round(100*P.above.target,1),"%",sep=""),
           col="black",cex=1.1,srt=45,pos=2,font=2)
    }
    if(P.between.thre.tar>0)
    {
      segments(YR[id.yr],B.target,YR[id.yr],B.threshold,col="yellow",lwd=8,lend="butt")  
      text(YR[id.yr],mean(c(B.target,B.threshold))*1.1,paste(round(100*P.between.thre.tar,1),"%",sep=""),
           col="black",cex=1.1,srt=45,pos=2,font=2)
    }
    if(P.between.lim.thre>0)
    {
      segments(YR[id.yr],B.threshold,YR[id.yr],B.limit,col="orange",lwd=8,lend="butt")
      text(YR[id.yr],mean(c(B.threshold,B.limit))*1.1,paste(round(100*P.between.lim.thre,1),"%",sep=""),
           col="black",cex=1.1,srt=45,font=2,pos=2)
    }
    if(P.below.limit>0)
    {
      segments(YR[id.yr],B.limit,YR[id.yr],LOW.100[id.yr],col="red",lwd=8,lend="butt")
      text(YR[id.yr],B.limit*0.8,paste(round(100*P.below.limit,1),"%",sep=""),
           col="black",cex=1.1,srt=45,pos=2,font=2)
    }
  }
  
  
  #All scenarios
  fn.fig("Biomass_relative",2000,2400)
  smart.par(n.plots=length(ktch_msy_scen),MAR=c(3,4,1,1),OMA=rep(.5,4),MGP=c(1,.6,0))
  for(sc in 1:length(ktch_msy_scen))
  {
    Ktch_MSY_Rel.bio=rbind(Ktch_MSY[[sc]]$bt.rel,Ktch_MSY[[sc]]$bt.rel.future)  
    
    #Percentile   
    fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES")
    
    # #Geometric mean
    # Ktch_MSY_rel_bt_mean=Ktch_MSY_rel_bt_lowSE=Ktch_MSY_rel_bt_upSE=nrow(Ktch_MSY_Rel.bio)
    # for(nr in 1:nrow(Ktch_MSY_Rel.bio))
    # {
    #   zz=Ktch_MSY_Rel.bio[nr,]
    #   zz[zz<0]=1e-6
    #   Ktch_MSY_rel_bt_mean[nr]=exp(mean(log(zz)))
    #   Ktch_MSY_rel_bt_upSE[nr]=exp(mean(log(zz)) + 1.96 * sd(log(zz)))
    #   Ktch_MSY_rel_bt_lowSE[nr]=exp(mean(log(zz)) - 1.96 * sd(log(zz)))
    # }
    # plot(Yrs,Ktch_MSY_rel_bt_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,1),xaxt='n',xlab="",ylab="",cex.axis=1.25)
    # segments(Yrs,Ktch_MSY_rel_bt_lowSE,Yrs,Ktch_MSY_rel_bt_upSE,col=CL)
    # points(Yrs[indx.ftur],Ktch_MSY_rel_bt_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
    # segments(Yrs[indx.ftur],Ktch_MSY_rel_bt_lowSE[indx.ftur],Yrs[indx.ftur],Ktch_MSY_rel_bt_upSE[indx.ftur],col="brown4")
    # abline(h=B.target,lwd=2,col='grey30',lty=2)
    # text(Yrs[3],B.target,"Target",pos=3,cex=1.1)
    # abline(h=B.threshold,lwd=2,col='grey30',lty=2)
    # text(Yrs[5],B.threshold,"Threshold",pos=3,cex=1.1)
    # abline(h=B.limit,lwd=2,col='grey30',lty=2)
    # text(Yrs[3],B.limit,"Limit",pos=3,cex=1.1)
    # axis(1,Yrs,labels=F,tck=-0.015)
    # axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
    
    
    legend("bottomleft",names(Ktch_MSY)[sc],bty='n',cex=1.5)
  }
  mtext("Relative biomass",2,cex=1.75,las=3,line=-1,outer=T)
  mtext("Financial year",1,line=-1,cex=1.75,outer=T)
  dev.off()
  
  #Base case only
  fn.fig("Biomass_relative_Base case",2000,2400)
  par(mar=c(3,3.5,1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
  sc=match("Base case",names(ktch_msy_scen))
  Ktch_MSY_Rel.bio=rbind(Ktch_MSY[[sc]]$bt.rel,Ktch_MSY[[sc]]$bt.rel.future)  
  
  #Percentile   
  fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES")
  
  #     #Geometric mean
  Ktch_MSY_rel_bt_mean=Ktch_MSY_rel_bt_lowSE=Ktch_MSY_rel_bt_upSE=nrow(Ktch_MSY_Rel.bio)
  for(nr in 1:nrow(Ktch_MSY_Rel.bio))
  {
    Ktch_MSY_rel_bt_mean[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])))
    Ktch_MSY_rel_bt_upSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) + 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
    Ktch_MSY_rel_bt_lowSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) - 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
  }
  # plot(Yrs,Ktch_MSY_rel_bt_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,1),xaxt='n',xlab="",ylab="",cex.axis=1.25)
  # segments(Yrs,Ktch_MSY_rel_bt_lowSE,Yrs,Ktch_MSY_rel_bt_upSE,col=CL)
  # points(Yrs[indx.ftur],Ktch_MSY_rel_bt_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
  # segments(Yrs[indx.ftur],Ktch_MSY_rel_bt_lowSE[indx.ftur],Yrs[indx.ftur],Ktch_MSY_rel_bt_upSE[indx.ftur],col="brown4")
  # abline(h=B.target,lwd=2,col='grey30',lty=2)
  # text(Yrs[3],B.target,"Target",pos=3,cex=1.25)
  # abline(h=B.threshold,lwd=2,col='grey30',lty=2)
  # text(Yrs[5],B.threshold,"Threshold",pos=3,cex=1.25)
  # abline(h=B.limit,lwd=2,col='grey30',lty=2)
  # text(Yrs[3],B.limit,"Limit",pos=3,cex=1.25)
  # axis(1,Yrs,labels=F,tck=-0.015)
  # axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
  legend("bottomleft",c("50%","75%","100%"),fill=COLS,bty='n',title="Percentile",cex=1.25)
  mtext("Relative biomass",2,cex=1.75,las=3,line=-1,outer=T)
  mtext("Financial year",1,line=-1,cex=1.75,outer=T)
  dev.off()
  
  write.csv(cbind(Yrs=Yrs,geom.mean=Ktch_MSY_rel_bt_mean,
                  lowSE=Ktch_MSY_rel_bt_lowSE,UpSE=Ktch_MSY_rel_bt_upSE),"rel.biom.base.case.csv",row.names=F)
  
  
  #show example of all runs
  rain.col=colorRampPalette(c("black", "grey","white"))(ncol(Ktch_MSY_Rel.bio))
  fn.fig("Biomass_relative_all_runs",2000,2400)
  par(mar=c(3,3.5,1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
  plot(Yrs,Ktch_MSY_Rel.bio[,1],ylim=c(0,1),xaxt='n',xlab="",ylab="",type='l')
  for(nn in 2:ncol(Ktch_MSY_Rel.bio)) lines(Yrs,Ktch_MSY_Rel.bio[,nn],col=sample(rain.col,1))
  abline(h=B.target,lwd=3,col='green',lty=2)
  text(Yrs[3],B.target,"Target",pos=3,cex=1.25)
  abline(h=B.threshold,lwd=3,col='yellow',lty=2)
  text(Yrs[5],B.threshold,"Threshold",pos=3,cex=1.25)
  abline(h=B.limit,lwd=3,col='orange',lty=2)
  text(Yrs[3],B.limit,"Limit",pos=3,cex=1.25)
  axis(1,Yrs,labels=F,tck=-0.015)
  axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
  mtext("Relative biomass",2,cex=1.75,las=3,line=-1,outer=T)
  mtext("Financial year",1,line=-1,cex=1.75,outer=T)
  dev.off()
  
  
  #Current depletion of total biomass
  
  #Base case only
  fn.fig("Posterior_current_year_depletion_Base case",1200,2400)
  par(mfcol=c(3,1),mar=c(3.5,3.5,.1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
  sc=match("Base case",names(ktch_msy_scen))
  Ktch_MSY_Rel.bio=rbind(Ktch_MSY[[sc]]$bt.rel,Ktch_MSY[[sc]]$bt.rel.future)
  Ktch_MSY_current_yr=Ktch_MSY_Rel.bio[match(Current,Yrs),]
  density.fun2(what=Ktch_MSY_current_yr,B.ref=B.target,CEX=1.5) 
  legend('bottomright',"Target",bty='n',cex=2)
  density.fun2(what=Ktch_MSY_current_yr,B.ref=B.threshold,CEX=1.5) 
  legend('bottomright',"Threshold",bty='n',cex=2)
  density.fun2(what=Ktch_MSY_current_yr,B.ref=B.limit,CEX=1.5) 
  legend('bottomright',"Limit",bty='n',cex=2)
  mtext("Probability",2,cex=1.75,las=3,line=-1.75,outer=T)
  mtext(paste(Current,"relative biomass"),1,line=-1,cex=1.75,outer=T)
  dev.off()
  
  
  #export probabilities for Weight of Evidence
  write.csv(Probs.ref.point(what=Ktch_MSY_current_yr),"Consequence_likelihood_total.csv",row.names=F)
  
  #export relative biomass runs
  write.csv(Ktch_MSY_Rel.bio,"Base case/Ktch_MSY_Rel.bio.csv",row.names=F)
  
  #Future depletion of total biomass
  
  #Base case only
  Ktch_MSY_future=Ktch_MSY_Rel.bio[length(Yrs),]
  
  fn.fig("Posterior_future_year_depletion_Base case",1200,2400)
  par(mfcol=c(3,1),mar=c(3.5,3.5,.1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
  density.fun2(what=Ktch_MSY_future,B.ref=B.target,CEX=1.5) 
  legend('bottomright',"Target",bty='n',cex=2)
  density.fun2(what=Ktch_MSY_future,B.ref=B.threshold,CEX=1.5) 
  legend('bottomright',"Threshold",bty='n',cex=2)
  density.fun2(what=Ktch_MSY_future,B.ref=B.limit,CEX=1.5) 
  legend('bottomright',"Limit",bty='n',cex=2)
  mtext("Probability",2,cex=1.75,las=3,line=-1.75,outer=T)
  mtext(paste(Yrs[length(Yrs)],"relative biomass"),1,line=-1,cex=1.75,outer=T)
  dev.off()
  
  #export probabilities for Weight of Evidence
  write.csv(Probs.ref.point(what=Ktch_MSY_future),"Consequence_likelihood_total_future.csv",row.names=F)
  
  
  
  
  #Fishing mortality trend
  Yrs=Yrs[-match(yr.future,Yrs)]             
  
  #All scenarios
  fn.fig("Fishing_mortality",2000,2400)
  smart.par(n.plots=length(ktch_msy_scen),MAR=c(3,5,1,1),OMA=rep(.5,4),MGP=c(1,.6,0))
  for(sc in 1:length(ktch_msy_scen))
  {
    Fish.mort=Ktch_MSY[[sc]]$Fish.mort
    
    #Percentiles
    fn.plot.percentile(DAT=Fish.mort,YR=Yrs,ADD.prob="NO")
    
    #Geometric mean
    # Fish.mort_mean=Fish.mort_lowSE=Fish.mort_upSE=nrow(Fish.mort)
    # for(nr in 1:nrow(Fish.mort))
    # {
    #   Fish.mort_mean[nr]=exp(mean(log(Fish.mort[nr,])))
    #   Fish.mort_upSE[nr]=exp(mean(log(Fish.mort[nr,])) + 1.96 * sd(log(Fish.mort[nr,])))
    #   Fish.mort_lowSE[nr]=exp(mean(log(Fish.mort[nr,])) - 1.96 * sd(log(Fish.mort[nr,])))
    # }
    # plot(Yrs,Fish.mort_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,1),xaxt='n',xlab="",ylab="",cex.axis=1.25)
    # segments(Yrs,Fish.mort_lowSE,Yrs,Fish.mort_upSE,col=CL)
    # points(Yrs[indx.ftur],Fish.mort_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
    # segments(Yrs[indx.ftur],Fish.mort_lowSE[indx.ftur],Yrs[indx.ftur],Fish.mort_upSE[indx.ftur],col="brown4")
    # axis(1,Yrs,labels=F,tck=-0.015)
    # axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
    legend("topleft",names(Ktch_MSY)[sc],bty='n',cex=1.5)
  }
  mtext(expression(paste(plain("Fishing mortality (year") ^ plain("-1"),")",sep="")),
        2,cex=1.75,las=3,line=-2.5,outer=T)   
  mtext("Financial year",1,line=-1,cex=1.75,outer=T)
  dev.off()
  
  
  #Base case only
  fn.fig("Fishing_mortality_Base case",2000,2400)
  par(mar=c(3,5,1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
  sc=match("Base case",names(ktch_msy_scen))
  Fish.mort=Ktch_MSY[[sc]]$Fish.mort
  
  #Percentiles
  fn.plot.percentile(DAT=Fish.mort,YR=Yrs,ADD.prob="NO")
  
  #   #Get geometric mean
  # Fish.mort_mean=Fish.mort_lowSE=Fish.mort_upSE=nrow(Fish.mort)
  # for(nr in 1:nrow(Fish.mort))
  # {
  #   Fish.mort_mean[nr]=exp(mean(log(Fish.mort[nr,])))
  #   Fish.mort_upSE[nr]=exp(mean(log(Fish.mort[nr,])) + 1.96 * sd(log(Fish.mort[nr,])))
  #   Fish.mort_lowSE[nr]=exp(mean(log(Fish.mort[nr,])) - 1.96 * sd(log(Fish.mort[nr,])))
  # }
  # plot(Yrs,Fish.mort_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,max(Fish.mort_upSE)),xaxt='n',xlab="",ylab="",cex.axis=1.25)
  # segments(Yrs,Fish.mort_lowSE,Yrs,Fish.mort_upSE,col=CL)
  # points(Yrs[indx.ftur],Fish.mort_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
  # segments(Yrs[indx.ftur],Fish.mort_lowSE[indx.ftur],Yrs[indx.ftur],Fish.mort_upSE[indx.ftur],col="brown4")
  # axis(1,Yrs,labels=F,tck=-0.015)
  # axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
  
  mtext(expression(paste(plain("Fishing mortality (year") ^ plain("-1"),")",sep="")),
        2,cex=1.75,las=3,line=-2,outer=T)
  mtext("Financial year",1,line=-1,cex=1.75,outer=T)
  legend("topleft",c("50%","75%","100%"),fill=COLS,bty='n',title="Percentile",cex=1.25)
  dev.off()
}




# Section H: DROPPED CODE -------------------------------------------------
#R2admb way of running scripts
# source("C:/Matias/Analyses/SOURCE_SCRIPTS/setupADMB.R") 
# system.time(for(i in 1:nrow(Scenarios))
# {
#   d=Scenarios[i,]
#   setPath(d$Model)
#   
#   Dat=Inputs[[i]]
#   Init.pars=as.list(Pin.pars[[i]])
#   
#   
#   clean_admb(Spec)
#   fn.create.bat.file("compile.bat")   #equivalent to setup_admb()
#   compile_admb(Spec,verbose=T)
#   MOD <- do_admb(Spec,data=Dat,params=Init.pars)
#   Store.Models[[i]]=MOD
#   Store.Reports[[i]]=reptoRlist(paste(Spec,".rep",sep=""))  #store Report file
#   
# })  


# #6. Asymptotic error for base base      #not tested yet    does this do extra things than #5  ???
# i=match(ID.base.Model,Scenarios$Model)
# PARS=Store.std[[i]]
# setPath(paste("2_Outputs/Model_outputs/",Scenarios[i,]$Model,sep=""))
# 
# Tot_biom_and_SE=c("TB_WC","TB_ZN1","TB_ZN2")
# Tot_virgin_Biom_and_SE=c("Bo_WC","Bo_ZN1","Bo_ZN2")
# 
# Biom_and_SE=c("FMB_WC","FMB_ZN1","FMB_ZN2")
# Virgin_Biom_and_SE=c("Bo_mature_WC","Bo_mature_ZN1","Bo_mature_ZN2")
# F_fem_and_SE=c("F_fem_WC","F_fem_ZN1","F_fem_ZN2")
# F_male_and_SE=c("F_male_WC","F_male_ZN1","F_male_ZN2")
# 
# fn.asymptotic.SE=function(what,what2,NAME,Add.Re.Pnt,YMAX)
# {
#   dat=subset(PARS,name%in%what,select=c(value,std.dev))
#   if(!is.na(what2))
#   {
#     dat2=subset(PARS,name%in%what2,select=c(value,std.dev))
#     Cv=dat$std.dev/dat$value
#     dat$value=dat$value/dat2$value
#     dat$std.dev=Cv*dat$value
#   }
#   
#   plot(Yrs,dat$value,type='l',lwd=3,col=CL,ylab="",xlab="",xaxs="i",ylim=c(0,YMAX))
#   lines(Yrs,dat$value+1.96*dat$std.dev,lwd=1.5)
#   lines(Yrs,dat$value-1.96*dat$std.dev,lwd=1.5)
#   legend("top",NAME,bty='n',cex=2.5)
#   
#   if(Add.Re.Pnt=="YES")
#   {
#     N.YRs=length(Yrs)
#     Year.Vec <- c(Yrs, tail(Yrs, 1), rev(Yrs), Yrs[1]) 
#     B.target.Vec <- c(rep(0,N.YRs), B.target, rep(B.target,N.YRs), 0)
#     polygon(Year.Vec, B.target.Vec, col = rgb(0.01, 0.5, 0.3,alpha =0.1), border = "transparent")
#   }
#   
#   
# }
# fn.Ymax=function(what,what2)
# {
#   d=subset(PARS,name%in%what,select=c(value,std.dev))
#   e=subset(PARS,name%in%what2,select=c(value,std.dev))
#   x=d$value/e$value
#   y=d$std.dev/e$std.dev
#   return(max(c(1,max(x)*1.1)))
# }
# 
# 
# #Relative total Biomass
# fn.fig("Total relative biomass by zone",2000,2600)
# par(mfcol=c(3,1),mai=c(.1,.1,.2,.1),oma=c(4,4,1,.1),las=1,mgp=c(1,.6,0),cex.axis=1.5)
# for(i in 1:3) fn.asymptotic.SE(what=Tot_biom_and_SE[i],Tot_virgin_Biom_and_SE[i],
#                                NAME=as.character(Areas.zones$zone[i]),Add.Re.Pnt="YES",
#                                YMAX=fn.Ymax(what=Tot_biom_and_SE[i],what2=Tot_virgin_Biom_and_SE[i]))
# mtext("Relative total biomass",2,outer=T,las=3,line=1.9,cex=2)
# mtext("Financial year",1,outer=T,line=2.25,cex=2)
# dev.off()
# 
# 
# #Relative female mature biomass
# fn.fig("Female mature relative biomass by zone",2000,2600)
# par(mfcol=c(3,1),mai=c(.1,.1,.2,.1),oma=c(4,4,1,.1),las=1,mgp=c(1,.6,0),cex.axis=1.5)
# for(i in 1:3) fn.asymptotic.SE(what=Biom_and_SE[i],Virgin_Biom_and_SE[i],
#                                NAME=as.character(Areas.zones$zone[i]),Add.Re.Pnt="YES",
#                                YMAX=fn.Ymax(what=Biom_and_SE[i],what2=Virgin_Biom_and_SE[i]))
# mtext("Relative female mature biomass",2,outer=T,las=3,line=1.9,cex=2)
# mtext("Financial year",1,outer=T,line=2.25,cex=2)
# dev.off()
# 
# 
# #Fishing mortality
# fn.fig("fishing mortality by zone",2600,2600)
# par(mfcol=c(3,2),mai=c(.1,.1,.2,.3),oma=c(4,5,1,.1),las=1,mgp=c(1,.6,0),cex.axis=1.5)
# 
# #females
# YMAX=max(c(max(subset(PARS,name%in%F_fem_and_SE,select=value)+1.1*subset(PARS,name%in%F_fem_and_SE,select=std.dev)),
#            max(subset(PARS,name%in%F_male_and_SE,select=value))))
# YMAX=0.85
# for(i in 1:3)
# {
#   fn.asymptotic.SE(what=F_fem_and_SE[i],NA,
#                    NAME=as.character(Areas.zones$zone[i]),Add.Re.Pnt="NO",YMAX=YMAX)
#   if(i==1)mtext("Females",3,line=0,cex=2)
# }
# #males
# for(i in 1:3)
# {
#   fn.asymptotic.SE(what=F_male_and_SE[i],NA,
#                    NAME=as.character(Areas.zones$zone[i]),Add.Re.Pnt="NO",YMAX=YMAX)
#   if(i==1)mtext("Males",3,line=0,cex=2)
#   
# }
# mtext(expression(paste(plain("Fishing mortality (years") ^ plain("-1"),")",sep="")),2,outer=T,las=3,line=1.9,cex=2)
# mtext("Financial year",1,outer=T,line=2.25,cex=2)
# dev.off()

