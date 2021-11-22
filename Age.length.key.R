library(readxl)
library(tidyverse)
library(fishmethods)
library(ggpubr)
library(FSA)
library(fitdistrplus)
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')


# Data section ------------------------------------------------------------
setwd(handl_OneDrive("Data/Age and growth"))
#Dusky.age <- read_excel("Dusky_fromRory.xlsx", sheet = "data_age_consensus")  #useless
#Dusky.size <- read_excel("Dusky_fromRory.xlsx", sheet = "data_size")

Dusky=read_excel("Dusky_fromColin.xlsx", sheet = "Sheet1")
Sandbar=read.csv('Sandbar_fromRory.csv')
Whiskery=read.csv('Whiskery.csv')
Gummy=read.csv('Terry/Gummy_Terry.csv')

Gummy.size.birth.TL=c(30,36) #Total length
Whiskery.size.birth.TL=c(22,27) #Total length

Bin.size=5

# Manipulations ------------------------------------------------------------
Dusky=Dusky%>%
  rename(Age='new final')%>%
  mutate(Sex=ifelse(SEX=='F','Female',ifelse(SEX=='M',"Male",NA)))%>%
  filter(!is.na(Age))%>%
  dplyr::select(FL,Sex,Age)

Sandbar=Sandbar%>%
  mutate(SEX=tolower(SEX),
         Age=All.reader.consensus.count,
         drop=ifelse(Age<7 & FL>125,'Yes','No'))%>%
  filter(drop=='No' & FL>10)%>%
  mutate(Sex=ifelse(SEX=='f','Female',ifelse(SEX=='m',"Male",NA)))%>%
  dplyr::select(FL,Sex,Age)

Whiskery=Whiskery%>%
  mutate(Age=Counts)%>%
  filter(Source=="Simpfendorfer et al 2000")%>%
  dplyr::select(FL,Sex,Age)

Gummy=Gummy%>%
  mutate(Age=Counts,
         FL=ifelse(is.na(FL),(TL-3.8499)/1.0952846,FL))%>%
  dplyr::select(FL,Sex,Age)

#put in list
Data=list(Dusky=Dusky,Sandbar=Sandbar,Whiskery=Whiskery,Gummy=Gummy)
names(Data)=paste(names(Data),'shark')

# Export original age-length data ------------------------------------------------------------
for(s in 1:length(Data))
{
  d=Data[[s]]%>%
    filter(Sex%in%c('Male','Female'))
  write.csv(d,handl_OneDrive(paste(
    paste("Analyses/Data_outs/",names(Data)[s],sep=''),
    paste(names(Data)[s],"age_length.csv",sep='.'),
    sep='/')),row.names = F)
}

# Add size at birth for gummy and whiskery ------------------------------------------------------------
Gummy.size.birth.FL=ceiling((Gummy.size.birth.TL-3.8499)/1.0952846)
Whiskery.size.birth.FL=ceiling((Whiskery.size.birth.TL-7.818082)/1.0700731)

add.Gummy=round(runif(10,Gummy.size.birth.FL[1],Gummy.size.birth.FL[2]))
add.Whiskery=round(runif(10,Whiskery.size.birth.FL[1],Whiskery.size.birth.FL[2]))

Gummy1=Gummy[1:length(add.Gummy),]%>%
  mutate(Sex=NA,
         FL=add.Gummy,
         Age=0)
Whiskery1=Whiskery[1:length(add.Whiskery),]%>%
  mutate(Sex=NA,
         FL=add.Whiskery,
         Age=0)

Gummy=rbind(Gummy,Gummy1)  
Whiskery=rbind(Whiskery,Whiskery1) 

#put in list
Data=list(Dusky=Dusky,Sandbar=Sandbar,Whiskery=Whiskery,Gummy=Gummy)
names(Data)=paste(names(Data),'shark')

# Plot age len ------------------------------------------------------------
fn.plt.age.len=function(d,TIT)
{
  d%>%
    ggplot(aes(Age,FL,colour=Sex))+
    geom_point(size=3)+
    ylim(0,max(d$FL))+
    ggtitle(TIT)
}
p=Data
for(s in 1:length(Data))p[[s]]=fn.plt.age.len(d=Data[[s]],TIT=names(Data)[s])
ggarrange(plotlist=p, ncol = 2, nrow = 2)


# Create age-length-key ------------------------------------------------------------
age.len.ki=function(Age,Size,Bin)
{
  Numbers.at.age=with(pinfish,alk(age=round(Age,0),size=Size,binsize=Bin))
  Proportions.at.age=with(pinfish,alk(age=round(Age,0),size=Size,binsize=Bin,type=2))
  Total.n.at.size=Numbers.at.age%>%dplyr::select(len,nl)
  Numbers.at.age=Numbers.at.age%>%dplyr::select(-nl)
  Proportions.at.age=Proportions.at.age%>%dplyr::select(-nl)
  
  return(list(Numbers.at.age=Numbers.at.age,
              Proportions.at.age=Proportions.at.age,
              Total.n.at.size=Total.n.at.size))
}
Age.LK=Data
for(s in 1:length(Data))Age.LK[[s]]=age.len.ki(Age=Data[[s]]$Age,
                                               Size=Data[[s]]$FL,
                                               Bin=Bin.size)

#plot
for(s in 1:length(Data))
{
  Size=Bin.size*round(Data[[s]]$FL/Bin.size)
  Age=round(Data[[s]]$Age)
  plot.key <- prop.table(xtabs(~Size+Age), margin=1)
  plot.key=alkPlot(plot.key,"bubble",col=col2rgbt("red",0.5),main=names(Data)[s])
}




# Export age-length-key ------------------------------------------------------------
for(s in 1:length(Age.LK))
{
  d=Age.LK[[s]]$Proportions.at.age
  write.csv(d,handl_OneDrive(paste(
    paste("Analyses/Data_outs/",names(Data)[s],sep=''),
    paste(names(Data)[s],"age_length_key.csv",sep='.'),
    sep='/')),row.names = F)
}


# Previous approach ------------------------------------------------------------
get.prop.at.age.from.length=function(age,mn.len,SD,N,int,Obs.len,min.obs)
{
  #create sample of random lengths  
  Key=data.frame(age=rep(age,each=N),
                 Mean.len=rep(mn.len,each=N),
                 SD=rep(SD,each=N))%>%
    mutate(len=round(rnorm(n(),Mean.len,SD)),
           LCat=lencat(len,startcat=int*floor(min(len)/int),w=int))
  
  #create age-length key
  raw <- table(Key$LCat, Key$age)
  WR.key <- prop.table(raw, margin=1)
  
  #Get age from size sample and selectivity at age
  if(length(Obs.len)>=min.obs)
  {
    #age from size
    WR1.len=data.frame(age=as.integer(NA),len=Obs.len)
    suppressWarnings(WR1.len <- alkIndivAge(WR.key,age~len,data=WR1.len))
    
    #selectivity at age
    n.class=length(age)
    #gamma
    fg=fitdist(WR1.len$age[WR1.len$age>0], "gamma")
    Sel=curve(dgamma(x,shape = fg$estimate[1], rate = fg$estimate[2] ),from=age[1],to=age[length(age)],n=n.class)
    
    #lognormal
    #fln <- fitdist(WR1.len$age, "lnorm")
    #Sel=curve(dlnorm(x,meanlog = fln$estimate[1], sdlog = fln$estimate[2] ),from=age[1],to=age[length(age)],n=n.class)
    
    Sel$y=Sel$y/max(Sel$y)
    
    return(list(dat=Key,age.len.key=WR.key,pred.age=WR1.len,Selectivity=Sel))
  }
  
}

# #Age length key for commercial shark species
# 
# #note: this combines data from different years (i.e. assumes no recruitment pulses and constant mortality)
# 
# #DATA
# setwd("C:/Matias/Data/Age and growth")
# whiskery=read.csv("Whiskery.csv")
# 
# 
# #Convert TL to FL McAuley unpublished (TL to FL)
# a.w=1.0044
# b.w=13.171
# TL_to_FL=function(TL,a,b) FL=(TL-b)/a
# whiskery$FL=with(whiskery,ifelse(is.na(FL)& !is.na(TL),TL_to_FL(TL,a.w,b.w),FL))
# names(whiskery)[match("Counts",names(whiskery))]="age"
# 
# 
# #Age length key
# fn.age.key=function(dat)
# {
#   # add length categories
#   WR.age.mod <- lencat(~FL,data=dat,startcat=min(dat$FL),w=10)
#   # create age-length key
#   raw <- table(WR.age.mod$LCat, WR.age.mod$age)
#   WR.key <- round(prop.table(raw, margin=1),2)
#   return(WR.key)
# }
# 
# Whiskery.key=fn.age.key(Subset(whiskery, !is.na(age)))
# 
# 
# 
# #2. Calculate age composition from length composition sample
# 
# Dummy.Whis.FL=data.frame(FL=round(runif(1000,min(whiskery$FL),max(whiskery$FL))))
# Dummy.Whis.FL <- lencat(~FL,data=Dummy.Whis.FL,startcat=min(whiskery$FL),w=10)
# Tbl.size=table(Dummy.Whis.FL$LCat)
# 
# Age.matrix=c(Tbl.size)*Whiskery.key
# Age.comp=round(colSums(Age.matrix))


