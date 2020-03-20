library(FSA)
library(fitdistrplus)
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


