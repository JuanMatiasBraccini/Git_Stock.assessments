
#Age length key for commercial shark species

#note: this combines data from different years (i.e. assumes no recruitment pulses and constant mortality)


library(FSA)

#DATA
setwd("C:/Matias/Data/Age and growth")
whiskery=read.csv("Whiskery.csv")


#Convert TL to FL McAuley unpublished (TL to FL)
a.w=1.0044
b.w=13.171
TL_to_FL=function(TL,a,b) FL=(TL-b)/a
whiskery$FL=with(whiskery,ifelse(is.na(FL)& !is.na(TL),TL_to_FL(TL,a.w,b.w),FL))
names(whiskery)[match("Counts",names(whiskery))]="age"


#Age length key
fn.age.key=function(dat)
{
  # add length categories
  WR.age.mod <- lencat(~FL,data=dat,startcat=min(dat$FL),w=10)
  # create age-length key
  raw <- table(WR.age.mod$LCat, WR.age.mod$age)
  WR.key <- round(prop.table(raw, margin=1),2)
  return(WR.key)
}

Whiskery.key=fn.age.key(Subset(whiskery, !is.na(age)))



#2. Calculate age composition from length composition sample

Dummy.Whis.FL=data.frame(FL=round(runif(1000,min(whiskery$FL),max(whiskery$FL))))
Dummy.Whis.FL <- lencat(~FL,data=Dummy.Whis.FL,startcat=min(whiskery$FL),w=10)
Tbl.size=table(Dummy.Whis.FL$LCat)

Age.matrix=c(Tbl.size)*Whiskery.key
Age.comp=round(colSums(Age.matrix))


