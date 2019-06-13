library(fishmethods)
library(lubridate)
data(tanaka)

#---- Estimation of Mortality using Times-At-Large Data from Tagging

# relyr= a vector of release year (or cohort) for individual times-at-large observations.
# tal= a vector of individual times-at-large observations.
# N= a vector of number of releases for each release year (or cohort). Each individual
#   observation from a release year should have the same N value.
#   method 1 = McGarvey et al., 2 = Gulland. Default is all (i.e., c(1,2)).
# np= the number of periods over which to combine data to make period estimates of
#   mortality. Set np=0 to estimate mortality for each release year.
# stper= vector of year values representing the beginning year of each period over which
#   to estimate mortality. The first year in c() must always be the first release year.
# nboot= the number of resamples for the Gulland method

#Mortality=mort.al(relyr = tanaka$relyr, tal = tanaka$tal, N = tanaka$N)




Dat=read.csv("C:/Matias/Analyses/Mark_recapture_analysis/Tagging.data.csv",stringsAsFactors=T)
Dat=subset(Dat,select=c(SHEET_NO,Species,Sex,Day.rel,Mn.rel,Yr.rel,Day.rec,Mn.rec,Yr.rec))

SPECS=unique(Dat$Species)
STORE.mortality=vector('list',length(SPECS))
names(STORE.mortality)=SPECS
fn.Mrt=function(sp)
{
  a=subset(Dat,Species==sp)
  a$N=1
  Yr.rel.N=aggregate(N~Yr.rel,a,sum)
  a=a[,-match("N",names(a))]
  
  Recs=subset(a,!is.na(Yr.rec))
  Recs=subset(a,!is.na(Mn.rec))
  Recs=subset(a,!is.na(Day.rec))

  Recs$Rel.date=with(Recs,as.POSIXct(paste(Yr.rel,"-",Mn.rel,"-",Day.rel,sep="")))
  Recs$Rec.date=with(Recs,as.POSIXct(paste(Yr.rec,"-",Mn.rec,"-",Day.rec,sep="")))
  Recs$tal=as.numeric(difftime(Recs$Rec.date, Recs$Rel.date, units="days"))/365
  Recs=subset(Recs,tal>0)
  
  Recs=merge(Recs,Yr.rel.N,by="Yr.rel")
  
  Mortality=mort.al(relyr = Recs$Yr.rel, tal = Recs$tal, N = Recs$N)
  return(Mortality)
}
for(i in 1:length(SPECS))STORE.mortality[[i]]=fn.Mrt(sp=SPECS[i])




#---- Estimate of Population Size from a Single Mark-Recapture Experiment
mrN.single(M=948,C=421,R=167)
mrN.single(M=9480,C=421,R=167)
