set.seed=666

fn.mod=function(N.shk.start,Yrs.mdld,Number.tagged.yr,Pop.grwz.r8t)
{
  YRs=1:Yrs.mdld
  N.shk.end=N.shk.start+round(N.shk.start*Pop.grwz.r8t*Yrs.mdld)
  N=vector('list',length=Yrs.mdld)
  Tagged=Recaptured=N
  
  #year 1
  N[[1]]=1:N.shk.start
  Tagged[[1]]=sample(N[[1]],round(runif(1,Number.tagged.yr[1],Number.tagged.yr[2])))
  names(N[[1]])=rep(1,length(N[[1]]))
  names(Tagged[[1]])=rep(1,length(Tagged[[1]]))
  
  #other years
  for(y in YRs[-1])
  {
    #grow population
    lastN=N[[y-1]][length(N[[y-1]])]
    N[[y]]=c(N[[y-1]],(lastN+1):(lastN+lastN*Pop.grwz.r8t))
    names(N[[y]])=rep(y,length(N[[y]]))
    
    #tagging event
    Tag=sample(N[[y]],round(runif(1,Number.tagged.yr[1],Number.tagged.yr[2])))
    rec.tag=match(unlist(Tagged),Tag)  
    if(length(rec.tag[!is.na(rec.tag)])>0)
    {
      Recaptured[[y]]=Tag[which(Tag%in%unlist(Tagged))]
      names(Recaptured[[y]])=rep(y,length(Recaptured[[y]]))
    }
    
    Tagged[[y]]=Tag[which(!Tag%in%unlist(Tagged))]  
    names(Tagged[[y]])=rep(y,length(Tagged[[y]]))
  }

  #Put as data frame
  Dummy.dat=data.frame(Year=names(unlist(Tagged)),Shark.ID=unlist(Tagged))
  Dummy.dat$Recapture=0
  dd=data.frame(Year=names(unlist(Recaptured)),Shark.ID=unlist(Recaptured))
  dd$Recapture=1
  Dummy.dat=subset(Dummy.dat,!Shark.ID%in%dd$Shark.ID)
  Dummy.dat=rbind(Dummy.dat,dd)
  
  
  Dummy.dat$Year=2008+as.numeric(as.character(Dummy.dat$Year))
  Dummy.dat$DDate=as.Date(paste(Dummy.dat$Year,'/01/01',sep=''))
  Dummy.dat$DDate1=as.Date(paste(Dummy.dat$Year+1,'/01/01',sep=''))
  Dummy.dat$Date=as.Date('1999/01/01')
  for(d in 1:nrow(Dummy.dat))Dummy.dat$Date[d]=sample(seq(Dummy.dat$DDate[d], Dummy.dat$DDate1[d], by="day"), 1)
  Dummy.dat=Dummy.dat[order(Dummy.dat$Date),match(c("Date","Shark.ID","Recapture"),names(Dummy.dat))]
  
  
  #plot population growth
  Pop.grow=data.frame(Yr=rep(NA,length(N)),N=rep(NA,length(N)))
  for(p in 1:nrow(Pop.grow))
  {
    Pop.grow$Yr[p]=p
    Pop.grow$N[p]=N[[p]][length(N[[p]])]
  }
  plot(Pop.grow$Yr,Pop.grow$N,type='b',ylab="Number of individuals",xlab="Year")
  return(Dummy.dat)
}

D10000=fn.mod(N.shk.start=10000,Yrs.mdld=10,Number.tagged.yr=c(20,30),Pop.grwz.r8t=0.05)
D1000=fn.mod(N.shk.start=1000,Yrs.mdld=10,Number.tagged.yr=c(10,20),Pop.grwz.r8t=0.05)
D100=fn.mod(N.shk.start=100,Yrs.mdld=10,Number.tagged.yr=c(1,10),Pop.grwz.r8t=0.05)
D100_10percent_20yrs=fn.mod(N.shk.start=100,Yrs.mdld=20,Number.tagged.yr=c(1,10),Pop.grwz.r8t=0.1)
sort(table(D100_10percent_20yrs$Shark.ID))

sort(table(D100$Shark.ID))

sum(D100$Recapture)
sum(D1000$Recapture)
sum(D10000$Recapture)

write.csv(D10000,"C:/Matias/Analyses/Population dynamics/White shark/Mark_recapture/dummy1.csv",row.names=F)
write.csv(D1000,"C:/Matias/Analyses/Population dynamics/White shark/Mark_recapture/dummy2.csv",row.names=F)
write.csv(D100,"C:/Matias/Analyses/Population dynamics/White shark/Mark_recapture/dummy3.csv",row.names=F)

write.csv(D100_10percent_20yrs[!duplicated(D100_10percent_20yrs$Shark.ID),],"C:/Matias/Analyses/Population dynamics/White shark/Mark_recapture/dummy4.csv",row.names=F)


