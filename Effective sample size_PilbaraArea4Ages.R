
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

# dat=read.csv("C:\\Ainslie\\Finfish\\Pilbara\\2015\\Peter\\PilbaraArea4Ages.csv", header=FALSE)
# dat=read.csv("C:\\Ainslie\\Finfish\\Pilbara\\2015\\Peter\\PilbaraArea4Ages2.csv", header=FALSE)
dat=read.csv(handl_OneDrive("Analyses/Effective_sample_size/PilbaraArea4Ages2.csv"), header=FALSE)
dat

ndat=as.matrix(dat[, -c(1:3)])

xx=seq(1, nrow(dat), 2)

obs=ndat[xx,]
pred=ndat[xx+1,]
Nassumed=dat[xx, 2]
years=dat[xx, 1]

##cjheck prop sum to 1
obs=obs/rowSums(obs)
pred=pred/rowSums(pred)

str(obs)
ages=1:25

##plot pred versus obs proportions
par(mfrow=c(3,2))
for (ii in 1:dim(obs)[1]){
  plot(ages, obs[ii,], "o", ylim=c(0, 0.25), xlab="age", ylab="Prop")
  lines(ages, pred[ii,], "o", col=2)
  mtext(years[ii], side=3, line=-1.2, adj=0.99)
  mtext(paste("N=", Nassumed[ii], sep=""), side=3, line=-2.2, adj=0.99)
}


##Effective Sample Size

## (1) mcallister and ianelli
ii=5
sum(pred[ii, ]*(1-pred[ii,]))/sum((obs[ii, ]-pred[ii, ])^2)
effN_MI=rep(NA, length(years))
par(mfrow=c(3,2))
for (ii in 1:dim(obs)[1]){
  plot(ages, obs[ii,], "o", ylim=c(0, 0.25), xlab="age", ylab="Prop")
  lines(ages, pred[ii,], "o", col=2)
  mtext(years[ii], side=3, line=-1.2, adj=0.99, cex=0.8)
  mtext(paste("N=", Nassumed[ii], sep=""), side=3, line=-2.2, adj=0.99, cex=0.8)
  effN=sum(pred[ii, ]*(1-pred[ii,]))/sum((obs[ii, ]-pred[ii, ])^2)
  mtext(paste("effN(MI)=", round(effN, 1), sep=""), side=3, line=-3.2, adj=0.99, cex=0.8)
  effN_MI[ii]=effN
}
mean(effN_MI)
effN_MI/Nassumed


xx=c(1, 3, 5, 7, 9)

obs=ndat[xx,]
pred=ndat[xx+1,]
Nassumed=dat[xx, 2]
years=dat[xx, 1]

##cjheck prop sum to 1
obs=obs/rowSums(obs)
pred=pred/rowSums(pred)

Nassumed

## francis TA1.8
My <- cbind(Obs=apply(obs,1,function(x)sum(ages*x)),
            Exp=apply(pred,1,function(x)sum(ages*x)))
Ry <- My[,'Obs']-My[,'Exp']
Sy <- sqrt(apply(pred,1,function(x)sum(x*ages^2))-My[,'Exp']^2)

var(Ry*sqrt(Nassumed)/Sy)
wj <- 1/var(Ry*sqrt(Nassumed)/Sy)
wj

1/wj

Nassumed*wj


par(mfrow=c(3,2))
for (ii in 1:dim(obs)[1]){
  plot(ages, obs[ii,], "o", ylim=c(0, 0.25), xlab="age", ylab="Prop")
  lines(ages, pred[ii,], "o", col=2)
  mtext(years[ii], side=3, line=-1.2, adj=0.99, cex=0.8)
  mtext(paste("N=", Nassumed[ii], sep=""), side=3, line=-2.2, adj=0.99, cex=0.8)
  effN=Nassumed[ii]*wj
  mtext(paste("effN(TA1.8)=", round(effN, 1), sep=""), side=3, line=-3.2, adj=0.99, cex=0.8)
}


effective.n=function(obs,pred,bin,Nassumed)
{
  My <- cbind(Obs=apply(obs,1,function(x)sum(bin*x)),
              Exp=apply(pred,1,function(x)sum(bin*x)))
  Ry <- My[,'Obs']-My[,'Exp']
  Sy <- sqrt(apply(pred,1,function(x)sum(x*bin^2))-My[,'Exp']^2)
  wj <- 1/var(Ry*sqrt(Nassumed)/Sy)
  return(wj)
}
