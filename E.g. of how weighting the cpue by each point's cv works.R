y=1:10
y.hat=c(4,5,1,6,y[5:10]*runif(6,.9,1.1))
std=rep(0.1,length(y))
std.CV=c(4,3,2,1,rep(.1,6))

par(mfcol=c(2,1),mai=c(.6,.6,.2,.1))

plot(y,y.hat)
lines(y,y,col=2)
segments(y,y.hat,y,(y.hat+std))
segments(y,y.hat,y,(y.hat-std))

plot(y,y.hat)
lines(y,y,col=2)
segments(y,y.hat,y,(y.hat+std.CV))
segments(y,y.hat,y,(y.hat-std.CV))


epsilon=(log(y)-log(y.hat))^2
#fn.neg.log.like=function(SD,EPSLN) log(SD)+(EPSLN/(2.*SD*SD))
#fn.neg.log.like=function(SD,OBS,PRED) (((log(OBS)-log(PRED))^2)/(2*SD))+(log(2*pi*SD))/2
fn.neg.log.like=function(SD,OBS,PRED) log(SD)+0.5*(log(OBS/PRED)/SD)^2


nobs=length(y)
LL=LL.CV=rep(NA,nobs)
for(i in 1:nobs)
{
  LL[i]=fn.neg.log.like(std[i],y[i],y.hat[i])
  LL.CV[i]=fn.neg.log.like(std.CV[i],y[i],y.hat[i])
}
A=data.frame(std.CV=std.CV,LL=round(LL,3),LL.largeCV=round(LL.CV,3))
A
sum(LL)
sum(LL.CV)
#Hence, sum of LL.CV is lower because the model is not forced to pass thru each point as the error bars are 
# touching the red line (prediction). For LL, all SD are very small so the prediction is forced to pass thru
# the points, hence giving a higher neg log like value
