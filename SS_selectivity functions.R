# selectivity functions

#' Calculate values for length logistic selectivity
#'
#' @param len A vector of lengths or ages
#' @param a The inflection point
#' @param b The 95% width
#' @return The selectivity curve as a vector
logistic1.fn <- function(len, a, b) {
  neglog19 <- -1 * log(19)
  denom <- 1. + exp(neglog19 * (len - a) / b)
  sel <- 1 / denom
  return(sel)
}

#' Calculate values for double normal selectivity
#'
#' @param x A vector of lengths or ages
#' @param a The peak
#' @param b The top
#' @param c The ascending width
#' @param d The Descending width
#' @param e The initial value
#' @param f The final value
#' @param use_e_999 Is -999 used for the initial value?
#' @param use_f_999 Is -999 used for the final value?
#' @return The double normal selectivity curve given the parameters as a vector
doubleNorm24.fn <- function(x, a, b, c, d, e, f, use_e_999, use_f_999) {
  # TODO: check if function not handling f < -1000 correctly (and -999 vals)
  if (use_e_999) {
    e <- -999
  }
  if (use_f_999) {
    f <- -999
  }
  if (e == 0) { # Avoid errors on the bounds
    e <- 1 - 0.999955 # an input that results in approx -10
  }
  if (e == 1) {
    e <- 0.999955 # an input that results in approx 10
  }
  if (e > 0) {
    e <- log(e / (1 - e)) # transform input to logit
  }
  
  if (f == 0) { # Avoid errors on the bounds
    f <- 1 - 0.999955 # an input that results in approx -10
  }
  if (f == 1) {
    f <- 0.999955 # an input that results in approx 10
  }
  if (f > 0) {
    f <- log(f / (1 - f)) # transform input to logit
  }
  sel <- rep(NA, length(x))
  startbin <- 1
  peak <- a
  upselex <- exp(c)
  downselex <- exp(d)
  final <- f
  if (e < -1000) {
    j1 <- -1001 - round(e)
    sel[1:j1] <- 1e-06
  }
  if (e >= -1000) {
    j1 <- startbin - 1
    if (e > -999) {
      point1 <- 1 / (1 + exp(-e))
      t1min <- exp(-(x[startbin] - peak)^2 / upselex)
    }
  }
  if (f < -1000) {
    j2 <- -1000 - round(f)
  }
  if (f >= -1000) {
    j2 <- length(x)
  }
  bin_width <- x[2] - x[1]
  peak2 <- peak + bin_width + (0.99 * x[j2] - peak - bin_width) / (1 +
                                                                     exp(-b))
  if (f > -999) {
    point2 <- 1 / (1 + exp(-final))
    t2min <- exp(-(x[j2] - peak2)^2 / downselex)
  }
  t1 <- x - peak
  t2 <- x - peak2
  join1 <- 1 / (1 + exp(-(20 / (1 + abs(t1))) * t1))
  join2 <- 1 / (1 + exp(-(20 / (1 + abs(t2))) * t2))
  if (e > -999) {
    asc <- point1 + (1 - point1) * (exp(-t1^2 / upselex) -
                                      t1min) / (1 - t1min)
  }
  if (e <= -999) {
    asc <- exp(-t1^2 / upselex)
  }
  if (f > -999) {
    dsc <- 1 + (point2 - 1) * (exp(-t2^2 / downselex) -
                                 1) / (t2min - 1)
  }
  if (f <= -999) {
    dsc <- exp(-(t2)^2 / downselex)
  }
  idx.seq <- (j1 + 1):j2
  sel[idx.seq] <- asc[idx.seq] * (1 - join1[idx.seq]) + join1[idx.seq] * (1 -
                                                                            join2[idx.seq] + dsc[idx.seq] * join2[idx.seq])
  if (startbin > 1 && e >= -1000) {
    sel[1:startbin] <- (x[1:startbin] / x[startbin])^2 *
      sel[startbin]
  }
  if (j2 < length(x)) {
    sel[(j2 + 1):length(x)] <- sel[j2]
  }
  return(sel)
}

wrapper.fn=function(x,a.dn,b.dn,c.dn,d.dn,e.dn,f.dn,a.log,b.log,out.logis=TRUE,out.DN=TRUE,XLAB='size')
{
  Sel.dn=doubleNorm24.fn(x,a=a.dn,b=b.dn, c=c.dn, d=d.dn, e=e.dn, f=f.dn,use_e_999=FALSE, use_f_999=FALSE)
  Sel.log=logistic1.fn(len=x,a=a.log, b=b.log)
  Col.DN='black'
  Col.logis='red'
  Lgn.DN=paste0('double normal (p1=',round(a.dn,1),', p2=',round(b.dn,1),', p3=',round(c.dn,1),
                ', p4=',round(d.dn,1),', p5=',round(e.dn,1),', p6=',round(f.dn,1),')')
  Lgn.Logis=paste0('logistic(p1=',round(a.log,1),', p2=',round(b.log,1),')')
  if(!out.DN)
  {
    Col.DN='transparent'
    Lgn.DN=''
  }
    
  if(!out.logis)
  {
    Col.logis='transparent'
    Lgn.Logis=''
  }
  par(xpd=T)
  plot(x,Sel.dn,ylab='Selectivity',xlab=XLAB,type='l',col=Col.DN)
  lines(x,Sel.log,col=Col.logis)
  legend('topleft',c(Lgn.DN,Lgn.Logis), inset=0.05,lty=1,col=c(Col.DN,Col.logis),bty='n',cex=0.75)
  
}
#wrapper.fn(x=40:160,a.dn=70,b.dn=-7,c.dn=6,d.dn=7,e.dn=-999,f.dn=-999,a.log=325,b.log=20,out.logis=FALSE)



ret.fn=function(L,p1,p2,p3,p4,p5,p6,p7)
{
  #logistic retention
  if(is.null(p5))
  {
    Retention=(p3/(1+exp(-(L-(p1+p4))/p2)))
    Tit=paste0("p1=",p1," p2=",p2," p3=",p3," p4=",p4)
  }
  #dome-shaped retention
  if(!is.null(p5))
  {
    Retention=(p3/(1+exp(-(L-(p1+p4))/p2)))*(1-(1/(1+exp(-(L-(p5+p7))/p6))))
    Tit=paste0("p1=",p1," p2=",p2," p3=",p3," p4=",p4," p5=",p5," p6=",p6," p7=",p7)
  }
  plot(L,Retention,type='l',main=Tit,ylim=c(0,1))
}
#ret.fn(L=0:200,p1=102,p2=2,p3=0.8,p4=0,p5=NULL,p6=NULL,p7=NULL)
#ret.fn(L=0:200,p1=87,p2=1,p3=1,p4=0,p5=NULL,p6=NULL,p7=NULL)
#ret.fn(L=90:400,p1=80,p2=1,p3=1,p4=0,p5=190,p6=1,p7=0)

mort.fn=function(L,p1,p2,p3,p4)
{
  Discard.mortality=1-((1-p3)/(1+exp((-(L-(p1+p4)))/p2)))
  Tit=paste0("p1=",p1," p2=",p2," p3=",p3," p4=",p4)
  plot(L,Discard.mortality,type='l',main=Tit,ylim=c(0,1))
}
#mort.fn(L=0:200,p1=20,p2=1.5,p3=.05,p4=0)  
#mort.fn(L=0:200,p1=5,p2=0,p3=.1,p4=0)  


doubleNorm24.fn_fit <- function(x, a, b, c, d) {

  e=-999
  f=-999
  if (e == 0) { # Avoid errors on the bounds
    e <- 1 - 0.999955 # an input that results in approx -10
  }
  if (e == 1) {
    e <- 0.999955 # an input that results in approx 10
  }
  if (e > 0) {
    e <- log(e / (1 - e)) # transform input to logit
  }
  
  if (f == 0) { # Avoid errors on the bounds
    f <- 1 - 0.999955 # an input that results in approx -10
  }
  if (f == 1) {
    f <- 0.999955 # an input that results in approx 10
  }
  if (f > 0) {
    f <- log(f / (1 - f)) # transform input to logit
  }
  sel <- rep(NA, length(x))
  startbin <- 1
  peak <- a
  upselex <- exp(c)
  downselex <- exp(d)
  final <- f
  if (e < -1000) {
    j1 <- -1001 - round(e)
    sel[1:j1] <- 1e-06
  }
  if (e >= -1000) {
    j1 <- startbin - 1
    if (e > -999) {
      point1 <- 1 / (1 + exp(-e))
      t1min <- exp(-(x[startbin] - peak)^2 / upselex)
    }
  }
  if (f < -1000) {
    j2 <- -1000 - round(f)
  }
  if (f >= -1000) {
    j2 <- length(x)
  }
  bin_width <- x[2] - x[1]
  peak2 <- peak + bin_width + (0.99 * x[j2] - peak - bin_width) / (1 +
                                                                     exp(-b))
  if (f > -999) {
    point2 <- 1 / (1 + exp(-final))
    t2min <- exp(-(x[j2] - peak2)^2 / downselex)
  }
  t1 <- x - peak
  t2 <- x - peak2
  join1 <- 1 / (1 + exp(-(20 / (1 + abs(t1))) * t1))
  join2 <- 1 / (1 + exp(-(20 / (1 + abs(t2))) * t2))
  if (e > -999) {
    asc <- point1 + (1 - point1) * (exp(-t1^2 / upselex) -
                                      t1min) / (1 - t1min)
  }
  if (e <= -999) {
    asc <- exp(-t1^2 / upselex)
  }
  if (f > -999) {
    dsc <- 1 + (point2 - 1) * (exp(-t2^2 / downselex) -
                                 1) / (t2min - 1)
  }
  if (f <= -999) {
    dsc <- exp(-(t2)^2 / downselex)
  }
  idx.seq <- (j1 + 1):j2
  sel[idx.seq] <- asc[idx.seq] * (1 - join1[idx.seq]) + join1[idx.seq] * (1 -
                                                                            join2[idx.seq] + dsc[idx.seq] * join2[idx.seq])
  if (startbin > 1 && e >= -1000) {
    sel[1:startbin] <- (x[1:startbin] / x[startbin])^2 *
      sel[startbin]
  }
  if (j2 < length(x)) {
    sel[(j2 + 1):length(x)] <- sel[j2]
  }
  return(sel)
}


fun.compare.sel.obs.size.comp=function(TL,SEL,SEl.sens_NSF,size.comps,Flts)
{
  smart.par(length(Flts),c(1,1.5,1,1),c(1,1,1,1),c(3, 1, 0))
  for(x in 1:length(Flts))
  {
    FLiT=Flts[x]
    if(FLiT=="Pilbara_Trawl") FLiT='Other'
    Sel=SEL%>%filter(Fleet==FLiT)
    Sel.ori=SS_selectivity_init_pars%>%filter(Species==names(List.sp)[i])
    a.log=1e5
    b.log=0
    OUT.LOGIS=FALSE
    OUT.DOUBLE.N=TRUE
    if(all(is.na(Sel[,c('P_3', 'P_4', 'P_5', 'P_6')])))
    {
      a.log=Sel$P_1
      b.log=Sel$P_2
      Sel[,c('P_1', 'P_2','P_3', 'P_4', 'P_5', 'P_6')]=c(10,100,-10,-10,0,0) 
      OUT.LOGIS=TRUE
      OUT.DOUBLE.N=FALSE
      
    }
    with(Sel,wrapper.fn(x=TL,a.dn=P_1,b.dn=P_2,c.dn=P_3,d.dn=P_4,e.dn=P_5,f.dn=P_6,
                        a.log=a.log,b.log=b.log,out.logis=OUT.LOGIS,out.DN=OUT.DOUBLE.N,XLAB='TL (cm)'))
    if(FLiT%in%c('Northern.shark','Survey') & !is.null(SEl.sens_NSF))
    {
      par=SEl.sens_NSF%>%filter(Fleet==FLiT)
      lines(TL,logistic1.fn(TL,par$P_1,par$P_2),lty=2,col='red')
      text(min(TL)+1,.05,'alternative',col='red',pos=4)
    }
    dd=size.comps%>%
      filter(Fleet==Flts[x])
    if(nrow(dd)>2)
    {
      dd$bin=TL.bins.cm*floor(dd$TL/TL.bins.cm)
      y=table(dd$bin)
      y=y/max(y)
      points(as.numeric(names(y)),y,type='h',col="forestgreen",lwd=1.5)
      d <- density(dd$TL)
      # lines(d$x,d$y/max(d$y),col="forestgreen",lwd=2)
    }
    #legend('bottomleft','size composition',lwd=2,col="forestgreen",bty='n')
    text(quantile(TL,.5),.5,Flts[x],cex=3,col='grey80')
  }
}

fun.check.mean.weight=function(TL,a,b,Mean.weight,SD.Mean.weight,Sel,NM)
{
  PP=data.frame(TL=TL,TWT=a*TL^b)
  with(Sel,wrapper.fn(x=TL,a.dn=P_1,b.dn=P_2,c.dn=P_3,d.dn=P_4,e.dn=P_5,f.dn=P_6,a.log=1e5,b.log=0,out.logis=FALSE))
  MeanTL=PP$TL[which.min(abs(PP$TWT-Mean.weight))]
  SDMeanTL.low=PP$TL[which.min(abs(PP$TWT-(Mean.weight-SD.Mean.weight)))]
  SDMeanTL.high=PP$TL[which.min(abs(PP$TWT-(Mean.weight+SD.Mean.weight)))]
  abline(v=MeanTL,lwd=2,col=2)
  abline(v=SDMeanTL.low,lty=3,lwd=1.5,col='grey30')
  abline(v=SDMeanTL.high,lty=3,lwd=1.5,col='grey30')
  text(PP$TL[which.min(abs(PP$TWT-Mean.weight))],0,"TL at mean meanbodywt input for SS",pos=4,srt=90,col=2)
  mtext(NM)
}