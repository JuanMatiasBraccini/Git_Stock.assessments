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

wrapper.fn=function(x,a.dn,b.dn,c.dn,d.dn,e.dn,f.dn,a.log,b.log,out.logis=TRUE)
{
  Sel.dn=doubleNorm24.fn(x,a=a.dn,b=b.dn, c=c.dn, d=d.dn, e=e.dn, f=f.dn,use_e_999=FALSE, use_f_999=FALSE)
  Sel.log=logistic1.fn(len=x,a=a.log, b=b.log)
  
  par(xpd=T)
  plot(x,Sel.dn,ylab='Selectivity',xlab='size',type='l')
  if(out.logis)
  {
    lines(x,Sel.log,col=2)
    legend('topleft',c(paste('double normal (a=',a.dn,', b=',b.dn,', c=',c.dn,', d=',d.dn,', e=',e.dn,', f=',f.dn,')',sep=''),
                       paste('logistic(a=',a.log,', b=',b.log,')',sep='')), inset=c(0,-0.1),
           lty=1,col=1:2,bty='n')
  }
  if(!out.logis)
  {
    legend('topleft',paste('double normal (a=',a.dn,', b=',b.dn,', c=',c.dn,', d=',d.dn,', e=',e.dn,', f=',f.dn,')',sep=''),
                       inset=c(0,-0.1),
           lty=1,col=1:2,bty='n')
  }
  
}
#wrapper.fn(x=85:300,a.dn=130,b.dn=-11,c.dn=9,d.dn=8,e.dn=-999,f.dn=-999,a.log=325,b.log=20)
#wrapper.fn(x=90:500,a.dn=150,b.dn=-1,c.dn=7,d.dn=8,e.dn=-999,f.dn=-999,a.log=325,b.log=10)
wrapper.fn(x=80:400,a.dn=110,b.dn=-11,c.dn=5,d.dn=10,e.dn=0.000001,f.dn=1,a.log=325,b.log=20,out.logis=FALSE)
