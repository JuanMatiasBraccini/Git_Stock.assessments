fn.compare.empirical.vs.SS.selectivities=function(Empiric,x,dn1,dn2,dn.init)
{
   DN_1=doubleNorm24.fn(x,a=dn1[1],b=dn1[2], c=dn1[3], d=dn1[4], e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE)
   DN_2=doubleNorm24.fn(x,a=dn2[1],b=dn2[2], c=dn2[3], d=dn2[4], e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE)
   DN_init=doubleNorm24.fn(x,a=dn.init[1],b=dn.init[2], c=dn.init[3], d=dn.init[4], e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE)
   
   par(mfcol=c(2,1),oma=c(.2,2,1,1),mai=c(.1,.1,.5,.1))
   plot(x,DN_1,ylab='Selectivity',type='l',lwd=3,xlab='TL',main='female',col=4)
   lines(x,DN_2,col=2,lwd=3)
   lines(x,DN_init,col=3,lwd=3)
   KL=1:(ncol(Empiric)-1)
   for(pp in 2:ncol(Empiric))
   {
      lines(Empiric$TL.mm,Empiric[,pp],col=KL[pp-1],lty=2)
   }
   DD=Size.compo.SS.format%>%filter(Sex==1)
   DD.f=DD[,c("year",colnames(Size.compo.SS.format)[grep('f',colnames(Size.compo.SS.format))])]%>%
      gather(Size,N,-year)%>%
      mutate(Size=as.numeric(gsub("\\D", "", Size)))%>%
      group_by(year) %>%
      mutate(n = max(N)) %>%
      mutate(freq = N / n)
   yrs=sort(unique(DD.f$year))
   GC=gray.colors(length(yrs))
   for(y in 1:length(yrs))
   {
      aa=DD.f%>%filter(year==yrs[y])
      lines(aa$Size,aa$freq,col=GC[y],lwd=2)
   }
   
   plot(x,DN_1,ylab='Selectivity',type='l',lwd=3,xlab='TL',main='male',col=4)
   lines(x,DN_2,col=2,lwd=3)
   lines(x,DN_init,col=3,lwd=3)
   KL=1:(ncol(Empiric)-1)
   for(pp in 2:ncol(Empiric))
   {
      lines(Empiric$TL.mm,Empiric[,pp],col=KL[pp-1],lty=2)
   }
   legend('topleft',c('SS_Southern.1','SS_Southern.2','Init',names(Empiric)[-1],'Observed'),
          bty='n',col=c(4,2,3,KL,1),lwd=c(3,3,3,rep(1,length(KL)),2),lty=c(1,1,1,rep(2,length(KL)),1))
   
   DD=Size.compo.SS.format%>%filter(Sex==2)
   DD.m=DD[,c("year",colnames(Size.compo.SS.format)[grep('m',colnames(Size.compo.SS.format))])]%>%
      dplyr::select(-Nsamp)%>%
      gather(Size,N,-year)%>%
      mutate(Size=as.numeric(gsub("\\D", "", Size)))%>%
      group_by(year) %>%
      mutate(n = max(N)) %>%
      mutate(freq = N / n)
   yrs=sort(unique(DD.m$year))
   GC=gray.colors(length(yrs))
   for(y in 1:length(yrs))
   {
      aa=DD.m%>%filter(year==yrs[y])
      lines(aa$Size,aa$freq,col=GC[y],lwd=2)
   }
   
   
   #Changes thru time   
   yrs=sort(unique(DD.f$year))
   smart.par(length(yrs),MAR=c(2,3,1,1),OMA=c(2.5,1,.05,2.5),MGP=c(1.8,.5,0))
   for(y in 1:length(yrs))
   {
      aa=DD.f%>%filter(year==yrs[y])
      plot(x,DN_1,ylab='Selectivity',type='l',lwd=2,xlab='TL',main=paste('female',yrs[y]),col=4)
      lines(x,DN_2,col=2,lwd=2)
      lines(x,DN_init,col=3,lwd=2)
      lines(aa$Size,aa$freq,col=1,lwd=3)
   }
   legend('topleft',c('SS_Southern.1','SS_Southern.2','Init'),
          bty='n',col=c(4,2,3),lwd=c(3,3,3),lty=c(1,1,1))
   
   yrs=sort(unique(DD.m$year))
   smart.par(length(yrs),MAR=c(2,3,1,1),OMA=c(2.5,1,.05,2.5),MGP=c(1.8,.5,0))
   for(y in 1:length(yrs))
   {
      aa=DD.m%>%filter(year==yrs[y])
      plot(x,DN_1,ylab='Selectivity',type='l',lwd=2,xlab='TL',main=paste('male',yrs[y]),col=4)
      lines(x,DN_2,col=2,lwd=2)
      lines(x,DN_init,col=3,lwd=2)
      lines(aa$Size,aa$freq,col=1,lwd=3)
   }
   legend('topleft',c('SS_Southern.1','SS_Southern.2','Init'),
          bty='n',col=c(4,2,3),lwd=c(3,3,3),lty=c(1,1,1))
}
pdf(paste(this.wd,"Compare empirical selectivity vs SS3 selectivity.pdf",sep='/'))
fn.compare.empirical.vs.SS.selectivities(Empiric=Species.data$`whiskery shark`$gillnet.selectivity,
                                         x=60:160,
                                         dn1=c(a=129.44,b=-4, c=5.286, d=3.62),
                                         dn2=c(a=128.389,b=-4, c=5.33, d=3.04),
                                         dn.init=c(a=129.58,b=-4, c=5.5, d=4.51))
dev.off()

x=60:160
DN_1=doubleNorm24.fn(x,a=130.18,b=-7.91, c=5.31, d=3.63, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE)
plot(x,DN_1)
lines(x,doubleNorm24.fn(x,a=130.18,b=-4, c=5.31, d=3.63, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE))
lines(x,doubleNorm24.fn(x,a=130.18,b=-4, c=5.5, d=4.5, e=-999, f=-999,use_e_999=FALSE, use_f_999=FALSE),col=3)
