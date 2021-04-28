 #gummy                #Walker 2010
theta1=186
theta2=36695
mesh= 6.5 #(in inches, = 16.5 cm)
alphabeta.g=theta1*mesh
beta.g=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
alpha.g=alphabeta.g/beta.g
 
 #whiskery                #Simpfendorfer & Unsworth 1998
alpha.w=64.01339
beta.w=18.53164
alphabeta.w=alpha.w*beta.w

mid.FL.fem=60:180
sel.g=((mid.FL.fem*10/alphabeta.g)^alpha.g)*(exp(alpha.g-(mid.FL.fem*10/beta.g)))
sel.w=((mid.FL.fem*10/alphabeta.w)^alpha.w)*(exp(alpha.w-(mid.FL.fem*10/beta.w)))

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

setwd(handl_OneDrive('Analyses/Population dynamics/Visualise data'))
tiff(file=paste("Selectivty.length.Whi.Gum.tiff"),width = 2400, height = 2400,
     units = "px", res = 300,compression = "lzw")
par(mai=c(1,1,.1,.1),oma=c(.1,.1,.1,.1),las=1,xpd=T,mgp=c(3,.65,0))
plot(mid.FL.fem,sel.w,type='l',ylab='Relative selectivity',xlab="Fork length (cm)",las=1,
     cex.lab=2,cex.axis=1.5,col=2,lwd=4)
lines(mid.FL.fem,sel.g,col=3,lwd=4)
legend('topright',c('Whiskery shark','Gummy shark'),bty='n',lty=1,lwd=4,col=c(2,3),cex=1.5)
dev.off()