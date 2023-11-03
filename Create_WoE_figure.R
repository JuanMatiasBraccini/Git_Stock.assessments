if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Population dynamics/Consequence_likelihood_plot.R'))



setwd(handl_OneDrive('Analyses/Population dynamics/WoE'))

#Figure 4. WOE risk assessment

  #Fill in cons-like matrix

#Outcomes from 2017 WoE assessment
WoE.year=2017
Species=c("Whiskery shark","Gummy shark","Dusky shark","Sandbar shark")
List=vector('list',length(Species))
names(List)=Species
List$"Whiskery shark"=data.frame(Consequence=1:4,Likelihood=c(4,2,2,1))
List$"Gummy shark"=data.frame(Consequence=1:4,Likelihood=c(4,2,2,0))
List$"Dusky shark"=data.frame(Consequence=1:4,Likelihood=c(2,3,2,1))
List$"Sandbar shark"=data.frame(Consequence=1:4,Likelihood=c(3,3,2,1))

  #Plot matrix
fn.fig(paste(WoE.year,".Risk_based_WoE",sep=''),2400,2400) 
smart.par(n.plots=length(List),MAR=c(1,1,.1,.1),OMA=c(3.5,6,.1,.1),MGP=c(2.5,.7,0))
for(i in 1:length(List))
{
  fun.cons.like.mat(TAB=List[[i]],CX=10,Species=Species[i])
  
  if(i %in%c(1,3))axis(2,1:4,c(expression(B[Cur] > B[Tar]),
                               expression(paste("B"["Thr"],"<","B"["Cur"],"<","B"["Tar"])),
                               expression(paste("B"["Lim"],"<","B"["Cur"],"<","B"["Thr"])),
                               expression(B[Cur] < B[Lim])),cex.axis=1.18)
  if(i %in%3:4)
  {
    axis(1,1:4,c("Remote","Unlikely","Possible","Likely"),cex.axis=1.2)
    axis(1,1:4,c("(<5%)","(5-20%)","(20-50%)","(>50%)"),cex.axis=1.18,line=1,col="transparent")
    
  }
  
  #axis(1,1:4,c("<5","5 - 20","20 - 50",expression(phantom(x) >=50)))
}
mtext("Consequence",2,4.5,outer=T,las=3,cex=1.5)
mtext("Likelihood",1,2,outer=T,cex=1.5)
dev.off()
