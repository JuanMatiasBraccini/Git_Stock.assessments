#library(ggplot2)
library(reshape)
library(tidyverse)
library(plyr)
#library(diags)
#library(dplyr)
library(gam)
#library(GGally)
library(corrplot)
#library(devtools)
#devtools::install_github("flr/FLCore", force=T)
library(FLCore)
#devtools::install_github("flr/diags", force=T)
#library(diags)
library(weights)
library(gridExtra)

theme_set(theme_bw(14))
options(digits=3)

#dirMy="C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/Population dynamics/Other people's code/Joel Rice"


#Dat=read.csv(paste(dirMy,"S_ALT_BSH.csv",sep='/')) 


fn.cpue.corr=function(Dat,AREA,WD,R=NULL)
{
  flts = names(Dat)[-1]
  u=Dat[-1,]
  names(u)[1]="year"
  
  unit=c(rep("CPUE",length(flts)))
  area=c(rep(AREA,length(flts)))   # this name has to be repeated in the plotting calls i.e. South_Atlantic_Ocean
  
  u= subset(melt(u,id="year"),!is.na(value))
  u=transform(u,area=area[variable],
              name=flts[variable],
              unit=unit[variable])[,-2]
  #u=ddply(u,.(name,area), transform, value2=diags:::stdz(value)) # DC 4/13/2017
  u=ddply(u,.(name,area), transform, value2=weights:::stdz(value)) # DC 4/13/2017
  
  scale<-function(x,y,...){ #DC Notes 4/13/2017 begin Function Call
    args=list(...)
    
    if (length(args)==0) group=rep(1,length(x)) else group=args[[1]]  
    
    gm=gam(y~lo(x)+group,data=data.frame(x=x,y=y,group=group))
    
    res=data.frame(hat =predict(gm),
                   y     =gm$y,
                   x     =x,
                   group =group,
                   scl   =c(0,coefficients(gm)[-(1:2)])[as.numeric(as.factor(group))]
    )
    res$y  =res$y  -res$scl
    res$hat=res$hat-res$scl
    
    if (length(args)==1) names(res)[4]=names(args)[1]  
    
    res[,-5]
  } #DC Notes 4/13/2017 end Function Call
  
  my_density <- function(data,mapping,...){
    ggplot(data=data,mapping=mapping)+
      geom_density(...,lwd=1)}
  
  my_bar <- function(data,mapping,...){
    ggplot(data=data,mapping=mapping)+
      geom_bar(...)}
  
  my_smooth <- function(data,mapping,...){
    ggplot(data=data,mapping=mapping)+
      geom_point(...,size=.5)+
      geom_smooth(...,method="lm",se=FALSE)}
  
  my_ccf<-function(cpue){
    cc=mdply(expand.grid(a=names(cpue),b=names(cpue)),
             function(a,b){
               #print(paste(a,b))
               res=model.frame(mcf(FLQuants(cpue[c(a,b)])))
               res=subset(res,!is.na(res[,7])&!is.na(res[,8]))
               
               if (dim(res)[1]>10){
                 res=data.frame(lag=-10:10,data=ccf(res[,7],res[,8],plot=F,
                                                    lag.max=10)$acf)
                 return(res)}else{return(NULL)}}
    )}
  
  dat=ddply(u,.(area), function(x) scale(x$year,x$value2,name=as.character(x$name)))
  
  #1. Time series of residuals from the smooth fit to CPUE indices. 
  p1=ggplot(merge(dat,data.frame(name=flts,unit=unit)))+
    geom_line( aes(x,hat),data=dat[,c("x","hat","area")],col="grey60")+
    geom_line( aes(x,hat))+
    geom_line( aes(x,y,  col=unit))+
    geom_point(aes(x,y,  col=unit))+
    facet_grid(name~area,scale="free",space="free_x")+
    theme_bw(14)+
    theme(legend.position="bottom")+
    labs(title='residuals from the smooth fit to CPUE indices',
         caption='X-axis is time, Y-axis are the scaled indices')
  
  #2. Pairwise scatter plots for CPUE indices.
  mat=cast(subset(u,area==AREA),year~name,value="value2")
  names(mat)=gsub(" ", "_",names(mat))
  mat=mat[,-1]%>%data.frame()
  mat=mat[!rowSums(is.na(mat))>1,]
  
  p2=NULL
  if(ncol(mat)>2)
  {
    p2=ggpairs(mat,columns=1:2,
               upper=list(continuous=wrap("cor", na.rm = TRUE,size=4, hjust=0.5)),
               lower=list(continuous = wrap(my_smooth, na.rm = TRUE)),
               diag=list(continuous="bar"))+          
      theme(legend.position="bottom")
  }

  
  
  #3. Correlation matrix for CPUE indices
  cr=cor(mat,use="pairwise.complete.obs")
  dimnames(cr)=list(gsub("_"," ",names(mat)),gsub("_"," ",names(mat)))
  cr[is.na(cr)]=0
  p3=NULL         

  
  #4. Cross-correlations between CPUE indices
  p4=NULL
  do.p4='no'
  if(do.p4!='no')
  {
    cpue=FLQuants(dlply(subset(u,area==AREA),.(name), with, 
                        as.FLQuant(data.frame(year=year,data=value2))))
    cc=my_ccf(cpue)
    p4=ggplot(cc)+
      geom_linerange(aes(x=lag,ymin=0,ymax=data))+
      facet_grid(a~b)+
      geom_vline(aes(xintercept=0))+
      theme_bw(14)+
      labs(title='Cross-correlations between CPUE indices',
           caption='identify lagged correlations (e.g., due to year-class effects). 
                   X-axis is lag number, and y-axis is cross-correlation')
    
  }
   
  if(!is.null(R))
  {
    dd=vector('list',length(flts))
    for(pp in 1:length(dd))
    {
      dd[[pp]]=Dat[,match(c('Year',flts[pp]),names(Dat))]%>%
              mutate(Lagged=lag(get(flts[pp])),
                     Delta=get(flts[pp])-Lagged,
                     Increment=ifelse(Delta>0,(get(flts[pp])*100/Lagged)-100,NA),
                     Mean=flts[pp])%>%
        dplyr::select(Year,Mean,Increment)
    }
    dd=do.call(rbind,dd)
    
    p5=Dat%>%gather(Mean,Fleet,-Year)%>%full_join(dd,by=c("Year","Mean"))%>%
      mutate(LBL=ifelse(!is.na(Increment) & Increment>2*R,paste(round(Increment),'%'),''))%>%
      ggplot(aes(Year,Fleet))+
      geom_point(col=2,size=3,alpha=.5)+geom_line(col=2,alpha=.5)+
      geom_text_repel(aes(label=LBL))+
      facet_wrap(~Mean)+ylab("Relative cpue")+
      labs(title='Percentage increase between years',
           caption=paste('Population growth rate (r) =', 100*round(R,2),'% per annum'))
    
    D.p6=Dat%>%gather(Mean,Fleet,-Year)%>%full_join(dd,by=c("Year","Mean"))
    my_formula=y ~ x
    p6=D.p6%>%
      ggplot(aes(Year,Fleet)) +
      geom_point(shape = 21, size = 5)+facet_wrap(~Mean)+
      geom_smooth(method = "lm", data = D.p6%>%filter(!is.na(Increment)),
                  se = F, fullrange = TRUE,colour="red")+
      stat_poly_eq(data = D.p6%>%filter(!is.na(Increment)),
                   aes(label =  paste(stat(eq.label),stat(adj.rr.label),stat(p.value.label),
                                      sep = "*\", \"*")),
                   formula = my_formula, parse = TRUE,
                   label.y = "bottom", label.x = "right", size = 4)
  }
  
  #Export figures
  p=list(p1,p2,p3,p4)
  if(exists('p5')) p=list(p1,p2,p3,p4,p5,p6)
  pdf(paste(WD, "CPUE correlations.pdf",sep='/'))
  for(g in 1:length(p))
  {
    if(g!=3) print(p[[g]])
    if(g==3)  corrplot(cr,diag=F,order="hclust",addrect=2,main='CPUE correlation matrix')
  }
    
  dev.off()


  
}




