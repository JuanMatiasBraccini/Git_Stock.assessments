#note: derive F from catch curve and gear selectivity
#      assume start of year in January (coincides with birth of most species)  
#      size is TL in mm
mm.conv=10 # convert total length to mm  
fn.source("Length_based.catch.curve.R")


#24.1 TDGDLF 
size.catch.curve_TDGDLF=vector('list',N.sp)
names(size.catch.curve_TDGDLF)=Keep.species
system.time({for(l in 1: N.sp)  
{
  # Calculate F
  Outfile='TDGDLF' 
  this.size.comp=paste('Size_composition',c('West.6.5','West.7','Zone1.6.5','Zone1.7','Zone2.6.5','Zone2.7'),sep="_")
  outfile=paste(Outfile,'_histogram',sep='')
  
  #get size composition
  iid=Species.data[[l]][fn.extract.dat(this.size.comp,names(Species.data[[l]]))]
  if(length(iid)>0)
  {
    print(paste("Size-based catch curve with dome-shape selectivity for --",names(Species.data)[l],"---",Outfile))
    
    dummy=do.call(rbind,iid)
    if(grepl("TDGDLF",outfile))
    {
      dummy=dummy%>%
        mutate(dummy=sub(".*Size_composition_", "", rownames(dummy)),
               Mesh=word(dummy,2,sep = "\\."),
               Mesh=ifelse(Mesh=='6','6.5',Mesh),
               Mesh=factor(Mesh,levels=c('6.5','7')),
               Zone=word(dummy,1,sep = "\\."))
      
    }
    dummy=dummy%>%filter(year<=as.numeric(substr(Last.yr.ktch,1,4)))
    
    N.min=dummy%>%
      group_by(year)%>%
      tally()%>%
      filter(n>=Min.annual.obs.ktch)%>%
      mutate(Keep=year)
    if(nrow(N.min)>0)
    {
      dummy=dummy%>%
        mutate(Keep=year)%>%
        filter(Keep%in%N.min$Keep)%>%
        mutate(TL=mm.conv*FL*List.sp[[l]]$a_FL.to.TL+List.sp[[l]]$b_FL.to.TL)     
      
      #1. Plot observed size frequency by year and mesh for years with minimum sample size
      if(grepl("TDGDLF",outfile))
      {
        p=dummy%>%
          ggplot( aes(x=TL/mm.conv, color=Mesh, fill=Mesh)) +
          geom_histogram(alpha=0.6, binwidth = TL.bins.cm)
        WHERE="top"
      }else
      {
        p=dummy%>%
          ggplot( aes(x=TL/mm.conv,color=year, fill=year)) +
          geom_histogram(alpha=0.6, binwidth = TL.bins.cm)
        WHERE="none"
      }
      p=p+
        facet_wrap(~year,scales='free_y')+
        xlab("Total length (cm)")+ylab("Count")+
        theme(legend.position=WHERE,
              legend.title = element_blank(),
              legend.text=element_text(size=14),
              strip.text.x = element_text(size = 12),
              axis.text=element_text(size=12),
              axis.title=element_text(size=16))
      print(p)
      ggsave(paste(handl_OneDrive("Analyses/Population dynamics/1."),
                   capitalize(List.sp[[l]]$Name),"/",AssessYr,
                   "/1_Inputs/Visualise data/Size.comp.",outfile,".tiff",sep=''),
             width = 12,height = 10,compression = "lzw")
      
      #2. Calculate F
      #initial value of F   
      params = log(0.2)   
      
      # age
      MaxAge = ceiling(mean(List.sp[[l]]$Max.age.F))
      
      # ages
      Ages = 0:MaxAge
      
      # Growth parameters
      Linf = mm.conv*with(List.sp[[l]],Growth.F$FL_inf*a_FL.to.TL+b_FL.to.TL)  #total length in mm 
      vbK = List.sp[[l]]$Growth.F$k
      Lo =  mm.conv*with(List.sp[[l]],Lzero*a_FL.to.TL+b_FL.to.TL)   #total length in mm   
      #tzero = 0              
      CVLenAtAge = 0.1      
      
      
      # length structure of model
      MaxLen = mm.conv*10*round(List.sp[[l]]$TLmax/10)
      LenInc = mm.conv*TL.bins.cm
      
      # length bins
      min.Lo=min(Lo-Lo*CVLenAtAge,min(dummy$TL))
      lbnd = seq(LenInc*floor((min.Lo)/LenInc),MaxLen - LenInc, LenInc)
      ubnd = lbnd + LenInc
      midpt = lbnd + (LenInc/2)
      nLenCl = length(midpt)
      
      # natural mortality  
      NatMort= mean(colMeans(store.species.M[[l]],na.rm=T))
      
      # gillnet selectivity (6.5 and 7 inch combined)
      SelAtLength=Selectivity.at.totalength[[l]]%>%               
        mutate(TL=TL*mm.conv)%>%
        filter( TL%in%midpt)%>%
        pull(Sel.combined)
      #plot(midpt,SelAtLength)
      
      #Execute model functions           
      MeanSizeAtAge = CalcMeanSizeAtAge(Lo,Linf, vbK)
      RecLenDist = CalcSizeDistOfRecruits(GrowthCurveResults=MeanSizeAtAge, CVLenAtAge)
      #plot(midpt,RecLenDist)
      
      LTM = CalcLTM(Linf, vbK, CVLenAtAge, midpt)   
      #image(1:nrow(LTM),1:ncol(LTM),as.matrix(LTM))
      
      #Fit model for each year with length data using parallel processing
      nyrs=dummy%>%
        group_by(year)%>%
        tally()%>%pull(year)
      add.dummy=data.frame(bin=midpt)
      
      cl <- makeCluster(detectCores()-1)
      registerDoParallel(cl)
      clusterEvalQ(cl, .libPaths('C:/~/R/win-library/4.0'))  #added bit to point where doparallel is 
      F.at.year=foreach(q=1:length(nyrs),.packages=c('dplyr','doParallel'),.export=c('ObjectiveFunc')) %dopar%
        {
          x=dummy%>%
            filter(year==nyrs[q])%>%
            mutate(bin=LenInc*floor(TL/LenInc)+LenInc/2)%>%
            group_by(bin)%>%
            tally()%>%
            full_join(add.dummy,by='bin')%>%
            arrange(bin)%>%
            mutate(n=ifelse(is.na(n),0,n))%>%
            filter(bin%in%midpt)
          ObsCatchFreqAtLen=x%>%pull(n)
          nlmb <- nlminb(params, ObjectiveFunc, gradient = NULL, hessian = TRUE)
          
          
          #Calculate uncertainty for parameter estimates
          #note: get variance-covariance matrix from fitted model
          hess.out = optimHess(nlmb$par, ObjectiveFunc)
          vcov.Params = solve(hess.out)
          ses = sqrt(diag(vcov.Params)) # asymptotic standard errors of parameter estimates
          EstFMort = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
          Table.check=data.frame(year=nyrs[q],
                                 objective.fun=nlmb$objective,
                                 convergence=nlmb$convergence,
                                 par.low95=EstFMort[2],
                                 par=EstFMort[1],
                                 par.up95=EstFMort[3])   
          
          #estimated
          ExpCatchAtLen = GetExpCatchAtLen(nlmb$par)
          
          return(list(Table=Table.check,
                      Obs=ObsCatchFreqAtLen/sum(ObsCatchFreqAtLen),
                      Pred=unlist(ExpCatchAtLen),
                      Size=midpt))
          rm(ExpCatchAtLen,ObsCatchFreqAtLen)
        }
      stopCluster(cl)
      names(F.at.year)=nyrs
      size.catch.curve_TDGDLF[[l]]=F.at.year
      rm(dummy,SelAtLength,NatMort,midpt,Ages,MeanSizeAtAge,RecLenDist,LTM,nyrs)
      
      
      #Export F  
      PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
                 capitalize(List.sp[[l]]$Name),"/",AssessYr,"/Size_based.Catch.curve",sep='')
      if(!file.exists(file.path(PATH))) dir.create(file.path(PATH))   
      out=do.call(rbind,subListExtract(size.catch.curve_TDGDLF[[l]],"Table"))%>%
        rename(Low95=par.low95,
               Mean=par,
               Up95=par.up95)
      write.csv(out,paste(PATH,paste(unlist(strsplit(Outfile, split='_', fixed=TRUE))[1],"F.csv",sep="."),sep='/'),row.names = F) 
      
    }  #end of "if  nrow(N.min)>0 "statement
    
  }  #end of "if length(iid)>0" statement
}  # end l  
})    #takes 10 mins

#obs vs pred  
for(l in 1: N.sp)
{
  dummy=size.catch.curve_TDGDLF[[l]]
  if(!is.null(dummy))
  {
    PATH=paste(handl_OneDrive("Analyses/Population dynamics/1."),
               capitalize(List.sp[[l]]$Name),"/",AssessYr,"/Size_based.Catch.curve",sep='')
    
    fn.fig(paste(PATH,"/TDGDLF_F.fit",sep=""),2400,2400)
    smart.par(n.plots=length(dummy),MAR=c(1.2,2,2.5,1.25),OMA=c(2,2,.2,2.1),MGP=c(1,.5,0))
    for(i in 1:length(dummy))
    {
      with(dummy[[i]],{
        plot(Size/mm.conv, Obs, xlab="", ylab="",pch=19, col="orange",cex=1.25)
        lines(Size/mm.conv, Pred,lwd=1.25)
        mtext(Table$year, cex=1.25,3)
      })
    }
    legend("topright",c("Observed","Predicted"),pch=c(19,NA),lty=c(NA,1),
           col=c("orange","black"),bty='n',cex=1.25)
    mtext("Total length (cm)",1,outer=T,cex=1.25,line=.5)
    mtext("Probability",2,outer=T,las=3,cex=1.25,line=.5)
    dev.off()
    
  }
  rm(dummy)
}

#27.2 NSF
#note: Not possible. Unknown selectivity for longline & only enough data for sandbar