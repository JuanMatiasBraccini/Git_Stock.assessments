#####################Original functions (glm approach) #############
require(msm)
gillnetfit=function(data,meshsizes,type="norm.loc",rel=NULL,
                    plots=c(T,T),plotlens=NULL,plotlens_age=NULL,details=F)
{
  if(sum(sort(meshsizes)==meshsizes)!=length(meshsizes))
    stop("Mesh sizes must be ascending order")
  lens=rep(data[,1],ncol(data[,-1]))
  msizes=rep(meshsizes,rep(nrow(data),ncol(data[,-1])))
  msize1=msizes[1]
  dat=as.vector(data[,-1])
  var1=lens*msizes; var2=msizes^2; var3=(lens/msizes)^2 
  var4=lens/msizes; var5=-log(msizes); var6=log(msizes/msizes[1])
  var7=var6*log(lens) - 0.5*var6*var6; var8=lens*lens
  var9=msizes/lens
  if(is.null(plotlens)) plotlens=data[,1]
  if(is.null(rel)) os=0
  else os=rep(log(rel),rep(nrow(data),ncol(data[,-1])))
  switch(type,
         "norm.loc"={
           if(missing(rel))
             fit=glm(dat ~ -1 + var1 + var2 + as.factor(lens),family=poisson)
           else
             fit=glm(dat ~ -1 + var1 + var2 + as.factor(lens) + offset(os),family=poisson)
           x=coef(fit)[c("var1","var2")]
           varx=summary(fit)$cov.unscaled[1:2,1:2]
           k=-2*x[2]/x[1]; sigma=sqrt(-2*x[2]/(x[1]^2))
           vartemp=deltamethod(list(~-2*x2/x1,~sqrt(-2*x2/(x1^2))),x,varx,ses=F)
           pars=c(k,sigma,k*msizes[1],sigma)
           form1=as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
           varpars=deltamethod(list(~x1,~x2,form1,~x2),c(k,sigma),vartemp,ses=F)
           gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars))) 
           rownames(gear.pars)=c("k","sigma","mode(mesh1)","std_dev(all meshes)") },
         "norm.sca"={
           if(missing(rel))
             fit=glm(dat ~ -1 + var3 + var4 + as.factor(lens),family=poisson)
           else
             fit=glm(dat ~ -1 + var3 + var4 + as.factor(lens) + offset(os),family=poisson)
           x=coef(fit)[c("var3","var4")]
           varx=summary(fit)$cov.unscaled[1:2,1:2]
           k1=-x[2]/(2*x[1]); k2=-1/(2*x[1])
           vartemp=deltamethod(list(~-x2/(2*x1),~-1/(2*x1)),x,varx,ses=F)
           pars=c(k1,k2,k1*msizes[1],sqrt(k2*msizes[1]^2))
           form1=as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
           form2=as.formula(sprintf("~sqrt(x2*%f^2)",msize1)) #Deltamethod quirk
           varpars=deltamethod(list(~x1,~x2,form1,form2),c(k1,k2),vartemp,ses=F)
           gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars))) 
           rownames(gear.pars)=c("k1","k2","mode(mesh1)","std_dev(mesh1)") },
         "gamma"={
           if(missing(rel))
             fit=glm(dat ~ -1 + var4 + var5 + as.factor(lens),family=poisson)
           else
             fit=glm(dat ~ -1 + var4 + var5 + as.factor(lens) + offset(os),family=poisson)
           x=coef(fit)[c("var4","var5")] 
           varx=summary(fit)$cov.unscaled[1:2,1:2]
           alpha=x[2]+1; k=-1/x[1]
           vartemp=deltamethod(list(~x2+1,~-1/x1),x,varx,ses=F)
           pars=c(alpha,k,(alpha-1)*k*msizes[1],sqrt(alpha*(k*msizes[1])^2))
           form1=as.formula(sprintf("~(x1-1)*x2*%f",msize1)) #Deltamethod quirk
           form2=as.formula(sprintf("~sqrt(x1*(x2*%f)^2)",msize1)) #Deltamethod quirk
           varpars=deltamethod(list(~x1,~x2,form1,form2),c(alpha,k),vartemp,ses=F)
           gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))  
           rownames(gear.pars)=c("alpha","k","mode(mesh1)","std_dev(mesh1)")  },
         "lognorm"={
           if(missing(rel))
             fit=glm(dat ~ -1 + var6 + var7 + as.factor(lens),family=poisson)
           else
             fit=glm(dat ~ -1 + var6 + var7 + as.factor(lens) + offset(os),family=poisson)
           x=coef(fit)[c("var6","var7")] 
           varx=summary(fit)$cov.unscaled[1:2,1:2]
           mu1=-(x[1]-1)/x[2]; sigma=sqrt(1/x[2])
           vartemp=deltamethod(list(~-(x1-1)/x2,~sqrt(1/x2)),x,varx,ses=F)
           pars=c(mu1,sigma,exp(mu1-sigma^2),sqrt(exp(2*mu1+sigma^2)*(exp(sigma^2)-1)))
           varpars=deltamethod(list(~x1,~x2,~exp(x1-x2^2),
                                    ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),c(mu1,sigma),vartemp,ses=F)
           gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))  
           rownames(gear.pars)=c("mu1(mode log-scale, mesh1)","sigma(std_dev log scale)",
                                 "mode(mesh1)","std_dev(mesh1)")  },
         stop(paste("\n",type, "not recognised, possible curve types are ", 
                    "\"norm.loc\", \"norm.sca\", \"gamma\", and \"lognorm\"")))
  rselect=rcurves(type,meshsizes,rel,pars,plotlens)
  rselect_age=rcurves(type,meshsizes,rel,pars,plotlens_age)
  
  devres=matrix(resid(fit,type="deviance"),nrow(data),ncol(data[,-1]))
  if(plots[1]) plot.curves(type,plotlens,rselect)
  if(plots[2]) plot.resids(devres,meshsizes,data[,1])
  g.o.f=c(deviance(fit),sum(resid(fit,type="pearson")^2),fit$df.res,fit$null)
  names(g.o.f)=c("model_dev","Pearson chi-sq","dof","null_dev")
  fit.type=paste(paste(type,ifelse(is.null(rel),"",": with unequal mesh efficiencies")))
  
  if(details==F)
  {
    return(list(fit.type=fit.type,gear.pars=gear.pars,fit.stats=g.o.f))
  }else 
  {
    return(list(fit.type=fit.type,gear.pars=gear.pars,fit.stats=g.o.f,
                devres=devres,rselect=rselect,rselect_age=rselect_age,
                type=type,meshsizes=meshsizes,plotlens=plotlens,lens=data[,1]))
  }
}

#Calculate the relative selection curves
rcurves=function(type,meshsizes,rel,pars,plotlens) {
  lens=rep(plotlens,length(meshsizes))
  relsizes=meshsizes/meshsizes[1]
  if(is.null(rel)) releff=1
  else releff=rep(rel,rep(length(plotlens),length(meshsizes)))
  msizes=rep(meshsizes,rep(length(plotlens),length(meshsizes)))
  relsizes=rep(relsizes,rep(length(plotlens),length(relsizes)))
  switch(type,
         "norm.loc"={ k=pars[1]; mean1=pars[3]; sigma=pars[2]
         rselect=exp(-((lens-k*msizes)^2)/(2*sigma^2)) },
         "norm.sca"={ k1=pars[1]; k2=pars[2]
         rselect=exp(-((lens-k1*msizes)^2)/(2*k2*msizes^2)) },
         "gamma"={ alpha=pars[1]; k=pars[2]
         rselect=(lens/((alpha-1)*k*msizes))^(alpha-1)*exp(alpha-1-lens/(k*msizes)) },
         "lognorm"={ mu1=pars[1]; sigma=pars[2]
         rselect=
           (1/lens)*exp(mu1+log(relsizes)-sigma^2/2-
                          (log(lens)-mu1-log(relsizes))^2/(2*sigma^2)) })   
  rselect=releff*rselect/max(releff)
  rselect=matrix(rselect,ncol=length(meshsizes))
  return(rselect)
}

#Plot the relative selection curves
plot.curves=function(type,plotlens,rselect,cOL)
{
  plot.title=switch(type,
                    "norm.loc"="Normal (common spread)",
                    "norm.sca"="Normal",
                    "gamma"="Gamma",
                    "lognorm"="Log-normal")
  plot.title=paste(plot.title,"retention curve")
  matplot(plotlens,rselect,type="l",lty=1,las=1,ylim=c(0,1),lwd=2,col=cOL,
          xlab="Length (cm)",ylab="Relative retention",main=plot.title)
}

#Plot the deviance residuals matrix
plot.resids=function(residuals,msizes,lens,cex=1,title="Deviance residuals",...)
{
    if(missing(lens)) lens=1:nrow(residuals)
    if(missing(msizes)) msizes=1:ncol(residuals)
    plot(c(min(lens),max(lens)),range(msizes),xlab="",ylab="",
         ylim=range(msizes)+(cex/25)*c(-1,1)*(max(msizes)-min(msizes)),
         type="n",main=title,...)
    for(i in 1:nrow(residuals))
      for(j in 1:ncol(residuals))
        points(lens[i],msizes[j],pch=ifelse(residuals[i,j]>0,16,1),
               cex=3*abs(residuals[i,j])*cex/(abs(max(residuals))))
}


##################### Next generation functions (optim approach) ##################
NetFit=function(Data,Meshsize,x0,rtype="norm.loc",rel.power=NULL) {
  if(sum(sort(Meshsize)==Meshsize)!=length(Meshsize))
    stop("Mesh size must be ascending order")
  if(is.null(rel.power)) rel.power=rep(1,length(Meshsize))
  Counts=Data[,-1]
  if(ncol(Counts)!=length(Meshsize))
    stop("Number of mesh sizes should be ",ncol(Counts))
  CountPropns=Counts/apply(Counts,1,sum,na.rm=TRUE)
  fullfit.l=sum(Counts*log(CountPropns),na.rm=TRUE)
  r=selncurves(rtype) #Get selection curve function
  fit=optim(x0,nllhood,Data=Data,Meshsize=Meshsize,r=r,rel.power=rel.power,
            hessian=T,control=list(trace=F))
  cat("Parameters=",fit$par,",    Deviance=",2*(fullfit.l+fit$value),"\n")
  invisible(c(fit,deviance=deviance,rtype=rtype,rel.power=list(rel.power),
              Meshsize=list(Meshsize),Data=list(Data))) }

nllhood=function(theta,Data,Meshsize,r,rel.power) {
  lens=Data[,1]; Counts=Data[,-1]
  rmatrix=outer(lens,Meshsize,r,theta)
  rmatrix[is.na(Counts)]=NA #No fitted retention for missing meshsizes
  rmatrix=t(t(rmatrix)*rel.power)
  phi=rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
  nll=-sum(Counts*log(phi),na.rm=TRUE)
  return(nll) }

Estimates=function(fit) {
  x=fit$par; varx=solve(fit$hess) 
  names=c("Mode(mesh1)","Std dev.(mesh1)")
  switch(fit$rtype,
         "norm.loc"={ pars=x; varpars=varx },
         "norm.sca"={ pars=x; varpars=varx },
         "lognorm"={
           pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)))
           varpars=deltamethod(list(~exp(x1-x2^2),
                                    ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),x,varx,ses=F)},
         "binorm.sca"={
           pars=c(x[1:4],exp(x[5])/(1+exp(x[5])))
           names=c("Mode1(mesh1)","Std dev.1(mesh1)",
                   "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
           varpars=deltamethod(list(~x1,~x2,~x3,~x4,~exp(x5)/(1+exp(x5))),
                               x,varx,ses=F)},
         "bilognorm"={
           pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)),
                  exp(x[3]-x[4]^2),sqrt(exp(2*x[3]+x[4]^2)*(exp(x[4]^2)-1)),
                  exp(x[5])/(1+exp(x[5])))
           names=c("Mode1(mesh1)","Std dev.1(mesh1)",
                   "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
           varpars=deltamethod(
             list(~exp(x1-x2^2),~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1)),
                  ~exp(x3-x4^2),~sqrt(exp(2*x3+x4^2)*(exp(x4^2)-1)),
                  ~exp(x5)/(1+exp(x5))),x,varx,ses=F)},  
         "tt.logistic"={
           pars=c(-x[1]/x[2],2*(log(3))/x[2],exp(x[3])/(1+exp(x[3])))
           names=c("L50","SR","p")
           varpars=deltamethod(list(~-x1/x2,~2*log(3)/x2,~exp(x3)/(1+exp(x3))),
                               x,varx,ses=F)},
         stop(paste("\n",fit$rtype, "not recognised, possible curve types are \n", 
                    "\"norm.loc\", \"norm.sca\", \"lognorm\" \n", 
                    "\"binorm.sca\", \"bilognorm\", and \"tt.logistic\"")) 
  )#End of switch
  estimates=cbind(pars,sqrt(diag(varpars)))
  colnames(estimates)=c("par","s.e.")
  rownames(estimates)=names
  return(estimates) }

PlotCurves=function(fit,Meshsize=NULL,plotlens=NULL,standardize=TRUE,...) {
  r=selncurves(fit$rtype) #Get selection curve function
  if(is.null(plotlens)) plotlens=fit$Data[,1]
  if(is.null(Meshsize)) Meshsize=fit$Meshsize
  plot.title=switch(fit$rtype,
                    "norm.loc"="Normal (common spread)",
                    "norm.sca"="Normal",
                    "lognorm"="Lognormal",
                    "binorm.sca"="Bi-normal",
                    "bilognorm"="Bi-lognormal",
                    "tt.logistic"="Control and logistic","")
  rmatrix=outer(plotlens,Meshsize,r,fit$par)
  rmatrix=t(t(rmatrix)*fit$rel.power)
  if(standardize) rmatrix=rmatrix/max(rmatrix)
  matplot(plotlens,rmatrix,type="l",las=1,ylim=c(0,1),
          xlab="Length (cm)",ylab="Relative retention",...) 
  #abline(h=seq(0,1,0.25),lty=3)
  lenrmatrix=cbind(plotlens,rmatrix)
  colnames(lenrmatrix)=c("Length",Meshsize)
  invisible(lenrmatrix) }

Summary=function(fit,label="Deviance residuals",
                 xlabel="Length (cm)",ylabel="Mesh size (cm)",cex=1) {
  r=selncurves(fit$rtype) #Get selection curve function
  lens=fit$Data[,1]; nlens=length(lens)
  Meshsize=fit$Meshsize; nmeshes=length(Meshsize)
  O=fit$Data[,-1]; #Matrix of observed counts
  rmatrix=outer(lens,Meshsize,r,fit$par)
  rmatrix[is.na(O)]=NA #No fitted retention for missing meshsizes
  rmatrix=t(t(rmatrix)*fit$rel.power)
  phi=rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
  E=apply(O,1,sum,na.rm=TRUE)*phi #Matrix of expected counts
  Pearson.resids=(O-E)/sqrt(E)
  Pearson.chisq=sum(Pearson.resids^2,na.rm=TRUE)
  wk=O*log(O/E); wk[is.na(wk)]=0
  Dev.resids=sign(O-E)*sqrt(2*(E-O+wk))
  Deviance=sum(Dev.resids^2,na.rm=TRUE)
  full.l=sum(-O+O*log(O),na.rm=TRUE)
  null.E=matrix(apply(O,1,mean,na.rm=TRUE),nrow=nlens,ncol=nmeshes)
  null.l=sum(-null.E+O*log(null.E),na.rm=TRUE)
  model.l=sum(-E+O*log(E),na.rm=TRUE)
  NonZeroDat=O[apply(O,1,sum,na.rm=TRUE)>0,]
  d.o.f.=nrow(NonZeroDat)*(nmeshes-1)-length(fit$par)-sum(is.na(NonZeroDat))
  out=rbind(null.l,model.l,full.l,Deviance,Pearson.chisq,d.o.f.)
  AreLensUnique=(length(lens)==length(unique(lens)))
  if(nmeshes>2&AreLensUnique) {
    plot(1,1,xlim=range(lens),xlab=xlabel,ylab=ylabel,
         ylim=range(Meshsize)+(cex/50)*c(-1,1)*(max(Meshsize)-min(Meshsize)),
         yaxt="n",type="n",main=label)
    axis(2,Meshsize,Meshsize,las=1)   
    for(i in 1:nlens)
      for(j in 1:nmeshes)
        points(lens[i],Meshsize[j],pch=ifelse(Dev.resids[i,j]>0,16,1),
               cex=3*abs(Dev.resids[i,j])*cex/(abs(max(Dev.resids)))) }
  else
    if(nmeshes==2) {
      Dev.resids.len=sign(Dev.resids[,2])*sqrt(apply(Dev.resids^2,1,sum))
      plot(lens,Dev.resids.len,type=ifelse(AreLensUnique,"h","p"),las=1,
           main=label,xlab=xlabel,ylab=ylabel,cex=cex)
      abline(h=0) }
  return(out)
}

