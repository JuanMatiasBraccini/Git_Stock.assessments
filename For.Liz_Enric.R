library(popbio)

#natural mortality (Jensen, Pauly,Hoenig, Theo)
M.fun=function(Amax,age.mat,Linf,k,Aver.T)
{
  m.Jensen.2=1.65/age.mat
  m.Jensen.2=rep(m.Jensen.2,length(age))
  
  #Pauly (1980)  
  m.Pauly=10^(-0.0066-0.279*log10(Linf)+0.6543*log10(k)+0.4634*log10(Aver.T))
  m.Pauly=rep(m.Pauly,length(age))
  
  #Hoenig (1983), combined teleost and cetaceans    
  m.Hoenig=exp(1.44-0.982*log(Amax))      
  m.Hoenig=rep(m.Hoenig,length(age))
  
  #Then et al (2015)
  m.Then.1=4.899*Amax^(-0.916)
  m.Then.1=rep(m.Then.1,length(age))
  
  
  #STEP 2. get mean at age
  nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Hoenig,m.Then.1)  

  
  return(rowMeans(nat.mort))
  apply(nat.mort, 1, function(x) weighted.mean(x, c(1,1,1.5,1.5)))
}
  
#r
Leslie=function(M,age.mat,Meanfec,CyclE)
{  
  #survivorship
  S=exp(-M)         
  
  #proportion surviving
  lx=rep(NA,length(age))
  lx[1]=1.0
  for (i in 2:(length(age)))lx[i]=lx[i-1]*S[i]
  
  #reproductive schedules   
  MF=c(rep(0,(age.mat-1)),Meanfec[age.mat:length(Meanfec)])
  mx=MF*sexratio/CyclE
  
  #probability of surviving (for birth-pulse, post-breeding census)
  px=vector(length=length(lx))
  for(i in 2:length(lx)) px[i-1]=(lx[i])/(lx[i-1])
  
  #fertility  (for birth-pulse, post-breeding census)
  bx=mx*px
  
  #projection matrix
  PX=px
  PX=PX[-length(PX)]
  BX=bx
  n=length(BX)
  Data=matrix(0,nrow=n,ncol=n)
  diag(Data)[-nrow(Data)]=PX
  Data=rbind(matrix(0,nrow=1,ncol=n),Data)
  Data=Data[-(n+1),]
  Data[1,]=BX
  rownames(Data)=colnames(Data)=(first.age+1):n
  
  #solve projection matrix
  LAMBDA=lambda(Data)
  r=log(LAMBDA)  
  
  
  return(r)  
}

#h
Stipns=function(M,age.mat,Meanfec,CyclE)
{  
  #survivorship
  surv=exp(-M)
  
  #fecundity  
  fecundity=Meanfec*sexratio/CyclE
  
  #maturity
  maturity=ifelse(age>=age.mat,1,0)   #knife edge
  
  # maximum age is plus group
  phi.o=0.0
  cum.survive=1.0
  z=0.0
  for (i in 2:(max.age)  )
  {
    z=M[i] + F.mult*Sel[i]
    z.ts=(M[i]+F.mult*Sel[i])*spawn.time
    phi.o=phi.o+cum.survive*fecundity[i]*maturity[i]*exp(-z.ts)
    cum.survive=cum.survive*exp(-z )
  }
  #plus group  
  z= M[max.age+1] + F.mult*Sel[max.age+1]
  z.ts=(M[max.age+1]+F.mult*Sel[max.age+1])*spawn.time
  phi.o=phi.o + fecundity[max.age+1]*maturity[max.age+1]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
  
  #maximum lifetime reproductive rate at low density
  alpha=phi.o*surv[1]
  
  #steepness
  h=alpha/(4+alpha)
  
  return(h)  
}


# Milk shark
first.age=0
max.age=12
age=first.age:max.age
age.mat=2
LINF=101
K=0.63
Temp=25
M.sim=M.fun(max.age,age.mat,LINF,K,Temp)
Age.mat.sim=2
Cycle.sim=1
sexratio=.5
Meanfec.sim=rep(7,length(age))
F.mult=0
Sel=rep(1,length(age))  #not used as F.mult set to 0
spawn.time=0

Leslie(M=M.sim,age.mat=Age.mat.sim,Meanfec=Meanfec.sim,CyclE=Cycle.sim)
Stipns(M=M.sim,age.mat=Age.mat.sim,Meanfec=Meanfec.sim,CyclE=Cycle.sim)


# Scalloped hammerhead
first.age=0
max.age=40
age=first.age:max.age
age.mat=14
LINF=301
K=0.07
Temp=27.3
M.sim=M.fun(max.age,age.mat,LINF,K,Temp)

Age.mat.sim=14
Cycle.sim=2
Meanfec.sim=rep(27,length(age))
Sel=rep(1,length(age))  #not used as F.mult set to 0
Leslie(M=M.sim,age.mat=Age.mat.sim,Meanfec=Meanfec.sim,CyclE=Cycle.sim)
Stipns(M=M.sim,age.mat=Age.mat.sim,Meanfec=Meanfec.sim,CyclE=Cycle.sim)