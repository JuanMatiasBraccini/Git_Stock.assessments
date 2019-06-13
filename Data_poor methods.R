
#1. Length-based Beverton-Holt Nonequilibrium Z Estimator

library(fishmethods)

#from: Gedamke, T. and J. M. Hoenig. 2006. Estimating mortality from mean length data in nonequilibrium
#     situations, with application to the assessment of goosefish. Trans. Am. Fish. Soc. 135:476-487

Z=bhnoneq(year=goosefish$year,mlen=goosefish$mlen, ss=goosefish$ss,
          K=0.108,Linf=126,Lc=30,nbreaks=1,styrs=c(1982),stZ=c(0.1,0.3),
          stsigma=20)



#2. Data-limited Methods  Toolkit (Carruthers & Hordyk)

library(DLMtool)
library(snowfall) # load package for parallel computing
sfInit(parallel=TRUE,cpus=2) # initiate the cluster with two cpus

#Example application to real fisheries data
mydata<-new('DLM_data') # create a new DLM data object and define:
mydata@Year<-2001:2010 # years
mydata@Cat<-matrix((11:20)*10*runif(10,0.5,1.5),nrow=1) # make up some annual catches
mydata@Ind<-matrix(seq(1.1,0.9,length.out=10)*runif(10,0.5,1.5),nrow=1)
mydata@Mort<-0.2 # instantaneous natural mortality rate
mydata@Abun<-1000 # current abundance estimate (biomass)
mydata@FMSY_M<-0.5 # guess of the ratio of FMSY to natural mortality rate
mydata@vbLinf<-200 # maximum length
mydata@vbK<-0.2 # von B growth coefficient k
mydata@LFC<-50 # length at first capture
mydata<-TAC(mydata) # calculate quotas
plot(mydata) # plot them
mydata<-Sense(mydata,"Fratio") # conduct a sensitivity analysis for one of the methods
sfStop()




