#  -------------- 1. SRAplus   --------------
#install.packages("remotes")
#remotes::install_github("james-thorson/FishLife@master")
remotes::install_github("DanOvando/sraplus")




#  -------------- CCSRA   --------------
#Catch curve stock reduction analysis with size composition in last year
#asumes logistic selectivity
#remotes::install_github("James-Thorson/CCSRA")
library(CCSRA)

#convert length to age composition
convert_length_to_age_samples(K, Linf, L0, Lcv, Lbin_mat, LengthComp_lt,
                              checkforbugs = TRUE)



library(fishmethods)

# 1 ----------------Mean length Beverton-Holt Nonequilibrium Z Estimator------------------------------

#from: Gedamke, T. and J. M. Hoenig. 2006. Estimating mortality from mean length data in nonequilibrium
#     situations, with application to the assessment of goosefish. Trans. Am. Fish. Soc. 135:476-487
data(goosefish)
Z=bhnoneq(year=goosefish$year,mlen=goosefish$mlen, ss=goosefish$ss,
          K=0.108,Linf=126,Lc=30,nbreaks=2,styrs=c(1982,1990),stZ=c(0.1,0.2,0.3))


# 2 ----------------Mortality Estimators from tagging------------------------------
## Data come from Appendix Table A2 and model structure from model (a) in
## Table 3.2 of Jiang (2005)
## Example takes a bit of time to run
data(Jiang)
model1<-irm_cr(relyrs = Jiang$relyrs, recapyrs = Jiang$recapyrs,
               N = Jiang$N, recapharv = Jiang$recapharv, recaprel = Jiang$recaprel,
               hlambda = Jiang$hlambda, rlambda = Jiang$rlambda, hphi = Jiang$hphi,
               rphi = Jiang$rphi, hmrate = Jiang$hmrate, Fyr = Jiang$Fyr,
               FAyr = Jiang$FAyr, Myr = Jiang$Myr, initial = c(0.1,0.05,0.1),
               lower = c(0.0001,0.0001,0.0001), upper=c(5,5,5),maxiter=10000)

# Data come from Table 4 and model structure from Table 5 under "year-specific F,
# constant M" in Hoenig et al. (1998)
data(Hoenig)
model2<-irm_h(relyrs = Hoenig$relyrs, recapyrs = Hoenig$recapyrs,
              N = Hoenig$N, recapharv = Hoenig$recapharv,lambda = Hoenig$lambda,
              phi = Hoenig$phi, Fyr = Hoenig$Fyr, Myr = Hoenig$Myr, initial = c(0.1,0.1),
              lower = c(0.0001,0.0001),upper = c(5,5), maxiter = 10000)


#Mortality using time series at large
data(tanaka)
model3=mort.al(relyr = tanaka$relyr, tal = tanaka$tal, N = tanaka$N)

#averaging the mortality models
tag_model_avg(model1,model2)



# -------------- 3. DBSRA Stock REduction (Dick and MAcCall 2011)   --------------
data(cowcod)
a=dbsra(year =cowcod$year, catch = cowcod$catch, catchCV = NULL, 
        catargs = list(dist="none",low=0,up=Inf,unit="MT"),
        agemat=11, k = list(low=100,up=15000,tol=0.01,permax=1000), 
        b1k = list(dist="none",low=0.01,up=0.99,mean=1,sd=0.1),
        btk = list(dist="beta",low=0.01,up=0.99,mean=0.1,sd=0.1,refyr=2009),
        fmsym = list(dist="lnorm",low=0.1,up=2,mean=-0.223,sd=0.2),
        bmsyk = list(dist="beta",low=0.05,up=0.95,mean=0.4,sd=0.05),
        M = list(dist="lnorm",low=0.001,up=1,mean=-2.90,sd=0.4),
        nsims = 10000)


library(TropFishR)

#  --------------Empirical M -------------
M_empirical(Linf = 80, K_l = 0.5, temp = 25, tmax = 30,
             method = c("Pauly_Linf","Hoenig","Then_growth","Then_tmax"))

# 3 ----------------Estimate Z based on a method derived by Beverton and Holt (1956)---------
#using length or age composition data

  # based on length-frequency data
data(synLFQ2)
Z_BevertonHolt(synLFQ2, catch_columns = 2, Lprime_tprime = 47.5)

  # based on age composition data
data(synCAA1)
Z_BevertonHolt(synCAA1, catch_columns = 3, Lprime_tprime = 2.5)


#4----------------Catch curve----
  # Variable paramter system (with catch vector) based on length frequency data
data(goatfish)
output <- catchCurve(goatfish)
summary(output$linear_mod)


  # Catch Curve with estimation of selection ogive
data(synLFQ3)
output <- catchCurve(synLFQ3, calc_ogive = TRUE)
summary(output$linear_mod_sel)

  # the same with predefined selection for regression line:
output <- catchCurve(synLFQ3, calc_ogive = TRUE, reg_int = c(9,21))
plot(output, plot_selec = TRUE)


#5----------------Selectivity----
# create list with selectivity information
select.list <- list(selecType = 'knife_edge',
                    Lc = 34, L75 = 37, tc = 5, meshSizes = c(60,80),
                    select_p1 = 2.7977, select_p2 = 0.1175)


# create vector with mid lengths
Lt <- seq(5, 50, 0.01)

# knife edge selectivity
sel_ke <- select_ogive(select.list, Lt)

# trawl ogive selectivity
select.list$selecType = "trawl_ogive"
sel_to <- select_ogive(select.list, Lt)

plot(Lt, sel_ke, type = 'l')
lines(Lt, sel_to, col = 'blue')

# Gillnet selectivity ("lognormal" and "normal_fixed")
select.list$selecType <- "lognormal"
sel_log <- select_ogive(select.list, Lt)

select.list$selecType <- "normal_fixed"
select.list$select_p1 <- 0.2
select.list$select_p2 <- 1.5
sel_nf <- select_ogive(select.list, Lt)

plot(Lt, sel_log, type = 'l')
lines(Lt, sel_nf, col = 'blue')



#Estimate selectivity from different experimental nets

  # Gillnet selectivity
data(tilapia)
out <- select(param = tilapia)
plot(out)


data(haddock)

output <- select_Millar(haddock, x0 = c(-10,0.3,0),
                        rtype = "tt.logistic")

plot(output, plotlens=seq(25,35,0.1), deviance_plot = FALSE)
legend("topleft",c("Control","Experimental"), lty=1:2, col=1:2)


# Gillnet
data(gillnet)

# Using inital estimates from old method
select_Millar(gillnet, x0 = NULL, rtype = "norm.loc")$value
select_Millar(gillnet, x0 = NULL, rtype = "norm.sca")$value
select_Millar(gillnet, x0 = NULL, rtype = "lognorm")$value


# 6 ----------LBSPR for estimating SPR based on size data and life history-------
library(LBSPR)

#Create objects
MyPars <- new("LB_pars")
slotNames(MyPars) #check what's in the object

#Populate the object
MyPars@Linf <- 100 
MyPars@L50 <- 66 
MyPars@L95 <- 70
MyPars@MK <- 1.5 

MyPars@SL50 <- 50 
MyPars@SL95 <- 65
MyPars@SPR <- 0.4
MyPars@BinWidth <- 5

#Run simulation
MySim <- LBSPRsim(MyPars)
MyPars@BinMax <- 150
MyPars@BinMin <- 0

#Outputs
MySim@SPR
MySim@FM  #ratio of F/M. the value for fishing mortality refers to the 
          # highest level of F experienced by any single size class

#Simulate specifying F/M instead of SPR
MyPars@SPR <- numeric() # remove value for SPR 
MyPars@FM <- 1 # set value for FM
MySim <- LBSPRsim(MyPars)
round(MySim@SPR, 2) # SPR at F/M = 1 

# Change the life history parameters
MyPars@MK <- 2.0 
MySim <- LBSPRsim(MyPars)
round(MySim@SPR, 2) # SPR 

MyPars@MK <- 0.5
MySim <- LBSPRsim(MyPars)
round(MySim@SPR, 2) # SPR 

MyPars@Linf <- 120
MySim <- LBSPRsim(MyPars)
round(MySim@SPR, 2) # SPR 

#Change selectivity parameters
MyPars@MK <- 1.5 
MyPars@SL50 <- 10
MyPars@SL95 <- 15 
MyPars@FM <- 1 
MySim <- LBSPRsim(MyPars)
round(MySim@SPR, 2) # SPR 

MyPars@SL50 <- 80
MyPars@SL95 <- 85 
MySim <- LBSPRsim(MyPars)
round(MySim@SPR, 2) # SPR 


#fitting empirical data
#need 2 objects: pars, and length data
# The length data can be either raw data or length frequency data
MyLengths <- new("LB_lengths") #The data type (freq or raw) must be 
                          # specified in the call to the new function
datdir <- DataDir()
list.files(datdir, pattern=".csv")

MyPars <- new("LB_pars")
# Note that only the life history parameters need to be specified for 
#            the estimation model. 
# The exploitation parameters will be estimated
MyPars@Species <- "MySpecies"
MyPars@Linf <- 100 
MyPars@L50 <- 66 
MyPars@L95 <- 70
MyPars@MK <- 1.5 
MyPars@L_units <- "mm"

#A ength freq data with multiple years
Len1 <- new("LB_lengths", LB_pars=MyPars, file=paste0(datdir, "/LFreq_MultiYr.csv"), 
            dataType="freq")
plotSize(Len1)
Len2 <- new("LB_lengths", LB_pars=MyPars, file=paste0(datdir, "/LFreq_MultiYrHead.csv"), 
            dataType="freq", header=TRUE)

#A raw length data set with multiple years
Len3 <- new("LB_lengths", LB_pars=MyPars, file=paste0(datdir, "/LRaw_MultiYr.csv"), 
            dataType="raw")
plotSize(Len3)

#Fit the model
myFit1 <- LBSPRfit(MyPars, Len1)
myFit3 <- LBSPRfit(MyPars, Len3)

myFit1@Ests
plotEsts(myFit1)
#By default the plotting function adds the smoother line to the estimated points



# 7 ----------------Data-limited Methods  Toolkit (Carruthers & Hordyk)------------------------------
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




