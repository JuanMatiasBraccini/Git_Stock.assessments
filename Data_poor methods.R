library(fishmethods)

# 1 ----------------Length-based Beverton-Holt Nonequilibrium Z Estimator------------------------------

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




library(TropFishR)
# 3 ----------------Estimate Z based on a method derived by Beverton and Holt (1956)---------
#using mide length or age composition data

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


# X ----------------Data-limited Methods  Toolkit (Carruthers & Hordyk)------------------------------
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




