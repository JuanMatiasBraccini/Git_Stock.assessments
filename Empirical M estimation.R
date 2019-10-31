library(fishmethods)
M.empirical(Linf=30.1,Kl=0.31,TC=24,method=c(1))

library(TropFishR)
M_empirical(Linf = 80, K_l = 0.5, temp = 25, tmax = 30,
            method = c("Pauly_Linf","Hoenig"))