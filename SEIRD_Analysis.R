## --------------------------------------------------------------------------------------------------
##
## Implementation of SEIRD model for sex- and age-specific forecasting of COVID-19 epidemic
##
## Purpose of script: 
##
## Author: Dr. Achim Doerre
## ORCID: https://orcid.org/0000-0001-9297-3675
##
## Date Created: 6 Oct 2020
##
## Copyright (c) Achim Doerre, 2020
## Email: achim.doerre@gmail.com
##
## --------------------------------------------------------------------------------------------------
##
## Notes:
## This script has been built for use in the original article
## 'Age- and Sex-Specific Modelling of the COVID-19 Epidemic'
## by Achim Doerre (ORCID: https://orcid.org/0000-0001-9297-3675 )
##  & Gabriele Doblhammer (ORCID: https://orcid.org/0000-0001-7746-0652 )
##
## The script has been edited and simplified for illustration of the method and for
## potentially adjusting it to similar setups and settings.
##   
##
## --------------------------------------------------------------------------------------------------



#set your working directory here
setwd("...")

#auxiliary functions and the core of the numerical implementation
source("general_functions.R")
source("Model_core.R")




#number of random iterations
NR = 100

#forecasting horizon (in days)
daysmax = 78
#discrete time step size (for Forward Euler Method; 1/4 means quarter-days)
td = 1/4




#group-specific fatality rates
dI = 1/100*rep(c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3), 2)
#group-specific hospitalisation fractions
ha = 1/100*rep(c(0.1, 0.3, 1.2, 3.2, 4.9, 10.2, 16.6, 24.3, 27.3), 2)
#group-specific intensive care fractions
ICUp = 1/100*rep(c(5.0, 5.0, 5.0, 5.0, 6.3, 12.2, 27.4, 43.2, 70.9), 2)



#color settings for different (explicit/derivative) compartments
shadecolS = rgb(0.75, 0.75, 0.75)
linecolS = rgb(0.2, 0.2, 0.2)

shadecolE = rgb(0.75, 0.75, 0.75)
linecolE = rgb(0.2, 0.2, 0.2)

shadecolI = rgb(1, 217/255, 179/255)
linecolI = rgb(1, 128/255, 0)

shadecolR = rgb(179/255, 1, 179/255)
linecolR = rgb(0.0, 153/255, 0.0)

linecolH = rgb(0.0, 0.0, 1.0)
shadecolH = rgb(179/255, 209/255, 255/255)

linecolICU = rgb(102/255, 0, 204/255)
shadecolICU = rgb(217/255, 179/255, 1)

shadecolD = rgb(1, 153/255, 153/255)
linecolD = rgb(204/255, 0, 0)




#number of discrete time steps
tmax = daysmax/td



#group names
groups = c("f0-9", "f10-19", "f20-29", "f30-39", "f40-49", "f50-59", "f60-69", "f70-79", "f80+",
           "m0-9", "m10-19", "m20-29", "m30-39", "m40-49", "m50-59", "m60-69", "m70-79", "m80+")
ngroups = length(groups)



#matrix of social contact rates
#here: 9 age groups x 2 sex = 18 groups, therefore 18x18 matrix
betaM = as.matrix(read.table("ContactMatrixMS18.txt", header=T))


#determining the initial cumulative and active compartment settings;
#example here: Germany on Aug 15 2020
Icum0 = c(4035, 9907, 16121, 16206, 19523, 17558, 9399, 8471, 13996,
          4361, 10342, 16740, 15988, 18589, 16718, 10369, 9345, 7915)

Rcum0 = c(3527, 9131, 15077, 15405, 18758, 16869, 8724, 7862, 10751,
          3846, 9430, 15460, 15016, 17671, 15891, 9090, 8192, 5056)

Dcum0 = c(1, 1, 3, 25, 48, 43, 481, 434, 3090,
          0, 4, 7, 68, 130, 117, 1085, 977, 2742)


#distribution of infections
Idist = Icum0/sum(Icum0)

#number of active cases at initial day
Iactive0 = Icum0 - Rcum0 - Dcum0



#dz is a factor to account for hypothetical missing cases due to under-reporting
#for no under-reporting, set
#dz = 0
#as a simple rule, we choose
dz = Iactive0/10



#initial settings for the compartments
#
#for S0, set the current population distribution across considered groups;
#here, we set the current German population minus individuals with infection history
S0 = 1000*c(3786, 3662, 4587, 5317, 5014, 6633, 5524, 4040, 3632,
            3974, 3875, 4986, 5623, 5091, 6698, 5217, 3416, 2274) - Icum0
#for the other compartments, we choose
E0 = dz
I0 = Iactive0
R0 = Rcum0
D0 = Dcum0





#wmin is the estimated secondary attack rate
#(separately estimated to be 13.2%)
wmin = 0.132*td


#mitigation factor to account for lockdown effect on contact reduction
#adjust accordingly with respect to number of groups and contact rate matrix
#
#mitfactor = 1 corresponds to no lockdown effect
mitfactor = matrix(1.0, ngroups, ngroups)


scenario = 1


#Scenario 1: general strict lockdown measures
if (scenario == 1){
  mitfactor = matrix(0.2, ngroups, ngroups)
  mitfactor1 = mitfactor
}

#Scenario 2: general strict lockdown measures, except offices reopen
#redfact is factor to account for newly increased contact rates
if (scenario == 2){
  redfact = 0.04
  mitfactor = matrix(0.2, ngroups, ngroups)
  mitfactor[4:6, 4:6] = 0.2 + redfact
  mitfactor[4:6, 7] = 0.2 + redfact/2
  mitfactor[7, 4:6] = 0.2 + redfact/2
  
  mitfactor[13:15, 13:15] = 0.2 + redfact
  mitfactor[13:15, 16] = 0.2 + redfact/2
  mitfactor[16, 13:15] = 0.2 + redfact/2
  
  mitfactor[4:6, 13:15] = 0.2 + redfact
  mitfactor[4:6, 16] = 0.2 + redfact/2
  mitfactor[7, 13:15] = 0.2 + redfact/2
  
  mitfactor[13:15, 4:6] = 0.2 + redfact
  mitfactor[13:15, 7] = 0.2 + redfact/2
  mitfactor[16, 4:6] = 0.2 + redfact/2
  mitfactor2 = mitfactor
}


#Scenario 2: general strict lockdown measures, except offices and schools reopen
#redfact is factor to account for newly increased contact rates
if (scenario == 3){
  redfact = 0.04
  mitfactor = matrix(0.2, ngroups, ngroups)
  mitfactor[2:3, 2:3] = 0.2 + redfact
  mitfactor[2:3, 4:6] = 0.2 + redfact/2
  mitfactor[4:6, 2:3] = 0.2 + redfact/2
  mitfactor[4:6, 4:6] = 0.2 + redfact
  mitfactor[4:6, 7] = 0.2 + redfact/2
  mitfactor[7, 4:6] = 0.2 + redfact/2
  
  mitfactor[11:12, 11:12] = 0.2 + redfact
  mitfactor[11:12, 13:15] = 0.2 + redfact/2
  mitfactor[11:15, 11:12] = 0.2 + redfact/2
  mitfactor[13:15, 13:15] = 0.2 + redfact
  mitfactor[13:15, 16] = 0.2 + redfact/2
  mitfactor[16, 13:15] = 0.2 + redfact/2
  
  mitfactor[2:3, 11:12] = 0.2 + redfact
  mitfactor[2:3, 13:15] = 0.2 + redfact/2
  mitfactor[4:6, 11:12] = 0.2 + redfact/2
  mitfactor[4:6, 13:15] = 0.2 + redfact
  mitfactor[4:6, 16] = 0.2 + redfact/2
  mitfactor[7, 13:15] = 0.2 + redfact/2
  
  mitfactor[11:12, 2:3] = 0.2 + redfact
  mitfactor[11:12, 4:6] = 0.2 + redfact/2
  mitfactor[13:15, 2:3] = 0.2 + redfact/2
  mitfactor[13:15, 4:6] = 0.2 + redfact
  mitfactor[13:15, 7] = 0.2 + redfact/2
  mitfactor[16, 4:6] = 0.2 + redfact/2
  mitfactor3 = mitfactor
}


colnames(mitfactor) = groups
rownames(mitfactor) = groups
round(mitfactor, 2)



result = ModelU_core(daysmax = daysmax, td=1/4, NR=NR,
                     S0, E0, I0, R0, D0,
                     betaM = betaM,
                     epsilonrange = c(0.19, 0.21),
                     gammarange = c(1/9, 1/7),
                     manrange = c(0.65, 0.75),
                     rrange = c(0.3, 0.6),
                     w = wmin,
                     mitfactor = mitfactor)


#results of the iterations
names(result)
tail(result$forecast)
tail(result$M)





#set logscale = 1 for displaying forecastings on a log scale
logscale = 0



#storing forecasting results in separate variables for convenience
#variables ending with 'R' contain the randomness
#variables ending with 'M' represent the average outcome with respect to the random repetitions
SR = result$SR
ER = result$ER
IR = result$IR
RR = result$RR
HR = result$HR
ICUR = result$ICUR
DR = result$DR

I1R = result$I1R
I2R = result$I2R
I3R = result$I3R
I4R = result$I4R
I5R = result$I5R
I6R = result$I6R
I7R = result$I7R
I8R = result$I8R
I9R = result$I9R
I10R = result$I10R
I11R = result$I11R
I12R = result$I12R
I13R = result$I13R
I14R = result$I14R
I15R = result$I15R
I16R = result$I16R
I17R = result$I17R
I18R = result$I18R


I1M = colMeans(I1R)
I2M = colMeans(I2R)
I3M = colMeans(I3R)
I4M = colMeans(I4R)
I5M = colMeans(I5R)
I6M = colMeans(I6R)
I7M = colMeans(I7R)
I8M = colMeans(I8R)
I9M = colMeans(I9R)
I10M = colMeans(I10R)
I11M = colMeans(I11R)
I12M = colMeans(I12R)
I13M = colMeans(I13R)
I14M = colMeans(I14R)
I15M = colMeans(I15R)
I16M = colMeans(I16R)
I17M = colMeans(I17R)
I18M = colMeans(I18R)


#approximate derivation of forecasted deaths
Ddist = Dcum0/sum(Dcum0)

D1R = Ddist[1]*result$DR
D2R = Ddist[2]*result$DR
D3R = Ddist[3]*result$DR
D4R = Ddist[4]*result$DR
D5R = Ddist[5]*result$DR
D6R = Ddist[6]*result$DR
D7R = Ddist[7]*result$DR
D8R = Ddist[8]*result$DR
D9R = Ddist[9]*result$DR
D10R = Ddist[10]*result$DR
D11R = Ddist[11]*result$DR
D12R = Ddist[12]*result$DR
D13R = Ddist[13]*result$DR
D14R = Ddist[14]*result$DR
D15R = Ddist[15]*result$DR
D16R = Ddist[16]*result$DR
D17R = Ddist[17]*result$DR
D18R = Ddist[18]*result$DR



D1M = colMeans(D1R)
D2M = colMeans(D2R)
D3M = colMeans(D3R)
D4M = colMeans(D4R)
D5M = colMeans(D5R)
D6M = colMeans(D6R)
D7M = colMeans(D7R)
D8M = colMeans(D8R)
D9M = colMeans(D9R)
D10M = colMeans(D10R)
D11M = colMeans(D11R)
D12M = colMeans(D12R)
D13M = colMeans(D13R)
D14M = colMeans(D14R)
D15M = colMeans(D15R)
D16M = colMeans(D16R)
D17M = colMeans(D17R)
D18M = colMeans(D18R)







#Forecasting plots
par(mar=c(2,2,2,2), new=FALSE)
options(scipen=5)
propmin = 0.7
propmax = 1.3

#Panel 1 (total)
IS = IR
DS = DR
#pdf("../figures/S1P1.pdf", width=116, height=77)
par(mar=c(2,2,2,2), new=FALSE)
forecast_band2(IS, linecolI, shadecolI, DS, linecolD, shadecolD, "Total")


x = 1:ncol(IS)
y = (I10M+I11M+I12M+I13M+I14M+I15M+I16M+I17M+I18M)/(I1M+I2M+I3M+I4M+I5M+I6M+I7M+I8M+I9M)
par(new = TRUE)
plot(x, y, typ="l", lty=2, axes = FALSE, xlab = "", ylab = "", ylim=c(propmin, propmax))
axis(side = 4, at = seq(propmin, propmax, 0.1))



#forecast for working population (ages 30-64)
IS = I4R+I5R+I6R+0.5*I7R + I13R+I14R+I15R+0.5*I16R
DS = D4R+D5R+D6R+0.5*D7R + D13R+D14R+D15R+0.5*D16R
forecast_band2(IS, linecolI, shadecolI, DS, linecolD, shadecolD, "Ages 30-64")


x = 1:ncol(IS)
y = (I13M+I14M+I15M+0.5*I16M)/(I4M+I5M+I6M+0.5*I7M)
par(new = TRUE)
plot(x, y, typ="l", lty=2, axes = FALSE, xlab = "", ylab = "", ylim=c(propmin, propmax))
axis(side = 4, at = seq(propmin, propmax, 0.1))

#END OF SCRIPT