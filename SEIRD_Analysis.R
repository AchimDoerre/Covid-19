## --------------------------------------------------------------------------------------------------
##
## Implementation of SEIRD model for sex-/gender- and age-specific assessment of COVID-19 epidemic
##
## Purpose of script: 
##
## Author: Dr. Achim Doerre
## ORCID: https://orcid.org/0000-0001-9297-3675
##
## Date Created: 25 Oct 2021
##
## Copyright (c) Achim Doerre, 2021
## Email: achim.doerre@gmail.com
##
## --------------------------------------------------------------------------------------------------
##
## Notes:
## This script has been built for use in the original article
## 'The influence of gender on COVID-19 infections and mortality in Germany: insights from age- and gender-specific
## modeling of contact rates, infections, and deaths in the early phase of the pandemic'
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
source("RCode/general_functions.R")
source("RCode/general_settings.R")
source("RCode/Model_core.R")


scenario = 1

#number of random iterations
NR = 1000

#forecasting horizon (in days)
daysmax = 75

#discrete time step size (for Forward Euler Method; 1/4 means quarter-days)
td = 1/4

#number of discrete time steps
tmax = daysmax/td

#alpha is the proportion of remaining contact intensity during strict lockdown
alpha = 0.2



#group-specific fatality rates
dI = 1/100*rep(c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3), 2)
#group-specific hospitalisation fractions (per infected)
ha = 1/100*rep(c(0.1, 0.3, 1.2, 3.2, 4.9, 10.2, 16.6, 24.3, 27.3), 2)
#group-specific intensive care fractions (per hospitalised)
ICUp = 1/100*rep(c(5.0, 5.0, 5.0, 5.0, 6.3, 12.2, 27.4, 43.2, 70.9), 2)




#group names
groups = c("f0-9", "f10-19", "f20-29", "f30-39", "f40-49", "f50-59", "f60-69", "f70-79", "f80+",
           "m0-9", "m10-19", "m20-29", "m30-39", "m40-49", "m50-59", "m60-69", "m70-79", "m80+")
ngroups = length(groups)


#matrix of social contact rates
#here: 9 age groups x 2 sex = 18 groups, therefore 18x18 matrix
CM = as.matrix(read.table("Data/ContactRates.txt", header=T))


#determining the initial cumulative and active compartment settings;
#choice here: Germany on Jan 5 2021
Icum0 = c(44149, 93731, 142086, 132923, 152165, 136865, 71863, 64785, 112766,
          47722, 95147, 140806, 124385, 136078, 122399, 70696, 63733, 56349)

Rcum0 = c(38039, 79380, 119445, 110116, 124666, 112128, 54750, 49356, 68119,
          41155, 81728, 120678, 105040, 113495, 102083, 52476, 47307, 30571)

Dcum0 = c(10, 8, 13, 84, 159, 143, 1778, 1603, 13653,
          5, 12, 20, 197, 378, 340, 3499, 3153, 11472)

#observed distribution of infected
Idist0 = Icum0/sum(Icum0)

#7-day incidence rate to start simulation from
base_7day_incidence = 7

#cumulative number of infections at beginning of forecast
base_Icum = 3000000

#correction of reported numbers to adjust for starting day of forecast
Idiff = base_Icum - sum(Icum0)
r = sum(Rcum0)/sum(Dcum0 + Rcum0)
Rdist0 = Rcum0/sum(Rcum0)
Ddist0 = Dcum0/sum(Dcum0)
Rcum0 = Rcum0 + r*Idiff*Rdist0
Dcum0 = Dcum0 + (1-r)*Idiff*Ddist0

#Iactive0 = Icum0 - Rcum0 - Dcum0
#sum(Iactive0)

P0 = 1000*c(3786, 3662, 4587, 5317, 5014, 6633, 5524, 4040, 3646,
            3974, 3875, 4986, 5623, 5091, 6698, 5217, 3416, 2277)

#S0: number of people in compartment S at beginning of simulation
S0 = P0 - Icum0

#I0: number of people in compartment I at beginning of simulation
#(=distribution of infected people, adjusted for base 7-day incidence rate)
I0 = 83000000*base_7day_incidence/100000 * Idist0

#E0: number of people in compartment E at beginning of simulation
#(= unknown number of exposed people)
#(approximated by number of infections in last 1/gamma = 5 days)
E0 = Idist0*sum(5/7*I0)

#R0 and D0: number of people in compartments R and D at beginning of simulation
#(equal to cumulative numbers due to being absorbing states)
R0 = Rcum0
D0 = Dcum0



#ICU0: number of people in ICU at beginning of simulation
ICU0 = (ICUp*ha)/sum(ICUp*ha)*2000

#Iactive0: number of active cases at initial day
Iactive0 = Icum0 - Rcum0 - Dcum0


#dz is a factor to account for hypothetical missing cases due to under-reporting
#for no under-reporting, set dz = 0;
#as a simple rule, we suggest:
dz = Iactive0/10


#days: vector of day numbers
days = 1:daysmax


#wmin: secondary attack rate (determined externally and adjusted for time step size)
wmin = 0.132*td


#mitfactor (mitigation factor): variable to account for lockdown effect on contact reduction
#adjust accordingly with respect to number of groups and contact rate matrix
#mitfactor = matrix(1.0, ngroups, ngroups) corresponds to no lockdown effect (i.e., scenario=0)
#note that mitfactor defines the mitigation with respect to all *pairs* of population groups
mitfactor = matrix(1.0, ngroups, ngroups)


#general compliance drift rate (=1 when no temporal effect of compliance change)
general_compliance_drift = rep(1, length(days)/td)



#- Scenario Definitions ---------------------------------------------------------------------------
#note that the scenario definitions are tailored to 9x2=18 groups (9 age categories, 2 genders)
#and would require appropriate adjustment in case of other models or datasets
#
#redfact: reduction factor with respect to lockdown impact (due to easing measures)

#scenario 1 (general lockdown measures)
if (scenario == 1){
  mitfactor = matrix(alpha, 18, 18)
}

#scenario 2 (offices reopen)
if (scenario == 2){
  redfact = 0.05
  mitfactor = matrix(alpha, 18, 18)
  mitfactor[4:6, 4:6] = alpha + redfact
  mitfactor[4:6, 7] = alpha + redfact/2
  mitfactor[7, 4:6] = alpha + redfact/2
  
  mitfactor[13:15, 13:15] = alpha + redfact
  mitfactor[13:15, 16] = alpha + redfact/2
  mitfactor[16, 13:15] = alpha + redfact/2
  
  mitfactor[4:6, 13:15] = alpha + redfact
  mitfactor[4:6, 16] = alpha + redfact/2
  mitfactor[7, 13:15] = alpha + redfact/2
  
  mitfactor[13:15, 4:6] = alpha + redfact
  mitfactor[13:15, 7] = alpha + redfact/2
  mitfactor[16, 4:6] = alpha + redfact/2
}



#scenario 3 (offices and schools reopen)
if (scenario == 3){
  redfact = 0.05
  mitfactor = matrix(alpha, 18, 18)
  mitfactor[2:3, 2:3] = alpha + redfact
  mitfactor[2:3, 4:6] = alpha + redfact/2
  mitfactor[4:6, 2:3] = alpha + redfact/2
  mitfactor[4:6, 4:6] = alpha + redfact
  mitfactor[4:6, 7] = alpha + redfact/2
  mitfactor[7, 4:6] = alpha + redfact/2
  
  mitfactor[11:12, 11:12] = alpha + redfact
  mitfactor[11:12, 13:15] = alpha + redfact/2
  mitfactor[11:15, 11:12] = alpha + redfact/2
  mitfactor[13:15, 13:15] = alpha + redfact
  mitfactor[13:15, 16] = alpha + redfact/2
  mitfactor[16, 13:15] = alpha + redfact/2
  
  mitfactor[2:3, 11:12] = alpha + redfact
  mitfactor[2:3, 13:15] = alpha + redfact/2
  mitfactor[4:6, 11:12] = alpha + redfact/2
  mitfactor[4:6, 13:15] = alpha + redfact
  mitfactor[4:6, 16] = alpha + redfact/2
  mitfactor[7, 13:15] = alpha + redfact/2
  
  mitfactor[11:12, 2:3] = alpha + redfact
  mitfactor[11:12, 4:6] = alpha + redfact/2
  mitfactor[13:15, 2:3] = alpha + redfact/2
  mitfactor[13:15, 4:6] = alpha + redfact
  mitfactor[13:15, 7] = alpha + redfact/2
  mitfactor[16, 4:6] = alpha + redfact/2
}


#scenario 4 (offices and schools reopen, female contacts same as males)
if (scenario == 4){
  redfact = 0.05
  mitfactor = matrix(alpha, 18, 18)
  mitfactor[2:3, 2:3] = alpha + redfact
  mitfactor[2:3, 4:6] = alpha + redfact/2
  mitfactor[4:6, 2:3] = alpha + redfact/2
  mitfactor[4:6, 4:6] = alpha + redfact
  mitfactor[4:6, 7] = alpha + redfact/2
  mitfactor[7, 4:6] = alpha + redfact/2
  
  mitfactor[11:12, 11:12] = alpha + redfact
  mitfactor[11:12, 13:15] = alpha + redfact/2
  mitfactor[11:15, 11:12] = alpha + redfact/2
  mitfactor[13:15, 13:15] = alpha + redfact
  mitfactor[13:15, 16] = alpha + redfact/2
  mitfactor[16, 13:15] = alpha + redfact/2
  
  mitfactor[2:3, 11:12] = alpha + redfact
  mitfactor[2:3, 13:15] = alpha + redfact/2
  mitfactor[4:6, 11:12] = alpha + redfact/2
  mitfactor[4:6, 13:15] = alpha + redfact
  mitfactor[4:6, 16] = alpha + redfact/2
  mitfactor[7, 13:15] = alpha + redfact/2
  
  mitfactor[11:12, 2:3] = alpha + redfact
  mitfactor[11:12, 4:6] = alpha + redfact/2
  mitfactor[13:15, 2:3] = alpha + redfact/2
  mitfactor[13:15, 4:6] = alpha + redfact
  mitfactor[13:15, 7] = alpha + redfact/2
  mitfactor[16, 4:6] = alpha + redfact/2
  
  CM[1:9,1:9] = CM[10:18,10:18]
  CM[1:9,10:18] = CM[1:9,10:18]
}


colnames(mitfactor) = groups
rownames(mitfactor) = groups



#==========================================================================================
# Core calculation - central procedure
#==========================================================================================

result = Model_core(days, td=td, NR=NR,
                     Icum0 = Icum0,
                     ICU0 = ICU0,
                     S0, E0, I0, R0, D0,
                     CM = CM,
                     w = wmin,
                     mitfactor = mitfactor,
                     general_compliance_drift = general_compliance_drift,
                     plot_overall_course = TRUE)


#==========================================================================================


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





#Forecast Plots
propmin = 0.7
propmax = 1.3

par(mar=c(4,2,2,2), new=FALSE)
options(scipen=5)

forecast_band(IR, linecolI, shadecolI, "Total Infected")

#gender ratio
x = 1:ncol(IR)
y = (I10M+I11M+I12M+I13M+I14M+I15M+I16M+I17M+I18M)/(I1M+I2M+I3M+I4M+I5M+I6M+I7M+I8M+I9M)
par(new = TRUE)
plot(x, y, typ="l", lty=1, axes = FALSE, xlab = "", ylab = "", ylim=c(propmin, propmax), lwd=2)
axis(side = 4, at = seq(propmin, propmax, 0.1))
abline(h=1, lty=2)


#END OF SCRIPT