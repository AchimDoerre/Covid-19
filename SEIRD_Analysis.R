## --------------------------------------------------------------------------------------------------
##
## Implementation of SEIRD model for sex- and age-specific forecasting of COVID-19 epidemic
##
## Purpose of script: 
##
## Author: Achim Doerre
## ORCID: https://orcid.org/0000-0001-9297-3675
##
## Date Created: 1 February 2021
##
## Copyright (c) Achim Doerre, 2021
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
source("general_settings.R")
source("Model_core.R")


#-------------------------------------------------------------------
#General Settings
#-------------------------------------------------------------------

#number of random iterations
NR = 100

#forecasting horizon (in days)
daysmax = 75

#discrete time step size (for Forward Euler Method; 1/4 means quarter-days)
td = 1/4



#contact reduction at base lockdown (alpha = 0.2 means 80% general reduction of contacts)
alpha = 0.2


#set the day, on which basis the distribution of infected people is determined
#(choose from {200430, 201017, 201216})
choose_IncDist_reference_day = "200430"

#choose download date to base Icum0, Rcum0 and Dcum0 on
#(choose from {201216, 210106})
download_date = "210106"

#7-day incidence rate to start simulation from
base_7day_incidence = 7

#cumulative number of infections at beginning of forecast (hypothetical/supposed)
base_Icum = 3000000


scenario = 3


generate_parameters = TRUE


#-------------------------------------------------------------------
# Determination of random parameter combinations
#-------------------------------------------------------------------


if (generate_parameters){
  set.seed(1)
  #set ranges of potential true parameters
  epsilonrange = c(0.19, 0.21)
  gammarange = c(1/9, 1/7)
  manrange = c(0.65, 0.75)
  rrange = c(0.4, 0.5)
  
  PR = matrix(0, NR, 4)
  for (j in 1:NR){
    gammanr = runif(1, min(gammarange), max(gammarange))
    epsilonnr = runif(1, min(epsilonrange), max(epsilonrange))
    mannr = runif(1, min(manrange), max(manrange))
    rnr = runif(1, min(rrange), max(rrange))
    PR[j,] = c(gammanr, epsilonnr, mannr, rnr)
  }
  colnames(PR) = c("gammanr", "epsilonnr", "mannr", "rnr")
  write.table(PR, "Model_parameters.txt", row.names=FALSE)
}





#-------------------------------------------------------------------
# Fundamental settings
#-------------------------------------------------------------------

#Matrix of distribution of infected persons with respect to groups
IncDist = as.matrix(read.table("IcrdistM.txt", header=TRUE))
Idist0 = IncDist[choose_IncDist_reference_day,]

#number of simulated days
tmax = daysmax/td

#group names and number of groups
groups = c("f0-9", "f10-19", "f20-29", "f30-39", "f40-49", "f50-59", "f60-69", "f70-79", "f80+",
           "m0-9", "m10-19", "m20-29", "m30-39", "m40-49", "m50-59", "m60-69", "m70-79", "m80+")
ngroups = length(groups)

#contact matrix (from de Kassteele et al.; is adjusted for Covid-19 below)
CM = as.matrix(read.table("ContactRates.txt", header=T))



#-----------------------------------------------------------------------------
# Determining the initial Compartment settings
#-----------------------------------------------------------------------------

IcumM = as.matrix(read.table("GER_Icum_210105.txt", header=T))
Icum0 = IcumM[nrow(IcumM), 3:20]

RcumM = as.matrix(read.table("GER_Rcum_210105.txt", header=T))
Rcum0 = RcumM[nrow(RcumM), 3:20]

DcumM = as.matrix(read.table("GER_Dcum_210105.txt", header=T))
Dcum0 = DcumM[nrow(DcumM), 3:20]


# ad-hoc correction of reported numbers to adjust for delayed starting day of forecast ------
Idiff = base_Icum - sum(Icum0)
r = sum(Rcum0)/sum(Dcum0 + Rcum0)
Rdist0 = Rcum0/sum(Rcum0)
Ddist0 = Dcum0/sum(Dcum0)
Rcum0 = Rcum0 + r*Idiff*Rdist0
Dcum0 = Dcum0 + (1-r)*Idiff*Ddist0
Icum0 = Icum0 + Idiff*Idist0


#load age distribution of considered country (e.g., Germany)
A = read.table("agedistributionGermany10.txt", header=TRUE)
P0 = 1000*c(A[,"female"], A[,"male"])

#number of people in compartment S at beginning of simulation
S0 = P0 - Icum0

#number of people in compartment I at beginning of simulation
#(=distribution of infected people, adjusted for base 7-day incidence rate)
I0 = 83000000*base_7day_incidence/100000 * Idist0

#number of people in compartment E at beginning of simulation
#(= unknown number of exposed people)
#(approximated by number of infections in last 1/gamma = 5 days)
E0 = Idist0*sum(5/7*I0)

#number of people in compartments R and D at beginning of simulation
#(equal to cumulative numbers due to being absorbing states)
R0 = Rcum0
D0 = Dcum0




#-----------------------------------------------------------------------------
# Further parameter settings
#-----------------------------------------------------------------------------

#vector of day numbers
days = 1:daysmax


#secondary attack rate (determined externally)
wmin = 0.1324848*td


#group-combination-specific mitigation factor (=1 when no reduction)
mitfactor = matrix(1, 18, 18)


#general compliance drift rate (=1 when no temporal effect of compliance change)
general_compliance_drift = rep(1, length(days)/td)



#- scenario definitions ----------------

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


#scenario 4 (offices and schools reopen, female contacts as males)
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
mitfactor




#==========================================================================================
# Core calculation - central procedure
#==========================================================================================

result = Model_core(days, td=td, NR=NR,
                    Icum0 = Icum0,
                    S0, E0, I0, R0, D0,
                    CM = CM,
                    w = wmin,
                    mitfactor = mitfactor,
                    general_compliance_drift = general_compliance_drift,
                    plot_overall_course = TRUE)


#==========================================================================================





#------------------------------------------------------------------------------------------
# Extraction of Results
#------------------------------------------------------------------------------------------


#essential forecast
result$forecast

#to extract specific results, check names of output list:
names(result)

#explanation:
#Xm (Sm, Em, ...): mean number of people in compartment X at each time step 1, ..., tmax
#XR (SR, ER, ...): number of people in compartment X for each simulation 1, ..., NR (rows)
#                  and at each time step 1, ..., tmax (columns)
#IaR (I1R, I2R, I3R, ...): number of infectious people in respective population groups
#                          for each simulation 1, ..., NR (rows) and at each time step 1, ..., tmax (columns)
#DaR (D1R, D2R, D3R, ...): number of ceased people in respective population groups
#                          for each simulation 1, ..., NR (rows) and at each time step 1, ..., tmax (columns)
#forecast: essential forecast for each time step 1, ..., tmax; including cumulative infections, active infections,
#          number of ceased people
#M: essential information on outcome across all simulation runs (index of peak day of active infections,
#   maximum number of ICU patients, number of cumulative deaths)
