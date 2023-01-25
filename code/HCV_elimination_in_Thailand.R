##############################
# HEP C AGE STRUCTURED MODEL #
##############################

# Cleanup - CAUTION clears R environment ####
rm(list=ls())

# Setup - load packages and define plotting functions ####

#setwd('~/MGH/thailand/dissertation/submission/code')

library(pacman)
p_load(deSolve, tidyverse, doParallel, manipulate, readxl, gridExtra, grid, scales, tictoc, Hmisc, viridis, rlist)

scaleFUN <- function(x) sprintf("%.2f", x)
scaleFUN2 <- function(x) sprintf("%.4f", x)

age_group_vector<-c('0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 -39',
                    '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69', '70 - 74', '75 - 79',
                    '80 - 84', '85 - 89', '90 - 94', '95 - 99', '100 and over')

get_only_legend <- function(plot) {
  # get tabular interpretation of plot
  plot_table <- ggplot_gtable(ggplot_build(plot)) 
  #  Mark only legend in plot
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
  # extract legend
  legend <- plot_table$grobs[[legend_plot]]
  # return legend
  return(legend) 
}

`%nin%` = Negate(`%in%`)

# Read Excel data ####


# Note for those getting code from github:
# Please create an R Project in the same directory as:
# - a directory called 'data' containing the excel file HCV_screening_data.xlsx
# - this model file

# HCV deaths data
hcv_deaths_data <- as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="HCV_deaths", range="A1:D9", col_names=TRUE)) #%>% 

# Read mortality and population structure data

scenario_base <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="mortality_scenarios", range="C3:W3", col_names=FALSE)))
scenario_1 <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="mortality_scenarios", range="C4:W4", col_names=FALSE)))
scenario_2 <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="mortality_scenarios", range="C5:W5", col_names=FALSE)))
scenario_3 <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="mortality_scenarios", range="C6:W6", col_names=FALSE)))

# Contact matrix data
contact0 <- as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="sexual_contact_matrix", range="D3:X23", col_names=FALSE))

# Initial conditions values
init_condits <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="initial_conditions_new", range="B5:V31", col_names=FALSE)))

# Birth rate values
birth.approx <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="birthrate", range="A1:B38", col_names=TRUE)))

# Reading mortality baseline (empty from 2023 onwards)
mortality_base_empty <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="mortality_rates", range="B3:W40", col_names=TRUE)))

# Population structure data
age_struc <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="Population_data_proportion", range="A2:V20", col_names=TRUE)))

# Population and birth rate data and projection - UN
pop_data <- cbind(c(2004:2021),as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="Population_data_absolute", range="W3:W20", col_names=FALSE))))
pop_proj <- cbind(c(2022:2040),as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="Population_data_absolute", range="W21:W39", col_names=FALSE))))

pop_data1 <- cbind(c(2004:2040),as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="Population_data_absolute", range="W3:X39", col_names=FALSE))))

# Population structure data and projection
age_struc_proportion_data <- cbind(2004:2021,as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="Population_data_absolute", range="B2:V20", col_names=TRUE))))
age_struc_proportion_proj <- cbind(2004:2040,as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="Population_data_absolute", range="B2:V39", col_names=TRUE))))

# Prevalence data
data_prev <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="prevalence_data_by_age", range="B1:G22", col_names=TRUE)))

# Simulations
simulations <- as_tibble(as.data.frame(read_excel("data/HCV_screening_data.xlsx", sheet="all_simulations", range="A1:F200", col_names=TRUE)))

# Aging process ####

age_groups <-  c(seq(1,21,length.out=21))
groups <- length(age_groups)
da <- rep(5,groups) # difference between consecutive age groups in Years

aging <- diag(-1/da)
aging[row(aging)-col(aging)==1] <- 1/head(da,-1) # time in years to enter and leave each age group 


# Contact matrix ####

contact0[contact0==0] <- min(contact0[contact0!=0])/50

beta <- matrix(NA,nrow=groups,ncol=groups)

beta_multiplier <- 0.025

for(j in 1:groups){
  for(i in 1:groups){
    beta[i,j] <- beta_multiplier*contact0[i,j]
  }
}


# Screening by age group ####

scr_start <- 19 # screening programme start year: 19=2023
scr_dur <- 7 # screening programme duration (years)

# Baseline screening coverage

cov_mean <- 0.06184121 # From model fit
cov_sd <- 0.01220783 # From model fit
cov_lower <- cov_mean - 1.96*cov_sd
cov_upper <- cov_mean + 1.96*cov_sd

groupBASE <- rep(cov_mean,groups)
groupBASE_lower <- rep(cov_lower,groups)
groupBASE_upper <- rep(cov_upper,groups)


# Targeted screening strategies 

coverageBASE <- cov_mean # ~6.2%
coverage1 <- 0.15
coverage2 <- 0.25
coverage3 <- 0.5

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

# Set up empty indices for storing results ####

Sindex <- 1:groups
F0index <- (groups+1):(2*groups)
F1index <- (2*groups+1):(3*groups)
F2index <- (3*groups+1):(4*groups)
F3index <- (4*groups+1):(5*groups)
C1index <- (5*groups+1):(6*groups)
C2index <- (6*groups+1):(7*groups)
C3index <- (7*groups+1):(8*groups)
C4index <- (8*groups+1):(9*groups)

HCCAindex <- (9*groups+1):(10*groups)
HCCBindex <- (10*groups+1):(11*groups)
HCCCindex <- (11*groups+1):(12*groups)
HCCDindex <- (12*groups+1):(13*groups)

F0cureindex <- (13*groups+1):(14*groups)
F1cureindex <- (14*groups+1):(15*groups)
F2cureindex <- (15*groups+1):(16*groups)
F3cureindex <- (16*groups+1):(17*groups)
C1cureindex <- (17*groups+1):(18*groups)
C2cureindex <- (18*groups+1):(19*groups)
C3cureindex <- (19*groups+1):(20*groups)
C4cureindex <- (20*groups+1):(21*groups)

Dindex <- (21*groups+1):(22*groups)
dthC14index <- (22*groups+1):(23*groups)
dthHCCindex <- (23*groups+1):(24*groups)

translivindex <- (24*groups+1):(25*groups)

aliveindex <- (1:21*groups)
infectindex <- (groups+1):(13*groups)
totalHCCindex <- (9*groups+1):(13*groups)
totalHCVindex <- (groups+1):(9*groups)

newdeathindex <- (25*groups+1):(26*groups)
incHCCindex <- (26*groups+1):(27*groups)

CIncindex <- (27*groups+1):(28*groups)
CScrindex <- (28*groups+1):(29*groups)

infectindex <- (groups+1):(13*groups)

totalHCCindex <- (9*groups+1):(13*groups)
totalHCVindex <- (groups+1):(9*groups)


# Parameters ####

start_year <- 2004
# Compartment flow rates
f0f1 <- 0.117
f1f2 <- 0.085
f2f3 <- 0.12
f3c1 <- 0.116
c1c2 <- 0.044
c2c3 <- 0.044
c3c4 <- 0.076
c1bA <- 0.0068
c1bB <- 0.0099
c1bC <- 0.0029
c1bD <- 0.0068
c2bA <- 0.0068
c2bB <- 0.0099
c2bC <- 0.0029
c2bD <- 0.0068
c3bD <- 0.0664
c4bD <- 0.0664
# Compartment mortality rates
dthc1 <- 0.01 
dthc2 <- 0.01
dthc3 <- 0.2
dthc4 <- 0.57
dthbA <- 1/(36/12) # 3 years to die from HCCA
dthbB <- 1/(16/12) # 1.3 years to die from HCCB
dthbC <- 1/(6/12) # 6 months to die from HCCC
dthbD <- 1/(3/12) # 3 months to die from HCCD
dthtrn <- 1/(240/12)
tranc4 <- 0.0015 # liver transplant rates
tranbA <- 0.0015
tranbB <- 0.0015
trt_start <- 15 # DAA treatment starting in 2019
sens <- 0.985 # sensitivity - combined

pF0scr<-1 # proportion of compartment screened
pF1scr<-1
pF2scr<-1
pF3scr<-1
pC1scr<-1
pC2scr<-1
pC3scr<-1
pC4scr<-1

std_cureF0<-0.7 # efficacy of standard treatment
std_cureF1<-0.7
std_cureF2<-0.7
std_cureF3<-0.7
std_cureC1<-0.7
new_cureF0<-0.985
new_cureF1<-0.985 # efficacy of new treatment - proportion cured
new_cureF2<-0.985
new_cureF3<-0.985
new_cureC1<-0.985
new_cureC2<-0.985
new_cureC3<-0.985
new_cureC4<-0.985

# set up parameters
parameters<- c(
  start_year,
  f0f1,
  f1f2,
  f2f3,
  f3c1,
  c1c2,
  c2c3,
  c3c4,
  c1bA,
  c1bB,
  c1bC,
  c1bD,
  c2bA,
  c2bB,
  c2bC,
  c2bD,
  c3bD,
  c4bD,
  dthc1,
  dthc2,
  dthc3,
  dthc4,
  dthbA,
  dthbB,
  dthbC,
  dthbD,
  dthtrn,
  tranc4,
  tranbA,
  tranbB,
  trt_start,
  sens,
  pF0scr,
  pF1scr,
  pF2scr,
  pF3scr,
  pC1scr,
  pC2scr,
  pC3scr,
  pC4scr,
  std_cureF0,
  std_cureF1,
  std_cureF2,
  std_cureF3,
  std_cureC1,
  new_cureF0,
  new_cureF1,
  new_cureF2,
  new_cureF3,
  new_cureC1,
  new_cureC2,
  new_cureC3,
  new_cureC4)

# Age dependent transition vectors. Currently uniform across age groups
f0f1_vec <- rep(f0f1,groups)
f3c1_vec <- rep(f3c1,groups)

# Set up time ####

simu.time <- seq(0, 36, by=1) 

# Initial conditions ####

initS <- as.numeric(init_condits[1,]) # S row

initF0 <- as.numeric(init_condits[2,]) # F0 row
initF1 <- as.numeric(init_condits[3,])
initF2 <- as.numeric(init_condits[4,])
initF3 <- as.numeric(init_condits[5,])

initC1 <- as.numeric(init_condits[6,])
initC2 <- as.numeric(init_condits[7,])
initC3 <- as.numeric(init_condits[8,])
initC4 <- as.numeric(init_condits[9,])

initHCCA <- as.numeric(init_condits[10,])
initHCCB <- as.numeric(init_condits[11,])
initHCCC <- as.numeric(init_condits[12,])
initHCCD <- as.numeric(init_condits[13,])

initF0cure <- as.numeric(init_condits[14,])
initF1cure <- as.numeric(init_condits[15,])
initF2cure <- as.numeric(init_condits[16,])
initF3cure <- as.numeric(init_condits[17,])

initC1cure <- as.numeric(init_condits[18,])
initC2cure <- as.numeric(init_condits[19,])
initC3cure <- as.numeric(init_condits[20,])
initC4cure <- as.numeric(init_condits[21,])

# Set up initial death
initD <- as.numeric(init_condits[22,])
initdthC14 <- as.numeric(init_condits[23,])
initdthHCC <- as.numeric(init_condits[24,])

initDHCC <- as.numeric(init_condits[25,])
initDC14 <- as.numeric(init_condits[26,])
inittransliv <- as.numeric(init_condits[27,])

initCInc <- rep(0,groups)
initCScr <- rep(0,groups)

# Set up initial values

init <- c(S=initS,F0=initF0,F1=initF1,F2=initF2,F3=initF3,C1=initC1,C2=initC2,C3=initC3,C4=initC4,
          HCCA=initHCCA,HCCB=initHCCB,HCCC=initHCCC,HCCD=initHCCD,
          F0cure=initF0cure,F1cure=initF1cure,F2cure=initF2cure,F3cure=initF3cure,
          C1cure=initC1cure,C2cure=initC2cure,C3cure=initC3cure,C4cure=initC4cure,
          D=initD,dthC14=initdthC14,dthHCC=initdthHCC,transliv=inittransliv, CInc=initCInc, CScr=initCScr)



# Define model ####
hepC.mod<- function(t, y, parameters)
{
  with(as.list(c(y, parameters)),
       {
         
         S <- y[Sindex]
         F0 <- y[F0index]
         F1 <- y[F1index]
         F2 <- y[F2index]
         F3 <- y[F3index]
         C1 <- y[C1index]
         C2 <- y[C2index]
         C3 <- y[C3index]
         C4 <- y[C4index]
         
         HCCA <- y[HCCAindex]
         HCCB <- y[HCCBindex]
         HCCC <- y[HCCCindex]
         HCCD <- y[HCCDindex]
         
         F0cure <- y[F0cureindex]
         F1cure <- y[F1cureindex]
         F2cure <- y[F2cureindex]
         F3cure <- y[F3cureindex]
         C1cure <- y[C1cureindex]
         C2cure <- y[C2cureindex]
         C3cure <- y[C3cureindex]
         C4cure <- y[C4cureindex]
         
         D <- y[Dindex]
         dthC14 <- y[dthC14index]
         dthHCC <- y[dthHCCindex]
         
         transliv <- y[translivindex]
         
         alive <- y[aliveindex]
         infect <- y[infectindex]
         totalHCC <- y[totalHCCindex]
         totalHCV <- y[totalHCVindex]
         
         newdeath <- y[newdeathindex]
         incHCC <- y[incHCCindex]
         
         dCInc <- y[CIncindex]
         
         dCScr <- y[CScrindex]
         
         
         # Screening
         
         if((t>scr_start)&(t<scr_start+scr_dur)){scr<-scr_group}else{scr<-groupBASELINE}
         
         alive <- S+F0+F1+F2+F3+C1+C2+C3+C4+HCCA+HCCB+HCCC+HCCD+F0cure+F1cure+F2cure+F3cure+C1cure+C2cure+C3cure+C4cure
         pop <- sum(alive)
         
         infect <- (F0+F1+F2+F3+C1+C2+C3+C4+HCCA+HCCB+HCCC+HCCD)
         lambda <- beta%*%infect/pop
         totalHCC <- HCCA+HCCB+HCCC+HCCD
         totalHCV  <- F0+F1+F2+F3+C1+C2+C3+C4
         prevalence <- 100*(infect/pop)
         
         
         flowin <- birth.func(t)*pop
         birth <- c(flowin, rep(0,groups-1))
         
         cureF0 <- ifelse(t<=trt_start, std_cureF0, new_cureF0) # efficacy of old and new treatment
         cureF1 <- ifelse(t<=trt_start, std_cureF1, new_cureF1)
         cureF2 <- ifelse(t<=trt_start, std_cureF2, new_cureF2)
         cureF3 <- ifelse(t<=trt_start, std_cureF3, new_cureF3)
         cureC1 <- ifelse(t<=trt_start, std_cureC1, new_cureC1)
         cureC2 <- ifelse(t<=trt_start, 0, new_cureC2)
         cureC3 <- ifelse(t<=trt_start, 0, new_cureC3)
         cureC4 <- ifelse(t<=trt_start, 0, new_cureC4)
         
         natdeath <- c(mortality_func1(t),mortality_func2(t),mortality_func3(t),mortality_func4(t),mortality_func5(t),mortality_func6(t),mortality_func7(t),mortality_func8(t),mortality_func9(t),mortality_func10(t),
                       mortality_func11(t),mortality_func12(t),mortality_func13(t),mortality_func14(t),mortality_func15(t),mortality_func16(t),mortality_func17(t),mortality_func18(t),mortality_func19(t),mortality_func20(t),mortality_func21(t))
         
         # ODEs
         dS <- birth -lambda*S -natdeath*S + aging %*% S
         dF0 <- ifelse(F0>=0,-f0f1_vec*F0 + lambda*S - scr*pF0scr*sens*cureF0*F0 -natdeath*F0 + aging %*% F0,lambda*S + aging %*% F0)
         dF1 <- ifelse(F1>=0,f0f1_vec*F0 -f1f2*F1 - scr*pF1scr*sens*cureF1*F1 -natdeath*F1 + aging %*% F1,rep(0,groups))
         dF2 <- ifelse(F2>=0,f1f2*F1 -f2f3*F2 - scr*pF2scr*sens*cureF2*F2  -natdeath*F2 + aging %*% F2,rep(0,groups))
         dF3 <- ifelse(F3>=0,f2f3*F2 -f3c1_vec*F3 - scr*pF3scr*sens*cureF3*F3 -natdeath*F3 + aging %*% F3,rep(0,groups))
         dC1 <- ifelse(C1>=0,f3c1_vec*F3 -dthc1*C1 -c1c2*C1 -scr*pC1scr*sens*cureC1*C1 -(c1bA + c1bB + c1bC + c1bD)*C1 -natdeath*C1 + aging %*% C1,rep(0,groups))
         dC2 <- ifelse(C2>=0,c1c2*C1 -dthc2*C2 -c2c3*C2 - scr*pC2scr*sens*cureC2*C2 -(c2bA + c2bB + c2bC + c2bD)*C2 -natdeath*C2 + aging %*% C2,rep(0,groups))
         dC3 <- ifelse(C3>=0,c2c3*C2 -dthc3*C3 -c3c4*C3 - scr*pC3scr*sens*cureC3*C3 -c3bD*C3 -natdeath*C3 + aging %*% C3,rep(0,groups))
         dC4 <- ifelse(C4>=0,c3c4*C3 -dthc4*C4 - scr*pC4scr*sens*cureC4*C4 - c4bD*C4 - tranc4*C4 -natdeath*C4 + aging %*% C4,rep(0,groups))
         
         dHCCA <- c1bA*(C1+C1cure) + c2bA*(C2+C2cure) -dthbA*HCCA -tranbA*HCCA -natdeath*HCCA + aging %*% HCCA
         dHCCB <- c1bB*(C1+C1cure) + c2bB*(C2+C2cure) -dthbB*HCCB -tranbB*HCCB -natdeath*HCCB + aging %*% HCCB
         dHCCC <- c1bC*(C1+C1cure) + c2bC*(C2+C2cure) -dthbC*HCCC -natdeath*HCCC + aging %*% HCCC
         dHCCD <- c1bD*(C1+C1cure) + c2bD*(C2+C2cure) + c3bD*(C3+C3cure) + c4bD*(C4+C4cure) -dthbD*HCCD -natdeath*HCCD + aging %*% HCCD
         
         dF0cure <- ifelse(F0cure>=0, scr*pF0scr*sens*cureF0*F0 - natdeath*F0cure + aging %*% F0cure, rep(0,groups))
         dF1cure <- ifelse(F1cure>=0, scr*pF1scr*sens*cureF1*F1 - natdeath*F1cure + aging %*% F1cure, rep(0,groups))
         dF2cure <- ifelse(F2cure>=0, scr*pF2scr*sens*cureF2*F2 - natdeath*F2cure + aging %*% F2cure, rep(0,groups))
         dF3cure <- ifelse(F3cure>=0, scr*pF3scr*sens*cureF3*F3 - natdeath*F3cure + aging %*% F3cure, rep(0,groups))
         dC1cure <- ifelse(C1cure>=0, scr*pC1scr*sens*cureC1*C1 - natdeath*C1cure -(c1bA + c1bB + c1bC + c1bD)*C1cure + aging %*% C1cure, rep(0,groups))
         dC2cure <- ifelse(C2cure>=0, scr*pC2scr*sens*cureC2*C2 - natdeath*C2cure -(c2bA + c2bB + c2bC + c2bD)*C2cure + aging %*% C2cure, rep(0,groups))
         dC3cure <- ifelse(C3cure>=0, scr*pC3scr*sens*cureC3*C3 - natdeath*C3cure - c3bD*C3cure + aging %*% C3cure, rep(0,groups))
         dC4cure <- ifelse(C4cure>=0, scr*pC4scr*sens*cureC4*C4 - natdeath*C4cure - c4bD*C4cure + aging %*% C4cure, rep(0,groups))
         
         dD <- dthc1*C1 + dthc2*C2 + dthc3*C3 + dthc4*C4 + dthbA*HCCA + dthbB*HCCB + dthbC*HCCC + dthbD*HCCD
         ddthC14 <- dthc1*C1 + dthc2*C2 + dthc3*C3 + dthc4*C4
         ddthHCC <- dthbA*HCCA + dthbB*HCCB + dthbC*HCCC + dthbD*HCCD
         
         dtransliv <- tranbA*HCCA + tranbB*HCCB + tranc4*C4
         dCInc <- lambda * S
         
         # Cumulative number of individuals screened/treated
         dCScr <- scr*sens*(pF0scr*cureF0*F0 + 
                              pF1scr*cureF1*F1 + 
                              pF2scr*cureF2*F2 + 
                              pF3scr*cureF3*F3 + 
                              pC1scr*cureC1*C1 +
                              pC2scr*cureC2*C2 + 
                              pC3scr*cureC3*C3 + 
                              pC4scr*cureC4*C4)
         
         newdeath <- ifelse(C1+C2+C3+totalHCC>0, dthc1*C1 + dthc2*C2 + dthc3*C3 + dthc4*C4 + dthbA*HCCA + dthbB*HCCB + dthbC*HCCC + dthbD*HCCD,0)
         incHCC <- c1bA*(C1+C1cure)+c2bA*(C2+C2cure)+c1bB*(C1+C1cure)+
           c2bB*(C2+C2cure)+c1bC*(C1+C1cure)+c2bC*(C2+C2cure)+c1bD*(C1+C1cure)+
           c2bD*(C2+C2cure)+c3bD*(C3+C3cure)+c4bD*(C4+C4cure)
         
         list(c(dS,dF0,dF1,dF2,dF3,dC1,dC2,dC3,dC4,
                dHCCA,dHCCB,dHCCC,dHCCD,
                dF0cure,dF1cure,dF2cure,dF3cure,dC1cure,dC2cure,dC3cure,dC4cure,
                dD,ddthC14,ddthHCC,dtransliv, dCInc, dCScr),
              incidenceHCC=incHCC,newdeath=newdeath,prevalence=prevalence,population=pop,total.infection=infect,
              totalHCC=totalHCC,totalHCV=totalHCV)
       }
  )
}

# Mutate data function ####

mutate_data <- function(ode_output){
  
  df1<-as_tibble(as.data.frame(ode_output)) %>% 
    mutate(
      
      # Totals across all age groups
      S=(S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+S11+S12+S3+S14+S15+S16+S17+S18+S19+S20+S21),
      
      # Total Screened
      CScr=(CScr1+CScr2+CScr3+CScr4+CScr5+CScr6+CScr7+CScr8+CScr9+CScr10+CScr11+CScr12+CScr3+CScr14+CScr15+CScr16+CScr17+CScr18+CScr19+CScr20+CScr21),
      
      # Total Cases
      CInc=(CInc1+CInc2+CInc3+CInc4+CInc5+CInc6+CInc7+CInc8+CInc9+CInc10+CInc11+CInc12+CInc3+CInc14+CInc15+CInc16+CInc17+CInc18+CInc19+CInc20+CInc21),
      
      # Cured
      F0cure=(F0cure1+F0cure2+F0cure3+F0cure4+F0cure5+F0cure6+F0cure7+F0cure8+F0cure9+F0cure10+F0cure11+
                F0cure12+F0cure13+F0cure14+F0cure15+F0cure16+F0cure17+F0cure18+F0cure19+F0cure20+F0cure21),
      F1cure=(F1cure1+F1cure2+F1cure3+F1cure4+F1cure5+F1cure6+F1cure7+F1cure8+F1cure9+F1cure10+F1cure11+
                F1cure12+F1cure13+F1cure14+F1cure15+F1cure16+F1cure17+F1cure18+F1cure19+F1cure20+F1cure21),
      F2cure=(F2cure1+F2cure2+F2cure3+F2cure4+F2cure3+F2cure6+F2cure7+F2cure8+F2cure9+F2cure10+F2cure11+
                F2cure12+F2cure13+F2cure14+F2cure15+F2cure16+F2cure17+F2cure18+F2cure19+F2cure20+F2cure21),
      F3cure=(F3cure1+F3cure2+F3cure3+F3cure4+F3cure5+F3cure6+F3cure7+F3cure8+F3cure9+F3cure10+F3cure11+
                F3cure12+F3cure13+F3cure14+F3cure15+F3cure16+F3cure17+F3cure18+F3cure19+F3cure20+F3cure21),
      C1cure=(C1cure1+C1cure2+C1cure3+C1cure4+C1cure5+C1cure6+C1cure7+C1cure8+C1cure9+C1cure10+C1cure11+
                C1cure12+C1cure13+C1cure14+C1cure15+C1cure16+C1cure17+C1cure18+C1cure19+C1cure20+C1cure21),
      C2cure=(C2cure1+C2cure2+C2cure3+C2cure4+C2cure5+C2cure6+C2cure7+C2cure8+C2cure9+C2cure10+C2cure11+
                C2cure12+C2cure13+C2cure14+C2cure15+C2cure16+C2cure17+C2cure18+C2cure19+C2cure20+C2cure21),
      C3cure=(C3cure1+C3cure2+C3cure3+C3cure4+C3cure5+C3cure6+C3cure7+C3cure8+C3cure9+C3cure10+C3cure11+
                C3cure12+C3cure13+C3cure14+C3cure15+C3cure16+C3cure17+C3cure18+C3cure19+C3cure20+C3cure21),
      C4cure=(C4cure1+C4cure2+C4cure3+C4cure4+C4cure5+C4cure6+C4cure7+C4cure8+C4cure9+C4cure10+C4cure11+
                C4cure12+C4cure13+C4cure14+C4cure15+C4cure16+C4cure17+C4cure18+C4cure19+C4cure20+C4cure21),
      
      cured = F0cure+F1cure+F2cure+F3cure+C1cure+C2cure+C3cure+C4cure,
      
      # Fibrosis, cirrhosis and HCC totals (sum of all age groups)
      F0=(F01+F02+F03+F04+F05+F06+F07+F08+F09+F010+F011+F012+F013+F014+F015+F016+F017+F018+F019+F020+F021),
      F1=(F11+F12+F13+F14+F15+F16+F17+F18+F19+F110+F111+F112+F113+F114+F115+F116+F117+F118+F119+F120+F121),
      F2=(F21+F22+F23+F24+F23+F26+F27+F28+F29+F210+F211+F212+F213+F214+F215+F216+F217+F218+F219+F220+F221),
      F3=(F31+F32+F33+F34+F35+F36+F37+F38+F39+F310+F311+F312+F313+F314+F315+F316+F317+F318+F319+F320+F321),
      C1=(C11+C12+C13+C14+C15+C16+C17+C18+C19+C110+C111+C112+C113+C114+C115+C116+C117+C118+C119+C120+C121),
      C2=(C21+C22+C23+C24+C25+C26+C27+C28+C29+C210+C211+C212+C213+C214+C215+C216+C217+C218+C219+C220+C221),
      C3=(C31+C32+C33+C34+C35+C36+C37+C38+C39+C310+C311+C312+C313+C314+C315+C316+C317+C318+C319+C320+C321),
      C4=(C41+C42+C43+C44+C45+C46+C47+C48+C49+C410+C411+C412+C413+C414+C415+C416+C417+C418+C419+C420+C421),
      HCCA=(HCCA1+HCCA2+HCCA3+HCCA4+HCCA5+HCCA6+HCCA7+HCCA8+HCCA9+HCCA10+HCCA11+HCCA12+HCCA13+HCCA14+HCCA15+HCCA16+HCCA17+HCCA18+HCCA19+HCCA20+HCCA21),
      HCCB=(HCCB1+HCCB2+HCCB3+HCCB4+HCCB5+HCCB6+HCCB7+HCCB8+HCCB9+HCCB10+HCCB11+HCCB12+HCCB13+HCCB14+HCCB15+HCCB16+HCCB17+HCCB18+HCCB19+HCCB20+HCCB21),
      HCCC=(HCCC1+HCCC2+HCCC3+HCCC4+HCCC5+HCCC6+HCCC7+HCCC8+HCCC9+HCCC10+HCCC11+HCCC12+HCCC13+HCCC14+HCCC15+HCCC16+HCCC17+HCCC18+HCCC19+HCCC20+HCCC21),
      HCCD=(HCCD1+HCCD2+HCCD3+HCCD4+HCCD5+HCCD6+HCCD7+HCCD8+HCCD9+HCCD10+HCCD11+HCCD12+HCCD13+HCCD14+HCCD15+HCCD16+HCCD17+HCCD18+HCCD19+HCCD20+HCCD21),
      
      # Incidence and mortality
      CInc=(CInc1+CInc2+CInc3+CInc4+CInc5+CInc6+CInc7+CInc8+CInc9+CInc10+CInc11+CInc12+CInc13+CInc14+CInc15+CInc16+CInc17+CInc18+CInc19+CInc20+CInc21),
      D=D1+D2+D3+D4+D5+D6+D7+D8+D9+D10+D11+D12+D13+D14+D15+D16+D17+D18+D19+D20+D21,
      # Yearly mortality
      Deaths = c(0, diff(D)),
      # Yearly incidence
      Inc = c(0, diff(CInc)),
      
      Inc5 = c(0, diff(CInc5)),
      Inc6 = c(0, diff(CInc6)),
      Inc7 = c(0, diff(CInc7)),
      Inc8 = c(0, diff(CInc8)),
      Inc9 = c(0, diff(CInc9)),
      Inc10 = c(0, diff(CInc10)),
      Inc11 = c(0, diff(CInc11)),
      Inc12 = c(0, diff(CInc12)),
      Inc13 = c(0, diff(CInc13)),
      Inc14 = c(0, diff(CInc14)),
      
      # Total infections, HCC and HCV (sum across all age groups)
      infect = (F0+F1+F2+F3+C1+C2+C3+C4+HCCA+HCCB+HCCC+HCCD),
      totalHCC = (HCCA+HCCB+HCCC+HCCD),
      totalHCV  = (F0+F1+F2+F3+C1+C2+C3+C4),
      
      # Total population of entire system
      total = population,
      
      # Prevalance of the above (%)
      infectprev = 100*(infect/total),
      HCCprev = 100*(totalHCC/total),
      HCVprev = 100*(totalHCV/total),
      
      # Total infections for each individual age group
      infect1 = (F01+F11+F21+F31+C11+C21+C31+C41+HCCA1+HCCB1+HCCC1+HCCD1),
      infect2 = (F02+F12+F22+F32+C12+C22+C32+C42+HCCA2+HCCB2+HCCC2+HCCD2),
      infect3 = (F03+F13+F23+F33+C13+C23+C33+C43+HCCA3+HCCB3+HCCC3+HCCD3),
      infect4 = (F04+F14+F24+F34+C14+C24+C34+C44+HCCA4+HCCB4+HCCC4+HCCD4),
      infect5 = (F05+F15+F25+F35+C15+C25+C35+C45+HCCA5+HCCB5+HCCC5+HCCD5),
      infect6 = (F06+F16+F26+F36+C16+C26+C36+C46+HCCA6+HCCB6+HCCC6+HCCD6),
      infect7 = (F07+F17+F27+F37+C17+C27+C37+C47+HCCA7+HCCB7+HCCC7+HCCD7),
      infect8 = (F08+F18+F28+F38+C18+C28+C38+C48+HCCA8+HCCB8+HCCC8+HCCD8),
      infect9 = (F09+F19+F29+F39+C19+C29+C39+C49+HCCA9+HCCB9+HCCC9+HCCD9),
      infect10 = (F010+F110+F210+F310+C110+C210+C310+C410+HCCA10+HCCB10+HCCC10+HCCD10),
      infect11 = (F011+F111+F211+F311+C111+C211+C311+C411+HCCA11+HCCB11+HCCC11+HCCD11),
      infect12 = (F012+F112+F212+F312+C112+C212+C312+C412+HCCA12+HCCB12+HCCC12+HCCD12),
      infect13 = (F013+F113+F213+F313+C113+C213+C313+C413+HCCA13+HCCB13+HCCC13+HCCD13),
      infect14 = (F014+F114+F214+F314+C114+C214+C314+C414+HCCA14+HCCB14+HCCC14+HCCD14),
      infect15 = (F015+F115+F215+F315+C115+C215+C315+C415+HCCA15+HCCB15+HCCC15+HCCD15),
      infect16 = (F016+F116+F216+F316+C116+C216+C316+C416+HCCA16+HCCB16+HCCC16+HCCD16),
      infect17 = (F017+F117+F217+F317+C117+C217+C317+C417+HCCA17+HCCB17+HCCC17+HCCD17),
      infect18 = (F018+F118+F218+F318+C118+C218+C318+C418+HCCA18+HCCB18+HCCC18+HCCD18),
      infect19 = (F019+F119+F219+F319+C119+C219+C319+C419+HCCA19+HCCB19+HCCC19+HCCD19),
      infect20 = (F020+F120+F220+F320+C120+C220+C320+C420+HCCA20+HCCB20+HCCC20+HCCD20),
      infect21 = (F021+F121+F221+F321+C121+C221+C321+C421+HCCA21+HCCB21+HCCC21+HCCD21),
      
      # Total proportion in each age group (sum of all compartments per age group)
      group1 = (S1+F01+F11+F21+F31+C11+C21+C31+C41+HCCA1+HCCB1+HCCC1+HCCD1+F0cure1+F1cure1+F2cure1+F3cure1+C1cure1+C2cure1+C3cure1+C4cure1), # which compartments contribute? not the cumulative ones but do I need std and new cured or will this result in overcounting?
      group2 = (S2+F02+F12+F22+F32+C12+C22+C32+C42+HCCA2+HCCB2+HCCC2+HCCD2+F0cure2+F1cure2+F2cure2+F3cure2+C1cure2+C2cure2+C3cure2+C4cure2),
      group3 = (S3+F03+F13+F23+F33+C13+C23+C33+C43+HCCA3+HCCB3+HCCC3+HCCD3+F0cure3+F1cure3+F2cure3+F3cure3+C1cure3+C2cure3+C3cure3+C4cure3),
      group4 = (S4+F04+F14+F24+F34+C14+C24+C34+C44+HCCA4+HCCB4+HCCC4+HCCD4+F0cure4+F1cure4+F2cure4+F3cure4+C1cure4+C2cure4+C3cure4+C4cure4),
      group5 = (S5+F05+F15+F25+F35+C15+C25+C35+C45+HCCA5+HCCB5+HCCC5+HCCD5+F0cure5+F1cure5+F2cure5+F3cure5+C1cure5+C2cure5+C3cure5+C4cure5),
      group6 = (S6+F06+F16+F26+F36+C16+C26+C36+C46+HCCA6+HCCB6+HCCC6+HCCD6+F0cure6+F1cure6+F2cure6+F3cure6+C1cure6+C2cure6+C3cure6+C4cure6),
      group7 = (S7+F07+F17+F27+F37+C17+C27+C37+C47+HCCA7+HCCB7+HCCC7+HCCD7+F0cure7+F1cure7+F2cure7+F3cure7+C1cure7+C2cure7+C3cure7+C4cure7),
      group8 = (S8+F08+F18+F28+F38+C18+C28+C38+C48+HCCA8+HCCB8+HCCC8+HCCD8+F0cure8+F1cure8+F2cure8+F3cure8+C1cure8+C2cure8+C3cure8+C4cure8),
      group9 = (S9+F09+F19+F29+F39+C19+C29+C39+C49+HCCA9+HCCB9+HCCC9+HCCD9+F0cure9+F1cure9+F2cure9+F3cure9+C1cure9+C2cure9+C3cure9+C4cure9),
      group10 = (S10+F010+F110+F210+F310+C110+C210+C310+C410+HCCA10+HCCB10+HCCC10+HCCD10+F0cure10+F1cure10+F2cure10+F3cure10+C1cure10+C2cure10+C3cure10+C4cure10), 
      group11 = (S11+F011+F111+F211+F311+C111+C211+C311+C411+HCCA11+HCCB11+HCCC11+HCCD11+F0cure11+F1cure11+F2cure11+F3cure11+C1cure11+C2cure11+C3cure11+C4cure11),
      group12 = (S12+F012+F112+F212+F312+C112+C212+C312+C412+HCCA12+HCCB12+HCCC12+HCCD12+F0cure12+F1cure12+F2cure12+F3cure12+C1cure12+C2cure12+C3cure12+C4cure12),
      group13 = (S13+F013+F113+F213+F313+C113+C213+C313+C413+HCCA13+HCCB13+HCCC13+HCCD13+F0cure13+F1cure13+F2cure13+F3cure13+C1cure13+C2cure13+C3cure13+C4cure13),
      group14 = (S14+F014+F114+F214+F314+C114+C214+C314+C414+HCCA14+HCCB14+HCCC14+HCCD14+F0cure14+F1cure14+F2cure14+F3cure14+C1cure14+C2cure14+C3cure14+C4cure14),
      group15 = (S15+F015+F115+F215+F315+C115+C215+C315+C415+HCCA15+HCCB15+HCCC15+HCCD15+F0cure15+F1cure15+F2cure15+F3cure15+C1cure15+C2cure15+C3cure15+C4cure15),
      group16 = (S16+F016+F116+F216+F316+C116+C216+C316+C416+HCCA16+HCCB16+HCCC16+HCCD16+F0cure16+F1cure16+F2cure16+F3cure16+C1cure16+C2cure16+C3cure16+C4cure16),
      group17 = (S17+F017+F117+F217+F317+C117+C217+C317+C417+HCCA17+HCCB17+HCCC17+HCCD17+F0cure17+F1cure17+F2cure17+F3cure17+C1cure17+C2cure17+C3cure17+C4cure17),
      group18 = (S18+F018+F118+F218+F318+C118+C218+C318+C418+HCCA18+HCCB18+HCCC18+HCCD18+F0cure18+F1cure18+F2cure18+F3cure18+C1cure18+C2cure18+C3cure18+C4cure18),
      group19 = (S19+F019+F119+F219+F319+C119+C219+C319+C419+HCCA19+HCCB19+HCCC19+HCCD19+F0cure19+F1cure19+F2cure19+F3cure19+C1cure19+C2cure19+C3cure19+C4cure19),
      group20 = (S20+F020+F120+F220+F320+C120+C220+C320+C420+HCCA20+HCCB20+HCCC20+HCCD20+F0cure20+F1cure20+F2cure20+F3cure20+C1cure20+C2cure20+C3cure20+C4cure20),
      group21 = (S21+F021+F121+F221+F321+C121+C221+C321+C421+HCCA21+HCCB21+HCCC21+HCCD21+F0cure21+F1cure21+F2cure21+F3cure21+C1cure21+C2cure21+C3cure21+C4cure21),
      
      # Prevalence by age group (%)
      prev1 = (infect1 * 100) / (group1),
      prev2 = (infect2 * 100) / (group2),
      prev3 = (infect3 * 100) / (group3),
      prev4 = (infect4 * 100) / (group4),
      prev5 = (infect5 * 100) / (group5),
      prev6 = (infect6 * 100) / (group6),
      prev7 = (infect7 * 100) / (group7),
      prev8 = (infect8 * 100) / (group8),
      prev9 = (infect9 * 100) / (group9),
      prev10 = (infect10 * 100) / (group10),
      prev11 = (infect11 * 100) / (group11),
      prev12 = (infect12 * 100) / (group12),
      prev13 = (infect13 * 100) / (group13),
      prev14 = (infect14 * 100) / (group14),
      prev15 = (infect15 * 100) / (group15),
      prev16 = (infect16 * 100) / (group16),
      prev17 = (infect17 * 100) / (group17),
      prev18 = (infect18 * 100) / (group18),
      prev19 = (infect19 * 100) / (group19),
      prev20 = (infect20 * 100) / (group20),
      prev21 = (infect21 * 100) / (group21),
      
      # Aggregated prevalence by age group from puclished data
      prev1_2 = 100*(infect1+infect2)/(group1+group2), # 0-9
      prev3_4 = 100*(infect3+infect4)/(group3+group4), # 10-19
      prev5_6 = 100*(infect5+infect6)/(group5+group6), # 20-29
      prev7_8 = 100*(infect7+infect8)/(group7+group8), # 30-39
      prev9_10 = 100*(infect9+infect10)/(group9+group10), # 40-49
      prev11_21 = 100*(infect11+infect12+infect13+infect14+infect15+infect16+infect17+infect18+infect19+infect20+infect21)/(group11+group12+group13+group14+group15+group16+group17+group18+group19+group20+group21) # Over 50
      
    ) %>% 
    pivot_longer(names_to = "variable", cols = !1) 
  return(df1)
}


# Collect results function ####

Collect_Results <- function(dataframe){
  
  inc_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("Inc")))$value
  inc_2015 <- as.numeric((dataframe %>% filter(variable %in% c("Inc")) %>% filter(Year %in% c("2015")))[,3])
  inc_target <- 0.1*inc_2015
  inc_target_diff <- inc_2030 - inc_target
  
  mort_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("Deaths")))$value
  mort_2015 <- as.numeric((dataframe %>% filter(variable %in% c("Deaths")) %>% filter(Year %in% c("2015")))[,3])
  mort_target <- 0.35*mort_2015
  mort_target_diff <- mort_2030 - mort_target
  
  CScr_2023 <- (dataframe %>% filter(Year %in% c("2023")) %>% filter(variable %in% c("CScr")))$value
  CScr_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("CScr")))$value
  total_screened <- CScr_2030 - CScr_2023
  
  CInc_2023 <- (dataframe %>% filter(Year %in% c("2023")) %>% filter(variable %in% c("CInc")))$value
  CInc_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("CInc")))$value
  total_cases <- CInc_2030 - CInc_2023
  
  CD_2023 <- (dataframe %>% filter(Year %in% c("2023")) %>% filter(variable %in% c("D")))$value
  CD_2030 <- (dataframe %>% filter(Year %in% c("2030")) %>% filter(variable %in% c("D")))$value
  total_deaths <- CD_2030 - CD_2023
  
  inc_elim_year <- 0
  inc_only <- dataframe %>% filter(variable %in% c("Inc"))
  for(i in 2:nrow(inc_only)){
    if(inc_only[i,]$value <= inc_target){
      inc_elim_year <- inc_only[i,]$Year
      break
    }
  }
  
  if(inc_elim_year==2040 | inc_elim_year==0){
    inc_elim_year <- "Beyond simulation" }
  
  mort_elim_year <- 0
  mort_only <- dataframe %>% filter(variable %in% c("Deaths"))
  for(i in 2:nrow(mort_only)){
    if(mort_only[i,]$value <= mort_target){
      mort_elim_year <- mort_only[i,]$Year
      break
    }
  }
  
  if(mort_elim_year==0){
    mort_elim_year <- "Beyond simulation"
  }
  
  
  
  if(inc_target_diff<=0){
    Inc_Elim = "Yes"}
  else{Inc_Elim = "No"}
  
  if(mort_target_diff<=0){
    Mort_Elim = "Yes"}
  else{Mort_Elim = "No"}
  
  
  return(c("Incidence 2030" = as.numeric(inc_2030), "Incidence Difference to Target" = as.numeric(inc_target_diff), "Deaths 2030" = as.numeric(mort_2030), "Deaths Difference to Target" = as.numeric(mort_target_diff), "Total Cases" = as.numeric(total_cases), "Cases Averted" = NA, "Total Deaths" = total_deaths, "Deaths Averted" = NA, "Total Screened" =  as.numeric(total_screened), "Extra Screened" = NA, "Year Incidence Elimination Reached" = inc_elim_year, "Year Deaths Elimination Reached" = mort_elim_year))
  
}

# Manipulate dataframe functions ####

annotate_dataframes <- function(df, population_scenario, baseline_coverage, brm, scr_group, scr_cov){
  
  df$population_scenario <- rep(population_scenario, nrow(df))
  df$baseline_coverage <- rep(baseline_coverage, nrow(df))
  df$brm <- rep(brm, nrow(df))
  df$scr_group <- rep(scr_group, nrow(df))
  df$scr_cov <- rep(scr_cov, nrow(df))
  
  return(df)
}

get_details <- function(x){
  
  dataframe <- as.character(x[1])
  population_scenario <- as.character(x[2])
  baseline_coverage <- as.character(x[3])
  brm <- as.character(x[4])
  scr_group <- as.character(x[5])
  scr_cov <- as.character(x[6])
  
  return(list(dataframe, population_scenario, baseline_coverage, brm, scr_group, scr_cov))
  
}

modify_dataframe <- function(df, char){
  
  line <- simulations[which(simulations==char),]
  details <- get_details(line)
  new_dataframe <- annotate_dataframes(df, details[[2]], details[[3]], details[[4]], details[[5]], details[[6]])
  
  return(new_dataframe)
}

combine_dataframes <- function(...){
  
  new <- list()
  for(i in 1:length(...)){
    char <- names(...)[i]
    new[[i]] <- modify_dataframe(...[[i]], char)
  }
  new_dataframe <- 
    return(new)
}

# Run all simulations ####

# Baseline population scenario

mort_scenario <- scenario_base
mortality_base <- mortality_base_empty

for(i in 19:37){
  mortality_base[i,2:22] <- mortality_base[i-1,2:22]*mort_scenario
}

mortality_approx <- mortality_base

mortality_func1 <- approxfun(mortality_approx$t,mortality_approx$rate1,method="linear")
mortality_func2 <- approxfun(mortality_approx$t,mortality_approx$rate2,method="linear")
mortality_func3 <- approxfun(mortality_approx$t,mortality_approx$rate3,method="linear")
mortality_func4 <- approxfun(mortality_approx$t,mortality_approx$rate4,method="linear")
mortality_func5 <- approxfun(mortality_approx$t,mortality_approx$rate5,method="linear")
mortality_func6 <- approxfun(mortality_approx$t,mortality_approx$rate6,method="linear")
mortality_func7 <- approxfun(mortality_approx$t,mortality_approx$rate7,method="linear")
mortality_func8 <- approxfun(mortality_approx$t,mortality_approx$rate8,method="linear")
mortality_func9 <- approxfun(mortality_approx$t,mortality_approx$rate9,method="linear")
mortality_func10 <- approxfun(mortality_approx$t,mortality_approx$rate10,method="linear")
mortality_func11 <- approxfun(mortality_approx$t,mortality_approx$rate11,method="linear")
mortality_func12 <- approxfun(mortality_approx$t,mortality_approx$rate12,method="linear")
mortality_func13 <- approxfun(mortality_approx$t,mortality_approx$rate13,method="linear")
mortality_func14 <- approxfun(mortality_approx$t,mortality_approx$rate14,method="linear")
mortality_func15 <- approxfun(mortality_approx$t,mortality_approx$rate15,method="linear")
mortality_func16 <- approxfun(mortality_approx$t,mortality_approx$rate16,method="linear")
mortality_func17 <- approxfun(mortality_approx$t,mortality_approx$rate17,method="linear")
mortality_func18 <- approxfun(mortality_approx$t,mortality_approx$rate18,method="linear")
mortality_func19 <- approxfun(mortality_approx$t,mortality_approx$rate19,method="linear")
mortality_func20 <- approxfun(mortality_approx$t,mortality_approx$rate20,method="linear")
mortality_func21 <- approxfun(mortality_approx$t,mortality_approx$rate21,method="linear")

birth.approx_mean <- birth.approx
birthrate_multiplier <- 1.09
birth.approx_mean$birth <- birthrate_multiplier * birth.approx_mean$birth
birth.func <- approxfun(birth.approx_mean$t,birth.approx_mean$birth,method="linear")

groupBASELINE <- groupBASE_lower
coverageBASE <- cov_lower

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE_lower
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df1 <- mutate_data(out)
results1 <- Collect_Results(df1)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df2 <- mutate_data(out)
results2 <- Collect_Results(df2)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df3 <- mutate_data(out)
results3 <- Collect_Results(df3)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df4 <- mutate_data(out)
results4 <- Collect_Results(df4)

df5 <- df1
results5 <- results1

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df6 <- mutate_data(out)
results6 <- Collect_Results(df6)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df7 <- mutate_data(out)
results7 <- Collect_Results(df7)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df8 <- mutate_data(out)
results8 <- Collect_Results(df8)

df9 <- df1
results9 <- results1

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df10 <- mutate_data(out)
results10 <- Collect_Results(df10)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df11 <- mutate_data(out)
results11 <- Collect_Results(df11)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df12 <- mutate_data(out)
results12 <- Collect_Results(df12)

df13 <- df1
results13 <- results1

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df14 <- mutate_data(out)
results14 <- Collect_Results(df14)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df15 <- mutate_data(out)
results15 <- Collect_Results(df15)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df16 <- mutate_data(out)
results16 <- Collect_Results(df16)

groupBASELINE <- groupBASE
coverageBASE <- cov_mean

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df17 <- mutate_data(out)
results17 <- Collect_Results(df17)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df18 <- mutate_data(out)
results18 <- Collect_Results(df18)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df19 <- mutate_data(out)
results19 <- Collect_Results(df19)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df20 <- mutate_data(out)
results20 <- Collect_Results(df20)

df21 <- df17
results21 <- results17

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df22 <- mutate_data(out)
results22 <- Collect_Results(df22)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df23 <- mutate_data(out)
results23 <- Collect_Results(df23)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df24 <- mutate_data(out)
results24 <- Collect_Results(df24)

df25 <- df17
results25 <- results17

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df26 <- mutate_data(out)
results26 <- Collect_Results(df26)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df27 <- mutate_data(out)
results27 <- Collect_Results(df27)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df28 <- mutate_data(out)
results28 <- Collect_Results(df28)

df29 <- df17
results29 <- results17

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df30 <- mutate_data(out)
results30 <- Collect_Results(df30)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df31 <- mutate_data(out)
results31 <- Collect_Results(df31)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df32 <- mutate_data(out)
results32 <- Collect_Results(df32)

groupBASELINE <- groupBASE_upper
coverageBASE <- cov_upper

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df33 <- mutate_data(out)
results33 <- Collect_Results(df33)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df34 <- mutate_data(out)
results34 <- Collect_Results(df34)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df35 <- mutate_data(out)
results35 <- Collect_Results(df35)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df36 <- mutate_data(out)
results36 <- Collect_Results(df36)

df37 <- df33
results37 <- results33

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df38 <- mutate_data(out)
results38 <- Collect_Results(df38)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df39 <- mutate_data(out)
results39 <- Collect_Results(df39)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df40 <- mutate_data(out)
results40 <- Collect_Results(df40)

df41 <- df33
results41 <- results33

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df42 <- mutate_data(out)
results42 <- Collect_Results(df42)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df43 <- mutate_data(out)
results43 <- Collect_Results(df43)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df44 <- mutate_data(out)
results44 <- Collect_Results(df44)

df45 <- df33
results45 <- results33

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df46 <- mutate_data(out)
results46 <- Collect_Results(df46)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df47 <- mutate_data(out)
results47 <- Collect_Results(df47)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df48 <- mutate_data(out)
results48 <- Collect_Results(df48)

# Decline population scenario

mort_scenario <- scenario_1
mortality_base <- mortality_base_empty

for(i in 19:37){
  mortality_base[i,2:22] <- mortality_base[i-1,2:22]*mort_scenario
}

mortality_approx <- mortality_base

mortality_func1 <- approxfun(mortality_approx$t,mortality_approx$rate1,method="linear")
mortality_func2 <- approxfun(mortality_approx$t,mortality_approx$rate2,method="linear")
mortality_func3 <- approxfun(mortality_approx$t,mortality_approx$rate3,method="linear")
mortality_func4 <- approxfun(mortality_approx$t,mortality_approx$rate4,method="linear")
mortality_func5 <- approxfun(mortality_approx$t,mortality_approx$rate5,method="linear")
mortality_func6 <- approxfun(mortality_approx$t,mortality_approx$rate6,method="linear")
mortality_func7 <- approxfun(mortality_approx$t,mortality_approx$rate7,method="linear")
mortality_func8 <- approxfun(mortality_approx$t,mortality_approx$rate8,method="linear")
mortality_func9 <- approxfun(mortality_approx$t,mortality_approx$rate9,method="linear")
mortality_func10 <- approxfun(mortality_approx$t,mortality_approx$rate10,method="linear")
mortality_func11 <- approxfun(mortality_approx$t,mortality_approx$rate11,method="linear")
mortality_func12 <- approxfun(mortality_approx$t,mortality_approx$rate12,method="linear")
mortality_func13 <- approxfun(mortality_approx$t,mortality_approx$rate13,method="linear")
mortality_func14 <- approxfun(mortality_approx$t,mortality_approx$rate14,method="linear")
mortality_func15 <- approxfun(mortality_approx$t,mortality_approx$rate15,method="linear")
mortality_func16 <- approxfun(mortality_approx$t,mortality_approx$rate16,method="linear")
mortality_func17 <- approxfun(mortality_approx$t,mortality_approx$rate17,method="linear")
mortality_func18 <- approxfun(mortality_approx$t,mortality_approx$rate18,method="linear")
mortality_func19 <- approxfun(mortality_approx$t,mortality_approx$rate19,method="linear")
mortality_func20 <- approxfun(mortality_approx$t,mortality_approx$rate20,method="linear")
mortality_func21 <- approxfun(mortality_approx$t,mortality_approx$rate21,method="linear")

birth.approx_mean <- birth.approx
birthrate_multiplier <- 0.95
birth.approx_mean$birth <- birthrate_multiplier * birth.approx_mean$birth
birth.func <- approxfun(birth.approx_mean$t,birth.approx_mean$birth,method="linear")

groupBASELINE <- groupBASE_lower
coverageBASE <- cov_lower

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE_lower
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df49 <- mutate_data(out)
results49 <- Collect_Results(df49)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df50 <- mutate_data(out)
results50 <- Collect_Results(df50)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df51 <- mutate_data(out)
results51 <- Collect_Results(df51)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df52 <- mutate_data(out)
results52 <- Collect_Results(df52)

df53 <- df49
results53 <- results49

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df54 <- mutate_data(out)
results54 <- Collect_Results(df54)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df55 <- mutate_data(out)
results55 <- Collect_Results(df55)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df56 <- mutate_data(out)
results56 <- Collect_Results(df56)

df57 <- df49
results57 <- results49

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df58 <- mutate_data(out)
results58 <- Collect_Results(df58)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df59 <- mutate_data(out)
results59 <- Collect_Results(df59)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df60 <- mutate_data(out)
results60 <- Collect_Results(df60)

df61 <- df49
results61 <- results49

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df62 <- mutate_data(out)
results62 <- Collect_Results(df62)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df63 <- mutate_data(out)
results63 <- Collect_Results(df63)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df64 <- mutate_data(out)
results64 <- Collect_Results(df64)

groupBASELINE <- groupBASE
coverageBASE <- cov_mean

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df65 <- mutate_data(out)
results65 <- Collect_Results(df65)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df66 <- mutate_data(out)
results66 <- Collect_Results(df66)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df67 <- mutate_data(out)
results67 <- Collect_Results(df67)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df68 <- mutate_data(out)
results68 <- Collect_Results(df68)

df69 <- df65
results69 <- results65

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df70 <- mutate_data(out)
results70 <- Collect_Results(df70)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df71 <- mutate_data(out)
results71 <- Collect_Results(df71)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df72 <- mutate_data(out)
results72 <- Collect_Results(df72)

df73 <- df65
results73 <- results65

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df74 <- mutate_data(out)
results74 <- Collect_Results(df74)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df75 <- mutate_data(out)
results75 <- Collect_Results(df75)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df76 <- mutate_data(out)
results76 <- Collect_Results(df76)

df77 <- df65
results77 <- results65

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df78 <- mutate_data(out)
results78 <- Collect_Results(df78)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df79 <- mutate_data(out)
results79 <- Collect_Results(df79)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df80 <- mutate_data(out)
results80 <- Collect_Results(df80)

groupBASELINE <- groupBASE_upper
coverageBASE <- cov_upper

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE_upper
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df81 <- mutate_data(out)
results81 <- Collect_Results(df81)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df82 <- mutate_data(out)
results82 <- Collect_Results(df82)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df83 <- mutate_data(out)
results83 <- Collect_Results(df83)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df84 <- mutate_data(out)
results84 <- Collect_Results(df84)

df85 <- df81
results85 <- results81

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df86 <- mutate_data(out)
results86 <- Collect_Results(df86)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df87 <- mutate_data(out)
results87 <- Collect_Results(df87)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df88 <- mutate_data(out)
results88 <- Collect_Results(df88)

df89 <- df81
results89 <- results81

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df90 <- mutate_data(out)
results90 <- Collect_Results(df90)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df91 <- mutate_data(out)
results91 <- Collect_Results(df91)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df92 <- mutate_data(out)
results92 <- Collect_Results(df92)

df93 <- df81
results93 <- results81

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df94 <- mutate_data(out)
results94 <- Collect_Results(df94)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df95 <- mutate_data(out)
results95 <- Collect_Results(df95)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df96 <- mutate_data(out)
results96 <- Collect_Results(df96)

# Growth population scenario

mort_scenario <- scenario_2
mortality_base <- mortality_base_empty

for(i in 19:37){
  mortality_base[i,2:22] <- mortality_base[i-1,2:22]*mort_scenario
}

mortality_approx <- mortality_base

mortality_func1 <- approxfun(mortality_approx$t,mortality_approx$rate1,method="linear")
mortality_func2 <- approxfun(mortality_approx$t,mortality_approx$rate2,method="linear")
mortality_func3 <- approxfun(mortality_approx$t,mortality_approx$rate3,method="linear")
mortality_func4 <- approxfun(mortality_approx$t,mortality_approx$rate4,method="linear")
mortality_func5 <- approxfun(mortality_approx$t,mortality_approx$rate5,method="linear")
mortality_func6 <- approxfun(mortality_approx$t,mortality_approx$rate6,method="linear")
mortality_func7 <- approxfun(mortality_approx$t,mortality_approx$rate7,method="linear")
mortality_func8 <- approxfun(mortality_approx$t,mortality_approx$rate8,method="linear")
mortality_func9 <- approxfun(mortality_approx$t,mortality_approx$rate9,method="linear")
mortality_func10 <- approxfun(mortality_approx$t,mortality_approx$rate10,method="linear")
mortality_func11 <- approxfun(mortality_approx$t,mortality_approx$rate11,method="linear")
mortality_func12 <- approxfun(mortality_approx$t,mortality_approx$rate12,method="linear")
mortality_func13 <- approxfun(mortality_approx$t,mortality_approx$rate13,method="linear")
mortality_func14 <- approxfun(mortality_approx$t,mortality_approx$rate14,method="linear")
mortality_func15 <- approxfun(mortality_approx$t,mortality_approx$rate15,method="linear")
mortality_func16 <- approxfun(mortality_approx$t,mortality_approx$rate16,method="linear")
mortality_func17 <- approxfun(mortality_approx$t,mortality_approx$rate17,method="linear")
mortality_func18 <- approxfun(mortality_approx$t,mortality_approx$rate18,method="linear")
mortality_func19 <- approxfun(mortality_approx$t,mortality_approx$rate19,method="linear")
mortality_func20 <- approxfun(mortality_approx$t,mortality_approx$rate20,method="linear")
mortality_func21 <- approxfun(mortality_approx$t,mortality_approx$rate21,method="linear")

birth.approx_mean <- birth.approx
birthrate_multiplier <- 1.11
birth.approx_mean$birth <- birthrate_multiplier * birth.approx_mean$birth
birth.func <- approxfun(birth.approx_mean$t,birth.approx_mean$birth,method="linear")

groupBASELINE <- groupBASE_lower
coverageBASE <- cov_lower

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE_lower
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df97 <- mutate_data(out)
results97 <- Collect_Results(df97)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df98 <- mutate_data(out)
results98 <- Collect_Results(df98)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df99 <- mutate_data(out)
results99 <- Collect_Results(df99)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df100 <- mutate_data(out)
results100 <- Collect_Results(df100)

df101 <- df97
results101 <- results97

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df102 <- mutate_data(out)
results102 <- Collect_Results(df102)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df103 <- mutate_data(out)
results103 <- Collect_Results(df103)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df104 <- mutate_data(out)
results104 <- Collect_Results(df104)

df105 <- df97
results105 <- results97

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df106 <- mutate_data(out)
results106 <- Collect_Results(df106)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df107 <- mutate_data(out)
results107 <- Collect_Results(df107)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df108 <- mutate_data(out)
results108 <- Collect_Results(df108)

df109 <- df97
results109 <- results97

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df110 <- mutate_data(out)
results110 <- Collect_Results(df110)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df111 <- mutate_data(out)
results111 <- Collect_Results(df111)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df112 <- mutate_data(out)
results112 <- Collect_Results(df112)

groupBASELINE <- groupBASE
coverageBASE <- cov_mean

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df113 <- mutate_data(out)
results113 <- Collect_Results(df113)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df114 <- mutate_data(out)
results114 <- Collect_Results(df114)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df115 <- mutate_data(out)
results115 <- Collect_Results(df115)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df116 <- mutate_data(out)
results116 <- Collect_Results(df116)

df117 <- df113
results117 <- results113

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df118 <- mutate_data(out)
results118 <- Collect_Results(df118)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df119 <- mutate_data(out)
results119 <- Collect_Results(df119)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df120 <- mutate_data(out)
results120 <- Collect_Results(df120)

df121 <- df113
results121 <- results113

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df122 <- mutate_data(out)
results122 <- Collect_Results(df122)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df123 <- mutate_data(out)
results123 <- Collect_Results(df123)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df124 <- mutate_data(out)
results124 <- Collect_Results(df124)

df125 <- df113
results125 <- results113

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df126 <- mutate_data(out)
results126 <- Collect_Results(df126)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df127 <- mutate_data(out)
results127 <- Collect_Results(df127)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df128 <- mutate_data(out)
results128 <- Collect_Results(df128)

groupBASELINE <- groupBASE_upper
coverageBASE <- cov_upper

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE_upper
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df129 <- mutate_data(out)
results129 <- Collect_Results(df129)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df130 <- mutate_data(out)
results130 <- Collect_Results(df130)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df131 <- mutate_data(out)
results131 <- Collect_Results(df131)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df132 <- mutate_data(out)
results132 <- Collect_Results(df132)

df133 <- df129
results133 <- results129

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df134 <- mutate_data(out)
results134 <- Collect_Results(df134)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df135 <- mutate_data(out)
results135 <- Collect_Results(df135)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df136 <- mutate_data(out)
results136 <- Collect_Results(df136)

df137 <- df129
results137 <- results129

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df138 <- mutate_data(out)
results138 <- Collect_Results(df138)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df139 <- mutate_data(out)
results139 <- Collect_Results(df139)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df140 <- mutate_data(out)
results140 <- Collect_Results(df140)

df141 <- df129
results141 <- results129

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df142 <- mutate_data(out)
results142 <- Collect_Results(df142)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df143 <- mutate_data(out)
results143 <- Collect_Results(df143)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df144 <- mutate_data(out)
results144 <- Collect_Results(df144)

# Plateau population scenario

mort_scenario <- scenario_3
mortality_base <- mortality_base_empty

for(i in 19:37){
  mortality_base[i,2:22] <- mortality_base[i-1,2:22]*mort_scenario
}

mortality_approx <- mortality_base

mortality_func1 <- approxfun(mortality_approx$t,mortality_approx$rate1,method="linear")
mortality_func2 <- approxfun(mortality_approx$t,mortality_approx$rate2,method="linear")
mortality_func3 <- approxfun(mortality_approx$t,mortality_approx$rate3,method="linear")
mortality_func4 <- approxfun(mortality_approx$t,mortality_approx$rate4,method="linear")
mortality_func5 <- approxfun(mortality_approx$t,mortality_approx$rate5,method="linear")
mortality_func6 <- approxfun(mortality_approx$t,mortality_approx$rate6,method="linear")
mortality_func7 <- approxfun(mortality_approx$t,mortality_approx$rate7,method="linear")
mortality_func8 <- approxfun(mortality_approx$t,mortality_approx$rate8,method="linear")
mortality_func9 <- approxfun(mortality_approx$t,mortality_approx$rate9,method="linear")
mortality_func10 <- approxfun(mortality_approx$t,mortality_approx$rate10,method="linear")
mortality_func11 <- approxfun(mortality_approx$t,mortality_approx$rate11,method="linear")
mortality_func12 <- approxfun(mortality_approx$t,mortality_approx$rate12,method="linear")
mortality_func13 <- approxfun(mortality_approx$t,mortality_approx$rate13,method="linear")
mortality_func14 <- approxfun(mortality_approx$t,mortality_approx$rate14,method="linear")
mortality_func15 <- approxfun(mortality_approx$t,mortality_approx$rate15,method="linear")
mortality_func16 <- approxfun(mortality_approx$t,mortality_approx$rate16,method="linear")
mortality_func17 <- approxfun(mortality_approx$t,mortality_approx$rate17,method="linear")
mortality_func18 <- approxfun(mortality_approx$t,mortality_approx$rate18,method="linear")
mortality_func19 <- approxfun(mortality_approx$t,mortality_approx$rate19,method="linear")
mortality_func20 <- approxfun(mortality_approx$t,mortality_approx$rate20,method="linear")
mortality_func21 <- approxfun(mortality_approx$t,mortality_approx$rate21,method="linear")

birth.approx_mean <- birth.approx
birthrate_multiplier <- 0.98
birth.approx_mean$birth <- birthrate_multiplier * birth.approx_mean$birth
birth.func <- approxfun(birth.approx_mean$t,birth.approx_mean$birth,method="linear")

groupBASELINE <- groupBASE_lower
coverageBASE <- cov_lower

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE_lower
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df145 <- mutate_data(out)
results145 <- Collect_Results(df145)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df146 <- mutate_data(out)
results146 <- Collect_Results(df146)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df147 <- mutate_data(out)
results147 <- Collect_Results(df147)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df148 <- mutate_data(out)
results148 <- Collect_Results(df148)

df149 <- df145
results149 <- results145

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df150 <- mutate_data(out)
results150 <- Collect_Results(df150)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df151 <- mutate_data(out)
results151 <- Collect_Results(df151)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df152 <- mutate_data(out)
results152 <- Collect_Results(df152)

df153 <- df145
results153 <- results145

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df154 <- mutate_data(out)
results154 <- Collect_Results(df154)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df155 <- mutate_data(out)
results155 <- Collect_Results(df155)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df156 <- mutate_data(out)
results156 <- Collect_Results(df156)

df157 <- df145
results157 <- results145

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df158 <- mutate_data(out)
results158 <- Collect_Results(df158)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df159 <- mutate_data(out)
results159 <- Collect_Results(df159)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df160 <- mutate_data(out)
results160 <- Collect_Results(df160)

groupBASELINE <- groupBASE
coverageBASE <- cov_mean

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df161 <- mutate_data(out)
results161 <- Collect_Results(df161)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df162 <- mutate_data(out)
results162 <- Collect_Results(df162)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df163 <- mutate_data(out)
results163 <- Collect_Results(df163)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df164 <- mutate_data(out)
results164 <- Collect_Results(df164)

df165 <- df161
results165 <- results161

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df166 <- mutate_data(out)
results166 <- Collect_Results(df166)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df167 <- mutate_data(out)
results167 <- Collect_Results(df167)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df168 <- mutate_data(out)
results168 <- Collect_Results(df168)

df169 <- df161
results169 <- results161

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df170 <- mutate_data(out)
results170 <- Collect_Results(df170)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df171 <- mutate_data(out)
results171 <- Collect_Results(df171)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df172 <- mutate_data(out)
results172 <- Collect_Results(df172)

df173 <- df161
results173 <- results161

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df174 <- mutate_data(out)
results174 <- Collect_Results(df174)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df175 <- mutate_data(out)
results175 <- Collect_Results(df175)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df176 <- mutate_data(out)
results176 <- Collect_Results(df176)

groupBASELINE <- groupBASE_upper
coverageBASE <- cov_upper

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 15% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 25% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 15% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 25% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 15% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 25% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 15% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 25% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 50% per year

scr_group <- groupBASE_upper
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df177 <- mutate_data(out)
results177 <- Collect_Results(df177)

scr_group <- groupA1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df178 <- mutate_data(out)
results178 <- Collect_Results(df178)

scr_group <- groupA2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df179 <- mutate_data(out)
results179 <- Collect_Results(df179)

scr_group <- groupA3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df180 <- mutate_data(out)
results180 <- Collect_Results(df180)

df181 <- df177
results181 <- results177

scr_group <- groupB1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df182 <- mutate_data(out)
results182 <- Collect_Results(df182)

scr_group <- groupB2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df183 <- mutate_data(out)
results183 <- Collect_Results(df183)

scr_group <- groupB3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df184 <- mutate_data(out)
results184 <- Collect_Results(df184)

df185 <- df177
results185 <- results177

scr_group <- groupC1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df186 <- mutate_data(out)
results186 <- Collect_Results(df186)

scr_group <- groupC2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df187 <- mutate_data(out)
results187 <- Collect_Results(df187)

scr_group <- groupC3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df188 <- mutate_data(out)
results188 <- Collect_Results(df188)

df189 <- df177
results189 <- results177

scr_group <- groupD1
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df190 <- mutate_data(out)
results190 <- Collect_Results(df190)

scr_group <- groupD2
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df191 <- mutate_data(out)
results191 <- Collect_Results(df191)

scr_group <- groupD3
out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
out <- cbind(out[,1]+start_year,out)
colnames(out)[1] <- "Year"
out <- out[,colnames(out)!="time"]
df192 <- mutate_data(out)
results192 <- Collect_Results(df192)

# Collect results ####

full_results <- rbind(results1,results2,results3,results4,results5,results6,results7,results8,results9,results10,results11,results12,results13,results14,results15,results16,results17,results18,results19,results20,results21,results22,results23,results24,results25,results26,results27,results28,results29,results30,results31,results32,results33,results34,results35,results36,results37,results38,results39,results40,results41,results42,results43,results44,results45,results46,results47,results48,results49,results50,results51,results52,results53,results54,results55,results56,results57,results58,results59,results60,results61,results62,results63,results64,results65,results66,results67,results68,results69,results70,results71,results72,results73,results74,results75,results76,results77,results78,results79,results80,results81,results82,results83,results84,results85,results86,results87,results88,results89,results90,results91,results92,results93,results94,results95,results96,results97,results98,results99,results100,results101,results102,results103,results104,results105,results106,results107,results108,results109,results110,results111,results112,results113,results114,results115,results116,results117,results118,results119,results120,results121,results122,results123,results124,results125,results126,results127,results128,results129,results130,results131,results132,results133,results134,results135,results136,results137,results138,results139,results140,results141,results142,results143,results144,results145,results146,results147,results148,results149,results150,results151,results152,results153,results154,results155,results156,results157,results158,results159,results160,results161,results162,results163,results164,results165,results166,results167,results168,results169,results170,results171,results172,results173,results174,results175,results176,results177,results178,results179,results180,results181,results182,results183,results184,results185,results186,results187,results188,results189,results190,results191,results192)

for(j in 1:nrow(full_results)){
  for(i in 1:5){
    full_results[j,i] <- round(as.numeric(full_results[j,i]),0)
    full_results[j,7] <- round(as.numeric(full_results[j,7]),0)
    full_results[j,9] <- round(as.numeric(full_results[j,9]),0)
  }
}

for(i in 1:nrow(full_results)){
  full_results[i,6] <- as.numeric(full_results[17,5]) - as.numeric(full_results[i,5])
  full_results[i,8] <- as.numeric(full_results[17,7]) - as.numeric(full_results[i,7])
  full_results[i,10] <- as.numeric(full_results[i,9]) - as.numeric(full_results[17,9])
}

#view(full_results)

# Copy results to clipboard to paste in Excel sheet
# write.table(full_results, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)

# Main figures ####
# Figure 1 - population structure over time (2004 - 2021 data) ####

age_struc_long <- age_struc %>% pivot_longer(names_to = "age_group", cols = !1)

age_struc_long[age_struc_long == "group1"] <- "0 - 4"
age_struc_long[age_struc_long == "group2"] <- "5 - 9"
age_struc_long[age_struc_long == "group3"] <- "10 - 14"
age_struc_long[age_struc_long == "group4"] <- "15 - 19"
age_struc_long[age_struc_long == "group5"] <- "20 - 24"
age_struc_long[age_struc_long == "group6"] <- "25 - 29"
age_struc_long[age_struc_long == "group7"] <- "30 - 34"
age_struc_long[age_struc_long == "group8"] <- "35 - 39"
age_struc_long[age_struc_long == "group9"] <- "40 - 44"
age_struc_long[age_struc_long == "group10"] <- "45 - 49"
age_struc_long[age_struc_long == "group11"] <- "50 - 54"
age_struc_long[age_struc_long == "group12"] <- "55 - 59"
age_struc_long[age_struc_long == "group13"] <- "60 - 64"
age_struc_long[age_struc_long == "group14"] <- "65 - 69"
age_struc_long[age_struc_long == "group15"] <- "70 - 74"
age_struc_long[age_struc_long == "group16"] <- "75 - 79"
age_struc_long[age_struc_long == "group17"] <- "80 - 84"
age_struc_long[age_struc_long == "group18"] <- "85 - 89"
age_struc_long[age_struc_long == "group19"] <- "90 - 94"
age_struc_long[age_struc_long == "group20"] <- "95 - 99"
age_struc_long[age_struc_long == "group21"] <- "Over 100"

age_struc_long$age_group <- factor(age_struc_long$age_group, levels = c('0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 - 39', '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69', '70 - 74', '75 - 79', '80 - 84', '85 - 89', '90 - 94', '95 - 99', 'Over 100'))

age_struc_long %>% 
  ggplot(aes(x=age_group,y=value,fill=as.factor(Year))) +
  geom_bar(stat="identity",position="dodge2", width=0.99) +
  labs(title="Population Structure Over Time", x="Age Group",y="Proportion of Population")+
  theme_bw(base_size = 18) +
  theme(
    strip.text.x = element_text(
      size = 14, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16)
  )+
  guides(fill=guide_legend(nrow=10,byrow=TRUE,title="Year"))+
  theme(plot.title = element_text(face = "italic", size=15))+
  scale_fill_viridis_d(direction=-1,option="D")


# Figure 3 - Total population of Thailand by scenario ####

dataframes_fig3 <- list(df1, df17, df33, df49, df65, df81, df97, df113, df129, df145, df161, df177)
names(dataframes_fig3) <- list('df1', 'df17', 'df33', 'df49', 'df65', 'df81', 'df97', 'df113', 'df129', 'df145', 'df161', 'df177')

dataframe_fig3 <- list.rbind(combine_dataframes(dataframes_fig3))

pop_totals <- dataframe_fig3 %>% filter(variable %in% c("total")) 

# Remove scr_group and scr_cov as they're the same for all data included in this figure
pop_totals <- pop_totals[-7]
pop_totals <- pop_totals[-7]

names(pop_data1)[names(pop_data1) == colnames(pop_data1)[1]] <- "Year"
names(pop_data1)[names(pop_data1) == colnames(pop_data1)[2]] <- "value"
names(pop_data1)[names(pop_data1) == colnames(pop_data1)[3]] <- "type"

pop_totals <- pivot_wider(pop_totals, names_from = baseline_coverage)

pop0 <- pop_totals[pop_totals$population_scenario=="baseline",]
pop1 <- pop_totals[pop_totals$population_scenario=="decline",]
pop2 <- pop_totals[pop_totals$population_scenario=="growth",]
pop3 <- pop_totals[pop_totals$population_scenario=="plateau",]

pop_totals %>% 
  
  ggplot() +
  
  geom_line(size = 1,
            aes(x = Year, y = mean, colour = population_scenario)) +
  
  geom_ribbon(data=pop0,
              mapping=aes(x = Year,ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "deeppink",
              alpha = 0.4
  ) +
  geom_ribbon(data=pop1,
              mapping=aes(x = Year,ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "lightgreen",
              alpha = 0.4
  ) +
  
  geom_ribbon(data=pop2,
              mapping=aes(x = Year,ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "turquoise2",
              alpha = 0.4
  ) +
  geom_ribbon(data=pop3,
              mapping=aes(x = Year,ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "darkorchid1",
              alpha = 0.4
  ) +
  
  geom_point(size = 2,
             colour = "gray43",
             data = pop_data1,
             aes(x = Year, y = value, shape = type)) +
  theme_bw(base_size = 18)+
  theme(
    strip.text.x = element_text(
      size = 14, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16)
  ) +
  scale_shape_discrete(labels=c('Total Population Data', 'Total Population Projection'))+
  scale_colour_discrete(labels=c('Baseline','Decline','Growth','Plateau'))+
  scale_y_continuous(labels=scales::comma) +
  theme(plot.title = element_text(face = "italic", size=15))+
  labs(title = "Population of Thailand", x="Year", y =("Size of Population"), colour="Model Output", shape="Data Point")


# Figure 4 - Age group populations of Thailand by scenario ####

pop_groups <- dataframe_fig3 %>% filter(variable %in% c("group1","group2","group3","group4","group5","group6","group7","group8","group9","group10","group11","group12","group13","group14","group15","group16","group17","group18","group19","group20","group21")) 

pop_groups[pop_groups == "group1"] <- "0 - 4"
pop_groups[pop_groups == "group2"] <- "5 - 9"
pop_groups[pop_groups == "group3"] <- "10 - 14"
pop_groups[pop_groups == "group4"] <- "15 - 19"
pop_groups[pop_groups == "group5"] <- "20 - 24"
pop_groups[pop_groups == "group6"] <- "25 - 29"
pop_groups[pop_groups == "group7"] <- "30 - 34"
pop_groups[pop_groups == "group8"] <- "35 - 39"
pop_groups[pop_groups == "group9"] <- "40 - 44"
pop_groups[pop_groups == "group10"] <- "45 - 49"
pop_groups[pop_groups == "group11"] <- "50 - 54"
pop_groups[pop_groups == "group12"] <- "55 - 59"
pop_groups[pop_groups == "group13"] <- "60 - 64"
pop_groups[pop_groups == "group14"] <- "65 - 69"
pop_groups[pop_groups == "group15"] <- "70 - 74"
pop_groups[pop_groups == "group16"] <- "75 - 79"
pop_groups[pop_groups == "group17"] <- "80 - 84"
pop_groups[pop_groups == "group18"] <- "85 - 89"
pop_groups[pop_groups == "group19"] <- "90 - 94"
pop_groups[pop_groups == "group20"] <- "95 - 99"
pop_groups[pop_groups == "group21"] <- "Over 100"

pop_groups$variable <- factor(pop_groups$variable, levels = c('0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 - 39', '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69', '70 - 74', '75 - 79', '80 - 84', '85 - 89', '90 - 94', '95 - 99', 'Over 100'))

names(age_struc_proportion_data)[names(age_struc_proportion_data) == colnames(age_struc_proportion_data)[1]] <- "Year"
UN_data_groups <- age_struc_proportion_data %>% pivot_longer(names_to = "variable", cols = !1)
UN_data_groups$type <- rep("data",nrow(UN_data_groups))
names(age_struc_proportion_proj)[names(age_struc_proportion_proj) == colnames(age_struc_proportion_proj)[1]] <- "Year"
age_struc_proportion_proj <- age_struc_proportion_proj %>% filter(Year %in% c("2022","2023","2024","2025","2026","2027","2028","2029","2030","2031","2032","2033","2034","2035","2036","2037","2038","2039","2040"))
proj_groups <- age_struc_proportion_proj %>% pivot_longer(names_to = "variable", cols = !1)
proj_groups$type <- rep("projection",nrow(proj_groups))
data_groups <- rbind(UN_data_groups,proj_groups)

data_groups[data_groups == "group1"] <- "0 - 4"
data_groups[data_groups == "group2"] <- "5 - 9"
data_groups[data_groups == "group3"] <- "10 - 14"
data_groups[data_groups == "group4"] <- "15 - 19"
data_groups[data_groups == "group5"] <- "20 - 24"
data_groups[data_groups == "group6"] <- "25 - 29"
data_groups[data_groups == "group7"] <- "30 - 34"
data_groups[data_groups == "group8"] <- "35 - 39"
data_groups[data_groups == "group9"] <- "40 - 44"
data_groups[data_groups == "group10"] <- "45 - 49"
data_groups[data_groups == "group11"] <- "50 - 54"
data_groups[data_groups == "group12"] <- "55 - 59"
data_groups[data_groups == "group13"] <- "60 - 64"
data_groups[data_groups == "group14"] <- "65 - 69"
data_groups[data_groups == "group15"] <- "70 - 74"
data_groups[data_groups == "group16"] <- "75 - 79"
data_groups[data_groups == "group17"] <- "80 - 84"
data_groups[data_groups == "group18"] <- "85 - 89"
data_groups[data_groups == "group19"] <- "90 - 94"
data_groups[data_groups == "group20"] <- "95 - 99"
data_groups[data_groups == "group21"] <- "Over 100"

data_groups$variable <- factor(data_groups$variable, levels = c('0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 - 39', '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69', '70 - 74', '75 - 79', '80 - 84', '85 - 89', '90 - 94', '95 - 99', 'Over 100'))

# Remove scr_group and scr_cov as they're the same for all data included in this figure
pop_groups <- pop_groups[-7]
pop_groups <- pop_groups[-7]

pop_groups <- pivot_wider(pop_groups, names_from = baseline_coverage)

group0 <- pop_groups[pop_groups$population_scenario=="baseline",]
group1 <- pop_groups[pop_groups$population_scenario=="decline",]
group2 <- pop_groups[pop_groups$population_scenario=="growth",]
group3 <- pop_groups[pop_groups$population_scenario=="plateau",]

pop_groups %>% 
  
  ggplot() +
  
  geom_line(size = 1,
            aes(x = Year, y = mean, colour = population_scenario)) +
  
  geom_ribbon(data=group0,
              mapping=aes(x = Year,ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "deeppink",
              alpha = 0.2
  ) +
  geom_ribbon(data=group1,
              mapping=aes(x = Year,ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "lightgreen",
              alpha = 0.2
  ) +
  
  geom_ribbon(data=group2,
              mapping=aes(x = Year,ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "turquoise2",
              alpha = 0.2
  ) +
  geom_ribbon(data=group3,
              mapping=aes(x = Year,ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "darkorchid1",
              alpha = 0.2
  ) +
  
  geom_point(size = 1.5,
             colour = "gray43",
             data = data_groups,
             aes(x = Year, y = value, shape = type)) +
  
  facet_wrap(~variable) +
  scale_shape_discrete(labels=c('Total Population Data', 'Total Population Projection'))+
  scale_colour_discrete(labels=c('Baseline','Decline','Growth','Plateau'))+
  scale_y_continuous(labels=scales::comma) +
  theme_bw(base_size = 18)+
  theme(
    strip.text.x = element_text(
      size = 14, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16)
  )+
  theme(plot.title = element_text(face = "italic", size=15))+
  labs(title = "Model vs. Data by Age Group", x="Year", y =("Size of Population"), shape="Data Point", colour="Model Output")

# Figure 5 - Prevalence by age group, model vs. data####

dataframes_fig5 <- list(df1, df17, df33)
names(dataframes_fig5) <- list("df1", "df17", "df33")
dataframe_fig5 <- list.rbind(combine_dataframes(dataframes_fig5))

model_prev <- dataframe_fig5 %>% filter(variable %in% c("prev1_2","prev3_4","prev5_6","prev7_8","prev9_10","prev11_21", "infectprev"))
model_prev <- model_prev[-4]
model_prev <- model_prev[-6]
model_prev <- model_prev[-6]
model_prev <- model_prev[-5]

model_prev <- pivot_wider(model_prev, names_from = baseline_coverage)

model_prev[model_prev == "infectprev"] <- "Total"
model_prev[model_prev == "prev1_2"] <- "0 - 9"
model_prev[model_prev == "prev3_4"] <- "10 - 19"
model_prev[model_prev == "prev5_6"] <- "20 - 29"
model_prev[model_prev == "prev7_8"] <- "30 - 39"
model_prev[model_prev == "prev9_10"] <- "40 - 49"
model_prev[model_prev == "prev11_21"] <- "Over 50"

model_prev$variable <- factor(model_prev$variable, levels = c('0 - 9','10 - 19','20 - 29','30 - 39','40 - 49','Over 50','Total'))

data_prev$lower_95 <- data_prev$value-data_prev$CI95
data_prev$upper_95 <- data_prev$value+data_prev$CI95
data_prev <- select(data_prev, -4,-5,-6)
data_prev$type <- rep("data", nrow(data_prev))
data_prev <- subset(data_prev, value != 0)

data_prev[data_prev == "infectprev"] <- "Total"
data_prev[data_prev == "prev1_2"] <- "0 - 9"
data_prev[data_prev == "prev3_4"] <- "10 - 19"
data_prev[data_prev == "prev5_6"] <- "20 - 29"
data_prev[data_prev == "prev7_8"] <- "30 - 39"
data_prev[data_prev == "prev9_10"] <- "40 - 49"
data_prev[data_prev == "prev11_21"] <- "Over 50"

data_prev$variable <- factor(data_prev$variable, levels = c('0 - 9','10 - 19','20 - 29','30 - 39','40 - 49','Over 50','Total'))


model_prev %>% 
  ggplot() +
  
  geom_line(aes(x = Year, y = mean)) +
  
  geom_ribbon(
    mapping=aes(x=Year,
                ymin = lower,
                ymax = upper),
    fill = "blue",
    alpha = 0.3
  ) +
  
  geom_point(data=data_prev,
             aes(x=Year, y=value)) +
  
  geom_errorbar(data_prev,
                mapping = aes(x = Year, ymin = lower_95, ymax= upper_95)) +
  theme_bw(base_size = 18)+
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16),
    legend.position = "bottom"
  )+
  
  facet_wrap(~variable) +
  scale_y_continuous(labels=scaleFUN, limits=c(0,4)) +
  theme(plot.title = element_text(face = "italic", size=15))+
  labs(title="Prevalence by Age Group: Model vs. Data", x="Year",y="Thailand Prevalence (%)") 

# Figure 6 - Main result: Incidence and mortality vs WHO Targets ####

dataframes_fig6 <- list(df1, df4, df17, df20, df33, df36)
names(dataframes_fig6) <- list('df1', 'df4', 'df17', 'df20', 'df33', 'df36')

dataframe_fig6 <- list.rbind(combine_dataframes(dataframes_fig6))

model_incidence <- dataframe_fig6 %>% filter(variable %in% c("Inc"))
model_incidence <- model_incidence[-4]
model_incidence <- model_incidence[-5]

model_incidence <- pivot_wider(model_incidence, names_from = baseline_coverage)

inc_2015_lower <- as.numeric((df1 %>% filter(variable %in% c("Inc")) %>% filter(Year %in% c("2015")))[,3])
inc_target_lower <- 0.1*inc_2015_lower

inc_2015_mean <- as.numeric((df17 %>% filter(variable %in% c("Inc")) %>% filter(Year %in% c("2015")))[,3])
inc_target_mean <- 0.1*inc_2015_mean

inc_2015_upper <- as.numeric((df33 %>% filter(variable %in% c("Inc")) %>% filter(Year %in% c("2015")))[,3])
inc_target_upper <- 0.1*inc_2015_upper

model_mortality <- dataframe_fig6 %>% filter(variable %in% c("Deaths"))
model_mortality <- model_mortality[-4]
model_mortality <- model_mortality[-5]

model_mortality <- pivot_wider(model_mortality, names_from = baseline_coverage)

mort_2015_lower <- as.numeric((df1 %>% filter(variable %in% c("Deaths")) %>% filter(Year %in% c("2015")))[,3])
mort_target_lower <- 0.35*mort_2015_lower

mort_2015_mean <- as.numeric((df17 %>% filter(variable %in% c("Deaths")) %>% filter(Year %in% c("2015")))[,3])
mort_target_mean <- 0.35*mort_2015_mean

mort_2015_upper <- as.numeric((df33 %>% filter(variable %in% c("Deaths")) %>% filter(Year %in% c("2015")))[,3])
mort_target_upper <- 0.35*mort_2015_upper

inc_baseline <- model_incidence %>% filter(scr_group %in% c("30+")) %>% filter(scr_cov %in% c("baseline"))
inc_max_scr <- model_incidence %>% filter(scr_group %in% c("30+")) %>% filter(scr_cov %in% c("0.5"))
mort_baseline <- model_mortality %>% filter(scr_group %in% c("30+")) %>% filter(scr_cov %in% c("baseline"))
mort_max_scr <- model_mortality %>% filter(scr_group %in% c("30+")) %>% filter(scr_cov %in% c("0.5"))

fig6a <- inc_baseline %>% filter(Year>2018) %>% 
  ggplot(aes(x = Year, y = mean)) +
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "deeppink",
              alpha = 0.2
  ) +
  geom_line(size = 1.3)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  geom_vline(data = data.frame(xint=2030), aes(xintercept = xint), linetype = "dashed") +
     annotate("rect", xmin = 2018, xmax = 2042, ymin = inc_target_lower, ymax = inc_target_upper,
              alpha = .7,
              fill="lightgrey") +
     annotate("segment", x = 2018, xend = 2042, y = inc_target_mean, yend = inc_target_mean,
              colour = "grey30", linetype=2)+
     annotate("segment", x = 2030, xend = 2030, y = 0, yend = 10000,
              colour = "grey30", linetype=2)+
     annotate("text", x = 2025, y = inc_target_mean-500,family = "", fontface = 3, size=5, label = "2030 Elimination Target (95% CI)") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "italic", size=15))+
  coord_cartesian(xlim = c(2020,2040), ylim = c(0, 12000)) +
  labs(title = "Incidence with Baseline Screening: 6.2% Across All Ages", x="Year", y ="New Cases")

fig6b <- inc_max_scr %>% filter(Year>2018) %>% 
  ggplot(aes(x = Year, y = mean)) +
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "darkorchid1",
              alpha = 0.2
  ) +
  geom_line(size = 1.3)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  geom_vline(data = data.frame(xint=2030), aes(xintercept = xint), linetype = "dashed") +
  annotate("rect", xmin = 2018, xmax = 2042, ymin = inc_target_lower, ymax = inc_target_upper,
           alpha = .7,
           fill="lightgrey") +
  annotate("segment", x = 2018, xend = 2042, y = inc_target_mean, yend = inc_target_mean,
           colour = "grey30", linetype=2)+
  annotate("segment", x = 2030, xend = 2030, y = 0, yend = 10000,
           colour = "grey30", linetype=2)+
  annotate("text", x = 2025, y = inc_target_mean-500,family = "", fontface = 3, size=5, label = "2030 Elimination Target (95% CI)") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "italic", size=15))+
  coord_cartesian(xlim = c(2020,2040), ylim = c(0, 12000)) +
  annotate(geom = "point", x = 2023, y = 7400, colour = "darkorchid1", size = 5) + 
  annotate(geom = "point", x = 2023, y = 7400, size=3) + 
  annotate(geom = "text", x = 2024, y = 7500, label = "Start of Screening Programme", hjust = "left",fontface = 3, size=5) +
  labs(title = "Incidence with Maximum Screening: 50% Across Ages 30+", x="Year", y ="") 

fig6c <- mort_baseline %>% filter(Year>2018) %>% 
  ggplot(aes(x = Year, y = mean)) +
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "turquoise2",
              alpha = 0.2
  ) +
  geom_line(size = 1.3)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  geom_vline(data = data.frame(xint=2030), aes(xintercept = xint), linetype = "dashed") +
  annotate("rect", xmin = 2018, xmax = 2042, ymin = mort_target_lower, ymax = mort_target_upper,
           alpha = .7,
           fill="lightgrey") +
  annotate("segment", x = 2018, xend = 2042, y = mort_target_mean, yend = mort_target_mean,
           colour = "grey30", linetype=2)+
  annotate("segment", x = 2030, xend = 2030, y = 0, yend = 10000,
           colour = "grey30", linetype=2)+
  annotate("text", x = 2025, y = mort_target_mean-500,family = "", fontface = 3, size=5, label = "2030 Elimination Target (95% CI)") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "italic", size=15))+
  coord_cartesian(xlim = c(2020,2040), ylim = c(0, 12000)) +
  labs(title = "Mortality with Baseline Screening: 6.2% Across All Ages", x="Year", y ="HCV-related Deaths")

fig6d <- mort_max_scr %>% filter(Year>2018) %>% 
  ggplot(aes(x = Year, y = mean)) +
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "lightgreen",
              alpha = 0.2
  ) +
  geom_line(size = 1.3)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  geom_vline(data = data.frame(xint=2030), aes(xintercept = xint), linetype = "dashed") +
  annotate("rect", xmin = 2018, xmax = 2042, ymin = mort_target_lower, ymax = mort_target_upper,
           alpha = .7,
           fill="lightgrey") +
  annotate("segment", x = 2018, xend = 2042, y = mort_target_mean, yend = mort_target_mean,
           colour = "grey30", linetype=2)+
  annotate("segment", x = 2030, xend = 2030, y = 0, yend = 10000,
           colour = "grey30", linetype=2)+
  annotate("text", x = 2025, y = mort_target_mean-500,family = "", fontface = 3, size=5, label = "2030 Elimination Target (95% CI)") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "italic", size=15))+
  coord_cartesian(xlim = c(2020,2040), ylim = c(0, 12000)) +
  annotate(geom = "point", x = 2023, y = 7900, colour = "lightgreen", size = 5) + 
  annotate(geom = "point", x = 2023, y = 7900, size=3) + 
  annotate(geom = "text", x = 2024, y = 8500, label = "Start of Screening Programme", hjust = "left",fontface = 3, size=5) +
  labs(title = "Mortality with Maximum Screening: 50% Across Ages 30+", x="Year", y ="")

grid.arrange(fig6a, fig6b, fig6c, fig6d, ncol=2, nrow=2,
             top = textGrob("Yearly Incidence and Mortality vs WHO Targets",gp=gpar(fontsize=20,font=3)))

# Supplementary Figure 1 - birth rate and population data and projection ####

pop_data <- cbind(pop_data,rep("data",nrow(pop_data)))
pop_proj <- cbind(pop_proj,rep("proj",nrow(pop_proj)))
names(pop_data)[names(pop_data) == colnames(pop_data)[1]] <- "Year"
names(pop_data)[names(pop_data) == colnames(pop_data)[2]] <- "value"
names(pop_data)[names(pop_data) == colnames(pop_data)[3]] <- "type"
names(pop_proj)[names(pop_proj) == colnames(pop_proj)[1]] <- "Year"
names(pop_proj)[names(pop_proj) == colnames(pop_proj)[2]] <- "value"
names(pop_proj)[names(pop_proj) == colnames(pop_proj)[3]] <- "type"
pop_data$type2 <- rep("pop", nrow(pop_data))
pop_proj$type2 <- rep("pop", nrow(pop_proj))

birth.func <- approxfun(birth.approx$t,birth.approx$birth,method="linear")
birthrate_data <- birth.func(0:18)
birthrate_proj <- birth.func(19:36)
birthrate_data <- cbind(c(2004:2022),birthrate_data,rep("data",length(birthrate_data)),rep("birth",length(birthrate_data)))
birthrate_proj <- cbind(c(2023:2040),birthrate_proj,rep("proj",length(birthrate_proj)),rep("birth",length(birthrate_proj)))
colnames(birthrate_data) <- c("Year", "value", "type", "type2")
colnames(birthrate_proj) <- c("Year"," value", "type", "type2")


birthrate_compare <- as.data.frame(rbind(birthrate_data,birthrate_proj))
birthrate_compare[,1] <- as.numeric(birthrate_compare[,1])
birthrate_compare[,2] <- as.numeric(birthrate_compare[,2])

pop_and_birth <- rbind(pop_data,pop_proj,birthrate_compare)

ggplot(pop_and_birth, aes(x = Year, y = value, colour = type)) +
  geom_point(size=2) +
  theme_minimal(base_size=14) +
  theme(legend.title=element_blank())+
  scale_y_continuous(labels=scales::comma) +
  facet_wrap(~type2, scales = "free_y", nrow = 2, 
             strip.position = "left", 
             labeller = as_labeller(c(pop = "Population", birth = "Birth Rate (Births per Person per Year)") ) )  +
  ylab(NULL) +
  theme_bw(base_size = 18) +
  theme(
    strip.text.x = element_text(
      size = 14, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16)
  )+
  scale_shape_discrete(labels=c('Birth Rate','Total Population'))+
  scale_colour_discrete(labels=c('Data', 'Projection')) +
  theme(plot.title = element_text(face = "italic", size=18))+
  labs(title = "Birth Rate and Population", x="Year", y =(""), colour="Data Type")

# Supplementary Figure 2 - beta matrix heat map ####

beta_matrix <- as.data.frame(cbind(c(1:21),beta))

beta_matrix  <- beta_matrix  %>% pivot_longer(names_to = "group", cols = !1)
beta_matrix  <- cbind(beta_matrix [,1],rep(1:21,21),beta_matrix [,3])
names(beta_matrix )[names(beta_matrix ) == colnames(beta_matrix )[1]] <- "X"
names(beta_matrix )[names(beta_matrix ) == colnames(beta_matrix )[2]] <- "Y"
names(beta_matrix )[names(beta_matrix ) == colnames(beta_matrix )[3]] <- "Z"
beta_matrix  <- beta_matrix  %>% mutate(X=5*X-1, Y=5*Y-1)

ggplot(beta_matrix , aes(Y, X, fill= Z)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, direction=1) +
  labs(title = "Beta Matrix: Derived from Sexual Contact and HPV in Laos", x="Age", y =("Age"), fill = "Transmission Coefficient") +
  theme_bw(base_size = 18)+
  theme(plot.title = element_text(face = "italic", size=15))

# Supplementary Figure 3 - Sensitivity analysis: population scenario on 2030 incidence ####

dataframes_suppfig3 <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, 
                            df16, df17, df18, df19, df20, df21, df22, df23, df24, df25, df26, df27, df28, df29, 
                            df30, df31, df32, df33, df34, df35, df36, df37, df38, df39, df40, df41, df42, df43, 
                            df44, df45, df46, df47, df48, df49, df50, df51, df52, df53, df54, df55, df56, df57, 
                            df58, df59, df60, df61, df62, df63, df64, df65, df66, df67, df68, df69, df70, df71, 
                            df72, df73, df74, df75, df76, df77, df78, df79, df80, df81, df82, df83, df84, df85, 
                            df86, df87, df88, df89, df90, df91, df92, df93, df94, df95, df96, df97, df98, df99, 
                            df100, df101, df102, df103, df104, df105, df106, df107, df108, df109, df110, df111, 
                            df112, df113, df114, df115, df116, df117, df118, df119, df120, df121, df122, df123, 
                            df124, df125, df126, df127, df128, df129, df130, df131, df132, df133, df134, df135, 
                            df136, df137, df138, df139, df140, df141, df142, df143, df144, df145, df146, df147, 
                            df148, df149, df150, df151, df152, df153, df154, df155, df156, df157, df158, df159, 
                            df160, df161, df162, df163, df164, df165, df166, df167, df168, df169, df170, df171, 
                            df172, df173, df174, df175, df176, df177, df178, df179, df180, df181, df182, df183, 
                            df184, df185, df186, df187, df188, df189, df190, df191, df192)

names(dataframes_suppfig3) <- list('df1', 'df2', 'df3', 'df4', 'df5', 'df6', 'df7', 'df8', 'df9', 'df10', 'df11',
                                   'df12', 'df13', 'df14', 'df15', 'df16', 'df17', 'df18', 'df19', 'df20', 'df21', 
                                   'df22', 'df23', 'df24', 'df25', 'df26', 'df27', 'df28', 'df29', 'df30', 'df31', 
                                   'df32', 'df33', 'df34', 'df35', 'df36', 'df37', 'df38', 'df39', 'df40', 'df41', 
                                   'df42', 'df43', 'df44', 'df45', 'df46', 'df47', 'df48', 'df49', 'df50', 'df51', 
                                   'df52', 'df53', 'df54', 'df55', 'df56', 'df57', 'df58', 'df59', 'df60', 'df61', 
                                   'df62', 'df63', 'df64', 'df65', 'df66', 'df67', 'df68', 'df69', 'df70', 'df71', 
                                   'df72', 'df73', 'df74', 'df75', 'df76', 'df77', 'df78', 'df79', 'df80', 'df81', 
                                   'df82', 'df83', 'df84', 'df85', 'df86', 'df87', 'df88', 'df89', 'df90', 'df91', 
                                   'df92', 'df93', 'df94', 'df95', 'df96', 'df97', 'df98', 'df99', 'df100', 'df101', 
                                   'df102', 'df103', 'df104', 'df105', 'df106', 'df107', 'df108', 'df109', 'df110', 
                                   'df111', 'df112', 'df113', 'df114', 'df115', 'df116', 'df117', 'df118', 'df119', 
                                   'df120', 'df121', 'df122', 'df123', 'df124', 'df125', 'df126', 'df127', 'df128', 
                                   'df129', 'df130', 'df131', 'df132', 'df133', 'df134', 'df135', 'df136', 'df137', 
                                   'df138', 'df139', 'df140', 'df141', 'df142', 'df143', 'df144', 'df145', 'df146', 
                                   'df147', 'df148', 'df149', 'df150', 'df151', 'df152', 'df153', 'df154', 'df155', 
                                   'df156', 'df157', 'df158', 'df159', 'df160', 'df161', 'df162', 'df163', 'df164', 
                                   'df165', 'df166', 'df167', 'df168', 'df169', 'df170', 'df171', 'df172', 'df173', 
                                   'df174', 'df175', 'df176', 'df177', 'df178', 'df179', 'df180', 'df181', 'df182', 
                                   'df183', 'df184', 'df185', 'df186', 'df187', 'df188', 'df189', 'df190', 'df191', 
                                   'df192')

dataframe_suppfig3 <- list.rbind(combine_dataframes(dataframes_suppfig3))

incidence_2030 <- dataframe_suppfig3 %>% filter(variable %in% c("Inc")) %>% filter(Year %in% c("2030"))

i1 <- incidence_2030 %>% filter(population_scenario %in% c("baseline"))
i1 <- pivot_wider(i1, names_from = baseline_coverage)
i2 <- incidence_2030 %>% filter(population_scenario %in% c("decline"))
i2 <- pivot_wider(i2, names_from = baseline_coverage)
i3 <- incidence_2030 %>% filter(population_scenario %in% c("growth"))
i3 <- pivot_wider(i3, names_from = baseline_coverage)
i4 <- incidence_2030 %>% filter(population_scenario %in% c("plateau"))
i4 <- pivot_wider(i4, names_from = baseline_coverage)
incidence_2030 <- rbind(i1,i2,i3,i4)

incidence_2030$population_scenario[incidence_2030$population_scenario == "baseline"] <- "Baseline"
incidence_2030$population_scenario[incidence_2030$population_scenario == "decline"] <- "Decline"
incidence_2030$population_scenario[incidence_2030$population_scenario == "growth"] <- "Growth"
incidence_2030$population_scenario[incidence_2030$population_scenario == "plateau"] <- "Plateau"
incidence_2030$scr_cov <- factor(incidence_2030$scr_cov, levels = c('baseline', '0.15', '0.25', '0.5'))

incidence_2030 %>% 
  ggplot(aes(x = population_scenario,
             y = mean,
             fill = factor(scr_cov))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))+
  facet_wrap(~scr_group) +
  theme_bw(base_size = 18) +
  annotate("text",
           x = 2,
           y = inc_target_mean +200,
           label = "2030 Incidence Target",
           family = "", fontface = 3, size=4) +
  annotate("segment", x = 0, xend = 5, y = inc_target_mean, yend = inc_target_mean,
           colour = "grey30", linetype=2)+
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16),
    legend.position = "bottom"
  ) +
  scale_fill_discrete(labels=c('6.2% (Baseline)', '15%', '25%', '50%')) +
  theme(plot.title = element_text(face = "italic", size=15))+
  labs(title="Sensitivity Analysis: Effect of Population Scenario on 2030 Incidence", x="Population Scenario", y =("New Cases in 2030"), fill="Yearly Screening Coverage")


# Supplementary figure 4 - Sensitivity analysis: population scenario on 2030 Mortality ####

mortality_2030 <- dataframe_suppfig3 %>% filter(variable %in% c("Deaths")) %>% filter(Year %in% c("2030"))

n1 <- mortality_2030 %>% filter(population_scenario %in% c("baseline"))
n1 <- pivot_wider(n1, names_from = baseline_coverage)
n2 <- mortality_2030 %>% filter(population_scenario %in% c("decline"))
n2 <- pivot_wider(n2, names_from = baseline_coverage)
n3 <- mortality_2030 %>% filter(population_scenario %in% c("growth"))
n3 <- pivot_wider(n3, names_from = baseline_coverage)
n4 <- mortality_2030 %>% filter(population_scenario %in% c("plateau"))
n4 <- pivot_wider(n4, names_from = baseline_coverage)
mortality_2030 <- rbind(n1,n2,n3,n4)

mortality_2030$population_scenario[mortality_2030$population_scenario == "baseline"] <- "Baseline"
mortality_2030$population_scenario[mortality_2030$population_scenario == "decline"] <- "Decline"
mortality_2030$population_scenario[mortality_2030$population_scenario == "growth"] <- "Growth"
mortality_2030$population_scenario[mortality_2030$population_scenario == "plateau"] <- "Plateau"
mortality_2030$scr_cov <- factor(mortality_2030$scr_cov, levels = c('baseline', '0.15', '0.25', '0.5'))

mortality_2030 %>% 
  ggplot(aes(x = population_scenario,
             y = mean,
             fill = factor(scr_cov))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))+
  facet_wrap(~scr_group) +
  theme_bw(base_size = 18) +
  annotate("text",
           x = 2,
           y = mort_target_mean +500,
           label = "2030 Mortality Target",
           family = "", fontface = 3, size=4) +
  annotate("segment", x = 0, xend = 5, y = mort_target_mean, yend = mort_target_mean,
           colour = "grey30", linetype=2)+
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16),
    legend.position = "bottom"
  ) +
  scale_fill_discrete(labels=c('6.2% (Baseline)', '15%', '25%', '50%')) +
  theme(plot.title = element_text(face = "italic", size=15))+
  labs(title="Sensitivity Analysis: Effect of Population Scenario on 2030 Mortality", x="Population Scenario", y =("HCV Related Deaths in 2030"), fill="Yearly Screening Coverage")


# Supplementary figure 5 - Sensitivity analysis: population scenario on 2030 Prevalence ####

prev_2030 <- dataframe_suppfig3 %>% filter(variable %in% c("prev1_2","prev3_4","prev5_6","prev7_8","prev9_10","prev11_21", "infectprev")) %>% filter(Year %in% c("2030")) %>% 
  filter(scr_group %in% c("30+")) %>% filter(scr_cov %in% c("baseline","0.5"))

p1 <- prev_2030 %>% filter(population_scenario %in% c("baseline"))
p1 <- pivot_wider(p1, names_from = baseline_coverage)
p2 <- prev_2030 %>% filter(population_scenario %in% c("decline"))
p2 <- pivot_wider(p2, names_from = baseline_coverage)
p3 <- prev_2030 %>% filter(population_scenario %in% c("growth"))
p3 <- pivot_wider(p3, names_from = baseline_coverage)
p4 <- prev_2030 %>% filter(population_scenario %in% c("plateau"))
p4 <- pivot_wider(p4, names_from = baseline_coverage)
prev_2030 <- rbind(p1,p2,p3,p4)

prev_2030$scr_cov <- factor(prev_2030$scr_cov, levels = c('baseline', '0.5'))
prev_2030$variable <- factor(prev_2030$variable, levels = c('prev1_2' ,'prev3_4', 'prev5_6', 'prev7_8', 'prev9_10', 'prev11_21', 'infectprev'))

prev_2030$population_scenario[prev_2030$population_scenario == "baseline"] <- "Baseline"
prev_2030$population_scenario[prev_2030$population_scenario == "decline"] <- "Decline"
prev_2030$population_scenario[prev_2030$population_scenario == "growth"] <- "Growth"
prev_2030$population_scenario[prev_2030$population_scenario == "plateau"] <- "Plateau"

prev_2030 %>% 
  ggplot(aes(x = population_scenario,
             y = mean,
             fill = factor(variable))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))+
  facet_wrap(~scr_cov, labeller = labeller(scr_cov = 
                                             c("baseline" = "Baseline Coverage: 6.2% per Year Across All Age Groups",
                                               "0.5" = "Heaviest Coverage: 50% per Year Across Ages 30+")
  ))  +
  theme_bw(base_size = 18) +
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16),
    legend.position = "bottom"
  )+
  scale_fill_discrete(labels=c('0 - 9', '10 - 19', '20 - 29', '30 - 39', '40 - 49', 'Over 50', 'Total')) +
  theme(plot.title = element_text(face = "italic", size=15))+
  scale_y_continuous(labels=scaleFUN, limits=c(0,1.3)) +
  labs(title="Sensitivity Analysis: Effect of Population Scenario on 2030 Prevalence", x="Population Scenario", y =("Prevalence in Thailand in 2030 (%)"), fill="Age Group")


# Supplementary figure 6 - Main result (incidence vs WHO goals) for all screening groups and coverages ####

dataframes_suppfig6 <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, 
                            df16, df17, df18, df19, df20, df21, df22, df23, df24, df25, df26, df27, df28, df29, 
                            df30, df31, df32, df33, df34, df35, df36, df37, df38, df39, df40, df41, df42, df43, 
                            df44, df45, df46, df47, df48)
names(dataframes_suppfig6) <- list('df1', 'df2', 'df3', 'df4', 'df5', 'df6', 'df7', 'df8', 'df9', 'df10', 
                                   'df11', 'df12', 'df13', 'df14', 'df15', 'df16', 'df17', 'df18', 'df19', 'df20', 
                                   'df21', 'df22', 'df23', 'df24', 'df25', 'df26', 'df27', 'df28', 'df29', 'df30', 
                                   'df31', 'df32', 'df33', 'df34', 'df35', 'df36', 'df37', 'df38', 'df39', 'df40', 
                                   'df41', 'df42', 'df43', 'df44', 'df45', 'df46', 'df47', 'df48')

dataframe_suppfig6 <- list.rbind(combine_dataframes(dataframes_suppfig6))

model_incidence_all <- dataframe_suppfig6 %>% filter(variable %in% c("Inc"))
model_incidence_all <- model_incidence_all[-4]
model_incidence_all <- model_incidence_all[-5]
model_incidence_all <- pivot_wider(model_incidence_all, names_from = baseline_coverage)

model_incidence_all$scr_cov <- as.character(model_incidence_all$scr_cov)
model_incidence_all$scr_group <- as.character(model_incidence_all$scr_group)
model_incidence_all[model_incidence_all == "baseline"] <- "6.2% (Baseline)"
model_incidence_all[model_incidence_all == "0.15"] <- "15%"
model_incidence_all[model_incidence_all == "0.25"] <- "25%"
model_incidence_all[model_incidence_all == "0.5"] <- "50%"

model_incidence_all$scr_cov <- factor(model_incidence_all$scr_cov, levels=c('6.2% (Baseline)','15%','25%','50%'))

model_incidence_all %>% 
  filter(Year > 2020) %>% 
  ggplot(aes(x = Year, y = mean, col=as.factor(scr_group)))+
  geom_ribbon(mapping=aes(x = Year,ymin = lower, ymax = upper, fill=as.factor(scr_group)),
              inherit.aes = FALSE,
              alpha = 0.2
  ) +
  geom_line(size=1.3) +
  facet_wrap(~scr_cov+scr_group)+
  annotate("rect", xmin = 2018, xmax = 2040, ymin = inc_target_lower, ymax = inc_target_upper,
           alpha = .7,
           fill="lightgrey") +
  annotate("segment", x = 2018, xend = 2040, y = inc_target_mean, yend = inc_target_mean,
           colour = "grey30", linetype=2)+
  annotate("segment", x = 2030, xend = 2030, y = 0, yend = 10000,
           colour = "grey30", linetype=2)+
  annotate("text", x = 2025, y = inc_target_mean-500, family = "", fontface = 3, size=4, label = "2030 Elimination Target (95% CI)") +
  theme_bw(base_size = 14) +
  theme(
    strip.text.x = element_text(
      size = 12, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16),
    legend.position = "none"
  )+
  theme(plot.title = element_text(face = "italic", size=15))+
  coord_cartesian(xlim = c(2020,2040), ylim = c(0, 12000)) +
  labs(title = "Yearly HCV Incidence by Target Age Group and Screening Coverage", x="Year", y =("New Cases"))

# Supplementary figure 7 - Main result (mortality on WHO goals) for all screening groups and coverages ####

model_mortality_all <- dataframe_suppfig6 %>% filter(variable %in% c("Deaths"))
model_mortality_all <- model_mortality_all[-4]
model_mortality_all <- model_mortality_all[-5]
model_mortality_all <- pivot_wider(model_mortality_all, names_from = baseline_coverage)

model_mortality_all$scr_cov <- as.character(model_mortality_all$scr_cov)
model_mortality_all$scr_group <- as.character(model_mortality_all$scr_group)
model_mortality_all[model_mortality_all == "baseline"] <- "6.2% (Baseline)"
model_mortality_all[model_mortality_all == "0.15"] <- "15%"
model_mortality_all[model_mortality_all == "0.25"] <- "25%"
model_mortality_all[model_mortality_all == "0.5"] <- "50%"

model_mortality_all$scr_cov <- factor(model_mortality_all$scr_cov, levels=c('6.2% (Baseline)','15%','25%','50%'))

model_mortality_all %>% 
  filter(Year > 2020) %>% 
  ggplot(aes(x = Year, y = mean, col=as.factor(scr_group)))+
  geom_ribbon(mapping=aes(x = Year,ymin = lower, ymax = upper, fill=as.factor(scr_group)),
              inherit.aes = FALSE,
              alpha = 0.2
  ) +
  geom_line(size=1.3) +
  facet_wrap(~scr_cov+scr_group)+
  annotate("rect", xmin = 2018, xmax = 2040, ymin = mort_target_lower, ymax = mort_target_upper,
           alpha = .7,
           fill="lightgrey") +
  annotate("segment", x = 2018, xend = 2040, y = mort_target_mean, yend = mort_target_mean,
           colour = "grey30", linetype=2)+
  annotate("segment", x = 2030, xend = 2030, y = 0, yend = 10000,
           colour = "grey30", linetype=2)+
  annotate("text", x = 2025, y = mort_target_mean-500, family = "", fontface = 3, size=4, label = "2030 Elimination Target (95% CI)") +
  theme_bw(base_size = 14) +
  theme(
    strip.text.x = element_text(
      size = 12, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16),
    legend.position = "none"
  )+
  theme(plot.title = element_text(face = "italic", size=15))+
  coord_cartesian(xlim = c(2020,2040), ylim = c(0, 10000)) +
  labs(title = "Yearly HCV-Related Mortality by Target Age Group and Screening Coverage", x="Year", y =("HCV-Related Deaths"))



# Supplementary figure 8 - Yearly HCV mortality model vs. "data" ####

dataframes_suppfig8 <- list(df1,df17,df33,df49,df65,df81,df97,df113,df129,df145,df161,df177)
names(dataframes_suppfig8) <- list('df1','df17','df33','df49','df65','df81','df97','df113','df129','df145','df161','df177')

dataframe_suppfig8 <- list.rbind(combine_dataframes(dataframes_suppfig8))
mortality_suppfig8 <- dataframe_suppfig8 %>% filter(variable %in% c("Deaths"))
mortality_suppfig8 <- mortality_suppfig8[-6]

m1<- mortality_suppfig8 %>% filter(population_scenario %in% c("baseline"))
m2<- mortality_suppfig8 %>% filter(population_scenario %in% c("growth"))
m3<- mortality_suppfig8 %>% filter(population_scenario %in% c("plateau"))
m4<- mortality_suppfig8 %>% filter(population_scenario %in% c("decline"))
m1 <- pivot_wider(m1, names_from = baseline_coverage)
m2 <- pivot_wider(m2, names_from = baseline_coverage)
m3 <- pivot_wider(m3, names_from = baseline_coverage)
m4 <- pivot_wider(m4, names_from = baseline_coverage)
mort_suppfig8 <- rbind(m1,m2,m3,m4)

mort_suppfig8$population_scenario[mort_suppfig8$population_scenario == "baseline"] <- "Baseline"
mort_suppfig8$population_scenario[mort_suppfig8$population_scenario == "decline"] <- "Decline"
mort_suppfig8$population_scenario[mort_suppfig8$population_scenario == "growth"] <- "Growth"
mort_suppfig8$population_scenario[mort_suppfig8$population_scenario == "plateau"] <- "Plateau"

mort_suppfig8 %>% 
  filter(Year > 2008) %>% 
  ggplot(aes(x = Year)) +
  geom_ribbon(mapping=aes(x = Year,ymin = lower, ymax = upper, fill=population_scenario),
              inherit.aes = FALSE,
              alpha = 0.2
  ) +
  geom_line(aes(col=population_scenario, y = mean), size=1.3)+
  geom_point(data=hcv_deaths_data,
             aes(x=Year, y=value)) +
  
  geom_errorbar(hcv_deaths_data,
                mapping = aes(x = Year, ymin = lower_95, ymax= upper_95)) +
  facet_wrap(~population_scenario) +
  annotate("rect", xmin = 2010, xmax = 2040, ymin = mort_target_lower, ymax = mort_target_upper,
           alpha = .7,
           fill="lightgrey") +
  annotate("segment", x = 2010, xend = 2040, y = mort_target_mean, yend = mort_target_mean,
           colour = "grey30", linetype=2)+
  annotate("segment", x = 2030, xend = 2030, y = 0, yend = 16000,
           colour = "grey30", linetype=2)+
  annotate("text", x = 2025, y = mort_target_mean-500, family = "", fontface = 3, size=4, label = "2030 Elimination Target (95% CI)") +
  theme_bw(base_size = 14) +
  theme(
    strip.text.x = element_text(
      size = 12, color = "black", face = "italic"
    ),
    legend.text=element_text(size=16),
    legend.position = "None"
  )+
  theme(plot.title = element_text(face = "italic", size=15))+
  labs(title = "Yearly HCV-Related Mortality by Population Scenario", x="Year", y =("HCV-Related Deaths"))