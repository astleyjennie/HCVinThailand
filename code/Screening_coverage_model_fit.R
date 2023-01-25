# Model fit for baseline coverage

##############################
# HEP C AGE STRUCTURED MODEL #
##############################

# Cleanup - CAUTION clears R environment ####
rm(list=ls())

# Setup - load packages and define plotting functions ####

#setwd('C:/Users/jenni/Documents/MGH/thailand/dissertation/submission/code')

library(pacman)
p_load(deSolve, tidyverse, doParallel, manipulate, readxl, gridExtra, grid, scales, tictoc, Hmisc, viridis)

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

cov_mean <- 0.06 # 0.5%
cov_sd <- 0.02 # 0.25% standard dev
cov_lower <- cov_mean + 1.96*cov_sd
cov_upper <- cov_mean - 1.96*cov_sd

groupBASE <- rep(cov_mean,groups)
groupBASE_lower <- rep(cov_lower,groups)
groupBASE_upper <- rep(cov_upper,groups)

# Targeted screening strategies 

coverageBASE <- cov_mean
coverage1 <- 0.1
coverage2 <- 0.15
coverage3 <- 0.3
coverage4 <- 0.4
coverage5 <- 0.5

# Group A: over 30s
groupA1 <- c(rep(coverageBASE, 6), rep(coverage1, 15)) # 30+ at 10% per year
groupA2 <- c(rep(coverageBASE, 6), rep(coverage2, 15)) # 30+ at 20% per year
groupA3 <- c(rep(coverageBASE, 6), rep(coverage3, 15)) # 30+ at 30% per year
groupA4 <- c(rep(coverageBASE, 6), rep(coverage4, 15)) # 30+ at 40% per year
groupA5 <- c(rep(coverageBASE, 6), rep(coverage5, 15)) # 30+ at 50% per year

# Group B: over 40s
groupB1 <- c(rep(coverageBASE, 8), rep(coverage1, 13)) # 40+ at 10% per year
groupB2 <- c(rep(coverageBASE, 8), rep(coverage2, 13)) # 40+ at 20% per year
groupB3 <- c(rep(coverageBASE, 8), rep(coverage3, 13)) # 40+ at 30% per year
groupB4 <- c(rep(coverageBASE, 8), rep(coverage4, 13)) # 40+ at 40% per year
groupB5 <- c(rep(coverageBASE, 8), rep(coverage5, 13)) # 40+ at 50% per year

# Group C: over 50s
groupC1 <- c(rep(coverageBASE, 10), rep(coverage1, 11)) # 50+ at 10% per year
groupC2 <- c(rep(coverageBASE, 10), rep(coverage2, 11)) # 50+ at 20% per year
groupC3 <- c(rep(coverageBASE, 10), rep(coverage3, 11)) # 50+ at 30% per year
groupC4 <- c(rep(coverageBASE, 10), rep(coverage4, 11)) # 50+ at 40% per year
groupC5 <- c(rep(coverageBASE, 10), rep(coverage5, 11)) # 50+ at 50% per year

# Group D: over 60s
groupD1 <- c(rep(coverageBASE, 12), rep(coverage1, 9)) # 60+ at 10% per year
groupD2 <- c(rep(coverageBASE, 12), rep(coverage2, 9)) # 60+ at 20% per year
groupD3 <- c(rep(coverageBASE, 12), rep(coverage3, 9)) # 60+ at 30% per year
groupD4 <- c(rep(coverageBASE, 12), rep(coverage4, 9)) # 60+ at 40% per year
groupD5 <- c(rep(coverageBASE, 12), rep(coverage5, 9)) # 60+ at 50% per year


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

#set up initial death
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



# Birth rate function ####

# Birth rate multiplier (brm) with uncertainty

brm_mean <- 1.09
brm_sd <- 0


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
      
      # Total population of entire system (should fit Thai population data)
      total = population,
      
      # Total infections, HCC and HCV (sum across all age groups)
      infect = (F0+F1+F2+F3+C1+C2+C3+C4+HCCA+HCCB+HCCC+HCCD),
      totalHCC = (HCCA+HCCB+HCCC+HCCD),
      totalHCV  = (F0+F1+F2+F3+C1+C2+C3+C4),
      
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
      prev9_10 = 100*(infect9+infect10)/(group9+group10), #40-49
      prev11_21 = 100*(infect11+infect12+infect13+infect14+infect15+infect16+infect17+infect18+infect19+infect20+infect21)/(group11+group12+group13+group14+group15+group16+group17+group18+group19+group20+group21) # Over 50
      
    ) %>% 
    pivot_longer(names_to = "variable", cols = !1) 
  return(df1)
}


# Define model fit ####

getNLL <- function(params){
  
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
  birthrate_multiplier <- brm_mean
  birth.approx_mean$birth <- birthrate_multiplier * birth.approx_mean$birth
  birth.func <- approxfun(birth.approx_mean$t,birth.approx_mean$birth,method="linear")

  groupBASE <- rep(params[1],groups)
  scr_group <- groupBASE
  mortality_base1 <- mortality_base
  
  for(i in 19:37){
    mortality_base1[i,2:22] <- mortality_base1[i-1,2:22]*mort_scenario
  }
  
  mortality_approx <- mortality_base1
  
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
           
           if((t>scr_start)&(t<scr_start+scr_dur)){scr<-scr_group}else{scr<-groupBASE}
           
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
  


  
  
  
  out <- ode(method="rk4", y = init, times = simu.time, func=hepC.mod, parms = parameters, hini=0.1)
  out <- cbind(out[,1]+start_year,out)
  colnames(out)[1] <- "Year"
  out <- out[,colnames(out)!="time"]
  
  df1_base_scr_base <- mutate_data(out)
  
  model_prev_all  <- df1_base_scr_base %>% filter(variable %in% c("prev1_2","prev3_4","prev5_6","prev7_8","prev9_10","prev11_21", "infectprev"))
  model_prev_2014 <- model_prev_all %>% filter(Year %in% c("2014"))
  model_prev_2022 <- model_prev_all %>% filter(Year %in% c("2022"))
  
  data_prev_2014 <- data_prev %>% filter(Year %in% c("2014"))
  data_prev_2022 <- data_prev %>% filter(Year %in% c("2022")) %>% filter(variable %in% c("infectprev"))
  
  model_prev <- rbind(model_prev_2014,model_prev_2022)
  data_prev_both <- rbind(data_prev_2014,data_prev_2022)

  NLL <- -sum(dnorm(x=data_prev_both$value,
                  mean=model_prev$value,
                  sd=0.5,
                  #size=20,
                  log = TRUE))
  return(NLL)
}

#



# Get fit results ####

start.vector <- c(0.01) # baseline screening coverage

tic("modelfit")
fit1 <- optim(start.vector, getNLL, hessian=T, method="Brent", lower=0, upper=0.5)
toc()

fit1$hessian
fit_sd<-sqrt(diag(solve(fit1$hessian)))

fit1$par
fit_sd
5
