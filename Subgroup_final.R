############  Microsimulation Alzheimer's Disease in Czechia  ##########
#########             Hana M. Broulikova           #####################
# This code is based on the work of the DARTH group (darthworkgroup.com)
## Jalal H, et al. An Overview of R in Health Decision Sciences. 
# Med. Decis. Making. 2017; 37(3): 735-746. 
## Krijkamp EM, et al. Microsimulation modeling for health decision sciences 
# using R: a tutorial. Med. Decis. Making. 2018; 38(3): 400-422. 


rm(list =ls())                                              # clear memory (removes all the variables from the workspace)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the folder where the file is saved 

#### 01 Load packages ####
library(plyr) 
library(readr)
library(tidyverse)
library(truncnorm)
library(readxl)
#library(patchwork)

#### 02 Load Functions ####
# This is the package of functions provided during the coures.
source("functions/Functions.R")

#### 03 Input Model Parameters ####
set.seed(1)                                                 # set the seed 

#### Model structure 
n.t   <- 40                                                 # time horizon, 40 cycles to achieve lifetime horizon
n.i   <- 10000                                             # number of simulated individuals
v.n   <- c("Mild","Moderate","Severe","Death")              # the model states names
n.s   <- length(v.n)                                        # the number of health states

d.r   <- 0.03                                               # discount rate of 3% per cycle
v.dwe <- v.dwc <- 1 / ((1 + d.r) ^ (0:n.t))                 # discount weight 


group   <- c("Mean", "Lopez")                               # decline schemes considered in this model
mmse    <- 28                                               # intital MMSE at the start of the model
#tmmse   <-  runif(n.i,11, 21)                               # MMSE of diagnosis and treatment initiation (SA tmmse19, SAimmse 20)
#immse   <-  runif(n.i,21, 28)   
#distr   <- read.csv("sampling.csv") %>% 
#  janitor::clean_names() %>% 
#  mutate(age = onset)

# load age at disease onset distribution
p.Treatment <- 0.5                                         # indicates the probability of receiving treatment without intervention
p.Treatment_i <- 1
p.Slow  <- 0.39                                             # probability of being slow progressor regardless treatment 
p.Fs    <- 0.6                                              # probability of switching from fast to slow progression with treatment
p.Sf    <- 0                                                # probability of switching from slow to fast with treatment
hr      <- 2.1                                              # hazard ratio for AD patients


#### Annual probability of death (Czech mortality tables increased by HR for AD)
p.Death <- read.csv("ut.csv") %>%                         # load Czech mortality tables and adjust by HR
  mutate(Male = 1 - exp(-(- log(1 - qxM))*hr), Female = 1 - exp(-(- log(1 - qxZ))*hr)) %>% 
  select(age, Male, Female) %>% 
  gather(key = "sex" , value = "p.Mortality", 2:3) 


#### Cost inputs (unit costs discounted to the value of 2017)
c.Mild.i      <- 951*12*((1+0.03)^3)                        # annual costs of indirect care for mild cognitive deficit
c.Moderate.i  <- 1872*12*((1+0.03)^3)                       # annual costs of indirect care for moderate cognitive deficit
c.Severe.i    <- 2155*12*((1+0.03)^3)                       # annual costs of indirect care for severe cognitive deficit
c.Outpatient  <- 97*((1+0.03)^3)                            # annual costs of outpatient care (for treated patients)
c.Mild.d      <- 93*2*((1+0.03)^5)                          # annual costs of drugs for mild cognitive deficit (for treated patients)
c.Moderate.d  <- 93*2*((1+0.03)^5)                          # annual costs of drugs for moderate cognitive deficit (for treated patients)
c.Severe.d    <- 289*2*((1+0.03)^5)                         # annual costs of drugs for severe cognitive deficit (for treated patients)
c.Screening   <- (97*((1+0.03)^3) + (26 + 18)/((1+0.03)^2))*31 # costs of screening (for treated patients in the intervention branch)
c.Death       <- 0                                          # annual costs for dead patients 

#### Static characteristics of the cohort
v.MMSE<- rep(mmse, n.i)                                     # a vector with the intitial mmse score 
#rg    <- runif(n.i)                                        # random number for gender sampling
#v.Age <- sample(x = distr$onset, prob = distr$freq_age,  size = n.i, replace = TRUE)          
# sample from age distribution an initial age for every individual
#df.X  <- data.frame(ID = 1:n.i, age = v.Age, rg = rg, MMSE = v.MMSE) # binding patient characteristics into the dataframe
#df.X  <- merge(df.X, distr, by = "age") %>% 
#  mutate(sex = ifelse(rg<freq_gen, "Male", "Female")) %>% 
#  select(age, sex, MMSE)                                    # merging probablity of being man to the dataframe, creating gender


#### 04 Define Simulation Functions ####
#### Probability function
# The Probs function that updates the transition probabilities of every cycle 
Probs <- function(M_it, df.X, MMSE, t) { 
  # Arguments:
  # M_it: health state occupied by individual i at cycle t (character variable)
  # df.X: dataframe used to produce mortality na mmse vectors 
  # MMSE: score from Mini mental state examination usdd to measure cognitive performace in state t
  # t:    current cycle 
  # Returns: 
  # transition probabilities for that cycle
  
  m.p.it          <- matrix(0, nrow = n.s, ncol = n.i)      # create matrix of state transition probabilities
  rownames(m.p.it) <-  v.n                                  # give the state names to the rows
  
  # lookup baseline probability of dying and MMSE score based on individual characteristics
  population <- join(df.X, p.Death, by = c("age", "sex")) 
  
  p.Mild.d          <- population[M_it == "Mild","p.Mortality"]
  p.Moderate.d      <- population[M_it == "Moderate","p.Mortality"]
  p.Severe.d        <- population[M_it == "Severe","p.Mortality"]
  
  p.Mild.mmse      <- population[M_it == "Mild","MMSE"]
  p.Moderate.mmse  <- population[M_it == "Moderate","MMSE"]
  p.Severe.mmse    <- population[M_it == "Severe","MMSE"]
  
  # update the v.p with the appropriate probabilities (order of probabilities: Mild, Moderate, Severe, Death) 
  m.p.it[,M_it == "Mild"]     <- rbind(ifelse(ifelse(p.Mild.mmse >= 21, 1, 0) >0, 1-p.Mild.d, 0),
                                       ifelse(ifelse(p.Mild.mmse < 21 & p.Mild.mmse >= 11, 1, 0) >0, 1-p.Mild.d, 0), 
                                       ifelse(ifelse(p.Mild.mmse < 11, 1, 0) >0, 1-p.Mild.d, 0),
                                       p.Mild.d)
  m.p.it[,M_it == "Moderate"] <- rbind(ifelse(ifelse(p.Moderate.mmse >= 21, 1, 0) >0, 1-p.Moderate.d, 0),
                                       ifelse(ifelse(p.Moderate.mmse < 21 & p.Moderate.mmse >= 11, 1, 0) >0, 1-p.Moderate.d, 0), 
                                       ifelse(ifelse(p.Moderate.mmse < 11, 1, 0) >0, 1-p.Moderate.d, 0),
                                       p.Moderate.d)
  m.p.it[,M_it == "Severe"]   <- rbind(ifelse(ifelse(p.Severe.mmse >= 21, 1, 0) >0, 1-p.Severe.d, 0),
                                       ifelse(ifelse(p.Severe.mmse < 21 & p.Severe.mmse >= 11, 1, 0) >0, 1-p.Severe.d, 0),
                                       ifelse(ifelse(p.Severe.mmse < 11, 1, 0) >0, 1-p.Severe.d, 0),
                                       p.Severe.d)
  m.p.it[,M_it == "Death"]      <- rbind(0, 0, 0, 1)
  
  
  
  
  return(t(m.p.it))        		                                 # return the transition probabilities
}       
#### Costs function
# The Costs function estimates the costs at every cycle.
Costs <- function (M_it, Streatment, Int= FALSE) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Streatment: has the treatment been initiated? 
  # Int: is the individual in the intervention branch? (default is FALSE)
  
  c.it <- 0                                                                                                     # by default the cost for everyone is zero 
  c.it[M_it == "Mild" & Streatment == 0]      <- c.Mild.i                                                       # update the cost if mild untreated
  c.it[M_it == "Mild" & Streatment == 1]      <- c.Mild.i + c.Mild.d + c.Outpatient + c.Screening*Int           # update the cost if mild detected
  c.it[M_it == "Mild" & Streatment > 1]       <- c.Mild.i + c.Mild.d + c.Outpatient                             # update the cost if mild treated
  c.it[M_it == "Moderate" & Streatment == 0]  <- c.Moderate.i                                                   # update the cost if moderate untreated
  c.it[M_it == "Moderate" & Streatment == 1]  <- c.Moderate.i + c.Moderate.d + c.Outpatient + c.Screening*Int   # update the cost if moderate detected
  c.it[M_it == "Moderate" & Streatment > 1]   <- c.Moderate.i + c.Moderate.d + c.Outpatient                     # update the cost if moderate treated
  c.it[M_it == "Severe" & Streatment == 0]    <- c.Severe.i                                                     # update the cost if severe untreated
  c.it[M_it == "Severe" & Streatment == 1]    <- c.Severe.i + c.Severe.d + c.Outpatient + c.Screening*Int       # update the cost if severe detected
  c.it[M_it == "Severe" & Streatment > 1]     <- c.Severe.i + c.Severe.d + c.Outpatient                         # update the cost if severe treated
  c.it[M_it == "Death"]  <- c.Death                                                                             # update the cost if dead
  
  return(c.it)        		                                                                                      # return the costs
}

#### Partial cost functions
# Informal care
Costs_i <- function (M_it) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Streatment: has the treatment been initiated? 
  # Int: is the individual in the intervention branch? (default is FALSE)
  
  c_inf <- 0                                                                # by default the cost for everyone is zero 
  c_inf[M_it == "Mild"]      <- c.Mild.i                                    # update the cost if mild 
  c_inf[M_it == "Moderate"]  <- c.Moderate.i                                # update the cost if mild
  c_inf[M_it == "Severe"]    <- c.Severe.i                                  # update the cost if severe 
  c_inf[M_it == "Death"]     <- c.Death                                     # update the cost if dead
  
  return(c_inf)        		                                                  # return the costs
}

# Outpatient care
Costs_o <- function (M_it, Streatment) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Streatment: has the treatment been initiated? 
  # Int: is the individual in the intervention branch? (default is FALSE)
  
  c_out <- 0                                                                # by default the cost for everyone is zero 
  c_out[M_it == "Mild" & Streatment == 0]      <- 0                         # update the cost if mild untreated
  c_out[M_it == "Mild" & Streatment >= 1]      <- c.Outpatient              # update the cost if mild detected
  c_out[M_it == "Moderate" & Streatment == 0]  <- 0                         # update the cost if moderate untreated
  c_out[M_it == "Moderate" & Streatment >= 1]  <- c.Outpatient              # update the cost if moderate detected
  c_out[M_it == "Severe" & Streatment == 0]    <- 0                         # update the cost if severe untreated
  c_out[M_it == "Severe" & Streatment >= 1]    <- c.Outpatient              # update the cost if severe detected
  c_out[M_it == "Death"]                       <- c.Death                   # update the cost if dead
  
  
  return(c_out)        		                                                  # return the costs
}

# Medication
Costs_m <- function (M_it, Streatment) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Streatment: has the treatment been initiated? 
  # Int: is the individual in the intervention branch? (default is FALSE)
  
  c_med <- 0                                                                # by default the cost for everyone is zero 
  c_med[M_it == "Mild" & Streatment == 0]      <- 0                         # update the cost if mild untreated
  c_med[M_it == "Mild" & Streatment >= 1]      <- c.Mild.d                  # update the cost if mild treated
  c_med[M_it == "Moderate" & Streatment == 0]  <- 0                         # update the cost if moderate untreated
  c_med[M_it == "Moderate" & Streatment >= 1]  <- c.Moderate.d              # update the cost if moderate treated
  c_med[M_it == "Severe" & Streatment == 0]    <- 0                         # update the cost if severe untreated
  c_med[M_it == "Severe" & Streatment >= 1]    <- c.Severe.d                # update the cost if severe treated
  c_med[M_it == "Death"]                       <- c.Death                   # update the cost if dead
  
  return(c_med)        		                                                  # return the costs
}

# Screening
Costs_s <- function (M_it, Streatment, Int= FALSE) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Streatment: has the treatment been initiated? 
  # Int: is the individual in the intervention branch? (default is FALSE)
  
  c_s <- 0                                                                  # by default the cost for everyone is zero 
  c_s[M_it == "Mild" & Streatment == 0]      <- 0                           # update the cost if mild untreated
  c_s[M_it == "Mild" & Streatment == 1]      <- c.Screening*Int             # update the cost if mild detected
  c_s[M_it == "Mild" & Streatment > 1]       <- 0                           # update the cost if mild treated
  c_s[M_it == "Moderate" & Streatment == 0]  <- 0                           # update the cost if moderate untreated
  c_s[M_it == "Moderate" & Streatment == 1]  <- c.Screening*Int             # update the cost if moderate detected
  c_s[M_it == "Moderate" & Streatment > 1]   <- 0                           # update the cost if moderate treated
  c_s[M_it == "Severe" & Streatment == 0]    <- 0                           # update the cost if severe untreated
  c_s[M_it == "Severe" & Streatment == 1]    <- c.Screening*Int             # update the cost if severe detected
  c_s[M_it == "Severe" & Streatment > 1]     <- 0                           # update the cost if severe treated
  c_s[M_it == "Death"]                       <- c.Death                     # update the cost if dead
  
  return(c_s)                                                               # return the costs     	
}

# Fraction treated
# Screening
Fraction <- function (M_it, Streatment) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Streatment: has the treatment been initiated? 
  # Int: is the individual in the intervention branch? (default is FALSE)
  
  fr <- 0                                                                  # by default the cost for everyone is zero 
  fr[M_it == "Mild" & Streatment == 0]      <- 0                           # update the cost if mild untreated
  fr[M_it == "Mild" & Streatment == 1]      <- 1                           # update the cost if mild detected
  fr[M_it == "Mild" & Streatment > 1]       <- 0                           # update the cost if mild treated
  fr[M_it == "Moderate" & Streatment == 0]  <- 0                           # update the cost if moderate untreated
  fr[M_it == "Moderate" & Streatment == 1]  <- 1                           # update the cost if moderate detected
  fr[M_it == "Moderate" & Streatment > 1]   <- 0                           # update the cost if moderate treated
  fr[M_it == "Severe" & Streatment == 0]    <- 0                           # update the cost if severe untreated
  fr[M_it == "Severe" & Streatment == 1]    <- 1                           # update the cost if severe detected
  fr[M_it == "Severe" & Streatment > 1]     <- 0                           # update the cost if severe treated
  fr[M_it == "Death"]                       <- c.Death                     # update the cost if dead
  
  return(fr)                                                               # return the costs     	
}


#### Function updating MMSE
# The upMMSE function updates decline scheme, calculates annual decline and substracts annual decline from cognitive score in previous cycle
upMMSE = function(MMSE, v.Scheme, v.Switch, v.Streatment){   
  # MMSE: cognitive score in previous cycle
  # vStreatment: has the treatment been initiated? 
  # Int: is the individual in the intervention branch? (default is FALSE) 
  
  # assign decline scheme and treatment availability
  decline <-  ifelse(v.Scheme == "Mean" & v.Streatment > 0, "treated_mean",
                     ifelse(v.Scheme == "Mean" & v.Streatment == 0, "untreated_mean",
                            ifelse(v.Scheme == "Lopez.s" & v.Streatment > 0, "treated_lopez_slow",
                                   ifelse (v.Scheme == "Lopez.s" & v.Streatment == 0, "untreated_lopez_slow",
                                           ifelse(v.Scheme == "Lopez.f" & v.Streatment > 0 & v.Switch == 0, "treated_lopez_fast",
                                                  ifelse(v.Scheme == "Lopez.f" & v.Streatment > 0 & v.Switch == 1, "treated_lopez_slow",
                                                         "untreated_lopez_fast"))))))
  
  # calculate number of patients declining according to each scheme
  tm    <- sum(decline == "treated_mean")              
  utm   <- sum(decline == "untreated_mean") 
  tls   <- sum(decline == "treated_lopez_slow")
  utls  <- sum(decline == "untreated_lopez_slow")
  tlf   <- sum(decline == "treated_lopez_fast")
  utlf  <- sum(decline == "untreated_lopez_fast")
  
  # draw annual cognitive decline 
  d.MMSE = MMSE
  d.MMSE[decline == "treated_mean"] = rtruncnorm(tm, a = 0, b = Inf, mean = 1.5, sd = 1.5)
  d.MMSE[decline == "untreated_mean"] = rtruncnorm(utm, a = 0, b = Inf, mean = 3.5, sd = 1.5)
  d.MMSE[decline == "treated_lopez_slow"] = runif(tls, min = -1, max = 2)
  d.MMSE[decline == "untreated_lopez_slow"] = runif(utls, min = -1, max = 2)
  d.MMSE[decline == "treated_lopez_fast"] = runif(tlf, min = 3, max = 5)
  d.MMSE[decline == "untreated_lopez_fast"] = runif(utlf, min = 3, max = 6.8) 
  
  MMSE <- ifelse((MMSE - d.MMSE) > 0, (MMSE - d.MMSE), 0)                                     # calculate current MMSE score
  return(MMSE)                                                                                # return the MMSE
}
#### 05 Microsimulation function ####
MicroSim <- function(n.i, df.X , Int = FALSE, seed = 1) {
  # Arguments:  
  # n.i:     number of individuals
  # df.X     data frame with individual data 
  # Age      age of the individuals
  # Sex      sex of the indivuduals 
  # Int:     is this the individual receiving intervention? (default is FALSE)
  # seed:    defauls is 1
  
  set.seed(seed)                                                                                  # set the seed
  n.s <- length(v.n)                                                                              # the number of health states
  
  #Static characteristics (additional to age and gender)  
  ra <- runif(n.i)                                                                                # vector of random numbers used for the decision about treatment availability
  rb <- runif(n.i)                                                                                # vector of random numbers used for the decision about initial progression pace
  rc <- runif(n.i)                                                                                # vector of random numbers used for switching from fast to slow disease progression if treated
  v.Group       <- sample(x = group, prob = c(0.5, 0.5), size = n.i, replace = TRUE)              # randomly assigned decline scheme
  v.Treatment   <- ifelse(ra<= ifelse(Int == 0, p.Treatment, p.Treatment_i), 1, 0)                # vector indicating treatment availability for a patient
  v.Scheme      <- ifelse(v.Group == "Mean", "Mean",                                              # vecotr indicating slow or fast disease progression for Lopez
                          ifelse(rb<=p.Slow, "Lopez.s", "Lopez.f"))
  v.Switch      <- ifelse(v.Group == "Mean", 1,                                                   # vecotr indicating switching from fast to slow progression if treated 
                          ifelse(v.Treatment == 1 & rc<=p.Fs, 1, 0))
  
  # Dynamic characteristics 
  v.M_Init            <- rep("Mild", n.i)                                                         # a vector with the initial health state for all individuals 
  v.Streatment_Init   <- rep(0, n.i)                                                              # a vector with the initial treatment status at the start of the model 
  
  # create three matrices called m.M, m.C, m.E, m.Streatment, m.MMSE
  # number of rows is equal to the n.i, the number of columns is equal to n.t  (the initial state and all the n.t cycles)
  # m.M is used to store the health state information over time for every individual
  # m.C is used to store the costs information over time for evey individual
  # m.E is used to store the effects information over time for every individual
  # m.Streatment is used to track time since treatment initiation 
  # m.MMSE is used to store cognitive scores over time for every individual
  
  m.F <- m.c.informal <- m.c.outpatient <- m.c.screening <- m.c.medication <- m.M <- m.C <- m.E <- m.Streatment <- m.MMSE <- matrix(nrow = n.i, ncol = n.t + 1, 
                                                                                                                                    dimnames = list(paste("ind"  , 1:n.i, sep = " "), 
                                                                                                                                                    paste("cycle", 0:n.t, sep = " ")))  
  
  m.M [, 1]         <- v.M_Init                                                                   # initial health state at cycle 0 for individual i
  m.Streatment[,1]  <- v.Streatment_Init                                                          # initialize time since treatment initition for individual i
  m.MMSE[, 1]       <- v.MMSE                                                                     # initial mmse at cycle 0 for individual i
  m.C[, 1]          <- Costs(m.M[, 1], m.Streatment[,1], Int)                                     # calculate costs per individual during cycle 0
  #m.E[, 1]         <- Effs (m.M[, 1], df.X)                                                      # calculate QALYs per individual during cycle 0
  m.E[, 1]          <- 0.408+0.01*df.X$MMSE
  m.c.informal[, 1] <- Costs_i(m.M[, 1])                                                          # calculate informal costs per individual during cycle 0
  m.c.outpatient[, 1]<- Costs_o(m.M[, 1], m.Streatment[,1])                                       # calculate outpatient costs per individual during cycle 0
  m.c.screening[, 1]<- Costs_s(m.M[, 1], m.Streatment[,1], Int)                                   # calculate screening costs per individual during cycle 0
  m.c.medication[, 1]<- Costs_m(m.M[, 1], m.Streatment[,1])                                       # calculate medication costs per individual during cycle 0
  m.F[, 1]<- Fraction(m.M[, 1], m.Streatment[,1])                                                 # calculate outpatient costs per individual during cycle 0
  
  # open a loop for time running cycles 1 to n.t 
  for (t in 1:n.t) { 
    
    v.p             <- Probs(m.M[, t], df.X, df.X$MMSE, t)                                        # calculate the transition probabilities for the cycle based on health state t
    m.M[, t + 1]    <- samplev(v.p, 1)                                                            # sample the current health state and store that state in matrix m.M 
    df.X$MMSE       <- upMMSE(df.X$MMSE, v.Scheme, v.Switch, m.Streatment[,t])                    # calculate mmse score during the cycle t + 1
    m.MMSE[, t + 1] <-  df.X$MMSE                                                                 # store MMSE score to the matrix
    m.C[, t + 1]    <- Costs(m.M[, t + 1], m.Streatment[,t], Int)                                 # calculate costs per individual during cycle t + 1
    m.c.informal[, t + 1]    <- Costs_i(m.M[, t + 1])                                             # calculate informal costs per individual during cycle t + 1
    m.c.outpatient[, t + 1]  <- Costs_o(m.M[, t + 1], m.Streatment[,t])                           # calculate outpatient costs per individual during cycle t + 1
    m.c.screening[, t + 1]   <- Costs_s(m.M[, t + 1], m.Streatment[,t], Int)                      # calculate screening costs per individual during cycle t + 1
    m.c.medication[, t + 1]  <- Costs_m(m.M[, t + 1], m.Streatment[,t])                           # calculate medication costs per individual during cycle t + 1
    m.F[, t + 1]  <- Fraction(m.M[, t + 1], m.Streatment[,t])                                     # calculate medication costs per individual during cycle t + 1
    m.E[, t + 1]    <- ifelse(m.M[, t + 1] != "Death", 0.408+0.01*df.X$MMSE,0)                    # calculate QALYs per individual during cycle t + 1
    # update time of treatment 
    m.Streatment[, t + 1]            <- ifelse((Int == 0 & v.Treatment==1 & df.X$MMSE <= tmmse & m.M[, t + 1] != "Death") 
                                               |(Int == 1 & v.Treatment==1 & df.X$MMSE <= immse & m.M[, t + 1] != "Death"), m.Streatment[,t] + 1,
                                               ifelse(m.Streatment[,t] > 0 & m.M[, t + 1] != "Death", m.Streatment[,t] + 1,
                                                      m.Streatment[,t]))
    # update age
    df.X$age[m.M[,t + 1] != "Death"] <- df.X$age[m.M[, t + 1] != "Death"] + 1
  } # close the loop for the time points
  
    # display simulation progress
  if(t/(n.t/10) == round(t/(n.t/10), 0)) {                                                        # display progress every 10%
    cat('\r', paste(t/n.t * 100, "% done", sep = " "))
  }
  
  tc <- m.C %*% v.dwc
  te <- m.E %*% v.dwe
  ti <- m.c.informal %*% v.dwc
  to <- m.c.outpatient %*% v.dwc
  ts <- m.c.screening %*% v.dwc
  tm <- m.c.medication %*% v.dwc
  
  tc_hat <- mean(tc)
  te_hat <- mean(te)
  results <<- list(m.Streatment = m.Streatment, m.MMSE = m.MMSE, m.M = m.M, m.C = m.C, m.E = m.E, tc = tc , te = te, tc_hat = tc_hat, te_hat = te_hat, m.c.informal = m.c.informal, m.c.outpatient = m.c.outpatient, m.c.screening = m.c.screening, m.c.medication = m.c.medication, ti = ti, to = to, ts = ts, tm = tm, m.F = m.F) # store the results from the simulation in a list  
  return(results)                                                                                 # return the results
}                                     

#### Extensisons

age_for = seq(from = 70, 
              to = 90, 
              by = 10)
sex_for = c("Male", "Female")
tmmse_for   <-  seq(from = 16, 
                    to = 16,
                    by = 1)
#tmmse = 16
immse_for   <-  seq(from = 16, 
                    to = 28,
                    by = 1) 


Results_w = matrix(NA,
                   nrow = length(age_for) * length(sex_for) * length(tmmse_for) * length(immse_for),
                   ncol = 20,
                   dimnames = list(NULL,
                                   c("age", "gender", "immse", "tmmse", 
                                     "Informal_no_int", "Informal_int", 
                                     "Medication_no_int","Medication_int", 
                                     "Outpatient_no_int", "Outpatient_int",
                                     "Screening_no_int", "Screening_int", 
                                     "Total_no_int","Total_int",
                                     "QALY_no_int" , "QALY_int", 
                                     "Treated_no_int", "Treated_int",
                                     "NB_no_int", "NB_int"))) # případně dodat  dobu dožití
# ---- For cycles ----
i = 1
#for (age in age_for) {
# for (sex in sex_for) { 
for (tmmse in tmmse_for) {
  for (immse in immse_for) {
    for (age in age_for) {
      for (sex in sex_for) {
        
        df.X  <- data.frame(ID = 1:n.i, age = age, sex = sex, MMSE = v.MMSE) %>% 
          select(age, sex, MMSE) 
        
        # ---- Simulation Early Assesment with/without Treatment ----
        outcomes_no_int  <- MicroSim(n.i, df.X, Int = FALSE, seed = 1)
        outcomes_int     <- MicroSim(n.i, df.X, Int = TRUE, seed = 1)
      }
      
      # save results into the matrix
      Res = c(age, sex, immse, tmmse,
              sum(outcomes_no_int[["ti"]])/n.i, sum(outcomes_int[["ti"]])/n.i, #informal
              sum(outcomes_no_int[["tm"]])/n.i, sum(outcomes_int[["tm"]])/n.i, #medication
              sum(outcomes_no_int[["to"]])/n.i, sum(outcomes_int[["to"]])/n.i, #outpatient
              sum(outcomes_no_int[["ts"]])/n.i, sum(outcomes_int[["ts"]])/n.i, #screening 
              NA, NA, #Total costs
              sum(outcomes_no_int[["te"]])/n.i, sum(outcomes_int[["te"]])/n.i, #qalys
              sum(outcomes_no_int[["m.F"]])/n.i, sum(outcomes_int[["m.F"]])/n.i, # fraction of treated
              NA, NA) # Net benefit
      
      Results_w[i,] = Res
      
      i = i + 1
      print(i)
    }
  }
}

age_for = seq(from = 70, 
              to = 90, 
              by = 10)
sex_for = c("Female", "Male")
tmmse_for   <-  seq(from = 16, 
                    to = 16,
                    by = 1)
#tmmse = 16
immse_for   <-  seq(from = 16, 
                    to = 28,
                    by = 1) 


Results_m = matrix(NA,
                   nrow = length(age_for) * length(sex_for) * length(tmmse_for) * length(immse_for),
                   ncol = 20,
                   dimnames = list(NULL,
                                   c("age", "gender", "immse", "tmmse", 
                                     "Informal_no_int", "Informal_int", 
                                     "Medication_no_int","Medication_int", 
                                     "Outpatient_no_int", "Outpatient_int",
                                     "Screening_no_int", "Screening_int", 
                                     "Total_no_int","Total_int",
                                     "QALY_no_int" , "QALY_int", 
                                     "Treated_no_int", "Treated_int",
                                     "NB_no_int", "NB_int"))) # případně dodat  dobu dožití
# ---- For cykly ----
i = 1
#for (age in age_for) {
# for (sex in sex_for) { 
for (tmmse in tmmse_for) {
  for (immse in immse_for) {
    for (age in age_for) {
      for (sex in sex_for) {
        
        df.X  <- data.frame(ID = 1:n.i, age = age, sex = sex, MMSE = v.MMSE) %>% 
          select(age, sex, MMSE) 
        
        # ---- Simulation Early Assesment with/without Treatment ----
        outcomes_no_int  <- MicroSim(n.i, df.X, Int = FALSE, seed = 1)
        outcomes_int     <- MicroSim(n.i, df.X, Int = TRUE, seed = 1)
      }
      
      #zapsání výsledků do matice
      Res = c(age, sex, immse, tmmse,
              sum(outcomes_no_int[["ti"]])/n.i, sum(outcomes_int[["ti"]])/n.i, #informal
              sum(outcomes_no_int[["tm"]])/n.i, sum(outcomes_int[["tm"]])/n.i, #medication
              sum(outcomes_no_int[["to"]])/n.i, sum(outcomes_int[["to"]])/n.i, #outpatient
              sum(outcomes_no_int[["ts"]])/n.i, sum(outcomes_int[["ts"]])/n.i, #screening 
              NA, NA, #Total costs
              sum(outcomes_no_int[["te"]])/n.i, sum(outcomes_int[["te"]])/n.i, #qalys
              sum(outcomes_no_int[["m.F"]])/n.i, sum(outcomes_int[["m.F"]])/n.i, # fraction of treated
              NA, NA) # Net benefit
      
      Results_m[i,] = Res
      
      i = i + 1
      print(i)
    }
  }
}





Results <- rbind(Results_m, Results_w) %>% 
  as.data.frame() 
#write.csv2(Results, "Results_nb.csv")

#### Graphical analysis
library(readxl)
Results_plots <- read_excel("Results_15_9.xlsx", 
                        col_types = c("skip", "numeric", "text", 
                        "numeric", "numeric", "numeric", 
                        "numeric", "numeric", "numeric", 
                        "numeric", "numeric", "numeric", 
                        "numeric", "numeric", "numeric", 
                        "numeric", "numeric", "numeric", 
                        "numeric", "numeric", "numeric"))

Results_plots <-   mutate(Results_plots, Total_no_int = Informal_no_int + Medication_no_int + Outpatient_no_int + Screening_no_int, 
                    Total_int = Informal_int + Medication_int + Outpatient_int + Screening_int,
                    k = 45000,
                    NB_no_int = k*QALY_no_int - Total_no_int,
                    NB_int = k*QALY_int - Total_int,
                    NB = NB_int -NB_no_int) 


Results_plots <- mutate(Results_plots, Age = age, 
                        Gender = ifelse(gender == "Male", 1, 0), 
                        MMSE = immse, 
                        Drugs = (Medication_int - Medication_no_int)*(-1),
                        Home = (Informal_int - Informal_no_int)*(-1),
                        Outpatient = (Outpatient_int - Outpatient_no_int)*(-1),
                        Healthcare = ((Medication_int + Outpatient_int + Screening_int) - (Medication_no_int + Outpatient_no_int + Screening_no_int))*(-1),
                        Total = (Total_int - Total_no_int)*(-1)) %>% select(Gender, Age, MMSE, Home, Healthcare, Total)

results_long = gather(Results_plots, key = Category, value = Cost, Home:Total) %>% 
  mutate(Category_gender = 
           ifelse (Category == "Home" & Gender == 1, "Home_men",
                   ifelse(Category == "Healthcare" & Gender == 1, "Healthcare_men", 
                          ifelse(Category == "Total" & Gender == 1, "Total_men",
                                 ifelse (Category == "Home" & Gender == 0, "Home_women",
                                         ifelse(Category == "Healthcare" & Gender == 0, "Healthcare_women", 
                                                ifelse(Category == "Total" & Gender == 0, "Total_women", "missing"
                                                )))))))
euro_format <- function(largest_with_cents = 100000) {
  function(x) {
    x <- round_any(x, 0.01)
    if (max(x, na.rm = TRUE) < largest_with_cents &
        !all(x == floor(x), na.rm = TRUE)) {
      nsmall <- 2L
    } else {
      x <- round_any(x, 1)
      nsmall <- 0L
    }
    str_c("€", format(x, nsmall = nsmall, trim = TRUE, big.mark = ",", scientific = FALSE, digits=1L))
  }
}

plot_70 <- ggplot(subset(results_long, Age==70), aes(x = MMSE, y = Cost, color = Category_gender,  linetype = Category)) + 
  geom_line() + 
  scale_color_manual(values=c("blue", "red", "blue", "red", "blue", "red")) +
  scale_linetype_manual(values=c("dashed", "dotted", "solid"), name = "Payer Perspective", labels=c("Insurer", "Patient", "Society")) +
  labs(title="Potential Savings for a Patient Contracting AD at the Age of 70") +
  xlab("MMSE at Treatment Initiation") +
  ylab("Potential Savings") +
  scale_y_continuous(labels = euro_format()) +
  scale_x_reverse(breaks = c(28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)) +
  guides(colour=FALSE) +
  theme(legend.position="none")
plot_70

plot_80 <- ggplot(subset(results_long, Age==80), aes(x = MMSE, y = Cost, color = Category_gender,  linetype = Category)) + 
  geom_line() + 
  scale_color_manual(values=c("blue", "red", "blue", "red", "blue", "red")) +
  scale_linetype_manual(values=c("dashed", "dotted", "solid"), name = "Payer Perspective", labels=c("Insurer", "Family", "Society")) +
  labs(title="Potential Savings for a Patient Contracting AD at the Age of 80") +
  xlab("MMSE at Treatment Initiation") +
  ylab("Potential Savings") +
  scale_y_continuous(labels = euro_format()) +
  scale_x_reverse(breaks = c(28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)) +
  guides(colour=FALSE) +
  theme(legend.position="none")
plot_80

plot_90 <- ggplot(subset(results_long, Age==90), aes(x = MMSE, y = Cost, color = Category_gender,  linetype = Category)) + 
  geom_line() + 
  scale_color_manual(values=c("blue", "red", "blue", "red", "blue", "red")) +
  scale_linetype_manual(values=c("dashed", "dotted", "solid"), name = "Payer Perspective", labels=c("Insurer", "Patient", "Society")) +
  labs(title="Potential Savings for a Patient Contracting AD at the Age of 90") +
  xlab("MMSE at Treatment Initiation") +
  ylab("Potential Savings") +
  scale_y_continuous(labels = euro_format()) +
  scale_x_reverse(breaks = c(28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16)) +
  guides(colour=FALSE) +
  theme(legend.position="none")
plot_90














