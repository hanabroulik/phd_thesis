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

#### 02 Load Functions ####
# These functions are a product of the DARTH group. 
source("functions/Functions.R") 

#### 03 Input Model Parameters ####
set.seed(1)                                                 # set the seed  

#### Model structure 
n.t   <- 40                                                 # time horizon, 40 cycles to achieve lifetime horizon
n.i   <- 100000                                             # number of simulated individuals
v.n   <- c("Mild","Moderate","Severe","Death")              # the model states names
n.s   <- length(v.n)                                        # the number of health states

d.r   <- 0.03                                               # discount rate of 3% per cycle
v.dwe <- v.dwc <- 1 / ((1 + d.r) ^ (0:n.t))                 # discount weight 


group   <- c("Mean", "Lopez")                               # decline schemes considered in this model
mmse    <- 28                                               # intital MMSE at the start of the model
tmmse   <-  runif(n.i,11, 21)                               # MMSE of diagnosis and treatment initiation (SA tmmse19, SAimmse 20)
immse   <-  runif(n.i,21, 28)   
distr   <- read.csv("sampling.csv") %>% 
  janitor::clean_names() %>% 
  mutate(age = onset)

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
rg    <- runif(n.i)                                         # random number for gender sampling
v.Age <- sample(x = distr$onset, prob = distr$freq_age,  size = n.i, replace = TRUE)          
# sample from age distribution an initial age for every individual
df.X  <- data.frame(ID = 1:n.i, age = v.Age, rg = rg, MMSE = v.MMSE) # binding patient characteristics into the dataframe
df.X  <- merge(df.X, distr, by = "age") %>% 
  mutate(sex = ifelse(rg<freq_gen, "Male", "Female")) %>% 
  select(age, sex, MMSE)                                    # merging probablity of being man to the dataframe, creating gender


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
  return(MMSE)                                                                                                # return the MMSE
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
  
  # Static characteristics (additional to age and gender)  
  ra <- runif(n.i)                                                                                # vector of random numbers used for the decision about treatment availability
  rb <- runif(n.i)                                                                                # vector of random numbers used for the decision about initial progression pace
  rc <- runif(n.i)                                                                                # vector of random numbers used for switching from fast to slow disease progression if treated
  v.Group       <- sample(x = group, prob = c(0.5, 0.5), size = n.i, replace = TRUE)              # randomly assigned decline scheme
  v.Treatment   <- ifelse(ra<= ifelse(Int == 0, p.Treatment, p.Treatment_i), 1, 0)                                     # vector indicating treatment availability for a patient
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
  # close the loop for the time points
  }
  
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

#### 06 Run microsimulation ####
# the outcomes are of the simulation are stored in the variables `outcomes_no_int` and `outcomes_int`.
outcomes_no_int  <- MicroSim(n.i, df.X, Int = FALSE, seed = 1)
outcomes_int     <- MicroSim(n.i, df.X, Int = TRUE, seed = 1)

#### 07 Outcomes ####
#### Without intervention
m.M <- outcomes_no_int$m.M

# get the distribution function 
# count the number of individuals in each health state at each cycle
m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
m.TR <- m.TR / n.i                                               # calculate the proportion of individuals 
colnames(m.TR) <- v.n                                            # name the rows of the matrix
rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")               # name the columns of the matrix

# plot trace of first health state
plot_CAU <- plot(0:n.t, m.TR[, 1], type = "l", main = "", 
                 ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
# add a line for each additional state
for (n.s in 2:length(v.n)) {
  lines(0:n.t, m.TR[, n.s], col = n.s)                            # adds a line to current plot
}

title(main = c(" "," ",""))
legend("right", v.n, col=1:4,                                     # add a legend to current plot
       lty = rep(1, 3), bty = "n", cex = 0.65)

#### With intervention
m.M <- outcomes_int$m.M

# get the distribution function 
m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
m.TR <- m.TR / n.i                                                # calculate the proportion of individuals 
colnames(m.TR) <- v.n                                             # name the rows of the matrix
rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")                # name the columns of the matrix

# plot trace of first health state
plot_I <- plot(0:n.t, m.TR[, 1], type = "l", main = "", 
               ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
# add a line for each additional state
for (n.s in 2:length(v.n)) {
  lines(0:n.t, m.TR[, n.s], col = n.s)                            # adds a line to current plot
}

#title(main = c(" "," ",""))
legend("right", v.n, col=1:4,                                     # add a legend to current plot
       lty = rep(1, 3), bty = "n", cex = 0.65)


#### Cost effectiveness
# Table 2.13
# store the mean costs and MSCE of each strategy in a new variable C (vector costs)
v.C <- c(outcomes_no_int$tc_hat, outcomes_int$tc_hat) 
se.C <- c(sd(outcomes_no_int$tc), sd(outcomes_int$tc))/sqrt(n.i)
# store the mean QALYs and MSCE of each strategy in a new variable E (vector effects)
v.E <- c(outcomes_no_int$te_hat, outcomes_int$te_hat)
se.E <- c(sd(outcomes_no_int$te), sd(outcomes_int$te))/sqrt(n.i)
v.tC <- cbind(outcomes_no_int$tc, outcomes_int$tc) 
delta.C <- v.C[2] - v.C[1]                                         # calculate incremental costs
names(delta.C) <- "Incremental costs"                              # name the variable delta.C
delta.E <- v.E[2] - v.E[1]                                         # calculate incremental QALYs
names(delta.E) <- "QALYs gained"                                   # name the variable delta.E
se.delta.C <- sd(outcomes_int$tc - outcomes_no_int$tc)/sqrt(n.i)   # Monte Carlo squared error of incremental costs
se.delta.E <- sd(outcomes_int$te - outcomes_no_int$te)/sqrt(n.i)   # Monte Carlo squared error of incremental QALYs

ICER <- delta.C / delta.E                                          # calculate the ICER
names(ICER)    <- "ICER"                                           # name the variable ICER
results <- c(delta.C, delta.E, ICER)                               # store the values in a new variable

# Create full incremental cost-effectiveness analysis table
v.C <- round(v.C, 0)                                               # round to 0 decimals
v.E <- round(v.E, 2)                                               # round to 2 decimals

strategies <- c("Care as usual", "Early detection")                # store the strategy names
delta.C <- c("", as.character(round(delta.C, 0)))                  # round to 0 decimals
delta.E <- c("", as.character(round(delta.E, 2)))                  # round to 2 decimals
ICER    <- c("", as.character(round(ICER, 0)))                     # round to 0 decimals

table_micro <- cbind(strategies, v.C, v.E, delta.C, delta.E, ICER) # bind the variables 
table_micro <- as.data.frame(table_micro)                          # create a table 
table_micro                                                        # print the table with CE-analysis results


#### Results
# create matrix of interest
# with intervention
checkMMSE_i <- outcomes_int[["m.MMSE"]]
checkE_i <- outcomes_int[["m.E"]]
checkM_i <- outcomes_int[["m.M"]]
checkC_i <- outcomes_int[["m.C"]]
checkInformal_i <- outcomes_int[["m.c.informal"]]
checkOut_i <- outcomes_int[["m.c.outpatient"]]
checkMed_i <- outcomes_int[["m.c.medication"]]
checkScr_i <- outcomes_int[["m.c.screening"]]
ti_i <- outcomes_int[["ti"]]
ts_i <- outcomes_int[["ts"]]
to_i <- outcomes_int[["to"]]
tm_i <- outcomes_int[["tm"]]
checkT_i <- outcomes_int[["m.Streatment"]]
checkF_i <- outcomes_int[["m.F"]]

# without intervention
checkMMSE_cau <- outcomes_no_int[["m.MMSE"]]
checkE_cau <- outcomes_no_int[["m.E"]]
checkM_cau <- outcomes_no_int[["m.M"]]
checkC_cau <- outcomes_no_int[["m.C"]]
checkInformal_cau <- outcomes_no_int[["m.c.informal"]]
checkOut_cau <- outcomes_no_int[["m.c.outpatient"]]
checkMed_cau <- outcomes_no_int[["m.c.medication"]]
checkScr_cau <- outcomes_no_int[["m.c.screening"]]
ti_cau <- outcomes_no_int[["ti"]]
ts_cau <- outcomes_no_int[["ts"]]
to_cau <- outcomes_no_int[["to"]]
tm_cau <- outcomes_no_int[["tm"]]
checkT_cau <- outcomes_no_int[["m.Streatment"]]
checkF_cau <- outcomes_no_int[["m.F"]]


#### Results
# Table 2.14
table(checkM_cau)
table(checkM_i)
table(checkT_i)
table(checkScr_i)
table(checkScr_cau)
table(checkF_i)
table(checkF_cau)


# Difference in cycles lived in Death, Mild, Moderate, Severe
3360832-3359560
370258-425291
144364-198023
224546-117126
# Total cycles lived CAU, I and difference
370398+143947+224823
425103+194651+120686
(370398+143947+224823)-(425103+194651+120686)

# Cost categories (Table 2.15)
informal_cau <- sum(ti_cau)/100000
informal_i <- sum(ti_i)/100000
screening_cau <- sum(ts_cau)/100000
screening_i <- sum(ts_i)/100000
outpatient_cau <- sum(to_cau)/100000
outpatient_i <- sum(to_i)/100000
medication_cau <- sum(tm_cau)/100000
medication_i <- sum(tm_i)/100000

cau <- informal_cau + screening_cau + outpatient_cau + medication_cau
i <- informal_i + screening_i + outpatient_i + medication_i
cau_hc <- screening_cau + outpatient_cau + medication_cau
i_hc <- screening_i + outpatient_i + medication_i
cau - i
hc <- cau_hc - i_hc
informal <- informal_cau - informal_i
outpatient <- outpatient_cau - outpatient_i
screeing <- screening_cau - screening_i
medication <- medication_cau - medication_i

# Per diagnosed
sum(checkF_i)
sum(checkF_cau)
pd_informal_cau <- sum(ti_cau)/100000
pd_informal_i <- sum(ti_i)/100000
pd_screening_cau <- sum(ts_cau)/sum(checkF_cau)
pd_screening_i <- sum(ts_i)/sum(checkF_i)
pd_outpatient_cau <- sum(to_cau)/sum(checkF_cau)
pd_outpatient_i <- sum(to_i)/sum(checkF_i)
pd_medication_cau <- sum(tm_cau)/sum(checkF_cau)
pd_medication_i <- sum(tm_i)/sum(checkF_i)
pd_outpatient <- outpatient_cau - outpatient_i
pd_screeing <- screening_cau - screening_i
pd_medication <- medication_cau - medication_i

cau_hc <- (screening_cau + outpatient_cau + medication_cau)*100000/sum(checkF_cau)
i_hc <- (screening_i + outpatient_i + medication_i)*100000/sum(checkF_i)
hc <- cau_hc - i_hc

# Net benefit (Table 2.13)
k = 45000
benefit_cau = k*outcomes_no_int$te_hat - outcomes_no_int$tc_hat
benefit_i = k*outcomes_int$te_hat - outcomes_int$tc_hat
nb = benefit_i - benefit_cau


