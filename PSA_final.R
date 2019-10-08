############  Microsimulation Alzheimer's Disease in Czechia  ##########
#########             Hana M. Broulikova           #####################
# This code is based on the work of the DARTH group (darthworkgroup.com)
## Jalal H, et al. An Overview of R in Health Decision Sciences. 
# Med. Decis. Making. 2017; 37(3): 735-746. 
## Krijkamp EM, et al. Microsimulation modeling for health decision sciences 
# using R: a tutorial. Med. Decis. Making. 2018; 38(3): 400-422. 
###########################################################################################
# set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set working directory 

# First run the Excersice Markov Sick Sicker model. The next step is to define the model input
source("functions/Functions.R")
source("functions/PSA_functions.R") # load custom made functions to create: 
# Cost-effectiveness plane, 
# Cost-effectiveness acceptability curves/frontiers and
# Expected value of Perfect Information
install.packages(c(
  #                    "reshape2",
  #                    "ellipse",
  #                    "plyr",
  #                    "ggplot2",
  #                    "scales",
  #                    "BCEA"
  "matrixStats"))

# install.packages("rpgm", dependencies = TRUE)
library(plyr) 
library(ggplot2)
library(ellipse)
library(BCEA)
library(tidyverse)
library(truncnorm)
##################################### Model input #########################################
source("functions/Functions.R")
set.seed(1)                                                 # set the seed  
# Model input
#### Model structure
v.Strategies <- c("CAU", "Intervention") # strategy names 
k = 45000                                                   # willingness to pay
n.sim <- 1000     
n.t   <- 40                                                 # time horizon, 40 cycles to achieve lifetime horizon
n.i   <- 10000                                              # number of simulated individuals
v.n   <- c("Mild","Moderate","Severe","Death")              # the model states names
n.s   <- length(v.n)                                        # the number of health states
d.r   <- 0.03                                               # discount rate of 3% per cycle
v.dwe <- v.dwc <- 1 / ((1 + d.r) ^ (0:n.t))                 # discount weight 
group   <- c("Mean", "Lopez")                               # decline schemes considered in this model
mmse    <- 28                                               # intital MMSE at the start of the model
tmmse   <-  runif(n.sim,11, 21)                             # MMSE of diagnosis and treatment initiation 
immse   <-  runif(n.sim,21, 28)   
distr   <- read.csv("sampling.csv") %>% 
  janitor::clean_names() %>% 
  mutate(age = onset)

# load age at disease onset distribution
p.Treatment <- runif(n.sim, 0.25, 1)                        # indicates the probability of receiving treatment without intervention
p.Treatment_i <- 1
p.Slow  <- 0.39                                             # probability of being slow progressor regardless treatment 
p.Fs    <- 0.6                                              # probability of switching from fast to slow progression with treatment
p.Sf    <- 0                                                # probability of switching from slow to fast with treatment
hr      <- 2.1                                              # hazard ratio for AD patients


#### Annual probability of death (Czech mortality tables increased by HR for AD)
p.Death <- read.csv("ut.csv") %>%                           # load Czech mortality tables and adjust by HR
  mutate(Male = 1 - exp(-(- log(1 - qxM))*hr), Female = 1 - exp(-(- log(1 - qxZ))*hr)) %>% 
  select(age, Male, Female) %>% 
  gather(key = "sex" , value = "p.Mortality", 2:3) 


#### Cost inputs (unit costs discounted to the value of 2017)
c.Mild.i      <- rgamma(n.sim, shape = 0.86,  scale = 17797)  # annual costs of indirect care for mild cognitive deficit
c.Moderate.i  <- rgamma(n.sim, shape = 4.89,  scale = 5212)   # annual costs of indirect care for moderate cognitive deficit
c.Severe.i    <- rgamma(n.sim, shape = 2.51,  scale = 12092)  # annual costs of indirect care for severe cognitive deficit
c.Outpatient  <- 97*((1+0.03)^3)                              # annual costs of outpatient care (for treated patients)
c.Mild.d      <- 93*2*((1+0.03)^5)                            # annual costs of drugs for mild cognitive deficit (for treated patients)
c.Moderate.d  <- 93*2*((1+0.03)^5)                            # annual costs of drugs for moderate cognitive deficit (for treated patients)
c.Severe.d    <- 289*2*((1+0.03)^5)                           # annual costs of drugs for severe cognitive deficit (for treated patients)
c.Screening   <- (97*((1+0.03)^3) + (26 + 18)/((1+0.03)^2))*31# costs of screening (for treated patients in the intervention branch)
c.Death       <- 0                                            # annual costs for dead patients 

#### Static characteristics of the cohort
v.MMSE<- rep(mmse, n.i)                                               # a vector with the intitial mmse score 
rg    <- runif(n.i)                                                   # random number for gender sampling
v.Age <- sample(x = distr$onset, prob = distr$freq_age,  size = n.i, replace = TRUE)          
# sample from age distribution an initial age for every individual
df.X  <- data.frame(ID = 1:n.i, age = v.Age, rg = rg, MMSE = v.MMSE)  # binding patient characteristics into the dataframe
df.X  <- merge(df.X, distr, by = "age") %>% 
  mutate(sex = ifelse(rg<freq_gen, "Male", "Female")) %>% 
  select(age, sex, MMSE)                                              # merging probablity of being man to the dataframe, creating gender

#input.psa <- cbind(p.Treatment, c.Mild.i, c.Moderate.i, c.Severe.i)

############################### MCM Model  ###########################

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
Costs <- function (M_it, Streatment, Int= FALSE, i) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Streatment: has the treatment been initiated? 
  # Int: is the individual in the intervention branch? (default is FALSE)
  
  c.it <- 0                                                                                                     # by default the cost for everyone is zero 
  c.it[M_it == "Mild" & Streatment == 0]      <- c.Mild.i[i]                                                    # update the cost if mild untreated
  c.it[M_it == "Mild" & Streatment == 1]      <- c.Mild.i[i] + c.Mild.d+ c.Outpatient + c.Screening*Int         # update the cost if mild detected
  c.it[M_it == "Mild" & Streatment > 1]       <- c.Mild.i[i] + c.Mild.d + c.Outpatient                          # update the cost if mild treated
  c.it[M_it == "Moderate" & Streatment == 0]  <- c.Moderate.i[i]                                                # update the cost if moderate untreated
  c.it[M_it == "Moderate" & Streatment == 1]  <- c.Moderate.i[i] + c.Moderate.d + c.Outpatient + c.Screening*Int# update the cost if moderate detected
  c.it[M_it == "Moderate" & Streatment > 1]   <- c.Moderate.i[i] + c.Moderate.d + c.Outpatient                  # update the cost if moderate treated
  c.it[M_it == "Severe" & Streatment == 0]    <- c.Severe.i[i]                                                  # update the cost if severe untreated
  c.it[M_it == "Severe" & Streatment == 1]    <- c.Severe.i[i] + c.Severe.d + c.Outpatient + c.Screening*Int    # update the cost if severe detected
  c.it[M_it == "Severe" & Streatment > 1]     <- c.Severe.i[i] + c.Severe.d + c.Outpatient                      # update the cost if severe treated
  c.it[M_it == "Death"]  <- c.Death                                                                             # update the cost if dead
  
  return(c.it)        		                                                                                      # return the costs
}

#### Partial cost functions
# Informal care
Costs_i <- function (M_it,i) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Streatment: has the treatment been initiated? 
  # Int: is the individual in the intervention branch? (default is FALSE)
  
  c_inf <- 0                                                                # by default the cost for everyone is zero 
  c_inf[M_it == "Mild"]      <- c.Mild.i[i]                                 # update the cost if mild 
  c_inf[M_it == "Moderate"]  <- c.Moderate.i[i]                             # update the cost if mild
  c_inf[M_it == "Severe"]    <- c.Severe.i[i]                               # update the cost if severe 
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
MicroSim <- function(n.i, df.X , Int = FALSE, seed = 1, i) {
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
  v.Treatment   <- ifelse(ra<= ifelse(Int == 0, p.Treatment[i], p.Treatment_i), 1, 0)                                     # vector indicating treatment availability for a patient
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
  m.C[, 1]          <- Costs(m.M[, 1], m.Streatment[,1], Int, i)                                  # calculate costs per individual during cycle 0
  m.E[, 1]          <- 0.408+0.01*df.X$MMSE
  m.c.informal[, 1] <- Costs_i(m.M[, 1],i)                                                        # calculate informal costs per individual during cycle 0
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
    m.C[, t + 1]    <- Costs(m.M[, t + 1], m.Streatment[,t], Int, i)                              # calculate costs per individual during cycle t + 1
    m.c.informal[, t + 1]    <- Costs_i(m.M[, t + 1], i)                                          # calculate informal costs per individual during cycle t + 1
    m.c.outpatient[, t + 1]  <- Costs_o(m.M[, t + 1], m.Streatment[,t])                           # calculate outpatient costs per individual during cycle t + 1
    m.c.screening[, t + 1]   <- Costs_s(m.M[, t + 1], m.Streatment[,t], Int)                      # calculate screening costs per individual during cycle t + 1
    m.c.medication[, t + 1]  <- Costs_m(m.M[, t + 1], m.Streatment[,t])                           # calculate medication costs per individual during cycle t + 1
    m.F[, t + 1]  <- Fraction(m.M[, t + 1], m.Streatment[,t])                                     # calculate medication costs per individual during cycle t + 1
    #m.E[, t + 1]    <- Effs (m.M[, t + 1], df.X)                                                 # calculate QALYs per individual during cycle t + 1
    m.E[, t + 1]    <- ifelse(m.M[, t + 1] != "Death", 0.408+0.01*df.X$MMSE,0)                    # calculate QALYs per individual during cycle t + 1
    # update time of treatment 
    m.Streatment[, t + 1]            <- ifelse((Int == 0 & v.Treatment==1 & df.X$MMSE <= tmmse[i] & m.M[, t + 1] != "Death") 
                                               |(Int == 1 & v.Treatment==1 & df.X$MMSE <= immse[i] & m.M[, t + 1] != "Death"), m.Streatment[,t] + 1,
                                               ifelse(m.Streatment[,t] > 0 & m.M[, t + 1] != "Death", m.Streatment[,t] + 1,
                                                      m.Streatment[,t]))
    # update age
    df.X$age[m.M[,t + 1] != "Death"] <- df.X$age[m.M[, t + 1] != "Death"] + 1
  } # close the loop for the time points
  
  # Display simulation progress
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

# initialize variables that will be stored in the probabilistic simulation loop
te_no_trt <- te_trt <- tc_no_trt <- tc_trt <- vector("numeric", n.sim)

# start the clock and run the model
p <- Sys.time()  # store system time 

for (i in 1:n.sim) {
  #### 06 Run microsimulation ####
  # the outcomes are of the simulation are stored in the variables `outcomes_no_int` and `outcomes_int`.
  outcomes_no_int  <- MicroSim(n.i, df.X, Int = FALSE, seed = 1, i)
  outcomes_int     <- MicroSim(n.i, df.X, Int = TRUE, seed = 1, i)
  
  #store
  
  te_no_trt[i] <- outcomes_no_int$te_hat # discount costs by multiplying the utility vector for no treatment with discount weights for effects  
  te_trt[i]    <- outcomes_int$te_hat # discount costs by multiplying the utility vector for no treatment with discount weights for effects  
  
  tc_no_trt[i] <- outcomes_no_int$tc_hat # discount costs by multiplying the cost vector for no treatment with discount weights for costs  
  tc_trt[i]    <- outcomes_int$tc_hat # discount costs by multiplying the cost vector for treatment with discount weights for costs  
}


Sys.time() - p # calculate time to run the analusis by extracting p from the current system time  

par(mfrow= c(2,2))
hist(tc_no_trt, main = "Costs CAU",xlab = "Costs")    # histogram of total discounted cost for control 
hist(te_no_trt, main = "QALY CAU", xlab = "QALY")     # histogram of total discounted QALYs for control 

hist(tc_trt, main = "Costs TT", xlab = "Costs")       # histogram of total discounted cost for treatment
hist(te_trt, main = "QALY TT", xlab = "QALY")         # histogram of total discounted QALYs for treatment
par(mfrow = c(1,1))
v.wtp <- c(1, seq(5000, 200000, length.out = 31)) # create vector with willingness to pay values

te <- data.frame(te_no_trt, te_trt)          # create data frame of effectiveness for both strategies
tc <- data.frame(tc_no_trt, tc_trt)          # create data frame of costs for both strategies
colnames(te) <- colnames(tc) <- v.Strategies # name the columns

scatterE(strategies = v.Strategies, m.c = tc, m.e = te)             # create cost-effectiveness plane
ceaf(v.wtp = v.wtp, strategies = v.Strategies, m.e = te, m.c = tc)  # create cost-effectiveness acceptability curve
evpiE(v.wtp = v.wtp, m.e = te, m.c = tc)                            # create expected value of information plot

m <- bcea(as.matrix(te), as.matrix(tc), ref= 2, Kmax = 200000)      # fit a bayesian cost-effectivenesss analysis object
ceplane.plot(m, wtp = 45000)
ceac.plot(m)


evi.plot(m)
evppi.res <- evppi(parameter = 1:4,input = input.psa , he = m, method = "gam", plot= T) #calculate expected value of partial perfect information
plot.evppi(evppi.res)


#output
psa_params <- as.data.frame(cbind(p.Treatment, tmmse, immse, c.Mild.i, c.Moderate.i, c.Severe.i))
psa_results <- as.data.frame(cbind(tc_no_trt, tc_trt, te_no_trt, te_trt))
psa_results <- mutate(psa_results, delta_c = tc_trt - tc_no_trt, delta_e = te_trt - te_no_trt, icer = delta_c/delta_e, benefit_cau = k*te_no_trt, benefit_i = k*te_trt, net_benefit =  k*delta_e - delta_c, 
                      oportunity = ifelse(net_benefit < 0, - net_benefit, 0))
sum(psa_results$oportunity)/1000
summary(tc_no_trt)
summary(tc_trt)
summary(te_no_trt)
summary(te_trt)
sd(tc_no_trt)
sd(tc_trt)
sd(te_no_trt)
sd(te_trt)
mean(tc_no_trt) - mean(tc_trt)
mean(te_no_trt) - mean(te_trt)
temp <- filter(psa_results, delta_c <= 0)
temp <- filter(psa_results, delta_e < 0)
temp <- filter(psa_results, te_trt > 3.784)
temp1 <- filter(psa_results, icer < 45000) 
temp2 <- filter(psa_results, net_benefit < 0) 
temp3 <- filter(psa_results, net_benefit > 0 & icer < 45000) 
temp <- filter(psa_results, icer < 0)
mean(temp2$oportunity)

install.packages("moments")
library(moments)
skewness(tc_trt)
skewness(tc_no_trt)
skewness(te_trt)
skewness(te_no_trt)
kurtosis(tc_trt)
kurtosis(tc_no_trt)
kurtosis(te_trt)
kurtosis(te_no_trt)











