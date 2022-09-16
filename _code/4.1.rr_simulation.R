
#### SIMULATING RR FUNCTIONS ####
#### SIMULATING RR FUNCTIONS ####

library(parallel)
library(MASS)

set.seed(2626)
# set up 1,000 simulations - must change files simulating parameters if we want more than 1,000 simulations
nnn <- 1000 

setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs")

# make list of RRs
source("3.rr_function.R")
source("3.rr_function_split.R")

#RRslist
#RRslist_split


### RR_SIM_FILE ###
### RR_SIM_FILE ###
# create list by sex
RR_SIM <- function(dis) { 
	list(disease = dis$disease, 
		RRCurrent = dis$RRCurrent, 
		betaCurrent = dis$betaCurrent, 
		beta_sim = mvrnorm(nnn, mu = dis$betaCurrent, Sigma = dis$covBetaCurrent),  #one set = one row
		lnRRFormer = dis$lnRRFormer,
		lnRRFormer_sim = rnorm(nnn, mean = dis$lnRRFormer, sd = sqrt(dis$varLnRRFormer)), 
		lowerlim = dis$lowerlim,
		midlim = dis$midlim,
		upperlim = dis$upperlim
		)
	}

RR_CANCER_SIM 		<- list(lapply(RRslist_m, RR_SIM), lapply(RRslist_f, RR_SIM))


RR_SIM_SPLIT <- function(dis) { 
  list(disease = dis$disease, 
       RRCurrent = dis$RRCurrent, 
       betaCurrent = dis$betaCurrent, 
       beta_sim = mvrnorm(nnn, mu = dis$betaCurrent, Sigma = dis$covBetaCurrent),  #one set = one row
       #       lnRRFormer = dis$lnRRFormer,
       #       lnRRFormer_sim = rnorm(nnn, mean = dis$lnRRFormer, sd = sqrt(dis$varLnRRFormer)), 
       lowerlim = dis$lowerlim,
       midlim = dis$midlim,
       upperlim = dis$upperlim,
       toplim = dis$toplim
  )
}

RR_CANCER_SIM_SPLIT 		<- list(lapply(RRslist_split_m, RR_SIM_SPLIT), lapply(RRslist_split_f, RR_SIM_SPLIT))

source("3.rr_function_split10.R")
RR_SIM_SPLIT10 <- function(dis) { 
  list(disease = dis$disease, 
       RRCurrent = dis$RRCurrent, 
       betaCurrent = dis$betaCurrent, 
       beta_sim = mvrnorm(nnn, mu = dis$betaCurrent, Sigma = dis$covBetaCurrent),  #one set = one row
       #       lnRRFormer = dis$lnRRFormer,
       #       lnRRFormer_sim = rnorm(nnn, mean = dis$lnRRFormer, sd = sqrt(dis$varLnRRFormer)), 
       lowerlim = dis$lowerlim,
       midlim = dis$midlim,
       upperlim = dis$upperlim,
       toplim = dis$toplim
  )
}

RR_CANCER_SIM_SPLIT10 		<- list(lapply(RRslist_split_m, RR_SIM_SPLIT10), lapply(RRslist_split_f, RR_SIM_SPLIT10))

setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/simulations")

save(RR_CANCER_SIM, file = "SIM_CANCER_RR_HR.RData")
save(RR_CANCER_SIM_SPLIT, file = "SIM_CANCER_RR_split_HR.RData")
save(RR_CANCER_SIM_SPLIT10, file = "SIM_CANCER_RR_split10_HR.RData")
rm(RRslist_m, RRslist_f, RRslist_split_m, RRslist_split_f, RR_CANCER_SIM, RR_CANCER_SIM_SPLIT)
