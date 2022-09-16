###############################################################
### This script produces PEs for PAFs and produces          ###
### files of the adjusted mean consumption values and        ###
### AAF point estimates for each age, sex and country       ###
### for the relevant cancers.                               ###
###                                                         ###
### This script uses 3 source files:                        ###
### 2.aaf_function.R - contains the Mean and AAF functions  ###
### 3.rr_function.R - contains the RR functions             ###
### 3.rr_function_split.R - contains the RR functions       ###
### split by consumption level                              ###
###############################################################

# Clear everything in the global environment
rm(list = ls())

# run the following libraries
library(parallel)
library(MASS)
library(Hmisc)
library(foreach)
library(doParallel)
library(tidyverse)


YEAR_INCIDENCE <- 2020 # insert year of cancer incidence data
YEAR_ALCOHOL 	<- YEAR_INCIDENCE - 10
SEX_LIST 		<- c("men", "women")
SEX_NUM_LIST 	<- c(1:2)
AGE_NUM_LIST 	<- c(1:6)

mycores <- 1 ## should be 1 for windows machines
myverbose <- TRUE
MyTestStatus <- "no"

# import alcohol PCA data
setwd("//Inti/cin/Users/RumgayH/alcohol project/IARC Cancer project/GLOBAL BURDEN OF CANCERS ATTRIBUTABLE TO ALCOHOL/DATA_2010")
drkData 			<- read.csv("DATA_INPUT_APC_LANCET_CAPPED_2019 - v2.csv")

#drkData %>% filter(year == YEAR_ALCOHOL) %>% dplyr::select(iso3a) %>% unique() %>% View() # check which countries have alc data - 189
drkData$BINGE_A <- 60  # introduce bingers variable with value 60
drkData$FORMER_DRINKERS 			<- drkData$FORMER_DRINKERS_CAPED
drkData$LIFETIME_ABSTAINERS 		<- drkData$LIFETIME_ABSTAINERS_CAPED
drkData$PCA <- drkData$PCA_SEX   # set PCA value to PCA sex value


# set wd to run source code for PAF calculations
setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs")
source("2.aaf_function.R")        # load PAF functions
source("3.rr_function.R") # run source code for RR function

# limit data to selected years and create new datasets
drkData_2010 <- drkData[drkData$year == YEAR_ALCOHOL,]

## calculate adjusted mean consumption by strata for selected year
meanCon_2010 <- computeMean(drkData_2010)
meanCon_2010_M <- subset(meanCon_2010, SEX == 1)
meanCon_2010_F <- subset(meanCon_2010, SEX == 2)

# create list of cancer years datasets
ALC_CANCER <- list(meanCon_2010_M , meanCon_2010_F)


## CANCER PAF ##
## CANCER PAF ##
# calculate PAF by sex
# calculateAAF1 uses new calculation which matches split drinkers (matches cases in split calculation exactly)
# calculateAAF2 uses formula for AAF3 but includes former drinkers in denom
# calculateAAF3 uses formula for AAF3 but includes former drinkers in numerator and denom - total incl former

# calculate PAF for current drinkers 
AAF_CANCER <- do.call(rbind, lapply(SEX_NUM_LIST, function(i1) {
  if (i1 == 1) {
    aaf_out <-
      calculateAAF1(data = ALC_CANCER[[1]],
                   disease = RRslist_m, 
                   mc.cores = 1)
  }
  if (i1 == 2) {
    aaf_out <-
      calculateAAF1(data = ALC_CANCER[[2]],
                   disease = RRslist_f, 
                   mc.cores = 1)
  }
  
  aaf_out$YEAR <- YEAR_ALCOHOL
  
  return(aaf_out)
}))


# calculate PAF by sex for former drinkers
AAF_CANCER_former <- do.call(rbind, lapply(SEX_NUM_LIST, function(i1) {
  if (i1 == 1) {
    aaf_out <-
      calculateAAFformer(data = ALC_CANCER[[1]],
                   disease = RRslist_m, 
                   mc.cores = 1)
  }
  if (i1 == 2) {
    aaf_out <-
      calculateAAFformer(data = ALC_CANCER[[2]],
                   disease = RRslist_f, 
                   mc.cores = 1)
  }
  
  aaf_out$YEAR <- YEAR_ALCOHOL
  
  return(aaf_out)
}))

# calculate PAF by sex for current drinkers including former drinkers in denom
AAF_CANCER_formercd <- do.call(rbind, lapply(SEX_NUM_LIST, function(i1) {
  if (i1 == 1) {
    aaf_out <-
      calculateAAF2(data = ALC_CANCER[[1]],
                          disease = RRslist_m, 
                          mc.cores = 1)
  }
  if (i1 == 2) {
    aaf_out <-
      calculateAAF2(data = ALC_CANCER[[2]],
                          disease = RRslist_f, 
                          mc.cores = 1)
  }
  
  aaf_out$YEAR <- YEAR_ALCOHOL
  
  return(aaf_out)
}))

setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs")
source("3.rr_function_split.R")  # run source code for RR function
## CANCER PAF ##
## CANCER PAF ##
# calculate PAF by sex by consumption level
AAF_CANCER_split <- do.call(rbind, lapply(SEX_NUM_LIST, function(i1) {
  if (i1 == 1) {
    aaf_out <-
      calculateAAF_split(data = ALC_CANCER[[1]],
                   disease = RRslist_split_m, 
                   mc.cores = 1)
  }
  if (i1 == 2) {
    aaf_out <-
      calculateAAF_split(data = ALC_CANCER[[2]],
                   disease = RRslist_split_f, 
                   mc.cores = 1)
  }
  
  aaf_out$YEAR <- YEAR_ALCOHOL
  
  return(aaf_out)
}))

AAF_CANCER_split %>% 
  bind_rows(AAF_CANCER_former %>% 
              mutate(DISEASE = paste0(DISEASE,"_former"))) %>% 
  bind_rows(AAF_CANCER_formercd %>% 
              mutate(DISEASE = paste0(DISEASE,"_formercd")))-> AAF_CANCER_split

#setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs")
source("3.rr_function_split10.R")  # run source code for RR function for consumption levels of 10g
## CANCER PAF ##
## CANCER PAF ##
# calculate PAF by sex
AAF_CANCER_split10 <- do.call(rbind, lapply(SEX_NUM_LIST, function(i1) {
  if (i1 == 1) {
    aaf_out <-
      calculateAAF_split(data = ALC_CANCER[[1]],
                         disease = RRslist_split_m, 
                         mc.cores = 1)
  }
  if (i1 == 2) {
    aaf_out <-
      calculateAAF_split(data = ALC_CANCER[[2]],
                         disease = RRslist_split_f, 
                         mc.cores = 1)
  }
  
  aaf_out$YEAR <- YEAR_ALCOHOL
  
  return(aaf_out)
}))

#### ---- WRITE OUTPUT FILES ----####
setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_results")
write.csv(AAF_CANCER, "AAF_CANCER_01.02.21.csv", row.names = FALSE)
write.csv(AAF_CANCER_split, "AAF_CANCER_split_01.02.21.csv", row.names = FALSE)
write.csv(AAF_CANCER_split10, "AAF_CANCER_split10_07.04.21.csv", row.names = FALSE)

# remove variables created that won't be used again
rm(list=ls())


