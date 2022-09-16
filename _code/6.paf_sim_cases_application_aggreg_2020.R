######################################################################
### This script calculates PAFs of cases using the point estimates ###
### produced in the script 4.2.aaf_simulation.R and Globocan 2020 data     ###
######################################################################

# Clear everything in the global environment
rm(list = ls())


# Run the following libraries
library(tidyverse)
library(data.table)
library(Rcan)
library(parallel)
library(MASS)

source("~/functions.R")
numCores = 16

#setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_results")
#AAF_CANCER 	<- as.data.table(read.csv("AAF_CANCER_11.08.20.csv"))
#AAF_CANCERf <- AAF_CANCER[SEX==1 & AGE_CATEGORY==5,]
#AAF_CANCERf %>% rename(ISO3 = REGION,
#                       disease = DISEASE) -> AAF_CANCERf
#paf <- PAF_DATA_1_5
#for (CAN in unique(PAF_DATA_1_5$disease)){
#  for (REG in unique(PAF_DATA_1_5$ISO3)){
#
#    cat("sex:", toString(CAN),", reg:",toString(REG), "\n") # print sex and country number
#
#    paf[disease==CAN & ISO3==REG, LCI := quantile(paf[disease==CAN & ISO3==REG,3:1002], probs = 0.025, na.rm = TRUE)]
#    paf[disease==CAN & ISO3==REG, UCI := quantile(paf[disease==CAN & ISO3==REG,3:1002], probs = 0.975, na.rm = TRUE)]
#  }
#}
#comp <- merge(paf[,c(1,2,1005:1006)],AAF_CANCERf,by=c("ISO3","disease"), all=T)
#comp %>% dplyr::select(ISO3, disease, LCI, AAF, UCI) %>% filter(UCI<AAF) %>% View()

####--- Import data ---####
# import Globocan data - need to update to 2020 data
# import Globocan data
setwd("~/AlcPAF_os/alc_data")
DATA_GLOB	<- as.data.table(read.csv("Globocan2020.csv"))
glob_inc18	<- as.data.table(read.csv("Globocan2018.csv"))
# Globocan incidence = type 0, mort = type 1
# Globocan persons = sex 0, men = sex 1, women = sex 2
# and apply 2018 HCC and SCC proportions
#setwd("//Inti/cin/DataShare/Globocan2018")
#DATA_GLOB	<- as.data.table(read.csv("Globocan.csv"))
# Globocan incidence = type 0, mort = type 1
# Globocan persons = sex 0, men = sex 1, women = sex 2

# import HCC and oesophageal SCC cases
HCCcases <-
  as.data.table(
    read.csv(
      "~/AlcPAF_os/alc_data/cases_method1_18.05.20.csv"
    )
  )
OSCCcases <-
  as.data.table(read.csv(
    "~/AlcPAF_os/alc_data/OSCC_MA.csv"
  ))

# create list of regions and countries to match
setwd("~/AlcPAF_os/alc_results")
countrymatch <-
  as.data.table(read.csv("pafs_inc_country_WCRF_02.02.21.csv"))
countrymatch %>%
  dplyr::select(country_code, subregion_label_2) %>%
  unique() -> subregmatch
countrymatch %>%
  dplyr::select(country_code, hdigrp) %>%
  unique() -> hdimatch
countrymatch %>%
  dplyr::select(country_code, who_region_label) %>%
  unique() -> whomatch
countrymatch %>%
  dplyr::select(country_code, region_label) %>%
  unique() -> contmatch
countrymatch %>%
  dplyr::select(country_code:ISO3, Development:hdigrp) %>%
  unique() -> countrymatch

# load age structure csv
setwd("~/AlcPAF_os/alc_data")
AGE_STRUCT <-
  as.data.table(read.csv("PAF_INC_STRUCTURE_AGE_HR.csv"))

# import populations data from Globocan for calculating rates
pop_glob <- read.csv("pops_globocan.csv")
# rename age variable to same as pafs dataset
names(pop_glob)[names(pop_glob) == "age"] <- "age_num"


####--- Reshape Globocan data ---####
# filter Globocan incidence data
glob_inc <-
  DATA_GLOB[(type == 0), ,][order(country_code, sex, age),][cancer_code %in% c(1, 3, 5, 6, 7, 8, 9, 11, 13, 14, 20, 39, 40),] #[age>=6,] # select incidence data ages 25+
glob_inc <-
  glob_inc[country_code < 900, ][sex != 0, ] # filter to countries only and remove persons cases from Globocan data
glob_inc %>% dplyr::select(-country_label, -type, -py) -> glob_inc

# filter Globocan 2018 data for OSCC and HCC
glob_inc18 <-
  glob_inc18[(type == 0), ,][order(country_code, sex, age),][cancer_code %in% c(6, 11),]
glob_inc18 <-
  glob_inc18[country_code < 900, ][sex != 0, ] # filter to countries only and remove persons cases from Globocan data
glob_inc18 %>% dplyr::select(-country_label, -type, -py) -> glob_inc18

# reshape HCC cases to same format as Globocan data
HCCcases %>%
  filter(count.name == "HCC") %>%
  mutate(count.name = "Liver_HCC") %>%
  dplyr::select(COUNTRY_CODE, SEX_LABEL, AGEG, count.name, count.value) %>%
  rename(country_code = COUNTRY_CODE,
         cancer_label = count.name,
         cases = count.value) %>%
  mutate(
    cancer_code = 11.1,
    # to match RRs list
    sex = case_when(SEX_LABEL == "Men" ~ 1,
                    SEX_LABEL == "Women" ~ 2),
    age = case_when(
      AGEG == "0-4" ~ 1,
      AGEG == "5-9" ~ 2,
      AGEG == "10-14" ~ 3,
      AGEG == "15-19" ~ 4,
      AGEG == "20-24" ~ 5,
      AGEG == "25-29" ~ 6,
      AGEG == "30-34" ~ 7,
      AGEG == "35-39" ~ 8,
      AGEG == "40-44" ~ 9,
      AGEG == "45-49" ~ 10,
      AGEG == "50-54" ~ 11,
      AGEG == "55-59" ~ 12,
      AGEG == "60-64" ~ 13,
      AGEG == "65-69" ~ 14,
      AGEG == "70-74" ~ 15,
      AGEG == "75-79" ~ 16,
      AGEG == "80-84" ~ 17,
      AGEG == "85+" ~ 18
    )
  ) %>%
  dplyr::select(-SEX_LABEL,-AGEG) %>%
  group_by(sex, country_code) %>%
  mutate(total = sum(cases)) %>%
  bind_rows(
    OSCCcases %>%
      filter(un_code != 1000) %>%
      rename(country_code = un_code,
             cases = N) %>%
      mutate(
        cancer_code = 6.1,
        # to match RRs list
        sex = case_when(sex == "Males" ~ 1,
                        sex == "Females" ~ 2)
      ) %>%
      group_by(sex, country_code) %>%
      mutate(
        total = sum(cases),
        cancer_label = paste0(cancer_label, "_SCC")
      ) %>%
      dplyr::select(country_code, cancer_label, cancer_code, sex, age, cases, total)
  ) -> HCCcases

# add HCC cases to globocan dataset - incidence
glob_inc18 <- bind_rows(glob_inc18, HCCcases)

# need to calculate distribution of subtype of liver cancer and oesophageal cancer to apply to cases
# proportion of SCC <64, 65+ - group cases <64 (age 1-13) and 65+ (age 14-18)
# proportion of HCC - group all ages
glob_inc18 %>%
  group_by(country_code, age, sex) %>%
  mutate(prop = case_when(
    cancer_code == 11.1 ~ cases[cancer_code == 11],
    cancer_code == 6.1 ~ cases[cancer_code == 6],
    TRUE ~ cases
  )) %>%
  filter(cancer_code %in% c(11.1, 6.1)) %>%
  mutate(grp = case_when(cancer_code == 6.1 &
                           age > 13 ~ 2, #split age groups for oesoph
                         TRUE ~ 1)) %>%
  group_by(country_code, grp, sex, cancer_code) %>%
  mutate(prop = sum(cases) / sum(prop)) %>%
  ungroup() %>%
  dplyr::select(-grp) -> glob_inc18

# apply 2018 proportions of HCC and OSCC to 2020 data
glob_inc18 %>%
  bind_rows(glob_inc %>%
              mutate(cancer_code = as.numeric(cancer_code))) %>%
  group_by(country_code, age, sex) %>%
  mutate(cases = case_when(
    cancer_code == 11.1 ~ cases[cancer_code == 11] * prop,
    cancer_code == 6.1 ~ cases[cancer_code == 6] *
      prop,
    TRUE ~ cases
  )) %>%
  group_by(country_code, sex) %>%
  mutate(
    total = case_when(
      cancer_code == 11.1 ~ sum(cases[cancer_code == 11.1]),
      cancer_code == 6.1 ~ sum(cases[cancer_code ==
                                       6.1]),
      TRUE ~ total
    ),
    cases = case_when(is.na(cases) ~ 0, TRUE ~ cases),
    total = case_when(is.na(total) ~ 0, TRUE ~ total)
  ) %>%
  dplyr::select(-prop) -> glob_inc
rm(glob_inc18, DATA_GLOB)

names(glob_inc)[names(glob_inc) == "age"] <-
  "age_num" # rename age variable
glob_inc <-
  merge(AGE_STRUCT, glob_inc, by = "age_num") # merge PAF age variable into inc data


####--- Match alcohol data and Globocan countries ---####
setwd("~/AlcPAF_os/alc_data")
regionsmatch 	<-
  as.data.table(read.csv("countries_region_matching.csv"))
regionsmatch %>%
  dplyr::select(M49Code, ISO3) %>%
  unique() %>%
  rename(country_code = M49Code) -> regionsmatch
# merge regions match with globocan data
glob_inc <- merge(glob_inc, regionsmatch, by = "country_code")
glob_inc$paf_age_cancer <- as.numeric(glob_inc$paf_age_cancer)

####--- Match cancer names ---####

setwd("~/AlcPAF_os/alc_data")
RRs	<- as.data.table(read.csv("WCRF_RRs3.csv"))
RRsplit	<- as.data.table(read.csv("WCRF_RRs_split3.csv"))
RRsplit10	<- as.data.table(read.csv("WCRF_RRs_split10.csv"))
# filter to non-regional RRs
RRs <- RRs[Regional == 0]
RRsplit <- RRsplit[Regional == 0]
RRsplit10 <- RRsplit10[Regional==0]
# separate disease name and Globocan number
Gcode <- RRs[sex == 1, c("disease", "cancer_code", "Main_analysis")]
Gcodes <-
  RRsplit[sex == 1, c("disease", "cancer_code", "Main_analysis")]
Gcodes10 <- RRsplit10[sex==1, c("disease", "cancer_code", "Main_analysis")] 
# add nonlin names and Globocan codes
Gcode <-
  rbind(Gcode, data.frame(
    c("Oesophagus_SCC_nonlin", "Colon_lin", "Rectum_lin"),
    c(6.1, 8, 9),
    c(1, 1, 1)
  ), use.names = FALSE)
Gcodes <-
  rbind(Gcodes, data.frame(
    c(
      "Oesophagus_SCC_nonlin_light",
      "Oesophagus_SCC_nonlin_mod",
      "Oesophagus_SCC_nonlin_heavy",
      "Colon_lin_light",
      "Colon_lin_mod",
      "Colon_lin_heavy",
      "Rectum_lin_light",
      "Rectum_lin_mod",
      "Rectum_lin_heavy"
    ),
    c(6.1, 6.1, 6.1, 8, 8, 8, 9, 9, 9),
    c(1, 1, 1, 1, 1, 1, 1, 1, 1)
  ), use.names = FALSE)
Gcodes10 <- rbind(Gcodes10, data.frame(c("Oesophagus_SCC_nonlin_1", "Oesophagus_SCC_nonlin_2", "Oesophagus_SCC_nonlin_3", "Oesophagus_SCC_nonlin_4", "Oesophagus_SCC_nonlin_5",
                                         "Oesophagus_SCC_nonlin_6", "Oesophagus_SCC_nonlin_7", "Oesophagus_SCC_nonlin_8", "Oesophagus_SCC_nonlin_9", "Oesophagus_SCC_nonlin_10",
                                         "Oesophagus_SCC_nonlin_11", "Oesophagus_SCC_nonlin_12", "Oesophagus_SCC_nonlin_13", "Oesophagus_SCC_nonlin_14", "Oesophagus_SCC_nonlin_15"), 
                                       c(6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1), 
                                       c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)), use.names=FALSE)


####---- Import country simulations ####

# load simulations data
year <- 2020 # set year
# change working directory
setwd("~/AlcPAF_os/simulations")
# load simulated point estimates introducing age and sex variable
load(file = paste0("CANCER_1_1_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <- 1
PAF_DATA_1_1 <- aaf_out_PAF
load(file = paste0("CANCER_1_2_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <- 1
PAF_DATA_1_2 <- aaf_out_PAF
load(file = paste0("CANCER_1_3_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <- 1
PAF_DATA_1_3 <- aaf_out_PAF
load(file = paste0("CANCER_1_4_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <- 1
PAF_DATA_1_4 <- aaf_out_PAF
load(file = paste0("CANCER_1_5_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <- 1
PAF_DATA_1_5 <- aaf_out_PAF
load(file = paste0("CANCER_1_6_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <- 1
PAF_DATA_1_6 <- aaf_out_PAF

load(file = paste0("CANCER_2_1_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <- 2
PAF_DATA_2_1 <- aaf_out_PAF
load(file = paste0("CANCER_2_2_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <- 2
PAF_DATA_2_2 <- aaf_out_PAF
load(file = paste0("CANCER_2_3_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <- 2
PAF_DATA_2_3 <- aaf_out_PAF
load(file = paste0("CANCER_2_4_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <- 2
PAF_DATA_2_4 <- aaf_out_PAF
load(file = paste0("CANCER_2_5_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <- 2
PAF_DATA_2_5 <- aaf_out_PAF
load(file = paste0("CANCER_2_6_", (year - 10), "new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <- 2
PAF_DATA_2_6 <- aaf_out_PAF

setwd("~/AlcPAF_os/alc_results")

sim <-
  list(
    rbind(
      PAF_DATA_1_1,
      PAF_DATA_1_2,
      PAF_DATA_1_3,
      PAF_DATA_1_4,
      PAF_DATA_1_5,
      PAF_DATA_1_6
    ),
    rbind(
      PAF_DATA_2_1,
      PAF_DATA_2_2,
      PAF_DATA_2_3,
      PAF_DATA_2_4,
      PAF_DATA_2_5,
      PAF_DATA_2_6
    )
  )


rm(
  PAF_DATA_1_1,
  PAF_DATA_1_2,
  PAF_DATA_1_3,
  PAF_DATA_1_4,
  PAF_DATA_1_5,
  PAF_DATA_1_6,
  PAF_DATA_2_1,
  PAF_DATA_2_2,
  PAF_DATA_2_3,
  PAF_DATA_2_4,
  PAF_DATA_2_5,
  PAF_DATA_2_6,
  aaf_out_PAF
)


# sim includes ISO3, disease name, paf_age_cancer, sex and 1000 simulated PEs
# add Globocan code to dt using disease names to match with globocan incidence data
# sim_t <- merge(sim, Gcode, by = "disease") # regional oesoph aren't included in RRs list so removed
sim_t <- lapply(1:2, function(SEX) {
  t <- merge(sim[[SEX]], Gcode, by = "disease")
  t %>%
    bind_rows(t %>% 
                filter(disease=="Pharynx") %>% # add rows for hypopharynx
                mutate(cancer_code=5)) %>% 
    rename_at(vars(3:1002), function(x)
      paste0("sim", x)) %>%
    filter(!str_detect(disease, "_lin")) -> t
})
rm(sim)


# merge simulations and globocan incidence datasets by sex, age, country and cancer type
sim_t. <- lapply(1:2, function(SEX) {
  t <-
    as.data.table(merge(
      glob_inc[sex == SEX, ],
      sim_t[[SEX]],
      by = c("cancer_code", "paf_age_cancer", "ISO3", "sex"),
      all = T,
      allow.cartesian = TRUE
    ))
  
  t2 <- lapply(unique(t$ISO3), function(REG) {
    t3 <- t[ISO3 == REG, ]
  })
  
})

# calculate attributable cases per sex, country and simulation number, replace AAF with attributable cases
sim_t2 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t2 <- as.data.table(sim_t.[[SEX]][[REG]])
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and country number
    
    for (i in 12:1011) {
      t2[, as.numeric(i)] <- t2[, ..i] * t2$cases
    }
    
    # remove premenop results in older ages
    t2 <- t2[(cancer_code == 20 & age_num < 6), disease := "Breast"]
    t2 <-
      t2[(
        cancer_code == 20 &
          age_num > 10 &
          country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)
      ), disease := "Breast_postmen"]
    t2 <-
      t2[(
        cancer_code == 20 &
          age_num %in% c(6:10) &
          country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)
      ), disease := "Breast_premen"]
    t2 <-
      t2[!(cancer_code == 20 &
             disease == "Breast_premen" & age_num > 10), ]
    # remove postmenop results in younger ages
    t2 <-
      t2[!(cancer_code == 20 &
             disease == "Breast_postmen" & age_num <= 10), ]
    # assing 0 to male breast cancer cases
    t2[cancer_code == 20 & sex == 1, c(12:1011)] <- 0
    
    t2
  })
  
})
rm(sim_t, sim_t.)

# add all cancer totals for Main analysis and sensitivity (include stomach & panc)
# also assign HCC and SCC totals to liver and oesophagus
sim_t3m <- lapply(1:199, function(REG) {
  t2 <- sim_t2[[1]][[REG]]
  cat("sex:", toString(1), ", reg:", toString(REG), "\n") # print sex and country number
  
  # summarise multiple columns with case_when creating total attributable cases for all ages per cancer type
  t2 %>%
    group_by(age_num) %>%
    mutate_at(vars(12:1011), funs(case_when(
      cancer_code %in% c(39, 40) ~ sum(., na.rm = T),
      TRUE ~ .
    ))) %>%
    bind_rows(
      t2 %>%
        group_by(age_num) %>%
        mutate_at(vars(12:1011), funs(case_when(
          cancer_code %in% c(39, 40) ~ sum(.[Main_analysis == 1], na.rm = T)
        ))) %>%
        mutate(cancer_label = case_when(
          cancer_code %in% c(39, 40) ~ paste0(cancer_label, "_M")
        )) %>% # add label to main all cancer results
        filter(cancer_code %in% c(39, 40))
    ) %>%
    group_by(age_num) %>%
    mutate_at(vars(12:1011), funs(
      case_when(
        cancer_code == 11 ~ sum(.[cancer_code == 11.1], na.rm = T),
        cancer_code == 6 ~ sum(.[cancer_code ==
                                   6.1], na.rm = T),
        TRUE ~ .
      )
    )) -> t3
  
  t3 %>% 
    bind_rows(t3 %>% 
                filter(cancer_code %in% c(3,5)) %>% 
                            group_by(age_num) %>% 
                            mutate(cancer_label = "Pharynx",
                                   disease= "Pharynx",
                                   cancer_code = 41,
                                   cases = sum(cases),
                                   total = sum(total)) %>%  
                mutate_at(vars(12:1011), funs(sum(., na.rm = T))) %>% 
                            unique()) -> t3
  
  t4 <- as.data.table(t3)
  
  t4
  
})
sim_t3f <- lapply(1:199, function(REG) {
  t2 <- sim_t2[[2]][[REG]]
  cat("sex:", toString(2), ", reg:", toString(REG), "\n") # print sex and country number
  
  # summarise multiple columns with case_when creating total attributable cases for all ages per cancer type
  t2 %>%
    group_by(age_num) %>%
    mutate_at(vars(12:1011), funs(case_when(
      cancer_code %in% c(39, 40) ~ sum(., na.rm = T),
      TRUE ~ .
    ))) %>%
    bind_rows(
      t2 %>%
        group_by(age_num) %>%
        mutate_at(vars(12:1011), funs(case_when(
          cancer_code %in% c(39, 40) ~ sum(.[Main_analysis == 1], na.rm = T)
        ))) %>%
        mutate(cancer_label = case_when(
          cancer_code %in% c(39, 40) ~ paste0(cancer_label, "_M")
        )) %>% # add label to main all cancer results
        filter(cancer_code %in% c(39, 40))
    ) %>%
    group_by(age_num) %>%
    mutate_at(vars(12:1011), funs(
      case_when(
        cancer_code == 11 ~ sum(.[cancer_code == 11.1], na.rm = T),
        cancer_code == 6 ~ sum(.[cancer_code ==
                                   6.1], na.rm = T),
        TRUE ~ .
      )
    )) -> t3
  
  t3 %>% 
    bind_rows(t3 %>% 
                filter(cancer_code %in% c(3,5)) %>% 
                group_by(age_num) %>% 
                mutate(cancer_label = "Pharynx",
                       disease= "Pharynx",
                       cancer_code = 41,
                       cases = sum(cases),
                       total = sum(total)) %>%  
                mutate_at(vars(12:1011), funs(sum(., na.rm = T))) %>% 
                unique()) -> t3
  
  t4 <- as.data.table(t3)
  
  t4
  
})
sim_t3 <- list(sim_t3m, sim_t3f)
setwd("~/AlcPAF_os/simulations")
save(sim_t3, file = "sim_t3all.RData")
rm(sim_t3, sim_t2)
# replace aafs in missing countries with regional average aafs
SIM_SUBREG <- function(sim) {
  setwd("~/AlcPAF_os/simulations")
  simname <- paste0((sim), ".RData") # should load as sim_t3
  load(simname)
  
  # calculate AAF per simulation per cancer type for subregion to use for missing countries
  sim_subreg <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_t3[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      left_join(subregmatch) %>%
      group_by(subregion_label_2, cancer_label, age_num) %>%
      mutate_at(vars(12:1011), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
      mutate(cases = sum(cases, na.rm = T)) %>%
      ungroup() %>%
      dplyr::select(subregion_label_2,
                    sex,
                    age_num,
                    cancer_code,
                    cancer_label,
                    cases,
                    sim1:sim1000) %>%
      unique() %>%
      mutate_at(vars(7:1006), funs(. / cases))  -> t4 # calc aaf per age group
    
    t4 <- as.data.table(t4)
    
    t4
    
  })
  setwd("~/AlcPAF_os/simulations")
  save(sim_subreg, file = "sim_subreg.RData")
  
}
SIM_SUBREG(sim = "sim_t3all")

SIM_MISSING <- function(sim) { # need to re-calculate totals for missing countries
  setwd("~/AlcPAF_os/simulations")
  simname <- paste0((sim), ".RData")
  load(simname)
  load("sim_t3all.RData")
  sim_missing <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_t3[[SEX]][[REG]]#[[1]][[145]]
      cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and country number
      
      t2 <- as.data.table(t1)
      
      if (nrow(t2) > 0) {
        if (t2$country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)) {
          t2 %>%
            left_join(subregmatch) %>%
            dplyr::select(-sim1:-sim1000) %>%
            left_join(
              sim_subreg[[SEX]] %>% 
                dplyr::select(
                  subregion_label_2,
                  age_num,
                  cancer_label,
                  sim1:sim1000
                )
            ) %>% 
            dplyr::select(-subregion_label_2) %>%
            mutate_at(vars(13:1012), funs(. * cases))  %>% # calc aa cases per age
            group_by(age_num) %>%
            mutate_at(vars(13:1012), funs(
              case_when(
                cancer_label %in% c(
                  "All cancers but non-melanoma skin cancer_M",
                  "All cancers_M"
                ) ~ sum(.[cancer_code %in% c(1, 3, 5, 6.1, 8, 9, 11.1, 14, 20)], na.rm = T),
                cancer_label %in% c(
                  "All cancers but non-melanoma skin cancer",
                  "All cancers"
                ) ~ sum(.[cancer_code %in% c(1, 3, 5, 6.1, 8, 9, 11.1, 14, 20, 7, 13)], na.rm =
                          T),
                cancer_label == "Liver and intrahepatic bile ducts" ~ sum(.[cancer_code ==
                                                                              11.1], na.rm = T),
                cancer_label == "Oesophagus" ~ sum(.[cancer_code ==
                                                       6.1], na.rm = T),
                cancer_label == "Pharynx" ~ sum(.[cancer_code %in% c(3,5)], na.rm = T),
                TRUE ~ .
              )
            )) %>% 
            dplyr::select(-Main_analysis, -disease) -> t3
          
          t3
          
        } else {
          t2 %>% 
            dplyr::select(-Main_analysis, -disease) -> t2
          t2
        }
      } else {
        t2
      }
      
    })
    
  })
  setwd("~/AlcPAF_os/simulations")
  save(sim_missing, file = "sim_missing.RData")
  
}
SIM_MISSING(sim = "sim_subreg")
setwd("~/AlcPAF_os/simulations")
load("sim_missing.RData")
# calculate PAF per simulation per cancer type per country and sex including missing countries
sim_t4 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t2 <- sim_missing[[SEX]][[REG]]
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and country number
    
    t2 %>%
      group_by(cancer_label) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>%  # sum cases for all age groups
      dplyr::select(-paf_age_cancer,
                    -age_num,
                    -age,
                    -cases) %>% 
      unique() %>%
      mutate_at(vars(7:1006), funs(. / total)) -> t3
    
    t4 <- as.data.table(t3)
    
    t4
  })
  
})

# calculate 95% CIs of PAF per cancer type, country and sex
sim_t5 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t3 <- sim_t4[[SEX]][[REG]]
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and country number
    
    if (nrow(t3) > 0) {
      for (CAN in unique(t3$cancer_label)) {
        t3[cancer_label == CAN, LCI := quantile(t3[cancer_label == CAN, 7:1006], probs = 0.025, na.rm = TRUE)]
        t3[cancer_label == CAN, UCI := quantile(t3[cancer_label == CAN, 7:1006], probs = 0.975, na.rm = TRUE)]
        t3[cancer_label == CAN, mdn := quantile(t3[cancer_label == CAN, 7:1006], probs = 0.5, na.rm = TRUE)]
      }
      
    } else{
      t3 <-
        data.table(
          cancer_code = NA,
          ISO3 = NA,
          country_code = NA,
          sex = NA,
          cancer_label = NA,
          LCI = NA,
          UCI = NA,
          mdn = NA,
          total = NA
        )
    }
    
    t4 <-
      t3[, c(
        "cancer_code",
        "ISO3",
        "country_code",
        "sex",
        "cancer_label",
        "LCI",
        "UCI",
        "mdn",
        "total"
      )]
    t4
    
  })
  t1 <- do.call(rbind.data.frame, t)
})

# calculate cases per simulation per cancer type, country, sex and age
sim_t6 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t2 <- sim_missing[[SEX]][[REG]]
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and country number
    
    t2 %>%
      dplyr::select(-paf_age_cancer,-age,-cases) %>%
      unique() -> t3
    
    t4 <- as.data.table(t3)
    
    if (nrow(t4) > 0) {
      for (CAN in unique(t4$cancer_label)) {
        for (AGE in unique(t4$age_num)) {
          t4[(cancer_label == CAN &
                age_num == AGE), LCIage := quantile(t4[(cancer_label == CAN &
                                                          age_num == AGE), 8:1007], probs = 0.025, na.rm = TRUE)]
          t4[(cancer_label == CAN &
                age_num == AGE), UCIage := quantile(t4[(cancer_label == CAN &
                                                          age_num == AGE), 8:1007], probs = 0.975, na.rm = TRUE)]
          t4[(cancer_label == CAN &
                age_num == AGE), mdnage := quantile(t4[(cancer_label == CAN &
                                                          age_num == AGE), 8:1007], probs = 0.5, na.rm = TRUE)]
          
        }
      }
      
    } else {
      t4 <-
        data.table(
          cancer_code = NA,
          ISO3 = NA,
          country_code = NA,
          sex = NA,
          cancer_label = NA,
          age_num = NA,
          LCIage = NA,
          UCIage = NA,
          mdnage = NA,
          total = NA
        )
    }
    
    t4 <-
      t4[, c(
        "cancer_code",
        "ISO3",
        "country_code",
        "sex",
        "age_num",
        "cancer_label",
        "LCIage",
        "UCIage",
        "mdnage",
        "total"
      )]
    t4
  })
  t1 <- do.call(rbind.data.frame, t)
})

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
pafssim <-
  rbind(sim_t5[[1]], sim_t5[[2]])  # get 95% CIs of PAF per country
pafssima <-
  rbind(sim_t6[[1]], sim_t6[[2]]) # get 95% CIs of cases per age per country
pafssim %>%
  filter(!is.na(country_code)) %>%
  mutate(LCIc = LCI * total,
         UCIc = UCI * total,
         mdnc = mdn * total) %>%
  left_join(pafssima) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, ISO3, country_code, age_num) %>%
      mutate(
        LCIc = sum(LCIc, na.rm = T),
        UCIc = sum(UCIc, na.rm = T),
        mdnc = sum(mdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        LCI = LCIc / total,
        UCI = UCIc / total,
        mdn = mdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import other paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafs_inc_country_WCRF_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(
    pafssimt2 %>%
      dplyr::select(
        cancer_label,
        country_code,
        sex,
        age_num,
        LCI,
        UCI,
        mdn,
        LCIc,
        UCIc,
        mdnc,
        LCIage,
        UCIage,
        mdnage
      )
  ) %>%
  rename(
    totLCI = LCI,
    totUCI = UCI,
    totmdn = mdn,
    totLCIc = LCIc,
    totUCIc = UCIc,
    totmdnc = mdnc
  ) -> grpCI

# calculate ASRs for UCI and LCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) %>%
  left_join(pop_glob) -> grpCI

grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "py",
        group_by = c("country_code", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(country_code, sex, cancer_label, LCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "py",
        group_by = c("country_code", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(country_code, sex, cancer_label, UCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "aacases",
        var_py = "py",
        group_by = c("country_code", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(country_code, sex, cancer_label, asr)
  ) %>%
  dplyr::select(
    -aacases,
    -cases,
    -age,
    -paf_age_cancer,
    -age_num,
    -aaf,
    -py,
    -LCIage,
    -UCIage,
    -mdnage,
    -disease,
    -Main_analysis
  ) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_11.04.21.csv", row.names = FALSE)

rm(sim_subreg, sim_t, sim_t2, sim_t3, sim_t4, sim_t5, sim_t6)

###---- subregion CIs ----

# calculate aa cases per cancer type for subregion by age
sim_subreg2 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    left_join(subregmatch) %>%
    group_by(subregion_label_2, cancer_label, age_num) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
    mutate(cases = sum(cases, na.rm = T)) %>%
    ungroup() %>%
    dplyr::select(subregion_label_2,
                  sex,
                  age_num,
                  cancer_code,
                  cancer_label,
                  cases,
                  sim1:sim1000) %>%
    unique() -> t4
  
  t4 <- as.data.table(t4)
  
  for (REG in unique(t4$subregion_label_2)) {
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and subreg
    
    for (CAN in unique(t4$cancer_label)) {
      for (AGE in unique(t4$age_num)) {
        # calc 95%CIs for aa cases by age per subreg
        
        t4[(subregion_label_2 == REG &
              cancer_label == CAN &
              age_num == AGE), LCIage := quantile(t4[(subregion_label_2 == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
        t4[(subregion_label_2 == REG &
              cancer_label == CAN &
              age_num == AGE), UCIage := quantile(t4[(subregion_label_2 == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
        t4[(subregion_label_2 == REG &
              cancer_label == CAN &
              age_num == AGE), mdnage := quantile(t4[(subregion_label_2 == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
        
      }
    }
  }
  
  t5 <-
    t4[, c(
      "subregion_label_2",
      "age_num",
      "cancer_code",
      "sex",
      "cancer_label",
      "LCIage",
      "UCIage",
      "mdnage",
      "cases"
    )]
  t5
  
})

# calculate total aa cases and PAF per cancer type by subregion
sim_subreg3 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    left_join(subregmatch) %>%
    group_by(subregion_label_2, cancer_label) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each subreg
    group_by(subregion_label_2, cancer_label, age_num) %>%
    mutate(total = sum(total)) %>%
    ungroup() %>%
    dplyr::select(subregion_label_2,
                  sex,
                  cancer_code,
                  cancer_label,
                  total,
                  sim1:sim1000) %>%
    unique() %>%
    mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per subreg
  
  t4 <- as.data.table(t4)
  
  for (REG in unique(t4$subregion_label_2)) {
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and subreg
    
    for (CAN in unique(t4$cancer_label)) {
      t4[(subregion_label_2 == REG &
            cancer_label == CAN), totLCI := quantile(t4[(subregion_label_2 == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
      t4[(subregion_label_2 == REG &
            cancer_label == CAN), totUCI := quantile(t4[(subregion_label_2 == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
      t4[(subregion_label_2 == REG &
            cancer_label == CAN), totmdn := quantile(t4[(subregion_label_2 == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
      
    }
  }
  
  t5 <-
    t4[, c("subregion_label_2",
           "sex",
           "cancer_label",
           "total",
           "totLCI",
           "totUCI",
           "totmdn")]
  t5
  
})

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
pafssima <-
  rbind(sim_subreg2[[1]], sim_subreg2[[2]]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_subreg3[[1]], sim_subreg3[[2]])  # get 95% CIs of PAF per subreg

pafssim %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, subregion_label_2, age_num) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import subreg paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafs_subregion2_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(pafssimt2 %>%
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI
grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("subregion_label_2", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(subregion_label_2, sex, cancer_label, LCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("subregion_label_2", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(subregion_label_2, sex, cancer_label, UCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("subregion_label_2", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(subregion_label_2, sex, cancer_label, asr)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_subreg_11.04.21.csv", row.names = FALSE)

rm(sim_subreg2, sim_subreg3)

###---- HDI CIs ----

# calculate aa cases per cancer type for subregion by age
sim_hdi <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    left_join(hdimatch) %>%
    group_by(hdigrp, cancer_label, age_num) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
    mutate(cases = sum(cases, na.rm = T)) %>%
    ungroup() %>%
    dplyr::select(hdigrp,
                  sex,
                  age_num,
                  cancer_code,
                  cancer_label,
                  cases,
                  sim1:sim1000) %>%
    unique() -> t4
  
  t4 <- as.data.table(t4)
  
  for (REG in unique(t4$hdigrp)) {
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and subreg
    
    for (CAN in unique(t4$cancer_label)) {
      for (AGE in unique(t4$age_num)) {
        # calc 95%CIs for aa cases by age per subreg
        
        t4[(hdigrp == REG &
              cancer_label == CAN &
              age_num == AGE), LCIage := quantile(t4[(hdigrp == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
        t4[(hdigrp == REG &
              cancer_label == CAN &
              age_num == AGE), UCIage := quantile(t4[(hdigrp == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
        t4[(hdigrp == REG &
              cancer_label == CAN &
              age_num == AGE), mdnage := quantile(t4[(hdigrp == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
        
      }
    }
  }
  
  t5 <-
    t4[, c(
      "hdigrp",
      "age_num",
      "cancer_code",
      "sex",
      "cancer_label",
      "LCIage",
      "UCIage",
      "mdnage",
      "cases"
    )]
  t5
  
})

# calculate total aa cases and PAF per cancer type by subregion
sim_hdi2 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    left_join(hdimatch) %>%
    group_by(hdigrp, cancer_label) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each HDI
    group_by(hdigrp, cancer_label, age_num) %>%
    mutate(total = sum(total)) %>%
    ungroup() %>%
    dplyr::select(hdigrp, sex, cancer_code, cancer_label, total, sim1:sim1000) %>%
    unique() %>%
    mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per HDI
  
  t4 <- as.data.table(t4)
  
  for (REG in unique(t4$hdigrp)) {
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and HDI
    
    for (CAN in unique(t4$cancer_label)) {
      t4[(hdigrp == REG &
            cancer_label == CAN), totLCI := quantile(t4[(hdigrp == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
      t4[(hdigrp == REG &
            cancer_label == CAN), totUCI := quantile(t4[(hdigrp == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
      t4[(hdigrp == REG &
            cancer_label == CAN), totmdn := quantile(t4[(hdigrp == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
      
    }
  }
  
  t5 <-
    t4[, c("hdigrp",
           "sex",
           "cancer_label",
           "total",
           "totLCI",
           "totUCI",
           "totmdn")]
  t5
  
})

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
pafssima <-
  rbind(sim_hdi[[1]], sim_hdi[[2]]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_hdi2[[1]], sim_hdi2[[2]])  # get 95% CIs of PAF per subreg

pafssim %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, hdigrp, age_num) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import subreg paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafs_hdi_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(pafssimt2 %>%
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI
grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("hdigrp", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(hdigrp, sex, cancer_label, LCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("hdigrp", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(hdigrp, sex, cancer_label, UCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("hdigrp", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(hdigrp, sex, cancer_label, asr)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_hdi_11.04.21.csv", row.names = FALSE)

rm(sim_hdi, sim_hdi2)


###---- WHO CIs ----

# calculate aa cases per cancer type for subregion by age
sim_who <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    left_join(whomatch) %>%
    group_by(who_region_label, cancer_label, age_num) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
    mutate(cases = sum(cases, na.rm = T)) %>%
    ungroup() %>%
    dplyr::select(who_region_label,
                  sex,
                  age_num,
                  cancer_code,
                  cancer_label,
                  cases,
                  sim1:sim1000) %>%
    unique() -> t4
  
  t4 <- as.data.table(t4)
  
  for (REG in unique(t4$who_region_label)) {
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and subreg
    
    for (CAN in unique(t4$cancer_label)) {
      for (AGE in unique(t4$age_num)) {
        # calc 95%CIs for aa cases by age per subreg
        
        t4[(who_region_label == REG &
              cancer_label == CAN &
              age_num == AGE), LCIage := quantile(t4[(who_region_label == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
        t4[(who_region_label == REG &
              cancer_label == CAN &
              age_num == AGE), UCIage := quantile(t4[(who_region_label == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
        t4[(who_region_label == REG &
              cancer_label == CAN &
              age_num == AGE), mdnage := quantile(t4[(who_region_label == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
        
      }
    }
  }
  
  t5 <-
    t4[, c(
      "who_region_label",
      "age_num",
      "cancer_code",
      "sex",
      "cancer_label",
      "LCIage",
      "UCIage",
      "mdnage",
      "cases"
    )]
  t5
  
})

# calculate total aa cases and PAF per cancer type by subregion
sim_who2 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    left_join(whomatch) %>%
    group_by(who_region_label, cancer_label) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each HDI
    group_by(who_region_label, cancer_label, age_num) %>%
    mutate(total = sum(total)) %>%
    ungroup() %>%
    dplyr::select(who_region_label, sex, cancer_code, cancer_label, total, sim1:sim1000) %>%
    unique() %>%
    mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per HDI
  
  t4 <- as.data.table(t4)
  
  for (REG in unique(t4$who_region_label)) {
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and HDI
    
    for (CAN in unique(t4$cancer_label)) {
      t4[(who_region_label == REG &
            cancer_label == CAN), totLCI := quantile(t4[(who_region_label == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
      t4[(who_region_label == REG &
            cancer_label == CAN), totUCI := quantile(t4[(who_region_label == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
      t4[(who_region_label == REG &
            cancer_label == CAN), totmdn := quantile(t4[(who_region_label == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
      
    }
  }
  
  t5 <-
    t4[, c("who_region_label",
           "sex",
           "cancer_label",
           "total",
           "totLCI",
           "totUCI",
           "totmdn")]
  t5
  
})

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
pafssima <-
  rbind(sim_who[[1]], sim_who[[2]]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_who2[[1]], sim_who2[[2]])  # get 95% CIs of PAF per subreg

pafssim %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, who_region_label, age_num) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import subreg paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafs_who_reg_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(pafssimt2 %>%
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI
grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("who_region_label", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(who_region_label, sex, cancer_label, LCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("who_region_label", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(who_region_label, sex, cancer_label, UCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("who_region_label", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(who_region_label, sex, cancer_label, asr)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_who_11.04.21.csv", row.names = FALSE)

rm(sim_who, sim_who2)


###---- Continent CIs ----

# calculate aa cases per cancer type for subregion by age
sim_cont <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    left_join(contmatch) %>%
    group_by(region_label, cancer_label, age_num) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
    mutate(cases = sum(cases, na.rm = T)) %>%
    ungroup() %>%
    dplyr::select(region_label,
                  sex,
                  age_num,
                  cancer_code,
                  cancer_label,
                  cases,
                  sim1:sim1000) %>%
    unique() -> t4
  
  t4 <- as.data.table(t4)
  
  for (REG in unique(t4$region_label)) {
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and subreg
    
    for (CAN in unique(t4$cancer_label)) {
      for (AGE in unique(t4$age_num)) {
        # calc 95%CIs for aa cases by age per subreg
        
        t4[(region_label == REG &
              cancer_label == CAN &
              age_num == AGE), LCIage := quantile(t4[(region_label == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
        t4[(region_label == REG &
              cancer_label == CAN &
              age_num == AGE), UCIage := quantile(t4[(region_label == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
        t4[(region_label == REG &
              cancer_label == CAN &
              age_num == AGE), mdnage := quantile(t4[(region_label == REG &
                                                        cancer_label == CAN &
                                                        age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
        
      }
    }
  }
  
  t5 <-
    t4[, c(
      "region_label",
      "age_num",
      "cancer_code",
      "sex",
      "cancer_label",
      "LCIage",
      "UCIage",
      "mdnage",
      "cases"
    )]
  t5
  
})

# calculate total aa cases and PAF per cancer type by subregion
sim_cont2 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    left_join(contmatch) %>%
    group_by(region_label, cancer_label) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each HDI
    group_by(region_label, cancer_label, age_num) %>%
    mutate(total = sum(total)) %>%
    ungroup() %>%
    dplyr::select(region_label, sex, cancer_code, cancer_label, total, sim1:sim1000) %>%
    unique() %>%
    mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per HDI
  
  t4 <- as.data.table(t4)
  
  for (REG in unique(t4$region_label)) {
    cat("sex:", toString(SEX), ", reg:", toString(REG), "\n") # print sex and HDI
    
    for (CAN in unique(t4$cancer_label)) {
      t4[(region_label == REG &
            cancer_label == CAN), totLCI := quantile(t4[(region_label == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
      t4[(region_label == REG &
            cancer_label == CAN), totUCI := quantile(t4[(region_label == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
      t4[(region_label == REG &
            cancer_label == CAN), totmdn := quantile(t4[(region_label == REG &
                                                           cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
      
    }
  }
  
  t5 <-
    t4[, c("region_label",
           "sex",
           "cancer_label",
           "total",
           "totLCI",
           "totUCI",
           "totmdn")]
  t5
  
})

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
pafssima <-
  rbind(sim_cont[[1]], sim_cont[[2]]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_cont2[[1]], sim_cont2[[2]])  # get 95% CIs of PAF per subreg

pafssim %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, region_label, age_num) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import subreg paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafs_region_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(pafssimt2 %>%
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI
grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("region_label", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(region_label, sex, cancer_label, LCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("region_label", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(region_label, sex, cancer_label, UCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("region_label", "sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(region_label, sex, cancer_label, asr)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_continent_11.04.21.csv", row.names = FALSE)

rm(sim_cont, sim_cont2)

###---- world CIs ----

# calculate aa cases per cancer type for subregion by age
sim_w <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    mutate(worldgrp = "world") %>%
    group_by(cancer_label, age_num) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
    mutate(cases = sum(cases, na.rm = T)) %>%
    ungroup() %>%
    dplyr::select(worldgrp,
                  sex,
                  age_num,
                  cancer_code,
                  cancer_label,
                  cases,
                  sim1:sim1000) %>%
    unique() -> t4
  
  t4 <- as.data.table(t4)
  
  for (CAN in unique(t4$cancer_label)) {
    cat("sex:", toString(SEX), "; ", toString(CAN), "\n") # print sex and cancer
    
    for (AGE in unique(t4$age_num)) {
      # calc 95%CIs for aa cases by age
      
      t4[(cancer_label == CAN &
            age_num == AGE), LCIage := quantile(t4[(cancer_label == CAN &
                                                      age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
      t4[(cancer_label == CAN &
            age_num == AGE), UCIage := quantile(t4[(cancer_label == CAN &
                                                      age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
      t4[(cancer_label == CAN &
            age_num == AGE), mdnage := quantile(t4[(cancer_label == CAN &
                                                      age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
      
    }
  }
  
  t5 <-
    t4[, c(
      "worldgrp",
      "age_num",
      "cancer_code",
      "sex",
      "cancer_label",
      "LCIage",
      "UCIage",
      "mdnage",
      "cases"
    )]
  t5
  
})

# calculate total aa cases and PAF per cancer type by subregion
sim_w2 <- lapply(1:2, function(SEX) {
  t <- lapply(1:199, function(REG) {
    t1 <- sim_missing[[SEX]][[REG]]
    
    t2 <- as.data.table(t1)
    t2
  })
  
  t3 <- do.call(rbind.data.frame, t)
  
  t3 %>%
    mutate(worldgrp = "world") %>%
    group_by(cancer_label) %>%
    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each HDI
    group_by(worldgrp, cancer_label, age_num) %>%
    mutate(total = sum(total)) %>%
    ungroup() %>%
    dplyr::select(worldgrp, sex, cancer_code, cancer_label, total, sim1:sim1000) %>%
    unique() %>%
    mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per HDI
  
  t4 <- as.data.table(t4)
  
  for (CAN in unique(t4$cancer_label)) {
    cat("sex:", toString(SEX), "; ", toString(CAN), "\n") # print sex and cancer
    
    t4[(cancer_label == CAN), totLCI := quantile(t4[(cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
    t4[(cancer_label == CAN), totUCI := quantile(t4[(cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
    t4[(cancer_label == CAN), totmdn := quantile(t4[(cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
    
  }
  
  t5 <-
    t4[, c("worldgrp",
           "sex",
           "cancer_label",
           "total",
           "totLCI",
           "totUCI",
           "totmdn")]
  t5
  
})

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
pafssima <-
  rbind(sim_w[[1]], sim_w[[2]]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_w2[[1]], sim_w2[[2]])  # get 95% CIs of PAF per subreg

pafssim %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, age_num) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import world paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafs_world_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(pafssimt2 %>%
              dplyr::select(-worldgrp,-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI
grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(sex, cancer_label, LCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(sex, cancer_label, UCIasr)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("sex", "cancer_label"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(sex, cancer_label, asr)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_world_11.04.21.csv", row.names = FALSE)

rm(sim_w, sim_w2)

####---- Simulations and calculating CIs for levels calcs ####


####---- Import country simulations by level

# load simulations data
year <- 2020 # set year
# change working directory
setwd("~/AlcPAF_os/simulations")
# load simulated point estimates introducing age, sex and level variable
load(file = paste0("CANCER_1_1_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <-
  1
PAF_DATA_1_1 <-
  aaf_out_PAF
PAF_DATA_1_1 <-
  PAF_DATA_1_1[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_1_2_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <-
  1
PAF_DATA_1_2 <-
  aaf_out_PAF
PAF_DATA_1_2 <-
  PAF_DATA_1_2[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_1_3_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <-
  1
PAF_DATA_1_3 <-
  aaf_out_PAF
PAF_DATA_1_3 <-
  PAF_DATA_1_3[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_1_4_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <-
  1
PAF_DATA_1_4 <-
  aaf_out_PAF
PAF_DATA_1_4 <-
  PAF_DATA_1_4[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_1_5_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <-
  1
PAF_DATA_1_5 <-
  aaf_out_PAF
PAF_DATA_1_5 <-
  PAF_DATA_1_5[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_1_6_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <-
  1
PAF_DATA_1_6 <-
  aaf_out_PAF
PAF_DATA_1_6 <-
  PAF_DATA_1_6[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]

load(file = paste0("CANCER_2_1_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <-
  2
PAF_DATA_2_1 <-
  aaf_out_PAF
PAF_DATA_2_1 <-
  PAF_DATA_2_1[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_2_2_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <-
  2
PAF_DATA_2_2 <-
  aaf_out_PAF
PAF_DATA_2_2 <-
  PAF_DATA_2_2[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_2_3_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <-
  2
PAF_DATA_2_3 <-
  aaf_out_PAF
PAF_DATA_2_3 <-
  PAF_DATA_2_3[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_2_4_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <-
  2
PAF_DATA_2_4 <-
  aaf_out_PAF
PAF_DATA_2_4 <-
  PAF_DATA_2_4[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_2_5_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <-
  2
PAF_DATA_2_5 <-
  aaf_out_PAF
PAF_DATA_2_5 <-
  PAF_DATA_2_5[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]
load(file = paste0("CANCER_2_6_", (year - 10), "split_new.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <-
  2
PAF_DATA_2_6 <-
  aaf_out_PAF
PAF_DATA_2_6 <-
  PAF_DATA_2_6[str_detect(disease, "_light"), level := 1][str_detect(disease, "_mod"), level :=
                                                            2][str_detect(disease, "_heavy"), level := 3]

load(file = paste0("CANCER_1_1_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <-
  1
PAF_DATA_1_1_F <-
  aaf_out_PAF
PAF_DATA_1_1_F$disease <-
  paste0(PAF_DATA_1_1_F$disease, "_former")
PAF_DATA_1_1_F <- PAF_DATA_1_1_F[, level := 4]
load(file = paste0("CANCER_1_2_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <-
  1
PAF_DATA_1_2_F <-
  aaf_out_PAF
PAF_DATA_1_2_F$disease <-
  paste0(PAF_DATA_1_2_F$disease, "_former")
PAF_DATA_1_2_F <- PAF_DATA_1_2_F[, level := 4]
load(file = paste0("CANCER_1_3_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <-
  1
PAF_DATA_1_3_F <-
  aaf_out_PAF
PAF_DATA_1_3_F$disease <-
  paste0(PAF_DATA_1_3_F$disease, "_former")
PAF_DATA_1_3_F <- PAF_DATA_1_3_F[, level := 4]
load(file = paste0("CANCER_1_4_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <-
  1
PAF_DATA_1_4_F <-
  aaf_out_PAF
PAF_DATA_1_4_F$disease <-
  paste0(PAF_DATA_1_4_F$disease, "_former")
PAF_DATA_1_4_F <- PAF_DATA_1_4_F[, level := 4]
load(file = paste0("CANCER_1_5_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <-
  1
PAF_DATA_1_5_F <-
  aaf_out_PAF
PAF_DATA_1_5_F$disease <-
  paste0(PAF_DATA_1_5_F$disease, "_former")
PAF_DATA_1_5_F <- PAF_DATA_1_5_F[, level := 4]
load(file = paste0("CANCER_1_6_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <-
  1
PAF_DATA_1_6_F <-
  aaf_out_PAF
PAF_DATA_1_6_F$disease <-
  paste0(PAF_DATA_1_6_F$disease, "_former")
PAF_DATA_1_6_F <- PAF_DATA_1_6_F[, level := 4]

load(file = paste0("CANCER_2_1_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <-
  2
PAF_DATA_2_1_F <-
  aaf_out_PAF
PAF_DATA_2_1_F$disease <-
  paste0(PAF_DATA_2_1_F$disease, "_former")
PAF_DATA_2_1_F <- PAF_DATA_2_1_F[, level := 4]
load(file = paste0("CANCER_2_2_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <-
  2
PAF_DATA_2_2_F <-
  aaf_out_PAF
PAF_DATA_2_2_F$disease <-
  paste0(PAF_DATA_2_2_F$disease, "_former")
PAF_DATA_2_2_F <- PAF_DATA_2_2_F[, level := 4]
load(file = paste0("CANCER_2_3_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <-
  2
PAF_DATA_2_3_F <-
  aaf_out_PAF
PAF_DATA_2_3_F$disease <-
  paste0(PAF_DATA_2_3_F$disease, "_former")
PAF_DATA_2_3_F <- PAF_DATA_2_3_F[, level := 4]
load(file = paste0("CANCER_2_4_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <-
  2
PAF_DATA_2_4_F <-
  aaf_out_PAF
PAF_DATA_2_4_F$disease <-
  paste0(PAF_DATA_2_4_F$disease, "_former")
PAF_DATA_2_4_F <- PAF_DATA_2_4_F[, level := 4]
load(file = paste0("CANCER_2_5_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <-
  2
PAF_DATA_2_5_F <-
  aaf_out_PAF
PAF_DATA_2_5_F$disease <-
  paste0(PAF_DATA_2_5_F$disease, "_former")
PAF_DATA_2_5_F <- PAF_DATA_2_5_F[, level := 4]
load(file = paste0("CANCER_2_6_", (year - 10), "former.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <-
  2
PAF_DATA_2_6_F <-
  aaf_out_PAF
PAF_DATA_2_6_F$disease <-
  paste0(PAF_DATA_2_6_F$disease, "_former")
PAF_DATA_2_6_F <- PAF_DATA_2_6_F[, level := 4]

load(file = paste0("CANCER_1_1_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <-
  1
PAF_DATA_1_1_Fcd <-
  aaf_out_PAF
PAF_DATA_1_1_Fcd$disease <-
  paste0(PAF_DATA_1_1_Fcd$disease, "_formercd")
PAF_DATA_1_1_Fcd <- PAF_DATA_1_1_Fcd[, level := 5]
load(file = paste0("CANCER_1_2_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <-
  1
PAF_DATA_1_2_Fcd <-
  aaf_out_PAF
PAF_DATA_1_2_Fcd$disease <-
  paste0(PAF_DATA_1_2_Fcd$disease, "_formercd")
PAF_DATA_1_2_Fcd <- PAF_DATA_1_2_Fcd[, level := 5]
load(file = paste0("CANCER_1_3_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <-
  1
PAF_DATA_1_3_Fcd <-
  aaf_out_PAF
PAF_DATA_1_3_Fcd$disease <-
  paste0(PAF_DATA_1_3_Fcd$disease, "_formercd")
PAF_DATA_1_3_Fcd <- PAF_DATA_1_3_Fcd[, level := 5]
load(file = paste0("CANCER_1_4_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <-
  1
PAF_DATA_1_4_Fcd <-
  aaf_out_PAF
PAF_DATA_1_4_Fcd$disease <-
  paste0(PAF_DATA_1_4_Fcd$disease, "_formercd")
PAF_DATA_1_4_Fcd <- PAF_DATA_1_4_Fcd[, level := 5]
load(file = paste0("CANCER_1_5_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <-
  1
PAF_DATA_1_5_Fcd <-
  aaf_out_PAF
PAF_DATA_1_5_Fcd$disease <-
  paste0(PAF_DATA_1_5_Fcd$disease, "_formercd")
PAF_DATA_1_5_Fcd <- PAF_DATA_1_5_Fcd[, level := 5]
load(file = paste0("CANCER_1_6_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <-
  1
PAF_DATA_1_6_Fcd <-
  aaf_out_PAF
PAF_DATA_1_6_Fcd$disease <-
  paste0(PAF_DATA_1_6_Fcd$disease, "_formercd")
PAF_DATA_1_6_Fcd <- PAF_DATA_1_6_Fcd[, level := 5]

load(file = paste0("CANCER_2_1_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <-
  2
PAF_DATA_2_1_Fcd <-
  aaf_out_PAF
PAF_DATA_2_1_Fcd$disease <-
  paste0(PAF_DATA_2_1_Fcd$disease, "_formercd")
PAF_DATA_2_1_Fcd <- PAF_DATA_2_1_Fcd[, level := 5]
load(file = paste0("CANCER_2_2_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <-
  2
PAF_DATA_2_2_Fcd <-
  aaf_out_PAF
PAF_DATA_2_2_Fcd$disease <-
  paste0(PAF_DATA_2_2_Fcd$disease, "_formercd")
PAF_DATA_2_2_Fcd <- PAF_DATA_2_2_Fcd[, level := 5]
load(file = paste0("CANCER_2_3_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <-
  2
PAF_DATA_2_3_Fcd <-
  aaf_out_PAF
PAF_DATA_2_3_Fcd$disease <-
  paste0(PAF_DATA_2_3_Fcd$disease, "_formercd")
PAF_DATA_2_3_Fcd <- PAF_DATA_2_3_Fcd[, level := 5]
load(file = paste0("CANCER_2_4_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <-
  2
PAF_DATA_2_4_Fcd <-
  aaf_out_PAF
PAF_DATA_2_4_Fcd$disease <-
  paste0(PAF_DATA_2_4_Fcd$disease, "_formercd")
PAF_DATA_2_4_Fcd <- PAF_DATA_2_4_Fcd[, level := 5]
load(file = paste0("CANCER_2_5_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <-
  2
PAF_DATA_2_5_Fcd <-
  aaf_out_PAF
PAF_DATA_2_5_Fcd$disease <-
  paste0(PAF_DATA_2_5_Fcd$disease, "_formercd")
PAF_DATA_2_5_Fcd <- PAF_DATA_2_5_Fcd[, level := 5]
load(file = paste0("CANCER_2_6_", (year - 10), "formercd.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <-
  2
PAF_DATA_2_6_Fcd <-
  aaf_out_PAF
PAF_DATA_2_6_Fcd$disease <-
  paste0(PAF_DATA_2_6_Fcd$disease, "_formercd")
PAF_DATA_2_6_Fcd <- PAF_DATA_2_6_Fcd[, level := 5]


# add sims for former drinkers to sims for current drinkers for each sex/country/cancer/age
sim <-
  list(
    rbind(
      PAF_DATA_1_1,
      PAF_DATA_1_2,
      PAF_DATA_1_3,
      PAF_DATA_1_4,
      PAF_DATA_1_5,
      PAF_DATA_1_6,
      PAF_DATA_1_1_F,
      PAF_DATA_1_2_F,
      PAF_DATA_1_3_F,
      PAF_DATA_1_4_F,
      PAF_DATA_1_5_F,
      PAF_DATA_1_6_F,
      PAF_DATA_1_1_Fcd,
      PAF_DATA_1_2_Fcd,
      PAF_DATA_1_3_Fcd,
      PAF_DATA_1_4_Fcd,
      PAF_DATA_1_5_Fcd,
      PAF_DATA_1_6_Fcd
    ),
    rbind(
      PAF_DATA_2_1,
      PAF_DATA_2_2,
      PAF_DATA_2_3,
      PAF_DATA_2_4,
      PAF_DATA_2_5,
      PAF_DATA_2_6,
      PAF_DATA_2_1_F,
      PAF_DATA_2_2_F,
      PAF_DATA_2_3_F,
      PAF_DATA_2_4_F,
      PAF_DATA_2_5_F,
      PAF_DATA_2_6_F,
      PAF_DATA_2_1_Fcd,
      PAF_DATA_2_2_Fcd,
      PAF_DATA_2_3_Fcd,
      PAF_DATA_2_4_Fcd,
      PAF_DATA_2_5_Fcd,
      PAF_DATA_2_6_Fcd
    )
  )


rm(
  PAF_DATA_1_1,
  PAF_DATA_1_2,
  PAF_DATA_1_3,
  PAF_DATA_1_4,
  PAF_DATA_1_5,
  PAF_DATA_1_6,
  PAF_DATA_1_1_F,
  PAF_DATA_1_2_F,
  PAF_DATA_1_3_F,
  PAF_DATA_1_4_F,
  PAF_DATA_1_5_F,
  PAF_DATA_1_6_F,
  PAF_DATA_1_1_Fcd,
  PAF_DATA_1_2_Fcd,
  PAF_DATA_1_3_Fcd,
  PAF_DATA_1_4_Fcd,
  PAF_DATA_1_5_Fcd,
  PAF_DATA_1_6_Fcd,
  PAF_DATA_2_1,
  PAF_DATA_2_2,
  PAF_DATA_2_3,
  PAF_DATA_2_4,
  PAF_DATA_2_5,
  PAF_DATA_2_6,
  PAF_DATA_2_1_F,
  PAF_DATA_2_2_F,
  PAF_DATA_2_3_F,
  PAF_DATA_2_4_F,
  PAF_DATA_2_5_F,
  PAF_DATA_2_6_F,
  aaf_out_PAF,
  PAF_DATA_2_1_Fcd,
  PAF_DATA_2_2_Fcd,
  PAF_DATA_2_3_Fcd,
  PAF_DATA_2_4_Fcd,
  PAF_DATA_2_5_Fcd,
  PAF_DATA_2_6_Fcd
)


sim_t <- lapply(1:5, function(CAT) {
  lapply(1:2, function(SEX) {
    t <- merge(sim[[SEX]], Gcodes, by = "disease")
    t %>%
      bind_rows(t %>% 
                  filter(str_detect(disease,"Pharynx")) %>% # add rows for hypopharynx
                  mutate(cancer_code=5)) %>% 
      rename_at(vars(3:1002), function(x)
        paste0("sim", x)) %>%
      filter(!str_detect(disease, "_lin")) %>%
      filter(level == CAT) %>% 
      dplyr::select(-level) -> t
  })
})
rm(sim)

# merge simulations and globocan incidence datasets by sex, age, country and cancer type
sim_t. <- lapply(1:5, function(CAT) {
  t3 <- lapply(1:2, function(SEX) {
    t <-
      as.data.table(merge(
        glob_inc[sex == SEX, ],
        sim_t[[CAT]][[SEX]],
        by = c("cancer_code", "paf_age_cancer", "ISO3", "sex"),
        all = T,
        allow.cartesian = TRUE
      ))
    
#    t %>%
#      bind_rows(t %>% filter(cancer_label %in% c("Oesophagus","Liver and intrahepatic bile ducts","All cancers","All cancers but non-melanoma skin cancer") |
#                             #  is.na(disease) | # adds rows for ages 0-5 which have no disease
#                               (
#                                 cancer_code == 20 &
#                                   country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)
#                               ))) -> t
    
    t2 <- lapply(unique(t$ISO3), function(REG) {
      cat("categ:",
          toString(CAT),
          "sex:",
          toString(SEX),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      t3 <- t[ISO3 == REG, ]
    })
  })
})

for (CAT in 1:5) {
  dt <- sim_t.[[CAT]]
  filename <- paste0("simsplit_", CAT, ".RData")
  save(dt, file = filename)
}
rm(sim_t, sim_t.)


# calculate attributable cases per sex, country and simulation number, replace AAF with attributable cases
SIM_COUNTRY <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <- paste0("simsplit_", CAT, ".RData") # should load as dt
  load(simname)
  
  sim_t2 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t2 <- as.data.table(dt[[SEX]][[REG]]) 
      cat("categ:",
          toString(CAT),
          "sex:",
          toString(SEX),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (i in 12:1011) {
        t2[, as.numeric(i)] <- t2[, ..i] * t2$cases
        
      }
      
      # remove premenop results in older ages
      t2 <-
        t2[(cancer_code == 20 & age_num < 6), disease := "Breast"]
      t2 <-
        t2[(
          cancer_code == 20 &
            age_num > 10 &
            country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)
        ), disease := "Breast_postmen"]
      t2 <-
        t2[(
          cancer_code == 20 &
            age_num %in% c(6:10) &
            country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)
        ), disease := "Breast_premen"]
      t2 <-
        t2[!(cancer_code == 20 &
               str_detect(disease,"Breast_premen") & age_num > 10), ]
      # remove postmenop results in younger ages
      t2 <-
        t2[!(cancer_code == 20 &
               str_detect(disease,"Breast_postmen") & age_num <= 10), ]
      # make male breast cancer 0
      t2[cancer_code == 20 & sex == 1, c(12:1011)] <- 0
      
      t2 %>%
        group_by(age_num) %>%
        mutate_at(vars(12:1011), funs(case_when(
          cancer_code %in% c(39, 40) ~ sum(., na.rm = T),
          TRUE ~ .
        ))) %>%
        bind_rows(
          t2 %>%
            group_by(age_num) %>%
            mutate_at(vars(12:1011), funs(
              case_when(cancer_code %in% c(39, 40) ~ sum(.[Main_analysis == 1], na.rm =
                                                           T))
            )) %>%
            mutate(cancer_label = case_when(
              cancer_code %in% c(39, 40) ~ paste0(cancer_label, "_M")
            )) %>% # add label to main all cancer results
            filter(cancer_code %in% c(39, 40))
        ) %>%
        group_by(age_num) %>%
        mutate_at(vars(12:1011), funs(
          case_when(
            cancer_code == 11 ~ sum(.[cancer_code == 11.1], na.rm = T),
            cancer_code == 6 ~ sum(.[cancer_code ==
                                       6.1], na.rm = T),
            TRUE ~ .
          )
        )) %>% 
        dplyr::select(-Main_analysis, -disease) -> t3
      
      t3 %>% 
        bind_rows(t3 %>% 
                    filter(cancer_code %in% c(3,5)) %>% 
                    group_by(age_num) %>% 
                    mutate(cancer_label = "Pharynx",
                           cancer_code = 41,
                           cases = sum(cases),
                           total = sum(total)) %>%  
                    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% 
                    unique()) -> t3
      
      t3 <- as.data.table(t3)
      
      t3
    })
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_t2_", CAT, ".RData")
  save(sim_t2, file = filename)
}
for (category in 1:5) {
  SIM_COUNTRY(CAT = category)
}


# calculate AAF per simulation per cancer type for subregion to use for missing countries
SIM_SUBREG <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_t2_", CAT, ".RData") # should load as sim_t2
  load(simname)
  
  sim_subreg <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_t2[[SEX]][[REG]]
      cat("categ:",
          toString(CAT),
          "sex:",
          toString(SEX),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>% 
      left_join(subregmatch) %>%
      group_by(subregion_label_2, cancer_label, age_num) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
      mutate(cases = sum(cases, na.rm = T)) %>%
      ungroup() %>%
      dplyr::select(
        subregion_label_2,
        sex,
        age_num,
        cancer_code,
        cancer_label,
        cases,
        sim1:sim1000
      ) %>%
      unique() %>%
      mutate_at(vars(7:1006), funs(. / cases))  -> t4 # calc aaf per age group
    
    t4 <- as.data.table(t4)
    
    t4
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_subreg_", CAT, ".RData")
  save(sim_subreg, file = filename)
}
for (category in 1:5) {
  SIM_SUBREG(CAT = category)
}

# replace aafs in missing countries with regional average aafs
SIM_MISSING <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_subreg_", CAT, ".RData") # should load as sim_subreg
  load(simname)
  simname <-
    paste0("sim_t2_", CAT, ".RData") # should load as sim_t2
  load(simname)
  
  sim_missing <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_t2[[SEX]][[REG]]
      cat("categ:",
          toString(CAT),
          "sex:",
          toString(SEX),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      t2 <- as.data.table(t1)
      
      if (nrow(t2) > 0) {
        if (t2$country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)) {
          t2 %>%
            left_join(subregmatch) %>%
            dplyr::select(-sim1:-sim1000) %>% 
            left_join(
              sim_subreg[[SEX]] %>%
                dplyr::select(
                  subregion_label_2,
                  age_num,
                  cancer_label,
                  sim1:sim1000
                )
            ) %>%
            dplyr::select(-subregion_label_2) %>% 
            mutate_at(vars(11:1010), funs(. * cases))  %>%  # calc aa cases per age
            group_by(age_num) %>%
            mutate_at(vars(11:1010), funs(
              case_when(
                cancer_label %in% c(
                  "All cancers but non-melanoma skin cancer_M",
                  "All cancers_M"
                ) ~ sum(.[cancer_code %in% c(1, 3, 5, 6.1, 8, 9, 11.1, 14, 20)], na.rm = T),
                cancer_label %in% c(
                  "All cancers but non-melanoma skin cancer",
                  "All cancers"
                ) ~ sum(.[cancer_code %in% c(1, 3, 5, 6.1, 8, 9, 11.1, 14, 20, 7, 13)], na.rm =
                          T),
                cancer_label == "Liver and intrahepatic bile ducts" ~ sum(.[cancer_code ==
                                                                              11.1], na.rm = T),
                cancer_label == "Oesophagus" ~ sum(.[cancer_code ==
                                                       6.1], na.rm = T),
                cancer_label == "Pharynx" ~ sum(.[cancer_code %in% c(3,5)], na.rm = T),
                TRUE ~ .
              )
            )) -> t3
          
          t3
          
        } else {
          t2
        }
      } else {
        t2
      }
      
    })
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_missing_", CAT, ".RData")
  save(sim_missing, file = filename)
}
for (category in 1:5) {
  SIM_MISSING(CAT = category)
}

###---- Subregion CIs ----
# calculate AAF per simulation per cancer type for subregion to use for missing countries
SIM_SUBREG2 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for subregion by age
  sim_subreg2 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      left_join(subregmatch) %>%
      group_by(subregion_label_2, cancer_label, age_num) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
      mutate(cases = sum(cases, na.rm = T)) %>%
      ungroup() %>%
      dplyr::select(
        subregion_label_2,
        sex,
        age_num,
        cancer_code,
        cancer_label,
        cases,
        sim1:sim1000
      ) %>%
      unique() -> t4
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$subregion_label_2)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        for (AGE in unique(t4$age_num)) {
          # calc 95%CIs for aa cases by age per subreg
          
          
          t4[(subregion_label_2 == REG &
                cancer_label == CAN &
                age_num == AGE), LCIage := quantile(t4[(subregion_label_2 == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
          t4[(subregion_label_2 == REG &
                cancer_label == CAN &
                age_num == AGE), UCIage := quantile(t4[(subregion_label_2 == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
          t4[(subregion_label_2 == REG &
                cancer_label == CAN &
                age_num == AGE), mdnage := quantile(t4[(subregion_label_2 == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
          
        }
      }
    }
    
    t5 <-
      t4[, c(
        "subregion_label_2",
        "age_num",
        "cancer_code",
        "sex",
        "cancer_label",
        "LCIage",
        "UCIage",
        "mdnage",
        "cases"
      )]
    t5
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_subreg2_", CAT, ".RData")
  save(sim_subreg2, file = filename)
}
for (category in 1:5) {
  SIM_SUBREG2(CAT = category)
}

# calculate total aa cases and PAF per cancer type by subregion
# calculate AAF per simulation per cancer type for subregion to use for missing countries
SIM_SUBREG3 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for subregion by age
  sim_subreg3 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      left_join(subregmatch) %>% 
      group_by(subregion_label_2, cancer_label) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each subreg
      group_by(subregion_label_2, cancer_label, age_num) %>%
      mutate(total = sum(total)) %>%
      ungroup() %>%
      dplyr::select(subregion_label_2,
                    sex,
                    cancer_code,
                    cancer_label,
                    total,
                    sim1:sim1000) %>%
      unique() %>%
      mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per subreg
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$subregion_label_2)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        t4[(subregion_label_2 == REG &
              cancer_label == CAN), totLCI := quantile(t4[(subregion_label_2 == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
        t4[(subregion_label_2 == REG &
              cancer_label == CAN), totUCI := quantile(t4[(subregion_label_2 == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
        t4[(subregion_label_2 == REG &
              cancer_label == CAN), totmdn := quantile(t4[(subregion_label_2 == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
        
        
      }
    }
    
    t5 <-
      t4[, c(
        "subregion_label_2",
        "sex",
        "cancer_label",
        "total",
        "totLCI",
        "totUCI",
        "totmdn"
      )]
    t5
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_subreg3_", CAT, ".RData")
  save(sim_subreg3, file = filename)
}
for (category in 1:5) {
  SIM_SUBREG3(CAT = category)
}

# for each level create persons totals - load sim_subreg2 and sim_subreg3

setwd("~/AlcPAF_os/simulations")
load("sim_subreg2_1.RData"); sim_subreg2_1 <- sim_subreg2
load("sim_subreg2_2.RData"); sim_subreg2_2 <- sim_subreg2
load("sim_subreg2_3.RData"); sim_subreg2_3 <- sim_subreg2
load("sim_subreg2_4.RData"); sim_subreg2_4 <- sim_subreg2
load("sim_subreg2_5.RData"); sim_subreg2_5 <- sim_subreg2

load("sim_subreg3_1.RData"); sim_subreg3_1 <- sim_subreg3
load("sim_subreg3_2.RData"); sim_subreg3_2 <- sim_subreg3
load("sim_subreg3_3.RData"); sim_subreg3_3 <- sim_subreg3
load("sim_subreg3_4.RData"); sim_subreg3_4 <- sim_subreg3
load("sim_subreg3_5.RData"); sim_subreg3_5 <- sim_subreg3

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
# add level variable back in
pafssima <-
  rbind(sim_subreg2_1[[1]][,level:=1], sim_subreg2_1[[2]][,level:=1],
        sim_subreg2_2[[1]][,level:=2], sim_subreg2_2[[2]][,level:=2],
        sim_subreg2_3[[1]][,level:=3], sim_subreg2_3[[2]][,level:=3],
        sim_subreg2_4[[1]][,level:=4], sim_subreg2_4[[2]][,level:=4],
        sim_subreg2_5[[1]][,level:=5], sim_subreg2_5[[2]][,level:=5]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_subreg3_1[[1]][,level:=1], sim_subreg3_1[[2]][,level:=1],
        sim_subreg3_2[[1]][,level:=2], sim_subreg3_2[[2]][,level:=2],
        sim_subreg3_3[[1]][,level:=3], sim_subreg3_3[[2]][,level:=3],
        sim_subreg3_4[[1]][,level:=4], sim_subreg3_4[[2]][,level:=4],
        sim_subreg3_5[[1]][,level:=5], sim_subreg3_5[[2]][,level:=5])  # get 95% CIs of PAF per subreg

pafssim %>%
  group_by(cancer_label, subregion_label_2, level, sex) %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, subregion_label_2, age_num, level) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import subreg paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafsplit_subregion2_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(pafssimt2 %>%
              filter(!is.na(subregion_label_2)) %>% 
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI

grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("subregion_label_2", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(subregion_label_2, sex, cancer_label, LCIasr,level)
  ) %>% 
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("subregion_label_2", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(subregion_label_2, sex, cancer_label, UCIasr,level)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("subregion_label_2", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(subregion_label_2, sex, cancer_label, asr,level)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_subregsplit_11.04.21.csv", row.names =
            FALSE)

rm(sim_subreg2_1,sim_subreg2_2,sim_subreg2_3,sim_subreg2_4,sim_subreg2_5,
   sim_subreg3_1,sim_subreg3_2,sim_subreg3_3,sim_subreg3_4,sim_subreg3_5,
   sim_subreg2, sim_subreg3)

###---- HDI CIs ----

SIM_HDI2 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for hdi by age
  sim_hdi2 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      left_join(hdimatch) %>%
      group_by(hdigrp, cancer_label, age_num) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
      mutate(cases = sum(cases, na.rm = T)) %>%
      ungroup() %>%
      dplyr::select(
        hdigrp,
        sex,
        age_num,
        cancer_code,
        cancer_label,
        cases,
        sim1:sim1000
      ) %>%
      unique() -> t4
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$hdigrp)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        for (AGE in unique(t4$age_num)) {
          # calc 95%CIs for aa cases by age per subreg
          
          
          t4[(hdigrp == REG &
                cancer_label == CAN &
                age_num == AGE), LCIage := quantile(t4[(hdigrp == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
          t4[(hdigrp == REG &
                cancer_label == CAN &
                age_num == AGE), UCIage := quantile(t4[(hdigrp == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
          t4[(hdigrp == REG &
                cancer_label == CAN &
                age_num == AGE), mdnage := quantile(t4[(hdigrp == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
          
        }
      }
    }
    
    t5 <-
      t4[, c(
        "hdigrp",
        "age_num",
        "cancer_code",
        "sex",
        "cancer_label",
        "LCIage",
        "UCIage",
        "mdnage",
        "cases"
      )]
    t5
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_hdi2_", CAT, ".RData")
  save(sim_hdi2, file = filename)
}
for (category in 1:5) {
  SIM_HDI2(CAT = category)
}

# calculate total aa cases and PAF per cancer type by hdi
# calculate AAF per simulation per cancer type for hdi to use for missing countries
SIM_HDI3 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for hdi by age
  sim_hdi3 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      left_join(hdimatch) %>% 
      group_by(hdigrp, cancer_label) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each subreg
      group_by(hdigrp, cancer_label, age_num) %>%
      mutate(total = sum(total)) %>%
      ungroup() %>%
      dplyr::select(hdigrp,
                    sex,
                    cancer_code,
                    cancer_label,
                    total,
                    sim1:sim1000) %>%
      unique() %>%
      mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per subreg
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$hdigrp)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        t4[(hdigrp == REG &
              cancer_label == CAN), totLCI := quantile(t4[(hdigrp == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
        t4[(hdigrp == REG &
              cancer_label == CAN), totUCI := quantile(t4[(hdigrp == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
        t4[(hdigrp == REG &
              cancer_label == CAN), totmdn := quantile(t4[(hdigrp == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
        
        
      }
    }
    
    t5 <-
      t4[, c(
        "hdigrp",
        "sex",
        "cancer_label",
        "total",
        "totLCI",
        "totUCI",
        "totmdn"
      )]
    t5
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_hdi3_", CAT, ".RData")
  save(sim_hdi3, file = filename)
}
for (category in 1:5) {
  SIM_HDI3(CAT = category)
}

# for each level create persons totals - load sim_subreg2 and sim_subreg3

setwd("~/AlcPAF_os/simulations")
load("sim_hdi2_1.RData"); sim_hdi2_1 <- sim_hdi2
load("sim_hdi2_2.RData"); sim_hdi2_2 <- sim_hdi2
load("sim_hdi2_3.RData"); sim_hdi2_3 <- sim_hdi2
load("sim_hdi2_4.RData"); sim_hdi2_4 <- sim_hdi2
load("sim_hdi2_5.RData"); sim_hdi2_5 <- sim_hdi2

load("sim_hdi3_1.RData"); sim_hdi3_1 <- sim_hdi3
load("sim_hdi3_2.RData"); sim_hdi3_2 <- sim_hdi3
load("sim_hdi3_3.RData"); sim_hdi3_3 <- sim_hdi3
load("sim_hdi3_4.RData"); sim_hdi3_4 <- sim_hdi3
load("sim_hdi3_5.RData"); sim_hdi3_5 <- sim_hdi3

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
# add level variable back in
pafssima <-
  rbind(sim_hdi2_1[[1]][,level:=1], sim_hdi2_1[[2]][,level:=1],
        sim_hdi2_2[[1]][,level:=2], sim_hdi2_2[[2]][,level:=2],
        sim_hdi2_3[[1]][,level:=3], sim_hdi2_3[[2]][,level:=3],
        sim_hdi2_4[[1]][,level:=4], sim_hdi2_4[[2]][,level:=4],
        sim_hdi2_5[[1]][,level:=5], sim_hdi2_5[[2]][,level:=5]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_hdi3_1[[1]][,level:=1], sim_hdi3_1[[2]][,level:=1],
        sim_hdi3_2[[1]][,level:=2], sim_hdi3_2[[2]][,level:=2],
        sim_hdi3_3[[1]][,level:=3], sim_hdi3_3[[2]][,level:=3],
        sim_hdi3_4[[1]][,level:=4], sim_hdi3_4[[2]][,level:=4],
        sim_hdi3_5[[1]][,level:=5], sim_hdi3_5[[2]][,level:=5])  # get 95% CIs of PAF per subreg

pafssim %>%
  group_by(cancer_label, hdigrp, level, sex) %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, hdigrp, age_num, level) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import subreg paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafsplit_hdi_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(pafssimt2 %>%
              filter(!is.na(hdigrp)) %>% 
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI

grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("hdigrp", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(hdigrp, sex, cancer_label, LCIasr,level)
  ) %>% 
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("hdigrp", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(hdigrp, sex, cancer_label, UCIasr,level)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("hdigrp", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(hdigrp, sex, cancer_label, asr,level)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_hdisplit_11.04.21.csv", row.names =
            FALSE)

rm(sim_hdi2_1,sim_hdi2_2,sim_hdi2_3,sim_hdi2_4,sim_hdi2_5,
   sim_hdi3_1,sim_hdi3_2,sim_hdi3_3,sim_hdi3_4,sim_hdi3_5,
   sim_hdi2, sim_hdi3)

###---- WHO region CIs ----

SIM_WHO2 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for who by age
  sim_who2 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      left_join(whomatch) %>%
      group_by(who_region_label, cancer_label, age_num) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
      mutate(cases = sum(cases, na.rm = T)) %>%
      ungroup() %>%
      dplyr::select(
        who_region_label,
        sex,
        age_num,
        cancer_code,
        cancer_label,
        cases,
        sim1:sim1000
      ) %>%
      unique() -> t4
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$who_region_label)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        for (AGE in unique(t4$age_num)) {
          # calc 95%CIs for aa cases by age per subreg
          
          
          t4[(who_region_label == REG &
                cancer_label == CAN &
                age_num == AGE), LCIage := quantile(t4[(who_region_label == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
          t4[(who_region_label == REG &
                cancer_label == CAN &
                age_num == AGE), UCIage := quantile(t4[(who_region_label == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
          t4[(who_region_label == REG &
                cancer_label == CAN &
                age_num == AGE), mdnage := quantile(t4[(who_region_label == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
          
        }
      }
    }
    
    t5 <-
      t4[, c(
        "who_region_label",
        "age_num",
        "cancer_code",
        "sex",
        "cancer_label",
        "LCIage",
        "UCIage",
        "mdnage",
        "cases"
      )]
    t5
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_who2_", CAT, ".RData")
  save(sim_who2, file = filename)
}
for (category in 1:5) {
  SIM_WHO2(CAT = category)
}

# calculate total aa cases and PAF per cancer type by who
# calculate AAF per simulation per cancer type for who to use for missing countries
SIM_WHO3 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for who by age
  sim_who3 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      left_join(whomatch) %>% 
      group_by(who_region_label, cancer_label) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each subreg
      group_by(who_region_label, cancer_label, age_num) %>%
      mutate(total = sum(total)) %>%
      ungroup() %>%
      dplyr::select(who_region_label,
                    sex,
                    cancer_code,
                    cancer_label,
                    total,
                    sim1:sim1000) %>%
      unique() %>%
      mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per subreg
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$who_region_label)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        t4[(who_region_label == REG &
              cancer_label == CAN), totLCI := quantile(t4[(who_region_label == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
        t4[(who_region_label == REG &
              cancer_label == CAN), totUCI := quantile(t4[(who_region_label == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
        t4[(who_region_label == REG &
              cancer_label == CAN), totmdn := quantile(t4[(who_region_label == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
        
        
      }
    }
    
    t5 <-
      t4[, c(
        "who_region_label",
        "sex",
        "cancer_label",
        "total",
        "totLCI",
        "totUCI",
        "totmdn"
      )]
    t5
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_who3_", CAT, ".RData")
  save(sim_who3, file = filename)
}
for (category in 1:5) {
  SIM_WHO3(CAT = category)
}

# for each level create persons totals - load sim_subreg2 and sim_subreg3

setwd("~/AlcPAF_os/simulations")
load("sim_who2_1.RData"); sim_who2_1 <- sim_who2
load("sim_who2_2.RData"); sim_who2_2 <- sim_who2
load("sim_who2_3.RData"); sim_who2_3 <- sim_who2
load("sim_who2_4.RData"); sim_who2_4 <- sim_who2
load("sim_who2_5.RData"); sim_who2_5 <- sim_who2

load("sim_who3_1.RData"); sim_who3_1 <- sim_who3
load("sim_who3_2.RData"); sim_who3_2 <- sim_who3
load("sim_who3_3.RData"); sim_who3_3 <- sim_who3
load("sim_who3_4.RData"); sim_who3_4 <- sim_who3
load("sim_who3_5.RData"); sim_who3_5 <- sim_who3

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
# add level variable back in
pafssima <-
  rbind(sim_who2_1[[1]][,level:=1], sim_who2_1[[2]][,level:=1],
        sim_who2_2[[1]][,level:=2], sim_who2_2[[2]][,level:=2],
        sim_who2_3[[1]][,level:=3], sim_who2_3[[2]][,level:=3],
        sim_who2_4[[1]][,level:=4], sim_who2_4[[2]][,level:=4],
        sim_who2_5[[1]][,level:=5], sim_who2_5[[2]][,level:=5]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_who3_1[[1]][,level:=1], sim_who3_1[[2]][,level:=1],
        sim_who3_2[[1]][,level:=2], sim_who3_2[[2]][,level:=2],
        sim_who3_3[[1]][,level:=3], sim_who3_3[[2]][,level:=3],
        sim_who3_4[[1]][,level:=4], sim_who3_4[[2]][,level:=4],
        sim_who3_5[[1]][,level:=5], sim_who3_5[[2]][,level:=5])  # get 95% CIs of PAF per subreg

pafssim %>%
  group_by(cancer_label, who_region_label, level, sex) %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, who_region_label, age_num, level) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import subreg paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafsplit_who_reg_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(pafssimt2 %>%
              filter(!is.na(who_region_label)) %>% 
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI

grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("who_region_label", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(who_region_label, sex, cancer_label, LCIasr,level)
  ) %>% 
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("who_region_label", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(who_region_label, sex, cancer_label, UCIasr,level)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("who_region_label", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(who_region_label, sex, cancer_label, asr,level)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_whosplit_11.04.21.csv", row.names =
            FALSE)

rm(sim_who2_1,sim_who2_2,sim_who2_3,sim_who2_4,sim_who2_5,
   sim_who3_1,sim_who3_2,sim_who3_3,sim_who3_4,sim_who3_5,
   sim_who2, sim_who3)


###---- Continent CIs ----
SIM_CONT2 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for contion by age
  sim_cont2 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      left_join(contmatch) %>%
      group_by(region_label, cancer_label, age_num) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & cont
      mutate(cases = sum(cases, na.rm = T)) %>%
      ungroup() %>%
      dplyr::select(
        region_label,
        sex,
        age_num,
        cancer_code,
        cancer_label,
        cases,
        sim1:sim1000
      ) %>%
      unique() -> t4
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$region_label)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        for (AGE in unique(t4$age_num)) {
          # calc 95%CIs for aa cases by age per cont
          
          
          t4[(region_label == REG &
                cancer_label == CAN &
                age_num == AGE), LCIage := quantile(t4[(region_label == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
          t4[(region_label == REG &
                cancer_label == CAN &
                age_num == AGE), UCIage := quantile(t4[(region_label == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
          t4[(region_label == REG &
                cancer_label == CAN &
                age_num == AGE), mdnage := quantile(t4[(region_label == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
          
        }
      }
    }
    
    t5 <-
      t4[, c(
        "region_label",
        "age_num",
        "cancer_code",
        "sex",
        "cancer_label",
        "LCIage",
        "UCIage",
        "mdnage",
        "cases"
      )]
    t5
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_cont2_", CAT, ".RData")
  save(sim_cont2, file = filename)
}
for (category in 1:5) {
  SIM_CONT2(CAT = category)
}

# calculate total aa cases and PAF per cancer type by contion
# calculate AAF per simulation per cancer type for contion to use for missing countries
SIM_CONT3 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for contion by age
  sim_cont3 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      left_join(contmatch) %>% 
      group_by(region_label, cancer_label) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each cont
      group_by(region_label, cancer_label, age_num) %>%
      mutate(total = sum(total)) %>%
      ungroup() %>%
      dplyr::select(region_label,
                    sex,
                    cancer_code,
                    cancer_label,
                    total,
                    sim1:sim1000) %>%
      unique() %>%
      mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per cont
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$region_label)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        t4[(region_label == REG &
              cancer_label == CAN), totLCI := quantile(t4[(region_label == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
        t4[(region_label == REG &
              cancer_label == CAN), totUCI := quantile(t4[(region_label == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
        t4[(region_label == REG &
              cancer_label == CAN), totmdn := quantile(t4[(region_label == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
        
        
      }
    }
    
    t5 <-
      t4[, c(
        "region_label",
        "sex",
        "cancer_label",
        "total",
        "totLCI",
        "totUCI",
        "totmdn"
      )]
    t5
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_cont3_", CAT, ".RData")
  save(sim_cont3, file = filename)
}
for (category in 1:5) {
  SIM_CONT3(CAT = category)
}

# for each level create persons totals - load sim_cont2 and sim_cont3

setwd("~/AlcPAF_os/simulations")
load("sim_cont2_1.RData"); sim_cont2_1 <- sim_cont2
load("sim_cont2_2.RData"); sim_cont2_2 <- sim_cont2
load("sim_cont2_3.RData"); sim_cont2_3 <- sim_cont2
load("sim_cont2_4.RData"); sim_cont2_4 <- sim_cont2
load("sim_cont2_5.RData"); sim_cont2_5 <- sim_cont2

load("sim_cont3_1.RData"); sim_cont3_1 <- sim_cont3
load("sim_cont3_2.RData"); sim_cont3_2 <- sim_cont3
load("sim_cont3_3.RData"); sim_cont3_3 <- sim_cont3
load("sim_cont3_4.RData"); sim_cont3_4 <- sim_cont3
load("sim_cont3_5.RData"); sim_cont3_5 <- sim_cont3

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
# add level variable back in
pafssima <-
  rbind(sim_cont2_1[[1]][,level:=1], sim_cont2_1[[2]][,level:=1],
        sim_cont2_2[[1]][,level:=2], sim_cont2_2[[2]][,level:=2],
        sim_cont2_3[[1]][,level:=3], sim_cont2_3[[2]][,level:=3],
        sim_cont2_4[[1]][,level:=4], sim_cont2_4[[2]][,level:=4],
        sim_cont2_5[[1]][,level:=5], sim_cont2_5[[2]][,level:=5]) # get 95% CIs of cases per age per cont
pafssim <-
  rbind(sim_cont3_1[[1]][,level:=1], sim_cont3_1[[2]][,level:=1],
        sim_cont3_2[[1]][,level:=2], sim_cont3_2[[2]][,level:=2],
        sim_cont3_3[[1]][,level:=3], sim_cont3_3[[2]][,level:=3],
        sim_cont3_4[[1]][,level:=4], sim_cont3_4[[2]][,level:=4],
        sim_cont3_5[[1]][,level:=5], sim_cont3_5[[2]][,level:=5])  # get 95% CIs of PAF per cont

pafssim %>%
  group_by(cancer_label, region_label, level, sex) %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, region_label, age_num, level) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import cont paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafsplit_region_11.04.21.csv")) # 2020 data
pafs %>%
  full_join(pafssimt2 %>%
              filter(!is.na(region_label)) %>% 
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI

grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("region_label", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(region_label, sex, cancer_label, LCIasr,level)
  ) %>% 
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("region_label", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(region_label, sex, cancer_label, UCIasr,level)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("region_label", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(region_label, sex, cancer_label, asr,level)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_contsplit_11.04.21.csv", row.names =
            FALSE)

rm(sim_cont2_1,sim_cont2_2,sim_cont2_3,sim_cont2_4,sim_cont2_5,
   sim_cont3_1,sim_cont3_2,sim_cont3_3,sim_cont3_4,sim_cont3_5,
   sim_cont2, sim_cont3)

###---- world CIs ----

SIM_WORLD2 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for world by age
  sim_w2 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      mutate(worldgrp = "world") %>%
      group_by(worldgrp, cancer_label, age_num) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
      mutate(cases = sum(cases, na.rm = T)) %>%
      ungroup() %>%
      dplyr::select(
        worldgrp,
        sex,
        age_num,
        cancer_code,
        cancer_label,
        cases,
        sim1:sim1000
      ) %>%
      unique() -> t4
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$worldgrp)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        for (AGE in unique(t4$age_num)) {
          # calc 95%CIs for aa cases by age per subreg
          
          
          t4[(worldgrp == REG &
                cancer_label == CAN &
                age_num == AGE), LCIage := quantile(t4[(worldgrp == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
          t4[(worldgrp == REG &
                cancer_label == CAN &
                age_num == AGE), UCIage := quantile(t4[(worldgrp == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
          t4[(worldgrp == REG &
                cancer_label == CAN &
                age_num == AGE), mdnage := quantile(t4[(worldgrp == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
          
        }
      }
    }
    
    t5 <-
      t4[, c(
        "worldgrp",
        "age_num",
        "cancer_code",
        "sex",
        "cancer_label",
        "LCIage",
        "UCIage",
        "mdnage",
        "cases"
      )]
    t5
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_w2_", CAT, ".RData")
  save(sim_w2, file = filename)
}
for (category in 1:5) {
  SIM_WORLD2(CAT = category)
}

# calculate total aa cases and PAF per cancer type by world
# calculate AAF per simulation per cancer type for world to use for missing countries
SIM_WORLD3 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for world by age
  sim_w3 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      mutate(worldgrp = "world") %>%
      group_by(worldgrp, cancer_label) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each subreg
      group_by(worldgrp, cancer_label, age_num) %>%
      mutate(total = sum(total)) %>%
      ungroup() %>%
      dplyr::select(worldgrp,
                    sex,
                    cancer_code,
                    cancer_label,
                    total,
                    sim1:sim1000) %>%
      unique() %>%
      mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per subreg
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$worldgrp)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        t4[(worldgrp == REG &
              cancer_label == CAN), totLCI := quantile(t4[(worldgrp == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
        t4[(worldgrp == REG &
              cancer_label == CAN), totUCI := quantile(t4[(worldgrp == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
        t4[(worldgrp == REG &
              cancer_label == CAN), totmdn := quantile(t4[(worldgrp == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
        
        
      }
    }
    
    t5 <-
      t4[, c(
        "worldgrp",
        "sex",
        "cancer_label",
        "total",
        "totLCI",
        "totUCI",
        "totmdn"
      )]
    t5
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_w3_", CAT, ".RData")
  save(sim_w3, file = filename)
}
for (category in 1:5) {
  SIM_WORLD3(CAT = category)
}

# for each level create persons totals - load sim_subreg2 and sim_subreg3

setwd("~/AlcPAF_os/simulations")
load("sim_w2_1.RData"); sim_w2_1 <- sim_w2
load("sim_w2_2.RData"); sim_w2_2 <- sim_w2
load("sim_w2_3.RData"); sim_w2_3 <- sim_w2
load("sim_w2_4.RData"); sim_w2_4 <- sim_w2
load("sim_w2_5.RData"); sim_w2_5 <- sim_w2

load("sim_w3_1.RData"); sim_w3_1 <- sim_w3
load("sim_w3_2.RData"); sim_w3_2 <- sim_w3
load("sim_w3_3.RData"); sim_w3_3 <- sim_w3
load("sim_w3_4.RData"); sim_w3_4 <- sim_w3
load("sim_w3_5.RData"); sim_w3_5 <- sim_w3

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
# add level variable back in
pafssima <-
  rbind(sim_w2_1[[1]][,level:=1], sim_w2_1[[2]][,level:=1],
        sim_w2_2[[1]][,level:=2], sim_w2_2[[2]][,level:=2],
        sim_w2_3[[1]][,level:=3], sim_w2_3[[2]][,level:=3],
        sim_w2_4[[1]][,level:=4], sim_w2_4[[2]][,level:=4],
        sim_w2_5[[1]][,level:=5], sim_w2_5[[2]][,level:=5]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_w3_1[[1]][,level:=1], sim_w3_1[[2]][,level:=1],
        sim_w3_2[[1]][,level:=2], sim_w3_2[[2]][,level:=2],
        sim_w3_3[[1]][,level:=3], sim_w3_3[[2]][,level:=3],
        sim_w3_4[[1]][,level:=4], sim_w3_4[[2]][,level:=4],
        sim_w3_5[[1]][,level:=5], sim_w3_5[[2]][,level:=5])  # get 95% CIs of PAF per subreg

pafssim %>%
  group_by(cancer_label, worldgrp, level, sex) %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, worldgrp, age_num, level) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import subreg paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafsplit_world_11.04.21.csv")) # 2020 data
pafs %>%
  rename(worldgrp = world) %>% 
  full_join(pafssimt2 %>%
              filter(!is.na(worldgrp)) %>% 
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI

grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("worldgrp", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(worldgrp, sex, cancer_label, LCIasr,level)
  ) %>% 
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("worldgrp", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(worldgrp, sex, cancer_label, UCIasr,level)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("worldgrp", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(worldgrp, sex, cancer_label, asr,level)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique() -> grpCI

write.csv(grpCI, "totalpafsCIasr_worldsplit_11.04.21.csv", row.names =
            FALSE)

rm(sim_w2_1,sim_w2_2,sim_w2_3,sim_w2_4,sim_w2_5,
   sim_w3_1,sim_w3_2,sim_w3_3,sim_w3_4,sim_w3_5,
   sim_w2, sim_w3)


####---- Simulations and calculating CIs for 10 g level calcs ####


####---- Import country simulations by level

# load simulations data
year <- 2020 # set year
# change working directory
setwd("~/AlcPAF_os/simulations")
# load simulated point estimates introducing age, sex and level variable
load(file = paste0("CANCER_1_1_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <-
  1
PAF_DATA_1_1 <-
  aaf_out_PAF
PAF_DATA_1_1 <-  
  PAF_DATA_1_1[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                   11][str_detect(disease,"_12"), level :=
                                                                        12][str_detect(disease,"_13"), level :=
                                                                             13][str_detect(disease,"_14"), level :=
                                                                                  14][str_detect(disease,"_15"), level :=15]

load(file = paste0("CANCER_1_2_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <-
  1
PAF_DATA_1_2 <-
  aaf_out_PAF
PAF_DATA_1_2 <-
  PAF_DATA_1_2[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]
load(file = paste0("CANCER_1_3_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <-
  1
PAF_DATA_1_3 <-
  aaf_out_PAF
PAF_DATA_1_3 <-
  PAF_DATA_1_3[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]
load(file = paste0("CANCER_1_4_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <-
  1
PAF_DATA_1_4 <-
  aaf_out_PAF
PAF_DATA_1_4 <-
  PAF_DATA_1_4[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]
load(file = paste0("CANCER_1_5_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <-
  1
PAF_DATA_1_5 <-
  aaf_out_PAF
PAF_DATA_1_5 <-
  PAF_DATA_1_5[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]
load(file = paste0("CANCER_1_6_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <-
  1
PAF_DATA_1_6 <-
  aaf_out_PAF
PAF_DATA_1_6 <-
  PAF_DATA_1_6[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]

load(file = paste0("CANCER_2_1_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  1
aaf_out_PAF$sex <-
  2
PAF_DATA_2_1 <-
  aaf_out_PAF
PAF_DATA_2_1 <-
  PAF_DATA_2_1[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]
load(file = paste0("CANCER_2_2_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  2
aaf_out_PAF$sex <-
  2
PAF_DATA_2_2 <-
  aaf_out_PAF
PAF_DATA_2_2 <-
  PAF_DATA_2_2[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]
load(file = paste0("CANCER_2_3_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  3
aaf_out_PAF$sex <-
  2
PAF_DATA_2_3 <-
  aaf_out_PAF
PAF_DATA_2_3 <-
  PAF_DATA_2_3[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]
load(file = paste0("CANCER_2_4_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  4
aaf_out_PAF$sex <-
  2
PAF_DATA_2_4 <-
  aaf_out_PAF
PAF_DATA_2_4 <-
  PAF_DATA_2_4[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]
load(file = paste0("CANCER_2_5_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  5
aaf_out_PAF$sex <-
  2
PAF_DATA_2_5 <-
  aaf_out_PAF
PAF_DATA_2_5 <-
  PAF_DATA_2_5[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]
load(file = paste0("CANCER_2_6_", (year - 10), "split10.RData"))
aaf_out_PAF <-
  as.data.table(aaf_out_PAF)
colnames(aaf_out_PAF)[c(1, 2)] <-
  c("ISO3", "disease")
aaf_out_PAF$paf_age_cancer <-
  6
aaf_out_PAF$sex <-
  2
PAF_DATA_2_6 <-
  aaf_out_PAF
PAF_DATA_2_6 <-
  PAF_DATA_2_6[substr(disease, nchar(disease)-1,nchar(disease))=="_1", level := 
                 1][str_detect(disease,"_2") , level := 
                      2][str_detect(disease,"_3"), level :=
                           3][str_detect(disease,"_4"), level :=
                                4][str_detect(disease,"_5"), level :=
                                     5][str_detect(disease,"_6"), level :=
                                          6][str_detect(disease,"_7"), level :=
                                               7][str_detect(disease,"_8"), level :=
                                                    8][str_detect(disease,"_9"), level :=
                                                         9][str_detect(disease,"_10"), level :=
                                                              10][str_detect(disease,"_11"), level :=
                                                                    11][str_detect(disease,"_12"), level :=
                                                                          12][str_detect(disease,"_13"), level :=
                                                                                13][str_detect(disease,"_14"), level :=
                                                                                      14][str_detect(disease,"_15"), level :=15]




# add sims for former drinkers to sims for current drinkers for each sex/country/cancer/age
sim <-
  list(
    rbind(
      PAF_DATA_1_1,
      PAF_DATA_1_2,
      PAF_DATA_1_3,
      PAF_DATA_1_4,
      PAF_DATA_1_5,
      PAF_DATA_1_6
    ),
    rbind(
      PAF_DATA_2_1,
      PAF_DATA_2_2,
      PAF_DATA_2_3,
      PAF_DATA_2_4,
      PAF_DATA_2_5,
      PAF_DATA_2_6
    )
  )


rm(
  PAF_DATA_1_1,
  PAF_DATA_1_2,
  PAF_DATA_1_3,
  PAF_DATA_1_4,
  PAF_DATA_1_5,
  PAF_DATA_1_6,
  PAF_DATA_2_1,
  PAF_DATA_2_2,
  PAF_DATA_2_3,
  PAF_DATA_2_4,
  PAF_DATA_2_5,
  PAF_DATA_2_6
)


sim_t <- lapply(1:15, function(CAT) {
  lapply(1:2, function(SEX) {
    t <- merge(sim[[SEX]], Gcodes10, by = "disease")
    cat("categ:",toString(CAT),"sex:",toString(SEX),"\n")
    t %>%
      bind_rows(t %>% 
                  filter(str_detect(disease,"Pharynx")) %>% # add rows for hypopharynx
                  mutate(cancer_code=5)) %>% 
      rename_at(vars(3:1002), function(x)
        paste0("sim", x)) %>%
      filter(!str_detect(disease, "_lin")) %>%
      filter(level == CAT) %>% 
      dplyr::select(-level) -> t
  })
})
rm(sim)

# merge simulations and globocan incidence datasets by sex, age, country and cancer type
sim_t. <- lapply(1:15, function(CAT) {
  t3 <- lapply(1:2, function(SEX) {
    t <-
      as.data.table(merge(
        glob_inc[sex == SEX, ],
        sim_t[[CAT]][[SEX]],
        by = c("cancer_code", "paf_age_cancer", "ISO3", "sex"),
        all = T,
        allow.cartesian = TRUE
      ))
    

    t2 <- lapply(unique(t$ISO3), function(REG) {
      cat("categ:",
          toString(CAT),
          "sex:",
          toString(SEX),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      t3 <- t[ISO3 == REG, ]
    })
  })
})

for (CAT in 1:15) {
  dt <- sim_t.[[CAT]]
  filename <- paste0("simsplit10_", CAT, ".RData")
  save(dt, file = filename)
}
rm(sim_t, sim_t.)


# calculate attributable cases per sex, country and simulation number, replace AAF with attributable cases
SIM_COUNTRY <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <- paste0("simsplit10_", CAT, ".RData") # should load as dt
  load(simname)
  
  sim_t2 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t2 <- as.data.table(dt[[SEX]][[REG]]) 
      cat("categ:",
          toString(CAT),
          "sex:",
          toString(SEX),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (i in 12:1011) {
        t2[, as.numeric(i)] <- t2[, ..i] * t2$cases
        
      }
      
      # remove premenop results in older ages
      t2 <-
        t2[(cancer_code == 20 & age_num < 6), disease := "Breast"]
      t2 <-
        t2[(
          cancer_code == 20 &
            age_num > 10 &
            country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)
        ), disease := "Breast_postmen"]
      t2 <-
        t2[(
          cancer_code == 20 &
            age_num %in% c(6:10) &
            country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)
        ), disease := "Breast_premen"]
      t2 <-
        t2[!(cancer_code == 20 &
               str_detect(disease,"Breast_premen") & age_num > 10), ]
      # remove postmenop results in younger ages
      t2 <-
        t2[!(cancer_code == 20 &
               str_detect(disease,"Breast_postmen") & age_num <= 10), ]
      # make male breast cancer 0
      t2[cancer_code == 20 & sex == 1, c(12:1011)] <- 0
      
      t2 %>%
        group_by(age_num) %>%
        mutate_at(vars(12:1011), funs(case_when(
          cancer_code %in% c(39, 40) ~ sum(., na.rm = T),
          TRUE ~ .
        ))) %>%
        bind_rows(
          t2 %>%
            group_by(age_num) %>%
            mutate_at(vars(12:1011), funs(
              case_when(cancer_code %in% c(39, 40) ~ sum(.[Main_analysis == 1], na.rm =
                                                           T))
            )) %>%
            mutate(cancer_label = case_when(
              cancer_code %in% c(39, 40) ~ paste0(cancer_label, "_M")
            )) %>% # add label to main all cancer results
            filter(cancer_code %in% c(39, 40))
        ) %>%
        group_by(age_num) %>%
        mutate_at(vars(12:1011), funs(
          case_when(
            cancer_code == 11 ~ sum(.[cancer_code == 11.1], na.rm = T),
            cancer_code == 6 ~ sum(.[cancer_code ==
                                       6.1], na.rm = T),
            TRUE ~ .
          )
        )) %>% 
        dplyr::select(-Main_analysis, -disease) -> t3
      
      t3 %>% 
        bind_rows(t3 %>% 
                    filter(cancer_code %in% c(3,5)) %>% 
                    group_by(age_num) %>% 
                    mutate(cancer_label = "Pharynx",
                           cancer_code = 41,
                           cases = sum(cases),
                           total = sum(total)) %>%  
                    mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% 
                    unique()) -> t3
      
      t3 <- as.data.table(t3)
      
      t3
    })
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_t210_", CAT, ".RData")
  save(sim_t2, file = filename)
}
for (category in 1:15) {
  SIM_COUNTRY(CAT = category)
}


# calculate AAF per simulation per cancer type for subregion to use for missing countries
SIM_SUBREG <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_t210_", CAT, ".RData") # should load as sim_t2
  load(simname)
  
  sim_subreg <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_t2[[SEX]][[REG]]
      cat("categ:",
          toString(CAT),
          "sex:",
          toString(SEX),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>% 
      left_join(subregmatch) %>%
      group_by(subregion_label_2, cancer_label, age_num) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
      mutate(cases = sum(cases, na.rm = T)) %>%
      ungroup() %>%
      dplyr::select(
        subregion_label_2,
        sex,
        age_num,
        cancer_code,
        cancer_label,
        cases,
        sim1:sim1000
      ) %>%
      unique() %>%
      mutate_at(vars(7:1006), funs(. / cases))  -> t4 # calc aaf per age group
    
    t4 <- as.data.table(t4)
    
    t4
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_subreg10_", CAT, ".RData")
  save(sim_subreg, file = filename)
}
for (category in 1:15) {
  SIM_SUBREG(CAT = category)
}

# replace aafs in missing countries with regional average aafs
SIM_MISSING <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_subreg10_", CAT, ".RData") # should load as sim_subreg
  load(simname)
  simname <-
    paste0("sim_t210_", CAT, ".RData") # should load as sim_t2
  load(simname)
  
  sim_missing <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_t2[[SEX]][[REG]]
      cat("categ:",
          toString(CAT),
          "sex:",
          toString(SEX),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      t2 <- as.data.table(t1)
      
      if (nrow(t2) > 0) {
        if (t2$country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)) {
          t2 %>%
            left_join(subregmatch) %>%
            dplyr::select(-sim1:-sim1000) %>% 
            left_join(
              sim_subreg[[SEX]] %>%
                dplyr::select(
                  subregion_label_2,
                  age_num,
                  cancer_label,
                  sim1:sim1000
                )
            ) %>%
            dplyr::select(-subregion_label_2) %>% 
            mutate_at(vars(11:1010), funs(. * cases))  %>%  # calc aa cases per age
            group_by(age_num) %>%
            mutate_at(vars(11:1010), funs(
              case_when(
                cancer_label %in% c(
                  "All cancers but non-melanoma skin cancer_M",
                  "All cancers_M"
                ) ~ sum(.[cancer_code %in% c(1, 3, 5, 6.1, 8, 9, 11.1, 14, 20)], na.rm = T),
                cancer_label %in% c(
                  "All cancers but non-melanoma skin cancer",
                  "All cancers"
                ) ~ sum(.[cancer_code %in% c(1, 3, 5, 6.1, 8, 9, 11.1, 14, 20, 7, 13)], na.rm =
                          T),
                cancer_label == "Liver and intrahepatic bile ducts" ~ sum(.[cancer_code ==
                                                                              11.1], na.rm = T),
                cancer_label == "Oesophagus" ~ sum(.[cancer_code ==
                                                       6.1], na.rm = T),
                cancer_label == "Pharynx" ~ sum(.[cancer_code %in% c(3,5)], na.rm = T),
                TRUE ~ .
              )
            )) -> t3
          
          t3
          
        } else {
          t2
        }
      } else {
        t2
      }
      
    })
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_missing10_", CAT, ".RData")
  save(sim_missing, file = filename)
}
for (category in 1:15) {
  SIM_MISSING(CAT = category)
}



SIM_WORLD2 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing10_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for world by age
  sim_w2 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      mutate(worldgrp = "world") %>%
      group_by(worldgrp, cancer_label, age_num) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each age group & subreg
      mutate(cases = sum(cases, na.rm = T)) %>%
      ungroup() %>%
      dplyr::select(
        worldgrp,
        sex,
        age_num,
        cancer_code,
        cancer_label,
        cases,
        sim1:sim1000
      ) %>%
      unique() -> t4
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$worldgrp)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        for (AGE in unique(t4$age_num)) {
          # calc 95%CIs for aa cases by age per subreg
          
          
          t4[(worldgrp == REG &
                cancer_label == CAN &
                age_num == AGE), LCIage := quantile(t4[(worldgrp == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.025, na.rm = TRUE)]
          t4[(worldgrp == REG &
                cancer_label == CAN &
                age_num == AGE), UCIage := quantile(t4[(worldgrp == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.975, na.rm = TRUE)]
          t4[(worldgrp == REG &
                cancer_label == CAN &
                age_num == AGE), mdnage := quantile(t4[(worldgrp == REG &
                                                          cancer_label == CAN &
                                                          age_num == AGE), 7:1006], probs = 0.5, na.rm = TRUE)]
          
        }
      }
    }
    
    t5 <-
      t4[, c(
        "worldgrp",
        "age_num",
        "cancer_code",
        "sex",
        "cancer_label",
        "LCIage",
        "UCIage",
        "mdnage",
        "cases"
      )]
    t5
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_w210_", CAT, ".RData")
  save(sim_w2, file = filename)
}
for (category in 1:15) {
  SIM_WORLD2(CAT = category)
}

# calculate total aa cases and PAF per cancer type by world
# calculate AAF per simulation per cancer type for world to use for missing countries
SIM_WORLD3 <- function(CAT) {
  setwd("~/AlcPAF_os/simulations")
  simname <-
    paste0("sim_missing10_", CAT, ".RData") # should load as sim_missing
  load(simname)
  
  # calculate aa cases per cancer type for world by age
  sim_w3 <- lapply(1:2, function(SEX) {
    t <- lapply(1:199, function(REG) {
      t1 <- sim_missing[[SEX]][[REG]]
      
      t2 <- as.data.table(t1)
      t2
    })
    
    t3 <- do.call(rbind.data.frame, t)
    
    t3 %>%
      mutate(worldgrp = "world") %>%
      group_by(worldgrp, cancer_label) %>%
      mutate_at(vars(11:1010), funs(sum(., na.rm = T))) %>% # sum aa cases for each subreg
      group_by(worldgrp, cancer_label, age_num) %>%
      mutate(total = sum(total)) %>%
      ungroup() %>%
      dplyr::select(worldgrp,
                    sex,
                    cancer_code,
                    cancer_label,
                    total,
                    sim1:sim1000) %>%
      unique() %>%
      mutate_at(vars(6:1005), funs(. / total)) -> t4 # calc PAF per subreg
    
    t4 <- as.data.table(t4)
    
    for (REG in unique(t4$worldgrp)) {
      cat("sex:",
          toString(SEX),
          "categ:",
          toString(CAT),
          ", reg:",
          toString(REG),
          "\n") # print sex and country number
      
      for (CAN in unique(t4$cancer_label)) {
        t4[(worldgrp == REG &
              cancer_label == CAN), totLCI := quantile(t4[(worldgrp == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.025, na.rm = TRUE)]
        t4[(worldgrp == REG &
              cancer_label == CAN), totUCI := quantile(t4[(worldgrp == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.975, na.rm = TRUE)]
        t4[(worldgrp == REG &
              cancer_label == CAN), totmdn := quantile(t4[(worldgrp == REG &
                                                             cancer_label == CAN), 6:1005], probs = 0.5, na.rm = TRUE)]
        
        
      }
    }
    
    t5 <-
      t4[, c(
        "worldgrp",
        "sex",
        "cancer_label",
        "total",
        "totLCI",
        "totUCI",
        "totmdn"
      )]
    t5
    
  })
  setwd("~/AlcPAF_os/simulations")
  filename <- paste0("sim_w310_", CAT, ".RData")
  save(sim_w3, file = filename)
}
for (category in 1:15) {
  SIM_WORLD3(CAT = category)
}

# for each level create persons totals - load sim_subreg2 and sim_subreg3

setwd("~/AlcPAF_os/simulations") #update with 10gs
load("sim_w210_1.RData"); sim_w2_1 <- sim_w2
load("sim_w210_2.RData"); sim_w2_2 <- sim_w2
load("sim_w210_3.RData"); sim_w2_3 <- sim_w2
load("sim_w210_4.RData"); sim_w2_4 <- sim_w2
load("sim_w210_5.RData"); sim_w2_5 <- sim_w2
load("sim_w210_6.RData"); sim_w2_6 <- sim_w2
load("sim_w210_7.RData"); sim_w2_7 <- sim_w2
load("sim_w210_8.RData"); sim_w2_8 <- sim_w2
load("sim_w210_9.RData"); sim_w2_9 <- sim_w2
load("sim_w210_10.RData"); sim_w2_10 <- sim_w2
load("sim_w210_11.RData"); sim_w2_11 <- sim_w2
load("sim_w210_12.RData"); sim_w2_12 <- sim_w2
load("sim_w210_13.RData"); sim_w2_13 <- sim_w2
load("sim_w210_14.RData"); sim_w2_14 <- sim_w2
load("sim_w210_15.RData"); sim_w2_15 <- sim_w2

load("sim_w310_1.RData"); sim_w3_1 <- sim_w3
load("sim_w310_2.RData"); sim_w3_2 <- sim_w3
load("sim_w310_3.RData"); sim_w3_3 <- sim_w3
load("sim_w310_4.RData"); sim_w3_4 <- sim_w3
load("sim_w310_5.RData"); sim_w3_5 <- sim_w3
load("sim_w310_6.RData"); sim_w3_6 <- sim_w3
load("sim_w310_7.RData"); sim_w3_7 <- sim_w3
load("sim_w310_8.RData"); sim_w3_8 <- sim_w3
load("sim_w310_9.RData"); sim_w3_9 <- sim_w3
load("sim_w310_10.RData"); sim_w3_10 <- sim_w3
load("sim_w310_11.RData"); sim_w3_11 <- sim_w3
load("sim_w310_12.RData"); sim_w3_12 <- sim_w3
load("sim_w310_13.RData"); sim_w3_13 <- sim_w3
load("sim_w310_14.RData"); sim_w3_14 <- sim_w3
load("sim_w310_15.RData"); sim_w3_15 <- sim_w3

# persons totals
# create attributable cases of LCI, UCI and median
# bind lists or create new list of sum of m and f lists
# add level variable back in
pafssima <-
  rbind(sim_w2_1[[1]][,level:=1], sim_w2_1[[2]][,level:=1],
        sim_w2_2[[1]][,level:=2], sim_w2_2[[2]][,level:=2],
        sim_w2_3[[1]][,level:=3], sim_w2_3[[2]][,level:=3],
        sim_w2_4[[1]][,level:=4], sim_w2_4[[2]][,level:=4],
        sim_w2_5[[1]][,level:=5], sim_w2_5[[2]][,level:=5],
        sim_w2_6[[1]][,level:=6], sim_w2_6[[2]][,level:=6],
        sim_w2_7[[1]][,level:=7], sim_w2_7[[2]][,level:=7],
        sim_w2_8[[1]][,level:=8], sim_w2_8[[2]][,level:=8],
        sim_w2_9[[1]][,level:=9], sim_w2_9[[2]][,level:=9],
        sim_w2_10[[1]][,level:=10], sim_w2_10[[2]][,level:=10],
        sim_w2_11[[1]][,level:=11], sim_w2_11[[2]][,level:=11],
        sim_w2_12[[1]][,level:=12], sim_w2_12[[2]][,level:=12],
        sim_w2_13[[1]][,level:=13], sim_w2_13[[2]][,level:=13],
        sim_w2_14[[1]][,level:=14], sim_w2_14[[2]][,level:=14],
        sim_w2_15[[1]][,level:=15], sim_w2_15[[2]][,level:=15]) # get 95% CIs of cases per age per subreg
pafssim <-
  rbind(sim_w3_1[[1]][,level:=1], sim_w3_1[[2]][,level:=1],
        sim_w3_2[[1]][,level:=2], sim_w3_2[[2]][,level:=2],
        sim_w3_3[[1]][,level:=3], sim_w3_3[[2]][,level:=3],
        sim_w3_4[[1]][,level:=4], sim_w3_4[[2]][,level:=4],
        sim_w3_5[[1]][,level:=5], sim_w3_5[[2]][,level:=5],
        sim_w3_6[[1]][,level:=6], sim_w3_6[[2]][,level:=6],
        sim_w3_7[[1]][,level:=7], sim_w3_7[[2]][,level:=7],
        sim_w3_8[[1]][,level:=8], sim_w3_8[[2]][,level:=8],
        sim_w3_9[[1]][,level:=9], sim_w3_9[[2]][,level:=9],
        sim_w3_10[[1]][,level:=10], sim_w3_10[[2]][,level:=10],
        sim_w3_11[[1]][,level:=11], sim_w3_11[[2]][,level:=11],
        sim_w3_12[[1]][,level:=12], sim_w3_12[[2]][,level:=12],
        sim_w3_13[[1]][,level:=13], sim_w3_13[[2]][,level:=13],
        sim_w3_14[[1]][,level:=14], sim_w3_14[[2]][,level:=14],
        sim_w3_15[[1]][,level:=15], sim_w3_15[[2]][,level:=15])  # get 95% CIs of PAF per subreg

pafssim %>%
  dplyr::filter(!is.na(total)) %>% 
  group_by(cancer_label, worldgrp, level, sex) %>%
  mutate(
    totLCIc = totLCI * total,
    totUCIc = totUCI * total,
    totmdnc = totmdn * total
  ) %>%
  left_join(pafssima %>% dplyr::filter(!is.na(LCIage)) %>% dplyr::select(-cases)) -> pafssimt

# create persons totals and merge with regions data
pafssimt %>%
  bind_rows(
    pafssimt %>%
      group_by(cancer_label, worldgrp, age_num, level) %>%
      mutate(
        totLCIc = sum(totLCIc, na.rm = T),
        totUCIc = sum(totUCIc, na.rm = T),
        totmdnc = sum(totmdnc, na.rm = T),
        LCIage = sum(LCIage, na.rm = T),
        UCIage = sum(UCIage, na.rm = T),
        mdnage = sum(mdnage, na.rm = T),
        total = sum(total, na.rm = T),
        totLCI = totLCIc / total,
        totUCI = totUCIc / total,
        totmdn = totmdnc / total,
        sex = 0
      ) %>%
      unique()
  ) ->  pafssimt2

# import paf estimate
setwd("~/AlcPAF_os/alc_results")
pafs 	<-
  as.data.table(read.csv("pafsplit10_world_11.04.21.csv")) # 2020 data
pafs %>%
  filter(!is.na(age_num)) %>% 
  rename(worldgrp = world) %>% 
  full_join(pafssimt2 %>%
              filter(!is.na(worldgrp)) %>% 
              dplyr::select(-cancer_code,-total)) -> grpCI
grpCI %>%
  mutate(
    UCIage = case_when(is.na(UCIage) ~ 0,
                       TRUE ~ UCIage),
    LCIage = case_when(is.na(LCIage) ~ 0,
                       TRUE ~ LCIage)
  ) -> grpCI

grpCI %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "LCIage",
        var_py = "grppyage",
        group_by = c("worldgrp", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "LCIasr"
      ) %>%
      dplyr::select(worldgrp, sex, cancer_label, LCIasr,level)
  ) %>% 
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "UCIage",
        var_py = "grppyage",
        group_by = c("worldgrp", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "UCIasr"
      ) %>%
      dplyr::select(worldgrp, sex, cancer_label, UCIasr,level)
  ) %>%
  left_join(
    grpCI %>%
      HR_csu_asr(
        var_age = "age_num",
        var_cases = "grpaacasesage",
        var_py = "grppyage",
        group_by = c("worldgrp", "sex", "cancer_label","level"),
        var_age_group = "cancer_label",
        var_asr = "asr"
      ) %>%
      dplyr::select(worldgrp, sex, cancer_label, asr,level)
  ) %>%
  dplyr::select(-age_num,-grpaacasesage,-grppyage,-LCIage,-UCIage,-mdnage) %>%
  unique()  -> grpCI

write.csv(grpCI, "totalpafsCIasr_worldsplit10_26.04.21.csv", row.names =
            FALSE)

rm(sim_w2_1,sim_w2_2,sim_w2_3,sim_w2_4,sim_w2_5,
   sim_w2_6,sim_w2_7,sim_w2_8,sim_w2_9,sim_w2_10,
   sim_w2_11,sim_w2_12,sim_w2_13,sim_w2_14,sim_w2_15,
   sim_w3_1,sim_w3_2,sim_w3_3,sim_w3_4,sim_w3_5,
   sim_w3_6,sim_w3_7,sim_w3_8,sim_w3_9,sim_w3_10,
   sim_w3_11,sim_w3_12,sim_w3_13,sim_w3_14,sim_w3_15,
   sim_w2, sim_w3)

