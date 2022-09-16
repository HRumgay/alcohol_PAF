######################################################################
### This script calculates PAFs of cases using the point estimates ###
### produced in the script 1.paf_calc.R and Globocan 2018 data     ###
######################################################################


# Run the following libraries
library(tidyverse)
library(data.table)
library(Rcan)

# Clear everything in the global environment
rm(list = ls())


####--- Import data ---####
# import Globocan data
setwd("//Inti/cin/DataShare/Globocan2020")
DATA_GLOB	<- as.data.table(read.csv("Globocan.csv"))
# import Globocan data
setwd("//Inti/cin/DataShare/Globocan2018")
glob_inc18	<- as.data.table(read.csv("Globocan.csv"))
# Globocan incidence = type 0, mort = type 1
# Globocan persons = sex 0, men = sex 1, women = sex 2

# import HCC and oesophageal SCC cases
HCCcases <- as.data.table(read.csv("//Inti/cin/Users/RumgayH/RProjects/Liver_subtypes/results/cases_method1_18.05.20.csv"))
OSCCcases <- as.data.table(read.csv("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_data/OSCC_MA.csv"))

# create list of regions and countries to match from HCC data
HCCcases %>% 
  dplyr::select(REGION_LABEL:COUNTRY_LABEL, COUNTRY_CODE, UNTERM_COUNTRY_LABEL:INCOME_GROUP) %>% 
  unique() -> countrymatch

# import AAFs results
setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_results")
AAF_CANCER 	<- as.data.table(read.csv("AAF_CANCER_01.02.21.csv"))
AAF_CANCER_split 	<- as.data.table(read.csv("AAF_CANCER_split_01.02.21.csv"))
AAF_CANCER_split10 	<- as.data.table(read.csv("AAF_CANCER_split10_07.04.21.csv"))

# load age structure csv
setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_data")
AGE_STRUCT <- as.data.table(read.csv("PAF_INC_STRUCTURE_AGE_HR.csv"))

# import populations data from Globocan for calculating rates
#DATA_GLOB %>% 
#  select(country_code,sex,age,py) %>% unique()->pop_glob
#write.csv(pop_glob,"pops_globocan.csv", row.names=FALSE)
pop_glob <- read.csv("pops_globocan.csv")
# rename age variable to same as pafs dataset
names(pop_glob)[names(pop_glob) == "age"] <- "age_num"

####--- Reshape Globocan data ---####
# filter Globocan incidence data
glob_inc <- DATA_GLOB[(type == 0), , ][order(country_code, sex, age), ][cancer_code %in% c(1,3,5,6,7,8,9,11,13,14,20,39,40), ] 
glob_inc <- glob_inc[country_code<900,][sex != 0,] # filter to countries only and remove persons cases from Globocan data
glob_inc %>% dplyr::select(-country_label,-type,-py) -> glob_inc

# filter Globocan mortality data
glob_mort <- DATA_GLOB[(type == 1), , ][order(country_code, sex, age), ][cancer_code %in% c(1,3,5,6,7,8,9,11,13,14,20,39,40), ] 
glob_mort <- glob_mort[country_code<900,][sex != 0,] # filter to countries only and remove persons cases from Globocan data
glob_mort %>% dplyr::select(-country_label,-type,-py) -> glob_mort

# filter Globocan 2018 data for OSCC and HCC
glob_inc18 <- glob_inc18[(type == 0), , ][order(country_code, sex, age), ][cancer_code %in% c(6,11), ] 
glob_inc18 <- glob_inc18[country_code<900,][sex != 0,] # filter to countries only and remove persons cases from Globocan data
glob_inc18 %>% dplyr::select(-country_label,-type,-py) -> glob_inc18

# reshape HCC cases to same format as Globocan data
HCCcases %>% 
  filter(count.name == "HCC") %>% 
  mutate(count.name = "Liver_HCC") %>% 
  dplyr::select(COUNTRY_CODE,SEX_LABEL,AGEG,count.name,count.value) %>% 
  rename(country_code = COUNTRY_CODE,
         cancer_label = count.name,
         cases = count.value) %>% 
  mutate(cancer_code = 11.1, # to match RRs list
         sex = case_when(SEX_LABEL == "Men" ~ 1,
                         SEX_LABEL == "Women" ~ 2),
         age = case_when(AGEG == "0-4" ~ 1,
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
                         AGEG == "85+" ~ 18)) %>% 
  dplyr::select(-SEX_LABEL, -AGEG) %>% 
  group_by(sex, country_code) %>%
  mutate(total = sum(cases)) %>% 
  bind_rows(OSCCcases %>% 
              filter(un_code != 1000) %>% 
              rename(country_code = un_code,
                     cases = N) %>% 
              mutate(cancer_code = 6.1, # to match RRs list
                     sex = case_when(sex == "Males" ~ 1,
                                     sex == "Females" ~ 2)) %>% 
              group_by(sex, country_code) %>%
              mutate(total = sum(cases),
                     cancer_label = paste0(cancer_label,"_SCC")) %>% 
              dplyr::select(country_code, cancer_label, cancer_code, sex, age, cases, total)) -> HCCcases

# add HCC cases to globocan dataset - incidence
glob_inc18 <- bind_rows(glob_inc18, HCCcases)

# need to calculate distribution of subtype of liver cancer and oesophageal cancer to apply to mortality
# proportion of SCC <64, 65+ - group cases <64 (age 1-13) and 65+ (age 14-18)
# proportion of HCC - group all ages
glob_inc18 %>% 
  group_by(country_code, age, sex) %>% 
  mutate(prop = case_when(cancer_code==11.1 ~ cases[cancer_code==11],
                          cancer_code==6.1 ~ cases[cancer_code==6],
                          TRUE ~ cases)) %>%
  filter(cancer_code %in% c(11.1, 6.1)) %>% 
  mutate(grp = case_when(cancer_code==6.1 & age > 13 ~ 2, #split age groups for oesoph
                         TRUE ~ 1)) %>% 
  group_by(country_code, grp, sex, cancer_code) %>% 
  mutate(prop = sum(cases)/sum(prop)) %>% 
  ungroup() %>% 
  dplyr::select(-grp) -> glob_inc18

# apply 2018 proportions of HCC and OSCC to 2020 data
glob_inc18 %>% 
  bind_rows(glob_inc %>% 
              mutate(cancer_code = as.numeric(cancer_code))) %>% 
  group_by(country_code, age, sex) %>%
  mutate(cases = case_when(cancer_code==11.1 ~ cases[cancer_code==11]*prop,
                           cancer_code==6.1 ~ cases[cancer_code==6]*prop,
                           TRUE ~ cases)) %>%
  group_by(country_code, sex) %>%
  mutate(total = case_when(cancer_code==11.1 ~ sum(cases[cancer_code==11.1]),
                           cancer_code==6.1 ~ sum(cases[cancer_code==6.1]),
                           TRUE ~ total),
         cases = case_when(is.na(cases) ~ 0, TRUE ~ cases),
         total = case_when(is.na(total) ~ 0, TRUE ~ total)) %>%
  dplyr::select(-prop) -> glob_inc
glob_inc18 %>% 
  bind_rows(glob_mort %>% 
              mutate(cancer_code = as.numeric(cancer_code))) %>% 
  group_by(country_code, age, sex) %>%
  mutate(cases = case_when(cancer_code==11.1 ~ cases[cancer_code==11]*prop,
                           cancer_code==6.1 ~ cases[cancer_code==6]*prop,
                           TRUE ~ cases)) %>%
  group_by(country_code, sex) %>%
  mutate(total = case_when(cancer_code==11.1 ~ sum(cases[cancer_code==11.1]),
                           cancer_code==6.1 ~ sum(cases[cancer_code==6.1]),
                           TRUE ~ total),
         cases = case_when(is.na(cases) ~ 0, TRUE ~ cases),
         total = case_when(is.na(total) ~ 0, TRUE ~ total)) %>%
  dplyr::select(-prop) -> glob_mort
rm(glob_inc18)

names(glob_inc)[names(glob_inc) == "age"] <- "age_num" # rename age variable
names(glob_mort)[names(glob_mort) == "age"] <- "age_num" # rename age variable


####--- Match alcohol data and Globocan countries ---####

setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_data")
regionsmatch 	<- as.data.table(read.csv("countries_region_matching.csv"))
regionsmatch %>% 
  dplyr::select(M49Code, ISO3) %>% 
  unique() %>% 
  rename(country_code = M49Code) -> regionsmatch
# merge regions match with globocan data
glob_inc <- merge(glob_inc, regionsmatch, by = "country_code")
glob_mort <- merge(glob_mort, regionsmatch, by = "country_code")

####--- Match cancer names ---####

setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_data")
RRs	<- as.data.table(read.csv("WCRF_RRs3.csv"))
RRsplit	<- as.data.table(read.csv("WCRF_RRs_split4.csv"))
RRsplit10	<- as.data.table(read.csv("WCRF_RRs_split10.csv"))
# filter to non-regional RRs
RRs <- RRs[Regional==0]
RRsplit <- RRsplit[Regional==0]
RRsplit10 <- RRsplit10[Regional==0]
# separate disease name and Globocan number
Gcode <- RRs[sex==1, c("disease", "cancer_code", "Main_analysis")]
Gcodes <- RRsplit[sex==1, c("disease", "cancer_code", "Main_analysis")] 
Gcodes10 <- RRsplit10[sex==1, c("disease", "cancer_code", "Main_analysis")] 
# add nonlin names and Globocan codes
Gcode <- rbind(Gcode, data.frame(c("Oesophagus_SCC_nonlin","Colon_lin","Rectum_lin"), c(6.1,8,9), c(1,1,1)), use.names=FALSE)
Gcodes <- rbind(Gcodes, data.frame(c("Oesophagus_SCC_nonlin_low", "Oesophagus_SCC_nonlin_light", "Oesophagus_SCC_nonlin_mod", "Oesophagus_SCC_nonlin_heavy"), 
                                   c(6.1,6.1,6.1,6.1), 
                                   c(1,1,1,1)), use.names=FALSE)
Gcodes10 <- rbind(Gcodes10, data.frame(c("Oesophagus_SCC_nonlin_1", "Oesophagus_SCC_nonlin_2", "Oesophagus_SCC_nonlin_3", "Oesophagus_SCC_nonlin_4", "Oesophagus_SCC_nonlin_5",
                                         "Oesophagus_SCC_nonlin_6", "Oesophagus_SCC_nonlin_7", "Oesophagus_SCC_nonlin_8", "Oesophagus_SCC_nonlin_9", "Oesophagus_SCC_nonlin_10",
                                         "Oesophagus_SCC_nonlin_11", "Oesophagus_SCC_nonlin_12", "Oesophagus_SCC_nonlin_13", "Oesophagus_SCC_nonlin_14", "Oesophagus_SCC_nonlin_15"), 
                                       c(6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1,6.1), 
                                       c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)), use.names=FALSE)

# change variable names to lower case to be able to match
names(AAF_CANCER) <- tolower(names(AAF_CANCER))
names(AAF_CANCER_split) <- tolower(names(AAF_CANCER_split))
names(AAF_CANCER_split10) <- tolower(names(AAF_CANCER_split10))
# match cancer name and Globocan code with PAFs data
AAF_CANCER <- merge(AAF_CANCER, Gcode, by = "disease") # assign Globocan cancer code to AAFs
AAF_CANCER_split <- merge(AAF_CANCER_split, Gcodes, by = "disease") # assign Globocan cancer code to AAFs
AAF_CANCER_split10 <- merge(AAF_CANCER_split10, Gcodes10, by = "disease") # assign Globocan cancer code to AAFs

####--- Merge alcohol and cancer data ---####

glob_inc <- merge(AGE_STRUCT, glob_inc, by = "age_num") # merge PAF age variable into inc data
glob_mort <- merge(AGE_STRUCT, glob_mort, by = "age_num") # merge PAF age variable into inc data
names(AAF_CANCER)[names(AAF_CANCER) == "region"] <- "ISO3" # rename country var
names(AAF_CANCER)[names(AAF_CANCER) == "age_category"] <- "paf_age_cancer"
names(AAF_CANCER_split)[names(AAF_CANCER_split) == "region"] <- "ISO3" # rename country var
names(AAF_CANCER_split)[names(AAF_CANCER_split) == "age_category"] <- "paf_age_cancer"
names(AAF_CANCER_split10)[names(AAF_CANCER_split10) == "region"] <- "ISO3" # rename country var
names(AAF_CANCER_split10)[names(AAF_CANCER_split10) == "age_category"] <- "paf_age_cancer"

# create PAFs dataset
aafs <- AAF_CANCER[, c("ISO3", "sex", "paf_age_cancer", "disease","cancer_code", "aaf", "Main_analysis")]
# create PAFs dataset for different levels
aafsplit <- AAF_CANCER_split[, c("ISO3", "sex", "paf_age_cancer", "disease","cancer_code", "aaf", "Main_analysis")]
aafsplit10 <- AAF_CANCER_split10[, c("ISO3", "sex", "paf_age_cancer", "disease","cancer_code", "aaf", "Main_analysis")]

# add grouping for hypopharynx
aafs %>% 
  bind_rows(aafs %>% 
              filter(str_detect(disease,"Pharynx")) %>% 
              mutate(cancer_code = 5))  -> aafs
aafsplit %>% 
  bind_rows(aafsplit %>% 
              filter(str_detect(disease,"Pharynx")) %>% 
              mutate(cancer_code = 5))  -> aafsplit
aafsplit10 %>% 
  bind_rows(aafsplit10 %>% 
              filter(str_detect(disease,"Pharynx")) %>% 
              mutate(cancer_code = 5))  -> aafsplit10

# calculate attributable cases and PAF for each age group, sex and cancer per country
# for incidence and mortality
for (dat in c("glob_inc", "glob_mort")) {
  
  glob <- get(dat)
  
  # merge AAFs data with incidence data
  pafs <- merge(glob, aafs, by = c("ISO3", "sex", "paf_age_cancer", "cancer_code"), all = TRUE)
  # merge AAFs data with incidence data
  pafsplit <- merge(glob, aafsplit, by = c("ISO3", "sex", "paf_age_cancer", "cancer_code"),allow.cartesian=TRUE, all = TRUE)
  
  
  ####--- Create age-specific alc attributable cases 
  
  # calculate aa cases per age group per cancer type
  pafs <- pafs[, aacases:= aaf*cases][age_num<6, aacases := 0][age_num<6, aaf := 0]
  pafsplit <- pafsplit[, aacases:= aaf*cases][age_num<6, aacases := 0][age_num<6, aaf := 0]
  
  # remove premenop results in older ages
  pafs <- pafs[(cancer_code==20 & age_num <6),disease:= "Breast"]
  pafs <- pafs[(cancer_code==20 & age_num >10 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)),disease:= "Breast_postmen"]
  pafs <- pafs[(cancer_code==20 & age_num %in% c(6:10) & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)),disease:= "Breast_premen"]
  pafsplit <- pafsplit[(cancer_code==20 & age_num <6),disease:= "Breast"]
  pafsplit <- pafsplit[(cancer_code==20 & age_num >10 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)),disease:= "Breast_postmen"]
  pafsplit <- pafsplit[(cancer_code==20 & age_num %in% c(6:10) & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)),disease:= "Breast_premen"]
  pafs <- pafs[!(cancer_code==20 & disease == "Breast_premen" & age_num > 10),]
  pafsplit <- pafsplit[!(cancer_code==20 & str_detect(disease,"Breast_premen") & age_num > 10),] 
  # remove postmenop results in younger ages
  pafs <- pafs[!(cancer_code==20 & disease == "Breast_postmen" & age_num <= 10),]
  pafsplit <- pafsplit[!(cancer_code==20 & str_detect(disease,"Breast_postmen") & age_num <= 10),] 
  # remove NAs for male breast cancer, replace with 0
  pafs <- pafs[cancer_code==20 & sex==1, aacases:=0]
  pafsplit <- pafsplit[cancer_code==20 & sex==1, aacases:=0]
  
  
  # create different pafs datasets for Bagnardi CRC and WCRF CRC estimates
  pafs %>% 
    filter(is.na(disease) | !str_detect(disease, "_lin")) -> pafs2     # keep WCRF CRC estimates
  pafsplit %>% 
    filter(is.na(disease) | !str_detect(disease, "_lin")) -> pafsplit2 # keep WCRF CRC estimates
  
#  pafs %>% 
#    filter(is.na(disease) | !str_detect(disease, "above20")) -> pafs  # keep Bagnardi CRC estimates
#  pafsplit %>% 
#    filter(is.na(disease) | !str_detect(disease, "above20")) -> pafsplit # keep Bagnardi CRC estimates

  
    
  #--- following was used to calc aa cases for Bagnardi CRC estimates - can probably remove this section if sticking to WCRF
  
  # bind second set of all cancer total excluding stomach and pancreas for main analysis results
  # sum total aa cases for all cancers and assign HCC aa cases to total liver aa cases
#  pafs %>% 
#    group_by(country_code, age_num, sex) %>% 
#    mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases, na.rm=T), 
#                               TRUE ~ aacases)) %>%
#    bind_rows(pafs %>% 
#                group_by(country_code, age_num, sex) %>% 
#                mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases[Main_analysis==1], na.rm=T)),
#                       cancer_label = case_when(cancer_code %in% c(39,40) ~ paste0(cancer_label,"_M"))) %>% # add label to main all cancer results
#                filter(cancer_code %in% c(39,40))) %>% 
#    mutate(aacases = case_when(cancer_code == 11 ~ sum(aacases[cancer_code==11.1], na.rm=T),
#                               cancer_code == 6 ~ sum(aacases[cancer_code==6.1], na.rm=T),
#                               TRUE ~ aacases),
#           aaf = case_when(is.na(aaf) ~ aacases/cases, # fill in gaps for AAF
#                           TRUE ~ aaf)) -> pafsmf
#  pafsplit %>% 
#    mutate(level = case_when(str_detect(disease,"light") ~ 1,
#                             str_detect(disease,"mod") ~ 2,
#                             str_detect(disease,"heavy") ~ 3,
#                             str_detect(disease,"former") ~ 4,
#                             TRUE ~ 1)) %>% 
#    bind_rows(pafsplit %>% filter(is.na(disease)) %>% 
#                mutate(level = 2)) %>% 
#    bind_rows(pafsplit %>% filter(is.na(disease)) %>% 
#                mutate(level = 3)) %>% 
#    bind_rows(pafsplit %>% filter(is.na(disease)) %>% 
#                mutate(level = 4)) %>% 
#    group_by(country_code, age_num, sex, level) %>% 
#    mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases, na.rm=T), 
#                               TRUE ~ aacases)) -> pafsplit
#  
#  pafsplit %>% 
#    bind_rows(pafsplit %>% 
#                group_by(country_code, age_num, sex, level) %>% 
#                mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases[Main_analysis==1], na.rm=T)),
#                       cancer_label = case_when(cancer_code %in% c(39,40) ~ paste0(cancer_label,"_M"))) %>% # add label to main all cancer results
#                filter(cancer_code %in% c(39,40))) %>% 
#    group_by(country_code, age_num, sex, level) %>%
#    mutate(aacases = case_when(cancer_code == 11 ~ sum(aacases[cancer_code==11.1], na.rm=T),
#                               cancer_code == 6 ~ sum(aacases[cancer_code==6.1], na.rm=T),
#                               TRUE ~ aacases),
#           aaf = case_when(is.na(aaf) ~ aacases/cases, # fill in gaps for AAF
#                           TRUE ~ aaf)) -> pafsplit
#  
#  # create persons totals
#  pafsmf %>% 
#    group_by(country_code,age_num, cancer_label) %>% 
#    mutate(sex = 0,
#           cases = sum(cases, na.rm=T),
#           total = sum(total, na.rm=T),
#           aacases = sum(aacases, na.rm=T),
#           aaf = aacases/cases) %>% unique() -> pafsp
#  # bind mf and persons
#  pafsmfp <- rbind(pafsmf, pafsp)
#  
#  # create persons totals
#  pafsplit %>% 
#    group_by(country_code, ISO3, age_num, cancer_label, level) %>% 
#    mutate(sex = 0,
#           cases = sum(cases, na.rm=T),
#           total = sum(total, na.rm=T),
#           aacases = sum(aacases, na.rm=T),
#           aaf = aacases/cases) %>% unique() -> pafsplitp
#  # bind mf and persons
#  pafsplit <- rbind(pafsplit, pafsplitp)
#  
  
  
  ####--- Calculate PAFs for Bagnardi CRC estimates
  # create pafs per sex, cancer type, and country
#  pafsmfp %>%
#    group_by(country_code, sex, cancer_label) %>% 
#    mutate(totaacases = sum(aacases),
#           totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
#                              TRUE ~ 0)) -> pafsmfp #check if India has duplicated totaacases for males _M?
#  # create pafs per sex, cancer type, and country
#  pafsplit %>%
#    group_by(country_code, sex, cancer_label, level) %>% 
#    mutate(totaacases = sum(aacases),
#           totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
#                              TRUE ~ 0)) -> pafsplit
#  rm(pafsmf, pafsp, pafsplitp)
#  
#  names(countrymatch) <- tolower(names(countrymatch))
#  pafsmfp <- merge(countrymatch, pafsmfp, by = "country_code")
#  pafsplit <- merge(countrymatch, pafsplit, by = "country_code")
#  
#  
#  pafsmfp %>% 
#    bind_rows(pafsmfp %>% 
#                filter(cancer_label %in% c("Colon", "Rectum")) %>% 
#                group_by(country_label, sex, age_num) %>% 
#                mutate(cancer_label = "Colorectum",
#                       disease= "Colorectum",
#                       cancer_code = NA,
#                       cases = sum(cases),
#                       total = sum(total),
#                       aacases = sum(aacases),
#                       aaf = aacases/cases,
#                       totaacases = sum(totaacases),
#                       totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
#                                          TRUE ~ 0)) %>% 
#                unique()) %>% 
#    bind_rows(pafsmfp %>% 
#                filter(cancer_label %in% c("Lip, oral cavity", "Oropharynx", "Larynx")) %>% 
#                group_by(country_label, sex, age_num) %>% 
#                mutate(cancer_label = "Head & Neck",
#                       disease = "Head & Neck",
#                       cancer_code = NA,
#                       cases = sum(cases),
#                       total = sum(total),
#                       aacases = sum(aacases),
#                       aaf = aacases/cases,
#                       totaacases = sum(totaacases),
#                       totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
#                                          TRUE ~ 0)) %>% 
#                unique()) -> pafsmfp
#  
#  pafsplit %>% 
#    bind_rows(pafsplit %>% 
#                filter(cancer_label %in% c("Colon", "Rectum")) %>% 
#                group_by(country_label, sex, age_num, level) %>% 
#                mutate(cancer_label = "Colorectum",
#                       disease= "Colorectum",
#                       cancer_code = NA,
#                       cases = sum(cases),
#                       total = sum(total),
#                       aacases = sum(aacases),
#                       aaf = aacases/cases,
#                       totaacases = sum(totaacases),
#                       totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
#                                          TRUE ~ 0)) %>% 
#                unique()) %>% 
#    bind_rows(pafsplit %>% 
#                filter(cancer_label %in% c("Lip, oral cavity", "Oropharynx", "Larynx")) %>% 
#                group_by(country_label, sex, age_num, level) %>% 
#                mutate(cancer_label = "Head & Neck",
#                       disease = "Head & Neck",
#                       cancer_code = NA,
#                       cases = sum(cases),
#                       total = sum(total),
#                       aacases = sum(aacases),
#                       aaf = aacases/cases,
#                       totaacases = sum(totaacases),
#                       totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
#                                          TRUE ~ 0)) %>% 
#                unique()) -> pafsplit
#  
#  
#  setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_results")
#  if (dat == "glob_inc"){
#    write.csv(pafsmfp, "pafs_inc_country_03.11.20.csv", row.names=FALSE)
#    write.csv(pafsplit, "pafs_inc_country_split_03.11..20.csv", row.names=FALSE)
#  } else if (dat == "glob_mort") {
#    write.csv(pafsmfp, "pafs_mort_country_03.11.20.csv", row.names=FALSE)
#    write.csv(pafsplit, "pafs_mort_country_split_03.11.20.csv", row.names=FALSE)
#  }
#  
  
  #--- calc aa cases for WCRF CRC estimates 
  # bind second set of all cancer total excluding stomach and pancreas for main analysis results
  # sum total aa cases for all cancers and assign HCC aa cases to total liver aa cases
  pafs2 %>% 
    group_by(country_code, age_num, sex) %>% 
    mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases, na.rm=T), 
                               TRUE ~ aacases)) %>%
    bind_rows(pafs2 %>% 
                group_by(country_code, age_num, sex) %>% 
                mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases[Main_analysis==1], na.rm=T)),
                       cancer_label = case_when(cancer_code %in% c(39,40) ~ paste0(cancer_label,"_M"))) %>% # add label to main all cancer results
                filter(cancer_code %in% c(39,40))) %>% 
    mutate(aacases = case_when(cancer_code == 11 ~ sum(aacases[cancer_code==11.1], na.rm=T),
                               cancer_code == 6 ~ sum(aacases[cancer_code==6.1], na.rm=T),
                               TRUE ~ aacases),
           aaf = case_when(is.na(aaf) ~ aacases/cases, # fill in gaps for AAF
                           TRUE ~ aaf)) -> pafsmf2
  pafsplit2 %>% 
    mutate(level = case_when(str_detect(disease,"light") ~ 1,
                             str_detect(disease,"mod") ~ 2,
                             str_detect(disease,"heavy") ~ 3,
                             substr(disease,nchar(as.character(disease))-5,nchar(as.character(disease)))=="former" ~ 4,
                             substr(disease,nchar(as.character(disease))-7,nchar(as.character(disease)))=="formercd" ~ 5,
                             str_detect(disease,"low") ~ 6,
                             TRUE ~ 1)) %>% 
    bind_rows(pafsplit2 %>% filter(is.na(disease) | (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 2)) %>% 
    bind_rows(pafsplit2 %>% filter(is.na(disease) | (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 3)) %>% 
    bind_rows(pafsplit2 %>% filter(is.na(disease) | (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 4)) %>% 
    bind_rows(pafsplit2 %>% filter(is.na(disease) | (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 5)) %>% 
    bind_rows(pafsplit2 %>% filter(is.na(disease) | (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 6)) %>% 
    bind_rows(pafsplit2 %>% filter(cancer_code %in% c(7,8,9,13), str_detect(disease,"mod"), age_num>5) %>% 
                mutate(level = 1,
                       aaf=0,
                       aacases=0)) %>% 
    group_by(country_code, age_num, sex, level) %>% 
    mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases, na.rm=T), 
                               TRUE ~ aacases)) -> pafsplit2

  pafsplit2 %>% 
    bind_rows(pafsplit2 %>% 
                group_by(country_code, age_num, sex, level) %>% 
                mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases[Main_analysis==1], na.rm=T)),
                       cancer_label = case_when(cancer_code %in% c(39,40) ~ paste0(cancer_label,"_M"))) %>% # add label to main all cancer results
                filter(cancer_code %in% c(39,40))) %>% 
    group_by(country_code, age_num, sex, level) %>%
    mutate(aacases = case_when(cancer_code == 11 ~ sum(aacases[cancer_code==11.1], na.rm=T),
                               cancer_code == 6 ~ sum(aacases[cancer_code==6.1], na.rm=T),
                               TRUE ~ aacases),
           aaf = case_when(is.na(aaf) ~ aacases/cases, # fill in gaps for AAF
                           TRUE ~ aaf)) -> pafsplit2
  
  # create persons totals
  pafsmf2 %>% 
    group_by(country_code,age_num, cancer_label) %>% 
    mutate(sex = 0,
           cases = sum(cases, na.rm=T),
           total = sum(total, na.rm=T),
           aacases = sum(aacases, na.rm=T),
           aaf = aacases/cases) %>% unique() -> pafsp2
  # bind mf and persons
  pafsmfp2 <- rbind(pafsmf2, pafsp2)
  
  # create persons totals
  pafsplit2 %>% 
    group_by(country_code, ISO3, age_num, cancer_label, level) %>% 
    mutate(sex = 0,
           cases = sum(cases, na.rm=T),
           total = sum(total, na.rm=T),
           aacases = sum(aacases, na.rm=T),
           aaf = aacases/cases) %>% unique() -> pafsplitp2
  # bind mf and persons
  pafsplit2 <- rbind(pafsplit2, pafsplitp2)

  ####--- Calculate PAFs for WCRF estimates 
  # create pafs per sex, cancer type, and country
  pafsmfp2 %>%
    group_by(country_code, sex, cancer_label) %>% 
    mutate(totaacases = sum(aacases),
           totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
                              TRUE ~ 0)) -> pafsmfp2
  # create pafs per sex, cancer type, and country
  pafsplit2 %>%
    group_by(country_code, sex, cancer_label, level) %>% 
    mutate(totaacases = sum(aacases),
           totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
                              TRUE ~ 0)) -> pafsplit2
  rm(pafsmf2, pafsp2, pafsplitp2)
  
  names(countrymatch) <- tolower(names(countrymatch))
  pafsmfp2 <- merge(countrymatch, pafsmfp2, by = "country_code")
  pafsplit2 <- merge(countrymatch, pafsplit2, by = "country_code")
  
  
    pafsmfp2 %>% 
      bind_rows(pafsmfp2 %>% 
                  filter(cancer_code %in% c(3,5)) %>% 
                  group_by(country_label, sex, age_num) %>% 
                  mutate(cancer_label = "Pharynx",
                         disease= "Pharynx",
                         cancer_code = 41,
                         cases = sum(cases),
                         total = sum(total),
                         aacases = sum(aacases),
                         aaf = aacases/cases,
                         totaacases = sum(totaacases),
                         totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
                                            TRUE ~ 0)) %>% 
                  unique()) -> pafsmfp2
    
    pafsplit2 %>% 
      bind_rows(pafsplit2 %>% 
                  filter(cancer_code %in% c(3,5)) %>%  
                  group_by(country_label, sex, age_num, level) %>% 
                  mutate(cancer_label = "Pharynx",
                         disease= "Pharynx",
                         cancer_code = 41,
                         cases = sum(cases),
                         total = sum(total),
                         aacases = sum(aacases),
                         aaf = aacases/cases,
                         totaacases = sum(totaacases),
                         totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
                                            TRUE ~ 0)) %>% 
                  unique()) -> pafsplit2
  
  # write to csv pafs for all ages combined and age specific attributable cases
  # per country with regional groupings for WCRF CRC RRs
  setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_results")
  if (dat == "glob_inc"){
    write.csv(pafsmfp2, "pafs_inc_country_WCRF.csv", row.names=FALSE)
    write.csv(pafsplit2, "pafs_inc_country_split_WCRF.csv", row.names=FALSE)
  } else if (dat == "glob_mort") {
    write.csv(pafsmfp2, "pafs_mort_country_WCRF.csv", row.names=FALSE)
    write.csv(pafsplit2, "pafs_mort_country_split_WCRF.csv", row.names=FALSE)
  }
 
}


####--- 10g increase - attributable cases and PAF for each age group, sex and cancer per country ----
# for incidence and mortality
for (dat in c("glob_inc", "glob_mort")) {
  
  glob <- get(dat)
  
  # merge AAFs data with incidence data
  pafsplit <- merge(glob, aafsplit10, by = c("ISO3", "sex", "paf_age_cancer", "cancer_code"),allow.cartesian=TRUE, all = TRUE)
  
  
  ####--- Create age-specific alc attributable cases 
  
  # calculate aa cases per age group per cancer type
  pafsplit <- pafsplit[, aacases:= aaf*cases][age_num<6, aacases := 0][age_num<6, aaf := 0]
  
  # remove premenop results in older ages
  pafsplit <- pafsplit[(cancer_code==20 & age_num <6),disease:= "Breast"] 
  pafsplit <- pafsplit[(cancer_code==20 & age_num >10 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)),disease:= "Breast_postmen"]
  pafsplit <- pafsplit[(cancer_code==20 & age_num %in% c(6:10) & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728)),disease:= "Breast_premen"]
  pafsplit <- pafsplit[!(cancer_code==20 & str_detect(disease,"Breast_premen") & age_num > 10),] 
  # remove postmenop results in younger ages
  pafsplit <- pafsplit[!(cancer_code==20 & str_detect(disease,"Breast_postmen") & age_num <= 10),] 
  # remove NAs for male breast cancer, replace with 0
  pafsplit <- pafsplit[cancer_code==20 & sex==1, aacases:=0]
  
  #--- calc aa cases 
  # bind second set of all cancer total excluding stomach and pancreas for main analysis results
  # sum total aa cases for all cancers and assign HCC aa cases to total liver aa cases
  
  pafsplit %>% 
    mutate(level = case_when(substr(disease,nchar(as.character(disease))-1,nchar(as.character(disease)))=="_1" ~ 1,
                             str_detect(disease,"_2") ~ 2,
                             str_detect(disease,"_3") ~ 3,
                             str_detect(disease,"_4") ~ 4,
                             str_detect(disease,"_5") ~ 5,
                             str_detect(disease,"_6") ~ 6,
                             str_detect(disease,"_7") ~ 7,
                             str_detect(disease,"_8") ~ 8,
                             str_detect(disease,"_9") ~ 9,
                             str_detect(disease,"_10") ~ 10,
                             str_detect(disease,"_11") ~ 11,
                             str_detect(disease,"_12") ~ 12,
                             str_detect(disease,"_13") ~ 13,
                             str_detect(disease,"_14") ~ 14,
                             str_detect(disease,"_15") ~ 15,
                             TRUE ~ 1)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 2)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 3)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 4)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 5)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 6)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 7)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 8)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 9)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 10)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 11)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 12)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 13)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 14)) %>% 
    bind_rows(pafsplit %>% filter(is.na(disease)| (cancer_code==20 & country_code %in% c(254, 258, 275, 312, 316, 474, 540, 630, 638, 728))| (cancer_code==20 & age_num<6)) %>% 
                mutate(level = 15)) %>% 
    bind_rows(pafsplit %>% filter(cancer_code %in% c(7,8,9,13), str_detect(disease,"_13"), age_num>5) %>% 
                mutate(level = 1,
                       aaf=0,
                       aacases=0)) %>% 
    group_by(country_code, age_num, sex, level) %>% 
    mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases, na.rm=T), 
                               TRUE ~ aacases)) -> pafsplit
  
  pafsplit %>% 
    bind_rows(pafsplit %>% 
                group_by(country_code, age_num, sex, level) %>% 
                mutate(aacases = case_when(cancer_code %in% c(39,40) ~ sum(aacases[Main_analysis==1], na.rm=T)),
                       cancer_label = case_when(cancer_code %in% c(39,40) ~ paste0(cancer_label,"_M"))) %>% # add label to main all cancer results
                filter(cancer_code %in% c(39,40))) %>% 
    group_by(country_code, age_num, sex, level) %>%
    mutate(aacases = case_when(cancer_code == 11 ~ sum(aacases[cancer_code==11.1], na.rm=T),
                               cancer_code == 6 ~ sum(aacases[cancer_code==6.1], na.rm=T),
                               TRUE ~ aacases),
           aaf = case_when(is.na(aaf) ~ aacases/cases, # fill in gaps for AAF
                           TRUE ~ aaf)) -> pafsplit
  
  # create persons totals for drinking category data
  pafsplit %>% 
    bind_rows(pafsplit %>% group_by(country_code, ISO3, age_num, cancer_label, level) %>% 
                mutate(sex = 0,
                       cases = sum(cases, na.rm=T),
                       total = sum(total, na.rm=T),
                       aacases = sum(aacases, na.rm=T),
                       aaf = aacases/cases) %>% unique()) -> pafsplit

  ####--- Calculate PAFs  
  # create pafs per sex, cancer type, and country for drinking categories
  pafsplit %>%
    group_by(country_code, sex, cancer_label, level) %>% 
    mutate(totaacases = sum(aacases),
           totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
                              TRUE ~ 0)) %>% 
    dplyr::select(-disease, -Main_analysis, -paf_age_cancer) %>% 
    left_join(countrymatch) -> pafsplit
  
  pafsplit %>% 
    bind_rows(pafsplit %>% 
                filter(cancer_code %in% c(3,5)) %>%  
                group_by(country_label, sex, age_num, level) %>% 
                mutate(cancer_label = "Pharynx",
                       disease= "Pharynx",
                       cancer_code = 41,
                       cases = sum(cases),
                       total = sum(total),
                       aacases = sum(aacases),
                       aaf = aacases/cases,
                       totaacases = sum(totaacases),
                       totpaf = case_when(!is.na(totaacases/total) ~ totaacases/total,
                                          TRUE ~ 0)) %>% 
                unique()) -> pafsplit

  
  # write pafs and ASRs to csv for all ages combined and age specific attributable cases
  # per country
  setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_results")
  
  # save country data for WHO Europe countries
  if (dat == "glob_inc"){
    write.csv(pafsplit, "pafs_inc_country_split10.csv", row.names=FALSE)
  } else if (dat == "glob_mort") {
    write.csv(pafsplit, "pafs_mort_country_split10.csv", row.names=FALSE)
  }
  
}

rm(AAF_CANCER,AAF_CANCER_split, aafs, aafsplit, HCCcases, OSCCcases, pafs, DATA_GLOB, pafsplit, pafs2, pafsmfp2, pafsplit2)

#escc proportions for WHOEurope tax modeling study
glob_inc18 %>% 
  filter(cancer_code == 6.1, age == 10) %>% 
  select(-age,-cancer_code,-cases,-total) %>% 
  left_join(countrymatch %>% 
              select(country_code, web_country_label, who_region_label)) %>% 
  filter(who_region_label=="EURO") %>% 
  select(-who_region_label)-> osccEur
write.csv(osccEur, "//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/WHOEurope/osccEur.csv", row.names=FALSE)

