library(tidyverse)
library(data.table)

# load RRs data
setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/alc_data")
RRraw <- read.csv("WCRF_RRs3.csv")

# calculate beta values
RRraw %>% 
  mutate(betaCurrent = log(RR)/Gramsalc, # calculate beta values 
         covbetaCurrent = ((log(UCI)-log(LCI))/(2*1.96*Gramsalc))^2) -> RRbeta  # calculate covariance from CIs
rm(RRraw)

# tests to see which method to use for variance
#Var1 <- ((RRbeta$betaCurrent-log(RRbeta$LCI)/RRbeta$Gramsalc)/1.96)^2
#Var2 <- ((log(RRbeta$UCI)/RRbeta$Gramsalc-RRbeta$betaCurrent)/1.96)^2
#Var3 <- ((log(RRbeta$UCI)-log(RRbeta$LCI))/(2*1.96*RRbeta$Gramsalc))^2

# convert to data table
RRbeta_t <- data.table(RRbeta)

# create list of cancer types with linear trend
listlin <- c(as.character(unique(RRbeta_t$disease[RRbeta_t$Linear==1])))
# create list of RR functions
for (gender in c(1,2)){
  RRbeta <- RRbeta_t[sex==gender,]
  if (gender == 1){
    RRslistlin_m <- lapply(listlin, function(name){
      list(disease = name,
           RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x*log(x)*beta[4])}, # relative risk function for the (current) drinkers
           betaCurrent = c(0,RRbeta$betaCurrent[RRbeta$disease==name],0,0), # coefficient for the RR function of (current) drinkers
           lowerlim = RRbeta$Alc_start[RRbeta$disease==name], # where alcohol distribution should start from
           midlim = RRbeta$Alc_middle[RRbeta$disease==name],
           upperlim = RRbeta$Alc_end[RRbeta$disease==name],
           covBetaCurrent = matrix(c(0,0,0,0,0, RRbeta$covbetaCurrent[RRbeta$disease==name],0,0,0,0,0,0,0, 0,0,0),4,4),
           lnRRFormer = RRbeta$lnRRFormer[RRbeta$disease==name],
           varLnRRFormer = RRbeta$varLnRRFormer[RRbeta$disease==name]) 
    })
  } else if (gender == 2){
    
    RRslistlin_f <- lapply(listlin, function(name){
      list(disease = name,
           RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x*log(x)*beta[4])}, # relative risk function for the (current) drinkers
           betaCurrent = c(0,RRbeta$betaCurrent[RRbeta$disease==name],0,0), # coefficient for the RR function of (current) drinkers
           lowerlim = RRbeta$Alc_start[RRbeta$disease==name], # where alcohol distribution should start from
           midlim = RRbeta$Alc_middle[RRbeta$disease==name],
           upperlim = RRbeta$Alc_end[RRbeta$disease==name],
           covBetaCurrent = matrix(c(0,0,0,0,0, RRbeta$covbetaCurrent[RRbeta$disease==name],0,0,0,0,0,0,0, 0,0,0),4,4),
           lnRRFormer = RRbeta$lnRRFormer[RRbeta$disease==name],
           varLnRRFormer = RRbeta$varLnRRFormer[RRbeta$disease==name]) 
    }) 
  }
}

rm(RRbeta, RRbeta_t)

## non-linear oesoph SCC
Oesophagus_SCC_nonlin = list(disease = "Oesophagus_SCC_nonlin",
                             RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x*log(x)*beta[4])},
                             betaCurrent = c(0,0.05593,0,-0.00789),
                             lowerlim = 0.1,
                             midlim = 0.1,
                             upperlim = 150,
                             covBetaCurrent = matrix(c(0,0,0,0, 0,0.000065,0,-0.00001, 0,0,0,0, 0,-0.00001,0,0.00000264),4,4),
                             lnRRFormer = 1.16,
                             varLnRRFormer = 0.243480229040442^2)

Colon_lin = list(disease = "Colon_lin",
                 RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x^3*beta[4])},
                 betaCurrent = c(0, 0.006279, 0, 0),
                 lowerlim = 0.1,
                 midlim = 0.1,
                 upperlim = 150,
                 covBetaCurrent = matrix(c(0,0,0,0,0, 0.000000907,0,0,0,0,0,0,0, 0,0,0),4,4),
                 lnRRFormer = 2.19,
                 varLnRRFormer = 0.0465106^2)

Rectum_lin = list(disease = "Rectum_lin",
                  RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x^3*beta[4])},
                  betaCurrent = c(0, 0.006279, 0, 0),
                  lowerlim = 0.1,
                  midlim = 0.1,
                  upperlim = 150,
                  covBetaCurrent = matrix(c(0,0,0,0,0, 0.000000907,0,0,0,0,0,0,0, 0,0,0),4,4),
                  lnRRFormer = 2.19,
                  varLnRRFormer = 0.0465106^2)

RRslistnl_m <- list(Oesophagus_SCC_nonlin, Colon_lin, Rectum_lin)

Oesophagus_SCC_nonlin = list(disease = "Oesophagus_SCC_nonlin",
                             RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x*log(x)*beta[4])},
                             betaCurrent = c(0,0.05593,0,-0.00789),
                             lowerlim = 0.1,
                             midlim = 0.1,
                             upperlim = 150,
                             covBetaCurrent = matrix(c(0,0,0,0, 0,0.000065,0,-0.00001, 0,0,0,0, 0,-0.00001,0,0.00000264),4,4),
                             lnRRFormer = 1.16,
                             varLnRRFormer = 0.243480229040442^2)

Colon_lin = list(disease = "Colon_lin",
                 RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x^3*beta[4])},
                 betaCurrent = c(0, 0.006279, 0, 0),
                 lowerlim = 0.1,
                 midlim = 0.1,
                 upperlim = 150,
                 covBetaCurrent = matrix(c(0,0,0,0,0, 0.000000907,0,0,0,0,0,0,0, 0,0,0),4,4),
                 lnRRFormer = 1.05,
                 varLnRRFormer = 0.145968002587317^2)

Rectum_lin = list(disease = "Rectum_lin",
                  RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x^3*beta[4])},
                  betaCurrent = c(0, 0.006279, 0, 0),
                  lowerlim = 0.1,
                  midlim = 0.1,
                  upperlim = 150,
                  covBetaCurrent = matrix(c(0,0,0,0,0, 0.000000907,0,0,0,0,0,0,0, 0,0,0),4,4),
                  lnRRFormer = 1.05,
                  varLnRRFormer = 0.145968002587317^2)

RRslistnl_f <- list(Oesophagus_SCC_nonlin, Colon_lin, Rectum_lin)

# create list of linear and non-linear RRs
RRslist_m <- c(RRslistlin_m, RRslistnl_m)
RRslist_f <- c(RRslistlin_f, RRslistnl_f)
rm(RRslistlin_m, RRslistlin_f, RRslistnl_m, RRslistnl_f, Oesophagus_SCC_nonlin, Colon_lin, Rectum_lin)


setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs")

## EXAMPLE OF NON-LIN FUNCTION from Kevin's code
# copies of oesoph and CRC RRs including covariance matrix
####### Oesophagus SCC Cancer ###########
## male
#Oesophagus_SCC_cancer_male = list(disease = "Oesophagus_SCC_Cancer",
#                                  RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x*log(x)*beta[4])},
#                                  betaCurrent = c(0,0.05593,0,-0.00789),
#                                  covBetaCurrent = matrix(c(0,0,0,0, 0,0.000065,0,-0.00001, 0,0,0,0, 0,-0.00001,0,0.00000264),4,4),
#                                  lnRRFormer = log(1.16),
#                                  varLnRRFormer = 0.243480229040442^2)
## female
#Oesophagus_SCC_cancer_female = list(disease = "Oesophagus_SCC_Cancer",
#                                    RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x*log(x)*beta[4])},
#                                    betaCurrent = c(0,0.05593,0,-0.00789),
#                                    covBetaCurrent = matrix(c(0,0,0,0, 0,0.000065,0,-0.00001, 0,0,0,0, 0,-0.00001,0,0.00000264),4,4),
#                                    lnRRFormer = log(1.16),
#                                    varLnRRFormer = 0.243480229040442^2)
####### Colorectal cancer #######
## male
#colorectalcancer_male = list(disease = "Colorectal_Cancer",
#                             RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x^3*beta[4])},
#                             betaCurrent = c(0, 0.006279, 0, 0),
#                             covBetaCurrent = matrix(c(0,0,0,0,0, 0.000000907,0,0,0,0,0,0,0, 0,0,0),4,4),
#                             lnRRFormer = log(2.19),
#                             varLnRRFormer = 0.0465106^2)
## female
#colorectalcancer_female = list(disease = "Colorectal_Cancer",
#                               RRCurrent = function(x, beta) {exp(1*beta[1] + x*beta[2] + x^2*beta[3] + x^3*beta[4])},
#                               betaCurrent =c(0, 0.006279, 0, 0),
#                               covBetaCurrent = matrix(c(0,0,0,0,0, 0.000000907,0,0,0,0,0,0,0, 0,0,0),4,4),
#                               lnRRFormer = log(1.05),
#                               varLnRRFormer = 0.145968002587317^2)

