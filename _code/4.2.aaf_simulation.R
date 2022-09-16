################################################################
#### The program produces PE and simulated estiamtes PAFS   ####
#### each value is outputted into a seperate file   	    	####
################################################################

library(parallel)
library(MASS)

rm(list = ls())
nnn = 1000 # number of simulations
#year_sim = 2010


# added if/else where normalising constant = 0, replace with 1

##---- RUN THIS - SIMULATIONS CODE -----
PROC_SIM_FILES <- function(year_sim){
  
  #### READING IN DATA ####
  #### READING IN DATA ####
  setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/simulations")
  load("SIM_CANCER_RR_HR.RData") # RR_CANCER_SIM
  
  setwd("//Inti/cin/Users/RumgayH/alcohol project/IARC Cancer project/GLOBAL BURDEN OF CANCERS ATTRIBUTABLE TO ALCOHOL/DATA_2010")
  drkData 			<- read.csv("DATA_INPUT_APC_LANCET_CAPPED_2019 - v2.csv")
  load(paste("DRINKING_STATUS_SIM_", year_sim, ".RData", sep="")); DRINKING_STATUS_SIM <- DRINKING_STATUS_C_SIM
  load(paste("K_SIM_", year_sim, ".RData", sep="")); K_SIM <- K_SIM
  load(paste("THETA_SIM_", year_sim, ".RData", sep="")); THETA_SIM <- THETA_SIM
  load(paste("NORM_CONS_SIM_", year_sim, ".RData", sep="")); NORM_CONS_SIM <- NORM_CONS_SIM
  
  DRINKING_STATUS_SIM_T <- DRINKING_STATUS_SIM
  K_SIM_T <- K_SIM
  THETA_SIM_T <- THETA_SIM
  NORM_CONS_SIM_T <- NORM_CONS_SIM
  RR_TEMP_T <- RR_CANCER_SIM
  YEAR <- year_sim
  BATCH <- "CANCER"
  
  REGIONS_LIST 	 	<- drkData$ISO_3[drkData$SEX == 1 & drkData$AGE_CATEGORY == 1 & drkData$year == year_sim]
  
  SIM_TEMP 				<- as.data.frame( drkData$ISO_3[drkData$SEX == 1 & drkData$AGE_CATEGORY == 1] ) 
  colnames(SIM_TEMP)		<- "ISO_3"
  SIM_TEMP$RUSSIA			<- NA
  SIM_TEMP$AGE_CATEGORY 	<- SIM_TEMP$SEX <- NA
  COL_A 					<- length(SIM_TEMP[1,])+1
  COL_B 					<- length(SIM_TEMP[1,])+nnn

  
  ### PREDEFINING FUNCTIONS ####
  ### PREDEFINING FUNCTIONS ####
  
  ## risk of drinkers
  #drkfun <- function(x, pabs, pform, nc, k, theta)
  #{
  #  (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta) * dis$RRCurrent(x, beta = dis$betaCurrent)
  #}

  #### SIMULATIONS ####
  #### SIMULATIONS ####
  
  PAF_NORMAL <- function(SEX, AGE){
    
    RR_TEMP			 <- RR_TEMP_T[[SEX]]
    pabs			 <- DRINKING_STATUS_SIM_T[[1]][[SEX]][[AGE]]
    pform			 <- DRINKING_STATUS_SIM_T[[2]][[SEX]][[AGE]]
    cd				 <- DRINKING_STATUS_SIM_T[[3]][[SEX]][[AGE]]
    # pabs[, (COL_A : COL_B) ] + pform[, (COL_A : COL_B) ] + cd[, (COL_A : COL_B) ]
    
    k				 <- K_SIM_T[[SEX]]
    theta			 <- THETA_SIM_T[[SEX]][[AGE]]
    NORM_CONS_SIM	 <- NORM_CONS_SIM_T[[SEX]][[AGE]]
    
    REGIONS_LIST 	  <- unique(pabs$ISO_3)
    
    aaf_out_C <- lapply(REGIONS_LIST, function(z1c){
      
      pabs_B			 <- pabs[pabs$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      pform_B			 <- pform[pform$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      k_B				 <- k[k$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      theta_B			 <- theta[theta$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      NORM_CONS_SIM_B	 <- NORM_CONS_SIM[NORM_CONS_SIM$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      #NORM_CONS_SIM_B	 <- NORM_CONS_SIM_T[,5:1004]
      
      aaf_out_B <- lapply(RR_TEMP, function(dis){
        
        cat("Region:", toString(z1c), "; disease:", dis$disease, "\n")
        
        aaf <- lapply(1:nnn, function(nnn2) {
          
          #if (NORM_CONS_SIM_B[nnn2][[1]]==0){
          #  normConst <- 1
          #} else {
          #  normConst <- NORM_CONS_SIM_B[nnn2][[1]]
          #}
          
#          ## risk of drinkers
          drkfun2 <- function(x, pabs, pform, nc, k, theta)
          {
            (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta)  * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
          } # set shape of gamma distribution and risk at each alcohol level
          
          
          drk1 <- mapply(function(pabs, pform, nc, k, theta)
          {
            integrate(drkfun2, lower = dis$lowerlim, upper = dis$midlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta, stop.on.error = FALSE)$value
          }, # compute population risk at each drinking level considering the proportion of drinkers
          pabs = pabs_B[nnn2][[1]], pform = pform_B[nnn2][[1]], nc = NORM_CONS_SIM_B[nnn2][[1]], k = k_B[nnn2][[1]], theta = theta_B[nnn2][[1]]
          ) # define abstainers, former drinkers, normalisation constant, k and theta variables
          
          drk2 <- mapply(function(pabs, pform, nc, k, theta)
          {
            integrate(drkfun2, lower = dis$midlim, upper = dis$upperlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta, stop.on.error = FALSE)$value
          }, # compute population risk at each drinking level considering the proportion of drinkers
          pabs = pabs_B[nnn2][[1]], pform = pform_B[nnn2][[1]], nc = NORM_CONS_SIM_B[nnn2][[1]], k = k_B[nnn2][[1]], theta = theta_B[nnn2][[1]]
          ) # define abstainers, former drinkers, normalisation constant, k and theta variables
          
          
          ## APPLYING FUNCTION
#          drk1 <- integrate(function(x){
#            (1 - pabs_B[nnn2][[1]] - pform_B[nnn2][[1]])/ nc * dgamma(x, shape = k_B[nnn2][[1]],scale = theta_B[nnn2][[1]]) * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
#          }, lower = dis$lowerlim, upper = dis$midlim, stop.on.error = FALSE)$value
#         
#          drk2 <- integrate(function(x){
#            (1 - pabs_B[nnn2][[1]] - pform_B[nnn2][[1]])/ nc * dgamma(x, shape = k_B[nnn2][[1]],scale = theta_B[nnn2][[1]]) * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
#          }, lower = dis$midlim, upper = dis$upperlim, stop.on.error = FALSE)$value
          
          ## AAF
          
          aaf <- (drk2) / (drk1 + drk2 + 1)
          
          aaf    
          
        })
        
        aaf_2 	<- do.call(cbind, aaf)
        aaf_out <- cbind.data.frame(toString(z1c), toString(dis$disease),(aaf_2)); aaf_out
        
      })
      
      aaf_out_B <- do.call(rbind, aaf_out_B); aaf_out_B} 	)
    
    aaf_out_PAF <- do.call(rbind, aaf_out_C)
    
    setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/simulations")
    
    PAF_name <- paste(paste(BATCH, SEX, AGE, YEAR, sep="_"), "new.RData", sep="")
    save(aaf_out_PAF, file = toString(PAF_name))	
  }
  
  PAF_NORMAL(SEX = 1, AGE = 1); PAF_NORMAL(SEX = 1, AGE = 2); PAF_NORMAL(SEX = 1, AGE = 3);
  PAF_NORMAL(SEX = 1, AGE = 4); PAF_NORMAL(SEX = 1, AGE = 5); PAF_NORMAL(SEX = 1, AGE = 6);
  
  PAF_NORMAL(SEX = 2, AGE = 1); PAF_NORMAL(SEX = 2, AGE = 2); PAF_NORMAL(SEX = 2, AGE = 3);
  PAF_NORMAL(SEX = 2, AGE = 4); PAF_NORMAL(SEX = 2, AGE = 5); PAF_NORMAL(SEX = 2, AGE = 6);
  
}

PROC_SIM_FILES(year_sim = 2010)

##---- RUN THIS - SIMULATIONS CODE for levels -  ----
PROC_SIM_FILES <- function(year_sim){
  
  #### READING IN DATA ####
  #### READING IN DATA ####
  setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/simulations")
  load("SIM_CANCER_RR_split_HR.RData") # RR_CANCER_SIM_SPLIT
  
  setwd("//Inti/cin/Users/RumgayH/alcohol project/IARC Cancer project/GLOBAL BURDEN OF CANCERS ATTRIBUTABLE TO ALCOHOL/DATA_2010")
  drkData 			<- read.csv("DATA_INPUT_APC_LANCET_CAPPED_2019 - v2.csv")
  load(paste("DRINKING_STATUS_SIM_", year_sim, ".RData", sep="")); DRINKING_STATUS_SIM <- DRINKING_STATUS_C_SIM
  load(paste("K_SIM_", year_sim, ".RData", sep="")); K_SIM <- K_SIM
  load(paste("THETA_SIM_", year_sim, ".RData", sep="")); THETA_SIM <- THETA_SIM
  load(paste("NORM_CONS_SIM_", year_sim, ".RData", sep="")); NORM_CONS_SIM <- NORM_CONS_SIM
  
  DRINKING_STATUS_SIM_T <- DRINKING_STATUS_SIM
  K_SIM_T <- K_SIM
  THETA_SIM_T <- THETA_SIM
  NORM_CONS_SIM_T <- NORM_CONS_SIM
  RR_TEMP_T <- RR_CANCER_SIM_SPLIT
  YEAR <- year_sim
  BATCH <- "CANCER"
  
  REGIONS_LIST 	 	<- drkData$ISO_3[drkData$SEX == 1 & drkData$AGE_CATEGORY == 1 & drkData$year == 2010]
  
  SIM_TEMP 				<- as.data.frame( drkData$ISO_3[drkData$SEX == 1 & drkData$AGE_CATEGORY == 1] ) 
  colnames(SIM_TEMP)		<- "ISO_3"
  SIM_TEMP$RUSSIA			<- NA
  SIM_TEMP$AGE_CATEGORY 	<- SIM_TEMP$SEX <- NA
  COL_A 					<- length(SIM_TEMP[1,])+1
  COL_B 					<- length(SIM_TEMP[1,])+nnn
  
  ### PREDEFINING FUNCTIONS ####
  ### PREDEFINING FUNCTIONS ####
  
  ## risk of drinkers
  #drkfun <- function(x, pabs, pform, nc, k, theta)
  #{
  #  (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta) * dis$RRCurrent(x, beta = dis$betaCurrent)
  #}
  
  #### SIMULATIONS ####
  #### SIMULATIONS ####
  
  PAF_NORMAL <- function(SEX, AGE){
    
    RR_TEMP			 <- RR_TEMP_T[[SEX]]
    pabs			 <- DRINKING_STATUS_SIM_T[[1]][[SEX]][[AGE]]
    pform			 <- DRINKING_STATUS_SIM_T[[2]][[SEX]][[AGE]]
    cd				 <- DRINKING_STATUS_SIM_T[[3]][[SEX]][[AGE]]
    # pabs[, (COL_A : COL_B) ] + pform[, (COL_A : COL_B) ] + cd[, (COL_A : COL_B) ]
    
    k				 <- K_SIM_T[[SEX]]
    theta			 <- THETA_SIM_T[[SEX]][[AGE]]
    NORM_CONS_SIM	 <- NORM_CONS_SIM_T[[SEX]][[AGE]]
    
    REGIONS_LIST 	  <- unique(pabs$ISO_3)
    
    aaf_out_C <- lapply(REGIONS_LIST, function(z1c){
      
      pabs_B			 <- pabs[pabs$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      pform_B			 <- pform[pform$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      k_B				 <- k[k$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      theta_B			 <- theta[theta$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      NORM_CONS_SIM_B	 <- NORM_CONS_SIM[NORM_CONS_SIM$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      
      
      aaf_out_B <- lapply(RR_TEMP, function(dis){
        
        cat("Region:", toString(z1c), "; disease:", dis$disease, "\n")
        
        aaf <- lapply(1:nnn, function(nnn2) {
          
          #if (NORM_CONS_SIM_B[nnn2][[1]]==0){
          #  normConst <- 1
          #} else {
          #  normConst <- NORM_CONS_SIM_B[nnn2][[1]]
          #}
          
          #          ## risk of drinkers
          drkfun2 <- function(x, pabs, pform, nc, k, theta)
          {
            (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta)  * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
          } # set shape of gamma distribution and risk at each alcohol level
          
          
          drk1 <- mapply(function(pabs, pform, nc, k, theta)
          {
            integrate(drkfun2, lower = dis$lowerlim, upper = dis$midlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta, stop.on.error = FALSE)$value
          }, # compute population risk at each drinking level considering the proportion of drinkers
          pabs = pabs_B[nnn2][[1]], pform = pform_B[nnn2][[1]], nc = NORM_CONS_SIM_B[nnn2][[1]], k = k_B[nnn2][[1]], theta = theta_B[nnn2][[1]]
          ) # define abstainers, former drinkers, normalisation constant, k and theta variables
          
          drk2 <- mapply(function(pabs, pform, nc, k, theta)
          {
            integrate(drkfun2, lower = dis$midlim, upper = dis$upperlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta, stop.on.error = FALSE)$value
          }, # compute population risk at each drinking level considering the proportion of drinkers
          pabs = pabs_B[nnn2][[1]], pform = pform_B[nnn2][[1]], nc = NORM_CONS_SIM_B[nnn2][[1]], k = k_B[nnn2][[1]], theta = theta_B[nnn2][[1]]
          ) # define abstainers, former drinkers, normalisation constant, k and theta variables
          
          drk3 <- mapply(function(pabs, pform, nc, k, theta)
          {
            integrate(drkfun2, lower = dis$upperlim, upper = dis$toplim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta, stop.on.error = FALSE)$value
          }, # compute population risk at each drinking level considering the proportion of drinkers
          pabs = pabs_B[nnn2][[1]], pform = pform_B[nnn2][[1]], nc = NORM_CONS_SIM_B[nnn2][[1]], k = k_B[nnn2][[1]], theta = theta_B[nnn2][[1]]
          ) # define abstainers, former drinkers, normalisation constant, k and theta variables
          
        ## APPLYING FUNCTION
#          drk1 <- integrate(function(x){
#            (1 - pabs_B[nnn2][[1]] - pform_B[nnn2][[1]])/ ifelse(NORM_CONS_SIM_B[nnn2][[1]]==0,1,NORM_CONS_SIM_B[nnn2][[1]]) * dgamma(x, shape = k_B[nnn2][[1]],scale = theta_B[nnn2][[1]]) * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
#          }, lower = dis$lowerlim, upper = dis$midlim, stop.on.error = FALSE)$value
#          
#          drk2 <- integrate(function(x){
#            (1 - pabs_B[nnn2][[1]] - pform_B[nnn2][[1]])/ ifelse(NORM_CONS_SIM_B[nnn2][[1]]==0,1,NORM_CONS_SIM_B[nnn2][[1]]) * dgamma(x, shape = k_B[nnn2][[1]],scale = theta_B[nnn2][[1]]) * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
#          }, lower = dis$midlim, upper = dis$upperlim, stop.on.error = FALSE)$value
#          
#          drk3 <- integrate(function(x){
#            (1 - pabs_B[nnn2][[1]] - pform_B[nnn2][[1]])/ ifelse(NORM_CONS_SIM_B[nnn2][[1]]==0,1,NORM_CONS_SIM_B[nnn2][[1]]) * dgamma(x, shape = k_B[nnn2][[1]],scale = theta_B[nnn2][[1]]) * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
#          }, lower = dis$upperlim, upper = dis$toplim, stop.on.error = FALSE)$value
          
          
          ## AAF
          
          aaf <- (drk2) / (drk1 + drk2 + drk3 + 1)
          
          aaf
        })
        
        aaf_2 	<- do.call(cbind, aaf)
        aaf_out <- cbind.data.frame(toString(z1c), toString(dis$disease),(aaf_2)); aaf_out
        
      })
      
      aaf_out_B <- do.call(rbind, aaf_out_B); aaf_out_B} 	)
    
    aaf_out_PAF <- do.call(rbind, aaf_out_C)
    
    setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/simulations")
    
    PAF_name <- paste(paste(BATCH, SEX, AGE, YEAR, sep="_"), "split_new.RData", sep="")
    save(aaf_out_PAF, file = toString(PAF_name))	
  }
  
  PAF_NORMAL(SEX = 1, AGE = 1); PAF_NORMAL(SEX = 1, AGE = 2); PAF_NORMAL(SEX = 1, AGE = 3);
  PAF_NORMAL(SEX = 1, AGE = 4); PAF_NORMAL(SEX = 1, AGE = 5); PAF_NORMAL(SEX = 1, AGE = 6);
  
  PAF_NORMAL(SEX = 2, AGE = 1); PAF_NORMAL(SEX = 2, AGE = 2); PAF_NORMAL(SEX = 2, AGE = 3);
  PAF_NORMAL(SEX = 2, AGE = 4); PAF_NORMAL(SEX = 2, AGE = 5); PAF_NORMAL(SEX = 2, AGE = 6);
  
}

PROC_SIM_FILES(year_sim = 2010)

##--- RUN THIS - SIMULATIONS CODE FORMER DRINKERS ----

PROC_SIM_FILES <- function(year_sim){
  
  #### READING IN DATA ####
  #### READING IN DATA ####
  setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/simulations")
  load("SIM_CANCER_RR_HR.RData") # RR_CANCER_SIM
  
  setwd("//Inti/cin/Users/RumgayH/alcohol project/IARC Cancer project/GLOBAL BURDEN OF CANCERS ATTRIBUTABLE TO ALCOHOL/DATA_2010")
  drkData 			<- read.csv("DATA_INPUT_APC_LANCET_CAPPED_2019 - v2.csv")
  load(paste("DRINKING_STATUS_SIM_", year_sim, ".RData", sep="")); DRINKING_STATUS_SIM <- DRINKING_STATUS_C_SIM
  load(paste("K_SIM_", year_sim, ".RData", sep="")); K_SIM <- K_SIM
  load(paste("THETA_SIM_", year_sim, ".RData", sep="")); THETA_SIM <- THETA_SIM
  load(paste("NORM_CONS_SIM_", year_sim, ".RData", sep="")); NORM_CONS_SIM <- NORM_CONS_SIM
  
  DRINKING_STATUS_SIM_T <- DRINKING_STATUS_SIM
  K_SIM_T <- K_SIM
  THETA_SIM_T <- THETA_SIM
  NORM_CONS_SIM_T <- NORM_CONS_SIM
  RR_TEMP_T <- RR_CANCER_SIM
  YEAR <- year_sim
  BATCH <- "CANCER"
  
  REGIONS_LIST 	 	<- drkData$ISO_3[drkData$SEX == 1 & drkData$AGE_CATEGORY == 1 & drkData$year == 2010]
  
  SIM_TEMP 				<- as.data.frame( drkData$ISO_3[drkData$SEX == 1 & drkData$AGE_CATEGORY == 1] ) 
  colnames(SIM_TEMP)		<- "ISO_3"
  SIM_TEMP$RUSSIA			<- NA
  SIM_TEMP$AGE_CATEGORY 	<- SIM_TEMP$SEX <- NA
  COL_A 					<- length(SIM_TEMP[1,])+1
  COL_B 					<- length(SIM_TEMP[1,])+nnn
  
  
  ### PREDEFINING FUNCTIONS ####
  ### PREDEFINING FUNCTIONS ####
  
  ## risk of drinkers
  #drkfun <- function(x, pabs, pform, nc, k, theta)
  #{
  #  (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta) * dis$RRCurrent(x, beta = dis$betaCurrent)
  #}
  
  #### SIMULATIONS ####
  #### SIMULATIONS ####
  
  PAF_NORMAL <- function(SEX, AGE){
    
    RR_TEMP			 <- RR_TEMP_T[[SEX]]
    pabs			 <- DRINKING_STATUS_SIM_T[[1]][[SEX]][[AGE]]
    pform			 <- DRINKING_STATUS_SIM_T[[2]][[SEX]][[AGE]]
    cd				 <- DRINKING_STATUS_SIM_T[[3]][[SEX]][[AGE]]
    # pabs[, (COL_A : COL_B) ] + pform[, (COL_A : COL_B) ] + cd[, (COL_A : COL_B) ]
    
    k				 <- K_SIM_T[[SEX]]
    theta			 <- THETA_SIM_T[[SEX]][[AGE]]
    NORM_CONS_SIM	 <- NORM_CONS_SIM_T[[SEX]][[AGE]]
    
    REGIONS_LIST 	  <- unique(pabs$ISO_3)
    
    aaf_out_C <- lapply(REGIONS_LIST, function(z1c){
      
      pabs_B			 <- pabs[pabs$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      pform_B			 <- pform[pform$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      k_B				 <- k[k$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      theta_B			 <- theta[theta$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      NORM_CONS_SIM_B	 <- NORM_CONS_SIM[NORM_CONS_SIM$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      
      
      aaf_out_B <- lapply(RR_TEMP, function(dis){
        
        cat("Region:", toString(z1c), "; disease:", dis$disease, "\n")
        
        aaf <- lapply(1:nnn, function(nnn2) {
          
#          if (NORM_CONS_SIM_B[nnn2][[1]]==0){
#            normConst <- 1
#          } else {
#            normConst <- NORM_CONS_SIM_B[nnn2][[1]]
#          }
          
          #          ## risk of drinkers
          drkfun2 <- function(x, pabs, pform, nc, k, theta)
          {
            (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta)  * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
          } # set shape of gamma distribution and risk at each alcohol level
          
          
          drk1 <- mapply(function(pabs, pform, nc, k, theta)
          {
            integrate(drkfun2, lower = dis$lowerlim, upper = dis$midlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta, stop.on.error = FALSE)$value
          }, # compute population risk at each drinking level considering the proportion of drinkers
          pabs = pabs_B[nnn2][[1]], pform = pform_B[nnn2][[1]], nc = NORM_CONS_SIM_B[nnn2][[1]], k = k_B[nnn2][[1]], theta = theta_B[nnn2][[1]]
          ) # define abstainers, former drinkers, normalisation constant, k and theta variables
          
          drk2 <- mapply(function(pabs, pform, nc, k, theta)
          {
            integrate(drkfun2, lower = dis$midlim, upper = dis$upperlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta, stop.on.error = FALSE)$value
          }, # compute population risk at each drinking level considering the proportion of drinkers
          pabs = pabs_B[nnn2][[1]], pform = pform_B[nnn2][[1]], nc = NORM_CONS_SIM_B[nnn2][[1]], k = k_B[nnn2][[1]], theta = theta_B[nnn2][[1]]
          ) # define abstainers, former drinkers, normalisation constant, k and theta variables
          
          ## APPLYING FUNCTION
#          drk1 <- integrate(function(x){
#            (1 - pabs_B[nnn2][[1]] - pform_B[nnn2][[1]])/ ifelse(NORM_CONS_SIM_B[nnn2][[1]]==0,1,NORM_CONS_SIM_B[nnn2][[1]]) * dgamma(x, shape = k_B[nnn2][[1]],scale = theta_B[nnn2][[1]]) * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
#          }, lower = dis$lowerlim, upper = dis$midlim, stop.on.error = FALSE)$value
#          
#          drk2 <- integrate(function(x){
#            (1 - pabs_B[nnn2][[1]] - pform_B[nnn2][[1]])/ ifelse(NORM_CONS_SIM_B[nnn2][[1]]==0,1,NORM_CONS_SIM_B[nnn2][[1]]) * dgamma(x, shape = k_B[nnn2][[1]],scale = theta_B[nnn2][[1]]) * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
#          }, lower = dis$midlim, upper = dis$upperlim, stop.on.error = FALSE)$value
          
          ## AAF
          
          aaf <- (pform_B[nnn2][[1]]*(dis$lnRRFormer_sim[nnn2]-1)) / (pform_B[nnn2][[1]]*(dis$lnRRFormer_sim[nnn2]-1) + drk1 + drk2 + 1)
          
          aaf    
          
        })
        
        aaf_2 	<- do.call(cbind, aaf)
        aaf_out <- cbind.data.frame(toString(z1c), toString(dis$disease),(aaf_2)); aaf_out
        
      })
      
      aaf_out_B <- do.call(rbind, aaf_out_B); aaf_out_B} 	)
    
    aaf_out_PAF <- do.call(rbind, aaf_out_C)
    
    setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/simulations")
    
    PAF_name <- paste(paste(BATCH, SEX, AGE, YEAR, sep="_"), "former.RData", sep="")
    save(aaf_out_PAF, file = toString(PAF_name))	
  }
  
  PAF_NORMAL(SEX = 1, AGE = 1); PAF_NORMAL(SEX = 1, AGE = 2); PAF_NORMAL(SEX = 1, AGE = 3);
  PAF_NORMAL(SEX = 1, AGE = 4); PAF_NORMAL(SEX = 1, AGE = 5); PAF_NORMAL(SEX = 1, AGE = 6);
  
  PAF_NORMAL(SEX = 2, AGE = 1); PAF_NORMAL(SEX = 2, AGE = 2); PAF_NORMAL(SEX = 2, AGE = 3);
  PAF_NORMAL(SEX = 2, AGE = 4); PAF_NORMAL(SEX = 2, AGE = 5); PAF_NORMAL(SEX = 2, AGE = 6);
  
}

PROC_SIM_FILES(year_sim = 2010)


##--- RUN THIS - SIMULATIONS CODE FOR CURRENT DRINKERS INCLUDING FORMER DRINKERS IN DENOM ----

PROC_SIM_FILES <- function(year_sim){
  
  #### READING IN DATA ####
  #### READING IN DATA ####
  setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/simulations")
  load("SIM_CANCER_RR_HR.RData") # RR_CANCER_SIM
  
  setwd("//Inti/cin/Users/RumgayH/alcohol project/IARC Cancer project/GLOBAL BURDEN OF CANCERS ATTRIBUTABLE TO ALCOHOL/DATA_2010")
  drkData 			<- read.csv("DATA_INPUT_APC_LANCET_CAPPED_2019 - v2.csv")
  load(paste("DRINKING_STATUS_SIM_", year_sim, ".RData", sep="")); DRINKING_STATUS_SIM <- DRINKING_STATUS_C_SIM
  load(paste("K_SIM_", year_sim, ".RData", sep="")); K_SIM <- K_SIM
  load(paste("THETA_SIM_", year_sim, ".RData", sep="")); THETA_SIM <- THETA_SIM
  load(paste("NORM_CONS_SIM_", year_sim, ".RData", sep="")); NORM_CONS_SIM <- NORM_CONS_SIM
  
  DRINKING_STATUS_SIM_T <- DRINKING_STATUS_SIM
  K_SIM_T <- K_SIM
  THETA_SIM_T <- THETA_SIM
  NORM_CONS_SIM_T <- NORM_CONS_SIM
  RR_TEMP_T <- RR_CANCER_SIM
  YEAR <- year_sim
  BATCH <- "CANCER"
  
  REGIONS_LIST 	 	<- drkData$ISO_3[drkData$SEX == 1 & drkData$AGE_CATEGORY == 1 & drkData$year == 2010]
  
  SIM_TEMP 				<- as.data.frame( drkData$ISO_3[drkData$SEX == 1 & drkData$AGE_CATEGORY == 1] ) 
  colnames(SIM_TEMP)		<- "ISO_3"
  SIM_TEMP$RUSSIA			<- NA
  SIM_TEMP$AGE_CATEGORY 	<- SIM_TEMP$SEX <- NA
  COL_A 					<- length(SIM_TEMP[1,])+1
  COL_B 					<- length(SIM_TEMP[1,])+nnn
  
  
  ### PREDEFINING FUNCTIONS ####
  ### PREDEFINING FUNCTIONS ####
  
  ## risk of drinkers
  #drkfun <- function(x, pabs, pform, nc, k, theta)
  #{
  #  (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta) * dis$RRCurrent(x, beta = dis$betaCurrent)
  #}
  
  #### SIMULATIONS ####
  #### SIMULATIONS ####
  
  PAF_NORMAL <- function(SEX, AGE){
    
    RR_TEMP			 <- RR_TEMP_T[[SEX]]
    pabs			 <- DRINKING_STATUS_SIM_T[[1]][[SEX]][[AGE]]
    pform			 <- DRINKING_STATUS_SIM_T[[2]][[SEX]][[AGE]]
    cd				 <- DRINKING_STATUS_SIM_T[[3]][[SEX]][[AGE]]
    # pabs[, (COL_A : COL_B) ] + pform[, (COL_A : COL_B) ] + cd[, (COL_A : COL_B) ]
    
    k				 <- K_SIM_T[[SEX]]
    theta			 <- THETA_SIM_T[[SEX]][[AGE]]
    NORM_CONS_SIM	 <- NORM_CONS_SIM_T[[SEX]][[AGE]]
    
    REGIONS_LIST 	  <- unique(pabs$ISO_3)
    
    aaf_out_C <- lapply(REGIONS_LIST, function(z1c){
      
      pabs_B			 <- pabs[pabs$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      pform_B			 <- pform[pform$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      k_B				 <- k[k$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      theta_B			 <- theta[theta$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      NORM_CONS_SIM_B	 <- NORM_CONS_SIM[NORM_CONS_SIM$ISO_3 %in% toString(z1c), (COL_A : COL_B) ]
      
      
      aaf_out_B <- lapply(RR_TEMP, function(dis){
        
        cat("Region:", toString(z1c), "; disease:", dis$disease, "\n")
        
        aaf <- lapply(1:nnn, function(nnn2) {
          
          #if (NORM_CONS_SIM_B[nnn2][[1]]==0){
          #  normConst <- 1
          #} else {
          #  normConst <- NORM_CONS_SIM_B[nnn2][[1]]
          #}
          
          #          ## risk of drinkers
          drkfun2 <- function(x, pabs, pform, nc, k, theta)
          {
            (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta)  * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
          } # set shape of gamma distribution and risk at each alcohol level
          
          
          drk1 <- mapply(function(pabs, pform, nc, k, theta)
          {
            integrate(drkfun2, lower = dis$lowerlim, upper = dis$midlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta, stop.on.error = FALSE)$value
          }, # compute population risk at each drinking level considering the proportion of drinkers
          pabs = pabs_B[nnn2][[1]], pform = pform_B[nnn2][[1]], nc = normConst, k = k_B[nnn2][[1]], theta = theta_B[nnn2][[1]]
          ) # define abstainers, former drinkers, normalisation constant, k and theta variables
          
          drk2 <- mapply(function(pabs, pform, nc, k, theta)
          {
            integrate(drkfun2, lower = dis$midlim, upper = dis$upperlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta, stop.on.error = FALSE)$value
          }, # compute population risk at each drinking level considering the proportion of drinkers
          pabs = pabs_B[nnn2][[1]], pform = pform_B[nnn2][[1]], nc = normConst, k = k_B[nnn2][[1]], theta = theta_B[nnn2][[1]]
          ) # define abstainers, former drinkers, normalisation constant, k and theta variables
          
          ## APPLYING FUNCTION
#          drk1 <- integrate(function(x){
#            (1 - pabs_B[nnn2][[1]] - pform_B[nnn2][[1]])/ ifelse(NORM_CONS_SIM_B[nnn2][[1]]==0,1,NORM_CONS_SIM_B[nnn2][[1]]) * dgamma(x, shape = k_B[nnn2][[1]],scale = theta_B[nnn2][[1]]) * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
#          }, lower = dis$lowerlim, upper = dis$midlim, stop.on.error = FALSE)$value
#          
#          drk2 <- integrate(function(x){
#            (1 - pabs_B[nnn2][[1]] - pform_B[nnn2][[1]])/ ifelse(NORM_CONS_SIM_B[nnn2][[1]]==0,1,NORM_CONS_SIM_B[nnn2][[1]]) * dgamma(x, shape = k_B[nnn2][[1]],scale = theta_B[nnn2][[1]]) * (dis$RRCurrent(x, beta = dis$beta_sim[nnn2,])-1)
#          }, lower = dis$midlim, upper = dis$upperlim, stop.on.error = FALSE)$value
          
          ## AAF
          
          aaf <- (drk2) / (pform_B[nnn2][[1]]*(dis$lnRRFormer_sim[nnn2]-1) + drk1 + drk2 + 1)
          
          aaf    
          
        })
        
        aaf_2 	<- do.call(cbind, aaf)
        aaf_out <- cbind.data.frame(toString(z1c), toString(dis$disease),(aaf_2)); aaf_out
        
      })
      
      aaf_out_B <- do.call(rbind, aaf_out_B); aaf_out_B} 	)
    
    aaf_out_PAF <- do.call(rbind, aaf_out_C)
    
    setwd("//Inti/cin/Users/RumgayH/RProjects/AlcPAFs/simulations")
    
    PAF_name <- paste(paste(BATCH, SEX, AGE, YEAR, sep="_"), "formercd.RData", sep="")
    save(aaf_out_PAF, file = toString(PAF_name))	
  }
  
  PAF_NORMAL(SEX = 1, AGE = 1); PAF_NORMAL(SEX = 1, AGE = 2); PAF_NORMAL(SEX = 1, AGE = 3);
  PAF_NORMAL(SEX = 1, AGE = 4); PAF_NORMAL(SEX = 1, AGE = 5); PAF_NORMAL(SEX = 1, AGE = 6);
  
  PAF_NORMAL(SEX = 2, AGE = 1); PAF_NORMAL(SEX = 2, AGE = 2); PAF_NORMAL(SEX = 2, AGE = 3);
  PAF_NORMAL(SEX = 2, AGE = 4); PAF_NORMAL(SEX = 2, AGE = 5); PAF_NORMAL(SEX = 2, AGE = 6);
  
}

PROC_SIM_FILES(year_sim = 2010)

