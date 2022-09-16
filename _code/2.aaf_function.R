## Source file to calculate adjusted mean alcohol consumption and alcohol-attributable fractions

## computeMeanStratum: computes the mean alcohol consumption for each age group
## in a given region and for a given sex. adjusting for spillage and unrecorded
## alcohol consumption (0.8). 
## popsize: Vector of X. Population size in the X age groups.
## relcoef: Vector of X. Relative coefficients for the X age groups.
## pabs: Vector of X. Proportion of abstainers for the X age groups.
## pformer: Vector of X. Proportion of former drinkers for the X age groups.
## pca: Single value or vector of X identical values. Per capita consumption of alcohol.
computeMeanStratum <- function(popsize, relcoef, pabs, pformer, pca, adjustPCA = 0.8)
{
  ## calculate the adjusted pca
  pca <- pca * 1000 * 0.789 * adjustPCA / 365   
  ## number of drinkers in each age group
  drk <- popsize * (1 - pabs - pformer)
  ## alcohol consumption over all age groups
  pcad <- pca * sum(popsize) / sum(drk)
  ## mean consumption per age group
  mu <- relcoef * pcad * sum(drk) / sum(relcoef * drk)
  return(mu)
}

############################################################################
#### the following function, computeMean, will compute the mean alcohol ####
#### consumption and standard deviation for all the inputs using        ####
#### the function computeMeanStratum                                    ####
############################################################################

## compute per gender - age category - region:
## mean consumption, SD, and the resulting estimates of the parameters
## for the gamma distribution (used to model alcohol consumption)
computeMean <- function(data, gender = NULL, age = NULL, adjustPCA = 0.8)
{
  ## handle gender subsetting
  if(!is.null(gender))
  {
    ## convert string input to numeric
    if(is.character(gender)) gender <- ifelse(gender == "male", 1, 2)
    ## take relevant gender subset,
    ## do not take age subset because mean of one age group is relative to the other groups, thus need all groups for calculation
    data <- subset(data, SEX %in% gender)
  }

  ## set up return argument
  ret <- subset(data, select = c(REGION, SEX, AGE_CATEGORY,
                        #PCA, VAR_PCA, POPULATION, RELATIVE_COEFFICIENTS, 
                        LIFETIME_ABSTAINERS, FORMER_DRINKERS))
  ret$MEAN <- NA
  
  
  ## generate index distinguishing all groups for which
  ## a separate mean should be calculated, namely all the
  ## possible region and sex combinations (labeled "pattern")
  index <- with(data, paste(REGION, SEX, sep = ":"))
  patterns <- unique(index)

  ## calculate mean for every region/sex/age combination
  for (i in patterns)
    {
      pat <- strsplit(i, ":")[[1]]
      dat <- data[data$REGION == pat[1] & data$SEX == pat[2],]
      ret$MEAN[ret$REGION == pat[1] & ret$SEX == pat[2]] <- computeMeanStratum(
      popsize = dat$POPULATION, relcoef = dat$RELATIVE_COEFFICIENT,
      pabs = dat$LIFETIME_ABSTAINERS, pformer = dat$FORMER_DRINKERS,
      pca = dat$PCA, adjustPCA = adjustPCA)
    }
 
  ## calculate sd (different factor for men and women)
  ret$SD <- 1.171 * ret$MEAN
  ret$SD[ret$SEX == 2] <- 1.258 * ret$MEAN[ret$SEX == 2]
   
  ## calculate estimates for the parameters of a gamma distribution
  ret$K <- ret$MEAN^2 / ret$SD^2 ## README: this is always (1/1.171)^2 or 1/1.258^2
  ret$THETA <- ret$SD^2 / ret$MEAN 
  
  ## take subset on age category
  if(!is.null(age)) ret <- subset(ret, AGE_CATEGORY %in% age)

  ## return calculated means
  return(ret)
}



#####################################################################
#### This function, calculates AAF, computes the AAFs and can use  ####
#### a defined number of cores (only available on linux so far)  ####
#####################################################################

## calculation of AAFs:
## input necessary:
## k and theta - from dataset returned of computeMean()
## prop of abstainers and former drinker per stratum  - from dataset returned of computeMean()
## region, sex, age category - from dataset returned of computeMean()
## number of diseases, relative risk function, coefficients, log relative risks for former drinkers - from list(s)-object
## README:
## attention: expects data structure like this:
## Oral_cavity = list(disease = name,
## RRCurrent = function(x, beta) {exp(x*beta)}, # relative risk function for the (current) drinkers
## betaCurrent = RRbeta$betaCurrent[RRbeta$Cancer_name==name]) # coefficient for the RR function of (current) drinkers

# calculate AAF for current drinkers
calculateAAF1 <- function(data, disease, mc.cores = 1, ...)
{
  ## for parallel computing
  require("parallel")
  
  ## define function calculating the AAF per disease
  ## this is to be applied to each element of the disease argument
  ## definion is done here to avoid having to pass on the data argument formally (lexical scoping)
  calAAFdisease <- function(dis){
    
    ## normalizing constant for gamma distribution
    normConst <- mapply(function(shape, scale) 
    {
      integrate(dgamma, lower = 0.1, upper = 150, shape = shape, scale = scale)$value
    },
    shape = data$K, scale = data$THETA
    )
    
    ## risk of drinkers
    drkfun2 <- function(x, pabs, pform, nc, k, theta)
    {
      (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta)  * (dis$RRCurrent(x, beta = dis$betaCurrent)-1)
    } # set shape of gamma distribution and risk at each alcohol level
    
    
    drk1 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2, lower = dis$lowerlim, upper = dis$midlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    drk2 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2, lower = dis$midlim, upper = dis$upperlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    
    ## AAF - calculate AAF per age group
    aaf <- (drk2) / (drk1 + drk2 + 1)      
    
    return(aaf)
    
  }
  
  ## run in parallel (on unix machines)
  aaf <- mclapply(disease, calAAFdisease, mc.cores = mc.cores)
  
  ## set up return object
  ret <- data.frame()
  for (i in 1:length(disease)) ret <- rbind(ret, data)
  ret$DISEASE <- rep(sapply(disease, function(x) x$disease), each = nrow(data))
  ret$AAF <- unlist(aaf)
  
  return(ret)
}


# calculate AAF for former drinkers only
calculateAAFformer <- function(data, disease, mc.cores = 1, ...)
{
  ## for parallel computing
  require("parallel")
  
  ## define function calculating the AAF per disease
  ## this is to be applied to each element of the disease argument
  ## definion is done here to avoid having to pass on the data argument formally (lexical scoping)
  calAAFdisease <- function(dis){
    
    ## normalizing constant for gamma distribution
    normConst <- mapply(function(shape, scale) 
    {
      integrate(dgamma, lower = 0.1, upper = 150, shape = shape, scale = scale)$value
    },
    shape = data$K, scale = data$THETA
    )
    
    ## risk of drinkers
    drkfun2 <- function(x, pabs, pform, nc, k, theta)
    {
      (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta)  * (dis$RRCurrent(x, beta = dis$betaCurrent)-1)
    } # set shape of gamma distribution and risk at each alcohol level
    
    
    drk1 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2, lower = dis$lowerlim, upper = dis$midlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    drk2 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2, lower = dis$midlim, upper = dis$upperlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    
    ## AAF - calculate AAF per age group
    aaf <- (data$FORMER_DRINKERS*(dis$lnRRFormer-1)) / (data$FORMER_DRINKERS*(dis$lnRRFormer-1) + drk1 + drk2 + 1)    
    
    return(aaf)
    
  }
  
  ## run in parallel (on unix machines)
  aaf <- mclapply(disease, calAAFdisease, mc.cores = mc.cores)
  
  ## set up return object
  ret <- data.frame()
  for (i in 1:length(disease)) ret <- rbind(ret, data)
  ret$DISEASE <- rep(sapply(disease, function(x) x$disease), each = nrow(data))
  ret$AAF <- unlist(aaf)
  
  return(ret)
}


# calculate AAF for current drinkers including former in denominator
calculateAAF2 <- function(data, disease, mc.cores = 1, ...)
{
  ## for parallel computing
  require("parallel")
  
  ## define function calculating the AAF per disease
  ## this is to be applied to each element of the disease argument
  ## definion is done here to avoid having to pass on the data argument formally (lexical scoping)
  calAAFdisease <- function(dis){
    
    ## normalizing constant for gamma distribution
    normConst <- mapply(function(shape, scale) 
    {
      integrate(dgamma, lower = 0.1, upper = 150, shape = shape, scale = scale)$value
    },
    shape = data$K, scale = data$THETA
    )
    
    ## risk of drinkers
    drkfun2 <- function(x, pabs, pform, nc, k, theta)
    {
      (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta)  * (dis$RRCurrent(x, beta = dis$betaCurrent)-1)
    } # set shape of gamma distribution and risk at each alcohol level
    
    
    drk1 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2, lower = dis$lowerlim, upper = dis$midlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    drk2 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2, lower = dis$midlim, upper = dis$upperlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    
    ## AAF - calculate AAF per age group
    aaf <- (drk2) / (data$FORMER_DRINKERS*(dis$lnRRFormer-1) + drk1 + drk2 + 1)    
    
    return(aaf)
    
  }
  
  ## run in parallel (on unix machines)
  aaf <- mclapply(disease, calAAFdisease, mc.cores = mc.cores)
  
  ## set up return object
  ret <- data.frame()
  for (i in 1:length(disease)) ret <- rbind(ret, data)
  ret$DISEASE <- rep(sapply(disease, function(x) x$disease), each = nrow(data))
  ret$AAF <- unlist(aaf)
  
  return(ret)
}


# calculate former + current drinkers together
calculateAAF3 <- function(data, disease, mc.cores = 1, ...)
{
  ## for parallel computing
  require("parallel")
  
  ## define function calculating the AAF per disease
  ## this is to be applied to each element of the disease argument
  ## definion is done here to avoid having to pass on the data argument formally (lexical scoping)
  calAAFdisease <- function(dis){
    
    ## normalizing constant for gamma distribution
    normConst <- mapply(function(shape, scale) 
    {
      integrate(dgamma, lower = 0.1, upper = 150, shape = shape, scale = scale)$value
    },
    shape = data$K, scale = data$THETA
    )
    
    ## risk of drinkers
    drkfun2 <- function(x, pabs, pform, nc, k, theta)
    {
      (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta)  * (dis$RRCurrent(x, beta = dis$betaCurrent)-1)
    } # set shape of gamma distribution and risk at each alcohol level
    
    
    drk1 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2, lower = dis$lowerlim, upper = dis$midlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    drk2 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2, lower = dis$midlim, upper = dis$upperlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    
    ## AAF - calculate AAF per age group
    aaf <- (data$FORMER_DRINKERS*(dis$lnRRFormer-1) + drk2) / (data$FORMER_DRINKERS*(dis$lnRRFormer-1) + drk1 + drk2 + 1)    
    
    return(aaf)
    
  }
  
  ## run in parallel (on unix machines)
  aaf <- mclapply(disease, calAAFdisease, mc.cores = mc.cores)
  
  ## set up return object
  ret <- data.frame()
  for (i in 1:length(disease)) ret <- rbind(ret, data)
  ret$DISEASE <- rep(sapply(disease, function(x) x$disease), each = nrow(data))
  ret$AAF <- unlist(aaf)
  
  return(ret)
}

# calculat AAF for current drinkers only by consumption level
calculateAAF_split <- function(data, disease, mc.cores = 1, ...)
{
  ## for parallel computing
  require("parallel")
  
  ## define function calculating the AAF per disease
  ## this is to be applied to each element of the disease argument
  ## definion is done here to avoid having to pass on the data argument formally (lexical scoping)
  calAAFdisease <- function(dis){
    
    ## normalizing constant for gamma distribution
    normConst <- mapply(function(shape, scale) 
    {
      integrate(dgamma, lower = 0.1, upper = 150, shape = shape, scale = scale)$value
    },
    shape = data$K, scale = data$THETA
    )
    
    ## risk of drinkers
    
    drkfun2b <- function(x, pabs, pform, nc, k, theta)
    {
      (1 - pabs - pform)/ nc * dgamma(x, shape = k, scale = theta) * (dis$RRCurrent(x, beta = dis$betaCurrent)-1)
    } # set shape of gamma distribution and risk at each alcohol level
    
    
    drk1 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2b, lower = dis$lowerlim, upper = dis$midlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    drk2 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2b, lower = dis$midlim, upper = dis$upperlim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    drk3 <- mapply(function(pabs, pform, nc, k, theta)
    {
      integrate(drkfun2b, lower = dis$upperlim, upper = dis$toplim, pabs = pabs, pform = pform, nc = nc,k = k, theta = theta)$value
    }, # compute population risk at each drinking level considering the proportion of drinkers
    pabs = data$LIFETIME_ABSTAINERS, pform = data$FORMER_DRINKERS, nc = normConst, k = data$K, theta = data$THETA
    ) # define abstainers, former drinkers, normalisation constant, k and theta variables
    
    ## AAF - calculate AAF per age group
    aaf <- drk2 / (drk1 + drk2 + drk3 + 1)    
    
    return(aaf)
    
  }
  
  ## run in parallel (on unix machines)
  aaf <- mclapply(disease, calAAFdisease, mc.cores = mc.cores)
  
  ## set up return object
  ret <- data.frame()
  for (i in 1:length(disease)) ret <- rbind(ret, data)
  ret$DISEASE <- rep(sapply(disease, function(x) x$disease), each = nrow(data))
  ret$AAF <- unlist(aaf)
  
  return(ret)
}

