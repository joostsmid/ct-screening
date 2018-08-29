# Here we set up a series of procedures to calculate the rate at which the respondents 
# males (or females), a particular age a, engage in sexual partnerships with members of the 
# opposite sex of age a' for the first time in their life. This is the rate of the 
# beginning of sexual activity. Also ages younget than 16 are allowed but computed differently
# from ages>=16. For Ages >= 16 the fractions of virgins are lost by considering each age group
# separately, and dteremining the fraction of people that age that have never had sex.
# For ages younger than 16, we determine ages at which people had first sex. This number is more
# prone to historical bais becuase data reportted about longer ago, by older people, is considered
# instead of fractions in 2011. To minimize this bias: only data of respondents max 25y of age is used.
# method used for smoothing: monotonic b-spline using smooth.monotone function
###############################################################################

### load libraries 
library(spatial)
library(sn)
library(RColorBrewer)
library(weights)
library(fda)

### function that generates probability for still being virgin as function of age for females and males
prob.virgin <- function(nat, age.classes, use.weights = TRUE, verbose = F){
  
  sex.names <- c("M","F") # gender
  sex.code <-c(1,2) #coding used in NATSAL for gender
  
  n.age.classes <- length(age.classes)-1
  age.class.groups <- age.classes[1:n.age.classes]
  
  v.prob <- to.tensor(NA, dims = list(sex=c("M","F") , age=paste0('age', age.class.groups)) )
  
  ############################ EXPLANATIONS OF THE VARIABLE NAMES IN THE DATASET ####################################################
  ### afsex is age at 1st het. sex. intercourse if aged 13+ (code 96 if sexual intercourse has not happened yet)
  ### dage is respondent's age at interview, in years 
  ### rsex is respondent's sex (1 = male, 2 = female)
  ### total_wt is the respondent's weight to produced a statistically well-balanced data set: (some respondent's are over or underrepresented, and this needs to be adjusted)
  ###################################################################################################################################
  
  for(g in sex.code){
    
    ### gender condition
    gender.cond <- nat$rsex %in% sex.code[g]
    
    ################################################################################################################################
    ### calculations for ages >= 16
    ################################################################################################################################
    
    ### loop through all respondent's ages above 16
    youngest_natsal <- which(age.classes>=16)[1]
    for(i in youngest_natsal:n.age.classes) {
      
      from <- age.classes[i]
      to <- age.classes[i+1]
      
      if(verbose){
        cat("looping through: gender ", sex.names[g], " and the age group from ", from, " to", to, '\n')
      }
      
      ### respondent's age should be within the predefined age band
      ageband.cond <- nat$dage >= from & nat$dage< to
      
      stillvirgin.cond <- nat$afsex == 96 | nat$afsex == 90 #(abuse (code 90) is counted as virginity)
      
      nomorevirgin.cond <- nat$afsex < 90 #no answers (code 97 and 99) are not taken into account
      
      stillvirgin.ind <- which(gender.cond & ageband.cond & stillvirgin.cond)
      nomorevirgin.ind <- which(gender.cond & ageband.cond & nomorevirgin.cond)
      
      #compute (weighted) proportion of virgins in that age band
      
      if(use.weights){
        weightedno.virgins <- sum(nat$total_wt[stillvirgin.ind])
        weightedtotal <- weightedno.virgins + sum(nat$total_wt[nomorevirgin.ind])
        v.prob[g,i] <- weightedno.virgins / weightedtotal # nat$total_wt is the respondent's weight to produce a statistically well-balanced data set
      }
      else{
        v.prob[g,i] <-  length(stillvirgin.ind) / sum(length(stillvirgin.ind),length(nomorevirgin.ind))
      }
    }
    
    ################################################################################################################################
    ### calculations for ages below 16
    ################################################################################################################################
    
    if (youngest_natsal>1){ # if there is a lower limit age band below age 16
      
      # new JS130217:
      # first compute proportion virgin of age 16 - age.classes[youngest_natsal]
      
      if (age.classes[youngest_natsal]>16){
        
        from <- 16
        to <- age.classes[youngest_natsal]
        
        ### respondent's age should be within the predefined age band
        ageband.cond <- nat$dage >= from & nat$dage< to
        
        stillvirgin.cond <- nat$afsex == 96 | nat$afsex == 90 #(abuse (code 90) is counted as virginity)
        
        nomorevirgin.cond <- nat$afsex < 90 #no answers (code 97 and 99) are not taken into account
        
        stillvirgin.ind <- which(gender.cond & ageband.cond & stillvirgin.cond)
        nomorevirgin.ind <- which(gender.cond & ageband.cond & nomorevirgin.cond)
        
        #compute (weighted) proportion of virgins in that age band
        
        if(use.weights){
          weightedno.virgins <- sum(nat$total_wt[stillvirgin.ind])
          weightedtotal <- weightedno.virgins + sum(nat$total_wt[nomorevirgin.ind])
          v.prob_16 <- weightedno.virgins / weightedtotal # nat$total_wt is the respondent's weight to produce a statistically well-balanced data set
        }
        else{
          # proportion virgins in [16, youngestnatsal], a subinterval of age.classes
          v.prob_16 <-  length(stillvirgin.ind) / sum(length(stillvirgin.ind),length(nomorevirgin.ind)) 
        }
        
      } else {
        v.prob_16<-as.numeric(v.prob[g,youngest_natsal])
      }
      
      ### Distribute all ages at which participants (with dage<=25) lost virginity before the age of age.classes[youngest_natsal]
      ### participants data are restricted to participant younger than 25, to minimize historical bias
      available.cond <- nat$afsex < 16 #= age.classes[youngest_natsal]
      age.cond <- nat$dage<=25 #to avoid too much historical bias: only include data of participants aged 16-25
      agesfirstsex <- nat$afsex[gender.cond & available.cond & age.cond] #vector of ages of first sex for participants <= 25y, who had first sex under 16
      #aggregate ages in vector agesfirstsex: replace by ages in lower bound in (aggregated) age classes
      
      if(use.weights){
        weights <- nat$total_wt[gender.cond & available.cond & age.cond]
        ### need to normalize the weights (because we are looking at a subset)
        weights <- weights/sum(weights)
        # weighted frequency table of ages first sex (expressed as percentages) of those that have lost virginity before 16
        pct.agesfirstsex <- wpct(agesfirstsex, weights)
      } else{
        # unweighted frequency table of ages first sex (expressed as percentages) of those that have lost virginity before 16
        pct.agesfirstsex <- table(agesfirstsex)/sum(table(agesfirstsex))
      }
      
      # P(virgin loss in ages 13-15) = P(virgin loss in ages 13-15 | no more virgin at age 16) * P(no more virgin at age 16) 
      #                              = P(virgin loss in ages 13-15 | no more virgin at age 16) * (1 - v.prob_16)
      # we scale the probabilities per age before 16 of having had first sex to match the probability of being no longer virgin before age 16
      prob.virgin.before16 <- v.prob_16 + (1-pct.agesfirstsex)*(1-v.prob_16)
      
      for(i in 1:(youngest_natsal-1)) {
        
        from <- age.classes[i]
        to <- age.classes[i+1]
        
        if (to!=age.classes[youngest_natsal]){
          #take mean over ages in that interval to account for the ages at which most people loose virginity within an age interval
          v.prob[g,i] <- mean(prob.virgin.before16[names(prob.virgin.before16)>= from & names(prob.virgin.before16)< to]) 
        } else {
          v.prob[g,i] <- mean(c(prob.virgin.before16[names(prob.virgin.before16)>= from & names(prob.virgin.before16)< to],
                                rep(v.prob_16,length(16:age.classes[youngest_natsal])-1)))
        }
        
      }
      
    }
    
  }
  
  return(v.prob)
  
}

prob.virgin.smooth <- function(v.prob){
  
  sex.names <- c("M","F") # gender
  sex.code <-c(1,2) #coding used in NATSAL for gender
  
  n.age.classes <- length(age.classes)-1
  age.class.groups <- age.classes[1:n.age.classes]
  
  ################################################################################################################################
  ### smoothing
  ################################################################################################################################
    
    v.prob.smooth <- to.tensor(NA, dims = list(sex=c("M","F") , age=paste0('age', age.class.groups)) )
    
    for(g in sex.code){
      
      #  Range of observations
      rng <- c(age.class.groups[1],tail(age.class.groups,n=1))
      #  First set up a basis for monotone smooth
      #  We use b-spline basis functions of order 6
      #  Knots are positioned at the ages of observation.
      norder <- 6
      nage   <- n.age.classes
      nbasis <- n.age.classes + norder - 2
      wbasis <- create.bspline.basis(rng, nbasis, norder, age.class.groups)
      #  starting values for coefficient
      cvec0 <- matrix(0,nbasis,1)
      Wfd0  <- fd(cvec0, wbasis)
      #  set up functional parameter object
      Lfdobj    <- 3          #  penalize curvature of acceleration
      lambda    <- 10^(-0.5)  #  smoothing parameter
      growfdPar <- fdPar(Wfd0, Lfdobj, lambda)
      
      result <- smooth.monotone(age.class.groups, as.numeric(v.prob[g,]), growfdPar)
      #  Extract the functional data object and regression coefficients
      Wfd  <- result$Wfdobj
      beta <- result$beta
      
      x.smooth <- beta[1] + beta[2]*eval.monfd(age.class.groups, Wfd)
      
      v.prob.smooth[g,] <- pmin(as.numeric(v.prob[g,1]), x.smooth) # virginity cannot be above 100%
      
    }

  return(v.prob.smooth)
  
}

