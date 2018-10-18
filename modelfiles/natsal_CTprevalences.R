#By this function, the chlamydia prevalence among sexually actives is calculated (only those are included in the urine study)
#The chlamydia prevalence in the low risk category and in the high risk category are calculated, by
#combing NATSAL data on number of partners for an individual aged a, and the chlamydia status of that individual.
#individuals are assigned probabilistically in the low and high risk group.
#uncertainty is considered using bootstrapping
#input: fold is the parameter indicating the proportion of people in the low resp. high risk group

#infered prelance per age and activity group can be used in a Bayesian framework to better infer the parameters
#of the STI model, in particular the assortativity index epsilon. Epsilon can influence the prevalnce per risk group,
#in relation to the overall prevalence (Fig 17 winter school Althaus), so putting observations of the prevalence
#stratified per activity class as well as observations of the prevalence combined over the activity classes may
#lead to a better understanding of epsilon

get.natsal.CTprevalences <- function(nat,aggregatebyAC=T, n.boot = 1000, verbose=F){

  # use bootstrapping only if individual prevalences per AC have to be calculated
  if (aggregatebyAC){
    # dimensions: sex * age
    numerator <- array(0, c(2,n.age.classes))
    denominator <- array(0, c(2,n.age.classes))
    prev <- array(0, c(2,n.age.classes))
  } else {

    # dimensions: sex * activityclasses * age * nboot
    prev <- array(0, c(2,2,n.age.classes,n.boot))
    numerator <- array(0, c(2,2,n.age.classes,n.boot))
    denominator <- array(0, c(2,2,n.age.classes,n.boot))
  }
  
  
  ### loop through all genders
  for(g in 1:2){
    
    ### loop through all age classes
    for(i in 1:(n.age.classes)) {
      
      if(verbose){
        cat("looping through: gender ", sex.names[g], " and age class: ", age.classes[i], '\n')
      }
      
      ############################ EXPLANATIONS OF THE VARIABLE NAMES IN THE DATASET ####################################################
      ### hetnonew is the number of of new het. sex partners, last year (99 = missing, -1 = not applicable)
      ### dage is respondent's age at interview, in years 
      ### rsex is respondent's sex (1 = male, 2 = female)
      ### total_wt is the respondent's weight to produced a statistically well-balanced data set: (some respondent's are over or underrepresented, and this needs to be adjusted)
      ### urine_wt is the weight for the subset of urine samples
      ### ct_posconfirmed is chlamydia test results from urine sample (0: not detected; 1: detected)
      ###################################################################################################################################
      
      ### filter the data set, choose only a subset that does not contain misinformation
      
      ### respondent's age should be within the predefined age band
      ageband.cond <- nat$dage >= age.classes[i] & nat$dage < age.classes[i+1]
      
      ### respondent's gender should be male 
      gender.cond <- nat$rsex %in% g
      
      ct.cond <- nat$ct_posconfirmed=="negative" | nat$ct_posconfirmed=="positive"
      
      urine_wt.cond <- !is.na(nat$urine_wt)
      
      ### respondent should be sexually active
      eversex.cond <- nat$hetlife != 0
      
      if (aggregatebyAC){
        ### combine all conditions
        all.conds <- which(ageband.cond & gender.cond & ct.cond & urine_wt.cond & eversex.cond)
      } else {
        ### applicable number of new partners
        hetnonew.cond <- nat$hetnonew >= 0 & nat$hetnonew != 99 & nat$hetnonew < 995
        ### combine all conditions
        all.conds <- which(hetnonew.cond & ageband.cond & gender.cond & ct.cond & urine_wt.cond & eversex.cond)
      }
      
      ###################################################################################################################################
      ### SUBSET OF DATA
      ###################################################################################################################################
      
      ### extract subset that matches specified conditions
      data.CT <- nat$ct_posconfirmed[all.conds] 
      
      ### extract the weights from this subdata
      weight <- nat$urine_wt[all.conds]
      
      ### need to normalize the weights (because we are looking at a subset)
      weight <- weight/mean(weight)
      
      ### to be used for bootstrapping, only when aggregatebyAC=T
      if (!aggregatebyAC){
        data.hetnonew <- nat$hetnonew[all.conds] ### data.hetnonew contains a vector of the number of new heterosexual partners per year
        ndata.hetnonew <- length(data.hetnonew)
      }
      
      ###################################################################################################################################
      ### CALCULATIONS
      ###################################################################################################################################
      
      if (aggregatebyAC){ #bootstrap sample if individual prevalence per AC have to be obtained
        if (length(all.conds)>0){
          numerator[g,i] <- sum( weight[which(data.CT=="positive")] )
          denominator[g,i] <- sum(weight)
          prev[g,i] <- numerator[g,i] / denominator[g,i] #prevalence
        }
        
      } else {
        
        ### bootstraps
        for(k in 1:n.boot) {
          
          if (length(all.conds)>0){
            highrisk.id <- sample(1:ndata.hetnonew,size=floor(fold[g,2]*ndata.hetnonew),prob=data.hetnonew,replace=F)
            lowrisk.id <- setdiff(1:ndata.hetnonew, highrisk.id)
            
            if (length(weight[lowrisk.id])>0) {
              
              numerator[g,1,i,k] <- sum(weight[intersect(which(data.CT=="positive"),lowrisk.id)])
              denominator[g,1,i,k] <- sum(weight[lowrisk.id])
              prev[g,1,i,k] <- numerator[g,1,i,k] / denominator[g,1,i,k] #prevalence low risk
              
            }
            
            if (length(weight[highrisk.id])>0) {
              numerator[g,2,i,k] <- sum(weight[intersect(which(data.CT=="positive"),highrisk.id)])
              denominator[g,2,i,k] <- sum(weight[highrisk.id])
              prev[g,2,i,k] <- numerator[g,2,i,k]/ denominator[g,2,i,k] #prevalence low risk
            }
            
          }
          
        }
      }
    }
  }
  return(list(numerator=numerator,denominator=denominator,prev=prev))
}

#############################################################################################################
### write epidata for mcmc
#############################################################################################################

# #make natsal 2 epidata
# library(foreign)
# load("C:/Users/smid/Dropbox/Natsal/Natsal-2/natsal2.RData")
# natsal2 <- natsal2
# natsal2$ct_posconfirmed <- natsal2$c_result
# natsal2$ct_posconfirmed[natsal2$ct_posconfirmed==0] <- "negative"
# natsal2$ct_posconfirmed[natsal2$ct_posconfirmed==1] <- "positive"
# natsal2.prevs.all  <- get.natsal.CTprevalences(natsal2,aggregatebyAC=F)
# natsal2.meannumerator.prevs <- apply(natsal2.prevs.all$numerator,MARGIN=c(1,2,3),mean)
# natsal2.meandenominator.prevs <- apply(natsal2.prevs.all$denominator,MARGIN=c(1,2,3),mean)
# 
# #make natsal 3 epidata
# load("C:/Users/smid/Dropbox/Natsal/Natsal-3/natsal3.RData")
# nat2 <- read.dta("C:/Users/smid/Dropbox/Natsal/Natsal-3/NatSal3_new12.dta")
# natsal3$ct_posconfirmed <- nat2$ct_posconfirmed[match(natsal3$sin2, nat2$sin2)]
# natsal3$urine_wt <- nat2$urine_wt[match(natsal3$sin2, nat2$sin2)]
# natsal3.prevs.all  <- get.natsal.CTprevalences(natsal3,aggregatebyAC=F)
# natsal3.meannumerator.prevs <- apply(natsal3.prevs.all$numerator,MARGIN=c(1,2,3),mean)
# natsal3.meandenominator.prevs <- apply(natsal3.prevs.all$denominator,MARGIN=c(1,2,3),mean)

#By this function, the M. genitalium prevalence is calculated

get.natsal.MGprevalences <- function(nat,aggregatebyAC=T, n.boot = 1000, verbose=F){
  
  # use bootstrapping only if individual prevalences per AC have to be calculated
  if (aggregatebyAC){
    # dimensions: sex * age
    numerator <- array(0, c(2,n.age.classes))
    denominator <- array(0, c(2,n.age.classes))
    prev <- array(0, c(2,n.age.classes))
  } else {
    
    # dimensions: sex * activityclasses * age * nboot
    prev <- array(0, c(2,2,n.age.classes,n.boot))
    numerator <- array(0, c(2,2,n.age.classes,n.boot))
    denominator <- array(0, c(2,2,n.age.classes,n.boot))
  }
  
  
  ### loop through all genders
  for(g in 1:2){
    
    ### loop through all age classes
    for(i in 1:(n.age.classes)) {
      
      if(verbose){
        cat("looping through: gender ", sex.names[g], " and age class: ", age.classes[i], '\n')
      }
      
      ############################ EXPLANATIONS OF THE VARIABLE NAMES IN THE DATASET ####################################################
      ### hetnonew is the number of of new het. sex partners, last year (99 = missing, -1 = not applicable)
      ### dage is respondent's age at interview, in years 
      ### rsex is respondent's sex (1 = male, 2 = female)
      ### total_wt is the respondent's weight to produced a statistically well-balanced data set: (some respondent's are over or underrepresented, and this needs to be adjusted)
      ### urine_wt is the weight for the subset of urine samples
      ### ct_posconfirmed is chlamydia test results from urine sample (0: not detected; 1: detected)
      ###################################################################################################################################
      
      ### filter the data set, choose only a subset that does not contain misinformation
      
      ### respondent's age should be within the predefined age band
      ageband.cond <- nat$dage >= age.classes[i] & nat$dage < age.classes[i+1]
      
      ### respondent's gender should be male 
      gender.cond <- nat$rsex %in% g
      
      mg.cond <- nat$mg_posconfirmed=="negative" | nat$mg_posconfirmed=="positive"
      
      urine_wt.cond <- !is.na(nat$urine_wt)
      
      ### respondent should be sexually active
      eversex.cond <- nat$hetlife != 0
      
      if (aggregatebyAC){
        ### combine all conditions
        all.conds <- which(ageband.cond & gender.cond & mg.cond & urine_wt.cond & eversex.cond)
      } else {
        ### applicable number of new partners
        hetnonew.cond <- nat$hetnonew >= 0 & nat$hetnonew != 99 & nat$hetnonew < 995
        ### combine all conditions
        all.conds <- which(hetnonew.cond & ageband.cond & gender.cond & mg.cond & urine_wt.cond & eversex.cond)
      }
      
      ###################################################################################################################################
      ### SUBSET OF DATA
      ###################################################################################################################################
      
      ### extract subset that matches specified conditions
      data.CT <- nat$mg_posconfirmed[all.conds] 
      
      ### extract the weights from this subdata
      weight <- nat$urine_wt[all.conds]
      
      ### need to normalize the weights (because we are looking at a subset)
      weight <- weight/mean(weight)
      
      ### to be used for bootstrapping, only when aggregatebyAC=T
      if (!aggregatebyAC){
        data.hetnonew <- nat$hetnonew[all.conds] ### data.hetnonew contains a vector of the number of new heterosexual partners per year
        ndata.hetnonew <- length(data.hetnonew)
      }
      
      ###################################################################################################################################
      ### CALCULATIONS
      ###################################################################################################################################
      
      if (aggregatebyAC){ #bootstrap sample if individual prevalence per AC have to be obtained
        if (length(all.conds)>0){
          numerator[g,i] <- sum( weight[which(data.CT=="positive")] )
          denominator[g,i] <- sum(weight)
          prev[g,i] <- numerator[g,i] / denominator[g,i] #prevalence
        }
        
      } else {
        
        ### bootstraps
        for(k in 1:n.boot) {
          
          if (length(all.conds)>0){
            highrisk.id <- sample(1:ndata.hetnonew,size=floor(fold[g,2]*ndata.hetnonew),prob=data.hetnonew,replace=F)
            lowrisk.id <- setdiff(1:ndata.hetnonew, highrisk.id)
            
            if (length(weight[lowrisk.id])>0) {
              
              numerator[g,1,i,k] <- sum(weight[intersect(which(data.CT=="positive"),lowrisk.id)])
              denominator[g,1,i,k] <- sum(weight[lowrisk.id])
              prev[g,1,i,k] <- numerator[g,1,i,k] / denominator[g,1,i,k] #prevalence low risk
              
            }
            
            if (length(weight[highrisk.id])>0) {
              numerator[g,2,i,k] <- sum(weight[intersect(which(data.CT=="positive"),highrisk.id)])
              denominator[g,2,i,k] <- sum(weight[highrisk.id])
              prev[g,2,i,k] <- numerator[g,2,i,k]/ denominator[g,2,i,k] #prevalence low risk
            }
            
          }
          
        }
      }
    }
  }
  return(list(numerator=numerator,denominator=denominator,prev=prev))
}

