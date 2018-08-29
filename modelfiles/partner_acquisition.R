# Here we set up a series of procedures to calculate the average number of times in which 
# males (or females), of a particular age a an a particular activity class j, 
# engage in sexual contacts with members of the opposite sex. To do this, we look 
# at the distribution of sex acts per year in cohorts of the same age, and 
# fit the sum of nJ poisson distributions to this empirical distribution, obtaining 
# parameter estimates for each. 
# The resulting tensors are the used to calculate the 
# parameter ct in the genericSTImod.R script. 
#
# Authors: Victor Garcia, Joost Smid and Christian Althaus, University of Bern, September 2015
###############################################################################

# The function new.contact.nmbs.nat2_3combined() combines partner acquisition data from natsal 2 and 3
new.contact.nmbs.nat2_3combined <- function(nJ, age.classes){
  
  #produce warning when first age class smaller than 13: then inconsistencies occur in calculation acquisition rates becuase proportion of lifetime lived in
  # sexual activity before and including 13 should only include age 13 and not ages before 13. Disproportionally high number of sex acts attributed to that age group
  # in line "prop.acts.beforeyoungestnatsal <- prop.nomorevirg.beforeyoungestnatsal * diff(age.classes)[1:(youngest_natsal-1)]"
  if(age.classes[1]<13) cat("Warning: First age class smaller than 13: inconsistencies occur in calculation acquisition rates for ages smaller than 13.","\n")
  
  sex.names <- c("M","F") # gender
  sex.code <-c(1, 2) #coding used in NATSAL for gender
  n.age.classes <- length(age.classes)-1
  dimnames1 <- list(sex=c("M","F"), j=paste0('j',1:nJ), age=paste0('age', age.classes[1:(n.age.classes)]))  # JS090316
  nElements <-  nJ*(n.age.classes)*length(sex.names) #JS130116
  
  ### define output tensor
  ct.mean <- to.tensor(rep(0,nElements), dims = dimnames1)
  
  #create matrix of proportions of people in risk categories, all but the highest
  #The proportion of people in the highest risk category is 1-this number
  #different rows refer to different genders
  dimnames2 <- list(sex=c("M","F"), j=paste0('j',1:nJ) )  # JS120416
  nElements2 <-  nJ*length(sex.names) #JS120416
  risk.cat <- to.tensor(rep(0,nElements2), dims = dimnames2) #JS120416
  
  youngest_natsal <- which(age.classes >= 16)[1] #youngest age in age.classes from which information is available from natsal
  
  ########################################################################################################################
  ### Calculate the frequencies of individuals within activity classes across all age groups
  ########################################################################################################################
  
  ### loop through all genders
  for(g in 1:length(sex.names)){
    
    sex <- sex.code[g]
    
    ############################ EXPLANATIONS OF THE VARIABLE NAMES IN THE DATASET ####################################################
    ### hetnonew is the number of of new het. sex partners, last year (99 = missing, -1 = not applicable)
    ### dage is respondent's age at interview, in years 
    ### rsex is respondent's sex (1 = male, 2 = female)
    ### total_wt is the respondent's weight to produced a statistically well-balanced data set: (some respondent's are over or underrepresented, and this needs to be adjusted)
    ###################################################################################################################################
    
    ### filter the data set, choose only a subset that does not contain misinformation
    ### applicable number of new partners (already conditioned on ever having had sex, otherwise this variable would take value -1)
    hetnonew.cond <- nat2$hetnonew >= 0 & nat2$hetnonew != 99 & nat2$hetnonew < 995
    
    ### respondent's age should be within the predefined age band
    ageband.cond.allages <- nat2$dage >= age.classes[youngest_natsal] & nat2$dage <= age.classes[n.age.classes]
    
    ### respondent's gender should be male 
    gender.cond <- nat2$rsex %in% sex
    
    ### respondent should be sexually active (new 130816)
    sex.active.cond <- nat2$afsex != 90 & nat2$afsex != 96 & nat2$afsex != 97 & nat2$afsex != 99
    
    ### combine all conditions
    all.conds.allages <- which(hetnonew.cond & ageband.cond.allages & gender.cond & sex.active.cond)
    
    ### extract subset that matches specified conditions
    data.allages2 <- nat2$hetnonew[all.conds.allages] ### data contains a vector of the number of new heterosexual partners per year
    
    ### extract the weights from this subdata
    weight.allages2 <- nat2$total_wt[all.conds.allages]
    
    ### Natsal-3 ###
    
    ### filter the data set, choose only a subset that does not contain misinformation
    ### applicable number of new partners (already conditioned on ever having had sex, otherwise this variable would take value -1)
    hetnonew.cond <- nat3$hetnonew >= 0 & nat3$hetnonew != 99 & nat3$hetnonew < 995
    
    ### respondent's age should be within the predefined age band
    ageband.cond.allages <- nat3$dage >= age.classes[youngest_natsal] & nat3$dage <= age.classes[n.age.classes]
    
    ### respondent's gender should be male 
    gender.cond <- nat3$rsex %in% sex
    
    ### respondent should be sexually active (new 130816)
    sex.active.cond <- nat3$afsex != 90 & nat3$afsex != 96 & nat3$afsex != 97 & nat3$afsex != 99
    
    ### combine all conditions
    all.conds.allages <- which(hetnonew.cond & ageband.cond.allages & gender.cond & sex.active.cond)
    
    ### extract subset that matches specified conditions
    data.allages3 <- nat3$hetnonew[all.conds.allages] ### data contains a vector of the number of new heterosexual partners per year
    
    ### extract the weights from this subdata
    weight.allages3 <- nat3$total_wt[all.conds.allages]
    
    ### Combine ###
    data.allages <- c(data.allages2,data.allages3)
    weight.allages <- c(weight.allages2,weight.allages3)
    
    ### need to normalize the weights (because we are looking at a subset)
    weight.allages <- weight.allages/mean(weight.allages)
    
    ###################################################################################################################################
    ### Maximum Likelihood estimation: proportions in each risk category
    ###################################################################################################################################
    
    if(nJ == 1){
      
      risk.cat[g] <- 1
      
    }
    
    if(nJ == 2){
      
      f.allages <- function(x,a,m1,m2){
        # a: proportion in high risk category
        # m1: expected number of partners in high risk category
        # m1: expected number of partners in low risk category
        return(a*dpois(x,m1)+(1-a)*dpois(x,m2))
      }
      
      ### the log likelihood function, constrain to mathematically reasonable parameter combinations
      nll.allages <- function(a,m1,m2) {
        
        if( a < 0 | a > 1 | m1 < 0 | m2 < 0 ){
          res <- NA
        }
        
        else{
          res <- -sum(weight.allages*log(f.allages(data.allages,a,m1,m2)))
        }
        
        return(res)
      }
      
      ### run estimation process with the stats4 package, mle function
      try(est3.allages <- bbmle::mle2(nll.allages,start=list(a=0.9,m1=0.5,m2=5), method = "SANN"), silent = TRUE)
      
      ### extract the fractions of the individuals in the age classes, and keep them fixed
      fractions <- est3.allages@coef[1]
      
      risk.cat[g,1] <- fractions #JS120416
      risk.cat[g,2] <- 1-fractions #JS120416
      
    }
    
    if(nJ == 3){
      
      f.allages <- function(x,a,b,m1,m2,m3){
        
        return(a*dpois(x,m1)+b*dpois(x,m2)+(1-a-b)*dpois(x,m3))
      }
      
      ### the log likelihood function
      nll.allages <- function(a,b,m1,m2,m3) {
        
        ### constrain to mathematically reasonable parameter combinations
        if( a < 0 | a > 1 |  b < 0 | b > 1 | (a + b ) > 1 |  m1 < 0 | m2 < 0 | m3 < 0 ){
          res <- NA
        }
        
        else{
          res <- -sum(weight.allages*log(f.allages(data.allages,a,b,m1,m2,m3)))
        }
        
        return(res)
      }
      
      ### run estimation process with the stats4 package, mle function
      try(est3.allages <- bbmle::mle2(nll.allages,start=list(a=0.9,b=0.05,m1=0.5,m2=5,m3=10), method = "SANN"), silent = TRUE )
      
      ### extract the fractions of the individuals in the age classes, and keep them fixed
      fractions <- est3.allages@coef[1:2]
      
      #risk.cat[g,] <- fractions #JS120416
      
      risk.cat[g,1:2] <- fractions #JS120416
      risk.cat[g,3] <- 1-sum(fractions) #JS120416
      
    }
    
    ########################################################################################################################
    ### Proceed to calculate estimates per age class, for ages >= 16
    ########################################################################################################################
    
    ### loop through all respondent's ages above 16
    for(i in youngest_natsal:(n.age.classes)) {
      
      # sex <- sex.code[g]  
      from <- age.classes[i]
      to <- age.classes[i+1]
      
      
      ### Natsal-2 ###
      
      ### applicable number of new partners (already conditioned on ever having had sex, otherwise this variable would take value -1)
      hetnonew.cond <- nat2$hetnonew >= 0 & nat2$hetnonew != 99 & nat2$hetnonew < 995
      
      ### respondent's age should be within the predefined age band
      ageband.cond <- nat2$dage >= from & nat2$dage < to
      
      ### respondent's gender should be male 
      gender.cond <- nat2$rsex %in% sex
      
      ### combine all conditions
      all.conds <- which(hetnonew.cond & ageband.cond & gender.cond)
      
      ### extract subset that matches specified conditions
      data2 <- nat2$hetnonew[all.conds] ### data contains a vector of the number of new heterosexual partners per year
      
      ### extract the weights from this subdata
      weight2 <- nat2$total_wt[all.conds]
      
      ### Natsal-3 ###
      
      ### applicable number of new partners (already conditioned on ever having had sex, otherwise this variable would take value -1)
      hetnonew.cond <- nat3$hetnonew >= 0 & nat3$hetnonew != 99 & nat3$hetnonew < 995
      
      ### respondent's age should be within the predefined age band
      ageband.cond <- nat3$dage >= from & nat3$dage < to
      
      ### respondent's gender should be male 
      gender.cond <- nat3$rsex %in% sex
      
      ### combine all conditions
      all.conds <- which(hetnonew.cond & ageband.cond & gender.cond)
      
      ### extract subset that matches specified conditions
      data3 <- nat3$hetnonew[all.conds] ### data contains a vector of the number of new heterosexual partners per year
      
      ### extract the weights from this subdata
      weight3 <- nat3$total_wt[all.conds]
      
      ### Combine ###
      
      data <- c(data2,data3)
      weight <- c(weight2,weight3)
      
      ### need to normalize the weights (because we are looking at a subset)
      weight <- weight/mean(weight)
      
      
      if(length(data) == 0){
        next;
        cat("no data to proceed with calculation; jumping to next case", '\n')
      }
      
      ###################################################################################################################################
      ### Maximum Likelihood estimation of the number of new partners per year for each activity class
      ###################################################################################################################################
      
      if(nJ == 1){
        
        ### combination of two poisson distributions
        f <- function(x,m1){
          
          return(dpois(x,m1))
        }
        
        ### the log likelihood function
        nll <- function(m1) {
          
          ### constrain to mathematically reasonable parameter combinations
          if(  m1 < 0){
            res <- NA
          }
          
          else{
            res <- -sum(weight*log(f(data,m1)))
          }
          return(res)
        }
        
        
        ### run estimation process with the stats4 package, mle function
        try(est3 <- bbmle::mle2(nll,start=list(m1=0.5), method = "SANN"), silent = TRUE )
        
      }
      
      if(nJ == 2){
        
        ### combination of two poisson distributions
        f <- function(x,m1,m2){
          
          return(fractions[1]*dpois(x,m1) + (1-fractions[1])*dpois(x,m2))
        }
        
        ### the log likelihood function
        nll <- function(m1,m2) {
          
          ### constrain to mathematically reasonable parameter combinations
          if(  m1 < 0 | m2 < 0 ){
            res <- NA
          }
          
          else{
            res <- -sum(weight*log(f(data,m1,m2)))
          }
          
          return(res)
        }
        
        
        ### RUN ESTIMATION PROCESS 
        ### run estimation process with nlm function
        try(est3 <- bbmle::mle2(nll,start=list(m1=0.5,m2=5), method = "SANN"), silent = TRUE )
        
      }
      
      ##### MLE of THREE poisson distributions #####
      if(nJ == 3){
        ### combination of three poisson distributions
        f <- function(x,m1,m2,m3){
          
          return(fractions[1]*dpois(x,m1)+fractions[2]*dpois(x,m2)+(1-fractions[1]-fractions[2])*dpois(x,m3))
        }
        
        ### the log likelihood function
        nll <- function(m1,m2,m3) {
          
          ### constrain to mathematically reasonable parameter combinations
          if(  m1 < 0 | m2 < 0 | m3 < 0 ){
            res <- NA
          }
          
          else{
            res <- -sum(weight*log(f(data,m1,m2,m3)))
          }
          return(res)
        }
        
        ### RUN ESTIMATION PROCESS 
        ### run estimation process with the stats4 package, mle function
        try(est3 <- bbmle::mle2(nll,start=list(m1=0.5,m2=5,m3=10), method = "SANN"), silent = TRUE )
        
      }
      
      ### check whether ther has been an error: if yes, enter an NA value
      if(class(est3) == "try-error" | is.null(est3) ){
        
        for(j in 1:nJ){
          ct.mean[sex.names[g],j, paste0('age',age.classes[i])] <- NA #JS220716
        }
      }
      else{ ### if no, enter the values of the estimates
        
        ### filling up the entries of the tensor in the j-th dimension
        for(j in 1:nJ){
          ct.mean[sex.names[g],j, paste0('age',age.classes[i])] <- est3@coef[j]
        }
        
      }
      
      ### reset estimation 
      est3 <- NULL
      
    }
    
    ########################################################################################################################
    ### Proceed to calculate estimates per age class, for ages < 16 (NEW 220816)
    ########################################################################################################################
    
    if (youngest_natsal>1){
      
      ### Natsal-2 ###
      
      ### applicable number of new partners (already conditioned on ever having had sex, otherwise this variable would take value -1)
      hetnonew.cond <- nat2$hetnonew >= 0 & nat2$hetnonew != 99 & nat2$hetnonew < 995
      
      ### applicable number of lifetime partners
      hetlife.cond <- nat2$hetlife >= 0 & nat2$hetlife != 99 & nat2$hetlife < 995
      
      age16.cond <- nat2$dage==16
      
      ### respondent's gender should be male 
      gender.cond <- nat2$rsex %in% sex
      
      ### respondent should be sexually active (new 130816)
      sex.active.cond <- nat2$afsex != 90 & nat2$afsex != 96 & nat2$afsex != 97 & nat2$afsex != 99
      
      all16.cond <- which(hetnonew.cond & hetlife.cond & age16.cond & gender.cond & sex.active.cond)
      
      partnersbefore16_2 <- nat2$hetlife[all16.cond] - nat2$hetnonew[all16.cond] # vector of number of partners before 16
      
      ### extract the weights from this subdata
      weight2 <- nat2$total_wt[all16.cond]
      
      
      
      ### Natsal-3 ###
      
      ### applicable number of new partners (already conditioned on ever having had sex, otherwise this variable would take value -1)
      hetnonew.cond <- nat3$hetnonew >= 0 & nat3$hetnonew != 99 & nat3$hetnonew < 995
      
      ### applicable number of lifetime partners
      hetlife.cond <- nat3$hetlife >= 0 & nat3$hetlife != 99 & nat3$hetlife < 995
      
      age16.cond <- nat3$dage==16
      
      ### respondent's gender should be male 
      gender.cond <- nat3$rsex %in% sex
      
      ### respondent should be sexually active (new 130816)
      sex.active.cond <- nat3$afsex != 90 & nat3$afsex != 96 & nat3$afsex != 97 & nat3$afsex != 99
      
      all16.cond <- which(hetnonew.cond & hetlife.cond & age16.cond & gender.cond & sex.active.cond)
      
      partnersbefore16_3 <- nat3$hetlife[all16.cond] - nat3$hetnonew[all16.cond] # vector of number of partners before 16
      
      ### extract the weights from this subdata
      weight3 <- nat3$total_wt[all16.cond]
      
      ### Combine ###
      
      weight <- c(weight2,weight3)
      partnersbefore16 <- c(partnersbefore16_2, partnersbefore16_3)
      
      ### need to normalize the weights (because we are looking at a subset)
      weight <- weight/mean(weight)
      
      ###################################################################################################################################
      ### Maximum Likelihood estimation of total number of new partners for ages below 16
      ###################################################################################################################################
      
      if(nJ == 1){
        
        ### combination of two poisson distributions
        f <- function(x,m1){
          
          return(dpois(x,m1))
        }
        
        ### the log likelihood function
        nll <- function(m1) {
          
          ### constrain to mathematically reasonable parameter combinations
          if(  m1 < 0){
            res <- NA
          }
          
          else{
            res <- -sum(weight*log(f(partnersbefore16,m1)))
          }
          return(res)
        }
        
        
        ### run estimation process with the stats4 package, mle function
        try(mle.est.below16 <- bbmle::mle2(nll,start=list(m1=0.5), method = "SANN"), silent = TRUE )
        
      }
      
      if(nJ == 2){
        
        ### combination of two poisson distributions
        f <- function(x,m1,m2){
          
          return(fractions[1]*dpois(x,m1) + (1-fractions[1])*dpois(x,m2))
        }
        
        ### the log likelihood function
        nll <- function(m1,m2) {
          
          ### constrain to mathematically reasonable parameter combinations
          if(  m1 < 0 | m2 < 0 ){
            res <- NA
          }
          
          else{
            res <- -sum(weight*log(f(partnersbefore16,m1,m2)))
          }
          
          return(res)
        }
        
        
        ### RUN ESTIMATION PROCESS 
        ### run estimation process with nlm function
        try(mle.est.below16 <- bbmle::mle2(nll,start=list(m1=0.5,m2=5), method = "SANN"), silent = TRUE )
        
      }
      
      ##### MLE of THREE poisson distributions #####
      if(nJ == 3){
        ### combination of three poisson distributions
        f <- function(x,m1,m2,m3){
          
          return(fractions[1]*dpois(x,m1)+fractions[2]*dpois(x,m2)+(1-fractions[1]-fractions[2])*dpois(x,m3))
        }
        
        ### the log likelihood function
        nll <- function(m1,m2,m3) {
          
          ### constrain to mathematically reasonable parameter combinations
          if(  m1 < 0 | m2 < 0 | m3 < 0 ){
            res <- NA
          }
          
          else{
            res <- -sum(weight*log(f(partnersbefore16,m1,m2,m3)))
          }
          return(res)
        }
        
        ### RUN ESTIMATION PROCESS 
        ### run estimation process with the stats4 package, mle function
        try(mle.est.below16 <- bbmle::mle2(nll,start=list(m1=0.5,m2=5,m3=10), method = "SANN"), silent = TRUE )
        
      }
      
      
      #proportion of repsondents that is not virgin anymore, before the youngest age in age.classes from which information is available from natsal (youngest_natsal)
      prop.nomorevirg.beforeyoungestnatsal <- as.numeric(1-vl[g,1:(youngest_natsal-1)]) #JS220816
      # proportion of the sex acts that is done when being in the different age groups before youngest_natsal (scaled with # years of being in that age band)
      prop.acts.beforeyoungestnatsal <- prop.nomorevirg.beforeyoungestnatsal * diff(age.classes)[1:(youngest_natsal-1)]
      prop.acts.beforeyoungestnatsal <- prop.acts.beforeyoungestnatsal/sum(prop.acts.beforeyoungestnatsal) #normalize
      
      ### check whether ther has been an error: if yes, enter an NA value
      if(class(mle.est.below16) == "try-error" | is.null(mle.est.below16) ){
        
        for(j in 1:nJ){
          ct.mean[sex.names[g],j, paste0('age',age.classes[1:(youngest_natsal-1)])] <- NA #JS220716
        }
      } else{ ### if no, enter the values of the estimates
        
        #MEAN VALUES
        ### filling up the entries of the tensor in the j-th dimension
        for(j in 1:nJ){
          # by multiplying the total number of sex acts before 16 (mle estimate) with the proportion of the sex acts before 16,
          # we get the number of sex acts in the different age classes before 16
          ct.mean[sex.names[g],j, paste0('age',age.classes[1:(youngest_natsal-1)])] <- mle.est.below16@coef[j] * prop.acts.beforeyoungestnatsal
        }
        
      }
      
      ### reset estimation 
      mle.est.below16 <- NULL
      
    }
    
  }

  return(list(ct.mean=ct.mean, risk.cat=risk.cat))

}

##### We replace non existing entries by replacing them by averages of the adjacent values. 
replace.na.by.averaging <- function(tensor.object){
  
  ### extract the information about the tensor object
  genders <- gsub("(.*)", '\\1' ,dimnames(tensor.object)$sex)
  activity.classes <- as.numeric(gsub("j(.*)", '\\1' ,dimnames(tensor.object)$j))
  ages <- as.numeric(gsub("age(.*)", '\\1' ,dimnames(tensor.object)$age))

  ### loop through genders
  for(g in 1:length(genders)){
    
    ### loop through activity classes
    for(j in 1:length(activity.classes)){
      
       ### age values 
      x <- ages
      
      ### entry values for spline
      ### the response variable are the number of contacts
      y <- tensor.object[g, j, ]
      
      na.ind <- which(is.na(y))
      
      ### proceed to approximate the missing values by adjacent values

      ### this only works if the na's are very sparse!!!!!!!!!!!!
      if(all(na.ind != 1) & all(na.ind != length(y)) ){
        y[na.ind] <- sapply(1:length(na.ind), function(w) mean(c(y[na.ind[w]-1],y[na.ind[w]+1])) )
      }
      else{
        y <- replace(y, na.ind, 0)
      }
      
      tensor.object[g, j, ] <- y
      
    }
  }
  
  return(tensor.object)
}


smooth.est.age <- function(tensor.object, df){
  
  ### extract the information about the tensor object
  genders <- gsub("(.*)", '\\1' ,dimnames(tensor.object)$sex)
  activity.classes <- as.numeric(gsub("j(.*)", '\\1' ,dimnames(tensor.object)$j))
  ages <- as.numeric(gsub("age(.*)", '\\1' ,dimnames(tensor.object)$age))
  #youngest_natsal <- which(ages >= 16)[1] #youngest age in age.classes from which information is available from natsal
  
  ### loop through genders
  for(g in 1:length(genders)){
    
    ### loop through activity classes
    for(j in 1:length(activity.classes)){
      
      ### UTILIZATION OF A SPLINE FOR SMOOTHING
      x <- ages
      
      ### y values for spline
      y <- tensor.object[g, j, ]
      
      #na.ind <- which(is.na(tensor.object[g, j, youngest_natsal:length(ages)]))
      na.ind <- which(is.na(tensor.object[g, j, ]))
      
      
      ### correct x and y values if necessary
      if(length(na.ind) > 0){
        y <- y[-na.ind]
        x <- x[-na.ind]
      }
      
      #compute spline on log-transformed data
      x.trans <- log(x-x[1]+1)
      spline.res <- smooth.spline(x.trans,y, df = df)
      smoother <- predict(spline.res, x.trans)

      tensor.object[g, j, ] <- smoother$y
      
    }
  }

  return(tensor.object)
}



