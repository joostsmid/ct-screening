# Here we set up a series of procedures to calculate the probabilities with which 
# males (or females), a particular age a, engage in sexual partnerships with members of the 
# opposite sex of age a'. The resulting matrices are the used to calculate the 
# parameter rho in the genericSTImod.R script. 
# Calculations for age >=16 are based on reported partner ages of most recent, second most recent and third most recent partners
# Calculations for age <16 are based on reported partner ages of first sexual partners. For these, only data from respondents
# younger than 20 are used to prevent historical bias. For first partners, we assume that relationships were heterosexual.
###############################################################################

### load libraries 
library(spatial)
library(sn)
library(RColorBrewer)

#####################################################################################
### function that generates mixing matrix for females and males
#####################################################################################

create.mixing.matrices <-  function(nat, age.classes, verbose = TRUE){
  
  ############################ EXPLANATIONS OF THE VARIABLE NAMES IN THE DATASET ####################################################
  ### rafsmr is respondent's age in completed years at 1st sex with most recent partner
  ### r1ptage is partner's age on 1st occasion
  ### dage is respondent's age at interview, in years 
  ### rsex is respondent's sex (1 = male, 2 = female)
  ### agefsp is (estimated) partner age at moment of first sex
  ### total_wt is the respondent's weight to produced a statistically well-balanced data set: (some respondent's are over or
  ### underrepresented, and this needs to be adjusted)
  ###################################################################################################################################
  
  mostrecentparticipantage <- nat$rafsmr
  mostrecentpartnerage <- nat$r1ptage
  mostrecentparticipantage2 <- nat$rafs2mr
  mostrecentpartnerage2 <- nat$r1ptage2
  mostrecentparticipantage3 <- nat$rafs3mr
  mostrecentpartnerage3 <- nat$r1ptage3
  
  n.age.classes <- length(age.classes)-1
  sex.names <- c("M","F") # gender
  sex.code <-c(1,2) #coding used in NATSAL for gender
  partnersex.code <- c(2,1)
  
  mm.m <- matrix(data=0,nrow=n.age.classes,ncol=n.age.classes)
  mm.f <- matrix(data=0,nrow=n.age.classes,ncol=n.age.classes)
  
  ### filter the data set: no-nonsense conditions for partner age
  no.nonsense.cond1 <- mostrecentpartnerage > 0 & mostrecentpartnerage < 91 & !is.na(nat$total_wt)
  no.nonsense.cond2 <- mostrecentpartnerage2 > 0 & mostrecentpartnerage2 < 91 & !is.na(nat$total_wt)
  no.nonsense.cond3 <- mostrecentpartnerage3 > 0 & mostrecentpartnerage3 < 91 & !is.na(nat$total_wt)
  
  ### loop through all genders
  for(g in sex.code){
    
    ### filter the data set: only heterosexual relations
    gender.cond <- nat$rsex == sex.code[g]
    genderr1.cond <-  nat$r1sex == partnersex.code[g]
    genderr2.cond <-  nat$r1sex2 == partnersex.code[g]
    genderr3.cond <-  nat$r1sex3 == partnersex.code[g]
    
    for(i in 1:(n.age.classes)) {
      
      ### filter the data set: is the age of the respondent (when had sex) in the age band?
      ageband.resp.cond1 <- mostrecentparticipantage >= age.classes[i] & mostrecentparticipantage < age.classes[i+1]
      ageband.resp.cond2 <- mostrecentparticipantage2 >= age.classes[i] & mostrecentparticipantage2 < age.classes[i+1]
      ageband.resp.cond3 <- mostrecentparticipantage3 >= age.classes[i] & mostrecentparticipantage3 < age.classes[i+1]
      
      ### loop through all age classes of partner
      for(j in 1:n.age.classes) {
        
        #cat(g,age.classes[i],age.classes[i+1], '\n') #print to screen
        
        ### filter the data set: is the age of the partner in the age band?
        ageband.part.cond1 <- mostrecentpartnerage >= age.classes[j] & mostrecentpartnerage < age.classes[j+1]
        ageband.part.cond2 <- mostrecentpartnerage2 >= age.classes[j] & mostrecentpartnerage2 < age.classes[j+1]
        ageband.part.cond3 <- mostrecentpartnerage3 >= age.classes[j] & mostrecentpartnerage3 < age.classes[j+1]
        
        ### combine all conditions for specific gender and age band
        all.conds1 <- which(gender.cond & genderr1.cond & ageband.resp.cond1 & ageband.part.cond1 & no.nonsense.cond1)
        all.conds2 <- which(gender.cond & genderr2.cond & ageband.resp.cond2 & ageband.part.cond2 & no.nonsense.cond2)
        all.conds3 <- which(gender.cond & genderr3.cond & ageband.resp.cond3 & ageband.part.cond3 & no.nonsense.cond3)
        
        ### extract subset that matches specified conditions
        summed.weight <- sum(nat$total_wt[all.conds1]) + sum(nat$total_wt[all.conds2]) + sum(nat$total_wt[all.conds3]) 
        
        if (g==1){
          mm.m[i,j] <- summed.weight
        }
        else{
          mm.f[i,j] <- summed.weight
        }
        
      }
    }
    
  }
  
  rownames(mm.m) <- paste0('age', age.classes[-(n.age.classes+1)])
  colnames(mm.m) <- paste0('age', age.classes[-(n.age.classes+1)])
  rownames(mm.f) <- paste0('age', age.classes[-(n.age.classes+1)])
  colnames(mm.f) <- paste0('age', age.classes[-(n.age.classes+1)])
  
  return(list(mm.m, mm.f))
  
}

#####################################################################################
### function that generates mixing matrix for females and males but extends age partner to a desired (extended) range
### used for making plot of mixing matrix
#####################################################################################

create.mixing.matrices.extended <-  function(nat, age.classes, age.classes.partner, verbose = TRUE){
  
  ############################ EXPLANATIONS OF THE VARIABLE NAMES IN THE DATASET ####################################################
  ### rafsmr is respondent's age in completed years at 1st sex with most recent partner
  ### r1ptage is partner's age on 1st occasion
  ### dage is respondent's age at interview, in years 
  ### rsex is respondent's sex (1 = male, 2 = female)
  ### total_wt is the respondent's weight to produced a statistically well-balanced data set: (some respondent's are over or underrepresented, and this needs to be adjusted)
  ###################################################################################################################################
  
  mostrecentparticipantage <- nat$rafsmr
  mostrecentpartnerage <- nat$r1ptage
  mostrecentparticipantage2 <- nat$rafs2mr
  mostrecentpartnerage2 <- nat$r1ptage2
  mostrecentparticipantage3 <- nat$rafs3mr
  mostrecentpartnerage3 <- nat$r1ptage3
  
  n.age.classes <- length(age.classes)-1
  n.age.classes.partner <- length(age.classes.partner)-1
  sex.names <- c("M","F") # gender
  sex.code <-c(1,2) #coding used in NATSAL for gender
  partnersex.code <- c(2,1)
  
  mm.m <- matrix(data=0,nrow=n.age.classes,ncol=n.age.classes.partner)
  mm.f <- matrix(data=0,nrow=n.age.classes,ncol=n.age.classes.partner)
  
  ### filter the data set: no-nonsense conditions for partner age
  no.nonsense.cond1 <- mostrecentpartnerage > 0 & mostrecentpartnerage < 91 & !is.na(nat$total_wt)
  no.nonsense.cond2 <- mostrecentpartnerage2 > 0 & mostrecentpartnerage2 < 91 & !is.na(nat$total_wt)
  no.nonsense.cond3 <- mostrecentpartnerage3 > 0 & mostrecentpartnerage3 < 91 & !is.na(nat$total_wt)
  
  ### loop through all genders
  for(g in sex.code){
    
    ### filter the data set: only heterosexual relations
    gender.cond <- nat$rsex == sex.code[g]
    genderr1.cond <-  nat$r1sex == partnersex.code[g]
    genderr2.cond <-  nat$r1sex2 == partnersex.code[g]
    genderr3.cond <-  nat$r1sex3 == partnersex.code[g]
    
    for(i in 1:(n.age.classes)) {
      
      ### filter the data set: is the age of the respondent (when had sex) in the age band?
      ageband.resp.cond1 <- mostrecentparticipantage >= age.classes[i] & mostrecentparticipantage < age.classes[i+1]
      ageband.resp.cond2 <- mostrecentparticipantage2 >= age.classes[i] & mostrecentparticipantage2 < age.classes[i+1]
      ageband.resp.cond3 <- mostrecentparticipantage3 >= age.classes[i] & mostrecentparticipantage3 < age.classes[i+1]
      
      ### loop through all age classes of partner
      for(j in 1:n.age.classes.partner) {
        
        #cat(g,age.classes[i],age.classes[i+1], '\n') #print to screen
        
        ### filter the data set: is the age of the partner in the age band?
        ageband.part.cond1 <- mostrecentpartnerage >= age.classes.partner[j] & mostrecentpartnerage < age.classes.partner[j+1]
        ageband.part.cond2 <- mostrecentpartnerage2 >= age.classes.partner[j] & mostrecentpartnerage2 < age.classes.partner[j+1]
        ageband.part.cond3 <- mostrecentpartnerage3 >= age.classes.partner[j] & mostrecentpartnerage3 < age.classes.partner[j+1]
        
        ### combine all conditions for specific gender and age band
        all.conds1 <- which(gender.cond & genderr1.cond & ageband.resp.cond1 & ageband.part.cond1 & no.nonsense.cond1)
        all.conds2 <- which(gender.cond & genderr2.cond & ageband.resp.cond2 & ageband.part.cond2 & no.nonsense.cond2)
        all.conds3 <- which(gender.cond & genderr3.cond & ageband.resp.cond3 & ageband.part.cond3 & no.nonsense.cond3)
        
        ### extract subset that matches specified conditions
        summed.weight <- sum(nat$total_wt[all.conds1]) + sum(nat$total_wt[all.conds2]) + sum(nat$total_wt[all.conds3]) 
        
        if (g==1){
          mm.m[i,j] <- summed.weight
        }
        else{
          mm.f[i,j] <- summed.weight
        }
        
      }
    }
    
  }
  
  rownames(mm.m) <- paste0('age', age.classes[-(n.age.classes+1)])
  colnames(mm.m) <- paste0('age', age.classes.partner[-(n.age.classes.partner+1)])
  rownames(mm.f) <- paste0('age', age.classes[-(n.age.classes+1)])
  colnames(mm.f) <- paste0('age', age.classes.partner[-(n.age.classes.partner+1)])
  
  return(list(mm.m, mm.f))
  
}

### test 
# mm <- create.mixing.matrices(nat,age.classes)[[1]]

# JS010716: deze functie heb ik aangepast.
#####################################################################################
### function that smoothes age normalized mixing matrix by different methods
### Takes a normalized matrix (raw data about which individuals (male) mixed with which (female), categorized by age)
### and returns a smoothened matrix. 
### The first couple of methods smooth each row in a sequential fashion. 
### Two method smoothes the entire surface at once. 
#####################################################################################
smooth.mixing.matrix <- function(mixing.matrix, method = "spline", age.classes, verbose = TRUE){
  
  ### mixing matrix should have equal number of rows as it has columns
  ### it should be enumerated from 1:n.age.classes, the maximum number of rows, 
  ### which should correspond to the ages of the respondent or partner. 
  
  iMaxAge <- tail(age.classes,n=1)
  n.age.classes <- length(age.classes)-1
  age.class.groups <- age.classes[1:n.age.classes]
  
  ### define matrix with smoothened values
  smooth.mat <- matrix(0, nrow = n.age.classes, ncol = n.age.classes)
  
  ### loop through rows
  for(i in 1:n.age.classes){
    
    ### define the predictor and response variables
    y <- mixing.matrix[i,]
    #x <- 1:n.age.classes
    x <- age.class.groups
    
    if(method == "spline"){
      ### must not weight the values at the end of the distribution, because otherwise smoothing spline is zero
      spline.func <- smooth.spline( x = x, y = y, w = dnorm(1:n.age.classes, mean  = i, sd = 10))
      predictions <- predict(spline.func, x)
      smooth.mat[i, ] <- predictions$y
    }
    
    if(method == "nadaraya-watson"){
      
      smooth.func <- ksmooth(x,y, kernel = "normal", bandwidth = 5,n.points = length(x))
      predictions <- smooth.func 
      smooth.mat[i, ] <- predictions$y
    }
    
    if(method == "skew-normal"){
      
      ### run optimization
      smooth.func <- optim(par = c(xi = 15, omega = 1, alpha = 0), function(pars) RSS.skewnorm(pars, x = x, y = y) )
      predictions <- list(x = age.class.groups, y = dsn(x, smooth.func$par[1], smooth.func$par[2], smooth.func$par[3]) )
      smooth.mat[i, ] <- predictions$y
    }
    
  }
  
  if(method == "smooth.2d"){
    
    ### cannot work out how surf.ls works, due to lack of documenation
    #trend.surface <- surf.ls(np = 2, x = 1:n.age.classes, y = 1:n.age.classes, z = mixing.matrix)
    #smooth.mat <- predict(trend.surface, x = 1:n.age.classes, y = 1:n.age.classes)
    
    rownames(mixing.matrix) <- age.class.groups
    colnames(mixing.matrix) <- age.class.groups
    
    
    smooth.md <- as.data.frame(as.table(mixing.matrix))
    names(smooth.md) <- c("respondent.age", "partner.age", "frequency")
    
    
    smooth.md[,1] <- as.numeric(smooth.md[,1])
    smooth.md[,2] <- as.numeric(smooth.md[,2])
    
    tr.srf <- smooth.2d(Y = smooth.md$frequency, x = smooth.md[,1:2], nrow=n.age.classes, ncol=n.age.classes, surface = FALSE)
    #image(tr.srf)
    
    smooth.mat <- tr.srf
  }
  
  if(method == "sn"){
    
    ###mimi <- matrix(0, nrow = n.age.classes, ncol = n.age.classes) #JS, for Pjoint with sum(Pjoint)=1
    
    ### RSS function for a 2 dim aging function
    sn.RSS <- function(mixing.matrix, par, verbose = TRUE){
      
      asym = par[1]
      lin.fac = par[2] 
      alpha = par[3]
      
      RSS <- 0
      
      for(i in 1:n.age.classes){
        ###mimi[i,] <-mixing.matrix[i,]/sum(mixing.matrix[i,]) #JS, for Pjoint with sum(Pjoint)=1
        ###RSS <- RSS + sum((mimi[i,] - dsn(1:n.age.classes, xi = i + asym*i, omega = lin.fac*i, alpha = alpha))^2) #JS, for Pjoint with sum(Pjoint)=1
        RSS <- RSS + sum((mixing.matrix[i,] - dsn(1:n.age.classes, xi = i + asym*i, omega = lin.fac*i, alpha = alpha))^2)
      }
      
      if(verbose){
        cat("RSS is: ", RSS, '\n' )
      }
      
      return(RSS) 
    }
    
    res.optim <- optim(par = c(asym = 0.001, lin.fac = 3, alpha = 1), function(x) sn.RSS(mixing.matrix, x) )
    
    
    ### fill smoothing matrix
    for(i in 1:n.age.classes){
      smooth.mat[i,] <- dsn(1:n.age.classes, xi = i + res.optim$par[1]*i, omega = res.optim$par[2]*i, alpha = res.optim$par[3])
    }
    ###smooth.mat<-smooth.mat/sum(smooth.mat) #JS, for Pjoint with sum(Pjoint)=1
  }
  
  return(smooth.mat)
}


#### Here we make a 3D plot of the matrix 
plot.3D.matrix <- function(mix.matrix, respondent = 16:60, partner = 16:60, normalize = TRUE, verbose = TRUE){
  
  if(normalize){
    mixing.matrix <- normalize.mat(mix.matrix)
  }
  
  
  #### define variables 
  x <- respondent
  y <- partner
  
  ### response variable
  z <- mix.matrix[respondent, partner]
  
  ### colors 
  n.col = 100
  color <- heat.colors(n.col)
  
  nrz <- nrow(z)
  ncz <- ncol(z)
  # Compute the z-value at the facet centres
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  # Recode facet z-values into color indices
  facetcol <- cut(zfacet, n.col)
  
  
  persp(x, 
        y, 
        z = z, 
        #col = "lightgreen", 
        #col = surf.colors(as.matrix(z), n.col), 
        col = color[facetcol], 
        #shade = 0.4,
        theta = 30, 
        phi = 30,
        xlab = "Respondent Age [years]", 
        ylab = "Partner Age [years]", 
        zlab = "Conditional Probability Density", 
        ticktype = "detailed")
  
}


### test 
#plot.3D.matrix(smooth.mat)



### This function draws a contour plot of the data gathered from NATSAL 
### and converted into an age mixing matrix 
plot.contour.agemix <- function(norm.res, female = TRUE, iMinAge = 16, iMaxAge = 60, with.averages = TRUE){
  
  ### define the boundaries for respondents and partners
  respondent <- iMinAge:iMaxAge
  partner <- iMinAge:iMaxAge
  
  {if(class(norm.res) == "list"){
    ### extract information of norm.res
    mf <- norm.res[[1]]
    mf1 <- norm.res[[2]]
    mf2 <- norm.res[[3]]
  }
    else{
      mf <- norm.res
      with.averages = FALSE
    }}
  
  #### labels 
  if(female){
    main.lab = "Male partners of female respondents"
    x.lab = "age of female respondent [years]"
    y.lab = "age of male partner [years]"
  }
  else{
    main.lab = "Female partners of male respondents"
    x.lab = "age of male respondent [years]"
    y.lab = "age of female partner [years]"
  }
  
  ### contour plot
  if(with.averages){
    ### with the mean and median lines
    filled.contour(respondent,partner,mf[respondent,partner],
                   col=gray(19:1/19),
                   plot.title = title(main= main.lab,xlab= x.lab, ylab = y.lab),
                   key.title = title(main="density"), 
                   plot.axes={axis(1);axis(2);abline(0,1,lwd=2);lines(respondent,mf1[respondent],lwd=2,lty=3);lines(respondent,mf2[respondent],lwd=2,lty=2)})
  }
  else{
    ### without the lines
    filled.contour(respondent, partner, mf[respondent, partner], col = gray(19:1/19), 
                   plot.title = title(main=main.lab,xlab=x.lab,ylab=y.lab),
                   key.title = title(main="density"),
                   plot.axes={axis(1);axis(2);abline(0,1,lwd=2)})
  }
  
}

### test 
#plot.contour.agemix(res)
#plot.contour.agemix(smooth.mat)  

### This function draws a contour plot of the data gathered from NATSAL 
### and converted into an age mixing matrix 
plot.image.agemix <- function(norm.res, female = TRUE, iMinAge = 16, iMaxAge = 60, with.averages = TRUE, cex = 1.8, line.dis = 3){
  
  
  ### define the boundaries for respondents and partners
  respondent <- iMinAge:iMaxAge
  partner <- iMinAge:iMaxAge
  
  ### extract information of norm.res
  {if(class(norm.res) == "list"){
    ### extract information of norm.res
    mf <- norm.res[[1]]
    mf1 <- norm.res[[2]]
    mf2 <- norm.res[[3]]
  }
    else{
      mf <- norm.res
      with.averages = FALSE
    }}
  
  #### labels 
  if(female){
    main.lab = "Male partners of female respondents"
    x.lab = "age of female respondent [years]"
    y.lab = "age of male partner [years]"
  }
  else{
    main.lab = "Female partners of male respondents"
    x.lab = "age of male respondent [years]"
    y.lab = "age of female partner [years]"
  }
  
  ### image plot
  image(respondent,partner,mf[respondent,partner],col=gray(20:1/20),
        main= main.lab,
        xlab= "", #x.lab,
        ylab= "", #y.lab,
        cex.axis = cex,
        cex.lab = cex)
  abline(0,1,lwd=2)
  
  if(with.averages){
    ### do not understand these lines
    lines(respondent,mf1[respondent],lwd=2,lty=3)
    lines(respondent,mf2[respondent],lwd=2,lty=2)
  }
  
  
  mtext(x.lab, side = 1, line = line.dis, cex = cex)
  mtext(y.lab, side = 2, line = line.dis, cex = cex)
}

### test 
#plot.image.agemix(norm.res = normalize.mat(f.mixing.matrix))