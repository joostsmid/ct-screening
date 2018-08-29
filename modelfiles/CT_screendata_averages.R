### This function tabulates screening data from N. Chandra (submitted to Eurosurveillance, 2016) with NATSAL-3 data
### to infer proportions of screening uptake among sexually actives from 2000 to 2012 in age agroups 15-19, 20-24,
### 25-24 and 35-44 subdivided for the age classes used in the model.

get.CT.screendata <- function(age.classes){
  
  screencategories.bounds <- c(16,20,25,35,45)
  screencategories <- c(16,20,25,35)
  years_screening <- 2000:2012
  sex <- c("M","F")
  
  ##########################################################################################
  ### read in NCSP data (from N. Chandra et al (2016, submitted to Eurosurveillance))
  #### denominator: total # people in stratum, numbers obtained from Office of National Statistics
  ##########################################################################################
  
  upt_NCSP <- read.table(file.path(path_data,"screeningdata_chandra.txt"), header=T, row.names=1) / 100 # per capita number of tests
  
  # minimum and maximum estimates given, expected value defined as average of min and max
  upt_NCSP_mean <- to.tensor(NA, dims = list(sex=sex, screencatage=paste0('screencatage', screencategories),screening_y = years_screening))
  upt_NCSP_mean["M",1,] <- (upt_NCSP$M_min_1519 + upt_NCSP$M_max_1519)/2
  upt_NCSP_mean["M",2,] <- (upt_NCSP$M_min_2024 + upt_NCSP$M_max_2024)/2
  upt_NCSP_mean["M",3,] <- (upt_NCSP$M_min_2534 + upt_NCSP$M_max_2534)/2
  upt_NCSP_mean["M",4,] <- (upt_NCSP$M_min_3544 + upt_NCSP$M_max_3544)/2
  upt_NCSP_mean["F",1,] <- (upt_NCSP$F_min_1519 + upt_NCSP$F_max_1519)/2
  upt_NCSP_mean["F",2,] <- (upt_NCSP$F_min_2024 + upt_NCSP$F_max_2024)/2
  upt_NCSP_mean["F",3,] <- (upt_NCSP$F_min_2534 + upt_NCSP$F_max_2534)/2
  upt_NCSP_mean["F",4,] <- (upt_NCSP$F_min_3544 + upt_NCSP$F_max_3544)/2
  
  
  ##########################################################################################
  ### put the inferred uptake values in age groups which are of interest for the rest of the model
  ### zero screening assumed for age classes older than 45 (calculate natsal_T_A(nat,c(45,60),nJ)$prev
  ### This gives 0.24% for males>45 in 2010, 0.20% for females>45 in 2010)
  ##########################################################################################
  
  n.age.classes <- length(age.classes)-1
  age.class.groups <- age.classes[1:n.age.classes]
  
  years.in.categories <- to.tensor(NA, dims = list(age=paste0('age', age.classes[1:(n.age.classes)]) ,screencatage=paste0('screencatage', screencategories)))
  for (k in 1:n.age.classes){
    for (y in 1:4){
      years.in.categories[k,y] <- length(intersect( seq(age.classes[k],age.classes[k+1]-1), seq(screencategories.bounds[y],screencategories.bounds[y+1]-1) ))
    }
  }
  norm.years.in.categories <- years.in.categories / margin.tensor(years.in.categories, i=2) # fraction of time that individuals in age.classes[k] are in the screencatage[j] class
  norm.years.in.categories[is.na(norm.years.in.categories)] <- 0 #possible divisions by zero are changed to zeros
  
  upt <- margin.tensor(norm.years.in.categories * upt_NCSP_mean, i=2)
  upt <- reorder.tensor(upt, c(2,1,3))
  
  return(upt)
  
}