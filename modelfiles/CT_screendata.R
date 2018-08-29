### This function tabulates test data from N. Chandra (submitted to Eurosurveillance, 2016): Number of tests pp
### in age agroups 15-19, 20-24,25-24 and 35-44 subdivided for the age classes used in the model

get.CT.screendata <- function(age.classes){
  
  agecategories.bounds <- c(16,20,25,35,45)
  agecategories <- c(16,20,25,35)
  years_screening <- 2000:2012
  sex <- c("M","F")
  minmax <- c("min","max")
  
  ##########################################################################################
  ### read in NCSP data (from N. Chandra et al (2016, submitted to Eurosurveillance))
  #### denominator: total # people in stratum, numbers obtained from Office of National Statistics
  ##########################################################################################
  
  upt_chandra_data <- read.table(file.path(path_data,"screeningdata_chandra.txt"), header=T, row.names=1) / 100 # per capita number of tests
  
  # minimum and maximum estimates given, expected value defined as average of min and max
  upt_chandra <- to.tensor(NA, dims = list(sex=sex, screencatage=paste0('screencatage', agecategories),minmax=minmax,screening_y = years_screening))
  upt_chandra["M",1,1,] <- upt_chandra_data$M_min_1519
  upt_chandra["M",1,2,] <- upt_chandra_data$M_max_1519
  upt_chandra["M",2,1,] <- upt_chandra_data$M_min_2024
  upt_chandra["M",2,2,] <- upt_chandra_data$M_max_2024
  upt_chandra["M",3,1,] <- upt_chandra_data$M_min_2534
  upt_chandra["M",3,2,] <- upt_chandra_data$M_max_2534
  upt_chandra["M",4,1,] <- upt_chandra_data$M_min_3544
  upt_chandra["M",4,2,] <- upt_chandra_data$M_max_3544
  upt_chandra["F",1,1,] <- upt_chandra_data$F_min_1519
  upt_chandra["F",1,2,] <- upt_chandra_data$F_max_1519
  upt_chandra["F",2,1,] <- upt_chandra_data$F_min_2024
  upt_chandra["F",2,2,] <- upt_chandra_data$F_max_2024
  upt_chandra["F",3,1,] <- upt_chandra_data$F_min_2534
  upt_chandra["F",3,2,] <- upt_chandra_data$F_max_2534
  upt_chandra["F",4,1,] <- upt_chandra_data$F_min_3544
  upt_chandra["F",4,2,] <- upt_chandra_data$F_max_3544
  
  ##########################################################################################
  ### put the inferred uptake values in age groups which are of interest for the rest of the model
  ### zero screening assumed for age classes older than 45 (calculate natsal_T_A(nat,c(45,60),nJ)$prev
  ### This gives 0.24% for males>45 in 2010, 0.20% for females>45 in 2010)
  ##########################################################################################
  
  n.age.classes <- length(age.classes)-1
  age.class.groups <- age.classes[1:n.age.classes]
  
  years.in.categories <- to.tensor(NA, dims = list(age=paste0('age', age.classes[1:(n.age.classes)]) ,screencatage=paste0('screencatage', agecategories)))
  for (k in 1:n.age.classes){
    for (y in 1:4){
      years.in.categories[k,y] <- length(intersect( seq(age.classes[k],age.classes[k+1]-1), seq(agecategories.bounds[y],agecategories.bounds[y+1]-1) ))
    }
  }
  norm.years.in.categories <- years.in.categories / margin.tensor(years.in.categories, i=2) # fraction of time that individuals in age.classes[k] are in the screencatage[j] class
  norm.years.in.categories[is.na(norm.years.in.categories)] <- 0 #possible divisions by zero are changed to zeros
  
  upt <- margin.tensor(norm.years.in.categories * upt_chandra, i=2)
  upt <- reorder.tensor(upt, c(2,1,3,4))
  
  return(upt)
  
}