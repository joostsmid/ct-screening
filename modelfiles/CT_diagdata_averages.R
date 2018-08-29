### This function tabulates diagnosis data from N. Chandra (submitted to Eurosurveillance, 2016) 
### in age agroups 15-19, 20-24, 25-24 and 35-44 subdivided for the age classes used in the model.

get.CT.diagdata <- function(age.classes){
  
  agecategories.bounds <- c(15,20,25,35,45)
  agecategories <- c(15,20,25,35)
  years <- 2000:2012
  sex <- c("M","F")
  
  ##########################################################################################
  ### read in NCSP data (from N. Chandra et al (2016, submitted to Eurosurveillance))
  #### denominator: total # people in stratum, numbers obtained from Office of National Statistics
  ##########################################################################################
  
  diag_chandra <- read.table(file.path(path_data,"diagnosisdata_chandra.txt"), header=T, row.names=1)
  
  # minimum and maximum estimates given, expected value defined as average of min and max
  diag_chandra_mean <- to.tensor(NA, dims = list(sex=sex, screencatage=paste0('screencatage', agecategories),screening_y = years))
  diag_chandra_mean["M",1,] <- (diag_chandra$M_min_1519 + diag_chandra$M_max_1519)/2
  diag_chandra_mean["M",2,] <- (diag_chandra$M_min_2024 + diag_chandra$M_max_2024)/2
  diag_chandra_mean["M",3,] <- (diag_chandra$M_min_2534 + diag_chandra$M_max_2534)/2
  diag_chandra_mean["M",4,] <- (diag_chandra$M_min_3544 + diag_chandra$M_max_3544)/2
  diag_chandra_mean["F",1,] <- (diag_chandra$F_min_1519 + diag_chandra$F_max_1519)/2
  diag_chandra_mean["F",2,] <- (diag_chandra$F_min_2024 + diag_chandra$F_max_2024)/2
  diag_chandra_mean["F",3,] <- (diag_chandra$F_min_2534 + diag_chandra$F_max_2534)/2
  diag_chandra_mean["F",4,] <- (diag_chandra$F_min_3544 + diag_chandra$F_max_3544)/2
  
  ##########################################################################################
  ### put the inferred uptake values in age groups which are of interest for the rest of the model
  ### zero screening assumed for age classes older than 45 (calculate natsal_T_A(nat,c(45,60),nJ)$prev
  ### This gives 0.24% for males>45 in 2010, 0.20% for females>45 in 2010)
  ##########################################################################################
  
  n.age.classes <- length(age.classes)-1
  age.class.groups <- age.classes[1:n.age.classes]
  
  years.in.categories <- to.tensor(NA, dims = list(age=paste0('age', age.classes[1:(n.age.classes)]) ,screencatage=paste0('diagcatage', agecategories)))
  for (k in 1:n.age.classes){
    for (y in 1:4){
      years.in.categories[k,y] <- length(intersect( seq(age.classes[k],age.classes[k+1]-1), seq(agecategories.bounds[y],agecategories.bounds[y+1]-1) ))
    }
  }
  norm.years.in.categories <- years.in.categories / margin.tensor(years.in.categories, i=2) # fraction of time that individuals in age.classes[k] are in the screencatage[j] class
  norm.years.in.categories[is.na(norm.years.in.categories)] <- 0 #possible divisions by zero are changed to zeros
  
  diag <- margin.tensor(norm.years.in.categories * diag_chandra_mean, i=2)
  diag <- reorder.tensor(diag, c(2,1,3))
  
  return(diag)
  
}