# Here we define functions that allow to convert cpp based vectors (vectors as passed on 
# to STImod.cpp function) into tensors and back. 
#
# Authors: Victor Garcia and Christian Althaus, University of Bern, October 2015
###############################################################################

#function to convert variable (object) name into String
varname2string <- function(v1) {
  deparse(substitute(v1))
}

# varname2string(I)

#### Forward conversion: TENSOR to CPP-VECTOR
################################## function to convert tensors into useful vector represenation for cpp passing ##########################################
#### how to convert from a tensor to a vector such that the .cpp file can read stuff in correctly
# tens.to.vec.cpp <- function(tens, tens_name = "test", with.names = TRUE, verbose = FALSE){
tens.to.vec.cpp <- function(tens, with.names = TRUE, verbose = FALSE){ #changed JS130116
  
  my.s.vec <- numeric()
  my.s.counter <- 1
  tens_name <- as.character(substitute(tens))

  ### for one dimensional tensors
  if(length(dim(tens)) == 1){
    
    for(x1 in 1:dim(tens)[1]){
      my.s.vec[my.s.counter] <- tens[x1]
      
      if(with.names){
        names(my.s.vec)[my.s.counter] <- paste0(c(tens_name, dimnames(tens)[[1]][x1]), collapse = "_")
      }
        
      if(verbose){
        cat("looping through x1:", x1, '\n')
      }
      
      my.s.counter = my.s.counter + 1
    }
    
  }
  
  ### for two dimensional tensors
  if(length(dim(tens)) == 2){
    
    for(x1 in 1:dim(tens)[1]){
      for(x2 in 1:dim(tens)[2]){
        my.s.vec[my.s.counter] <- tens[x1,x2]
        
        if(with.names){
          names(my.s.vec)[my.s.counter] <- paste0(c(tens_name,dimnames(tens)[[1]][x1],dimnames(tens)[[2]][x2]), collapse = "_")
        }
        
        if(verbose){
          cat("looping through x1:", x1, " x2: ", x2, '\n')
        }
        
        my.s.counter = my.s.counter +1
      }
    }
    
  }
  
  ### For three dimensional tensors
  if(length(dim(tens)) == 3){
    
    for(x1 in 1:dim(tens)[1]){
      for(x2 in 1:dim(tens)[2]){
        for(x3 in 1:dim(tens)[3]){
          my.s.vec[my.s.counter] <- tens[x1,x2,x3]
          
          if(with.names){
            names(my.s.vec)[my.s.counter] <- paste0(c(tens_name,dimnames(tens)[[1]][x1],dimnames(tens)[[2]][x2],dimnames(tens)[[3]][x3]), collapse = "_")
          }
          
          if(verbose){
            cat("looping through x1:", x1, " x2: ", x2, " x3: ", x3, '\n')
          }
          
          my.s.counter = my.s.counter +1
        }
      }
    }
    
  }
  
  ### four dimensional tensor
  if(length(dim(tens)) == 4){
    
    for(x1 in 1:dim(tens)[1]){
      for(x2 in 1:dim(tens)[2]){
        for(x3 in 1:dim(tens)[3]){
          for(x4 in 1:dim(tens)[4]){
            my.s.vec[my.s.counter] <- tens[x1,x2,x3, x4]
            
            if(with.names){
              names(my.s.vec)[my.s.counter] <- paste0(c(tens_name,dimnames(tens)[[1]][x1],dimnames(tens)[[2]][x2],dimnames(tens)[[3]][x3],dimnames(tens)[[4]][x4]), collapse = "_")
            }
            
            if(verbose){
              cat("looping through x1:", x1, " x2: ", x2, " x3: ", x3, " x4: ", x4,'\n')
            }
            
            my.s.counter = my.s.counter +1
          }
        }
      }
    }
    
  }
  
  ### five dimensional tensor
  if(length(dim(tens)) == 5){
    
    for(x1 in 1:dim(tens)[1]){
      for(x2 in 1:dim(tens)[2]){
        for(x3 in 1:dim(tens)[3]){
          for(x4 in 1:dim(tens)[4]){
            for(x5 in 1:dim(tens)[5]){
              my.s.vec[my.s.counter] <- tens[x1,x2,x3, x4, x5]
              my.s.counter = my.s.counter +1
            }
          }
        }
      }
    }
    
  }
  
  
  ### six-dimensional tensor
  if(length(dim(tens)) == 6){
    
    for(x1 in 1:dim(tens)[1]){
      for(x2 in 1:dim(tens)[2]){
        for(x3 in 1:dim(tens)[3]){
          for(x4 in 1:dim(tens)[4]){
            for(x5 in 1:dim(tens)[5]){
              for(x6 in 1:dim(tens)[6]){
                my.s.vec[my.s.counter] <- tens[x1,x2,x3, x4, x5, x6]
                my.s.counter = my.s.counter +1
              }
            }
          }
        }
      }
    }
    
  }
  
  ### six-dimensional tensor
  if(length(dim(tens)) == 7){
    
    for(x1 in 1:dim(tens)[1]){
      for(x2 in 1:dim(tens)[2]){
        for(x3 in 1:dim(tens)[3]){
          for(x4 in 1:dim(tens)[4]){
            for(x5 in 1:dim(tens)[5]){
              for(x6 in 1:dim(tens)[6]){
                for(x7 in 1:dim(tens)[7]){
                  my.s.vec[my.s.counter] <- tens[x1, x2, x3, x4, x5, x6, x7]
                  my.s.counter = my.s.counter +1
                }
              }
            }
          }
        }
      }
    }
    
  }
  
  return(my.s.vec)
}



#### conversion of an initial vector to a series of tensors representing the different compartments 
initvec.to.tensor.cpp <- function(init.cpp, par.cond = list(age_classes, nmb_act_classes),  verbose = FALSE, plot.progress = TRUE){

  ################################## convert the init.cpp vector into tensors ##################################
  nmb_age_classes <- length(par.cond$age_classes)
  
  if(verbose){
    cat("converting the vector if initial.cpp conditions into a series of tensor class objects", '\n')
  }
  
  ##############################################################################################################################
  ### auxiliary parameter sex/gender
  sex <- c("M","F") ### gender
  sex_length <- 2
  
  ##############################################################################################################################
  ### Define a vector of virgins, structured by age, sex
  dimnames1 <- list(sex = sex, age = paste0('age', par.cond$age_classes))
  nElements <-  nmb_age_classes*length(sex)
  
  U <- to.tensor(rep(0, nElements), dims  = dimnames1)
  
  ##############################################################################################################################
  ### Define a vector of susceptibles, structured by age, sex
  nJ <- par.cond$nmb_act_classes ### number of activity classes
  dimnames1 <- list(sex=sex, j=paste0('j',1:nJ), age=paste0('age', par.cond$age_classes))
  nElements <-  nJ*nmb_age_classes*length(sex)
  
  S <- to.tensor(rep(0, nElements), dims  = dimnames1)
  
  ##############################################################################################################################
  ### Define a vector of infecteds I_A, structured by age, sex, activity class
  
  I_A <- to.tensor(rep(0, nElements), dims  = dimnames1)
  ##############################################################################################################################
  ### Define a vector of infecteds I_S, structured by age, sex, activity class
  
  I_S <- to.tensor(rep(0, nElements), dims  = dimnames1)
  
  ##############################################################################################################################
  ### Define a vector of recovereds R, structured by age, sex, activity class
  ### gender
  R <- to.tensor(rep(0, nElements), dims  = dimnames1)
  
  ##############################################################################################################################
  ### Define a vector of infecteds I_A, structured by age, sex, activity class
  
  I_A2 <- to.tensor(rep(0, nElements), dims  = dimnames1)
  ##############################################################################################################################
  ### Define a vector of infecteds I_S, structured by age, sex, activity class
  
  I_S2 <- to.tensor(rep(0, nElements), dims  = dimnames1)
  
  
  ##############################################################################################################################
  ### Define a vector of recovereds vaccinated, structured by age, sex, activity class
  ### gender
  D <- to.tensor(rep(0, nElements), dims  = dimnames1)
  
  
  ###################### Fill up the tensors #################################################################################
  tensor_entry_numbers = sex_length * par.cond$nmb_act_classes * nmb_age_classes;
  
  ## define new counter
  counter = 1
  
  ## read in U: virgins
  for(g in 1:sex_length){
    for (j in 1:nmb_age_classes) {
      U [g,j] = as.numeric(init.cpp[counter]);
      counter = counter + 1;
    }
  } 
  
  ## read in the initial conditions passed on by the user call
  for(g in 1:sex_length){
    for (j in 1:par.cond$nmb_act_classes) {
      for (a in 1:nmb_age_classes) {
        
        if(verbose){
          cat("looping through x1:", g, " x2: ", j, " x3: ", a, '\n')
        }
        
        S [g,j,a] = as.numeric(init.cpp[counter]);
        I_A [g,j,a] = as.numeric(init.cpp[counter + 1*tensor_entry_numbers]);
        I_S [g,j,a] = as.numeric(init.cpp[counter + 2*tensor_entry_numbers]);
        R [g,j,a] = as.numeric(init.cpp[counter + 3*tensor_entry_numbers]);
        I_A2 [g,j,a] = as.numeric(init.cpp[counter + 4*tensor_entry_numbers]);
        I_S2 [g,j,a] = as.numeric(init.cpp[counter + 5*tensor_entry_numbers]);
        D [g,j,a] = as.numeric(init.cpp[counter + 6*tensor_entry_numbers]);
        counter = counter + 1;
      }
    }
  } 
  
  
  if(verbose){
    
    cat("Conversion into tensors finalized", '\n')
    
    if(counter - 1 + 6*tensor_entry_numbers != length(init.cpp)){
      
      cat("last entry of init.cpp used: ", counter -1, '\n')
      cat("length of cpp-based vector to be converted: ", length(init.cpp), '\n')
      cat("ERROR COMMITED DURING CONVERSION", '\n')
    }
    
  }
  
  return(list(U = U, S = S, I_A = I_A, I_S = I_S, R = R, I_A2 = I_A2, I_S2 = I_S2, D = D))
  
}