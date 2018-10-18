

### source relevant scripts
source(file.path(path_in, "loadlibraries.R"))
source(file.path(path_in, "my_mcmcMH2.R")) #use this in case mcmcMH is not working: there seems to be a problem with the version in the fitR package
source(file.path(path_in, "generate_cpp.R")) # functions to calculate the differentials using a .cpp file and create an R compatible function thereof
source(file.path(path_in, "cppvectorconv.R")) # script to convert tensors into vectors and back
# source(file.path(path_in, "load_chlamydia_params_natsal.R")) # load model parameters using a Natsal dataset
source(file.path(path_in, "load_chlamydia_params.R")) # load saved model parameters
# source(file.path(path_in, "natsal_CTprevalences.R")) # load function to compute prevalences using a Natsal dataset

#############################################################################################################
### compile c++ model
#############################################################################################################

generate.cpp.code(parms = parameters)

#system(paste("R CMD SHLIB ", filename, sep = "")) # compile for R
system("R CMD SHLIB STIreduced.cpp") # compile for R

if (R.Version()$os=="mingw32"){
  try(dyn.load("STIreduced.dll"))
} else {
  try(dyn.load("STIreduced.so"))
}

#############################################################################################################
### prevalence epidata (from Natsal 2 and 3)
#############################################################################################################

#make natsal 2 epidata
# natsal2.prevs.all  <- get.natsal.CTprevalences(nat2) # compute prevalence using a Natsal-2 dataset
load(file.path(path_data,"prevs_nat2.Rda")) # load saved prevalence from Natsal-2 data
natsal2.meannumerator.prevs <- apply(natsal2.prevs.all$numerator,MARGIN=c(1,2),mean)
natsal2.meandenominator.prevs <- apply(natsal2.prevs.all$denominator,MARGIN=c(1,2),mean)

#make natsal 3 epidata
# natsal3.prevs.all  <- get.natsal.CTprevalences(nat3) # compute prevalence using a Natsal-3 dataset
load(file.path(path_data,"prevs_nat3.Rda")) # load saved prevalence from Natsal-3 data
natsal3.meannumerator.prevs <- apply(natsal3.prevs.all$numerator,MARGIN=c(1,2),mean)
natsal3.meandenominator.prevs <- apply(natsal3.prevs.all$denominator,MARGIN=c(1,2),mean)

#combine data
dimnames.epidata <- list(sex=c("M","F"), age=paste0('age', age.class.groups), quotient=c("num","denom"))
epidata_2000 <- to.tensor(c(round(as.vector(natsal2.meannumerator.prevs)),as.vector(natsal2.meandenominator.prevs)), dims = dimnames.epidata)
epidata_2011 <- to.tensor(c(round(as.vector(natsal3.meannumerator.prevs)),as.vector(natsal3.meandenominator.prevs)), dims = dimnames.epidata)

### diagnosis epidata
diagdata <- read.table(file.path(path_data,diagdatafile), header=T, row.names=1)
agecategories_diag <- c(15,20,25,35)
diag <- to.tensor(NA, dims = list(sex=sex, age=paste0('age', agecategories_diag),screening_y = as.numeric(rownames(diagdata))))
for (a in 1:length(agecategories_diag)){
  for (y in 1:dim(diagdata)[1]){
    diag[1,a,y] <- diagdata[y,a]
    diag[2,a,y] <- diagdata[y,a+4]
  }
}

#we have epidata at time burnintime + 0 (we have screening data from 2000, NATSAL-2 data is of 2000) and at time burnintime +11 (NATSAL-3; 2011)
epidata <- list(time_0 = epidata_2000,time_11 = epidata_2011, epidata_diag=diag)

##############################################################################
# BAYESIAN ESTIMATION OF beta USING MH ALGORITHM
##############################################################################

#statistics prior
beta.prior.min <- 0
beta.prior.max <- 1
epsilon.prior.min <- 0
epsilon.prior.max <- 1
gamma.prior.min <- 0
gamma.prior.max <- 10
kappa.prior.min <- 0
kappa.prior.max <- 1
omega.A.prior.min <- 0
omega.A.prior.max <- 10
omega.S.prior.min <- 0
omega.S.prior.max <- 10
fsymp.M.prior.min <- 0
fsymp.M.prior.max <- 1
fsymp.F.prior.min <- 0
fsymp.F.prior.max <- 1
treat.prior.min <- 0
treat.prior.max <- 20
eta1.prior.min <- 1
eta1.prior.max <- 10
eta2.prior.min <- -2
eta2.prior.max <- 1
r.prior.min <- 1
r.prior.max <- 1000

#used for the priors
beta.prior.shape1 <- 6
beta.prior.shape2 <- 4
epsilon.prior.shape1 <- 8
epsilon.prior.shape2 <- 2
gamma.prior.mean <- 365/433 #Althaus, C. L., Heijne, J. C. M., Roellin, A. & Low, N. Transmission dynamics of Chlamydia trachomatis affect the impact of screening programmes. Epidemics-Neth 2, 123-131, doi:10.1016/j.epidem.2010.04.002 (2010
gamma.prior.sd <- ((365/433-365/447) + (365/420-365/433) / 2) / 1.96 #Althaus, C. L., Heijne, J. C. M., Roellin, A. & Low, N. 
omega.A.prior.mean <- 0.921*(1-0.194) #Batteiger, JID 2010:201
omega.A.prior.sd <- (1-0.194)*((0.921-0.899) + (0.96-0.921) / 2) / 1.96 #Batteiger, JID 2010:201
omega.S.prior.mean <- 0.921*(1-0.194) #Batteiger, JID 2010:201
omega.S.prior.sd <- (1-0.194)*((0.921-0.899) + (0.96-0.921) / 2) / 1.96 #Batteiger, JID 2010:201
treat.prior.mean <- 365/33 # Bethan Davies, Sarah-Jane Anderson, [...], and Helen Ward, How robust are the natural history parameters used in chlamydia transmission dynamic models? A systematic review
treat.prior.sd <- ((365/33-365/40) + (365/30-365/33) / 2) / 1.96

n.unknown <- length(testmodel) # number of parameters that need to be inferred

#define an model of type mcmc.model
STI_MCMC_model <- list(name="generic STI model",

                       simulate = function (theta, times) # SIMULATE THE MODEL OUTCOMES
                       {

                         #sample from parameters
                         parameters$beta["M"] <- theta[["beta"]]
                         parameters$beta["F"] <- theta[["beta"]]
                         parameters$epsilon <- theta[["epsilon"]]
                         parameters$gamma["M"]  <- theta[["gamma"]]
                         parameters$gamma["F"]  <- theta[["gamma"]]
                         parameters$kappa["M"] <- theta[["kappa"]]
                         parameters$kappa["F"] <- theta[["kappa"]]
                         omega.A <- theta[["omega.A"]]
                         omega.S <- theta[["omega.S"]]
                         parameters$omega <- c(omega.A,omega.S)
                         fsymp.M <- theta[["fsymp.M"]]
                         fsymp.F <- theta[["fsymp.F"]]
                         parameters$fsymp <- c(fsymp.M,fsymp.F)
                         parameters$treat <- theta[["treat"]]
                         eta1 <- theta[["eta1"]]
                         eta2 <- 10^theta[["eta2"]]
                         parameters$eta <- c(eta1,eta2)

                         parameters.cpp <- list(beta = parameters$beta,
                                                epsilon = parameters$epsilon,
                                                gamma = parameters$gamma,
                                                kappa = parameters$kappa,
                                                omega = parameters$omega,
                                                fsymp = parameters$fsymp,
                                                treat = parameters$treat,
                                                eta = parameters$eta)

                         parameters.cpp <- unlist(parameters.cpp)

                         trajectory <- data.frame(ode( y = init.cpp, times=times,
                                                       func = function(x,y,z) SIR.cpp(x,y,z),
                                                       parms = parameters.cpp, method = "bdf"))

                         names(trajectory) <- c("time",names(init.cpp))

                         # compute number of new diagnoses per year instead of cumulative number
                         for (k in grep("D_", names(trajectory))){
                           D_name <- names(trajectory)[k]
                           mutatetext <- paste0("mutate(trajectory,", D_name, "=c(0,diff(", D_name, ")))")
                           trajectory <- eval(parse(text=mutatetext))
                         }
                         # compute number of new cases per year instead of cumulative number
                         for (k in grep("Inc_", names(trajectory))){
                           Inc_name <- names(trajectory)[k]
                           mutatetext <- paste0("mutate(trajectory,", Inc_name, "=c(0,diff(", Inc_name, ")))")
                           trajectory <- eval(parse(text=mutatetext))
                         }

                         return(trajectory)

                       },

                       dprior=function (theta, log = FALSE) # COMPUTE THE DENSITIES OF A VALUE IN THE PRIOR DISTRIBUTION
                       {

                         #log.prior.beta <- dunif(theta[["beta"]], min = beta.prior.min, max = beta.prior.max, log = TRUE)
                         log.prior.beta <- dbeta(theta[["beta"]], shape1 = beta.prior.shape1, shape2 = beta.prior.shape2, log = TRUE)
                         #log.prior.epsilon <- dunif(theta[["epsilon"]], min = epsilon.prior.min, max = epsilon.prior.max, log = TRUE)
                         log.prior.epsilon <- dbeta(theta[["epsilon"]], shape1 = epsilon.prior.shape1, shape2 = epsilon.prior.shape2, log = TRUE)
                         log.prior.gamma <- dnorm(theta[["gamma"]], mean = gamma.prior.mean, sd = gamma.prior.sd, log = TRUE)
                         log.prior.kappa <- dunif(theta[["kappa"]], min = kappa.prior.min, max = kappa.prior.max, log = TRUE)
                         log.prior.omega.A <- dnorm(theta[["omega.A"]], mean = omega.A.prior.mean, sd = omega.A.prior.sd, log = TRUE)
                         log.prior.omega.S <- dnorm(theta[["omega.S"]], mean = omega.S.prior.mean, sd = omega.S.prior.sd, log = TRUE)
                         log.prior.fsymp.M <- dunif(theta[["fsymp.M"]], min = fsymp.M.prior.min, max = fsymp.M.prior.max, log = TRUE)
                         log.prior.fsymp.F <- dunif(theta[["fsymp.F"]], min = fsymp.F.prior.min, max = fsymp.F.prior.max, log = TRUE)
                         log.prior.treat <- dnorm(theta[["treat"]], mean = treat.prior.mean, sd = treat.prior.sd, log = TRUE)
                         log.prior.eta1 <- dunif(theta[["eta1"]], min = eta1.prior.min, max = eta1.prior.max, log = TRUE)
                         log.prior.eta2 <- dunif(theta[["eta2"]], min = eta2.prior.min, max = eta2.prior.max, log = TRUE)

                         log.prior.r <- dgamma(theta[["r"]], shape = 100^2, scale=100/5000, log = TRUE)
                         #mu=100;v=5000;x<-seq(0,1000,1);plot(x,dgamma(x,mu^2/v,mu/v), type="l")

                         log.sum <- sum(log.prior.beta) + sum(log.prior.epsilon) + sum(log.prior.gamma) + sum(log.prior.kappa ) +
                           sum(log.prior.omega.A) +sum(log.prior.omega.S) +  sum(log.prior.fsymp.M) + sum(log.prior.fsymp.F) +
                           sum(log.prior.treat) + sum(log.prior.eta1) + sum(log.prior.eta2) + log.prior.r
                         if (log.sum==-Inf) cat("wrong priors given to mcmc algorithm; priors:",theta, "\n")
                         return(ifelse(log, log.sum, exp(log.sum)))
                       },

                       dPointObsPrev = function (data.infecteds, data.agegroupsizes, model.prev, theta, log = FALSE)  # COMPUTE THE LOG LIKELIHOOD OF OBSERVING PREVALENCE DATA
                       {
                         return(dbinom(x = data.infecteds, size = data.agegroupsizes, prob = model.prev, log = log)) #point log likelihood
                       },


                       dPointObsDiag = function (data.diag, model.diag, theta, log = FALSE)  # COMPUTE THE LOG LIKELIHOOD OF OBSERVING DIAGNOSIS DATA
                       {
                         #return(dnbinom(x = round(data.diag) , mu = model.diag * 100000, size= 1/theta[["r"]], log = log)) #size param negbin unscaled
                         return(dnbinom(x = round(data.diag) , mu = model.diag * 100000, size= model.diag * 100000 / (theta[["r"]]-1), log = log)) #size param negbin scaled
                       }

)


#other functions: function that calculates the likelihood of a trajectory given data
my_dTrajObs <- function (fitmodel, theta, data)
{
  if (is.null(dimnames(data[[1]])$j)){
    nJ_data <- 1
  } else {
    nJ_data <- length(dimnames(data[[1]])$j)
  }

  age<-paste0('age', age.class.groups)

  times <- c(0, (parameters$burnintime - parameters$testperiod_in_burnin):(parameters$burnintime+dim(diagdata)[1]))
  # we should start at burnintime - testperiod_in_burnin because diagnoses have to be mutated from that point on

  traj <- fitmodel$simulate(theta, times)
  lppd <- c()
  par.cond <- list(age_classes =age.class.groups, nmb_act_classes = nJ)

  ###################################################
  # compute density for prevalence data
  ###################################################

  years.prevdata <- c(2000, 2011) - 2000 + parameters$burnintime #2000, 2011
  for (k in 1:2) { # loop for different time points (for us: NATSAL-2 and NATSAL-3)

    # model output
    traj.natsal <- traj[traj$time==years.prevdata[k], -1]

    tensors <- initvec.to.tensor.cpp(traj.natsal, par.cond = par.cond)

    U.mod <- tensors$U
    S.mod <- tensors$S
    I_A.mod <- tensors$I_A
    I_S.mod <- tensors$I_S
    R.mod <- tensors$R
    I_A.mod2 <- tensors$I_A2
    I_S.mod2 <- tensors$I_S2

    if (nJ_data==1){

      S.mod <- margin.tensor(S.mod,i=2)
      I_A.mod <- margin.tensor(I_A.mod,i=2)
      I_S.mod <- margin.tensor(I_S.mod,i=2)
      R.mod <- margin.tensor(R.mod,i=2)
      I_A.mod2 <- margin.tensor(I_A.mod2,i=2)
      I_S.mod2 <- margin.tensor(I_S.mod2,i=2)


    }

    prevalence <- (I_A.mod + I_S.mod + I_A.mod2 + I_S.mod2) / (S.mod + I_A.mod + I_S.mod + R.mod + I_A.mod2 + I_S.mod2) # model prevalence

    # data
    data.time <- data[[k]]

    lppd_k <- fitmodel$dPointObsPrev(data.infecteds = as.vector(data.time[1:(2*nJ_data*n.age.classes)]),
                                     data.agegroupsizes = as.vector(data.time[(2*nJ_data*n.age.classes+1):(4*nJ_data*n.age.classes)]),
                                     model.prev = as.vector(prevalence),
                                     theta = theta, log = F)

    lppd <- c(lppd, lppd_k)

  }

  ###################################################
  # compute density for diagnoses data
  ###################################################

  age.catdata <- c(15,20,25,35,45) # age categories data
  years.in.categories <- to.tensor(NA, dims = list(age=paste0('age', age.classes[1:(n.age.classes)]) ,screencatage=paste0('agedata', age.catdata[-length(age.catdata)])))
  for (k in 1:n.age.classes){
    for (y in 1:(length(age.catdata)-1)){
      years.in.categories[k,y] <- length(intersect( seq(age.classes[k],age.classes[k+1]-1), seq(age.catdata[y],age.catdata[y+1]-1) ))
    }
  }

  norm.years.in.categories <- years.in.categories / margin.tensor(years.in.categories, i=1)
  norm.years.in.categories[is.na(norm.years.in.categories)] <- 0 #possible divisions by zero are changed to zeros

  size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) #JS040416

  for (y in 1:(dim(data[[3]])[3])) { # loop over the years of which data about diagnosis is available

    #model output
    traj.year <- traj[y+1+testperiod_in_burnin, -1]

    tensors <- initvec.to.tensor.cpp(traj.year, par.cond = par.cond)
    D.mod <- tensors$D # number of diagnoses, predicted by the model

    if (nJ_data==1){
      D.mod <- margin.tensor(D.mod,i=2)
    }

    D.mod <- D.mod / size.age.classes # per capita number of diagnoses, predicted by the model
    D.mod <- D.mod * norm.years.in.categories # per capita number of diagnoses/year in age classes model and age classes in data
    D.mod <- margin.tensor(D.mod, i=2) # per capita number of diagnoses in age classes in data (weighted by time spent in data age groups)

    # data
    data.year <- data[[3]][,,y]

    lppd_k <- fitmodel$dPointObsDiag(data.diag = as.vector(data.year),
                                     model.diag = as.vector(D.mod),
                                     theta = theta, log = F)

    lppd <- c(lppd, lppd_k)

  }

  return(lppd=lppd)

}

#calculate the posterior
my_logPosterior <- function(fitmodel, theta, data) {

  trajobs <- my_dTrajObs(fitmodel, theta, data)

  log.prior <- fitmodel$dprior(theta, log=T)
  log.likelihood <- sum(log(trajobs))
  log.posterior <- log.prior+log.likelihood


  return(log.density=log.posterior)
}

#define wrapper
my_logPosterior_epidata <- function(Pars) {
  log.posterior <- my_logPosterior(fitmodel = STI_MCMC_model,
                                  theta = Pars,
                                  data = epidata)
  return(log.density=log.posterior)
}