# The set of important parameters for the running of the SIR model genericSTImod is defined here. 
# The script uses other scripts to load information from databases, and to prepare the information 
# to be used in auxiliary functions and the SIR model. 
# Here we specify the parameters for Chlamydia modeling. 
#
###############################################################################

n.age.classes <- length(age.classes)-1
age.class.groups <- age.classes[1:n.age.classes]

##############################################################################################################################
### probability of still being virgin at age a, depending on age and gender
##############################################################################################################################


load(file.path(path_data,"vl.Rda"))

##############################################################################################################################
### acquisition rate of new sexual partners per year, activity class and age
##############################################################################################################################

load(file.path(path_data,"ct.Rda")) # partner change rates
load(file.path(path_data,"fold.Rda")) # proportion in low and high activity class

##############################################################################################################################
### age related sexual contact preference matrix
##############################################################################################################################

load(file.path(path_data,"mm.Rda"))
load(file.path(path_data,"mf.Rda"))

rho_age <- to.tensor(rep(0, n.age.classes^2*2), dims = list(sex=c("M","F"), age=paste0('age', age.class.groups), ageprime=paste0('age', age.class.groups) ))

### fill tensor up
for(a in 1:n.age.classes){
  rho_age["F", a, ] <- mf[a,]
  rho_age["M", a, ] <- mm[a,]
}

##############################################################################################################################
### aging rate
##############################################################################################################################

ar <- 1/c(diff(age.classes))

##############################################################################################################################
### transmission probability per infected-susceptible partnership
##############################################################################################################################

beta <- to.tensor(rep(0, 2), dims = list(sex=c("M","F")))
beta["M"] <- 0.7
beta["F"] <- 0.7

##############################################################################################################################
### adjustment for contact rates among genders
##############################################################################################################################

theta <- to.tensor(c(0.5,0.5), dims = list(sex=c("M","F")))

##############################################################################################################################
### assortativity measure (mixing coefficient by class of sexual activity)
##############################################################################################################################

epsilon <- 0.9

##############################################################################################################################
### shifting rate of activitity classes (rate by gender, of switching from one sexual activity classe to another)
##############################################################################################################################

### set default to zero
sw <- to.tensor(rep(0,nJ^2*n.age.classes*2),
                dims = list(sex=c("M","F"), j=paste0('j',1:nJ),jprime=paste0('j',1:nJ), age=paste0('age', age.class.groups)))

for(g in 1:2){
  for(a in 1:n.age.classes){
    switchmat <- matrix(rep(fold[g,],nJ),nJ,nJ)
    diag(switchmat) <- 0
    sw[g,,,a] <- switchmat
  }
}

##############################################################################################################################
### clearance rate of infection
##############################################################################################################################

gamma <- to.tensor(rep(0, 2), dims = list(sex=c("M","F")))
gamma["M"] <- 1
gamma["F"] <- 1

##############################################################################################################################
### factor by which susceptibility to reinfection is reduced
##############################################################################################################################

kappa <- to.tensor(rep(0, length(c("M","F"))), dims = list(sex=c("M","F")) )
kappa["M"] <- .9
kappa["F"] <- .9

##############################################################################################################################
### testing rate
##############################################################################################################################

#also see natsal_screendata_plots.R: there are no screen data available in NATSAL 2. We assume that at that time, nobody was
#screened (prescreen.prop) and that screening was gradually put in place over a period of 9 years to arrive at the screening uptake seen in NATSAL3.
#Here we retrieve the data from NATSAL 3 and distinguish between gender, acticvity class and ages
#it is assumed that after the last year of screening for which data are available, screening rates stay at that level

### tests epidata
testdata <- read.table(file.path(path_data,testdatafile), header=T, row.names=1)
testspp <- to.tensor(NA, dims = list(sex=sex, age=paste0('age', age.class.groups),screening_y = as.numeric(rownames(testdata))))

for (y in 1:dim(testdata)[1]){
  testspp[1,1,y] <- testdata[y,1]/100
  testspp[1,2,y] <- testdata[y,1]/100
  testspp[1,3,y] <- testdata[y,2]/100
  testspp[1,4,y] <- testdata[y,3]/100
  testspp[1,5,y] <- testdata[y,4]/100
  
  testspp[2,1,y] <- testdata[y,5]/100
  testspp[2,2,y] <- testdata[y,5]/100
  testspp[2,3,y] <- testdata[y,6]/100
  testspp[2,4,y] <- testdata[y,7]/100
  testspp[2,5,y] <- testdata[y,8]/100
}

# account for the different size of the age groups: number of tests done per modelled compartment (marginalized over activity classes)
size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) #JS040416
chi_all <- testspp * size.age.classes

##############################################################################################################################
### average factor by which infected/exposed people receive more testing pp than susceptibles/recovereds  (exponential relation)
##############################################################################################################################

eta <- c(3,0)

##############################################################################################################################
### efficacy of treatment (asymptomatic and symptomatic infections), estimable
##############################################################################################################################

# 0.921 Batteiger, JID 2010:201
# This term also accounts for reinfection (resulting in lower treatment efficacy: if treated they don't become susceptible but
# stay in the infected class. We assume that they stay in the same infected class (asympt/sympt) because they are infected by
# the same Ct subtype.

omega <- c(0.921*(1-0.194), 0.921*(1-0.194)) # for asymptomatic and symptomatic infections, respectively


##############################################################################################################################
### fraction symptomatic
##############################################################################################################################

fsymp  <- to.tensor(c(0.1,0.1), dims = list(sex=c("M","F")) )

##############################################################################################################################
### treatment rate (symptomatic cases)
##############################################################################################################################

treat <- 365/33 # Davies et al, 2014

##############################################################################################################################
### time to run system into steady state, and period in which screening is assumed to have increased leinearly before 2000
##############################################################################################################################

burnintime <- 100
testperiod_in_burnin <- 10

##############################################################################################################################
### initials
##############################################################################################################################

size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) #JS040416

U <- to.tensor(rep(0,n.age.classes*2), dims = list(sex=c("M","F"), age=paste0('age', age.classes[1:n.age.classes])) )
U["M",] <- vl[1,1:(length(age.classes)-1)]
U["F",] <- vl[2,1:(length(age.classes)-1)]

initprev <- to.tensor(0.05, dims = list(j=paste0('j',1:nJ))) # initial prevalence of 5% assumed

S <- reorder.tensor( (1-U)*fold*(1-initprev), c("sex","j","age"))
I_A <- reorder.tensor( (1-U)*fold*(1-fsymp)*initprev, c("sex","j","age"))
I_S <- reorder.tensor( (1-U)*fold*fsymp*initprev, c("sex","j","age"))
R <- to.tensor(rep(0, nJ*n.age.classes*length(c("M","F"))), dims  = list(sex=c("M","F"), j=paste0('j',1:nJ), age=paste0('age', age.class.groups)))
I_A2 <- to.tensor(rep(0, nJ*n.age.classes*length(c("M","F"))), dims  = list(sex=c("M","F"), j=paste0('j',1:nJ), age=paste0('age', age.class.groups)))
I_S2 <- to.tensor(rep(0, nJ*n.age.classes*length(c("M","F"))), dims  = list(sex=c("M","F"), j=paste0('j',1:nJ), age=paste0('age', age.class.groups)))
D <- to.tensor(rep(0, nJ*n.age.classes*length(c("M","F"))), dims  = list(sex=c("M","F"), j=paste0('j',1:nJ), age=paste0('age', age.class.groups)))
Inc <- to.tensor(rep(0, nJ*n.age.classes*length(c("M","F"))), dims  = list(sex=c("M","F"), j=paste0('j',1:nJ), age=paste0('age', age.class.groups)))


U <- U*size.age.classes
S <- S*size.age.classes
I_A <- I_A*size.age.classes
I_S <- I_S*size.age.classes
R <- R*size.age.classes
I_A2 <- I_A2*size.age.classes
I_S2 <- I_S2*size.age.classes
D <- D*size.age.classes

init.cpp <- c(tens.to.vec.cpp(U), tens.to.vec.cpp(S), tens.to.vec.cpp(I_A), tens.to.vec.cpp(I_S),tens.to.vec.cpp(R),
              tens.to.vec.cpp(I_A2), tens.to.vec.cpp(I_S2), tens.to.vec.cpp(D), tens.to.vec.cpp(Inc))

parameters <- list(sex = c("M","F"), nJ = nJ, n.age.classes = n.age.classes, age.classes = age.classes, ar = ar, vl = vl, fold = fold,
                   beta = beta, ct = ct, rho_age = rho_age, theta = theta,
                   epsilon = epsilon, sw = sw, gamma = gamma,  kappa = kappa,
                   chi_all=chi_all, eta=eta, omega=omega, fsymp=fsymp, treat=treat,
                   burnintime = burnintime, testperiod_in_burnin = testperiod_in_burnin)