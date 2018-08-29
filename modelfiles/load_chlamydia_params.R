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
### screening rate
##############################################################################################################################

source(file.path(path_in,"CT_screendata_averages.R"))

# per capita number of tests in the entire population
chi_all <- get.CT.screendata(age.classes)

# account for the different size of the age groups: number of tests done per modelled compartment (marginalized over activity classes)
size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) 
chi_all <- chi_all * size.age.classes
chi_all <- chi_all[,,1:(modelled_NCSP_years+1)]

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
### fraction of the total Ct treatments (per year), that are treatments because of partner notification.
##############################################################################################################################

notif <- 0

##############################################################################################################################
### fraction symptomatic
##############################################################################################################################

fsymp  <- to.tensor(c(0.1,0.1), dims = list(sex=c("M","F")) )

##############################################################################################################################
### treatment rate (symptomatic cases)
##############################################################################################################################

treat <- 365/33 # Davies et al, 2014

##############################################################################################################################
### years of testing before years_screening (intervention is assumed to occur at t=burnintime)
### the testing rate is assumed to increase linearly in this period from zero to chi(0,g,j,a)
##############################################################################################################################

testperiod_in_burnin <- 10

##############################################################################################################################
### initials
##############################################################################################################################

source(file.path(path_in,"CT_diagdata_averages.R"))
epidata_diag <- get.CT.diagdata(age.classes)

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

U <- U*size.age.classes
S <- S*size.age.classes
I_A <- I_A*size.age.classes
I_S <- I_S*size.age.classes
R <- R*size.age.classes
I_A2 <- I_A2*size.age.classes
I_S2 <- I_S2*size.age.classes
D <- D*size.age.classes

init.cpp <- c(tens.to.vec.cpp(U), tens.to.vec.cpp(S), tens.to.vec.cpp(I_A), tens.to.vec.cpp(I_S),tens.to.vec.cpp(R),
              tens.to.vec.cpp(I_A2), tens.to.vec.cpp(I_S2), tens.to.vec.cpp(D))

parameters <- list(sex = c("M","F"), nJ = nJ, n.age.classes = n.age.classes, age.classes = age.classes, ar = ar, vl = vl, fold = fold,
                   beta = beta, ct = ct, rho_age = rho_age, theta = theta,
                   epsilon = epsilon, sw = sw, gamma = gamma,  kappa = kappa,
                   chi_all=chi_all, eta=eta, omega=omega, notif=notif, fsymp=fsymp, treat=treat,
                   burnintime = burnintime, testperiod_in_burnin = testperiod_in_burnin)
