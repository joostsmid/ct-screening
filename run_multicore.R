# Args from array in bash
args=(commandArgs(TRUE))
args=as.numeric(unlist(args))

# different start values of the parameters
model.run=1:5

# different models (1-4) and sensitivity analyses for these models (5-28)
model.type=rep(list(c("bet_","eps_","kap_","eta1_","eta2_","fsyM_","fsyF_"),
                    c("bet_","eps_","eta1_","eta2_","fsyM_","fsyF_"),
                    c("bet_","eps_","kap_","eta1_","fsyM_","fsyF_"),
                    c("bet_","eps_","eta1_","fsyM_","fsyF_")),7)

# initial conditions, for each model
model.type.inittheta <-
  #mcmc_reinf_bl
  list(c(beta = .85,epsilon=0.8,gamma=365/433,kappa=.7,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = 0,r=200),
       c(beta = 0.64,epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 5.9,eta2 = 0.06,r=200),
       c(beta = .89,epsilon=0.8,gamma=365/433,kappa=.8,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = -2,r=200),
       c(beta = 0.64, epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 3.5,eta2 = -2,r=200),
       
       #mcmc_reinf_gamma90
       c(beta = .85,epsilon=0.8,gamma=0.9*365/433,kappa=.7,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = 0,r=200),
       c(beta = 0.64,epsilon=0.8,gamma=0.9*365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 5.9,eta2 = 0.06,r=200),
       c(beta = .89,epsilon=0.8,gamma=0.9*365/433,kappa=.8,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = -2,r=200),
       c(beta = 0.64, epsilon=0.8,gamma=0.9*365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 3.5,eta2 = -2,r=200),
       #mcmc_reinf_gamma110
       c(beta = .85,epsilon=0.8,gamma=1.1*365/433,kappa=.7,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = 0,r=200),
       c(beta = 0.64,epsilon=0.8,gamma=1.1*365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 5.9,eta2 = 0.06,r=200),
       c(beta = .89,epsilon=0.8,gamma=1.1*365/433,kappa=.8,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = -2,r=200),
       c(beta = 0.64, epsilon=0.8,gamma=1.1*365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 3.5,eta2 = -2,r=200),
       
       #mcmc_reinf_omegar90
       c(beta = .85,epsilon=0.8,gamma=365/433,kappa=.7,omega.A=0.9*0.921*(1-0.194),omega.S=0.9*0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = 0,r=200),
       c(beta = 0.64,epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.9*0.921*(1-0.194),omega.S=0.9*0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 5.9,eta2 = 0.06,r=200),
       c(beta = .89,epsilon=0.8,gamma=365/433,kappa=.8,omega.A=0.9*0.921*(1-0.194),omega.S=0.9*0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = -2,r=200),
       c(beta = 0.64, epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.9*0.921*(1-0.194),omega.S=0.9*0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 3.5,eta2 = -2,r=200),
       #mcmc_reinf_omegar110
       c(beta = .85,epsilon=0.8,gamma=365/433,kappa=.7,omega.A=1.1*0.921*(1-0.194),omega.S=1.1*0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = 0,r=200),
       c(beta = 0.64,epsilon=0.8,gamma=365/433,kappa=0,omega.A=1.1*0.921*(1-0.194),omega.S=1.1*0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 5.9,eta2 = 0.06,r=200),
       c(beta = .89,epsilon=0.8,gamma=365/433,kappa=.8,omega.A=1.1*0.921*(1-0.194),omega.S=1.1*0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = -2,r=200),
       c(beta = 0.64, epsilon=0.8,gamma=365/433,kappa=0,omega.A=1.1*0.921*(1-0.194),omega.S=1.1*0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 3.5,eta2 = -2,r=200),
       
       #mcmc_reinf_chiS90
       c(beta = .85,epsilon=0.8,gamma=365/433,kappa=.7,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=0.9*365/33,eta1 = 2.5,eta2 = 0,r=200),
       c(beta = 0.64,epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=0.9*365/33,eta1 = 5.9,eta2 = 0.06,r=200),
       c(beta = .89,epsilon=0.8,gamma=365/433,kappa=.8,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=0.9*365/33,eta1 = 2.5,eta2 = -2,r=200),
       c(beta = 0.64, epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=0.9*365/33,eta1 = 3.5,eta2 = -2,r=200),
       #mcmc_reinf_chiS110
       c(beta = .85,epsilon=0.8,gamma=365/433,kappa=.7,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=1.1*365/33,eta1 = 2.5,eta2 = 0,r=200),
       c(beta = 0.64,epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=1.1*365/33,eta1 = 5.9,eta2 = 0.06,r=200),
       c(beta = .89,epsilon=0.8,gamma=365/433,kappa=.8,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=1.1*365/33,eta1 = 2.5,eta2 = -2,r=200),
       c(beta = 0.64, epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),notif=0,fsymp.M=0.1,fsymp.F=0.1,treat=1.1*365/33,eta1 = 3.5,eta2 = -2,r=200))

# seed for random number generator per model type, model run
seed.rnm.all <- 1:(length(model.run)*length(model.type))

testmodel <- model.type[[args[1]]]
my.init.theta <- model.type.inittheta[[args[1]]]
nrun <- model.run[args[2]]
seed.rnm <- seed.rnm.all[(args[1]-1)*length(model.run)+args[2]]

path_in <- "modelfiles"
path_data <- "data"
path_results <- "results"

# number of mcmc iterations
n.mcmc.iterations <- 20000

# define basic parameters
sex = c("M","F")
age.classes <- c(15,18,20,25,35,45)
nJ <- 2

# Parameters of full (saturated) model
fullmodel<- c("bet_","eps_","gam_","kap_","omA_","omS_","not_","fsyM_","fsyF_","tre_","eta1_","eta2_","r_")

modelled_NCSP_years <- 11 #2003-2011
burnintime <- 100 #burn period until NCSP
max.simtime <-  burnintime + modelled_NCSP_years

source(file.path(path_in, "define_mcmc.R")) #define model, load mcmc functions and epidata
source(file.path(path_in, "run_mcmc.R")) # run model
