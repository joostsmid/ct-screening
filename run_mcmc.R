# Args from array in bash
args=(commandArgs(TRUE))
args=as.numeric(unlist(args))

# different start values of the parameters
model.run=1:5

# different model types...
model.type=list(c("bet_","eps_","kap_","eta1_","eta2_","fsyM_","fsyF_","gam_","omA_","omS_","tre_"),
                c("bet_","eps_","eta1_","eta2_","fsyM_","fsyF_","gam_","omA_","omS_","tre_"),
                c("bet_","eps_","kap_","eta1_","fsyM_","fsyF_","gam_","omA_","omS_","tre_"),
                c("bet_","eps_","eta1_","fsyM_","fsyF_","gam_","omA_","omS_","tre_"))

# ...with different initial values to start the mcmc chains with (per mcmc chain we will run a slight variation of these inital values)
model.type.inittheta <-
  #mcmc_reinf_bl
  list(c(beta = .85,epsilon=0.8,gamma=365/433,kappa=.7,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = 0,r=200),
       c(beta = 0.64,epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 5.9,eta2 = 0.06,r=200),
       c(beta = .89,epsilon=0.8,gamma=365/433,kappa=.8,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 2.5,eta2 = -2,r=200),
       c(beta = 0.64, epsilon=0.8,gamma=365/433,kappa=0,omega.A=0.921*(1-0.194),omega.S=0.921*(1-0.194),fsymp.M=0.1,fsymp.F=0.1,treat=365/33,eta1 = 3.5,eta2 = -2,r=200))

# data file (maximum/minimum/midpoint data, as in Chandra et al, Eurosurveillance 2016)
testdatafile <- c("testdata_mean.txt","testdata_min.txt","testdata_max.txt")[args[1]]
diagdatafile <- c("diagnosisdata_mean.txt","diagnosisdata_min.txt","diagnosisdata_max.txt")[args[1]]

# seed for random number generator per model type, model run
seed.rnm.all <- 1:(length(model.run)*length(model.type)*3)

# run one model.type..
testmodel <- model.type[[args[2]]]
my.init.theta <- model.type.inittheta[[args[2]]]
# ..starting with a slight variation of initial value to start the mcmc chains with
nrun <- model.run[args[3]]
seed.rnm <- seed.rnm.all[(args[1]-1)*length(model.run)*4+(args[2]-1)*length(model.run)+args[3]]

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
fullmodel<- c("bet_","eps_","gam_","kap_","omA_","omS_","fsyM_","fsyF_","tre_","eta1_","eta2_","r_")

source(file.path(path_in, "define_mcmc.R")) #define model, load mcmc functions and epidata
source(file.path(path_in, "run_mcmc.R")) # run model
