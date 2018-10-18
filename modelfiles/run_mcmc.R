#############################################################################
### RUN MCMC
#############################################################################

limits.lower <- c(beta = beta.prior.min, 
                  epsilon = epsilon.prior.min, 
                  gamma = gamma.prior.min, 
                  kappa = kappa.prior.min,
                  omega.A = omega.A.prior.min,
                  omega.S = omega.S.prior.min,
                  fsymp.M = fsymp.M.prior.min,
                  fsymp.F = fsymp.F.prior.min,
                  treat = treat.prior.min,
                  eta1 = eta1.prior.min,
                  eta2 = eta2.prior.min,
                  r = r.prior.min)

limits.upper <- c(beta = beta.prior.max, 
                  epsilon = epsilon.prior.max, 
                  gamma = gamma.prior.max, 
                  kappa = kappa.prior.max,
                  omega.A = omega.A.prior.max,
                  omega.S = omega.S.prior.max,
                  fsymp.M = fsymp.M.prior.max,
                  fsymp.F = fsymp.F.prior.max,
                  treat = treat.prior.max,
                  eta1 = eta1.prior.max,
                  eta2 = eta2.prior.max,
                  r = r.prior.max)

# choose sd proposal distribution
my.proposal.sd <- c(beta = 0.02,
                    epsilon = 0.1,
                    gamma = 0.05,
                    kappa = 0.1,
                    omega.A = 0.02,
                    omega.S = 0.02,
                    fsymp.M = 0.02,
                    fsymp.F = 0.02,
                    treat = 0.8,
                    eta1 = 0.7,
                    eta2 = 0.1,
                    r=15)

# set sd of parameters that should not be inferred to zero
# if you want a parameter to be not inferred, choose sd=0. Then parameter retains initial value
my.proposal.sd[!(fullmodel %in% testmodel)] <- 0

my.init.theta2 <- my.init.theta
minv <- limits.lower[fullmodel %in% testmodel]
maxv <- limits.upper[fullmodel %in% testmodel]
set.seed(seed.rnm);randv <- rnorm(rep(1,length(testmodel)), mean=my.init.theta[fullmodel %in% testmodel],sd=my.init.theta[fullmodel %in% testmodel]/100)

my.init.theta2[fullmodel %in% testmodel] <- randv

trace <- my_mcmcMH(target = my_logPosterior_epidata,
                   init.theta = my.init.theta2,
                   proposal.sd = my.proposal.sd,
                   n.iterations =  n.mcmc.iterations,
                   adapt.size.start = NULL,
                   adapt.shape.start = NULL,
                   adapt.size.cooling= NULL,
                   limits = list(lower = limits.lower, upper= limits.upper))

trace$model.type <- testmodel
trace$model.run <- nrun
trace$init.theta <- my.init.theta2

save(trace, file=file.path(path_results,"trace.RData" ))