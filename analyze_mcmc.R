# Args from array in bash
args=(commandArgs(TRUE))
args=as.numeric(unlist(args))

# mcmcresultsdir <- "C:/Users/smid/temp genericSTI/scripts_final/Rdata files/"###############################################
# mcmcresultsdir <- "../genericSTI/results/scripts_final_Rdata/"

# different start values of the parameters
modelruns=1:5
runspermodel <- tail(modelruns, n=1)

# different model types
model.type=rep(list(c("bet_","eps_","kap_","eta1_","eta2_","fsyM_","fsyF_"),
                    c("bet_","eps_","eta1_","eta2_","fsyM_","fsyF_"),
                    c("bet_","eps_","kap_","eta1_","fsyM_","fsyF_"),
                    c("bet_","eps_","eta1_","fsyM_","fsyF_")),7)

testmodel <- model.type[[args]]
testmodel.number <- args

fullmodel<- c("bet_","eps_","gam_","kap_","omA_","omS_","not_","fsyM_","fsyF_","tre_","eta1_","eta2_","r_")

sex = c("M","F")
age.classes <- c(15,18,20,25,35,45)
nJ <- 2

modelled_NCSP_years <- 11 #2003-2011
burnintime <- 100 #burn period until NCSP
max.simtime <-  burnintime + modelled_NCSP_years

path_in <- "modelfiles"
path_data <- "data"
path_results <- "results"
path_figures <- "figures"

n.mcmc.iterations <- 20000
nsimsplot <- 100

source(file.path(path_in, "define_mcmc.R"))
source(file.path(path_in, "plot_functions.R"))

mycol <- brewer.pal(9,"Set1")
mycol <- c(mycol[9],mycol[1:8])

##############################################################################
# mcmc traces
##############################################################################

col1M <- "blue"
col2M <- "green4"
col1F <- "red"
col2F <- "orange"

mcmclst <- list()

for (nrun in modelruns){
  
  load(file.path(path_results,paste0("trace", (testmodel.number-1)*runspermodel+nrun, ".RData")))
  tracetrace <- as.matrix(trace$trace)
  allvarnames <- colnames(tracetrace)
  tracetrace <- as.data.frame(matrix(unlist(tracetrace), ncol = dim(tracetrace)[2], byrow = FALSE))
  names(tracetrace) <- allvarnames
  
  assign(paste0("tracetrace.params",nrun), tracetrace[,c((1:length(fullmodel))[fullmodel %in% testmodel],dim(tracetrace)[2])])
  assign(paste0("trace",nrun), mcmc(eval(parse(text = paste0("tracetrace.params",nrun)))))
  mcmclst[[nrun]] <- eval(parse(text = paste0("trace",nrun)))
  
}

# combine traces as mcmc.list object
trace.mcmc.params <- mcmc.list(mcmclst)
nall.mcmc.iterations <- dim(trace.mcmc.params[[1]])[1]

#acceptance rate, effective sampling size
accrate <- 1 - rejectionRate(trace.mcmc.params)
effectiveSize <- effectiveSize(trace.mcmc.params)

# remove burnin
trace.mcmc.params.burned <- burnAndThin(trace.mcmc.params, burn = round(nall.mcmc.iterations*0.4))

#do thinning
trace.mcmc.params.burned.thinned.xy <- burnAndThin(trace.mcmc.params, 
                                                   burn = round(nall.mcmc.iterations*0.4), 
                                                   thin = 20)

trace.mcmc.params.burned.thinned <- burnAndThin(trace.mcmc.params, 
                                                burn = round(nall.mcmc.iterations*0.4), 
                                                thin = round(runspermodel*nall.mcmc.iterations*0.6/nsimsplot))

summary.trace.mcmc.params.burned <- summary(trace.mcmc.params.burned) 

pdf(file.path(path_figures,"xyplot.pdf"),width=5,height=7)
xyplot(trace.mcmc.params.burned.thinned.xy)
dev.off()

pdf(file.path(path_figures,"densityplot.pdf"),width=5,height=7)
densityplot(trace.mcmc.params.burned.thinned.xy)
dev.off()

##############################################################################
# compute Gelman-Rubin diagnostic
##############################################################################

sj2 <- 1/(dim(trace.mcmc.params.burned[[1]])[1]-1) 

n.chains <- runspermodel
n.params <- dim(trace.mcmc.params.burned[[1]])[2]
n.it <- dim(trace.mcmc.params.burned[[1]])[1]

Rhat <- rep(NA, n.params-1)
for (k in 1:(n.params-1)){ # loop over number of parameters
  
  sj2 <- rep(NA, n.chains)
  thetabarj <- rep(NA, n.chains)
  for (j in 1:n.chains){ # loop over number of chains
    thetabarj[j] <- mean(trace.mcmc.params.burned[[j]][,k])
    sj2[j] <- 1/(n.it-1) * sum((trace.mcmc.params.burned[[j]][,k]-thetabarj[j] )^2)
  }
  
  thetabarbar <- 1/n.chains * sum(thetabarj)
  W <- 1/n.chains * sum(sj2)
  B <- n.it / (n.chains-1) * sum((thetabarj[j] - thetabarbar)^2)
  
  varthetahat <- (1-1/n.it) * W + 1/n.it * B
  Rhat[k] <- sqrt(varthetahat / W)
  
}

##############################################################################
# COMPUTE DIC
##############################################################################

trace.combined <- ldply(trace.mcmc.params.burned.thinned)

# matrix for log pointwise predictive density
lppd.M <- matrix(NA, nrow=dim(trace.combined)[1], ncol=prod(dim(epidata[[1]])[1:2]) + prod(dim(epidata[[2]])[1:2]) + prod(dim(epidata[[3]])))
for (k in 1:dim(trace.combined)[1]){
  theta.k.allpars <- unlist(c(trace.combined[k,1:(dim(trace.combined)[2]-1)], tracetrace[1,which(!(fullmodel %in% testmodel))]))
  lppd.M[k,] <- my_dTrajObs(STI_MCMC_model, theta.k.allpars , epidata)
}

theta.bar <- colMeans(trace.combined[-dim(trace.combined)[2]])
theta.bar.allpars <- unlist(c(theta.bar, tracetrace[1,which(!(fullmodel %in% testmodel))]))
log.like.theta.bar <- sum(log(my_dTrajObs(STI_MCMC_model, theta.bar.allpars , epidata)))
D.theta.bar <- -2 * log.like.theta.bar
#p.D <- var(-2 * trace.combined$log.density)/2
p.D <- var(-2 * rowSums(log(lppd.M)))/2
DIC <- D.theta.bar + 2 * p.D
# Spiegelhalter 2002: 9.2.4. What is an important difference in DIC?
# Burnham and Anderson (1998) suggested models receiving AIC within 1-2 of the 'best' deserve
# consideration, and 3-7 have considerably less support: these rules of thumb appear to work
# reasonably well for DIC.

##############################################################################
# COMPUTE WAIC
##############################################################################

# https://link.springer.com/content/pdf/10.1007%2Fs11222-013-9416-2.pdf
pW <- 2* sum( log(colMeans(lppd.M))-colMeans(log(lppd.M)))
pW.alt <- sum(apply(log(lppd.M), MARGIN=2, var))

lppd <- sum( log(colMeans(lppd.M)))

WAIC <- -2 * (lppd - pW)

##############################################################################
# COMPARE MCMC RESULT TO PREVALENCE DATA
##############################################################################

myyears <- 2000:2011

# function to create posterior estimates for mean prevalence and numbers of diagnoses pp (plus CI)
simulate_prevs_diag <- function(trace.mcmc.pars, fitmodel, years=myyears, CI=.95){
  # function to plot CT prevalences in certain year
  
  plottedsims <- dim(trace.mcmc.pars)[1]
  times <- c(0, (burnintime - testperiod_in_burnin):max.simtime)
  prevalenceT <- to.tensor(NA, dims = list(sex=c("M","F"), age=paste0('age', age.class.groups), screening_y = years, nsim = 1:plottedsims))
  diagT <- to.tensor(NA, dims = list(sex=c("M","F"), age=paste0('age', age.class.groups), screening_y = years, nsim = 1:plottedsims))
  
  
  for (s in 1:plottedsims){
    theta <- trace.mcmc.pars[s, 1:length(testmodel)]
    theta.allpars <- unlist(c(theta, tracetrace[1,which(!(fullmodel %in% testmodel))]))
    
    simulation <- fitmodel$simulate(theta.allpars, times)
    
    for (y in 1:length(years)){
      
      # model output
      #traj.last <- traj[traj$time==max.simtime, -1]
      simulation.y <- simulation[which(times+2000-burnintime==years[y]), -1]
      
      tensors <- initvec.to.tensor.cpp(simulation.y, par.cond = list(age_classes =age.class.groups, nmb_act_classes = nJ))
      
      U.mod <- tensors$U
      S.mod <- tensors$S
      I_A.mod <- tensors$I_A
      I_S.mod <- tensors$I_S
      R.mod <- tensors$R
      I_A.mod2 <- tensors$I_A2
      I_S.mod2 <- tensors$I_S2
      D.mod <- tensors$D
      
      #marginalize across activity classes
      S.mod <- margin.tensor(S.mod,i=2)
      I_A.mod <- margin.tensor(I_A.mod,i=2)
      I_S.mod <- margin.tensor(I_S.mod,i=2)
      R.mod <- margin.tensor(R.mod,i=2)
      I_A.mod2 <- margin.tensor(I_A.mod2,i=2)
      I_S.mod2 <- margin.tensor(I_S.mod2,i=2)
      D.mod <- margin.tensor(D.mod,i=2)
      
      prevalence <- (I_A.mod + I_S.mod + I_A.mod2 + I_S.mod2) / (S.mod + I_A.mod + I_S.mod + R.mod + I_A.mod2 + I_S.mod2) # model prevalence
      
      prevalenceT[,,y,s] <- prevalence
      diagT[,,y,s] <- D.mod
      
    }
    
  }
  
  sorted_prevalenceT <- prevalenceT
  sorted_diagT <- diagT
  for (g in 1:dim(sorted_prevalenceT)[1]){
    for (a in 1:dim(sorted_prevalenceT)[2]){
      for (y in 1:dim(sorted_prevalenceT)[3]){
        sorted_prevalenceT[g,a,y,] <- sort(prevalenceT[g,a,y,])
        sorted_diagT[g,a,y,] <- sort(diagT[g,a,y,])
      }
    }
  }
  
  mean_prevalenceT <- apply(prevalenceT, MARGIN=c(1,2,3), FUN=mean)
  low_prevalenceT <- sorted_prevalenceT[,,,ceiling(((1-CI)/2)*dim(trace.mcmc.pars)[1])]
  upp_prevalenceT <- sorted_prevalenceT[,,,floor((1-(1-CI)/2)*dim(trace.mcmc.pars)[1])]
  
  mean_diagT <- apply(diagT, MARGIN=c(1,2,3), FUN=mean)
  low_diagT <- sorted_diagT[,,,ceiling(((1-CI)/2)*dim(trace.mcmc.pars)[1])]
  upp_diagT <- sorted_diagT[,,,floor((1-(1-CI)/2)*dim(trace.mcmc.pars)[1])]
  
  return(list(prevalenceT=prevalenceT, mean_prevalenceT=mean_prevalenceT, low_prevalenceT=low_prevalenceT, upp_prevalenceT=upp_prevalenceT,
              diagT=diagT, mean_diagT=mean_diagT, low_diagT=low_diagT, upp_diagT=upp_diagT))
  
}

postsims <- simulate_prevs_diag(trace.combined, STI_MCMC_model, years=2000:2011)

### (comparison 2000-2011)
plotdataM.2000 <- data.frame(ages=age.class.groups,
                             agesto=age.classes[2:(n.age.classes+1)],
                             prev.mean=summary(mcmc(data.matrix(t(postsims$prevalenceT[1,,1,]))))$statistics[,"Mean"],
                             prev.low=summary(mcmc(data.matrix(t(postsims$prevalenceT[1,,1,]))))$quantiles[,"2.5%"],
                             prev.high=summary(mcmc(data.matrix(t(postsims$prevalenceT[1,,1,]))))$quantiles[,"97.5%"],
                             prevdata.mean=binom.confint(as.numeric(epidata$time_0[1,,1]),as.numeric(epidata$time_0[1,,2]), conf.level=.95,method="wilson")$mean,
                             prevdata.low=binom.confint(as.numeric(epidata$time_0[1,,1]),as.numeric(epidata$time_0[1,,2]), conf.level=.95,method="wilson")$lower,
                             prevdata.upp=binom.confint(as.numeric(epidata$time_0[1,,1]),as.numeric(epidata$time_0[1,,2]), conf.level=.95,method="wilson")$upper)

plotdataF.2000 <- data.frame(ages=age.class.groups,
                             agesto=age.classes[2:(n.age.classes+1)],
                             prev.mean=summary(mcmc(data.matrix(t(postsims$prevalenceT[2,,1,]))))$statistics[,"Mean"],
                             prev.low=summary(mcmc(data.matrix(t(postsims$prevalenceT[2,,1,]))))$quantiles[,"2.5%"],
                             prev.high=summary(mcmc(data.matrix(t(postsims$prevalenceT[2,,1,]))))$quantiles[,"97.5%"],
                             prevdata.mean=binom.confint(as.numeric(epidata$time_0[2,,1]),as.numeric(epidata$time_0[2,,2]), conf.level=.95,method="wilson")$mean,
                             prevdata.low=binom.confint(as.numeric(epidata$time_0[2,,1]),as.numeric(epidata$time_0[2,,2]), conf.level=.95,method="wilson")$lower,
                             prevdata.upp=binom.confint(as.numeric(epidata$time_0[2,,1]),as.numeric(epidata$time_0[2,,2]), conf.level=.95,method="wilson")$upper)

plotdataM.2011 <- data.frame(ages=age.class.groups,
                             agesto=age.classes[2:(n.age.classes+1)],
                             prev.mean=summary(mcmc(data.matrix(t(postsims$prevalenceT[1,,12,]))))$statistics[,"Mean"],
                             prev.low=summary(mcmc(data.matrix(t(postsims$prevalenceT[1,,12,]))))$quantiles[,"2.5%"],
                             prev.high=summary(mcmc(data.matrix(t(postsims$prevalenceT[1,,12,]))))$quantiles[,"97.5%"],
                             prevdata.mean=binom.confint(as.numeric(epidata$time_11[1,,1]),as.numeric(epidata$time_11[1,,2]), conf.level=.95,method="wilson")$mean,
                             prevdata.low=binom.confint(as.numeric(epidata$time_11[1,,1]),as.numeric(epidata$time_11[1,,2]), conf.level=.95,method="wilson")$lower,
                             prevdata.upp=binom.confint(as.numeric(epidata$time_11[1,,1]),as.numeric(epidata$time_11[1,,2]), conf.level=.95,method="wilson")$upper)

plotdataF.2011 <- data.frame(ages=age.class.groups,
                             agesto=age.classes[2:(n.age.classes+1)],
                             prev.mean=summary(mcmc(data.matrix(t(postsims$prevalenceT[2,,12,]))))$statistics[,"Mean"],
                             prev.low=summary(mcmc(data.matrix(t(postsims$prevalenceT[2,,12,]))))$quantiles[,"2.5%"],
                             prev.high=summary(mcmc(data.matrix(t(postsims$prevalenceT[2,,12,]))))$quantiles[,"97.5%"],
                             prevdata.mean=binom.confint(as.numeric(epidata$time_11[2,,1]),as.numeric(epidata$time_11[2,,2]), conf.level=.95,method="wilson")$mean,
                             prevdata.low=binom.confint(as.numeric(epidata$time_11[2,,1]),as.numeric(epidata$time_11[2,,2]), conf.level=.95,method="wilson")$lower,
                             prevdata.upp=binom.confint(as.numeric(epidata$time_11[2,,1]),as.numeric(epidata$time_11[2,,2]), conf.level=.95,method="wilson")$upper)

plotprevs.M.2000 <- ggplot() +
  geom_segment(data = plotdataM.2000, aes(x = ages, xend = agesto, y = prev.mean, yend = prev.mean), size = 1, lineend = "butt") +
  geom_rect(data = plotdataM.2000, aes(xmin=ages, xmax=agesto, ymin=prev.low, ymax=prev.high),alpha=.3) + 
  geom_point(data=plotdataM.2000, mapping=aes(x=(ages+agesto)/2, y=prevdata.mean), size=3, shape=19, alpha = 1) +
  geom_errorbar(data=plotdataM.2000, mapping=aes(x=(ages+agesto)/2, ymin=prevdata.low, ymax=prevdata.upp), width=0.2, size=.1, alpha = 1) + 
  scale_y_continuous(breaks = seq(0,0.07, by=0.01)) +
  labs(x = "Ages") +
  labs(y = "Prevalence") +
  annotate("text", x=40, y=0.069, label= "Men 2000", size=4) +
  theme_bw(base_family = "Times") +
  coord_cartesian(ylim=c(0,0.07)) +
  scale_x_continuous(breaks=(age.class.groups+age.classes[-1])/2, labels=paste0(age.class.groups,"-",(age.classes[-1]-1))) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size=8))

plotprevs.F.2000 <- ggplot() +
  geom_segment(data = plotdataF.2000, aes(x = ages, xend = agesto, y = prev.mean, yend = prev.mean), size = 1, lineend = "butt") +
  geom_rect(data = plotdataF.2000, aes(xmin=ages, xmax=agesto, ymin=prev.low, ymax=prev.high), alpha=.3) + 
  geom_point(data=plotdataF.2000, mapping=aes(x=(ages+agesto)/2, y=prevdata.mean), size=3, shape=19, alpha = 1) +
  geom_errorbar(data=plotdataF.2000, mapping=aes(x=(ages+agesto)/2, ymin=prevdata.low, ymax=prevdata.upp), width=0.2, size=.1, alpha = 1) + 
  scale_y_continuous(breaks = seq(0,0.07, by=0.01)) +
  labs(x = "Ages") +
  labs(y = "Prevalence") +
  annotate("text", x=38, y=0.069, label= "Women 2000", size=4) +
  theme_bw(base_family = "Times") +
  coord_cartesian(ylim=c(0,0.07)) +
  scale_x_continuous(breaks=(age.class.groups+age.classes[-1])/2, labels=paste0(age.class.groups,"-",(age.classes[-1]-1))) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size=8))

plotprevs.M.2011 <- ggplot() +
  geom_segment(data = plotdataM.2011, aes(x = ages, xend = agesto, y = prev.mean, yend = prev.mean), size = 1, lineend = "butt") +
  geom_rect(data = plotdataM.2011, aes(xmin=ages, xmax=agesto, ymin=prev.low, ymax=prev.high),alpha=.3) + 
  geom_point(data=plotdataM.2011, mapping=aes(x=(ages+agesto)/2, y=prevdata.mean), size=3, shape=19, alpha = 1) +
  geom_errorbar(data=plotdataM.2011, mapping=aes(x=(ages+agesto)/2, ymin=prevdata.low, ymax=prevdata.upp), width=0.2, size=.1, alpha = 1) + 
  scale_y_continuous(breaks = seq(0,0.07, by=0.01)) +
  labs(x = "Ages") +
  labs(y = "Prevalence") +
  annotate("text", x=40, y=0.069, label= "Men 2011", size=4) +
  theme_bw(base_family = "Times") +
  coord_cartesian(ylim=c(0,0.07)) +
  scale_x_continuous(breaks=(age.class.groups+age.classes[-1])/2, labels=paste0(age.class.groups,"-",(age.classes[-1]-1))) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size=8))

plotprevs.F.2011 <- ggplot() +
  geom_segment(data = plotdataF.2011, aes(x = ages, xend = agesto, y = prev.mean, yend = prev.mean), size = 1, lineend = "butt") +
  geom_rect(data = plotdataF.2011, aes(xmin=ages, xmax=agesto, ymin=prev.low, ymax=prev.high), alpha=.3) + 
  geom_point(data=plotdataF.2011, mapping=aes(x=(ages+agesto)/2, y=prevdata.mean), size=3, shape=19, alpha = 1) +
  geom_errorbar(data=plotdataF.2011, mapping=aes(x=(ages+agesto)/2, ymin=prevdata.low, ymax=prevdata.upp), width=0.2, size=.1, alpha = 1) + 
  scale_y_continuous(breaks = seq(0,0.07, by=0.01)) +
  labs(x = "Ages") +
  labs(y = "Prevalence") +
  annotate("text", x=38, y=0.069, label= "Women 2011", size=4) +
  theme_bw(base_family = "Times") +
  coord_cartesian(ylim=c(0,0.07)) +
  scale_x_continuous(breaks=(age.class.groups+age.classes[-1])/2, labels=paste0(age.class.groups,"-",(age.classes[-1]-1))) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size=8))

pdf(file.path(path_figures,"epiprevplot.pdf") ,width=7,height=7)
multiplot(plotlist=list(plotprevs.M.2000,plotprevs.F.2000,
                        plotprevs.M.2011,plotprevs.F.2011), layout=matrix(c(1,2,3,4), nrow=2, byrow=TRUE))
dev.off()

# change of prevalence in time for different sex and age groups
gendernames <- c("M","F")
df_prev <- data.frame(years=myyears)

for (g in 1:2){
  
  df_prev[paste0("mean_",gendernames[g],"_1519")] <- as.vector(3/5*postsims$mean_prevalenceT[g,1,] + 2/5*postsims$mean_prevalenceT[g,2,])
  df_prev[paste0("mean_",gendernames[g],"_2024")] <- as.vector(postsims$mean_prevalenceT[g,3,])
  df_prev[paste0("mean_",gendernames[g],"_2534")] <- as.vector(postsims$mean_prevalenceT[g,4,])
  df_prev[paste0("mean_",gendernames[g],"_3544")] <- as.vector(postsims$mean_prevalenceT[g,5,])
  
  df_prev[paste0("low_",gendernames[g],"_1519")] <- as.vector(3/5*postsims$low_prevalenceT[g,1,] + 2/5*postsims$low_prevalenceT[g,2,])
  df_prev[paste0("low_",gendernames[g],"_2024")] <- as.vector(postsims$low_prevalenceT[g,3,])
  df_prev[paste0("low_",gendernames[g],"_2534")] <- as.vector(postsims$low_prevalenceT[g,4,])
  df_prev[paste0("low_",gendernames[g],"_3544")] <- as.vector(postsims$low_prevalenceT[g,5,])
  
  df_prev[paste0("upp_",gendernames[g],"_1519")] <- as.vector(3/5*postsims$upp_prevalenceT[g,1,] + 2/5*postsims$upp_prevalenceT[g,2,])
  df_prev[paste0("upp_",gendernames[g],"_2024")] <- as.vector(postsims$upp_prevalenceT[g,3,])
  df_prev[paste0("upp_",gendernames[g],"_2534")] <- as.vector(postsims$upp_prevalenceT[g,4,])
  df_prev[paste0("upp_",gendernames[g],"_3544")] <- as.vector(postsims$upp_prevalenceT[g,5,])
  
}

prevtimeplot_M <- ggplot() +
  geom_line(data=df_prev, mapping=aes(x=years, y=mean_M_1519), colour=mycol[1]) +
  geom_ribbon(data=df_prev, aes(x=years, ymin=low_M_1519, ymax=upp_M_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_prev, mapping=aes(x=years, y=mean_M_2024), colour=mycol[2]) +
  geom_ribbon(data=df_prev, aes(x=years, ymin=low_M_2024, ymax=upp_M_2024), fill=mycol[2], alpha = .3) +
  geom_line(data=df_prev, mapping=aes(x=years, y=mean_M_2534), colour=mycol[3]) +
  geom_ribbon(data=df_prev, aes(x=years, ymin=low_M_2534, ymax=upp_M_2534), fill=mycol[3], alpha = .3) +
  geom_line(data=df_prev, mapping=aes(x=years, y=mean_M_3544), colour=mycol[4]) +
  geom_ribbon(data=df_prev, aes(x=years, ymin=low_M_3544, ymax=upp_M_3544), fill=mycol[4], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_prev$years), max(df_prev$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,0.06, by=0.01)) +
  coord_cartesian(ylim=c(0,0.06)) +
  labs(x = "Years") +
  labs(y = "Prevalence") +
  ggtitle("Men") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), text=element_text(size=8))


prevtimeplot_F <- ggplot() +
  geom_line(data=df_prev, mapping=aes(x=years, y=mean_F_1519), colour=mycol[1]) +
  geom_ribbon(data=df_prev, aes(x=years, ymin=low_F_1519, ymax=upp_F_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_prev, mapping=aes(x=years, y=mean_F_2024), colour=mycol[2]) +
  geom_ribbon(data=df_prev, aes(x=years, ymin=low_F_2024, ymax=upp_F_2024), fill=mycol[2], alpha = .3) +
  geom_line(data=df_prev, mapping=aes(x=years, y=mean_F_2534), colour=mycol[3]) +
  geom_ribbon(data=df_prev, aes(x=years, ymin=low_F_2534, ymax=upp_F_2534), fill=mycol[3], alpha = .3) +
  geom_line(data=df_prev, mapping=aes(x=years, y=mean_F_3544), colour=mycol[4]) +
  geom_ribbon(data=df_prev, aes(x=years, ymin=low_F_3544, ymax=upp_F_3544), fill=mycol[4], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_prev$years), max(df_prev$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,0.06, by=0.01)) +
  coord_cartesian(ylim=c(0,0.06)) +
  labs(x = "Years") +
  labs(y = "Prevalence") +
  ggtitle("Women") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

pdf(file.path(path_figures,"prevtimeplot.pdf") ,width=7,height=4)
multiplot(plotlist=list(prevtimeplot_M,prevtimeplot_F), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

####################################################################
# diagnosis rates and diagnosis data
####################################################################

age.catdata <- c(15,20,25,35,45)
years.in.categories <- to.tensor(NA, dims = list(age=paste0('age', age.classes[1:(n.age.classes)]) ,screencatage=paste0('diagcatage', age.catdata[1:(length(age.catdata)-1)])))
for (k in 1:n.age.classes){
  for (y in 1:(length(age.catdata)-1)){
    years.in.categories[k,y] <- length(intersect( seq(age.classes[k],age.classes[k+1]-1), seq(age.catdata[y],age.catdata[y+1]-1) ))
  }
}
norm.years.in.categories <- years.in.categories / margin.tensor(years.in.categories, i=1) # fraction of time that individuals in age.classes[k] are in the screencatage[j] class

# test data
source(file.path(path_in,"CT_screendata.R"))
uptake.minmax <- get.CT.screendata(age.catdata)
names_uptake <- dimnames(uptake.minmax)
names_uptake$screening_y <- names_uptake$screening_y[1:12]
uptake.minmax <- uptake.minmax[,,,1:12]
dimnames(uptake.minmax) <- names_uptake

# diagnoses
source(file.path(path_in,"CT_diagdata.R"))
epidata_diag <- get.CT.diagdata(age.catdata)
names_diag <- dimnames(epidata_diag)
names_diag$screening_y <- names_diag$screening_y[1:12]
epidata_diag <- epidata_diag[,,,1:12]
dimnames(epidata_diag) <- names_diag

mean_diagT <- to.tensor(as.vector(postsims$mean_diagT), dims= list(sex=c("M","F"), age=paste0('age', age.classes[1:(n.age.classes)]), screening_y=myyears ))
size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) #JS040416
mean_diagT <- mean_diagT / size.age.classes * norm.years.in.categories
mean_diagT <- margin.tensor(mean_diagT, i=2)

low_diagT <- to.tensor(as.vector(postsims$low_diagT), dims= list(sex=c("M","F"), age=paste0('age', age.classes[1:(n.age.classes)]), screening_y=myyears ))
size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) #JS040416
low_diagT <- low_diagT / size.age.classes * norm.years.in.categories
low_diagT <- margin.tensor(low_diagT, i=2)

upp_diagT <- to.tensor(as.vector(postsims$upp_diagT), dims= list(sex=c("M","F"), age=paste0('age', age.classes[1:(n.age.classes)]), screening_y=myyears ))
size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) #JS040416
upp_diagT <- upp_diagT / size.age.classes * norm.years.in.categories
upp_diagT <- margin.tensor(upp_diagT, i=2)

gencode<-1
diagdataM <- data.frame(years=rep(myyears,4),
                        shift=c(rep(-.15,length(myyears)),
                                rep(-.05,length(myyears)),
                                rep(.05,length(myyears)),
                                rep(.15,length(myyears))),
                        testdata.upt.min=c(as.vector(uptake.minmax[gencode,1,1,]),
                                           as.vector(uptake.minmax[gencode,2,1,]),
                                           as.vector(uptake.minmax[gencode,3,1,]),
                                           as.vector(uptake.minmax[gencode,4,1,])),
                        testdata.upt.max=c(as.vector(uptake.minmax[gencode,1,2,]),
                                           as.vector(uptake.minmax[gencode,2,2,]),
                                           as.vector(uptake.minmax[gencode,3,2,]),
                                           as.vector(uptake.minmax[gencode,4,2,])),
                        
                        diagdata.upt.min=c(as.vector(epidata_diag[gencode,1,1,])/100000,
                                           as.vector(epidata_diag[gencode,2,1,])/100000,
                                           as.vector(epidata_diag[gencode,3,1,])/100000,
                                           as.vector(epidata_diag[gencode,4,1,])/100000),
                        diagdata.upt.max=c(as.vector(epidata_diag[gencode,1,2,])/100000,
                                           as.vector(epidata_diag[gencode,2,2,])/100000,
                                           as.vector(epidata_diag[gencode,3,2,])/100000,
                                           as.vector(epidata_diag[gencode,4,2,])/100000),
                        
                        diagmodel.upt.mean=c(as.vector(mean_diagT[gencode,,1]),
                                             as.vector(mean_diagT[gencode,,2]),
                                             as.vector(mean_diagT[gencode,,3]),
                                             as.vector(mean_diagT[gencode,,4])),
                        
                        diagmodel.upt.low=c(as.vector(low_diagT[gencode,,1]),
                                            as.vector(low_diagT[gencode,,2]),
                                            as.vector(low_diagT[gencode,,3]),
                                            as.vector(low_diagT[gencode,,4])),
                        
                        diagmodel.upt.upp=c(as.vector(upp_diagT[gencode,,1]),
                                            as.vector(upp_diagT[gencode,,2]),
                                            as.vector(upp_diagT[gencode,,3]),
                                            as.vector(upp_diagT[gencode,,4])),
                        
                        grp=c(rep("15-19",length(myyears)),
                              rep("20-24",length(myyears)),
                              rep("25-34",length(myyears)),
                              rep("35-44",length(myyears))))

gencode<-2
diagdataF <- data.frame(years=rep(myyears,4),
                        shift=c(rep(-.15,length(myyears)),
                                rep(-.05,length(myyears)),
                                rep(.05,length(myyears)),
                                rep(.15,length(myyears))),
                        testdata.upt.min=c(as.vector(uptake.minmax[gencode,1,1,]),
                                           as.vector(uptake.minmax[gencode,2,1,]),
                                           as.vector(uptake.minmax[gencode,3,1,]),
                                           as.vector(uptake.minmax[gencode,4,1,])),
                        testdata.upt.max=c(as.vector(uptake.minmax[gencode,1,2,]),
                                           as.vector(uptake.minmax[gencode,2,2,]),
                                           as.vector(uptake.minmax[gencode,3,2,]),
                                           as.vector(uptake.minmax[gencode,4,2,])),
                        
                        diagdata.upt.min=c(as.vector(epidata_diag[gencode,1,1,])/100000,
                                           as.vector(epidata_diag[gencode,2,1,])/100000,
                                           as.vector(epidata_diag[gencode,3,1,])/100000,
                                           as.vector(epidata_diag[gencode,4,1,])/100000),
                        diagdata.upt.max=c(as.vector(epidata_diag[gencode,1,2,])/100000,
                                           as.vector(epidata_diag[gencode,2,2,])/100000,
                                           as.vector(epidata_diag[gencode,3,2,])/100000,
                                           as.vector(epidata_diag[gencode,4,2,])/100000),
                        
                        diagmodel.upt.mean=c(as.vector(mean_diagT[gencode,,1]),
                                             as.vector(mean_diagT[gencode,,2]),
                                             as.vector(mean_diagT[gencode,,3]),
                                             as.vector(mean_diagT[gencode,,4])),
                        
                        diagmodel.upt.low=c(as.vector(low_diagT[gencode,,1]),
                                            as.vector(low_diagT[gencode,,2]),
                                            as.vector(low_diagT[gencode,,3]),
                                            as.vector(low_diagT[gencode,,4])),
                        
                        diagmodel.upt.upp=c(as.vector(upp_diagT[gencode,,1]),
                                            as.vector(upp_diagT[gencode,,2]),
                                            as.vector(upp_diagT[gencode,,3]),
                                            as.vector(upp_diagT[gencode,,4])),
                        
                        grp=c(rep("15-19",length(myyears)),
                              rep("20-24",length(myyears)),
                              rep("25-34",length(myyears)),
                              rep("35-44",length(myyears))))

mycol <- brewer.pal(9,"Set1")
mycol <- c(mycol[9],mycol[1:8])

# for women from 2000: only point estimates of data available. No min and max. For plotting: adapt this:
diagdataF$testdata.upt.max[diagdataF$testdata.upt.min==diagdataF$testdata.upt.max] <- 
  diagdataF$testdata.upt.max[diagdataF$testdata.upt.min==diagdataF$testdata.upt.max] + .005
diagdataF$diagdata.upt.max[which(diagdataF$diagdata.upt.min==diagdataF$diagdata.upt.max)] <- 
  diagdataF$diagdata.upt.max[which(diagdataF$diagdata.upt.min==diagdataF$diagdata.upt.max)] + .0005



testplotM <- ggplot() +
  geom_segment(data = diagdataM, aes(x = years+shift, xend = years+shift, y = testdata.upt.min, yend = testdata.upt.max, col=grp), size = 1, lineend = "butt") +
  #geom_point(data=diagdataM, aes(x=years+shift, y=(testdata.upt.min+testdata.upt.max)/2, col=grp), size=3, shape=19, alpha = 1) +
  scale_x_continuous(breaks = seq(min(diagdataM$years), max(diagdataM$years), by = 2)) +
  scale_color_manual("", values = mycol[1:4]) +
  labs(x = "Years") +
  labs(y = "Yearly number of tests pp") +
  scale_y_continuous(breaks = seq(0,0.4, by=0.05)) +
  coord_cartesian(ylim=c(0,0.4)) +
  ggtitle("Men") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.justification=c(0.01,.99), legend.position=c(0.01,.99),text=element_text(size=8))


diagplotM <- ggplot() +
  geom_segment(data = diagdataM, aes(x = years+shift, xend = years+shift, y = diagdata.upt.min, yend = diagdata.upt.max, col=grp), size = 1, lineend = "butt") + 
  #geom_point(data=diagdataM, aes(x=years+shift, y=(diagdata.upt.min+diagdata.upt.max)/2, col=grp), size=3, shape=19, alpha = 1) +
  scale_x_continuous(breaks = seq(min(diagdataM$years), max(diagdataM$years), by = 2)) +
  geom_line(data = diagdataM, aes(x = years+shift, y = diagmodel.upt.mean, col=grp), size = 1) + 
  geom_ribbon(data=diagdataM, aes(x=years+shift, ymin=diagmodel.upt.low, ymax=diagmodel.upt.upp, fill=grp), alpha = .3)  +
  scale_color_manual("", values = mycol[1:4]) +
  scale_fill_manual("", values = mycol[1:4]) +
  labs(x = "Years") +
  labs(y = "Yearly number of diagnoses pp") +
  ggtitle("Men") +
  theme_bw(base_family = "Times") +
  scale_y_continuous(breaks = seq(0,0.04, by=0.005)) +
  coord_cartesian(ylim=c(0,0.04)) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.justification=c(0.01,.99), legend.position=c(0.01,.99),text=element_text(size=8))

testplotF <- ggplot() +
  geom_segment(data = diagdataF, aes(x = years+shift, xend = years+shift, y = testdata.upt.min, yend = testdata.upt.max, col=grp), size = 1, lineend = "butt") +
  #geom_point(data=diagdataF, aes(x=years+shift, y=(testdata.upt.min+testdata.upt.max)/2, col=grp), size=3, shape=19, alpha = 1) +
  scale_x_continuous(breaks = seq(min(diagdataF$years), max(diagdataF$years), by = 2)) +
  scale_color_manual("", values = mycol[1:4]) +
  labs(x = "Years") +
  labs(y = "Yearly number of tests pp") +
  ggtitle("Women") +
  theme_bw(base_family = "Times") +
  scale_y_continuous(breaks = seq(0,0.4, by=0.05)) +
  coord_cartesian(ylim=c(0,0.4)) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.justification=c(0.01,.99), legend.position=c(0.01,.99),text=element_text(size=8))

diagplotF <- ggplot() +
  geom_segment(data = diagdataF, aes(x = years+shift, xend = years+shift, y = diagdata.upt.min, yend = diagdata.upt.max, col=grp), size = 1, lineend = "butt") + 
  #geom_point(data=diagdataF, aes(x=years+shift, y=(diagdata.upt.min+diagdata.upt.max)/2, col=grp), size=1, shape=15, alpha = 1) +
  scale_x_continuous(breaks = seq(min(diagdataF$years), max(diagdataF$years), by = 2)) +
  geom_line(data = diagdataF, aes(x = years+shift, y = diagmodel.upt.mean, col=grp), size = 1) + 
  geom_ribbon(data=diagdataF, aes(x=years+shift, ymin=diagmodel.upt.low, ymax=diagmodel.upt.upp, fill=grp), alpha = .3)  +
  scale_color_manual("", values = mycol[1:4]) +
  scale_fill_manual("", values = mycol[1:4]) +
  labs(x = "Years") +
  labs(y = "Yearly number of diagnoses pp") +
  ggtitle("Women") +
  theme_bw(base_family = "Times") +
  scale_y_continuous(breaks = seq(0,0.04, by=0.005)) +
  coord_cartesian(ylim=c(0,0.04)) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.justification=c(0.01,.99), legend.position=c(0.01,.99),text=element_text(size=8))

pdf(file.path(path_figures,"testsplot.pdf") ,width=7,height=4)
multiplot(plotlist=list(testplotM,testplotF), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

pdf(file.path(path_figures,"diagplot.pdf" ),width=7,height=4)
multiplot(plotlist=list(diagplotM,diagplotF), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

positivityrateplotM <- ggplot() +
  geom_segment(data = diagdataM, aes(x = years+shift, xend = years+shift, 
                                     y = diagdata.upt.min/((testdata.upt.min+testdata.upt.max)/2), 
                                     yend = diagdata.upt.max/((testdata.upt.min+testdata.upt.max)/2),
                                     col=grp), size = 1, lineend = "butt") + 
  geom_point(data=diagdataM, aes(x=years+shift, y=((diagdata.upt.min+diagdata.upt.max)/2)/((testdata.upt.min+testdata.upt.max)/2), col=grp), size=3, shape=19, alpha = 1) +
  geom_line(data = diagdataM, aes(x = years+shift, y = diagmodel.upt.mean/((testdata.upt.min+testdata.upt.max)/2), col=grp), size = 1) + 
  geom_ribbon(data=diagdataM, aes(x=years+shift,
                                  ymin=diagmodel.upt.low/((testdata.upt.min+testdata.upt.max)/2), 
                                  ymax=diagmodel.upt.upp/((testdata.upt.min+testdata.upt.max)/2), 
                                  fill=grp), alpha = .3)  +
  scale_x_continuous(breaks = seq(min(diagdataM$years), max(diagdataM$years), by = 2)) +
  scale_color_manual("", values = mycol[1:4]) +
  scale_fill_manual("", values = mycol[1:4]) +
  labs(x = "Years") +
  labs(y = "Positivity rate") +
  theme_bw(base_family = "Times") +
  coord_cartesian(ylim=c(0,0.4)) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.justification=c(0.99,.99), legend.position=c(0.99,.99),text=element_text(size=8))

positivityrateplotF <- ggplot() +
  geom_segment(data = diagdataF, aes(x = years+shift, xend = years+shift, 
                                     y = diagdata.upt.min/((testdata.upt.min+testdata.upt.max)/2), 
                                     yend = diagdata.upt.max/((testdata.upt.min+testdata.upt.max)/2),
                                     col=grp), size = 1, lineend = "butt") + 
  geom_point(data=diagdataF, aes(x=years+shift, y=((diagdata.upt.min+diagdata.upt.max)/2)/((testdata.upt.min+testdata.upt.max)/2), col=grp), size=3, shape=19, alpha = 1) +
  geom_line(data = diagdataF, aes(x = years+shift, y = diagmodel.upt.mean/((testdata.upt.min+testdata.upt.max)/2), col=grp), size = 1) + 
  geom_ribbon(data=diagdataF, aes(x=years+shift,
                                  ymin=diagmodel.upt.low/((testdata.upt.min+testdata.upt.max)/2), 
                                  ymax=diagmodel.upt.upp/((testdata.upt.min+testdata.upt.max)/2), 
                                  fill=grp), alpha = .3)  +
  scale_x_continuous(breaks = seq(min(diagdataF$years), max(diagdataF$years), by = 2)) +
  scale_color_manual("", values = mycol[1:4]) +
  scale_fill_manual("", values = mycol[1:4]) +
  labs(x = "Years") +
  labs(y = "Positivity rate") +
  theme_bw(base_family = "Times") +
  coord_cartesian(ylim=c(0,0.4)) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.justification=c(0.99,.99), legend.position=c(0.99,.99),text=element_text(size=8))


pdf(file.path(path_figures,"posrateplot.pdf" ),width=7,height=4)
multiplot(plotlist=list(positivityrateplotM,positivityrateplotF), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

####################################################################
# plot change of value eta in time (for 16-24 and 25-44 separately) and plot prevalences
####################################################################

compute_rates_time <- function(trace.mcmc.pars,fitmodel, years){
  
  source(file.path(path_in,"CT_screendata_simple.R"))
  chi_all <- get.CT.screendata(age.classes)
  
  # account for the different size of the age groups: number of tests done per modelled compartment (marginalized over activity classes)
  size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) #JS040416
  chi_all <- chi_all * size.age.classes
  
  age.class.groups <- age.classes[1:n.age.classes]
  plottedsims <- dim(trace.mcmc.pars)[1]
  
  eta_time <- to.tensor(NA, dims = list(sex=c("M","F"), age=paste0('age', age.class.groups), screening_y = years, nsim = 1:plottedsims))
  screenrate_inf <- to.tensor(NA, dims = list(sex=c("M","F"), age=paste0('age', age.class.groups), screening_y = years, nsim = 1:plottedsims)) # screening rate in infected people
  screenrate_tot <- to.tensor(NA, dims = list(sex=c("M","F"), age=paste0('age', age.class.groups), screening_y = years, nsim = 1:plottedsims)) # screening rate in all people
  
  times <- c(0, (burnintime - testperiod_in_burnin):max.simtime)
  
  
  for (s in 1:plottedsims){
    theta <- trace.mcmc.pars[s, 1:length(testmodel)]
    theta.allpars <- unlist(c(theta, tracetrace[1,which(!(fullmodel %in% testmodel))]))
    
    simulation <- fitmodel$simulate(theta.allpars, times)
    
    
    for (y in 1:length(years)){
      
      tensors <- initvec.to.tensor.cpp(simulation[which(times+2000-burnintime==years[y]),-1], par.cond = list(age_classes =age.classes[1:n.age.classes], nmb_act_classes = nJ))
      
      # compute screen_all
      sympttest_all <- theta.allpars["treat"] * (margin.tensor(tensors$I_S,i=2) + margin.tensor(tensors$I_S2,i=2)) # all tests for treatment in compartment
      alltest_all <- chi_all[,,which(dimnames(chi_all)$screening_y==as.character(years[y]))] # all tests in compartment
      screen_all <- alltest_all - sympttest_all # all screening tests in compartment
      
      eta_exp <- 1 + (theta.allpars["eta1"]-1) * exp(-10^theta.allpars["eta2"]*screen_all)
      
      eta_time[,,y,s] <- eta_exp
      
      
      screenrate_inf[,,y,s] <- eta_exp * screen_all /
        (margin.tensor(tensors$S,i=2) + margin.tensor(tensors$R,i=2) + eta_exp * (margin.tensor(tensors$I_A,i=2) + margin.tensor(tensors$I_A2,i=2)))
      
      screenrate_tot[,,y,s] <- screen_all /
        (margin.tensor(tensors$S,i=2) + margin.tensor(tensors$R,i=2) + margin.tensor(tensors$I_A,i=2) + margin.tensor(tensors$I_A2,i=2))
      
    }
    
  }
  
  sorted_eta_time <- eta_time
  for (g in 1:dim(sorted_eta_time)[1]){
    for (a in 1:dim(sorted_eta_time)[2]){
      for (y in 1:dim(sorted_eta_time)[3]){
        sorted_eta_time[g,a,y,] <- sort(eta_time[g,a,y,])
      }
    }
  }
  
  sorted_screenrate_inf <- screenrate_inf
  for (g in 1:dim(sorted_screenrate_inf)[1]){
    for (a in 1:dim(sorted_screenrate_inf)[2]){
      for (y in 1:dim(sorted_screenrate_inf)[3]){
        sorted_screenrate_inf[g,a,y,] <- sort(screenrate_inf[g,a,y,])
      }
    }
  }
  
  sorted_screenrate_tot <- screenrate_tot
  for (g in 1:dim(sorted_screenrate_tot)[1]){
    for (a in 1:dim(sorted_screenrate_tot)[2]){
      for (y in 1:dim(sorted_screenrate_tot)[3]){
        sorted_screenrate_tot[g,a,y,] <- sort(screenrate_tot[g,a,y,])
      }
    }
  }
  
  mean_eta_time <- apply(eta_time, MARGIN=c(1,2,3), FUN=mean)
  low_eta_time <- sorted_eta_time[,,,ceiling(0.025*dim(trace.mcmc.pars)[1])]
  upp_eta_time <- sorted_eta_time[,,,floor(0.975*dim(trace.mcmc.pars)[1])]
  
  mean_screenrate_inf <- apply(screenrate_inf, MARGIN=c(1,2,3), FUN=mean)
  low_screenrate_inf <- sorted_screenrate_inf[,,,ceiling(0.025*dim(trace.mcmc.pars)[1])]
  upp_screenrate_inf <- sorted_screenrate_inf[,,,floor(0.975*dim(trace.mcmc.pars)[1])]
  
  mean_screenrate_tot <- apply(screenrate_tot, MARGIN=c(1,2,3), FUN=mean)
  low_screenrate_tot <- sorted_screenrate_tot[,,,ceiling(0.025*dim(trace.mcmc.pars)[1])]
  upp_screenrate_tot <- sorted_screenrate_tot[,,,floor(0.975*dim(trace.mcmc.pars)[1])]
  
  return(list(eta_time=eta_time, mean_eta_time=mean_eta_time, low_eta_time=low_eta_time, upp_eta_time=upp_eta_time,
              screenrate_inf=screenrate_inf, mean_screenrate_inf=mean_screenrate_inf, low_screenrate_inf=low_screenrate_inf, upp_screenrate_inf=upp_screenrate_inf,
              screenrate_tot=screenrate_tot, mean_screenrate_tot=mean_screenrate_tot, low_screenrate_tot=low_screenrate_tot, upp_screenrate_tot=upp_screenrate_tot))
  
}

rates_time <- compute_rates_time(trace.combined, STI_MCMC_model, years=myyears)

# change of eta in time for different sex and age groups
gendernames <- c("M","F")
df_eta <- data.frame(years=myyears)

for (g in 1:2){
  
  df_eta[paste0("mean_",gendernames[g],"_1519")] <- as.vector(3/5*rates_time$mean_eta_time[g,1,] + 2/5*rates_time$mean_eta_time[g,2,])
  df_eta[paste0("mean_",gendernames[g],"_2024")] <- as.vector(rates_time$mean_eta_time[g,3,])
  df_eta[paste0("mean_",gendernames[g],"_2534")] <- as.vector(rates_time$mean_eta_time[g,4,])
  df_eta[paste0("mean_",gendernames[g],"_3544")] <- as.vector(rates_time$mean_eta_time[g,5,])
  
  df_eta[paste0("low_",gendernames[g],"_1519")] <- as.vector(3/5*rates_time$low_eta_time[g,1,] + 2/5*rates_time$low_eta_time[g,2,])
  df_eta[paste0("low_",gendernames[g],"_2024")] <- as.vector(rates_time$low_eta_time[g,3,])
  df_eta[paste0("low_",gendernames[g],"_2534")] <- as.vector(rates_time$low_eta_time[g,4,])
  df_eta[paste0("low_",gendernames[g],"_3544")] <- as.vector(rates_time$low_eta_time[g,5,])
  
  df_eta[paste0("upp_",gendernames[g],"_1519")] <- as.vector(3/5*rates_time$upp_eta_time[g,1,] + 2/5*rates_time$upp_eta_time[g,2,])
  df_eta[paste0("upp_",gendernames[g],"_2024")] <- as.vector(rates_time$upp_eta_time[g,3,])
  df_eta[paste0("upp_",gendernames[g],"_2534")] <- as.vector(rates_time$upp_eta_time[g,4,])
  df_eta[paste0("upp_",gendernames[g],"_3544")] <- as.vector(rates_time$upp_eta_time[g,5,])
  
}

etaM_timeplot_1519_2024 <- ggplot() +
  geom_line(data=df_eta, mapping=aes(x=years, y=mean_M_1519), colour=mycol[1]) +
  geom_ribbon(data=df_eta, aes(x=years, ymin=low_M_1519, ymax=upp_M_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_eta, mapping=aes(x=years, y=mean_M_2024), colour=mycol[2]) +
  geom_ribbon(data=df_eta, aes(x=years, ymin=low_M_2024, ymax=upp_M_2024), fill=mycol[2], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_eta$years), max(df_eta$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,8, by=1)) +
  coord_cartesian(ylim=c(0,8)) +
  labs(x = "Years") +
  labs(y = expression("Differential screening uptake ("*eta*" )")) +
  ggtitle("Men") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

etaF_timeplot_1519_2024 <- ggplot() +
  geom_line(data=df_eta, mapping=aes(x=years, y=mean_F_1519), colour=mycol[1]) +
  geom_ribbon(data=df_eta, aes(x=years, ymin=low_F_1519, ymax=upp_F_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_eta, mapping=aes(x=years, y=mean_F_2024), colour=mycol[2]) +
  geom_ribbon(data=df_eta, aes(x=years, ymin=low_F_2024, ymax=upp_F_2024), fill=mycol[2], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_eta$years), max(df_eta$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,8, by=1)) +
  coord_cartesian(ylim=c(0,8)) +
  labs(x = "Years") +
  labs(y = expression("Differential screening uptake ("*eta*" )")) +
  ggtitle("Women") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

pdf(file.path(path_figures,"eta_timeplot_1519_2024.pdf" ),width=7,height=4)
multiplot(plotlist=list(etaM_timeplot_1519_2024,etaF_timeplot_1519_2024), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

etaM_timeplot <- ggplot() +
  geom_line(data=df_eta, mapping=aes(x=years, y=4/9*mean_M_1519+5/9*mean_M_2024), colour=mycol[1]) +
  geom_ribbon(data=df_eta, aes(x=years, ymin=4/9*low_M_1519+5/9*low_M_2024, ymax=4/9*upp_M_1519+5/9*upp_M_2024), fill="black", alpha = .1) +
  scale_x_continuous(breaks = seq(min(df_eta$years), max(df_eta$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,8, by=1)) +
  coord_cartesian(ylim=c(0,8)) +
  labs(x = "Years") +
  labs(y = "Differential screening uptake") +
  ggtitle("Men (16-24 years old)") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

etaF_timeplot <- ggplot() +
  geom_line(data=df_eta, mapping=aes(x=years, y=4/9*mean_F_1519+5/9*mean_F_2024), colour=mycol[1]) +
  geom_ribbon(data=df_eta, aes(x=years, ymin=4/9*low_F_1519+5/9*low_F_2024, ymax=4/9*upp_F_1519+5/9*upp_F_2024), fill="black", alpha = .1) +
  scale_x_continuous(breaks = seq(min(df_eta$years), max(df_eta$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,8, by=1)) +
  coord_cartesian(ylim=c(0,8)) +
  labs(x = "Years") +
  labs(y = "Differential screening uptake") +
  ggtitle("Women (16-24 years old)") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

pdf(file.path(path_figures,"eta_timeplot.pdf" ),width=7,height=4)
multiplot(plotlist=list(etaM_timeplot,etaF_timeplot), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

# screen rate in all
gendernames <- c("M","F")
df_screenrate_tot <- data.frame(years=myyears)

for (g in 1:2){
  
  df_screenrate_tot[paste0("mean_",gendernames[g],"_1519")] <- as.vector(3/5*rates_time$mean_screenrate_tot[g,1,] + 2/5*rates_time$mean_screenrate_tot[g,2,])
  df_screenrate_tot[paste0("mean_",gendernames[g],"_2024")] <- as.vector(rates_time$mean_screenrate_tot[g,3,])
  df_screenrate_tot[paste0("mean_",gendernames[g],"_2534")] <- as.vector(rates_time$mean_screenrate_tot[g,4,])
  df_screenrate_tot[paste0("mean_",gendernames[g],"_3544")] <- as.vector(rates_time$mean_screenrate_tot[g,5,])
  
  df_screenrate_tot[paste0("low_",gendernames[g],"_1519")] <- as.vector(3/5*rates_time$low_screenrate_tot[g,1,] + 2/5*rates_time$low_screenrate_tot[g,2,])
  df_screenrate_tot[paste0("low_",gendernames[g],"_2024")] <- as.vector(rates_time$low_screenrate_tot[g,3,])
  df_screenrate_tot[paste0("low_",gendernames[g],"_2534")] <- as.vector(rates_time$low_screenrate_tot[g,4,])
  df_screenrate_tot[paste0("low_",gendernames[g],"_3544")] <- as.vector(rates_time$low_screenrate_tot[g,5,])
  
  df_screenrate_tot[paste0("upp_",gendernames[g],"_1519")] <- as.vector(3/5*rates_time$upp_screenrate_tot[g,1,] + 2/5*rates_time$upp_screenrate_tot[g,2,])
  df_screenrate_tot[paste0("upp_",gendernames[g],"_2024")] <- as.vector(rates_time$upp_screenrate_tot[g,3,])
  df_screenrate_tot[paste0("upp_",gendernames[g],"_2534")] <- as.vector(rates_time$upp_screenrate_tot[g,4,])
  df_screenrate_tot[paste0("upp_",gendernames[g],"_3544")] <- as.vector(rates_time$upp_screenrate_tot[g,5,])
  
}

screenrateM_tot_timeplot <- ggplot() +
  geom_line(data=df_screenrate_tot, mapping=aes(x=years, y=mean_M_1519), colour=mycol[1]) +
  geom_ribbon(data=df_screenrate_tot, aes(x=years, ymin=low_M_1519, ymax=upp_M_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_screenrate_tot, mapping=aes(x=years, y=mean_M_2024), colour=mycol[2]) +
  geom_ribbon(data=df_screenrate_tot, aes(x=years, ymin=low_M_2024, ymax=upp_M_2024), fill=mycol[2], alpha = .3) +
  geom_line(data=df_screenrate_tot, mapping=aes(x=years, y=mean_M_2534), colour=mycol[3]) +
  geom_ribbon(data=df_screenrate_tot, aes(x=years, ymin=low_M_2534, ymax=upp_M_2534), fill=mycol[3], alpha = .3) +
  geom_line(data=df_screenrate_tot, mapping=aes(x=years, y=mean_M_3544), colour=mycol[4]) +
  geom_ribbon(data=df_screenrate_tot, aes(x=years, ymin=low_M_3544, ymax=upp_M_3544), fill=mycol[4], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_screenrate_tot$years), max(df_screenrate_tot$years), by = 2)) +
  labs(x = "Years") +
  labs(y = "Screening rates per year") +
  ggtitle("Men") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

screenrateF_tot_timeplot <- ggplot() +
  geom_line(data=df_screenrate_tot, mapping=aes(x=years, y=mean_F_1519), colour=mycol[1]) +
  geom_ribbon(data=df_screenrate_tot, aes(x=years, ymin=low_F_1519, ymax=upp_F_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_screenrate_tot, mapping=aes(x=years, y=mean_F_2024), colour=mycol[2]) +
  geom_ribbon(data=df_screenrate_tot, aes(x=years, ymin=low_F_2024, ymax=upp_F_2024), fill=mycol[2], alpha = .3) +
  geom_line(data=df_screenrate_tot, mapping=aes(x=years, y=mean_F_2534), colour=mycol[3]) +
  geom_ribbon(data=df_screenrate_tot, aes(x=years, ymin=low_F_2534, ymax=upp_F_2534), fill=mycol[3], alpha = .3) +
  geom_line(data=df_screenrate_tot, mapping=aes(x=years, y=mean_F_3544), colour=mycol[4]) +
  geom_ribbon(data=df_screenrate_tot, aes(x=years, ymin=low_F_3544, ymax=upp_F_3544), fill=mycol[4], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_screenrate_tot$years), max(df_screenrate_tot$years), by = 2)) +
  labs(x = "Years") +
  labs(y = "Screening rates per year") +
  ggtitle("Women") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

pdf(file.path(path_figures,"screenrate_tot_timeplot.pdf") ,width=7,height=4)
multiplot(plotlist=list(screenrateM_tot_timeplot,screenrateF_tot_timeplot), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

# screen rate in infecteds (rate xi in model)
gendernames <- c("M","F")
df_screenrate_inf <- data.frame(years=myyears)

for (g in 1:2){
  
  df_screenrate_inf[paste0("mean_",gendernames[g],"_1519")] <- as.vector(3/5*rates_time$mean_screenrate_inf[g,1,] + 2/5*rates_time$mean_screenrate_inf[g,2,])
  df_screenrate_inf[paste0("mean_",gendernames[g],"_2024")] <- as.vector(rates_time$mean_screenrate_inf[g,3,])
  df_screenrate_inf[paste0("mean_",gendernames[g],"_2534")] <- as.vector(rates_time$mean_screenrate_inf[g,4,])
  df_screenrate_inf[paste0("mean_",gendernames[g],"_3544")] <- as.vector(rates_time$mean_screenrate_inf[g,5,])
  
  df_screenrate_inf[paste0("low_",gendernames[g],"_1519")] <- as.vector(3/5*rates_time$low_screenrate_inf[g,1,] + 2/5*rates_time$low_screenrate_inf[g,2,])
  df_screenrate_inf[paste0("low_",gendernames[g],"_2024")] <- as.vector(rates_time$low_screenrate_inf[g,3,])
  df_screenrate_inf[paste0("low_",gendernames[g],"_2534")] <- as.vector(rates_time$low_screenrate_inf[g,4,])
  df_screenrate_inf[paste0("low_",gendernames[g],"_3544")] <- as.vector(rates_time$low_screenrate_inf[g,5,])
  
  df_screenrate_inf[paste0("upp_",gendernames[g],"_1519")] <- as.vector(3/5*rates_time$upp_screenrate_inf[g,1,] + 2/5*rates_time$upp_screenrate_inf[g,2,])
  df_screenrate_inf[paste0("upp_",gendernames[g],"_2024")] <- as.vector(rates_time$upp_screenrate_inf[g,3,])
  df_screenrate_inf[paste0("upp_",gendernames[g],"_2534")] <- as.vector(rates_time$upp_screenrate_inf[g,4,])
  df_screenrate_inf[paste0("upp_",gendernames[g],"_3544")] <- as.vector(rates_time$upp_screenrate_inf[g,5,])
  
}

screenrateM_inf_timeplot <- ggplot() +
  geom_line(data=df_screenrate_inf, mapping=aes(x=years, y=mean_M_1519), colour=mycol[1]) +
  geom_ribbon(data=df_screenrate_inf, aes(x=years, ymin=low_M_1519, ymax=upp_M_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_screenrate_inf, mapping=aes(x=years, y=mean_M_2024), colour=mycol[2]) +
  geom_ribbon(data=df_screenrate_inf, aes(x=years, ymin=low_M_2024, ymax=upp_M_2024), fill=mycol[2], alpha = .3) +
  geom_line(data=df_screenrate_inf, mapping=aes(x=years, y=mean_M_2534), colour=mycol[3]) +
  geom_ribbon(data=df_screenrate_inf, aes(x=years, ymin=low_M_2534, ymax=upp_M_2534), fill=mycol[3], alpha = .3) +
  geom_line(data=df_screenrate_inf, mapping=aes(x=years, y=mean_M_3544), colour=mycol[4]) +
  geom_ribbon(data=df_screenrate_inf, aes(x=years, ymin=low_M_3544, ymax=upp_M_3544), fill=mycol[4], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_screenrate_inf$years), max(df_screenrate_inf$years), by = 2)) +
  labs(x = "Years") +
  labs(y = expression("Screening rate asymptomatically infected men ("*chi^A*" )")) +
  ggtitle("Men") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

screenrateF_inf_timeplot <- ggplot() +
  geom_line(data=df_screenrate_inf, mapping=aes(x=years, y=mean_F_1519), colour=mycol[1]) +
  geom_ribbon(data=df_screenrate_inf, aes(x=years, ymin=low_F_1519, ymax=upp_F_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_screenrate_inf, mapping=aes(x=years, y=mean_F_2024), colour=mycol[2]) +
  geom_ribbon(data=df_screenrate_inf, aes(x=years, ymin=low_F_2024, ymax=upp_F_2024), fill=mycol[2], alpha = .3) +
  geom_line(data=df_screenrate_inf, mapping=aes(x=years, y=mean_F_2534), colour=mycol[3]) +
  geom_ribbon(data=df_screenrate_inf, aes(x=years, ymin=low_F_2534, ymax=upp_F_2534), fill=mycol[3], alpha = .3) +
  geom_line(data=df_screenrate_inf, mapping=aes(x=years, y=mean_F_3544), colour=mycol[4]) +
  geom_ribbon(data=df_screenrate_inf, aes(x=years, ymin=low_F_3544, ymax=upp_F_3544), fill=mycol[4], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_screenrate_inf$years), max(df_screenrate_inf$years), by = 2)) +
  labs(x = "Years") +
  labs(y = expression("Screening rate asymptomatically infected women ("*chi^A*" )")) +
  ggtitle("Women") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

pdf(file.path(path_figures,"screenrate_inf_timeplot.pdf") ,width=7,height=4)
multiplot(plotlist=list(screenrateM_inf_timeplot,screenrateF_inf_timeplot), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

###########################################################################################################
### ANALYSIS OF WHAT WOULD HAVE HAPPENED WITH PREVALENCE IF ETA2 had been constant
###########################################################################################################

trace.combined.etaconstant <- trace.combined
trace.combined.etaconstant["eta2"] <- -6
postsims_eta2constant <- simulate_prevs_diag(trace.combined.etaconstant, STI_MCMC_model, years=myyears)

# change of prevalence in time for different sex and age groups
gendernames <- c("M","F")
df_prev_eta2constant <- data.frame(years=myyears)

for (g in 1:2){
  
  df_prev_eta2constant[paste0("mean_",gendernames[g],"_1519")] <- as.vector(3/5*postsims_eta2constant$mean_prevalenceT[g,1,] + 2/5*postsims_eta2constant$mean_prevalenceT[g,2,])
  df_prev_eta2constant[paste0("mean_",gendernames[g],"_2024")] <- as.vector(postsims_eta2constant$mean_prevalenceT[g,3,])
  df_prev_eta2constant[paste0("mean_",gendernames[g],"_2534")] <- as.vector(postsims_eta2constant$mean_prevalenceT[g,4,])
  df_prev_eta2constant[paste0("mean_",gendernames[g],"_3544")] <- as.vector(postsims_eta2constant$mean_prevalenceT[g,5,])
  
  df_prev_eta2constant[paste0("low_",gendernames[g],"_1519")] <- as.vector(3/5*postsims_eta2constant$low_prevalenceT[g,1,] + 2/5*postsims_eta2constant$low_prevalenceT[g,2,])
  df_prev_eta2constant[paste0("low_",gendernames[g],"_2024")] <- as.vector(postsims_eta2constant$low_prevalenceT[g,3,])
  df_prev_eta2constant[paste0("low_",gendernames[g],"_2534")] <- as.vector(postsims_eta2constant$low_prevalenceT[g,4,])
  df_prev_eta2constant[paste0("low_",gendernames[g],"_3544")] <- as.vector(postsims_eta2constant$low_prevalenceT[g,5,])
  
  df_prev_eta2constant[paste0("upp_",gendernames[g],"_1519")] <- as.vector(3/5*postsims_eta2constant$upp_prevalenceT[g,1,] + 2/5*postsims_eta2constant$upp_prevalenceT[g,2,])
  df_prev_eta2constant[paste0("upp_",gendernames[g],"_2024")] <- as.vector(postsims_eta2constant$upp_prevalenceT[g,3,])
  df_prev_eta2constant[paste0("upp_",gendernames[g],"_2534")] <- as.vector(postsims_eta2constant$upp_prevalenceT[g,4,])
  df_prev_eta2constant[paste0("upp_",gendernames[g],"_3544")] <- as.vector(postsims_eta2constant$upp_prevalenceT[g,5,])
  
}

prevtime_eta2constant_plot_M <- ggplot() +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=mean_M_1519), colour=mycol[1]) +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=low_M_1519, ymax=upp_M_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=mean_M_2024), colour=mycol[2]) +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=low_M_2024, ymax=upp_M_2024), fill=mycol[2], alpha = .3) +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=mean_M_2534), colour=mycol[3]) +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=low_M_2534, ymax=upp_M_2534), fill=mycol[3], alpha = .3) +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=mean_M_3544), colour=mycol[4]) +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=low_M_3544, ymax=upp_M_3544), fill=mycol[4], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_prev_eta2constant$years), max(df_prev_eta2constant$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,0.06, by=0.01)) +
  coord_cartesian(ylim=c(0,0.06)) +
  labs(x = "Years") +
  labs(y = "Prevalence") +
  ggtitle("Men") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))


prevtime_eta2constant_plot_F <- ggplot() +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=mean_F_1519), colour=mycol[1]) +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=low_F_1519, ymax=upp_F_1519), fill=mycol[1], alpha = .3) +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=mean_F_2024), colour=mycol[2]) +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=low_F_2024, ymax=upp_F_2024), fill=mycol[2], alpha = .3) +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=mean_F_2534), colour=mycol[3]) +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=low_F_2534, ymax=upp_F_2534), fill=mycol[3], alpha = .3) +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=mean_F_3544), colour=mycol[4]) +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=low_F_3544, ymax=upp_F_3544), fill=mycol[4], alpha = .3) +
  scale_x_continuous(breaks = seq(min(df_prev_eta2constant$years), max(df_prev_eta2constant$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,0.06, by=0.01)) +
  coord_cartesian(ylim=c(0,0.06)) +
  labs(x = "Years") +
  labs(y = "Prevalence") +
  ggtitle("Women") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

pdf(file.path(path_figures,"prevtime_eta2constant_plot.pdf") ,width=7,height=4)
multiplot(plotlist=list(prevtime_eta2constant_plot_M,prevtime_eta2constant_plot_F), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

# plot fitted prevalence and prevalence in hypothetical scenario in same plot, 16-24 years old
prevtime_eta2var_eta2constant_plot_M <- ggplot() +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=4/9*mean_M_1519+5/9*mean_M_2024), colour=mycol[1],linetype = "dashed") +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=4/9*low_M_1519+5/9*low_M_2024), colour=mycol[1],linetype = "dashed") +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=4/9*upp_M_1519+5/9*upp_M_2024), colour=mycol[1],linetype = "dashed") +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=4/9*low_M_1519+5/9*low_M_2024, ymax=4/9*upp_M_1519+5/9*upp_M_2024), fill="black", alpha = .1) +
  geom_line(data=df_prev, mapping=aes(x=years, y=4/9*mean_M_1519+5/9*mean_M_2024), colour=mycol[1],linetype = "solid") +
  geom_line(data=df_prev, mapping=aes(x=years, y=4/9*low_M_1519+5/9*low_M_2024), colour=mycol[1],linetype = "solid") +
  geom_line(data=df_prev, mapping=aes(x=years, y=4/9*upp_M_1519+5/9*upp_M_2024), colour=mycol[1],linetype = "solid") +
  geom_ribbon(data=df_prev, aes(x=years, ymin=4/9*low_M_1519+5/9*low_M_2024, ymax=4/9*upp_M_1519+5/9*upp_M_2024), fill="black", alpha = .1) +
  scale_x_continuous(breaks = seq(min(df_prev_eta2constant$years), max(df_prev_eta2constant$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,0.06, by=0.01)) +
  coord_cartesian(ylim=c(0,0.06)) +
  annotate("text", x=2010, y=0.025, label= "best-fit", size=3) +
  annotate("text", x=2009.5, y=0.015, label= "hypothetical", size=3) +
  labs(x = "Years") +
  labs(y = "Prevalence") +
  ggtitle("Men (16-24 years old)") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

prevtime_eta2var_eta2constant_plot_F <- ggplot() +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=4/9*mean_F_1519+5/9*mean_F_2024), colour=mycol[1],linetype = "dashed") +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=4/9*low_F_1519+5/9*low_F_2024), colour=mycol[1],linetype = "dashed") +
  geom_line(data=df_prev_eta2constant, mapping=aes(x=years, y=4/9*upp_F_1519+5/9*upp_F_2024), colour=mycol[1],linetype = "dashed") +
  geom_ribbon(data=df_prev_eta2constant, aes(x=years, ymin=4/9*low_F_1519+5/9*low_F_2024, ymax=4/9*upp_F_1519+5/9*upp_F_2024), fill="black", alpha = .1) +
  geom_line(data=df_prev, mapping=aes(x=years, y=4/9*mean_F_1519+5/9*mean_F_2024), colour=mycol[1],linetype = "solid") +
  geom_line(data=df_prev, mapping=aes(x=years, y=4/9*low_F_1519+5/9*low_F_2024), colour=mycol[1],linetype = "solid") +
  geom_line(data=df_prev, mapping=aes(x=years, y=4/9*upp_F_1519+5/9*upp_F_2024), colour=mycol[1],linetype = "solid") +
  geom_ribbon(data=df_prev, aes(x=years, ymin=4/9*low_F_1519+5/9*low_F_2024, ymax=4/9*upp_F_1519+5/9*upp_F_2024), fill="black", alpha = .1) +
  scale_x_continuous(breaks = seq(min(df_prev_eta2constant$years), max(df_prev_eta2constant$years), by = 2)) +
  scale_y_continuous(breaks = seq(0,0.06, by=0.01)) +
  coord_cartesian(ylim=c(0,0.06)) +
  annotate("text", x=2010, y=0.025, label= "best-fit", size=3) +
  annotate("text", x=2009.5, y=0.015, label= "hypothetical", size=3) +
  labs(x = "Years") +
  labs(y = "Prevalence") +
  ggtitle("Women (16-24 years old)") +
  theme_bw(base_family = "Times") +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=8))

pdf(file.path(path_figures,"prevtime_eta2var_eta2constant_plot.pdf" ),width=7,height=4)
multiplot(plotlist=list(prevtime_eta2var_eta2constant_plot_M,prevtime_eta2var_eta2constant_plot_F), layout=matrix(c(1,2), nrow=1, byrow=TRUE))
dev.off()

##############################################################################
# return non-pdf output
##############################################################################

mcmcdata <- list(accrate=accrate,
                   effectiveSize=effectiveSize,
                   summary.trace.mcmc.params.burned=summary.trace.mcmc.params.burned,
                   Rhat=Rhat,
                   D.theta.bar=D.theta.bar,
                   p.D=p.D,
                   DIC=DIC,
                   pW=pW,
                   pW.alt=pW.alt,
                   WAIC=WAIC)
save(mcmcdata, file=file.path(path_figures,"mcmcdata.RData" ))

plotsummarydata <- list(plotdataM.2000=plotdataM.2000,plotdataF.2000=plotdataF.2000,plotdataM.2011=plotdataM.2011,plotdataF.2011=plotdataF.2011,
                        df_prev=df_prev,diagdataM=diagdataM,diagdataF=diagdataF,df_eta=df_eta,
                        df_screenrate_tot=df_screenrate_tot,df_screenrate_inf=df_screenrate_inf,
                        df_prev_eta2constant=df_prev_eta2constant)
save(plotsummarydata, file=file.path(path_figures,"plotsummarydata.RData" ))

