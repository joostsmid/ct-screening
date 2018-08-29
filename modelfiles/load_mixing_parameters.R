# SMOOTHED AGE MIXING DISTRIBUTIONS, USING SKEWNORMAL FIT
# Parameter values are from article Smid, J. H., Garcia, V., Low, N., Mercer, C. H. & Althaus, C. 
# Age difference between heterosexual partners: implications for the spread of Chlamydia trachomatis. Epidemiology, doi:10.1016/j.epidem.2018.03.004 (2018).
# https://doi.org/10.1016/j.epidem.2018.03.004 (Supplementary materials, Figure S1)

##############################################################################################################################
### age related sexual contact preference matrix 
##############################################################################################################################

rho_age <- to.tensor(rep(0, n.age.classes^2*2), dims = list(sex=c("M","F"), age=paste0('age', age.class.groups), ageprime = paste0('age', age.class.groups) ))

### construct mixing matrix

#converts central moments into parameters sn distribution: xi, omega, alpha
moment.to.parsSN <- function(mu,sd,sk){
  absdelta <- sqrt( (pi/2) * abs(sk)^(2/3) / (((4-pi)/2)^(2/3) + abs(sk)^(2/3)) )
  if (sk<0){absdelta <- -absdelta}
  alpha <- absdelta / sqrt(1-absdelta^2)
  omega <- sqrt(sd^2 / (1 - 2 * absdelta^2 / pi) )
  xi <- mu - omega * absdelta * sqrt(2/pi)
  return(c(xi,omega,alpha))
}

# Read in fitted parameters (central moments as function of participant age)
intercept.mu.m <- 6.146126
slope.mu.m <- 0.6812418
intercept.sd.m <- -3.709598
slope.sd.m <- 0.4032263
slope2.sd.m <- -0.003054981
steepness.sk.m <- -0.1226717
midpoint.sk.m <- 38.93827
 
intercept.mu.f <- 6.907002
slope.mu.f <- 0.7968689
intercept.sd.f <- -4.344449
slope.sd.f <- 0.5368711
slope2.sd.f <- -0.005689269
steepness.sk.f <- -0.1308938
midpoint.sk.f <- 37.01026
  
mm <- matrix(NA,n.age.classes,n.age.classes)
mf <- matrix(NA,n.age.classes,n.age.classes)

#extended partner ages: only used for plotting
n.age.classes.partner <- length(age.classes.partner)-1
mm.extended <- matrix(NA,n.age.classes,n.age.classes.partner)
mf.extended <- matrix(NA,n.age.classes,n.age.classes.partner)

for (k in 1:n.age.classes){ #loop over age classes of participants
  sn.backtrans.par.m <- moment.to.parsSN(mu=intercept.mu.m + slope.mu.m * age.class.groups[k],
                                       sd=intercept.sd.m + slope.sd.m * age.class.groups[k] + slope2.sd.m * age.class.groups[k]^2,
                                       sk=-2/(1+exp(steepness.sk.m * (age.class.groups[k] - midpoint.sk.m )))+1)
  
  if (any(is.na(sn.backtrans.par.m))) {cat("error in (back)-transformation skew-negative distribution")}
  else {
    xm <- dsn(age.class.groups, xi= sn.backtrans.par.m[1], omega= sn.backtrans.par.m[2],alpha=sn.backtrans.par.m[3])
    mm[k,] <- xm / sum(xm) #normalize
    
    xm.extended <- dsn(age.classes.partner[1:n.age.classes.partner], xi= sn.backtrans.par.m[1], omega= sn.backtrans.par.m[2],alpha=sn.backtrans.par.m[3])  #JS210716
    mm.extended[k,] <- xm.extended / sum(xm.extended) #normalize
  }
  
  sn.backtrans.par.f <- moment.to.parsSN(mu=intercept.mu.f + slope.mu.f * age.class.groups[k],
                                       sd=intercept.sd.f + slope.sd.f * age.class.groups[k] + slope2.sd.f * age.class.groups[k]^2,#sk=min(.995, intercept.sk.f + slope.sk.f * age.class.groups[k]))
                                       sk=-2/(1+exp(steepness.sk.f * (age.class.groups[k] - midpoint.sk.f )))+1)

  if (any(is.na(sn.backtrans.par.f))) {cat("error in (back)-transformation skew-negative distribution")}
  
  else {
    xf <- dsn(age.class.groups, xi= sn.backtrans.par.f[1], omega= sn.backtrans.par.f[2], alpha=sn.backtrans.par.f[3])
    mf[k,] <- xf / sum(xf) #normalize
    
    
    xf.extended <- dsn(age.classes.partner[1:n.age.classes.partner], xi= sn.backtrans.par.f[1], omega= sn.backtrans.par.f[2],alpha=sn.backtrans.par.f[3])  #JS210716
    mf.extended[k,] <- xf.extended / sum(xf.extended) #normalize
  }
}

mm[mm<1E-7]<-0 #otherwise problems may occur in solving odes 
mf[mf<1E-7]<-0

mm <- sweep(mm, 1, rowSums(mm), FUN="/") #normalize across rows
mf <- sweep(mf, 1, rowSums(mf), FUN="/") #normalize across rows

### fill tensor up
for(a in 1:n.age.classes){
  rho_age["F", a, ] <- mf[a,]
  rho_age["M", a, ] <- mm[a,]
}

##############################################################################################################################
### adjustment for contact rates among genders
##############################################################################################################################

if (exists("parameters")){
  sex <- c("M","F") ### gender
  dimnames1 <- list(sex=sex)
  nElements <- length(sex)
  theta <- to.tensor(c(0.5,0.5), dims = dimnames1)
  
  parameters$rho_age <- rho_age
  parameters$theta <- theta
}