# This function is a modification of the existing mcmcMH from fitR package
# It can be used to take an existing mcmc chain as input and extend it with another n.iterations iterations
my_mcmcMH <- function (target, init.theta = NULL, proposal.sd = NULL, n.iterations, 
          covmat = NULL, limits = list(lower = NULL, upper = NULL), 
          adapt.size.start = NULL, adapt.size.cooling = 0.99, adapt.shape.start = NULL, 
          print.info.every = n.iterations/100, verbose = FALSE, max.scaling.sd = 50,
          acceptance.rate.weight = NULL, acceptance.window = NULL, existingtrace=NULL) 
{
  
  if (is.null(init.theta) & is.null(existingtrace)){
    stop("need to input either initial theta or existing trace", call. = FALSE)
  } else if (!is.null(init.theta) & is.null(existingtrace)){
    theta.current <- init.theta
    theta.mean <- theta.current
  } else if (is.null(init.theta) & !is.null(existingtrace)){
    theta.current <- unlist(existingtrace$trace[dim(existingtrace$trace)[1],]) # includes log.density
    theta.current <- theta.current[-length(theta.current)]# remove log.density
    
    tracemat <- matrix(as.vector(unlist(existingtrace$trace)), 
                       nrow=dim(existingtrace$trace)[1], ncol=dim(existingtrace$trace)[2])
    colnames(tracemat) <- colnames(existingtrace$trace)
    theta.mean <- colMeans(tracemat)
    theta.mean <- theta.mean[-length(theta.mean)]# remove log.density
  }
  
  theta.propose <- theta.current
  lower.proposal <- limits$lower
  upper.proposal <- limits$upper
  theta.names <- names(theta.current)
  
  if (is.null(lower.proposal)) {
    lower.proposal <- init.theta
    lower.proposal[] <- -Inf
  } else {
    lower.proposal <- lower.proposal[theta.names]
  }
  if (is.null(upper.proposal)) {
    upper.proposal <- init.theta
    upper.proposal[] <- Inf
  } else {
    upper.proposal <- upper.proposal[theta.names]
  }
  
  covmat.proposal <- covmat
  
  if (!is.null(proposal.sd) && is.null(names(proposal.sd))) {
    names(proposal.sd) <- theta.names
  }
  
  if (is.null(existingtrace)){#JS280717
    
    if (is.null(covmat.proposal)) {
      if (is.null(proposal.sd)) {
        proposal.sd <- init.theta/10
      }
      covmat.proposal <- matrix(diag(proposal.sd[theta.names]^2, 
                                     nrow = length(theta.names)), nrow = length(theta.names), 
                                dimnames = list(theta.names, theta.names))
    } else {
      covmat.proposal <- covmat.proposal[theta.names, theta.names]
    }
    
  } else { #JS280717
    if (is.null(covmat.proposal)) { #if chain has to be continued using all information of existing chain #JS280717
      covmat.proposal <- existingtrace$covmat.proposal #JS280717
    } else { #if chain has to be continued using another proposal covariance matrix (e.g. of another chain) #JS280717
      covmat.proposal <- covmat.proposal[theta.names, theta.names] #JS280717
    } #JS280717
    
  } #JS280717
  
  adapting.size <- FALSE
  adapting.shape <- FALSE
  
  if (is.null(existingtrace)){ #JS15052017

    covmat.empirical <- covmat.proposal
    covmat.empirical[, ] <- 0
    
    acceptance.rate <- 0
    first.iteration <- 1
  
  } else {
    
    covmat.empirical <- existingtrace$covmat.empirical
      
    acceptance.rate <- existingtrace$acceptance.rate
    first.iteration <- dim(tracemat)[1]+1
  
  }
  
  covmat.proposal.init <- covmat.proposal
  theta.estimated.names <- names(which(diag(covmat.proposal) >  0))
  
  target.theta.current <- target(theta.current)
  if (class(target.theta.current) == "numeric") {
    target.theta.current <- list(log.density = target.theta.current, 
                                 trace = theta.current)
  }
  
  if (!is.null(print.info.every)) {
    message("Init: ", printNamedVector(theta.current[theta.estimated.names]), 
            ", target: ", target.theta.current[["log.density"]])
  }
  
  if (is.null(existingtrace)){
    trace <- data.frame(t(c(target.theta.current[["trace"]], target.theta.current["log.density"])))
  } else {
    trace <- rbind(existingtrace$trace, c(target.theta.current[["trace"]], 
                            target.theta.current[["log.density"]])) 
  }
  
  if (!is.null(acceptance.window)) {
    acceptance.series <- c()
  }
  scaling.sd <- 1
  scaling.multiplier <- 1
  
  if (is.null(print.info.every)) {
    print.info.every <- n.iterations + 1
  }

  start_iteration_time <- Sys.time()
  
  for (i.iteration in first.iteration:(first.iteration+n.iterations-1)) { #JS15052017
    if (!is.null(adapt.size.start) && i.iteration >= adapt.size.start && 
        (is.null(adapt.shape.start) || acceptance.rate * 
         i.iteration < adapt.shape.start)) {
      if (!adapting.size) {
        message("\n---> Start adapting size of covariance matrix")
        adapting.size <- TRUE
      }
      scaling.multiplier <- exp(adapt.size.cooling^(i.iteration - 
                                                      adapt.size.start) * (acceptance.rate - 0.234))
      scaling.sd <- scaling.sd * scaling.multiplier
      scaling.sd <- min(c(scaling.sd, max.scaling.sd))
      covmat.proposal.new <- scaling.sd^2 * covmat.proposal.init
      if (!(any(diag(covmat.proposal.new)[theta.estimated.names] < .Machine$double.eps))) {
        covmat.proposal <- covmat.proposal.new
      }
    }
    else if (!is.null(adapt.shape.start) && acceptance.rate * i.iteration >= adapt.shape.start) {
      if (!adapting.shape) {
        message("\n---> Start adapting shape of covariance matrix")
        adapting.shape <- TRUE
      }
      scaling.sd <- 2.38/sqrt(length(theta.estimated.names))
      covmat.proposal <- scaling.sd^2 * covmat.empirical
    }
    if (i.iteration%%ceiling(print.info.every) == 0) {
      state.mcmc <- target.theta.current$trace
      message("Iteration: ", i.iteration, "/", n.iterations+first.iteration-1, 
              ", acceptance rate: ", sprintf("%.3f", acceptance.rate), 
              appendLF = FALSE)
      if (!is.null(adapt.size.start) || !is.null(adapt.shape.start)) {
        message(", scaling.sd: ", sprintf("%.3f", scaling.sd), 
                ", scaling.multiplier: ", sprintf("%.3f", scaling.multiplier), 
                appendLF = FALSE)
      }
      message(", state: ", printNamedVector(state.mcmc))
      message(", logdensity: ", target.theta.current$log.density)
    }
    if (any(diag(covmat.proposal)[theta.estimated.names] < 
            .Machine$double.eps)) {
      print(covmat.proposal[theta.estimated.names, theta.estimated.names])
      stop("non-positive definite covmat", call. = FALSE)
    }
    if (length(theta.estimated.names) > 0) {
      theta.propose[theta.estimated.names] <- as.vector(rtmvnorm(1, 
                                                                 mean = theta.current[theta.estimated.names], 
                                                                 sigma = covmat.proposal[theta.estimated.names, 
                                                                                         theta.estimated.names], lower = lower.proposal[theta.estimated.names], 
                                                                 upper = upper.proposal[theta.estimated.names]))
    }
    
    #modified 130417
    try(target.theta.propose <- target(theta.propose), silent = TRUE)
    if(class(target.theta.propose) == "try-error") browser()
    
    
    #target.theta.propose <- target(theta.propose)
    if (class(target.theta.propose) == "numeric") {
      target.theta.propose <- list(log.density = target.theta.propose, 
                                   trace = theta.propose)
    }
    if (!is.finite(target.theta.propose$log.density)) {
      log.acceptance <- -Inf
    }
    else {
      log.acceptance <- target.theta.propose$log.density - 
        target.theta.current$log.density
      log.acceptance <- log.acceptance + dtmvnorm(x = theta.current[theta.estimated.names], 
                                                  mean = theta.propose[theta.estimated.names], 
                                                  sigma = covmat.proposal[theta.estimated.names, 
                                                                          theta.estimated.names], lower = lower.proposal[theta.estimated.names], 
                                                  upper = upper.proposal[theta.estimated.names], 
                                                  log = TRUE)
      log.acceptance <- log.acceptance - dtmvnorm(x = theta.propose[theta.estimated.names], 
                                                  mean = theta.current[theta.estimated.names], 
                                                  sigma = covmat.proposal[theta.estimated.names, 
                                                                          theta.estimated.names], lower = lower.proposal[theta.estimated.names], 
                                                  upper = upper.proposal[theta.estimated.names], 
                                                  log = TRUE)
    }
    if (verbose) {
      message("Propose: ", theta.propose[theta.estimated.names], 
              ", target: ", target.theta.propose[["log.density"]], 
              ", acc prob: ", exp(log.acceptance), ", ", appendLF = FALSE)
    }
    if (is.accepted <- (log(runif(1)) < log.acceptance)) {
      theta.current <- theta.propose
      target.theta.current <- target.theta.propose
      if (verbose) {
        message("accepted")
      }
    }
    else if (verbose) {
      message("rejected")
    }

    trace <- rbind(trace, c(target.theta.current[["trace"]], 
                            target.theta.current[["log.density"]])) #WAS THIS A MISTAKE?? CHANGED 29/11/2016 JS
    if (i.iteration == 1) {
      acceptance.rate <- is.accepted
    }
    else {
      if (is.null(acceptance.rate.weight)) {
        if (is.null(acceptance.window)) {
          acceptance.rate <- acceptance.rate + (is.accepted - 
                                                  acceptance.rate)/i.iteration
        }
        else {
          acceptance.series <- c(is.accepted, acceptance.series)
          if (length(acceptance.series) > acceptance.window) {
            acceptance.series <- acceptance.series[-length(acceptance.series)]
          }
          acceptance.rate <- mean(acceptance.series)
        }
      }
      else {
        acceptance.rate <- acceptance.rate * (1 - acceptance.rate.weight) + 
          is.accepted * acceptance.rate.weight
      }
    }
    tmp <- updateCovmat(covmat.empirical, theta.mean, theta.current, 
                        i.iteration)
    covmat.empirical <- tmp$covmat
    theta.mean <- tmp$theta.mean
  }
  return(list(trace = trace, acceptance.rate = acceptance.rate, 
              covmat.empirical = covmat.empirical,
              covmat.proposal = covmat.proposal))
}
