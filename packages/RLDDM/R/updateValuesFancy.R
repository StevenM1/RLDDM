updateValuesFancy <- function(outcomes, values, eta1, eta2) {
  
  nTrials <- nrow(outcomes)
  nChoices <- ncol(outcomes)
  if(length(eta1) == 1) {
    eta1 <- rep(eta1, nTrials)
  }
  if(length(eta2) == 1) {
    eta2 <- rep(eta2, nTrials)
  }
  
  # initialize here in R, because we need these values
  VV <- PE <- matrix(nrow=nTrials, ncol=nChoices)
  
  out = .C('doubleUpdateValuesCMC',
           nTrials=nTrials,
           nChoices=nChoices,
           values=as.double(values),
           VV=as.double(VV),
           PE=as.double(PE),
           outcomes=as.double(outcomes),
           eta1=as.double(eta1),
           eta2=as.double(eta2), NAOK=TRUE)
  VV <- matrix(out$VV, nrow=nTrials, ncol=nChoices)
  PE <- matrix(out$PE, nrow=nTrials, ncol=nChoices)
  list(VV=VV, PE=PE)
}



# test
compareCvsR <- function(choice, outcomes, values, eta1, eta2) {
  library(zoo)

  nt <- length(choice)  # number of trials
  nc <- 2  # number of choice options (left/right)
  R.PE <- numeric(nt) # PE per trial
  values <- c(.5, .5)

  # C
  updated <- updateValuesFancy(outcomes=outcomes, values=values, eta1, eta2)
  C.VV <- updated$VV
  C.PE <- updated$PE

  # R
  # make sure values are reset
  outcome = apply(outcome, 1, function(x) x[!is.na(x)])
  
  values <- c(.5, .5)
  R.VV <- matrix(NA, nrow=nt, ncol=2)
  
  # Update all values (can we vectorize this? Hopefully)
  for(trial in 1:nt) {
    c_ = choice[trial] # c_ = choice
    o_ = outcome[trial] # o_ = outcome
    
    # Keep track of current value, used later for drift rates
    R.VV[trial, c_] = values[c_]
    R.VV[trial, -c_] = values[-c_]
    
    # calculate PE
    dv = o_-values[c_]  # prediction error = outcome (reward) - predicted value
    R.PE[trial] <- dv  # keep track of this
    
    # update values
    values[c_] = values[c_] + ifelse(dv>0, eta1, eta2)*dv
  }
  

  # forward fill NAs
  all(R.VV==C.VV)
}

# 
# test
# nTrials <- 20
# notChosen <- sample(c(1, 2), size=nTrials, replace=TRUE)
# eta1 <- .3
# eta2 <- .4
# values <- c(.5, .5)
# outcome <- matrix(rnorm(nTrials*2, 0, 3), nrow=nTrials)
# outcome[cbind(1:nTrials, as.numeric(notChosen))] <- NA

# tmp = updateValuesFancy(outcomes=outcome, values=values, eta1, eta2)
# 

#
# #debug(compareCvsR)
# compareCvsR(ifelse(notChosen==1, 2, 1), outcome, values, eta1, eta2)
# 
