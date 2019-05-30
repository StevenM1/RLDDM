updateValues <- function(choice, outcome, values, PE, eta1, eta2) {

  # for output
  VV1 <- VV2 <- numeric(length(choice))

  out = .C('updateValuesC',
           nTrials=as.integer(length(choice)),
           value1=as.double(values[1]),
           value2=as.double(values[2]),
           VV1=as.double(VV1),
           VV2=as.double(VV2),
           choice=as.integer(choice),
           outcome=as.integer(outcome),
           PE=as.double(PE),
           eta1=as.double(eta1),
           eta2=as.double(eta2))
  VV <- cbind(out$VV1, out$VV2)
  list(VV=VV, PE=out$PE)
}



# test

compareCvsR <- function(choice, outcome, values, eta1, eta2) {
  library(zoo)
  
  nt <- length(choice)  # number of trials
  nc <- 2  # number of choice options (left/right)
  PE.C <- PE.R <- numeric(nt) # PE per trial
  values <- c(.5, .5)
  
  # C
  updated <- updateValues(choice, outcome, values, PE.C, eta1, eta2)
  C.VV <- updated$VV
  C.PE <- updated$PE
  
  # R
  # make sure values are reset
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
    PE.R[trial] <- dv  # keep track of this
    
    # update values
    values[c_] = values[c_] + ifelse(dv>0, eta1, eta2)*dv
  }
  
  # forward fill NAs
  all(R.VV==C.VV)
}

# test
# choice <- sample(c(1, 2), size=20, replace=TRUE)
# eta1 <- .3
# eta2 <- .4
# values <- c(.5, .5)
# outcome <- sample(c(1, 0), size=20, replace=TRUE)
# PE <- numeric(20)
# updateValues(choice=choice, outcome=outcome, values=values, PE=PE, eta1=eta1, eta2=eta2)
# 
# #debug(compareCvsR)
# compareCvsR(choice, outcome, values, eta1, eta2)
