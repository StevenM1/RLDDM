updateValuesDouble <- function(choice, outcome1, outcome2, values, PE1, PE2, eta1, eta2) {
  
  # for output
  VV1 <- VV2 <- numeric(length(choice))
  
  out = .C('doubleUpdateValuesC',
           nTrials=as.integer(length(choice)),
           value1=as.double(values[1]),
           value2=as.double(values[2]),
           VV1=as.double(VV1),
           VV2=as.double(VV2),
           choice=as.integer(choice),
           outcome1=as.integer(outcome1),
           outcome2=as.integer(outcome2),
           PE1=as.double(PE1),
           PE2=as.double(PE2),
           eta1=as.double(eta1),
           eta2=as.double(eta2))
  VV <- cbind(out$VV1, out$VV2)
  list(VV=VV, PE1=out$PE1, PE2=out$PE2)
}



# test
compareCvsR <- function(choice, outcome1, outcome2, values, eta1, eta2) {
  library(zoo)
  
  nt <- length(choice)  # number of trials
  nc <- 2  # number of choice options (left/right)
  PE1.C <- PE2.C <- PE1.R <- PE2.R <- numeric(nt) # PE per trial
  values <- c(.5, .5)
  
  # C
  updated <- updateValuesDouble(choice, outcome1, outcome2, values, PE1.C, PE2.C, eta1, eta2)
  C.VV <- updated$VV
  C.PE <- updated$PE
  
  # R
  # make sure values are reset
  values <- c(.5, .5)
  R.VV <- matrix(NA, nrow=nt, ncol=2)
  
  # Update all values (can we vectorize this? Hopefully)
  for(trial in 1:nt) {
    o1_ = outcome1[trial] # o_ = outcome
    o2_ = outcome2[trial]
    
    # Keep track of current value, used later for drift rates
    R.VV[trial, 1] = values[1]
    R.VV[trial, 2] = values[2]
    
    # calculate PE
    dv1 = o1_ - values[1]  # prediction error = outcome (reward) - predicted value
    dv2 = o2_ - values[2]
    PE1.R[trial] <- dv1  # keep track of this
    PE2.R[trial] <- dv2  # keep track of this
    
    # update values
    values[1] = values[1] + ifelse(dv1>0, eta1, eta2)*dv1
    values[2] = values[2] + ifelse(dv2>0, eta1, eta2)*dv2
  }
  
  # forward fill NAs
  all(R.VV==C.VV)
}

# test
# choice <- sample(c(1, 2), size=20, replace=TRUE)
# eta1 <- .3
# eta2 <- .4
# values <- c(.5, .5)
# outcome1 <- sample(c(1, 0), size=20, replace=TRUE)
# outcome2 <- sample(c(1, 0), size=20, replace=TRUE)
# PE1 <- PE2 <- numeric(20)
# updateValuesDouble(choice=choice, outcome1=outcome1, outcome2=outcome2, values=values, PE1=PE1, PE2=PE2, eta1=eta1, eta2=eta2)
# 
# #debug(compareCvsR)
# compareCvsR(choice, outcome1=outcome1, outcome2=outcome2, values, eta1, eta2)

