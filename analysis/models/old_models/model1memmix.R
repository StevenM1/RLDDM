### Model 1
# a ~ Cue
# # mixture distribution of 2 DDMs, with a "memory mix-up" mixture proportion
# we assume that on a subset of the trials, the value differnece is accumulated properly
# but the stimuli are mis-identified

### modelSpec is a list containing:
# 1. The parameters to fit, and the factors they depend on
# 2. constants in the model
# 3. The factors from (1), and their levels
modelSpec = list('variablePars'=list('a' = 'condition',
                                     'mix' = 1,
                                     't0' = 1,
                                     'eta1' = 1, 
                                     'm'=1),
                 'constants'=c('sv'=0, 'sz'=0, 'st0'=0, 'z'=0.5, 's'=1, 'eta2'=-Inf),
                 'condition'=c('SPD', 'ACC'),
                 'learningRule'='Qlearning')

obj <- objRLDDMMultiCond

### transformLearningRate is a function transforming
### "global" parameters to trial-by-trial values, dependent 
### on the condition
transformLearningRate <- function(pars, condition) {
  # "Declare"
  eta1 <- eta2 <- rep(pars[['eta1']], length(condition))
  return(list(eta1=eta1, eta2=eta2))
}

### the following function gets trial-by-trial DDM pars
transformDDMPars <- function(pars, condition, delta_ev) {
  ### Gets trial-by-trial DDM parameters ###
  nTrials = length(condition)
  a <- v <- t0 <- z <- sv <- sz <- s <- rep(NA, nTrials)
  
  # all current models have no variability in sz, sv, s, t0
  t0 = rep(pars[['t0']], nTrials)
  z = rep(pars[['z']], nTrials)
  sv <- rep(pars[['sv']], nTrials)
  sz <- rep(pars[['sz']], nTrials)
  s <- rep(pars[['s']], nTrials)
  st0 <- rep(pars[['st0']], nTrials)
  
  # all models assume a linear relation between delta_ev and v
  v = delta_ev*pars[['m']]
  
  # a differs by condition
  a[condition=='SPD'] <- pars[['a.SPD']]
  a[condition=='ACC'] <- pars[['a.ACC']]
  
  # rescale z from [0, 1] to [0, a]
  z = z*a
  
  # repeat all pars by 2
  t0 <- rep(t0, 2)
  a <- rep(a, 2)
  z <- rep(z, 2)
  sz <- rep(sz, 2)
  sv <- rep(sv, 2)
  s <- rep(s, 2)
  st0 <- rep(s, 2)
  v <- c(v, -v)

  return(list(t0=t0, a=a, v=v, z=z, sz=sz, sv=sv, s=s, st0=st0))
}

obj <- function(pars,
                rt, choice, condition, outcomes, parNames, learningRule="SARSA",
                values=c(0.5, 0.5), VVchoiceIdx, constants=NULL, 
                returnType='full', min.like=-1e20,
                backend = 'c') {
  ### Transcoded from Hanneke den Ouden's Matlab code (2015, see below)
  
  ## Learning model is Rescorla-Wagner:
  # vA <- vA + eta*(r-vA)   -> vA = value of option A, and is updated with predcition error (reward-vA) modulated by learning rate eta
  
  ## Choice function is DDM:
  # a = Threshold
  # v = drift rate = difference in EV, multiplied by a linear scaling factor
  # t0 = non-decision time
  
  # name parameters, add constants to vector
  names(pars) <- parNames
  if(!is.null(constants)) for(i in 1:length(constants)) pars[[names(constants[i])]] <- constants[[i]]
  etas <- transformLearningRate(condition=condition, 
                                pars=pars)
  
  # Update values trial-by-trial
  updated <- updateValuesFancy(outcomes, values, etas$eta1, etas$eta2, learningRule=learningRule)
  VV <- updated$VV
  PE <- updated$PE
  
  # On the basis of which values is chosen? Here, we check in 'outcomes' for which alternatives got feedback, and use this to index in VV
  VV_choice = matrix(VV[VVchoiceIdx], ncol=2, byrow=TRUE)
  EV_diff = VV_choice[,1] - VV_choice[,2]  # difference in expected value over trials
  
  # DDM choice function
  ddmPars <- transformDDMPars(pars, condition, delta_ev=EV_diff)
  like <- ddiffusion(rt=rep(rt, 2), 
                     response=rep(choice, rt), 
                     a=ddmPars[['a']], 
                     v=ddmPars[['v']], 
                     t0=ddmPars[['t0']],
                     sz=ddmPars[['sz']],
                     z=ddmPars[['z']],
                     sv=ddmPars[['sz']],
                     s=ddmPars[['s']],
                     st0=ddmPars[['st0']])*c(pars[['mix']], 1-pars[['mix']])
  nt <- length(rt)
  like <- like[1:nt] + like[(nt:1):(nt*2)]
  LL <- sum(log(like))
  if(any(is.na(LL)) | any(is.infinite(LL)) | any(ddmPars[['s']] < 0)) {
    LL <- min.like
  }
  
  if(returnType=='full') {
    # pdiffusion requires the RTs to be ordered, so lets order
    PE_choice = matrix(PE[VVchoiceIdx], ncol=2, byrow=TRUE)
    return(list(VV=VV_choice, VV_full=VV, LL=LL, PE_full=PE, PE=PE_choice))
  } else {
    return(-LL)
  }
}


defaultBounds <- function() {
  list('a'=c(1e-3, 5),
       'v'=c(1e-3, 5),
       'm'=c(1e-3, 250),
       'z'=c(.2, .8),
       't0'=c(0, .5),
       'eta1'=c(0, 1),
       'eta2'=c(0, 1),
       'sz'=c(0, 5),
       'sv'=c(0, 5),
       's'=c(0, 5),
       'vmax'=c(0, 100),
       'k'=c(0, 50),
       'st0'=c(0, .5),
       'beta'=c(0, 50),
       'mix'=c(0.5, 1))
}
