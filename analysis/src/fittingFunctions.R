# Functions required for fitting RL-DDM
prepareForFitting <- function(dat) {
  df <- dat[!is.na(dat$rt),] # remove nonresponse trials
  df <- df[df$rt>.1,] # remove extremely fast responses
  
  # "Condition" here is the stimulus set under consideration. This should be a numerical (1, 2, 3, 4 ...)
  df$condition <- df$stimulus_set-min(df$stimulus_set)
  df$condition <- match(df$condition, unique(df$condition))-1
  
  # Get response "correctness" - i.e., whether a choice corresponds to the optimal choice
  df$response.corr <- df$choice <- ifelse(df$choiceIsHighP, 2, 1) # code choice in terms of upper bound (2) or lower bound (1)
  df$reward <- df$outcome / 100  # re-code rewards to 0/1, assuming an identity value function
  
  # define bins per stimulus; useful for later checking model fit
  df$trialN_this_stim <- NA
  for(lvl in unique(df$condition)) {
    df$trialN_this_stim[df$condition==lvl] <- seq(1, sum(df$condition==lvl))
  }
  df$bin <- as.numeric(cut(df$trialN_this_stim, 5))
  
  # Set-up outcome matrix
  outcomes <- matrix(NA, nrow=nrow(df), ncol=length(unique(df$condition))*2)
  for(row in 1:nrow(df)) {
    cond = df$condition[row]
    outcomes[row,(cond)*2+ifelse(df$choice[row]==1, 2, 1)] <- df$reward[row]
  }
  
  # Set-up values vector
  values <- rep(0.5, ncol(outcomes))
  
  # On the basis of which alternatives is chosen? 
  # Make a matrix of nTrials x 2; first column = trialN, second column = choice number
  VVchoiceIdx <- matrix(FALSE, nrow=nrow(df), ncol=ncol(outcomes))
  for(tr in 1:nrow(df)) {
    condition <- df[tr, 'condition']
    VVchoiceIdx[tr, ((condition)*2+1):((condition)*2+2)] <- TRUE
  }
  VVchoiceIdx <- which(t(VVchoiceIdx), arr.ind = TRUE)[,2:1]  # Gives for every trial each column that is chosen
  choice <- df$choice
  return(list(choice=choice, VVchoiceIdx=VVchoiceIdx, outcomes=outcomes, df=df, values=values))
}

obj <- function(pars, rt, choice, condition, outcomes, parNames, learningRule='SARSA',
                choiceFunction='SM', 
                values=c(0.5, 0.5), VVchoiceIdx, constants=NULL, 
                returnType='full', min.like=-1e20,
                backend = 'c') {

  # name parameters, add constants to vector
  names(pars) <- parNames
  if(!is.null(constants)) for(i in 1:length(constants)) pars[[names(constants[i])]] <- constants[[i]]
  etas <- transformLearningRate(condition=condition, pars=pars)
  
  # Update values trial-by-trial
  updated <- updateValuesFancy(outcomes, values, etas$eta1, etas$eta2, learningRule=learningRule)
  VV <- updated$VV
  PE <- updated$PE
  
  # On the basis of which values is chosen? Here, we check in 'outcomes' for which alternatives got feedback, 
  # and use this to index in VV
  VV_choice = matrix(VV[VVchoiceIdx], ncol=2, byrow=TRUE)
  
  if(choiceFunction == 'SM') {
    # Softmax
    choicePars <- transformChoicePars(pars, condition)
    PP <- exp(VV_choice*choicePars[['beta']])
    PP <- PP/apply(PP, 1, sum)
    ## note that choice == 2 is "correct", which corresponds to the *first* column in the PP matrix
    ## Sorry for that confusion
    LL <- sum(log(PP[choice==1,2])) + sum(log(PP[choice==2,1]))	
  } else if(choiceFunction == 'DDM') {
    # DDM choice function
    choicePars <- transformChoicePars(pars, condition, VV_choice)
    like <- ddiffusion(rt=rt, response=choice, 
                       a=choicePars[['a']], 
                       v=choicePars[['v']], 
                       t0=choicePars[['t0']],
                       sz=choicePars[['sz']],
                       z=choicePars[['z']],
                       sv=choicePars[['sv']],
                       s=choicePars[['s']],
                       st0=choicePars[['st0']])
    LL <- sum(log(like))
  } else if(choiceFunction == 'RW') {
    # Racing Walds
    choicePars <- transformChoicePars(pars, condition, VV_choice)
    like <- dWaldRace(rt=rt, 
              # In the data, choice == 2 in data corresponds to correct answer (nb: this is what rdiffusion/ddiffusion assume)
              response=choice,  
              v=list(choicePars[['v2']], choicePars[['v1']]),
              B=choicePars[['B']], 
              A=choicePars[['A']], 
              t0=choicePars[['t0']], 
              st0=0,  # NOT (YET) IMPLEMENTED!
              silent=TRUE)
    LL <- sum(log(like))
  } else if(choiceFunction == 'LBA') {
    ### NB not tested!
    choicePars <- transformChoicePars(pars, condition, VV_choice)
    like <- dLBA(rt=rt, response=choice, 
                 A=choicePars[['A']], 
                 b=choicePars[['b']], 
                 t0=choicePars[['t0']],
                 mean_v=choicePars[['mean_v']],
                 sd_v=choicePars[['sd_v']], 
                 silent=TRUE)
    LL <- sum(log(like))
  }

  if(any(is.na(LL)) | any(is.infinite(LL))) {
    LL <- min.like
  }
  
  if(returnType=='full') {
    PE_choice = matrix(PE[VVchoiceIdx], ncol=2, byrow=TRUE)
    return(list(VV=VV_choice, VV_full=VV, LL=LL, PE_full=PE, PE=PE_choice))
  } else {
    return(-LL)
  }
}

# Model set-up ------------------------------------------------------------
defaultBounds <- function() {
  list('a'=c(1e-3, 5),
       'v'=c(1e-3, 5),
       'm'=c(1e-3, 50),
       'z'=c(.2, .8),
       't0'=c(0, .5),
       'eta1'=c(0, 1),
       'eta2'=c(0, 1),
       'sz'=c(0, 5),
       'sv'=c(0, 2),
       's'=c(0, 5),
       'vmax'=c(0, 100),
       'k'=c(0, 50),
       'st0'=c(0, .5),
       'beta'=c(0, 50))
}

model <- function(modelSetup, bounds=defaultBounds()) {
  # make p.vector
  p.vector <- c()
  for(parName in names(modelSetup$variablePars)) {
    valueInList = modelSetup$variablePars[[parName]]
    if(valueInList == 1) {
      # intercept only
      p.vector <- c(p.vector, parName)
    } else {
      # dependent on condition
      p.vector <- c(p.vector, paste(parName, modelSetup[[valueInList]], sep='.'))
    }
  }
  
  # get bounds
  lower <- upper <- c()
  for(par in p.vector) {
    parameterCharacter = strsplit(par, '.', fixed=TRUE)[[1]][1]
    lower <- c(lower, bounds[[parameterCharacter]][1])
    upper <- c(upper, bounds[[parameterCharacter]][2])
  }
  
  # print('p.vector:')
  # print(allPars)
  # print('Constants:')
  # print(model$constants)
  return(list(p.vector=p.vector, 
              lowerBounds=lower,
              upperBounds=upper,
              constants=modelSetup$constants,
              learningRule=modelSetup$learningRule,
              choiceFunction=modelSetup$choiceFunction))
}

# 
# # Objective function ------------------------------------------------------
# objRLDDMMultiCond <- function(pars,
#                               rt, choice, condition, outcomes, parNames, learningRule="SARSA",
#                               values=c(0.5, 0.5), VVchoiceIdx, constants=NULL, 
#                               returnType='full', min.like=-1e20,
#                               backend = 'c') {
#   ### Transcoded from Hanneke den Ouden's Matlab code (2015, see below)
#   
#   ## Learning model is Rescorla-Wagner:
#   # vA <- vA + eta*(r-vA)   -> vA = value of option A, and is updated with predcition error (reward-vA) modulated by learning rate eta
#   
#   ## Choice function is DDM:
#   # a = Threshold
#   # v = drift rate = difference in EV, multiplied by a linear scaling factor
#   # t0 = non-decision time
#   
#   # name parameters, add constants to vector
#   names(pars) <- parNames
#   if(!is.null(constants)) for(i in 1:length(constants)) pars[[names(constants[i])]] <- constants[[i]]
#   etas <- transformLearningRate(condition=condition, 
#                                 pars=pars)
#   
#   # Update values trial-by-trial
#   updated <- updateValuesFancy(outcomes, values, etas$eta1, etas$eta2, learningRule=learningRule)
#   VV <- updated$VV
#   PE <- updated$PE
#   
#   # On the basis of which values is chosen? Here, we check in 'outcomes' for which alternatives got feedback, and use this to index in VV
#   VV_choice = matrix(VV[VVchoiceIdx], ncol=2, byrow=TRUE)
#   EV_diff = VV_choice[,1] - VV_choice[,2]  # difference in expected value over trials
#   
#   # DDM choice function
#   ddmPars <- transformDDMPars(pars, condition, delta_ev=EV_diff)
#   like <- ddiffusion(rt=rt, response=choice, 
#                      a=ddmPars[['a']], 
#                      v=ddmPars[['v']], 
#                      t0=ddmPars[['t0']],
#                      sz=ddmPars[['sz']],
#                      z=ddmPars[['z']],
#                      sv=ddmPars[['sv']],
#                      s=ddmPars[['s']],
#                      st0=ddmPars[['st0']])
#   LL <- sum(log(like))
#   if(any(is.na(LL)) | any(is.infinite(LL)) | any(ddmPars[['s']] < 0)) {
#     LL <- min.like
#   }
#   
#   if(returnType=='full') {
#     # pdiffusion requires the RTs to be ordered, so lets order
#     idx <- order(rt)
#     PP <- pdiffusion(rt[order(rt)], response=rep(2, length(choice)), 
#                      a=ddmPars[['a']][order(rt)], 
#                      v=ddmPars[['v']][order(rt)], 
#                      t0=ddmPars[['t0']][order(rt)],
#                      sz=ddmPars[['sz']][order(rt)],
#                      z=ddmPars[['z']][order(rt)],
#                      sv=ddmPars[['sz']][order(rt)],
#                      s=ddmPars[['s']][order(rt)])
#     # and reverse order
#     PP <- PP[order(order(rt))]
#     PE_choice = matrix(PE[VVchoiceIdx], ncol=2, byrow=TRUE)
#     return(list(VV=VV_choice, VV_full=VV, LL=LL, PE_full=PE, PE=PE_choice))
#   } else {
#     return(-LL)
#   }
# }
# 
# objRLRWMultiCond <- function(pars,
#                              rt, choice, condition, outcomes, parNames, learningRule="SARSA",
#                              values=c(0.5, 0.5), VVchoiceIdx, constants=NULL, 
#                              returnType='full', min.like=-1e20,
#                              backend = 'c') {
#   ### Transcoded from Hanneke den Ouden's Matlab code (2015, see below)
#   
#   ## Learning model is Rescorla-Wagner:
#   # vA <- vA + eta*(r-vA)   -> vA = value of option A, and is updated with predcition error (reward-vA) modulated by learning rate eta
#   
#   ## Choice function is DDM:
#   # a = Threshold
#   # v = drift rate = difference in EV, multiplied by a linear scaling factor
#   # t0 = non-decision time
#   
#   # name parameters, add constants to vector
#   names(pars) <- parNames
#   if(!is.null(constants)) for(i in 1:length(constants)) pars[[names(constants[i])]] <- constants[[i]]
#   etas <- transformLearningRate(condition=condition, 
#                                 pars=pars)
#   if('startValue' %in% names(pars)) {
#     startValue = pars[['startValue']]
#     idx1 <- seq(1, ncol(outcomes), 2)
#     idx2 <- idx1 + 1
#     values[idx1] <- values[idx1] + startValue/2
#     values[idx2] <- values[idx2] - startValue/2
#   }
#   # Update values trial-by-trial
#   updated <- updateValuesFancy(outcomes, values, etas$eta1, etas$eta2, learningRule=learningRule)
#   VV <- updated$VV
#   PE <- updated$PE
#   
#   # On the basis of which values is chosen? Here, we check in 'outcomes' for which alternatives got feedback, and use this to index in VV
#   VV_choice = matrix(VV[VVchoiceIdx], ncol=2, byrow=TRUE)
# #  EV_diff = VV_choice[,1] - VV_choice[,2]  # difference in expected value over trials
#   # VV_choice[,1] corresponds to the "correct" answer
#   # VV_choice[,2] to the "incorrect" answer
#   
#   # DDM choice function
#   choicePars <- transformChoicePars(pars, condition, VV_choice)
#   like <- dWaldRace(rt=rt, 
#                     # In the data, choice == 2 in data corresponds to correct answer (nb: this is what rdiffusion/ddiffusion assume)
#                     # so flip 1 and 2 to make sure accumulator 1 corresponds to choice 1 = correct
#                     response=choice,#ifelse(choice==2, 1, 2),  
#                     v=list(choicePars[['v1']], choicePars[['v2']]),
#                     B=choicePars[['B']], 
#                     A=choicePars[['A']], 
#                     t0=choicePars[['t0']], 
#                     st0=0,  # NOT YET IMPLEMENTED! 
#                     silent=TRUE)
#   
#   LL <- sum(log(like))
#   if(any(is.na(LL)) | any(is.infinite(LL))) {
#     LL <- min.like
#   }
#   
#   if(returnType=='full') {
#     # pdiffusion requires the RTs to be ordered, so lets order
#     PE_choice = matrix(PE[VVchoiceIdx], ncol=2, byrow=TRUE)
#     return(list(VV=VV_choice, VV_full=VV, LL=LL, PE_full=PE, PE=PE_choice))
#   } else {
#     return(-LL)
#   }
# }
# 
# objSoftmaxMultiCond <- function(pars,
#                               rt, choice, condition, outcomes, parNames, learningRule="Qlearning",
#                               values=c(0.5, 0.5), VVchoiceIdx, constants=NULL, 
#                               returnType='full', min.like=-1e20,
#                               backend = 'c') {
#   ### Transcoded from Hanneke den Ouden's Matlab code (2015, see below)
#   
#   ## Learning model is Rescorla-Wagner:
#   # vA <- vA + eta*(r-vA)   -> vA = value of option A, and is updated with predcition error (reward-vA) modulated by learning rate eta
#   
#   ## Choice function is DDM:
#   # a = Threshold
#   # v = drift rate = difference in EV, multiplied by a linear scaling factor
#   # t0 = non-decision time
#   
#   # name parameters, add constants to vector
#   names(pars) <- parNames
#   if(!is.null(constants)) for(i in 1:length(constants)) pars[[names(constants[i])]] <- constants[[i]]
#   etas <- transformLearningRate(condition=condition, 
#                                 pars=pars)
#   
#   # Update values trial-by-trial
#   updated <- updateValuesFancy(outcomes, values, etas$eta1, etas$eta2, learningRule=learningRule)
#   VV <- updated$VV
#   PE <- updated$PE
#   
#   # On the basis of which values is chosen? Here, we check in 'outcomes' for which alternatives got feedback, and use this to index in VV
#   VV_choice = matrix(VV[VVchoiceIdx], ncol=2, byrow=TRUE)
#   EV_diff = VV_choice[,1] - VV_choice[,2]  # difference in expected value over trials
#   
#   # DDM choice function
# #  ddmPars <- transformChoicePars(pars, condition, delta_ev=EV_diff)
#   choicePars <- transformChoicePars(pars, condition)
#   PP <- exp(VV_choice*choicePars[['beta']])
#   PP <- PP/apply(PP, 1, sum)
#   LL <- sum(log(PP[choice==1,1])) + sum(log(PP[choice==2,2]))	 ## note that choice == 2 is "correct", which corresponds to the *first* column in the PP matrix
#   
#   if(any(is.na(LL)) | any(is.infinite(LL))) {
#     LL <- min.like
#   }
#   
#   if(returnType=='full') {
#     PE_choice = matrix(PE[VVchoiceIdx], ncol=2, byrow=TRUE)
#     return(list(VV=VV_choice, VV_full=VV, LL=LL, PE_full=PE, PE=PE_choice))
#   } else {
#     return(-LL)
#   }
# }
# objDDMMultiCond <- function(pars,
#                             rt, choice, condition, outcomes=NULL, parNames,
#                             values=NULL, VVchoiceIdx=NULL, constants=NULL, 
#                             returnType='full', min.like=-1e20,
#                             backend = 'c') {
#   ### Transcoded from Hanneke den Ouden's Matlab code (2015, see below)
#   
#   ## Learning model is Rescorla-Wagner:
#   # vA <- vA + eta*(r-vA)   -> vA = value of option A, and is updated with predcition error (reward-vA) modulated by learning rate eta
#   
#   ## Choice function is DDM:
#   # a = Threshold
#   # v = drift rate = difference in EV, multiplied by a linear scaling factor
#   # t0 = non-decision time
#   
#   # name parameters, add constants to vector
#   names(pars) <- parNames
#   if(!is.null(constants)) for(i in 1:length(constants)) pars[[names(constants[i])]] <- constants[[i]]
#   # etas <- transformLearningRate(condition=condition, 
#   #                               pars=pars)
#   
#   # Update values trial-by-trial
#   # updated <- updateValuesFancy(outcomes, values, etas$eta1, etas$eta2)
#   # VV <- updated$VV
#   # PE <- updated$PE
#   
#   # On the basis of which values is chosen? Here, we check in 'outcomes' for which alternatives got feedback, and use this to index in VV
#   # VV_choice = matrix(VV[VVchoiceIdx], ncol=2, byrow=TRUE)
#   # EV_diff = VV_choice[,1] - VV_choice[,2]  # difference in expected value over trials
#   
#   # DDM choice function
#   ddmPars <- transformDDMPars(pars, condition)
#   like <- ddiffusion(rt=rt, response=choice, 
#                      a=ddmPars[['a']], 
#                      v=ddmPars[['v']], 
#                      t0=ddmPars[['t0']],
#                      sz=ddmPars[['sz']],
#                      z=ddmPars[['z']],
#                      sv=ddmPars[['sz']],
#                      s=ddmPars[['s']],
#                      st0=ddmPars[['st0']])
#   LL <- sum(log(like))
#   if(any(is.na(LL)) | any(is.infinite(LL)) | any(ddmPars[['s']] < 0)) {
#     LL <- min.like
#   }
#   
#   if(returnType=='full') {
#     # pdiffusion requires the RTs to be ordered, so lets order
#     idx <- order(rt)
#     PP <- pdiffusion(rt[order(rt)], response=rep(2, length(choice)), 
#                      a=ddmPars[['a']][order(rt)], 
#                      v=ddmPars[['v']][order(rt)], 
#                      t0=ddmPars[['t0']][order(rt)],
#                      sz=ddmPars[['sz']][order(rt)],
#                      z=ddmPars[['z']][order(rt)],
#                      sv=ddmPars[['sz']][order(rt)],
#                      s=ddmPars[['s']][order(rt)])
#     # and reverse order
#     PP <- PP[order(order(rt))]
#     # PE_choice = matrix(PE[VVchoiceIdx], ncol=2, byrow=TRUE)
#     return(list(VV=NULL, VV_full=NULL, LL=LL, PE_full=NULL, PE=NULL))
#   } else {
#     return(-LL)
#   }
# }
# 
# 
# objLBAMultiCond <- function(pars,
#                             rt, choice, condition, outcomes=NULL, parNames,
#                             values=NULL, VVchoiceIdx=NULL, constants=NULL, 
#                             returnType='full', min.like=-1e20,
#                             backend = 'c') {
#   ### Transcoded from Hanneke den Ouden's Matlab code (2015, see below)
#   
#   ## Learning model is Rescorla-Wagner:
#   # vA <- vA + eta*(r-vA)   -> vA = value of option A, and is updated with predcition error (reward-vA) modulated by learning rate eta
#   
#   ## Choice function is DDM:
#   # a = Threshold
#   # v = drift rate = difference in EV, multiplied by a linear scaling factor
#   # t0 = non-decision time
#   
#   # name parameters, add constants to vector
#   names(pars) <- parNames
#   if(!is.null(constants)) for(i in 1:length(constants)) pars[[names(constants[i])]] <- constants[[i]]
#   # etas <- transformLearningRate(condition=condition, 
#   #                               pars=pars)
#   
#   # Update values trial-by-trial
#   # updated <- updateValuesFancy(outcomes, values, etas$eta1, etas$eta2)
#   # VV <- updated$VV
#   # PE <- updated$PE
#   
#   # On the basis of which values is chosen? Here, we check in 'outcomes' for which alternatives got feedback, and use this to index in VV
#   # VV_choice = matrix(VV[VVchoiceIdx], ncol=2, byrow=TRUE)
#   # EV_diff = VV_choice[,1] - VV_choice[,2]  # difference in expected value over trials
#   
#   # DDM choice function
#   lbaPars <- transformLBAPars(pars, condition)
#   like <- dLBA(rt=rt, response=choice, 
#                A=lbaPars[['A']], 
#                b=lbaPars[['b']], 
#                t0=lbaPars[['t0']],
#                mean_v=lbaPars[['mean_v']],
#                sd_v=lbaPars[['sd_v']])
#   LL <- sum(log(like))
#   if(any(is.na(LL)) | any(is.infinite(LL))) {
#     LL <- min.like
#   }
#   
#   if(returnType=='full') {
#     # pdiffusion requires the RTs to be ordered, so lets order
#     # idx <- order(rt)
#     # PP <- pLBA(rt[order(rt)], response=rep(2, length(choice)), 
#     #            A=lbaPars[['A']][order(rt)], 
#     #            b=lbaPars[['b']][order(rt)], 
#     #            t0=lbaPars[['t0']][order(rt)],
#     #            mean_v=lbaPars[['mean_v']][order(rt),],
#     #            sd_v=lbaPars[['sd_v']][order(rt),])
#     # # and reverse order
#     # PP <- PP[order(order(rt))]
#     # PE_choice = matrix(PE[VVchoiceIdx], ncol=2, byrow=TRUE)
#     return(list(VV=NULL, VV_full=NULL, LL=LL, PE_full=NULL, PE=NULL))
#   } else {
#     return(-LL)
#   }
# }
# 
# objRLLBAMultiCond <- function(pars,
#                             rt, choice, condition, outcomes=NULL, parNames,
#                             values=NULL, VVchoiceIdx=NULL, constants=NULL, 
#                             returnType='full', min.like=-1e20,
#                             backend = 'c') {
#   ### Transcoded from Hanneke den Ouden's Matlab code (2015, see below)
#   
#   ## Learning model is Rescorla-Wagner:
#   # vA <- vA + eta*(r-vA)   -> vA = value of option A, and is updated with predcition error (reward-vA) modulated by learning rate eta
#   
#   ## Choice function is DDM:
#   # a = Threshold
#   # v = drift rate = difference in EV, multiplied by a linear scaling factor
#   # t0 = non-decision time
#   
#   # name parameters, add constants to vector
#   names(pars) <- parNames
#   if(!is.null(constants)) for(i in 1:length(constants)) pars[[names(constants[i])]] <- constants[[i]]
#   etas <- transformLearningRate(condition=condition,
#                                 pars=pars)
# 
#   # Update values trial-by-trial
#   updated <- updateValuesFancy(outcomes, values, etas$eta1, etas$eta2)
#   VV <- updated$VV
#   PE <- updated$PE
#   
#   # On the basis of which values is chosen? Here, we check in 'outcomes' for which alternatives got feedback, and use this to index in VV
#   VV_choice = matrix(VV[VVchoiceIdx], ncol=2, byrow=TRUE)
#   # EV_diff = VV_choice[,1] - VV_choice[,2]  # difference in expected value over trials
#   
#   # DDM choice function
#   lbaPars <- transformLBAPars(pars, condition, ev1=VV_choice[,1], ev2=VV_choice[,2])
#   like <- dLBA(rt=rt, response=choice, 
#                A=lbaPars[['A']], 
#                b=lbaPars[['b']], 
#                t0=lbaPars[['t0']],
#                mean_v=lbaPars[['mean_v']],
#                sd_v=lbaPars[['sd_v']], silent=TRUE)
#   LL <- sum(log(like))
#   if(any(is.na(LL)) | any(is.infinite(LL))) {
#     LL <- min.like
#   }
#   
#   if(returnType=='full') {
#     # pdiffusion requires the RTs to be ordered, so lets order
#     # idx <- order(rt)
#     # PP <- pdiffusion(rt[order(rt)], response=rep(2, length(choice)), 
#     #                  a=ddmPars[['a']][order(rt)], 
#     #                  v=ddmPars[['v']][order(rt)], 
#     #                  t0=ddmPars[['t0']][order(rt)],
#     #                  sz=ddmPars[['sz']][order(rt)],
#     #                  z=ddmPars[['z']][order(rt)],
#     #                  sv=ddmPars[['sz']][order(rt)],
#     #                  s=ddmPars[['s']][order(rt)])
#     # # and reverse order
#     # PP <- PP[order(order(rt))]
#     PE_choice = matrix(PE[VVchoiceIdx], ncol=2, byrow=TRUE)
#     return(list(VV=VV_choice, VV_full=VV, LL=LL, PE_full=PE, PE=PE_choice))
#   } else {
#     return(-LL)
#   }
# }

