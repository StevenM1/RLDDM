rm(list=ls())
source('./src/racingWaldDists.R')

##### Simple parameter recovery -----------------------------------------------
# ### No trial-by-trial estimation
# obj <- function(pars, rt, response, parNames, constants=NULL, min.like=-1e20) {
#   names(pars) <- parNames
#   if(!is.null(constants)) for(i in 1:length(constants)) pars[[names(constants[i])]] <- constants[[i]]
# 
#   LL <- dWaldRace(rt=rt, response=response,
#                   v=c(pars[['v1']], pars[['v2']]),
#                   B=pars[['B']],
#                   A=pars[['A']], t0=pars[['t0']], st0=0, silent=TRUE)
#   LL <- sum(log(LL))
#   if(any(is.na(LL)) | any(is.infinite(LL))) {
#     LL <- min.like
#   }
#   -LL
# }
# 
# library(DEoptim)
# lower <- c(v1=0, v2=0, B=1e-3, t0=1e-3, A=1e-3)
# upper <- c(v1=5, v2=5, B=5, t0=1, A=5)
# constants <- c()
# 
# nSets <- 100
# results <- data.frame(v1=rep(NA, nSets), v2=NA, B=NA, t0=NA, A=NA,
#                       v1_fit=NA, v2_fit=NA, B_fit=NA, t0_fit=NA, A_fit=NA)
# for(i in 1:nSets) {
#   print(i)
#   while(TRUE) {
#     v1 <- runif(1, 0, 5)
#     v2 <- runif(1, 0, 5)
#     t0 <- runif(1, 0, .5)
#     B <- runif(1, 1e-4, 5)
#     A <- runif(1, 1e-4, 5)
#     dat <- rWaldRace(n=500, v=c(v1, v2), B=B, A=A, t0=t0)
# 
#     if(mean(dat$RT) > 3) next
#     if(mean(dat$RT) < .2) next
#     if(mean(dat$R==1) < .05) next
#     if(mean(dat$R==2) < .05) next
#     break
#   }
#   results[i, c('v1', 'v2', 'B', 't0', 'A')] <- c(v1, v2, B, t0, A)
#   #debug(n1Wald)
#   outDEoptim <- DEoptim(fn=obj, lower=lower, upper=upper,
#                         control=DEoptim.control(itermax=500, steptol=100, parallelType=1,
#                                                 parVar=c('dWaldRace', 'check_i_arguments', 'check_n1_arguments',
#                                                          'dWald', 'pWald', 'n1Wald')),
#                         rt=dat$RT,
#                         response=dat$R,
#                         parNames=names(lower),
#                         constants=constants)
#   bestPars <- outDEoptim$optim$bestmem
#   names(bestPars) <- names(lower)
# 
#   results[i, c('v1', 'v2', 'B', 't0', 'A',
#                'v1_fit', 'v2_fit', 'B_fit', 't0_fit', 'A_fit')] <- c(v1, v2, B, t0, A,
#                                                                      bestPars[['v1']],
#                                                                      bestPars[['v2']],
#                                                                      bestPars[['B']],
#                                                                      bestPars[['t0']],
#                                                                      bestPars[['A']])
#   par(mfrow=c(2,3))
#   for(param in c('v1', 'v2', 'B', 't0', 'A')) {
#     plot(results[,param], results[,paste0(param, '_fit')], main=param,
#          xlim=c(0, 5), ylim=c(0, 5), xlab=paste('True', param), ylab=paste('Fitted', param))
#     abline(a=0, b=1)
#   }
#   
#   # plot fit
#   d1 <- density(dat$RT[dat$R==1])
#   d2 <- density(dat$RT[dat$R==2])
#   d1$y <- d1$y*mean(dat$R==1)
#   d2$y <- d2$y*mean(dat$R==2)
#   plot(d1$x, d1$y, lwd=2, type='l')
#   lines(d2$x, d2$y, col=2, lwd=2)
#   
#   xs <- seq(0, 10, .05)
#   lines(xs, dWaldRace(xs, response=1, A=bestPars[['A']], B=bestPars[['B']], v=c(bestPars[['v1']], bestPars[['v2']]), t0=bestPars[['t0']]),
#         lty=2)
#   lines(xs, dWaldRace(xs, response=2, A=bestPars[['A']], B=bestPars[['B']], v=c(bestPars[['v1']], bestPars[['v2']]), t0=bestPars[['t0']]),
#         col=2, lty=2)
# }
# #
# par(mfrow=c(2,3))
# for(param in c('v1', 'v2', 'B', 't0', 'A')) {
#   plot(results[,param], results[,paste0(param, '_fit')], main=param, xlab=paste('True', param), ylab=paste('Fitted', param))
# }


##### Parameter recovery of trial-by-trial parameters ----------
makeStimuli <- function(nTrialsPerStimSet, stimRewardProbs=list(c(.7, .3), c(.7, .3), c(.7, .3))) {
  
  df <- NULL
  stimSet <- 1
  for(stimRewardProb in stimRewardProbs) {
    thisDf <- data.frame(trial=1:nTrialsPerStimSet,
                         rewardProb1=stimRewardProb[1],
                         rewardProb2=stimRewardProb[2],
                         stimulus_set=stimSet)
    stimSet <- stimSet + 1
    df <- rbind(df, thisDf)
  }
  # shuffle
  df <- df[sample(1:nrow(df)),]
  df$trial <- 1:nrow(df)
  df
}

simRLRW <- function(nTrials, B, t0, eta1, eta2=eta1, b0, b1, bU=0, bB=0, A=0, stimuli=NULL, learningRule='SARSA') {
  # Simulates reinforcement learning with Racing Walds as a decision model
  if(is.null(stimuli)) stimuli <- makeStimuli(nTrials/3)
  
  stimuli$choice=NA
  stimuli$rt <- NA
  nStimSets <- length(unique(stimuli$stimulus_set))
  
  # pre-allocate some memory
  B_tr = B - (1:nTrials)*bB
  vs <- c()
  values <- matrix(NA, nrow=nTrials+1, ncol=2*nStimSets)
  values[1,] <- 0.5
  rewards <- matrix(NA, nrow=nTrials, ncol=1)
  PEs <- c()
  
  # loop
  for(tr in 1:nTrials) {
    i = stimuli[tr,'stimulus_set']
    stimulus2 <- c(i:(i+1)+1*(i-1))[2]
    stimulus1 <- c(i:(i+1)+1*(i-1))[1]
    v1 <- b0 + bU*tr + b1*(values[tr,stimulus1] - values[tr,stimulus2])
    v2 <- b0 + bU*tr + b1*(values[tr,stimulus2] - values[tr,stimulus1])
    v <- c(v1, v2)
    vs <- rbind(vs, v)
    
    # simulate diffusion
    rtChoice <- rWaldRaceSM(n=1, B=B_tr[tr], v=v, t0=t0, A=A)
    stimuli[tr, 'choice'] <- choice <- as.numeric(rtChoice$R)[1]
    stimuli[tr, 'rt'] <- rtChoice$RT[1]
    
    # Reward?
    rewards[tr] <- ifelse(runif(1) < stimuli[tr, c('rewardProb1', 'rewardProb2')][choice], 1, 0)
    
    # update values
    chosenStim <- c(i:(i+1)+1*(i-1))[-choice]
    if(learningRule == 'SARSA') {
      PE <- as.numeric(rewards[tr] - values[tr, chosenStim])  # PE = outcome - expected value of chosen action
    } else if(learningRule == 'Qlearning') {
      PE <- as.numeric(rewards[tr] - max(values[tr, c(stimulus1, stimulus2)]))  # PE = outcome - max expected value
    }
    values[tr+1, chosenStim] <- values[tr, chosenStim] + ifelse(PE>0, PE*eta1, PE*eta2)  # TD learning
    values[tr+1, -chosenStim] <- values[tr, -chosenStim]  # value of not-chosen stimulus remains the same
    PEs <- c(PEs, PE)
  }
  
  values <- data.frame(values)
  colnames(values) <- c('value1', 'value2')
  stimuli <- cbind(stimuli, values[1:nTrials,])
  stimuli$reward <- rewards
  stimuli$PE <- PEs
  stimuli$v <- vs
  stimuli$B <- B_tr
  stimuli
}

simRLRWSingleSub <- function(modelSetup) {
  # Simulate, check out data
  while(TRUE) {
    truePars <- c()
    for(pNum in 1:length(modelSetup$p.vector)) {
      truePars <- c(truePars, runif(1, modelSetup$lowerBounds[pNum], modelSetup$upperBounds[pNum]))
    }
    names(truePars) <- modelSetup$p.vector
    truePars[['b1']] <- rbeta(1, 1.2, 6)*100
    truePars[['eta1']] <- rbeta(1, 1.2, 6)
    
    if('bB' %in% modelSetup$p.vector) {
      if(truePars[['B']]-truePars[['bB']]*nTrials <= 0.001) next
    }
    
    simDat <- simRLRW(nTrials=300, 
                      B=truePars[['B']], t0=truePars[['t0']],
                      eta1=truePars[['eta1']], b0=truePars[['b0']],
                      b1=-truePars[['b1']], 
                      bU=ifelse('bU' %in% modelSetup$p.vector, truePars[['bU']], 0),
                      bB=ifelse('bB' %in% modelSetup$p.vector, truePars[['bB']], 0),
                      learningRule = modelSetup$learningRule)
    simDat$cue <- 1
    simDat$choiceIsHighP <- simDat$choice == 1
    if(mean(simDat$rt) > 2) next
    if(mean(simDat$rt) < .3) next
    if(mean(simDat$choiceIsHighP) < .55) next
    if(mean(simDat$choiceIsHighP) > .90) next
    return(list(simDat=simDat, truePars=truePars))
  }
}

# Function to fit a single sub
fitSingleSub <- function(pp, dat, modelSetup, outputFn, truePars,
                         overwritePreviousFit=FALSE,
                         nCores=1) {
  print(pp)
  # ugly
  thisSubDat <- dat
  if(nrow(thisSubDat) == 0) {
    stop("Data empty!")
  }
  # Prep data for fitting
  d = prepareForFitting(thisSubDat)
  df=d$df
  outcomes=d$outcomes
  VVchoiceIdx=d$VVchoiceIdx
  values=d$values
  choice=d$choice
  # NB: choice==2 is "Correct", choice==1 is "Wrong". This is the "ddiffusion convention" in rtdists
  
  # test obj function (no errors, NA, inf...)
  obj(pars=c(modelSetup$upperBounds-modelSetup$lowerBounds)/2,
      rt=df$rt,
      choice=choice,
      condition=df$cue,
      values=values,
      outcomes=outcomes,
      VVchoiceIdx=VVchoiceIdx,
      returnType='negLL',
      parNames=modelSetup$p.vector,
      constants=modelSetup$constants,
      learningRule=modelSetup$learningRule,
      choiceFunction=modelSetup$choiceFunction)
  
  if(file.exists(outputFn) & !overwritePreviousFit) {
    load(outputFn)
    if(outDEoptim$optim$iter < 100) {
      # Continue from previous fit
      optimize = TRUE
      initialPop = outDEoptim$member$pop
    } else {
      # start new fit? should be done now
      optimize = FALSE
      initialPop = NULL
    }
  } else {
    optimize = TRUE
    initialPop = NULL
  }
  if(overwritePreviousFit) {
    optimize = TRUE
    initialPop = NULL
  }
  if(optimize) {
    # Fit ---------------------------------------------------------------------
    if(nCores > 1) {
      # Manual assignment of nr of cores for parallel
      cl <- parallel::makeCluster(nCores)
      packFn <- function(packages) {
        for (i in packages) library(i, character.only = TRUE, lib.loc='/home/stevenm/rpackages')
      }
      parallel::clusterCall(cl, packFn, c('RLDDM', 'rtdists', 'SuppDists'))
      parallel::clusterExport(cl, c('df', 
                                    'obj', 
                                    'outcomes', 
                                    'transformChoicePars', 
                                    'transformLearningRate',
                                    'n1Wald', 'rWald', 'pWald', 'dWald', 'dWaldRace',
                                    'check_i_arguments', 'check_n1_arguments'))
    } else {
      cl <- NULL
    }
    outDEoptim <- DEoptim(obj, 
                          lower=modelSetup$lowerBounds, 
                          upper=modelSetup$upperBounds, 
                          control=DEoptim.control(itermax=5000, 
                                                  parallelType = ifelse(nCores>1, 1, 0), 
                                                  cluster=cl,
                                                  step=500, 
                                                  initialpop = initialPop,
                                                  packages = c('RLDDM', 'rtdists')
                          ), 
                          parNames=modelSetup$p.vector, 
                          constants=modelSetup$constants, 
                          learningRule=modelSetup$learningRule,
                          choice=choice,
                          rt=df$rt,
                          condition=df$cue,
                          values=values, 
                          outcomes=outcomes, 
                          VVchoiceIdx=VVchoiceIdx,
                          returnType='negLL',
                          choiceFunction=modelSetup$choiceFunction)
    bestPars <- outDEoptim$optim$bestmem
    
    # save everything
    save(df, outDEoptim, bestPars, truePars, modelSetup, transformChoicePars, transformLearningRate, file=outputFn)
    
    # plot
    # nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
    # par(mar=c(5, 4, 2, 2) + 0.1)
    # tmp = plotSingleSub()
  }
}

#
doSingleParRec <- function(i, modelSetup, parRecDir, modelN, overwrite) {
  out <- simRLRWSingleSub(modelSetup = modelSetup)
  truePars <- out$truePars
  simDat <- out$simDat
  simDat$outcome = simDat$reward*100
  outputFn <- file.path(parRecDir, paste0('model-', modelN, '_dataset-', i, '.rdat'))
  
  fitSingleSub(i, simDat, modelSetup, outputFn, truePars=truePars, overwritePreviousFit = overwrite)
}

# Libraries, directories ---------------------------------------------------------------
source('./setup_locale.R') # gets dataDir, workDir
library(rtdists)
library(DEoptim)
library(RLDDM)
source(file.path(workDir, 'analysis/src/fittingFunctions.R'))
# for fitting
parRecDir <- file.path(workDir, 'parRec')
nDataSets <- 150
modelType = 'RLRW'
nCores <- 20
overwrite = TRUE
nTrials <- 300
maxIter <- 5000

# for plotting
source(file.path(workDir, 'analysis/src/plottingFunctions.R'))
library(snowfall)

modelNs <- c('LR-SARSA_CF-RW-AccAdvantages_SAT-NULL')
for(modelN in modelNs) {
  source(file.path(workDir, 'analysis', 'models', paste0(modelN, '.R')))
  modelSetup <- model(modelSpec)
  
  #
  if(nCores > 1) {
    library(snowfall)
    sfInit(cpus = nCores, parallel = TRUE)
    sfExportAll()
    for(package.name in c('rtdists', 'RLDDM', 'SuppDists', 'DEoptim')) sfLibrary(package.name, lib.loc='/home/stevenm/rpackages', character.only = TRUE)
    sfLapply(1:nDataSets, function(x) doSingleParRec(i=x, 
                                                     modelSetup=modelSetup, 
                                                     parRecDir=parRecDir,
                                                     modelN=modelN,
                                                     overwrite=overwrite))
    sfStop()
  } else {
    lapply(1:nDataSets, function(x) doSingleParRec(i=x, 
                                                   modelSetup=modelSetup, 
                                                   parRecDir=parRecDir,
                                                   modelN=modelN,
                                                   overwrite=overwrite))
  }
}


# plot overall
library(stringr)
modelN <- 'LR-SARSA_CF-RW-AccAdvantages_SAT-NULL'
source(file.path(workDir, 'analysis', 'models', paste0(modelN, '.R')))
plotSingleSubs <- FALSE
dirs <- list.files(parRecDir, pattern=paste0('model-', modelN, '_dataset'))
dirs
parDf <- NULL
for(dir in dirs) {
  if('truePars' %in% ls())  rm(truePars)
  load(file.path(parRecDir, dir))
  if(!'truePars' %in% ls()) next
  names(bestPars) <- paste0(modelSetup$p.vector, '.fit')
  pp <- sub('.*?dataset-(.*?).rdat', '\\1', dir)
  
  if(plotSingleSubs) {
    nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
    par(mar=c(5, 4, 2, 2) + 0.1)
    # prep data
    d = prepareForFitting(df)
    df=d$df
    outcomes=d$outcomes
    VVchoiceIdx=d$VVchoiceIdx
    values=d$values
    choice=d$choice
    
    # plot, get quantiles
    out = plotSingleSub()
  }
  parDf <- rbind(parDf, c('dataset'=pp, bestPars, truePars))
}
parDf <- data.frame(parDf)
for(col in colnames(parDf)) {
  parDf[,col] <- as.numeric(as.character(parDf[,col]))
}

par(mfrow=c(2,3))
for(pNum in 1:length(modelSetup$p.vector)) {
  pname = modelSetup$p.vector[pNum]
  x = parDf[,pname]
  y = parDf[,paste0(pname, '.fit')]
  if(pname == 'b1') {
    xlim <- ylim <- c(-2, 100)
  } else {
    xlim = ylim = c(modelSetup$lowerBounds[pNum], modelSetup$upperBounds[pNum])
  }
  plot(x, y, main=pname, xlab='True', ylab='Fit', xlim=xlim, ylim=ylim)
  abline(a=0, b=1, col='grey')
  abline(reg=lm(y~x), lty=2)
  corcoef <- round(cor(x, y), 2)
  rmse <- round(sqrt(mean((y-x)^2)),3)
  legend('bottomright', paste0('r = ', corcoef, '\nrmse = ', rmse), bty='n')
}





# # Simulate, check out data
# while(TRUE) {
#   truePars <- c()
#   for(pNum in 1:length(modelSetup$p.vector)) {
#     truePars <- c(truePars, runif(1, modelSetup$lowerBounds[pNum], modelSetup$upperBounds[pNum]))
#   }
#   names(truePars) <- modelSetup$p.vector
#   truePars[['b1']] <- rbeta(1, 1.2, 6)*100
#   truePars[['eta1']] <- rbeta(1, 1.2, 6)
#   
#   if('bB' %in% modelSetup$p.vector) {
#     if(truePars[['B']]-truePars[['bB']]*nTrials <= 0.001) next
#   }
#   
#   simDat <- simRLRW(nTrials=300, 
#                     B=truePars[['B']], t0=truePars[['t0']],
#                     eta1=truePars[['eta1']], b0=truePars[['b0']],
#                     b1=-truePars[['b1']], 
#                     bU=ifelse('bU' %in% modelSetup$p.vector, truePars[['bU']], 0),
#                     bB=ifelse('bB' %in% modelSetup$p.vector, truePars[['bB']], 0))
#   simDat$cue <- 1
#   simDat$choiceIsHighP <- simDat$choice == 1
#   if(mean(simDat$rt) > 2) next
#   if(mean(simDat$rt) < .3) next
#   if(mean(simDat$choiceIsHighP) < .55) next
#   if(mean(simDat$choiceIsHighP) > .90) next
#   break
# }

# 
# for(i in 1:nDataSets) {
#   
#   modelFn <- file.path(parRecDir, paste0('model-', modelN, '_dataset-', i, '.rdat'))
#   if(file.exists(modelFn) & !overwrite) {
#     load(modelFn)
#     if(outDEoptim$optim$iter < maxIter) next
#   }
#   
#   
#   out <- simRLRWSingleSub(modelSetup)
#   simDat <- out$simDat
#   truePars <- out$truePars

#      
#     # Prepare data for fitting
#     simDat$outcome = simDat$reward*100
#     thisPpDat <- simDat
#     d = prepareForFitting(thisPpDat)
#     df=d$df
#     outcomes=d$outcomes
#     VVchoiceIdx=d$VVchoiceIdx
#     values=d$values
#     choice=d$choice
# #    choice = ifelse(choice==2, 1, 2)
#     
#     # test obj function (no errors, NA, inf...)
#     obj(pars=truePars,
#         rt=df$rt,
#         choice=choice,
#         condition=df$cue,
#         values=values,
#         outcomes=outcomes,
#         VVchoiceIdx=VVchoiceIdx,
#         returnType='negLL',
#         parNames=names(truePars), 
#         constants=modelSetup$constants,
#         choiceFunction=modelSetup$choiceFunction,
#         learningRule=modelSetup$learningRule)
#     
#     optimize = TRUE
#     initialPop = NULL
#     if(optimize) {
#       # Fit ---------------------------------------------------------------------
#       if(nCores > 1) {
#         # Manual assignment of nr of cores for parallel
#         cl <- parallel::makeCluster(nCores)
#         packFn <- function(packages) {
#           for (i in packages) library(i, character.only = TRUE, lib.loc='/home/stevenm/rpackages')
#         }
#         parallel::clusterCall(cl, packFn, c('RLDDM', 'rtdists'))
#         parallel::clusterExport(cl, c('df',
#                                       'obj',
#                                       'outcomes',
#                                       'transformChoicePars',
#                                       'transformLearningRate',
#                                       'n1Wald', 'rWald', 'pWald', 'dWaldRace',
#                                       'dWald', 'check_i_arguments', 'check_n1_arguments'))
#       } else {
#         cl <- NULL
#       }
#       outDEoptim <- DEoptim(obj,
#                             lower=modelSetup$lowerBounds,
#                             upper=modelSetup$upperBounds,
#                             control=DEoptim.control(itermax=maxIter,
#                                                     parallelType = ifelse(nCores>1, 1, 0),
#                                                     cluster=cl,
#                                                     step=100,
#                                                     initialpop = initialPop),
#                             parNames=modelSetup$p.vector,
#                             constants=modelSetup$constants,
#                             choice=choice,
#                             rt=df$rt,
#                             condition=df$cue,
#                             values=values,
#                             outcomes=outcomes,
#                             VVchoiceIdx=VVchoiceIdx,
#                             returnType='negLL')
#       bestPars <- outDEoptim$optim$bestmem
#       
#       # save everything
#       save(d, truePars, outDEoptim, bestPars, modelSetup, file=modelFn)
#     }
#     
#     # plot
#     nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
#     par(mar=c(5, 4, 2, 2) + 0.1)
#     pp <- i
#     tmp = plotSingleSub()
#   }
# }
