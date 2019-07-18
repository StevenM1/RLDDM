rm(list=ls())  # fresh start

# Libraries, directories ---------------------------------------------------------------
source('./setup_locale.R') # gets dataDir, workDir
library(rtdists)
library(DEoptim)
library(RLDDM)
library(moments)

# Model to fit -------------------------------------------------------------------
source(file.path(workDir, 'analysis/src/fittingFunctions.R'))
source(file.path(workDir, 'analysis/src/plottingFunctions.R'))
source(file.path(workDir, 'analysis/src/racingWaldDists.R'))

# Optimization options ---------------------------------------------------------
nCores <- 15
overwritePreviousFit = FALSE
# fit 2 (DDM/RW) x 2 (drift x threshold effect) x 2 (Q-learning/SARSA) = 8 models
# each to Miniblocks and Trialwise
modelsToFit <- c('LR-Q_CF-DDM_SAT-a_st0', 'LR-Q_CF-DDM_SAT-m_st0',
                 'LR-SARSA_CF-DDM_SAT-a_st0', 'LR-SARSA_CF-DDM_SAT-m_st0',
                 'LR-Q_CF-RW-AccAdvantages_SAT-B', 'LR-Q_CF-RW-AccAdvantages_SAT-b0',
                 'LR-SARSA_CF-RW-AccAdvantages_SAT-B', 'LR-SARSA_CF-RW-AccAdvantages_SAT-b0')
exps <- c('exp2')
blocks <- c('Miniblocks', 'Trialwise')

# Function to fit a single sub
fitSingleSubWrapper <- function(pp, dat, modelSetup, modelN, resDir, block,
                                overwritePreviousFit) {
  ## wraps for loading previous fit & saving everything
  thisSubDat <- dat[dat$pp==pp&dat$Block==block,]
  outputFn <- file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_', block, '.rdat'))
  if(file.exists(outputFn) & !overwritePreviousFit) {
    load(outputFn)
    if(!'outDEoptim' %in% ls()) outDEoptim <- NULL
  } else {
    outDEoptim <- NULL
  }
  
  outDEoptim <- fitSingleSub(dat, modelSetup, outDEoptim=outDEoptim, nCores=1)
  bestPars <- outDEoptim$optim$bestmem
  names(bestPars) <- modelSetup$p.vector
  save(df, outDEoptim, bestPars, modelSetup, transformChoicePars, transformLearningRate, file=outputFn)
}

fitSingleSub <- function(dat, modelSetup, outDEoptim=NULL, nCores=1) {
  if(nrow(dat) == 0) {
    stop("Data empty!")
  }
  # Prep data for fitting
  d = prepareForFitting(dat)
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
  
  # Set-up output file name, check if previous fit exists
  if(!is.null(outDEoptim)) {
    if(outDEoptim$optim$iter < 100) {
      # Continue from previous fit
      optimize = TRUE
      initialPop = outDEoptim$member$pop
    } else {
      # start new fit
      optimize = FALSE
      initialPop = NULL
    }
  } else {
    # start new fit
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
  }
  return(outDEoptim)
}

# Start fitting loop ------------------------------------------------------
for(exp in exps) {
  # Load data
  load(file.path(dataDir, paste0('data_', exp, '.Rdata')))
  # very minor changes to data
  dat$cue <- as.character(dat$cue)
  dat[is.na(dat$cue), 'cue'] <- 'NEU'
  dat$cue <- as.factor(dat$cue)
  ppsToFit <- unique(dat$pp)
  
  for(block in blocks) {
    for(modelN in modelsToFit) {
      # set-up output path
      modelType <- sub(".*?CF-(.*?)_.*", "\\1", modelN)
      resDir <- file.path(workDir, 'fits', 'Barbara', exp, paste0('model-', modelType))
      dir.create(resDir, showWarnings = FALSE, recursive=TRUE)
      
      # load model specification, and set-up  
      source(file.path(workDir, 'analysis', 'models', paste0(modelN, '.R')))
      modelSetup <- model(modelSpec)
      
      if(nCores > 1) {
        library(snowfall)
        sfInit(cpus = nCores, parallel = TRUE)
        sfExportAll()
        for(package.name in c('rtdists', 'RLDDM', 'SuppDists', 'DEoptim')) sfLibrary(package.name, lib.loc='/home/stevenm/rpackages', character.only = TRUE)
        sfLapply(ppsToFit, function(x) fitSingleSubWrapper(pp = x, 
                                                           dat = dat, 
                                                           block = block, 
                                                           modelSetup = modelSetup,
                                                           modelN = modelN, 
                                                           resDir = resDir,
                                                           overwritePreviousFit = overwritePreviousFit))
        sfStop()
      } else {
        lapply(ppsToFit, function(x) fitSingleSubWrapper(pp=x, 
                                                         dat=dat, 
                                                         block=block, 
                                                         modelSetup=modelSetup,
                                                         modelN=modelN, resDir = resDir,
                                                         overwritePreviousFit=overwritePreviousFit))
      }
    }
  }
}



# 
# # Loop over participants
# for(pp in ppsToFit) {
#   print(pp)
#   # Prepare data for fitting
#   thisPpDat <- dat[dat$pp==pp&as.character(dat$Block)==block,]
#   if(nrow(thisPpDat) == 0) next
#   d = prepareForFitting(thisPpDat)
#   df=d$df
#   outcomes=d$outcomes
#   VVchoiceIdx=d$VVchoiceIdx
#   values=d$values
#   choice=d$choice
#   choice = ifelse(choice==2, 1, 2)
#   
#   # test obj function (no errors, NA, inf...)
#   obj(pars=c(modelSetup$upperBounds-modelSetup$lowerBounds)/2,
#             rt=df$rt,
#             choice=choice,
#             condition=df$cue,
#             values=values,
#             outcomes=outcomes,
#             VVchoiceIdx=VVchoiceIdx,
#             returnType='negLL',
#             parNames=modelSetup$p.vector,
#             constants=modelSetup$constants,
#             learningRule=modelSetup$learningRule,
#             choiceFunction=modelSetup$choiceFunction)
#   
#   # Set-up output file name, check if previous fit exists
#   modelFn = file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_', block, '.rdat'))
#   if(file.exists(modelFn) & !overwritePreviousFit) {
#     load(modelFn)
#     if(outDEoptim$optim$iter < 100) {
#       # Continue from previous fit
#       optimize = TRUE
#       initialPop = outDEoptim$member$pop
#     } else {
#       # start new fit? should be done now
#       optimize = FALSE
#       initialPop = NULL
#     }
#   } else {
#     optimize = TRUE
#     initialPop = NULL
#   }
#   if(overwritePreviousFit) {
#     optimize = TRUE
#     initialPop = NULL
#   }
#   if(optimize) {
#     # Fit ---------------------------------------------------------------------
#     if(nCores > 1) {
#       # Manual assignment of nr of cores for parallel
#       cl <- parallel::makeCluster(nCores)
#       packFn <- function(packages) {
#         for (i in packages) library(i, character.only = TRUE, lib.loc='/home/stevenm/rpackages')
#       }
#       parallel::clusterCall(cl, packFn, c('RLDDM', 'rtdists', 'SuppDists'))
#       parallel::clusterExport(cl, c('df', 
#                                     'obj', 
#                                     'outcomes', 
#                                     'transformChoicePars', 
#                                     'transformLearningRate',
#                                     'n1Wald', 'rWald', 'pWald', 'dWald', 'dWaldRace',
#                                     'check_i_arguments', 'check_n1_arguments'))
#     } else {
#       cl <- NULL
#     }
#     outDEoptim <- DEoptim(obj, 
#                           lower=modelSetup$lowerBounds, 
#                           upper=modelSetup$upperBounds, 
#                           control=DEoptim.control(itermax=5000, 
#                                                   parallelType = ifelse(nCores>1, 1, 0), 
#                                                   cluster=cl,
#                                                   step=500, 
#                                                   initialpop = initialPop,
#                                                   packages = c('RLDDM', 'rtdists')
#                                                   ), 
#                           parNames=modelSetup$p.vector, 
#                           constants=modelSetup$constants, 
#                           learningRule=modelSetup$learningRule,
#                           choice=choice,
#                           rt=df$rt,
#                           condition=df$cue,
#                           values=values, 
#                           outcomes=outcomes, 
#                           VVchoiceIdx=VVchoiceIdx,
#                           returnType='negLL',
#                           choiceFunction=modelSetup$choiceFunction)
#     bestPars <- outDEoptim$optim$bestmem
#     
#     # save everything
#     save(df, outDEoptim, bestPars, modelSetup, transformChoicePars, transformLearningRate, file=modelFn)
#     
#     # plot
#     # nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
#     # par(mar=c(5, 4, 2, 2) + 0.1)
#     # tmp = plotSingleSub()
#   }
# }