rm(list=ls())  # fresh start

# Libraries, directories ---------------------------------------------------------------
.libPaths(c(.libPaths(), '/home/stevenm/rpackages')) # for tux server
library(rtdists)
library(DEoptim)
library(RLDDM)
dataDir <- '~/surfdrive/data/learningTask/Barbara_preprocessed'
workDir <- '/Users/steven/Sync/PhDprojects/RLDDM'

# Model to fit -------------------------------------------------------------------
modelType = 'RLDDM'
source(file.path(workDir, 'analysis/fittingFunctions.R'))
exp <- 'exp3'
resDir <- file.path(workDir, 'fits', 'Barbara', exp, paste0('model-', modelType))
dir.create(resDir, showWarnings = FALSE, recursive=TRUE)

# Optimization options ---------------------------------------------------------
nCores <- 2
# Load data
load(file.path(dataDir, paste0('data_', exp, '.Rdata')))
# very minor changes to data
dat$cue <- as.character(dat$cue)
dat[is.na(dat$cue), 'cue'] <- 'NEU'
dat$cue <- as.factor(dat$cue)
ppsToFit <- unique(dat$pp)
modelsToFit <- 1:4


# Start fitting loop ------------------------------------------------------
for(modelN in modelsToFit) {
  
  # load model specification, and set-up  
  source(file.path(workDir, 'analysis', 'models', paste0('model', modelN, '.R')))
  modelSetup <- model(modelSpec)
  
  # Loop over participants
  for(pp in ppsToFit) {
    # Prepare data for fitting
    d = prepareForFitting(dat[dat$pp==pp&as.character(dat$Block)=='Miniblocks',]) 
    df=d$df
    outcomes=d$outcomes
    VVchoiceIdx=d$VVchoiceIdx
    values=d$values
    choice=d$choice

    # test obj function (no errors, NA, inf...)
    # objRLDDMMultiCond(pars=c(modelSetup$upperBounds-modelSetup$lowerBounds)/2,
    #                   rt=df$rt,
    #                   choice=choice,
    #                   condition=df$cue,
    #                   values=values,
    #                   outcomes=outcomes,
    #                   VVchoiceIdx=VVchoiceIdx,
    #                   returnType='negLL',
    #                   parNames=modelSetup$p.vector, constants=modelSetup$constants)

    # Set-up output file name, check if previous fit exists
    modelFn = file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_Miniblocks.rdat'))
    if(file.exists(modelFn)) {
      load(modelFn)
      if(outDEoptim$optim$iter < 100) {
        # Continue from previous fit
        optimize = TRUE
        start_pop = outDEoptim$member$pop
      } else {
        # start new fit
        optimize = FALSE
        start_pop = NULL
      }
    } else {
      optimize = TRUE
      start_pop = NULL
    }
    if(optimize) {
      # Fit ---------------------------------------------------------------------
      if(nCores > 1) {
        # Manual assignment of nr of cores for parallel
        cl <- parallel::makeCluster(nCores)
        packFn <- function(packages) {
          for (i in packages) library(i, character.only = TRUE)
        }
        parallel::clusterCall(cl, packFn, c('RLDDM', 'rtdists'))
        parallel::clusterExport(cl, c('df', 'objRLDDMMultiCond', 
                                      'outcomes', 
                                      'transformDDMPars', 
                                      'transformLearningRate'))
      } else {
        cl <- NULL
      }
      outDEoptim <- DEoptim(objRLDDMMultiCond, 
                            lower=modelSetup$lowerBounds, 
                            upper=modelSetup$upperBounds, 
                            control=DEoptim.control(itermax=1000, 
                                                    parallelType = ifelse(nCores>1, 1, 0), 
                                                    cluster=cl,
                                                    reltol=100, 
                                                    initialpop = start_pop,
                                                    parVar=c('df', 'objRLDDMMultiCond', 'outcomes', 
                                                             'transformDDMPars', 'transformLearningRate'),
                                                    packages = c('RLDDM', 'rtdists')), 
                            parNames=modelSetup$p.vector, 
                            constants=modelSetup$constants, 
                            choice=choice,
                            rt=df$rt,
                            condition=df$cue,
                            values=values, 
                            outcomes=outcomes, 
                            VVchoiceIdx=VVchoiceIdx,
                            returnType='negLL')
      bestPars <- outDEoptim$optim$bestmem
      
      # save everything
      save(df, outDEoptim, bestPars, modelSetup, transformDDMPars, transformLearningRate, file=modelFn)
    }
  }
}

