rm(list=ls())  # fresh start

# Libraries ---------------------------------------------------------------
.libPaths(c(.libPaths(), '/home/stevenm/rpackages'))
library(rtdists)
library(DEoptim)
library(RLDDM)

# Directories -------------------------------------------------------------
dataDir <- '~/surfdrive/Testfolder/exp2'
prepDataDir <- '~/surfdrive/data/learningTask/Barbara_preprocessed'
workDir <- '/Users/steven/Sync/PhDprojects/RLDDMproj'
source(file.path(workDir, 'analysis/fittingFunctions.R'))


# Model to fit -------------------------------------------------------------------
modelType = 'RLDDM'
modelN <- 1
exp <- 'exp3'
resDir <- file.path(workDir, 'fits', 'Barbara', exp, paste0('model-', modelType))
dir.create(resDir, showWarnings = FALSE, recursive=TRUE)

getModelSpec <- function(modelN) {
  if(modelN == 1) {
    modelSpec = list('variablePars'=list('a' = 'condition',
                                         'm' = 1,
                                         't0' = 1,
                                         'eta1' = 1),
                     'constants'=c('sv'=0, 'sz'=0, 'z'=0.5, 's'=1, 'eta2'=-Inf),
                     'condition'=c('SPD', 'ACC'))
  } else if(modelN == 2) {
    modelSpec = list('variablePars'=list('a' = 1,
                                         'm' = 1,
                                         't0' = 1,
                                         'eta1' = 'condition'),
                     'constants'=c('sv'=0, 'sz'=0, 'z'=0.5, 's'=1, 'eta2'=-Inf),
                     'condition'=c('SPD', 'ACC'))
  } else if(modelN == 3) {
    modelSpec = list('variablePars'=list('a' = 1,
                                         'm' = 1,
                                         't0' = 1,
                                         'eta1' = 1),
                     'constants'=c('sv'=0, 'sz'=0, 'z'=0.5, 's'=1, 'eta2'=-Inf),
                     'condition'=c('SPD', 'ACC'))
  } else if(modelN == 4) {
    modelSpec = list('variablePars'=list('a' = 'condition',
                                         'm' = 1,
                                         't0' = 1,
                                         'eta1' = 'condition'),
                     'constants'=c('sv'=0, 'sz'=0, 'z'=0.5, 's'=1, 'eta2'=-Inf),
                     'condition'=c('SPD', 'ACC'))
  }
  modelSpec
}

# Load data
load(file.path(prepDataDir, paste0('data_', exp, '.Rdata')))

# Fitting routine
nCores <- 2
allPps <- unique(dat$pp)
#pp <- allPps[1]
for(modelN in 1:4) {
  for(pp in allPps) {
    d = prepareForFitting(dat[dat$pp==pp&as.character(dat$Block)=='Miniblocks',]) 
    df=d$df
    outcomes=d$outcomes
    VVchoiceIdx=d$VVchoiceIdx
    values=d$values
    choice=d$choice
    etaMax <- 1
    
    df$cue <- as.character(df$cue)
    df[is.na(df$cue), 'cue'] <- 'NEU'
    df$cue <- as.factor(df$cue)
    
    modelSetup = setupModel(getModelSpec(modelN=modelN))
    
    # test obj function
    objRLDDMMultiCond(pars=c(modelSetup$upperBounds-modelSetup$lowerBounds)/2,
                      rt=df$rt,
                      choice=choice,
                      condition=df$cue,
                      values=values,
                      outcomes=outcomes,
                      VVchoiceIdx=VVchoiceIdx,
                      modelN=modelN,
                      returnType='negLL',
                      parNames=modelSetup$p.vector, constants=modelSetup$constants)

    modelFn = file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_Miniblocks.rdat'))
    if(file.exists(modelFn)) {
      load(modelFn)
      if(outDEoptim$optim$iter < 100) {
        optimize = TRUE
        start_pop = outDEoptim$member$pop
      } else {
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
        cl <- parallel::makeCluster(12)
        packFn <- function(packages) {
          for (i in packages) library(i, character.only = TRUE)
        }
        parallel::clusterCall(cl, packFn, c('RLDDM', 'rtdists'))
        parallel::clusterExport(cl, c('df', 'objRLDDMMultiCond', 
                                      'outcomes', 'transformDDMPars', 
                                      'transformLearningRate'))
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
                            modelN=modelN,
                            returnType='negLL')
      bestPars <- outDEoptim$optim$bestmem
      # stopCluster(cl)
      save(df, outDEoptim, bestPars, file=modelFn)
    }
    
    # Plot learning -----------------------------------------------------------
    # if(savePlots) pdf(file=file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_fit.pdf')), width=15, height=15)
    
    # par(mfrow=c(1,1), las=1)
    # plotPars <- list(bestPars)[[1]]
    # model <- objRLDDMMultiCond(plotPars, 
    #                            rt=df$rt,
    #                            parNames=modelSetup$parNames,
    #                            fixedPars=modelSetup$fixedPars,
    #                            choice=choice,
    #                            condition=df$cue,
    #                            values=values, 
    #                            outcomes=outcomes, 
    #                            VVchoiceIdx=VVchoiceIdx,
    #                            modelN=modelN,
    #                            returnType='full')
    # plotVV <- model$VV
    # par(mfcol=c(3,2), mar=c(3,4,3,2)+.1, lwd=2, cex=1.05)
    
    # palette(c('cornflowerblue', 'orchid'))
    # for(i in 1:length(unique(df$condition))) {
    #   idx <- df$condition==unique(df$condition)[i]
    #   # Value difference
    #   plotTrace(model$VV[idx,1]-model$VV[idx,2], 
    #             add=FALSE, ylab='Value difference', xlab='Trial n')
    #   abline(h=0, lty=2, col='grey')
    #   points(1:sum(idx), ifelse(df$choice[idx]==2, par()$usr[4], par()$usr[3]), pch=4, xpd=TRUE)
    #   
    #   # Prediction errors
    #   pes <- apply(model$PE, 1, sum, na.rm=TRUE)[idx]
    #   plotTrace(pes, 
    #             add=FALSE, ylab='Prediction errors', xlab='Trial n')
    #   plotTrace(pes)#*ifelse(pes>0, bestPars[['eta1']], bestPars[['eta2']]), add=TRUE, lwd=2, lty=2, col=2)
    #   abline(h=0, lty=2, col='grey')
    #   points(1:sum(idx), ifelse(df$choice[idx]==2, par()$usr[4], par()$usr[3]), pch=4, xpd=TRUE)
    #   
    #   # Probability of choosing upper bound
    #   plotTrace(model$PP[idx,1], ylim=c(0, 1),
    #             add=FALSE, ylab='Probability of "correct"', xlab='Trial n')
    #   abline(h=0.5, lty=2, col='grey')
    #   points(1:sum(idx), ifelse(df$choice[idx]==2, par()$usr[4], par()$usr[3]), pch=4, xpd=TRUE)
    #   plotRunningWindowChoice2(df$choice[idx]-1, binSize=10, lty=2, col='darkred', lwd=2)
    #   
    #   title(formatStr(c('eta1', 'eta2', 'beta'), plotPars, n.decimals=2), outer=TRUE, line=-1)
    # }
    # if(savePlots) dev.off()
  }
}


# Collect all parameters --------------------------------------------------
## Get all parameters, mean RT and mean accuracy into single dataframe
fit <- data.frame(expand.grid(pp=allPps, modelN=1:4, cue=c('SPD', 'ACC')),
                  eta1=NA, eta2=NA, a=NA, t0=NA, m=NA, mrt=NA, accObjective=NA, accSubjective=NA, 
                  negLL=NA, nPars=NA, nObs=NA, BIC=NA, minBIC=NA, wBIC=NA, AIC=NA, minAIC=NA, wAIC=NA)
# allData <- NULL
# keepCols <- c('TrialNumber', "stim_right", 'stim_left', 'condition', 'p_win_correct', 'p_win_incorrect',
#               'rt', 'choiceSymbol', 'choiceIsHighP', 'outcome', 'pp', 'task', 'EVchoice', 'eta1', 'beta',
#               'PE', 'update')# 'eta2', 'beta')

for(modelN in 1:4) {
  for(pp in allPps) {
    # load data
    d = prepareForFitting(dat[dat$pp==pp&as.character(dat$Block)=='Miniblocks',]) 
    #    fnBase = d$fnBase
    df=d$df
    outcomes=d$outcomes
    VVchoiceIdx=d$VVchoiceIdx
    values=d$values
    choice=d$choice
    etaMax <- 1
    
    df$cue <- as.character(df$cue)
    df[is.na(df$cue), 'cue'] <- 'NEU'
    df$cue <- as.factor(df$cue)
    
    # get parameter names & bounds for current model
    modelSetup = setupModel(getModelSpec(modelN=modelN))
    
    # load model
    modelFn = file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_Miniblocks.rdat'))
    load(modelFn)
    names(bestPars) <- modelSetup$parNames
    model <- objRLDDMMultiCond(bestPars,
                               rt=df$rt,
                               parNames=modelSetup$parNames,
                               constants=modelSetup$constants,
                               choice=choice,
                               condition=df$cue,
                               values=values,
                               outcomes=outcomes,
                               VVchoiceIdx=VVchoiceIdx,
                               modelN=modelN,
                               returnType='full')
    df$EVchoice <- model$VV[,1]-model$VV[,2]
    df$modelChoice <- ifelse(df$EVchoice > 0, 2, 1)
    # "subjective" accuracy: choice in line with highest EV-stimulus
    df$accuracySubjective <- df$modelChoice == df$choice
    # "objective" accuracy: Choice is optimal
    df$accuracyObjective <- df$choice == 2
    df$PE <- apply(model$PE, 1, sum, na.rm=TRUE)
    
    # some metrics
    mrt <- aggregate(rt ~ cue, df, mean)
    accObjective <- aggregate(accuracyObjective ~ cue, df, mean)
    accSubjective <- aggregate(accuracySubjective ~ cue, df, mean)
    idx1 <- fit$pp==pp & fit$modelN==modelN
    fit[idx1, 'nObs'] <- nrow(df)
    for(cue in mrt$cue) {
      fit[idx1&fit$cue==cue, 'mrt'] <- mrt[mrt$cue==cue, 'rt']
      fit[idx1&fit$cue==cue, 'accObjective'] <- accObjective[accObjective$cue==cue, 'accuracyObjective']
      fit[idx1&fit$cue==cue, 'accSubjective'] <- accSubjective[accSubjective$cue==cue, 'accuracySubjective']
    }
    
    fit[idx1, 'negLL'] <- -model$LL
    fit[idx1, 't0'] <- bestPars[['t0']]
    fit[idx1, 'm'] <- bestPars[['m']]
    fit[idx1, 'nPars'] <- length(bestPars)
    if(modelN == 1) {
      fit[idx1, 'eta1'] <- bestPars[['eta1']]
      fit[idx1&fit$cue=='SPD', 'a'] <- bestPars[['a.SPD']]
      fit[idx1&fit$cue=='ACC', 'a'] <- bestPars[['a.ACC']]
      #      fit[idx1&fit$cue=='NEU', 'a'] <- bestPars[['a.NEU']]
    } else if(modelN == 2) {
      fit[idx1, 'a'] <- bestPars[['a']]
      fit[idx1&fit$cue=='SPD', 'eta1'] <- bestPars[['eta1.SPD']]
      fit[idx1&fit$cue=='ACC', 'eta1'] <- bestPars[['eta1.ACC']]
      #      fit[idx1&fit$cue=='NEU', 'eta1'] <- bestPars[['eta1.NEU']]
    } else if(modelN == 3) {
      fit[idx1, 'a'] <- bestPars[['a']]
      fit[idx1, 'eta1'] <- bestPars[['eta1']]
    } else if(modelN == 4) {
      fit[idx1&fit$cue=='SPD', 'a'] <- bestPars[['a.SPD']]
      fit[idx1&fit$cue=='ACC', 'a'] <- bestPars[['a.ACC']]
      #      fit[idx1&fit$cue=='NEU', 'a'] <- bestPars[['a.NEU']]
      fit[idx1&fit$cue=='SPD', 'eta1'] <- bestPars[['eta1.SPD']]
      fit[idx1&fit$cue=='ACC', 'eta1'] <- bestPars[['eta1.ACC']]
      #      fit[idx1&fit$cue=='NEU', 'eta1'] <- bestPars[['eta1.NEU']]
    }
    fit[idx1, 'BIC'] <- 2*fit[idx1, 'negLL'] + log(fit[idx1, 'nObs'])*fit[idx1, 'nPars']
    fit[idx1, 'AIC'] <- 2*fit[idx1, 'negLL'] + 2*fit[idx1, 'nPars']
    # fit[idx1, 'minBIC'] <- fit[idx1, 'BIC']-min(fit[idx1, 'BIC'])
  }
}





###
### only if fit on single block
#fit[fit$modelN %in% c(1,2), 'nPars'] <- fit[fit$modelN %in% c(1,2), 'nPars']-1
#fit[fit$modelN==4, 'nPars'] <- fit[fit$modelN==4, 'nPars']-2

fitWide = reshape(data=fit,#[fit$modelN!=3,],
                  v.names=c('eta1', 'eta2', 'a', 't0', 'm', 'mrt', 'accObjective', 'accSubjective'), 
                  direction='wide', idvar = c('pp', 'modelN'), timevar=c('cue'))
for(pp in fitWide$pp) {
  idx <- fitWide$pp==pp
  fitWide[idx, 'minBIC'] = fitWide[idx, 'BIC'] - min(fitWide[idx, 'BIC'])
  fitWide[idx, 'minAIC'] = fitWide[idx, 'AIC'] - min(fitWide[idx, 'AIC'])
}

fitWide[fitWide$minBIC==0, 'modelN']
fitWide[fitWide$minAIC==0, 'modelN']

winningModels <- data.frame(pp=unique(fitWide$pp), BIC=NA, AIC=NA, nObs=NA)
for(pp in fitWide$pp) {
  winningModels[winningModels$pp==pp, 'BIC'] = fitWide[fitWide$pp==pp & fitWide$minBIC==0, 'modelN']
  winningModels[winningModels$pp==pp, 'AIC'] = fitWide[fitWide$pp==pp & fitWide$minAIC==0, 'modelN']
  winningModels[winningModels$pp==pp, 'nObs'] = fitWide[fitWide$pp==pp & fitWide$minBIC==0, 'nObs']
}
table(winningModels$BIC)
table(winningModels$BIC)/sum(table(winningModels$BIC))


#table(winningModels$AIC)



####
fitWide$accObjDiff <- fitWide$accObjective.ACC-fitWide$accObjective.SPD
fitWide$accSubjDiff <- fitWide$accSubjective.ACC-fitWide$accSubjective.SPD
fitWide$mrtDiff <- fitWide$mrt.ACC-fitWide$mrt.SPD

mean(fitWide[fitWide$modelN==1, 'accObjDiff'])
mean(fitWide[fitWide$modelN==1, 'accSubjDiff'])