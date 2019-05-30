rm(list=ls())  # fresh start

# Libraries ---------------------------------------------------------------
.libPaths(c(.libPaths(), '/home/stevenm/rpackages')) # for tux server
library(rtdists)
library(DEoptim)
library(RLDDM)

# Directories -------------------------------------------------------------
dataDir <- '~/surfdrive/Testfolder/exp2'
prepDataDir <- '~/surfdrive/data/learningTask/Barbara_preprocessed'
workDir <- '/Users/steven/Sync/PhDprojects/RLDDM'
source(file.path(workDir, 'analysis/fittingFunctions.R'))

# What to plot? -----------------------------------------------------------
exp <- 'exp2'
resDir <- file.path(workDir, 'fits', 'Barbara', exp, paste0('model-', modelType))

# single participant
modelN <- 1
load(file.path(prepDataDir, paste0('data_', exp, '.Rdata')))

# Plot data ---------------------------------------------------------------
# allPps <- unique(dat$pp)
# pp <- allPps[2]
# df = dat[dat$pp==pp & dat$Block == 'Miniblocks',]
# df$choiceIsHighP[is.na(df$rt)] = FALSE
# # Data
# probs <- seq(0.1, .9, .2)
# par(mfrow=c(1,1), lwd=2)
# plot(0, 0, type='n', xlim=c(min(df$rt, na.rm=TRUE), max(df$rt, na.rm = TRUE)), ylim=c(0, 1), ylab='defective CDF', xlab='RT')
# for(cue in c('SPD', 'ACC')) {
#   choiceProb = mean(df$choiceIsHighP[df$Cue==cue])
#   qps1 = quantile(df$rt[df$Cue==cue & df$choiceIsHighP==TRUE], probs, na.rm=TRUE)
#   qps2 = quantile(df$rt[df$Cue==cue & df$choiceIsHighP==FALSE], probs, na.rm=TRUE)
#   lines(qps1, probs*(choiceProb), type='b', lty=ifelse(cue=='SPD', 2, 1), col=ifelse(cue=='SPD', 2, 1))
#   lines(qps2, probs*(1-choiceProb), type='b', lty=ifelse(cue=='SPD', 2, 1), col=ifelse(cue=='SPD', 2, 1))
# }
# legend('topright', c('SPD', 'ACC'), lty=c(2, 1), lwd=c(2,2), col=c(2, 1), bty='n', cex=1.2)
# 

# Plot RTs and learning -----------------------------------------------------------
def.par <- par(no.readonly = TRUE) # save default, for resetting layout
all.sims.qps <- all.dat.qps <- list()
for(pp in unique(dat$pp)) {
  # Plot Model --------------------------------------------------------------
  load(file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_Miniblocks.rdat')))
  if(savePlots) pdf(file=file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_Miniblocks_cdf.pdf')))
  
  # Layout
  nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), respect = TRUE)
  layout.show(nf)
  par(mar=c(5, 4, 2, 2) + 0.1)

  # prepare data
  d = prepareForFitting(dat[dat$pp==pp&as.character(dat$Block)=='Miniblocks',]) 
  df=d$df
  outcomes=d$outcomes
  VVchoiceIdx=d$VVchoiceIdx
  values=d$values
  choice=d$choice
  etaMax <- 1
  
  # Set-up
  modelSetup <- setupModel(modelN=modelN)
  names(bestPars) <- modelSetup$parNames
  model <- objRLDDMMultiCond(bestPars, 
                             rt=df$rt,
                             parNames=modelSetup$parNames,
                             fixedPars=modelSetup$fixedPars,
                             choice=choice,
                             condition=df$cue,
                             values=values, 
                             outcomes=outcomes, 
                             VVchoiceIdx=VVchoiceIdx,
                             modelN=modelN,
                             returnType='full')
  at0m = getPars(condition = df$cue, modelN = modelN, pars=bestPars)
  pars.by.trial <- list(v=at0m[['m']]*(model$VV[,1]-model$VV[,2]),
                        a=at0m[['a']], t0=at0m[['t0']], z=at0m[['a']]*.5, s=1)
  cat('simulating', sep='')
  quants <- c(.1, .3, .5, .7, .9)
  
  library(snowfall)
  sfInit(parallel=TRUE, cpus=4)
  sfLibrary(rtdists)
  sfExportAll()
  sims.qps <- sfLapply(1:100, function(x) {
    cat('.', sep='')
    simdat <- rdiffusion(n=nrow(df), 
                         a=pars.by.trial[['a']], 
                         v=pars.by.trial[['v']], 
                         z=pars.by.trial[['z']], s=pars.by.trial[['s']], t0=pars.by.trial[['t0']])
    qps1 <- quantile(simdat[simdat$response=='upper','rt'], quants)
    qps2 <- quantile(simdat[simdat$response=='lower','rt'], quants)
    simdat$trialN <- 1:nrow(simdat)
    simdat$bin <- cut(simdat$trialN, 10, labels=FALSE)
    rtByBin <- aggregate(rt~bin, simdat, mean)
    accByBin <- aggregate((response=='upper')~bin, simdat, mean)

    list('r2'=list(prop=mean(simdat$response=='upper'),
                   qps=qps1,
                   def.probability=quants*mean(simdat$response=='upper')),
         'r1'=list(prop=mean(simdat$response=='lower'),
                   qps=qps2,
                   def.probability=quants*mean(simdat$response=='lower')),
         'rtByBin'=rtByBin[,2], 
         'accByBin'=accByBin[,2])
  })
  sfStop()
  
  qps1 <- quantile(df$rt[choice==2], quants)
  qps2 <- quantile(df$rt[choice==1], quants)
  df$trialN <- 1:nrow(df)
  df$bin <- cut(df$trialN, 10, labels=FALSE)
  rtByBin <- aggregate(rt~bin, df, mean)
  accByBin <- aggregate(choiceIsHighP~bin, df, mean)
  
  dat.qps <- list('r2'=list(prop=mean(choice==2),
                            qps=qps1, 
                            def.probability=quants*mean(choice==2)),
                  'r1'=list(prop=mean(choice==1),
                            qps=qps2,
                            def.probability=quants*mean(choice==1)))
  xlim.max <- max(max(sapply(sims.qps, function(x) sapply(c('r1', 'r2'), function(y) max(x[[y]]$qps)))),
                  max(sapply(dat.qps, function(x) max(x$qps))))
  ylim.max <- max(max(sapply(sims.qps, function(x) sapply(c('r1', 'r2'), function(y) max(x[[y]]$prop)))),
                  max(sapply(dat.qps, function(x) max(x$prop))))
  dat.qps$rtByBin = rtByBin[,2]
  dat.qps$accByBin = accByBin[,2]
  all.dat.qps[[pp]] = dat.qps
  
  plot(0, 0, type='n', xlim=c(0, xlim.max), ylim=c(0, ylim.max), xlab='RT', ylab='Def. prob.')
  abline(v=seq(0, 10, .25), lty='dotted', col='grey')
  abline(h=seq(0, 1, .1), lty='dotted', col='grey')
  tmp = lapply(sims.qps, function(x) lapply(1:2, function(y) plotCDF(qps=x[[y]], add=TRUE, type='l', lty=1, col='grey')))
  sims.mean <- list('r1'=list(qps=apply(do.call(rbind, lapply(sims.qps, function(x) x$r1$qps)), 2, mean),
                              def.prob=apply(do.call(rbind, lapply(sims.qps, function(x) x$r1$def.probability)), 2, mean),
                              prop=mean(sapply(sims.qps, function(x) x$r1$prop))),
                    'r2'=list(qps=apply(do.call(rbind, lapply(sims.qps, function(x) x$r2$qps)), 2, mean),
                              def.prob=apply(do.call(rbind, lapply(sims.qps, function(x) x$r2$def.probability)), 2, mean),
                              prop=mean(sapply(sims.qps, function(x) x$r2$prop))),
                    'rtByBin'=apply(sapply(sims.qps, function(x) x$rtByBin), 1, mean),
                    'accByBin'=apply(sapply(sims.qps, function(x) x$accByBin), 1, mean))
  all.sims.qps[[pp]] = sims.mean
  plotCDF(qps=dat.qps$r1, add=TRUE, type='b', lty=3, lwd=3, col=1)
  plotCDF(qps=dat.qps$r2, add=TRUE, type='b', lty=3, lwd=3, col=2)
  plotCDF(sims.mean$r1, add=TRUE, type='b', pch=4, lwd=2, col=1)
  plotCDF(sims.mean$r2, add=TRUE, type='b', pch=4, lwd=2, col=2)
  legend('bottomright', c('Data', 'Model'), lty=c(3, 1), lwd=c(3,2), pch=c(1, 4), bty='n')
  title(main=paste('Sub', pp), sub=paste(names(bestPars), '=', c(round(bestPars[1],4), round(bestPars[2:length(bestPars)],2)), collapse='  '))
  
  plot(x=1:10, y=sims.mean$rtByBin, type='l', ylim=c(0.3, .9), xlab='Trial (binned)', ylab='Mean RT')
  lines(x=1:10, y=dat.qps$rtByBin, type='b')
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')
  
  plot(x=1:10, y=sims.mean$accByBin, type='l', ylim=c(0.3, 1), xlab='Trial (binned)', ylab='Accuracy')
  lines(x=1:10, y=dat.qps$accByBin, type='b')
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')
  if(savePlots) dev.off()
}

# overall


# plot mean of data
all.sims.qps.means = list('r1' = list(qps=apply(t(sapply(all.sims.qps, function(x) x$r1$qps)), 2, mean),
                                      def.prob=apply(t(sapply(all.sims.qps, function(x) x$r1$def.prob)), 2, mean)),
                          'r2' = list(qps=apply(t(sapply(all.sims.qps, function(x) x$r2$qps)), 2, mean),
                                      def.prob=apply(t(sapply(all.sims.qps, function(x) x$r2$def.prob)), 2, mean)),
                          'rtByBin' = apply(sapply(all.sims.qps, function(x) x$rtByBin), 1, mean),
                          'accByBin' = apply(sapply(all.sims.qps, function(x) x$accByBin), 1, mean))
all.dat.qps.means = list('r1' = list(qps=apply(t(sapply(all.dat.qps, function(x) x$r1$qps)), 2, mean),
                                     def.prob=apply(t(sapply(all.dat.qps, function(x) x$r1$def.prob)), 2, mean)),
                         'r2' = list(qps=apply(t(sapply(all.dat.qps, function(x) x$r2$qps)), 2, mean),
                                     def.prob=apply(t(sapply(all.dat.qps, function(x) x$r2$def.prob)), 2, mean)),
                         'meanRtByBin' = apply(sapply(all.dat.qps, function(x) x$rtByBin), 1, mean),
                         'sdRtByBin' = apply(sapply(all.dat.qps, function(x) x$rtByBin), 1, sd),
                         'meanAccByBin' = apply(sapply(all.dat.qps, function(x) x$accByBin), 1, mean),
                         'sdAccByBin' = apply(sapply(all.dat.qps, function(x) x$accByBin), 1, sd))
if(savePlots) pdf(file=file.path(resDir, paste0('sub-ALL_model-', modelN, '_Miniblocks_cdf.pdf')))
plot(0, 0, type='n', xlim=c(0, xlim.max), ylim=c(0, ylim.max), xlab='RT', ylab='Def. prob.')
abline(v=seq(0, 10, .25), lty='dotted', col='grey')
abline(h=seq(0, 1, .1), lty='dotted', col='grey')
plotCDF(qps=dat.qps$r1, add=TRUE, type='b', lty=3, lwd=3, col=1)
plotCDF(qps=dat.qps$r2, add=TRUE, type='b', lty=3, lwd=3, col=2)
plotCDF(sims.mean$r1, add=TRUE, type='b', pch=4, lwd=2, col=1)
plotCDF(sims.mean$r2, add=TRUE, type='b', pch=4, lwd=2, col=2)
legend('bottomright', c('Data', 'Model'), lty=c(3, 1), lwd=c(3,2), pch=c(1, 4), bty='n')
title(main=paste('Across subjects'))

# rt by bin
plot(x=1:10, y=all.sims.qps.means$rtByBin, type='l', ylim=c(0.5, .8), xlab='Trial (binned)', ylab='Mean RT')
lines(x=1:10, y=all.dat.qps.means$meanRtByBin, type='b')
arrows(x0=1:10, x1=1:10, 
       y0=all.dat.qps.means$meanRtByBin-(all.dat.qps.means$sdRtByBin)/sqrt(length(allPps)), 
       y1=all.dat.qps.means$meanRtByBin+(all.dat.qps.means$sdRtByBin)/sqrt(length(allPps)), 
       angle = 90, code = 3, length = .05)
legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')


plot(x=1:10, y=all.sims.qps.means$accByBin, type='l', ylim=c(0.3, 1), xlab='Trial (binned)', ylab='Accuracy')
lines(x=1:10, y=all.dat.qps.means$meanAccByBin, type='b')
arrows(x0=1:10, x1=1:10, 
       y0=all.dat.qps.means$meanAccByBin-(all.dat.qps.means$sdAccByBin)/sqrt(length(allPps)), 
       y1=all.dat.qps.means$meanAccByBin+(all.dat.qps.means$sdAccByBin)/sqrt(length(allPps)), 
       angle = 90, code = 3, length = .05)
legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')

if(savePlots) dev.off()






###### old stuff below
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
