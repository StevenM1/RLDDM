rm(list=ls())  # fresh start

# Libraries, directories ---------------------------------------------------------------
source('./setup_locale.R') # gets dataDir, workDir
library(rtdists)
library(DEoptim)
library(RLDDM)
library(moments)
source(file.path(workDir, 'analysis/src/fittingFunctions.R'))
source(file.path(workDir, 'analysis/src/racingWaldDists.R'))
source(file.path(workDir, 'analysis/src/plottingFunctions.R'))

# What to plot? -----------------------------------------------------------
modelN <- 'LR-Q_CF-RW-AccAdvantages_SAT-NULL' # 'LR-Q_CF-DDM_SAT-NULL'  #'LR-Q_CF-RW-AccAdvantages_SAT-NULL'
modelType <- sub(".*?CF-(.*?)_.*", "\\1", modelN)
exp <- 'exp2'
block <- 'Calibration'
savePlots <- TRUE
resDir <- file.path(workDir, 'fits', 'Barbara', exp, paste0('model-', modelType))
nSim = 50

# load data
load(file.path(dataDir, paste0('data_', exp, '.Rdata')))
allPps <- unique(dat$pp)

# Plot RTs and learning -----------------------------------------------------------
if(savePlots) pdf(file=file.path(resDir, paste0('model-', modelN, '_', block, '_overall_cdf_bysub.pdf')))
def.par <- par(no.readonly = TRUE) # save default, for resetting layout
all.sims.qps <- all.dat.qps <- list()
for(pp in allPps) {
  # Plot Model --------------------------------------------------------------
  fitFn = file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_', block, '.rdat'))
  if(!file.exists(fitFn)) {print(paste0('sub ', pp, ' doesnt exist...')); next}
  load(fitFn)

  # Layout
  nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
  #  layout.show(nf)
  par(mar=c(5, 4, 2, 2) + 0.1)
  
  # prepare data
  d = prepareForFitting(dat[dat$pp==pp&as.character(dat$Block)==block,])
  df=d$df
  outcomes=d$outcomes
  VVchoiceIdx=d$VVchoiceIdx
  values=d$values
  choice=d$choice

  # Set-up
  source(file.path(workDir, 'analysis', 'models', paste0(modelN, '.R')))
  
  # plot, get quantiles
  out = plotSingleSub()
  all.dat.qps[[pp]] = out$all.dat.qps
  all.sims.qps[[pp]] = out$all.sims.qps
}

# overall
all.sims.qps.means <- getMeanSimQps(all.sims.qps) 
all.dat.qps.means <- getMeanDatQps(all.dat.qps)

# plot
nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
par(mar=c(5, 4, 2, 2) + 0.1)
plot(0, 0, type='n', xlim=c(0, 1.2), ylim=c(0, 1), xlab='RT', ylab='Def. prob.')
abline(v=seq(0, 10, .25), lty='dotted', col='grey')
abline(h=seq(0, 1, .1), lty='dotted', col='grey')
plotCDF(qps=all.dat.qps.means$r1, add=TRUE, type='b', lty=3, lwd=3, col=1)
plotCDF(qps=all.dat.qps.means$r2, add=TRUE, type='b', lty=3, lwd=3, col=2)
plotCDF(all.sims.qps.means$r1, add=TRUE, type='b', pch=4, lwd=2, col=1)
plotCDF(all.sims.qps.means$r2, add=TRUE, type='b', pch=4, lwd=2, col=2)
legend('bottomright', c('Data', 'Model'), lty=c(3, 1), lwd=c(3,2), pch=c(1, 4), bty='n')
title(main=paste('Across subjects'))

# rt by bin
plot(x=1:10, y=all.sims.qps.means$rtByBin, type='l', ylim=c(0.4, .8), xlab='Trial (binned)', ylab='Mean RT')
lines(x=1:10, y=all.dat.qps.means$meanRtByBin, type='b')
arrows(x0=1:10, x1=1:10, 
       y0=all.dat.qps.means$meanRtByBin-(all.dat.qps.means$sdRtByBin)/sqrt(length(allPps)), 
       y1=all.dat.qps.means$meanRtByBin+(all.dat.qps.means$sdRtByBin)/sqrt(length(allPps)), 
       angle = 90, code = 3, length = .05)
legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')

# acc by bin
plot(x=1:10, y=all.sims.qps.means$accByBin, type='l', ylim=c(0.3, 1), xlab='Trial (binned)', ylab='Accuracy')
lines(x=1:10, y=all.dat.qps.means$meanAccByBin, type='b')
arrows(x0=1:10, x1=1:10, 
       y0=all.dat.qps.means$meanAccByBin-(all.dat.qps.means$sdAccByBin)/sqrt(length(allPps)), 
       y1=all.dat.qps.means$meanAccByBin+(all.dat.qps.means$sdAccByBin)/sqrt(length(allPps)), 
       angle = 90, code = 3, length = .05)
legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')

#skew by bin
plot(x=1:10, y=all.sims.qps.means$skewByBin, type='l', ylim=c(0, 3), xlab='Trial (binned)', ylab='Skewness')
lines(x=1:10, y=all.dat.qps.means$meanSkewByBin, type='b')
arrows(x0=1:10, x1=1:10, 
       y0=all.dat.qps.means$meanSkewByBin-(all.dat.qps.means$sdSkewByBin)/sqrt(length(allPps)), 
       y1=all.dat.qps.means$meanSkewByBin+(all.dat.qps.means$sdSkewByBin)/sqrt(length(allPps)), 
       angle = 90, code = 3, length = .05)
legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')

#IQR by bin
plot(x=1:10, y=all.sims.qps.means$iqrByBin, type='l', ylim=c(0, 1), xlab='Trial (binned)', ylab='IQR')
lines(x=1:10, y=all.dat.qps.means$meanIqrByBin, type='b')
arrows(x0=1:10, x1=1:10, 
       y0=all.dat.qps.means$meanIqrByBin-(all.dat.qps.means$sdIqrByBin)/sqrt(length(allPps)), 
       y1=all.dat.qps.means$meanIqrByBin+(all.dat.qps.means$sdIqrByBin)/sqrt(length(allPps)), 
       angle = 90, code = 3, length = .05)
legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')

if(savePlots) dev.off()


# End overall plot --------------------------------------------------------






# Plot per condition ------------------------------------------------------
load(file.path(dataDir, paste0('data_', exp, '.Rdata')))
allPps <- unique(dat$pp)

# Plot RTs and learning -----------------------------------------------------------
def.par <- par(no.readonly = TRUE) # save default, for resetting layout
all.sims.qps <- all.dat.qps <- list()
for(pp in allPps) {
  # Plot Model --------------------------------------------------------------
  load(file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_', block, '.rdat')))
  if(savePlots) pdf(file=file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_', block, '_bycue_cdf.pdf')))

  # Layout
  nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
  #  layout.show(nf)
  par(mar=c(5, 4, 2, 2) + 0.1)

  # prepare data
  d = prepareForFitting(dat[dat$pp==pp&as.character(dat$Block)==block,])
  df=d$df
  outcomes=d$outcomes
  VVchoiceIdx=d$VVchoiceIdx
  values=d$values
  choice=d$choice
  etaMax <- 1

  # Set-up
  modelSetup <- model(modelSpec)
  names(bestPars) <- modelSetup$p.vector
  fitModel <- obj(bestPars, 
                  rt=df$rt,
                  parNames=modelSetup$p.vector,
                  constants=modelSetup$constants,
                  choice=choice,
                  condition=df$cue,
                  values=values,
                  outcomes=outcomes, 
                  VVchoiceIdx=VVchoiceIdx,
                  returnType='full',
                  learningRule = modelSetup$learningRule)
  
  for(i in 1:length(modelSetup$constants)) bestPars[[names(modelSetup$constants[i])]] <- modelSetup$constants[[i]]
  choicePars = transformChoicePars(condition = df$cue, pars=bestPars, fitModel$VV)
  
  sims.qps <- lapply(1:nSim, function(x) {
    cat('.', sep='')
    simdat <- rWaldRaceSM(n=nrow(df), 
                          A=choicePars[['A']], 
                          B=choicePars[['B']], 
                          v=list(choicePars[['v1']], choicePars[['v2']]), 
                          t0=choicePars[['t0']])
    simdat$rt <- simdat$RT
    simdat$response <- ifelse(simdat$R == 1, 'upper', 'lower')
    simdatAll <- simdat
    out <- list()
    for(cond in c('SPD', 'ACC')) {
      simdat <- simdatAll[df$cue==cond,]
      qps1 <- quantile(simdat[simdat$response=='upper','rt'], quants)
      qps2 <- quantile(simdat[simdat$response=='lower','rt'], quants)
      simdat$trialN <- 1:nrow(simdat)
      simdat$bin <- cut(simdat$trialN, 10, labels=FALSE)
      rtByBin <- aggregate(rt~bin, simdat, mean)
      accByBin <- aggregate((response=='upper')~bin, simdat, mean)
      skewByBin <- aggregate(rt~bin, simdat, skewness)
      iqrByBin <- aggregate(rt~bin, simdat, IQR)

      out[[cond]] = list('r2'=list(prop=mean(simdat$response=='upper'),
                                   qps=qps1,
                                   def.probability=quants*mean(simdat$response=='upper')),
                         'r1'=list(prop=mean(simdat$response=='lower'),
                                   qps=qps2,
                                   def.probability=quants*mean(simdat$response=='lower')),
                         'rtByBin'=rtByBin[,2],
                         'accByBin'=accByBin[,2],
                         'skewByBin'=skewByBin[,2],
                         'iqrByBin'=iqrByBin[,2]
      )
    }
    out
  })

  dat.qps <- list()
  for(cue in c('SPD', 'ACC')) {
    dfThisCue <- df[df$cue==cue,]
    idx <- df$cue == cue
    qps1 <- quantile(dfThisCue$rt[choice[idx]==2], quants)
    qps2 <- quantile(dfThisCue$rt[choice[idx]==1], quants)
    dfThisCue$trialN <- 1:nrow(dfThisCue)
    dfThisCue$bin <- cut(dfThisCue$trialN, 10, labels=FALSE)
    rtByBin <- aggregate(rt~bin, dfThisCue, mean)
    accByBin <- aggregate(choiceIsHighP~bin, dfThisCue, mean)
    skewByBin <- aggregate(rt~bin, dfThisCue, skewness)
    iqrByBin <- aggregate(rt~bin, dfThisCue, IQR)

    dat.qps[[cue]] <- list('r2'=list(prop=mean(choice[idx]==2),
                                     qps=qps1,
                                     def.probability=quants*mean(choice[idx]==2)),
                           'r1'=list(prop=mean(choice[idx]==1),
                                     qps=qps2,
                                     def.probability=quants*mean(choice[idx]==1)),
                           'rtByBin' = rtByBin[,2],
                           'accByBin' = accByBin[,2],
                           'skewByBin' = skewByBin[,2],
                           'iqrByBin' = iqrByBin[,2]
    )
  }
  all.dat.qps[[pp]] = dat.qps
  all.sims.qps[[pp]] = list()

  # start plot
  xlim.max = 1.5
  ylim.max = 1
  plot(0, 0, type='n', xlim=c(0, xlim.max), ylim=c(0, ylim.max), xlab='RT', ylab='Def. prob.')
  # expand, loop
  for(cond in c('SPD', 'ACC')) {
    sim.qps.thisCond = lapply(sims.qps, function(x) x[[cond]])
    dat.qps.thisCond = dat.qps[[cond]] #lapply(dat.qps, function(x) x[[cond]])

    # to fix
    # xlim.max <- max(max(sapply(sims.qps, function(x) sapply(c('r1', 'r2'), function(y) max(x[[y]]$qps)))),
    #                 max(sapply(dat.qps, function(x) max(x$qps))))
    #  ylim.max <- max(max(sapply(sims.qps, function(x) sapply(c('r1', 'r2'), function(y) max(x[[y]]$prop)))),
    #                  max(sapply(dat.qps, function(x) max(x$prop))))
    #    ylim.max = 1
    #    if(is.na(xlim.max)) {
    #      warning('xlim.max could not be determined')
    #    xlim.max = 1.5
    #    }
    abline(v=seq(0, 10, .25), lty='dotted', col='grey')
    abline(h=seq(0, 1, .1), lty='dotted', col='grey')
    abline(v=choicePars[['t0']][1], lty=2, col='black')
    abline(v=choicePars[['t0']][1]+choicePars[['st0']][1], lty=2, col='black')

    # plot individual sims
    tmp = lapply(sim.qps.thisCond, function(x) lapply(1:2, function(y) plotCDF(qps=x[[y]], add=TRUE, type='l', lty=1, col='grey')))
    sims.mean <- list('r1'=list(qps=apply(do.call(rbind, lapply(sim.qps.thisCond, function(x) x$r1$qps)), 2, mean),
                                def.prob=apply(do.call(rbind, lapply(sim.qps.thisCond, function(x) x$r1$def.probability)), 2, mean),
                                prop=mean(sapply(sim.qps.thisCond, function(x) x$r1$prop))),
                      'r2'=list(qps=apply(do.call(rbind, lapply(sim.qps.thisCond, function(x) x$r2$qps)), 2, mean),
                                def.prob=apply(do.call(rbind, lapply(sim.qps.thisCond, function(x) x$r2$def.probability)), 2, mean),
                                prop=mean(sapply(sim.qps.thisCond, function(x) x$r2$prop))),
                      'rtByBin'=apply(sapply(sim.qps.thisCond, function(x) x$rtByBin), 1, mean),
                      'accByBin'=apply(sapply(sim.qps.thisCond, function(x) x$accByBin), 1, mean),
                      'skewByBin'=apply(sapply(sim.qps.thisCond, function(x) x$skewByBin), 1, mean),
                      'iqrByBin'=apply(sapply(sim.qps.thisCond, function(x) x$iqrByBin), 1, mean)
    )
    all.sims.qps[[pp]][[cond]] = sims.mean
    plotCDF(sims.mean$r1, add=TRUE, type='b', pch=4, lwd=2, col=1)
    plotCDF(sims.mean$r2, add=TRUE, type='b', pch=4, lwd=2, col=2)
    plotCDF(qps=dat.qps.thisCond$r1, add=TRUE, type='b', lty=3, lwd=3, col=1)
    plotCDF(qps=dat.qps.thisCond$r2, add=TRUE, type='b', lty=3, lwd=3, col=2)
  }
  legend('bottomright', c('Data', 'Model'), lty=c(3, 1), lwd=c(3,2), pch=c(1, 4), bty='n')
  title(main=paste('Sub', pp),
        sub=paste(names(bestPars), '=', c(round(bestPars[1],4), round(bestPars[2:length(bestPars)],2)), collapse='  '))

  plot(x=1:10, y=sims.mean$rtByBin, type='l', ylim=c(0.3, .9), xlab='Trial (binned)', ylab='Mean RT')
  lines(x=1:10, y=dat.qps$rtByBin, type='b')
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')

  plot(x=1:10, y=sims.mean$accByBin, type='l', ylim=c(0.3, 1), xlab='Trial (binned)', ylab='Accuracy')
  lines(x=1:10, y=dat.qps$accByBin, type='b')

  plot(x=1:10, y=sims.mean$skewByBin, type='l', ylim=c(0, 2.5), xlab='Trial (binned)', ylab='Skewness')
  lines(x=1:10, y=dat.qps$skewByBin, type='b')

  plot(x=1:10, y=sims.mean$iqrByBin, type='l', ylim=c(0.3, 1), xlab='Trial (binned)', ylab='IQR')
  lines(x=1:10, y=dat.qps$iqrByBin, type='b')
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')
  if(savePlots) dev.off()
}

# overall
if(savePlots) pdf(file=file.path(resDir, paste0('sub-ALL_model-', modelN, '_', block, '_bycue_cdf.pdf')))
nf <- layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
for(cond in c('SPD', 'ACC')) {
  par(mar=c(5, 4, 2, 2) + 0.1)
  plot(0, 0, type='n', xlim=c(0, xlim.max), ylim=c(0, ylim.max), xlab='RT', ylab='Def. prob.')
  all.sims.qps.thisCond <- lapply(all.sims.qps, function(x) x[[cond]])
  all.dat.qps.thisCond <-lapply(all.dat.qps, function(x) x[[cond]])
  # plot mean of data
  all.sims.qps.means = list('r1' = list(qps=apply(t(sapply(all.sims.qps.thisCond, function(x) x$r1$qps)), 2, mean, na.rm=TRUE),
                                        def.prob=apply(t(sapply(all.sims.qps.thisCond, function(x) x$r1$def.prob)), 2, mean, na.rm=TRUE)),
                            'r2' = list(qps=apply(t(sapply(all.sims.qps.thisCond, function(x) x$r2$qps)), 2, mean, na.rm=TRUE),
                                        def.prob=apply(t(sapply(all.sims.qps.thisCond, function(x) x$r2$def.prob)), 2, mean, na.rm=TRUE)),
                            'rtByBin' = apply(sapply(all.sims.qps.thisCond, function(x) x$rtByBin), 1, mean, na.rm=TRUE),
                            'accByBin' = apply(sapply(all.sims.qps.thisCond, function(x) x$accByBin), 1, mean, na.rm=TRUE),
                            'skewByBin' = apply(sapply(all.sims.qps.thisCond, function(x) x$skewByBin), 1, mean, na.rm=TRUE),
                            'iqrByBin' = apply(sapply(all.sims.qps.thisCond, function(x) x$iqrByBin), 1, mean, na.rm=TRUE)
  )
  all.dat.qps.means = list('r1' = list(qps=apply(t(sapply(all.dat.qps.thisCond, function(x) x$r1$qps)), 2, mean, na.rm=TRUE),
                                       def.prob=apply(t(sapply(all.dat.qps.thisCond, function(x) x$r1$def.prob)), 2, mean, na.rm=TRUE)),
                           'r2' = list(qps=apply(t(sapply(all.dat.qps.thisCond, function(x) x$r2$qps)), 2, mean, na.rm=TRUE),
                                       def.prob=apply(t(sapply(all.dat.qps.thisCond, function(x) x$r2$def.prob)), 2, mean, na.rm=TRUE)),
                           'meanRtByBin' = apply(sapply(all.dat.qps.thisCond, function(x) x$rtByBin), 1, mean, na.rm=TRUE),
                           'sdRtByBin' = apply(sapply(all.dat.qps.thisCond, function(x) x$rtByBin), 1, sd, na.rm=TRUE),
                           'meanAccByBin' = apply(sapply(all.dat.qps.thisCond, function(x) x$accByBin), 1, mean, na.rm=TRUE),
                           'sdAccByBin' = apply(sapply(all.dat.qps.thisCond, function(x) x$accByBin), 1, sd, na.rm=TRUE),
                           'meanSkewByBin' = apply(sapply(all.dat.qps.thisCond, function(x) x$skewByBin), 1, mean, na.rm=TRUE),
                           'sdSkewByBin' = apply(sapply(all.dat.qps.thisCond, function(x) x$skewByBin), 1, sd, na.rm=TRUE),
                           'meanIqrByBin' = apply(sapply(all.dat.qps.thisCond, function(x) x$iqrByBin), 1, mean, na.rm=TRUE),
                           'sdIqrByBin' = apply(sapply(all.dat.qps.thisCond, function(x) x$iqrByBin), 1, sd, na.rm=TRUE)
  )
  abline(v=seq(0, 10, .25), lty='dotted', col='grey')
  abline(h=seq(0, 1, .1), lty='dotted', col='grey')
  plotCDF(qps=all.dat.qps.means$r1, add=TRUE, type='b', lty=3, lwd=3, col=1)
  plotCDF(qps=all.dat.qps.means$r2, add=TRUE, type='b', lty=3, lwd=3, col=2)
  plotCDF(all.sims.qps.means$r1, add=TRUE, type='b', pch=4, lwd=2, col=1)
  plotCDF(all.sims.qps.means$r2, add=TRUE, type='b', pch=4, lwd=2, col=2)
  legend('bottomright', c('Data', 'Model'), lty=c(3, 1), lwd=c(3,2), pch=c(1, 4), bty='n')
  title(main=paste('Across subjects'))

  # rt by bin
  plot(x=1:10, y=all.sims.qps.means$rtByBin, type='l', ylim=c(0.4, .7), xlab='Trial (binned)', ylab='Mean RT')
  lines(x=1:10, y=all.dat.qps.means$meanRtByBin, type='b')
  arrows(x0=1:10, x1=1:10,
         y0=all.dat.qps.means$meanRtByBin-(all.dat.qps.means$sdRtByBin)/sqrt(length(allPps)),
         y1=all.dat.qps.means$meanRtByBin+(all.dat.qps.means$sdRtByBin)/sqrt(length(allPps)),
         angle = 90, code = 3, length = .05)
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')

  # accuracy per bin
  plot(x=1:10, y=all.sims.qps.means$accByBin, type='l', ylim=c(0.3, 1), xlab='Trial (binned)', ylab='Accuracy')
  lines(x=1:10, y=all.dat.qps.means$meanAccByBin, type='b')
  arrows(x0=1:10, x1=1:10,
         y0=all.dat.qps.means$meanAccByBin-(all.dat.qps.means$sdAccByBin)/sqrt(length(allPps)),
         y1=all.dat.qps.means$meanAccByBin+(all.dat.qps.means$sdAccByBin)/sqrt(length(allPps)),
         angle = 90, code = 3, length = .05)
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')

  # skewness per bin
  plot(x=1:10, y=all.sims.qps.means$skewByBin, type='l', ylim=c(0.3, 1), xlab='Trial (binned)', ylab='Skewness')
  lines(x=1:10, y=all.dat.qps.means$meanSkewByBin, type='b')
  arrows(x0=1:10, x1=1:10,
         y0=all.dat.qps.means$meanSkewByBin-(all.dat.qps.means$sdSkewByBin)/sqrt(length(allPps)),
         y1=all.dat.qps.means$meanSkewByBin+(all.dat.qps.means$sdSkewByBin)/sqrt(length(allPps)),
         angle = 90, code = 3, length = .05)
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')

  # iqr per bin
  plot(x=1:10, y=all.sims.qps.means$iqrByBin, type='l', ylim=c(0.3, 1), xlab='Trial (binned)', ylab='IQR')
  lines(x=1:10, y=all.dat.qps.means$meanIqrByBin, type='b')
  arrows(x0=1:10, x1=1:10,
         y0=all.dat.qps.means$meanIqrByBin-(all.dat.qps.means$sdIqrByBin)/sqrt(length(allPps)),
         y1=all.dat.qps.means$meanIqrByBin+(all.dat.qps.means$sdIqrByBin)/sqrt(length(allPps)),
         angle = 90, code = 3, length = .05)
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')
  title(cond, outer=TRUE)
}
# if(savePlots) dev.off()
