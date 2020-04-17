# replicate behavioral data by sumulating with individual fitted parameters
expModelRep = function(modelName){
  set.seed(123)
  # create output directories
  dir.create("figures/expModelRep/")
  dir.create(sprintf("figures/expModelRep/%s",modelName))
  dir.create("genData/expModelRep")
  
  # load experiment parameters
  load('expParas.RData')
  nBlock = 2
  
  # load packages and sub functions 
  library("tidyverse")
  source("expSchematics.R")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  source("subFxs/analysisFxs.R") 
  source("subFxs/modelRep.R")
  
  # num of repetitions 
  nRep = 10
  
  # load empirical data 
  allData = loadAllData()
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$id
  nSub = length(ids)

  ## check fit
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  
  ################ compare AUCs and CIPs from empirical and replicated data ################
  ## empirical data 
  source("MFAnalysis.R")
  MFResults = MFAnalysis(isTrct = T)
  blockStats = MFResults[['blockStats']]
  blockStats = blockStats[blockStats$blockNum <= 2, ]
  muWTWEmp = blockStats$muWTW
  stdWTWEmp = blockStats$stdWTW
  
  ## rep
  repOutputs =  modelRep(trialData, ids, nRep, T, modelName)
  plotData = data.frame(mu =  repOutputs$muWTWRep_mu, std = repOutputs$stdWTWRep_mu,
                        empMu = muWTWEmp, empStd = stdWTWEmp,
                          passCheck = rep(passCheck, each = 2), 
                        condition = blockStats$condition) %>% filter(passCheck)
  save(repOutputs, file = sprintf("genData/expModelRep/%s_trct.RData", modelName))
  
  ## plot to compare AUCs
    plotData %>%
    ggplot(aes(empMu, mu)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21) + 
    geom_abline(slope = 1, intercept = 0)  + 
    ylab("Model-generated (s)") + xlab("Observed (s)") +
      ggtitle(sprintf("AUC, n = %d", sum(passCheck))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5))+
    scale_color_manual(values = conditionColors) +
      theme(legend.position = "none")
  fileName = sprintf("figures/expModelRep/%s/muWTW_muWTWRep.eps", modelName)
  fileName = sprintf("figures/expModelRep/%s/muWTW_muWTWRep.png", modelName) 
  ggsave(filename = fileName,  width = 4, height = 4)
  
  ## plot to compare CIPs
    plotData %>%
    ggplot(aes(empStd, std, shape = condition)) + 
    geom_point(size = 4, aes(color = condition), stroke = 1, shape = 21)  + 
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(bold(paste("Model-generated (s"^2,")")))) +
      xlab(expression(bold(paste("Observed (s"^2,")")))) + ggtitle(sprintf("CIP, n = %d", sum(passCheck))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + 
    scale_x_continuous(breaks = c(0, 5, 15), limits = c(0, 16)) + 
    scale_y_continuous(breaks = c(0, 5, 15), limits = c(0, 16)) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none")
  fileName = sprintf("figures/expModelRep/%s/stdWTW_stdWTWRep.eps", modelName)
  fileName = sprintf("figures/expModelRep/%s/stdWTW_stdWTWRep.png", modelName)
  ggsave(filename = fileName,  width = 4, height = 4)
  
  ################ compare WTW time courses from empirical and replicated data ################
  # rep again 
  repOutputs =  modelRep(trialData, ids, nRep, F, modelName)
  save(repOutputs, file = sprintf("genData/expModelRep/%s.RData", modelName))
  timeWTW_ = repOutputs$timeWTW_
  

  # plot to check time course of learning
 greenData = data.frame(
    xmin = c(0, blockSec),
    xmax = c(blockSec, blockSec * 2),
    cbal = c(1, 2)
  )
  purpleData = data.frame(
    xmin = c(blockSec, 0),
    xmax =  c(blockSec * 2, blockSec),
    cbal = c(1, 2)
  )
  plotData = data.frame(
    wtw = as.vector(timeWTW_),
    t = rep(seq(0, blockSec * 2-1, by = 1), nSub),
    cbal = rep(blockStats$cbal[blockStats$condition == "HP"], each = length(tGrid) * 2)
  ) %>% group_by(cbal, t) %>%
    dplyr::summarise(mu = mean(wtw, na.rm = F),
              se = sd(wtw, na.rm = F) / sqrt(length(wtw)),
              min = mu - se,
              max = mu + se) 
  plotData$mu[plotData$t == blockSec] = NA
  plotData %>% ggplot(aes(t, mu)) +
    geom_rect(data =  greenData, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 20),
              inherit.aes = F, fill = conditionColorBacks[1]) + 
    geom_rect(data =  purpleData, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 20),
              inherit.aes = F, fill = conditionColorBacks[2]) +
    geom_line() +
    facet_grid(cbal~.) +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.5) +
    facet_grid(cbal~.) + 
    scale_color_manual(values = conditionColors) +
    scale_fill_manual(values = conditionColorBacks) + myTheme +
    theme(legend.position = "none") +
    ylab("WTW (s)") + xlab("Task time (min)")  +
    scale_x_continuous(breaks = 0:2 * 10 * 60 , labels = 0:2 * 10 * 60 / 60) + 
    theme(legend.position = "none") + 
    scale_y_continuous(breaks = c(0, 10, 20), labels = c(0, 10, 20))
  fileName = sprintf("figures/expModelRep/%s/timecourse.eps", modelName)
  ggsave(filename = fileName,  width = 5, height = 3)
  fileName = sprintf("figures/expModelRep/%s/timecourse.png", modelName)
  ggsave(filename = fileName,  width = 5, height = 3)
}

