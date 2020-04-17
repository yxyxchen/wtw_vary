expModelFit = function(modelName, isFirstFit){
  # generate output directories
  dir.create("genData")
  dir.create("genData/expModelFit")
  dir.create("stanWarnings")
  
  # load experiment parameters
  load("expParas.RData")
  
  # load sub-functions and packages
  library("dplyr"); library("tidyr")
  source("subFxs/loadFxs.R")
  source("subFxs/helpFxs.R")
  source('subFxs/modelFitGroup.R')
  source("expSchematics.R")
  
  # prepare inputs
  allData = loadAllData()
  hdrData = allData$hdrData
  ids = unique(hdrData$id)
  trialData = allData$trialData
  outputDir = sprintf("genData/expModelFit/%s", modelName)
  config = list(
    nChain = 4,
    nIter = 8000,
    adapt_delta = 0.99,
    max_treedepth = 11,
    warningFile = sprintf("stanWarnings/exp_%s.txt", modelName)
  )
  
  # if it is the first time to fit the model, fit all participants
  # otherwise, check model fitting results and refit those that fail any of the following criteria
  ## no divergent transitions 
  ## Rhat < 1.01 
  ## Effective Sample Size (ESS) > nChain * 100
  if(!isFirstFit){
    ids = names(trialData)
    paraNames = getParaNames(modelName)
    expPara = loadExpPara(paraNames, outputDir)
    passCheck = checkFit(paraNames, expPara)
    trialData = trialData[!passCheck]
    
    # increase the num of Iterations 
    config = list(
      nChain = 4,
      nIter = 8000,
      adapt_delta = 0.99,
      max_treedepth = 11,
      warningFile = sprintf("stanWarnings/exp_refit_%s.txt", modelName)
    )
  }
  
  # fit the model for all participants
  modelFitGroup(modelName, trialData, config, outputDir, T)
}

