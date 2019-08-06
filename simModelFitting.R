getSimTrialData = function(){
  dir.create("genData")
  dir.create("genData/simulation")
  library('ggplot2')
  library('plyr')
  library('dplyr')
  library('tidyr')
  load("wtwSettings.RData")
  source('subFxs/repetitionFxs.R') # called by simulate 
  source("subFxs/helpFxs.R") # getParas
  source("subFxs/plotThemes.R")
  source("subFxs/loadFxs.R") # load scheduledWait from empirical data
  source("subFxs/analysisFxs.R") 
  
  # modelName 
  modelName = "RL2"
  repFun = getRepFun(modelName)
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$ID          
  nSub = length(ids)      
  
  # load expPara
  paraNames = getParaNames(modelName)
  parentDir ="genData/expModelFitting"; dirName = sprintf("%s/%sdb",parentDir, modelName)
  expPara = loadExpPara(paraNames, dirName)
  
  set.seed(123)
  simTrialData = list()
  for(sIdx in 1 : nSub){
    id = ids[sIdx]
    paras = as.double(expPara[expPara$id == id, 1 : length(paraNames)])
    # prepare input
    thisTrialData = trialData[[id]] # here we id instead of sIdx
    # excluded some trials
    excluedTrials1 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                             thisTrialData$condition == conditions[1])
    excluedTrials2 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                             thisTrialData$condition == conditions[2])
    excluedTrials = c(excluedTrials1, excluedTrials2)
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials & 
                                    thisTrialData$blockNum <= 2,]
    condition = thisTrialData$condition
    scheduledWait = thisTrialData$scheduledWait
    blockNum = thisTrialData$blockNum
    scheduledReward = thisTrialData$trialEarnings
    scheduledReward[scheduledReward == 0] = ifelse(scheduledWait[scheduledReward == 0] < 6.6195,
                                                   ifelse(condition[scheduledReward == 0] == "Rising", min(tokenValue), max(tokenValue)),
                                                   ifelse(condition[scheduledReward == 0]  == "Rising", max(tokenValue), min(tokenValue)))
    scheduledReward = sapply(scheduledReward, function(x) which(tokenValue == x))
    id = ids[sIdx]
    simTrialData[[id]] = repFun(paras, condition, scheduledWait, scheduledReward, blockNum)
  }
  hdrData$ID = hdrData$ID
  save(simTrialData, hdrData, file = "genData/simulation/simTrialData.RData")
}
# differences from the cluster version can be found by searching # needed in local computers
expModelFitting = function(modelName){
  # create outfiles
  dir.create("genData") 
  dir.create("genData/simModelFitting") 
  dir.create(sprintf("genData/simModelFitting/%s", modelName))
  
  #load libraries
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');
  library("loo")
  library("coda") 
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getparaNames
  source("subFxs/analysisFxs.R") # for block2session
  load("wtwSettings.RData")
  
  #  set the environment for Rstan
  library('rstan')
  options(warn=-1, message =-1) # run without this for one participant to chec everything
  Sys.setenv(USE_CXX14=1) # needed in local computeres
  rstan_options(auto_write = TRUE) 
  
  # compile the stan model 
  dir.create(sprintf("genData/simModelFitting/%s", modelName))
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))

  # load simData
  load("genData/simulation/simTrialData.RData")
  idList = hdrData$ID
  n = length(idList)             
  
  # determine paraNames
  paraNames = getParaNames(modelName)
  if(paraNames == "wrong model name"){
    print(paraNames)
    break
  }

  # loop over participants 
  library("doMC")
  library("foreach")
  # nCore = as.numeric(Sys.getenv("NSLOTS")) # needed for cluster
  # if(is.na(nCore)) nCore = 1 # needed for cluster
  nCore = parallel::detectCores() -1 # only for the local computer
  registerDoMC(nCore)
  
  foreach(i = 1 : n) %dopar% {
    thisID = idList[[i]]
    thisTrialData = simTrialData[[thisID]]
    cond = thisTrialData$condition
    scheduledWait = thisTrialData$scheduledWait
    # thisTrialData = block2session(thisTrialData) not needed, since we only use timeWaited and trialEarnings
    fileName = sprintf("genData/simModelFitting/%s/s%s", modelName, thisID)
    modelFitting(thisTrialData, fileName, paraNames, model, modelName)
  }
}
