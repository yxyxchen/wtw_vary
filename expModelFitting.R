# differences from the cluster version can be found by searching # needed in local computers
expModelFitting = function(modelName){
  # create outfiles
  dir.create("genData") 
  dir.create("genData/expModelFitting") 
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
  
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
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))

  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
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
    thisTrialData = trialData[[thisID]]
    # excluded some trials
    excluedTrials1 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                             thisTrialData$condition == conditions[1])
    excluedTrials2 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                             thisTrialData$condition == conditions[2])
    excluedTrials = c(excluedTrials1, excluedTrials2)
    thisTrialData = thisTrialData[(!(1 : nrow(thisTrialData)) %in% excluedTrials) & thisTrialData$blockNum <= 2,]
    cond = thisTrialData$condition
    scheduledWait = thisTrialData$scheduledWait
    # thisTrialData = block2session(thisTrialData) not needed, since we only use timeWaited and trialEarnings
    fileName = sprintf("genData/expModelFitting/%s/s%s", modelName, thisID)
    modelFitting(thisTrialData, fileName, paraNames, model, modelName)
  }
}
