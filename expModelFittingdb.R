# this script fits the RL model for each participant
# using Rstan
# here I change all my modelFitting function into the risk version
# while in stan, I have different expMofelfitting and modelFitting scripts for different things 
# current algorithm might, you know increase nFit yet didn't really refit
expModelFitting = function(modelName){
  #load libraries
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');
  library("stringr")
  library("loo")
  library("coda") 
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getparaNames
  load("wtwSettings.RData")
  source("subFxs/analysisFxs.R")
  
  #  set the environment for Rstan
  library('rstan')
  options(warn=-1, message =-1) # run without this for one participant to chec everything
  Sys.setenv(USE_CXX14=1) # needed in local computeres
  rstan_options(auto_write = TRUE) 
  
  # loop over participants 
  library("doMC")
  library("foreach")
  # nCore = as.numeric(Sys.getenv("NSLOTS")) # needed for cluster
  # if(is.na(nCore)) nCore = 1 # needed for cluster
  nCore = parallel::detectCores() -1 # only for the local computer
  registerDoMC(nCore)
  
  # compile the stan model 
  model = stan_model(file = sprintf("stanModels/%s.stan", paste(modelName, "db", sep = "")))
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$ID          
  nSub = length(ids)                    
  
  originalFile = sprintf("genData/expModelFitting/%s", modelName)
  dbFile = sprintf("genData/expModelFitting/%sdb", modelName)
  if(!file.exists(dbFile)){
    dir.create(dbFile)
    allFiles = list.files(path = originalFile)
    nFile = length(allFiles)
    if(nFile == nSub * 3){
      lapply(1 : nFile, function(i) file.copy(sprintf("%s/%s", originalFile, allFiles[i]),
                                              sprintf("%s/%s", dbFile, allFiles[i])))
      print("creat the debug folder")
    }else{
      print("Wrong number of files in the original folder!")
      break
    }
  }
  
  # determine paraNames
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  if(paraNames == "wrong model name"){
    print(paraNames)
    break
  }
  
  # enter the refit stage
  nLoop = 1
  while(nLoop < 15){
    # determine excID
    expPara = loadExpPara(paraNames,
                          sprintf("genData/expModelFitting/%sdb", modelName))
    useID = getUseID(expPara, paraNames)
    excID = ids[!ids %in% useID]
    
    # loop over excID
    n = length(excID)
    if(n > 0){
      text = sprintf("%s, Start to refit %d participants", modelName, length(excID))
      print(text)
      foreach(i = 1 : n) %dopar% {
        thisID = excID[[i]]
        text = sprintf("refit s%d", thisID)
        print(text)
        # update nFits and converge
        fitFile = sprintf("genData/expModelFitting/%sdb/afit_s%s.RData", modelName, thisID)
        if(file.exists(fitFile)){
          load(fitFile); nFit = nFit  + 1; save(nFit, file = fitFile)
        }else{
          nFit = 2; save(nFit, file = fitFile)
        }
        
        # prepare
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
        # determine fileName
        fileName = sprintf("genData/expModelFitting/%sdb/s%s", modelName, thisID)
        # refit
        # load upper and lower
        tempt = read.csv(sprintf("genData/expModelFitting/%sdb/s%s_summary.txt", modelName, thisID),
                         header = F)
        low= tempt[1:nPara,4]
        up = tempt[1 : nPara,8]
        converge = modelFittingdb(thisTrialData, fileName, paraNames, model, modelName, nPara, low, up)
      }# loop over participants
      nLoop = nLoop + 1
    }else{
      break
    }
  }
  # evaluate useID again
  expPara = loadExpPara(paraNames,
                        sprintf("genData/expModelFitting/%sdb", modelName))
  useID = getUseID(expPara, paraNames)
  print("finish!")
  print(modelName)
  print(length(useID))
}
