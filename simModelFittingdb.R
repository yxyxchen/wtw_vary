expModelFitting = function(decodeModel){
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
  rstan_options(auto_write = TRUE) 
  
  # loop over participants 
  library("doMC")
  library("foreach")
  nCore = as.numeric(Sys.getenv("NSLOTS")) # needed for cluster
  if(is.na(nCore)) nCore = 1 # needed for cluster
  # nCore = parallel::detectCores() -1 # only for the local computer
  registerDoMC(nCore)
  
  # compile the stan model 
  model = stan_model(file = sprintf("stanModels/%s.stan", paste(decodeModel, "db", sep = "")))
  
  # load simData
  load(sprintf("genData/simulation/%s.RData", encodeModel))
  ids = hdrData$ID
  nSub = length(ids)                    
  
  originalFile = sprintf("genData/simModelFitting/%s/%s", encodeModel, decodeModel)
  dbFile = sprintf("genData/simModelFitting/%s/%sdb", encodeModel, decodeModel)
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
  paraNames = getParaNames(decodeModel)
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
                          sprintf("genData/simModelFitting/%sdb", decodeModel))
    useID = getUseID(expPara, paraNames)
    excID = ids[!ids %in% useID]
    
    # loop over excID
    n = length(excID)
    if(n > 0){
      text = sprintf("Start to refit %d participants", length(excID))
      print(text)
      foreach(i = 1 : n) %dopar% {
        thisID = excID[[i]]
        text = sprintf("refit s%s", thisID)
        print(text)
        # update nFits and converge
        fitFile = sprintf("genData/simModelFitting/%s/%sdb/afit_s%s.RData", encodeModel, decodeModel, thisID)
        if(file.exists(fitFile)){
          load(fitFile); nFit = nFit  + 1; save(nFit, file = fitFile)
        }else{
          nFit = 2; save(nFit, file = fitFile)
        }
        
        # prepare
        thisTrialData = simTrialData[[thisID]]
        cond = thisTrialData$condition
        scheduledWait = thisTrialData$scheduledWait
        # determine fileName
        fileName = sprintf("genData/simModelFitting/%s/%sdb/s%s", encodeModel, decodeModel, thisID)
        # refit
        # load upper and lower
        tempt = read.csv(sprintf("genData/simModelFitting/%s/%sdb/s%s_summary.txt", encodeModel, decodeModel, thisID),
                         header = F)
        low= tempt[1:nPara,4]
        up = tempt[1 : nPara,8]
        converge = modelFittingdb(thisTrialData, fileName, paraNames, model, decodeModel, nPara, low, up)
      }# loop over participants
      nLoop = nLoop + 1
    }else{
      break
    }
  }
  # evaluate useID again
  expPara = loadExpPara(paraNames,
                        sprintf("genData/simModelFitting/%s/%sdb", encodeModel, decodeModel))
  useID = getUseID(expPara, paraNames)
  print("finish")
  print(decodeModel)
  print(length(useID))
}
