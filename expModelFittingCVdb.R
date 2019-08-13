# if there is an error, especially when the errors occurs after modelFitting
# nFits will not match summary.txt
# in expModelFittingCV.R, ids in useID, expPara, converges are 1,2,3,
# others, like all ids in the filename, are like original ids
expModelFitting = function(modelName){
  
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library("stringr")
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getparaNames
  load("wtwSettings.RData")
  
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
  
  # parameters
  nFold = 10
  
  # load expData
  # one sub's data, two conditions together, are in one element in trialData
  # hdrData only has 42 entries
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  nSub = nrow(hdrData)
  ids = hdrData$ID# id encoded in trialData
  # initialize outputs
  
  # for a specific model 
  
  # detect the debug folder
  originalFile = sprintf("genData/expModelFittingCV/%s", modelName)
  dbFile = sprintf("genData/expModelFittingCV/%sdb", modelName)
  if(!file.exists(dbFile)){
    dir.create(dbFile)
    allFiles = list.files(path = originalFile)
    nFile = length(allFiles)
    if(nFile == (nSub * nFold)){
      lapply(1 : nFile, function(i) file.copy(sprintf("%s/%s", originalFile, allFiles[i]),
                                              sprintf("%s/%s", dbFile, allFiles[i])))
      print("creat the debug folder")
    }else{
      print("Wrong number of files in the original folder!")
      break
    }
  }
  
  # loop over models
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # load cvPara
  # enter the refit stage
  nLoop = 1
  while(nLoop < 2){
    cvPara = loadCVPara(paraNames, sprintf("genData/expModelFittingCV/%sdb", modelName),
                        "*_summary.txt")
    idsCV = cvPara$id
    useID = getUseID(cvPara, paraNames)
    excID = idsCV[!idsCV %in% useID]
    
    # refit the mode
    if(length(excID) > 0){
      text = sprintf("%s, Start to refit %d participants", modelName, length(excID))
      print(text)
      # compile the debug version of the model
      model = stan_model(file = sprintf("stanModels/%sdb.stan", modelName))
      foreach(i = 1 : length(excID)) %dopar% {
        # extract sIdx and fIdx from the id encoded in cvPara
        id = excID[i]
        junk = str_locate(id, "s[0-9]+")
        sIdx = substr(id, (junk[1] + 1), junk[2]) # use to load fit.RData and trialData
        junk = str_locate(id, "f[0-9]+")
        fIdx =  as.double(substr(id, (junk[1] + 1), junk[2]))
        text = sprintf("reFit %s", id)
        print(text)
        # update nFits and converge
        fitFile = sprintf("genData/expModelFittingCV/%sdb/afit_%s.RData", modelName, id)
        if(file.exists(fitFile)){
          load(fitFile)
          nFit = nFit  + 1
          save(nFit, file = fitFile)
        }else{
          nFit = 2
          save(nFit, file = fitFile)
        }
        # prepare data
        thisTrialData = trialData[[sIdx]]
        # excluded some trials
        excluedTrials1 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                                 thisTrialData$condition == conditions[1])
        excluedTrials2 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                                 thisTrialData$condition == conditions[2])
        excluedTrials = c(excluedTrials1, excluedTrials2)
        thisTrialData = thisTrialData[(!(1 : nrow(thisTrialData)) %in% excluedTrials) & thisTrialData$blockNum <= 2,]
        cond = thisTrialData$condition
        scheduledWait = thisTrialData$scheduledWait
        # select the training set
        load(sprintf("genData/expModelFittingCV/split/s%s.RData", sIdx))
        select = c(1:5, as.vector(partTable[-fIdx,]))
        thisTrialData = thisTrialData[(1 : nrow(thisTrialData)) %in% select,]
        fileName = sprintf("genData/expModelFittingCV/%sdb/%s", modelName,
                           id)
        # refit
        # load upper and lower
        tempt = read.csv(sprintf("genData/expModelFittingCV/%sdb/%s_summary.txt", modelName,
                                 id),header = F)
        low= tempt[1:nPara,4]
        up = tempt[1 : nPara,8]
        lowLimits = sapply(1 : nPara, function(i) ifelse(paraNames[i] == "gamma", 0.7,ifelse(paraNames[i] == "tau", 0.1, 0)))
        upMatchTable = c(0.3, 5, 22, 1, 0.3, 65);
        names(upMatchTable) = c("phi", "nega", "tau", "gamma", "beta", "prior")
        upLimits = sapply(1 : nPara, function(i){
          paraName = paraNames[i]
          as.double(upMatchTable[paraName])})
        if(any(abs(low - up) < 0.03)){
          closeIdxs = which(abs(low - up) < 0.01)
          low[closeIdxs] = pmax(lowLimits[closeIdxs], low[closeIdxs] - 0.01)
          up[closeIdxs] = pmin(upLimits[closeIdxs], up[closeIdxs] + 0.01)
        }
        converge = modelFittingCVdb(thisTrialData, fileName, paraNames, model, modelName, nPara, low, up)
      }# loop over participants  
      nLoop = nLoop + 1
    }else{
      break
    }
  }
  # evaluate useID again
  cvPara = loadCVPara(paraNames, sprintf("genData/expModelFittingCV/%sdb", modelName),"*_summary.txt")
  useID = getUseID(cvPara, paraNames)
  print("finish!")
  print(modelName)
  print(length(useID))
}# end of the function




