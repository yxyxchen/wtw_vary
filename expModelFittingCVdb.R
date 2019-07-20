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
  idsCV = 1 : (nSub * nFold) # id encoded in cvPara
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
    if(nFile == length(idsCV)){
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
  while(nLoop < 15){
    cvPara = loadCVPara(paraNames, sprintf("genData/expModelFittingCV/%sdb", modelName),
                        "*_summary.txt")
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
        sIdx = ceiling(excID[i] / nFold)  # ceiling groups 1-10 together yet floor + 1 groups 0-9 together
        fIdx = excID[i] - (sIdx-1) * nFold
        text = sprintf("reFit s%d_f%d", ids[sIdx], fIdx)
        print(text)
        # update nFits and converge
        fitFile = sprintf("genData/expModelFittingCV/%sdb/afit_s%d_f%d.RData", modelName, ids[sIdx], fIdx)
        if(file.exists(fitFile)){
          load(fitFile)
          nFit = nFit  + 1
          save(nFit, file = fitFile)
        }else{
          nFit = 2
          save(nFit, file = fitFile)
        }
        # prepare data
        thisTrialData = trialData[[ids[sIdx]]]
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
        load(sprintf("genData/expModelFittingCV/split/s%d.RData", ids[sIdx]))
        select = as.vector(partTable[-fIdx,])
        thisTrialData = thisTrialData[(1 : nrow(thisTrialData)) %in% select,]
        fileName = sprintf("genData/expModelFittingCV/%sdb/s%d_f%d", modelName,
                           ids[sIdx], fIdx)
        # refit
        # load upper and lower
        tempt = read.csv(sprintf("genData/expModelFittingCV/%sdb/s%d_f%d_summary.txt", modelName,
                                 ids[sIdx], fIdx),header = F)
        low= tempt[1:nPara,4]
        up = tempt[1 : nPara,8]
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




