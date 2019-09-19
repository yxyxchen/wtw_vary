# model-free analysis

# inputs:
# isTrct : logical variable determining whether the last portion in each block is truncated 

# outputs (summarised stats for each participant and each condition, 42 * 2):
# sumStats = {
  # id : [84x1 id]
  # condition : [84x1 fac]
  # nExcl : [84x1 int] # total number of excluded trials 
  # muWTWs : [84x1 num] # average willingness to wait (WTW), measured by area under the Kaplan-Meier survival curve
  # stdWTWs : [84x1 num] # standard deviation of WTW, measured in Kaplan-Meier survival analysis
  # totalEarnings_s :  [84x1 num] 
# }
# timeWTW_ : list(84x1) # wtw timecourse, each element is a vector
# trialWTW_ : list(84x1) # trial-wise WTW, each element is a vector
# survCurve_ : list(84x1) # Kaplan-Meier survival curve, each element is a vector

MFAnalysis = function(isTrct){
  # load libraries
  source('subFxs/loadFxs.R') 
  source('subFxs/analysisFxs.R') 
  library('dplyr')
  library("tidyr")
  
  # create the output directory
  dir.create("genData")
  dir.create("genData/MFAnalysis")
  
  # load experiment parameters
  load("expParas.RData")
  
  # load exp data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id 
  nSub = length(ids)                    # n
  cat('Analyzing data for',nSub,'subjects.\n')
  
  # initialize output variables 
  nExcls = numeric(length = nSub * nBlock)
  muWTWs = numeric(length = nSub * nBlock) 
  stdWTWs = numeric(length = nSub * nBlock) 
  totalEarnings_s =  numeric(length = nSub * nBlock) 
  conditions = numeric(length = nSub * nBlock) 
  timeWTW_ = vector(mode = "list", length = nSub * nBlock) 
  trialWTW_ = vector(mode = "list", length = nSub * nBlock) 
  survCurve_ = vector(mode = "list", length = nSub * nBlock) 
  
  # loop over inidviduals
  for (sIdx in 1 : nSub) {
    # loop over blocks
    for(bkIdx in 1 : nBlock){
      # load trialData 
      id = ids[sIdx]
      thisTrialData = trialData[[id]]
      # index for elements in trialData
      noIdx = (sIdx - 1) * nBlock + bkIdx # 
      # extract (and truncate) trialData for this block
      thisTrialData = thisTrialData %>% filter(thisTrialData$blockNum == bkIdx)
      if(isTrct){
        trctLine = blockSec - max(tMaxs)
        # truncate trials completed after tractline in each block
        nExcls[noIdx] = sum(thisTrialData$trialStartTime > trctLine)
        thisTrialData = thisTrialData %>% filter(trialStartTime <=  trctLine )
      }else{
        nExcls[noIdx] = 0
      }
      
      # determine condition
      conditions[noIdx] = unique(thisTrialData$condition)
      
        
      # calcualte totalEarnings
      totalEarnings_s[noIdx] =  sum(thisTrialData$trialEarnings)
      
      # survival analysis
      kmscResults = kmsc(thisTrialData, min(tMaxs), F, kmGrid)
      
      # 
      muWTWs[noIdx] = kmscResults[['auc']]
      survCurve_[[noIdx]] = kmscResults$kmOnGrid
      stdWTWs[[noIdx]] = kmscResults$stdWTW
      
      # WTW timecourse
      wtwtsResults = wtwTS(thisTrialData, tGrid, min(tMaxs), F)
      timeWTW_[[noIdx]] = wtwtsResults$timeWTW
      trialWTW_[[noIdx]] = wtwtsResults$trialWTW
    }
      
  }
  # return outputs
  sumStats = data.frame(
    id = rep(ids, each = 2),
    condition = conditions,
    nExcl = nExcls,
    totalEarnings = totalEarnings_s,
    muWTW = muWTWs,
    stdWTW = stdWTWs
  )
  outputs = list(
    sumStats = sumStats,
    survCurve_ = survCurve_,
    trialWTW_ = trialWTW_,
    timeWTW_ = timeWTW_ 
  )
  return(outputs)
}
