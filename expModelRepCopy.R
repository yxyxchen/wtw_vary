# replicate behavioral data by sumulating with individual fitted parameters
expModelRep = function(modelName){
  # create output directories
  dir.create("figures/expModelRep/")
  dir.create(sprintf("figures/expModelRep/%s",modelName))
  
  # load experiment parameters
  load('expParas.RData')
  
  # load packages and sub functions 
  library("ggplot2") 
  library("dplyr")
  library("tidyr")
  source("subFxs/plotThemes.R")
  source("subFxs/helpFxs.R") 
  source("subFxs/loadFxs.R") # 
  
  
  # normative analysis
  normResults = expSchematics(0, iti, F)
  optimRewardRates = normResults$optimRewardRates
  optimWaitThresholds = normResults$optimWaitThresholds
  presentValues = normResults$presentValues
  
  # get the generative model 
  source(sprintf("subFxs/gnrModels/%s.R", modelName))
  gnrModel = get(modelName)
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  # num of repetitions 
  nRep = 10
  
  # load empirical data 
  allData = loadAllData()
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$id
  nSub = length(ids)

  # initialize outputs
  repTrialData = vector(length = nSub * nRep, mode ='list')
  repNo = matrix(1 : (nSub * nRep), nrow = nRep, ncol = nSub)

  # loop over participants
  for(sIdx in 1 : nSub){
    # prepare empirical data 
    id = ids[[sIdx]]
    thisTrialData = trialData[[id]] 
    # excluded trials at the end of blocks 
    excluedTrials = which(thisTrialData$trialStartTime > (blockSec - max(tMaxs)))
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials &
                                    thisTrialData$blockNum <= 2,]
    # scheduled rewards on non-rewards trials are not recorded and not used in parameter
    # estimation.However, given the scheduled delay, we can infer whether the scheduled 
    # reward is more likely to be 8 or -1. Since scheduled rewards on non-rewards were
    # unknown for participants and didn't influence their choices, we can use the 
    # most likely values as if they were the real scheduled rewards. 
    
    scheduledWait = thisTrialData$scheduledWait
    condition = thisTrialData$condition
    trialEarnings = thisTrialData$trialEarnings
    inferredWait = ifelse(scheduledWait < 6.6195,
           ifelse(condition == "HP", min(tokenValue), max(tokenValue)),
           ifelse(condition == "HP", max(tokenValue), min(tokenValue)))
    scheduledReward = vector(length = length(trialEarnings))
    scheduledReward[trialEarnings == 0] = inferredWait[trialEarnings == 0]
    scheduledReward[trialEarnings != 0] = trialEarnings[trialEarnings != 0]
    
    # load individually fitted paramters 
    fitSummary = read.table(sprintf("genData/expModelFit/%s/s%s_summary.txt",  modelName, id),sep = ",", row.names = NULL)
    paras =  fitSummary[1 : nPara,1]
    # simulate nRep times
    for(rIdx in 1 : nRep){
      tempt = gnrModel(paras, condition, scheduledWait, scheduledReward)
      repTrialData[[repNo[rIdx, sIdx]]] = tempt
    }
  }
  
  ## WTW from empirical data 
  source("MFAnalysis.R")
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  muWTWEmp = sumStats$muWTW
  stdWTWEmp = sumStats$stdWTW
  ## initialize WTW from replicated data 
  muWTWRep_ = matrix(NA, nrow = nRep , ncol = nSub * 2)
  stdWTWRep_ = matrix(NA, nrow = nRep, ncol = nSub * 2)
  ## loop over participants 
  for(sIdx in 1 : nSub){
    id = ids[[sIdx]]
    cbal = hdrData$cbal[sIdx]
    ## loop over conditions
    for(bkIdx in 1 : 2){
      if(hdrData$cbal[sIdx] == 1){
        condition = ifelse(bkIdx == 1, "HP", "LP")
      }else{
        condition = ifelse(bkIdx == 1, "LP", "HP")
      }
      colIdx = sIdx * 2 - 2 + bkIdx # the column index to store results 
      ## loop over repetitions
      for(rIdx in 1 : nRep){
        thisRepTrialData = repTrialData[[repNo[rIdx, sIdx]]]
        thisRepTrialData = as.data.frame(thisRepTrialData[!names(thisRepTrialData) %in% c("Qwaits_", "Viti_")])
        thisRepTrialData = thisRepTrialData[thisRepTrialData$condition == condition,]
        kmscResults = kmsc(thisRepTrialData, min(tMaxs), F, kmGrid)
        muWTWRep_[rIdx,colIdx] = kmscResults$auc
        stdWTWRep_[rIdx, colIdx] = kmscResults$stdWTW
      }
    }
  }
  ## summarise WTW across simulations for replicated data 
  muWTWRep_mu = apply(muWTWRep_, MARGIN = 2, mean) # mean of average willingness to wait
  muWTWRep_std = apply(muWTWRep_, MARGIN = 2, sd) # std of average willingess to wait
  stdWTWRep_mu = apply(stdWTWRep_, MARGIN = 2, mean) # mean of std willingness to wait
  stdWTWRep_std = apply(stdWTWRep_, MARGIN = 2, sd) # std of std willingess to wait
  ## check whether parameter estimates are reliable 
  expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  
  ## save Results
  plotData = data.frame(mu =  muWTWRep_mu, std = stdWTWRep_mu,
                        empMu = muWTWEmp, empStd = stdWTWEmp,
                        passCheck = rep(passCheck, each = 2), 
                        condition = sumStats$condition) %>% filter(passCheck == T)
  save(plotData, file = "genData/varyRep.RData")
  cor.test(plotData$mu, plotData$empMu, method = "spearman")
  ## plot to compare average willingess to wait
    plotData %>%
    ggplot(aes(empMu, mu, shape = condition)) + 
    geom_point(size = 2, color = themeColor, fill = "#ccece6", stroke = 1) + 
    geom_abline(slope = 1, intercept = 0) +
    scale_shape_manual(values = c(21, 23)) + 
    ylab("Model-generated (s)") + xlab("Observed (s)") + ggtitle(sprintf("Average WTW, n = %d", sum(passCheck))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_x_continuous(breaks = c(0, 30), limits = c(-1,31)) + 
    scale_y_continuous(breaks = c(0, 30), limits = c(-1,31))
    
  fileName = sprintf("figures/expModelRep/%s/muWTW_muWTWRep.eps", modelName) 
  ggsave(filename = fileName,  width = 4, height = 4)

  ## plot to compare std willingess to wait
    plotData %>%
    ggplot(aes(empStd, mu, shape = condition)) + 
    geom_point(size = 2, color = themeColor, fill = "#ccece6", stroke = 1)  + 
    geom_abline(slope = 1, intercept = 0) +
    ylab(expression(bold(paste("Model-generated (s"^2,")")))) + xlab(expression(bold(paste("Observed (s"^2,")")))) + ggtitle(sprintf("Std WTW, n = %d", sum(passCheck))) +
    myTheme + theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    scale_shape_manual(values = c(21, 23)) + 
    scale_x_continuous(breaks = c(0, 15)) + 
    scale_y_continuous(breaks = c(0, 15))
  fileName = sprintf("figures/expModelRep/%s/stdWTW_stdWTWRep.eps", modelName)
  ggsave(filename = fileName,  width = 4, height = 4)
}
