expModelCmp = function(){
  # libraries and scripts
  library("ggplot2")
  library("dplyr")
  library("tidyr")
  source("subFxs/helpFxs.R")
  source("subFxs/loadFxs.R")
  source("subFxs/plotThemes.R")
  load("expParas.RData")
  
  # load data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id    
  nSub = length(ids) 
  
  # check fit
  modelNames = c("QL1", "QL2", "RL1", "RL2_v2")
  nModel = length(modelNames)
  passCheck_ = matrix(NA, nrow = nSub, ncol = nModel)
  for(i in 1 : nModel){
    modelName = modelNames[i]
    paraNames = getParaNames(modelName)
    expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))
    passCheck_[,i] = checkFit(paraNames, expPara)
  }
  
  # extract leave-one-out results
  logEvidence_ = matrix(NA, nSub, nModel)
  pWaic_ = matrix(NA, nSub, nModel)
  for(m in 1 : nModel){
    modelName = modelNames[m]
    for(sIdx in 1 : nSub ){
      id = ids[sIdx]
      fileName = sprintf("genData/expModelFit/%s/s%s_waic.RData", modelName, id)
      load(fileName)
      logEvidence_[sIdx, m] = LOO$elpd_loo 
      pWaic_[sIdx, m] = LOO$p_loo
    }
  }
  
  # compare model with one learning rates and two seprate learning rates
  library("ggpubr")
  bestNums = sapply(1 : 2, function(i) sum(apply(
    logEvidence_[passCheck_[,1] & passCheck_[,2],1:2],MARGIN = 1, FUN = function(x) which.max(x) == i)))
  data.frame(model = modelNames[1:2], bestNums = bestNums) %>%
    ggplot(aes(x="", y=bestNums, fill=model)) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0) + ylab("") + xlab("") +
    ggtitle(sprintf("Best described (n = %d)", sum(bestNums)))+ 
    myTheme
  dir.create("figures/expModelCmp")
  ggsave("figures/expModelCmp/loo_QL1_QL2.eps", width = 4, height = 3.5)
  
  bestNums = sapply(1 : 2, function(i) sum(apply(
    logEvidence_[passCheck_[,3] & passCheck_[,4],3:4],MARGIN = 1, FUN = function(x) which.max(x) == i)))
  data.frame(model = modelNames[3:4], bestNums = bestNums) %>%
    ggplot(aes(x="", y=bestNums, fill=model)) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0) + ylab("") + xlab("") +
    ggtitle(sprintf("Best described (n = %d)", sum(bestNums)))+ 
    myTheme
  dir.create("figures/expModelCmp")
  ggsave("figures/expModelCmp/loo_RL1_RL2.eps", width = 4, height = 3.5)
  
  # compare RL2 and QL2
  bestNums = sapply(1 : 2, function(i) sum(apply(
    logEvidence_[passCheck_[,2] & passCheck_[,4],c(2,4)],MARGIN = 1, FUN = function(x) which.max(x) == i)))
  data.frame(model = modelNames[c(2,4)], bestNums = bestNums) %>%
    ggplot(aes(x="", y=bestNums, fill=model)) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0) + ylab("") + xlab("") +
    ggtitle(sprintf("Best described (n = %d)", sum(bestNums)))+ 
    myTheme
  dir.create("figures/expModelCmp")
  ggsave("figures/expModelCmp/loo_QL2_RL2.eps", width = 4, height = 3.5)
  
  # sumStats
  MFResults = MFAnalysis(isTrct = T)
  sumStats = MFResults[['sumStats']]
  # save outputs 
  outputs = data.frame(
    logEvidence_,
    passCheck1 = passCheck_[,1],
    passCheck2 = passCheck_[,2],
    passCheck3 = passCheck_[,3],
    passCheck4 = passCheck_[,4],
    muWTW = sumStats$muWTW)
  dir.create("genData/expModelCmp")
  write.table(file = "genData/expModelCmp/loo.csv", outputs,
              sep = ",", col.names = F, row.names = F)
  
}


