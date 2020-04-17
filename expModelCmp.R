expModelCmp = function(){
  # libraries and scripts
  library("ggplot2")
  library("dplyr")
  library("tidyr")
  source("subFxs/helpFxs.R")
  source("subFxs/loadFxs.R")
  source("subFxs/plotThemes.R")
  load("expParas.RData")
  dir.create("figures/expModelCmp")
  
  # load data
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  ids = hdrData$id             
  nSub = length(ids) 
  
  # check fit
  modelNames = c("QL1", "QL2", "RL1", "RL2", "optim_noise", "optim_noise_bias")
  nModel = length(modelNames)
  passCheck_ = matrix(NA, nrow = nSub, ncol = nModel)
  for(i in 1 : nModel){
    modelName = modelNames[i]
    paraNames = getParaNames(modelName)
    expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))
    passCheck_[,i] = checkFit(paraNames, expPara)
  }
  
  # extract leave-one-out results
  waic_ = matrix(NA, nSub, nModel)
  for(m in 1 : nModel){
    modelName = modelNames[m]
    for(sIdx in 1 : nSub ){
      id = ids[sIdx]
      fileName = sprintf("genData/expModelFit/%s/s%s_waic.RData", modelName, id)
      load(fileName)
      waic_[sIdx, m] = WAIC$waic
    }
  }
  
  # 
  outputTable = cbind(waic_, passCheck_)
  write.table(outputTable, "genData/waic.csv", col.names = F, sep = ",")
  
  # plot WAIC
  waic_[!passCheck_] = NA
  data.frame(
    modelName = c("Q1", "Q2", "R1", "R2", "ON", "ONB"),
    mu = apply(waic_, 2, mean, na.rm = T),
    se = apply(waic_, 2, sd, na.rm = T) / sqrt(apply(waic_, 2, function(x) sum(!is.na(x))))
  ) %>% mutate(modelName = factor(modelName, levels = (modelName[order(-mu)]))) %>%
    ggplot(aes(modelName, mu)) + geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mu - se, ymax = mu + se), width = 0.4) +
    myTheme + xlab("") + ylab("WAIC") 
  fileName = "figures/expModelCmp/cmp.eps"
  fileName = "figures/expModelCmp/cmp.png"
  ggsave(filename = fileName,  width = 4, height = 3)
}


