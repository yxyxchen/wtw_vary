expModelRepitation = function(modelName){
  isTrun = T
  library("ggplot2") 
  library("ggpubr")
  library("dplyr")
  library("tidyr")
  source("subFxs/plotThemes.R")
  load("wtwSettings.RData") 
  source("subFxs/helpFxs.R") # getPars
  source("subFxs/loadFxs.R") # load  expPara
  source("subFxs/taskFxs.R") # drawSamples
  source("subFxs/repetitionFxs.R") # getRepFunction
  source("subFxs/analysisFxs.R") # kmsc, trialPlot
  
  # load blockData
  nBlock = 2
  nComb = 10
  load("genData/expDataAnalysis/blockData.RData")
  load("genData/expDataAnalysis/kmOnGridBlock.RData")
  # load trialData since we need scheduledWait 
  allData = loadAllData()
  hdrData = allData$hdrData 
  trialData = allData$trialData       
  ids = hdrData$ID 
  nSub = length(ids)
  
  # re-simulate data
  dir.create("figures/expModelRepitation")
  dir.create(sprintf("figures/expModelRepitation/%s",modelName))
  thisRep = modelRepitation(modelName, summaryData, trialData, nComb) # set seeds indise

  # initialize 
  expPara = thisRep$expPara
  repTrialData = thisRep$repTrialData
  
  paraNames = getParaNames(modelName)
  useID = getUseID(expPara, paraNames)
  repNo = thisRep$repNo
  nSub =(length(useID))
  AUCRep_ = matrix(NA, nrow = nComb , ncol = nSub * nBlock)
  stdWdRep_ = matrix(NA, nrow = nComb, ncol = nSub * nBlock)
  kmOnGridRep_ = vector(mode = "list", length = nSub * nBlock)
  ks_ = matrix(NA, nrow = nComb , ncol = nSub * nBlock)
  dist_ = matrix(NA, nrow = nComb , ncol = nSub * nBlock)
  plotKMSC = F
  
  for(sIdx in 1 : nSub){
    # prepare inputs
    id = useID[[sIdx]]
    label = sprintf("sub%s", id)
    kmOnGridMatrix = matrix(NA, nrow = length(kmGrid), ncol = nComb)
    for(cIdx in 1 : nComb){
      thisRepTrialData = repTrialData[[repNo[cIdx, which(expPara$id == id)]]]
      for(bkIdx in 1 : nBlock){
        noIdx = sIdx * nBlock - nBlock + bkIdx
        thisRepTrialData[c("Qwaits", "reRates","Vitis", "targets", "deltas")] = NULL
        thisRepTrialData = as.data.frame(thisRepTrialData)
        kmscResults = kmsc(thisRepTrialData[thisRepTrialData$blockNum == bkIdx,], min(tMaxs), label ,plotKMSC, kmGrid)
        AUCRep_[cIdx,noIdx] = kmscResults[['auc']]
        stdWdRep_[cIdx, noIdx] = kmscResults$stdWd
        kmOnGridMatrix[,cIdx] = kmscResults$kmOnGrid
        junk = ks.test(kmscResults$kmOnGrid,
                       kmOnGrid_[[which(ids == id) * 2 - 2 +  bkIdx]])
        ks_[cIdx, noIdx] = as.numeric(junk$statistic)
        dist_[cIdx, noIdx] = sum((thisRepTrialData$timeWaited - thisRepTrialData$timeWaited)^2)
      }
    }
    kmOnGridRep_[[noIdx]] = kmOnGridMatrix
  }
  
  # save ks_
  dir.create("genData/expModelRepitation")
  dir.create(sprintf("genData/expModelRepitation/%s", modelName))
  save("dist_", "ks_", file = sprintf("genData/expModelRepitation/%s/compare.RData", modelName))  

  # compare emipirical and reproduced AUC
  muAUCRep = apply(AUCRep_, MARGIN = 2, mean);stdAUCRep = apply(AUCRep_, MARGIN = 2, sd)
  minAUCRep = muAUCRep - stdAUCRep;maxAUCRep = muAUCRep + stdAUCRep
  muStdWdRep = apply(stdWdRep_, MARGIN = 2, mean);stdStdWdRep = apply(stdWdRep_, MARGIN = 2, sd)
  minStdWdRep = muStdWdRep - stdStdWdRep;maxStdWdRep = muStdWdRep + stdStdWdRep
  data.frame(muAUCRep, minAUCRep, maxAUCRep,muStdWdRep, minStdWdRep, maxStdWdRep,
             AUC = blockData$AUC[blockData$id %in% useID], stdWD = blockData$stdWd[blockData$id %in% useID],
             condition = blockData$condition[blockData$id %in% useID],
             blockNum = blockData$blockNum[blockData$id %in% useID]) %>%
    ggplot(aes(AUC, muAUCRep)) +  geom_errorbar(aes(ymin = minAUCRep, ymax = maxAUCRep), color = "grey") +
    geom_point(size = 2) + facet_grid(~blockNum) + 
    geom_abline(slope = 1, intercept = 0) + saveTheme + xlim(c(-4, 34)) + ylim(c(-4, 34)) +
    ylab("Model-generated (s)") + xlab("Observed (s)") + ggtitle(sprintf("Average WTW, n = %d", length(useID))) +
    myThemeBig + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  fileName = sprintf("figures/expModelRepitation/%s/AUC_AUCRep.png", modelName)
  ggsave(filename = fileName,  width = 6, height = 4)
  
  data.frame(muAUCRep, minAUCRep, maxAUCRep,muStdWdRep, minStdWdRep, maxStdWdRep,
             AUC = blockData$AUC[blockData$id %in% useID], stdWD = blockData$stdWd[blockData$id %in% useID],
             condition = blockData$condition[blockData$id %in% useID],
             blockNum = blockData$blockNum[blockData$id %in% useID]) %>%
    ggplot(aes(AUC, muAUCRep)) +  geom_errorbar(aes(ymin = minAUCRep, ymax = maxAUCRep), color = "grey") +
    geom_point(size = 2) + facet_grid(~condition) + 
    geom_abline(slope = 1, intercept = 0) + saveTheme + xlim(c(-4, 34)) + ylim(c(-4, 34)) +
    ylab("Model-generated (s)") + xlab("Observed (s)") + ggtitle(sprintf("Average WTW, n = %d", length(useID))) +
    myThemeBig + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  fileName = sprintf("figures/expModelRepitation/%s/AUC_AUCRep2.png", modelName)
  ggsave(filename = fileName,  width = 6, height = 4)
  
  # try the difference
  data.frame(diff = muAUCRep[blockData$condition == "Falling"] -
               muAUCRep[blockData$condition == "Rising"],
             cbal = factor(blockData$cbal[blockData$condition == "Falling"],
                           labels =c("Rise-Fall", "Fall-Rise")))%>%
    ggplot(aes(cbal, diff)) + geom_boxplot() +
    stat_compare_means(comparisons = list(c("Rise-Fall", "Fall-Rise")),
                       aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                       bracket.size = 1, size = 6) + ylim(c(-30, 30)) +
    xlab("") + ylab("Fall - Rise (s)") + myTheme
  fileName = sprintf("figures/expModelRepitation/%s/optimism.png", modelName)
  ggsave(filename = fileName,  width = 6, height = 4)
  
  # Icompare emipirical and reproduced stdWtW
  data.frame(muAUCRep, minAUCRep, maxAUCRep,muStdWdRep, minStdWdRep, maxStdWdRep,
             AUC = blockData$AUC[blockData$id %in% useID], stdWd = blockData$stdWd[blockData$id %in% useID],
             condition = blockData$condition[blockData$id %in% useID],
             blockNum = blockData$blockNum[blockData$id %in% useID]) %>%
    ggplot(aes(stdWd, muStdWdRep)) + geom_point() + geom_errorbar(aes(ymin = minStdWdRep, ymax = maxStdWdRep), color = "grey") +
    geom_point(size = 2) + facet_grid(~blockNum) + 
    geom_abline(slope = 1, intercept = 0) + saveTheme  +
    ylab(expression(bold(paste("Model-generated (s"^2,")")))) +
    xlab(expression(bold(paste("Observed (s"^2,")")))) +ggtitle(sprintf("Std WTW, n = %d", length(useID)))+
    myThemeBig + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  fileName = sprintf("figures/expModelRepitation/%s/std_stdRep.png", modelName)
  ggsave(filename = fileName,  width = 6, height = 4)
  data.frame(muAUCRep, minAUCRep, maxAUCRep,muStdWdRep, minStdWdRep, maxStdWdRep,
             AUC = blockData$AUC[blockData$id %in% useID], stdWd = blockData$stdWd[blockData$id %in% useID],
             condition = blockData$condition[blockData$id %in% useID],
             blockNum = blockData$blockNum[blockData$id %in% useID]) %>%
    ggplot(aes(stdWd, muStdWdRep)) + geom_point() + geom_errorbar(aes(ymin = minStdWdRep, ymax = maxStdWdRep), color = "grey") +
    geom_point(size = 2) + facet_grid(~condition) + 
    geom_abline(slope = 1, intercept = 0) + saveTheme  +
    ylab(expression(bold(paste("Model-generated (s"^2,")")))) +
    xlab(expression(bold(paste("Observed (s"^2,")")))) +ggtitle(sprintf("Std WTW, n = %d", length(useID)))+
    myThemeBig + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  fileName = sprintf("figures/expModelRepitation/%s/std_stdRep2.png", modelName)
  ggsave(filename = fileName,  width = 6, height = 4)
  
  # # compare emipircal and reproduced trialPlot, for one participant 
  # id = 9
  # sIdx = which(useID  == id)
  # cond = unique(blockData$condition[blockData$id == id])
  # label = sprintf("Sub %d, %s", id, cond)
  # if(isTrun){
  #   junk = lastTrunc(expTrialData[[id]])
  # }
  # junk = block2session(junk)
  # trialPlots(junk, "Observed Data") 
  # ggsave(sprintf("figures/expModelRepitation/%s/actual_data_%d.png", modelName, id),
  #        width = 5, height = 4)
  # 
  # with(thisRep,{
  #   id = 1
  #   sIdx = which(useID  == id)
  #   tempt = repTrialData[[repNo[1,sIdx]]]
  #   tempt$timeWaited =  matrix(unlist(lapply(1:nComb, function(i) repTrialData[[repNo[i,sIdx]]]$timeWaited)), ncol = nComb) %>%
  #     apply(MARGIN  = 1, FUN = mean) 
  #   tempt = within(tempt, sapply(1 : length(timeWaited), function(i) ifelse(timeWaited[i] >= scheduledWait[i], tokenValue, 0)))
  #   tempt$blockNum = junk$blockNum
  #   trialPlots(tempt,"Model-generated Data")
  #   ggsave(sprintf("figures/expModelRepitation/%s/sim_data__%d.png", modelName, id),
  #          width = 5, height = 4)
  # })
  # 
  # 
  # # survival curve prediction
  # idList = hdrData$ID[hdrData$stress == "no stress"]
  # for(sIdx in 1 : nSub){
  #   thisID = idList[sIdx]
  #   if(thisID %in% useID){
  #     para = as.double(expPara[sIdx, 1 : length(paraNames)])
  #     label = sprintf('Subject %s, %s, %s, LL = %.1f',thisID, hdrData$condition[sIdx], hdrData$stress[sIdx], expPara$LL_all[sIdx])
  #     label = paste(label, paste(round(para, 3), collapse = "", seq = " "))
  #     # prepara data 
  #     cond = hdrData$condition[hdrData$ID == thisID]
  #     kmOnGrid = kmOnGrid_[[which(hdrData$ID == thisID)]]
  #     tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  #     kmOnGridRep = apply(kmOnGridRep_[[which(useID== thisID)]], MARGIN = 1, mean)
  #     junk = data.frame(time = kmGrid, exp = kmOnGrid, rep= kmOnGridRep)
  #     plotData = gather(junk, source, survival_rate, -time)
  #     p = ggplot(plotData, aes(time, survival_rate, color = source)) + geom_line(size = 2) + ggtitle(label) + displayTheme
  #     p = p + ylim(c(-0.1,1.1))
  #     print(p)
  #     readline("continue")
  #     fileName = sprintf("%s_s%d.png", modelName, thisID)
  #     #ggsave(fileName, width = 4, height =4)
  #   }
  # }
  # 
}

