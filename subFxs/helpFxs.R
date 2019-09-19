library("stringr")
getParaNames = function(modelName){
  if(modelName == "QL1") paraNames = c("phi", "tau", "gamma", "prior")
  else if(modelName == "QL2") paraNames = c("phi_pos", "phi_neg", "tau", "gamma", "prior")
  else if(modelName == "RL1") paraNames = c("phi", "tau", "prior", "beta")
  else if(modelName =="RL2") paraNames = c("phi_pos", "phi_neg", "tau", "prior", "beta")
  else if(modelName == "BL") paraNames = c("pwait")
  return(paraNames)
}

checkFit = function(paraNames, expPara){
  ids = expPara$id
  # detect participants with high Rhats 
  RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(paraNames)] # columns recording Rhats
  high_Rhat_ids = ids[apply(expPara[,RhatCols] >= 1.01, MARGIN = 1, sum) > 0]
  
  # detect participants with low ESSs
  ESSCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(paraNames)]# columns recording ESSs
  low_ESS_ids = ids[apply(expPara[,ESSCols] < (4 * 100), MARGIN = 1, sum) > 0]
  
  # detect divergent transitions
  dt_ids = ids[expPara$nDt > 0]
  
  # identify participants satisifying all three criteria:
  passCheck = !ids %in% unique(c(dt_ids, high_Rhat_ids, low_ESS_ids))
  
  return(passCheck)
}


# check R-hat and ESS of parameter estimations
## 
getParaComb = function(paraTable){
  paraNames = names(paraTable)
  nPara = length(paraTable)
  nValue = unlist(lapply(1 : nPara, function(i) length(paraTable[[i]])))
  
  output = matrix(NA, nrow = prod(nValue), ncol = nPara)
  for(pIdx in 1 : nPara){
    eachRep = ifelse(pIdx >= nPara, 1, prod(nValue[(pIdx+1) : nPara])) # repetition number for each element in the seq
    seqRep = ifelse(pIdx <= 1, 1, prod(nValue[1 : (pIdx - 1)])) # repetition number for the seq
    output[,pIdx]= rep(rep(paraTable[[pIdx]], each = eachRep), seqRep)
  }
  colnames(output) = paraNames
  return(output)
}



softMax = function(values, tau){
  probAction1 = 1 / (1 + exp((values[2] -values[1]) * tau))
  return(probAction1)
}

getWD = function(pWaits, stepDuration){
  nStep = length(pWaits)
  survProbs = cumprod(pWaits) # survProbs[i] is the prob that survive after the ith steps
  pAlwaysWait = survProbs[nStep]
  quitProbs = c(1, survProbs[1:(nStep-1)]) - survProbs # quitProbs[i], prob that quit before the ith steps
  waitDurations = c(1 : nStep * stepDuration, nStep * stepDuration)
  probWd = c(quitProbs, pAlwaysWait)
  muWd = sum(waitDurations * probWd) 
  stdWd = sqrt(sum((waitDurations - muWd)^2 * probWd))
  outputs = list(waitDurations = waitDurations,
                 probs = probWd,
                 mu = muWd,
                 std = stdWd)
  return(outputs)
}