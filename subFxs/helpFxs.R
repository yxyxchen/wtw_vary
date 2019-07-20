library("stringr")
getParaNames = function(modelName){
  if(modelName == "QL1") paraNames = c("phi", "tau", "gamma", "prior")
  else if(modelName == "QL2") paraNames = c("phi", "phiP", "tau", "gamma", "prior")
  else if(modelName == "RL1") paraNames = c("phi", "tau", "prior", "beta")
  else if(modelName =="RL2") paraNames = c("phi", "phiP", "tau", "prior", "beta", "betaP")
  else if(modelName == "baseline") paraNames = c("pWait")
  else return("wrong model name")
}

getUseID = function(expPara, paraNames){
  paraNames = c(paraNames, "LL_all")
  idList = expPara$id
  RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(paraNames)]
  EffeCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(paraNames)]
  if(length(RhatCols) == 1){
    useID = idList[expPara[,RhatCols] < 1.1 & 
                     expPara[,EffeCols] >100]
  }else{
    useID = idList[apply(expPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(paraNames) & 
                     apply(expPara[,EffeCols] >100, MARGIN = 1, sum) == length(paraNames)]
  }
  return(useID)
}
  
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