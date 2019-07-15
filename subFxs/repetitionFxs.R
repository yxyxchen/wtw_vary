# select modelFun by modelName
getRepModelFun = function(modelName){
  if(modelName == "para4"){
    repModelFun = para4
  }else if(modelName %in% c("PR", "PR_5")){
    repModelFun = PR
  }else if(modelName == "PRNC"){
    repModelFun = PRNC
  }else if(modelName == "baseline"){
    repModelFun = baseLine
  }else if(modelName == "uniPrior"){
    repModelFun = uniPrior
  }else if(modelName == "uniPriorNC"){
    repModelFun = uniPriorNC
  }else if(modelName == "hyper"){
    repModelFun = hyper
  }else if(modelName == "PRbs"){
    repModelFun = PRbs
  }else if(modelName == "PRbsNC"){
    repModelFun = PRbsNC
  }else if(modelName == "MVT"){
    repModelFun = MVT
  }else if(modelName %in% c("Rlearn", "Rlearndb")){
    repModelFun = Rlearn
  }else if(modelName == "RlearnL"){
    repModelFun = RlearnL
  }else if(modelName == "reduce_gamma"){
    repModelFun = reduce_gamma
  }else{
    return("wrong model name!")
  }
  return(repModelFun)
}

modelRepitation = function(modelName, summaryData, expTrialData,  nComb){
  paras = getParas(modelName)
  parentDir ="genData/expModelFitting"
  dirName = sprintf("%s/%sdb",parentDir, modelName)
  tempt = loadExpPara(paras, dirName)
  useID = as.factor(tempt$id)
  expPara = tempt
  
  # simulate nRep ztimes for each participants, with different parameter samples
  repModelFun = getRepModelFun(modelName)
  nSub = length(useID)
  repTrialData = vector(length = nSub * nComb, mode ='list')
  repNo = matrix(1 : (nSub * nComb), nrow = nComb, ncol = nSub)
  set.seed(231)
  for(sIdx in 1 : nSub){
    id = useID[[sIdx]]
    # load para samples
    paraSamples = read.table(sprintf("%s/%sdb/s%s.txt", parentDir, modelName, id),sep = ",", row.names = NULL)
    # load behavioral inputs
    thisExpTrialData = expTrialData[[id]] # here we useID
    # truncate 
    # excluded some trials
    excluedTrialsHP = which(thisExpTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                              thisExpTrialData$condition == "HP")
    excluedTrialsLP = which(thisExpTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                              thisExpTrialData$condition == "LP")
    excluedTrials = c(excluedTrialsHP, excluedTrialsLP)
    thisExpTrialData = thisExpTrialData[!(1 : nrow(thisExpTrialData)) %in% excluedTrials,]
    # prepare data
    cond = thisExpTrialData$condition
    scheduledWait = thisExpTrialData$scheduledWait
    scheduledReward = thisExpTrialData$trialEarnings
    scheduledReward[scheduledReward == 0] = ifelse(scheduledWait[scheduledReward == 0] < 6.6195,
                                                   ifelse(cond[scheduledReward == 0] == "Rising", min(tokenValue), max(tokenValue)),
                                                   ifelse(cond[scheduledReward == 0]  == "Rising", max(tokenValue), min(tokenValue)))
    scheduledReward = sapply(scheduledReward, function(x) which(tokenValue == x))
    # simulate
    for(cbIdx in 1 : nComb){
      paraSample = as.double(paraSamples[sample(1 : nrow(paraSamples), 1), 1 : length(paras)])
      tempt = repModelFun(paraSample, cond, scheduledWait, scheduledReward)
      repTrialData[[repNo[cbIdx, sIdx]]] = tempt
    }
  }
  outputs = list(expPara = expPara, useID = useID, repTrialData = repTrialData, repNo = repNo)
  return(outputs)
}



reduce_gamma = function(paras, cond, scheduledWait, scheduledReward){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; zeroPoint = paras[4]
  gamma = 1;
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= tMaxs[2]
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9)  * subOptimalRatio
  
  Qquit = wIni; Viti = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    thisScheduledWait = scheduledWait[tIdx]
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue[scheduledReward[tIdx]], 0) 
      # dertime whether St+1 is the terminal state
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      if(getReward){
        Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }else{
        Viti = Viti + phiP*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }
      
      # update Qquit by counterfactual learning
      if(getReward){
        Qquit = Qquit + phi*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
      }else{
        Qquit = Qquit + phiP*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
      }
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis,
    "cond" = cond
  )
  return(outputs)
}


Rlearn = function(paras, cond, scheduledWait, scheduledReward){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; zeroPoint = paras[4]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= tMaxs[2]
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration  * subOptimalRatio
  
  Qquit = 0; Viti = 0; reRate = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  reRates = vector(length = nTrial); reRates[1] = reRate
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    thisScheduledWait = scheduledWait[tIdx]
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue[scheduledReward[tIdx]], 0) 
      # dertime whether St+1 is the terminal state
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) nextReward - reRate * (T-t) + Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      delta = (returns[1] - reRate * (iti / stepDuration) - Viti)
      if(getReward){
        Viti = Viti + phi * delta
      }else{
        Viti = Viti + phiP* delta
      }
      
      # update Qquit by counterfactual learning
      if(getReward){
        Qquit = Qquit + phi*(returns[1] - reRate * (iti / stepDuration + 1) - Qquit)
      }else{
        Qquit = Qquit + phiP*(returns[1] - reRate * (iti / stepDuration + 1) - Qquit)
      }
      
      # update reRate 
      if(getReward){
        reRate = reRate + phi * delta
      }else{
        reRate = reRate + phiP * delta
      }      
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
      reRates[tIdx + 1] = reRate
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas,
    "Vitis" = Vitis, "reRates" = reRates, "cond" = cond
  )
  return(outputs)
}


RlearnL = function(paras, cond, scheduledWait, scheduledReward){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; zeroPoint = paras[4]
  beta = paras[5]; betaP = paras[6]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= tMaxs[2]
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration  * subOptimalRatio
  
  Qquit = 0; Viti = 0; reRate = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  reRates = vector(length = nTrial); reRates[1] = reRate
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    thisScheduledWait = scheduledWait[tIdx]
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue[scheduledReward[tIdx]], 0) 
      # dertime whether St+1 is the terminal state
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) nextReward - reRate * (T-t) + Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      delta = (returns[1] - reRate * (iti / stepDuration) - Viti)
      if(getReward){
        Viti = Viti + phi * delta
      }else{
        Viti = Viti + phiP* delta
      }
      
      # update Qquit by counterfactual learning
      if(getReward){
        Qquit = Qquit + phi*(returns[1] - reRate * (iti / stepDuration + 1) - Qquit)
      }else{
        Qquit = Qquit + phiP*(returns[1] - reRate * (iti / stepDuration + 1) - Qquit)
      }
      
      # update reRate 
      if(getReward){
        reRate = reRate + beta * delta
      }else{
        reRate = reRate + betaP * delta
      }      
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
      reRates[tIdx + 1] = reRate
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas,
    "Vitis" = Vitis, "reRates" = reRates, "cond" = cond
  )
  return(outputs)
}
################### PRNC
PRNC = function(paras, cond, scheduledWait, scheduledReward){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; gamma = paras[4]; zeroPoint = paras[5]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= tMaxs[2]
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration  * subOptimalRatio
  
  Qquit = wIni; Viti = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    thisScheduledWait = scheduledWait[tIdx]
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue[scheduledReward[tIdx]], 0) 
      # dertime whether St+1 is the terminal state
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        Qquit = Qquit + phiP *(returns[T-1] - Qquit)
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      if(getReward){
        Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }else{
        Viti = Viti + phiP*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis, "cond" = cond
  )
  return(outputs)
}




#### PRbs
PRbsNC = function(paras, cond, scheduledWait, scheduledReward){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; gamma = paras[4]; zeroPoint = paras[5]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= tMaxs[2]
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates) ) * stepDuration / (1 - 0.9)  * subOptimalRatio
  
  Qquit = wIni; Viti = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    thisScheduledWait = scheduledWait[tIdx]
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue[scheduledReward[tIdx]], 0) 
      # dertime whether St+1 is the terminal state
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      if(getReward){
        Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }else{
        Viti = Viti + phiP*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }
      
      # update Qquit by counterfactual learning
      if(tIdx > 1){
        if(trialEarnings[tIdx - 1] == 0){
          if(getReward){
            Qquit = Qquit + phi*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
          }else{
            Qquit = Qquit + phiP*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
          }
        }
      }
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis, "cond" = cond
  )
  return(outputs)
}

PRbs = function(paras, cond, scheduledWait, scheduledReward){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; gamma = paras[4]; zeroPoint = paras[5]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= tMaxs[2]
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9)  * subOptimalRatio
  
  Qquit = wIni; Viti = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial); timeWaited = rep(0, nTrial); sellTime = rep(0, nTrial); elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    thisScheduledWait = scheduledWait[tIdx]
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue[scheduledReward[tIdx]], 0) 
      # dertime whether St+1 is the terminal state
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = nextReward
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
    
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      if(getReward){
        Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }else{
        Viti = Viti + phiP*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }
      
      # update Qquit by counterfactual learning
      if(getReward){
        Qquit = Qquit + phi*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
      }else{
        Qquit = Qquit + phiP*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
      }
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial, "trialEarnings" = trialEarnings, "timeWaited" = timeWaited,
    "sellTime" = sellTime, "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis, "cond" = cond
  )
  return(outputs)
}