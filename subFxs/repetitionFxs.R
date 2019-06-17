# select modelFun by modelName
getRepModelFun = function(modelName){
  if(modelName == "para4"){
    repModelFun = para4
  }else if(modelName %in% c("PR", "PR_5")){
    repModelFun = PR
  }else if(modelName == "curiosityTrialSp"){
    repModelFun = curiosityTrialSp
  }else if(modelName == "baseline"){
    repModelFun = baseLine
  }else{
    return("wrong model name!")
  }
  return(repModelFun)
}

PR = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; gamma = paras[4]; zeroPoint = paras[5]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  subOptimalRatio = 0.9
  QHPApOptim = 5 / 6 * stepDuration / (1 - 0.9) * subOptimalRatio
  QLPApOptim = 0.93 * stepDuration / (1 - 0.9) * subOptimalRatio
  wIni = (QHPApOptim + QLPApOptim)/ 2 
  
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
      nextReward = ifelse(getReward, tokenValue, 0) 
      # dertime whether St+1 is the terminal state
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
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
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis
  )
  return(outputs)
}

################ monte ######################
curiosityTrialSp = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  gamma = paras[3]
  zeroPoint = paras[4]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  # in addition, I use the gamma from 0.5s stepDuration, just hope the Q is similiar to the asympototic value in this RL
  # finally, we use / (1 - gamma) instead of the gamma / (1 - gamma), it assumes the results always happen as the begging 
  # so it is a upper
  # here we use 0.9 as the discount rate for one stepDuration
  QHPApOptim = 5 / 6 * stepDuration / (1 - 0.9) 
  QLPApOptim = 0.93 * stepDuration / (1 - 0.9) 
  wIni = (QHPApOptim + QLPApOptim)/ 2
  
  # Qwait = rep(wIni, nTimeStep)
  # since the participants start the trial with , we assume max(Qwait0) = wini * 0.8
  # again, we assume it has a slope 
  Qquit = wIni * 0.9
  Viti = wIni * 0.9
  #Qwait = rep(wIni*0.93, nTimeStep)
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  #Qwait = rep(wIni, nTimeStep)
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial);
  Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial);
  Qquits[1] = Qquit
  Vitis = vector(length = nTrial);
  Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # initialize elapsed time
  elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine 
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
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether St+1 is the terminal state
      # if the trial terminates, track terminal timestep index T, trialEarnings, timeWaited, sellTime and elapsedTime
      # otherwise, continue
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
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
      
      # update action values for each timestep t
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      # returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      # when the agent always wait and get the reward, update Qwait[1:(T-1)]
      # otherwise, update Qquit and Qwait[1 : (T-2)]      
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)]
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        Qquit = Qquit + phi*(returns[T-1] - Qquit)
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      
      # update Viti
      Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      
      # update Qquit by counterfactual learning
      Qquit = Qquit + phi*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits,
    "Qquits" = Qquits,
    "Gs" = Gs,
    "deltas" = deltas,
    "Vitis" = Vitis
  )
  return(outputs)
}


################ monte ######################
para4 = function(paras, cond, scheduledWait){
  # parse para
  phi = paras[1]
  tau = paras[2]
  gamma = paras[3]
  zeroPoint = paras[4]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize actionValues
  # here we use the optimal reward rates from the normative analysis in Lempert 2018
  # it is more accurate then the one I calcualte in wtwSettings.R
  # in addition, I use the gamma from 0.5s stepDuration, just hope the Q is similiar to the asympototic value in this RL
  # finally, we use / (1 - gamma) instead of the gamma / (1 - gamma), it assumes the results always happen as the begging 
  # so it is a upper
  # here we use 0.9 as the discount rate for one stepDuration
  QHPApOptim = 5 / 6 * stepDuration / (1 - 0.9) 
  QLPApOptim = 0.93 * stepDuration / (1 - 0.9) 
  wIni = (QHPApOptim + QLPApOptim)/ 2
  
  # Qwait = rep(wIni, nTimeStep)
  # since the participants start the trial with , we assume max(Qwait0) = wini * 0.8
  # again, we assume it has a slope 
  Qquit = wIni * 0.9
  Viti = wIni * 0.9
  #Qwait = rep(wIni*0.93, nTimeStep)
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  #Qwait = rep(wIni, nTimeStep)
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial);
  Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial);
  Qquits[1] = Qquit
  Vitis = vector(length = nTrial);
  Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # initialize elapsed time
  elapsedTime = 0
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine 
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
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether St+1 is the terminal state
      # if the trial terminates, track terminal timestep index T, trialEarnings, timeWaited, sellTime and elapsedTime
      # otherwise, continue
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
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
      
      # update action values for each timestep t
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      # returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      # when the agent always wait and get the reward, update Qwait[1:(T-1)]
      # otherwise, update Qquit and Qwait[1 : (T-2)]      
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)]
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
      }else{
        Qquit = Qquit + phi*(returns[T-1] - Qquit)
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phi*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      
      # update Viti
      Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      
      # update Qquit by counterfactual learning
      Qquit = Qquit + phi*(gamma^(iti / stepDuration + 1) * returns[1] - Qquit)
      
      # record updated values
      Qwaits[,tIdx + 1] = Qwait
      Qquits[tIdx + 1] = Qquit
      Vitis[tIdx + 1] = Viti
    }# end of the value update section
    
  } # end of the trial loop
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = scheduledWait,
    "Qwaits" = Qwaits,
    "Qquits" = Qquits,
    "Gs" = Gs,
    "deltas" = deltas,
    "Vitis" = Vitis
  )
  return(outputs)
}

baseline = function(paras, cond, scheduledWait){
  waitRate = paras[1]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  nTimeStep = tMax / stepDuration
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # initialize elapsed time
  elapsedTime = 0
  
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # determine 
    thisScheduledWait = scheduledWait[tIdx]
    # loop for each timestep t and determine At
    t = 1
    while(t <= nTimeStep){
      # determine At
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # observe St+1 and Rt+1
      rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether St+1 is the terminal state
      # if the trial terminates, track terminal timestep index T, trialEarnings, timeWaited, sellTime and elapsedTime
      # otherwise, continue
      nextStateTerminal = (getReward || action == "quit")
      if(nextStateTerminal){
        T = t+1
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
        sellTime[tIdx] = elapsedTime + timeWaited[tIdx] 
        elapsedTime = elapsedTime + timeWaited[tIdx] + iti
        break
      }else{
        t = t + 1
      }
    }# end of the action selection section
  }
}