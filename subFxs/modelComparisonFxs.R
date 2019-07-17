# since reward can be -1, we change the defination of getReward
# we also change the rule of using phiP and betaP
# besides, we change tMax and cond
# select modelFun by modelName
getLogLikFun = function(modelName){
  if(modelName == "para4"){
    logLikFun = para4
  }else if(modelName == "baseline"){
    logLikFun = baseLine
  }else if(modelName == "uniPrior"){
    logLikFun = uniPrior
  }else if(modelName == "uniPriorNC"){
    logLikFun = uniPriorNC
  }else if(modelName == "hyper"){
    logLikFun = hyper
  }else if(modelName == "PRbs"){
    logLikFun = PRbs
  }else if(modelName == "PRbsNC"){
    logLikFun = PRbsNC
  }else if(modelName == "MVT"){
    logLikFun = MVT
  }else if(modelName %in% c("Rlearn", "Rlearndb")){
    logLikFun = Rlearn
  }else if(modelName == "RlearnL"){
    logLikFun = RlearnL
  }else if(modelName == "reduce_gamma"){
    logLikFun = reduce_gamma
  }else{
    return("wrong model name!")
  }
  return(logLikFun)
}


PRbs = function(paras, cond, trialEarnings, timeWaited){
  # parse para
  phi = thisParas[1]; phiP = thisParas[2]; tau = thisParas[3]; gamma = thisParas[4]; zeroPoint = thisParas[5]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  
  # initialize actionValues
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration  / (1 - 0.9) * subOptimalRatio
  
  Qquit = wIni; Viti = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  lik_ = matrix(nrow = nTimeStep, ncol = nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    nextReward = trialEarnings[tIdx]
    getReward = ifelse(nextReward != 0, T, F)
    T = Ts[tIdx]
    # calculate logLik
    lik_[,tIdx] =  sapply(1 : nTimeStep, function(i) 1 / sum(1  + exp((Qquit- Qwait[i])* tau)))
    
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        if(nextReward > 0){
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }else{
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phiP*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      if(nextReward > 0){
        Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }else{
        Viti = Viti + phiP*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }
      
      # update Qquit by counterfactual learning
      if(nextReward > 0){
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
    "lik_" = lik_,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis,
    "cond" = cond
  )
  return(outputs)
}

PRbsNC = function(thisParas, cond, trialEarnings, timeWaited){
  # parse para
  phi = thisParas[1]; phiP = thisParas[2]; tau = thisParas[3]; gamma = thisParas[4]; zeroPoint = thisParas[5]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  # initialize actionValues
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9) * subOptimalRatio
  
  Qquit = wIni; Viti = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  lik_ = matrix(nrow = nTimeStep, ncol = nTrial)
    
  # loop over trials
  for(tIdx in 1 : nTrial) {
    nextReward = trialEarnings[tIdx]
    getReward = ifelse(nextReward != 0, T, F)
    T = Ts[tIdx]
    # loop for each timestep t and determine At
    lik_[,tIdx] =  sapply(1 : nTimeStep, function(i) 1 / sum(1  + exp(Qquit - Qwait[i])* tau))
    
    # update values 
    if(tIdx < nTrial){
      nextReward = trialEarnings[tIdx]
      getReward = ifelse(nextReward != 0, T, F)
      T = Ts[tIdx]
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        if(nextReward > 0){
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }else{
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phiP*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      if(nextReward > 0){
        Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }else{
        Viti = Viti + phiP*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }
      
      # update Qquit by counterfactual learning
      if(tIdx > 1){
        if(trialEarnings[tIdx - 1] == 0){
          if(nextReward > 0){
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
    "lik_" = lik_,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis, "cond" = cond
  )
  return(outputs)
}

Rlearn = function(thisParas, cond, trialEarnings, timeWaited){
  # parse para
  phi = thisParas[1]; phiP = thisParas[2]; tau = thisParas[3]; zeroPoint = thisParas[4]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  
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
  lik_ = matrix(nrow = nTimeStep, ncol = nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    nextReward = trialEarnings[tIdx]
    getReward = ifelse(nextReward != 0, T, F)
    T = Ts[tIdx]
    # loop for each timestep t and determine At
    lik_[,tIdx] =  sapply(1 : nTimeStep, function(i) 1 / sum(1  + exp(Qquit - Qwait[i])* tau))
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) nextReward - reRate * (T-t) + Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        if(nextReward > 0){
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }else{
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phiP*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      delta = (returns[1] - reRate * (iti / stepDuration) - Viti)
      if(nextReward > 0){
        Viti = Viti + phi * delta
      }else{
        Viti = Viti + phiP* delta
      }
      
      # update Qquit by counterfactual learning
      if(nextReward > 0){
        Qquit = Qquit + phi*(returns[1] - reRate * (iti / stepDuration + 1) - Qquit)
      }else{
        Qquit = Qquit + phiP*(returns[1] - reRate * (iti / stepDuration + 1) - Qquit)
      }
      
      # update reRate 
      if(nextReward > 0){
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
    "lik_" = lik_,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas,
    "Vitis" = Vitis, "reRates" = reRates, "cond" = cond
  )
  return(outputs)
}


RlearnL = function(thisParas, cond, trialEarnings, timeWaited){
  # parse para
  phi = thisParas[1]; phiP = thisParas[2]; tau = thisParas[3]; zeroPoint = thisParas[4]
  beta = thisParas[5]; betaP = thisParas[6]
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
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
  lik_ = matrix(nrow = nTimeStep, ncol = nTrial)
  
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    nextReward = trialEarnings[tIdx]
    getReward = ifelse(nextReward != 0, T, F)
    T = Ts[tIdx]
    # loop for each timestep t and determine At
    lik_[,tIdx] =  sapply(1 : nTimeStep, function(i) 1 / sum(1  + exp(Qquit - Qwait[i])* tau))
    
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) nextReward - reRate * (T-t) + Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        if(nextReward > 0){
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }else{
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phiP*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      delta = (returns[1] - reRate * (iti / stepDuration) - Viti)
      if(nextReward > 0){
        Viti = Viti + phi * delta
      }else{
        Viti = Viti + phiP* delta
      }
      
      # update Qquit by counterfactual learning
      if(nextReward > 0){
        Qquit = Qquit + phi*(returns[1] - reRate * (iti / stepDuration + 1) - Qquit)
      }else{
        Qquit = Qquit + phiP*(returns[1] - reRate * (iti / stepDuration + 1) - Qquit)
      }
      
      # update reRate 
      if(nextReward > 0){
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
    "lik_" = lik_,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas,
    "Vitis" = Vitis, "reRates" = reRates, "cond" = cond
  )
  return(outputs)
}

reduce_gamma = function(thisParas, cond, trialEarnings, timeWaited){
  # parse para
  phi = thisParas[1]; phiP = thisParas[2]; tau = thisParas[3]; zeroPoint = thisParas[4]
  gamma = 1;
  
  # determine number of trials and nTimeSteps 
  nTrial = length(scheduledWait)
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  # initialize actionValues
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration  / (1 - 0.9) * subOptimalRatio
  
  Qquit = wIni; Viti = wIni 
  Qwait = zeroPoint*0.1 - 0.1*(0 : (nTimeStep - 1)) + Qquit
  
  # initialize varibles for recording action values
  Qwaits = matrix(NA, nTimeStep, nTrial); Qwaits[,1] = Qwait
  Qquits = vector(length = nTrial); Qquits[1] = Qquit
  Vitis = vector(length = nTrial); Vitis[1] = Viti
  deltas = matrix(NA, nTimeStep, nTrial)
  Gs = matrix(NA, nTimeStep, nTrial)
  
  # initialize outputs 
  lik_ = matrix(nrow = nTimeStep, ncol = nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    nextReward = trialEarnings[tIdx]
    getReward = ifelse(nextReward != 0, T, F)
    T = Ts[tIdx]
    # loop for each timestep t and determine At
    lik_[,tIdx] =  sapply(1 : nTimeStep, function(i) 1 / sum(1  + exp(Qquit - Qwait[i])* tau))
  
    
    # update values 
    if(tIdx < nTrial){
      returns = sapply(1 : (T-1), function(t) gamma^(T-t-1) *nextReward + gamma^(T-t) * Viti)
      if(getReward){
        Gs[1 : (T-1), tIdx] = returns[1 : (T-1)];
        deltas[1 : (T-1), tIdx] = returns[1 : (T-1)] - Qwait[1 : (T-1)]
        if(nextReward > 0){
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }else{
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phiP*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }
      }else{
        if(T > 2){
          Gs[1 : (T-2), tIdx] = returns[1 : (T-2)]
          deltas[1 : (T-2), tIdx] = returns[1 : (T-2)] - Qwait[1 : (T-2)]
          Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
        }
      }
      # update Viti
      if(nextReward > 0){
        Viti = Viti + phi*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }else{
        Viti = Viti + phiP*(gamma^(iti / stepDuration) * returns[1] - Viti)
      }
      
      # update Qquit by counterfactual learning
      if(nextReward > 0){
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
    "lik_" = lik_,
    "Qwaits" = Qwaits, "Qquits" = Qquits, "Gs" = Gs, "deltas" = deltas, "Vitis" = Vitis,
    "cond" = cond
  )
  return(outputs)
}
