# we don't need lp__ which is a sum loglikelyhood scaled by a constant, something like a dispersion 
# we don't need lp__ for model comparison 
# we save LL_all in both the summary and the all samples data since some times out of memory will change it a lot
# we use log_like to calculate WAIC and looStat
# but we don't save log_like
modelFitting = function(thisTrialData, fileName, paras, model, modelName){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni
  # since the participants' initial strategies are unlikely optimal
  # we multiple the optimal opportunity cost by subOptimalRatio
  subOptimalRatio = 0.9 
  if(modelName %in% c("baseline", "MVT", "Rlearn", "RlearnL")){
    wIni = mean(as.double(optimRewardRates)) * stepDuration  * subOptimalRatio
  }else if(any(paras %in% c("gamma", "k")) || modelName == "reduce_gamma"){
    wIni =  mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9)  * subOptimalRatio
  }else{
    print("wrong model name!")
    break
  }
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings != loseValue] = scheduledWait[trialEarnings != loseValue]
  tMax = tMaxs[2] # the second change
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(tMax = tMax,
                    wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    # real data
                    timeWaited = timeWaited,
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    #
                    iti = iti,
                    stepDuration = stepDuration,
                    tokenValue = tokenValue)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  # extract parameters
  extractedPara = fit %>%
    rstan::extract(permuted = F, pars = c(paras, "LL_all"))
  # save sampling sequences
  tempt = extractedPara %>%
    adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
    dplyr::select(-chains) 
  write.table(matrix(unlist(tempt), ncol = length(paras) + 1), file = sprintf("%s.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE) 
  # calculate and save WAIC
  log_lik = extract_log_lik(fit) # quit time consuming
  WAIC = waic(log_lik)
  looStat = loo(log_lik)
  save("WAIC", "looStat", file = sprintf("%s_waic.RData", fileName))
  fitSummary <- summary(fit,pars = c(paras, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paras) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
}

modelFittingCV = function(thisTrialData, fileName, paras, model, modelName){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni
  # since the participants' initial strategies are unlikely optimal
  # we multiple the optimal opportunity cost by subOptimalRatio
  subOptimalRatio = 0.9 
  if(modelName %in% c("baseline", "MVT", "Rlearn", "RlearnL")){
    wIni = mean(as.double(optimRewardRates)) * stepDuration  * subOptimalRatio
  }else if(any(paras %in% c("gamma", "k")) || modelName == "reduce_gamma"){
    wIni =  mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9)  * subOptimalRatio
  }else{
    print("wrong model name!")
    break
  }
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings != loseValue] = scheduledWait[trialEarnings != loseValue]
  tMax = tMaxs[2] # the second change
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(tMax = tMax,
                    wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    # real data
                    timeWaited = timeWaited,
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    #
                    iti = iti,
                    stepDuration = stepDuration,
                    tokenValue = tokenValue)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  fitSummary <- summary(fit,pars = c(paras, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paras) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
}

modelFittingdb = function(thisTrialData, fileName, paras, model, modelName,nPara, low, up){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni
  # since the participants' initial strategies are unlikely optimal
  # we multiple the optimal opportunity cost by subOptimalRatio
  subOptimalRatio = 0.9 
  if(modelName %in% c("baseline", "MVT", "Rlearn", "RlearnL")){
    wIni = mean(as.double(optimRewardRates)) * stepDuration  * subOptimalRatio
  }else if(any(paras %in% c("gamma", "k")) || modelName == "reduce_gamma"){
    wIni =  mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9)  * subOptimalRatio
  }else{
    print("wrong model name!")
    break
  }
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings != loseValue] = scheduledWait[trialEarnings != loseValue]
  tMax = tMaxs[2] # the second change
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(tMax = tMax,
                    wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    nPara = nPara,
                    # real data
                    timeWaited = timeWaited,
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    low = low,
                    up = up,
                    #
                    iti = iti,
                    stepDuration = stepDuration,
                    tokenValue = tokenValue)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  # extract parameters
  extractedPara = fit %>%
    rstan::extract(permuted = F, pars = c(paras, "LL_all"))
  # save sampling sequences
  tempt = extractedPara %>%
    adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
    dplyr::select(-chains) 
  write.table(matrix(unlist(tempt), ncol = length(paras) + 1), file = sprintf("%s.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE) 
  # calculate and save WAIC
  log_lik = extract_log_lik(fit) # quit time consuming
  WAIC = waic(log_lik)
  looStat = loo(log_lik)
  save("WAIC", "looStat", file = sprintf("%s_waic.RData", fileName))
  fitSummary <- summary(fit,pars = c(paras, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paras) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)

  # detmerine converge
  converge = all(fitSummary[,"Rhat"] < 1.1) & all(fitSummary[, "n_eff"] >100)
  return(converge)  
}


modelFittingCVdb = function(thisTrialData, fileName, paras, model, modelName,nPara, low, up){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni
  # since the participants' initial strategies are unlikely optimal
  # we multiple the optimal opportunity cost by subOptimalRatio
  subOptimalRatio = 0.9 
  if(modelName %in% c("baseline", "MVT", "Rlearn", "RlearnL")){
    wIni = mean(as.double(optimRewardRates)) * stepDuration  * subOptimalRatio
  }else if(any(paras %in% c("gamma", "k")) || modelName == "reduce_gamma"){
    wIni =  mean(as.double(optimRewardRates)) * stepDuration / (1 - 0.9)  * subOptimalRatio
  }else{
    print("wrong model name!")
    break
  }
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings != loseValue] = scheduledWait[trialEarnings != loseValue]
  tMax = tMaxs[2] # the second change
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(tMax = tMax,
                    wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    nPara = nPara,
                    # real data
                    timeWaited = timeWaited,
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts,
                    low = low,
                    up = up,
                    #
                    iti = iti,
                    stepDuration = stepDuration,
                    tokenValue = tokenValue)
  fit = sampling(object = model, data = data_list, cores = 1, chains = nChain,
                 iter = nIter) 
  fitSummary <- summary(fit,pars = c(paras, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paras) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
  
  # detmerine converge
  converge = all(fitSummary[,"Rhat"] < 1.1) & all(fitSummary[, "n_eff"] >100)
  return(converge)  
}