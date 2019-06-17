# we don't need lp__ which is a sum loglikelyhood scaled by a constant, something like a dispersion 
# we don't need lp__ for model comparison 
# we save LL_all in both the summary and the all samples data since some times out of memory will change it a lot
# we use log_like to calculate WAIC and looStat
# but we don't save log_like
modelFitting = function(thisTrialData, fileName, paras, model){
  #
  load("wtwSettings.RData")
  # simulation parameters
  nChain = 4
  nIter = 5000
  
  # determine wIni
  # since the participants' initial strategies are unlikely optimal
  # we multiple the optimal opportunity cost by subOptimalRatio
  subOptimalRatio = 0.9 
  if(any(paras == "phiR")){
    wIni = optimRewardRates[[2]] * stepDuration  * subOptimalRatio
  }else if(any(paras == "gamma") || modelName == "baseline"){
    wIni = optimRewardRates[[2]] * stepDuration / (1 - 0.9)  * subOptimalRatio
  }else{
    print("wrong model name!")
    break
  }
  
  # prepare input
  timeWaited = thisTrialData$timeWaited
  scheduledWait = thisTrialData$scheduledWait
  trialEarnings = thisTrialData$trialEarnings
  timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
  cond = unique(thisTrialData$condition)
  tMax = max(tMaxs)
  nTimeSteps = tMax / stepDuration
  Ts = round(ceiling(timeWaited / stepDuration) + 1)
  data_list <- list(tMax = tMax,
                    wIni = wIni,
                    nTimeSteps = nTimeSteps,
                    timeWaited = timeWaited,
                    N = length(timeWaited),
                    trialEarnings = trialEarnings,
                    Ts = Ts)
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
  fitSumary <- summary(fit,pars = c(paras, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSumary, nrow = length(paras) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
            col.names = F, row.names=FALSE)
}

