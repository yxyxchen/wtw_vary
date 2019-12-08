library("stringr")
getParaNames = function(modelName){
  if(modelName %in% c("QL1", "QL1_prime")) paraNames = c("phi", "tau", "gamma", "prior")
  else if(modelName %in% c("QL2", "QL2_prime")) paraNames = c("phi_pos", "phi_neg", "tau", "gamma", "prior")
  else if(modelName %in% c("RL1", "RL1_prime")) paraNames = c("phi", "tau", "prior", "beta")
  else if(modelName %in% c("RL2", "RL2_prime")) paraNames = c("phi_pos", "phi_neg", "tau", "prior", "beta")
  else if(modelName %in% c("BL")) paraNames = c("pWait")
  return(paraNames)
}

checkFit = function(paraNames, expPara){
  ids = expPara$id
  # detect participants with high Rhats 
  RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(paraNames)] # columns recording Rhats
  if(length(RhatCols) > 1){
    high_Rhat_ids = ids[apply(expPara[,RhatCols] >= 1.01, MARGIN = 1, sum) > 0]
  }else{
    high_Rhat_ids = ids[expPara[,RhatCols] >= 1.01 ]
  }
  
  
  # detect participants with low ESSs
  ESSCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(paraNames)]# columns recording ESSs
  if(length(ESSCols) > 1){
    low_ESS_ids = ids[apply(expPara[,ESSCols] < (4 * 100), MARGIN = 1, sum) > 0]
  }else{
    low_ESS_ids = ids[expPara[,ESSCols] < (4 * 100)]
  }

  
  # detect divergent transitions
  dt_ids = ids[expPara$nDt > 0]
  
  # identify participants satisifying all three criteria:
  passCheck = !ids %in% unique(c(dt_ids, high_Rhat_ids, low_ESS_ids))
  
  return(passCheck)
}


