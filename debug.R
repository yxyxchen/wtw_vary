vars = paste0("log_lik[", 1 : 910, "]", sep = "")
a = fit %>%
  rstan::extract(permuted = F, pars = c(vars)) %>%
  adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
  dplyr::select(-chains) 

b = a[1,]

repTrialL = sapply(1 : length(trialEarnings), function(i){
  if(trialEarnings[i] != 0){
    junk = log(lik_[1 : max(Ts[i]-1, 1), i])
    junk[is.infinite(junk)] = -10000
    sum(junk)
  }else{
    junk = c(log(lik_[1:max(Ts[i] - 2,1), i]), log(1-lik_[Ts[i] - 1, i]))
    junk[is.infinite(junk)] = -10000
    sum(junk)
  }
})

fitTrialL = vector(length = nTrial)
no = 1
for(t in 1 : nTrial){
  T = Ts[t]
  fitTrialL[t] = sum(b[no : (no + T - 2)])
  no = no + T - 1
}

compare = cbind(fitTrialL, repTrialL)
wrongTrials = which(round(fitTrialL, 3) != round(repTrialL, 3))
fitTrialL[1] == repTrialL[1]

