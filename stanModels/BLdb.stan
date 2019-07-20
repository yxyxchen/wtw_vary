data {
  // depending on the task
  real stepDuration;
  real iti;
  
  // depending on the condition
  real wIni;
  int nTimeSteps; // nTimeSteps = tMax / stepDuration
  
  // depending on each subject
  int N; // number of trials
  vector[N] trialEarnings;
  int Ts[N]; // terminal time step index 
  
  // debug
  int nPara;
  vector[nPara] low;
  vector[nPara] up;
}
transformed data {
  int totalSteps = sum(Ts) - N;
}
parameters {
  real<lower = low[1], upper = up[1]> pWait;
}
model {
  pWait ~ uniform(low[1], up[1]);
  // calculate the likelihood 
  for(tIdx in 1 : N){
    int action;
    for(i in 1 :(Ts[tIdx] - 1)){
    if(trialEarnings[tIdx] == 0 && i == (Ts[tIdx] - 1)){
      action = 0; // quit
    }else{
      action = 1; // wait
    }
      target += bernoulli_lpmf(action | pWait);// theta defines the prob of 1
    } 
  }
}
generated quantities {
// initialize log_lik
  vector[totalSteps] log_lik = rep_vector(0, totalSteps);
  real LL_all;
  int no = 1;
  // loop over trials
  for(tIdx in 1 : N){
    int action;
    for(i in 1 : (Ts[tIdx] - 1)){
      if(trialEarnings[tIdx] == 0 && i == (Ts[tIdx] - 1)){
        action = 0; // quit
      }else{
        action = 1; // wait
      }
      log_lik[no] = bernoulli_lpmf(action | pWait);
      no = no + 1;
    }
  }// end of the loop
  LL_all =sum(log_lik);
}



