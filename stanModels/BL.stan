data {
  // experiment parameters
  real stepSec;// duration of one step
  int nStepMax; // max num of steps in one trial
  real iti;// iti duration, unit = second
  
  // initialization parameters
  real reRateIni; 
  
  // experiment data
  int N; // number of trials
  int Rs[N]; // reward in each trial
  int Ts[N]; // terminal state in each trial
}
transformed data {
  // total number of steps in all trials
  // Notiably, T - 1 gives the number of steps in a trial
  int nStepTotal = sum(Ts) - N;
}
parameters {
  // parameters:
  // pWait : probability of waiting at each stepSep
  real<lower = -0.5, upper = 0.5> raw_pWait;
}
transformed parameters{
  // transfer paras
  real pWait = (raw_pWait + 0.5) ; // pWait ~ unif(0, 1)
}
model {
  raw_pWait ~ uniform(-0.5, 0.5);
  
  // calculate the likelihood 
  for(tIdx in 1 : N){
    int action;
    for(i in 1 :(Ts[tIdx] - 1)){
    if(Rs[tIdx] == 0 && i == (Ts[tIdx] - 1)){
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
  vector[nStepTotal] log_lik = rep_vector(0, nStepTotal);
  real LL_all;
  int no = 1;
  // loop over trials
  for(tIdx in 1 : N){
    int action;
    for(i in 1 : (Ts[tIdx] - 1)){
      if(Rs[tIdx] == 0 && i == (Ts[tIdx] - 1)){
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



