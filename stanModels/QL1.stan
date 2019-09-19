data {
  // experiment parameters
  real stepSec;// duration of one step
  int nStepMax; // max num of steps in one trial
  real iti;// iti duration
  
  // initialization parameters
  real VitiIni; 
  
  // experiment data
  int N; // number of trials
  int Rs[N]; // reward in each trial
  int Ts[N]; // terminal state in each trial
}
transformed data {
  // total number of steps in all trials
  // Notiably, T - 1 gives the number of steps in a trial
  int nStepTotal = sum(Ts) - N;
  
  // the number of waiting steps in each trial
  int nWait_s[N];
  for(i in 1 : N){
    // if R > 0, the agent always chooses to wait.
    // while if R = 0, the agent quits in the last step
    if(Rs[i] != 0){
      nWait_s[i] = Ts[i] - 1;
    }else{
      nWait_s[i] = Ts[i] - 2;
    }
  }
}
parameters {
  // parameters:
  // phi : learning rate
  // tau : action consistency, namely the soft-max temperature parameter
  // gamma: discount factor
  // prior: prior belief parameter
  
  // for computational efficiency,we sample raw parameters from unif(-0.5, 0.5)
  // which are later transformed into actual parameters
  real<lower = -0.5, upper = 0.5> raw_phi;
  real<lower = -0.5, upper = 0.5> raw_tau;
  real<lower = -0.5, upper = 0.5> raw_gamma;
  real<lower = -0.5, upper = 0.5> raw_prior;
}
transformed parameters{
  // scale raw parameters into real parameters
  real phi = (raw_phi + 0.5) * 0.3; // phi ~ unif(0, 0.3)
  real tau = (raw_tau + 0.5) * 21.9 + 0.1; // tau ~ unif(0.1, 22)
  real gamma = (raw_gamma + 0.5) * 0.3 + 0.7; // gamma ~ unif(0.7, 1)
  real prior = (raw_prior + 0.5) * 65; // prior ~ unif(0, 65)
  
  // declare variables 
  // // state value of the ITI state
  real Viti; 
  // // action value of waiting in each step after ITI
  vector[nStepMax] Qwaits; 
  // // variables to record action values 
  matrix[nStepMax, N] Qwaits_ = rep_matrix(0, nStepMax, N);
  vector[N] Viti_ = rep_vector(0, N);
  // // the reward signal for updating action values at the end of each trial 
  real rwdSignal;  
  
  // initialize action values 
  //// the initial value of the ITI state 
  Viti = VitiIni; 
  // the initial waiting value delines with elapsed time 
  // and the prior parameter determines at which step it falls below Viti
  for(i in 1 : nStepMax){
    Qwaits[i] = (prior - i) * 0.1 + Viti;
  }
  
  // record initial action values
  Qwaits_[,1] = Qwaits;
  Viti_[1] = Viti;
 
  //loop over trials
  for(tIdx in 1 : (N - 1)){
    int T = Ts[tIdx]; // current terminal state
    int R = Rs[tIdx]; // current reward
    
    //calculate the reward signal for updating action value 
    // which equals R plus the discounted value of the successor state. Noticably, 
    // the successor state at the end of trial is always the iti state before the next trial
    rwdSignal = R + gamma * Viti;
    
    // update Qwaits and Viti towards the discounted reward signals 
    for(t in 1 : nWait_s[tIdx]) {
      real discRwdsignal = rwdSignal * gamma^(T - t -1);// the discounted reward signal 
      Qwaits[t] = Qwaits[t] + phi * (discRwdsignal - Qwaits[t]);
    }
    Viti = Viti + phi * (gamma ^ (T -2 + iti / stepSec) * rwdSignal - Viti);
    
    // save action values
    Qwaits_[,tIdx+1] = Qwaits;
    Viti_[tIdx+1] = Viti;
  }
}
model {
  // delcare variables 
  int action; 
  vector[2] actionValues; 
  
  // distributions for raw parameters
  raw_phi ~ uniform(-0.5, 0.5);
  raw_tau ~ uniform(-0.5, 0.5);
  raw_gamma ~ uniform(-0.5, 0.5);
  raw_prior ~ uniform(-0.5, 0.5);
  
  // loop over trials
  for(tIdx in 1 : N){
    int T = Ts[tIdx]; // current terminal state
    int R = Rs[tIdx]; // current reward
    // loop over steps
    for(t in 1 : (T - 1)){
      // determine the action
      // the agent wait in every steps in rewarded trials
      // and wait except for the last step in non-rewarded trials
      if(R == 0 && t == (T-1)){
        action = 2; // proceed to the next ITI
      }else{
        action = 1; // wait
      }
      // calculate the likelihood using the soft-max function
      actionValues[1] = Qwaits_[t, tIdx] * tau;
      actionValues[2] = Viti_[tIdx] * tau;
      target += categorical_logit_lpmf(action | actionValues);
    } 
  }
}
generated quantities {
 // generate action-wise log likelihood and total log likelihood
 
 // initialize variables
  vector[2] actionValues;
  int action;
  vector[nStepTotal] log_lik = rep_vector(0, nStepTotal);
  real LL_all; // total log likelihood
  int no = 1; // action index
  
  // loop over trials
  for(tIdx in 1 : N){
    int T = Ts[tIdx]; // current terminal state
    int R = Rs[tIdx]; // current reward
    // loop over steps
    for(t in 1 : (T - 1)){
      if(R == 0 && t == (T-1)){
        action = 2; // proceed to the next ITI
      }else{
        action = 1; // wait
      }
      // calculate the likelihood using the soft-max function
      actionValues[1] = Qwaits_[t, tIdx] * tau;
      actionValues[2] = Viti_[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | actionValues);
      no = no + 1;
    }
  }
  // calculate total log likelihood
  LL_all =sum(log_lik);
}
