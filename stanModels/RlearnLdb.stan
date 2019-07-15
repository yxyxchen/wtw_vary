data {
  // depending on the condition
  real wIni;
  int tMax;
  int nTimeSteps; // nTimeSteps = tMax / stepDuration
  int nPara;
  
  // depending on each subject
  int N; // number of trials
  vector[N] timeWaited;
  vector[N] trialEarnings;
  int Ts[N]; // terminal time step index 
  real stepDuration;
  real iti;
  real tokenValue;
  vector[nPara] low;
  vector[nPara] up; 
}
transformed data {
  int totalSteps = sum(Ts) - N;
}
parameters {
  real<lower = low[1], upper = up[1]> phi;
  real<lower = low[2], upper = up[2]> phiP; 
  real<lower = low[3], upper = up[3]> tau;
  real<lower = low[4], upper = up[4]> zeroPoint; 
  real<lower = low[5], upper = up[5]> beta;
  real<lower = low[6], upper = up[6]> betaP;   
}
transformed parameters{
  // initialize action values 
  // especially for this version use 0.9, original 1
  real Qquit = 0;
  real Viti = 0;
  real reRate = wIni;
  vector[nTimeSteps] Qwait;
    // initialize variables to record action values 
  matrix[nTimeSteps, N] Qwaits = rep_matrix(0, nTimeSteps, N);
  vector[N] Qquits = rep_vector(0, N);

  // initialize caching variables
  real G1;
  real delta;
  // fill values
  for(i in 1 : nTimeSteps){
    Qwait[i] = zeroPoint*0.1 - 0.1*(i - 1) + Qquit;
  }
  
  // fill the first element of Qwaits, Quits and Vitis 
  Qwaits[,1] = Qwait;
  Qquits[1] = Qquit;
  
  //loop over trial
  for(tIdx in 1 : (N -1)){
    // determine the termial timestep T 
    int T = Ts[tIdx];
    real RT = trialEarnings[tIdx];
    
    // update action values for rewarded trials
    if(RT > 0){
      for(t in 1 : (T - 1)){
        real G = RT - reRate * (T - t) + Viti;
        Qwait[t] = Qwait[t] + phi * (G - Qwait[t]);
      }
    }else{
      if(T > 2){
        for(t in 1 : (T-2)){
          real G =  RT  - reRate * (T - t) + Viti;
          Qwait[t] = Qwait[t] + phiP * (G - Qwait[t]);    
        }
      }
    }
    // update Qquit by counterfactual thiking
    G1 =  RT  - reRate*(T - 1) + Viti;
    if(RT > 0){
      Qquit = Qquit + phi * (G1 - reRate * (iti / stepDuration + 1) - Qquit);
    }else{
      Qquit = Qquit + phiP * (G1 - reRate * (iti / stepDuration + 1) - Qquit);
    }
    
    // update Viti
    delta = (G1 - reRate * (iti / stepDuration) - Viti);
    if(RT > 0){
       Viti = Viti + phi * delta;
    }else{
       Viti = Viti + phiP * delta;
    }
   
    // update reRate 
    if(RT > 0){
      reRate = reRate +  beta * delta;
    }else{
      reRate = reRate + betaP * delta;
    }
    
    // save action values
    Qwaits[,tIdx+1] = Qwait;
    Qquits[tIdx+1] = Qquit;
  }// end of the loop
}
model {
  phi ~ uniform(low[1], up[1]);
  phiP ~ uniform(low[2], up[2]);
  tau ~ uniform(low[3], up[3]);
  zeroPoint ~ uniform(low[4], up[4]);
  beta ~ uniform(low[5], up[5]);
  betaP ~ uniform(low[6], up[6]);  
  // calculate the likelihood 
  for(tIdx in 1 : N){
    int action;
    vector[2] values;
    int T = Ts[tIdx];
    for(i in 1 : (T - 1)){
      if(trialEarnings[tIdx] == 0 && i == (T-1)){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = Qwaits[i, tIdx] * tau;
      values[2] = Qquits[tIdx] * tau;
      //action ~ categorical_logit(values);
      target += categorical_logit_lpmf(action | values);
    } 
  }
}
generated quantities {
// initialize log_lik
  vector[totalSteps] log_lik = rep_vector(0, totalSteps);
  vector[2] values;
  real LL_all;
  int no = 1;
  // loop over trials
  for(tIdx in 1 : N){
    int action;
    int T = Ts[tIdx];
    for(i in 1 : (T - 1)){
      if(trialEarnings[tIdx] == 0 && i == (T-1)){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = Qwaits[i, tIdx] * tau;
      values[2] = Qquits[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | values);
      no = no + 1;
    }
  }// end of the loop
  LL_all =sum(log_lik);
}
