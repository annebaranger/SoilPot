// input data
data {
  int<lower=1> N; // number of observations
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  real prior_t_fsm;
  real prior_t_hsm;
  vector[N] fsm; // predictor : frost safety margins
  vector[N] hsm;
}

// the parameters to be estimated from the data
parameters {
  real <lower=-3,upper=1.5> plateau_fsm; 
  real <lower=-3,upper=1.5> plateau_hsm; 
  real<lower=0> b_fsm; // slope of the first segment, associated to predictor, before threshold
  real<lower=0> b_hsm;
  real <lower=min(fsm),upper=max(fsm)> threshold_fsm; //point where regression changes
  real <lower=min(hsm),upper=max(hsm)> threshold_hsm; 
}

// The model to be estimated.
// presence are distributed according to Bernoulli law
// logit used as link-function
model {
  //logical to catch changing point
  vector[N] fsm_app; 
  vector[N] hsm_app;
  for (i in 1:N) {
    if (fsm[i] < threshold_fsm) {
      fsm_app[i] = 0;
    } else {
      fsm_app[i] = 1;
    }
  }
  
    for (i in 1:N) {
    if (hsm[i] < threshold_hsm) {
      hsm_app[i] = 0;
    } else {
      hsm_app[i] = 1;
    }
  }
  
  //priors
  b_fsm~normal(0,1);
  b_hsm~normal(0,1);
  threshold_fsm~normal(prior_t_fsm,1);
  threshold_hsm~normal(prior_t_hsm,1);

  // How the data are expected to have been generated:
  // for (i in 1:N) {
  //   presence[i] ~ bernoulli_logit((1-fsm_app[i]) * (b_fsm * (fsm[i] - threshold) + plateau[species[i]])
  //                                 + fsm_app[i] * plateau);
  // }
  presence ~ bernoulli_logit( (1-fsm_app) .* (b_fsm * (fsm - threshold_fsm) + plateau_fsm) 
                              + fsm_app * plateau_fsm
                              + (1-hsm_app) .* (b_hsm * (hsm - threshold_hsm) + plateau_hsm) 
                              + hsm_app * plateau_hsm
  );
}
