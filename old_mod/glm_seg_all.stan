// input data
data {
  int<lower=1> N; // number of observations
  int<lower=1> S; //number of species
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  int <lower=1,upper=S> species [N];
  real prior_threshold_fsm;
  real prior_threshold_hsm;
  vector[N] fsm; // predictor : frost safety margins
  vector[N] hsm;
}

// the parameters to be estimated from the data
parameters {
  real <lower=-3,upper=1.5> plateau_int_fsm; 
  real <lower=-3,upper=1.5> plateau_int_hsm; 
  vector <lower=-3,upper=1.5> [S] plateau_fsm; // plateau of the second segment, or threshold value
  vector <lower=-3,upper=1.5> [S] plateau_hsm;
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
  threshold_fsm~normal(prior_threshold_fsm,1);
  threshold_hsm~normal(prior_threshold_hsm,1);
  plateau_fsm~normal(plateau_int_fsm,1);
  plateau_hsm~normal(plateau_int_hsm,1);

  // How the data are expected to have been generated:
  // for (i in 1:N) {
  //   presence[i] ~ bernoulli_logit((1-fsm_app[i]) * (b_fsm * (fsm[i] - threshold) + plateau[species[i]])
  //                                 + fsm_app[i] * plateau);
  // }
  presence ~ bernoulli_logit( (1-fsm_app) .* (b_fsm * (fsm - threshold_fsm) + plateau_fsm[species]) 
                              + fsm_app .* plateau_fsm[species]
                              + (1-hsm_app) .* (b_hsm * (hsm - threshold_hsm) + plateau_hsm[species]) 
                              + hsm_app .* plateau_hsm[species]
  );
}

// predictions
generated quantities {
  int presence_pred [N];
  int hsm_app;
  int fsm_app;
  for(i in 1:N) {
    hsm_app=hsm[i]<threshold_hsm?0:1;
    fsm_app=fsm[i]<threshold_fsm?0:1;
    presence_pred[i] = bernoulli_rng(inv_logit(
                              (1-fsm_app) * (b_fsm * (fsm[i] - threshold_fsm) + plateau_fsm[species[i]]) 
                              + fsm_app * plateau_fsm[species[i]]
                              + (1-hsm_app) * (b_hsm * (hsm[i] - threshold_hsm) + plateau_hsm[species[i]]) 
                              + hsm_app * plateau_hsm[species[i]]
                              )
                              );
  }
}

