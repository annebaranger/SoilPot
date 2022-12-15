// input data
data {
  int<lower=1> N; // number of observations
  int<lower=1> S; //number of species
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  int <lower=1,upper=S> species [N];
  real prior_breakpoint_fsm;
  real prior_breakpoint_hsm;
  vector[N] fsm; // predictor : frost safety margins
  vector[N] hsm;
}

// the parameters to be estimated from the data
parameters {
  real <lower=-3,upper=1.5> plateau_int_fsm; 
  real <lower=-3,upper=1.5> plateau_int_hsm; 
  vector <lower=-3,upper=1.5> [S] plateau_fsm; // plateau of the second segment, or breakpoint value
  vector <lower=-3,upper=1.5> [S] plateau_hsm;
  real<lower=0> b_fsm; // slope of the first segment, associated to predictor, before breakpoint
  real<lower=0> b_hsm;
  real <lower=min(fsm),upper=max(fsm)> breakpoint_fsm; //point where regression changes
  real <lower=min(hsm),upper=max(hsm)> breakpoint_hsm; 
}

// Functions of estimated parameters.
transformed parameters{
  vector[N] fsm_app; // apparent fsm to catch if predictor is before or after breakpoint
  vector[N] hsm_app;

  for (i in 1:N) {
    if (fsm[i] < breakpoint_fsm) {
      fsm_app[i] = 0;
    } else {
      fsm_app[i] = 1;
    }
  }
  
    for (i in 1:N) {
    if (hsm[i] < breakpoint_hsm) {
      hsm_app[i] = 0;
    } else {
      hsm_app[i] = 1;
    }
  }
}

// The model to be estimated.
// presence are distributed according to Bernoulli law
// logit used as link-function
model {
  //priors
  b_fsm~normal(0,1);
  breakpoint_fsm~normal(prior_breakpoint_fsm,1);
  breakpoint_hsm~normal(prior_breakpoint_hsm,1);
  plateau_fsm~normal(plateau_int_fsm,1);
  plateau_hsm~normal(plateau_int_hsm,1);

  // How the data are expected to have been generated:
  // for (i in 1:N) {
  //   presence[i] ~ bernoulli_logit((1-fsm_app[i]) * (b_fsm * (fsm[i] - breakpoint) + plateau[species[i]])
  //                                 + fsm_app[i] * plateau);
  // }
  presence ~ bernoulli_logit( (1-fsm_app) .* (b_fsm * (fsm - breakpoint_fsm) + plateau_fsm[species]) 
                              + fsm_app .* plateau_fsm[species]
                              + (1-hsm_app) .* (b_hsm * (hsm - breakpoint_hsm) + plateau_hsm[species]) 
                              + hsm_app .* plateau_hsm[species]
  );
}
