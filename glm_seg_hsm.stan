// input data
data {
  int<lower=1> N; // number of observations
  int<lower=1> S; //number of species
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  int <lower=1,upper=S> species [N];
  real prior_breakpoint;
  vector[N] fsm; // predictor : frost safety margins
}

// the parameters to be estimated from the data
parameters {
  real <lower=-3,upper=1.5> plateau_int; 
  vector <lower=-3,upper=1.5> [S] plateau; // plateau of the second segment, or breakpoint value
  real<lower=0> b_fsm; // slope of the first segment, associated to predictor, before breakpoint
  real <lower=min(fsm),upper=max(fsm)> breakpoint; //point where regression changes
}

// Functions of estimated parameters.
transformed parameters{
  vector[N] fsm_app; // apparent fsm to catch if predictor is before or after breakpoint
  
  // fsm_app=fsm>breakpoint;
  for (i in 1:N) {
    if (fsm[i] < breakpoint) {
      fsm_app[i] = 0;
    } else {
      fsm_app[i] = 1;
    }
  }
}

// The model to be estimated.
// presence are distributed according to Bernoulli law
// logit used as link-function
model {
  //priors
  b_fsm~normal(0,1);
  breakpoint~normal(prior_breakpoint,1);
  plateau~normal(plateau_int,1);

  // How the data are expected to have been generated:
  // for (i in 1:N) {
  //   presence[i] ~ bernoulli_logit((1-fsm_app[i]) * (b_fsm * (fsm[i] - breakpoint) + plateau[species[i]])
  //                                 + fsm_app[i] * plateau);
  // }
  presence ~ bernoulli_logit( (1-fsm_app) .* (b_fsm * (fsm - breakpoint) + plateau[species]) + fsm_app .* plateau[species]);
}
