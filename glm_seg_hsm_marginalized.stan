// input data
data {
  int<lower=1> N; // number of observations
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  vector[N] fsm; // predictor : frost safety margins
  
  int <lower=1> breakpoint_n;
  real breakpoints [breakpoint_n];
}

transformed data {
  real lambda;
  lambda = -log(breakpoint_n);
}
// the parameters to be estimated from the data
parameters {
  real <lower=-3,upper=1.5> plateau; // plateau of the second segment, or breakpoint value
  real<lower=0> b_fsm; // slope of the first segment, associated to predictor, before breakpoint
}

transformed parameters {
  vector [breakpoint_n] lp;
  lp=rep_vector(lambda,breakpoint_n);
  
  for (bp in 1:breakpoint_n)
    for (n in 1:N)
      lp[bp]=lp[bp]+bernoulli_logit_lpmf(presence[n] | fsm[n] < breakpoints[bp] ? (b_fsm * (fsm[n] - breakpoints[bp]) + plateau) : plateau );
}

// The model to be estimated.
// presence are distributed according to Bernoulli law
// logit used as link-function
model {
  //priors
  b_fsm~normal(0,1);
  target += log_sum_exp(lp);
}
