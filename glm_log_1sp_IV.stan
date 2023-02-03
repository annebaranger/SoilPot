// input data
data {
  int<lower=1> N; // number of observation
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  real prior_K;
  vector[N] fsm; // predictor : frost safety margins
  vector[N] hsm; // predictor : hydraulic safety margins
}

// the parameters to be estimated from the data
parameters {
  real <lower=0,upper=1> K_int;
  real <lower=0> r_fsm;
  real <lower=0> r_hsm;
  real <lower=min(fsm),upper=max(fsm)> t_fsm; //point where regression changes
  real <lower=min(hsm),upper=max(hsm)> t_hsm; //point where regression changes
}

transformed parameters {
  vector <lower=0,upper=1> [N] proba;
  vector <lower=0,upper=1> [N] K_vect;
  real <lower=0> prior_r_fsm;
  real <lower=0> prior_r_hsm;
  K_vect=rep_vector(K_int,N);
  proba = K_vect ./
                (
                  (1 + exp(-r_fsm *(fsm - t_fsm))).*
                  (1 + exp(-r_hsm *(hsm - t_hsm)))
                  );
  prior_r_fsm =  (2*K_int)/(max(fsm)-min(fsm)) ;
  prior_r_hsm =  (2*K_int)/(max(hsm)-min(hsm)) ;
}

model {
  //priors
  K_int~normal(prior_K,0.1);
  t_fsm~normal(0,0.1);
  t_hsm~normal(0,0.1);
  r_fsm~normal(prior_r_fsm,0.1);
  r_hsm~normal(prior_r_hsm,0.1);

  // How the data are expected to have been generated:
  presence ~ bernoulli(proba);
}
