// input data
data {
  int<lower=1> N; // number of observation
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  real prior_K;
  vector[N] fsm; // predictor : frost safety margins
}

// the parameters to be estimated from the data
parameters {
  real <lower=0,upper=1> K_int;
  real <lower=0> r_fsm;
  real <lower=min(fsm),upper=max(fsm)> t_fsm; //point where regression changes
}

transformed parameters {
  vector <lower=0,upper=1> [N] proba;
  vector <lower=0,upper=1> [N] K_vect;
  K_vect=rep_vector(K_int,N);
  proba = K_vect ./
                (
                  (1 + exp(-r_fsm *(fsm - t_fsm)))
                  );

}

model {
  //priors
  K_int~normal(prior_K,0.1);
  t_fsm~normal(0,0.1);

  // How the data are expected to have been generated:
  presence ~ bernoulli(proba);
}
