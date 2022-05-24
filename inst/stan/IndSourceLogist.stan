data {
  int<lower=1> m;                     // observations
  int p;                              // columns (#covariates + intercept) in model matrix
  int<lower=0, upper=1> y[m];         // dv ; binary
  matrix[m, p] X;                     // model design matrix
  real<lower=0> tau;                  // mu prior hyperparameter
  real<lower=0> a_omega;              // omega prior shape hyperparameter
  real<lower=0> b_omega;              // omega prior scale hyperparameter
  vector[p] eta;                      // betas mean
  matrix[p, p] I;                     // betas covariance matrix structure ; benchmark identity divided by omega
}

parameters {
  row_vector[p-1] beta;               // regression coefficients
  real mu;                            // overall effect
  real<lower=0> omega;                // betas prior hyperparameter
}

transformed parameters {
   vector[m] linpred;                 // linear predictor
   row_vector[p] betamu;              // coefficients+mu vector

   betamu = append_col(mu, beta);     // appending mu and betas
   linpred = X * betamu';             // linear predictor
}

model {
  omega ~ gamma(a_omega, b_omega);    //
  mu ~ normal(eta[1], 1/sqrt(tau));   //

   for(i in 1:(p-1)){                 //
      beta[i] ~ normal(eta[i+1], (1/sqrt(omega))*I[i,i]);
  }

  y ~ bernoulli_logit(linpred);       // fitting model
}

