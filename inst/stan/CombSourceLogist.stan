functions{
  /* Log Likelihood in-vivo data  */
  real ll_bern_logitPr(vector y, vector lp) {
    vector[num_elements(y)] pr;
    real sum_log_pr;
    real aux;
    for (i in 1:num_elements(y)){
      aux = exp(lp[i]) / (1+exp(lp[i]));
      pr[i] = y[i] * log(aux) + (1-y[i]) * log(1-aux);
    }
    sum_log_pr = sum(pr);
    return sum_log_pr;
  }
}

data {
  int<lower=0> p_v;           // number of parameters ; covariates + intercept in in-vivo data
  int<lower=0> p_b;           // number of common parameters ; common covariates from in-silico/vivo data
  int<lower=0> p_vo;          // number of parameters ; covariates  only in in-vivo data

  real<lower=0> tau;          // mu prior hyperparameter
  real<lower=0> a_omega;      // omega prior shape hyperparameter
  real<lower=0> b_omega;      // omega prior scale hyperparameter
  vector[p_vo] eta;           // prior betas_(v-c) mean ; only in-vivo
  matrix[p_vo, p_vo] I;       // prior betas_(v-c) covariance matrix structure ; benchmark identity divided by omega - only in-vivo
  row_vector[p_b] mu_s_b;      // prior betas_c mean Student's t ; from in-silico fitted model
  matrix[p_b, p_b] Sigma_s_b;  // prior betas_c covariance matrix Student's t ; from in-silico fitted model
  real<lower=1> nu_s_b;        // prior betas_c df Student's t ; from in-silico fitted model

  int<lower=0> m_v;           // number of in-vivo patients
  vector[m_v] y_v;            // dv in-vivo ; binary
  matrix[m_v, p_v] X_v;       // design matrix in-vivo
  real<lower=0> alpha;        // power prior parameter alpha
}

parameters {
  row_vector[p_v-1] beta_v;   // regression coefficients
  real mu;                    // overall effect
  real<lower=0> omega;        // scalar parameter of covariance matrix
}

transformed parameters {
   vector[m_v] lp_v;          // Linear predictor in-vivo data
   matrix[p_vo, p_vo] Omega;  // Covariance matrix for betas_(v-c) ; only in-vivo
   row_vector[p_v] betamu_v;  // coefficients+mu ; in-vivo
   row_vector[p_b] beta_b;    // coefficients ; common features
   row_vector[p_vo] beta_vo;  // coefficients ; only in in-vivo

   betamu_v = append_col(mu, beta_v);
   lp_v = X_v * betamu_v';

   beta_b = beta_v[1:p_b];                // Getting common betas
   beta_vo = beta_v[(p_b+1):(p_b+p_vo)];  // Getting only in-vivo betas

   Omega = 1/sqrt(omega) * I;
}

model{
  omega ~ gamma(a_omega, b_omega);
  mu ~ normal(0, 1/sqrt(tau));

  /* posterior distribution for betas in in-vivo data only */
  /* log-likelihood ; using in-vivo data */
  target += ll_bern_logitPr(y_v, lp_v);
  /* power prior ; using in-silico data */
  target += alpha * multi_student_t_lpdf(beta_b | nu_s_b, mu_s_b, Sigma_s_b);
  /* prior for betas in in-vivo data only */
  if( p_vo!=0 ) { target += multi_normal_lpdf(beta_vo | eta, Omega); }
}

