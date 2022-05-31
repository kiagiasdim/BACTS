# **B**ayesian **A**ugmented **C**linical **T**rialS - **BACTS**

This library implements the methods described in [Kiagias et al (2021)](https://www.frontiersin.org/article/10.3389/fmedt.2021.719380) for combining information from RCT for TB therapeutic vaccination ([STriTuVaD](https://www.strituvad.eu)) and computer simualtions from [UISS-TB](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03762-5). 

To install run
``` R
devtools::install_github("kiagiasdim/BACTS")
```
You will need [Stan](https://mc-stan.org/users/interfaces/rstan) and the [fitHeavyTail](https://cran.r-project.org/web/packages/fitHeavyTail/) package to run `BACTS`

The library provides code to fit hierarchical Bayesian logistic models for dichotomous end points, with and without random effects.  We use the former to fit the *in silico* data and the latter for the *in vivo* trials.  

A second suite of routines combine the posterior distributions from these two sources of information into an in silico Augmented Clinial Trial.

## Combining *in silico* and *in vivo* experiments

The vignette below illustrates how to fit the logist models to *in silico* and *in vivo* sources and to combine them in a **B**ayesian **A**ugmented **C**linical **T**rial.

```R
# Loading simulated UISS-TB dataset
data("UISS_TB_data")

library(rstan)
library(fitHeavyTail)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Initial fitting of independent models for in vivo and in silico trials - treatment group.
## These models will be used to obtain a power prior for the common betas of the trials.
## We will then apply the combined model (augmented trial).

# Fitting logistic RE model on in-silico patients only - treatment group
fit_is_t <- IndSourceLogistRE_stan( m = nrow(X_is_t), p = ncol(X_is_t), y = y_is_t, X = X_is_t,
                                    tau = 2, eta = rep(0, ncol(X_is_t)), I = diag(ncol(X_is_t)),
                                    npatID = length(table(ID_is_t)), patID = ID_is_t,
                                    a_omega_u = a_omega_u_t, b_omega_u = b_omega_u_t,
                                    a_omega = 2, b_omega = .6, iter = 1e4, warmup = 1e3, chains = 3,
                                    pars = c("mu", "beta"), control = list(adapt_delta = .99) )

# Fitting logistic model on in-vivo patients only - treatment group
fit_iv_t <- IndSourceLogist_stan( m = nrow(X_iv_t), p = ncol(X_iv_t), y = y_iv_t, X = X_iv_t,
                                  tau = 2, eta = rep(0, ncol(X_iv_t)), I = diag(ncol(X_iv_t)),
                                  a_omega = 2, b_omega = .6, iter = 1e4, warmup = 1e3, chains = 3,
                                  pars = c("mu", "beta"), control = list(adapt_delta = .99) )

# Extracting output from fitting
out_is_t <- rstan::extract(fit_is_t)
out_iv_t <- rstan::extract(fit_iv_t)

# Isolating overall effect 'mu' and 'betas' from output - in vivo trial
PostS_iv_t <- cbind( out_iv_t$mu, out_iv_t$beta )

# Isolating overall effect 'mu' and 'betas' from output - in silico trial
# Note: We get only the common 'betas' with in vivo trial
PostS_is_t <- cbind( out_is_t$mu, out_is_t$beta[, 1:(ncol(PostS_iv_t)-1)] )

# Renaming columns
colnames(PostS_iv_t) <- colnames(PostS_is_t) <- c( "mu", paste("beta", seq(1:(ncol(X_iv_t)-1)), sep = "_") )

# Obtaining power prior 'alpha' using 'alpha_update' function
alpha_t <- alpha_update( PostS_iv_mu = PostS_iv_t[, "mu"], PostS_is_mu = PostS_is_t[, "mu"], M = ncol(X_is_t) )

# Getting power prior for common betas for in silico trial
# using a multivariate t-student distribution
betasPrior_b_t <- fit_mvt( out_is_t$beta[, 1:(ncol(PostS_iv_t)-1)], nu = "MLE-diag")

# Fitting combined model for in vivo patients (augmented trial)
fit_comb_t <- CombSourceLogist_stan( m_v = nrow(X_iv_t), y_v = y_iv_t, X_v = X_iv_t,
                                     p_v = ncol(X_iv_t), p_b = ncol(X_iv_t) - 1, p_vo = 0,
                                     mu_s_b = unname(betasPrior_b_t$mu),
                                     Sigma_s_b = unname(betasPrior_b_t$cov), nu_s_b = betasPrior_b_t$nu,
                                     alpha = alpha_t, iter = 1e4, warmup = 1e3, chains = 3,
                                     pars = c("mu", "beta_v"), control = list(adapt_delta = .99) )

# Extracting output from fitting
out_comb_t <- rstan::extract(fit_comb_t)

# Isolating overall effect 'mu' and 'betas' from output - augmented trial
PostS_comb_t <- cbind( out_comb_t$mu, out_comb_t$beta_v )

# Preview results
head(PostS_comb_t)
```
