#' Augmented trial fit with STAN - combing information from in-silico and in-vivo data
#'
#' @description This STAN function fits a logistic model to the patients data using random intercepts.
#' @export
#' @param m_v Integer ; number of in-vivo patients (>=1)
#' @param y_v Vector ; binary dependent variable (size 'm_v')
#' @param X_v Matrix ; model design matrix (dimension 'm_v x p_v')
#' @param p_v Integer ; number of in-vivo parameters (including intercept) (>=2)
#' @param p_b Integer ; number of common parameters from in-vivo and in-silico trials (without intercept) (>=2)
#' @param p_vo Integer ; number of in-vivo only parameters (without intercept) (>=2)
#' @param tau Real ; intercept prior precision
#' @param eta Numeric vector ; only in-vivo betas prior mean (size 'p_vo')
#' @param I Matrix ; betas covariance matrix (dimension of 'p x p')
#' @param a_omega Real ; omega prior shape hyperparameter (only in-vivo data)
#' @param b_omega Real ; omega prior scale hyperparameter (only in-vivo data)
#' @param mu_s_b Numeric vector ; prior - t-student - mean (size 'p_b') for common betas (from in-vivo and in-silico data)
#' @param Sigma_s_b Matrix ; prior - t-student - covariance matrix (dimension of 'p_b x p_b') for common betas (from in-vivo and in-silico data)
#' @param nu_s_b Real ; rior - t-student - degrees of freedom (>0)
#' @param alpha Real ; power of power-prior, indicating the strength of similarity from in-vivo and in-silico data
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains, etc ...).
#' @details The main scope of this model is to take into consideration the potential variability arising from UISS-TB when applied for patients with identical vector of features. Specifically, using UISS-TB numerous times for a specific patients profile will result in slighty differents results due to its stochastic nature.
#' Therefore, fitting this model shall be considered for in-silico data only. Running the model for either dataset (in-silico or in-vivo) corresponds to drawing inference from an independent source, without combining any information from the trials.
#' @return An object of class `stanfit` returned by `rstan::sampling`.
#' @examples
#' # Loading simulated UISS-TB dataset
#' data("UISS_TB_data")
#'
#' # Cores for stan + rstan options
#' library(rstan)
#' library(fitHeavyTail)
#' # options(mc.cores = parallel::detectCores())
#' # rstan_options(auto_write = TRUE)
#'
#' ## Initial fitting of independent models for in vivo and in silico trials - treatment group.
#' ## These models will be used to obtain a power prior for the common betas of the trials.
#' ## We will then apply the combined model (augmented trial).
#'
#' # Fitting logistic RE model on in-silico patients only - treatment group
#' fit_is_t <- IndSourceLogistRE_stan( m = nrow(X_is_t), p = ncol(X_is_t), y = y_is_t, X = X_is_t,
#'                                     tau = 2, eta = rep(0, ncol(X_is_t)), I = diag(ncol(X_is_t)),
#'                                     npatID = length(table(ID_is_t)), patID = ID_is_t,
#'                                     a_omega_u = a_omega_u_t, b_omega_u = b_omega_u_t,
#'                                     a_omega = 2, b_omega = .6, iter = 1e4, warmup = 1e3, chains = 3,
#'                                     pars = c("mu", "beta"), control = list(adapt_delta = .99) )
#'
#' # Fitting logistic model on in-vivo patients only - treatment group
#' fit_iv_t <- IndSourceLogist_stan( m = nrow(X_iv_t), p = ncol(X_iv_t), y = y_iv_t, X = X_iv_t,
#'                                   tau = 2, eta = rep(0, ncol(X_iv_t)), I = diag(ncol(X_iv_t)),
#'                                   a_omega = 2, b_omega = .6, iter = 1e4, warmup = 1e3, chains = 3,
#'                                   pars = c("mu", "beta"), control = list(adapt_delta = .99) )
#'
#' # Extracting output from fitting
#' out_is_t <- rstan::extract(fit_is_t)
#' out_iv_t <- rstan::extract(fit_iv_t)
#'
#' # Isolating overall effect 'mu' and 'betas' from output - in vivo trial
#' PostS_iv_t <- cbind( out_iv_t$mu, out_iv_t$beta )
#'
#' # Isolating overall effect 'mu' and 'betas' from output - in silico trial
#' # Note: We get only the common 'betas' with in vivo trial
#' PostS_is_t <- cbind( out_is_t$mu, out_is_t$beta[, 1:(ncol(PostS_iv_t)-1)] )
#'
#' # Renaming columns
#' colnames(PostS_iv_t) <- colnames(PostS_is_t) <- c( "mu", paste("beta", seq(1:(ncol(X_iv_t)-1)), sep = "_") )
#'
#' # Obtaining power prior 'alpha' using 'alpha_update' function
#' alpha_t <- alpha_update( PostS_iv_mu = PostS_iv_t[, "mu"], PostS_is_mu = PostS_is_t[, "mu"], M = ncol(X_is_t) )
#'
#' # Getting power prior for common betas for in silico trial
#' # using a multivariate t-student distribution
#' betasPrior_b_t <- fit_mvt( out_is_t$beta[, 1:(ncol(PostS_iv_t)-1)], nu = "MLE-diag")
#'
#' # Fitting combined model for in vivo patients (augmented trial)
#' fit_comb_t <- CombSourceLogist_stan( m_v = nrow(X_iv_t), y_v = y_iv_t, X_v = X_iv_t,
#'                                      p_v = ncol(X_iv_t), p_b = ncol(X_iv_t) - 1, p_vo = 0,
#'                                      mu_s_b = unname(betasPrior_b_t$mu),
#'                                      Sigma_s_b = unname(betasPrior_b_t$cov), nu_s_b = betasPrior_b_t$nu,
#'                                      alpha = alpha_t, iter = 1e4, warmup = 1e3, chains = 3,
#'                                      pars = c("mu", "beta_v"), control = list(adapt_delta = .99) )
#'
#' # Extracting output from fitting
#' out_comb_t <- rstan::extract(fit_comb_t)
#'
#' # Isolating overall effect 'mu' and 'betas' from output - augmented trial
#' PostS_comb_t <- cbind( out_comb_t$mu, out_comb_t$beta_v )
#'
#' # Preview results
#' head(PostS_comb_t)
#'
CombSourceLogist_stan <- function( m_v = m_v, y_v = y_v, X_v = X_v,
                                   p_v = p_v, p_b = p_b, p_vo = p_vo,
                                   tau = 2, eta = rep(0, p_vo), I = diag(p_vo), a_omega = 2, b_omega = .6,
                                   mu_s_b = mu_s_b, Sigma_s_b = Sigma_s_b, nu_s_b = nu_s_b,
                                   alpha = alpha, ... )
{
  standata <- list( m_v = m_v, y_v = y_v, X_v = X_v, p_v = p_v, p_b = p_b, p_vo = p_vo,
                    tau = tau, eta = eta, I = I, a_omega = a_omega, b_omega = b_omega,
                    mu_s_b = mu_s_b, Sigma_s_b = Sigma_s_b, nu_s_b = nu_s_b, alpha = alpha )
  out <- rstan::sampling( stanmodels$CombSourceLogist, data = standata, ... )
  return( out )
}


