#' Logistic fit with STAN
#'
#' @description This STAN function fits a logistic model to the patients data.
#' @export
#' @param m Integer ; number of patients (>=1)
#' @param p Integer ; number of parameters (including intercept) (>=2)
#' @param y Vector ; binary dependent variable (size 'm')
#' @param X Matrix ; model design matrix (dimension 'm x p')
#' @param tau Real ; intercept prior precision
#' @param eta Numeric vector ; intercept and betas prior mean (1st element of vector corresponds to intercept's mean)
#' @param I Matrix ; betas covariance matrix (dimension 'p x p')
#' @param a_omega Real ; omega prior shape hyperparameter
#' @param b_omega Real ; omega prior scale hyperparameter
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains, etc ...).
#' @details Fitting this logistic model can be applied for either in-silico or in-vivo data, from either the control or treatment group. Running the model for either dataset (in-silico or in-vivo) corresponds to drawing inference from an independent source, without combining any information from the trials.
#' @return An object of class `stanfit` returned by `rstan::sampling`.
#' @examples
#' # Loading simulated UISS-TB dataset
#' data("UISS_TB_data")
#'
#' # Cores for stan + rstan options
#' library(rstan)
#' # options(mc.cores = parallel::detectCores())
#' # rstan_options(auto_write = TRUE)
#'
#' # Fitting model on in-vivo patients only - control group
#' fit_iv_c <- IndSourceLogist_stan( m = nrow(X_iv_c), p = ncol(X_iv_c), y = y_iv_c, X = X_iv_c,
#'                                   tau = 2, eta = rep(0, ncol(X_iv_c)), I = diag(ncol(X_iv_c)),
#'                                   a_omega = 2, b_omega = .6, iter = 1e4, warmup = 1e3, chains = 3,
#'                                   pars = c("mu", "beta"), control = list(adapt_delta = .99) )
#'
#' # Extracting output from fitting
#' out_iv_c <- rstan::extract(fit_iv_c)
#'
#' # Isolating overall effect 'mu' and 'betas'
#' PostS_iv_c <- cbind( out_iv_c$mu, out_iv_c$beta )
#' colnames(PostS_iv_c) <- c( "mu", paste("beta", seq(1:(ncol(X_iv_c)-1)), sep = "_") )
#'
#' # Preview of posterior samples
#' head(PostS_iv_c)
#'
#'
#'
#'
IndSourceLogist_stan <- function( m = m, p = p, y = y, X = X, tau = tau, eta = eta,
                                  I = I, a_omega = a_omega, b_omega = b_omega, ... )
{
  standata <- list( m = m, p = p, y = y, X = X, tau = tau, eta = eta, I = I, a_omega = a_omega, b_omega = b_omega )
  out <- rstan::sampling( stanmodels$IndSourceLogist, data = standata, ... )
  return( out )
}




