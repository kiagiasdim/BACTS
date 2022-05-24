#' Logistic fit with STAN using Random intercepts (RE)
#'
#' @description This STAN function fits a logistic model to the patients data using random intercepts.
#' @export
#' @param m Integer ; number of in-silico patients (>=1)
#' @param p Integer ; number of parameters (including intercept) (>=2)
#' @param y Vector ; binary dependent variable (size 'm')
#' @param X Matrix ; model design matrix (dimension 'm x p')
#' @param tau Real ; intercept prior precision
#' @param eta Numeric vector ; intercept and betas prior mean (1st element of vector corresponds to intercept's mean)
#' @param I Matrix ; betas covariance matrix (dimension 'p x p')
#' @param a_omega Real ; omega prior shape hyperparameter
#' @param b_omega Real ; omega prior scale hyperparameter
#' @param npatID Integer ; number of unique in-silico patients profiles
#' @param patID Vector of integers ; in-silico patients ID labels. This vector consists of integers that represent the unique patients profiles included in data.
#' @param a_omega_u Real ; omega_u prior shape hyperparameter
#' @param b_omega_u Real ; omega_u prior rate hyperparameter
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains, etc ...).
#' @details The main scope of this model is to take into consideration the potential variability arising from UISS-TB when applied for patients with identical vector of features. Specifically, using UISS-TB numerous times for a specific patients profile will result in slighty differents results due to its stochastic nature.
#' Therefore, fitting this model shall be considered for in-silico data only. Running the model for either dataset (in-silico or in-vivo) corresponds to drawing inference from an independent source, without combining any information from the trials.
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @examples
#' # Loading simulated UISS-TB dataset
#' data("UISS_TB_data")
#'
#' # Cores for stan + rstan options
#' library(rstan)
#' # options(mc.cores = parallel::detectCores())
#' # rstan_options(auto_write = TRUE)
#'
#' # Fitting model on in-silico patients only - treatment group
#' fit_is_t <- IndSourceLogistRE_stan( m = nrow(X_is_t), p = ncol(X_is_t), y = y_is_t, X = X_is_t,
#'                                     tau = 2, eta = rep(0, ncol(X_is_t)), I = diag(ncol(X_is_t)),
#'                                     npatID = length(table(ID_is_t)), patID = ID_is_t,
#'                                     a_omega_u = a_omega_u_t, b_omega_u = b_omega_u_t,
#'                                     a_omega = 2, b_omega = .6, iter = 1e4, warmup = 1e3, chains = 3,
#'                                     pars = c("mu", "beta"), control = list(adapt_delta = .99) )
#'
#'
#' # Extracting output from fitting
#' out_is_t <- rstan::extract(fit_is_t)
#'
#' # Isolating overall effect 'mu' and 'betas'
#' PostS_is_t <- cbind( out_is_t$mu, out_is_t$beta )
#' colnames(PostS_is_t) <- c( "mu", paste("beta", seq(1:(ncol(X_is_t)-1)), sep = "_") )
#'
#' # Preview of posterior samples
#' head(PostS_is_t)
#'
#'
#'
#'
IndSourceLogistRE_stan <- function( m = m, p = p, y = y, X = X, tau = tau, eta = eta,
                                    I = I, a_omega = a_omega, b_omega = b_omega,
                                    npatID = npatID, patID = patID,
                                    a_omega_u = a_omega_u, b_omega_u = b_omega_u, ... )
{
  standata <- list( m = m, p = p, y = y, X = X, tau = tau, eta = eta, I = I, a_omega = a_omega, b_omega = b_omega,
                    npatID = npatID, patID = patID, a_omega_u = a_omega_u, b_omega_u = b_omega_u )
  out <- rstan::sampling( stanmodels$IndSourceLogistRE, data = standata, ... )
  return( out )
}





