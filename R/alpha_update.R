#' Update of the power prior parameter 'alpha'
#'
#' @description This function obtains the power prior parameter 'alpha' that describes the agreement level of the posterior samples from the in-silico and in-vivo trials, using the weight function 'hprob' to penalise their similarity probability. The parameter 'alpha' is expressed as the ratio of the effective size of the in-silico trial and the size of the virtual cohort, with the latter being a function of 'hprob' and the maximum number of in-silico patients allowed for the augmented trial.
#' @export
#' @param PostS_iv_mu vector ; posterior sample of the overall effect 'mu' from the in-vivo trial.
#' @param PostS_is_mu vector ; posterior sample of the overall effect 'mu' from the in-silico trial.
#' @param M integer ; size of the in-silico patients cohort.
#' @param Mmax integer ; maximum number of in-silico patients allowed in the augmented trial.
#' @examples
#' # Loading simulated UISS-TB dataset
#' data("UISS_TB_data")
#'
#' # Cores for stan + rstan options
#' library(rstan)
#' # options(mc.cores = parallel::detectCores())
#' # rstan_options(auto_write = TRUE)
#'
#' ## Initial fitting of independent models for in vivo and in silico trials - control group.
#' ## These models will be used to obtain a power prior for the common betas of the trials.
#' ## We will then apply the combined model (augmented trial).
#'
#' # Fitting logistic RE model on in-silico patients only - control group
#' fit_is_c <- IndSourceLogistRE_stan( m = nrow(X_is_c), p = ncol(X_is_c), y = y_is_c, X = X_is_c,
#'                                     tau = 2, eta = rep(0, ncol(X_is_c)), I = diag(ncol(X_is_c)),
#'                                     npatID = length(table(ID_is_c)), patID = ID_is_c,                                                                 a_omega_u = a_omega_u_c, b_omega_u = b_omega_u_c,
#'                                     a_omega = 2, b_omega = .6, iter = 1e4, warmup = 1e3, chains = 3,
#'                                     pars = c("mu", "beta"), control = list(adapt_delta = .99) )
#'
#' # Fitting logistic model on in-vivo patients only - control group
#' fit_iv_c <- IndSourceLogist_stan( m = nrow(X_iv_c), p = ncol(X_iv_c), y = y_iv_c, X = X_iv_c,
#'                                   tau = 2, eta = rep(0, ncol(X_iv_c)), I = diag(ncol(X_iv_c)),
#'                                   a_omega = 2, b_omega = .6, iter = 1e4, warmup = 1e3, chains = 3,
#'                                   pars = c("mu", "beta"), control = list(adapt_delta = .99) )
#'
#' # Extracting output from fitting
#' out_is_c <- rstan::extract(fit_is_c)
#' out_iv_c <- rstan::extract(fit_iv_c)
#'
#' # Isolating overall effect 'mu' and 'betas' from output - in vivo trial
#' PostS_iv_c <- cbind( out_iv_c$mu, out_iv_c$beta )
#'
#' # Isolating overall effect 'mu' and 'betas' from output - in silico trial
#' # Note: We get only the common 'betas' with in vivo trial
#' PostS_is_c <- cbind( out_is_c$mu, out_is_c$beta[, 1:(ncol(PostS_iv_c)-1)] )
#'
#' # Renaming columns
#' colnames(PostS_iv_c) <- colnames(PostS_is_c) <- c( "mu", paste("beta", seq(1:(ncol(X_iv_c)-1)), sep = "_") )
#'
#' # Obtaining power prior 'alpha' using 'alpha_update' function
#' alpha_c <- alpha_update( PostS_iv_mu = PostS_iv_c[, "mu"], PostS_is_mu = PostS_is_c[, "mu"], M = ncol(X_is_c) )
#' alpha_c
#'
alpha_update <- function( PostS_iv_mu, PostS_is_mu, M, Mmax = M ){

  # Calculating probability on endpoints of in-silico and in-vivo trials
  Pr_phi.iv_LESS_phi.is <- sum( PostS_iv_mu <= PostS_is_mu ) / length(PostS_is_mu)

  # Getting effective size using penatly function
  m <- hprob( p = Pr_phi.iv_LESS_phi.is ) * Mmax

  # Derived 'alpha'
  alpha <- m / M

  return(alpha)

}

