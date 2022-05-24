\name{UISS_TB_data}
\alias{UISS_TB_data}
\docType{data}
\title{
Simulated UISS-TB dataset
}
\description{
This dataset has been created using the Universal Immune System Simulator for Tuberculosis vaccination (UISS-TB) to represent both outcomes from the in-silico and in-vivo trials.
}
\usage{data("UISS_TB_data")}
\format{
The dataset consists of individual matrices and vectors corresponding to the input values for implementing the individual and combined models for both control and treatment groups.
  \describe{
    \item{\code{X_is_c}}{design matrix from in-silico trial - control group}
    \item{\code{X_is_t}}{design matrix from in-silico trial - treatment group}
    \item{\code{X_iv_c}}{design matrix from in-vivo trial - control group}
    \item{\code{X_iv_t}}{design matrix from in-vivo trial - treatment group}
    \item{\code{y_is_c}}{dependent variable (0 or 1) from in-silico trial - control group}
    \item{\code{y_is_t}}{dependent variable (0 or 1) from in-silico trial - treatment group}
    \item{\code{y_iv_c}}{dependent variable (0 or 1) from in-vivo trial - control group}
    \item{\code{y_iv_t}}{dependent variable (0 or 1) from in-vivo trial - treatment group}
    \item{\code{ID_is_c}}{vector of integers represting the unique profiles of patients from in-silico trial - control group. Here every patient has a different initial vector of characteristics.}
    \item{\code{ID_is_t}}{vector of integers represting the unique profiles of patients from in-silico trial - treatment group. Here every patient has a different initial vector of characteristics.}
    \item{\code{a_omega_u_t}}{numeric value ; fitted estimate for hyperparameter of random intercepts prior shape - treatment group }
    \item{\code{b_omega_u_t}}{numeric value ; fitted estimate for hyperparameter of random intercepts prior rate - treatment group }
    \item{\code{a_omega_u_c}}{numeric value ; fitted estimate for hyperparameter of random intercepts prior shape - control group }
    \item{\code{b_omega_u_c}}{numeric value ; fitted estimate for hyperparameter of random intercepts prior rate - control group }
  }
}
\details{
...
}
\source{
UISS-TB (github link here)
}
\references{
UISS-TB paper
}


