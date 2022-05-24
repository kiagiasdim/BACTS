\name{FrontiersData}
\alias{FrontiersData}
\docType{data}
\title{
Simulated UISS-TB dataset used in Frontiers paper
}
\description{
This dataset has been created using the Universal Immune System Simulator for Tuberculosis vaccination (UISS-TB) to represent both outcomes from the in-silico and in- vivo trials.
}
\usage{data("FrontiersData")}
\format{
The dataset consists of individual matrices and vectors corresponding to the input values for implementing the individual and combined models for both control and treatment groups.
  \describe{
    \item{\code{IS_c}}{list containing quantities for in-silico trial - control group}
    \item{\code{IS_t}}{list containing quantities for in-silico trial - treatment group}
    \item{\code{IV_c}}{list containing quantities for in-vivo trial - control group}
    \item{\code{IV_t}}{list containing quantities for in-vivo trial - treatment group}
    \item{\code{AG_c}}{list containing quantities for augmented trial - control group}
    \item{\code{AG_t}}{list containing quantities for augmented trial - treatment group}
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


