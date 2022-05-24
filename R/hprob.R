#' Weight function based on the similarity probability of two samples
#'
#' @description This function provides a weight with respect to the probability 'p' describing the similarity of two (posterior) samples. Its main feature is the symmetry around .5 and it is based on two parameters 'lamda' and 'k', controlling the scale and shape respectively, determining the penalty level as 'p' deviates for .5.
#' @export
#' @param k scalar (positive) ; parameter controlling the shape of the function.
#' @param lambda scalar (positive - between 0 and 1) ; parameter controlling the scale of the function.
#' @param p scalar (between 0 and 1) ; probability to get weight for - here this is related on how similar two posterior samples appear to be.
#' @details The choice of parameters 'lambda' and 'k' denote the scale and shape of the function respectively. Specifically, larger values of 'lambda' and/or 'k' tend to produce a faster decrease from the peak of the function, whereas, smaller values correspond to a smoother and slower decrease from the peak.
#' @examples
#' hprob(k = 1, lambda = .1, p = .1)
#'
hprob <- function(k = 1, lambda = .4, p){

  if ( lambda >=1 ) { stop("Invalid parameter 'lambda'") }
  if ( p < .5 & p >= 0 ) { f <- 1 - exp( -( p/lambda )^k )
  } else if ( p >= .5 & p <= 1 ) { f <- 1 - exp( -( (1-p)/lambda )^k )
  } else { stop("Invalid prob 'p'") }
  return( f/(1 - exp( -(.5/lambda)^k) ) )

}
