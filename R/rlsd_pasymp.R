#' Draw proportion of asymptomatic for LSD 
#'
#' This random number generator draws the proportion of asymptomatic cases for
#' LSD using literature-driven parameters and beta distributions. See [rbeta()]
#' for details on means and variances. . 
#'
#' @source 
#'
#' @export
#' @param n the number of values to draw
#' @param shape1 the first shape parameter of the beta distribution; defaults to
#' 100
#' @param shape2 the second shape parameter of the beta distribution; defaults to
#' 100
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @examples
#' 
#' hist(
#'   rlsd_pasymp(1e5), 
#'   main = "Proportion of asymptomatic cases", 
#'   xlab = "P (asymptomatic)"
#' )
#' 
rlsd_pasymp <- function(n, shape1 = 100, shape2 = 100) {
  rbeta(n, shape1, shape2)
}
