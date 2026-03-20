#' Draw proportion of asymptomatic for LSD 
#'
#' This random number generator draws the proportion of asymptomatic cases for
#' LSD using literature-driven parameters and beta distributions. See [rbeta()]
#' for details on means and variances. . 
#'
#' @source TBC
#'
#' @export
#' @param n the number of values to draw
#' @param mu the mean of the distribution; defaults to 0.5
#' @param sd the standard deviation of the distribution; defaults to 0.05
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @examples
#' 
#' hist(
#'   rlsd_pasymp(1e5), 
#'   main = "Proportion of asymptomatic cases", 
#'   xlab = "P (asymptomatic)"
#' )
#' 
rlsd_pasymp <- function(n, mu = 0.5, sd = 0.05) {
  params <- beta_musd_to_shapes(mu, sd)
  rbeta(n, params$shape1, params$shape2)
}
