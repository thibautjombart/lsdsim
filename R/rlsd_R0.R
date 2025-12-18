#' Draw LSD basic reproduction numbers
#'
#' This random number generator draws average values of the basic reproduction
#' number (R0) for LSD using literature-driven parameters and log-normal
#' distributions. R0 is defined as the average number of secondary infections
#' per infected cattle in a fully susceptible herd. Note that the value returned
#' and associated variability refer to a population average, and do not reflect
#' individual-level variation. The values used here refer to transmission via
#' vector Stomoxys calcitrans (see Table 4 in source for other vectors), based
#' on virus isolation in cattle blood. 
#'
#' @source Gubbins S. Using the basic reproduction number to assess the risk of
#'   transmission of lumpy skin disease virus by biting insects. Transbound
#'   Emerg Dis. 2019;66: 1873–1883.
#'
#' @export
#' @param n the number of values to draw
#' @param mu the population average, defaults to 9.7
#' @param sd the standard deviation of the lognormal distribution, defaults to
#'   0.1
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
rlsd_R0 <- function(n, mu = 9.7, sd = 0.05) {
  rlnorm(n, meanlog = log(mu), sdlog = sd)
}
