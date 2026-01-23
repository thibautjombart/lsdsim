#' Draw LSD basic reproduction numbers
#'
#' This random number generator draws average values of the basic reproduction
#' number (R0) for LSD using literature-driven parameters and log-normal
#' distributions. `rlsd_R0_I()` is for symptomatic individuals, while
#' `rlsd_R0_A` is for asymptomatic cases. R0 is defined as the average number of
#' secondary infections per infected cattle in a fully susceptible herd. Note
#' that the value returned and associated variability refer to a population
#' average, and do not reflect individual-level variation. The values used here
#' refer to transmission via vector Stomoxys calcitrans (see Table 4 in source
#' for other vectors), based on virus isolation in cattle blood.
#'
#' @source Gubbins S. Using the basic reproduction number to assess the risk of
#'   transmission of lumpy skin disease virus by biting insects. Transbound
#'   Emerg Dis. 2019;66: 1873–1883.
#'
#' @export
#' @param n the number of values to draw
#' @param mu the population average, defaults to 9.7 for symptomatic individuals,
#' and 1.2 for asymptomatic cases
#' @param sd the standard deviation of the lognormal distribution, defaults to
#'   0.1 for symptomatic individuals, and 0.1 for asymptomatic cases
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @rdname rlsd_R0
#' @examples
#' 
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(2, 1))
#' hist(rlsd_R0_I(1e5), main = "R0 for symptomatic cases", xlab = "R0")
#' hist(rlsd_R0_A(1e5), main = "R0 for asymptomatic cases", xlab = "R0")
#' par(opar)
#' 
rlsd_R0_I <- function(n, mu = 9.7, sd = 0.05) {
  rlnorm(n, meanlog = log(mu), sdlog = sd)
}

#' @rdname rlsd_R0
#' @export
rlsd_R0_A <- function(n, mu = 1.2, sd = 0.1) {
  rlnorm(n, meanlog = log(mu), sdlog = sd)
}
