#' Draw simulation parameters for LSD
#'
#' This function creates tables of input parameters for LSD transmission using
#' estimates from Gubbins et al. (2019), based on virus isolation in cattle
#' blood, assuming transmission via Stomoxys calcitrans. It uses sub-functions
#' [rlsd_latent], [rlsd_infec_period], and [rlsd_R0]. See these functions to
#' change default values. Note that the value returned and associated
#' variability refer to a population average, and do not reflect
#' individual-level variation.
#'
#' @source Gubbins S. Using the basic reproduction number to assess the risk of
#'   transmission of lumpy skin disease virus by biting insects. Transbound
#'   Emerg Dis. 2019;66: 1873–1883.
#'   
#' @returns A data.frame with one set of parameters per row, including: R0, the 
#' latent period, the infectious period, the transmission rate (beta), rate of 
#' incubation (sigma), and the clearance rate (gamma).
#'
#' @export
#' @param n the number of values to draw
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' 
#' @examples 
#' x <- rlsd_param(1e3)
#' head(x)
#' summary(x)
#' hist(x$R0, main = "Distribution of R0", xlab = "R0", border = "white")
#' hist(x$latent, main = "Latent period", 
#'      xlab = "Latent period (days)", 
#'      border = "white")
#' hist(x$infec_period, main = "Infectious period", 
#'      xlab = "Infectious period (days)", 
#'      border = "white")
#' 
rlsd_param <- function(n) {
  out <- data.frame(
    R0 = rlsd_R0(n),
    latent = rlsd_latent(n),
    infec_period = rlsd_infec_period(n)
  )
  
  out$beta <- out$R0 / out$infec_period
  out$sigma <- 1 / out$latent
  out$gamma <- 1 / out$infec_period
  out
}
