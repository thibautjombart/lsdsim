#' Draw simulation parameters for LSD
#'
#' This function creates tables of input parameters for LSD transmission using
#' estimates from Gubbins et al. (2019), based on virus isolation in cattle
#' blood, assuming transmission via Stomoxys calcitrans. It uses sub-functions
#' [rlsd_latent], [rlsd_infec_period], [rlsd_R0] and [rlsd_cfr]. See these
#' functions to change default values. Note that the value returned and
#' associated variability refer to a population average, and do not reflect
#' individual-level variation.
#'
#' @source Gubbins S. Using the basic reproduction number to assess the risk of
#'   transmission of lumpy skin disease virus by biting insects. Transbound
#'   Emerg Dis. 2019;66: 1873–1883.
#'
#'   Manić M, Stojiljković M, Petrović M, Nišavić J, Bacić D, Petrović T, et al.
#'   Epizootic features and control measures for lumpy skin disease in
#'   south-east Serbia in 2016. Transbound Emerg Dis. 2019;66: 2087–2099.
#'
#'   Abera Z, Degefu H, Gari G, Dibessa ZA. Review on epidemiology and economic
#'   importance of lumpy skin disease. Int J Basic Appl Virol. 2015. Available:
#'   https://www.academia.edu/download/76656747/2.pdf
#'
#'   Şevik M, Doğan M. Epidemiological and molecular studies on lumpy skin
#'   disease outbreaks in Turkey during 2014-2015. Transbound Emerg Dis.
#'   2017;64: 1268–1279.
#'
#' @returns A data.frame with one set of parameters per row, including: R0, the
#'   latent period, the infectious period, the transmission rate (beta), rate of
#'   incubation (sigma), and the clearance rate (gamma).
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
    R0_I = rlsd_R0_I(n),
    R0_A = rlsd_R0_A(n),
    latent = rlsd_latent(n),
    infec_period = rlsd_infec_period(n),
    cfr = rlsd_cfr(n), 
    pasymp = rlsd_pasymp(n) 
  )

  out$beta_I <- out$R0_I / out$infec_period
  out$beta_A <- out$R0_A / out$infec_period
  out$sigma <- 1 / out$latent
  out$gamma <- 1 / out$infec_period
  out
}
