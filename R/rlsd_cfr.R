#' Draw LSD case fatality ratio
#'
#' This random number generator draws average case fatality ratio (CFR) for LSD
#' using literature-driven parameters and beta distributions. The CFR is defined
#' as the proportion of infected cattles who die from the disease. Reported
#' values refer to virus isolation in cattle blood. Defaults correpond to a CFR
#' of 3%. See [rbeta()] for details on means and variances.
#'
#' @source Manić M, Stojiljković M, Petrović M, Nišavić J, Bacić D, Petrović T,
#'   et al. Epizootic features and control measures for lumpy skin disease in
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
#' @export
#' @param n the number of values to draw
#' @param mu the mean of the distribution; defaults to 0.03
#' @param sd the standard deviation of the distribution; defaults to 0.05
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @examples
#'
#' hist(rlsd_cfr(1e5), main = "CFR of LSD", xlab = "Proportion of deaths")
#' 
rlsd_cfr <- function(n, mu = 0.03, sd = 0.05) {
  params <- beta_musd_to_shapes(mu, sd)
  rbeta(n, params$shape1, params$shape2)
}

