#' Alternative beta parametrisation
#' 
#' This internal function generates shape parameters for beta distributions
#' using a specified mean and standard deviation
#' 
#' @noRd

beta_musd_to_shapes <- function(mu, sd) {
  list(
    shape1 = mu / (sd * sd), 
    shape2 = (1 - mu) / (sd * sd)  
  )
}
