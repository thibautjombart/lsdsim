#' Build a delta matrix from geographic coordinates
#' 
#' Here coordinates are provided as x and y in a 2-columns matrix
#' 
#' @export
#' @param coords a 2-columns matrix of xy coordinates
#' @param diffusion the value of the diffusion to be used, defined as the
#'   fraction of the force of infection going outside each patch towards the
#'   neighbours
#'   
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}

make_delta <- function(coords, diffusion = 0) {
  out <- 1 * as.matrix(dist(coords))^2 < 1.01 # identify neighbours
  diag(out) <- 0
  out <- prop.table(out, 1) * diffusion
  diag(out) <- 1 - diffusion
  dimnames(out) <- NULL
  out
}
