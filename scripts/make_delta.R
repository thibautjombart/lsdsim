#' Build a delta matrix from geographic coordinates
#' 
#' Here coordinates are provided as x and y in a 2-columns matrix
#' 
#' @param coords a 2-columns matrix of xy coordinates
#' @diffusion the value of the diffusion to be used, defined as the fraction of 
#' the force of infection going outside each patch towards the neighbours
make_delta <- function(coords, diffusion = 0) {
  if (ncol(coords) > 2) {
    coords <- coords[, c("x", "y")]
  }
  out <- prop.table(1 * (as.matrix(dist(coords)) < 1.01), 1) * diffusion 
  diag(out) <- 1 - diffusion
  out
}
