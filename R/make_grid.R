#' Make a 2-D grid
#' 
#' Builds a square, regular grid of pre-specified size
#' 
#' @export
#' @param size the size of the square grid
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}

make_grid <- function(size) {
  xy <- expand.grid(seq_len(size), seq_len(size))
  xy <- as.data.frame(xy)
  names(xy) <- c("x", "y")
  xy
}
