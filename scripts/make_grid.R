#' Make grid
#' 
#' Builds a regular grid of pre-specified size
#' 
make_grid <- function(size) {
  xy <- expand.grid(seq_len(size), seq_len(size))
  xy <- as.data.frame(xy)
  names(xy) <- c("x", "y")
  xy
}
