#' Build a delta matrix from geographic coordinates
#'
#' Spatial coordinates are provided as x and y in a 2-columns matrix.
#' Connectivity is based upon rook-relationships, whereby a location is
#' connected to locations immediately left, right, above, and under it. In
#' addition, cliques can be added to the matrix, i.e. highly-connected locations
#' (see details).
#'
#' @details Cliques are chosen at random from the existing locations, using a
#'   Binomial draw with probability `p_clique`. New connections are added so
#'   that cliques are connected to an average of `clique_connect` of all other
#'   locations.
#'
#' @export
#' @param coords a 2-columns matrix of xy coordinates
#' @param diffusion the value of the diffusion to be used, defined as the
#'   fraction of the force of infection going outside each patch towards the
#'   neighbours
#' @param p_clique the proportion of cliques (highly connected locations);
#'   defaults to zero (no cliques)
#' @param clique_connect the average connectivity of cliques; defaults to 0
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}

make_delta <- function(coords, 
                       diffusion = 0,
                       p_clique = 0,
                       clique_connect = 0
                       ){
  
  ## define neighours from rook connectivity
  out <- 1 * as.matrix(dist(coords))^2 < 1.01 # identify neighbours
  
  ## add cliques using hypergeometric distribution
  n_cliques <- stats::rbinom(1, nrow(out), p_clique)
  id_cliques <- sample.int(nrow(out), n_cliques)
  
  for (i in id_cliques) {
    to_replace <- !out[i, ]
    new_connections <- sample(
      c(TRUE, FALSE), 
      size = sum(to_replace), 
      prob = c(clique_connect, 1 - clique_connect), 
      replace = TRUE
    )
    out[i, to_replace] <- new_connections
  }
  
  ## zero the diagonal, standardise the matrix
  diag(out) <- 0
  out <- prop.table(out, 1) * diffusion
  diag(out) <- 1 - diffusion
  dimnames(out) <- NULL
  out
}
