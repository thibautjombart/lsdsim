#' Change connectivity matrix in quarantine
#' 
#' This function changes the delta matrix used in [lsdsim] to implement 
#' quanrantine, defined as a reduction in outward-going transmission. The 
#' strategy is to produce a complete new delta matrix where all populations are
#' in quarantine. Only specific rows of this matrix will then be used during the 
#' simulation as needed.
#'
#' @param delta the connectivity matrix used in [lsdsim]
#' @param efficacy the relative reduction in outward-going transmission 
#' 
add_quarantine <- function(delta, efficacy) {
  ## strategy: 
  ## - we isolate non-diagonal terms
  ## - we reduce them by a relative fraction equal to the efficacy
  ## - we find new diagonal terms to ensure row-standardisation
  out <- delta
  diag(out) <- 0
  out <- out * (1 - efficacy)
  new_diag <- 1 - rowSums(out)
  diag(out) <- new_diag
  out
}
