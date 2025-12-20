#' Change connectivity matrix in quarantine
#' 
#' This function changes the delta matrix used in [lsdsim] to implement 
#' quanrantine, defined as a reduction in outward-going transmission. The 
#' strategy is to produce a complete new delta matrix where all populations are
#' in quarantine. Only specific rows of this matrix will then be used during the 
#' simulation as needed.
#'
#' @export
#' @param delta the connectivity matrix used in [lsdsim]
#' @param efficacy_in the relative reduction in transmission within herds 
#' @param efficacy_out the relative reduction in outward-going transmission 
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' 
add_quarantine <- function(delta, efficacy_in, efficacy_out) {
  ## strategy:
  ## - we reduce all terms by a relative fraction equal to efficacy_out
  ## - we restore the diagonal terms of the original matrix, and reduce them by 
  ##   a relative fraction equal to efficacy_in
  ## 
  ## note that the resulting matrix is no longer row-standardized, but this is
  ## fine, we do want quarantine to reduce transmission, not merely impact where
  ## new cases go
  out <- delta
  out <- out * (1 - efficacy_out)
  diag(out) <- diag(delta) * (1 - efficacy_in)
  out
}
