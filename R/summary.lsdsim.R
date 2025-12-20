#' Summary for lsdsim objects
#' 
#' This function sums compartments across all populations in a [lsdsim] object, 
#' for each time step.
#' 
#' @export
#' @param object an [lsdsim] object

summary.lsdsim <- function(object) {
  x <- object
  n_steps <- nrow(x)
  comps <- c("S", "E", "I", "D", "C", "R", "V")
  comp_sums <- list()
  for(e in comps) {
    to_keep <- grep(sprintf("^[%s]", e), colnames(x))
    comp_sums[[e]] <- rowSums(x[, to_keep])
  }
  out <- cbind(time = seq_len(n_steps), data.frame(comp_sums))
  out
}