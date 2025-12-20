#' Summary for lsdsim objects
#' 
#' This function sums compartments across all populations in a [lsdsim] object, 
#' for each time step. It also reports the number of patches in 'naive', 
#' 'response' and 'ring' state at every time step.
#' 
#' @export
#' @param object an [lsdsim] object

summary.lsdsim <- function(object) {
  x <- object
  n_steps <- nrow(x)
  
  ## sum compartments over all populations
  comps <- c("S", "E", "I", "D", "C", "R", "V", "N")
  comp_sums <- list()
  for(e in comps) {
    to_keep <- grep(sprintf("^[%s]", e), colnames(x))
    comp_sums[[e]] <- rowSums(x[, to_keep])
  }
  
  ## total cases: including E, I, D, R
  ## AR is total cases / population size
  total_cases <- comp_sums[["E"]] + 
    comp_sums[["I"]] + 
    comp_sums[["D"]] + 
    comp_sums[["R"]]
  ar <- total_cases / comp_sums[["N"]]
  ar[comp_sums[["N"]] == 0] <- 0
  
  ## number of populations with at least one infectious individual (I):
  ## - n_pop_infec: currently
  ## - n_pop_ever_infec: ever infected since the beginning of the outbreak
  I <- x[, grep("^I", colnames(x))]
  n_pop_infec <- rowSums(I > 0)
  n_pop_ever_infec <- sum(colSums(I) > 0)
  
  ## get the breakdown of status 
  status <- x[, grep("status", names(x))]
  n_naive <- rowSums(status == "naive")
  n_response <- rowSums(status == "response")
  n_ring <- rowSums(status == "ring")
  
  ## pull results together
  out <- cbind(
    time = seq_len(n_steps), 
    data.frame(comp_sums),
    total_cases = total_cases,
    ar = ar,
    n_pop_infec = n_pop_infec, 
    n_pop_ever_infec = n_pop_ever_infec,
    n_naive = n_naive,
    n_response = n_response,
    n_ring = n_ring
  )
  out$total_naive <- sum(n_naive)
  out$total_response <- sum(n_response)
  out$total_ring <- sum(n_ring)
  out
}