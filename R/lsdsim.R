#' Function to simulate an LSD outbreak in a metapopulation
#'
#' This function simulates an outbreak of Lumpy skin disease in a meta-population
#' arranged on a regular grid of size `grid_size`, and a diffusion process 
#' defined by a `delta` matrix.
#'
lsdsim <- function(grid_size = 1,
                   time = 365,
                   beta = 0.1, # density-dependent infection rate
                   sigma = 1/7, # inv. latent period
                   gamma = 1/20, # inv. duration of infection
                   cfr = 0.1, # case fatality ratio
                   vacc_coverage = 0,
                   vacc_efficacy = 0.65,
                   delta = NULL, 
                   diffusion = 0,
                   ini_S = 0,
                   ini_E = 0,
                   ini_I = 0,
                   ini_C = 0,
                   ini_D = 0,
                   ini_R = 0,
                   ini_V = 0
                   ) {
  
  ## handle arguments
  ### recycle arguments as needed
  ### define defaults for NULL values
  n_pop <- grid_size^2
  ini_S <- rep(ini_S, length.out = n_pop)
  ini_E <- rep(ini_E, length.out = n_pop)
  ini_I <- rep(ini_I, length.out = n_pop)
  ini_C <- rep(ini_C, length.out = n_pop)
  ini_D <- rep(ini_D, length.out = n_pop)
  ini_R <- rep(ini_R, length.out = n_pop)
  ini_V <- rep(ini_V, length.out = n_pop)
  
  
  ## geographic structure
  xy <- make_grid(grid_size)
  delta <- make_delta(xy, diffusion)
  
  S <- E <- I <- C <- D <- R <- V <- matrix(0, nrow = time, ncol = n_pop)
  S[1, ] <- ini_S
  E[1, ] <- ini_E
  I[1, ] <- ini_I
  C[1, ] <- ini_C
  D[1, ] <- ini_D
  R[1, ] <- ini_R
  V[1, ] <- ini_V
  
  has_cases <- rep(FALSE, n_pop)
  in_response <- rep(FALSE, n_pop)
  days_with_cases <- rep(0, n_pop)
  
  p_S_V <- vacc_coverage * vacc_efficacy
  rate_S_V <- -log(1 - p_S_V)
  
  p_E_I <- 1 - exp(-sigma)
  
  p_I_out <- 1 - exp(-gamma)
  
  for (t in seq_len(time - 1)) {
    
    ## individuals leaving S...
    ## - to E (new infections)
    ## - to V (vaccination)
    ## (culling will be handled separately)
    
    rate_S_E <- as.vector(delta %*% (beta * I[t, ]))
    rate_S_out <- rate_S_E + rate_S_V
    p_S_out <- 1 - exp(-rate_S_out)
    n_S_out <- rbinom(n_pop, size = S[t, ], prob = p_S_out)
    rate_ratio <- rate_S_E / rate_S_out
    rate_ratio[!is.finite(rate_ratio)] <- 0
    n_S_E <- rbinom(n_pop, size = n_S_out, prob = rate_ratio)
    n_S_V <- n_S_out - n_S_E
    
    
    ## individuals leaving E
    n_E_I <- rbinom(n_pop, size = E[t, ], prob = p_E_I)
    
    ## individuals leaving I...
    ## - to recover (R)
    ## - to die (D)
  
    n_I_out <- rbinom(n_pop, size = I[t, ], prob = p_I_out)
    n_I_D <- rbinom(n_pop, size = n_I_out, prob = cfr)
    n_I_R <- n_I_out - n_I_D
    
    ## all changes
    S[t + 1, ] <- S[t, ] - n_S_E - n_S_V
    E[t + 1, ] <- E[t, ] + n_S_E - n_E_I
    I[t + 1, ] <- I[t, ] + n_E_I - n_I_D - n_I_R
    D[t + 1, ] <- D[t, ] + n_I_D
    R[t + 1, ] <- R[t, ] + n_I_R
    V[t + 1, ] <- V[t, ] + n_S_V
  }
  
  
  ## format output
  colnames(S) <- paste("S", seq_len(n_pop), sep = "_")
  colnames(E) <- paste("E", seq_len(n_pop), sep = "_")
  colnames(I) <- paste("I", seq_len(n_pop), sep = "_")
  colnames(C) <- paste("C", seq_len(n_pop), sep = "_")
  colnames(D) <- paste("D", seq_len(n_pop), sep = "_")
  colnames(R) <- paste("R", seq_len(n_pop), sep = "_")
  colnames(V) <- paste("V", seq_len(n_pop), sep = "_")
  
  cbind.data.frame(S, E, I, C, D, R, V)
}

