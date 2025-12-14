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
                   interv_delay = 1e30, # how many day after 1st case to start
                   interv_release = 1e30, # how many days after the last case to stop
                   interv_type = c("cull", "quarantine"),
                   rate_cull = 1e30, # immediate culling once response starts
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
  
  interv_type <- match.arg(interv_type)
  
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
  
  S <- E <- I <- C <- D <- R <- V <- status <- matrix(0, nrow = time, ncol = n_pop)
  S[1, ] <- ini_S
  E[1, ] <- ini_E
  I[1, ] <- ini_I
  C[1, ] <- ini_C
  D[1, ] <- ini_D
  R[1, ] <- ini_R
  V[1, ] <- ini_V
  
  has_cases <- rep(FALSE, n_pop)
  days_with_cases <- rep(0, n_pop)
  days_without_cases <- rep(0, n_pop)
  in_response <- rep(FALSE, n_pop)

  p_S_V <- vacc_coverage * vacc_efficacy
  rate_S_V <- -log(1 - p_S_V)
  
  
  for (t in seq_len(time - 1)) {
    
    ## state monitoring variables
    has_cases <- has_cases | (I[t, ] > 0)
    days_with_cases[has_cases] <- days_with_cases[has_cases] + 1
    days_without_cases[I[t, ] == 0] <- days_without_cases[I[t, ] == 0] + 1
    in_response <- (days_with_cases >= interv_delay) & (days_without_cases <= interv_release)
    status[t + 1, ] <- as.integer(in_response)
    
    ## individuals leaving S...
    ## - to E (new infections)
    ## - to V (vaccination)
    ## (culling will be handled separately)
    
    rate_S_E <- as.vector(delta %*% (beta * I[t, ]))
    if (interv_type == "cull") {
      rate_into_C <- c(0, rate_cull)[in_response + 1] # culling
    } else {
      rate_into_C <- 0
    }
    
    ## nested binomials are used to decide where individuals leaving S go
    rate_S_out <- rate_S_E + rate_S_V + rate_into_C
    p_S_out <- 1 - exp(-rate_S_out)
    n_S_out <- rbinom(n_pop, size = S[t, ], prob = p_S_out)
    
    S_E_ratio <- rate_S_E / rate_S_out
    S_E_ratio[!is.finite(S_E_ratio)] <- 0
    n_S_E <- rbinom(n_pop, size = n_S_out, prob = S_E_ratio)
    
    S_V_ratio <- rate_S_V / (rate_S_V + rate_into_C)
    S_V_ratio[!is.finite(S_V_ratio)] <- 0
    n_S_V <- rbinom(n_pop, size = n_S_out - n_S_E, prob = S_V_ratio)
    n_S_C <- n_S_out - n_S_E - n_S_V
    
    
    ## individuals leaving E
    rate_E_out <- sigma + rate_into_C
    p_E_out <- 1 - exp(-rate_E_out)
    n_E_out <- rbinom(n_pop, size = E[t, ], prob = p_E_out)
    E_I_ratio <- sigma / rate_E_out
    E_I_ratio[!is.finite(E_I_ratio)] <- 0
    n_E_I <- rbinom(n_pop, size = n_E_out, prob = E_I_ratio)
    n_E_C <- n_E_out - n_E_I
    
    ## individuals leaving I...
    ## - to recover (R)
    ## - to die from the disease (D)
    ## - to be culled (C)
    ##
    ## We first draw individuals going to R or D (i.e. not culled), then draw
    ## disease outcome from the CFR.
  
    rate_I_out <- gamma + rate_into_C
    p_I_out <- 1 - exp(-rate_I_out)
    n_I_out <- rbinom(n_pop, size = I[t, ], prob = p_I_out)
    I_RD_ratio <- gamma / rate_I_out
    I_RD_ratio[!is.finite(I_RD_ratio)] <- 0
    n_I_RD <- rbinom(n_pop, size = n_I_out, prob = I_RD_ratio)
    n_I_D <- rbinom(n_pop, size = n_I_RD, prob = cfr)
    n_I_R <- n_I_RD - n_I_D
    n_I_C <- n_I_out - n_I_RD
    
    
    ## Recovered individuals can be culled too
    n_R_C <- rbinom(n_pop, size = R[t, ], prob = 1 - exp(-rate_into_C))
    
    ## Vaccinated individuals can be culled too
    n_V_C <- rbinom(n_pop, size = V[t, ], prob = 1 - exp(-rate_into_C))
    
    ## all changes
    S[t + 1, ] <- S[t, ] - n_S_E - n_S_V - n_S_C
    E[t + 1, ] <- E[t, ] + n_S_E - n_E_I - n_E_C
    I[t + 1, ] <- I[t, ] + n_E_I - n_I_D - n_I_R - n_I_C
    C[t + 1, ] <- C[t, ] + n_S_C + n_E_C + n_I_C + n_R_C + n_V_C
    D[t + 1, ] <- D[t, ] + n_I_D
    R[t + 1, ] <- R[t, ] + n_I_R - n_R_C
    V[t + 1, ] <- V[t, ] + n_S_V - n_V_C
  }
  
  
  ## format output
  colnames(S) <- paste("S", seq_len(n_pop), sep = "_")
  colnames(E) <- paste("E", seq_len(n_pop), sep = "_")
  colnames(I) <- paste("I", seq_len(n_pop), sep = "_")
  colnames(C) <- paste("C", seq_len(n_pop), sep = "_")
  colnames(D) <- paste("D", seq_len(n_pop), sep = "_")
  colnames(R) <- paste("R", seq_len(n_pop), sep = "_")
  colnames(V) <- paste("V", seq_len(n_pop), sep = "_")
  colnames(status) <- paste("status", seq_len(n_pop), sep = "_")
  
  cbind.data.frame(S, E, I, C, D, R, V, status)
}

