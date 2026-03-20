#' Function to simulate an LSD outbreak in a metapopulation
#'
#' This function simulates an outbreak of Lumpy skin disease in a meta-population
#' arranged on a regular grid of size `grid_size`, and a diffusion process 
#' defined by a `delta` matrix.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @export
#'
#' @param grid_size The size of the regular grid to be used. The number of
#'   populations will be this number squared. Defaults to 1.
#'
#' @param time An `integer` indicating the duration of the simulation, in
#'   days. Defaults to 365.
#'
#' @param beta_I The rate of infection for symptomatic cases (I); defaults to
#'   0.1.
#'
#' @param beta_A The rate of infection for asymptomatic cases (A); defaults to
#'   0.01.
#'
#' @param sigma The inverse of the latency period; defaults to 1/7.
#'
#' @param gamma The inverse of the duration of illness; defaults to 1/20.
#'
#' @param cfr The case fatality ratio; defaults to 0.1 (i.e. 10%).
#'
#' @param pasymp The average proportion of asymptomatic cases; defaults to 0.5
#'   (i.e. 50%).
#'
#' @param mass_culling A `logical` indicating if mass culling should be used in
#'   the response. Only the farm where the response takes place is
#'   affected. Defaults to FALSE.
#' 
#' @param select_culling A `logical` indicating if selective culling should be
#'   used in the response. Only the farm where the response takes place is
#'   affected. Defaults to FALSE.
#' 
#' @param vaccination A `logical` indicating if vaccination should be used in
#'   the response. The farm where the response takes place and its neighbours
#'   are affected. Defaults to FALSE.
#' 
#' @param quarantine A `logical` indicating if quarantine should be used in
#'   the response. The farm where the response takes place and its neighbours
#'   are affected. Defaults to FALSE.
#' 
#' @param insecticide A `logical` indicating if insecticide should be used in
#'   the response. The farm where the response takes place and its neighbours
#'   are affected. Defaults to FALSE.
#'
#' @param rate_cull The rate at which culling takes place once the response has
#'   been triggered. Defaults to 1e30, i.e. all animals are killed in 1 day.
#'
#' @param vacc_coverage The vaccine coverage, i.e. the proportion of individuals
#'   getting vaccinated once the response has been triggered. Defaults to 0 - no
#'   vaccination.
#'
#' @param vacc_efficacy The vaccine efficacy, i.e. the proportion of vaccinated
#'   individuals who effectively become immune. Defaults to 0.65 (i.e. 65%).
#'
#' @param quarant_efficacy_in The relative reduction of transmission inside the
#'   farm due to quarantine. Defaults to 0.2 (i.e. 20%).
#'
#' @param quarant_efficacy_out The relative reduction of transmission towards
#'   neighbouring farms due to quarantine. Defaults to 0.9 (i.e. 90%).
#'
#' @param insect_efficacy The relative reduction of transmission due to
#'   insecticide spray. Defaults to  0.5 (i.e. 50%).
#'
#' @param interv_delay The delay to intervention, defined as the number of days
#'   after the first symptomatic case for the intervention to start. Defaults to
#'   1e30, i.e. no intervention.
#'
#' @param interv_release The number of days after the last symptomatic case
#'   after which the intervention stops. Defaults to 28 days.
#'
#' @param delta The connectivity matrix defining the proportions of the forces
#'   of infection going from populations (in rows), to populations (in
#'   columns). Defaults to `NULL`, in which case the matrix is generated
#'   automatically using rook connectivity and the `diffusion` argument.
#'
#' @param diffusion The proportion of the force of infection directed towards
#'   neighbouring populations, used the build the `delta` matrix if it is set to
#'   `NULL`. Defaults to 0, i.e. no connectivity (`delta` is the identity
#'   matrix).
#'
#' @param ini_S The initial number of susceptible individuals, recycled as
#'   needed to match the total number of populations. Defaults to 0.
#' 
#' @param ini_E The initial number of exposed individuals, recycled as needed to
#'   match the total number of populations. Defaults to 0.
#' 
#' @param ini_A The initial number of asymptomatic cases, recycled as needed to
#'   match the total number of populations. Defaults to 0.
#' 
#' @param ini_I The initial number of symptomatic cases, recycled as needed to
#'   match the total number of populations. Defaults to 0.
#' 
#' @param ini_C The initial number of culled animals, recycled as needed to
#'   match the total number of populations. Defaults to 0.
#' 
#' @param ini_D The initial number of individuals who died from the disease,
#'   recycled as needed to match the total number of populations. Defaults to 0.
#' 
#' @param ini_R The initial number of recovered individuals, recycled as needed
#'   to match the total number of populations. Defaults to 0.
#' 
#' @param ini_V The initial number of individuals immunized through vaccination,
#'   recycled as needed to match the total number of populations. Defaults to 0.
#' 


lsdsim <- function(
    grid_size = 1,
    time = 365,
    beta_I = 0.1, # infection rate for I
    beta_A = 0.01, # infection rate for A
    sigma = 1/7, # inv. latent period
    gamma = 1/20, # inv. duration of infection
    cfr = 0.1, # case fatality ratio
    pasymp = 0.5, # proportion of asymptomatic cases
    mass_culling = FALSE, # mass culling in affected populations?
    select_culling = FALSE, # culling of infected (I) individuals only
    vaccination = FALSE, # vaccinate affected population and neighbours?
    quarantine = FALSE, # quarantine in affected pop and neighbours?
    insecticide = FALSE, # use insecticide in affected pop and neighbours?
    rate_cull = 1e30, # immediate culling once response starts
    vacc_coverage = 0, # prop of individuals getting vaccinated
    vacc_efficacy = 0.65, # prop of vaccinated individuals getting protection
    quarant_efficacy_in = 0.2, # 20% transmission reduction within herd
    quarant_efficacy_out = 0.9, # 90% outward transmission reduction
    insect_efficacy = 0.5, # % reduction of transmission due to insecticide
    interv_delay = 1e30, # how many day after 1st case to start
    interv_release = 28, # how many days after the last case to stop
    delta = NULL, 
    diffusion = 0,
    ini_S = 0,
    ini_E = 0,
    ini_A = 0,
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
  ini_A <- rep(ini_A, length.out = n_pop)
  ini_I <- rep(ini_I, length.out = n_pop)
  ini_C <- rep(ini_C, length.out = n_pop)
  ini_D <- rep(ini_D, length.out = n_pop)
  ini_R <- rep(ini_R, length.out = n_pop)
  ini_V <- rep(ini_V, length.out = n_pop)
  ini_N <- ini_S + ini_E + ini_A + ini_I + ini_R + ini_V
  
  ## geographic structure
  xy <- make_grid(grid_size)
  if (is.null(delta)) {
    delta <- make_delta(xy, diffusion)
  } else {
    stopifnot(nrow(delta) == ncol(delta))
    stopifnot(nrow(delta) == n_pop)
  }
  
  ## note: new_I will be used to track the incidence of new cases
  S <- E <- A <- I <- new_I <- C <- D <- R <- V <- N <- 
    matrix(0, nrow = time, ncol = n_pop)
  S[1, ] <- ini_S
  E[1, ] <- ini_E
  A[1, ] <- ini_A
  I[1, ] <- new_I[1, ] <- ini_I
  C[1, ] <- ini_C
  D[1, ] <- ini_D
  R[1, ] <- ini_R
  V[1, ] <- ini_V
  N[1, ] <- ini_N
  status <- matrix("naive", nrow = time, ncol = n_pop)
  
  has_an_outbreak <- rep(FALSE, n_pop)
  days_in_outbreak <- rep(0, n_pop)
  days_without_cases <- rep(0, n_pop)
  in_response <- rep(FALSE, n_pop)

  rate_vacc <- -log(1 - (vacc_coverage * vacc_efficacy))
  
  delta_quarant <- add_quarantine(
    delta, 
    quarant_efficacy_in,
    quarant_efficacy_out
  )
  
  delta_ori <- delta # copy of the original delta matrix
  temp <- delta
  diag(temp) <- 0
  list_neighbours <- apply(temp > 0, 1, which) # list of neighbours for each pop
  
  for (t in seq_len(time - 1)) {
    
    ## state monitoring variables
    ##
    ## - in_response: logical indicating populations with ongoing outbreak
    ## response
    ## - was_in_response: same, for the previous time steps
    ## - id_pop_in_response: indices of populations in response
    ## - id_pop_in_ring: indices of neighbours to populations in response
    ## - id_pop_response_and_ring: indices of pop in response and ring
    ##
    ## The process for monitoring outbreaks and response states is as follows:
    ##
    ## - has_an_outbreak: indicates pops which has had at least one case; once 
    ##   TRUE, it stays TRUE until the end of an outbreak response; at any time
    ##   step, it includes all previously TRUE, and pops with at least one I; 
    ##   resets to FALSE once an outbreak response ends
    ## - days_in_outbreak: counts the number of days since the first case; this 
    ##   is used to know when to trigger a new response; it is reset to zero 
    ##   once an outbreak response ends
    ## - days_without_cases: counts the number of days with no case in the pop;
    ##   it is reset to 0 once there is a least one case
    
    
    ## monitor pops with an ongoing outbreak
    has_an_outbreak <-  has_an_outbreak | (I[t, ] > 0)
    days_in_outbreak[has_an_outbreak] <- days_in_outbreak[has_an_outbreak] + 1
    
    ## monitor pops with a period of no cases
    days_without_cases[I[t, ] == 0] <- days_without_cases[I[t, ] == 0] + 1
    days_without_cases[I[t, ] > 0] <- 0
    
    ## trigger new responses
    was_in_response <- in_response
    in_response <- (days_in_outbreak >= interv_delay) & (days_without_cases <= interv_release)
    
    ## stop responses after x days with no case, and reset counters
    response_stopped <- was_in_response & !in_response
    has_an_outbreak[response_stopped] <- FALSE
    days_in_outbreak[response_stopped] <- 0
    
    ## monitor IDs of pops in response, and the associated rings
    id_pop_in_response <- which(in_response)
    id_pop_in_ring <-  unique(unlist(list_neighbours[in_response]))
    id_pop_response_and_ring <- unique(c(id_pop_in_response, id_pop_in_ring))
    status[t + 1, id_pop_in_ring] <- "ring"
    status[t + 1, id_pop_in_response] <- "response"
    
    ## response-associated processes
    ##
    ## - mass_culling: mass culling of affected population 
    ## - select_culling: culling of I only once in response mode
    ## - vaccination: vaccination of affected pop and neighbours
    ## - quarantine: quarantine of affected pop and neighbours
    ## - insecticide: insecticide administered in affected pop and neighbours
    if (mass_culling) {
      rate_into_C <- rep(0, n_pop)
      rate_into_C[in_response] <- rate_cull
    } else {
      rate_into_C <- 0
    }
    
    if (select_culling) {
      rate_I_C <- rep(0, n_pop)
      rate_I_C[in_response] <- rate_cull
    } else {
      rate_I_C <- 0
    }
    
    if (vaccination) {
      rate_S_V <- rep(0, n_pop)
      rate_S_V[id_pop_response_and_ring] <- rate_vacc
    } else {
      rate_S_V <- 0
    }
    
    if (quarantine) {
      delta <- delta_ori
      delta[id_pop_response_and_ring, ] <- delta_quarant[id_pop_response_and_ring, ]
    }
   
    if (insecticide) {
       current_beta_I <- rep(beta_I, n_pop)
       current_beta_I[id_pop_response_and_ring] <- 
         current_beta_I[id_pop_response_and_ring] * (1 - insect_efficacy)
       
       current_beta_A <- rep(beta_A, n_pop)
       current_beta_A[id_pop_response_and_ring] <- 
         current_beta_A[id_pop_response_and_ring] * (1 - insect_efficacy)
    } else {
      current_beta_I <- beta_I
      current_beta_A <- beta_A
    }
    
    
    ## individuals leaving S...
    ## - to be infected (E)
    ## - to be vaccinated (V)
    ## - to be culled (C)
    ##
    ## nested binomials are used to decide where individuals leaving S go
    rate_S_E <- as.vector(
      delta %*% ((current_beta_I * I[t, ] + current_beta_A * A[t, ]) / N[t, ])
    )
    
    rate_S_E[is.na(rate_S_E)] <- 0 # for populations where N = 0
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
    
    
    ## individuals leaving E...
    ## - to become asymptomatic (A)
    ## - to become infectious, symptomatic (I)
    ## - to be culled (C)
    ##
    ## We first draw the number of individuals leaving E, then decide which ones
    ## go to A or I, then draw the ones going to A, I, and deduce C from the 
    ## rest.
    rate_E_out <- sigma + rate_into_C
    p_E_out <- 1 - exp(-rate_E_out)
    n_E_out <- rbinom(n_pop, size = E[t, ], prob = p_E_out)
    E_AI_ratio <- sigma / rate_E_out
    E_AI_ratio[!is.finite(E_AI_ratio)] <- 0
    n_E_AI <- rbinom(n_pop, size = n_E_out, prob = E_AI_ratio)
    n_E_A <- rbinom(n_pop, size = n_E_AI, prob = pasymp)
    n_E_I <- n_E_AI - n_E_A
    n_E_C <- n_E_out - n_E_AI
   
    
    ## individuals leaving A...
    ## - to recover (R)
    ## - to be culled (C), either though mass culling of selective culling
    ##
    
    rate_A_out <- gamma + rate_into_C
    p_A_out <- 1 - exp(-rate_A_out)
    n_A_out <- rbinom(n_pop, size = A[t, ], prob = p_A_out)
    A_R_ratio <- gamma / rate_A_out
    A_R_ratio[!is.finite(A_R_ratio)] <- 0
    n_A_R <- rbinom(n_pop, size = n_A_out, prob = A_R_ratio)
    n_A_C <- n_A_out - n_A_R
    
     
    ## individuals leaving I...
    ## - to recover (R)
    ## - to die from the disease (D)
    ## - to be culled (C), either though mass culling of selective culling
    ##
    ## We first draw individuals going to R or D (i.e. not culled), then draw
    ## disease outcome from the CFR.
  
    rate_I_out <- gamma + rate_into_C + rate_I_C
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
    E[t + 1, ] <- E[t, ] + n_S_E - n_E_A - n_E_I - n_E_C
    A[t + 1, ] <- A[t, ] + n_E_A - n_A_R - n_A_C
    I[t + 1, ] <- I[t, ] + n_E_I - n_I_D - n_I_R - n_I_C
    C[t + 1, ] <- C[t, ] + n_S_C + n_E_C + n_A_C + n_I_C + n_R_C + n_V_C
    D[t + 1, ] <- D[t, ] + n_I_D
    R[t + 1, ] <- R[t, ] + n_A_R + n_I_R - n_R_C
    V[t + 1, ] <- V[t, ] + n_S_V - n_V_C
    N[t + 1, ] <- S[t + 1, ] + E[t + 1, ] + A[t + 1, ] + I[t + 1, ] + R[t + 1, ] + V[t + 1, ]
    new_I[t + 1, ] <- n_E_I
  }
  
  
  ## format output
  colnames(S) <- paste("S", seq_len(n_pop), sep = "_")
  colnames(E) <- paste("E", seq_len(n_pop), sep = "_")
  colnames(A) <- paste("A", seq_len(n_pop), sep = "_")
  colnames(I) <- paste("I", seq_len(n_pop), sep = "_")
  colnames(C) <- paste("C", seq_len(n_pop), sep = "_")
  colnames(D) <- paste("D", seq_len(n_pop), sep = "_")
  colnames(R) <- paste("R", seq_len(n_pop), sep = "_")
  colnames(V) <- paste("V", seq_len(n_pop), sep = "_")
  colnames(N) <- paste("N", seq_len(n_pop), sep = "_")
  colnames(new_I) <- paste("new_I", seq_len(n_pop), sep = "_")
  colnames(status) <- paste("status", seq_len(n_pop), sep = "_")
  
  out <- cbind.data.frame(S, E, A, I, C, D, R, V, N, new_I, status)
  class(out) <- c("lsdsim", class(out))
  out
}

