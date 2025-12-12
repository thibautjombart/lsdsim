#   #############
#  ###############
# #################
# # SEIRDV  MODEL #
# ################
#  ##############
#   ############


#######################################################
# Core equations for transitions between compartments #
#######################################################

## Note the use of [i] syntax for indexing populations
##
## Compartments include:
## - S: susceptible
## - E: exposed (not infectious yet)
## - I: infectious
## - D: dead
## - R: recovered
## - V: vaccinated
##
## Transmission is frequency-dependent.
##
## Intervention monitoring:
## - 'status' has n_populations values being
##     + 0: naive (no cases so far)
##     + 1: in response mode
## - intervention kicks in 'interv_delays' days after the first I
## - there are 2 types of intervention:
##     + interv_type = 1: mass culling of infected herds
##     + interv_type = 2: quarantine reducing diffusion by a factor of 
## `quarant_efficacy`

# Core equations for transitions between compartments

## Note the use of [i] syntax for indexing populations

update(S[]) <- S[i] - n_S_E[i] - n_S_V[i] - n_S_D[i]
update(E[]) <- E[i] + n_S_E[i] - n_E_I[i] - n_E_D[i]
update(I[]) <- I[i] + n_E_I[i] - n_I_R[i] - n_I_D[i]
update(R[]) <- R[i] + n_I_R[i] - n_R_D[i]
update(D[]) <- D[i] + n_S_D[i] + n_E_D[i] + n_I_D[i] + n_R_D[i] + n_V_D[i]
update(V[]) <- V[i] + n_S_V[i] - n_V_D[i]

## N is the living population, i.e. people able to have contacts with S
N[] <- S[i] + E[i] + I[i] + R[i] + V[i]

output(status[]) <- status[i]



##################
# Initial states #
##################

## compartments
initial(S[]) <- ini_S[i]
initial(E[]) <- ini_E[i]
initial(I[]) <- ini_I[i]
initial(R[]) <- ini_R[i]
initial(D[]) <- ini_D[i]
initial(V[]) <- ini_V[i]

## active: 0 = no outbreak going on ; 1 = outbreak going on
active[] <- if (as.integer(step) == 1) 0 else active[i]

## days_active: number of days since 1st case of the outbreak
days_active[] <- if (as.integer(step) == 1) 0 else days_active[i]

## days_no_case: number of days with 0 case
days_no_case[] <- if (as.integer(step) == 1) 0 else days_no_case[i]

## status: 0 = naive population (no response); 1 = response
status[] <- if (as.integer(step) == 1) 0 else status[i]


#####################################
# Update on intervention monitoring #
#####################################
## counting days of activity (i.e. since 1st case appeared)
active[] <- if (I[i] > 0) 1 else active[i]
days_active[] <- if (active[i] == 1) days_active[i] + 1 else 0

## promotion to response mode
status[] <- if (days_active[i] > interv_delay) 1 else status[i]

## counting days with 0 case
days_no_case[] <- if (I[i] == 0) days_no_case[i] + 1 else 0

## stop response where cases has stopped for long enough
status[] <- if (days_no_case[i] > interv_release) 0 else status[i]
active[] <- if (days_no_case[i] > interv_release) 0 else active[i]


###########################
# Parameters declarations #
###########################

ini_S[] <- user(0)
ini_E[] <- user(0)
ini_I[] <- user(0)
ini_R[] <- user(0)
ini_D[] <- user(0)
ini_V[] <- user(0)

beta <- user(0)
sigma <- user(0)
gamma <- user(0)
mu <- user(0)
vacc_coverage <- user(0)
vacc_efficacy <- user(0)
interv_delay <- user(0)
quarant_efficacy <- user(0)
interv_type <- user(1)
interv_release <- user(0)

delta[,] <- user(0) # population x population matrix



###################
# Item dimensions #
###################

n_populations <- user(1)

dim(ini_S) <- n_populations
dim(ini_E) <- n_populations
dim(ini_I) <- n_populations
dim(ini_R) <- n_populations
dim(ini_D) <- n_populations
dim(ini_V) <- n_populations

dim(S) <- n_populations
dim(E) <- n_populations
dim(I) <- n_populations
dim(R) <- n_populations
dim(D) <- n_populations
dim(V) <- n_populations
dim(N) <- n_populations

dim(rate_S_V) <- n_populations
dim(p_S_V) <- n_populations
dim(rate_S_out) <- n_populations
dim(delta) <- c(n_populations, n_populations)
dim(lambda_prod) <- c(n_populations, n_populations)
dim(lambda) <- n_populations
dim(p_S_out) <- n_populations
dim(n_S_out) <- n_populations
dim(n_S_E) <- n_populations
dim(n_S_V) <- n_populations
dim(n_S_D) <- n_populations


dim(n_E_I) <- n_populations
dim(n_E_D) <- n_populations

dim(n_I_out) <- n_populations
dim(n_I_R) <- n_populations
dim(n_I_D) <- n_populations

dim(n_R_D) <- n_populations
dim(n_V_D) <- n_populations

dim(active) <- n_populations
dim(status) <- n_populations
dim(days_active) <- n_populations
dim(days_no_case) <- n_populations

dim(beta_vec) <- n_populations

####################################################
# Calculations of the stochastic changes of states #
####################################################

##############################
# Transitions: S->E and S->V #
##############################
#
# The code below calculates individual rates of S->E
#
# 'delta' is a population x population dispersal matrix, where rows are sources
# and columns recipients of the FOI; each row sums to 1
#
# Individuals can leave S to go to either E or V. The rates are summed to
# calculate how many people leave S, and nested binomials are used to decide
# where they go.

## coverage depends on response status

p_S_V[] <- vacc_coverage * vacc_efficacy # prop of vaccinated; could exceed 1
p_S_V[] <- if (p_S_V[i] > 1) 1 else p_S_V[i] # avoid over-shooting again, should be redundant
rate_S_V[] <- -log(1 - p_S_V[i]) # corresponding rate

## infectivity depends on response status
#beta_vec[] <- if (status[i] == 1) beta * (1 - interv_efficacy) else beta
#lambda_prod[, ] <- beta_vec[i] * delta[i, j] * I[i] / N[i]
lambda_prod[, ] <- beta * delta[i, j] * I[i] / N[i]
lambda[] <- sum(lambda_prod[, i]) # sum by column
rate_S_out[] <- lambda[i] + rate_S_V[i] # rate at which individuals leave S

## draw out of S to E and V
p_S_out[] <- 1 - exp(-rate_S_out[i]) # indiv proba of leaving S
n_S_out[] <- if (S[i] > 0) rbinom(S[i], p_S_out[i]) else 0
n_S_E[] <- if (n_S_out[i] > 0) rbinom(n_S_out[i], lambda[i] / rate_S_out[i]) else 0
n_S_V[] <- n_S_out[i] - n_S_E[i]



#####################
# Transitions: E->I #
#####################

n_E_I[] <- if (E[i] > 0) rbinom(E[i], 1 - exp(-sigma)) else 0



##############################
# Transitions: I->R and I->D #
##############################

# Here we first compute how many individuals leave I (n_I_out), and then decide
# where they go using nested binomial distributions.
#
# 'gamma' is the rate at which individuals leave I
# 'mu' is the CFR

n_I_out[] <- if (I[i] > 0) rbinom(I[i], 1 - exp(-gamma)) else 0
n_I_R[] <- if (n_I_out[i] > 0) rbinom(n_I_out[i], 1 - mu) else 0
n_I_D[] <- n_I_out[i] - n_I_R[i]



