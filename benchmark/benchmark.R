
# Some small benchmarking on a standard desktop...
#
# Expectations
# ------------
#
# The algorithm should scale linearly with time and badly (^2) with the number 
# of patches, due to the delta matrix. Future implementations might want to try
# and use sparse matrices.
# 
# Interpreter calls can be condensed and there might be some time gain to be 
# made by compiling the R code. We leave this to the user.
# 
# 
# What we do
# ----------
# 
# 1. compare compiled version to normal implementation
# 2. evaluate scaling with time
# 3. evaluate scaling of the algorithm for increasing number of populations
# 4. compare "no intervention" vs "full interventions"

devtools::load_all()

# 1. compiled version vs normal implementation
lsdsim_compiled <- compiler::cmpfun(lsdsim)
compiled_vec <- rep(c(FALSE, TRUE), each = 20)
vec_fun <- c(lsdsim, lsdsim_compiled)[compiled_vec + 1]
exp_1_time <- lapply(
  vec_fun, 
  function(f) system.time(
    f(grid_size = 10, time = 365)
  )
)
exp_1_res <- cbind.data.frame(
  compiled = compiled_vec, 
  Reduce(rbind, exp_1_time)
)
dotchart(
  exp_1_res$elapse, pch = 20, 
  main = "lsddim: scaling with/without compilation", 
  xlab = "Runtime (s)",
  ylab = "Simulation (red = compiled code)",
  col = compiled_vec + 1
)


# 2. scaling with time
time_vec <- seq(50, to = 1000, by = 50)
exp_2_time <- lapply(
  time_vec, 
  function(t) system.time(lsdsim_compiled(4, time = t))
)
exp_2_res <- cbind.data.frame(
  time = time_vec, 
  Reduce(rbind, exp_2_time)
)
plot(
  exp_2_res$time, exp_2_res$elapse, pch = 20, 
  main = "lsddim: scaling with time", 
  xlab = "Duration of the simulation (days)", 
  ylab = "Runtime (s)"
)


# 3. scaling with grid size (nb of populations ^2)
grid_size <- c(1:5, seq(10, 60, by = 10))
exp_3_time <- lapply(
  grid_size, 
  function(k) system.time(lsdsim_compiled(grid_size = k, time = 365))
)
exp_3_res <- cbind.data.frame(
  grid_size = grid_size, 
  Reduce(rbind, exp_3_time)
)
plot(
  exp_3_res$grid_size, exp_3_res$elapse, pch = 20, 
  main = "lsddim: scaling with grid size", 
  xlab = "Grid size (sqrt(nb. populations))", 
  ylab = "Runtime (s)"
)


# 4. scaling with/without intervention
interv_vec <- rep(c(FALSE, TRUE), each = 20)
exp_4_time <- lapply(
  interv_vec, 
  function(e) system.time(
    lsdsim_compiled(grid_size = 10, time = 300, 
           vaccination = e, 
           culling = e, 
           quarantine = e, 
           insecticide = e, 
           interv_delay = 1
    )
  )
)
exp_4_res <- cbind.data.frame(
  interv = interv_vec, 
  Reduce(rbind, exp_4_time)
)
dotchart(
  exp_4_res$elapse, pch = 20, 
  main = "lsddim: scaling with intervention", 
  xlab = "Runtime (s)",
  ylab = "Simulation (red = intervention)",
  col = interv_vec + 1
)


