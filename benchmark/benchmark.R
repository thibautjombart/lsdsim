
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
# 1. evaluate scaling with time
# 2. evaluate scaling of the algorithm for increasing number of populations
# 3. compare "no intervention" vs "full interventions"

devtools::load_all()


# 1. scaling with time
time_vec <- seq(50, to = 1000, by = 50)
exp_1_time <- lapply(
  time_vec, 
  function(t) system.time(lsdsim(4, time = t))
)
exp_1_res <- cbind.data.frame(time = time_vec, Reduce(rbind, exp_1_time))
plot(
  exp_1_res$time, exp_1_res$elapse, pch = 20, 
  main = "lsddim: scaling with time", 
  xlab = "Duration of the simulation (days)", 
  ylab = "Runtime (s)"
)


# 2. scaling with grid size (nb of populations ^2)
grid_size <- 1:30
exp_2_time <- lapply(
  grid_size, 
  function(k) system.time(lsdsim(grid_size = k, time = 20))
)
exp_2_res <- cbind.data.frame(grid_size = grid_size, Reduce(rbind, exp_2_time))
plot(
  exp_2_res$grid_size, exp_2_res$elapse, pch = 20, 
  main = "lsddim: scaling with grid size", 
  xlab = "Grid size (sqrt(nb. populations))", 
  ylab = "Runtime (s)"
)


# 3. scaling with/without intervention
interv_vec <- rep(c(FALSE, TRUE), each = 20)
exp_3_time <- lapply(
  interv_vec, 
  function(e) system.time(
    lsdsim(grid_size = 1, time = 50, 
           vaccination = e, 
           culling = e, 
           quarantine = e, 
           insecticide = e, 
           interv_delay = 1
    )
  )
)
exp_3_res <- cbind.data.frame(interv = interv_vec, Reduce(rbind, exp_3_time))
dotchart(
  exp_3_res$elapse, pch = 20, 
  main = "lsddim: scaling with grid size", 
  xlab = "Runtime (s)",
  ylab = "Simulation (red = intervention)",
  col = interv_vec + 1
)
