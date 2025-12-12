test_that(
  "lsdsim output has the right dimensions", {
    res <- lsdsim()
    expect_equal(dim(res), c(365, 7))
    
    res <- lsdsim(time = 1)
    expect_equal(dim(res), c(1, 7))
    
    res <- lsdsim(time = 10, grid_size = 5)
    expect_equal(dim(res), c(10, 7*5*5))
  }
)



test_that(
  "lsdsim output has the right column names", {
    res <- lsdsim(time = 10, grid_size = 2)
    exp_names <- lapply(
      c("S", "E", "I", "C", "D", "R", "V"), 
      paste, 1:4, sep = "_")
    exp_names <- unlist(exp_names)
    expect_equal(colnames(res), exp_names)
  }
)



test_that(
  "no case when no initial infection", {
    res <- lsdsim(time = 10, grid_size = 2, ini_I = 0)
    expect_true(all(unlist(res) == 0))
  }
)


test_that(
  "no further infection when beta = 0", {
    res <- lsdsim(time = 10, grid_size = 2, 
                  ini_S = 5000, ini_I = 100, 
                  beta = 0, # no infection
                  gamma = 1e30 # immediately leaving I
                  )
    I <- res[, grep("I_", colnames(res))]
    
    expect_true(all(I[1,] == 100))
    expect_true(all(I[-1,] == 0))
    
  }
)


test_that(
  "all infected when beta is high", {
    res <- lsdsim(time = 10, grid_size = 3, 
                  ini_S = 5000, ini_I = 1, 
                  beta = 1e30, # super high infection
                  sigma = 1e30, # super fast E->I
                  gamma = 0 # no leaving I
    )
    E <- res[, grep("E_", colnames(res))]
    I <- res[, grep("I_", colnames(res))]
    
    expect_true(all(E[2,] == 5000))
    expect_true(all(E[-2,] == 0))
    expect_true(all(I[-(1:2),] == 5001))
    
  }
)


test_that(
  "infections do not spread when no diffusion", {
    res <- lsdsim(time = 10, grid_size = 3, 
                  ini_S = 5000, ini_I = c(1, rep(0, 8)), 
                  beta = 1e30, # super high infection
                  sigma = 1e30, # super fast E->I
                  gamma = 0 # no leaving I
    )
    I <- res[, grep("I_", colnames(res))]
 
    expect_true(all(I[3:10, 1] == 5001))
    expect_true(all(I[, -1] == 0))
    
  }
)


test_that(
  "infections spread with diffusion", {
    res <- lsdsim(time = 10, grid_size = 3, 
                  ini_S = 5000, ini_I = c(1, rep(0, 8)), 
                  beta = 1e30, # super high infection
                  sigma = 1e30, # super fast E->I
                  gamma = 0, # no leaving I
                  diffusion = 0.1 # 10% diffusion
    )
    I <- res[, grep("I_", colnames(res))]
   
    expect_true(all(I[10, ] >= 5000))
    
  }
)


test_that(
  "deaths are within expected ratios", {
    res <- lsdsim(time = 10, grid_size = 10, 
                  ini_S = 1e5, 
                  ini_I = c(10, rep(0, 8)), 
                  beta = 1e30, # super high infection
                  sigma = 1e30, # super fast E->I
                  gamma = 1e30, # super fast I->...
                  cfr = 0.4, # 40% mortality
                  diffusion = 0.1 # 10% diffusion
    )
    R <- res[, grep("R_", colnames(res))]
    D <- res[, grep("D_", colnames(res))]
    res <- as.numeric(na.omit(D/(R+D)))
    expect_true(all(res > 0.39))
    expect_true(all(res < 0.41))
  }
)



