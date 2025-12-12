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
  "lsdsim output has no case when no initial infection", {
    res <- lsdsim(time = 10, grid_size = 2, ini_I = 0)
    expect_true(all(unlist(res) == 0))
  }
)


test_that(
  "lsdsim output no further infection when beta = 0", {
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
