test_that(
  "lsdsim output has the right dimensions", {
    res <- lsdsim()
    expect_equal(dim(res), c(365, 10))
    
    res <- lsdsim(time = 1)
    expect_equal(dim(res), c(1, 10))
    
    res <- lsdsim(time = 10, grid_size = 5)
    expect_equal(dim(res), c(10, 10*5*5))
  }
)



test_that(
  "lsdsim output has the right column names", {
    res <- lsdsim(time = 10, grid_size = 2)
    exp_names <- lapply(
      c("S", "E", "I", "C", "D", "R", "V", "N", "new_I","status"), 
      paste, 1:4, sep = "_")
    exp_names <- unlist(exp_names)
    expect_equal(colnames(res), exp_names)
  }
)


test_that(
  "all compartments are zero", {
    res <- lsdsim(time = 10, grid_size = 2, ini_I = 0)
    expect_true(all(unlist(res[, 1:36]) == 0))
    expect_true(all(res[, grep("status", names(res))] == "naive"))
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
    E <- res[, grep("^E_", colnames(res))]
    I <- res[, grep("^I_", colnames(res))]
    
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
    I <- res[, grep("^I_", colnames(res))]
 
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
    I <- res[, grep("^I_", colnames(res))]
   
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


test_that(
  "response triggers when it should", {
    res <- lsdsim(time = 20, grid_size = 3, 
                  ini_S = 1e5, 
                  ini_I = c(1, rep(0, 9)),
                  interv_delay = 14, # response 14 days after 1st case
                  interv_release = 1e6 # response stays on
    )
    
    status <- res[, grep("status_", colnames(res))]

    ## expectation: 
    ## response in patch 1 on day 15, stays on
    ## no response anywhere else
    expect_equal(status[, 1], rep(c("naive", "response"), c(14, 6)))
    expect_true(all(status[, -1] == "naive"))
  }
)


test_that(
  "response stops when it should", {
    res <- lsdsim(time = 20, grid_size = 3, 
                  ini_S = 1e5, 
                  ini_I = c(1, rep(0, 9)),
                  gamma = 1e30, # immediately leaving I
                  beta = 0, # no transmission beyond first case
                  interv_delay = 4, # response 4 days after 1st case
                  interv_release = 10 # response stops 10 days after last case
    )
    
    status <- res[, grep("status_", colnames(res))]
    
    ## expectation: intervention from day 5 to 12
    expect_true(all(status[1:4, 1] == "naive"))
    expect_true(all(status[5:12, 1] == "response"))
    expect_true(all(status[13:20, 1] == "naive"))
  }
)


test_that(
  "status output is correct", 
  {
    ## expected result:
    ##
    ## - all naive except for time steps 3:5
    ## - patches 1, 5, 8 go into response on time steps 3:5
    ## - patches 2, 6, 9, 4, 7, 12 go into ring on time steps 3:5
    ##
    ## To check the spatial structure, use:
    ##  plot(make_grid(4), pch = "")
    ##  text(make_grid(4), labels = 1:16, col = ini_I + 1)
    
    ini_I <- rep(0, 16)
    ini_I[c(1, 5, 8)] <- 1
    res <- lsdsim(time = 10, 
                  grid_size = 4, 
                  ini_S = 100, 
                  ini_I = ini_I,
                  gamma = 1e30, # fast leaving I
                  beta = 0, # no transmission beyond first case
                  interv_delay = 2, # response 4 days after 1st case
                  interv_release = 3, # response stops 10 days after last case)
                  diffusion = 0.01
    )
    status <- res[, grep("status_", colnames(res))]
    expected_response <- ini_I > 0
    expected_ring <- c(2, 6, 9, 4, 7, 12)
    expect_true(all(status[-(3:5), ] == "naive"))
    expect_true(all(status[3:5, expected_response] == "response"))
    expect_true(all(status[3:5, expected_ring] == "ring"))
  }
)


test_that(
  "mass culling works as expected", {
    res <- lsdsim(time = 20, grid_size = 3, 
                  ini_S = 1000, 
                  ini_E = c(10, rep(0, 8)),
                  ini_I = c(8, rep(0, 8)),
                  ini_R = c(5, rep(0, 8)),
                  ini_V = c(123, rep(0, 8)),
                  sigma = 0, # no leaving E
                  gamma = 0, # no leaving I
                  beta = 0, # no transmission beyond first case
                  mass_culling = TRUE,
                  interv_delay = 4, # response 4 days after 1st case
                  interv_release = 10 # response stops 10 days after last case
    )
    
    C <- res[, grep("C_", colnames(res))]
    
    ## expectation: intervention from day 5 to 12
    ## mass culling takes place on day 5 in patch 1
    ## 
    expect_true(all(C[1:4, 1] == 0))
    expect_true(all(C[5:20, 1] == 1146))
    expect_equal(res[, "S_1"], rep(c(1000, 0), c(4, 16)))
    expect_equal(res[, "E_1"], rep(c(10, 0), c(4, 16)))
    expect_equal(res[, "I_1"], rep(c(8, 0), c(4, 16)))
    expect_equal(res[, "R_1"], rep(c(5, 0), c(4, 16)))
    expect_equal(res[, "V_1"], rep(c(123, 0), c(4, 16)))
  }
)


test_that(
  "selective culling works as expected", {
    res <- lsdsim(time = 20, grid_size = 3, 
                  ini_S = 1000, 
                  ini_E = c(10, rep(0, 8)),
                  ini_I = c(8, rep(0, 8)),
                  ini_R = c(5, rep(0, 8)),
                  ini_V = c(123, rep(0, 8)),
                  sigma = 0, # no leaving E
                  gamma = 0, # no leaving I
                  beta = 0, # no transmission beyond first case
                  select_culling = TRUE,
                  rate_cull = 1000, # immediate kill of all infected cattles once in response
                  interv_delay = 4, # response 4 days after 1st case
                  interv_release = 10 # response stops 10 days after last case
    )
    
    C <- res[, grep("C_", colnames(res))]
    
    ## expectation: intervention from day 5 to 12
    ## selective culling takes place on day 5 in patch 1
    ## only I are culled (C) the rest stays as is
    expect_true(all(C[1:4, 1] == 0))
    expect_true(all(C[5:20, 1] == 8))
    expect_equal(res[, "S_1"], rep(1000, 20))
    expect_equal(res[, "E_1"], rep(10, 20))
    expect_equal(res[, "I_1"], rep(c(8, 0), c(4, 16)))
    expect_equal(res[, "R_1"], rep(5, 20))
    expect_equal(res[, "V_1"], rep(123, 20))
  }
)



test_that(
  "vaccination works as expected", {
    res <- lsdsim(time = 10, grid_size = 10, 
                  ini_S = 1e5, 
                  ini_V = 0,
                  ini_I = 1,
                  interv_delay = 1, 
                  beta = 0,
                  vaccination = TRUE,
                  vacc_coverage = 0.1, # 10% vaccination per day
                  vacc_efficacy = 0.5 # 50% efficacy
    )
    
    S <- res[, grep("S_", colnames(res))]
    V <- res[, grep("V_", colnames(res))]
    
    ## expectation: about 5% of individuals vaccinated
    res <- V[2, ] / S[1, ] 
    expect_true(all(res > 0.04))
    expect_true(all(res < 0.06))
  }
)



test_that(
  "quarantine works as expected", 
  {
    ## expectation: 
    ## infection spreads everywhere, last step is all I
    set.seed(1)
    res <- lsdsim(grid_size = 3, time = 20, 
                  ini_S = 1e4,
                  ini_I = c(1, rep(0, 8)),
                  interv_delay = 1, # intervention 1 day after 1st case
                  quarantine = TRUE, 
                  quarant_efficacy_in = 0, # no effect of quarantine isinde farm
                  quarant_efficacy_out = 0, # no reduction in outwards transmission
                  sigma = 1e30, # fast E->I
                  gamma = 0, # no leaving I
                  beta = 10,
                  diffusion = 0.1 # 10% diffusion of infection
    )
    
    I <- res[, grep("^I_", colnames(res))]
    expect_true(all(I[20,] >= 1e4))
    
    ## expectation: 
    ## infection does not spread due to quarantine stopping all outwards 
    ## transmissions
    res <- lsdsim(grid_size = 3, time = 20, 
                  ini_S = 1e4,
                  ini_I = c(1, rep(0, 8)),
                  interv_delay = 1, # intervention 1 day after 1st case
                  quarantine = TRUE,
                  quarant_efficacy_out = 1, # no outwards transmission
                  sigma = 1/14, # E->I in about 14 days
                  gamma = 1/7, # disease lasts about 7 days
                  beta = 10,
                  diffusion = 0.01 # 1% diffusion of infection
    )
    
    I <- res[, grep("^I_", colnames(res))]
    expect_true(all(I[, 1] > 0))
    expect_true(all(I[, -1] == 0))
    
    
    
    ## expectation: 
    ## infection does not spread due to quarantine stopping all transmissions
    ## inside the farms; patch 1 stays at 1 case
    res <- lsdsim(grid_size = 3, time = 20, 
                  ini_S = 1e4,
                  ini_I = c(1, rep(0, 8)),
                  interv_delay = 1, # intervention 1 day after 1st case
                  quarantine = TRUE,
                  quarant_efficacy_in = 1, # no outwards transmission
                  quarant_efficacy_out = 0, # no outwards transmission
                  sigma = 1/14, # E->I in about 14 days
                  gamma = 0, # no leaving I
                  beta = 10,
                  diffusion = 0.01 # 1% diffusion of infection
    )
    
    I <- res[, grep("^I_", colnames(res))]
    expect_equal(I[, 1], rep(1, 20))
    expect_true(all(I[, -1] == 0))
    
    }
)



test_that(
  "insecticide works as expected", 
  {
    ## expectation: 
    ## infection stops because insecticide has perfect efficacy
    ## only 1 infected in pop 1, no infection elsewhere
    res <- lsdsim(grid_size = 3, time = 20, 
                  ini_S = 1e4,
                  ini_I = c(1, rep(0, 8)),
                  interv_delay = 1, # intervention 1 day after 1st case
                  insecticide = TRUE, 
                  insect_efficacy = 1,
                  sigma = 1e30, # fast E->I
                  gamma = 0, # no leaving I
                  beta = 1,
                  diffusion = 0.1 # 10% diffusion of infection
    )
    
    I <- res[, grep("^I_", colnames(res))]
    expect_true(all(I[, 1] == 1))
    expect_true(all(I[, -1] == 0))
  }
)
