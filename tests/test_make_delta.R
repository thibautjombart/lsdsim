test_that(
  "make_delta produces correct outputs", 
  {
    xy <- make_grid(4)
    
    ## diagonal matrix when no diffusion
    res <- make_delta(xy, 0)
    expect_equal(res, diag(1, 16))
    
    ## diagonal is 1 - diffusion
    res <- make_delta(xy, 0.123)
    expect_equal(diag(res), rep(1 - 0.123, 16))
    
    ## matrix is row-standardizes
    expect_equal(rowSums(res), rep(1, 16))
  }
)