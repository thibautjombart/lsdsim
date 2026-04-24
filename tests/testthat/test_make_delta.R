test_that(
  "make_delta works well without cliques", 
  {
    ## diagonal delta with no diffusion
    delta <- make_delta(make_grid(3), 0)
    expect_equal(delta, diag(1, 9))
    
    ## delta for 3x3 grid with 10% diffusion
    xy <- make_grid(3)
    res <- make_delta(xy, 0.1)
    expect_equal(diag(res), rep(0.9, 9))
    
    ## non-neighbours are all further than 1 unit
    expect_true(all(as.matrix(dist(xy))[res == 0] > 1))

    ## neighbours are all 1 unit away from each other at most
    expect_true(all(as.matrix(dist(xy))[res > 0] <= 1))
  }
)