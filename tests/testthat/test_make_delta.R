test_that(
  "make_delta works well without cliques", 
  {
    ## diagonal delta with no diffusion
    delta <- make_delta(make_grid(3), 0)
    expect_equal(delta, diag(1, 9))
    
    ## diagonal delta with no diffusion, even with all cliques
    delta <- make_delta(make_grid(3), 0, 
                        p_clique = 1, clique_connect = 1)
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


test_that(
  "make_delta works with cliques", 
  {
    ## corner case: all cliques
    res <- make_delta(make_grid(5), diffusion = 0.3, 
                      p_clique = 1, clique_connect = 1)
    expect_equal(rowSums(res > 0), rep(25, 25))
    
    ## 20% cliques, fully connected, 50% diffusion
    ## 
    ## we check that between 15% and 25% of location are fully connected
    res <- make_delta(make_grid(30), diffusion = 0.5, 
                        p_clique = 0.2, 
                        clique_connect = 1)
    expect_equal(diag(res), rep(0.5, 900))
    p_cliques <- mean(rowSums(res > 0) == 900)
    expect_true(p_cliques > 0.15)
    expect_true(p_cliques < 0.25)
    
  }
)
