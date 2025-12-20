test_that(
  "add_quarantine produces expected results", 
  {
    delta <- make_delta(make_grid(2), 0.05)
    
    ## no change when efficacy is 0
    expect_equal(delta, add_quarantine(delta, 0, 0))
    
    ## off-diagonal terms decrease by the right amount
    off_terms <- delta > 0
    diag(off_terms) <- FALSE
    res <- add_quarantine(delta, 0, 0.2)
    ratio <- res[off_terms] / delta[off_terms]
    expect_equal(ratio, rep(0.8, 8))
    
    ## diagonal terms decrease by the right amount
    res <- add_quarantine(delta, 0.3, 0)
    ratio <- diag(res) / diag(delta)
    expect_equal(ratio, rep(0.7, 4))
    
  }
)