# Change connectivity matrix in quarantine

This function changes the delta matrix used in
[lsdsim](https://thibautjombart.github.io/lsdsim/reference/lsdsim.md) to
implement quanrantine, defined as a reduction in outward-going
transmission. The strategy is to produce a complete new delta matrix
where all populations are in quarantine. Only specific rows of this
matrix will then be used during the simulation as needed.

## Usage

``` r
add_quarantine(delta, efficacy_in, efficacy_out)
```

## Arguments

- delta:

  the connectivity matrix used in
  [lsdsim](https://thibautjombart.github.io/lsdsim/reference/lsdsim.md)

- efficacy_in:

  the relative reduction in transmission within herds

- efficacy_out:

  the relative reduction in outward-going transmission

## Author

Thibaut Jombart <thibautjombart@gmail.com>
