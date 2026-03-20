# Build a delta matrix from geographic coordinates

Here coordinates are provided as x and y in a 2-columns matrix

## Usage

``` r
make_delta(coords, diffusion = 0)
```

## Arguments

- coords:

  a 2-columns matrix of xy coordinates

- diffusion:

  the value of the diffusion to be used, defined as the fraction of the
  force of infection going outside each patch towards the neighbours

## Author

Thibaut Jombart <thibautjombart@gmail.com>
