# Draw proportion of asymptomatic for LSD

This random number generator draws the proportion of asymptomatic cases
for LSD using literature-driven parameters and beta distributions. See
[`rbeta()`](https://rdrr.io/r/stats/Beta.html) for details on means and
variances. .

## Usage

``` r
rlsd_pasymp(n, mu = 0.5, sd = 0.05)
```

## Source

TBC

## Arguments

- n:

  the number of values to draw

- mu:

  the mean of the distribution; defaults to 0.5

- sd:

  the standard deviation of the distribution; defaults to 0.05

## Author

Thibaut Jombart <thibautjombart@gmail.com>

## Examples

``` r
hist(
  rlsd_pasymp(1e5), 
  main = "Proportion of asymptomatic cases", 
  xlab = "P (asymptomatic)"
)

```
