# Draw LSD latent period

This random number generator draws average latent period for LSD using
literature-driven parameters and log-normal distributions. The latent
period is defined as the time interval between infection of a cattle,
and this cattle becoming infectious. Reported values refer to virus
isolation in cattle blood. Note that the value returned and associated
variability is a population average, and does not reflect
individual-level variation.

## Usage

``` r
rlsd_latent(n, mu = 8.1, sd = 0.05)
```

## Source

Gubbins S. Using the basic reproduction number to assess the risk of
transmission of lumpy skin disease virus by biting insects. Transbound
Emerg Dis. 2019;66: 1873–1883.

## Arguments

- n:

  the number of values to draw

- mu:

  the population average, defaults to 8.1 days

- sd:

  the standard deviation of the lognormal distribution, defaults to 0.05

## Author

Thibaut Jombart <thibautjombart@gmail.com>
