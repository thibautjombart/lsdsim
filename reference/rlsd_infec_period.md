# Draw LSD infectious period

This random number generator draws average infectious period for LSD
using literature-driven parameters and log-normal distributions. The
infectious period is defined as the time interval between the onset of
infectiousness in an infected cattle and recovery or death. Reported
values refer to virus isolation in cattle blood. Note that the value
returned and associated variability is a population average, and does
not reflect individual-level variation.

## Usage

``` r
rlsd_infec_period(n, mu = 8.8, sd = 0.05)
```

## Source

Gubbins S. Using the basic reproduction number to assess the risk of
transmission of lumpy skin disease virus by biting insects. Transbound
Emerg Dis. 2019;66: 1873–1883.

## Arguments

- n:

  the number of values to draw

- mu:

  the population average, defaults to 8.8 days

- sd:

  the standard deviation of the lognormal distribution, defaults to 0.05

## Author

Thibaut Jombart <thibautjombart@gmail.com>
