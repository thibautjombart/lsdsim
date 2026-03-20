# Function to simulate an LSD outbreak in a metapopulation

This function simulates an outbreak of Lumpy skin disease in a
meta-population arranged on a regular grid of size `grid_size`, and a
diffusion process defined by a `delta` matrix.

## Usage

``` r
lsdsim(
  grid_size = 1,
  time = 365,
  beta_I = 0.1,
  beta_A = 0.01,
  sigma = 1/7,
  gamma = 1/20,
  cfr = 0.1,
  pasymp = 0.5,
  mass_culling = FALSE,
  select_culling = FALSE,
  vaccination = FALSE,
  quarantine = FALSE,
  insecticide = FALSE,
  rate_cull = 1e+30,
  vacc_coverage = 0,
  vacc_efficacy = 0.65,
  quarant_efficacy_in = 0.2,
  quarant_efficacy_out = 0.9,
  insect_efficacy = 0.5,
  interv_delay = 1e+30,
  interv_release = 28,
  delta = NULL,
  diffusion = 0,
  ini_S = 0,
  ini_E = 0,
  ini_A = 0,
  ini_I = 0,
  ini_C = 0,
  ini_D = 0,
  ini_R = 0,
  ini_V = 0
)
```

## Author

Thibaut Jombart <thibautjombart@gmail.com>
