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

## Arguments

- grid_size:

  The size of the regular grid to be used. The number of populations
  will be this number squared. Defaults to 1.

- time:

  An `integer` indicating the duration of the simulation, in days.
  Defaults to 365.

- beta_I:

  The rate of infection for symptomatic cases (I); defaults to 0.1.

- beta_A:

  The rate of infection for asymptomatic cases (A); defaults to 0.01.

- sigma:

  The inverse of the latency period; defaults to 1/7.

- gamma:

  The inverse of the duration of illness; defaults to 1/20.

- cfr:

  The case fatality ratio; defaults to 0.1 (i.e. 10%).

- pasymp:

  The average proportion of asymptomatic cases; defaults to 0.5 (i.e.
  50%).

- mass_culling:

  A `logical` indicating if mass culling should be used in the response.
  Only the farm where the response takes place is affected. Defaults to
  FALSE.

- select_culling:

  A `logical` indicating if selective culling should be used in the
  response. Only the farm where the response takes place is affected.
  Defaults to FALSE.

- vaccination:

  A `logical` indicating if vaccination should be used in the response.
  The farm where the response takes place and its neighbours are
  affected. Defaults to FALSE.

- quarantine:

  A `logical` indicating if quarantine should be used in the response.
  The farm where the response takes place and its neighbours are
  affected. Defaults to FALSE.

- insecticide:

  A `logical` indicating if insecticide should be used in the response.
  The farm where the response takes place and its neighbours are
  affected. Defaults to FALSE.

- rate_cull:

  The rate at which culling takes place once the response has been
  triggered. Defaults to 1e30, i.e. all animals are killed in 1 day.

- vacc_coverage:

  The vaccine coverage, i.e. the proportion of individuals getting
  vaccinated once the response has been triggered. Defaults to 0 - no
  vaccination.

- vacc_efficacy:

  The vaccine efficacy, i.e. the proportion of vaccinated individuals
  who effectively become immune. Defaults to 0.65 (i.e. 65%).

- quarant_efficacy_in:

  The relative reduction of transmission inside the farm due to
  quarantine. Defaults to 0.2 (i.e. 20%).

- quarant_efficacy_out:

  The relative reduction of transmission towards neighbouring farms due
  to quarantine. Defaults to 0.9 (i.e. 90%).

- insect_efficacy:

  The relative reduction of transmission due to insecticide spray.
  Defaults to 0.5 (i.e. 50%).

- interv_delay:

  The delay to intervention, defined as the number of days after the
  first symptomatic case for the intervention to start. Defaults to
  1e30, i.e. no intervention.

- interv_release:

  The number of days after the last symptomatic case after which the
  intervention stops. Defaults to 28 days.

- delta:

  The connectivity matrix defining the proportions of the forces of
  infection going from populations (in rows), to populations (in
  columns). Defaults to `NULL`, in which case the matrix is generated
  automatically using rook connectivity and the `diffusion` argument.

- diffusion:

  The proportion of the force of infection directed towards neighbouring
  populations, used the build the `delta` matrix if it is set to `NULL`.
  Defaults to 0, i.e. no connectivity (`delta` is the identity matrix).

- ini_S:

  The initial number of susceptible individuals, recycled as needed to
  match the total number of populations. Defaults to 0.

- ini_E:

  The initial number of exposed individuals, recycled as needed to match
  the total number of populations. Defaults to 0.

- ini_A:

  The initial number of asymptomatic cases, recycled as needed to match
  the total number of populations. Defaults to 0.

- ini_I:

  The initial number of symptomatic cases, recycled as needed to match
  the total number of populations. Defaults to 0.

- ini_C:

  The initial number of culled animals, recycled as needed to match the
  total number of populations. Defaults to 0.

- ini_D:

  The initial number of individuals who died from the disease, recycled
  as needed to match the total number of populations. Defaults to 0.

- ini_R:

  The initial number of recovered individuals, recycled as needed to
  match the total number of populations. Defaults to 0.

- ini_V:

  The initial number of individuals immunized through vaccination,
  recycled as needed to match the total number of populations. Defaults
  to 0.

## Author

Thibaut Jombart <thibautjombart@gmail.com>
