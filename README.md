A stochastic meta-population model of LSD epidemics
================
Thibaut Jombart
2025-12-31

This repository contains the **R** package *lsdsim* which implements a
stochastic, meta-population compartmental model of Lumpy Skin Disease
(LSD), including the simulation of different types of interventions. It
supports the paper titled “*To cull or not to cull: a model-based
evaluation of response strategies against Lumpy Skin Disease outbreaks*”
by Jombart and Abbate (2025). This paper is currently under peer-review,
and being processed on bioRxiv. Meanwhile, a copy of the manuscript and
supporting material is provided on this repository (folder *paper*).

Direct links:

- [manuscript](https://github.com/thibautjombart/lsdsim/blob/main/paper/manuscript_1.0.pdf)
- [supporting
  material](https://github.com/thibautjombart/lsdsim/blob/main/paper/supplementary_material.pdf)

The repository also contains an an instance of
[*reportfactory*](https://www.reconverse.org/reportfactory/) to run
epidemic simulations in a reproducible way.

This repository is still under development. In particular, documentation
and testing will be improved in early 2026. If you intend on using its
content, do not hesitate to contact me at <thibautjombart@gmail.com>.

## Getting started

### Downloading the repository

You first need to clone or download the content of this repository, and
open an R session inside the main folder.

To clone the repository, use:

``` bash
git clone git@github.com:thibautjombart/lsdsim.git
```

If you have already cloned the repository, make sure you have the latest
version by opening a terminal in the main folder and typing:

``` bash
git pull
```

If you use [Rstudio](https://posit.co/download/rstudio-desktop/), the
easiest way to open a session is via opening the file
*lsdv_simulations.Rproj*.

To install *lsdsim* you will need to have *devtools* installed on your
computer and to be able to compile R packages from the source. You will
also need to have the *reportfactory* package installed. The latter will
be handy to handle dependencies of the Rmarkdown documents stored in
*report_sources/*.

### Installing dependencies

The following code should install all needed packages:

``` r
## install and load lsdsim
if (!required(devtools)) install.packages("devtools")
devtools::load_all()

## install reportfactory, install Rmd deps
if (!required(reportfactory)) install.packages("reportfactory")
reportfactory::install_deps()
```

## Generating results

### Published results

Results presented in the initial publications are stored as:

- *outputs/res_2025-12-29.qs*: R object storing summaries of 6,000
  simulations
- *outputs/simulations_final.html*: report of the 6,000 simulations

### Generating new results

The *reportfactory* can be used to generate a report running simulations
of epidemics with different types of responses, and storing the output
in a time-stamped folder in *outputs/*. The sources of the report itself
are in *report_sources/simulations.Rmd*. Users can modify this report as
they see fit, or add a new report to the *report_sources* folder to run
new simulations.

The *simulations.Rmd* report has the following parameters:

- n_rep: an `integer` indicating the number of replicates to generate;
  in our original study, we used 1,000, but beware that this will take
  hours to run on a standard desktop; defaults to 10
- run_sim: a `logical` indicating whether new simulations should be run;
  if `FALSE`, the latest simulation results stored inside `outputs` will
  be used to generate the report; defaults to `TRUE`

For instance, if you want to generate new results with 30 replicates,
you can type:

``` r
reportfactory::compile_reports(params = list(n_rep = 30, run_sim = TRUE))
```

``` r
>>> Compiling report: simulations
      - with parameters: n_rep = 30, run_sim = TRUE
All done!
```

## The *lsdsim* package

The *lsdsim* package was started around Christmas time to provide a
platform for testing different responses to LSD epidemics, in the wake
of mass culling imposed in France which led to protests and social
unrest. In a nutshell: it works, all main features are unit-tested, but
it is still under development. At this stage:

- unit-testing covers all core features including epidemic spread and
  interventions; full test coverage will be ensured later
- documentation is largely lacking, but is underway
- CRAN checks are not passing (mostly due to the lack of documentation)
  but the package does install on current R versions
- a list of all core functions of the package will be provided here
- continuous integration will be implemented via github actions
- a dedicated website will be provided using *pkgdown*
