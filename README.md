# Deadwood cavity-nesting Hymenoptera analyses

This repository contains the R code used for the study *"Moisture
limits the diversity of nesting bees, wasps, and parasitoids in lying and
standing deadwood"*. (Martini et al., submitted).

The study examines cavity-nesting Hymenoptera and their associated natural
enemies collected from experimental deadwood traps at the BEF-China platform.
The experiment comprised 192 traps in 64 plots, with three traps per plot:

- `G`: lying deadwood in contact with the ground;
- `A`: standing deadwood accessible to ants; and
- `nA`: standing deadwood with ant exclusion.

The workflow includes sampling-completeness analyses, generalized linear
mixed models, sensitivity analyses, community-composition analyses, beta
diversity, and piecewise structural equation models.

## Data availability

Data and R scripts are stored on GitHub + Zenodo.
Data is also available on
[Figshare](https://doi.org/10.6084/m9.figshare.31979481.v2):

- `log_data.csv` contains one row per deadwood trap. Its 192 rows include the
  experimental design, plot and environmental covariates, ant observations,
  and trap-level host and parasitoid responses. Variables with imputed values
  are supplied alongside their corresponding complete cases.
- `deadwd_insect_data.csv` contains emergence records for
  hosts and their parasitoids. It is used to construct the community matrices.
- `log_data_metadata.csv` contains metadata for `log_data.csv`.

## Project structure

```text
.
├── README.md
├── deadwd_script.R
├── deadwd_functions.R
├── path_analysis.R
├── log_data.csv
├── deadwd_insect_data.csv
├── log_data_metadata.csv
├── renv.lock
├── renv/
├── .Rprofile
```

The scripts have the following roles:

- `deadwd_script.R` runs the primary models, sampling-completeness analyses,
  sensitivity analyses, and community analyses.
- `deadwd_functions.R` contains helper functions and is sourced automatically.
- `path_analysis.R` runs the full and parsimonious piecewise SEM workflows.

## Reproducing the R environment

The analyses were run with **R 4.6.1**. Exact package versions and sources are
recorded in `renv.lock`.

Open R or RStudio in this directory. The supplied `.Rprofile` should activate
the project automatically. If necessary, activate it manually:

```r
source("renv/activate.R")
```

Restore and verify the environment:

```r
renv::restore(retry = FALSE)
renv::status()
```

A successful restoration reports that the project is in a consistent state.
Some packages may require an appropriate compiler or operating-system library.

## Running the analyses

Keep the three data files in the project directory. Run the two executable
scripts in this order:

```r
source("deadwd_script.R", echo = TRUE)
source("path_analysis.R", echo = TRUE)
```

Both scripts clear the workspace and automatically source `deadwd_functions.R`.

The three `iNEXT` calls in `deadwd_script.R` use 100 bootstrap replications for
a faster run. The reported analyses used 1,000 replications, as noted beside
those calls in the script.

## Outputs

- `deadwd_script_output/` contains the main and sensitivity model summaries
  and the RDS inputs used for supplementary tables.
- `path_analysis_output/` contains `sem_table_results.rds`.

The two optional PDF-rendering blocks in `deadwd_script.R` are commented out by
default. Uncomment them to create PDF versions of the model summaries. These
blocks require `rmarkdown` and LuaLaTeX; they are not required for the analyses
or numerical results.

Random-number-dependent procedures use explicit seeds where relevant. The
`renv` lockfile reproduces the R package environment but does not manage system
components such as compilers or LaTeX distributions.
