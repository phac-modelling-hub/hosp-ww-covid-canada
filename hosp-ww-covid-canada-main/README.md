# hosp-ww-covid-canada

Code and data accompanying publication by Yusuf &amp; Champredon

## Data

The following data files are required to create the model:

1.  `nml-database.RDS`: wastewater concentration data
2.  `ari_hosp.rds`: hospital admission data
3.  `allprop.rds`: provincial variant proportion data
4.  `counts-by-prov.rds`: provincial case count data

The following data files are for the sensitivity analyses:

1.  `ww-flow.rds`: wastewater flow data for available sites

These files are stored in the `data/` folder.

## Parameter values

The following data files are used to parameterize the model:

1.  `cityprov.csv`: list of cities and provinces to include in the model.
2.  `model-params.csv`: model parameter values.
3.  `sim-params.csv`: model parameter values for the simulated dataset.

These files are stored in the `prm/` folder.

## Model object creation

Scripts are stored in the `src/` folder.

The script `data.R` must be run first, as it generates the dataframe object used for the model development.

Note : 

 - that the `CCT` library is needed for the data scripts to run: `remotes::install_github("TheZetner/cct")`
 - `data.R` takes about 30 min to run.

## Model execution

The file `main-infere.R` executes the full analysis using the model object.

**Warning**: If the number of MCMC iterations (in `prm/model.param.csv`) is kept at the default 5e5, the computation time will be very long (likely longer than 12 hours). You may want to reduce the number of iterations (to, say, 1000) to test if the code runs correctly on your computer.  


### Appendix scripts

-   `inference-hmc.R` executes an Hamiltonian MCMC algorithm using a simulated dataset.

-   `sensitivity-variant.R` executes a sensitivity analysis of the variant dominance definition.

-   `autocorrelation.R` tests for autocorrelation of the data

-   `compare_norm_vs_not.R` compares the MCMC model output when using normalized vs non-normalized wastewater data. Non-normalized wastewater data is used in the manuscript.

-   `leave-one-out.R` executes a leave-one-out cross validation of the model.

-   `leave-one-out-analysis.R` analyses the results of the leave-one-out cross validation.

-   `residuals.R` executes a residual analysis of the model posterior fit.

## Utility scripts

-   `utils.R` contains utility functions used in the data processing and model execution scripts. Script is sourced in the other execution scripts.

-   `func-data.R` contains a function to digests the model data prior to model execution.

-   `run-tables.R` generates the LaTeX tables for the manuscript.

-   `run-figures.R` generates the figures for the manuscript using functions stored in `figures.R`.

-   `lm.R` contains the naive linear regression function.

-   `model-hmc.R` contains functions for the Hamiltonian MCMC model execution.

-   `plot.R` contains several plotting functions for diagnostics and manuscript preparation.

-   `var.R` contains functions for the variant dominance definition.

## Outputs

-   `output/` contains the model output objects, including the model object and the MCMC posterior samples.

-   `figs/` contains the figures of the data and model preparation, as well as model outputs.

-   Tables and figures for the manuscript are stored in `tables/` and `figs/` subfolders of `ms/`, respectively.
