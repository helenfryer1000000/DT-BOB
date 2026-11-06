# Usage

## Installing

**Best of both, top-m and exhaustive** need only R (3.6 or later). The package `glmnet` is optional; the search uses its lasso path for starting models:

```
Rscript -e 'install.packages("glmnet")'
```

`ncores > 1` uses the base R package `parallel` (forking, so not on Windows).

**The MCMC** needs R, Python, SLURM and the packages listed in `code/R_packages.txt` and `code/Python_packages.txt`. The easiest route is the Docker image: see [MCMC.md](MCMC.md). The image also installs `glmnet`, so every method runs inside it. In the container the repository is at `/doublethink-mcmc`: run the new methods from that folder, and set `workdir="/workspace"` in the config so the results land in the folder shared with your computer (the example configs use `workdir="output"`, which stays inside the container).

## Running

```
./doublethink.sh <config file>
```

The `method` line of the config chooses the method. `doublethink.sh` passes MCMC configs to `code/qsub_wrapper.sh` and the others to `bob/run-bob.R`, which can also be called directly:

```
Rscript bob/run-bob.R <config file>
```

Relative paths in the config are read from the folder you run in.

## Input files

The same files as the MCMC (see `example/` for working copies):

| Config entry | File |
|---|---|
| `input_filename` | CSV of the variables, one row per individual, with a header row. The first column is an unnamed row number (as written by R's `write.csv`); it is dropped. Must contain the ID column `idcol`. |
| `columns_filename` | CSV with no header: each row gives a column name of the input file and its type, `double`, `float`, `integer` or `factor`. Must list every column of the input file except the first. |
| `cases_filename` | CSV with columns `idcol` and `phenocol`: the outcome for each ID (0/1 for logistic regression). |
| `exclude_filename` | IDs to leave out, one per line (optional). |
| `exclude_columns` | Variables to leave out, one per line (optional). |
| `groups_filename` | Optional, best of both / top-m / exhaustive only. One test per line: comma-separated variable names tested together, for example `ApoB,LDL.D`. |

Data preparation is as for the MCMC: missing values of numeric variables are replaced by the column mean; negative factor levels are recoded as `-999`; factors with one level or more than `nlevmax` levels are dropped. Factors enter as treatment-contrast dummy variables (a factor with L levels has L − 1 coefficients) but are added or dropped as one variable.

## Config file

The config is read by both R and bash, so keep the `name=value` form with no spaces around `=`. Lines starting with `#` are comments.

### Common to all methods

| Entry | Default | Meaning |
|---|---|---|
| `outname` | required | Name of the run. Results go to `workdir/outname/`. |
| `workdir` | required | Folder for results. |
| `method` | `"mcmc"` | `"best_of_both"`, `"top_m"`, `"exhaustive"` or `"mcmc"`. |
| `input_filename`, `columns_filename`, `cases_filename`, `exclude_filename`, `exclude_columns` | | See above. |
| `idcol`, `phenocol` | `"eid"`, `"pheno"` | ID and outcome column names. |
| `glm_family` | `"binomial"` | `"gaussian"`, `"binomial"` (logistic) or `"poisson"` for the new methods. |
| `intercept` | `"true"` | Include an intercept in every model. |
| `mu_prior` | required | Prior odds of including each variable (μ). |
| `h_prior` | required | Prior effect-size scale (h): 1 is AIC-like, n is BIC-like. |
| `nlevmax` | 50 | Factors with more levels are dropped. |
| `controls_subsample`, `seed_subsample` | 0, 1 | Optional subsampling of controls (0 = keep all), as in the MCMC. |

### Best of both, top-m and exhaustive

| Entry | Default | Meaning |
|---|---|---|
| `top_m` | 3 | m for `method="top_m"`, from 2 to 10. `method="exhaustive"` also uses it for its exact top-m column. |
| `groups_filename` | none | Extra tests of groups of variables (see above). |
| `test_variables` | `"all"` | `"all"` tests every variable. `"best"` tests only the variables in the 10 best models found, which is quicker. For best of both every other variable has a statistic of at most 1 (see [METHODS.md](METHODS.md#the-three-statistics)), so this rarely loses a discovery, but check the p-value a statistic of 1 would give with your μ, h and n before relying on it. |
| `search_candidates` | 10 | Adds and swaps fitted exactly per step, after screening. Larger is slower and more thorough. |
| `search_restarts` | 20 | Random restarts of the search for the best model. |
| `search_perturb` | 3 | Random moves per restart. |
| `lasso_init` | `"true"` | Start from models on the lasso path (needs `glmnet`). |
| `lasso_nstarts` | 5 | Number of lasso starting models. |
| `max_model_size` | `Inf` | Largest number of variables in a model. |
| `elite_size` | 50 | Best models kept as extra starting points for each test. |
| `seed` | 1 | Random seed (the search is reproducible given the seed). |
| `ncores` | 1 | CPU cores used to run the tests in parallel. |
| `alpha` | 0.01 | FWER level for the "Significant (FWER)" column. |
| `pvalue_cutoff` | 0.02 | Adjusted p-values above this are shown as "-". |

`method="exhaustive"` ignores the search settings and stops if ν > 20.

With factor variables, the new methods need `intercept="true"`. A missing `exclude_filename` or `exclude_columns` file is an error; names in `exclude_columns` that are not in the data give a warning.

### MCMC only

`sourcedir`, `R_lib_loc`, `python_loc`, `chain_numbers`, `niter`, `burnin`, `thinning`, `seed_offset`, `flat_likelihood`: see `code/README.md` and [MCMC.md](MCMC.md).

## Output files (best of both, top-m, exhaustive)

All in `workdir/outname/`, tab-separated:

**`results.bob.variables.tsv`**: one row per variable, sorted by the statistic.

| Column | Meaning |
|---|---|
| Test | The variable. |
| Variables tested | k (1 here). |
| In best model | Whether the variable is in the best model overall. |
| log10 statistic | $\log_{10}$ of $T_1$ (best of both), $T_m$ (top-m) or $\mathrm{PO}_{\mathcal{A}:\mathcal{O}}$ (exhaustive). Above 0 means the "with" side wins. |
| log10 best of both (exact), log10 top-m (exact) | Exhaustive only: the exact $T_1$ and $T_m$ (m = `top_m`), for comparison with the full statistic. The p-value columns use the full statistic. |
| -log10 Adjusted p-value | FWER-adjusted p-value over all tests (main result); "-" if above `pvalue_cutoff`. |
| Significant (FWER) | `*` if the adjusted p-value is at most `alpha`. |
| -log10 Single-test p-value | p-value for this test alone, not adjusted for multiple testing. |
| Effect size in best model with, Standard error | The variable's coefficient and its standard error in the best model containing it (on the scale of the input data; "-" for factors and groups). |
| Best model with, Best model without | The two models compared by best of both (the top model on each side for top-m). |
| log10 PO best with, log10 PO best without | Their posterior odds against the grand null. |

**`results.bob.groups.tsv`**: the same for each group in `groups_filename`.

**`results.bob.top-models.tsv`** (search methods): for every test, the m best models found on each side ("with" and "without") and their $\log_{10}$ posterior odds against the grand null.

**`results.bob.best-models.tsv`**: the best models found overall, with their $\log_{10}$ posterior odds and size.

**`results.bob.best-model-coefficients.tsv`**: coefficients and standard errors of the best model overall (not written when the best model is the grand null model, as in the binary example).

**`results.bob.log`**: settings, the best model, the number of models scored and the run time. A copy of the config is saved beside it.

### Reading the results

* A variable is a **discovery** when its adjusted p-value is small: the "Significant (FWER)" column marks those at `alpha`. The family-wise error rate over all tests in the run, groups included, is then controlled at `alpha` (asymptotically).
* For best of both, only variables in the best model can have a positive "log10 statistic". A variable outside it can still matter jointly with correlated variables: test them as a group.
* When two variables are strongly correlated, often neither is significant alone but the pair is: this is the same situation as the MCMC's "post-hoc groups". Use `groups_filename` to test such sets.
* Compare "Best model with" and "Best model without" to see which variable takes over when the tested one is removed.

## Run time and memory

Measured on one core:

| Data | Method | Time |
|---|---|---|
| AMD example (n = 145, ν = 15) | best of both, top-m | about 2 seconds |
| AMD example (n = 145, ν = 15) | exhaustive, Gaussian / logistic | about 4 / 20 seconds |
| Simulated, n = 100,000, ν = 1000, Gaussian | top-3, all 1000 tests | about 3 minutes (plus 10 seconds for the best model) |
| Simulated, n = 20,000, ν = 200, logistic | best of both, all 200 tests | about 20 minutes (3 minutes for the best model) |

The Gaussian family needs only $X^\top X$ and $X^\top y$, so its cost hardly depends on n once these are formed. Logistic and Poisson models refit the GLM for each model scored, so their cost grows with n. Roughly, a logistic run costs n × ν × (number of tests), so at n = 100,000 and ν = 1000 testing all 1000 variables would take more than a day on one core. For logistic models at that scale use `test_variables="best"` (tests only the variables in the best models, typically a few dozen), `ncores`, or both. Memory is dominated by the data matrix (8 bytes × n × number of design columns) plus up to about 200 MB of cached fits for logistic and Poisson models.

## Troubleshooting

* **"glm_family = binomial needs a 0/1 phenotype"**: the outcome has other values; use `glm_family="gaussian"` or recode.
* **"No variables left to analyse"**: every column was excluded; check `exclude_columns` and `columns_filename`.
* **"Some IDs in input_filename have no phenotype"**: every ID in the input file must appear in `cases_filename` (or be excluded).
* **A model with collinear variables** scores minus infinity and is never chosen, as in the MCMC.
* **Results change with the seed**: the search has missed a model in one of the runs. The search can miss models but never overstates a model's posterior odds, so for each test the better run is the one whose "log10 PO best with" and "log10 PO best without" are higher. Rerun with larger `search_candidates` and `search_restarts` until the results stop changing.
