# Doublethink: MCMC, best of both and top-m

Doublethink tests which variables (risk factors) have a direct effect on an outcome. Each test is simultaneously Bayesian (a posterior odds) and frequentist (a p-value that controls the family-wise error rate). This repository offers three ways to compute the test:

| Method | `method=` in the config | What it computes | How |
|---|---|---|---|
| **Best of both** | `"best_of_both"` | For each tested variable: the best model that contains it against the best model that does not | Search: lasso starting points, then a greedy add/drop/swap search with random restarts |
| **Top-m** | `"top_m"` with `top_m=2`…`10` | The sum over the *m* best models with the variable against the sum over the *m* best without | The same search, then a check of the neighbours of the *m* best models on each side |
| **Doublethink MCMC** | `"mcmc"` (or no `method` line) | The full model-averaged posterior odds, estimated by Markov chain Monte Carlo | The original [doublethink-mcmc](https://github.com/danny-wilson/doublethink-mcmc) code, unchanged, in `code/` |

A fourth method, `"exhaustive"`, scores every one of the 2^ν models. It gives the exact Doublethink answer (and exact best-of-both and top-m answers) when the number of variables ν is about 20 or fewer. It is there for small problems and for checking the search.

The idea behind best of both and top-m: the full Doublethink statistic sums the posterior odds of all 2^ν models, which needs MCMC when ν is large. In large samples the sums are dominated by their largest terms, so the ratio of the largest terms, or of the few largest, has the same asymptotic false positive rate and can be found by optimisation instead of sampling. [docs/METHODS.md](docs/METHODS.md) gives the formulas and what is and is not guaranteed.

## Quick start (no Docker needed for the new methods)

You need R (version 3.6 or later). The R package `glmnet` is optional but recommended: the search uses its lasso path for starting models.

```
git clone <this repository>
cd <repository folder>
Rscript -e 'install.packages("glmnet")'      # optional

./doublethink.sh example/example_best_of_both.cfg     # best of both, Gaussian AMD example
./doublethink.sh example/example_top3.cfg             # top-3
./doublethink.sh example/binary_best_of_both.cfg      # best of both, logistic regression
./doublethink.sh example/example_groups.cfg           # top-3, with groups of variables tested together
./doublethink.sh example/example_exhaustive.cfg       # exact answer by enumeration (15 variables)
```

Run these from the repository folder: the example configs use paths relative to it. Each run takes a few seconds and writes its results to `output/<outname>/`. The main table is `results.bob.variables.tsv`:

```
column -s $'\t' -t output/example_top3/results.bob.variables.tsv | less -S
```

To run the original MCMC, follow [docs/MCMC.md](docs/MCMC.md) (Docker and SLURM).

## Which method to use

* For **posterior probabilities** (the Bayesian reading), use the MCMC, or `exhaustive` when ν ≤ 20. The search methods give a frequentist test with the same asymptotic level; they do not give a posterior probability (see [docs/METHODS.md](docs/METHODS.md#bayesian-reading)).
* For **p-values with family-wise error control** when ν is large (hundreds to thousands of variables), use the search. Top-m with m = 3 is a good default. In simulations, best of both (m = 1) disagrees with the full Doublethink test more often than top-3, and larger m agrees more closely (numbers in [docs/METHODS.md](docs/METHODS.md#how-closely-the-search-statistics-agree-with-the-full-test)). Where they disagree under the null, the search statistics nearly always fail to reject when the full test rejects, so they err on the conservative side.
* The search is deterministic given the `seed`, needs no convergence checks, and runs on one machine without SLURM.

## Documentation

| File | Contents |
|---|---|
| [docs/USAGE.md](docs/USAGE.md) | Input files, every config option, output files and how to read them, run times, troubleshooting |
| [docs/METHODS.md](docs/METHODS.md) | The model, the prior, the three statistics, p-values, FWER, the Bayesian reading, agreement with the full test |
| [docs/SEARCH.md](docs/SEARCH.md) | How the search works, and how it was checked |
| [docs/MCMC.md](docs/MCMC.md) | The original doublethink-mcmc instructions (Docker, SLURM, MCMC outputs) |
| `code/README.md` | The original developer notes on the MCMC code and its config file |

## Repository layout

```
doublethink.sh            one entry point for every method (reads 'method' from the config)
bob/                      best of both, top-m and exhaustive (new)
  run-bob.R               reads a config, runs the method, writes the results
  doublethink-bob-functions.R   likelihood, posterior odds, search, tests
code/                     the original Doublethink MCMC code (unchanged)
example/                  the AMD example data, the original MCMC configs and new configs
tests/                    checks of the search against exhaustive enumeration, and a timing run
docs/                     documentation
Dockerfile, slurm.conf, start_slurm.sh   container for the MCMC (the new methods also run in it)
```

## Checks

`Rscript tests/test_search_vs_exhaustive.R` compares the search with exhaustive enumeration on six data sets (the AMD example as Gaussian and logistic regressions, and simulated correlated designs with 12 to 16 variables, up to n = 10,000). For every single-variable test and for m = 1, 3 and 10 the search found the exact statistic. `Rscript tests/test_scale.R` times a Gaussian run with n = 100,000 and ν = 1000 and re-checks it with a much wider search. See [docs/SEARCH.md](docs/SEARCH.md#how-the-search-was-checked).

## Citation and credits

The MCMC code, the example data and the Docker set-up are from [doublethink-mcmc](https://github.com/danny-wilson/doublethink-mcmc) by Daniel J. Wilson and Nicolas Arning, released with

> Arning N, Fryer HR, Wilson DJ (2026). Identifying direct risk factors in UK Biobank via simultaneous Bayesian-frequentist model-averaged hypothesis testing using Doublethink. *PNAS* 123(1): e2514138122.

The theory is in Fryer HR, Arning N, Wilson DJ, *Doublethink: simultaneous Bayesian-frequentist model-averaged hypothesis testing*. The best-of-both and top-m code in `bob/` was added to this fork.

MIT licence: see [LICENSE](LICENSE).
