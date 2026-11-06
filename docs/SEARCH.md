# How the search works

Best of both and top-m need, for every test, the best model (or the *m* best models) on each side: with a tested variable, and without. With ν variables there are 2^ν models, so they cannot all be scored when ν is more than about 20. The search finds the best ones by optimisation. The objective is the exact log posterior odds of a model,

$$\log \mathrm{PO}_{s:\varnothing} = |s|\log\mu + \frac{\mathrm{df}_s}{2}\log\xi_n + (1-\xi_n)\,\ell_{s:\varnothing},$$

which is a log-likelihood with a penalty for each variable: an L0-penalised (best-subset) regression problem. The code is in `bob/doublethink-bob-functions.R`.

## Step 1: starting models from the lasso

The lasso (L1-penalised regression, R package `glmnet`) is fitted along its whole regularisation path. Each point on the path selects a set of variables. These sets are scored by $\log\mathrm{PO}$ and the best `lasso_nstarts` (default 5) become starting models, together with the empty model. The lasso is fast even with thousands of variables, and it points the search at the right part of model space. If `glmnet` is not installed, or `lasso_init="false"`, the search starts from the empty model only.

## Step 2: hill climbing with add, drop and swap moves

From each starting model the search repeatedly moves to the best **neighbour**, until no neighbour is better (a local maximum). The neighbours of a model $s$ are:

* **add:** $s$ plus one variable;
* **drop:** $s$ minus one variable;
* **swap:** $s$ with one variable replaced by another (this lets the search move between correlated variables, which add and drop alone do badly).

Scoring every add and swap exactly would need about $\nu\,|s|$ model fits per step. Instead, each candidate is first **screened** by a quick estimate of how much it would raise $\log\mathrm{PO}$:

* Gaussian family: the exact gain, computed for all candidates at once from $X^\top X$ and $X^\top y$ (no pass over the data);
* logistic and Poisson: the Rao score test for adding the variable to the current fit, also computed for all candidates at once.

Only the `search_candidates` best adds (default 10), and the `search_candidates` best swaps plus the best swap for each dropped variable, are then fitted exactly. Drops are always fitted exactly. Every fitted model is cached, so no model is fitted twice.

## Step 3: random restarts (iterated local search)

To escape local maxima, the best model found so far is perturbed by `search_perturb` random add or drop moves (default 3) and the climb is repeated, `search_restarts` times (default 20). A restart that finds a better model replaces the best model.

The result is the best model overall, $\hat s$, and a cache of every model scored on the way. The best `elite_size` models in the cache (default 50) are kept as extra starting points for step 4.

## Step 4: the two sides of each test

For a test of variables $\mathcal{V}$:

* **With side** (models containing a tested variable): the climb starts from $\hat s$ if it contains a tested variable, otherwise from $\hat s$ plus each tested variable, and also from the best elite model on this side. Moves that would remove the last tested variable are not allowed.
* **Without side** (models containing no tested variable): the climb starts from $\hat s$ with the tested variables removed, and from the best elite model on this side. Tested variables may not be added.

Because the best model with a variable and the best model without it are usually close to $\hat s$, these climbs are short.

For **top-m** (m ≥ 2) the search then scores the neighbours of each of the current *m* best models on each side, and repeats until the list of the *m* best stops changing (at most 5 rounds). In simulations, the runners-up nearly always lie within one or two add, drop or swap moves of the best model, which is why this works.

If any test finds a model better than $\hat s$, the search climbs from that model to a new $\hat s$ and reruns all the tests from it (at most 5 times). So every test is anchored on the same, best-known model.

## How the search was checked

**Against exhaustive enumeration.** `tests/test_search_vs_exhaustive.R` scores all 2^ν models and compares the exact statistic with the search's, for every single-variable test and for m = 1, 3 and 10:

| Data set | ν | n | Family | Tests exact (m = 1, 3, 10) |
|---|---|---|---|---|
| AMD example | 15 | 145 | Gaussian, no intercept | 15/15, 15/15, 15/15 |
| AMD example, binary outcome | 15 | 145 | logistic | 15/15, 15/15, 15/15 |
| AR(1) correlation 0.8, 4 effects | 16 | 10,000 | Gaussian | 16/16 each |
| AR(1) correlation 0.95, 4 effects | 16 | 10,000 | Gaussian | 16/16 each |
| AR(1) correlation 0.5, no effects | 16 | 10,000 | Gaussian | 16/16 each |
| AR(1) correlation 0.8, 3 effects | 12 | 5,000 | logistic | 12/12 each |

The search found the exact best model and the exact statistic (difference 0 in $\log T$) in every case.

**At scale.** `tests/test_scale.R` simulates n = 100,000 observations of ν = 1000 variables in blocks of 10 with AR(1) correlation 0.9 within each block, and 10 true effects. The default search took 8 seconds to find the best model and about 4 minutes for all 1000 top-3 tests on one core (Gaussian). A second search with four times as many candidates per move and 100 restarts found the same best model and the same statistic for every test whose statistic was within a factor $e^{5}$ of 1 or above.

**What is not guaranteed.** A local search can in principle miss the best model. The lasso starts, swap moves, restarts and the re-anchoring rule make this unlikely, and the checks above found no case of it, but there is no proof. For important results, rerun with a larger `search_candidates` and `search_restarts` (and a different `seed`) and check that the statistics do not change.
