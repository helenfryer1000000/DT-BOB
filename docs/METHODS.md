# Methods

This page defines what each method computes. The notation follows the Doublethink papers.

## The model and the prior

There are ν candidate variables. A **model** $s$ is a subset of them (equivalently a binary vector in $\{0,1\}^\nu$); $\varnothing$ is the grand null model with no variables. Each model is a generalised linear model (Gaussian, logistic or Poisson), with or without an intercept, as set in the config.

For each model:

* $|s|$ is the number of variables in $s$;
* $\mathrm{df}_s$ is its number of coefficients (a numeric variable has 1; a factor with $L$ levels has $L-1$);
* $\ell_{s:\varnothing}$ is the maximised log-likelihood of $s$ minus that of the grand null model (for the Gaussian family the variance is estimated by maximum likelihood, as in `glm`).

With $n$ observations, prior odds $\mu$ of including each variable (`mu_prior`) and effect-size scale $h$ (`h_prior`), let

$$\xi_n = \frac{h}{n+h}.$$

The posterior odds of model $s$ against the grand null are

$$\log \mathrm{PO}_{s:\varnothing} = |s|\log\mu + \frac{\mathrm{df}_s}{2}\log\xi_n + (1-\xi_n)\,\ell_{s:\varnothing}.$$

The first term is the prior odds. The other two are the log Bayes factor (Johnson's approximation, `calc.logBF` in `code/doublethink-functions.R`). This is exactly what the MCMC uses, so all four methods score models identically. (One exception: without an intercept, `glm` codes the first factor in a model with all its levels, so the new methods require `intercept="true"` when there are factors.) (Check: `method="exhaustive"` on the AMD example reproduces the "log10 Posterior odds" column of the exact results in [MCMC.md](MCMC.md) to two decimal places.)

## The tests

A test names a set $\mathcal{V}$ of $k$ **tested variables**. Its null hypothesis is that all of them have zero effect. By default there is one test per variable ($k=1$); a groups file adds tests of several variables together.

The models split into two sides:

$$\mathcal{A} = \{\text{models containing at least one tested variable}\},$$
$$\mathcal{O} = \{\text{models containing no tested variable}\}.$$

$\mathcal{O}$ has $2^{\nu-k}$ models.

### The three statistics

**Full Doublethink (MCMC, exhaustive).** The model-averaged posterior odds against the null:

$$\mathrm{PO}_{\mathcal{A}:\mathcal{O}} = \frac{\sum_{s\in\mathcal{A}}\mathrm{PO}_{s:\varnothing}}{\sum_{s\in\mathcal{O}}\mathrm{PO}_{s:\varnothing}}.$$

**Best of both (top-1).** The best model with a tested variable against the best model without:

$$T_1 = \frac{\max_{s\in\mathcal{A}}\mathrm{PO}_{s:\varnothing}}{\max_{s\in\mathcal{O}}\mathrm{PO}_{s:\varnothing}}.$$

**Top-m.** The sum of the $m$ largest posterior odds on each side:

$$T_m = \frac{\text{sum of the } m \text{ largest } \mathrm{PO}_{s:\varnothing},\ s\in\mathcal{A}}{\text{sum of the } m \text{ largest } \mathrm{PO}_{s:\varnothing},\ s\in\mathcal{O}}.$$

$T_1$ is best of both. If $m$ is at least the number of models on each side, $T_m = \mathrm{PO}_{\mathcal{A}:\mathcal{O}}$. The tables report $\log_{10}$ of the statistic in the column "log10 statistic".

A useful fact: let $\hat s$ be the best model overall. If a tested variable is not in $\hat s$, then $\hat s\in\mathcal{O}$ is also the best model without it, so $T_1\le 1$ and $\log_{10} T_1\le 0$. Only variables in the best model (or groups that meet it) can have $T_1>1$.

### P-values and the family-wise error rate

Each statistic $T$ (any of the three) is turned into a p-value by the Doublethink asymptotic formula:

$$p = \Pr\!\left(\chi^2_1 \ \ge\ 2\log\frac{T}{c\,\mu\sqrt{\xi_n}}\right),$$

with $p=1$ when $T\le c\,\mu\sqrt{\xi_n}$. Two choices of $c$ are reported:

* **"-log10 Adjusted p-value"**: $c=\nu$. This is the main result. Rejecting every test with adjusted $p\le\alpha$ controls the family-wise error rate (FWER) at $\alpha$ asymptotically, over all single-variable tests and groups together. This is the same adjustment as the MCMC tables use.
* **"-log10 Single-test p-value"**: $c=k$, the number of tested variables. It is valid for one pre-chosen test only.

Why the adjustment covers all tests at once: for all three statistics, shrinking the null (testing more variables) can only increase the statistic, because the "with" side gains models and the "without" side loses them. So each is a shortcut closed-testing procedure (Theorem 1 of the paper), and the FWER of the whole family is the false positive rate of one test with $c=\nu$.

The formula is asymptotic: the p-value is accurate when it is small and $n$ is large. It uses $\chi^2_1$ even for a factor with several coefficients, as the MCMC tables do. As in the MCMC tables, adjusted p-values above 0.02 (`pvalue_cutoff`) are shown as "-". The "Significant (FWER)" column marks adjusted $p\le$ `alpha` (default 0.01).

### What is guaranteed

* **Full Doublethink:** the asymptotic level and strong FWER control are proved in the paper, under its conditions (a fixed number of variables, true effects not too small relative to $\sqrt{\log n/n}$, and the regularity conditions of each theorem).
* **Best of both:** under the same conditions, $T_1$ has the same asymptotic false positive rate as the full statistic at the same threshold, and the same FWER control. The proof (in a note accompanying the revised paper, not yet published) replaces one lemma of the paper, calibration of the sum, by its counterpart for the maximum; every other step carries over unchanged.
* **Top-m:** asymptotically every $T_m$ has the same level. At finite $n$ the rigorous bound for $m\ge2$ is only about $m$ times the nominal level, because several models each below the threshold can add up to it; in simulations top-m sits close to nominal.
* **All methods:** the guarantees hold for the statistic itself. The search must find the right models; see [SEARCH.md](SEARCH.md) for how this was checked.

### Bayesian reading

The posterior probability that the null is true is

$$\Pr(\text{null}\mid\text{data}) = \frac{1}{1+\mathrm{PO}_{\mathcal{A}:\mathcal{O}}}.$$

This needs the full sums, which only the MCMC (or exhaustive enumeration) estimates. The search statistics give a lower bound on $\mathrm{PO}_{\mathcal{A}:\mathcal{O}}$, line by line:

$$\sum_{s\in\mathcal{A}}\mathrm{PO}_{s:\varnothing}\ \ge\ \text{sum of the } m \text{ largest on } \mathcal{A},$$
$$\sum_{s\in\mathcal{O}}\mathrm{PO}_{s:\varnothing}\ \le\ \frac{|\mathcal{O}|}{m}\times\text{sum of the } m \text{ largest on } \mathcal{O}\qquad\text{(every other term is at most the mean of the } m \text{ largest)},$$
$$\text{so}\qquad \mathrm{PO}_{\mathcal{A}:\mathcal{O}}\ \ge\ \frac{\min(m,\,2^{\nu-k})}{2^{\nu-k}}\,T_m .$$

(When $m\ge|\mathcal{O}|=2^{\nu-k}$ the "without" side is summed in full and the factor is 1.) So the Bayesian false discovery rate is at most $1/(1+\tau)$ if a test rejects only when $T_m\ge 2^{\nu-k}\tau/\min(m,2^{\nu-k})$ (for best of both, $2^{\nu-k}\tau$). The bound assumes the search found the true $m$ best models on the "without" side; missing one of them makes $T_m$ too large. With many variables this threshold is enormous, so in practice **read the search methods as frequentist tests** and use the MCMC for posterior probabilities.

## How closely the search statistics agree with the full test

These numbers come from simulations done for the accompanying methods work (the scripts are not part of this repository). All $2^\nu$ models were enumerated, so the full statistic $\mathrm{PO}_{\mathcal{A}:\mathcal{O}}$ was known exactly, and the search was not involved. They use a Gaussian model with known variance (this software estimates the variance) and $\tau=\mu=h=1$. To compare the statistics on your own data with ν ≤ 20, use `method="exhaustive"`, which reports all three. "Disagreement" is the number of data sets on which exactly one of the two tests rejects, as a share of the data sets on which the full test rejects.

**Equicorrelated designs, ν = 2 to 8, n = 100,000 (worst of 10 designs):**

| Statistic | Disagreement, no effects | Disagreement, one effect |
|---|---|---|
| Best of both (top-1) | 44% | 10% |
| Top-3 | 3.4% | 1.1% |

**The 15 AMD biomarkers' design (ν = 15, all tested together):**

| Statistic | No effects, n = 145 | No effects, n = 14,500 | One effect, n = 145 | One effect, n = 14,500 |
|---|---|---|---|---|
| Best of both (top-1) | 34% | 17% | 2.1% | 0.1% |
| Top-3 | 23% | 4.2% | 0.8% | 0.0% |
| Top-10 | 12.5% | 1.1% | 0.3% | 0.0% |

Under the null, almost all of the disagreement is the search statistic *not* rejecting when the full test does (in the AMD design, top-3 and top-10 never rejected when the full test did not). So the search statistics are conservative relative to the full test, and agreement improves with $m$ and with $n$.
