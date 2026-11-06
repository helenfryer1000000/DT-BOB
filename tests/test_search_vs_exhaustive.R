#!/usr/bin/env Rscript
############################################################################
# Checks that the search finds the exact top-1 and top-m statistics.
#
# For data sets small enough to enumerate all 2^nu models, it compares the
# search's statistic for every single-variable test with the exact value
# from exhaustive enumeration, for m = 1, 3 and 10.
#
# Usage (from the repository folder): Rscript tests/test_search_vs_exhaustive.R
############################################################################

this_dir = local({
	a = commandArgs(trailingOnly=FALSE)
	f = sub("^--file=", "", a[grep("^--file=", a)])
	if(length(f)) dirname(normalizePath(f)) else getwd()
})
repo = dirname(this_dir)
source(file.path(repo, "bob", "doublethink-bob-functions.R"))

check = function(label, y, design, family, intercept, mu, h, ms=c(1, 3, 10)) {
	E = make_engine(y, design, family=family, intercept=intercept, mu=mu, h=h)
	t0 = Sys.time()
	ex = exhaustive_scores(E)
	t_ex = as.numeric(difftime(Sys.time(), t0, units="secs"))
	# A fresh engine for the search so it cannot read the exhaustive cache
	E2 = make_engine(y, design, family=family, intercept=intercept, mu=mu, h=h)
	t0 = Sys.time()
	best = find_best_model(E2, y, design, n_restarts=20, verbose=FALSE)
	elite = cache_top(E2, 50)
	out = c()
	for(m in ms) {
		worst = 0; nbad = 0
		for(j in seq_len(E$nu)) {
			r = run_test(E2, j, best, elite, m=m)
			e = exhaustive_test(ex, j, m=m)
			d = abs(r$logT - e$logT)
			worst = max(worst, d); if(d > 1e-6) nbad = nbad + 1
		}
		out = c(out, sprintf("m=%d: %d/%d tests exact (max |error| in log T = %.2g)", m, E$nu - nbad, E$nu, worst))
	}
	t_s = as.numeric(difftime(Sys.time(), t0, units="secs"))
	best_exact = max(ex$lpo)
	cat(sprintf("%-38s nu=%d n=%d  best model found: %s | %s | search %.1fs, exhaustive %.1fs\n",
		label, E$nu, length(y), ifelse(abs(best$lpo - best_exact) < 1e-8, "yes", "NO"), paste(out, collapse="; "), t_s, t_ex))
	invisible(nbad)
}

## 1. The bundled AMD example (Gaussian, no intercept), and its binary version
cfgdir = file.path(repo, "example")
d = bob_load_data(file.path(cfgdir, "columns.csv"), file.path(cfgdir, "features.csv"), file.path(cfgdir, "outcomes.csv"), file.path(cfgdir, "exclude_rows.csv"))
excl = readLines(file.path(cfgdir, "exclude_columns.csv"))
design = bob_build_design(d$a, excl)
check("AMD example, gaussian", d$y, design, "gaussian", FALSE, 0.1, 1)
db = bob_load_data(file.path(cfgdir, "columns.csv"), file.path(cfgdir, "features.csv"), file.path(cfgdir, "binary_outcomes.csv"), file.path(cfgdir, "exclude_rows.csv"))
check("AMD example, binomial", db$y, design, "binomial", TRUE, 0.1, 1)

## 2. Simulated correlated designs with several real effects
sim = function(n, nu, rho, beta, family, seed) {
	set.seed(seed)
	R = rho^abs(outer(1:nu, 1:nu, "-"))        # AR(1) correlation: neighbours strongly correlated
	X = matrix(rnorm(n * nu), n, nu) %*% chol(R)
	colnames(X) = paste0("x", 1:nu)
	eta = drop(X %*% beta)
	y = if(family == "gaussian") eta + rnorm(n) else rbinom(n, 1, plogis(-1 + eta))
	list(y=y, design=list(X=X, blocks=setNames(as.list(1:nu), colnames(X)), varnames=colnames(X), df=setNames(rep(1L, nu), colnames(X))))
}
b = c(0.04, 0, 0, 0.03, 0, 0, 0, 0.05, 0, 0, 0, 0, 0.02, 0, 0, 0)
s = sim(10000, 16, 0.8, b, "gaussian", 1); check("AR(0.8), 4 effects, gaussian", s$y, s$design, "gaussian", TRUE, 0.1, 1)
s = sim(10000, 16, 0.95, b, "gaussian", 2); check("AR(0.95), 4 effects, gaussian", s$y, s$design, "gaussian", TRUE, 0.1, 1)
s = sim(10000, 16, 0.5, rep(0, 16), "gaussian", 3); check("AR(0.5), grand null, gaussian", s$y, s$design, "gaussian", TRUE, 0.1, 1)
s = sim(5000, 12, 0.8, 3 * b[1:12], "binomial", 4); check("AR(0.8), 3 effects, binomial", s$y, s$design, "binomial", TRUE, 0.1, 1)
