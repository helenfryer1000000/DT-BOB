#!/usr/bin/env Rscript
############################################################################
# Timing and consistency at a realistic scale, where exhaustive enumeration
# is impossible. Simulates a correlated design, runs the default search, then
# reruns every test with a much wider search (more candidates per move, more
# restarts) and reports how often the two agree.
#
# Usage: Rscript tests/test_scale.R [n] [nu] [family] [m]
#   defaults: n = 100000, nu = 1000, family = gaussian, m = 3
############################################################################

this_dir = local({
	a = commandArgs(trailingOnly=FALSE)
	f = sub("^--file=", "", a[grep("^--file=", a)])
	if(length(f)) dirname(normalizePath(f)) else getwd()
})
source(file.path(dirname(this_dir), "bob", "doublethink-bob-functions.R"))
args = commandArgs(trailingOnly=TRUE)
n = if(length(args) >= 1) as.numeric(args[1]) else 1e5
nu = if(length(args) >= 2) as.integer(args[2]) else 1000L
family = if(length(args) >= 3) args[3] else "gaussian"
m = if(length(args) >= 4) as.integer(args[4]) else 3L

set.seed(7)
# Blocks of 10 variables with AR(0.9) correlation inside each block
R = 0.9^abs(outer(1:10, 1:10, "-")); L = chol(R)
X = do.call(cbind, lapply(1:(nu / 10), function(b) matrix(rnorm(n * 10), n, 10) %*% L))
colnames(X) = paste0("x", 1:nu)
beta = rep(0, nu); eff = seq(5, nu, by=nu / 10); beta[eff] = 0.02
eta = drop(X %*% beta)
y = if(family == "gaussian") eta + rnorm(n) else rbinom(n, 1, plogis(-1 + 3 * eta))
design = list(X=X, blocks=setNames(as.list(1:nu), colnames(X)), varnames=colnames(X), df=setNames(rep(1L, nu), colnames(X)))
cat("n =", n, " nu =", nu, " family =", family, " m =", m, " true effects:", paste(colnames(X)[eff], collapse=" "), "\n")

E = make_engine(y, design, family=family, intercept=TRUE, mu=0.1, h=1)
t0 = Sys.time()
best = find_best_model(E, y, design, n_restarts=20, verbose=TRUE)
elite = cache_top(E, 50)
cat("Best model:", paste(colnames(X)[best$model], collapse=" "), "\n")
t1 = Sys.time()
res = lapply(1:nu, function(j) run_test(E, j, best, elite, m=m))
t2 = Sys.time()
cat(sprintf("Global search %.1fs, %d tests %.1fs, models scored %d\n",
	as.numeric(difftime(t1, t0, units="secs")), nu, as.numeric(difftime(t2, t1, units="secs")), length(ls(E$cache))))
logT = vapply(res, `[[`, 0, "logT")
p = doublethink_p(logT, nu, 0.1, E$xi)
cat("Variables with adjusted p < 0.01:", paste(colnames(X)[p < 0.01], collapse=" "), "\n")

# Wider search on every test whose statistic could matter (log T > -5)
chk = which(logT > -5)
best2 = find_best_model(E, y, design, ncand=40, n_restarts=100, verbose=FALSE)
elite2 = cache_top(E, 200)
d = vapply(chk, function(j) run_test(E, j, best2, elite2, m=m, ncand=40)$logT - logT[j], 0)
cat(sprintf("Wider search: same best model: %s; %d tests rechecked, %d changed by more than 1e-6 (max change in log T %.3g)\n",
	identical(best$model, best2$model), length(chk), sum(abs(d) > 1e-6), max(abs(d))))
