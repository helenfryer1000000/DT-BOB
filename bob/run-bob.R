#!/usr/bin/env Rscript
############################################################################
# run-bob.R: Doublethink by search ("best of both", top-m) or by exhaustive
# enumeration. Reads the same config file as the MCMC version.
#
# Usage: Rscript bob/run-bob.R config_file
#
# See README.md and docs/USAGE.md.
############################################################################

help = "Usage: Rscript run-bob.R config_file"
args = commandArgs(trailingOnly=TRUE)
if(length(args) != 1) { cat(help, "\n"); stop("Incorrect usage") }
cfg = normalizePath(args[1])
source(cfg)

# Locate this script's folder so the functions can be found wherever the repo is.
this_dir = local({
	a = commandArgs(trailingOnly=FALSE)
	f = sub("^--file=", "", a[grep("^--file=", a)])
	if(length(f)) dirname(normalizePath(f)) else getwd()
})
source(file.path(this_dir, "doublethink-bob-functions.R"))

## Settings and defaults ###################################################
default = function(name, value) if(!exists(name, envir=.GlobalEnv, inherits=FALSE)) assign(name, value, envir=.GlobalEnv)
default("method", "best_of_both")      # best_of_both | top_m | exhaustive
default("top_m", 3)                     # m for method="top_m" (2 to 10)
default("idcol", "eid")
default("phenocol", "pheno")
default("glm_family", "binomial")
default("intercept", "true")
default("controls_subsample", 0)
default("seed_subsample", 1)
default("seed", 1)
default("nlevmax", 50)
default("exclude_filename", "")
default("exclude_columns", "")
default("groups_filename", "")          # optional: one test per line, comma-separated variable names
default("test_variables", "all")        # all | best (only variables in the best model and the best runner-up models)
default("search_candidates", 10)        # adds/swaps scored exactly per move, after screening
default("search_restarts", 20)          # random restarts of the global search
default("search_perturb", 3)            # random moves per restart
default("lasso_init", "true")           # use the lasso path (glmnet) for starting models
default("lasso_nstarts", 5)
default("max_model_size", Inf)
default("elite_size", 50)               # best models kept from the global search, used as extra starts
default("ncores", 1)
default("alpha", 0.01)                  # FWER level for the "Significant" column (keep <= pvalue_cutoff)
default("pvalue_cutoff", 0.02)          # p-values above this are shown as '-', as in the MCMC tables

method = tolower(method)
if(!method %in% c("best_of_both", "top_m", "exhaustive")) stop("method must be best_of_both, top_m or exhaustive (use doublethink.sh for mcmc)")
m = if(method == "best_of_both") 1 else as.integer(top_m)
if(method == "top_m" && (m < 2 || m > 10)) stop("top_m must be between 2 and 10")
intercept = toupper(as.character(intercept)) == "TRUE"
lasso_init = toupper(as.character(lasso_init)) == "TRUE"
mu_prior = as.numeric(mu_prior); h_prior = as.numeric(h_prior)
if(!(mu_prior > 0)) stop("mu_prior must be positive")
if(!(h_prior > 0)) stop("h_prior must be positive")
set.seed(as.numeric(seed))

outdir = file.path(workdir, outname)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
invisible(file.copy(cfg, file.path(outdir, basename(cfg)), overwrite=TRUE))
logfile = file.path(outdir, "results.bob.log")
logcon = file(logfile, open="wt")
say = function(...) { msg = paste0(...); cat(msg, "\n"); cat(msg, "\n", file=logcon); flush(logcon) }
t0 = Sys.time()
say("Doublethink (", method, if(method == "top_m") paste0(", m = ", m), ") started ", format(t0))

## Data ####################################################################
data = bob_load_data(columns_filename, input_filename, cases_filename,
	if(nzchar(exclude_filename)) exclude_filename else NULL,
	idcol=idcol, phenocol=phenocol, controls_subsample=as.numeric(controls_subsample), seed_subsample=seed_subsample)
colnames.remove = if(nzchar(exclude_columns) && file.exists(exclude_columns)) { x = readLines(exclude_columns); x[nzchar(x)] } else character(0)
design = bob_build_design(data$a, colnames.remove, as.numeric(nlevmax))
say("n = ", length(data$y), " observations, nu = ", length(design$blocks), " variables (", ncol(design$X), " design columns)")
if(glm_family == "binomial" && !all(data$y %in% c(0, 1))) stop("glm_family = binomial needs a 0/1 phenotype")

E = make_engine(data$y, design, family=glm_family, intercept=intercept, mu=mu_prior, h=h_prior)
nu = E$nu
say("mu = ", mu_prior, ", h = ", h_prior, ", xi = h/(n+h) = ", signif(E$xi, 4), ", family = ", glm_family, ", intercept = ", intercept)

## The tests ###############################################################
tests = lapply(seq_len(nu), function(j) j)
names(tests) = design$varnames
groups = list()
if(nzchar(groups_filename)) {
	gl = readLines(groups_filename); gl = gl[nzchar(trimws(gl))]
	for(line in gl) {
		vars = trimws(strsplit(line, ",", fixed=TRUE)[[1]])
		idx = match(vars, design$varnames)
		if(any(is.na(idx))) stop("Group '", line, "' names variables that are not analysed: ", paste(vars[is.na(idx)], collapse=", "))
		groups[[paste(vars, collapse="+")]] = sort(unique(idx))
	}
	say(length(groups), " groups read from ", groups_filename)
}
all_tests = c(tests, groups)

## Exhaustive ##############################################################
if(method == "exhaustive") {
	say("Scoring all 2^", nu, " models")
	ex = exhaustive_scores(E)
	ord = order(-ex$lpo)[1:min(50, length(ex$lpo))]
	best = list(model=which(ex$M[ord[1], ] == 1), lpo=ex$lpo[ord[1]])
	elite = lapply(ord, function(i) list(model=which(ex$M[i, ] == 1), lpo=ex$lpo[i]))
	res = lapply(all_tests, function(t) {
		r = exhaustive_test(ex, t, m=1)
		inA = rowSums(ex$M[, t, drop=FALSE]) > 0
		ia = which(inA)[which.max(ex$lpo[inA])]; io = which(!inA)[which.max(ex$lpo[!inA])]
		list(tested=t, logT=r$logPO, logT1=r$logT1,
			topA=list(keys=model_key(which(ex$M[ia, ] == 1)), lpo=ex$lpo[ia]),
			topO=list(keys=model_key(which(ex$M[io, ] == 1)), lpo=ex$lpo[io]))
	})
} else {
## Search ##################################################################
	say("Searching for the best model (lasso starts: ", lasso_init && requireNamespace("glmnet", quietly=TRUE), ", restarts: ", search_restarts, ")")
	best = find_best_model(E, data$y, design, ncand=search_candidates, n_restarts=search_restarts,
		perturb_moves=search_perturb, use_lasso=lasso_init, nlasso=lasso_nstarts, max_size=max_model_size)
	say("Best model (log10 PO vs null = ", signif(best$lpo / log(10), 6), "): ", if(length(best$model)) paste(design$varnames[best$model], collapse=", ") else "(null model)")
	run_all = function(best) {
		elite = cache_top(E, elite_size)
		todo = all_tests
		if(test_variables == "best") {
			keep = unique(unlist(lapply(elite[1:min(10, length(elite))], `[[`, "model")))
			todo = c(tests[keep], groups)
		}
		say("Running ", length(todo), " tests with m = ", m)
		f = function(t) run_test(E, t, best, elite, m=m, ncand=search_candidates, max_size=max_model_size)
		if(as.integer(ncores) > 1) parallel::mclapply(todo, f, mc.cores=as.integer(ncores)) else lapply(todo, f)
	}
	# If a test finds a model better than the best model, restart with it:
	# the best model is the anchor for every test.
	for(round in 1:5) {
		res = run_all(best)
		imp = which(vapply(res, `[[`, TRUE, "improved"))
		if(length(imp) == 0) break
		cand = lapply(res[imp], function(r) {
			k = if(r$topA$lpo[1] >= r$topO$lpo[1]) r$topA$keys[1] else r$topO$keys[1]
			list(model=key2model(k), lpo=max(r$topA$lpo[1], r$topO$lpo[1]))
		})
		nb = cand[[which.max(vapply(cand, `[[`, 0, "lpo"))]]
		best = climb(E, nb$model, SIDE_FREE, search_candidates, max_model_size)
		say("A test found a better model; restarting the tests from it (log10 PO = ", signif(best$lpo / log(10), 6), ")")
	}
	elite = cache_top(E, elite_size)
}

## Output ##################################################################
vn = function(k) { s = key2model(k); if(length(s) == 0) "(null model)" else paste(design$varnames[s], collapse=" + ") }
fmt = function(x, d=2) ifelse(is.na(x), "-", ifelse(is.finite(x), formatC(round(x, d), format="f", digits=d), as.character(x)))
stat_name = switch(method, best_of_both="best of both (top-1)", top_m=paste0("top-", m), exhaustive="exact model-averaged posterior odds")

rows = lapply(names(res), function(nm) {
	r = res[[nm]]
	k = length(r$tested)
	logT = r$logT
	p_adj = doublethink_p(logT, nu, mu_prior, E$xi)
	p_unadj = doublethink_p(logT, k, mu_prior, E$xi)
	cf = if(k == 1) E$coef(key2model(r$topA$keys[1])) else NULL
	est = se = NA
	if(!is.null(cf)) { w = which(cf$variable == design$varnames[r$tested]); if(length(w) == 1) { est = cf$estimate[w]; se = cf$se[w] } }
	data.frame(
		"Test" = nm,
		"Variables tested" = k,
		"In best model" = ifelse(any(r$tested %in% best$model), "yes", "no"),
		"log10 statistic" = fmt(logT / log(10)),
		"-log10 Adjusted p-value" = ifelse(p_adj <= pvalue_cutoff, fmt(-log10(p_adj)), "-"),
		"Significant (FWER)" = ifelse(p_adj <= as.numeric(alpha), "*", ""),
		"-log10 Single-test p-value" = fmt(-log10(p_unadj)),
		"Effect size in best model with" = fmt(est),
		"Standard error" = fmt(se),
		"Best model with" = vn(r$topA$keys[1]),
		"Best model without" = vn(r$topO$keys[1]),
		"log10 PO best with" = fmt(r$topA$lpo[1] / log(10)),
		"log10 PO best without" = fmt(r$topO$lpo[1] / log(10)),
		check.names=FALSE, stringsAsFactors=FALSE)
})
tab = do.call(rbind, rows)
tab$`.ord` = vapply(res, `[[`, 0, "logT")
is_group = tab$Test %in% names(groups)
write_tab = function(d, f) { d$`.ord` = NULL; write.table(d, file.path(outdir, f), sep="\t", quote=FALSE, row.names=FALSE) }
write_tab(tab[!is_group, ][order(-tab$`.ord`[!is_group]), ], "results.bob.variables.tsv")
if(any(is_group)) write_tab(tab[is_group, ][order(-tab$`.ord`[is_group]), ], "results.bob.groups.tsv")

# The m best models on each side of every test (long format)
if(method != "exhaustive") {
	long = do.call(rbind, lapply(names(res), function(nm) {
		r = res[[nm]]
		rbind(
			data.frame(Test=nm, Side="with", Rank=seq_along(r$topA$keys), Model=vapply(r$topA$keys, vn, ""), log10PO=fmt(r$topA$lpo / log(10), 3), stringsAsFactors=FALSE),
			data.frame(Test=nm, Side="without", Rank=seq_along(r$topO$keys), Model=vapply(r$topO$keys, vn, ""), log10PO=fmt(r$topO$lpo / log(10), 3), stringsAsFactors=FALSE))
	}))
	write.table(long, file.path(outdir, "results.bob.top-models.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}
# Best models overall
bm = data.frame(Rank=seq_along(elite), "log10 PO vs null"=fmt(vapply(elite, `[[`, 0, "lpo") / log(10), 3),
	Size=vapply(elite, function(e) length(e$model), 0L),
	Model=vapply(elite, function(e) vn(model_key(e$model)), ""), check.names=FALSE)
write.table(bm, file.path(outdir, "results.bob.best-models.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# Coefficients of the best model
cf = E$coef(best$model)
if(!is.null(cf)) write.table(cf, file.path(outdir, "results.bob.best-model-coefficients.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

say("Statistic: ", stat_name)
say("Adjusted p-values: Pr(chi^2_1 >= 2 log(T / (nu mu sqrt(xi)))) with nu = ", nu, "; they control the FWER over all tests (closed testing)")
say("Models scored: ", length(ls(E$cache, all.names=TRUE)))
say("Finished in ", format(round(difftime(Sys.time(), t0, units="secs"), 1)))
close(logcon)
cat("\nTop results (", file.path(outdir, "results.bob.variables.tsv"), "):\n", sep="")
show = tab[!is_group, ][order(-tab$`.ord`[!is_group]), c("Test", "log10 statistic", "-log10 Adjusted p-value", "Significant (FWER)", "Effect size in best model with", "Standard error")]
print(head(show, 20), row.names=FALSE)
