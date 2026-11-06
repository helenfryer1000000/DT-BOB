############################################################################
# doublethink-bob-functions.R
#
# Functions for the "best of both" (top-1) and top-m search versions of
# Doublethink. They sit beside the original MCMC code in code/ and use the
# same model, prior and likelihood:
#
#   log PO(s : 0) = |s| log(mu) + (df_s / 2) log(xi) + (1 - xi) * loglik(s : 0)
#
# where s is a model (a set of variables), df_s its number of coefficients,
# xi = h / (n + h), mu the prior odds of including each variable, and
# loglik(s : 0) the maximised log-likelihood of s minus that of the grand
# null model. This is calc.logBF() plus the prior of lr.mcmc() in
# code/doublethink-functions.R.
#
# See docs/METHODS.md for the statistics and docs/SEARCH.md for the search.
############################################################################

## Small helpers ###########################################################

logsumexp = function(x) {
	x = x[is.finite(x)]
	if(length(x) == 0) return(-Inf)
	mx = max(x)
	mx + log(sum(exp(x - mx)))
}

# A model is an increasing integer vector of variable indices.
# Its key is a string: "m" followed by the indices, so "m" is the grand null model.
model_key = function(s) paste0("m", paste(sort.int(s), collapse=","))
key2model = function(k) if(k == "m") integer(0) else as.integer(strsplit(substring(k, 2), ",", fixed=TRUE)[[1]])
norm_model = function(s) sort.int(unique(as.integer(s)))

# Doublethink asymptotic p-value for a statistic with log value logT.
# Pr(chi^2_1 >= 2 log(T / (k mu sqrt(xi)))), equal to 1 when the log is negative.
# k = number of variables tested (unadjusted) or nu (adjusted for the FWER).
doublethink_p = function(logT, k, mu, xi, log.p=FALSE) {
	x = logT - log(k) - log(mu) - 0.5 * log(xi)
	q = 2 * pmax(x, 0)
	pchisq(q, df=1, lower.tail=FALSE, log.p=log.p)
}

## Data #####################################################################

# Read the data exactly as code/run-mcmc.R does (load.data), without the
# side effects that are only needed by the MCMC chains.
bob_load_data = function(columns_filename, input_filename, cases_filename, exclude_filename,
                         idcol="eid", phenocol="pheno", controls_subsample=0, seed_subsample=1) {
	col = rbind(cbind("row", "character"), read.csv(columns_filename, header=FALSE, stringsAsFactors=FALSE), stringsAsFactors=FALSE)
	colnames(col) = c("column", "class")
	col$class[col$class == "float"] = "double"
	a = read.csv(input_filename, check.names=FALSE, colClasses=as.character(col$class))
	if(!is.null(exclude_filename) && !file.exists(exclude_filename)) stop("exclude_filename not found: ", exclude_filename)
	exclude.eids = if(!is.null(exclude_filename)) scan(exclude_filename, quiet=TRUE) else numeric(0)
	a = a[!a[, idcol] %in% exclude.eids, ]
	isf = col$class == "factor"
	if(any(isf)) {
		wh = suppressWarnings(which(apply(a[, isf, drop=FALSE], 2, function(y) as.numeric(as.character(y))) < 0, arr.ind=TRUE))
		a[, isf][wh] = "-999"
		stopifnot(sum(is.na(a[, isf, drop=FALSE])) == 0)
		for(j in which(isf)) a[, j] = droplevels(as.factor(a[, j]))
	}
	for(j in which(!isf)) {
		if(!is.numeric(a[, j])) next
		IS.NA = is.na(a[, j])
		if(sum(IS.NA) > 0) a[IS.NA, j] = mean(a[, j], na.rm=TRUE)
	}
	case.eids = read.csv(cases_filename)
	case.eids = case.eids[!case.eids[, idcol] %in% exclude.eids, ]
	pheno = data.frame(case.eids[match(a[, idcol], case.eids[, idcol]), ])
	if(controls_subsample != 0) {
		controls = pheno[pheno[, phenocol] == 0, ]
		cases = pheno[pheno[, phenocol] == 1, ]
		saved.seed = if(exists(".Random.seed", envir=.GlobalEnv)) get(".Random.seed", envir=.GlobalEnv) else NULL
		set.seed(seed_subsample)
		keep = c(cases[, idcol], sample(controls[, idcol], controls_subsample, replace=FALSE))
		if(!is.null(saved.seed)) assign(".Random.seed", saved.seed, envir=.GlobalEnv)
		a = a[a[, idcol] %in% keep, ]
		pheno = pheno[pheno[, idcol] %in% keep, ]
	}
	stopifnot(all(a[, idcol] == pheno[, idcol]))
	if(any(is.na(pheno[, phenocol]))) stop("Some IDs in input_filename have no phenotype in cases_filename")
	rownames(a) = a[, idcol]
	a[, idcol] = NULL
	a[, match("", colnames(a))] = NULL
	list(a=a, y=as.numeric(pheno[, phenocol]))
}

# Turn the data frame of variables into a numeric design matrix X plus a list
# 'blocks' giving the columns of X that belong to each variable. A numeric
# variable is one column; a factor with L levels is L - 1 treatment-contrast
# columns, as glm() would code it. Variables with one level, with more than
# nlevmax levels, or named in colnames.remove are left out, as in lr.mcmc().
bob_build_design = function(a, colnames.remove=character(0), nlevmax=50, intercept=TRUE) {
	nlev = unlist(lapply(a, nlevels))
	rm = match(colnames.remove, colnames(a)); rm = rm[!is.na(rm)]
	use = setdiff(which(nlev != 1 & nlev <= nlevmax), rm)
	if(length(use) == 0) stop("No variables left to analyse")
	if(!intercept && any(vapply(use, function(j) is.factor(a[, j]), TRUE)))
		stop("Factors need intercept=\"true\" for best of both, top-m and exhaustive (without an intercept glm codes the first factor in each model differently, so a factor's contribution would depend on the model)")
	pieces = list(); blocks = list(); start = 0
	for(j in use) {
		if(is.factor(a[, j])) {
			M = model.matrix(~ f, data.frame(f=a[, j]))[, -1, drop=FALSE]
			colnames(M) = paste0(colnames(a)[j], ":", levels(a[, j])[-1])
		} else {
			M = matrix(as.numeric(a[, j]), ncol=1, dimnames=list(NULL, colnames(a)[j]))
		}
		pieces[[length(pieces) + 1]] = M
		blocks[[length(blocks) + 1]] = start + seq_len(ncol(M))
		start = start + ncol(M)
	}
	X = do.call(cbind, pieces)
	names(blocks) = colnames(a)[use]
	list(X=X, blocks=blocks, varnames=colnames(a)[use], df=vapply(blocks, length, 0L))
}

## The scoring engine ########################################################
# make_engine() returns an environment with
#   lpo(s)          exact log PO(s : 0), cached
#   add_gain(s, J)  approximate gain in log PO from adding each variable in J to s
#                   (exact for the Gaussian family; score test for GLMs)
#   coef(s)         coefficient table of model s
make_engine = function(y, design, family="gaussian", intercept=TRUE, mu=0.1, h=1, chunk=50) {
	E = new.env()
	X = design$X; blocks = design$blocks; df = design$df
	n = length(y); nu = length(blocks)
	E$n = n; E$nu = nu; E$df = df; E$varnames = design$varnames
	E$mu = mu; E$h = h; E$xi = h / (n + h)
	E$family = family; E$intercept = intercept
	E$cache = new.env(hash=TRUE, size=100000L)
	E$n_fits = 0
	logmu = log(mu); logxi = log(E$xi); one_minus_xi = 1 - E$xi
	single = all(df == 1)
	col_of = if(single) unlist(blocks) else NULL

	if(family == "gaussian") {
		# Everything is a function of X'X, X'y and y'y: no pass over the data is needed.
		if(intercept) {
			Xc = sweep(X, 2, colMeans(X)); yc = y - mean(y)
		} else {
			Xc = X; yc = y
		}
		G = crossprod(Xc); g = drop(crossprod(Xc, yc)); RSS0 = sum(yc^2); dG = diag(G)
		rm(Xc)
		if(RSS0 <= 0) stop("The phenotype has no variation")
		fitg = function(cols) {
			if(length(cols) == 0) return(list(ok=TRUE, R=NULL, b=numeric(0), RSS=RSS0))
			Gs = G[cols, cols, drop=FALSE]
			R = tryCatch(chol(Gs), error=function(e) NULL)
			if(is.null(R) || min(diag(R))^2 < 1e-10 * max(diag(Gs))) return(list(ok=FALSE))
			b = backsolve(R, forwardsolve(t(R), g[cols]))
			RSS = RSS0 - sum(g[cols] * b)
			if(RSS <= 1e-12 * RSS0) return(list(ok=FALSE))
			list(ok=TRUE, R=R, b=b, RSS=RSS)
		}
		loglik = function(s) {
			f = fitg(unlist(blocks[s], use.names=FALSE))
			if(!f$ok) return(-Inf)
			-0.5 * n * log(f$RSS / RSS0)
		}
		E$add_gain = function(s, J) {
			if(length(J) == 0) return(numeric(0))
			cols = unlist(blocks[s], use.names=FALSE)
			f = fitg(cols)
			if(!f$ok) return(rep(-Inf, length(J)))
			out = numeric(length(J))
			if(all(df[J] == 1)) {
				jc = unlist(blocks[J], use.names=FALSE)
				if(length(cols)) {
					W = forwardsolve(t(f$R), G[cols, jc, drop=FALSE])
					r = g[jc] - drop(crossprod(G[cols, jc, drop=FALSE], f$b))
					v = dG[jc] - colSums(W^2)
				} else {
					r = g[jc]; v = dG[jc]
				}
				frac = ifelse(v > 1e-10 * dG[jc], r^2 / (v * f$RSS), NA)
				dll = ifelse(is.na(frac) | frac >= 1, -Inf, -0.5 * n * log1p(-pmin(frac, 1 - 1e-15)))
				out = logmu + 0.5 * logxi + one_minus_xi * dll
			} else {
				for(i in seq_along(J)) {
					jc = blocks[[J[i]]]
					if(length(cols)) {
						W = forwardsolve(t(f$R), G[cols, jc, drop=FALSE])
						r = g[jc] - drop(crossprod(G[cols, jc, drop=FALSE], f$b))
						V = G[jc, jc, drop=FALSE] - crossprod(W)
					} else {
						r = g[jc]; V = G[jc, jc, drop=FALSE]
					}
					q = tryCatch(sum(r * solve(V, r)), error=function(e) NA)
					frac = q / f$RSS
					dll = if(is.na(frac) || frac >= 1) -Inf else -0.5 * n * log1p(-frac)
					out[i] = logmu + 0.5 * df[J[i]] * logxi + one_minus_xi * dll
				}
			}
			out
		}
		E$coef = function(s) {
			cols = unlist(blocks[s], use.names=FALSE)
			f = fitg(cols)
			if(!f$ok || length(cols) == 0) return(NULL)
			sigma2 = f$RSS / (n - length(cols) - intercept)
			Ri = backsolve(f$R, diag(length(cols)))
			se = sqrt(sigma2 * rowSums(Ri^2))
			data.frame(variable=rep(names(blocks)[s], df[s]), term=colnames(X)[cols], estimate=f$b, se=se, stringsAsFactors=FALSE)
		}
	} else if(family %in% c("binomial", "poisson")) {
		fam = get(family)()
		Z0 = if(intercept) matrix(1, n, 1) else matrix(0, n, 0)
		fitglm = function(cols) {
			Z = cbind(Z0, X[, cols, drop=FALSE])
			E$n_fits = E$n_fits + 1
			if(ncol(Z) == 0) {
				eta = rep(0, n); m = fam$linkinv(eta)
				return(list(ok=TRUE, dev=sum(fam$dev.resids(y, m, rep(1, n))), mu=m, w=fam$mu.eta(eta)^2 / fam$variance(m), Z=Z, coef=numeric(0), rank=0))
			}
			fit = suppressWarnings(glm.fit(Z, y, family=fam, intercept=intercept))
			if(fit$rank < ncol(Z) || any(is.na(fit$coefficients))) return(list(ok=FALSE))
			list(ok=TRUE, dev=fit$deviance, mu=fit$fitted.values, w=fit$weights, Z=Z, coef=fit$coefficients, fit=fit)
		}
		dev0 = fitglm(integer(0))$dev
		# Fitted means and weights of recently fitted models, so that screening
		# candidates against a model just scored needs no refit.
		recent = new.env(hash=TRUE); recent_keys = character(0); n_recent = max(20, min(1000, floor(1.25e7 / n)))  # about 200 MB
		remember = function(key, f) {
			assign(key, list(mu=f$mu, w=f$w), envir=recent)
			recent_keys <<- c(recent_keys, key)
			if(length(recent_keys) > n_recent) { rm(list=recent_keys[1], envir=recent); recent_keys <<- recent_keys[-1] }
		}
		loglik = function(s) {
			f = fitglm(unlist(blocks[s], use.names=FALSE))
			if(!f$ok) return(-Inf)
			remember(model_key(s), f)
			0.5 * (dev0 - f$dev)
		}
		# Rao score test for adding each variable in J to s. Half the score
		# statistic approximates the gain in log-likelihood.
		# The fit of the model last used for screening is kept, because the
		# search screens several candidate sets against the same model.
		last = list(key=NULL)
		E$add_gain = function(s, J) {
			if(length(J) == 0) return(numeric(0))
			key = model_key(s)
			if(!identical(last$key, key)) {
				f = recent[[key]]
				if(is.null(f)) {
					f = fitglm(unlist(blocks[s], use.names=FALSE))
					if(!f$ok) return(rep(-Inf, length(J)))
				}
				sw = sqrt(f$w)
				Z = cbind(Z0, X[, unlist(blocks[s], use.names=FALSE), drop=FALSE])
				Q = if(ncol(Z)) qr.Q(qr(Z * sw)) else matrix(0, n, 0)
				last <<- list(key=key, res=y - f$mu, sw=sw, Q=Q)
			}
			res = last$res; sw = last$sw; Q = last$Q
			out = numeric(length(J))
			one = df[J] == 1
			# Single-column variables: vectorised over chunks of columns
			Js = which(one)
			for(ii in split(Js, ceiling(seq_along(Js) / chunk))) {
				M = X[, unlist(blocks[J[ii]], use.names=FALSE), drop=FALSE]
				U = colSums(M * res)
				WM = M * sw
				if(ncol(Q)) WM = WM - Q %*% crossprod(Q, WM)
				V = colSums(WM^2)
				S = ifelse(V > 1e-10 * colSums((M * sw)^2), U^2 / V, NA)
				out[ii] = ifelse(is.na(S), -Inf, logmu + 0.5 * logxi + one_minus_xi * 0.5 * S)
			}
			# Factors: one block at a time
			for(i in which(!one)) {
				M = X[, blocks[[J[i]]], drop=FALSE]
				U = drop(crossprod(M, res))
				WM = M * sw
				if(ncol(Q)) WM = WM - Q %*% crossprod(Q, WM)
				S = tryCatch(sum(U * solve(crossprod(WM), U)), error=function(e) NA)
				out[i] = if(is.na(S)) -Inf else logmu + 0.5 * df[J[i]] * logxi + one_minus_xi * 0.5 * S
			}
			out
		}
		E$coef = function(s) {
			cols = unlist(blocks[s], use.names=FALSE)
			f = fitglm(cols)
			if(!f$ok || length(cols) == 0) return(NULL)
			tab = coef(summary.glm(structure(f$fit, class=c("glm", "lm"))))
			tab = tab[nrow(tab) - rev(seq_along(cols)) + 1, , drop=FALSE]
			data.frame(variable=rep(names(blocks)[s], df[s]), term=colnames(X)[cols], estimate=tab[, 1], se=tab[, 2], stringsAsFactors=FALSE)
		}
	} else stop("glm_family must be gaussian, binomial or poisson")

	E$lpo = function(s) {
		k = model_key(s)
		v = E$cache[[k]]
		if(!is.null(v)) return(v)
		v = length(s) * logmu + 0.5 * sum(df[s]) * logxi + one_minus_xi * loglik(s)
		assign(k, v, envir=E$cache)
		v
	}
	E
}

## Constraints ##############################################################
# A "side" is the set of models a search may visit:
#   forbid      variables that may not be in the model
#   require_any variables of which at least one must be in the model
side_ok = function(s, side) {
	if(length(side$forbid) && any(s %in% side$forbid)) return(FALSE)
	if(length(side$require_any) && !any(s %in% side$require_any)) return(FALSE)
	TRUE
}
SIDE_FREE = list(forbid=integer(0), require_any=integer(0))

## Neighbourhood and local search ############################################
# Neighbours of s on a side: add one variable, drop one variable, or swap one
# variable for another. Adds and swaps are screened by add_gain(): only the
# ncand most promising candidates are scored exactly. Returns a list of models.
neighbours = function(E, s, side, ncand=10, max_size=Inf) {
	out = list()
	cand = setdiff(seq_len(E$nu), c(s, side$forbid))
	top_adds = function(base, cand) {
		if(length(cand) == 0) return(integer(0))
		gn = E$add_gain(base, cand)
		ok = is.finite(gn)
		cand = cand[ok]; gn = gn[ok]
		cand[order(-gn)][seq_len(min(ncand, length(cand)))]
	}
	if(length(s) < max_size) for(j in top_adds(s, cand)) {
		t = norm_model(c(s, j)); if(side_ok(t, side)) out[[length(out) + 1]] = t
	}
	# Swaps: screen every (drop i, add j) pair by log PO(s - i) + gain(j),
	# then score exactly the ncand best pairs and the best pair for each i.
	sw = list()
	for(i in s) {
		d = s[s != i]
		if(side_ok(d, side)) out[[length(out) + 1]] = d
		cand_i = cand
		if(length(side$require_any) && !any(d %in% side$require_any)) cand_i = intersect(cand, side$require_any)
		if(length(cand_i) == 0) next
		gn = E$add_gain(d, cand_i)
		ok = is.finite(gn)
		if(!any(ok)) next
		sw[[length(sw) + 1]] = data.frame(i=i, j=cand_i[ok], score=E$lpo(d) + gn[ok])
	}
	if(length(sw)) {
		sw = do.call(rbind, sw)
		pick = unique(c(order(-sw$score)[seq_len(min(ncand, nrow(sw)))], vapply(split(seq_len(nrow(sw)), sw$i), function(r) r[which.max(sw$score[r])], 0L)))
		for(r in pick) {
			t = norm_model(c(s[s != sw$i[r]], sw$j[r])); if(side_ok(t, side)) out[[length(out) + 1]] = t
		}
	}
	out
}

# Steepest-ascent hill climb on log PO from s, staying on the side.
# Every model scored is recorded in 'seen' (an environment key -> log PO).
climb = function(E, s, side, ncand=10, max_size=Inf, seen=NULL, max_steps=1000) {
	s = norm_model(s)
	# A start larger than max_size is first cut down by the best drops.
	while(length(s) > max_size) {
		drops = lapply(s, function(i) s[s != i])
		drops = drops[vapply(drops, side_ok, TRUE, side=side)]
		if(length(drops) == 0) break
		s = drops[[which.max(vapply(drops, E$lpo, 0))]]
	}
	cur = E$lpo(s)
	if(!is.null(seen)) assign(model_key(s), cur, envir=seen)
	for(step in seq_len(max_steps)) {
		nb = neighbours(E, s, side, ncand, max_size)
		if(length(nb) == 0) break
		vals = vapply(nb, E$lpo, 0)
		if(!is.null(seen)) for(i in seq_along(nb)) assign(model_key(nb[[i]]), vals[i], envir=seen)
		i = which.max(vals)
		if(!(vals[i] > cur + 1e-10)) break
		s = nb[[i]]; cur = vals[i]
	}
	list(model=s, lpo=cur)
}

# A random add or drop move, used to perturb a model between climbs.
perturb = function(s, nu, side, nmoves=3, max_size=Inf) {
	for(r in seq_len(nmoves)) {
		pool = setdiff(seq_len(nu), side$forbid)
		if(length(s) > 0 && (runif(1) < 0.5 || length(s) >= max_size)) {
			t = s[-sample.int(length(s), 1)]
		} else {
			add = setdiff(pool, s)
			if(length(add) == 0) next
			t = norm_model(c(s, add[sample.int(length(add), 1)]))
		}
		if(side_ok(t, side)) s = t
	}
	s
}

## Starting models from the lasso #############################################
# The lasso (L1-penalised GLM, package glmnet) is a fast machine-learning
# screen: the supports along its regularisation path are good starting points
# for the L0 search on log PO. Returns a list of distinct models.
lasso_starts = function(y, design, family, intercept, nstarts=5, E=NULL) {
	if(!requireNamespace("glmnet", quietly=TRUE)) return(list())
	fit = tryCatch(glmnet::glmnet(design$X, y, family=family, intercept=intercept, nlambda=100, dfmax=min(ncol(design$X), 200)), error=function(e) NULL)
	if(is.null(fit)) return(list())
	B = as.matrix(fit$beta != 0)
	col2var = rep(seq_along(design$blocks), design$df)
	models = unique(lapply(seq_len(ncol(B)), function(l) norm_model(col2var[B[, l]])))
	if(!is.null(E)) {
		# Models grow along the path and log PO usually rises then falls, so
		# stop scoring after 10 models in a row that do not improve on the best.
		v = rep(-Inf, length(models)); bestv = -Inf; since = 0
		for(i in seq_along(models)) {
			v[i] = E$lpo(models[[i]])
			if(v[i] > bestv) { bestv = v[i]; since = 0 } else since = since + 1
			if(since >= 10) break
		}
		keep = is.finite(v)
		models = models[keep][order(-v[keep])]
	}
	models[seq_len(min(nstarts, length(models)))]
}

## Global search for the best model ###########################################
# Climbs from the empty model, from the lasso starts and, by iterated local
# search, from random perturbations of the best model found so far.
find_best_model = function(E, y, design, ncand=10, n_restarts=20, perturb_moves=3, use_lasso=TRUE, nlasso=5, max_size=Inf, verbose=TRUE) {
	starts = list(integer(0))
	if(use_lasso) {
		ls = lasso_starts(y, design, E$family, E$intercept, nstarts=nlasso, E=E)
		if(verbose) cat("Lasso path gave", length(ls), "starting models\n")
		starts = c(starts, ls)
	}
	best = list(model=integer(0), lpo=E$lpo(integer(0)))
	for(s0 in unique(starts)) {
		r = climb(E, s0, SIDE_FREE, ncand, max_size)
		if(r$lpo > best$lpo) best = r
	}
	if(verbose) cat("Best model after deterministic starts: log10 PO =", signif(best$lpo / log(10), 6), "\n")
	for(r in seq_len(n_restarts)) {
		s0 = perturb(best$model, E$nu, SIDE_FREE, perturb_moves, max_size)
		res = climb(E, s0, SIDE_FREE, ncand, max_size)
		if(res$lpo > best$lpo + 1e-10) {
			best = res
			if(verbose) cat("  restart", r, "improved: log10 PO =", signif(best$lpo / log(10), 6), "\n")
		}
	}
	best
}

# The 'n' best models held in the engine's cache.
cache_top = function(E, n=50, side=SIDE_FREE) {
	keys = ls(E$cache, all.names=TRUE)
	vals = unlist(mget(keys, envir=E$cache))
	o = order(-vals)
	out = list(); k = 0
	for(i in o) {
		s = key2model(keys[i])
		if(side_ok(s, side)) { k = k + 1; out[[k]] = list(model=s, lpo=vals[i]) }
		if(k >= n) break
	}
	out
}

## One test ###################################################################
# Tests the null hypothesis that every variable in 'tested' has a zero effect.
# Side A: models with at least one tested variable. Side O: models with none.
# Finds the m best models on each side and returns
#   log T_m = log(sum of the m best PO on side A) - log(sum of the m best PO on side O).
# With m = 1 this is best of both: the best model with a tested variable
# against the best model without.
run_test = function(E, tested, best, elite=list(), m=1, ncand=10, max_size=Inf, enrich_rounds=5) {
	tested = norm_model(tested)
	sideA = list(forbid=integer(0), require_any=tested)
	sideO = list(forbid=tested, require_any=integer(0))
	seen = new.env(hash=TRUE)
	# Starting points: the best model projected onto each side, and the
	# best elite model on each side.
	startsA = if(any(best$model %in% tested)) list(best$model) else lapply(tested, function(t) norm_model(c(best$model, t)))
	if(length(startsA) > 10) {
		gn = E$add_gain(best$model, tested)
		startsA = startsA[order(-gn)[1:10]]
	}
	startsO = list(setdiff(best$model, tested))
	for(e in elite) {
		if(side_ok(e$model, sideA)) { startsA[[length(startsA) + 1]] = e$model; break }
	}
	for(e in elite) {
		if(side_ok(e$model, sideO)) { startsO[[length(startsO) + 1]] = e$model; break }
	}
	for(e in elite) assign(model_key(e$model), e$lpo, envir=seen)
	for(s in unique(startsA)) climb(E, s, sideA, ncand, max_size, seen)
	for(s in unique(startsO)) climb(E, s, sideO, ncand, max_size, seen)

	top_side = function(side, mm) {
		keys = ls(seen, all.names=TRUE)
		vals = unlist(mget(keys, envir=seen))
		ok = vapply(keys, function(k) side_ok(key2model(k), side), TRUE)
		keys = keys[ok]; vals = vals[ok]
		o = order(-vals)[seq_len(min(mm, sum(ok)))]
		list(keys=keys[o], lpo=vals[o])
	}
	# For m > 1, score the neighbourhood of each of the current m best models
	# on each side until the list of m best stops changing.
	if(m > 1) for(side in list(sideA, sideO)) {
		prev = ""
		for(r in seq_len(enrich_rounds)) {
			tp = top_side(side, m)
			sig = paste(tp$keys, collapse="|")
			if(sig == prev) break
			prev = sig
			for(k in tp$keys) {
				nb = neighbours(E, key2model(k), side, ncand, max_size)
				for(t in nb) assign(model_key(t), E$lpo(t), envir=seen)
			}
		}
	}
	tA = top_side(sideA, m); tO = top_side(sideO, m)
	list(tested=tested,
		logT=logsumexp(tA$lpo) - logsumexp(tO$lpo),
		topA=tA, topO=tO,
		n_seen=length(ls(seen, all.names=TRUE)),
		improved=max(c(tA$lpo[1], tO$lpo[1])) > best$lpo + 1e-8)
}

## Exhaustive enumeration (small nu only) ######################################
# Scores all 2^nu models. Gives the exact Doublethink statistic PO_{A:O}
# (full model averaging), and exact top-1 and top-m statistics. Used for
# checking the search; feasible up to about nu = 20.
exhaustive_scores = function(E, max_nu=20) {
	nu = E$nu
	if(nu > max_nu) stop("Exhaustive enumeration is limited to nu <= ", max_nu)
	M = as.matrix(expand.grid(rep(list(0:1), nu)))
	lp = apply(M, 1, function(r) E$lpo(which(r == 1)))
	list(M=M, lpo=lp)
}
exhaustive_test = function(ex, tested, m=1) {
	inA = rowSums(ex$M[, tested, drop=FALSE]) > 0
	a = sort(ex$lpo[inA], decreasing=TRUE); o = sort(ex$lpo[!inA], decreasing=TRUE)
	list(logPO=logsumexp(a) - logsumexp(o),
		logT=logsumexp(a[seq_len(min(m, length(a)))]) - logsumexp(o[seq_len(min(m, length(o)))]),
		logT1=a[1] - o[1])
}
