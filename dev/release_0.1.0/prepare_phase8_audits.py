from pathlib import Path
import hashlib
import json

original = Path('C:/Users/peterman.73/OneDrive - The Ohio State University/Research/R_packages/terradish/review_work_20260924/package_audit')
target = Path('dev/release_0.1.0/phase8_audits')
target.mkdir(exist_ok=True)
files = ['core/helper.R', 'core/exp1b_wishart_clustered.R', 'core/exp7_optim.R',
         'core/exp9_rho0.R', 'scale/shift_check.R', 'scale/ksr_units.R',
         'flex/common.R', 'flex/exp2_spline_support.R', 'compute/common.R',
         'compute/e1_graph.R', 'compute/e4_misc.R', 'compute/e6_melip_crop.R',
         'ext/check1_bruteforce.R', 'ext/check2_behaviour.R', 'ext/check5_kernel_diag.R',
         'cv/setup.R', 'cv/e5.R', 'cv/e6.R']
manifest = []
for name in files:
    src = original / name
    s = src.read_text(encoding='utf-8')
    original_hash = hashlib.sha256(src.read_bytes()).hexdigest()
    s = s.replace('/tmp/claude-0/audit/', target.resolve().as_posix() + '/')
    s = s.replace('source("common.R")', f'source("{(target / "flex/common.R").resolve().as_posix()}")')
    s = s.replace('source("setup.R")', f'source("{(target / "cv/setup.R").resolve().as_posix()}")')
    s = s.replace('terradish:::.design_loglinear_model', 'terradish:::.loglinear_conductance_from_matrix')
    if name == 'compute/e1_graph.R':
        s = s[:s.index('# helper:')]
        s += '\nstopifnot(isTRUE(all.equal(values(scale_to_0_1(r))[, 1], values(scale_to_0_1(r))[, 2])))\n'
    if name == 'compute/e4_misc.R':
        s = s.splitlines()[0] + '\n' + s[s.index('cat("==== E4b'):]
    if name == 'compute/e6_melip_crop.R':
        s = s[:s.index('cat("==== nested')]
        s = s.replace('sqrt(diag(solve(-ff$mle$hessian)))', 'sqrt(diag(vcov(ff)))')
        s = s.replace('sqrt(diag(solve(-fcr$mle$hessian)))', 'sqrt(diag(vcov(fcr)))')
    if name == 'ext/check1_bruteforce.R':
        a, b = s.index('cat("\\n==== (c)'), s.index('cat("\\n==== (d)')
        s = s[:a] + s[b:]
        s = s[:s.index('cat("\\n==== (e)')]
        s = s.replace('zc <- z - mean(z)', 'zc <- as.numeric(scale(z))\nH <- diag(n) - 1/n\nL <- cbind(diag(n - 1), -1)')
        s = s.replace('"kernel == centered zz\'/mean(zc^2):", all.equal(K, tcrossprod(zc) / mean(zc^2)', '"kernel == centered absdiff:", all.equal(K, -0.5 * H %*% abs(outer(zc, zc, "-")) %*% H')
        s = s.replace('-ldwish(nu * S, nu, Sig)', '-ldwish(nu * L %*% S %*% t(L), nu, L %*% Sig %*% t(L))')
        s = s.replace('cat("objective differences (pkg):', 'stopifnot(max(abs(diff(obj) - diff(bf))) < 1e-8)\ncat("objective differences (pkg):', 1)
        s = s.replace('cat("objective differences (pkg):   ", round(diff(objg)', 'stopifnot(max(abs(diff(objg) - diff(bfg))) < 1e-8)\ncat("objective differences (pkg):   ", round(diff(objg)')
        s = s.replace('lambda_env', 'lambda_absdiff_env')
        # Direct MLPE pair construction below uses unscaled site values.
        s = s.replace('gm <- mlpe_covariates(site, transform = "absdiff")', 'gm <- mlpe_covariates(site, transform = "absdiff", scale = FALSE)')
        s = s.replace('pairwise_endpoint_covariates(site, transform = "absdiff")', 'pairwise_endpoint_covariates(site, transform = "absdiff", scale = FALSE)')
        s += '\nstopifnot(abs(om + llb) < 1e-8)\n'
    if name == 'ext/check2_behaviour.R':
        s = s[:s.index('# kernel vs drift')]
        s = s.replace('scale` is a no-op when normalize = TRUE', 'scale changes absdiff kernel magnitude')
        s = s.replace('rank-1 kernel: covariance by position on the axis, not similarity', 'absdiff kernel uses centered pairwise differences')
        s = s.replace('(z_i-z_j)^2/mean(zc^2)', 'abs difference after site standardization')
        s = s.replace('outer(zz, zz, "-")^2 / mean(zz^2)', 'abs(outer(as.numeric(scale(zz)), as.numeric(scale(zz)), "-"))')
        s = s.replace('anova uses chi2_1, not chi-bar', 'anova now uses chi-bar-square')
        s = s.replace('anova accepts non-nested measurement models', 'anova must reject non-nested measurement models')
        s = s.replace('wishart_covariates(env)', 'wishart_covariates(env, model = "wishart_covariance")')
        s = s.replace('lambda_var1', 'lambda_absdiff_var1')
        s += '\ncat("Corrected README argument order:"); print(dim(pairwise_endpoint_covariates(cov, co20)))\n'
    if name == 'ext/check5_kernel_diag.R':
        s = s.replace('# Does the rank-1 kernel', '# Does the new default absdiff kernel')
        s = s.replace('lambda_env', 'lambda_absdiff_env')
        s = s.replace('model = mdl)', 'model = mdl, transform = "absdiff")')
    if name == 'cv/e5.R':
        s = s.replace('z[, "a", drop = FALSE]', 'as.data.frame(z[, "a", drop = FALSE])')
        s = s.replace('z[, c("b","c")]', 'as.data.frame(z[, c("b","c")])')
        s = s.replace('wishart_covariates(as.data.frame(z[, "a", drop = FALSE]))', 'wishart_covariates(as.data.frame(z[, "a", drop = FALSE]), model = "wishart_covariance")')
        s = s.replace('wishart_covariates(as.data.frame(z[, c("b","c")]))', 'wishart_covariates(as.data.frame(z[, c("b","c")]), model = "wishart_covariance")')
    if name == 'cv/e6.R':
        s = s.replace('ck <- file.path(getwd(), "ck.rds"); unlink(ck)', 'ck <- tempfile("phase8-cv-checkpoint-", fileext = ".rds")')
        a = s.index('b <- terradish_cv_folds')
        s = s[:a] + '''message <- tryCatch({
  terradish_cv_folds(surface, list(m = S ~ x1), folds = rf,
    model = wishart_covariance, nu = 25, checkpoint = ck)
  "ERROR: changed response was accepted"
}, error = conditionMessage)
cat("Changed-response resume:", message, "\\n")
stopifnot(grepl("incompatible|match|changed", message, ignore.case = TRUE))
unlink(ck)
})
'''
    dst = target / name
    dst.parent.mkdir(exist_ok=True)
    dst.write_text(s, encoding='utf-8')
    manifest.append(dict(source=str(src), source_sha256=original_hash,
                         adapted=str(dst), adapted_sha256=hashlib.sha256(dst.read_bytes()).hexdigest()))
(target / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
