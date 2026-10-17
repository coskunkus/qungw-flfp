###############################################################################
## 06_lr_intervals.R
## Likelihood ratio (LR) confidence intervals for the regression coefficients
## (Section 2.3 and Section 4.3, Fig. 2), computed by the formulation of
## Doganaksoy (2021, Technometrics 63:127-135): the lower (upper) 100(1-a)% LR
## limit of beta_j is the smallest (largest) value of beta_j satisfying
##     l(psi) = kappa,   kappa = l(psi_hat) - chi2_{1;1-a} / 2,
## with all the parameters treated as optimisation variables. Each limit is
## one constrained optimisation (Rsolnp::solnp). As a check, the limits are
## also obtained in the conventional way, as the roots of the profile
## log-likelihood, and the two sets of limits are compared.
## Requires qungw_core.R (best, TAU, X, y, make_loglik).
###############################################################################

cat("\n================ 06_lr_intervals.R ================\n")
if (!requireNamespace("Rsolnp", quietly = TRUE)) install.packages("Rsolnp", repos = "https://cloud.r-project.org")

ll_sum <- function(p) sum(make_loglik(TAU)(p))
p_hat  <- coef(best); l_hat <- ll_sum(p_hat)
kappa  <- l_hat - qchisq(0.95, 1) / 2
H      <- (best$hessian + t(best$hessian)) / 2; se_t <- sqrt(diag(solve(-H)))
jj     <- 3:length(p_hat)                                   # beta_0, beta_1, beta_2

## Doganaksoy (2021): min / max beta_j subject to l(psi) = kappa
lr_limit <- function(j, side) {                             # side = -1 lower, +1 upper
  st <- p_hat; st[j] <- p_hat[j] + side * 1.96 * se_t[j]    # start near the Wald limit
  f  <- Rsolnp::solnp(st, fun = function(p) -side * p[j],
                      eqfun = function(p) ll_sum(p), eqB = kappa,
                      control = list(trace = 0, tol = 1e-10, delta = 1e-8))
  c(limit = unname(f$pars[j]), constraint = ll_sum(f$pars) - kappa, conv = f$convergence)
}

## conventional check: profile log-likelihood l(beta_j, lambda_tilde)
profile_ll <- function(j, b, start = p_hat) {
  free <- setdiff(seq_along(p_hat), j)
  g <- function(q) { p <- start; p[free] <- q; p[j] <- b; -ll_sum(p) }
  o <- optim(start[free], g, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
  p <- start; p[free] <- o$par; p[j] <- b; list(l = -o$value, par = p)
}
profile_root <- function(j, side) {
  h <- function(b) profile_ll(j, b)$l - kappa
  ends <- p_hat[j] + side * c(0, 4) * se_t[j]
  uniroot(h, sort(ends), tol = 1e-10)$root
}

LR <- t(sapply(jj, function(j) {
  lo <- lr_limit(j, -1); hi <- lr_limit(j, +1)
  c(est = unname(p_hat[j]), wald_lo = unname(p_hat[j] - 1.96 * se_t[j]), wald_hi = unname(p_hat[j] + 1.96 * se_t[j]),
    lr_lo = lo[["limit"]], lr_hi = hi[["limit"]],
    prof_lo = profile_root(j, -1), prof_hi = profile_root(j, +1),
    max_constraint_gap = max(abs(c(lo[["constraint"]], hi[["constraint"]]))),
    solnp_conv = max(lo[["conv"]], hi[["conv"]]))
}))
rownames(LR) <- c("Intercept", "Fertility", "Log GDP per capita")
cat(sprintf("kappa = %.4f (l_hat = %.4f)\n", kappa, l_hat)); print(round(LR, 5))
cat(sprintf("Largest difference between the two LR computations: %.2e\n",
            max(abs(LR[, c("lr_lo", "lr_hi")] - LR[, c("prof_lo", "prof_hi")]))))

## profile curves for Fig. 2
PROF <- lapply(seq_along(jj), function(k) { j <- jj[k]
  b <- seq(LR[k, "lr_lo"] - 0.35 * (LR[k, "lr_hi"] - LR[k, "lr_lo"]),
           LR[k, "lr_hi"] + 0.35 * (LR[k, "lr_hi"] - LR[k, "lr_lo"]), length.out = 61)
  data.frame(b = b, l = sapply(b, function(v) profile_ll(j, v)$l)) })
RES$lr <- list(kappa = kappa, l_hat = l_hat, table = LR, profiles = PROF)
save_rds(RES$lr, file.path(OUTDIR, "lr_intervals.rds"))
