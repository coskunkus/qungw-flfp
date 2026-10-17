###############################################################################
##  qungw_core.R  -  shared core of the QUNGW unit quantile regression
##  ---------------------------------------------------------------------------
##  Paper : "A New Unit Quantile Regression Based on the Unit NG-Weibull
##           Distribution, with an Application to Female Labour Force
##           Participation" (Computational Economics, submitted)
##
##  Distribution functions of the UNGW / QUNGW family, the admissibility
##  condition gamma(-log mu)^alpha < B(u) of Remark 1, the penalised
##  log-likelihood, the multi-start maxLik protocol of Section 2.4, and the
##  beta, Kumaraswamy and unit-Weibull competitors. Every other script sources
##  this file, so the application, the diagnostics and the simulation use one
##  implementation. Nothing is written to disk here.
###############################################################################

## ---------------------------------------------------------------------------
## 0. Packages
## ---------------------------------------------------------------------------
pkgs <- c("maxLik", "numDeriv", "lamW")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
invisible(lapply(pkgs, library, character.only = TRUE))

## Every random stream is seeded with the generator named explicitly, so that a
## different RNGkind() in the user's session cannot change the results.
seed_mt <- function(s) set.seed(s, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
seed_mt(20260815)

## Set this to the same folder used by get_data.R. Forward slashes on Windows.
## folders: this file lives in R/, the data in ../data, outputs go to ../results
if (!exists("ROOT")) ROOT <- normalizePath(file.path(getwd(), ".."))
source("io_utils.R")   # .rds and figure I/O safe for non-ASCII paths
WORKDIR <- file.path(ROOT, "data")
OUTDIR  <- file.path(ROOT, "results"); dir.create(OUTDIR, showWarnings = FALSE)
## Do NOT setwd() here. These scripts used to switch the working directory to
## WORKDIR, which meant that a later source("other_script.R") with a relative
## path resolved against C:/ck instead of the code folder and silently ran a
## stale copy. Every path below is absolute, built from WORKDIR, so the working
## directory is left exactly as the caller set it.

## BASELINE QUANTILE.
## u = 0.25 rather than the median: by Remark 1 the lambda > 0 restriction is
## vacuous for u <= exp(-1) = 0.3679, so the optimum is interior, the standard
## errors are valid, and the comparison against the unconstrained competitors
## is like-for-like. At u = 0.5 on this sample the largest log-likelihood lies on
## the boundary of the admissible region (Table 4) and none of that holds.
TAU <- 0.25           # the quantile modelled (u in the paper)
OUT <- OUTDIR         # tables and figures are written here

###############################################################################
##  REPORT
##  Everything printed below goes to the console AND to report_qungw.txt.
##  split = TRUE keeps the console live so progress is still visible.
###############################################################################

cat("=============================================================\n")
cat("QUNGW UNIT QUANTILE REGRESSION : FEMALE LABOUR FORCE PARTICIPATION\n")
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Quantile modelled (u): ", TAU, "\n")
cat("=============================================================\n")

###############################################################################
## 1. DATA
###############################################################################
##
##  This script does NOT touch the network. Run get_data.R once to produce
##  flfp_data.csv, then run this script as often as you like.

DATAFILE <- file.path(WORKDIR, "flfp_data.csv")
if (!file.exists(DATAFILE))
  stop("Data file not found: ", DATAFILE,
       "\n  Run get_data.R first.", call. = FALSE)

dat_all <- read.csv(DATAFILE, stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## COVARIATE SET
## ---------------------------------------------------------------------------
## The full seven-covariate set is badly collinear in this cross-section:
## VIFs of 6-7, cor(lgdppc, internet) = 0.88, cor(fert, dep) = 0.83. With
## n = 119 that leaves no covariate significant and flips the signs of educ_f,
## dep and internet. The three retained below are the least mutually
## correlated representatives of the three channels the literature identifies:
##   fert    demographic burden          (Bloom et al., 2009)
##   lgdppc  level of development        (Goldin, 1995; feminisation U)
##   serv_f  structural transformation   (Olivetti, 2013)
## Set XVARS_FULL <- TRUE to reproduce the seven-covariate specification as a
## robustness check.
## The two-covariate specification gives k = 5 for QUNGW and k = 4 for the
## beta, Kumaraswamy and unit-Weibull competitors (Section 4.1).
XVARS_FULL <- FALSE
xvars <- if (XVARS_FULL)
  c("fert", "lgdppc", "urb", "dep", "serv_f", "internet", "educ_f") else
  c("fert", "lgdppc")

## Listwise deletion is applied to the response and to the covariates actually
## used, and to nothing else. Deleting on unused columns throws away countries
## for no reason: on this data the seven-covariate set leaves 119 economies
## while the three-covariate set leaves considerably more, and the countries
## recovered are disproportionately those with low participation.
dat <- dat_all[complete.cases(dat_all[, c("y", xvars)]), ]
if (exists("EXCLUDE") && length(EXCLUDE)) dat <- dat[!dat$iso3c %in% EXCLUDE, ]

stopifnot(all(dat$y > 0 & dat$y < 1))      # no boundary observations
cat(sprintf("\n%d economies with an observed response; %d complete on %s\n",
            nrow(dat_all), nrow(dat), paste(xvars, collapse = " + ")))
cat(sprintf("n = %d countries, y in [%.4f, %.4f]\n",
            nrow(dat), min(dat$y), max(dat$y)))

## ---------------------------------------------------------------------------
## DATA AUDIT. Run this first and read it before trusting anything downstream.
## ---------------------------------------------------------------------------
cat("\n================ DATA AUDIT ================\n")
cat("Rows:", nrow(dat), "  Columns:", ncol(dat), "\n\n")
cat("Column types:\n"); print(sapply(dat, class))
cat("\nNA count per column:\n"); print(colSums(is.na(dat)))
cat("\nFirst 6 rows:\n"); print(head(dat))
cat("\nRange of every numeric column:\n")
print(round(t(sapply(dat[sapply(dat, is.numeric)],
                     function(v) c(min = min(v), max = max(v)))), 4))
cat("\nSanity checks:\n")
cat("  y strictly inside (0,1):        ", all(dat$y > 0 & dat$y < 1), "\n")
cat("  shares urb/dep/serv_f/internet in [0,1.5]: ",
    all(sapply(dat[, c("urb","dep","serv_f","internet")],
               function(v) all(v >= 0 & v <= 1.5))), "\n")
cat("  lgdppc plausible (5 to 13):     ",
    all(dat$lgdppc > 5 & dat$lgdppc < 13), "\n")
cat("  fertility plausible (0.5 to 8): ",
    all(dat$fert > 0.5 & dat$fert < 8), "\n")
cat("  duplicated countries:           ", anyDuplicated(dat$iso3c) > 0, "\n")
cat("\nLowest and highest 5 by y:\n")
o <- order(dat$y)
print(dat[c(head(o, 5), tail(o, 5)), c("country", "iso3c", "y")])
cat("============================================\n\n")

## Covariates are standardised so that the shape parameters are well scaled.
## Coefficients are on the standardised scale; divide by the reported sd to
## return to natural units.
Xraw  <- as.matrix(dat[, xvars])
X     <- cbind(1, scale(Xraw))
colnames(X) <- c("(Intercept)", xvars)
y     <- dat$y
n     <- length(y)
p     <- ncol(X)

## ---------------------------------------------------------------------------
## Collinearity check. A singular or near-singular Hessian in the optimiser is
## almost always caused here, not by the likelihood: if two covariates are
## nearly collinear the log-likelihood is flat in one direction and Newton-type
## steps fail with "singular matrix 'a' in solve".
## ---------------------------------------------------------------------------
kappa_X <- kappa(X, exact = TRUE)
vif <- sapply(xvars, function(v) {
  others <- setdiff(xvars, v)
  r2 <- summary(lm(Xraw[, v] ~ Xraw[, others]))$r.squared
  1 / (1 - r2)
})
cat("\nCondition number of the design matrix:", round(kappa_X, 1), "\n")
cat("Variance inflation factors:\n"); print(round(vif, 2))
if (kappa_X > 100 || max(vif) > 10)
  warning("Design matrix is ill-conditioned (kappa > 100 or VIF > 10). ",
          "Drop or combine the offending covariates before trusting the ",
          "standard errors.")
cat("\nPairwise correlations:\n")
print(round(cor(Xraw), 2))

## Table 1: descriptive statistics
desc <- t(apply(cbind(y = y, Xraw), 2, function(v)
  c(Mean = mean(v), SD = sd(v), Min = min(v),
    Median = median(v), Max = max(v))))
print(round(desc, 3))

###############################################################################
## 2. THE QUNGW DISTRIBUTION
###############################################################################
##
##  PARAMETER SPACE. After the quantile re-parameterisation, lambda is no
##  longer a free parameter but the function xi(mu, alpha, gamma) of Eq. (8):
##
##      theta_i = xi_i = (1/gamma) * w_i^(-1) * ( exp(-gamma*w_i) - log(u) - 1 )
##      with w_i = (-log mu_i)^alpha
##
##  The baseline NG-Weibull requires lambda > 0, and this is binding rather
##  than cosmetic: for lambda < 0 the density takes negative values near y = 0
##  and the cdf is not monotone (verified numerically). Since lambda > 0
##  requires exp(-gamma*w_i) > 1 + log(u), the constraint is empty when
##  u <= exp(-1) = 0.3679 and binds otherwise. At the median, u = 0.5, it
##  requires
##
##      gamma * (-log mu_i)^alpha < 1.1814   for every i,
##
##  which is a joint restriction on (alpha, gamma) and on the fitted quantiles.
##  It is not a restriction on the data: gamma can always shrink to accommodate
##  a small mu_i. In practice it means the optimiser must be started from small
##  gamma when the sample contains low values of the response, which is why the
##  starting grid below reaches down to gamma = 0.02. The likelihood applies a
##  smooth penalty rather than returning -Inf, so that a search step into the
##  infeasible region can be reversed instead of stalling.

## xi(mu, alpha, gamma), Eq. (8). No Lambert W function is needed here.
xi_fun <- function(mu, alpha, gamma, u) {
  lm  <- -log(mu)                          # positive because 0 < mu < 1
  lma <- lm^alpha
  ## expm1 keeps the expression accurate as gamma -> 0 (Proposition 2)
  (expm1(-gamma * lma) - log(u)) / (gamma * lma)
}

## cdf, Eq. (9)
## written as log F = expm1(-gamma w) - (w/v) {expm1(-gamma v) - log u}, with
## w = (-log y)^alpha and v = (-log mu)^alpha, which is stable as gamma -> 0,
## where F tends to the unit-Weibull c.d.f. u^(w/v) (Proposition 2)
pQUNGW <- function(q, alpha, gamma, mu, u = TAU) {
  w <- (-log(q))^alpha; v <- (-log(mu))^alpha
  exp(expm1(-gamma * w) - (w / v) * (expm1(-gamma * v) - log(u)))
}

## pdf, Eq. (10)
dQUNGW <- function(x, alpha, gamma, mu, u = TAU, log = FALSE) {
  gxi <- gamma * xi_fun(mu, alpha, gamma, u)
  ly <- -log(x); w <- ly^alpha; eg <- exp(-gamma * w)
  ld <- log(alpha) + (alpha - 1) * log(ly) - log(x) +
        expm1(-gamma * w) - gxi * w + log(gxi + gamma * eg)
  if (log) ld else exp(ld)
}

## quantile function, Eqs. (5)-(6), through lamW::lambertW0, with a numerical
## inversion of the cdf as fallback.
qQUNGW <- function(prob, alpha, gamma, mu, u = TAU) {
  ## closed-form inversion through the Lambert W function, Eqs. (5)-(6):
  ## Q = exp{-(xi(lambda)/(lambda*gamma))^(1/alpha)} with lambda = xi(mu, alpha, gamma)
  mu  <- rep_len(mu, length(prob))
  lam <- xi_fun(mu, alpha, gamma, u)
  xl  <- lam * lamW::lambertW0(exp((1 + log(prob)) / lam) / lam) - 1 - log(prob)
  q   <- exp(-(xl / (lam * gamma))^(1 / alpha))
  bad <- !is.finite(q)                      # numerical fall-back: invert the c.d.f.
  if (any(bad)) q[bad] <- vapply(which(bad), function(i)
    uniroot(function(t) pQUNGW(t, alpha, gamma, mu[i], u) - prob[i],
            c(.Machine$double.xmin, 1 - .Machine$double.eps), tol = 1e-12)$root, numeric(1))
  q
}

rQUNGW <- function(alpha, gamma, mu, u = TAU)
  qQUNGW(runif(length(mu)), alpha, gamma, mu, u)

## largest gamma*(-log mu)^alpha compatible with lambda > 0 at a given u
feasible_bound <- function(u) if (1 + log(u) <= 0) Inf else -log(1 + log(u))

## SMALLEST ATTAINABLE CONDITIONAL QUANTILE.
## Letting lambda -> 0+ in Eq. (8) gives gamma*(-log mu)^alpha -> B(u), so for
## fixed (alpha, gamma) the model can only place the uth conditional quantile
## in the interval ( mu_min , 1 ) with
##      mu_min = exp( -( B(u)/gamma )^(1/alpha) ),  B(u) = -log(1+log u).
## If mu_min exceeds the smallest observed response, the model is structurally
## unable to fit the lower tail, whatever the covariates do. This is checked
## and reported below.
mu_min_attainable <- function(alpha, gamma, u) {
  B <- feasible_bound(u)
  if (!is.finite(B)) return(0)
  exp(-(B / gamma)^(1 / alpha))
}

###############################################################################
## 3. LOG-LIKELIHOOD FOR maxLik
##    par = ( log(alpha), log(gamma), beta_0, ..., beta_p )
##    Returning a vector of per-observation contributions lets maxLik use BHHH
##    and lets sandwich standard errors be computed if wanted.
###############################################################################

## The likelihood builders take the design matrix and the response as
## arguments, defaulting to the full sample. Nothing writes to the global
## environment: the case-deletion loop below simply passes a subset. An
## earlier version mutated X, y and n with `<<-`, which fails outright when
## the script is sourced into a locked environment ("cannot change value of
## locked binding for 'n'") and is fragile even when it works.
make_loglik <- function(u, Xm = X, yv = y) function(par) {
  nn <- length(yv)
  alpha <- exp(par[1]); gamma <- exp(par[2]); beta <- par[-(1:2)]

  ## keep the search inside a numerically safe region. Very large alpha
  ## overflows (-log y)^alpha. As gamma -> 0 the model tends to the unit-Weibull
  ## quantile regression (Proposition 2); the expm1 form below evaluates the
  ## log-likelihood accurately in that limit.
  if (!is.finite(alpha) || alpha < 1e-3 || alpha > 50)  return(rep(-1e10, nn))
  if (!is.finite(gamma) || gamma < 1e-6 || gamma > 1e3) return(rep(-1e10, nn))

  eta <- as.vector(Xm %*% beta)
  mu  <- 1 / (1 + exp(-eta))                 # logit link, Eq. (13)

  if (!all(is.finite(mu)) || any(mu <= 1e-10) || any(mu >= 1 - 1e-10))
    return(rep(-1e10, nn))

  xi <- xi_fun(mu, alpha, gamma, u)
  ly <- -log(yv); w <- ly^alpha; eg <- exp(-gamma * w)
  if (any(!is.finite(w)) || any(!is.finite(eg))) return(rep(-1e10, nn))

  ## smooth penalty on the infeasible region instead of a hard wall
  bad <- !is.finite(xi) | xi <= 0
  if (any(bad)) {
    viol <- sum(pmax(0, -xi[is.finite(xi)])) + sum(!is.finite(xi))
    return(rep(-1e8 * (1 + viol) / nn, nn))
  }
  if (any(xi + eg <= 0)) return(rep(-1e10, nn))

  ## gamma * xi = {expm1(-gamma v) - log u} / v, so that log(gamma) + log(xi + eg)
  ## = log(gamma * xi + gamma * eg) has a finite limit as gamma -> 0
  gxi <- gamma * xi
  ll <- log(alpha) + (alpha - 1) * log(ly) - log(yv) +
        expm1(-gamma * w) - gxi * w + log(gxi + gamma * eg)

  if (any(!is.finite(ll))) return(rep(-1e10, nn))
  ll
}

###############################################################################
## 4. ESTIMATION : multi-start over starting values and maxLik algorithms
###############################################################################

## starting values for beta from a quasi-binomial GLM with a logit link
glm0  <- glm(y ~ X[, -1], family = quasibinomial(link = "logit"))
beta0 <- unname(coef(glm0))

## A run is accepted only if it converged, reached a finite log-likelihood,
## has a negative definite Hessian, and sits away from the boundary of the
## parameter space. As gamma -> 0 the model tends to the unit-Weibull quantile
## regression (Proposition 2), which is fitted separately as a competitor; a run
## with gamma below GAMMA_MIN is therefore close to that limit, where the
## information about gamma vanishes. 01_estimation.R reports the profile
## log-likelihood at gamma = GAMMA_MIN, which shows that the rule does not
## affect the estimate.
GAMMA_MIN <- 1e-3
ALPHA_RNG <- c(1e-2, 40)

## Two separate things have to be handled here.
##
## 1. maxLik's Newton-type methods (NR, BHHH) invert the Hessian at every step
##    via qr.solve. When the Hessian is singular at some iterate, maxNR calls
##    it inside try(silent = FALSE), so the message
##       Error in qr.solve(H[!fixed, !fixed, ...]) : singular matrix 'a' in solve
##    is PRINTED to stderr but is NOT propagated: the run continues and
##    returns an object with a non-zero convergence code. tryCatch therefore
##    cannot suppress it, because nothing is actually being thrown. The
##    message is noise, not a crash. It is silenced below by diverting the
##    message connection for the duration of the call.
##
## 2. Genuine errors and warnings are caught by tryCatch and turned into NULL,
##    so that one bad start can never abort the multi-start loop.
safe_maxLik <- function(llf, start, method) {
  tmp <- tempfile(); con <- file(tmp, open = "wt")
  sink(con, type = "message")
  on.exit({ sink(type = "message"); close(con); unlink(tmp) }, add = TRUE)
  tryCatch(
    suppressWarnings(
      maxLik(logLik = llf, start = start, method = method,
             control = list(iterlim = 1000, printLevel = 0))),
    error = function(e) NULL, warning = function(w) NULL)
}

## Returns a list(fit = <maxLik object or NULL>, ok = TRUE/FALSE, why = "...").
## The fitted object is returned even when rejected, so that a diagnosis never
## requires running the optimiser a second time.
fit_one <- function(start, method, llf) {
  out <- safe_maxLik(llf, start, method)
  bad <- function(w) list(fit = out, ok = FALSE, why = w)

  if (is.null(out))                                  return(bad("error"))
  if (is.null(out$code) || !(out$code %in% 0:2))     return(bad("no convergence"))
  if (!is.finite(logLik(out)) || logLik(out) < -1e7) return(bad("penalty region"))
  th <- coef(out)
  if (exp(th[2]) < GAMMA_MIN)                        return(bad("gamma at floor"))
  if (exp(th[1]) < ALPHA_RNG[1] || exp(th[1]) > ALPHA_RNG[2])
                                                     return(bad("alpha out of range"))
  H <- out$hessian
  if (is.null(H) || any(!is.finite(H)))              return(bad("bad Hessian"))
  ## The numerical Hessian is symmetric in theory but not to the last bit in
  ## floating point. eigen() then treats it as a general matrix and returns
  ## COMPLEX eigenvalues, on which max() fails with
  ##   "invalid 'type' (complex) of argument".
  ## Symmetrising and passing symmetric = TRUE forces the real symmetric
  ## routine and guarantees real eigenvalues.
  H  <- (H + t(H)) / 2
  ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || !all(is.finite(ev)))            return(bad("eigen failed"))
  if (max(ev) >= 0)                                  return(bad("not a maximum"))
  if (max(abs(ev)) / min(abs(ev)) > 1e12)            return(bad("singular Hessian"))
  out$hessian <- H
  list(fit = out, ok = TRUE, why = "ok")
}

## Hessian used for inference. The Hessian returned by maxLik is a one-step
## finite difference of a numerical gradient, and it amplifies differences of
## order 1e-9 in the optimum (which arise between platforms and between
## optimisers) into changes of about 0.1% in the standard errors. The
## Richardson-extrapolated Hessian of numDeriv is stable to about 1e-9 at the
## same points, so the standard errors, p-values and GD_i do not depend on the
## platform. maxLik's Hessian is still used for the acceptance rules above.
accurate_hessian <- function(llf, par, fallback = NULL) {
  H <- tryCatch(numDeriv::hessian(function(p) sum(llf(p)), par), error = function(e) NULL)
  if (is.null(H) || any(!is.finite(H))) return(fallback)
  H <- (H + t(H)) / 2; dimnames(H) <- list(names(par), names(par)); H }

grid    <- expand.grid(la = log(c(0.5, 1, 2, 4)),
                       lg = log(c(0.02, 0.1, 0.5, 2)))
## BFGS and Nelder-Mead do not invert a Hessian at each step, so they survive
## an ill-conditioned problem that makes NR and BHHH fail outright. NR and BHHH
## are kept for stage 2 because when they do work they are more accurate.
METHODS <- c("BFGS", "NM", "NR", "BHHH")

## Two-stage multi-start.
##   Stage 1 runs the cheap gradient method over the whole starting grid.
##   Stage 2 re-runs only the most promising starts with every method.
## The earlier version ran 16 starts x 4 methods, and re-ran every rejected
## fit a second time for diagnosis: 128 optimisations of a 10-parameter
## likelihood per quantile, which takes many minutes. This does 16 + 3 x 4 = 28
## and reuses the rejected objects instead of recomputing them.
estimate_qungw <- function(u, verbose = TRUE, Xm = X, yv = y, b0 = beta0) {
  llf <- make_loglik(u, Xm, yv)
  starts <- lapply(seq_len(nrow(grid)), function(g) {
    st <- c(logAlpha = grid$la[g], logGamma = grid$lg[g], b0)
    names(st)[-(1:2)] <- colnames(Xm); st
  })

  if (verbose) cat(sprintf("u = %.2f  stage 1 (BFGS, %d starts): ",
                           u, length(starts)))
  s1 <- lapply(starts, function(st) {
    r <- fit_one(st, "BFGS", llf)
    if (verbose) cat(if (r$ok) "+" else ".")
    r
  })
  if (verbose) cat("\n")

  ll_of <- function(r) if (is.null(r$fit)) -Inf else
    suppressWarnings(as.numeric(logLik(r$fit)))
  ord  <- order(vapply(s1, ll_of, numeric(1)), decreasing = TRUE)
  top  <- starts[head(ord, 3)]

  if (verbose) cat("           stage 2 (all methods, top 3 starts): ")
  s2 <- unlist(lapply(top, function(st)
    lapply(METHODS, function(m) {
      r <- fit_one(st, m, llf)
      if (verbose) cat(if (r$ok) "+" else ".")
      r
    })), recursive = FALSE)
  if (verbose) cat("\n")

  all_r <- c(s1, s2)
  ok    <- Filter(function(r) r$ok, all_r)

  if (verbose)
    cat(sprintf("           %d of %d runs accepted\n",
                length(ok), length(all_r)))

  if (!length(ok)) {
    why <- table(vapply(all_r, function(r) r$why, character(1)))
    cat("\nNo run was accepted. Reasons:\n"); print(why)
    cat("  'singular Hessian' or 'not a maximum' means the likelihood is flat\n")
    cat("  in some direction, usually near-collinear covariates: check the\n")
    cat("  VIF and correlation output printed above.\n")
    cat("  'gamma at floor' means the optimum is degenerate.\n")
    return(NULL)
  }

  best_i <- which.max(vapply(ok, ll_of, numeric(1)))
  f <- ok[[best_i]]$fit
  f$hessian <- accurate_hessian(llf, coef(f), f$hessian)
  f
}

best <- estimate_qungw(TAU)
if (is.null(best)) {
  stop("QUNGW: no run converged to an interior maximum. ",
       "Try extending the gamma grid downwards.")
}

cat(sprintf("\nQUNGW logLik = %.4f\n", logLik(best)))

## standard errors. The delta method moves alpha and gamma off the log scale,
## since d alpha / d log(alpha) = alpha.
est_t <- coef(best)
Hb    <- best$hessian; Hb <- (Hb + t(Hb)) / 2
V_t   <- solve(-Hb)
J     <- diag(c(exp(est_t[1]), exp(est_t[2]), rep(1, p)))
est   <- c(alpha = exp(est_t[1]), gamma = exp(est_t[2]), est_t[-(1:2)])
V     <- J %*% V_t %*% t(J)
se    <- sqrt(diag(V))
zval  <- est / se
pval  <- 2 * pnorm(-abs(zval))
## sandwich covariance of White (1982), A^{-1} B A^{-1}, with A the observed
## information and B the outer product of the per-observation scores; it remains
## valid when the model is misspecified
S_obs <- numDeriv::jacobian(make_loglik(TAU), est_t)
V_s   <- J %*% (V_t %*% crossprod(S_obs) %*% V_t) %*% t(J)
se_s  <- sqrt(diag(V_s)); p_s <- 2 * pnorm(-abs(est / se_s))

## 95% Wald intervals exactly as defined in the paper (Eq. ci): on the log
## scale for the positive shape parameters, mapped back by exponentiation;
## on the natural scale for the regression coefficients. These are also the
## intervals whose coverage the simulation study evaluates.
z975  <- qnorm(0.975)
se_t  <- sqrt(diag(V_t))                 # log-scale SEs for the shapes
ci_lo <- c(exp(est_t[1] - z975 * se_t[1]), exp(est_t[2] - z975 * se_t[2]),
           est[-(1:2)] - z975 * se[-(1:2)])
ci_hi <- c(exp(est_t[1] + z975 * se_t[1]), exp(est_t[2] + z975 * se_t[2]),
           est[-(1:2)] + z975 * se[-(1:2)])

qungw_tab <- data.frame(Estimate = est, SE = se, z = zval, p = pval,
                        ci95_lo = ci_lo, ci95_hi = ci_hi)
rownames(qungw_tab) <- c("alpha", "gamma", colnames(X))
cat(sprintf("\n--- QUNGW regression at u = %.2f ---\n", TAU))
print(round(qungw_tab, 4))

ll_q  <- as.numeric(logLik(best)); k_q <- length(est_t)
AIC_q <- -2 * ll_q + 2 * k_q
BIC_q <- -2 * ll_q + k_q * log(n)

## report how close the fit sits to the lambda > 0 boundary
alpha_h <- est[1]; gamma_h <- est[2]
mu_h    <- 1 / (1 + exp(-as.vector(X %*% est[-(1:2)])))
slack   <- gamma_h * (-log(mu_h))^alpha_h
mu_min  <- mu_min_attainable(alpha_h, gamma_h, TAU)
n_below <- sum(y < mu_min)

cat("\n--- Feasibility of the fitted model ---\n")
## NOTE. At top level R parses `if (...) expr` as complete as soon as the line
## ends, so a bare `else` on the next line is a syntax error. The braces are
## required here, not decorative.
if (is.finite(feasible_bound(TAU))) {
  cat(sprintf("max gamma*(-log mu)^alpha = %.4f of bound %.4f  (%.0f%%)\n",
              max(slack), feasible_bound(TAU),
              100 * max(slack) / feasible_bound(TAU)))
} else {
  cat(sprintf("u = %.2f <= exp(-1), so lambda > 0 holds for every mu and the\n",
              TAU))
  cat("restriction of Remark 1 is vacuous. No boundary issue.\n")
}
cat(sprintf("smallest attainable conditional quantile mu_min = %.4f\n", mu_min))
cat(sprintf("smallest observed response               min(y) = %.4f\n", min(y)))
cat(sprintf("observations below mu_min: %d of %d (%.0f%%)\n",
            n_below, n, 100 * n_below / n))
if (max(slack) / feasible_bound(TAU) > 0.95) {
  cat("\n  WARNING. The optimum sits on the lambda > 0 boundary, not in the\n")
  cat("  interior of the parameter space. The usual maximum likelihood\n")
  cat("  asymptotics do not apply there, so the standard errors, z values and\n")
  cat("  p values printed above are NOT valid, and the log-likelihood is a\n")
  cat("  constrained rather than an unconstrained maximum. Comparing it with\n")
  cat("  the other families on AIC or BIC is not a like-for-like comparison.\n")
  cat("  Use a quantile u <= exp(-1) = 0.3679, where the restriction is\n")
  cat("  vacuous, as the baseline instead.\n")
}

###############################################################################
## 5. COMPETING MODELS, same link and same covariates
###############################################################################

## 5a. Unit-Weibull quantile regression (Mazucheli et al., 2020)
##     F(y) = u ^ { ( log(y) / log(mu) )^b }
make_UW <- function(u, Xm = X, yv = y) function(par) {
  nn <- length(yv)
  b <- exp(par[1])
  mu <- 1 / (1 + exp(-as.vector(Xm %*% par[-1])))
  if (any(mu <= 1e-10 | mu >= 1 - 1e-10)) return(rep(-1e10, nn))
  r  <- log(yv) / log(mu)
  ll <- log(b) - log(yv) + log(-log(u)) - log(-log(mu)) +
        (b - 1) * log(r) + log(u) * r^b
  if (any(!is.finite(ll))) return(rep(-1e10, nn))
  ll
}

## 5b. Kumaraswamy quantile regression (Mitnik and Baek, 2013)
##     F(y) = 1 - (1 - y^a)^b, with a pinned so the uth quantile equals mu
make_KW <- function(u, Xm = X, yv = y) function(par) {
  nn <- length(yv)
  b <- exp(par[1])
  mu <- 1 / (1 + exp(-as.vector(Xm %*% par[-1])))
  if (any(mu <= 1e-10 | mu >= 1 - 1e-10)) return(rep(-1e10, nn))
  a <- log(1 - (1 - u)^(1 / b)) / log(mu)
  if (any(!is.finite(a) | a <= 0)) return(rep(-1e10, nn))
  ll <- log(a) + log(b) + (a - 1) * log(yv) + (b - 1) * log1p(-yv^a)
  if (any(!is.finite(ll))) return(rep(-1e10, nn))
  ll
}

## 5c. Beta regression (Ferrari and Cribari-Neto, 2004)
##     mean parameterisation, logit link, constant precision. It has no u,
##     so the same fit appears in every panel below; it is a conditional-mean
##     model and is included only as the standard benchmark.
make_BE <- function(u = NULL, Xm = X, yv = y) function(par) {
  nn <- length(yv)
  phi <- exp(par[1])
  mu <- 1 / (1 + exp(-as.vector(Xm %*% par[-1])))
  if (any(mu <= 1e-10 | mu >= 1 - 1e-10)) return(rep(-1e10, nn))
  ll <- lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
        (mu * phi - 1) * log(yv) + ((1 - mu) * phi - 1) * log1p(-yv)
  if (any(!is.finite(ll))) return(rep(-1e10, nn))
  ll
}

## Every competitor gets the same multi-start protection as the QUNGW model:
## a grid over its shape parameter crossed with two algorithms. A single
## default start can strand a rival at an inferior local maximum, and a gap
## computed against an under-optimised rival is not a gap.
fit_competitor <- function(maker, u, Xm = X, yv = y, b0 = beta0,
                           shape_grid = c(-2, -1, 0, 1, 2, 3)) {
  best <- NULL
  for (s0 in shape_grid) for (m in c("BFGS", "NM")) {
    r <- safe_maxLik(maker(u, Xm, yv), c(shape = s0, b0), m)
    if (!is.null(r) && is.finite(logLik(r)) &&
        (is.null(best) || logLik(r) > logLik(best))) best <- r
  }
  best
}

uw <- fit_competitor(make_UW, TAU)
kw <- fit_competitor(make_KW, TAU)
be <- fit_competitor(make_BE, TAU, shape_grid = c(1, 2, 3, 4, 5))  # log(phi)

ic <- function(m) {
  if (is.null(m)) return(c(logLik = NA, k = NA, AIC = NA, BIC = NA))
  k <- length(coef(m)); l <- as.numeric(logLik(m))
  c(logLik = l, k = k, AIC = -2 * l + 2 * k, BIC = -2 * l + k * log(n))
}
comp <- rbind(Beta           = ic(be),
              Kumaraswamy    = ic(kw),
              `Unit-Weibull` = ic(uw),
              QUNGW = c(logLik = ll_q, k = k_q, AIC = AIC_q, BIC = BIC_q))
cat("\n--- Model comparison ---\n"); print(round(comp, 3))
cat("\nCovariates:", paste(xvars, collapse = ", "),
    sprintf("  (k = %d for QUNGW, %d for the two-parameter families)\n",
            k_q, k_q - 1))
bestic <- rownames(comp)[which.min(comp[, "AIC"])]
cat("Lowest AIC:", bestic, " Lowest BIC:",
    rownames(comp)[which.min(comp[, "BIC"])], "\n")

