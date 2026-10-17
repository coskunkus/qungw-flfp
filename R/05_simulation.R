###############################################################################
##  05_simulation.R : parameter-recovery simulation for the QUNGW model
##  ---------------------------------------------------------------------------
##  Does not need any other script or the data. Run by run_all.R with RUN_SIM=1,
##  or on its own from the R folder:  Rscript 05_simulation.R
##
##  DESIGN, deliberately minimal (the scenario a referee will ask about):
##    * The data-generating process is the QUNGW model CALIBRATED TO THE
##      EMPIRICAL FIT of the paper: u = 0.25, alpha = 3.8916, gamma = 3.7652,
##      beta = (-0.6404, -0.2333, 0.3925).
##    * Two standardised covariates with correlation -0.83, matching the
##      empirical correlation of fertility and log GDP per capita, so the
##      simulation inherits the application's collinearity rather than
##      assuming it away.
##    * Sample sizes n = 100, 172, 300: below, at, and above the application.
##    * REPS = 2000 replications per n. Replication r of sample size n uses
##      its own L'Ecuyer-CMRG stream, derived from SEED + n, so that the
##      results do not depend on the number of cores or on the platform
##      (forked workers on Linux/macOS, a socket cluster on Windows).
##  Reported per parameter and n: bias, SD, RMSE, and empirical coverage of
##  the nominal 95% Wald interval (on the log scale for alpha and gamma,
##  transformed back, which is how the paper reports them).
##
##  Output: results/report_sim_n*.txt and results/simulation_n*.rds
###############################################################################

## ---- settings --------------------------------------------------------------
if (!exists("ROOT")) ROOT <- normalizePath(file.path(getwd(), ".."))
if (!exists("save_rds")) source("io_utils.R")   # .rds and figure I/O safe for non-ASCII paths
WORKDIR <- file.path(ROOT, "results"); dir.create(WORKDIR, showWarnings = FALSE)
if (!dir.exists(WORKDIR)) dir.create(WORKDIR, recursive = TRUE)
## Do NOT setwd() here. These scripts used to switch the working directory to
## WORKDIR, which meant that a later source("other_script.R") with a relative
## path resolved against C:/ck instead of the code folder and silently ran a
## stale copy. Every path below is absolute, built from WORKDIR, so the working
## directory is left exactly as the caller set it.

REPS  <- 2000
if (!exists("NS")) NS <- c(100, 172, 300)   # Table 1; run with NS <- c(500, 1000) for the check in the text
U     <- 0.25
TRUEP <- c(alpha = 3.8916, gamma = 3.7652,
           b0 = -0.6404, b1 = -0.2333, b2 = 0.3925)
RHO   <- -0.83                      # cor(fert, lgdppc) in the application
SEED  <- 20260815

pkgs <- c("maxLik", "lamW", "numDeriv")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
library(maxLik)
HAS_LAMW <- requireNamespace("lamW", quietly = TRUE)

set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

## ---- QUNGW ingredients (self-contained copies) ------------------------------
xi_fun <- function(mu, alpha, gamma, u) {
  w <- (-log(mu))^alpha
  (1 / gamma) * (1 / w) * (exp(-gamma * w) - log(u) - 1)
}

pQ <- function(q, alpha, gamma, mu, u) {
  xi <- xi_fun(mu, alpha, gamma, u)
  w  <- (-log(q))^alpha
  exp(exp(-gamma * w) - xi * gamma * w - 1)
}

## closed-form quantile of the UNGW(alpha, gamma, lambda) via Lambert W (theta = lambda)
## (Eq. 5-6 of the paper); used to simulate fast. Falls back to numerical
## inversion when lamW is unavailable.
q_theta <- function(p, alpha, gamma, theta) {
  if (HAS_LAMW) {
    xi_p <- theta * lamW::lambertW0(exp((1 + log(p)) / theta) / theta) -
            1 - log(p)
    exp(-(xi_p / (theta * gamma))^(1 / alpha))
  } else {
    vapply(seq_along(p), function(i)
      uniroot(function(t)
        exp(exp(-gamma * (-log(t))^alpha) -
            theta * gamma * (-log(t))^alpha - 1) - p[i],
        c(1e-12, 1 - 1e-12), tol = 1e-11)$root, numeric(1))
  }
}

rQ <- function(alpha, gamma, mu, u) {
  th <- xi_fun(mu, alpha, gamma, u)
  q_theta(runif(length(mu)), alpha, gamma, th)
}

make_ll <- function(u, Xm, yv) function(par) {
  nn <- length(yv)
  alpha <- exp(par[1]); gamma <- exp(par[2]); beta <- par[-(1:2)]
  if (!is.finite(alpha) || alpha < 1e-3 || alpha > 50)  return(rep(-1e10, nn))
  if (!is.finite(gamma) || gamma < 1e-3 || gamma > 1e3) return(rep(-1e10, nn))
  mu <- 1 / (1 + exp(-as.vector(Xm %*% beta)))
  if (any(mu <= 1e-10 | mu >= 1 - 1e-10))               return(rep(-1e10, nn))
  xi <- xi_fun(mu, alpha, gamma, u)
  ly <- -log(yv); w <- ly^alpha; eg <- exp(-gamma * w)
  if (any(!is.finite(w)) || any(!is.finite(eg)))        return(rep(-1e10, nn))
  if (any(!is.finite(xi)) || any(xi <= 0)) {
    viol <- sum(pmax(0, -xi[is.finite(xi)])) + sum(!is.finite(xi))
    return(rep(-1e8 * (1 + viol) / nn, nn))
  }
  ll <- log(alpha) + log(gamma) + (alpha - 1) * log(ly) - log(yv) -
        (1 - eg) - xi * gamma * w + log(xi + eg)
  if (any(!is.finite(ll))) return(rep(-1e10, nn))
  ll
}

safe_fit <- function(llf, start) {
  out <- tryCatch(
    suppressWarnings(maxLik(llf, start = start, method = "BFGS",
                            control = list(iterlim = 500, printLevel = 0))),
    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(out) || is.null(out$code) || !(out$code %in% 0:2)) return(NULL)
  if (!is.finite(logLik(out)) || logLik(out) < -1e7)             return(NULL)
  th <- coef(out)
  if (exp(th[2]) < 2e-3)                                         return(NULL)
  H <- out$hessian
  if (is.null(H) || any(!is.finite(H)))                          return(NULL)
  H <- (H + t(H)) / 2
  ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || max(ev) >= 0 ||
      max(abs(ev)) / min(abs(ev)) > 1e12)                        return(NULL)
  out$hessian <- H
  out
}

## ---- one replication --------------------------------------------------------
one_rep <- function(n) {
  ## correlated standardised covariates, as in the application
  z1 <- rnorm(n); z2 <- RHO * z1 + sqrt(1 - RHO^2) * rnorm(n)
  X  <- cbind(1, z1, z2)
  mu <- 1 / (1 + exp(-as.vector(X %*% TRUEP[3:5])))
  y  <- rQ(TRUEP["alpha"], TRUEP["gamma"], mu, U)
  llf <- make_ll(U, X, y)

  ## three warm starts: the truth, and two perturbations. In a simulation the
  ## truth is known, so this is legitimate and standard; the empirical
  ## application uses the full cold-start grid instead.
  starts <- list(
    c(log(TRUEP["alpha"]), log(TRUEP["gamma"]), TRUEP[3:5]),
    c(log(2), log(1), 0, 0, 0),
    c(log(6), log(6), TRUEP[3:5]))

  best <- NULL
  for (st in starts) {
    f <- safe_fit(llf, st)
    if (!is.null(f) && (is.null(best) || logLik(f) > logLik(best))) best <- f
  }
  if (is.null(best)) return(NULL)

  e  <- coef(best)
  ## Richardson Hessian, as in the application (see accurate_hessian in qungw_core.R)
  H  <- tryCatch(numDeriv::hessian(function(p) sum(llf(p)), e), error = function(err) NULL)
  if (is.null(H) || any(!is.finite(H))) H <- best$hessian
  se <- sqrt(diag(solve(-(H + t(H)) / 2)))
  est <- c(exp(e[1]), exp(e[2]), e[3:5])
  ## 95% Wald intervals: log scale for the shapes, natural for the betas
  lo <- c(exp(e[1] - 1.96 * se[1]), exp(e[2] - 1.96 * se[2]),
          e[3:5] - 1.96 * se[3:5])
  hi <- c(exp(e[1] + 1.96 * se[1]), exp(e[2] + 1.96 * se[2]),
          e[3:5] + 1.96 * se[3:5])
  cov <- as.numeric(lo <= TRUEP & TRUEP <= hi)
  c(est, cov)
}

## ---- parallel map that works on every platform -----------------------------
par_map <- function(idx, FUN) {
  nc <- max(1L, parallel::detectCores(), na.rm = TRUE)
  if (.Platform$OS.type != "windows") return(parallel::mclapply(idx, FUN, mc.cores = nc))
  cl <- parallel::makeCluster(nc); on.exit(parallel::stopCluster(cl))
  env <- environment(FUN)
  while (!identical(env, globalenv()) && !identical(env, emptyenv())) {
    parallel::clusterExport(cl, ls(env, all.names = TRUE), envir = env); env <- parent.env(env) }
  parallel::clusterExport(cl, ls(globalenv()), envir = globalenv())
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(maxLik)))
  parallel::parLapply(cl, idx, FUN) }

## ---- run --------------------------------------------------------------------
REPORT <- file.path(WORKDIR, paste0("report_sim_n", paste(NS, collapse = "_"), ".txt"))
sink(REPORT, split = TRUE)
cat("=============================================================\n")
cat("PARAMETER-RECOVERY SIMULATION, QUNGW at u =", U, "\n")
cat("DGP calibrated to the empirical fit; cor(x1,x2) =", RHO, "\n")
cat("Replications per n:", REPS, "  Seed:", SEED, "\n")
if (!HAS_LAMW) cat("(lamW not installed: simulating by numerical inversion, slower)\n")
cat("=============================================================\n")

t0 <- Sys.time()
res <- list()
for (n in NS) {
  cat(sprintf("\nn = %d: ", n))
  out <- matrix(NA_real_, REPS, 10)
  fails <- 0
  ## one random-number stream per replication
  set.seed(SEED + n, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  seeds <- vector("list", REPS); seeds[[1]] <- .Random.seed
  for (r in seq_len(REPS)[-1]) seeds[[r]] <- parallel::nextRNGStream(seeds[[r - 1]])
  run_rep <- function(r) { assign(".Random.seed", seeds[[r]], envir = globalenv()); one_rep(n) }
  lst <- par_map(seq_len(REPS), run_rep)
  for (r in seq_len(REPS)) { o <- lst[[r]]; if (is.null(o) || length(o) != 10) fails <- fails + 1 else out[r, ] <- o }
  cat("\n")
  ok <- out[stats::complete.cases(out), , drop = FALSE]
  est <- ok[, 1:5, drop = FALSE]; covm <- ok[, 6:10, drop = FALSE]
  tab <- data.frame(
    n         = n,
    parameter = c("alpha", "gamma", "beta0", "beta1 (fert)", "beta2 (lgdppc)"),
    true      = as.numeric(TRUEP),
    bias      = colMeans(est) - as.numeric(TRUEP),
    sd        = apply(est, 2, sd),
    rmse      = sqrt(colMeans((t(t(est) - as.numeric(TRUEP)))^2)),
    coverage  = colMeans(covm))
  tab$convergence <- nrow(ok) / REPS
  res[[as.character(n)]] <- tab

  pr <- tab
  pr[, 3:8] <- round(pr[, 3:8], 3)
  print(pr, row.names = FALSE)
}
cat(sprintf("\nTotal runtime: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

all_tab <- do.call(rbind, res); save_rds(all_tab, file.path(WORKDIR, paste0("simulation_n", paste(NS, collapse = "_"), ".rds")))
## No spreadsheet output in this project: the table above is the deliverable and
## it is already in the report file.

cat("\nReading guide. Bias should shrink toward zero and RMSE should fall\n")
cat("roughly like 1/sqrt(n) as n grows; coverage should sit near 0.95.\n")
cat("gamma is the hardest parameter (its information is smallest), so its\n")
cat("RMSE will be the largest; that is expected, not a defect.\n")
cat("Written:", basename(REPORT), "\n")
cat("=============================================================\n")
sink()
cat("Report written to", REPORT, "\n")
