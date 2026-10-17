###############################################################################
## 03_residuals.R
## Residual analysis of Section 4.6: Table 6 (KS on Cox-Snell, SW on quantile
## residuals, with parametric-bootstrap p-values) and the simulated envelopes of Fig. 4 (200 refitted samples,
## Atkinson, 1981). Requires qungw_core.R.
###############################################################################

cat("\n================ 03_residuals.R ================\n")
seed_mt(20261018)
cl <- function(F) pmin(pmax(F, 1e-12), 1 - 1e-12)
MODELS <- list(
  QUNGW = list(par = coef(best), mk = function(Xm, yv) make_loglik(TAU, Xm, yv),
    cdf  = function(p, yv) { m <- 1 / (1 + exp(-as.vector(X %*% p[-(1:2)]))); pQUNGW(yv, exp(p[1]), exp(p[2]), m, TAU) },
    rgen = function(p) { m <- 1 / (1 + exp(-as.vector(X %*% p[-(1:2)]))); qQUNGW(runif(n), exp(p[1]), exp(p[2]), m, TAU) }),
  Kumaraswamy = list(par = coef(kw), mk = function(Xm, yv) make_KW(TAU, Xm, yv),
    cdf  = function(p, yv) { b <- exp(p[1]); m <- 1 / (1 + exp(-as.vector(X %*% p[-1]))); a <- log(1 - (1 - TAU)^(1 / b)) / log(m); 1 - (1 - yv^a)^b },
    rgen = function(p) { b <- exp(p[1]); m <- 1 / (1 + exp(-as.vector(X %*% p[-1]))); a <- log(1 - (1 - TAU)^(1 / b)) / log(m); (1 - (1 - runif(n))^(1 / b))^(1 / a) }),
  Beta = list(par = coef(be), mk = function(Xm, yv) make_BE(TAU, Xm, yv),
    cdf  = function(p, yv) { ph <- exp(p[1]); m <- 1 / (1 + exp(-as.vector(X %*% p[-1]))); pbeta(yv, m * ph, (1 - m) * ph) },
    rgen = function(p) { ph <- exp(p[1]); m <- 1 / (1 + exp(-as.vector(X %*% p[-1]))); rbeta(n, m * ph, (1 - m) * ph) }),
  `Unit-Weibull` = list(par = coef(uw), mk = function(Xm, yv) make_UW(TAU, Xm, yv),
    cdf  = function(p, yv) { b <- exp(p[1]); m <- 1 / (1 + exp(-as.vector(X %*% p[-1]))); TAU^((log(yv) / log(m))^b) },
    rgen = function(p) { b <- exp(p[1]); m <- 1 / (1 + exp(-as.vector(X %*% p[-1]))); exp(log(m) * (log(runif(n)) / log(TAU))^(1 / b)) }))

## Table 6
RES$tab6 <- do.call(rbind, lapply(names(MODELS), function(k) { M <- MODELS[[k]]
  F <- cl(M$cdf(M$par, y)); e <- -log(1 - F); r <- qnorm(F)
  ks <- suppressWarnings(ks.test(e, "pexp", 1)); sw <- shapiro.test(r)
  data.frame(model = k, KS = unname(ks$statistic), KS_p = ks$p.value, CS_mean = mean(e),
             SW = unname(sw$statistic), SW_p = sw$p.value,
             skew = mean((r - mean(r))^3) / sd(r)^3, min_r = min(r)) }))
print(RES$tab6, digits = 3)

## simulated envelopes (Fig. 4)
## The same parametric bootstrap gives the envelopes and calibrates the KS and
## SW statistics of Table 6: the residuals of each sample generated from the
## fitted model are computed after refitting the model to that sample, so that
## the bootstrap null distribution includes the effect of estimation. The
## p-values are (r + 1) / (B + 1), with r the number of replicates at least as
## extreme as the observed statistic (large KS, small SW).
B_ENV <- 200
stat_of <- function(F) { F <- cl(F); r <- qnorm(F)
  c(KS = unname(suppressWarnings(ks.test(-log(1 - F), "pexp", 1))$statistic), SW = unname(shapiro.test(r)$statistic)) }
ENV <- lapply(MODELS, function(M) {
  r <- sort(qnorm(cl(M$cdf(M$par, y))))
  B <- lapply(seq_len(B_ENV), function(b) { ys <- pmin(pmax(M$rgen(M$par), 1e-10), 1 - 1e-10)
    f <- NULL; for (m in c("BFGS", "NM", "NR")) { f <- safe_maxLik(M$mk(X, ys), M$par, m)
      if (!is.null(f) && f$code %in% 0:2) break; f <- NULL }
    if (is.null(f)) return(NULL)
    Fs <- M$cdf(coef(f), ys); list(r = sort(qnorm(cl(Fs))), st = stat_of(Fs)) })
  ok <- !vapply(B, is.null, logical(1)); B <- B[ok]
  E  <- sapply(B, `[[`, "r"); ST <- t(sapply(B, `[[`, "st"))
  list(r = r, lo = apply(E, 1, quantile, 0.05), hi = apply(E, 1, quantile, 0.95), m = rowMeans(E),
       stats = ST, fail = sum(!ok)) })
save_rds(ENV, file.path(OUTDIR, "envelopes.rds"))
RES$tab6$KS_pb <- mapply(function(k, s) (sum(ENV[[k]]$stats[, "KS"] >= s) + 1) / (nrow(ENV[[k]]$stats) + 1), RES$tab6$model, RES$tab6$KS)
RES$tab6$SW_pb <- mapply(function(k, s) (sum(ENV[[k]]$stats[, "SW"] <= s) + 1) / (nrow(ENV[[k]]$stats) + 1), RES$tab6$model, RES$tab6$SW)
RES$env_B <- B_ENV; RES$env_fail <- sapply(ENV, `[[`, "fail")
print(RES$tab6[, c("model", "KS", "KS_p", "KS_pb", "SW", "SW_p", "SW_pb")], digits = 3); print(RES$env_fail)
thq <- qnorm(ppoints(n)); lt <- thq < -1.3
RES$env <- sapply(ENV, function(e) c(lower_tail_n = sum(lt), lower_tail_below = sum(e$r[lt] < e$lo[lt]),
                                     above_center = sum((e$r > e$hi)[thq > -1.6 & thq < -0.2])))
## economies whose QUNGW residual lies below the envelope in the lower tail
local({ o <- order(qnorm(cl(MODELS$QUNGW$cdf(MODELS$QUNGW$par, y)))); e <- ENV$QUNGW
  RES$env_q_below_iso <<- dat$iso3c[o][lt & e$r < e$lo] })
print(RES$env); print(RES$env_q_below_iso)
