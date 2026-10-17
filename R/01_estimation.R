###############################################################################
## 01_estimation.R
## Tables 2-4, the quantile process (Fig. 1), the choice of u, the robustness of
## the model comparison (Vuong, bootstrap, contributions, leave-one-out with
## refits) and the missing-data check of Section 4.1.
## Requires qungw_core.R to have been sourced (run_all.R does this).
###############################################################################

if (!exists("RES")) RES <- list()
cat("\n================ 01_estimation.R ================\n")

## ---------------------------------------------------------------- data, Table 2
RES$n_response  <- sum(!is.na(dat_all$y))
RES$n           <- n
RES$y_range     <- range(y)
RES$desc        <- round(desc, 3)
RES$kappa_X     <- kappa_X
RES$vif         <- unname(vif[1])
RES$cor_x       <- cor(Xraw)[1, 2]

## wider covariate set examined first and rejected (Section 4.1)
all7 <- c("fert", "lgdppc", "urb", "dep", "serv_f", "internet", "educ_f")
d7   <- dat_all[complete.cases(dat_all[, c("y", all7)]), all7]
vif7 <- sapply(all7, function(v) 1 / (1 - summary(lm(d7[[v]] ~ as.matrix(d7[, setdiff(all7, v)])))$r.squared))
cr7  <- cor(d7); diag(cr7) <- NA
RES$vif7_range <- range(vif7); RES$cor7_max <- max(abs(cr7), na.rm = TRUE); RES$n7 <- nrow(d7)

## ---------------------------------------------------------------- Table 3
RES$tab3 <- data.frame(est = est, se = se, z = zval, p = pval, se_s = se_s, p_s = p_s,
                       lo = c(exp(est_t[1:2] - 1.96 * sqrt(diag(V_t))[1:2]), est[-(1:2)] - 1.96 * se[-(1:2)]),
                       hi = c(exp(est_t[1:2] + 1.96 * sqrt(diag(V_t))[1:2]), est[-(1:2)] + 1.96 * se[-(1:2)]))
RES$comp <- comp
cat("\nTable 3\n"); print(round(RES$tab3, 4)); print(round(comp, 3))

## ---------------------------------------------------------------- Table 4
## For u > exp(-1) the restriction of Remark 1 can make the likelihood increase
## towards the boundary of the admissible set. If no run is accepted, or if the
## accepted run is beaten by a run on the boundary, the boundary value is
## reported (it is not an interior maximum).
boundary_sup <- function(u) {
  llf <- make_loglik(u); best_l <- -Inf; info <- NULL
  g2 <- expand.grid(la = log(c(0.5, 1, 2, 3, 4, 6)), lg = log(c(0.01, 0.05, 0.1, 0.5, 1, 2, 4)))
  for (k in seq_len(nrow(g2))) for (m in c("BFGS", "NM")) {
    f <- safe_maxLik(llf, c(g2$la[k], g2$lg[k], beta0), m); if (is.null(f)) next
    l <- sum(llf(coef(f))); if (!is.finite(l) || l < -1e6) next
    if (l > best_l) { best_l <- l; e <- coef(f); a <- exp(e[1]); g <- exp(e[2])
      mu <- 1 / (1 + exp(-as.vector(X %*% e[-(1:2)])))
      info <- c(alpha = a, gamma = g, used = 100 * max(g * (-log(mu))^a) / feasible_bound(u),
                mu_min = mu_min_attainable(a, g, u)) } }
  c(logLik = best_l, info)
}
tab4 <- list()
for (uu in c(0.10, 0.25, 0.50, 0.75)) {
  f  <- estimate_qungw(uu, verbose = FALSE)
  bs <- if (uu > exp(-1)) boundary_sup(uu) else NULL
  if (!is.null(f)) {
    e <- coef(f); a <- exp(e[1]); g <- exp(e[2]); mu <- 1 / (1 + exp(-as.vector(X %*% e[-(1:2)])))
    used <- if (is.finite(feasible_bound(uu))) 100 * max(g * (-log(mu))^a) / feasible_bound(uu) else 0
    q <- c(logLik = as.numeric(logLik(f)), used = used, interior = 1)
  } else q <- c(logLik = NA, used = NA, interior = 0)
  if (!is.null(bs) && (is.na(q["logLik"]) || bs["logLik"] > q["logLik"]))
    q <- c(logLik = unname(bs["logLik"]), used = unname(bs["used"]), interior = 0)
  fk <- fit_competitor(make_KW, uu); fu <- fit_competitor(make_UW, uu)
  tab4[[as.character(uu)]] <- data.frame(u = uu, model = c("QUNGW", "Kumaraswamy", "Unit-Weibull"),
    logLik = c(q["logLik"], logLik(fk), logLik(fu)), k = c(5, 4, 4),
    AIC = c(ifelse(q["interior"] == 1, -2 * q["logLik"] + 10, NA), -2 * as.numeric(logLik(fk)) + 8, -2 * as.numeric(logLik(fu)) + 8),
    used = c(q["used"], NA, NA), interior = c(q["interior"], 1, 1))
  if (uu == 0.50) RES$median_boundary <- bs
}
RES$tab4 <- do.call(rbind, tab4); rownames(RES$tab4) <- NULL
cat("\nTable 4\n"); print(RES$tab4, digits = 5)
## Afghanistan: largest gamma allowing its fitted median to equal its observed rate
mb <- RES$median_boundary
RES$gamma_max_AFG <- feasible_bound(0.5) / (-log(min(y)))^mb["alpha"]
RES$n_below_mumin <- sum(y < mb["mu_min"])

## accepted runs at u = 0.25 (choice section)
RES$runs_accepted_025 <- local({
  llf <- make_loglik(TAU)
  starts <- lapply(seq_len(nrow(grid)), function(g) { st <- c(logAlpha = grid$la[g], logGamma = grid$lg[g], beta0); st })
  s1 <- lapply(starts, function(st) fit_one(st, "BFGS", llf))
  ll1 <- vapply(s1, function(r) if (is.null(r$fit)) -Inf else as.numeric(logLik(r$fit)), 0)
  top <- starts[head(order(ll1, decreasing = TRUE), 3)]
  s2 <- unlist(lapply(top, function(st) lapply(METHODS, function(m) fit_one(st, m, llf))), recursive = FALSE)
  sum(vapply(c(s1, s2), function(r) r$ok, TRUE))
})

## ---------------------------------------------------------------- quantile process (Fig. 1)
taus <- seq(0.05, 0.35, by = 0.025)
RES$proc <- do.call(rbind, lapply(taus, function(tt) {
  f <- estimate_qungw(tt, verbose = FALSE); e <- coef(f)
  s <- sqrt(diag(solve(-((f$hessian + t(f$hessian)) / 2))))
  data.frame(u = tt, term = colnames(X)[-1], est = e[-(1:3)], se = s[-(1:3)])
}))
RES$proc$lo <- RES$proc$est - 1.96 * RES$proc$se
RES$proc$hi <- RES$proc$est + 1.96 * RES$proc$se

## ---------------------------------------------------------------- unit-Weibull limit (Proposition 2)
## numerical check of the limit: at gamma = 1e-4 the QUNGW log-likelihood,
## evaluated with alpha = shape and beta = coefficients of the fitted
## unit-Weibull model, equals the unit-Weibull log-likelihood
ll_tot <- function(p) sum(make_loglik(TAU)(p))
RES$uw_limit <- c(qungw = ll_tot(c(coef(uw)[1], log(1e-4), coef(uw)[-1])), uw = as.numeric(logLik(uw)))
## profile log-likelihood below the acceptance floor gamma = GAMMA_MIN, on a grid
## from GAMMA_MIN down to 1e-6; its largest value shows whether the rule can
## have excluded a better optimum
g_grid <- 10^seq(log10(GAMMA_MIN), -6, by = -0.5)
prof_g <- sapply(g_grid, function(g) -optim(c(coef(uw)[1], coef(uw)[-1]), function(q) -ll_tot(c(q[1], log(g), q[-1])),
                                            method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))$value)
RES$ll_gamma_floor <- max(prof_g); RES$gamma_grid <- data.frame(gamma = g_grid, profile = prof_g)
cat(sprintf("\nUnit-Weibull limit: QUNGW(gamma=1e-4) %.6f, unit-Weibull %.6f; largest profile for 1e-6 <= gamma <= %g: %.4f\n",
            RES$uw_limit[1], RES$uw_limit[2], GAMMA_MIN, RES$ll_gamma_floor))

## refit from a given start: BFGS, then Nelder-Mead and Newton-Raphson if BFGS fails
refit <- function(llf, start) {
  for (m in c("BFGS", "NM", "NR")) { f <- safe_maxLik(llf, start, m)
    if (!is.null(f) && f$code %in% 0:2 && is.finite(logLik(f))) return(f) }
  NULL }

## ---------------------------------------------------------------- robustness of the comparison
k_q <- length(coef(best)); k_k <- length(coef(kw))
d_i  <- make_loglik(TAU)(coef(best)) - make_KW(TAU)(coef(kw))
vuong_stat <- function(d, n_) sqrt(n_) * mean(d - (k_q - k_k) / n_) / sd(d - (k_q - k_k) / n_)
RES$gap <- sum(d_i); RES$vuong <- vuong_stat(d_i, n); RES$vuong_p <- 2 * pnorm(-abs(RES$vuong))
o5 <- order(d_i, decreasing = TRUE)[1:5]
RES$top5 <- data.frame(iso = dat$iso3c[o5], contrib = d_i[o5]); RES$top5_sum <- sum(d_i[o5])
RES$decile1 <- sum(d_i[y <= quantile(y, 0.10)])
## leave-one-out with both models refitted
loo <- t(sapply(seq_len(n), function(i) {
  Xi <- X[-i, , drop = FALSE]; yi <- y[-i]
  fq <- refit(make_loglik(TAU, Xi, yi), coef(best)); fk <- refit(make_KW(TAU, Xi, yi), coef(kw))
  if (is.null(fq) || is.null(fk)) return(c(NA, NA))
  dd <- make_loglik(TAU, Xi, yi)(coef(fq)) - make_KW(TAU, Xi, yi)(coef(fk))
  c(vuong_stat(dd, n - 1), sum(dd))
}))
stopifnot(all(is.finite(loo)))                            # every leave-one-out refit converged
RES$loo_vuong <- range(loo[, 1]); RES$loo_vuong_where <- dat$iso3c[c(which.min(loo[, 1]), which.max(loo[, 1]))]
RES$loo_gap <- range(loo[, 2])
## parametric bootstrap of the gap (B_BOOT replications under each fitted model)
seed_mt(20260815)
alpha_h <- exp(coef(best)[1]); gamma_h <- exp(coef(best)[2]); mu_h <- 1 / (1 + exp(-as.vector(X %*% coef(best)[-(1:2)])))
bk <- exp(coef(kw)[1]); muk <- 1 / (1 + exp(-as.vector(X %*% coef(kw)[-1]))); ak <- log(1 - (1 - TAU)^(1 / bk)) / log(muk)
refit_gap <- function(ys) {
  fq <- refit(make_loglik(TAU, X, ys), coef(best)); fk <- refit(make_KW(TAU, X, ys), coef(kw))
  if (is.null(fq) || is.null(fk)) return(NA_real_); as.numeric(logLik(fq)) - as.numeric(logLik(fk)) }
B_BOOT <- 1000
g_q <- replicate(B_BOOT, refit_gap(rQUNGW(alpha_h, gamma_h, mu_h, TAU)))
g_k <- replicate(B_BOOT, refit_gap((1 - (1 - runif(n))^(1 / bk))^(1 / ak)))
RES$boot_fail <- c(qungw = sum(!is.finite(g_q)), kum = sum(!is.finite(g_k)))
g_q <- g_q[is.finite(g_q)]; g_k <- g_k[is.finite(g_k)]
## Monte Carlo p-value (r + 1) / (B + 1), with r the number of replicates at least as large as the observed gap
RES$boot_kum <- c(median = median(g_k), q975 = unname(quantile(g_k, 0.975)), p = mean(g_k >= RES$gap),
                  r = sum(g_k >= RES$gap), B = length(g_k), p_mc = (sum(g_k >= RES$gap) + 1) / (length(g_k) + 1), max = max(g_k))
RES$boot_q   <- c(median = median(g_q), pctile = 100 * mean(g_q <= RES$gap), B = length(g_q))
cat(sprintf("\nVuong %.3f (p %.3f); LOO refit range %.3f-%.3f; bootstrap Kum median %.2f q975 %.2f p %.4f; QUNGW median %.2f pctile %.0f\n",
            RES$vuong, RES$vuong_p, RES$loo_vuong[1], RES$loo_vuong[2], RES$boot_kum[1], RES$boot_kum[2],
            RES$boot_kum[3], RES$boot_q[1], RES$boot_q[2]))

## ---------------------------------------------------------------- missing-data check
drop_idx <- !complete.cases(dat_all[, c("y", "fert", "lgdppc")]) & !is.na(dat_all$y)
D <- dat_all[drop_idx, ]; K <- dat_all[!drop_idx & !is.na(dat_all$y), ]
y_all <- dat_all$y[!is.na(dat_all$y)]
RES$miss <- list(n_excluded = nrow(D), min_excluded = min(D$y), iso_min = D$iso3c[which.min(D$y)],
                 in_lowest_decile = sum(D$y < quantile(y_all, 0.10)), expected = 0.1 * nrow(D),
                 wilcox_y = wilcox.test(D$y, K$y)$p.value, wilcox_fert = wilcox.test(D$fert, K$fert)$p.value)
## fertility is standardised with the mean and SD of the full sample in both
## fits, so that the two coefficients are on the same scale
fert_all <- dat_all$fert[!is.na(dat_all$y) & !is.na(dat_all$fert)]
fit1 <- function(d) { yv <- d$y; Xm <- cbind(1, (d$fert - mean(fert_all)) / sd(fert_all)); colnames(Xm) <- c("(Intercept)", "fert")
  b0 <- unname(coef(glm(yv ~ Xm[, -1], family = quasibinomial(link = "logit"))))
  q <- estimate_qungw(TAU, verbose = FALSE, Xm = Xm, yv = yv, b0 = b0); k <- fit_competitor(make_KW, TAU, Xm = Xm, yv = yv, b0 = b0)
  c(beta_fert = unname(coef(q)[4]), dAIC = (-2 * as.numeric(logLik(q)) + 8) - (-2 * as.numeric(logLik(k)) + 6)) }
RES$miss$full <- fit1(dat_all[!is.na(dat_all$y) & !is.na(dat_all$fert), ])
RES$miss$est  <- fit1(dat_all[complete.cases(dat_all[, c("y", "fert", "lgdppc")]), ])
str(RES$miss)
