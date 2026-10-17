###############################################################################
## 02_influence.R
## Influence diagnostics of Section 2.5 / 4.5 (Fig. 3, Table 5):
## exact case deletion (GD_i, LD_i), local influence of Cook (1986) under
## case-weight and covariate perturbation (C_i, d_max), and relative changes
## upon deleting the flagged economies. Requires qungw_core.R and 01_estimation.R.
###############################################################################

cat("\n================ 02_influence.R ================\n")
th <- coef(best); ll_full <- make_loglik(TAU)
H  <- (best$hessian + t(best$hessian)) / 2; Hinv <- solve(H)

## global influence: exact case deletion, alpha and gamma on the log scale
del <- t(sapply(seq_len(n), function(i) {
  f <- safe_maxLik(make_loglik(TAU, X[-i, , drop = FALSE], y[-i]), th, "BFGS")
  if (is.null(f) || !(f$code %in% 0:2)) f <- safe_maxLik(make_loglik(TAU, X[-i, , drop = FALSE], y[-i]), th, "NM")
  coef(f) }))
D  <- sweep(del, 2, th)
GD <- rowSums((D %*% (-H)) * D)
LD <- 2 * (sum(ll_full(th)) - apply(del, 1, function(p) sum(ll_full(p))))

## local influence: C_d = 2|d' Delta' L^-1 Delta d|
curv <- function(Delta) { B <- -t(Delta) %*% Hinv %*% Delta; ev <- eigen((B + t(B)) / 2, symmetric = TRUE)
  list(Cmax = 2 * ev$values[1], dmax = abs(ev$vectors[, 1]), Ci = 2 * abs(diag(B))) }
S  <- numDeriv::jacobian(ll_full, th)                    # score contributions (case weights)
LI <- list(case = curv(t(S)))
for (j in 2:3) {                                         # x_ij + omega_i, s_j = 1 (standardised)
  dfdw <- function(p) { e <- 1e-5; Xp <- X; Xm <- X; Xp[, j] <- X[, j] + e; Xm[, j] <- X[, j] - e
    (make_loglik(TAU, Xp, y)(p) - make_loglik(TAU, Xm, y)(p)) / (2 * e) }
  LI[[colnames(X)[j]]] <- curv(t(numDeriv::jacobian(dfdw, th)))
}
INFL <- list(lab = dat$iso3c, GD = GD, LD = LD, LI = LI)
save_rds(INFL, file.path(OUTDIR, "influence.rds"))
top <- function(v, k = 4) dat$iso3c[order(v, decreasing = TRUE)[1:k]]
RES$infl <- list(GD_top = setNames(sort(GD, decreasing = TRUE)[1:4], top(GD)),
                 LD_top = setNames(sort(LD, decreasing = TRUE)[1:4], top(LD)),
                 dmax_case_top = top(LI$case$dmax, 2), Ci_fert_top = top(LI$fert$Ci, 1),
                 Ci_lgdp_top = top(LI$lgdppc$Ci, 1), Cmax = sapply(LI, `[[`, "Cmax"))
print(RES$infl)

## impact of the flagged economies (Table 5)
p_of <- function(f) { Hh <- (f$hessian + t(f$hessian)) / 2; Vt <- solve(-Hh); e <- coef(f)
  J <- diag(c(exp(e[1:2]), rep(1, length(e) - 2))); est_ <- c(exp(e[1:2]), e[-(1:2)])
  se_ <- sqrt(diag(J %*% Vt %*% t(J))); list(est = est_, p = 2 * pnorm(-abs(est_ / se_))) }
full <- p_of(best)
sets <- list(None = character(0), Iran = "IRN", Iraq = "IRQ", `Puerto Rico` = "PRI",
             `Saudi Arabia` = "SAU", `All four` = c("IRN", "IRQ", "PRI", "SAU"))
RES$tab5 <- do.call(rbind, lapply(names(sets), function(k) {
  keep <- !dat$iso3c %in% sets[[k]]; Xk <- X[keep, , drop = FALSE]; yk <- y[keep]
  b0 <- unname(coef(glm(yk ~ Xk[, -1], family = quasibinomial(link = "logit"))))
  f  <- estimate_qungw(TAU, verbose = FALSE, Xm = Xk, yv = yk, b0 = b0)
  fk <- fit_competitor(make_KW, TAU, Xm = Xk, yv = yk, b0 = b0)
  pp <- p_of(f)
  data.frame(set = k, par = c("alpha", "gamma", "b0", "fert", "lgdp"), est = pp$est,
             RC = 100 * (full$est - pp$est) / full$est, p = pp$p,
             llQ = as.numeric(logLik(f)), llK = as.numeric(logLik(fk))) }))
print(RES$tab5, digits = 4)
