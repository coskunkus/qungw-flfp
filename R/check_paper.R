###############################################################################
## check_paper.R  -  compares every number printed in the manuscript with the
## value recomputed by run_all.R. The manuscript reads its numbers from
## results/values.tex and its table bodies from results/tab*.tex; these are
## substituted into the LaTeX source first, so that the check is made on the
## text as it is typeset. Each check formats the computed value with the
## rounding used in the paper and looks for it in that text. The script also
## checks that the manuscript uses every macro of values.tex and no other.
###############################################################################

if (!exists("ROOT")) ROOT <- normalizePath(file.path(getwd(), ".."))
if (!exists("save_rds")) source("io_utils.R")   # .rds and figure I/O safe for non-ASCII paths
RD  <- file.path(ROOT, "results")
RES <- read_rds(file.path(RD, "results.rds"))
TEX <- paste(readLines(file.path(ROOT, "latex", "flfp_ungw_compecon.tex"), warn = FALSE), collapse = " ")
## table bodies (\restab{...})
for (tb in regmatches(TEX, gregexpr("\\\\restab\\{[A-Za-z0-9_]+\\}", TEX))[[1]]) {
  fn <- sub("\\\\restab\\{(.*)\\}", "\\1", tb)
  TEX <- sub(tb, paste(readLines(file.path(RD, paste0(fn, ".tex")), warn = FALSE)[-1], collapse = " "), TEX, fixed = TRUE) }
## macros of values.tex (\newcommand{\vName}{value})
VAL <- readLines(file.path(RD, "values.tex"), warn = FALSE)[-1]
vn  <- sub("^\\\\newcommand\\{\\\\(v[A-Za-z]+)\\}.*$", "\\1", VAL)
vv  <- sub("^\\\\newcommand\\{\\\\v[A-Za-z]+\\}\\{(.*)\\}$", "\\1", VAL)
used <- unique(sub("^\\\\", "", regmatches(TEX, gregexpr("\\\\v[A-Z][A-Za-z]*", TEX))[[1]]))
undefined <- setdiff(used, vn); unused <- setdiff(vn, used)
for (k in order(-nchar(vn)))
  TEX <- gsub(paste0("\\\\", vn[k], "(?![A-Za-z])(\\{\\}| *)"), gsub("\\", "\\\\", vv[k], fixed = TRUE), TEX, perl = TRUE)
TEX <- gsub("\\s+", " ", TEX)
TEX <- gsub("$", "", TEX, fixed = TRUE)            # compare numbers without math delimiters
f <- function(x, d) formatC(x, format = "f", digits = d)
m <- function(x, d) f(x, d)
out <- list(); chk <- function(what, s) out[[length(out) + 1]] <<- data.frame(item = what, expected = s, found = grepl(s, TEX, fixed = TRUE))

## data and Table 2
chk("economies with response", as.character(RES$n_response)); chk("n", as.character(RES$n))
chk("min y", f(RES$y_range[1], 4)); chk("max y", f(RES$y_range[2], 4))
for (v in rownames(RES$desc)) for (s in colnames(RES$desc)) chk(paste("Table 2", v, s), f(RES$desc[v, s], 3))
chk("condition number", f(RES$kappa_X, 1)); chk("VIF", f(RES$vif, 2)); chk("cor(fert,lgdp)", m(RES$cor_x, 2))
chk("VIF7 low", f(RES$vif7_range[1], 1)); chk("VIF7 high", f(RES$vif7_range[2], 1)); chk("max |cor| 7", f(RES$cor7_max, 2))
## Table 3
t3 <- RES$tab3
for (i in 1:5) { chk(paste("T3 est", rownames(t3)[i]), m(t3$est[i], 4)); chk(paste("T3 se", rownames(t3)[i]), f(t3$se[i], 4))
  chk(paste("T3 CI lo", rownames(t3)[i]), m(t3$lo[i], 3)); chk(paste("T3 CI hi", rownames(t3)[i]), f(t3$hi[i], 3)) }
for (i in 3:5) chk(paste("T3 Wald", rownames(t3)[i]), m(t3$z[i], 3))
chk("T3 p fert", f(t3$p[4], 3))
lab3 <- c("\\alpha", "\\gamma", "Intercept", "Fertility", "Log GDP per capita")
for (i in 1:5) { cells <- c(f(t3$est[i], 4), f(t3$se[i], 4), if (i > 2) f(t3$z[i], 3) else "--",
    if (i > 2) (if (t3$p[i] < 0.001) "<0.001" else f(t3$p[i], 3)) else "--", sprintf("(%s,\\,%s)", f(t3$lo[i], 3), f(t3$hi[i], 3)),
    f(t3$se_s[i], 4), if (i > 2) (if (t3$p_s[i] < 0.001) "<0.001" else f(t3$p_s[i], 3)) else "--")
  chk(paste("T3 row", lab3[i]), paste(c(lab3[i], cells), collapse = " & ")) }
for (k in rownames(RES$comp)) { chk(paste("T3 logLik", k), f(RES$comp[k, "logLik"], 3)); chk(paste("T3 AIC", k), m(RES$comp[k, "AIC"], 2)); chk(paste("T3 BIC", k), m(RES$comp[k, "BIC"], 2)) }
## LR confidence intervals (Section 4.3, Fig. 2)
if (!is.null(RES$lr)) { L <- RES$lr$table
  for (i in 1:3) chk(paste("LR CI", rownames(L)[i]), sprintf("(%s,\\,%s)", f(L[i, "lr_lo"], 3), f(L[i, "lr_hi"], 3)))
  chk("kappa", paste0("\\kappa=", f(RES$lr$kappa, 2)))
  chk("LR vs Wald < 0.01", if (max(abs(L[, c("lr_lo", "lr_hi")] - L[, c("wald_lo", "wald_hi")])) < 0.01) "by less than 0.01" else "MISMATCH")
  chk("LR vs profile < 1e-10", if (max(abs(L[, c("lr_lo", "lr_hi")] - L[, c("prof_lo", "prof_hi")])) < 1e-10) "to within 10^{-10}" else "MISMATCH")
  chk("LR left of Wald", if (all(L[, "lr_lo"] < L[, "wald_lo"] & L[, "lr_hi"] < L[, "wald_hi"])) "lie slightly to the left" else "MISMATCH") }
## Table 4
t4 <- RES$tab4
for (i in seq_len(nrow(t4))) {
  chk(sprintf("T4 logLik u=%.2f %s", t4$u[i], t4$model[i]), f(t4$logLik[i], 2))
  if (!is.na(t4$AIC[i])) chk(sprintf("T4 AIC u=%.2f %s", t4$u[i], t4$model[i]), m(t4$AIC[i], 2)) }
## Table 4, row by row (log-likelihood, k, AIC)
for (i in seq_len(nrow(t4))) { ll <- f(t4$logLik[i], 2); if (t4$interior[i] == 0) ll <- paste0(ll, "^{a}")
  chk(sprintf("T4 row u=%.2f %s", t4$u[i], t4$model[i]), paste(c(t4$model[i], ll, t4$k[i], if (is.na(t4$AIC[i])) "--" else f(t4$AIC[i], 2)), collapse = " & ")) }
## Table 3 model-comparison rows
for (k in rownames(RES$comp)) chk(paste("T3 row", k), paste(c(k, f(RES$comp[k, "logLik"], 3), RES$comp[k, "k"], f(RES$comp[k, "AIC"], 2), f(RES$comp[k, "BIC"], 2)), collapse = " & "))
t4q <- RES$tab4[RES$tab4$model == "QUNGW", ]
chk("median: drop of the log-likelihood", paste0("is ", f(t4q$logLik[t4q$u == 0.25] - t4q$logLik[t4q$u == 0.50], 2), " below that at u=0.25"))
chk("median: Kumaraswamy invariant", if (diff(range(RES$tab4$logLik[RES$tab4$model == "Kumaraswamy"])) < 0.01) "changes by less than 0.01" else "MISMATCH")
chk("unit-Weibull limit value", paste0("model, ", f(RES$comp["Unit-Weibull", "logLik"], 2), ", is the limiting value"))
chk("profile below gamma floor", paste0("does not exceed ", f(RES$ll_gamma_floor, 2)))
chk("accepted runs at u=0.25", paste(RES$runs_accepted_025, "of the 28"))
## quantile process
pr <- RES$proc
chk("fert u=0.05", m(pr$est[pr$u == 0.05 & pr$term == "fert"], 3)); chk("fert u=0.35", m(pr$est[pr$u == 0.35 & pr$term == "fert"], 3))
chk("lgdp u=0.05", f(pr$est[pr$u == 0.05 & pr$term == "lgdppc"], 3)); chk("lgdp u=0.35", f(pr$est[pr$u == 0.35 & pr$term == "lgdppc"], 3))
chk("bands exclude zero (all)", if (all(pr$hi < 0 | pr$lo > 0)) "exclude zero at every level" else "MISMATCH")
## robustness
chk("gap", f(RES$gap, 2)); chk("Vuong", f(RES$vuong, 3)); chk("Vuong p", f(RES$vuong_p, 3))
chk("top5 sum", f(RES$top5_sum, 1)); chk("decile 1", f(RES$decile1, 1))
chk("LOO Vuong min", f(RES$loo_vuong[1], 3)); chk("LOO Vuong max", f(RES$loo_vuong[2], 3))
chk("LOO gap min", f(RES$loo_gap[1], 2)); chk("LOO gap max", f(RES$loo_gap[2], 2))
chk("boot Kum median", m(RES$boot_kum["median"], 1))
chk("boot Kum: largest replicate below the gap", if (RES$boot_kum["max"] < RES$gap) paste0("largest value, ", f(RES$boot_kum["max"], 1), ", is below the observed value ", f(RES$gap, 2)) else "MISMATCH")
chk("boot Kum MC p-value", paste0("p-value 1/(", RES$boot_kum["B"], "+1) is ", f((RES$boot_kum["r"] + 1) / (RES$boot_kum["B"] + 1), 3)))
chk("boot: all refits converged", if (all(RES$boot_fail == 0)) paste(RES$boot_kum["B"], "replications in all of which both") else "MISMATCH")
chk("boot QUNGW median", f(RES$boot_q["median"], 1)); chk("boot QUNGW pctile", paste0("only ", f(RES$boot_q["pctile"], 1), "\\% of the bootstrap"))
## missing data
chk("Yemen", f(RES$miss$min_excluded, 4)); chk("excluded in lowest decile", paste0(c("One","Two","Three","Four")[RES$miss$in_lowest_decile], " of the twelve"))
chk("expected", f(RES$miss$expected, 1)); chk("Wilcoxon y", f(RES$miss$wilcox_y, 2)); chk("Wilcoxon fert", f(RES$miss$wilcox_fert, 2))
chk("fert full", m(RES$miss$full["beta_fert"], 2)); chk("fert est", m(RES$miss$est["beta_fert"], 2))
chk("dAIC full", f(-RES$miss$full["dAIC"], 1)); chk("dAIC est", f(-RES$miss$est["dAIC"], 1))
## influence and Table 5
chk("GD Iran", f(RES$infl$GD_top[1], 3)); chk("LD Iran", f(RES$infl$LD_top[1], 3))
chk("GD 2nd", f(RES$infl$GD_top[2], 3)); chk("LD 2nd", f(RES$infl$LD_top[2], 3))
t5 <- RES$tab5
for (i in seq_len(nrow(t5))) { chk(paste("T5 est", t5$set[i], t5$par[i]), m(t5$est[i], 3))
  if (t5$set[i] != "None") chk(paste("T5 RC", t5$set[i], t5$par[i]), m(t5$RC[i], 1)) }
pf <- function(p) if (p < 0.001) "(<0.001)" else sprintf("(%s)", f(p, 3))
for (s in unique(t5$set)) { r <- t5[t5$set == s, ]
  chk(paste("T5 row", s), paste(c(s, f(r$est, 3), f(r$llQ[1], 2), f(r$llK[1], 2)), collapse = " & "))
  chk(paste("T5 p row", s), paste(c("--", "--", sapply(r$p[3:5], pf)), collapse = " & "))
  if (s != "None") chk(paste("T5 RC row", s), paste(sprintf("[%s]", f(r$RC, 1)), collapse = " & ")) }
for (s in unique(t5$set)) { chk(paste("T5 llQ", s), f(t5$llQ[t5$set == s][1], 2)); chk(paste("T5 llK", s), f(t5$llK[t5$set == s][1], 2)) }
## Table 6 and envelopes
t6 <- RES$tab6
for (i in 1:4) chk(paste("T6 row", t6$model[i]), paste(c(t6$model[i], f(t6$KS[i], 3), f(t6$KS_pb[i], 3), f(t6$CS_mean[i], 3), f(t6$SW[i], 3), f(t6$SW_pb[i], 3)), collapse = " & "))
chk("T6 bootstrap p max", paste0("no bootstrap p-value exceeds ", f(max(t6$KS_pb, t6$SW_pb), 3)))
for (i in 1:4) { chk(paste("T6 KS", t6$model[i]), f(t6$KS[i], 3)); chk(paste("T6 CS mean", t6$model[i]), f(t6$CS_mean[i], 3))
  chk(paste("T6 SW", t6$model[i]), f(t6$SW[i], 3)); chk(paste("T6 skew", t6$model[i]), m(t6$skew[i], 2)); chk(paste("T6 min", t6$model[i]), m(t6$min_r[i], 2)) }
chk("env lower tail n", paste("Of the", RES$env["lower_tail_n", 1], "residuals"))
chk("env Kum/Beta/UW below", sprintf("%d, %d and %d", RES$env["lower_tail_below", "Kumaraswamy"], RES$env["lower_tail_below", "Beta"], RES$env["lower_tail_below", "Unit-Weibull"]))
chk("env QUNGW below = 1", if (RES$env["lower_tail_below", "QUNGW"] == 1) "only that of Iran" else "MISMATCH")
## simulation (if run)
sf <- file.path(RD, "simulation_n100_172_300.rds")
if (file.exists(sf)) { st <- read_rds(sf)
  for (i in seq_len(nrow(st))) { chk(sprintf("T1 bias n=%d %s", st$n[i], st$parameter[i]), m(st$bias[i], 3))
    chk(sprintf("T1 sd n=%d %s", st$n[i], st$parameter[i]), f(st$sd[i], 3)); chk(sprintf("T1 rmse n=%d %s", st$n[i], st$parameter[i]), f(st$rmse[i], 3))
    chk(sprintf("T1 cov n=%d %s", st$n[i], st$parameter[i]), f(st$coverage[i], 3)) } }

chk("no undefined macro", if (length(undefined)) paste("UNDEFINED:", paste(undefined, collapse = " ")) else "")
chk("every macro of values.tex used", if (length(unused)) paste("UNUSED:", paste(unused, collapse = " ")) else "")
res <- do.call(rbind, out)
cat(sprintf("\nCHECK AGAINST THE MANUSCRIPT: %d of %d values found\n", sum(res$found), nrow(res)))
if (any(!res$found)) { cat("\nNot found in the manuscript:\n"); print(res[!res$found, ], row.names = FALSE) }
write.table(res, file.path(RD, "check_paper.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
