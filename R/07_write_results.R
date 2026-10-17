###############################################################################
## 07_write_results.R
## Writes every number that the manuscript reports into results/:
##   values.tex      one LaTeX macro per number quoted in the text
##   tab1_sim.tex    body of Table 1 (simulation)
##   tab2_desc.tex   body of Table 2 (descriptive statistics)
##   tab3_par.tex    body of Table 3, estimates
##   tab3_comp.tex   body of Table 3, model comparison
##   tab4_quant.tex  body of Table 4 (quantile levels)
##   tab5_impact.tex body of Table 5 (case deletion)
##   tab6_diag.tex   body of Table 6 (residual tests)
## The figures Fig1-Fig4 are written to the same folder by 04_figures.R. The
## manuscript reads all of these files with \input and \includegraphics, so it
## contains no result typed by hand.
## Run by run_all.R; it can also be run on its own from the R folder, since it
## only reads results/*.rds and data/flfp_data.csv.
###############################################################################

cat("\n================ 07_write_results.R ================\n")
if (!exists("ROOT")) ROOT <- normalizePath(file.path(getwd(), ".."))
if (!exists("save_rds")) source("io_utils.R")   # .rds and figure I/O safe for non-ASCII paths
RD  <- file.path(ROOT, "results")
RES <- read_rds(file.path(RD, "results.rds"))
sim_files <- file.path(RD, c("simulation_n100_172_300.rds", "simulation_n500_1000.rds"))
if (!all(file.exists(sim_files)))
  stop("The simulation results are missing from results/. Run the simulation once with\n",
       "  RUN_SIM=1 Rscript run_all.R   (Windows: set RUN_SIM=1, then Rscript run_all.R)")
SIM <- read_rds(file.path(RD, "simulation_n100_172_300.rds"))
SIML <- read_rds(file.path(RD, "simulation_n500_1000.rds"))
ENV <- read_rds(file.path(RD, "envelopes.rds"))
D   <- read.csv(file.path(ROOT, "data", "flfp_data.csv"), stringsAsFactors = FALSE)

f  <- function(x, d) formatC(unname(x), format = "f", digits = d)
fl <- function(x, d) f(floor(unname(x) * 10^d) / 10^d, d)        # rounded down ("at least")
fc <- function(x, d) f(ceiling(unname(x) * 10^d) / 10^d, d)      # rounded up ("below", "at most")
nm <- function(iso) sub(" \\(.*\\)$", "", sub(",.*$", "", D$country[match(iso, D$iso3c)]))  # "Iran, Islamic Rep." -> "Iran", "Puerto Rico (US)" -> "Puerto Rico"
join <- function(x, sep = " and ") if (length(x) < 2) x else paste(paste(x[-length(x)], collapse = ", "), x[length(x)], sep = sep)
word <- function(k) c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve")[k]
cap  <- function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))
pv   <- function(p) if (p < 0.001) "<0.001" else f(p, 3)
sci  <- function(p) { e <- floor(log10(p)); m <- round(p / 10^e); if (m == 10) { m <- 1; e <- e + 1 }
  sprintf("%d\\times10^{%d}", m, e) }
neg  <- function(s) ifelse(substr(s, 1, 1) == "-", paste0("$", s, "$"), s)  # minus sign in text cells

V <- character(0)
M <- function(name, value) { stopifnot(grepl("^v[A-Z][A-Za-z]*$", name), !name %in% names(V)); V[name] <<- as.character(value) }

## ----------------------------------------------------------------- data (Sec. 1, 4.1)
yf <- D[!is.na(D$y), ]; yf <- yf[order(yf$y), ]
lowest <- f(yf$y[1], 2); M("vIntroLow", lowest)
M("vIntroLowNames", join(nm(yf$iso3c[f(yf$y, 2) == lowest])))       # in increasing order of y
M("vIntroHigh", f(max(yf$y), 2)); M("vIntroHighName", nm(yf$iso3c[which.max(yf$y)]))
M("vNresp", RES$n_response); M("vN", RES$n)
smp <- D[complete.cases(D[, c("y", "fert", "lgdppc")]), ]
M("vYmin", f(RES$y_range[1], 4)); M("vYminName", nm(smp$iso3c[which.min(smp$y)]))
M("vYmax", f(RES$y_range[2], 4)); M("vYmaxName", nm(smp$iso3c[which.max(smp$y)]))
M("vVIFsevenLo", f(RES$vif7_range[1], 1)); M("vVIFsevenHi", f(RES$vif7_range[2], 1)); M("vCorSevenMax", f(RES$cor7_max, 2))
M("vCondX", f(RES$kappa_X, 1)); M("vVIF", f(RES$vif, 2)); M("vCorX", f(RES$cor_x, 2))
mi <- RES$miss
M("vNexclWord", word(mi$n_excluded)); M("vMissMinName", nm(mi$iso_min)); M("vMissMin", f(mi$min_excluded, 4))
M("vMissDecileWord", cap(word(mi$in_lowest_decile))); M("vMissExpected", f(mi$expected, 1))
M("vWilcoxY", f(mi$wilcox_y, 2)); M("vWilcoxFert", f(mi$wilcox_fert, 2))
M("vFertFull", f(mi$full["beta_fert"], 2)); M("vFertEst", f(mi$est["beta_fert"], 2))
M("vDAICfull", f(-mi$full["dAIC"], 1)); M("vDAICest", f(-mi$est["dAIC"], 1))
M("vFertSD", f(RES$desc["fert", "SD"], 2))

## ----------------------------------------------------------------- Sec. 4.2
mb <- RES$median_boundary
M("vRunsAccepted", RES$runs_accepted_025)
t4m <- RES$tab4
M("vLLmedDrop", f(t4m$logLik[t4m$u == 0.25 & t4m$model == "QUNGW"] - t4m$logLik[t4m$u == 0.50 & t4m$model == "QUNGW"], 2))
stopifnot(t4m$interior[t4m$u == 0.50 & t4m$model == "QUNGW"] == 0,                      # boundary at the median
          diff(range(t4m$logLik[t4m$model == "Kumaraswamy"])) < 0.01)                    # "by less than 0.01"
## unit-Weibull limit (Proposition 2) and the gamma floor of the acceptance rules
stopifnot(abs(diff(RES$uw_limit)) < 1e-4)
M("vLLuw", f(RES$uw_limit["uw"], 2)); M("vLLfloor", f(RES$ll_gamma_floor, 2))

## ----------------------------------------------------------------- Sec. 4.3
t3 <- RES$tab3; cp <- RES$comp; t4 <- RES$tab4
M("vPfert", f(t3$p[4], 3)); stopifnot(t3$p[4] >= 0.001)
M("vPfertS", f(t3$p_s[4], 3)); M("vPgdpS", if (t3$p_s[5] < 0.001) "<0.001" else f(t3$p_s[5], 3))
stopifnot(t3$p_s[4] > 0.05, t3$p_s[5] < 0.05)                                       # fertility not, income significant
M("vAlpha", f(t3$est[1], 4)); M("vGamma", f(t3$est[2], 4))
M("vBzero", f(t3$est[3], 4)); M("vBone", f(t3$est[4], 4)); M("vBtwo", f(t3$est[5], 4))
M("vLLQ", f(cp["QUNGW", "logLik"], 2)); M("vLLK", f(cp["Kumaraswamy", "logLik"], 2)); M("vGap", f(RES$gap, 2))
M("vLLQten", f(t4$logLik[t4$u == 0.10 & t4$model == "QUNGW"], 2))
M("vLLQfifty", f(t4$logLik[t4$u == 0.50 & t4$model == "QUNGW"], 2))
M("vLLQseventyfive", f(t4$logLik[t4$u == 0.75 & t4$model == "QUNGW"], 2))
pr <- RES$proc; pe <- function(u, k) pr$est[abs(pr$u - u) < 1e-9 & pr$term == k]
M("vProcFertLo", f(pe(0.05, "fert"), 3)); M("vProcFertHi", f(pe(0.35, "fert"), 3))
M("vProcGdpLo", f(pe(0.05, "lgdppc"), 3)); M("vProcGdpHi", f(pe(0.35, "lgdppc"), 3))
L <- RES$lr$table; ci <- function(i) sprintf("(%s,\\,%s)", f(L[i, "lr_lo"], 3), f(L[i, "lr_hi"], 3))
M("vLRint", ci(1)); M("vLRfert", ci(2)); M("vLRgdp", ci(3)); M("vKappa", f(RES$lr$kappa, 2))
stopifnot(max(abs(L[, c("lr_lo", "lr_hi")] - L[, c("prof_lo", "prof_hi")])) < 1e-10,   # "to within 10^{-10}"
          max(abs(L[, c("lr_lo", "lr_hi")] - L[, c("wald_lo", "wald_hi")])) < 0.01,    # "by less than 0.01"
          all(pr$hi < 0 | pr$lo > 0))                                                  # "exclude zero at every level"

## ----------------------------------------------------------------- Sec. 4.4
M("vVuong", f(RES$vuong, 3)); M("vVuongP", f(RES$vuong_p, 3))
M("vBootKmed", f(RES$boot_kum["median"], 1)); M("vBootQmed", f(RES$boot_q["median"], 1)); M("vBootQpct", f(RES$boot_q["pctile"], 1))
stopifnot(RES$boot_kum["p"] == 0, all(RES$boot_fail == 0))                              # "none of its values reaches", "every refit converged"
M("vBootB", RES$boot_kum["B"]); M("vBootKmax", f(RES$boot_kum["max"], 1)); M("vBootP", f(RES$boot_kum["p_mc"], 3))
M("vTopFiveNames", join(nm(RES$top5$iso))); M("vTopFiveSum", f(RES$top5_sum, 1)); M("vDecileOne", f(RES$decile1, 1))
M("vLooVuongLo", f(RES$loo_vuong[1], 3)); M("vLooVuongLoName", nm(RES$loo_vuong_where[1]))
M("vLooVuongHi", f(RES$loo_vuong[2], 3)); M("vLooVuongHiName", nm(RES$loo_vuong_where[2]))
M("vLooGapLo", f(RES$loo_gap[1], 2)); M("vLooGapHi", f(RES$loo_gap[2], 2))
stopifnot(RES$loo_gap[1] > 0)                                                          # "larger log-likelihood in every subsample"

## ----------------------------------------------------------------- Sec. 4.5
inf <- RES$infl; t5 <- RES$tab5
stopifnot(names(inf$GD_top)[1:2] == names(inf$LD_top)[1:2])
M("vInflFour", join(sort(nm(names(inf$GD_top)))))
M("vInflTopName", nm(names(inf$GD_top)[1])); M("vInflSecondName", nm(names(inf$GD_top)[2]))
M("vGDtop", f(inf$GD_top[1], 3)); M("vLDtop", f(inf$LD_top[1], 3)); M("vGDsecond", f(inf$GD_top[2], 3)); M("vLDsecond", f(inf$LD_top[2], 3))
stopifnot(inf$dmax_case_top[1] == names(inf$GD_top)[1])
M("vCfertTopName", nm(inf$Ci_fert_top)); M("vCgdpTopName", nm(inf$Ci_lgdp_top))
a <- D[D$iso3c == inf$Ci_fert_top, ]; b <- D[D$iso3c == inf$Ci_lgdp_top, ]
M("vCfertTopY", f(a$y, 3)); M("vCfertTopFert", f(a$fert, 2)); M("vCgdpTopY", f(b$y, 3)); M("vCgdpTopGdp", f(b$lgdppc, 2))
del <- t5[t5$set != "None", ]
M("vShapeRCmax", fc(max(abs(del$RC[del$par %in% c("alpha", "gamma")])), 0))
M("vFertRClo", f(min(del$RC[del$par == "fert"]), 1)); M("vFertRChi", f(max(del$RC[del$par == "fert"]), 1))
ns <- del[del$par == "fert" & del$p >= 0.05, ]
M("vFertNsNames", join(ns$set, " or ")); M("vFertNsP", join(f(ns$p, 3)))
stopifnot(all(del$p[del$par == "lgdp"] < 0.05))                                       # "income coefficient remains significant"
dl <- tapply(t5$llQ - t5$llK, t5$set, `[`, 1)
M("vLLDiffLo", f(min(dl), 1)); M("vLLDiffHi", f(max(dl), 1))

## ----------------------------------------------------------------- Sec. 4.6 (Fig. 4)
E <- RES$env
M("vEnvTailN", E["lower_tail_n", 1]); M("vEnvKum", E["lower_tail_below", "Kumaraswamy"])
M("vEnvBeta", E["lower_tail_below", "Beta"]); M("vEnvUW", E["lower_tail_below", "Unit-Weibull"])
stopifnot(E["lower_tail_below", "QUNGW"] == 1); M("vEnvQName", nm(RES$env_q_below_iso))
thq <- qnorm(ppoints(length(ENV[[1]]$r))); band <- thq > -1.6 & thq < 1.6
## range of theoretical quantiles over which EVERY model has residuals above
## (left of the centre) or below (right of the centre) its envelope: the
## intersection of the per-model ranges
up <- sapply(ENV, function(e) range(thq[band & thq < 0 & e$r > e$hi]))
dn <- sapply(ENV, function(e) range(thq[band & thq > 0 & e$r < e$lo]))
up <- c(max(up[1, ]), min(up[2, ])); dn <- c(max(dn[1, ]), min(dn[2, ]))
stopifnot(up[1] < up[2], dn[1] < dn[2])
M("vEnvUpLo", f(up[1], 1)); M("vEnvUpHi", f(up[2], 1)); M("vEnvDnLo", f(dn[1], 1)); M("vEnvDnHi", f(dn[2], 1))
M("vEnvB", RES$env_B); stopifnot(all(RES$env_fail == 0))
M("vDiagPmax", f(max(RES$tab6$KS_pb, RES$tab6$SW_pb), 3)); stopifnot(max(RES$tab6$KS_pb, RES$tab6$SW_pb) < 0.05)

## ----------------------------------------------------------------- Sec. 3 (Table 1)
ss <- function(n, k) SIM[SIM$n == n & SIM$parameter == k, ]
## the simulation design uses the estimates of Table 3 rounded to four decimals
stopifnot(f(SIM$true[match(c("alpha", "gamma", "beta0", "beta1 (fert)", "beta2 (lgdppc)"), SIM$parameter)], 4) == f(t3$est, 4),
          f(SIML$true[1:5], 4) == f(t3$est, 4))
shape <- c("alpha", "gamma"); coefs <- setdiff(unique(SIM$parameter), shape)
M("vSimCovBetaLo", f(min(SIM$coverage[SIM$parameter %in% coefs]), 3)); M("vSimCovBetaHi", f(max(SIM$coverage[SIM$parameter %in% coefs]), 3))
M("vSimBiasAlphaLo", f(ss(100, "alpha")$bias, 3)); M("vSimBiasAlphaHi", f(ss(300, "alpha")$bias, 3))
M("vSimBiasGammaLo", f(ss(100, "gamma")$bias, 3)); M("vSimBiasGammaHi", f(ss(300, "gamma")$bias, 3))
M("vSimRmseAlphaLo", f(ss(100, "alpha")$rmse, 3)); M("vSimRmseAlphaHi", f(ss(300, "alpha")$rmse, 3))
M("vSimCovAlphaSmall", f(ss(100, "alpha")$coverage, 3)); M("vSimCovGammaSmall", f(ss(100, "gamma")$coverage, 3))
big <- SIM$coverage[SIM$n >= 172 & SIM$parameter %in% shape]
M("vSimCovShapeLo", f(min(big), 3)); M("vSimCovShapeHi", f(max(big), 3))
M("vSimLargeCovLo", f(min(SIML$coverage), 3)); M("vSimLargeCovHi", f(max(SIML$coverage), 3))
M("vSimLargeConv", fl(min(SIML$convergence), 3))
cv <- sapply(c(100, 172, 300), function(n) SIM$convergence[SIM$n == n][1])
M("vSimConv", sprintf("%s, %s and %s", f(cv[1], 3), f(cv[2], 3), f(cv[3], 3)))

writeLines(c("%% Generated by R/07_write_results.R -- do not edit by hand.",
             sprintf("\\newcommand{\\%s}{%s}", names(V), V)), file.path(RD, "values.tex"))

## ----------------------------------------------------------------- table bodies
wr <- function(file, lines) writeLines(c("%% Generated by R/07_write_results.R -- do not edit by hand.", lines), file.path(RD, file))
lab <- c(alpha = "$\\alpha$", gamma = "$\\gamma$", beta0 = "$\\beta_0$", "beta1 (fert)" = "$\\beta_1$", "beta2 (lgdppc)" = "$\\beta_2$")

## Table 1
t1 <- unlist(lapply(c(100, 172, 300), function(n) { s <- SIM[SIM$n == n, ]
  c(sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", ifelse(seq_len(nrow(s)) == 1, n, ""), lab[s$parameter],
            neg(f(s$true, 3)), neg(f(s$bias, 3)), f(s$sd, 3), f(s$rmse, 3), f(s$coverage, 3)), "\\addlinespace") }))
wr("tab1_sim.tex", head(t1, -1))

## Table 2
dl2 <- c(y = "FLFP rate ($y$)", fert = "Fertility", lgdppc = "Log GDP per capita")
wr("tab2_desc.tex", sprintf("%s & %s \\\\", dl2[rownames(RES$desc)],
   apply(RES$desc, 1, function(r) paste(f(r, 3), collapse = " & "))))

## Table 3
l3 <- c("$\\alpha$", "$\\gamma$", "Intercept", "Fertility", "Log GDP per capita")
pp3 <- function(p) if (p < 0.001) "$<0.001$" else f(p, 3)
wr("tab3_par.tex", sapply(1:5, function(i) sprintf("%s & %s & %s & %s & %s & $(%s,\\,%s)$ & %s & %s \\\\", l3[i],
   neg(f(t3$est[i], 4)), f(t3$se[i], 4), if (i > 2) neg(f(t3$z[i], 3)) else "--",
   if (i > 2) pp3(t3$p[i]) else "--", f(t3$lo[i], 3), f(t3$hi[i], 3),
   f(t3$se_s[i], 4), if (i > 2) pp3(t3$p_s[i]) else "--")))
wr("tab3_comp.tex", sprintf("%s & %s & %d & $%s$ & $%s$ & \\\\", rownames(cp), f(cp[, "logLik"], 3), as.integer(cp[, "k"]), f(cp[, "AIC"], 2), f(cp[, "BIC"], 2)))

## Table 4
t4l <- unlist(lapply(unique(t4$u), function(u) { s <- t4[t4$u == u, ]
  c(sprintf("%s & %s & %s & %d & %s & %s \\\\", ifelse(seq_len(nrow(s)) == 1, f(u, 2), ""), s$model,
            ifelse(s$interior == 0, sprintf("$%s^{a}$", f(s$logLik, 2)), f(s$logLik, 2)), s$k,
            ifelse(is.na(s$AIC), "--", sprintf("$%s$", f(s$AIC, 2))), ifelse(is.na(s$used), "--", f(s$used, 0))), "\\addlinespace") }))
wr("tab4_quant.tex", head(t4l, -1))

## Table 5
pp <- function(p) if (p < 0.001) "($<0.001$)" else sprintf("(%s)", f(p, 3))
t5l <- unlist(lapply(unique(t5$set), function(s) { r <- t5[t5$set == s, ]
  c(sprintf("%s & %s & %s & %s \\\\", s, paste(neg(f(r$est, 3)), collapse = " & "), f(r$llQ[1], 2), f(r$llK[1], 2)),
    if (s != "None") sprintf("     & %s & & \\\\", paste(sprintf("[%s]", neg(f(r$RC, 1))), collapse = " & ")),
    sprintf("     & -- & -- & %s & & \\\\", paste(sapply(r$p[3:5], pp), collapse = " & ")), "\\addlinespace") }))
wr("tab5_impact.tex", head(t5l, -1))

## Table 6
t6 <- RES$tab6
wr("tab6_diag.tex", sprintf("%s & %s & %s & %s & %s & %s & $%s$ & $%s$ \\\\", t6$model, f(t6$KS, 3), f(t6$KS_pb, 3), f(t6$CS_mean, 3),
   f(t6$SW, 3), f(t6$SW_pb, 3), f(t6$skew, 2), f(t6$min_r, 2)))

cat(sprintf("Written to %s: values.tex (%d macros) and 7 table files\n", RD, length(V)))
