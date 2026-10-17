###############################################################################
## 04_figures.R  -  Fig. 1 (quantile process), Fig. 2 (influence index plots),
## Fig. 3 (quantile residuals with simulated envelopes), as PDF, EPS and PNG.
## Requires 01-03 to have been run in the same session (or their .rds files).
###############################################################################
cat("\n================ 04_figures.R ================\n")
proc <- RES$proc
env  <- read_rds(file.path(OUTDIR, "envelopes.rds"))
a    <- read_rds(file.path(OUTDIR, "influence.rds"))
n <- length(a$GD); idx <- seq_len(n)
pal <- list(blue = "#0057B8", blueL = "#9CC3EE", orange = "#E8590C", orangeL = "#F9B98F",
            green = "#00A651", greenL = "#A6E3BF", purple = "#7B2CBF", purpleL = "#D2B7F0",
            grid = "#ECECEC", ax = "#4A4A4A")      # saturated palette
base_par <- function(...) par(..., las = 1, tcl = -0.25, mgp = c(2.4, 0.6, 0), cex.axis = 0.85,
                              col.axis = pal$ax, col.lab = "#333333", fg = pal$ax, bty = "n")
gridlines <- function() invisible(NULL)        # no grid lines inside the panels
## ---------------- Fig 1 ----------------
fig1 <- function() {
  base_par(mfrow = c(1, 2), mar = c(4, 4.2, 2, 0.8))
  spec <- list(fert = c("Fertility", pal$blue, pal$blueL), lgdppc = c("Log GDP per capita", pal$orange, pal$orangeL))
  for (k in names(spec)) { d <- proc[proc$term == k, ]; s <- spec[[k]]
    yl <- range(c(0, d$lo, d$hi)); yl <- yl + c(-1, 1) * 0.05 * diff(yl)
    plot(d$u, d$est, type = "n", ylim = yl, xlab = expression("Quantile level " * italic(u)),
         ylab = "Coefficient", xaxt = "n", main = s[1], font.main = 1, cex.main = 1, col.main = "#333333")
    axis(1, at = seq(0.05, 0.35, 0.1)); gridlines()
    polygon(c(d$u, rev(d$u)), c(d$lo, rev(d$hi)), col = s[3], border = NA)
    abline(h = 0, lty = 3, col = "#888888"); abline(v = 0.25, lty = 2, col = "#AAAAAA")
    lines(d$u, d$est, lwd = 2, col = s[2]); points(d$u, d$est, pch = 21, bg = "white", col = s[2], cex = 0.85, lwd = 1.3) } }
## ---------------- Fig 2: profile log-likelihoods and LR limits ----------------
lr <- read_rds(file.path(OUTDIR, "lr_intervals.rds"))
fig_lr <- function() {
  base_par(mfrow = c(1, 3), mar = c(5.2, 4.4, 2, 0.8)); par(cex = 1)
  labs <- c("Intercept", "Fertility", "Log GDP per capita")
  xl <- list(expression(beta[0]), expression(beta[1]), expression(beta[2]))
  for (k in 1:3) { d <- lr$profiles[[k]]; lo <- lr$table[k, "lr_lo"]; hi <- lr$table[k, "lr_hi"]
    yl <- c(min(d$l), lr$l_hat + 0.12 * (lr$l_hat - min(d$l)))
    plot(d$b, d$l, type = "n", ylim = yl, xlab = "", ylab = "Profile log-likelihood", xaxt = "n",
         main = labs[k], font.main = 1, cex.main = 1, col.main = "#333333"); gridlines()
    axis(1, at = pretty(d$b, n = 6))
    lines(d$b, d$l, lwd = 1.8, col = "black")
    abline(h = lr$kappa, lty = 2, lwd = 1.3, col = pal$blue)
    text(par("usr")[1], lr$kappa, expression(kappa), adj = c(-0.4, -0.35), cex = 1.05, col = pal$blue)
    segments(c(lo, hi), lr$kappa, c(lo, hi), par("usr")[3], lty = 3, col = "black", lwd = 1.2)
    mtext(sprintf("%.3f", c(lo, hi)), side = 1, at = c(lo, hi), line = 1.9, cex = 0.72, col = "#333333")
    mtext(xl[[k]], side = 1, line = 3.6, cex = 0.8, col = "#333333") } }
## ---------------- Fig 3 ----------------
## one colour per row: case deletion (blue), case-weight perturbation (green),
## covariate perturbation (orange); bars and points below the benchmark are grey
fig3 <- function() {
  base_par(mfrow = c(3, 2), mar = c(3.6, 4.2, 1.8, 0.8)); par(cex = 0.9)
  P <- list(list(a$GD, "(a) Generalised Cook's distance", 2*mean(a$GD), expression(GD[i]), pal$blue),
            list(a$LD, "(b) Likelihood distance", 2*mean(a$LD), expression(LD[i]), pal$blue),
            list(a$LI$case$Ci, "(c) Case-weight perturbation", 2*mean(a$LI$case$Ci), expression(C[i]), pal$green),
            list(a$LI$case$dmax, "(d) Case-weight perturbation", 2/sqrt(n), expression("|"*d[max]*"|"), pal$green),
            list(a$LI$fert$Ci, "(e) Perturbation of fertility", 2*mean(a$LI$fert$Ci), expression(C[i]), pal$orange),
            list(a$LI$lgdppc$Ci, "(f) Perturbation of log GDP per capita", 2*mean(a$LI$lgdppc$Ci), expression(C[i]), pal$orange))
  for (q in P) { v <- q[[1]]
    plot(idx, v, type = "n", xlab = "Index", ylab = q[[4]], main = q[[2]], ylim = c(0, max(v) * 1.18),
         font.main = 1, cex.main = 1, col.main = "#333333"); gridlines()
    hi <- v >= q[[3]]
    segments(idx, 0, idx, v, col = ifelse(hi, q[[5]], "#C8C8C8"), lwd = ifelse(hi, 1.6, 1))
    abline(h = q[[3]], lty = 2, col = "#9A9A9A")
    o <- sort(order(v, decreasing = TRUE)[1:4]); dx <- rep(0, 4)
    for (m in 2:4) if (o[m] - o[m-1] < 12) { dx[m-1] <- dx[m-1] - 5; dx[m] <- dx[m] + 5 }
    points(o, v[o], pch = 19, cex = 0.7, col = ifelse(v[o] >= q[[3]], q[[5]], "#A8A8A8"))
    text(o + dx, v[o], a$lab[o], pos = 3, cex = 0.68, col = "#333333", xpd = TRUE) } }
## ---------------- Fig 4 ----------------
fig4 <- function() {
  base_par(mfrow = c(2, 2), mar = c(3.8, 4.2, 1.8, 0.8)); par(cex = 0.9)
  th <- qnorm(ppoints(n)); yl <- range(unlist(lapply(env, function(e) c(e$r, e$lo, e$hi))))
  for (k in c("QUNGW", "Kumaraswamy", "Beta", "Unit-Weibull")) { e <- env[[k]]
    plot(th, e$r, type = "n", ylim = yl, xlab = "Standard normal quantiles", ylab = "Quantile residuals",
         main = k, font.main = 1, cex.main = 1, col.main = "#333333"); gridlines()
    polygon(c(th, rev(th)), c(e$lo, rev(e$hi)), col = pal$blueL, border = NA)
    lines(th, e$m, lty = 2, col = "#7F7F7F")
    out <- e$r < e$lo | e$r > e$hi
    points(th[!out], e$r[!out], pch = 16, cex = 0.5, col = "black")
    points(th[out], e$r[out], pch = 21, cex = 0.6, col = "black", bg = "white", lwd = 1) } }
out <- function(fn, f, w, h) {
  ## cairo devices embed the fonts, as required by the publisher
  save_fig(paste0(fn, ".pdf"), cairo_pdf, f, width = w, height = h, family = "sans")
  save_fig(paste0(fn, ".eps"), cairo_ps, f, width = w, height = h, family = "sans", fallback_resolution = 600)
  save_fig(paste0(fn, ".png"), png, f, width = w, height = h, units = "in", res = 220) }
out(file.path(OUTDIR, "Fig1"), fig1, 7, 3.4); out(file.path(OUTDIR, "Fig2"), fig_lr, 7, 3.3)
out(file.path(OUTDIR, "Fig3"), fig3, 7, 8.4); out(file.path(OUTDIR, "Fig4"), fig4, 7, 6.2)
