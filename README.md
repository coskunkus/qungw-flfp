# qungw-flfp

R code for the paper

> Soydal, H., Acar, M. S., Yorgancılar Atatoprak, F. N., Karakaya, K., & Kuş, C.
> *A New Unit Quantile Regression Based on the Unit NG-Weibull Distribution, with
> an Application to Female Labour Force Participation.* Submitted to
> *Computational Economics*.

## Method in one paragraph

The unit NG-Weibull (UNGW) distribution is obtained from the NG-Weibull
distribution of Wang et al. (2021) through Y = exp(-X); its quantile function is
explicit through the Lambert W function. Re-parameterising it in terms of the
u-th quantile mu gives the QUNGW regression, logit(mu_i) = x_i beta. For
u <= exp(-1) the parameters vary freely, so the model is designed for the lower
quantiles; above exp(-1) the condition gamma(-log mu)^alpha < B(u) of Remark 1
must hold and the maximum likelihood estimate may lie on the boundary. The
parameters are identifiable (Proposition 3), and the unit-Weibull quantile regression
is the limit of the model as gamma -> 0 (Proposition 2). Estimation is by maxLik with a
two-stage multi-start protocol that checks this condition; likelihood ratio
confidence intervals of the coefficients are obtained by the constrained
formulation of Doganaksoy (2021).
Global (generalised Cook's distance, likelihood distance) and local (Cook, 1986)
influence diagnostics, sandwich (White, 1982) standard errors, quantile residuals with simulated envelopes and
bootstrap-calibrated KS and SW tests, and a Monte Carlo study complete the package.

## Files

| File | Content |
|---|---|
| `R/qungw_core.R` | shared core: distribution functions, log-likelihood, estimation protocol, competitors |
| `R/01_estimation.R` | Tables 2-4, quantile process, Vuong test, bootstrap, leave-one-out, missing-data check |
| `R/02_influence.R` | Fig. 3 and Table 5: case deletion and local influence |
| `R/03_residuals.R` | Table 6 and the simulated envelopes of Fig. 4 |
| `R/06_lr_intervals.R` | LR confidence intervals of the regression coefficients (Doganaksoy, 2021) and the profile log-likelihoods of Fig. 2 |
| `R/04_figures.R` | Figs 1-4 (pdf, eps, png) |
| `R/07_write_results.R` | writes `results/values.tex` (every number quoted in the text, as LaTeX macros) and the table bodies `results/tab*.tex` |
| `R/05_simulation.R` | Monte Carlo study (Table 1), 2000 replications |
| `R/run_all.R` | runs everything and writes `results/report.txt` |
| `R/check_paper.R` | substitutes `results/` into the manuscript source and compares every printed number with the recomputed value (needs the manuscript source in `latex/`, not included) |
| `R/io_utils.R` | `.rds` and figure input/output that also works when the folder path contains non-ASCII characters (Windows) |
| `R/get_data.R` | rebuilds `data/flfp_data.csv` from the World Bank API (not needed to reproduce) |
| `R/publish_github.R` | publishes this repository |
| `data/flfp_data.csv` | World Development Indicators, 2023: 184 economies, of which 172 report both covariates |
| `data/data_provenance.txt` | extraction date (20 August 2026), API endpoints and coverage of every series |
| `results/` | everything the manuscript reads: `values.tex`, `tab1_sim.tex`-`tab6_diag.tex`, `Fig1`-`Fig4`, plus the intermediate `.rds` files and the run reports |

## Run order

```sh
cd R
Rscript run_all.R             # tables and figures, a few minutes
RUN_SIM=1 Rscript run_all.R   # adds the simulation study
```

All results are written to `results/`. The manuscript contains no number typed
by hand: it reads `results/values.tex` and the table bodies with `\input`, and the
figures from the same folder, so recompiling it after `run_all.R` reproduces the
paper from the code. Without `RUN_SIM=1`, the simulation results already in
`results/` are used for Table 1.

## Packages

`maxLik`, `numDeriv`, `lamW`, `Rsolnp` (installed automatically if missing); `parallel`
(base R) for the simulation.

## Seeds

Bootstrap and simulation: 20260815. In the simulation, replication r of sample size n has
its own L'Ecuyer-CMRG stream derived from `SEED + n`, so that the results do not depend on
the number of cores or on the platform.
Envelopes: 20261018.

## Notation shared with the paper

| Code | Paper |
|---|---|
| `TAU` | quantile level u (0.25) |
| `logAlpha`, `logGamma` | log alpha, log gamma |
| `mu` | conditional quantile mu_i |
| `feasible_bound(u)` | B(u) in Remark 1 |
| `mu_min_attainable()` | mu_min in eq. (12) |
| `GD`, `LD` | GD_i, LD_i |
| `Cmax`, `dmax`, `Ci` | C_max, d_max, C_i |

## Licence

MIT
