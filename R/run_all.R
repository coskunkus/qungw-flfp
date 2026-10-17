###############################################################################
## run_all.R  -  reproduces every number, table and figure of the manuscript
##
##   cd R ; Rscript run_all.R            (a few minutes without simulation)
##   RUN_SIM=1 Rscript run_all.R         (adds the simulation, about 8 minutes)
##
## Folders: R/ (code), data/flfp_data.csv (input), results/ (all results:
## the numbers, tables and figures that the manuscript reads; see 07_write_results.R).
## The data file is produced by get_data.R from the World Bank API; it is kept
## so that later revisions of the database cannot change the results.
###############################################################################

ROOT <- normalizePath(file.path(getwd(), ".."))
t0 <- Sys.time()
RES <- list()
dir.create(file.path(ROOT, "results"), showWarnings = FALSE)
sink(file.path(ROOT, "results", "report.txt"), split = TRUE)
source("qungw_core.R")
source("01_estimation.R")
source("06_lr_intervals.R")
source("02_influence.R")
source("03_residuals.R")
source("04_figures.R")
save_rds(RES, file.path(ROOT, "results", "results.rds"))
if (identical(Sys.getenv("RUN_SIM"), "1")) {
  sys.source("05_simulation.R", envir = new.env(parent = globalenv()))
  e2 <- new.env(parent = globalenv()); assign("NS", c(500, 1000), envir = e2)
  sys.source("05_simulation.R", envir = e2)
}
## values.tex and the table files read by the manuscript
sys.source("07_write_results.R", envir = new.env(parent = globalenv()))
## software versions, for the record
cat("\n", R.version.string, "|", sessionInfo()$running, "\n")
cat("Packages:", paste(sapply(c("maxLik", "numDeriv", "lamW", "Rsolnp"), function(p) paste(p, packageVersion(p))), collapse = ", "), "\n")
cat(sprintf("\nTotal runtime: %.1f minutes\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sink()
## the comparison with the manuscript needs latex/ (not part of the repository)
if (file.exists(file.path(ROOT, "latex", "flfp_ungw_compecon.tex"))) source("check_paper.R")
