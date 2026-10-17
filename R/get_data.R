###############################################################################
##  get_data.R : builds the estimation sample and writes flfp_data.csv
##  ---------------------------------------------------------------------------
##  Run this ONCE. After that qungw_maxLik.R reads the CSV and never touches
##  the network.
##
##  WHY THE OLD VERSION WAS SLOW
##  ----------------------------
##  1. WDI(..., extra = TRUE) downloads and merges the full country metadata
##     table on every call. That is the main cost and it is unnecessary.
##  2. WDI() defaults to per_page = 32500 and then walks the pagination, with
##     retries on every timeout. For a single year the whole series fits in
##     one page of 400 records, so the pagination is wasted work.
##  3. WDI() also validates indicator codes against a cached catalogue, which
##     triggers a further download the first time it runs in a session.
##
##  This script calls the API directly: one request per indicator, 400 records
##  per request, no metadata download, no catalogue validation, no pagination.
##  Eight requests of roughly 50 KB each. On a normal connection it finishes in
##  a few seconds.
##
##  The API throttles rapid successive requests and reports the refusal as
##  HTTP 400 Bad Request rather than 429, which looks like an invalid
##  indicator code but is not. Each request is therefore retried with an
##  exponential backoff, and every series that arrives is cached under
##  wdi_cache/ so that re-running the script only fetches what is missing.
##
##  IF THE API IS UNREACHABLE (proxy, firewall, university network)
##  --------------------------------------------------------------
##  Set USE_API <- FALSE below and download the eight series by hand from
##      https://databank.worldbank.org/source/world-development-indicators
##  choosing all countries, the eight series codes listed in `codes`, and the
##  single year, then export as CSV to wdi_manual.csv. The script reshapes it.
###############################################################################

## ---- paths -----------------------------------------------------------------
## The CSV is written to data/ of the project, which is where qungw_core.R
## reads it. Run this script from the R folder.
if (!exists("ROOT")) ROOT <- normalizePath(file.path(getwd(), ".."))
WORKDIR <- file.path(ROOT, "data")

if (!dir.exists(WORKDIR)) dir.create(WORKDIR, recursive = TRUE)
## Do NOT setwd() here. These scripts used to switch the working directory to
## WORKDIR, which meant that a later source("other_script.R") with a relative
## path resolved against the data folder instead of the code folder and ran a
## stale copy. Every path below is absolute, built from WORKDIR, so the working
## directory is left exactly as the caller set it.
cat("Working directory:", getwd(), "\n\n")

USE_API <- TRUE
YEAR    <- 2023
OUTFILE <- file.path(WORKDIR, "flfp_data.csv")

codes <- c(
  flfp     = "SL.TLF.CACT.FE.ZS",   # DV: female LFP rate, % of female pop 15+
  fert     = "SP.DYN.TFRT.IN",      # total fertility rate
  urb      = "SP.URB.TOTL.IN.ZS",   # urban population, % of total
  dep      = "SP.POP.DPND",         # age dependency ratio
  gdppc    = "NY.GDP.PCAP.PP.KD",   # GDP per capita, PPP, constant 2021 int$
  serv_f   = "SL.SRV.EMPL.FE.ZS",   # employment in services, female
  internet = "IT.NET.USER.ZS",      # individuals using the internet, %
  educ_f   = "SE.SEC.ENRR.FE"       # secondary enrolment, female, gross
)

## Regional and income aggregates, hardcoded so that no metadata request is
## needed. Both the legacy codes (HIC, LIC, LMC, UMC) and the current ones
## (XD, XM, XN, XT) are listed.
AGGREGATES <- c(
  "AFE","AFW","ARB","CSS","CEB","EAR","EAS","EAP","TEA","EMU","ECS","ECA",
  "TEC","EUU","FCS","HPC","HIC","XD","IBD","IBT","IDB","IDX","IDA","INX",
  "LTE","LCN","LAC","TLA","LDC","LMY","LIC","XM","LMC","XN","MEA","MNA",
  "TMN","MIC","NAC","OED","OSS","PSS","PST","PRE","SST","SAS","TSA","SSF",
  "SSA","TSS","UMC","XT","WLD"
)

###############################################################################
## 1. Fetch
###############################################################################

if (USE_API) {

  if (!requireNamespace("jsonlite", quietly = TRUE))
    install.packages("jsonlite", repos = "https://cloud.r-project.org")

  options(timeout = 120)

  ## The World Bank API throttles rapid successive requests and answers with
  ## HTTP 400 rather than 429. The failure is transient: the same URL succeeds
  ## a second later. So each indicator is retried with a backoff, there is a
  ## short pause between indicators, and every successful response is cached
  ## to disk so that a re-run only fetches what is still missing.

  CACHE_DIR <- file.path(WORKDIR, "wdi_cache")
  if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR)

  fetch_one <- function(code, tries = 5) {
    cf <- file.path(CACHE_DIR, paste0(code, "_", YEAR, ".csv"))
    if (file.exists(cf)) {
      d <- read.csv(cf, stringsAsFactors = FALSE)
      cat(sprintf("  %-20s %3d records  (cached)\n", code, nrow(d)))
      return(d)
    }

    url <- sprintf(
      "https://api.worldbank.org/v2/country/all/indicator/%s?date=%d&format=json&per_page=400",
      code, YEAR)

    for (k in seq_len(tries)) {
      t0 <- Sys.time()
      js <- try(suppressWarnings(jsonlite::fromJSON(url, flatten = TRUE)),
                silent = TRUE)
      ok <- !inherits(js, "try-error") && length(js) >= 2 &&
            !is.null(js[[2]]) && is.data.frame(js[[2]]) && nrow(js[[2]]) > 0
      if (ok) {
        d <- js[[2]]
        out <- data.frame(iso3c   = d$countryiso3code,
                          country = d$country.value,
                          value   = d$value,
                          stringsAsFactors = FALSE)
        write.csv(out, cf, row.names = FALSE)
        cat(sprintf("  %-20s %3d records  %.1fs%s\n", code, nrow(out),
                    as.numeric(difftime(Sys.time(), t0, units = "secs")),
                    if (k > 1) sprintf("  (attempt %d)", k) else ""))
        return(out)
      }
      if (k < tries) {
        wait <- 2^(k - 1)          # 1, 2, 4, 8 seconds
        cat(sprintf("  %-20s attempt %d failed, retrying in %ds\n",
                    code, k, wait))
        Sys.sleep(wait)
      }
    }
    cat(sprintf("  %-20s FAILED after %d attempts\n", code, tries))
    NULL
  }

  cat("Fetching", length(codes), "indicators for", YEAR, "\n")
  t_all <- Sys.time()
  lst <- list()
  for (nm in names(codes)) {
    lst[[nm]] <- fetch_one(codes[[nm]])
    Sys.sleep(0.5)                 # be polite to the API
  }
  cat(sprintf("Total: %.1fs\n",
              as.numeric(difftime(Sys.time(), t_all, units = "secs"))))

  failed <- names(lst)[vapply(lst, is.null, logical(1))]
  if (length(failed))
    stop("These indicators could not be fetched: ",
         paste(codes[failed], collapse = ", "),
         ".\n  Everything that succeeded is cached in ", CACHE_DIR,
         "/, so simply run this script again and only the missing series ",
         "will be requested. If they keep failing, set USE_API <- FALSE and ",
         "use the manual DataBank export.")

  ## merge the eight series on iso3c
  wide <- lst[[1]][, c("iso3c", "country")]
  for (nm in names(lst)) {
    tmp <- lst[[nm]][, c("iso3c", "value")]
    names(tmp)[2] <- nm
    wide <- merge(wide, tmp, by = "iso3c", all = TRUE)
  }

} else {

  ## Manual DataBank export: rows are country x series, columns are years.
  raw <- read.csv(file.path(WORKDIR, "wdi_manual.csv"), stringsAsFactors = FALSE,
                  na.strings = c("", "..", "NA"))
  ycol <- grep(as.character(YEAR), names(raw), value = TRUE)[1]
  raw$value <- suppressWarnings(as.numeric(raw[[ycol]]))
  key <- setNames(names(codes), codes)
  raw$var <- key[raw$Series.Code]
  raw <- raw[!is.na(raw$var), ]
  wide <- reshape(raw[, c("Country.Code", "Country.Name", "var", "value")],
                  idvar = c("Country.Code", "Country.Name"),
                  timevar = "var", direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  names(wide)[1:2] <- c("iso3c", "country")
}

###############################################################################
## 2. Clean
###############################################################################

wide <- wide[!is.na(wide$iso3c) & wide$iso3c != "", ]
wide <- wide[!wide$iso3c %in% AGGREGATES, ]

dat <- data.frame(
  country  = wide$country,
  iso3c    = wide$iso3c,
  y        = wide$flfp / 100,          # to a proportion in (0,1)
  fert     = wide$fert,
  lgdppc   = log(wide$gdppc),
  urb      = wide$urb / 100,
  dep      = wide$dep / 100,
  serv_f   = wide$serv_f / 100,
  internet = wide$internet / 100,
  educ_f   = wide$educ_f / 100,
  stringsAsFactors = FALSE
)

## Which covariate costs the most observations? Reported so that the trade-off
## between sample size and covariate coverage is explicit rather than hidden
## inside na.omit().
n_before <- nrow(dat)
miss <- sort(colSums(is.na(dat[, -(1:2)])), decreasing = TRUE)
cat("\nMissing values by variable (out of", n_before, "economies):\n")
print(miss[miss > 0])

drop1 <- sapply(names(miss), function(v) {
  keep <- setdiff(names(dat)[-(1:2)], v)
  sum(complete.cases(dat[, keep]))
})
cat("\nSample size if each variable were dropped:\n")
print(sort(drop1, decreasing = TRUE))

## IMPORTANT. The CSV is written WITH the missing values left in place.
## Deleting rows here would apply listwise deletion across all eight columns,
## so a later analysis that uses only three covariates would still be paying
## for missingness in the five it does not use. The analysis script deletes on
## the columns it actually needs. Only rows with a missing RESPONSE are
## dropped, since those are useless to every specification.
dat <- dat[!is.na(dat$y), ]

stopifnot(all(dat$y > 0 & dat$y < 1))

cat(sprintf("\n%d economies with an observed response.\n", nrow(dat)))
cat(sprintf("y in [%.4f, %.4f], no observation on the boundary.\n",
            min(dat$y), max(dat$y)))
cat(sprintf("Complete cases on all eight variables: %d\n",
            sum(complete.cases(dat))))
cat("Sample size for some candidate covariate sets:\n")
for (set in list(c("fert", "lgdppc", "serv_f"),
                 c("fert", "lgdppc"),
                 c("fert", "lgdppc", "serv_f", "educ_f"),
                 c("fert", "lgdppc", "urb", "dep", "serv_f", "internet", "educ_f")))
  cat(sprintf("  %-55s n = %d\n", paste(set, collapse = " + "),
              sum(complete.cases(dat[, c("y", set)]))))

write.csv(dat, OUTFILE, row.names = FALSE)

## ---------------------------------------------------------------------------
## Data provenance record.
## The World Bank Indicators API serves a LIVE database: series are revised, so
## the same request made on two different days can return different numbers.
## A URL alone is therefore not a reproducible data citation. This block writes
## the extraction date, the exact endpoints and the per-indicator coverage to
## data_provenance.txt. Cite that date in the paper and deposit this file with
## the replication material.
## ---------------------------------------------------------------------------
PROV <- file.path(WORKDIR, "data_provenance.txt")
prov <- file(PROV, open = "wt")
writeLines(c(
  "DATA PROVENANCE RECORD",
  "----------------------------------------------------------------",
  paste("Extracted on          :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Reference year        :", YEAR),
  paste("Source                : World Bank World Development Indicators"),
  paste("Response series origin: ILO modelled estimates, ILOSTAT"),
  paste0("API endpoint pattern  : ",
         "https://api.worldbank.org/v2/country/all/indicator/<CODE>",
         "?date=", YEAR, "&format=json&per_page=400"),
  paste("Licence               : WDI is released under CC BY 4.0"),
  paste("R version             :", R.version.string),
  paste("jsonlite version      :",
        as.character(utils::packageVersion("jsonlite"))),
  "",
  "Indicators requested:"), prov)
## Two indicators are renamed when the analysis frame is built (flfp becomes the
## proportion y, gdppc becomes lgdppc on the log scale), so the counts have to be
## looked up under the analysis names, not the download names.
COLMAP <- c(flfp = "y", gdppc = "lgdppc")
## When a cached response is reused the numbers come from the day the cache file
## was written, not from today. Stamping the run time alone would misstate the
## vintage, so each indicator reports its own cache date and whether this run
## fetched it live.
CDIR <- if (exists("CACHE_DIR")) CACHE_DIR else file.path(WORKDIR, "wdi_cache")
cache_when <- function(code) {
  cf <- file.path(CDIR, paste0(code, "_", YEAR, ".csv"))
  if (file.exists(cf)) format(file.info(cf)$mtime, "%Y-%m-%d") else "live fetch"
}
for (nm in names(codes)) {
  col <- if (nm %in% names(COLMAP)) COLMAP[[nm]] else nm
  writeLines(sprintf("  %-9s %-20s -> %-8s n non-missing = %3d   obtained %s",
                     nm, codes[[nm]], col, sum(!is.na(dat[[col]])),
                     cache_when(codes[[nm]])), prov)
}
vint <- vapply(unname(unlist(codes)), cache_when, character(1))
vint <- vint[vint != "live fetch"]
if (length(vint))
  writeLines(c("", sprintf("DATA VINTAGE (cache dates): %s to %s",
                           min(vint), max(vint)),
               "This, not the run date above, is when the numbers were obtained",
               "from the World Bank. Cite it in the paper."), prov)
n_est <- sum(complete.cases(dat[, c("y", "fert", "lgdppc")]))
writeLines(c(
  "",
  sprintf("Economies with an observed response      : %d", nrow(dat)),
  sprintf("ESTIMATION SAMPLE, y + fert + lgdppc     : %d   <-- n reported in the paper",
          n_est),
  sprintf("Complete cases on all %d variables        : %d   (the wider set, rejected)",
          length(codes), sum(complete.cases(dat))),
  "",
  "Aggregates (regional and income groupings) were removed before analysis.",
  "----------------------------------------------------------------"), prov)
close(prov)
cat("Provenance record written to", PROV, "\n")
cat("Written to", OUTFILE, "\n")
