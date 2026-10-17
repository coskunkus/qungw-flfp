###############################################################################
## publish_github.R  -  publishes the code (R/) and the data (data/) with
## README.md, LICENSE, .gitignore and .gitattributes to GitHub. The results and
## the manuscript stay local; anyone who clones reproduces them with run_all.R.
##
## Before the first run: create an EMPTY repository named qungw-flfp on GitHub
## (no README, no licence) and set GH_USER below. git must be installed and
## able to push (GitHub Desktop or a personal access token).
##
## Run from the R folder:  Rscript publish_github.R
###############################################################################

GH_USER <- "coskunkus"
REPO    <- "qungw-flfp"
BRANCH  <- "main"
## The first run creates the history; later runs clone the repository, replace
## its files with the current ones and push a new commit, so that the history
## of earlier versions is kept. Nothing is force-pushed.
MESSAGE <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "Update code and results"

if (GH_USER == "YOUR_GITHUB_USERNAME") stop("Set GH_USER first.")
ROOT <- normalizePath(file.path(getwd(), ".."))
if (!file.exists("qungw_core.R")) stop("Run this script from the R folder.")

## A folder name with Turkish characters reaches git and, on Windows, even
## list.files() corrupted when it is given as an absolute path. The files are
## therefore listed and copied with paths RELATIVE to the R folder, which are
## ASCII, into an ASCII-only staging folder, and git is always called from
## inside that folder.
## a new folder on every run: on Windows, unlink() cannot remove the read-only
## object files that git leaves in .git, so a reused folder would be corrupt
STAGE <- tempfile(paste0(REPO, "_"))
URL <- sprintf("https://github.com/%s/%s.git", GH_USER, REPO)
old <- getwd()
## everything below runs inside tryCatch, so that the working directory is the
## R folder again afterwards, also when git fails
tryCatch({
## clone from inside tempdir() with a relative target, so that no path has to
## be quoted for the Windows shell; then come back to the R folder for copying
setwd(tempdir())
suppressWarnings(system2("git", c("clone", "-q", "-b", BRANCH, URL, basename(STAGE)), stdout = FALSE, stderr = FALSE))
fresh <- TRUE
if (dir.exists(STAGE)) { setwd(STAGE)
  fresh <- !identical(suppressWarnings(system2("git", c("rev-parse", "--is-inside-work-tree"), stdout = TRUE, stderr = FALSE)), "true") }
setwd(old)
if (fresh) { unlink(STAGE, recursive = TRUE); STAGE <- tempfile(paste0(REPO, "_")); dir.create(STAGE) }
cat(if (fresh) "Creating a new history.\n" else "Adding a commit to the existing history.\n")
## replace the tracked content with the current files
## the code and the data are published: R/, data/ and the repository files
for (d in c("R", "data", "results")) unlink(file.path(STAGE, d), recursive = TRUE)
dir.create(file.path(STAGE, "R")); dir.create(file.path(STAGE, "data"))
copy_set <- function(from, pattern, to) {
  fs <- list.files(from, pattern = pattern, full.names = TRUE)
  if (!length(fs)) stop("no files found in ", from)
  if (!all(file.copy(fs, file.path(STAGE, to), overwrite = TRUE))) stop("copying ", from, " failed")
  cat(sprintf("  %-8s %3d files\n", to, length(fs))) }
copy_set(".", "\\.R$", "R")
copy_set("../data", "\\.(csv|txt)$", "data")
for (f in c("README.md", "LICENSE", ".gitignore", ".gitattributes"))
  if (!file.copy(file.path("..", f), STAGE, overwrite = TRUE)) stop("copying ", f, " failed")

setwd(STAGE)
git <- function(...) { a <- c(...); cat("git", a, "\n")
  st <- system2("git", a); if (st != 0) stop("git ", a[1], " failed") }
if (fresh) { git("init", "-q"); git("checkout", "-q", "-b", BRANCH); git("remote", "add", "origin", URL) }
git("add", "-A")
if (!length(system2("git", c("status", "--porcelain"), stdout = TRUE))) {
  cat("Nothing has changed since the last publication.\n") } else {
  git("commit", "-q", "-m", shQuote(MESSAGE))
  git("push", "-u", "origin", BRANCH) }
cat(sprintf("\nPublished: https://github.com/%s/%s\n", GH_USER, REPO))
}, finally = setwd(old))
