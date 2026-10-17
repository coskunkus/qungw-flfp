###############################################################################
## io_utils.R  -  file input/output that works when the project folder has
## non-ASCII characters in its path (e.g. "G:/Drive'ım/yeni-çalışmalar/...").
## On Windows, saveRDS()/readRDS() with a file name go through gzfile(), and
## the graphics devices through the C library, which fail with "Illegal byte
## sequence" on such paths. file() connections and file.copy() do not, so
##   - .rds files are written and read through gzcon(file(...)), and
##   - figures are drawn into tempdir() and then copied to their place.
## The files produced are identical to those of saveRDS()/the devices.
###############################################################################

save_rds <- function(x, path) {
  con <- gzcon(file(path, "wb")); on.exit(close(con)); saveRDS(x, con); invisible(path) }

read_rds <- function(path) {
  con <- gzcon(file(path, "rb")); on.exit(close(con)); readRDS(con) }

## device(file, ...) is opened on a temporary file, draw() is called, and the
## result is copied to path
save_fig <- function(path, device, draw, ...) {
  tmp <- tempfile(fileext = paste0(".", tools::file_ext(path)))
  device(tmp, ...); on.exit(if (file.exists(tmp)) unlink(tmp))
  tryCatch(draw(), finally = dev.off())
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("could not write ", path)
  invisible(path) }
