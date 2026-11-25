#' download_tep_data.R
#' ----------------------------------------------------------------------
#' Downloads the five Tennessee Eastman Process files used by
#' tep/TEP_Pipeline.R (Russell, Chiang and Braatz benchmark data) into
#' tep/TEP_data/ (or into $DTMDP_TEP_DATA) and verifies their MD5 sums.
#'
#' The data are not redistributed in this repository. They are the
#' Braatz-group TEP simulation files (d00.dat = normal operation, training;
#' dXX_te.dat = test runs with fault IDV(XX) introduced at sample 161),
#' see L.H. Chiang, E.L. Russell, R.D. Braatz, "Fault Detection and
#' Diagnosis in Industrial Systems", Springer, 2001. The default source is
#' a public GitHub mirror of the original files; any copy with identical
#' MD5 sums can be used (place the files in tep/TEP_data/ by hand, or set
#' the TEP_BASE_URL environment variable to another mirror).
#'
#' Usage (repository root):  Rscript tep/download_tep_data.R
#' ----------------------------------------------------------------------

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })

base_url <- Sys.getenv("TEP_BASE_URL",
  "https://raw.githubusercontent.com/camaramm/tennessee-eastman-profBraatz/master")
dest_dir <- Sys.getenv("DTMDP_TEP_DATA", repo_path("tep", "TEP_data"))
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

# MD5 sums of the files used for the paper
md5 <- c(
  "d00.dat"    = "cd5043102a43f5917c3b7a853d16071b",   # normal operation, training (52 x 500, transposed)
  "d00_te.dat" = "0f5c530464b358b3b69b117a95d50081",   # normal operation, test (960 x 52)
  "d01_te.dat" = "b1511cd725076fe9e4b8f2d7093781e5",   # IDV(1)  A/C feed ratio step   -> DynamicsOff
  "d11_te.dat" = "5cb55340bbd61b771ee1d4385e85f3c9",   # IDV(11) reactor cooling water -> SensorNoisy
  "d13_te.dat" = "af68cf12e2a504315f8fa383e86c6b55"    # IDV(13) slow kinetics drift   -> Drift
)

ok <- TRUE
for (f in names(md5)) {
  dest <- file.path(dest_dir, f)
  if (file.exists(dest) && unname(tools::md5sum(dest)) == md5[[f]]) {
    cat(sprintf("  %-11s present, MD5 ok\n", f)); next
  }
  url <- paste0(base_url, "/", f)
  cat(sprintf("  %-11s downloading %s\n", f, url))
  res <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE),
                  error = function(e) { message("    download failed: ", conditionMessage(e)); 1L })
  if (!identical(as.integer(res), 0L) || !file.exists(dest)) { ok <- FALSE; next }
  got <- unname(tools::md5sum(dest))
  if (got != md5[[f]]) {
    message(sprintf("    MD5 mismatch for %s (got %s, expected %s)", f, got, md5[[f]]))
    ok <- FALSE
  }
}
if (ok) {
  cat("TEP data ready in", dest_dir, "\n")
} else {
  message("Some TEP files could not be downloaded or verified. Place the files listed ",
          "above in ", dest_dir, " manually (see README.md).")
  quit(save = "no", status = 1)
}
