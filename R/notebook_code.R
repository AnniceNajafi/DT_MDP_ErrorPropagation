#' notebook_code.R
#' Reads a figure script for display and execution in a notebook, without its
#' progress messages: every cat() call and the banner comments (# ====,
#' # ---- chunk: ...) are dropped. cat() only prints, so the notebook produces
#' the same figures as the script itself.
notebook_code <- function(path) {
  x <- readLines(path, warn = FALSE)
  pd <- utils::getParseData(parse(text = x, keep.source = TRUE))
  sym <- pd$parent[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "cat"]
  rng <- pd[match(pd$parent[match(sym, pd$id)], pd$id), c("line1", "col1", "line2", "col2")]
  emptied <- logical(length(x))
  for (i in order(rng$line1, rng$col1, decreasing = TRUE)) {
    r <- rng[i, ]
    x[r$line1] <- paste0(substr(x[r$line1], 1, r$col1 - 1),
                         substr(x[r$line2], r$col2 + 1, nchar(x[r$line2])))
    if (r$line2 > r$line1) x[(r$line1 + 1):r$line2] <- NA
    x[r$line1] <- sub(";\\s*$", "", x[r$line1])            # "a; cat(...)" -> "a"
    emptied[r$line1] <- !grepl("\\S", x[r$line1])
  }
  banner <- grepl("^\\s*#'?\\s*(={5,}|-{4} chunk:)", x)
  x <- x[!is.na(x) & !emptied & !banner]
  invisible(parse(text = x))                               # still valid R
  x
}
