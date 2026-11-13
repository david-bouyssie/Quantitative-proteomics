## Helpers shared by the test scripts (sourced by each tests/test_*.R)

# Repository root (the parent of the tests folder) and package library of the renv project
tests_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 1) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath("tests")
})
repo_root <- normalizePath(file.path(tests_dir, ".."))
renv_library <- file.path(repo_root, "env", "renv", "library", "R-4.2", "x86_64-w64-mingw32")
if (!dir.exists(renv_library)) stop("Package library not found: ", renv_library)
.libPaths(c(renv_library, .libPaths()))

example_dir <- file.path(repo_root, "Example", "Datasets", "QC-DA")

# Results of the checks
test_results <- data.frame(check = character(), passed = logical(), stringsAsFactors = FALSE)

# Records a check: prints PASS or FAIL with the description (and details when it fails)
check <- function(description, condition, details = NULL) {
  passed <- isTRUE(condition)
  cat(if (passed) "PASS " else "FAIL ", description, if (!passed && !is.null(details)) paste0("  [", details, "]") else "", "\n", sep = "")
  test_results[nrow(test_results) + 1, ] <<- list(description, passed)
  invisible(passed)
}

# Returns the error message of an expression, or NA if it runs without error
error_message <- function(expr) tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))

# Prints the summary and stops with a non-zero exit status if a check failed
finish <- function(name) {
  n_fail <- sum(!test_results$passed)
  cat(sprintf("\n%s: %d checks, %d passed, %d failed\n", name, nrow(test_results), sum(test_results$passed), n_fail))
  if (n_fail > 0) quit(status = 1, save = "no")
}

# Renders a report of Modules/Quantitative data analysis in a short temporary folder (Windows paths are limited to 260 characters):
# the knitting folder is <tmp>/k/a/b, so that the report finds the library at ../../Library as in the repository
render_report <- function(rmd, params, output_name) {
  tmp <- tempfile("qp_"); dir.create(file.path(tmp, "k", "a", "b"), recursive = TRUE)
  file.copy(file.path(repo_root, "Library"), file.path(tmp, "k"), recursive = TRUE)
  out_dir <- file.path(tmp, "out"); dir.create(out_dir)
  pandoc <- file.path(R.home(), "bin")
  if (file.exists(file.path(pandoc, "pandoc.exe")) || file.exists(file.path(pandoc, "pandoc"))) Sys.setenv(RSTUDIO_PANDOC = pandoc)
  # A copy of the module folder, so that the intermediate files of rmarkdown are not written into the repository
  module_dir <- file.path(tmp, "module"); dir.create(module_dir)
  file.copy(file.path(repo_root, "Modules", "Quantitative data analysis", rmd), module_dir)
  html <- file.path(out_dir, output_name)
  error <- error_message(suppressMessages(rmarkdown::render(file.path(module_dir, rmd), output_file = html, params = params, envir = new.env(),
                                                              quiet = TRUE, knit_root_dir = file.path(tmp, "k", "a", "b"))))
  list(error = error, html = if (file.exists(html)) paste(readLines(html, warn = FALSE, encoding = "UTF-8"), collapse = "\n") else "", dir = out_dir)
}

# Returns the text of a report without its HTML tags and embedded images
report_text <- function(html) gsub("<[^>]+>", " ", gsub("data:image/[^\"']+", "", html))
