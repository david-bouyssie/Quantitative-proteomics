## Tests of the QC report (Modules/Quantitative data analysis/QC.Rmd) on the example dataset
## Usage: Rscript --vanilla tests/test_qc.R (or tests/run_all.bat)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
cat("QC report tests\n\n")

quant <- file.path(example_dir, "Parsed_proteins_set.txt")
design <- file.path(example_dir, "exp_design.txt")
qc_parameters <- file.path(example_dir, "QC_parameters.txt")

# Writes a copy of the QC parameters with some values changed (NULL removes the parameter)
qc_params_with <- function(changes) {
  p <- readLines(qc_parameters, warn = FALSE)
  for (name in names(changes)) {
    p <- p[!startsWith(p, paste0(name, "="))]
    if (!is.null(changes[[name]])) p <- c(p, paste0(name, "=", changes[[name]]))
  }
  f <- tempfile(fileext = ".txt"); writeLines(p, f); f
}
qc_render <- function(quantif_file = quant, design_exp = design, parameters = qc_parameters) {
  render_report("QC.Rmd", list(quantif_file = quantif_file, design_exp = design_exp, parameters = parameters, origin_file = ""), "QC.html")
}

## 1. With the experimental design
r <- qc_render()
check("QC renders with the experimental design", is.na(r$error), r$error)
txt <- report_text(r$html)
check("imputation with the noise model, sd estimated within conditions",
      grepl("imputed with the noise model", txt) && grepl("measurements within conditions", txt))
check("noise level at the 5% percentile of each sample", grepl("the 5% percentile of its observed intensities", txt))
sd_value <- as.numeric(sub(".*standard deviation of ([0-9.]+) log2.*", "\\1", txt))
check("noise sd of the example close to 0.8 log2", isTRUE(abs(sd_value - 0.81) < 0.05), paste("sd =", sd_value))
check("dimension reduction section present", grepl("Dimension reduction", txt))

## 2. Without experimental design
r <- qc_render(design_exp = "")
check("QC renders without experimental design", is.na(r$error), r$error)
txt <- report_text(r$html)
check("without design, the noise sd is estimated over all the samples taken as one group", grepl("over all the samples taken as one group", txt))

## 3. Small dataset: the noise model cannot be estimated, the percentile model is used
small <- tempfile(fileext = ".txt"); writeLines(readLines(quant, n = 31), small)
r <- qc_render(quantif_file = small)
check("QC renders on a 30-protein dataset", is.na(r$error), r$error)
check("small dataset: fallback to the percentile model, with the reason", grepl("replaced by the 5% percentile.*noise model could not be estimated", report_text(r$html)))

## 4. Imputation.MNAR.percentile: absent (5 by default) and invalid
r <- qc_render(parameters = qc_params_with(list(Imputation.MNAR.percentile = NULL)))
check("QC renders without Imputation.MNAR.percentile (5 by default)", is.na(r$error) && grepl("the 5% percentile of its observed intensities", report_text(r$html)), r$error)
r <- qc_render(parameters = qc_params_with(list(Imputation.MNAR.percentile = "20")))
check("Imputation.MNAR.percentile = 20 stops with a clear message", grepl("Imputation.MNAR.percentile must be greater than 0 and at most 10", r$error), r$error)

## 5. Sample names over 15 characters: shortened to S1, S2... without design, to <Condition>_<index> with a design
long_quant <- read.delim(quant, check.names = FALSE, stringsAsFactors = FALSE)
names(long_quant) <- sub("^(Intensity|Identification_type)_(.*)$", "\\1_\\2_long_sample_name", names(long_quant))
long_quant_file <- tempfile(fileext = ".txt")
write.table(long_quant, long_quant_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
r <- qc_render(quantif_file = long_quant_file, design_exp = "")
txt <- report_text(r$html)
check("long sample names without design: QC renders, samples numbered S1 to S40", is.na(r$error) && grepl("10amol_1_long_sample_name\\s+S1\\s", txt) && grepl("\\sS40\\s", txt), r$error)

# Design with long replicate names and its rows shuffled (the conditions are interleaved)
long_design <- read.delim(design, stringsAsFactors = FALSE)
long_design$File <- paste0(long_design$File, "_long_sample_name")
long_design$repBio <- paste0("injection_", long_design$repBio)
set.seed(1); long_design <- long_design[sample(nrow(long_design)), ]
long_design_file <- tempfile(fileext = ".txt")
write.table(long_design, long_design_file, sep = "\t", quote = FALSE, row.names = FALSE)
r <- qc_render(quantif_file = long_quant_file, design_exp = long_design_file)
txt <- report_text(r$html)
short_names <- regmatches(txt, gregexpr("\\b[0-9]+(amol|fmol)_[0-9]+\\b", txt))[[1]]
check("long sample names with a shuffled design: QC renders, names <Condition>_<index> from 1 to 4 in each condition",
      is.na(r$error) && all(paste0(rep(unique(long_design$Condition), each = 4), "_", 1:4) %in% short_names) && !any(grepl("_[5-9]$", short_names)), r$error)

## 6. Launcher (Launchers/run_proteomics_stats.R), as called by the webserver
launcher_dir <- tempfile("qc_launcher_"); dir.create(launcher_dir)
# Original file, then experimental design as last argument
html <- file.path(launcher_dir, "QC.html")
l <- run_launcher("run_proteomics_stats.R", c("QC", quant, html, qc_parameters, quant, design))
txt <- if (file.exists(html)) report_text(paste(readLines(html, warn = FALSE, encoding = "UTF-8"), collapse = "\n")) else ""
check("launcher: QC renders with the original file and the experimental design", l$status == 0 && file.exists(html), tail(l$output, 3))
check("launcher: paths with forward slashes in the report, complete report (pandoc reads backslashes as LaTeX commands)",
      grepl(normalizePath(quant, winslash = "/"), txt, fixed = TRUE) && grepl("Dimension reduction", txt))
check("launcher: experimental design used", nzchar(txt) && !grepl("over all the samples taken as one group", txt))
# No original file ("") and an experimental design file that does not exist (webserver without a valid design): QC without design
html <- file.path(launcher_dir, "QC_no_design.html")
l <- run_launcher("run_proteomics_stats.R", c("QC", quant, html, qc_parameters, "", file.path(launcher_dir, "missing_exp_design.txt")))
txt <- if (file.exists(html)) report_text(paste(readLines(html, warn = FALSE, encoding = "UTF-8"), collapse = "\n")) else ""
check("launcher: QC renders without original file and with a missing experimental design file", l$status == 0 && grepl("over all the samples taken as one group", txt), tail(l$output, 3))
# Failure: non-zero exit status, the R error and the position in the report in the output (read by the webserver for the task message)
l <- run_launcher("run_proteomics_stats.R", c("QC", quant, file.path(launcher_dir, "QC_error.html"), qc_params_with(list(Imputation.MNAR.percentile = "20")), "", design))
check("launcher: a QC error gives a non-zero exit status, the error message and the report position",
      l$status != 0 && any(grepl("^Error", l$output)) && any(grepl("Imputation.MNAR.percentile must be greater than 0", l$output)) &&
        any(grepl("^Quitting from .*QC.Rmd", l$output)),
      paste(tail(l$output, 4), collapse = " | "))

finish("QC report tests")
