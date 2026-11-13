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

finish("QC report tests")
