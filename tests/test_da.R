## Tests of the DA report (Modules/Quantitative data analysis/DA.rmd) on the example dataset
## The p-values, q-values and significant flags of the output table are recomputed independently from the exported test values
## Usage: Rscript --vanilla tests/test_da.R (or tests/run_all.bat)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
cat("DA report tests\n\n")

quant <- file.path(example_dir, "Parsed_proteins_set.txt")
design <- file.path(example_dir, "exp_design.txt")
da_parameters <- file.path(example_dir, "Parameters.txt")

# Writes a copy of the DA parameters with some values changed (NULL removes the parameter)
da_params_with <- function(changes = list()) {
  p <- readLines(da_parameters, warn = FALSE)
  for (name in names(changes)) {
    p <- p[!startsWith(p, paste0(name, "="))]
    if (!is.null(changes[[name]])) p <- c(p, paste0(name, "=", changes[[name]]))
  }
  f <- tempfile(fileext = ".txt"); writeLines(p, f); f
}
# Renders DA and reads the first sheet of the output table (NULL if the render failed)
da_render <- function(changes = list(), quant_df = quant) {
  r <- render_report("DA.rmd", list(quant_df = quant_df, parameters = da_params_with(changes), design_exp = design, origin_file = ""),
                     "DA.html", output_dir_param = "output_dir")
  xlsx <- file.path(r$dir, "DA_output_table.xlsx")
  r$table <- if (is.na(r$error) && file.exists(xlsx)) openxlsx::read.xlsx(xlsx, sheet = 1, check.names = FALSE, sep.names = " ") else NULL
  r$sheets <- if (file.exists(xlsx)) openxlsx::getSheetNames(xlsx) else character(0)
  r
}
num <- function(x) suppressWarnings(as.numeric(x))
flag <- function(x) x %in% c(TRUE, "TRUE")

# Log intensities used by the test for a condition of a comparison, one column per replicate, in the order of the experimental design
test_values <- function(tab, condition, comparison) {
  reps <- read.delim(design)
  reps <- reps$repBio[reps$Condition == condition]
  as.matrix(sapply(reps, function(rep) num(tab[[paste0(rep, "_", condition, "_mean_imputed_log_", comparison)]])))
}
# Two-group p-values recomputed for every tested protein (NA rows of the export are kept NA)
recompute_pvalues <- function(tab, comparison, test, b_order = NULL) {
  conds <- strsplit(comparison, "/")[[1]]
  A <- test_values(tab, conds[1], comparison); B <- test_values(tab, conds[2], comparison)
  if (!is.null(b_order)) B <- B[, b_order]
  tested <- !is.na(num(tab[[paste0("pval_", comparison)]]))
  p <- rep(NA_real_, nrow(tab))
  p[tested] <- sapply(which(tested), function(i) test(A[i, ], B[i, ]))
  p
}
# Spike-ins (UPS1, expected up in the first condition) found, and background (yeast) hits, of a comparison
score <- function(tab, comparison) {
  sig <- flag(tab[[paste0("significant_", comparison)]]); ups <- grepl("_UPS$", tab$Accession)
  c(spike_ins = sum(sig & ups & num(tab[[paste0("log2_fc_", comparison)]]) > 0), ups = sum(ups & !is.na(num(tab[[paste0("pval_", comparison)]]))),
    yeast = sum(sig & !ups))
}
# Significance rule of the report: q-value at most the threshold and absolute log2 fold change at least the threshold
expected_significant <- function(tab, comparison, q, threshold = 0.05, fc = 0.58) {
  !is.na(q) & q <= threshold & abs(num(tab[[paste0("log2_fc_", comparison)]])) >= fc
}
max_diff <- function(a, b) max(abs(a - b), na.rm = TRUE)

comparisons <- c("50fmol/25fmol", "50fmol/500amol")

## 1. Example parameters: Student t-test on log10 values, BH
r <- da_render()
check("DA renders with the example parameters", is.na(r$error), r$error)
tab <- r$table
check("output workbook with the data, parameters, design and metadata sheets",
      identical(r$sheets, c("Parsed_proteins_set", "R_Params", "Exp_Design", "R_Metadata")), paste(r$sheets, collapse = ", "))
check("one Accession per row in the output table", !is.null(tab) && !anyDuplicated(tab$Accession))
for (cm in comparisons) {
  expected_cols <- paste0(c("pval_", "qval_", "significant_", "log2_fc_", "fc_", "mean_intensity_"), cm)
  check(paste0(cm, ": result columns present"), all(expected_cols %in% names(tab)), paste(setdiff(expected_cols, names(tab)), collapse = ", "))
  p_ref <- recompute_pvalues(tab, cm, function(a, b) t.test(a, b, var.equal = TRUE)$p.value)
  p <- num(tab[[paste0("pval_", cm)]])
  check(paste0(cm, ": p-values equal to an independent Student t-test on the exported log values"), max_diff(p, p_ref) < 1e-10, max_diff(p, p_ref))
  q_ref <- p; q_ref[!is.na(p)] <- p.adjust(p[!is.na(p)], "BH")
  check(paste0(cm, ": q-values equal to BH over the tested proteins"), max_diff(num(tab[[paste0("qval_", cm)]]), q_ref) < 1e-10)
  check(paste0(cm, ": significant = q-value <= 0.05 and |log2 FC| >= 0.58"),
        identical(flag(tab[[paste0("significant_", cm)]]), expected_significant(tab, cm, q_ref)))
  check(paste0(cm, ": log2_fc = log2(fc)"), max_diff(num(tab[[paste0("log2_fc_", cm)]]), log2(num(tab[[paste0("fc_", cm)]]))) < 1e-10)
  s <- score(tab, cm)
  check(paste0(cm, ": at least 46 of the 48 spike-ins found, at most 10 yeast hits"), s["spike_ins"] >= 46 && s["yeast"] <= 10,
        paste0(s["spike_ins"], "/", s["ups"], " spike-ins, ", s["yeast"], " yeast"))
}
check("filtering summary describes the observation filter per group", grepl("observations in every group", report_text(r$html)))

## 2. Paired t-test: pairs by position (default) and by name
# In the experimental design the 500amol replicates are listed as inj1, inj3, inj2, inj4: the two pairings differ for 50fmol/500amol
cm <- "50fmol/500amol"
paired_t <- function(a, b) t.test(a, b, paired = TRUE)$p.value
r <- da_render(list(Test.paired = "TRUE"))
check("paired t-test renders (pairs by position)", is.na(r$error), r$error)
check("paired by position: pairs listed in the report", grepl("by position): inj1 - inj1, inj2 - inj3, inj3 - inj2, inj4 - inj4", report_text(r$html), fixed = TRUE))
if (!is.null(r$table)) {
  # test_values() returns the replicates in the design order, so the i-th columns of A and B form the i-th pair
  p_ref <- recompute_pvalues(r$table, cm, paired_t)
  check("paired by position: p-values equal to an independent paired t-test", max_diff(num(r$table[[paste0("pval_", cm)]]), p_ref) < 1e-10)
}
r <- da_render(list(Test.paired = "TRUE", Test.paired.by = "name"))
check("paired t-test renders (pairs by name)", is.na(r$error), r$error)
if (!is.null(r$table)) {
  # 500amol columns reordered to inj1, inj2, inj3, inj4, as the 50fmol columns
  p_ref <- recompute_pvalues(r$table, cm, paired_t, b_order = c(1, 3, 2, 4))
  p_position <- recompute_pvalues(r$table, cm, paired_t)
  check("paired by name: p-values equal to an independent paired t-test of the same-name replicates, and different from the pairing by position",
        max_diff(num(r$table[[paste0("pval_", cm)]]), p_ref) < 1e-10 && max_diff(num(r$table[[paste0("pval_", cm)]]), p_position) > 1e-3)
}

## 3. ANOVA with Holm post-hoc and BH gating on three comparisons
comparisons3 <- c("50fmol/25fmol", "50fmol/500amol", "25fmol/500amol")
r <- da_render(list(Test.type = "anova", Test.posthoc.procedure = "holm", Test.anova.gating = "bh", Comparisons = paste(comparisons3, collapse = ";")))
check("ANOVA with Holm post-hoc and BH gating renders", is.na(r$error), r$error)
tab <- r$table
if (!is.null(tab)) {
  check("ANOVA: anova_pvalue and anova_qvalue exported", all(c("anova_pvalue", "anova_qvalue") %in% names(tab)))
  ap <- num(tab$anova_pvalue)
  aq_ref <- ap; aq_ref[!is.na(ap)] <- p.adjust(ap[!is.na(ap)], "BH")
  check("ANOVA: anova_qvalue = BH of anova_pvalue", max_diff(num(tab$anova_qvalue), aq_ref) < 1e-10)
  for (cm in comparisons3) {
    check(paste0("ANOVA ", cm, ": raw and post-hoc p-values exported"), all(paste0(c("pval_", "posthoc_pval_"), cm) %in% names(tab)))
    sig <- flag(tab[[paste0("significant_", cm)]])
    check(paste0("ANOVA ", cm, ": significant proteins found, all passing the gate (anova_qvalue <= 0.05)"),
          sum(sig) > 0 && all(num(tab$anova_qvalue)[sig] <= 0.05), paste(sum(sig), "significant"))
    pp <- num(tab[[paste0("posthoc_pval_", cm)]]); praw <- num(tab[[paste0("pval_", cm)]])
    check(paste0("ANOVA ", cm, ": post-hoc (Holm) p-values not lower than the raw p-values"), all(pp >= praw - 1e-12, na.rm = TRUE))
  }
}

## 4. Kruskal-Wallis with Bonferroni post-hoc
r <- da_render(list(Test.type = "kruskall", Test.posthoc.procedure = "bonferroni", Comparisons = paste(comparisons3, collapse = ";")))
check("Kruskal-Wallis with Bonferroni post-hoc renders", is.na(r$error), r$error)
check("Kruskal-Wallis: kruskal_pvalue and kruskal_qvalue exported, no anova columns",
      !is.null(r$table) && all(c("kruskal_pvalue", "kruskal_qvalue") %in% names(r$table)) && !any(c("anova_pvalue", "anova_qvalue") %in% names(r$table)))

## 5. limma and the k-nearest-neighbour MCAR imputation
r <- da_render(list(Test.type = "limma"))
check("limma renders", is.na(r$error), r$error)
check("limma: p-values for both comparisons", !is.null(r$table) && all(sapply(comparisons, function(cm) sum(!is.na(num(r$table[[paste0("pval_", cm)]]))) > 2000)))
r <- da_render(list(Imputation.MCAR.model = "knn"))
check("knn MCAR imputation renders", is.na(r$error), r$error)
if (!is.null(r$table)) {
  s <- score(r$table, "50fmol/25fmol")
  check("knn: at least 46 spike-ins and at most 10 yeast hits on 50fmol/25fmol", s["spike_ins"] >= 46 && s["yeast"] <= 10,
        paste0(s["spike_ins"], " spike-ins, ", s["yeast"], " yeast"))
}

## 6. Invalid parameters and inputs stop with a clear message
r <- da_render(list(Volcano.color.style = "signifigance"))
check("Volcano.color.style = signifigance stops with a clear message", grepl("Volcano.color.style must be significance, intensity or monochrome", r$error), r$error)
r <- da_render(list(Test.type = "kruskall", Test.posthoc.procedure = "tukey", Comparisons = paste(comparisons3, collapse = ";")))
check("kruskall with tukey post-hoc stops with a clear message", grepl("tukey requires Test.type=anova", r$error, fixed = TRUE), r$error)
r <- da_render(list(Test.paired = "TRUE", Test.paired.by = "order"))
check("Test.paired.by = order stops with a clear message", grepl("Test.paired.by must be position or name", r$error), r$error)
dup <- read.delim(quant, check.names = FALSE, stringsAsFactors = FALSE)
dup$Accession[2] <- dup$Accession[1]
dup_file <- tempfile(fileext = ".txt"); write.table(dup, dup_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
r <- da_render(quant_df = dup_file)
check("duplicated Accession stops with a clear message", grepl("duplicated Accession", r$error), r$error)

finish("DA report tests")
