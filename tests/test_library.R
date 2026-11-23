## Tests of the functions of Library/quantitativeProteomics.R on the example dataset
## Usage: Rscript --vanilla tests/test_library.R (or tests/run_all.bat)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
suppressMessages(library(matrixStats))
source(file.path(repo_root, "Library", "quantitativeProteomics.R"))
cat("Library tests\n\n")

data <- read.csv(file.path(example_dir, "Parsed_proteins_set.txt"), check.names = FALSE, sep = "\t")
design <- read.csv(file.path(example_dir, "exp_design.txt"), sep = "\t")
files <- design$File; cond <- design$Condition
raw <- data[c("Id", files)]
msms <- as.matrix(data[paste0("Identification_type", sub("Intensity", "", files))]) == "By MS/MS"

## 1. Normalization: median of ratios to the reference sample
mat <- numeric_sample_matrix(raw[-1])
reference <- choose_reference_sample(mat)
check("reference sample = the sample with the fewest missing values", reference == files[which.min(colSums(is.na(mat)))])
norm <- normalize(raw, reference)
factors <- apply(sweep(mat, 1, mat[, reference], "/"), 2, median, na.rm = TRUE)
check("normalize(): each sample divided by the median of its ratios to the reference",
      isTRUE(all.equal(unname(as.matrix(norm[paste0(files, "_norm")])), unname(sweep(mat, 2, factors, "/")))))
check("normalize(): factors returned in the norm_factors attribute", isTRUE(all.equal(attr(norm, "norm_factors"), factors)))
subset_rows <- seq_len(nrow(raw)) <= 300
norm_sub <- normalize(raw, reference, rows = subset_rows)
factors_sub <- apply(sweep(mat[subset_rows, ], 1, mat[subset_rows, reference], "/"), 2, median, na.rm = TRUE)
check("normalize(rows = ): factors computed on the given proteins, applied to all",
      isTRUE(all.equal(unname(as.matrix(norm_sub[paste0(files, "_norm")])), unname(sweep(mat, 2, factors_sub, "/")))))
check("normalize(): ratio and normalized columns named as expected",
      all(c(paste0(setdiff(files, reference), "_ratio_", reference), paste0(files, "_norm")) %in% colnames(norm)))

## Kept proteins (global filter of the example parameters) and normalized log2 values
L <- as.matrix(norm[paste0(files, "_norm")]); colnames(L) <- files
X <- log2(L)
keep <- Reduce(`|`, lapply(unique(cond), function(g) rowMeans(!is.na(X[, cond == g, drop = FALSE])) >= 0.5)) & rowSums(msms, na.rm = TRUE) >= 2
df <- data.frame(Id = data$Id, L, check.names = FALSE)

## 2. Noise model
z <- sapply(2:6, function(k) { set.seed(k); draw_exact_sd(k) })
check("draw_exact_sd(): mean 0, sd 1 and values within +/- 2", all(sapply(z, function(v) abs(mean(v)) < 1e-12 && abs(sd(v) - 1) < 1e-12 && all(abs(v) <= 2))))
check("draw_exact_sd(1) = 0", draw_exact_sd(1) == 0)
model <- estimate_noise_model(df, cond, 0.05, rows = keep)
check("estimate_noise_model(): noise sd of the example (mean of the row sds, 0.563 log2 from 6,564 protein/condition rows)",
      abs(model$noise_sd - 0.563334) < 1e-5 && model$n_rows == 6564, sprintf("sd %.6f, rows %d", model$noise_sd, model$n_rows))
check("estimate_noise_model(): noise level = 5% percentile of each sample",
      isTRUE(all.equal(unname(model$noise_level), unname(apply(X[keep, ], 2, quantile, 0.05, na.rm = TRUE)))))
check("estimate_noise_model(): percentile outside ]0, 0.10] stops", grepl("at most 10", error_message(estimate_noise_model(df, cond, 0.2, rows = keep))))
check("estimate_noise_model(): too few rows stops with a clear message", grepl("are needed", error_message(estimate_noise_model(df[1:10, ], cond, 0.05))))

set.seed(123); before <- runif(1)
set.seed(123); a <- impute_using_noise_model(df, cond, model, rows = keep, seed = 1); after <- runif(1)
b <- impute_using_noise_model(df, cond, model, rows = keep, seed = 1)
A <- as.matrix(a[-1])
check("impute_using_noise_model(): reproducible with the same seed", identical(a, b))
check("impute_using_noise_model(): state of the random generator restored", identical(before, after))
check("impute_using_noise_model(): observed values unchanged", identical(unname(A[!is.na(L)]), unname(L[!is.na(L)])))
check("impute_using_noise_model(): kept proteins fully imputed, other proteins unchanged",
      !anyNA(A[keep, ]) && identical(unname(A[!keep, ]), unname(L[!keep, ])))
full_rows <- unlist(lapply(unique(cond), function(g) { s <- which(cond == g); r <- which(keep & rowSums(is.na(X[, s])) == length(s))
  sapply(r, function(i) sd(log2(A[i, s]) - model$noise_level[s])) }))
check("impute_using_noise_model(): fully missing protein/condition rows get exactly the noise sd", all(abs(full_rows - model$noise_sd) < 1e-9))
skip <- is.na(L) & keep & (seq_len(nrow(L)) %% 2 == 0)
c_skip <- as.matrix(impute_using_noise_model(df, cond, model, rows = keep, seed = 1, skip = skip)[-1])
check("impute_using_noise_model(skip = ): skipped cells left missing", all(is.na(c_skip[skip])))
untouched <- rowSums(skip) == 0
check("impute_using_noise_model(): values of a protein do not depend on the other proteins (one seed per protein and condition)",
      identical(unname(c_skip[untouched, ]), unname(A[untouched, ])))
p <- as.matrix(impute_using_noise_percentile(df, 0.05, skip = skip)[-1])
check("impute_using_noise_percentile(): missing values = 5% percentile of the sample, skipped cells left missing",
      all(is.na(p[skip])) && all(sapply(seq_along(files), function(j) { m <- is.na(L[, j]) & !skip[, j]; all(p[m, j] == quantile(L[, j], 0.05, na.rm = TRUE)) })))

## 3. MCAR classification and knn imputation
Xk <- X[keep, ]; msms_k <- msms[keep, ]
cutoff <- estimate_mcar_intensity_cutoff(Xk, cond)
check("estimate_mcar_intensity_cutoff(): cutoff of the example (about 21.5 log2)", isTRUE(abs(cutoff$cutoff - 21.49) < 0.05), paste("cutoff", cutoff$cutoff))
check("estimate_mcar_intensity_cutoff(): not estimated (NA, with a message) on too few rows", is.na(estimate_mcar_intensity_cutoff(Xk[1:50, ], cond)$cutoff))
mcar_all <- classify_mcar(Xk, cond, msms_k, 3, 3)
mcar_cut <- classify_mcar(Xk, cond, msms_k, 3, 3, cutoff = cutoff$cutoff)
check("classify_mcar(): 178 MCAR values without cutoff, 16 with the cutoff (example)", sum(mcar_all) == 178 && sum(mcar_cut) == 16,
      paste(sum(mcar_all), sum(mcar_cut)))
check("classify_mcar(): only missing values, and the cutoff selects a subset", all(is.na(Xk[mcar_all])) && all(mcar_all[mcar_cut]))
check("classify_mcar(): min_obs < 2 stops", grepl("at least 2", error_message(classify_mcar(Xk, cond, msms_k, 1, 0))))

set.seed(7); Xm <- Xk; mk <- matrix(FALSE, nrow(Xk), ncol(Xk))
for (g in unique(cond)) { s <- which(cond == g); r <- which(rowSums(is.na(Xk[, s])) == 0); r <- r[runif(length(r)) < 0.1]
  for (i in r) mk[i, s[sample.int(length(s), sample(1:2, 1))]] <- TRUE }
Xm[mk] <- NA
knn <- impute_mcar_neighbours(Xm, cond, mk, 3, seed = 1)
check("impute_mcar_neighbours(): reproducible, observed values unchanged, no MCAR value left missing",
      identical(knn, impute_mcar_neighbours(Xm, cond, mk, 3, seed = 1)) && identical(knn$values[!mk], Xm[!mk]) && !anyNA(knn$values[mk]))
model_rows <- condition_rows(Xm, cond); model_rows <- model_rows[model_rows$k >= 3, ]; model_rows <- model_rows[order(model_rows$mean), ]
nb <- knn$n_neighbours
props <- do.call(rbind, lapply(unique(cond), function(g) { s <- which(cond == g)
  do.call(rbind, lapply(which(rowSums(mk[, s, drop = FALSE]) > 0), function(i) {
    obs <- Xm[i, s][!is.na(Xm[i, s])]; k <- length(obs); nm <- sum(mk[i, s]); m <- mean(obs)
    pos <- findInterval(m, model_rows$mean); from <- min(max(pos - nb %/% 2, 0), nrow(model_rows) - nb)
    target <- median(model_rows$sd[(from + 1):(from + nb)]) * sqrt(length(s) / k); full <- knn$values[i, s]
    reachable <- if (nm == 1) (k * target^2 - (k - 1) * sd(obs)^2) > 0 else ((k + nm - 1) * target^2 - (k - 1) * sd(obs)^2) > 0
    data.frame(nm = nm, ok = if (reachable) abs(sd(full) - target) < 1e-9 else all(abs(full[mk[i, s]] - m) < 1e-9), mean_kept = abs(mean(full) - m) < 1e-9)
  })) }))
check("impute_mcar_neighbours(): sd of the protein = target (or mean imputed when the target cannot be reached)", all(props$ok))
check("impute_mcar_neighbours(): with several missing values, the mean of the protein is unchanged", all(props$mean_kept[props$nm > 1]))
check("impute_mcar_neighbours(): min_occurrences < 2 stops", grepl("at least 2", error_message(impute_mcar_neighbours(Xm, cond, mk, 1))))

## 4. CV and median deviation
toy <- data.frame(a = c(10, 10, NA, 5), b = c(20, NA, NA, 5), c = c(30, 12, 4, 5))
cv <- compute_cv(toy)
check("compute_cv(): sd / mean of the observed values, in %", isTRUE(all.equal(cv[1:2], c(sd(c(10, 20, 30)) / 20, sd(c(10, 12)) / 11) * 100)))
check("compute_cv(): NA with fewer than 2 observed values, 0 for constant values", is.na(cv[3]) && cv[4] == 0)
dev <- compute_median_deviation(toy)
check("compute_median_deviation(): (X - median) / median per protein", isTRUE(all.equal(dev[1, ], data.frame(a = -0.5, b = 0, c = 0.5), check.attributes = FALSE)))

## 5. with_fixed_seed()
set.seed(42); before <- runif(1); set.seed(42); x <- with_fixed_seed(5, runif(1)); after <- runif(1)
check("with_fixed_seed(): fixed result and state of the random generator restored",
      identical(x, { set.seed(5); runif(1) }) && identical(before, after))

## 6. validate_inputs()
check("validate_inputs(): example inputs are valid",
      is.na(error_message(validate_inputs(data, design, comparisons = "50fmol/25fmol", required_columns = c("Accession", "Gene_name", "Specific_peptides")))))
bad <- data; bad$Id[2] <- bad$Id[1]
check("validate_inputs(): duplicated Id", grepl("duplicated Id", error_message(validate_inputs(bad, design))))
check("validate_inputs(): missing required column", grepl("missing column\\(s\\) Specific_peptides",
      error_message(validate_inputs(data[setdiff(names(data), "Specific_peptides")], design, required_columns = "Specific_peptides"))))
bad <- data; bad$Accession[2] <- bad$Accession[1]
check("validate_inputs(): duplicated Accession (when Accession is required)", grepl("duplicated Accession", error_message(validate_inputs(bad, design, required_columns = "Accession"))))
bad <- data; bad$Specific_peptides[1:2] <- c(NA, -1)
check("validate_inputs(): invalid Specific_peptides values", grepl("2 value\\(s\\) of Specific_peptides", error_message(validate_inputs(bad, design, required_columns = "Specific_peptides"))))
check("validate_inputs(): unknown condition in a comparison", grepl("does not match any condition", error_message(validate_inputs(data, design, comparisons = "50fmol/2fmol"))))
check("validate_inputs(): intensity column of the design missing from the table",
      grepl("intensity column\\(s\\) of the experimental design not found", error_message(validate_inputs(data[setdiff(names(data), files[1])], design))))

## 7. plot_sample_legend(): legend panel shared by several plots
suppressMessages(library(ggplot2))
samples <- sub("Intensity_", "", files)  # 40 names of up to 11 characters
legend_panel <- plot_sample_legend(samples, rainbow(length(samples)))
legend_df <- legend_panel$plot$data
check("plot_sample_legend(): 40 samples in 8 columns of 5 rows, filled by columns, height for 5 rows",
      max(legend_df$column) == 7 && max(legend_df$row) == 4 && legend_df$row[6] == 0 && legend_df$column[6] == 1 && legend_panel$height == 60 + 26 * 5)
check("plot_sample_legend(): fewer samples than columns, one row",
      nrow(plot_sample_legend(c("A", "B", "C"), c("red", "green", "blue"))$plot$data) == 3 &&
        max(plot_sample_legend(c("A", "B", "C"), c("red", "green", "blue"))$plot$data$row) == 0)
check("plot_sample_legend(): rendered as a PNG figure",
      grepl("^<div class='figure-png'><img width='1200' src='data:image/png;base64,", plot_figure_as_html(legend_panel$plot, "Samples_legend", 1200, legend_panel$height)))

## 8. hclust_segments() and plot_clustered_heatmap()
hc <- hclust(dist(c(1, 2, 10, 11, 30)))
seg <- hclust_segments(hc)
check("hclust_segments(): 3 segments per merge, leaves at positions 1..n and height 0, top merge at the largest height",
      nrow(seg) == 3 * 4 && setequal(seg$x[seg$y == 0], 1:5) && max(seg$yend) == max(hc$height))
values <- as.matrix(data[files[1:8]]); values <- log10(values[complete.cases(values), ][1:300, ])
check("plot_clustered_heatmap(): rendered as a PNG figure (300 proteins, 8 samples, z-scores), in both orientations",
      grepl("^<div class='figure-png'>", plot_figure_as_html(plot_clustered_heatmap(t(scale(t(values))), limits = c(-3, 3)), "heatmap", 800, 600)) &&
        grepl("^<div class='figure-png'>", plot_figure_as_html(plot_clustered_heatmap(values, proteins_in_rows = FALSE), "heatmap", 800, 600)))
check("heatmap_proteins_in_rows(): proteins in rows with thousands of proteins and 40 samples, not with 500 proteins and 300 samples",
      heatmap_proteins_in_rows(2898, 40) && !heatmap_proteins_in_rows(500, 300))

## 10. condition_colors(): Paul Tol "bright" up to 7 conditions, "Paired" beyond
check("condition_colors(): distinct colors, Tol bright up to 7 (blue, red, green, yellow for 4), Paired up to 12, interpolated beyond",
      identical(condition_colors(4), c("#4477AA", "#EE6677", "#228833", "#CCBB44")) && length(unique(condition_colors(7))) == 7 &&
        identical(condition_colors(8), RColorBrewer::brewer.pal(8, "Paired")) && identical(condition_colors(10), RColorBrewer::brewer.pal(10, "Paired")) &&
        length(unique(condition_colors(20))) == 20)

## 11. Text measures: widths measured on the PNG device, sizes fitted to them
w <- text_width_px(c("C1", "C2_SAX HP like_1"), 14)
check("text_width_px(): longer texts wider, width proportional to the size, plotmath labels measured, device state restored",
      w[2] > 3 * w[1] && abs(text_width_px("C2_SAX HP like_1", 28) / w[2] - 2) < 0.1 &&
        text_width_px(arrow_label("up", "higher in A"), 14, plotmath = TRUE) > text_width_px("higher in A", 14) && grDevices::dev.cur() == 1)
check("sample_axis_text_size(): at most 14 pt, smaller for many samples and for long names (at most 150 px long)",
      sample_axis_text_size(10) == 14 && sample_axis_text_size(80) < 6 + 1e-9 + 520 / 80 / 1.25 &&
        text_width_px(paste(rep("long name", 3), collapse = " "), sample_axis_text_size(10, labels = paste(rep("long name", 3), collapse = " "))) <= 151)
check("condition_legend_guide(): a ggplot guides object (legend wrapping in the figure width)",
      inherits(condition_legend_guide(c("C1_Kit HP", "C2_SAX HP like", "C3_SAX", "C4_NEAT"), 600), "Guides") ||
        is.list(condition_legend_guide(c("C1_Kit HP", "C2_SAX HP like"), 600)))

## 9. arrow_label(): plotmath labels with an arrow, whatever the characters of the text
check("arrow_label(): valid plotmath expressions, also with quotes and backslashes in the text",
      all(sapply(c(arrow_label("up", "higher in A"), arrow_label("right", 'higher in "B" \\ C', arrow_after = TRUE)),
                 function(l) !inherits(try(parse(text = l), silent = TRUE), "try-error"))) &&
        grepl("%up%", arrow_label("up", "A")) && grepl("%->%", arrow_label("right", "A", arrow_after = TRUE)))

## 10b. display_sample_names(): names of the samples in the figures (QC and DA)
check("display_sample_names(): names kept up to 15 characters, else <Condition>_<index> (numbered within the condition), or S<index>",
      identical(display_sample_names(c("inj1_A", "inj2_A")), c("inj1_A", "inj2_A")) &&
        identical(display_sample_names(c("OAEMB261105_PA_04_C1", "OAEMB261105_PA_08_C2", "OAEMB261105_PA_10_C1"), c("C1", "C2", "C1")), c("C1_1", "C2_1", "C1_2")) &&
        identical(display_sample_names(c("a_very_long_sample_name", "b")), c("S1", "S2")))

## 10. exclude_samples(): samples removed from the design and the table (all their columns), checks after the filtering
ex <- exclude_samples(data, design, c("10amol_1", "Intensity_50amol_3"), comparisons = "50amol/10amol")
check("exclude_samples(): design rows and all the columns of the samples removed, the others kept",
      nrow(ex$expDesign) == nrow(design) - 2 && !any(c("Intensity_10amol_1", "Identification_type_10amol_1", "Intensity_50amol_3", "Identification_type_50amol_3") %in% names(ex$data)) &&
        ncol(ex$data) == ncol(data) - 4 && "Intensity_10amol_2" %in% names(ex$data) &&
        is.na(error_message(validate_inputs(ex$data, ex$expDesign, comparisons = "50amol/10amol"))))
check("exclude_samples(): nothing to exclude leaves the inputs unchanged", identical(exclude_samples(data, design, character(0)), list(data = data, expDesign = design)))
check("exclude_samples(): unknown sample stops", grepl("not found in the experimental design: 10amol_9", error_message(exclude_samples(data, design, "10amol_9"))))
check("exclude_samples(): compared condition left with fewer than 2 replicates stops",
      grepl("fewer than 2 replicates: 10amol \\(1\\)", error_message(exclude_samples(data, design, paste0("10amol_", 1:3), comparisons = "50amol/10amol"))))

finish("Library tests")
