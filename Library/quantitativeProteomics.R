
#####################
# Inputs validation #
#####################

# This function checks that the experimental design and the quantification table are consistent and unambiguous,
# so that per-sample columns can be selected by their exact names. It stops with the list of all the problems found.
validate_inputs <- function(data, expDesign, comparisons = NULL, derived_suffixes = c("_mean", "_identification_type"), required_columns = character(0)){

  # Inputs :
  # data : quantification table (with an "Id" column, "Intensity_*" and "Identification_type_*" columns)
  # expDesign : experimental design (columns File, repBio, Condition)
  # comparisons : optional vector of comparisons "A/B"
  # derived_suffixes : suffixes appended to sample ids ("<repBio>_<Condition>") to create new columns
  # required_columns : other columns the quantification table must contain ("Specific_peptides" must also hold counts: numbers >= 0, no missing value)
  # Output :
  # invisible TRUE, or an error listing the problems

  problems <- c()

  # Experimental design
  required <- c("File", "repBio", "Condition")
  missing_design_cols <- setdiff(required, colnames(expDesign))
  if (length(missing_design_cols) > 0) {
    stop("Experimental design: missing column(s) ", paste(missing_design_cols, collapse = ", "), " (expected: File, repBio, Condition).")
  }
  for (col in required) {
    n_empty <- sum(is.na(expDesign[[col]]) | trimws(as.character(expDesign[[col]])) == "")
    if (n_empty > 0) problems <- c(problems, paste0("Experimental design: ", n_empty, " empty value(s) in column ", col))
  }
  files <- as.character(expDesign$File)
  if (anyDuplicated(files)) {
    problems <- c(problems, paste0("Experimental design: duplicated File name(s): ", paste(unique(files[duplicated(files)]), collapse = ", ")))
  }
  sample_ids <- paste0(expDesign$repBio, "_", expDesign$Condition)
  pairs_by_sample <- tapply(paste(expDesign$repBio, expDesign$Condition, sep = "\r"), sample_ids, function(x) length(unique(x)))
  if (any(pairs_by_sample > 1)) {
    problems <- c(problems, paste0("Experimental design: ambiguous sample id(s) <repBio>_<Condition> obtained from different repBio/Condition pairs: ",
                                   paste(names(pairs_by_sample)[pairs_by_sample > 1], collapse = ", ")))
  }

  # Quantification table
  if (anyDuplicated(colnames(data))) {
    problems <- c(problems, paste0("Quantification table: duplicated column name(s): ", paste(unique(colnames(data)[duplicated(colnames(data))]), collapse = ", ")))
  }
  if (!"Id" %in% colnames(data)) {
    problems <- c(problems, "Quantification table: missing 'Id' column")
  } else if (anyDuplicated(data$Id)) {
    problems <- c(problems, paste0("Quantification table: duplicated Id value(s): ", paste(head(unique(data$Id[duplicated(data$Id)]), 10), collapse = ", ")))
  }
  missing_files <- setdiff(files, colnames(data))
  if (length(missing_files) > 0) {
    problems <- c(problems, paste0("Quantification table: intensity column(s) of the experimental design not found: ", paste(missing_files, collapse = ", ")))
  }
  id_type_cols <- paste0("Identification_type", sub("Intensity", "", files))
  missing_id_types <- setdiff(id_type_cols, colnames(data))
  if (length(missing_id_types) > 0) {
    problems <- c(problems, paste0("Quantification table: identification type column(s) not found: ", paste(missing_id_types, collapse = ", ")))
  }
  missing_required <- setdiff(required_columns, colnames(data))
  if (length(missing_required) > 0) {
    problems <- c(problems, paste0("Quantification table: missing column(s) ", paste(missing_required, collapse = ", "),
                                   " (expected in the parsed table, see Parsers)"))
  }
  if ("Specific_peptides" %in% intersect(required_columns, colnames(data))) {
    counts <- suppressWarnings(as.numeric(as.character(data$Specific_peptides)))
    n_bad <- sum(is.na(counts) | counts < 0)
    if (n_bad > 0) problems <- c(problems, paste0("Quantification table: ", n_bad, " value(s) of Specific_peptides are missing, not numbers or negative"))
  }
  derived_cols <- c(outer(unique(sample_ids), derived_suffixes, paste0))
  colliding_cols <- intersect(derived_cols, colnames(data))
  if (length(colliding_cols) > 0) {
    problems <- c(problems, paste0("Quantification table: column(s) already named like the columns created by the analysis: ", paste(colliding_cols, collapse = ", ")))
  }

  # Comparisons
  if (length(comparisons) > 0) {
    conditions <- unique(as.character(expDesign$Condition))
    for (comparison in comparisons) {
      couple <- strsplit(comparison, "/")[[1]]
      if (length(couple) != 2) {
        problems <- c(problems, paste0("Invalid comparison \"", comparison, "\": expected two conditions separated by \"/\" (e.g. \"",
                                       paste(head(conditions, 2), collapse = "/"), "\")"))
        next
      }
      unknown <- setdiff(couple, conditions)
      if (length(unknown) > 0) {
        problems <- c(problems, paste0("Invalid comparison \"", comparison, "\": ", paste0("\"", unknown, "\"", collapse = " and "),
                                       if (length(unknown) > 1) " do" else " does",
                                       " not match any condition defined in the experimental design (existing conditions: ",
                                       paste(conditions, collapse = ", "), ")"))
      } else if (couple[1] == couple[2]) {
        problems <- c(problems, paste0("Invalid comparison \"", comparison, "\": a condition cannot be compared with itself"))
      }
    }
  }

  if (length(problems) > 0) {
    stop(paste(c("Invalid inputs:", problems), collapse = "\n - "), call. = FALSE)
  }

  # Intensity columns absent from the design are allowed (samples can be excluded on purpose), but reported
  unused_files <- setdiff(grep("^Intensity", colnames(data), value = TRUE), files)
  if (length(unused_files) > 0) {
    warning("Intensity column(s) not used because absent from the experimental design: ", paste(unused_files, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}

#################
# Normalization #
#################

# Normalization by the median of ratios to a reference sample: each sample is divided by the median of its ratios
# to the reference sample (equivalently, the median log2 ratio is subtracted on the log2 scale).
# normalize() is the entry point; the helpers below are the steps, kept separate for testing and QC.

# Converts the sample columns of a table to a numeric matrix (same column names)
numeric_sample_matrix <- function(df){
  mat <- do.call(cbind, lapply(df, function(x) if (is.numeric(x)) as.numeric(x) else as.numeric(as.character(unlist(x)))))
  colnames(mat) <- colnames(df)
  mat
}

# Returns the name of the reference sample: the (first) sample with the fewest missing values
choose_reference_sample <- function(mat){
  mat <- as.matrix(mat)
  colnames(mat)[which.min(colSums(is.na(mat)))]
}

# Returns the matrix of the ratios sample / reference, protein by protein (the reference column is 1 where it is observed)
compute_ratios <- function(mat, reference){
  if (!reference %in% colnames(mat)) stop("Normalization: reference sample '", reference, "' not found in the intensity columns.")
  mat / mat[, reference]
}

# Returns the normalization factor of each sample: the median of its ratios to the reference (1 for the reference).
# rows : optional subset of the proteins used to compute the factors (logical vector or row indices)
compute_norm_factors <- function(ratios, rows = NULL){
  if (!is.null(rows)) ratios <- ratios[rows, , drop = FALSE]
  factors <- apply(ratios, 2, median, na.rm = TRUE)
  bad <- names(factors)[is.na(factors) | !is.finite(factors) | factors <= 0]
  if (length(bad) > 0) {
    stop("Normalization: no usable ratio to the reference sample for ", paste(bad, collapse = ", "),
         " (no protein quantified in both samples among the proteins used for normalization).")
  }
  factors
}

# Divides each sample by its normalization factor
apply_norm_factors <- function(mat, factors){
  sweep(mat, 2, factors[colnames(mat)], "/")
}

normalize <- function(df, reference_sample, rows = NULL){

  # Inputs :
  # df : quantification matrix whom 1rst column contains proteins id
  # reference_sample : name of the column corresponding to the reference sample
  # rows : optional subset of the proteins used to compute the normalization factors (logical vector or row indices);
  #        the factors are applied to all the proteins
  # Output :
  # dataframe containing the ratio values for each sample ("<sample>_ratio_<reference>"), and the normalized
  # quantification matrix ("<sample>_norm"). The normalization factors are returned in the "norm_factors" attribute.

  mat <- numeric_sample_matrix(df[-1])
  ratios <- compute_ratios(mat, reference_sample)
  factors <- compute_norm_factors(ratios, rows)
  normalized <- apply_norm_factors(mat, factors)

  others <- colnames(mat) != reference_sample
  normalized_data <- data.frame(df[[1]], ratios[, others, drop = FALSE], normalized, check.names = FALSE)
  colnames(normalized_data) <- c("Id", paste0(colnames(mat)[others], "_ratio_", reference_sample), paste0(colnames(mat), "_norm"))
  attr(normalized_data, "norm_factors") <- factors

  return(normalized_data)

}



######
# CV #
######

# This function compute the median deviation for each sample of the quantification matrix.
compute_median_deviation <- function(intensity_matrix){
  
  # Input :
  # intensity_matrix : quantification matrix
  # Output :
  # matrix containing the median deviation for each intensity of each sample
  
  intensity_matrix$medians = rowMedians(data.matrix(intensity_matrix), na.rm=T)
  for(sample in colnames(intensity_matrix)[-ncol(intensity_matrix)]){
    intensity_matrix[sample]=(as.numeric(intensity_matrix[[sample]])-intensity_matrix$medians)/intensity_matrix$medians
  }
  
  return(intensity_matrix[,-which(names(intensity_matrix)=="medians")])
  
}

library(matrixStats)

# This function compute the median deviation for each sample of the quantification matrix. It requires the library "matrixStats".
compute_cv <- function(intensity_matrix){
  
  # Input :
  # intensity_matrix : quantification matrix
  # Output :
  # vector containing the CV for each protein
  
  medians = rowMedians(data.matrix(intensity_matrix), na.rm=T)
  sd = apply(intensity_matrix,1,sd)
  cv = sd/medians * 100
  
  return(cv)
  
}

##############
# Imputation #
##############
# Converts a table whose cells may contain the "MCAR" marker (cells deliberately left missing) to a numeric matrix.
# data.matrix() must not be used here: it turns character columns into factor codes.
numeric_intensity_matrix <- function(df){
  df[df=="MCAR"] <- NA
  as.matrix(as.data.frame(lapply(df, function(x) if (is.numeric(x)) as.numeric(x) else as.numeric(as.character(x))), check.names = FALSE))
}

# This function impute missing values with background noise using a percentile.
impute_using_noise_percentile <- function(df,per,skip = NULL){

  # Input :
  # df : quantification matrix whom 1rst column contains protein id
  # per : chosen percentile
  # skip : optional logical matrix (one column per sample column of df): missing values to leave missing (MCAR values not imputed)
  # Output :
  # quantification matrix with imputed intensities

  for(j in seq_along(colnames(df))[-1]){
    current_intensities = df[[j]]
    df_without_string = numeric_intensity_matrix(df[j])
    percentile = as.numeric(quantile(df_without_string,per, na.rm = T))
    to_impute = is.na(current_intensities)
    if (!is.null(skip)) to_impute = to_impute & !skip[, j - 1]
    current_intensities[to_impute] <- percentile
    df[[j]] = current_intensities
  }
  
  colnames(df)[2:ncol(df)] = paste0(colnames(df)[2:ncol(df)],"_imputed")
  return(df)
  
}

###################
# MCAR imputation #
###################

# Seed of the random draws of one protein (row) in one condition: the values imputed for a protein and a condition depend only
# on the seed and on this protein, not on the other proteins (so that runs with different settings can be compared protein by protein)
row_seed <- function(seed, row, group_index, n_groups){
  seed + (row - 1) * n_groups + group_index
}

# Runs code with a fixed seed of the random generator, and restores the state of the generator afterwards
with_fixed_seed <- function(seed, code){
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old_seed)) suppressWarnings(rm(".Random.seed", envir = globalenv())) else assign(".Random.seed", old_seed, envir = globalenv()))
  set.seed(seed)
  code
}

# Returns one row per protein and condition with at least one observed value: protein (row index), condition, mean, sd
# (NA with less than 2 values), number of observed values (k) and number of samples of the condition (n), on the log2 matrix X
condition_rows <- function(X, groups){
  do.call(rbind, lapply(unique(groups), function(g){
    Y <- X[, groups == g, drop = FALSE]
    k <- rowSums(!is.na(Y)); ok <- which(k >= 1)
    m <- rowMeans(Y[ok, , drop = FALSE], na.rm = TRUE)
    s <- sqrt(rowSums((Y[ok, , drop = FALSE] - m)^2, na.rm = TRUE) / (k[ok] - 1)); s[k[ok] < 2] <- NA
    data.frame(protein = ok, condition = g, mean = m, sd = s, k = k[ok], n = ncol(Y), row.names = NULL)
  }))
}

# Estimates the intensity above which missing values are considered as missing completely at random (MCAR).
# The missing-value rate of the protein-condition rows is computed in bins of their mean intensity (equal-count bins), and smoothed
# with a decreasing (isotonic) fit. The plateau is the median fitted rate of the upper half of the bins. The cutoff is the upper
# edge of the highest bin whose fitted rate exceeds the plateau by more than 'tol' (-Inf if none: no intensity dependence).
estimate_mcar_intensity_cutoff <- function(X, groups, n_bins = 25, tol = 0.01, min_rows = 500){

  # Inputs :
  # X : log2 intensities (proteins kept for the analysis), one column per sample
  # groups : condition of each column of X
  # Output :
  # list(cutoff (log2, NA if it cannot be estimated), plateau, bins (data frame), message)

  r <- condition_rows(X, as.character(groups))
  if (nrow(r) < min_rows) {
    return(list(cutoff = NA, plateau = NA, bins = NULL,
                message = paste0("too few protein/condition rows (", nrow(r), ", at least ", min_rows, " are needed) to estimate the intensity cutoff")))
  }
  breaks <- unique(quantile(r$mean, seq(0, 1, length.out = n_bins + 1), names = FALSE))
  r$bin <- cut(r$mean, breaks, include.lowest = TRUE, labels = FALSE)
  bins <- do.call(rbind, lapply(sort(unique(r$bin)), function(b){
    d <- r[r$bin == b, ]
    data.frame(bin = b, from = min(d$mean), to = max(d$mean), center = median(d$mean), rows = nrow(d), rate = sum(d$n - d$k) / sum(d$n))
  }))
  bins$fitted <- rev(isoreg(rev(bins$rate))$yf)
  upper <- bins$fitted[(nrow(bins) %/% 2 + 1):nrow(bins)]
  plateau <- median(upper)
  above <- which(bins$fitted > plateau + tol)
  cutoff <- if (length(above) > 0) bins$to[max(above)] else -Inf
  list(cutoff = cutoff, plateau = plateau, bins = bins, message = NULL)
}

# Classifies the missing values as MCAR: the missing values of a protein in a condition are MCAR if the protein has at least
# min_obs observed values and at least min_msms MS/MS identifications in this condition, and (if a cutoff is given) a mean of its
# observed values at least equal to the cutoff. The other missing values are MNAR.
classify_mcar <- function(X, groups, msms, min_obs, min_msms, cutoff = NA){

  # Inputs :
  # X : log2 intensities, one column per sample
  # groups : condition of each column of X
  # msms : logical matrix (same dimensions as X), TRUE when the protein is identified by MS/MS in the sample
  # min_obs, min_msms : thresholds in each condition (min_obs >= 2)
  # cutoff : log2 intensity cutoff (NA: not used)
  # Output :
  # logical matrix (same dimensions as X), TRUE for the MCAR missing values

  if (is.na(min_obs) || min_obs < 2) stop("Imputation.MCAR.threshold.obs must be at least 2.", call. = FALSE)
  groups <- as.character(groups)
  msms[is.na(msms)] <- FALSE
  mcar <- matrix(FALSE, nrow(X), ncol(X), dimnames = dimnames(X))
  for (g in unique(groups)) {
    s <- which(groups == g)
    k <- rowSums(!is.na(X[, s, drop = FALSE]))
    m <- suppressWarnings(rowMeans(X[, s, drop = FALSE], na.rm = TRUE))
    rows <- k >= min_obs & k < length(s) & rowSums(msms[, s, drop = FALSE]) >= min_msms
    if (!is.na(cutoff)) rows <- rows & !is.na(m) & m >= cutoff
    mcar[rows, s] <- is.na(X[rows, s, drop = FALSE])
  }
  mcar
}

# Imputes the MCAR missing values from the neighbours ("knn" model), on the log2 scale.
# - Neighbours of a protein in a condition: the sqrt(number of proteins) protein-condition rows nearest in mean intensity, among the rows
#   with at least min_occurrences observed values.
# - Target sd: the median sd of the neighbours, multiplied by sqrt(n/k) (n samples in the condition, k observed values), because the
#   test then counts the imputed values as measured.
# - One missing value: the value which brings the sd of the protein to the target (on a random side of its mean); its mean if the target
#   cannot be reached. Several missing values: drawn jointly so that the mean of the protein is unchanged and its sd equals the target.
impute_mcar_neighbours <- function(X, groups, mcar, min_occurrences, seed = 1, cutoff = 2){

  # Inputs :
  # X : log2 intensities (proteins kept for the analysis), one column per sample
  # groups : condition of each column of X
  # mcar : logical matrix of the MCAR missing values (see classify_mcar)
  # min_occurrences : minimum number of observed values of a neighbour row (at least 2)
  # Output :
  # list(values (X with the MCAR values imputed), n_neighbours, n_model_rows)

  if (is.na(min_occurrences) || min_occurrences < 2) stop("Imputation.knn.min.occurrences must be at least 2.", call. = FALSE)
  groups <- as.character(groups)
  model <- condition_rows(X, groups)
  model <- model[model$k >= min_occurrences, ]
  model <- model[order(model$mean), ]
  n_nb <- floor(sqrt(nrow(X)))
  if (nrow(model) < n_nb) {
    stop("MCAR imputation (knn): only ", nrow(model), " protein/condition rows with at least ", min_occurrences,
         " observed values are available as neighbours (at least ", n_nb, " are needed).", call. = FALSE)
  }
  neighbour_sd <- function(v){
    p <- findInterval(v, model$mean)
    from <- min(max(p - n_nb %/% 2, 0), nrow(model) - n_nb)
    median(model$sd[(from + 1):(from + n_nb)])
  }

  out <- X
  with_fixed_seed(seed, {
    for (g in unique(groups)) {
      s <- which(groups == g)
      for (i in which(rowSums(mcar[, s, drop = FALSE]) > 0)) {
        set.seed(row_seed(seed, i, match(g, unique(groups)), length(unique(groups))))
        miss <- mcar[i, s]
        obs <- X[i, s][!is.na(X[i, s])]
        k <- length(obs); nm <- sum(miss); m <- mean(obs); s_own <- sd(obs)
        target <- neighbour_sd(m) * sqrt(length(s) / k)
        if (nm == 1) {
          d2 <- (k * target^2 - (k - 1) * s_own^2) * (k + 1) / k
          values <- m + sample(c(-1, 1), 1) * sqrt(max(d2, 0))
        } else {
          need <- (k + nm - 1) * target^2 - (k - 1) * s_own^2
          values <- m + draw_exact_sd(nm, cutoff) * sqrt(max(need, 0) / (nm - 1))
        }
        out[i, s[miss]] <- values
      }
    }
  })
  list(values = out, n_neighbours = n_nb, n_model_rows = nrow(model))
}

###############
# Noise model #
###############

# Draws k standard values with a mean of exactly 0 and a sd of exactly 1, all within +/- cutoff (0 when k = 1).
# Single draws beyond the cutoff are redrawn, and the whole set is redrawn if a value is beyond the cutoff after standardization.
draw_exact_sd <- function(k, cutoff = 2){
  if (k == 1) return(0)
  repeat {
    z <- rnorm(k)
    while (any(far <- abs(z) > cutoff)) z[far] <- rnorm(sum(far))
    z <- z - mean(z)
    z <- z / sd(z)
    if (all(abs(z) <= cutoff)) return(z)
  }
}

# Returns the log2 intensities of a quantification matrix (linear intensities, 1rst column = protein id).
# Non-positive intensities and cells holding the "MCAR" marker are missing values.
log2_intensity_matrix <- function(df){
  linear <- numeric_intensity_matrix(df[-1]); colnames(linear) <- colnames(df)[-1]
  linear[!is.na(linear) & linear <= 0] <- NA
  log2(linear)
}

# This function estimates the background noise model used to impute missing values (MNAR), on the log2 scale.
# - Noise level of a sample: the given percentile of its observed values.
# - Noise sd: pooled within-condition sd of the protein-condition rows with at least 2 observed values whose mean lies
#   between the given percentile and the upper percentile (25%) of all the observed values.
estimate_noise_model <- function(df, groups, percentile, rows = NULL, upper_percentile = 0.25, min_rows = 100){

  # Inputs :
  # df : quantification matrix (linear intensities) whom 1rst column contains protein id
  # groups : condition of each sample column (same order as the columns of df, without the id column)
  # percentile : percentile of the noise level, as a fraction (0 < percentile <= 0.10)
  # rows : proteins used to estimate the model (logical vector, default: all)
  # Output :
  # list(noise_level (log2, per sample), noise_sd (log2), n_rows (protein-condition rows used for the sd), range (log2 bounds of these rows))

  if (is.na(percentile) || percentile <= 0 || percentile > 0.10) {
    stop("Imputation.MNAR.percentile must be greater than 0 and at most 10 (%) for the noise model.", call. = FALSE)
  }
  groups <- as.character(groups)
  if (length(groups) != ncol(df) - 1) stop("estimate_noise_model: 'groups' must give the condition of each sample column.")
  if (is.null(rows)) rows <- rep(TRUE, nrow(df))
  kept <- log2_intensity_matrix(df)[rows, , drop = FALSE]

  # Noise level of each sample
  noise_level <- apply(kept, 2, function(x) if (all(is.na(x))) NA else quantile(x, percentile, na.rm = TRUE, names = FALSE))
  if (any(is.na(noise_level))) {
    stop("Imputation: no observed value in sample(s) ", paste(colnames(kept)[is.na(noise_level)], collapse = ", "), " among the proteins kept for the analysis.", call. = FALSE)
  }

  # Noise sd: pooled within-condition sd of the low-intensity protein-condition rows
  stats <- do.call(rbind, lapply(unique(groups), function(g){
    Y <- kept[, groups == g, drop = FALSE]
    n <- rowSums(!is.na(Y)); ok <- n >= 2
    m <- rowMeans(Y[ok, , drop = FALSE], na.rm = TRUE)
    data.frame(mean = m, ss = rowSums((Y[ok, , drop = FALSE] - m)^2, na.rm = TRUE), df = n[ok] - 1)
  }))
  range <- quantile(kept[!is.na(kept)], c(percentile, upper_percentile), names = FALSE)
  in_range <- stats$mean >= range[1] & stats$mean <= range[2]
  if (sum(in_range) < min_rows) {
    stop("Imputation: only ", sum(in_range), " low-intensity protein/condition rows with at least 2 observed values are available to estimate the noise sd (at least ",
         min_rows, " are needed). Use the percentile model, or disable the imputation.", call. = FALSE)
  }

  list(noise_level = noise_level, noise_sd = sqrt(sum(stats$ss[in_range]) / sum(stats$df[in_range])), n_rows = sum(in_range), range = range)
}

# This function imputes missing values with background noise, using a model returned by estimate_noise_model().
# For each protein and condition with at least one missing value, a set of values is drawn for all the samples of the condition
# around their noise levels, with a sd of exactly the noise sd (see draw_exact_sd), and the observed values are kept.
impute_using_noise_model <- function(df, groups, noise_model, rows = NULL, seed = 1, cutoff = 2, skip = NULL){

  # Inputs :
  # df : quantification matrix (linear intensities) whom 1rst column contains protein id. Cells holding the "MCAR" marker are left missing.
  # groups : condition of each sample column (same order as the columns of df, without the id column)
  # noise_model : model returned by estimate_noise_model() for the same matrix
  # rows : proteins to impute (logical vector, default: all); the other proteins are returned unchanged
  # seed : seed of the random generator (the imputation is reproducible; the state of the generator is restored afterwards)
  # skip : optional logical matrix (one column per sample column of df): missing values to leave missing (MCAR values not imputed)
  # Output :
  # quantification matrix with imputed intensities ("_imputed" columns)

  groups <- as.character(groups)
  if (length(groups) != ncol(df) - 1) stop("impute_using_noise_model: 'groups' must give the condition of each sample column.")
  if (is.null(rows)) rows <- rep(TRUE, nrow(df))

  left_missing <- as.matrix(df[-1] == "MCAR"); left_missing[is.na(left_missing)] <- FALSE
  if (!is.null(skip)) left_missing <- left_missing | skip
  X <- log2_intensity_matrix(df)
  noise_level <- noise_model$noise_level[colnames(X)]

  # Imputation, with a fixed seed per protein and condition (the state of the random generator is restored afterwards)
  imputed <- X
  to_impute <- is.na(X) & !left_missing
  with_fixed_seed(seed, {
    for (g in unique(groups)) {
      s <- which(groups == g)
      for (i in which(rows & rowSums(to_impute[, s, drop = FALSE]) > 0)) {
        set.seed(row_seed(seed, i, match(g, unique(groups)), length(unique(groups))))
        full <- noise_level[s] + noise_model$noise_sd * draw_exact_sd(length(s), cutoff)
        miss <- to_impute[i, s]
        imputed[i, s[miss]] <- full[miss]
      }
    }
  })

  # Back to linear intensities; the observed values are returned exactly as given
  linear <- numeric_intensity_matrix(df[-1])
  out <- 2^imputed
  observed <- !is.na(linear) & linear > 0
  out[observed] <- linear[observed]
  table <- data.frame(df[[1]], out, check.names = FALSE)
  colnames(table) <- c(colnames(df)[1], paste0(colnames(df)[-1], "_imputed"))

  return(table)
}

###########
# Heatmap #
###########

# This function create a heatmap with a hierarchical clustering of proteins
create_heatmap <- function(intensities,imputed,samples_names,group_ids,htitle){
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(ggdendro)
  library(patchwork)
  
  # Input:
  # intensities : quantification matrix where the first two columns are identifiers
  #               Example: Accession, Label, sample1, sample2, ...
  #
  # imputed : boolean matrix used as a missing-value mask
  #           Expected structure:
  #           Id, Label, sample1_mean, sample2_mean, ...
  #
  # samples_names : list of sample names without "_mean"
  #                 Example: K13D_1_C3_K13D, K13D_2_C3_K13D, ...
  #
  # group_ids : list of biological conditions
  #             Example: C3_K13D, C4_K13ND
  #
  # Output:
  # heatmap plot
  
  
  ############################
  # 1. Basic safety checks
  ############################
  
  # The heatmap matrix needs at least:
  # column 1 = Accession or protein identifier
  # column 2 = Label or gene name
  # column 3+ = intensity columns
  if (ncol(intensities) < 3) {
    stop("create_heatmap(): intensities must contain at least 3 columns: identifier, label, and intensity values.")
  }
  
  # The imputed mask should have the same number of protein rows as the intensity matrix.
  if (nrow(intensities) != nrow(imputed)) {
    stop("create_heatmap(): intensities and imputed do not have the same number of rows.")
  }
  
  # In your current setup:
  # intensities[[1]] is Accession
  # imputed[[1]] is named Id, but contains the same Accession values.
  #
  # This check makes sure the protein order is the same before any long-format reshaping.
  if (!identical(as.character(intensities[[1]]), as.character(imputed[[1]]))) {
    stop("create_heatmap(): protein order differs between intensities and imputed mask.")
  }
  
  
  ############################
  # 2. Build numeric matrix for clustering
  ############################
  
  # This replaces:
  # inputProteins = intensities[,3:ncol(intensities)]
  #
  # We keep only the sample intensity columns.
  inputProteins <- intensities[, 3:ncol(intensities), drop = FALSE]
  
  # Make sure all intensity columns are numeric.
  # This avoids hidden character columns causing strange behavior in dist/hclust/ggplot.
  inputProteins <- data.frame(
    lapply(inputProteins, as.numeric),
    check.names = FALSE
  )
  
  # Use gene labels as rownames for clustering.
  # make.names(..., unique = TRUE) prevents duplicate rowname problems.
  rownames(inputProteins) <- make.names(intensities[[2]], unique = TRUE)
  
  # Cluster proteins based on the intensity matrix.
  distMatrixProteins <- dist(inputProteins)
  hClustProteins <- hclust(distMatrixProteins)
  
  
  ############################
  # 3. Convert intensity matrix to long format
  ############################
  
  # These are the heatmap sample columns.
  # Example:
  # K13D_1_C3_K13D, K13D_2_C3_K13D, ...
  intensity_cols <- colnames(intensities)[3:ncol(intensities)]
  
  # This replaces the old manual construction:
  #
  # proteinRank = c()
  # for (i in 1:(ncol(inputProteins))){
  #   proteinRank = c(proteinRank,1:nrow(inputProteins))
  # }
  # proteinId = rep(intensities[,1],ncol(inputProteins))
  # gene_name = rep(make.names(intensities[,2],unique=T),ncol(inputProteins))
  #
  # The old version assumed a specific row order after pivoting.
  # This version keeps the protein identity attached to each row before pivoting.
  intensity_long <- intensities %>%
    dplyr::mutate(
      proteinRank = dplyr::row_number(),
      proteinId = .data[[colnames(intensities)[1]]],
      Gene_name = make.names(.data[[colnames(intensities)[2]]], unique = TRUE)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(intensity_cols),
      names_to = "Sample",
      values_to = "log_intensity"
    ) %>%
    dplyr::mutate(
      log_intensity = as.numeric(log_intensity)
    )
  
  
  ############################
  # 4. Convert imputed mask to long format
  ############################
  
  # Your imputed mask columns should contain "mean".
  # Example:
  # K13D_1_C3_K13D_mean, K13D_2_C3_K13D_mean, ...
  mask_cols <- grep("mean$", colnames(imputed), value = TRUE)
  
  if (length(mask_cols) == 0) {
    stop("create_heatmap(): imputed mask has no columns ending in '_mean'.")
  }
  
  # Check that the mask sample columns correspond to the heatmap sample columns.
  # The heatmap uses:
  # K13D_1_C3_K13D
  #
  # The mask uses:
  # K13D_1_C3_K13D_mean
  #
  # After removing "_mean", they should be identical.
  if (!identical(intensity_cols, sub("_mean$", "", mask_cols))) {
    stop(
      "create_heatmap(): heatmap sample columns and imputed mask columns do not match after removing '_mean'.\n",
      "Heatmap columns: ", paste(intensity_cols, collapse = ", "), "\n",
      "Mask columns: ", paste(mask_cols, collapse = ", ")
    )
  }
  
  # This replaces:
  #
  # pivot_long_imputed <- imputed %>%
  #   pivot_longer(cols=contains("mean"),
  #                names_to="condition",
  #                values_to="is_imputed")
  #
  # The important addition is proteinRank.
  # That gives us a stable way to join the mask to the intensity values.
  imputed_long <- imputed %>%
    dplyr::mutate(
      proteinRank = dplyr::row_number()
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(mask_cols),
      names_to = "Sample",
      values_to = "imputed"
    ) %>%
    dplyr::mutate(
      Sample = sub("_mean$", "", Sample),
      imputed = as.logical(imputed)
    )
  
  
  ############################
  # 5. Join intensity values and imputation mask
  ############################
  
  # This replaces:
  #
  # heatmapDF = data.frame(
  #   proteinId,
  #   gene_name,
  #   as.character(pivot_long_intensitiesDF$condition),
  #   pivot_long_intensitiesDF$intensity,
  #   proteinRank,
  #   pivot_long_imputed$is_imputed
  # )
  #
  # The old version put columns together by position.
  # The new version joins them by proteinRank and Sample.
  # This avoids misalignment when pivot_longer ordering is not what we assumed.
  heatmapDF <- intensity_long %>%
    dplyr::left_join(
      imputed_long[, c("proteinRank", "Sample", "imputed"), drop = FALSE],
      by = c("proteinRank", "Sample")
    )
  
  # If this happens, the mask did not match some heatmap cells.
  # That would make blank spaces unreliable, so we stop immediately.
  if (any(is.na(heatmapDF$imputed))) {
    stop("create_heatmap(): some heatmap cells do not have a matching imputation-mask value.")
  }
  
  # This is the line that creates blank spaces.
  # Imputed cells are removed from the plotted dataframe.
  # geom_raster() then leaves those positions blank.
  heatmapDF <- subset(heatmapDF, imputed == FALSE)
  
  if (nrow(heatmapDF) == 0) {
    stop("create_heatmap(): no non-imputed heatmap cells remain after filtering.")
  }
  
  
  ############################
  # 6. Order samples by group
  ############################
  
  # This replaces the original orderedSamples block.
  #
  # Important change:
  # I keep ordered_sample_levels separate from heatmapDF$orderedSamples.
  # In the old function, orderedSamples was first a sample list, then later reused
  # as a factor the length of the full heatmap dataframe.
  ordered_sample_levels <- c()
  
  group_ids <- sort(group_ids)
  
  for (group in group_ids) {
    ordered_sample_levels <- c(
      ordered_sample_levels,
      sort(grep(paste0("_", group), samples_names, value = TRUE))
    )
  }
  
  ordered_sample_levels <- unique(ordered_sample_levels)
  
  # Safety fallback:
  # If any intensity column was not captured by the group-based ordering,
  # append it at the end instead of creating NA x positions.
  missing_ordered_samples <- setdiff(intensity_cols, ordered_sample_levels)
  
  if (length(missing_ordered_samples) > 0) {
    ordered_sample_levels <- c(ordered_sample_levels, missing_ordered_samples)
  }
  
  heatmapDF$orderedSamples <- factor(
    heatmapDF$Sample,
    levels = ordered_sample_levels
  )
  
  heatmapDF$x <- as.numeric(heatmapDF$orderedSamples)
  
  if (any(is.na(heatmapDF$x))) {
    stop("create_heatmap(): some samples could not be assigned an x-axis position.")
  }
  
  
  ############################
  # 7. Order proteins according to clustering
  ############################
  
  # This keeps the original idea:
  # proteins are displayed in the order from hclust().
  heatmapDF$orderedProteins <- factor(
    heatmapDF$proteinRank,
    levels = hClustProteins$order
  )
  
  heatmapDF$y <- as.numeric(heatmapDF$orderedProteins)
  
  if (any(is.na(heatmapDF$y))) {
    stop("create_heatmap(): some proteins could not be assigned a y-axis position.")
  }
  
  # Labels in clustered order.
  # This is safer than using unique(heatmapDF$Gene_name), because heatmapDF has
  # multiple rows per protein and some cells have been removed.
  y_axis_labels <- make.names(intensities[[2]], unique = TRUE)[hClustProteins$order]
  
  
  ############################
  # 8. Adjust color limits
  ############################
  
  # These were in the original function.
  # They are currently not used directly in scale_fill_gradientn(),
  # but I am keeping them in case you use them later.
  minRatio <- min(heatmapDF$log_intensity, na.rm = TRUE)
  maxRatio <- max(heatmapDF$log_intensity, na.rm = TRUE)
  maxabsRatio <- max(abs(heatmapDF$log_intensity), na.rm = TRUE)
  
  
  ############################
  # 9. Make heatmap plot
  ############################
  
  hcPlot <- ggplot(data = heatmapDF) +
    geom_raster(mapping = aes(x = x, y = y, fill = log_intensity)) +
    scale_fill_gradientn(colours = c("green", "yellow", "red")) +
    ggtitle("") +
    scale_x_continuous(
      name = "",
      breaks = seq_along(ordered_sample_levels),
      labels = ordered_sample_levels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = "",
      breaks = seq_along(hClustProteins$order),
      labels = y_axis_labels,
      expand = c(0, 0),
      position = "right"
    ) +
    labs(title = htitle) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 30),
      axis.title = element_text(size = 25),
      axis.text.x = element_text(size = 24, angle = 90),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    )
  
  
  ############################
  # 10. Make protein dendrogram
  ############################
  
  dendrogramProteins <- dendro_data(hClustProteins, type = "rectangle")
  dendrogramDataProteins <- segment(dendrogramProteins)
  
  dendrogramProteinsPlot <- ggplot() +
    geom_segment(
      data = dendrogramDataProteins,
      mapping = aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    coord_flip(clip = "off") +
    scale_y_reverse(expand = c(0.5, 0)) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  
  ############################
  # 11. Combine dendrogram and heatmap
  ############################
  
  result <- dendrogramProteinsPlot +
    plot_spacer() +
    hcPlot +
    plot_layout(
      widths = c(4, -1.1, 4.5)
    )
  
  # These were in your original function.
  result$heights[7] <- unit(0.2, "null")
  result$widths[5] <- unit(0.2, "null")
  
  return(result)
}

########################
#Heatmap w/o clustering#
########################

# This function create a heatmap with a hierarchical clustering of proteins
create_heatmap_2 <- function(intensities,imputed,samples_names,group_ids){
  
  # Input :
  # intensities : quantification matrix whom 1rst column contains protein id
  # imputed : boolean matrix whom 1rst column contains protein id
  # samples_names : list of sample names
  # group_ids : list of biological conditions
  # Output :
  # heatmap plot
  
  inputProteins = intensities[,3:ncol(intensities)]
  rownames(inputProteins) = make.names(intensities$Label,unique=T)
  distMatrixProteins <- dist(inputProteins)
  hClustProteins <- hclust(distMatrixProteins)
  
  # Create data frame
  proteinRank = c()
  for (i in 1:(ncol(inputProteins))){
    proteinRank = c(proteinRank,1:nrow(inputProteins))
  }
  proteinId = rep(intensities[,1],ncol(inputProteins))
  gene_name = rep(make.names(intensities[,2],unique=T),ncol(inputProteins))
  intensitiesDF = melt(inputProteins)
  imputed = melt(imputed,id.vars = 1)
  heatmapDF = cbind(proteinId,gene_name,intensitiesDF,proteinRank,imputed$value)
  colnames(heatmapDF) = c("Id","Gene_name","Sample","log_intensity","proteinRank","imputed")
  
  heatmapDF=subset(heatmapDF,heatmapDF$imputed==F)
  
  # Order sample by group
  orderedSamples = c()
  group_ids = sort(group_ids)
  for (group in group_ids){
    orderedSamples = c(orderedSamples,sort(grep(paste0("_",group),samples_names,value=T)))
  }
  heatmapDF$orderedSamples <- factor(heatmapDF$Sample, levels = unique(orderedSamples))
  heatmapDF$x <- as.numeric(heatmapDF$orderedSamples)
  
  # Order proteins according to clustering results
  heatmapDF$orderedProteins <- factor(heatmapDF$proteinRank, levels = hClustProteins$order)
  heatmapDF$y <- as.numeric(heatmapDF$orderedProteins)
  
  # Adjust colors
  minRatio <- min(heatmapDF$log_intensity)
  maxRatio <- max(heatmapDF$log_intensity)
  maxabsRatio <- max(abs(heatmapDF$log_intensity))
  
  # Transforming Heatmap DataFrame
  # heatmapDF <- t(heatmapDF)
  # write.csv(heatmapDF, "heatmapDF")
  
  # Make plot
  hcPlot <- ggplot(data = heatmapDF) +
    geom_raster(mapping = aes(x = y, y = x, fill = log_intensity)) +
    scale_fill_gradientn(colours=c("green","yellow","red")) +
    #scale_fill_scico(palette = "berlin") +
    ggtitle("")+
    scale_y_continuous(name = "", breaks = 1:length(orderedSamples), labels = orderedSamples, expand = c(0, 0), position = "right") +
    scale_x_continuous(name = "", breaks = unique(heatmapDF$y), labels = unique(heatmapDF$Gene_name), expand = c(0, 0)) +
    theme(#axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title=element_text(size=30),
      axis.title=element_text(size=25),
      axis.text.y=element_text(size=24),
      axis.text.x=element_text(angle = 90),
      legend.text=element_text(size=16),
      legend.title=element_text(size=18))
  
  return(hcPlot)
  
  # No dendrogram needed
  
  # Make grobs from plots
  # matrixGrob <- ggplotGrob(hcPlot)
  
  # Assemble
  # bottom
}

#######
# ROC #
#######

# This function compute the sensibility and the specificity varying ratio and pvalue thresholds
compute_ROC <- function(statistic_table,n,ratioVar,pvalVar,variants){
  
  # Input :
  # statistic_table : table containing results of statistical analysis and whom 1rst column corresponds to protein ids
  # n : number of tests
  # ratioVar : ratio thresholds to test
  # pvalVar : pval thresholds to test
  # variants : list of ids of expected variants
  # Output :
  # heatmap plot
  
  mat <- matrix(ncol = 4, nrow = (n*n))
  dimnames(mat)[[2]] <- c("FDP", "Sensitivity", "ratio_threshold", "pval_threshold")
  l <- 0
  for (i in 1:n) {
    Z <- ratioVar[i]
    for (j in 1:n) {
      p <- pvalVar[j]
      
      FP <- nrow(statistic_table[(abs(statistic_table$ratio) >= Z) & (statistic_table$pval <= p) & !(statistic_data$Accession %in% variants),])
      TP <- nrow(statistic_table[(abs(statistic_table$ratio) >= Z) & (statistic_table$pval <= p) & (statistic_data$Accession %in% variants),])
      TN <- nrow(statistic_table[((abs(statistic_table$ratio) < Z) | (statistic_table$pval > p)) & !(statistic_data$Accession %in% variants),])
      FN <- nrow(statistic_table[((abs(statistic_table$ratio) < Z) | (statistic_table$pval > p)) & (statistic_data$Accession %in% variants),])
      
      speci <- (FP/(FP+TP))*100
      sensi <- (TP/(TP+FN))*100
      mat[j+l,] <- c(speci, sensi, Z, p)
    }
    l <- n+l
  }
  gtab <- as.data.frame(mat)
  return(gtab)
  
}