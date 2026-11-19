library(readxl)

# Parameters set in launcher script (see "Lauchers" section)
# quantif_file = quantification file to parse (Proline export: xlsx, or tab-separated text)
# output_file = output table
# intensities_type = abundance (normalized by Proline) or raw_abundance
# sheet = sheet number if quantif_file is a xlsx file, used only if there is no "Protein sets" sheet
# ibaq = TRUE to convert the intensities to iBAQ (intensity divided by the number of observable peptides)

# Load data, keeping the column names as they are (sample names may contain spaces or other characters)
if(grepl("\\.xlsx?$", quantif_file, ignore.case = TRUE)){
  # The "Protein sets" sheet is found by its name, its position depends on the sheets selected in the export
  sheets = excel_sheets(quantif_file)
  protein_sets_sheet = which(tolower(gsub("[ _]", "", sheets)) == "proteinsets")
  if (length(protein_sets_sheet) == 1) {
    sheet = protein_sets_sheet
  } else if (is.na(sheet) || sheet < 1 || sheet > length(sheets)) {
    stop("Proline parser: no \"Protein sets\" sheet found (sheets: ", paste(sheets, collapse = ", "), ") and no valid sheet number given.")
  }
  cat("Proline export: sheet", sheet, "(", sheets[sheet], ")\n")
  proteinGroupsInput = as.data.frame(read_excel(quantif_file, sheet, col_names = TRUE), check.names = FALSE)
}else{
  # quote = "" and comment.char = "": quotes or "#" in protein descriptions must not shift or drop rows
  proteinGroupsInput <- read.delim(quantif_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
}
cat("Proline export:", nrow(proteinGroupsInput), "rows,", ncol(proteinGroupsInput), "columns\n")

# Returns the columns whose names start with the given prefix, as a data frame (column names without the prefix)
prefixed_columns <- function(prefix){
  cols <- names(proteinGroupsInput)[startsWith(names(proteinGroupsInput), prefix)]
  out <- proteinGroupsInput[cols]
  names(out) <- substring(cols, nchar(prefix) + 1)
  out
}

# Id: the protein set id
id_col <- if ("protein_set_id" %in% names(proteinGroupsInput)) "protein_set_id" else grep("id$", names(proteinGroupsInput), value = TRUE)[1]
if (is.na(id_col)) stop("Proline parser: no protein set id column (protein_set_id) found.")
Id = data.frame(Id = proteinGroupsInput[[id_col]])
if (anyDuplicated(Id$Id)) stop("Proline parser: the ", id_col, " column has duplicated values.")

if (!"accession" %in% names(proteinGroupsInput)) stop("Proline parser: no accession column found.")
Accession = data.frame(Accession = proteinGroupsInput[["accession"]])

# Gene names (missing values are allowed; the column is created empty if absent)
Gene_name = data.frame(Gene_name = if ("gene_name" %in% names(proteinGroupsInput)) proteinGroupsInput[["gene_name"]] else NA_character_)

# Intensities: "<intensities_type>_<sample>" columns; zeros are missing values
intensities = prefixed_columns(paste0(intensities_type, "_"))
if (ncol(intensities) == 0) {
  available = c("abundance", "raw_abundance")[sapply(c("abundance_", "raw_abundance_"), function(p) any(startsWith(names(proteinGroupsInput), p)))]
  stop("Proline parser: no '", intensities_type, "_<sample>' column found (intensity types available: ",
       if (length(available)) paste(available, collapse = ", ") else "none", "). Check the intensity type and the sheet number.")
}
samples = names(intensities)
intensities[] = lapply(intensities, function(x) { x = as.numeric(x); x[!is.na(x) & x <= 0] = NA; x })
names(intensities) = paste0("Intensity_", samples)

# iBAQ: intensities divided by the number of observable peptides. Proline fills the #observable_peptides column (-1 otherwise)
# only after the task retrieving the protein sequences from the FASTA files.
if (ibaq) {
  if (!"#observable_peptides" %in% names(proteinGroupsInput)) {
    stop("Proline parser: iBAQ conversion impossible, no #observable_peptides column found.")
  }
  observable = as.numeric(proteinGroupsInput[["#observable_peptides"]])
  valid = !is.na(observable) & observable > 0
  if (!any(valid)) {
    stop("Proline parser: iBAQ conversion impossible, the number of observable peptides is unknown for all the proteins (#observable_peptides = -1). ",
         "Retrieve the protein sequences from the FASTA files in Proline, then export the dataset again.")
  }
  if (any(!valid)) {
    cat("Warning: number of observable peptides unknown for", sum(!valid), "proteins (", paste(head(Accession$Accession[!valid], 5), collapse = ", "),
        if (sum(!valid) > 5) "..." , "), their iBAQ intensities are missing.\n")
  }
  intensities[] = lapply(intensities, function(x) x / ifelse(valid, observable, NA))
}

# Identification types from the "psm_count_<sample>" columns: at least one PSM = "By MS/MS", none = "By matching", empty = not observed
psm_counts = prefixed_columns("psm_count_")
missing_samples = setdiff(samples, names(psm_counts))
if (length(missing_samples) > 0) {
  stop("Proline parser: no psm_count column for sample(s) ", paste(missing_samples, collapse = ", "), ".")
}
identification_types = as.data.frame(lapply(psm_counts[samples], function(p) {
  p = as.numeric(p)
  ifelse(is.na(p), NA, ifelse(p > 0, "By MS/MS", "By matching"))
}), stringsAsFactors = FALSE)
names(identification_types) = paste0("Identification_type_", samples)

# Number of specific peptides (the DA report removes the proteins without specific peptide)
if (!"#specific_peptides" %in% names(proteinGroupsInput)) {
  stop("Proline parser: no #specific_peptides column found (it is part of the Proline protein sets export).")
}
specific_peptides = data.frame(Specific_peptides = as.numeric(proteinGroupsInput[["#specific_peptides"]]))

# Build dataframe
proteinGroupsOutput = cbind(Id, Accession, Gene_name, identification_types, specific_peptides, intensities)
write.table(proteinGroupsOutput, output_file, row.names = FALSE, sep = "\t")
cat("Parsed:", nrow(proteinGroupsOutput), "proteins,", length(samples), "samples (", intensities_type, if (ibaq) "converted to iBAQ", "intensities), written to", output_file, "\n")
