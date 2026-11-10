# Package library of the renv project (see Quantitative-proteomics/env/renv.lock), placed first so that its packages take precedence
# over the copies in the R system library. Paths are relative to R-Portable_with_QP (set by run_script.bat).
renv_library <- "./Quantitative-proteomics/env/renv/library/R-4.2/x86_64-w64-mingw32"
if (!dir.exists(renv_library)) {
  stop("Package library not found: ", normalizePath(renv_library, mustWork = FALSE))
}
.libPaths(c(renv_library, .libPaths()))

## Usage :
## Rscript --vanilla run_parsing.R quantification_file output quantification_file_source intensity_type sheet

#library(here)
#setwd(here())

args <- commandArgs(trailingOnly = TRUE)
quantif_file = args[1]
output_file = args[2]
source = args[3]
intensities_type = args[4]
sheet = as.numeric(args[5])

source(paste0("./Quantitative-proteomics/Parsers/parse_",source,".R"))
#source(paste0("./Parsers/parse_",source,".R"))
