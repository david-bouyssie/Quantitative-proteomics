.libPaths(
  c(
    .libPaths(),
    "./my_project/renv/library/R-3.6/x86_64-w64-mingw32",
    "./my_project/library"
  )
)

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

source(paste0("./my_project/Quantitative-proteomics/Parsers/parse_",source,".R"))
