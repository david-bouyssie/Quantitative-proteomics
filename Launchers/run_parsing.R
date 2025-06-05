.libPaths(
  c(
    .libPaths(),
    "./renv/library/R-4.2/x86_64-w64-mingw32"
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

source(paste0("./Quantitative-proteomics/Parsers/parse_",source,".R"))
#source(paste0("./Parsers/parse_",source,".R"))
