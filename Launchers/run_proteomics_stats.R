.libPaths(
  c(
    .libPaths(),
    "./my_project/renv/library/R-3.6/x86_64-w64-mingw32",
    "./my_project/library"
  )
)

### PATCH WINDOWS FONTS (see: https://stackoverflow.com/questions/55933524/r-can-not-find-fonts-to-be-used-in-plotting) ###
if (.Platform$OS.type == "windows") {
  library(showtext)
  font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
  showtext_auto()
}

library(stringr)

args <- commandArgs(trailingOnly = TRUE)
analysis = args[1]
quantif_file = normalizePath(args[2])
output_file = file.path(normalizePath(dirname(args[3])), basename(args[3]))
print(output_file)

if (analysis == "QC") {

  normalization = args[4]
  keep_empty_rows = args[5]
  
  rmarkdown::render(
    './my_project/Quantitative-proteomics/Modules/Quantitative data analysis/QC.Rmd',
    output_file = output_file,
    params = list(
      quantif_file = quantif_file,
      normalization = normalization,
      keep_empty_rows = keep_empty_rows
    )
  )

} else if (analysis=="DA") {

  output_table = file.path(normalizePath(dirname(args[4])), basename(args[4]))
  design_exp = normalizePath(args[5])
  parameters = normalizePath(args[6])
  
  rmarkdown::render(
    './my_project/Quantitative-proteomics/Modules/Quantitative data analysis/DA.Rmd',
    output_file = output_file,
    params = list(
      quantif_file = quantif_file,
      parameters = parameters,
      design_exp = design_exp,
      output_table = output_table
    )
  )

} else {
  print("Wrong argument")
}
