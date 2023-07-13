.libPaths(
  c(
    .libPaths(),
    "./Quantitative-proteomics/renv/library/R-4.2/x86_64-w64-mingw32"
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

if (analysis == "QC") {

  output_file = file.path(normalizePath(dirname(args[3])), basename(args[3]))
  parameters = normalizePath(args[4])
  exp_design = normalizePath(args[5])

  rmarkdown::render(
    './Quantitative-proteomics/Modules/Quantitative data analysis/QC.Rmd',
    output_file = output_file,
    params = list(
      quantif_file = quantif_file,
      design_exp = exp_design,
      parameters = parameters 
    )
  )

} else if (analysis=="DA") {

  output_dir = normalizePath(args[3])
  exp_design = normalizePath(args[4])
  parameters = normalizePath(args[5])
  output_file = file.path(output_dir, "DA.html")
  
  rmarkdown::render(
    './Quantitative-proteomics/Modules/Quantitative data analysis/DA.Rmd',
    output_file = output_file ,
    params = list(
      quant_df = quantif_file,
      parameters = parameters,
      design_exp = exp_design,
      output_dir = output_dir,
      proline_source_file = "none",
      proline_sheet = 4
    )
  )

} else {
  print("Wrong argument")
}
