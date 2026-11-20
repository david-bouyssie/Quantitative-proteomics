# Repository root (Quantitative-proteomics), found from the location of this launcher, so that it runs from any current directory.
# With Rscript, the script path is given by --file=. When sourced from an R session, the working directory must be the Launchers folder.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) == 1) {
  launcher_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else if (basename(getwd()) == "Launchers" && file.exists("run_proteomics_stats.R")) {
  launcher_dir <- getwd()
} else {
  stop("Run this launcher with Rscript, or source it with the Launchers folder as working directory.")
}
repo_root <- dirname(launcher_dir)

# Package library of the renv project (see env/renv.lock), placed first so that its packages take precedence
# over the copies in the R system library
renv_library <- file.path(repo_root, "env", "renv", "library", "R-4.2", "x86_64-w64-mingw32")
if (!dir.exists(renv_library)) {
  stop("Package library not found: ", renv_library)
}
.libPaths(c(renv_library, .libPaths()))
modules_dir <- file.path(repo_root, "Modules", "Quantitative data analysis")
# Input and output paths given as arguments are relative to the current directory

### PATCH WINDOWS FONTS (see: https://stackoverflow.com/questions/55933524/r-can-not-find-fonts-to-be-used-in-plotting) ###
if (.Platform$OS.type == "windows") {
  library(showtext)
  font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
  showtext_auto()
  # showtext draws the texts at its own resolution (96 dpi by default): the figures are designed at 72 pixels per inch (1 pt = 1 px,
  # see plot_figure_as_html in the library), otherwise their texts are a third larger than planned and the long ones are cut
  showtext_opts(dpi = 72)
}

library(stringr)

args <- commandArgs(trailingOnly = TRUE)
analysis = args[1]
quantif_file = normalizePath(args[2], winslash = "/")

if (analysis == "QC") {

  output_file = file.path(normalizePath(dirname(args[3]), winslash = "/"), basename(args[3]))
  parameters = normalizePath(args[4], winslash = "/")
  # Optional original quantification file (absent, or "" to give an experimental design without it), then optional experimental design
  if (length(args) >= 5 && nzchar(args[5])){
    origin_file = normalizePath(args[5], winslash = "/")
  }else{
    origin_file = ""
  }
  if (length(args) >= 6){
    exp_design = normalizePath(args[6], winslash = "/", mustWork = FALSE)  # QC without experimental design if the file does not exist
  }else{
    exp_design = ""
  }

  rmarkdown::render(
    file.path(modules_dir, "QC.Rmd"),
    output_file = output_file,
    params = list(
      quantif_file = quantif_file,
      design_exp = exp_design,
      parameters = parameters,
      origin_file = origin_file
    )
  )

} else if (analysis=="DA") {

  dir.create(args[3], recursive = TRUE, showWarnings = FALSE)  # the output folder is created if needed
  output_dir = normalizePath(args[3], winslash = "/")
  exp_design = normalizePath(args[4], winslash = "/")
  parameters = normalizePath(args[5], winslash = "/")
  if (length(args) >= 6){
    origin_file = normalizePath(args[6], winslash = "/")
  }else{
    origin_file = ""
  }
  output_file = file.path(output_dir, "DA.html")
  
  rmarkdown::render(
    file.path(modules_dir, "DA.rmd"),
    
    output_file = output_file ,
    params = list(
      quant_df = quantif_file,
      parameters = parameters,
      design_exp = exp_design,
      output_dir = output_dir,
      origin_file = origin_file
    )
  )

} else {
  print("Wrong argument")
}
