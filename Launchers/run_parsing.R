# Repository root (Quantitative-proteomics), found from the location of this launcher, so that it runs from any current directory.
# With Rscript, the script path is given by --file=. When sourced from an R session, the working directory must be the Launchers folder.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) == 1) {
  launcher_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else if (basename(getwd()) == "Launchers" && file.exists("run_parsing.R")) {
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

## Usage (input and output paths are relative to the current directory):
## Rscript --vanilla <path to>/Launchers/run_parsing.R quantification_file output quantification_file_source intensity_type sheet

args <- commandArgs(trailingOnly = TRUE)
quantif_file = args[1]
output_file = args[2]
source = args[3]
intensities_type = args[4]
sheet = as.numeric(args[5])

source(file.path(repo_root, "Parsers", paste0("parse_", source, ".R")))
