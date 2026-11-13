# Repository root (Quantitative-proteomics), found from the location of this launcher, so that it runs from any current directory.
# With Rscript, the script path is given by --file=. When sourced from an R session, the working directory must be the Launchers folder.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) == 1) {
  launcher_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else if (basename(getwd()) == "Launchers" && file.exists("run_cytoc.R")) {
  launcher_dir <- getwd()
} else {
  stop("Run this launcher with Rscript, or source it with the Launchers folder as working directory.")
}
repo_root <- dirname(launcher_dir)

# Package library of the renv project (see env/renv.lock), placed first so that its packages take precedence
# over the copies in the R system library. pandoc is looked for on the PATH by rmarkdown.
renv_library <- file.path(repo_root, "env", "renv", "library", "R-4.2", "x86_64-w64-mingw32")
if (!dir.exists(renv_library)) {
  stop("Package library not found: ", renv_library)
}
.libPaths(c(renv_library, .libPaths()))

## Usage (input and output paths are relative to the current directory):
## Rscript --vanilla <path to>/Launchers/run_cytoc.R [cytoc_xics file] [output report]

args <- commandArgs(trailingOnly = TRUE)
# The report is knitted in its own folder: input and output paths are made absolute
input_file = normalizePath(args[1])
output_file = file.path(normalizePath(dirname(args[2])), basename(args[2]))

rmarkdown::render(file.path(repo_root, "Modules", "Instrumental QC", "Descriptive_cytoc.rmd"),
                  output_file = output_file,
                  params=list(file_name=input_file))
