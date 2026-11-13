# Repository root (Quantitative-proteomics), found from the location of this launcher, so that it runs from any current directory.
# With Rscript, the script path is given by --file=. When sourced from an R session, the working directory must be the Launchers folder.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) == 1) {
  launcher_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else if (basename(getwd()) == "Launchers" && file.exists("run_ROC.R")) {
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
## Rscript --vanilla <path to>/Launchers/run_ROC.R [input folder] [variant list] [output report] [number of tests] [min ratio] [max ratio] [default ratio] [min pval] [max pval] [default pval]

args <- commandArgs(trailingOnly = TRUE)
# The report is knitted in its own folder: input and output paths are made absolute
repository = paste0(normalizePath(args[1]), "/")
variants = normalizePath(args[2])
output_file = file.path(normalizePath(dirname(args[3])), basename(args[3]))
nb_tests = as.numeric(args[4])
min_ratio = as.numeric(args[5])
max_ratio = as.numeric(args[6])
default_ratio_threshold = as.numeric(args[7])
min_pval = as.numeric(args[8])
max_pval = as.numeric(args[9])
default_pval_threshold = as.numeric(args[10])

rmarkdown::render(file.path(repo_root, "Modules", "Quantitative data analysis", "ROC.rmd"),
                  output_file = output_file,
                  params=list(
                    repository = repository,
                    variants = variants,
                    nb_tests = nb_tests,
                    min_ratio = min_ratio,
                    max_ratio = max_ratio,
                    default_ratio_threshold = default_ratio_threshold,
                    min_pval = min_pval,
                    max_pval = max_pval,
                    default_pval_threshold = default_pval_threshold))
