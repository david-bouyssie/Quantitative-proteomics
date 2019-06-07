
args <- commandArgs(trailingOnly = TRUE)
repository = paste0("../",args[1])
variants = paste0("../",args[2])
output_file = paste0("../",args[3])
nb_tests = as.numeric(args[4])
min_ratio = as.numeric(args[5])
max_ratio = as.numeric(args[6])
default_ratio_threshold = as.numeric(args[7])
min_pval = as.numeric(args[8])
max_pval = as.numeric(args[9])
default_pval_threshold = as.numeric(args[10])

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")

normalization = args[4]
rmarkdown::render('Rmd/ROC.Rmd',
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
