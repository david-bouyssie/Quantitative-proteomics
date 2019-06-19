### See:
### - https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/
### - https://stackoverflow.com/questions/28507693/call-rmarkdown-on-command-line-using-a-r-that-is-passed-a-file
### - https://stackoverflow.com/questions/32479130/passing-parameters-to-r-markdown

args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")
output_dir <- "results/"
rmarkdown::render('../Modules/Instrumental QC/Descriptive_cytoc.Rmd',
                  output_file = args[2],
                  params=list(file_name=args[1]))
