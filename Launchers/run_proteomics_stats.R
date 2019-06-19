library(stringr)

args <- commandArgs(trailingOnly = TRUE)
analysis = args[1]
quantif_file = args[2]
output_file = args[3]

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")

if(analysis=="QC"){

  normalization = args[4]
  keep_empty_rows = args[5]
  rmarkdown::render('../Modules/Quantitative data analysis/QC.Rmd',
                    output_file = output_file,
                    params=list(
                      quantif_file=quantif_file,
                      normalization=normalization,
                      keep_empty_rows=keep_empty_rows))


}else if(analysis=="DA"){

  output_table = args[4]
  print(quantif_file)
  design_exp = args[5]
  parameters = args[6]
  rmarkdown::render('../Modules/Quantitative data analysis/DA.Rmd',
                    output_file = output_file,
                    params=list(
                      quantif_file=quantif_file,
                      parameters=parameters,
                      design_exp=design_exp,
                      output_table=output_table))

}else{
  print("Wrong argument")
}
