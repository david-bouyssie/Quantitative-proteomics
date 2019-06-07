library(stringr)

args <- commandArgs(trailingOnly = TRUE)
analysis = args[1]
quantif_file = paste0("../",args[2])
output_file = paste0("../",args[3])

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")

if(analysis=="QC"){

  normalization = args[4]
  keep_empty_rows = args[5]
  rmarkdown::render('Rmd/QC.Rmd',
                    #output_dir="../Results",
                    output_file = output_file,
                    params=list(
                      quantif_file=quantif_file,
                      normalization=normalization,
                      keep_empty_rows=keep_empty_rows))


}else if(analysis=="DA"){

  output_table = paste0("../",args[4])
  design_exp = paste0("../",args[5])
  parameters = paste0("../",args[6])
  rmarkdown::render('Rmd/DA.Rmd',
                    #output_dir="../Results",
                    output_file = output_file,
                    params=list(
                      quantif_file=quantif_file,
                      parameters=parameters,
                      design_exp=design_exp,
                      output_table=output_table))

}else{
  print("Wrong argument")
}
