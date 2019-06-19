
:: ###################
:: # Parsing modules #
:: ###################

:: Command line format : Rscript --vanilla run_parsing.R [file to parse] [output] [parser option : Proline/Maxquant/ProteomeDiscoverer] [abundance type option : LFQ/iBAQ/raw_abundance/abundance] [if xlsx file, sheet number]
:: Example : 
::Rscript --vanilla run_parsing.R "path_to_file_to_parse/name_of_file_to_parse" "path_to_desired_output_directory/desired_name" Proline abundance 4


:: ###############################
:: # Quantitative data QC module #
:: ###############################

:: Command line format : Rscript --vanilla run_proteomics_stats.R QC [input (standardized format)] [output] [normalization option: TRUE/FALSE] [delete empty lines : TRUE/FALSE]
:: Example :
::Rscript --vanilla run_proteomics_stats.R QC "../../Example/Datasets/QC-DA/Parsed_proteins_set.txt" "../../Example/Reports/QC_example.html" TRUE FALSE


:: ################################
:: # Differential analysis module #
:: ################################

:: Command line format : Rscript --vanilla run_proteomics_stats.R DA [input (standardized format)] [report output] [table output] [experimental design file] [parameters file]
:: Example :
::Rscript --vanilla run_proteomics_stats.R DA "../../Example/Datasets/QC-DA/Parsed_proteins_set.txt" "../../Example/Reports/DA_example.html" "../../Example/Reports/DA_example.txt" "../../Example/Datasets/QC-DA/designExp.txt" "../../Example/Datasets/QC-DA/Parameters.txt"


:: ##############
:: # ROC module #
:: ##############

:: Command line format : Rscript --vanilla run_proteomics_stats.R DA [input folder] [variant list] [output report] [number of tests] [min ratio] [max ratio] [default ratio] [min pval] [max pval] [default pval]
:: Example :
::Rscript --vanilla run_ROC.R "../../Example/Datasets/ROC/" "../../Example/Datasets/ROC/variants/variants.csv" "../../Example/Reports/ROC_example.html" 10 0.5 10 0.58 0 0.05 0.05


:: ############
:: # CytoC QC #
:: ############

:: Command line format : Rscript --vanilla run_proteomics_stats.R DA [input] [output folder]
:: Example :
::Rscript --vanilla run_cytoc.R "../../Example/Datasets/CytoC/cytoc_xics_test.csv" "../../Example/Reports"