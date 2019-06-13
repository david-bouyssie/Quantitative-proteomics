@echo off
set PATH=%PATH%;C:\Program Files\R\R-3.5.1\bin


:: ###################
:: # Parsing modules #
:: ###################

:: Command line format : Rscript --vanilla run_parsing.R [file to parse] [output] [parser option : Proline/Maxquant/ProteomeDiscoverer] [abundance type option : LFQ/iBAQ/raw_abundance/abundance] [if xlsx file, sheet number]
:: Example : 
Rscript --vanilla run_parsing.R "../Data/Search Xtandem PME12modified Percolator validation_MRF_uniquepep_discardMC-pepnorm_16-05-2019_0941.xlsx" "../Data/Quanti/Parsed_data/Search Xtandem PME12modified Percolator validation_MRF_uniquepep_discardMC-pepnorm_16-05-2019_0941.txt" Proline abundance 4


:: ###############################
:: # Quantitative data QC module #
:: ###############################

:: Command line format : Rscript --vanilla run_proteomics_stats.R QC [input (standardized format)] [output] [normalization option: TRUE/FALSE] [delete empty lines : TRUE/FALSE]
:: Example :
Rscript --vanilla run_proteomics_stats.R QC ../Data/standardized_file.txt /Results/QC/test.html TRUE FALSE


:: ################################
:: # Differential analysis module #
:: ################################

:: Command line format : Rscript --vanilla run_proteomics_stats.R DA [input (standardized format)] [report output] [table output] [experimental design file] [parameters file]
:: Example :
Rscript --vanilla run_proteomics_stats.R DA ../Data/standardized_file.txt /Results/DA/test.html /Results/DA/test.txt ../Data/designExp.txt ../Data/Params_da.txt


:: ##############
:: # ROC module #
:: ##############

:: Command line format : Rscript --vanilla run_proteomics_stats.R DA [input folder] [variant list] [output report] [number of tests] [min ratio] [max ratio] [default ratio] [min pval] [max pval] [default pval]
:: Example :
Rscript --vanilla run_ROC.R ../Results/imputation/ ../Results/variants.csv ../Results/imputation/comparison.html 10 0.5 10 0.58 0 0.05 0.05


:: ############
:: # CytoC QC #
:: ############

:: Command line format : Rscript --vanilla run_proteomics_stats.R DA [input] [output folder]
:: Example :
Rscript --vanilla run_cytoc.R ../../Data/CytoC/cytoc_xics.csv ../Results/CytoC/
