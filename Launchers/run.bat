::@echo off
::set PATH=%PATH%;C:/Program Files/R/R-3.5.1/bin


:: ###################
:: # Parsing modules #
:: ###################

:: Command line format : Rscript --vanilla run_parsing.R [file to parse] [output] [parser option : Proline/Maxquant/ProteomeDiscoverer] [abundance type option : LFQ/iBAQ/raw_abundance/abundance] [if xlsx file, sheet number]
:: Example : 
::Rscript --vanilla run_parsing.R "D:/Utilisateurs/cludwig/stage/RStudio/Tests/Data/Quanti/Raw/Proline/Sandrine/180619_PCP Analysis from DEAE sample inj 1 SLB_18-06-2019_0931 (1).xlsx" "D:/Utilisateurs/cludwig/stage/RStudio/Tests/test.txt" Proline abundance 4


:: ###############################
:: # Quantitative data QC module #
:: ###############################

:: Command line format : Rscript --vanilla run_proteomics_stats.R QC [input (standardized format)] [output] [normalization option: TRUE/FALSE] [delete empty lines : TRUE/FALSE]
:: Example :
::Rscript --vanilla run_proteomics_stats.R QC "D:/Utilisateurs/cludwig/stage/RStudio/Tests/test.txt" "D:/Utilisateurs/cludwig/stage/RStudio/Tests/test.html" TRUE FALSE


:: ################################
:: # Differential analysis module #
:: ################################

:: Command line format : Rscript --vanilla run_proteomics_stats.R DA [input (standardized format)] [report output] [table output] [experimental design file] [parameters file]
:: Example :
::Rscript --vanilla run_proteomics_stats.R DA "D:/Utilisateurs/cludwig/stage/RStudio/Tests/Data/Quanti/Parsed_data/Human_yeast/all_mrf_withsamplename.txt" "D:/Utilisateurs/cludwig/stage/RStudio/Tests/test.html" "D:/Utilisateurs/cludwig/stage/RStudio/Tests/test.txt" "D:/Utilisateurs/cludwig/stage/RStudio/Tests/Data/Quanti/Parsed_data/Human_yeast/designExp.txt" "D:/Utilisateurs/cludwig/stage/RStudio/Tests/Data/Quanti/Parsed_data/Human_yeast/Params_da.txt"


:: ##############
:: # ROC module #
:: ##############

:: Command line format : Rscript --vanilla run_proteomics_stats.R DA [input folder] [variant list] [output report] [number of tests] [min ratio] [max ratio] [default ratio] [min pval] [max pval] [default pval]
:: Example :
::Rscript --vanilla run_ROC.R "D:/Utilisateurs/cludwig/stage/RStudio/Tests/Results/Human_yeast/ROC/imputation/gaussian" "D:/Utilisateurs/cludwig/stage/RStudio/Tests/Results/Human_yeast/ROC/variants.csv" "D:/Utilisateurs/cludwig/stage/RStudio/Tests/test.html" 10 0.5 10 0.58 0 0.05 0.05


:: ############
:: # CytoC QC #
:: ############

:: Command line format : Rscript --vanilla run_proteomics_stats.R DA [input] [output folder]
:: Example :
Rscript --vanilla run_cytoc.R "D:/Utilisateurs/cludwig/stage/RStudio/Tests/Data/CytoC/Rapport/cytoc_xics_test.csv" "D:/Utilisateurs/cludwig/stage/RStudio/Tests/"