@echo off
set PATH=%PATH%;C:\Program Files\R\R-3.5.1\bin


:: ###################
:: # Parsing modules #
:: ###################

:: The command line expected : Rscript --vanilla run_parsing.R [file to parse] [output] [parser option : Proline/Maxquant/ProteomeDiscoverer] [abundance type option : LFQ/iBAQ/raw_abundance/abundance] [if xlsx file, sheet number]
:: Example : 
Rscript --vanilla run_parsing.R "../Data/Search Xtandem PME12modified Percolator validation_MRF_uniquepep_discardMC-pepnorm_16-05-2019_0941.xlsx" "../Data/Quanti/Parsed_data/Search Xtandem PME12modified Percolator validation_MRF_uniquepep_discardMC-pepnorm_16-05-2019_0941.txt" Proline abundance 4


:: ###############################
:: # quantitative data QC module #
:: ###############################

:: The command line expected : Rscript --vanilla run_proteomics_stats.R QC [input] [output] [normalization : TRUE/FALSE] [delete empty lines : TRUE/FALSE]
:: Example :
Rscript --vanilla run_proteomics_stats.R QC ../Data/Quanti/Parsed_data/PBB_190418_BCG_IP_n1_MabA_Glu-C+Tryp_19-04-2019_Proline_abundance.txt ../Results/Pascaline/QC/test.html TRUE FALSE


:: ################################
:: # Differential analysis module #
:: ################################

:: The command line expected : Rscript --vanilla run_proteomics_stats.R DA [input] [report output] [data output] [experimental design file] [parameters file]
:: Example :
Rscript --vanilla run_proteomics_stats.R DA ../Data/Quanti/Parsed_data/Human_yeast/all_mrf_withsamplename.txt ../Results/Human_yeast/DA/111.html ../Results/Human_yeast/DA/111.txt ../Data/Quanti/Parsed_data/Human_yeast/designExp.txt ../Data/Quanti/Parsed_data/Human_yeast/Params_da.txt


:: #######
:: # ROC #
:: #######

:: The command line expected : Rscript --vanilla run_proteomics_stats.R DA [input folder] [variant list] [output] [number of tests] [min ratio] [max ratio] [chosen ratio] [min pval] [max pval] [chosen pval]
:: Example :
Rscript --vanilla run_ROC.R ../Results/Human_yeast/ROC/imputation/gaussian ../Results/Human_yeast/ROC/variants.csv ../Results/Human_yeast/ROC/imputation/gaussian/comparison.html 10 0.5 10 0.58 0 0.05 0.05


:: ############
:: # CytoC QC #
:: ############

:: The command line expected : Rscript --vanilla run_proteomics_stats.R DA [input] [output folder]:: Example :
Rscript --vanilla run_cytoc.R ../../Data/CytoC/Rapport/cytoc_xics_ref.csv ../Results/CytoC/
