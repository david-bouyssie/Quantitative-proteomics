
:: Command lines of the modules, run with the portable R after putting it and its pandoc on the PATH:
::   set "PATH=<path to>\R-Portable_with_QP\bin;%PATH%"
:: The launchers find the repository files from their own location, so they run from any folder;
:: input and output paths are relative to the current folder. The examples below are run from R-Portable_with_QP.

:: ###################
:: # Parsing modules #
:: ###################

:: Command line format : Rscript --vanilla Quantitative-proteomics\Launchers\run_parsing.R [file to parse] [output] [parser option : Proline/MaxQuant/ProteomeDiscoverer] [abundance type option : LFQ/iBAQ/raw_abundance/abundance] [if xlsx file, sheet number]
:: Example (abundance intensities of sheet 4, "Protein sets", of a Proline export) :
Rscript --vanilla Quantitative-proteomics\Launchers\run_parsing.R "<Proline export>.xlsx" parsed_example.tsv Proline abundance 4


:: ###############################
:: # Quantitative data QC module #
:: ###############################

:: Command line format : Rscript --vanilla Quantitative-proteomics\Launchers\run_proteomics_stats.R QC [input (standardized format)] [output report] [parameters] [experimental design (optional)] [original quantification file (optional)]
:: Example :
Rscript --vanilla Quantitative-proteomics\Launchers\run_proteomics_stats.R QC Quantitative-proteomics\Example\Datasets\QC-DA\Parsed_proteins_set.txt QC_example.html Quantitative-proteomics\Example\Datasets\QC-DA\QC_parameters.txt Quantitative-proteomics\Example\Datasets\QC-DA\exp_design.txt

:: ################################
:: # Differential analysis module #
:: ################################

:: Command line format : Rscript --vanilla Quantitative-proteomics\Launchers\run_proteomics_stats.R DA [input (standardized format)] [output folder, created if needed] [experimental design file] [parameters file] [original quantification file (optional)]
:: DA.html and DA_output_table.xlsx are written in the output folder
:: Example :
Rscript --vanilla Quantitative-proteomics\Launchers\run_proteomics_stats.R DA Quantitative-proteomics\Example\Datasets\QC-DA\Parsed_proteins_set.txt DA_example Quantitative-proteomics\Example\Datasets\QC-DA\exp_design.txt Quantitative-proteomics\Example\Datasets\QC-DA\Parameters.txt

:: ##############
:: # ROC module #
:: ##############

:: Command line format : Rscript --vanilla Quantitative-proteomics\Launchers\run_ROC.R [input folder] [variant list] [output report] [number of tests] [min ratio] [max ratio] [default ratio] [min pval] [max pval] [default pval]
:: Example :
Rscript --vanilla Quantitative-proteomics\Launchers\run_ROC.R Quantitative-proteomics\Example\Datasets\ROC Quantitative-proteomics\Example\Datasets\ROC\variants.csv ROC_example.html 10 0.5 10 0.58 0 0.05 0.05


:: ############
:: # CytoC QC #
:: ############

:: Command line format : Rscript --vanilla Quantitative-proteomics\Launchers\run_cytoc.R [cytoc_xics file] [output report]
:: Example :
Rscript --vanilla Quantitative-proteomics\Launchers\run_cytoc.R Quantitative-proteomics\Example\Datasets\CytoC\cytoc_xics_test.csv CytoC_example.html
