@echo off
set PATH=%PATH%;C:\Program Files\R\R-3.5.1\bin


::###############################################################################################################################################################
:: ###########
:: # Parsing #
:: ###########



:: ### CARINE ###
::Rscript --vanilla run_parsing.R ../Data/Quanti/Raw/Proline/Carine/proteinSet.txt ../Data/Quanti/Parsed_data/Carine/Proline/20180320_I2MC_Prats_QC20190404.txt Proline abundance

:: ### PASCALINE ###
::Rscript --vanilla run_parsing.R ../Data/Quanti/Raw/Proline/Pascaline/PBB_190418_BCG_IP_n1_MabA_Glu-C+Tryp_19-04-2019_Proline.xlsx ../Data/Quanti/Parsed_data/Pascaline/PBB_190418_BCG_IP_n1_MabA_Glu-C+Tryp_19-04-2019_Proline_abundance.txt Proline raw_abundance 4

:: ### HUMAN_YEAST ###
::Rscript --vanilla run_parsing.R ../Data/Quanti/Raw/Proline/Human_yeast/all_mrf_withsamplename.xlsx ../Data/Quanti/Parsed_data/Human_yeast/all_mrf_withsamplename.txt Proline abundance 4

:: ### DAVID ###
::Rscript --vanilla run_parsing.R "../Data/Search Xtandem PME12modified Percolator validation_MRF_uniquepep_discardMC-pepnorm_16-05-2019_0941.xlsx" "../Data/Quanti/Parsed_data/David/Search Xtandem PME12modified Percolator validation_MRF_uniquepep_discardMC-pepnorm_16-05-2019_0941.txt" Proline abundance 4



::###############################################################################################################################################################
:: ######
:: # QC #
:: ######


:: ### PASCALINE ###
::Rscript --vanilla run_proteomics_stats.R QC ../Data/Quanti/Parsed_data/Pascaline/PBB_190418_BCG_IP_n1_MabA_Glu-C+Tryp_19-04-2019_Proline_abundance.txt ../Results/Pascaline/QC/test.html TRUE FALSE

:: ### HUMAN_YEAST ###
::Rscript --vanilla run_proteomics_stats.R QC ../Data/Quanti/Parsed_data/Human_yeast/all_mrf_withsamplename.txt ../Results/Human_yeast/QC/all_mrf_withsamplename_new.html TRUE FALSE



::###############################################################################################################################################################
:: ######
:: # DA #
:: ######



::## CARINE ##
::Rscript --vanilla run_proteomics_stats.R DA ../Data/Quanti/Parsed_data/Carine/Proline/20180320_I2MC_Prats_QC20190404(2).txt ../Results/Carine/DA/test.html ../Results/Carine/DA/test.txt ../Data/Quanti/Parsed_data/Carine/Proline/designExp.txt ../Data/Quanti/Parsed_data/Carine/Proline/Params_da.txt

::## CAROLE ##
::Rscript --vanilla run_proteomics_stats.R DA ../Data/Quanti/Parsed_data/Carole/190409_Quanti_4_conditions_3R_2r_avec_MV_11-04-2019_1523.txt ../Results/Carole/DA/190409_Quanti_4_conditions_3R_2r_avec_MV_11-04-2019_1523.html ../Results/Carole/DA/190409_Quanti_4_conditions_3R_2r_avec_MV_11-04-2019_1523.txt ../Data/Quanti/Parsed_data/Carole/designExp.txt ../Data/Quanti/Parsed_data/Carole/Params_da.txt

::## ALEX ##
::Rscript --vanilla run_proteomics_stats.R DA ../Data/Quanti/Parsed_data/Alex/Quantification_Manip_DCIR_EGTA_Injection_1_+_2_061218_10-12-2018_1541.txt ../Results/Alex/DA/test_LOG_Student.html ../Results/Alex/DA/test_percentileCol_LOG_Student.txt ../Data/Quanti/Parsed_data/Alex/designExp.txt ../Data/Quanti/Parsed_data/Alex/Params_da.txt

::## DAVID ##
::Rscript --vanilla run_proteomics_stats.R DA "../Data/Quanti/Parsed_data/David/Search Mascot PME12 modified import validation Mascot quantification Proline_Medianbilogicalprofile_discardallMC_pepUnique_normpep_16-05-2019_0929.txt" "../Results/David/DA/Search Mascot PME12 modified import validation Mascot quantification Proline_Medianbilogicalprofile_discardallMC_pepUnique_normpep_16-05-2019_0929.html" "../Results/David/DA/Search Mascot PME12 modified import validation Mascot quantification Proline_Medianbilogicalprofile_discardallMC_pepUnique_normpep_16-05-2019_0929.txt" ../Data/Quanti/Parsed_data/David/designExp.txt ../Data/Quanti/Parsed_data/David/Params_da.txt

::Rscript --vanilla run_proteomics_stats.R DA "../Data/Quanti/Parsed_data/David/Search Mascot PME12 modified import validation quantification Proline_no post process.txt" "../Results/David/DA/Search Mascot PME12 modified import validation quantification Proline_no post process.html" "../Results/David/DA/Search Mascot PME12 modified import validation quantification Proline_no post process.txt" ../Data/Quanti/Parsed_data/David/designExp.txt ../Data/Quanti/Parsed_data/David/Params_da.txt

::Rscript --vanilla run_proteomics_stats.R DA "../Data/Quanti/Parsed_data/David/Search Xtandem PME12modified Percolator validation_Medianbilogicalprofile_discardallMC_pepUnique_normpep.txt" "../Results/David/DA/Search Xtandem PME12modified Percolator validation_Medianbilogicalprofile_discardallMC_pepUnique_normpep.html" "../Results/David/DA/Search Xtandem PME12modified Percolator validation_Medianbilogicalprofile_discardallMC_pepUnique_normpep.txt" ../Data/Quanti/Parsed_data/David/designExp.txt ../Data/Quanti/Parsed_data/David/Params_da.txt

::Rscript --vanilla run_proteomics_stats.R DA "../Data/Quanti/Parsed_data/David/Search Xtandem PME12modified Percolator validation_MRF_uniquepep_discardMC-pepnorm_16-05-2019_0941.txt" "../Results/David/DA/Search Xtandem PME12modified Percolator validation_MRF_uniquepep_discardMC-pepnorm_16-05-2019_0941.html" "../Results/David/DA/Search Xtandem PME12modified Percolator validation_MRF_uniquepep_discardMC-pepnorm_16-05-2019_0941.txt" ../Data/Quanti/Parsed_data/David/designExp.txt ../Data/Quanti/Parsed_data/David/Params_da.txt

::## HUMAN_YEAST ##
::Rscript --vanilla run_proteomics_stats.R DA ../Data/Quanti/Parsed_data/Human_yeast/all_mrf_withsamplename.txt ../Results/Human_yeast/DA/111.html ../Results/Human_yeast/DA/111.txt ../Data/Quanti/Parsed_data/Human_yeast/designExp.txt ../Data/Quanti/Parsed_data/Human_yeast/Params_da.txt


::###############################################################################################################################################################
:: #######
:: # ROC #
:: #######


Rscript --vanilla run_ROC.R ../Results/Human_yeast/ROC/imputation/gaussian ../Results/Human_yeast/ROC/variants.csv ../Results/Human_yeast/ROC/imputation/gaussian/comparison.html 10 0.5 10 0.58 0 0.05 0.05


::###############################################################################################################################################################
:: #########
:: # CytoC #
:: #########

::Rscript --vanilla run_cytoc.R ../../Data/CytoC/Rapport/cytoc_xics_test.csv ../Results/CytoC
::Rscript --vanilla run_cytoc.R ../../Data/CytoC/Rapport/cytoc_xics_ref.csv ../Results/CytoC/
