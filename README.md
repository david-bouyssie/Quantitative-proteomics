# Quantitative-proteomics

This repository contains Rmarkdown modules dedicated to statistical analysis of quantitative proteomics data. <br/>
They are divided up into 2 projects : 
  - Instrumental quality control
  - Quantitative data analysis 
<br/>
The 1rst project contains modules that allow the controlling of LC-MS devices stability through a serie of acquisitions by analyzing the generated signal for a reference protein (ex : CytoC). <br/>
The 2nd project is dedicated to quantitative proteomic data analysis from QC to statistical differential analysis.
<br/>
<br/>

All modules can be launched by command line thanks to the R scripts stored in the **"Launchers"** folder. The modules from the quantitative data analysis project make use of the functions implemented in the library stored in the **"Library"** folder. These modules work with a standardized input format that can be obtain from any quantification file thanks to the parsers stored in the **"Parsers"** folder.
