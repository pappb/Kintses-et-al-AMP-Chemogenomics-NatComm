# Kintses-et-al-AMP-Chemogenomics-NatComm
Code deposited for paper "Chemical-genetic profiling reveals limited cross-resistance between antimicrobial peptides with different modes of action", Kintses et al. Nature Communications


The programs are written in R.

 __script folder__
* _script/cg01_RemoveBackground.R_  
* _script/cg02_rescaleData.R_
* _script/cg03_plotsOfDataAfterBgRemoveAndRescale.R_ Generates figures to visualise data grouped by treatments. It is not needed for  further analisys.
* _script/cg04_plotRawBarplotsForEachLine.R_ Generates figures to visualise data grouped by genes. It is not needed for further analisys.
* _script/cg05_normalizationOfBimodalDist.R_ Normalises data
* _script/functions-for-chemogenomic-project.R_ contains some helper functions.
* _script/pdfTools.R_ contains some more helper functions for the normalisation script.


__data folder__
 * It contains the input files of the scripts 
 * _data/dataTable1-expression_values_of_sequencing.csv_ - this file contains the read numbers from the sequencing
 * _data/batch-structure.csv_ - this table describes which data columns come from the same assay
 * all the intermediate data tables will be placed here.
 
