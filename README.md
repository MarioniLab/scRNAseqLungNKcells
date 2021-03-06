# scRNAseqLungNKcells

This repo contains Rmarkdown files to analyze 10x Genomics single-cell RNAseq data in Schuijs et al. 2020 comparing NK cells sorted from the lungs of mice that had been treated with PBS or IL-33. 

## Preparing data and directories

These scripts use the output of CellRanger (v.3.1.0) and assume that the raw (unfiltered) CellRanger output is in a directory called "raw_data" within the repository.

To run the scripts, also create directories called "data" and "tables_plots" within the repository.

## Scripts

Generally useful functions are contained in 
NKcell_plotting_functions.R

Scripts were run in the following order on the unfiltered CellRanger output.

Pre-processing:
NKcellIL33Preprocess.Rmd

QC and dimensionality reduction:
NKcellIL33QCandDimRed.Rmd

Removal of a few aberrant macrophages:
NKcellIL33rmMP.Rmd

Differential abundance and expression analyses:
NKcellIL33AnalysisDADE.Rmd

Plotting results of various analyses for figures:
NKcellIL33Plots.Rmd
