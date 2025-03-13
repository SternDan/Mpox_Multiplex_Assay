# Mpox_Multiplex_Assay
# Description of Data and Analysis Files

The following data and analysis files are contained in the subfolders accompanying the manuscript 
"Differentiation between mpox infection and MVA immunization by a novel machine learning-supported 
serological multiplex assay" by Surtees et al.

## Installation and usage
Install R (version 4.3.0) and R Studio (2023-06-2). Install all necessary packages, which are listed in each R script before running the scripts. All necessary data files (input and output) are contained within the different folders.

## Figure_01_Method_Comparison
Compare results of novel multiplex assay with reference methods (ELISA, IFA, and NT) for a panel of selected sera, which have been quantified by the different methods. Method comparison is performed in the R script `Fig_1_Method_Comparison.R`.

Depends on the following input files:
- `input/dataInputELISAMultiplex.Rdata`
- `input/dataInputNT.Rdata`
- `input/dataInputQuant.Rdata`
- `input/metadata_IFA.Rdata`

The following output files are generated:
- Figure 1: Comparison between results obtained by the multiplex assay and ELISA, IFA, and NT as reference methods. Correlation plots and plot for comparison with Delta-VACV antigen are combined in Adobe Illustrator
- Figure S1: Plot of multiplex antigens in comparison to ELISA as reference method
- parampaBa.xlsx: Parameters for Passing-Bablock-Regression save as Excel-file

## Figure_02_Compare_Panels
Plot IgG and IgM results of the differen serogroups (pre-immune, MVA, MPXV) and panels (acute, epi) stratified by childhood immunisation against smallpox. Analysis and plotting is performed in the R script `Figure_2_Compare_Panels.R`.

Depends on the following input files:
- `input/dataInputComparePanels.Rdata`

The following output files are generated:
- Figure_2: Spider plots and plots of selected antigens from the acute and epi panel. Plot of ratios of selected homologue antigen pairs.
- Figure S2: Plot of antigens not contained in figure 2
- Figure S3: Plot of ratios for antigen pairs not contained in figure 2

## Figure_04_ML_Performance

