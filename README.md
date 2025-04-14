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
Reads assay performance of the comparison of different ML algorithms tested on different panels. Analysis and plotting is performed in the R script `Figure_4.R`.

Depends on the following input files:
- `input/dataIn.Rdata`
- `input/dataInMeta.Rdata`
- `input/dataInputAnalyte.Rdata`
- `input/dataInRep.Rdata`
- `input/statisticalDataCombined.Rdata`

The following output files are generated:
- Fig4ab.pdf: Plot of F1 scores for the comparison of the different algorithms and panels. Plot of specificity of GBC and LDA on the negative control panel
- Fig4_all_xgboost_lda.pdf: Circular plot of misclassifications on all data set, stratified by childhood smallpox vaccionation
- freqTableGBCLDA.xlsx: Frequency table underlying the circula plot in figure 4
- freqTableAll.Rdata: Frequence table for misclassificaton for other algorithms
- supTablePerformanceRev.Rdata: Data underlying Table S9 (Generated using SupportingTablesPerformanceRev.Rmd -> Generates Docx file)

## Figure_05_Validation
Analysis and generation of figures and table for the validation measurements. Analysis and plotting is performed in the R script `Figure_5_Validation.R`.

Depends on the following input files:
- `input/dataInputComparison.Rdata`
- `input/ensemblePrediction.Rdata`
- `input/heatmap_input.Rdata`

The following output files are generated:
- heatmap_IgG_white.pdf: Heatmap of validation panel together with plots of metadata for IgG data
- heatmap_IgM_white.pdf: Heatmap of validataion panel together with plots of metadata for IgM data
- plotConfusionVal.pdf: Confusion matrices for validation measurements depending on all data, curated data, and data filtered by ensemble prediction confidence
- plotCorrelationAllIgG.pdf: Correlation plot between IgG results, ensemble prediction confidence for pre-immune, MVA, or MPXV class
- plotMPXVMVAPre.pdf: Plot of classification for ensemble prediction
- STableEnsemblePrediction.xlsx: Table S10 with all data for the validaton panel

## Figure_S05_Bead_coupling
Generates plot for coupling control and batch-to-batch variabiliy using in the R script analyseCoupling.R.

Depends on the following input files:
- `input/dataInputBatch.Rdata`
- `input/dataInputPlotting.Rdata`

The following output files are generated:
- plotCoupling.png: Figure S4 containing the plots for the coupling control and the batch-to-batch variability

## Figure_S06_07_Exclusion_Training
Generates plots and population-based cut-off values used to exclude samples with hints on unrecognized infections from the epi panel from ML based training and testing. Not to confuse with more robust cut-off values, which have been generated to classify sera based on binary classifiers as described in Figure_S07_ROC.

Depends on the following input files: 
- `input/dataClustering_post.Rdata`
- `input/dataInWide.Rdata`

The following output files are generated (amoung others):
- SFig_plotCombinedyoung.png: Figure S5 Comparison between antigens and panels in young population
- SFig_plotROC.png: Figure S6 ROC curve

## Figure_S08_Classical_ROC
Generate plot for the performance of single antigens based on ROC analysis as well as threshold values. 

Depends on the following input files: 
- `input/dataInWide.Rdata`
- `input/plotImprovePerformance.Rdata`

The following output files are generated:
- plotFig6.png: Figure S8
- ROC_parameters_classic.xlsx: Table S11 with performance parameters for ROC analysis for single antigens
