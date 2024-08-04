# Mpox_Multiplex_Assay
# Description of Data and Analysis Files

The following data and analysis files are contained in the subfolders accompanying the manuscript 
"Differentiation between mpox infection and MVA immunization by a novel machine learning-supported 
serological multiplex assay" by Surtees et al.

## 0 Analyse Coupling
Raw data files in the "input" folder and R script `analyseCoupling_final.R` containing measurements for the reactivity of orthopoxvirus-specific antigens against different mono- and polyclonal but monospecific antibodies (a) or vaccinia immune globulin (VIG) against different coupling batches (b).

The following figures and/or tables are generated:
- Figure S15

## 1 Import Data
Data wrangling of raw data to generate clean data frames for subsequent use in other scripts ("3 Quantify Normalize Data"). Raw data files are in the "input" folder, and the R script `generate_dataInput_final.R` contains functions to import ELISA and Bio-Plex results. Two data frames containing the measured data from a panel of CPXV positive and negative samples, samples of an in-house MVA vaccination panel, and MPXV positive samples, and a negative control panel from MMR diagnostics are generated. `Description.html` contains more information on the data frames.

Input data folder is included as zip file. Unzip to "input" folder to be able to execute to R code

Generated data is contained in two files:
- `dataInput.Rdata`
- `dataInputNK.Rdata`

## 2 Import Metadata
Data wrangling of metadata to generate clean data frames for subsequent use in other scripts. Input metafiles are not made public due to the protection of sensitive data. If needed, data can be obtained for the purpose of peer review only.

Four output files are generated:
- `metadata_IFA.Rdata`: Contains sampleID_metadata, antibody isotype (IgG, IgM), IFA titer, and panel information ("CPXV_MVA") of serum from CPXV diagnostics and in-house MVA study.
- `metadata_MPXV_patients.Rdata`: Full dataset of MPXV metadata (not available).
- `metadata_MPXV_patients_only.Rdata`: CaseID_metadata, date of birth, vaccination_pox (only available for peer review, confidential).
- `metadata_MVA_time.Rdata`: Metadata for MVA study.

## 3 Quantify Normalized Data
Quantification of Bio-Plex raw data (MFI) using 4-parameter fits from VIG dilution series. Contains R script and functions to quantify the IgG and IgM results (also IgA, unpublished data) using the functionality of the `drLumi` package (archived version available here for installation: [drLumi](https://cran.r-project.org/src/contrib/Archive/drLumi/)). For some antigens, automatic fitting on the standard curve measured on the same plate was not possible due to lack of convergence. In this case, standard curves from plates measured on the same day which were highly similar were used for quantification (for details see comments in the R script). Due to the lack of a suitable standard, IgM data was normalized using the IgG standard to account for batch-to-batch variability despite the lack of a suitable standard.

Depends on two input files:
- `../1 Import Data/output/dataInput.Rdata`
- `../1 Import Data/output/dataInputNK.Rdata`

Two output files are generated:
- `dataInputQuant.Rdata` (also exported as an Excel file)
- `dataInputQuantNK.Rdata`

## 4 Method Comparison
Compare results of novel multiplex assay with reference methods (ELISA, IFA, and NT) for a panel of selected sera, which have been quantified by the different methods. Method comparison is performed in the R script `analysisMethodComparison_final.R`, with functions to determine the cut-off values and perform ROC analysis.

Depends on the following input files:
- `../2 Import Metadata/output/metadata_IFA.Rdata`
- `../3 Quantify Normalize Data/output/dataInputQuant.Rdata`
- `../3 Quantify Normalize Data/output/dataInputQuantNK.Rdata`
- `input/corrNT.Rdata` (Generated in the script "9 Correlation NT")
- `input/dataInputNT.R` (Generated as output in the script "9 Correlation NT")

A total of 30 output files is generated, leading to the following figures and tables in the final manuscript:
- Figure 1: Comparison between results obtained by the multiplex assay and ELISA, IFA, and NT as reference methods.
- Figure S9: ROC analysis for comparison between the IFA results (titer ≤ 1:80 negative) and the IgG (a) or IgM (b) binding response determined in the multiplex assay against the tested antigens.
- Figure S10: Correlation between binding to the recombinant protein antigens in the multiplex assay and results obtained by an in-house ELISA for IgG and IgM detection.
- Figure S11: Comparison between the antigens implemented in the multiplex assay with the IFA titers (a) or the NT titers (b) for IgM detection.
- Figure S12: Plot of Pearson’s or Spearman correlation coefficients between the antigens implemented in the multiplex assay and the ELISA for IgG (a) or IgM (b) detection (Pearson’s r), the IFA for IgG (c) and IgM (d), or the NT for IgG (e) or IgM (f) detection (Spearman rho rank correlation). Combined using Adobe Illustrator.
- Table S9: Cut-off values and 95% confidence intervals (CI) as determined by method comparison between IFA and the multiplex assay for IgG antibodies.
- Table S10: Cut-off values and 95% confidence intervals (CI) as determined by method comparison between IFA and the multiplex assay for IgM antibodies.
- Table S11: Performance parameters (median, lower and upper 95% confidence intervals) for method comparison between IFA and the multiplex assay for IgG antibodies.
- Table S12: Performance parameters (median, lower and upper 95% confidence intervals) for method comparison between IFA and the multiplex assay for IgM antibodies.
- Table S13: Pearson’s correlation coefficients between ELISA and multiplex results for the tested antigens and antibody isotypes.
- Table S14: Spearman rho rank correlation coefficients between IFA titers and multiplex results for the tested antigens and antibody isotypes.
- Table S15: Spearman rho rank correlation coefficients and Pearson’s r correlation coefficients between NT titers and multiplex results for the tested antigens and antibody isotypes.

## 7 Patient Panel Merge
Merges metadata for mpox panel and MVA panel with quantified data. The R file `mergeMetadata_final.R` is used. The sera count towards the "acute" panel of the manuscript.

Depends on the following input files:
- `../4 Method Comparison/output/dataInputQuantCat.Rdata`
- `../2 Import Metadata/output/metadata_MVA_time.Rdata`
- `../2 Import Metadata/output/metadata_MPXV_patients.Rdata`

The following output files containing the merged data files are generated:
- `dataInputMPXVmeta`, file = `output/dataInputMPXVmeta_all.Rdata`
- `dataInputMVAmeta`, file = `output/dataInputMVAmeta.Rdata`

## 9 Correlation NT
Method comparison between results from virus NT and multiplex assay. The R file `AnalyseCorrelation_final.R` contains the data and analysis for the comparison between the NT test and the multiplex assay.

Depends on the following input files:
- `../7 Patient Panel Merge/output/dataInputMPXVmeta.Rdata`
- `input/Ergebnisse_NT MPXV_für ZBS3.xlsx` (Not contained due to data safety reasons)
- `input/2022_10_21_Imvanex_Studie_anonym_Titer für Daniel.xlsx` (Not contained due to data safety reasons)

The following output files are generated:
- `../4 Method Comparison/input/dataInputNT.R` (Input files)
- `../4 Method Comparison/input/corrNT.Rdata` (Results of correlation analysis as data frames and plots)
- `output/sensNT.xlsx`
- `output/specNT.xlsx`
- `output/accNT.xlsx`
- `output/plotROCNT.png`

The following figures and/or tables are generated:
- Parts of Figure S1
- Parts of Figure S11
- Parts of Figure S12

## 11 Shiny App Spox
A Shiny app to analyze results of measurements for the epi-panel. It analyzes results of two plates simultaneously. Metadata is loaded from Word protocol files, MFI data are quantified using VIG as a standard (IgG data for both IgG and IgM). Predictions based on an early LDA model trained on acute data are also contained but not used for the final predictions. The folder contains input and output from one exemplarly run. Full 
input and output files can be provided, if needed. 

Depends on the following input files:
- `input_final/` (All files in the folder)
- `input_quant/` (All files in the folder)
- `input_nc/` (All files in the folder)
- `input_nc_quant/` (All files in the folder)
- `input_rep/` (All files in the folder)
- `input_rep_quant/` (All files in the folder)
- `MPOX_clinic2.xlsx` (Meta data mpox)
- `MMR Kontrollseren MPox Panel.xlsx` (Meta data NC panel)

Generates the following output files:
- `data/dataInSpox.R`
- `data/dataInQuantSpox.R`
- `data/dataInQuantNCCombined.R`
- `data/dataInSpoxRep.R`
- `data/dataInQuantSpoxRep.R`

## 14 Vergleich Spox
Merge results for quantification generated by the Shiny app. The R script `analyseSpox_final.R` contains the

 code to generate data frames for the results of the quantification and predictions of an early LDA model (legacy code) for the epidemiological panel measured (Spox), repetition measurements (rep), and a second negative panel (NC) which was included, as the first panel did not contain information about the age of the subjects.

## 15 Generate Dataframe ML
Combines data from different panels and metadata to generate a final CSV file as input for ML training and testing of different models. The R script `generateMLInput_final.R` combines data frames and metadata from different panels in one large CSV file for subsequent input for ML training and testing.

Depends on the following input files:
- `../14 Vergleich Spox/data/dataInQuantSpox.R` (dataInQuant -> Spox Dataframe IgG)
- `../14 Vergleich Spox/data/dataInSpox.R` (dataIn -> Spox Dataframe IgG for analysis)
- `../14 Vergleich Spox/data/dataInQuantNCCombined.R` (dataInQuantNCCombined -> New NC panel)
- `../14 Vergleich Spox/data/dataInNCCombined.R` (dataInNCCombined -> IgG analyzed and classified panel)
- `../14 Vergleich Spox/data/dataInQuantSpoxRep.R` (dataInQuantRep -> Repeated measurement of select samples)
- `../14 Vergleich Spox/data/dataInSpoxRep.R` (dataInRep -> IgG analyzed repeated measurements)
- `../4 Method Comparison/output/dataInputQuantCat.Rdata`
- `../4 Method Comparison/output/dataInputQuantCatNK.Rdata`
- `../4 Method Comparison/output/threshold.Rdata`
- `../2 Import Metadata/output/metadata_MVA_time.Rdata`
- `../2 Import Metadata/output/metadata_MPXV_patients.Rdata`
- `../7 Patient Panel Merge/output/dataInputMPXVmeta_all.Rdata` (dataInputMPXVmeta -> New: including all data)
- `../7 Patient Panel Merge/output/dataInputMVAmeta.Rdata` (dataInputMVAmeta)

Generates the following output files:
- `output/dataInputAll.csv` (All data)
- `output/dataInputComplete.csv` (Only data with complete metadata files)
- `output/dataPredictionSPox.csv` (Results only for the seroepidemiological panel)

## 16 Compare Panels
Combines plots from Three-Way ANOVA and analysis of stratified panels to generate figures and supplementary figures. The R script `CompareSeroconversion_final.R` generated the plot for an overview of the different panels.

Depends on the following input files:
- `../24 Spox Pre Pos Clustering/input/dataClustering.Rdata`
- `../24 Spox Pre Pos Clustering/output/plotCombinedAll.Rdata`
- `../25 Three-Way ANOVA Panel/output/boxplotsAnova.Rdata`

Generates the following output files:
- `output/dataInputPanelANOVA.Rdata` (Data frame with harmonized metadata for Three-Way ANOVA)
- `output/plotFig3New.png`
- `output/plotFig3New.pdf`
- `output/Sup_Fig_Ratios.png`

The following figures and/or tables are generated:
- Figure S14
- Figure 3 (Final figure labels corrected in Adobe Illustrator)

## 17 Compare Duplicate Measurements
Analyze results between quantified results for independent replications. The R file `CompareDuplicates.R` reads results for original and replicate measurements for a panel of n = 40 duplicate measurements.

## 23 Spox Mpox Pos
Combine and harmonize age group data for determination of unrecognized mpox infections in unvaccinated subjects in the epidemiological panel. The R file `AnalyseSpoxPos_final.R` contains code for the harmonization of the age groups, which have been coded differently in the different panels, to generate one harmonized data frame.

Depends on the following input files:
- `input/SPox_Proben.xlsx` (Not contained due to data protection)
- `input/predictionSpoxComplete.xlsx`
- `input/dataInputAll.csv`
- `input/MMR Kontrollseren MPox Panel.xlsx`

Generates the following output file:
- `output/dataClustering.Rdata`

## 24 Spox Pre Pos Clustering
Identify sera within the epidemiological panel that are potentially misclassified due to reactivity against different OPXV-specific antigens but unlikely to have received childhood vaccination based on age. The R file `ClusterPrePos_final.R` contains the code to generate population-based cut-off values to discern samples with unusually high reactivity against OPXV-specific antigens in age groups with low likelihood of having received smallpox vaccinations in early childhood.

Depends on the following input files:
- `input/dataClustering.Rdata`

Generates 13 output files, among those a list with sample IDs to exclude from ML training and testing based on the population-based cut-off values.

The following figures and/or tables are generated:
- Figure S16
- Figure S17

## 25 Three-Way ANOVA Panel
Performs Three-Way ANOVA to compare differences between the different panels. The R script `Three-Way-Anova_final.R` contains the code to perform the Three-Way ANOVA on the different panels tested.

Depends on the following input files:
- `../16 Compare Panels/output/dataInputPanelANOVA.Rdata`

Generates the following output files:
- `output/boxplotsAnova.Rdata` (Combined boxplot for generation of Figure 3) -> Not contained: file too big for upload: will be posted to figshare
- `output/tables.Rdata` (Generation of supporting tables) -> Not contained: file too big for upload: will be posted to figshare

The following figures and/or tables are generated:
- Figure 3
- Figure S13
- Figure S22
- Figure S23
- Table S16
- Table S17
- Table S18
- Table S19
- Table S20
- Table S21
- Table S22
- Table S23
- Table S24
- Table S25
- Table S26
- Table S27
