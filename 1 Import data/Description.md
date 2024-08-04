# Data import for Mpox-specific multiplex-assay

The R script *generate_dataInput.R* takes Excel-files as data input and generates a large dataframe containing all experimental data and assay-specific metadata to analyse the establishment and validation of a Mpox-specific multiplex-assay.
The final dataframe is called *dataInput* and is exported as *dataInput.Rdata*. It contains the following variables: 

1. **experiment**: *Character vector* describing the overall experiment of the data file
2. **filename**: *Character vector* containing the file name of the experimental data file
3. **date**: *Datetime vector* containing the date when the experiment was performed
4. **plate**: *Numeric vector* containing number of the plate that was measured in the experiment
5. **well**: *Character vector* naming the well position on each plate
6. **assaytype**: *Character vector* describing which assay was performed (ELISA or Multiplex)
7. **sample_type**: *Character vector* with information on type of date (unknown, standard, negative)
8. **sampleID_assay**: *Character vector* containing the sample ID used in the assay
9. **sampleID_metadata**: *Character vector* containing the sample ID of the metadata
10. **dilution**: *Numeric vector* with information on the dilution
11. **isotype**: *Character vector* with information on the antibody isotype
12. **antigen**: *Character vector* with information on the antigen
13. **data**: *Numeric vector* with data (either OD or MFI)
14. **panel**: *Character vector* with information on the serum panel

The following serum panels are included in this data set:
- "CPXV_MVA": sera positive for CPXV (as determined by diagnostic PCR) or sampled within an in house MVA study
- "MPXV": sera positive for MPXV (as determined by diagnostic PCR) and sampled during the
peak of the mpox outbreak in Berlin in 2022

The second datafile termed "dataInputNK" contains sera from an mpox-negative panel, sampled 
before the outbreak in 2019 without known history of mpox infection.
