addQC_function <- function(inputData){
  # Add flags for QC samples
  dataInputQCHSA <-
    inputData %>%
    filter(analyte == "HSA", dilution == 100) %>% 
    group_by(experiment, plate_assay, sampleID_assay, isotype, dilution) %>% 
    dplyr::summarise(highBG = if_else(data > 500, T, F),
              background = data)
  
  dataInputQCIgG <-
    inputData %>%
    filter(analyte == "ahIgG", dilution == 100, 
           isotype == "IgG") %>% 
    group_by(experiment, plate_assay, sampleID_assay, isotype, dilution) %>% 
    dplyr::summarise(lowIg = if_else(data < 10000, T, F),
              positiveIg = data)
  
  dataInputQCIgM <-
    inputData %>%
    filter(analyte == "ahIgM", dilution == 100,
           isotype == "IgM") %>% 
    group_by(experiment, plate_assay, sampleID_assay, isotype, dilution) %>% 
    dplyr::summarise(lowIg = if_else(data < 6000, T, F),
              positiveIg = data)
  
  dataInputQC <-
    rbind(dataInputQCIgG, dataInputQCIgM)
  
  # Normalize data
  dataInput <-
    inputData %>%
    dplyr::left_join(dataInputQC, by = c("experiment", "plate_assay","sampleID_assay", "isotype",
                                  "dilution")) %>%
    dplyr::left_join(dataInputQCHSA, by = c("experiment", "plate_assay","sampleID_assay", "isotype",
                                     "dilution"))
  return(dataInput)
}