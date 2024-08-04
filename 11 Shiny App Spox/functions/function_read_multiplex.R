## Define function for data input and merge with experimental metadata ----
read_plate_function <-
  function(filename1, plate_input){
    mfidata_plate_1 <- 
      readxl::read_xlsx(filename1, skip = 7) %>%
      filter(`...2` != "Well") %>%
      dplyr::rename(type = `...1`, well = `...2`) %>%
      mutate_at(vars(-well, -type), ~(as.numeric(sub(",", ".", ., fixed = TRUE)))) %>%
      pivot_longer(c(-well, -type), names_to = "analyte", values_to = "data") %>%
      dplyr::mutate(analyte = str_replace(analyte, " \\([^\\)]+\\)", ""),
             analyte = str_remove(analyte, " Lysat")) %>%
      dplyr::select(-type) %>%
      unique() %>%
      pivot_wider(names_from = analyte,
                  values_from = data) %>%
      dplyr::mutate(Delta = VACV - Hep2,
             Delta = Delta + abs(min(Delta))) %>%
      pivot_longer(cols = - well, names_to = "analyte", values_to = "data") %>% 
      dplyr::mutate(plate_assay = plate_input) %>% 
      filter(!is.na(data))
    return(mfidata_plate_1)
  }


read_metadata_function <- function(metafilename, dataFile1, dataFile2) {
  # Read metadata contained in metadata file
  #read_metadata <- function(dataIn){
  dataIn <- metafilename
  
  # Extract used bead regions
  beadRegions <- docx_extract_tbl(docxtractr::read_docx(dataIn),
                                  tbl_number = 1) %>% 
    dplyr::rename(region = Bead.region, antigen = starts_with("Bezeichnung"))
  
  # Extract plate layout plate 1
  plateLayout_plate_1 <- docx_extract_tbl(docxtractr::read_docx(dataIn),
                                          tbl_number = 2,
                                          preserve = TRUE) %>% 
    pivot_longer(c(-`X.`), names_to = "col", values_to = "sample") %>%
    dplyr::mutate(col = str_remove(col, "X"),
           well = paste(`X.`, col, sep = "")) %>%
    dplyr::select(well, sampleID_assay = sample) %>%
    separate(sampleID_assay, into = c("sampleID_assay", "dilution"), sep = "\n") %>%
    dplyr::mutate(dilution = if_else(is.na(dilution), parse_number(sampleID_assay),
                              as.numeric(dilution)),
           sampleID_assay = if_else(is.na(as.numeric(sampleID_assay)),
                                    gsub("[[:digit:]]", "", sampleID_assay),
                                    sampleID_assay),
           sampleID_assay = case_when(sampleID_assay == "Blank" ~ "Background0",
                                      sampleID_assay == "VIG " & dilution == 500 ~ "Standard1",
                                      sampleID_assay == "VIG " & dilution == 2000 ~ "Standard2",
                                      sampleID_assay == "VIG " & dilution == 8000 ~ "Standard3",
                                      sampleID_assay == "VIG " & dilution == 32000 ~ "Standard4",
                                      sampleID_assay == "VIG " & dilution == 128000 ~ "Standard5",
                                      TRUE ~ sampleID_assay),
           sampleID_assay = str_remove_all(sampleID_assay, " "),
           sample_type = case_when(sampleID_assay == "NK" ~ "negative",
                                   grepl("Background", sampleID_assay) ~ "standard",
                                   grepl("Standard", sampleID_assay) ~ "standard",
                                   TRUE ~ "unknown"),
           plate_assay = 1)
  
  # Extract plate layout plate 2
  plateLayout_plate_2 <- docx_extract_tbl(docxtractr::read_docx(dataIn),
                                          tbl_number = 3,
                                          preserve = TRUE) %>% 
    pivot_longer(c(-`X.`), names_to = "col", values_to = "sample") %>%
    dplyr::mutate(col = str_remove(col, "X"),
           well = paste(`X.`, col, sep = "")) %>%
    dplyr::select(well, sampleID_assay = sample) %>%
    separate(sampleID_assay, into = c("sampleID_assay", "dilution"), sep = "\n") %>%
    dplyr::mutate(dilution = if_else(is.na(dilution), parse_number(sampleID_assay),
                              as.numeric(dilution)),
           sampleID_assay = if_else(is.na(as.numeric(sampleID_assay)),
                                    gsub("[[:digit:]]", "", sampleID_assay),
                                    sampleID_assay),
           sampleID_assay = case_when(sampleID_assay == "Blank" ~ "Background0",
                                      sampleID_assay == "VIG " & dilution == 500 ~ "Standard1",
                                      sampleID_assay == "VIG " & dilution == 2000 ~ "Standard2",
                                      sampleID_assay == "VIG " & dilution == 8000 ~ "Standard3",
                                      sampleID_assay == "VIG " & dilution == 32000 ~ "Standard4",
                                      sampleID_assay == "VIG " & dilution == 128000 ~ "Standard5",
                                      TRUE ~ sampleID_assay),
           sampleID_assay = str_remove_all(sampleID_assay, " "),
           sample_type = case_when(sampleID_assay == "NK" ~ "negative",
                                   grepl("Background", sampleID_assay) ~ "standard",
                                   grepl("Standard", sampleID_assay) ~ "standard",
                                   TRUE ~ "unknown"),
           plate_assay = 2)
  
  # Extract isotype on each plate
  isotypes <- docx_extract_tbl(docxtractr::read_docx(dataIn),
                               tbl_number = 4,
                               preserve = TRUE) 
  isotypesClean <- isotypes[-1,]
  
  names(isotypesClean) <- as.character(isotypes[1,])
  names(isotypesClean)[1] <- "row"
  
  isotypes <-
    isotypesClean %>% 
    pivot_longer(c(-row), names_to = "col", values_to = "isotype") %>%
    dplyr::mutate(well = paste(row, col, sep = ""),
           isotype = str_remove(isotype, "Anti-")) %>%
    dplyr::select(well, isotype)
  
  # Extract patientIDs
  patientID_plate_1 <- docx_extract_tbl(docxtractr::read_docx(dataIn),
                                        tbl_number = 5)[,c(1:3)] %>% 
    dplyr::rename(sampleID_assay = Nummer, sampleID_meta = Probe,
           plate_assay = Platte) %>% 
    dplyr::mutate(sampleID_assay = sampleID_assay,
           plate_assay = as.numeric(str_remove(plate_assay, "Platte ")))
  
  patientID_plate_2 <- docx_extract_tbl(docxtractr::read_docx(dataIn),
                                        tbl_number = 6)[,c(1:3)] %>% 
    dplyr::rename(sampleID_assay = Nummer, sampleID_meta = Probe,
           plate_assay = Platte) %>% 
    dplyr::mutate(sampleID_assay = sampleID_assay,
           plate_assay = as.numeric(str_remove(plate_assay, "Platte ")))
  
  patientID <-
    rbind(patientID_plate_1, patientID_plate_2)
  
  # Extract experimental metadata
  metaDataAll <- docx_summary(officer::read_docx(dataIn))
  date <- metaDataAll$text[grepl("Versuchsdurch", metaDataAll$text)] %>% 
    str_remove("Versuchsdurchführung: ") %>% 
    as_datetime(format = "%d.%m.%Y")
  operator <- metaDataAll$text[grepl("Operator", metaDataAll$text)] %>% 
    str_remove("Operatoren: ")
  bead_lot <- metaDataAll$text[grepl("Beads Lot: ", metaDataAll$text)] %>% 
    str_remove("Beads Lot: ")
  assay_buffer_lot <- metaDataAll$text[grepl("Assay-Puffer: 1% BSA/PBS, ", metaDataAll$text)] %>% 
    str_remove("Assay-Puffer: 1% BSA/PBS, ")
  NK_serum <- metaDataAll$text[grepl("NK-Serum:", metaDataAll$text)] %>% 
    str_remove("NK-Serum: ")
  experiment <- metaDataAll$text[grepl("^Experiment", metaDataAll$text)] %>% 
    str_remove("Experiment ")
  
  # Combine all metadata in one dataframe to combine with experimental data
  metadata_joined <-
    rbind(plateLayout_plate_1, plateLayout_plate_2) %>% 
    left_join(isotypes, by = c("well")) %>% 
    left_join(patientID, by = c("sampleID_assay", "plate_assay")) %>% 
    dplyr::mutate(sampleID_meta = if_else(sampleID_assay == "NK",
                                   NK_serum,
                                   sampleID_meta),
           date = date[1],
           experiment = experiment,
           operator = operator,
           bead_lot = bead_lot, # Include information on bead lot in metadata
           file_data_plate_1 = dataFile1,
           file_data_plate_2 = dataFile2,
           file_meta = metafilename)
  return(metadata_joined)
}

join_data_function <- function(mfidata_plate_1, 
                               mfidata_plate_2,
                               metadata_joined) {
  mfidata <- 
    rbind(mfidata_plate_1, mfidata_plate_2)
  output <- 
    metadata_joined %>% 
    left_join(mfidata, by = c("well", "plate_assay"))
  return(output)
}
