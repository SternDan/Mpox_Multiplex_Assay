## Define function for data input and merge with experimental metadata ----
read_multiplex <- function(filename, metafilename){
  # Read mfi data, calculate and normalize Delta
  mfidata <- 
    import(filename, skip = 7) %>%
    filter(...2 != "Well") %>%
    rename(type = ...1, well = ...2) %>%
    mutate_at(vars(-well, -type), ~(as.numeric(sub(",", ".", ., fixed = TRUE)))) %>%
    pivot_longer(c(-well, -type), names_to = "analyte", values_to = "data") %>%
    mutate(analyte = str_replace(analyte, " \\([^\\)]+\\)", "")) %>%
    dplyr::select(-type) %>%
    unique() %>%
    pivot_wider(names_from = analyte,
                values_from = data) %>%
    mutate(Delta = VACV - Hep2,
           Delta = Delta + abs(min(Delta))) %>%
    pivot_longer(cols = - well, names_to = "analyte", values_to = "data")
  
  # Read date when plate was read
  date_character <-import(filename, skip = 1, col_names = FALSE, n_max = 1) %>% 
    mutate(date = str_remove(`...1`, "Acquisition Date: "),
           date = as_datetime(date, format = "%d-%b-%Y, %H:%M")) %>% 
    dplyr::select(-`...1`) %>% 
    pull()
  
  # Read metadata contained in metadata file
  metadata <- import(metafilename,
                     sheet = "metadata") %>% 
    select(experiment = Assay, plate = Plate, batch = Batch)
  
  layout <- import(metafilename,
                   sheet = "layout") %>%
    pivot_longer(c(-...1), names_to = "col", values_to = "sample") %>%
    mutate(well = paste(`...1`, col, sep = "")) %>%
    dplyr::select(well, sampleID_assay = sample) %>% 
    mutate(sample_type = case_when(sampleID_assay == "NK" ~ "negative",
                                   grepl("Background", sampleID_assay) ~ "standard",
                                   grepl("Standard", sampleID_assay) ~ "standard",
                                   TRUE ~ "unknown"))
  
  isotype <- import(metafilename,
                    sheet = "isotype") %>%
    pivot_longer(c(-...1), names_to = "col", values_to = "isotype") %>%
    mutate(well = paste(`...1`, col, sep = "")) %>%
    dplyr::select(well, isotype)
  
  dilution <- import(metafilename,
                     sheet = "dilution") %>%
    pivot_longer(c(-...1), names_to = "col", values_to = "dilInput") %>%
    mutate(well = paste(`...1`, col, sep = "")) %>%
    dplyr::select(well, dilInput)
  
  # Combine metadata in one file
  metadata_joined <-
    layout %>%
    left_join(isotype, by = c("well")) %>%
    left_join(dilution, by = c("well")) %>%
    mutate(date = date_character,
           experiment = metadata$experiment,
           plate = metadata$plate,
           batch = metadata$batch,
           filename = filename,
           assaytype = "Multiplex") %>% 
    dplyr::select(experiment, filename, date, plate, batch, well, assaytype, sampleID_assay,
           dilution = dilInput, isotype)
  
  output <- 
    metadata_joined %>% 
    left_join(mfidata, by = "well")
  return(output)
}
