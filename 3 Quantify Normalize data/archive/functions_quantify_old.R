## Quantifiy Standards

quantifyMultiplex <- function(ecdata, multiInputList){
  ## select only analytes, that are contained in the ec data file
  analytes <- c(as.character(unique(ecdata$analyte)))
  fitInput <- as.data.frame(dplyr::filter(multiInputList[[i]],
                                          multiInputList[[i]]$analyte%in%analytes)) %>% 
    rename(sample = sampleID_assay) %>% 
    mutate(plate = "plate_1") 
  
  ## Connect ec data with sample data
  datasets <- data_selection(x = fitInput,
                             ecfile = ecdata,
                             byvar.ecfile = c("sample", "analyte"),
                             backname = "Background0",
                             stanname = "Standard",
                             posname = "Controls",
                             fbatch = "plate")
  
  ## Fit to generate standardcurves
  allanalytes <- scluminex(plateid = "newplate", 
                           standard = datasets$plate_1$standard, 
                           background = datasets$plate_1$background,
                           bkg = "constraint", 
                           lfct = c("SSl4"), 
                           fmfi = "data",
                           verbose = TRUE)
  
  ## Determine limit of quantification
  cvAllAnalytes <- loq_cv(allanalytes, max.cv = 0.3)
  limitsOfQuantification <- summary(cvAllAnalytes)
  
  # Quantify concentrations
  concUnknown <- est_conc(allanalytes,
                          datasets$plate_1$unknowns,
                          fmfi = "data",
                          dilution = 1) %>% 
    select(experiment:isotype,
           analyte, data:dil.ub.conc) %>% 
    left_join(limitsOfQuantification, by = "analyte") %>% 
    mutate(remarks = case_when(log10.fitted.conc > ly & log10.fitted.conc < uy~ "ok",
                               log10.fitted.conc <= ly ~ "below_fit", 
                               log10.fitted.conc >= uy ~ "above_fit", 
                               is.na(dil.fitted.conc) & log10(data) <= ly ~ 
                                 "below",
                               is.na(dil.fitted.conc) & log10(data) >= uy ~ 
                                 "above"),
           dil.fitted.final = case_when(remarks == "ok" ~ dil.fitted.conc,
                                        remarks == "below_fit" ~ dil.fitted.conc,
                                        remarks == "above_fit" ~ dil.fitted.conc))
  
  concFinalDil <-
    concUnknown %>% 
    select(dilution_assay, analyte, sampleID_metadata, data, dil.fitted.conc, dil.fitted.final, remarks) %>% 
    mutate(quant_conc_dil = case_when(remarks == "ok" ~ dilution_assay/100*dil.fitted.final,
                                      grepl("above",remarks) & !(dilution_assay %in% c(1000, 6400)) ~ NA_real_,
                                      grepl("below",remarks) & dilution_assay != 100 ~ NA_real_,
                                      TRUE ~ dilution_assay/100*dil.fitted.final)) %>% 
    group_by(analyte, sampleID_metadata) %>% 
    summarise(quant_conc = mean(quant_conc_dil, na.rm = TRUE),
              quant_sd = sd(quant_conc_dil, na.rm = TRUE))
  
  concFinalSingleDil <-
    concUnknown %>% 
    select(dilution_assay, analyte, sampleID_metadata, data, dil.fitted.final, remarks) %>% 
    pivot_wider(names_from = c(dilution_assay), values_from = c(data, dil.fitted.final, remarks)) %>% 
    mutate(quant_single_dil = case_when(grepl("above", remarks_1000) ~ dil.fitted.final_1000*10,
                                        remarks_1000 == "ok" ~ dil.fitted.final_1000*10,
                                        remarks_100 == "ok" ~ dil.fitted.final_100,
                                        grepl("below", remarks_100) ~ dil.fitted.final_100),
           source_single_dil = case_when(grepl("above", remarks_1000) ~ "above",
                                         remarks_1000 == "ok" ~ "1000",
                                         remarks_100 == "ok" ~ "100",
                                         grepl("below", remarks_100) ~ "below")) %>% 
    select(analyte, sampleID_metadata, quant_single_dil,  source_single_dil)
  
  
  concUnknownDil <-
    concUnknown %>% 
    left_join(concFinalDil, by = c("analyte", "sampleID_metadata")) %>% 
    left_join(concFinalSingleDil, by = c("analyte", "sampleID_metadata")) %>% 
    select(-ec, -data.std, -plate, -dilution)
  
  outputListAnalysis <- list(standards = allanalytes,
                             dataOutput = concUnknownDil)
}




quantifyMultiplexReFit <- function(ecdata, multiInputList, allanalytesIn, i){
  ## select only analytes, that are contained in the ec data file
  analytes <- c(as.character(unique(ecdata$analyte)))
  fitInput <- as.data.frame(dplyr::filter(multiInputList[[i]],
                                          multiInputList[[i]]$analyte%in%analytes)) %>% 
    rename(sample = sampleID_assay) %>% 
    mutate(plate = "plate_1") 
  
  ## Connect ec data with sample data
  datasets <- data_selection(x = fitInput,
                             ecfile = ecdata,
                             byvar.ecfile = c("sample", "analyte"),
                             backname = "Background0",
                             stanname = "Standard",
                             posname = "Controls",
                             fbatch = "plate")
  
  ## Fit to generate standardcurves
  allanalytes <- allanalytesIn
  
  ## Determine limit of quantification
  cvAllAnalytes <- loq_cv(allanalytes, max.cv = 0.3)
  limitsOfQuantification <- summary(cvAllAnalytes)
  
  # Quantify concentrations
  concUnknown <- est_conc(allanalytes,
                          datasets$plate_1$unknowns,
                          fmfi = "data",
                          dilution = 1) %>% 
    select(experiment:isotype,
           analyte, data:dil.ub.conc) %>% 
    left_join(limitsOfQuantification, by = "analyte") %>% 
    mutate(remarks = case_when(log10.fitted.conc > ly & log10.fitted.conc < uy~ "ok",
                               log10.fitted.conc <= ly ~ "below_fit", 
                               log10.fitted.conc >= uy ~ "above_fit", 
                               is.na(dil.fitted.conc) & log10(data) <= ly ~ 
                                 "below",
                               is.na(dil.fitted.conc) & log10(data) >= uy ~ 
                                 "above"),
           dil.fitted.final = case_when(remarks == "ok" ~ dil.fitted.conc,
                                        remarks == "below_fit" ~ dil.fitted.conc,
                                        remarks == "above_fit" ~ dil.fitted.conc))
  
  concFinalDil <-
    concUnknown %>% 
    select(dilution_assay, analyte, sampleID_metadata, data, dil.fitted.conc, dil.fitted.final, remarks) %>% 
    mutate(quant_conc_dil = case_when(remarks == "ok" ~ dilution_assay/100*dil.fitted.final,
                                      grepl("above",remarks) & !(dilution_assay %in% c(1000, 6400)) ~ NA_real_,
                                      grepl("below",remarks) & dilution_assay != 100 ~ NA_real_,
                                      TRUE ~ dilution_assay/100*dil.fitted.final)) %>% 
    group_by(analyte, sampleID_metadata) %>% 
    summarise(quant_conc = mean(quant_conc_dil, na.rm = TRUE),
              quant_sd = sd(quant_conc_dil, na.rm = TRUE))
  
  concFinalSingleDil <-
    concUnknown %>% 
    select(dilution_assay, analyte, sampleID_metadata, data, dil.fitted.final, remarks) %>% 
    pivot_wider(names_from = c(dilution_assay), values_from = c(data, dil.fitted.final, remarks)) %>% 
    mutate(quant_single_dil = case_when(grepl("above", remarks_1000) ~ dil.fitted.final_1000*10,
                                        remarks_1000 == "ok" ~ dil.fitted.final_1000*10,
                                        remarks_100 == "ok" ~ dil.fitted.final_100,
                                        grepl("below", remarks_100) ~ dil.fitted.final_100),
           source_single_dil = case_when(grepl("above", remarks_1000) ~ "above",
                                         remarks_1000 == "ok" ~ "1000",
                                         remarks_100 == "ok" ~ "100",
                                         grepl("below", remarks_100) ~ "below")) %>% 
    select(analyte, sampleID_metadata, quant_single_dil,  source_single_dil)
  
  
  concUnknownDil <-
    concUnknown %>% 
    left_join(concFinalDil, by = c("analyte", "sampleID_metadata")) %>% 
    left_join(concFinalSingleDil, by = c("analyte", "sampleID_metadata")) %>% 
    select(-ec, -data.std, -plate, -dilution)
  
  outputListAnalysis <- list(standards = allanalytes,
                             dataOutput = concUnknownDil)
}


##### Dev Quantify ELISA
quantifyELISA <- function(ecdata, multiInputList){
  ## select only analytes, that are contained in the ec data file
  analytes <- c(as.character(unique(ecdata$analyte)))
  fitInput <- as.data.frame(dplyr::filter(multiInputList[[i]],
                                          multiInputList[[i]]$analyte%in%analytes)) %>% 
    rename(sample = sampleID_assay) %>% 
    mutate(plate = "plate_1") 
  
  ## Connect ec data with sample data
  datasets <- data_selection(x = fitInput,
                             ecfile = ecdata,
                             byvar.ecfile = c("sample", "analyte"),
                             backname = "Background0",
                             stanname = "Standard",
                             posname = "Controls",
                             fbatch = "plate")
  
  ## Fit to generate standardcurves
  allanalytes <- scluminex(plateid = "newplate", 
                           standard = datasets$plate_1$standard, 
                           background = datasets$plate_1$background,
                           bkg = "constraint", 
                           lfct = c("SSl4"), 
                           fmfi = "data",
                           verbose = TRUE)
  
  plot(allanalytes)
  
  ## Determine limit of quantification
  cvAllAnalytes <- loq_cv(allanalytes, max.cv = 0.3)
  limitsOfQuantification <- summary(cvAllAnalytes)
  
  # Quantify concentrations
  concUnknown <- est_conc(allanalytes,
                          datasets$plate_1$unknowns,
                          fmfi = "data",
                          dilution = 1) %>% 
    select(experiment:isotype,
           analyte, data:dil.ub.conc) %>% 
    left_join(limitsOfQuantification, by = "analyte") %>% 
    mutate(remarks = case_when(log10.fitted.conc > ly & log10.fitted.conc < uy~ "ok",
                               log10.fitted.conc <= ly ~ "below_fit", 
                               log10.fitted.conc >= uy ~ "above_fit", 
                               is.na(dil.fitted.conc) & log10(data) <= ly ~ 
                                 "below",
                               is.na(dil.fitted.conc) & log10(data) >= uy ~ 
                                 "above"),
           dil.fitted.final = case_when(remarks == "ok" ~ dil.fitted.conc,
                                        remarks == "below_fit" ~ dil.fitted.conc,
                                        remarks == "above_fit" ~ dil.fitted.conc))
  
  concFinalDil <-
    concUnknown %>% 
    select(dilution_assay, analyte, sampleID_metadata, data, dil.fitted.conc, dil.fitted.final, remarks) %>% 
    mutate(quant_conc_dil = case_when(remarks == "ok" ~ dilution_assay/100*dil.fitted.final,
                                      grepl("above",remarks) & !(dilution_assay %in% c(1000, 6400)) ~ NA_real_,
                                      grepl("below",remarks) & dilution_assay != 100 ~ NA_real_,
                                      TRUE ~ dilution_assay/100*dil.fitted.final)) %>% 
    group_by(analyte, sampleID_metadata) %>% 
    summarise(quant_conc = mean(quant_conc_dil, na.rm = TRUE),
              quant_sd = sd(quant_conc_dil, na.rm = TRUE))
  
  
  concUnknownDil <-
    concUnknown %>% 
    left_join(concFinalDil, by = c("analyte", "sampleID_metadata")) %>% 
    mutate(quant_single_dil = NA_real_,  
           source_single_dil = NA_real_) %>% 
    select(-ec, -data.std, -plate, -dilution)
  
  outputListAnalysis <- list(standards = allanalytes,
                             dataOutput = concUnknownDil)
}

quantifyMultiplexReFit <- function(ecdata, multiInputList, allanalytesIn, i){
  ## select only analytes, that are contained in the ec data file
  analytes <- c(as.character(unique(ecdata$analyte)))
  fitInput <- as.data.frame(dplyr::filter(multiInputList[[i]],
                                          multiInputList[[i]]$analyte%in%analytes)) %>% 
    rename(sample = sampleID_assay) %>% 
    mutate(plate = "plate_1") 
  
  ## Connect ec data with sample data
  datasets <- data_selection(x = fitInput,
                             ecfile = ecdata,
                             byvar.ecfile = c("sample", "analyte"),
                             backname = "Background0",
                             stanname = "Standard",
                             posname = "Controls",
                             fbatch = "plate")
  
  ## Fit to generate standardcurves
  allanalytes <- allanalytesIn
  
  ## Determine limit of quantification
  cvAllAnalytes <- loq_cv(allanalytes, max.cv = 0.3)
  limitsOfQuantification <- summary(cvAllAnalytes)
  
  # Quantify concentrations
  concUnknown <- est_conc(allanalytes,
                          datasets$plate_1$unknowns,
                          fmfi = "data",
                          dilution = 1) %>% 
    select(experiment:isotype,
           analyte, data:dil.ub.conc) %>% 
    left_join(limitsOfQuantification, by = "analyte") %>% 
    mutate(remarks = case_when(log10.fitted.conc > ly & log10.fitted.conc < uy~ "ok",
                               log10.fitted.conc <= ly ~ "below_fit", 
                               log10.fitted.conc >= uy ~ "above_fit", 
                               is.na(dil.fitted.conc) & log10(data) <= ly ~ 
                                 "below",
                               is.na(dil.fitted.conc) & log10(data) >= uy ~ 
                                 "above"),
           dil.fitted.final = case_when(remarks == "ok" ~ dil.fitted.conc,
                                        remarks == "below_fit" ~ dil.fitted.conc,
                                        remarks == "above_fit" ~ dil.fitted.conc))
  
  concFinalDil <-
    concUnknown %>% 
    select(dilution_assay, analyte, sampleID_metadata, data, dil.fitted.conc, dil.fitted.final, remarks) %>% 
    mutate(quant_conc_dil = case_when(remarks == "ok" ~ dilution_assay/100*dil.fitted.final,
                                      grepl("above",remarks) & !(dilution_assay %in% c(1000, 6400)) ~ NA_real_,
                                      grepl("below",remarks) & dilution_assay != 100 ~ NA_real_,
                                      TRUE ~ dilution_assay/100*dil.fitted.final)) %>% 
    group_by(analyte, sampleID_metadata) %>% 
    summarise(quant_conc = mean(quant_conc_dil, na.rm = TRUE),
              quant_sd = sd(quant_conc_dil, na.rm = TRUE))
  
  concFinalDil <-
    concUnknown %>% 
    select(dilution_assay, analyte, sampleID_metadata, data, dil.fitted.conc, dil.fitted.final, remarks) %>% 
    mutate(quant_conc_dil = case_when(remarks == "ok" ~ dilution_assay/100*dil.fitted.final,
                                      grepl("above",remarks) & !(dilution_assay %in% c(1000, 6400)) ~ NA_real_,
                                      grepl("below",remarks) & dilution_assay != 100 ~ NA_real_,
                                      TRUE ~ dilution_assay/100*dil.fitted.final)) %>% 
    group_by(analyte, sampleID_metadata) %>% 
    summarise(quant_conc = mean(quant_conc_dil, na.rm = TRUE),
              quant_sd = sd(quant_conc_dil, na.rm = TRUE))
  
  
  concUnknownDil <-
    concUnknown %>% 
    left_join(concFinalDil, by = c("analyte", "sampleID_metadata")) %>% 
    mutate(quant_single_dil = NA_real_,  
           source_single_dil = NA_real_) %>% 
    select(-ec, -data.std, -plate, -dilution)
  
  outputListAnalysis <- list(standards = allanalytes,
                             dataOutput = concUnknownDil)
}
