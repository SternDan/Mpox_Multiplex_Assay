## Quantifiy Standards
quantifyMultiplex <- function(ecdata, multiInputList){
  ## select only analytes, that are contained in the ec data file
  analytes <- c(as.character(unique(ecdata$analyte)))
  fitInput <- as.data.frame(dplyr::filter(multiInputList,
                                          multiInputList$analyte%in%analytes)) %>% 
    dplyr::rename(sample = sampleID_assay) %>% 
    dplyr::mutate(plate = "plate_1") 
  
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
                           bkg = "ignore", 
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
   # dplyr::select(experiment:isotype,
  #         analyte, sample, data:dil.ub.conc) %>% 
    dplyr::left_join(limitsOfQuantification, by = "analyte") %>% 
    dplyr::mutate(remarks = case_when(log10.fitted.conc > ly & log10.fitted.conc < uy~ "ok",
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
    dplyr::select(dilution_assay, analyte, sample, data, dil.fitted.conc, dil.fitted.final, remarks) %>% 
    dplyr::mutate(quant_conc_dil = case_when(remarks == "ok" ~ dilution_assay/100*dil.fitted.final,
                                      grepl("above",remarks) & !(dilution_assay %in% c(1000, 6400)) ~ NA_real_,
                                      grepl("below",remarks) & dilution_assay != 100 ~ NA_real_,
                                      TRUE ~ dilution_assay/100*dil.fitted.final)) %>% 
    group_by(analyte, sample) %>% 
    dplyr::summarise(quant_conc = mean(quant_conc_dil, na.rm = TRUE),
              quant_sd = sd(quant_conc_dil, na.rm = TRUE))
  
  concFinalSingleDil <-
    concUnknown %>% 
    dplyr::select(dilution_assay, analyte, sample, data, dil.fitted.final, remarks) %>% 
    pivot_wider(names_from = c(dilution_assay), values_from = c(data, dil.fitted.final, remarks), values_fn = NULL) %>% ## Neu 26.09.2023 values_fn zu Null
    dplyr::mutate(quant_single_dil = case_when(grepl("above", remarks_1000) ~ dil.fitted.final_1000*10,
                                        remarks_1000 == "ok" ~ dil.fitted.final_1000*10,
                                        remarks_100 == "ok" ~ dil.fitted.final_100,
                                        grepl("below", remarks_100) ~ dil.fitted.final_100),
           source_single_dil = case_when(grepl("above", remarks_1000) ~ "above",
                                         remarks_1000 == "ok" ~ "1000",
                                         remarks_100 == "ok" ~ "100",
                                         grepl("below", remarks_100) ~ "below")) %>% 
    dplyr::select(analyte, sample, quant_single_dil,  source_single_dil)
  
  
  concUnknownDil <-
    concUnknown %>% 
    dplyr::left_join(concFinalDil, by = c("analyte", "sample")) %>% 
    dplyr::left_join(concFinalSingleDil, by = c("analyte", "sample")) #%>% 
   # dplyr::select(-ec, -plate, -dilution)
  
  outputListAnalysis <- list(standards = allanalytes,
                             dataOutput = concUnknownDil)
}






quantifyMultiplexReFit <- function(ecdata, multiInputList, allanalytesIn, i){
  ## select only analytes, that are contained in the ec data file
  analytes <- c(as.character(unique(ecdata$analyte)))
  fitInput <- as.data.frame(dplyr::filter(multiInputList,
                                          multiInputList$analyte%in%analytes)) %>% 
    dplyr::rename(sample = sampleID_assay) %>% 
    dplyr::mutate(plate = "plate_1") 
  
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
   # dplyr::select(experiment:isotype,
  #         analyte, sample, data:dil.ub.conc) %>% 
    dplyr::left_join(limitsOfQuantification, by = "analyte") %>% 
    dplyr::mutate(remarks = case_when(log10.fitted.conc > ly & log10.fitted.conc < uy~ "ok",
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
    dplyr::select(dilution_assay, analyte, sample, data, dil.fitted.conc, dil.fitted.final, remarks) %>% 
    dplyr::mutate(quant_conc_dil = case_when(remarks == "ok" ~ dilution_assay/100*dil.fitted.final,
                                      grepl("above",remarks) & !(dilution_assay %in% c(1000, 6400)) ~ NA_real_,
                                      grepl("below",remarks) & dilution_assay != 100 ~ NA_real_,
                                      TRUE ~ dilution_assay/100*dil.fitted.final)) %>% 
    group_by(analyte, sample) %>% 
    dplyr::summarise(quant_conc = mean(quant_conc_dil, na.rm = TRUE),
              quant_sd = sd(quant_conc_dil, na.rm = TRUE))
  
  concFinalSingleDil <-
    concUnknown %>% 
    dplyr::select(dilution_assay, analyte, sample, data, dil.fitted.final, remarks) %>% 
    pivot_wider(names_from = c(dilution_assay), values_from = c(data, dil.fitted.final, remarks), values_fn = NULL) %>% 
    dplyr::mutate(quant_single_dil = case_when(grepl("above", remarks_1000) ~ dil.fitted.final_1000*10,
                                        remarks_1000 == "ok" ~ dil.fitted.final_1000*10,
                                        remarks_100 == "ok" ~ dil.fitted.final_100,
                                        grepl("below", remarks_100) ~ dil.fitted.final_100),
           source_single_dil = case_when(grepl("above", remarks_1000) ~ "above",
                                         remarks_1000 == "ok" ~ "1000",
                                         remarks_100 == "ok" ~ "100",
                                         grepl("below", remarks_100) ~ "below")) %>% 
    dplyr::select(analyte, sample, quant_single_dil,  source_single_dil)
  
  
  concUnknownDil <-
    concUnknown %>% 
    dplyr::left_join(concFinalDil, by = c("analyte", "sample")) %>% 
    dplyr::left_join(concFinalSingleDil, by = c("analyte", "sample")) #%>% 
    #dplyr::select(-ec, -plate, -dilution)
  
  outputListAnalysis <- list(standards = allanalytes,
                             dataOutput = concUnknownDil)
}

quantifyMultiplexNK <- function(ecdata, multiInputList){
  ## select only analytes, that are contained in the ec data file
  analytes <- c(as.character(unique(ecdata$analyte)))
  fitInput <- as.data.frame(dplyr::filter(multiInputList,
                                          multiInputList$analyte%in%analytes)) %>% 
    dplyr::rename(sample = sampleID_assay) %>% 
    dplyr::mutate(plate = "plate_1") 
  
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
                           bkg = "ignore", # Deviation from other analysis: here ignore, else constraint
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
    dplyr::select(experiment:isotype,
           analyte, data:dil.ub.conc) %>% 
    dplyr::left_join(limitsOfQuantification, by = "analyte") %>% 
    dplyr::mutate(remarks = case_when(log10.fitted.conc > ly & log10.fitted.conc < uy~ "ok",
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
    dplyr::select(dilution_assay, analyte, sample, data, dil.fitted.conc, dil.fitted.final, remarks) %>% 
    dplyr::mutate(quant_conc_dil = case_when(remarks == "ok" ~ dilution_assay/100*dil.fitted.final,
                                      grepl("above",remarks) & !(dilution_assay %in% c(1000, 6400)) ~ NA_real_,
                                      grepl("below",remarks) & dilution_assay != 100 ~ NA_real_,
                                      TRUE ~ dilution_assay/100*dil.fitted.final)) %>% 
    group_by(analyte, sample) %>% 
    dplyr::summarise(quant_conc = mean(quant_conc_dil, na.rm = TRUE),
              quant_sd = sd(quant_conc_dil, na.rm = TRUE))
  
  concFinalSingleDil <-
    concUnknown %>% 
    dplyr::select(dilution_assay, analyte, sample, data, dil.fitted.final, remarks) %>% 
    pivot_wider(names_from = c(dilution_assay), values_from = c(data, dil.fitted.final, remarks)) %>% 
    dplyr::mutate(quant_single_dil = dil.fitted.final_100,
           source_single_dil = case_when(remarks_100 == "ok" ~ "100",
                                         grepl("below", remarks_100) ~ "below")) %>% 
    dplyr::select(analyte, sample, quant_single_dil,  source_single_dil)
  
  
  concUnknownDil <-
    concUnknown %>% 
    dplyr::left_join(concFinalDil, by = c("analyte", "sample")) %>% 
    dplyr::left_join(concFinalSingleDil, by = c("analyte", "sample")) %>% 
    dplyr::select(-ec, -data.std, -plate, -dilution)
  
  outputListAnalysis <- list(standards = allanalytes,
                             dataOutput = concUnknownDil)
}

quantifyMultiplexReFitNK <- function(ecdata, multiInputList, allanalytesIn){
  ## select only analytes, that are contained in the ec data file
  analytes <- c(as.character(unique(ecdata$analyte)))
  fitInput <- as.data.frame(dplyr::filter(multiInputList,
                                          multiInputList$analyte%in%analytes)) %>% 
    dplyr::rename(sample = sampleID_assay) %>% 
    dplyr::mutate(plate = "plate_1") 
  
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
    dplyr::select(experiment:isotype,
           analyte, data:dil.ub.conc) %>% 
    dplyr::left_join(limitsOfQuantification, by = "analyte") %>% 
    dplyr::mutate(remarks = case_when(log10.fitted.conc > ly & log10.fitted.conc < uy~ "ok",
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
    dplyr::select(dilution_assay, analyte, sample, data, dil.fitted.conc, dil.fitted.final, remarks) %>% 
    dplyr::mutate(quant_conc_dil = case_when(remarks == "ok" ~ dilution_assay/100*dil.fitted.final,
                                      grepl("above",remarks) & !(dilution_assay %in% c(1000, 6400)) ~ NA_real_,
                                      grepl("below",remarks) & dilution_assay != 100 ~ NA_real_,
                                      TRUE ~ dilution_assay/100*dil.fitted.final)) %>% 
    group_by(analyte, sample) %>% 
    dplyr::summarise(quant_conc = mean(quant_conc_dil, na.rm = TRUE),
              quant_sd = sd(quant_conc_dil, na.rm = TRUE))
  
  concFinalSingleDil <-
    concUnknown %>% 
    dplyr::select(dilution_assay, analyte, sample, data, dil.fitted.final, remarks) %>% 
    pivot_wider(names_from = c(dilution_assay), values_from = c(data, dil.fitted.final, remarks)) %>% 
    dplyr::mutate(quant_single_dil = dil.fitted.final_100,
           source_single_dil = case_when(remarks_100 == "ok" ~ "100",
                                         grepl("below", remarks_100) ~ "below")) %>% 
    dplyr::select(analyte, sample, quant_single_dil,  source_single_dil)
  
  
  concUnknownDil <-
    concUnknown %>% 
    dplyr::left_join(concFinalDil, by = c("analyte", "sample")) %>% 
    dplyr::left_join(concFinalSingleDil, by = c("analyte", "sample")) %>% 
    dplyr::select(-ec, -data.std, -plate, -dilution)
  
  outputListAnalysis <- list(standards = allanalytes,
                             dataOutput = concUnknownDil)
}

generate_quantify_output_function <- function(ecdata, dataInput, threshold){
  quantInputIgG_plate_1 <-
    dataInput%>% 
    dplyr::rename(dilution_assay = dilution) %>% 
    filter(isotype == "IgG" & plate_assay == 1)
  
  quantInputIgG_plate_2 <-
    dataInput%>% 
    dplyr::rename(dilution_assay = dilution) %>% 
    filter(isotype == "IgG" & plate_assay == 2)
  
  quantInputIgM_plate_1 <-
    dataInput%>% 
    dplyr::rename(dilution_assay = dilution) %>% 
    filter(isotype == "IgM"& plate_assay == 1)
  
  quantInputIgM_plate_2 <-
    dataInput%>% 
    dplyr::rename(dilution_assay = dilution) %>% 
    filter(isotype == "IgM" & plate_assay == 2)
  
  quantifiedIgG_plate_1 <- quantifyMultiplex(ecdata = ecdata,
                                             multiInputList = quantInputIgG_plate_1)
  
  quantifiedIgG_plate_2 <- quantifyMultiplex(ecdata = ecdata,
                                             multiInputList = quantInputIgG_plate_2)
  
  quantifiedIgM_plate_1 <- quantifyMultiplexReFit(ecdata = ecdata,
                                                  multiInputList = quantInputIgM_plate_1,
                                                  allanalytesIn = quantifiedIgG_plate_1$standards,
                                                  1)
  
  quantifiedIgM_plate_2 <- quantifyMultiplexReFit(ecdata = ecdata,
                                                  multiInputList = quantInputIgM_plate_2,
                                                  allanalytesIn = quantifiedIgG_plate_1$standards,
                                                  1)
  
  
  standard_plate_1 <- quantifiedIgG_plate_1$standards  
  standard_plate_2 <- quantifiedIgG_plate_2$standards 
  
  standard_info_plate_1 <- do.call(rbind, lapply(standard_plate_1, '[[', 2)) %>% 
    add_column(plate_assay = 1, 
               date = unique(dataInput$date), 
               experiment = unique(dataInput$experiment), 
               operator = unique(dataInput$operator), 
               bead_lot = unique(dataInput$bead_lot))
  standard_info_plate_2 <- do.call(rbind, lapply(standard_plate_1, '[[', 2)) %>% 
    add_column(plate_assay = 2, 
               date = unique(dataInput$date), 
               experiment = unique(dataInput$experiment), 
               operator = unique(dataInput$operator), 
               bead_lot = unique(dataInput$bead_lot))
  standard_info <-
    bind_rows(standard_info_plate_1, standard_info_plate_2)
  
  dataInputQuant <- 
    rbind(quantifiedIgG_plate_1$dataOutput, 
          quantifiedIgG_plate_2$dataOutput, 
          quantifiedIgM_plate_1$dataOutput,
          quantifiedIgM_plate_2$dataOutput)  %>% 
    dplyr::mutate(analyte = factor(analyte, levels = c("A27L",
                                                "A29", 
                                                "L1R", 
                                                "M1",
                                                "D8L",
                                                "E8",
                                                "H3L",
                                                "A33R", 
                                                "A35R",
                                                "B5R", 
                                                "B6", 
                                                "A5L",
                                                "ATI-C", 
                                                "ATI-N",
                                                "VACV", 
                                                "Delta"), ordered = TRUE)) %>% 
    group_by(analyte, experiment, plate_assay) %>% 
    dplyr::mutate(dataIn = case_when(is.na(quant_conc) & grepl("above", remarks) ~ max(quant_conc, na.rm = TRUE),
                              is.na(quant_conc) & grepl("below", remarks) ~ min(quant_conc, na.rm = TRUE),
                              TRUE ~ quant_conc),
           dataIn = if_else(dataIn < 1, 1, dataIn)) %>% 
    ungroup()
  
  ##
  # 03) Determine serostatus for IgG and IgM compared to the delta results
  dataInputQuantCat <-
    dataInputQuant %>%  
    filter(dilution_assay == 100) %>% 
    dplyr::left_join(threshold, by = c("isotype" = "assay", "analyte" = "antigene")) %>% 
    dplyr::mutate(serostatus = if_else(dataIn > median, 1, 0),
           serostatus_cat = case_when(dataIn > high ~ "positive",
                                      dataIn > median & dataIn <= high ~ "borderline positive",
                                      dataIn > low & dataIn <= median ~ "borderline negative",
                                      dataIn <= low ~ "negative"),
           serostatus_cat = factor(serostatus_cat, levels = c("negative", "borderline negative",
                                                              "borderline positive", "positive"),
                                   ordered = TRUE)) 
  
  dataInputQuantCatDelta <-
    dataInputQuantCat %>% 
    filter(dilution_assay == 100 & analyte == "Delta") %>% 
    dplyr::select(experiment, plate_assay, isotype, sampleID_meta, dataIn,
                  serostatus, serostatus_cat)
  
  dataInput <-
    dataInputQuantCat %>% 
    dplyr::left_join(dataInputQuantCatDelta, by = c("experiment", "plate_assay", "isotype", "sampleID_meta"), suffix = c("", ".delta")) %>% 
    unique()
  
  outlist <- list()
  outlist[[1]] <- standard_plate_1
  outlist[[2]] <- standard_plate_2
  outlist[[3]] <- dataInput
  outlist[[4]] <- standard_info
  return(outlist)
}


