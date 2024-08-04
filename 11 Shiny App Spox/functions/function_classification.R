classify_sero_function <- function(dataInput,
                                    lDAPanelsAll,
                                    lDAPanelsDelta,
                                    lDAPanelsHigh,
                                    Combined_spec,
                                    ATI_N_spec,
                                    aucCombinedIgGPos){
                                
  dataInputLDAAll <- dataInput %>% 
    dplyr::filter(dilution_assay == 100 & isotype == "IgG") %>%
    dplyr::filter(analyte != "VACV") %>% 
    dplyr::filter(!(analyte %in% c("L1R", "M1"))) %>% # New 2023-08-21: Deselect L1R and M1
    dplyr::select(experiment, plate_assay, bead_lot, isotype, sampleID_meta,
                  analyte, serostatus, dataIn, serostatus.delta, dataIn.delta) 
   
   dataInputLDADelta <- dataInput %>% 
     dplyr::filter(dilution_assay == 100 & isotype == "IgG") %>% 
     dplyr::filter(assaytype == "Multiplex") %>% 
     dplyr::filter(analyte != "VACV") %>% 
     dplyr::filter(!(analyte %in% c("L1R", "M1"))) %>% # New 2023-08-21: Deselect L1R and M1
     dplyr::filter(serostatus.delta == 1) %>% 
     dplyr::select(experiment, plate_assay, bead_lot, isotype, sampleID_meta, 
                   analyte, serostatus, dataIn, serostatus.delta, dataIn.delta)
   
   dataInputLDAHigh <- dataInput %>% 
     dplyr::filter(dilution_assay == 100 & isotype == "IgG") %>% 
     dplyr::filter(assaytype == "Multiplex") %>% 
     dplyr::filter(analyte != "VACV") %>% 
     dplyr::filter(!(analyte %in% c("L1R", "M1"))) %>% # New 2023-08-21: Deselect L1R and M1
     dplyr::filter(log10(dataIn.delta) >= 2.08) %>% 
     dplyr::select(experiment, plate_assay, bead_lot, isotype, sampleID_meta, 
                   analyte, serostatus, dataIn, serostatus.delta, dataIn.delta) 
   
   
   dataInputLDAAllwide <-
     dataInputLDAAll %>% 
     dplyr::select(experiment, plate_assay, bead_lot, sampleID_meta, analyte, dataIn) %>% 
     mutate(dataIn = log10(dataIn)) %>% 
     pivot_wider(names_from = analyte, values_from = dataIn, values_fn = NULL)
   
   dataInputLDADeltawide <-
     dataInputLDADelta %>% 
     dplyr::select(experiment, plate_assay, bead_lot, sampleID_meta, analyte, dataIn) %>% 
     mutate(dataIn = log10(dataIn)) %>% 
     pivot_wider(names_from = analyte, values_from = dataIn, values_fn = NULL)
   
   dataInputLDAHighwide <-
     dataInputLDAHigh %>% 
     dplyr::select(experiment, plate_assay, bead_lot, sampleID_meta, analyte, dataIn) %>% 
     mutate(dataIn = log10(dataIn)) %>% 
     pivot_wider(names_from = analyte, values_from = dataIn, values_fn = NULL)
   
   
   predictAll <- predict(lDAPanelsAll , dataInputLDAAllwide[5:ncol(dataInputLDAAllwide)])
   predictDelta <- predict(lDAPanelsDelta , dataInputLDADeltawide[5:ncol(dataInputLDADeltawide)])
   predictHigh <- predict(lDAPanelsHigh , dataInputLDAHighwide[5:ncol(dataInputLDAHighwide)])
   
   
   dataInputGLMAllwide <-
     dataInputLDAAll %>% 
     dplyr::select(experiment, plate_assay, sampleID_meta, analyte, dataIn) %>% 
     mutate(dataIn = log10(dataIn),
            analyte = str_replace(analyte, "-", ".")) %>% 
     pivot_wider(names_from = analyte, values_from = dataIn, values_fn = NULL)
   
   #predictGLM <- predict(aucCombinedIgGPos, dataInputGLMAllwide)
   
  # predictionGLMDelta <- 
  #   data.frame(GLM_value =
  #                predict(aucCombinedIgGPos, dataInputGLMAllwide)) %>% 
  #   mutate(Combined_infected = (if_else(GLM_value > Combined_spec, "infected", "uninfected")))
   
   ATI_infected <-
     dataInputGLMAllwide %>% 
     mutate(ATI_infected = (if_else(ATI.N > ATI_N_spec, "infected", "uninfected"))) %>% 
     dplyr::select(sampleID_meta, ATI_infected)#,
   #ATI_infected = if_else(is.na(MPXV), ("Seronegative"),  ATI_infected))
   
   ##
   # Generate prediction for different groups
   # MPXV
   predictAllpost <- predictAll$posterior[,4] 
   predictDeltapost <- predictDelta$posterior[,4] 
   predictHighpost <- predictHigh$posterior[,4] 
   
   # Pre
   predictAllpostPre <- predictAll$posterior[,1] 
   predictDeltapostPre <- predictDelta$posterior[,1] 
   predictHighpostPre <- predictHigh$posterior[,1] 
   
   # MVA
   predictAllpostMVA <- predictAll$posterior[,2] 
   predictDeltapostMVA <- predictDelta$posterior[,2] 
   predictHighpostMVA <- predictHigh$posterior[,2] 
   
   # CPXV
   predictAllpostCPXV <- predictAll$posterior[,3] 
   predictDeltapostCPXV <- predictDelta$posterior[,3] 
   predictHighpostCPXV <- predictHigh$posterior[,3] 
   
   LDAAll <-
     dataInputLDAAllwide %>% 
     dplyr::select(experiment, plate_assay, sampleID_meta) %>% 
     add_column(MPXV_all = predictAllpost,
                Pre_all = predictAllpostPre,
                MVA_all = predictAllpostMVA,
                CPXV_all = predictAllpostCPXV)
   
   LDADelta <-
     dataInputLDADeltawide %>% 
     dplyr::select(experiment, plate_assay, sampleID_meta) %>% 
     add_column(MPXV_delta = predictDeltapost,
                Pre_delta = predictDeltapostPre,
                MVA_delta = predictDeltapostMVA,
                CPXV_delta = predictDeltapostCPXV)
   
   LDAHigh <-
     dataInputLDAHighwide %>% 
     dplyr::select(experiment, plate_assay, sampleID_meta) %>% 
     add_column(MPXV_high = predictHighpost,
                Pre_high = predictHighpostPre,
                MVA_high = predictHighpostMVA,
                CPXV_high = predictHighpostCPXV)
   
   
   
   unique(dataInput$analyte)
   
   ##
   # Generate final data output
   
   # Select 
   
   dataOutputIgG <-
     dataInput %>% 
     filter(isotype == "IgG") %>% 
     filter(analyte != "VACV") %>% 
     dplyr::select(experiment, plate_assay, bead_lot, date, operator, 
            sampleID_meta, isotype, highBG, lowIg,
            analyte, data, dataIn, serostatus, serostatus_cat, serostatus.delta, serostatus_cat.delta) %>% 
     pivot_wider(names_from = analyte, values_from = c(data, dataIn, serostatus, serostatus_cat)) %>% 
     dplyr::left_join(LDAAll, by = c("experiment", "sampleID_meta", "plate_assay")) %>% 
     dplyr::left_join(LDADelta, by = c("experiment", "sampleID_meta", "plate_assay")) %>% 
     dplyr::left_join(LDAHigh, by = c("experiment", "sampleID_meta", "plate_assay")) %>% 
     mutate(MPXV_all = round(MPXV_all, digits = 2),
            MPXV_delta = round(MPXV_delta, digits = 2),
            MPXV_high = round(MPXV_high, digits = 2),
            Pre_all = round(Pre_all, digits = 2),
            Pre_delta = round(Pre_delta, digits = 2),
            Pre_high = round(Pre_high, digits = 2),
            MVA_all = round(MVA_all, digits = 2),
            MVA_delta = round(MVA_delta, digits = 2),
            MVA_high = round(MVA_high, digits = 2),
            CPXV_all = round(CPXV_all, digits = 2),
            CPXV_delta = round(CPXV_delta, digits = 2),
            CPXV_high = round(CPXV_high, digits = 2),
            Final_all = case_when(MPXV_all > 0.95 ~ "MPXV",
                                  Pre_all > 0.95 ~ "Pre",
                                  MVA_all > 0.95 ~ "MVA",
                                  CPXV_all > 0.95 ~ "CPXV",
                                  TRUE ~ "Ambiguous"),
            Final_delta = case_when(MPXV_delta > 0.95 ~ "MPXV",
                                    Pre_delta > 0.95 ~ "Pre",
                                    MVA_delta > 0.95 ~ "MVA",
                                    CPXV_delta > 0.95 ~ "CPXV",
                                    is.na(MPXV_delta) ~ NA_character_,
                                    TRUE ~ "Ambiguous"),
            Final_high = case_when(MPXV_high > 0.95 ~ "MPXV",
                                   Pre_high > 0.95 ~ "Pre",
                                   MVA_high > 0.95 ~ "MVA",
                                   CPXV_high > 0.95 ~ "CPXV",
                                   is.na(MPXV_high) ~ NA_character_,
                                   TRUE ~ "Ambiguous"),
             Final_combined = case_when(!is.na(Final_high) 
                                         & Final_high != "Ambiguous" ~ Final_high,
                                        !is.na(Final_delta)
                                        & Final_delta != "Ambiguous" ~ Final_delta,
                                         Final_all != "Ambiguous" ~ Final_all,
                                         is.na(Final_delta)  ~ "seronegative"),
          # Final_combined = case_when(Final_high != "NA" 
          #                             & Final_high == Final_delta ~ Final_high,
          #                             Final_high != "NA" & Final_high != "Ambiguous"
          #                             & Final_high != Final_delta ~ Final_high,
          #                             Final_high != "NA" & Final_high == "Ambiguous"
          #                             & Final_high != Final_delta ~ Final_delta,
          #                             Final_high == "NA" & Final_delta != "NA" ~ Final_delta,
          #                             Final_delta == "NA" ~ "seronegative"),
          #  Confidence = case_when(Final_high != "NA" 
          #                         & Final_high == Final_delta ~ "Highest",
          #                         Final_high != "NA" & Final_high != "Ambiguous"
          #                         & Final_high != Final_delta ~ "High",
          #                         Final_high != "NA" & Final_high == "Ambiguous"
          #                         & Final_high != Final_delta ~ "Medium",
          #                         Final_high == "NA" & Final_delta != "NA" ~ "Medium",
          #                         Final_delta == "NA" ~ "seronegative"),
            ratio_A29_A27L = dataIn_A29/dataIn_A27L,
            ratio_E8_D8L = dataIn_E8/dataIn_D8L,
            ratio_A35R_A33R = dataIn_A35R/dataIn_A33R,
            ratio_B6_B5R = dataIn_B6/dataIn_B5R) %>% 
     dplyr::left_join(ATI_infected, by = "sampleID_meta") %>%
     unique() %>%  
    # add_column(predictionGLMDelta) %>% 
     mutate(#Combined_infected = if_else(is.na(MPXV_delta), ("seronegative"),  Combined_infected),
            ATI_infected = if_else(is.na(MPXV_delta), ("seronegative"),  ATI_infected),
            Final_combined = case_when(is.na(Final_combined) & MPXV_high > 0.75 ~ "MPXV", # New 2023-08-31
                                       is.na(Final_combined) & MVA_high > 0.75 ~ "MVA",
                                       TRUE ~ Final_combined))
  return(dataOutputIgG)
}