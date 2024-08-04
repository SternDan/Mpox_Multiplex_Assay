rocFunction <- function(data, predictor, response){
  preCol <-
    data %>%
    #   filter(population == "Kupferzell") %>%
    #filter(!is.na(NT_Charite)) %>%
    dplyr::select({{predictor}})
  
  respCol <-
    data %>%
    #   filter(population == "Kupferzell") %>%
    #filter(!is.na(NT_Charite)) %>%
    dplyr::select({{response}})
  
  roc(preCol[[1]], respCol[[1]])
}

rocFunctionELISA <- function(data, predictor, response){
  preCol <-
    data %>%
    #   filter(population == "Kupferzell") %>%
  #  filter(!is.na(NT_Charite)) %>%
    dplyr::select({{predictor}})
  
  respCol <-
    data %>%
    #   filter(population == "Kupferzell") %>%
#    filter(!is.na(NT_Charite)) %>%
    dplyr::select({{response}})
  
  roc(preCol[[1]], respCol[[1]])
}

# Define functions for ROC-Analysis
readRocParam <- function(input){
  dataFrame <- do.call(rbind.data.frame, get(input))# %>%
  # dplyr::select(input = c(`50%`))
  colnames(dataFrame) <- input
  return(dataFrame)
}

transformRocParamFunc <- function(param_list, assay_name){
  paramDataFrame <- do.call(cbind.data.frame, param_list)
  names(paramDataFrame)
  head <- tidyr::fill(as_tibble(names(paramDataFrame)), value)
  headFix <- paste(head$value, c("low", "median", "high"), sep = "-")
  colnames(paramDataFrame) <- headFix
  paramDataFrame_transpose <- t(paramDataFrame)
  antigene <- row.names(paramDataFrame_transpose)
  paramDataFrame_transpose <- as_tibble(paramDataFrame_transpose)
  paramDataFrame_transpose$antigene <- antigene
  paramDataFrame_transpose_sep <-
    paramDataFrame_transpose %>%
    tidyr::separate(antigene, c("antigene", "CI"), sep = "-") %>%
    mutate(assay = assay_name)
}


rocParamFunc <- function(param_list){
  paramDataFrame <- do.call(cbind.data.frame, param_list)
  names(paramDataFrame)
  head <- tidyr::fill(as_tibble(names(paramDataFrame)), value)
  headFix <- paste(head$value, c("low", "median", "high"), sep = "-")
  colnames(paramDataFrame) <- headFix
  return(paramDataFrame)
}


prepDataFrameFunc <- function(input, firstCol, remove){
  input %>%
    dplyr::select({{ firstCol }}, antigene, CI) %>%
    pivot_wider(names_from = c("CI"), values_from = c(firstCol)) %>%
    mutate(antigene = gsub(remove, "", antigene)) %>%
    mutate(antigene = gsub("_ratio", "", antigene)) %>%
    mutate(antigene = gsub("_norm", "", antigene)) %>%
    mutate(antigene = gsub("\\.", "-", antigene)) %>%
    mutate(antigene = gsub("_", " ", antigene)) %>%
    mutate(antigene_class = case_when(grepl("SARS.CoV.2", antigene) ~ "SARS-CoV-2",
                                      grepl("Combined", antigene) ~ "SARS-CoV-2",
                                      grepl("SARS-CoV S1", antigene) ~ "Highly pathogenic",
                                      grepl("MERS-CoV S1", antigene) ~ "Highly pathogenic",
                                      grepl("HCoV", antigene) ~ "Low pathogenic",
                                      grepl("ELISA", antigene) ~ "ELISA", 
                                      grepl("Roche", antigene) ~ "ELISA", 
                                      grepl("Abbott", antigene) ~ "ELISA", 
                                      grepl("DiaSorin", antigene) ~ "ELISA", 
                                      TRUE ~ "Control"))
}

plotSpecFunc <- function(input, ylabel){
  input %>%
    arrange(median) %>%
    mutate(antigene = factor(antigene, levels = antigene)) %>%
    ggplot(mapping = aes(x = (antigene), y = median, ymin = low, ymax = high, color = antigene_class)) +
    geom_pointrange(size = 1) +
    scale_y_continuous(name = ylabel)+
    xlab("Antigene")+
    labs(color = "")+
    coord_flip() +
    scale_color_manual(values=c(wes_palette(n=5, name="Darjeeling1")))+
    theme_bw()
}

rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
          "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")