####
# Analyse Multiplex-Assay: Coupling controls
# Date final version: 03.08.2024
# Analysis: Daniel Stern
####

#### Remarks ####
# For the comparison with the monospecific antibodies, results from 
# three experimentes are used
# 
# 1) Pox 09
# Due to a mixup of the detection antibodies, only results for mouse mAbs
# can be used from this experiment. mAbs were targeting A27, B5, His-Tag.
# 
# 2) Pox 11
# First panel measurement with human sera. Results for VIG are used from this
# experiment
# 
# 3) Pox 13
# Repetition of coupling control with monospecific polyclonal antibodies
# targeting D8, H3, A33, and L1
####

#### Description of the script ####
# 1) Load packages and clean environment
# 2) Load experimental data for testing of monospecific antibodies
# 3) Load metadata for testing of monospecific antibodies
# 4) Plot data for testing of monospecific antibodies
# 5) Load experimental data for batch comparison
# 6) Load metadata for batch comparison
# 7) Merge data
# 8) Fit standard curves for all antigen using the drlumi-package
# 9) Plot combined standard curves: Color by batch
####

#### Prepare workspace ####
# Clean environment
rm(list = ls(all.names = TRUE))

# Load packages
library(tidyverse)
library(rio)
library(ggthemes)
library(drLumi)
library(ggpubr)
library(drc)

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### Pox 09 ####
#Import data
dataInputCouplingmAbs <- import("input/Pox 09 Ladungskontrol Assay.xlsx",
                                range = c("B10:T27"),
                                col_names = FALSE) %>%
  mutate_at(vars(-`...1`), ~(as.numeric(sub(",", ".", ., fixed = TRUE))))

namesTable <- import("input/Pox 09 Ladungskontrol Assay.xlsx",
                     range = c("C8:T8"),
                     col_names = FALSE) %>%
  pivot_longer(cols = everything()) %>%
  mutate(value = str_replace(value, " \\([^\\)]+\\)", ""), 
         value = str_replace_all(value, " ", "_"),
         value = str_replace_all(value, "-", ".")) %>%
  pull(value)

names(dataInputCouplingmAbs) <- c("well", namesTable)

# Import plate payout mAbs only
plateLayoutmAbs <- import("input/Pox 09_plateLayout.xlsx")

# Combine Data 
dataInputCouplingmAbs <- 
  plateLayoutmAbs %>%
  left_join(dataInputCouplingmAbs, by = c("well"))


#### Pox 11 #### Import data
dataInputCouplingVIG <- import("input/Pox 11 Ladungskontroll Assay.xlsx",
                              range = c("B10:T43"),
                              col_names = FALSE) %>%
  mutate_at(vars(-`...1`), funs(as.numeric(sub(",", ".", ., fixed = TRUE))))

namesTableVIG <- import("input/Pox 11 Ladungskontroll Assay.xlsx",
                        range = c("C8:T8"),
                        col_names = FALSE) %>%
  pivot_longer(cols = everything()) %>%
  mutate(value = str_replace(value, " \\([^\\)]+\\)", ""), 
         value = str_replace_all(value, " ", "_"),
         value = str_replace_all(value, "-", ".")) %>%
  pull(value)

names(dataInputCouplingVIG) <- c("well", namesTableVIG)

# Import plate plate layout
plateLayoutVIG <- import("input/Pox 11_plateLayout.xlsx")

# Combine Data 
dataInputCouplingVIG <- 
  plateLayoutVIG %>%
  left_join(dataInputCouplingVIG, by = c("well"))

#### Pox 13 #### Import data
dataInputCouplingpAbs <- import("input/Pox 13 Ladungskontroll Assay.xlsx",
                                readxl = FALSE,
                                cols = c(2:20),
                                rows = c(10:26),
                                colNames = FALSE) %>%
  mutate_at(vars(-X1), funs(as.numeric(sub(",", ".", ., fixed = TRUE))))

namesTablepAbs <- import("input/Pox 13 Ladungskontroll Assay.xlsx",
                         readxl = FALSE,
                         cols = c(2:20),
                         rows = c(8),
                         colNames = FALSE) %>%
  pivot_longer(cols = everything()) %>%
  mutate(value = str_replace(value, " \\([^\\)]+\\)", ""), 
         value = str_replace_all(value, " ", "_"),
         value = str_replace_all(value, "-", ".")) %>%
  pull(value)

names(dataInputCouplingpAbs) <- c("well", namesTablepAbs)

# Import plate layout
plateLayoutpAbs <- import("input/Pox 13_plateLayout.xlsx")

# Combine data 
dataInputCouplingpAbs <- 
  plateLayoutpAbs %>%
  left_join(dataInputCouplingpAbs, by = c("well"))

#### Combine data and delete non needed data ####
dataInputCoupling <- 
  dataInputCouplingmAbs %>%
  rbind(dataInputCouplingVIG) %>%
  rbind(dataInputCouplingpAbs)

rm(dataInputCouplingmAbs, dataInputCouplingpAbs, dataInputCouplingVIG,
   plateLayoutmAbs, plateLayoutpAbs, plateLayoutVIG,
   namesTable, namesTablepAbs, namesTableVIG)


#### Prepare data frame for plotting ####
dataInputPlotting <- 
  dataInputCoupling %>%
  mutate(plotDil = case_when((type == "mAb" & conc_dil == 1) ~ 1,
                             (type == "mAb" & conc_dil == 0.2) ~ 2,
                             (type == "mAb" & conc_dil == 0.04) ~ 3,
                             (type == "pAb" & conc_dil == 100) ~ 1,
                             (type == "pAb" & conc_dil == 500) ~ 2,
                             (type == "pAb" & conc_dil == 2500) ~ 3,
                             (type == "pAb" & conc_dil == 10) ~ 1,
                             (type == "pAb" & conc_dil == 2) ~ 2,
                             (type == "pAb" & conc_dil == 0.4) ~ 3)) %>%
  filter(!is.na(plotDil)) %>% 
  pivot_longer(c(A29:Hep2), names_to = "Antigen", values_to = "MFI") %>%
  mutate(Antigen = case_when(Antigen == "A27L" ~ "A27",
                             Antigen == "A35R" ~ "A35",
                             Antigen == "A33R" ~ "A33",
                             Antigen == "B5R" ~ "B5",
                             Antigen == "L1R" ~ "L1",
                             Antigen == "D8L" ~ "D8",
                             Antigen == "H3L" ~ "H3",
                             Antigen == "A5L" ~ "A5",
                             Antigen == "ATI.C" ~ "ATI-C",
                             Antigen == "ATI.N" ~ "ATI-N",
                             TRUE ~ Antigen),
         Strain = case_when(Antigen %in% c("VACV", "A27", "A33", "B5", "L1",
                                           "D8", "H3") ~ "VACV",
                            Antigen %in% c("A29", "A35", "B6", "E8", "M1") ~ "MPXV",
                            Antigen %in% c("ATI-C", "ATI-N", "A5") ~ "CPXV",
                            TRUE ~ "Control"),
         Antigen = factor(Antigen, levels = c("A27",
                                              "A29", 
                                              "L1", 
                                              "M1",
                                              "D8",
                                              "E8",
                                              "H3",
                                              "A33", 
                                              "A35",
                                              "B5", 
                                              "B6", 
                                              "A5",
                                              "ATI-C", 
                                              "ATI-N",
                                              "VACV", 
                                              "Hep2"), 
                          ordered = TRUE),
         specificity = if_else(specificity == "Vaccinia", "VIG", specificity),
         specificity = factor(specificity, levels = c("A27", 
                                                      "L1",
                                                      "D8",
                                                      "H3",
                                                      "A33",
                                                      "B5",
                                                      "A5",
                                                      "ATI-C", 
                                                      "ATI-N",
                                                      "VIG",
                                                      "His"), ordered = TRUE),
         Strain = factor(Strain, levels = c("VACV",
                                            "MPXV",
                                            "CPXV",
                                            "Control"),
                         ordered = TRUE)) %>% 
  filter(specificity != "His")



plotCouplingControl <- 
  dataInputPlotting %>%
  ggplot(mapping = aes(x = plotDil, y = MFI, fill = Strain)) + 
  geom_bar(stat = "identity") +
  facet_grid(specificity ~ Antigen) +
  scale_fill_manual(values = colorBlindBlack8[2:8]) +
  scale_x_continuous(name = "Dilution") +
  scale_y_continuous(name = "Multiplex (mfi)") +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        # axis.title.x = element_blank(),
        #  axis.text.x = element_blank(),
        #  axis.line.x = element_blank(),
        plot.title = element_text(size = 10))

ggsave("output/plotCoupling.png", plotCouplingControl, width = 12, height = 8, dpi = 600)


#### Load experimental data for batch to batch comparison ####
dataInputPox17 <- import("input/Pox 17 Ladungskontroll Assay.xlsx",
                         readxl = FALSE,
                         cols = c(2:20),
                         rows = c(10:43),
                         colNames = FALSE) %>%
  mutate_at(vars(-X1), ~(as.numeric(sub(",", ".", ., fixed = TRUE))))

namesTablePox17 <- import("input/Pox 17 Ladungskontroll Assay.xlsx",
                          readxl = FALSE,
                          cols = c(2:20),
                          rows = c(8),
                          colNames = FALSE) %>%
  pivot_longer(cols = everything()) %>%
  mutate(value = str_replace(value, " \\([^\\)]+\\)", ""), 
         value = str_replace_all(value, " ", "_")) %>%
  pull(value)

names(dataInputPox17) <- c("well", namesTablePox17)


dataInputPox24 <- import("input/Pox 24 Kopplungskontrolle.xlsx",
                         readxl = FALSE,
                         cols = c(2:22),
                         rows = c(10:43),
                         colNames = FALSE) %>%
  mutate_at(vars(-X1), ~(as.numeric(sub(",", ".", ., fixed = TRUE))))

namesTablePox24 <- import("input/Pox 24 Kopplungskontrolle.xlsx",
                          readxl = FALSE,
                          cols = c(2:22),
                          rows = c(8),
                          colNames = FALSE) %>%
  pivot_longer(cols = everything()) %>%
  mutate(value = str_replace(value, " \\([^\\)]+\\)", ""), 
         value = str_replace_all(value, " ", "_")) %>%
  pull(value)

names(dataInputPox24) <- c("well", namesTablePox24)


dataInputPox24B <- import("input/Pox 24 B Kopplungskontrolle.xlsx",
                          readxl = FALSE,
                          cols = c(2:22),
                          rows = c(10:34),
                          colNames = FALSE) %>%
  mutate_at(vars(-X1), ~(as.numeric(sub(",", ".", ., fixed = TRUE))))

namesTablePox24B <- import("input/Pox 24 B Kopplungskontrolle.xlsx",
                           readxl = FALSE,
                           cols = c(2:22),
                           rows = c(8),
                           colNames = FALSE) %>%
  pivot_longer(cols = everything()) %>%
  mutate(value = str_replace(value, " \\([^\\)]+\\)", ""), 
         value = str_replace_all(value, " ", "_")) %>%
  pull(value)

names(dataInputPox24B) <- c("well", namesTablePox24B)


dataInputPox27 <- import("input/Pox 27 Ladungskontrolle 160123.xlsx",
                         readxl = FALSE,
                         cols = c(2:22),
                         rows = c(10:43),
                         colNames = FALSE) %>%
  mutate_at(vars(-X1), ~(as.numeric(sub(",", ".", ., fixed = TRUE))))

namesTablePox27 <- import("input/Pox 27 Ladungskontrolle 160123.xlsx",
                          readxl = FALSE,
                          cols = c(2:22),
                          rows = c(8),
                          colNames = FALSE) %>%
  pivot_longer(cols = everything()) %>%
  mutate(value = str_replace(value, " \\([^\\)]+\\)", ""), 
         value = str_replace_all(value, " ", "_")) %>%
  pull(value)

names(dataInputPox27) <- c("well", namesTablePox27)


#### Load metadata for batch to batch comparison ####
metadataInputPox17 <- import("input/Pox 17_plateLayout.xlsx")
metadataInputPox24 <- import("input/Pox 24_plateLayout.xlsx")
metadataInputPox24B <- import("input/Pox 24 B_plateLayout.xlsx")
metadataInputPox27 <- import("input/Pox 27_plateLayout.xlsx")

dataInputPox17combined <-
  dataInputPox17 %>% 
  left_join(metadataInputPox17, by = "well") %>% 
  dplyr::select(-ahIgG, -HSA) %>% 
  mutate(assay = "Pox17")

dataInputPox24combined <-
  dataInputPox24 %>% 
  left_join(metadataInputPox24, by = "well") %>% 
  dplyr::select(-ahIgG, -HSA, -ah_IgM, -ah_IgA) %>% 
  mutate(assay = "Pox24")

dataInputPox24Bcombined <-
  dataInputPox24B %>% 
  left_join(metadataInputPox24B, by = "well") %>% 
  mutate(assay = "Pox24B")

dataInputPox27combined <-
  dataInputPox27 %>% 
  left_join(metadataInputPox27, by = "well") %>% 
  dplyr::select(-ahIgG, -HSA, -IgM, -IgA) %>% 
  mutate(assay = "Pox27")

## Delete non-needed columns ##
rm(dataInputPox17, dataInputPox24B, dataInputPox24, dataInputPox27,
   metadataInputPox17, metadataInputPox24B, metadataInputPox24, metadataInputPox27)

dataInputPox17long <-
  dataInputPox17combined %>% 
  pivot_longer(cols = c(A29:Hep2), names_to = "antigen",
               values_to = "mfi")

dataInputPox24long <-
  dataInputPox24combined %>% 
  pivot_longer(cols = c(A29:Hep2), names_to = "antigen",
               values_to = "mfi")

dataInputPox24Blong <-
  dataInputPox24Bcombined %>% 
  pivot_longer(cols = c(E8), names_to = "antigen",
               values_to = "mfi")

dataInputPox27long <-
  dataInputPox27combined %>% 
  pivot_longer(cols = c(A29:Hep2), names_to = "antigen",
               values_to = "mfi")


dataInputBatch <-
  rbind(dataInputPox17long,
        dataInputPox24long,
        # dataInputPox24Blong,
        dataInputPox27long) %>% 
  filter(!is.na(mfi)) %>% 
  filter(!(assay == "Pox24" & batch == 3 & antigen == "E8")) %>% # Exclude bad coupling
  mutate(batch = case_when(batch == 1 ~ "Batch_09/22",
                           batch == 2 ~ "Batch_10/22",
                           batch == 3 ~ "Batch_12/22",
                           batch == 4 ~ "Batch_01/23"))
                           
rm(dataInputPox17combined,
   dataInputPox24combined,
   dataInputPox24Bcombined,
   dataInputPox27combined,
   dataInputPox17long,
   dataInputPox24long,
   dataInputPox24Blong,
   dataInputPox27long)


meandataInputBatch <-
  dataInputBatch %>% 
  group_by(antibody, antigen, concentration) %>% 
  summarize(meanMFI = mean(mfi),
            sdMFI = sd(mfi),
            nMFI = length(unique(mfi)))

meandataInputBatch %>% 
  ggplot(mapping = aes(x = log10(concentration), y = (meanMFI))) +
  geom_pointrange(aes(ymin = (meanMFI-sdMFI), ymax = (meanMFI+sdMFI))) +
  facet_wrap("antigen")

plotBatch <-
dataInputBatch %>% 
  mutate(antigen = factor(antigen, levels = c("A27L",
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
                                              "Hep2"), 
                          ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = log10(concentration), y = (mfi), color = as.factor(batch))) +
  geom_point() +
  geom_smooth(method = drm, 
              method.args = list(fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))), se = F) +
 # geom_smooth(color = "#0072bc", alpha = 0.2, method = drm, 
#              method.args = list(fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))), se = F) +
  facet_wrap("antigen") +
  scale_color_manual(name = "Batch", values = colorBlindBlack8[2:8]) +
  scale_x_continuous(name = "VIG (log10 µg/mL)") +
  scale_y_continuous(name = "Multiplex (mfi)") +
  theme_bw() +
  theme(legend.position = "top",
        # panel.grid = element_blank(),
        strip.background = element_blank())

plotCoupling <- 
ggarrange(plotCouplingControl, plotBatch, nrow = 2, align = "hv",
          labels = "auto")

ggsave("output/plotCoupling.png", width = 10, height = 12, dpi = 300)
