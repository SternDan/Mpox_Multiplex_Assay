####
# Analyse Multiplex-Assay: Coupling controls
# Date: 213.03.2025
# Analysis: Daniel Stern
####


#### Prepare workspace ####
# Clean environment
rm(list = ls(all.names = TRUE))

# Load packages
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(drc)

####
# Load data
load("input/dataInputBatch.Rdata")
load("input/dataInputPlotting.Rdata")


plotCouplingControl <- 
  dataInputPlotting %>%
  ggplot(mapping = aes(x = plotDil, y = MFI, fill = Strain)) + 
  geom_bar(stat = "identity") +
  facet_grid(specificity ~ Antigen) +
  scale_fill_manual(values = colorblind_pal()(8)[2:8]) +
  scale_x_continuous(name = "Dilution") +
  scale_y_continuous(name = "Multiplex (mfi)") +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size = 10))

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
  facet_wrap("antigen") +
  scale_color_manual(name = "Batch", values = colorblind_pal()(8)[2:8]) +
  scale_x_continuous(name = "VIG (log10 µg/mL)") +
  scale_y_continuous(name = "Multiplex (mfi)") +
  theme_bw() +
  theme(legend.position = "top",
        strip.background = element_blank())

plotCoupling <- 
  ggarrange(plotCouplingControl, plotBatch, nrow = 2, align = "hv",
            labels = "auto")

ggsave("output/plotCoupling.png", width = 10, height = 12, dpi = 300)
