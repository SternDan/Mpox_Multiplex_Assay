dataInputQuant %>% 
  filter(dil.fitted.final > 1) %>% 
  filter(dilution_assay == 1000) %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = log10(data), y = log10(quant_conc),
                       color = as.factor(batch))) +
  geom_point() +
  facet_wrap("analyte")


dataInputQuant %>% 
  filter(quant_conc > 1) %>% 
  #  filter(dilution_assay == 100) %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = (data.norm), y = log10(quant_conc),
                       color = as.factor(dilution_assay))) +
  geom_point() +
  facet_wrap("analyte")



dataInputQuant %>% 
  filter(quant_conc > 1) %>% 
  filter(dilution_assay == 1000) %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = log10(data), y = data.norm,
                       color = as.factor(panel))) +
  geom_point() +
  facet_wrap("analyte")


dataInputQuant %>% 
  filter(analyte == "Delta") %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  select(sampleID_metadata, assaytype, data.norm) %>% 
  pivot_wider(id_cols = sampleID_metadata, names_from = assaytype, values_from = data.norm, 
              values_fn = mean) %>% 
  ggplot(mapping= aes(x = (ELISA), y = (Multiplex))) +
  geom_point()


dataInputQuant %>% 
  filter(analyte == "Delta") %>% 
  filter(dilution_assay == 1000) %>% 
  filter(isotype == "IgG") %>% 
  select(sampleID_metadata, assaytype, data.norm) %>% 
  pivot_wider(id_cols = sampleID_metadata, names_from = assaytype, values_from = data.norm, 
              values_fn = mean) %>% 
  ggplot(mapping= aes(x = (ELISA), y = (Multiplex))) +
  geom_point()


dataInputQuant %>% 
  filter(analyte == "Delta") %>% 
  filter(dilution_assay == 100) %>% 
  filter(isotype == "IgG") %>% 
  select(sampleID_metadata, assaytype, quant_conc) %>% 
  pivot_wider(id_cols = sampleID_metadata, names_from = assaytype, values_from = quant_conc, 
              values_fn = mean) %>% 
  ggplot(mapping= aes(x = log10(ELISA), y = log10(Multiplex))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  geom_abline(slope = 1)+
  theme_bw()


dataInputQuant %>% 
  mutate(quant_conc = if_else(quant_conc<1, 1, quant_conc)) %>% 
  filter(dilution_assay == 100) %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = log10(data), y = log10(quant_conc),
                       color = as.factor(batch))) +
  geom_pointrange(aes(ymin = log10(quant_conc-quant_sd), 
                      ymax = log10(quant_conc+quant_sd))) +
  facet_wrap("analyte")


dataInputQuant %>% 
  mutate(dil.fitted.final = if_else(dil.fitted.final<1, 1, dil.fitted.final)) %>% 
  filter(dilution_assay == 1000) %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = log10(data), y = log10(dil.fitted.final),
                       color = as.factor(batch))) +
  geom_point() +
  facet_wrap("analyte")

dataInputQuant %>% 
  mutate(quant_single_dil = if_else(quant_single_dil<1, 1, quant_single_dil)) %>% 
  filter(dilution_assay == 1000) %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = log10(data), y = log10(quant_single_dil),
                       color = as.factor(batch))) +
  geom_point() +
  facet_wrap("analyte")


dataInputQuant %>% 
  mutate(quant_conc = if_else(quant_conc<1, 1, quant_conc)) %>% 
  filter(dilution_assay == 100) %>% 
  filter(assaytype == "Multiplex") %>% 
  filter(isotype == "IgG") %>% 
  ggplot(mapping = aes(x = log10(quant_single_dil), y = log10(quant_conc),
                       color = as.factor(source_single_dil))) +
  geom_pointrange(aes(ymin = log10(quant_conc-quant_sd), 
                      ymax = log10(quant_conc+quant_sd))) +
  facet_wrap("analyte")
