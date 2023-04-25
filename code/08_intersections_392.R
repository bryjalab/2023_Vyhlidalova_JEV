################################################################################
######################## 9/11 PATIENTS PROTEINS ################################
################################################################################

# NOTE
# 392 proteins were identified in at least 9 of 11 patients
# less restrictive filtering than 11/11 (intersection of intersections)

# Create output directory
dir.create(here('outputs', '08_intersections_392'))

# Libraries required:
library(here)
library(dplyr)
library(ggplot2)
library(gprofiler2)

# Load the data
load(here('outputs', '06_methods-intersection', '06_intersections.Rdata'))
load(here("outputs", "05_filter-out-B-samples", "05_data.filtered.Rdata"))
patients <- read.csv(here('outputs', '06_methods-intersection', '06_upsetplot_background-table.csv'))
patients$X <- NULL

# Filter the proteins which are in >= 9/11 patients
data_392 <- patients %>%
  rowwise() %>%
  mutate(sum_patients = sum(c_across(patient_1: patient_11))) %>%
  filter(sum_patients > 8) %>%
  select(!sum_patients)

# Grep the 392 proteins from filtered table
data_392 <- data.filtered[data.filtered$Suggested.Symbol %in% data_392$value, ]

############################# Gene ontology ####################################

gostres_392 <- gost(query = data_392$Suggested.Symbol,
                    organism = "hsapiens", ordered_query = FALSE,
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                    measure_underrepresentation = FALSE, evcodes = TRUE,
                    user_threshold = 0.05, correction_method = "g_SCS",
                    domain_scope = "annotated", custom_bg = NULL,
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)

GO_392 <- gostres_392$result

GO_392_csv <- apply(GO_392,2,as.character) # otherwise an error occurs
write.csv(GO_392_csv, here('outputs', '08_intersections_392', '08_392-proteins_gprofiler.csv'))

svg(here('outputs', '08_intersections_392', '08_392-proteins_GO-BP.svg'))
#pdf(here('outputs', '08_intersections_392', '08_392-proteins_GO-BP.pdf'))
GO_392 %>%
  filter(source == "GO:BP") %>%
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  ggplot(aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
  geom_bar(stat = "identity",  fill = "#CBE11E") +
  coord_flip() +
  labs(title = "GO: Biological process, top 10",
       x = "term",
       y = "-log10(FDR)") +
  theme_minimal()
dev.off()

svg(here('outputs', '08_intersections_392', '08_392-proteins_GO-CC.svg'))
#pdf(here('outputs', '08_intersections_392', '08_392-proteins_GO-CC.pdf'))
GO_392 %>%
  filter(source == "GO:CC") %>%
  arrange(p_value) %>%
  slice_head(n = 10) %>%
  ggplot(aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
  geom_bar(stat = "identity",  fill = "#CBE11E") +
  coord_flip() +
  labs(title = "GO: Cellular component, top 10",
       x = "term",
       y = "-log10(FDR)") +
  theme_minimal()
dev.off()

############################# Human Cell Map ###################################

HCM <- read.csv(here('outputs', '02_data-cleaning', 'updated_preys-latest.csv'))

HCM <- left_join(data_392, HCM %>% select(Suggested.Symbol, MMF.localization), by = c("Suggested.Symbol" = "Suggested.Symbol"))

svg(here('outputs', '08_intersections_392', '08_392-proteins_HCM.svg'))
pdf(here('outputs', '08_intersections_392', '08_392-proteins_HCM.pdf'))
HCM %>%
  group_by(MMF.localization) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = reorder(MMF.localization, n), y = n))+
  geom_bar(stat= "identity")+
  coord_flip() +
  theme_classic()
dev.off()

svg(here('outputs', '08_intersections_392', '08_392-proteins_HCM_stacked.svg'))
#pdf(here('outputs', '08_intersections_392', '08_392-proteins_HCM_stacked.pdf'))
HCM %>%
  group_by(MMF.localization) %>%
  summarise(n = n()) %>%
  mutate(x = '392_proteins') %>%
  filter(!MMF.localization == "NA") %>%
  mutate(compartment = factor(MMF.localization, levels = c("plasma membrane", "cell junction", "miscellaneous",
                                                           "nuclear outer membrane-ER membrane network",
                                                           "mitochondrial outer membrane, peroxisome",
                                                           "endosome, lysosome",
                                                           "early endosome, recycling endosome",
                                                           "actin cytoskeleton, cytosol", "mitochondrial matrix",
                                                           "Golgi apparatus", "ER membrane",
                                                           "ER lumen",
                                                           "mitochondrial inner membrane, mitochondrial intermembrane space",
                                                           "chromatin", "nucleoplasm", "nucleolus",
                                                           "cytoplasmic ribonucleoprotein granule" ))) %>%
  ggplot(aes(x = x, y = n , fill= compartment))+
  geom_bar(stat= "identity", position = "stack")+
  theme_classic()
dev.off()

# Export the csv table
data.export.392 <- patients[patients$value %in% data_392$Suggested.Symbol, ]
write.csv(data.export.392, here('outputs', '08_intersections_392', '08_intersections_392.csv'))
