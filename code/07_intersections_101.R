################################################################################
############### INTERSECTION of INTERSECTIONS = 101 proteins ###################
################################################################################

# NOTE
# 101 proteins are in the 'intersection of intersections', ie were identified by both
# UC and SEC in all patients

# Libraries required:
library(here)
library(dplyr)
library(ggplot2)
library(gprofiler2)

# Create output directory
dir.create(here('outputs', '07_intersections_101'))

# Load the data
load(here('outputs', '06_methods-intersection', '06_intersections.Rdata'))
load(here("outputs", "05_filter-out-B-samples", "05_data.filtered.Rdata"))

# Grep the 101 proteins from filtered table
data_101 <- data.filtered[data.filtered$Suggested.Symbol %in% intersection_of_intersections, ]

############################# Gene ontology ####################################
# actual version of database: e109_eg56_p17_1d3191d, database updated on 29/03/2023 

gostres_101 <- gost(query = data_101$Suggested.Symbol,
                    organism = "hsapiens", ordered_query = FALSE,
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                    measure_underrepresentation = FALSE, evcodes = TRUE,
                    user_threshold = 0.05, correction_method = "g_SCS",
                    domain_scope = "annotated", custom_bg = NULL,
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)

GO_101 <- gostres_101$result

GO_101_csv <- apply(GO_101,2,as.character) # otherwise an error occurs
write.csv(GO_101_csv, here('outputs', '07_intersections_101', '07_101-proteins_gprofiler.csv'))

svg(here('outputs', '07_intersections_101', '07_101-proteins_GO-BP.svg'))
#pdf(here('outputs', '07_intersections_101', '07_101-proteins_GO-BP.pdf'))
GO_101 %>%
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

svg(here('outputs', '07_intersections_101', '07_101-proteins_GO-CC.svg'))
#pdf(here('outputs', '07_intersections_101', '07_101-proteins_GO-CC.pdf'))
GO_101 %>%
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

HCM <- left_join(data_101, HCM %>% select(Suggested.Symbol, MMF.localization), by = c("Suggested.Symbol" = "Suggested.Symbol"))

svg(here('outputs', '07_intersections_101', '07_101-proteins_HCM.svg'))
#pdf(here('outputs', '07_intersections_101', '07_101-proteins_HCM.pdf'))
HCM %>%
  group_by(MMF.localization) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = reorder(MMF.localization, n), y = n))+
  geom_bar(stat= "identity", position = "stack")+
  coord_flip() +
  theme_classic()
dev.off()

svg(here('outputs', '07_intersections_101', '07_101-proteins_HCM_stacked.svg'))
#pdf(here('outputs', '07_intersections_101', '07_101-proteins_HCM_stacked.pdf'))
HCM %>%
  group_by(MMF.localization) %>%
  summarise(n = n()) %>%
  mutate(x = '101_proteins') %>%
  filter(!MMF.localization == "NA") %>%
  mutate(compartment = factor(MMF.localization, levels = c("plasma membrane", "cell junction", "miscellaneous",
                                                           "nuclear outer membrane-ER membrane network",
                                                           "mitochondrial outer membrane, peroxisome",
                                                           "endosome, lysosome",
                                                           "early endosome, recycling endosome",
                                                           "actin cytoskeleton, cytosol", "mitochondrial matrix",
                                                           "Golgi apparatus", "ER membrane"))) %>%
  ggplot(aes(x = x, y = n , fill= compartment))+
  geom_bar(stat= "identity", position = "stack")+
  theme_classic()
dev.off()

# Export the csv table
write.csv(data_101, here('outputs', '07_intersections_101', '07_intersections_101.csv'))
