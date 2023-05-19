################################################################################
##################### REPEATED MEASURE CORRELATION #############################
################################################################################

# Libraries required:
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rmcorr)
library(ggprism)


# Load data
SEC.df <- read.csv(here("outputs", "11_RNAseq-comparison_II", "11_RNAseq_heatmap-SEC_percent-table.csv"))
UC.df <- read.csv(here("outputs", "11_RNAseq-comparison_II", "11_RNAseq_heatmap-UC_percent-table.csv"))

# Transform data to long format
SEC_long <- SEC.df[, -1] %>%
  pivot_longer(cols = starts_with ("patient"), names_to = "patient", names_prefix = "patient_", values_to = "SEC") %>% 
  mutate(cell.type = recode(cell.type, "epithelial_cancer_cells" = "epithelial")) %>% 
  mutate(id = paste(cell.type, patient, sep = "_"))

UC_long <- UC.df[, -1] %>%
  pivot_longer(cols = starts_with ("patient"), names_to = "patient", names_prefix = "patient_", values_to = "UC") %>% 
  mutate(cell.type = recode(cell.type, "epithelial_cancer_cells" = "epithelial")) %>% 
  mutate(id = paste(cell.type, patient, sep = "_"))

data.long <- inner_join(x = SEC_long, y = UC_long[, 3:4], by = "id")
 

# Repeated measure correlation - SEC vs. UC
# Model
my.rmc.cells <- rmcorr(participant = cell.type , measure1 = SEC, measure2 = UC,
                       dataset = data.long)
# Plot with legend
svg(here("outputs", "11_RNASeq-comparison_II", "11_rmcorr_SECvsUC_legend.svg"))
pdf(here("outputs", "11_RNASeq-comparison_II", "11_rmcorr_SECvsUC_legend.pdf"))
ggplot(data.long, aes(x = SEC , y = UC, color = cell.type , group = cell.type )) + 
  geom_line(aes(x = SEC, y = my.rmc.cells$model$fitted.values), size = 1) + 
  geom_point() + 
  theme_prism() + 
  xlim(0, 100) +
  ylim(0, 100) +
  labs(x = "SEC (%)", y = "UC (%)") +
  annotate("text", x = 60, y = 14, label = paste("r", "=", round(my.rmc.cells$r, digits = 3), sep = " "), hjust = 0) +
  annotate("text", x = 60, y = 8, label = paste("95% CI:", round(my.rmc.cells$CI[1], digits = 3), "-", round(my.rmc.cells$CI[2], digits = 3)), hjust = 0) +
  annotate("text", x = 60, y = 2, label = paste("p = ", format(my.rmc.cells$p, scientific = TRUE, digits = 3)), hjust = 0) 
dev.off() 



