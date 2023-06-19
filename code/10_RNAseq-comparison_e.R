################################################################################
##################### REPEATED MEASURE CORRELATION #############################
################################################################################

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
data.long$patient <- as.factor(data.long$patient)

primary.cells <- data.frame(patient = c(1:10),
                            epithelial = c(62.13,	31.30,	90.34,	7.82,	98.85,	86.96,	3.79,	85.85,
                                           43.00,	17.56),
                            fibroblast = c(2.29,	6.15,	1.21,	6.41,	0.88,	11.88,	2.17,	3.88,	1.59,	2.67),
                            macrophages = c(35.58, 62.55, 8.45, 85.77, 0.27,	1.16,	94.03,	10.28,55.40,	79.77
                            ))

primary.cells.long <- pivot_longer(primary.cells, cols = epithelial:macrophages, names_to = "cell.type", values_to = "FC")
#primary.cells.long$patient <- paste("patient_", primary.cells.long$patient, collapse = NULL, sep = "")
primary.cells.long$patient <- as.factor(primary.cells.long$patient)

# Merge flow cytometry data and SEC/UC data
data.long <- left_join(data.long, primary.cells.long)
data.long$cell.type <- as.factor(data.long$cell.type)

data.long <- na.omit(data.long)

# Repeated measure correlation - SEC vs. FC
# Model
my.rmc.cells.SEC <- rmcorr(participant = cell.type, measure1 = FC, measure2 = SEC,
                           dataset = data.long)
# Plot
svg(here("outputs", "11_RNASeq-comparison_II", "11_rmcorr_SECvsFC.svg"))
pdf(here("outputs", "11_RNASeq-comparison_II", "11_rmcorr_SECvsFC.pdf"))
ggplot(data.long, aes(x = FC, y = SEC, color = cell.type, group = cell.type)) + 
  geom_line(aes(x = FC, y = my.rmc.cells.SEC$model$fitted.values), size = 1) + 
  geom_point() + 
  theme_prism() + 
  xlim(0,100) +
  ylim(0, 100) +
  labs(x = "FC (%)", y = "SEC (%)", title = "Repeated measure correlation: SEC") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 60, y = 14, label = paste("r", "=", round(my.rmc.cells.SEC$r, digits = 3), sep = " "), hjust = 0) +
  annotate("text", x = 60, y = 8, label = paste("95% CI:", round(my.rmc.cells.SEC$CI[1], digits = 3), "-", round(my.rmc.cells.SEC$CI[2], digits = 3)), hjust = 0) +
  annotate("text", x = 60, y = 2, label = paste("p = ", format(my.rmc.cells.SEC$p, scientific = TRUE, digits = 3)), hjust = 0) 
dev.off()

# Repeated measure correlation - UC vs. FC
# Model
my.rmc.cells.UC <- rmcorr(participant = cell.type, measure1 = FC, measure2 = UC,
                          dataset = data.long)
# Plot
svg(here("outputs", "11_RNASeq-comparison_II", "11_rmcorr_UCvsFC.svg"))
pdf(here("outputs", "11_RNASeq-comparison_II", "11_rmcorr_UCvsFC.pdf"))
ggplot(data.long, aes(x = FC, y = UC, color = cell.type, group = cell.type)) + 
  geom_line(aes(x = FC, y = my.rmc.cells.UC$model$fitted.values), size = 1) + 
  geom_point() + 
  theme_prism() + 
  xlim(0,100) +
  ylim(0, 100) +
  labs(x = "FC (%)", y = "UC (%)", title = "Repeated measure correlation: UC") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 60, y = 14, label = paste("r", "=", round(my.rmc.cells.UC$r, digits = 3), sep = " "), hjust = 0) +
  annotate("text", x = 60, y = 8, label = paste("95% CI:", round(my.rmc.cells.UC$CI[1], digits = 3), "-", round(my.rmc.cells.UC$CI[2], digits = 3)), hjust = 0) +
  annotate("text", x = 60, y = 2, label = paste("p = ", format(my.rmc.cells.UC$p, scientific = TRUE, digits = 3)), hjust = 0) 
dev.off()
