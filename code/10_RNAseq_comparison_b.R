################################################################################
##################### COMPARISON WITH IZAR II. #################################
################################################################################

# Libraries required:
library(here)
library(dplyr)
library(ComplexHeatmap)
library(tibble)
library(tidyr)
library(ggplot2)

# Create output directory
dir.create(here('outputs', '11_RNAseq-comparison_II'))

# Load the data:
RNAseq_markers <- read.csv(here('outputs', '02_data-cleaning', 'updated_Izar_cell_markers_updated_cropped.csv'))
RNAseq_markers$X <- NULL
RNAseq_markers$cell.type[RNAseq_markers$cell.type == "CAFs"] <- "fibroblast"
patients <- read.csv(here('outputs', '06_methods-intersection', '06_upsetplot_background-table.csv')) # 2418 proteins present in >=1 S&U
patients$X <- NULL
data.filtered <- read.csv(here('outputs', '05_filter-out-B-samples', '05_data_filtered-out-B-samples.csv'))

# Plot binary heatmap for Izar markers present in our data
RNAseq_subset_binary <- left_join(patients, RNAseq_markers, by = c("value" = "Suggested.Symbol"))
RNAseq_subset_binary <- na.omit(RNAseq_subset_binary)

RNAseq_subset_binary_mat <- as.matrix(RNAseq_subset_binary[, 2:12])
rownames(RNAseq_subset_binary_mat) <-RNAseq_subset_binary$value

category <- factor(RNAseq_subset_binary$cell.type)

row_ha <- rowAnnotation(category = category)

svg(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-binary.svg'))
#pdf(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq-heatmap-binary.pdf'))
Heatmap(RNAseq_subset_binary_mat,  col = c("white", "black"),
        right_annotation = row_ha,
        row_names_gp = gpar(fontsize = 5))
dev.off()

# SEC barplots
SEC1 <- RNAseq_subset_binary_mat
SEC2 <- data.filtered[data.filtered$Suggested.Symbol %in% rownames(RNAseq_subset_binary_mat), ]

SEC2 <- SEC2 %>% select(-starts_with("U"), -X, -name, -ID)
rownames(SEC2) <- SEC2$Suggested.Symbol
SEC2$Suggested.Symbol <- NULL

SEC1 <- SEC1[order(row.names(SEC1)), ]
SEC2 <- SEC2[order(row.names(SEC2)), ]

SEC_final <- SEC1

for (i in 1:nrow(SEC_final)) {
  for (j in 1: ncol(SEC_final)){
    if (SEC_final[i,j] == 1) {
      SEC_final[i,j] <- SEC2[i,j]
    }
  }
}

SEC_long <- as.data.frame(SEC_final)
SEC_long <- left_join(rownames_to_column(SEC_long), RNAseq_subset_binary %>% select(value, cell.type), by = c("rowname" = "value"))
SEC_long <- pivot_longer(SEC_long, cols = patient_1:patient_11, names_to = 'patient', values_to = 'intensity')

svg(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-SEC.svg'))
#pdf(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq-heatmap-SEC.pdf'))
SEC_long %>%
  group_by(patient) %>%
  mutate(intensity_norm = intensity/sum(intensity)) %>%
  group_by(patient, cell.type) %>%
  summarise(sum_category = sum(intensity_norm)) %>%
  mutate(patients_factor = factor(patient, levels = c("patient_1", "patient_2", "patient_3", "patient_4",
                                                      "patient_5", "patient_6", "patient_7", "patient_8",
                                                      "patient_9", "patient_10", "patient_11"))) %>%
  ggplot(., aes(x = patients_factor, y = sum_category, fill = cell.type))+
  geom_bar(stat = "identity")  +
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  labs(title = "RNAseq SEC")
dev.off()

# UC barplots

UC1 <- RNAseq_subset_binary_mat
UC2 <- data.filtered[data.filtered$Suggested.Symbol %in% rownames(RNAseq_subset_binary_mat), ]

UC2 <- UC2 %>% select(-c(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11), -X, -name, -ID)
rownames(UC2) <- UC2$Suggested.Symbol
UC2$Suggested.Symbol <- NULL

UC1 <- UC1[order(row.names(UC1)), ]
UC2 <- UC2[order(row.names(UC2)), ]

UC_final <- UC1

for (i in 1:nrow(UC_final)) {
  for (j in 1: ncol(UC_final)){
    if (UC_final[i,j] == 1) {
      UC_final[i,j] <- UC2[i,j]
    }
  }
}

UC_long <- as.data.frame(UC_final)
UC_long <- left_join(rownames_to_column(UC_long), RNAseq_subset_binary %>% select(value, cell.type), by = c("rowname" = "value"))
UC_long <- pivot_longer(UC_long, cols = patient_1:patient_11, names_to = 'patient', values_to = 'intensity')

svg(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-UC.svg'))
pdf(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq-heatmap-UC.pdf'))
UC_long %>%
  group_by(patient) %>%
  mutate(intensity_norm = intensity/sum(intensity)) %>%
  group_by(patient, cell.type) %>%
  summarise(sum_category = sum(intensity_norm)) %>%
  mutate(patients_factor = factor(patient, levels = c("patient_1", "patient_2", "patient_3", "patient_4",
                                                      "patient_5", "patient_6", "patient_7", "patient_8",
                                                      "patient_9", "patient_10", "patient_11"))) %>%
  ggplot(., aes(x = patients_factor, y = sum_category, fill = cell.type))+
  geom_bar(stat = "identity")  +
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  labs(title = "RNAseq UC")
dev.off()


# Create tables with percentages

SEC.df <- SEC_long %>%
  group_by(patient) %>%
  mutate(intensity_norm = intensity/sum(intensity)*100) %>%
  group_by(patient, cell.type) %>%
  summarise(sum_category = sum(intensity_norm)) %>%
  pivot_wider(names_from = "patient", values_from = "sum_category")

write.csv(SEC.df, here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-SEC_percent-table.csv'))

UC.df <- UC_long %>%
  group_by(patient) %>%
  mutate(intensity_norm = intensity/sum(intensity)*100) %>%
  group_by(patient, cell.type) %>%
  summarise(sum_category = sum(intensity_norm)) %>%
  pivot_wider(names_from = "patient", values_from = "sum_category")

write.csv(UC.df, here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-UC_percent-table.csv'))

