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
RNAseq_markers <- read.csv(here('outputs', '02_data-cleaning', 'updated_Izar_cell_markers_updated.csv'))

# Filter just the subset of proteins defined by Izar:
RNAseq_markers <- RNAseq_markers[RNAseq_markers$Gene.Name %in% c("CLDN3", "FOLR1", "ELF3", "CLDN4", "EPCAM", "TACSTD2", "KRTCAP3", "MMP7", "CLDN7", "SOX17",
                                                                 "PRSS22", "C1orf186", "TSPAN1", "CKB", "S100A14", "LCN2", "SMIM22", "KLF5", "TMEM139",
                                                                 "FXYD3", "PRSS8", "IGF2", "LYNX1", "WNT7A", "SPINT2", "KCNK15", "LYPD1", "SCNN1A", "MEDAG",
                                                                 "FGF7", "COL1A1", "DCN", "SERPINE1", "COL1A2", "VCAM1", "EMILIN1", "NID2", "BDKRB1", "COL3A1",
                                                                 "RGS4", "CDH11", "POSTN", "GPC3", "COL5A1", "CCDC80", "COL5A2", "LOX", "PROCR", "CALB2", "SERPINB2",
                                                                 "TDO2", "FILIP1L", "COL8A1", "PTGER3", "DPP4", "PLS3", "AIF1", "MS4A4A", "MS4A6A", "LY86", "LRRC25",
                                                                 "PILRA", "STAB1", "VSIG4", "CD14", "CD68", "MS4A7", "SLCO2B1", "LILRB2", "CYBB", "ASGR1", "C5AR1",
                                                                 "CD84", "CSF1R", "FCGR1A", "CD163", "ADAP2", "C1QC", "IGSF6", "SLC11A1", "SCIMP", "CD33",
                                                                 "RP11-290F20.3", "FPR1"), ]

RNAseq_markers$X <- NULL
RNAseq_markers$cell.type[RNAseq_markers$cell.type == "CAFs"] <- "fibroblast"


# Heatmaps ####
## SEC ####
# Load MS intensities and reduce them to binary form according to their presence/non-presence
data.filtered <- read.csv(here('outputs', '05_filter-out-B-samples', '05_data_filtered-out-B-samples.csv'))

data.filtered.binary <- data.filtered
data.filtered.binary[is.na(data.filtered.binary)] <- 0
data.filtered.binary <- data.filtered.binary %>%
  mutate(across(U1:S11, ~ ifelse(. > 0, 1, 0)))


# Select 41 cell-specific proteins present in >=1 S&U and add their presence info specific to SEC and UC
RNAseq_subset_binary <- left_join(patients, RNAseq_markers, by = c("value" = "Suggested.Symbol"))
RNAseq_subset_binary <- na.omit(RNAseq_subset_binary) %>% 
  mutate(cell.type = factor(cell.type, levels = c("epithelial_cancer_cells", "fibroblast", "macrophages"))) 
levels(RNAseq_subset_binary$cell.type) <- list("malignant" = "epithelial_cancer_cells", "fibroblast" = "fibroblast", "macrophage" = "macrophages")


RNAseq_subset_binary <- data.filtered.binary %>% 
  select(c(starts_with("S"), starts_with("U"))) %>% 
  left_join(x = RNAseq_subset_binary, y = ., by = c("value" = "Suggested.Symbol"))

# Plot binary heatmap for Izar markers present in our data specific to SEC
RNAseq_subset_binary_SEC <- RNAseq_subset_binary %>% 
  select(c("value", "cell.type", starts_with("S"))) 

RNAseq_subset_binary_SEC_mat <- as.matrix(RNAseq_subset_binary_SEC[, 3:13])
rownames(RNAseq_subset_binary_SEC_mat) <-RNAseq_subset_binary$value

category <- RNAseq_subset_binary_SEC$cell.type

row_ha <- rowAnnotation(category = category)

svg(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-binary_SEC.svg'))
pdf(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq-heatmap-binary_SEC.pdf'))
Heatmap(RNAseq_subset_binary_SEC_mat,  col = c("white", "black"),
        right_annotation = row_ha,
        split = RNAseq_subset_binary_SEC$cell.type,
        row_title = NULL,
        cluster_row_slices = FALSE,
        row_names_gp = gpar(fontsize = 5))
dev.off()

# Plot binary heatmap for Izar markers present in our data specific to UC
RNAseq_subset_binary_UC <- RNAseq_subset_binary %>% 
  select(c("value", "cell.type", starts_with("U"))) 

RNAseq_subset_binary_UC_mat <- as.matrix(RNAseq_subset_binary_UC[, 3:13])
rownames(RNAseq_subset_binary_UC_mat) <-RNAseq_subset_binary$value

category <- RNAseq_subset_binary_SEC$cell.type

row_ha <- rowAnnotation(category = category)

svg(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-binary_UC.svg'))
pdf(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq-heatmap-binary_UC.pdf'))
Heatmap(RNAseq_subset_binary_UC_mat,  col = c("white", "black"),
        right_annotation = row_ha,
        split = RNAseq_subset_binary_UC$cell.type,
        cluster_row_slices = FALSE,
        row_title = NULL,
        row_names_gp = gpar(fontsize = 5))
dev.off()


# SEC barplots
SEC1 <- RNAseq_subset_binary_SEC_mat
SEC2 <- data.filtered[data.filtered$Suggested.Symbol %in% rownames(RNAseq_subset_binary_SEC_mat), ]

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
SEC_long <- pivot_longer(SEC_long, cols = S1:S11, names_to = 'patient', values_to = 'intensity')
SEC_long$patient <- gsub(pattern = "S", replacement = "patient_", SEC_long$patient)

svg(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-SEC.svg'))
#pdf(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq-heatmap-SEC.pdf'))
SEC_long %>%
  group_by(patient, cell.type) %>%
  mutate(mean.intensity = mean(intensity)) %>%
  ungroup() %>%
  group_by(patient) %>%
  mutate(intensity.relative = mean.intensity/sum(intensity)) %>%
  mutate(patients_factor = factor(patient, levels = c("patient_1", "patient_2", "patient_3", "patient_4",
                                                      "patient_5", "patient_6", "patient_7", "patient_8",
                                                      "patient_9", "patient_10", "patient_11"))) %>%
  ggplot(., aes(x = patients_factor, y = intensity.relative, fill = cell.type))+
  geom_bar(stat = "identity")  +
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  labs(title = "RNAseq SEC")
dev.off()

# UC barplots

UC1 <- RNAseq_subset_binary_UC_mat
UC2 <- data.filtered[data.filtered$Suggested.Symbol %in% rownames(RNAseq_subset_binary_UC_mat), ]

UC2 <- UC2 %>% select(c("Suggested.Symbol", starts_with("U")))
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
UC_long <- pivot_longer(UC_long, cols = U1:U11, names_to = 'patient', values_to = 'intensity')
UC_long$patient <- gsub(pattern = "U", replacement = "patient_", UC_long$patient)

# Save SEC and UC long dataframes
save(SEC_long, UC_long, file = here('outputs', '11_RNAseq-comparison_II', '11_SEC-UC-long-df.RData'))

svg(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-UC.svg'))
#pdf(here('outputs', '11_RNAseq-comparison_II', '11_RNAseq-heatmap-UC.pdf'))
UC_long %>%
  group_by(patient, cell.type) %>%
  mutate(mean.intensity = mean(intensity)) %>%
  ungroup() %>%
  group_by(patient) %>%
  mutate(intensity.relative = mean.intensity/sum(intensity)) %>%
  mutate(patients_factor = factor(patient, levels = c("patient_1", "patient_2", "patient_3", "patient_4",
                                                      "patient_5", "patient_6", "patient_7", "patient_8",
                                                      "patient_9", "patient_10", "patient_11"))) %>%
  ggplot(., aes(x = patients_factor, y = intensity.relative, fill = cell.type))+
  geom_bar(stat = "identity")  +
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  labs(title = "RNAseq UC")
dev.off()


# Create tables with percentages

SEC.df <- SEC_long %>%
  group_by(patient, cell.type) %>%
  mutate(mean.intensity = mean(intensity)) %>%
  ungroup() %>%
  group_by(patient) %>%
  mutate(intensity.percent = mean.intensity/sum(intensity)*100) %>%
  group_by(patient, cell.type) %>%
  summarise(sum.percent = sum(intensity.percent)) %>%
  pivot_wider(names_from = "patient", values_from = "sum.percent")

SEC.df <- SEC.df[, c(1,2, 5:12, 3,4)]

write.csv(SEC.df, here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-SEC_percent-table.csv'))

UC.df <- UC_long %>%
  group_by(patient, cell.type) %>%
  mutate(mean.intensity = mean(intensity)) %>%
  ungroup() %>%
  group_by(patient) %>%
  mutate(intensity.percent = mean.intensity/sum(intensity)*100) %>%
  group_by(patient, cell.type) %>%
  summarise(sum.percent = sum(intensity.percent)) %>%
  pivot_wider(names_from = "patient", values_from = "sum.percent")

UC.df <- UC.df[, c(1,2, 5:12, 3,4)]

write.csv(UC.df, here('outputs', '11_RNAseq-comparison_II', '11_RNAseq_heatmap-UC_percent-table.csv'))
