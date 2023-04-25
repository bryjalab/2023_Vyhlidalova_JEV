################################################################################
###################### COMPARISON WITH IZAR ####################################
################################################################################

# Libraries required:
library(here)
library(dplyr)
library(ComplexHeatmap)

# Create output directory
dir.create(here('outputs', '10_RNAseq-comparison'))

# Load the data
load(here('outputs', '09_controls-addition', '09_patient-specific-data.Rdata'))
RNAseq <- read.csv(here('outputs', '02_data-cleaning', 'updated_RNAseq_S2_table.csv'))
RNAseq$X <- NULL

# Map our data onto Izar's RNAseq data: 8 proteins
RNAseq_8 <- RNAseq[RNAseq$Suggested.Symbol %in% patients.9_ctrl.0$value, ]

RNAseq_8.scaled <- as.matrix(RNAseq_8[, 2:19])
rownames(RNAseq_8.scaled) <- RNAseq_8$Suggested.Symbol
RNAseq_8.scaled <- t(scale(t(RNAseq_8.scaled), scale = TRUE, center = TRUE))

svg(here('outputs', '10_RNAseq-comparison', '10_RNAseq_heatmap-8.svg'))
#pdf(here('outputs', '10_RNAseq-comparison', '10_RNAseq-heatmap-8.pdf'))
Heatmap(RNAseq_8.scaled)
dev.off()

# Map our data onto Izar's RNAseq data: 157 proteins
RNAseq_157 <- RNAseq[RNAseq$Suggested.Symbol %in% patient.specific$Suggested.Symbol, ]

RNAseq_157.scaled <- as.matrix(RNAseq_157[, 2:19])
rownames(RNAseq_157.scaled) <- RNAseq_157$Suggested.Symbol
RNAseq_157.scaled <- t(scale(t(RNAseq_157.scaled), scale = TRUE, center = TRUE))

svg(here('outputs', '10_RNAseq-comparison', '10_RNAseq_heatmap-157.svg'))
#pdf(here('outputs', '10_RNAseq-comparison', '10_RNAseq-heatmap-157.pdf'))
Heatmap(RNAseq_157.scaled, row_names_gp = gpar(fontsize = 2))
dev.off()

