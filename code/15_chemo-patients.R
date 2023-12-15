#####################################################################################################################
##################################### Patients after chemotherapy addition ##########################################
#####################################################################################################################

# Libraries
library(here)
library(dplyr)
library(ComplexHeatmap)
library(UpSetR)

# Create a new directory for processed data:
dir.create(here('outputs', '15_chemo-patients'))

# Data input: chemo samples
data <- read.csv(here('outputs', '02_data-cleaning', 'updated_proteinGroups.csv'), stringsAsFactors = FALSE)

## Select columns containing intensities
data <- data %>%
  select(Suggested.Symbol, name, ID, starts_with("Intensity.")) 

## Rename the columns
colnames(data) <- gsub("Intensity.", "", colnames(data))

## Select the chemo samples
data.chemo <- data %>%
  select(Suggested.Symbol, ID, a71, a54, a55, a65)

data.chemo <- data.chemo %>%
  rename(CH1 = a54,
         CH2 = a55,
         CH3 = a65,
         CH4 = a71)

data.chemo <- data.chemo[, c(1, 2, 4:6, 3)]

rm(data)

#################################### MISEV mapping #######################################################

data.MISEV <- read.csv(here('outputs', '02_data-cleaning', 'updated_MISEV_protein_categories.csv'))

# Map MISEV data onto our protein data
data.merged <- left_join(data.MISEV, data.chemo, by = c("Suggested.Symbol" = "Suggested.Symbol"))

# Remove unnecessary columns
data.merged <- data.merged %>%
  select(-X, -Approved, -ID)

# Define columns order
col_order <- c("category",  "marker", "Suggested.Symbol", "CH1", "CH2", "CH3", "CH4")
data.merged <- data.merged[, col_order]

colnames(data.merged)[3] <- "Gene.name.updated"

# Export the table with all MISEV markers (including the ones not present in our data)
#write.csv(data.merged, here("outputs", "04_MISEV-markers-mapping", "04_MISEV-markers-mapped_raw-intensities_all-markers.csv"))

# Filter out NA values
data.merged <- na.omit(data.merged)

# Filter out rows with zeros only
data.merged <- data.merged %>%
  rowwise() %>%
  mutate(sums = sum(c_across(CH1:CH4))) %>%
  filter(sums != 0) %>%
  select(-sums)

# Create a binary matrix
data.merged.binary <- data.merged
data.merged.binary[,4:ncol(data.merged.binary)][data.merged.binary[,4:ncol(data.merged.binary)] >0] <- 1

# Prepare data for heatmap plotting
data.heatmap <- as.matrix(data.merged.binary %>%
                            select(c("CH1", "CH2", "CH3", "CH4")))

rownames(data.heatmap) <- data.merged.binary$Gene.name.updated

category <- factor(data.merged.binary$category)
levels(category) <- c("1 - Transmembrane in PM", "2 - Cytosolic",
                      "3 - Non-EV co-isolated", "4 - Transmembrane other",
                      "5 - Secreted")
row_ha <- rowAnnotation(category = category)

col_order <-  c("CH1", "CH2", "CH3", "CH4")

data.heatmap <- data.heatmap[, col_order]
data.heatmap2 <- as.data.frame(data.heatmap)
data.heatmap2$category <- data.merged.binary$category

# Plot the heatmap
svg(here('outputs', '15_chemo-patients', '15_MISEV-markers-heatmap.svg'))
#pdf(here('outputs', '15_chemo-patients', '15_MISEV-markers-heatmap.pdf'))
Heatmap(data.heatmap,
        col = c("white", "black"),
        row_names_gp = grid::gpar(fontsize = 4),
        column_names_gp = grid::gpar(fontsize = 9),
        name = "Protein present",
        right_annotation = row_ha,
#        top_annotation = col_ha,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        row_split = data.heatmap2$category
)
dev.off()

rm(data.heatmap, data.heatmap2, data.merged, data.merged.binary, data.MISEV, row_ha, category, col_order)

#################################### Addition to UC+ctrls #######################################################

# Load table (from 09_controls-addition) containing controls mapped on patients
patients <- read.csv(here('outputs', '09_controls-addition', '09_controls-mapped-on-patients.csv'))
patients$X <- NULL

# Map chemo patients onto UC-measured patients+controls
patients <- left_join(patients, data.chemo, by = c("value" = "Suggested.Symbol", "ID" = "ID"))

# Create binary matrix from chemo patients
patients[,20:23][patients[, 20:23] >0] <- 1

# Firstly, find proteins which are in all 11 patients, but no controls, or at least 9/11 patients, respectively
patients.11_ctrl.0 <- patients %>%
  rowwise() %>%
  mutate(sum_patients = sum(c_across(patient_1: patient_11))) %>%
  mutate(sum_ctrl = sum(c_across(C1:C5)>0)) %>%
  filter(sum_ctrl == 0) %>%
  filter(sum_patients == 11)

patients.11_ctrl.0.heatmap <- as.matrix(patients.11_ctrl.0[,4:23])
rownames(patients.11_ctrl.0.heatmap) <- patients.11_ctrl.0$value
pdf(here('outputs', '15_chemo-patients', '15_UC-ctrl_patients11.pdf'))
Heatmap(patients.11_ctrl.0.heatmap,
        cluster_columns = FALSE,
        col = c("white", "black"))
dev.off()

patients.9_ctrl.0 <- patients %>%
  rowwise() %>%
  mutate(sum_patients = sum(c_across(patient_1: patient_11))) %>%
  mutate(sum_ctrl = sum(c_across(C1:C5)>0)) %>%
  filter(sum_ctrl == 0) %>%
  filter(sum_patients > 8)

patients.9_ctrl.0.heatmap <- as.matrix(patients.9_ctrl.0[,4:23])
rownames(patients.9_ctrl.0.heatmap) <- patients.9_ctrl.0$value
pdf(here('outputs', '15_chemo-patients', '15_UC-ctrl_patients9.pdf'))
Heatmap(patients.9_ctrl.0.heatmap,
        cluster_columns = FALSE,
        col = c("white", "black"))
dev.off()

patients.157 <- read.csv(here('outputs', '09_controls-addition', '09_157-proteins_patient-specific.csv'))
patients.157$X <- NULL
patients.157 <- left_join(patients.157, data.chemo, by = c("Suggested.Symbol" = "Suggested.Symbol"))
# Create binary matrix from chemo patients
patients.157[,21:24][patients.157[, 21:24] >0] <- 1

patients.157.heatmap <- patients.157 %>% select(-sum_patients, -sum_ctrl, -ID)
patients.157.heatmap <- as.matrix(patients.157.heatmap[,2:21])
rownames(patients.157.heatmap) <- patients.157$Suggested.Symbol 

pdf(here('outputs', '15_chemo-patients', '15_UC-ctrl_patients157.pdf'))
Heatmap(patients.157.heatmap,
        cluster_columns = FALSE,
        row_names_gp = grid::gpar(fontsize = 4),
        column_names_gp = grid::gpar(fontsize = 9),
        col = c("white", "black"))
dev.off()
