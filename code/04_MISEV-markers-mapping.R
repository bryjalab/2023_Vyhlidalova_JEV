################################################################################
########################### MISEV MARKERS MAPPING ##############################
################################################################################

# Libraries required:
library(here)
library(dplyr)
library(ComplexHeatmap)

# Create a new directory for processed data:
dir.create(here('outputs', '04_MISEV-markers-mapping'))

# Input cleaned data from 03_data-processing.R
load(here("outputs", "03_data-processed", "03_data-processed.RData"))
data.MISEV <- read.csv(here('outputs', '02_data-cleaning', 'updated_MISEV_protein_categories.csv'))

# Map MISEV data onto our protein data
data.merged <- left_join(data.MISEV, data, by = c("Suggested.Symbol" = "Suggested.Symbol"))

data.merged <- data.merged %>%
  select(-c(C1:C5)) %>%
  select(-X, -Approved, -ID, - name)

col_order <- c("category",  "marker", "Suggested.Symbol", "U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11",
               "S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11",
               "B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11")
data.merged <- data.merged[, col_order]

colnames(data.merged)[3] <- "Gene.name.updated"

# Export the table with all MISEV markers (including the ones not present in our data)
write.csv(data.merged, here("outputs", "04_MISEV-markers-mapping", "04_MISEV-markers-mapped_raw-intensities_all-markers.csv"))

# Filter out NA values
data.merged <- na.omit(data.merged)
# Filter out rows with zeros only
data.merged <- data.merged %>%
  rowwise() %>%
  mutate(sums = sum(c_across(U1:B11))) %>%
  filter(sums != 0) %>%
  select(-sums)

# Create a binary matrix
data.merged.binary <- data.merged
data.merged.binary[,4:ncol(data.merged.binary)][data.merged.binary[,4:ncol(data.merged.binary)] >0] <- 1

# Prepare data for heatmap plotting
data.heatmap <- as.matrix(data.merged.binary %>%
                            select(c("U1", "S1", "B1",
                                     "U2", "S2", "B2",
                                     "U3", "S3", "B3",
                                     "U4", "S4", "B4",
                                     "U5", "S5", "B5",
                                     "U6", "S6", "B6",
                                     "U7", "S7", "B7",
                                     "U8", "S8", "B8",
                                     "U9", "S9", "B9",
                                     "U10", "S10", "B10",
                                     "U11", "S11", "B11")))

rownames(data.heatmap) <- data.merged.binary$Gene.name.updated

category <- factor(data.merged.binary$category)
levels(category) <- c("1 - Transmembrane in PM", "2 - Cytosolic",
                      "3 - Non-EV co-isolated", "4 - Transmembrane other",
                      "5 - Secreted")
row_ha <- rowAnnotation(category = category)

col_order <-  c("U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11",
                "S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11",
                "B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11")

data.heatmap <- data.heatmap[, col_order]
col_ha <- HeatmapAnnotation(method = rep(c("UC", "SEC", "Bulk of proteins"), each = 11),
                            col = list(method = c("UC" = "#1ECBE1", "SEC" = "#E11ECB", "Bulk of proteins" = "#000000")))

data.heatmap2 <- as.data.frame(data.heatmap)
data.heatmap2$category <- data.merged.binary$category

# Plot the heatmap
svg(here('outputs', '04_MISEV-markers-mapping', '04_MISEV-markers-heatmap.svg'))
#pdf(here('outputs', '04_MISEV-markers-mapping', '04_MISEV-markers-heatmap.pdf'))
Heatmap(data.heatmap,
        col = c("white", "black"),
        row_names_gp = grid::gpar(fontsize = 4),
        column_names_gp = grid::gpar(fontsize = 9),
        name = "Protein present",
        right_annotation = row_ha,
        top_annotation = col_ha,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        row_split = data.heatmap2$category
)
dev.off()

# Export the table behind the heatmap (binary)
write.csv(data.heatmap, here("outputs", "04_MISEV-markers-mapping", "04_MISEV-markers-mapped_binary-matrix.csv"))

