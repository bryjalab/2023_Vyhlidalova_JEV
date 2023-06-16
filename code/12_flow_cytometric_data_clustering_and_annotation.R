
# Libraries required
library(CATALYST)
library(flowCore)
library(SummarizedExperiment)
library(pheatmap)
library(ggplot2)
library(FlowSOM)
library(Seurat)
library(xlsx)

# Creating the outputs directory
dir.create(here('outputs'))

# Creating the directory for files output
dir.create(here('outputs','12_flow-cytometric-data-percentages-of-celltypes'))
dir.create(here('outputs','13_flow-cytometric-data-output-plots'))

# Loading of transformed files
fcs_files <- list.files(pattern = ".fcs$")

fs <- read.flowSet(fcs_files, 
                   transformation = FALSE, 
                   truncate_max_range = FALSE)

# Panel annotation

channel_name <- c("FJComp-APC-A",
                  "FJComp-APC-Fire 750-A",
                  "FJComp-APC-Fire 810-A",
                  "FJComp-Alexa Fluor 647-A",
                  "FJComp-Alexa Fluor 700-A",
  
                  "FJComp-BB515-A",
                  "FJComp-BB700-A",
                  "FJComp-BUV395-A",
                  "FJComp-BUV496-A",
                  "FJComp-BUV563-A",
  
                  "FJComp-BUV661-A",
                  "FJComp-BUV737-A",
                  "FJComp-BUV805-A",
                  "FJComp-BV421-A",
                  "FJComp-BV510-A",
  
                  "FJComp-BV570-A",
                  "FJComp-BV650-A",
                  "FJComp-BV711-A",
                  "FJComp-PE-A",
                  "FJComp-PE-CF594-A",
  
                  "FJComp-PE-Fire 810-A",
                  "FJComp-PerCP-A",
                  "FJComp-Super Bright 600-A")

marker <- c("CD204",
            "CD88",
            "CD8",
            "Lyve1",
            "EpCAM",
  
            "CD25",
            "CD95",
            "CD15",
            "CD9",
            "CD90",
  
            "CD169",
            "CD11b",
            "CD16",
            "CD14",
            "CD56",
  
            "CD3",
            "CD19",
            "CD206",
            "FOLR2",
            "CD163",
  
            "HLA-DR",
            "CD45",
            "VSIG4")

marker_class <- c("type",
                  "type",
                  "type",
                  "type",
                  "type",
                  
                  "type",
                  "type",
                  "type",
                  "type",
                  "type",
                  
                  "type",
                  "type",
                  "type",
                  "type",
                  "type",
                  
                  "type",
                  "type",
                  "type",
                  "type",
                  "type",
                  
                  "type",
                  "type",
                  "type")

panel <- data.frame(channel_name,marker,marker_class)


# Experimental design annotation

file_name <- c("P1_transformed.fcs",
               "P2_transformed.fcs",
               "P3_transformed.fcs",
               "P4_transformed.fcs",
               "P5_transformed.fcs",
               "P6_transformed.fcs",
               "P7_transformed.fcs",
               "P8_transformed.fcs",
               "P9_transformed.fcs",
               "P10_transformed.fcs",
               "P12_transformed.fcs",
               "P13_transformed.fcs",
               "P14_transformed.fcs")

sample_id <- c("01",
               "02",
               "03",
               "04",
               "05",
               "06",
               "07",
               "08",
               "09",
               "10",
               "12",
               "13",
               "14")

exp_info <- data.frame(file_name, sample_id)


# Dataset preparation:

ascites <- prepData(fs, 
                    panel, 
                    exp_info,
                    panel_cols = list(channel = "channel_name", antigen = "marker"),
                    md_cols = list(file = "file_name", id = "sample_id"),
                    FACS = TRUE)

# Dimensional reduction (tSNE)

ascites <- runDR(ascites, 
                 dr = "TSNE", 
                 cells = NULL, 
                 features = NULL,
                 assay = "counts",
                 seed = 1,
                 verbose = TRUE)

# FlowSOM clustering

ascites <- cluster(ascites, 
                   features = c("CD3",
                                "CD19",
                                "CD11b",
                                "CD15",
                                "CD16",
                                "CD56",
                                "CD14",
                                "CD163",
                                "CD88",
                                "EpCAM",
                                "FOLR2",
                                "CD8",
                                "CD90",
                                "CD45"), 
                   xdim = 10, 
                   ydim = 10, 
                   maxK = 35, 
                   verbose = TRUE, 
                   seed = 1)

plotDR(ascites, dr = "TSNE", color_by = "meta25")

# Filtering out of unknown populations

ascites_filtered <- filterSCE(ascites, 
                              k = 'meta25',
                              cluster_id %in% c(1,2,3,4,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25))

# Dimensional reduction no. 2

ascites_filtered <- runDR(ascites_filtered, 
                          dr = "TSNE", 
                          cells = NULL, 
                          features = c("CD3",
                                       "CD8",
                                       "CD19",
                                       "CD45",
                                       "CD11b",
                                       "CD56",
                                       "CD90",
                                       "CD16",
                                       "CD15",
                                       "EpCAM",
                                       "CD14",
                                       "HLA-DR",
                                       "FOLR2",
                                       "CD206"),
                           assay = "counts",
                           seed = 1,
                           verbose = TRUE)

# Re-clustering of the dataset

ascites_filtered <- cluster(ascites_filtered, 
                            features = c("CD3",
                                         "CD8",
                                         "CD19",
                                         "CD45",
                                         "CD11b",
                                         "CD56",
                                         "CD90",
                                         "CD16",
                                         "CD15",
                                         "EpCAM",
                                         "CD14",
                                         "HLA-DR",
                                         "FOLR2",
                                         "CD206"), 
                            xdim = 10, 
                            ydim = 10, 
                            maxK = 35, 
                            verbose = TRUE, 
                            seed = 1)

# Creating annotation table

old_cluster <- c("1",
                  "2",
                  "3",
                  "4",
                  "5",
                  "6",
                  "7",
                  "8",
                  "9",
                  "10",
                  
                  "11",
                  "12",
                  "13",
                  "14",
                  "15",
                  "16",
                  "17",
                  "18",
                  "19",
                  "20",
                  
                  "21",
                  "22",
                  "23",
                  "24",
                  "25",
                  "26",
                  "27",
                  "28",
                  "29",
                  "30")

new_cluster <- c("Epithelial cells",
                  "Epithelial cells",
                  "CD4 T cells",
                  "NK cells",
                  "Fibroblasts",
                  "NK cells",
                  "CD4 T cells",
                  "Epithelial cells",
                  "Epithelial cells",
                  "Fibroblasts",
                  
                  "Epithelial cells",
                  "Epithelial cells",
                  "Neutrophils",
                  "Fibroblasts",
                  "CD8 T cells",
                  "B cells",
                  "Neutrophils",
                  "Neutrophils",
                  "Epithelial cells",
                  "Macrophages",
                  
                  "Macrophages",
                  "Macrophages",
                  "Macrophages",
                  "Monocytes",
                  "Macrophages",
                  "Macrophages",
                  "Macrophages",
                  "Macrophages",
                  "Macrophages",
                  "CD8 T cells")

annotation_table <- data.frame(old_cluster, new_cluster)

# Dataset annotation

ascites_annotated <- mergeClusters(ascites_filtered, 
                                   k = "meta30",
                                   id = "Celltypes",
                                   table = annotation_table)

ascites_annotated <- filterSCE(ascites_annotated,
                               sample_id %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10"))

# Plot of t-SNE with annotated populations

plot1 <- plotDR(ascites_annotated, 
       dr = "TSNE", 
       color_by = "Celltypes", 
       k_pal = c('#f2746a','#7bae41','#039be5','#1A237E','#ce93d8','#ab47bc','#7b1fa2','#f9a31b','#717171'))

save(plot1, file = here('outputs','13_flow-cytometric-data-output-plots','ascites_all_populations_tSNE.svg'))
save(plot1, file = here('outputs','13_flow-cytometric-data-output-plots','ascites_all_populations_tSNE.pdf'))


# Plot with abundances of annotated populations

plot2 <- plotAbundances(ascites_annotated, 
               k = "Celltypes", 
               by = "sample_id", 
               group_by = NULL, 
               col_clust = FALSE, 
               k_pal = c('#f2746a','#7bae41','#039be5','#1A237E','#ce93d8','#ab47bc','#7b1fa2','#f9a31b','#717171'))

save(plot2, file = here('outputs','13_flow-cytometric-data-output-plots','percentages_all_populations.svg'))
save(plot2, file = here('outputs','13_flow-cytometric-data-output-plots','percentages_all_populations.pdf'))


# Extraction of counts of distinct celltypes across patients

ascites_annotated$Celltypes <- cluster_ids(ascites_annotated, "Celltypes")

ascites_seurat <- as.Seurat(ascites_annotated,
                            counts = "counts",
                            data = NULL)

annotated_populations <- table(ascites_seurat@meta.data$sample_id, 
                               group.by = ascites_seurat@meta.data$Celltypes)

# Export of the celltype counts across patients
write.csv(annotated_populations, file = here('outputs', '12_flow-cytometric-data-percentages-of-celltypes', 'ascites_annotated_populations_counts.csv'))

