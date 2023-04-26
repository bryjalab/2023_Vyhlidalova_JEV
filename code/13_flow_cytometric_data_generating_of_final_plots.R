
# Libraries required
library(ggplot2)
library(dplyr)
library(tidyr)
library(xlsx)

# Setting of the working directory

setwd("outputs/14_percentages-of-celltypes/")

# Plot with populations of epithelial cells, fibroblasts, macrophages and T cells
## Loading of the data

ep_fib_mac_Tcells <- read.xlsx('ascites_epithelial_cells_fibroblasts_macrophages_Tcells.xlsx', sheetName = 'long format')

## Re-typing and re-ordering of the patient IDs

ep_fib_mac_Tcells$Patient <- as.character(ep_fib_mac_Tcells$Patient)
ep_fib_mac_Tcells$Patient <- factor(ep_fib_mac_Tcells$Patient, 
                                    levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))

## Re-ordering of the celltypes

ep_fib_mac_Tcells$celltype <- factor(ep_fib_mac_Tcells$celltype, 
                                     levels = c('Epithelial', 'Fibroblasts', 'Macrophages', 'T cells'))

## Generating the plot

ggplot(ep_fib_mac_Tcells, 
       aes(x = Patient, y = percentage, fill = celltype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#f2746a','#7bae41','#1ebdbf','#a880b9','#f9a31b','#717171')) +
  ylab('Cell type summary')






# Plot with populations of epithelial cells, fibroblasts and macrophages
## Loading of the data

ep_fib_macs <- read.xlsx('ascites_epithelial_cells_fibroblasts_macrophages.xlsx', sheetName = 'long format')

## Re-typing and re-ordering of the patient IDs

ep_fib_macs$Patient <- as.character(ep_fib_macs$Patient)
ep_fib_macs$Patient <- factor(ep_fib_macs$Patient, 
                              levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))

## Re-ordering of the celltypes

ep_fib_macs$celltype <- factor(ep_fib_macs$celltype, 
                               levels = c('Epithelial cells', 'Fibroblasts', 'Macrophages'))

## Generating the plot

ggplot(ep_fib_macs, 
       aes(x = Patient, y = percentage, fill = celltype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#f2746a','#7bae41','#1ebdbf','#a880b9','#f9a31b','#717171')) +
  ylab('Cell type summary')
