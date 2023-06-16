
# Libraries required
library(ggplot2)
library(dplyr)
library(tidyr)
library(xlsx)

# Creating the outputs directory
dir.create(here('outputs'))

# Creating the directory for files output
dir.create(here('outputs','13_flow-cytometric-data-output-plots'))

# Setting of the directory for input download
setwd("outputs/12_flow-cytometric-data-percentages-of-celltypes/")

# Plot with populations of epithelial cells, fibroblasts, macrophages and T cells
## Loading of the data

ep_fib_mac_Tcells <- read.xlsx('ascites_main_celltypes_with_Tcells.xlsx', sheetName = 'long format')

## Re-typing and re-ordering of the patient IDs

ep_fib_mac_Tcells$Patient <- as.character(ep_fib_mac_Tcells$Patient)
ep_fib_mac_Tcells$Patient <- factor(ep_fib_mac_Tcells$Patient, 
                                    levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))

## Re-ordering of the celltypes

ep_fib_mac_Tcells$celltype <- factor(ep_fib_mac_Tcells$celltype, 
                                     levels = c('Epithelial', 'Fibroblasts', 'Macrophages', 'T cells'))

## Generating the plot

plot1 <- ggplot(ep_fib_mac_Tcells, 
       aes(x = Patient, y = percentage, fill = celltype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#f2746a','#7bae41','#1ebdbf','#a880b9')) +
  ylab('Cell type summary')

save(plot1, file = here('outputs', '13_flow-cytometric-data-output-plots','percentages_main_celltypes_with_Tcells.svg'))
save(plot1, file = here('outputs', '13_flow-cytometric-data-output-plots','percentages_main_celltypes_with_Tcells.pdf'))




# Plot with populations of epithelial cells, fibroblasts and macrophages
## Loading of the data

ep_fib_macs <- read.xlsx('ascites_main_celltypes.xlsx', sheetName = 'long format')

## Re-typing and re-ordering of the patient IDs

ep_fib_macs$Patient <- as.character(ep_fib_macs$Patient)
ep_fib_macs$Patient <- factor(ep_fib_macs$Patient, 
                              levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))

## Re-ordering of the celltypes

ep_fib_macs$celltype <- factor(ep_fib_macs$celltype, 
                               levels = c('Epithelial cells', 'Fibroblasts', 'Macrophages'))

## Generating the plot

plot2 <- ggplot(ep_fib_macs, 
       aes(x = Patient, y = percentage, fill = celltype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#f2746a','#7bae41','#1ebdbf')) +
  ylab('Cell type summary')


save(plot2, file = here('outputs', '13_flow-cytometric-data-output-plots','percentages_main_celltypes.svg'))
save(plot2, file = here('outputs', '13_flow-cytometric-data-output-plots','percentages_main_celltypes.pdf'))
