
# Libraries required
library(flowCore)

# Loading of preprocessed samples
P1 <- flowCore::read.FCS("P1_viable_cells.fcs", truncate_max_range = FALSE)
P2 <- flowCore::read.FCS("P2_viable_cells.fcs", truncate_max_range = FALSE)
P3 <- flowCore::read.FCS("P3_viable_cells.fcs", truncate_max_range = FALSE)
P4 <- flowCore::read.FCS("P4_viable_cells.fcs", truncate_max_range = FALSE)
P5 <- flowCore::read.FCS("P5_viable_cells.fcs", truncate_max_range = FALSE)
P6 <- flowCore::read.FCS("P6_viable_cells.fcs", truncate_max_range = FALSE)
P7 <- flowCore::read.FCS("P7_viable_cells.fcs", truncate_max_range = FALSE)
P8 <- flowCore::read.FCS("P8_viable_cells.fcs", truncate_max_range = FALSE)
P9 <- flowCore::read.FCS("P9_viable_cells.fcs", truncate_max_range = FALSE)
P10 <- flowCore::read.FCS("P10_viable_cells.fcs", truncate_max_range = FALSE)
P12 <- flowCore::read.FCS("P12_viable_cells.fcs", truncate_max_range = FALSE)
P13 <- flowCore::read.FCS("P13_viable_cells.fcs", truncate_max_range = FALSE)
P14 <- flowCore::read.FCS("P14_viable_cells.fcs", truncate_max_range = FALSE)

# Specification of the channels selected for logicle transformation
selected_channels <- c("FJComp-APC-A",
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


# File transformation via estimateLogicle

P1_transformList <- estimateLogicle(P1, channels = selected_channels)
P1_transformed <- flowCore::transform(P1, P1_transformList)

P2_transformList <- estimateLogicle(P2, channels = selected_channels)
P2_transformed <- flowCore::transform(P2, P2_transformList)

P3_transformList <- estimateLogicle(P3, channels = selected_channels)
P3_transformed <- flowCore::transform(P3, P3_transformList)

P4_transformList <- estimateLogicle(P4, channels = selected_channels)
P4_transformed <- flowCore::transform(P4, P4_transformList)

P5_transformList <- estimateLogicle(P5, channels = selected_channels)
P5_transformed <- flowCore::transform(P5, P5_transformList)

P6_transformList <- estimateLogicle(P6, channels = selected_channels)
P6_transformed <- flowCore::transform(P6, P4_transformList)

P7_transformList <- estimateLogicle(P7, channels = selected_channels)
P7_transformed <- flowCore::transform(P7, P7_transformList)

P8_transformList <- estimateLogicle(P8, channels = selected_channels)
P8_transformed <- flowCore::transform(P8, P8_transformList)

P9_transformList <- estimateLogicle(P9, channels = selected_channels)
P9_transformed <- flowCore::transform(P9, P9_transformList)

P10_transformList <- estimateLogicle(P10, channels = selected_channels)
P10_transformed <- flowCore::transform(P10, P10_transformList)

P12_transformList <- estimateLogicle(P12, channels = selected_channels)
P12_transformed <- flowCore::transform(P12, P12_transformList)

P13_transformList <- estimateLogicle(P13, channels = selected_channels)
P13_transformed <- flowCore::transform(P13, P13_transformList)

P14_transformList <- estimateLogicle(P14, channels = selected_channels)
P14_transformed <- flowCore::transform(P14, P14_transformList)


# Export of the transformed fcs files 
write.FCS(P1_transformed, filename = "P1_transformed.fcs")
write.FCS(P2_transformed, filename = "P2_transformed.fcs")
write.FCS(P3_transformed, filename = "P3_transformed.fcs")
write.FCS(P4_transformed, filename = "P4_transformed.fcs")
write.FCS(P5_transformed, filename = "P5_transformed.fcs")
write.FCS(P6_transformed, filename = "P6_transformed.fcs")
write.FCS(P7_transformed, filename = "P7_transformed.fcs")
write.FCS(P8_transformed, filename = "P8_transformed.fcs")
write.FCS(P9_transformed, filename = "P9_transformed.fcs")
write.FCS(P10_transformed, filename = "P10_transformed.fcs")
write.FCS(P12_transformed, filename = "P12_transformed.fcs")
write.FCS(P13_transformed, filename = "P13_transformed.fcs")
write.FCS(P14_transformed, filename = "P14_transformed.fcs")
