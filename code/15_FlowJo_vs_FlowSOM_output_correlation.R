#Libraries required
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rmcorr)
library(ggprism)
library(ggpubr)
library(lme4)
library(lmerTest)
library(effects)

#Creating the outputs directory
dir.create(here('outputs'))

#Creating the directory for files output
dir.create(here('outputs','15_FlowJo_FlowSOM_output_correlation'))

#Loading of the data
data <- read.xlsx(here('data','FlowJo_vs_FlowSOM.xlsx', sheetName = "rmcorr"))

#Ordering of the factor parameters
data$Celltype <- factor(data$Celltype, levels = c("Epithelial cells",
                                                  "Fibroblasts",
                                                  "Macrophages",
                                                  "Monocytes",
                                                  "B cells",
                                                  "CD4 T cells",
                                                  "CD8 T cells",
                                                  "Neutrophils",
                                                  "NK cells"))

data$Patient_ID <- factor(data$Patient_ID, levels = c("1",
                                                      "2",
                                                      "3",
                                                      "4",
                                                      "5",
                                                      "6",
                                                      "7",
                                                      "8",
                                                      "9",
                                                      "10"))


#Repeated measurement correlation
my.rmc.cells <- rmcorr(participant = Celltype, measure1 = FlowJo, measure2 = FlowSOM,
                       dataset = data)


plot <- ggplot(data, aes(x = FlowJo , y = FlowSOM, color = Celltype, group = Celltype)) + 
  geom_line(aes(x = FlowJo, y = my.rmc.cells$model$fitted.values), size = 0.9) + 
  scale_color_manual(values=c("#F2746A", "#2FB34A", "#6F93CC", "#2B3079", "#C694C3", "#9A51A0", "#783294", "#F8A31F", "#717171")) +
  geom_point() + 
  theme_prism() + 
  xlim(0,100) +
  ylim(0, 100) +
  labs(x = "FlowJo (%)", y = "FlowSOM (%)", title = "Repeated measure correlation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 60, y = 14, label = paste("r", "=", round(my.rmc.cells$r, digits = 3), sep = " "), hjust = 0) +
  annotate("text", x = 60, y = 8, label = paste("95% CI:", round(my.rmc.cells$CI[1], digits = 3), "-", round(my.rmc.cells$CI[2], digits = 3)), hjust = 0) +
  annotate("text", x = 60, y = 2, label = paste("p = ", format(my.rmc.cells$p, scientific = TRUE, digits = 3)), hjust = 0) 
save(plot, file = here('outputs','15_FlowJo_FlowSOM_output_correlation','FlowJo_FlowSOM_rmcorr.svg'))
save(plot, file = here('outputs','15_FlowJo_FlowSOM_output_correlation','FlowJo_FlowSOM_rmcorr.pdf'))


#Repeated measurement correlation on log transformed data
data$log10.FlowJo <- log10(data$FlowJo + 0.00001)
data$log10.FlowSOM <- log10(data$FlowSOM + 0.00001)

my.rmc.log <- rmcorr(participant = Celltype , measure1 = log10.FlowJo, measure2 = log10.FlowSOM,
                     dataset = data)

plot_transformed <- ggplot(data, aes(x = log10.FlowJo, y = log10.FlowSOM, color = Celltype , group = Celltype )) + 
  geom_line(aes(x = log10.FlowJo, y = my.rmc.log$model$fitted.values), size = 1) + 
  scale_color_manual(values=c("#F2746A", "#2FB34A", "#6F93CC", "#2B3079", "#C694C3", "#9A51A0", "#783294", "#F8A31F", "#717171")) +
  geom_point() + 
  theme_prism() + 
  xlim(-5, 2) +
  ylim(-5, 2) +
  labs(x = "log10.FlowJo", y = "log10.FlowSOM", title = "RMCorr - transformation log10(x + 0.0001)") +
  annotate("text", x = 0, y = -3, label = paste("r", "=", round(my.rmc.log$r, digits = 3), sep = " "), hjust = 0) +
  annotate("text", x = 0, y = -3.5, label = paste("95% CI:", round(my.rmc.log$CI[1], digits = 3), "-", round(my.rmc.log$CI[2], digits = 3)), hjust = 0) +
  annotate("text", x = 0, y = -4, label = paste("p = ", format(my.rmc.log$p, scientific = TRUE, digits = 3)), hjust = 0) 
save(plot_transformed, file = here('outputs','15_FlowJo_FlowSOM_output_correlation','FlowJo_FlowSOM_rmcorr_transformed.svg'))
save(plot_transformed, file = here('outputs','15_FlowJo_FlowSOM_output_correlation','FlowJo_FlowSOM_rmcorr_transformed.pdf'))

