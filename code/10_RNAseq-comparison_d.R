################################################################################
########################## ENRICHMENT COEFFICIENTS #############################
################################################################################

# Libraries required:
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)
library(rstatix)
library(ggpubr)

# Load the data
primary.cells <- data.frame(patient = c(1:10),
                            epithelial_cancer_cells = c(62.13,	31.30,	90.34,	7.82,	98.85,	86.96,	3.79,	85.85,
                                           43.00,	17.56),
                            fibroblast = c(2.29,	6.15,	1.21,	6.41,	0.88,	11.88,	2.17,	3.88,	1.59,	2.67),
                            macrophages = c(35.58, 62.55, 8.45, 85.77, 0.27,	1.16,	94.03,	10.28,55.40,	79.77
                            ))
primary.cells.long <- pivot_longer(primary.cells, cols = epithelial_cancer_cells:macrophages, names_to = "cell.type", values_to = "percent")
primary.cells.long$patient <- paste("patient_", primary.cells.long$patient, collapse = NULL, sep = "")

load(here('outputs', '11_RNAseq-comparison_II', '11_SEC-UC-long-df.RData'))

### SEC
SEC.merge <- SEC_long %>%
  group_by(patient) %>%
  mutate(intensity_norm = intensity/sum(intensity)*100) %>%
  group_by(patient, cell.type) %>%
  summarise(sum_category = sum(intensity_norm)) %>%
  left_join(., primary.cells.long, by = c("patient" = "patient", "cell.type" = "cell.type")) %>%
  filter(patient != "patient_11") %>%
  mutate(fold_enrichment = sum_category/percent) %>%
  ungroup()
SEC.merge$cell.type <- as.factor(SEC.merge$cell.type)

SEC.merge$fold_enrichment_log2 <- log2(SEC.merge$fold_enrichment)

### UC
UC.merge <- UC_long %>%
  group_by(patient) %>%
  mutate(intensity_norm = intensity/sum(intensity)*100) %>%
  group_by(patient, cell.type) %>%
  summarise(sum_category = sum(intensity_norm)) %>%
  left_join(., primary.cells.long, by = c("patient" = "patient", "cell.type" = "cell.type")) %>%
  filter(patient != "patient_11") %>%
  mutate(fold_enrichment = sum_category/percent) %>%
  ungroup()
UC.merge$cell.type <- as.factor(UC.merge$cell.type)

UC.merge$fold_enrichment_log2 <- log2(UC.merge$fold_enrichment)

# Compute medians
SEC.merge %>%
  group_by(cell.type) %>%
  summarise(median.FE = median(fold_enrichment),
            median.FE.log = median(fold_enrichment_log2))

UC.merge %>%
  group_by(cell.type) %>%
  summarise(median.FE = median(fold_enrichment),
            median.FE.log = median(fold_enrichment_log2))

# Compute means
SEC.merge %>%
  group_by(cell.type) %>%
  summarise(mean.FE = mean(fold_enrichment),
            mean.FE.log = mean(fold_enrichment_log2))

UC.merge %>%
  group_by(cell.type) %>%
  summarise(mean.FE = mean(fold_enrichment),
            mean.FE.log = mean(fold_enrichment_log2))


# Repeated measures ANOVA

# SEC
# Compute repeated measures anova
res.aov.SEC <- anova_test(data = SEC.merge, dv = fold_enrichment_log2 , wid = patient, within = cell.type)
get_anova_table(res.aov.SEC)
pwc.SEC <- SEC.merge %>%
  pairwise_t_test(
    fold_enrichment_log2 ~ cell.type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc.SEC

pwc.SEC <- pwc.SEC %>% add_xy_position(x = "cell.type")

svg(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_SEC_anova.svg"))
pdf(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_SEC_anova.pdf"))
q <- SEC.merge %>% 
  ggplot(aes(x = cell.type, y = fold_enrichment_log2)) +
  geom_point(size = 3, aes(col = cell.type, shape = cell.type)) +
  theme_prism(axis_text_angle = 45, base_size = 12)+
  stat_summary(fun="mean", geom="point", color="black", size=2,
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) }) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.4, yend=..y..),  color="black")+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.4, yend=..y..),  color="black")+
  scale_shape_prism()

q + add_pvalue(pwc.SEC, remove.bracket = TRUE, y.position = c(10, 11, 10))
dev.off()

# UC
res.aov.UC <- anova_test(data = UC.merge, dv = fold_enrichment_log2 , wid = patient, within = cell.type)
get_anova_table(res.aov.UC)
pwc.UC <- UC.merge %>%
  pairwise_t_test(
    fold_enrichment_log2 ~ cell.type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc.UC

pwc.UC <- pwc.UC %>% add_xy_position(x = "cell.type")

svg(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_UC_anova.svg"))
pdf(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_UC_anova.pdf"))
q <- UC.merge %>% 
  ggplot(aes(x = cell.type, y = fold_enrichment_log2)) +
  geom_point(size = 3, aes(col = cell.type, shape = cell.type)) +
  theme_prism(axis_text_angle = 45, base_size = 12)+
  stat_summary(fun="mean", geom="point", color="black", size=2,
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) }) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.4, yend=..y..),  color="black")+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.4, yend=..y..),  color="black")+
  scale_shape_prism()

q + add_pvalue(pwc.UC, remove.bracket = TRUE, y.position = c(10, 11, 10))
dev.off()

####################################################################################################
################ PRIMARY CELLS CONTAINING ALSO PERIPHERAL BLOOD CELLS ##############################
####################################################################################################

# Load the data
primary.cells.all <- read.csv2(here('data', 'flow_data_percentages_all_populations_final.csv'))
colnames(primary.cells.all)[2] <- "epithelial_cancer_cells"
colnames(primary.cells.all)[3] <- "fibroblast"
colnames(primary.cells.all)[4] <- "macrophages"

primary.cells.all.long <- pivot_longer(primary.cells.all, cols = epithelial_cancer_cells:NK_cells, names_to = "cell.type", values_to = "percent")
primary.cells.all.long$patient <- paste("patient_", primary.cells.all.long$Patient_ID, collapse = NULL, sep = "")

primary.cells.all.long$percent <- as.numeric(primary.cells.all.long$percent)

################################# Merge with SEC and UC by MS ################################

### SEC

SEC.merge.all <- SEC_long %>%
  group_by(patient) %>%
  mutate(intensity_norm = intensity/sum(intensity)*100) %>%
  group_by(patient, cell.type) %>%
  summarise(sum_category = sum(intensity_norm)) %>%
  left_join(., primary.cells.all.long, by = c("patient" = "patient", "cell.type" = "cell.type")) %>%
  filter(patient != "patient_11") %>%
  mutate(fold_enrichment = sum_category/percent) %>%
  ungroup()
SEC.merge.all$cell.type <- as.factor(SEC.merge$cell.type)

SEC.merge.all$fold_enrichment_log2 <- log2(SEC.merge.all$fold_enrichment)

res.aov.SEC.all <- anova_test(data = SEC.merge.all, dv = fold_enrichment_log2 , wid = patient, within = cell.type)
get_anova_table(res.aov.SEC.all)
pwc.SEC.all <- SEC.merge.all %>%
  pairwise_t_test(
    fold_enrichment_log2 ~ cell.type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc.SEC.all

pwc.SEC.all <- pwc.SEC.all %>% add_xy_position(x = "cell.type")

svg(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_SEC_anova_primaryCells.svg"))
pdf(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_SEC_anova_primaryCells.pdf"))
q <- SEC.merge.all %>% 
  ggplot(aes(x = cell.type, y = fold_enrichment_log2)) +
  geom_point(size = 3, aes(col = cell.type, shape = cell.type)) +
  theme_prism(axis_text_angle = 45, base_size = 12)+
  stat_summary(fun="mean", geom="point", color="black", size=2,
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) }) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.4, yend=..y..),  color="black")+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.4, yend=..y..),  color="black")+
  scale_shape_prism()

q + add_pvalue(pwc.SEC.all, remove.bracket = TRUE, y.position = c(10, 11, 10))
dev.off()

### UC

UC.merge.all <- UC_long %>%
  group_by(patient) %>%
  mutate(intensity_norm = intensity/sum(intensity)*100) %>%
  group_by(patient, cell.type) %>%
  summarise(sum_category = sum(intensity_norm)) %>%
  left_join(., primary.cells.all.long, by = c("patient" = "patient", "cell.type" = "cell.type")) %>%
  filter(patient != "patient_11") %>%
  mutate(fold_enrichment = sum_category/percent) %>%
  ungroup()
UC.merge.all$cell.type <- as.factor(UC.merge$cell.type)

UC.merge.all$fold_enrichment_log2 <- log2(UC.merge.all$fold_enrichment)

res.aov.UC.all <- anova_test(data = UC.merge.all, dv = fold_enrichment_log2 , wid = patient, within = cell.type)
get_anova_table(res.aov.UC.all)
pwc.UC.all <- UC.merge.all %>%
  pairwise_t_test(
    fold_enrichment_log2 ~ cell.type, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc.UC.all

pwc.UC.all <- pwc.UC.all %>% add_xy_position(x = "cell.type")

svg(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_UC_anova_primaryCells.svg"))
pdf(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_UC_anova_primaryCells.pdf"))
q <- UC.merge.all %>% 
  ggplot(aes(x = cell.type, y = fold_enrichment_log2)) +
  geom_point(size = 3, aes(col = cell.type, shape = cell.type)) +
  theme_prism(axis_text_angle = 45, base_size = 12)+
  stat_summary(fun="mean", geom="point", color="black", size=2,
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) }) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.4, yend=..y..),  color="black")+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.4, yend=..y..),  color="black")+
  scale_shape_prism()

q + add_pvalue(pwc.UC.all, remove.bracket = TRUE, y.position = c(10, 11, 10))
dev.off()
