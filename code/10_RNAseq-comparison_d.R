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
                            epithelial_cancer_cells = c(49.57,	16.69,	54.19,	2.17,	92.07,	46.04,	0.82,	46.80,	2.36,	2.18),
                            fibroblast = c(1.83,	3.28,	0.72,	1.78,	0.82,	6.29,	0.47,2.11,	0.09,	0.33),
                            macrophages = c(28.39,	33.36,	5.07,	23.85,	0.25,	0.61,	20.29	,5.60	,3.05	,9.88))

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

# SEC.pval <- SEC.merge %>%
#     emmeans_test(fold_enrichment_log2 ~ cell.type, p.adjust.method = "bonferroni") %>%
#     rstatix::add_xy_position()
# 
# # Plot with legend
# svg(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_SEC.svg"))
# pdf(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_SEC.pdf"))
# q <- SEC.merge %>% 
#   ggplot(aes(x = cell.type, y = fold_enrichment_log2)) +
#   geom_point(size = 3, aes(col = cell.type, shape = cell.type)) +
#   theme_prism(axis_text_angle = 45, base_size = 12)+
#   stat_summary(fun="mean", geom="point", color="black", size=2,
#                fun.min = function(z) { quantile(z,0.25) },
#                fun.max = function(z) { quantile(z,0.75) }) +
#   stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.4, yend=..y..),  color="black")+
#   stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.4, yend=..y..),  color="black")+
#   scale_shape_prism()
# 
# q + add_pvalue(SEC.pval, remove.bracket = TRUE, y.position = c(10, 11, 10))
# dev.off()

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

# UC.pval <- UC.merge %>%
#   emmeans_test(fold_enrichment_log2 ~ cell.type, p.adjust.method = "bonferroni") %>%
#   rstatix::add_xy_position()
# 
# # Plot with legend
# svg(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_UC.svg"))
# pdf(here("outputs", "11_RNASeq-comparison_II", "11_fold-enrichment_UC.pdf"))
# q <- UC.merge %>% 
#   ggplot(aes(x = cell.type, y = fold_enrichment_log2)) +
#   geom_point(size = 3, aes(col = cell.type, shape = cell.type)) +
#   theme_prism(axis_text_angle = 45, base_size = 12)+
#   stat_summary(fun="mean", geom="point", color="black", size=2,
#                fun.min = function(z) { quantile(z,0.25) },
#                fun.max = function(z) { quantile(z,0.75) }) +
#   stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.4, yend=..y..),  color="black")+
#   stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.4, yend=..y..),  color="black")+
#   scale_shape_prism()
# 
# q + add_pvalue(UC.pval, remove.bracket = TRUE, y.position = c(10, 11, 10))
# dev.off()

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
