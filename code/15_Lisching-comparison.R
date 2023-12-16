#####################################################################################################################
##################################### Comparison with Lischnig et al. 2022 ##########################################
#####################################################################################################################

# Libraries
library(here)
library(dplyr)
library(ggplot2)
library(HGNChelper)
library(tidyr)
library(rstatix)
library(ggpubr)

# Directory for outputs
dir.create(here('outputs', '15_Lisching-comparison'))

# Data input: Lischnig et al. (2022)
data.lischnig <- read.csv(here('data', 'Lischnig_2022_suppTableS1.csv'))
## Update the gene names to be comparable with others used in our manuscript
load(here('data', 'genenames_update_20220816.RData')) 
update_data.lischnig <- checkGeneSymbols(data.lischnig$Gene.name, species = "human", map = genenames_newest, unmapped.as.na = FALSE)

## Replace non-approved gene symbols
update_data.lischnig$Suggested.Symbol[update_data.lischnig$Suggested.Symbol == "EPRS1 /// QARS1"] <- "QARS1"
update_data.lischnig$Suggested.Symbol[update_data.lischnig$Suggested.Symbol == "SARS1 /// SARS2"] <- "SARS1"
update_data.lischnig$Suggested.Symbol[update_data.lischnig$Suggested.Symbol == "MPHOSPH6 /// PALS2"] <- "PALS2"
update_data.lischnig$Suggested.Symbol[update_data.lischnig$Suggested.Symbol == "C12orf75 /// FSTL1"] <- "C12orf75"
update_data.lischnig$Suggested.Symbol[update_data.lischnig$Suggested.Symbol == "MASP2 /// BRD8 /// C11orf58 /// KIFAP3"] <- "KIFAP3"

data.lischnig$Suggested.Symbol[duplicated(data.lischnig$Suggested.Symbol)]

## rowID columns
update_data.lischnig$rowID <- rownames(update_data.lischnig)
data.lischnig$rowID <- rownames(data.lischnig)

## Join the dataframes
data.lischnig <- left_join(data.lischnig, update_data.lischnig, by = c("Gene.name" = "x", "rowID" = "rowID"))

# Data input: our data, S+U+B fractions (ie supp table 2)
load(here('outputs', '03_data-processed', '03_data-processed.RData'))

#########################################################################################################
# Small vesicles -> having pvalue < 0.05 & Fold.Change(sEV/lEV) > 1.5
data.lischnig.small <- data.lischnig %>%
  filter(P.value < 0.05 & Fold.change.sEVs.lEVs > 1.5)

# Large vesicles -> having pvalue < 0.05 & Fold.Change(sEV/lEV) < -1.5
data.lischnig.large <- data.lischnig %>%
  filter(P.value < 0.05 & Fold.change.sEVs.lEVs < -1.5)

# Plot the number for large and small vesicles
data.lischnig$vesicle <- ""
data.lischnig$vesicle[data.lischnig$Suggested.Symbol %in% data.lischnig.small$Suggested.Symbol] <- "small"
data.lischnig$vesicle[data.lischnig$Suggested.Symbol %in% data.lischnig.large$Suggested.Symbol] <- "large"

# Map Lischnig data onto ours
data <- data %>%
  left_join(., data.lischnig %>% select(Suggested.Symbol, vesicle), by = c("Suggested.Symbol" = "Suggested.Symbol"))

# Replace NA by not.detected
data$vesicle[is.na(data$vesicle)] <- "not.detected"
data$vesicle[data$vesicle == ""] <- "not.significant"

# Transform to long table  
data.long <- pivot_longer(data, cols = c(U7:B6), names_to = "fraction", values_to = "intensity")
data.long$fraction.ID <- substr(data.long$fraction, 1, 1)

data.long <- data.long %>%
  mutate(fraction.ID = factor(fraction.ID, levels = c("C", "U", "S", "B"))) %>%
  mutate(vesicle  = factor(vesicle, levels = c("not.detected", "not.significant", "small", "large"))) %>%
  mutate(fraction = factor(fraction, levels = c("C1", "C2", "C3", "C4", "C5",
                                                "U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "U9", "U10", "U11",
                                                "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11",
                                                "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11")))
  
# Visualize: all fractions
pdf(here('outputs', '15_Lisching-comparison', '15_barplot-all-fractions.pdf'))
data.long %>%
  filter(intensity != 0) %>%
  group_by(fraction, vesicle, fraction.ID) %>%
  summarise(n_prot = n()) %>%
  ggplot(aes(x = fraction, y = n_prot, fill = fraction.ID))+
    geom_bar(stat = "identity") +
    facet_grid(. ~ vesicle)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Visualize: fractions summed
pdf(here('outputs', '15_Lisching-comparison', '15_barplot-fractions-summarized_a.pdf'))
data.long %>%
  filter(intensity != 0) %>%
  group_by(vesicle, fraction.ID) %>%
  summarise(n_prot = n()) %>%
  ggplot(aes(x = fraction.ID, y = n_prot, fill = fraction.ID))+
  geom_bar(stat = "identity") +
  facet_grid(. ~ vesicle)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(here('outputs', '15_Lisching-comparison', '15_barplot-fractions-summarized_b.pdf'))
data.long %>%
  filter(intensity != 0) %>%
  group_by(vesicle, fraction.ID) %>%
  summarise(n_prot = n()) %>%
  ggplot(aes(x = vesicle, y = n_prot, fill = vesicle))+
  geom_bar(stat = "identity") +
  facet_grid(. ~ fraction.ID)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(here('outputs', '15_Lisching-comparison', '15_barplot-fractions-summarized_c.pdf'))
data.long %>%
  filter(intensity != 0) %>%
  group_by(vesicle, fraction.ID) %>%
  summarise(n_prot = n()) %>%
  ggplot(aes(x = fraction.ID, y = n_prot, fill = vesicle))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(here('outputs', '15_Lisching-comparison', '15_barplot-fractions-percent.pdf'))
data.long %>%
  filter(intensity != 0) %>%
  group_by(fraction.ID) %>% 
  mutate(sum.prot = n()) %>%
  ungroup() %>%
  group_by(fraction.ID, vesicle) %>%
  mutate(sum.vesicle = n()) %>%
  summarise(percent = sum.vesicle/sum.prot*100) %>%
  distinct() %>%
  ggplot(aes(x = fraction.ID, y = percent, fill = vesicle))+
  geom_bar(stat = "identity", position = "stack")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Do statistics on the percent
stat <- data.frame(data.long %>%
  filter(intensity != 0) %>%
  group_by(fraction.ID) %>% 
  mutate(sum.prot = n()) %>%
  ungroup() %>%
  group_by(fraction.ID, vesicle) %>%
  mutate(sum.vesicle = n()) %>%
  summarise(percent = sum.vesicle/sum.prot*100) %>%
  distinct())

res.aov <- anova_test(data = stat, dv = percent, wid = fraction.ID, within = vesicle)
get_anova_table(res.aov)
pwc <- stat %>%
  pairwise_t_test(
    percent ~ vesicle, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

pwc <- pwc %>% add_xy_position(x = "vesicle")

