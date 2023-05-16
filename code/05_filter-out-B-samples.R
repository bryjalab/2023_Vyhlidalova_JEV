################################################################################
########################     Filter-out B samples     ##########################
################################################################################

# Libraries required:
library(here)
library(dplyr)
library(ggplot2)
library(tidyr)

source(here('code', "functions.R"))

# Create a new directory for processed data:
dir.create(here('outputs', '05_filter-out-B-samples'))

# Load the input data
load(here("outputs", "03_data-processed", "03_data-processed.RData"))

# Plot the number of proteins per each method
svg(here("outputs", "05_filter-out-B-samples", "05_barplot-protein-numbers-methods.svg"))
pdf(here("outputs", "05_filter-out-B-samples", "05_barplot-protein-numbers-methods.pdf"))
data %>%
  pivot_longer(cols = c(4:42), names_to = "Condition", values_to = "Intensity") %>%
  mutate(Method = case_when(grepl("S", Condition) ~ "S",
                            grepl("B", Condition) ~ "B",
                            TRUE ~ "U") ) %>%
  group_by(Condition, Method) %>%
  summarize(n = sum(Intensity > 0)) %>%
  filter(Condition %in% c("U1", "S1", "B1",
                          "U2", "S2", "B2",
                          "U3", "S3", "B3",
                          "U4", "S4", "B4",
                          "U5", "S5", "B5",
                          "U6", "S6", "B6",
                          "U7", "S7", "B7",
                          "U8", "S8", "B8",
                          "U9", "S9", "B9",
                          "U10", "S10", "B10",
                          "U11", "S11", "B11")) %>%
  mutate(Condition = factor(Condition, levels = c("U1", "S1", "B1",
                                                  "U2", "S2", "B2",
                                                  "U3", "S3", "B3",
                                                  "U4", "S4", "B4",
                                                  "U5", "S5", "B5",
                                                  "U6", "S6", "B6",
                                                  "U7", "S7", "B7",
                                                  "U8", "S8", "B8",
                                                  "U9", "S9", "B9",
                                                  "U10", "S10", "B10",
                                                  "U11", "S11", "B11"))) %>%
  ggplot(., aes(x = Condition, y = n, fill = Method))+
  geom_bar(stat= "identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("#848382", "#E11ECB", "#1ECBE1"))+
  labs(x = "Patient",
       y = "Number of proteins")
dev.off()

# Filter-out B samples
tmp <- data.frame(Suggested.Symbol = data$Suggested.Symbol)
data.filtered <- data.frame(Suggested.Symbol = data$Suggested.Symbol)

data.filtered <- filter_B(data, B1, U1, S1)
data.filtered <- filter_B(data, B2, U2, S2)
data.filtered <- filter_B(data, B3, U3, S3)
data.filtered <- filter_B(data, B4, U4, S4)
data.filtered <- filter_B(data, B5, U5, S5)
data.filtered <- filter_B(data, B6, U6, S6)
data.filtered <- filter_B(data, B7, U7, S7)
data.filtered <- filter_B(data, B8, U8, S8)
data.filtered <- filter_B(data, B9, U9, S9)
data.filtered <- filter_B(data, B10, U10, S10)
data.filtered <- filter_B(data, B11, U11, S11)

colnames(data.filtered)[1] <- "Suggested.Symbol"

data.filtered <- left_join(data.filtered, data %>% select(Suggested.Symbol, name, ID), by = c("Suggested.Symbol" = "Suggested.Symbol"))
data.filtered$Suggested.Symbol[duplicated(data.filtered$Suggested.Symbol)]

# Save filtered data
write.csv(data.filtered, here("outputs", "05_filter-out-B-samples", "05_data_filtered-out-B-samples.csv"))
save(data.filtered, file = here("outputs", "05_filter-out-B-samples", "05_data_filtered-out-B-samples.Rdata"))
save(data.filtered, file = here("outputs", "05_filter-out-B-samples", "05_data.filtered.Rdata"))
