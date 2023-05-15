################################################################################
############################     DLS Analysis     ##############################
################################################################################

# Libraries required:
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)

# Create outputs directory
dir.create(here('outputs'))

# Create directory for files output
dir.create(here('outputs', '01_DLS-data'))

# Data input: DLS analysis output, in long format
data <- read.delim(here('data', 'DSV_analysis_input.txt'))

# Create unique IDs for each sample
id <- paste(1:150, "id", sep = "_")
data$id <- rep(id, nrow(data)/150)

# Remove patients 3 and 4
data <- data %>%
  filter(! Patient %in% c("U3", "U4", "S3", "S4", "B3", "B4"))

# Summarize (mean) Number of particles and Size of the particles
data.summarized <- data %>%
  group_by(Patient, id) %>%
  summarise(mean_number = mean(Number),
            mean_size = mean(Size))

# Plot the graphs
for (i in levels(factor(data.summarized$Patient))) {
  file_name <- paste(i, ".svg", sep = "")
  gg <- data.summarized %>%
    filter(Patient == i) %>%
    ggplot(., aes(x = mean_size, y = mean_number, group = 1))+
    geom_line()+
    scale_x_log10() +
    theme_classic() +
    labs(title = i,
         x = "Size (d.nm)",
         y = "Number (%)") +
    scale_y_continuous(breaks=seq(0,14,by=2), limits = c(0,14))
  print(gg)
  ggsave(filename = here('outputs', '01_DLS-data', file_name))
}

# Save summarized DLS data
save(data.summarized, file = here("outputs", "01_DLS-data", "01_DLS-data_summarized.RData"))

