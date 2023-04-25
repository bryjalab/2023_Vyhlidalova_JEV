################################################################################
##################### PROTEOMIC DATA CLEANING AND PROCESSING ###################
################################################################################

# Libraries required
library(here)
library(DEP)
library(dplyr)

# Create a new directory for processed data:
dir.create(here('outputs', '03_data-processed'))

# Proteomics data input (proteinGroups.txt table from MaxQuant)
data <- read.csv(here('outputs', '02_data-cleaning', 'updated_proteinGroups.csv'), stringsAsFactors = FALSE)

# Select columns containing intensities
data <- data %>%
  select(Suggested.Symbol, name, ID, starts_with("Intensity."))

# Rename the columns
colnames(data) <- gsub("Intensity.", "", colnames(data))

# Remove samples with chemo patients and patient 88
data <- data %>%
  select(-c("a71", "a54", "a55", "a65", "a88", "a88A", "a88B"))

# Remove qEV-only samples: 59, 74, 93
data <- data %>%
  select(-c("a59A", "a59B", "a74A", "a74B", "a93A", "a93B"))

# Remove UC-only samples: 33, 51, 62, 85
data <- data %>%
  select(-c("a33", "a51", "a62", "a85"))

# Rename controls
data <- data %>%
  rename(C1 = a42,
         C2 = a43a,
         C3 = a43b,
         C4 = a44c,
         C5 = a48,
  )

# Rename ultracentrifugation samples
data <- data %>%
  rename(U1 = a57,
         U2 = a77,
         U3 = a80,
         U4 = a87,
         U5 = a90,
         U6 = a92,
         U7 = a102,
         U8 = a107,
         U9 = a109,
         U10 = a111,
         U11= a120
  )

# Rename qEV samples
data <- data %>%
  rename(S1 = a57A,
         S2 = a77A,
         S3 = a80A,
         S4 = a87A,
         S5 = a90A,
         S6 = a92A,
         S7 = a102A,
         S8 = a107A,
         S9 = a109A,
         S10 = a111A,
         S11= a120A,
         B1 = a57B,
         B2 = a77B,
         B3 = a80B,
         B4 = a87B,
         B5 = a90B,
         B6 = a92B,
         B7 = a102B,
         B8 = a107B,
         B9 = a109B,
         B10 = a111B,
         B11= a120B
  )



# Save the cleaned data
save(data, file = here("outputs", "03_data-processed", "03_data-processed.RData"))
write.csv(data, file = here("outputs", "03_data-processed", "03_data-processed.csv"))

