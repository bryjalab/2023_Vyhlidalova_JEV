################################################################################
###################### ADDITION OF CONTROLS ####################################
################################################################################

# Libraries required:
library(here)
library(dplyr)

# Create output directory
dir.create(here('outputs', '09_controls-addition'))

# Load the data
load(here('outputs', '06_methods-intersection', '06_intersections.Rdata'))
load(here("outputs", "05_filter-out-B-samples", "05_data.filtered.Rdata")) # after B-proteins removal
patients <- read.csv(here('outputs', '06_methods-intersection', '06_upsetplot_background-table.csv'))
load(here('outputs', '03_data-processed', '03_data-processed.Rdata'))

# Add controls to the patients table
patients.ctrl <- left_join(patients, data %>% select(Suggested.Symbol, name, ID, starts_with("C")), by = c("value" = "Suggested.Symbol"))
patients.ctrl$X <- NULL
patients.ctrl <- patients.ctrl[, c(1, 13, 14, 2:12, 15:19)] # reorder columns, so IDs and gene names are in the beginning
patients.ctrl[,15:19][patients.ctrl[,15:19] >0] <- 1 # create binary matrix also from controls

# Export tables:
## Controls from the original table
ctrls <- data %>% select(Suggested.Symbol, name, ID, starts_with("C"))
write.csv(ctrls, here('outputs', '09_controls-addition', '09_controls-filtered-from-PG-table.csv'))
write.csv(patients.ctrl, here('outputs', '09_controls-addition', '09_controls-mapped-on-patients.csv'))

# Firstly, find proteins which are in all 11 patients, but no controls, or at least 9/11 patients, respectively
patients.11_ctrl.0 <- patients.ctrl %>%
  rowwise() %>%
  mutate(sum_patients = sum(c_across(patient_1: patient_11))) %>%
  mutate(sum_ctrl = sum(c_across(C1:C5)>0)) %>%
  filter(sum_ctrl == 0) %>%
  filter(sum_patients == 11)

patients.9_ctrl.0 <- patients.ctrl %>%
  rowwise() %>%
  mutate(sum_patients = sum(c_across(patient_1: patient_11))) %>%
  mutate(sum_ctrl = sum(c_across(C1:C5)>0)) %>%
  filter(sum_ctrl == 0) %>%
  filter(sum_patients > 8)

# Find the proteins, which are in the intersection (S&U) and also in U, but not in controls
SU_U_merge <- mapply(c, intersections, UC_specific, SIMPLIFY=FALSE)

SU_U <- data.frame(Suggested.Symbol = data$Suggested.Symbol)

for (i in c(1:11)){
  patient.name <- names(SU_U_merge)[i]
  patient <- data.frame(Suggested.Symbol = SU_U_merge[i], protein = 1)
  colnames(patient) <- c("Suggested.Symbol", patient.name)
  SU_U <- SU_U %>%
    left_join(., patient, by = c("Suggested.Symbol" = "Suggested.Symbol"))
}

SU_U[is.na(SU_U)] <- 0

SU_U <- left_join(SU_U, data %>% select(Suggested.Symbol, starts_with("C")), by = c('Suggested.Symbol' = "Suggested.Symbol"))
SU_U[,13:17][SU_U[,13:17] >0] <- 1

# 207 proteins! S&U+U proteins
SU_U_list <- SU_U %>%
  rowwise() %>%
  mutate(sum_patients = sum(c_across(patient_1: patient_11))) %>%
  mutate(sum_ctrl = sum(c_across(C1:C5)>0)) %>%
  filter(sum_ctrl == 0) %>%
  filter(sum_patients > 8)

# 157 proteins, they have to be also in at least one intersection
patient.specific <- SU_U_list[SU_U_list$Suggested.Symbol %in% patients$value, ]

# Export tables
write.csv(patients.9_ctrl.0, here('outputs', '09_controls-addition', '09_8-proteins_patients9-controls0.csv'))
write.csv(patient.specific, here('outputs', '09_controls-addition', '09_157-proteins_patient-specific.csv'))

save(patient.specific, patients.9_ctrl.0, file = here('outputs', '09_controls-addition', '09_patient-specific-data.Rdata'))
