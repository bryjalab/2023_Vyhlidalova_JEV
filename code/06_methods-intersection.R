################################################################################
############################## METHODS INTERSECTION ############################
################################################################################

# Libraries required:
library(here)
library(dplyr)
library(reshape2)
library(UpSetR)
library(ggplot2)

# Create a new directory for processed data:
dir.create(here('outputs', '06_methods-intersection'))

# Load the data
load(here("outputs", "05_filter-out-B-samples", "05_data.filtered.Rdata"))

# Function to create intersections for each patient between UC and SEC methods
make_intersections <- function (data, U, S, number) {
  require(dplyr)
  require(lazyeval)
  
  U <- enquo(U)
  S <- enquo(S)
  
  pU <- data %>%
    filter(!is.na(!!U)) %>%
    select(Suggested.Symbol)
  
  pS <- data %>%
    filter(!is.na(!!S)) %>%
    select(Suggested.Symbol)
  
  patient_tmp <- list( pU, pS)
  names_list <- c(paste("U", number, sep = ""),
                  paste("S", number, sep = ""))
  names(patient_tmp) <- names_list
  
  patients <- append(patients, list(patient_tmp))
}

# Create the intersections
patients <- list()

patients <- make_intersections(data.filtered, U1, S1, 1)
patients <- make_intersections(data.filtered, U2, S2, 2)
patients <- make_intersections(data.filtered, U3, S3, 3)
patients <- make_intersections(data.filtered, U4, S4, 4)
patients <- make_intersections(data.filtered, U5, S5, 5)
patients <- make_intersections(data.filtered, U6, S6, 6)
patients <- make_intersections(data.filtered, U7, S7, 7)
patients <- make_intersections(data.filtered, U8, S8, 8)
patients <- make_intersections(data.filtered, U9, S9, 9)
patients <- make_intersections(data.filtered, U10, S10, 10)
patients <- make_intersections(data.filtered, U11, S11, 11)

# Find the intersections:
intersections <- list(
  patient_1 = intersect(patients [[1]]$U1$Suggested.Symbol, patients[[1]]$S1$Suggested.Symbol),
  patient_2 = intersect(patients [[2]]$U2$Suggested.Symbol, patients[[2]]$S2$Suggested.Symbol),
  patient_3 = intersect(patients [[3]]$U3$Suggested.Symbol, patients[[3]]$S3$Suggested.Symbol),
  patient_4 = intersect(patients [[4]]$U4$Suggested.Symbol, patients[[4]]$S4$Suggested.Symbol),
  patient_5 = intersect(patients [[5]]$U5$Suggested.Symbol, patients[[5]]$S5$Suggested.Symbol),
  patient_6 = intersect(patients [[6]]$U6$Suggested.Symbol, patients[[6]]$S6$Suggested.Symbol),
  patient_7 = intersect(patients [[7]]$U7$Suggested.Symbol, patients[[7]]$S7$Suggested.Symbol),
  patient_8 = intersect(patients [[8]]$U8$Suggested.Symbol, patients[[8]]$S8$Suggested.Symbol),
  patient_9 = intersect(patients [[9]]$U9$Suggested.Symbol, patients[[9]]$S9$Suggested.Symbol),
  patient_10 = intersect(patients [[10]]$U10$Suggested.Symbol, patients[[10]]$S10$Suggested.Symbol),
  patient_11 = intersect(patients [[11]]$U11$Suggested.Symbol, patients[[11]]$S11$Suggested.Symbol))

# Which proteins are UC specific?
UC_specific <- list(
  patient_1 = setdiff(patients [[1]]$U1$Suggested.Symbol, patients[[1]]$S1$Suggested.Symbol),
  patient_2 = setdiff(patients [[2]]$U2$Suggested.Symbol, patients[[2]]$S2$Suggested.Symbol),
  patient_3 = setdiff(patients [[3]]$U3$Suggested.Symbol, patients[[3]]$S3$Suggested.Symbol),
  patient_4 = setdiff(patients [[4]]$U4$Suggested.Symbol, patients[[4]]$S4$Suggested.Symbol),
  patient_5 = setdiff(patients [[5]]$U5$Suggested.Symbol, patients[[5]]$S5$Suggested.Symbol),
  patient_6 = setdiff(patients [[6]]$U6$Suggested.Symbol, patients[[6]]$S6$Suggested.Symbol),
  patient_7 = setdiff(patients [[7]]$U7$Suggested.Symbol, patients[[7]]$S7$Suggested.Symbol),
  patient_8 = setdiff(patients [[8]]$U8$Suggested.Symbol, patients[[8]]$S8$Suggested.Symbol),
  patient_9 = setdiff(patients [[9]]$U9$Suggested.Symbol, patients[[9]]$S9$Suggested.Symbol),
  patient_10 = setdiff(patients [[10]]$U10$Suggested.Symbol, patients[[10]]$S10$Suggested.Symbol),
  patient_11 = setdiff(patients [[11]]$U11$Suggested.Symbol, patients[[11]]$S11$Suggested.Symbol))

# Which proteins are qEV specific?
SEC_specific <- list(
  patient_1 = setdiff(patients [[1]]$S1$Suggested.Symbol, patients[[1]]$U1$Suggested.Symbol),
  patient_2 = setdiff(patients [[2]]$S2$Suggested.Symbol, patients[[2]]$U2$Suggested.Symbol),
  patient_3 = setdiff(patients [[3]]$S3$Suggested.Symbol, patients[[3]]$U3$Suggested.Symbol),
  patient_4 = setdiff(patients [[4]]$S4$Suggested.Symbol, patients[[4]]$U4$Suggested.Symbol),
  patient_5 = setdiff(patients [[5]]$S5$Suggested.Symbol, patients[[5]]$U5$Suggested.Symbol),
  patient_6 = setdiff(patients [[6]]$S6$Suggested.Symbol, patients[[6]]$U6$Suggested.Symbol),
  patient_7 = setdiff(patients [[7]]$S7$Suggested.Symbol, patients[[7]]$U7$Suggested.Symbol),
  patient_8 = setdiff(patients [[8]]$S8$Suggested.Symbol, patients[[8]]$U8$Suggested.Symbol),
  patient_9 = setdiff(patients [[9]]$S9$Suggested.Symbol, patients[[9]]$U9$Suggested.Symbol),
  patient_10 = setdiff(patients [[10]]$S10$Suggested.Symbol, patients[[10]]$U10$Suggested.Symbol),
  patient_11 = setdiff(patients [[11]]$S11$Suggested.Symbol, patients[[11]]$U11$Suggested.Symbol))

# Create a dataframe with number of proteins for each category
## UC specific, SEC specific, in the intersection
df.intersections <- data.frame(S_U = lengths(intersections),
                               U = lengths(UC_specific),
                               S = lengths(SEC_specific))
df.intersections <- df.intersections[, c(3,1,2)]

# Plot the numbers for UC, SEC and intersection
svg(here('outputs', '06_methods-intersection', '06_barplot-methods-specific.svg'))
#pdf(here('outputs', '06_methods-intersection', '06_barplot-methods-specific.pdf'))
df.intersections %>%
  mutate(patient = rownames(df.intersections)) %>%
  select(patient, U, S, S_U) %>%
  pivot_longer(cols = U:S_U, names_to = "method", values_to = "protein_nr") %>%
  mutate(method_fact = factor(method, levels = c("U", "S_U", "S"))) %>%
  mutate(patient_fact = factor(patient, levels = c("patient_1", "patient_2", "patient_3", "patient_4",
                                                   "patient_5", "patient_6", "patient_7", "patient_8",
                                                   "patient_9", "patient_10", "patient_11"))) %>%
  ggplot(aes(x = patient_fact, y = protein_nr, fill = method_fact)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#1ECBE1", "#CBE11E", "#E11ECB" )) +
  coord_flip()+
  theme_minimal()+
  labs(x = "Patients",
       y = "Number of proteins",
       fill = "Method")
dev.off()


# Which proteins were identified by both methods in all patients
# Intersection of the intersections
intersection_of_intersections <- Reduce(intersect, intersections)

# Which proteins were identified by at least one method
# Union of the intersections
union_of_intersections <- Reduce(union, intersections)

# Which proteins are present in the intersection?
patients_intersection <- melt(intersections)
patients_intersection$present <- 1   #add variable that protein is present in intersection for particular patient

# Transform the long table to wide table (proteins = rows, patients = columns)
patients_intersection <- pivot_wider(patients_intersection, names_from = L1, values_from = present)
patients_intersection[is.na(patients_intersection)] <- 0


# Plot the upset plot of intersections
# Short version (default) - order by degree
svg(here('outputs', '06_methods-intersection', '06_upsetplot_intersections_short_degree.svg'))
#pdf(here('outputs', '06_methods-intersection', '06_upsetplot_intersections_short_degree.pdf'))
upset(fromList(intersections), order.by = "degree", nsets = 11)
dev.off()

# Short version (default) - order by freq
svg(here('outputs', '06_methods-intersection', '06_upsetplot_intersections_short_freq.svg'))
#pdf(here('outputs', '06_methods-intersection', '06_upsetplot_intersections_short_freq.pdf'))
upset(fromList(intersections), order.by = "freq", nsets = 11)
dev.off()

# Long version (nintersects = NA) - saved manually 6500x500
svg(here('outputs', '06_methods-intersection', '06_upsetplot_intersections_all-intersections_freq.svg'))
upset(fromList(intersections), order.by = "freq", nsets = 11, nintersects = NA)
dev.off()

# Save the outputs
write.csv(df.intersections, here('outputs', '06_methods-intersection', '06_proteins-methods_numbers.csv'))
write.csv(patients_intersection, here('outputs', '06_methods-intersection', '06_upsetplot_background-table.csv'))

save(intersections, UC_specific, intersection_of_intersections, union_of_intersections,
     file = here('outputs', '06_methods-intersection', '06_intersections.Rdata'))