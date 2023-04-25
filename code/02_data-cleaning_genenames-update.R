################################################################################
################ INITIAL PROCESSING OF DATA AND GENE NAMES UPDATE ##############
################################################################################

# Libraries required
library(HGNChelper)
library(here)
library(dplyr)
library(tidyr)
library(DEP)

# Create a new directory for clean data:
dir.create(here('outputs', '02_data-cleaning'))

# Get the input files in 'data' directory:
filenames <- list.files(here('data'), pattern = ".txt")
names <- gsub(".txt", "", filenames)
for (i in 1:length(names)) assign(names[i], read.delim(here('data', filenames[i])))

# Get the most actual gene names:
#genenames_newest <- getCurrentHumanMap()

# Using the gene names from 2022-08-16
load(here('data', 'genenames_update_20220816.RData'))

# General info how gene names are updated:
## if there is updated (suggested) ID, that one is used
## if there is NA value in suggested, it is replaced by original (although old ID)
####### done by 'unmapped.as.na = FALSE' parameter
## if there are more suggested IDs, delimited by "///", the one with higher GeneCards score is used

##################  proteinGroups TABLE ###########################

d1 <- proteinGroups

# Table of cRAP contaminants
cRAP.proteins <- d1[grepl("cRAP", d1$Majority.protein.IDs), ]
cRAP.proteins <- cRAP.proteins %>% select(Protein.IDs, Majority.protein.IDs, Gene.names, Protein.names, Fasta.headers, Peptides,
                                          starts_with("Peptides."))
write.csv(cRAP.proteins, here('outputs', '02_data-cleaning', 'proteinGroups_cRAP-filtered.csv'))

d1 <- d1 %>%
  filter(Reverse != "+") %>%
  filter(!grepl("cRAP", Majority.protein.IDs)) %>%
  filter(Only.identified.by.site != "+") %>%
  filter(!grepl("keratin", Fasta.headers)) %>%
  filter(!grepl("Keratin", Fasta.headers))

# Make unique identifiers for gene IDs and protein IDs
d1 <- make_unique(d1, "Gene.names", "Protein.IDs", delim = ";")

# Update gene names
update_d1 <- checkGeneSymbols(d1$name, species = "human", map = genenames_newest, unmapped.as.na = FALSE)

# Replace non-approved gene symbols
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "EPRS1 /// QARS1"] <- "QARS1"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "SEPTIN2 /// SEPTIN6"] <- "SEPTIN2"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "C12orf75 /// FSTL1"] <- "C12orf75"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "LNPK /// NUSAP1"] <- "LNPK"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "SARS1 /// SARS2"] <- "SARS1"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "MPHOSPH6 /// PALS2"] <- "PALS2"

update_d1$Suggested.Symbol[duplicated(update_d1$Suggested.Symbol)] # check for the duplicates

d1 <- left_join(d1, update_d1, by = c("name" = "x"))

write.csv(d1, here('outputs', '02_data-cleaning', 'updated_proteinGroups.csv'))


#####################  RNAseq_S2_table ##############################

d2 <- RNAseq_S2_table

# Update gene names
update_d2 <- checkGeneSymbols(d2$Cell.type, species = "human", map = genenames_newest, unmapped.as.na = FALSE)

update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "SARS1 /// SARS2"] <- "SARS1"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "SEPTIN2 /// SEPTIN6"] <- "SEPTIN2"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "EPRS1 /// QARS1"] <- "QARS1"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "GPAT3 /// LPCAT1"] <- "GPAT3"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MTARC1 /// MARCHF1"] <- "MARCHF1"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MPHOSPH6 /// PALS2"] <- "PALS2"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "TAFAZZIN /// WWTR1"] <- "TAFAZZIN"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "CCN3 /// PLXNA1 /// RPL10"] <- "CCN3"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "C11orf98 /// LBHD1"] <- "C11orf98"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "B3GNT2 /// B4GAT1"] <- "B4GAT1"  
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "BHLHE40 /// CENPX"] <- "BHLHE40"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "KAT14 /// PET117"] <- "KAT14"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "IRF4 /// PWWP3A"] <- "PWWP3A"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MTARC2 /// MARCHF2"] <- "MARCHF2"

# Duplicates solving
update_d2$Suggested.Symbol[duplicated(update_d2$Suggested.Symbol)] # check for the duplicates

d2 <- left_join(d2, update_d2, by = c("Cell.type" = "x"))

d2$Suggested.Symbol[duplicated(d2$Suggested.Symbol)] # check for the duplicates
# For these cases delete the older version of the protein
# SRSF10 is there 4x, so deleted
d2 <- d2[ !(d2$Cell.type %in% c("HNRNPU-AS1",  "CASC7",  "LSMD1", "STRA13")), ]
d2 <- d2[ !(d2$Cell.type %in% c("SRSF10")), ]

write.csv(d2, here('outputs', '02_data-cleaning', 'updated_RNAseq_S2_table.csv'))

#####################  Human Cell Map ##############################

d3 <- `preys-latest`

# Update gene names
update_d3 <- checkGeneSymbols(d3$symbol, species = "human", map = genenames_newest, unmapped.as.na = FALSE)
nrow(update_d3[update_d3$Approved == FALSE, ])

update_d3$Suggested.Symbol[update_d3$Suggested.Symbol == "SEPTIN2 /// SEPTIN6"] <- "SEPTIN2"
update_d3$Suggested.Symbol[update_d3$Suggested.Symbol == "EPRS1 /// QARS1"] <- "QARS1"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MPHOSPH6 /// PALS2"] <- "PALS2"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MT-CO2 /// PTGS2"] <- "MT-CO2"

d3 <- left_join(d3, update_d3, by = c("symbol" = "x"))

d3$Suggested.Symbol[duplicated(d3$Suggested.Symbol)] # duplicates

d3 <- d3[ !(d3$symbol %in% c("C10orf12", "KIAA0754")), ]

write.csv(d3, here('outputs', '02_data-cleaning', 'updated_preys-latest.csv'))


#################  MISEV protein categories #########################

d4 <- MISEV2018_protein_categories

# Transform to long-format table
colnames(d4) <- c("cat1", "cat2", "cat3", "cat4", "cat5")
d4 <- pivot_longer(d4, cat1:cat5, names_to = "category", values_to = "marker")
d4 <- d4 %>%
  arrange(category)
d4 <- d4[!d4$marker == "",]

# Update gene names
update_d4 <- checkGeneSymbols(d4$marker, species = "human", map = genenames_newest, unmapped.as.na = FALSE)
nrow(update_d4[update_d4$Approved == FALSE, ])
#update_d4$Suggested.Symbol <- sapply(strsplit(update_d4$Suggested.Symbol," /// "), `[`, 1)

d4 <- left_join(d4, update_d4, by = c("marker" = "x"))
d4$Suggested.Symbol[duplicated(d4$Suggested.Symbol)]

d4 <- d4[ !(d4$marker %in% c("IL27", "ITGA2B")), ] # ITGA2B is there 4x, so deleted

write.csv(d4, here('outputs', '02_data-cleaning', 'updated_MISEV_protein_categories.csv'))

#################  Izar (RNAseq) markers #########################

d5 <- Izar_cell_markers_updated

# Update gene names
update_d5 <- checkGeneSymbols(d5$Gene.Name, species = "human", map = genenames_newest, unmapped.as.na = FALSE)
nrow(update_d5[update_d5$Approved == FALSE, ])

d5 <- left_join(d5, update_d5, by = c("Gene.Name" = "x"))
d5$Suggested.Symbol[duplicated(d5$Suggested.Symbol)]

write.csv(d5, here('outputs', '02_data-cleaning', 'updated_Izar_cell_markers_updated.csv'))
