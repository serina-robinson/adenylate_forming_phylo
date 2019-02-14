# Install packages
pacman::p_load("DECIPHER", "Biostrings", "readxl", "tidyverse", "gtools")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Read in the tbl
tbl <- read_excel("data/anl_training_set_20190213_current.xlsx")
tbl_nodup <- tbl[!duplicated(tbl$entry_name),]
write_csv(tbl_nodup, "data/anl_training_set_20190213_current_nodups.csv")

# Read in the updated table
table(tbl_new$substrate_group)
table(tbl_new$class)
tbl_new$class[tbl_new$class == "UNKNOWN"] <- NA
tbl_new$class[grep("fatty acid transport", tbl_new$protein_names)] <- "FAT"
tbl_new$substrate[grep("fatty acid transport", tbl_new$protein_names)] <- "long_chain"
tbl_new$substrate_group[grep("fatty acid transport", tbl_new$protein_names)] <- "long_chain"
tbl_new$substrate_group[tbl_new$substrate_group == "bile_acid"] <- "very_long_chain_bile"
tbl_new$substrate_group[tbl_new$substrate_group == "very_long_chain"] <- "very_long_chain_bile"
tbl_new$class[tbl_new$class == "ALA"] <- "PEPTIDE"
tbl_new$protein_names[tbl_new$class == "PEPTIDE"]
table(tbl_new$class)
table(tbl_new$confidence_score)
tbl_new$confidence_score[tbl_new$confidence_score == 1] <- 2
tbl_new$confidence_score[tbl_new$confidence_score == 5] <- 2

# write_csv(tbl_new, "data/anl_training_set_20190214_current.csv")

# Get some more luciferase sequences


# Trim the training set

