# Install packages
pacman::p_load("DECIPHER", "Biostrings", "readxl", "tidyverse", "gtools")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Read in the ANL set
train <- read_excel("data/anl_training_set_merged_current.xlsx")

# Read in the UniPROT download
uni <- fread("data/uniprot-yourlist_M201902136746803381A1F0E0DB47453E0216320D25CA225.tab", data.table = F) %>%
  janitor::clean_names() %>%
  dplyr::select(entry_name,  pathway, sequence, kinetics, function_cc)

uni$pathway

# Merge
merg <- train %>%
  rename(pdb_seq = sequence) %>%
  inner_join(., uni, by = "entry_name")
dim(merg)

write_csv(merg, "data/anl_training_set_20190213_current.csv")
