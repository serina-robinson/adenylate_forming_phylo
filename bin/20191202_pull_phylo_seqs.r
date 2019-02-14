# Install packages
pacman::p_load("DECIPHER", "Biostrings", "readxl")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Read in the NRPS Adomain seqs from SANDPUMA2
# https://bitbucket.org/chevrm/sandpuma2/raw/10b732ca7ce501c2c326db6a662100d76d08ba71/flat/sp2.adomains.faa
adom <- readAAStringSet("data/sandpuma2_adomains.faa")
# 1,093 sequences with lengths ranging from 254 to 762 aa

# Read in the experimentally characterized OleC (beta-latone synthetase) sequences
olec <- readAAStringSet("data/17_experimentally_characterized_OleC_seqs.fasta")

