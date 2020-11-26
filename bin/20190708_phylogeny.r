# Install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("ggtree")
# BiocManager::install("treeio")
pacman::p_load("DECIPHER", "plot3D", "plot3Drgl", "ade4", "seqinr",
               "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", "ggtree", "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Read in the full length sequences
rawdat <- read_excel("../mibig_training_set_build_test/data/combined_adenylate_forming_training_set_for_db_20191404.xlsx")
rawdat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "NRPS"), collapse = "|"), rawdat$functional_class),] # 658 observations
dim(rawdat) # 627 sequences
sqnams_tr <- paste0(rawdat$acc, "_", rawdat$organism, "_", rawdat$small_substrate_group, "_", rawdat$functional_class)
which_mib <- rawdat$data_source == "mibig"
for_coloring <- c(which_mib, rep(TRUE, 46))
for_coloring_final <- for_coloring[-grep("V9I5V9_9FABA", sqnams_tr)]


rawdat_sqs <- AAStringSet(rawdat$aa_seq)
names(rawdat_sqs) <- sqnams_tr
table(word(names(rawdat_sqs), -2, sep = "_"))

# Read in the NRPS sequences (one from each class)
nrps <- readAAStringSet("../mibig_training_set_build_test/data/sp2.adomainsD.faa")
monomers <- word(names(nrps), 2, sep = "\t")
monomers
names(nrps) <- paste0(names(nrps), "_", monomers, "_amino.acid_NRPS")
names(nrps) <- gsub("\\\t", "_", names(nrps))
names(nrps)

mondf <- data.frame(cbind(monomers, names(nrps)), stringsAsFactors = F) %>%
  dplyr::add_count(monomers) %>%
  dplyr::group_by(monomers) %>%
  dplyr::slice(1) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 2)
mondf # 49 observations with n > 2

nrps_to_include <- nrps[names(nrps) %in% mondf$V2]

# Combine everything and align
comb <- AAStringSet(c(rawdat_sqs, nrps_to_include))
length(comb) # 676 sequences
comb_dedup <- comb[!duplicated(comb)] # removed three duplicates
table(word(names(comb_dedup), sep = "_", -1))
table(word(names(comb_dedup), sep = "_", -2))

names(comb_dedup) <- gsub(paste0(c(" ", "\\(","\\)","\\+","\\:","\\?","\\,","\\/","\\;","\\[","\\]","\\-","\\="), collapse = "|"), "_", names(comb_dedup))
head(names(comb_dedup))
tail(names(comb_dedup))
writeXStringSet(comb_dedup, "output/671_seqs_for_phylogeny_unaligned_20191404.fa")

# Read in the HMMaligned results and remove all the lower case letters
hmmaligned <- readAAStringSet('data/671_seqs_HMMaligned_converted.fasta')
width(hmmaligned)
#hmmstr <- # sapply(hmmaligned, function(x) { splitseq(x, word = 1) })
hmmstr <- str_split(pattern = "", hmmaligned)

replace_lower <- function(x) {
  # Input is a vector of letters from a HMMaligned sequence
  # Output is a vector with lowercase letters replaced
  bool_lett <- x %in% letters
  x[bool_lett] <- "-"
  vec <- paste0(x, collapse = "")
  return(vec)
}

vecall <- lapply(hmmstr, replace_lower)
vecaa <- AAStringSet(unlist(vecall))
names(vecaa) <- names(hmmaligned)
#":", ",", ")", "(", ";", "]", "[", "'" 
names(vecaa) <- gsub("'|:| |,|\\(|\\)|\\[|\\]|;|#|-", "", names(vecaa))
names(vecaa)
writeXStringSet(vecaa, "output/671_seqs_HMMAligned_converted_uppercased.fasta")
