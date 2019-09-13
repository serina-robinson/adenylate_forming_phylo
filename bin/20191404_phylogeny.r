# Install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("ggtree")
# BiocManager::install("treeio")
pacman::p_load("DECIPHER", "plot3D", "plot3Drgl", "ade4",
               "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", "ggtree", "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Read in the full length sequences
rawdat <- read_excel("../mibig_training_set_build_test/data/combined_adenylate_forming_training_set_for_db_20191404.xlsx")
rawdat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "NRPS"), collapse = "|"), rawdat$functional_class),] # 658 observations
dim(rawdat) # 627 sequences
sqnams_tr <- paste0(rawdat$acc, "_", rawdat$organism, "_", rawdat$small_substrate_group, "_", rawdat$functional_class)
which_mib <- rawdat$data_source == "mibig"
p2$data$label[p2$data$isTip] # 670
length(which_mib)
670-625
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

# Now align 
# dec_aln <- AlignSeqs(comb_dedup)
BrowseSeqs(dec_aln)

# Trim the boundaries
bound_tr <- AAStringSet(substr(dec_aln, start = 260, stop = 4030))
# BrowseSeqs(bound_tr)
# writeXStringSet(bound_tr, "output/671_seqs_for_phylogeny_aligned_20191404.fa")
checkr <- readAAStringSet("output/671_seqs_for_phylogeny_aligned_20191404.fa")
BrowseSeqs(checkr)
names(checkr)
table(word(names(checkr), sep = "_", -1))
table(word(names(checkr), sep = "_", -2))
grep("\\)", names(checkr))

# Now submit to FastTree for treeing
# nwk_tr <- treeio::read.newick("output/671_seqs_for_phylogeny_20191404_fasttree.nwk")

# nwk_tr <- read.table("output/671_seqs_for_phylogeny_20191404_fasttree.nwk", stringsAsFactors = F)
head(nwk_tr)
# nwk_try2 <- as.phylo(nwk_tr[,1])

# my_nwk<-newick2phylog(x.tre = nwk_tr[,1])

# nwk_phylo <- as.phylo(my_nwk)
# write.tree(nwk_phylo, "output/671_sqs_nwk_conversion.nwk")

phylo_fin <- treeio::read.newick("output/671_sqs_nwk_conversion.nwk")
phylo_find <- drop.tip(phylo_fin, "V9I5V9_9FABA_Caragana_korshinskii_cinnamate_and_succinylbenzoate_derivatives_ARYL")
phylo_fin$tip.label
# phylo_fin$node.label <- gsub(phylo_fin$node.label, "X", "")        
#phylo_fin$node.label <- paste0(word(phylo_fin$node.label, sep = "_", 1), ".", word(phylo_fin$node.label, sep = "_", 2))
head(phylo_fin$node.label)

# Make a phylogenetic tree
pl <- ggtree(phylo_find, layout = "circular")
tofix <- pl$data$label[!pl$data$isTip]
tofix <- gsub("X", "", tofix)    
tofix <- paste0(word(tofix, sep = "_", 1), ".", word(tofix, sep = "_", 2))
pl$data$label[!pl$data$isTip] <- tofix
tofix
table(word(pl$data$label[pl$data$isTip], -1, sep = "_")) # looks good

# Append metadata
dd <- data.frame(label = pl$data$label,
                 grp = ifelse(pl$data$isTip, word(pl$data$label, -1, sep = "_"), "NONE"), 
                 size = ifelse(as.numeric(pl$data$label) > 0.75, 0.1, 0))
as.numeric(pl$data$label)
p2 <- pl %<+% dd
tail(pl$data$label)

# Set the color palette
pal1
pal1 <- colorRampPalette(colors=brewer.pal(8, "Set1"))(8)
pal1[pal1 == "#377EB8"] <- "#92D050"
pal1[pal1 == "#A65628"] <- "gray68"
pal1[pal1 == "#F781BF"] <- "#A65628"
pal1[pal1 == "#4DAF4A"] <- "#377EB8"
pal1[pal1 == "#FFFF33"] <- "goldenrod"



pal2 <- c(pal1, "#F781BF", "blue1", "darkorchid1", "navy", #"black", 
          "gray68", "plum1", "blue1",
          "deepskyblue", "gold", "darkorchid1", 
          "deeppink2", "lightslateblue",
          "lightblue2", "darkseagreen1")
palette(pal2)
pal2

pdf("output/671_FastTree_ANL_seqs_unlabeled_circular_nodes_labeled.pdf", width = 15, height = 15)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  aes(color = grp) +
  scale_color_manual(values = pal2) +
  #ggplot2::xlim(-0.1, NA) +
  geom_nodepoint(size = p2$data$size[!p2$data$isTip], color = "gray40") +
  geom_tiplab2(aes(label = label), size = 0.5) +
  # geom_tiplab2(aes(label = ifelse(grepl("Xanthomonas", label), label, ""), size = 0.5)) +
  theme(legend.position = "right") +
  geom_treescale(offset = 1, fontsize = 4, x=3, y=1) +
  geom_tippoint(alpha = 0.7, size = ifelse(for_coloring_final, -1, 2), color = "gray40") +#aes(x = x+0.2), alpha = 0.7, size = ifelse(for_coloring_final, -1, 1.5), color = "gray40")
  
# geom_nodelab(label = p2$data$node[!p2$data$isTip], color = "black", size = 2) +
 geom_hilight_encircle(node = 727, fill = "#92D050", alpha = 0.2, linetype = 0)
  # geom_tiplab2(aes(label=label), size = 1)
  # ggplot2::xlim(-0.2, NA)
 #  geom_tiplab(aes(label=label), size=12, hjust=-.2)
#   geom_tiplab2(aes(label=label), size=12, hjust=-.2)
ptree
#scale_color_manual(pal)
# geom_text2(aes(subset=!isTip,label=label),size=14,hjust=-.2)
dev.off()



