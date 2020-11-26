pacman::p_load("DECIPHER", "plot3D", "plot3Drgl", "ade4", "seqinr",
               "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", 
               "ggtree", "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Read in the phyogenetic tree
#phylo_fin <- treeio::read.newick("data/RAxML_bestTree.671_seqs_ANL_20190708_bootstrapped.nwk")
# phylo2 <- treeio::read.newick("data/RAxML_bootstrap.671_seqs_ANL_20190708_bootstrapped.nwk")
phylo3 <- treeio::read.newick("data/RAxML_bipartitions.671_seqs_ANL_20190708_bootstrapped.nwk")
  
# Make a phylogenetic tree
pl <- ggtree(phylo3, layout = "circular")

# Test draw
pdf("output/671_RaxML_raw.pdf", width = 15, height = 15)
pl
dev.off()

tofix <- pl$data$label[!pl$data$isTip]
table(word(pl$data$label[pl$data$isTip], -1, sep = "_")) # looks good

pl$data$label

# Append metadata
dd <- data.frame(label = pl$data$label,
                 grp = ifelse(pl$data$isTip, word(pl$data$label, -1, sep = "_"), "NONE"), 
                 size = ifelse(as.numeric(pl$data$label) > 75, 0.5, -1))
as.numeric(pl$data$label)
p2 <- pl %<+% dd


# Detect biosynthetic gene clusters
rawdat <- read_excel("../mibig_training_set_build_test/data/combined_adenylate_forming_training_set_for_db_20191404.xlsx")
rawdat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "NRPS"), collapse = "|"), rawdat$functional_class),] # 658 observations
dim(rawdat) # 627 sequences
sqnams_tr <- paste0(rawdat$acc, "_", rawdat$organism, "_", rawdat$small_substrate_group, "_", rawdat$functional_class)
which_mib <- rawdat$data_source == "mibig"
nnodes <- length(p2$data$label[!p2$data$isTip])
671-625
for_coloring_final <- c(which_mib, rep(TRUE, 46))

# for_coloring_final <- for_coloring[-grep("V9I5V9_9FABA", sqnams_tr)]


# Set the color palette
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

pdf("output/671_RaxML_bootstrapped_labeled.pdf", width = 12, height = 12)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  aes(color = grp) +
  scale_color_manual(values = pal2) +
  #ggplot2::xlim(-0.1, NA) +
  geom_nodepoint(size = p2$data$size[!p2$data$isTip], color = "gray40") +
  geom_tiplab2(aes(label = label), size = 0.5) +
  # geom_tiplab2(aes(label = ifelse(grepl("Xanthomonas", label), label, ""), size = 0.5)) +
  theme(legend.position = "right") +
  geom_treescale(offset = 1, fontsize = 4, x=3, y=1) #+
  #ggplot2::xlim(0.5, NA) 
  #geom_tippoint(alpha = 0.7, size = ifelse(for_coloring_final, -1, 2), color = "gray40")#aes(x = x+0.2), alpha = 0.7, size = ifelse(for_coloring_final, -1, 1.5), color = "gray40")
 
  # geom_nodelab(label = p2$data$node[!p2$data$isTip], color = "black", size = 2) +
  # geom_hilight_encircle(node = 727, fill = "#92D050", alpha = 0.2, linetype = 0)
# geom_tiplab2(aes(label=label), size = 1)
# 
#  geom_tiplab(aes(label=label), size=12, hjust=-.2)
#   geom_tiplab2(aes(label=label), size=12, hjust=-.2)
ptree
#scale_color_manual(pal)
# geom_text2(aes(subset=!isTip,label=label),size=14,hjust=-.2)
dev.off()

