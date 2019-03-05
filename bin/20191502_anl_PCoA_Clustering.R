# Install packages
pacman::p_load("DECIPHER", "plot3D", "plot3Drgl",
               "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", "ggtree", "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtree", version = "3.8")

# BiocManager::install("treeio", version = "3.8")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Align sequences with DECIPHER

# Trim
# trimal -in 338_aligned_seqs_DECIPHER.fasta -gt 0.8 -out 338_aligned_seqs_DECIPHER_trimal.fasta -htmlout 338_aligned_seqs_upper_trimal.html

# FastTree
# FastTree <338_aligned_seqs_upper_trimal.fasta> 338_aligned_seqs_fasttree.nwk

# Read in the newick tree for ggtree

# nwk_tr <- treeio::read.newick("output/338_aligned_seqs_DECIPHER_Fasttree.nwk")

nwk_ape <- ape::read.tree("output/338_aligned_seqs_DECIPHER_Fasttree.nwk")
todrop <- nwk_ape$tip.label[grep("HOLDOUTTEST|PEPTIDE|OTHER|unknown|131|89|153|164", nwk_ape$tip.label)]
nwk_drop <- drop.tip(nwk_ape, todrop)

# Make a phylogenetic tree
pl <- ggtree(nwk_drop, layout = "circular", size = 2)

# Append metadata
dd <- data.frame(label = pl$data$label,
                 grp = word(pl$data$label, 5, sep = "_"))
p2 <- pl %<+% dd

# Set the color palette
pal<-colorRampPalette(brewer.pal(8,"Set2"))(length(unique(p2$data$grp)))
# pal <- c(pal, "black")


p2$data$grp <- as.character(p2$data$grp)
p2$data$grp[is.na(p2$data$grp)] <- "black"
p2$data$grp <- as.factor(p2$data$grp)
grpcol <- pal[p2$data$grp]
grpcol[grpcol == "#C5A07A"] <- "black"
grpcol
pal
pdf("output/20191802_FastTree_ANL_seqs.pdf", width = 30, height = 30)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  scale_color_manual(values = pal) +
  aes(color = grpcol) +
  ggplot2::xlim(0, NA) +
  # geom_nodelab2(aes(label = node), size = 12) +
  # geom_tiplab2(size = 12) +
  geom_hilight(node = 485, fill = "#FCD738", alpha = 0.2) +
  geom_hilight(node = 385, fill = "#B1C968", alpha = 0.2) +
  geom_hilight(node = 474, fill = "#BE93C6", alpha = 0.2) +
  geom_hilight(node = 371, fill = "#D2BD9F", alpha = 0.2) +
  geom_hilight(node = 557, fill = "#DB98AE", alpha = 0.2) +
  geom_hilight(node = 436, fill = "#DD927E", alpha = 0.2)
 # geom_tiplab(aes(label=label), color = grpcol, size=12, hjust=-.2)
#   geom_tiplab2(aes(label=label), size=12, hjust=-.2)
ptree
#scale_color_manual(pal)
# geom_text2(aes(subset=!isTip,label=label),size=14,hjust=-.2)
dev.off()

show_col(pal)
cbind(pal, txt)
txt <- c("ARYL", "CAR", "LACS", "VLACSBILE", "FAAL", "LUCIFERASE", "SACS", "FAT", "BLS", "MMCS", "MACS")

pdf(file="FastTree_ANL_legend.pdf",bg="white",width = 20)
jpeg(file = "FastTree_ANL_legend.jpg", bg = "white")
plot.new()
legend("center", pch='.', ncol=1,
       legend= txt, 
       fill = pal[1:length(pal)-1], bty="n",
       text.col="black")
dev.off()


