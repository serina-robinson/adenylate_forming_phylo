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

# nwk_tr <- treeio::read.newick("data/338_aligned_seqs_fasttree.nwk")
# 
# # Make a phylogenetic tree
# pl <- ggtree(nwk_tr)
# 
# # Append metadata
# dd <- data.frame(label = pl$data$label,
#                  grp = word(pl$data$label, 5, sep = "_"))
# p2 <- pl %<+% dd
# 
# # Set the color palette
# pal<-colorRampPalette(brewer.pal(12,"Paired"))(length(unique(p2$data$grp)))
# pal <- c(pal, "black")


# pdf("output/338_FastTree_ANL_seqs_nodes_labeled.pdf", width = 100, height = 100)
# par(mar=c(0.001,0.001,0.001,0.001))
# ptree <- p2 +
#   aes(color = grp) +
#   scale_color_manual(values = pal) +
#   ggplot2::xlim(-5, NA) +
#   geom_tiplab(aes(label=label), size=12, hjust=-.2)
# #   geom_tiplab2(aes(label=label), size=12, hjust=-.2)
# ptree
# #scale_color_manual(pal)
# # geom_text2(aes(subset=!isTip,label=label),size=14,hjust=-.2)
# dev.off()


# DistML
# aa338 <- readAAStringSet("output/cdhit_90perc_clusters.fasta")
aa459 <- readAAStringSet("output/anl_459_training_set_seqs.fasta")
names(aa459)

# dec_aln <- AlignSeqs(aa459)
writeXStringSet(dec_aln, "output/459_aligned_seqs_DECIPHER.fasta")
phy<-read.phyDat("output/459_aligned_seqs_DECIPHER.fasta", format="fasta",type="AA")


# Make a distance matrix
dm<-dist.ml(phy)
rownames(dm)
 #rownames(dm) <- names(phy)
# mat <- as.matrix(dm)
# write.csv(mat, paste0("data/",numseqs,"_distance_matrix.csv"), quote = FALSE, row.names = T)
# write.table(mat, paste0("data/",numseqs,"_distance_matrix.tsv"), quote = FALSE, row.names = T, sep = "\t")

# Multidimensional scaling using cmdscale
mds<-cmdscale(dm,eig=TRUE,k=3)
x<-mds$points[,1]
y<-mds$points[,2]
z<-mds$points[,3]

# Calculate percent explained by each principal component
pc1 <- mds$GOF[1]
pc2 <- (mds$GOF[2]-mds$GOF[1])
pc1
pc2

# Make x- and y-axis labels
xlab<-paste0("PC 1 (", round( pc1 * 100, 2), "% tot. explained var.)")
ylab<-paste0("PC 2 (", round( pc2 * 100, 2), "% tot. explained var.)")

numseqs <- length(x)
numseqs

# Make a 2D plot
dat <- data.frame(cbind(x, y))
head(dat)

# pdf(paste0("output/", numseqs,"_PCoA.pdf"))
# par(mar=c(0.01, 0.01, 0.01, 0.01))
# ggplot(data = dat, aes(x=x, y=y)) + 
#   geom_point() +
#   # geom_point(size=sz,colour=clusts, alpha = ap, shape = shap2)+
#   scale_shape(solid = TRUE)+
#   labs(x=xlab,y=ylab) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 14),
#         legend.key = element_rect(fill = "white"),
#         legend.background = element_rect(fill = "white"),
#         legend.position = c(0.14, 0.80),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.title=element_text(size=20, face="bold", hjust=0)
#   )
# dev.off()

# Plot 3D
# scatter3D(x,y,z)
# plotrgl()
# rgl::plot3d(x,y,z)

############# 
# Do seqeunces cluster by their substrate?


# Color points for MDS 
source("src/color_pcoa_points.r")
sub_tofind <- unique(word(rownames(dat), 5, sep = "_"))
sub_tofind
tocol <- sub_tofind[!sub_tofind %in% c("PEPTIDE", "HOLDOUTTEST", "OTHER")]
tocol
pal <- palette(colorRampPalette(colors=brewer.pal(8,"Accent"))(length(tocol)))
pal
cdat <- color_pcoa_points(tocol, dat, pal)


# Remove certain points
rownames(cdat)
word(rownames(cdat), 5, sep = "_")
trdat <- cdat[!word(rownames(cdat), 5, sep = "_") == "OTHER",]
trdat <- trdat[!word(rownames(trdat), 5, sep = "_") == "HOLDOUTTEST",]
trdat <- trdat[!word(rownames(trdat), 5, sep = "_") == "PEPTIDE",]
head(trdat)
numseqs <- length(trdat$x)
numseqs

# Make a colored 2D plot
pdf(paste0("output/", numseqs,"_PCoA_colored.pdf"))
par(mar=c(0.01, 0.01, 0.01, 0.01))
ggplot(data = trdat, aes(x=x, y=y)) + 
  geom_text(label = trdat$nms) +
  # geom_text(x = cdat$x, y = cdat$y, label = cdat$labl, check_overlap = T) +
  geom_point(fill = trdat$colrs, size = trdat$sz, 
             alpha = trdat$alph, shape = 21) +
  scale_color_manual(unique(trdat$colrs)) +
  scale_shape(solid = TRUE) +
  labs(x=xlab,y=ylab) +
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = "bottom",
        legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=20, face="bold", hjust=0))
dev.off()

with(cdat, plot3d(x = x, y = y, z = z,
                  col = colrs, type = 's', size = 1))




pdf(file="426_ANL_legend.pdf",bg="white",width = 20)
plot.new()
legend("center", pch='.', ncol=1,
       legend= unique(trdat$nms), fill = unique(trdat$colrs), bty="n",
       text.col="black")
dev.off()

length(tocol)
tocol

pllist <- list()
for(i in 1:length(tocol)) {
  trdat_one <- trdat[trdat$nms == tocol[i],]
 
  pllist[[i]] <- ggplot(data = trdat_one, aes(x=x, y=y)) + 
    geom_text(label = rownames(trdat_one)) +
    # geom_text(x = cdat$x, y = cdat$y, label = cdat$labl, check_overlap = T) +
    geom_point(fill = trdat_one$colrs, size = trdat_one$sz, 
               alpha = trdat_one$alph, shape = 21) +
    scale_color_manual(unique(trdat_one$colrs)) +
    scale_shape(solid = TRUE) +
    labs(x=xlab,y=ylab) +
    theme_bw() + 
    theme(axis.text = element_text(size = 14),
          legend.key = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          legend.position = "bottom",
          legend.box = "horizontal",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title=element_text(size=20, face="bold", hjust=0))
}

pdf(paste0("output/Plot_grid_PCoA_colored.pdf"), width = 35, height=  20)
par(mar=c(0.01, 0.01, 0.01, 0.01))
plot_grid(plotlist=pllist)
dev.off()
