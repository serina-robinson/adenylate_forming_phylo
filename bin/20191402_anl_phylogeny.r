# Install packages
pacman::p_load("DECIPHER", "plot3D", "plot3Drgl",
               "phangorn", "RColorBrewer", "phylobase", "treeio", "ggtree", "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtree", version = "3.8")

BiocManager::install("treeio", version = "3.8")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Align sequences with HMMAlign
# hmmalign -o 338_aligned_seqs.fasta --outformat A2M ../AMP-binding.hmm 338_anl_seqs_for_treeing_cdhit90.fasta

# Reformat alignment
# t_coffee -other_pg seq_reformat -in 338_aligned_seqs.fasta -action +upper > output/338_aligned_seqs_upper.fasta

# Trim
# trimal -in 338_aligned_seqs_upper.fasta -gt 0.8 -out 338_aligned_seqs_upper_trimal.fasta -htmlout 338_aligned_seqs_upper_trimal.html

# FastTree
# FastTree <338_aligned_seqs_upper_trimal.fasta> 338_aligned_seqs_fasttree.nwk

# Read in the newick tree for ggtree

nwk_tr <- treeio::read.newick("data/338_aligned_seqs_fasttree.nwk")

# Make a phylogenetic tree
pl <- ggtree(nwk_tr)

# Append metadata
dd <- data.frame(label = pl$data$label,
                 grp = word(pl$data$label, 5, sep = "_"))
p2 <- pl %<+% dd

# Set the color palette
pal<-colorRampPalette(brewer.pal(12,"Paired"))(length(unique(p2$data$grp)))
pal <- c(pal, "black")


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
aa338 <- readAAStringSet("output/cdhit_90perc_clusters.fasta")

dec_aln <- AlignSeqs(aa338)
writeXStringSet(dec_aln, "output/338_aligned_seqs_DECIPHER.fasta")
phy<-read.phyDat("output/338_aligned_seqs_DECIPHER.fasta", format="fasta",type="AA")


# Make a distance matrix
dm<-dist.ml(phy)
rownames(dm)
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



pdf(paste0("output/", numseqs,"_PCoA.pdf"))
par(mar=c(0.01, 0.01, 0.01, 0.01))
ggplot(data = dat, aes(x=x, y=y)) + 
  geom_point() +
  # geom_point(size=sz,colour=clusts, alpha = ap, shape = shap2)+
  scale_shape(solid = TRUE)+
  labs(x=xlab,y=ylab) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.14, 0.80),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=20, face="bold", hjust=0)
  )
dev.off()

# Plot 3D
# scatter3D(x,y,z)
# plotrgl()
# rgl::plot3d(x,y,z)

############# 
# Do seqeunces cluster by their substrate?


# Color points for MDS 
source("src/color_pcoa_points.r")

# Remove certain points
rownames(cdat)
word(rownames(cdat), 5, sep = "_")
trdat <- cdat[!word(rownames(cdat), 5, sep = "_") == "OTHER",]
trdat <- trdat[!word(rownames(trdat), 5, sep = "_") == "HOLDOUTTEST",]
head(trdat)
numseqs <- length(trdat$x)

sub_tofind <- unique(word(rownames(trdat), 5, sep = "_"))
sub_tofind
tocol <- sub_tofind[!is.na(sub_tofind)]
tocol
pal <- palette(colorRampPalette(colors=brewer.pal(8,"Accent"))(length(tocol)))

cdat <- color_pcoa_points(tocol, dat, pal)
head(cdat)

numseqs
trdat$nms
# Make a colored 2D plot
pdf(paste0("output/", numseqs,"_PCoA_colored.pdf"))
par(mar=c(0.01, 0.01, 0.01, 0.01))
ggplot(data = trdat, aes(x=x, y=y)) + 
  geom_text(label = trdat$nms) +
  # geom_text(x = cdat$x, y = cdat$y, label = cdat$labl, check_overlap = T) +
  geom_point(fill = trdat$colrs, size = trdat$sz, 
             alpha = trdat$alph, shape = 21) +
  scale_color_manual(pal) +
  scale_shape(solid = TRUE) +
  labs(x=xlab,y=ylab) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.80, 0.80),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=20, face="bold", hjust=0)
  )
dev.off()

with(cdat, plot3d(x = x, y = y, z = z,
                  col = colrs, type = 's', size = 1))

tocol

pdf(file="324_ANL_legend.pdf",bg="white",width = 20)
plot.new()
legend("center", pch='.', ncol=1,
       legend= tocol, fill = pal, bty="n",
       text.col="black")
dev.off()
