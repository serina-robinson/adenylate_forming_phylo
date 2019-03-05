# Install packages
pacman::p_load("DECIPHER", "Biostrings", "data.table", "readxl", "tidyverse", "gtools", "rentrez")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Read in the tbl
# tbl <- read_excel("data/anl_training_set_20190213_current.xlsx")
# tbl_nodup <- tbl[!duplicated(tbl$entry_name),]
# write_csv(tbl_nodup, "data/anl_training_set_20190213_current_nodups.csv")
# 
# # Read in the updated table
# table(tbl_new$substrate_group)
# table(tbl_new$class)
# tbl_new$class[tbl_new$class == "UNKNOWN"] <- NA
# tbl_new$class[grep("fatty acid transport", tbl_new$protein_names)] <- "FAT"
# tbl_new$substrate[grep("fatty acid transport", tbl_new$protein_names)] <- "long_chain"
# tbl_new$substrate_group[grep("fatty acid transport", tbl_new$protein_names)] <- "long_chain"
# tbl_new$substrate_group[tbl_new$substrate_group == "bile_acid"] <- "very_long_chain_bile"
# tbl_new$substrate_group[tbl_new$substrate_group == "very_long_chain"] <- "very_long_chain_bile"
# tbl_new$class[tbl_new$class == "ALA"] <- "PEPTIDE"
# tbl_new$protein_names[tbl_new$class == "PEPTIDE"]
# table(tbl_new$class)
# table(tbl_new$confidence_score)
# tbl_new$confidence_score[tbl_new$confidence_score == 1] <- 2
# tbl_new$confidence_score[tbl_new$confidence_score == 5] <- 2

# write_csv(tbl_new, "data/anl_training_set_20190214_current.csv")

# Trim the training set
train <- read_excel("data/anl_training_set_20190214_current.xlsx") %>%
  dplyr::select(entry, gene_names, entry_name, organism, length, structure_title,
                protein_names, substrate, substrate_group, class,
                function_cc, kinetics, pdb_id.x, 
                structure_title, pub_med_id.x, 
                doi, year, month, day, jabbrv, title, abstract, ec_numbers,
                sequence) %>%
  dplyr:: rename(pub_med_id = pub_med_id.x, 
         pdb_id = pdb_id.x,
         functional_class = class,
         aa_seq = sequence) %>%
  dplyr::filter(!functional_class == "NA")

train$functional_class[train$functional_class == "PDB_TEST"] <- "HOLDOUT_TEST"

# 386 sequences

# Add in the OleC sequences
svtn <- readAAStringSet("data/17_experimentally_characterized_OleC_seqs.fasta")
sxty <- readAAStringSet("data/68_sukovich_oleCs.fasta")
orfs <- readAAStringSet("data/lstC_orf1.fasta")
orfs <- orfs[-grep("toxytricini", names(orfs))] # Remove virginiae
comb <- c(orfs, svtn, sxty)
comb_aa <- AAStringSet(comb[!duplicated(comb)])
length(comb_aa) # 74 sequences

# Extract the organism name
steno <- names(comb_aa)[names(comb_aa) == "AFC01244.1 Sequence 6 from patent US 8110093"]
names(comb_aa)[names(comb_aa) == "AFC01244.1 Sequence 6 from patent US 8110093"] <- paste0(steno, "[Stenotrophomonas maltophilia sequence 6 from patent US 8110093]")
orgnams <- unlist(str_extract_all(names(comb_aa), "\\[[^()]+\\]"))
orgnams_clnd <- substring(orgnams, 2, nchar(orgnams)-1)
names(comb_aa)
names(comb_aa) <- trimws(gsub("OleC", "", names(comb_aa)))


# substrates <- rep("median_beta_hydroxy_long", length(orgnams_clnd))
orgnams_clnd
length(orgnams_clnd)

bltibb <- data.frame(organism = orgnams_clnd,
                 length = width(comb_aa),
                 protein_names = names(comb_aa),
                 substrate = "median_beta_hydroxy_acid_long",
                 substrate_group = "median_beta_hydroxy_acid_long",
                 functional_class = "BLS",
                 ec_numbers = "EC 6.1.3.1",
                 confidence_score = 1,
                 aa_seq = as.character(comb_aa), stringsAsFactors = F)

bltibb$functional_class[grep("Nocardia", bltibb$organism)] <- "HOLDOUT_TEST"

train_aug <- train %>%
  full_join(., bltibb)
write_csv(train_aug, "data/anl_training_set_with_olec_20190214.csv")

table(train_aug$substrate_group)

# Lucs from Gimenez et al. squid paper
# Candidate luciferases from bioluminescent creatures: BAE80728
# Enzymes that produce light from non-luminescent creatures: NP_651221

# Add in the squid coelanterazine sequences
squid <- read_excel("data/Watasenia_squid_luciferase.xlsx") %>%
  mutate(confidence_score = 2)

colnames(squid)
colnames(train_aug)
train_aug2 <- train_aug %>%
  full_join(., squid)
tran
write_csv(train_aug2, "data/anl_training_set_with_olec_and_squid_20190214.csv")

# Add in the cyanobacterial FAAL sequences?

  
# Add in the luciferases
blast_hits <- fread("data/Luciferase_PSI_BLAST_HitTable.csv", data.table = F)
accs <- blast_hits$V2[blast_hits$V13 > 40] # 524 sequences
# write.table(accs, "data/luciferase_hits_forty_percent_accs.txt", row.names = F, col.names = F, quote = F)
# accs_prot <- entrez_fetch(db="protein", id = accs[1:200], rettype = "fasta", ap_key="826357f5ff17c7ec62e583909071e94f9d08")
fa <- readAAStringSet("data/top_500_luciferase_seqs.fasta")
luc_seqs <- fa[grep("lucifer", names(fa))]
luc_seqs <- luc_seqs[!duplicated(luc_seqs)]
luc_no_syn <- luc_seqs[-grep("synthetic|construct|Cloning|vector|Aedes|partial|Trichoplusia|Harpegnathos|Red-bioluminescence|replicon", names(luc_seqs))]
# names(luc_no_syn)[names(luc_no_syn]
orgnams <- unlist(str_extract_all(names(luc_no_syn), "\\[[^()]+\\]"))
str_extract_all(names(luc_no_syn), "\\[[^()]+\\]")
orgnams
orgnams_clnd <- substring(orgnams, 2, nchar(orgnams)-1)
luc_one_per_genus <- AAStringSet(luc_no_syn[!duplicated(orgnams_clnd)])
orgnams_nodups <- orgnams_clnd[!duplicated(orgnams_clnd)]

names(luc_one_per_genus)
length(orgnams_nodups)

luc_df <- data.frame(organism = orgnams_nodups,
                     length = width(luc_one_per_genus),
                     protein_names = names(luc_one_per_genus),
                     substrate = "luciferin",
                     substrate_group = "luciferin",
                     functional_class = "LUCIFERASE",
                     ec_numbers = "EC 1.13.12.",
                     confidence_score = 1,
                     aa_seq = as.character(luc_one_per_genus), stringsAsFactors = F)

head(luc_df)

train_aug3 <- train_aug2 %>%
  full_join(., luc_df)

# write_csv(train_aug3, "data/anl_training_set_20190214_olec_squid_luciferase.csv")

#write_csv(df2[duplicated(df2$aa_seq),], "data/df2_seq_duplicates.csv")

# Cross reference with Swiss-Prot luciferases
# luc <- fread("data/uniprot-luciferase-filtered-reviewed_yes_300.txt", data.table = F) %>%
#   janitor::clean_names()
# luc_aa <- AAStringSet(luc$sequence)
# names(luc_aa) <- paste0(luc$entry_name, "_", luc$protein_names)
# # luc_aa <- luc_aa[!duplicated(luc_aa)]
# # luc_comb <- AAStringSet(c(luc_one_per_genus, luc_aa))
# # luc_dup <- names(luc_comb)[duplicated(luc_comb)]
# # luc_dup
# tok <- luc[grep("Luciferin 4-monooxygenase", luc$protein_names),]
# tok2 <- tok[!tok$entry_name %in% train_aug$entry_name,]
# luc_df[grep("lateralis|pensylvanica", luc_df$organism),]

# Add in the JGI sequences (These will be test holdout sequences)
# As will Orf1, LstC, and the Nocardia sequence

jgi <- readAAStringSet("data/5_JGI_OleC_seqs.fasta")
orgnams <- unlist(str_extract_all(names(jgi), "\\[[^()]+\\]"))
orgnams
orgnams_clnd <- substring(orgnams, 2, nchar(orgnams)-1)
orgnams_clnd

jgi_tibb <- data.frame(organism = orgnams_clnd,
                     length = width(comb_aa),
                     protein_names = names(comb_aa),
                     substrate = "median_beta_hydroxy_acid_long",
                     substrate_group = "median_beta_hydroxy_acid_long",
                     functional_class = "BLS",
                     ec_numbers = "6.1.3.1",
                     confidence_score = 1,
                     aa_seq = as.character(comb_aa), stringsAsFactors = F)


# Fix the duplicates problem

# CmiS6 was duplicated
df2 <- read_csv("data/anl_training_set_20190214_olec_squid_luciferase.csv")
cmis6 <- readAAStringSet("data/cmis6.fasta")
df2$aa_seq[grep("cmiS6", df2$gene_names)] <- as.character(cmis6)

df_nodups <- df2[!duplicated(df2$aa_seq),]
dim(df_nodups) # 511 final sequences
table(df_nodups$substrate)

# Remove any sequences that are shorter than 200 amino acids or longer than 1500
df3 <- df_nodups[nchar(df_nodups$aa_seq) < 1500,]
df4 <- df3[nchar(df3$aa_seq) > 200,]
df_nofrag <- df4[-grep("Fragment", df4$protein_names),] #459 seqs
table(df_nofrag$substrate_group)
# df_nofrag$protein_names[is.na(df_nofrag$organism)]

# write_csv(df_nofrag, "data/final_anl_training_set_20190214.csv")

df_nofrag$organism[grep("Wat", df_nofrag$organism)] <- "Watasenia scintillans"


df_nofrag <- df_nofrag[-grep("Bifunctional", df_nofrag$protein_names),] 
df_nofrag$functional_class[df_nofrag$protein_names == "Uncharacterized protein, isoform E (EC 6.2.1.3)" & df_nofrag$organism == "Drosophila melanogaster (Fruit fly)"] <- "HOLDOUT_TEST"
faals <- df_nofrag[df_nofrag$functional_class == "FAAL",]
faal2 <- faals[grep("CoA", faals$protein_names),]
df_nofrag$functional_class[df_nofrag$entry_name %in% faal2$entry_name] <- "LACS"

# write_csv(df_nofrag, "output/anl_training_set_updated_20190215.csv") 

df_for_phylo <- df_nofrag %>%
  mutate(substrate_group = str_replace_all(substrate_group, pattern = "_|' '", replacement =  "")) %>%
  mutate(org_short = word(organism, 1,2, sep=" ")) %>%
  mutate(org_short = str_replace_all(org_short, pattern = " ", replacement =  "_")) %>%
  mutate(org_short = str_replace_all(org_short, pattern = "\\.", replacement =  "")) %>%
  mutate(substrate = str_replace_all(substrate, pattern = "_|' '", replacement =  "")) %>%
  mutate(substrate = str_replace_all(substrate, pattern = "[[:punct:]]", replacement =  "")) %>%
  mutate(substrate = str_replace_all(substrate, pattern = " ", replacement =  "")) %>%
  mutate(functional_class = str_replace_all(functional_class, pattern = "_|' '", replacement =  "")) %>%
 #  mutate(entry_name = str_replace_all(entry_name, pattern = "_", replacement =  "")) %>%
  mutate(sqnams = paste0(1:nrow(df_nofrag), "_", org_short, "_", substrate, "_", functional_class)) 

write_csv(df_for_phylo, "output/anl_training_set_updated_20190215_fixnams.csv")

# Sequences to exclude from training because unclear
torem <- c(131, 89, 153, 164)
df_for_phylo[46,] #bifuncitonal AAs protein



train459_aa <- AAStringSet(df_for_phylo$aa_seq)

names(train459_aa) <- df_for_phylo$sqnams


# Add in sequences from the 4-chlorobenzoate
chloro <- fread("data/uniprot-_4+coumarate+coa+ligase_+existence__evidence+at+transcript+level--.tab", data.table = F)
head(chloro)
table(chloro$Organism)


# summary(width(train459_aa))
# writeXStringSet(train459_aa, "output/anl_459_training_set_seqs.fasta")


# Now try adding in the carnitine/crotonobetaine sequences
# carn <- fread("data/uniprot-_crotonobetaine+carnitine+coa_+ligase-filtered-reviewed_yes.tab", data.table = F)
# carn$`Protein existence`

# Read in CD-HIT results at 95% sequence id

cdhit95 <- readAAStringSet("output/cdhit_95perc_clusters.fasta")
grps <- word(names(cdhit95), 4, sep="_")
names(cdhit95)[grep("coelenterazine", names(cdhit95))] <- "Watasenia_scintillans_coelenterazinedisulfate_OTHER"
table(grps)
cdhit97 <- readAAStringSet("output/cdhit_97perc_clusters.fasta")
grps <- word(names(cdhit97), 4, sep = "_")
table(grps)
orgs <- word(names(cdhit97), 1, sep = "_")
orgs
names(cdhit97)[grep("coelenterazine", names(cdhit97))] <- "Watasenia_scintillans_coelenterazinedisulfate_OTHER"
cdhit90 <- readAAStringSet("output/cdhit_90perc_clusters.fasta")
grps <- word(names(cdhit90), 4, sep = "_")
table(grps)

names(cdhit90)[grep("coelenterazine", names(cdhit90))] <- "Watasenia_scintillans_coelenterazinedisulfate_OTHER"


# We'll do it at 90% for treeing?
length(cdhit90) # 338 sequences
writeXStringSet(cdhit90, "data/338_anl_seqs_for_treeing_cdhit90.fasta")


