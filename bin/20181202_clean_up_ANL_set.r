# Install packages
pacman::p_load("DECIPHER", "Biostrings", "readxl", "tidyverse", "gtools")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/adenylate_forming_phylo/")

# Read in the ANL training set
anl <- read_excel("data/20191202_ANL_training_set.xlsx") %>%
  dplyr::filter(substrate != "NRPS") 
dim(anl)

pmids <- anl$pub_med_id.x %>%
  str_split(., pattern = ";") %>%
  map(., 1) %>%
  unlist() %>%
  unique()

pmid_no_missing <- pmids[pmids!="NA"]
write.table(pmid_no_missing, "output/all_pmids.txt", quote = F, col.names = F, row.names = F)



# Read in the PubMed result 
pubmed_raw <- read_csv("data/pubmed_result.csv")
pubmed_gs <- gsub("create date:", "", pubmed_raw$Properties)
pubmed_spl <- str_split(pubmed_gs, "/")
yr <- sapply(pubmed_spl, head, 1)
dy <- substr(sapply(pubmed_spl, tail, 1), 1, 2)
mon <- pubmed_spl %>%
  map_chr(c(2)) %>%
  as.character()


pubmed_res <- pubmed_raw %>%
  rename(title = Title,
         jabbrv = ShortDetails,
         pmid_first = EntrezUID) %>%
  mutate(year = yr,
         month = mon,
         day = dy) %>%
  dplyr::select(pmid_first, year, month, day, jabbrv, title)


# Merge with the existing ANL set 
anl_no_pmid <- anl %>%
  dplyr::filter(title == "NA") %>%
  mutate(pmid_first  = as.numeric(pmid_first))
#  dplyr::filter(grepl(paste0(pubmed_res$EntrezUID, collapse = "|"), pmid_first))

pubmed_res$pmid_first %in% anl_no_pmid$pmid_first

anl_merge <- anl_no_pmid %>%
  dplyr::select(-jabbrv, -title, -year, -day, -month) %>%
  dplyr::left_join(., pubmed_res, by = "pmid_first")

# Merge with the original set 
anl_pmid_known <- anl %>%
  dplyr::filter(title != "NA") %>%
  mutate(pmid_first = as.numeric(pmid_first))


final <- rbind(anl_pmid_known, anl_merge)

table(final$abstract)
write_csv(final, "data/anl_training_set_merged.csv")

# # anl_final <- anl_pmid_known %>%
#   full_join(., anl_merge)
# 
# 
# 
# write_csv(anl_final, "output/anl_merge_pubmed.csv")
