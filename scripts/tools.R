#'
#' Match Yao 2014 supple_table3 - 67 DUX4 robust genes up-regulated in DUX4-target-positive 
#' FSHD biopsy samples, FSHD myotubes and DUX4-transduced cells
#' 1. match to GRCh38 (hg38)  Ensembl v. 88
#' 2. match to GRCh38 (hg38) Gencode v.36
#'
#'

library(tidyverse)
library(DESeq2)
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
# Originally we have 67 genes but remove: 
#   1. KHDC1L with id ENSG00000243501 was removed, this in Gencode 40 is called KHDC1-AS1
#   2. four have no identifier in GRCh38
# So we have 62 total left

targets <- readxl::read_xls(file.path(pkg_dir, "extdata", "ddu251supp_table3.xls"), 
                            skip = 1, sheet=1) %>%
  dplyr::select(ENSEMBL, gene.name) %>%
  dplyr::rename(gene_name = gene.name) %>%
  dplyr::distinct(gene_name, .keep_all=TRUE)

# match to ensembl v88; sanitized dds
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
load(file.path(pkg_dir, "data", "sanitized.dds.rda"))

annotation <- as.data.frame(rowData(sanitized.dds)) %>%
  dplyr::select(gene_id, gene_name) 
  #dplyr::mutate(ENSEMBL = sapply(str_split(gene_id, "\\."), "[[", 1))

# tools
get_id_by_name = function(g_name) {
  if (g_name=="no identifier") return(NA)
  annotation %>% dplyr::filter(gene_name %in% g_name) %>% pull(gene_id) 
}

# 1. find matched name
# 2. fix the names that are matched - search on-line information 
match_name = targets %>% left_join(annotation, by="gene_name", 
                                   suffix=c(".hg19", ".hg38")) %>%
    dplyr::mutate(match_name = !is.na(gene_id)) %>%
    dplyr::mutate(gene_name.hg38 = case_when(
    gene_name == "WI2-2994D6.2" ~ "PRAMEF25",
    gene_name == "PRAMEF3" ~ "PRAMEF13",
    gene_name == "RP13-221M14.1" ~ "HNRNPCL4",
    gene_name == "WI2-3308P17.2" ~ "PRAMEF11",
    gene_name == "RP13-221M14.5" ~ "PRAMEF9",
    gene_name == "RP13-221M14.3" ~ "no identifier",
    gene_name == "TRIM49DP" ~ "TRIM49D1",
    gene_name == "TRIM49L1" ~ "TRIM49D2",
    gene_name == "WI2-2994D6.1" ~ "no identifier", 
    gene_name == "XX-FW84067D5.1" ~ "no identifier",
    gene_name == "AC010606.1" ~ "no identifier",
    TRUE ~ gene_name
  )) %>% # 3. fixed gene_id and ENSEMBL.hg38
  dplyr::mutate(gene_id = if_else(is.na(gene_id), 
                                  map_chr(gene_name.hg38, get_id_by_name), 
                                  gene_id)) %>%
  dplyr::rename(ENSEMBL.hg19 = ENSEMBL, gene_name.hg19 = gene_name,
                gene_id.hg38 = gene_id) %>%
  dplyr::mutate(ENSEMBL.hg38 = sapply(str_split(gene_id.hg38, "\\."), "[[", 1)) 


######### search on-line ########
# WI2-2994D6.2: PRAMEF25; ENSG00000229571
# PRAMEF3: PRAMEF13; ENSG00000279169
# XX-FW84067D5.1: no current identifier
# RP13-221M14.1: HNRNPCL4; ENSG00000179412
# WI2-3308P17.2: PRAMEF11; ENSG00000239810
# RP13-221M14.5: PRAMEF9 ENSG00000204505
# RP13-221M14.3: no identifier -- but by gencode 40 we have PRAMEF35P ENSG00000236179
# WI2-2994D6.1: no current identifier
# TRIM49DP: TRIM49D1 ENSG00000223417
# TRIM49L1: TRIM49D2 ENSG00000233802
# AC010606.1: no current identifier, it is lincRNA, probably not MBD3L2B/CTB-25J19.1 ENSG00000196589


