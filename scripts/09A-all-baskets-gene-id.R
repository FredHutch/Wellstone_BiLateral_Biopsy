# make data frame to include the gencode id of all baskets

library(tidyverse)
library(DESeq2)

# get annotation from gencode36
load(file.path(pkg_dir, "data", "dds.rda"))
anno_gencode36 <- anno_gencode36 <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id") %>%
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", ""))

#
# Inflamm, ECM, complement, and IG baskets :gene_id, ens_id, gene_name in gencode v35
#
load(file.path(pkg_dir, "data", "FSHD_signatures_baskets.rda"))
load(file.path(pkg_dir, "data", "candidates.rda"))
load(file.path(pkg_dir, "data", "complement_df.rda")) # 21 genes from four complement activation pathway
FSHD_signatures_baskets$inflamm <- FSHD_signatures_baskets$inflamm %>% 
  add_row(FSHD_signatures_baskets$stress %>% dplyr::filter(gene_name == "CDKN1A")) %>%
  dplyr::filter(! gene_name %in% c("SPP1", "THBS1", "DMBT1", "TREM2"))

FSHD_signatures_baskets$ecm <- FSHD_signatures_baskets$ecm %>% 
  dplyr::filter(!gene_name == "SPP1")

# these is the basket id of ensembl88
baskets_id <- list(`DUX4-M` = candidates %>% dplyr::filter(`basket-M`) %>% pull(gene_id),
                   `DUX4-M6` = candidates %>% dplyr::filter(`basket-M6`) %>% pull(gene_id),
                   `DUX4-M12` =  candidates %>% dplyr::filter(`basket-M12`) %>% pull(gene_id),
                   ECM  = FSHD_signatures_baskets$ecm$gene_id,
                   Inflamm = FSHD_signatures_baskets$inflamm$gene_id,
                   Complement = complement_df %>% dplyr::filter(gene_name %in% c("C1QA", "C1QB", "C1QC", "C1R", "C1S", "C3")) %>% pull(gencode_id),
                   IG = complement_df %>% dplyr::filter(gene_name %in% c("IGHG2", "IGHG1", "IGHG3", "IGHG4", "IGKC", "FCGR2B")) %>% pull(gencode_id))

# exclude

all_baskets <- map(baskets_id, function(id) {
    data.frame(gene_id_v88 = id) %>%
      dplyr::mutate(ens_id =str_replace(gene_id_v88, "\\..*", "")) %>%
      dplyr::left_join(anno_gencode36, by="ens_id") %>%
      rename(gencode_v35 = "gene_id")
}) 

save(all_baskets, file = file.path(pkg_dir, "data", "all_baskets.rda"))
