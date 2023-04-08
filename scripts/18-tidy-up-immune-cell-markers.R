
library(dplyr)
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
load(file.path(pkg_dir, "data", "bilat_dds.rda"))
load(file.path(pkg_dir, "data", "longitudinal_dds.rda"))

anno_gencode35 <- as.data.frame(rowData(bilat_dds)) %>%
  rownames_to_column(var="gencode35_id") %>% # BiLat study using Gencode 36
  dplyr::mutate(ens_id=str_replace(gencode35_id, "\\..*", "")) %>%
  dplyr::distinct(gene_name, .keep_all = TRUE)

anno_ens88 <- as.data.frame(rowData(longitudinal_dds)) %>%
  rownames_to_column(var="ens88_id") %>% # longitudinal study ens v88
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", "")) %>%
  dplyr::distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::select(ens88_id, ens_id, gene_name)


#
# b cell markers from estipona
# 
MHC_II <- tibble(`Gene` = paste0(c("HLA-DM", "HLA-DO", "HLA-DP", "HLA-DQ", "HLA-DR"), 
                                       c("A", "B")),
                     `Marker Type` = "MHCII",
                     `Protein Type` = "Receptor", 
                     `Localization`="Cell Membrane", note="add on by CJ")
# note: Estipona only select HLA-DOA                     
IL12 <- tibble(Gene = c("IL12A", "IL27"),
                   `Marker Type`= "Breg",
                   `Protein Type` = "Cytokine",
                   `Localization` =	"Secreted", note="add on by CJ")

b_cell_markers <- readxl::read_xlsx(path=file.path(pkg_dir, "extdata", 
                                             "B-cell-markers-Estipona-2020.xlsx")) %>%
  dplyr::filter(!`Protein Type` == "Immunoglobulin") %>% # exclude immunoglobulin (we know IgG are expressed)
  dplyr::select(Gene, `Marker Type`, `Protein Type`, `Localization`) %>% 
  tibble::add_column(note="Estipona 2020") %>%
  dplyr::bind_rows(MHC_II) %>%
  dplyr::bind_rows(IL12) %>%
  dplyr::rename(gene_name = Gene) %>%
  left_join(anno_ens88, by="gene_name") %>%
  left_join(dplyr::select(anno_gencode35, gene_name, gencode35_id), 
            by="gene_name") %>%
  drop_na(ens88_id, gencode35_id)

save(b_cell_markers, file=file.path(pkg_dir, "data", "b_cell_markers.rda"))