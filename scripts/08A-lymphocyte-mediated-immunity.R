# This script takes both the bilateral and longitudinal studies and
# takes a closer look at the lymphocyte and B-cell mediated immunities that
# invode the immunse response mediated by circulating immunoglobulin 
# - GO:0002455 (12/41) humoral immune response mediated by circulating immunoglobulin 
# - GO:0006958 (11/30) complement activation classical pathway
# - GO:0006956 (14/58) complement activation
# - GO:0006957 (4/12) alternative pathway
#


# load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggrepel))
library(corrr)
library(GO.db)
require(goseq)
library(org.Hs.eg.db)

pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"

#
# suppl. Table 5
#
high_sig <- readxl::read_xlsx(file.path(pkg_dir, "extdata", 
    "suppl_table_5_Candidate_Biomarkers_and_Enriched_Go.xlsx"),
    sheet=1, skip=3) %>%
  dplyr::select(gene_name, gencode_id) %>%
  dplyr::mutate(ens_id = str_replace(gencode_id, "\\..*", ""))


#
# tidy santized.dds: remove muscle-low FSHD sampels
#
load(file.path(pkg_dir, "data", "sanitized.rlg.rda"))
load(file.path(pkg_dir, "data", "cluster_df.rda"))

all(cluster_df$sample_name == colnames(sanitized.rlg)) # checked
sanitized.rlg$cluster <- as.character(cluster_df$new_cluster_name)
sanitized.rlg$cluster[sanitized.rlg$pheno_type == "Control"] <- "Control"
sanitized.rlg$cluster <- factor(sanitized.rlg$cluster,
                                levels=c("Control", "Mild", "Moderate", "IG-High", "High", "Muscle-Low"))

muscle_low <- cluster_df %>% dplyr::filter(new_cluster_name == "Muscle-Low") %>% 
  pull(sample_name) %>% as.character(.)

sub_rlg <- sanitized.rlg[, !sanitized.dds$cluster=="Muscle-Low"]

annotation_ens88 <- as.data.frame(rowData(sanitized.dds)) %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::rename(ens88_id = gene_id) %>%
  dplyr::mutate(ens_id = str_replace(ens88_id, "\\..*", ""))

#
# (lymphocyte activation subset): lymphocyte migration, activation and proliferation
#
lymphocyte_go <- 
  data.frame(GO_ID = c("GO:0046649", "GO:0046651", "GO:2000401"),
             NAME = c("lymphocyte activation", "lymphocyte proliferation",
                      "regulation of lymphocyte migration"))

gene2cat <- goseq::getgo(high_sig$ens_id, "hg38", "ensGene", fetch.cats = "GO:BP")
names(gene2cat) <- high_sig$ens_id
cat2gene <- goseq:::reversemapping(gene2cat)[lymphocyte_go$GO_ID]
lymphocyte_act <- map_dfr(cat2gene, function(cat_gene) {
  high_sig %>% dplyr::filter(ens_id %in% cat_gene)
}, .id="GO_ID") 
lymphocyte_act

lymphocyte_act %>% distinct(gene_name, .keep_all=TRUE) # 61 genes


#
# T cell activation
#
tcell_go <- 
  data.frame(GO_ID = c("GO:0042110", "GO:0042098", "GO:0042129"),
             NAME = c("T cell activation", "T cell proliferation", 
                      "Regulation of T cell activation"))

gene2cat <- goseq::getgo(high_sig$ens_id, "hg38", "ensGene", fetch.cats = "GO:BP")
names(gene2cat) <- high_sig$ens_id
cat2gene <- goseq:::reversemapping(gene2cat)[tcell_go$GO_ID]
tcell_prolif <- map_dfr(cat2gene, function(cat_gene) {
  high_sig %>% dplyr::filter(ens_id %in% cat_gene)
}, .id="GO_ID") 

tcell_prolif %>% dplyr::filter(GO_ID=="GO:0042098") %>% as.data.frame()

tcell_prolif %>% distinct(gene_name, .keep_all=TRUE)  %>% as.data.frame()

lymphocyte_act <- lymphocyte_act %>% bind_rows(tcell_prolif) %>% 
  distinct(gene_name, .keep_all=TRUE)  

#
# lymphocyte and B-cell mediated immune response genes
# 
bcell_go <- 
  data.frame(GO_ID = c("GO:0002449", "GO:0019724", "GO:0002455"),
             NAME = c("lympocyte mediated immunity", 
                       "B-cell mediated immuse response",
                       "humoral immune response mediated by circulating immunoglobulin"))
                      
# get sig (high vs. control) genes in the categores
gene2cat <- goseq::getgo(high_sig$ens_id, "hg38", "ensGene", fetch.cats = "GO:BP")
names(gene2cat) <- high_sig$ens_id
cat2gene <- goseq:::reversemapping(gene2cat)[bcell_go$GO_ID]
bcell_df <- map_dfr(cat2gene, function(cat_gene) {
  high_sig %>% dplyr::filter(ens_id %in% cat_gene)
}, .id="GO_ID") 
bcell_df %>% distinct(gene_name, .keep_all=TRUE) # 36

write_xlsx(x=list(`lymphocyte-and-T-call-activation` = lymphocyte_act,
                  `B-cell-mediated-immunity` = bcell_df %>% distinct(gene_name, .keep_all=TRUE)), 
           path=file.path(pkg_dir, "stats", 
           "differentailly-expressed-High-vs-controls-lymphyocyte-complement.xlsx" ))

#
# Estipona 2020 B cell markers
#

# note: the table has no HLA-DPA and HLA-DQB
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

b_markers <- readxl::read_xlsx(path=file.path(pkg_dir, "extdata", 
                                             "B-cell-markers-Estipona-2020.xlsx")) %>%
  dplyr::filter(!`Protein Type` == "Immunoglobulin") %>% # exclude immunoglobulin (we know IgG are expressed)
  dplyr::select(Gene, `Marker Type`, `Protein Type`, `Localization`) %>% 
  tibble::add_column(note="Estipona 2020") %>%
  dplyr::bind_rows(MHC_II) %>%
  dplyr::bind_rows(IL12) %>%
  dplyr::rename(gene_name = Gene) %>%
  left_join(annotation_ens88, by="gene_name") %>%
  drop_na(ens88_id)

# HLA-DR*n-> replaced by HLA-DRA; FCRL4 unknown ID
# IL35*: IL-12α + IL-27β(IL12A; IL27RA, IL27); replaced by IL12A and IL27


#
# DE in high samples
#
de_bcell <- b_markers %>% dplyr::filter(ens88_id %in% high_sig$ens88_id)

b_markers %>% inner_join(high_sig, by="ens88_id")
#
# viz de b cell markers (5)
#
library(pheatmap)
annotation_col <- colData(sub_rlg) %>% as.data.frame() %>%
  dplyr::select(pheno_type, cluster) %>%
  rownames_to_column(var="sample_name")
data <- assay(sub_rlg[de_bcell$ens88_id]) %>% as.data.frame() %>%
  rownames_to_column(var="ens88_id") %>%
  gather(key=sample_name, value=rlog, -ens88_id) %>%
  left_join(de_bcell, by="ens88_id") %>%
  left_join(annotation_col, by= "sample_name")

ggplot(data, aes(x=cluster, y=rlog)) +
  geom_boxplot(width=0.3, outlier.shape=NA) +
  geom_jitter(width=0.1) +
  facet_wrap(~gene_name, nrow=5, scale="free_y") +
  theme_classic() +
  labs(title="B cell markers DE in High group", y="expression (rlog)", x="") +
  theme(plot.title=element_text(hjust = 0.5, size=10))
ggsave(file.path(pkg_dir, "figures", "B-cell-markers-DE.pdf"))  



#
# viz DE MS4A family by boxplot
#

#
# viz LILRB family
#