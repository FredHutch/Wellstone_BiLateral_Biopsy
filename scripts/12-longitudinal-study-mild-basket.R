# longitudinal study: mild vs. control markers

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(corrr)
library(GO.db)
library(readxl)
library(writexl)
library(org.Hs.eg.db)

pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
fig_dir <- file.path(pkg_dir, "figures", "longitudinal_cohort")

load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "data", "cluster_df.rda"))
load(file.path(pkg_dir, "data", "all_baskets.rda"))
load(file.path(pkg_dir, "data", "candidates.rda"))
all(cluster_df$sample_name == colnames(sanitized.dds)) # checked
sanitized.dds$cluster <- as.character(cluster_df$new_cluster_name)
sanitized.dds$cluster[sanitized.dds$pheno_type == "Control"] <- "Control"
sanitized.dds$cluster <- factor(sanitized.dds$cluster,
                                levels=c("Control", "Mild", "Moderate", "IG-High", "High", "Muscle-Low"))

anno_ens88 <- as.data.frame(rowData(sanitized.dds)) %>%
  dplyr::select(gene_id, gene_name)

col_data <- as.data.frame(colData(sanitized.dds)) %>%
  dplyr::select(sample_name, pheno_type, cluster, fat_fraction, STIR_rating, visit) %>%
  dplyr::mutate(STIR_status = if_else(STIR_rating > 0, "STIR+", "STIR-")) %>%
  dplyr::mutate(STIR_status = factor(STIR_status, levels=c("STIR-", "STIR+")))

# color stir
library(wesanderson)
pal <- wes_palette("Darjeeling1", n=5)
stir_pal <- c("STIR-" = pal[5], "STIR+"= pal[1])
cntr_pal <- pal[2]

#
# Early markers!!!  see https://fredhutch.github.io/RWellstone_FSHD_muscle_biopsy/mild-fshd.html
# objective: two baskets: one choosed objectively whether or not already exists in other baskets
#            one choosed to be new candidate that is not in any other basket; we want to see if any
#            of the genes of baskets more sensitive to 
#

# here are the candidates orederd by random forest from https://fredhutch.github.io/RWellstone_FSHD_muscle_biopsy/mild-fshd.html

# top_rf are a subset of 164 candidates and top 12 ordered by raondom forest
top_rf <- c("FGF18", "CAPG", "CDKN1A", "FOSB", 
            "GDF15", "GLIPR2", "LINC00152", 
            "RUNX1", "COMP", "HDC", "PRAMEF4", 
            "ACTA1", "C3AR1", "TIMP1", "CLCLC1")

mild_markers <- readxl::read_xlsx(path=file.path(pkg_dir,"extdata", "suppl_table_7_mild_FSHD_164_potential_markers.xlsx"),
                                  skip=4, sheet=1) %>%
  dplyr::select(-contains("#")) %>%                                  
  dplyr::filter(AUC >= 0.9, gene_name %in% top_rf) %>%
  left_join(anno_ens88, by="gene_name")                     

# let's not spend too much time on this

#
# Viz: TPM by cluster 
#
data <- assays(sanitized.dds[mild_markers$gene_id])[["TPM"]] %>%
  as.data.frame() %>%
  rownames_to_column(var="gene_id") %>%
  left_join(anno_ens88, by="gene_id") %>%
  gather(key=sample_name, value=TPM, -gene_name, -gene_id) %>%
  dplyr::mutate(sample_name = if_else(sample_name=="32-0002b1", "32-0002b", sample_name)) %>%
  left_join(col_data, by="sample_name") 

ggplot(data, aes(x=cluster, y=TPM)) +
  geom_boxplot(width=0.3, outlier.shape=NA) +
  geom_jitter(width=0.3, size=1, alpha=0.5) +
  facet_wrap(~gene_name, scale="free_y") +
  scale_y_continuous(trans='log10') +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file=file.path(fig_dir, "mild-markers-by-cluster.pdf")  )

#
# Viz: TPM vs. fat_fraction
#
data %>% dplyr::filter(!is.na(fat_fraction)) %>%
ggplot(aes(x=fat_fraction, y=TPM)) +
  geom_point(size=1, aes(color=STIR_status)) +
  facet_wrap(~gene_name, scale="free_y", nrow=4) +
  scale_y_continuous(trans='log10') +  
  theme_bw() +
  scale_color_manual(values=stir_pal) +
  theme(legend.position="none")
ggsave(file=file.path(fig_dir, "mild_markers-TPM-vs-fat-fraction.pdf"))  

# what if not include muscle-low?
data %>% dplyr::filter(!is.na(fat_fraction), !cluster == "Muscle-Low") %>%
ggplot(aes(x=fat_fraction, y=TPM)) +
  geom_point(size=1, aes(color=STIR_status)) +
  facet_wrap(~gene_name, scale="free_y", nrow=4) +
  scale_y_continuous(trans='log10') +  
  theme_bw() +
  scale_color_manual(values=stir_pal) +
  theme(legend.position="none")
ggsave(file=file.path(fig_dir, "mild-markers-TPM-vs-fat-fraction-exclude-muscle-low.pdf"))  

#
# Viz: TPM vs STIR_status
#
data %>% dplyr::filter(!is.na(STIR_status)) %>%
ggplot(aes(x=STIR_status, y=TPM)) +
  geom_boxplot(width=0.3, aes(color=STIR_status)) +
  facet_wrap(~gene_name, scale="free_y", nrow=4) +
  scale_y_continuous(trans='log10') +  
  theme_bw() +
  scale_color_manual(values=stir_pal) +
  theme(legend.position="none")
ggsave(file=file.path(fig_dir, "mild-markers-TPM-by-STIR-status.pdf"))  

data %>% dplyr::filter(!is.na(STIR_status), !cluster == "Muscle-Low") %>%
ggplot(aes(x=STIR_status, y=TPM)) +
  geom_boxplot(width=0.3, aes(color=STIR_status)) +
  facet_wrap(~gene_name, scale="free_y", nrow=4) +
  scale_y_continuous(trans='log10') +  
  theme_bw() +
  scale_color_manual(values=stir_pal) +
  theme(legend.position="none")
ggsave(file=file.path(fig_dir, "mild-markers-TPM-by-STIR-status-exclude-muscle-low.pdf"))

#
# how about IG vs. fat fraction?
#
ig_data <- assays(sanitized.dds[all_baskets$IG$gene_id_v88])[["TPM"]] %>%
  as.data.frame() %>%
  rownames_to_column(var="gene_id") %>%
  left_join(anno_ens88, by="gene_id") %>%
  gather(key=sample_name, value=TPM, -gene_name, -gene_id) %>%
  dplyr::mutate(sample_name = if_else(sample_name=="32-0002b1", "32-0002b", sample_name)) %>%
  left_join(col_data, by="sample_name") 

ig_data %>% dplyr::filter(!is.na(STIR_status), !is.na(TPM)) %>%
  ggplot(aes(x=fat_fraction, y=TPM)) +
    geom_point(size=1, aes(color=STIR_status)) +
    facet_wrap(~gene_name, scale="free_y", nrow=4) +
    scale_y_continuous(trans='log10') +  
    theme_bw() +
    scale_color_manual(values=stir_pal) +
    theme(legend.position="none")
ggsave(file=file.path(fig_dir, "ig-markers-TPM-by-STIR-status.pdf"))

# boxplot by visits
ig_data %>% dplyr::filter(!is.na(STIR_status), !is.na(TPM)) %>%
ggplot(aes(x=visit, y=TPM)) +
  geom_boxplot(width=0.3, aes(fill=visit)) +
  facet_wrap(~gene_name, scale="free_y", nrow=4) +
  scale_y_continuous(trans='log10') +  
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme(legend.position="none")
ggsave(file=file.path(fig_dir, "ig-markers-TPM-by-visits.pdf"))

#
# how about a better visulation between visit I and II 
#

##################################### below are the parallel analysis for Bilat ##################
# 
# bilat: TPM (per gene) vs fat_infilt
#
load(file.path(pkg_dir, "data", "dds.rda"))
load(file.path(pkg_dir, "data", "comprehensive_df.rda"))

anno_gencode35 <- anno_gencode35 <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gencode_id") 

tmp_markers <- mild_markers %>% dplyr::select(gene_id, gene_name) %>% 
  dplyr::mutate(gene_name=if_else(gene_name=="LINC00152", "CYTOR", gene_name)) %>%
  left_join(anno_gencode36, by="gene_name") 


tpm_bilat <- assays(dds[tmp_markers$gencode_id])[["TPM"]] %>%
  as.data.frame() %>%
  rownames_to_column(var="gencode_id") %>%
  left_join(anno_gencode35, by="gencode_id") %>%
  dplyr::select(-gencode_id, -gene_type) %>%
  gather(key=sample_name, value=TPM, -gene_name) %>%
  dplyr::mutate(sample_id = str_replace(sample_name, "_.*", "") )%>%
  dplyr::mutate(sample_id = str_replace(sample_id, "b$|b1$|-1$", "")) %>%
  dplyr::left_join(comprehensive_df, by="sample_id") %>%
  dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R"))

tpm_bilat %>% dplyr::select(gene_name, TPM, Fat_Infilt_Percent) %>%
  group_by(gene_name) %>%
  corrr::correlate(.[-1])[2,2]
# is there any linear relaship fat [0, 0.2]
tpm_bilat %>% dplyr::filter(!is.na(Fat_Infilt_Percent)) %>% 
  #dplyr::filter(Fat_Infilt_Percent < 0.2) %>%
  ggplot(aes(x=Fat_Infilt_Percent, y=TPM)) +
    geom_point(size=1, alpha=0.5, aes(color=STIR_status)) +
    #geom_smooth(method="lm") +
    facet_wrap(~gene_name, scale="free_y", nrow=4) +
    scale_y_continuous(trans='log10') +  
    theme_bw() +
    scale_color_manual(values=stir_pal) +
    theme(legend.position="none")
ggsave(file=file.path(fig_dir, "bilat-mild-markers-TPM-vs-fat-infilt-percent.pdf"))  

tpm_bilat %>% dplyr::filter(!is.na(Fat_Infilt_Percent)) %>% 
  dplyr::filter(Fat_Infilt_Percent < 0.2) %>%
  ggplot(aes(x=Fat_Infilt_Percent, y=TPM)) +
    geom_point(aes(color=STIR_status)) +
    geom_smooth(method="lm") +
    facet_wrap(~gene_name, scale="free_y", nrow=4) +
    scale_y_continuous(trans='log10') +  
    theme_bw() +
    scale_color_manual(values=stir_pal) +
    theme(legend.position="none")
ggsave(file=file.path(fig_dir, "bilat-mild-markers-TPM-vs-fat-infilt-percent-low.pdf"))  


ig_bilat <- assays(dds[all_baskets$IG$gencode_v35])[["TPM"]] %>%
  as.data.frame() %>%
  rownames_to_column(var="gencode_id") %>%
  left_join(anno_gencode35, by="gencode_id") %>%
  dplyr::select(-gencode_id, -gene_type) %>%
  gather(key=sample_name, value=TPM, -gene_name) %>%
  dplyr::mutate(sample_id = str_replace(sample_name, "_.*", "") )%>%
  dplyr::mutate(sample_id = str_replace(sample_id, "b$|b1$|-1$", "")) %>%
  dplyr::left_join(comprehensive_df, by="sample_id") %>%
  dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R"))

ig_bilat %>% dplyr::filter(!is.na(Fat_Infilt_Percent)) %>% 
  #dplyr::filter(Fat_Infilt_Percent < 0.2) %>%
  ggplot(aes(x=Fat_Infilt_Percent, y=TPM)) +
    geom_point(size=1, alpha=0.5, aes(color=STIR_status)) +
    #geom_smooth(method="lm") +
    facet_wrap(~gene_name, scale="free_y", nrow=4) +
    scale_y_continuous(trans='log10') +  
    theme_bw() +
    scale_color_manual(values=stir_pal) +
    theme(legend.position="none")
ggsave(file=file.path(fig_dir, "bilat-ig-basket-TPM-vs-fat-infilt-percent.pdf"))    

ig_bilat %>% dplyr::filter(!is.na(Fat_Infilt_Percent)) %>% 
  dplyr::filter(Fat_Infilt_Percent < 0.2) %>%
  ggplot(aes(x=Fat_Infilt_Percent, y=TPM)) +
    geom_point(aes(color=STIR_status)) +
    geom_smooth(method="lm") +
    facet_wrap(~gene_name, scale="free_y", nrow=4) +
    scale_y_continuous(trans='log10') +  
    theme_bw() +
    scale_color_manual(values=stir_pal) +
    theme(legend.position="none")
ggsave(file=file.path(fig_dir, "bilat-ig-basket-TPM-vs-fat-infilt-percent-low.pdf"))  