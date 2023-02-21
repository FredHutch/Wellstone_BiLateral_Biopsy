# This script uses the longitudinal study data and 
# takes a close look at complement activation and humoral immune 
# response mediated by circulating immunoglobulin 
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
#library(hg38.HomoSapiens.Gencode.v24)

# parameters
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
fig_dir <- file.path(pkg_dir, "figures", "grant-fig")

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
load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "data", "cluster_df.rda"))

all(cluster_df$sample_name == colnames(sanitized.dds)) # checked
sanitized.dds$cluster <- as.character(cluster_df$new_cluster_name)
sanitized.dds$cluster[sanitized.dds$pheno_type == "Control"] <- "Control"
sanitized.dds$cluster <- factor(sanitized.dds$cluster,
                                levels=c("Control", "Mild", "Moderate", "IG-High", "High", "Muscle-Low"))

muscle_low <- cluster_df %>% dplyr::filter(new_cluster_name == "Muscle-Low") %>% 
  pull(sample_name) %>% as.character(.)

sub_dds <- sanitized.dds[, !sanitized.dds$cluster=="Muscle-Low"]

annotation <- as.data.frame(rowData(sanitized.dds)) %>%
  dplyr::select(gene_id, gene_name)

#
# extract genes from enriched GO terms of interest: total 21 gene
#  
complement_go <- 
  data.frame(id = c("GO:0002455", "GO:0006958", "GO:0006956", "GO:0006957"),
              name = c("humoral immune response mediated by circulating immunoglobulin",
                       "complement activation classical pathway",
                       "complement activation", # this woulde include some circulating immunoglobulin
                       "complement activation alternative pathway"))

gene2cat <- goseq::getgo(high_sig$ens_id, "hg38", "ensGene", fetch.cats = "GO:BP")
names(gene2cat) <- high_sig$ens_id
cat2gene <- goseq:::reversemapping(gene2cat)[complement_go$id]
complement_df <- map_dfr(cat2gene, function(cat_gene) {
  high_sig %>% dplyr::filter(ens_id %in% cat_gene)
}) %>%
  dplyr::distinct(ens_id, .keep_all=TRUE)
save(complement_df, file=file.path(pkg_dir, "data", "complement_df.rda"))

#
# boxplot of log10(TPM+1) of complement_df by cluster
#
col_data <- as.data.frame(colData(sub_dds)) %>%
  dplyr::select(cluster, pheno_type, visit) %>%
  rownames_to_column(var="sample_id")

data <- t(assays(sub_dds)[["TPM"]][complement_df$gencode_id, ]) %>% 
  as.data.frame() %>% rownames_to_column(var="sample_id") %>%
  gather(key=gencode_id, value=TPM, -sample_id) %>%
  dplyr::mutate(log10TPM = log10(TPM+1)) %>%
  dplyr::left_join(complement_df, by="gencode_id") %>%
  dplyr::left_join(col_data, by="sample_id") %>%
  dplyr::mutate(subject = str_replace(sample_id, "b$|b1$|-1$", "")) 

ggplot(data, aes(x=cluster, y=log10TPM)) +
  geom_boxplot(outlier.shape = NA, width=0.7) +
  geom_jitter(width=0.2, size=0.7) + 
  labs(y=bquote(log[10]~"(TPM+1)"), title="Complement") +
  theme_bw() +
  facet_wrap(~gene_name, nrow=3, scales="free_y") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=0.5, size=10),
        axis.title.y = element_text(size=9),
        legend.justification=c(0,1), legend.position=c(0, 1)) 
ggsave(file=file.path(pkg_dir, "figures", "complement-log10TPM-boxplot-by-cluster.pdf"))


ggplot(data, aes(x=cluster, y=TPM)) +
  geom_boxplot(outlier.shape = NA, width=0.7) +
  geom_jitter(width=0.2, size=0.7) + 
  labs(y="TPM", title="Complement") +
  theme_bw() +
  facet_wrap(~gene_name, nrow=3, scales="free_y") +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=0.5, size=10),
        axis.title.y = element_text(size=9),
        legend.justification=c(0,1), legend.position=c(0, 1)) 
ggsave(file=file.path(pkg_dir, "figures", "complement-TPM-boxplot-by-cluster.pdf"))

#
# Correlation between two visit?
#


#
# correlation between two visits
#

# remove control; add Visit column
# get pairs

data_corr <- data %>% dplyr::filter(pheno_type == "FSHD") %>%
  dplyr::select(gene_name, visit, log10TPM, subject) %>%
  spread(key=visit, value=log10TPM) %>%
  drop_na("I", "II")

cor_visits <- data_corr %>%
  group_by(gene_name) %>%
  summarise(cor = cor(I, II)) %>%
  dplyr::mutate(cor = paste0("cor=", format(cor, digits=2)))

ggplot(data_corr, aes(x=`I`, y=`II`)) +
  geom_point(size=0.7) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ gene_name, nrow=4, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_visits, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x=bquote("Visit I:"~log[10]~"(TPM+1)"), 
       y=bquote("Visit II:" ~ log[10]~"(TPM+1)"), title="Complement activation") 
ggsave(file=file.path(pkg_dir, "figures", "complement-corr-log10TPM.pdf"))

  

#
# senstive in Mild/Moderate FSHD ?
#

m_dds <- sanitized.dds[, sanitized.dds$cluster %in% c("Control", "Mild", "Moderate")]
DESeq2::design(m_dds) <- ~ pheno_type
m_dds <- estimateDispersions(m_dds)
m_dds <- DESeq(m_dds)
res <- as.data.frame(results(m_dds, alpha=0.05, lfcThreshold=1)) %>%
  rownames_to_column(var="gene_id") %>%
  dplyr::select(gene_id, log2FoldChange, padj) %>%
  dplyr::mutate(up_sig = padj < 0.05 & log2FoldChange > 0) 

tpm <- assays(sanitized.dds)[["TPM"]]
com_avg_TPM <- 
  data.frame(avg_TPM_controls = 
                        rowMeans(tpm[, sanitized.dds$cluster == "Control"]),
             avg_TPM_M = 
                        rowMeans(tpm[, sanitized.dds$cluster %in% c("Mild", "Moderate")]),
             avg_TPM_H = 
                        rowMeans(tpm[, sanitized.dds$cluster %in% c("IG-High", "High")])) %>%
  rownames_to_column(var="gene_id")

com_candidates <- res %>% dplyr::filter(gene_id %in% complement_df$gencode_id) %>%
  dplyr::left_join(annotation, by="gene_id") %>% 
  left_join(com_avg_TPM, by="gene_id") %>%
  relocate(gene_name, .after=gene_id)

com_candidates  

# are the comlement candidates DE in Modereate?
moderate_sig <- readxl::read_xlsx(file.path(pkg_dir, "extdata", 
    "suppl_table_5_Candidate_Biomarkers_and_Enriched_Go.xlsx"),
    sheet=5, skip=3) %>%
  dplyr::select(gene_name, gencode_id) %>%
  dplyr::mutate(ens_id = str_replace(gencode_id, "\\..*", ""))
moderate_sig %>% dplyr::filter(gencode_id %in% complement_df$gencode_id)  

m2_dds <- sanitized.dds[, sanitized.dds$cluster %in% c("Control", "Moderate")]
DESeq2::design(m2_dds) <- ~ pheno_type
m2_dds <- estimateDispersions(m2_dds)
m2_dds <- DESeq(m2_dds)
as.data.frame(results(m_dds, alpha=0.05, lfcThreshold=1)) %>%
  rownames_to_column(var="gene_id") %>%
  dplyr::select(gene_id, log2FoldChange, padj) %>%
  dplyr::mutate(up_sig = padj < 0.05 & log2FoldChange > 0) %>%
  dplyr::filter(gene_id %in% complement_df$gencode_id)  


#
# complement in BiLat: between bilateral TAs
#
load(file.path(pkg_dir, "data", "dds.rda"))
gencode_v35 <-  as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id_v35") %>%
  dplyr::select(gene_id_v35) %>%
  dplyr::mutate(ens_id = str_replace(gene_id_v35, "\\..*", ""))

new_complement_df <- complement_df %>%
  left_join(gencode_v35, by="ens_id")
bilat_TPM <- assays(dds)[["TPM"]][new_complement$gene_id_v35, ]

# join the longitudinal
bl_avg_TPM <- data.frame(avg_TPM_BL = rowMeans(bilat_TPM)) %>%
  rownames_to_column(var="gene_id_v35") %>%
  dplyr::mutate(ens_id =  str_replace(gene_id_v35, "\\..*", ""))

com_candidates %>% 
  dplyr::mutate(ens_id = str_replace(gene_id, "\\..*", ""))  %>%
  left_join(bl_avg_TPM, by="ens_id") %>%
  dplyr::select(-gene_id_v35, -ens_id, -gene_id)
  
# per-gene correlation between bilat TAs
col_data <- as.data.frame(colData(dds)) %>% dplyr::select(sample_name, Subject, location)
data  <- t(bilat_TPM) %>% as.data.frame() %>%
  rownames_to_column(var="sample_name") %>%
  left_join(col_data, by="sample_name") %>%
  tidyr::gather(key=gene_id, value=TPM, -sample_name, -Subject, -location) %>%
  dplyr::select(-sample_name) %>%
  tidyr::spread(key=location, value=TPM) %>%
  tidyr::drop_na(L, R) %>%
  rename(gene_id_v35=gene_id) %>%
  left_join(new_complement_df, by="gene_id_v35") %>%
  dplyr::mutate(log10L = log10(L+1), log10R=log10(R+1)) %>%
  dplyr::filter(!Subject %in% c("13-0009", "13-0007"))

cor_bilat <- data %>%
  group_by(gene_name) %>%
  summarise(cor = cor(log10L, log10R)) %>%
  dplyr::mutate(cor = paste0("cor=", format(cor, digits=2)))
# just plot IG genes

ig <- c("IGKC", "IGHG4", "IGHG2", "IGHG1", "IGHG3", "IGHV3-23")
data %>% 
  dplyr::filter(gene_name %in% ig) %>%
  ggplot(aes(x=log10L, y=log10R)) +
    geom_point() +
    geom_text_repel(aes(label=Subject), size=1.8, show.legend=FALSE) +
    facet_wrap(~gene_name, nrow=2, scales="free_y") +
    theme_bw() +
    geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
    geom_text(data = dplyr::filter(cor_bilat, gene_name %in% ig), aes(label = cor), 
              x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
    labs(x=bquote("L:"~log[10]~"(TPM+1)"), 
         y=bquote("R:" ~ log[10]~"(TPM+1)"), title="IG") 
    
ggsave(file.path(pkg_dir, "figures", "complement-bilat-corr-ig-onlye.pdf"))

data %>% 
  ggplot(aes(x=log10L, y=log10R)) +
    geom_point(size=0.7) +
    facet_wrap(~gene_name, nrow=4, scales="free") +
    geom_smooth(method="lm", se=FALSE) +
    theme_bw() +
    geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
    geom_text(data = cor_bilat, aes(label = cor), 
              x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
    labs(x=bquote("L:"~log[10]~"(TPM+1)"), 
         y=bquote("R:" ~ log[10]~"(TPM+1)"), title="complement") 
    
ggsave(file.path(pkg_dir, "figures", "complement-bilat-corr.pdf"))

data %>% 
  dplyr::filter(gene_name %in% c("C3", "VSIG4", "C4A", "C6")) %>%
  ggplot(aes(x=log10L, y=log10R)) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_text_repel(aes(label=Subject), size=1.8, show.legend=FALSE) +
    facet_wrap(~gene_name, nrow=2, scales="free") +
    theme_bw() +
    geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
    geom_text(data = dplyr::filter(cor_bilat, gene_name %in% c("C3", "VSIG4", "C4A", "C6")), aes(label = cor), 
              x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
    labs(x=bquote("L:"~log[10]~"(TPM+1)"), 
         y=bquote("R:" ~ log[10]~"(TPM+1)"), title="IG") 
    
ggsave(file.path(pkg_dir, "figures", "test.pdf"))

