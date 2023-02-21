#'
#' This script performs bilateral correlation analysis and makes
#' Figure 6 for the manuscript
#'

# Bilat correlation for basket genes
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggrepel))

#
# load dataset
#
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"

load(file.path(pkg_dir, "data", "dds.rda"))
load(file.path(pkg_dir, "data", "all_baskets.rda"))
load(file.path(pkg_dir, "data", "comprehensive_df.rda"))
col_data <- as.data.frame(colData(dds)) %>%
  dplyr::mutate(sample_id = str_replace(sample_name, "[b]*_.*", "")) 

#
# tidy-up some data frames
#
anno_gencode35 <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id") %>% # BiLat study using Gencode 36
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", ""))

#
# basket gene-by-gene bilateral correlation using TPM 
#
bilat_tpm <- map_dfr(names(all_baskets), function(basket_name) {
    id <- all_baskets[[basket_name]]$gencode_v35
    tpm <- assays(dds[id])[["TPM"]] %>% t(.) %>% 
      as.data.frame() %>%
      rownames_to_column(var="sample_id") %>%
      dplyr::mutate(sample_id = str_replace(sample_id, "[b]*_.*", "")) %>%
      gather(key=gene_id, value=TPM, -sample_id) %>% 
      left_join(anno_gencode35, by="gene_id") %>%
      dplyr::select(-ens_id, -gene_id, -gene_type) %>%
      left_join(dplyr::select(col_data, sample_id, location, Subject), by="sample_id") %>%
      dplyr::select(-sample_id) %>%
      spread(key=location, value=TPM) %>% 
      dplyr::filter(!is.na(L), !is.na(R)) %>%
      dplyr::mutate(basket=basket_name)
}) %>% dplyr::mutate(basket=factor(basket, levels=names(all_baskets))) 
#  dplyr::filter(!Subject %in% c("13-0007", "13-0009"))
  

gene_cor <- bilat_tpm %>% 
  dplyr::mutate(L_log = log10(L+1), R_log=log10(R+1)) %>%
  group_by(basket, gene_name) %>% 
  summarise(correlation=cor(L, R), cor_log = cor(L_log, R_log)) 

gene_cor %>% dplyr::filter(!is.na(correlation), !basket %in% c("DUX4-M", "DUX4-M12")) %>%
  ggplot(aes(x=basket, y=correlation)) +
    geom_boxplot(width=0.5, outlier.shape=NA) +
    geom_point(alpha=0.6, color="grey50") +
    theme_minimal() +
    labs(x="", y="Pearson by TPM") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file=file.path(pkg_dir, "figures", "bilat-basket-comparison-cor-by-gene.pdf"), width=2, height=3)  

#
# scatter plot per 
#

# IG
ig_cor <- gene_cor %>% dplyr::filter(basket=="IG")
bilat_tpm %>% dplyr::filter(basket=="IG") %>%
  ggplot(aes(x=L, y=R)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~gene_name, scales="free") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = ig_cor, aes(label=format(correlation, digit=2)), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
  theme_bw() +
  labs(x="L (TPM)", y="R (TPM)", title="IG basket")
ggsave(file=file.path(pkg_dir, "figures", "bilat-basket-comparison-IG.pdf"), width=6, height=4)


# Complement
comp_cor <- gene_cor %>% dplyr::filter(basket=="Complement")
bilat_tpm %>% dplyr::filter(basket=="Complement") %>%
  ggplot(aes(x=L, y=R)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~gene_name, scales="free") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = comp_cor, aes(label=format(correlation, digit=2)), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
  theme_bw() +
  labs(x="L (TPM)", y="R (TPM)", title="Complement basket")
ggsave(file=file.path(pkg_dir, "figures", "bilat-basket-comparison-Complement.pdf"), width=6, height=4)

# ECM
ECM_cor <- gene_cor %>% dplyr::filter(basket=="ECM")
bilat_tpm %>% dplyr::filter(basket=="ECM") %>%
  ggplot(aes(x=L, y=R)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~gene_name, scales="free") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = ECM_cor, aes(label=format(correlation, digit=2)), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
  theme_bw() +
  labs(x="L (TPM)", y="R (TPM)", title="ECM basket")
ggsave(file=file.path(pkg_dir, "figures", "bilat-basket-comparison-ECM.pdf"), width=6, height=4)

# Inflamm
inflamm_cor <- gene_cor %>% dplyr::filter(basket=="Inflamm")
bilat_tpm %>% dplyr::filter(basket=="Inflamm") %>%
  ggplot(aes(x=L, y=R)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~gene_name, scales="free") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = inflamm_cor, aes(label=format(correlation, digit=2)), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
  theme_bw() +
  labs(x="L (TPM)", y="R (TPM)", title="Inflamm basket")
ggsave(file=file.path(pkg_dir, "figures", "bilat-basket-comparison-inflamm.pdf"), width=6, height=4)

# DUX4-M6
dux4m6_cor <- gene_cor %>% dplyr::filter(basket=="DUX4-M6")
bilat_tpm %>% dplyr::filter(basket=="DUX4-M6") %>%
  ggplot(aes(x=L, y=R)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~gene_name, scales="free") +
  #scale_y_continuous(trans='log10') +
  #scale_x_continuous(trans='log10') +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = dux4m6_cor, aes(label=format(correlation, digit=2)), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
  theme_bw() +
  labs(x="L (TPM)", y="R (TPM)", title="DUX4-M6 basket")
ggsave(file=file.path(pkg_dir, "figures", "bilat-basket-comparison-DUX4-M6.pdf"), width=6, height=4)

#
#
#

bk_scores <- bilat_tpm %>% dplyr::mutate(L_log = log10(L+1), R_log=log10(R+1)) %>%
  group_by(basket, Subject) %>% 
  summarise(L_logscore=sum(L_log), R_logscore= sum(R_log), L_score=mean(L), R_score=mean(R)) %>%
  dplyr::filter(! basket %in% c("DUX4-M", "DUX4-M12")) %>%
  dplyr::mutate(basket = factor(basket, levels=c("DUX4-M6", "Inflamm", "ECM", "Complement", "IG"))) 

bk_cor <- bk_scores %>%
  group_by(basket) %>% summarise(cor_log=cor(L_logscore, R_logscore), 
                                 cor_tpm = cor(L_score, R_score))

bk_scores %>% 
ggplot(aes(x=L_logscore, y=R_logscore)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~basket, scale="free") +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = bk_cor, aes(label=format(cor_log, digit=2)), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
  theme_bw() +
    scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  labs(x="L (sum log10TPM)", y="R (sum log10TPM)", title="baskets bilateral correlation using sum log10")
ggsave(file=file.path(pkg_dir, "figures", "bilat-basket-comparison-log10TPM-score.pdf"), width=6, height=4)

bk_scores %>% 
ggplot(aes(x=L_score, y=R_score)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~basket, scale="free") +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = bk_cor, aes(label=format(cor_tpm, digit=2)), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
  theme_bw() +
    scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  labs(x="L (TPM)", y="R (TPM)", title="baskets bilateral correlation using average TPM")
ggsave(file=file.path(pkg_dir, "figures", "bilat-basket-comparison-average-TPM-score.pdf"), width=6, height=4)