# this script generates figures for the grand
#


# load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggrepel))
# define parameters and load data sets
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
fig_dir <- file.path(pkg_dir, "figures", "grant-fig")
load(file.path(pkg_dir, "data", "dds.rda"))
load(file.path(pkg_dir, "data", "mri.rda"))
load(file.path(pkg_dir, "data", "DUX4_positive.rda"))
load(file.path(pkg_dir, "data", "FSHD_signatures_baskets.rda"))


#
# tidy annotation
#
anno_gencode36 <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id") %>% # BiLat study using Gencode 36
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", ""))

#
# tidy DUX4_positive
#
DUX4_positive <- DUX4_positive %>%
  dplyr::filter(! sample_id %in% c("13-0009R", "13-0007R")) %>%
  dplyr::mutate(`DUX4+ (M6)` = if_else(`DUX4+ (M6)`, "DUX4+", "DUX4-")) %>%
  dplyr::mutate(`DUX4+ (M6)` = factor(`DUX4+ (M6)`, levels=c("DUX4+", "DUX4-")))

#
# tidy muscle strenght
#
muscle_strength <- get(load(file.path(pkg_dir, "data", "muscle_strength.rda"))) %>%
  dplyr::rename(Subject=`Record ID`) %>%
  gather(key=`location`, value=`Foot Dorsiflexors`, -Subject) %>%
  dplyr::mutate(sample_id = paste0(Subject, 
                                   str_sub(location, start=1L, end=1L)))

#
# filter FSHD signature baskets
# the goal is to select the biomarkers that are reliably measureable
# 1. SPP1 and THBS1 in Inflamm basket are overlapping with ECM
# 2. after correlation comparison, we decided to remove DMBT1 and TREM2 from inflamm
#   and SPP1 from ECM baskets; so each ECM and inflamm basket obtain six biomarkers 
load(file.path(pkg_dir, "data", "FSHD_signatures_baskets.rda"))
FSHD_signatures_baskets$inflamm <- FSHD_signatures_baskets$inflamm %>% 
  add_row(FSHD_signatures_baskets$stress %>% dplyr::filter(gene_name == "CDKN1A")) %>%
  dplyr::filter(! gene_name %in% c("SPP1", "THBS1", "DMBT1", "TREM2"))

FSHD_signatures_baskets$ecm <- FSHD_signatures_baskets$ecm %>% 
  dplyr::filter(!gene_name == "SPP1")

#
# tidy up the ECM and inflamm baskets
# 
baskets_TPM <- map_dfc(FSHD_signatures_baskets[1:2], function(basket) {
  tmp <- basket %>% dplyr::select(gene_id, gene_name) %>%
    dplyr::mutate(ens_id = str_replace(gene_id, "\\..*", "")) %>%
    dplyr::left_join(anno_gencode36 %>% dplyr::select(gene_id, gene_name, ens_id),
                     by="ens_id", suffix=c(".ens88", ".gencode36")) 
  assays(dds[tmp$gene_id.gencode36])[["TPM"]] %>% as.data.frame() %>%
      summarise(across(where(is.numeric), mean)) %>% t(.) %>% as.data.frame() %>%
      rownames_to_column(var="sample_name") %>%
      dplyr::rename(avgTPM = V1)
}) %>% dplyr::rename(sample_name = `sample_name...1`,
              ECM = `avgTPM...2`,
              Inflamm = `avgTPM...4`) %>%
  dplyr::select(sample_name, ECM, Inflamm) %>%
  dplyr::mutate(sample_id = str_replace(sample_name, "_.*", "")) %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "b$", ""))

#
# merge all variables
#
merge_df <- mri %>% left_join(DUX4_positive, by = "sample_id") %>%
  left_join(baskets_TPM, by = "sample_id") %>%
  dplyr::filter(!is.na(`basket-M6`)) %>%
  dplyr::rename(`DUX4-M6` = `basket-M6`) %>%
  left_join(muscle_strength, by = "sample_id")

#
# Historical controls - basket avg TPM
#
load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "data", "candidates.rda"))
controls_dds <- sanitized.dds[, sanitized.dds$pheno_type == "Control"]
baskets_id <- list(`DUX4-M6` = candidates %>% dplyr::filter(`basket-M6`) %>% pull(gene_id),
                   ECM  = FSHD_signatures_baskets$ecm$gene_id,
                   Inflamm = FSHD_signatures_baskets$inflamm$gene_id)

basket_cntr <- map_dfr(baskets_id, function(gene_id) {
    c(mean=colMeans(assays(controls_dds[gene_id])[["TPM"]]) %>% mean(.),
      sd = colSds(assays(controls_dds[gene_id])[["TPM"]]) %>% sd(.))
}) %>% 
  add_column(basket=names(baskets_id)) %>%
  dplyr::mutate(mean_minus_sd = mean - sd,
                mean_plus_sd = mean + sd) %>%
  dplyr::mutate(mean_minus_sd = if_else(mean_minus_sd < 0, 0, mean_minus_sd))                

control_basket <- map_dfr(names(baskets_id), function(name) {
  assays(controls_dds[baskets_id[[name]]])[["TPM"]] %>% as.data.frame() %>%
      summarise(across(where(is.numeric), mean)) %>% t(.) %>% as.data.frame() %>%
      dplyr::rename(TPM = V1) %>%
      rownames_to_column(var="sample_id") %>%
      add_column(STIR_status = "Control", .after="sample_id") %>%
      add_column(basket=name, .after="STIR_status")
}) 
 
#
# color
#
library(wesanderson)
pal <- wes_palette("Darjeeling1", n=5)
stir_pal <- c("STIR-" = pal[5], "STIR+"= pal[1])
cntr_pal <- pal[2]

#####################################################
#
# (1) DUX4 M6 basket vs Fat_Infilt_Percent (color by STIR)
#
#####################################################
ggplot(merge_df, aes(x=Fat_Infilt_Percent, y=`DUX4-M6`)) +
  geom_point(aes(color=STIR_status)) +
  theme_minimal() +
  geom_hline(yintercept=basket_cntr %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean),
             color="gray10") +
  #geom_hline(yintercept=basket_cntr %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean_minus_sd),
  #           color="gray10", linetype="dashed") +
  geom_hline(yintercept=basket_cntr %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean_plus_sd),
             color="gray10", linetype="dashed") +
  scale_color_manual(values = stir_pal) +           
  scale_y_continuous(trans='log10') +  
  labs(x="Fat infiltration percent", y="DUX4 M6 (TPM)") 
ggsave(file=file.path(fig_dir, "DUX4-score-vs-fat-infilt-percent.pdf"),
       width=6, height=4)

# also similar plot with M12, M32, and inflam, ig and ecm baskets


#####################################################
#
# (2) DUX4 M6 basket vs muscle strength (color by STIR)
#
#####################################################
merge_df %>% dplyr::filter(!is.na(`Foot Dorsiflexors`)) %>%
 ggplot(aes(x=`Foot Dorsiflexors`, y=`DUX4-M6`)) +
   geom_point(aes(color=STIR_status)) +
   theme_minimal() +
   geom_smooth(method="lm", se=FALSE) +
   geom_hline(yintercept=basket_cntr %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean),
              color="gray10") +
   geom_hline(yintercept=basket_cntr %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean_plus_sd),
              color="gray10", linetype="dashed") +
   scale_color_manual(values = stir_pal) +           
   scale_y_continuous(trans='log10') +  
   labs(y="DUX4 M6 (TPM)") 
ggsave(file=file.path(fig_dir, "DUX4-score-vs-muscle-strength.pdf"),
       width=6, height=4)   
       

#####################################################
#
# (3) STIR vs DUX4-M6, ECM, and inflamm baskets
#     color by STIR status; add *** referring p-value of F-statistics
#     note: all pvalue < 2e-4
#     add control boxplots
#####################################################
#
#summary(lm(Inflamm ~ STIR_status, data=merge_df)) #2e-4
#summary(lm(ECM ~ STIR_status, data=merge_df)) #6e-6

tmp_tidy <- merge_df %>% 
  dplyr::select(sample_id, ECM, Inflamm, `DUX4-M6`, STIR_status) %>%
  dplyr::filter(!is.na(STIR_status)) %>%
  gather(key=basket, value=TPM, -sample_id, -STIR_status,) %>%
  add_row(control_basket) %>%
  dplyr::mutate(STIR_status = factor(STIR_status, levels=c("Control", "STIR-", "STIR+")),
                basket = factor(basket, levels=c("DUX4-M6", "Inflamm", "ECM")))

  ggplot(tmp_tidy, aes(x=STIR_status, y=TPM)) +
    geom_boxplot(width=0.7, outlier.shape=NA, fill="grey75", alpha=0.5) +
    geom_jitter(width = 0.3, size=0.5, alpha=0.5) +
    facet_wrap(~ basket, scales="free_y", nrow=1) +
    theme_bw() +
    labs(x="", y="basket scores (TPM)") +
    scale_y_continuous(trans='log10') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggsave(file=file.path(fig_dir, "STIR-vs-all-baskets.pdf"), width=4.5, height=3 )


#####################################################
#
# (4) basekt correlation by ggpair
#
#####################################################
library(GGally)
merge_df %>% dplyr::select(`DUX4-M6`, Inflamm, ECM) %>%
  ggpairs() + 
     scale_y_continuous(trans='log10') +
     scale_x_continuous(trans='log10') +
     theme_bw() 
ggsave(file=file.path(fig_dir, "basekt-m6-other-baskets.pdf"), width=6, height=6)  
