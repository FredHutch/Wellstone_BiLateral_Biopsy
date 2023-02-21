# longitudinal study: the vairability of the ECM and inflamm basket genes between two visits
# 

# load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggrepel))
library("corrr")

# parameters
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
fig_dir <- file.path(pkg_dir, "figures", "grant-fig")


#
# tidy santized.dds: remove muscle-low FSHD sampels
#
load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "data", "cluster_df.rda"))
all(cluster_df$sample_name == colnames(sanitized.dds)) # checked
muscle_low <- cluster_df %>% dplyr::filter(new_cluster_name == "Muscle-Low") %>% 
  pull(sample_name) %>% as.character(.)
design(sanitized.dds) < ~ pheno_type
sanitized.dds <- sanitized.dds[, ! colnames(sanitized.dds) %in% muscle_low]
sanitized.rlg <- rlog(sanitized.dds)

#
# filter FSHD signature baskets
# the goal is to select the biomarkers that are reliably measureable
# 1. SPP1 and THBS1 in Inflamm basket are overlapping with ECM
# 2. after correlation comparison, we decided to remove DMBT1 and TREM2 from inflamm
#   and SPP1 from ECM baskets; so each ECM and inflamm basket obtain six biomarkers 
#
load(file.path(pkg_dir, "data", "FSHD_signatures_baskets.rda"))
FSHD_signatures_baskets$inflamm <- FSHD_signatures_baskets$inflamm %>% 
  add_row(FSHD_signatures_baskets$stress %>% dplyr::filter(gene_name == "CDKN1A")) %>%
  dplyr::filter(! gene_name %in% c("SPP1", "THBS1", "DMBT1", "TREM2"))

FSHD_signatures_baskets$ecm <- FSHD_signatures_baskets$ecm %>% 
  dplyr::filter(!gene_name == "SPP1")
  
#
# tidy datasets
#
anno_ens <- as.data.frame(rowData(sanitized.dds)) %>%
  dplyr::select(gene_id, gene_name) %>% # BiLat study using Gencode 36
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", ""))

baskets_id <- list(`DUX4-M6` = candidates %>% dplyr::filter(`basket-M6`) %>% pull(gene_id),
                   `DUX4-M12` =  candidates %>% dplyr::filter(`basket-M12`) %>% pull(gene_id),
                   ECM  = FSHD_signatures_baskets$ecm$gene_id,
                   Inflamm = FSHD_signatures_baskets$inflamm$gene_id)

#
# Year I & II pairs: sample_id = rownames of the dds instance
#
visit <- colData(sanitized.dds[, sanitized.dds$pheno_type == "FSHD"]) %>% as.data.frame() %>%
  dplyr::select(visit) %>%
  rownames_to_column(var="sample_id") 

# sanity check and get subjects of paired
sanity_check <- visit %>%
  dplyr::mutate(subject = str_replace(sample_id, "b$|b1$|-1$", ""))
# how many pairs: 26 subjects
sum(table(sanity_check$subject) == 2)
# how many subject total: 35 subjects
length(unique(sanity_check$subject)  )
# get paried subjects (26)
paired_subject <- names(which(table(sanity_check$subject) == 2))


#
# control basket rlog
#
controls_rlg <- sanitized.rlg[, sanitized.rlg$pheno_type=="Control"]
# get average rlg for each basket genes
baskets_cntr_rlg <- map(baskets_id, function(id) {
  rowMeans(assay(controls_rlg[id]))
}) 

#
# get log2FoldChagne = sample rlg - average controls rlog
#
baskets_rlgFC <- map(names(baskets_id), function(basket_name){
    id <- baskets_id[[basket_name]]
    rlog <- assay(sanitized.rlg)[id, sanitized.rlg$pheno_type == "FSHD"] 
    rlogFC <- rlog - baskets_cntr_rlg[[basket_name]] 
    t(rlogFC) %>% as.data.frame() %>%
      rownames_to_column(var="sample_id") %>%
      dplyr::mutate(subject = str_replace(sample_id, "b$|b1$|-1$", "")) %>%
      gather(key=gene_id, value=rlogFC, -sample_id, -subject) %>%
      left_join(visit, by="sample_id") %>%
      dplyr::select(-sample_id) %>%
      dplyr::filter(subject %in% paired_subject) %>%
      spread(key=visit, value=rlogFC) %>%
      left_join(anno_ens, by="gene_id") 
}) 
names(baskets_rlgFC) <- names(baskets_id)


#
# get tpm for the basket genes
#
baskets_log10tpm <- map(baskets_id, function(id) {
  tpm <- assays(sanitized.dds)[["TPM"]][id, sanitized.dds$pheno_type == "FSHD" ] 
  tpm <- log10(tpm + 1) %>%
    t(.) %>% as.data.frame() %>%
    rownames_to_column(var="sample_id") %>%
    dplyr::mutate(subject = str_replace(sample_id, "b$|b1$|-1$", "")) %>%
    gather(key=gene_id, value=TPM, -sample_id, -subject) %>%
    left_join(visit, by="sample_id") %>%
    dplyr::select(-sample_id) %>%
    dplyr::filter(subject %in% paired_subject) %>%
    spread(key=visit, value=TPM) %>%
    left_join(anno_ens, by="gene_id") 
}) 

baskets_tpm <- map(baskets_id, function(id) {
  tpm <- assays(sanitized.dds)[["TPM"]][id, sanitized.dds$pheno_type == "FSHD" ] %>%
    t(.) %>% as.data.frame() %>%
    rownames_to_column(var="sample_id") %>%
    dplyr::mutate(subject = str_replace(sample_id, "b$|b1$|-1$", "")) %>%
    gather(key=gene_id, value=TPM, -sample_id, -subject) %>%
    left_join(visit, by="sample_id") %>%
    dplyr::select(-sample_id) %>%
    dplyr::filter(subject %in% paired_subject) %>%
    spread(key=visit, value=TPM) %>%
    left_join(anno_ens, by="gene_id") 
})


# corr between visit I and II: tpm and log10tpm
cor_visits <- map(baskets_log10tpm, function(x) {
  x %>% group_by(gene_name) %>%
  summarise(cor = cor(I, II)) %>%
  dplyr::mutate(cor = paste0("cor=", format(cor, digits=2)))
})
  
cor_visits_TPM <- map(baskets_tpm, function(x) {
  x %>% group_by(gene_name) %>%
  summarise(cor = cor(I, II)) %>%
  dplyr::mutate(cor = paste0("cor=", format(cor, digits=2)))
})

cor_visits_rlogFC<- map(baskets_rlgFC, function(x) {
  x %>% group_by(gene_name) %>%
  summarise(cor = cor(I, II)) %>%
  dplyr::mutate(cor = paste0("cor=", format(cor, digits=2)))
})


#
# ggplots - scatter and cor
#

# ECM
ggplot(baskets_log10tpm$ECM, aes(x=`I`, y=`II`)) +
  geom_point() +
  #geom_text_repel(aes(label=subject), size=1.8, show.legend=FALSE) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ gene_name, nrow=2, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_visits$ECM, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I (log10(TPM+1))", y="Visit II (log10(TPM+1))", title="ECM") 
ggsave(file=file.path(fig_dir, "basket-ECM-two-visits-log10TPM.pdf"))



# DUX4-M6
ggplot(baskets_log10tpm$`DUX4-M6`, aes(x=`I`, y=`II`)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ gene_name, nrow=2, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_visits$`DUX4-M6`, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I (log10(TPM+1))", y="Visit II (log10(TPM+1))", title="DUX4 M6")  
ggsave(file=file.path(fig_dir, "basket-DUX4-M6-two-visits-log10TPM.pdf"))

# DUX4-M12
ggplot(baskets_log10tpm$`DUX4-M12`, aes(x=`I`, y=`II`)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ gene_name, nrow=3, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_visits$`DUX4-M12`, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I (log10(TPM+1))", y="Visit II (log10(TPM+1))", title="DUX4 M12") 
ggsave(file=file.path(fig_dir, "basket-DUX4-M12-two-visits-log10TPM.pdf"))

# Inflamm
ggplot(baskets_log10tpm$Inflamm, aes(x=`I`, y=`II`)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ gene_name, nrow=2, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_visits$Inflamm, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I (log10(TPM+1))", y="Visit II (log10(TPM+1))", title="Inflamm") 
ggsave(file=file.path(fig_dir, "basket-Inflamm-two-visits-log10TPM.pdf"))

#
# scatter plot of baskets_rlogFC
#

# ECM
ggplot(baskets_rlgFC$ECM, aes(x=`I`, y=`II`)) +
  geom_point() +
  #geom_text_repel(aes(label=subject), size=1.8, show.legend=FALSE) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ gene_name, nrow=2, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_visits_rlogFC$ECM, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I (rlogFC)", y="Visit II (rlogFC)", title="ECM") 
ggsave(file=file.path(fig_dir, "basket-ECM-two-visits-rlogFC.pdf"))


# DUX4-M6
ggplot(baskets_rlgFC$`DUX4-M6`, aes(x=`I`, y=`II`)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ gene_name, nrow=2, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_visits_rlogFC$`DUX4-M6`, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I (rlogFC)", y="Visit II (rlogFC)", title="DUX4 M6")  
ggsave(file=file.path(fig_dir, "basket-DUX4-M6-two-visits-rlogFC.pdf"))

# DUX4-M12
ggplot(baskets_rlgFC$`DUX4-M12`, aes(x=`I`, y=`II`)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ gene_name, nrow=3, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_visits_rlogFC$`DUX4-M12`, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I (rlogFC))", y="Visit II (rlogFC)", title="DUX4 M12") 
ggsave(file=file.path(fig_dir, "basket-DUX4-M12-two-visits-rlogFC.pdf"))

# Inflamm
ggplot(baskets_rlgFC$Inflamm, aes(x=`I`, y=`II`)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ gene_name, nrow=2, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_visits_rlogFC$Inflamm, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I (rlogFC)", y="Visit II (rlogFC)", title="Inflamm") 
ggsave(file=file.path(fig_dir, "basket-Inflamm-two-visits-rlogFC.pdf"))

# composite score

############################################
# part II: basket composite score
############################################


#
# composite score by average TPM of basket genes: average over each genes' TPM z-score
#
baskets_avglog10TPM <- map_dfr(names(baskets_id), function(basket_name) {
    # per-gene z-score of TPM
    id <- baskets_id[[basket_name]]
    tpm <- assays(sanitized.dds)[["TPM"]][id, sanitized.dds$pheno_type == "FSHD" ] 

    # per-sample average TPM -> scaled by log10
    data.frame(avg_TPM=log10(colMeans(tpm)+1)) %>%
      rownames_to_column(var="sample_id") %>%
      dplyr::mutate(subject = str_replace(sample_id, "b$|b1$|-1$", "")) %>%
      left_join(visit, by="sample_id") %>%
      dplyr::select(-sample_id) %>%
      dplyr::filter(subject %in% paired_subject) %>%
      spread(key=visit, value=avg_TPM) %>%
      add_column(basket=basket_name)
}) %>%
   dplyr::mutate(basket = factor(basket, levels=c("DUX4-M6", "DUX4-M12", "Inflamm", "ECM")))

cor_basket_avglog10TPM <- baskets_avglog10TPM %>% group_by(basket) %>%
  summarize(cor = cor(I, II)) %>%
  dplyr::mutate(cor = paste0("cor=", format(cor, digits=2)))

ggplot(baskets_avglog10TPM, aes(x=`I`, y=`II`)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ basket, nrow=2, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_basket_avglog10TPM, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I: log10(avg. TPM+1)", y="Visit II: log10(avg. TPM+1)")
ggsave(file=file.path(fig_dir, "basket-avglog10TPM-corr.pdf"))   

#
# composite score by rlogFC
#

baskets_avg_rlogFC <- map_dfr(names(baskets_rlgFC), function(basket_name) {
  avg_rlogFC <- baskets_rlgFC[[basket_name]] %>% group_by(subject) %>%
    summarise(I=mean(I), II=mean(II)) %>%
    add_column(basket=basket_name)
}) %>% 
  dplyr::mutate(basket = factor(basket, levels=c("DUX4-M6", "DUX4-M12", "Inflamm", "ECM")))

cor_basket_avg_rlogFC <- baskets_avg_rlogFC %>% 
  group_by(basket) %>%
  summarize(cor = cor(I, II)) %>%
  dplyr::mutate(cor = paste0("cor=", format(cor, digits=2)))

ggplot(baskets_avg_rlogFC, aes(x=`I`, y=`II`)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap( ~ basket, nrow=2, scales="free") +
  theme_bw() +
  geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
  geom_text(data = cor_basket_avg_rlogFC, aes(label = cor), 
            x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="blue", size=2) +
  labs(x="Visit I: avg. rlogFC", y="Visit II: avg rlogFC")
ggsave(file=file.path(fig_dir, "basket-avg-rlogFC-corr.pdf"))   