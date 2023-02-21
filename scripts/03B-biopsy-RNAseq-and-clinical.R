# 03B-biopsy-markers-vs-clinical.R
# 
# 1) inter and intrapersonal variation in clinial/MRI scores


library(DESeq2)
library(tidyverse)
library(ggrepel)


pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
load(file.path(pkg_dir, "data", "dds.rda"))
load(file.path(pkg_dir, "data", "methyl_clinic.rda"))
load(file.path(pkg_dir, "data", "t1.rda"))

annotation <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id") %>%
  dplyr::mutate(ensembl = sapply(str_split(gene_id, "[.]"), "[", 1))

#
# 1) MRI: STIR and FAT - trajectory plot; only show label ones with changes
#
stir_diff <- methyl_clinic %>% dplyr::select(Subject, location, STIR) %>%
  tidyr::drop_na(STIR, location) %>% 
  dplyr::mutate(STIR = as.numeric(STIR)) %>%
  dplyr::filter(!grepl("2", location)) %>%
  tidyr::spread(key=location, value=STIR) %>%
  dplyr::mutate(delta= abs(R-L)) %>% 
  dplyr::filter(delta > 0) %>% pull(Subject)
# "01-0030" "13-0003" "13-0004" "32-0022" "32-0030"  

methyl_clinic %>% tidyr::drop_na(STIR, location) %>%
  dplyr::filter(!grepl("2", location)) %>%
  dplyr::mutate(label_change = ifelse(Subject %in% stir_diff, Subject, NA)) %>%
  ggplot(aes(x=location, y=STIR)) +
    geom_line(aes(group=Subject), color="gray50", show.legend=FALSE) +
    geom_text_repel(aes(label=label_change), segment.curvature = -1e-20, nudge_x=0.2, nudge_y=0.2, size=2, show.legend=FALSE) +
    geom_jitter(width=0.08, height=0.08, size=1, color="steelblue") +
    theme_bw() 

ggsave(file.path(pkg_dir, "figures", "MRI-paired-trajectory.pdf"), width=3, height=4)  
# not many have STIR ans FAT changes


#
# 2) FAT: FAT - trajectory plot; only show label ones with changes
#
fat_diff <- methyl_clinic %>% dplyr::select(Subject, location, FAT) %>%
  tidyr::drop_na(FAT, location) %>% 
  #dplyr::mutate(FAT = as.numeric(STIR)) %>%
  dplyr::filter(!grepl("2", location)) %>%
  tidyr::spread(key=location, value=FAT) %>%
  dplyr::mutate(delta= abs(as.numeric(R)-as.numeric(L))) %>% 
  dplyr::filter(delta > 0) %>% pull(Subject)
#  "13-0003" "32-0020" "32-0022" "32-0025" "32-0027" "32-0028" "32-0029"  

methyl_clinic %>% tidyr::drop_na(FAT, location) %>%
  dplyr::filter(!grepl("2", location)) %>%
  dplyr::mutate(label_change = ifelse(Subject %in% fat_diff, Subject, NA)) %>%
  ggplot(aes(x=location, y=FAT)) +
    geom_line(aes(group=Subject), color="gray50", show.legend=FALSE) +
    geom_text_repel(aes(label=label_change), segment.curvature = -1e-20, nudge_x=0.2, nudge_y=0.2, size=2, show.legend=FALSE) +
    geom_jitter(width=0.08, height=0.08, size=1, color="steelblue") +
    theme_bw() 

ggsave(file.path(pkg_dir, "figures", "FAT-paired-trajectory.pdf"), width=3, height=4)  

# T1 fraction counterlateral pairs
t1_diff <- t1 %>% dplyr::select(Subject, location, FAT_FRACTION) %>%
  tidyr::drop_na(FAT_FRACTION, location) %>%
  tidyr::spread(key=location, value=FAT_FRACTION) %>%
  dplyr::mutate(delta= abs(as.numeric(R)-as.numeric(L))) %>% 
  dplyr::filter(delta > 0.1) %>% pull(Subject)

# [1] "13-0006" "13-0008" "13-0010" "32-0020" "32-0025" "32-0027" "32-0028" "32-0029"  
t1 %>% tidyr::drop_na(FAT_FRACTION, location) %>%
  dplyr::mutate(label_change = ifelse(Subject %in% t1_diff, Subject, NA)) %>%
  ggplot(aes(x=location, y=FAT_FRACTION, group=Subject)) +
    geom_line(color="gray50", show.legend=FALSE) +
    geom_text_repel(aes(label=label_change), segment.curvature = -1e-20, nudge_x=0.05, nudge_y=0.05, size=2, show.legend=FALSE) +
    theme_bw()
ggsave(file.path(pkg_dir, "figures", "FAT-FRACTION-paired-trajectory.pdf"), width=3, height=4)  

  
#
# 3) Methyl: dynamic between  counterlateral pairs
#

methyl_clinic %>% tidyr::drop_na(`BSS Value`) %>%
  dplyr::filter(!grepl("2", location)) %>%
  ggplot(aes(x=location, y=`BSS Value`)) +
    geom_line(aes(group=Subject), color="gray50", show.legend=FALSE) +
    #geom_text_repel(aes(label=label_change), segment.curvature = -1e-20, nudge_x=0.2, nudge_y=0.2, size=2, show.legend=FALSE) +
    #geom_jitter(width=0.08, height=0.08, size=1, color="steelblue") +
    geom_point(size=1, color="steelblue") +
    theme_bw() 
ggsave(file.path(pkg_dir, "figures", "BSS-paired-trajectory.pdf"), width=3, height=4)  

#
# 4) BSS vs. Fat vs. STIR vs. pathology
#    density plot for each
#
methyl_clinic %>% tidyr::drop_na(`BSS Value`, `STIR`) %>% 
  ggplot(aes(x=`BSS Value`)) +
    geom_density() +
    theme_bw()
ggsave(file.path(pkg_dir, "figures", "BSS_density.pdf"), width=3, height = 2.5)    

#
# 4a) BSS vs. Fat
#
x = methyl_clinic %>% tidyr::drop_na(`BSS Value`, `FAT`)
pval = anova(lm(`BSS Value` ~ `FAT`, x)) # 0.9

methyl_clinic %>% tidyr::drop_na(`BSS Value`, `FAT`, location) %>% 
  dplyr::mutate(location = str_replace(location, "2", "")) %>%
  ggplot(aes(x=FAT, y=`BSS Value`)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    facet_wrap(~location) +
    theme_bw()
ggsave(file.path(pkg_dir, "figures", "BSS-vs-FAT-boxplot.pdf"), width=5, height=3)    

x=methyl_clinic %>% inner_join(t1, by="sample_id") %>%  tidyr::drop_na(`BSS Value`, `FAT_FRACTION`, location.x)
lm(`BSS Value` ~ FAT_FRACTION, x) # p-value = 0.7, cor=-0.04
methyl_clinic %>% inner_join(t1, by="sample_id") %>%
  tidyr::drop_na(`BSS Value`, `FAT_FRACTION`, location.x) %>%
  dplyr::mutate(location = str_replace(location.x, "2", "")) %>%
  ggplot(aes(x=FAT_FRACTION, y=`BSS Value`)) +
    geom_point() +
    facet_wrap(~location.x) +
    theme_bw()
ggsave(file.path(pkg_dir, "figures", "BSS-vs-FAT_FRACTION-scatter.pdf"), width=5, height=3)    

#
# 4b) BSS vs. STIR 
#

x = methyl_clinic %>% tidyr::drop_na(`BSS Value`, STIR, location) %>%
  dplyr::mutate(location = str_replace(location, "2", "")) %>% 
  dplyr::filter(location == "L")
pval = anova(lm(`BSS Value` ~ STIR, x)) # 0.1

x = methyl_clinic %>% tidyr::drop_na(`BSS Value`, STIR, location) %>%
  dplyr::mutate(location = str_replace(location, "2", "")) %>% 
  dplyr::filter(location == "R")
pval = anova(lm(`BSS Value` ~ STIR, x)) # 0.05

methyl_clinic %>% tidyr::drop_na(`BSS Value`, STIR, location) %>% 
  dplyr::mutate(location = str_replace(location, "2", "")) %>%
  ggplot(aes(x=STIR, y=`BSS Value`)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    facet_wrap(~location) +
    theme_bw()
ggsave(file.path(pkg_dir, "figures", "BSS-vs-STIR-boxplot.pdf"), width=5, height=3)    

#
# 4c) BSS vs. pathology
#
methyl_clinic %>% tidyr::drop_na(`BSS Value`, `Pathology Score`, location) %>% 
  dplyr::mutate(location = str_replace(location, "2", "")) %>%
  ggplot(aes(x=`Pathology Score`, y=`BSS Value`)) +
    geom_point() +
    facet_wrap(~location) +
    theme_bw()
ggsave(file.path(pkg_dir, "figures", "BSS-vs-Pathology-scatter.pdf"), width=5, height=3) 
summary(lm(`Pathology Score` ~ `BSS Value`, methyl_clinic)) # p-value - 0.26

#
# 5) RNA-seq score: TPM and rlog
#
  
load(file.path(pkg_dir, "data", "FSHD_markers.rda"))
load(file.path(pkg_dir, "data", "four_markers.rda"))
load(file.path(pkg_dir, "data", "rlog.rda"))


# 5a RNA_score (sum of TPM)
markers_4 <- annotation %>% dplyr::filter(gene_name %in% four_markers$gene_name)
fshd_markers <- annotation %>% dplyr::filter(ensembl %in% FSHD_markers$ensembl)
fshd_markers <- fshd_markers[fshd_markers$gene_id %in% rownames(rlog),] # total 52 markers that have some expression

# which FSHD_markers are not included in "annotation"? TPM here are in log10 scale
rna_score <- data.frame(tpm_4_sum = assays(dds)[["TPM"]][markers_4$gene_id, ] %>% colSums(),
                        tpmlog_4 = log10(assays(dds)[["TPM"]][markers_4$gene_id, ] + 1) %>% colSums(),
                        tpmlog_fshd = log10(assays(dds)[["TPM"]][fshd_markers$gene_id, ] +1) %>% colSums(),
                        rlog_4 = assay(rlog)[markers_4$gene_id, ] %>% colSums(),
                        rlog_fshd = assay(rlog)[fshd_markers$gene_id, ] %>% colSums() ) %>%
  rownames_to_column(var = "RNA_sample_id") %>%                        
  dplyr::mutate(sample_id = sapply(str_split(RNA_sample_id, "_"), "[", 1), .after="RNA_sample_id")   %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "b", "")) %>%
  dplyr::mutate(Subject = str_replace(sample_id, "[L|R]", "")) %>%
  dplyr::mutate(DUX4_targeted = tpm_4_sum > 1)             

save(rna_score, file=file.path(pkg_dir, "data", "rna_score.rda"))

# TPM_4 density plot
ggplot(rna_score, aes(x=tpm_4_sum)) +
  geom_density() +
  theme_minimal()
ggsave(file.path(pkg_dir, "figures",  "test.pdf"), width=4, height=3)  

#
# 6) rna score vs. all other metrics
#
fig_dir <- file.path(pkg_dir, "figures")
merge_df <- rna_score %>%
  dplyr::inner_join(methyl_clinic, by="sample_id") %>%
  dplyr::inner_join(dplyr::select(t1, -location, -Subject), by="sample_id") %>%
  dplyr::mutate(STIR_status =if_else(STIR=="0", "STIR-", "STIR+"))

# 6a) ran score vs. STIR
merge_df %>% tidyr::drop_na(STIR) %>%
  ggplot(aes(x=STIR, y=tpmlog_4)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    facet_wrap(~location) +
    labs(y="sum of log10(TPM)") +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-vs-stir-tpm4.pdf"), width=5, height=3)    

merge_df %>% tidyr::drop_na(STIR) %>%
  ggplot(aes(x=STIR, y=tpm_fshd)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    facet_wrap(~location) +
    labs(y="sum of log10(TPM) of 52 markers") +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-vs-stir-tpm52.pdf"), width=5, height=3)   

merge_df %>% tidyr::drop_na(STIR_status) %>%
  ggplot(aes(x=STIR_status, y=tpm_fshd)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    facet_wrap(~location) +
    labs(y="sum of log10(TPM) of 52 markers") +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-vs-stir-status-tpm52.pdf"), width=5, height=3)   

merge_df %>% tidyr::drop_na(STIR) %>%
  ggplot(aes(x=STIR, y=rlog_4)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    facet_wrap(~location) +
    labs(y="RNA score (rlog)") +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-vs-stir-rlog4.pdf"), width=5, height=3)  

x = merge_df %>% tidyr::drop_na(STIR)
anova(lm(rlog_fshd ~ STIR_status, x)) # 1.7e-5

# 6b) ran score vs. FAT
merge_df %>% tidyr::drop_na(FAT) %>%
  ggplot(aes(x=FAT, y=rlog_4)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    facet_wrap(~location) +
    labs(y="RNA score (rlog)") +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-vs-FAT-rlog4.pdf"), width=5, height=3)  
x = merge_df %>% tidyr::drop_na(FAT)
anova(lm(rlog_fshd ~ FAT, x))

merge_df %>% tidyr::drop_na(FAT) %>%
  ggplot(aes(x=FAT, y=tpm_4)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    facet_wrap(~location) +
    labs(y="RNA score (TPM 52)") +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-vs-FAT-tpm4.pdf"), width=5, height=3)  


# 6c) ran score vs. FAT_FRACTION
merge_df %>% tidyr::drop_na(FAT_FRACTION) %>%
  ggplot(aes(x=FAT_FRACTION, y=tpm_4)) +
    geom_point() +
    geom_smooth(method="loess") +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    labs(y="RNA score (tpm_4)") +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-vs-FAT_FRACTION-tpm4.pdf"), width=5, height=3)  

merge_df %>% tidyr::drop_na(FAT_FRACTION) %>%
  ggplot(aes(x=FAT_FRACTION, y=tpm_fshd)) +
    geom_point() +
    geom_smooth(method="loess") +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    labs(y="RNA score (tpm_52)") +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-vs-FAT_FRACTION-tpm52.pdf"), width=5, height=3)

x = merge_df %>% tidyr::drop_na(FAT_FRACTION)
anova(lm(tpm_fshd ~ FAT_FRACTION, x))

#
# how about tpm_fhsd vs. fat_content? fat_content vs. FAT_FRACTION?
#
load(file.path(pkg_dir, "data", "blood_fat_muscle_content.rda"))
x = blood_fat_muscle_content %>% 
  dplyr::rename(RNA_sample_id = sample_name) %>%                        
  left_join(merge_df, by="RNA_sample_id")
cor(x$FAT_FRACTION, x$fat)

x %>% tidyr::drop_na(FAT_FRACTION, fat) %>%
  ggplot(aes(x=FAT_FRACTION, y=fat)) +
    geom_point() +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    labs(y="fat content (log10TPM)") +
    theme_bw()
ggsave(file.path(fig_dir, "rna_fat-vs-FAT_FRACTION.pdf"), width=5, height=3)

# FAT vs. FAT_FRACTION
x %>% tidyr::drop_na(FAT_FRACTION, FAT) %>%
  ggplot(aes(y=FAT_FRACTION, x=FAT)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter() +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    labs(y="FAT_FRACTION") +
    theme_bw()
ggsave(file.path(fig_dir, "FAT-vs-FAT_FRACTION.pdf"), width=5, height=3)

# FAT vs. fat content
x %>% tidyr::drop_na(fat, FAT) %>%
  ggplot(aes(y=fat, x=FAT)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter() +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    labs(y="fat content (RNA TPM") +
    theme_bw()
ggsave(file.path(fig_dir, "FAT-vs-fat-content-tpm.pdf"), width=5, height=3)

# fat content vs. RNA score: do they have any relationship?
x %>% tidyr::drop_na(fat, tpm_fshd) %>%
  ggplot(aes(x=fat, y=tpm_fshd)) +
    geom_point() +
    geom_smooth(lm="loess") +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    labs(y="RNA score (TPM 52)", x="fat content (FASN, LEP, SCD)") +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-vs-rna_fat.pdf"), width=5, height=3)

# FAT vs

# 6c) ran score vs. BSS
merge_df %>% tidyr::drop_na(`BSS Value`) %>%
  ggplot(aes(x=`BSS Value`, y=rlog_4)) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE) +
    labs(y="RNA score (rlog 4)") +
    facet_wrap(~location) +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-bs-BSS-rlog4.pdf"), width=5, height=3)
summary(lm(`BSS Value` ~ rlog_4, merge_df)) # p-value = 0.5


# 6d) ran score vs. pathology
merge_df %>% tidyr::drop_na(`Pathology Score`, location) %>%
  ggplot(aes(x=`Pathology Score`, y=tpm_4)) +
    geom_point() +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    geom_smooth(method="lm", se=FALSE) +
    labs(y="RNA score (TPM 4)") +
    facet_wrap(~location) +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-bs-pathology-tpm4.pdf"), width=5, height=3)

x=merge_df %>% tidyr::drop_na(`Pathology Score`)
cor(x$`Pathology Score`, x$rlog_4) # corr = 0.5

# 6e) rna score vs. inflammation score
merge_df %>% tidyr::drop_na(`Inflammation Score`, location) %>%
  dplyr::mutate(`Inflammation Score` = factor(`Inflammation Score`)) %>%
  ggplot(aes(x=`Inflammation Score`, y=rlog_4)) +
    geom_boxplot(oulier.shape=NA) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE) +
    labs(y="RNA score (rlog 4)") +
    facet_wrap(~location) +
    theme_bw()
ggsave(file.path(fig_dir, "rnascore-bs-inflammation-rlog.pdf"), width=5, height=3)

x=merge_df %>% tidyr::drop_na(`Inflammation Score`)
cor(x$`Inflammation Score`, x$rlog_4) # corr = 0.44

#
# PCA for all matrices
#