#
# This script examines the relationship of RNA score (four markers) and mri data. 
# use logistic regression to predict whether the muscle is DUX4-targeted (DUX4+) based on MRI (STIR + FAT)
#  

library(DESeq2)
library(tidyverse)
library(ggrepel)


pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biops"
load(file.path(pkg_dir, "data", "dds.rda"))
load(file.path(pkg_dir, "data", "methyl_clinic.rda"))
load(file.path(pkg_dir, "data", "mri.rda"))
load(file.path(pkg_dir, "data", "rna_score.rda"))
load(file.path(pkg_dir, "data", "blood_fat_muscle_content.rda"))

blood_fat_muscle_content <- blood_fat_muscle_content %>%
  dplyr::rename(RNA_sample_id = sample_name)
annotation <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id") %>%
  dplyr::mutate(ensembl = sapply(str_split(gene_id, "[.]"), "[", 1))

merge_df <- rna_score %>% left_join(blood_fat_muscle_content, by="RNA_sample_id") %>%
  inner_join(select(mri, -Subject), by="sample_id") %>%
  dplyr::mutate(STIR_RATING = factor(STIR_RATING)) %>%
  dplyr::mutate(DUX4_targeted = if_else(DUX4_targeted, 1, 0)) %>%
  dplyr::filter(!sample_id %in% c("13-0007R", "13-0009R"))

# 1. use four biomarkers here

#
# rna_score vs. STIR-/+ and GLM
#

ggplot(merge_df, aes(x=STIR_status, y=tpmlog_4)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=DUX4_targeted), size=1, width=0.3) +
  theme_bw()
ggsave(file.path(pkg_dir, "figures", "test.pdf"), width=3, height=3)  

library(aod)
bistir_logit <- glm(DUX4_targeted ~ STIR_status, data = merge_df, family = "binomial")
summary(bistir_logit)
exp(coef(bistir_logit))
aod::wald.test(b = coef(bistir_logit), Sigma = vcov(bistir_logit), Terms=2)


stir_logit <- glm(DUX4_targeted ~ STIR_RATING, data = merge_df, family = "binomial")
summary(stir_logit)
aod::wald.test(b = coef(stir_logit), Sigma = vcov(stir_logit), Terms=c(2:5))
exp(coef(stir_logit))

# using STIR-/+ is better to perdict DUX4-targeted muscle
new_df = data.frame(STIR_status=factor(c("STIR-", "STIR+")))
predict(bistir_logit, newdata = new_df, type = "response")

new_df = data.frame(STIR_RATING=factor(c(0:4)))
predict(stir_logit, newdata = new_df, type = "response")

#
# rna_score vs. Fat_Infil_percent
#
merge_df %>% dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R")) %>%
ggplot(aes(x=Fat_Infilt_Percent, y=tpmlog_4)) +
  geom_point(size=1, aes(color=DUX4_targeted)) +
  geom_smooth(method="loess") +
  theme_minimal() +
  theme(legend.position="none") +
  geom_text_repel(aes(label=sample_id), size=1.8, show.legend=FALSE) 

ggsave(file.path(pkg_dir, "figures", "test.pdf"), width=4, height=3)  

merge_df %>% dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R")) %>%
ggplot(aes(x=FAT_FRACTION, y=tpmlog_4)) +
  geom_point(size=1, aes(color=DUX4_targeted)) +
  geom_smooth(method="loess") +
  theme_minimal() +
  theme(legend.position="none") +
  geom_text_repel(aes(label=sample_id), size=1.8, show.legend=FALSE) 

ggsave(file.path(pkg_dir, "figures", "test.2.pdf"), width=4, height=3)  

merge_df %>% dplyr::filter(Subject=="13-0010")

#
# logistics: 
#
tmp <- merge_df %>% dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R"))
fat_logit <- glm(DUX4_targeted ~ Fat_Infilt_Percent, data = tmp, family = "binomial")
summary(fat_logit)
new_df <- data.frame(Fat_Infilt_Percent = seq(0, 1, by=0.01))
new_df$prob = predict(fat_logit, newdata = new_df, type = "response")
ggplot(new_df, aes(x=Fat_Infilt_Percent, y=prob)) +
  geom_line() +
  theme_minimal()
ggsave(file.path(pkg_dir, "figures", "test.pdf"), width=4, height=3)

tmp <- merge_df %>% dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R"))
fat_logit2 <- glm(DUX4_targeted ~ FAT_FRACTION, data = tmp, family = "binomial")
summary(fat_logit2)
new_df <- data.frame(FAT_FRACTION = seq(0, 1, by=0.01))
new_df$prob = predict(fat_logit2, newdata = new_df, type = "response")
ggplot(new_df, aes(x=FAT_FRACTION, y=prob)) +
  geom_line() +
  theme_minimal()
ggsave(file.path(pkg_dir, "figures", "test.2.pdf"), width=4, height=3)


# logistic regresssion:  + STIR-/+
tmp <- merge_df %>% dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R"))
fat_stir_logit <- glm(DUX4_targeted ~ Fat_Infilt_Percent + STIR_status, data = tmp, family = "binomial")
summary(fat_stir_logit)
new_df <- data.frame(Fat_Infilt_Percent = rep(seq(0, 1, by=0.01), 2), 
                     STIR_status = rep(c("STIR-", "STIR+"), each=length(seq(0, 1, by=0.01))))
new_df$prob = predict(fat_stir_logit, newdata = new_df, type = "response")
ggplot(new_df, aes(x=Fat_Infilt_Percent, y=prob, group=STIR_status, color=STIR_status)) +
  geom_line() +
  theme(legend.position="top") +
  theme_minimal()
ggsave(file.path(pkg_dir, "figures", "test.2.pdf"), width=4, height=3)

                    
