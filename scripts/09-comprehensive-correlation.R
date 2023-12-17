# comprephensive correlation between variables

# This script suppports sections of Association of DUX4 scores to histopathology and clinical data

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggrepel))

#
# load data
#
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"

load(file.path(pkg_dir, "data", "mri.rda"))
load(file.path(pkg_dir, "data", "DUX4_positive.rda"))
load(file.path(pkg_dir, "data", "dds.rda"))

anno_gencode36 <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id") %>% # BiLat study using Gencode 36
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", ""))

# merge mri + dux4 score
merge_df <- mri %>% inner_join(DUX4_positive, by="sample_id") %>%
  dplyr::filter(!is.na(`Cumulative Score`), !is.na(`basket-M6`)) %>%
  dplyr::mutate(`DUX4+ (M6)` = if_else(`DUX4+ (M6)`, "DUX4+", "DUX4-"))


#
#  pathology vs. dux4
#
cor(merge_df$`Cumulative Score`, merge_df$`basket-M6`)
summary(lm(merge_df$`Cumulative Score` ~ merge_df$`basket-M6`))

merge_df %>% dplyr::mutate(`DUX4+ (M6)` = factor(`DUX4+ (M6)`, levels=c("DUX4-", "DUX4+"))) %>%
  group_by(`DUX4+ (M6)`) %>%
  summarise(mean=mean(`Cumulative Score`), sd=sd(`Cumulative Score`))
#
# CSS
#
css <- readxl::read_xlsx(file.path(pkg_dir, "extdata",
                                  "simple_Wellstone_CSS_age_gender_strength_2aug22.xlsx")) %>%
  dplyr::rename(Subject = `Record ID`) %>%
  dplyr::select(Subject, CSS)

css_dux4 <- DUX4_positive %>%
  dplyr::mutate(`DUX4+ (M6)` = if_else(`DUX4+ (M6)`, "DUX4+", "DUX4-")) %>% 
  dplyr::mutate(`DUX4+ (M6)` = factor(`DUX4+ (M6)`, levels=c("DUX4-", "DUX4+"))) %>%
  dplyr::mutate(Subject = str_replace(sample_id, "[L|R]", "")) %>%
  dplyr::left_join(css, by="Subject") 


cor(css_dux4$`CSS`, css_dux4$`basket-M6`) 
css_dux4 %>% 
  group_by(`DUX4+ (M6)`) %>%
  summarise(mean=mean(`CSS`), sd=sd(`CSS`))

#
# TA muscle strength
#
muscle_strength <- get(load(file.path(pkg_dir, "data", "muscle_strength.rda"))) %>%
  dplyr::rename(Subject=`Record ID`) %>%
  gather(key=`location`, value=`Foot Dorsiflexors`, -Subject) %>%
  dplyr::mutate(sample_id = paste0(Subject, 
                                   str_sub(location, start=1L, end=1L)))

muscle_strength %>% summarize(mean=mean(`Foot Dorsiflexors`),
                              sd = sd(`Foot Dorsiflexors`))
strength_dux4 <- DUX4_positive %>%
  dplyr::mutate(`DUX4+ (M6)` = if_else(`DUX4+ (M6)`, "DUX4+", "DUX4-")) %>% 
  dplyr::mutate(`DUX4+ (M6)` = factor(`DUX4+ (M6)`, levels=c("DUX4-", "DUX4+"))) %>%
  #dplyr::mutate(Subject = str_replace(sample_id, "[L|R]", "")) %>%
  dplyr::left_join(muscle_strength, by="sample_id") %>%
  dplyr::filter(!is.na(`Foot Dorsiflexors`))

# correlation: 
cor(strength_dux4$`Foot Dorsiflexors`, strength_dux4$`basket-M6`) 
strength_dux4 %>% group_by(`DUX4+ (M6)`) %>%
    summarize(mean=mean(`Foot Dorsiflexors`), sd=sd(`Foot Dorsiflexors`))


ggplot(strength_dux4, aes(x=`Foot Dorsiflexors`, y=`basket-M6`)) +
    geom_point() +
    theme_minimal() +
    geom_smooth(method="lm", se=FALSE) +
    scale_y_continuous(trans='log10') +
    theme(legend.position="bottom")
ggsave(file=file.path(pkg_dir, "figures", "muscle-strength-vs-dux4-score.pdf"),
       width=3.5, height=2.5)

pval <- summary(lm(`Foot Dorsiflexors` ~ `DUX4+ (M6)`, data=strength_dux4))
ggplot(strength_dux4, aes(x=`DUX4+ (M6)`, y=`Foot Dorsiflexors`)) +
  geom_boxplot(width=0.6) +
  geom_text(aes(label="p-val=0.0004", x=1.7, y=40), size=2) +
  #geom_jitter(width=0.3) +
  theme_minimal() +
  labs(x="")

ggsave(file=file.path(pkg_dir, "figures", "muscle-strength-vs-DUX4-status.pdf"),
       width=1.5, height=2)       

#
# Inflamm, ECM, complement, and IG basket score :gene_id, ens_id, gene_name in gencode v35
#

load(file.path(pkg_dir, "data", "all_baskets.rda"))

# note: later remove fat samples: "13-0009R", "13-0007R"
all_baskets_TPM <- map_dfr(names(all_baskets), function(name) {
  id <- all_baskets[[name]]$gencode_v35    
  assays(dds[id])[["TPM"]] %>% as.data.frame() %>%
      summarise(across(where(is.numeric), mean)) %>% t(.) %>% as.data.frame() %>%
      rownames_to_column(var="sample_name") %>%
      add_column(basket = name) %>%
      dplyr::rename(TPM = V1)
}) %>% 
  dplyr::mutate(basket = factor(basket, levels=c("DUX4-M", "DUX4-M6", "DUX4-M12", "Inflamm", "ECM", "Complement", "IG")),
                sample_id = str_replace(sample_name, "[b]*_.*", ""), .before="sample_name") %>%
  dplyr::select(-sample_name)                


# control basket
load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "data", "candidates.rda"))
controls_dds <- sanitized.dds[, sanitized.dds$pheno_type == "Control"]
control_baskets <- map_dfr(names(all_baskets), function(name) {
  id <- all_baskets[[name]]$gene_id_v88
  assays(controls_dds[id])[["TPM"]] %>% as.data.frame() %>%
      summarise(across(where(is.numeric), mean)) %>% t(.) %>% as.data.frame() %>%
      dplyr::rename(TPM = V1) %>%
      rownames_to_column(var="sample_id") %>%
      add_column(basket=name) 
}) %>%
  dplyr::mutate(basket = factor(basket, levels=c("DUX4-M", "DUX4-M6", "DUX4-M12", "Inflamm", "ECM", "Complement", "IG")))
 
save(control_baskets, file=file.path(pkg_dir, "data", "control_baskets.rda"))


#
# combine muscle strength, CSS, age, geneder, mri, bss, baskets (DUX4, inflamm, ECM, Complement, and  IG)
#       
spread_baskets <- all_baskets_TPM %>%
  spread(key=basket, value=TPM)

comprehensive_df <- mri %>% 
  #dplyr::filter(!is.na(`Cumulative Score`), !is.na(`basket-M6`)) %>%
  dplyr::full_join(css, by="Subject") %>%
  dplyr::full_join(muscle_strength %>% dplyr::select(-location, -Subject), by="sample_id") %>%
  dplyr::full_join(spread_baskets, by="sample_id") %>%
  dplyr::mutate(`DUX4+ (M6)` = if_else(`DUX4-M6` > 0.3, "DUX4+", "DUX4-")) %>%
  dplyr::mutate(`DUX4+ (M6)` = factor(`DUX4+ (M6)`, levels=c("DUX4-", "DUX4+")))


# A. muscle strength vs. STIR
comprehensive_df %>% dplyr::filter(!is.na(`Foot Dorsiflexors`), !is.na(STIR_status)) %>%
  group_by(STIR_status) %>%
  summarize(mean=mean(`Foot Dorsiflexors`), sd=sd(`Foot Dorsiflexors`))

comprehensive_df %>% dplyr::filter(!is.na(`Foot Dorsiflexors`), !is.na(STIR_status)) %>%
  ggplot(aes(x=STIR_status, y=`Foot Dorsiflexors`)) +
    geom_boxplot(width=0.6) +
    theme_minimal() +
    labs(x="")
ggsave(file=file.path(pkg_dir, "figures", "muscle-strength-vs-STIR-status.pdf"),
       width=1.5, height=2)     

# baskets vs. STIR : remore 13-0009R and 13-0007R
data <- all_baskets_TPM %>% 
  left_join(dplyr::select(mri, sample_id, STIR_status), by="sample_id") %>%
  add_row(control_baskets %>% add_column(STIR_status = "Control")) %>%
  dplyr::mutate(STIR_status = factor(STIR_status, levels=c("Control", "STIR-", "STIR+"))) %>%
  dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R")) %>%
  dplyr::filter(!basket == "DUX4-M12") %>%
  ggplot(aes(x=STIR_status, y=TPM)) +
    geom_boxplot(width=0.7, outlier.shape=NA, fill="grey75", alpha=0.5) +
    geom_jitter(width = 0.3, size=0.5, alpha=0.5) +
    facet_wrap(~ basket, scales="free_y", nrow=1) +
    theme_bw() +
    labs(x="", y="basket scores (TPM)") +
    scale_y_continuous(trans='log10') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggsave(file=file.path(pkg_dir, "figures", "STIR-vs-all-baskets.pdf"), width=6.2, height=2.8)  

# basekts scatter plots - exclude 13-0009R and 13-0007R
library(GGally)
spread_baskets %>% dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R")) %>%
  dplyr::select(`DUX4-M6`, Inflamm, ECM, Complement, IG) %>%
  ggpairs() + 
     scale_y_continuous(trans='log10') +
     scale_x_continuous(trans='log10') +
     theme_bw() 
ggsave(file=file.path(pkg_dir, "figures", "baskets-all-scatter-corr.pdf"), width=6, height=6)  

# comprehensive correlation (pearson) - figure; cannot do STIR_status
comprehensive_df %>% 
  dplyr::filter(!is.na(`Foot Dorsiflexors`), !is.na(STIR_status), !is.na(`DUX4-M6`)) %>%
  dplyr::select(`I. Variability in Fiber`, `II. Extent of Central Nucleation`, 
                `III. Necrosis/Regeneration`, `IV. Interstitial Fiobrsis`, 
                `Cumulative Score`, Fat_Infilt_Percent, FAT_FRACTION, STIR_RATING, CSS,
                `Foot Dorsiflexors`, 
                `DUX4-M6`, Inflamm, ECM, Complement, IG) %>%
  GGally::ggcorr(label = TRUE, label_round=1, label_size = 3, size=2, hjust=0.9, legend.position="none", layout.exp=3)

ggsave(file=file.path(pkg_dir, "figures", "comprehensive-correlation.pdf"), width=5, height=4)

#
# bilaterals comparisons: fat infilt, STIR_rating, pathology score, baskets, 
#
library(corrr)
# A. Fat infilt
variables <- c("I. Variability in Fiber", "II. Extent of Central Nucleation", 
               "III. Necrosis/Regeneration", "IV. Interstitial Fiobrsis",
               "Cumulative Score", "Fat_Infilt_Percent", "FAT_FRACTION", "STIR_RATING", 
               "Foot Dorsiflexors", "DUX4-M6", "Inflamm", "ECM", "Complement", "IG")
location_df <- map_dfr(variables, function(variable) {
        comprehensive_df %>% dplyr::select(Subject, !!sym(variable), location) %>%
          dplyr::filter(!is.na(!!sym(variable))) %>%
          tidyr::spread(key=location, value=!!sym(variable)) %>%
          dplyr::filter(!is.na(L), !is.na(R)) %>%
          add_column(var = variable)
      }) %>%
  dplyr::filter(!(var %in% c("DUX4-M6", "Inflamm", "ECM", "Complement", "IG") & Subject %in% c("13-0009", "13-0007"))) %>%
  dplyr::mutate(var = factor(var, levels = variables))
table(location_df$var)       

bilat_corr <- location_df %>%  group_by(var) %>%
  summarise(cor = cor(L, R)) %>%
  dplyr::mutate(cor = paste0("cor=", format(cor, digits=2)))


location_df %>% 
  dplyr::filter(var %in% variables[1:9]) %>%
  ggplot(aes(x=L, y=R)) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE) +
    facet_wrap( ~ var, nrow=3, scales="free") +
    geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
    geom_text(data = dplyr::filter(bilat_corr, var %in% variables[1:9]), aes(label = cor), 
              x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
    theme_bw()
ggsave(file=file.path(pkg_dir, "figures", "bilat-clinical-comparison.pdf"), width=5, height=6)

location_df %>% 
  dplyr::filter(var %in% variables[10:14]) %>%
  ggplot(aes(x=L, y=R)) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE) +
    facet_wrap( ~ var, nrow=2, scales="free") +
     scale_y_continuous(trans='log10') +
    scale_x_continuous(trans='log10') +
    geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
    geom_text(data = dplyr::filter(bilat_corr, var %in% variables[10:14]), aes(label = cor), 
              x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
    theme_bw()

ggsave(file=file.path(pkg_dir, "figures", "bilat-basket-comparison.pdf"), width=5, height=4)


#
# BSS
#
load(file.path(pkg_dir, "data", "methyl_clinic.rda"))
methyl_clinic <- methyl_clinic %>% 
  dplyr::mutate(haplotype = paste0(haplotype_1, "/", haplotype_2))
# subject 01-0022: use BSSA Q1 value
idx <- which(methyl_clinic$Subject == "01-0022")
methyl_clinic$`BSS Value`[idx] <- methyl_clinic$`BSSA Q1`[idx] 

# scatter plot between L and R
bss_df <- methyl_clinic %>% dplyr::select(Subject, location, `BSS Value`) %>%
  dplyr::filter(!is.na(`BSS Value`)) %>%
  dplyr::mutate(`BSS Value` = as.numeric(`BSS Value`)) %>%
  tidyr::spread(key=location, value=`BSS Value`) %>%
          dplyr::filter(!is.na(L), !is.na(R))
bss_cor <- cor(bss_df$L, bss_df$R)       
ggplot(bss_df, aes(x=L, y=R)) +
    geom_point(size=1) +
    geom_smooth(method="lm", se=FALSE) +
    geom_abline(slope=1, intercept=0, color="gray50", alpha=0.5, linetype="dashed") +
    geom_text(aes(label = format(bss_cor, digit=2)), 
              x=Inf, y=-Inf, hjust=1.1, vjust=-1.2, color="grey25", size=3) +
    theme_minimal()
ggsave(file=file.path(pkg_dir, "figures", "bilat-bss-comparison.pdf"), width=2, height=1.8)

# bss vs. DUX4 (remove 13-0009R/13-0007R); remove 01-0022?
merge_df <- comprehensive_df %>% dplyr::select(sample_id, `DUX4-M6`, STIR_status, `Cumulative Score`) %>%
  dplyr::filter(!is.na(`DUX4-M6`)) %>%
  dplyr::filter(!sample_id %in% c("13-0009R", "13-0007R")) %>%
  left_join(methyl_clinic, by="sample_id") %>%
  dplyr::filter(!Subject == "01-0022") %>%
  dplyr::filter(!is.na(`BSS Value`)) %>%
  dplyr::mutate(`BSS Value` = as.numeric(`BSS Value`)) 
  

cor(merge_df$`BSS Value`, merge_df$`DUX4-M6`)


ggplot(merge_df, aes(x=`BSS Value`, y=`DUX4-M6`)) +
    geom_point(size=1) +
    geom_smooth(method="lm", se=FALSE) +
    geom_text(aes(label = "0.19"), 
              x=Inf, y=Inf, hjust=1.1, vjust=1.2, color="grey25", size=3) +
    theme_minimal() 
#    scale_y_continuous(trans='log10') 
ggsave(file=file.path(pkg_dir, "figures", "BSS-vs-DUX4-M6.pdf"), width=2, height=1.8)

# bss vs. pathology
tmp <- merge_df %>% dplyr::filter(!is.na(`Cumulative Score`))
cor(tmp$`BSS Value`, tmp$`DUX4-M6`)
merge_df %>% dplyr::filter(!is.na(`Cumulative Score`)) %>%
  ggplot(aes(x=`BSS Value`, y=`Cumulative Score`)) +
    geom_point(size=1) +
    geom_smooth(method="lm", se=FALSE) +
    geom_text(aes(label = "0.17"), 
              x=Inf, y=Inf, hjust=1.1, vjust=1.2, color="grey25", size=3) +
    theme_minimal() 

ggsave(file=file.path(pkg_dir, "figures", "BSS-vs-pathology.pdf"), width=2, height=1.8)