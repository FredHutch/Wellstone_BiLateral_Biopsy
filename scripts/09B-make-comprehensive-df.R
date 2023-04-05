
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

load(file.path(pkg_dir, "data", "mri.rda"))
load(file.path(pkg_dir, "data", "DUX4_positive.rda"))
load(file.path(pkg_dir, "data", "dds.rda"))
load(file.path(pkg_dir, "data", "all_baskets.rda"))
load(file.path(pkg_dir, "data", "blood_fat_muscle_content.rda"))

#
# tidy-up some data frames
#
anno_gencode36 <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id") %>% # BiLat study using Gencode 36
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", ""))


# css / age/sex
css <- readxl::read_xlsx(file.path(pkg_dir, "extdata",
                                  "simple_Wellstone_CSS_age_gender_strength_2aug22.xlsx")) %>%
  dplyr::rename(Subject = `Record ID`) %>%
  dplyr::select(Subject, CSS, `Age at this visit`, Sex)

#
# muscle strength
#
muscle_strength <- get(load(file.path(pkg_dir, "data", "muscle_strength.rda"))) %>%
  dplyr::rename(Subject=`Record ID`) %>%
  gather(key=`location`, value=`Foot Dorsiflexors`, -Subject) %>%
  dplyr::mutate(sample_id = paste0(Subject, 
                                   str_sub(location, start=1L, end=1L)))
#
# blood_fat_muscle_content (sum(log10(TPM+1)))
#
blood_fat_muscle_content <- blood_fat_muscle_content %>%
  dplyr::select(-RNA_sample_id) %>%
  rename(sample_id = sample_name)

# baskets score; not going to include
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

#
# baskets score in log10 scale: logSum = sum of log10(basket genes)
#
all_baskets_logTPM <- sapply(names(all_baskets), function(name) {
  id <- all_baskets[[name]]$gencode_v35 
  tpm_score <- colSums(log10(assays(dds[id])[["TPM"]]+1))
}) %>% as.data.frame() %>%
  rename_with(~paste0(.x, "-logSum")) %>%
  rownames_to_column(var="sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "[b]*_.*", "")) 


#
# complement scores
#

# 1. tidy up the data frame to be merge to comprehensive_df
complement_scoring <- readxl::read_xlsx(file.path(pkg_dir, "extdata", 
                                          "FINAL.Revised.BilatStudy.Path.MACstaining data.FINAL.xlsx"),
                                sheet=1) %>%
  dplyr::select(`Record ID`, `Complement Scoring (Left)`, `Complement Scoring (Right)`) %>%
  tidyr::gather(key=Location, value=`Complement Scoring`, -`Record ID`)   %>%
  dplyr::mutate(Location = if_else(str_detect(Location, "Left"), "L", "R")) %>%
  dplyr::mutate(sample_id = paste0(`Record ID`, Location)) %>%
  dplyr::select(sample_id, `Complement Scoring`)    


#
# comprehensive covariates data frame (not include BSS yet) to include:
# 1. muscle strength 
# 2. basket score (logTPM: logSum over basket genes)
# 3. blood, fat, muscle content
# 4. complement scoring (using the final scoring)
#

#spread_baskets <- all_baskets_TPM %>%
#  spread(key=basket, value=TPM)

comprehensive_df <- mri %>% 
  dplyr::full_join(css, by="Subject") %>%
  dplyr::full_join(muscle_strength %>% dplyr::select(-location, -Subject), by="sample_id") %>%
  #dplyr::full_join(spread_baskets, by="sample_id") %>%
  #dplyr::mutate(`DUX4+ (M6)` = if_else(`DUX4-M6` > 0.3, "DUX4+", "DUX4-")) %>%
  #dplyr::mutate(`DUX4+ (M6)` = factor(`DUX4+ (M6)`, levels=c("DUX4-", "DUX4+"))) %>%
  dplyr::full_join(all_baskets_logTPM, by="sample_id") %>%
  dplyr::full_join(blood_fat_muscle_content, by="sample_id") %>%
  dplyr::full_join(complement_scoring, by="sample_id") %>%
  dplyr::relocate(`Complement Scoring`, .after="Inflammation Score") %>%
  dplyr::relocate(`CSS`, .after=`Complement Scoring`) %>%
  dplyr::relocate(`Foot Dorsiflexors`, .after=`CSS`) %>%
  dplyr::relocate(`Age at this visit`, .after=`sample_id`) %>%
  dplyr::relocate(`Sex`, .after=`Age at this visit`)

#
# add class - that classified by ML with longitudinal training data
#
load(file.path(pkg_dir, "data", "bilat_MLpredict.rda"))
comprehensive_df <- comprehensive_df %>%
  left_join(dplyr::select(bilat_MLpredict, sample_id, class), by="sample_id")

save(comprehensive_df, file = file.path(pkg_dir, "data", "comprehensive_df.rda"))
writexl::write_xlsx(comprehensive_df, path=file.path(pkg_dir, "stats", "comprehensive_metadata.xlsx"))



#
# control baskets
#

load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
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

control_baskets_logsum <- map_dfr(names(all_baskets), function(name) {
  id <- all_baskets[[name]]$gene_id_v88
  log10(assays(controls_dds[id])[["TPM"]] + 1) %>%
      as.data.frame() %>%
      summarise(across(where(is.numeric), sum)) %>% t(.) %>% as.data.frame() %>%
      dplyr::rename(log10TPM = V1) %>%
      rownames_to_column(var="sample_id") %>%
      add_column(basket=name) 
}) %>%
  dplyr::mutate(basket = factor(basket, levels=c("DUX4-M", "DUX4-M6", "DUX4-M12", "Inflamm", "ECM", "Complement", "IG")))
 
save(control_baskets_logsum, file=file.path(pkg_dir, "data", "control_baskets_logsum.rda"))
