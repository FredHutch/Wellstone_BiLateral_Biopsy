# 00-data-cleaning.R: clean up clinical/MRI/Methyl data and make a composite dataset
#
# 1. Generate bogus IDs.
# 2. Convert xlsx and csv sheet to data.frame instance.
# 3. dataset (/data): 
#    - bogus_id.rda: bogus ids
#    - methyl_clinic.rda: bi-sulfi methylation, MIR, and pathology scores
#
#
# note: I add the last two columns of Wellstone BiLat_Muscle Biopsy_Jones MethylSeq Data.xlsx
#       and bumped the version to .v2

library(tidyverse)
library(readxl)
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
ext_dir <- file.path(pkg_dir, "extdata")
bs_methyl_filename <- file.path(ext_dir, "Wellstone BiLat_Muscle Biopsy_Jones MethylSeq Data.v2.xlsx")
histology_filename <- file.path(ext_dir, "BiLAT Biopsy Data 1.14.2022.csv")
# muscle strength: Left/Right Food Dorsiflexors - EX and EY by volume 
baseline_filename <- file.path(ext_dir, "BILLATStudy_Baseline_29March2022.xlsx")
mri_filename <- file.path(ext_dir, "T1_STIR_ratings_bilat_sample_26july21.xlsx")
#summary_filename <- file.path(ext_dir, "WELLSTONE_summary_data_28mar22_send.xlsx")
summary_filename <- file.path(ext_dir, "WELLSTONE_summary_data_30may22_sendC.xlsx")

#
# bogus subject ID
#
bogus_id <- readxl::read_xlsx(path=bs_methyl_filename,
                              sheet=1, skip=1, range=readxl::cell_cols(1)) %>%
  dplyr::rename(Subject = `Subject #`)  %>%                     
  dplyr::select(Subject) %>% tidyr::drop_na() %>%
  dplyr::distinct(Subject) %>%
  dplyr::mutate(bogus_subject = c(paste0("00", 1:9), paste0("0", 10:length(Subject)))) %>%
  dplyr::mutate(bogus_subject= paste0("F", bogus_subject))
save(bogus_id, file=file.path(pkg_dir, "data", "bogus_id.rda"))

#
# Methylation data: data cleaning
#
bs_methyl <- read_xlsx(path=bs_methyl_filename,
                     sheet=1, skip=1) %>%
  dplyr::rename(Subject = `Subject #`)  %>%
  dplyr::select(Subject, Sex, `Predicted Haplotype`, `Muscle Biopsy`, ADF, FAT, STIR, 
               `BSSA Q1`, `BSSA Q2`, `BSSL Q1`, `BSSL Q2`, `BSSX Ave`,
               `FSHD1 FSHD2 Healthy`, `BSS Value`, `use BSSL`) %>%
  dplyr::mutate(location = if_else(str_detect(`Muscle Biopsy`, "R"), "R", "L")) %>%
  dplyr::mutate(location = if_else(str_detect(`Muscle Biopsy`, "[*][*]"), 
                                   paste0(location, "2"), location)) %>%     
  dplyr::mutate(`use BSSL` = if_else(`use BSSL` == "T", TRUE, FALSE)) %>%
  dplyr::mutate(`use BSSL` = if_else(is.na(`use BSSL`), FALSE, `use BSSL`))   %>%
  dplyr::mutate(ADF = factor(ADF, levels=c("1", "2", "3-", "3", "4-", "4", "4+", "5-", "5"))) %>%
  dplyr::mutate(FAT = factor(FAT), STIR=factor(STIR))  

# (1.) clean up 
#      - a. fill up missing components
#      - b. add haplotype
#      - c. remove 32-0029R and make 32-0029 R** to be 32-0029
#      - d. remove other duplicates *-*L2 or *-*R2

# (a) remove empty row
bs_methyl <- bs_methyl[-c(25, 48), ]
# (a) fill up missing components
fill_missing_values_based_on_subject <- function(df, col_name) {
    df[seq(from=2, to=64, by=2), col_name] <- df[seq(from=1, to=63, by=2), col_name, drop=TRUE]
    df[65, col_name] <- df[64, col_name, drop=TRUE]
    df[c(67:69), col_name] <- df[66, col_name, drop=TRUE]
    df[71, col_name] <- df[70, col_name, drop=TRUE]
    df
}

bs_methyl <- fill_missing_values_based_on_subject(bs_methyl, "Subject")
bs_methyl <- fill_missing_values_based_on_subject(bs_methyl, "Predicted Haplotype")
bs_methyl <- fill_missing_values_based_on_subject(bs_methyl, "Sex")
bs_methyl <- fill_missing_values_based_on_subject(bs_methyl, "FSHD1 FSHD2 Healthy")


#
# (c) haplotype
#
bs_methyl <- bs_methyl %>%
  dplyr::mutate(haplotype_chr4 = sapply(str_split(`Predicted Haplotype`, ","), "[[", 1)) %>%
  dplyr::mutate(haplotype_1 = sapply(str_split(haplotype_chr4, "/"), "[[", 1),
                haplotype_2 = sapply(str_split(haplotype_chr4, "/"), "[[", 2)) %>%
  dplyr::mutate(haplotype_1 = sapply(str_split(haplotype_1, "[(]"), "[[", 1),
                          haplotype_2 = sapply(str_split(haplotype_2, "[(]"), "[[", 1)) %>%
  dplyr::mutate(haplotype_1 = factor(haplotype_1), 
                haplotype_2 = factor(haplotype_2))    

# (d/e) 
bs_methyl <- bs_methyl %>%
  dplyr::filter(!is.na(`BSS Value`)) %>%
  dplyr::mutate(location = if_else(Subject=="32-0029" & location=="R2", "R", location)) %>%
  dplyr::mutate(sample_id = paste0(Subject, location)) %>%
  dplyr::filter(!location %in% c("R2", "L2")) %>%
  dplyr::relocate(location, sample_id, .after="Subject")

bs_methyl %>% dplyr::filter(Subject=="32-0029") %>% as.data.frame()  

save(bs_methyl, file=file.path(pkg_dir, "data", "bs_methyl.rda"))

# (3.) add bogus ID
#bs_methyl <- bs_methyl %>%
#  left_join(bogus_id, by="Subject") %>%
#  dplyr::mutate(bogus_sample_id = paste0(bogus_subject, location), sample_id = paste0(Subject, location)) %>%
#  dplyr::relocate(sample_id, bogus_subject, bogus_sample_id, .after="Subject")


### NOTE: I don't think we need histology/score/inflam scores from shitology_filename file!

#
# Histology: pathology and inflammation scores appended to bs_methyl
#
histology <- read_csv(file=histology_filename) %>%
  rename(Subject = "Record ID")

# note: Cumulative Score (Right)...24 and Cumulative Score (left) columns are the sum score of 
# variability in Firber, extent of central nucleation, necrosis and interstitial fiobrosis
save(histology, file=file.path(pkg_dir, "data", "histology.rda"))

# append pathology and inflammation score to bs_methyl
score <- histology %>% 
  dplyr::select(Subject, `Cumulative Score (Right)...24`, `Cumulative Score (left)`) %>%
  rename(R=`Cumulative Score (Right)...24`, L=`Cumulative Score (left)`) %>%
  gather(key=location, value=`Pathology Score`, -Subject) %>%
  #left_join(bogus_id, by="Subject") %>%
  dplyr::mutate(sample_id = paste0(Subject, location)) %>% 
  dplyr::select(-Subject, -location)

inflam <- histology %>%
  dplyr::select(Subject, `Inflammation Score (right)`, `Inflammation Score (left)`) %>%
  rename(R=`Inflammation Score (right)`, L=`Inflammation Score (left)`) %>%
  gather(key=location, value=`Inflammation Score`, -Subject) %>%
  dplyr::mutate(sample_id = paste0(Subject, location)) %>% 
  dplyr::select(-Subject, -location)

methyl_clinic <- bs_methyl %>%
  left_join(score, by="sample_id") %>%
  left_join(inflam, by="sample_id")

save(methyl_clinic, file=file.path(pkg_dir, "data", "methyl_clinic.rda"))  

# note that 32-0006 doesn't have Right biopsy

#
# summary data: this one has 13-0006R
#
mri <- read_xlsx(path=summary_filename, sheet=1) %>%
  dplyr::mutate(Site = if_else(`Site...1` == 1, "01", as.character(`Site...1`)),
                subj = if_else(subj < 10, paste0("000", subj), paste0("00", subj))) %>%
  dplyr::mutate(Subject = paste0(Site, "-", subj)) %>%
  dplyr::mutate(location = if_else(`SIDE (1=right)`==1, "R", "L")) %>%
  dplyr::mutate(sample_id = paste0(Subject, location)) %>%
  dplyr::rename(FAT_FRACTION = `_FAT_FRACTION...36`, `Cumulative Score` = `Cumulative Score...13`) %>%
  dplyr::select(Subject, sample_id, location, 
                `I. Variability in Fiber`,`II. Extent of Central Nucleation`,
                `III. Necrosis/Regeneration`, `IV. Interstitial Fiobrsis`,
                `Cumulative Score`,`Inflammation Score`,
                DIXON_VOL, FAT_MEAN, FAT_SD, 
                WAT_MEAN, WAT_SD, FAT_RATING, STIR_RATING, FAT_FRACTION, Fat_Infilt_Percent,
                `_UTE_MEAN_DIFF`) %>%
  dplyr::mutate(STIR_status =if_else(`STIR_RATING` == "0", "STIR-", "STIR+"))

ggplot(mri, aes(x=Fat_Infilt_Percent, y=FAT_FRACTION)) +
  geom_point(size=1, color="steelblue") +
  theme_bw()
ggsave(file.path(pkg_dir, "figures", "fat_fraction-vs-fat_infil_TA.pdf"), width=5, height=3)  

save(mri, file = file.path(pkg_dir, "data", "mri.rda"))

muscle_strength <- read_xlsx(path=baseline_filename, sheet=1) %>% 
  dplyr::select(`Record ID`, contains("Dorsiflexors")) %>%
  dplyr::filter(!is.na(`Record ID`)) %>%
  dplyr::mutate(`Right Foot Dorsiflexors` = as.numeric(`Right Foot Dorsiflexors`),
                `Left Foot Dorsiflexors`  = as.numeric(`Left Foot Dorsiflexors`)) %>%
  dplyr::filter(!is.na(`Right Foot Dorsiflexors`), !is.na(`Left Foot Dorsiflexors`))  

save(muscle_strength, file=file.path(pkg_dir, "data", "muscle_strength.rda"))               