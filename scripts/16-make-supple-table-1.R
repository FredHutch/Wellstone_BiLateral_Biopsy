# make suppl table: from Yao2014-robust-66-DUX4-targets.xlsx
library(tidyverse)
library(readxl)
library(writexl)

pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
load(file.path(pkg_dir, "data", "comprehensive_df.rda"))
load(file.path(pkg_dir, "data", "bs_methyl.rda"))
rin <- readxl::read_xlsx(path=file.path(pkg_dir, "manuscript", 
                                        "comprehensive_df-RINN.xlsx")) %>%
  dplyr::select(sample_id, RINN)

# comprehensive-data-frame.xlsx to include comprehensive_df and BSS values
file_name <- file.path(pkg_dir, "manuscript", "suppl-tables", 
                       "suppl-table-3-comprehensive-data-frame.xlsx")

out <- comprehensive_df %>%
  dplyr::select(-`DUX4-M`, -`DUX4-M6`, -`DUX4-M12`, -`Inflamm`,
               -ECM, `Complement`, `DUX4+ (M6)`) %>%
  dplyr::left_join(rin, by="sample_id") %>%
  dplyr::rename(RIN=RINN) %>%
  dplyr::relocate(RIN, .after=location)

writexl::write_xlsx(list(comprehensive_df = out, bisulfite_seq=bs_methyl),
                    path=file_name)
  
