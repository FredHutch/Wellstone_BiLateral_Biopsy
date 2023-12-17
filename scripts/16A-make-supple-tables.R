# make suppl table: 
# 1. DUX4 baskets statistics
# 2. all other basket genes

library(tidyverse)
library(readxl)
library(writexl)

pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
load(file.path(pkg_dir, "data", "comprehensive_df.rda"))
#load(file.path(pkg_dir, "data", "bs_methyl.rda"))

#####
## Supple Table 2
## Basket genes
#####


#####
## Suppl. Table 1
## comprehensive_df.rda including all clinical data, pathological scores,  basket scores
#####

# it is made by 09b-make-comprehensive-df.R ->  already
file_name = file.path(pkg_dir, "stats", "comprehensive_metadata.xlsx")

###
## Suppl. Table 2
## DUX4 baskets and others
###
load(file.path(pkg_dir, "data", "DUX4_baskets.rda"))
load(file.path(pkg_dir, "data", "all_baskets.rda"))
writexl::write_xlsx(list(`DUX4 basekets statistics` = DUX4_baskets,
                         `DUX4-M` = all_baskets[["DUX4-M"]],
                         `DUX4-M6` = all_baskets[["DUX4-M6"]],
                         `DUX4-M12` = all_baskets[["DUX4-M12"]],
                         `ECM` = all_baskets[["ECM"]],
                         `Inflamm` = all_baskets[["Inflamm"]],
                         `Complement` = all_baskets[["Complement"]],
                         `IG` = all_baskets[["IG"]]),
                    path = file.path(pkg_dir, "manuscript", 
                                     "suppl-tables/suppl-table-2-DUX4-and-other-baskets.xlsx"))
  
