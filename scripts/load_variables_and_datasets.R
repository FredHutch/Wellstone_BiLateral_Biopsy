# load longitudinal.dds and bilateral.dds and associated metdata

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(kableExtra)
suppressPackageStartupMessages(library(writexl))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(purrr))
library(wesanderson)
library(latex2exp)

pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
load(file.path(pkg_dir, "data", "bilat_dds.rda"))
load(file.path(pkg_dir, "data", "longitudinal_dds.rda"))
load(file.path(pkg_dir, "data", "bilat_MLpredict.rda"))
load(file.path(pkg_dir, "data", "comprehensive_df.rda"))
source(file.path(pkg_dir, "scripts", "viz_tools.R"))

# tidy annotation from two datasets
anno_gencode35 <- as.data.frame(rowData(bilat_dds)) %>%
  rownames_to_column(var="gencode35_id") %>% # BiLat study using Gencode 36
  dplyr::mutate(ens_id=str_replace(gencode35_id, "\\..*", "")) %>%
  dplyr::distinct(gene_name, .keep_all = TRUE)

anno_ens88 <- as.data.frame(rowData(longitudinal_dds)) %>%
  rownames_to_column(var="ens88_id") %>% # longitudinal study ens v88
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", "")) %>%
  dplyr::distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::select(ens88_id, ens_id, gene_name)
# insert class to bilat_dds
col_data <- colData(bilat_dds) %>% as.data.frame() %>%
  left_join(bilat_MLpredict, by="sample_id")
bilat_dds$class <- col_data$class

longitudinal_dds$STIR_status <- if_else(longitudinal_dds$STIR_rating > 0, "STIR+", "STIR-")
