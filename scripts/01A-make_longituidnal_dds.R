# convert sanity.dds to longitudinal_dds
# main added feature: (1) add class to metadata, and (2) perform DESeq based on cluster

library(DESeq2)

pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "data", "cluster_df.rda"))

all(cluster_df$sample_name == colnames(sanitized.dds)) # checked

sanitized.dds$cluster <- as.character(cluster_df$new_cluster_name)
sanitized.dds$cluster[sanitized.dds$pheno_type == "Control"] <- "Control"
sanitized.dds$cluster <- factor(sanitized.dds$cluster,
                                levels=c("Control", "Mild", "Moderate", "IG-High", "High", "Muscle-Low"))

#
# perform DESeq based on cluster; note that sizeFactor and dispersion have been pre-defined
#
longitudinal_dds <- sanitized.dds
design(longitudinal_dds) <- ~ cluster
# sanity check
colnames(longitudinal_dds)==longitudinal_dds$sample_id
# conform colnames and sample_id
colnames(longitudinal_dds) <- longitudinal_dds$sample_id

save(longitudinal_dds, file=file.path(pkg_dir, "data", "longitudinal_dds.rda"))