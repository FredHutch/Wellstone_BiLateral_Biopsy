# This script checks 01-0034 and 01-0026 and H4 and M6 baskets

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(corrr)
library(GO.db)
require(goseq)
library(org.Hs.eg.db)

library(viridis)
library(hrbrthemes)

pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
fig_dir <- file.path(pkg_dir, "figures", "longitudinal_cohort")

load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "data", "cluster_df.rda"))
load(file.path(pkg_dir, "data", "all_baskets.rda"))
load(file.path(pkg_dir, "data", "candidates.rda"))
all(cluster_df$sample_name == colnames(sanitized.dds)) # checked
sanitized.dds$cluster <- as.character(cluster_df$new_cluster_name)
sanitized.dds$cluster[sanitized.dds$pheno_type == "Control"] <- "Control"
sanitized.dds$cluster <- factor(sanitized.dds$cluster,
                                levels=c("Control", "Mild", "Moderate", "IG-High", "High", "Muscle-Low"))

anno_ens88 <- as.data.frame(rowData(sanitized.dds)) %>%
  dplyr::select(gene_id, gene_name)


#
# color for stir_pal
#
library(wesanderson)
pal <- wes_palette("Darjeeling1", n=5)
stir_pal <- c("STIR-" = pal[5], "STIR+"= pal[1])
cntr_pal <- pal[2]
# make scatter plot of fat fraction vs. DUX4-H4, DUX4-M6, colored by STIR+
# make scatter plot for DUX4-M6, DUX4-H4 and other baskets and get the correlation        

#
# paired
#
col_data <- as.data.frame(colData(sanitized.dds))
paired_subject <- names(which(table(col_data$patient_id) == 2))

#
# fat/blood/muscle content = sum of log10(TPM +1) of marker genes in each content
#
markers_ens88 <- tibble(marker_type=c(rep("blood", 3), rep("fat", 3), rep("muscle", 3)),
                  gene_name=c("HBA1", "HBA2", "HBB",
                              "FASN", "LEP", "SCD",
                              "ACTA1", "TNNT3", "MYH1"))  %>%
  dplyr::mutate(marker_type = factor(marker_type)) %>%
  left_join(anno_ens88, by="gene_name")

markertype_tpm_longi <- sapply(levels(markers_ens88$marker_type), function(type) {
   id <- markers_ens88 %>% dplyr::filter(marker_type == type) %>%
      pull(gene_id)
   sub <- sanitized.dds[id]
   tpm_score <- colSums(log10(assays(sub)[["TPM"]]+1))
}) %>% as.data.frame() %>%
  rownames_to_column(var="sample_name") %>%
  dplyr::mutate(sample_name = if_else(sample_name=="32-0002b1", "32-0002b", sample_name))

#
# all basket score  and make comprehensive data.frame (longi_comprehensive_df)
#

# include DUX4-H4
all_baskets$`DUX4-H4` <- candidates %>% dplyr::filter(`basket-H4`) %>%
  rename(gene_id_v88 = gene_id) %>% dplyr::select(gene_id_v88, gene_name)
 
all_baskets_TPM <- map_dfr(names(all_baskets), function(name) {
  id <- all_baskets[[name]]$gene_id_v88    
  assays(sanitized.dds[id])[["TPM"]] %>% as.data.frame() %>%
      summarise(across(where(is.numeric), mean)) %>% t(.) %>% as.data.frame() %>%
      rownames_to_column(var="sample_name") %>%
      dplyr::mutate(sample_name = if_else(sample_name=="32-0002b1", "32-0002b", sample_name)) %>%
      add_column(basket = name) %>%
      dplyr::rename(TPM = V1)
}) %>% 
  dplyr::mutate(basket = factor(basket, levels=c("DUX4-H4", "DUX4-M", "DUX4-M6", "DUX4-M12", "Inflamm", "ECM", "Complement", "IG")))


spread_baskets <- all_baskets_TPM %>%
  spread(key=basket, value=TPM)

longi_comprehensive_df <- col_data %>%
  left_join(spread_baskets, by="sample_name") %>%
  dplyr::select(-lib_size, -raw_dir, -file_bam, -sizeFactor, -paired, -paired_end, -read_length, -frag_length) %>%
  dplyr::mutate(STIR_status = if_else(STIR_rating > 0, "STIR+", "STIR-")) %>%
  dplyr::mutate(STIR_status = factor(STIR_status, levels=c("STIR-", "STIR+"))) %>%
  left_join(markertype_tpm_longi, by="sample_name") %>%
  dplyr::arrange(patient_id)

writexl::write_xlsx(longi_comprehensive_df, path=file.path(pkg_dir, "stats", "longitudinal-cohort-comprehensive-df.xlsx"))  

#
# control baskets
#
load(file.path(pkg_dir, "data", "control_baskets.rda"))

basket_cntr <- control_baskets %>% group_by(basket) %>%
  summarise(mean=mean(TPM), sd = sd(TPM)) %>%
  dplyr::mutate(mean_plus_sd = mean+sd, mean_minus_sd = mean-sd)

#
# fat fraction vs DUX4-H4 and M6; color by STIR and shape by whether it is muscle-low
#
longi_comprehensive_df %>%
  dplyr::filter(pheno_type == "FSHD") %>%
  dplyr::filter(!is.na(fat_fraction), !is.na(`DUX4-M6`)) %>%
  dplyr::mutate(muscle_low = if_else(cluster == "Muscle-Low", TRUE, FALSE)) %>%
  ggplot(aes(x=fat_fraction, y=`DUX4-M6`)) +
    geom_point(aes(color=STIR_status, shape=muscle_low)) +
    theme_minimal() +
    scale_color_manual(values = stir_pal) +           
    scale_y_continuous(trans='log10') +  
    theme(legend.position="none") +
    labs(y="DUX4-M6 (TPM)") +
    geom_text_repel(aes(label=sample_id), size=1.8, show.legend=FALSE) +
    geom_hline(yintercept=basket_cntr %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean),
             color="gray50") +
    geom_hline(yintercept=basket_cntr %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean_plus_sd),
             color="gray50", linetype="dashed") 
ggsave(file=file.path(fig_dir, "fat-fraction-vs-DUX4-M6.pdf"), height=3, width = 4)


############### between visits ########################
#
# bar chart of DUX4-H4 and DUX4-M6 arranged by first year value; paired samples only
#
order_M6 <- longi_comprehensive_df %>%
  dplyr::filter(patient_id %in% paired_subject, visit=="I") %>%
  dplyr::arrange(`DUX4-M6`) %>% pull(patient_id)

longi_comprehensive_df %>%
  dplyr::filter(patient_id %in% paired_subject ) %>%
  dplyr::mutate(patient_id = factor(patient_id, levels=order_M6)) %>%
  dplyr::mutate(log10_DUX4M6 = log10(`DUX4-M6` +1)) %>%
  ggplot(aes(x=patient_id, y=log10_DUX4M6, fill=visit)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_minimal() +
    labs(y="log10(DUX4-M6 TPM +1)", x="") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(aes(label = Pathology.Score), size=2.5,  vjust = -0.2, hjust=0.5, position = position_dodge(.9))
ggsave(file=file.path(fig_dir, "DUX4-M6-bar-two-visit-log10.pdf"), width=6, height=3)

# Stephen mentioned about the super baskets, but their magnitude might not be compariable ...
longi_comprehensive_df %>%
  dplyr::filter(patient_id %in% paired_subject) %>%
  dplyr::select(sample_name, patient_id, visit, `DUX4-M6`,  Inflamm, ECM, Complement, IG ) %>%
  gather(key=basket, value=TPM, -sample_name, -patient_id, -visit) %>%
  dplyr::mutate(basket=factor(basket, levels=c("DUX4-M6", "Inflamm", "ECM", "Complement", "IG"))) %>%
  dplyr::mutate(log10_TPM = log10(TPM +1)) %>%
  ggplot(aes(x=basket, y=log10_TPM, fill=visit)) +
    geom_bar(width = 0.8, stat="identity", position=position_dodge(), alpha=0.7) +
    facet_wrap(~patient_id) +
    theme_minimal() +
    theme(legend.position=c(0.8, 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values=c("#999999", "#E69F00")) +
    labs(x="")
ggsave(file=file.path(fig_dir, "baskets-bar-two-visit-log10-baskets.pdf"))

longi_comprehensive_df %>%
  dplyr::filter(patient_id %in% paired_subject) %>%
  dplyr::select(sample_name, patient_id, visit, `DUX4-M6`,  Inflamm, ECM, Complement, IG ) %>%
  gather(key=basket, value=TPM, -sample_name, -patient_id, -visit) %>%
  dplyr::mutate(basket=factor(basket, levels=c("DUX4-M6", "Inflamm", "ECM", "Complement", "IG"))) %>%
  dplyr::mutate(log10_TPM = log10(TPM +1)) %>%
  ggplot(aes(x=basket, y=TPM, fill=visit)) +
    geom_bar(width=0.8, stat="identity", position=position_dodge(), alpha=0.7) +
    facet_wrap(~patient_id, scale="free_y") +
    theme_minimal() +
    theme(legend.position="none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values=c("#999999", "#E69F00")) +
    labs(x="")
ggsave(file=file.path(fig_dir, "baskets-bar-two-visit-TPM-baskets.pdf"))    

#################### blood/fat/muscle ###################
#
# boxplots of blood/fat/muscle content
#
longi_data <- longi_comprehensive_df %>% dplyr::select(sample_name, cluster, blood, fat, muscle) %>%
  gather(key=cell_type, value=TPM, -sample_name, -cluster) %>%
  dplyr::mutate(cohort = if_else(cluster=="Control", "Historical control", "Longitudinal")) %>%
  dplyr::mutate(muscle_low = if_else(cluster=="Muscle-Low", "muscle-low", "others")) %>%
  dplyr::select(-cluster)

ggplot(longi_data, aes(x=TPM)) +
  geom_density() +
  #geom_boxplot(width=0.7) +
  facet_wrap(~cell_type) +
  theme_bw() +
  geom_jitter(aes(x=TPM, y=0), size=1, height=0.01) +
  #geom_text_repel(data=mark_label, aes(y=0, label=sample_name), 
  #                size=2.5, fontface = "bold", color="red3", 
  #                min.segment.length = unit(0, 'lines'), 
  #                nudge_y = .08) +
  labs(x="sum of log10(TPM+1)", y="density")
ggsave(file.path(fig_dir, "longitudinal-biopsy-fat-blood-muscle-density.pdf"), height=2, width=5) 

longi_data %>% 
  ggplot(aes(y=TPM, x=cell_type)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0.3, size=1, aes(color=`muscle_low`), alpha=0.7) +
  theme_bw() +
  labs(y="sum of log10(TPM+1)", x="") +
  theme(legend.position="none") +
  scale_color_manual(values=c("#E69F00", "#999999"))

ggsave(file.path(fig_dir, "longitudinal-biopsy-fat-blood-muscle-boxplot.pdf"), height=3, width=2) 
 
# bilat + longitudinal; group by cohort
load(file.path(pkg_dir, "data", "comprehensive_df.rda"))

bilat_data <- comprehensive_df %>% dplyr::select(sample_id, blood, fat, muscle) %>%
    gather(key=cell_type, value=TPM, -sample_id) %>%
    dplyr::filter(!is.na(TPM)) %>%
    add_column(cohort = "Bilateral") %>%
    dplyr::mutate(muscle_low = if_else(sample_id %in% c("13-0007R", "13-0009R"), "muscle-low", "others")) %>%
    rename(sample_name=sample_id)

comb_data <- bilat_data  %>% bind_rows(longi_data) %>%
  dplyr::mutate(cohort = factor(cohort, level=c("Historical control", "Longitudinal", "Bilateral")))

comb_data %>% 
  ggplot(aes(y=TPM, x=cell_type)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0.3, size=1, aes(color=`muscle_low`), alpha=0.7) +
  theme_bw() +
  facet_wrap(~cohort, nrow=1) +
  labs(y="sum of log10(TPM+1)", x="") +
  theme(legend.position="none") +
  scale_color_manual(values=c("#E69F00", "#999999"))
ggsave(file.path(fig_dir, "longitudinal-bilat-fat-blood-muscle-boxplot.pdf"), height=3, width=4) 

# group by contnet
ggplot(comb_data, aes(x=cohort, y=TPM)) +
  geom_boxplot(width=0.5, aes(fill=cohort), alpha=0.5) +
  theme_bw() +
  facet_wrap(~cell_type, scale="free", nrow=3) +
  labs(x="", y="sum of log10(TPM+1)") +
  #annotate("text", x=2, y=2, label="Muscle-Low", size=2, color="red3") +
  coord_flip()  +
  theme(legend.position="none")
  #geom_text_repel(data=mark_label, aes(label=sample_name), 
  #                size=2, color="red3", 
  #                min.segment.length = unit(0, 'lines'), 
  #                nudge_y = .08, nudge_x=0.02) 

ggsave(file.path(fig_dir, "longitudinal-bilat-fat-blood-muscle-boxplot-by-content.pdf"), height=4, width=4) 

ggplot(comb_data, aes(x=cohort, y=TPM)) +
  geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0.3, size=0.6, aes(color=muscle_low), alpha=0.7) +
  theme_bw() +
  facet_wrap(~cell_type, scale="free", nrow=3) +
  labs(x="", y="sum of log10(TPM+1)") +
  #annotate("text", x=2, y=2, label="Muscle-Low", size=2, color="red3") +
  coord_flip()  +
  scale_color_manual(values=c("#E69F00", "#999999")) +
  theme(legend.position="none")
  #geom_text_repel(data=mark_label, aes(label=sample_name), 
  #                size=2, color="red3", 
  #                min.segment.length = unit(0, 'lines'), 
  #                nudge_y = .08, nudge_x=0.02) 

ggsave(file.path(fig_dir, "longitudinal-bilat-fat-blood-muscle-boxplot-point-by-content.pdf"), height=4, width=4) 


#
# bloody samples in bilat 
#

# boxplot of baskets for muscle-low, bloody (blood content > 12), and others
comprehensive_df %>% 
  dplyr::mutate(group=case_when(sample_id %in% c("13-0007R", "13-0009R") ~ "muscle-low",
                                sample_id %in% c("32-0020R", "32-0029R") ~ "bloody",
                                TRUE ~ "others")) %>%
  dplyr::select(group, `DUX4-M6`, Inflamm, ECM, Complement, IG) %>%
  dplyr::filter(!is.na(`DUX4-M6`)) %>%
  gather(key=basket, value=TPM, -group) %>%
  ggplot(aes(x=group, y=TPM)) +
    geom_boxplot(width=0.5, outlier.shape=NA, aes(fill=group), show.legend=FALSE) +
    geom_jitter(width=0.2, color="grey50", alpha=0.7) +
    facet_wrap(~basket, scale="free_y") +
    theme_bw() +
    scale_y_continuous(trans='log10') 

ggsave(file=file.path(pkg_dir, "figures", "baskets-group-by-cell-characteristics.pdf"))