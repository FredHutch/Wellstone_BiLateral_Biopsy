library(DESeq2)
library(tidyverse)
library(ggrepel)
library(corrr)
library(latex2exp)

pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
fig_dir <- file.path(pkg_dir, "manuscript", "figures")
load(file.path(pkg_dir, "data", "dds.rda"))
load(file.path(pkg_dir, "data", "mri.rda"))
load(file.path(pkg_dir, "data", "blood_fat_muscle_content.rda"))
load(file.path(pkg_dir, "data", "comprehensive_df.rda"))
load(file.path(pkg_dir, "data", "control_baskets.rda"))
load(file.path(pkg_dir, "data", "control_baskets_logsum.rda"))

anno_gencode36 <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id") %>% # BiLat study using Gencode 36
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", ""))

#
# longitudinal study
#
load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "data", "cluster_df.rda"))
sanitized.dds$cluster <- as.character(cluster_df$new_cluster_name)
sanitized.dds$cluster[sanitized.dds$pheno_type == "Control"] <- "Control"
sanitized.dds$cluster <- factor(sanitized.dds$cluster,
                                levels=c("Control", "Mild", "Moderate", "IG-High", "High", "Muscle-Low"))

anno_ens88 <- as.data.frame(rowData(sanitized.dds)) %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::mutate(ens_id=str_replace(gene_id, "\\..*", ""))

markers_ens88 <- tibble(marker_type=c(rep("blood", 3), rep("fat", 3),
                                      rep("muscle", 3)),
                  gene_name=c("HBA1", "HBA2", "HBB",
                              "FASN", "LEP", "SCD",
                              "ACTA1", "TNNT3", "MYH1"))  %>%
  dplyr::mutate(marker_type = factor(marker_type)) %>%
  left_join(anno_ens88, by="gene_name")
 
  
#
# color for stir_pal
#
library(wesanderson)
pal <- wes_palette("Darjeeling1", n=5)
stir_pal <- c("STIR-" = pal[5], "STIR+"= pal[1])
cntr_pal <- pal[2]

#
# control basket
#
basket_cntr <- control_baskets %>% group_by(basket) %>%
  summarise(mean=mean(TPM), sd = sd(TPM)) %>%
  dplyr::mutate(mean_plus_sd = mean+sd, mean_minus_sd = mean-sd)

basket_cntr_logsum <- control_baskets_logsum %>% 
  group_by(basket) %>%
  summarise(mean=mean(log10TPM), sd = sd(log10TPM)) %>%
  dplyr::mutate(mean_plus_sd = mean+sd, mean_minus_sd = mean-sd)


#
# Fig 1:; cell type content and PCA dottype by STIR status
#

# (a) cell type content density plot
data <- blood_fat_muscle_content  %>% 
  gather(key=cell_type, value=TPM, -sample_name, -RNA_sample_id)  %>%
  dplyr::mutate(cell_type = factor(cell_type))
mark_label <- data %>% 
  dplyr::filter(sample_name %in% c("13-0009R", "13-0007R", "32-0020R", "32-0029R")) %>%
  dplyr::filter(!(sample_name %in% c("13-0009R", "13-0007R") & cell_type == "blood")) %>%
  dplyr::filter(!(sample_name %in% c("32-0020R", "32-0029R") & 
                                       cell_type %in% c("fat", "muscle"))) 
# qunatile for each of the cell type content
df_quantile <- data %>% group_by(cell_type) %>% 
  summarise(q0 = quantile(TPM, probs=seq(0, 1, 0.01))["0%"], 
            q3 = quantile(TPM, probs=seq(0, 1, 0.01))["3%"],
            q97 = quantile(TPM, probs=seq(0, 1, 0.01))["97%"],
            q100 = quantile(TPM, probs=seq(0, 1, 0.01))["100%"]) %>% as.data.frame()
df_rect <- data.frame(cell_type = df_quantile$cell_type,
                      x1 = c(0, 0, df_quantile[3, 3]),
                      x2 = c(df_quantile[1:2, 4], Inf))
ggplot(data, aes(x=TPM)) +
  geom_density() +
  facet_wrap(~cell_type, scales="free_x") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        strip.background =element_rect(fill="transparent")) +
  geom_jitter(aes(x=TPM, y=0), size=1, height=0.01) +
  geom_text_repel(data=mark_label, aes(y=0, label=sample_name), 
                  size=2.5, color="gray20", 
                  min.segment.length = unit(0, 'lines'), 
                  nudge_y = .1) +
  geom_rect(data = df_rect, fill="grey50", alpha=0.2,
            aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf), inherit.aes = FALSE) +                  
  labs(x=latex2exp::TeX("RNA-seq cell type content ($\\log_{10} TPM$$)"), y="density") 

ggsave(file.path(fig_dir, "biopsy-fat-blood-muscle-density.pdf"), 
       height=1.8, width=4.8)

# (b) PCA colored by STIR status
library(ggrepel)
load(file.path(pkg_dir, "data", "rlog.rda"))
rlog <- rlog[, !colnames(rlog) == "32-0028L_S4"]
col_data <- as.data.frame(colData(rlog)) %>%
  dplyr::mutate(sample_id = str_replace(sample_name, "_.*", "") )%>%
  dplyr::mutate(sample_id = str_replace(sample_id, "b$|b1$|-1$", "")) %>%
  left_join(dplyr::select(mri, sample_id, STIR_status), by="sample_id") %>%
  dplyr::mutate(STIR_status = factor(STIR_status, levels=c("STIR-", "STIR+")))
rlog$STIR_status <- col_data$STIR_status
rlog$sample_id <- col_data$sample_id

data <- plotPCA(rlog, intgroup=c("location", "STIR_status", "sample_id"), returnData=TRUE)

percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, shape=STIR_status)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title="PCA") +
  geom_text_repel(aes(label=sample_id), size=1.8, 
                  show.legend=FALSE) +
  theme_minimal() +
  #guides(fill=guide_legend(ncol=2)) +
  #scale_color_manual(values = stir_pal) +           
  theme(#legend.position="bottom", 
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.1, 0.1),
        panel.grid.minor.x = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=10))
ggsave(file.path(fig_dir, "biopsy-PCA.pdf"), width=3.5, height=3)

# (c) longitudinal and bilat combined

# longitudinal data
markertype_tpm_longi <- sapply(levels(markers_ens88$marker_type), function(type) {
   id <- markers_ens88 %>% dplyr::filter(marker_type == type) %>%
      pull(gene_id)
   sub <- sanitized.dds[id]
   tpm_score <- colSums(log10(assays(sub)[["TPM"]]+1))
}) %>% as.data.frame() %>%
  rownames_to_column(var="sample_name") %>%
  add_column(pheno_type = sanitized.dds$pheno_type,
             classes = sanitized.dds$cluster)

tmp_longi <- markertype_tpm_longi %>% 
  dplyr::mutate(study = if_else(pheno_type=="FSHD", "Longitudinal", "Historical control")) %>%
  dplyr::mutate(group = if_else(classes == "Muscle-Low", "Muscle-Low", "others")) %>%
  dplyr::select(-pheno_type, -classes)

tmp_rbind <- blood_fat_muscle_content %>%
  dplyr::select(-RNA_sample_id) %>%
  dplyr::relocate(sample_name, .before=blood) %>%
  dplyr::mutate(study = "Bilat") %>%
  dplyr::mutate(group = if_else(sample_name %in% c("13-0007R", "13-0009R"), "Muscle-Low", "others")) %>%
  bind_rows(tmp_longi) %>%
  gather(key=type, value=TPM, -sample_name, -study, -group) %>%
  dplyr::filter(type=="muscle") %>%
  dplyr::mutate(study = factor(study, levels=c("Historical control", "Longitudinal", "Bilat")))

mark_label <- tmp_rbind %>% 
  dplyr::filter(group == "Muscle-Low" & study == "Bilat")

ggplot(tmp_rbind, aes(x=study, y=TPM)) +
  geom_boxplot(width=0.5) +
  theme_minimal() +
  labs(x="", y=latex2exp::TeX("RNA-seq muscle content ($\\log_{10}TPM$)")) +
  annotate("text", x=2.5, y=3.75, hjust=0.5, vjust=0.5,
           label="Muscle-Low", size=2, color="grey25") +
  #facet_wrap(~type, nrow=3, scale="free_y") +
  coord_flip() +
  geom_text_repel(data=mark_label, aes(label=sample_name), 
                  size=2, color="grey25", 
                  min.segment.length = unit(0, 'lines'), 
                  nudge_y = .08, nudge_x=0.02) +
  labs(title="Bilat and longitudinal studies", 
       y=latex2exp::TeX("RNA-seq muscle cell content $\\log_{10} TPM") )+
  theme(plot.title=element_text(hjust=0.5, size=10), 
        axis.title.x=element_text(size=10))
  
ggsave(file.path(fig_dir, "longitudinal-bilat-muscle-boxplot.pdf"),
        width=4, height=1.5)

#
# Fig 2: (a) and (b)
#

# prepare data for (a) and (b)
fat_df  <- blood_fat_muscle_content %>%
  dplyr::rename(sample_id = sample_name) %>%
  left_join(dplyr::select(mri, sample_id, Fat_Infilt_Percent, FAT_FRACTION, STIR_status), by="sample_id")

fat_df %>% 
  dplyr::select(Fat_Infilt_Percent, FAT_FRACTION, fat) %>%
  corrr::correlate() 
  
fat_df %>% dplyr::filter(!sample_id %in% c("13-0007R",  "13-0009R")) %>%
  ggplot(aes(x=Fat_Infilt_Percent, y=fat)) +
    geom_smooth(method="lm", se=FALSE, color="grey75", alpha=0.5) +
    geom_point(size=1.2) +
    theme_bw() +
    geom_text_repel(aes(label=sample_id), size=2, show.legend=FALSE) +
    geom_text(aes(x=Inf, y=-Inf, label="cor = 0.53"), size=2.5, vjust=-1, hjust=1.2) +
    annotate("text", x=Inf, y=-Inf, label="cor = 0.53", size=2.5, vjust=-1, hjust=1.2) +
  labs(x="whole muscle fat percent", y="RNA-seq fat content")
ggsave(file.path(pkg_dir, "figures", "Fat-Infilt-vs-fat-content.pdf"),
       height = 2.5, width = 3.5)

# another try: figure to include fat fraction vs fat content
tmp <- comprehensive_df %>% 
  dplyr::filter(!is.na(fat), !is.na(Fat_Infilt_Percent)) %>%
  dplyr::filter(!sample_id %in% c("13-0007R",  "13-0009R")) %>%
  dplyr::select(Fat_Infilt_Percent, FAT_FRACTION, fat) %>%
  rename(`whole muscle fat percent`=Fat_Infilt_Percent, 
         `regional fat fraction`=FAT_FRACTION)

corr_df <- tmp %>% corrr::correlate()  %>%
  dplyr::select(term, fat) %>%
  dplyr::filter(term != "fat") %>%
  dplyr::mutate(term = factor(term, levels=c("whole muscle fat percent", "regional fat fraction")))

# Figure 2(a)
ggplot(tmp, aes(x=`whole muscle fat percent`, y=`regional fat fraction`)) +
  geom_point() +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  geom_hline(yintercept=0.1, color="grey50", alpha=0.5, linetype="dashed") 

ggsave(file.path(fig_dir, "fat-infilt-percent-vs-fat-faction.pdf"),
       height=2, width = 3)
# (b)
tmp %>%
  gather(key=term, value=percentage, -fat) %>%
  dplyr::mutate(term = factor(term)) %>% 
  ggplot(aes(x=percentage, y=fat)) + 
  geom_smooth(method="lm", se=FALSE, color="grey75", alpha=0.5) +
    geom_point(size=1.2) +
    facet_wrap(~term, nrow=1) +
    theme_bw() +
    geom_text(data=corr_df, aes(x=Inf, y=-Inf, label=paste0("Pearson=", format(fat, digit=2))), 
              size=2.5, vjust=-1, hjust=1.2) +
    labs(x="", y="RNA-seq fat content")

ggsave(file.path(fig_dir, "T1-vs-fat-marker-content.pdf"),
       height = 2, width = 4.5)
#
# Fig 3
#
basket_cntr_logsum <- control_baskets_logsum %>% group_by(basket) %>%
  summarise(mean=mean(log10TPM), sd = sd(log10TPM), n=length(basket)) %>%
  dplyr::mutate(lower_ci = mean - 1.96*(sd/n), upper_ci=mean+1.96*(sd/n),
                mean_plus_sd = mean+sd, mean_minus_sd = mean-sd)

# (a) whole muscel fat percent vs STIR-/+
pval = summary(lm(comprehensive_df$Fat_Infilt_Percent ~ comprehensive_df$STIR_status))$coefficients[2, "Pr(>|t|)"] 
outliers = comprehensive_df %>% dplyr::filter(STIR_status == "STIR-", Fat_Infilt_Percent > 0.7)
annot = data.frame(pval=pval)

ggplot(comprehensive_df, aes(x=STIR_status, y=Fat_Infilt_Percent)) +
  geom_boxplot(width=0.6) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(), 
        axis.title.y=element_text(size=9),
        axis.text.x = element_text(size=8)) +
  labs(x="", y="whole muscle fat percent") +
  annotate("text", color="grey25", x=1.5, y=0.94, label=TeX("p=2e-8"), size=1.8) +
  geom_segment(aes(x = 1, y = 0.91, xend = 2, yend = 0.91), size=0.2) +
  #geom_text(data=annot, size=2.5, fontface="italic", aes(x=1.5, y=0.9, 
  #          label=TeX("$p=2e-8$")) )+
  geom_text_repel(data=outliers, aes(label=sample_id), color="gray20",
                  nudge_x=0.2, size=1.5, show.legend=FALSE) 
                    
ggsave(file=file.path(fig_dir, "Fat-infilt-vs-STIR-status.pdf"), 
       width=1.5, height=2)

# (b) DUX4-M6 baskets (logSum score) vs fat infiltrate percent
gg <- comprehensive_df %>% dplyr::filter(!is.na(`DUX4-M6`), !is.na(Fat_Infilt_Percent)) %>%
  dplyr::filter(!sample_id %in% c("13-0007R",  "13-0009R")) %>%
  ggplot(aes(x=Fat_Infilt_Percent, y=`DUX4-M6-logSum`)) +
    geom_point(aes(color=STIR_status), alpha=0.8) +
    theme_minimal() +
    scale_color_manual(values = stir_pal) +           
    labs(x="whole muscle fat percent", y=TeX("DUX4 M6 score: $\\sum \\log_{10}(TPM_i+1)$")) +
    theme(legend.title = element_blank(),
          legend.position = "none",
          legend.box.background = element_rect(colour = "grey50"))
gg +  
  #geom_smooth(method="loess", se=FALSE) +
  geom_hline(yintercept=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean),
             color="gray75", linetype="dashed") +
  annotate("text", label="historical control mean", x=0.85, 
             y=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean), 
             color="gray50", hjust=1, vjust=-0.5)             
ggsave(file=file.path(fig_dir, "fat-infilt-percent-vs-DUX4-M6-logsum.pdf"),
       width=5, height=3.5)

gg +  
  geom_rect(aes(xmin=-Inf, xmax=Inf, 
                ymin=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(lower_ci), 
                ymax=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(upper_ci)), 
                fill="grey92", alpha=0.08) +
  annotate("text", label="historical control 95% CI", x=0.85, 
             y=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean), 
             color="gray50", hjust=1, vjust=0.5) +             
  scale_y_continuous(trans='log10')   
ggsave(file=file.path(fig_dir, "fat-infilt-percent-vs-DUX4-M6-logsum-scale-y.pdf"),
       width=5, height=3.5)    

#
# supplementary figure 3a-b: 
# logSum for all baskets, exclude M6 bud include DUX4-M (rename to M34)
#
lvls <- c( "DUX4-M12", "DUX4-M34", "Inflamm",
             "ECM", "Complement", "IG")

cntr <- basket_cntr_logsum %>% dplyr::filter(basket!="DUX4-M6") %>%
  dplyr::mutate(basket = as.character(basket)) %>%
  dplyr::mutate(basket = if_else(basket=="DUX4-M", "DUX4-M34", basket)) %>% 
  dplyr::mutate(basket = factor(basket, levels=lvls))

gg <- comprehensive_df %>% dplyr::filter(!is.na(`DUX4-M6`), !is.na(Fat_Infilt_Percent)) %>%
  dplyr::filter(!sample_id %in% c("13-0007R",  "13-0009R")) %>%
  dplyr::select(Fat_Infilt_Percent, `DUX4-M-logSum`, `DUX4-M12-logSum`, `ECM-logSum`,
               `Inflamm-logSum`, `Complement-logSum`, `IG-logSum`, sample_id, STIR_status) %>%
  gather(key=basket, value=logSum, -Fat_Infilt_Percent, -sample_id, -STIR_status) %>%
  dplyr::mutate(basket=str_replace(basket, "-logSum", "")) %>%
  dplyr::mutate(basket = if_else(basket == "DUX4-M", "DUX4-M34",basket)) %>% 
  dplyr::mutate(basket=factor(basket, levels=lvls)) %>%
  ggplot(aes(x=Fat_Infilt_Percent, y=logSum)) +
    geom_point(aes(color=STIR_status), alpha=0.7) +
    theme_bw() +
    facet_wrap(~basket, nrow=2, scale="free_y") +
    #geom_hline(aes(yintercept = mean_plus_sd), data=a, linetype="dashed", color="gray75") +
    #geom_rect(data=cntr, aes(xmin=-Inf, xmax=Inf, ymin=lower_ci, ymax=upper_ci),
    #          fill="gray92", alpha=0.08, inherit.aes=FALSE) +
    scale_color_manual(values = stir_pal) +           
    labs(title="basket scores vs. whole muscle fat percent",
    x="whole muscle fat percent", y=TeX("basket: $\\sum \\log_{10}(TPM_i+1)$")) +
    theme(legend.position="none", panel.grid.minor = element_blank(),
          plot.title=element_text(hjust=0.5),
          strip.background =element_rect(fill="transparent") )

gg + 
  geom_hline(aes(yintercept = mean), data = cntr, color="gray75", linetype="dashed") +
  geom_smooth(method="loess", se=FALSE, linewidth=0.7)
ggsave(file=file.path(fig_dir, "fat-infilt-percent-vs-baskets-logsum.pdf"),
       height=4, width=7)

gg + 
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=lower_ci, ymax=upper_ci),
            data=cntr, fill="grey25", alpha=0.15, inherit.aes=FALSE) +
  scale_y_continuous(trans='log10')            

ggsave(file=file.path(fig_dir, "fat-infilt-percent-vs-baskets-logsum-scaled-y.pdf"),
       height=4, width=7)         

# suppl figure 3c? where is it?

#
# Figure 4: DUX4-M6 vs other clinical variables
#       

# (c) DUX4-M6 vs. cumulative score

gg <- comprehensive_df %>% dplyr::filter(!is.na(`Cumulative Score`), !is.na(`DUX4-M6-logSum`)) %>%
  dplyr::filter(!sample_id %in% c("13-0007R",  "13-0009R")) %>%
  rename(`histopathologic score`=`Cumulative Score`) %>%
  ggplot(aes(x=`histopathologic score`, y=`DUX4-M6-logSum`)) +
    geom_point(aes(color=STIR_status)) +
    geom_smooth(method="lm", se=FALSE, color="grey75", alpha=0.5, linewidth=0.7) +
    theme_minimal() +
    scale_color_manual(values = stir_pal) +           
    labs(x="", y="DUX4-M6 score",
         title="Histopathologic score vs. DUX4 score") +
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=12))      

gg +   
  geom_hline(yintercept=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean),
             color="gray75", linetype="dashed") +
  annotate("text", label="historical control mean", x=12.5, 
           y=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean), 
           color="gray50", hjust=1, vjust=-0.5) 
ggsave(file=file.path(fig_dir, "histo-vs-DUX4-M6-logsum.pdf"),
       width=4.5, height=3) 

gg +     
  geom_rect(aes(xmin=-Inf, xmax=Inf, 
                ymin=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(lower_ci), 
                ymax=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(upper_ci)), 
                fill="grey92", alpha=0.08) + 
  annotate("text", label="historical control 95% CI", x=12.5, 
           y=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean), 
          color="gray50", hjust=1, vjust=0.5) +
  scale_y_continuous(trans='log10', breaks=c(0.01, 0.1, 1, 3))
ggsave(file=file.path(fig_dir, "histo-vs-DUX4-M6-logsum-scale-y.pdf"),
       width=4.5, height=3) 


# (b) DUX4 vs. CSS
gg <- comprehensive_df %>% dplyr::filter(!is.na(`CSS`), !is.na(`DUX4-M6-logSum`)) %>%
  dplyr::filter(!sample_id %in% c("13-0007R",  "13-0009R")) %>%
  ggplot(aes(x=`CSS`, y=`DUX4-M6-logSum`)) +
    geom_point(aes(color=STIR_status)) +
    geom_smooth(method="lm", se=FALSE, color="grey75", alpha=0.5, linewidth=0.7) +
    theme_minimal() +
    theme(legend.position="none", plot.title=element_text(hjust=0.5, size=12)) +
    scale_color_manual(values = stir_pal) +           
    labs(x="", y="DUX4-M6 score", 
         title="CSS vs. DUX4 score")
gg +   
  geom_hline(yintercept=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean),
             color="gray75", linetype="dashed") +
  annotate("text", label="historical control mean", x=8, 
           y=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean), 
           color="gray50", hjust=1, vjust=-0.5) 
ggsave(file=file.path(fig_dir, "CSS-vs-DUX4-M6-logsum.pdf"),
       width=4.5, height=3)              

gg +     
  geom_rect(aes(xmin=-Inf, xmax=Inf, 
                ymin=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(lower_ci), 
                ymax=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(upper_ci)), 
                fill="grey92", alpha=0.08) + 
  annotate("text", label="historical control 95% CI", x=8, 
           y=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean), 
          color="gray50", hjust=1, vjust=0.5) +
  scale_y_continuous(trans='log10', breaks=c(0.01, 0.1, 1, 3))
ggsave(file=file.path(fig_dir, "CSS-vs-DUX4-M6-logsum-scale-y.pdf"),
       width=4.5, height=3)            

# (a) DUX4 vs. foot dorsiflexors
gg <- comprehensive_df %>% dplyr::filter(!is.na(`Foot Dorsiflexors`), !is.na(`DUX4-M6-logSum`)) %>%
  dplyr::filter(!sample_id %in% c("13-0007R",  "13-0009R")) %>%
  ggplot(aes(x=`Foot Dorsiflexors`, y=`DUX4-M6-logSum`)) +
    geom_point(aes(color=STIR_status)) +
    geom_smooth(method="lm", se=FALSE, color="grey75", alpha=0.5, linewidth=0.57) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust=0.5, size=12)) +
    scale_color_manual(values = stir_pal) +           
    labs(x="", y="DUX4-M6 score", 
         title="Foot dorsiflexors vs. DUX4 score")

gg +   
  geom_hline(yintercept=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean),
             color="gray75", linetype="dashed") +
  annotate("text", label="historical control mean", x=43, 
           y=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean), 
           color="gray50", hjust=1, vjust=1) 
ggsave(file=file.path(fig_dir, "Foot-Dorsiflexors-vs-DUX4-M6-logsum.pdf"),
       width=4.5, height=3)              

gg +     
  geom_rect(aes(xmin=-Inf, xmax=Inf, 
                ymin=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(lower_ci), 
                ymax=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(upper_ci)), 
                fill="grey92", alpha=0.08) + 
  annotate("text", label="historical control 95% CI", x=43, 
           y=basket_cntr_logsum %>% dplyr::filter(basket=="DUX4-M6") %>% pull(mean), 
          color="gray50", hjust=1, vjust=0.5) +
  scale_y_continuous(trans='log10', breaks=c(0.01, 0.1, 1, 3))
ggsave(file=file.path(fig_dir, "Foot-Dorsiflexors-vs-DUX4-M6-logsum-scale-y.pdf"),
       width=4.5, height=3)  

# 

# comprehensive correlation: supplemental figure 3d
comprehensive_df %>% 
  dplyr::filter(!is.na(`Foot Dorsiflexors`), !is.na(STIR_status), !is.na(`DUX4-M6`)) %>%
  dplyr::select(`I. Variability in Fiber`, `II. Extent of Central Nucleation`, 
                `III. Necrosis/Regeneration`, `IV. Interstitial Fiobrsis`, 
                `Cumulative Score`, Fat_Infilt_Percent, FAT_FRACTION, STIR_RATING, CSS,
                `Foot Dorsiflexors`, 
                `DUX4-M6-logSum`, `Inflamm-logSum`, `ECM-logSum`, `Complement-logSum`, `IG-logSum`) %>%
  GGally::ggcorr(label = TRUE, label_round=1, label_size = 3, size=2, hjust=0.9, legend.position="none", layout.exp=3)

ggsave(file=file.path(fig_dir, "comprehensive-correlation-basket-logsum.pdf"), width=6, height=4)




#
# figure 5: baskets vs STIR+/- and between basket correlation with scatter plots
# those were in 06-figures-for-grant.R
#


# (a) baskets vs STIR+/-
add_control <- control_baskets_logsum %>%
  dplyr::filter(basket %in% c("DUX4-M6", "ECM", "Inflamm", "Complement", "IG")) %>%
  tibble::add_column(STIR_status="Control", .after="sample_id") %>%
  dplyr::relocate(basket, .before=log10TPM)

gg <- comprehensive_df %>% 
  dplyr::filter(!sample_id %in% c("13-0007R",  "13-0009R")) %>%
  dplyr::select(sample_id, `ECM-logSum`, `Inflamm-logSum`, `Complement-logSum`, 
                `IG-logSum`, `DUX4-M6-logSum`, STIR_status) %>%
  dplyr::filter(!is.na(STIR_status)) %>%
  gather(key=basket, value=log10TPM, -sample_id, -STIR_status) %>%
  dplyr::mutate(basket = str_replace(basket, "-logSum", "")) %>%
  add_row(add_control) %>%
  dplyr::mutate(STIR_status = factor(STIR_status, levels=c("Control", "STIR-", "STIR+")),
                basket = factor(basket, levels=c("DUX4-M6", "Inflamm", "ECM", "Complement", "IG"))) %>%
  ggplot(aes(x=STIR_status, y=log10TPM)) +
    geom_boxplot(width=0.7, outlier.shape=NA, fill="grey75", alpha=0.5) +
    geom_jitter(width = 0.3, size=0.5, alpha=0.5) +
    facet_wrap(~ basket, scales="free_y", nrow=1) +
    theme_bw() +
    labs(x="", y="basket scores") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8, hjust=1),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill="transparent"), 
          strip.text.x = element_text(size = 8))
gg
ggsave(file=file.path(fig_dir, "STIR-vs-all-baskets.pdf"), width=5.2, height=2.5 )

gg + scale_y_continuous(trans='log10')
ggsave(file=file.path(fig_dir, "STIR-vs-all-baskets-scale-y.pdf"), width=5.4, height=2.5 )


# (b) correlation between baskets with scatter plot
library(GGally)
gg <- comprehensive_df %>% 
  dplyr::filter(!sample_id %in% c("13-0007R",  "13-0009R")) %>%
  dplyr::select(`DUX4-M6-logSum`, `Inflamm-logSum`, `ECM-logSum`, `Complement-logSum`, 
                `IG-logSum`) %>%
  dplyr::rename_with(~str_replace(.x, "-logSum", "")) %>%              
  ggpairs(aes(ize=0.7, alpha=0.7)) + 
       theme_bw() +
     theme(panel.grid.minor = element_blank(), 
           strip.text.y = element_text(size = 8),
           strip.text.x = element_text(size = 8),
           strip.background = element_rect(fill="transparent") )
gg     
ggsave(file=file.path(fig_dir, "basekt-m6-other-baskets.pdf"), width=4.5, height=4.5)  

gg  + scale_y_continuous(trans='log10')
ggsave(file=file.path(fig_dir, "basekt-m6-other-baskets-scale-y.pdf"), width=4.5, height=4.5)  


#
# Figure 6: gitbook.2 -> 07-bilateral-comparison.Rmd
#


#
# Figure 7: 03D-BSS.R; need to migrate to 08-BSS.Rmd
#