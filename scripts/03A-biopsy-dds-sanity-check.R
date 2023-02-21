# 03A-biopsy-quality-control.R
#
# this quality control following this gitpage: https://fredhutch.github.io/RWellstone_FSHD_muscle_biopsy/biopsy-qc.html
# check:
# 1) blood/muscle/fat content (checked) - heatmap and density; found some outliers
# 2) PCA sample distance?
# 4) correlation among samples (TPM)?
# 3) heatmap of FSHD markers and top with clinical data 

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(pheatmap)

# define and load parameters
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
load(file.path(pkg_dir, "data", "dds.rda"))

annotation <- as.data.frame(rowData(dds)) %>%
  rownames_to_column(var="gene_id") 

#
# A) blood/muscle/fat content
#
markers <- tibble(cell_type=c(rep("blood", 3), rep("fat", 3), rep("muscle", 3)),
                  gene_name=c("HBA1", "HBA2", "HBB",
                              "FASN", "LEP", "SCD",
                              "ACTA1", "TNNT3", "MYH1"))  %>%
  dplyr::mutate(cell_type = factor(cell_type))

tmp_id <- annotation %>% dplyr::filter(gene_name %in% markers$gene_name)  %>%
  select(gene_id, gene_name)
  
markers <- markers %>% dplyr::left_join(tmp_id, by="gene_name")

# score for each cell_type: avg of log10(TPM+1)
celltype_tpm <- sapply(levels(markers$cell_type), function(type) {
   id <- markers %>% dplyr::filter(cell_type == type) %>%
      pull(gene_id)
   sub <- dds[id]
   tpm_score <- colSums(log10(assays(sub)[["TPM"]]+1))
   #tpm_score <- colMeans(log10(assays(sub)[["TPM"]]+1))
}) %>% as.data.frame() %>%
  rownames_to_column(var="sample_name")

blood_fat_muscle_content <- celltype_tpm
save(blood_fat_muscle_content, file=file.path(pkg_dir, "data", "blood_fat_muscle_content.rda"))

#
quantile(markertype_tpm$muscle, probs=seq(0, 1, 0.01))
quantile(markertype_tpm$blood, probs=seq(0, 1, 0.05))
quantile(markertype_tpm$fat, probs=seq(0, 1, 0.05))

# density: visualization of the distribution:
celltype_tpm %>% gather(key=cell_type, value=TPM, -sample_name) %>%
 ggplot(aes(x=TPM)) +
   geom_density() +
   facet_wrap(~cell_type) +
   theme_bw()
ggsave(file.path(pkg_dir, "figures", "biopsy-fat-blood-muscle-densty.pdf"), height=3, width=6)

# heatmap of the nine genes
library(pheatmap)
sub <- dds[markers$gene_id]
mat <- log10(assays(sub)[["TPM"]]+1)
rownames(mat) <- markers$gene_name
anno_row <- markers %>% column_to_rownames(var="gene_name") %>%
  dplyr::select(cell_type)
pheatmap(mat, annotation_row = anno_row, fontsize_col=5.5, fontsize_row=7, cluster_row=FALSE,
        file=file.path(pkg_dir, "figures", "biopsy-fat-blood-muscle-heatmap.pdf"), width=8, height=4)

# which samples have relatively more blood or fat content? Which has less muscle content
# normal range: blood < 4, fat < 1.5, muscle > 2.5
celltype_tpm %>% dplyr::filter(blood > 10)
celltype_tpm %>% dplyr::filter(fat > 4)
celltype_tpm %>% dplyr::filter(muscle < 2.5)

celltype_tpm %>% summarize(blood_mean=mean(blood), blood_sd = sd(blood),
                           fat_mean=mean(fat), fat_sd = sd(fat), 
                           muscle_mean=mean(muscle), muscle_sd=sd(muscle))

# conclusion: 
#   (a) 13-0007R is muscle-less.
#   (b) 32-0020R ad 32-0029R are bloody
#   (c) 13-0007R and 13-0009R are fattier

#
# correlation among samples
#


#
# B) PCA: can we find outliers (possibly the ones that are muscle-less, blooder and fatter)
#  use rlog data and DESeq2::plotPCA()     
#
library(ggrepel)
load(file.path(pkg_dir, "data", "rlog.rda"))
rlog <- rlog[, !colnames(rlog) == "32-0028L_S4"]
data <- plotPCA(rlog,
                intgroup=c("location"), returnData=TRUE)
                
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=location)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title="PCs on the BiLat RNA-seq muscle biopsies") +
  geom_text_repel(aes(label=rownames(data)), size=1.8, show.legend=FALSE) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(#legend.position="bottom",
        legend.title = element_blank(),
        legend.position = c(0.1, 0.1),
        legend.box.background = element_rect(colour = "black"),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=10))
ggsave(file.path(pkg_dir, "figures", "biopsy-PCA.pdf"), width=5, height=4)


# colored by DUX4+ (need to load rna_score)

data$`DUX4+` = rna_score$DUX4_targeted
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=`DUX4+`)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title="PCs on the BiLat RNA-seq muscle biopsies") +
  geom_text_repel(aes(label=rownames(data)), size=1.8, show.legend=FALSE) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(#legend.position="bottom",
        #legend.title = element_blank(),
        legend.position = c(0.1, 0.2),
        legend.box.background = element_rect(colour = "black"),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=10))
ggsave(file.path(pkg_dir, "figures", "biopsy-PCA.2.pdf"), width=5, height=4)

# conclusion:
#   (1) 32-0020L seems an outlier; find out its MRI/histology data
#   (2) 13-0009R and 13-0007R on the top side are fattier samples
#   (3) 32-0020R, a bloody sample, is away from the cluster
#   (4) the PC1 is driven by immunoglobulin H, L, K, and J
# ??? check the first principal component loading variable -> find out why 



#
# explore pricipal components and loading variables
#

# first pricipal component
.getPCA <- function(object, ntop=500) {
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
}
pca <- .getPCA(rlog)

# what drives PC1: higher value -> higher IG expression
loading <- as.data.frame(pca$rotation) %>%
  rownames_to_column(var="gene_id") %>% dplyr::arrange(desc(PC1))
# JCHAIN (immuglobulin J polypepetide), IGKC, IGHG1 IGHM, IGKV3-11 IGHV1-69-2
annotation %>% dplyr::filter(gene_id %in% loading$gene_id[1:20])

# what drived PC2: lower value -> higher IG and higher muscle expression
loading <- as.data.frame(pca$rotation) %>%
  rownames_to_column(var="gene_id") %>% dplyr::arrange(PC2)
annotation %>% dplyr::filter(gene_id %in% loading$gene_id[1:30])


#
# D) heatmap on DUX4/extracellular/inflamm/cell-cycle/immunoglobulin marker genes; use rlog
#
library(pheatmap)
load(file.path(pkg_dir, "data", "rlog.rda"))
markers_list <- list(dux4 = c("LEUTX", "KHDC1L", "PRAMEF2", "TRIM43"),
                     extracellular_matrix = c("PLA2G2A", "COL19A1", "COMP", "COL1A1"), # remove SFPR2
                     inflamm = c("CCL18", "CCL13" ,"C6", "C7"), # remove CCL19
                     cell_cycle = c("CCNA1", "CDKN1A", "CDKN2A"), # newly added
                     immunoglobulin = c("IGHA1", "IGHG4", "IGHGP"))
marker_id <- annotation %>% dplyr::filter(gene_name %in% unlist(markers_list))  %>%
  select(gene_id, gene_name)                     

marker_df <- unlist(markers_list) %>% as.data.frame() %>% rownames_to_column(var="type") %>%
  rename(gene_name=".") %>%
  dplyr::mutate(type = stringr::str_sub(type, 1L, -2L)) %>%
  dplyr::left_join(marker_id, by="gene_name")

data <- assay(rlog[marker_df$gene_id])
rownames(data) <- marker_df$gene_name
anno_row <- marker_df %>% dplyr::select(gene_name, type) %>% column_to_rownames(var="gene_name")
pheatmap(data, annotation_row=anno_row,
         fontsize_col=5.5, fontsize_row=7, cluster_row=FALSE, 
         file=file.path(pkg_dir, "figures", "biopsy-special-markers-heatmap.pdf"), width=8, height=4)

# how about markers on FSHD markers?        
load(file.path(pkg_dir, "data", "FSHD_markers.rda"))
sel <- intersect(rownames(rlog), FSHD_markers$gene_id)
data <- assay(rlog[sel])
rownames(data) <- rowData(rlog[sel])$gene_name
pheatmap(data, fontsize_col=5.5, fontsize_row=5.5,
         file=file.path(pkg_dir, "figures", "biopsy-FSHD-markers-heatmap.pdf"), width=6, height=4)

# conclusion:
#  1) 32-0020L and 32-0020R are high in DUX4 target expression
#  2) the outlier (PCA) 32-0022L is high in DUX4 target expression


#
# D) correlation among samples distance
#
sample_dists <- dist(t(assay(rlog)))
library("RColorBrewer")
sample_dist_matrix <- as.matrix(sample_dists)
#rownames(sample_dist_matrix) <- paste(log$condition, vsd$type, sep="-")
colnames(sample_dist_matrix) <- rownames(sample_dist_matrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sample_dist_matrix, fontsize_col=5.5, fontsize_row=5.5,
         clustering_distance_rows=sample_dists,
         clustering_distance_cols=sample_dists,
         col=colors, 
         file=file.path(pkg_dir, "figures", "biopsy-sample-distance-heatmap.pdf"), width=7, height=7)

# conclusion:
#  1) most of the paired-samples have higher correlation (cluster analysis by sample distance)
#  2) there is always a reason why some paried-samples aren't "paired": 
#     - 13-0009R and 13-0007R are not correlated with their paired because of their fat content.
#     - 32-0022L has high DUX4 target expression where 32-0022R doesn't not; they have different STIR and FAT MRI status
#     - 32-0027R: a bit high in some FSHD markers / 32-0027L: low in FSHD markers; this pair has different MRI and RNA profile