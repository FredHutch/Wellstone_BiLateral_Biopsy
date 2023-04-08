# make_RSE.R
# This script (1) sort/index the bam file; (2) get gene counts and make RES instance
# R/4.1.2

library(GenomicFeatures)
library(GenomicAlignments)
library(hg38.HomoSapiens.Gencode.v35)
data(gene.anno)
library(tidyverse)
library(Rsamtools)
library(DESeq2)
library(BiocParallel)
library(ggrepel)
bp_param <- MulticoreParam(workers=12L)
register(bp_param)

pkg_dir <- file.path("/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.BiLat.FSHD.biopsy")
bam_dir <- file.path("/fh/fast/tapscott_s/user/sbennett/Wellstone_bx_08-19-21/bams/bams")
scratch_dir <- file.path("/fh/scratch/delete90/tapscott_s/hg38.BiLat.FSHD.biopsy/sort_index_bam")

setwd(file.path(pkg_dir, "scripts"))
#
# sort and index: stort in 
#
bam_files <- list.files(bam_dir, full.name=TRUE, pattern="\\.bam$")
for (file in bam_files) {
    cmt <- sprintf("sbatch -n4 ./sort_index.sh %s", file)
    system(cmt)
}
  

#
# sample info
#
sort_files <- list.files(scratch_dir, full.name=TRUE, pattern="\\.bam$")
# BAM info: library size
lib_size <- bplapply(sort_files, countBam, BPPARAM=bp_param)
lib_size_df <- do.call(rbind, lib_size)[, c("file", "records")] %>%
  dplyr::rename(sample_name = file, lib_size = records) %>%
dplyr::mutate(sample_name = str_replace(sample_name, "_R1_001.bam", ""))

# sample info data.frame
si <- data.frame(bam_file=sort_files) %>%
  dplyr::mutate(sample_name = str_replace(basename(bam_file), "_R1_001.bam", "")) %>%
  dplyr::mutate(Subject = str_replace_all(sapply(str_split(sample_name, "_"), "[[", 1), "L|R|b", "")) %>%
  dplyr::mutate(location = str_sub(str_replace(sapply(str_split(sample_name, "_"), "[[", 1), "b", ""), -1)) %>%
  dplyr::mutate(b_flag = str_detect(sample_name, "b")) %>%
  dplyr::left_join(lib_size_df, by="sample_name") %>%
  add_column(read_length = 100)
rownames(si) <- si$sample_name
si <- as(si, "DataFrame")

# check if the bam files are single or paired-ends
library(Rsamtools)
quickBamFlagSummary(si$bam_file[1], main.groups.only=TRUE)

#
# gene counts: this should not include 32-0028L
#
features <- exonsBy(hg38.HomoSapiens.Gencode.v35, by="gene")
features <- keepStandardChromosomes(features,
                                    species = "Homo_sapiens",
                                    pruning.mode="fine")
rse <- GenomicAlignments::summarizeOverlaps(features=features,
                                            reads = Rsamtools::BamFileList(si$bam_file),
                                            mode = "IntersectionStrict",
                                            ignore.strand=TRUE, singleEnd=TRUE,
                                            inter.feature=TRUE, BPPARAM=bp_param)
colnames(rse) <- str_replace(colnames(rse), "_R1_001.bam", "")
SummarizedExperiment::colData(rse) <- si[colnames(rse), ]
SummarizedExperiment::rowData(rse)$gene_name <- gene.anno[rownames(rse), "gene_name"]
SummarizedExperiment::rowData(rse)$gene_type <- gene.anno[rownames(rse), "gene_type"]

# annotation
rowData(rse)$gene_name <- gene.anno[rownames(rse), "gene_name"]
rowData(rse)$gene_type <- gene.anno[rownames(rse), "gene_type"]

# TPM and RPKM
source(file.path(pkg_dir, "scripts", "getTPMperTx.R"))
assays(rse)$TPM <- getTPMperTx(rse)
source(file.path(pkg_dir, "scripts", "getScaledCounts.R"))
assays(rse)$RPKM <- getScaledCountsPerTx(rse)

library(DESeq2)
dds <- DESeqDataSet(rse, design = ~ location)
dds <- DESeq2::estimateSizeFactors(dds)
save(dds, file=file.path(pkg_dir, "data", "dds.rda"))

# tidy up dds to correct the colnames and sync with sample_id: bilat_dds
bilat_dds <- dds
bilat_dds <- bilat_dds[rowSums(counts(dds)) >= 30 ]
bilat_dds <- DESeq2::estimateSizeFactors(bilat_dds)
bilat_dds$sizeFactors
colnames(bilat_dds) <- bilat_dds$sample_id <- str_replace(bilat_dds$sample_name, "[b]*_.*", "")
save(bilat_dds, file=file.path(pkg_dir, "data", "bilat_dds.rda"))


sub <- dds[rowSums(counts(dds)) >= 30]
rlog <- rlog(sub)
save(rlog, file=file.path(pkg_dir, "data", "rlog.rda"))

#
# rlog with 32-0028L
#
library(DESeq2)
dds <- DESeqDataSet(rse, design = ~ location)
dds <- DESeq2::estimateSizeFactors(dds)
sub <- dds[rowSums(counts(dds)) >= 30]
rlog <- rlog(sub)

tmp_dir <- file.path("/fh/scratch/delete90/tapscott_s/hg38.BiLat.FSHD.biopsy")
save(rlog, file=file.path(tmp_dir, "data", "rlog_w_outlier.rda"))
save(dds, file=file.path(tmp_dir, "data", "dds_w_outlier.rda"))
# PCA with the outlier
data <- DESeq2::plotPCA(rlog, intgroup=c("location", "sample_name"), returnData=TRUE)

percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=location)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title="PCA: BiLat RNA-seq with 32-0028L") +
  geom_text_repel(aes(label=sample_name), size=1.8, show.legend=FALSE) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(#legend.position="bottom",
        legend.title = element_blank(),
        legend.position = c(0.9, 0.9),
        legend.box.background = element_rect(colour = "black"),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=10))
ggsave(file.path(tmp_dir, "figures", "biopsy-PCA-w-32-0028L.pdf"), width=5, height=4)


#
# QC: exclude 32-0028L, use loading variables to prove it is a fatty tissue
#
#tmp_dir <- file.path("/fh/scratch/delete90/tapscott_s/hg38.BiLat.FSHD.biopsy")
tmp_dir <- pkg_dir
load(file.path(tmp_dir, "data", "rlog_w_outlier.rda"))
load(file.path(tmp_dir, "data", "dds_w_outlier.rda"))
anno_ens88 <- as.data.frame(rowData(rlog)) %>%
  rownames_to_column(var="gene_id") %>%
  dplyr::select(gene_id, gene_name, gene_type)


# (a) get PCA tool
getPCA <- function(object, ntop=500) {
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
}

pca <- getPCA(rlog, ntop=500)
pcs <- pca$x
loading_var <- as.data.frame(pca$rotation)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
names(percentVar) <- colnames(pca$x)

# (b) get firt 40 first loading variables 
loading_var <- as.data.frame(pca$rotation) %>%
  rownames_to_column(var="gene_id") %>%
  dplyr::select(gene_id, PC1, PC2) %>%
  left_join(anno_ens88, by="gene_id") %>%
  arrange(desc(PC1))

loading_var[1:45, ]

# viz by heatmap; top 20
library(pheatmap)
data <- assay(rlog)[loading_var$gene_id[1:50],]
rownames(data) <- loading_var$gene_name[1:50]
pheatmap(data,
 cellheight=6, fontsize = 6, silent=TRUE, scale="row",
         filename=file.path(tmp_dir, "figures", "32-0028L-heatmap-top-50-loading-variables-zscore.pdf"))

data <- assay(rlog)[loading_var$gene_id[1:80],]
rownames(data) <- loading_var$gene_name[1:80]
pheatmap(data,  show_rownames=FALSE,
 fontsize = 6, silent=TRUE, scale="row",
         filename=file.path(tmp_dir, "figures", "32-0028L-heatmap-top-80-loading-variables-zscore.pdf"))
#
# check their FAT/muscle content; this should not be the reason
#

assays(dds)$TPM <- getTPMperTx(dds) # function from the scripts folder
markers <- tibble(cell_type=c(rep("blood", 3), rep("fat", 3), rep("muscle", 3)),
                  gene_name=c("HBA1", "HBA2", "HBB",
                              "FASN", "LEP", "SCD",
                              "ACTA1", "TNNT3", "MYH1"))  %>%
  dplyr::mutate(cell_type = factor(cell_type)) %>%
  left_join(anno_ens88, by="gene_name")


# (1) heatmap using rlog
library(pheatmap)
data <- assay(rlog[markers$gene_id, ])
rownames(data) <- markers$gene_name
pheatmap(data, 
         cellheight=12, fontsize = 6, silent=TRUE, scale="row", 
         filename=file.path(tmp_dir, "figures", "32-0028L-heatmap-fat-blood-muscel-content.pdf"),
         width=6, height=6)

# (2) boxplot using content score

celltype_tpm <- sapply(levels(markers$cell_type), function(type) {
   id <- markers %>% dplyr::filter(cell_type == type) %>%
      pull(gene_id)
   sub <- dds[id]
   tpm_score <- colSums(log10(assays(sub)[["TPM"]]+1))
   #tpm_score <- colMeans(log10(assays(sub)[["TPM"]]+1))
}) %>% as.data.frame() %>%
  rownames_to_column(var="sample_name")

celltype_tpm %>%
  gather(key=content, value=score, -sample_name) %>%
  dplyr::mutate(content=factor(content)) %>%
  dplyr::mutate(colored = if_else(sample_name=="32-0028L_S4", FALSE, TRUE)) %>%
  ggplot(aes(x=content, y=score)) +
    geom_boxplot(width=0.5, outlier.shape=NA) +
    geom_jitter(width=0.3, size=1, aes(color=colored, alpha=0.7)) +
    scale_color_manual(values=c("red", "black")) +
    theme_bw() +
    #facet_wrap(~content, scale="free_y") +
    theme(legend.position="none") 

ggsave(file.path(pkg_dir, "figures", "32-0028L-boxplot-fat-blood-muscle-conent.pdf"),
       width=2, height=2)    


