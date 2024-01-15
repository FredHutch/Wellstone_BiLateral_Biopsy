--- 
title: "Wellstone bilateral tibialis anterior FSHD biopsied muscle study"
author: "Chao-Jen Wong"
date: "2024-01-13"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
cover-image: images/cover.jpeg
url: https://fredHutch.github.io/Wellstone_BiLateral_Biopsy
description: |
  This book supports computational transparency and reproducibility for our publication titled "Validation study of the association of MRI and FSHD gene signature reveals markers of whole muscle and systemic disease progression".
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---

# About

The objective of this book is to advocate for computational transparency and reproducibility of our FSHD bilateral cohort study, titled ["Validation study of the association of MRI and FSHD gene signature reveals markers of whole muscle and systemic disease progression."](https://www.biorxiv.org/content/biorxiv/early/2023/02/20/2023.02.20.529303.full.pdf) This book includes detailed explanations of our decision-making processes grounded in bioinformatics analysis, machine learning, and statistical inferences along with reproducible codes (in R) that generate results and figures on the fly. 

In this FSHD bilateral cohort study, we recruited 34 FSHD patients and obtained their left and right tibialis anteria (TA) muscle biopsies. On the biopsies, we performed MRI characterization, histopathological scoring, RNA extraction for RNA sequencing, and DNA extraction for bisulfite sequencing of the 4qA permissive allele in the region of exon 2-3. Through the metrics yield by these processes, we studied the following topics:

1. Transcript-based assessment of muscle cell type content of the FSHD biopsied muscle
2. FSHD molecular signatures, including DUX4, extracurricular matrix, inflammatory, complement activation, and immunoglobulins
3. ML classification based on the FSHD molecular signatures
4. Identification of six baskets containing genes that exhibit specific signatures associated with FSHD
5. Association of MRI characteristics with signatures of DUX4 expression 
6. Bilateral comparison analysis
7. Verification of immune cell infiltrate signatures and immune cell type proportions


## Software

- R (4.2.2) and packages from the Bioconductor (3.17) and Tidyverse projects. Most used packages include `GenomicAlignments`, `GenomicFeatures`, `DESeq2`, `dplyr`, `ggplot2`, `caret`, and genomic-related BioC packages for annotation.
- `fastqc`, `Subread`, `samTools`, and `cutadapt` for RNA-seq data preprocessing, and `PLIER` for immune cell type proportions.

## Repos structure
The section outlines the structure of our repository, which is available at [here](https://github.com/FredHutch/Wellstone_BiLateral_Biopsy). 

```
\data: 
    |- longitudinal_dds.rda: a DESeqDataSet instance include the 
       longitudinal RNA-seq gene counts, TPM, and published metadata 
       including ML-based classification labels, MRI, histopathology 
       scores, and clinical data
    |- bilateral_dds.rda: a DESeqDataSet instance; bilateral RNA-seq gene
       counts, TPM, and metadata including ML-based classification labels 
       and published MRI and clinical data
    |- comprehensive_df.rda: a data.frame instance obtaining the biletarl
       cohort's clinical, MRI, pathology, FSHD molecular signature scores,
       and DNA methylation data
    |- all_baskets.rda: basekts of genes representing FSHD disease 
       signatures (DUX4-M, DUX4-M6, DUX4-M12, ECM, Inflamm, Complement,
       and IG baskets)
    |- bilat_MLpredict.rda: categorization of the bilateral biopsy
       samples by machine learning models (random forest and 
       KNN) trained by the longitudinal cohort and the signature
       genes in the baskets
       cohort 
       

\docs: *.html, the pages rendered by *.Rmd files in the \gitbook
       directory
\gitbook: *.Rmd, source code of this book
\scripts: *.R; un-orgnaized R code of our initial data 
          exploration and bioinformatics analysis
\extdata: *.xlsx; extra files and supplemental tables
```

## Annotation 
The gene annotation for our RNA-seq analysis is based on Gencode version 35. To facilitate gene counts using Bioconductor packages, we customized a Bioconductor `TxDB` package named `hg38.HomoSapiens.Gencode.v35`. See Appendix \@ref(appendixA-buildTxDB) on how to build Bioconductor `TxDB` or `EsbDB` annotation packages.


## RNA-seq preprocessing and gene counts
The pre-preprocessing of RNA-seq data, we filtered out unqualified raw reads, trimmed the Illumina universal adapters using _Trimmomatic_, and aligned the remaining reads to GRCh38.p13 with _Rsubread_.

The code chunk below provides an example of how we used `GenomicAlignments`, `BiocParallel`, and a customized `TxDB` package (`hg38.HomoSapiens.Gencode.v35`) to count reads and generate a `RangedSummarizedExperiment` object. We used `GenomicAlignments::summarizeOveralsp()` to count the reads. With `Intersectionstrict` mode that we only count reads that are completely contained within the range of exons, and ignore any ambiguous reads that straddle different gene features. 

- `sort_files` parameter below indicates the location of the sorted, indexed BAM files used in our analysis


```r
# sort_files gives the location of our bam files
sort_files <- list.files(scratch_dir, full.name=TRUE, pattern="\\.bam$")

library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(BiocParallel)
library(hg38.HomoSapiens.Gencode.v35)
bp_param <- BiocParallel::MulticoreParam(workers=12L)
register(bp_param)

features <- GenomicFeatures::exonsBy(hg38.HomoSapiens.Gencode.v35, by="gene")
features <- keepStandardChromosomes(features,
                                    species = "Homo_sapiens",
                                    pruning.mode="fine")
rse <- 
  GenomicAlignments::summarizeOverlaps(features=features,
                                       reads =
                                         Rsamtools::BamFileList(sort_file),
                                       mode = "IntersectionStrict",
                                       ignore.strand=TRUE, 
                                       singleEnd=TRUE,
                                       inter.feature=TRUE, 
                                       BPPARAM=bp_param)
```

