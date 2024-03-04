# Wellstone_BiLateral_Biopsy

This repository is created to support the computational transparency and reproducibility of our published paper, [_Regional and bilateral MRI and gene signatures in facioscapulohumeral dystrophy: implications for clinical trial design and mechanisms of disease progression_](https://doi.org/10.1093/hmg/ddae007). It includes the preprocessed RNA-seq data, metadata, histology scores, MRI characteristics, scripts, and a [GitBook](https://fredhutch.github.io/Wellstone_BiLateral_Biopsy/index.html) that provides a narrative on our analytical processes. This [GitBook](https://fredhutch.github.io/Wellstone_BiLateral_Biopsy/index.html) includes reproducible R code, enabling the generation of figures and tables directly from the code. 

__Note:__ For replicating the analysis results, the reader should substitute the `pkg_dir` parameter with the path to the location where the clone of this repository is saved.

## Structure
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
