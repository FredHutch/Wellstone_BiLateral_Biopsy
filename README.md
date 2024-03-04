# Wellstone_BiLateral_Biopsy

This repository is created to support the computational transparency and reproducibility of our published paper,_Validation study of the association of MRI and FSHD gene signature reveals markers of whole muscle and systemic disease progression_. It includes the preprocessed RNA-seq data, metadata, clinical scores, analysis results and a gitbook with details of our anlaysis and reproducible R code that generates figures and tables on the fly. 

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
