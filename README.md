# Wellstone_BiLateral_Biopsy

This repository is created to support the computational transparency and reproducibility of our published paper,_Validation study of the association of MRI and FSHD gene signature reveals markers of whole muscle and systemic disease progression_. It includes the preprocessed RNA-seq data, metadata, clinical scores, analysis results and a gitbook with details of our anlaysis and reproducible R code that generates figures and tables on the fly. 

## Structure
```
\data: clinial data, preprocessed RNA-seq datasets from the longigutinal and bilateral 
stuides and their metadata

    |- sanitized.dds.rda
    |- cluster_df.rda
    |- dds.rda
    |- comprehensive_df.rda

\docs: folder hosts the gitbook (*.html); all figures were generated on the fly of 
      the *.Rmd files
\gitbook: orignal *.Rmd that makes the gitbook
\scripts: un-orgnaized R code of our initial data exploration and analysis
\extdata: supplemental tables 

```
