# Immune cell infiltrates signature {#immune-cell-infiltrates}


Please note that in the longitudinal study, we employed K-means clustering to identify five distinct clusters labeled as Mild, Moderate, IG-High, High, and Muscle-Low. Through this study, we have found that T cell proliferation, migration, and B cell-mediated immune responses are more prominent in severely affected muscles, specifically those labelled as "IG-high" and "High".  The objective of this chapter is as follows:

1. Identify the enriched GO terms (biological processes) directly linked to T/B cells and circulating immunoglobulins in the longitudinal study (High or IG0high vs. control).
2. Determine the differentially expressed genes (High vs. control) associated with the enriched GO terms (from step 1), and use heatmap visualization to observe the trend of an increase in gene expression across all longitudinal and bilateral samples, associated with the severity of FSHD (characterised by STIR status, fat fraction/infiltration percent) and/or complement scoring.
3. Identify immune cell markers that are differentially expressed (High vs. control) using the IRIS database (Abbs 2009) as a source.




```r
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
fig_dir <- file.path(pkg_dir, "figures", "immune-infiltration")
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

# color pal
col <- c("#ccece6", "#006d2c") # control vs. FSHD?
stir_pal <- c("STIR+" = "#006d2c", "STIR-" = "#ccece6")
complement_pal = c("3" = "#006d2c", "2" = "#66c2a4", "1" = "#ccece6")
```

## Enriched biological processes directly associated with T/B cells
The code chunk below (1) obtains the enriched GO terms in more affected muscles from the longitudinal study, supplementary table 5, and (2) retains the enriched GO terms containing the terms including `T cell`, `B cell`, `humoral`, `immunoglobulin`, and `complement`. Total 20 terms are extracted. 


```r
high_enriched <- readxl::read_xlsx(file.path(pkg_dir, "extdata", 
    "suppl_table_5_Candidate_Biomarkers_and_Enriched_Go.xlsx"),
    sheet=2, skip=2)

ig_enriched <- readxl::read_xlsx(file.path(pkg_dir, "extdata", 
    "suppl_table_5_Candidate_Biomarkers_and_Enriched_Go.xlsx"),
    sheet=4, skip=2)

moderate_enriched <- readxl::read_xlsx(file.path(pkg_dir, "extdata", 
    "suppl_table_5_Candidate_Biomarkers_and_Enriched_Go.xlsx"),
    sheet=6, skip=2)

# get the lymphocyte related GO ID
infiltrates_go <- high_enriched %>%
  dplyr::filter(str_detect(term, "T cell|B cell|humoral|immunoglobulin|complement"))

infiltrates_go %>% 
  dplyr::select(category, over_represented_pvalue, term) %>%
  dplyr::rename(pvalue = over_represented_pvalue) %>%
  dplyr::mutate(pvalue = format(pvalue, scientific = TRUE)) %>%
  dplyr::arrange(category) %>%
  kbl(caption="Enriched biological process GO terms (n=20)  related to T/B cells, humoral immune response, and immunoglobulin domains. Identified using the longitudinal IG-High and High samples (vs. controls).") %>%
  kable_paper("hover", full_width = F)
```

<table class=" lightable-paper lightable-hover" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption>(\#tab:sup-table-5)Enriched biological process GO terms (n=20)  related to T/B cells, humoral immune response, and immunoglobulin domains. Identified using the longitudinal IG-High and High samples (vs. controls).</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> category </th>
   <th style="text-align:left;"> pvalue </th>
   <th style="text-align:left;"> term </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0002455 </td>
   <td style="text-align:left;"> 4.391414e-10 </td>
   <td style="text-align:left;"> humoral immune response mediated by circulating immunoglobulin </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0002460 </td>
   <td style="text-align:left;"> 6.665762e-09 </td>
   <td style="text-align:left;"> adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0002920 </td>
   <td style="text-align:left;"> 2.865839e-11 </td>
   <td style="text-align:left;"> regulation of humoral immune response </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006956 </td>
   <td style="text-align:left;"> 4.086619e-10 </td>
   <td style="text-align:left;"> complement activation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006957 </td>
   <td style="text-align:left;"> 1.965806e-04 </td>
   <td style="text-align:left;"> complement activation, alternative pathway </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006958 </td>
   <td style="text-align:left;"> 1.619492e-10 </td>
   <td style="text-align:left;"> complement activation, classical pathway </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006959 </td>
   <td style="text-align:left;"> 1.158839e-15 </td>
   <td style="text-align:left;"> humoral immune response </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0010818 </td>
   <td style="text-align:left;"> 2.450256e-09 </td>
   <td style="text-align:left;"> T cell chemotaxis </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0016064 </td>
   <td style="text-align:left;"> 1.175922e-08 </td>
   <td style="text-align:left;"> immunoglobulin mediated immune response </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0019724 </td>
   <td style="text-align:left;"> 1.380039e-08 </td>
   <td style="text-align:left;"> B cell mediated immunity </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0030449 </td>
   <td style="text-align:left;"> 7.071256e-11 </td>
   <td style="text-align:left;"> regulation of complement activation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0031295 </td>
   <td style="text-align:left;"> 1.531739e-04 </td>
   <td style="text-align:left;"> T cell costimulation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0042098 </td>
   <td style="text-align:left;"> 5.735434e-06 </td>
   <td style="text-align:left;"> T cell proliferation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0042102 </td>
   <td style="text-align:left;"> 4.877625e-05 </td>
   <td style="text-align:left;"> positive regulation of T cell proliferation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0042110 </td>
   <td style="text-align:left;"> 6.083250e-07 </td>
   <td style="text-align:left;"> T cell activation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0042129 </td>
   <td style="text-align:left;"> 1.273162e-05 </td>
   <td style="text-align:left;"> regulation of T cell proliferation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0050863 </td>
   <td style="text-align:left;"> 3.541698e-08 </td>
   <td style="text-align:left;"> regulation of T cell activation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0050870 </td>
   <td style="text-align:left;"> 1.090819e-06 </td>
   <td style="text-align:left;"> positive regulation of T cell activation </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0072678 </td>
   <td style="text-align:left;"> 1.727781e-11 </td>
   <td style="text-align:left;"> T cell migration </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:2000404 </td>
   <td style="text-align:left;"> 1.335915e-05 </td>
   <td style="text-align:left;"> regulation of T cell migration </td>
  </tr>
</tbody>
</table>

Display the p-values yielded by GSEA:

```r
df <- infiltrates_go %>%
  dplyr::select(category, term, over_represented_pvalue) %>%
  arrange(term) %>% 
  left_join(dplyr::select(ig_enriched, category, over_represented_pvalue),
            by="category", suffix=c(".High", ".IG-High")) %>%
  tidyr::gather(key=class, value=pvalue, -category, -term) %>%
  dplyr::mutate(class = str_replace(class, 
                                    "over_represented_pvalue.", ""),
                term = str_trunc(term, 45),
                class = factor(class, levels=c("IG-High", "High"))) %>%
  dplyr::mutate(log10pvalue = -log10(pvalue)) 

cat_levels <- df %>% dplyr::filter(class == "High") %>%
  dplyr::arrange(category) %>% pull(term)

df %>% dplyr::mutate(term=factor(term, levels=cat_levels)) %>%
  #dplyr::mutate(term = factor(term,levels=cat_levels)) %>%
  ggplot(aes(x=class, y=term)) +
    geom_point(aes(size=log10pvalue), color="red1", alpha=0.7,
               show.legend = FALSE) +
    theme_minimal() +
  labs(title="Enriched GO terms related to T/B cells and immunoglobulins") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust=1, size=9.5),
        axis.title.y = element_blank())
```

<div class="figure">
<img src="08-immune-cell-infiltrates-signature_files/figure-html/infiltrates-go-p-values-1.png" alt="Enriched GO terms related to T and B cells and immunoglobulin domains. Red dots present the p-value of enriched GO terms in IG or High samples." width="672" />
<p class="caption">(\#fig:infiltrates-go-p-values)Enriched GO terms related to T and B cells and immunoglobulin domains. Red dots present the p-value of enriched GO terms in IG or High samples.</p>
</div>

```r

ggsave(file.path(fig_dir, "longitudinal-enriched-immunity-GO.pdf"), 
       width=3.8, height=4.2)  
```

## Significant genes associated with the T/B-associated enriched GO terms
The code chunk below identifies significantly expressed genes in High or IG-High (compared to the controls) that are associated with the enriched GO terms obtained previously.

__Results:__ Total 95 B/T cell immunity response significant markers are identified. 


```r
high_sig <- readxl::read_xlsx(file.path(pkg_dir, "extdata", 
    "suppl_table_5_Candidate_Biomarkers_and_Enriched_Go.xlsx"),
    sheet=1, skip=3) %>%
  dplyr::select(gene_name, gencode_id) %>%
  dplyr::mutate(ens_id = str_replace(gencode_id, "\\..*", ""))

ig_sig <- readxl::read_xlsx(file.path(pkg_dir, "extdata", 
    "suppl_table_5_Candidate_Biomarkers_and_Enriched_Go.xlsx"),
    sheet=3, skip=3) %>%
  dplyr::select(gene_name, gencode_id) %>%
  dplyr::mutate(ens_id = str_replace(gencode_id, "\\..*", ""))
sig <- bind_rows(high_sig, ig_sig) %>%
  dplyr::distinct(gene_name, .keep_all=TRUE)

## keep the genes that involved in the infiltrates_go
gene2cat <- goseq::getgo(sig$ens_id, "hg38", "ensGene", fetch.cats = "GO:BP")
names(gene2cat) <- sig$ens_id
cat2gene <- goseq:::reversemapping(gene2cat)[infiltrates_go$category]
infiltrates_act <- map_dfr(cat2gene, function(cat_gene) {
  sig %>% dplyr::filter(ens_id %in% cat_gene)
}, .id="category") %>% 
  arrange(category) %>%
  left_join(dplyr::select(infiltrates_go, category, term), by="category") %>%
  dplyr::rename(GOID = category) %>%
  distinct(gene_name, .keep_all=TRUE) %>%
  dplyr::rename(gene_id = gencode_id) %>%
  #arrange(gene_name) %>%
  dplyr::mutate(category = case_when(
    str_detect(term, "T cell") ~ "T cell migration / activation",
    str_detect(term, "humoral|complement") ~ "B cell mediated / humoral immune response",
    str_detect(term, "adaptive") ~ "adaptive immune response (IG domains)"
  )) %>% 
  arrange(category, gene_name)
```

### Heatmap: Longitudinal study
We have employed a heatmap to visualize the trend of increased expression levels for the 95 immune signature genes previously identified to be associated with enriched T/B and IG Gene Ontology (GO) terms. The order of the samples in the heatmap is arranged as Control, Mild, Moderate, IG-High, High, and Muscle-low. 


```r
cluster_color <- c(Control="#ff7f00", 
                    Mild="#a65628", 
                    Moderate="#f781bf", 
                    `IG-High`="#984ea3", High="#e41a1c", 
                   `Muscle-Low`="#377eb8")
category_color <- wes_palette("Darjeeling2", n=4)[c(1,3,4)]
names(category_color) <- levels(factor(infiltrates_act$category))
col <- c("#ccece6", "#006d2c")
stir_pal <- c("STIR+" = "#006d2c", "STIR-" = "#ccece6")
complement_pal = c("3" = "#006d2c", "2" = "#66c2a4", "1" = "#ccece6")
# column annotation: class, fat fraction, STIR rating, and histophatology

annotation_col <- colData(longitudinal_dds) %>% as.data.frame() %>%
  dplyr::select(cluster, fat_fraction, STIR_status, 
                #STIR_rating,
                Pathology.Score) %>%
  dplyr::rename(class=cluster, `STIR+/-` = STIR_status,
                `fat fraction`=fat_fraction,
                `histopathology` = Pathology.Score) %>%
  dplyr::mutate(`STIR+/-` = factor(`STIR+/-`))
# column annotation color
ann_cor <-  list(class = cluster_color, category = category_color,
                 `STIR+/-` = stir_pal,
                 `fat fraction` = col, 
                 `histopathology` = col)
markers_heatmp_group_by_types_class(markers = infiltrates_act, 
                                    dds = longitudinal_dds, 
                                    factor = "cluster", 
                                    annotation_col = annotation_col,
                                    ann_cor = ann_cor, 
                                    annotation_legend=TRUE)
```

<div class="figure">
<img src="08-immune-cell-infiltrates-signature_files/figure-html/heatmap-immune-infiltrate-longi-1.png" alt="Heatmap demonstrating the increased expression trend for 95 previously identifed T/B and IG domain related genes on more affected muscle. Color presented the row-wised z-score of expression levels." width="768" />
<p class="caption">(\#fig:heatmap-immune-infiltrate-longi)Heatmap demonstrating the increased expression trend for 95 previously identifed T/B and IG domain related genes on more affected muscle. Color presented the row-wised z-score of expression levels.</p>
</div>



### Heatmap:  bilateral study
The code chunk below generates a heatmap depicting the expression levels of the 95 immune signature genes from the bilateral samples. Noted that the samples are ordered based on their classification into Control-like and Moderate+, determined by ML-based classification, and Muscle-Low based on the low expression levels of muscle markers.


```r
annotation_col <- comprehensive_df %>% 
  drop_na(class) %>%
  dplyr::select(sample_id, class, Fat_Infilt_Percent, STIR_status,
                `Cumulative Score`, `Complement Scoring`) %>%
  dplyr::rename(`fat percent` = Fat_Infilt_Percent, 
                `STIR+/-` = STIR_status,
                `histopathology` = `Cumulative Score`) %>%
  column_to_rownames(var="sample_id") %>%
  dplyr::mutate(`Complement Scoring` = factor(`Complement Scoring`, 
                                              levels=c("3", "2", "1")))
                
class_color <- c(`Control-like`="#a65628", 
                 `Moderate+`="#984ea3", 
                 `Muscle-Low`="#377eb8")
ann_cor <- list(class = class_color, category = category_color,
                `fat percent` = col, 
                `STIR+/-` = stir_pal,
                `histopathology` = col, 
                `Complement Scoring` = complement_pal)

markers <- infiltrates_act %>%
  left_join(dplyr::select(anno_gencode35, ens_id, gencode35_id),
            by="ens_id") %>%
  dplyr::select(-gene_id) %>%
  dplyr::rename(gene_id = gencode35_id) %>%
  drop_na(gene_id)

markers_heatmp_group_by_types_class(markers = markers, 
                                    dds = bilat_dds, 
                                    factor = "class", 
                                    annotation_legend=TRUE,
                                    annotation_col = annotation_col,
                                    border_color = "transparent",
                                    ann_cor = ann_cor)
```

<img src="08-immune-cell-infiltrates-signature_files/figure-html/heatmap-immune-infiltrate-bilat-1.png" width="768" />


## Immune cell infiltrate signature
Using a generic cell-type marker dataset from the immune response in silico (IRIS) (Abbas 2009) and longitudinal cohort, we detected 63 immune cell markers whose expression levels are significantly elevated in more affected muscles. With additional significant immunity-related markers, including such as cytokine receptors in B cells (IL-4R, IL-21R), STAT6, Leukocyte immunoglobulin-like receptor B family (LILRB4, LILRB5), T cell surface antigens (CD2 and CD4), T and B cell interaction and migration (CCL19), and STAT6, we demonstrated a immune cell infiltrate signatures in both longitudinal and bilateral samples.
                                   
Code chunk below identified 63 immune cell markers as a subset of the IRIS database (Abbas 2009).

```r
library(PLIER)
data(bloodCellMarkersIRISDMAP)
# about 860 markers
blood_iris <- as.data.frame(bloodCellMarkersIRISDMAP) %>%
  rownames_to_column(var="gene_name") %>%
  dplyr::select(gene_name, starts_with("IRIS")) %>%
  dplyr::filter(rowSums(across(where(is.numeric))) > 0) 

# 63 are DE in High, longitudinal cohort
sig_blood_iris <- high_sig %>%
  inner_join(blood_iris, by="gene_name") 
category <- names(sig_blood_iris)[4:ncol(sig_blood_iris)]
sig_blood_iris$category <- 
  apply(sig_blood_iris[, 4:ncol(sig_blood_iris)], 1, 
        function(x) {
    tmp <- if_else(unlist(as.vector(x)) == 1, TRUE, FALSE)
    return(category[tmp][1])
       })

sig_blood_iris  <- sig_blood_iris %>%
  dplyr::select(gene_name, gencode_id, ens_id, category) %>%
  dplyr::rename(gene_id = gencode_id) %>%
  dplyr::mutate(category = str_replace(category, "IRIS_", "")) %>%
  dplyr::mutate(category = str_replace(category, "\\-.*", "")) %>%
  arrange(category, gene_name) %>%
  dplyr::mutate(category = factor(category, 
                                  levels = c("Bcell", "CD4Tcell", "CD8Tcell",
                                             "DendriticCell",
                                             "Monocyte", "Neutrophil",
                                             "NKcell", "PlasmaCell")))

# modification: do not consider those "other" markers
# add others: IL-4R, IL-21R, LILRB4, LILRB5, CD2, CD4, CCL19
#others <- data.frame(gene_name = c("IL4R", "IL21R", "LILRB4",
#                                   "LILRB5",
#                                   "CD2", "CD4", "CCL19", "STAT6")) %>%
#  left_join(anno_ens88, by="gene_name") %>% 
#  dplyr::rename(gene_id = ens88_id) %>%
#  add_column(category="Others") %>%
#  arrange(gene_name)
```

### Heamap: the longitudinal study
The following heatmap illustrates that there is a tendency of higher expression levels in muscles that are more affected.


```r
cluster_color <- c(Control="#ff7f00", 
                    Mild="#a65628", 
                    Moderate="#f781bf", 
                    `IG-High`="#984ea3", High="#e41a1c", 
                   `Muscle-Low`="#377eb8")


annotation_col <- colData(longitudinal_dds) %>% as.data.frame() %>%
  dplyr::select(cluster, fat_fraction, STIR_status, 
                Pathology.Score) %>%
  dplyr::rename(class=cluster, 
                `fat fraction`=fat_fraction,
                `STIR+/-` = STIR_status,
                `histopathology` = Pathology.Score) %>%
  dplyr::mutate(`STIR+/-` = factor(`STIR+/-`))

category_color <- wes_palette("Darjeeling1", n=8, typ="continuous")[1:8]
names(category_color) <- levels(sig_blood_iris$category)

ann_cor <-  list(class = cluster_color, 
                 category = category_color,
                 `fat fraction` = col, 
                 `STIR+/-` = stir_pal,
                 `histopathology` = col)
# heatmap
markers_heatmp_group_by_types_class(markers = sig_blood_iris, 
                                    dds = longitudinal_dds, 
                                    factor = "cluster", 
                                    border_color = "transparent",
                                    annotation_col = annotation_col,
                                    ann_cor = ann_cor)
```

<div class="figure">
<img src="08-immune-cell-infiltrates-signature_files/figure-html/longitudinal-63-sig-immune-cell-markers-heatmap-1.png" alt="Heapmap demonstrating the expression levels (by row-wise z-score of log10(TPM++1) in the longitudinal study." width="672" />
<p class="caption">(\#fig:longitudinal-63-sig-immune-cell-markers-heatmap)Heapmap demonstrating the expression levels (by row-wise z-score of log10(TPM++1) in the longitudinal study.</p>
</div>


### Heatmap: Bilateral study
The heatmap for the bilateral study below incorporates the complement scoring.

```r
annotation_col <- comprehensive_df %>% 
  drop_na(class) %>%
  dplyr::select(sample_id, class, Fat_Infilt_Percent, STIR_status,
                `Cumulative Score`, `Complement Scoring`) %>%
  dplyr::rename(`fat percent` = Fat_Infilt_Percent, 
                `STIR+/-` = STIR_status,
                `histopathology` = `Cumulative Score`) %>%
  column_to_rownames(var="sample_id") %>%
  dplyr::mutate(`Complement Scoring` = factor(`Complement Scoring`))

class_color <- c(`Control-like`="#a65628", 
                    `Moderate+`="#984ea3", 
                    `Muscle-Low`="#377eb8")
ann_cor <- list(class = class_color, category = category_color,
                `fat percent` = col, 
                `STIR+/-` = stir_pal,
                `histopathology` = col, 
                `Complement Scoring` = complement_pal)

# update sig_blood_iris gene_id to gencode35_id
tmp <- sig_blood_iris %>%
  left_join(dplyr::select(anno_gencode35, ens_id, gencode35_id), 
            by="ens_id") %>%
  dplyr::select(-gene_id) %>%
  dplyr::rename(gene_id = gencode35_id)

markers_heatmp_group_by_types_class(markers = tmp,
                                    dds = bilat_dds, 
                                    factor = "class", 
                                    border_color = "transparent",
                                    annotation_col = annotation_col,
                                    ann_cor = ann_cor )
```

<img src="08-immune-cell-infiltrates-signature_files/figure-html/bilat-63-sig-immune-cell-markers-heatmap-1.png" width="672" />


### DE immune marker in the Bilat corhort
We compared the Control-like and Moderate+ samples in the Bilat cohort to determine which the immune cell-type genes in the IRIS dataset are differentially expressed. 


```r
sub_dds <- bilat_dds[, bilat_dds$class %in% c("Control-like", "Moderate+")]
sub_dds$class <- factor(sub_dds$class)
design(sub_dds) <- ~ class
sub_dds <- DESeq(sub_dds)
```


```r
immune_signature_gene <- sig_blood_iris %>%
  left_join(dplyr::select(anno_gencode35, ens_id, gencode35_id), 
            by="ens_id") %>%
  dplyr::select(-gene_id) %>%
  dplyr::rename(gene_id = gencode35_id)

res <- results(sub_dds, lfcThreshold = 0.5) %>% 
  as.data.frame() %>%
  rownames_to_column(var="gene_id") %>%
  dplyr::filter(padj < 0.05, log2FoldChange > 0) %>%
  inner_join(immune_signature_gene, by="gene_id")

res %>%
  dplyr::select(gene_id, gene_name, log2FoldChange, category) %>%
  kbl(caption='44 IRIS immune markers that are significently expressed in Modereate+ samples relative to the control-like samples.') %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:sig-blood-iris-in-bilat)44 IRIS immune markers that are significently expressed in Modereate+ samples relative to the control-like samples.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> gene_id </th>
   <th style="text-align:left;"> gene_name </th>
   <th style="text-align:right;"> log2FoldChange </th>
   <th style="text-align:left;"> category </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000002933.9 </td>
   <td style="text-align:left;"> TMEM176A </td>
   <td style="text-align:right;"> 1.545149 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000038427.16 </td>
   <td style="text-align:left;"> VCAN </td>
   <td style="text-align:right;"> 1.254348 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000038945.15 </td>
   <td style="text-align:left;"> MSR1 </td>
   <td style="text-align:right;"> 1.631995 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000060558.4 </td>
   <td style="text-align:left;"> GNA15 </td>
   <td style="text-align:right;"> 1.732585 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000078081.8 </td>
   <td style="text-align:left;"> LAMP3 </td>
   <td style="text-align:right;"> 3.208762 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000090104.12 </td>
   <td style="text-align:left;"> RGS1 </td>
   <td style="text-align:right;"> 3.254608 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000090659.18 </td>
   <td style="text-align:left;"> CD209 </td>
   <td style="text-align:right;"> 1.687079 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000100979.15 </td>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 1.539781 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000100985.7 </td>
   <td style="text-align:left;"> MMP9 </td>
   <td style="text-align:right;"> 2.886832 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000104951.16 </td>
   <td style="text-align:left;"> IL4I1 </td>
   <td style="text-align:right;"> 3.070091 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000105967.16 </td>
   <td style="text-align:left;"> TFEC </td>
   <td style="text-align:right;"> 1.733737 </td>
   <td style="text-align:left;"> Bcell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000108846.16 </td>
   <td style="text-align:left;"> ABCC3 </td>
   <td style="text-align:right;"> 2.173910 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000110077.14 </td>
   <td style="text-align:left;"> MS4A6A </td>
   <td style="text-align:right;"> 1.591415 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000110079.18 </td>
   <td style="text-align:left;"> MS4A4A </td>
   <td style="text-align:right;"> 1.647453 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000114013.16 </td>
   <td style="text-align:left;"> CD86 </td>
   <td style="text-align:right;"> 1.609474 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000118785.14 </td>
   <td style="text-align:left;"> SPP1 </td>
   <td style="text-align:right;"> 4.957166 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000124491.16 </td>
   <td style="text-align:left;"> F13A1 </td>
   <td style="text-align:right;"> 1.777443 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000124772.12 </td>
   <td style="text-align:left;"> CPNE5 </td>
   <td style="text-align:right;"> 2.192108 </td>
   <td style="text-align:left;"> Bcell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000125730.17 </td>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 1.689339 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000132205.11 </td>
   <td style="text-align:left;"> EMILIN2 </td>
   <td style="text-align:right;"> 1.440005 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000132514.13 </td>
   <td style="text-align:left;"> CLEC10A </td>
   <td style="text-align:right;"> 1.486597 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000133048.13 </td>
   <td style="text-align:left;"> CHI3L1 </td>
   <td style="text-align:right;"> 2.736904 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000137441.8 </td>
   <td style="text-align:left;"> FGFBP2 </td>
   <td style="text-align:right;"> 2.177403 </td>
   <td style="text-align:left;"> NKcell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000137801.11 </td>
   <td style="text-align:left;"> THBS1 </td>
   <td style="text-align:right;"> 1.995249 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000139572.4 </td>
   <td style="text-align:left;"> GPR84 </td>
   <td style="text-align:right;"> 2.554428 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000143320.9 </td>
   <td style="text-align:left;"> CRABP2 </td>
   <td style="text-align:right;"> 1.922802 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000143387.13 </td>
   <td style="text-align:left;"> CTSK </td>
   <td style="text-align:right;"> 1.456861 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000155659.15 </td>
   <td style="text-align:left;"> VSIG4 </td>
   <td style="text-align:right;"> 1.435979 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000157227.13 </td>
   <td style="text-align:left;"> MMP14 </td>
   <td style="text-align:right;"> 1.302515 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000161921.16 </td>
   <td style="text-align:left;"> CXCL16 </td>
   <td style="text-align:right;"> 1.362402 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000163599.17 </td>
   <td style="text-align:left;"> CTLA4 </td>
   <td style="text-align:right;"> 2.343052 </td>
   <td style="text-align:left;"> CD4Tcell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000166278.15 </td>
   <td style="text-align:left;"> C2 </td>
   <td style="text-align:right;"> 2.486130 </td>
   <td style="text-align:left;"> PlasmaCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000166927.13 </td>
   <td style="text-align:left;"> MS4A7 </td>
   <td style="text-align:right;"> 1.716174 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000169245.6 </td>
   <td style="text-align:left;"> CXCL10 </td>
   <td style="text-align:right;"> 2.637787 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000170458.14 </td>
   <td style="text-align:left;"> CD14 </td>
   <td style="text-align:right;"> 1.377451 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000171812.13 </td>
   <td style="text-align:left;"> COL8A2 </td>
   <td style="text-align:right;"> 1.870228 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000171848.15 </td>
   <td style="text-align:left;"> RRM2 </td>
   <td style="text-align:right;"> 2.263203 </td>
   <td style="text-align:left;"> CD4Tcell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000173369.17 </td>
   <td style="text-align:left;"> C1QB </td>
   <td style="text-align:right;"> 1.939756 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000173372.17 </td>
   <td style="text-align:left;"> C1QA </td>
   <td style="text-align:right;"> 1.684487 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000177575.12 </td>
   <td style="text-align:left;"> CD163 </td>
   <td style="text-align:right;"> 1.736284 </td>
   <td style="text-align:left;"> Monocyte </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000181458.10 </td>
   <td style="text-align:left;"> TMEM45A </td>
   <td style="text-align:right;"> 1.878089 </td>
   <td style="text-align:left;"> PlasmaCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000182578.14 </td>
   <td style="text-align:left;"> CSF1R </td>
   <td style="text-align:right;"> 1.345719 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000187474.5 </td>
   <td style="text-align:left;"> FPR3 </td>
   <td style="text-align:right;"> 1.862982 </td>
   <td style="text-align:left;"> DendriticCell </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000198794.12 </td>
   <td style="text-align:left;"> SCAMP5 </td>
   <td style="text-align:right;"> 1.762450 </td>
   <td style="text-align:left;"> PlasmaCell </td>
  </tr>
</tbody>
</table>

## Complement scoring vs. immune cell infiltrates signature
The earlier mentioned 63 immune cell markers contribute to characterizing the signature of immune cell infiltration. The following code calculates the Spearman correlation coefficient between complement scoring and immune cell signature (average $\log_{10}$TPM of 63 immune cell markers), which yields a value of 0.456.


```r
tpm <- assays(bilat_dds[immune_signature_gene$gene_id])[["TPM"]]

signature <- colMeans(log10(tpm+1)) %>% as.data.frame() %>%
  rownames_to_column(var="sample_id") 
names(signature)[2] <- "avg_tpm"

tmp_df <- comprehensive_df %>% 
  dplyr::select(sample_id, `Complement Scoring`, `STIR_status`,
                `STIR_RATING`) %>%
  left_join(signature, by="sample_id") %>%
  drop_na(`Complement Scoring`, `avg_tpm`) 
  
# Spearman correlation between complement and immune infiltrate signatures
tmp_df %>% 
  drop_na(`Complement Scoring`, `avg_tpm`) %>%
  dplyr::select(`Complement Scoring`, `avg_tpm`) %>%
  corrr::correlate(., method="spearman")
#> # A tibble: 2 × 3
#>   term               `Complement Scoring` avg_tpm
#>   <chr>                             <dbl>   <dbl>
#> 1 Complement Scoring               NA       0.454
#> 2 avg_tpm                           0.454  NA
```



```r
tmp_df %>% 
  drop_na(`Complement Scoring`, `avg_tpm`) %>%
  dplyr::mutate(`Complement Scoring`= factor(`Complement Scoring`)) %>%
  ggplot(aes(x=`Complement Scoring`, y=avg_tpm)) +
    geom_boxplot(width=0.5, outlier.shape = NA) +
    geom_dotplot(aes(fill=`Complement Scoring`),
                 binaxis='y', stackdir='center',
                 show.legend = FALSE, alpha=0.8,
                 stackratio=1.5, dotsize=1.5) +
    theme_minimal()
```

<div class="figure" style="text-align: center">
<img src="08-immune-cell-infiltrates-signature_files/figure-html/63-blood-markers-correlated-with-complement-scoring-1.png" alt="Dotplot of immune cell infiltrates signature (presented by average TPM) groups by complement scoring 1, 2, and 3." width="672" />
<p class="caption">(\#fig:63-blood-markers-correlated-with-complement-scoring)Dotplot of immune cell infiltrates signature (presented by average TPM) groups by complement scoring 1, 2, and 3.</p>
</div>

## STIR status vs. immune cell infiltrates signature
The elevated immune cell signatures is associated with with STIR status: Wilcox test p-value = $2.9e-6$. 


```r
tmp_df %>% 
  drop_na(`STIR_status`, `avg_tpm`) %>%
  ggplot(aes(x=`STIR_status`, y=avg_tpm)) +
    geom_boxplot(width=0.5, outlier.shape = NA) +
    geom_dotplot(aes(fill=`STIR_status`),
                 binaxis='y', stackdir='center',
                 show.legend = FALSE,
                 stackratio=1.5, dotsize=1.5, alpha=0.7) +
    theme_minimal()
```

<img src="08-immune-cell-infiltrates-signature_files/figure-html/STIR-vs-immune-cell-infilrates-signatures-1.png" width="672" style="display: block; margin: auto;" />


```r

# wilcox.test
x <- tmp_df %>% 
  drop_na(`STIR_status`, `avg_tpm`) %>%
  dplyr::select(sample_id, `STIR_status`, `avg_tpm`) %>%
  spread(key=STIR_status, value=avg_tpm) 
wilcox.test(x$`STIR-`, x$`STIR+` )
#> 
#> 	Wilcoxon rank sum exact test
#> 
#> data:  x$`STIR-` and x$`STIR+`
#> W = 134, p-value = 2.902e-06
#> alternative hypothesis: true location shift is not equal to 0
```




## Bilateral comparison analysis for the immune cell infiltrate markers

### 63 significant immune cell markers
In the preceding section, we uncovered a distinctive signature of immune cell infiltrates defined by 63 immune cell markers within the longitudinal study. Now, calculate the expression correlation of these markers between the left and right muscles in the bilateral study.

The code below displays the tools to compute the correlation between L and R muscles on the immune cell infiltrate markers.

```r
col_data <- colData(bilat_dds) %>% as.data.frame()
# convert gene_id to Gencode35
markers <- sig_blood_iris %>%
  left_join(dplyr::select(anno_gencode35, ens_id, gencode35_id), 
            by="ens_id") %>%
  dplyr::select(-gene_id) %>%
  dplyr::rename(gene_id = gencode35_id)

.markers_bilateral_correlation <- function(markers, dds) {
  tpm <- assays(dds[markers$gene_id])[["TPM"]] %>% t(.) %>% 
      as.data.frame() %>%
      rownames_to_column(var="sample_id") %>%
      gather(key=gene_id, value=TPM, -sample_id) %>% 
      left_join(dplyr::select(col_data, sample_id, location,
                              Subject), by="sample_id") %>%
      dplyr::select(-sample_id) %>%
      spread(key=location, value=TPM) %>% 
      drop_na(L, R) 
  gene_cor <- tpm %>% 
    dplyr::mutate(L_log = log10(L+1), R_log=log10(R+1)) %>%
    group_by(gene_id) %>% 
    summarise(Pearson_by_TPM=cor(L, R), 
              Pearson_by_log10TPM = cor(L_log, R_log)) %>%
    dplyr::left_join(markers, by="gene_id") %>%
    arrange(desc(Pearson_by_log10TPM)) %>%
    dplyr::relocate(gene_name, .after=gene_id)
  
  return(gene_cor)
}

.avg_TPM_by_class <- function(markers, dds) {
   avg_TPM <- map_dfc(levels(dds$class), function(x) {
     assays(dds[markers$gene_id])[["TPM"]][, dds$class == x] %>% 
       rowMeans(.) 
   }) %>%
     tibble::add_column(gene_id = markers$gene_id) %>%
     dplyr::rename_with(~paste0(levels(dds$class), " avg TPM"), 
                        starts_with("..."))
}
```

Display the correlation and average TPM grouped by ML-based classes.

```r
gene_cor_63 <- .markers_bilateral_correlation(markers=markers,
                                           dds=bilat_dds) %>%
  dplyr::left_join(.avg_TPM_by_class(markers = markers, dds=bilat_dds), 
                   by="gene_id") %>%
  arrange(category) 

gene_cor_63 %>% 
  dplyr::select(-Pearson_by_TPM) %>%
  kbl(caption="Sixty-three significant immune cell markers correlation analysis between left and right biopsied muscles. The Peason correlation coefficients are compuated using the expression levels presented by Log10(TPM+1).") %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:display-63-sig-immune-cell-tpm)Sixty-three significant immune cell markers correlation analysis between left and right biopsied muscles. The Peason correlation coefficients are compuated using the expression levels presented by Log10(TPM+1).</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> gene_id </th>
   <th style="text-align:left;"> gene_name </th>
   <th style="text-align:right;"> Pearson_by_log10TPM </th>
   <th style="text-align:left;"> ens_id </th>
   <th style="text-align:left;"> category </th>
   <th style="text-align:right;"> Control-like avg TPM </th>
   <th style="text-align:right;"> Moderate+ avg TPM </th>
   <th style="text-align:right;"> Muscle-Low avg TPM </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000124772.12 </td>
   <td style="text-align:left;"> CPNE5 </td>
   <td style="text-align:right;"> 0.8373333 </td>
   <td style="text-align:left;"> ENSG00000124772 </td>
   <td style="text-align:left;"> Bcell </td>
   <td style="text-align:right;"> 0.0443273 </td>
   <td style="text-align:right;"> 0.3835624 </td>
   <td style="text-align:right;"> 1.0004216 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000105967.16 </td>
   <td style="text-align:left;"> TFEC </td>
   <td style="text-align:right;"> 0.8176975 </td>
   <td style="text-align:left;"> ENSG00000105967 </td>
   <td style="text-align:left;"> Bcell </td>
   <td style="text-align:right;"> 0.0234566 </td>
   <td style="text-align:right;"> 0.1235314 </td>
   <td style="text-align:right;"> 0.1023323 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000138180.16 </td>
   <td style="text-align:left;"> CEP55 </td>
   <td style="text-align:right;"> 0.7993619 </td>
   <td style="text-align:left;"> ENSG00000138180 </td>
   <td style="text-align:left;"> CD4Tcell </td>
   <td style="text-align:right;"> 0.0348918 </td>
   <td style="text-align:right;"> 0.1666287 </td>
   <td style="text-align:right;"> 0.2252888 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000163599.17 </td>
   <td style="text-align:left;"> CTLA4 </td>
   <td style="text-align:right;"> 0.7880751 </td>
   <td style="text-align:left;"> ENSG00000163599 </td>
   <td style="text-align:left;"> CD4Tcell </td>
   <td style="text-align:right;"> 0.0315704 </td>
   <td style="text-align:right;"> 0.2644233 </td>
   <td style="text-align:right;"> 0.1854032 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000171848.15 </td>
   <td style="text-align:left;"> RRM2 </td>
   <td style="text-align:right;"> 0.6609791 </td>
   <td style="text-align:left;"> ENSG00000171848 </td>
   <td style="text-align:left;"> CD4Tcell </td>
   <td style="text-align:right;"> 0.0090876 </td>
   <td style="text-align:right;"> 0.0662736 </td>
   <td style="text-align:right;"> 0.1102069 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000174807.4 </td>
   <td style="text-align:left;"> CD248 </td>
   <td style="text-align:right;"> 0.6261514 </td>
   <td style="text-align:left;"> ENSG00000174807 </td>
   <td style="text-align:left;"> CD8Tcell </td>
   <td style="text-align:right;"> 14.3514643 </td>
   <td style="text-align:right;"> 48.2904687 </td>
   <td style="text-align:right;"> 253.0578657 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000114013.16 </td>
   <td style="text-align:left;"> CD86 </td>
   <td style="text-align:right;"> 0.8269231 </td>
   <td style="text-align:left;"> ENSG00000114013 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.0663682 </td>
   <td style="text-align:right;"> 0.3421849 </td>
   <td style="text-align:right;"> 0.4462872 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000132514.13 </td>
   <td style="text-align:left;"> CLEC10A </td>
   <td style="text-align:right;"> 0.8181479 </td>
   <td style="text-align:left;"> ENSG00000132514 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.3866177 </td>
   <td style="text-align:right;"> 1.7920233 </td>
   <td style="text-align:right;"> 3.3907985 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000173369.17 </td>
   <td style="text-align:left;"> C1QB </td>
   <td style="text-align:right;"> 0.8099608 </td>
   <td style="text-align:left;"> ENSG00000173369 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 2.0357169 </td>
   <td style="text-align:right;"> 13.7189810 </td>
   <td style="text-align:right;"> 28.4657562 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000110077.14 </td>
   <td style="text-align:left;"> MS4A6A </td>
   <td style="text-align:right;"> 0.8063214 </td>
   <td style="text-align:left;"> ENSG00000110077 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.4436241 </td>
   <td style="text-align:right;"> 2.1606915 </td>
   <td style="text-align:right;"> 3.6416739 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000166927.13 </td>
   <td style="text-align:left;"> MS4A7 </td>
   <td style="text-align:right;"> 0.7842506 </td>
   <td style="text-align:left;"> ENSG00000166927 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.2090296 </td>
   <td style="text-align:right;"> 1.1084212 </td>
   <td style="text-align:right;"> 1.9491638 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000124491.16 </td>
   <td style="text-align:left;"> F13A1 </td>
   <td style="text-align:right;"> 0.7788406 </td>
   <td style="text-align:left;"> ENSG00000124491 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 2.0078292 </td>
   <td style="text-align:right;"> 12.2302448 </td>
   <td style="text-align:right;"> 21.8912499 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000182578.14 </td>
   <td style="text-align:left;"> CSF1R </td>
   <td style="text-align:right;"> 0.7750512 </td>
   <td style="text-align:left;"> ENSG00000182578 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.7874655 </td>
   <td style="text-align:right;"> 3.3187330 </td>
   <td style="text-align:right;"> 7.1733571 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000282608.2 </td>
   <td style="text-align:left;"> ADORA3 </td>
   <td style="text-align:right;"> 0.7736084 </td>
   <td style="text-align:left;"> ENSG00000282608 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.0386997 </td>
   <td style="text-align:right;"> 0.1584448 </td>
   <td style="text-align:right;"> 0.2628639 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000173372.17 </td>
   <td style="text-align:left;"> C1QA </td>
   <td style="text-align:right;"> 0.7445602 </td>
   <td style="text-align:left;"> ENSG00000173372 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 3.4317578 </td>
   <td style="text-align:right;"> 19.1178534 </td>
   <td style="text-align:right;"> 52.6300218 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000161921.16 </td>
   <td style="text-align:left;"> CXCL16 </td>
   <td style="text-align:right;"> 0.7398670 </td>
   <td style="text-align:left;"> ENSG00000161921 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.3268772 </td>
   <td style="text-align:right;"> 1.4080810 </td>
   <td style="text-align:right;"> 5.0737114 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000010327.10 </td>
   <td style="text-align:left;"> STAB1 </td>
   <td style="text-align:right;"> 0.7266172 </td>
   <td style="text-align:left;"> ENSG00000010327 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 1.5300961 </td>
   <td style="text-align:right;"> 5.3923377 </td>
   <td style="text-align:right;"> 17.3660509 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000078081.8 </td>
   <td style="text-align:left;"> LAMP3 </td>
   <td style="text-align:right;"> 0.7185273 </td>
   <td style="text-align:left;"> ENSG00000078081 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.0108336 </td>
   <td style="text-align:right;"> 0.1663221 </td>
   <td style="text-align:right;"> 0.2174447 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000090104.12 </td>
   <td style="text-align:left;"> RGS1 </td>
   <td style="text-align:right;"> 0.7173044 </td>
   <td style="text-align:left;"> ENSG00000090104 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.0130089 </td>
   <td style="text-align:right;"> 0.2018628 </td>
   <td style="text-align:right;"> 0.1915478 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000184060.11 </td>
   <td style="text-align:left;"> ADAP2 </td>
   <td style="text-align:right;"> 0.7103740 </td>
   <td style="text-align:left;"> ENSG00000184060 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.3838474 </td>
   <td style="text-align:right;"> 1.1441190 </td>
   <td style="text-align:right;"> 3.2110364 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000104951.16 </td>
   <td style="text-align:left;"> IL4I1 </td>
   <td style="text-align:right;"> 0.6975342 </td>
   <td style="text-align:left;"> ENSG00000104951 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.0220052 </td>
   <td style="text-align:right;"> 0.3013386 </td>
   <td style="text-align:right;"> 0.1689496 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000120708.17 </td>
   <td style="text-align:left;"> TGFBI </td>
   <td style="text-align:right;"> 0.6748693 </td>
   <td style="text-align:left;"> ENSG00000120708 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 2.4659110 </td>
   <td style="text-align:right;"> 8.9473151 </td>
   <td style="text-align:right;"> 33.0267668 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000002933.9 </td>
   <td style="text-align:left;"> TMEM176A </td>
   <td style="text-align:right;"> 0.6728172 </td>
   <td style="text-align:left;"> ENSG00000002933 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.3762751 </td>
   <td style="text-align:right;"> 1.8857893 </td>
   <td style="text-align:right;"> 11.2110257 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000169245.6 </td>
   <td style="text-align:left;"> CXCL10 </td>
   <td style="text-align:right;"> 0.6550995 </td>
   <td style="text-align:left;"> ENSG00000169245 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.5894382 </td>
   <td style="text-align:right;"> 5.9803700 </td>
   <td style="text-align:right;"> 2.1547531 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000137801.11 </td>
   <td style="text-align:left;"> THBS1 </td>
   <td style="text-align:right;"> 0.5769849 </td>
   <td style="text-align:left;"> ENSG00000137801 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.8099008 </td>
   <td style="text-align:right;"> 5.0986959 </td>
   <td style="text-align:right;"> 17.6310690 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000132205.11 </td>
   <td style="text-align:left;"> EMILIN2 </td>
   <td style="text-align:right;"> 0.5297056 </td>
   <td style="text-align:left;"> ENSG00000132205 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.3997817 </td>
   <td style="text-align:right;"> 1.9355804 </td>
   <td style="text-align:right;"> 10.0447402 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000108846.16 </td>
   <td style="text-align:left;"> ABCC3 </td>
   <td style="text-align:right;"> 0.5169606 </td>
   <td style="text-align:left;"> ENSG00000108846 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.0253774 </td>
   <td style="text-align:right;"> 0.2010368 </td>
   <td style="text-align:right;"> 0.4090524 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000088827.12 </td>
   <td style="text-align:left;"> SIGLEC1 </td>
   <td style="text-align:right;"> 0.5044042 </td>
   <td style="text-align:left;"> ENSG00000088827 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 1.0464540 </td>
   <td style="text-align:right;"> 4.0403822 </td>
   <td style="text-align:right;"> 10.0460345 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000166145.14 </td>
   <td style="text-align:left;"> SPINT1 </td>
   <td style="text-align:right;"> 0.5003076 </td>
   <td style="text-align:left;"> ENSG00000166145 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.0426139 </td>
   <td style="text-align:right;"> 0.1894087 </td>
   <td style="text-align:right;"> 0.7413172 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000187474.5 </td>
   <td style="text-align:left;"> FPR3 </td>
   <td style="text-align:right;"> 0.4933989 </td>
   <td style="text-align:left;"> ENSG00000187474 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.3026753 </td>
   <td style="text-align:right;"> 1.7773001 </td>
   <td style="text-align:right;"> 1.9219966 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000090659.18 </td>
   <td style="text-align:left;"> CD209 </td>
   <td style="text-align:right;"> 0.4608114 </td>
   <td style="text-align:left;"> ENSG00000090659 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.2322114 </td>
   <td style="text-align:right;"> 1.2533845 </td>
   <td style="text-align:right;"> 3.3028487 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000176014.13 </td>
   <td style="text-align:left;"> TUBB6 </td>
   <td style="text-align:right;"> 0.4602100 </td>
   <td style="text-align:left;"> ENSG00000176014 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 3.5972900 </td>
   <td style="text-align:right;"> 11.0828601 </td>
   <td style="text-align:right;"> 48.3484150 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000118785.14 </td>
   <td style="text-align:left;"> SPP1 </td>
   <td style="text-align:right;"> 0.4479528 </td>
   <td style="text-align:left;"> ENSG00000118785 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.0207772 </td>
   <td style="text-align:right;"> 1.1074521 </td>
   <td style="text-align:right;"> 1.8615908 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000011422.12 </td>
   <td style="text-align:left;"> PLAUR </td>
   <td style="text-align:right;"> 0.4301508 </td>
   <td style="text-align:left;"> ENSG00000011422 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 0.1783956 </td>
   <td style="text-align:right;"> 0.5793819 </td>
   <td style="text-align:right;"> 2.2774659 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000141753.7 </td>
   <td style="text-align:left;"> IGFBP4 </td>
   <td style="text-align:right;"> 0.3779018 </td>
   <td style="text-align:left;"> ENSG00000141753 </td>
   <td style="text-align:left;"> DendriticCell </td>
   <td style="text-align:right;"> 47.2044767 </td>
   <td style="text-align:right;"> 167.3634231 </td>
   <td style="text-align:right;"> 1281.8420877 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000038945.15 </td>
   <td style="text-align:left;"> MSR1 </td>
   <td style="text-align:right;"> 0.8538713 </td>
   <td style="text-align:left;"> ENSG00000038945 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.1165807 </td>
   <td style="text-align:right;"> 0.5913336 </td>
   <td style="text-align:right;"> 0.8065051 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000177575.12 </td>
   <td style="text-align:left;"> CD163 </td>
   <td style="text-align:right;"> 0.8078951 </td>
   <td style="text-align:left;"> ENSG00000177575 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.3361453 </td>
   <td style="text-align:right;"> 1.8946880 </td>
   <td style="text-align:right;"> 2.0304719 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000170458.14 </td>
   <td style="text-align:left;"> CD14 </td>
   <td style="text-align:right;"> 0.7842175 </td>
   <td style="text-align:left;"> ENSG00000170458 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 1.7751882 </td>
   <td style="text-align:right;"> 7.7468300 </td>
   <td style="text-align:right;"> 21.4222643 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000161944.16 </td>
   <td style="text-align:left;"> ASGR2 </td>
   <td style="text-align:right;"> 0.7740608 </td>
   <td style="text-align:left;"> ENSG00000161944 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.0212327 </td>
   <td style="text-align:right;"> 0.1023889 </td>
   <td style="text-align:right;"> 0.1496403 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000110079.18 </td>
   <td style="text-align:left;"> MS4A4A </td>
   <td style="text-align:right;"> 0.7206380 </td>
   <td style="text-align:left;"> ENSG00000110079 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.1930113 </td>
   <td style="text-align:right;"> 1.0180924 </td>
   <td style="text-align:right;"> 2.7351409 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000143320.9 </td>
   <td style="text-align:left;"> CRABP2 </td>
   <td style="text-align:right;"> 0.7144369 </td>
   <td style="text-align:left;"> ENSG00000143320 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.6114227 </td>
   <td style="text-align:right;"> 3.8166488 </td>
   <td style="text-align:right;"> 13.9456375 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000019169.10 </td>
   <td style="text-align:left;"> MARCO </td>
   <td style="text-align:right;"> 0.7105930 </td>
   <td style="text-align:left;"> ENSG00000019169 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.1199304 </td>
   <td style="text-align:right;"> 0.6352074 </td>
   <td style="text-align:right;"> 1.9900735 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000157227.13 </td>
   <td style="text-align:left;"> MMP14 </td>
   <td style="text-align:right;"> 0.6813758 </td>
   <td style="text-align:left;"> ENSG00000157227 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 2.4939224 </td>
   <td style="text-align:right;"> 10.9206905 </td>
   <td style="text-align:right;"> 68.0308145 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000143387.13 </td>
   <td style="text-align:left;"> CTSK </td>
   <td style="text-align:right;"> 0.6812790 </td>
   <td style="text-align:left;"> ENSG00000143387 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 2.5631582 </td>
   <td style="text-align:right;"> 12.3782990 </td>
   <td style="text-align:right;"> 70.4075115 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000038427.16 </td>
   <td style="text-align:left;"> VCAN </td>
   <td style="text-align:right;"> 0.6795530 </td>
   <td style="text-align:left;"> ENSG00000038427 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.5364233 </td>
   <td style="text-align:right;"> 2.1562978 </td>
   <td style="text-align:right;"> 3.3628331 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000133048.13 </td>
   <td style="text-align:left;"> CHI3L1 </td>
   <td style="text-align:right;"> 0.6346830 </td>
   <td style="text-align:left;"> ENSG00000133048 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.1946220 </td>
   <td style="text-align:right;"> 2.1621032 </td>
   <td style="text-align:right;"> 2.7397122 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000165140.11 </td>
   <td style="text-align:left;"> FBP1 </td>
   <td style="text-align:right;"> 0.6256957 </td>
   <td style="text-align:left;"> ENSG00000165140 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.1543315 </td>
   <td style="text-align:right;"> 0.6023780 </td>
   <td style="text-align:right;"> 3.7040663 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000141505.12 </td>
   <td style="text-align:left;"> ASGR1 </td>
   <td style="text-align:right;"> 0.5509758 </td>
   <td style="text-align:left;"> ENSG00000141505 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.0351094 </td>
   <td style="text-align:right;"> 0.1354246 </td>
   <td style="text-align:right;"> 0.6025999 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000171812.13 </td>
   <td style="text-align:left;"> COL8A2 </td>
   <td style="text-align:right;"> 0.5351123 </td>
   <td style="text-align:left;"> ENSG00000171812 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.1500623 </td>
   <td style="text-align:right;"> 0.9979313 </td>
   <td style="text-align:right;"> 3.0772092 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000060558.4 </td>
   <td style="text-align:left;"> GNA15 </td>
   <td style="text-align:right;"> 0.5164559 </td>
   <td style="text-align:left;"> ENSG00000060558 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.1154703 </td>
   <td style="text-align:right;"> 0.6176169 </td>
   <td style="text-align:right;"> 1.4995404 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000155659.15 </td>
   <td style="text-align:left;"> VSIG4 </td>
   <td style="text-align:right;"> 0.5133328 </td>
   <td style="text-align:left;"> ENSG00000155659 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.5261769 </td>
   <td style="text-align:right;"> 2.2279206 </td>
   <td style="text-align:right;"> 4.3380216 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000100979.15 </td>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 0.5012153 </td>
   <td style="text-align:left;"> ENSG00000100979 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 1.8178334 </td>
   <td style="text-align:right;"> 9.1650567 </td>
   <td style="text-align:right;"> 44.0062775 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000100985.7 </td>
   <td style="text-align:left;"> MMP9 </td>
   <td style="text-align:right;"> 0.4984437 </td>
   <td style="text-align:left;"> ENSG00000100985 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.1003691 </td>
   <td style="text-align:right;"> 1.5203490 </td>
   <td style="text-align:right;"> 8.2657241 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000139572.4 </td>
   <td style="text-align:left;"> GPR84 </td>
   <td style="text-align:right;"> 0.4751397 </td>
   <td style="text-align:left;"> ENSG00000139572 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.0077004 </td>
   <td style="text-align:right;"> 0.0773290 </td>
   <td style="text-align:right;"> 0.0618142 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000125730.17 </td>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.4263273 </td>
   <td style="text-align:left;"> ENSG00000125730 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 4.4394449 </td>
   <td style="text-align:right;"> 25.6278376 </td>
   <td style="text-align:right;"> 67.5924876 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000134668.12 </td>
   <td style="text-align:left;"> SPOCD1 </td>
   <td style="text-align:right;"> 0.3293381 </td>
   <td style="text-align:left;"> ENSG00000134668 </td>
   <td style="text-align:left;"> Monocyte </td>
   <td style="text-align:right;"> 0.0101972 </td>
   <td style="text-align:right;"> 0.0482261 </td>
   <td style="text-align:right;"> 0.3785615 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000100365.16 </td>
   <td style="text-align:left;"> NCF4 </td>
   <td style="text-align:right;"> 0.5626164 </td>
   <td style="text-align:left;"> ENSG00000100365 </td>
   <td style="text-align:left;"> Neutrophil </td>
   <td style="text-align:right;"> 0.3574015 </td>
   <td style="text-align:right;"> 1.3059616 </td>
   <td style="text-align:right;"> 3.1460980 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000137441.8 </td>
   <td style="text-align:left;"> FGFBP2 </td>
   <td style="text-align:right;"> 0.7745707 </td>
   <td style="text-align:left;"> ENSG00000137441 </td>
   <td style="text-align:left;"> NKcell </td>
   <td style="text-align:right;"> 0.4225620 </td>
   <td style="text-align:right;"> 3.5214963 </td>
   <td style="text-align:right;"> 12.7187471 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000166278.15 </td>
   <td style="text-align:left;"> C2 </td>
   <td style="text-align:right;"> 0.7007918 </td>
   <td style="text-align:left;"> ENSG00000166278 </td>
   <td style="text-align:left;"> PlasmaCell </td>
   <td style="text-align:right;"> 0.0039973 </td>
   <td style="text-align:right;"> 0.0403557 </td>
   <td style="text-align:right;"> 0.1117194 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000110492.15 </td>
   <td style="text-align:left;"> MDK </td>
   <td style="text-align:right;"> 0.7003678 </td>
   <td style="text-align:left;"> ENSG00000110492 </td>
   <td style="text-align:left;"> PlasmaCell </td>
   <td style="text-align:right;"> 0.5089112 </td>
   <td style="text-align:right;"> 2.5043028 </td>
   <td style="text-align:right;"> 13.2510295 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000181458.10 </td>
   <td style="text-align:left;"> TMEM45A </td>
   <td style="text-align:right;"> 0.7003172 </td>
   <td style="text-align:left;"> ENSG00000181458 </td>
   <td style="text-align:left;"> PlasmaCell </td>
   <td style="text-align:right;"> 0.0529060 </td>
   <td style="text-align:right;"> 0.3722313 </td>
   <td style="text-align:right;"> 1.0593557 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000198794.12 </td>
   <td style="text-align:left;"> SCAMP5 </td>
   <td style="text-align:right;"> 0.6259889 </td>
   <td style="text-align:left;"> ENSG00000198794 </td>
   <td style="text-align:left;"> PlasmaCell </td>
   <td style="text-align:right;"> 0.0735910 </td>
   <td style="text-align:right;"> 0.4475840 </td>
   <td style="text-align:right;"> 2.3682200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000142552.8 </td>
   <td style="text-align:left;"> RCN3 </td>
   <td style="text-align:right;"> 0.5367903 </td>
   <td style="text-align:left;"> ENSG00000142552 </td>
   <td style="text-align:left;"> PlasmaCell </td>
   <td style="text-align:right;"> 3.1691866 </td>
   <td style="text-align:right;"> 10.1816596 </td>
   <td style="text-align:right;"> 48.4596901 </td>
  </tr>
</tbody>
</table>

__Scatter plot__

```r
x <- gene_cor_63 %>%
  dplyr::filter(category %in% c("Bcell", "PlasmaCell"))

tpm <- assays(bilat_dds[x$gene_id])[["TPM"]] %>% t(.) %>% 
      as.data.frame() %>%
      rownames_to_column(var="sample_id") %>%
      gather(key=gene_id, value=TPM, -sample_id) %>% 
      left_join(dplyr::select(col_data, sample_id, location,
                              Subject), by="sample_id") %>%
      dplyr::select(-sample_id) %>%
      spread(key=location, value=TPM) %>% 
      drop_na(L, R) %>%
      dplyr::left_join(x, by="gene_id") %>%
      dplyr::mutate(log10L = log10(L+1)) %>%
      dplyr::mutate(log10R = log10(R+1))

ggplot(tpm, aes(x=log10L, y=log10R)) +
    geom_point(size=1, alpha=0.7) +
    geom_smooth(method="lm", linewidth=0.5, se=FALSE, alpha=0.5) +
    facet_wrap(~gene_name, nrow = 2, scales="free") +
    theme_classic() +
    theme(plot.title=element_text(hjust=0.5, size=10)) +
    labs(x="L (log10(TPM+1))", y="R (log10(TPM+1))")
```

<img src="08-immune-cell-infiltrates-signature_files/figure-html/bcell-plasma-markers-scatter-plot-TPM-1.png" width="672" />

```r
ggsave(file.path(pkg_dir, "figures", "immune-infiltration",
                 "bilat-sig-Bcell-Plasma-in-Longitudinal.pdf"),
       width=6, height=3)
```

### Top loading variables by the PLIER package

In the previous chapter, we employed The PLIER package to estimate the contributions of immune cell types to the FSHD transcriptome. The top Z-loading variables play a crucial role in determining the relative proportions of these immune cell type contributions. Below, we list these top loading variables and their bilateral correlation between L and R biopsies. Please note that they may not necessarily exhibit strong expression levels.


```r
load(file.path(pkg_dir, "data", "plier_topZ_markers.rda"))
markers <- plier_topZ_markers$bilat_markers
gene_cor_topZ <- .markers_bilateral_correlation(markers=markers,
                                           dds=bilat_dds) %>%
    dplyr::left_join(.avg_TPM_by_class(markers = markers, dds=bilat_dds), 
                   by="gene_id")

gene_cor_topZ %>% 
  dplyr::select(-Pearson_by_TPM) %>%
  kbl(caption="Identified by the PLIER package, listed below are the 18 top Z markers contributed to cell type proportions in bilateral study. The Peason correlation coefficients are compuated using the expression levels presented by Log10(TPM+1).") %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:plier-topz-markers-bilateral)Identified by the PLIER package, listed below are the 18 top Z markers contributed to cell type proportions in bilateral study. The Peason correlation coefficients are compuated using the expression levels presented by Log10(TPM+1).</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> gene_id </th>
   <th style="text-align:left;"> gene_name </th>
   <th style="text-align:right;"> Pearson_by_log10TPM </th>
   <th style="text-align:left;"> LV_index </th>
   <th style="text-align:left;"> cell_type </th>
   <th style="text-align:left;"> ens_id </th>
   <th style="text-align:right;"> Control-like avg TPM </th>
   <th style="text-align:right;"> Moderate+ avg TPM </th>
   <th style="text-align:right;"> Muscle-Low avg TPM </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000110777.12 </td>
   <td style="text-align:left;"> POU2AF1 </td>
   <td style="text-align:right;"> 0.9220746 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000110777 </td>
   <td style="text-align:right;"> 0.0312330 </td>
   <td style="text-align:right;"> 0.3724736 </td>
   <td style="text-align:right;"> 0.1057589 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000136573.14 </td>
   <td style="text-align:left;"> BLK </td>
   <td style="text-align:right;"> 0.8141694 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000136573 </td>
   <td style="text-align:right;"> 0.0108037 </td>
   <td style="text-align:right;"> 0.1151240 </td>
   <td style="text-align:right;"> 0.0270017 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000110077.14 </td>
   <td style="text-align:left;"> MS4A6A </td>
   <td style="text-align:right;"> 0.8063214 </td>
   <td style="text-align:left;"> 20,IRIS_DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> ENSG00000110077 </td>
   <td style="text-align:right;"> 0.4436241 </td>
   <td style="text-align:right;"> 2.1606915 </td>
   <td style="text-align:right;"> 3.6416739 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000177455.15 </td>
   <td style="text-align:left;"> CD19 </td>
   <td style="text-align:right;"> 0.7963163 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000177455 </td>
   <td style="text-align:right;"> 0.0199848 </td>
   <td style="text-align:right;"> 0.1477995 </td>
   <td style="text-align:right;"> 0.0461397 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000166927.13 </td>
   <td style="text-align:left;"> MS4A7 </td>
   <td style="text-align:right;"> 0.7842506 </td>
   <td style="text-align:left;"> 20,IRIS_DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> ENSG00000166927 </td>
   <td style="text-align:right;"> 0.2090296 </td>
   <td style="text-align:right;"> 1.1084212 </td>
   <td style="text-align:right;"> 1.9491638 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000132704.16 </td>
   <td style="text-align:left;"> FCRL2 </td>
   <td style="text-align:right;"> 0.7786763 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000132704 </td>
   <td style="text-align:right;"> 0.0216208 </td>
   <td style="text-align:right;"> 0.1253008 </td>
   <td style="text-align:right;"> 0.0336529 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000282608.2 </td>
   <td style="text-align:left;"> ADORA3 </td>
   <td style="text-align:right;"> 0.7736084 </td>
   <td style="text-align:left;"> 20,IRIS_DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> ENSG00000282608 </td>
   <td style="text-align:right;"> 0.0386997 </td>
   <td style="text-align:right;"> 0.1584448 </td>
   <td style="text-align:right;"> 0.2628639 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000156738.18 </td>
   <td style="text-align:left;"> MS4A1 </td>
   <td style="text-align:right;"> 0.6979562 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000156738 </td>
   <td style="text-align:right;"> 0.0939018 </td>
   <td style="text-align:right;"> 0.5551574 </td>
   <td style="text-align:right;"> 0.1685982 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000119866.22 </td>
   <td style="text-align:left;"> BCL11A </td>
   <td style="text-align:right;"> 0.6909547 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000119866 </td>
   <td style="text-align:right;"> 0.0025620 </td>
   <td style="text-align:right;"> 0.0148748 </td>
   <td style="text-align:right;"> 0.0093517 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000174123.11 </td>
   <td style="text-align:left;"> TLR10 </td>
   <td style="text-align:right;"> 0.6724404 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000174123 </td>
   <td style="text-align:right;"> 0.0186041 </td>
   <td style="text-align:right;"> 0.0667165 </td>
   <td style="text-align:right;"> 0.0748956 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000147454.14 </td>
   <td style="text-align:left;"> SLC25A37 </td>
   <td style="text-align:right;"> 0.6618809 </td>
   <td style="text-align:left;"> 22,IRIS_Neutrophil-Resting </td>
   <td style="text-align:left;"> Neutrophil-Resting </td>
   <td style="text-align:left;"> ENSG00000147454 </td>
   <td style="text-align:right;"> 6.7495459 </td>
   <td style="text-align:right;"> 9.1052736 </td>
   <td style="text-align:right;"> 12.1514379 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000100473.18 </td>
   <td style="text-align:left;"> COCH </td>
   <td style="text-align:right;"> 0.6310987 </td>
   <td style="text-align:left;"> 22,IRIS_Neutrophil-Resting </td>
   <td style="text-align:left;"> Neutrophil-Resting </td>
   <td style="text-align:left;"> ENSG00000100473 </td>
   <td style="text-align:right;"> 0.1073993 </td>
   <td style="text-align:right;"> 0.4540094 </td>
   <td style="text-align:right;"> 0.9624094 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000128218.8 </td>
   <td style="text-align:left;"> VPREB3 </td>
   <td style="text-align:right;"> 0.5787516 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000128218 </td>
   <td style="text-align:right;"> 0.3698683 </td>
   <td style="text-align:right;"> 0.8011012 </td>
   <td style="text-align:right;"> 0.5077516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000117322.18 </td>
   <td style="text-align:left;"> CR2 </td>
   <td style="text-align:right;"> 0.5075859 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000117322 </td>
   <td style="text-align:right;"> 0.0166510 </td>
   <td style="text-align:right;"> 0.0830997 </td>
   <td style="text-align:right;"> 0.0050237 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000187474.5 </td>
   <td style="text-align:left;"> FPR3 </td>
   <td style="text-align:right;"> 0.4933989 </td>
   <td style="text-align:left;"> 20,IRIS_DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> ENSG00000187474 </td>
   <td style="text-align:right;"> 0.3026753 </td>
   <td style="text-align:right;"> 1.7773001 </td>
   <td style="text-align:right;"> 1.9219966 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000158477.7 </td>
   <td style="text-align:left;"> CD1A </td>
   <td style="text-align:right;"> 0.4457300 </td>
   <td style="text-align:left;"> 20,IRIS_DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> DendriticCell-LPSstimulated </td>
   <td style="text-align:left;"> ENSG00000158477 </td>
   <td style="text-align:right;"> 0.0171722 </td>
   <td style="text-align:right;"> 0.0480016 </td>
   <td style="text-align:right;"> 0.1425001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000244734.4 </td>
   <td style="text-align:left;"> HBB </td>
   <td style="text-align:right;"> 0.4277861 </td>
   <td style="text-align:left;"> 22,IRIS_Neutrophil-Resting </td>
   <td style="text-align:left;"> Neutrophil-Resting </td>
   <td style="text-align:left;"> ENSG00000244734 </td>
   <td style="text-align:right;"> 468.3524397 </td>
   <td style="text-align:right;"> 1290.3398086 </td>
   <td style="text-align:right;"> 720.0911969 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000163534.15 </td>
   <td style="text-align:left;"> FCRL1 </td>
   <td style="text-align:right;"> 0.3398749 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000163534 </td>
   <td style="text-align:right;"> 0.0200861 </td>
   <td style="text-align:right;"> 0.0731329 </td>
   <td style="text-align:right;"> 0.0216181 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000143546.10 </td>
   <td style="text-align:left;"> S100A8 </td>
   <td style="text-align:right;"> 0.3275266 </td>
   <td style="text-align:left;"> 22,IRIS_Neutrophil-Resting </td>
   <td style="text-align:left;"> Neutrophil-Resting </td>
   <td style="text-align:left;"> ENSG00000143546 </td>
   <td style="text-align:right;"> 4.3832407 </td>
   <td style="text-align:right;"> 7.8536112 </td>
   <td style="text-align:right;"> 3.5999534 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000196092.14 </td>
   <td style="text-align:left;"> PAX5 </td>
   <td style="text-align:right;"> 0.3205516 </td>
   <td style="text-align:left;"> 2,IRIS_Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> Bcell-Memory_IgG_IgA </td>
   <td style="text-align:left;"> ENSG00000196092 </td>
   <td style="text-align:right;"> 0.0535907 </td>
   <td style="text-align:right;"> 0.0706754 </td>
   <td style="text-align:right;"> 0.0370464 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000140932.10 </td>
   <td style="text-align:left;"> CMTM2 </td>
   <td style="text-align:right;"> 0.2400213 </td>
   <td style="text-align:left;"> 22,IRIS_Neutrophil-Resting </td>
   <td style="text-align:left;"> Neutrophil-Resting </td>
   <td style="text-align:left;"> ENSG00000140932 </td>
   <td style="text-align:right;"> 0.1152954 </td>
   <td style="text-align:right;"> 0.2188711 </td>
   <td style="text-align:right;"> 0.2500067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000163221.9 </td>
   <td style="text-align:left;"> S100A12 </td>
   <td style="text-align:right;"> 0.2350988 </td>
   <td style="text-align:left;"> 22,IRIS_Neutrophil-Resting </td>
   <td style="text-align:left;"> Neutrophil-Resting </td>
   <td style="text-align:left;"> ENSG00000163221 </td>
   <td style="text-align:right;"> 0.6872278 </td>
   <td style="text-align:right;"> 1.6299728 </td>
   <td style="text-align:right;"> 0.6324510 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000100721.11 </td>
   <td style="text-align:left;"> TCL1A </td>
   <td style="text-align:right;"> 0.2112394 </td>
   <td style="text-align:left;"> 22,IRIS_Neutrophil-Resting </td>
   <td style="text-align:left;"> Neutrophil-Resting </td>
   <td style="text-align:left;"> ENSG00000100721 </td>
   <td style="text-align:right;"> 0.0296388 </td>
   <td style="text-align:right;"> 0.0680252 </td>
   <td style="text-align:right;"> 0.0049452 </td>
  </tr>
</tbody>
</table>


```r
write_xlsx(list(`63 immune-cell-markers`= gene_cor_63,
                `23 top loading variables (PLIER)`=gene_cor_topZ),
           path=file.path(pkg_dir, "stats",
          "Bilateral-correlation-63-immune-cell-markers-and-23-topZ-by-PLIER.xlsx")) 
```

Check the scatter plot: TPM on x and y coordinates are displayed in log10 scale.


<img src="08-immune-cell-infiltrates-signature_files/figure-html/some-sig-markers-bcell-1.png" width="672" />

 
