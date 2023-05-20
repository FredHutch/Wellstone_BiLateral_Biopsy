# viz tools

per_gene_boxplot_groupby_cluster <- 
  function(dds, gene_info, factor, title=NULL,
           subtitle=NULL) {
    # gene_info: data.frame with columns gene_id and gene_name 
    # start with TPM
    sub <- dds[gene_info$gene_id, ]
    col_data <- colData(sub) %>% as.data.frame() %>%
      dplyr::select(sample_id, factor)

    data <- assays(sub)[["TPM"]] %>% t(.) %>%
      as.data.frame() %>%
      rownames_to_column(var="sample_id") %>%
      gather(key=gene_id, value=TPM, -sample_id) %>%
      left_join(gene_info,
                by="gene_id") %>%
      left_join(col_data, by="sample_id") %>%
      dplyr::mutate(log10TPM = log10(TPM+1))

  ## boxplot
  ggplot(data, aes(x=get(factor), y=TPM)) +
    geom_boxplot(width=0.5, outlier.shape=NA) +
    geom_jitter(width=0.3, size=0.3, alpha=0.5, color="grey50")+
    facet_wrap(~gene_name, scale="free_y") +
    theme_classic() +
    #scale_y_continuous(trans='log10') +
    labs(title=title, subtitle=subtitle) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x= element_blank())
}

# show the re
topZ_heatmp_group_by_types_class <- 
  # row gaps based on cell types
  # column gaps based on cluster
  function(markers_topZ, dds, factor, 
           annotation_col, ann_cor, ...) {
    # markers_topZ is either bilat_markers or longi_markers including the columns: gene_id, gene_name, cell_type
    class <- colData(dds)[, factor]
    tpm <- lapply(levels(class), function(x) {
      s <- colnames(dds)[class == x]
      data <- assays(dds[markers_topZ$gene_id, s])[["TPM"]] 
      data <- log10(data+1)
      # re-arrange by the colmeans (look heatmap)
      Colv <- colMeans(data, na.rm = TRUE)
      data <- data[, order(Colv)]
    })

    tpm <- do.call(cbind, tpm)
    rownames(tpm) <- markers_topZ$gene_name
    zscore_tpm <- (tpm - rowMeans(tpm)) / rowSds(tpm)
    gaps_col <- cumsum(table(class))
    
    #annotation_col <- colData(dds)[colnames(tpm), factor, drop=FALSE] %>% as.data.frame()
    #colnames(annotation_col) <- "class"
    
    # rows
    annotation_row <- markers_topZ[, "cell_type", drop=FALSE]
    rownames(annotation_row) <- markers_topZ$gene_name
    gaps_row <- cumsum(table(markers_topZ$cell_type))
   
    breaks <- c(seq(-2, 2, length.out=98), c(3,4))
    pheatmap(zscore_tpm,
             breaks = breaks,
             gaps_col = gaps_col,
             gaps_row = gaps_row,
             show_colnames = FALSE,
             fontsize_row=7,
             cellheight=8,
             annotation_col = annotation_col,
             annotation_row = annotation_row,
             scale="none",
             annotation_colors=ann_cor,
             angle_col = 90,
             cluster_rows=FALSE, 
             cluster_cols=FALSE, ...)
}

markers_heatmp_group_by_types_class <- 
  # markers: data.frame with gene_id and gene_name
  # row gaps based on cell types
  # column gaps based on cluster
  function(markers, dds, factor, #file_name,
           annotation_col, ann_cor, ...) {
    class <- colData(dds)[, factor]
    tpm <- lapply(levels(class), function(x) {
      s <- colnames(dds)[class == x]
      data <- assays(dds[markers$gene_id, s])[["TPM"]] 
      data <- log10(data+1)
      # re-arrange by the colmeans (look heatmap)
      Colv <- colMeans(data, na.rm = TRUE)
      data <- data[, order(Colv)]
    })
    tpm <- do.call(cbind, tpm)
    rownames(tpm) <- markers$gene_name
    zscore_tpm <- (tpm - rowMeans(tpm)) / rowSds(tpm)
    gaps_col <- cumsum(table(class))

    # column annotation
    annotation_col <- annotation_col[colnames(tpm), ]

    # rows
    annotation_row <- markers[, "category", drop=FALSE] %>% as.data.frame()
    rownames(annotation_row) <- markers$gene_name
    #gaps_row <- cumsum(table(markers_topZ$cell_type))
   
    breaks <- c(seq(-2, 2, length.out=98), 3, 4)
    pheatmap(zscore_tpm,
             breaks = breaks,
             gaps_col = gaps_col,
             #gaps_row = gaps_row,
             show_colnames = FALSE,
             fontsize_row=6,
             cellheight=6.5,
             annotation_col = annotation_col,
             annotation_row = annotation_row,
             scale="none",
             annotation_colors=ann_cor,
             angle_col = 90,
             #annotation_legend = FALSE,
             cluster_rows=FALSE, cluster_cols=FALSE,  ...)
}

