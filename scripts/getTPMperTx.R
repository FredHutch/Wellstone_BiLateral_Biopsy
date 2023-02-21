#'
#' SE must have colData with read_length and lib_size columns
#' 
getTPMperTx <- function(SE,
                        features = rowRanges(SE),
                        counts = assays(SE)$counts,
                        read_length = colData(SE)$read_length) {
    # the input should be a RangeSummarized instance
    if (is.null(read_length)) 
      stop("read_length is essential")
    if (is.null(features)) 
      stop("features (GRangeList) is essential")

    n_features <- length(features)
    n_samples <- ncol(SE)
    
    L <- rep(sum(width(features)), n_samples)
    FL <- sum(width(features))
    #FL <- matrix(L, nrow = n_features, ncol = n_samples)
    X <- sweep(counts, 2, read_length, FUN="*") 
    tx <- X / FL
    T <- colSums(tx)
    TPM <- sweep(tx, 2, T * 1e-6, FUN="/")
}
