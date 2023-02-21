# 14-basket-ML-training-model.R

#
# load data and libraries
#
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(writexl))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
# sanitized.dds is the longitudinal set using Ensembl v88 annotation
load(file.path(pkg_dir, "data", "all_baskets.rda"))
load(file.path(pkg_dir, "data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "data", "cluster_df.rda"))
annotation <- as.data.frame(rowData(sanitized.dds)) %>%
    dplyr::select(gene_id, gene_name)
set.seed(123)

#
# tidy sanitized.dds$class
#
sanitized.dds$class <- cluster_df %>%
  dplyr::mutate(class = as.character(new_cluster_name)) %>%
  dplyr::mutate(class = if_else( RNA_cluster == "A_Cntr", "Control",
                                class)) %>%
  dplyr::mutate(class = factor(class, levels=c("Control",
                                               levels(new_cluster_name)))) %>%
  pull(class)                                               


#
#  (A) build training model using longitudinal study - control vs. moderate+IG+High
# . (similar to what I did on appendix A)

# 1. get gene id
baskets_genes <- do.call(rbind, all_baskets[c("DUX4-M6", "ECM", "Inflamm", "Complement", "IG")])

# 2. get log10TPM matrix: log10(TPM + 1) of control and moderate+IG+High from sanitized.dds
.predictor_df <- function(predictor, class, dds) {
  sub <- dds[predictor, dds$class %in% class]
  data <- log10(assays(sub)[["TPM"]] +1) %>% t(.) %>%
    as.data.frame() %>%
    add_column(pheno_type = sub$pheno_type) %>%
    rename_with(~ str_replace(.x, "\\..*", ""), starts_with("ENSG"))
}

# 3. ML random forest model fitting: use caret
library(caret)
library(randomForest)
# 1. tools
.rf_knn_fit <- function(predictor_df) {
    # parameter tunning
    fit_ctrl <- trainControl(method = "LOOCV",
                             number = 10,
                             #repeats = 10,
                             classProbs = TRUE)
                             #summaryFunction = twoClassSummary)
    # model - rf
    rf_fit <- train(pheno_type ~ .,
                    data = predictor_df,
                    method = "rf",
                    tuneLength = 35, 
                    preProcess = "pca",
                    trControl = fit_ctrl,
                    metric = "Accuracy",
                    na.action = "na.omit")
    # knn
    knn_fit <- train(pheno_type ~ .,
                    data = predictor_df,
                    method = "knn",
                    trControl = fit_ctrl,
                    metric = "Accuracy",
                    na.action = "na.omit")

    fits <- list(knn=knn_fit, rf=rf_fit) # return list
    return(fits)
}

# tidy the fit accuracy, predictor_des: description of the predictor, mode: control vs. moderate+IG+High
.tidy_fit_accuracy <- function(fits, predictor_description="", model=NULL) {
    data.frame(rf = max(fits$rf$results[, "Accuracy"]),
               knn = max(fits$knn$results[, "Accuracy"]),
               des = predictor_description,
               model=model) 
}


# 2. model fitting: predictor_df (baskets_genes) -> modeling fitting (basket_fit) -> performance
class <- c("Control", "Moderate", "IG-High", "High")
predictor_df <- .predictor_df(predictor = baskets_genes$gene_id_v88, 
                              class = class, dds = sanitized.dds)
basket_fit <- .rf_knn_fit(predictor_df)
.tidy_fit_accuracy(basket_fit, predictor_description="Baskets", model="Moderate+ vs. Controls")


#
# (C) prediction bilateral using the training model based on the longitudinal study
#
load(file.path(pkg_dir, "data", "dds.rda")) # bilate dds
dds$sample_id <- str_replace(dds$sample_name, "[b]*_.*", "")
dds$class <- dds$pheno_type <- "FSHD"
sub <- dds[, ! dds$sample_id %in% c("13-0007R", "13-0009R")]
# prediction
predictor_df_bilat <- .predictor_df(predictor = baskets_genes$gencode_v35, 
                                    class = "FSHD", dds = sub)
predict_rfFit <- predict(basket_fit$rf, newdata=predictor_df_bilat) %>% 
  as.data.frame() %>% dplyr::rename(randomForest.fit = ".") %>%
  dplyr::mutate(sample_id = str_replace(rownames(predictor_df_bilat), "[b]*_.*", ""))
predict_KNNFit <- predict(basket_fit$knn, newdata=predictor_df_bilat) %>% 
  as.data.frame() %>%  dplyr::rename(KNN.fit = ".") %>%
  dplyr::mutate(sample_id = str_replace(rownames(predictor_df_bilat), "[b]*_.*", ""))

#
# (D) gives the classes of bilateral_class (data.frame=sample_name, sample_id, FSHD/control-like)
#
bilat_MLpredict <- data.frame(sample_id = dds$sample_id) %>%
  dplyr::left_join(predict_rfFit, by="sample_id") %>% 
  dplyr::left_join(predict_KNNFit, by="sample_id") %>%
  dplyr::mutate(class = randomForest.fit) %>%
  dplyr::mutate(class = case_when(
    class == "Control" ~ "Control-like",
    class == "FSHD" ~ "Moderate+",
    is.na(class) ~ "Muscle-Low"
  )) %>%
  dplyr::mutate(class=factor(class))
save(bilat_MLpredict, file=file.path(pkg_dir, "data", "bilat_MLpredict.rda"))  