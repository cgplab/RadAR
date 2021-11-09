# perform hierarchical clustering on rdr obj ---------------------------------------------------
#' Perform hierarchical clustering of the radiomic dataset
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}}).
#' @param which_data (character) Which data use for the computation of the correlation coefficients. It can be one of the following: "normal",
#' "scaled", "normalized".
#' @param method_dist_row (character) Which method use to computed distance matrix between features (rows). Print available
#' methods by \code{\link{print_distance_methods}}
#' @param method_dist_col (character) Which method use to computed distance matrix between samples (columns). Print available
#' methods by \code{\link{print_distance_methods}}
#' @param method_hcl_row (character) Which method use to generate hierarchical clustering for features (rows). Print available
#' methods by \code{\link{print_hcl_methods}}
#' @param method_hcl_col (character) Which method use to generate hierarchical clustering for samples (columns). Print available
#' methods by \code{\link{print_hcl_methods}}
#'
#' @return An updated rdr (a RadAR object)
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
do_hierarchical_clustering <- function(rdr = NULL,
                                       which_data = "scaled",
                                       method_dist_row = "euclidean",
                                       method_dist_col = "correlation.pearson",
                                       method_hcl_row = "ward.D",
                                       method_hcl_col = "ward.D"

)
{
  methods_hcl <- c("ward.D",
                   "ward.D2",
                   "single",
                   "complete",
                   "average" ,
                   "mcquitty",
                   "median",
                   "centroid")
  methods_dist <- c("euclidean",
                    "maximum",
                    "manhattan",
                    "canberra",
                    "binary",
                    "minkowski",
                    "correlation.pearson",
                    "correlation.spearman",
                    "correlation.kendall"
  )

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(which_data %in% c("normal", "scaled", "normalized"), msg = "[RadAR] Error: Invalid data type")
  assertthat::assert_that(method_dist_row %in% methods_dist & length(method_dist_row) == 1, msg = "[RadAR] Error: Invalid method for calculating distance")
  assertthat::assert_that(method_dist_col %in% methods_dist & length(method_dist_col) == 1, msg = "[RadAR] Error: Invalid method for calculating distance")
  assertthat::assert_that(method_hcl_row %in% methods_hcl & length(method_hcl_row) == 1, msg = "[RadAR] Error: Invalid method for performing hierarcgical clustering")
  assertthat::assert_that(method_hcl_col %in% methods_hcl & length(method_hcl_col) == 1, msg = "[RadAR] Error: Invalid method for performing hierarcgical clustering")

  if (which_data == "normal") {
    data <- assays(rdr)$values
  }
  if (which_data == "scaled") {
    data <- assays(rdr)$scaled_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet scaled")
  }
  if (which_data == "normalized") {
    data <- assays(rdr)$norm_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet normalized")
  }

  ## calc distance
  dist0 <- unlist(strsplit(method_dist_col, ".", fixed = T))
  if (dist0[1] != "correlation") {
    dist_col <- dist(t(data), method = method_dist_col)
  } else  {
    dist_col <- as.dist(1 - cor(data, method = dist0[2], use = "na"))
  }
  dist0 <- unlist(strsplit(method_dist_row, ".", fixed = T))
  if (dist0[1] != "correlation") {
    dist_row <- dist(data, method = method_dist_row)
  } else  {
    dist_row <- as.dist(1 - cor(t(data), method = dist0[2], use = "na"))
  }
  hcl_col <- hclust(dist_col, method = method_hcl_col)
  hcl_row <- hclust(dist_row, method = method_hcl_row)
  hcl <- list (hcl_col = hcl_col,
               hcl_row = hcl_row)

  metadata(rdr)$hcl <- hcl
  metadata(rdr)$which_data_hcl <- which_data
  return(rdr)
}

# calc differential radiomics -------------------------------------------------------------------
#' Perform differential radiomic analysis
#'
#' This function implements two non-parametric statistical tests (Wilcoxon-Mann-Whitney and AUC)
#' to perform differential analysis of radiomic features. In case of Wilcoxon-Mann-Whitney p-values can be adjusted using
#' different methods.
#'
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}}).
#' @param conditions (numeric, character) Vector of labels for each sample in rdr. Should have the same
#' length of \code{ncol(rdr)}. Required.
#' @param which_data (character) Which data use for plot. It can be one of the following: "normal", "scaled", "normalized".
#' @param method (character) Which statistical methods use for differential analysis.
#' It can be "wilcox" (Wilcoxon-Mann-Whitney), "auc" (Area Under the Curve) or "kruskal-wallis".
#' @param adjust_pvalues_by (character) Which method use to correct p-values from wilcox test. Print available methods by
#' \code{\link{p.adjust.methods}}.
#' @param thr_pvalue (numeric) P-value threshold to identify statistically significant features for wilcox or kruskal-wallis.
#' It should be in the range (0, 1].
#' @param thr_auc (numeric) AUC threshold to identify statistically significant features.
#' It should be in the range (0.5,1].
#'
#' @return An updated rdr (a RadAR object)
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
calc_differential_radiomics <- function(rdr = NULL,
                                        conditions = NULL,
                                        which_data = "scaled",
                                        method = "wilcox",
                                        adjust_pvalues_by = "BH",
                                        thr_pvalue = NA,
                                        thr_auc = NA)

{

  assertthat::assert_that(length(rdr) > 0,
                          msg = "[RadAR] Error: rdr object required")
  if (!is.na(thr_pvalue)) {
    assertthat::assert_that(thr_pvalue > 0,
                            msg = "[RadAR] Error: thr_pvalue should be positive")
  }
  if (!is.na(thr_auc)) {
    assertthat::assert_that(thr_auc > 0.5 & thr_auc <= 1,
                            msg = "[RadAR] Error: thr_auc should be in the range (0.5,1]")
  }
  assertthat::assert_that(length(conditions) == ncol(rdr),
                          msg = "[RadAR] Error: conditions need to be specified and should be a vector of length ncol(rdr)")
  assertthat::assert_that(which_data %in% c("normal", "scaled", "normalized"),
                          msg = "[RadAR] Error: Invalid data type")
  assertthat::assert_that(method %in% c("wilcox", "auc", "kruskal-wallis"),
                          msg = "[RadAR] Error: invalid method to perform differential statistics")

  if (which_data == "normal") {
    data <- assays(rdr)$values
  }
  if (which_data == "scaled") {
    data <- assays(rdr)$scaled_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: feature values have not been yet scaled")
  }
  if (which_data == "normalized") {
    data <- assays(rdr)$norm_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: feature values have not been yet normalized")
  }

  if (method %in% c("wilcox", "auc")) {

    assertthat::assert_that(length(unique(conditions)) == 2 |
                              (length(unique(conditions)) == 3 & any(is.na(unique(conditions)))),
                            msg = "[RadAR] Error: two conditions should be indicated (test vs ctrl)")

    unique_conditions <- unique(conditions)
    if (any(is.na(unique_conditions))) {
      unique_conditions <- unique_conditions[-which(is.na(unique_conditions))]
    }

    median_cond1 <- apply(data[, which(conditions == unique_conditions[1])], 1, median, na.rm = T)
    median_cond2 <- apply(data[, which(conditions == unique_conditions[2])], 1, median, na.rm = T)

    if (method == "wilcox") {
      assertthat::assert_that(adjust_pvalues_by %in% p.adjust.methods,
                              msg = "[RadAR] Invalid method for adjusting p-values")
      pvalues <- apply(data, 1, calc_wilcox, conditions)
      adjusted_pvalues <- p.adjust(pvalues, method = adjust_pvalues_by)
      rowData(rdr)$wilcox_test_pvalue <- adjusted_pvalues
      metadata(rdr)$which_data_wilcox <- which_data
      metadata(rdr)$conditions_wilcox <- conditions

      if (!is.na(thr_pvalue)) {
        ix_features <- which(rowData(rdr)$wilcox_test_pvalue <= thr_pvalue  )
        if (length(ix_features) == 0) {
          message("Warning: [RadAR] No statistically significant features found")
        }
        flag_condition <- rep("", nrow(rdr))
        flag_condition[which(median_cond1 > median_cond2 &
                               adjusted_pvalues <= thr_pvalue)] <- paste("Higher in", unique_conditions[1])
        flag_condition[which(median_cond2 > median_cond1 &
                               adjusted_pvalues <= thr_pvalue)] <- paste("Higher in", unique_conditions[2])
        rowData(rdr)$wilcox_test_description <- flag_condition
      }
    }

    if (method == "auc") {
      auc <- apply(data, 1, calc_auc, conditions)
      rowData(rdr)$auc_value <- auc
      metadata(rdr)$which_data_auc <- which_data
      metadata(rdr)$conditions_auc <- conditions
      if (!is.na(thr_auc)) {
        ix_features <- which(rowData(rdr)$auc_value >= thr_auc | rowData(rdr)$auc_value <= (1-thr_auc) )
        if (length(ix_features) == 0) {
          message("Warning: [RadAR] No statistically significant features found")
        }
        flag_condition <- rep("", nrow(rdr))
        flag_condition[which(auc >= thr_auc)] <- paste("Higher in", unique_conditions[1])
        flag_condition[which(auc <= (1-thr_auc) )] <- paste("Higher in", unique_conditions[2])
        rowData(rdr)$auc_description <- flag_condition
      }
    }
  }

  if (method == "kruskal-wallis") {
    assertthat::assert_that(length(unique(conditions)) > 2 |
                              (length(unique(conditions)) > 3 & any(is.na(unique(conditions)))),
                            msg = "[RadAR] Error: at least three conditions should be indicated (test vs ctrl)")

    assertthat::assert_that(adjust_pvalues_by %in% p.adjust.methods,
                            msg = "[RadAR] Invalid method for adjusting p-values")

    pvalues <- apply(data, 1, calc_kruskal_wallis, conditions)
    adjusted_pvalues <- p.adjust(pvalues, method = adjust_pvalues_by)
    rowData(rdr)$kruskal_wallis_test_pvalue <- adjusted_pvalues
    metadata(rdr)$which_data_kruskal_wallis <- which_data
    metadata(rdr)$conditions_kruskal_wallis <- conditions

    if (!is.na(thr_pvalue)) {
      ix_features <- which(rowData(rdr)$kruskal_wallis_test_pvalue <= thr_pvalue  )
      if (length(ix_features) == 0) {
        message("Warning: [RadAR] No statistically significant features found")
      }
    }
  }


  return (rdr)

}


calc_wilcox <- function(x, conditions) {
  y <- wilcox.test(x ~ conditions)
  return (y$p.value)
}

calc_kruskal_wallis <- function(x, conditions) {
  y <- kruskal.test(x ~ conditions)
  return (y$p.value)
}

calc_auc <- function(x, conditions) {
  y <- wilcox.test(x ~ conditions)
  return (y$statistic/prod(table(conditions)))
}

#' Select top radiomic features according to a given statistics
#'
#' This function create a summary table of most significant features based on previously computed statistics
#' (wilcox, AUC, concordance index or cox regression)
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}}).
#' @param which_statistics (character) Select top features based on one of the following pre-computed statistics:
#' "wilcox", "AUC", "concordance" or "cox".
#' @param thr_pvalue (numeric) P-value threshold to identify statistically significant features from wilcox.
#' It should be in the range (0, 1].
#' @param thr_auc (numeric) AUC threshold to identify statistically significant features from AUC
#' It should be in the range (0.5, 1].
#' @param thr_concordance (numeric) Threshold to identify statistically significant features from concordance analysis.
#' It should be in the range (0, .5].
#' @param thr_cox_zvalue (numeric) Z threshold (Wald statistics) to identify statistically significant features from cox regression analysis.
#' It should be in the range (0, inf).
#' @param write_to (character) If specified, filename to output top statistically significant features (tab-delimited).
#'
#' @return A table of class \code{\link{tibble}} reporting top features.
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
select_top_features  <- function(rdr = NULL,
                                 which_statistics = "wilcox",
                                 thr_pvalue = 5e-2,
                                 thr_auc = .8,
                                 thr_concordance = 0.10,
                                 thr_cox_zvalue = 3,
                                 write_to = NULL
)

{

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(thr_pvalue > 0, msg = "[RadAR] Error: thr_pvalue_wilcox should be positive")
  assertthat::assert_that(thr_cox_zvalue > 0, msg = "[RadAR] Error: thr_cox_zvalue should be positive")
  assertthat::assert_that(thr_auc > 0.5 & thr_auc <= 1, msg = "[RadAR] Error: auc_thr should be in the range (0.5,1]")
  assertthat::assert_that(thr_concordance > 0 & thr_concordance < 0.5, msg = "[RadAR] Error: thr_concordance should be in the range (0, .5)")
  assertthat::assert_that(length(which_statistics) == 1, msg = "[RadAR] Error: Only one statistics should be indicated")
  assertthat::assert_that(which_statistics %in% c("wilcox", "auc", "concordance", "cox", "kruskal-wallis"),
                          msg = "[RadAR] Error: Invalid statistcs")

  if (which_statistics == "kruskal-wallis") {
    which_statistics <- "kruskal_wallis"
  }
  which_data_wilcox <- metadata(rdr)$which_data_wilcox
  which_data_auc <- metadata(rdr)$which_data_auc
  which_data_cox <- metadata(rdr)$which_data_cox
  which_data_concordance <- metadata(rdr)$which_data_concordance
  which_data_kruskal_wallis <- metadata(rdr)$which_data_kruskal_wallis

  conditions_wilcox <- metadata(rdr)$conditions_wilcox
  conditions_auc <- metadata(rdr)$conditions_auc
  conditions_kruskal_wallis <- metadata(rdr)$conditions_wilcox

  which_data <- get(paste0("which_data", "_", which_statistics))

  if (which_statistics == "wilcox") {
    assertthat::assert_that(any("wilcox_test_pvalue" %in% colnames(rowData(rdr))),
                            msg = "[RadAR] Error: Wilcox-Mann-Whitney test has not been yet computed")
    thr_pvalue_wilcox <- thr_pvalue
  }
  if (which_statistics == "auc") {
    assertthat::assert_that(any("auc_value" %in% colnames(rowData(rdr))),
                            msg = "[RadAR] Error: AUC test has not been yet computed")
  }
  if (which_statistics == "kruskal_wallis") {
    assertthat::assert_that(any("kruskal_wallis_test_pvalue" %in% colnames(rowData(rdr))),
                            msg = "[RadAR] Error: Kruskal-Wallis test has not been yet computed")

  }

  if (which_data == "normal") {
    data <- assays(rdr)$values
  }
  if (which_data == "scaled") {
    data <- assays(rdr)$scaled_values
    assertthat::assert_that(!is.null(data),
                            msg = "[RadAR] Error: Feature values have not been scaled yet")
  }
  if (which_data == "normalized") {
    data <- assays(rdr)$norm_values
    assertthat::assert_that(!is.null(data),
                            msg = "[RadAR] Error: Feature values have not been normalized yet")
  }

  if (which_statistics %in% c("wilcox", "auc")) {
    conditions <- get(paste0("conditions", "_", which_statistics))
    unique_conditions <- unique(conditions)
    if (any(is.na(unique_conditions))) {
      unique_conditions <- unique_conditions[-which(is.na(unique_conditions))]
    }
    median_cond1 <- apply(data[, which(conditions == unique_conditions[1])], 1, median, na.rm = T)
    median_cond2 <- apply(data[, which(conditions == unique_conditions[2])], 1, median, na.rm = T)
  }

  if (which_statistics == "wilcox") {

    ix_features <- which(rowData(rdr)$wilcox_test_pvalue <= thr_pvalue_wilcox)
    assertthat::assert_that(length(ix_features) > 0, msg = "[RadAR] Error: No statistically significant features found")
    flag_condition <- rep("", nrow(rdr))
    flag_condition[which(median_cond1 > median_cond2 &
                           rowData(rdr)$wilcox_test_pvalue <= thr_pvalue_wilcox)] <- paste("Higher in", unique_conditions[1])
    flag_condition[which(median_cond2 > median_cond1 &
                           rowData(rdr)$wilcox_test_pvalue <= thr_pvalue_wilcox)] <- paste("Higher in", unique_conditions[2])

    df_top_features <- data.frame(feature_name = rowData(rdr)$feature_name[ix_features],
                                  image_type = rowData(rdr)$image_type[ix_features],
                                  feature_type = rowData(rdr)$feature_type[ix_features],
                                  feature_description = rowData(rdr)$feature_description[ix_features],
                                  statistics = rowData(rdr)$wilcox_test_pvalue[ix_features],
                                  statistics_description = flag_condition[ix_features],
                                  median_cond1[ix_features],
                                  median_cond2[ix_features]
    )
    colnames(df_top_features)[c(7, 8)] <- c(paste0("median_in_", unique_conditions[1]),
                                            paste0("median_in_", unique_conditions[2]))
    rownames(df_top_features) <- rownames(rdr)[ix_features]
    df_top_features <- df_top_features[order(df_top_features$statistics), ]
  }

  if (which_statistics == "auc") {

    ix_features <- which(rowData(rdr)$auc_value >= thr_auc | rowData(rdr)$auc_value <= (1-thr_auc) )
    assertthat::assert_that(length(ix_features) > 0, msg = "[RadAR] Error: No statistically significant features found")

    flag_condition <- rep("", nrow(rdr))
    flag_condition[which(median_cond1 > median_cond2 &
                           abs(rowData(rdr)$auc_value-0.5) > (thr_auc-0.5))] <- paste("Higher in", unique_conditions[1])
    flag_condition[which(median_cond2 > median_cond1 &
                           abs(rowData(rdr)$auc_value-0.5) > (thr_auc-0.5))] <- paste("Higher in", unique_conditions[2])

    df_top_features <- data.frame(feature_name = rowData(rdr)$feature_name[ix_features],
                                  image_type = rowData(rdr)$image_type[ix_features],
                                  feature_type = rowData(rdr)$feature_type[ix_features],
                                  feature_description = rowData(rdr)$feature_description[ix_features],
                                  statistics = rowData(rdr)$auc_value[ix_features],
                                  statistics_description = flag_condition[ix_features],
                                  median_cond1[ix_features],
                                  median_cond2[ix_features]
    )
    colnames(df_top_features)[c(7, 8)] <- c(paste0("median_in_", unique_conditions[1]),
                                            paste0("median_in_", unique_conditions[2]))
    rownames(df_top_features) <- rownames(rdr)[ix_features]
    df_top_features <- df_top_features[order(abs(df_top_features$statistics-0.5), decreasing = T), ]
  }

  if (which_statistics == "kruskal_wallis") {

    ix_features <- which(rowData(rdr)$kruskal_wallis_test_pvalue <= thr_pvalue)
    assertthat::assert_that(length(ix_features) > 0,
                            msg = "[RadAR] Error: No statistically significant features found")

    df_top_features <- data.frame(feature_name = rowData(rdr)$feature_name[ix_features],
                                  image_type = rowData(rdr)$image_type[ix_features],
                                  feature_type = rowData(rdr)$feature_type[ix_features],
                                  feature_description = rowData(rdr)$feature_description[ix_features],
                                  statistics = rowData(rdr)$kruskal_wallis_test_pvalue[ix_features]
    )
    rownames(df_top_features) <- rownames(rdr)[ix_features]
    df_top_features <- df_top_features[order(df_top_features$statistics), ]
  }
  if (which_statistics == "concordance") {

    ix_features <- which(abs(rowData(rdr)$concordance_index - 0.5) >= thr_concordance  )
    assertthat::assert_that(length(ix_features) > 0, msg = "[RadAR] No statistically significant features found")

    flag_condition <- rep("", nrow(rdr))
    flag_condition[intersect(ix_features, which(rowData(rdr)$concordance_index - .5 < 0))] <- "worse increasing"
    flag_condition[intersect(ix_features, which(rowData(rdr)$concordance_index - .5 > 0))] <- "worse decreasing"

    df_top_features <- data.frame(feature_name = rowData(rdr)$feature_name[ix_features],
                                  image_type = rowData(rdr)$image_type[ix_features],
                                  feature_type = rowData(rdr)$feature_type[ix_features],
                                  feature_description = rowData(rdr)$feature_description[ix_features],
                                  statistics = rowData(rdr)$concordance_index[ix_features],
    )
    rownames(df_top_features) <- rownames(rdr)[ix_features]
    df_top_features <- df_top_features[order(abs(df_top_features$statistics-0.5), decreasing = T), ]
  }

  if (which_statistics == "cox") {

    ix_features <- which(abs(rowData(rdr)$cox_regression_statistics) >= thr_cox_zvalue)
    assertthat::assert_that(length(ix_features) > 0, msg = "[RadAR] No statistically significant features found")

    flag_condition <- rep("", nrow(rdr))
    flag_condition[intersect(ix_features, which(rowData(rdr)$cox_regression_statistics >= thr_cox_zvalue))] <- "worse increasing"
    flag_condition[intersect(ix_features, which(rowData(rdr)$cox_regression_statistics <= -thr_cox_zvalue))] <- "worse decreasing"
    rowData(rdr)$cox_regression_description <- flag_condition

    df_top_features <- data.frame(feature_name = rowData(rdr)$feature_name[ix_features],
                                  image_type = rowData(rdr)$image_type[ix_features],
                                  feature_type = rowData(rdr)$feature_type[ix_features],
                                  feature_description = rowData(rdr)$feature_description[ix_features],
                                  statistics = rowData(rdr)$cox_regression_statistics[ix_features],
                                  statistics_description = flag_condition[ix_features]
    )
    rownames(df_top_features) <- rownames(rdr)[ix_features]
    df_top_features <- df_top_features[order(abs(df_top_features$statistics), decreasing = T), ]
  }
  if (!is.null(write_to)) {
    write.table(df_top_features, file = write_to, quote = F, row.names = F, sep ="\t")
  }
  df_top_features <- tibble::as_tibble(df_top_features)
  return (df_top_features)

}

#' Compute concordance index for radiomic features
#'
#' This function calculates concordance index of radiomic features.
#' Response should be a survival object obtained by \code{\link{Surv}}
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}}).
#' @param surv_obj An object of class \code{\link{Surv}}.
#' @param thr_concordance (numeric) Threshold to identify statistically significant features from concordance analysis.
#' It should be in the range (0, .5].
#' @param which_data (character) Which data use to compute concordance index. It can be one of the following: "normal",
#' "scaled", "normalized".
#'
#' @return An updated rdr (a RadAR object)
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
calc_concordance_index <- function(rdr = NULL,
                                   surv_obj = NULL,
                                   thr_concordance = NA,
                                   which_data = "scaled"
)
{

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(is.Surv(surv_obj) > 0, msg = "[RadAR] Error: surv_obj object required")
  assertthat::assert_that(which_data %in% c("normal", "scaled", "normalized"), msg = "[RadAR] Error: Invalid data type")

  if (!is.na(thr_concordance)) {
    assertthat::assert_that(thr_concordance > 0 & thr_concordance < 0.5, msg = "[RadAR] Error: thr_concordance should be in the range (0, .5)")
  }

  if (which_data == "normal") {
    data <- assays(rdr)$values
  }
  if (which_data == "scaled") {
    data <- assays(rdr)$scaled_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet scaled")
  }
  if (which_data == "normalized") {
    data <- assays(rdr)$norm_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet normalized")
  }

  concordance_index <- rep(NA, nrow(rdr))
  for (i in 1: nrow(rdr)) {
    concordance_index[i] <- concordance(surv_obj ~ data[i, ])$concordance
    #    concordance_index[i] <- concordance(surv_obj ~ data[i, ], data = data.frame(rdr_original@assays$data$values[i, ]))$concordance
  }
  rowData(rdr)$concordance_index <- concordance_index
  metadata(rdr)$which_data_concordance <- which_data

  if (!is.na(thr_concordance)) {
    ix_features <- which(abs(rowData(rdr)$concordance_index - .5) >= thr_concordance)
    if (length(ix_features) == 0) {
      message("Warning: [RadAR] No statistically significant features found")
    }
    flag_condition <- rep("", nrow(rdr))
    flag_condition[intersect(ix_features, which(rowData(rdr)$concordance_index - .5 < 0))] <- "worse increasing"
    flag_condition[intersect(ix_features, which(rowData(rdr)$concordance_index - .5 > 0))] <- "worse decreasing"
    rowData(rdr)$concordance_description <- flag_condition
  }
  return(rdr)

}


#' Fit a Cox regression model using radiomic features
#'
#' This function fits a Cox Proportional-Hazards Model using radiomic features as independent variables.
#' Response should be an object obtained by \code{\link{Surv}}
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}}).
#' @param surv_obj An object of class \code{\link{Surv}}.
#' @param thr_cox_zvalue (numeric) Z threshold (Wald statistics) to identify statistically significant features from cox regression analysis.
#' It should be in the range (0, inf).
#' @param which_data (character) Which data use to compute concordance index. It can be one of the following: "normal",
#' "scaled", "normalized".
#' @param multiple_regression (character). Cox regression can be "multiple" (surv ~ all features) or "single" (surv ~ each feature).
#'
#' @return An updated rdr (a RadAR object)
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
calc_cox_regression <- function(rdr = NULL,
                                surv_obj = NULL,
                                thr_cox_zvalue = NA,
                                which_data = "scaled",
                                multiple_regression = F
)
{

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(is.Surv(surv_obj) > 0, msg = "[RadAR] Error: surv_obj object required")
  assertthat::assert_that(which_data %in% c("normal", "scaled", "normalized"), msg = "[RadAR] Error: Invalid data type")
  assertthat::assert_that(is.logical(multiple_regression), msg = "[RadAR] Error: multiple_regression should be logical")

  if (!is.na(thr_cox_zvalue)) {
    assertthat::assert_that(thr_cox_zvalue > 0, msg = "[RadAR] Error: thr_cox_zvalue should be positive")
  }

  if (which_data == "normal") {
    data <- assays(rdr)$values
  }
  if (which_data == "scaled") {
    data <- assays(rdr)$scaled_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet scaled")
  }
  if (which_data == "normalized") {
    data <- assays(rdr)$norm_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet normalized")
  }

  if (!multiple_regression) {
    cox_regression_statistics <- rep(NA, nrow(rdr))
    cox_hazard_ratio <- rep(NA, nrow(rdr))
    for (i in 1: nrow(rdr)) {
      tmp <- summary(coxph(surv_obj ~ data[i, ]))$coefficients
      cox_regression_statistics[i] <- tmp[1, 4]
      cox_hazard_ratio[i] <- tmp[1, 2]
    }
  }
  if (multiple_regression) {
    tmp <- summary(coxph(surv_obj ~ t(data)))$coefficients
    cox_regression_statistics <- tmp[, 4]
    cox_hazard_ratio <-  tmp[, 2]
  }

  rowData(rdr)$cox_regression_statistics <- cox_regression_statistics
  rowData(rdr)$cox_hazard_ratio <- cox_hazard_ratio
  metadata(rdr)$which_data_cox <- which_data
  metadata(rdr)$which_cox_regression <- ifelse(multiple_regression, "multiple", "single")

  if (!is.na(thr_cox_zvalue)) {
    ix_features <- which(abs(rowData(rdr)$cox_regression_statistics) >= thr_cox_zvalue)
    if (length(ix_features) == 0) {
      message("Warning: [RadAR] No statistically significant features found")
    }
    flag_condition <- rep("", nrow(rdr))
    flag_condition[intersect(ix_features, which(rowData(rdr)$cox_regression_statistics >= thr_cox_zvalue))] <- "worse increasing"
    flag_condition[intersect(ix_features, which(rowData(rdr)$cox_regression_statistics <= -thr_cox_zvalue))] <- "worse decreasing"
    rowData(rdr)$cox_regression_description <- flag_condition
  }
  return(rdr)
}

#' Test radiomic signature by Cox regression model
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}}).
#' @param surv_obj An object of class \code{\link{Surv}}.
#' @param signature (character) vector of radiomic features defining the signature.
#' @param signature_rdr  RadAR object (class \code{\link{SummarizedExperiment}}) obtained by
#' \code{\link{do_feature_selection}}. Required to test a model generated using
#'  glmnet-cox or glmnet-binom.
#' @param which_data (character) Which data use to compute concordance index. It can be one of the following: "normal",
#' "scaled", "normalized".
#'
#' @return In case of mRMR, hcl and pca methods, the function returns the results of the cox regression model, including the concordance index.
#' In case of glmnet-cox method, the function returns a list of two elements, including 1) the prediction of the glmnet model on the new data
#' and 2) the computation of concordance index.
#'
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
test_radiomic_signature <- function(rdr = NULL,
                                    surv_obj = NULL,
                                    signature = NULL,
                                    signature_rdr = NULL,
                                    which_data = "scaled"
)
{

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(is.Surv(surv_obj) > 0, msg = "[RadAR] Error: surv_obj object required")
  assertthat::assert_that(length(signature) > 0 | length(signature_rdr) > 0, msg = "[RadAR] Error: signature object required")
  assertthat::assert_that(which_data %in% c("normal", "scaled", "normalized"), msg = "[RadAR] Error: Invalid data type")

  if (length(signature_rdr) > 0) {
    signature <- rownames(signature_rdr)
  }

  if (!all(signature %in% rownames(rdr))) {
    ix <- which(signature %in% rownames(rdr) == F)
    assertthat::assert_that(1<0, msg = paste0("[RadAR] Error:", toString(signature[ix], " not in rdr")))
  }

  if (which_data == "normal") {
    data <- assays(rdr)$values
  }
  if (which_data == "scaled") {
    data <- assays(rdr)$scaled_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet scaled")
  }
  if (which_data == "normalized") {
    data <- assays(rdr)$norm_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet normalized")
  }

  data_sign <- data[signature, ]

  if (length(signature_rdr) > 0) {
    if (length(metadata(signature_rdr)$glmnet_model) == 1) {
      res <- summary(coxph(surv_obj ~ t(data_sign)))
    } else {
      if (metadata(signature_rdr)$glmnet_model_lambda == "min") {
        my_s <- "lambda.min"
      } else {
        my_s <- "lambda.1se"
      }
      cvfit <- metadata(signature_rdr)$glmnet_model
      res <- list()
      res$pred <- predict(cvfit,
                          newx = t(data),
                          s = my_s,
                          type = "response")

      res$concordance <- Cindex(pred = res$pred,
                       y = surv_obj)

    }
  } else {
    res <- summary(coxph(surv_obj ~ t(data_sign)))
  }

  return(res)
}

#' Perform feature selection based on different methods
#'
#' This function implements different methods to perform feature selection of radiomic datasets.
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}}).
#' @param n_features (numeric) Number of features to be selected. Required.
#' @param select_by (character) Which criteria use to select informative radiomic features within clusters
#' of similar (i.e., redundant) features. It can be one of the following: "variability", "random", "concordance".
#' @param method (character) Which  method use to identify redundant features. It can be one of the following:
#' "mRMR" (minimum-redundancy-maximum-relevance),"hcl" (hierarchical clustering of correlation matrix), "pca" (K-means applied to Principal Component Analysis),
#' "glmnet-cox" (generalized linear model via penalized maximum likelihood (glmnet) fitting cox regression model),
#' "glmnet-binonial" (glmnet fitting binomial regression model),
#' Using mRMR, this function works as a wrapper to \code{\link{mRMR}} package.
#' Using glmnet-*, this function works as a wrapper to \code{\link{glmnet}} package.
#'
#' @param surv_obj An object of class \code{\link{Surv}}. Required if select_by is "concordance".
#' @param which_data (character) Which data use to compute concordance index. It can be one of the following: "normal",
#' "scaled", "normalized".
#' @param corr_measure (character) Which method use to calculate correlation. It can be one of the following:
#' "pearson", "kendall", "spearman".
#' @param min_features_per_group (numeric) Minimum number of features for each cluster.
#' @param thr_pca_cum_prop (numeric) Threshold to select number of components based on
#' cumulative proportion of explained variance criterion.
#' @param response (numeric) A response variable, required if any of mRMR or glmnet-binomial or
#' methods are used.
#' @param lambda (character) In glmnet, it controls the overall strength of the penalty.
#' Possible values are "min" or "1se" (1 standard deviation). For more details see \code{\link{glmnet}}
#' @param alpha (numeric) In glmnet, it controls elastic-net penalty.
#' Typical values are 0 (ridge) or 1 (lasso). For more details see \code{\link{glmnet}}.
#'
#' @return A list including two elements:
#' `rdr`: the updated (reduced) rdr (a RadAR object)
#' `signature`: the radiomic features included in the signature
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
do_feature_selection <- function(rdr = NULL,
                                 n_features = NULL,
                                 select_by = NULL,
                                 method =  "hcl",
                                 surv_obj = NULL,
                                 which_data = "scaled",
                                 corr_measure = "pearson",
                                 min_features_per_group = 5,
                                 thr_pca_cum_prop = .8,
                                 response = NULL,
                                 lambda = "min",
                                 alpha = NULL
)
{
  set.seed(1)
  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] rdr object required")
  assertthat::assert_that(which_data %in% c("normal", "scaled", "normalized"), msg = "[RadAR] Invalid data type")
  assertthat::assert_that(corr_measure %in% c("pearson", "kendall", "spearman"), msg = "[RadAR] Invalid corr_measure setting")
  assertthat::assert_that(length(method) > 0 , msg = "[RadAR] method required")
  assertthat::assert_that(method %in% c("hcl", "pca", "mRMR",
                                        "glmnet-cox", "glmnet-binomial") ,
                          msg = "[RadAR] Invalid method")

  if (method %in% c("hcl", "pca", "mRMR")) {
    assertthat::assert_that(min_features_per_group > 0 & min_features_per_group < (nrow(rdr)-n_features) ,
                            msg = "[RadAR] min_features_per_group should be in the range [1, (total features - n_features)]")
    assertthat::assert_that(!is.null(n_features), msg = "[RadAR] n_features required for methods hcl, pca or mRMR")

  }

  if (method %in% c("mRMR", "glmnet-cox", "glmnet-binomial") == F) {
    assertthat::assert_that(length(select_by) > 0 , msg = "[RadAR] select_by required")
    assertthat::assert_that(select_by %in% c("variability", "random", "concordance"), msg = "[RadAR] Invalid method for feature prioritization")
  }
  if (method == "mRMR") {
    assertthat::assert_that(length(response) > 0, msg = "[RadAR] Error: For mRMR a response variable is required")
    assertthat::assert_that(length(response) == ncol(rdr), msg = "[RadAR] Error: Response variable should have same length of ncol(rdr)")
    assertthat::assert_that(length(unique(response)) == 2  |
                              (length(unique(response)) == 3 & any(is.na(unique(response)))),
                            msg = "[RadAR] Error: Response variable should have two states")
    select_by <- "none"
  }
  if (method %in% c("glmnet-cox", "glmnet-binomial")) {
    if (length(alpha) > 0) {
      assertthat::assert_that(alpha >= 0 & alpha <= 1,
                              msg = "[RadAR] Invalid value for alpha. It should be in the range [0, 1].")
    }
    assertthat::assert_that(lambda %in% c("min", "1se"),
                            msg = "[RadAR] Invalid value for lamba. It should min or 1se.")
    select_by <- "none"
  }
  if (select_by == "interpretability" ) {
    assertthat::assert_that(length(rowData(rdr)$interpretability_score) > 0, msg = "[RadAR] Error: Interpretability of features is not available")
  }
  if (select_by == "concordance") {
    assertthat::assert_that(length(surv_obj) > 0, msg = "[RadAR] Error: For concordance a surv_obj is required")
  }
  if (method == "glmnet-cox") {
    assertthat::assert_that(length(surv_obj) > 0, msg = "[RadAR] Error: For glmnet-cox a surv_obj is required")
  }
  if (method %in% c("glmnet-binomial")) {
    assertthat::assert_that(length(response) > 0, msg = "[RadAR] Error: For glmnet-binomial a response variable is required")
    assertthat::assert_that(length(response) == ncol(rdr), msg = "[RadAR] Error: Response variable should have same length of ncol(rdr)")
    assertthat::assert_that(length(unique(response)) == 2 |
                              (length(unique(response)) == 3 & any(is.na(unique(response)))),
                            msg = "[RadAR] Error: Response variable should have two states")
  }

  if (which_data == "normal") {
    data <- assays(rdr)$values
  }
  if (which_data == "scaled") {
    data <- assays(rdr)$scaled_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet scaled")
  }
  if (which_data == "normalized") {
    data <- assays(rdr)$norm_values
    assertthat::assert_that(!is.null(data), msg = "[RadAR] Error: Feature values have not been yet normalized")
  }

  if (method %in% c("mRMR", "glmnet-binomial")) {
    if(any(is.na(response))) {
      ixna <- which(is.na(response))
      response <- response[-ixna]
      data <- data[, -ixna]
    }
  }
  if (method == "glmnet-cox" | select_by == "concordance") {
    if(any(is.na(surv_obj[, 1]))) {
      ixna <- which(is.na(surv_obj[, 1]))
      surv_obj <- surv_obj[-ixna, ]
      data <- data[, -ixna]
    }
  }


  if (method == "hcl") {
    dist_row <- as.dist( 1- cor(t(data), use = "na", method = corr_measure))
    hcl <- hclust(dist_row, method = "ward.D")
    clusts <- cutree(hcl, k = n_features)
  }

  if (method == "pca") {
    prc <- prcomp(data)
    summ_prc <- summary(prc)
    n_components <- which(summ_prc$importance[3, ] > thr_pca_cum_prop)[1]-1
    kmeans_pca <- kmeans(prc$x[, 1:n_components], centers = n_features)
    clusts <- kmeans_pca$cluster
  }

  if (method %in% c("pca", "hcl")) {
    assertthat::assert_that(all(table(clusts) > min_features_per_group),
                            msg = "[RadAR] Error: one or more groups with too few features. Decrease min_features_per_group or n_features")

    if (select_by == "variability") {
      q1 <- apply(data, 1, quantile, .25, na.rm = T)
      q2 <- apply(data, 1, quantile, .5, na.rm = T)
      q3 <- apply(data, 1, quantile, .75, na.rm = T)
      variability <- abs((q3-q1)/q2)
      sel_features <- rep("", n_features)
      for (i in 1: n_features) {
        sel_features[i] <- names(which.max(variability[which(clusts == i)]))
      }
    }

    if (select_by == "random") {
      sel_features <- rep("", n_features)
      for (i in 1: n_features) {
        sel_features[i] <- names(sample(which(clusts == i), size = 1))
      }
    }

    if (select_by == "interpretability") {
      sel_features <- rep("", n_features)
      for (i in 1: n_features) {
        sel_features[i] <- names(which.max(rowData(rdr)$interpretability_score[which(clusts == i)]))
      }
    }
    if (select_by == "concordance") {
      concordance_index <- rep(NA, nrow(rdr))
      for (i in 1: nrow(rdr)) {
        concordance_index[i] <- concordance(surv_obj ~ data[i, ])$concordance
      }
      names(concordance_index) <- rownames(rdr)
      sel_features <- rep("", n_features)
      for (i in 1: n_features) {
        sel_features[i] <- names(which.max(abs(concordance_index[which(clusts == i)]-0.5) ))
      }
    }
  }
  if (method == "mRMR") {
    mrmr_in_data <- t(rbind(data, response))
    mrmr_data <- mRMRe::mRMR.data(data = data.frame(mrmr_in_data))
    mrmr_res <- mRMRe::mRMR.classic("mRMRe.Filter", data = mrmr_data,
                                    target_indices = ncol(mrmr_in_data),
                                    feature_count = n_features)
    ix_features <- as.numeric(unlist(solutions(mrmr_res)))
    sel_features <- rownames(rdr)[ix_features]
  }

  if (method == "glmnet-cox") {
    if (lambda == "min") {
      my_s <- "lambda.min"
    } else {
      my_s <- "lambda.1se"
    }
    if (length(alpha) > 0) {
      cvfit <- cv.glmnet(x = t(data),
                         y = surv_obj, family = "cox",
                         alpha = alpha)
    } else {
      cvfit <- cv.glmnet(x = t(data),
                         y = surv_obj,
                         family = "cox")
    }
    coef.min <- coef(cvfit, s = my_s)
    cf <- as.matrix(coef(cvfit, s = my_s))
    res <- cf[which(cf!=0), ]
    if (any(names(res) == "(Intercept)")) {
      res <- res[-which(names(res) == "(Intercept)")]
    }
    if (length(res) > 0) {
      message(paste("[RadAR]", "n=", length(res), "features selected with glmnet"))
    } else {
      assertthat::assert_that(1 < 0,
                              msg = "[RadAR] Error: No features selected with glmnet. Try changing lambda and/or alpha parameters.")
    }
    ix_features <- which(rownames(rdr) %in% names(res))
    sel_features <- names(res)
  }

  if (method == "glmnet-binomial") {
    if (lambda == "min") {
      my_s <- "lambda.min"
    } else {
      my_s <- "lambda.1se"
    }
    if (length(alpha) > 0) {
      cvfit <- cv.glmnet(x = t(data),
                         y = response,
                         family = "binomial",
                         alpha = alpha)
    } else {
      cvfit <- cv.glmnet(x = t(data),
                         y = response,
                         family = "binomial")
    }
    coef.min <- coef(cvfit, s = my_s)
    cf <- as.matrix(coef(cvfit, s = my_s))
    res <- cf[which(cf!=0), ]
    if (any(names(res) == "(Intercept)")) {
      res <- res[-which(names(res) == "(Intercept)")]
    }
    if (length(res) > 0) {
      message(paste("[RadAR]", "n=", length(res), "features selected with glmnet"))
    } else {
      assertthat::assert_that(1 < 0,
                              msg = "[RadAR] Error: No features selected with glmnet. Try changing lambda and/or alpha parameters.")
    }
    ix_features <- which(rownames(rdr) %in% names(res))
    sel_features <- names(res)
  }


  message(paste("[RadAR]", "Selected features are:", toString(sel_features)))
  rdr <- rdr[sel_features, ]
  metadata(rdr)$feature_selection_method <- method
  metadata(rdr)$feature_selection_selectby <- select_by
  metadata(rdr)$feature_selection_which_data <- which_data
  if (method %in% c("glmnet-cox", "glmnet-binomial")) {
    metadata(rdr)$glmnet_model <- cvfit
    metadata(rdr)$glmnet_model_lambda <- lambda
  } else {
    metadata(rdr)$glmnet_model <- F
    metadata(rdr)$glmnet_model_lambda <- NULL
  }
  out <- list()
  out$rdr <- rdr
  out$signature <- rownames(rdr)
  return(out)

}








