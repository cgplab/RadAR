

# Detect and replace outliers based on IQR method.

#' Find and replace radiomic feature outliers based on IQR.
#'
#' Feature values below Q1-1.5 IQR or above Q3+1.5 IQR are
#' identified as outliers.
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}})
#' @param replace_with (character) Replace outliers with one of the following: "NA", "mean", "median"
#'
#' @return An updated rdr (a RadAR object)
#' @export
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @examples
find_feature_outliers <- function(rdr = NULL,
                                  replace_with="NA")
{

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] rdr object required")
  assertthat::assert_that(replace_with %in% c("NA", "mean", "median"), msg = "[RadAR] Invalid replacement rule for outliers")
  method <- "IQR"
  data <- rdr@assays$data$values
  mean_data <- apply(data, 1, mean, na.rm = T)
  median_data <- apply(data, 1, median, na.rm = T)
  tot_outl <- c()
  if (method == "IQR") {
    for ( i in 1: nrow(data)) {
      outl <- boxplot.stats(data[i, ])$out
      if (replace_with == "NA") {
        data[i, names(outl)] <- NA
      }
      if (replace_with == "mean") {
        data[i, names(outl)] <- mean_data[i]
      }
      if (replace_with == "median") {
        data[i, names(outl)] <- median_data[i]
      }
      tot_outl <- c( tot_outl, outl)
    }
    message(paste("[RadAR]", length(tot_outl), "outliers found."))
  }
  assay(rdr) <- data
  return(rdr)
}

#' Filter out oulier patients
#'
#' Identify feature outliers based oin IQR and exclude patients having a number of feature outliers
#' greater than or equal to a predifined thereshold.
#' This function is usefult to exclude heterogeneous patients and make the dataset more homogenous.
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}})
#' @param fraction_outliers (numeric) Minimum fraction of feature outliers to exclude a patient.
#'
#' @return An updated, eventually reduced, rdr (a RadAR object)
#' @export
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @examples
find_outliers <- function(rdr = NULL,
                          fraction_outliers=NULL)
{

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] rdr object required")
  assertthat::assert_that(length(fraction_outliers) > 0, msg = "[RadAR] fraction_outliers parameter should be defined")
  assertthat::assert_that(fraction_outliers > 0 & fraction_outliers <= 1, msg = "[RadAR] fraction_outliers parameter should be in the range (0, 1]")
  method <- "IQR"
  data <- rdr@assays$data$values
  outl <- vector("list", length = nrow(data))
  outl_samples <- c()
  if (method == "IQR") {
    for ( i in 1: nrow(data)) {
      outl[[i]] <- boxplot.stats(data[i, ])$out
      outl_samples <- c(outl_samples, names(outl[[i]]))
    }
  }
  stats_outl_samples <- sort(table(outl_samples), decreasing = T)/nrow(data)
  outl_samples <- names(which(stats_outl_samples >= fraction_outliers))
  nooutl_samples <- names(which(stats_outl_samples < fraction_outliers))
  message(paste("[RadAR]", length(outl_samples), "samples were removed as those have fraction_outliers>=", fraction_outliers,". \n Samples include the following:", toString(outl_samples)))
  rdr_out <- rdr[,nooutl_samples]
  return(rdr_out)
}

# Normalize feature values, only quntile norm implemented so far --------------------------------------
#' Apply normalization to feature values
#'
#' Feature values are quantile normalized by limma package.
#' After quantile normalization, the distributions of feature values are identical in terms of
#' statistical properties.
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}})
#' @param method Only quantile normalization has been implemented so far. After quantile normalization,
#' the distributions of feature values are identical in terms of statistical properties
#' @param which_data (character) Which data should be normalized. It can be one of "normal" or "scaled".
#'
#' @return An updated rdr (a RadAR object)
#' @export
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @examples
normalize_feature_values <- function(rdr = NULL,
                          method = "quantile",
                          which_data = "scaled")
{

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] rdr object required")
  assertthat::assert_that(method %in% c("quantile"),
                          msg = "[RadAR] Invalid method for values normalization")
  assertthat::assert_that(which_data %in% c("normal", "scaled"),
                          msg = "[RadAR] Invalid value for variable which_data")
  if (which_data == "normal") {
    data <- rdr@assays$data$values
  }
  if (which_data == "scaled") {
    data <- rdr@assays$data$scaled_values
  }
  data <- t(limma::normalizeQuantiles(t(data)))
  rdr@assays$data$norm_values <- data
  return(rdr)
}

# Scale feature values --------------------------------------------------------------------
#' This function implements different scaling strategies to feature values
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}})
#' @param method (character) Which method use for scaling feature values. It can be one of "minmax" or "median".
#' Using min-max, all features will take values in the range [0,1].
#' Using median, all features will be median subctracted.
#'
#' @return An updated rdr (a RadAR object)
#' @export
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @examples
scale_feature_values <- function(rdr = NULL,
                                 method = "minmax")
{

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(method %in% c("minmax", "median"),
                          msg = "[RadAR] Error: Invalid method for scaling feature values")
  data <- rdr@assays$data$values
  if (method == "minmax") {
    range_features <- apply(data, 1, range, na.rm = T)
    data <- (data - range_features[1,])/ ( range_features[2,] -  range_features[1,])
  }
  if (method == "median") {
    median_data <- apply(data, 1, median, na.rm = T)
    data <- data - median_data
  }
  rdr@assays$data$scaled_values <- data
  return(rdr)
}

#' Filter a RadAR object by image type(s)
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}})
#' @param image_type (character) Only features corresponding to image_type will be retained.
#' Available values for image_type depend on the software used for feature extraction.
#' More than one image_type can be passed. Use \code{print_image_type} to visualize available image_type.
#'
#' @return An updated rdr (a RadAR object)
#' @export
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @examples
filter_by_image_type <- function(rdr = NULL,
                                 image_type = NULL) {

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(length(image_type) > 0, msg = "[RadAR] Error: image_type values required")
  assertthat::assert_that(all(image_type %in% rowData(rdr)$image_type),
                          msg = "[RadAR] Error: image_type doesn't match available image types")
  ix_image_types <- which(rowData(rdr)$image_type %in% image_type)
  rdr <- rdr[ix_image_types, ]
  return(rdr)
}


#' Filter a RadAR object by feature type(s)
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}})
#' @param image_type (character) Only features corresponding to feature_type will be retained.
#' Available values for feature_type depend on the software used for feature extraction.
#' More than one feature_type can be passed. Use \code{print_feature_type} to visualize available feature_type.
#'
#' @return An updated rdr (a RadAR object)
#' @export
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @examples
#'
filter_by_feature_type <- function(rdr = NULL,
                                   feature_type = NULL) {

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(length(feature_type) > 0, msg = "[RadAR] Error: feature_type values required")
  assertthat::assert_that(all(feature_type %in% rowData(rdr)$feature_type),
                          msg = "[RadAR] Error: feature_type doesn't match available feature types")
  ix_feature_types <- which(rowData(rdr)$feature_type %in% feature_type)
  rdr <- rdr[ix_feature_types, ]
  return(rdr)
}

# Find clusters   --------------------------------------------------------------------
#' Flag samples and/or features by cluster membership
#'
#' This function assigns cluster membership to patients/ROIs and/or features
#' based on previous hierarchical clustering (\code{\link{do_hierarchical_clustering}})
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiment}})
#' @param n_clusters_cols (numeric) Number of clusters for patients/ROIs (columns)
#' @param n_clusters_features (numeric) Number of clusters for features (rows)
#'
#' @return An updated rdr (a RadAR object)
#' @export
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @examples
find_clusters <- function(rdr = NULL,
                          n_clusters_cols = NULL,
                          n_clusters_features = NULL)
{
  if (is.null(n_clusters_cols)) {
    n_clusters_cols <- -1
  }
  if (is.null(n_clusters_features)) {
    n_clusters_features <- -1
  }
  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(n_clusters_cols > 0 | n_clusters_features > 0,
                          msg = "[RadAR] Error: n_clusters_cols and/or n_clusters_features need to be specified")
  assertthat::assert_that(length(metadata(rdr)$hcl) > 0, msg = "[RadAR] Error: hierachical clustering has not been yet computed")
  hcl <- metadata(rdr)$hcl
  if (n_clusters_cols > 0) {
    clusts <- cutree(hcl$hcl_col, k = n_clusters_cols)
    colData(rdr)$hcl_cols <- clusts
  }
  if (n_clusters_features > 0) {
    clusts <- cutree(hcl$hcl_row, k = n_clusters_features)
    rowData(rdr)$hcl_features <- clusts
  }
  return(rdr)
}




