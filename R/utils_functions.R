#' Print available image types
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiments}})
#'
#' @return A screen message showing available image types.
#' @export
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @examples
print_image_type <- function(rdr = NULL) {

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  available <- toString(unique(rowData(rdr)$image_type))
  available <- ifelse(available == "", "none available", available)
  message(paste0("[RadAR] ", "Available image types are\n", gsub(", ", "\n", available)))

}

#' Print available feature types
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiments}})
#'
#' @return A screen message showing available feature types.
#' @export
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @examples
print_feature_type <- function(rdr = NULL) {

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  available <- toString(unique(rowData(rdr)$feature_type))
  available <- ifelse(available == "", "none available", available)
  message(paste0("[RadAR] ", "Available feature types are\n", gsub(", ", "\n", available)))

}

#' Print available methods for hierarchical clustering
#'
#' @return A screen message with available methods.
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
print_hcl_methods <- function() {

  methods_hcl <- c("ward.D",
                   "ward.D2",
                   "single",
                   "complete",
                   "average" ,
                   "mcquitty",
                   "median",
                   "centroid")

  message(paste0("[RadAR] ", "Available hiearchical clustering methods are\n", gsub(", ", "\n", toString(methods_hcl))))

}

#' Print available methods for computing distance
#'
#' @return A screen message with available methods.
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
print_distance_methods <- function() {

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

  message(paste0("[RadAR] ", "Available distance methods are\n", gsub(", ", "\n", toString(methods_dist))))

}

#' Retrieve the description of a list of input features
#'
#' @param rdr A RadAR object (class \code{\link{SummarizedExperiments}})
#' @param feature_names (character) names of input features.
#'
#' @return A table of class \code{\link{tibble}} reporting feature names and dictionary.
#' @author Matteo Benelli (\email{matteo.benelli@uslcentro.toscana.it})
#' @export
#'
#' @examples
#'
features_to_dictionary <- function(rdr = NULL,
                                   feature_names = NULL) {

  assertthat::assert_that(length(rdr) > 0, msg = "[RadAR] Error: rdr object required")
  assertthat::assert_that(length(feature_names) > 0, msg = "[RadAR] Error: feature_names required")
  assertthat::assert_that(all(feature_names %in% rowData(rdr)$feature_name),
                              msg = paste("[RadAR] Error: the following features are not in rdr:",
                                          toString(feature_names[which(feature_names %in% rowData(rdr)$feature_name == F)])))

  out <- cbind(feature_names,
               rowData(rdr)[feature_names, ]$feature_type,
               rowData(rdr)[feature_names, ]$feature_description
               )
  colnames(out) <- c("feature_name", "feature_type", "description")
  out <- tibble::as_tibble(out)
  return(out)
}




