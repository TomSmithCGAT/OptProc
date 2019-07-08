#' Plot the missing values.
#'
#' @param obj Input
#' @param verbose Tally the missing values
#' @param ... further arguments passed to gplots::heatmap.2
#' @seealso \code{\link{heatmap.2}} which this function wraps
#' @return NA
#' @export
plotMissing <- function(obj, verbose=TRUE, ...){
  tmp_obj <- obj
  exprs(tmp_obj)[!is.na(exprs(tmp_obj))] <- 1
  exprs(tmp_obj)[is.na(exprs(tmp_obj))] <- 0

  missing <- exprs(tmp_obj)
  missing <- missing[rowSums(missing==0)>0,] # identify features without missing values

  if(verbose){
    cat(sprintf("Out of %s total features, %s (%s%%) have missing values\n",
                length(rownames(exprs(tmp_obj))), length(rownames(missing)),
                round(100*length(rownames(missing))/length(rownames(exprs(tmp_obj))),3)))

    print(table(rowSums(missing==0)))
  }

  if(length(rownames(missing))>1){
    missing_twice <- missing[rowSums(missing==0)>1,]

    if(verbose){
      cat(sprintf("And %s features have more than one missing value\n", length(rownames(missing_twice))))
    }

    colnames(missing) <- pData(tmp_obj)$Sample_name

    gplots::heatmap.2(missing, col = c("lightgray", "black"),
              scale = "none", dendrogram = "none",
              trace = "none", keysize = 0.5, key = FALSE,Colv = F, labRow=F,
              ...)
  }

  return(NA)
}


#' Simulate missing values with mimimal assumptions
#'
#' Assumes that PSMs from the same peptide should have the same
#' relative abundances, and that the missing values are missing at random (MAR)
#' with peptide intensity the main cause of missingness.
#'
#' Identifies groups of PSMs with at least one with missing values (PSMm) and
#' one without (PSMnm) from the same peptide. PSMnm are then randomly paired
#' with a PSMm, intensity values scaled to the PSMm and missing values
#' (\code{NA}) added at the same positions as missing values in PSMm.
#'
#' Returns two MSnSets: "truth" and "missing". Both have the same dimensions as
#' the input MSnSet. Truth contains the scaled PSMs. Missing contains the scaled
#' PSMs with missing values added.
#'
#'
#' @param obj Input MSnSet
#' @param n Number of PSMs to impute missing values in
#' @param id_column Column to group PSM by
#' @param verbose Output description of missing values to console
#' @return list("truth":Input data set with ground truths for missing values,
#'              "missing":Input data set with NA for missing values)
#' @export
#'
addMissing <- function(obj, n=Inf, id_column, verbose=TRUE){

  # identify ids with missing values
  missing <- obj[rowSums(is.na(obj))>0,]
  id_missing <- fData(missing)[[id_column]]

  # identify ids without missing values
  not_missing <- obj[rowSums(is.na(obj))==0,]
  id_not_missing <- fData(not_missing)[[id_column]]

  # identify ids with both missing and no missing
  id_both <- sample(intersect(id_missing, id_not_missing))

  # how many features without missing values for these ids
  n_features <- nrow(not_missing[fData(not_missing)[[id_column]] %in% id_both,])

  if(verbose){
    message(
      sprintf("%s features where missing values can be simulated.
These comprise %s unique meta-features\n", n_features, length(id_both)))
  }

  # create new array matrices to hold the truth and missing quantification values
  truth_e <- missing_e <- exprs(obj)


  missing_n <- 0
  for(id in id_both){

    if(missing_n == n) break

    .e <- exprs(obj)[fData(obj)[[id_column]]==id,]
    .e_m <- .e[rowSums(is.na(.e))>0,,drop=FALSE]
    .e_nm <- .e[rowSums(is.na(.e))==0,,drop=FALSE]

    feature_nm <- rownames(.e_nm)
    feature_m <- rownames(.e_m)

    rand_m <- sample(feature_m,
                     size=min(length(feature_nm), length(feature_m)),
                     replace=FALSE)

    feature_nm <- sample(feature_nm)[1:length(rand_m)]


    for(ix in seq_along(feature_nm)){

      if(missing_n == n) break

      id_nm <- feature_nm[ix]
      id_m <- rand_m[ix]

      ave_e_m <- sum(.e_m[id_m,], na.rm=TRUE)
      ave_e_nm <- sum(.e_nm[id_nm,])
      intensity_ratio <- ave_e_nm/ave_e_m

      reduced_values <- .e_nm[id_nm,]/intensity_ratio
      reduced_values_add_m <- reduced_values
      reduced_values_add_m[is.na(.e_m[id_m,])] <- NA

      truth_e[id_nm, ] <- reduced_values
      missing_e[id_nm, ] <- reduced_values_add_m

      missing_n <- missing_n + 1
    }
  }

  truth <- missing <- obj
  exprs(truth) <- truth_e
  exprs(missing) <- missing_e

  fData(missing)$sim_missing <- (
    rowSums(is.na(truth_e))==0 &
      rowSums(is.na(missing_e))>0)

  fData(truth)$sim_missing <- (
    rowSums(is.na(truth_e))==0 &
      rowSums(is.na(missing_e))>0)

  if(verbose){
    message(sprintf("%s/%s features have had missing values simulated",
                    sum(fData(missing)$sim_missing),
                    nrow(missing)))
  }

  invisible(list("truth"=truth, "missing"=missing))
}
