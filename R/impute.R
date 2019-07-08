#' Impute missing values using compositional KNN imputation
#'
#' @param obj Input
#' @param method impKNNa=compositional KNN or
#'               impCoda=compositional KNN with iterative improvement in estimates
#' @param k nearest neighbours
#' @param sum_first Sum normalise the array data before imputation, and convert
#'   back to non-sum normalised afterwards
#' @param verbose Log the start and end time for imputation
#' @seealso \code{\link{robCompositions::impKNNa}} and
#' \code{\link{robCompositions::impCoda}} which this function wraps
#' @return Input data with missing values imputed
#' @export
robComp <- function(obj, method="comp-knn", k, sum_first=FALSE, verbose=FALSE){

  if(verbose){
    cat(sprintf("Starting Compositions KNN imputation at: %s\n", Sys.time()))
  }

  if(sum_first){
    obj_e <- obj %>% normalise("sum") %>% exprs()
  } else {
    obj_e <- obj %>% exprs()
  }

  if(method=="comp-knn"){
    eImp <- robCompositions::impKNNa(obj_e, k=k)
  } else if(method=="it-comp-knn"){
    eImp <- robCompositions::impCoda(obj_e, method='ltsReg', k=5)
  } else { stop("Unknown method")}

  if(verbose){
    cat(sprintf("Finished imputation at: %s\n", Sys.time()))
  }

  .imputed <- obj
  exprs(.imputed) <- eImp$xImp
  exprs(.imputed) <- exprs(.imputed) * rowSums(exprs(obj), na.rm=TRUE)

  invisible(.imputed)
}

#' Impute missing values using KNN imputation with prior sum normalisation
#'
#' @param obj Input
#' @param k nearest neighbours
#' @param verbose Log the start and end time for imputation
#' @return Input data with missing values imputed
#' @export
SnKnn <- function(obj, k=5, verbose=FALSE){

  sum_intensity <- rowSums(exprs(obj), na.rm=TRUE)

  if(!args$verbose){
    sink("/dev/null")
  }

  sn_knn <- obj %>%
    normalise(method="sum") %>%
    MSnbase::impute(method="knn", k=k)

  if(!args$verbose){
    sink()
  }

  exprs(sn_knn) <- exprs(sn_knn) * sum_intensity

  return(sn_knn)

}

#' Extract the truth and imputed values
#'
#' @param obj MSnSet with imputed values
#' @param missing_obj MSnSet with missing values
#' @param truth_obj MSnSet with ground truths
#' @return data.frame() with truth and imputed columns
#' @export
getTruthImputed <- function(obj, missing_obj, truth_obj){

  missing <- missing_obj[fData(missing_obj)$sim_missing,]
  truth <- truth_obj[fData(missing_obj)$sim_missing,]
  obj_imp_only <- obj[fData(missing_obj)$sim_missing,]

  missing_ix <- is.na(exprs(missing))
  missing_truth <- exprs(truth)[missing_ix]
  missing_imputed <- exprs(obj_imp_only)[missing_ix]

  return(data.frame("truth"=missing_truth, "imputed"=missing_imputed))
}

#' Obtain the Root mean squared error (RMSE)
#' @param obj data.frame() with truth and imputed columns
#' @return RMSE
#' @export
#' @export
#' .data <- data.frame("truth"=c(1,2,3,4,5), "imputed"=c(1,3,2,4,6))
#' print(getRMSE(.data))
getRMSE <- function(obj){
  RMSE <- sqrt(mean((obj$truth - obj$imputed)^2))
  invisible(RMSE)
}


#' Plot the truth vs imputed
#' @param obj data.frame() with truth and imputed columns
#' @return scatter plot for truth vs imputed
#' @export
plotTruthImputed <- function(obj){
  p <- ggplot(obj, aes(log2(truth), log2(imputed))) +
    geom_point(size=0.25, alpha=0.25) +
    my_theme +
    theme(text=element_text(size=14),
          plot.title=element_text(size=10, hjust=0.5)) +
    geom_abline(slope=1, linetype="dashed", colour=cbPalette[7])

  return(p)
}
