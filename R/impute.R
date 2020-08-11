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
    eImp <- robCompositions::impCoda(obj_e, method='ltsReg', k=k)
  } else { stop("Unknown method")}

  if(verbose){
    cat(sprintf("Finished imputation at: %s\n", Sys.time()))
  }

  .imputed <- obj
  exprs(.imputed) <- as.matrix(eImp$xImp)

  if(sum_first){
    exprs(.imputed) <- exprs(.imputed) * rowSums(exprs(obj), na.rm=TRUE)
  }

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

  if(!verbose){
    sink("/dev/null")
  }

  sn_knn <- obj %>%
    normalise(method="sum") %>%
    MSnbase::impute(method="knn", k=k)

  if(!verbose){
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

  missing_obj <- missing_obj[rownames(obj),]
  truth_obj <- truth_obj[rownames(obj),]

  missing <- missing_obj[fData(missing_obj)$sim_missing,]
  truth <- truth_obj[fData(missing_obj)$sim_missing,]
  obj_imp_only <- obj[fData(missing_obj)$sim_missing,]

  missing_ix <- is.na(exprs(missing))

  missing_truth <- exprs(truth)[missing_ix]

  missing_imputed <- exprs(obj_imp_only)[missing_ix]

  return(data.frame("truth"=missing_truth, "imputed"=missing_imputed))
}

#' Obtain the Root mean squared error (RMSE)
#' @param truth ground truths
#' @param imputed imputed values
#' @return RMSE
#' @export
#' truth <- c(1,2,3,4,5)
#' imputed <- c(1,3,2,4,6))
#' print(getRMSE(truth, imputed))
getRMSE <- function(truth, imputed){
  RMSE <- sqrt(mean((truth - imputed)^2))
  invisible(RMSE)
}


#' Obtain the Root median squared error (RMedSE)
#' @param truth ground truths
#' @param imputed imputed values
#' @return RMedSE
#' @export
#' truth <- c(1,2,3,4,5)
#' imputed <- c(1,3,2,4,6)
#' print(getRMedSE(truth, imputed))
getRMedSE <- function(truth, imputed){
  RMedSE <- sqrt(median((truth - imputed)^2))
  invisible(RMedSE)
}

#' Obtain the Mean Absolute Relative Difference (MARD)
#' @param truth ground truths
#' @param imputed imputed values
#' @return MARD
#' @export
#' truth <- c(1,2,3,4,5)
#' imputed <- c(1,3,2,4,6)
#' print(getMARD(truth, imputed))
getMARD <- function(truth, imputed){
  MARD <- mean((abs(truth - imputed)/truth))
  invisible(MARD)
}


#' Obtain the Mean Absolute Log2 Quotient Error (MAQE)
#' @param truth ground truths
#' @param imputed imputed values
#' @return MAQE
#' @export
#' truth <- c(1,2,3,4,5)
#' imputed <- c(1,3,2,4,6)
#' print(getMAQE(truth, imputed))
getMAQE <- function(truth, imputed){
  MAQE <- mean(abs(log2(truth/imputed)))
  invisible(MAQE)
}

#' Plot the truth vs imputed
#' @param obj data.frame() with truth and imputed columns
#' @return scatter plot for truth vs imputed
#' @export
plotTruthImputed <- function(obj){

  min_value = min(c(log2(obj$truth), log2(obj$imputed)), na.rm=TRUE)

  max_value = max(c(log2(obj$truth), log2(obj$imputed)), na.rm=TRUE)

  p <- ggplot(obj, aes(log2(truth), log2(imputed))) +
    geom_point(size=0.25, alpha=0.25) +
    ggtitle(round(getMAQE(obj$truth, obj$imputed), 1)) +
    my_theme +
    theme(text=element_text(size=20),
          plot.title=element_text(size=20, hjust=0.5)) +
    geom_abline(slope=1, linetype="dashed", colour=cbPalette[7]) +
    xlim(min_value, max_value) +
    ylim(min_value, max_value)

  return(p)
}

#' Impute missing values using either \code{\link{MSnbase::impute}} or
#' \code{\link{RobCompositions}}
#' @param obj MSnSet with missing values
#' @param method imputation method, Choose from knn, MinDet, MinProb,
#' sn-knn, it-comp-knn or comp-knn
#' @return input MSnSet with missing values imputed
#' @seealso \code{\link{robCompositions::impKNNa}} for comp-knn and
#' \code{\link{robCompositions::impCoda}} for it-comp-knn methods
#' @export
imputeOptProc <- function(obj, method, k, verbose=FALSE){
  if(verbose) write("Imputing missing values", stdout())
  msnbase_methods <- c("knn", "MinDet", "MinProb")
  allowed_methods <- c(msnbase_methods, "sn-knn","it-comp-knn", "comp-knn")

  if(!method %in% allowed_methods){
    stop(sprintf("method must be in %s", paste(allowed_methods, collapse=",")))
  }

  if(method %in% msnbase_methods){
    if(!verbose){
      sink("/dev/null")
    }

    if(method=="knn"){
      .imputed <- obj %>% MSnbase::impute(method="knn", k=k)
    } else if(method=='MinProb') {
      .imputed <- obj %>% MSnbase::impute(method='MinProb', tune.sigma=0.01)
    } else {
      .imputed <- obj %>% MSnbase::impute(method=method)
    }

    if(!verbose){
      sink()
    }

  } else if(method=="sn-knn"){
    .imputed <- obj %>%
      SnKnn(k=k, verbose=verbose)
  } else {
    .imputed <- obj %>%
      robComp(method=method, k=k, sum_first=TRUE, verbose=TRUE)
  }

  return(.imputed)
}
