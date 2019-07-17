#' Combine TMT sets.
#' 
#' Note that only intersecting features are returned
#'
#' @param x MSnSet containing quanification values for TMT set 1
#' @param y MSnSet containing quanification values for TMT set 2
#' @return Combined MsnSet
#' @export
combineTMTsets <- function(x, y){
  library(MSnbase)
  
  if(any(rownames(x) %in% rownames(y))){
    x <- updateSampleNames(x, "TMT_1")
    y <- updateSampleNames(y, "TMT_2")
  }

  # on the reasonable assumption that the features names will be overlapping
  y <- updateFvarLabels(y, "TMT_2")
  
  intersecting_peptides <- intersect(rownames(x), rownames(y))
  
  combined_quant <- MSnbase::combine(x[intersecting_peptides,],
                                     y[intersecting_peptides,])
  
  return(combined_quant)
}

#' Adjust quantification values using bridging channels
#' 
#' Adjust quantification values using bridging channels to
#' normalise across isobaric tag multiplexes.
#' Note that the bridge channels are not normalised.
#'
#' @param values vector of quantification values
#' @param plex1_ix index for quantification values in multiplex 1
#' @param plex2_ix index for quantification values in multiplex 2
#' @param bridge1_ix index for bridge sample(s) in multiplex 1
#' @param bridge2_ix index for bridge sample(s) in multiplex 1
#' @return bridge normalised values
#' @examples
#' bn(c(1,1,2,2,2,4), plex1_ix=1:2, plex2_ix=4:5, bridge1_ix=3, bridge2_ix=6)
bn <- function(values, plex1_ix, plex2_ix, bridge1_ix, bridge2_ix){
  bridge_1 <- values[bridge1_ix]
  bridge_2 <- values[bridge2_ix]
  if(any(is.na(c(bridge_1, bridge_2)))){
    values[plex1_ix] <- NA
    values[plex2_ix] <- NA
  }
  
  else{
    
    bridge_ratios <- bridge_1/bridge_2
    print(bridge_ratios)
    bridge_mean_ratio = mean(bridge_ratios)
    print(bridge_mean_ratio)
    print(c(plex2_ix, bridge2_ix))
    values[c(plex2_ix, bridge2_ix)] <- (
      values[c(plex2_ix, bridge2_ix)] * bridge_mean_ratio)
  }
  return(values)
}

#' Adjust quantification values in MSnSet
#' 
#' Adjust quantification values in MSnSet using bridging channels to
#' normalise across isobaric tag multiplexes.
#' Note that the bridge channels are not normalised.
#'
#' @param values vector of quantification values
#' @param plex1_ix index for quantification values in multiplex 1
#' @param plex2_ix index for quantification values in multiplex 2
#' @param bridge1_ix index for bridge sample(s) in multiplex 1
#' @param bridge2_ix index for bridge sample(s) in multiplex 1
#' @param keep_bridge How many values to retain for each bridge sample (0, 1, or 2)
#' @return bridge normalised values
#' @export
bridgeNormalise <- function(obj, plex1_ix=1:10, plex2_ix=12:21,
                            bridge1_ix=11, bridge2_ix=22, 
                            keep_bridge=0){
  return_obj <- obj
  exprs(return_obj) <- exprs(return_obj) %>%
    apply(1, FUN=function(x) bn(x, plex1_ix, plex2_ix, bridge1_ix, bridge2_ix)) %>%
    t()
  
  if(keep_bridge==0){
    invisible(return_obj[,c(plex1_ix, plex2_ix)])
  } else if(keep_bridge==1){
    invisible(return_obj[,c(plex1_ix, bridge1_ix, plex2_ix)])
  } else if(keep_bridge==2){
    invisible(return_obj)
  } else{ stop("keep_bridge must be 0, 1 or 2")}
}

#' Average the bridge sample(s) in a vector of quantification values
#' 
#' Note, redundant bridge sample(s) are removed
#'
#' @param values vector of quantification values
#' @param plex1_ix index for quantification values in multiplex 1
#' @param plex2_ix index for quantification values in multiplex 2
#' @param bridge1_ix index for bridge sample(s) in multiplex 1
#' @param bridge2_ix index for bridge sample(s) in multiplex 1
#' @return bridge normalised values
#' @examples
#' q <- c(1,1,NA,3,2,2,4,6)
#' MeanBridges(q, plex1_ix=1:2, plex2_ix=5:6, bridge1_ix=3:4, bridge2_ix=7:8)
mb <- function(values, plex1_ix, plex2_ix, bridge1_ix, bridge2_ix){
  bridge_1 <- values[bridge1_ix]
  bridge_2 <- values[bridge2_ix]
  
  bridge_means <- colMeans(rbind(bridge_1, bridge_2), na.rm=TRUE)
  
  values[bridge1_ix] <- bridge_means
  return(values[c(plex1_ix, bridge1_ix, plex2_ix)])
}


#' Average the bridge sample(s) for each feature in an MSnSet
#' 
#' Note, redundant bridge sample(s) are removed
#'
#' @param values vector of quantification values
#' @param plex1_ix index for quantification values in multiplex 1
#' @param plex2_ix index for quantification values in multiplex 2
#' @param bridge1_ix index for bridge sample(s) in multiplex 1
#' @param bridge2_ix index for bridge sample(s) in multiplex 1
#' @return bridge normalised values
#' @export
MeanBridges <- function(obj, plex1_ix, plex2_ix, bridge1_ix, bridge2_ix){
  return_obj <- obj[,c(plex1_ix, bridge1_ix, plex2_ix)]
  exprs(return_obj) <- exprs(obj) %>%
    apply(1, FUN=function(x) mb(x, plex1_ix, plex2_ix, bridge1_ix, bridge2_ix)) %>%
    t()
  
  invisible(return_obj)
}