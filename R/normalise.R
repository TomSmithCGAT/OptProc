#' Normalise data using sum +/- VSN
#'
#' @param obj Input
#' @param method sum = row-wise sum normalisation, or
#'               vsn-sum = VSN then row-wise sum normalisation, or
#'               sum-vsn = row-wise sum normalisation then VSN
#' @return Normalised Input data 
#' @export
spatial_proteomics_normalise <- function(obj, method="sum"){
  
  if(method=="sum"){
    .norm <- obj %>% normalise('sum')
  } else if(method=="vsn-sum"){
    .norm <- obj %>% normalise('vsn') %>% normalise('sum')
  } else if(method=="sum-vsn"){
    .norm <- obj %>% normalise('sum') %>% normalise('vsn')
  } else { stop("Unknown method")}
  
  invisible(.norm)
}