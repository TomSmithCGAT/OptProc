#' Plot the missing values.
#'
#' @param params Output from \code{svmOptimisation}
#' @return \code{data.frame()} with per organelle F1
#' @export
getOrganelleF1 <- function(params){

  # code below modified from Olly's
  times <- params@design[['times']]
  f1scoreperO <- matrix(0, times, nrow(params@cmMatrices[[1]]))
  
  for(i in 1:times){
    
    conf <- params@cmMatrices[[i]]
    
    f1perO <- MLInterfaces::F1(conf, naAs0 = TRUE)
    
    f1scoreperO[i, ] <- f1perO
  }
  
  OrganelleF1 <- data.frame(f1scoreperO)
  colnames(OrganelleF1) <- names(params@datasize$data.markers)
  
  return(OrganelleF1)
}
