#' Obtain the QSep distances for an MSnSet with marker proteins annotated
#' @param obj MSnSet input
#' @param ... further arguments passed to QSep
#' @return \code{data.frame()} with qsep distances
#' @export
getQsepDistances <- function(obj,  ...){
  suppressWarnings(suppressPackageStartupMessages(library("pRoloc")))
  qsep_dist <- qsep(QSep(obj, ...))
  clustering <- hclust(as.dist(qsep_dist))
  
  qsep_df <-  cbind(expand.grid(
    "other" = rownames(qsep_dist),
    "reference" = colnames(qsep_dist),
    stringsAsFactors=FALSE),
    "qsep" = as.vector(qsep_dist)) %>%
    data.frame() %>%
    mutate(other=factor(other, levels=clustering$labels[clustering$order]),
           reference=factor(reference, levels=clustering$labels[clustering$order]))

  invisible(qsep_df)
}

#' Plot the QSep distances obtained from \code{getQsepDistances}
#' @param qsep_df Output from \code{getQsepDistances}
#' @return \code{ggplot2} plot
#' @export
plotQsepDistances <- function(qsep_df){
    
    p <- ggplot(qsep_df, aes(reference, other, fill=qsep)) +
      geom_tile() +
      my_theme +
      theme(text=element_text(size=10),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
      xlab("Reference") +
      ylab("") +
      scale_fill_continuous(low="white", high="steelblue", name="QSep\ndistance")
    
    invisible(p)
}

