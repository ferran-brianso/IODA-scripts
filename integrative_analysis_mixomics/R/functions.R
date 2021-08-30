#' @title Functions to run the "IODA (Mixomics part)" in targets' mode
#' @description .......
#' @return ..........
#' @param .......
#' @examples
#' library(ggplot2)
#' library(tidyverse)
#' library(stats)


load_and_check <- function(inFile, type = "G") {
  ## expected input: path to raw data file in .Rda format 
  ##                 which should have gene.data and prot.data objects!!
  ## RETURNS: gene or protein data (checked in order to remove duplicates)

  load(file=inFile, verbose = TRUE)
  
  if (type=="G"){
    in.data <- t(gene.data)
  }else{
    in.data <- t(prot.data)
  }
  #dim(in.data)
  #sum(is.na(in.data))
  ## check if there are some duplicated gene names
  dups <- rownames(in.data)[which(duplicated(rownames(in.data)))]
  if (length(dups)>0){
    in.data[which(rownames(in.data) %in% dups),]
    in.data <- in.data[!duplicated(rownames(in.data)),]
  }
  #dim(in.data)
  return(in.data)
}


plot_corr_matrix <- function(X, Y, resultsDir){ 
  ## 
  ## 
  require(mixOmics)
  ## Save plot out into a file
  png(file=file.path(resultsDir, "corrMatrix.png"), 
      width=800, height=800, res=120)
  imgCor(X, Y)
  dev.off()
}