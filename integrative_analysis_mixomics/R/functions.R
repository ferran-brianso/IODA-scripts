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


## A PARTIR D'AQUI, DEFINIR LES FUNCIONS QUE EXECUTIN 
## L'ANALISI INTEGRATIVA SEGONS EL METODE (I EL PAQUET)
## QUE ES VULGUI APLICAR:

perform_rCCA <- function(X, Y, resultsDir, scoresFile, rccFile){
  ##
  ##
  require(mixOmics)
  ## Run rCCA given the provided parameters:
  load(file=file.path(resultsDir, scoresFile))
  lambda1 <- cv.score$opt.lambda1
  lambda2 <- cv.score$opt.lambda2
  result <- rcc(X, Y, ncomp = 3, lambda1 = lambda1, lambda2 = lambda2)
  #head(result$cor)
  #head(result$loadings$X)
  #head(result$loadings$Y)
  save(result, file=file.path(resultsDir, rccFile))
  return(result)
}

plot_indiv <- function(rccaResult, resultsDir){
  ##
  ##
  require(mixOmics)
  ## Get sample labels from original data and plot individual samples
  data("breast.TCGA")
  #dim(breast.TCGA$data.train$mrna)
  
  col.groups <- as.numeric(breast.TCGA$data.train$subtype)
  #col.groups

  plotIndiv(rccaResult, comp = 1:2, col = col.groups, cex=3)
  ## Save individual samples plot into a file
  png(file=file.path(resultsDir, "indivSamples.png"), 
      width=800, height=800, res=120)
    plotIndiv(rccaResult, comp = 1:2, col = col.groups, cex=3)
  dev.off()  
}



plot_corrCirc <- function(rccaResult, resultsDir, 
                          comp = 1:2, cutOff=0.6, cex = c(3, 3)){
  ##
  ##
  require(mixOmics)
  ## PLor correlation circle plot at a given correlation cutoff

  plotVar(rccaResult, comp = 1:2, cutoff=cutOff, cex = c(3, 3))
  ## Save corr circle plot into a file
  png(file=file.path(resultsDir, paste0("corrCircle.",cutOff,".png")), 
      width=800, height=800, res=120)
    plotVar(rccaResult, comp = 1:2, cutoff=cutOff, cex = c(3, 3))
  dev.off()
}

