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
  ##                 which should have data.1 (e.g. genes) and data.2 (e.g. prots) objects!!
  ## RETURNS: data1 or data2 loaded and checked (in order to remove duplicates)

  load(file=inFile, verbose = TRUE)
  
  if (type=="G"){
    in.data <- data.1 #t(data.1)
  }else{
    in.data <- data.2 #t(data.2)
  }
  #dim(in.data)
  #sum(is.na(in.data))
  ## check if there are some duplicated gene names (rows)
  dups <- rownames(in.data)[which(duplicated(rownames(in.data)))]
  if (length(dups)>0){
    in.data[which(rownames(in.data) %in% dups),]
    in.data <- in.data[!duplicated(rownames(in.data)),]
  }
  # ## and check if there are some duplicated sample names (columns)
  # dups2 <- colnames(in.data)[which(duplicated(colnames(in.data)))]
  # if (length(dups2)>0){
  #   in.data[which(colnames(in.data) %in% dups2),]
  #   in.data <- in.data[!duplicated(colnames(in.data)),]
  # }
  
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
  
  ## Plot correlation circle at a given correlation cutoff
  plotVar(rccaResult, comp = 1:2, cutoff=cutOff, cex = c(3, 3))
  ## Save corr circle plot into a file
  png(file=file.path(resultsDir, paste0("corrCircle.",cutOff,".png")), 
      width=800, height=800, res=120)
    plotVar(rccaResult, comp = 1:2, cutoff=cutOff, cex = c(3, 3))
  dev.off()
}



tag_features <- function(rccaResult, x.tag = "x.", y.tag = "y."){
  ##
  ##
  
  ## Prepare rCCA result in order to avoid having duplicated vertex names in the resulting network
  ## That means adding a  'x' or 'y' to each of the feature names in the rcc result object
  colnames(rccaResult$X) <- paste0(x.tag, colnames(rccaResult$X))
  colnames(rccaResult$X)
  rccaResult$names$colnames$X <- colnames(rccaResult$X)

  colnames(rccaResult$Y) <- paste0(y.tag, colnames(rccaResult$Y))
  colnames(rccaResult$Y)
  rccaResult$names$colnames$Y <- colnames(rccaResult$Y)
  return(rccaResult)
}



get_corrNetwork <- function(rccaResult, resultsDir, netw.threshold = 0.5){
  ##
  ##
  require(mixOmics)
  require(igraph)
  
  ## Create the network object
  net <- network(rccaResult, comp = 1:3, interactive = FALSE, cutoff=netw.threshold)
  dev.off()

  ## Save relevance network into a png file
  png(file=file.path(resultsDir, 
                     paste0("relNetwork", as.character(netw.threshold), ".png")),
      width=800, height=800, res=120)
    network(rccaResult, comp = 1:3, interactive = FALSE, cutoff=netw.threshold)
  dev.off()
  
  ## Export it to Cytoscape-readable .graphml format
  write.graph(net$gR, 
              file=file.path(resultsDir, 
                             paste0("relNetwork", 
                                   as.character(netw.threshold),".graphml")), 
              format = "graphml")
  
  ## And obtain and save the specific pairs having corr value over the threshold
  net.out <- net$M
  net.out[abs(net.out)<netw.threshold] <- NA #drop cases below the threshold
  net.out <- as.data.frame(as.table(net.out)) #turn into a 3-column table
  net.out <- na.omit(net.out)
  colnames(net.out) <- c("X-feat", "Y-feat", "CorrValue")
  net.out <- net.out[order(-abs(net.out$CorrValue)),]    #sort by highest abs corr value
  write.csv2(net.out, file=file.path(resultsDir, 
                                paste("relNetworkValues", 
                                      as.character(netw.threshold),".csv", 
                                      sep="")),
             row.names = TRUE)
  return(net)
}



plot_cim <- function(rccaResult, resultsDir, x.lab = "X", y.lab = "Y"){
  ##
  ##
  require(mixOmics)
  
  ## Create and export CIM heatmap
  #cim(rccaResult, comp = 1:3, ylab = y.lab, xlab = x.lab, margins = c(5, 6))
  ## Save cim heatmap into a file
  png(file=file.path(resultsDir, "cimHeatMap.png"), 
      width=800, height=800, res=120)
    cim(rccaResult, comp = 1:3, ylab = x.lab, xlab = y.lab, ## x AND y LABELS INTENTIONALLY REVERSED
        margins = c(5, 6))
  dev.off()
}
