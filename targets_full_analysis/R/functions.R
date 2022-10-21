#' @title Functions to run the "Annotate and Expand Omics Data Matrix" script as a Targets' Pipeline
#' @description .......
#' @return ..........
#' @param data Data frame, preprocessed expression matrix as df.
#' @examples
#' library(ggplot2)
#' library(tidyverse)
#' library(stats)

library(topGO)
library(org.Hs.eg.db)
library(ggplot2)
library(stats)
library(mixOmics)
library(FactoMineR)
library(factoextra)


set_dframe <- function(my.data, resultsDir, outTag) {
  ## expected input: samples in rows (1st column is sample name!!)
  ## feature ids in the remaining columns
  ## RETURNS: a data frame with FEATURES in ROWS and SAMPLES in COLUMNS
  my.frame <- as.data.frame(my.data)
  rownames(my.frame) <- my.frame$sample
  my.frame <- my.frame[ ,-1]
  my.matrix <- t(as.matrix(my.frame))
  save(my.matrix, file = file.path(resultsDir, paste0("input.matrix", outTag, ".Rda")))
  return(t(my.frame))
}


get_categ_matrix <- function(onto="BP", resultsDir, N=10, df, afile=NA, outTag = "") {
  # if afile is NA, annotate the data considering onto and N
  if (is.na(afile)){
    # onto: should be one of "BP"(default), "MF" or "CC"
    # N: number to filter out those GO categs not having N+ elements from our gene list
    # df: should be an expression matrix with gene symbols in rows and samples in columns
    ## RETURNS: the given data frame with biological categories added, with cbind, as extra columns
    annot.list <- annFUN.org(onto, feasibleGenes = rownames(df), 
                             mapping = "org.Hs.eg.db", ID = "symbol")

    # filtro casos en que tinguem N+ elements en la categoria GO
    GOgtN <- annot.list[which(lengths(annot.list)>=N)]
    num.categs <- length(GOgtN)
  
    # from now on, samples will be in colnames and gene symbols in rownames
    num.feats <- length(rownames(df))
    num.samples <- length(colnames(df))
  
    # inicialitzo la matriu d'anotacions
    categ.matrix <- matrix(0, nrow = num.feats, ncol = num.categs)
    rownames(categ.matrix) <- rownames(df)
    colnames(categ.matrix) <- names(GOgtN)
  
    # assigno les anotacions trobades
    for (j in 1:num.categs){
      categ.matrix[which(rownames(categ.matrix) %in% GOgtN[[j]]), j] <- 1
    }
  # otherwise, load annotations from corresponding input  file
  }else{
    categ.matrix <- read.csv(afile, header = TRUE, sep = ",")
    categ.matrix <- data.matrix(categ.matrix[ ,-1]) # first column should be exactly the same feature names 
                                                    # as provided for the expression data 
  }
  
  # guardo i retorno la matriu de categories
  save(categ.matrix, file = file.path(resultsDir, paste0("categ.matrix", outTag, ".Rda")))
  return(categ.matrix)
}  


###################################################
## expand_annot_matrix call example
###################################################
#expd_matrix <- expand_annot_matrix(
#  annot.matrix, ## matrix including expr values and binary annotations to bio categs
#  1:num.samples, ## range of columns containing sample-gene values (1:N in the demo)
#  (num.samples+1):(num.samples+num.categs), ## range of columns with categs (the rest of the cols)
#  method = "mean"
#)
expand_annot_matrix <- function(x, resultsDir, s.cols, c.cols, method="mean", outTag = "") {
  # function that, given a matrix including expr values for genes (rows) and samples + binary 
  ## annotations to bio categs (in columns), returns the same matrix with additional 
  ## rows containing the "mean", or "sum", values of those genes included in the categories
  for (j in 1:length(c.cols)){
    if (method=="sum") {
      vals <- apply(x[,s.cols], 2, function(y){sum(y[which(x[,c.cols[j]]==1)])})
    }else{
      vals <- apply(x[,s.cols], 2, function(y){mean(y[which(x[,c.cols[j]]==1)])})
    }
    #print(vals)
    if (j==1){
      expanded.matrix <- rbind(x, c(vals,rep(NA, length(c.cols))))
      expanded.matrix[dim(expanded.matrix)[1],c.cols[j]] <- length(which(x[,c.cols[j]]==1))
    }else{
      expanded.matrix <- rbind(expanded.matrix, c(vals,rep(NA, length(c.cols))))
      expanded.matrix[dim(expanded.matrix)[1],c.cols[j]] <- length(which(x[,c.cols[j]]==1))
    }
    rownames(expanded.matrix)[dim(expanded.matrix)[1]] <- colnames(x)[length(s.cols)+j]
  }
  expanded.matrix <- expanded.matrix[ , s.cols]
  
  ## guardo i retorno la matriu expandida
  save(expanded.matrix, file = file.path(resultsDir, paste0("expanded.matrix", outTag, ".Rda")))
  return(expanded.matrix)
}


get_exclusive_expmat <- function(x, y, resultsDir, s.cols, c.cols1, c.cols2, outTag = "") {
  # function that, given a pair of matrices, including expr values for gene-like features (rows) and samples + binary 
  ## annotations to bio categs (in columns, s.cols + n.cols), returns a single matrix with rows containing the "mean",
  ## values of those features included in annotated categories coming exclusively from the first matrix (x), 
  ## not present in the second (y)
  for (j in 1:length(c.cols1)){
    vals <- apply(x[,s.cols], 2, function(f){mean(f[which(x[,c.cols1[j]]==1)])})
    #print(vals)
    if (j==1){
      expanded.matrix1 <- rbind(x, c(vals,rep(NA, length(c.cols1))))
      expanded.matrix1[dim(expanded.matrix1)[1],c.cols1[j]] <- length(which(x[,c.cols1[j]]==1))
    }else{
      expanded.matrix1 <- rbind(expanded.matrix1, c(vals,rep(NA, length(c.cols1))))
      expanded.matrix1[dim(expanded.matrix1)[1],c.cols1[j]] <- length(which(x[,c.cols1[j]]==1))
    }
    rownames(expanded.matrix1)[dim(expanded.matrix1)[1]] <- colnames(x)[length(s.cols)+j]
  }
  expanded.matrix1 <- expanded.matrix1[ , s.cols]

  for (j in 1:length(c.cols2)){
    vals <- apply(y[,s.cols], 2, function(f){mean(f[which(y[,c.cols2[j]]==1)])})
    #print(vals)
    if (j==1){
      expanded.matrix2 <- rbind(y, c(vals,rep(NA, length(c.cols2))))
      expanded.matrix2[dim(expanded.matrix2)[1],c.cols2[j]] <- length(which(y[,c.cols2[j]]==1))
    }else{
      expanded.matrix2 <- rbind(expanded.matrix2, c(vals,rep(NA, length(c.cols2))))
      expanded.matrix2[dim(expanded.matrix2)[1],c.cols2[j]] <- length(which(y[,c.cols2[j]]==1))
    }
    rownames(expanded.matrix2)[dim(expanded.matrix2)[1]] <- colnames(y)[length(s.cols)+j]
  }
  expanded.matrix2 <- expanded.matrix2[ , s.cols]
  
  intersect.matrix <- expanded.matrix1[intersect(rownames(expanded.matrix1)[(nrow(x)+1):(nrow(x)+length(c.cols1))], 
                                                 rownames(expanded.matrix2)[(nrow(y)+1):(nrow(y)+length(c.cols2))]), ]
  
  exclusive_expmat <- expanded.matrix1[-which(rownames(expanded.matrix1) %in% rownames(intersect.matrix)), ]
  exclusive_expmat <- exclusive_expmat[(nrow(x)+1):(nrow(exclusive_expmat)), ]

  ## guardo i retorno les matrius amb anotacions exclusives
  save(exclusive_expmat, file = file.path(resultsDir, paste0("exclusive.matrix", outTag, ".Rda")))

  return(exclusive_expmat)
}


get_intersect_expmat <- function(x, outTag1 = "1", y, outTag2 = "2", resultsDir, s.cols, c.cols1, c.cols2, wt1 = 0.5, wt2 = 0.5) {
  # function that, given a pair of matrices, including expr values for gene-like features (rows) and samples + binary 
  ## annotations to bio categs (in columns, s.cols + n.cols), returns a single matrix with rows containing the weighted
  ## values of those features included in common annotated categories (present in both x and y annot.matrices)
  for (j in 1:length(c.cols1)){
    vals <- apply(x[,s.cols], 2, function(f){mean(f[which(x[,c.cols1[j]]==1)])})
    #print(vals)
    if (j==1){
      expanded.matrix1 <- rbind(x, c(vals,rep(NA, length(c.cols1))))
      expanded.matrix1[dim(expanded.matrix1)[1],c.cols1[j]] <- length(which(x[,c.cols1[j]]==1))
    }else{
      expanded.matrix1 <- rbind(expanded.matrix1, c(vals,rep(NA, length(c.cols1))))
      expanded.matrix1[dim(expanded.matrix1)[1],c.cols1[j]] <- length(which(x[,c.cols1[j]]==1))
    }
    rownames(expanded.matrix1)[dim(expanded.matrix1)[1]] <- colnames(x)[length(s.cols)+j]
  }
  expanded.matrix1 <- expanded.matrix1[ , s.cols]
  
  for (j in 1:length(c.cols2)){
    vals <- apply(y[,s.cols], 2, function(f){mean(f[which(y[,c.cols2[j]]==1)])})
    #print(vals)
    if (j==1){
      expanded.matrix2 <- rbind(y, c(vals,rep(NA, length(c.cols2))))
      expanded.matrix2[dim(expanded.matrix2)[1],c.cols2[j]] <- length(which(y[,c.cols2[j]]==1))
    }else{
      expanded.matrix2 <- rbind(expanded.matrix2, c(vals,rep(NA, length(c.cols2))))
      expanded.matrix2[dim(expanded.matrix2)[1],c.cols2[j]] <- length(which(y[,c.cols2[j]]==1))
    }
    rownames(expanded.matrix2)[dim(expanded.matrix2)[1]] <- colnames(y)[length(s.cols)+j]
  }
  expanded.matrix2 <- expanded.matrix2[ , s.cols]
  
  intersect.matrix1 <- expanded.matrix1[intersect(rownames(expanded.matrix1)[(nrow(x)+1):(nrow(x)+length(c.cols1))], 
                                                 rownames(expanded.matrix2)[(nrow(y)+1):(nrow(y)+length(c.cols2))]), ]
  intersect.matrix2 <- expanded.matrix2[intersect(rownames(expanded.matrix2)[(nrow(y)+1):(nrow(y)+length(c.cols2))], 
                                                  rownames(expanded.matrix1)[(nrow(x)+1):(nrow(x)+length(c.cols1))]), ]
  
  intersect.matrix <- intersect.matrix1*wt1 + intersect.matrix2*wt2
    
  ## guardo i retorno les matrius expandides
  save(intersect.matrix, file = file.path(resultsDir, paste0("intersect.matrix.Rda")))
  save(intersect.matrix1, file = file.path(resultsDir, paste0("intersect.matrix", outTag1, ".Rda")))
  save(intersect.matrix2, file = file.path(resultsDir, paste0("intersect.matrix", outTag2, ".Rda")))
  
  return(intersect.matrix)
}


get_groupInfo <- function() {
  ## expected input: none
  ## RETURNS: a list with the cancer subtypes (as class factor list)
  ## Note: will be able to load the factor in future versions
  require(mixOmics)
  data("breast.TCGA")
  brca.subtype <- breast.TCGA$data.train$subtype
  return(brca.subtype)
}

get_basicMfaResults <- function(basic.data, groups, group.names, group.types) {
  require(FactoMineR)
  mfa.res <- MFA(base=basic.data, #jointBasicDataSc,
                 group=groups, #blocksBasic,
                 name.group = group.names, #p.BlockNamesMFA[1:2],
                 type=group.types, #p.BlockTypesMFA[1:2],
                 graph=FALSE)
  return(mfa.res)
}

get_expMfaResults <- function(expd.data, groups, group.names, group.types, supplGroups) {
  require(FactoMineR)
  mfa.res <- MFA(base=expd.data, #jointBasicDataSc,
                 group=groups, #blocksBasic,
                 name.group = group.names, #p.BlockNamesMFA[1:2],
                 type=group.types, #p.BlockTypesMFA[1:2],
                 num.group.sup=supplGroups,
                 graph=FALSE)
  return(mfa.res)
}
