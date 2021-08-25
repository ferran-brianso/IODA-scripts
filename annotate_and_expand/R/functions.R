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


set_dframe <- function(my.data) {
  ## expected input: samples in rows (1st column is sample name!!) and features in columns
  ## feature ids in the remaining columns
  ## RETURNS: a data frame with FEATURES in ROWS and SAMPLES in COLUMNS
  my.frame <- as.data.frame(my.data)
  rownames(my.frame) <- my.frame$sample
  my.frame <- my.frame[ ,-1]
  return(t(my.frame))
}


# create_plot <- function(my.matrix) {
#   my.matrix %>% 
#     as.data.frame() %>% 
#     rownames_to_column("feature") %>% 
#     pivot_longer(-c(feature), names_to = "sample", values_to = "values") %>% 
#     ggplot(aes(x=sample, y=feature, fill=values)) + 
#       geom_raster() +
#       scale_fill_viridis_c()
# }


get_categ_matrix <- function(onto="BP", N=10, df, afile=NA) {
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
    categ.matrix <- data.matrix(categ.matrix[ ,-1])
  }
  
  # retorno la matriu de categories
  return(categ.matrix)
}  


expand_annot_matrix <- function(x, s.cols, c.cols, method="mean"){
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
      em <- rbind(x, c(vals,rep(NA, length(c.cols))))
      em[dim(em)[1],c.cols[j]] <- length(which(x[,c.cols[j]]==1))
    }else{
      em <- rbind(em, c(vals,rep(NA, length(c.cols))))
      em[dim(em)[1],c.cols[j]] <- length(which(x[,c.cols[j]]==1))
    }
    rownames(em)[dim(em)[1]] <- colnames(x)[length(s.cols)+j]
  }
  em <- em[ , s.cols]
  return(em)
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
