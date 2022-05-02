############################################################
##### Created by Ferran Brianso
##### ferran.brianso_at_gmail.com
############################################################


############################################################
##### To run this targets script:
## 1- load 'targets' and 'tarchetypes' libs
library(targets)
library(tarchetypes)

## 2- load functions, params and general options
source("R/functions.R")
source("R/params.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", "dplyr", "readr", "tidyr"))

## 3- run the following comands IN CONSOLE
#tar_make()       # to run the script
#tar_glimpse()    # to see a simple view of the workflow
#tar_visnetwork() # to see more details of the workflow
############################################################


############################################################
##### DO NOT RUN ANY OF THE REMAINING LINES
############################################################
list(
  ## 01- Lectura de dades d'entrada
  tar_target(
    data_file1,
    p.inFile1,
    format = "file"
  ),
  tar_target(
    data_file2,
    p.inFile2,
    format = "file"
  ),
  tar_target(
    in_data1,
    read_csv(data_file1, col_types = cols())
  ),
  tar_target(
    in_data2,
    read_csv(data_file2, col_types = cols())
  ),
  
  
  ## 02- Pas de les dades a data frames amb estructura feats x samples
  tar_target(
    dframe1,
    set_dframe(in_data1)
  ),

  tar_target(
    dframe2,
    set_dframe(in_data2)
  ),

  # ## Per crear el heat map amb ggplot2 de les dades inicials
  # tar_target(
  #   in_hmplot, 
  #   #create_plot(data.matrix(dframe))
  #   create_plot(data.matrix(dframe[1:30, 1:16]))
  # ),


  ## 03- Creacio de la taula de categories anotades (amb valors 1/0)
  tar_target(
    categ_matrix1, 
    get_categ_matrix(onto=p.GO.onto, resultsDir=p.resultsDir, N=p.GO.miN, df=dframe1, 
                     afile=p.annotFile, outTag = p.outTag1) 
  ),

  tar_target(
    categ_matrix2, 
    get_categ_matrix(onto=p.GO.onto, resultsDir=p.resultsDir, N=p.GO.miN, df=dframe2, 
                     afile=p.annotFile, outTag = p.outTag2) 
  ),

  
  tar_target(
    categ_intersect1,
    categ_matrix1[, intersect(colnames(categ_matrix1), colnames(categ_matrix2))]
  ),

  tar_target(
    categ_intersect2,
    categ_matrix2[, intersect(colnames(categ_matrix1), colnames(categ_matrix2))]
  ),

  tar_target(
    categ_intersect,
    rbind(categ_intersect1, categ_intersect2)
  ),
    
  ## 04- Recompte de gens anotats per categoria (amb valor 1)
  tar_target(
    categ_sums1, 
    colSums(categ_matrix1) 
  ),
    
  tar_target(
    categ_sums2, 
    colSums(categ_matrix2) 
  ),
  
  ## 05- Unio del data frame amb les dades i el de les categories anotades  
  tar_target(
    annot_matrix1, 
    ## annot.matrix will be the transposition of the original matrix 
    ## with some GO categs (with 1/0 values) added as extra columns
    cbind(data.matrix(dframe1), categ_matrix1)
  ),
  
  tar_target(
    annot_matrix2, 
    ## annot.matrix will be the transposition of the original matrix 
    ## with some GO categs (with 1/0 values) added as extra columns
    cbind(data.matrix(dframe2), categ_matrix2)
  ),


  ## 06- Expansio de files afegint les mitjanes dels valors de les categories anotades
  tar_target(
    expd_matrix1,
    expand_annot_matrix(annot_matrix1, p.resultsDir,
                        1:ncol(dframe1), # range of cols containing samples
                        (ncol(dframe1)+1):(ncol(dframe1)+ncol(categ_matrix1)), # idem for categs
                        method="mean", outTag = p.outTag1) 
  ),
  
  tar_target(
    expd_matrix2,
    expand_annot_matrix(annot_matrix2, p.resultsDir,
                        1:ncol(dframe2), # range of cols containing samples
                        (ncol(dframe2)+1):(ncol(dframe2)+ncol(categ_matrix2)), # idem for categs
                        method="mean", outTag = p.outTag2) 
  ),

  tar_target(
    exc_exp_mat1,
    get_exclusive_expmat(annot_matrix1, annot_matrix2, p.resultsDir,
                        s.cols = 1:ncol(dframe1), # range of cols containing samples
                        c.cols1 = (ncol(dframe1)+1):(ncol(dframe1)+ncol(categ_matrix1)), # idem for categs
                        c.cols2 = (ncol(dframe2)+1):(ncol(dframe2)+ncol(categ_matrix2)), # idem for categs
                        outTag = p.outTag1) 
  ),
  
  tar_target(
    exc_exp_mat2,
    get_exclusive_expmat(annot_matrix2, annot_matrix1, p.resultsDir,
                           s.cols = 1:ncol(dframe2), # range of cols containing samples
                           c.cols1 = (ncol(dframe2)+1):(ncol(dframe2)+ncol(categ_matrix2)), # idem for categs
                           c.cols2 = (ncol(dframe1)+1):(ncol(dframe1)+ncol(categ_matrix1)), # idem for categs
                           outTag = p.outTag2) 
  ),
  
  tar_target(
    int_exp_mat,
    get_intersect_expmat(annot_matrix1, annot_matrix2, p.resultsDir,
                         s.cols = 1:ncol(dframe1), # range of cols containing samples
                         c.cols1 = (ncol(dframe1)+1):(ncol(dframe1)+ncol(categ_matrix1)), # idem for categs
                         c.cols2 = (ncol(dframe2)+1):(ncol(dframe2)+ncol(categ_matrix2)), # idem for categs
                         wt1 = 0.25, wt2 = 0.75) 
  ),
  
  
  
  # ## Per crear el heat map amb ggplot2 de les dades expandides
  # tar_target(
  #   expd_hmplot, 
  #   #create_plot(expd_matrix)
  #   create_plot(expd_matrix[c(1:24,(nrow(annot_matrix)+1):nrow(expd_matrix)), 1:16])
  # ),
  

  ## 07- Creacio de l'informe en HTML
  tar_render(report, "report.Rmd")
)

############################################################
############################################################