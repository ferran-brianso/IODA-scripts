############################################################
##### Created by Ferran Brianso
##### ferran.brianso_at_gmail.com
############################################################


############################################################
##### To run this targets script:

## 0- Set working directory to Source File Location
# setwd("~/2021-IODA/IODA-scripts/targets_full_analysis")

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
  ## 1.01- Lectura de dades d'entrada
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
  
  
  ## 1.02- Pas de les dades a data frames amb estructura feats x samples
  tar_target(
    dframe1,
    set_dframe(in_data1, resultsDir=p.resultsDir, outTag=p.outTag1)
  ),

  tar_target(
    dframe2,
    set_dframe(in_data2, resultsDir=p.resultsDir, outTag=p.outTag2)
  ),

  # ## Per crear el heat map amb ggplot2 de les dades inicials
  # tar_target(
  #   in_hmplot, 
  #   #create_plot(data.matrix(dframe))
  #   create_plot(data.matrix(dframe[1:30, 1:16]))
  # ),


  ## 1.03- Creacio de la taula de categories anotades (amb valors 1/0)
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
    
  ## 1.04- Recompte de gens anotats per categoria (amb valor 1)
  tar_target(
    categ_sums1, 
    colSums(categ_matrix1) 
  ),
    
  tar_target(
    categ_sums2, 
    colSums(categ_matrix2) 
  ),
  
  ## 1.05- Unio del data frame amb les dades i el de les categories anotades  
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


  ## 1.06- Expansio de files afegint les mitjanes dels valors de les categories anotades
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
    get_intersect_expmat(annot_matrix1, outTag1 = p.outTag1, 
                         annot_matrix2, outTag2 = p.outTag2, 
                         p.resultsDir,
                         s.cols = 1:ncol(dframe1), # range of cols containing samples
                         c.cols1 = (ncol(dframe1)+1):(ncol(dframe1)+ncol(categ_matrix1)), # idem for categs
                         c.cols2 = (ncol(dframe2)+1):(ncol(dframe2)+ncol(categ_matrix2)), # idem for categs
                         wt1 = p.weight1, wt2 = p.weight2) 
  ),
  
  
  # ## Per crear el heat map amb ggplot2 de les dades expandides
  # tar_target(
  #   expd_hmplot, 
  #   #create_plot(expd_matrix)
  #   create_plot(expd_matrix[c(1:24,(nrow(annot_matrix)+1):nrow(expd_matrix)), 1:16])
  # ),
  
  ## 1.07- Creacio de l'informe en HTML part1
  #tar_render(report1, "report_part1.Rmd"),

  
  
  ############################################
  ############################################
  ### PART 2 Processat del MFA 
  ### FALTA CONNECTAR AQUESTA PART AMB L'ANTERIOR SENSE CARREGAR FITXERS EXTERNS
  
  ## 2.01- to apply MFA we need a unique matrix with obtained by merging all datasets 
  ## through its common dimension, which is set to be "rows", transposing the original matrices
  tar_target(
    jointBasicData,
    cbind(t(as.matrix(dframe1)), t(as.matrix(dframe2)))
  ),
  
  tar_target(
    jointBasicDataSc,
    cbind(scale(t(as.matrix(dframe1))), scale(t(as.matrix(dframe2))))
  ),

  tar_target(
    jointExpData,
    cbind(t(as.matrix(dframe1)), t(as.matrix(dframe2)),
          t(exc_exp_mat1), t(exc_exp_mat2),
          t(int_exp_mat))
  ),
  
  tar_target(
    jointExpDataSc,
    cbind(scale(t(as.matrix(dframe1))), scale(t(as.matrix(dframe2))),
          scale(t(exc_exp_mat1)), scale(t(exc_exp_mat2)),
          scale(t(int_exp_mat)))
  ),

  ##numSamples <- ncol(annotGenes)  <= exc_exp_mat1
  ##numGenes <- nrow(matGenes)  <= dframe1
  ##numProts <- nrow(matProts)  <= dframe2
  ##numGOsInGenes <- nrow(annotGenes) <= exc_exp_mat1
  ##numGOsInProts <- nrow(annotProts) <= exc_exp_mat2
  ##numGOsInCommon <- nrow(annotCommon) <= int_exp_mat
  
  ## MFA data blocks are defined implicitly by the (number of) columns 
  ## of the data matrix they are made of

  tar_target(
    blocksBasic,
    c(nrow(dframe1), nrow(dframe2))
  ),
  
  tar_target(
    blocksExp,
    c(nrow(dframe1), nrow(dframe2), 
      nrow(exc_exp_mat1), nrow(exc_exp_mat2), 
      nrow(int_exp_mat))
  ),
  
  
  ## 2.02- Carrega la info de groups de mostres (cancer subtype en BRCA data)
  tar_target(
    load_grinfo,
    get_groupInfo()
  ),
    
  ## 2.XX- Creacio de l'informe en HTML part2
  tar_render(report2, "report_part2.Rmd")
)





############################################################
############################################################