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
#tar_glimpse()    # to see a simple workflow
#tar_visnetwork() # to see a more complete workflow
############################################################


############################################################
##### DO NOT RUN ANY OF THE REMAINING LINES
############################################################
list(
  ## 01- Lectura de dades d'entrada
  tar_target(
    data_file,
    p.inFile,
    format = "file"
  ),
  tar_target(
    in_data,
    read_csv(data_file, col_types = cols())
  ),
  
  
  ## 02- Pas de les dades a data frame amb estructura feats x samples
  tar_target(
    dframe,
    set_dframe(in_data)
  ),

  # ## Per crear el heat map amb ggplot2 de les dades inicials
  # tar_target(
  #   in_hmplot, 
  #   #create_plot(data.matrix(dframe))
  #   create_plot(data.matrix(dframe[1:30, 1:16]))
  # ),


  ## 03- Creacio de la taula de categories anotades (amb valors 1/0)
  tar_target(
    categ_matrix, 
    get_categ_matrix(onto=p.GO.onto, N=p.GO.miN, df=dframe, afile=p.annotFile) 
  ),
  
  ## 04- Recompte de gens anotats per categoria (amb valor 1)
  tar_target(
    categ_sums, 
    colSums(categ_matrix) 
  ),
    

  ## 05- Unio del data frame amb les dades i el de les categories anotades  
  tar_target(
    annot_matrix, 
    ## annot.matrix will be the transposition of the original matrix 
    ## with some GO categs (with 1/0 values) added as extra columns
    cbind(data.matrix(dframe), categ_matrix)
  ),
  
  ## 06- Expansio de files afegint les mitjanes dels valors de les categories anotades
  tar_target(
    expd_matrix,
    expand_annot_matrix(annot_matrix, 
                        1:ncol(dframe), # range of cols containing samples
                        (ncol(dframe)+1):(ncol(dframe)+ncol(categ_matrix)), # idem for categs
                        method="mean")
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