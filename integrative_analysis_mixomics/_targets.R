############################################################
##### Created by Ferran Brianso
##### ferran.brianso_at_gmail.com
############################################################
#####
##### targets_.R file for Integrative Analysis with Mixomics
#####
############################################################

############################################################
##### To run this targets script:
## 1- load 'targets' and 'tarchetypes' libs
library(targets)
library(tarchetypes)

## 2- set the working dir to the source file directory

## 3- load functions, params and general options
source("R/functions.R")
source("R/params.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", "biglm", "dplyr", "readr", "tidyr"))

## 4- run the following comands IN CONSOLE
#tar_make()       # to run the script
#tar_glimpse()    # to see a simple view of the workflow
#tar_visnetwork() # to see the same with functions and params
############################################################


############################################################
##### DO NOT RUN ANY OF THE LINES BELOW
############################################################
list(
  ## 01- Carrega dades d'entrada
  tar_target(
    gene_data,
    load_and_check(p.inFile, type = "G")
  ),
  tar_target(
    prot_data,
    load_and_check(p.inFile, type = "P")
  ),
  
  ## 02- Prepara les dades
  tar_target(
    gene_sc,
    scale(gene_data)
  ),
  tar_target(
    prot_sc,
    scale(prot_data)
  ),
  tar_target(
    X,
    t(gene_sc)
  ),
  tar_target(
    Y,
    t(prot_sc)
  ),
  
  ## 03- Executa l'analisi integrativa
  tar_target(
    rCCA,
    perform_rCCA(X, Y, p.resultsDir, p.scoresFile, p.rccFile)
  ),

  
  ## 00- Creacio de l'informe en HTML
  tar_render(report, "report.Rmd")
)

############################################################
############################################################