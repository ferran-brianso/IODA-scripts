############################################################
##### Created by Ferran Briansó
##### ferran.brianso_at_gmail.com
############################################################
#####
##### files for Annotation and Expansion of Omics Data
#####
############################################################

############################################################
##### To run this targets script:
## 0- allocate all files and folders in you working dir
##    and open _targets.R file in Rstudio

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
