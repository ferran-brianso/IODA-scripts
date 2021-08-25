############################################################
############################################################
##### Work in progress
############################################################
##### Created by Ferran Brianso
##### ferran.brianso_at_gmail.com
############################################################
############################################################


############################################################
############################################################
## Params to run the "Annotate and Expand Omics Data Matrix" 
## script as a Targets' Pipeline
############################################################
############################################################


############################################################
#### Params for loading the main raw data
## related targets: 
##       data_file

p.inFile="data/mrna.csv"    # csv file with raw data from a single omics source


############################################################
## Params for the biological annotation
## related targets: 
##       categ_matrix

p.GO.onto="BP" # onto: should be one of "BP"(default), "MF" or "CC"
p.GO.miN=8      # N: number to filter out those GO categs not having N+ elements from our gene list

p.annotFile="data/annotations.csv" # should be NA if not available!!!

if(!is.na(p.annotFile)){
  p.GO.onto <<- paste("None.",  "Annotations loaded from file:",  p.annotFile)
  p.GO.miN <<- "NA"
}


############################################################