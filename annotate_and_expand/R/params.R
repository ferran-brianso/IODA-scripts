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

p.inFile1="data/mrna.csv"    # csv file with raw data from a first omics source (mandatory)
p.inFile2="data/prots.csv"    # csv file with raw data from a second omics source (optional, or NA)


############################################################
#### Params for exporting results
## related targets: 
##

p.resultsDir = "results"  ## results folder
p.outTag1 = ".genes"       ## tag to be added to the outputs created
p.outTag2 = ".prots"       ## tag to be added to the outputs created


############################################################
## Params for the biological annotation
## related targets: 
##       categ_matrix

p.GO.onto="BP" # onto: should be one of "BP"(default), "MF" or "CC"
p.GO.miN=8      # N: number to filter out those GO categs not having N+ elements from our gene list

p.annotFile=NA #"data/annotations.csv" # should be NA if not available!!!

if(!is.na(p.annotFile)){
  p.GO.onto <<- paste("None.",  "Annotations loaded from file:",  p.annotFile)
  p.GO.miN <<- "NA"
}


############################################################