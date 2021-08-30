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
## Params to run the "IODA (Mixomics part)" in targets' mode
############################################################
############################################################


############################################################
#### Params for loading the main raw data
## related targets: 
##       in_data

p.inFile = "data/raw.data.Rda"  # Rda file with raw data from genes and prots (both with gene symbols)

############################################################


############################################################
#### Params regarding mixomics-derived plots
## related targets: 
##       X
##       Y

p.export = TRUE  # TRUE/FALSE
p.resultsDir = "results/mixomics"
p.scoresFile = "cv.score.Rdata"
p.rccFile = "rccResult.Rdata"

p.cutOffs = c(0.5, 0.6, 0.7)

############################################################