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

## HAUREM DE DECIDIR EN QUINS ALTRES FORMATS ES PODEN LLEGIR LES DADES, 
## JA QUE PER ARA SOLS LES LLEGIM D'UN Rda

############################################################


############################################################
#### Params regarding mixomics-derived plots
## related targets: 
##       X
##       Y

p.export = TRUE  # TRUE/FALSE
p.resultsDir = "results/mixomics"
p.scoresFile = "cv.score.Rda"
p.rccFile = "rccResult.Rda"

p.circ.cutOffs = c(0.55, 0.65, 0.75)
p.netw.threshold = p.circ.cutOffs[1] #0.5

p.x.tag = "g."
p.y.tag = "p."

p.x.lab = "Genes"
p.y.lab = "Prots"

############################################################