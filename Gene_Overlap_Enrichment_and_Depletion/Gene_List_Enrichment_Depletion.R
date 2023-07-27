rm(list = ls())
library(rio)
source('~/onlybiohpc/TUTORIALS/Gene_Overlap_Enrichment_and_Depletion/Tutorial_Functions.R')

# Set the working directory to wherever you have all the tutorial materials
setwd('~/onlybiohpc/TUTORIALS/')

# Load major cell type markers
ctmarks_adult = rio::import('DATASETS/Human_AdultBrain_CellType_Markers.xlsx')
ctmarks_fetal = rio::import('DATASETS/Human_FetalBrain_CellType_Markers.xlsx')

# Create list of cell type marker genes
ctmarks_adultL = split(ctmarks_adult, ctmarks_adult$CellType)
ctmarks_adultL = lapply(ctmarks_adultL, function(x){x$Gene})

ctmarks_fetalL = split(ctmarks_fetal, ctmarks_fetal$CellType)
ctmarks_fetalL = lapply(ctmarks_fetalL, function(x){x$Gene})

# Load gene list of interest (SFARI genes in this example)
sfari = rio::import('DATASETS/SFARI-Gene_genes_07-17-2023release_07-26-2023export.csv')

# Create list of sfari genes split by their score
sfariL = split(sfari, sfari$`gene-score`)
names(sfariL) = paste0('Score_', names(sfariL))
sfariL = lapply(sfariL, function(x){x$`gene-symbol`})

## ENRICHMENT WITH ADULT ##
# Set the background for enrichment test. We will use all genes tested for cell type marker analysis.
matForDEG = readRDS('DATASETS/ADULT_pseudobulk_perBroadCellType_normalized.RDS')
bcg = nrow(matForDEG)

geneOvEnrDep(gnL1 = sfariL, gnL2 = ctmarks_adultL, bcg, plot = T, hg = 10, wd = 14, fn = 'Adult_Cortical_CellTypeMarkers_SFARI_EnrichmentDepletion')


## ENRICHMENT WITH FETAL ##
# Set the background for enrichment test. We will use all genes tested for cell type marker analysis.
matForDEG = readRDS('DATASETS/FETAL_pseudobulk_perBroadCellType_normalized.RDS')
bcg = nrow(matForDEG)

geneOvEnrDep(gnL1 = sfariL, gnL2 = ctmarks_fetalL, bcg, plot = T, hg = 10, wd = 14, fn = 'Fetal_Cortical_CellTypeMarkers_SFARI_EnrichmentDepletion')












