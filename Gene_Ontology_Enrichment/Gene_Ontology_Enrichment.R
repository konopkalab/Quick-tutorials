rm(list = ls())
library(rio)
library(ggpubr)

# Set the working directory to wherever you have all the tutorial materials
setwd('~/onlybiohpc/TUTORIALS/')
source('Gene_Ontology_Enrichment/Tutorial_Functions.R')

# Load major cell type markers
ctmarks = rio::import('DATASETS/Human_AdultBrain_CellType_Markers.xlsx')

# Set the background for enrichment test. We will use all genes tested for cell type marker analysis.
matForDEG = readRDS('DATASETS/ADULT_pseudobulk_perBroadCellType_normalized.RDS')
bcg = rownames(matForDEG)

# Loop through each cell type marker and perform GO enrichment
ctypes = unique(ctmarks$CellType)
resL = list()
for(i in 1:length(ctypes)){

	# Extract marker genes per cell type
	gnstotest = ctmarks[ctmarks$CellType == ctypes[i], 'Gene']

	# Run GO enrichment. Do not apply any statistical cutoff yet.
	res = GOenrich(gnstotest, bcg, pCut = 1, qCut = 1, species = 'human')
	
	# Add cell type info
	res$CellType = ctypes[i]

	resL[[i]] = res
	print(paste0('GO performed in: ', ctypes[i]))
}

# Combine all cell type results
resdf = do.call(rbind, resL)

# Keep significant GO enrichments (FDR < 0.05 and fold change > 1.3)
resdf_sign = resdf[resdf$OddsRatio > 1.3 & resdf$qvalue < 0.05,]

# Save the results
rio::export(resdf, 'GO_Enrichment_CellTypeMarkers_All.xlsx')
rio::export(resdf_sign, 'GO_Enrichment_CellTypeMarkers_Significant.xlsx')


## Plot some of the GO enrichments. Include statistics from the other cell types for comparison

# Selected GO terms
terms = c('myelination', 'extracellular matrix structural constituent', 'regulation of trans-synaptic signaling', 'immune receptor activity')

# Extract selected GO terms
plotenr = resdf[resdf$Description %in% terms,]

# For plotting
plotenr$log10FDR = -log10(plotenr$qvalue)
plotenr$is_sign = ifelse(plotenr$qvalue < 0.05, '*', '')
plotenr$gonameY = sapply(plotenr$Description, function(x){wrapper(x, width = 20)})

pdf('CellTypeMarkers_Gene_Ontology', width = 8, height = 6)
ggplot(plotenr, aes(x = CellType, y = gonameY)) +
  geom_point(aes(size = log10FDR, fill = OddsRatio), shape = 21, colour = "black") +
  scale_fill_gradient2(midpoint = 1, low = 'blue', high = 'red') +
  theme_classic() + xlab('') + ylab('') +
  theme(text = element_text(size=20, face = 'bold'), legend.pos = 'right') +
  rotate_x_text(45) +
  geom_text(aes(label = is_sign), vjust = 0.8, colour = "blue", size = 12 ) +
  scale_size_continuous(range = c(7,20))
dev.off()


















