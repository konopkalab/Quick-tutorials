GOenrich = function(gns, uni, pCut = 0.05, qCut = 0.2, species = 'human'){

	require(clusterProfiler)

	if(species != 'human' & species != 'mouse'){
		stop('Please select either mouse or human')
	}
	if(species == 'mouse'){
		require(org.Mm.eg.db)
		orgdb = org.Mm.eg.db
	}
	if(species == 'human'){
		require(org.Hs.eg.db)
		orgdb = org.Hs.eg.db
	}
	
	# Change from symbol to entrezid
	gene_df = bitr(gns, fromType = "SYMBOL",
		toType = c("ENTREZID"),
		OrgDb = orgdb)

	uni_df = bitr(uni, fromType = "SYMBOL",
		toType = c("ENTREZID"),
		OrgDb = orgdb)
	
	# Run all gene ontology enrichments
	goterms = c('MF', 'BP', 'CC')
	egos = list()
	for(i in 1:length(goterms)){
		ego = enrichGO(gene          = gene_df$ENTREZID,
				universe      = uni_df$ENTREZID,
				OrgDb         = orgdb,
				ont           = goterms[i],
				pAdjustMethod = "BH",
				pvalueCutoff  = pCut,
				qvalueCutoff  = qCut)
		ego = setReadable(ego, OrgDb = orgdb)
		ego = as.data.frame(ego)
		if(nrow(ego) == 0){next}

		ego$GO = goterms[i]
		egos[[i]] = ego
	}
	
	# Combine
	egodf = do.call(rbind, egos)

	# Calculate overlap amounts and odds ratios
	egodf$ObsOv = sapply(1:nrow(egodf), function(x){gsub('/.*', '', egodf[x, 'GeneRatio']) %>% as.numeric()})
	egodf$ObsAll = sapply(1:nrow(egodf), function(x){gsub('.*/', '', egodf[x, 'GeneRatio']) %>% as.numeric()})
	egodf$BcgOv = sapply(1:nrow(egodf), function(x){gsub('/.*', '', egodf[x, 'BgRatio']) %>% as.numeric()})
	egodf$BcgAll = sapply(1:nrow(egodf), function(x){gsub('.*/', '', egodf[x, 'BgRatio']) %>% as.numeric()})
	egodf$OddsRatio = (egodf$ObsOv / egodf$ObsAll) / (egodf$BcgOv / egodf$BcgAll)

	return(egodf)
}


wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}


