geneOvEnr = function(gnL1, gnL2, bcg, plot = T, fn = 'Gene_Enrichment_Plot', wd = 20, hg = 10){

	require(GeneOverlap)
	require(reshape2)
	require(ggpubr)

	# Run enrichment
	resgom = newGOM(gnL1, gnL2, genome.size = bcg)
	pvalmat = getMatrix(resgom, name="pval")
	oddsmat = getMatrix(resgom, name="odds.ratio")
	fdrmat = apply(pvalmat, 2, function(x){p.adjust(x, method = 'fdr')})

	# Long format
	oddsrM = melt(oddsmat)
	pvalsM = melt(fdrmat)

	# Combine metrics
	toplot = cbind(oddsrM, FDR = pvalsM$value)

	# FDR for scaling
	toplot$log10FDR = -log10(toplot$FDR)

	# FDR for labeling
	toplot$sign_label = formatC(toplot$FDR, format = "e", digits = 1)

	# Asterisk for labeling
	toplot$sign_label_2 = ifelse(toplot$FDR < 0.05 &
			toplot$value > median(as.numeric(oddsmat)), '*', '')

	if(plot == T){

	pdf(paste0(fn, '.pdf'), width = wd, height = hg)
	print( ggscatter(toplot, x = 'Var1', y = 'Var2', color = 'value', size = 'log10FDR') +
	geom_label(data = toplot, aes(label = sign_label), color="black",
		label.size = NA, fill = alpha(c("white"),0), fontface = 'bold', size = 6) + 
	labs(x="", y="") +
	scale_size_continuous(range = c(2,25)) +
	scale_color_gradient2(midpoint = 1, low = 'blue', high = 'red') +
	theme_classic() +
	theme(text = element_text(size=20, face = 'bold')) +
	geom_text(aes(label = sign_label_2), vjust = 1.2, colour = "darkgreen", fontface = 'bold', size = 15 ) +
	rotate_x_text(45) )
	dev.off()

	}

	if(plot == F){
		return(toplot)
	}
}

