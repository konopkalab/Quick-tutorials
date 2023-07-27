geneOvEnrDep = function(gnL1, gnL2, bcg, plot = T, fn = 'Gene_Enrichment_Plot', wd = 20, hg = 10, xlabel = '', ylabel = ''){

	require(ggpubr)
	require(dplyr)

	# Fisher's exact test enrichment
	resL = list()
	for(i in 1:length(gnL2)){

		resL[[i]] = lapply(1:length(gnL1), function(x){
			A = gnL1[[x]]
			B = gnL2[[i]]
			q1 = bcg - length(union(A,B))
			q2 = length(setdiff(A,B))
			q3 = length(setdiff(B,A))
			q4 = length(intersect(A,B))
			mat = matrix(c(q1,q2,q3,q4), nrow=2)
			res = fisher.test(mat)
			or = res$estimate
			pval = res$p.value
			tmp = data.frame(var1 = names(gnL1)[x], var2 = names(gnL2)[i], OddsRatio = or, pval = pval)
			tmp
		}) %>% do.call(rbind, .)

		resL[[i]]$FDR = p.adjust(resL[[i]]$pval, method = 'fdr')
	}

	toplot = do.call(rbind, resL)
	rownames(toplot) = NULL

	# FDR for scaling
	toplot$log10FDR = -log10(toplot$FDR)

	# FDR for labeling
	toplot$sign_label = formatC(toplot$FDR, format = "e", digits = 1)

	# Asterisk for labeling
	toplot$sign_label_2 = ifelse(toplot$FDR < 0.05, '*', '')

	# Plot or return the data frame
	if(plot == T){

	pdf(paste0(fn, '.pdf'), width = wd, height = hg)
		print( ggscatter(toplot, x = 'var1', y = 'var2', color = 'OddsRatio', size = 'log10FDR') +
		geom_label(data = toplot, aes(label = sign_label), color="black",
			label.size = NA, fill = alpha(c("white"),0), fontface = 'bold', size = 6) + 
		labs(x=xlabel, y=ylabel) +
		scale_size_continuous(range = c(5,25)) +
		scale_color_gradient2(midpoint = 1, low = 'blue', high = 'red') +
		theme_classic() +
		theme(text = element_text(size=20, face = 'bold')) +
		geom_text(aes(label = sign_label_2), vjust = 1.2, colour = "darkgreen", fontface = 'bold', size = 20 ) +
		rotate_x_text(45) )
	dev.off()

	}

	if(plot == F){
		return(toplot)
	}
}

