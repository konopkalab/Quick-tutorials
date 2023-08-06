
plotVenn = function(gnList, fn, margins = 0.1){

	require(VennDiagram)

	# Check if the list is greater than 3 gene sets.
	if(length(gnList) > 3){
		stop("Your input has >3 gene lists. Use an upset plot.")
	} else if(length(gnList) == 2){

		# Extract gene sets
		set1 = gnList[[1]]
		set2 = gnList[[2]]

		# Extract names of the gene sets
		nm1 = names(gnList)[1]
		nm2 = names(gnList)[2]

		# Plot double venn diagram
		pdf(paste0(fn, '.pdf'), width = 8)
		print( draw.pairwise.venn(area1 = length(set1),
					area2 = length(set2),
					cross.area = length(intersect(set1, set2)),
					fill = c('orange', 'darkgreen'), col = c('orange', 'darkgreen'),
					category = c(nm1, nm2), fontface = 'bold', cex = 2,
					cat.cex = 2, cat.fontface = 'bold', margin = rep(margins, 4)) )
		dev.off()

	} else if(length(gnList) == 3){

		# Extract gene sets
		set1 = gnList[[1]]
		set2 = gnList[[2]]
		set3 = gnList[[3]]

		# Extract names of the gene sets
		nm1 = names(gnList)[1]
		nm2 = names(gnList)[2]
		nm3 = names(gnList)[3]

		# Plot triple venn diagram
		pdf(paste0(fn, '.pdf'), width = 8)
		print( draw.triple.venn(area1 = length(set1),
					area2 = length(set2),
					area3 = length(set3),
					n12 = length(intersect(set1, set2)),
					n13 = length(intersect(set1, set3)),
					n23 = length(intersect(set2, set3)),
					n123 = length(Reduce(intersect, list(set1, set2, set3))),
					fill = c('orange', 'darkgreen', 'lightblue'),
					col = c('orange', 'darkgreen', 'lightblue'),
					category = c(nm1, nm2, nm3), fontface = 'bold', cex = 2,
					cat.cex = 2, cat.fontface = 'bold', margin = rep(margins, 4)) )
		dev.off()
	}

}



