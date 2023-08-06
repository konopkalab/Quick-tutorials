rm(list = ls())
library(ComplexUpset)
library(dplyr)
library(ggpubr)
source('~/onlybiohpc/TUTORIALS/Visualizations_of_gene_overlaps/Tutorial_Functions.R')

####
## PREPARE LIST OF GENE SETS
####

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

####
## OVERLAP BETWEEN TWO GENE SETS -- VENN DIAGRAM
####

# A list of genesets. Elements of the list will be overlapped. Let's use GABAergic markers from both adult and fetal
markersL = list(ctmarks_adultL[["Inhibitory"]], ctmarks_fetalL[["GABAergic"]])
names(markersL) = c('Inh_Adult', 'Inh_Fetal')

# Set file name
fn = 'Adult_Fetal_GABAergic_Overlap'

# This function utilizes VennDiagram package to plot pairwise or triple venn diagrams
plotVenn(gnL = markersL, fn = fn, margins = 0.1)


####
## OVERLAP BETWEEN MULTIPLE GENE SETS -- UPSET PLOTS
####


# Create a boolean matrix from the gene list
allgns = Reduce(union, ctmarks_adultL)
logmat = sapply(names(ctmarks_adultL), function(x){allgns %in% ctmarks_adultL[[x]]}) %>% as.data.frame
rownames(logmat) = allgns

# Cell types to compare
ctypes = colnames(logmat)

# Display all non-zero intersections
pdf('Adult_CellType_Marker_Overlaps.pdf', width = 10, height = 8)
ComplexUpset::upset(logmat, ctypes, name='Cell_Types', width_ratio=0.1, sort_sets = F,
set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=90))),
themes=upset_modify_themes(
        list(
            'intersections_matrix'=theme(text=element_text(size=30)),
            'Intersection size'=theme(text=element_text(size=30)),
            'overall_sizes'=theme(text=element_text(size=30)),
            'default'=theme(text=element_text(size=30))
        )))
dev.off()


# Let's do the same with fetal markers
allgns = Reduce(union, ctmarks_fetalL)
logmat = sapply(names(ctmarks_fetalL), function(x){allgns %in% ctmarks_fetalL[[x]]}) %>% as.data.frame
rownames(logmat) = allgns
ctypes = colnames(logmat)

pdf('Fetal_CellType_Marker_Overlaps.pdf', width = 10, height = 8)
ComplexUpset::upset(logmat, ctypes, name='Cell_Types', width_ratio=0.1, sort_sets = F,
set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=90))),
themes=upset_modify_themes(
        list(
            'intersections_matrix'=theme(text=element_text(size=30)),
            'Intersection size'=theme(text=element_text(size=30)),
            'overall_sizes'=theme(text=element_text(size=30)),
            'default'=theme(text=element_text(size=30))
        )))
dev.off()



