##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Mar  9 09:28:52 2022
##########################################################################
##########################################################################

##########################################
# NPPA discovery 
##########################################
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
annot = readRDS(paste0(annotDir,
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
annot$manual = NA

jj = which(annot$geneID == 'AMEX60DD051098')
annot$manual[jj] = 'NPPA'
annot$gene.symbol.toUse[jj] = annot$manual[jj]

saveRDS(annot, file = paste0(annotDir,
                             'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse_manual_v1.rds'))
