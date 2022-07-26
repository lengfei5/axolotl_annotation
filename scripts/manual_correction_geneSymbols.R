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



##########################################
# collect all genes  
##########################################
aa = rtracklayer::import(paste0(annotDir, 'AmexT_v47_Hox.patch.gtf'))
aa = aa[which(aa$type == 'gene')]
aa = as.data.frame(aa)
aa = aa[, c(1:5, 10:11)]

annot = readRDS(paste0(annotDir,
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse_manual_v1.rds'))
annot = annot[which(!is.na(annot$gene.symbol.toUse)), ]

mm = match(aa$gene_id, annot$geneID)

aa$gene_name[which(is.na(mm))] = NA
aa$gene_name[which(!is.na(mm))] = annot$gene.symbol.toUse[mm[which(!is.na(mm))]]

#aa$gene_name = gsub("\\s*\\([^\\)]+\\)", "", aa$gene_name[2])

write.table(aa, file = paste0(resDir, '/amexT_v47_hoxPatch.txt'), sep = '\t', col.names = TRUE, row.names = FALSE, 
            quote = FALSE)
