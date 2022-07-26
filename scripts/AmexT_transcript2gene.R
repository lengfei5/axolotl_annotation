##########################################################################
##########################################################################
# Project:
# Script purpose: mapping transcripts to genes
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jan 11 16:55:06 2022
##########################################################################
##########################################################################
inputDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/Transcriptomics/'
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'

aa = read.delim(file = paste0(inputDir, 'AmexT_v47_FULL_transcripts.all.txt'), header = FALSE)
aa$V1 = gsub('>', '', aa$V1)

annot = readRDS(file = paste0(annotDir, 'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID.rds'))
aa$gene = annot$geneID[match(aa$V1, annot$transcritID)]

write.table(aa, file = paste0(inputDir, 'AmexT_v47_transcripts_genes_t2g.txt'), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')





