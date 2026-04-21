##########################################################################
##########################################################################
# Project:
# Script purpose: process the annotation file gtf of UKY_AmexF1_1 GCF_040938575.1_UKY_AmexF1_1
# the gtf and genome files were downloaded from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_040938575.1/
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Apr 21 10:36:44 2026
##########################################################################
##########################################################################
rm(list = ls())

library('rtracklayer')
library(GenomicRanges)
library('GenomicFeatures')


inputDir = '/groups/tanaka/People/current/jiwang/Genomes/axolotl_new/UKY_AmexF1_1/annotations/RefSeq/'

annotDir = '/groups/tanaka/People/current/jiwang/Genomes/axolotl_new/UKY_AmexF1_1/annotations/RefSeq/'

##########################################
# import the original gtf file and collect protein coding genes  
##########################################
gtf = paste0(annotDir, '/genomic.gtf')
annot = import(gtf, format = "gtf")
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")

kk = which(annot$gene_biotype == 'protein_coding')

genes = data.frame(annot[kk], stringsAsFactors = FALSE)

ggs = data.frame(chrom = genes$seqnames,
                 start = genes$start, 
                 end = genes$end, 
                 names = genes$gene, 
                 scores = 0, 
                 strand = genes$strand, 
                 stringsAsFactors = FALSE)

write.table(ggs, file = paste0(annotDir, 'UKY_AmexF1_1_proteincoding_genes.bed'), 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

rm(ggs)



##########################################
# make transcripts and gene mapping files 
##########################################
kk = which(annot$transcript_biotype == 'mRNA' & annot$type == 'transcript')

transcripts = data.frame(annot[kk], stringsAsFactors = FALSE)
transcripts = data.frame(chrom = transcripts$seqnames, 
                         start = transcripts$start, 
                         end = transcripts$end, 
                         strand = transcripts$strand, 
                         width = transcripts$width, 
                         source = transcripts$source,
                         type = transcripts$type, 
                         gene = transcripts$gene, 
                         gene_id = transcripts$gene_id, 
                         gene_biotype = transcripts$gene_biotype, 
                         transcript_biotype = transcripts$transcript_biotype, 
                         product = transcripts$product, 
                         model_evidence = transcripts$model_evidence
                         )

mm = match(genes$gene, transcripts$gene)

write.table(transcripts, file = paste0(annotDir, 'UKY_AmexF1_1_proteincoding_transcript.to.gene_mapping.txt'), 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)



ggs$names = paste0(ggs$names, '.0')


aa = read.delim(file = paste0(inputDir, 'AmexT_v47_FULL_transcripts.all.txt'), header = FALSE)
aa$V1 = gsub('>', '', aa$V1)

annot = readRDS(file = paste0(annotDir, 'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID.rds'))
aa$gene = annot$geneID[match(aa$V1, annot$transcritID)]

write.table(aa, file = paste0(inputDir, 'AmexT_v47_transcripts_genes_t2g.txt'), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')