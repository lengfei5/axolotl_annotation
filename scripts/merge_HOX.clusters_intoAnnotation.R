##########################################################################
##########################################################################
# Project:
# Script purpose: add HOX A, B, C, D into annotation
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Apr  6 11:13:11 2022
##########################################################################
##########################################################################
rm(list = ls())

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
resDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/axolotl_annotation/results'
########################################################
########################################################
# Section : first double check the origial full annotation is the same as the one from UCSC
# 
########################################################
########################################################
library('rtracklayer')
library(GenomicRanges)
library('GenomicFeatures')

annot = import(paste0(annotDir, 'AmexT_v47.release.gtf'))
ucsc = import(paste0(annotDir, 'UCSC_downloaded_14032022/ax6_v47_scafford.gtf'))

seqnames(annot)
aa = annot[grep('^chr', seqnames(annot))]

transcript1 = sapply(aa$transcript_id, function(x) {x = unlist(strsplit(as.character(x), '[|]')); x[length(x)]})
transcript2 = sapply(ucsc$transcript_id, function(x) {x = unlist(strsplit(as.character(x), '[|]')); x[length(x)]})
transcript1 = unique(transcript1)
transcript2 = unique(transcript2)
transcript1 = transcript1[!is.na(transcript1)]
transcript2 = transcript2[!is.na(transcript2)]

cat(length(transcript1), ' transcript IDs from released version \n')
cat(length(transcript2), ' transcript IDs from UCSC version \n')

## check the overlap 
mm = match(transcript1, transcript2)
length(intersect(transcript1, transcript2))
length(which(is.na(mm))) ## all transcript ID from release version is covered in the UCSC version

mm = match(transcript2, transcript1)
length(which(is.na(mm)))

kk = which(is.na(mm)) ## 9 transcripts found in UCSC but not in released version, all from patch.gtf in UCSC
transcript2[kk]  # "Col4a5", "Col4a5_dup1", "Hoxd1", "Hoxd4", "Hoxd8", "Hoxd9","Hoxd10", "Hoxd11", "Hoxd13"


########################################################
########################################################
# Section : add HOX A, B, C D clusters into both gtf files
# 
########################################################
########################################################
hox2 = rtracklayer::import('../toSave/HoxD_cluster.gtf')

hox1 = read.delim(file = paste0('../toSave/HoxABC_chr_cleaned.bed'), sep = '\t', header = FALSE, stringsAsFactors = FALSE)
hox1$V4 = gsub('HOx', 'Hox', hox1$V4)

hox1$V6[grep('HoxB|HoxC', hox1$V4)] = '+'
genes = hox1$V4

hox1 = makeGRangesFromDataFrame(hox1, seqnames.field=c("V1"),
                                  start.field="V2", end.field="V3", strand.field="V6", keep.extra.columns = FALSE)

hox1$source = 'ambMex60DD_annotation_patch.hoxABC'
hox1$type = 'exon'
hox1$score = 1000
hox1$phase = NA
hox1$gene_id = toupper(sapply(genes, function(x) unlist(strsplit(as.character(x), '_'))[1]))
hox1$transcript_id = hox1$gene_id

hoxs = c(hox1, hox2)
hoxs$gene_id = toupper(hoxs$gene_id)
hoxs$transcript_id = toupper(hoxs$transcript_id)


##########################################
# define the HOX A, B, C, D regions and remove everything from that regions and add hox cluster in that regions
# first for released version
##########################################
# HOXA chr2p:873,124,843-873,486,933
# HOXB chr13q:321,499,384-324,268,993
# HOXC chr3q:411,306,369-414,279,019
# HOXD chr9q:426,118,724-427,689,905
gr <- GRanges(
  seqnames = Rle(c("chr2p", "chr13q", "chr3q", "chr9q")),
  ranges = IRanges(c(873124843, 321499384, 411306369, 426118724), end = c(873486933, 324268993, 414279019, 427689905), 
                   names = c('HOXA', 'HOXB', 'HOXC', 'HOXD')),
  strand = Rle(strand(c("-", "+", "+", "-")), c(1, 1, 1, 1)))

gr

jj = findOverlaps(gr, aa, type = 'any', ignore.strand=TRUE)
#jj = findOverlaps(hoxs, aa, type = 'any')

xx = aa[jj@to, ]
export(xx, con = paste0(resDir, '/annotation_discarded_from_AmexT_v47_hoxPatch.gtf'), format = 'gtf')

xx = xx[grep('HOX', xx$transcript_id, invert = TRUE)]
xx = xx[grep('HOX', xx$gene_name, invert = TRUE)]

ggs = unique(xx$gene_id)
transcripts = unique(xx$transcript_id)

#ggs = paste0(ggs[, 1], '_',  ggs[, 2])
#ggs = ggs[grep('HOX', ggs, invert = TRUE)]

xx = aa[-jj@to]
xx[grep('^HOX', xx$gene_name)]
xx[grep('^HOX', xx$transcript_id)]

hoxs_a = hoxs
hoxs_a$gene_name = hoxs_a$gene_id
xx = c(hoxs_a, xx)
xx <- sortSeqlevels(xx)
xx <- sort(xx, ignore.strand=TRUE)

export(xx, con = paste0(resDir, '/AmexT_v48.Hox.gtf'), format = 'gtf')

##########################################
# remove all annotation in Hox clusters and add new hox genes
##########################################
jj = findOverlaps(gr, ucsc, type = 'any', ignore.strand=TRUE)

xx = ucsc[jj@to, ]
export(xx, con = paste0(resDir, '/annotation_discarded_from_AmexT_v47_UCSC_hoxPatch.gtf'), format = 'gtf')

xx = xx[grep('Hox', xx$transcript_id, invert = TRUE)]

ggs = unique(xx$gene_id)
transcripts = unique(xx$transcript_id)

#ggs = paste0(ggs[, 1], '_',  ggs[, 2])
#ggs = ggs[grep('HOX', ggs, invert = TRUE)]

xx = ucsc[-jj@to]

xx = c(xx, hoxs)
xx <- sortSeqlevels(xx)
xx <- sort(xx, ignore.strand=TRUE)

export(xx, con = paste0(resDir, '/AmexT_v48.Hox_UCSC.gtf'), format = 'gtf')
