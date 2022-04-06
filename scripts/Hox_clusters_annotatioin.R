##########################################################################
##########################################################################
# Project:
# Script purpose: take the manually correction from Akane and liftover from contigs to scafford
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Apr  5 11:02:22 2022
##########################################################################
##########################################################################
dataDir = '/Volumes/groups/tanaka/Collaborations/Jingkui-Akane/amHoxABCseqences' 
resDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/axolotl_annotation/results'

##########################################
# HoxA cluster
##########################################
cluster = 'HoxAs'
aa = openxlsx::read.xlsx(paste0(dataDir, '/Am', cluster, '_Predicted_AK.xlsx'))
kk = which(colnames(aa) == "Contigs" | colnames(aa) ==  "Comment" )
aa = aa[, c(1, kk[1]:kk[2])]

aa = aa[which(!is.na(aa$Contigs)), ]

res = c()
for(n in 2:nrow(aa))
{
  # n = 6
  chr = aa[n, 2]
  strand = aa[n, 7]
  coords = as.character(aa[n, 3:6])
  coords = coords[!is.na(coords)]
  name = aa[n, 1]
  name = gsub('[*]', '', name)
  
  for(ii in 1:length(coords))
  {
    test = unlist(strsplit(as.character(gsub(':', '-', coords[ii])), '-'))
    if(length(test) == 3){
      res = rbind(res, c(paste0('C00', test[1]), test[2:3], strand, name) )
    }
    if(length(test) == 2){
      res = rbind(res, c(chr, test, strand, name))
    }
  }
}

res = data.frame(res, stringsAsFactors = FALSE)
res$score = 0
res$name = paste0(res$X5, '_', res[,2], '_', res[, 3])
res = res[, c(1:3, 7, 6, 4)]

res = res[match(unique(res$name), res$name), ]
jj = which(res$X2 > res$X3)
starts = res$X2[jj]
res$X2[jj] = res$X3[jj]
res$X3[jj] = starts

res = res[order(as.numeric(res$X2)), ]

write.table(res, file = paste0(resDir, '/Clusters_', cluster, '_contigs_coordinate.bed'), sep = '\t', quote = FALSE, row.names = FALSE,
            col.names = FALSE)


###### double-check after blasting the axolotl gene sequence
bst = read.table(file = paste0(resDir, '/AmHoxAs_Predicted_AK.outfmt6'), sep = '\t', header = FALSE)
bst = bst[which(bst$V2 == 'chr2p'), ]
bst = data.frame(bst, stringsAsFactors = FALSE)
jj = which(bst$V9 > bst$v10)
starts = bst$V9[jj]
bst$V9[jj] = bst$V10[jj]
bst$V10[jj] = starts

bst = bst[order(-bst$V9), ]
bst = bst[, -c(3:8)]

#### clean the liftovered HOXA clusters
res = read.table(file = paste0(resDir, '/HoxA_cluster_chr.bed'), sep = '\t', header = FALSE)

# discard the third exon of HOXA11
res = res[which(res$V4 != 'HoxA11_103792_104289'), ]
# res = res[which(res$V4 != 'HoxA3_105291_105417'), ]
write.table(res, file = paste0(resDir, '/HoxA_cluster_chr_cleaned.bed'), sep = '\t', col.names = FALSE, row.names = FALSE,
            quote = FALSE)

##########################################
# HoxB cluster 
##########################################
cluster = 'HoxBs'

aa = openxlsx::read.xlsx(paste0(dataDir, '/Am', cluster, '_Predicted_AK.xlsx'))
kk = which(colnames(aa) == "Contigs" | colnames(aa) ==  "Comment" )
aa = aa[, c(1, kk[1]:kk[2])]

aa = aa[which(!is.na(aa$Contigs)), ]

res = c()
for(n in 2:nrow(aa))
{
  # n = 6
  chr = aa[n, 2]
  strand = aa[n, 7]
  coords = as.character(aa[n, 3:6])
  coords = coords[!is.na(coords)]
  name = aa[n, 1]
  name = gsub('[*]', '', name)
  
  for(ii in 1:length(coords))
  {
    test = unlist(strsplit(as.character(gsub(':', '-', coords[ii])), '-'))
    if(length(test) == 3){
      res = rbind(res, c(paste0('C00', test[1]), test[2:3], strand, name) )
    }
    if(length(test) == 2){
      res = rbind(res, c(chr, test, strand, name))
    }
  }
}

res = data.frame(res, stringsAsFactors = FALSE)
res$score = 0
res$name = paste0(res$X5, '_', res[,2], '_', res[, 3])
res = res[, c(1:3, 7, 6, 4)]

res = res[match(unique(res$name), res$name), ]
jj = which(res$X2 > res$X3)
starts = res$X2[jj]
res$X2[jj] = res$X3[jj]
res$X3[jj] = starts

res = res[order(as.numeric(res$X2)), ]

write.table(res, file = paste0(resDir, '/Clusters_', cluster, '_contigs_coordinate.bed'), sep = '\t', quote = FALSE, row.names = FALSE,
            col.names = FALSE)


### clean the liftovered annotation
res = read.table(file = paste0(resDir, '/HoxB_cluster_chr.bed'), sep = '\t', header = FALSE)

res = res[order(res$V2), ]
# discard the third exon of HOXA11
#res = res[which(res$V4 != 'HoxA11_103792_104289'), ]
#res = res[which(res$V4 != 'HoxA3_105291_105417'), ]

write.table(res, file = paste0(resDir, '/HoxB_cluster_chr_cleaned.bed'), sep = '\t', col.names = FALSE, row.names = FALSE,
            quote = FALSE)


##########################################
# HOXC 
##########################################
cluster = 'HoxCs'

aa = openxlsx::read.xlsx(paste0(dataDir, '/Am', cluster, '_Predicted_AK.xlsx'))
kk = which(colnames(aa) == "Contigs" | colnames(aa) ==  "Comment" )
aa = aa[, c(1, kk[1]:kk[2])]

aa = aa[which(!is.na(aa$Contigs)), ]

res = c()
for(n in 2:nrow(aa))
{
  # n = 6
  chr = aa[n, 2]
  strand = aa[n, 7]
  coords = as.character(aa[n, 3:6])
  coords = coords[!is.na(coords)]
  name = aa[n, 1]
  name = gsub('[*]', '', name)
  
  for(ii in 1:length(coords))
  {
    test = unlist(strsplit(as.character(gsub(':', '-', coords[ii])), '-'))
    if(length(test) == 3){
      res = rbind(res, c(paste0('C00', test[1]), test[2:3], strand, name) )
    }
    if(length(test) == 2){
      res = rbind(res, c(chr, test, strand, name))
    }
  }
}

res = data.frame(res, stringsAsFactors = FALSE)
res$score = 0
res$name = paste0(res$X5, '_', res[,2], '_', res[, 3])
res = res[, c(1:3, 7, 6, 4)]

res = res[match(unique(res$name), res$name), ]
jj = which(res$X2 > res$X3)
starts = res$X2[jj]
res$X2[jj] = res$X3[jj]
res$X3[jj] = starts

res = res[order(as.numeric(res$X2)), ]

write.table(res, file = paste0(resDir, '/Clusters_', cluster, '_contigs_coordinate.bed'), sep = '\t', quote = FALSE, row.names = FALSE,
            col.names = FALSE)

### clean the liftovered annotation
res = read.table(file = paste0(resDir, '/HoxC_cluster_chr.bed'), sep = '\t', header = FALSE)

res = res[order(-res$V2), ]
# discard the third exon of HOXA11
#res = res[which(res$V4 != 'HoxA11_103792_104289'), ]
#res = res[which(res$V4 != 'HoxA3_105291_105417'), ]

write.table(res, file = paste0(resDir, '/HoxC_cluster_chr_cleaned.bed'), sep = '\t', col.names = FALSE, row.names = FALSE,
            quote = FALSE)

