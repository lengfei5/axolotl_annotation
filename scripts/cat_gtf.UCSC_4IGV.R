##########################################################################
##########################################################################
# Project: axolotl gene annotation
# Script purpose: concatenate gtf files downloaded from UCSC for IGV gene annotation
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Mar 25 12:24:21 2022
##########################################################################
##########################################################################
dataDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/UCSC_downloaded_14032022'

files = list.files(path = dataDir, pattern = '*.gtf', full.names = TRUE)

gtfs = c()
for(n in 1:length(files))
{
  # n = 1
  cat(n, ' -- ', files[n], '\n')
  tmp = read.delim(files[n], sep = '\t', header = FALSE)
  gtfs = rbind(gtfs, tmp)
   
}

gtfs = data.frame(gtfs, stringsAsFactors = FALSE)

## remove contigs 
xx = gtfs[grep('^chr', gtfs$V1), ]
xx$V1 = as.character(xx$V1)
gtfs = xx

gtfs = gtfs[order(gtfs$V1,  gtfs$V4), ]


xx$transcript = sapply(xx$V9, function(x) {
  temp = gsub('[|]', ' ', unlist(strsplit(as.character(x), ';'))[2])
  temp = unlist(strsplit(as.character(temp), ' '))
  temp[length(temp)]
})

write.table(gtfs, file = paste0(dataDir, '/ax6_v47_scafford.gtf'), sep = '\t', row.names = FALSE, col.names = FALSE,
            quote = FALSE)

