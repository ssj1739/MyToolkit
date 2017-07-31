main <- function(args){
  #gene <- as.character(args[0])
  gene <- "WDR77"
  load("~/Google Drive/Lab Files/WDR77-PRMT5/data/ucsc_annot.RData")
  tissue <- get_tissue(args[1])
  df <- read_delim(args[1], delim='\t', escape_double = F, trim_ws = T, progress = F)
  #print(nrow(r))
  isocount <- countIsos(df, gene, tissue, ucsc_annot)
  r <- cbind(isocount, tissue)
  f <- file(paste0("~/Documents/",tissue,"_",gene,"_isoform_results.txt"),'w')
  write.table(r, file=f, sep='\t', row.names = F, col.names = F)
  close(f)
  #print(args)
  #print(args[1])
}

get_tissue <- function(file){
  s =strsplit(file, split = '/', fixed=T)
  n=length(s[[1]])
  t=strsplit(s[[1]][n], split='.', fixed=T)
  return(t[[1]][1])
}

trim_ids <- function(id){
  return(strsplit(id, split = '.', fixed=T)[[1]][1])
}

mergeFiles <- function(dir, file_patt){
  f.list <- dir(path = dir, pattern=file_patt, full.names=T)
  file.create(file_patt)
  for(f in f.list){
    file.append(file1 = file_patt, file2=f)
  }
}

## Functionalizing the above analysis
countIsos <- function(df, gene, tissue, ucsc_annot, scale=F){
  # feed in path to dataset containing isoforms from TCGA
  # feed in gene name of interest to search for all isoforms
  # feed in object containing ucsc_annotations with at least one column called "ucsc_ids" that contains trimmed ids.
  isoform_data <- as.data.frame(df)
  isoform_ids <- sapply(isoform_data[,1], trim_ids) # process ids for matching with ucsc_annotations
  found <- isoform_ids %in% ucsc_annot$ucsc_ids # conduct matching
  #browser()
  geneNames <- as.data.frame(cbind(ucsc_annot[match(isoform_ids, ucsc_annot$ucsc_ids),], isoform_ids, found)) # subset matched genes from annotation
  
  # identify isoforms of genes of interest
  n.gene <- which(geneNames$hg19.ensemblToGeneName.value==gene)
  gene_isoforms <- isoform_data[n.gene,] #
  
  # Extract scaled or raw counts
  s=2
  if(scale){s=3}
  iso_counts <- as.data.frame(gene_isoforms[,seq(s, ncol(gene_isoforms), by=2)], row.names = gene_isoforms$`Hybridization REF`)
  iso_counts_t <- data.frame(apply(t(iso_counts),2,as.numeric), Sample=colnames(iso_counts))
  
  # Plot, save fig
  pdf(paste0("~/Google Drive/Lab Files/WDR77-PRMT5/figs/",tissue,"_",gene,".pdf"))
  barplot((colMeans(iso_counts_t[,1:ncol(iso_counts_t)-1])),
          col=1:ncol(iso_counts_t)*10,
          main=paste0("Log Scaled Expression of ",gene, " Isoforms\n in ",tissue, " samples - TCGA")
  )
  dev.off()
  return(iso_counts_t)
}  

library("dplyr", quietly = T, warn.conflicts = F)
library("readr", quietly = T, warn.conflicts = F)

main(commandArgs(TRUE))




