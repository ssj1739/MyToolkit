### My toolkit of semi-useful functions that I use often


# FILE IO FUNCTIONS -------------------------------------------------------


#'@name read_gct
#'@author Sidharth Jain
#'@usage read_gct(file=path/to/file.gct)
#'@details Uses readr package to read in gct files
#'@param file Path to a gct file
#'@return the tibble object containing all information from the .gct file
#'
read_gct <- function(file){
  library(readr)
  gct <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE, skip = 2)
}

read_raw_tcga <- function(file, scaled=T){
  TCGA <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # Getting rownames (genes) and colnames (samples) ### Following CCLE format
  samples <- unique(sapply(colnames(TCGA)[-1], strsplit2, split="_", n=1))
  unknown_genes <- grepl("?", TCGA$`Hybridization REF`, fixed=T)
  genes <- sapply(TCGA$`Hybridization REF`, strsplit2, split="|", fixed=T, n=1)
  genes[unknown_genes] <- TCGA$`Hybridization REF`[unknown_genes]
  
  # Getting numeric data to fill matrix
  ii = grep("*_1", colnames(TCGA))
  if(scaled){
    dat = as.data.frame(TCGA[-1,ii])
  }else{
    dat = as.data.frame(TCGA[-1,(ii-1)])
  }
  
  colnames(dat) <- samples
  rownames(dat) <- deduplicate(genes[-1])
  
  return(dat)
}

read_norm_tcga <- function(file){
  if(!file.exists(file)){
    stop("File not found")
  }
  # if(grepl('tar.gz', file)){
  #   dir=strsplit2(file,split='.', fixed=T)
  #   dir=paste(dir[1:(length(dir)-2)], collapse = '.')
  #   
  #   unzipped_files <- untar(tarfile = file,compressed = T, list=T)
  #   untar(tarfile = file, compressed = T, verbose=T, exdir = dir)
  #  
  #   file=paste0(dir,'/',unzipped_files[grepl('data.txt', unzipped_files)])
  # }
  TCGA <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE)
  
  samples <- colnames(TCGA)[-1]
  unknown_genes <- grepl("?", TCGA$`Hybridization REF`, fixed=T)
  genes <- sapply(TCGA$`Hybridization REF`, strsplit2, split="|", fixed=T, n=1)
  genes[unknown_genes] <- TCGA$`Hybridization REF`[unknown_genes]
  
  dat = as.data.frame(TCGA[-1,-1])
  colnames(dat) <- samples
  rownames(dat) <- deduplicate(genes[-1])
  
  return(dat)
}

read_clinical <- function(file){
  library(readr)
  clinical <- as.data.frame(t(read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE)))
  colnames(clinical) <- sapply(clinical[1,], as.character)
  clinical$Tumor_Sample_Barcode <- rownames(clinical)
  clinical <- as.data.frame(clinical[-1,])
  return(clinical)
}

write_gct <- function(X, file){
  library(readr)
  write_lines(c("#1.2", paste0(nrow(X), "\t", ncol(X)-2)), path = file)
  write_delim(x=X, path=file, delim="\t", append=T, col_names = T)
}

loadConstants <- function(){
  load(file = '~/Data/constants.RData', envir = .GlobalEnv)
}

loadData <- function(type, gene=NULL, cell_line=NULL, clean_features=T, clean_celllines=T, convert_id=F){
  library("taigr")
  library('stringr')
  data <- switch(type,
                 CRISPR = load.from.taiga(data.name='avana-broad-17q4-2d4b', data.file='gene_effect', quiet = T),
                 Avana = load.from.taiga(data.name='avana-broad-17q4-2d4b', data.file='gene_effect', quiet = T),
                 Broad_RNAi = t(load.from.taiga(data.name='achilles-demeter2-f87d', quiet = T)),
                 Achilles = t(load.from.taiga(data.name='achilles-demeter2-f87d', quiet = T)),
                 DRIVE = t(load.from.taiga(data.name='drive-ataris-scores-f72c', data.file='expanded_ATARiS_data', quiet = T)),
                 Novartis_RNAi = t(load.from.taiga(data.name='drive-ataris-scores-f72c', data.file='expanded_ATARiS_data', quiet = T)),
                 expression = load.from.taiga(data.name='ccle-rnaseq-expression-genes', quiet = T),
                 mutation = load.from.taiga(data.name='pooled-mutation-6481', quiet = T),
                 mut = load.from.taiga(data.name='pooled-mutation-6481', quiet = T),
                 mut_damaging = load.from.taiga(data.name='damaging-729f', quiet = T),
                 cn = load.from.taiga(data.name='gene-level-cn-87aa', data.file='gene_CN_SNP_priority', quiet = T),
                 copynumber = load.from.taiga(data.name='gene-level-cn-87aa', data.file='gene_CN_SNP_priority', quiet = T),
                 RNAseq = load.from.taiga(data.name='ccle-rnaseq-expression-genes', quiet = T),
                 DRIVE_RSA = t(load.from.taiga(data.name='drive-rsa-scores-5668', quiet = T)),
                 drive_rsa = t(load.from.taiga(data.name='drive-rsa-scores-5668', quiet = T)),
                 DRIVE_rsa = t(load.from.taiga(data.name='drive-rsa-scores-5668', quiet = T)),
                 chromatin = load.from.taiga(data.name='ccle-global-chromatin-prof-f49b', quiet = T),
                 emt = load.from.taiga(data.name='ccle-emt-score-b0da', quiet = T),
                 EMT = load.from.taiga(data.name='ccle-emt-score-b0da', quiet = T),
                 GSEA = t(load.from.taiga(data.name='ssgsea-enrichment-scores-for-msigdb-h-using-ccle-rnaseq-expression', quiet = T)),
                 gsea_hallmarks = t(load.from.taiga(data.name='ssgsea-enrichment-scores-for-msigdb-h-using-ccle-rnaseq-expression', quiet = T)),
                 methylation = load.from.taiga(data.name='rrbs-4b29', quiet = T),
                 rrbs = load.from.taiga(data.name='rrbs-4b29', quiet = T),
                 metabolomics = load.from.taiga(data.name='metabolomics-cd0c', quiet = T),
                 metabolome = load.from.taiga(data.name='metabolomics-cd0c', quiet = T),
                 demeter2 = t(load.from.taiga(data.name='demeter2-combined-dc9c', data.file='gene_means_proc', quiet = T)),
                 Combined_DEMETER2 = t(load.from.taiga(data.name='demeter2-combined-dc9c', data.file='gene_means_proc', quiet = T)),
                 Combined = t(load.from.taiga(data.name='demeter2-combined-dc9c', data.file='gene_means_proc', quiet = T)),
                 Combined_RNAi = t(load.from.taiga(data.name='demeter2-combined-dc9c', data.file='gene_means_proc', quiet = T)),
                 demeter2_combined =  t(load.from.taiga(data.name='demeter2-combined-dc9c', data.file='gene_means_proc', quiet = T)),
                 Achilles_DEMETER2 = t(load.from.taiga(data.name='demeter2-achilles-5386', data.file='gene_means_proc', quiet = T)),
                 DRIVE_DEMETER2 = t(load.from.taiga(data.name='demeter2-drive-0591', data.file='gene_means_proc')),
                 Position = load.from.taiga(data.name='gene-cytogenic-position-2b92'),
                 GDSC = load.from.taiga(data.name='gdsc-ic50-ccle-aligned', data.file='data'),
                 IC50 = load.from.taiga(data.name='gdsc-ic50-ccle-aligned', data.file='data'),
                 miRNA = load.from.taiga(data.name='mirna-expression-2c5f', data.version=2)
  )
  force(data)
  if(convert_id || length(grep("DEMETER2", ignore.case = T, type))==1 || length(grep("combined", ignore.case = T, type))==1 ){
    colnames(data) <- as.character(sapply(colnames(data), strsplit2, n=1, split="[&, ]"))
  }
  if(clean_features){
    colnames(data) <- sapply(colnames(data), strsplit2, split=" ", n=1)
  }
  if(clean_celllines){
    rownames(data) <- sapply(rownames(data), str_to_upper)
  }
  if(!is.null(cell_line)){
    if(!any(grep(cell_line, rownames(data))>0)){
      warning(paste0("Cell line(s) ", cell_line, "is not in dataset ", type, ".  Returned NULL instead."))
      return(NULL)
    }
    if(length(cell_line) > 1){
      data <- data[match(cell_line, rownames(data)),]
    }else{
      data <- data[grep(cell_line, rownames(data)),]
    }
    if(length(data)==0 & !is.null(data)){
      warning("No cell lines match query.  Returned NULL instead.")
      return(NULL)
    }
  }
  if(!is.null(gene)){
    if(!gene %in% colnames(data)){
      warning(paste0("Gene ", gene, " is not in dataset ", type, ".  Returned NULL instead."))
      return(NULL)
    }
    data <- data[,colnames(data) %in% gene]
    return(data)
  }
  return(as.data.frame(data))
}

cell_line_status <- function(cell_line, type="RNAseq", gene=NULL){
  data=loadData(type = type)
  ii = grepl(pattern = cell_line, rownames(data), ignore.case = T)
  if(sum(ii)==0){
    stop("Cell line not found!")
  }
  if(!gene %in% colnames(data)){
    stop("Invalid gene!")
  }
  return(data[ii, gene])
}

pullCBPdata <- function(genes, profile_type=c("EXPR", "COPY", "MUT")[1]){
  library(cgdsr)
  library(tidyverse)
  library(magrittr)
  
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  cancerstudies <- getCancerStudies(mycgds)[,1]
  ii = sapply(cancerstudies, function(x){
    prof <- getGeneticProfiles(mycgds, x)
    return(prof[which(grepl(profile_type, prof$genetic_alteration_type)),1])
  })
  studies <- names(ii[sapply(ii, length) > 0])
  
  # Remove all tcga-pub studies, only use supplemented TCGA data
  studies <- studies[!grepl('tcga_pub', studies)]
  
  # Retrieve data for corresponding studies and cases in those studies.
  df_all <- lapply(studies, function(x){
    caselist = getCaseLists(mycgds, x)[1,1]
    df_list <- lapply(ii[[x]], function(y){
      df = getProfileData(x=mycgds, genes = genes, geneticProfiles = y, caseList = caselist)
      df <- as.data.frame(apply(df, 2, as.numeric))
      df$Sample <- rownames(df)
      df$Profile <- y
      return(as_tibble(df))
    })
    df_bound = bind_rows(df_list)
    df_bound$Study = x
    return(df_bound)
  })
  
  df_all %<>% 
    bind_rows() %<>%
    mutate("Type" = sapply(Study, strsplit2, split="_", n=1))
  
  return(df_all)
}

write_fasta <- function(file = "", seq, name, append=FALSE){
  require(readr)
  write_lines(path = file, x=c(paste(">", name), seq, "\n"), append = append)
}

read_mappings <- function(file, sep = "\t", IDsep = ",") {
  a <- read.delim(file = file, header = FALSE,
                  quote = "", sep = sep, colClasses = "character")
  
  ## a bit of preprocesing to get it in a nicer form
  map <- a[, 2]
  names(map) <- gsub(" ", "", a[, 1]) ## trim the spaces
  
  ## split the IDs
  return(lapply(map, function(x) gsub(" ", "", strsplit(x, split = IDsep)[[1]])))
}


# DATA WRANGLING FUNCTIONS ------------------------------------------------
#'@name deduplicate
#'@author Sidharth Jain
#'@description Deduplicates a vector - any values X in the vector that are the nth duplicate are replaced by X{sep}n
#'@example x <- c('A', 'A', 'B', 'B', 'C', 'C', 'C')
#'deduplicate(x) #[1] "A"   "A_0" "B"   "B_0" "C"   "C_0" "C_1"
#'
#'@
deduplicate <- function(x, sep="_"){
# Deduplicates a vector - any values X in the vector that are the nth duplicate are replaced by X{sep}n
# User should take note of any values that are already called X{sep}n - dangerous!
  if(!anyDuplicated(x))
    return(x)
  
  count = 0
  while(anyDuplicated(x)){
     first=which(duplicated(x))
     x[first] <- sapply(x[first], function(x) strsplit2(x, split=sep, n=1))
     x[first] <- sapply(x[first], function(x) paste0(x, sep, count))
     count = count + 1
  }
  return(x)
  return(x)
}

interpolate_matrix <- function(x){
  # Fills in NAs in a matrix with average of row mean and column mean.  
  cm=colMeans(x, na.rm=T)
  rm=rowMeans(x, na.rm=T)
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(is.na(x[i,j])){
        n=mean(rm[i], cm[j])
        x[i,j] = n
      }
    }
    if(i %% (nrow(x)/10) == 0){
      print(paste0(i/(nrow(x)) * 100, "% complete"))
    }
  }
  return(x)
}

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

reorder2 <- function(X, x.1=NULL, by){
  if(is.null(x.1) & is.null(dim(X))){
    return(X[order(match(X, by))])
  }
  
  return(X[order(match(x.1, by)),])
}

strsplit2 <- function(x, n, ...){
  # cleans name and returns nth argument, quick and easy for gene or cell line names
  strsplit(x, ...)[[1]][n]
}

clean_name <- function(x, ...){
  strsplit2(x, n=1, ...)
}

head2 <- function(x, n=5){
  # Only shows first n rows and n columns in 2D matrix/dataframe
  n=ifelse(nrow(x)>n,n,nrow(x))
  n=ifelse(ncol(x)>n,n,ncol(x))
  return(x[1:n,1:n])
}

top <- function(x, n=10, decreasing=T){
  sort(x, decreasing = decreasing)[1:n]
}

tsne.dist <- function(tsne.obj, i=NULL, gene_names=NULL){
  ndim=ncol(tsne.obj$Y)
  dat=tsne.obj$Y
  # Select int (interested rows) if mentioned in i argument
  if(length(i) >= 1){
    if(is.factor(i))
      i=as.character(i)
    if(is.character(i)){
      int=match(i, gene_names)
    }else{
      int=i
    }
  }else{
    int=1:nrow(dat)
    i=int
  }
  res=matrix(data = 0, nrow = length(int), ncol=nrow(dat))
  for(ii in 1:length(int)){
    x0=dat[int[ii],]
    if(all(is.na(x0))){
      res[ii,] <- NA
    }else{
      for(j in 1:nrow(dat)){
        x=dat[j,1:ndim]
        res[ii,j]=sqrt(sum((x0-x)^2))
      }
    }
  }
  colnames(res) <- gene_names
  return(data.frame(res, row.names = i))
}


intersect_recursive <- function(set_list){
  if(length(set_list)==2){
    return(intersect(unlist(set_list[1]), unlist(set_list[2])))
  }else{
    return(intersect(unlist(set_list[1]), intersect_recursive(set_list[-1])))
  }
}

union_recursive <- function(set_list){
  if(length(set_list)==2){
    return(union(unlist(set_list[1]), unlist(set_list[2])))
  }else{
    return(union(unlist(set_list[1]), union_recursive(set_list[-1])))
  }
}

convert_ids <- function(old_ids){
  new_mappings <- load.from.taiga(data.name='gpp-shrna-mapping-8759', data.version=1, data.file='CP1175_20171102_19mer')
  new_ids <- new_mappings$`Gene Symbol`[match(old_ids, new_mappings$`Gene ID`)]
  return(new_ids)
}

get_tissue <- function(cl){
  s=strsplit(cl, split = "_")[[1]][-1]
  return(paste(s, collapse=" "))
}

extract_info_tcga_sample_ids <- function(ids){
  id_breakdown = as.data.frame(t(sapply(ids, function(x) strsplit(x, split='-')[[1]])))
  colnames(id_breakdown) <- c("Project", "TSS", "Participant", "Sample/Vial", "Portion/Analyte", "Plate", "Center")
  id_breakdown$Phenotype <- "Normal"
  id_breakdown$Phenotype[grep('0[1-9]', id_breakdown$`Sample/Vial`)] <- "Tumor"
  return(id_breakdown)
}

`%,%` <- function(a, b){
  paste0(a, b)
}

# PLOTTING FUNCTIONS ------------------------------------------------------
dPlot <- function(geneX, geneY, data.set.x="CRISPR", data.set.y="Expression", geneContext=NULL, context=NULL){
  # Creates a plot of dependency vs anything, including expression, copy number, mutation, etc.
  
  # Load in data
  library("taigr")
  library("readr")
  library("ggplot2")
  
  data.x <- loadData(gene=geneX, type = data.set.x)
  data.y <- loadData(gene=geneY, type = data.set.y)
  
  cls <- intersect(names(data.x), names(data.y))
  tissue <- sapply(cls, get_tissue)
  
  if(!is.null(context)){
    annotation=context
  }else{
    annotation=tissue
  }
  
  df <- data.frame(cell_lines=cls, 
                   geneX, 
                   X=data.x[names(data.x) %in% cls], 
                   geneY, 
                   Y=data.y[names(data.y) %in% cls],
                   annotation = annotation[cls])
  g <- ggplot(data=df, aes(x=X, y=Y, color=annotation)) +
    geom_point() +
    labs(x=data.set.x, y=data.set.y) +
    theme(legend.position = 'none')
  
  library(plotly)
  
  return(ggplotly(g))
}

qplot3d <- function(mat3d, features=rownames(mat3d), axis.labels=colnames(mat3d), ...){
  # Quickly plots a 3 dimensional matrix
  # Does not add any colors, etc. by default other than those already included in plotly
  library(plotly)
  p <- plot_ly(data = as.data.frame(mat3d), x = ~mat3d[,1], y = ~mat3d[,2], z = ~mat3d[,3], ...) %>%
    #add_text() %>%
    add_markers(hovertext=features) %>%
    layout(
      scene = list(
        xaxis = list(title=colnames(mat3d)[1]),
        yaxis = list(title=colnames(mat3d)[2]),
        zaxis = list(title=colnames(mat3d)[3])
      )
    )
  return(p)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}

waterfallplot <- function(gene, cls=NULL, data.set="CRISPR", annotate.by="tissue", annotate.gene=NULL, label_bars=F, static=F, legend=T){
  library("taigr")
  library("ggplot2")
  library("plotly")
  stopifnot(length(gene)>0)
  
  gene_effect <- loadData(gene = gene, type=data.set)
  cls = names(gene_effect)
  tissues <- sapply(cls, function(x){ 
    s=strsplit(x, split = "_")[[1]][-1]
    return(paste(s, collapse=" "))
  })
  
  annotate.gene <- ifelse(is.null(annotate.gene), gene, annotate.gene)
  gene_annot <- NA
  if(annotate.by %in% c("mutation", "mut")){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    gene_annot <- ifelse(ccle_binary[match(paste0(annotate.gene, "_", "MUT"), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]==1, 1, 0) 
    
    
    #ccle_mut_del <- load.from.taiga(data.name='ccle-mut-data-binary-matrix', data.version=1)
    #ccle_mut_mis <- load.from.taiga(data.name='ccle-mis-mut-binary-matrix', data.version=1)
    # a <- ccle_mut_del[match( cls, rownames(ccle_mut_del)), match(paste0(annotate.gene, "_", "DEL"), colnames(ccle_mut_del))]
    # b <- ccle_mut_mis[match( cls, rownames(ccle_mut_mis)), match(gene, colnames(ccle_mut_mis))]
    # gene_annot <- as.logical(a+b)
  }else if(annotate.by %in% c("copy number", "cn", "cnv") ){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    
    gene_annot <- ifelse(ccle_binary[match(paste0(annotate.gene, "_", "AMP"), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]==1, 1, 0)
    # ccle_copynumber <- t(load.from.taiga(data.name='ccle-copy-number-variants-hgnc-mapped', data.version=4))
    # colnames(ccle_copynumber) <- sapply(colnames(ccle_copynumber), strsplit2, split=" ", n=1)
    # gene_annot <- round(ccle_copynumber[match(cls, rownames(ccle_copynumber)), match(gene, colnames(ccle_copynumber))])
  }else if(annotate.by %in% c("deletion", "del")){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    
    gene_annot <- ifelse(ccle_binary[match(paste0(annotate.gene, "_", "DEL"), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]==1, 1, 0)
    
  }else{
    gene_annot <- tissues
  }
  
  if(all(is.na(gene_annot))){
    warning("No known annotations for " %,% annotate.by %,% ".  Instead annotating by tissue.")
    gene_annot <- tissues
  }
  
  df <- data.frame(
    names=cls, 
    clean_names = sapply(cls, strsplit2, split="_", n=1),
    scores = gene_effect,
    context=as.factor(as.character(gene_annot))
  )
  
  g <- ggplot(df, aes(x=reorder(names, -scores), y=scores, fill=context))
  
  g <- g + geom_bar(stat = 'identity', na.rm = T, aes(text=sprintf("Cell line: %s<br>Score: %s<br>Context: %s", names, scores, context))) +
    scale_fill_discrete(drop=F) +
    labs(x="Cell line", y="Dependency") +
    theme(axis.text.x=element_text(angle=45, hjust=0, size=4), legend.text = element_text(size = 7)) +
    ggtitle(paste0("Waterfall Plot - ", gene, " in ", data.set))
  
  if(!legend)
    g <- g+theme(legend.position = 'none')
  if(!label_bars)
    g <- g+theme(axis.text.x = element_blank())
  if(static)
    return(g)
  return(ggplotly(g))  
}


drive_waterfallplot <- function(gene, cls=NULL, annotate.by='tissue', annotate.gene=NULL, static=F, label_bars=F, color_by_discrete=T){
  library("taigr")
  library("ggplot2")
  library("plotly")
  stopifnot(length(gene)>0)
  
  gene_effect <- loadData(gene = gene, type="DRIVE")
  gene_rsa <- loadData(gene = gene, type="DRIVE_RSA")
  cls = names(gene_effect)
  
  gene_annot.by_tissue=sapply(cls, function(x){ # By default, annotate by tissue/lineage
    s=strsplit(x, split = "_")[[1]][-1]
    return(paste(s, collapse=" "))
  })
  gene_annot <- gene_annot.by_tissue
  
  if(is.null(annotate.gene)){
    annotate.gene <- gene
  }
  
  
  if(annotate.by %in% c("mutation", "mut")){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    gene_annot <- ccle_binary[grep(glob2rx(paste0(annotate.gene, "*_MUT")), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
    
    
    #ccle_mut_del <- load.from.taiga(data.name='ccle-mut-data-binary-matrix', data.version=1)
    #ccle_mut_mis <- load.from.taiga(data.name='ccle-mis-mut-binary-matrix', data.version=1)
    # a <- ccle_mut_del[match( cls, rownames(ccle_mut_del)), match(paste0(annotate.gene, "_", "DEL"), colnames(ccle_mut_del))]
    # b <- ccle_mut_mis[match( cls, rownames(ccle_mut_mis)), match(gene, colnames(ccle_mut_mis))]
    # gene_annot <- as.logical(a+b)
  }else if(annotate.by %in% c("copy number", "cn", "cnv") ){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    
    gene_annot <- ccle_binary[grep(glob2rx(paste0(annotate.gene, "*_AMP")), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
    # ccle_copynumber <- t(load.from.taiga(data.name='ccle-copy-number-variants-hgnc-mapped', data.version=4))
    # colnames(ccle_copynumber) <- sapply(colnames(ccle_copynumber), strsplit2, split=" ", n=1)
    # gene_annot <- round(ccle_copynumber[match(cls, rownames(ccle_copynumber)), match(gene, colnames(ccle_copynumber))])
  }else if(annotate.by %in% c("deletion", "del")){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    
    gene_annot <- ccle_binary[grep(glob2rx(paste0(annotate.gene, "*_DEL")), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
    
  }else if(annotate.by %in% c("RNAseq", "expr", "expression")){
    expr <- loadData("RNAseq", gene=annotate.gene)
    gene_annot <- expr[match(cls, names(expr), nomatch=NA)]
    color_by_discrete=F
  }else if(annotate.by %in% c('any', 'all')){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    gene_annot <- ccle_binary[grep(paste0(annotate.gene,"[[:punct:]]"), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch=1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
  }
  
  if(all(is.na(gene_annot))){
    warning("No known annotations for " %,% annotate.by %,% ".  Instead annotating by tissue.")
    gene_annot <- as.factor(as.character(gene_annot.by_tissue))
  }
  df <- data.frame(
    names=cls, 
    clean_names = sapply(cls, strsplit2, split="_", n=1),
    ataris_scores = gene_effect,
    rsa_scores = gene_rsa,
    context=as.factor(as.character(gene_annot))
  )
  
  g1 <- ggplot(df, aes(x=reorder(names, -ataris_scores), y=ataris_scores, fill=context)) +
    geom_bar(stat = 'identity', 
             na.rm = T, 
             aes(text=sprintf("Cell line: %s<br>Score: %s<br>Context: %s", names, ataris_scores, context))) +
    labs(x="", y="ATARiS Score") +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=4), legend.position = 'none')
  
  g2 <- ggplot(df, aes(x=reorder(names, -ataris_scores), y=rsa_scores, fill=context)) +
    geom_bar(stat = 'identity', 
             na.rm = T, 
             aes(text=sprintf("Cell line: %s<br>Score: %s<br>Context: %s", names, rsa_scores, context))) +
    labs(x="", y="RSA Score") +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=4), legend.position='none')
  
  p <- lapply(list(g1, g2), function(x){
    if(label_bars){
      x <- x+geom_text(aes(label=df$clean_names, size=0.5, hjust=1, vjust=1, angle=45))
    }
    if(color_by_discrete){
      x <- x+scale_fill_discrete(drop=F)
    }
    return(x)
  })
  
  if(static){
    return(multiplot(p))
  }
  
  p <- lapply(p, ggplotly, tooltip='text')
  
  return(subplot(p[1:2], titleY=TRUE, nrows=2, shareX=T))
}

comparison_waterfallplot <- function(gene, annotate.by="tissue", annotate.gene=NULL, annotate.custom=NULL, label_bars=F, reorder.by="DRIVE", static=F, color_by_discrete=T, legend=T){
  library("taigr")
  library('readr')
  library("ggplot2")
  library("plotly")
  stopifnot(length(gene)>0)
  
  gene_effect_avana <- loadData(type = "Avana")
  gene_effect_achilles <- loadData(type="Achilles")
  gene_effect_drive <- loadData(type="DRIVE")
  gene_effect_combined <- loadData(type="Combined")
  
  which.plots <- c(Avana=gene %in% colnames(gene_effect_avana), Achilles=gene %in% colnames(gene_effect_achilles), DRIVE=gene %in% colnames(gene_effect_drive), Combined=gene %in% colnames(gene_effect_combined))
  
  
  cls <- switch(reorder.by,
                DRIVE=rownames(gene_effect_drive),
                Achilles=rownames(gene_effect_achilles),
                Avana=rownames(gene_effect_avana),
                CRISPR=rownames(gene_effect_avana),
                Combined=rownames(gene_effect_combined))
  
  nplots <- sum(which.plots)
  
  gene_annot.by_tissue=sapply(cls, function(x){ # By default, annotate by tissue/lineage
    s=strsplit(x, split = "_")[[1]][-1]
    return(paste(s, collapse=" "))
  })
  gene_annot <- gene_annot.by_tissue
  
  if(is.null(annotate.gene)){
    annotate.gene <- gene
  }
  
  
  if(annotate.by %in% c("mutation", "mut")){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    gene_annot <- ccle_binary[grep(glob2rx(paste0(annotate.gene, "*_MUT")), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
    
    
    #ccle_mut_del <- load.from.taiga(data.name='ccle-mut-data-binary-matrix', data.version=1)
    #ccle_mut_mis <- load.from.taiga(data.name='ccle-mis-mut-binary-matrix', data.version=1)
    # a <- ccle_mut_del[match( cls, rownames(ccle_mut_del)), match(paste0(annotate.gene, "_", "DEL"), colnames(ccle_mut_del))]
    # b <- ccle_mut_mis[match( cls, rownames(ccle_mut_mis)), match(gene, colnames(ccle_mut_mis))]
    # gene_annot <- as.logical(a+b)
  }else if(annotate.by %in% c("copy number", "cn", "cnv") ){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    
    gene_annot <- ccle_binary[grep(glob2rx(paste0(annotate.gene, "*_AMP")), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
    # ccle_copynumber <- t(load.from.taiga(data.name='ccle-copy-number-variants-hgnc-mapped', data.version=4))
    # colnames(ccle_copynumber) <- sapply(colnames(ccle_copynumber), strsplit2, split=" ", n=1)
    # gene_annot <- round(ccle_copynumber[match(cls, rownames(ccle_copynumber)), match(gene, colnames(ccle_copynumber))])
  }else if(annotate.by %in% c("deletion", "del")){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    
    gene_annot <- ccle_binary[grep(glob2rx(paste0(annotate.gene, "*_DEL")), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
    
  }else if(annotate.by %in% c("RNAseq", "expr", "expression")){
    expr <- loadData("RNAseq", gene=annotate.gene)
    gene_annot <- expr[match(cls, names(expr), nomatch=NA)]
    color_by_discrete=F
  }else if(annotate.by %in% c('any', 'all')){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    gene_annot <- ccle_binary[grep(paste0(annotate.gene,"[[:punct:]]"), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch=1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
  }else if(annotate.by %in% c('DRIVE', 'Achilles', 'CRISPR')){
    dep = loadData(type=annotate.by, gene=annotate.gene)
    gene_annot <- dep[match(cls, names(dep), nomatch = NA)]
  }
  
  if(all(is.na(gene_annot))){
    warning("No known annotations for " %,% annotate.by %,% ".  Instead annotating by tissue.")
    gene_annot <- as.factor(as.character(gene_annot.by_tissue))
  }
  
  drive_scores = gene_effect_drive[match(cls, rownames(gene_effect_drive)),match(gene, colnames(gene_effect_drive))]
  avana_scores=gene_effect_avana[match(cls, rownames(gene_effect_avana)),match(gene, colnames(gene_effect_avana))]
  achilles_scores=gene_effect_achilles[match(cls, rownames(gene_effect_achilles)),match(gene, colnames(gene_effect_achilles))]
  combined_scores=gene_effect_combined[match(cls, rownames(gene_effect_combined)),match(gene, colnames(gene_effect_combined))]
  
  df <- data.frame(
    names=cls, 
    clean_names = sapply(cls, strsplit2, split="_", n=1),
    drive_scores=ifelse(is.null(drive_scores), NA, drive_scores), 
    avana_scores=ifelse(is.null(avana_scores), NA, avana_scores),
    achilles_scores=ifelse(is.null(achilles_scores), NA, achilles_scores),
    combined_scores=ifelse(is.null(combined_scores), NA, combined_scores),
    context=gene_annot
  )
  
  df$reorder.by <-  switch(reorder.by, 
                           DRIVE=df$drive_scores,
                           Avana=df$avana_scores,
                           CRISPR=df$avana_scores,
                           Achilles=df$achilles_scores,
                           Combined=df$combined_scores)
  
  iplots <- which(apply(df[,c('drive_scores', 'avana_scores', "achilles_scores", "combined_scores")], 2, function(x) !all(is.na(x))))
  
  g1 <- ggplot(df, aes(x=reorder(names, -reorder.by), y=drive_scores, fill=context)) +
    geom_bar(stat = 'identity', 
             na.rm = T, 
             aes(text=sprintf("Cell line: %s<br>Score: %s<br>Context: %s", names, drive_scores, context))) +
    labs(x="", y="ATARiS Score\n") +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=4), legend.position = 'none')
  
  g2 <- ggplot(df, aes(x=reorder(names, -reorder.by), y=avana_scores, fill=context)) +
    geom_bar(stat = 'identity', 
             na.rm = T, 
             aes(text=sprintf("Cell line: %s<br>Score: %s<br>Context: %s", names, avana_scores, context))) +
    labs(x="", y="CERES Score<br>") +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=4), legend.position='none')
  
  g3 <- ggplot(df, aes(x=reorder(names, -reorder.by), y=achilles_scores, fill=context)) +
    geom_bar(stat = 'identity', 
             na.rm = T, 
             aes(text=sprintf("Cell line: %s<br>Score: %s<br>Context: %s", names, achilles_scores, context))) +
    labs(x="", y="Achilles D1 Score<br>") +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=4), legend.position='none')
  
  g4 <- ggplot(df, aes(x=reorder(names, -reorder.by), y=combined_scores, fill=context)) +
    geom_bar(stat = 'identity', 
             na.rm = T, 
             aes(text=sprintf("Cell line: %s<br>Score: %s<br>Context: %s", names, combined_scores, context))) +
    labs(x="", y="Combined D2 Score\n", title=paste0("Dependency Comparison - ", gene)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=4), legend.text = element_text(size = 7), legend.title=element_blank())
  
  p <- lapply(list(g1, g2, g3, g4), function(x){
    if(label_bars){
      x <- x+geom_text(aes(label=df$clean_names, size=0.5, hjust=1, vjust=1, angle=45))
    }
    if(color_by_discrete){
      x <- x+scale_fill_discrete(drop=F)
    }
    if(!legend){
      x <- x+theme(legend.position = 'none')
    }
    return(x)
  })
  
  if(static){
    return(multiplot(p))
  }
  
  p <- lapply(p, ggplotly, tooltip='text')
  
  return(subplot(p[iplots], titleY=TRUE, nrows=nplots, shareX=T))
}

comparison_waterfallplot_expression <- function(gene, ...){
  plotlist <- list(dPlot("AGO2", gene, "RNAseq", "DRIVE"),
                   dPlot("AGO2", gene, "RNAseq", "CRISPR"),
                   dPlot("AGO2", gene, "RNAseq", "Achilles"))
  p <- subplot(plotlist, nrows = 3)
  subplot(comparison_waterfallplot(gene, ...), p, titleY = T)
}

comparison_dataframe <- function(gene, annotate.by="tissue", annotate.gene=NULL, reorder.by="DRIVE"){
  library("taigr")
  library('readr')
  stopifnot(length(gene)>0)
  
  gene_effect_avana <- loadData(type = "Avana", gene=gene)
  gene_effect_achilles <- loadData(type="Achilles", gene=gene)
  gene_effect_drive <- loadData(type="DRIVE", gene=gene)
  gene_effect_combined <- loadData(type="Combined", gene=gene)
  
  cls <- switch(reorder.by,
                DRIVE=names(gene_effect_drive),
                Achilles=names(gene_effect_achilles),
                Avana=names(gene_effect_avana),
                CRISPR=names(gene_effect_avana),
                Combined=names(gene_effect_combined))
  
  #nplots <- gene %in% colnames(gene_effect_avana) + gene %in% colnames(gene_effect_achilles) + gene %in% colnames(gene_effect_drive) + gene %in% colnames(gene_effect_combined)
  
  gene_annot.by_tissue=sapply(cls, function(x){ # By default, annotate by tissue/lineage
    s=strsplit(x, split = "_")[[1]][-1]
    return(paste(s, collapse=" "))
  })
  gene_annot <- gene_annot.by_tissue
  
  if(is.null(annotate.gene)){
    annotate.gene <- gene
  }
  
  
  if(annotate.by %in% c("mutation", "mut")){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    gene_annot <- ccle_binary[grep(glob2rx(paste0(annotate.gene, "*_MUT")), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
    
    
    #ccle_mut_del <- load.from.taiga(data.name='ccle-mut-data-binary-matrix', data.version=1)
    #ccle_mut_mis <- load.from.taiga(data.name='ccle-mis-mut-binary-matrix', data.version=1)
    # a <- ccle_mut_del[match( cls, rownames(ccle_mut_del)), match(paste0(annotate.gene, "_", "DEL"), colnames(ccle_mut_del))]
    # b <- ccle_mut_mis[match( cls, rownames(ccle_mut_mis)), match(gene, colnames(ccle_mut_mis))]
    # gene_annot <- as.logical(a+b)
  }else if(annotate.by %in% c("copy number", "cn", "cnv") ){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), trim_ws = TRUE, skip = 2))
    
    gene_annot <- ccle_binary[grep(glob2rx(paste0(annotate.gene, "*_AMP")), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
    # ccle_copynumber <- t(load.from.taiga(data.name='ccle-copy-number-variants-hgnc-mapped', data.version=4))
    # colnames(ccle_copynumber) <- sapply(colnames(ccle_copynumber), strsplit2, split=" ", n=1)
    # gene_annot <- round(ccle_copynumber[match(cls, rownames(ccle_copynumber)), match(gene, colnames(ccle_copynumber))])
  }else if(annotate.by %in% c("deletion", "del")){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    
    gene_annot <- ccle_binary[grep(glob2rx(paste0(annotate.gene, "*_DEL")), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch = 1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
    
  }else if(annotate.by %in% c("RNAseq", "expr", "expression")){
    expr <- loadData("RNAseq", gene=annotate.gene)
    gene_annot <- expr[match(cls, names(expr), nomatch=NA)]
    color_by_discrete=F
  }else if(annotate.by %in% c('any', 'all')){
    library(readr)
    ccle_binary <- as.data.frame(read_delim("~/Data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct", 
                                            "\t", escape_double = FALSE, col_types = cols(Description = col_skip()), 
                                            trim_ws = TRUE, skip = 2))
    gene_annot <- ccle_binary[grep(paste0(annotate.gene,"[[:punct:]]"), ccle_binary$Name),match(cls, colnames(ccle_binary), nomatch=1)]
    gene_annot <- apply(gene_annot, 2, function(x) any(as.logical(as.numeric(x))))
  }else if(annotate.by %in% c('DRIVE', 'Achilles', 'CRISPR')){
    dep = loadData(type=annotate.by, gene=annotate.gene)
    gene_annot <- dep[match(cls, names(dep), nomatch = NA)]
  }
  
  if(all(is.na(gene_annot))){
    warning("No known annotations for " %,% annotate.by %,% ".  Instead annotating by tissue.")
    gene_annot <- as.factor(as.character(gene_annot.by_tissue))
  }
  
  df <- data.frame(
    names=cls, 
    clean_names = sapply(cls, strsplit2, split="_", n=1),
    drive_scores=gene_effect_drive[match(cls, names(gene_effect_drive))], 
    avana_scores=gene_effect_avana[match(cls, names(gene_effect_avana))],
    achilles_scores=gene_effect_achilles[match(cls, names(gene_effect_achilles))],
    combined_scores=gene_effect_combined[match(cls, names(gene_effect_combined))],
    context=gene_annot
  )
  
  return(df)
}

network_graph <- function(gene, cor_mat="DRIVE", r_cutoff=0.15, k=5, save=NULL){
  # Given a correlation matrix of the genome and a gene, this function should create a network around that gene with k neighbors
  # First, validate inputs:
  type=NULL
  # get correlation matrix
  if(is.character(cor_mat)){
    type=cor_mat
    cor_mat <- switch(cor_mat,
                      DRIVE=load.from.taiga(data.name='drive-spearman-4707'),
                      Achilles=load.from.taiga(data.name='achilles-spearman-1a1a'),
                      CRISPR=load.from.taiga(data.name='crispr-spearman-254d'),
                      Combined=load.from.taiga(data.name='combined-pearson-6b02')
    )  
  }
  
  if(!gene %in% colnames(cor_mat)){
    warning(paste0("Gene: ", gene, " not found in correlation matrix."))
    return(NULL)
  }
  
  # Filter to only display k neighbors in network
  if(sum(cor_mat[,colnames(cor_mat)==gene]!=0) <= 2*k){
    warning("k-value is too high or r_cutoff is too low.")
  }
  i.poscor = order(cor_mat[,colnames(cor_mat)==gene], decreasing = T)[-1][1:k]
  i.negcor = order(cor_mat[,colnames(cor_mat)==gene], decreasing = F)[1:k]
  
  # Create node and edge objects for graph display:
  genes = c(colnames(cor_mat)[c(i.poscor, i.negcor)], gene)
  
  # generate iGraph object to create edge list
  library(igraph)
  network=graph_from_adjacency_matrix(cor_mat[genes, genes], weighted=T, mode='undirected', diag = F)
  links <- igraph::as_data_frame(network)

  links <- links[abs(links$weight) >= r_cutoff,]
  # d=as.dist(distances(network, weights = abs(E(network)$weight)))
  # plot(hclust(d))
  
  # Create node list manually:
  library(taigr)
  propeller_list <- readRDS("~/Data/propeller_list.RDS")
  cgc <- read.csv(file = "~/Data/cancer_gene_census_v2.csv")$Gene.Symbol
  stringdb <- load.from.taiga(data.name='stringdb-neighboring-genes-ensp-converted-tohgnc-symbol', data.version=1)
  
  is.propeller <- genes %in% propeller_list
  is.cancer_gene <- genes %in% cgc
  
  if(!any(is.propeller)){
    warning("No propellers were found in any immediate neighbors")
    # if(as.logical(as.integer(readline(prompt = "Display current graph (0) or search neighbors for propeller correlates (1)?")))){
    #   deep_network_graph(cor_mat=type, gene=gene, save=save)
    # }
  }
  
  validated.with.gene <- sapply(genes, function(x){
    df <- stringdb[stringdb$gene1==x & stringdb$gene2==gene,]
    return(nrow(df)>0)
  })
  
  nodes <- data.frame(id=genes, propeller_status = as.numeric(is.propeller)+1, cgc_status = as.numeric(is.cancer_gene)+1, validated=validated.with.gene)
  
  # Generate interactive network graph with visNetwork library
  library(visNetwork)
  #nodes$size <- 20*abs(cor_mat[gene, as.character(nodes$id)])
  nodes$label <- nodes$id
  nodes$title <- sapply(nodes$validated, function(x) paste0("Validated interaction with ", gene, " in StringDB: ", x))
  nodes$shape <- c("dot", "circularImage")[nodes$propeller_status]
  nodes$image <- c("", "https://d30y9cdsu7xlg0.cloudfront.net/png/77747-200.png")[nodes$propeller_status]
  nodes$color.background <- ifelse(nodes$validated, "#D8B365", "white")
  nodes$color.highlight.background <- "yellow"
  nodes$color.border <- c("black", "red")[nodes$cgc_status]
  nodes$color.highlight.border <- nodes$color.border
  nodes$borderWidth <- c(1, 4)[nodes$cgc_status]
  nodes$font.size <- 20
  nodes$font.face <- 'bold'
  nodes$font.color <- 'white'
  
  links$width <- sqrt(abs(links$weight))*10
  links$smooth=F
  links$color <- c('#FC8D59', '#67A9CF')[as.numeric(links$weight>0) + 1]
  links$title=round(links$weight, digits=2)
  links$font.size=7
  
  vn = visNetwork(nodes=nodes, edges=links, main=list(text=paste0(gene, " - ", type), style="font-family:Georgia, Times New Roman, Times, serif;font-size:24px;text-align:center;color:white;"), background = 'black') %>%
    visIgraphLayout()
  
  if(!is.null(save)){
    if(is.character(save)){
      visSave(vn, file=paste0(save,"/",gene,"_",type,"_network.html"))
    }else{
      visSave(vn, file=paste0("~/Downloads/tmp/",gene,"_",type,"_network.html"))
    }
  }
  
  return(vn)
}

#' Create a deep network graph, recursively searching a correlation matrix until a particular query is found
#' 
#' @param query: A particular gene to query as a root node for a network graph.
#' @param type: Which correlation matrix to use (DRIVE, Achilles, CRISPR)
#' @param searchFor: A gene or list of genes to be linked to query.  For our purposes, usually propellers.
#' @param cor_mat: A matrix of correlations.  If null, type will provide matrix.
#' @param k: number of neighbors to search around query node (actually returns 2*k neighbors, positive and negatively correlated).
#' @param recurse_cutoff: Number of recursions from root node.
#' @param r_cutoff: Correlation coefficient to use as cutoff - no links below this cutoff should be displayed.
deep_network_graph <- function(query, type, cor_mat=NULL, searchFor, k=5, recurse_cutoff=2, r_cutoff=0.25){
  # input validation
  if(is.null(cor_mat)){
    cor_mat <- switch(type,
                      DRIVE=readRDS("~/Documents/wdr-data/correlations/drive_cor_spearman.RDS"),
                      Achilles=readRDS("~/Documents/wdr-data/correlations/achilles_cor_spearman.RDS"),
                      CRISPR=readRDS("~/Documents/wdr-data/correlations/crispr_cor_spearman.RDS")
    ) 
  }

  if(!query %in% colnames(cor_mat)){
    warning("Query does not exist in correlation matrix!")
    return(NULL)
  }
  if(!any(searchFor %in% colnames(cor_mat))){
    warning("No searched items exist in correlation matrix!")
    return(NULL)
  }
  
  # Extract correlates, filter, and sort
  correlates <- cor_mat[,query]
  top_k <- c(sort(correlates)[1:k], sort(correlates, decreasing = T)[1:k])
  
  # base case:
  if(any(names(top_k) %in% searchFor)){
    library(igraph)
    cor_mat_filt=cor_mat
    cor_mat_filt[abs(cor_mat) < r_cutoff] <- 0
    network=graph_from_adjacency_matrix(cor_mat_filt[names(top_k), names(top_k)], weighted=T, mode='undirected', diag = F)
    links <- igraph::as_data_frame(network)
    links[abs(links$weight) < r_cutoff,'weight'] <- 0
    links$type <- type
    return(links)
    
  # recursive case:
  }else{
    if(recurse_cutoff > 0){
      library(dplyr)
      library(igraph)
      a <- lapply(names(top_k), function(x){
        n_l <- deep_network_graph(query=x, type=type, cor_mat=cor_mat, searchFor = searchFor, k=k, recurse_cutoff = recurse_cutoff - 1, r_cutoff = r_cutoff)
        return(bind_rows(n_l))
      })
      df <- bind_rows(a)
      network=graph_from_adjacency_matrix(cor_mat[names(top_k), names(top_k)], weighted=T, mode='undirected', diag = F)
      links <- igraph::as_data_frame(network)
      links[abs(links$weight) < r_cutoff,'weight'] <- 0
      links$type <- type
      return(rbind(links, df))
    }else{
      return(NULL)
    }
  }
}

find_propellers_graph <- function(gene, save=NULL, searchFor=NULL, annot=NULL, ...){
  # Define query family - default to propellers
  if(is.null(searchFor)){
    searchFor <- readRDS("~/Data/propeller_list.RDS")
  }
  
  # Define links data frame
  links_DRIVE <- deep_network_graph(query=gene, type="DRIVE", searchFor=searchFor, ...)
  links_CRISPR <- deep_network_graph(query=gene, type="CRISPR", searchFor=searchFor, ...)
  links_Achilles <- deep_network_graph(query=gene, type="Achilles", searchFor=searchFor, ...)
  
  links = rbind(links_DRIVE, links_CRISPR, links_Achilles)
  
  if(!nrow(links)>1){
    stop("No propellers found at current recurse level.  Either increase recurse_cutoff, increase k, or decrease r_cutoff.")
  }
  
  # Calculate distances
  library(igraph)
  g <- graph_from_data_frame(links, directed=F, vertices = NULL)
  dist <- shortest_paths(g, from = gene, weights = NA)
  degrees_of_separation <- sapply(dist$vpath, length)
  
  # Create nodes data frame
  library("taigr")
  genes = unique(c(links$from, links$to))
  
  cgc <- read.csv(file = "~/Documents/wdr-data/cancer_gene_census.csv")$Gene.Symbol
  stringdb <- load.from.taiga(data.name='stringdb-neighboring-genes-ensp-converted-tohgnc-symbol', data.version=1)
  
  is.propeller <- genes %in% searchFor
  is.cancer_gene <- genes %in% cgc
  validated.with.gene <- sapply(genes, function(x){
    df <- stringdb[stringdb$gene1==x & stringdb$gene2==gene,]
    return(nrow(df)>0)
  })
  
  nodes <- data.frame(id=genes, propeller_status = as.numeric(is.propeller)+1, cgc_status = as.numeric(is.cancer_gene)+1, validated=validated.with.gene, degrees_of_separation=degrees_of_separation)
  
  # Generate interactive network graph with visNetwork library
  library(visNetwork)
  #nodes$size <- 20*abs(cor_mat[gene, as.character(nodes$id)])
  nodes$label <- nodes$id
  nodes$title <- sapply(nodes$degrees_of_separation, function(x) paste0(x-1, " degree(s) of separation from ", gene))
  nodes$shape <- c("dot", "circularImage")[nodes$propeller_status]
  nodes$image <- c("", "https://upload.wikimedia.org/wikipedia/commons/thumb/2/2e/1erj_7bladed_beta_propeller.png/220px-1erj_7bladed_beta_propeller.png")[nodes$propeller_status]
  nodes$brokenImage <- c("", "https://d30y9cdsu7xlg0.cloudfront.net/png/77747-200.png")[nodes$propeller_status]
  nodes$color.background <- ifelse(nodes$validated, "#D8B365", "white")
  nodes$color.background <- ifelse(nodes$id %in% annot, "lightgreen", nodes$color.background)
  nodes$color.highlight.background <- "yellow"
  nodes$color.border <- c("black", "red")[nodes$cgc_status]
  nodes$color.highlight.border <- nodes$color.border
  nodes$borderWidth <- c(1, 4)[nodes$cgc_status]
  nodes$font.size <- 20
  nodes$font.face <- 'bold'
  nodes$font.color <- 'white'
  
  links$width <- sqrt(abs(links$weight))*10
  links$length <- abs(links$weight)*100
  links$smooth=F
  links$hidden <- ifelse(links$weight==0, T, F)
  links$color <- c('#FC8D59', '#67A9CF')[as.numeric(links$weight>0) + 1]
  links$title=paste0(links$type, ": ", round(links$weight, digits=2))
  links$font.size=7
  
  vn = visNetwork(nodes=nodes, edges=links, main=list(text=paste0(gene), style="font-family:Georgia, Times New Roman, Times, serif;font-size:24px;text-align:center;color:white;"), background = 'black') %>%
    visIgraphLayout()
  
  if(!is.null(save)){
    visSave(vn, file=paste0(save,"/",gene,"_allData_deep_network.html"), background = 'black')
  }
  
  return(vn)
}

volcano_plot <- function(lm_df, static=F, label=NULL){
  library(ggplot2)
  library(magrittr)
  
  plot_df <- data.frame(
    Gene=lm_df$Gene,
    neg_log_pvalue = -log(lm_df$p.value),
    log_fc = lm_df$EffectSize
  )
  
  plot_df %<>%
    mutate(Significant = -log(lm_df$q.value) > 3 & abs(lm_df$EffectSize) > 0.25)
  
  gg = ggplot(data=plot_df, aes(x=log_fc, y=neg_log_pvalue, fill=Significant, text=Gene)) +
    geom_point()
  
  if(!is.null(label) & length(label)==nrow(plot_df)){
    gg + 
      geom_text(aes(label=label))
    print("Added labels")
  }else{
    warning("No labels!")
  }
  
  if(static){
    return(gg)
  }
  
  plotly::ggplotly(gg)
}


# ANALYSIS FUNCTIONS ------------------------------------------------------



#' Estimate linear-model stats for a matrix of data using limma with empirical Bayes moderated t-stats for p-values
#' Borrowed from CDS team's cdsr package.
#'
#' @param mat: Nxp data matrix with N cell lines and p genes
#' @param vec: N vector of independent variables. Can be two-group labels as factors, bools, or can be numeric
#' @param covars: Optional Nxk matrix of covariates
#' @param weights: Optional N vector of precision weights for each data point
#' @param target_type: Name of the column variable in the data (default 'Gene')
#'
#' @return: data frame of stats
#' @export
#'
#' @examples
#' CRISPR = load.from.taiga(data.name='avana-2-0-1-d98f', 
#' data.version=1, 
#' data.file='ceres_gene_effects',
#' transpose = T)
#' is_panc <- load.from.taiga(data.name = 'ccle-lines-lineages') %>% .[, 'pancreas']
#' ulines <- intersect(rownames(CRISPR), names(is_panc))
#' lim_res <- run_lm_stats_limma(CRISPR[ulines,], is_panc[ulines])
#' @export run_lm_stats_limma
run_lm_stats_limma <- function(mat, vec, covars = NULL, weights = NULL, target_type = 'Gene') {
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  
  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  } else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata,]))
  }
  #there must be more than one unique value of the independent variable
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  #if using covariates add them as additional predictors to the model
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata,, drop = FALSE]
    combined[['pred']] <- pred
    form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  } else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata,])
    } else{
      weights <- weights[udata]
    }
  }
  fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
  fit <- limma::eBayes(fit)
  targ_coef <- grep('pred', colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf) 
  
  if (colnames(results)[1] == 'ID') {
    colnames(results)[1] <- target_type
  } else {
    results %<>% rownames_to_column(var = target_type)  
  } 
  results$min_samples <- min_samples[results[[target_type]]]            
  
  two_to_one_sided <- function(two_sided_p, stat, test_dir) {
    #helper function for converting two-sided p-values to one-sided p-values
    one_sided_p <- two_sided_p / 2
    if (test_dir == 'right') {
      one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
    } else {
      one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
    }
    return(one_sided_p)
  }
  results %<>% set_colnames(revalue(colnames(.), c('logFC' = 'EffectSize', 'AveExpr' = 'Avg', 't' = 't_stat', 'B' = 'log_odds',
                                                   'P.Value' = 'p.value', 'adj.P.Val' = 'q.value', 'min_samples' = 'min_samples'))) %>% na.omit()
  results %<>% dplyr::mutate(p.left = two_to_one_sided(p.value, EffectSize, 'left'),
                             p.right = two_to_one_sided(p.value, EffectSize, 'right'),
                             q.left = two_to_one_sided(q.value, EffectSize, 'left'),
                             q.right = two_to_one_sided(q.value, EffectSize,'right'))
  return(results)
}

find_correlations_from_dependency <- function(gene=NULL, k=10, r_cutoff=0.1, from=c("DRIVE", "Achilles", "CRISPR", "Combined"), to=c("mutation", "copy number", "expression", 'dependency', 'methylation', 'chromatin', 'metabolome', 'gsea_hallmarks', 'bioplex'), annotate=T, manual_dependency=NULL, method='spearman', updateProgress=NULL){
  return_df <- data.frame()
  dependency_data <- NULL
  if(!is.null(manual_dependency))
    dependency_data <- list("Manual"=manual_dependency)
  if(!is.null(gene))
    dependency_data <- sapply(from, loadData, gene=gene)
  
  stopifnot(!is.null(dependency_data))
  library("taigr")
  library('dplyr')
  
  progress <- function(updateProgress, status, value){
    if(!is.function(updateProgress)){
      return()
    }else{
      updateProgress(detail = paste0("Completing ", status), value=value)
    }
  }
  
  if('dependency' %in% to){
    progress(updateProgress, status="dependency calculations", value=1/8)
    dependency_else <- sapply(from, loadData)
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(dependency_else[[x]]))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], dependency_else[[x]][cls,], use='pairwise.complete.obs', method=method)
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        from = x,
        to = x
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
  }
  
  if('mutation' %in% to){
    progress(updateProgress, status="mutation calculations", value=2/8)
    
    pooled_mut <- load.from.taiga(data.name='pooled-mutation-6481', data.version=1)
    colnames(pooled_mut) <- sapply(colnames(pooled_mut), strsplit2, split=" ", n=1)
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(pooled_mut))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], pooled_mut[cls,], use='complete', method=method)
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        from = x,
        to = "Mutation"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
  }
  
  if('copy number' %in% to){
    progress(updateProgress, status="copynumber calculations", value=3/8)
    
    gene.CN.SNP.priority <- load.from.taiga(data.name='gene-level-cn-87aa', data.version=3, data.file='gene_CN_SNP_priority')
    rownames(gene.CN.SNP.priority) <- sapply(rownames(gene.CN.SNP.priority), strsplit2, split="snp_", n=2)
    colnames(gene.CN.SNP.priority) <- sapply(colnames(gene.CN.SNP.priority), strsplit2, split=" ", n=1)
    
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(gene.CN.SNP.priority))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], gene.CN.SNP.priority[cls,], use='pairwise', method=method)
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        from = x,
        to = "Copy Number"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
  }
  
  if('expression' %in% to){
    progress(updateProgress, status="expression calculations", value=4/8)
    
    RNAseq <- loadData(type='RNAseq')
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(RNAseq))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], RNAseq[cls,], use='complete', method=method)
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        from = x,
        to = "Expression"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
  }
  
  if('methylation' %in% to){
    progress(updateProgress, status="methylation calculations", value=5/8)
    
    RBBS <- loadData(type='methylation')
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(RBBS))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], RBBS[cls,], use='pairwise.complete.obs', method=method)
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        from = x,
        to = "Methylation"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
  }
  
  if('chromatin' %in% to){
    progress(updateProgress, status="chromatin calculations", value=6/8)
    
    chromatin <- loadData("chromatin")
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(chromatin))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], chromatin[cls,], use='complete', method=method)
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        from = x,
        to = "Chromatin"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
  }
  
  if('metabolome' %in% to){
    progress(updateProgress, status="metabolomics calculations", value=7/8)
    
    metabolome <- loadData("metabolome")
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(metabolome))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], metabolome[cls,], use='complete', method=method)
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        from = x,
        to = "Metabolome"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
  }
  
  if('gsea_hallmarks' %in% to){
    progress(updateProgress, status="dependency calculations", value=8/8)
    
    gsea_hallmarks <- t(load.from.taiga(data.name='ssgsea-enrichment-scores-for-msigdb-h-using-ccle-rnaseq-expression', data.version=1))
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(gsea_hallmarks))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], gsea_hallmarks[cls,], use='complete', method=method)
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        from = x,
        to = "GSEA Hallmark"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
  }
  
  
  if("bioplex" %in% to){
    progress(updateProgress, status="Bioplex integrations", value=8/8)
    
    bioplex <- load.from.taiga(data.name='bioplex-ppi-e127', data.version=1)
    df_1 <- bioplex[bioplex$SymbolA==gene,]
    df_2 <- bioplex[bioplex$SymbolB==gene,]
    
    p <- data.frame(
      genes = c(df_1$SymbolB, df_2$SymbolA),
      scores = c(df_1$`p(Interaction)`, df_2$`p(Interaction)`),
      from = "BioPlex",
      to = "BioPlex"
    )
    
    return_df <- rbind(return_df, p)
  }
  ### Additional filtering and annotation
  # Filter by correlation coefficient
  progress(updateProgress, status="final annotations", value=8/8)
  
  return_df <- return_df[abs(return_df$scores)>r_cutoff,]
  
  # Duplicate correlations are marked TRUE
  ii = which(return_df$genes %in% return_df$genes[duplicated(return_df$genes)])
  return_df$starred <- FALSE
  return_df$starred[ii] <- TRUE
  
  # CGCs are marked as TRUE
  cgc <- as.character(read.csv('~/Data/cancer_gene_census_v2.csv')$Gene.Symbol)
  return_df$isCGC <- ifelse(return_df$genes %in% cgc, TRUE, FALSE)
  
  # Propellers are marked TRUE
  propeller_list <- readRDS('~/Data/propeller_list.RDS')
  return_df$isPropeller <- ifelse(return_df$genes %in% propeller_list, TRUE, FALSE)
  
  # Additional Annotation
  return_df$annot <- annotate_genelist(return_df$genes, type = 'position')
  
  return(return_df)
}

find_lm_from_dependency <- function(gene, k=10, r_cutoff=0.1, from=c("DRIVE", "Achilles", "CRISPR", "Combined"), to=c("mutation", "copy number", "expression", 'dependency', 'methylation', 'chromatin', 'metabolome', 'gsea_hallmarks'), annotate=T, updateProgress=NULL){
  return_df <- data.frame()
  dependency_data <- sapply(from, loadData, gene=gene)
  library("taigr")
  library('dplyr')
  
  progress <- function(updateProgress, status, value){
    if(!is.function(updateProgress)){
      return()
    }else{
      updateProgress(detail = paste0("Completing ", status), value=value)
    }
  }
  
  if('dependency' %in% to){
    dependency_else <- sapply(from, loadData)
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(dependency_else[[x]]))
      if(length(cls)==0)
        return(NULL)
      lm_df = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=dependency_else[[x]][cls,])
      ii <- order(lm_df$log_odds, decreasing = T)[1:k]
      data.frame(
        genes = lm_df$Gene[ii],
        scores = lm_df$EffectSize[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = x
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    progress(updateProgress, status="dependency calculations", value=1/8)
  }
  
  if('Achilles' %in% to){
    Achilles <- loadData(type='Achilles')
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(Achilles))
      if(length(cls)==0)
        return(NULL)
      lm_df = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=Achilles[cls,])
      ii <- order(lm_df$log_odds, decreasing = T)[1:k]
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = "Achilles"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
  }
  
  if('CRISPR' %in% to){
    CRISPR <- loadData(type='CRISPR')
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(CRISPR))
      if(length(cls)==0)
        return(NULL)
      cc = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=CRISPR[cls,])
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = "CRISPR"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    progress(updateProgress, status="dependency calculations", value=2/8)
  }
  
  if('mutation' %in% to){
    pooled_mut <- load.from.taiga(data.name='pooled-mutation-6481', data.version=1)
    colnames(pooled_mut) <- sapply(colnames(pooled_mut), strsplit2, split=" ", n=1)
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(pooled_mut))
      if(length(cls)==0)
        return(NULL)
      lm_df = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=pooled_mut[cls,])
      ii <- order(lm_df$log_odds, decreasing = T)[1:k]
      data.frame(
        genes = lm_df$Gene[ii],
        scores = lm_df$EffectSize[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = "Mutation"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    progress(updateProgress, status="mutation calculations", value=1/8)
  }
  
  if('copy number' %in% to){
    gene.CN.SNP.priority <- loadData(type="cn")
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(gene.CN.SNP.priority))
      if(length(cls)==0)
        return(NULL)
      lm_df = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=gene.CN.SNP.priority[cls,])
      ii <- order(lm_df$log_odds, decreasing = T)[1:k]
      data.frame(
        genes = lm_df$Gene[ii],
        scores = lm_df$EffectSize[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = "Copy Number"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    progress(updateProgress, status="copynumber calculations", value=3/8)
  }
  
  if('expression' %in% to){
    RNAseq <- loadData(type='RNAseq')
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(RNAseq))
      if(length(cls)==0)
        return(NULL)
      lm_df = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=RNAseq[cls,])
      ii <- order(lm_df$log_odds, decreasing = T)[1:k]
      data.frame(
        genes = lm_df$Gene[ii],
        scores = lm_df$EffectSize[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = "Expression"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    progress(updateProgress, status="expression calculations", value=4/8)
  }
  
  if('methylation' %in% to){
    RRBS <- loadData(type='methylation')
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(RRBS))
      if(length(cls)==0)
        return(NULL)
      lm_df = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=RRBS[cls,])
      ii <- order(lm_df$log_odds, decreasing = T)[1:k]
      data.frame(
        genes = lm_df$Gene[ii],
        scores = lm_df$EffectSize[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = "Methylation"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    progress(updateProgress, status="methylation calculations", value=4.5/8)
  }
  
  if('chromatin' %in% to){
    chromatin <- loadData("chromatin")
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(chromatin))
      if(length(cls)==0)
        return(NULL)
      lm_df = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=chromatin[cls,])
      ii <- order(lm_df$log_odds, decreasing = T)[1:k]
      data.frame(
        genes = lm_df$Gene[ii],
        scores = lm_df$EffectSize[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = "Chromatin"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    progress(updateProgress, status="chromatin calculations", value=5/8)
  }
  
  if('metabolome' %in% to){
    metabolome <- loadData("metabolome")
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(metabolome))
      if(length(cls)==0)
        return(NULL)
      lm_df = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=metabolome[cls,])
      ii <- order(lm_df$log_odds, decreasing = T)[1:k]
      data.frame(
        genes = lm_df$Gene[ii],
        scores = lm_df$EffectSize[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = "Metabolome"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    progress(updateProgress, status="metabolomics calculations", value=6/8)
  }
  
  if('gsea_hallmarks' %in% to){
    gsea_hallmarks <- t(load.from.taiga(data.name='ssgsea-enrichment-scores-for-msigdb-h-using-ccle-rnaseq-expression', data.version=1))
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(gsea_hallmarks))
      if(length(cls)==0)
        return(NULL)
      lm_df = run_lm_stats_limma(vec=dependency_data[[x]][cls], mat=gsea_hallmarks[cls,])
      ii <- order(lm_df$log_odds, decreasing = T)[1:k]
      data.frame(
        genes = lm_df$Gene[ii],
        scores = lm_df$EffectSize[ii],
        p.value = lm_df$p.value[ii],
        from = x,
        to = "GSEA Hallmark"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    progress(updateProgress, status="GSEA hallmarks calculations", value=7/8)
  }
  
  ### Additional filtering and annotation
  # Filter by correlation coefficient
  return_df <- return_df[abs(return_df$scores)>r_cutoff,]
  
  # Duplicate correlations are marked TRUE
  ii = which(return_df$genes %in% return_df$genes[duplicated(return_df$genes)])
  return_df$starred <- FALSE
  return_df$starred[ii] <- TRUE
  
  # CGCs are marked as TRUE
  cgc <- as.character(read.csv('~/Data/cancer_gene_census_v2.csv')$Gene.Symbol)
  return_df$isCGC <- ifelse(return_df$genes %in% cgc, TRUE, FALSE)
  
  # Propellers are marked TRUE
  propeller_list <- readRDS("~/Data/propeller_list.RDS")
  return_df$isPropeller <- ifelse(return_df$genes %in% propeller_list, TRUE, FALSE)
  
  # Additional Annotation
  progress(updateProgress, status="final calculations", value=8/8)
  return_df$annot <- annotate_genelist(return_df$genes, type = 'position')
  
  return(return_df)
}

annotate_genelist <- function(genes, type="dependency"){
  if(type %in% c("dependency", 'expression')){
    library(topGO)
    # library(biomaRt)
    # ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    # ensembl_genes_df <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 
    #                                        'chromosome_name', 'start_position', 'end_position', 'band',
    #                                        'entrezgene'), 
    #                           mart = ensembl)
    # annFUN.org(whichOnto = "BP", mapping = 'org.Hs.eg.db', ID = 'symbol')
    
    geneList <- factor(as.integer(genelist %in% genes))
    names(geneList) <- genelist
    GOdata <- new("topGOdata", ontology="BP", allGenes=geneList, annot = annFUN.org, mapping='org.Hs.eg.db', ID='symbol')
    resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    GenTable(GOdata, classic=resultsFis)
  }
  if(type %in% c("copy number", "position", "chromosome")){
    Position <- loadData("Position")
    ii = match(genes, Position$gene)
    
    bands <- sapply(ii, function(x){
      if(is.na(Position$chr[x])){
        return(NA)
      }
      if(is.na(Position$band[x])){
        return(Position$chr[x])
      }
      return(paste0(Position$chr[x], Position$band[x]))
    })
    
    return(bands)
  }
}

# CONSTANTS ---------------------------------------------------------------
# Accessible by loadConstants function
# CRISPR <- loadData("CRISPR")
# COMBINED <- loadData("Combined")
# ACHILLES <- loadData("Achilles")
# DRIVE <- loadData("DRIVE")
# 
# ALLGENES <- union_recursive(list(colnames(CRISPR), colnames(COMBINED), colnames(ACHILLES), colnames(DRIVE)))
# ALLCELLS <- union_recursive(list(rownames(CRISPR), rownames(COMBINED), rownames(ACHILLES), rownames(DRIVE)))
# COMMONGENES <- intersect_recursive(list(colnames(CRISPR), colnames(COMBINED), colnames(ACHILLES), colnames(DRIVE)))
# COMMONCELLS <- intersect_recursive(list(rownames(CRISPR), rownames(COMBINED), rownames(ACHILLES), rownames(DRIVE)))
# 
# save(CRISPR, COMBINED, ACHILLES, DRIVE, ALLGENES, ALLCELLS, COMMONGENES, COMMONCELLS, file = '~/constants.RData')
