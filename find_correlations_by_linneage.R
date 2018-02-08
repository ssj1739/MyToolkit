find_correlations_by_linneage <- function(gene=NULL, k=10, r_cutoff=0.1, from=c("DRIVE", "Achilles", "CRISPR", "Combined"), to=c("mutation", "copy number", "expression", 'dependency', 'methylation', 'chromatin', 'metabolome', 'gsea_hallmarks', 'bioplex'), annotate=T, linneage=NULL, method='spearman', updateProgress=NULL){
  return_df <- data.frame()
  dependency_data <- NULL
  if(!is.null(linneage)){
    dependency_data <- sapply(from, loadData, gene=gene)
    dependency_data = lapply(dependency_data, function(x){
      ii = grep(linneage, names(x))
      return(x[ii])
    })
  }
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
    if(!is.null(linneage)){
      dependency_else <- sapply(dependency_else, function(x){
        ii = grep(linneage, rownames(x))
        return(x[ii,])
      })
    }
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
    
    pooled_mut <- loadData("mutation")
    if(!is.null(linneage)){
      ii = grep(linneage, rownames(pooled_mut))
      pooled_mut <- pooled_mut[ii,]
    }
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
    
    copynumber <- loadData("copynumber")
    if(!is.null(linneage)){
      ii = grep(linneage, rownames(copynumber))
      copynumber <- copynumber[ii,]
    }
    
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(copynumber))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], copynumber[cls,], use='pairwise', method=method)
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
    if(!is.null(linneage)){
      ii = grep(linneage, rownames(RNAseq))
      RNAseq <- RNAseq[ii,]
    }
    
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
    
    RRBS <- loadData(type='methylation')
    if(!is.null(linneage)){
      ii = grep(linneage, rownames(RRBS))
      RRBS <- RRBS[ii,]
    }
    
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(names(dependency_data[[x]]), rownames(RRBS))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls], RRBS[cls,], use='pairwise.complete.obs', method=method)
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
    if(!is.null(linneage)){
      ii = grep(linneage, rownames(chromatin))
      chromatin <- chromatin[ii,]
    }
    
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
    if(!is.null(linneage)){
      ii = grep(linneage, rownames(metabolome))
      metabolome <- metabolome[ii,]
    }
    
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
    
    gsea_hallmarks <- loadData("GSEA")
    if(!is.null(linneage)){
      ii = grep(linneage, rownames(gsea_hallmarks))
      gsea_hallmarks <- gsea_hallmarks[ii,]
    }
    
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
    if(nrow(df_1) + nrow(df_2) == 0){ # gene not bait or prey node
      return_df <- return_df
    }else{
      p <- data.frame(
        genes = c(df_1$SymbolB, df_2$SymbolA),
        scores = c(df_1$`p(Interaction)`, df_2$`p(Interaction)`),
        from = "BioPlex",
        to = "BioPlex"
      )
      return_df <- rbind(return_df, p)
    }
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
  cgc <- as.character(read.csv('data/cancer_gene_census_v2.csv')$Gene.Symbol)
  return_df$isCGC <- ifelse(return_df$genes %in% cgc, TRUE, FALSE)
  
  # Propellers are marked TRUE
  propeller_list <- readRDS('data/propeller_list.RDS')
  return_df$isPropeller <- ifelse(return_df$genes %in% propeller_list, TRUE, FALSE)
  
  # Additional Annotation
  return_df$annot <- annotate_genelist(return_df$genes, type = 'position')
  
  return(return_df)
}