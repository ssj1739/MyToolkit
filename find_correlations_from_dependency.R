#' Correlating Dependencies
#' @name find_correlations_from_dependency
#' @author Sidharth Jain
#' @description Calculates all correlations between the dependencies for a given gene to all other 
#' features in other datasets, including copy number, mutation, expression, etc.
#' @param gene - character of HUGO gene symbol
#' @param k - integer of nearest correlates to return (will return 2*k, positive and negative)
#' @param r_cutoff - cutoff for r.squared value
#' @param from - character vector of dependency datasets to extract from (see details)
#' @param to - character vector of datasets to correlate to correlate to (see details)
#' @param annotate - logical, annotate by chromosomal position?
#' @param manual_dependency - optional numeric vector of dependency
#' @param method - character to be passed to "cor", usually spearman or pearson (or kendall)
#' @param updateProgress - only include if using in Shiny app context to update progress bar
#' @export
find_correlations_from_dependency <- function(gene=NULL, k=10, r_cutoff=0.1, from=c("DRIVE", "Achilles", "CRISPR", "Combined"), to=c("mutation", "copy number", "expression", 'dependency', 'methylation', 'chromatin', 'metabolome', 'gsea_hallmarks', 'bioplex', 'protein'), annotate=T, manual_dependency=NULL, method='spearman', updateProgress=NULL){
  require('taigr')
  library("taigr")
  library('dplyr')
  
  return_df <- data.frame()
  dependency_data <- NULL
  if(!is.null(manual_dependency))
    dependency_data <- list("Manual"=manual_dependency)
  if(!is.null(gene))
    dependency_data <- lapply(from, loadData, gene=gene)
  names(dependency_data) <- from
  
  stopifnot(!is.null(dependency_data))

  
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
      cls <- intersect(rownames(dependency_data[[x]]), rownames(dependency_else[[x]]))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls,], dependency_else[[x]][cls,], use='pairwise.complete.obs', method=method)
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
    colnames(pooled_mut) <- sapply(colnames(pooled_mut), strsplit2, split=" ", n=1)
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(rownames(dependency_data[[x]]), rownames(pooled_mut))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls,], pooled_mut[cls,], use='complete', method=method)
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
    
    gene.CN.SNP.priority <- loadData('copynumber')
    
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(rownames(dependency_data[[x]]), rownames(gene.CN.SNP.priority))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls,], gene.CN.SNP.priority[cls,], use='pairwise', method=method)
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
      cls <- intersect(rownames(dependency_data[[x]]), rownames(RNAseq))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls,], RNAseq[cls,], use='complete', method=method)
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
      cls <- intersect(rownames(dependency_data[[x]]), rownames(RBBS))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls,], RBBS[cls,], use='pairwise.complete.obs', method=method)
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
      cls <- intersect(rownames(dependency_data[[x]]), rownames(chromatin))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls,], chromatin[cls,], use='complete', method=method)
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
      cls <- intersect(rownames(dependency_data[[x]]), rownames(metabolome))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls,], metabolome[cls,], use='complete', method=method)
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
    progress(updateProgress, status="gsea calculations", value=8/8)
    
    gsea_hallmarks <- t(load.from.taiga(data.name='ssgsea-enrichment-scores-for-msigdb-h-using-ccle-rnaseq-expression', data.version=1))
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(rownames(dependency_data[[x]]), rownames(gsea_hallmarks))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls,], gsea_hallmarks[cls,], use='complete', method=method)
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
  
  if('protein' %in% to){
    progress(updateProgress, status="protein calculations", value=8/8)
    
    protein <- load.from.taiga(data.name='ms-protein-a729')
    p <- lapply(names(dependency_data), function(x){
      cls <- intersect(rownames(dependency_data[[x]]), rownames(protein))
      if(length(cls)==0)
        return(NULL)
      cc = cor(dependency_data[[x]][cls,], protein[cls,], use='pairwise.complete', method=method)
      ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
      data.frame(
        genes = colnames(cc)[ii],
        scores = cc[ii],
        from = x,
        to = "MS Protein"
      )
    })
    return_df <- rbind(return_df, bind_rows(p))
    
    # RPPA <- loadData("RPPA")
    # p <- lapply(names(dependency_data), function(x){
    #   cls <- intersect(rownames(dependency_data[[x]]), rownames(RPPA))
    #   if(length(cls)==0)
    #     return(NULL)
    #   cc = cor(dependency_data[[x]][cls,], RPPA[cls,], use='pairwise.complete', method=method)
    #   ii <- c(order(cc, decreasing = T)[1:k], order(cc, decreasing=F)[k:1])
    #   data.frame(
    #     genes = colnames(cc)[ii],
    #     scores = cc[ii],
    #     from = x,
    #     to = "RPPA Protein"
    #   )
    # })
    # return_df <- rbind(return_df, bind_rows(p))
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