#' Convert gene-set collection object to sparse matrix
#'
#' @param gsc gene set collection object
#' @param gene_list vector of gene names, specifying the set of genes included in analysis
#'
#' @return: sparse binary matrix of gene set memberships
#' @export
#'
#' @examples
gsc_to_sparseMat <- function(gsc, gene_list) {
    gsc_M <- Matrix::Matrix(0, nrow = length(gene_list), ncol = length(gsc), sparse = T)
    rownames(gsc_M) <- gene_list
    colnames(gsc_M) <- names(gsc)
    for (ii in seq_along(gsc)) {
        cur_in <- gene_list %in% GSEABase::geneIds(gsc[[ii]])
        gsc_M[cur_in, ii] <- rep(1, sum(cur_in))
    }
    gsc_M
}

#' Run npGSEA analysis, testing for gene set enrichment of association with a vector across a collection of gene sets
#'
#' @param gsc: Gene set collection object
#' @param X: Matrix of gene-level data. Rows are cell line, columns are genes
#' @param y: Vector of data to compare with dep_mat. Can be categorical or continuous
#' @param covars: optional matrix of covariates to control for. This matrix is linearly projected out of both X and y before further analysis
#' @param approx: type of approximation to use [default 'norm', other options 'beta' and 'chiSq'. See npGSEA docs]
#' @param min_set_size: minimum size of gene sets to consider [default 3]
#' @param max_set_size: Max size of gene sets to consider [default 300]
#' @param normalize_X: whether to normalize scale of data in each column of the X matrix to have unit variance [default FALSE]
#' @param n_top_genes: Optional number of top contributing genes from each set to return [default 0]
#' @param top_dir: If using n_top_genes, which 'direction' of test statistic should genes be ranked by [default 'left']
#' @param nproc: How many processes to use (default 0, no multithreading). Only helps when using covariates
#'
#' @return: data frame of results. Below assumes you're using the normal approximation. See docs if you're using a different model
#' \itemize{
#' \item gene_set: name of gene set
#' \item enrich_stat: sum of gene-level test statistics over genes in set
#' \item norm_enrich_stat: mean of gene-level test stats over genes in set
#' \item set_size: number of genes in set
#' \item p_left/p_right/p_twoSided: P-values for testing whether gene stats are decreased/increased/either for genes in set
#' \item q_left/q_right/q_twoSided: FDR_adjusted p-values
#' }
#'
#'
#' @export
#'
#' @examples
run_npGSEA <- function(gsc, X, y, covars = NULL, approx = 'norm',
                       min_set_size = 3, max_set_size = 300,
                       normalize_X = FALSE, n_top_genes = 0, top_dir = 'left',
                       nproc = 0) {
    library(npGSEAm)
    is_used <- !is.na(y)
    if (!is.null(covars)) {
       if (is.vector(covars)) {
         covars <- data.frame(matrix(covars, ncol = 1))
         colnames(covars) <- 'cov'
       }
        is_used[rowSums(is.na(covars)) > 0] <- FALSE
        covars %<>% filter(is_used)
        covars <- model.matrix(as.formula(paste0('~0 + ', paste(colnames(covars), collapse = '+'))), covars)
    }
    X <- X[is_used,]
    Y <- y[is_used]

    #project covar matrix out of each column of X, and from Y before running analysis
    if (!is.null(covars)) {
      if (nproc > 0) {
        library(plyr)
        library(doMC)
        doMC::registerDoMC(cores=nproc)
        X <- plyr::aaply(X, 2, function(col) {
          ncol <- rep(NA, length(col))
          ncol[!is.na(col)] <- resid(lm(col ~ covars))
          return(ncol)
        }, .parallel=TRUE) %>% t()
      } else {
        X <- apply(X, 2, function(col) {
          ncol <- rep(NA, length(col))
          ncol[!is.na(col)] <- resid(lm(col ~ covars))
          return(ncol)
        })
      }
      Y <- resid(lm(Y ~ covars))
    }
    used_genes <- colnames(X)
    gsc_M <- gsc_to_sparseMat(gsc, used_genes)
    n_used_genes <- Matrix::colSums(gsc_M)
    use_gscs <- n_used_genes >= min_set_size & n_used_genes <= max_set_size
    cur_gsc <- gsc[use_gscs]
    gsc_M <- gsc_M[, use_gscs]

    res <- npGSEAm::npGSEA(t(X), Y, cur_gsc, covars = NULL,
                           approx = approx, scaleXY=TRUE, uniVarX=normalize_X)

    results <- data.frame(gene_set = names(cur_gsc),
                          enrich_stat = unlist(stat(res)),
                          p_left = unlist(pLeft(res)),
                          p_twoSided = unlist(pTwoSided(res)),
                          p_right = unlist(pRight(res)),
                          set_size = n_used_genes[use_gscs]) %>%
        mutate(norm_enrich_stat = enrich_stat / set_size,
               q_left = p.adjust(p_left, method = 'BH'),
               q_right = p.adjust(p_right, method = 'BH'),
               q_twoSided = p.adjust(p_twoSided, method = 'BH'))

    if (n_top_genes > 0) {
        print('computing top contributing genes')
        all_betas <- betaHats(res)
        top_genes <- lapply(all_betas, function(x) {
            x %<>%
                as.data.frame() %>%
                rownames_to_column(var = 'Gene')
            if (top_dir == 'left') {
                x %<>% dplyr::arrange(V1)
            } else if (top_dir == 'right') {
                x %<>% dplyr::arrange(desc(V1))
            } else {stop('invalid top_dir')}
            x %>% head(n_top_genes) %>% .[['Gene']]
        })
        results$top_genes = I(top_genes)
    }
    results
}


#' Helper function for computing gene level stats
#'
#' @param X Samples by genes data matrix
#' @param y Phenotype vector
#' @param stat_type Type of stat c(t_stat, cor, mod_t_stat, pval, mod_pval, regr_coef)]
#' @param stat_trans Type of transformation c(square, abs, rank)
#'
#' @return
#'
#' @examples
get_gene_stats <- function(X, y, stat_type = 't_stat', stat_trans = 'none') {
  stopifnot(stat_type %in% c('t_stat', 'cor', 'pval', 'mod_pval', 'mod_t_stat', 'regr_coef'))
  dir <- NULL
  if (stat_type == 'mod_t_stat') {
    stats <- run_lm_stats_limma(mat = X, vec = y) %>% .[['t_stat']]
    dir <- sign(stats)
  } else if (stat_type == 'regr_coef') {
    stats <- run_lm_stats_limma(mat = X, vec = y) %>% .[['EffectSize']]
    dir <- sign(stats)
  } else if (stat_type == 'cor') {
    stats <- corr_test(X, y, adjust = 'none') %>% .[['r']]
    dir <- sign(stats)
  } else if (stat_type == 't_stat') {
    stats <- corr_test(X, y, adjust = 'none') %>% .[['t']][,1]
    dir <- sign(stats)
  } else if (stat_type == 'pval') {
    stats_mat <- corr_test(X, y, adjust = 'none')
    stats <- stats_mat$p[, 1]
    dir <- sign(stats_mat$r)
  } else if (stat_type == 'mod_pval') {
    stats_mat <- run_lm_stats_limma(mat = X, vec = y)
    stats <- stats_mat$t_stat
    dir <- sign(stats_mat$p.value)
  }
  stopifnot(stat_trans %in% c('none', 'square', 'abs', 'rank'))
  if (stat_trans == 'square') {
    stats <- stats^2
  } else if (stat_trans == 'abs') {
    stats <- abs(stats)
  } else if (stat_trans == 'rank') {
    stats <- order(stats, decreasing = TRUE)
  }
  names(stats) <- colnames(X)
  names(dir) <- colnames(X)
  used <- which(!is.na(stats))
  stats <- stats[used]
  dir <- dir[used]
  return(list(stat = stats, dir = dir))
}

#' Run GSEA using fgsea package, either with gene or label permutation null
#'
#' @param gsc gene-set collection object from GSEAbase
#' @param X matrix of data with rows as cell lines and columns as genes (same gene names as in gsc object)
#' @param y vector of phenotypes.
#' @param perm_type type of permutations to generate null ['label', or 'gene']
#' @param nperm number of permutations
#' @param min_set_size minimum gene set size
#' @param max_set_size max gene set size
#' @param stat_type if using gene-permutation, stat to use. Options ['t_stat', 'mod_t_stat', 'cor', 'regr_coef', 'pval', 'mod_pval']
#' @param stat_trans if using gene-permutation, stat transformation to apply Options ['none', 'square', 'abs', 'rank']
#' @param gseaParam: GSEA parameter value (power)
#' @param nproc: number of processes to run in parallel
#' @param gene_stat: vector of pre-computed gene-level stats. Only works for 'gene-permutation' testing. IF provided X and y will be ignored
#'
#' @return
#' @export
#'
#' @examples
run_fGSEA <- function(gsc, X = NULL, y = NULL, perm_type = 'label', nperm = 10000,
                       min_set_size = 3, max_set_size = 300,
                      stat_type = 't_stat', stat_trans = 'none',
                      gseaParam = 1, nproc = 0, gene_stat = NULL) {
  if (any(grepl('package:piano', search()))) {
    detach("package:piano", unload=TRUE)
  }
  library(fgsea)
  stopifnot(perm_type %in% c('label', 'gene'))
  if (perm_type == 'label') {
    stopifnot(is.matrix(X) & !is.null(y))
    print(sprintf('Running sample-permutation testing with %d perms', nperm))
    used_samples <- which(!is.na(y))
    used_genes <- which(apply(X[used_samples,], 2, var, na.rm=T) > 0)
    fgseaRes <- fgsea::fgseaL(pathways = geneIds(gsc),
                      mat = t(X[used_samples, used_genes]),
                      labels = y[used_samples],
                      minSize=min_set_size,
                      maxSize=max_set_size,
                      nperm=nperm,
                      gseaParam = gseaParam,
                      nproc = nproc)
  } else if (perm_type == 'gene') {
    print(sprintf('Running gene-permutation testing with %d perms', nperm))
    if (is.null(gene_stat)) {
      gene_stat <- get_gene_stats(X, y, stat_type = stat_type, stat_trans = stat_trans)$stat
    }
    fgseaRes <- fgsea::fgsea(pathways = geneIds(gsc),
                      stats = gene_stat,
                      minSize=min_set_size,
                      maxSize=max_set_size,
                      nperm=nperm,
                      gseaParam = gseaParam,
                      nproc = nproc)
  }
  return(fgseaRes)
}



#' Run gene set analysis using piano package
#'
#' @param gsc gene set collection as gene-set collection object
#' @param X matrix of data with rows as cell lines and columns as genes
#' @param y phenotype vector
#' @param gene_stat Optional: pre-computed gene-level stats. If given, X and y will be ignored
#' @param stat_dir Optional: pre-computed gene-stat directions
#' @param method GSA method to use. Can be ['fisher', 'stouffer', 'reporter', 'tailStrength', 'mean', 'median', 'page']
#' @param stat_type Type of statistic to use for testing enrichment. Can be ['t_stat', 'mod_t_stat', 'cor', 'cor_rank', 'pval', 'mod_pval']. For methods 'fisher', 'stouffer', 'reporter', and 'tailStrength', must be p
#' @param stat_trans Type of transformation applied to gene stat Can be ['none', 'square', 'rank', 'abs'].
#' @param nperm Number of permutations to run (gene-permutation)
#' @param nproc Number of threads to use (requires snowfall library)
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
run_piano <- function(gsc, X = NULL, y = NULL, gene_stat = NULL, stat_dir = NULL, method, stat_type, stat_trans = 'none', nperm = 10000, nproc = 1, verbose = TRUE) {
  library(piano)
  library(GSEABase)

  if (nproc > 1) {
    library(snowfall)
  }
  stopifnot(method %in% c('mean', 'median', 'page', 'fisher', 'stouffer', 'reporter', 'tailStrength'))
  if (method %in% c('fisher', 'stouffer', 'reporter', 'tailStrength')) {
    stopifnot(stat_type %in% c('pval', 'mod_pval'))
  }
  if (method %in% c('page', 'fisher', 'stouffer', 'reporter')) {
    signifMethod <- 'nullDist'
  } else {
    signifMethod <- 'geneSampling'
  }
  gsc_obj <- GSEABase::geneIds(gsc)
  gs_names <- do.call(c, sapply(names(gsc_obj), function(name) {
    rep(name, length(gsc_obj[[name]]))
  }))
  gsc <- data.frame(gene = unlist(gsc_obj), gene_set = gs_names)

  if (is.null(gene_stat)) {
    stats <- get_gene_stats(X, y, stat_type = stat_type, stat_trans = stat_trans)
    gene_stat <- stats$stat
    stat_dir <- stats$dir
  }
  res <- runGSA(gene_stat,
                gsc = loadGSC(gsc),
                geneSetStat = method,
                directions = stat_dir,
                signifMethod = signifMethod,
                nPerm = nperm,
                ncpus = nproc,
                verbose = verbose)
  GSAsummaryTable(res)
}


#' Run set overlap test using piano 'runGSAhyper' function
#'
#' @param hit_genes list of genes scoring as hits
#' @param all_genes list of all genes tested
#' @param gsc gene set collection as gene-set collection object
#' @param gsSizeLim 2x1 vector of size limits on the gene sets (defaults to c(1,Inf))
#' @param adjMethod Method for p-valuea adjustment (from available in p.adjust, default = 'BH')
#'
#' @return Table of stats for all gene sets in collection
#' @export
#'
#' @examples
run_GSAhyper <- function(hit_genes, all_genes, gsc, gsSizeLim = c(1, Inf), adjMethod = 'BH') {
  library(piano)
  library(GSEABase)
  library(tibble)
  gsc_obj <- GSEABase::geneIds(gsc)
  gs_names <- do.call(c, sapply(names(gsc_obj), function(name) {
    rep(name, length(gsc_obj[[name]]))
  }))
  gsc <- data.frame(gene = unlist(gsc_obj), gene_set = gs_names)

  res <- runGSAhyper(genes = hit_genes,
                     universe = all_genes,
                     gsc = loadGSC(gsc),
                     gsSizeLim = gsSizeLim,
                     adjMethod = adjMethod) %>%
    .[['resTab']] %>%
    data.frame(check.names=FALSE) %>%
    tibble::rownames_to_column(var = 'gene_set') %>%
    mutate(
      set_size = `Significant (in gene set)` + `Non-significant (in gene set)`,
      prob_in_set = `Significant (in gene set)` / (`Significant (in gene set)` + `Non-significant (in gene set)`),
      prob_not_in_set = `Significant (not in gene set)` / (`Significant (not in gene set)` + `Non-significant (not in gene set)`),
      odds_ratio = (prob_in_set / (1 - prob_in_set)) / (prob_not_in_set / (1 - prob_not_in_set))
    ) %>%
    dplyr::select(-prob_in_set, -prob_not_in_set)
  return(res)
}



#' Make network visualization of gene set stats
#'
#' @param gene_stats dataframe containing at minimum columns 'gene_set' and 'stat'
#' @param gsc Original gene set collection object
#' @param all_genes List of all genes used in original analysis
#' @param top_N number of top gene sets to visualize in the graph (ranked by stat)
#' @param one_per_cluster Whether to label just the top gene set per cluster (default TRUE)
#' @param min_cluster_size Minimum cluster size to include in the final plot (default 2)
#' @param vertex.label.cex Size of vertex labels (default 0.5)
#' @param set_cor_thresh Minimum Jaccard similarity between sets to include as edge
#'
#' @return igraph plot
#' @export
#'
#' @examples
make_GSA_network <- function(gene_set_stats,
                             gsc,
                             all_genes,
                             top_N = 100,
                             one_per_cluster = TRUE,
                             min_cluster_size = 2,
                             vertex.label.cex = 0.5,
                             set_cor_thresh = 0.1) {
  library(philentropy)
  library(igraph)
  node_df <- gene_set_stats %>%
    arrange(desc(stat)) %>%
    head(top_N) %>%
    mutate(gene_set = as.vector(gene_set))
  node_df[['id']] <- node_df$gene_set
  node_df[['size']] <- node_df$stat / sd(node_df$stat)

  gsc_M <- cdsr::gsc_to_sparseMat(gsc[as.vector(node_df$gene_set)], all_genes)

  C_mat <- philentropy::distance(t(as.matrix(gsc_M[, node_df$gene_set])), method = 'jaccard')
  colnames(C_mat) <- node_df$gene_set
  rownames(C_mat) <- node_df$gene_set
  edge_df <- melt(1-C_mat) %>%
    set_colnames(c('from', 'to', 'cor')) %>%
    filter(from != to,
           cor >= set_cor_thresh)

  #make initial graph to get neighborhoods
  graph <- igraph::graph_from_data_frame(d = edge_df,
                                         vertices = node_df,
                                         directed = FALSE)
  clust <- igraph::cluster_edge_betweenness(graph, E(graph)$cor, directed = TRUE)
  node_df$clust <- clust$membership
  clust_sizes <- sapply(groups(clust), length)
  good_clusts <- which(clust_sizes >= min_cluster_size)
  node_df %<>% filter(clust %in% good_clusts)
  edge_df %<>% filter(from %in% node_df$gene_set, to %in% node_df$gene_set)

  if (nrow(node_df) == 0) {
    return(NULL)
  }

  #only retain labels for best node per group
  node_labels <- node_df %>%
    group_by(clust) %>%
    dplyr::slice(which.max(stat)) %>%
    ungroup()

  graph <- igraph::graph_from_data_frame(d = edge_df,
                                         vertices = node_df,
                                         directed = FALSE)

  clust_to_num <- seq(length(good_clusts))
  names(clust_to_num) <- good_clusts

  if (length(good_clusts) <= 12) {
    node_df$color <- RColorBrewer::brewer.pal(length(good_clusts), "Set3") %>%
      magrittr::extract(clust_to_num[as.character(node_df$clust)])
  } else {
    node_df$color <- rainbow(length(good_clusts)) %>%
      magrittr::extract(clust_to_num[as.character(node_df$clust)])
  }

  if (one_per_cluster) {
    node_df %<>% mutate(node_label = ifelse(gene_set %in% node_labels$gene_set, gene_set, NA))
  } else {
    node_df %<>% mutate(node_label = gene_set)
  }

  vertex_size <- node_df[["size"]] %>%
    as.numeric() %>%
    abs() %>%
    magrittr::divide_by(max(.)) %>%
    magrittr::multiply_by(20)

  igraph::plot.igraph(graph,
                      vertex.label.dist = 0,
                      vertex.label = node_df$node_label,
                      vertex.label.color = 'black',
                      vertex.size = vertex_size,
                      vertex.alpha = 0.75,
                      vertex.color = adjustcolor(node_df$color, alpha.f=0.75),
                      vertex.shape = 'circle',
                      edge.arrow.size = 0.1,
                      edge.curved = 0.1,
                      layout = layout_nicely,
                      vertex.label.cex = vertex.label.cex,
                      edge.width=10*edge_df$cor)

  }


# compute gene set analysis ------------
geneset_matsum <- function(mx, mx.msigdb){ mx %*% mx.msigdb }
geneset_matmean <- function(matsum, geneset.size){
  one_over_n <- matrix(1/geneset.size, nrow(matsum), ncol(matsum), byrow = T)
  matsum * one_over_n
}

geneset_matvar <- function(mx.sumsqr, mx.mean, geneset.size){
  n <- matrix(geneset.size, nrow(mx.sumsqr), ncol(mx.sumsqr), byrow = T)
  one_over_nminus1 <- matrix(1/(geneset.size - 1), nrow(mx.sumsqr), ncol(mx.sumsqr), byrow = T)

  (mx.sumsqr - mx.mean^2 * n) * one_over_nminus1
}

geneset_t <- function(mx, mx.msigdb){

  geneset.size <- colSums(mx.msigdb)


  sum_within <- geneset_matsum(mx, mx.msigdb) %>% data.matrix()
  mean_within <- geneset_matmean(sum_within, geneset.size)
  sum_bkgd <- matrix(rowSums(mx), nrow = nrow(mx), ncol = ncol(mx.msigdb)) - sum_within
  mean_bkgd <- geneset_matmean(sum_bkgd, ncol(mx) - geneset.size) %>% data.matrix()

  sumsqrt_within <- geneset_matsum(mx^2, mx.msigdb) %>% data.matrix()
  var_within <- geneset_matvar(sumsqrt_within, mean_within, geneset.size) %>% data.matrix()
  sumsqrt_bkgd <- matrix(rowSums(mx^2), nrow = nrow(mx), ncol = ncol(mx.msigdb)) - sumsqrt_within
  var_bkgd <- geneset_matvar(sumsqrt_bkgd, mean_bkgd, ncol(mx) - geneset.size) %>% data.matrix()

  delta_mean <- mean_within - mean_bkgd
  pooled_var <- var_within/geneset.size + var_bkgd/(ncol(mx) - geneset.size)

  melt(mean_within, varnames = c('dataset', 'geneset.name'), value.name = 'mean.corr') %>%
    mutate(t.statistic = c(delta_mean/sqrt(pooled_var)))
}


t_test_msigdb <- function(cormat, mx.msigdb, min.geneset.size = 10){
  cormat %<>% set_colnames(gsub('.*_(.*)', '\\1', colnames(.)))

  common_genes <- intersect(colnames(cormat), rownames(mx.msigdb))
  cormat <- cormat[, common_genes]

  mx.msigdb <- mx.msigdb[common_genes, ]

  # exclude gene sets less than 10 genes or more than 500 genes
  mx.msigdb <- mx.msigdb[, colSums(mx.msigdb) <= 500 & colSums(mx.msigdb) >= min.geneset.size]

  # compute t statistic using matrix multiplications
  geneset_t(cormat, mx.msigdb)
}

find_geneset_category <- function(x){
  x <- ifelse(grepl('Pertb\\d+', x), 'MSigDB Perturbation', x)
  x <- ifelse(grepl('GO\\d+', x), 'MSigDB GO', x)
  x <- ifelse(grepl('Pthwy\\d+', x), 'MSigDB Pathway', x)
  x <- ifelse(grepl('Hallmk\\d+', x), 'MSigDB Hallmark', x)
  x <- ifelse(grepl('target_', x), 'REP target gene', x)
  x <- ifelse(grepl('MOA_', x), 'REP MOA', x)
  x
}


get_ttest_table <- function(cor_table, data_msigdb, data_geneset, feature){

  cormat <- cor_table %>%
    reshape2::acast(Target ~ dataset, value.var = 'Pearson.corr', fun.aggregate = median) %>%
    t()


  t_test_msigdb(cormat, data_msigdb) %>%
    left_join(cor_table %>% select(dataset, global.mean, global.sd) %>% unique(), by = 'dataset') %>%
    left_join(data_geneset, by = 'geneset.name') %>%
    mutate(feature = feature,
           mean.Z = (mean.corr - global.mean)/global.sd,
           category = find_geneset_category(geneset.name)) %>%
    select(feature, dataset, category, geneset, mean.corr, mean.Z, t.statistic, geneset.name) %>%
    filter(!grepl('NIKOLSKY|CHARAFE', geneset)) %>%
    arrange(desc(abs(mean.Z)))
}

