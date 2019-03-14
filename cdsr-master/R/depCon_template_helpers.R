get_mat_vec_cors <- function(mat, vec, mat_SE = NULL) {
  if (is.null(mat_SE)) {
    return(cor(mat, vec, use = 'pairwise.complete.obs')[,1])
  } else {
    library(weights)
    CL_weights = 1/rowMeans(mat_SE^2, na.rm=T)
    return(wtd.cors(mat, vec, weight = CL_weights)[,1])
  }
}

two_to_one_sided <- function(two_sided_p, stat, test_dir) {
  one_sided_p <- two_sided_p / 2
  if (test_dir == 'right') {
    one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
  } else {
    one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
  }
  return(one_sided_p)
}

get_limma_stats <- function(dep_mat, CL_vec, covars = NULL, dep_SE_mat = NULL, target_type = 'Gene', test_dir = 'left') {
  common_CLs <- intersect(names(CL_vec), rownames(dep_mat))
  if (!is.null(covars)) {
    common_CLs %<>% intersect(names(covars))
  }
  if (is.null(dep_SE_mat)) {
    weights = NULL
  } else {
    weights <- 1/dep_SE_mat[common_CLs,]^2
  }
  lim_res <- cdsr::run_lm_stats_limma(mat = dep_mat[common_CLs,],
                                      vec = CL_vec[common_CLs],
                                      covars = covars[common_CLs],
                                      weights = weights,
                                      target_type = 'Gene')

  if (test_dir == 'both') {
    lim_res %<>% mutate(p_val = p.value, q_val = q.value)
  } else {
    lim_res %<>% mutate(p_val = two_to_one_sided(p.value, EffectSize, test_dir),
                        q_val = p.adjust(p_val, method = 'BH'))
  }

  pear_cors <- get_mat_vec_cors(dep_mat[common_CLs,], CL_vec[common_CLs], dep_SE_mat[common_CLs,])
  lim_res %<>% left_join(data.frame(Gene = names(pear_cors), Pear_cor = pear_cors), by = 'Gene') %>%
    dplyr::select(Gene, EffectSize, Pear_cor, p_val, q_val)
  if (target_type != 'Gene') {
      names(lim_res)[1] <- target_type
  }

  return(lim_res)
}


make_dep_clusters <- function(mat, stats_df, max_genes = 50, two_sided = FALSE,
                              var_type = 'Gene', key_var = 'brd_id',
                              cor_p_thresh = 0.01, min_cor = 0.25) {
  top_genes <- filter(stats_df, pref_dep) %>%
    head(max_genes)
  if (var_type == 'cmpd_name' & key_var == 'brd_id') {
    top_genes %<>% mutate(brd_id = make.names(str_sub(brd_id, 1,13)))
  }
  if (nrow(top_genes) > 2) {
    node_df <- top_genes
    node_df[['id']] <- node_df[[var_type]]
    node_df[['label']] <- node_df[[var_type]]
    node_df[['size']] <- node_df[['MeanDiff']]
    if (two_sided) {
      node_df %<>% mutate(size = abs(size))
    }
    if (var_type == 'cmpd_name') {
      C_mat <- psych::corr.test(mat[, top_genes[[key_var]]], use = 'pairwise.complete.obs')
    } else {
      C_mat <- psych::corr.test(mat[, top_genes[[var_type]]], use = 'pairwise.complete.obs')
    }
    C_mat$r[C_mat$p > cor_p_thresh] <- 0
    edge_df <- melt(C_mat$r) %>%
      set_colnames(c('from', 'to', 'cor')) %>%
      filter(from != to,
             abs(cor) >= min_cor)

    invisible(make_a_network(node_df, edge_df, directed=F))
  } else {
    return(NULL)
  }
}

get_gsea_results <- function(dep_mat, CL_var, gsc, covars, test_dir, target_type = 'Gene') {
  if (target_type == 'broad_id') {
    colnames(dep_mat) <- make.names(str_sub(colnames(dep_mat), 1, 13))
  } else if (target_type == 'DRUG ID') {
    colnames(dep_mat) <- paste0('ID', colnames(dep_mat))
  } else if (target_type == 'master_cpd_id') {
    colnames(dep_mat) <- paste0('X', colnames(dep_mat))
  }
  gsea_res <- cdsr::run_npGSEA(gsc,
                              dep_mat,
                               CL_var,
                               covars = covars,
                               min_set_size = 3,
                               max_set_size = Inf,
                               normalize_X = FALSE,
                               n_top_genes = 5,
                               top_dir = test_dir)
  if (test_dir == 'both') {
    gsea_res %<>% dplyr::mutate(p_val = p_twoSided, q_val = q_twoSided)
  } else if (test_dir == 'left') {
    gsea_res %<>% dplyr::mutate(p_val = p_left, q_val = q_left)
  } else if (test_dir == 'right') {
    gsea_res %<>% dplyr::mutate(p_val = p_right, q_val = q_right)
  }
  gsea_res %<>%
    dplyr::select(gene_set, p_val, q_val, set_size, top_genes, norm_enrich_stat) %>%
    arrange(p_val) %>%
    mutate(p_val = cdsr:::format_sci_not(p_val, digits = 2),
           q_val = cdsr:::format_sci_not(q_val, digits = 2))
  return(gsea_res)
}

make_network_target_vis <- function(top_gene_stats,
                                  mat,
                                  min_cluster_size = 1,
                                  vertex.label.cex = 1,
                                  target_type = 'Gene') {
  library(igraph)

  cor_p_thresh <- 0.05
  min_cluster_size <- 1
  cor_thresh_quantile <- 0.9

  node_df <- top_gene_stats
  node_df[['id']] <- node_df[[target_type]]
  node_df[['label']] <- node_df[[target_type]]
  node_df[['size']] <- abs(node_df[['EffectSize']])

  if (target_type == 'broad_id') {
    C_mat <- psych::corr.test(mat[, top_gene_stats[[target_type]]], use = 'pairwise.complete.obs')
    node_df[['label']] <- node_df$pert_iname
  } else if (target_type == 'DRUG ID') {
    C_mat <- psych::corr.test(mat[, top_gene_stats[[target_type]]], use = 'pairwise.complete.obs')
    node_df[['label']] <- node_df$`DRUG NAME`
  } else if (target_type == 'master_cpd_id') {
    C_mat <- psych::corr.test(mat[, top_gene_stats[[target_type]]], use = 'pairwise.complete.obs')
    node_df[['label']] <- node_df$cpd_name
  } else {
    C_mat <- psych::corr.test(mat[, top_gene_stats[[target_type]]], use = 'pairwise.complete.obs')
  }
  min_cor <- quantile(C_mat$r, cor_thresh_quantile, na.rm=T)

  C_mat$r[C_mat$p > cor_p_thresh] <- 0
  edge_df <- melt(C_mat$r) %>%
    set_colnames(c('from', 'to', 'cor')) %>%
    mutate(cor = abs(cor)) %>%
    filter(from != to, cor >= min_cor)

  #make initial graph to get neighborhoods
  graph <- igraph::graph_from_data_frame(d = edge_df,
                                         vertices = node_df,
                                         directed = FALSE)
  clust <- igraph::cluster_edge_betweenness(graph, E(graph)$cor, directed = TRUE)
  node_df$clust <- clust$membership
  clust_sizes <- sapply(groups(clust), length)
  good_clusts <- which(clust_sizes >= min_cluster_size)
  node_df %<>% filter(clust %in% good_clusts)
  edge_df %<>% filter(from %in% node_df[[target_type]], to %in% node_df[[target_type]])

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
                      edge.width=20*edge_df$cor)
}


create_link_gene_page <- function(val) {
  sprintf('<a href="https://cds.team/depmap/gene/%s">%s</a>', val, val)
}
create_link_taiga_page <- function(data_set_name, data_version) {
  sprintf('<a href="https://cds.team/taiga/dataset/%s/%s">%s</a>', data_set_name, data_version, data_set_name)
}
create_link_compound_page <- function(brd_id, name) {
  sprintf('<a href="http://datasci-dep:3838/compoundReport/?targetID=%s">%s</a>', brd_id, name)
}


#combine top model info, and best feature info
proc_gene_biomarker_data <- function(df) {
  join_all(list(
    df %>%
      group_by(gene) %>%
      top_n(n = 1, wt = overall_auc) %>%
      ungroup() %>%
      dplyr::select(gene, best_model = model, overall_auc, overall_pearson),
    df %>%
      filter(model == 'KitchenSink') %>%
      dplyr::select(gene, best_overall_feature = feature0),
    df %>%
      filter(model == 'RelatedFeatures') %>%
      dplyr::select(gene, best_related_feature = feature0)
  ), by = 'gene', type = 'full') %>%
    mutate(best_overall_feature = str_replace_all(best_overall_feature, ' \\(.+\\)', ''),
           best_related_feature = str_replace_all(best_related_feature, ' \\(.+\\)', '')) %>%
    dplyr::rename(Gene = gene)
}

get_df_90_rank <- function(dep_df, target_name) {
  #calculate 90th percentile of gene dep ranks per cell line for each gene/cmpd
  df <- dep_df %>%
    melt() %>%
    set_colnames(c('CCLE_ID', 'target', 'GS')) %>%
    dplyr::group_by(CCLE_ID) %>%
    mutate(gene_rank = percent_rank(GS)) %>%
    ungroup() %>%
    dplyr::group_by(target) %>%
    dplyr::summarise(deprank_90perc = quantile(gene_rank, probs = 0.9, na.rm=T)) %>%
    dplyr::rename_(.dots = setNames('target', target_name))
  df[[target_name]] <- as.character(df[[target_name]])
  return(df)
}

format_sci_not <- function(vec, digits = 2) {
  exp <- floor(log10(vec))
  rnums <- round(vec/10^exp, digits-1)
  return(rnums*10^exp)
}

make_pref_dep_table <- function(dep_stats, target_info, target_type, test_dir, max_pref_deps, min_show_top_deps, q_thresh, normLRT_thresh = NULL) {
  if (!is.null(dep_stats)) {
    if (test_dir == 'both') {
      dep_stats %<>% mutate(esize = abs(EffectSize))
    } else if (test_dir == 'left') {
      dep_stats %<>% mutate(esize = -EffectSize)
    } else {
      dep_stats %<>% mutate(esize = EffectSize)
    }
    dep_stats %<>% arrange(dplyr::desc(esize))

    pref_deps <- dep_stats %>%
      filter(q_val < q_thresh)

    n_pref_deps <- nrow(pref_deps)
    if (n_pref_deps < min_show_top_deps) {
      pref_deps <- dep_stats %>%
        head(min_show_top_deps)
    }
    pref_deps %<>% head(max_pref_deps)

    pref_deps %<>%
      left_join(target_info, by = target_type) %>%
      mutate(p_val = format_sci_not(p_val, digits = 2),
             q_val = format_sci_not(q_val, digits = 2)) %>%
      dplyr::select(-esize)

    if (!is.null(normLRT_thresh)) {
      pref_deps %<>% filter(normLRT > normLRT_thresh)
      #move normLRT to second col position
      pref_deps <- pref_deps[, c(target_type, 'normLRT', setdiff(colnames(pref_deps), c(target_type, 'normLRT')))]
    }
    if (target_type == 'broad_id') {
      pref_deps %<>% dplyr::mutate(broad_id = pert_iname) %>%
        dplyr::rename(Compound = broad_id)
    } else if (target_type == 'DRUG ID') {
      pref_deps %<>% dplyr::mutate(`DRUG ID` = `DRUG NAME`) %>%
        dplyr::rename(Compound = `DRUG ID`)
    } else if (target_type == 'master_cpd_id') {
      pref_deps %<>% dplyr::mutate(master_cpd_id = cpd_name) %>%
        dplyr::rename(Compound = master_cpd_id)
    }
    if (target_type == 'Gene') {
        pref_deps %<>% mutate(Gene = create_link_gene_page(Gene))
    }
    # } else if (target_type == 'Compound') {
    #     pref_deps %<>% mutate(Compoud = create_link_gene_page(Compoud))
    # }
    n_rows_used <- nrow(pref_deps)
    tab <- DT::datatable(pref_deps, filter = 'top', rownames = FALSE, escape = FALSE,
                  options = list(lengthMenu = c(5, 10, 25, 50, 100), pageLength = 25)) %>%
      DT::formatRound(intersect(colnames(pref_deps),
                                c('EffectSize', 'Pear_cor', 'overall_auc', 'overall_pearson', 'deprank_90perc',
                                  'in_group_mean_dep', 'normLRT', 'in_group_frac_FDR_deps', 'frac_FDR_deps', 'in_group_pref_dep_frac')), 3)
  }
  return(list(tab = tab, n_pref_deps = n_pref_deps, n_rows_used = n_rows_used))
}

make_volcano_plot <- function(pref_dep_stats, test_dir, n_targ_label = 20, target_type = 'Gene', q_thresh = 0.2) {
  if (test_dir == 'both') {
    pref_dep_stats %<>% mutate(esize = abs(EffectSize))
  } else if (test_dir == 'left') {
    pref_dep_stats %<>% mutate(esize = -EffectSize)
  } else {
    pref_dep_stats %<>% mutate(esize = EffectSize)
  }
  if (target_type %in% c('broad_id', 'DRUG ID', 'master_cpd_id')) {
    cur_target_type <- 'Compound'
  } else {
    cur_target_type <- target_type
  }
  ggplot(pref_dep_stats, aes(EffectSize, -log10(p_val), color = q_val < q_thresh)) +
    geom_point(alpha = 0.75) +
    geom_text_repel(data = pref_dep_stats %>%
                      filter(q_val < q_thresh) %>%
                      arrange(dplyr::desc(esize)) %>%
                      head(n_targ_label),
                    aes_string('label' = cur_target_type),
                    size = 3) +
    guides(color = FALSE) +
    theme_Publication() +
    scale_colour_Publication()
}


theme_Publication <- function(base_size=12, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.text = element_text(size = rel(1.2)),
            legend.direction = "horizontal",
            legend.key.size= unit(0.3, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(5,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#fdb462","#386cb0","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#fdb462","#386cb0","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

