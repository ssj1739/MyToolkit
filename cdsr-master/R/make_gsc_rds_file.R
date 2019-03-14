
# get MSIGDB gene set data
# [LW] I put the following code inside a function to avoid auto-loading libraries and generate conflicts in NAMESPACE
make_MSigDB_gs_data <- function(){
  library(readr)
  library(magrittr)
  library(tibble)
  library(reshape2)
  library(GSEABase)
  
  gsc_dir <- '~/CPDS/data/MSigDB'
  
  gene_set_info <- list(
    hallmark = list(
      file = 'h.all.v6.0.symbols.gmt',
      col_type = 'h',
      key_type = 'Gene'),
    positional = list(
      file = 'c1.all.v6.0.symbols.gmt',
      col_type = 'c1',
      key_type = 'Gene'),
    oncogenes = list(
      file = 'c6.all.v6.0.symbols.gmt',
      col_type = 'c6',
      key_type = 'Gene'),
    GO_biological_process = list(
      file = 'c5.bp.v6.0.symbols.gmt',
      col_type = 'c5',
      key_type = 'Gene'),
    GO_cellular_component = list(
      file = 'c5.cc.v6.0.symbols.gmt',
      col_type = 'c5',
      key_type = 'Gene'),
    GO_molecular_function = list(
      file = 'c5.mf.v6.0.symbols.gmt',
      col_type = 'c5',
      key_type = 'Gene'),
    canonical = list(
      file = 'c2.cp.v6.0.symbols.gmt',
      col_type = 'c2',
      key_type = 'Gene'),
    reactome = list(
      file = 'c2.cp.reactome.v6.0.symbols.gmt',
      col_type = 'c2',
      key_type = 'Gene'),
    biocarta = list(
      file = 'c2.cp.biocarta.v6.0.symbols.gmt',
      col_type = 'c2',
      key_type = 'Gene'),
    KEGG = list(
      file = 'c2.cp.kegg.v6.0.symbols.gmt',
      col_type = 'c2',
      key_type = 'Gene'),
    cancer_modules = list(
      file = 'c4.cm.v6.0.symbols.gmt',
      col_type = 'c2',
      key_type = 'Gene'),
    CORUM = list(
      file = 'corum.all.symbols.dedup.gmt',
      col_type = 'other',
      key_type = 'Gene'),
    Transc_facs = list(
      file = 'c3.tft.v6.0.symbols.gmt',
      col_type = 'c3',
      key_type = 'Gene'),
    Rep_moa = list(
      file = 'rep_moa.gmt',
      col_type = 'other',
      key_type = 'Broad_ID'),
    Rep_target = list(
      file = 'rep_target.gmt',
      col_type = 'other',
      key_type = 'Broad_ID'),
    GDSC_target = list(
      file = 'gdsc_target.gmt',
      col_type = 'other',
      key_type = 'DRUG_ID'
    ),
    GDSC_target_pathway = list(
      file = 'gdsc_target_pathway.gmt',
      col_type = 'other',
      key_type = 'DRUG_ID'
    ),
    CTD2_target = list(
      file = 'ctd2_target.gmt',
      col_type = 'other',
      key_type = 'CPD_ID'
    ),
    CTD2_moa = list(
      file = 'ctd2_moa.gmt',
      col_type = 'other',
      key_type = 'CPD_ID'
    )
  )
  
  gene_set_info
}

make_gsc_data <- function(gene_set_info, gsc_dir) {
  gsc_data <- lapply(gene_set_info, function(gs) {
    if (gs$col_type == 'other') {
      collect <- ComputedCollection()
    } else {
      collect <- BroadCollection(category = gs$col_type)
    }
    getGmt(file.path(gsc_dir, gs$file),
           collectionType = collect,
           geneIdType = SymbolIdentifier())
  })

  #make chrome-arm positional gsc
  all_chroms <- paste0('chr', c(seq(22), c('X', 'Y')))
  all_arms <- apply(expand.grid(all_chroms, c('p', 'q')),1, paste0, collapse = '')
  gsc_data$chr_arms <- lapply(all_arms %>% set_names(all_arms), function(arm_name) {
    print(arm_name)
    used_gscs <- grepl(arm_name, names(gsc_data$positional))
    if (sum(used_gscs) > 0) {
      temp_set <- Reduce(GSEABase::union, gsc_data$positional[used_gscs])
      setName(temp_set) <- arm_name
      return(temp_set)
    }
  })
  gsc_data$chr_arms <- gsc_data$chr_arms[sapply(gsc_data$chr_arms, function(x) !is.null(x))]
  gsc_data$chr_arms %<>% GeneSetCollection()
  gene_set_info$chr_arms <- list(
    file = NA,
    col_type = 'c1',
    key_type = 'Gene')

  attr(gsc_data, 'key_type') <- lapply(gene_set_info, function(x) x$key_type)
  write_rds(gsc_data, file.path(gsc_dir, 'gsc_data.rds'))
}
