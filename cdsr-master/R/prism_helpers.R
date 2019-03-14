#' @export load.from.taiga.aws
load.from.taiga.aws <- function(data.name, data.version){
  path = '~/.taiga/'
  filename <- paste(path, data.name, '_', data.version, '.Rdata', sep = '')
  load(filename)
  data
}

#'  viability transformation
#'  @export log_viab100
log_viab100 <- function(x, cap = 10){
  x[x < cap] <- cap
  log2(x/100)
} %>% Vectorize()

#'  @export log_viab
log_viab <- function(x, cap = 0.1){
  x[x < cap] <- cap
  log2(x)
} %>% Vectorize()

#'  @export cap
cap <- function(x, max = 100, min = 0){
  x <- ifelse(x > max, max, x)
  x <- ifelse(x < min, min, x)
  x
}

#'  @export saveToCSV
saveToCSV <- function(data.name, data.version){
  x <- load.from.taiga.aws(data.name = data.name, data.version = data.version)
  write.csv(x, paste0("data/toTaiga/", data.name, '.csv'), row.names = T)
}

#'  @export saveToTaiga
saveToTaiga <- function(x, data.name, data.version){
  data <- x
  save(data, file = paste0("~/Documents/Broad/PRISM/DepModeling/data/toTaiga/", data.name, '_', data.version, '.Rdata'))
}

#'  @export get_binary_mut_mat
get_binary_mut_mat <- function(df){
  mx <- df %>% mutate(value = 1) %>%
    reshape2::dcast(CCLE_ID ~ Hugo_Symbol, value.var = 'value', fun.aggregate = sum) %>%
    magrittr::set_rownames(.[['CCLE_ID']]) %>% select(-CCLE_ID) %>%
    data.matrix()
  mx[mx > 1] <- 1
  mx
}

#' take first column as rownames and convert to matrix
#'  @export to_mat
to_mat <- function(df){
  df %>% magrittr::set_rownames(df[, 1]) %>% .[, -1, drop = F] %>% data.matrix()
}

#' take rownames as the first column and convert to data frame
#'  @export to_df
to_df <- function(mx, colname1){
  cbind(data.frame(as.character(rownames(mx))) %>% magrittr::set_colnames(colname1),
        as_data_frame(mx))
}

#' convert pert_dose to standard dose dilution
#'  @export to_df
get_effective_dose <- function(pert_dose, effective_dose){
  pert_effective_dose <- pert_dose
  idx_cpd <- which(pert_dose > 0)

  pert_effective_dose[idx_cpd] <-
    kmeans(log(pert_dose[idx_cpd]), centers=log(effective_dose)) %>%
    magrittr::extract2("cluster") %>%
    magrittr::extract(effective_dose, .)
  pert_effective_dose
}

#' Ari selectivity score for PRISM profile
#'  @export ari_score
#'  @param x viability, must be in 0 ~ 100 format, i.e. capped from 0 to 100
ari_score <- function(x){

  #set anything below 0 to zero and anything above 100 to 100
  x <- x[!is.na(x)]
  x[x < 0] <- 0
  x[x > 100] <- 100
  x <- sort(x)

  # the first 5% quantile
  x_s1 <- x[1: floor(length(x)*0.05)]
  x_s1 <- sum(100 - 4*x_s1, na.rm = T)

  # the 20% to 40% quantile
  x_s2 <- x[ceiling(length(x)*0.2) : floor(length(x)*0.4)]
  x_s2 <- sum(x_s2 - 100, na.rm = T)

  # 40% to 100% quantile
  x_s3 <- x[ceiling(length(x)*0.4) : length(x)]
  x_s3 <- sum(x_s3 - 100, na.rm = T)

  (8*x_s1 + 3*x_s2 + 4*x_s3)/100

}

#' Compute correlation matrix for two matrices
#' @export get_cormat
#' @param mx1 first input matrix, cannot be NA
#' @param mx2 second input matrix, can be NA. When mx2 is NA, it will compute self correlation for mx1
#' @param Z.score one of strings "by mx1", by mx2", or "global", it will transform each correlation value to Z scores
get_cormat <- function(mx1, mx2, Z.score = NA){
  stopifnot(class(mx1) == 'matrix')

  if(class(mx1) == 'logical') {mx2 <- mx1} # self correlation if mx2 is not provided

  rows <- intersect(rownames(mx1), rownames(mx2))
  cormat <- suppressWarnings(cor(mx1[rows, ], mx2[rows, ], use = 'p'))

  if(is.na(Z.score)){
    return(cormat)
  }
}


#' convert dose to character
#'  @export make_dose
make_dose <- function(x){
  ifelse(x >= 1, paste0(round(x, 2), 'uM'),
         ifelse(x > 0, paste0(round(1000*x, 0), 'nM'), 'DMSO'))
}

#' make diagonal value NA for self correlation matrix
#' @export set_diag_na
#' @param cormat square matrix with diagnal values equal to one
set_diag_na <- function(mx){
  stopifnot(class(mx) == 'matrix')
  stopifnot(nrow(mx) == ncol(mx))
  stopifnot(all(diag(mx) == 1))

  mx[mx == 1] <- NA
  mx
}

#' run combat for batch correction
#' @export run_combat
#' @param df_data data frame containing columns of "pert_type", "pert_mfc_id", "pert_effective_dose", "det_plate", "cell_line_id",  and "pool_id"
#' @param value.var colnames of the viability column
run_combat <- function(df_data, value.var, output.name = 'norm_viab_combat'){
  df_combat <- df_data %>%
    filter(pert_type == 'trt_cp') %>%
    plyr::ddply(plyr::.(pert_mfc_id, pert_well), function(x){
      print(paste(x$pert_mfc_id[1], x$pert_effective_dose[1]))
      dat_mat <- x %>%
        mutate(condition_id=paste(pert_well, det_plate, pool_id, sep=":")) %>%
        reshape2::dcast(cell_line_id + condition_id ~ pert_mfc_id,
                        value.var=value.var)
      batch <- dat_mat$condition_id
      cell_ids <- dat_mat$cell_line_id
      bat_mat <- try(sva::ComBat(dat_mat %>%
                                   select(-cell_line_id, -condition_id) %>%
                                   as.matrix() %>%
                                   t() %>%
                                   rbind(dummy=rnorm(ncol(.))),
                                 batch=batch))
      if(class(bat_mat) != "try-error"){
        data.frame(cell_line_id=cell_ids,
                   condition_id=batch,
                   t(bat_mat[!grepl("dummy", row.names(bat_mat)), , drop=F])) %>%
          magrittr::set_colnames(gsub("[.]", "-", colnames(.))) %>%
          reshape2::melt(id.vars=c("cell_line_id", "condition_id"),
                         variable.name="pert_mfc_id",
                         value.name = output.name) %>%
          tidyr::separate("condition_id",
                          c("pert_well", "det_plate", 'pool_id'),
                          sep=":")
      } else {
        NULL
      }

    }, .parallel = T)


}


#' general function to compute FCPC
#' removing the top2 doses
#' @export compute_FCPC
#' @param df_mfi data frame containing columns of mfi, det_plate, pool_id, cell_line_id
compute_FCPC <- function(df_mfi){
  df_FCPC <- df_mfi %>%
    # Normalize the MFI (log version) values to plate control
    group_by(det_plate, pool_id, cell_line_id) %>%
    mutate(log_fc=log2(mfi / median(mfi[pert_type != "trt_poscon" &
                                                      round(pert_effective_dose, 0) < 1], na.rm=T))) %>%
    ungroup()
}


#' general function to compute ssmd
#' @export compute_ssmd
compute_ssmd <- function(df_FCPC){
  df_ssmd <- df_FCPC %>%

    # Perform SSMD QC
    group_by(det_plate, pool_id, cell_line_id) %>%
    summarise(ssmd=(median(log_fc[pert_type=="ctl_vehicle"], na.rm=T) -
                      median(log_fc[pert_type=="trt_poscon"], na.rm=T)) /
                sqrt(mad(log_fc[pert_type=="ctl_vehicle"], na.rm=T)^2 +
                       mad(log_fc[pert_type=="trt_poscon"], na.rm=T)^2),
              plate_ssmd=(median(log_fc[pert_type != "trt_poscon" & round(pert_effective_dose, 0) < 5], na.rm=T) -
                            median(log_fc[pert_type=="trt_poscon"], na.rm=T)) /
                sqrt(mad(log_fc[pert_type != "trt_poscon" &
                                  round(pert_effective_dose, 0) < 5], na.rm=T)^2 +
                       mad(log_fc[pert_type=="trt_poscon"], na.rm=T)^2)) %>%
    ungroup() %>%
    mutate(ssmd=pmax(ssmd, plate_ssmd))
}

#' a ggplot theme that removes grid, background color, and top and right borders
#' @export theme_custom
theme_custom <- function(){
  theme_bw()  %+replace%
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}

#' get the largest absolute value with sign
#' @export get_best
get_best <- function(x){
  ifelse(max(x, na.rm = T) > -min(x, na.rm = T), max(x, na.rm = T), min(x, na.rm = T))
}

#' query dbPRISM for compound assay data
#' @export dbPRISM_get_assaydata
dbPRISM_get_assaydata <- function(Broad_ID, project_ID){
  # match project ID
  project_list <- c('CalicoT1', 'Repurposing', 'MTS_Rep', 'CalicoT2A', 'REP_PR500')
  print(paste('Valid project IDs are', paste(project_list, collapse = ', ')))
  stopifnot(project_ID %in% project_list)
  
  require(RPostgreSQL)
  require(dbplyr)
  
  prismdb <- src_postgres(dbname = 'prism',
                          host = 'prism.ceiatqqqcshl.us-east-1.rds.amazonaws.com',
                          port = 5432,
                          user = 'prismadmin',
                          password = '?Ct7)4exFz')
  
  
  projectID_list <- list(CalicoT1 = 'HTS_CalicoT1dose_1119pert_576cl_2016', 
                         Repurposing = 'HTS_REPdose_1400pert_576cl_2016', 
                         MTS_Rep = 'MTS_REP_539pert_562cl_2017Q2', 
                         CalicoT2A = 'HTS_CalicoT2A_2884pert_489cl_2017Q2',
                         REP_PR500 = 'HTS_REPdose_redetection_1400pert_489cl_2017Q3')
  
  projectID <- projectID_list[[project_ID]]
  query <- paste0('SELECT assay_cell_line_id, perturbation_id,rep_idx, effective_dose, viability FROM assay_data WHERE perturbation_id =',
                  '\'', Broad_ID, '\'', ' AND project_id = ', '\'', projectID, '\'')
  
  assaydata_colnames <- c('cell_line_id', 'pert_mfc_id', 'replicate_idx', 'pert_effective_dose', 'viability')
  print(query)
  assaydata_postgres <- tbl(prismdb, sql(query)) %>% 
    as.data.frame() %>% 
    set_colnames(assaydata_colnames)
  
  # convert viability to 0~100 scale for consistency
  # if(project == 'MTS_Rep'){assaydata_postgres %<>% mutate(viability = 100*viability)}
  assaydata_postgres
}


#' process gctx format PRISM data
#' @export process_gctx_prism
#' @param file.name a character string of the path of the input file
#' @param exlucde.wells well to exclude, a vector of platewell_id, i.e. det_plate:pert_well
process_gctx_prism <- function(file.name, exlucde.wells){
  require(rhdf5)
  
  gct <- suppressMessages(roller::parse.gctx(file.name))
  raw_data <- melt(gct@mat) %>% 
    rename(platewell_id = Var2, cell_line_id = Var1) %>%
    separate(platewell_id, into = c('det_plate', 'pert_well'), sep = ':', remove = F) %>%
    left_join(cl_meta_pr500) %>%
    left_join(cpd_meta) %>%
    # exlude wells
    mutate(value = ifelse(platewell_id %in% exlucde.wells, NA, value)) %>%
    # two compound has lost dose information, exclude mislabled doses
    filter(!(pert_type == 'trt_cp' & pert_effective_dose < 0)) %>%
    select(det_plate, pert_well, cell_line_id, pool_id, pert_mfc_id, pert_type, pert_effective_dose, value) 
}

#' run combat for batch correction
#' @export rev_factor
#' @param x a vector, can be numeric, character or factor
#' @return  a vector of factor with the default order of factor reversed
rev_factor <- function(x){
  if(class(x) != 'factor'){
    x <- as.factor(x)
  }
  
  factor(x, levels = rev(levels(x)))
}
