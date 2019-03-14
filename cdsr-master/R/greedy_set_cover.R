binarize_matrix <- function(mat, threshold){
  
  mat[is.na(mat)] <- 0
  mat[mat > threshold] <- 0
  mat[mat <= threshold] <- 1
  
  return(mat)
  
}

boot_greedy_set_cover <- function(binary_matrix, n_samples=1000){
  
  require(magrittr)
  require(boot)
  require(cowplot)
  
  boot_obj <- boot::boot(binary_matrix, statistic = greedy_set_cover, R = n_samples, stype="i")
  
  boot_df <- boot_obj$t %>%
                magrittr::set_colnames(colnames(binary_matrix)) %>%
                magrittr::set_rownames(paste0("sample", 1:nrow(.))) %>%
                cdsr::melt_matrix(row_name="sample_id", col_name="feature", value_name = "member") %>%
                dplyr::filter(member == 1) %>%
                dplyr::group_by(sample_id) %>%
                dplyr::mutate(set_size=n()) %>%
                dplyr::mutate(set_id=sort(feature %>% as.numeric()) %>% paste(collapse="|")) %>%
                dplyr::ungroup() %>%
                dplyr::group_by(set_id) %>%
                dplyr::mutate(set_frequency=n()/n_samples) %>%
                dplyr::ungroup() %>%
                dplyr::group_by(feature) %>%
                dplyr::mutate(feature_frequency=n()/n_samples) %>%
                dplyr::ungroup()
  
  return(list(boot_obj=boot_obj, 
              boot_df=boot_df))
  
}


greedy_set_cover <- function(binary_matrix, indices){
  
  require(magrittr)
  
  minimal_set <- rep(0, ncol(binary_matrix))
  
  if(any(rowSums(binary_matrix) == 0)){
    warning("No feasible solutions...")
    return(minimal_set)
  }
  
  while(nrow(binary_matrix) > 0){
    max_covering_column <- which(colSums(binary_matrix) == max(colSums(binary_matrix))) %>%
                              magrittr::extract(., sample(1:length(.), 1))
    minimal_set[max_covering_column] <- 1
    binary_matrix <- binary_matrix[ifelse(binary_matrix[, max_covering_column] == 1, F, T), -max_covering_column, drop=F]
  }
  
  return(minimal_set)
  
}