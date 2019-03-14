
Rsquared_rm_na <- function(col){
  ifelse(is.na(col), 0, col) %>% round(4)
}

#' parsing an ATLANTIS summary table
#' 
#' @importFrom magrittr "%>%"
#' @param df The ATLANTIS summary table download from Tundra 
#' @param match.confounder logic of whether to match the confounder-only model for each targetID
#' @param cutoff a numeric value between -1 and 1. When match.confounder is TRUE, the corrected Rsquared will be calculated and value below or equal to this cutoff will be filtered out
#' @return if match.confounder is TRUE, A data frame in which each row contains the ATLANTIS model metrics of Rsquared and topFeatureas well as the model with confounder only as feature. If match.confounder is FALSE, no confounder models will be matched and the output also include a column of report url 
#' @author Li Wang
#'
#' @details  For generating the main table in the project dashboard, use match.confounder = T, cutoff = 0. For generating the table in the Biomarker tab in the compound dashboard, use match.confounder = F

#' @export parse_atlantis_table

parse_atlantis_table <- function(df, match.confounder = T, cutoff = 0){
  
  df %<>% dplyr::select(targetID, predictiveFeatures, Rsquared, topFeature, report) %>%
    dplyr::mutate(targetID = gsub('^GS_', '', targetID), 
           predictiveFeatures = gsub('\', \'', '_', predictiveFeatures),
           predictiveFeatures = gsub('\\[\'|\'\\]|PCAL_|PR500_|BAYER_REP_', '', predictiveFeatures),
           Rsquared = Rsquared_rm_na(Rsquared))
  if(match.confounder){
    stopifnot(cutoff <= 1 & cutoff >= -1)
    df %>%
      dplyr::left_join(filter(., predictiveFeatures == 'Confounders') %>%
                  rename(Rsquared_CONF = Rsquared, topFeature_CONF = topFeature) %>%
                  mutate(Rsquared_CONF = Rsquared_rm_na(Rsquared_CONF)) %>%
                  select(-predictiveFeatures, -report)) %>%
      dplyr::select(-report) %>%
      dplyr::mutate(Rsquared_corrected = Rsquared - Rsquared_CONF) %>%
      dplyr::filter(Rsquared_corrected > cutoff)
  }else{
    df
  }
}

