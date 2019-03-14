library(RSellersLab)
library(dplyr)
library(tidyverse)
library(e1071)

combined = load_data("Combined", gene="WSB2")

tissues = sapply(rownames(combined), get_tissue)

df = data.frame(
  rnai = combined,
  tissue = tissues
)

df <- df %>%
  group_by(tissue) %>%
  mutate(skew_by_tissue = skewness(WSB2, na.rm=T)) %>%
  mutate(n_cells = length(WSB2))
  ungroup()

ggplot(data=subset(df,n_cells > 1), aes(x=WSB2, y=reorder(tissue, -skew_by_tissue))) + 
  geom_point(color="orange", size=3, alpha=0.5) +
  theme_minimal() +
  theme(line=element_blank())
