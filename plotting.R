source("~/Documents/MyToolkit/Rtoolkit.R")
library(tidyverse)
library(ggplot2)

DRIVE = loadData("DRIVE")
CN = loadData("cn")

cls = intersect(rownames(DRIVE), rownames(CN))
colors = c(rgb(red = 252,green = 0,blue = 76, maxColorValue = 255), rgb(red=44, green=108, blue=171, maxColorValue = 255))


ggplot(data=DRIVE[cls,], aes(x=reorder(rownames(DRIVE[cls,]), -PRMT5), y=PRMT5, fill=ifelse(CN[cls,"MTAP"] < -1, "Deleted", "WT"))) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual(values = colors, name="MTAP Status") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = c(0.15, 0.1), strip.background = element_blank(), 
        panel.grid = element_blank(), panel.background = element_blank(), plot.title=element_text(face = 'bold', size = 16)
        ) +
  labs(x="Cell Line", y="Dependency Score", title="PRMT5 - shRNA (Project DRIVE, Novartis)") +
  ggsave(filename = "PRMT5_DRIVE_waterfall_plot-SJ-02_23_17.png")

# Redo with CRISPR

CRISPR = loadData("CRISPR")

cls = intersect(rownames(CRISPR), rownames(CN))

ggplot(data=CRISPR[cls,], aes(x=reorder(rownames(CRISPR[cls,]), -PRMT5), y=PRMT5, fill=ifelse(CN[cls,"MTAP"] < -1, "Deleted", "WT")))+
  geom_bar(stat='identity', width=1) +
  scale_fill_manual(values=colors, name="MTAP Status") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = c(0.15, 0.1), strip.background = element_blank(),
        panel.grid = element_blank(), panel.background = element_blank(), plot.title=element_text(face = 'bold', size = 16)
        ) +
  labs(x="Cell Line", y="Dependency Score", title="PRMT5 - CRISPR (Achilles, Broad)") +
  ggsave(filename = "PRMT5_CRISPR_waterfall_plot-SJ-02_23_17.png")
