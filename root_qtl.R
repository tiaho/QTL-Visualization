# root_qtl.R

# pick a threshold for the log scores
library(ggplot2)

setwd("~/Desktop/Lab/Maloof Lab/QTL_Visualization/")

qtl <- read.csv("chr_pos_lod_name.csv")

ggplot(qtl) +
  geom_line(aes(pos, lod), size = 1.5) +
  facet_grid(~ chr, scales = "free_x", space = "free_x") + 
  geom_rect(aes(fill = as.factor(chr)),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.01) +
  geom_hline(yintercept = 0.50, color = "red", size = 1) +
  theme(legend.position = "none") +
  geom_histogram(aes(pos), binwidth = 0.1)
