# root_qtl.R

library(ggplot2)
library(grid)

setwd("~/Desktop/Lab/Maloof Lab/QTL_Visualization/")

qtl <- read.csv("chr_pos_lod_name.csv")

qtl$background.color <- "1"
qtl$background.color[(qtl$chr %% 2) == 0] <- "0"

ggplot(qtl) +
  facet_grid(~ chr, scales = "free_x", space = "free_x") + 
  geom_rect(aes(fill = background.color, xmin = -Inf, xmax = Inf, ymin = -0.02, ymax = Inf)) +
  geom_line(aes(pos, lod), size = 1) +
  geom_hline(yintercept = 0.50, color = "red", size = 1) +
  geom_segment(aes(x = pos, xend = pos, y = -0.005, yend = -0.015)) +
  scale_x_continuous(breaks = seq(0, 1000, by = 50)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.02, max(qtl$lod * 1.02))) +
  scale_fill_manual(values = c("0" = "gray90", "1" = "white")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90),
        axis.line=element_line(),
        panel.margin = unit(0, "cm")) +
  ggtitle("LOD Curves for QTLs") +
  xlab("Genetic Position in cM") +
  ylab("LOD Score")


# if num of chromosome is odd, 1 == gray90. otherwise, 1 == white