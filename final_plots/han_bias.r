library(tidyverse)
library(ggplot2)
library(gridExtra)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

plots <- list()
dir <- '203.han_g1k'; name <- 'IGSR'
dir <- '204.han_bgi'; name <- 'BGI'

args <- parser$parse_args(
  c('-i', sprintf('~/Data/hg38/jvc/runs/%s/NCF1/v100/extra', dir), '-o', '/tmp/',
    '-w', '25', '-r'))

(plots[[paste0(name, '_a')]] <- ggplot(depth_v.av) +
    geom_line(aes(middle, copy_num), size=2, data=windows_v.av, color='gray60') +
    geom_line(aes(middle, norm_cn1, group=sample), alpha=0.4, size=0.2, color='blue') +
    geom_line(aes(middle, copy_num), size=2, data=windows_v.av, alpha=.2) +
    scale_x_continuous(x_label, labels=fmt_pos) +
    scale_y_continuous('Normalized read depth', breaks=0:40,
                       labels=function(y) ifelse(y < 8, as.character(y), '')) +
    coord_cartesian(ylim=c(3.1, 7.9)) +
    theme_light() +
    theme(axis.title = element_blank()))

(plots[[paste0(name, '_b')]] <- ggplot(depth_v.av) +
    geom_line(aes(middle, copy_num), size=2, data=windows_v.av, color='gray60') +
    geom_line(aes(middle, corr_cn1, group=sample), alpha=0.4, size=0.2, color='blue') +
    geom_line(aes(middle, copy_num), size=2, data=windows_v.av, alpha=.2) +
    scale_x_continuous(x_label, labels=fmt_pos) +
    scale_y_continuous('Corrected norm. read depth', breaks=0:40,
                       labels=function(y) ifelse(y < 8, as.character(y), '')) +
    theme_light() +
    theme(axis.title = element_blank()))

# depth_g1k <- depth_v.av
# depth_bgi <- depth_v.av

g <- plot_grid(plots[['IGSR_a']], plots[['IGSR_b']],
          plots[['BGI_a']], plots[['BGI_b']], labels=LETTERS,
          label_x = -0.02)
g2 <- add_axis(g, 'Position (chr7)', 'Normalized read depth')
ggsave('~/Tmp/1.png', g2, width=10, height=8, scale=.7, dpi=450)

min(windows$start)
max(windows$end)
