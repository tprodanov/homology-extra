library(ggplot2)
library(tidyverse)
library(ggrepel)
library(grid)
library(gridExtra)

robustness <- data.frame()
for (name in c('han', 'children', 'IBS', 'CHB')) {
  tmp <- read_delim(
    sprintf('~/Data/hg38/jvc/summaries/pair/%s.csv', name),
    '\t', comment='#', col_types = cols())
  cat(sprintf('%s: %d samples\n', name, length(unique(filter(tmp, type == 'sample')$name))))
  robustness <- rbind(robustness,
      tmp %>% filter(type == 'gene' & refCN != 'all') %>% add_column(dataset = name))
}

robustness$label <- c(
  'han' = '83 Han Chinese samples',
  'children' = '129 IGSR relatives',
  'IBS' = '107 IBS samples',
  'CHB' = '103 CHB samples'
  )[robustness$dataset]
robustness$label <- factor(robustness$label, levels=unique(robustness$label))

plots <- list()
i = 1
for (curr_label in levels(robustness$label)) {
  curr <- filter(robustness, label == curr_label)
  plots[[i]] <- ggplot(curr) +
    geom_point(aes(obs_agCN, agCN_mean), alpha=.5) +
    geom_text_repel(aes(obs_agCN, agCN_mean, label=name), size=2.9,
                    data=filter(curr, agCN_mean < 90),
                    nudge_x = 0, nudge_y = -5, seed = 5, segment.alpha=.3) +
    labs(subtitle = curr_label) +
    scale_x_continuous('Average agCN', breaks=seq(0, 20, 2),
                       limits=c(
                         min(robustness$obs_agCN, na.rm=T),
                         max(robustness$obs_agCN, na.rm=T)), expand=expansion(add=0.5)) +
    scale_y_continuous('Average concordance between runs (%)',
                       limits=c(ifelse(i < 3, 0, 60), 100),
                       breaks=seq(0, 100, ifelse(i < 3, 20, 10))) +
    theme_light() +
    theme(axis.title = element_blank())
  i <- i + 1
}

merged <- do.call(plot_grid, c(plots,
                     list(labels=LETTERS, rel_heights=c(1, 0.6))))
g <- add_axis(merged, 'Average agCN', 'Average concordance between runs (%)')
ggsave('~/Tmp/1.png', g, width=10, height=8, scale=0.6, dpi=450)
