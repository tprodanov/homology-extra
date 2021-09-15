library(tidyverse)
library(ggplot2)
library(cowplot)
library(gridExtra)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

full <- data.frame()
dir <- '~/Data/hg38/jvc/summaries/pair/'
for (filename in list.files(dir)) {
  if (grepl('merged', filename)) {
    next;
  }
  tmp <- read_delim(file.path(dir, filename), '\t', comment='#',
                    col_types = cols())
  tmp$dataset <- sub('.csv', '', filename)
  full <- rbind(full, tmp)
}
rm(tmp)

genes <- filter(full, type == 'gene' & refCN != 'all')

calc_cumulative <- function(df, colname) {
  df <- df[!is.na(df[[colname]]),]
  dataset <- df$dataset
  col <- df[[colname]]
  
  # thresholds <- c(0, round(col, 2))
  # thresholds <- sort(unique(thresholds))
  thresholds <- seq(0, 100, 0.5)

  res <- data.frame()
  for (threshold in thresholds) {
    res <- rbind(res,
        data.frame(dataset = dataset, is_concordant = col >= threshold) %>%
          count(dataset, is_concordant) %>% add_column(threshold = threshold))
  }
  max_values <- aggregate(n ~ dataset, res, max)
  max_values2 <- max_values$n
  names(max_values2) <- max_values$dataset
  print(max_values2)
  
  res$perc <- 100 * res$n / max_values2[res$dataset]
  pivot_wider(res, names_from='is_concordant', values_from=c('n', 'perc'), values_fill=0) %>%
    rename(conc = n_TRUE, disc = n_FALSE, conc_perc = perc_TRUE, disc_perc = perc_FALSE)
}

plots <- list()

# Han samples ====

genes_han <- filter(genes, dataset == 'han')
counts_han <- rbind(
  calc_cumulative(genes_han, 'agCN_mean') %>%
    add_column(label='agCN: All loci'),
  calc_cumulative(filter(genes_han, obs_agCN < 5), 'agCN_mean') %>%
    add_column(label='agCN: Mean agCN < 5'),
  calc_cumulative(filter(genes_han, obs_agCN < 7), 'agCN_mean') %>%
    add_column(label='agCN: Mean agCN < 7'),
  calc_cumulative(genes_han, 'psCN_mean') %>%
    add_column(label='psCN: All loci')
)

# Two model parameters ====

curr_genes <- filter(genes, dataset %in% c('children', 'IBS', 'CHB'))

counts_2m <- calc_cumulative(curr_genes, 'agCN_mean')
full_names <- c(
  'children' = '129 EUR relatives',
  'CHB' = '103 Han Chinese samples',
  'IBS' = '107 Iberian samples'
)
counts_2m$label <- factor(full_names[counts_2m$dataset], levels=full_names)

# Subsampling ====

curr_genes <- filter(genes, dataset %in% c('IBS', 'IBS10x', 'IBS20x'))

counts_sub <- calc_cumulative(curr_genes, 'agCN_mean')
full_names <- c(
  'IBS' = 'Original: 30x',
  'IBS20x' = 'Subsampled to 20x',
  'IBS10x' = 'Subsampled to 10x'
)
counts_sub$label <- factor(full_names[counts_sub$dataset], levels=full_names)

# Plots ====

(plots[['han']] <- ggplot(counts_han) +
  geom_line(aes(threshold, conc_perc, color=label), size=1) +
  scale_x_continuous('Concordance between runs',
                     breaks=seq(0, 100, 5), minor_breaks=0:100,
                     expand=expansion(add=.5)) +
  scale_y_continuous('Percentage of loci') +
  coord_cartesian(xlim = c(80, 100)) +
  scale_color_discrete('') +
  theme_light() +
  ggtitle('83 Han Chinese samples: independent IGSR & BGI datasets') +
  theme(
    axis.title = element_blank(),
    legend.position = c(0.04, 0.01),
    legend.justification = c('left', 'bottom'),
    legend.background = element_rect(fill='transparent'),
    legend.key = element_rect(fill='transparent'),
    plot.title = element_text(size=11)))

(plots[['2m']] <- ggplot(counts_2m) +
    geom_line(aes(threshold, conc_perc, color=label), size=1) +
    scale_x_continuous('Concordance between runs',
                       breaks=seq(0, 100, 5), minor_breaks=0:100,
                       expand=expansion(add=.5)) +
    scale_y_continuous('Percentage of loci') +
    coord_cartesian(xlim = c(80, 100), ylim=c(90, 100)) +
    scale_color_discrete('') +
    theme_light() +
    ggtitle('Parascopy analysis with two independent sets of model parameters') +
    theme(
      axis.title = element_blank(),
      legend.position = c(0.04, 0.1),
      legend.justification = c('left', 'bottom'),
      legend.background = element_rect(fill='transparent'),
      legend.key = element_rect(fill='transparent'),
      plot.title = element_text(size=11)))

(plots[['subs']] <- ggplot(counts_sub) +
    geom_line(aes(threshold, conc_perc, color=label), size=1) +
    scale_x_continuous('Concordance between runs',
                       breaks=seq(0, 100, 5), minor_breaks=0:100,
                       expand=expansion(add=.5)) +
    scale_y_continuous('Percentage of loci') +
    coord_cartesian(xlim = c(80, 100), ylim=c(90, 100)) +
    scale_color_discrete('') +
    theme_light() +
    ggtitle('107 Iberian samples: two independent sets of model parameters') +
    theme(
      axis.title = element_blank(),
      legend.position = c(0.04, 0.1),
      legend.justification = c('left', 'bottom'),
      legend.background = element_rect(fill='transparent'),
      legend.key = element_rect(fill='transparent'),
      plot.title = element_text(size=11)))

(g <- add_axis(
  plot_grid(plots[['han']],
          plots[['2m']],
          plots[['subs']],
          ncol=1, labels=LETTERS),
  'Concordance threshold', 'Percentage of loci'))
ggsave('~/Tmp/1.png', g, width=10, height=10, scale=.6, dpi=450)

ggplot(counts_han) +
    geom_line(aes(threshold, conc_perc, color=label), size=1) +
    scale_x_continuous('Concordance between runs',
                       breaks=seq(0, 100, 5), minor_breaks=0:100,
                       expand=expansion(add=.5)) +
    scale_y_continuous('Percentage of loci') +
    coord_cartesian(xlim = c(80, 100)) +
    scale_color_discrete('') +
    theme_light() +
    # ggtitle('83 Han Chinese samples: independent IGSR & BGI datasets') +
    theme(
      legend.position = c(0.04, 0.01),
      legend.justification = c('left', 'bottom'),
      legend.background = element_rect(fill='transparent'),
      legend.key = element_rect(fill='transparent'),
      plot.title = element_text(size=11))
ggsave('~/Tmp/1.png', width=6, height=2.5, scale=.8, dpi=600)

han <- read_delim('~/Data/hg38/jvc/summaries/pair/han.csv', '\t', comment = '#')
han <- filter(han, type == 'gene')
han <- group_by(han, name) %>% slice_tail(n = 1) %>% ungroup()
filter(han, is.na(agCN_mean))

filter(han, agCN_mean >= 99) %>% nrow
filter(han, agCN_mean < 90) %>% nrow
filter(han, agCN_mean < 90)
filter(han, agCN_mean < 90 & obs_agCN > 7) %>% nrow


ibs <- read_delim('~/Data/hg38/jvc/summaries/pair/IBS.csv', '\t', comment = '#') %>%
  filter(type == 'gene') %>% group_by(name) %>% slice_tail(n = 1) %>% ungroup()
chb <- read_delim('~/Data/hg38/jvc/summaries/pair/CHB.csv', '\t', comment = '#') %>%
  filter(type == 'gene') %>% group_by(name) %>% slice_tail(n = 1) %>% ungroup()

filter(ibs, agCN_mean >= 99) %>% nrow
filter(chb, agCN_mean >= 99) %>% nrow
filter(ibs, agCN_mean < 99)
filter(chb, agCN_mean < 99)

filter(ibs, is.na(agCN_mean))
filter(chb, is.na(agCN_mean))
