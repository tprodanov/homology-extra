#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pdf(NULL)

suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(zoo))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(viridis))

if (length(args) == 0) {
  data_dir <- '~/Data/hg38/jvc/runs/024.split/C4A-a/extra'
  plots_prefix <- '~/Tmp/'
  important_samples <- c('HG00136', 'HG00105', 'HG00717', 'HG00378', 'HG00740',
    'HG00110', 'HG01334', 'HG00336', 'HG00380', 'HG01260', 'HG00620',
    'HG00620', 'HG01174', 'HG01495')
  important_samples <- c()
} else {
  data_dir <- trimws(args[1], which='right', whitespace='/')
  plots_dir <- trimws(args[2], which='right', whitespace='/')
  if (length(args) > 2) {
    important_samples <- args[3:length(args)]
  } else {
    important_samples <- c()
  }
  
  dir_split <- unlist(strsplit(data_dir, '/', fixed=T))
  dir_split <- dir_split[dir_split != 'extra']
  start <- match('.', dir_split)
  start <- ifelse(is.na(start), 1, start + 1)
  plots_prefix <- paste(dir_split[start:length(dir_split)], collapse='.')
  plots_prefix <- sprintf('%s/%s.', plots_dir, plots_prefix)
}
has_important <- length(important_samples) > 0
cat(sprintf('Reading depth from %s/*.csv\n', data_dir))
cat(sprintf('Saving plots to %s*.png\n', plots_prefix))

get.mav <- function(values, n) {
  rollapply(values, width=n, FUN=mean, align="right")
}

aggregate.mav <- function(df, ploidies, columns, min_mav=1, max_mav=30) {
  if (!('sample' %in% colnames(df))) {
    df$sample <- 'NULL'
  }
  n_windows <- length(ploidies)
  n_samples <- length(unique(df$sample))
  
  if (n_windows * n_samples != nrow(df)) {
    cat(sprintf('Ploidies: n_windows = %d\n', n_windows))
    cat(sprintf('Df:       n_samples = %d\n', n_samples))
    cat(sprintf('n_windows * n_samples = %d\n', n_windows * n_samples))
    cat(sprintf('Dataframe length      = %d\n', nrow(df)))
    stop(1)
  }

  const_ends <- cumsum(rle(ploidies)$lengths)
  const_ranges <- data.frame(start=c(0, head(const_ends, -1)) + 1, end=const_ends)
  const_ranges$mva <- pmax(
    min_mav,
    pmin(max_mav, floor((const_ranges$end - const_ranges$start) / 2)))

  if (nrow(const_ranges) > 1) {
    regions.mav <- function(values) {
      do.call(c, apply(const_ranges, 1,
                       function(row) get.mav(values[row[1]:row[2]], row[3])))
    }
  } else {
    regions.mav <- function(values) get.mav(values, const_ranges[1, 3])
  }

  new_n_windows <- length(regions.mav(1:n_windows))
  res <- data.frame(window_ix=rep(1:new_n_windows, n_samples))
  for (col in columns) {
    if (any(is.na(df[col]))) {
      cat(sprintf('Column "%s" contains NAs\n', col))
      stop(1)
    }
    tmp <- aggregate(as.formula(sprintf('%s ~ sample', col)), df, regions.mav)
    tmp <- cbind(data.frame(sample=tmp$sample), as.data.frame(tmp[[col]]))
    tmp <- pivot_longer(tmp, starts_with('V'), names_prefix='V',
                         names_to='window', values_to=col)
    
    if (ncol(res) == 1) {
      res$sample <- tmp$sample
    }
    res <- cbind(res, tmp[col])
  }
  res
}

load <- function(name, comment.char='#', ...) {
  filename <- sprintf('%s/%s%s', data_dir, name, ifelse(grepl('\\.', name), '', '.csv'))
  if (!file.exists(filename)) {
    cat(sprintf('Cannot load csv file "%s"\n', filename))
    ignore <- file.remove(tmp_file)
    stop(1)
  }
  read.csv(filename, sep='\t', comment.char=comment.char, ...)
}

windows <- load('windows.bed', comment.char='')
depth <- load('depth')
viterbi <- load('viterbi_states')
transitions <- load('hmm_params')

# Transform input data frames

transitions <- group_by(transitions, window_ix) %>% slice_tail(n = 1)
windows$multiplier <- transitions[match(windows$window_ix,
                                        transitions$window_ix),]$multiplier
depth$norm_ploidy <- depth$depth1 / depth$bg_depth1 * 2

windows_v <- filter(windows, !is.na(multiplier))
windows_v$viterbi_ix <- 1:nrow(windows_v) - 1

viterbi$iteration <- factor(viterbi$iteration,
  levels=c(sprintf('%s%s', rep(1:100, 2), rep(c('a', 'b'), each=100)), 'final'))
viterbi$viterbi_ix <- windows_v[match(viterbi$window_ix,
    windows_v$window_ix),]$viterbi_ix
viterbi <- pivot_longer(viterbi, starts_with('HG'),
                        names_to='sample', values_to='pred_ploidy')
viterbi_last <- group_by(viterbi, sample, window_ix) %>% slice_tail(n = 1)

depth$norm_ploidy_offset <- depth$norm_ploidy / windows[depth$window_ix + 1,]$multiplier
depth_v <- filter(depth, !is.na(norm_ploidy_offset))
depth_v <- left_join(depth_v, viterbi_last, by=c('window_ix', 'sample'))

# Moving average (without removing offsets)

depth.mav <- aggregate.mav(depth, windows$ploidy, 'norm_ploidy')
ploidy.mav <- aggregate.mav(windows, windows$ploidy, 'ploidy')
(g_all <- ggplot(depth.mav, aes(window_ix, norm_ploidy)) +
    geom_line(aes(y=ploidy), data=ploidy.mav, color='black', size=1) +
    geom_line(aes(group=sample), alpha=0.1,
              color=ifelse(has_important, 'gray30', 'blue')) +
    geom_line(aes(color=sample),
              data=subset(depth.mav, sample %in% important_samples)) +
    scale_y_continuous('Normalized ploidy', breaks=-10:20) +
    scale_x_continuous('Moving window') +
    scale_color_discrete('Sample') +
    ggtitle('(a) Normalized ploidy (all windows)') +
    theme_bw())
ggsave(sprintf('%snorm_ploidy_a.png', plots_prefix), width=8, height=5)

# Moving average (only good windows)

depth_v.mav <- aggregate.mav(depth_v, windows_v$ploidy,
                             c('norm_ploidy', 'norm_ploidy_offset', 'pred_ploidy'))
ploidy_v.mav <- aggregate.mav(windows_v, windows_v$ploidy, 'ploidy')
(g_v <- ggplot(depth_v.mav, aes(window_ix, norm_ploidy)) +
  geom_line(aes(y=ploidy), data=ploidy_v.mav, size=1, color='black') +
  geom_line(aes(group=sample), alpha=.1,
            color=ifelse(has_important, 'gray30', 'blue')) +
    geom_line(aes(color=sample),
              data=subset(depth_v.mav, sample %in% important_samples)) +
  scale_y_continuous('Normalized ploidy', breaks=-10:20) +
  scale_x_continuous('Moving window') +
  ggtitle('(b) Normalized ploidy (good windows)') +
  theme_bw())
ggsave(sprintf('%snorm_ploidy_b.png', plots_prefix), width=8, height=5)

# Moving average (only good windows, removing offsets)

(g_v_offset <- ggplot(depth_v.mav, aes(window_ix, norm_ploidy_offset)) +
    geom_line(aes(y=ploidy), data=ploidy_v.mav, size=1, color='black') +
    geom_line(aes(group=sample), alpha=.1,
              color=ifelse(has_important, 'gray30', 'blue')) +
    geom_line(aes(color=sample),
              data=subset(depth_v.mav, sample %in% important_samples)) +
    scale_y_continuous('Normalized ploidy', breaks=-10:20) +
    scale_x_continuous('Moving window') +
    ggtitle('(c) Normalized ploidy (good windows, remove offset)') +
    theme_bw())
ggsave(sprintf('%snorm_ploidy_c.png', plots_prefix), width=8, height=5)

# Norm ploidy, colored by predicted ploidy

depth_v.mav <- mutate(depth_v.mav,
    pred_ploidy_int = ifelse(pred_ploidy == round(pred_ploidy), pred_ploidy, NA))
(g_v_offset <- ggplot(depth_v.mav, aes(window_ix, norm_ploidy_offset)) +
    geom_line(aes(y=ploidy), data=ploidy_v.mav, size=1, color='black') +
    geom_line(aes(group=sample, color=factor(pred_ploidy_int)), alpha=.6) +
    scale_y_continuous('Normalized ploidy', breaks=-10:20) +
    scale_x_continuous('Moving window') +
    ggtitle('(d) Normalized ploidy (good windows, remove offset)') +
    guides(color = guide_legend(override.aes = list(alpha=1, size=2))) +
    scale_color_viridis('Predicted\nploidy', discrete=T, na.value='gray60') +
    theme_bw())
ggsave(sprintf('%snorm_ploidy_d.png', plots_prefix), width=8, height=5)

# Viterbi paths (geom_tile)

viterbi_matrix <- select(viterbi_last, sample, window_ix, pred_ploidy) %>%
  pivot_wider(names_from = 'sample', values_from = 'pred_ploidy') %>%
  ungroup %>%
  select(!window_ix) %>%
  as.matrix
sample_dist <- dist(t(viterbi_matrix))
sample_clust <- hclust(sample_dist)

viterbi_first_last <- rbind(
  filter(viterbi, iteration == '1a'),
  viterbi_last %>% mutate(iteration = 'final'))
viterbi_first_last$sample2 <- factor(viterbi_first_last$sample,
                                     levels=sample_clust$labels[sample_clust$order])

ggplot(viterbi_first_last) +
  geom_tile(aes(x=viterbi_ix, y=sample2, fill=factor(pred_ploidy))) +
  scale_x_continuous('100bp window') +
  scale_y_discrete('Samples', labels=NULL, breaks=NULL) +
  scale_fill_viridis('Predicted\nploidy', discrete=T) +
  facet_wrap(~ iteration, ncol=1) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank())
ggsave(sprintf('%sviterbi_tile.png', plots_prefix), width=8, height=8)

# Viterbi paths (geom_line)

transitions <- unlist(aggregate(pred_ploidy ~ sample,
    viterbi_last, function(x) which(diff(x) != 0))$pred_ploidy)
n_samples <- length(unique(viterbi$sample))
transitions <- table(transitions)
transitions <- as.numeric(names(transitions)[transitions >= 0.1 * n_samples])
max_viterbi_ix <- max(windows_v$viterbi_ix)
text_pos <- c(0, round((head(transitions, -1) + tail(transitions, -1)) / 2),
                 max_viterbi_ix)
text_pos <- unique(round(text_pos / 3) * 3)
text_pos <- sapply(text_pos, function(x) min(x, max_viterbi_ix))

viterbi_counts <- filter(viterbi_first_last, viterbi_ix %in% text_pos) %>%
  aggregate(sample ~ iteration + viterbi_ix + pred_ploidy, ., length)

ggplot(viterbi_first_last, aes(viterbi_ix, pred_ploidy)) +
  geom_line(aes(y=ploidy), data=windows_v, size=5, color='gray70') +
  geom_line(aes(group=sample), alpha=.2, color='blue') +
  geom_text(aes(y = pred_ploidy + 0.3, label = sample), data=viterbi_counts, size=3) +
  scale_x_continuous('100bp window') +
  scale_y_continuous('Predicted ploidy', breaks=0:20, minor_breaks=NULL) +
  scale_color_discrete('Sample') +
  facet_wrap(~ iteration, ncol=1) +
  theme_bw()
ggsave(sprintf('%sviterbi_line.png', plots_prefix), width=8, height=5)

if (has_important) {
  ggplot(viterbi_first_last, aes(viterbi_ix, pred_ploidy)) +
    geom_line(aes(y=ploidy), data=windows_v, size=5, color='gray70') +
    geom_line(aes(group=sample), alpha=0.1, color='gray30') +
    geom_line(aes(color=sample), size=2,
              data=subset(viterbi_first_last, sample %in% important_samples)) +
    scale_x_continuous('100bp window') +
    scale_y_continuous('Predicted ploidy', breaks=0:20, minor_breaks=NULL) +
    scale_color_discrete('Sample') +
    facet_wrap(~ iteration, ncol=1) +
    theme_bw()
  ggsave(sprintf('%sviterbi_important.png', plots_prefix), width=8, height=5)
}

n_important <- length(important_samples)
if (has_important && n_important < 40) {
  ggplot(filter(depth_v, sample %in% important_samples)) +
    geom_vline(xintercept=transitions, color='gray80', linetype='dashed') +
    geom_line(aes(viterbi_ix, ploidy), data=windows_v, size=2, color='gray80') +
    geom_line(aes(viterbi_ix, norm_ploidy_offset, group=sample,
                  color = factor(pred_ploidy))) +
    facet_wrap(~ sample, ncol=3, scales='free_y') +
    scale_x_continuous('100bp window') +
    scale_y_continuous('Normalized ploidy', breaks=0:20, minor_breaks=NULL) +
    scale_color_discrete('Predicted\nploidy') +
    guides(color = guide_legend(override.aes = list(alpha=1, size=2))) +
    theme_bw() +
    theme(panel.grid.major.y=element_line(color='gray80'))
  
  width <- 5 + max(depth_v$viterbi_ix) / 100 * 1.5 * min(n_important, 3)
  height <- max(5, ceiling(length(important_samples) / 3) * 2)
  ggsave(sprintf('%snorm_ploidy_important.png', plots_prefix),
         width=width, height=height, limitsize=F)
}
