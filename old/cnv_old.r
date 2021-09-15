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

get.mav <- function(values, n) {
  rollapply(values, width=n, FUN=mean, align="right")
}

aggregate.mav <- function(df, ploidies, columns, min_mav=1, max_mav=30) {
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

load <- function(name, ...) {
  filename <- sprintf('%s/%s.csv', data_dir, name)
  if (!file.exists(filename)) {
    cat(sprintf('Cannot load csv file "%s"\n', filename))
    ignore <- file.remove(tmp_file)
    stop(1)
  }
  read.csv(filename, sep='\t', comment.char='#', ...)
}

assign_groups <- function(df, groups) {
  df$group <- groups[match(df$sample, groups$sample), 'path']
  df$shift <- groups[match(df$sample, groups$sample), 'shift']
  n_windows <- sum(df$sample == df$sample[1])
  groups_table <- table(df$group) / n_windows
  groups_table <- ifelse(groups_table >= 5,
                         sprintf('Path %s', names(groups_table)), 'Other')
  df$group <- groups_table[as.character(df$group)]
  df$group <- reorder(df$group, df$group, FUN=function(x) -length(x))
  df
}

if (length(args) == 0) {
  data_dir <- '~/Data/hg38/jvc/runs/024.split/C4A-e'
  plots_prefix <- '~/Tmp/'
} else {
  data_dir <- args[1]
  plots_prefix <- args[2]
}
cat(sprintf('Reading depth from %s/*.csv\n', data_dir))
cat(sprintf('Saving plots to %s.*.png\n', plots_prefix))

windows <- load('windows')
depth <- load('depth')
viterbi <- load('viterbi')
groups <- load('viterbi_paths')
groups <- subset(groups, region == 1)

depth$norm_ploidy <- depth$depth1 / depth$bg_depth1 * 2
depth$norm_ploidy_offset <- depth$norm_ploidy - windows$offset[depth$window_ix + 1]

MAX_VITERBI <- 239
windows <- subset(windows, viterbi_ix <= MAX_VITERBI)
viterbi <- subset(viterbi, viterbi_ix <= MAX_VITERBI)
MAX_WINDOW <- tail(windows, 1)$window_ix
depth <- subset(depth, window_ix <= MAX_WINDOW)

viterbi$iteration <- factor(viterbi$iteration, levels=unique(viterbi$iteration))
row.names(groups) <- groups$sample

depth2 <- subset(depth, !is.na(norm_ploidy_offset))
depth2$norm_ploidy_wo_ploidy <- depth2$norm_ploidy - windows$ploidy[depth2$window_ix + 1]
means <- aggregate(norm_ploidy ~ sample, depth, function(x)
  c('mean' = mean(x), 'var' = var(x)))
means <- do.call(data.frame, means) %>%
  rename(mean=norm_ploidy.mean, var=norm_ploidy.var)
means <- assign_groups(means, groups)
means$group_shift <- sprintf('%s(%+d)', means$group, means$shift)

(g6 <- ggplot(means) +
    geom_point(aes(mean, var), size=1, color='#FF4444') +
    scale_x_continuous('Mean normalized ploidy', breaks=seq(-10, 20, 0.5)) +
    scale_y_continuous('Variation') +
    theme_bw() +
    theme(panel.grid.major.x = element_line(color='gray70')))
ggsave('~/Work/presentations/november/mean_ploidy.png', width=8, height=3)

(g6 <- ggplot(means) +
    geom_point(aes(mean, var, color=group), size=1) +
    scale_x_continuous('Mean normalized ploidy', breaks=seq(-10, 20, 0.5)) +
    scale_y_continuous('Variation') +
    theme_bw() +
    theme(panel.grid.major.x = element_line(color='gray70')))

(g6 <- ggplot(means) +
  geom_point(aes(mean, var, color='A'), size=2) +
  scale_x_continuous('Mean offset', breaks=seq(-10, 20, 0.5)) +
  scale_y_continuous('Variation') +
  scale_color_discrete('Path(shift)') +
  theme_bw())
(g7 <- ggplot(means) +
    geom_histogram(aes(mean, fill=group_shift), binwidth=0.02) +
    scale_x_continuous('Mean offset', breaks=seq(-10, 20, 0.5)) +
    scale_y_continuous('Count') +
    scale_fill_discrete('Path(shift)') +
    theme_bw())
(g8 <- ggplot(means) +
    geom_vline(xintercept=0, color='red', alpha=.9) +
    geom_histogram(aes(mean - round(mean), fill=group_shift), binwidth=0.01) +
    scale_x_continuous('Distance to nearest integer', breaks=seq(-0.5, 0.5, 0.1)) +
    scale_y_continuous('Count') +
    coord_cartesian(xlim=c(-0.5, 0.5)) +
    theme_bw())

plot_grid(g6, g7, g8, ncol=1, rel_heights=c(0.5, 0.25, 0.25))
ggsave(sprintf('%s.4.png', plots_prefix), width=10, height=10, dpi=600)

samples9 <- sprintf('HG00%03d', c(99, 97, 102, 113, 118, 130))
windows9 <- windows[windows$window_ix >= 200 & windows$window_ix <= 386
                & !is.na(windows$viterbi_ix),]
depth9 <- depth[depth$window_ix >= 200 & depth$window_ix <= 386
                & depth$sample %in% samples9 & !is.na(depth$norm_ploidy_offset),]

depth.mav <- aggregate.mav(depth, windows$ploidy, 'norm_ploidy_offset')
windows9$sample <- 'NULL'
ploidy.mav <- aggregate.mav(windows9, windows9$ploidy, 'ploidy', max_mav=15)

depth.mav <- aggregate.mav(depth, rep(4, MAX_WINDOW + 1), 'norm_ploidy')

(g1 <- ggplot(depth.mav) +
    geom_line(aes(window_ix, norm_ploidy, group=sample), alpha=0.1, color='blue') +
    scale_y_continuous('Normalized ploidy', breaks=-10:20) +
    scale_x_continuous('Moving window') +
    theme_bw())
ggsave('~/Work/presentations/november/norm_depth.png', width=8, height=4)

ggplot(subset(viterbi, iteration == 'rev')) +
  geom_line(aes(viterbi_ix, state, group=sample), alpha=.1, color='blue')

(g1 <- ggplot(depth.mav) +
    geom_line(aes(window_ix, norm_ploidy_offset, group=sample, color=sample), alpha=1) +
    scale_y_continuous('Normalized ploidy', breaks=-10:20) +
    scale_x_continuous('Moving window') +
    theme_bw())

(g1 <- ggplot(depth.mav) +
    geom_line(aes(window_ix, norm_ploidy, group=sample), color='blue', alpha=1) +
    scale_y_continuous('Normalized ploidy', breaks=-10:20) +
    scale_x_continuous('Moving window') +
    theme_bw())

(g1 <- ggplot(depth.mav) +
  geom_line(aes(window_ix, ploidy), data=ploidy.mav, size=2, color='black') +
  geom_line(aes(window_ix, norm_ploidy, group=sample), color='blue', alpha=.1) +
  scale_y_continuous('Normalized ploidy', breaks=-10:20) +
  scale_x_continuous('Moving window') +
  ggtitle('(a) Normalized ploidy (all windows)') +
  theme_bw())

windows2 <- subset(windows, !is.na(viterbi_ix))
depth2.mav <- aggregate.mav(depth2, windows2$ploidy,
  c('norm_ploidy', 'norm_ploidy_offset')) #%>% assign_groups(groups)
windows2$sample <- 'NULL'
ploidy2.mav <- aggregate.mav(windows2, windows2$ploidy, 'ploidy')

(g2 <- ggplot(depth2.mav) +
  geom_line(aes(window_ix, ploidy), data=ploidy2.mav, size=2, color='black') +
  geom_line(aes(window_ix, norm_ploidy, group=sample), color='blue', alpha=.1) +
  scale_y_continuous('Normalized ploidy', breaks=-10:20) +
  scale_x_continuous('Moving window') +
  ggtitle('(b) Normalized ploidy (good windows)') +
  theme_bw())

(g3 <- ggplot(depth2.mav) +
    geom_line(aes(window_ix, norm_ploidy_offset, group=sample), color='blue', alpha=.1) +
    scale_y_continuous('Shifted normalized ploidy', breaks=-10:20) +
    scale_x_continuous('Moving window') +
    theme_bw())
ggsave('~/Work/presentations/november/shifted_norm_depth.png', width=8, height=4)

viterbi_x <- subset(viterbi, iteration == '1a')
viterbi_x <- subset(viterbi, iteration == 'rev')
viterbi_matrix <- pivot_wider(
  viterbi_x, names_from = 'sample', values_from = 'state')[, c(-1, -2)] %>%
  as.matrix
sample_dist <- dist(t(viterbi_matrix + 4))
sample_clust <- hclust(sample_dist)
viterbi_x$sample <- factor(viterbi_x$sample, levels=sample_clust$labels[sample_clust$order])

ggplot(subset(viterbi_x, viterbi_ix >= 200)) +
  geom_tile(aes(x=viterbi_ix, y=sample, fill=factor(state + 4))) +
  scale_x_continuous('100bp window') +
  scale_y_discrete('Samples', labels=NULL, breaks=NULL) +
  scale_fill_viridis('Pooled\nploidy', discrete=T) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank())
ggsave('~/Work/presentations/november/viterbi_1.png', width=16, height=8, scale=.515)
ggsave('~/Work/presentations/november/viterbi_1_zoom.png', width=16, height=8, scale=.515)
ggsave('~/Work/presentations/november/viterbi_rev_zoom.png', width=16, height=8, scale=.515)

ggplot(viterbi_x) +
  geom_line(aes(x=viterbi_ix, y=state + 4, group=sample), color='blue', alpha=.2) +
  scale_x_continuous('100bp window') +
  scale_y_continuous('Pooled ploidy') +
  theme_bw()
ggsave('~/Work/presentations/november/viterbi_1b.png', width=8, height=4)
ggsave('~/Work/presentations/november/viterbi_revb.png', width=8, height=4)


(g3 <- ggplot(depth2.mav) +
  geom_line(aes(window_ix, ploidy), data=ploidy2.mav, size=2, color='black') +
  geom_line(aes(window_ix, norm_ploidy_offset, group=sample), color='blue', alpha=.1) +
  scale_y_continuous('Shifted normalized ploidy', breaks=-10:20) +
  scale_x_continuous('Moving window') +
  ggtitle('(c) Shifted normalized ploidy (good windows)') +
  theme_bw())
plot_grid(g1, g2, g3, ncol=1)
ggsave(sprintf('%s.1.png', plots_prefix), width=10, height=8, dpi=600)

(g4 <- ggplot(depth2.mav) +
  geom_line(aes(window_ix, ploidy), data=ploidy2.mav, size=.5, color='black') +
  geom_line(aes(window_ix, norm_ploidy_offset, group=sample,
                color=sprintf('%+d', shift)), alpha=.5) +
  scale_y_continuous('Normalized ploidy - offset', breaks=-10:20) +
  scale_x_continuous('Moving window') +
  facet_wrap(~ group, ncol=1) +
  scale_color_discrete('Shift') +
  theme_bw() +
  theme(strip.text.x=element_text(margin=margin(b=0, t=0))))
ggsave(sprintf('%s.2.png', plots_prefix), g4, width=10, height=8, dpi=600)

(g5 <- ggplot(viterbi) +
  geom_line(aes(viterbi_ix, state, group=sample), color='blue', alpha=.5) +
  facet_wrap(~ iteration, ncol=3) +
  scale_x_continuous('Window') +
  scale_y_continuous('Hidden state') +
  theme_bw() +
  theme(strip.text.x=element_text(margin=margin(b=0, t=0))))
ggsave(sprintf('%s.3.png', plots_prefix), g5, width=10, height=10, dpi=600)
cat(sprintf('Success [%s]\n', data_dir))

# VITERBI_IX <- 0
# WINDOW_IX <- subset(windows, viterbi_ix == VITERBI_IX)$window_ix
# depth_w <- subset(depth, window_ix == WINDOW_IX)
# levels(viterbi$iteration)
# viterbi_w <- subset(viterbi, viterbi_ix == VITERBI_IX & iteration == '9a')
# samples_state <- viterbi_w$state
# names(samples_state) <- as.character(viterbi_w$sample)
# depth_w$state <- samples_state[depth_w$sample]
# 
# ggplot(depth_w) +
#   geom_histogram(aes(norm_ploidy, fill=factor(state)),
#                  position='identity', binwidth=0.1, alpha=.4) +
#   geom_vline(xintercept=subset(windows, window_ix == WINDOW_IX)$shift + 1:6) +
#   scale_x_continuous(breaks=0:10)
# 
# ggsave('~/Downloads/2.png', width=10, height=15)
