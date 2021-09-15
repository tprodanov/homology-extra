library(ggplot2)
library(tidyverse)

depth <- read.csv('~/Data/hg38/jvc/runs/201.han_g1k/ANAPC1/extra/depth.csv', sep='\t')
all_samples <- unique(depth$sample)

data_dir <- '~/Data/hg38/jvc/runs'
datasets <- c('201.han_g1k' = '1000 genomes',  '202.han_hc' = 'Han high cov')
region <- 'ANAPC1'
rundir <- ''
samples <- c('HG00428')
# samples <- all_samples

data <- data.frame()
for (i in seq_along(datasets)) {
  curr_dir <- sprintf('%s/%s/%s/%s/extra', data_dir, names(datasets)[i], region, rundir)
  depth <- read.csv(sprintf('%s/depth.csv', curr_dir), sep='\t')
  depth <- filter(depth, sample %in% samples)
  
  windows <- read.csv(sprintf('%s/windows.bed', curr_dir), sep='\t')
  colnames(windows)[1] <- 'chrom'
  hmm_params <- read.csv(sprintf('%s/hmm_params.csv', curr_dir), sep='\t')
  hmm_params <- group_by(hmm_params, window_ix) %>% slice_tail(n = 1)
  windows$multiplier <- hmm_params[match(windows$window_ix,
                                   hmm_params$window_ix),]$multiplier
  
  viterbi <- read.csv(sprintf('%s/viterbi_states.csv', curr_dir), sep='\t')
  viterbi <- select(viterbi, 'window_ix' | all_of(samples))
  viterbi <- pivot_longer(viterbi, all_of(samples),
                          names_to='sample', values_to='pred_ploidy')
  viterbi_last <- group_by(viterbi, sample, window_ix) %>% slice_tail(n = 1)

  joined <- left_join(depth, windows, by='window_ix')
  joined <- left_join(joined, viterbi_last, by=c('window_ix', 'sample'))
  joined$dataset <- datasets[i]
  data <- rbind(data, joined)
}

data <- data %>%
  mutate(norm_ploidy = depth1 / bg_depth1 * 2) %>%
  mutate(norm_ploidy_offset = norm_ploidy / multiplier)

ggplot(data) +
  geom_line(aes(window_ix, ploidy), color='gray70', size=3) +
  geom_line(aes(window_ix, norm_ploidy, color=dataset)) +
  facet_wrap(~ sample) +
  theme_bw() +
  scale_y_continuous(breaks=0:20)

ggplot(filter(data, in_viterbi)) +
  geom_line(aes(window_ix, ploidy), color='gray70', size=3) +
  geom_line(aes(window_ix, norm_ploidy, color=dataset)) +
  facet_wrap(~ sample) +
  theme_bw() +
  scale_y_continuous(breaks=0:20)

ggplot(filter(data, in_viterbi)) +
  geom_line(aes(window_ix, ploidy), color='gray70', size=3) +
  geom_line(aes(window_ix, pred_ploidy, color=dataset), size=1.5) +
  geom_point(aes(window_ix, norm_ploidy, color=dataset)) +
  facet_wrap(~ region_group, scales='free') +
  scale_y_continuous('Normalized ploidy', breaks=0:20) +
  scale_x_continuous('Window') +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color='gray70'))
ggsave('~/Tmp/1.png', width=14, height=7)

ggplot(filter(data, in_viterbi)) +
  geom_tile(aes(region_group, ploidy), fill='gray60', height=0.2) +
  geom_boxplot(aes(region_group, norm_ploidy, fill=dataset)) +
  scale_y_continuous('Normalized ploidy', breaks=0:20) +
  scale_x_discrete('Region group') +
  scale_fill_discrete('Dataset') +
  theme_bw() +
  theme(panel.grid.major=element_line(color='gray75'))
ggsave('~/Tmp/2.png', width=14, height=7)

ggplot(filter(data, in_viterbi)) +
  geom_tile(aes(region_group, ploidy), fill='gray60', height=0.2) +
  geom_boxplot(aes(region_group, norm_ploidy_offset, fill=dataset)) +
  scale_y_continuous('Normalized ploidy (with offsets)', breaks=0:20) +
  scale_x_discrete('Region group') +
  scale_fill_discrete('Dataset') +
  theme_bw() +
  theme(panel.grid.major=element_line(color='gray75'))
ggsave('~/Tmp/3.png', width=14, height=7)



aggr_depth <- aggregate(norm_ploidy ~ sample + region_group + dataset,
                        data, FUN=function(x) c(mean = mean(x), sd = sd(x)))
aggr_depth$mean <- aggr_depth$norm_ploidy[,1]
aggr_depth$sd <- aggr_depth$norm_ploidy[,2]
aggr_depth$norm_ploidy <- NULL

ggplot(filter(aggr_depth, region_group == '04-06')) +
  geom_point(aes(mean, sd)) +
  facet_wrap(~ dataset) +
  scale_x_continuous('Mean normalized ploidy', breaks=0:20) +
  scale_y_continuous('Standard deviation') +
  theme_bw() + 
  theme(panel.grid.major.x = element_line(color = 'gray70'))

ggplot(aggr_depth, aes(dataset, mean)) +
  geom_boxplot(aes(fill=dataset)) +
  geom_jitter(alpha=.5) +
  facet_wrap(~ region_group, scales='free_y') +
  scale_y_continuous('Mean normalized ploidy', breaks=0:20) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = 'gray70'))
