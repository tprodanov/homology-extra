library(ggplot2)
library(tidyverse)

RELIABLE_THRESHOLD <- 0.95

transform_f_values <- function(f_values, likelihoods, min_samples=30) {
  f_values <- filter(f_values, !is.na(copy1) & n_samples >= min_samples)
  f_values$mean_fval <- select(f_values, starts_with('copy')) %>%
    apply(1, mean, na.rm=T)
  f_values$min_fval <- select(f_values, starts_with('copy')) %>%
    apply(1, min, na.rm=T)
  f_values$reliable <- f_values$min_fval >= RELIABLE_THRESHOLD
  f_values
}

load_f_values <- function(upper_dir, region_name, dataset_names, dir_suffix='extra') {
  all_f_values <- data.frame()
  for (i in seq_along(dataset_names)) {
    dir <- sprintf('%s/%s/%s/%s', upper_dir, names(dataset_names)[i],
                   region_name, dir_suffix)
    if (!file.exists(sprintf('%s/em_likelihoods.csv', dir))) {
      cat(sprintf('WARN: Directory %s does not contain necessary files.\n', dir))
      next
    }
    f_values <- read.csv(sprintf('%s/em_f_values.csv', dir),
                         sep='\t', comment.char='#')
    f_values <- transform_f_values(f_values)
    if (nrow(f_values) == 0) {
      next
    }
    
    f_values$dataset <- dataset_names[[i]]
    all_f_values <- rbind(all_f_values, f_values)
  }
  all_f_values$dataset <- factor(all_f_values$dataset, levels=dataset_names)
  all_f_values
}

dataset_names <- c('244.AMR' = 'Admixed-American',
                   '243.AFR' = 'African',
                   '241.EAS' = 'East-Asian',
                   '242.EUR' = 'European',
                   '245.SAS' = 'South-Asian')

upper_dir <- '~/Data/hg38/jvc/runs'
plots_dir <- '~/Data/hg38/jvc/plots/population_comparison/r009/f_val'
region_name <- 'SMN1'

all_f_values <- load_f_values(upper_dir, region_name, dataset_names, 'r009/extra')
f_values <- filter(all_f_values, region_group == '02-01') %>%
  select(!c('copy3', 'copy4'))
f_values$pos <- as.numeric(sub('chr5:', '', f_values$psv))

f_values2 <- select(f_values, !'reliable') %>%
  pivot_longer('info_content' | starts_with('copy') | ends_with('fval'),
               names_to='stat', values_to='value')

f_values2$stat <- as.character(f_values2$stat)
f_values2 <- mutate(f_values2,
                    stat = case_when(stat == 'info_content' ~ 'Information content',
                                     stat == 'mean_fval' ~ 'Mean f-value',
                                     stat == 'min_fval' ~ 'Minimal f-value',
                                     stat == 'copy1' ~ 'Copy 1',
                                     stat == 'copy2' ~ 'Copy 2',
                                     stat == 'copy3' ~ 'Copy 3',
                                     stat == 'copy4' ~ 'Copy 4',
                                     stat == 'copy5' ~ 'Copy 5',
                                     T ~ stat))

smn_pos <- c(70925030, 70952347)
f_values2$in_smn <- with(f_values2, smn_pos[1] <= pos & pos <= smn_pos[2])
f_values2$pos_f <- factor(f_values2$pos)
uniq_pos <- levels(f_values2$pos_f)
f_values2$ix <- as.numeric(f_values2$pos_f)

smn_caller_psvs <- c(70950493, 70950966, 70951392, 70951463,
                     70951897, 70951946, 70952094, 70952209)
x_colors <- ifelse(uniq_pos %in% as.character(smn_caller_psvs),
                   'red', 'black')

colors <- ggthemes::tableau_color_pal()(10)[c(1, 2, 3, 5, 7)]
smn_start_ix <- min(filter(f_values2, pos >= smn_pos[1])$ix)
start_ix <- max(0, smn_start_ix - 10)
smn_end_ix <- max(filter(f_values2, pos <= smn_pos[2])$ix)
end_ix <- min(length(uniq_pos), smn_end_ix + 10)

del_pos <- 70948120
del_ix <- min(filter(f_values2, pos >= del_pos)$ix)

f_values_filt <- filter(f_values2, startsWith(stat, 'Copy') & ix >= start_ix & ix <= end_ix)
rel_count <- filter(f_values, reliable) %>% count(pos)
blue_vlines <- filter(rel_count, n >= length(dataset_names) - 1)$pos
blue_vlines <- match(as.character(blue_vlines), uniq_pos)

f_values_filt$label <- ifelse(f_values_filt$stat == 'Copy 1', 'SMN1', 'SMN2')

uniq_pos2 <- sprintf('%s%s', ifelse(uniq_pos %in% smn_caller_psvs, '🞷 ', ''),
                     as.numeric(uniq_pos) - 70e6)

boundaries <- data.frame(pos=c(smn_start_ix - 0.5, smn_end_ix + 0.5, del_ix - 0.5),
      linetype=c(rep('Gene', 2), 'Deletion'))
boundaries$linetype <- factor(boundaries$linetype,
                              levels=unique(boundaries$linetype))

(g_smn_b <- ggplot(f_values_filt) +
  #geom_vline(xintercept=c(smn_start_ix - 0.5, smn_end_ix + 0.5),
             #linetype='dashed', color='gray50') +
  #geom_vline(xintercept=del_ix,
             #linetype='23', color='gray50') +
  geom_segment(aes(x=pos, xend=pos, y=-Inf, yend=Inf, linetype=linetype),
               data=boundaries, color='gray30') +
  geom_vline(xintercept = blue_vlines, size=3, color='blue', alpha=.1) +
  geom_point(aes(ix, value, color=dataset),
             position=position_dodge(width=0.5)) +
  geom_hline(yintercept = RELIABLE_THRESHOLD) +
  scale_x_continuous('PSV position (-70 Mb)',
                     breaks=start_ix:end_ix, minor_breaks=NULL,
                     labels=function(i) uniq_pos2[i],
                     expand = c(0.01, 0.01)) +
  scale_y_continuous('Frequency of the reference allele') +
  scale_color_manual('', values=colors) +
  scale_linetype_manual('', values=c('dashed', '13')) +
  facet_wrap(~ label, ncol=1) +
  guides(color = guide_legend(override.aes = list(alpha=1, size=2.5))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=7,
                                   color=x_colors[start_ix:end_ix]),
        legend.position='top',
        legend.margin=margin(-5, 0, -10, 0),
        legend.key.size = unit(15, 'pt')))
ggsave('~/Tmp/1.png', width=10, height=6)

f_values$ix <- match(as.character(f_values$pos), uniq_pos)
ggplot(filter(f_values, ix >= start_ix & ix <= end_ix)) +
  #annotate(geom='segment', x=0.95, xend=1, y=0.95, yend=0.95) +
  #annotate(geom='segment', x=0.95, xend=0.95, y=0.95, yend=1) +
  #annotate(geom='segment', x=0.8, xend=1, y=0.8, yend=0.8) +
  #annotate(geom='segment', x=0.8, xend=0.8, y=0.8, yend=1) +
  annotate(geom='rect', xmin=0.8, xmax=1, ymin=0.8, ymax=1, fill='gold', alpha=.3) +
  annotate(geom='rect', xmin=0.95, xmax=1, ymin=0.95, ymax=1, fill='lawngreen', alpha=1) +
  geom_point(aes(copy1, copy2, color=!(pos %in% smn_caller_psvs)), alpha=.5) +
  facet_wrap(~ dataset) +
  scale_color_manual('PSV is used in\nSMNCopyNumberCaller', values=c('red', 'black'),
                     labels=c('Yes', 'No')) +
  theme_bw() +
  theme(legend.position=c(0.92, 0.11), legend.justification=c('right', 'bottom'))
ggsave('~/Tmp/2.png', width=10, height=6)

length(uniq_pos)
smn_end_ix - smn_start_ix + 1

serf1a_pos <- c(70900687, 70918530)
length(unique(filter(f_values, pos >= serf1a_pos[1] & pos <= serf1a_pos[2])$pos))
end_ix - start_ix + 1

filter(f_values, ix >= start_ix & ix <= end_ix & reliable) %>% count(dataset)
filter(f_values, ix >= start_ix & ix <= end_ix & copy1 >= 0.8 & copy2 >= 0.8) %>%
  count(dataset)

filter(f_values, ix >= start_ix & ix <= end_ix & reliable) %>%
  count(psv) %>% filter(n >= 4) %>% nrow
filter(f_values, ix >= start_ix & ix <= end_ix & copy1 >= 0.8 & copy2 >= 0.8) %>%
  count(psv) %>% filter(n >= 4) %>% nrow
