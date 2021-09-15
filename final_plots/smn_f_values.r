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

dataset_names <- c('203.han_g1k' = '1000 genomes',
                   '204.han_bgi' = 'High coverage BGI')


upper_dir <- '~/Data/hg38/jvc/runs'
plots_dir <- '~/Data/hg38/jvc/plots/population_comparison/v100/f_val'
region_name <- 'SMN1'

all_f_values <- load_f_values(upper_dir, region_name, dataset_names, 'v100/extra')
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

# =============
# uniq_pos <- setdiff(uniq_pos, c('70952674', '70954868'))
# =============

f_values2$ix <- match(f_values2$pos_f, uniq_pos)
smn_caller_psvs <- c(70950493, 70950966, 70951392, 70951463,
                     70951897, 70951946, 70952094, 70952209)
x_colors <- ifelse(uniq_pos %in% as.character(smn_caller_psvs),
                   'red', 'black')

# colors <- ggthemes::tableau_color_pal()(10)[c(1, 2, 3, 5, 7)]
colors <- ggthemes::tableau_color_pal()(10)[c(1, 2, 3, 5)]
smn_start_ix <- min(which(as.numeric(uniq_pos) >= smn_pos[1]))
start_ix <- max(0, smn_start_ix - 10)
smn_end_ix <- max(which(as.numeric(uniq_pos) <= smn_pos[2]))
end_ix <- min(length(uniq_pos), smn_end_ix + 10)

del_pos <- 70948287
del_ix <- min(which(as.numeric(uniq_pos) >= del_pos))

f_values_filt <- filter(f_values2, startsWith(stat, 'Copy') & ix >= start_ix & ix <= end_ix)
rel_count <- filter(f_values, reliable) %>% count(pos)
blue_vlines <- filter(rel_count, n >= length(dataset_names) - 1)$pos
blue_vlines <- match(as.character(blue_vlines), uniq_pos)
# present_in_all <- (count(f_values, pos) %>% filter(n == length(dataset_names)))$pos

f_values_filt$label <- ifelse(f_values_filt$stat == 'Copy 1', 'SMN1', 'SMN2')

# '＊'
# '✳'
# '🞷'
uniq_pos2 <- sprintf('%s%s', ifelse(uniq_pos %in% smn_caller_psvs, '✳', ''),
                     format(as.numeric(uniq_pos) - 70.9e6))

boundaries <- data.frame(pos=c(smn_start_ix - 0.5, smn_end_ix + 0.5, del_ix - 0.5),
      linetype=c(rep('Gene', 2), 'Deletion'))
boundaries$linetype <- factor(boundaries$linetype,
                              levels=unique(boundaries$linetype))

(g_smn_b <- ggplot(filter(f_values_filt, dataset != 'Admixed-American')) +
  geom_segment(aes(x=pos, xend=pos, y=-Inf, yend=Inf, linetype=linetype),
               data=boundaries, color='gray30') +
  #geom_vline(xintercept = blue_vlines, size=3, color='blue', alpha=.1) +
  geom_hline(yintercept = RELIABLE_THRESHOLD) +
  geom_point(aes(ix, value, color=dataset), size=1.1,
             position=position_dodge(width=0.5)) +
  scale_x_continuous('PSV position (starting at 70.9 Mb)',
                     breaks=start_ix:end_ix, minor_breaks=NULL,
                     labels=function(i) uniq_pos2[i],
                     expand = c(0.01, 0.01)) +
  # scale_y_continuous('Frequency of the ref. allele') + # TODO: RETURN
  scale_y_continuous('Frequency of the ref. allele      ') +
  scale_color_manual('', values=colors) +
  # scale_linetype_manual('Boundary', values=c('dashed', '13')) + # TODO: RETURN
  scale_linetype_manual('', values=c('dashed', '13')) +
  facet_wrap(~ label, ncol=1, strip.position = 'right') +
  guides(color = guide_legend(override.aes = list(alpha=1, size=2.5))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8,
                                   color=x_colors[start_ix:end_ix]),
        #axis.title.y = element_text(hjust=1),
        legend.position='top',
        legend.margin=margin(-5, 0, -10, 5),
        legend.title = element_text(size=9),
        legend.key.height = unit(10, 'pt'),
        # legend.key.width = unit(28, 'pt'),
        legend.spacing.x = unit(1, 'pt'), # TODO: REMOVE
        plot.margin = margin(10, 5, 5, 5),
        strip.text = element_text(margin=margin(0, 2, 0, 2))))
ggsave('~/Tmp/1.png', width=10, height=6, scale=.6, dpi=450)
ggsave('~/Tmp/1.png', width=10, height=4, scale=.55, dpi=600)

length(unique(filter(f_values_filt, dataset != 'Admixed-American')$pos))

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

f_values_wide <- filter(f_values2, T) %>% select(dataset, pos, stat, value) %>%
  pivot_wider(names_from='stat', values_from='value')
f_values_wide$min_fval <- pmin(f_values_wide[['Copy 1']], f_values_wide[['Copy 2']])
f_values_wide$reliable <- f_values_wide$min_fval >= 0.95
filter(f_values_wide, reliable) %>% count(dataset)

f_values_wide2 <- select(f_values_wide, pos, dataset, min_fval) %>%
  pivot_wider(names_from='dataset', values_from='min_fval')
names(f_values_wide2) <- c('pos', 'g1k', 'bgi')

f_values_wide3 <- filter(f_values_wide2, !is.na(g1k) & !is.na(bgi))

ggplot(f_values_wide2) +
  geom_point(aes(g1k, bgi))

nrow(f_values_wide3)
sum(with(f_values_wide3, pos >= smn_pos[1] & pos <= smn_pos[2]))

sum(with(f_values_wide3, g1k >= 0.95))
sum(with(f_values_wide3, bgi >= 0.95))
sum(with(f_values_wide3, (g1k >= 0.95) != (bgi >= 0.95)))

sum(with(f_values_wide3, g1k >= 0.8))
sum(with(f_values_wide3, bgi >= 0.8))
filter(f_values_wide3, (g1k >= 0.8) != (bgi >= 0.8))

with(f_values_wide3, cor(g1k, bgi))
with(filter(f_values_wide3, g1k >= 0.6 | bgi >= 0.6), cor(g1k, bgi))
