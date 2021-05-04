library(cowplot)

samples <- c('HG00233', 'HG00103', 'HG00107', 'HG00137', 'HG00130', 'NA12874')
samples <- c('HG00233', 'HG00103', 'HG00107', 'HG00137', 'HG00130', 'NA12006')
sample_labels <- LETTERS[1:length(samples)]
names(sample_labels) <- samples
min_ix <- 292
max_ix <- 561

fdepth_v <- filter(depth_v, sample %in% samples &
                     window_ix >= min_ix & window_ix <= max_ix) %>%
  mutate(sample = factor(sample, levels=samples))
fdepth_v$middle <- windows_v[match(fdepth_v$window_ix, windows_v$window_ix),]$middle
fwindows_v <- filter(windows_v, window_ix >= min_ix & window_ix <= max_ix)
fdepth_v.av <- aggregate.average(fdepth_v, fwindows_v, 'corr_cn1', n=5)

f_states <- filter(hmm_states, sample %in% samples &
                     window_ix >= min_ix & window_ix <= max_ix &
                     (iteration == '1' | iteration == 'v')) %>%
  mutate(sample = factor(sample, levels=samples))

hmm_wide <- filter(hmm_states, iteration == 'v' &
                     (window_ix == min_ix | window_ix == max_ix)) %>%
  select(c('sample', 'window_ix', 'pred_cn')) %>%
  pivot_wider(names_from='window_ix', values_from='pred_cn')
colnames(hmm_wide) <- c('sample', 'start', 'end')
filter(hmm_wide, start != end)

hmm_wide <- filter(hmm_states, window_ix == max_ix &
                     (iteration == '1' | iteration == 'v')) %>%
  select(c('sample', 'iteration', 'pred_cn')) %>%
  pivot_wider(names_from='iteration', values_from='pred_cn')
colnames(hmm_wide) <- c('sample', 'start', 'end')
filter(hmm_wide, round(start) != end)

fdepth_v.av$label <- sample_labels[fdepth_v.av$sample]
colors <- palette.colors(n=6, 'Dark2')
del_pos <- 70948120
(g1 <- ggplot(fdepth_v.av) +
  geom_vline(xintercept=del_pos, linetype='dashed', color='gray60') +
  geom_hline(yintercept=c(3, 4, 5), color='gray60') +
  #geom_line(aes(middle, pred_cn, linetype='First'),
            #data=filter(f_states, iteration == '1')) +
  #geom_line(aes(middle, pred_cn, linetype='Last'),
            #data=filter(f_states, iteration == 'v')) +
  geom_point(aes(middle, corr_cn1, color=label), alpha=1, size=1.3) +
  facet_wrap(~ label, strip.position='right', ncol=1) +
  scale_x_continuous('Position (chr5, kb)', minor_breaks=NULL,
                     breaks=70e6 + c(930, 940, 950) * 1e3,
      labels=function(x) format(x / 1e3, big.mark=',', digits=0, scientific=F)) +
  scale_y_continuous('Normalized pooled read depth', breaks=1:10, minor_breaks=NULL) +
  scale_color_manual(values=colors) +
  scale_linetype_manual('HMM iteration', values=c('twodash', 'solid')) +
  guides(color=F) +
  theme_bw() +
  theme(legend.position='top', legend.margin=margin(-2, 0, -5, 0),
        legend.key.width = unit(20, 'pt'),
        strip.text.y = element_text(angle=0, margin = margin(0, 2, 0, 2, 'pt')),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)))
ggsave('~/Tmp/1.png', g1, width=4, height=4)

f_states$label <- sample_labels[f_states$sample]
f_states$iter_label <- ifelse(f_states$iteration == '1',
                              'First HMM iteration', 'Last HMM iteration')

annotation <- rbind(
  data.frame(middle=70940000,
             pred_cn=rep(c(5, 4, 3), 2),
             label=rep(c('A', 'B D E F', 'C'), 2),
             iter_label = rep(c('First HMM iteration', 'Last HMM iteration'), each=3)),
  data.frame(middle=70951370,
             pred_cn=c(5, 4, 3),
             label=c('A', 'B', 'C D E F'),
             iter_label = 'Last HMM iteration'),
  filter(f_states,
         (window_ix == 550 & label %in% c('A', 'B', 'C')) |
         (window_ix == 522 & label == 'D') |
         (window_ix == 528 & label == 'E') |
         (window_ix == 534 & label == 'F')) %>%
    filter(iteration == '1') %>%
    select(c('middle', 'pred_cn', 'label', 'iter_label'))
)
annotation$nudge_x <- with(annotation,
    case_when(label == 'D' | label == 'E' ~ 200,
              label == 'F' ~ 800,
              T ~ 0))
annotation$nudge_y <- with(annotation, 
    case_when(label == 'D' | label == 'E' ~ 0.2,
              label == 'F' ~ -0.2,
              T ~ 0.2))

(g2 <- ggplot(f_states) +
  geom_vline(xintercept=del_pos, linetype='dashed', color='gray60') +
  geom_line(aes(middle, pred_cn, color=label), size=1, alpha=.8) +
  geom_text(aes(middle, pred_cn, label=label),
                  data=annotation, size=3,
                  nudge_x=annotation$nudge_x,
                  nudge_y=annotation$nudge_y) +
  facet_wrap(~ iter_label, ncol=1) +
  coord_cartesian(xlim=c(70935000, NA)) +
  scale_x_continuous('Position (chr5, kb)', minor_breaks=NULL,
    labels=function(x) format(x / 1e3, big.mark=',', digits=0, scientific=F)) +
  scale_y_continuous('Copy number estimate', breaks=1:10, minor_breaks=NULL,
                     limits=c(2.8, 5.3)) +
  scale_color_manual(values=colors) +
  guides(color=F) +
  theme_bw() +
  theme(strip.text.x = element_text(margin = margin(1, 0, 1, 0, 'pt'))))
ggsave('~/Tmp/2.png', g2, width=5, height=3, scale=0.85)

plot_grid(g1, g2, ncol=1, rel_heights=c(1, 0.65))# , labels=c('C', 'D'))
ggsave('~/Tmp/3.png', width=6, height=9, scale=.8)





data_dir <- '~/Data/hg38/jvc/runs/245.SAS/SMN1/r009/extra/'
psv_obs <- read.csv(sprintf('%s/psv_observations.csv', data_dir),
                    sep='\t', comment.char='#')
f_values <- read.csv(sprintf('%s/em_f_values.csv', data_dir),
                     sep='\t', comment.char='#')
f_values <- filter(f_values, region_group == '02-01')

em_lik <- read.csv(sprintf('%s/em_likelihoods.csv', data_dir),
                   sep='\t', comment.char='#')
em_lik <- em_lik %>% filter(region_group == '02-01') %>%
  group_by(cluster) %>% slice_tail(n=1) %>% ungroup()
ix <- which.max(em_lik$likelihood)
best_cluster <- em_lik$cluster[ix]
best_iter <- em_lik$iteration[ix]

sample_gts <- read.csv(sprintf('%s/em_sample_gts.csv', data_dir),
                       sep='\t', comment.char='#')
sample_gts <- filter(sample_gts, region_group == '02-01' &
                       cluster == best_cluster & iteration == best_iter)
# Remove samples with aggr copy number != 4.
sample_gts <- sample_gts[, !is.na(sample_gts[1,])]
best_gts <- sample_gts[, 6:ncol(sample_gts)] %>% apply(2, which.max)
best_gts <- sample_gts$genotype[best_gts]
names(best_gts) <- colnames(sample_gts)[6:ncol(sample_gts)]
psv_obs$sample_gt <- best_gts[psv_obs$sample]

# psvs <- c(70924721, 70927608, 70950340, 70952209)
# psvs <- c(70909380, 70921845, 70928300, 70943791, 70950966, 70952209)
psvs <- c(70909380, 70928300, 70952209)
fpsvs_obs <- filter(psv_obs, pos %in% psvs & !is.na(sample_gt)) %>%
  mutate(alt_cov = as.numeric(alt_cov),
         ref_fraq = ref_cov / (ref_cov + alt_cov))

fpsvs_obs$psv <- sprintf('chr5:%d', fpsvs_obs$pos)
ixs <- match(fpsvs_obs$psv, f_values$psv)
fpsvs_obs$label <- with(f_values[ixs,],
    sprintf('%s\nf1 = %.3f,  f2 = %.3f\ninformation content = %.3f',
            psv, copy1, copy2, info_content))

binwidth <- 0.05
ylim <- 125
col_height <- sum(with(fpsvs_obs, pos == psvs[1] & ref_fraq >= 1 - binwidth))
annotate_df <- data.frame(label = fpsvs_obs$label[1],
                          x=1 - 2 * binwidth, y=ylim - 10,
                          text=sprintf('ᐱ\n%d', col_height))

fill_colors <- ggthemes::tableau_color_pal()(10)
ggplot(fpsvs_obs) +
  #geom_vline(xintercept=c(0, 1), color='gray80') +
  geom_histogram(aes(ref_fraq, fill=sample_gt), center=0.5,
                 color='black', size=.1, binwidth=binwidth, position="stack") +
  geom_text(aes(x, y, label=text), data=annotate_df, size=3) +
  facet_wrap(~ label) +
  scale_y_continuous('Number of samples') +
  scale_x_continuous('Fraction of the allele corresponding to the copy 1') +
  coord_cartesian(ylim=c(0, 125)) +
  scale_fill_manual('Sample\ngenotype', values=fill_colors) +
  theme_bw()
ggsave('~/Tmp/3.png', width=10, height=3)

fpsvs_obs2 <- filter(fpsvs_obs, pos == psvs[3])
fpsvs_obs3 <- rbind(
  mutate(fpsvs_obs2, flag = 'B',
         grid = 'All samples'),
  mutate(fpsvs_obs2, flag = ifelse(sample_gt == '2,2', 'B', 'A'),
         grid = 'Samples with genotype 2,2')
)

ggplot(fpsvs_obs3) +
  geom_histogram(aes(ref_fraq, fill=flag), center=0.5,
                 color='black', size=.1, binwidth=binwidth, position="stack") +
  facet_wrap(~ grid) +
  scale_y_continuous('Number of samples') +
  scale_x_continuous('Fraction of the allele corresponding to the copy 1',
                     breaks=0:4 / 4) +
  scale_fill_manual(values=c('#d6e0eb', '#4E79A7'), guide=F) +
  theme_bw() +
  theme(strip.text.x = element_text(margin = margin(1, 0, 1, 0, 'pt')),
        axis.title.x = element_text(size=9, margin=margin(0, 0, -5, 0)),
        axis.title.y = element_text(size=9, margin=margin(0, -2, 0, -5)))
ggsave('~/Tmp/3.png', width=6, height=2, scale=.9)
#bacbde
###

samples <- c('HG03963', 'HG03895', 'HG03611', 'NA21129')
samples <- 'HG03895'
psv_obs$psv <- sprintf('chr5:%d', psv_obs$pos)
sample_obs <- filter(psv_obs, sample %in% samples & use & psv %in% f_values$psv) %>%
  mutate(alt_cov = as.numeric(alt_cov),
         ref_fraq = ref_cov / (ref_cov + alt_cov))
sample_obs$sample <- factor(sample_obs$sample, levels=samples)
sample_obs$f1 <- f_values[match(sample_obs$psv, f_values$psv),]$copy1
sample_obs$f2 <- f_values[match(sample_obs$psv, f_values$psv),]$copy2
sample_obs$f_case <- with(sample_obs, case_when(
  f1 >= 0.95 & f2 >= 0.95 ~ 'Reliable',
  f1 >= 0.8 & f2 >= 0.8 ~ 'Semi-reliable',
  f1 < 0.8 & f2 >= 0.8 ~ 'Low f1',
  f2 < 0.8 & f1 >= 0.8 ~ 'Low f2',
  T ~ 'Low f1 and f2'
))
sample_obs$f_case <- factor(sample_obs$f_case,
    levels=rev(c('Reliable', 'Semi-reliable', 'Low f1', 'Low f2', 'Low f1 and f2')))

labels <- c(
  #'Sample genotype = 2,2',
  'Sample genotype = 2,2'
  #'Sample genotype = 2,2.    Has gene conversion',
  #'Sample genotype = 3,1'
)
names(labels) <- samples

scales::show_col(fill_colors)
scales::show_col(c('red2', 'gold2', 'darkorange','deepskyblue', 'blue3'))
f_val_colors <- fill_colors[c(3, 6, 2, 4, 1)]
ggplot(sample_obs) +
  geom_histogram(aes(ref_fraq, fill=f_case), binwidth=binwidth,
                 color='black', size=.1) +
  facet_wrap(~ sample, labeller = as_labeller(labels)) +
  scale_x_continuous('Fraction of the allele corresponding to the copy 1') +
  scale_y_continuous('Number of PSVs') +
  scale_fill_manual('PSV type',
    values=f_val_colors,
    guide = guide_legend(reverse = TRUE)) +
  theme_bw() 
ggsave('~/Tmp/4.png', width=10, height=6)

filter(sample_obs, sample == samples[3] & f_case == 'Semi-reliable' & ref_fraq < 0.35)

sample_obs2 <- rbind(
  mutate(sample_obs, flag='B', grid='All PSVs'),
  mutate(sample_obs, flag=ifelse(f_case == 'Reliable', 'B', 'A'), grid='Reliable PSVs'))
         
ggplot(sample_obs2) +
  geom_histogram(aes(ref_fraq, fill=flag), binwidth=binwidth,
                 size=.1, color='black') +
  facet_wrap(~ grid) +
  scale_x_continuous('Fraction of the allele corresponding to the copy 1') +
  scale_y_continuous('Number of PSVs') +
  scale_fill_manual(values=c('#d6e0eb', '#4E79A7'), guide=F) +
  #scale_fill_manual(values=c('#e0ccdb', '#8f567f'), guide=F) +
  theme_bw() +
  theme(strip.text.x = element_text(margin = margin(1, 0, 1, 0, 'pt')),
        axis.title.x = element_text(size=9, margin=margin(0, 0, -5, 0)),
        axis.title.y = element_text(size=9, margin=margin(0, 0, 0, -5)))
ggsave('~/Tmp/4.png', width=6, height=2, scale=.9)
