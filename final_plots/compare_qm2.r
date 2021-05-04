library(cowplot)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

data_dir <- '~/Data/hg38/jvc/comparisons/populations/r009'
plot_dir <- '~/Data/hg38/jvc/plots/population_comparison/r009/experimental'

abline_color <- 'gray80'
abline_size <- 1

# ====== ABCC6 ======

gene <- 'ABCC6'
method <- 'QuicK-mer2'; method_short <- 'qm2'
y_label <- sprintf('%s CN estimate', method)

df <- load(gene, method_short, keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 3)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2', 'paralog3'))

obs <- rbind(
  get_observations(df, 'Aggregate', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'ABCC6', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'ABCC6p1', 'paralog3_noise', 'b_paralog3', T),
  get_observations(df, 'ABCC6p2', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(filter(obs, as.numeric(obs$label) == 1)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  coord_fixed(ratio=.5) +
  facet_grid( ~ label) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        plot.margin = unit(c(3, 3, 3, 3), "pt")))
(g2 <- ggplot(filter(obs, as.numeric(obs$label) > 1)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  facet_grid( ~ label, scales='free', space='free') +
  theme_bw() +
  theme(plot.margin = unit(c(3, 3, 3, 3), "pt")))
plot_grid(g1, g2, ncol=1)
ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table


# ====== ABCC6 ======

gene <- 'GTF2I'
method <- 'QuicK-mer2'; method_short <- 'qm2'
y_label <- sprintf('%s CN estimate', method)

df <- load(gene, method_short, keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 3)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2', 'paralog3'))

obs <- rbind(
  get_observations(df, 'Aggregate', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'GTF2I', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'GTF2Ip1', 'paralog3_noise', 'b_paralog3', T),
  get_observations(df, 'GTF2Ip4', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(filter(obs, as.numeric(obs$label) == 1)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous(y_label, breaks=0:20) +
    coord_fixed(ratio=.5) +
    facet_grid( ~ label) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          plot.margin = unit(c(3, 3, 3, 3), "pt")))
(g2 <- ggplot(filter(obs, as.numeric(obs$label) > 1)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous(y_label, breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    theme_bw() +
    theme(plot.margin = unit(c(3, 3, 3, 3), "pt")))
plot_grid(g1, g2, ncol=1)
ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

# ====== CEL ======

gene <- 'CEL'
method <- 'QuicK-mer2'; method_short <- 'qm2'
y_label <- sprintf('%s CN estimate', method)

df <- load(gene, method_short, keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, paralog_filter == 'PASS' & paralog_qual1 >= 30)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'CEL + CELP', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'CEL', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'CELP', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

ggplot(filter(obs, as.numeric(obs$label) > 0)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous(y_label, breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    theme_bw() +
    theme(plot.margin = unit(c(3, 3, 3, 3), "pt"))
ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

# ====== C4A ======

gene <- 'C4A'
method <- 'QuicK-mer2'; method_short <- 'qm2'
y_label <- sprintf('%s CN estimate', method)

df <- load(gene, method_short, keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, paralog_filter == 'PASS' & paralog_qual1 >= 30)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'C4A + C4B', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'C4A', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'C4B', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

ggplot(filter(obs, as.numeric(obs$label) > 0)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  facet_grid( ~ label, scales='free', space='free') +
  theme_bw() +
  theme(plot.margin = unit(c(3, 3, 3, 3), "pt"))
ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

