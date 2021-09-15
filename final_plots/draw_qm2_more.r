library(cowplot)
library(grid)
library(gridExtra)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

data_dir <- '~/Data/hg38/jvc/comparisons/genes/'

abline_color <- 'gray80'
abline_size <- 1

plots <- list()
QUAL <- 20

# ====== SMN1 ======

gene <- 'SMN1'

df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- filter(df, region == 'SMN1_exon7')
df <- extend_paralog(df, 2)
df <- filter(df, agCN_filter == 'PASS' & agCN_qual >= QUAL &
               psCN_filter == 'PASS' & psCN_qual1 >= QUAL)
df <- add_noise(df, c('agCN', 'psCN1', 'psCN2'))

obs <- rbind(
  get_observations(df, 'SMN1 + SMN2', 'agCN_noise', 'b_copy_num',
                   region == 'SMN1_exon7'),
  get_observations(df, 'SMN1', 'psCN1_noise', 'b_paralog1',
                   region == 'SMN1_exon7'),
  get_observations(df, 'SMN2', 'psCN2_noise', 'b_paralog2',
                   region == 'SMN1_exon7')
)
obs$label <- factor(obs$label, levels=unique(obs$label))

(g1 <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid(. ~ label, scales='free', space='free') +
    theme_bw() +
    theme(axis.title = element_blank()))
(plots[[gene]] <- add_axis(g1,
         x_name = 'Parascopy CN estimate',
         y_name = 'QuicK-mer2 CN estimate',
         title = 'SMN1 (exons 7-8): 13 / 17 reliable PSVs'))
# cat em_f_values.csv | sed 's/chr5://' | awk '$2 >= 70948120 && $2 <= 70953015' | awk '$6 >= 0.95 && $7 >= 0.95' | wc -l

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

cn_match(df, region == 'SMN1_middle')
cn_match(df, region == 'SMN1_end')
par_cn_match(df, 2, region == 'SMN1_end')

# ====== ABCC6 ======

gene <- 'ABCC6'
df <- load(gene, 'qm2-matrix', keep_paralog=T, keep_b=T, keep_qual=T)
df <- load(gene, 'qm2-bed', keep_paralog=T, keep_b=T, keep_qual=T)
df <- filter(df, region == 'exon2')

df <- extend_paralog(df, 3)
df <- filter(df, agCN_qual >= QUAL & agCN_filter == 'PASS' &
      psCN_filter == 'PASS' & psCN_qual1 >= QUAL & psCN_qual2 >= QUAL & psCN_qual3 >= QUAL)
df <- add_noise(df, c('agCN', 'psCN1', 'psCN2', 'psCN3'))

obs <- rbind(
  get_observations(df, 'Aggregate', 'agCN_noise', 'b_copy_num', T),
  get_observations(df, 'ABCC6', 'psCN1_noise', 'b_paralog1', T),
  get_observations(df, 'ABCC6p1', 'psCN3_noise', 'b_paralog3', T),
  get_observations(df, 'ABCC6p2', 'psCN2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

# filter(df, region == 'exon2' & copy_num != round(b_copy_num)) %>% as.data.frame

(g1 <- ggplot(filter(obs, as.numeric(obs$label) == 1)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
  coord_fixed(ratio=.5) +
  facet_grid( ~ label) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(3, 3, 3, 3), "pt")))
(g2 <- ggplot(filter(obs, as.numeric(obs$label) > 1)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
  facet_grid( ~ label, scales='free', space='free') +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(3, 3, 3, 3), "pt")))

# cat em_f_values.csv | sed 's/chr.*://' | awk '$2 <= 16223494' | awk '$6 >= 0.95 && $7 >= 0.95 && (NF == 7 || $8 >= 0.95)' | wc -l
(plots[[gene]] <- add_axis(plot_grid(g1, g2, ncol=1),
                     x_name = 'Parascopy CN estimate',
                     y_name = 'QuicK-mer2 CN estimate',
                     title = 'ABCC6: 355 / 399 reliable PSVs'))

obs <- rbind(
  get_observations(df, 'ABCC6 + p1 + p2', 'agCN_noise', 'b_copy_num', T),
  get_observations(df, 'ABCC6', 'psCN1_noise', 'b_paralog1', T),
  get_observations(df, 'ABCC6p1', 'psCN3_noise', 'b_paralog3', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, space='free', scales='free') +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.margin = unit(c(3, 3, 3, 3), "pt")))
(plots[[gene]] <- add_axis(g1,
                           x_name = 'Parascopy CN estimate',
                           y_name = 'QuicK-mer2 CN estimate',
                           title = 'ABCC6: 355 / 399 reliable PSVs'))

# ====== GTF2I ======

gene <- 'GTF2I'

df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 3)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual &
               paralog_qual3 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2', 'paralog3'))

obs <- rbind(
  get_observations(df, 'Aggregate CN', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'GTF2I', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'GTF2Ip1', 'paralog3_noise', 'b_paralog3', T),
  get_observations(df, 'GTF2Ip4', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(filter(obs, as.numeric(obs$label) == 1)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    coord_fixed(ratio=.5) +
    facet_grid( ~ label) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.margin = unit(c(3, 3, 3, 3), "pt")))
(g2 <- ggplot(filter(obs, as.numeric(obs$label) > 1)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.margin = unit(c(3, 3, 3, 3), "pt")))
(plots[[gene]] <- add_axis(plot_grid(g1, g2, ncol=1),
     x_name = 'Parascopy CN estimate',
     y_name = 'QuicK-mer2 CN estimate',
     title = 'GTF2I: 25 / 116 reliable PSVs'))

# ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

# ====== CEL ======

gene <- 'CEL'

df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'CEL + CELP', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'CEL', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'CELP', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

g1 <- ggplot(filter(obs, as.numeric(obs$label) > 0)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      plot.margin = unit(c(3, 3, 3, 3), "pt"))
(plots[[gene]] <- add_axis(g1,
   x_name = 'Parascopy CN estimate',
   y_name = 'QuicK-mer2 CN estimate',
   title = 'CEL: 46 / 61 reliable PSVs'))

# ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

# ====== C4A ======

gene <- 'C4A'

df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'C4A + C4B', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'C4A', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'C4B', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

g1 <- ggplot(filter(obs, as.numeric(obs$label) > 0)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('', breaks=0:20) +
  facet_grid( ~ label, scales='free', space='free') +
  theme_bw() +
  theme(axis.title = element_blank(),
        plot.margin = unit(c(3, 3, 3, 3), "pt"))

(plots[[gene]] <- add_axis(g1,
     x_name = 'Parascopy CN estimate',
     y_name = 'QuicK-mer2 CN estimate',
     title = 'C4A: 0 / 17 reliable PSVs (7 / 50 in the vicinity)'))
# ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

# ====== NCF1 ======

gene <- 'NCF1'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 3)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual &
               paralog_qual3 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2', 'paralog3'))

obs <- rbind(
  get_observations(df, 'NCF1 + NCF1B + NCF1C', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'NCF1', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'NCF1B', 'paralog2_noise', 'b_paralog2', T),
  get_observations(df, 'NCF1C', 'paralog3_noise', 'b_paralog3', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(filter(obs, as.numeric(obs$label) == 1)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('', breaks=0:20) +
    coord_fixed(ratio=.5) +
    facet_grid( ~ label) +
    theme_bw() +
    theme(axis.title=element_blank(),
          plot.margin = unit(c(3, 3, 3, 3), "pt")))
(g2 <- ggplot(filter(obs, as.numeric(obs$label) > 1)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    theme_bw() +
    theme(axis.title = element_blank(),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))
(plots[[gene]] <- add_axis(plot_grid(g1, g2, ncol=1),
       x_name = 'Parascopy CN estimate',
       y_name = 'QuicK-mer2 CN estimate',
       title = 'NCF1: 43 / 92 reliable PSVs'))
# ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

# ====== SRGAP2 ======

gene <- 'SRGAP2'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 4)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual &
               paralog_qual3 >= qual & paralog_qual4 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2', 'paralog3', 'paralog4'))

obs <- rbind(
  get_observations(df, 'Aggregate', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'SRGAP2', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'SRGAP2B', 'paralog2_noise', 'b_paralog2', T),
  get_observations(df, 'SRGAP2C', 'paralog3_noise', 'b_paralog3', T),
  get_observations(df, 'SRGAP2D', 'paralog4_noise', 'b_paralog4', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(filter(obs, as.numeric(obs$label) == 1)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('', breaks=0:20) +
    coord_fixed(ratio=.5) +
    facet_grid( ~ label) +
    theme_bw() +
    theme(axis.title=element_blank(),
          plot.margin = unit(c(3, 3, 3, 3), "pt")))
(g2 <- ggplot(filter(obs, as.numeric(obs$label) > 1)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('', breaks=0:20) +
    facet_grid( ~ label) +
    theme_bw() +
    theme(axis.title = element_blank(),
          plot.margin = unit(c(3, 3, 3, 3), "pt")))
(plots[[gene]] <- add_axis(plot_grid(g1, g2, ncol=1),
                           x_name = 'Parascopy CN estimate',
                           y_name = 'QuicK-mer2 CN estimate',
                           title = 'SRGAP2: 1458 / 1916 reliable PSVs'))

cn_match(df)
par_cn_match(df, 4)
filter(df, copy_num != round(b_copy_num))$sample

# ====== RHCE ======

gene <- 'RHCE'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'RHCE + RHD', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'RHCE', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'RHD', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

g1 <- ggplot(filter(obs, as.numeric(obs$label) > 0)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('', breaks=0:20) +
  facet_grid( ~ label, scales='free', space='free') +
  theme_bw() +
  theme(axis.title = element_blank(),
    plot.margin = unit(c(3, 3, 3, 3), "pt"))

(plots[[gene]] <- add_axis(g1,
                           x_name = 'Parascopy CN estimate',
                           y_name = 'QuicK-mer2 CN estimate',
                           title = 'RHCE: 898 / 1097 reliable PSVs'))

cn_match(df)
par_cn_match(df, 2)

# ====== AMY1C ======

gene <- 'AMY1C'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- filter(df, copy_num_qual >= qual)

obs <- rbind(
  get_observations(df, 'AMY1C', 'copy_num', 'b_copy_num', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

g1 <- ggplot(filter(obs, as.numeric(obs$label) > 0)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=seq(0, 20, 2)) +
  scale_y_continuous('', breaks=seq(0, 20, 2)) +
  # facet_grid( ~ label, scales='free', space='free') +
  theme_bw() +
  theme(axis.title = element_blank(),
        plot.margin = unit(c(3, 3, 3, 3), "pt"))

(plots[[gene]] <- add_axis(g1,
                           x_name = 'Parascopy CN estimate',
                           y_name = 'QuicK-mer2 CN estimate',
      title = sprintf('AMY1C: correlation coef. = %.3f', cor(df$copy_num, df$b_copy_num))))

# ====== NPY4R ======

gene <- 'NPY4R'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- filter(df, copy_num_qual >= qual)
df <- add_noise(df, 'copy_num')

obs <- rbind(
  get_observations(df, 'NPY4R', 'copy_num_noise', 'b_copy_num', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(filter(obs, as.numeric(obs$label) > 0)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('', breaks=0:20) +
  # facet_grid( ~ label, scales='free', space='free') +
  theme_bw() +
  theme(axis.title = element_blank(),
        plot.margin = unit(c(3, 3, 3, 3), "pt")))

(plots[[gene]] <- add_axis(g1,
                           x_name = 'Parascopy CN estimate',
                           y_name = 'QuicK-mer2 CN estimate',
                           title = 'NPY4R: 0 / 11 reliable PSVs (3 / 34 in the vicinity)'))

# ====== FCGR3A ======

gene <- 'FCGR3A'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'FCGR3A + FCGR3B', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'FCGR3A', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'FCGR3B', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(filter(obs, as.numeric(obs$label) > 0)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
  facet_grid( ~ label, scales='free', space='free') +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    plot.margin = unit(c(3, 3, 3, 3), "pt")))
(plots[[gene]] <- add_axis(g1,
                           x_name = 'Parascopy CN estimate',
                           y_name = 'QuicK-mer2 CN estimate',
                           title = 'FCGR3A: 122 / 179 reliable PSVs'))

# ====== APOBEC3A ======

gene <- 'APOBEC3A'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'APOBEC3A + APOBEC3B', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'APOBEC3A', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'APOBEC3B', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(filter(obs, as.numeric(obs$label) > 0)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))
(plots[[gene]] <- add_axis(g1,
                           x_name = 'Parascopy CN estimate',
                           y_name = 'QuicK-mer2 CN estimate',
                           title = 'APOBEC3A: 44 / 51 reliable PSVs'))

# ====== HYDIN ======

gene <- 'HYDIN'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'HYDIN + HYDIN2', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'HYDIN', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'HYDIN2', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(g1 <- ggplot(filter(obs, as.numeric(obs$label) > 0)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))
(plots[[gene]] <- add_axis(g1,
                           x_name = 'Parascopy CN estimate',
                           y_name = 'QuicK-mer2 CN estimate',
                           title = 'HYDIN: 1678 / 2499 reliable PSVs'))

dir <- '~/Data/hg38/jvc/plots/qm2'
for (name in names(plots)) {
  cat(sprintf('%s\n', name))
  ggsave(file.path(dir, sprintf('%s.png', name)), plots[[name]],
         width=6, height=4)
}

# ====== ZP3 ======

gene <- 'ZP3'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'ZP3 + POMZP3', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'ZP3', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'POMZP3', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

g1 <- ggplot(filter(obs, as.numeric(obs$label) > 0)) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
  facet_grid( ~ label, scales='free', space='free') +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    plot.margin = unit(c(3, 3, 3, 3), "pt"))
(plots[[gene]] <- add_axis(g1,
                           x_name = 'Parascopy CN estimate',
                           y_name = 'QuicK-mer2 CN estimate',
                           title = 'ZP3: 18 / 46 reliable PSVs'))

# ====================
# Final
# ====================

point_color <- 'black'

# -------------------- ABCC6

gene <- 'ABCC6'
df <- load(gene, 'qm2-bed', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 3)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual &
               paralog_qual3 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2', 'paralog3'))

obs <- rbind(
  get_observations(df, 'ABCC6 + p1 + p2', 'copy_num_noise', 'b_copy_num', region == 'exon2'),
  get_observations(df, 'ABCC6', 'paralog1_noise', 'b_paralog1', region == 'exon2'),
  get_observations(df, 'ABCC6p1', 'paralog3_noise', 'b_paralog3', region == 'exon2'),
  get_observations(df, 'ABCC6p2', 'paralog2_noise', 'b_paralog2', region == 'exon2')
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(plots[[gene]] <- ggplot(filter(obs, as.numeric(obs$label) < 4)) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3, color=point_color) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20, expand=expansion(add=.3)) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    ggtitle('ABCC6: 355 / 399 reliable PSVs') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=.5, size=11, margin = margin(3, 0, 3, 0)),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))

# -------------------- FCGR3A

gene <- 'FCGR3A'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'FCGR3A + FCGR3B', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'FCGR3A', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'FCGR3B', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(plots[[gene]] <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3, color=point_color) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20, expand=expansion(add=.3)) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    ggtitle('FCGR3A: 122 / 179 reliable PSVs') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=.5, size=11, margin = margin(3, 0, 3, 0)),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))

# -------------------- SMN1

gene <- 'SMN1'
df <- load(gene, 'qm2-matrix', keep_paralog=T, keep_b=T, keep_qual=T)

keep_samples <- (df %>% select(sample, region, copy_num) %>%
   pivot_wider(id_cols='sample', names_from='region', values_from='copy_num') %>%
   filter(SMN1_exon2 == SMN1_exon7))$sample

df <- filter(df, sample %in% keep_samples)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'SMN1 + SMN2', 'copy_num_noise', 'b_copy_num',
                   region == 'SMN1_exon7'),
  get_observations(df, 'SMN1', 'paralog1_noise', 'b_paralog1',
                   region == 'SMN1_exon7'),
  get_observations(df, 'SMN2', 'paralog2_noise', 'b_paralog2',
                   region == 'SMN1_exon7')
)
obs$label <- factor(obs$label, levels=unique(obs$label))

(plots[[gene]] <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3, color=point_color) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20, expand=expansion(add=.3)) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    ggtitle('SMN1: 13 / 24 reliable PSVs') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=.5, size=11, margin = margin(3, 0, 3, 0)),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))

# -------------------- ZP3

gene <- 'ZP3'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, copy_num_qual >= qual & copy_num_filter == 'PASS' &
               paralog_filter == 'PASS' & paralog_qual1 >= qual & paralog_qual2 >= qual)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'ZP3 + POMZP3', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'ZP3', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'POMZP3', 'paralog2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(plots[[gene]] <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3, color=point_color) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20, expand=expansion(add=.3)) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    ggtitle('ZP3: 18 / 46 reliable PSVs') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=.5, size=11, margin = margin(3, 0, 3, 0)),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))

# --------- Merge

plot_grid(plots[['ABCC6']], plots[['FCGR3A']],
          plots[['SMN1']], plots[['ZP3']],
          nrow = 2, labels = LETTERS)
ggsave('~/Tmp/1.png', width=10, height=9, scale=.7, dpi=450)
