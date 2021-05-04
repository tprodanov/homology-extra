library(cowplot)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

data_dir <- '~/Data/hg38/jvc/comparisons/populations/r009'
plot_dir <- '~/Data/hg38/jvc/plots/population_comparison/r009/experimental'

abline_color <- 'gray80'
abline_size <- 1

# ====== SMN1 ======

# ------ SMNCopyNumberCaller ------

gene <- 'SMN1'
method <- 'SMNCopyNumberCaller'; method_short <- 'caller'
y_label <- sprintf('%s CN estimate', method)
df <- load(gene, method_short, keep_qual=T, keep_b=T)
df <- extend_paralog(df, 2)
df$method <- sub('caller_', '', df$method)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2', 'b_paralog1', 'b_paralog2'))

obs <- rbind(
  get_observations(df, 'SMN1 + SMN2 exons 1-6', 'copy_num_noise', 'b_copy_num',
                   method == '16'),
  get_observations(df, 'SMN1 exons 7-8', 'paralog1_noise', 'b_paralog1_noise',
                   method == '78' & paralog_qual1 >= 30),
  get_observations(df, 'SMN2 exons 7-8', 'paralog2_noise', 'b_paralog2_noise',
                   method == '78' & paralog_qual1 >= 30)
)

(g_smn_c <- ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=0.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(bquote(bold(.(method)) ~ 'CN estimate'),
                     breaks=0:20) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw())
ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=12, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

filter(df, method == '78' & copy_num == round(b_copy_num)) %>%
  mutate(eq = paralog1 == b_paralog1) %>% count(eq)

filter(df, method == '78' & !is.na(paralog1) & is.na(b_paralog1)) %>% nrow
filter(df, method == '78' & !is.na(paralog1) & is.na(b_paralog1) &
       paralog_qual1 >= 30 & paralog_filter == 'PASS') %>% nrow
filter(df, method == '78' & is.na(paralog1) & !is.na(b_paralog1)) %>% nrow

# ------ MLPA ------

gene <- 'SMN1'
method <- 'MLPA'; method_short <- 'mlpa'
y_label <- sprintf('%s CN estimate', method)

df <- load(gene, method_short, keep_paralog=F, keep_b=T, keep_qual=T)
# For some reasons, some samples have two values in MLPA.
df <- group_by(df, sample, method) %>% slice_head(n=1) %>% ungroup
df$method <- sub('mlpa_', '', df$method)
df <- add_noise(df, c('copy_num'))

obs <- rbind(
  get_observations(df, 'SMN1 + SMN2 exons 1-6', 'copy_num_noise', 'b_copy_num',
                   method == '16'),
  get_observations(df, 'SMN1 + SMN2 exons 7-8', 'copy_num_noise', 'b_copy_num',
                   method == '78'))

(g_smn_d <- ggplot(filter(obs, population != 'BEB' & population != 'TSI')) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(bquote(bold(.(method)) ~ 'CN estimate'),
                     breaks=0:20) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw())
ggsave(sprintf('%s/%s.%s.v1.png', plot_dir, gene, method_short), width=10, height=6)

obs$pop2 <- factor(obs$population, levels=c('BEB', 'TSI', 'Other')) %>%
  replace_na('Other')
ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs, color=pop2), alpha=.5) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  scale_color_brewer(palette='Set1', 'Population') +
  facet_grid(. ~ label, scales='free', space='free') +
  guides(color = guide_legend(override.aes = list(alpha=1, size=2.5))) +
  theme_bw() +
  theme(legend.position=c(0.995, 0.005),
        legend.justification = c('right', 'bottom'))
ggsave(sprintf('%s/%s.%s.v2.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  filter(population != 'BEB' & population != 'TSI') %>%
  select(label, eq) %>% table

plot_grid(g_smn_a, g_smn_b, g_smn_c, g_smn_d, labels=LETTERS)
ggsave('~/Tmp/1.png', width=15, height=8)

# ------ QuicK-mer2 ------

gene <- 'SMN1'
method <- 'QuicK-mer2'; method_short <- 'qm2'
y_label <- sprintf('%s CN estimate', method)

keep_samples <- local({
  caller_df <- load(gene, 'caller', keep_paralog=F)
  caller_df <- select(caller_df, c('sample', 'method', 'copy_num')) %>%
    pivot_wider(names_from='method', values_from='copy_num')
  filter(caller_df, caller_16 == caller_78)$sample
})

df <- load(gene, method_short, keep_paralog=T, keep_b=T, keep_qual=T)
df <- filter(df, sample %in% keep_samples)
df <- extend_paralog(df, 2)
df$method <- sub('qm2_', '', df$method)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'))

obs <- rbind(
  get_observations(df, 'SMN1 + SMN2 Exons 1-6', 'copy_num_noise', 'b_copy_num', T),
  get_observations(df, 'SMN1 Exons 1-6', 'paralog1_noise', 'b_paralog1', T),
  get_observations(df, 'SMN2 Exons 1-6', 'paralog2_noise', 'b_paralog2', T)
)

(g3 <- ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(bquote(bold(.(method)) ~ 'CN estimate'),
                     breaks=0:20) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw())
ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=10, height=6)

mutate(obs, eq = round(a_obs) == round(b_obs)) %>%
  select(label, eq) %>% table

plot_grid(g0, g1, g2, g3, labels=LETTERS)
ggsave('~/Tmp/1.png', width=15, height=8)

# ====== AMY1C ======

gene <- 'AMY1C'
method <- 'qPCR'; method_short <- 'qpcr'
y_label <- sprintf('CN estimate', method)

df <- load(gene, method_short, keep_paralog=F, keep_b=T, keep_qual=T)
obs <- rbind(
  # get_observations(df, 'qPCR', 'copy_num', 'b_copy_num', method == 'qPCR'),
  get_observations(df, 'qPCR', 'copy_num', 'b_copy_num', method == 'qPCR_14'),
  get_observations(df, 'PRT', 'copy_num', 'b_copy_num', method == 'PRT'),
  get_observations(df, 'WGS', 'copy_num', 'b_copy_num', method == 'g1k')
)
obs <- filter(obs, !is.na(a_obs) & !is.na(b_obs))

obs_stat <- group_by(obs, label) %>%
  summarize(cor=cor(a_obs, b_obs),
            delta=sum(abs(a_obs - b_obs)) / length(a_obs),
            delta2=sum(abs(b_obs / a_obs - 1)) / length(a_obs))
obs_stat$text <- with(obs_stat, sprintf(' cor = %.3f\n Δ = %.3f', cor, delta))

ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.3) +
  geom_text(aes(-Inf, Inf, label=text), data=obs_stat, hjust=0, vjust=1.1) +
  scale_x_continuous('Parascopy CN estimate', breaks=seq(0, 20, 2)) +
  scale_y_continuous(y_label, breaks=seq(0, 20, 2)) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=11, height=6)

# ====== NPY4R ======

gene <- 'NPY4R'
method <- ''; method_short <- 'comp'
y_label <- 'CN estimate'

df <- load(gene, method_short, keep_paralog=F, keep_b=T, keep_qual=T)
df <- add_noise(df, 'copy_num', 0.015)

obs <- rbind(
  get_observations(df, 'FREEC', 'copy_num_noise', 'b_copy_num', method == 'FREEC'),
  get_observations(df, 'CNVnator', 'copy_num_noise', 'b_copy_num', method == 'CNVnator'),
  get_observations(df, 'ddPCR', 'copy_num_noise', 'b_copy_num', method == 'ddPCR')
)

ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.6) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=10, height=6)

# ====== FCGR3A ======

gene <- 'FCGR3A'
method <- ''; method_short <- 'comp'
y_label <- 'CN estimate'

df <- load(gene, method_short, keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df$b_copy_num <- as.numeric(df$b_copy_num)
df <- add_noise(df, c('copy_num', 'b_copy_num',
                      'paralog1', 'paralog2', 'b_paralog1', 'b_paralog2'))

obs <- rbind(
  get_observations(df, 'TaqMan', 'copy_num_noise', 'b_copy_num_noise',
                   method == 'TaqMan'),
  get_observations(df, 'PRT-REDVR', 'copy_num_noise', 'b_copy_num_noise',
                   method == 'PRT_REDVR'),
  get_observations(df, 'SYBR Green', 'copy_num_noise', 'b_copy_num_noise',
                   method == 'SYBR_Green')
)

ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.4) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw() +
  ggtitle('FCGR3A + FCGR3B')
ggsave(sprintf('%s/%s.a.png', plot_dir, gene), width=10, height=6)

df_filt <- filter(df, paralog_filter == 'PASS' & paralog_qual1 >= 30)
obs <- rbind(
  get_observations(df_filt, 'TaqMan', 'paralog2_noise', 'b_paralog2_noise',
                   method == 'TaqMan'),
  get_observations(df_filt, 'PRT-REDVR', 'paralog2_noise', 'b_paralog2_noise',
                   method == 'PRT_REDVR'),
  get_observations(df_filt, 'STR', 'paralog2_noise', 'b_paralog2_noise',
                   method == 'STR')
)

ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.4) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw() +
  ggtitle('FCGR3B')
ggsave(sprintf('%s/%s.b.png', plot_dir, gene), width=10, height=6)

all_match_cn <- filter(df, !is.na(b_copy_num)) %>%
  aggregate(b_copy_num ~ sample, ., function(x) length(unique(x)))
all_match_cn <- filter(all_match_cn, b_copy_num == 1)$sample

all_match_pcn <- filter(df, !is.na(b_paralog2)) %>%
  aggregate(b_paralog2 ~ sample, ., function(x) length(unique(x)))
all_match_pcn <- filter(all_match_pcn, b_paralog2 == 1)$sample
all_match_pcn <- intersect(all_match_cn, all_match_pcn)

obs <- rbind(
  get_observations(df, 'FCGR3A + FCGR3B', 'copy_num_noise', 'b_copy_num_noise',
                   method == 'TaqMan' & sample %in% all_match_cn),
  get_observations(df_filt, 'FCGR3B', 'paralog2_noise', 'b_paralog2_noise',
                   method == 'TaqMan' & sample %in% all_match_pcn))
ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.6) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw() +
  ggtitle('Three methods agree')
ggsave(sprintf('%s/%s.c.png', plot_dir, gene), width=8, height=6)

# ====== RHD/RHCE ======

gene <- 'RHCE'
method <- ''; method_short <- 'comp'
y_label <- 'CN estimate'

df <- load(gene, method_short, keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- add_noise(df, c('copy_num', 'b_copy_num',
                      'paralog1', 'paralog2', 'b_paralog1', 'b_paralog2'))

obs <- rbind(
  get_observations(df, 'RHD + RHCE', 'copy_num_noise', 'b_copy_num',
                   method == 'WGS'),
  get_observations(df, 'RHD', 'paralog2_noise', 'b_paralog2',
                   method == 'WGS'),
  get_observations(df, 'RHCE', 'paralog1_noise', 'b_paralog1',
                   method == 'WGS')
)
obs$label <- factor(obs$label, levels = unique(obs$label))

# WGS
ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.4) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw()
ggsave(sprintf('%s/%s.wgs.png', plot_dir, gene), width=10, height=6)

obs <- rbind(
  get_observations(df, 'RHD + RHCE', 'copy_num_noise', 'b_copy_num_noise',
                   method == 'MIP'),
  get_observations(df, 'RHD', 'paralog2_noise', 'b_paralog2_noise',
                   method == 'MIP'),
  get_observations(df, 'RHCE', 'paralog1_noise', 'b_paralog1_noise',
                   method == 'MIP')
)
obs$label <- factor(obs$label, levels = unique(obs$label))

# MIP
ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.4) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous(y_label, breaks=0:20) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw()
ggsave(sprintf('%s/%s.mip.png', plot_dir, gene), width=10, height=6)

# ====== SRGAP2 ======

gene <- 'SRGAP2'
method <- ''; method_short <- 'comp'
y_label <- 'CN estimate'

df <- load(gene, method_short, keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 4)
df <- add_noise(df, c('copy_num', 'b_copy_num',
                      'paralog1', 'paralog2', 'paralog3', 'paralog4',
                      'b_paralog1', 'b_paralog2', 'b_paralog3', 'b_paralog4'))

obs <- rbind(
  get_observations(df, 'Aggregate CN', 'copy_num_noise', 'b_copy_num',
                   method == 'WGS'),
  get_observations(df, 'SRGAP2A', 'paralog1_noise', 'b_paralog1',
                   method == 'WGS'),
  get_observations(df, 'SRGAP2B', 'paralog2_noise', 'b_paralog2',
                   method == 'WGS'),
  get_observations(df, 'SRGAP2C', 'paralog3_noise', 'b_paralog3',
                   method == 'WGS'),
  get_observations(df, 'SRGAP2D', 'paralog4_noise', 'b_paralog4',
                   method == 'WGS')
)
obs$label2 <- ifelse(grepl('SRGAP', obs$label), 'Paralog-specific CN', 'Total CN')
obs$label2 <- factor(obs$label2, unique(obs$label2))

ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs, color = label), alpha=.8) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('WGS CN estimate', breaks=0:20) +
  facet_wrap(. ~ label2, scales='free') +
  scale_color_manual('', values = c('gray20', RColorBrewer::brewer.pal(4, 'Set1'))) +
  theme_bw()
ggsave(sprintf('%s/%s.wgs.png', plot_dir, gene), width=10, height=6)

# ------ MIP ------

obs <- rbind(
  get_observations(df, 'Aggregate CN', 'copy_num_noise', 'b_copy_num_noise',
                   method == 'MIP'),
  get_observations(df, 'SRGAP2A', 'paralog1_noise', 'b_paralog1_noise',
                   method == 'MIP'),
  get_observations(df, 'SRGAP2B', 'paralog2_noise', 'b_paralog2_noise',
                   method == 'MIP'),
  get_observations(df, 'SRGAP2C', 'paralog3_noise', 'b_paralog3_noise',
                   method == 'MIP'),
  get_observations(df, 'SRGAP2D', 'paralog4_noise', 'b_paralog4_noise',
                   method == 'MIP')
)
obs$label2 <- ifelse(grepl('SRGAP', obs$label), 'Paralog-specific CN', 'Total CN')
obs$label2 <- factor(obs$label2, unique(obs$label2))

ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs, color = label), alpha=.8) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('MIP CN estimate', breaks=0:20) +
  facet_wrap(. ~ label2, scales='free') +
  scale_color_manual('', values = c('gray20', RColorBrewer::brewer.pal(4, 'Set1'))) +
  theme_bw()
ggsave(sprintf('%s/%s.mip.png', plot_dir, gene), width=10, height=6)

# ------ FISH ------

obs <- rbind(
  get_observations(df, 'Aggregate CN', 'copy_num_noise', 'b_copy_num_noise',
                   method == 'FISH'),
  get_observations(df, 'SRGAP2A', 'paralog1_noise', 'b_paralog1_noise',
                   method == 'FISH'),
  get_observations(df, 'SRGAP2B', 'paralog2_noise', 'b_paralog2_noise',
                   method == 'FISH'),
  get_observations(df, 'SRGAP2C', 'paralog3_noise', 'b_paralog3_noise',
                   method == 'FISH'),
  get_observations(df, 'SRGAP2D', 'paralog4_noise', 'b_paralog4_noise',
                   method == 'FISH')
)
obs$label2 <- ifelse(grepl('SRGAP', obs$label), 'Paralog-specific CN', 'Total CN')
obs$label2 <- factor(obs$label2, unique(obs$label2))

ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs, color = label), alpha=.8) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  scale_y_continuous('FISH CN estimate', breaks=0:20) +
  facet_wrap(. ~ label2, scales='free') +
  scale_color_manual('', values = c('gray20', RColorBrewer::brewer.pal(4, 'Set1'))) +
  theme_bw()
ggsave(sprintf('%s/%s.fish.png', plot_dir, gene), width=10, height=6)
