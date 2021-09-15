library(cowplot)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

data_dir <- '~/Data/hg38/jvc/comparisons/genes'
plot_dir <- '~/Data/hg38/jvc/plots/population_comparison/v100/experimental'

abline_color <- 'gray80'
abline_size <- 1

QUAL <- 20

# ====== SMN1 ======

# ------ SMNCopyNumberCaller ------

gene <- 'SMN1'
method <- 'SMNCopyNumberCaller'; method_short <- 'caller'
y_label <- sprintf('%s CN estimate', method)
df <- load(gene, method_short, keep_qual=T, keep_b=T)
df <- extend_paralog(df, 2)
df$method <- sub('caller_', '', df$method)
df <- add_noise(df, c('agCN', 'psCN1', 'psCN2', 'b_paralog1', 'b_paralog2'))

obs <- rbind(
  get_observations(df, 'SMN1 + SMN2\nexons 1-6', 'agCN_noise', 'b_copy_num',
                   method == '16' & agCN_qual >= 30),
  get_observations(df, 'SMN1\nexons 7-8', 'psCN1_noise', 'b_paralog1_noise',
                   method == '78' & psCN_filter == 'PASS' & psCN_qual1 >= 30),
  get_observations(df, 'SMN2\nexons 7-8', 'psCN2_noise', 'b_paralog2_noise',
                   method == '78' & psCN_filter == 'PASS' & psCN_qual1 >= 30)
)
obs$label <- factor(obs$label, levels=unique(obs$label))

(g_smn_c <- ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=0.6, size=0.5) +
  scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
  # scale_y_continuous(bquote(bold(.(method)) ~ 'CN estimate'),
  scale_y_continuous('SMNCopyNumberCaller CN estimate',
                     breaks=0:20, limits=c(NA, max(obs$b_obs, na.rm=T) + 0.5)) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw() +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(margin=margin(2, 0, 2, 0))))
filter(obs, round(a_obs) != round(b_obs))

ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=12, height=6)

cn_match(df, method == '16')
cn_match(df, method == '78')
par_cn_match(df, 2, method == '78')

# ------ QM2 ------

gene <- 'SMN1'
method <- 'QuicK-mer2'; method_short <- 'qm2'
y_label <- sprintf('%s CN estimate', method)
df <- load(gene, method_short, keep_qual=T, keep_b=T)
df <- extend_paralog(df, 2)
df$region <- sub('SMN1_', '', df$region)
df <- add_noise(df, c('agCN', 'psCN1', 'psCN2'))

length(unique(df$sample))

obs <- rbind(
  get_observations(df, 'SMN1 + SMN2\nexons 1-6', 'agCN_noise', 'b_copy_num',
                   region == 'exon2' & agCN_qual >= QUAL),
  get_observations(df, 'SMN1\nexons 7-8', 'psCN1_noise', 'b_paralog1',
                   region == 'exon7' & psCN_filter == 'PASS' & psCN_qual1 >= QUAL),
  get_observations(df, 'SMN2\nexons 7-8', 'psCN2_noise', 'b_paralog2',
                   region == 'exon7' & psCN_filter == 'PASS' & psCN_qual1 >= QUAL)
)
obs$label <- factor(obs$label, levels=unique(obs$label))

(g_smn_e <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=0.6, size=0.5) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous(y_label,
                       breaks=0:20, limits=c(NA, max(obs$b_obs, na.rm=T) + 0.5)) +
    facet_grid(. ~ label, scales='free', space='free') +
    theme_bw() +
    theme(axis.title.y = element_text(size=10),
          strip.text = element_text(margin=margin(2, 0, 2, 0))))
filter(obs, round(a_obs) != round(b_obs))

plot_grid(g_smn_c, g_smn_e, ncol=1, labels=LETTERS)
ggsave('~/Tmp/1.png', width=10, height=10, scale=.7, dpi=450)

ggsave(sprintf('%s/%s.%s.png', plot_dir, gene, method_short), width=12, height=6)

cn_match(df, method == '16')
cn_match(df, method == '78')
par_cn_match(df, 2, method == '78')


# ------ MLPA ------

gene <- 'SMN1'
method <- 'MLPA'; method_short <- 'mlpa'
y_label <- sprintf('%s CN estimate', method)

# df <- load(gene, 'old/mlpa', keep_paralog=F, keep_b=T, keep_qual=T)
df <- load(gene, method_short, keep_paralog=F, keep_b=T, keep_qual=T)
# For some reasons, some samples have two values in MLPA.
df <- group_by(df, sample, method) %>% slice_head(n=1) %>% ungroup
df$method <- sub('mlpa_', '', df$method)
df <- add_noise(df, c('agCN'), 0.02)

obs <- rbind(
  get_observations(df, 'Exons 1-6', 'agCN_noise', 'b_copy_num',
                   method == '16' & agCN_qual >= QUAL),
  get_observations(df, 'Exons 7-8', 'agCN_noise', 'b_copy_num',
                   method == '78' & agCN_qual >= QUAL))
obs$label <- factor(obs$label, levels=unique(obs$label))

(g_smn_d <- ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=0.6, size=0.5) +
  scale_x_continuous('Parascopy agCN estimate', breaks=0:20) +
  # scale_y_continuous(bquote(bold(.(method)) ~ 'CN estimate'),
  scale_y_continuous('MLPA agCN estimate', breaks=0:20,
    expand=expansion(add=0.7)) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw() +
  theme(strip.text = element_text(margin=margin(2, 0, 2, 0))))
# sum(with(obs, !is.na(a_obs) & !is.na(b_obs)))
ggsave('~/Tmp/1.png', width=6, height=3, scale=.85, dpi=600)

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

plot_grid(g_smn_a, g_smn_b, g_smn_c, g_smn_d, labels=LETTERS,
          rel_heights=c(0.6, 0.4), label_y=1.01)
ggsave('~/Tmp/1.png', width=15, height=8)

plot_grid(g_smn_a, g_smn_c, g_smn_b, g_smn_d, labels=c('A', 'C', 'B', 'D'),
          rel_widths=c(0.7, 0.3), label_y=1.01)
ggsave('~/Tmp/1.png', width=15, height=8, scale=0.9)

#####
plot_grid(plot_grid(g_smn_a, g_smn_c, g_smn_d, labels=c('A', 'B', 'C'), nrow=1,
                    rel_widths=c(1, 0.9, 0.9)),
          g_smn_b, labels=c('', 'D'), rel_heights=c(0.5, 0.5), nrow=2)
ggsave('~/Tmp/1.png', width=15, height=9, scale=0.7)
#####

plot_grid(plot_grid(g_smn_a, g_smn_c, g_smn_d, labels=c('A', 'B', 'C'),
                    ncol=1, rel_heights = c(1.3, 0.8, 0.8), label_y=1.015),
          g_smn_b, labels=c('', 'D'), rel_widths=c(0.45, 0.55), ncol=2, label_y=1.007)
ggsave('~/Tmp/1.png', width=15, height=8, scale=1)

cn_match(df, method == '16' & population != 'BEB' & population != 'TSI')
cn_match(df, method == '78' & population != 'BEB' & population != 'TSI')

plot_grid(plot_grid(g_viterbi, g_smn_d, labels=c('A', 'B'), nrow=1,
                    rel_widths=c(1.2, 1)),
          g_smn_b, labels=c('', 'C'), rel_heights=c(0.5, 0.5), nrow=2,
          label_y = 1.05)
# ggsave('~/Tmp/1.png', width=15, height=9, scale=0.7)
ggsave('~/Tmp/1.png', width=15, height=9, scale=0.55, dpi=600)

# ====== AMY1C ======

gene <- 'AMY1C'
method <- 'qPCR'; method_short <- 'qpcr'
y_label <- sprintf('CN estimate', method)

df <- load(gene, method_short, keep_paralog=F, keep_b=T, keep_qual=T)
df <- filter(df, agCN_qual >= QUAL)
obs <- rbind(
  # get_observations(df, 'qPCR', 'agCN', 'b_copy_num', method == 'qPCR'),
  get_observations(df, 'qPCR', 'agCN', 'b_copy_num', method == 'qPCR_14'),
  get_observations(df, 'PRT', 'agCN', 'b_copy_num', method == 'PRT'),
  get_observations(df, 'WGS', 'agCN', 'b_copy_num', method == 'g1k')
)
obs <- filter(obs, !is.na(a_obs) & !is.na(b_obs))

obs_stat <- group_by(obs, label) %>%
  summarize(cor=cor(a_obs, b_obs),
            delta=sum(abs(a_obs - b_obs)) / length(a_obs),
            delta2=sum(abs(b_obs / a_obs - 1)) / length(a_obs))
obs_stat$text <- with(obs_stat, sprintf(' r  = %.3f\n Δ = %.3f', cor, delta))

ggplot(obs) +
  geom_abline(color = abline_color, size = abline_size) +
  geom_point(aes(a_obs, b_obs), alpha=.5, size=0.9) +
  geom_text(aes(-Inf, Inf, label=text), data=obs_stat, hjust=0, vjust=1.1) +
  scale_x_continuous('Parascopy agCN estimate', breaks=seq(0, 20, 2)) +
  scale_y_continuous('agCN estimate', breaks=seq(0, 20, 2)) +
  facet_grid(. ~ label, scales='free', space='free') +
  theme_bw()
ggsave('~/Tmp/1.png', width=11, height=6, scale=.6, dpi=450)

length(unique(df$sample))
table(unique(df[c('sample', 'superpopulation')])$superpopulation)

filter(df, agCN == 18)

# cn_match(df, method == 'g1k')

# ====== NPY4R ======

gene <- 'NPY4R'
method <- ''; method_short <- 'comp'
y_label <- 'CN estimate'

df <- load(gene, method_short, keep_paralog=F, keep_b=T, keep_qual=T)
df <- add_noise(df, 'agCN', 0.015)

obs <- rbind(
  get_observations(df, 'FREEC', 'agCN_noise', 'b_copy_num', method == 'FREEC'),
  get_observations(df, 'CNVnator', 'agCN_noise', 'b_copy_num', method == 'CNVnator'),
  get_observations(df, 'ddPCR', 'agCN_noise', 'b_copy_num', method == 'ddPCR')
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
df <- add_noise(df, c('agCN', 'b_copy_num',
                      'psCN1', 'psCN2', 'b_paralog1', 'b_paralog2'))

# cn_match(df, m)

obs <- rbind(
  get_observations(df, 'TaqMan', 'agCN_noise', 'b_copy_num_noise',
                   method == 'TaqMan'),
  get_observations(df, 'PRT-REDVR', 'agCN_noise', 'b_copy_num_noise',
                   method == 'PRT_REDVR'),
  get_observations(df, 'SYBR Green', 'agCN_noise', 'b_copy_num_noise',
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

df_filt <- filter(df, psCN_filter == 'PASS' & psCN_qual1 >= 30)
obs <- rbind(
  get_observations(df_filt, 'TaqMan', 'psCN2_noise', 'b_paralog2_noise',
                   method == 'TaqMan'),
  get_observations(df_filt, 'PRT-REDVR', 'psCN2_noise', 'b_paralog2_noise',
                   method == 'PRT_REDVR'),
  get_observations(df_filt, 'STR', 'psCN2_noise', 'b_paralog2_noise',
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
  get_observations(df, 'FCGR3A + FCGR3B', 'agCN_noise', 'b_copy_num_noise',
                   method == 'TaqMan' & sample %in% all_match_cn),
  get_observations(df_filt, 'FCGR3B', 'psCN2_noise', 'b_paralog2_noise',
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
df <- add_noise(df, c('agCN', 'b_copy_num',
                      'psCN1', 'psCN2', 'b_paralog1', 'b_paralog2'))

obs <- rbind(
  get_observations(df, 'RHD + RHCE', 'agCN_noise', 'b_copy_num',
                   method == 'WGS'),
  get_observations(df, 'RHD', 'psCN2_noise', 'b_paralog2',
                   method == 'WGS'),
  get_observations(df, 'RHCE', 'psCN1_noise', 'b_paralog1',
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
  get_observations(df, 'RHD + RHCE', 'agCN_noise', 'b_copy_num_noise',
                   method == 'MIP'),
  get_observations(df, 'RHD', 'psCN2_noise', 'b_paralog2_noise',
                   method == 'MIP'),
  get_observations(df, 'RHCE', 'psCN1_noise', 'b_paralog1_noise',
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

cn_match(df, method == 'MIP')
par_cn_match(df, 2, method == 'MIP')

cn_match(df, method == 'WGS')
par_cn_match(df, 2, method == 'WGS')

# ====== SRGAP2 ======

gene <- 'SRGAP2'
method <- ''; method_short <- 'comp'
y_label <- 'CN estimate'

df <- load(gene, method_short, keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 4)
df <- add_noise(df, c('agCN', 'b_copy_num',
                      'psCN1', 'psCN2', 'psCN3', 'psCN4',
                      'b_paralog1', 'b_paralog2', 'b_paralog3', 'b_paralog4'))

obs <- rbind(
  get_observations(df, 'Aggregate CN', 'agCN_noise', 'b_copy_num',
                   method == 'WGS'),
  get_observations(df, 'SRGAP2A', 'psCN1_noise', 'b_paralog1',
                   method == 'WGS'),
  get_observations(df, 'SRGAP2B', 'psCN2_noise', 'b_paralog2',
                   method == 'WGS'),
  get_observations(df, 'SRGAP2C', 'psCN3_noise', 'b_paralog3',
                   method == 'WGS'),
  get_observations(df, 'SRGAP2D', 'psCN4_noise', 'b_paralog4',
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

cn_match(df, method == 'WGS')
par_cn_match(df, 4, method == 'WGS')

# ------ MIP ------

obs <- rbind(
  get_observations(df, 'Aggregate CN', 'agCN_noise', 'b_copy_num_noise',
                   method == 'MIP'),
  get_observations(df, 'SRGAP2A', 'psCN1_noise', 'b_paralog1_noise',
                   method == 'MIP'),
  get_observations(df, 'SRGAP2B', 'psCN2_noise', 'b_paralog2_noise',
                   method == 'MIP'),
  get_observations(df, 'SRGAP2C', 'psCN3_noise', 'b_paralog3_noise',
                   method == 'MIP'),
  get_observations(df, 'SRGAP2D', 'psCN4_noise', 'b_paralog4_noise',
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

cn_match(df, method == 'MIP')
par_cn_match(df, 4, method == 'MIP')

# ------ FISH ------

obs <- rbind(
  get_observations(df, 'Aggregate CN', 'agCN_noise', 'b_copy_num_noise',
                   method == 'FISH'),
  get_observations(df, 'SRGAP2A', 'psCN1_noise', 'b_paralog1_noise',
                   method == 'FISH'),
  get_observations(df, 'SRGAP2B', 'psCN2_noise', 'b_paralog2_noise',
                   method == 'FISH'),
  get_observations(df, 'SRGAP2C', 'psCN3_noise', 'b_paralog3_noise',
                   method == 'FISH'),
  get_observations(df, 'SRGAP2D', 'psCN4_noise', 'b_paralog4_noise',
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

cn_match(df, method == 'FISH')
par_cn_match(df, 4, method == 'FISH')
