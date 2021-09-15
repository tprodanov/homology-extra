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
               psCN_filter == 'PASS' & psCN_qual1 >= QUAL & psCN_qual2 >= QUAL)
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

(plots[[gene]] <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    ggtitle('SMN1 (exons 7-8): 13 / 17 reliable PSVs') +
    facet_grid(. ~ label, scales='free', space='free') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=.5, size=11, margin = margin(3, 0, 3, 0)),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))

# cat em_f_values.csv | sed 's/chr5://' | awk '$2 >= 70948120 && $2 <= 70953015' | awk '$6 >= 0.95 && $7 >= 0.95' | wc -l

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
  get_observations(df, 'ABCC6 + p1 + p2', 'agCN_noise', 'b_copy_num', T),
  get_observations(df, 'ABCC6', 'psCN1_noise', 'b_paralog1', T),
  get_observations(df, 'ABCC6p1', 'psCN3_noise', 'b_paralog3', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(plots[[gene]] <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20,
                       expand=expansion(add=0.2)) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, space='free', scales='free') +
    ggtitle('ABCC6: 355 / 399 reliable PSVs') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=.5, size=11, margin = margin(3, 0, 3, 0)),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))
# cat em_f_values.csv | sed 's/chr.*://' | awk '$2 <= 16223494' | awk '$6 >= 0.95 && $7 >= 0.95 && (NF == 7 || $8 >= 0.95)' | wc -l

# -------------------- FCGR3A

gene <- 'FCGR3A'
df <- load(gene, 'qm2', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, agCN_filter == 'PASS' & agCN_qual >= QUAL &
               psCN_filter == 'PASS' & psCN_qual1 >= QUAL & psCN_qual2 >= QUAL)
df <- add_noise(df, c('agCN', 'psCN1', 'psCN2'))

obs <- rbind(
  get_observations(df, 'FCGR3A + FCGR3B', 'agCN_noise', 'b_copy_num', T),
  get_observations(df, 'FCGR3A', 'psCN1_noise', 'b_paralog1', T),
  get_observations(df, 'FCGR3B', 'psCN2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(plots[[gene]] <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20, expand=expansion(add=.3)) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    ggtitle('FCGR3A: 122 / 179 reliable PSVs') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=.5, size=11, margin = margin(3, 0, 3, 0)),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))

# -------------------- ZP3

gene <- 'ZP3'
df <- load(gene, 'qm2-matrix', keep_paralog=T, keep_b=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, agCN_filter == 'PASS' & agCN_qual >= QUAL &
               psCN_filter == 'PASS' & psCN_qual1 >= QUAL & psCN_qual2 >= QUAL)
df <- add_noise(df, c('agCN', 'psCN1', 'psCN2'))

obs <- rbind(
  get_observations(df, 'ZP3 + POMZP3', 'agCN_noise', 'b_copy_num', T),
  get_observations(df, 'ZP3', 'psCN1_noise', 'b_paralog1', T),
  get_observations(df, 'POMZP3', 'psCN2_noise', 'b_paralog2', T)
)
obs$label <- factor(obs$label, levels = unique(obs$label))

(plots[[gene]] <- ggplot(obs) +
    geom_abline(color = abline_color, size = abline_size) +
    geom_point(aes(a_obs, b_obs), alpha=.3) +
    scale_x_continuous('Parascopy CN estimate', breaks=0:20, expand=expansion(add=.3)) +
    scale_y_continuous('QuicK-mer2 CN estimate', breaks=0:20) +
    facet_grid( ~ label, scales='free', space='free') +
    ggtitle('ZP3: 18 / 46 reliable PSVs') +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=.5, size=11, margin = margin(3, 0, 3, 0)),
      plot.margin = unit(c(3, 3, 3, 3), "pt")))
# cat em_f_values.csv | sed 's/chr.*://' | awk '$2 <= 76442071' | awk '$6 >= 0.95 && $7 >= 0.95' | wc -l

# --------- Merge

plot_grid(plots[['ABCC6']], plots[['FCGR3A']],
          plots[['SMN1']], plots[['ZP3']],
          nrow = 2, labels = LETTERS)
ggsave('~/Tmp/1.png', width=10, height=9, scale=.7, dpi=450)
