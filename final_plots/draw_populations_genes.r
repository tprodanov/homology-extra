library(ggh4x)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

data_dir <- '~/Data/hg38/jvc/comparisons/populations/r009'
plot_dir <- '~/Data/hg38/jvc/plots/population_comparison/r009/genes'

colors <- ggthemes::tableau_color_pal()(10)[c(1, 2, 3, 5, 7)]
border_color <- 'black'
border_size <- 0.3
bar_width <- 0.7

str_samples <- function(count) {
  sprintf('%d sample%s', count, ifelse(count %% 10 == 1 & count != 11, '', 's'))
}

count_events <- function(df, col) {
  pop_counts <- table(df$superpopulation)
  event_counts <- rev(sort(table(df[[col]])))
  for (event in names(event_counts)) {
    cat(sprintf('%s have %s:\n    ', str_samples(event_counts[event]), event))
    curr_counts <- table(df[df[[col]] == event,]$superpopulation)
    for (i in seq_along(curr_counts)) {
      pop <- names(curr_counts)[i]
      cat(sprintf('%s%s: %d / %d', ifelse(i > 1, ',  ', ''), pop,
                  curr_counts[i], pop_counts[pop]))
    }
    cat('\n')
  }
}

# ====== SMN1 ======

gene <- 'SMN1'
df <- load(gene, 'from_all')
df <- filter(df, grepl(gene, region))
df <- extend_paralog(df, 2)
df$region <- sub('SMN1_', '', df$region)
df_wide <- to_wide(df)
df_wide <- filter_df(df_wide, !is.na(paralog1_end) & !is.na(paralog2_end))

df_wide <- mutate(df_wide,
                  cn1 = copy_num_middle,
                  pcn1 = paralog1_end,
                  pcn2 = paralog2_end) %>%
  mutate(cn1=ifelse(cn1 >= 5, '≥ 5', as.character(cn1)),
         pcn1=as.character(pcn1), pcn2=as.character(pcn2))

region_names <- c('cn1' = 'SMN1 + SMN2 exons 1-6',
                  'pcn1' = 'SMN1 exons 7-8',
                  'pcn2' = 'SMN2 exons 7-8')
cn_counts <- to_counts(df_wide, region_names, c('cn1', 'pcn1', 'pcn2'))
cn_counts <- mutate(cn_counts,
                    value=factor(value, levels=c(as.character(0:4), '≥ 5')))
ref_cns <- ref_cns_df(cn_counts, c(4, 2, 2))

labels1 <- c('SMN1 + SMN2', 'SMN1', 'SMN2')
labels2 <- c('Exons 1-6', 'Exons 7-8', 'Exons 7-8')
cn_counts <- mutate(cn_counts,
                    label1 = labels1[as.numeric(region)],
                    label2 = labels2[as.numeric(region)])
ref_cns <- mutate(ref_cns,
                    label1 = labels1[as.numeric(region)],
                    label2 = labels2[as.numeric(region)])

(g_smn_a <- ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=0.95, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0),
           width=0.8,
           color=border_color, size=border_size) +
  facet_nested(label2 + label1 ~ .) +
  scale_x_discrete('Parascopy CN estimate') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw() +
  theme(legend.position='top', legend.margin = margin(-2, 0, -5, 0),
        legend.key.size = unit(10, 'pt')))
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

mutate(df_wide, 'cn_all' = sprintf('%d,%d,%d', cn1, cn2, cn2_delta)) %>%
  select(c('superpopulation', 'cn_all')) %>%
  table() %>% t()

# ====== SERF1A ======

gene <- 'SERF1A'
df <- load('SMN1', 'from_all', keep_paralog=F)
df <- filter(df, grepl(gene, region))
df$region <- sub('SERF1A_', '', df$region)
df_wide <- to_wide(df)
df_wide$n_del <- df_wide$middle - df_wide$end
df_wide <- filter_df(df_wide, n_del >= 0)

region_names <- c('middle' = 'Exon 1',
                  'end' = 'Exon 3')
cn_counts <- to_counts(df_wide, region_names, c('middle', 'end'))
ref_cns <- ref_cns_df(cn_counts, c(4, 4))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s_v1.png', plot_dir, gene), width=8, height=6)

region_names <- c('middle' = 'Exon 1',
                  'n_del' = 'Exon 3 deletion')
cn_counts <- to_counts(df_wide, region_names, c('middle', 'n_del'))
ref_cns <- ref_cns_df(cn_counts, c(4, 0))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s_v2.png', plot_dir, gene), width=8, height=6)

# ====== C4A =======

gene <- 'C4A'
df <- load(gene, 'from_pops', keep_paralog=F)
df_wide <- to_wide(df)
df_wide$n_del <- df_wide$end - df_wide$middle
df_wide <- filter_df(df_wide, n_del >= 0 & start == end)

region_names <- c('start' = 'Exon 1',
                  #'middle' = 'Intron 9',
                  #'end' = 'Exon 22',
                  'n_del' = 'Intron 9 deletion')
cn_counts <- to_counts(df_wide, region_names, c('start', 'n_del'))
ref_cns <- ref_cns_df(cn_counts, c(4, 0))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== FAM185A ======

gene <- 'FAM185A'
df <- load(gene, 'from_pops', keep_paralog=F)
df_wide <- to_wide(df)
df_wide$n_del <- df_wide$after - df_wide$del
df_wide <- filter_df(df_wide, n_del >= 0 & before == after)

region_names <- c('before' = 'FAM185A',
                  'n_del' = 'Deletion')
cn_counts <- to_counts(df_wide, region_names, c('before', 'n_del'))
ref_cns <- ref_cns_df(cn_counts, c(4, 0))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== NBPF4 ======

gene <- 'NBPF4'
df <- load(gene, 'from_pops', keep_paralog=F)
df_wide <- to_wide(df)

region_names <- c('start' = 'Exons 10-14',
                  'middle' = 'Exon 5',
                  'end' = 'Exon 1')
cn_counts <- to_counts(df_wide, region_names, c('start', 'middle', 'end'))
ref_cns <- ref_cns_df(cn_counts, c(6, 6, 8))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== EIF3C ======

gene <- 'EIF3C'
df <- load(gene, 'from_pops', keep_paralog=F)
df_wide <- to_wide(df)
df_wide <- filter_df(df_wide, !is.na(start))

region_names <- c('start' = 'Exons 1',
                  'middle' = 'Intron 1',
                  'end' = 'Exon 9')
cn_counts <- to_counts(df_wide, region_names, c('start', 'middle', 'end'))
ref_cns <- ref_cns_df(cn_counts, c(6, 4, 4))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== SRGAP2 ======

gene <- 'SRGAP2'
df <- load(gene, 'from_pops', keep_paralog=T)
df <- filter(df, region == 'd' | region == 'e')
df_wide <- to_wide(df)
df_wide <- mutate(df_wide,
    paralog_copy_num_d = reorder_copies(paralog_copy_num_d, c(1, 3, 2)),
    paralog_copy_num_e = reorder_copies(paralog_copy_num_e, c(1, 4, 2, 3)))

region_names <- c('d' = '3-copy region',
                  'e' = '4-copy region')
cn_counts <- to_counts(df_wide, region_names, matches('^paralog.*[de]'),
                       names_prefix='paralog_copy_num_')
ref_cns <- ref_cns_df(cn_counts, c('2,2,2', '2,2,2,2'))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_wrap(region ~ ., scales='free', strip.position='right', ncol=1) +
  scale_x_discrete('Paralog-specific copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s_v1.png', plot_dir, gene), width=11, height=6)

sum_perc <- aggregate(perc ~ value, cn_counts, sum)
sum_perc <- mutate(sum_perc,
                   value2 = ifelse(perc >= 2, as.character(value), 'Other'))
sum_perc$value2 <- factor(sum_perc$value2)

df_wide$clip_paralog_copy_num_d <- (sum_perc[
  match(df_wide$paralog_copy_num_d, sum_perc$value),]$value2)
df_wide$clip_paralog_copy_num_e <- (sum_perc[
  match(df_wide$paralog_copy_num_e, sum_perc$value),]$value2)
cn_counts2 <- to_counts(df_wide, region_names, matches('^clip_paralog.*[de]'),
                       names_prefix='clip_paralog_copy_num_')

ggplot(cn_counts2) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_wrap(region ~ ., scales='free', strip.position='right', ncol=1) +
  scale_x_discrete('Paralog-specific copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s_v2.png', plot_dir, gene), width=8, height=6)

# ====== STRC ======

gene <- 'STRC'
df <- load(gene, 'from_pops', keep_paralog=T)
df <- filter(df, region == 'gene_start') %>% mutate(region = '0')
df <- extend_paralog(df, 2)
df_wide <- to_wide(df)
df_wide <- rename(df_wide, 'cn1' = 'paralog1_0', 'cn2' = 'paralog2_0')
df_wide <- filter_df(df_wide, !is.na(cn1) & !is.na(cn2))

region_names <- c('cn1' = 'STRC', 'cn2' = 'pSTRC')
cn_counts <- to_counts(df_wide, region_names, starts_with('cn'))
ref_cns <- ref_cns_df(cn_counts, c(2, 2))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== NEB ======

gene <- 'NEB'
df <- load(gene, 'from_pops', keep_paralog=T)
df <- extend_paralog(df, 3)
df_wide <- to_wide(df)
df_wide <- rename(df_wide,
    'cn1' = 'paralog1_middle', 'cn2' = 'paralog2_middle', 'cn3' = 'paralog3_middle')
df_wide <- filter_df(df_wide, !is.na(cn1) & !is.na(cn2) & !is.na(cn3))

region_names <- c('cn1' = 'Copy 1', 'cn2' = 'Copy 2', 'cn3' = 'Copy 3')
cn_counts <- to_counts(df_wide, region_names, starts_with('cn'))
ref_cns <- ref_cns_df(cn_counts, c(2, 2, 2))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== RHCE ======

gene <- 'RHCE'
df <- load(gene, 'from_pops', keep_paralog=T)
df <- extend_paralog(df, 2)
df_wide <- to_wide(df)
df_wide <- rename(df_wide,
                  'cn1' = 'paralog1_0', 'cn2' = 'paralog2_0')
df_wide <- filter_df(df_wide, !is.na(cn1) & !is.na(cn2))

region_names <- c('cn1' = 'RHCE', 'cn2' = 'RHD')
cn_counts <- to_counts(df_wide, region_names, starts_with('cn'))
ref_cns <- ref_cns_df(cn_counts, c(2, 2))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== NPY4R ======

gene <- 'NPY4R'
df <- load(gene, 'from_pops', keep_paralog=F)
df$region <- 'cn'
df_wide <- to_wide(df)

region_names <- c('cn' = gene)
cn_counts <- to_counts(df_wide, region_names, starts_with('cn'))
ref_cns <- ref_cns_df(cn_counts, c(4))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== AMY1C ======

gene <- 'AMY1C'
df <- load(gene, 'from_pops', keep_paralog=F)
df$region <- 'cn'
df_wide <- to_wide(df)
df_wide <- mutate(df_wide,
                  cn_char = ifelse(cn >= 13, '≥ 13', as.character(cn)))

region_names <- c('cn_char' = gene)
cn_counts <- to_counts(df_wide, region_names, 'cn_char')
cn_counts$value <- factor(cn_counts$value, levels=c(as.character(1:12), '≥ 13'))
ref_cns <- ref_cns_df(cn_counts, c(6))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== FCGR3A ======

gene <- 'FCGR3A'
df <- load(gene, 'from_pops', keep_paralog=T)
df <- extend_paralog(df, 2)
df_wide <- to_wide(df)
df_wide <- rename(df_wide,
                  'cn1' = 'paralog1_0', 'cn2' = 'paralog2_0')
df_wide <- filter_df(df_wide, !is.na(cn1) & !is.na(cn2))

region_names <- c('cn1' = 'FCGR3A', 'cn2' = 'FCGR3B')
cn_counts <- to_counts(df_wide, region_names, starts_with('cn'))
ref_cns <- ref_cns_df(cn_counts, c(2, 2))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== OTOA ======

gene <- 'OTOA'
df <- load(gene, 'from_pops', keep_paralog=F)
df$region <- 'cn'
df_wide <- to_wide(df)

region_names <- c('cn' = gene)
cn_counts <- to_counts(df_wide, region_names, 'cn')
ref_cns <- ref_cns_df(cn_counts, c(4))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== ABCC6 ======

gene <- 'ABCC6'
df <- load(gene, 'from_pops', keep_paralog=F)
df_wide <- to_wide(df)

region_names <- c('02-01' = '2-copy region',
                  '03-01' = '3-copy region')
cn_counts <- to_counts(df_wide, region_names, c('02-01', '03-01'))
ref_cns <- ref_cns_df(cn_counts, c(4, 6))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

# ====== PMS2 ======

gene <- 'PMS2'
df <- load(gene, 'from_pops', keep_paralog=F)
df_wide <- to_wide(df)

region_names <- c('exon15' = 'Exon 15', 'exon14' = 'Exon 14')
cn_counts <- to_counts(df_wide, region_names, starts_with('exon'))
ref_cns <- ref_cns_df(cn_counts, c(4, 4))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  coord_cartesian(ylim=c(0, 10)) +
  scale_fill_manual('', values = colors) +
  theme_bw() +
  theme(axis.text.y=element_text(face=c(rep('plain', 4), 'bold')))
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

cn_counts
filter(cn_counts, key == 'exon14') %>% aggregate(n ~ superpopulation, ., sum)

# ====== HYDIN ======

gene <- 'HYDIN'
df <- load(gene, 'from_pops', keep_paralog=T)
df <- extend_paralog(df, 2)
df_wide <- to_wide(df)

region_names <- c('paralog1_exon18' = 'Copy 1', 'paralog2_exon18' = 'Copy 2')
cn_counts <- to_counts(df_wide, region_names, starts_with('paralog'))
ref_cns <- ref_cns_df(cn_counts, c(2, 2))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  coord_cartesian(ylim=c(0, 5)) +
  scale_fill_manual('', values = colors) +
  theme_bw() +
  theme(axis.text.y=element_text(face=c(rep('plain', 5), 'bold')))
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

df_wide %>% filter(paralog1_exon18 != 2 | paralog2_exon18 != 2) %>%
  count(superpopulation, paralog1_exon18, paralog2_exon18)
cn_counts

# ====== CEL ======

gene <- 'CEL'
df <- load(gene, 'from_pops', keep_paralog=T, keep_qual=T)
df <- extend_paralog(df, 2)
df <- filter(df, paralog_filter == 'PASS' & paralog_qual1 >= 30)

df_wide <- to_wide(df)
df_wide <- rename(df_wide, 'cn1' = 'paralog1_exon9', 'cn2' = 'paralog2_exon9')

region_names <- c('cn1' = 'CEL', 'cn2' = 'CELP')
cn_counts <- to_counts(df_wide, region_names, starts_with('cn'))
ref_cns <- ref_cns_df(cn_counts, c(2, 2))

ggplot(cn_counts) +
  geom_blank(aes(value)) +
  geom_tile(aes(cn, 0), width=bar_width + 0.1, height=Inf, data=ref_cns, alpha=.2) +
  geom_bar(aes(value, perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0), width=bar_width,
           color=border_color, size=border_size) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_manual('', values = colors) +
  theme_bw()
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

df_wide$event <- sprintf('%d,%d', df_wide$cn1, df_wide$cn2)
count_events(df_wide, 'event')
