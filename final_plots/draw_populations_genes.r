library(ggh4x)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))

data_dir <- '~/Data/hg38/jvc/comparisons/populations/r009'
plot_dir <- '~/Data/hg38/jvc/plots/population_comparison/r009/genes'

# colors <- ggthemes::tableau_color_pal()(10)[c(1, 2, 3, 5, 7)]
colors <- ggthemes::tableau_color_pal()(10)[c(1, 2, 3, 5)]
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

# labels1 <- c('SMN1 + 2', 'SMN1', 'SMN2')
labels1 <- c('SMN1 + 2\nexons 1-6', 'SMN1\nexons 7-8', 'SMN2\nexons 7-8')
labels1 <- factor(labels1, levels=labels1)
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
  # facet_nested(label2 + label1 ~ .) +
  facet_wrap(~ label1, ncol=1, strip.position='right') +
  scale_x_discrete('Parascopy CN estimate') +
  scale_y_continuous('Percentage of samples') +
  scale_fill_manual('', values = colors) +
  theme_bw() +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position='top', legend.margin = margin(-2, 0, -7, 0),
        legend.spacing.y = unit(3, 'pt'),
        #legend.text = element_text(margin=margin(-2, 0, -2, 0)),
        legend.key.size = unit(10, 'pt'),
        strip.text = element_text(margin = margin(0, 2, 0, 2))))
ggsave(sprintf('%s/%s.png', plot_dir, gene), width=8, height=6)

mutate(df_wide, 'cn_all' = sprintf('%d,%d,%d', cn1, cn2, cn2_delta)) %>%
  select(c('superpopulation', 'cn_all')) %>%
  table() %>% t()


df <- load(gene, 'from_all') %>% filter(region == 'SMN1_middle') %>% extend_paralog(2)
df <- filter(df, !is.na(paralog1))

paralog_counts <- count(df, copy_num, paralog1, paralog2)
paralog_counts$perc <- paralog_counts$n / nrow(df) * 100
annot <- data.frame()
for (cn in unique(paralog_counts$copy_num)) {
  label <- sprintf('Total: %s', (sum(filter(paralog_counts, copy_num == cn)$n)))
  annot <- rbind(annot, data.frame(copy_num=cn,
                                   perc=sum(filter(paralog_counts, copy_num == cn)$perc),
                                   label=label))
}

ggplot(paralog_counts) +
  geom_bar(aes(copy_num, perc, fill=factor(paralog1)),
           stat='identity', color='black', size=.2) +
  geom_text(data=annot, aes(copy_num, perc + .5, label=label),
            vjust=0, size=3.5) +
  scale_x_continuous('SMN1 + SMN2 copy number (exons 1-6)',
                     breaks=1:10, minor_breaks=NULL) +
  scale_y_continuous('Percentage of samples') +
  scale_fill_manual('SMN1 copy number',
                    values=RColorBrewer::brewer.pal(5, 'RdPu')[-1]) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(0, 0, 0, 0)),
        legend.position = 'top',
        legend.margin = margin(-2, 0, -7, 0),
        legend.key.height = unit(10, 'pt'))
ggsave('~/Tmp/1.png', width=8, height=6, scale=.65)

####################

annot2 <- data.frame()
for (cn in unique(paralog_counts$paralog1)) {
  label <- sprintf('Total: %s', (sum(filter(paralog_counts, paralog1 == cn)$n)))
  annot2 <- rbind(annot2, data.frame(paralog1=cn,
                                   perc=sum(filter(paralog_counts, paralog1 == cn)$perc),
                                   label=label))
}

(g_smn1 <- ggplot(paralog_counts) +
  geom_bar(aes(paralog1, perc, fill=factor(paralog2)),
           stat='identity', color='black', size=.2) +
  geom_text(data=annot2, aes(paralog1, perc + .5, label=label),
            vjust=0, size=3.5) +
  scale_x_continuous('SMN1 copy number',
                     breaks=1:10, minor_breaks=NULL) +
  scale_y_continuous('Percentage of samples') +
  scale_fill_manual('SMN2 copy number (exons 1-6)',
                    values=RColorBrewer::brewer.pal(5, 'RdPu')) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(0, 0, 0, 0)),
        legend.position = 'top',
        legend.margin = margin(-2, 0, -7, 0),
        legend.key.height = unit(10, 'pt')))
ggsave('~/Tmp/1.png', width=8, height=6, scale=.65)

ggplot(paralog_counts) +
  geom_tile(aes(x=factor(paralog1), y=factor(paralog2), fill=n), color='black') +
  geom_text(aes(x=factor(paralog1), y=factor(paralog2), label=n)) +
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, 'RdPu')[-8:-9],
                       trans='log', breaks=c(1, 10, 100, 1000)) +
  scale_x_discrete('SMN1 copy number', expand=c(0, 0)) +
  scale_y_discrete('SMN2 copy number', expand=c(0, 0)) +
  guides(fill=F) +
  theme_bw() +
  theme(panel.background=element_rect(fill='gray80'),
        panel.grid = element_blank())
ggsave('~/Tmp/1.png', width=8, height=6, scale=.65)

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
df <- load(gene, 'from_pops', keep_paralog=T, keep_qual=T)
df <- extend_paralog(df, 2)
# df <- filter(df, region != 'start')

df_wide <- to_wide(df)
filter(df_wide, copy_num_qual_start >= 30 & copy_num_qual_end >= 30) %>% nrow
filter(df_wide, copy_num_qual_start >= 30 & copy_num_qual_end >= 30) %>%
  filter(copy_num_start != copy_num_end) %>% nrow
filter(df_wide, copy_num_start != copy_num_end) %>% nrow

df_wide <- filter(df_wide, copy_num_qual_middle >= 30 & copy_num_qual_end >= 30 &
                    paralog_filter_end == 'PASS' & paralog_qual1_end >= 30) %>%
  mutate(short=copy_num_end - copy_num_middle,
         long=copy_num_middle)

region_names <- c('paralog1_end' = 'any C4A', #  'C4AS + C4AL',
                  'paralog2_end' = 'any C4B', #  'C4BS + C4BL',
                  'long' = 'any long', # 'C4AL + C4AL',
                  'short' = 'any short') # 'C4AS + C4BS')
cn_counts <- to_counts(df_wide, region_names,
        c('paralog1_end', 'paralog2_end', 'short', 'long'))
ref_cns <- ref_cns_df(cn_counts, c(2, 2, NA, NA))

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

df_wide$event <- with(df_wide,
                      sprintf('%s,%s,%s,%s', paralog1_end, paralog2_end, long, short))
rev(sort(table(df_wide$event)))

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

df <- load(gene, 'from_pops') %>% filter(region == 'gene_start') %>% extend_paralog(2)
df <- filter(df, !is.na(paralog1))

paralog_counts <- count(df, copy_num, paralog1, paralog2)
paralog_counts$perc <- paralog_counts$n / nrow(df) * 100

annot2 <- data.frame()
for (cn in unique(paralog_counts$paralog1)) {
  label <- sprintf('Total: %s', (sum(filter(paralog_counts, paralog1 == cn)$n)))
  annot2 <- rbind(annot2, data.frame(paralog1=cn,
                                     perc=sum(filter(paralog_counts, paralog1 == cn)$perc),
                                     label=label))
}

(g_strc <- ggplot(paralog_counts) +
  geom_bar(aes(paralog1, perc, fill=factor(paralog2)),
           stat='identity', color='black', size=.2) +
  geom_text(data=annot2, aes(paralog1, perc + .5, label=label),
            vjust=0, size=3.5) +
  scale_x_continuous('STRC copy number',
                     breaks=1:10, minor_breaks=NULL) +
  scale_y_continuous('Percentage of samples') +
  scale_fill_manual('pSTRC copy number',
                    values=RColorBrewer::brewer.pal(9, 'RdPu')[c(1, 4, 6, 9)]) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(0, 0, 0, 0)),
        legend.position = 'top',
        legend.margin = margin(-2, 0, -7, 0),
        legend.key.height = unit(10, 'pt')))
ggsave('~/Tmp/1.png', width=8, height=6, scale=.65)

ggplot(paralog_counts) +
  geom_tile(aes(x=factor(paralog1), y=factor(paralog2), fill=n), color='black') +
  geom_text(aes(x=factor(paralog1), y=factor(paralog2), label=n)) +
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, 'RdPu')[-8:-9],
                       trans='log', breaks=c(1, 10, 100, 1000)) +
  scale_x_discrete('STRC copy number', expand=c(0, 0)) +
  scale_y_discrete('pSTRC copy number', expand=c(0, 0)) +
  guides(fill=F) +
  theme_bw() +
  theme(panel.background=element_rect(fill='gray80'),
        panel.grid = element_blank())
ggsave('~/Tmp/1.png', width=8, height=6, scale=.65)

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

df2 <- mutate(df, superpopulation = 'AFR')
cn_counts <- to_counts(df2, c('copy_num' = 'CN'), 'copy_num') %>%
  mutate(superpopulation = NA)
ggplot(cn_counts) +
  geom_bar(aes(value, n), stat='identity') +
  theme_bw()

df <- load(gene, 'from_pops', keep_paralog=T)
# (paralog_table <- table(df$paralog_copy_num) / nrow(df) * 100)
# df <- mutate(df, paralog_group = case_when(
#   grepl('?', paralog_copy_num, fixed=T) ~ 'Other | Unknown',
#   copy_num < 4 | copy_num > 7 ~ paralog_copy_num,
#   paralog_table[paralog_copy_num] > 0.3 ~ paralog_copy_num,
#   T ~ 'Other | Unknown'
# ))
# df$paralog_group <- factor(df$paralog_group,
#     levels=rev(c('Other | Unknown', 'Unknown', 'Other', names(sort(paralog_table)))))

paralog_counts <- count(df, copy_num, paralog_copy_num)
paralog_counts$perc <- paralog_counts$n / nrow(df) * 100
annot <- data.frame()
for (cn in unique(paralog_counts$copy_num)) {
  label <- sprintf('Total: %s', (sum(filter(paralog_counts, copy_num == cn)$n)))
  for (i in which(paralog_counts$copy_num == cn)) {
    if (paralog_counts[i,]$n > 50) {
      label <- sprintf('%s\n%s: %s', label, paralog_counts[i,]$paralog_copy_num,
                       (paralog_counts[i,]$n))
    }
  }
  annot <- rbind(annot, data.frame(copy_num=cn,
                                   perc=sum(filter(paralog_counts, copy_num == cn)$perc),
                                   label=label))
}

aggr_counts <- count(df, copy_num)
aggr_counts$perc <- aggr_counts$n / nrow(df) * 100
breaks <- c(0, 1, 5, 10, 50, 100)
(g_neb <- ggplot(aggr_counts) +
  geom_bar(aes(copy_num, log(perc + 1), fill=factor(copy_num)),
           stat='identity', color='black', size=.2) +
  geom_text(data=annot, aes(copy_num, log(pmin(perc, 150) + 1) + .05, label=label),
            vjust=0, size=3.5) +
  scale_x_continuous('Aggregate copy number (NEB)', breaks=1:10, minor_breaks=NULL) +
  coord_cartesian(xlim=c(1.8, 10)) +
  scale_y_continuous('Percentage of samples (log)',
      breaks=log(1 + breaks), labels=breaks, minor_breaks=NULL,
      limits=c(0, 5), expand=expansion(add=.1)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(9, 'YlGnBu')[-1:-3]) +
  guides(fill=F) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(0, 0, 0, 0))))
ggsave('~/Tmp/1.png', width=10, height=6, scale=.65)

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

df <- load(gene, 'from_pops', keep_paralog=F) %>% filter(region == 'exon14')
aggr_counts <- count(df, copy_num)
aggr_counts$perc <- aggr_counts$n / nrow(df) * 100
annot <- data.frame()
for (cn in unique(aggr_counts$copy_num)) {
  label <- sprintf('Total: %s', (sum(filter(aggr_counts, copy_num == cn)$n)))
  annot <- rbind(annot, data.frame(copy_num=cn,
                                   perc=sum(filter(aggr_counts, copy_num == cn)$perc),
                                   label=label))
}
breaks <- c(0, 1, 5, 10, 50, 100)
(g_pms2 <- ggplot(aggr_counts) +
  geom_bar(aes(copy_num, log(perc + 1), fill=factor(copy_num)),
           stat='identity', color='black', size=.2) +
  geom_text(data=annot, aes(copy_num, log(pmin(perc, 150) + 1) + .05, label=label),
            vjust=0, size=3.5) +
  scale_x_continuous('Aggregate copy number (PMS2 exon 14)', breaks=1:10, minor_breaks=NULL) +
  #coord_cartesian(xlim=c(1.8, 10)) +
  scale_y_continuous('Percentage of samples (log)',
                     breaks=log(1 + breaks), labels=breaks, minor_breaks=NULL,
                     limits=c(0, 4.8), expand=expansion(add=.1)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(9, 'YlGnBu')[-1:-3]) +
  guides(fill=F) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(0, 0, 0, 0))))
ggsave('~/Tmp/1.png', width=8, height=6, scale=.65)

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

# ====== NCF1 ======

gene <- 'NCF1'
df <- load(gene, 'from_pops', keep_paralog=T, keep_qual=T)
df_wide <- df

region_names <- c('paralog_copy_num' = 'NCF1')
cn_counts <- to_counts(df_wide, region_names, 'paralog_copy_num')
ref_cns <- ref_cns_df(cn_counts, '2,2,2')

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

df_wide2 <- df_wide
df_wide2$paralog_copy_num <- (sum_perc[
  match(df_wide$paralog_copy_num, sum_perc$value),]$value2)
cn_counts2 <- to_counts(df_wide2, region_names, 'paralog_copy_num')

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


# =========================================

plot_grid(g_smn1, g_neb, g_strc, g_pms2, ncol=2, labels=LETTERS)
ggsave('~/Tmp/1.png', width=15, height=8, scale=.7)
