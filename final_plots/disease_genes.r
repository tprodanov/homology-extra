library(cowplot)

wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))
data_dir <- '~/Data/hg38/jvc/comparisons/genes/'
colors <- ggthemes::tableau_color_pal()(10)[c(1, 2, 3, 5)]
log_breaks <- c(0, 1, 5, 10, 50, 100)
QUAL <- 20

# ====== SMN1 ======

df <- load('SMN1', 'from_pops', keep_qual = T) %>% filter(grepl('SMN1_', region))
df <- extend_paralog(df, 2)

df$region <- sub('SMN1_', '', df$region)
df_wide <- to_wide(df)
df_wide <- filter_df(df_wide,
      agCN_qual_exon2 >= QUAL & agCN_qual_exon7 >= QUAL
      & psCN_filter_exon7 == 'PASS' & psCN_qual1_exon7 >= QUAL)

df_wide <- mutate(df_wide,
                  cn1 = agCN_exon2,
                  pcn1 = psCN1_exon7,
                  pcn2 = psCN2_exon7) %>%
  mutate(cn1=ifelse(cn1 >= 5, '≥ 5', as.character(cn1)),
         pcn1=as.character(pcn1), pcn2=as.character(pcn2))

region_names <- c('cn1' = 'SMN1 + 2\nexons 1-6',
                  'pcn1' = 'SMN1\nexons 7-8',
                  'pcn2' = 'SMN2\nexons 7-8')
cn_counts <- to_counts(df_wide, region_names, c('cn1', 'pcn1', 'pcn2'))
cn_counts <- mutate(cn_counts,
                    value=factor(value, levels=c(as.character(0:4), '≥ 5')))
# ref_cns <- ref_cns_df(cn_counts, c(4, 2, 2))

cn_counts <- filter(cn_counts, pop_name != 'Admixed-American')
(g_smn <- ggplot(cn_counts) +
    # geom_blank(aes(value)) +
    # geom_tile(aes(cn, 0), width=0.95, height=Inf, data=ref_cns, alpha=.2) +
    geom_bar(aes(value, perc, fill=pop_name), stat='identity',
             position=position_dodge2(preserve = 'single', padding = 0),
             width=0.8,
             color='black', size=0.3) +
    # facet_nested(label2 + label1 ~ .) +
    facet_wrap(~ region, ncol=1, strip.position='right') +
    scale_x_discrete('Copy number') +
    scale_y_continuous('Samples (%)', breaks=seq(0, 90, 30), minor_breaks=NULL) +
    scale_fill_manual('', values = colors[-1]) +
    theme_bw() +
    theme(# axis.line = element_line(color='black', size=.3),
          axis.title.y = element_text(hjust=0.5, margin=margin(0, 0, 0, 0)),
          axis.title.x = element_text(margin=margin(2, 0, 0, 0)),
          legend.title = element_blank(),
          legend.position = c(0.005, 1.05),
          legend.justification = c('left', 'top'),
          legend.key.height = unit(10, 'pt'),
          legend.spacing.x = unit(3, 'pt'),
          legend.background = element_rect(color=NA, fill='transparent'),
          strip.text = element_text(margin=margin(0, 2, 0, 2))))

sum(with(df_wide, superpopulation != 'AFR' & psCN1_exon7 >= 2)) /
  sum(with(df_wide, superpopulation != 'AFR'))

# ====== STRC ======

df <- load('STRC', 'from_pops', keep_qual = T) %>%
  filter(region == 'gene_start') %>% extend_paralog(2)
df <- filter(df, agCN_qual >= QUAL
             & !is.na(psCN1) & psCN_filter == 'PASS' & psCN_qual1 >= QUAL)

paralog_counts <- count(df, agCN, psCN1, psCN2)
paralog_counts$perc <- paralog_counts$n / nrow(df) * 100
# 
# annot <- data.frame()
# for (cn in unique(paralog_counts$agCN)) {
#   sum_n <- sum(filter(paralog_counts, agCN == cn)$n)
#   label <- sprintf('%s', sum_n)
#   annot <- rbind(annot, data.frame(agCN=cn, n=sum_n, label=label))
# }
# annot$perc <- 100 * annot$n / nrow(df)

# ref_cns <- ref_cns_df(cn_counts, c(4, 2, 2))
(g_strc <- ggplot(paralog_counts) +
    geom_bar(aes(agCN, log(1 + perc), fill=factor(psCN1)),
             stat='identity', color='black', size=.2, width=.95,
             position=position_dodge2(preserve = 'single', padding=0)) +
    scale_x_continuous('Aggregate copy number (STRC + pSTRC)',
                       breaks=1:10, minor_breaks=NULL,
                       expand=expansion(add=.1)) +
    scale_y_continuous('Samples (%, log-scale)', expand=expansion(add=.1),
                       breaks=log(1 + log_breaks), minor_breaks=NULL,
                       labels = log_breaks) +
    scale_fill_manual('STRC\npsCN',
                      values=RColorBrewer::brewer.pal(9, 'RdPu')[c(2, 5, 7, 9)]) +
    guides(fill = guide_legend(nrow=5, title.position = 'left',
                               label.position='left', hjust=1)) +
    theme_bw() +
    theme(axis.line = element_line(color='black', size=.3),
          panel.border = element_blank(),
          axis.title.x = element_text(margin=margin(2, 0, 0, 0)),
          axis.title.y = element_text(hjust=0.5, margin=margin(0, 2, 0, 0)),
          legend.position = c(0.995, 1.03),
          legend.title = element_text(hjust=0, vjust=0.9, margin=margin(0, 5, 0, 0)),
          legend.justification = c('right', 'top'),
          legend.background = element_rect(color=NA, fill='transparent'),
          legend.key.height = unit(10, 'pt'),
          legend.spacing.x = unit(3, 'pt')))

# ====== NEB ======

df <- load('NEB', 'from_pops', keep_paralog=T, keep_qual=T) %>% extend_paralog(3)
df <- filter(df, agCN_qual >= QUAL & psCN_filter == 'PASS'
         & psCN_qual1 >= QUAL & psCN_qual2 >= QUAL & psCN_qual3 >= QUAL)
df <- mutate(df, psCN = sprintf('%s,%s,%s', psCN1, psCN2, psCN3))

paralog_counts <- count(df, agCN, psCN)
paralog_counts$perc <- paralog_counts$n / nrow(df) * 100
# annot <- data.frame()
# for (cn in unique(paralog_counts$agCN)) {
#   label <- sprintf('Total: %s', (sum(filter(paralog_counts, agCN == cn)$n)))
#   for (i in which(paralog_counts$agCN == cn)) {
#     if (paralog_counts[i,]$n > 50) {
#       label <- sprintf('%s\n%s: %s', label, paralog_counts[i,]$psCN,
#                        (paralog_counts[i,]$n))
#     }
#   }
#   annot <- rbind(annot, data.frame(agCN=cn,
#                                    perc=sum(filter(paralog_counts, agCN == cn)$perc),
#                                    label=label))
# }

aggr_counts <- count(df, agCN)
aggr_counts$perc <- aggr_counts$n / nrow(df) * 100
breaks <- c(0, 1, 5, 10, 50, 100)
(g_neb <- ggplot(aggr_counts) +
    geom_bar(aes(agCN, log(perc + 1)),
             fill = RColorBrewer::brewer.pal(9, 'RdPu')[5],
             stat='identity', color='black', size=.2) +
    scale_x_continuous('Aggregate copy number (NEB, 3-copy intergenic repeat)',
                       breaks=1:10, minor_breaks=NULL) +
    coord_cartesian(xlim=c(1.8, 10)) +
    scale_y_continuous('Samples (%, log-scale)',
                       breaks=log(1 + breaks), labels=breaks, minor_breaks=NULL,
                       expand=expansion(add=.1)) +
    # scale_fill_manual(values=RColorBrewer::brewer.pal(9, 'YlGnBu')[-1:-3]) +
    # scale_fill_manual(values=RColorBrewer::brewer.pal(9, 'RdPu')[2:7]) +
    guides(fill=F) +
    theme_bw() +
    theme(axis.line = element_line(color='black', size=.3),
          panel.border = element_blank(),
          axis.title.y = element_text(hjust=0.5, margin=margin(0, 0, 0, 0)),
          axis.title.x = element_text(margin=margin(2, 0, 0, 0))))

# ====== PMS2 ======

df <- load('PMS2', 'from_pops', keep_paralog=F, keep_qual=T) %>% filter(region == 'exon14')
df <- filter(df, agCN_qual >= QUAL)
cn_counts <- to_counts(df, c('exon14' = 'Exon 14'), 'agCN')

cn_counts

(g_pms2 <- ggplot(filter(cn_counts, pop_name != 'Admixed-American')) +
  geom_bar(aes(value, log(1 + perc), fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single', padding = 0),
           width=1, color='black', size=.2) +
  scale_x_discrete('Aggregate copy number (PMS2 + PMS2CL)') +
  coord_cartesian(xlim=c(1.35, 3.55)) +
  scale_y_continuous('Samples (%, log-scale)',
                     breaks=log(1 + log_breaks), minor_breaks=NULL,
                     labels = log_breaks,
                     limits = c(0, log(101)),
                     expand=expansion(add=.1)) +
  scale_fill_manual('', values = colors) +
  guides(fill = guide_legend(nrow=5, label.position='left', hjust=1)) +
  theme_bw() +
  theme(axis.line = element_line(color='black', size=.3),
        panel.border = element_blank(),
        axis.title.y = element_text(hjust=0.5, margin=margin(0, 0, 0, 0)),
        axis.title.x = element_text(margin=margin(2, 0, 0, 0)),
        legend.title = element_blank(),
        legend.position = c(0.995, 1.05),
        legend.justification = c('right', 'top'),
        legend.key.height = unit(10, 'pt'),
        legend.spacing.x = unit(3, 'pt'),
        legend.background = element_rect(color=NA, fill='transparent')))

# ====== OTOA ======

df <- load('OTOA', 'from_pops', keep_paralog=F, keep_qual=T) %>% mutate(region='R')
df <- filter(df, agCN_qual >= QUAL)
df$agCN2 <- with(df, factor(case_when(
  # agCN == 4 ~ '4 (✳)',
  agCN >= 8 ~ '≥ 8',
  T ~ as.character(agCN),
), levels=c('4', '5', '6', '7', '≥ 8')))

cn_counts <- to_counts(df, c('R' = 'R'), 'agCN2')
ref_cns <- ref_cns_df(cn_counts, 4)
(g_otoa <- ggplot(filter(cn_counts, pop_name != 'Admixed-American'), aes(value)) +
    geom_blank() +
    geom_tile(aes(cn, 0), width=0.95, height=Inf, data=ref_cns, alpha=.2) +
    geom_bar(aes(value, perc, fill=pop_name), stat='identity',
             position=position_dodge2(preserve = 'single', padding = 0),
             width=0.8, color='black', size=.2) +
    scale_x_discrete('Aggregate copy number (OTOA + OTOAP1)',
                     expand=expansion(add=0.2)) +
    scale_y_continuous('Samples (%)',
      breaks=seq(0, 100, 20), expand=expansion(add=1)) +
    scale_fill_manual('', values = colors) +
    guides(fill = guide_legend(nrow=5, label.position='left', hjust=1)) +
    theme_bw() +
    theme(axis.line = element_line(color='black', size=.3),
          panel.border = element_blank(),
          axis.title.y = element_text(hjust=0.5, margin=margin(0, 2, 0, 0)),
          axis.title.x = element_text(margin=margin(2, 0, 0, 0)),
          legend.title = element_blank(),
          legend.position = c(0.995, 1.05),
          legend.justification = c('right', 'top'),
          legend.key.height = unit(10, 'pt'),
          legend.spacing.x = unit(3, 'pt'),
          legend.background = element_rect(color=NA, fill='transparent')))
ggsave('~/Tmp/1.png', width=8, height=4, scale=.7, dpi=500)
ggsave('~/Tmp/1.png', width=8, height=3.5, scale=.5, dpi=500)

# ====== OTOA (3 copies) ======

df <- load('OTOA', '3c', keep_paralog=T, keep_qual=T) %>% mutate(region='R') %>%
  extend_paralog(3)
df <- filter(df, agCN_qual >= QUAL & psCN_filter == 'PASS'
             & psCN_qual1 >= QUAL & psCN_qual2 >= QUAL & psCN_qual3 >= QUAL)

paralog_counts <- count(df, agCN, psCN1)
paralog_counts$perc <- paralog_counts$n / nrow(df) * 100

(g_otoa <- ggplot(paralog_counts) +
    annotate('rect', xmin = 4 - 0.3, xmax = 4 + 0.3, ymin = -Inf, ymax = Inf,
             alpha = .2) +
    geom_bar(aes(agCN, log(1 + perc), fill=factor(psCN1)),
             stat='identity', color='black', size=.2, width=.85,
             position=position_dodge2(preserve = 'single', padding=0)) +
    scale_x_continuous('Aggregate copy number (OTOA + OTOAP1)',
                       breaks=1:10, minor_breaks=NULL,
                       expand=expansion(add=.1)) +
    scale_y_continuous('Samples (%, log-scale)', expand=expansion(add=.1),
                       breaks=log(1 + log_breaks), minor_breaks=NULL,
                       labels = log_breaks) +
    scale_fill_manual('OTOA\npsCN',
                      values=RColorBrewer::brewer.pal(9, 'RdPu')[c(2, 5, 8)]) +
    guides(fill = guide_legend(nrow=5, title.position = 'left',
                               label.position='left', hjust=1)) +
    theme_bw() +
    theme(axis.line = element_line(color='black', size=.3),
          panel.border = element_blank(),
          axis.title.x = element_text(margin=margin(2, 0, 0, 0)),
          axis.title.y = element_text(hjust=0.5, margin=margin(0, 2, 0, 0)),
          legend.position = c(0.995, 1.03),
          legend.title = element_text(hjust=0, vjust=0.9, margin=margin(0, 5, 0, 0)),
          legend.justification = c('right', 'top'),
          legend.background = element_rect(color=NA, fill='transparent'),
          legend.key.height = unit(10, 'pt'),
          legend.spacing.x = unit(3, 'pt')))

# ====== NCF1 =======

df <- load('NCF1', 'from_pops', keep_qual = T) %>% extend_paralog(3)
df <- filter(df, agCN_qual >= QUAL & psCN_filter == 'PASS' & psCN_qual1 >= QUAL)
paralog_counts <- count(df, agCN, psCN1) %>% mutate(perc = 100 * n / nrow(df))

(g_ncf1 <- ggplot(paralog_counts) +
    geom_bar(aes(agCN, log(1 + perc), fill=factor(psCN1)),
             stat='identity', color='black', size=.2, width=1,
             position=position_dodge2(preserve = 'single', padding=0)) +
    scale_x_continuous('Aggregate copy number (NCF1 + NCF1B + NCF1C)',
                       breaks=1:10, minor_breaks=NULL,
                       expand=expansion(add=.05)) +
    scale_y_continuous('Samples (%, log-scale)', expand=expansion(add=.1),
                       breaks=log(1 + log_breaks), minor_breaks=NULL,
                       labels = log_breaks) +
    scale_fill_manual('NCF1\npsCN',
                      values=RColorBrewer::brewer.pal(9, 'RdPu')[c(2, 5, 7, 9)]) +
    guides(fill = guide_legend(nrow=5, title.position = 'left',
                               label.position='left', hjust=1)) +
    theme_bw() +
    theme(axis.line = element_line(color='black', size=.3),
          panel.border = element_blank(),
          axis.title.x = element_text(margin=margin(2, 0, 0, 0)),
          axis.title.y = element_text(hjust=0.5, margin=margin(0, 2, 0, 0)),
          legend.position = c(0.995, 1.03),
          legend.title = element_text(hjust=0, vjust=0.9, margin=margin(0, 5, 0, 0)),
          legend.justification = c('right', 'top'),
          legend.background = element_rect(color=NA, fill='transparent'),
          legend.key.height = unit(10, 'pt')))

sum(with(df, psCN1 >= 2)) /
  nrow(df)

# ====== Plot grid ======

# plot_grid(g_strc, g_neb, g_pms2, g_otoa, ncol=2, labels=LETTERS,
#           label_x=-0.01, label_y=1.02)
# ggsave('~/Tmp/1.png', width=15, height=8, scale=.6, dpi=600)

plot_grid(g_smn, g_pms2, g_ncf1, g_strc, g_neb, g_otoa, ncol=2, labels=LETTERS,
          label_x=-0.01, label_y=1.02)
ggsave('~/Tmp/1.png', width=15, height=12, scale=.6, dpi=600)

