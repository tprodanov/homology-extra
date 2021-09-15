library(ggplot2)
library(tidyverse)

args <- parser$parse_args(c('-i', '~/Data/hg38/jvc/runs/242.EUR/SMN1/v100/extra/',
                            '-o', '~/Tmp/'))

smn_annot <- read_delim('~/Data/hg38/genome/annot/annotation.SMN1.gff3', '\t', col_names=F)
names(smn_annot) <- c('chrom', 'db', 'type', 'start', 'end', 'score', 'strand',
                      'phase', 'attr')

start <- filter(smn_annot, type == 'gene')$start
end <- filter(smn_annot, type == 'gene')$end
hmm_filter <- filter(hmm_states_last, middle >= start & middle <= end)

exons <- filter(smn_annot, type == 'exon' & grepl('SMN1-202', attr)) %>%
  select(db, start, end)
exons$ix <- 1:nrow(exons)
exons_wide <- pivot_wider(exons, ix, values_from=c('start', 'end'))
names(exons_wide) <- c('ix', 'start', 'end')
exons_wide$name <- sprintf('Exon %s', c(1, '2a', '2b', 3:8))
exons_wide$nudge_x <- 0
exons_wide$nudge_x[4:5] <- c(-100, 100)
# exons_wide$name <- sprintf('Exon %d', exons_wide$ix)

hmm_filter$cn <- hmm_filter$pred_cn
left <- hmm_filter %>% select(sample, ix, middle, cn)
right <- hmm_filter %>% select(sample, ix, middle, cn) %>% mutate(ix = ix - 1)
hmm_embed <- inner_join(left, right, by=c('sample', 'ix'), suffix=c('.1', '.2'))
hmm_embed <- mutate(hmm_embed, diff = cn.2 - cn.1) %>%
  mutate(diff_f = factor(diff, levels=c('-1', '-2')))
hmm_segments <- count(hmm_embed,
                      ix, middle.1, cn.1, middle.2, cn.2, diff, diff_f)
hmm_segments <- mutate(hmm_segments,
    has_arrow = diff < 0 | ix == 280 | ix == 415)

hmm_segments0 <- filter(hmm_segments, diff == 0)
hmm_segments1 <- filter(hmm_segments, diff == -1)
hmm_segments2 <- filter(hmm_segments, diff == -2)
arrow <- arrow(length=unit(5, 'pt'))

del_pos <- 70948287

(g_viterbi <- ggplot(hmm_segments0) +
    geom_blank(aes(color=diff_f)) +
    geom_vline(xintercept = del_pos, linetype='dashed') +
    annotate('text', x = del_pos - 200, y = 0.52, label='Deletion',
              angle=90, size=2.8, hjust=0, vjust=0) +
    # Exons.
    geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), data=exons_wide,
              fill='black', alpha=.2) +
    geom_text(aes(x = (start + end) / 2, y = 0.52, label=name), data=exons_wide,
              angle=90, size=2.8, hjust=0,
              nudge_x=exons_wide$nudge_x) +
    # Segments with diff = 0.
    geom_segment(aes(x = middle.1, xend = middle.2, y = cn.1, yend = cn.2)) +
    geom_segment(aes(x = middle.1, xend = middle.2, y = cn.1, yend = cn.2),
                 data=filter(hmm_segments0, has_arrow), arrow = arrow) +
    geom_text(aes(x = middle.1, y = (cn.1 + cn.2) / 2, label = n),
              data = filter(hmm_segments0, has_arrow), nudge_y = 0.4) +
    # Segments with diff = -1.
    geom_curve(aes(x = middle.1, xend = middle.2, y = cn.1, yend = cn.2, color = diff_f),
               data = hmm_segments1, arrow = arrow, curvature = -0.2, alpha=.6,
               size = 0.8) +
    geom_text(aes(x = middle.1, y = (cn.1 + cn.2) / 2, label = n, color = diff_f),
               data = hmm_segments1, show.legend=F, nudge_x = 1200) +
    # Segments with diff = -2.
    geom_curve(aes(x = middle.1, xend = middle.2, y = cn.1, yend = cn.2, color = diff_f),
               data = hmm_segments2, arrow = arrow, curvature = 0.2, alpha=.6,
               size = 0.8) +
    geom_text(aes(x = middle.1, y = (cn.1 + cn.2) / 2, label = n, color = diff_f),
              data = hmm_segments2, show.legend=F, nudge_x = -700, nudge_y = 0.3) +
    # Scales.
    scale_color_manual('agCN transition', values = c('blue', 'red')) +
    scale_x_continuous('Position (chr5, kb)',
                       labels=function(x) format(x / 1000, big.mark=','),
                       expand=expansion(mult = 0.02)) +
    scale_y_continuous('Aggregate copy number', breaks=2:6, expand=c(0, 0)) +
    coord_cartesian(ylim=c(0.5, 6.7)) +
    guides(label = F) +
    theme_bw() +
    theme(strip.text = element_text(margin = margin(2, 0, 2, 0)),
          legend.position = 'top',
          legend.margin = margin(-5, 0, -10, 0),
          legend.background = element_rect(fill='transparent'),
          panel.grid = element_blank()))
ggsave('~/Tmp/1.png', width=5, height=3, scale=.85, dpi=600)


length(unique(hmm_filter$sample))
