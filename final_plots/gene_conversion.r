library(ggplot2)
library(tidyverse)

setwd('~/Data/hg38/jvc/runs')
plots_dir <- '~/Data/hg38/jvc/plots/gene_conversion/'

# RHCE Exons
gene <- 'RHCE'
transcript <- 'RHCE-204'
nudge_x <- rep(-750, 12)
nudge_x[c(2, 12)] <- -1000
exon_names <- paste('Exon', 1:12)

# SMN1 Exons
gene <- 'SMN1'
transcript <- 'SMN1-202'
nudge_x <- rep(-750, 9)
# nudge_x[c(2, 12)] <- -1000
exon_names <- paste('Exon', c('1', '2a', '2b', '3', '4', '5', '6', '7', '8'))

# Load annotation

annotation <- read_delim(sprintf('~/Data/hg38/genome/annot/annotation.%s.gff3', gene),
                         '\t', comment = '#', col_names = F)
gene_start <- min(annotation$X4)
gene_end <- max(annotation$X5)
exons <- filter(annotation, X3 == 'exon' & grepl(transcript, X9, fixed=T))
exons$name <- exon_names

# Select sample

dir <- '241.EAS/RHCE/r009c/extra'
dir <- '241.EAS/SMN1/r009/extra'

fval <- read_delim(file.path(dir, 'em_f_values.csv'), '\t', comment = '#')
fval <- mutate(fval,
               pos = as.numeric(sapply(strsplit(fval$psv, ':'), `[`, 2)),
               copy1 = as.numeric(copy1),
               copy2 = as.numeric(copy2),
               copy3 = as.numeric(copy3),
               copy4 = as.numeric(copy4)) %>%
  mutate(min_fval = pmin(copy1, copy2, copy3, copy4, na.rm = T))

gene_conv <- read_delim(file.path(dir, 'gene_conversion.bed'), '\t')
gene_conv <- gene_conv[order(-gene_conv$qual),]
filter(gene_conv, start <= gene_end & end >= gene_start & qual >= 20) %>% nrow
filter(gene_conv, qual >= 20) %>% nrow

sample <- 'NA18555'
sample_conv <- gene_conv[gene_conv$sample == sample,]

psv_obs_all <- read_delim(file.path(dir, 'psv_observations.csv'), '\t', comment='#')
psv_obs <- psv_obs_all[psv_obs_all$sample == sample,]
psv_obs <- mutate(psv_obs, total = ref_cov + alt_cov) %>%
  mutate(frac = ref_cov / total)
psv_obs <- left_join(psv_obs, fval, by='pos') %>%
  mutate(ftype = factor(case_when(
    is.na(min_fval) ~ 'Not used',
    min_fval >= 0.95 ~ 'Reliable',
    min_fval >= 0.8 ~ 'Semi-reliable',
    T ~ 'Unreliable'), levels=c('Reliable', 'Semi-reliable', 'Unreliable', 'Not used')))

colors <- ggthemes::tableau_color_pal()(10)

ggplot(psv_obs) +
  #geom_rect(data=sample_conv,
            #aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf),
            #fill='yellow', alpha=.1) +

  geom_rect(aes(xmin=X4, xmax=X5, ymin=-Inf, ymax=Inf), data=exons, alpha=0.4) +
  geom_text(aes(X4 + nudge_x, -0.25, label=name), data=exons, hjust=0, angle=90, size=3) +
  geom_vline(xintercept = c(gene_start, gene_end)) +
  # annotate('text', gene_end - 500, -0.2, label='Gene start', hjust=0, angle=90) +
  
  geom_point(aes(pos, frac, color = ftype), alpha=.7) +
  # coord_cartesian(xlim = c(gene_start - 1000, gene_end + 2000)) +
  
  scale_x_continuous(sprintf('Position (%s, kb)', annotation$X1[1]),
                     labels=function(x) format(x / 1000, big.mark = ','),
                     expand=expansion(mult=0.02),
                     breaks=(10000 * round(gene_start / 10000):round(gene_end / 10000)),
                     minor_breaks=NULL) +
  scale_y_continuous('Fraction of reads with Copy 1 allele', breaks=seq(0, 1, .25)) +
  scale_color_manual('PSV type   ', values=c(colors[c(5, 1, 3)], 'gray40')) +
  ggtitle(sprintf('%s: %s', gene, sample)) +
  guides(color = guide_legend(override.aes = list(size=2, alpha=1))) +
  theme_light() +
  theme(panel.grid.major.x = element_line(color=NA),
        legend.position = 'bottom',
        legend.margin = margin(-10, 0, -2, 0),
        legend.key.width = unit(1, 'pt'),
        legend.spacing = unit(1, 'pt'),
        legend.text = element_text(margin=margin(0, 5, 0, 0)))
# ggsave('~/Tmp/1.png', width=6, height=4, scale=1, dpi=450)
ggsave(sprintf('%s/%s.%s.png', plots_dir, gene, sample), width=6, height=4, scale=1, dpi=450)
