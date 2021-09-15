library(tidyverse)
library(ggplot2)
library(cowplot)

QUAL <- 20

all_f_val <- read.csv('~/Data/hg38/jvc/comparisons/f_val/f_val.genes.bed',
           sep='\t', comment.char='#', header=F,
           col.names=c('chrom', 'start', 'end', 'filename',
                       'group', 'n_samples', 'in_em',
                       'info_content', 'copy1', 'copy2', 'copy3', 'copy4',
                       'gene_chrom', 'gene_start', 'gene_end', 'gene', 'gene_type'))
all_f_val$min_fval <- apply(select(all_f_val, starts_with('copy')), 1, min, na.rm=T)
all_f_val$min_fval <- with(all_f_val, ifelse(is.finite(min_fval), min_fval, 0))
all_f_val$psv_type <- with(all_f_val, factor(case_when(
  min_fval >= 0.95 ~ 'Reliable',
  min_fval >= 0.8 ~ 'Semi-reliable',
  T ~ 'Unreliable'), levels=c('Unreliable', 'Semi-reliable', 'Reliable')))

all_f_val <- mutate(all_f_val,
    dataset=sapply(strsplit(filename, '/', fixed=T), `[`, 2),
    region=sapply(strsplit(filename, '/', fixed=T), `[`, 3),
    run=sapply(strsplit(filename, '/', fixed=T), `[`, 4),
    pos=sprintf('%s:%s', chrom, start)) %>%
  mutate(dataset=sapply(strsplit(dataset, '.', fixed=T), `[`, 2))

# Sort the dataset.
# all_f_val <- all_f_val[with(all_f_val, order(dataset, region, run)),]

# Remove cases when the same gene was analyzed several times.
# f_val <- group_by(all_f_val, dataset, pos) %>% slice_tail(n=1) %>% ungroup()
# filter(f_val, gene == 'C4A') %>% select(run) %>% unique

f_val <- all_f_val
pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')
f_val$pop <- pop_names[f_val$dataset]
f_val <- filter(f_val, !is.na(pop) & !is.na(gene))
f_val2 <- filter(f_val, gene_type == 'protein_coding')

type_counts <- count(f_val2, pop, gene, gene_type, psv_type)
total_counts <- count(f_val2, pop, gene, name='total')
type_counts <- left_join(type_counts, total_counts)
type_counts$perc <- with(type_counts, n / total * 100)
# At least 10 PSVs.
# type_counts <- filter(type_counts, total >= 10)

perc_sums <- aggregate(perc ~ gene + psv_type, type_counts, sum) %>%
  pivot_wider(names_from='psv_type', values_from='perc')
names(perc_sums)[3] <- 'Semi'
perc_sums <- mutate(perc_sums,
    Unreliable = replace_na(Unreliable, 0),
    Semi = replace_na(Semi, 0), Reliable = replace_na(Reliable, 0)) %>%
  mutate(sum_obs = Unreliable + Semi + Reliable) %>%
  mutate(value = (0.001 * Unreliable + 0.01 * Semi + Reliable) / sum_obs)
gene_ord <- perc_sums$gene[order(perc_sums$value)]
type_counts$gene2 <- factor(type_counts$gene, levels=gene_ord)

gene_cov <- read_delim('~/Data/hg38/jvc/abstract/primary_genes.bed',
                       '\t', comment = '#',
  col_names = c('chrom', 'start', 'end',
                'gene_chrom', 'gene_start', 'gene_end', 'gene', 'gene_type'))
gene_cov <- mutate(gene_cov,
    cov = pmin(end, gene_end) - pmax(start, gene_start)) %>%
  mutate(perc = cov / (gene_end - gene_start) * 100)
gene_cov$n_psvs <- type_counts[match(gene_cov$gene, type_counts$gene),]$total
gene_cov$n_psvs[is.na(gene_cov$n_psvs)] <- 0
gene_cov$psv_density <- gene_cov$n_psvs / gene_cov$cov * 1000
gene_cov$gene2 <- factor(gene_cov$gene, levels=gene_ord)

disease_assoc <- read_delim('~/Data/hg38/genome/curated_gene_disease_associations.tsv',
                            '\t')
# disease_assoc <- filter(disease_assoc, diseaseType == 'disease')

gene_subset <- read_delim('~/Data/hg38/jvc/regions/gene_centers2.bed',
                          '\t', col_names = F)$X4

# gene_subset <- unique(type_counts$gene)
# gene_subset <- gene_subset[gene_subset %in% disease_assoc$geneSymbol]

# gene_subset <- gene_subset[gene_subset %in% filter(perc_sums, sum_obs == 500)$gene]
# Reorder gene subset.
gene_subset <- intersect(gene_ord, gene_subset)
# gene_subset <- setdiff(gene_subset, 'KRT86')

type_counts$ix <- match(type_counts$gene, gene_subset)
(g1 <- ggplot(filter(type_counts, gene %in% gene_subset & pop != 'Admixed-American')) +
  geom_bar(aes(ix, perc, fill=psv_type),
           width=1, stat='identity', color='black', size=.1) +
  geom_hline(yintercept=c(25, 50, 75), linetype='solid', alpha=.3, size=.3) +
  scale_x_continuous(expand=c(0, 0), position='top',
                     breaks=1:length(gene_subset) + 0.2, labels=gene_subset) +
  scale_y_continuous('Percentage of PSVs', expand=c(0, 0), breaks=c(0, 50, 100)) +
  scale_fill_manual('PSV type  ',
                    values=RColorBrewer::brewer.pal(11, 'RdYlGn')[c(3, 6, 9)]) +
  facet_wrap(~ pop, ncol=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=7, vjust=1, hjust=0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=8),
        panel.background = element_rect(fill='gray40'),
        plot.background = element_rect(size=0, fill = 'transparent'),
        panel.grid = element_blank(),
        strip.text = element_text(margin=margin(2, 0, 2, 0)),
        plot.margin = margin(5, 2, 0, 2),
        legend.position = 'top',
        legend.key.height = unit(10, 'pt'),
        legend.key.width = unit(13, 'pt'),
        legend.margin = margin(0, 0, -5, 0)))

gene_cov <- mutate(gene_cov, n_psvs_group = case_when(
  n_psvs == 0 ~ '0 ',
  n_psvs <= 5 ~ '1-5 ',
  n_psvs <= 10 ~ '6-10 ',
  n_psvs <= 20 ~ '11-20 ',
  n_psvs <= 50 ~ '21-50 ',
  n_psvs <= 100 ~ '51-100 ',
  T ~ '> 100 '
))
gene_cov$n_psvs_group <- factor(gene_cov$n_psvs_group,
    levels=unique(gene_cov[order(gene_cov$n_psvs),]$n_psvs_group))
# n_colors <- length(unique(filter(gene_cov, gene %in% gene_subset)$n_psvs_group))
colors <- RColorBrewer::brewer.pal(11, 'RdYlBu')[c(2, 3, 4, 9, 10, 11)]

gene_cov$ix <- match(gene_cov$gene, gene_subset)
(g2 <- ggplot(filter(gene_cov, gene %in% gene_subset)) +
  geom_bar(aes(ix, psv_density, fill=n_psvs_group), stat='identity') +
  scale_x_continuous(expand=c(0, 0),
                     breaks=1:length(gene_subset), labels=gene_subset) +
  scale_y_continuous('PSVs / 1kb', expand=expansion(mult=.02),
                     breaks=seq(0, 100, 10)) +
  scale_fill_manual('Total number of PSVs ', values=colors) +
  guides(fill=guide_legend(nrow=1)) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle=90, size=7, vjust=0.5, hjust=1,
                                   #margin=margin(1.4, 0, 0, 0)),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(margin=margin(0, -8, 0, 0)),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        #panel.grid.minor.x = element_line(color='gray50'),
        strip.text = element_text(margin=margin(2, 0, 2, 0)),
        plot.margin = margin(3, 2, 2, 2),
        plot.background = element_rect(size=0, fill = 'transparent'),
        legend.position = 'top',
        legend.key.height = unit(10, 'pt'),
        legend.key.width = unit(13, 'pt'),
        legend.spacing.x = unit(3, 'pt'),
        legend.margin = margin(0, 0, -10, 0)))

# do.call(paste, as.list(gene_subset))

# ====== Load CN distribution ======

exons <- read_delim('~/Data/hg38/jvc/regions/gene_exons/genes_exons.isec.bed',
           '\t', col_names = c('chrom', 'start', 'end', 'gene'))

summaries <- data.frame()
for (filename in sprintf('~/Data/hg38/jvc/comparisons/children/trios/%s.bed',
                         c('241.EAS', '242.EUR', '243.AFR', '244.AMR', '245.SAS'))) {
  tmp <- read.csv(filename, sep = '\t', comment.char = '#', header = F)
  summaries <- rbind(summaries, tmp)
}
names(summaries) <- c('chrom', 'start', 'end', 'locus', 'sample', 'reg2',
                      'filter', 'copy_num', 'qual',
                      'psCN_filter', 'psCN', 'psCN_qual', 'gene', 'exon')
summaries$ref_cn <- 2 * (stringr::str_count(summaries$reg2, 'chr') + 1)
table(summaries$ref_cn)
# summaries2 <- filter(summaries, ref_cn > 2)

# positions <- data.frame()
# copy_num_distr <- data.frame()
# for (gene in unique(exons$gene)) {
#   curr_exons <- exons[exons$gene == gene,]
#   exon_chrom = curr_exons[1,]$chrom
#   best <- 0
#   cn_obs <- NULL
#   best_pos <- NULL
#   best_exon <- NULL
#   for (exon_ix in 1:nrow(curr_exons)) {
#     pos <- with(curr_exons[exon_ix,], (start + end) %/% 2)
#     tmp <- filter(summaries2, chrom == exon_chrom & start <= pos & pos < end) %>%
#       select(sample, copy_num, qual, ref_cn) %>% add_column(gene = gene)
#     curr <- sum(tmp$qual >= QUAL)
#     # cat(sprintf('%s exon %d: %d\n', gene, exon_ix, curr))
#     if (curr > best) {
#       best <- curr
#       cn_obs <- tmp
#       best_pos <- pos
#       best_exon <- exon_ix
#     }
#   }
#   # cat(sprintf('%s best: %d\n', gene, best))
#   copy_num_distr <- rbind(copy_num_distr, cn_obs)
#   positions <- rbind(positions,
#                      data.frame(gene=gene, chrom=exon_chrom, pos=best_pos, exon=best_exon))
# }
# positions$start <- positions$pos - 1
# positions$end <- positions$pos
# positions <- positions[c('chrom', 'start', 'end', 'gene', 'exon')]
# write_delim(positions, '~/Data/hg38/jvc/regions/gene_centers.unsort.bed', '\t',
#             col_names = F)

# positions <- read_delim('~/Data/hg38/jvc/regions/gene_centers2.bed', '\t', comment = '#',
#                         col_names = F)
# names(positions) <- c('chrom', 'start', 'end', 'gene', 'exon')
# table(positions$gene)[table(positions$gene) > 1]
# positions <- filter(positions, case_when(
#   gene == 'NCF1' ~ exon == 7,
#   gene == 'SMN1' ~ exon == 2,
#   T ~ T
# ))
# 
# copy_num_distr <- data.frame()
# for (i in 1:nrow(positions)) {
#   pos_chrom <- positions[i,]$chrom
#   pos_start <- positions[i,]$start
#   pos_end <- positions[i,]$end
# 
#   subs <- filter(summaries2, chrom == pos_chrom & start <= pos_end & end >= pos_start)
#   subs$gene <- positions[i,]$gene
#   subs$exon <- positions[i,]$exon
#   copy_num_distr <- rbind(copy_num_distr, subs)
# }

sum(summaries$qual >= QUAL)
copy_num_distr <- summaries
copy_num_distr$copy_num_f <- factor(as.character(copy_num_distr$copy_num),
                                       levels=as.character(0:20))
fcopy_num_distr <- filter(copy_num_distr, qual >= QUAL)
fcopy_num_distr$ix <- match(fcopy_num_distr$gene, gene_subset)
copy_num_counts <- count(fcopy_num_distr, ix, gene, copy_num)
ref_cn <- count(fcopy_num_distr, ix, gene, ref_cn)
copy_num_counts$perc <- 100 * copy_num_counts$n /
  ref_cn[match(copy_num_counts$gene, ref_cn$gene),]$n

copy_num_counts <- mutate(copy_num_counts, perc_group = case_when(
  perc <= 1 ~ '*',
  perc < 10 ~ '1-10 ',
  perc < 25 ~ '10-25 ',
  perc < 50 ~ '25-50 ',
  perc < 75 ~ '50-75 ',
  T ~ '75-100'
))

(g_cn <- ggplot(filter(copy_num_counts, perc_group != '*')) +
  geom_tile(aes(ix, copy_num, fill=perc_group), color='black', size=.1) +
  geom_point(aes(ix, ref_cn), data=ref_cn, color='white', size=1) +
  scale_fill_manual('Percentage of samples ',
    values = RColorBrewer::brewer.pal(9, 'YlGn')[c(2, 3, 5, 7, 9)]) +
  #scale_fill_gradientn('Percentage of samples',
    #colors = RColorBrewer::brewer.pal(9, 'YlGn')[c(-1, -9)],
                       #limits=c(0, 100)) +
  scale_x_continuous(expand=c(0, 0),
                     breaks=1:length(gene_subset), labels=gene_subset) +
  scale_y_continuous('agCN', expand=expansion(add=0.2),
                     breaks=seq(0, 20, 2)) +
  #guides(fill=guide_legend(nrow=1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=7, vjust=0.5, hjust=1,
                                   margin=margin(1.4, 0, 0, 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(margin=margin(0, -8, 0, 0)),
        panel.grid.major.x = element_blank(),
        #panel.grid.minor.x = element_line(color='gray50'),
        strip.text = element_text(margin=margin(2, 0, 2, 0)),
        plot.margin = margin(0, 2, 2, 2),
        plot.background = element_rect(size=0, fill = 'transparent'),
        legend.position = c(0.999, 0.9),
        legend.justification = c('right', 'top'),
        legend.direction = 'horizontal',
        legend.key.height = unit(10, 'pt'),
        legend.key.width = unit(13, 'pt'),
        legend.spacing.x = unit(3, 'pt'),
        legend.margin = margin(-5, 0, 0, 0)))

plot_grid(g1, g2, g_cn, ncol=1, align='v', rel_heights = c(0.9, 0.3, 0.5),
          labels=LETTERS, label_x = -0.005, label_y = c(1, 1.07, 1))
ggsave('~/Tmp/1.png', width=12, height=8, scale=0.75, dpi=600)
