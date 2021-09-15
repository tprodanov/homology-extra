library(ggplot2)

comparison <- read.csv('~/Data/hg38/jvc/tools/QuickMer2_1KG/comparison/000.bed',
         sep='\t', comment.char='#', header=F)
names(comparison) <- c('chrom', 'start', 'end', 'region', 'sample',
  'other_regions', 'ploidy_filter', 'ploidy', 'ploidy_qual', 'paralog_filter',
  'paralog_ploidy', 'paralog_qual', 'info', 'qm2_ploidy', 'qm2_paralog',
  'qm2_entries', 'pooled_err', 'paralog_err', 'paralog_int_err')

pos1 <- 70940000
comp1 <- subset(comparison, start <= pos1 & pos1 < end
                & ploidy_filter == 'PASS' & paralog_filter == 'PASS')

pos2 <- 70950000
comp2 <- subset(comparison, start <= pos2 & pos2 < end
                & ploidy_filter == 'PASS' & paralog_filter == 'PASS')

ggplot(comp1) +
  geom_histogram(aes(pooled_err), binwidth=.05) +
  theme_bw()

ggplot(comp1) +
  geom_bar(aes(ploidy - round(qm2_ploidy)), width=0.1) +
  theme_bw()

ggplot(comp1) +
  geom_histogram(aes(paralog_err), binwidth=0.2) +
  theme_bw()

ggplot(comp1) +
  geom_bar(aes(paralog_int_err), width=0.2) +
  theme_bw()

View(comp1[comp1$paralog_int_err > 1,])
