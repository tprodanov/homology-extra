library(ggplot2)
library(tidyverse)

filename <- '~/Data/hg38/jvc/comparisons/han/r009/concordance.csv'
out <- '~/Data/hg38/jvc/plots/han/r009/concordance.png'
filename <- '~/Data/hg38/jvc/comparisons/children/r009/concordance.csv'
out <- '~/Data/hg38/jvc/plots/246.children/r009c/concordance.png'

conc <- read_delim(filename, '\t', comment = '#')
conc_aggr <- aggregate(cbind(concordant, total) ~ dataset1 + dataset2 + gene, conc, sum)
conc_aggr <- mutate(conc_aggr, perc = concordant / total * 100)

length(unique(conc$sample))
nrow(conc_aggr)
sum(conc_aggr$perc >= 90)
sum(conc_aggr$perc >= 95)
mean(conc_aggr$perc)
sum(conc_aggr$total) / 1e8

conc_aggr$gene_ord <- reorder(conc_aggr$gene, conc_aggr$perc)
ggplot(conc_aggr) +
  geom_bar(aes(gene_ord, perc), stat='identity', fill='#9999ff', color='gray30', size=0.2) +
  geom_hline(yintercept = 90) +
  scale_x_discrete('Gene') +
  scale_y_continuous('Concordance (%)', breaks=seq(0, 100, 10)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=5))
ggsave(out, width=11, height=4)

thresholds <- c(0, 50, 80, 90, 95, 99, 100)
for (i in rev(seq_along(thresholds))) {
  lower <- thresholds[i]
  upper <- ifelse(lower == 100, 101, thresholds[i + 1])
  cat(sprintf('in [%d, %d): %d genes\n', lower, upper,
              sum(with(conc_aggr, perc >= lower & perc < upper))))
  cat(sprintf('    %s\n=====\n',
              do.call(paste, as.list(filter(conc_aggr, perc >= lower & perc < upper)$gene))))
}

comp <- read.csv('~/Data/hg38/jvc/comparisons/children/r009/comparison.csv',
         sep='\t', comment.char='#')
QUAL <- 30
fcomp <- filter(comp, ploidy1 == ploidy2 &
                  grepl(',', paralog_ploidy1),
         ploidy_filter1 == 'PASS' & ploidy_filter2 == 'PASS' &
         ploidy_qual1 >= QUAL & ploidy_qual2 >= QUAL &
         paralog_filter1 == 'PASS' & paralog_filter2 == 'PASS' &
         !grepl('\\?', paralog_ploidy1) & !grepl('\\?', paralog_ploidy2))
View(filter(fcomp, paralog_ploidy1 != paralog_ploidy2))

sum(fcomp$length)
sum(filter(fcomp, paralog_ploidy1 != paralog_ploidy2)$length)
100 * sum(filter(fcomp, paralog_ploidy1 != paralog_ploidy2)$length) / sum(fcomp$length)

sep_column <- function(vec, ix=1, sep=',') {
  suppressWarnings(as.numeric(sapply(strsplit(vec, sep), `[`, ix)))
}

extend_paralog <- function(df, copy_num) {
  for (j in 1:2) {
    for (i in 1:copy_num) {
      df[[sprintf('paralog_%d_%s', i, letters[j])]] <- 
        sep_column(df[[sprintf('paralog_ploidy%d', j)]], i)
      df[[sprintf('qual_%d_%s', i, letters[j])]] <-
        sep_column(df[[sprintf('paralog_qual%d', j)]], i)
    }
  }
  df
}

fcomp <- filter(comp, ploidy1 == ploidy2 &
                ploidy_filter1 == 'PASS' & ploidy_filter2 == 'PASS' &
                ploidy_qual1 >= QUAL & ploidy_qual2 >= QUAL &
                paralog_filter1 == 'PASS' & paralog_filter2 == 'PASS' &
                grepl(',', paralog_ploidy1))

fcomp2 <- fcomp %>% mutate(
  paralog_ploidy1 = sub(',[?0-9]*', '', paralog_ploidy1),
  paralog_ploidy2 = sub(',[?0-9]*', '', paralog_ploidy2)
)

fcomp3 <- extend_paralog(fcomp2, 5)
fcomp3_wide <- select(fcomp3, c('start', 'end', 'length', 'region', 'sample'),
                      matches('_[ab]$')) %>%
  pivot_longer(cols=matches('_[ab]$'),
               names_to=c('key', 'copy', 'dataset'),
               names_sep='_',
               values_to='value') %>%
  pivot_wider(names_from=c('key', 'dataset'))

fcomp3_wide <- filter(fcomp3_wide,
    !is.na(paralog_a) & !is.na(paralog_b) & !is.na(qual_a) & !is.na(qual_b) &
      qual_a >= QUAL & qual_b >= QUAL)
table(fcomp3_wide$copy)

sum(fcomp3_wide$length) / 1e6
sum(filter(fcomp3_wide, paralog_a != paralog_b)$length)
100 * sum(filter(fcomp3_wide, paralog_a != paralog_b)$length) / sum(fcomp3_wide$length)


region <- 'SMN1'

f_values_a <- read.csv(sprintf(
  '~/Data/hg38/jvc/runs/203.han_g1k/%s/r009/extra/em_f_values.csv', region),
                       sep='\t', comment='#')
f_values_b <- read.csv(sprintf(
  '~/Data/hg38/jvc/runs/204.han_bgi/%s/r009/extra/em_f_values.csv', region),
                       sep='\t', comment='#')

f_values_a$min_fval <- suppressWarnings(select(f_values_a, starts_with('copy')) %>%
                                          apply(1, min, na.rm=T))
f_values_b$min_fval <- suppressWarnings(select(f_values_b, starts_with('copy')) %>%
                                          apply(1, min, na.rm=T))
f_values_a$min_fval[is.infinite(f_values_a$min_fval)] <- NA
f_values_b$min_fval[is.infinite(f_values_b$min_fval)] <- NA
f_values_a$f_case <- with(f_values_a, case_when(
  is.na(min_fval) | min_fval < 0.8 ~ 'Unreliable',
  min_fval < 0.95 ~ 'Semi-reliable',
  T ~ 'Reliable'))
f_values_b$f_case <- with(f_values_b, case_when(
  is.na(min_fval) | min_fval < 0.8 ~ 'Unreliable',
  min_fval < 0.95 ~ 'Semi-reliable',
  T ~ 'Reliable'))

f_values <- full_join(f_values_a, f_values_b, c('region_group', 'psv'),
                      suffix=c('_a', '_b'))
f_values$pos <- as.numeric(sub('chr5:', '', f_values$psv))

ggplot(f_values) +
  annotate(geom='rect', xmin=0, xmax=1, ymin=0.95, ymax=1, fill='gold', alpha=.3) +
  annotate(geom='rect', xmin=0.95, xmax=1, ymin=0, ymax=1, fill='gold', alpha=.3) +
  annotate(geom='rect', xmin=0.95, xmax=1, ymin=0.95, ymax=1, fill='lawngreen', alpha=1) +
  geom_point(aes(min_fval_a, min_fval_b)) +
  geom_abline(slope=1) +
  scale_x_continuous('1000 genomes: min(f-values)') +
  scale_y_continuous('BGI samples: min(f-values)') +
  theme_bw()
ggsave('~/Tmp/1.png', width=6, height=6)

select(f_values, starts_with('f_case')) %>% table

smn_pos <- c(70925030, 70952347)
filter(f_values, pos >= smn_pos[1] & pos <= smn_pos[2]) %>%
  select(starts_with('f_case')) %>% table

f_values$ix <- 1:nrow(f_values)
ggplot(f_values, aes(ix)) +
  geom_point(aes(y=min_fval_a, color='a')) +
  geom_point(aes(y=min_fval_b, color='b'))
