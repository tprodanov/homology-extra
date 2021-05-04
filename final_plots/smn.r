library(tidyverse)
library(ggplot2)

setwd('~/Data/hg38/jvc/comparisons/populations/r009/')

perc_msg <- function(msg, x, y) {
  x <- ifelse(is.numeric(x), x, nrow(x))
  y <- ifelse(is.numeric(y), y, nrow(y))
  cat(sprintf('%-50s %d / %d (%.2f%%)\n', msg, x, y, ifelse(y > 0, x / y * 100, 0.0)))
}

load <- function(gene, filename, keep_paralog=T) {
  df <- read.csv(sprintf('%s/%s.csv', gene, filename), sep='\t', comment='#') %>%
    as_tibble
  df <- select(df, !matches('b_|dist'))
  if (!keep_paralog) {
    df <- select(df, !matches('paralog'))
  }
  df %>% rename('region' = 'method')
}

# ====== SMN1 ======

# ------ SMN1 Barplots ------

gene <- 'SMN1'
df <- load(gene, 'from_all')


smn_data <- read.csv('SMN1/from_all_with_none.csv', sep = '\t', comment = '#') %>% as_tibble
smn_data <- select(smn_data, !matches('b_|dist'))
smn_data <- filter(smn_data, grepl('SMN1', method))
smn_data$paralog1 <- as.numeric(sapply(strsplit(smn_data$paralog_copy_num, ','),
                                       `[`, 1))
smn_data$paralog2 <- as.numeric(sapply(strsplit(smn_data$paralog_copy_num, ','),
                                       `[`, 2))
smn_data_wide <- smn_data %>%
  mutate(exons = ifelse(method == 'SMN1_middle', '16', '78')) %>%
  select(c('sample', 'copy_num', 'exons') | matches('population|paralog[12]')) %>%
  pivot_wider(names_from = 'exons', values_from = c('copy_num', 'paralog1', 'paralog2'))

smn_data_filt <- filter(smn_data_wide, !is.na(paralog1_78))
cat(sprintf('Removed %d samples with unknown paralog copy number\n',
            nrow(smn_data_wide) - nrow(smn_data_filt)))

smn_data_filt <- mutate(smn_data_filt,
  cn1 = paralog1_78,
  cn2 = paralog2_78,
  cn2_delta = copy_num_16 - paralog1_78 - paralog2_78)
sample_contribution <- 100 / table(smn_data_filt$superpopulation)

cn_counts <- smn_data_filt %>%
  select(c('sample') | ends_with('population') | starts_with('cn')) %>%
  pivot_longer(starts_with('cn'), names_to='key') %>%
  count(superpopulation, key, value)
cn_counts$perc <- cn_counts$n * as.numeric(sample_contribution[cn_counts$superpopulation])

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')
cn_counts$pop_name <- pop_names[cn_counts$superpopulation]

region_names <- c('cn1' = 'SMN1',
                  'cn2' = 'SMN2',
                  'cn2_delta' = 'SMN2Δ7-8')
cn_counts$region <- region_names[cn_counts$key]

ggplot(cn_counts) +
  geom_bar(aes(factor(value), perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single'), width=.6) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_brewer('', palette = 'Dark2') +
  theme_bw()

# ------ SERF1A Barplots ------

serf1a_data <- read.csv('SMN1/with_none.csv', sep = '\t', comment = '#') %>% as_tibble
serf1a_data <- select(serf1a_data, !matches('paralog|b_|dist'))
serf1a_data <- filter(serf1a_data, grepl('SERF1A', method))

serf1a_wide <- serf1a_data %>%
  mutate(region = ifelse(method == 'SERF1A_middle', 'middle', 'end')) %>%
  select(c('sample', 'copy_num', 'region') | ends_with('population')) %>%
  pivot_wider(names_from = 'region', values_from = 'copy_num')
serf1a_wide$n_del <- with(serf1a_wide, middle - end)

cat(sprintf('There are %d strange sample(s)\n', sum(serf1a_wide$n_del < 0)))
serf1a_wide <- filter(serf1a_wide, n_del >= 0 & middle <= 6)
sample_contribution <- 100 / table(serf1a_wide$superpopulation)

cn_counts <- serf1a_wide %>%
  select(c('sample', 'middle', 'n_del') | ends_with('population')) %>%
  pivot_longer(c('middle', 'n_del'), names_to='key') %>%
  count(superpopulation, key, value)
cn_counts$perc <- cn_counts$n * as.numeric(sample_contribution[cn_counts$superpopulation])

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')
cn_counts$pop_name <- pop_names[cn_counts$superpopulation]

region_names <- c('middle' = 'SERF1A',
                  'n_del' = 'Deletion')
cn_counts$region <- factor(region_names[cn_counts$key],
                           levels=region_names)

ggplot(cn_counts) +
  geom_bar(aes(factor(value), perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single'), width=.6) +
  #facet_grid(. ~ region, scales='free_x', space='free_x') +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_brewer('', palette = 'Dark2') +
  theme_bw()

# ------ C4A barplots ------

c4a_data <- read.csv('C4A/with_none.csv', sep = '\t', comment = '#') %>% as_tibble
c4a_data <- select(c4a_data, !matches('paralog|b_|dist'))

c4a_wide <- c4a_data %>%
  rename('region' = 'method') %>%
  select(c('sample', 'copy_num', 'region') | ends_with('population')) %>%
  pivot_wider(names_from = 'region', values_from = 'copy_num')
c4a_wide$n_del <- with(c4a_wide, C4A_end - C4A_middle)

cat(sprintf('There are %d sample(s) with insertion\n', sum(c4a_wide$n_del < 0)))
cat(sprintf('There are %d sample(s) with different copy num before and after deletion\n',
            sum(c4a_wide$C4A_start != c4a_wide$C4A_end)))
c4a_wide <- filter(c4a_wide, n_del >= 0 & C4A_start == C4A_end)
sample_contribution <- 100 / table(c4a_wide$superpopulation)

cn_counts <- c4a_wide %>%
  select(c('sample', 'C4A_start', 'n_del') | ends_with('population')) %>%
  pivot_longer(c('C4A_start', 'n_del'), names_to='key') %>%
  count(superpopulation, key, value)
cn_counts$perc <- cn_counts$n * as.numeric(sample_contribution[cn_counts$superpopulation])

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')
cn_counts$pop_name <- pop_names[cn_counts$superpopulation]

region_names <- c('C4A_start' = 'C4A',
                  'n_del' = 'Deletion')
cn_counts$region <- factor(region_names[cn_counts$key],
                           levels=region_names)

ggplot(cn_counts) +
  geom_bar(aes(factor(value), perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single'), width=.6) +
  #facet_grid(. ~ region, scales='free_x', space='free_x') +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_brewer('', palette = 'Dark2') +
  theme_bw()

# ------ FAM185A barplots ------

fam <- read.csv('FAM185A/with_none.csv', sep = '\t', comment = '#') %>% as_tibble
fam <- select(fam, !matches('paralog|b_|dist'))

fam_wide <- fam %>%
  rename('region' = 'method') %>%
  select(c('sample', 'copy_num', 'region') | ends_with('population')) %>%
  pivot_wider(names_from = 'region', values_from = 'copy_num')
fam_wide$n_del <- with(fam_wide, after - del)

cat(sprintf('There are %d sample(s) with insertion\n', sum(fam_wide$n_del < 0)))
cat(sprintf('There are %d sample(s) with different copy num before and after deletion\n',
            sum(fam_wide$after != fam_wide$before)))
fam_wide <- filter(fam_wide, n_del >= 0 & before == after)
sample_contribution <- 100 / table(fam_wide$superpopulation)

cn_counts <- fam_wide %>%
  select(c('sample', 'before', 'n_del') | ends_with('population')) %>%
  pivot_longer(c('before', 'n_del'), names_to='key') %>%
  count(superpopulation, key, value)
cn_counts$perc <- cn_counts$n * as.numeric(sample_contribution[cn_counts$superpopulation])

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')
cn_counts$pop_name <- pop_names[cn_counts$superpopulation]

region_names <- c('before' = 'FAM185A',
                  'n_del' = 'Deletion')
cn_counts$region <- factor(region_names[cn_counts$key],
                           levels=region_names)

ggplot(cn_counts) +
  geom_bar(aes(factor(value), perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single'), width=.6) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_brewer('', palette = 'Dark2') +
  theme_bw()

# ------ NBPF4 barplots ------

nbpf <- read.csv('NBPF4/with_none.csv', sep = '\t', comment = '#') %>% as_tibble
nbpf <- select(nbpf, !matches('paralog|b_|dist'))

nbpf_wide <- nbpf %>%
  rename('region' = 'method') %>%
  select(c('sample', 'copy_num', 'region') | ends_with('population')) %>%
  pivot_wider(names_from = 'region', values_from = 'copy_num')

sample_contribution <- 100 / table(nbpf_wide$superpopulation)

cn_counts <- nbpf_wide %>%
  select(c('sample', 'start', 'middle', 'end') | ends_with('population')) %>%
  pivot_longer(c('start', 'middle', 'end'), names_to='key') %>%
  count(superpopulation, key, value)
cn_counts$perc <- cn_counts$n * as.numeric(sample_contribution[cn_counts$superpopulation])

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')
cn_counts$pop_name <- pop_names[cn_counts$superpopulation]

region_names <- c('start' = 'exons 10-14',
                  'middle' = '~ exon 5',
                  'end' = '~ exon 1')
cn_counts$region <- factor(region_names[cn_counts$key],
                           levels=region_names)

ggplot(cn_counts) +
  geom_bar(aes(factor(value), perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single'), width=.6) +
  facet_grid(region ~ ., scales='free_x') +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_brewer('', palette = 'Dark2') +
  theme_bw()

# ------ EIF3C barplots ------

eif3c <- read.csv('EIF3C/with_none.csv', sep = '\t', comment = '#') %>% as_tibble
eif3c <- select(eif3c, !matches('paralog|b_|dist'))

eif3c_wide <- eif3c %>%
  rename('region' = 'method') %>%
  select(c('sample', 'copy_num', 'region') | ends_with('population')) %>%
  pivot_wider(names_from = 'region', values_from = 'copy_num')

sample_contribution <- 100 / table(nbpf_wide$superpopulation)

cn_counts <- eif3c_wide %>%
  select(c('sample', 'start', 'middle', 'end') | ends_with('population')) %>%
  pivot_longer(c('start', 'middle', 'end'), names_to='key') %>%
  count(superpopulation, key, value)
cn_counts$perc <- cn_counts$n * as.numeric(sample_contribution[cn_counts$superpopulation])

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')
cn_counts$pop_name <- pop_names[cn_counts$superpopulation]

region_names <- c('start' = 'Exon 1',
                  'middle' = 'Intron 1',
                  'end' = 'Exon 9')
cn_counts$region <- factor(region_names[cn_counts$key],
                           levels=region_names)

ggplot(cn_counts) +
  geom_bar(aes(factor(value), perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single'), width=.6) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_brewer('', palette = 'Dark2') +
  theme_bw()

# ------ SRGAP2 barplots ------

srgap <- read.csv('SRGAP2/with_none.csv', sep = '\t', comment = '#') %>% as_tibble
srgap <- select(srgap, !matches('b_|dist'))
srgap$paralog1 <- as.numeric(sapply(strsplit(srgap$paralog_copy_num, ','),
                                    `[`, 1))
srgap$paralog2 <- as.numeric(sapply(strsplit(srgap$paralog_copy_num, ','),
                                    `[`, 2))
srgap$paralog3 <- as.numeric(sapply(strsplit(srgap$paralog_copy_num, ','),
                                    `[`, 3))
srgap$paralog4 <- as.numeric(sapply(strsplit(srgap$paralog_copy_num, ','),
                                    `[`, 4))
srgap_wide <- srgap %>%
  rename(region = method) %>%
  select(c('sample', 'copy_num', 'region') | matches('population|paralog[1-4]')) %>%
  pivot_wider(names_from = 'region',
              values_from = 'copy_num' | matches('paralog[1-4]')) %>%
  select(!c('paralog4_c', 'paralog4_d'))

srgap_filt <- filter(srgap_wide, copy_num_a == copy_num_b & copy_num_c == copy_num_d)
cat(sprintf('Removed %d samples with varying copy number\n',
            nrow(srgap_wide) - nrow(srgap_filt)))

#smn_data_filt <- mutate(smn_data_filt,
                        #cn1 = paralog1_78,
                        #cn2 = paralog2_78,
                        #cn2_delta = copy_num_16 - paralog1_78 - paralog2_78)
sample_contribution <- 100 / table(srgap_filt$superpopulation)

cn_counts <- srgap_filt %>%
  pivot_longer(starts_with('copy_num'), names_to='key') %>%
  count(superpopulation, key, value)
cn_counts$perc <- cn_counts$n * as.numeric(sample_contribution[cn_counts$superpopulation])

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')
cn_counts$pop_name <- pop_names[cn_counts$superpopulation]

region_names <- c('cn1' = 'SMN1',
                  'cn2' = 'SMN2',
                  'cn2_delta' = 'SMN2Δ7-8')
cn_counts$region <- region_names[cn_counts$key]

ggplot(cn_counts) +
  geom_bar(aes(factor(value), perc, fill=pop_name), stat='identity',
           position=position_dodge2(preserve = 'single'), width=.6) +
  facet_grid(region ~ .) +
  scale_x_discrete('Copy number') +
  scale_y_continuous('Percentage') +
  scale_fill_brewer('', palette = 'Dark2') +
  theme_bw()

# ------ SMN caller ------

smn_caller <- read.csv('smn1_caller.csv', sep = '\t', comment = '#') %>% as_tibble

compare_smn_caller <- function(smn_caller, qual=30) {
  smn_caller$paralog_qual_1 <- suppressWarnings(as.numeric(
    sapply(strsplit(smn_caller$paralog_qual, ','), `[`, 1)))
  
  df_16 <- filter(smn_caller, method == 'caller_16')
  df_78 <- filter(smn_caller, method == 'caller_78')
  df_join <- inner_join(df_16, df_78, by=c('sample', 'population'))

  cn_present <- filter(df_join,
    copy_num.x != '*' & copy_num.y != '*' & b_copy_num.x != '*' & b_copy_num.y != '*')
  cn_eq = filter(cn_present,
    copy_num.x == round(b_copy_num.x) & copy_num.y == round(b_copy_num.y))
  perc_msg('Copy number concordance:', cn_eq, cn_present)

  cn_present_q <- filter(cn_present,
    copy_num_filter.x == 'PASS' & copy_num_filter.y == 'PASS' &
    copy_num_qual.x >= qual & copy_num_qual.y >= qual)
  cn_eq_q <- filter(cn_present_q,
    copy_num.x == round(b_copy_num.x) & copy_num.y == round(b_copy_num.y))
  perc_msg('Copy number concordance (high quality):', cn_eq_q, cn_present_q)
  
  par_present <- filter(df_78, paralog_copy_num != '?,?' & b_paralog != '*')
  par_eq <- filter(par_present, paralog_copy_num == b_paralog)
  perc_msg('Paralog copy number concordance:', par_eq, par_present)
  
  par_present_q <- filter(par_present, paralog_filter == 'PASS' & paralog_qual_1 >= qual)
  par_eq_q <- filter(par_present_q, paralog_copy_num == b_paralog)
  perc_msg('Paralog copy number concordance (high quality):', par_eq_q, par_present_q)
  
  extra_pc <- filter(df_78, paralog_copy_num != '?,?' & b_paralog == '*')
  extra_pc_q <- filter(extra_pc, paralog_filter == 'PASS' & paralog_qual_1 >= qual)
  extra_smnc <- filter(df_78, paralog_copy_num == '?,?' & b_paralog != '*')
  perc_msg('Missing in SMNCaller', extra_pc, df_78)
  #perc_msg('Found only using Parascopy (high quality)', extra_pc_q, df_78)
  perc_msg('Missing in Parascope', extra_smnc, df_78)
}

compare_smn_caller(smn_caller)

# ------ MLPA ------

smn1_mlpa <- read.csv('smn1_mlpa.csv', sep = '\t', comment = '#') %>% as_tibble

draw_mlpa <- function(smn1_mlpa, smn_caller, qual=30) {
  mlpa_16 <- filter(smn1_mlpa, method == 'mlpa_16')
  mlpa_78 <- filter(smn1_mlpa, method == 'mlpa_78')
  mlpa_joined <- full_join(mlpa_16, mlpa_78, by=c('sample', 'population'))

  cn_present <- mlpa_joined
  cn_eq <- filter(cn_present,
    copy_num.x == round(b_copy_num.x) & copy_num.y == round(b_copy_num.y))
  perc_msg('Copy number concordance:', cn_eq, cn_present)

  cn_wo_beb_tsi <- filter(cn_present, population != 'BEB' & population != 'TSI')
  cn_eq_wo_beb_tsi <- filter(cn_wo_beb_tsi,
    copy_num.x == round(b_copy_num.x) & copy_num.y == round(b_copy_num.y))
  perc_msg('Copy number concordance (w/o BEB & TSI):', cn_eq_wo_beb_tsi, cn_wo_beb_tsi)

  mlpa_smnc_16 <- left_join(mlpa_16, filter(smn_caller, method == 'caller_16'),
    by = c('sample', 'population'))
  mlpa_smnc_78 <- left_join(mlpa_78, filter(smn_caller, method == 'caller_78'),
    by = c('sample', 'population'))
  mlpa_smnc <- rbind(add_column(mlpa_smnc_16, region = 'Exons 1-6'),
                     add_column(mlpa_smnc_78, region = 'Exons 7-8'))

  ggplot(mlpa_smnc) +
    geom_abline(slope=1, color='gray80') +
    geom_point(aes(b_copy_num.y, b_copy_num.x, color=as.character(copy_num.x)), alpha=.8) +
    theme_bw() +
    facet_wrap(~ region) +
    scale_x_continuous('SMNCopyNumberCaller copy number estimate') +
    scale_y_continuous('MLPA copy number estimate') +
    scale_colour_brewer('Parascopy\nCN estimate', palette = "Dark2") +
    guides(color = guide_legend(override.aes = list(alpha=1, size=2.5)))
}

draw_mlpa(smn1_mlpa, smn_caller)

# ------ Quick-Mer2 ------

smn1_qm2 <- read.csv('smn1_qm2.csv', sep = '\t', comment = '#') %>% as_tibble

draw_qm2 <- function(smn1_qm2, smn_caller, qual=30) {
  wo_del <- (full_join(filter(smn_caller, method == 'caller_16'),
                       filter(smn_caller, method == 'caller_78'), by='sample') %>%
    filter(copy_num.x == copy_num.y))$sample
  qm2_wodel <- filter(smn1_qm2, sample %in% wo_del)
  
  cn_present <- qm2_wodel
  cn_eq <- filter(cn_present, copy_num == round(b_copy_num))
  perc_msg('Total copy number concordance:', cn_eq, cn_present)

  cn_present_q <- filter(qm2_wodel, copy_num_filter == 'PASS' & copy_num_qual >= qual)
  cn_eq_q <- filter(cn_present_q, copy_num == round(b_copy_num))
  perc_msg('Total copy number concordance (high quality):', cn_eq_q, cn_present_q)

  par_present <- filter(qm2_wodel, paralog_copy_num != '?,?' & b_paralog != '*')
  par_present$b_paralog_round <- lapply(strsplit(par_present$b_paralog, ','),
    function(x) {
      x <- round(as.numeric(x))
      sprintf('%d,%d', x[1], x[2])
  }) %>% unlist
  par_eq <- filter(par_present, paralog_copy_num == b_paralog_round)
  perc_msg('Paralog copy number concordance:', par_eq, par_present)
  
  qm2_smnc <- left_join(par_present, filter(smn_caller, method == 'caller_16'),
                        by=c('sample', 'population'))
  qm2_smnc_1 <- add_column(qm2_smnc, region = 'SMN1')
  qm2_smnc_1$qm2_cn <- as.numeric(sapply(strsplit(qm2_smnc_1$b_paralog.x, ','), `[`, 1))
  qm2_smnc_1$smnc_cn <- suppressWarnings(as.numeric(sapply(
    strsplit(qm2_smnc_1$b_paralog.y, ','), `[`, 1)))
  qm2_smnc_1$pc_cn <- as.numeric(sapply(
    strsplit(qm2_smnc_1$paralog_copy_num.x, ','), `[`, 1))

  qm2_smnc_2 <- add_column(qm2_smnc, region = 'SMN2')
  qm2_smnc_2$qm2_cn <- as.numeric(sapply(strsplit(qm2_smnc_2$b_paralog.x, ','), `[`, 2))
  qm2_smnc_2$smnc_cn <- suppressWarnings(as.numeric(sapply(
    strsplit(qm2_smnc_2$b_paralog.y, ','), `[`, 2)))
  qm2_smnc_2$pc_cn <- as.numeric(sapply(
    strsplit(qm2_smnc_2$paralog_copy_num.x, ','), `[`, 2))

  qm2_smnc_both <- rbind(qm2_smnc_1, qm2_smnc_2)
  ggplot(qm2_smnc_both) +
    geom_abline(slope=1, color='gray80') +
    geom_point(aes(smnc_cn, qm2_cn, color=as.character(pc_cn)), alpha=.8) +
    facet_wrap(~ region) +
    theme_bw() +
    scale_x_continuous('SMNCopyNumberCaller copy number estimate') +
    scale_y_continuous('QuicK-mer2 copy number estimate') +
    scale_colour_brewer('Parascopy\nCN estimate', palette = "Dark2") +
    guides(color = guide_legend(override.aes = list(alpha=1, size=2.5)))
}

draw_qm2(smn1_qm2, smn_caller)

# ====== AMY1C ======

amy1c <- read.csv('amy1c_qpcr.csv', sep = '\t', comment = '#') %>% as_tibble
amy1c <- select(amy1c, !contains('paralog')) %>% select(!total_dist)
aggregate(b_copy_num ~ method, amy1c, function(x) sum(!is.na(x)))

amy1c_wide <- pivot_wider(amy1c, names_from = 'method', values_from = 'b_copy_num')
ggplot(amy1c_wide) +
  geom_abline(slope=1, color='gray80') +
  geom_point(aes(PRT, qPCR, color=factor(copy_num)), alpha=.8) +
  theme_bw() +
  scale_x_continuous('PRT copy number estimate') +
  scale_y_continuous('qPCR copy number estimate') +
  # scale_colour_brewer('Parascopy\nCN estimate', palette = "Dark2") +
  guides(color = guide_legend(override.aes = list(alpha=1, size=2.5)))

ggplot(amy1c_wide) +
  geom_abline(slope=1, color='gray80') +
  geom_boxplot(aes(copy_num, qPCR, group=copy_num), alpha=.8) +
  geom_jitter(aes(copy_num, qPCR, group=copy_num), alpha=.8, width=.2) +
  theme_bw() +
  scale_x_continuous('Parascopy copy number estimate') +
  scale_y_continuous('qPCR copy number estimate') +
  # scale_colour_brewer('Parascopy\nCN estimate', palette = "Dark2") +
  guides(color = guide_legend(override.aes = list(alpha=1, size=2.5)))

# ====== NPY4R ======

npy4r <- read.csv('npy4r.csv', sep = '\t', comment = '#') %>% as_tibble
npy4r <- select(npy4r, !contains('paralog')) %>% select(!total_dist)
aggregate(b_copy_num ~ method, npy4r, function(x) sum(!is.na(x)))

npy4r_wide <- pivot_wider(npy4r, names_from = 'method', values_from = 'b_copy_num')
ggplot(npy4r_wide) +
  geom_abline(slope=1, color='gray80') +
  geom_point(aes(CNVnator, ddPCR, color=factor(copy_num)), alpha=.8) +
  theme_bw() +
  scale_x_continuous('CNVnator copy number estimate') +
  scale_y_continuous('ddPCR copy number estimate') +
  scale_colour_brewer('Parascopy\nCN estimate', palette = "Dark2") +
  guides(color = guide_legend(override.aes = list(alpha=1, size=2.5)))

ggplot(npy4r_wide) +
  geom_abline(slope=1, color='gray80') +
  geom_point(aes(CNVnator, FREEC, color=factor(copy_num)), alpha=.8) +
  theme_bw() +
  scale_x_continuous('CNVnator copy number estimate') +
  scale_y_continuous('FREEC copy number estimate') +
  scale_colour_brewer('Parascopy\nCN estimate', palette = "Dark2") +
  guides(color = guide_legend(override.aes = list(alpha=1, size=2.5)))
