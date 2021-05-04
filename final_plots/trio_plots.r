wdir <- '~/Code/homology-extra/final_plots/'
source(sprintf('%s/common.r', wdir))
plot_dir <- '~/Data/hg38/jvc/plots/population_comparison/r009/trio/'

QUAL = 30
pedigree <- read.csv('~/Data/hg38/data/g1k/g1k.ped', sep='\t')

# === GTF2I ===

df <- (rbind(
  read.csv('~/Data/hg38/jvc/comparisons/children/r009/GTF2I/children.csv',
         sep='\t', comment.char='#'),
  read.csv('~/Data/hg38/jvc/comparisons/populations/r009/GTF2I/from_pops.csv',
           sep='\t', comment.char='#')) %>%
  select(!matches('b_|dist')) %>% filter(region != 'gene_start'))
df <- extend_paralog(df, 3)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2', 'paralog3'))
df <- mutate(df,
               cn_passes=copy_num_qual >= QUAL & copy_num_filter == 'PASS',
               par1_passes=paralog_qual1 >= QUAL & paralog_filter == 'PASS',
               par2_passes=paralog_qual2 >= QUAL & paralog_filter == 'PASS',
               par3_passes=paralog_qual3 >= QUAL & paralog_filter == 'PASS')

df_fam <- combine_families(df, pedigree)
df_fam$label <- sprintf('Individual CN = %d', df_fam$copy_num_child)

ggplot(df_fam) +
  geom_count(aes(copy_num_noise_mother, copy_num_noise_father),
             alpha=.5, size=1.5) +
  facet_wrap(~ label) +
  theme_bw() +
  scale_x_continuous('Maternal CN') +
  scale_y_continuous('Paternal CN') +
  ggtitle('GTF2I: Aggregate copy number')

# ====== APOBEC3A ======

df <- (rbind(
  read.csv('~/Data/hg38/jvc/comparisons/children/r009/APOBEC3A/children.csv',
           sep='\t', comment.char='#'),
  read.csv('~/Data/hg38/jvc/comparisons/populations/r009/APOBEC3A/from_pops.csv',
           sep='\t', comment.char='#')) %>%
  select(!matches('b_|dist')) %>% filter(region != 'gene_start'))

df <- extend_paralog(df, 2)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'), sd=0.03)
df <- mutate(df,
             cn_passes=copy_num_qual >= QUAL & copy_num_filter == 'PASS',
             par1_passes=paralog_qual1 >= QUAL & paralog_filter == 'PASS')

df_fam <- combine_families(df, pedigree)
df_fam$label <- sprintf('Individual CN = %d', df_fam$copy_num_child)

ggplot(df_fam) +
  geom_count(aes(copy_num_noise_mother, copy_num_noise_father),
             alpha=.6, size=1) +
  facet_wrap(~ label) +
  theme_bw() +
  scale_x_continuous('Maternal CN') +
  scale_y_continuous('Paternal CN') +
  ggtitle('APOBEC3A + APOBEC3B')
ggsave(sprintf('%s/APOBEC3A.aggr.png', plot_dir), width=11, height=6)

df_fam$par_label <- sprintf('Individual CN = %d', df_fam$paralog1_child)
df_filt <- filter(df_fam,
    paralog_filter_child == 'PASS' & paralog_filter_father == 'PASS' &
      paralog_filter_mother == 'PASS' & paralog_qual1_child >= QUAL &
      paralog_qual1_father >= QUAL & paralog_qual1_mother >= QUAL)


ggplot(df_filt) +
  geom_count(aes(paralog1_noise_mother, paralog1_noise_father),
             alpha=.6, size=1) +
  facet_wrap(~ par_label) +
  theme_bw() +
  scale_x_continuous('Maternal CN') +
  scale_y_continuous('Paternal CN') +
  ggtitle('APOBEC3A')
ggsave(sprintf('%s/APOBEC3A.copy1.png', plot_dir), width=11, height=6)

# ====== STRC ======

df <- (rbind(
  read.csv('~/Data/hg38/jvc/comparisons/children/r009/STRC/children.csv',
           sep='\t', comment.char='#'),
  read.csv('~/Data/hg38/jvc/comparisons/populations/r009/STRC/from_pops.csv',
           sep='\t', comment.char='#')) %>%
  select(!matches('b_|dist')) %>% filter(region == 'gene_start'))
df <- extend_paralog(df, 2)
df <- add_noise(df, c('copy_num', 'paralog1', 'paralog2'), sd=0.03)
df <- mutate(df,
             cn_passes=copy_num_qual >= QUAL & copy_num_filter == 'PASS',
             par1_passes=paralog_qual1 >= QUAL & paralog_filter == 'PASS')

df_fam <- combine_families(df, pedigree)
df_fam$label <- sprintf('Individual CN = %d', df_fam$copy_num_child)
df_fam <- mutate(df_fam,
                 min_cn_qual = pmin(copy_num_qual_child, copy_num_qual_father,
                                    copy_num_qual_mother)) %>%
  mutate(cn_qual_case = factor(case_when(
    min_cn_qual < 30 ~ 'Low',
    min_cn_qual < 200 ~ 'Interm.',
    min_cn_qual >= 200 ~ 'High'),
    levels=c('Low', 'Interm.', 'High')
  ))

ggplot(filter(df_fam, min_cn_qual >= 30)) +
  geom_count(aes(copy_num_noise_mother, copy_num_noise_father),
             alpha=.6, size=2) +
  facet_wrap(~ label) +
  theme_bw() +
  scale_x_continuous('Maternal CN') +
  scale_y_continuous('Paternal CN') +
  ggtitle('STRC: Aggregate copy number')
ggsave(sprintf('%s/STRC.aggr.png', plot_dir), width=11, height=6)

(mutate(df_fam, cn_prof = sprintf(
  '%d,%d,%d', copy_num_child, copy_num_mother, copy_num_father)) %>%
  count(cn_prof))
nrow(df_fam)

df_fam$par_label <- sprintf('Individual CN = %d', df_fam$paralog1_child)
df_fam <- mutate(df_fam,
                 min_pcn_qual = pmin(paralog_qual1_child, paralog_qual1_father,
                                     paralog_qual1_mother),
                 par_passes_all = paralog_filter_child == 'PASS' &
                   paralog_filter_mother == 'PASS' & paralog_filter_father == 'PASS') %>%
  mutate(pcn_qual_case = factor(case_when(
    min_pcn_qual < 30 | !par_passes_all ~ 'Low',
    min_pcn_qual < 200 ~ 'Interm.',
    min_pcn_qual >= 200 ~ 'High'),
    levels=c('Low', 'Interm.', 'High')
  ))

df_filt <- filter(df_fam, min_cn_qual >= 30 & min_pcn_qual >= 30 & par_passes_all)
ggplot(df_filt) +
  geom_count(aes(paralog1_noise_mother, paralog1_noise_father),
             alpha=.6, size=2) +
  facet_wrap(~ par_label) +
  theme_bw() +
  scale_x_continuous('Maternal CN') +
  scale_y_continuous('Paternal CN') +
  ggtitle('STRC')
ggsave(sprintf('%s/STRC.copy1.png', plot_dir), width=11, height=6)

(mutate(df_filt, pcn_prof = sprintf(
  '%d,%d,%d', paralog1_child, paralog1_mother, paralog1_father)) %>%
    count(pcn_prof))
nrow(df_filt)
