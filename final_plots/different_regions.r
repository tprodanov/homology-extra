library(ggplot2)
library(tidyverse)


data_dir <- '~/Data/hg38/jvc/comparisons/populations/r009'
plot_dir <- '~/Data/hg38/jvc/plots/population_comparison/r009/genes'

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')

perc_msg <- function(msg, x, y) {
  x <- ifelse(is.numeric(x), x, nrow(x))
  y <- ifelse(is.numeric(y), y, nrow(y))
  cat(sprintf('%-50s %d / %d (%.2f%%)\n', msg, x, y, ifelse(y > 0, x / y * 100, 0.0)))
}

load <- function(gene, filename, keep_paralog=T, keep_only_values=T) {
  df <- read.csv(sprintf('%s/%s/%s.csv', data_dir, gene, filename),
                 sep='\t', comment='#') %>%
    as_tibble
  df <- select(df, !matches('b_|dist'))
  if (!keep_paralog) {
    df <- select(df, !matches('paralog'))
  }
  if (keep_only_values) {
    df <- select(df, !matches('qual|filter'))
  }
  df
}

sep_column <- function(vec, ix=1, sep=',') {
  suppressWarnings(as.numeric(sapply(strsplit(vec, sep), `[`, ix)))
}

extend_paralog <- function(df, copy_num) {
  for (i in 1:copy_num) {
    df[[paste0('paralog', i)]] <- sep_column(df$paralog_copy_num, i)
    df[[paste0('paralog_qual', i)]] <- sep_column(df$paralog_qual, i)
  }
  df$paralog_copy_num <- NULL
  df$paralog_qual <- NULL
  df
}

to_wide <- function(df) {
  pivot_wider(df, names_from = 'region', values_from = matches('copy_num|paralog'))
}

QUAL <- 30

# ====== STRC ======

gene <- 'STRC'
df <- load(gene, 'from_pops', keep_only_values=F)
df <- filter(df, grepl('02-', region))
df$region <- recode(df$region, '02-01' = 'a', '02-02' = 'b', '02-03' = 'c')
df <- extend_paralog(df, 2)
df <- mutate(df,
  cn_passes = copy_num_filter == 'PASS' & copy_num_qual >= QUAL,
  par_passes1 = paralog_filter == 'PASS' & paralog_qual1 >= QUAL)
df_wide <- pivot_wider(df, names_from = 'region',
                       values_from = matches('copy_num|paralog|passes'))

df_cn_qual <- filter(df_wide, cn_passes_a & cn_passes_b & cn_passes_c)
perc_msg('High quality CN in all regions:', df_cn_qual, df_wide)
match_ab <- with(df_cn_qual, copy_num_a == copy_num_b)
match_ac <- with(df_cn_qual, copy_num_a == copy_num_c)
match_bc <- with(df_cn_qual, copy_num_b == copy_num_c)
correct_cn <- df_cn_qual[match_ab & match_ac & match_bc,]
do.call(paste, as.list(df_cn_qual[!(match_ab & match_ac & match_bc),]$sample))

perc_msg('CN match in regions A & B:', sum(match_ab), df_cn_qual)
perc_msg('CN match in all regions:', correct_cn, df_cn_qual)

df_par_qual <- filter(correct_cn, par_passes1_a & par_passes1_b & par_passes1_c &
                        !is.na(paralog1_a) & !is.na(paralog1_b) & !is.na(paralog1_c))
perc_msg('High quality PCN in all regions:', df_par_qual, correct_cn)
pmatch_ab <- with(df_par_qual, paralog1_a == paralog1_b)
pmatch_ac <- with(df_par_qual, paralog1_a == paralog1_c)
pmatch_bc <- with(df_par_qual, paralog1_b == paralog1_c)
correct_pcn <- df_par_qual[pmatch_ab & pmatch_ac & pmatch_bc,]
perc_msg('PCN match in regions A & B:', sum(pmatch_ab), df_par_qual)
perc_msg('PCN match in regions A & C:', sum(pmatch_ac), df_par_qual)
perc_msg('PCN match in regions B & C:', sum(pmatch_bc), df_par_qual)
perc_msg('PCN match in all regions:', correct_pcn, df_par_qual)

View(filter(df_par_qual, paralog1_a != paralog1_b))

do.call(paste, as.list(df_par_qual[!pmatch_bc,]$sample))
