library(tidyverse)
library(ggplot2)

pop_names <- c('AFR' = 'African',
               'EUR' = 'European',
               'SAS' = 'South-Asian',
               'EAS' = 'East-Asian',
               'AMR' = 'Admixed-American')

load <- function(gene, filename, keep_paralog=T, keep_qual=F, keep_b=F) {
  df <- read.csv(sprintf('%s/%s/%s.csv', data_dir, gene, filename),
                 sep='\t', comment='#') %>%
    as_tibble
  if (!keep_b) {
    df <- select(df, !matches('b_|dist'))
  }
  if (!keep_paralog) {
    df <- select(df, !matches('paralog'))
  }
  if (!keep_qual) {
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
    if ('paralog_qual' %in% colnames(df)) {
      df[[paste0('paralog_qual', i)]] <- sep_column(df$paralog_qual, i)
    }
    if ('b_paralog' %in% colnames(df)) {
      df[[paste0('b_paralog', i)]] <- sep_column(df$b_paralog, i)
    }
  }
  df$paralog_copy_num <- NULL
  df$paralog_qual <- NULL
  df$b_paralog <- NULL
  df
}

to_wide <- function(df) {
  pivot_wider(df, names_from = 'region', values_from = matches('copy_num|paralog'))
}

filter_df <- function(df, ...) {
  df_filt <- filter(df, ...)
  cat(sprintf('Removed %d samples', nrow(df) - nrow(df_filt)))
  df_filt
}

to_counts <- function(df, region_names, ..., names_prefix=NULL) {
  sample_contribution <- 100 / table(df$superpopulation)
  cn_counts <- select(df, c('sample', 'population', 'superpopulation') | all_of(...)) %>%
    pivot_longer(!c('sample', 'population', 'superpopulation'),
                 names_to='key', names_prefix) %>%
    count(superpopulation, key, value)
  cn_counts$perc <- cn_counts$n * as.numeric(sample_contribution[cn_counts$superpopulation])
  cn_counts$pop_name <- pop_names[cn_counts$superpopulation]
  cn_counts$region <- factor(region_names[cn_counts$key], levels=region_names)
  names(cn_counts$region) <- NULL
  cn_counts$value <- factor(cn_counts$value)
  cn_counts
}

ref_cns_df <- function(cn_counts, values) {
  regions <- levels(cn_counts$region)
  data.frame(region=factor(regions, levels=regions),
             cn=factor(values))
}

reorder_copies <- function(vec, ixs) {
  sapply(strsplit(vec, ','),
         function(row) do.call(paste, c(as.list(row[ixs]), sep=',')))
}

add_noise <- function(df, columns, sd=0.05) {
  n <- nrow(df)
  max_dev <- sd * 3
  for (col in columns) {
    df[[paste0(col, '_noise')]] <- df[[col]] +
      pmin(max_dev, pmax(-max_dev, rnorm(n, 0, sd)))
  }
  df
}

get_observations <- function(df, label, col1, col2, ...) {
  df_filt <- filter(df, ...)
  res <- select(df_filt, c('sample', 'population', 'superpopulation'))
  res$a_obs <- df_filt[[col1]]
  res$b_obs <- df_filt[[col2]]
  res$dist <- abs(res$a_obs - res$b_obs)
  res$label <- label
  res
}

combine_families <- function(df, pedigree) {
  df$status <- NA
  df$family <- NA
  ixs <- match(df$sample, pedigree$Individual.ID)
  jxs <- !is.na(ixs)
  df$family[jxs] <- pedigree[ixs[jxs],]$Individual.ID
  df$status[jxs] <- 'child'
  
  ixs <- match(df$sample, pedigree$Paternal.ID)
  jxs <- !is.na(ixs)
  df$family[jxs] <- pedigree[ixs[jxs],]$Individual.ID
  df$status[jxs] <- 'father'
  
  ixs <- match(df$sample, pedigree$Maternal.ID)
  jxs <- !is.na(ixs)
  df$family[jxs] <- pedigree[ixs[jxs],]$Individual.ID
  df$status[jxs] <- 'mother'
  
  df_wide <- df %>%
    select(!'sample') %>%
    pivot_wider(c('population', 'superpopulation', 'family', 'region'),
                names_from='status', values_from=matches('copy_num|paralog'))
  df_wide <- filter(df_wide, !is.na(copy_num_child) & !is.na(copy_num_father) &
                      !is.na(copy_num_mother))
  df_wide
}
